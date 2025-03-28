/*
Copyright (C) 2016-2023 Deep Genomics Inc. All Rights Reserved.
*/
#include "sam_line_parser.h"

#include "genome.h"
#include "gk_assert.h"
#include "strutil.h"
#include "util.h"
#include <cstdlib>

using namespace std;

BEGIN_NAMESPACE_GK

///////////////////////////////////////////////////////////////
// PARSE CIGAR STRING
//
// - REFERENCE GENOME vs PRIVATE GENOME
//   Junction positions are on the REFERENCE genome, NOT the
//   private genome, for obvious reasons.
//
// - COUNTING INSERTIONS/DELETIONS
//   When counting matches on the PRIVATE genome, one must
//   include 'I' insertions but exclude 'D' deletions.
//   When counting matches on the REFERENCE genome, it's
//   the opposite. So, we must sum the matches two separate
//   ways in order to calculate both the junction start/end
//   (reference genome) and the overhangs (alignment matches only).
//
//   For example, a read at start=2 had CIGAR string 1M1D1M2N1M
//   matches the reference genome as follows:
//
//      0 1 2 3 4 5 6 7 8 9     <-- reference genome coordinates
//          M D M - - M
//
//   so the start of the junction is 2+1+1+1 = 5, but the
//   the left overhang is only 2, with the D not being counted,
//   since that base wasn't there in the private genome.
//
//   Similarly for insertions, with CIGAR string 1M1I1M2N1M:
//
//      0 1 2 3 4 5 6 7 8 9
//          M M - - M
//           ^
//           I
//
//   the start of the junction is 2+1+1 = 4, but the left
//   overhang is actually 2, since
//   that base was there only in the private genome.
//
// - MULTI-SPLIT MAPPINGS
//   If a CIGAR string spans multiple junctions, then a junction
//   entry is created for each split. That means the read contributes
//   to multiple junctions. The left/right overhang for each junction
//   will the sum of all matches before/after that particular junction,
//   leaving out any gaps caused by splits (Ns) before/after the junction
//   in question.
//
//   For example, the CIGAR string 2M2N1M1N1M corresponds to true alignment
//
//      0 1 2 3 4 5 6 7 8 9
//          M M - - M - M
//
//   but for the junction table the overhangs of this read in its respective
//   junctions will be as if it were two reads aligned as follows:
//
//      0 1 2 3 4 5 6 7 8 9
//          M M - - M M
//              M M M - M
//
//   Obviously the read sequence doesn't actually match the reference genome
//   at the above positions, but this scheme facilitates counting of relative
//   junction usage (assuming usage of the two junctions is independent)
//   and also is allows positional bootstrap to be applied without modification
//   (since PCR duplicates will pile up in both junctions, with no positional
//    bias exon start/end)
//
//  - CIGAR types
//      consumes query  consumes reference  description
//  M       yes             yes             alignment match (sequence match or unmatch)
//  I       yes             no              insertion to reference
//  D       no              yes             deletion from reference
//  N       no              yes             skipped region from reference
//  S       yes             no              soft clipping (clipped sequence in SEQ)
//
//////////////////////////////////////////////////////////////////////
static void parse_next_cigar(string_view& cigar, int& len, char& code)
{
	while (!empty(cigar)) {
		auto pos = cigar.find_first_not_of("0123456789");
		if (pos != string_view::npos) {
			code = cigar[pos];
			len  = as_int(cigar.substr(0, pos));
			GK_CHECK(len != 0, value, "Could not convert '{}' to integer in CIGAR string", cigar);
			cigar.remove_prefix(pos + 1); // advance cigar pointer to next start
			return;
		}
		cigar.remove_prefix(size(cigar));
	}
}

void sam_line_parser::detect_strand_with_library(const char* format)
{
	// see https://salmon.readthedocs.io/en/stable/library_type.html
	if (endswith(format, "SF")) {
		_first_read_direction = direction::forward;
	} else if (endswith(format, "SR")) {
		_first_read_direction = direction::reverse;
	} else if (endswith(format, "U")) {
		_first_read_direction = direction::unknown;
	} else {
		GK_THROW(value, "Unknown library format string {}", format);
	}
}

void sam_line_parser::process_file(line_reader& lr)
{
	sam_cols cols;
	constexpr int max_cigar_codes = 200; // 53 block count for AF150400 EST
	cigar_op codes[max_cigar_codes];

	for (; !lr.done(); ++lr) {
		auto line = lr.line();

		// Skip empty lines, comment lines, any line not starting with "chr", or any line for the mitochondrial DNA
		if (empty(line) || line[0] == '@') {
			continue;
		}

		// Split tab-separated level. The SAM fields at this level are
		// [0] QNAME
		// [1] FLAG
		// [2] RNAME
		// [3] POSITION (1-based, leftmost mapping)
		// [4] MAPQ
		// [5] CIGAR
		// [6] RNEXT
		// [7] PNEXT
		// [8] TLEN
		// [9] SEQ
		// [10] QUAL
		// [11] Optional fields
		GK_CHECK(split_view(line, '\t', cols, num_cols) >= num_cols - 1, value,
				 "Expected at least {} tab-separated columns", num_cols - 1);

		// Check SAM flags
		unsigned flags = as_int(cols[1]);
		if ((flags & 0x0900) != 0) {
			continue; // A secondary alignment of multi-alignment read
		}
		if (!_include_duplicates && (flags & 0x0400) != 0) {
			continue; // A PCR duplicate
		}
		if ((flags & 0x0004) != 0) {
			continue; // Unmapped
		}
		strand_t segment_strand = (flags & 0x0010) != 0 ? neg_strand : pos_strand;

		const auto chrom = _chrom_names.as_chrom(cols[2]);

		// Parse cigar into a parallel array of integers and codes ('S', 'M', 'N', etc)
		int m      = 0; // number of cigar string components (exclude `N`s)
		int n      = 0; // number of splits (Ns) in cigar string
		auto cigar = cols[5];
		GK_CHECK(!empty(cigar), value, "Empty CIGAR string");

		bool allow_indels = cols[9] != "=";

		do {
			GK_CHECK(m < max_cigar_codes, value, "Too many parts in CIGAR string: '{}'", cigar);
			char code  = 0;
			int len   = 0;
			parse_next_cigar(cigar, len, code);
			GK_CHECK(code == 'S' || code == 'M' || code == 'N' || code == 'I' || code == 'D', value,
					 "Unrecognized code '{}' in CIGAR string", code);
			GK_CHECK(allow_indels || (code != 'I' && code != 'D'), value,
					 "Indel '{}' specified, but SEQ is '='", code);

			if (code == 'N') {
				n++;
			}
			codes[m].code   = code;
			codes[m].length = len;
			++m;
		} while (!empty(cigar));

		auto pos = as_pos(cols[3]) - 1;
		strand_t strand = pos_strand;
		auto strand_set = infer_strand(strand, segment_strand, flags);
		strand = infer_strand(strand, strand_set, chrom, pos, codes, m, n);
		process_line(chrom, pos, strand, segment_strand, codes, m, n, cols);
	}
}

bool sam_line_parser::infer_strand(strand_t& out, strand_t segment_strand, unsigned int flags) const
{
	if (_first_read_direction == direction::unknown) {
		return false;
	}

	constexpr unsigned multi_segment    = 0x1;
	constexpr unsigned first_segment    = 0x40;
	constexpr unsigned last_segment     = 0x80;
	constexpr unsigned terminal_segment = 0xC0;
	GK_CHECK((flags & multi_segment) == multi_segment
				 && ((flags & terminal_segment) == first_segment || (flags & terminal_segment) == last_segment),
			 value, "Library format specified paired-end, but single end read encountered");

	// read 1 + forward -> segment_strand
	// read 1 + reverse -> !segment_strand
	// read 2 + forward -> !segment_strand
	// read 2 + reverse -> segment_strand
	auto strand = segment_strand;
	if ((_first_read_direction == direction::reverse && (flags & first_segment) != 0)
		|| (_first_read_direction == direction::forward && (flags & first_segment) == 0)) {
		strand = segment_strand == pos_strand ? neg_strand : pos_strand;
	}
	out = strand;
	return true;
}

strand_t sam_line_parser::infer_strand(strand_t strand, bool strand_set, chrom_t chrom, pos_t pos,
									   const cigar_op* codes, int m, int n) const
{
	if (_introns == nullptr) {
		return strand;
	}
	GK_CHECK(n > 0, value, "Annotation strand inference requires junctions. No N CIGAR op found.");

	for (int i_code = 0, i_n = 0; i_code < m && i_n < n; ++i_code) {
		auto [length, op] = codes[i_code];
		if (op == 'N') {
			auto interval    = interval_t::from_dna0(chrom, pos, pos + length, strand, refg());
			auto overlapping = _introns->find_overlapping(interval);
			if (begin(overlapping) != end(overlapping)) {
				strand_set = true;
			} else {
				GK_CHECK(_first_read_direction == direction::unknown, value,
						 "Annotation strand inference: cannot validate library format");
			}
			if (_first_read_direction == direction::unknown || !strand_set) {
				interval    = interval.as_opp_strand();
				overlapping = _introns->find_overlapping(interval);
				if (begin(overlapping) != end(overlapping)) {
					if (!strand_set) {
						strand = interval.strand;
						strand_set = true;
					} else {
						GK_CHECK(_first_read_direction != direction::unknown, value,
								 "Annotation strand inference: intron can ambiguously align to both strands");
					}
				}
			}
			GK_CHECK(strand_set, value, "Annotation strand inference: cannot align read to an intron");
			++i_n;
		}
		if (op == 'M' || op == 'D' || op == 'N') {
			pos += length;
		}
	}

	return strand;
}

sam_line_parser::sam_line_parser(const genome_t& genome)
: _chrom_names{genome.chrom_names()}
, _interval_filter{[&](interval_t i) {
	GK_CHECK(i.refg == _chrom_names.refg(), value, "Cannot filter {} for {}", i, _chrom_names.refg_name());
}}
, _introns{std::empty(genome.anno()) ? nullptr : &genome.anno().intrs()}
{
}

END_NAMESPACE_GK
