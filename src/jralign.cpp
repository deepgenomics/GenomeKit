/*
Copyright (C) 2016-2023 Deep Genomics Inc. All Rights Reserved.
*/
#include "jralign.h"

#include "format.h"
#include "genome.h"
#include "genome_dna.h"
#include "strutil.h"
#include <cstring>
#include <map>
#include <utility>

using namespace std;

BEGIN_NAMESPACE_GK

const unsigned short c_jralign_sig = 0x0dec;
const unsigned short c_jralign_ver = 0x0003;
// versions:
//   0001: initial format
//   0002: include_variants
//   0003: Restructured to support all chromosomes and more species
//   <--- INSERT VERSION CHANGE SUMMARIES HERE

////////////////////////////////////////////////
// jraligns_t
////////////////////////////////////////////////

jraligns_t::jraligns_t(index_t                src, const jraligns_table& table) { unpack_from(table[src], table); }
jraligns_t::jraligns_t(const packed_jraligns& src, const jraligns_table& table) { unpack_from(src,        table); }

INLINE void jraligns_t::unpack_from(const packed_jraligns& src, const jraligns_table& table)
{
	this->as_interval() = src.as_interval();
	_has_variants = src.has_variants;
	_num_reads = src.num_reads;
	_overhangs = (const unsigned char*)(table.aux() + src.aux);
	_strands   = (const unsigned char*)(_overhangs + 2*num_reads()); // 2*reads is OK since clamped to 2^31
	_indices_begin    = has_variants() ? (const unsigned char*)(_strands + (num_reads() + 7) / 8) : nullptr;
	_variants_indices = has_variants() ? (const index_t*)(_indices_begin + (num_reads() + 1)) : nullptr;
}

jralign_t jraligns_t::operator[](unsigned i) const
{
	GK_CHECK(i < num_reads(), index, "Invalid read index {}", i);

	// 2*i+1 is OK since clamped to 2^31-1
	return jralign_t(_overhangs[2*i + 0],
					 _overhangs[2*i + 1],
					 _strands[i / 8] & (1 << i%8) ? pos_strand : neg_strand,
					 _indices_begin ? (
						 (_indices_begin[i+1] >= _indices_begin[i]) ?
						 (_indices_begin[i+1] - _indices_begin[i]) : (0xff - _indices_begin[i] + _indices_begin[i+1])
					) : 0,
					_indices_begin ? (_variants_indices + _indices_begin[i]) : nullptr);
}

////////////////////////////////////////////////////////////////
// junction_read_alignments
////////////////////////////////////////////////////////////////

junction_read_alignments::junction_read_alignments()
: _juncs(*this)
{
}

const jraligns_table& junction_read_alignments::juncs() const
{
	ensure_open();
	return _juncs;
}
const variant_table& junction_read_alignments::variants() const
{
	ensure_open();
	return _variants;
}

void junction_read_alignments::set_source(string sourcefile)
{
	GK_CHECK(!is_open(), runtime, "Cannot set source when file already open.");
	_sourcefile = std::move(sourcefile);
}

void junction_read_alignments::open_on_demand() const { const_cast<junction_read_alignments*>(this)->open(); }
void junction_read_alignments::close() { _fmap.close(); }

void junction_read_alignments::open()
{
	GK_CHECK(!is_open(), runtime, "jraligns_table::open() already opened");
	GK_CHECK(!_sourcefile.empty(), value, "No file was specified");

	// Memory map the source file
	_fmap.open(_sourcefile);

	// Read the GKJUNCREADS file signature
	unsigned short sig, ver;
	_fmap.read(sig);
	_fmap.read(ver);
	GK_CHECK(sig == c_jralign_sig, file, "Expected valid JRALIGN file signature {:x} but found {:x}.", c_jralign_sig, sig);
	GK_CHECK(ver == c_jralign_ver, file, "Expected JRALIGN file version {:x} but found {:x}.", c_jralign_ver, ver);

	// Read number of total reads
	_fmap.read(_num_reads);

	// Read the junction table itself
	_juncs.load(_fmap);

	// Read the variant table optionally
	if (_fmap.read<unsigned>() != 0)
		_variants.load(_fmap);

	// <-- INSERT NEW TABLES HERE
	_fmap.read_checkpoint(0x85420333);
}

int junction_read_alignments::jralign_version() { return c_jralign_ver; }

junction_read_alignments::builder::builder(const char* outfile, const genome_t& genome)
: sam_line_parser{genome}
, _file{outfile, "w"}
, _dna(&genome.dna())
{
}

void junction_read_alignments::builder::add(const char* infile)
{
	_verbose = getenv("GENOMEKIT_QUIET") == nullptr;
	if (_verbose)
		print("Loading {}...", infile);

	zline_reader lr(infile);
	try {
		process_file(lr);
	}
	GK_RETHROW("In SAM file: {}:{}", infile, lr.line_num());

	if (_verbose)
		print("\rLoaded  {}... ({} reads, {} junctions)    \n", infile, _num_reads_loaded,
			  _juncs.size());
}

void junction_read_alignments::builder::collect_variant(const b_variant_t& v)
{
	auto loc = _ra_variants.lower_bound(v);
	if (loc != _ra_variants.end() && !(v < loc->first))
		_ra_variants_indices.push_back(loc->second);
	else {
		auto index = (index_t)_ra_variants_dup.size();
		_ra_variants_dup.push_back(v);
		_ra_variants.insert(loc, make_pair(v, index));
		_ra_variants_indices.push_back(index);
	}
}

void junction_read_alignments::builder::process_line(chrom_t chrom, pos_t pos, strand_t strand, strand_t segment_strand,
													 const cigar_op* codes, int m, int n, sam_cols cols)
{
	// Find variants for this read
	auto seq = cols[9];
	bool has_variant = seq != "*" && seq != "=";
	if (include_variants() && has_variant) {
		_ra_variants_indices.clear();

		for (int i = 0; i < m; ++i) {
			auto [len, code] = codes[i];

			// Find SNV by comparing SEQ and reference genome
			// specifically, dnastr for each DNA sequence in alignment matches (M)
			if (code == 'M') {
				interval_t interval(chrom, pos, pos + len - 1, pos_strand, refg());
				dnastr ref_seq = (*_dna)(interval);
				GK_ASSERT(ref_seq.size() == len);
				for (int j = 0; j < ref_seq.size(); ++j) {
					if (ref_seq[j] == seq[j])
						continue;

					b_variant_t v;
					v.i   = interval_t(chrom, pos + j, pos + j, pos_strand, refg());
					v.ref = ref_seq[j];
					v.alt = seq[j];

					collect_variant(v);
				}
			}
			// Find Insertion variants
			if (code == 'I') {
				b_variant_t v;
				// end = pos - 1 for empty ref under inclusive-end convention
				v.i   = interval_t(chrom, pos, pos - 1, pos_strand, refg());
				v.ref = "";
				v.alt = seq.substr(0, len);
				collect_variant(v);
			}

			// Find Deletion variants
			if (code == 'D') {
				interval_t interval(chrom, pos, pos + len - 1, pos_strand, refg());
				dnastr ref_seq = (*_dna)(interval);
				GK_ASSERT(ref_seq.size() == len);

				b_variant_t v;
				v.i   = interval;
				v.ref.assign(ref_seq.data(), ref_seq.size());
				v.alt = "";
				collect_variant(v);
			}

			// M,I,S consumes query (i.e. seq); M,D,N consumes reference (i.e. pos)
			if (code == 'M' || code == 'I' || code == 'S')
				seq.remove_prefix(len);
			if (code == 'M' || code == 'D' || code == 'N')
				pos += len;
		}
		GK_CHECK(empty(seq), value,
				 "Expect length of CIGAR that consumes query ({}) to correspond to SEQ length ({})",
				 size(cols[9]) - size(seq), size(cols[9]));
		GK_CHECK(_ra_variants_indices.size() <= 255, value,
				 "Expect number of variants for one read to be less than 256 but found {}",
				 _ra_variants_indices.size());
	}

	// For each split i, contribute a read alignment to the
	// junction corresponding to that split.
	int match_start = as_int(cols[3]) - 1;
	for (int i = 0; i < n; ++i) {
		// Calculate junction start/end by summing up the M, D, and N
		// lengths until we encounter the N that corresponds to junction i.
		int start = match_start;
		int j     = 0;
		for (int countdown = i + 1; j < m; ++j) {
			auto [length, c] = codes[j];
			if (c == 'N')
				if (--countdown == 0)
					break;
			if (c == 'M' || c == 'D' || c == 'N')
				start += length;
		}
		GK_ASSERT(j < m); // j is now the index of the code, which is 'N', for junction i
		int end = start + codes[j].length;

		// Calculate left/right overhangs by summing all S, M, and I lengths.
		// before and after the junction's code entry.
		int left = 0;
		for (int k = 0; k < j; ++k) {
			auto [length, c] = codes[k];
			if (c == 'M')
				left += length;
		}
		int right = 0;
		for (int k = j + 1; k < m; ++k) {
			auto [length, c] = codes[k];
			if (c == 'M')
				right += length;
		}
		if (_overhang_error == error_handling::error) {
			GK_CHECK(left <= 255 && right <= 255, value,
						"Overhangs >255 not supported for overhangs ({}, {})", left, right);
		} else if (_overhang_error == error_handling::clamp) {
			left = min(left, 255);
			right = min(right, 255);
		}

		// If a match on either side of the intron is less than min_overhang, drop it.
		if (left < _min_overhang || right < _min_overhang)
			continue;

		// Define the interval and ensure that it satisfies exclude/allow
		auto interval = interval_t::from_dna0(chrom, start, end, strand, refg());
		if (!get_interval_filter().filter(interval))
			continue;

		// Find the junction corresponding to this exact interval
		jralign_entry& junc = _juncs[interval];
		size_t n            = junc.overhangs.size() / 2;
		junc.overhangs.push_back(int_cast<unsigned char>(left));
		junc.overhangs.push_back(int_cast<unsigned char>(right));
		if (n % 8 == 0)
			junc.strands.push_back(0);
		if (segment_strand == pos_strand)
			junc.strands[n / 8] |= 1 << (n % 8); // orientation

		if (include_variants()) {
			// Collects variant into this junction entry
			if (!empty(_ra_variants_indices)) {
				junc.variants_indices.insert(junc.variants_indices.end(), _ra_variants_indices.begin(),
											 _ra_variants_indices.end());
			}

			// Keep track of the total number of variants collected by this junction so far
			// This value act as index into which subsequent variants are being collected
			// Note raligns_variants_indices.size() < 256 checked at runtime
			junc.indices_begin.push_back(junc.indices_begin.back()
										 + static_cast<unsigned char>(_ra_variants_indices.size()));
		}
	}

	// Count the number of reads loaded, which may be less than the number of alignments created
	_num_reads_loaded++;

	if (_verbose && _num_reads_loaded % 234152 == 0)
		print("\rLoading... ({} reads, {} junctions)", _num_reads_loaded,
			  _juncs.size());
}

void junction_read_alignments::builder::finalize()
{
	long long num_reads_final = 0;
	GK_ASSERT(_ra_variants.size() == _ra_variants_dup.size());

	jraligns_table::builder juncs_builder(is_stranded());
	variant_table::builder  variants_builder(false);

	// Build jraligns_table
	for (const auto& [interval, src] : _juncs) {
		// Ignore junctions with too few reads across them
		auto num_reads = int_cast<unsigned>(src.overhangs.size() / 2);
		if (num_reads < _min_reads)
			continue;

		// Sets most significant bits of num_reads to indicate
		// any reads of this junction has variants
		GK_CHECK(!(num_reads & (1 << 31)), value, "Exceeded 2^31 reads limit on a single junction");

		// Fill record for this junction that will be dumped to disk
		packed_jraligns junc;
		junc.as_interval() = interval; // Copy the interval fields
		junc.has_variants = src.variants_indices.empty() ? 0 : 1;
		junc.num_reads = num_reads;
		junc.aux = int_cast<offset_t>(juncs_builder.curr_aux());
		juncs_builder.add_elem(junc);  // Add it to the table

		// Fill the aux data that the record points to.
		juncs_builder.add_aux(src.overhangs);
		juncs_builder.add_aux(src.strands);
		if (!src.variants_indices.empty()) {
			GK_ASSERT(include_variants());
			juncs_builder.add_aux(src.indices_begin);
			juncs_builder.add_aux(src.variants_indices);
		}

		num_reads_final += num_reads;
	}
	// Build variant_table
	for (const auto& built_variant : _ra_variants_dup) {
		packed_variant v;
		v.as_interval() = built_variant.i;
		v.aux = int_cast<offset64_t>(variants_builder.curr_aux());

		variants_builder.add_elem(v);
		variants_builder.add_aux(built_variant.ref);
		variants_builder.add_aux(built_variant.alt);
	}

	// Dump signature
	_file.write(c_jralign_sig);
	_file.write(c_jralign_ver);

	// Dump the total number of reads in the file
	_file.write(num_reads_final);

	// Dump the read alignment tables
	juncs_builder.dump(_file);

	// Optionally dump the variant table
	if (include_variants()) {
		_file.write<unsigned>(1);
		variants_builder.dump(_file);
	} else {
		_file.write<unsigned>(0);
	}
	// <-- INSERT NEW TABLES HERE

	_file.write_checkpoint(0x85420333);
	_file.close();
}

bool operator<(const junction_read_alignments::builder::b_variant_t& lhs, const junction_read_alignments::builder::b_variant_t& rhs) {
	if (lhs.i != rhs.i)     return lhs.i < rhs.i;
	if (auto cmp = lhs.ref.compare(rhs.ref); cmp != 0) return cmp < 0;
	if (auto cmp = lhs.alt.compare(rhs.alt); cmp != 0) return cmp < 0;
	return false;
}

END_NAMESPACE_GK
