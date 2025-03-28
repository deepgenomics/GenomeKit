/*
Copyright (C) 2016-2023 Deep Genomics Inc. All Rights Reserved.
*/
#include "ralign.h"

#include "format.h"
#include "genome.h"
#include "genome_dna.h"
#include "strutil.h"
#include <cstring>
#include <iterator>
#include <map>
#include <utility>

using namespace std;

BEGIN_NAMESPACE_GK

const unsigned short c_ralign_sig = 0x0efd;
const unsigned short c_ralign_ver = 0x0002;
// versions:
//   0001: initial format
//   0002: Restructured to support all chromosomes and more species
//   <--- INSERT VERSION CHANGE SUMMARIES HERE

/////////////////////////////////////////////////////////////////
//      packed_junction / junction_t / junction_table
/////////////////////////////////////////////////////////////////

junction_t::junction_t(index_t                src, const read_alignments& jralign){ unpack_from(jralign.junctions()[src], jralign); }
junction_t::junction_t(const packed_junction& src, const read_alignments& jralign){ unpack_from(src, jralign); }

void junction_t::unpack_from(const packed_junction& src, const read_alignments& jralign)
{
	this->as_interval() = src.as_interval();
	num_aligns          = src.num_aligns;
	aligns              = reinterpret_cast<const index_t*>(jralign.junctions().aux() + src.aux);
}

/////////////////////////////////////////////////////////////////
//      packed_align / align_t / align_table
/////////////////////////////////////////////////////////////////

align_t::align_t(index_t             src, const read_alignments& jralign){ unpack_from(jralign.alignments()[src], jralign); }
align_t::align_t(const packed_align& src, const read_alignments& jralign){ unpack_from(src, jralign); }

void align_t::unpack_from(const packed_align& src, const read_alignments& jralign)
{
	this->as_interval() = src.as_interval();
	num_matches    = src.num_matches;
	matches.a      = &jralign.matches()[src.match0];
	matches.b      = matches.a + num_matches;
}

/////////////////////////////////////////////////////////////////
//      packed_align_match / align_match_t / align_match_table
/////////////////////////////////////////////////////////////////

align_match_t::align_match_t(index_t                  src, const read_alignments& jralign){ unpack_from(jralign.matches()[src], jralign); }
align_match_t::align_match_t(const packed_align_match& src, const read_alignments& jralign){ unpack_from(src, jralign); }

void align_match_t::unpack_from(const packed_align_match& src, const read_alignments& jralign)
{
	this->as_interval() = src.as_interval();
	num_variants        = src.num_variants;
	variants            = reinterpret_cast<const index_t*>(jralign.matches().aux() + src.aux);
}


////////////////////////////////////////////////////////////////
// read_alignments
////////////////////////////////////////////////////////////////

read_alignments::read_alignments()
: _junctions(*this)
, _alignments(*this)
, _matches(*this)
{
}

const junction_table&    read_alignments::junctions() const
{
	ensure_open();
	return _junctions;
}
const align_table&       read_alignments::alignments() const
{
	ensure_open();
	return _alignments;
}
const align_match_table& read_alignments::matches() const
{
	ensure_open();
	return _matches;
}
const variant_table&     read_alignments::variants() const
{
	ensure_open();
	return _variants;
}

void read_alignments::set_source(string sourcefile)
{
	GK_CHECK2(!is_open(), runtime, "Cannot set source when file already open.");
	_sourcefile = std::move(sourcefile);
}

void read_alignments::open_on_demand() const { const_cast<read_alignments*>(this)->open(); }
void read_alignments::close() { _fmap.close(); }

void read_alignments::open()
{
	GK_CHECK2(!is_open(), runtime, "jraligns_table::open() already opened");
	GK_CHECK2(!_sourcefile.empty(), value, "No file was specified");

	// Memory map the source file
	_fmap.open(_sourcefile);

	// Read the GKJUNCREADS file signature
	unsigned short sig, ver;
	_fmap.read(sig);
	_fmap.read(ver);
	GK_CHECK(sig == c_ralign_sig, file, "Expected valid RALIGN file signature {:x} but found {:x}.", c_ralign_sig, sig);
	GK_CHECK(ver == c_ralign_ver, file, "Expected RALIGN file version {:x} but found {:x}.", c_ralign_ver, ver);

	// Read the tables
	_junctions.load(_fmap);
	_alignments.load(_fmap);
	_matches.load(_fmap);
	_variants.load(_fmap);
	// <-- INSERT NEW TABLES HERE

	_fmap.read_checkpoint(0x85420222);
}

int read_alignments::ralign_version() { return c_ralign_ver; }

read_alignments::builder::builder(const char* outfile, const genome_t& genome)
: sam_line_parser{genome}
, _file{outfile, "w"}
, _dna{&genome.dna()}
{
}

void read_alignments::builder::collect_variant(b_align_match_t& match, const c_variant_t& v)
{
	auto loc = _ra_variants.lower_bound(v);
	if (loc != _ra_variants.end() && !(v < loc->first)) {
		index_t idx = loc->second;
		match.variants.push_back(idx);  // records index to previously added variant
	}
	else {
		auto idx = (index_t)_ra_variants_dup.size();
		_ra_variants.insert(loc, make_pair(v, idx));
		match.variants.push_back(idx);
		_ra_variants_dup.push_back(v);
	}
}

void read_alignments::builder::add(const char* infile)
{
	_verbose = getenv("GENOMEKIT_QUIET") == nullptr;
	if (_verbose)
		print("Loading {}...", infile);

	zline_reader lr(infile);
	try {
		process_file(lr);
	}
	GK_RETHROW2("In SAM file: {}:{}", infile, lr.line_num());

	if (_verbose)
		print("\rLoaded  {}... ({} reads, {} junctions, {} variants)    \n", infile, _num_reads_loaded,
			  _juncs.size(), _ra_variants_dup.size());
}

void read_alignments::builder::process_line(chrom_t chrom, pos_t pos, strand_t strand, strand_t segment_strand,
											const cigar_op* codes, int m, int n, sam_cols cols)
{
	int match_start = pos;
	int match_end   = pos;
	auto seq        = cols[9];
	b_alignment_t  read;
	b_align_match_t a_match;

	bool has_seq     = seq != "*";
	bool has_variant = has_seq && seq != "=";

	if (!is_stranded()) {
		segment_strand = pos_strand;  // use pos_strand to represent unstranded searches
	}

	// Find variants for this read alignment
	for (int i = 0; i < m; ++i) {
		auto [len, code] = codes[i];

		if (has_variant) {
			// Find SNV by comparing SEQ and reference genome
			if (code == 'M') {
				interval_t interval(chrom, pos, pos + len - 1, pos_strand, refg());
				dnastr     ref_seq = (*_dna)(interval);
				GK_ASSERT(ref_seq.size() == len);
				for (int j = 0; j < ref_seq.size(); ++j) {
					if (ref_seq[j] == seq[j]) continue;

					c_variant_t v;
					v.i   = interval_t(chrom, pos + j, pos + j, pos_strand, refg());
					v.ref = ref_seq[j];
					v.alt = seq[j];

					collect_variant(a_match, v);
				}
			}

			// Find Insertion variants
			if (code == 'I') {
				c_variant_t v;
				// end = pos - 1 for empty ref under inclusive-end convention
				v.i   = interval_t(chrom, pos, pos - 1, pos_strand, refg());
				v.ref = "";
				v.alt = seq.substr(0, len);
				collect_variant(a_match, v);
			}

			// Find Deletion variants
			if (code == 'D') {
				interval_t interval(chrom, pos, pos + len - 1, pos_strand, refg());
				dnastr     ref_seq = (*_dna)(interval);
				GK_ASSERT(ref_seq.size() == len);

				c_variant_t v;
				v.i   = interval;
				v.ref.assign(ref_seq.data(), ref_seq.size());
				v.alt = "";
				collect_variant(a_match, v);
			}
		}

		// Record the matching region before N
		if (code == 'N') {
			GK_CHECK2(a_match.variants.size() <= 255, value,
				"Expect number of variants for one read alignment match to be less than 256 but found {}",
				a_match.variants.size());
			// assumption this isn't necessary, but revisit if we remove jralign
			GK_CHECK2(match_start < match_end, value,
						"No alignment match upstream of intron. Consider filtering out the junction-only "
						"reads or use jralign.");

			a_match.i = interval_t::from_dna0(chrom, match_start, match_end, segment_strand, refg());
			read.matches.push_back(std::move(a_match));
			a_match = b_align_match_t{};

			match_start = match_end + len;
		}

		// M,I,S consumes query (i.e. seq); M,D,N consumes reference (i.e. pos)
		if (code == 'M' || code == 'I' || code == 'S') seq.remove_prefix(len);
		if (code == 'M' || code == 'D' || code == 'N') {
			pos += len;
			match_end += len;
		}
	}
	GK_CHECK2(!has_seq || empty(seq), value,
			 "Expect length of CIGAR that consumes query ({}) to correspond to SEQ length ({})",
			 size(cols[9]) - size(seq), size(cols[9]));

	// Remember to collect the last matching region
	GK_CHECK2(a_match.variants.size() <= 255, value,
			"Expect number of variants for one read alignment match to be less than 256 but found {}",
			a_match.variants.size());
	GK_CHECK2(match_start < match_end, value,
				"No ending alignment match. Consider filtering out the junction-only "
				"reads or use jralign.");
	a_match.i = interval_t::from_dna0(chrom, match_start, match_end, segment_strand, refg());
	read.matches.push_back(std::move(a_match));

	// Collects data for this read alignment
	read.i = interval_t::from_dna0(chrom, as_int(cols[3]) - 1, match_end, segment_strand, refg());
	_aligns.push_back(std::move(read));

	// For each split i, contribute a read alignment to the
	// junction corresponding to that split.
	int original_match_start = as_int(cols[3])-1;
	for (int i = 0; i < n; ++i) {

		// Calculate junction start/end by summing up the M, D, and N
		// lengths until we encounter the N that corresponds to junction i.
		int start = original_match_start;
		int j = 0;
		for (int countdown = i+1; j < m; ++j) {
			auto [length, c] = codes[j];
			if (c == 'N')
				if (--countdown == 0)
					break;
			if (c == 'M' || c == 'D' || c == 'N')
				start += length;
		}
		GK_ASSERT(j < m);  // j is now the index of the code, which is 'N', for junction i
		int end = start + codes[j].length;

		// Define the interval and ensure that it satisfies exclude/allow
		auto interval = interval_t::from_dna0(chrom, start, end, strand, refg());
		if (!get_interval_filter().filter(interval))
			continue;

		// Find the junction corresponding to this exact interval
		// (always positive strand for now - TODO: determine likely strand from gene annotations?)
		b_junction_t& junc = _juncs[interval];
		junc.alignments.push_back(int_cast<index_t>(_aligns.size()-1));
	}

	// Count the number of reads loaded, which may be less than the number of alignments created
	_num_reads_loaded++;

	if (_verbose && _num_reads_loaded % 234152 == 0)
		print("\rLoading... ({} reads, {} junctions, {} variants)", _num_reads_loaded, _juncs.size(), _ra_variants_dup.size());
}

void read_alignments::builder::finalize()
{
	GK_ASSERT(_ra_variants.size() == _ra_variants_dup.size());

	junction_table::builder     juncs_builder(is_stranded());
	align_table::builder        aligns_builder(is_stranded());
	align_match_table::builder  matches_builder(is_stranded());
	variant_table::builder      variants_builder(false);

	// Build junction_table
	for (const auto& [interval, src] : _juncs) {

		// NOTE: The reason the junction table is built after-the-fact and not inline with the parsing
		//       is that historically there was a `min_alignments_per_junction` option that would filter
		//       out junctions that did not have enough read support. However, that's not the ideal way
		//       to filter out spurious junctions or alignments, so this feature was removed, and the
		//       loop here in finalize() was retained in case we want to implement smarter filtering
		//       option in the future (since can dramatically, DRAMATICALLY reduce file size if done right).
		//       However, that may be done as part of the RNA-seq pipeline, so not clear if build_ralign
		//       should implement any smart filtering.

		// Fill record for this junction that will be dumped to disk
		packed_junction junction;
		junction.as_interval()  = interval; // Copy the interval fields
		junction.num_aligns     = int_cast<index_t>(src.alignments.size());
		junction.aux            = int_cast<offset_t>(juncs_builder.curr_aux());

		// Fill aux data with record points to.
		juncs_builder.add_elem(junction);
		juncs_builder.add_aux(src.alignments);
	}

	// Build align_table and align_match_table
	for (const auto& src : _aligns) {
		// Fill record for this read alignment that will be dumped to disk
		packed_align alignment;
		alignment.as_interval() = src.i;
		alignment.num_matches   = int_cast<unsigned char>(src.matches.size());
		alignment.match0        = -1;

		// Fills record for an array of alignment matches
		GK_ASSERT(src.matches.size());     // At least 1 alignment, which is the read itself
		for (const auto& match_src : src.matches) {
			packed_align_match match;
			match.as_interval() = match_src.i;
			match.num_variants  = int_cast<unsigned char>(match_src.variants.size());
			match.aux           = int_cast<offset_t>(matches_builder.curr_aux());

			// Fill aux data with record points to.
			index_t idx = matches_builder.add_elem(match);
			matches_builder.add_aux(match_src.variants);

			alignment.match0 = (alignment.match0 == -1) ? idx : alignment.match0;
		}

		// Fill aux data with record points to.
		aligns_builder.add_elem(alignment);
	}

	// Build variant_table
	for (const auto& src : _ra_variants_dup) {
		// Fill record for this read alignment that will be dumped to disk
		packed_variant variant;
		variant.as_interval()   = src.i;
		variant.aux             = int_cast<offset64_t>(variants_builder.curr_aux());

		// Fill aux data with record points to.
		variants_builder.add_elem(variant);
		variants_builder.add_aux(src.ref);
		variants_builder.add_aux(src.alt);
	}

	// Dump signature
	_file.write(c_ralign_sig);
	_file.write(c_ralign_ver);

	// Dump 4 tables
	juncs_builder.dump(_file);
	aligns_builder.dump(_file);
	matches_builder.dump(_file);
	variants_builder.dump(_file);
	// <-- INSERT NEW TABLES HERE

	_file.write_checkpoint(0x85420222);
	_file.close();
}

bool operator<(const read_alignments::builder::c_variant_t& lhs,
               const read_alignments::builder::c_variant_t& rhs) {
	if (lhs.i != rhs.i)     return lhs.i < rhs.i;
	if (auto cmp = lhs.ref.compare(rhs.ref); cmp != 0) return cmp < 0;
	if (auto cmp = lhs.alt.compare(rhs.alt); cmp != 0) return cmp < 0;
	return false;
}

END_NAMESPACE_GK
