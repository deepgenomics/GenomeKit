/*
Copyright (C) 2016-2023 Deep Genomics Inc. All Rights Reserved.
*/
#include "genome_anno.h"

#include "genome.h"
#include "strutil.h"
#include <algorithm>
#include <filesystem>
#include <fstream>
#include <iterator>
#include <map>
#include <numeric>
#include <optional>
#include <set>
#include <string_view>
#include <unordered_map>
#include <utility>

using namespace std;

BEGIN_NAMESPACE_GK

constexpr unsigned short c_dganno_sig = 0xedba;
// NOTE: If you bump this, then `python -m genome_kit build --anno --test-anno`
//       Also `git rm` the older versioned test annos
constexpr unsigned short c_dganno_ver = 0x0016;
// versions:
//   0001: initial format
//   0002: add max interval size to interval_table indices
//   0003: add intron table
//   0004: add cds table
//   0005: only store index of first child element, since indices are contiguous
//   0006: add is_first/is_last to tran/exon/intr/cds; add cds list to transcripts
//   0007: add "translated_processed_pseudogene" to biotypes
//   0008: add per-transcript status field to packed_tran ("transcript_status")
//   0009: tables write 64-bit sizes to support aux pool > 2GB
//   000a: pos5 index secondarily sorted by pos3
//   000b: add version string
//         add attributes:
//            tran.product     -- "breast cancer type 1 susceptibility protein isoform 2"
//            tran.ccds_id     -- "CCDS53759.1"
//            tran.protein_id  -- "ENSP00000350283.3"
//            exon.id          -- "ENSE00002281982.1"
//         remove attributes:
//            tran.is_first (unused)
//            tran.is_last  (unused)
//   000c: more biotypes from Ensembl
//   000d: UCSC refseq rebuilt from original 2017-06-25 dated sources
//   000e: NCBI refseq off-by-1 end issue
//   000f: use new checkpoint value for stranded tables
//   0010: NCBI refseq v109 uses type==pseudogene for some genes; more biotypes
//   0011: biotypes and tsl for non-gk standard Ensembl annotations
//   0012: NCBI use refseq sources (which are archived) instead of all(GCF), but some indices are swapped
//   0013: biotypes from Gencode.vM31
//   0014: biotypes from Gencode.v41
//   0015: support for utr
//   0016: Restructured to support all chromosomes and more species
//   <--- INSERT VERSION CHANGE SUMMARIES HERE

namespace {
constexpr elem_num          transcript_exon_index_max   = 32767;
// besides cds, these maximums also apply to utr5 and utr3
constexpr elem_num          cds_exon_index_max          = 4094;
constexpr elem_num          cds_exon_index_invalid      = 4095;
constexpr const char* const cds_exon_index_invalid_cstr = "4095";

// The Ensembl v35+ lift37 pipeline had a bug where 0-based indexing was used in the exon_number
// attribute instead of 1-based for the original h38 file.
constexpr auto gencode_v35_lift37_0based_bug = true;
}  // namespace

// swap attr<->field order to avoid accidentally using wrong get_attr
static string_view get_attr_alt(string_view field, const vector<string_view>& attrs, const char* alternate_field)
{
	auto value = get_attr(attrs, field, "");
	if (std::empty(value))
		value = get_attr(attrs, alternate_field);
	return value;
}

/////////////////////////////////////////////////////////////

// Convert between exon number and string

static elem_num as_elem_num(string_view s)
{
	return int_cast<elem_num>(as_int(s) - 1);
}

static elem_num as_elem_num(string_view s, elem_num max_value)
{
	const elem_num num = as_elem_num(s);
	GK_CHECK(num <= max_value, runtime,
			 "Overflow detected when parsing elem_num: {} is greater than {}.", s, max_value);
	return num;
}

/////////////////////////////////////////////////////////////

// Convert between evidence level and string

// If you change the list, the indexes associated with types will change.
constexpr const char* g_elevel_to_cstr[num_elevel] = {
	"1",
	"2",
	"3",
};

static elevel_t as_elevel(string_view s)
{
	auto start = cbegin(g_elevel_to_cstr);
	auto stop  = cend(g_elevel_to_cstr);
	auto it    = find(start, stop, s);
	GK_CHECK(it != stop, value, "Unrecognized evidence level '{}'.", s);
	return elevel_t(distance(start, it));
}

static bool iequal(string_view x, string_view y)
{
	return equal(cbegin(x), cend(x), cbegin(y), cend(y), [](auto u, auto v) { return lower(u) == lower(v); });
}

static elevel_t ensembl_source_as_elevel(string_view data_source)
{
	// if not gencode, use the ensembl data source to map to confidence level as per FAQ
	// https://www.gencodegenes.org/faq.html
	if (iequal("ensembl", data_source)) {
		return elevel_t{ 2 };
	} else if (iequal("havana", data_source) || iequal("ensembl_havana", data_source)) {
		return elevel_t{ 1 };
	}
	return invalid_elevel;
}

static elevel_t ncbi_source_as_elevel(string_view data_source)
{
	// for refseq, follow https://www.ncbi.nlm.nih.gov/genome/annotation_euk/process/#products
	if (iequal("Gnomon", data_source)) {
		return elevel_t{ 2 };
	} else if (iequal("Curated Genomic", data_source) || iequal("Curated Genomic%2CGnomon", data_source)
			   || iequal("RefSeq", data_source) || iequal("BestRefSeq", data_source)
			   || iequal("BestRefSeq%2CGnomon", data_source)) {
		return elevel_t{ 1 };
	}
	return invalid_elevel;
}

static elevel_t as_elevel(string_view field, const vector<string_view>& attrs, string_view data_source)
{
	auto level   = get_attr(attrs, field, "");
	elevel_t ret = invalid_elevel;
	if (!std::empty(level))
		ret = as_elevel(level);
	else {
		if (ret = ensembl_source_as_elevel(data_source); ret == invalid_elevel)
			ret = ncbi_source_as_elevel(data_source);
	}

	return ret;
}

/////////////////////////////////////////////////////////////

// Convert between tsl and string

constexpr const char* g_tsl_as_cstr[num_tsl] = { "1", "2", "3", "4", "5", "NA" };

static tsl_t as_tsl(string_view s)
{
	auto start = cbegin(g_tsl_as_cstr);
	auto stop  = cend(g_tsl_as_cstr);
	auto it    = find(start, stop, s.substr(0, s.find(" "))); // handle "1 (assigned to previous version x)"
	GK_CHECK(it != stop, value, "Unrecognized transcript_support_level '{}'.", s);
	return tsl_t(distance(start, it));
}

/////////////////////////////////////////////////////////////

// Convert between biotype and string

// This list must stay in sorted (case-sensitive) order for the sake of binary search.
// If you change the list, the indexes associated with types will change.
// validate via http://sequenceontology.org/browser/obob.cgi
constexpr const char* g_biotype_to_cstr[num_biotype] = {
	"3prime_overlapping_ncRNA",
	"C_region",
	"C_region_pseudogene",
	"D_segment",
	"D_segment_pseudogene",
	"IG_C_gene",
	"IG_C_pseudogene",
	"IG_D_gene",
	"IG_D_pseudogene",
	"IG_J_gene",
	"IG_J_pseudogene",
	"IG_LV_gene",
	"IG_V_gene",
	"IG_V_pseudogene",
	"IG_pseudogene",
	"J_segment",
	"J_segment_pseudogene",
	"Mt_rRNA",
	"Mt_tRNA",
	"RNase_MRP_RNA",
	"RNase_P_RNA",
	"SRP_RNA",
	"TEC",
	"TR_C_gene",
	"TR_D_gene",
	"TR_J_gene",
	"TR_J_pseudogene",
	"TR_V_gene",
	"TR_V_pseudogene",
	"V_segment",
	"V_segment_pseudogene",
	"Y_RNA",
	"antisense",
	"antisense_RNA",
	"artifact",
	"bidirectional_promoter_lncRNA",
	"guide_RNA",
	"lincRNA",
	"lnc_RNA",
	"macro_lncRNA",
	"miRNA",
	"misc_RNA",
	"ncRNA",
	"ncRNA_pseudogene",
	"non_coding",
	"non_stop_decay",
	"nonsense_mediated_decay",
	"piRNA",
	"polymorphic_pseudogene",
	"pre_miRNA",
	"processed_pseudogene",
	"processed_transcript",
	"protein_coding",
	"protein_coding_CDS_not_defined",
	"protein_coding_LoF",
	"pseudogene",
	"rRNA",
	"rRNA_pseudogene",
	"retained_intron",
	"ribozyme",
	"sRNA",
	"scRNA",
	"scaRNA",
	"sense_intronic",
	"sense_overlapping",
	"snRNA",
	"snoRNA",
	"tRNA",
	"tRNA_pseudogene",
	"telomerase_RNA",
	"transcribed_processed_pseudogene",
	"transcribed_unitary_pseudogene",
	"transcribed_unprocessed_pseudogene",
	"translated_processed_pseudogene",
	"translated_unprocessed_pseudogene",
	"unitary_pseudogene",
	"unprocessed_pseudogene",
	"vaultRNA",
};

static biotype_t as_biotype(string_view s)
{
	auto start = cbegin(g_biotype_to_cstr);
	auto stop  = cend(g_biotype_to_cstr);
	struct CheckSorted {
		CheckSorted(decltype(start) start, decltype(stop) stop)
		{
			GK_ASSERT(is_sorted(start, stop, [](auto lhs, auto rhs) { return strcmp(lhs, rhs) < 0; }));
		}
	};
	[[maybe_unused]] static const CheckSorted check_sorted(start, stop);

	auto r = lower_bound(start, stop, s);
	GK_CHECK((r != stop) && s == *r, value, "Unrecognized biotype string '{}'.", s);
	return biotype_t(distance(start, r));
}

const char* biotype_as_cstr(biotype_t type)
{
	GK_DBASSERT(type < num_biotype);
	return g_biotype_to_cstr[type];
}

static biotype_t as_biotype_compatible(string_view s)
{
	static const unordered_map<string_view, string_view> biotype_aliases{
		{ "3prime_overlapping_ncrna", "3prime_overlapping_ncRNA" },
		{ "lincrna", "lincRNA" },
		{ "lncRNA", "lnc_RNA" },
		{ "other", "misc_RNA" },
		{ "transcribed_pseudogene", "transcribed_unprocessed_pseudogene" },
		{ "unknown_likely_coding", "misc_RNA" },
		{ "vault_RNA", "vaultRNA" },
	};

	if (auto it = biotype_aliases.find(s); it != cend(biotype_aliases))
		return as_biotype(it->second);
	return as_biotype(s);
}

/////////////////////////////////////////////////////////////

using packed_exons_t = vector<packed_exon>;
using packed_cdss_t  = vector<packed_cds>;
using packed_utrs_t  = vector<packed_utr>;

struct gff_record {
	explicit gff_record(string_view line, const chrom_names_t* chrom_names) : chrom_names(chrom_names)
	{
		using namespace std;
		GK_CHECK(split_view(line, '\t', cols, (int)size(cols)) == (int)size(cols), file,
				 "Expected {} tab-separated columns.", size(cols));
	}

	chrom_t        seqid() const noexcept { return chrom_names->as_chrom(cols[0]); }
	string_view    source() const noexcept { return cols[1]; }
	string_view    type() const noexcept { return cols[2]; }
	uint8_t        phase() const noexcept { return int_cast<uint8_t>(as_int(cols[7])); }
	void           set_attributes(vector<string_view>& out) const noexcept { split_view(cols[8], ';', out); }

	interval_t as_interval() const noexcept
	{
		return interval_t::from_dna0(seqid(), as_pos(cols[3]) - 1, as_pos(cols[4]), as_strand(cols[6]), chrom_names->refg());
	}

	string_view cols[9];
	const chrom_names_t* chrom_names;
};

class genome_anno::builder {
public:
	builder(const genome_t& genome): chrom_names{genome.chrom_names()}
	{
		verbose = getenv("GENOMEKIT_QUIET") == nullptr;
	}

	refg_t refg() const { return chrom_names.refg(); }

	void load_gencode(const char* infile)
	{
		if (verbose)
			print("Loading {}...\n", infile);

		// First figure out the reference genome by looking for aliases in the header lines
		//#description: evidence-based annotation of the human genome, version 25 (Ensembl 85), mapped to GRCh37 with gencode-backmap
		zline_reader lr(infile);

		const int num_cols = 9;
		string_view cols[num_cols];
		vector<string_view> attrs;
		map<string, offset_t, std::less<>> exon_id_offset;
		// handle unordered rows (strand-reversed or inconsistent CDS-interleaving within ensembl) in second pass
		packed_exons_t exons_to_link;
		packed_cdss_t cdss_to_link;
		packed_utrs_t utr5s_to_link;
		packed_utrs_t utr3s_to_link;
		packed_utrs_t generic_utrs_to_link;
		string last_tran_id;
		string last_gene_id;
		set<string, less<>> unknown_gff_types;

		static constexpr string_view transcript_prefix = "transcript:";

		for (; !lr.done(); ++lr) {
			auto line = lr.line();

			// Skip empty lines, comment lines, any line not starting with 'c' for "chr" or chromosome number
			if (line[0] == '#' || line[0] == 0)
				continue;

			// Split tab-separated level. The fields at this level are
			// [0] chr
			// [1] data source
			// [2] element type
			// [3] start position (1-based index, inclusive)
			// [4] end position (1-based index, inclusive)
			// [5] score (if applicable)
			// [6] strand
			// [7] phase (if applicable)
			// [8] attributes (a semicolon-separated string)
			try {
				GK_CHECK(split_view(line, '\t', cols, num_cols) == num_cols, value,
						 "Expected at least {} tab-separated columns", num_cols);

				// Pull out fields at tab-separated level
				const auto type = cols[2];

				const auto parsed_interval = [&]() {
					// convert from 1-based closed to 0-based half-open
					const auto chrom = chrom_names.as_chrom(cols[0]);
					return interval_t::from_dna0(chrom, as_pos(cols[3]) - 1, as_pos(cols[4]), as_strand(cols[6]), refg());
				};

				static constexpr string_view gene_prefix       = "gene:";

				if (type == "gene" || startswith(cols[8], "ID=gene:")) { // ensembl occasionally misnames type
					split_view(cols[8], ';', attrs);

					auto id   = get_attr(attrs, "ID");
					auto name = get_attr(attrs, "gene_name", "");
					auto type = get_attr_alt("gene_type", attrs, "biotype");

					id = skip_prefix(id, gene_prefix);
					if (std::empty(name))
						name = get_attr_alt("Name", attrs, "gene_id");

					// Fill the gene entry fields
					packed_gene gene = {parsed_interval()};
					// Gene-specific fields
					gene.level  = as_elevel("level", attrs, cols[1]);
					gene.type   = as_biotype_compatible(type);
					gene.num_trans = 0;            // Start as 0. Incremented as transcripts are read in. Child transcripts always defined after their parent gene.
					gene.tran0 = invalid_index;
					// Gene-specific variable-length fields to the aux able
					gene.aux = (offset_t)genes.curr_aux();
					genes.add_aux(id);
					genes.add_aux(name);

					// Add the entry and allow transcripts to look up its index by gene ID
					genes.add_elem(gene);
					last_gene_id = id;

				} else if (type == "transcript" || startswith(cols[8], "ID=transcript:")) {
					// process any unlinked children from previous transcript
					link(exons_to_link, cdss_to_link, utr5s_to_link, utr3s_to_link, generic_utrs_to_link);

					split_view(cols[8], ';', attrs);

					auto gene_id = get_attr(attrs, "Parent");
					auto id      = get_attr(attrs, "ID");
					auto type    = get_attr_alt("transcript_type", attrs, "biotype");
					auto tsl     = get_attr(attrs, "transcript_support_level", "NA");
					auto ccds_id = get_attr(attrs, "ccdsid", "");
					auto protein_id = get_attr(attrs, "protein_id", "");

					gene_id = skip_prefix(gene_id, gene_prefix);
					id      = skip_prefix(id, transcript_prefix);

					// If this fails: use a more heavy-weight map recording all prior genes
					GK_CHECK(!genes.empty() && last_gene_id == gene_id, key,
							 "Preceding gene '{}' is not the parent of transcript '{}'", last_gene_id,
							 id);
					index_t gene_index = genes.size() - 1;

					// Fill the transcript entry fields
					packed_tran tran = {parsed_interval()};
					// Transcript-specific fields
					tran.level  = as_elevel("level", attrs, cols[1]);
					tran.tsl    = as_tsl(tsl);
					tran.type   = as_biotype_compatible(type);
					tran.gene = gene_index;
					// Transcript-specific variable-sized fields in aux table
					tran.aux = (offset_t)trans.curr_aux();
					trans.add_aux(id);
					trans.add_aux(ccds_id);
					trans.add_aux(protein_id);
					trans.add_aux("");   // No 'product' field in GENCODE. Empty string ensures None is returned.

					// Add the entry to the transcript table
					index_t index = trans.add_elem(tran);
					last_tran_id  = id;

					// Increment the parent gene's transcript counter
					packed_gene& parent = genes[tran.gene];
					parent.num_trans++;
					if (parent.tran0 == invalid_index)
						parent.tran0 = index;
					GK_ASSERT(index == parent.tran0 + parent.num_trans-1);

				} else if (type == "exon") {
					split_view(cols[8], ';', attrs);
					auto interval = parsed_interval();

					auto tran_id  = get_attr(attrs, "Parent");
					auto num      = get_attr(attrs, "exon_number", "");
					auto id       = get_attr(attrs, "exon_id");
					if (std::empty(num)) num = get_attr(attrs, "rank");

					tran_id = skip_prefix(tran_id, transcript_prefix);

					// If this fails: use a more heavy-weight map recording all prior transcripts
					GK_CHECK(!trans.empty() && last_tran_id == tran_id, key,
							 "Preceding transcript '{}' is not the parent of exon '{}'", last_tran_id,
							 id);
					index_t tran_index = trans.size() - 1;

					// Fill the exon entry fields
					packed_exon exon = {interval};
					// Exon-specific fields
					exon.is_last = false;
					exon.index  = as_elem_num(num, transcript_exon_index_max);
					exon.tran = tran_index;
					exon.cds  = invalid_index; // May be set to a CDS index later
					exon.utr5  = invalid_index;
					exon.utr3  = invalid_index;
					// Exon-specific variable-sized fields in aux table
					auto id_offset = exon_id_offset.lower_bound(id);
					if (id_offset == exon_id_offset.end() || id_offset->first != id) {
						// This exon ID string hasn't appeared before, so add it to the pool
						exon.id = (offset_t)exons.curr_aux();
						exons.add_aux(id);
						exon_id_offset.insert(id_offset, pair<string, offset_t>(id, exon.id));
					} else {
						// This exon ID string has appeared before, so re-use it from the existing string pool
						exon.id = id_offset->second;
					}
					exons_to_link.push_back(exon);
				} else if (type == "CDS") {
					auto cds = process_exonic_component<packed_cds>(cols, last_tran_id, parsed_interval);
					phase_t phase = (pos_t)as_int(cols[7]);
					cds.phase    = (unsigned short)phase;
					cdss_to_link.push_back(cds);
				} else if (type == "five_prime_UTR") {
					const auto utr5 = process_exonic_component<packed_utr>(cols, last_tran_id, parsed_interval);
					utr5s_to_link.push_back(utr5);
				} else if (type == "three_prime_UTR") {
					const auto utr3 = process_exonic_component<packed_utr>(cols, last_tran_id, parsed_interval);
					utr3s_to_link.push_back(utr3);
				} else if (type == "UTR") {
					auto throwaway = process_exonic_component<packed_utr>(cols, last_tran_id, parsed_interval);
					generic_utrs_to_link.push_back(throwaway);
				} else if (auto loc = unknown_gff_types.lower_bound(type);
						   loc == cend(unknown_gff_types) || *loc != type) {
					unknown_gff_types.emplace_hint(loc, type);
					print("Unknown GFF type={} on line {}.\n", type, lr.line_num());
				}
				// <--- INSERT NEW ANNOTATION TYPES HERE
			}
			GK_RETHROW("In GFF3 file: {}:{}", infile, lr.line_num());

			if (verbose && lr.line_num() % 8000 == 1)
				print("  {} genes, {} transcripts, {} exons...\r", genes.size(), trans.size(), exons.size() + exons_to_link.size());

			//if (lr.line_num() > 10000) break; // UNCOMMENT FOR SPEED
		}
		link(exons_to_link, cdss_to_link, utr5s_to_link, utr3s_to_link, generic_utrs_to_link); // process any unlinked children from final transcript

		if (verbose)
			print("  {} genes, {} transcripts, {} exons...\n", genes.size(), trans.size(), exons.size());
	}

	void link(const packed_exon& exon)
	{
		// Add the entry to the exon table
		index_t index = exons.add_elem(exon);

		// Make note of first exon's index in the parent struct
		packed_tran& parent = trans[exon.tran];
		parent.num_exons++;
		GK_CHECK(exon.index <= parent.num_exons, index, "Exon number {} is too large for exon {}.", exon.index, exon);
		if (parent.exon0 == invalid_index)
			parent.exon0 = index;
		GK_ASSERT(index == parent.exon0 + parent.num_exons-1);

		// If this isn't the first exon, create a corresponding upstream intron record
		if (exon.index > 0) {
			// Find 3p end of previous exon and 5p end of current exon
			pos_t prev_pos3 = exons[parent.exon0 + exon.index - 1].pos3;
			pos_t curr_pos5 = exon.pos5;

			// 1-based closed to 0-based closed
			packed_intr intr = {exon.with_pos(prev_pos3, curr_pos5).expand(-1, -1)};
			GK_ASSERT(intr.size() >= 0);
			// Intron-specific fields
			intr.is_last = false;
			intr.index  = exon.index-1;  // intron 0 precedes exon 1
			intr.tran = exon.tran;

			// Add the entry to the intron table
			index_t index = intrs.add_elem(intr);

			// Make note of first intron's index in the parent struct
			if (parent.intr0 == invalid_index)
				parent.intr0 = index;
			GK_ASSERT(index == parent.intr0 + parent.num_exons-2);
		}
	}

	template <class T> void link_exonic_element(
			typename genome_anno_table<T>::builder &table,
			const T& element,
			index_t packed_exon::*at_index_ptr,
			elem_num packed_tran::*el_num_ptr,
			index_t packed_tran::*start_ptr,
			const char* type_name) {
		// Add the entry
		index_t index = table.add_elem(element);
		// Make note of element's index in the parent exon and parent transcript
		packed_tran& parent_tran = trans[element.tran];
		packed_exon& parent_exon = exons[parent_tran.exon0 + element.exon_index];
		auto el_idx = parent_exon.*at_index_ptr;
		GK_CHECK(el_idx == invalid_index, runtime, "Expected unique {0} per exon for {0} {1}.",
				 type_name, element);
		parent_exon.*at_index_ptr = index;
		(parent_tran.*el_num_ptr)++;

		if ((parent_tran.*start_ptr) == invalid_index)
			(parent_tran.*start_ptr) = index;
		GK_ASSERT(index == (parent_tran.*start_ptr) + (parent_tran.*el_num_ptr) - 1);
	}

	template <class T, class OrderF, class It, class PredicateF>
	void link_exonic_elements(
			vector<T>& elements_to_link,
			OrderF& upstream_ordering,
			It& exon_start,
			It& exon_stop,
			PredicateF &needs_remap,
			typename genome_anno_table<T>::builder &table,
			index_t packed_exon::*at_index_ptr,
			elem_num packed_tran::*el_num_ptr,
			index_t packed_tran::*start_ptr,
			const char* type_name) {
		if (elements_to_link.empty()) {
			return;
		}
		sort(begin(elements_to_link), end(elements_to_link), upstream_ordering);

		auto exon_it = exon_start;
		for_each(begin(elements_to_link), end(elements_to_link), [&](auto& element) {
			if (needs_remap(type_name, element)) {
				exon_it = find_if(exon_it, exon_stop, [&](auto exon) { return exon.contains(element); });
				GK_CHECK(exon_it != exon_stop, runtime, "Could not map {} \"{}\" to a parent exon", type_name, element);
				element.exon_index = exon_it->index;
			}
			link_exonic_element<T>(table, element, at_index_ptr, el_num_ptr, start_ptr, type_name);
		});
		elements_to_link.clear();
	}

	// Handle ensembl GFF3s when cdss may be inconsistently interleaved with
	// exons or negative strands may be in downstream ordering.
	// This relies on collecting all the exons/cdss in a transcript and then
	// processing when fully known (we get to a new transcript).
	void link(packed_exons_t& exons_to_link, packed_cdss_t& cdss_to_link, packed_utrs_t& utr5s_to_link, packed_utrs_t& utr3s_to_link, packed_utrs_t& generic_utrs_to_link)
	{
		static constexpr auto upstream_ordering = [](auto& lhs, auto& rhs) { return lhs.upstream_of(rhs); };
		sort(begin(exons_to_link), end(exons_to_link), upstream_ordering);
		for_each(cbegin(exons_to_link), cend(exons_to_link), [this](const auto& x) { link(x); });
		exons_to_link.clear();

		if (cdss_to_link.empty() && utr5s_to_link.empty() && utr3s_to_link.empty() && generic_utrs_to_link.empty()) {
			return;
		}

		sort(begin(cdss_to_link), end(cdss_to_link), upstream_ordering);

		for (packed_utr& utr : generic_utrs_to_link) {
			const packed_cds first_cds = cdss_to_link.front();
			const packed_cds last_cds = cdss_to_link.back();
			GK_CHECK(utr.upstream_of(first_cds) || utr.dnstream_of(last_cds), value,
				"Found utr {} in middle of cds", utr);
			if (utr.upstream_of(first_cds)) {
				utr5s_to_link.push_back(utr);
			}
			if (utr.dnstream_of(last_cds)) {
				utr3s_to_link.push_back(utr);
			}
		}
		generic_utrs_to_link.clear();

		const packed_tran& parent_tran = (!cdss_to_link.empty()) ? trans[cdss_to_link[0].tran] :
				(!utr5s_to_link.empty()) ? trans[utr5s_to_link[0].tran] :
				trans[utr3s_to_link[0].tran];
		const auto         exon_start  = next(exons.begin(), parent_tran.exon0);
		const auto         exon_stop   = next(exon_start, parent_tran.num_exons);

		const auto needs_remap = [exon_start](auto type_name, auto& exon_child) {
			bool remap = exon_child.exon_index == cds_exon_index_invalid;
			if constexpr (gencode_v35_lift37_0based_bug)
				remap = remap || !next(exon_start, exon_child.exon_index)->contains(exon_child);
			else if (!remap)
				GK_CHECK(next(exon_start, exon_child.exon_index)->contains(exon_child), runtime,
						 "{} not contained in parent exon: \"{}\"", type_name, exon_child);
			return remap;
		};

		link_exonic_elements<packed_cds>(cdss_to_link, upstream_ordering, exon_start, exon_stop, needs_remap, cdss,
										   &packed_exon::cds, &packed_tran::num_cdss, &packed_tran::cds0, "CDS");
		link_exonic_elements<packed_utr>(utr5s_to_link, upstream_ordering, exon_start, exon_stop, needs_remap, utr5s,
										   &packed_exon::utr5, &packed_tran::num_utr5s, &packed_tran::utr50, "UTR5");
		link_exonic_elements<packed_utr>(utr3s_to_link, upstream_ordering, exon_start, exon_stop, needs_remap, utr3s,
											&packed_exon::utr3, &packed_tran::num_utr3s, &packed_tran::utr30, "UTR3");
	}

	//////////////////////////////////////////////////////////////////

	void load_ucsc_refseq(const char* ucsc_db_dir)
	{
		if (verbose)
			print("Loading UCSC RefSeq from {}...\n", ucsc_db_dir);

		// Since RefSeq exons don't have IDs, make exon.id share the same empty string (None)
		auto empty_exon_id = (offset_t)exons.curr_aux();
		exons.add_aux("");

		biotype_t type_protein_coding = as_biotype("protein_coding");
		biotype_t type_non_coding     = as_biotype("non_coding");

		///////////////////////////////////////////////////////////
		// Load refLink.
		// refLink contains entries for several species, not just human, so we
		// can't populate the `genes` table directly with refLink. Instead, we
		// create an over-complete table indexed by transcript ID (e.g. "NM_007300").
		// Later, when we load refGene (which is human-only), we'll use the
		// transcript IDs encountered there to know which subset of refLink belongs
		// in the `genes` table.
		map<string, reflink_entry> reflink;
		string reflink_file = format("{}/refLink.txt.gz", ucsc_db_dir);
		if (!is_file(reflink_file))
			reflink_file = format("{}/refLink.txt", ucsc_db_dir);   // mini refLink.txt for testing is not zipped
		for (zline_reader lr(reflink_file); !lr.done(); ++lr) {
			// Split tab-separated refLink fields.
			// The columns at this level correspond to the UCSC table schema:
			// [ 0] string name;         // Name displayed in UI
			// [ 1] string product;      // Name of protein product
			// [ 2] string mrnaAcc;      // mRNA accession
			// [ 3] string protAcc;      // protein accession
			// [ 4] uint geneName;       // pointer to geneName table
			// [ 5] uint prodName;       // pointer to prodName table
			// [ 6] uint locusLinkId;    // Entrez ID (formerly LocusLink ID)
			// [ 7] uint omimId;         // OMIM ID
			static constexpr int num_cols = 8;
			string_view cols[num_cols];
			try {
				GK_CHECK(split_view(lr.line(), '\t', cols, num_cols) == num_cols, value,
					"Expected {} tab-separated columns", num_cols);

				reflink_entry& entry = reflink[string(cols[2])];
				GK_CHECK(entry.name.empty(), value, "Duplicate transcript ID \"{}\" found in refLink",
						cols[2]);
				entry.name          = cols[0];
				entry.product = cols[1];
				entry.transcript_id = cols[2];
				entry.protein_id = cols[3];
				entry.entrez_id = cols[6];
			}
			GK_RETHROW("In refLink file: {}:{}", reflink_file, lr.line_num());
		}

		///////////////////////////////////////////////////////////////////////
		// Load refGene, pre-grouping the lines (transcripts) by the corresponding
		// entrez_id, so that transcripts belonging to the same gene will be
		// created in consecutive order.
		using refgene_group_t = multiset<refgene_entry, compare_refgene_entry>;
		map<string, refgene_group_t> refgene;  // { entrez_id : {refgene entries with that entrez_id} }
		string refgene_file = format("{}/refGene.txt.gz", ucsc_db_dir);
		if (!is_file(refgene_file))
			refgene_file = format("{}/refGene.txt", ucsc_db_dir);  // mini refGene.txt for testing is not zipped
		for (zline_reader lr(refgene_file); !lr.done(); ++lr) {
			// Split tab-separated refGene transcript fields.
			// The columns at this level correspond to the UCSC table schema:
			// [ 0] uint   bin;          // No description
			// [ 1] string name;         // Name of gene (actually the transcript_id from GTF)
			// [ 2] string chrom;        // Reference sequence chromosome or scaffold
			// [ 3] char   strand;       // + or - for strand
			// [ 4] uint   txStart;      // Transcription start position (or end position for minus strand item)
			// [ 5] uint   txEnd;        // Transcription end position (or start position for minus strand item)
			// [ 6] uint   cdsStart;     // Coding region start (or end position for minus strand item)
			// [ 7] uint   cdsEnd;       // Coding region end (or start position for minus strand item)
			// [ 8] uint   exonCount;    // Number of exons
			// [ 9] uint[exonCount] exonStarts;  // Exon start positions (or end positions for minus strand item)
			// [10] uint[exonCount] exonEnds;    // Exon end positions (or start positions for minus strand item)
			// [11] uint score;          // Score
			// [12] string name2;        // Alternate name (e.g. gene_id from GTF)
			// [13] string cdsStartStat; // enum('none','unk','incmpl','cmpl')
			// [14] string cdsEndStat;   // enum('none','unk','incmpl','cmpl')
			// [15] int[exonCount] exonFrames;   // Exon frame {0,1,2}, or -1 if no frame for exon
			const int num_cols = 16;
			string_view cols[num_cols];
			auto line = lr.line();
			try {
				GK_CHECK(split_view(line, '\t', cols, num_cols) == num_cols, value,
						 "Expected {} tab-separated columns on line beginning \"{}\"", num_cols, line);

				// Pull out fields at tab-separated level
				refgene_entry e;
				e.tran_id   = cols[1];
				e.chrom     = chrom_names.as_chrom(cols[2]);
				e.strand    = as_strand(cols[3]);
				e.start     = as_pos(cols[4]);
				e.end       = as_pos(cols[5]);
				e.cds_start = as_pos(cols[6]);  // already 0-based, inclusive
				e.cds_end   = as_pos(cols[7]);  // already 0-based, exclusive
				e.num_exons = as_int(cols[8]);
				e.exon_starts = cols[9];
				e.exon_ends   = cols[10];
				GK_CHECK(e.num_exons > 0, value, "Expected positive exon_count but found {}", e.num_exons);
				GK_CHECK(e.tran_id.length() > 3 && e.tran_id[2] == '_', value, "Invalid transcript ID {}", e.tran_id);

				// Set the entrez_id member of this refgene entry, for sake of ordering/grouping
				auto reflink_entry = reflink.find(e.tran_id);
				GK_CHECK(reflink_entry != reflink.end(), value, "Could not find Entrez ID for {}", e.tran_id);
				e.entrez_id = reflink_entry->second.entrez_id;

				// Insert the refGene entry under the entrez_id key
				refgene[e.entrez_id].insert(e);
			}
			GK_RETHROW("In refGene file: {}:{}", refgene_file, lr.line_num());
		}

		///////////////////////////////////////////////////////////////////////
		// Process the refGene table with transcript-, exon-, and cds-level attributes
		// For each transcript, we'll ensure that its parent gene is in the `genes`
		// table, and then build out its exons, introns, and cds elements.
		vector<string_view> exon_starts;
		vector<string_view> exon_ends;
		// Structures to look up element indices by their ID string while loading
		map<string, index_t> tran_by_id;
		for (const auto& [ignore, group] : refgene) {
			// We'll need to know the entrez_id and the gene name of the gene group we're working with.
			// Look them up from the transcript
			auto gene_reflink_it = reflink.find(group.begin()->tran_id);
			GK_CHECK(gene_reflink_it != reflink.end(), key, "Could not find transcript '{}' in refLink table", group.begin()->tran_id);
			const string& gene_name = gene_reflink_it->second.name;
			const string& gene_id = gene_reflink_it->second.entrez_id;

			//////////////////////////////////////////////////
			// PHASE 1: Separate the transcripts into candidate gene clusters, where
			// groups of transcripts with large gaps are separated.
			// `clusters` will be a set of non-overlapping intervals, where each cluster
			// interval is constructed from a mutually-exclusive set of transcript intervals.
			vector<interval_t> clusters;
			for (const auto& j : group) {
				// Check if interval j can be added to the most recent cluster. This loop only works
				// because the group elements are sorted by (chrom, strand, start)
				if (!clusters.empty()) {
					interval_t& c = clusters.back();
					if (c.strand == j.strand && c.chrom == j.chrom && j.start < c.end() + 200000) {
						// Expand this cluster to include transcript interval j.
						GK_ASSERT(c.start() <= j.start);
						c = interval_t::from_dna0(j.chrom, c.start(), max(c.end(), j.end), c.strand, refg());
						continue;
					}
				}

				// Otherwise, start a new cluster
				clusters.push_back(interval_t::from_dna0(j.chrom, j.start, j.end, j.strand, refg()));
			}

			//////////////////////////////////////////////////
			// PHASE 2: For each cluster, create a gene entry and load it with
			// transcripts falling entirely within.
			for (auto& cluster : clusters) {
				// Look up a reflink row for this transcript group, which all share the same entrez_id
				// Create the gene entry.
				index_t gene_index = genes.add_elem();
				packed_gene& gene = genes[gene_index];
				gene.as_interval() = cluster;
				// Gene-specific fields
				gene.level  = invalid_elevel;   // Not specified by UCSC RefSeq
				gene.type   = invalid_biotype;  // Not specified by UCSC RefSeq
				gene.num_trans = 0;             // Start as 0. Incremented as transcripts are read in. Child transcripts always defined after their parent gene.
				gene.tran0 = invalid_index;
				// Gene-specific variable-length fields to the aux able
				gene.aux = (offset_t)genes.curr_aux();
				genes.add_aux(gene_id);
				genes.add_aux(gene_name);

				//////////////////////////////////////////////////
				// PHASE 4: Add each transcript contained within the final gene_interval
				for (const auto& e : group) {
					// If this particular transcript's interval is not within the pre-determined
					// gene interval (i.e. a gene cluster, possibly an explicitly-preferred one),
					// then we discard this transcript since it is probably a duplicate positioned
					// at some irrelevant locus due to UCSC's crappy BLAST method.
					interval_t interval = interval_t::from_dna0(e.chrom, e.start, e.end, e.strand, refg());
					if (!interval.within(gene))
						continue;

					// Pull out fields into local variables, for convenience
					const string& tran_id = e.tran_id;
					int num_exons = e.num_exons;
					strand_t strand = e.strand;

					// Determine biotype as non_coding or protein_coding
					biotype_t type = invalid_biotype;
					if      (tran_id[1] == 'M') type = type_protein_coding;
					else if (tran_id[1] == 'R') type = type_non_coding;
					else GK_THROW(value, "Invalid transcript ID {}", tran_id);

					// Ensure the gene biotype corresponds to the transcript biotypes, with priority None < non_coding < protein_coding
					if (type == type_protein_coding || (type == type_non_coding && gene.type != type_protein_coding))
						gene.type = type;

					auto tran_reflink_it = reflink.find(tran_id);
					GK_CHECK(tran_reflink_it != reflink.end(), key, "Could not find transcript '{}' in refLink table", tran_id);

					//////////////////////////////////////////////////////////
					// Add the transcript entry
					index_t tran_index = trans.add_elem();
					packed_tran& tran = trans[tran_index];
					tran_by_id[tran_id] = tran_index;
					tran.as_interval() = interval;
					// Transcript-specific fields
					// level/tsl not specified by UCSC RefSeq
					tran.type   = type;
					tran.gene = gene_index;
					// Transcript-specific variable-sized fields in aux table
					tran.aux = (offset_t)trans.curr_aux();
					trans.add_aux(tran_id);
					trans.add_aux("");   // 'ccds_id' not specified by UCSC RefSeq.
					trans.add_aux(tran_reflink_it->second.protein_id);
					trans.add_aux(tran_reflink_it->second.product);

					// Increment the parent gene's transcript counter
					gene.num_trans++;
					if (gene.tran0 == invalid_index)
						gene.tran0 = tran_index;
					GK_ASSERT(tran_index == gene.tran0 + gene.num_trans-1);

					///////////////////////////////////////////////////
					// parse the comma-separated exon arrays.
					// These always end in an extra command, e.g. 4 exon starts could be string "2000,2500,3000,3500,",
					exon_starts.resize(num_exons);
					exon_ends.resize(num_exons);
					// trim the trailing ","
					GK_CHECK(split_view(string_view(e.exon_starts).substr(0, e.exon_starts.size() - 1), ',',
										&exon_starts[0], num_exons)
								 == num_exons,
							 value, "Expected {} exon starts for {}", num_exons, tran_id);
					GK_CHECK(split_view(string_view(e.exon_ends).substr(0, e.exon_ends.size() - 1), ',', &exon_ends[0],
										num_exons)
								 == num_exons,
							 value, "Expected {} exon ends for {}", num_exons, tran_id);
					if (strand == neg_strand) {
						std::reverse(exon_starts.begin(), exon_starts.end());  // Make sure we add the exons in sense-strand order.
						std::reverse(exon_ends.begin(),   exon_ends.end());
					}
					int cds_len = 0;
					for (int i = 0; i < num_exons; ++i) {
						pos_t exon_start = as_pos(exon_starts[i]);
						pos_t exon_end   = as_pos(exon_ends[i]);

						///////////////////////////////////////////////////////
						// Add exon entry to the `exons` table.
						index_t exon_index = exons.add_elem(
							{ interval_t::from_dna0(tran.chrom, exon_start, exon_end, tran.strand, tran.refg) });
						packed_exon& exon = exons[exon_index];
						// Exon-specific fields
						exon.is_last = false;
						exon.index  = int_cast<elem_num>(i);
						exon.tran = tran_index;
						exon.cds  = invalid_index; // May be set to a CDS index later
						exon.utr5 = invalid_index;
						exon.utr3 = invalid_index;
						// Exon-specific variable-sized fields in aux table
						exon.id = empty_exon_id;

						// Make note of first exon's index in the parent struct
						tran.num_exons++;
						if (tran.exon0 == invalid_index)
							tran.exon0 = exon_index;
						GK_ASSERT(exon_index == tran.exon0 + tran.num_exons-1);

						////////////////////////////////////////////////////////////////
						// Add upstream intron entry to the `intrs` table, if applicable
						if (exon.index > 0) {
							// Find 3p end of previous exon and 5p end of current exon
							pos_t prev_pos3 = exons[tran.exon0 + exon.index-1].pos3;
							pos_t curr_pos5 = exon.pos5;

							// 1-based closed to 0-based closed
							index_t intr_index = intrs.add_elem({ tran.with_pos(prev_pos3, curr_pos5).expand(-1, -1) });
							packed_intr& intr = intrs[intr_index];
							// Intron-specific fields
							intr.is_last = false;
							intr.index  = exon.index-1;  // intron 0 precedes exon 1
							intr.tran = tran_index;

							// Make note of first intron's index in the parent struct
							if (tran.intr0 == invalid_index)
								tran.intr0 = intr_index;
							GK_ASSERT(intr_index == tran.intr0 + tran.num_exons-2);
						}

						if (e.cds_start < e.cds_end && exon_end >= e.cds_start && exon_start <= e.cds_end) {
							// Figure out CDS start/end within this exon
							pos_t exon_cds_start = max(exon_start, e.cds_start);
							pos_t exon_cds_end   = min(exon_end, e.cds_end);

							// Calculate the phase of this exon, i.e. how many nt would be trimmed to
							// put the coding sequence in frame 0
							uint8_t phase = int_cast<uint8_t>(3 - umod(cds_len, 3));
							if (phase == 3)
								phase = 0;
							cds_len += exon_cds_end-exon_cds_start;

							////////////////////////////////////////////////////////////////
							// Add cds entry to the `cdss` table, if applicable
							index_t cds_index = cdss.add_elem({{ interval_t::from_dna0(
								tran.chrom, exon_cds_start, exon_cds_end, tran.strand, tran.refg) }});
							packed_cds& cds   = cdss[cds_index];
							// CDS-specific fields
							cds.is_first = false;
							cds.is_last  = false;
							cds.phase    = phase;
							cds.exon_index = exon.index;
							cds.tran  = tran_index;

							// Make note of CDS's index in the parent exon and parent transcript
							exon.cds = cds_index;
							tran.num_cdss++;
							if (tran.cds0 == invalid_index)
								tran.cds0 = cds_index;
							GK_ASSERT(cds_index == tran.cds0 + tran.num_cdss-1);
						}
					}
				}

				GK_CHECK(gene.num_trans > 0, value, "Gene {} was assigned 0 transcripts -- logic bug!", gene_name);
			}
		}

		add_utrs();
	}

	void load_ncbi_refseq(const std::string& infile)
	{
		using namespace std;

		if (verbose)
			print("Loading {}...\n", infile);

		// NCBI doesn't provide exon ids, pool the empty string
		const auto exon_empty_id = int_cast<offset_t>(exons.add_aux(""));

		// NCBI RefSeq <v109 have interleaving (wrt genes/trans/exon/cds)
		// so store the records and only resolve when we have a clean break
		// (changing chromosome or ### directive)
		struct gene_record {
			packed_gene packed;
			string      id;
			string      name;
			string      gff_id;
		};
		struct tran_record {
			packed_tran packed;
			string      id;
			string      gff_id;
		};
		struct cds_record {
			packed_cds packed;
			string     product;
			string     protein_id;
		};
		vector<gene_record>                        gene_records;
		unordered_map<string, vector<tran_record>> tran_by_parent;
		unordered_map<string, vector<packed_exon>> exon_by_parent;
		unordered_map<string, vector<cds_record>> cds_by_parent;

		auto resolve_forward_references = [&]() {
			string protein_id; // propagate the CDS protein_id up to the transcript
			// use the CDS product (eg, "breast cancer type 1 susceptibility protein isoform 2") rather than the
			// transcript product (eg, "BRCA1, DNA repair associated, transcript variant 2")
			string product;
			// NCBI links exons and CDS to either the gene or transcript;
			// for now, only include features with a proper transcript_id.
			for (auto gene : gene_records) {
				auto i_gene = genes.add_elem();
				for (auto tran : tran_by_parent[gene.gff_id]) {
					auto i_tran = trans.add_elem();

					for (auto exon : exon_by_parent[tran.gff_id]) {
						exon.index = tran.packed.num_exons;
						exon.tran  = i_tran;
						exon.id    = exon_empty_id;

						auto i_exon = exons.add_elem(exon);
						if (tran.packed.num_exons == 0) {
							tran.packed.exon0 = i_exon;
						} else {
							auto i_intr = intrs.add_elem({ exon.with_pos(exons[i_exon - 1].pos3, exon.pos5).expand(-1, -1) });
							auto& intr  = intrs[i_intr];
							intr.index  = int_cast<elem_num>(exon.index - 1);
							intr.tran   = i_tran;
							if (tran.packed.num_exons == 1)
								tran.packed.intr0 = i_intr;
						}
						++tran.packed.num_exons;
					}

					protein_id.clear();
					product.clear();

					auto  one_cds_per_exon = true;
					auto  it_exon          = next(cbegin(exons), tran.packed.exon0);
					auto  end_exon         = next(it_exon, tran.packed.num_exons);
					auto& tran_cdss        = cds_by_parent[tran.gff_id];
					for (auto& cds : tran_cdss) {
						if (protein_id.empty())
							protein_id = cds.protein_id;
						if (product.empty())
							product = cds.product;

						cds.packed.tran = i_tran;
						it_exon = find_if(it_exon, end_exon, [&cds](auto exon) { return exon.contains(cds.packed); });
						if (it_exon != end_exon) {
							cds.packed.exon_index = it_exon->index;
							++it_exon;
						} else {
							one_cds_per_exon = false;
							break;
						}
					}

					if (one_cds_per_exon) {
						// only support simple CDS representations for now: Entrez IDs 51686, 23089, 4947, 4946
						// have programmed frameshifts but aren't implicated in any Mendialian diseases.
						// NOTE: GENCODE and UCSC decided to split exons over the frameshift gap rather than
						// exclusively over spliced regions.
						for (const auto& cds : tran_cdss) {
							auto i_cds = cdss.add_elem(cds.packed);
							if (tran.packed.num_cdss == 0)
								tran.packed.cds0 = i_cds;
							++tran.packed.num_cdss;

							GK_ASSERT(cds.packed.exon_index != cds_exon_index_invalid); // only one exon per CDS

							auto& exon = exons[tran.packed.exon0 + cds.packed.exon_index];
							GK_ASSERT(exon.cds == invalid_index); // only one CDS per exon
							exon.cds = i_cds;
						}
					}

					tran.packed.type = gene.packed.type;
					tran.packed.gene = i_gene;
					tran.packed.aux  = int_cast<offset_t>(trans.add_aux(tran.id));
					trans.add_aux("");
					trans.add_aux(protein_id);
					trans.add_aux(product);
					trans[i_tran] = tran.packed;

					auto num_tran = gene.packed.num_trans;
					if (num_tran == 0)
						gene.packed.tran0 = i_tran;
					gene.packed.num_trans = num_tran + 1;
				}

				gene.packed.aux = int_cast<offset_t>(genes.add_aux(gene.id));
				genes.add_aux(gene.name);
				genes[i_gene] = gene.packed;
			}

			cds_by_parent.clear();
			exon_by_parent.clear();
			tran_by_parent.clear();
			gene_records.clear();
		};

		zline_reader lr{infile};
		vector<string_view> attrs;
		optional<chrom_t>   last_chrom;
		vector<string_view> dbxrefs;
		set<string, less<>> unknown_gff_types;
		for (; !lr.done(); ++lr) {
			auto line = lr.line();
			if (line == "###") {
				last_chrom.reset();  // mark to resolve forward references
				continue;
			}
			if (line[0] == '#')
				continue;

			try {
				gff_record record{lr.line(), &chrom_names};
				auto       chrom = record.seqid();
				if (!last_chrom || chrom != *last_chrom) {
					resolve_forward_references();
					last_chrom = chrom;
				}

				record.set_attributes(attrs);
				auto gff_id = get_attr(attrs, "ID", "");
				if (record.type() == "gene" || startswith(gff_id, "gene")) {
					split_view(get_attr(attrs, "Dbxref"), ',', dbxrefs);
					gene_records.push_back(
						{ { record.as_interval(), ncbi_source_as_elevel(record.source()),
							as_biotype_compatible(get_attr(attrs, "gene_biotype")) },
						  string(get_attr(dbxrefs, "GeneID", nullptr, ':')),
						  string(get_attr(attrs, "gene")),
						  string(gff_id) });
				} else if (startswith(gff_id, "rna")) {
					// ignore any unstable (across versions) transcripts for now.
					if (auto transcript_id = get_attr(attrs, "transcript_id", ""); !std::empty(transcript_id)) {
						auto parent_id = get_attr(attrs, "Parent");
						GK_ASSERT(!contains(parent_id, ",")); // multiple parents not implemented
						tran_record tran{
							{}, // zero init the bit pad in the bitfield
							string(transcript_id),
							string(gff_id),
						};
						tran.packed = { record.as_interval() };
						tran.packed.level = ncbi_source_as_elevel(record.source());
						tran_by_parent[string(parent_id)].push_back(tran);
					}
				} else if (record.type() == "exon") {
					auto parent_id = get_attr(attrs, "Parent");
					GK_ASSERT(!contains(parent_id, ",")); // multiple parents not implemented
					exon_by_parent[string(parent_id)].push_back({ record.as_interval(), 0,
																  transcript_exon_index_max, invalid_index,
																  invalid_index, invalid_index,
																  invalid_index, offset_t(-1) });
				} else if (record.type() == "CDS") {
					auto parent_id = get_attr(attrs, "Parent");
					GK_ASSERT(!contains(parent_id, ",")); // multiple parents not implemented
					cds_by_parent[string(parent_id)].push_back(
						{ {{ record.as_interval(), 0, 0, record.phase(), cds_exon_index_invalid, invalid_index}},
						  string(get_attr(attrs, "product", "")),
						  string(get_attr(attrs, "protein_id", "")) });
				} else if (auto loc = unknown_gff_types.lower_bound(record.type());
						   loc == cend(unknown_gff_types) || *loc != record.type()) {
					unknown_gff_types.emplace_hint(loc, record.type());
					print("Unknown GFF type={} on line {}.\n", record.type(),
							lr.line_num());
				}
			}
			GK_RETHROW("In GFF3 file: {}:{}", infile, lr.line_num());

			if (verbose && lr.line_num() % 8000 == 1) {
				auto sum_value = [](auto init, const auto& x) { return init + x.second.size(); };
				print("  {} genes, {} transcripts, {} exons...\n", genes.size() + gene_records.size(),
					  accumulate(cbegin(tran_by_parent), cend(tran_by_parent), trans.size(), sum_value),
					  accumulate(cbegin(exon_by_parent), cend(exon_by_parent), exons.size(), sum_value));
			}
		}
		resolve_forward_references();
		add_utrs();

		if (verbose)
			print("  {} genes, {} transcripts, {} exons...\n", genes.size(), trans.size(), exons.size());
	}

	std::vector<std::string> dump(const char* outfile)
	{
		// First, sanity-check the tables and ensure aux data is allocated
		finalize_tables();

		// Dump the DGANNO file using the builder objects we've constructed from the GFF3 file
		auto        versioned_filename = format("{}.v{}.dganno", outfile, genome_anno::binary_version());
		binary_file out(versioned_filename.c_str(), "wb");

		// Write the signature
		out.write(c_dganno_sig);
		out.write(c_dganno_ver);

		// The primary access patterns to annotations are:
		// 1. linear search over features to filter them
		// 2. random access a feature since it's a known feature of interest
		// 3. binary search over the feature indices via an interval
		//
		// Because of the way the packed/unpacked structs are arranged, loading
		// a feature will likely expand all linked feature fields.
		// This means there's not much benefit in lazy loading single
		// feature tables, and we can lay out all the feature tables to be
		// loaded in order, all at once.
		//
		// For binary searches, the search will be done on a single feature,
		// so there's no benefit to outlining all feature's indices together
		// instead of keeping them adjacent to their element data.

		// Write each table in a specific order (see genome_anno::open)
		genes.dump(out);
		trans.dump(out);
		exons.dump(out);
		intrs.dump(out);
		cdss.dump(out);
		utr5s.dump(out);
		utr3s.dump(out);
		// <--- INSERT NEW TABLE DUMPS HERE

		auto cfg_filename = format("{}.cfg", outfile);
		std::ofstream(cfg_filename) << "refg=" << chrom_names.refg_name() << std::endl;

		if (verbose)
			print("Wrote {}\n", versioned_filename);

		return {versioned_filename, cfg_filename};
	}

private:
	struct reflink_entry {
		string name;
		string product;
		string transcript_id;
		string protein_id;
		string entrez_id;
	};

	struct refgene_entry {
		string     entrez_id;
		string     tran_id;
		chrom_t    chrom;
		strand_t   strand;
		pos_t      start;
		pos_t      end;
		pos_t      cds_start;
		pos_t      cds_end;
		int        num_exons;
		string     exon_starts;
		string     exon_ends;
	};

	struct compare_refgene_entry {
		// When adding UCSC refgene transcripts to a parent gene entry,
		// sorting them by (chrom,strand,start) helps to detect flakey
		// transcripts that are mapped to multiple disconnected sites.
		// Sorting them also by tran_id helps keep order of transcripts
		// consistent across multiple versions of gene, since UCSC transcript
		// IDs are not unique (same IDs can appear in two 'copies' of gene).
		bool operator()(const refgene_entry& a, const refgene_entry& b) const
		{
			if (a.chrom  != b.chrom)  return a.chrom < b.chrom;
			if (a.strand != b.strand) return a.strand < b.strand;
			if (a.start  != b.start)  return a.start < b.start;
			return a.tran_id.compare(b.tran_id) < 0;
		}
	};

	void finalize_tables()
	{
		//////////////////////////////////////////////////////////////////////////
		// Store auxiliary indices and strings for the gene and transcript entries.
		//////////////////////////////////////////////////////////////////////////

		// Per-transcript post-processing
		for (const auto& tran : trans) {
			if (tran.num_exons > 0) exons[tran.exon0 + tran.num_exons-1].is_last = true;
			if (tran.num_exons > 1) intrs[tran.intr0 + tran.num_exons-2].is_last = true;
			if (tran.num_cdss  > 0) {
				cdss[tran.cds0                  ].is_first = true;
				cdss[tran.cds0 + tran.num_cdss-1].is_last  = true;
			}
			if (tran.num_utr5s  > 0) {
				utr5s[tran.utr50                     ].is_first = true;
				utr5s[tran.utr50 + tran.num_utr5s - 1].is_last  = true;
			}
			if (tran.num_utr3s  > 0) {
				utr3s[tran.utr30                     ].is_first = true;
				utr3s[tran.utr30 + tran.num_utr3s - 1].is_last  = true;
			}
		}

		// <--- INSERT NEW AUXILIARY TABLE CASES HERE
	}

	// Add UTRs to transcripts where they are not provided as part of the annotation
	void add_utrs() {
		for (int tran_idx = 0; tran_idx < trans.size(); tran_idx++) {
			const packed_tran &tran = trans[tran_idx];
			if (verbose) {
				print("transcript: {}\n", tran);
			}
			add_utrs_to_single_transcript(tran, tran_idx);
		}
	}

	void add_utrs_to_single_transcript(const packed_tran& tran, const int tran_idx)
	{
		if (tran.num_cdss == 0) {
			if (verbose) {
				print("non-coding transcript {}\n", tran);
			}
			return;
		}

		for (index_t i = tran.exon0; i <= tran.exon0 + cdss[tran.cds0].exon_index; ++i) {
			auto exon = exons[i];
			auto interval = exon.as_interval();
			if (exon.cds != invalid_index) {
				auto cds = cdss[exon.cds];
				interval = interval.with_pos(interval.end5().start(), cds.end5().start()).shrink(0, 1);
			}
			if (!interval.empty()) {
				link_exonic_element<packed_utr>(utr5s, { {interval, false, false, 0, exon.index, tran_idx} },
												&packed_exon::utr5, &packed_tran::num_utr5s, &packed_tran::utr50, "UTR5");
			}
		}
		for (index_t i = tran.exon0 + cdss[tran.cds0 + tran.num_cdss - 1].exon_index; i < tran.exon0 + tran.num_exons; ++i) {
			auto exon    = exons[i];
			auto interval = exon.as_interval();
			if (exon.cds != invalid_index) {
				auto cds = cdss[exon.cds];
				interval = interval.with_pos(cds.end3().start(), interval.end3().start()).shrink(1, 0);
			}
			if (!interval.empty()) {
				link_exonic_element<packed_utr>(utr3s, { {interval, false, false, 0, exon.index, tran_idx} },
												&packed_exon::utr3, &packed_tran::num_utr3s, &packed_tran::utr30, "UTR3");
			}
		}
	}

	template <class T> T process_exonic_component(string_view *cols, const string &last_tran_id, const auto parsed_interval) {
		static constexpr string_view transcript_prefix = "transcript:";
		vector<string_view> attrs;
		split_view(cols[8], ';', attrs);
		auto tran_id = get_attr(attrs, "Parent");
		// only available in gencode
		auto exon_index = get_attr(attrs, "exon_number", cds_exon_index_invalid_cstr);
		if constexpr (gencode_v35_lift37_0based_bug) {
			/* can't rely on this being present to increment to 1-based, in case the first exon is purely UTR */
			if (exon_index == "0")
				exon_index = cds_exon_index_invalid_cstr;
		}

		tran_id = skip_prefix(tran_id, transcript_prefix);

		// If this fails: use a more heavy-weight map recording all prior transcripts
		auto type_name = cols[1];
		GK_CHECK(!trans.empty() && last_tran_id == tran_id, key,
				 "Preceding transcript '{}' is not the parent of {} '{}'", last_tran_id,
				 type_name, get_attr(attrs, "ID"));
		index_t tran_index = trans.size() - 1;

		// Fill the element entry fields
		T res = {{parsed_interval()}};
		res.is_first = false;
		res.is_last  = false;
		res.tran  = tran_index;
		res.exon_index = exon_index != cds_exon_index_invalid_cstr
						 ? as_elem_num(exon_index, cds_exon_index_max)
						 : cds_exon_index_invalid;
		return res;
	}

	bool verbose;
	chrom_names_t chrom_names;

	// Cross-referenced tables that comprise the final annotations.
	gene_table::builder genes;
	tran_table::builder trans;
	exon_table::builder exons;
	intr_table::builder intrs;
	cds_table::builder  cdss;
	utr_table::builder utr5s;
	utr_table::builder utr3s;
	// <--- INSERT NEW ANNOTATION BUILDER INSTANCES HERE

	// <--- INSERT NEW LOAD HELPER INSTANCES HERE
};

int genome_anno::binary_version() { return c_dganno_ver; }

std::vector<std::string> genome_anno::build_gencode(const char* infile, const char* outfile, const genome_t& genome)
{
	builder b{genome};
	b.load_gencode(infile);
	return b.dump(outfile);
}

std::vector<std::string> genome_anno::build_ucsc_refseq(const char* ucsc_db_dir, const char* outfile, const genome_t& genome)
{
	builder b{genome};
	b.load_ucsc_refseq(ucsc_db_dir);
	return b.dump(outfile);
}

std::vector<std::string> genome_anno::build_ncbi_refseq(const char* infile, const char* outfile, const genome_t& genome)
{
	builder b{genome};
	b.load_ncbi_refseq(infile);
	return b.dump(outfile);
}

void genome_anno::open()
{
	GK_CHECK(!is_open(), runtime, "genome_anno::open() already opened");
	GK_CHECK(!_sourcefile.empty(), value, "Genome is not a registered annotation");

	// Resolve the source file, possibly downloading it
	try {
		auto resolved = resolve_datafile_path(_sourcefile);
		if (std::filesystem::exists(resolved)) {
			_sourcefile = std::move(resolved);
		}
	} catch (const value_error& e) {
		print("{}\n", e.what());
	}

	// Memory map the source file
	_fmap.open(_sourcefile);

	// Read the DGANNO file signature
	unsigned short sig;
	unsigned short ver;
	_fmap.read(sig);
	_fmap.read(ver);
	GK_CHECK(sig == c_dganno_sig, file, "Expected valid DGANNO file signature {:x} but found {:x}.", c_dganno_sig, sig);
	GK_CHECK(ver == c_dganno_ver, file, "Expected DGANNO file version {:x} but found {:x}.", c_dganno_ver, ver);

	// Read and map each table in a specific order (see genome_anno::build_gencode)
	_genes.load(_fmap);
	_trans.load(_fmap);
	_exons.load(_fmap);
	_intrs.load(_fmap);
	_cdss.load(_fmap);
	_utr5s.load(_fmap);
	_utr3s.load(_fmap);
	// <--- INSERT NEW TABLE LOADS HERE
}

END_NAMESPACE_GK
