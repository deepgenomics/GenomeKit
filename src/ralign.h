/*
Copyright (C) 2016-2023 Deep Genomics Inc. All Rights Reserved.
*/
#pragma once
#ifndef __GENOME_KIT_RALIGN_H__
#define __GENOME_KIT_RALIGN_H__

#include "file.h"
#include "sam_line_parser.h"
#include "table.h"
#include "util.h"
#include "variant_table.h"
#include <map>
#include <string>
#include <vector>

BEGIN_NAMESPACE_GK
using std::vector;
using std::map;
using std::string;

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//      Descriptions
//
//  1. `packed_junction` / `junction_t` represents a skipped region when reads are aligned to a reference genome.
//      Existence of a junction is supported by one or more read alignments.
//  2. `packed_align` / `align_t` represents a single read alignment (i.e. a row in the BAM file).
//      A read alignment consists of potentially alternating sequence of matching and non-matching subsequences
//          -- matching (CIGAR: consecutive M,D,I) : matching region with potential variants (possibly SNVs for M, indels for I or D)
//          -- nonmatching (CIGAR: one N) : skipped region which we call junction
//      Each read alignment is uniquely stored and can be referenced from any junctions that it supports.
//      Each read alignment owns a range of alignment matches.
//      Each read alignment may own a sequence of variants for all of its matching regions
//  3. `packed_align_match` / `align_match_t` represents a matching subsequence of a read alignment
//      A alignment match is affiliated with a single read alignment, and can refer to the owning read alignment
//
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

// forward declare & typedefs
class align_match_table;
class align_table;
class genome_dna;
class genome_t;
class junction_table;
class read_alignments;

struct packed_align_match;
using align_match_range = range_t<const packed_align_match *>;

/////////////////////////////////////////////////////////////////
//      packed_junction / junction_t / junction_table
//
// packed_junction memory layout: offset name (size+padding)
//  0 pos5         (4)
//  4 pos3         (4)
//  8 refg         (8)
// 16 chrom        (4)
// 20 strand       (1)
// 21 num_aligns   (4)
// 25 aux          (4)
// 29 total
//
// Each packed_junction aux pool entry has a dynamic layout:
//
//    struct {
//        index_t*      aligns[n];     // indicies of supporting read alignments, where n = num_aligns
//    };
/////////////////////////////////////////////////////////////////

#pragma pack(push, 1)
struct packed_junction : public interval_t {
	index_t     num_aligns;
	offset_t    aux;
};
#pragma pack(pop)

struct junction_t : public interval_t {
	index_t         num_aligns;         // number of read alignments supporting this junction
	const index_t*  aligns;             // pointer start of indices array in aux mempool

	junction_t(index_t                src, const read_alignments& jralign);
	junction_t(const packed_junction& src, const read_alignments& jralign);
	void unpack_from(const packed_junction& src, const read_alignments& jralign);
};

class junction_table: public interval_table<packed_junction> { NOCOPY(junction_table)
public:
	junction_table(read_alignments& raligns) : raligns(raligns) { }
	read_alignments& raligns;
};

/////////////////////////////////////////////////////////////////
//      packed_align / align_t / align_table
//
// packed_align memory layout: offset name (size+padding)
//  0 pos5              (4)
//  4 pos3              (4)
//  8 refg         		(8)
// 16 chrom        		(4)
// 20 strand       		(1)
// 21 num_matches       (1)
// 22 match0            (4)
// 26 total
//
/////////////////////////////////////////////////////////////////

#pragma pack(push, 1)
struct packed_align : public interval_t {
	unsigned char   num_matches;   // Number of matches (consecutive M,I,D) in this read alignment
	index_t         match0;        // Index into read_alignment::alignments_matches of the first match that belongs to this read alignment
};
#pragma pack(pop)

struct align_t : public interval_t {
	unsigned char     num_matches;
	align_match_range matches;

	align_t(index_t                src, const read_alignments& jralign);
	align_t(const packed_align&    src, const read_alignments& jralign);
	void unpack_from(const packed_align& src, const read_alignments& jralign);
};

class align_table: public interval_table<packed_align> { NOCOPY(align_table)
public:
	align_table(read_alignments& raligns) : raligns(raligns) { }
	read_alignments& raligns;
};


/////////////////////////////////////////////////////////////////
//      packed_align_match / align_match_t / align_match_table
//
// packed_align_match memory layout: offset name (size+padding)
//  0 pos5         (4)
//  4 pos3         (4)
//  8 refg         (8)
// 16 chrom        (4)
// 20 strand       (1)
// 21 num_variants (1)
// 22 aux          (4)
// 26 total
//
// Each packed_align_match aux pool entry has a dynamic layout:
//
//    struct {
//        index_t*      variants[n];     // indices of owning variants, where n = num_variants
//    };
/////////////////////////////////////////////////////////////////
#pragma pack(push, 1)
struct packed_align_match : public interval_t {
	unsigned char num_variants;
	offset_t      aux;
};
#pragma pack(pop)

struct align_match_t : public interval_t {
	unsigned char   num_variants;   // Number of variants owned by the read alignment
	const index_t*  variants;       // Points to start of an array of index (that addresses variant_table) in aux mem pool

	align_match_t(index_t                     src, const read_alignments& jralign);
	align_match_t(const packed_align_match&   src, const read_alignments& jralign);
	void unpack_from(const packed_align_match& src, const read_alignments& jralign);
};

class align_match_table: public interval_table<packed_align_match> { NOCOPY(align_match_table)
public:
	align_match_table(read_alignments& raligns) : raligns(raligns) { }
	read_alignments& raligns;
};

/////////////////////////////////////////////////////////////////
// read_alignments  (junction + aligned reads + alignment matches + variants)
/////////////////////////////////////////////////////////////////

class read_alignments {
public:
	read_alignments();

	const junction_table&      junctions() const;
	const align_table&         alignments() const;
	const align_match_table&   matches() const;
	const variant_table&       variants() const;

	// Sets the original source file (e.g. GFF3 file), but does
	// not actually open it until open() is called.
	void set_source(string sourcefile);
	void open();
	void close();
	INLINE bool is_open() const { return _fmap.is_open(); }
	INLINE void ensure_open() const { if (!is_open()) open_on_demand(); }  // Fast check if file / indices already loaded
	INLINE const string& source() const { return _sourcefile; }
	static int ralign_version();

	// Takes SAM files as input, and outputs a RALIGN file.
	class builder : public sam_line_parser {
	public:
		builder(const char* outfile, const genome_t& genome);
		void add(const char* infile);
		void finalize();

		// variants
		struct c_variant_t {
			interval_t i;
			string ref;
			string alt;
		};

		// junction
		struct b_junction_t { NOCOPY(b_junction_t)
			b_junction_t() = default;

			vector<index_t> alignments;
		};
		// read alignment match
		struct b_align_match_t { NOCOPY(b_align_match_t)
			b_align_match_t() = default;
			b_align_match_t(b_align_match_t&&) = default;
			b_align_match_t& operator=(b_align_match_t&&) = default;

			interval_t i;
			vector<index_t> variants;
		};
		// read alignment
		struct b_alignment_t { NOCOPY(b_alignment_t)
			b_alignment_t() = default;
			b_alignment_t(b_alignment_t&&) = default;
			b_alignment_t& operator=(b_alignment_t&&) = default;

			interval_t         i;
			vector<b_align_match_t>  matches;
		};

		void collect_variant(b_align_match_t& read, const c_variant_t& v);

	protected:
		void process_line(chrom_t chrom, pos_t pos, strand_t strand, strand_t segment_strand, const cigar_op* codes,
						  int m, int n, sam_cols cols) override;

	private:
		binary_file _file;
		long long   _num_reads_loaded{};

		map<interval_t, b_junction_t>   _juncs;
		vector<b_alignment_t>           _aligns;
		// A set of variants for lookup of their index in `_ra_variants_dup`
		using variant_to_index = map<c_variant_t, index_t>;
		const genome_dna*   _dna;
		variant_to_index    _ra_variants;
		vector<c_variant_t> _ra_variants_dup;

		bool _verbose{};
	};

private:
	void open_on_demand() const; // threadsafe version (const)
	mmap_file _fmap;
	string    _sourcefile;
	junction_table    _junctions;
	align_table       _alignments;
	align_match_table _matches;
	variant_table     _variants;
};

bool operator<(const read_alignments::builder::c_variant_t& lhs,
               const read_alignments::builder::c_variant_t& rhs);

END_NAMESPACE_GK

#endif // __GENOME_KIT_RALIGN_H__
