/*
Copyright (C) 2016-2023 Deep Genomics Inc. All Rights Reserved.
*/
#pragma once
#ifndef __GENOME_KIT_JRALIGN_H__
#define __GENOME_KIT_JRALIGN_H__

#include "file.h"
#include "sam_line_parser.h"
#include "table.h"
#include "variant_table.h"
#include <map>
#include <string>
#include <vector>

BEGIN_NAMESPACE_GK

using std::vector;
using std::map;
using std::string;

class genome_dna;
class genome_t;
class jraligns_t;
class jraligns_table;
class junction_read_alignments;

/////////////////////////////////////////////////////////////////
// junction read alignment
/////////////////////////////////////////////////////////////////

struct jralign_t {
public:
	unsigned char left;   // Upstream overhang on reference strand
	unsigned char right;  // Downstream overhang on reference strand
	strand_t strand;  // Strand the read mapped to (pos_strand for reference strand, neg_strand for reverse strand)
	unsigned char num_variants; // Number of variants associated with this read
	const index_t* variants_indices; // Points to start of an array of indices (that addresses variant_table) in aux mem pool

	jralign_t() = default;
	INLINE jralign_t(unsigned char left, unsigned char right, strand_t strand, unsigned char num_variants, const index_t* variants_indices)
		: left(left), right(right), strand(strand), num_variants(num_variants), variants_indices(variants_indices) { }
	INLINE int len() const { return (int)left + (int)right; } // length of the read
};

inline bool operator==(const jralign_t& a, const jralign_t& b) { return a.left == b.left && a.right == b.right && a.strand == b.strand; }
inline bool operator!=(const jralign_t& a, const jralign_t& b) { return !(a == b); }


/////////////////////////////////////////////////////////////////
// jraligns_t (junction read alignments)
/////////////////////////////////////////////////////////////////

#pragma pack(push, 1)
struct packed_jraligns: public interval_t {
	uint8_t has_variants;  // Non-zero if this junction holds any variants
	unsigned  num_reads;   // Number of reads across this junction
	offset_t aux;          // Offset into junc_reads_table::auxpool for auxiliary data
};
#pragma pack(pop)
// packed_junc_reads memory layout: offset name (size+padding)
//  0 pos5         (4)
//  4 pos3         (4)
//  8 refg         (8)
// 16 chrom        (4)
// 20 strand       (1)
// 21 has_variants (1)
// 22 num_reads    (4)
// 26 aux          (4)
// 30 total
//
// Each packed_jraligns aux pool entry has a dynamic layout:
//
//    struct {
//        unsigned char overhangs[n][2];     // [i][0/1] = 5'/3' overhang of read i, respectively, where n = num_reads
//        unsigned char strands[(n+7)/8];    // bit [i] = alignment strand of read i (0 = aligned to reverse, 1 aligned to reference)
//        unsigned char indices_begin[n+1];  // index to variants_indices for read i
//        index_t*      variants_indices;    // [indices_begin[i]] = array of index_t for read i referencing `_variant_table`
//                                           // indices_begin[i+1] - indices_begin[i] = length of array for read i
//    };
//

class jraligns_t: public interval_t {
public:
	jraligns_t(index_t                src, const jraligns_table& table);
	jraligns_t(const packed_jraligns& src, const jraligns_table& table);

	jralign_t operator[](unsigned i) const;
	INLINE unsigned num_reads() const { return _num_reads; }
	INLINE bool has_variants() const { return _has_variants; }

private:
	bool _has_variants;
	unsigned _num_reads; // Number of reads across this junction
	const unsigned char* _overhangs;
	const unsigned char* _strands;
	const unsigned char* _indices_begin;          // 0 if no variants, set in ctor
	const index_t*       _variants_indices;       // 0 if no variants, set in ctor

	void unpack_from(const packed_jraligns& src, const jraligns_table& table);
};

/////////////////////////////////////////////////////////////////
// junction reads table
/////////////////////////////////////////////////////////////////

class jraligns_table: public interval_table<packed_jraligns> {
public:
	jraligns_table(junction_read_alignments& jraligns) : jraligns(jraligns) { }
	junction_read_alignments& jraligns;
};

/////////////////////////////////////////////////////////////////
// junction read alignments table (junction reads + variant_table)
/////////////////////////////////////////////////////////////////

class junction_read_alignments {
public:
	junction_read_alignments();

	const jraligns_table& juncs() const;
	const variant_table& variants() const;

	// Sets the original source file (e.g. GFF3 file), but does
	// not actually open it until open() is called.
	void set_source(string sourcefile);
	void open();
	void close();
	INLINE bool is_open() const { return _fmap.is_open(); }
	INLINE void ensure_open() const { if (!is_open()) open_on_demand(); }  // Fast check if file / indices already loaded
	INLINE const string& source() const { return _sourcefile; }
	INLINE long long num_reads() const { ensure_open(); return _num_reads; }
	static int jralign_version();

	// Takes SAM files as input, and outputs a RALIGN file.
	class builder : public sam_line_parser {
	public:
		enum class error_handling { error, clamp };

		builder(const char* outfile, const genome_t& genome);
		INLINE void set_min_reads(int value)    { _min_reads = (unsigned)value; }
		INLINE void set_min_overhang(int value) { _min_overhang = value; }
		INLINE void set_include_variants(bool value) { _include_variants = value; }
		INLINE bool include_variants() { return _include_variants; }
		INLINE void set_overhang_error(error_handling value) { _overhang_error = value; }
		void add(const char* infile);
		void finalize();

		struct b_variant_t {
			interval_t i;
			string ref;
			string alt;
		};

		void collect_variant(const b_variant_t& v);

	protected:
		void process_line(chrom_t chrom, pos_t pos, strand_t strand, strand_t segment_strand, const cigar_op* codes,
						  int m, int n, sam_cols cols) override;

	private:
		struct jralign_entry { NOCOPY(jralign_entry)
			// read specific information
			vector<unsigned char> overhangs;
			vector<unsigned char> strands;

			// optional: variant specific information
			vector<unsigned char> indices_begin;        // len = num_reads + 1
			vector<index_t> variants_indices;           // len = num_variants

			// default constructor
			jralign_entry() : indices_begin(1, 0) {}
		};

		binary_file _file;
		bool        _include_variants{};
		unsigned    _min_reads{};
		int         _min_overhang{};
		error_handling _overhang_error{};
		long long _num_reads_loaded {};
		map<interval_t, jralign_entry> _juncs;

		// Order b_variant_t as key for fast lookup,
		// and an index as mapped_type keeping track of the order in which keys are added
		using variant_to_index = map<b_variant_t, index_t>;
		const genome_dna* _dna;
		variant_to_index _ra_variants;
		vector<b_variant_t> _ra_variants_dup;
		vector<index_t> _ra_variants_indices;

		bool _verbose{};
	};

private:
	void open_on_demand() const; // threadsafe version (const)

	mmap_file _fmap;
	string    _sourcefile;
	jraligns_table _juncs;
	variant_table  _variants;
	long long      _num_reads{};
};

bool operator<(const junction_read_alignments::builder::b_variant_t& lhs,
               const junction_read_alignments::builder::b_variant_t& rhs);

END_NAMESPACE_GK

#endif // __GENOME_KIT_JRALIGN_H__
