/*
Copyright (C) 2016-2023 Deep Genomics Inc. All Rights Reserved.
*/
#pragma once
#ifndef __GENOME_KIT_GENOME_ANNO_H__
#define __GENOME_KIT_GENOME_ANNO_H__

#include "table.h"
#include "file.h"
#include <vector>

BEGIN_NAMESPACE_GK
using std::vector;

using phase_t   = unsigned char;        		// Phase 0,1,2 of CDS elements
using tsl_t     = unsigned char;        		// GENCODE transcript support level (transcript_support_level=1,...,5,NA)
using elevel_t  = unsigned char;        		// GENCODE gene/transcript evidence level (level=1,2,3)
using biotype_t = unsigned char;        		// GENCODE gene/transcript biotype
using elem_num  = unsigned short;       		// Index of an element within its transcript
static constexpr tsl_t     num_tsl = 6;     	// [tsl1,tsl2,tsl3,tsl4,tsl5,tslNA] is 6 unique GENCODE identifiers mapped to [0..5] or 6=None
static constexpr elevel_t  num_elevel = 3;  	// [1,2,3] is 3 unique GENCODE identifiers mapped to [0..2] or 3=None
static constexpr biotype_t num_biotype = 78;	// 78 unique GENCODE identifiers mapped to [0..77]
static constexpr tsl_t     invalid_tsl = num_tsl;
static constexpr elevel_t  invalid_elevel = num_elevel;
static constexpr biotype_t invalid_biotype = num_biotype;

// fwd declare
class genome_anno;
class genome_t;
struct packed_gene;
struct packed_tran;
struct packed_exon;
struct packed_intr;
struct packed_cds;
struct packed_utr;
using tran_range = range_t<const packed_tran *>;
using exon_range = range_t<const packed_exon *>;
using intr_range = range_t<const packed_intr *>;
using cds_range = range_t<const packed_cds *>;
using utr_range = range_t<const packed_utr *>;

/////////////////////////////////////////////////////////////////
// genome_anno_table
/////////////////////////////////////////////////////////////////

// Basically an interval_table but with a reference to the parent genome_anno
// object so that a one table (e.g. genes) can contain indices identifying elements
// of another table (e.g. transcripts) and there's convenient way to access any of
// the other tables given only a reference to one of the tables.
template <typename T>
class genome_anno_table : public interval_table<T> {
public:
	genome_anno_table(genome_anno& anno): anno(anno) {}
	const genome_anno& anno;
	friend class genome_anno;
};

// This macro is used by genome annotation types to declare
// convenient 'unpacking' constructors on their unpacked row type.
// Specifically, these constructors facilitate easy conversion from
// packed_xxxx to an unpacked xxxx_t structure which
// contains more convenient form of access to the auxiliary members.
// Here xxxx can be gene, tran, exon, etc.
// Assumes naming convention of
//    xxxx_table      for the table type
//    packed_xxxx     for the table row type
//    xxxx_t          for the unpacked row type
#define GK_DECL_UNPACKED_ANNO_CTORS(name) \
	using packed_value = packed_##name; \
	using packed_table = genome_anno_table<packed_value>; \
	name##_t(const packed_value& src, const packed_table& table); \
	void unpack_from(const packed_value& src, const packed_table& table)

/////////////////////////////////////////////////////////////////
// genes
/////////////////////////////////////////////////////////////////

struct packed_gene: public interval_t {
	elevel_t   level;       // Evidence level
	biotype_t  type;        // Type of gene
	std::byte  pad{};
	index_t    num_trans;   // Number of transcripts in this gene
	index_t    tran0;       // Index into genome_anno::trans of first transcript that belong to this gene
	offset_t   aux;         // Offset into gene_table::auxpool for auxiliary data
};
// packed_gene memory layout: offset name (size+padding)
//  0 pos5         (4)
//  4 pos3         (4)
//  8 refg         (8)
// 16 chrom        (4)
// 20 strand       (1)
// 21 level        (1)
// 22 type         (1)
// 23 pad		   (1)
// 24 num_trans    (4)
// 28 tran0        (4)
// 32 aux          (4)
// 36 total
//
// Each packed_gene aux pool entry has a dynamic layout:
//    struct {
//        char id[n];       // n = strlen+1 (for NULL)
//        char name[n];     // n = strlen+1 (for NULL)
//    };

struct gene_t: public interval_t {
	elevel_t   level;       // Evidence level of gene
	biotype_t  type;        // Type of gene
	index_t    num_trans;   // Number of transcripts in this gene
	tran_range trans;       // Pointers into genome_anno::trans of transcripts belonging to this gene
	const char* id;         // ID of this gene
	const char* name;       // Name of this gene
	GK_DECL_UNPACKED_ANNO_CTORS(gene);
};

/////////////////////////////////////////////////////////////////
// transcripts
/////////////////////////////////////////////////////////////////

struct packed_tran : public interval_t {
	elevel_t level { invalid_elevel };        // Evidence level of transcript
	tsl_t    tsl {invalid_tsl};               // Transcript support level of this transcript
	biotype_t type { invalid_biotype };       // Type of transcript
	elem_num num_cdss { 0 };                  // Number of CDSs in this transcript
	elem_num num_utr5s { 0 };                 // Number of UTR5s in this transcript
	elem_num num_utr3s { 0 };                 // Number of UTR3s in this transcript
	elem_num num_exons { 0 };                 // Number of exons in this transcript
	index_t gene { invalid_index };  // Index into genome_anno::genes of the gene that this transcript belongs to
	index_t exon0 { invalid_index }; // Index into genome_anno::exons of the first exon that belongs to this transcript
	index_t intr0 { invalid_index }; // Index into genome_anno::intrs of the first intron that belongs to this transcript
	index_t cds0 { invalid_index }; // Index into genome_anno::cdss of the first CDS that belongs to this transcript
	index_t utr50 {invalid_index }; // Index into genome_anno::utr5s of the first UTR5 that belongs to this transcript
	index_t utr30 {invalid_index }; // Index into genome_anno::utr3s of the first UTR3 that belongs to this transcript
	offset_t aux { offset_t(-1) };  // Offset into tran_table::auxpool for auxiliary data
};

// packed_tran memory layout: offset name (size+padding)
//  0 pos5         (4)
//  4 pos3         (4)
//  8 refg         (8)
// 16 chrom        (4)
// 20 strand       (1)
// 21 level        (1)
// 22 tsl          (1)
// 23 type         (1)
// 24 num_cdss     (2)
// 26 num_utr5s    (2)
// 28 num_utr3s    (2)
// 30 num_exons    (2)
// 32 gene         (4)
// 36 exon0        (4)
// 40 intr0        (4)
// 44 cds0         (4)
// 48 utr50        (4)
// 52 utr30        (4)
// 56 aux          (4)
// 60 total
//
// Each packed_tran aux pool entry has a dynamic layout:
//    struct {
//        char id[n];          // n = strlen+1 (for NULL)
//        char ccds_id[n];     // n = strlen+1 (for NULL)
//        char protein_id[n];  // n = strlen+1 (for NULL)
//        char product[n];     // n = strlen+1 (for NULL)
//    };

struct tran_t: public interval_t {
	elevel_t   level;       // Evidence level of transcript
	biotype_t  type;        // Type of parent gene
	tsl_t      tsl;         // Transcript support level of this transcript
	elem_num   num_exons;   // Number of exons in this transcript
	elem_num   num_cdss;    // Number of CDSs in this transcript
	elem_num   num_utr5s;   // Number of UTR5s in this transcript
	elem_num   num_utr3s;   // Number of UTR3s in this transcript
	const packed_gene* gene;// Pointer  into genome_anno::genes of the gene that this transcript belongs to
	exon_range exons;       // Pointers into genome_anno::exons of exons belonging to this transcript
	intr_range intrs;       // Pointers into genome_anno::intrs of introns belonging to this transcript
	cds_range  cdss;        // Pointers into genome_anno::cdss of CDSs belonging to this transcript
	utr_range utr5s;       // Pointers into genome_anno::utr5s of UTR5s belonging to this transcript
	utr_range utr3s;       // Pointers into genome_anno::utr3s of UTR3s belonging to this transcript
	const char* id;         // ID of this transcript
	const char* ccds_id;    // CCDSID of this transcript, if applicable
	const char* protein_id; // Protein of this transcript, if applicable
	const char* product;    // Description of transcript's product, if applicable
	GK_DECL_UNPACKED_ANNO_CTORS(tran);
};

/////////////////////////////////////////////////////////////////
// exons
/////////////////////////////////////////////////////////////////

struct packed_exon: public interval_t {
	uint8_t  is_last;          // Last in parent transcripts's list of exons?
	elem_num index;            // The number of this exon within its transcript
	index_t  tran;             // Index into genome_anno::trans of the transcript that this exon belongs to
	index_t  cds;              // Index into genome_anno::cdss of the CDS that belongs to this exon, if any
	index_t  utr5;             // Index into genome_anno::utr5s of the UTR5 that belongs to this exon, if any
	index_t  utr3;             // Index into genome_anno::utr3s of the UTR3 that belongs to this exon, if any
	offset_t id;               // Offset into exon_table::auxpool for this exon's ID, which may be shared by other exons
};
// packed_exon memory layout: offset name (size+padding)
//  0 pos5         (4)
//  4 pos3         (4)
//  8 refg         (8)
// 16 chrom        (4)
// 20 strand       (1)
// 21 is_last      (1)
// 22 index        (2)
// 24 tran         (4)
// 28 cds          (4)
// 32 utr5         (4)
// 36 utr3         (4)
// 40 id           (4)
// 44 total
//
// Each packed_exon aux pool entry has a dynamic layout:
//    struct {
//        char id[n];       // n = strlen+1 (for NULL)
//    };

struct exon_t: public interval_t {
	elem_num   index;       // The number of this exon within its transcript
	elevel_t   level;       // Evidence level of parent transcript
	biotype_t  type;        // Type of parent gene
	tsl_t      tsl;         // Transcript support level of the parent transcript
	const packed_tran* tran;// Pointer into genome_anno::trans of the transcript that this exon belongs to
	const packed_cds*  cds; // Pointer into genome_anno::cdss of the CDS that belongs to this exon, if any
	const packed_utr* utr5;// Pointer into genome_anno::utr5s of the UTR5 that belongs to this exon, if any
	const packed_utr* utr3;// Pointer into genome_anno::utr3s of the UTR3 that belongs to this exon, if any
	const char* id;         // ID of this exon
	GK_DECL_UNPACKED_ANNO_CTORS(exon);
};

INLINE const packed_exon* get_prev(const packed_exon* exon) { return exon->index==0 ? nullptr : exon-1; }
INLINE const packed_exon* get_next(const packed_exon* exon) { return exon->is_last ? nullptr : exon+1; }
       const packed_intr* get_prev_intr(const packed_exon* exon, const genome_anno& anno);
       const packed_intr* get_next_intr(const packed_exon* exon, const genome_anno& anno);
       const packed_intr* get_prev_intr(const packed_exon* exon, const genome_t& genome);
       const packed_intr* get_next_intr(const packed_exon* exon, const genome_t& genome);


/////////////////////////////////////////////////////////////////
// introns
/////////////////////////////////////////////////////////////////

struct packed_intr: public interval_t {
	uint8_t  is_last;          // Last in parent transcripts's list of introns?
	elem_num index;            // The number of this intron within its transcript
	index_t  tran;             // Index into genome_anno::trans of the transcript that this intr belongs to
};
// packed_intr memory layout: offset name (size+padding)
//  0 pos5         (4)
//  4 pos3         (4)
//  8 refg         (8)
// 16 chrom        (4)
// 20 strand       (1)
// 21 is_last      (1)
// 22 index        (2)
// 24 tran         (4)
// 28 total

struct intr_t: public interval_t {
	elem_num   index;       // The number of this intron within its transcript
	const packed_tran* tran;// Pointer into genome_anno::trans of the transcript that this intron belongs to
	elevel_t   level;       // Evidence level of parent transcript
	biotype_t  type;        // Type of parent gene
	tsl_t      tsl;         // Transcript support level of the parent transcript
	GK_DECL_UNPACKED_ANNO_CTORS(intr);
};

INLINE const packed_intr* get_prev(const packed_intr* intr) { return intr->index==0 ? nullptr : intr-1; }
INLINE const packed_intr* get_next(const packed_intr* intr) { return intr->is_last ? nullptr : intr+1; }
       const packed_exon* get_prev_exon(const packed_intr* intr, const genome_anno& anno);
       const packed_exon* get_next_exon(const packed_intr* intr, const genome_anno& anno);
       const packed_exon* get_prev_exon(const packed_intr* intr, const genome_t& genome);
       const packed_exon* get_next_exon(const packed_intr* intr, const genome_t& genome);

/////////////////////////////////////////////////////////////////
// cds
/////////////////////////////////////////////////////////////////

struct packed_exonic_element: public interval_t {
	uint8_t is_first: 1;  // First in parent transcripts's list of CDS elements?
	uint8_t is_last: 1;   // Last in parent transcripts's list of CDS elements?
	uint8_t phase: 6;     // Phase 0, 1, or 2 modulo reading frame, only used for CDS
	elem_num exon_index;  // The 0-based index of the exon that this element belongs to within its parent transcript (overflows at 4096, but human genes top out at <400 exons)
	index_t  tran;        // Index into genome_anno::trans of the transcript that this element belongs to
};

struct packed_cds: public packed_exonic_element {};

// packed_cds/packed_utr memory layout: offset name (size+padding)
//  0 pos5             			(4)
//  4 pos3             			(4)
//  8 refg         				(8)
// 16 chrom        				(4)
// 20 strand       				(1)
// 21 is_first/is_last/phase	(1)
// 22 exon_index 				(2)
// 24 tran             			(4)
// 28 total

struct cds_t: public interval_t {
	phase_t    phase;       // 0,1,2
	tsl_t      tsl;         // Transcript support level of the parent transcript
	elevel_t   level;       // Evidence level of parent transcript
	const packed_tran* tran;// Pointer into genome_anno::trans of the transcript that this CDS belongs to
	const packed_exon* exon;// Pointer into genome_anno::exons of the exon that this CDS belongs to
	GK_DECL_UNPACKED_ANNO_CTORS(cds);
};

INLINE const packed_cds* get_prev(const packed_cds* cds) { return cds->is_first ? nullptr : cds-1; }
INLINE const packed_cds* get_next(const packed_cds* cds) { return cds->is_last  ? nullptr : cds+1; }

/////////////////////////////////////////////////////////////////
// utr5/3
/////////////////////////////////////////////////////////////////

struct packed_utr : public packed_exonic_element {};

struct utr_t: public interval_t {
	const packed_tran* tran; /* Pointer into genome_anno::trans of the transcript that this UTR belongs to */
	const packed_exon* exon; /* Pointer into genome_anno::exons of the exon that this UTR belongs to */
	GK_DECL_UNPACKED_ANNO_CTORS(utr);
};
INLINE const packed_utr* get_prev(const packed_utr* utr) { return utr->is_first ? nullptr : utr - 1; }
INLINE const packed_utr* get_next(const packed_utr* utr) { return utr->is_last ? nullptr : utr + 1; }

using gene_table = genome_anno_table<packed_gene>;
using tran_table = genome_anno_table<packed_tran>;
using exon_table = genome_anno_table<packed_exon>;
using intr_table = genome_anno_table<packed_intr>;
using cds_table = genome_anno_table<packed_cds>;
using utr_table = genome_anno_table<packed_utr>;

/////////////////////////////////////////////////////////////////
// genome_anno
/////////////////////////////////////////////////////////////////

class genome_anno {
public:
	genome_anno();

	const gene_table& genes() const;
	const tran_table& trans() const;
	const exon_table& exons() const;
	const intr_table& intrs() const;
	const cds_table&  cdss() const;
	const utr_table&  utr5s() const;
	const utr_table&  utr3s() const;

	bool empty() const noexcept;

	// Sets the original source file (e.g. GFF3 file), but does
	// not actually open it until open() is called.
	void set_source(string sourcefile);
	void open();
	void close();
	INLINE bool is_open() const { return _fmap.is_open(); }
	INLINE void ensure_open() const { if (!is_open()) open_on_demand(); }  // Fast check if file / indices already loaded

	INLINE const string& source() const { return _sourcefile; }
	static int  binary_version();
	static std::vector<std::string> build_gencode(const char* infile, const char* outfile, const genome_t& genome);
	static std::vector<std::string> build_ucsc_refseq(const char* ucsc_db_dir, const char* outfile, const genome_t& genome);
	static std::vector<std::string> build_ncbi_refseq(const char* infile, const char* outfile, const genome_t& genome);

private:
	void open_on_demand() const; // threadsafe version (const)
	class builder;

	gene_table _genes;
	tran_table _trans;
	exon_table _exons;
	intr_table _intrs;
	cds_table  _cdss;
	utr_table  _utr5s;
	utr_table  _utr3s;
	mmap_file _fmap;
	string _sourcefile;
};

//////////////////////////////////////////////////////////////////

// Operator== is needed to ensure interval_t comparisons aren't implicitly used to compare annotation objects.
// Packed annotation objects are unique based on either their aux offset or their index+parent
INLINE bool operator==(const packed_gene& x, const packed_gene& y) { return x.aux == y.aux; }
INLINE bool operator==(const packed_tran& x, const packed_tran& y) { return x.aux == y.aux; }
INLINE bool operator==(const packed_exon& x, const packed_exon& y) { return x.tran == y.tran && x.index == y.index; }
INLINE bool operator==(const packed_intr& x, const packed_intr& y) { return x.tran == y.tran && x.index == y.index; }
INLINE bool operator==(const packed_cds&  x, const packed_cds&  y) { return x.tran == y.tran && x.exon_index == y.exon_index; }
INLINE bool operator!=(const packed_gene& x, const packed_gene& y) { return !(x == y); }
INLINE bool operator!=(const packed_tran& x, const packed_tran& y) { return !(x == y); }
INLINE bool operator!=(const packed_exon& x, const packed_exon& y) { return !(x == y); }
INLINE bool operator!=(const packed_intr& x, const packed_intr& y) { return !(x == y); }
INLINE bool operator!=(const packed_cds&  x, const packed_cds&  y) { return !(x == y); }

INLINE size_t hash(const packed_gene& g) { return g.aux; }
INLINE size_t hash(const packed_tran& t) { return t.aux; }
INLINE size_t hash(const packed_exon& e)
{
	size_t x = 23;
	x ^= ((size_t)e.tran+23)*102564;
	x ^= ((size_t)e.index+562123)*73;
	return x;
}
INLINE size_t hash(const packed_intr& e)
{
	size_t x = 812997;
	x ^= ((size_t)e.tran+23)*102564;
	x ^= ((size_t)e.index+562123)*73;
	return x;
}

INLINE size_t hash(const packed_exonic_element& e)
{
	size_t x = 812997;
	x ^= ((size_t)e.tran+23)*102564;
	x ^= ((size_t)e.exon_index+562751)*91;
	return x;
} // TODO: replace all these with standard byte hashing algorithm

const char* biotype_as_cstr(biotype_t type);

string default_anno_sourcefile(string_view config, string_view data_dir = default_data_directory);

END_NAMESPACE_GK

#endif // __GENOME_KIT_GENOME_ANNO_H__
