/*
Copyright (C) 2016-2023 Deep Genomics Inc. All Rights Reserved.
*/
#pragma once
#ifndef __GENOME_KIT_JRDIST_H__
#define __GENOME_KIT_JRDIST_H__

#include "jralign.h"
#include <functional>
#include <map>
#include <string>
#include <vector>

BEGIN_NAMESPACE_GK

using std::vector;
using std::map;
using std::string;

class jrdist_t;
class jrdist_table;

/////////////////////////////////////////////////////////////////
// junction read distribution
/////////////////////////////////////////////////////////////////

struct packed_jrdist: public interval_t {
	std::byte      pad;
	unsigned short num_counts; // Total number of read counts on this junction
	unsigned       num_reads;  // Total number of reads counted (sum of all read counts)
	offset_t  aux;             // Offset into junc_readdist_table::auxpool for auxiliary data
};
// packed_jrdist memory layout: offset name (size+padding)
//  0 pos5         (4)
//  4 pos3         (4)
//  8 refg         (8)
// 16 chrom        (4)
// 20 strand       (1)
// 21 pad          (1)
// 22 num_counts   (2)
// 24 num_reads    (4)
// 28 aux          (4)
// 32 total
//
// Each packed_jrdist aux pool entry has a dynamic layout:
//
//    struct {
//        unsigned char count_size;          // size of each count (1, 2, or 4 bytes)
//        ...                                // enough padding bytes to align counts to count_size bytes; padding bytes filled with count_size
//        uint8/16/32   counts[n];           // [i] = read count for count number i; size depends on count_size
//        unsigned char ushifts[n];          // [i] = shift for count entry i
//        unsigned char flags[(2*n+3)/4];    // bit [2*i] = strand of count entry i, and bit [2*i+1] = sign of shift for count entry i
//    };
//

struct jrcount_t {
	unsigned count;   // Count at this position; whether it's starts-only or starts-or-ends depends on the function that produced it
	pos_t    shift;   // Shift relative to junction pos5 or pos3 for this count. Always non-zero (-1 = last upstream position, +1 first downstream position)
	strand_t strand;  // Strand the read mapped to (pos_strand for reference strand, neg_strand for reverse strand)
	INLINE jrcount_t(unsigned count, pos_t shift, strand_t strand): count(count), shift(shift), strand(strand) { }
	INLINE bool is_start() const { return (shift < 0) == (strand == pos_strand); }
	INLINE bool is_end()   const { return (shift > 0) == (strand == pos_strand); }
};

inline bool operator==(const jrcount_t& a, const jrcount_t& b) { return a.count == b.count && a.shift == b.shift && a.strand == b.strand; }
inline bool operator!=(const jrcount_t& a, const jrcount_t& b) { return !(a == b); }

class jrdist_t: public interval_t {
public:
	INLINE int      num_counts() const { return _num_counts; }
	INLINE unsigned num_reads()  const { return _num_reads; }
	const jrcount_t operator[](int i) const;

	jrdist_t(index_t              src, const jrdist_table& table);
	jrdist_t(const packed_jrdist& src, const jrdist_table& table);

private:
	int                  _num_counts;
	unsigned             _num_reads;
	const unsigned char* _counts;  // Not necessarily stored as char, could be short or int
	const unsigned char* _ushifts;
	const unsigned char* _flags;

	void unpack_from(const packed_jrdist& src, const jrdist_table& table);
};

/////////////////////////////////////////////////////////////////
// junction read distribution table
/////////////////////////////////////////////////////////////////

class jrdist_table: public interval_table<packed_jrdist> { };

/////////////////////////////////////////////////////////////////
// read distribution table (junction reads + body reads)
/////////////////////////////////////////////////////////////////

class read_distributions {
public:
	const jrdist_table& juncs() const;

	// Sets the original source file (e.g. GFF3 file), but does
	// not actually open it until open() is called.
	void set_source(string sourcefile);
	void open();
	void close();
	INLINE bool is_open() const { return _fmap.is_open(); }
	INLINE void ensure_open() const { if (!is_open()) open_on_demand(); }  // Fast check if file / indices already loaded
	INLINE const string& source() const { return _sourcefile; }
	INLINE long long num_reads() const { ensure_open(); return _num_reads; }

	static int rdist_version();

	// Takes JRALIGN files as input, and outputs an JRDIST file.
	class builder {
	public:
		using filterfn
			= std::function<bool(const junction_read_alignments&, int junction, unsigned int read)>;

		builder(const char* outfile);
		INLINE void set_min_reads(int value) { _min_reads = (unsigned)value; }
		INLINE void set_min_overhang(int value) { _min_overhang = value; }
		void exclude(const interval_t& interval);
		void allow(const interval_t& interval);
		void filter(filterfn filter) { _filter = filter; }
		void add(const char* infile);
		void finalize();

	private:

		struct jrdist_collector;

		class jrdist_entry { NOCOPY(jrdist_entry)
		public:
			jrdist_entry() = default;
			~jrdist_entry();

			unsigned char mode()       const { return _mode; }
			unsigned char count_size() const { GK_ASSERT(_mode != list_mode); return 1 << (_mode - count8_mode); }
			void add(const jralign_t read);
			void collect(jrdist_collector& collector);
			void clear();

		private:
			// pair<interval_t, jrdist_entry> will be 12+max_list_size*3 bytes, and we want that size
			// to be a multiple of 8 since std::map node alignment is 8; so max_list_size = 12 is a good number
			enum { max_list_size = 12 };  // must be >= 4 based on assumptions made by add(), and <= 255 due to _size type
			enum { list_mode, count8_mode, count16_mode, count32_mode, num_modes };

			unsigned char _mode{list_mode}; // 0 (as_list mode), 1 (byte count), 2 (short count), or 4 (unsigned count); count size upgraded on demand
			unsigned char _size{};          // If mode == 0, size of list; otherwise size of each direction in the counts array (-_size+1 and _size-1 are valid)
			jralign_t _list[max_list_size]; // When mode == 0, explicitly stored read alignments

			// The _counts member is stored in memory overlapping _list, because these two modes (count and list)
			// are mutually exclusive, and because it allows us to store pair<interval_t, jrdist_entry> in a
			// more compact struct layout. The offset of 6 ensures the _counts pointer is aligned to 8 bytes,
			// but since it sits in the middle of the _list array, special care must be taken in ensure_counts_mode.
			INLINE const void*& _counts() const { return *(const void**)((char*)_list+4); }
			INLINE       void*& _counts()       { return *(      void**)((char*)_list+4); }

			template <typename T> void inc_count(unsigned char overhang, strand_t strand, int direction);
			template <typename T> void collect_impl(jrdist_collector& collector);
			void ensure_counts_mode();

			static void (jrdist_entry::*inc_count_fn[num_modes])(unsigned char overhang, strand_t strand, int direction);
		};

		binary_file _file;
		std::optional<refg_t> _refg;
		unsigned    _min_reads{1};
		int         _min_overhang{1};
		vector<interval_t> _exclude;
		vector<interval_t> _allow;
		map<interval_t, jrdist_entry> _juncs;
		filterfn _filter{};

		enum class stranded { unknown, no, yes } _stranded{};
	};

private:
	void open_on_demand() const; // threadsafe version (const)

	mmap_file _fmap;
	string    _sourcefile;
	jrdist_table _juncs;
	long long    _num_reads{};
};

END_NAMESPACE_GK

#endif // __GENOME_KIT_JRDIST_H__
