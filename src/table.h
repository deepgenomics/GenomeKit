/*
Copyright (C) 2016-2023 Deep Genomics Inc. All Rights Reserved.
*/
#pragma once
#ifndef __GENOME_KIT_TABLE_H__
#define __GENOME_KIT_TABLE_H__

#include "array_view.h"
#include "file.h"
#include "interval.h"
#include "util.h"
#include <algorithm>
#include <iterator>
#include <vector>

BEGIN_NAMESPACE_GK
using std::vector;
using std::upper_bound;
using std::lower_bound;
using std::max;

using index_t    = int;
using indices_t  = vector<index_t>;
using index_iter = const index_t*;

const index_t invalid_index = (index_t)0x80000000; // largest negative int32; increases likelihood of crash if dereferenced

/////////////////////////////////////////////////////////////////
// Macros to help define comparison functors
/////////////////////////////////////////////////////////////////

// TODO: use get_pos3/get_pos5 etc to turn these macros instantiations into a single template
// Compare an element's field to a specific value
#define DEFN_COMPARE_FIELD_TO_VALUE(field, type) \
	template <typename T> \
	struct compare_##field##_to_value { \
		struct lower { \
			INLINE lower(const T* items): items(items) { } \
			INLINE bool operator()(index_t i, type field) const { return items[i].field < field; } \
			const T* const items; \
		}; \
		struct upper { \
			INLINE upper(const T* items): items(items) { } \
			INLINE bool operator()(type field, index_t i) const { return field < items[i].field; } \
			const T* const items; \
		}; \
	}

// Compare two elements by a specific field
#define DEFN_COMPARE_FIELDS(field) \
	template <typename T> \
	struct compare_##field { \
		INLINE compare_##field(const T* items): items(items) { } \
		INLINE bool operator()(index_t a, index_t b) const { return items[a].field < items[b].field; } \
		const T* const items; \
	}

template <typename T>
struct compare_pos5_to_value {
	struct lower {
		INLINE lower(const T* items): items(items) { }
		INLINE bool operator()(index_t i, pos_t pos5) const { return items[i].pos5 < pos5; }
		const T* const items;
	};
	struct upper {
		INLINE upper(const T* items): items(items) { }
		INLINE bool operator()(pos_t pos5, index_t i) const {
			return pos5 < items[i].pos5;
		}
		const T* const items;
	};
};

//DEFN_COMPARE_FIELD_TO_VALUE(pos5, pos_t);
DEFN_COMPARE_FIELD_TO_VALUE(pos3, pos_t);

DEFN_COMPARE_FIELDS(pos5);
DEFN_COMPARE_FIELDS(pos3);

///////////////////////////////////////////////////////////////////////////
// table<T,I> - interface/implementation details common to tables
//              with compact, memory mapped, indexed representations
///////////////////////////////////////////////////////////////////////////

template <typename T, typename I>
class table { NOCOPY(table)
public:
	using value_t = T;

	table() = default;

	GK_DECLARE_ARRAY_TYPES(T, int)
	GK_DECLARE_ARRAY_CONST_METHODS(_elems, (int)_elems.size())
	GK_DECLARE_ARRAY_CONST_ITERATOR(_elems.data())

	INLINE int index_of(const T& elem) const {
		std::ptrdiff_t index = &elem - _elems.data();
		GK_CHECK(index >= 0 && (size_t)index < _elems.size(), value, "Table does not contain the given element");
		return int_cast<int>(index);
	}

	const char* aux() const { return _aux.data(); }
	void load(mmap_file& file);
	bool valid() const { return _loaded_mmap != nullptr && _loaded_mmap->is_open(); }
	bool stranded() const { return _stranded; }

	// cursor_t
	//   Cursor for iterating over results of a find operation that returns
	//   a contiguous range of indices, such as find, find_5p_within, find_3p_within.
	//   A cursor will iterate over those indices and dereference the index.
	//
	class cursor_t {
	public:
		index_iter cur;  // Pointer to index of current element
		const T* elems;  // Pointer to actual elements; needed to dereference the cur index
		INLINE cursor_t(index_iter cur, const T* elems): cur(cur), elems(elems) { }
		INLINE cursor_t& operator++() { ++cur; return *this; }
		INLINE cursor_t operator++(int) { cursor_t c = *this; ++(*this); return c; }
		const T& operator*()  const { return *(elems + *cur); }
		const T* operator->() const { return  (elems + *cur); }
		INLINE bool operator==(const cursor_t& rhs) const { GK_DBASSERT(elems == rhs.elems); return cur == rhs.cur; }
		INLINE bool operator!=(const cursor_t& rhs) const { return !(*this == rhs); }
	};

	// cursor_range
	//   A begin/end range to iterate over using cursor_t.
	//   Makes for nice C++11-style for-loops over the results.
	//
	class cursor_range {
	public:
		using cursor_t = typename table::cursor_t;
		index_iter first;
		index_iter last;
		const T* elems;
		INLINE cursor_range(index_iter begin, index_iter end, const T* elems): first(begin), last(end), elems(elems) { }
		INLINE cursor_t begin() const { return cursor_t(first, elems); }
		INLINE cursor_t end()   const { return cursor_t(last,  elems); }
		INLINE size_t   size()  const { return (size_t)(last - first); }
	};

	// table<T,I>::builder
	//    Utility class for building a table.
	//
	// Typical use case is as follows:
	//   - add each row one by one, possibly adding auxiliary
	//     data for each row.
	//   - build index over the rows
	//   - dump the structure to disk so that it can
	//     be easily mapped by load() later on.
	//
	class builder {
	public:
		builder(bool stranded=true): _stranded(stranded) { }
		GK_DECLARE_ARRAY_TYPES(T, int)
		GK_DECLARE_ARRAY_METHODS(_elems, (int)_elems.size())
		GK_DECLARE_ARRAY_CONST_ITERATOR(&_elems[0])

		index_t add_elem();
		index_t add_elem(const T& elem);

		template <typename X>
		offset64_t add_aux(const X* first, const X* last);
		template <typename X>
		offset64_t add_aux(const vector<X>& v);
		offset64_t add_aux(std::string_view str);
		void       align_aux(int bytes); // add some arbitrary aux bytes until offset is multiple of 'bytes'
		INLINE offset64_t curr_aux() { return (offset64_t)_aux.size(); }

		void dump(binary_file& out) const;

	private:
		vector<T>     _elems;
		vector<char>  _aux;
		bool          _stranded;
	};

protected:
	friend class builder;
	using indices_t = I;

	// The three inputs are as follows:
	//   1) A query interval 'c'.
	//   2) A comparison functor template 'class P<T>' defined by
	//      DEFN_COMPARE_FIELD_TO_VALUE, which is responsible for
	//      comparing the field of interest to a particular value
	//      (e.g. a position from the query interval); the particular
	//      field, and therefor functor, is known at compile time.
	//   3) An array of indices 'idx' into elems[...] that is pre-sorted
	//      according to the field of interest, to facilitate binary
	//      search.
	// The result is a range (a,b) where iterating over that range
	// returns all the elements that fall within the query interval.
	//
	template <template <typename> class P>
	cursor_range find_by_field(const interval_t& c, const array_view<index_t>& idx) const;

	// Data members
	array_view<T>    _elems;
	array_view<char> _aux;
	mmap_file*       _loaded_mmap{};
	indices_t        _idx;
	bool             _stranded{true};
};

/////////////////////////////////////////////////////////////////
// interval_table
/////////////////////////////////////////////////////////////////

// TODO change the way array_view works so that even the size is
// embedded in the memory mapped address space, and nothing is explicitly
// loaded into local structs of table, coord_idx, or interval_idx.
// This doesn't change the file format, but does change the struct
// definitions on the client side.

template <typename T>
class interval_idx { NOCOPY(interval_idx)
public:
	interval_idx() = default;

	pos_t                      max_interval_size(chrom_t chrom, strand_t strand) const;
	const array_view<index_t>& by_pos5(chrom_t chrom, strand_t strand) const;
	const array_view<index_t>& by_pos3(chrom_t chrom, strand_t strand) const;

	static void dump(binary_file& out, const vector<T>& rows);
	void load(mmap_file& in);

private:
	struct sorted_indices {
		pos_t               max_interval_size{};
		array_view<index_t> by_pos5;
		array_view<index_t> by_pos3;
	};

	const sorted_indices& get_indices(chrom_t chrom, strand_t strand) const;

	chrom_map_t<strand_t, sorted_indices> indices_by_stranded_key;
};

template <typename T>
class interval_table: public table<T, interval_idx<T> > {
	using base_class = table<T, interval_idx<T> >;
public:
	using cursor_range = typename base_class::cursor_range;

	// filtered_cursor<P>
	//   Iterates over results of a find operation that needs certain elements
	//   skipped according to predicate P, such as find_within or find_overlapping.
	//   Currently implemented in lieu of a real interval_tree structure.
	//
	template <typename P>
	class filtered_cursor {
	public:
		index_iter cur;  // Pointer to index of current interval; valid if < last
		index_iter last; // Terminator that this iterator must not walk past; alternative design is to add sentinel to elems instead
		const T* elems;  // Pointer to actual elements; needed to filter intervals based on pos5, and to dereference the cur index
		pos_t pos5;      // The 5p end of the original query interval; needed to filter intervals
		INLINE filtered_cursor(index_iter cur, index_iter end, const T* elems, pos_t pos5): cur(cur), last(end), elems(elems), pos5(pos5) { if (cur < last) advance_until_valid(); }
		INLINE filtered_cursor& operator++() { GK_DBASSERT(cur < last); if (++cur < last) advance_until_valid(); return *this; }
		INLINE filtered_cursor operator++(int) { filtered_cursor c = *this; ++(*this); return c; }
		const T& operator*()  const { GK_DBASSERT(cur < last); return *(elems + *cur); }
		const T* operator->() const { GK_DBASSERT(cur < last); return  (elems + *cur); }
		INLINE bool operator==(const filtered_cursor& rhs) const { GK_DBASSERT(elems == rhs.elems); return cur == rhs.cur; }
		INLINE bool operator!=(const filtered_cursor& rhs) const { return !(*this == rhs); }
	private:
		void advance_until_valid(); // If the current position is not in the query interval, advance
	};

	// filtered_cursor_range<P>
	//   A begin/end range to iterate over using filtered_cursor<P>.
	//   Makes for nice C++11-style for-loops over the results.
	//
	template <typename P>
	class filtered_cursor_range {
	public:
		using cursor_t = filtered_cursor<P>;
		index_iter first;  // Pointer to index of the first interval; valid if < last
		index_iter last;   // Terminator that this iterator must not walk past; alternative design is to add sentinel to elems instead
		const T* elems;    // Pointer to actual elements; needed to filter intervals based on pos5, and to dereference the cur index
		pos_t pos5;
		INLINE filtered_cursor_range(index_iter first, index_iter last, const T* elems, pos_t pos5): first(first), last(last), elems(elems), pos5(pos5) { }
		INLINE filtered_cursor_range(typename interval_table::cursor_range range, pos_t pos5): first(range.first), last(range.last), elems(range.elems), pos5(pos5) { }
		INLINE cursor_t begin() const { return cursor_t(first, last, elems, pos5); }
		INLINE cursor_t end()   const { return cursor_t(last,  last, elems, pos5); }
	};

	// When filtering a "within" query, must compare the candidate interval 5p end with
	// the query interval 5p end. When filtering an "overlap" query, must compare the
	// candidate interval 3p end with the query interval 5p end. This is a result of
	// doing interval queries with simple binary search + filtering, rather than with
	// an interval_tree data structure.
	using cursor_within            = filtered_cursor<get_pos5<T>>;
	using cursor_overlapping       = filtered_cursor<get_pos3<T>>;
	using cursor_range_within      = filtered_cursor_range<get_pos5<T>>;
	using cursor_range_overlapping = filtered_cursor_range<get_pos3<T>>;

	cursor_range_within find_within(const interval_t& i) const;
	cursor_range_overlapping find_overlapping(const interval_t& i) const;
	cursor_range find_5p_within(const interval_t& i) const;
	cursor_range find_3p_within(const interval_t& i) const;
	cursor_range find_5p_aligned(const interval_t& i) const;
	cursor_range find_3p_aligned(const interval_t& i) const;
	cursor_range find_exact(const interval_t& i) const;

private:
	// Private versions of find_5/3p_within don't assert !_single_stranded, and are used as subroutines by the other find_ methods.
	cursor_range _find_5p_within(const interval_t& i) const;
	cursor_range _find_3p_within(const interval_t& i) const;
};

////////////////////////////////////////////////////////
////////////////////////////////////////////////////////
////////////////////////////////////////////////////////

template <typename T, typename I>
void table<T,I>::load(mmap_file& in) {
	GK_ASSERT(_loaded_mmap == nullptr); // disallow multiple loads to simplify validation of mapped memory
	_loaded_mmap = &in;

	// Read the rows themselves
	in.read_until_align32();
	in.read_checkpoint(0x85420001);
	_elems.read(in);

	// Read the auxiliary memory pool
	in.read_checkpoint(0x85420002);
	_aux.read(in);

	// Read the indices
	_stranded = in.read<uint8_t>() != 0;
	in.read_until_align32();
	in.read_checkpoint(0x85420003);
	_idx.load(in);

	// Final checkpoint to make sure everything's ok
	in.read_checkpoint(0x85420990);
}

template <typename T, typename I> template <template <typename> class P>
typename table<T,I>::cursor_range table<T,I>::find_by_field(const interval_t& i, const array_view<index_t>& idx) const
{
	GK_CHECK(this->valid(), file, "IntervalTable has been invalidated by close or with statement on source object.");
	if (!_elems.empty()) {
		// TODO: data_dir injected as a context
		const auto& refgs = get_refg_registry();
		GK_CHECK(i.refg == _elems[0].refg, value, "Reference genome '{}' does not match query coordinate '{}'.",
				 refgs.refg_as_sv(_elems[0].refg), refgs.refg_as_sv(i.refg));
	}
	pos_t a = i.start();
	pos_t b = i.end() - 1;

	// Binary search using ordering 'idx' to identify rows that fall within query interval,
	// where P<T>::lower will test "elems[i].field < field" for each index i encountered during binary
	// search, and P<T>::upper will compare "field < elems[i].field".
	// Here, 'field' can be pos, pos5, pos3, etc, where 'idx' should be ordered according to that field.
	index_iter first = idx.begin(), last = idx.end();
	first = lower_bound(first, last, a, typename P<T>::lower(_elems.data())); // First row >= a, or b
	last  = upper_bound(first, last, b, typename P<T>::upper(_elems.data())); // First row >  b, or b
	return cursor_range(first, last, _elems.data());
}

////////////////////////////////////////////////////////

template <typename T, typename I>
index_t table<T,I>::builder::add_elem()
{
	auto i = int_cast<index_t>(_elems.size());
	_elems.resize(i+1);
	return i;
}

template <typename T, typename I>
index_t table<T,I>::builder::add_elem(const T& row)
{
	auto i = int_cast<index_t>(_elems.size());
	_elems.push_back(row);
	return i;
}

template <typename T, typename I> template <typename X>
offset64_t table<T,I>::builder::add_aux(const X* first, const X* last)
{
	offset64_t offset = curr_aux();
	_aux.insert(_aux.end(), (char*)first, (char*)last); // Copy the element indices into an array
	return offset;
}

template <typename T, typename I> template <typename X>
offset64_t table<T,I>::builder::add_aux(const vector<X>& v)
{
	if (v.empty())
		return curr_aux();
	return add_aux(&v[0], &v[0] + v.size());
}

template <typename T, typename I>
offset64_t table<T,I>::builder::add_aux(std::string_view str)
{
	offset64_t offset = curr_aux();
	// Append the string plus NULL terminator to the array
	_aux.insert(_aux.end(), std::cbegin(str), std::cend(str));
	_aux.push_back(0);
	return offset;
}

template <typename T, typename I>
void table<T,I>::builder::align_aux(int bytes)
{
	offset64_t offset = curr_aux();
	while (offset++ % bytes)
		_aux.push_back(0);
}

template <typename T, typename I>
void table<T,I>::builder::dump(binary_file& out) const {

	// Write the rows themselves, specifically all fixed-width data not stored in aux.
	// IMPORTANT: This kind of simple serialization assumes T is POD type
	//            with all data either internal or stored with relative offsets
	//            -- members in type T cannot be pointers, as the stored address
	//            will be invalid the next time they're read in.
	out.write_until_align(4);
	out.write_checkpoint(0x85420001);
	out.write_array(_elems);

	// Write auxiliary memory pool
	out.write_checkpoint(0x85420002);
	out.write_array(_aux);

	// Write the indices for this table
	out.write<uint8_t>(_stranded ? 1 : 0);
	out.write_until_align(4);
	out.write_checkpoint(0x85420003);
	I::dump(out, _elems);

	// Final checkpoint to make sure everything's ok
	out.write_checkpoint(0x85420990);
}

////////////////////////////////////////////////////////
////////////////////////////////////////////////////////
////////////////////////////////////////////////////////

template <typename T>
pos_t interval_idx<T>::max_interval_size(chrom_t chrom, strand_t strand) const
{
	return get_indices(chrom, strand).max_interval_size;
}
template <typename T>
const array_view<index_t>& interval_idx<T>::by_pos5(chrom_t chrom, strand_t strand) const
{
	return get_indices(chrom, strand).by_pos5;
}
template <typename T>
const array_view<index_t>& interval_idx<T>::by_pos3(chrom_t chrom, strand_t strand) const
{
	return get_indices(chrom, strand).by_pos3;
}

template <typename T>
void interval_idx<T>::dump(binary_file& out, const vector<T>& rows)
{
	struct to_dump {
		pos_t           max_interval_size{};
		vector<index_t> by_pos3;
	};
	chrom_map_t<strand_t, to_dump> indices;

	// Add index of each element to the corresponding indices array
	for (index_t i = 0; i < (index_t)rows.size(); ++i) {
		const T& r = rows[i];
		auto& data = indices[{r.chrom, r.strand}];
		data.by_pos3.push_back(i);
		data.max_interval_size = max(data.max_interval_size, r.size());
	}

	out.write(uint64_t{size(indices)});

	// sort the chromosomes so files are stable between hash_map implementations
	std::vector<typename decltype(indices)::key_type> keys;
	keys.reserve(std::size(indices));
	std::transform(std::begin(indices), std::end(indices), std::back_inserter(keys), [](auto& kv) { return kv.first; });
	std::sort(std::begin(keys), std::end(keys));

	// Sort all the indices by their respective position values, irrespective of strand
	compare_pos5<T> cmp5(rows.empty() ? nullptr : &rows[0]);
	compare_pos3<T> cmp3(rows.empty() ? nullptr : &rows[0]);
	vector<index_t> by_pos5;

	for (auto chrom_strand : keys) {
		auto [chrom, strand] = chrom_strand;
		out.write(chrom);
		out.write(strand);

		auto& data = indices[chrom_strand];
		// Write out the maximum interval sizes
		out.write(data.max_interval_size);

		stable_sort(data.by_pos3, cmp3); // avoid spurious binary diffs
		by_pos5 = data.by_pos3;  // Copy indices sorted by pos3
		stable_sort(by_pos5, cmp5);  // Stable sort by pos5, so that secondary sort is pos3
		out.write_array(by_pos5);
		out.write_array(data.by_pos3);
	}
}

template <typename T>
void interval_idx<T>::load(mmap_file& in)
{
	auto num_indices = in.read<uint64_t>();
	for (uint64_t i = 0; i < num_indices; ++i) {
		auto chrom  = in.read<chrom_t>();
		auto strand = in.read<strand_t>();

		auto& data             = indices_by_stranded_key[{chrom, strand}];
		data.max_interval_size = in.read<pos_t>();
		data.by_pos5.read(in);
		data.by_pos3.read(in);
	}
}

template <typename T>
const typename interval_idx<T>::sorted_indices& interval_idx<T>::get_indices(chrom_t chrom, strand_t strand) const
{
	static constexpr sorted_indices null_indices{};
	return find_or(indices_by_stranded_key, typename decltype(indices_by_stranded_key)::key_type{chrom, strand},
				   null_indices);
}

////////////////////////////////////////////////////////
////////////////////////////////////////////////////////
////////////////////////////////////////////////////////

template <typename T>
typename interval_table<T>::cursor_range_within interval_table<T>::find_within(const interval_t& i) const
{
	if (!this->_stranded && i.is_neg_strand())
		return find_within(i.as_pos_strand());

	// Get all intervals with pos3 appearing in closed interval [i.pos5, i.pos3]
	// Rely on the interval_cursors to filter out elements with pos5 < i.pos5 precedes
	return cursor_range_within(_find_3p_within(i), i.pos5);
}

template <typename T>
typename interval_table<T>::cursor_range_overlapping interval_table<T>::find_overlapping(const interval_t& i) const
{
	if (!this->_stranded && i.is_neg_strand())
		return find_overlapping(i.as_pos_strand());

	// Get all intervals with pos5 appearing in closed interval [i.pos5-max_size, i.pos3]
	// Rely on the interval_cursors to filter out the extra intervals with pos3 < i.pos5
	// but which were included despite not overlapping with the query interval.
	pos_t max_size = this->_idx.max_interval_size(i.chrom, i.strand);
	interval_t j = i.expand(max_size, 0);
	return cursor_range_overlapping(_find_5p_within(j), i.pos5);
}

template <typename T>
typename interval_table<T>::cursor_range interval_table<T>::_find_5p_within(const interval_t& i) const
{
	return base_class::template find_by_field<compare_pos5_to_value>(i, this->_idx.by_pos5(i.chrom, i.strand));
}

template <typename T>
typename interval_table<T>::cursor_range interval_table<T>::_find_3p_within(const interval_t& i) const
{
	return base_class::template find_by_field<compare_pos3_to_value>(i, this->_idx.by_pos3(i.chrom, i.strand));
}

template <typename T>
typename interval_table<T>::cursor_range interval_table<T>::find_5p_within(const interval_t& i) const
{
	GK_CHECK(this->_stranded || i.is_pos_strand(), value, "Cannot call find_5p_within on negative strand for unstranded table");
	return _find_5p_within(i);
}

template <typename T>
typename interval_table<T>::cursor_range interval_table<T>::find_3p_within(const interval_t& i) const
{
	GK_CHECK(this->_stranded || i.is_pos_strand(), value, "Cannot call find_3p_within on negative strand for unstranded table");
	return _find_3p_within(i);
}

template <typename T>
typename interval_table<T>::cursor_range interval_table<T>::find_5p_aligned(const interval_t& i) const
{
	GK_CHECK(this->_stranded || i.is_pos_strand(), value, "Cannot call find_5p_aligned on negative strand for unstranded table");

	// Case-by-case proof of correctness, using the (pos5, pos3) representation of intervals
	//
	// POSITIVE STRAND
	//
	// Interval A:
	//
	//   ...AAAA...
	//   0123456789
	//      [  ]        <-- (3, 6)
	//      []]]]]]     <-- (3, *)  => use (i.pos5, i.pos5)
	//
	// Interval A.end5() (empty):
	//
	//   ...AAAA...
	//   0123456789
	//     ][           <-- (3, 2)
	//      []]]]]]     <-- (3, *)  => use (i.pos5, i.pos5)
	//
	// Interval A.end3() (empty):
	//
	//   ...AAAA...
	//   0123456789
	//         ][       <-- (7, 6)
	//          []]     <-- (7, *)  => use (i.pos5, i.pos5)
	//
	// NEGATIVE STRAND
	//
	// Interval A:
	//
	//   ...AAAA...
	//   0123456789
	//      [  ]        <-- (6, 3)
	//   [[[[[[]        <-- (6, *)  => use (i.pos5, i.pos5)
	//
	// Interval A.end5() (empty):
	//
	//   ...AAAA...
	//   0123456789
	//         ][       <-- (6, 7)
	//   [[[[[[]        <-- (6, *)  => use (i.pos5, i.pos5)
	//
	// Interval A.end3() (empty):
	//
	//   ...AAAA...
	//   0123456789
	//     ][           <-- (2, 3)
	//   [[]            <-- (2, *)  => use (i.pos5, i.pos5)
	//
	return _find_5p_within(i.end5());
}

template <typename T>
typename interval_table<T>::cursor_range interval_table<T>::find_3p_aligned(const interval_t& i) const
{
	GK_CHECK(this->_stranded || i.is_pos_strand(), value, "Cannot call find_3p_aligned on negative strand for unstranded table");

	// Case-by-case proof of correctness, using the (pos5, pos3) representation of intervals
	//
	// POSITIVE STRAND
	//
	// Interval A:
	//
	//   ...AAAA...
	//   0123456789
	//      [  ]        <-- (3, 6)
	//   [[[[[[]        <-- (*, 6)  => use (i.pos3, i.pos3)
	//
	// Interval A.end5() (empty):
	//
	//   ...AAAA...
	//   0123456789
	//     ][           <-- (3, 2)
	//   [[]            <-- (*, 2)  => use (i.pos3, i.pos3)
	//
	// Interval A.end3() (empty):
	//
	//   ...AAAA...
	//   0123456789
	//         ][       <-- (7, 6)
	//   [[[[[[]        <-- (*, 6)  => use (i.pos3, i.pos3)
	//
	// NEGATIVE STRAND
	//
	// Interval A:
	//
	//   ...AAAA...
	//   0123456789
	//      [  ]        <-- (6, 3)
	//      []]]]]]     <-- (*, 3)  => use (i.pos3, i.pos3)
	//
	// Interval A.end5() (empty):
	//
	//   ...AAAA...
	//   0123456789
	//         ][       <-- (6, 7)
	//          []]     <-- (*, 7)  => use (i.pos3, i.pos3)
	//
	// Interval A.end3() (empty):
	//
	//   ...AAAA...
	//   0123456789
	//     ][           <-- (2, 3)
	//      []]]]]]     <-- (*, 3)  => use (i.pos3, i.pos3)
	//
	return _find_3p_within(i.end3());
}

template <typename T>
typename interval_table<T>::cursor_range interval_table<T>::find_exact(const interval_t& i) const
{
	if (!this->_stranded && i.is_neg_strand())
		return find_exact(i.as_pos_strand());

	GK_CHECK(this->valid(), file, "IntervalTable has been invalidated by close or with statement on source object.");
	if (!this->_elems.empty()) {
		// TODO: data_dir injected as a context
		const auto& refgs = get_refg_registry();
		GK_CHECK(i.refg == this->_elems[0].refg, value, "Reference genome '{}' does not match query coordinate '{}'.",
				 refgs.refg_as_sv(this->_elems[0].refg), refgs.refg_as_sv(i.refg));
	}

	// The pos5 index is deliberately sorted secondarily by pos3.
	// This facilitates a fast two-stage binary search:
	//   1) first find the [first, last) of all intervals with pos5 matching i.pos5
	//   2) then narrow down the [first, last) to be the sub-range with pos3 matching i.pos3
	const array_view<index_t>& idx = this->_idx.by_pos5(i.chrom, i.strand);
	index_iter first = idx.begin(), last = idx.end();
	first = lower_bound(first, last, i.pos5, typename compare_pos5_to_value<T>::lower(this->_elems.data()));
	last  = upper_bound(first, last, i.pos5, typename compare_pos5_to_value<T>::upper(this->_elems.data()));
	first = lower_bound(first, last, i.pos3, typename compare_pos3_to_value<T>::lower(this->_elems.data()));
	last  = upper_bound(first, last, i.pos3, typename compare_pos3_to_value<T>::upper(this->_elems.data()));
	return cursor_range(first, last, this->_elems.data());
}

template <typename T> template <typename P>
INLINE void interval_table<T>::filtered_cursor<P>::advance_until_valid()
{
	P get_pos;  // Returns either pos3 or pos5
	GK_DBASSERT(cur < last);
	const T* row = elems + (*cur);
	if (row->is_pos_strand()) {
		while (get_pos(*row) < pos5 && ++cur < last)
			row = elems + (*cur);
	} else {
		while (get_pos(*row) > pos5 && ++cur < last)
			row = elems + (*cur);
	}
}

END_NAMESPACE_GK

#endif // __GENOME_KIT_TABLE_H__
