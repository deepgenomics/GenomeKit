/*
Copyright (C) 2016-2023 Deep Genomics Inc. All Rights Reserved.
*/
#pragma once
#ifndef __GENOME_KIT_INTERVAL_H__
#define __GENOME_KIT_INTERVAL_H__

#include "chrom.h"
#include "defines.h"
#include "gk_assert.h"
#include "refg.h"
#include "util.h"

#include <climits>
#include <concepts>
#include <format>
#include <functional>
#include <string>
#include <string_view>
#include <unordered_map>
#include <vector>

BEGIN_NAMESPACE_GK
using std::max;
using std::min;
using std::string;

using pos_t = int32_t;          // Position within a coordinate system
enum class strand_t : uint8_t { // Strand index [+,-]
	neg_strand,
	pos_strand,
	num_strand // [+,-] is two unique identifiers mapped to [0,1]
};

using offset_t = uint32_t;   // Byte offset type. A pointer stored as 4-byte relative offset rather than as 8-byte
							 // absolute pointer. For easy serialization and memory savings.
using offset64_t = uint64_t; // Larger offset

static constexpr auto neg_strand = strand_t::neg_strand;
static constexpr auto pos_strand = strand_t::pos_strand;
static constexpr auto num_strand = as_ordinal(strand_t::num_strand);

/////////////////////////////////////////////////////////////////
// conversion routines
/////////////////////////////////////////////////////////////////

// Convert between pos_t and string
pos_t as_pos(const std::string_view s);

// Convert between '+'/'-' and strand_t
INLINE strand_t as_strand(char c)            { GK_CHECK(c == '+' || c == '-', value, "Expected strand to be '+' or '-' but found '{}'.", c); return c=='+' ? pos_strand : neg_strand; }
INLINE strand_t as_strand(std::string_view s)
{
	GK_CHECK(s.size() == 1, value, "Expected strand string \"{}\" to be \"+\" or \"-\".", s);
	return as_strand(s[0]);
}
INLINE bool is_valid_strand(strand_t strand) { return as_ordinal(strand) < num_strand; }
INLINE char strand_as_char(strand_t strand)  { GK_DBASSERT(is_valid_strand(strand)); return strand == pos_strand ? '+' : '-'; }

INLINE pos_t shift(pos_t pos, strand_t strand, pos_t amount) { return strand == pos_strand ? pos+amount : pos-amount; }
INLINE bool is_upstream(pos_t a, pos_t b, strand_t strand)   { return strand == pos_strand ? a < b : b < a; }
INLINE bool is_dnstream(pos_t a, pos_t b, strand_t strand)   { return strand == pos_strand ? a > b : b > a; }

template <class T>
struct chrom_key_t {
	chrom_t  chrom{};
	T other{};

	auto operator<=>(const chrom_key_t&) const = default;
};
template <class T>
struct chrom_key_hash_t {
	std::size_t operator()(chrom_key_t<T> k) const
	{
		if constexpr (sizeof(k.chrom) + sizeof(k.other) <= sizeof(std::size_t)) {
			return (scast<std::size_t>(k.chrom) << sizeof(k.other) * CHAR_BIT) | scast<std::size_t>(k.other);
		} else {
			static_assert(sizeof(k.chrom) <= sizeof(std::size_t) && sizeof(k.other) <= sizeof(std::size_t));
			return scast<std::size_t>(k.chrom) ^ scast<std::size_t>(k.other);
		}
	}
};

template <class T, class V>
using chrom_map_t = std::unordered_map<chrom_key_t<T>, V, chrom_key_hash_t<T>>;

/////////////////////////////////////////////////////////////////
// coordinate structs
/////////////////////////////////////////////////////////////////

#pragma pack(push, 1)  // Ensures subclasses don't waste 2 bytes on alignment padding
struct interval_t {
	pos_t pos5;        // 4 bytes - 0-based
	pos_t pos3;        // 4 bytes - 0-based inclusive
	refg_t refg;	   // 8 bytes - fnv1a_hash64 of name
	chrom_t chrom;     // 4 bytes - fnv1a_hash32 of name
	strand_t strand;   // 1 byte

	// TODO: with current constructor overloads, C++ user code might accidentally do interval_t(chrom, pos5, pos3, strand) and forget refg, but it would still default to hg19; that's not good
	INLINE static interval_t from_dna0(chrom_t chrom, pos_t start, pos_t end, strand_t strand, refg_t refg) { return interval_t(chrom, strand == pos_strand ? start : end-1, strand == pos_strand ? end-1 : start, strand, refg); }
	       interval_t() = default;
	INLINE interval_t(chrom_t chrom, pos_t pos5, pos_t pos3, strand_t strand, refg_t refg): pos5(pos5), pos3(pos3), refg(refg), chrom(chrom), strand(strand) { GK_DBASSERT(is_valid_strand(strand)); }
	INLINE pos_t   start()  const { return is_pos_strand() ? pos5   : pos3;   } // 0-based
	INLINE pos_t   end()    const { return is_pos_strand() ? pos3+1 : pos5+1; } // 0-based exclusive
	INLINE pos_t   size()   const { return is_pos_strand() ? pos3-pos5+1 : pos5-pos3+1; }
	INLINE bool    empty()  const { return size() == 0; }
	INLINE interval_t end5() const { return with_pos(pos5, pos5); }
	INLINE interval_t end3() const { return with_pos(pos3, pos3); }
	INLINE interval_t shift(pos_t amount) const { return with_pos(gk::shift(pos5, strand, amount), gk::shift(pos3, strand, amount)); }
	INLINE interval_t expand(pos_t upstream, pos_t dnstream) const { return with_pos(gk::shift(pos5, strand, -upstream), gk::shift(pos3, strand,  dnstream)); }
	INLINE interval_t shrink(pos_t upstream, pos_t dnstream) const { return with_pos(gk::shift(pos5, strand,  upstream), gk::shift(pos3, strand, -dnstream)); }
	INLINE interval_t slice(pos_t first, pos_t last)         const { return with_pos(gk::shift(pos5, strand,  first),    gk::shift(pos5, strand,  last-1)); }
	INLINE bool upstream_of(const interval_t& i) const { return (is_pos_strand() ? pos3 < i.pos5 : pos3 > i.pos5) && same_strand(i); }
	INLINE bool dnstream_of(const interval_t& i) const { return (is_pos_strand() ? pos5 > i.pos3 : pos5 < i.pos3) && same_strand(i); }
	INLINE bool contains(const interval_t& i) const { return !is_upstream(i.pos5, pos5, strand) && !is_dnstream(i.pos3, pos3, strand) && same_strand(i); }
	INLINE bool within(const interval_t& i)   const { return i.contains(*this); }
	INLINE bool overlaps(const interval_t& i) const { return !is_dnstream(pos5, i.pos3, strand) && !is_upstream(pos3, i.pos5, strand) && same_strand(i); }
	INLINE bool       is_pos_strand()         const { return strand == pos_strand; }
	INLINE bool       is_neg_strand()         const { return strand == neg_strand; }
	INLINE interval_t as_pos_strand()         const { return with_strand(pos_strand); }
	INLINE interval_t as_neg_strand()         const { return with_strand(neg_strand); }
	INLINE interval_t as_opp_strand()         const { return with_strand(is_pos_strand() ? neg_strand : pos_strand); }
	INLINE const interval_t& as_interval() const { return *this; }
	INLINE       interval_t& as_interval()       { return *this; }
	string as_str() const; // deliberately not virtual, we do NOT want VMT on interval_t

	INLINE interval_t with_pos(pos_t new_pos5, pos_t new_pos3) const { return interval_t{chrom, new_pos5, new_pos3, strand, refg}; }
	INLINE interval_t with_strand(strand_t new_strand) const { return strand == new_strand ? *this : interval_t{chrom, pos3, pos5, new_strand, refg}; }
	INLINE bool same_chrom(const interval_t& i) const { return chrom == i.chrom && refg == i.refg; }
	INLINE bool same_strand(const interval_t& i) const { return same_chrom(i) && strand == i.strand; }
};
#pragma pack(pop)

INLINE bool operator==(const interval_t& x, const interval_t& y) { return x.pos5 == y.pos5 && x.pos3 == y.pos3 && x.same_strand(y); }
INLINE bool operator!=(const interval_t& x, const interval_t& y) { return !(x == y); }
INLINE bool operator< (const interval_t& x, const interval_t& y) {
	return x.refg != y.refg       ? as_ordinal(x.refg) < as_ordinal(y.refg)
		   : x.chrom != y.chrom   ? as_ordinal(x.chrom) < as_ordinal(y.chrom)
		   : x.pos5 != y.pos5     ? x.pos5 < y.pos5
		   : x.strand != y.strand ? as_ordinal(x.strand) < as_ordinal(y.strand)
								  : x.pos3 < y.pos3;
}
INLINE bool operator> (const interval_t& x, const interval_t& y) { return y < x; }
INLINE bool operator<=(const interval_t& x, const interval_t& y) { return x<y || x==y; }
INLINE bool operator>=(const interval_t& x, const interval_t& y) { return x>y || x==y; }

INLINE size_t hash(const interval_t& c)
{
	size_t x = 0;
	x ^= ((size_t)c.pos5+23)*102564;  // TODO use hashing scheme with better properties
	x ^= ((size_t)c.pos3+562123)*73;
	x ^= (((size_t)c.strand) << 6);
	x ^= chrom_key_hash_t<refg_t>()({c.chrom, c.refg});
	return x;
}

// UNPACK_INTERVAL macros create local variable copies of interval fields.
// Makes for convenient access to the fields, plus frees the optimizer
// from worrying about aliasing. If any of these fields are not used in
// the code that follows, they'll be optimized away so there's no penalty.

#define UNPACK_INTERVAL(c) \
	[[maybe_unused]] pos_t pos5 = c.pos5; \
	[[maybe_unused]] pos_t pos3 = c.pos3; \
	[[maybe_unused]] chrom_t chrom = c.chrom; \
	[[maybe_unused]] refg_t refg = c.refg; \
	[[maybe_unused]] strand_t strand = c.strand;


// Handy functors for extracting a field known at compile time.
template <typename T> struct get_pos3 {
	INLINE pos_t operator()(const T& elem) { return elem.pos3; }
};

template <typename T> struct get_pos5 {
	INLINE pos_t operator()(const T& elem) { return elem.pos5; }
};

class interval_filter {
public:
	using interval_fn_t = std::function<void(interval_t)>;

	explicit interval_filter(interval_fn_t validate_interval = {});
	void allow(interval_t i);
	void exclude(interval_t i);
	bool filter(interval_t i) const;
	void validate() const; // revalidate the intervals

private:
	std::function<void(interval_t)> _validate_interval;
	std::vector<interval_t>         _allow;
	std::vector<interval_t>         _exclude;
};

END_NAMESPACE_GK

template <typename T>
	requires std::derived_from<T, gk::interval_t>
struct std::formatter<T> : std::formatter<std::string> {
	// T instead of interval_t since as_str is not virtual
	template <typename FormatCtx>
	auto format(const T& x, FormatCtx& ctx) const
	{
		return std::formatter<std::string>::format(x.as_str(), ctx);
	}
};

template <>
struct std::formatter<gk::strand_t> : std::formatter<std::uint8_t> {
	template <typename FormatCtx>
	auto format(gk::strand_t x, FormatCtx& ctx) const
	{
		return std::formatter<std::uint8_t>::format(static_cast<std::uint8_t>(x), ctx);
	}
};

#endif // __GENOME_KIT_INTERVAL_H__
