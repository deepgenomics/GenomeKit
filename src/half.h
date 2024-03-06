/*
Copyright (C) 2016-2023 Deep Genomics Inc. All Rights Reserved.
*/
#pragma once
#ifndef __GENOME_KIT_HALF_H__
#define __GENOME_KIT_HALF_H__

#include "defines.h"
#include <cmath>
#include <fmt/format.h>
#include <limits>

BEGIN_NAMESPACE_GK

// Bitwise operations on half are adapted from Numpy halffloat.c
union half_t {
	uint16_t h{};
	half_t() = default;
	INLINE half_t(uint16_t h): h(h) { }
	       half_t(float    x);  // error if overflows -- see clamp_for_half() and overflow_for_half()
	       half_t(double   x);  // error if overflows
	       half_t(int      x);  // error if overflows
	       half_t(unsigned x);  // error if overflows
};

float clamp_for_half(float f);    // clamps numeric values to range [-65504, +65504], so they can be safely converted by to half_t
float overflow_for_half(float f); // overflows numeric values that are beyond [-65504, +65504] to be -inf or +inf, so they can be safely converted to half_t

// Adapted from Numpy halffloat.c
INLINE uint16_t operator&(half_t   x, unsigned y) { return x.h & y; }
INLINE uint16_t operator&(unsigned x, half_t   y) { return x & y.h; }
INLINE bool is_nan(half_t h)    { return ((h & 0x7c00u) == 0x7c00u) && ((h & 0x03ffu) != 0); }
INLINE bool is_zero(half_t h)   { return  (h & 0x7fff) == 0; }
INLINE bool is_inf(half_t h)    { return ((h & 0x7fffu) == 0x7c00u); }
INLINE bool is_finite(half_t h) { return ((h & 0x7c00u) != 0x7c00u); }
INLINE bool sign_bit(half_t h)  { return  (h & 0x8000u) != 0; }
INLINE bool lt_nonan(half_t x, half_t y);
INLINE bool le_nonan(half_t x, half_t y);
INLINE bool eq_nonan(half_t x, half_t y) { return x.h == y.h || ((x.h | y.h) & 0x7fff) == 0; }
INLINE bool operator==(half_t x, half_t y) { return (!is_nan(x) && !is_nan(y)) && eq_nonan(x, y); }
INLINE bool operator!=(half_t x, half_t y) { return !(x == y); }
INLINE bool operator< (half_t x, half_t y) { return (!is_nan(x) && !is_nan(y)) && lt_nonan(x, y); }
INLINE bool operator<=(half_t x, half_t y) { return (!is_nan(x) && !is_nan(y)) && le_nonan(x, y); }
INLINE bool operator> (half_t x, half_t y) { return y <  x; }
INLINE bool operator>=(half_t x, half_t y) { return y <= x; }
INLINE half_t min(half_t x, half_t y) { return x < y ? x : y; }
INLINE half_t max(half_t x, half_t y) { return x > y ? x : y; }
INLINE half_t nanh() { return (uint16_t)0x7fffu; }

INLINE bool is_nan(float f)     { return std::isnan(f); }
INLINE bool is_inf(float  f)    { return std::isinf(f); }
INLINE bool is_finite(float  f) { return std::isfinite(f); }
constexpr float nanf()          { return std::numeric_limits<float>::quiet_NaN(); }

// Adapted from Numpy halffloat.c
uint32_t _as_float_special(half_t);
INLINE float as_float(float  f) { return f; }
INLINE float as_float(half_t h) // inline the normal case for speed since it's just a few bitwise ops
{
	union { float ret; uint32_t retbits; } conv;
	uint16_t h_exp;
	uint32_t f_sgn;
	h_exp = (h & 0x7c00u);
	f_sgn = ((uint32_t)h.h & 0x8000u) << 16;
	if (h_exp != 0x0000u && h_exp != 0x7c00u)
		conv.retbits = f_sgn + (((uint32_t)(h & 0x7fffu) + 0x1c000u) << 13); // normal
	else
		conv.retbits = _as_float_special(h); // don't inline the inf, NaN, or subnormal cases
	return conv.ret;
}

// Adapted from Numpy halffloat.c
INLINE bool lt_nonan(half_t x, half_t y)
{
    if (x & 0x8000u) return (y & 0x8000u) ?         (x & 0x7fffu) > (y & 0x7fffu) : (x.h != 0x8000u) || (y.h != 0x0000u);
    else             return (y & 0x8000u) ? false : (x & 0x7fffu) < (y & 0x7fffu);
}

// Adapted from Numpy halffloat.c
INLINE bool le_nonan(half_t x, half_t y)
{
    if (x & 0x8000u) return (y & 0x8000u) ? (x & 0x7fffu) >= (y & 0x7fffu) : true;
    else             return (y & 0x8000u) ? (x.h == 0x0000u) && (y.h == 0x8000u) : (x & 0x7fffu) <= (y & 0x7fffu);
}

END_NAMESPACE_GK

template <>
struct fmt::formatter<gk::half_t> : fmt::formatter<float> {
	template <typename FormatCtx>
	auto format(gk::half_t x, FormatCtx& ctx) const
	{
		return fmt::formatter<gk::half_t>::format(as_float(x), ctx);
	}
};

#endif // __GENOME_KIT_HALF_H__
