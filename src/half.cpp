/*
Copyright (C) 2016-2023 Deep Genomics Inc. All Rights Reserved.
*/
#include "half.h"
#include "gk_assert.h"

BEGIN_NAMESPACE_GK

float clamp_for_half(float f)
{
	if (f >  65504.0f) return  65504.0f;
	if (f < -65504.0f) return -65504.0f;
	return f;
}

float overflow_for_half(float f)
{
	union { float f; uint32_t fbits; } r;
	if (f >  65504.0f) { r.fbits = 0x7f800000; return r.f; } // +inf
	if (f < -65504.0f) { r.fbits = 0xff800000; return r.f; } // -inf
	return f;
}

// Adapted from Numpy halffloat.c
INLINE uint16_t float_as_halfbits(float x)
{
	union { float f; int32_t fbits; } conv;
	conv.f = x;
	uint32_t f = conv.fbits;
	uint32_t f_exp, f_sig;
	uint16_t h_sgn, h_exp, h_sig;
	h_sgn = (uint16_t) ((f & 0x80000000u) >> 16);
	f_exp = (f & 0x7f800000u);
	if (f_exp >= 0x47800000u) {
		if (f_exp == 0x7f800000u) { // Inf or NaN
			f_sig = (f & 0x007fffffu);
			if (f_sig != 0) {       // NaN - propagate the flag in the significand...
				auto ret = (uint16_t)(0x7c00u + (f_sig >> 13));
				if (ret == 0x7c00u)
					ret++;
				return h_sgn + ret;
			} else
				return (uint16_t)(h_sgn + 0x7c00u); // Signed inf
		} else {
			GK_THROW(value, "Overflow to inf detected when converting float32 value {} to float16. Try clamping or scaling values to a range no larger than [-65504, +65504]", x);
			//return (uint16_t)(h_sgn + 0x7c00u); // Overflow to signed inf
		}
	}
	if (f_exp <= 0x38000000u) { // Exponent underflow converts to a subnormal half or signed zero
		if (f_exp < 0x33000000u)
			return h_sgn;
		f_exp >>= 23;
		f_sig = (0x00800000u + (f&0x007fffffu));
		f_sig >>= (113 - f_exp);
		f_sig += 0x00001000u;
		h_sig = (uint16_t)(f_sig >> 13);
		return (uint16_t)(h_sgn + h_sig);
	}
	// Regular case (no overflow or underflow)
	h_exp = (uint16_t)((f_exp - 0x38000000u) >> 13);
	f_sig = (f&0x007fffffu);
	f_sig += 0x00001000u;
	h_sig = (uint16_t)(f_sig >> 13);
	return h_sgn + h_exp + h_sig;
}

half_t::half_t(float    x): h(float_as_halfbits(x)) { }
half_t::half_t(double   x): h(float_as_halfbits((float)x)) { }
half_t::half_t(int      x): h(float_as_halfbits((float)x)) { }
half_t::half_t(unsigned x): h(float_as_halfbits((float)x)) { }

// Adapted from Numpy halffloat.c
uint32_t _as_float_special(half_t h)
{
	uint16_t h_exp, h_sig;
	uint32_t f_sgn, f_exp, f_sig;
	h_exp = (h & 0x7c00u);
	f_sgn = ((uint32_t)h.h & 0x8000u) << 16;
	if (h_exp == 0x7c00u)
		return f_sgn + 0x7f800000u + (((uint32_t)(h & 0x03ffu)) << 13); // inf or NaN

	// 0 or subnormal
	h_sig = (h & 0x03ffu);
	if (h_sig == 0)
		return f_sgn; // Signed zero
	h_sig <<= 1;
	while ((h_sig & 0x0400u) == 0) {
		h_sig <<= 1;
		h_exp++;
	}
	f_exp = ((uint32_t)(127 - 15 - h_exp)) << 23;
	f_sig = ((uint32_t)(h_sig & 0x03ffu)) << 13;
	return f_sgn + f_exp + f_sig;
}

END_NAMESPACE_GK
