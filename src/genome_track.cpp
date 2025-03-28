/*
Copyright (C) 2016-2023 Deep Genomics Inc. All Rights Reserved.
*/
#include "genome_track.h"
#include <algorithm>
#include <numeric>
#include <optional>
#include <utility>
#include <variant>

BEGIN_NAMESPACE_GK

const char* genome_track::etype_as_cstr[genome_track::num_etype] = {
	"m0",
	"u1",
	"u2",
	"u3",
	"u4",
	"u5",
	"u6",
	"u8",
	"i8",
	"f2",
	"f3",
	"f4",
	"f5",
	"f6",
	"f8",
	"f16",
	"f32"
};

// Default decoded types for various encodings.
genome_track::dtype_t genome_track::etype_default_dtype[genome_track::num_etype] = {
	genome_track::bool_,    // m0
	genome_track::bool_,    // u1
	genome_track::uint8,    // u2
	genome_track::uint8,    // u3
	genome_track::uint8,    // u4
	genome_track::uint8,    // u5
	genome_track::uint8,    // u6
	genome_track::uint8,    // u8
	genome_track::int8,     // i8
	genome_track::float16,  // f2
	genome_track::float16,  // f3
	genome_track::float16,  // f4
	genome_track::float16,  // f5
	genome_track::float16,  // f6
	genome_track::float16,  // f8
	genome_track::float16,  // f16
	genome_track::float32,  // f32
};

const char* genome_track::dtype_as_cstr[as_ordinal(genome_track::num_dtype)] = {
	"bool",
	"uint8",
	"int8",
	"float16",
	"float32",
};

int genome_track::dtype_size[as_ordinal(genome_track::num_dtype)] = {
	sizeof(bool),    // bool_
	sizeof(uint8_t), // uint8
	sizeof(int8_t),  // int8
	sizeof(half_t),  // float16
	sizeof(float),   // float32
};

genome_track::dtype_t genome_track::as_dtype(const char* s)
{
	for (int i = 0; i < as_ordinal(genome_track::num_dtype); ++i)
		if (!strcmp(s, dtype_as_cstr[i]))
			return (genome_track::dtype_t)i;
	GK_THROW(value, "Unrecognized dtype '{}'", s);
}

genome_track::etype_t genome_track::as_etype(const char* s)
{
	for (int i = 0; i < as_ordinal(genome_track::num_etype); ++i)
		if (!strcmp(s, etype_as_cstr[i]))
			return (genome_track::etype_t)i;
	GK_THROW(value, "Unrecognized etype '{}'", s);
}

////////////////////////////////////////////////////////////

void genome_track::set_source(string sourcefile)
{
	// Try not to put anything here that would throw in a newly constructed object.
	// We want genome_t constructor to not throw, even if it's initialized outside
	GK_CHECK2(!_fmap.is_open(), runtime, "Cannot set source when file already open.");
	_sourcefile = std::move(sourcefile);
}

void genome_track::close()
{
	_fmap.close();
	_sourcefile.clear();
}

void genome_track::open_on_demand() const
{
	// TODO: acquire lock here and check _file for NULL again

	// Once we're inside the lock it's safe to use ncthis
	auto* ncthis = const_cast<genome_track*>(this);
	ncthis->open();

	// TODO: release lock here (implicitly, when falls out of scope)
}

void genome_track::operator()(const interval_t& c, bool*    dst, int stride) const { (*this)(c, dst, bool_,   stride); }
void genome_track::operator()(const interval_t& c, uint8_t* dst, int stride) const { (*this)(c, dst, uint8,   stride); }
void genome_track::operator()(const interval_t& c, int8_t*  dst, int stride) const { (*this)(c, dst, int8,    stride); }
void genome_track::operator()(const interval_t& c, half_t*  dst, int stride) const { (*this)(c, dst, float16, stride); }
void genome_track::operator()(const interval_t& c, float*   dst, int stride) const { (*this)(c, dst, float32, stride); }
void genome_track::operator()(const interval_t& c, void*    dst, dtype_t dtype, int stride) const
{
	ensure_open();
	GK_CHECK2(refg() == c.refg, value, "Reference genome mismatch");

	if (stride == 0)
		stride = dim();
	GK_CHECK(stride > 0, value, "Negative strides not supported: stride={}", stride);
	GK_CHECK2(stride >= dim(), value, "Stride is too small: stride={}, dim={}", stride, dim());
	const int layout = as_ordinal(stride == dim() ? encoding::layout_t::contiguous : encoding::layout_t::noncontiguous);

	// Get callbacks that are specialized to decode and default fill for this dtype and strand direction.
	encoding::decode_fn decode = _encoding.decoders[as_ordinal(dtype)][layout][as_ordinal(c.strand)];
	encoding::dfill_fn  dfill  = _encoding.dfillers[as_ordinal(dtype)][layout][as_ordinal(c.strand)];
	GK_CHECK(decode, type, "Cannot decode as {} from encoded type {}", dtype_as_cstr[as_ordinal(dtype)], etype_as_cstr[_encoding.etype]);
	GK_DBASSERT(dfill);

	static constexpr track_index_t null_index{};
	static constexpr track_info_t  null_info{.index = &null_index};
	const auto ct =
		find_or(_tracks, decltype(_tracks)::key_type{c.chrom, stranded() ? c.strand : pos_strand}, null_info);

	const void* RESTRICT          data = ct.data;
	const track_index_t* RESTRICT idx  = ct.index;

	// The index structure is stored in three parallel arrays, sorted by pos3.
	//
	//   int      jumps[num_jumps];          // sorted, increasing order
	//   pos_t    ends[num_blocks];          // sorted, increasing order
	//   pos_t    starts[num_blocks];        // sorted, increasing order
	//   uintNN_t data_offsets[num_blocks];  // relative to data segment for this (chrom, strand)
	//                                       // (NN = 32 or 64 depending on is_32bit_index)
	//
	// The initial binary search is done on 'ends' only, so having it as a separate array
	// makes better use of cache.
	const auto num_blocks = idx->num_blocks;
	const auto jumps      = rcast<const int32_t*>(idx + 1);
	const auto ends       = rcast<const pos_t*>(jumps + idx->num_jumps);
	const auto starts     = rcast<const pos_t*>(ends + num_blocks);

	std::optional<std::variant<const uint32_t*, const uint64_t*>> data_offsets;
	if (_encoding.etype != m0) {
		const auto next_field = starts + num_blocks;
		if (idx->is_32bit)
			data_offsets = rcast<const uint32_t*>(next_field);
		else
			data_offsets = rcast<const uint64_t*>(next_field);
	}

	// To understand how the decoding loop works, consider the following 10bp chromosome
	// with two intervals filled with data labeled UVXYZ and the rest default (.).
	//
	// The 'virtual' chromosome looks like this
	//
	//     0123456789ABCDEF   <-- genomic coordinates (hexadecimal to show as 1 digit)
	//     ....UV..WXYZ....   <-- indexed blocks (4,6) and (9,12)
	//
	// But the data is actually stored in two separate 'src' blocks, like this:
	//
	//     01 0123            <-- index of each value in respective src block
	//     UV WXYZ            <-- encoded values in respective src block
	//
	//     48                 <-- starts array {4, 8}  for src blocks
	//     6C                 <-- ends   array {6, 12} for src blocks
	//
	// On the positive strand, decoding (2,14) should output a 'dst' block:
	//
	//     0123456789AB       <-- index of each value in dst[]
	//     ..UV..WXYZ..       <-- decoded values
	//
	// On the negative strand, the decoding process is identical but the data is
	// stored in reverse order in the output, starting from the last position.
	//
	//     0123456789AB
	//     ..ZYXW..VU..
	//
	// To explain how decoding works, we'll unroll the state of the loop variables
	// as they go build these outputs, step by step, in the comments below.
	//
	// ----------------------------------------------
	// EXAMPLE 1: decoding (2,14) on positive strand
	// ----------------------------------------------
	//
	// STATE 1: entering main loop, after initial value of `i` is determined
	//          by calling upper_bound on `ends` array {6,12}
	//
	//     genomic track          data        index  output array
	//       b           a                    i      d
	//     0123456789ABCDEF       01 0123     48     0123456789AB
	//     ....UV..WXYZ....       UV WXYZ     6C     ????????????
	//
	// STATE 2: determine `a` for next default range [b,a) and fill it into index `d`
	//       b a                              i      d
	//     0123456789ABCDEF       01 0123     48     0123456789AB
	//     ....UV..WXYZ....       UV WXYZ     6C     ..??????????
	//
	// STATE 3: determine `b` for next data range [a,b) and decode it from index `s` into index `d`
	//         a b                s           i        d
	//     0123456789ABCDEF       01 0123     48     0123456789AB
	//     ....UV..WXYZ....       UV WXYZ     6C     ..UV????????
	//
	// STATE 4: determine `a` for next default range [b,a) and fill it into index `d`
	//           b a                           i         d
	//     0123456789ABCDEF       01 0123     48     0123456789AB
	//     ....UV..WXYZ....       UV WXYZ     6C     ..UV..??????
	//
	// STATE 5: determine `b` for next data range [a,b) and decode it from index `s` into index `d`
	//             a   b             s         i           d
	//     0123456789ABCDEF       01 0123     48     0123456789AB
	//     ....UV..WXYZ....       UV WXYZ     6C     ..UV..WXYZ??
	//
	// STATE 6: determine `a` for final default range [b,a) and fill it into index `d`
	//                 b a                      i              d
	//     0123456789ABCDEF       01 0123     48     0123456789AB
	//     ....UV..WXYZ....       UV WXYZ     6C     ..UV..WXYZ..
	//
	// DONE!
	//
	// ----------------------------------------------
	// EXAMPLE 2: decoding (5,10) on positive strand
	// ----------------------------------------------
	// Here the output should be a dst array containing
	//
	//     01234
	//     V..WX
	//
	// STATE 1: entering main loop, after initial value of `i` is determined by
	//          calling upper_bound on `ends` array {6,12}
	//
	//     genomic track          data        index  output array
	//          b    a                        i      d
	//     0123456789ABCDEF       01 0123     48     01234
	//     ....UV..WXYZ....       UV WXYZ     6C     ?????
	//
	// STATE 2: determine `a` for next default range [b,a) and do nothing since it's empty
	//          a
	//          b                             i      d
	//     0123456789ABCDEF       01 0123     48     01234
	//     ....UV..WXYZ....       UV WXYZ     6C     ?????
	//
	// STATE 3: determine `b` for next data range [a,b) and decode it from index `s` into index `d`
	//          ab                 s          i      d
	//     0123456789ABCDEF       01 0123     48     01234
	//     ....UV..WXYZ....       UV WXYZ     6C     V????
	//
	// STATE 4: determine `a` for next default range [b,a) and fill it into index `d`
	//           b a                           i      d
	//     0123456789ABCDEF       01 0123     48     01234
	//     ....UV..WXYZ....       UV WXYZ     6C     V..??
	//
	// STATE 5: determine `b` for next data range [a,b) and decode it from index `s` into index `d`
	//             a b               s         i        d
	//     0123456789ABCDEF       01 0123     48     01234
	//     ....UV..WXYZ....       UV WXYZ     6C     V..WX
	//
	// STATE 6: determine `a` for final default range [b,a) and do nothing since it's empty
	//               a
	//               b                          i         d
	//     0123456789ABCDEF       01 0123     48     01234
	//     ....UV..WXYZ....       UV WXYZ     6C     V..WX
	//
	// DONE!
	//
	// ----------------------------------------------
	// EXAMPLE 3: decoding (2,14) on negative strand
	// ----------------------------------------------
	//
	// All the states are the same, except dst index `d` goes in the opposite direction:
	//
	// STATE 1: entering main loop, after initial `i` is determined by upper_bound on ends[]={6,12}
	//       b           a                    i                 d
	//     0123456789ABCDEF       01 0123     48     0123456789AB
	//     ....UV..WXYZ....       UV WXYZ     6C     ????????????
	//
	// {steps 2-5 same as before, but with `d` decremented}
	//
	// STATE 6: determine `a` for final default range [b,a) and fill it into index `d`
	//                 b a                      i     d
	//     0123456789ABCDEF       01 0123     48     0123456789AB
	//     ....UV..WXYZ....       UV WXYZ     6C     ..ZYXW..UV..
	//
	// DONE!
	//
	// TODO: optimize this code for the case where the query enterval falls entirely inside
	// one of the data blocks, since that's the most common case, one would expect; after all,
	// tracks will tend to be relatively dense.
	//
	const void* dict = (dtype == float16) ? rcast<const void*>(_encoding.dict.h)
	                                      : rcast<const void*>(_encoding.dict.f);
	const any_t& fill = _encoding.default_value;
	pos_t cs = c.start();
	pos_t ce = c.end();
	int dim = _dim;

	if (cs == ce)
		return;

	GK_CHECK(cs >= 0, value, "Interval start {} is invalid", cs);

	if (_res > 1) {
		// If track resolution is >1bp, then move to coarse coordinate system on which data is stored.
		// This will result in only a subset of 'dst' being written to, but we'll repeat those decoded
		// values into the rest of the 'dst' array after the main decoding loop is complete.
		cs = udivdn(cs, _res);  // Round down
		ce = udivup(ce, _res);  // Round up
	}

	// Determine the index i of the leftmost data block to be decoded.
	// Use the jump array to quickly narrow down the scope of the binary search to a jump_size region.
	int i;
	{
		pos_t lo = udivdn(cs, (pos_t)jump_size);
		if (lo >= idx->num_jumps) {
			i = num_blocks;  // Fall through to a complete default_value fill
		} else {
			// Otherwise, do a binary search.
			const pos_t* first = ends + jumps[lo];
			const pos_t* last;
			pos_t hi = lo + 1;
			if (hi >= idx->num_jumps)
				last = ends + num_blocks;
			else
				last = ends + jumps[hi];
			i = (int)(std::upper_bound(first, last, cs) - ends);
		}
	}

	// Initialize b so that final decode() dfills [b,a) if loop doesn't iterate (num_blocks == 0)
	pos_t b = cs;

	// Start outputting to dst[0] for positive strand, dst[n-1] for negative strand
	int d = c.is_pos_strand() ? 0 : ce-cs-1;

	// ----------- MAIN DECODING LOOP ---------------
	for (; i < num_blocks; ++i) {
		// Fill [b,a) with default value.
		int s;
		pos_t a = starts[i];
		if (a <= b) {     // Data block overlaps start of query interval? (only relevant for first iteration of loop)
			s = b-a;      // If so, remember how much of the src block's encoded values we should skip.
			a = b;        // The default fill interval is empty, so do nothing.
		} else {
			s = 0;        // No overlap, so the decoding step will decode from the start of the data block.
			if (a >= ce)  // If this data block is beyond the end of the query interval, truncate at end.
				break;    // (Let the final dfill handle this case)
			d += dfill(dst, fill, a-b, dim, d, stride); // Fill with default
		}

		// Get pointer to the encoded data block (src)
		const void* const src = !data_offsets
									? nullptr
									: std::visit(
										  [&](auto offsets) {
											  using T  = std::remove_pointer_t<decltype(offsets)>;
											  T offset = offsets[i];
											  // If upper bit is set, remaining bits encode data inplace
											  // (fractional_decode will mask the marker+unused bits)
											  return (offset & high_bit<T>())
														 ? scast<const void*>(offsets + i)
														 : scast<const void*>(scast<const char*>(data) + offset);
										  },
										  *data_offsets);

		// Fill [a,b) with src decoded data.
		b = ends[i];
		if (b >= ce) {
			decode(dst, src, dict, ce-a, dim, d, s, stride);
			goto done; // Skip the final dfill
		} else {
			d += decode(dst, src, dict, b-a, dim, d, s, stride);
		}
	}

	// Fill the last [b,a) with default value where a == ce.
	if (ce != b)
		dfill(dst, fill, ce-b, dim, d, stride);

done:

	// REPEATING DECODED DATA UNTIL AT FULL-RESOLUTION
	//
	// If track resolution is >1bp, then all we've done so far is decode the data values
	// at a coarse scale, filling up to half of the 'dst' array with values. We still need
	// to move and repeat those valuesto fill the dst array completely.
	//
	// For example, suppose the track below has distinct values X, Y, Z specified on two.
	// separate intervals:
	//
	//   0  1  2  3  4  5  6  7  8  9  10 11 12 13 14 15  <-- genomic coordinate
	//   .  .  .  X  X  X  Y  Y  Y  .  .  .  Z  Z  Z  .   <-- decoded values
	//
	// If we encode this track at 3bp resolution, what is actually stored in the index are
	// the following two (start, end) intervals and corresponding data blocks:
	//
	//     (1, 3) : [X, Y]   # data for genomic interval (3, 9)
	//     (4, 5) : [Z]      # data for genomic interval (12, 15)
	//
	// Note the start/end in the index are divided by 3 compared to their genomic coordinate.
	//
	// ----------------------------------------------
	// EXAMPLE 1: decoding (2,13) on positive strand
	// ----------------------------------------------
	//
	// Given genomic interval (2,13), the coarse interval that must be decoded
	// is (0,5), which is computed as 2/3=0 and divup(13, 3)=5 (round up).
	// After the main decoding loop finishes decoding (0,5) from the coarse track,
	// the state we're starting from is:
	//
	//     (coarse)   (coarse)  (coarse)   output array (size 11)
	//     track      data      index           d
	//     012345     01 0      14         0123456789A
	//     .XY.Z.     XY Z      35         .XY.Z??????
	//
	// where the last value decoded was 'Z'. The final state of output index d=5.
	//
	// The coarse track, coarse data, and coarse index are no longer needed going forward.
	// Only the output array (dst) is needed.
	//
	// STATE 1: immediately after main decoding loop, initialize source index s=ce-cs (5)
	//          and then new output index d=11.
	//
	//          s     d
	//     0123456789A
	//     .XY.Z??????
	//
	// STATE 2: advance `s` by -1 and advance `d` by `phase = (13 % res)` = -1,
	//          copying values into each `d` along the way; here 13 comes from the genomic
	//          coordinate corresponding to `d` when expanding begins.
	//
	//         s     d
	//     0123456789A
	//     .XY.Z?????Z
	//
	// STATE 3: advance `s` by -1 and advance `d` by `-resolution`, copying values along the way.
	//
	//        s   d
	//     0123456789A
	//     .XY.Z??...Z
	//
	// STATE 4: again
	//
	//       s d
	//     0123456789A
	//     .XY.YYY...Z
	//
	// STATE 5: again
	//      d
	//      s
	//     0123456789A
	//     .XXXYYY...Z
	//
	// STATE 6: again, but stopping as soon as d==0; no new data is copied in this case.
	//     d
	//     s
	//     0123456789A
	//     .XXXYYY...Z
	//
	// The final output (dst) array has now been decoded at full-resolution.
	//
	// ----------------------------------------------
	// EXAMPLE 2: decoding (4,12) on negative strand
	// ----------------------------------------------
	//
	// Given genomic interval (4,12), the coarse interval that must be decoded
	// is (1,4), which is computed as 4/3=1 and divup(12, 3)=4 (round up).
	// The main decoding loop handles the strand orientation, and the
	// subsequent moving/repeating of values is not affected by strand.
	// The state we're starting from is:
	//
	//     (coarse)   (coarse)  (coarse)   output array (size 8)
	//     track      data      index     d
	//     012345     01 0      14         01234567
	//     .XY.Z.     XY Z      35         .YX?????
	//
	// where the last value decoded was '.' (default_value).
	// The final state of output index d=-1.
	//
	// The coarse track, coarse data, and coarse index are no longer needed going forward.
	// Only the output array (dst) is needed.
	//
	// STATE 1: immediately after main decoding loop, initialize source index s=ce-cs (3)
	//          and then new output index d=8.
	//
	//        s    d
	//     01234567
	//     .YX?????
	//
	// STATE 2: advance `s` by -1 and advance `d` by `phase = res - (4 % res) - 1` = 2,
	//          so copy two values; here 4 comes from the genomic coordinate corresponding
	//          to `d` when expanding begins, and we use res-(4%res) because it's reverse strand.
	//
	//       s   d
	//     01234567
	//     .YX???XX
	//
	// STATE 3: advance `s` by -1 and advance `d` by `-resolution`, copying values along the way.
	//
	//      s d
	//     01234567
	//     .YXYYYXX
	//
	// STATE 4: again, but then stop because d == 0
	//     d
	//     s
	//     01234567
	//     ...YYYXX
	//
	// The final output (dst) array has now been decoded at full-resolution.
	//
	if (_res > 1) {
		int phase;
		if (c.is_pos_strand()) {
			phase = umod(c.end(), _res);           // c.end % resolution
		} else {
			phase = _res - umod(c.start(), _res);  // c.start % resolution
			if (phase == _res)
				phase = 0;
		}
		_encoding.expanders[as_ordinal(dtype)][layout](dst, c.size(), dim, ce-cs, _res, phase, stride);
	}
}

genome_track::etype_t genome_track::etype() const { ensure_open(); return _fmap.as_ptr<header_t>(0)->etype; }
genome_track::dtype_t genome_track::dtype() const { return etype_default_dtype[etype()]; }
bool genome_track::supports_dtype(dtype_t dtype) const { ensure_open(); return _encoding.supports_dtype(dtype); }
bool genome_track::empty() const noexcept
{
	ensure_open();
	return std::empty(_tracks);
}

std::vector<interval_t> genome_track::intervals() const
{
	ensure_open();

	std::vector<interval_t> ret(std::accumulate(std::begin(_tracks), std::end(_tracks), 0,
												[](auto num, auto kv) { return num + kv.second.index->num_blocks; }));

	size_t i_itv = 0;
	for (auto [k, v] : _tracks) {
		const auto num_blocks = v.index->num_blocks;
		const auto jumps      = rcast<const int32_t*>(v.index + 1);
		const auto ends       = rcast<const pos_t*>(jumps + v.index->num_jumps);
		const auto starts     = rcast<const pos_t*>(ends + num_blocks);

		for (std::decay_t<decltype(num_blocks)> i_block{}; i_block < num_blocks; ++i_block) {
			ret[i_itv] = interval_t::from_dna0(k.chrom, starts[i_block] * _res, ends[i_block] * _res, k.other, _refg);
			++i_itv;
		}
	}
	return ret;
}

///////////////////////////////////////////////////////////////

template <typename D, typename S>
void generic_fdict_init(D* dst, const S* src, int size, int& nval)
{
	// Validate the 'src' dict before copying it over.
	// Rules for dicts:
	//   The -inf and +inf values are forbidden and will raise an error.
	//   The nan value is allowed, but it must appear at the end of the dict.
	//   For example, in a dict of size 256, the nan (if any) must appear
	//   in slot [255]. The two valid dicts of size 256 are therefore:
	//
	//       { ... 256 finite values ... } or
	//       { ... 255 finite values ..., nan }
	//
	//   Quantization to the dict is done by std::lower_bound (binary search)
	//   followed by rounding to the nearest value in the dictionary.
	//   All finite values must therefore appear in non-increasing order.
	//
	for (int i = 0; i < size; ++i) {
		GK_CHECK2(!is_inf(src[i]), value, "Dictionary cannot contain inf");
		GK_CHECK2(!is_nan(src[i]) || i == size-1, value, "Only final entry of dictionary can be nan");
		if (i > 0 && !is_nan(src[i]))
			GK_CHECK2(src[i-1] <= src[i], value, "Dictionary must be in non-decreasing order of value");
	}
	for (int i = 0; i < size; ++i)
		dst[i] = D(as_float(src[i]));

	nval = is_nan(src[size-1]) ? size-1 : size;
}

genome_track::float_dict::float_dict()
{
	for (int i = 0; i < 256; ++i) {
		h[i] = nanh();
		f[i] = nanf();
	}
}

void genome_track::float_dict::init(const half_t* dict, int size) { generic_fdict_init(h, dict, size, nval); generic_fdict_init(f, dict, size, nval); }
void genome_track::float_dict::init(const float*  dict, int size) { generic_fdict_init(h, dict, size, nval); generic_fdict_init(f, dict, size, nval); }

uint8_t genome_track::float_dict::encode(half_t x) const { return encode(as_float(x)); }
uint8_t genome_track::float_dict::encode(float  x) const
{
	GK_CHECK2(nval > 0, runtime, "Dictionary uninitialized. (Forgot to call set_dict?)");

	// It's probably faster to do binary search on 32-bit float (hardware support)
	// using f[...] than it is to do binary search on half-float (emulation) using hdict.
	// Use hdict only for decoding, not encoding.
	const float* dict = this->f;
	const float* first = dict;
	const float* last  = dict+nval;
	if (is_nan(x)) {
		bool has_nan = (nval & 1) != 0; // If nval is odd, last slot must be nan
		GK_CHECK2(has_nan, value, "Cannot encode nan to a dict with no nan entry");
		return nval;  // n is the index the the nan in this case
	}
	GK_CHECK2(!is_inf(x), value, "Can only encode finite values into a dict");
	GK_CHECK(x >= *first,    value, "Value {} was less than smallest dictionary value {} (wrong dict? use set_clamping?)", x, *first);
	GK_CHECK(x <= *(last-1), value, "Value {} was larger than largest dictionary value {} (wrong dict? use set_clamping?)", x, *(last-1));

	// Binary search.
	// By using lower_bound instead of upper_bound, it guarantees that
	// if x==*(last-1) then we'll get back i==last-1 and not last
	const float* i = std::lower_bound(first, last, x);
	GK_DBASSERT(i < last);

	// Round to nearest element between i-1 and i.
	// Check 'i > first' to ensure we don't risk *(i-1) evaluating to nan.
	if (i > first && x - *(i-1) < *i - x)
		--i;
	return (uint8_t)(i - dict); // Return index of dict element.
}

END_NAMESPACE_GK
