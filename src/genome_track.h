/*
Copyright (C) 2016-2023 Deep Genomics Inc. All Rights Reserved.
*/
#pragma once
#ifndef __GENOME_KIT_GENOME_TRACK__
#define __GENOME_KIT_GENOME_TRACK__

#include "file.h"
#include "half.h"
#include "interval.h"
#include "util.h"
#include <concepts>
#include <cfloat>
#include <climits>
#include <map>
#include <optional>
#include <string>
#include <vector>

BEGIN_NAMESPACE_GK

using std::string;
using std::map;

class genome_t;

//////////////////////////////////////////////////////////

namespace detail {
	// The same value represented as multiple types.
	// Used for default_value, max_value, min_value.
	// 8 bytes total.
	struct any_t {
		uint8_t u{};       // decoded uint8_t value
		int8_t  i{};       // decoded int8_t value
		half_t  h{nanh()}; // decoded half_t value
		float   f{nanf()}; // decoded float value
		any_t() = default;
		any_t(bool    x): u((uint8_t)x), i((int8_t)x), h((uint8_t)x), f((uint8_t)x) { }
		any_t(uint8_t x): u(x), i((int8_t)x), h(x), f(x) { }
		any_t(int8_t  x): u((uint8_t)x), i(x), h(x), f(x) { }
		any_t(half_t  x): u((uint8_t)gk::as_float(x)), i((int8_t)gk::as_float(x)), h(x), f(gk::as_float(x)) { }
		any_t(float   x): u((uint8_t)x), i((int8_t)x), h(x), f(x) { }
		any_t(uint8_t u, int8_t i, half_t h, float f): u(u), i(i), h(h), f(f) { }
		template <typename T> T as() const;
	};

	template <> INLINE bool    any_t::as<bool>()    const { return u != 0; }
	template <> INLINE uint8_t any_t::as<uint8_t>() const { return u; }
	template <> INLINE int8_t  any_t::as<int8_t>()  const { return i; }
	template <> INLINE half_t  any_t::as<half_t>()  const { return h; }
	template <> INLINE float   any_t::as<float>()   const { return f; }


	template <int dir> struct fractional_store { };

	// When storing a decoded value in forward direction, we simply store using [i]
	template <>
	struct fractional_store<1> {
		template <typename T>
		INLINE static void apply(T* RESTRICT& dst, int i, int& k, int dim, T value) { dst[i] = value; }
	};

	// When storing a decoded value in reverse direction, we must step dst backwards each time
	// we've filled the last position in
	template <>
	struct fractional_store<-1> {
		template <typename T>
		INLINE static void apply(T* RESTRICT& dst, int i, int& k, int dim, T value)
		{
			*(dst++) = value;

			// If we've just written the last dimension of this position,
			// then advance to the first dimension of the next  position
			// to be written to, which is backwards, not forwards.
			if (--k == 0) {
				dst -= 2*dim;
				k = dim;
			}
		}
	};

	template <typename T> struct _sizeof       { enum { size = sizeof(T) }; };
	template <>           struct _sizeof<void> { enum { size = 0 }; };
}

//////////////////////////////////////////////////////////

class genome_track {
public:
	// Encoding formats.
	enum etype_t : int8_t { // DECODED TYPE(S)            DECODED RANGE       ENCODED TYPE
		m0, 				// float16/32, uint8, bool    {1}                 0-bit (constant 1)
		u1, 				// float16/32, uint8, bool    {0..1}              1-bit
		u2, 				// float16/32, uint8          {0..3}              2-bit
		u3, 				// float16/32, uint8          {0..7}              3-bit
		u4, 				// float16/32, uint8          {0..15}             4-bit
		u5, 				// float16/32, uint8          {0..31}             5-bit
		u6, 				// float16/32, uint8          {0..63}             6-bit
		u8, 				// float16/32, uint8          {0..255}            8-bit
		i8, 				// float16/32, int8           {-128..127}         8-bit
		f2, 				// float16/32                 defaults to [0,1]   2-bit index into arbitrary float lookup table (possible non-linearly spaced)
		f3, 				// float16/32                 defaults to [0,1]   3-bit index into arbitrary float lookup table (possible non-linearly spaced)
		f4, 				// float16/32                 defaults to [0,1]   4-bit index into arbitrary float lookup table (possible non-linearly spaced)
		f5, 				// float16/32                 defaults to [0,1]   5-bit index into arbitrary float lookup table (possible non-linearly spaced)
		f6, 				// float16/32                 defaults to [0,1]   6-bit index into arbitrary float lookup table (possible non-linearly spaced)
		f8, 				// float16/32                 defaults to [0,1]   8-bit index into arbitrary float lookup table (possible non-linearly spaced)
		f16,				// float16/32                 arbitrary           16-bit arbitrary half-precision (IEEE 754)
		f32,				// float32                    arbitrary           32-bit arbitrary full-precision (IEEE 754)
		// NOTE: if you alter etype_t, be sure to update etype_as_cstr, etype_default_dtype and be sure to use GK_BEGIN/END_ENCODING
		num_etype
	};
	// Decoding formats.
	enum dtype_t {
		bool_,   // bool
		uint8,   // uint8_t
		int8,    // int8_t
		float16, // half_t
		float32, // float
		// NOTE: if you alter dtype_t, be sure to update GK_BEGIN_ENCODING, dtype_size, dtype_as_cstr, as_dtype,
		//       g_reverse_track_data, any_t, py_dtypes, genome_track::builder::set_data, PyGenomeTrackBuilder::set_data
		num_dtype
	};

	static const char* etype_as_cstr[num_etype];
	static dtype_t     etype_default_dtype[num_etype];
	static const char* dtype_as_cstr[num_dtype];
	static int         dtype_size[num_dtype];
	static dtype_t     as_dtype(const char* name);
	static etype_t     as_etype(const char* name);

	void set_source(string sourcefile);
	void open();
	void close();
	INLINE bool is_open()     const { return _fmap.is_open(); }
	INLINE void ensure_open() const { if (!is_open()) open_on_demand(); }  // Fast check if file / indices already loaded
	INLINE const string& source() const { return _sourcefile; }

	// Decode a genomic interval into pre-allocated destination array `dst`.
	// The destination array must contain c.size() * dim elements.
	//
	// TODO: add stride to the destination writes, so that user can
	// create numpy array for mini-batch and fill data elements directly
	// regardless of how the features are strided in mini-batch memory.
	void operator()(const interval_t& c, bool*    dst) const;
	void operator()(const interval_t& c, int8_t*  dst) const;
	void operator()(const interval_t& c, uint8_t* dst) const;
	void operator()(const interval_t& c, half_t*  dst) const;
	void operator()(const interval_t& c, float*   dst) const;
	void operator()(const interval_t& c, void*    dst, dtype_t dtype) const;

	etype_t etype() const;
	dtype_t dtype() const;
	bool    supports_dtype(dtype_t dtype) const;
	bool    empty() const noexcept;
	INLINE refg_t refg()  const { ensure_open(); return _refg; }
	INLINE int    dim()   const { ensure_open(); return _dim; }
	INLINE int    res()   const { ensure_open(); return _res; }
	INLINE bool   stranded()   const { ensure_open(); return _stranded; }

	static int gtrack_version();

private:

	using any_t = detail::any_t;

	// Fractional encodings packed into this data type.
	// Warning: changing this hasn't been tested.
	// Warning: making this 64-bit may not work directly on big-endian,
	//          due to the way data is encoded into 32-bit offset bits.
	using dword_t = uint32_t;

	// Stored at each index_offset of a GTRACK file is the following header,
	// followed by parallel starts/ends/data_offsets arrays comprising the index itself.
	struct track_index_t {
		int32_t num_blocks{};  // Number of blocks within a particular chromosome-strand track.
		int32_t num_jumps{};   // Number of block indices stored in the jumps array.
		uint8_t is_32bit{};    // Is the index's data_offsets array 32-bit
		uint8_t reserved[3]{};
	}; // ... parallel arrays after here.

	enum { jump_size = 8192 };

	// The in-memory version of header_track_info_t.
	// When a header is loaded, the offsets therein are converted
	// into pointers for simpler access into the memory mapped file.
	struct track_info_t {
		const void*          data{};  // Pointer to the start of data for a particular chromosome-strand track.
		const track_index_t* index{}; // Pointer offset to the per-chromosome, per-strand track data index
	};

	// When file is on disk, this is the structure containing the info for initializing track_info_t
	struct header_track_info_t {
		uint64_t data_offset{};   // File offset to the start of data for a particular chromosome-strand track.
		uint64_t index_offset{};  // File offset to the per-chromosome, per-strand track data index
		chrom_t  chrom{};
		strand_t strand{};
	};

	// GTRACK dictionary (lookup table)
	// Lookup table for fast decoding of values.
	// float value f[i] is always equal to half_t value h[i].
	struct float_dict {
		int32_t nval{};  // Number of non-NaN entries, so h[0:nval] and f[0:nval] are numeric.
		half_t  h[256];  //  512 byte dict for decoding to half_t
		float   f[256];  // 1024 byte dict for decoding to float

		float_dict();

		// Initialize the dict with 'size' (16 or 256) elements
		void init(const half_t* h, int size);
		void init(const float*  f, int size);

		// Encode a half_t/float value into a dictionary index
		uint8_t encode(half_t x) const;
		uint8_t encode(float  x) const;

		// Returns h or f depending on T=half_t or T=float
		template <typename T> INLINE const T* get() const;

        INLINE float min() const;
        INLINE float max() const;
	};

	// GTRACK file header
	struct header_t {
		// File signature
		uint16_t sig{};    		// GTRACK signature (0x70ac)
		uint16_t ver{};    		// GTRACK version number

		// Track properties
		int32_t  dim{};          // Dimensionality of the track
		int32_t  res{};          // Resolution of the track (1pb, 25bp, etc)
		any_t    default_value;  // Default value, used to fill all positions that do not contain data.
		any_t    max_value;      // Maximum observed track datum
		any_t    min_value;      // Minimum observed track datum
		uint8_t  stranded{};     // Whether the track has separate data for each strand
		etype_t  etype{};        // Encoding type of the track data
		uint8_t  reserved[2]{};  // Reserved bytes; serves to align next entry at 8-byte
		refg_t   refg{};         // Reference genome to which this track is aligned

		// Dictionary for decoding track data by lookup
		float_dict dict;

		uint64_t num_tracks;

		// On disk, the header is followed by locations of data and index within file,
		// for each chromosome and strand.
		// header_track_info_t tracks[num_tracks]
	};

	// Fast way to ensure the file has actually been
	// opened and mapped to memory, before trying to
	// access any data / fields from the file.
	// Useful for opening the file on-demand within
	// public API functions, allowing the user to not
	// explicitly open the file ahead of time.
	void open_on_demand() const;

	friend class builder;

	// An 'encoding' object stores the pointers and fields
	// necessary to do fast encoding/decoding for a particular
	// etype and track dimension.
	//
	// NOTE: because of the weird way the encoding object does
	//       polymorphism, subclasses CANNOT contain new data
	//       members, and must have empty default constructors.
	//       For example, see how encoders are initialized by
	//       genome_track::builder constructor.
	//
	struct encoding {

		// Function signatures for encoding, decoding, and default filling.
		//
		//   encode_fn:
		//      dst: Buffer to encode into. Must be num_encoded_bytes(size, dim) bytes.
		//      src: Buffer to encode from. Must contain size*dim data values.
		//      dict: Dictionary for encoding.
		//      size: number of positions in src.
		//      dim: number of values per position.
		//
		//   decode_fn:
		//      dst: Buffer to decode into.
		//      src: Buffer to decode from.
		//      dict: Dictionary for decoding.
		//      size: Number of positions to decode.
		//      dim: Number of values per position.
		//      d: Position (within dst) to start decoding into, i.e. start at dst[dim*d]
		//      s: Position (within src) to start decoding from
		//
		//   dfill_fn:
		//      dst: Buffer to fill.
		//      fill: Value to fill with.
		//      size: Number of positions to fill.
		//      dim: Number of values per position.
		//      d: Position (within dst) to start filling into, i.e. start at dst[dim*d]
		//
		//   expand_fn:
		//      dst: Buffer with coarse data, ready to be expanded to fill the whole buffer.
		//      size: Number of positions to fill in the whole dst buffer.
		//      dim: Number of values per position.
		//      s: Number of coarse values stored in 'dst', i.e. the index after the rightmost decoded value.
		//      res: Resolution of the track, i.e. number of times to repeat each decoded value.
		//      phase: Phase (modulo resolution) of dst's end position with respect
		//             to genomic coordinates.
		//
		using encode_fn = void (*)(void *__restrict, const void *__restrict, const float_dict &, int, int);
		using decode_fn = int (*)(void *__restrict, const void *__restrict, const void *__restrict, int, int, int, int);
		using dfill_fn = int (*)(void *__restrict, const any_t &, int, int, int);
		using expand_fn = void (*)(void *__restrict, int, int, int, int, int);

		// Function pointers to implement encode, decode, default fill, and expand.
		// Each pointer is specialized by etype, dtype, and possibly dim and strand.
		// The pointers are initialized by init().
		encode_fn encoders[num_dtype];
		decode_fn decoders[num_dtype][2];  // [neg_strand] = reverse, [pos_strand] = forward
		dfill_fn  dfillers[num_dtype][2];  // [neg_strand] = reverse, [pos_strand] = forward
		expand_fn expanders[num_dtype];

		// Parameters of this particular encoding.
		int bits_per_encoded_datum;  // 1-8 or 16
		int bytes_per_encoded_word;  // 1 or 4
		etype_t etype;
		any_t default_value;
		any_t range_min;
		any_t range_max;
		float_dict dict;

		void init(etype_t etype, int dim, int res, any_t default_value);
		void init_dict();

		INLINE size_t num_required_bits(int size, int dim) const
		{
			return bits_per_encoded_datum * size * dim;
		}

		INLINE size_t num_encoded_bytes(int size, int dim) const
		{
			if (bits_per_encoded_datum == 0)
				return 0;
			// A datum encoding is not allowed to cross the word boundary to simplify the encoders/decoders
			const int num_per_word = 8*bytes_per_encoded_word/bits_per_encoded_datum; // Round down.
			return size_t(bytes_per_encoded_word)*divup(size*dim, num_per_word); // Round up to the nearest word size.
		}

		// Encodes by applying encoder::apply repeatedly to each src element and storing the encoded values in dst.
		// (We don't bother with all that unrolling / specialization nonsense for encoding, doesn't have to be fast.)
		template <typename encoder>
		INLINE static void generic_encode(      typename encoder::dst_t* RESTRICT dst,
		                                  const typename encoder::src_t* RESTRICT src,
		                                  const float_dict& dict,
		                                  int size, int dim)
		{
			// generic_decode works when e->bits_per_encoded_datum us a multiples of 8, because
			// then it's safe to simply advance the pointers (no fractional dst bytes to worry about).
			for (int i = 0; i < size; ++i)
				for (int j = 0; j < dim; ++j)
					dst[i*dim+j] = encoder::apply(src[i*dim+j], dict);
		}

		// Decodes by applying decoder::apply repeatedly to each src element and storing the result in dst.
		// Important to INLINE this function so that constant propagation of `dim` unrolls the inner loop.
		template <typename decoder, int dir, int unroll>
		INLINE static int generic_decode(       typename decoder::dst_t* RESTRICT dst,
		                                  const typename decoder::src_t* RESTRICT src,
		                                  const typename decoder::dst_t* RESTRICT dict,
		                                  int size, int dim, int d, int s)
		{
			// generic_decode works when e->bits_per_encoded_datum us a multiples of 8, because
			// then it's safe to simply advance the pointers (no fractional src bytes to worry about).
			GK_ASSERT(size > 0);
			GK_ASSERT(dim > 0);
			GK_ASSERT(dim % unroll == 0);  // Ensure that unroll divides dim
			dst += d*dim;
			src += s*dim;
			for (int i = 0; i < size; ++i)
				for (int j = 0; j < dim; j += unroll)
					for (int k = 0; k < unroll; ++k)
						dst[i*dim*dir+j+k] = decoder::apply(src[i*dim+j+k], dict);
			return size*dir;
		}

		template <typename dst_t, int dir, int unroll>
		INLINE static int default_fill(dst_t* RESTRICT dst, const any_t& fill, int size, int dim, int d)
		{
			GK_ASSERT(size > 0);
			GK_ASSERT(dim > 0);
			GK_ASSERT(dim % unroll == 0);  // Ensure that unroll divides dim
			dst += d*dim;
			dst_t x = fill.as<dst_t>();
			for (int i = 0; i < size; ++i)
				for (int j = 0; j < dim; j += unroll)
					for (int k = 0; k < unroll; ++k)
						dst[i*dim*dir+j+k] = x;
			return size*dir;
		}

		// Define some specializations where constant-propagation of 'dim' can reach the inner loops of the generic implementation
		template <typename decoder, int dir, int unroll, int const_dim>
		static int generic_decode_dim(      typename decoder::dst_t* RESTRICT dst,
		                              const typename decoder::src_t* RESTRICT src,
		                              const typename decoder::dst_t* RESTRICT dict,
		                              int size, int dim, int d, int s)
		{
			GK_DBASSERT(dim == const_dim);
			return generic_decode<decoder, dir, unroll>(dst, src, dict, size, const_dim, d, s); // inline generic_decode with 'dim' constant-propagated into the inner loop
		}

		template <typename dst_t, int dir, int unroll, int const_dim>
		INLINE static int default_fill_dim(dst_t* RESTRICT dst, const any_t& fill, int size, int dim, int d)
		{
			GK_DBASSERT(dim == const_dim);
			return default_fill<dst_t, dir, unroll>(dst, fill, size, const_dim, d); // inline default_fill with 'dim' constant-propagated into the inner loop
		}

		template <typename decoder, int dir>
		static decode_fn specialized_decode_fn_dir(int dim)
		{
			// Specialize based on *exact* value of 'dim', i.e. inner loop completely unrolled at compile time.
			// Clang/GCC/MSVC isn't able/willing to unroll this for an extra 6.34% boost
			switch (dim) {
			case 1:  return (decode_fn)generic_decode_dim<decoder, dir, 1, 1>; // unroll = dim = 1
			case 2:  return (decode_fn)generic_decode_dim<decoder, dir, 2, 2>; // unroll = dim = 2
			case 3:  return (decode_fn)generic_decode_dim<decoder, dir, 3, 3>; // unroll = dim = 3
			case 4:  return (decode_fn)generic_decode_dim<decoder, dir, 4, 4>; // unroll = dim = 4
			}

			// Specialize based on *divisor* of 'dim', i.e. inner loop can be partially unrolled at compile time.
			if (dim % 5 == 0) return (decode_fn)generic_decode<decoder, dir, 5>; // unroll = 5
			if (dim % 4 == 0) return (decode_fn)generic_decode<decoder, dir, 4>; // unroll = 4
			if (dim % 3 == 0) return (decode_fn)generic_decode<decoder, dir, 3>; // unroll = 3
			if (dim % 2 == 0) return (decode_fn)generic_decode<decoder, dir, 2>; // unroll = 2
			else              return (decode_fn)generic_decode<decoder, dir, 1>; // unroll = 1, i.e. generic version
		}

		template <typename decoder>
		static decode_fn specialized_decode_fn(int dim, int dir)
		{
			return dir == 1 ? specialized_decode_fn_dir<decoder,  1>(dim)
			                : specialized_decode_fn_dir<decoder, -1>(dim);
		}

		template <typename dst_t, int dir>
		static dfill_fn specialized_default_fill_fn_dir(int dim)
		{
			// Specialize based on *exact* value of 'dim', i.e. inner loop completely unrolled at compile time.
			switch (dim) {
			case 1:  return (dfill_fn)default_fill_dim<dst_t, dir, 1, 1>; // unroll = dim = 1
			case 2:  return (dfill_fn)default_fill_dim<dst_t, dir, 2, 2>; // unroll = dim = 2
			case 3:  return (dfill_fn)default_fill_dim<dst_t, dir, 3, 3>; // unroll = dim = 3
			case 4:  return (dfill_fn)default_fill_dim<dst_t, dir, 4, 4>; // unroll = dim = 4
			}

			// Specialize based on *divisor* of 'dim', i.e. inner loop can be partially unrolled at compile time.
			if (dim % 5 == 0) return (dfill_fn)default_fill<dst_t, dir, 5>; // unroll = 5
			if (dim % 4 == 0) return (dfill_fn)default_fill<dst_t, dir, 4>; // unroll = 4
			if (dim % 3 == 0) return (dfill_fn)default_fill<dst_t, dir, 3>; // unroll = 3
			if (dim % 2 == 0) return (dfill_fn)default_fill<dst_t, dir, 2>; // unroll = 2
			else              return (dfill_fn)default_fill<dst_t, dir, 1>; // unroll = 1, i.e. generic version
		}

		template <typename dst_t>
		static dfill_fn specialized_default_fill_fn(int dim, int dir)
		{
			return dir == 1 ? specialized_default_fill_fn_dir<dst_t,  1>(dim)
			                : specialized_default_fill_fn_dir<dst_t, -1>(dim);
		}

		// See the documentation at end of genome_track::operator() for explanation of how this works.
		// Important to INLINE this function so that constant propagation of `dim` optimizes away the 'dim' loops.
		template <typename T, int unroll>
		INLINE static void generic_expand(T* RESTRICT dst, int size, int dim, int s, int res, int phase)
		{
			// restrict only tells the compiler the pointer will not alias
			// another restrict pointer, so to tell the compiler that we're
			// copying from non-overlapping regions, we need to pass dst as
			// separate pointers.

			#pragma GCC diagnostic push
			#pragma GCC diagnostic ignored "-Wrestrict"
			generic_expand_inner<T, unroll>(dst, dst, size, dim, s, res, phase);
			#pragma GCC diagnostic pop
		}
		template <typename T, int unroll>
		INLINE static void generic_expand_inner(T* RESTRICT dst, T* RESTRICT src, int size, int dim, int s, int res,
												int phase)
		{
			GK_ASSERT(size > 0);
			GK_ASSERT(res > 1);
			GK_ASSERT(dim > 0);
			GK_ASSERT(res % unroll == 0);  // Ensure that unroll divides res
			GK_ASSERT(phase >= 0 && phase < res);

			// Start moving items to the end of the array, so our new destination 'd' will be the end.
			int d = size;

			// First fill any residual slots at the end of the array, to bring us into phase
			// with the track resolution. For this non-performance-critical loop, innermost
			// loop is over dim, which makes the loop stopping criterion much simpler.
			if (phase) {
				--s;      // Advance `s` by 1 position, only once.
				do {
					// Advance `d` by 1 position until we have reached a phase of 0.
					// If we reach `s` before the phase, then the query interval resided
					// entirely within one span (i.e. all the same value).
					if (--d == 0)  // <0 to handle case when size==0
						return;

					// Otherwise copy all dimensions for this particular destination index `d`.
					for (int j = dim-1; j >= 0; --j)  // High-to-low access pattern, for consistency.
						dst[d*dim+j] = src[s*dim+j];

				} while (--phase);
			}

			GK_ASSERT(s >= 0);
			GK_ASSERT(s <= d);

			// Main performance-critical loop, where all expanded runs are full-length and
			// the inner loop is unrolled.
			for (;;) {
				s -= 1;    // Advance `s` by 1 position.
				d -= res;  // Advance `d` by res positions.

				// If we've stepped past the source data, it means we've reached the last chunk
				// of a query interval that was out of phase with the track resolution, and so we
				// can't unroll the loop anymore.
				if (d < 0) {
					d += res;      // Push 'd' back to where it was, after 's'.
					GK_DBASSERT(s <= 0);
					GK_DBASSERT(d >= 0);
					if (d > 1) {
						while (--d)    // Advance `d` by one position at a time until we hit `s`.
							for (int j = dim-1; j >= 0; --j)
								dst[d*dim+j] = src[j];
					}
					return;
				}

				GK_DBASSERT(s >= 0);

				// Otherwise, we're expanding somewhere safe in the middle of the data,
				// and can unroll the 'res' loop. Outer loop is over dimensions so that
				// we can read each value once and write it multiple times.
				for (int j = dim-1; j >= 0; --j) {
					const T x = src[s*dim+j];  // Read the value for dimension j that we need to repeat.

					// Inner loop is over resolution. Go from large addresses to smaller
					// addresses to keep the memory access pattern in a consistent direction.
					for (int r = res-1; r >= 0; r -= unroll)
						for (int k = 0; k < unroll; ++k)
							dst[(d+r-k)*dim+j] = x;
				}
			}
		}

		// Define a specialization so that constant-propagation of 'dim=1' can optimize away the 'dim' loop of generic_expand().
		template <typename T, int unroll, int const_dim>
		static void generic_expand_dim(T* RESTRICT dst, int size, int dim, int s, int res, int phase)
		{
			GK_DBASSERT(dim == const_dim);
			generic_expand<T, unroll>(dst, size, const_dim, s, res, phase); // inline generic_expand with 'dim' constant-propagated into the inner loop
		}

		template <typename T, int unroll>
		static expand_fn specialized_expand_fn_res(int dim)
		{
			// If 1-dimensional track, then specialize on dim to optimize away the dim loop.
			// Since dim isn't inner-most loop, don't bother specializing on other values.
			return (dim == 1) ? (expand_fn)generic_expand_dim<T, unroll, 1>  // dim = 1 at compile time
			                  : (expand_fn)generic_expand<T, unroll>;        // generic dim
		}

		template <typename T>
		static expand_fn specialized_expand_fn(int dim, int res)
		{
			if (res == 0) return nullptr; // Possible if not set yet by genome_track::builder

			// Specialize based on *divisor* of 'res', i.e. inner loop can be partially unrolled at compile time.
			if (res % 5 == 0) return specialized_expand_fn_res<T, 5>(dim); // unroll = 5
			if (res % 4 == 0) return specialized_expand_fn_res<T, 4>(dim); // unroll = 4
			if (res % 3 == 0) return specialized_expand_fn_res<T, 3>(dim); // unroll = 3
			if (res % 2 == 0) return specialized_expand_fn_res<T, 2>(dim); // unroll = 2
			else              return specialized_expand_fn_res<T, 1>(dim); // unroll = 1, i.e. generic resolution
		}

		template <typename decoder, int dir>
		static int decode_m0(      typename decoder::dst_t* RESTRICT dst,
		                     const typename decoder::src_t* RESTRICT src,
		                     const typename decoder::dst_t* RESTRICT dict,
		                     int size, int dim, int d, int s)
		{
			dst += d;
			for (int i = 0; i < size; ++i)
				dst[i*dir] = 1;  // Fill the dst with 1s, ignoring 'src' and which is void* for m0 masks
			return size*dir;
		}

		template <typename decoder>
		static decode_fn specialized_decode_m0_fn(int dim, int dir)
		{
			GK_CHECK(dim == 1, value, "Masks can only be 1 dimensional");
			return dir == 1 ? (decode_fn)decode_m0<decoder, 1> : (decode_fn)decode_m0<decoder, -1>;
		}

		// Encoder for etypes that use 1, 2, or 4 bits.
		template <typename encoder>
		INLINE static void fractional_encode(      typename encoder::dst_t* RESTRICT dst,
		                                     const typename encoder::src_t* RESTRICT src,
		                                     const float_dict& dict,
		                                     int size, int dim)
		{
			const int nbits = encoder::nbits;
			const int num_per_dword = 8*sizeof(dword_t)/nbits; // Round down.

			// Pre-initialize all dwords to zero. We're guaranteed to have this many dwords to work with.
			for (int i = 0; i < divup(size*dim, num_per_dword); ++i)
				dst[i] = 0;

			// Store the encoded bits into each dword using bitwise OR.
			// Bits are stored in dword bit order, with least significant bits
			// representing lower positions / dimensions.
			// For example, a 4-bit encoded type would pack 8 values into 32-bits as follows
			//
			//  [  upper 16 bits  ] [  lower 16 bits  ]
			//  abcd ef98 7654 3210 abcd ef98 7654 3210  <-- bit index, least-significant on right
			//  7777 6666 5555 4444 3333 2222 1111 0000  <-- index of element associated with each bit
			//
			// This facilitates fast decoding via incremental small shifts, since the low-order
			// bits are the first ones that need to be decoded.
			//
			for (int i = 0; i < size; ++i)
				for (int j = 0; j < dim; ++j)
					dst[(i*dim+j)/num_per_dword] |= encoder::apply(src[i*dim+j], dict) << (nbits*((i*dim+j) % num_per_dword));
		}


		template <typename decoder, int dir>
		INLINE static int fractional_decode(      typename decoder::dst_t* RESTRICT dst,
		                                    const typename decoder::src_t* RESTRICT src,
		                                    const typename decoder::dst_t* RESTRICT dict,
		                                    int size, int dim, int d, int s)
		{
			const int nbits = decoder::nbits;
			const int num_per_dword = 8*sizeof(dword_t)/nbits; // Round down.
			const unsigned mask = (1 << nbits)-1;

			int m = size*dim;
			int a = dim*s;
			int b = a + m;
			int first_dword = udivdn(a, num_per_dword);  // round down
			int last_dword  = udivup(b, num_per_dword);  // round up
			int num_dwords  = last_dword - first_dword;

			// Advance pointers to starting positions
			dst += d*dim;
			src += first_dword;

			if (num_dwords <= 1) {

				// All values decoded from a single dword
				dword_t dword = *src;
				int ra = a - num_per_dword*first_dword;
				dword >>= ra*nbits;
				for (int i = 0, k = dim; i < size*dim; ++i, dword >>= nbits)
					detail::fractional_store<dir>::apply(dst, i, k, dim, decoder::apply(dword & mask, dict));

			} else {

				// Handle low-order bits in first dword, if any
				int i = 0, k = dim;
				int ra = umod(a, num_per_dword);
				if (ra > 0) {
					dword_t dword = *src++;
					dword >>= ra*nbits;
					for (; i < num_per_dword-ra; ++i, dword >>= nbits)
						detail::fractional_store<dir>::apply(dst, i, k, dim, decoder::apply(dword & mask, dict));
				}

				// Handle chunks that span entire dword
				int rb = umod(b, num_per_dword);
				int n  = m - rb;
				for (; i < n; i += num_per_dword) {
					dword_t dword = *src++;
					for (int j = 0; j < num_per_dword; ++j, dword >>= nbits)
						detail::fractional_store<dir>::apply(dst, i+j, k, dim, decoder::apply(dword & mask, dict));
				}
				// Handle any remaining high-order bits in final dword, if any
				if (i < m) {
					dword_t dword = *src;
					for (; i < m; ++i, dword >>= nbits)
						detail::fractional_store<dir>::apply(dst, i, k, dim, decoder::apply(dword & mask, dict));
				}
			}

			return size*dir;
		}

		// Define some specializations where constant-propagation of 'dim' can reach the inner loops of the generic implementation
		template <typename decoder, int dir, int const_dim>
		static int fractional_decode_dim(      typename decoder::dst_t* RESTRICT dst,
		                                 const typename decoder::src_t* RESTRICT src,
		                                 const typename decoder::dst_t* RESTRICT dict,
		                                 int size, int dim, int d, int s)
		{
			GK_DBASSERT(dim == const_dim);
			return fractional_decode<decoder, dir>(dst, src, dict, size, const_dim, d, s); // inline with 'dim' constant-propagated into the inner loop
		}

		template <typename decoder, int dir>
		static decode_fn specialized_fractional_decode_fn_dir(int dim)
		{
			// Specialize based on *exact* value of 'dim'.
			// This is important for fractional decoder on reverse strand because
			// its inner-most loop uses 'dim' when advancing the dst pointer.
			switch (dim) {
			case 1:  return (decode_fn)fractional_decode_dim<decoder, dir, 1>; // dim = 1
			case 2:  return (decode_fn)fractional_decode_dim<decoder, dir, 2>; // dim = 2
			case 3:  return (decode_fn)fractional_decode_dim<decoder, dir, 3>; // dim = 3
			case 4:  return (decode_fn)fractional_decode_dim<decoder, dir, 4>; // dim = 4
			}

			// No specialization, generic version
			return (decode_fn)fractional_decode<decoder, dir>;
		}

		template <typename decoder>
		static decode_fn specialized_fractional_decode_fn(int dim, int dir)
		{
			return dir == 1 ? specialized_fractional_decode_fn_dir<decoder,  1>(dim)
			                : specialized_fractional_decode_fn_dir<decoder, -1>(dim);
		}


		// TODO: provide straight memcpy specialization for non-reverse case where transcoding is a no-op, e.g. for half_t -> half_t or i8 -> i8
		// TODO: provide reverse memcpy specialization for reverse case where transcoding is a no-op AND dim=1

	}; // encoding


	#define GK_BEGIN_ENCODING(etype, _bits_per_encoded_datum, _encoded_type, _range_min, _range_max) \
		struct etype##_encoding: public encoding { \
			using encoded_type = _encoded_type; \
			enum { c_bits_per_encoded_datum = _bits_per_encoded_datum }; \
			void init(int dim, int res) \
			{ \
				bits_per_encoded_datum = _bits_per_encoded_datum; \
				bytes_per_encoded_word = detail::_sizeof<encoded_type>::size; \
				range_min = any_t _range_min; \
				range_max = any_t _range_max; \
				encoders[bool_  ] = etype##_encoding::bool__encoder::as_encode_fn(dim); \
				encoders[uint8  ] = etype##_encoding::uint8_encoder::as_encode_fn(dim); \
				encoders[int8   ] = etype##_encoding::int8_encoder::as_encode_fn(dim); \
				encoders[float16] = etype##_encoding::float16_encoder::as_encode_fn(dim); \
				encoders[float32] = etype##_encoding::float32_encoder::as_encode_fn(dim); \
				decoders[bool_  ][as_ordinal(neg_strand)] = etype##_encoding::bool__decoder::as_decode_fn(dim, -1); \
				decoders[bool_  ][as_ordinal(pos_strand)] = etype##_encoding::bool__decoder::as_decode_fn(dim,  1); \
				decoders[uint8  ][as_ordinal(neg_strand)] = etype##_encoding::uint8_decoder::as_decode_fn(dim, -1); \
				decoders[uint8  ][as_ordinal(pos_strand)] = etype##_encoding::uint8_decoder::as_decode_fn(dim,  1); \
				decoders[int8   ][as_ordinal(neg_strand)] = etype##_encoding::int8_decoder::as_decode_fn(dim, -1); \
				decoders[int8   ][as_ordinal(pos_strand)] = etype##_encoding::int8_decoder::as_decode_fn(dim,  1); \
				decoders[float16][as_ordinal(neg_strand)] = etype##_encoding::float16_decoder::as_decode_fn(dim, -1); \
				decoders[float16][as_ordinal(pos_strand)] = etype##_encoding::float16_decoder::as_decode_fn(dim,  1); \
				decoders[float32][as_ordinal(neg_strand)] = etype##_encoding::float32_decoder::as_decode_fn(dim, -1); \
				decoders[float32][as_ordinal(pos_strand)] = etype##_encoding::float32_decoder::as_decode_fn(dim,  1); \
				dfillers[bool_  ][as_ordinal(neg_strand)] = specialized_default_fill_fn<bool   >(dim, -1); \
				dfillers[bool_  ][as_ordinal(pos_strand)] = specialized_default_fill_fn<bool   >(dim,  1); \
				dfillers[uint8  ][as_ordinal(neg_strand)] = specialized_default_fill_fn<uint8_t>(dim, -1); \
				dfillers[uint8  ][as_ordinal(pos_strand)] = specialized_default_fill_fn<uint8_t>(dim,  1); \
				dfillers[int8   ][as_ordinal(neg_strand)] = specialized_default_fill_fn<int8_t >(dim, -1); \
				dfillers[int8   ][as_ordinal(pos_strand)] = specialized_default_fill_fn<int8_t >(dim,  1); \
				dfillers[float16][as_ordinal(neg_strand)] = specialized_default_fill_fn<half_t >(dim, -1); \
				dfillers[float16][as_ordinal(pos_strand)] = specialized_default_fill_fn<half_t >(dim,  1); \
				dfillers[float32][as_ordinal(neg_strand)] = specialized_default_fill_fn<float  >(dim, -1); \
				dfillers[float32][as_ordinal(pos_strand)] = specialized_default_fill_fn<float  >(dim,  1); \
				expanders[bool_  ] = specialized_expand_fn<bool   >(dim, res); \
				expanders[uint8  ] = specialized_expand_fn<uint8_t>(dim, res); \
				expanders[int8   ] = specialized_expand_fn<int8_t >(dim, res); \
				expanders[float16] = specialized_expand_fn<half_t >(dim, res); \
				expanders[float32] = specialized_expand_fn<float  >(dim, res); \
				init_dict(); /* static dispatch to base class (no-op) or subclass specialization */ \
			}

	#define GK_END_ENCODING \
		};

	#define _GK_DTYPE_bool_    bool
	#define _GK_DTYPE_uint8    uint8_t
	#define _GK_DTYPE_int8     int8_t
	#define _GK_DTYPE_float16  half_t
	#define _GK_DTYPE_float32  float

	#define GK_BEGIN_ENCODER(_decoded_type) \
		struct _decoded_type##_encoder { \
			enum { nbits = c_bits_per_encoded_datum }; \
			using encoder_type = _decoded_type##_encoder; \
			using dst_t        = encoded_type; \
			using src_t        = _GK_DTYPE_##_decoded_type;

	#define GK_END_ENCODER \
		};

	#define GK_BEGIN_DECODER(_decoded_type) \
		struct _decoded_type##_decoder { \
			enum { nbits = c_bits_per_encoded_datum }; \
			using decoder_type = _decoded_type##_decoder; \
			using dst_t        = _GK_DTYPE_##_decoded_type; \
			using src_t        = encoded_type;

	#define GK_END_DECODER \
		};

	#define GK_NULL_ENCODER(_decoded_type) \
			GK_BEGIN_ENCODER(_decoded_type) \
				INLINE static encode_fn as_encode_fn(int dim) { return nullptr; } \
			GK_END_ENCODER

	#define GK_NULL_DECODER(_decoded_type) \
			GK_BEGIN_DECODER(_decoded_type) \
				INLINE static decode_fn as_decode_fn(int dim, int dir) { return nullptr; } \
			GK_END_DECODER

	#define GK_NULL_TRANSCODER(_decoded_type) \
			GK_NULL_ENCODER(_decoded_type) \
			GK_NULL_DECODER(_decoded_type)

	#define GK_ENCODER_APPLY_RANGED(x, lo, hi) \
			dst_t y = (dst_t)x; \
			GK_CHECK((in_range(x, lo, hi) && src_t(y) == x), value, "Value {} can't be encoded, must be integral value in range [{},{}]", x, lo, hi); \
			return y

	#define GK_ENCODER_APPLY \
		INLINE static dst_t apply(src_t x, const float_dict& dict)

	#define GK_DECODER_APPLY \
		INLINE static dst_t apply(src_t y, const dst_t* RESTRICT dict)

	#define GK_ENCODER_FN \
		INLINE static encode_fn as_encode_fn(int dim)

	#define GK_DECODER_FN \
		INLINE static decode_fn as_decode_fn(int dim, int dir)

	GK_BEGIN_ENCODING(m0, 0, void, (1, 1, 1, 1), (1, 1, 1, 1))
		// No encoders for m0 -- treated as special case where anything in the genome-wide index gets decoded to 1
		GK_NULL_ENCODER(bool_)
		GK_BEGIN_DECODER(bool_)
			GK_DECODER_FN    { return specialized_decode_m0_fn<uint8_decoder>(dim, dir); }
		GK_END_DECODER
		GK_NULL_ENCODER(uint8)
		GK_BEGIN_DECODER(uint8)
			GK_DECODER_FN    { return specialized_decode_m0_fn<uint8_decoder>(dim, dir); }
		GK_END_DECODER
		GK_NULL_ENCODER(int8)
		GK_BEGIN_DECODER(int8)
			GK_DECODER_FN    { return specialized_decode_m0_fn<uint8_decoder>(dim, dir); }  // deliberately share decoder to reduce binary size
		GK_END_DECODER
		GK_NULL_ENCODER(float16)
		GK_BEGIN_DECODER(float16)
			GK_DECODER_FN    { return specialized_decode_m0_fn<decoder_type>(dim, dir); }
		GK_END_DECODER
		GK_NULL_ENCODER(float32)
		GK_BEGIN_DECODER(float32)
			GK_DECODER_FN    { return specialized_decode_m0_fn<decoder_type>(dim, dir); }
		GK_END_DECODER
	GK_END_ENCODING

	// All the encodings have an encoder from float16/32, to support loading of WIG files.

	GK_BEGIN_ENCODING(u1, 1, dword_t, (0, 0, 0, 0), (1, 1, 1, 1))
		GK_BEGIN_ENCODER(uint8)
			GK_ENCODER_APPLY { GK_ENCODER_APPLY_RANGED(x, 0, 1); }
			GK_ENCODER_FN    { return (encode_fn)fractional_encode<encoder_type>; }
		GK_END_ENCODER
		GK_BEGIN_DECODER(uint8)
			GK_DECODER_APPLY { return (dst_t)y; }
			GK_DECODER_FN    { return specialized_fractional_decode_fn<decoder_type>(dim, dir); }
		GK_END_DECODER
		GK_NULL_TRANSCODER(int8)
		GK_BEGIN_ENCODER(bool_)
			GK_ENCODER_APPLY { return x ? 1 : 0; }
			GK_ENCODER_FN    { return (encode_fn)fractional_encode<encoder_type>; }
		GK_END_ENCODER
		GK_BEGIN_DECODER(bool_)
			GK_DECODER_FN    { return uint8_decoder::as_decode_fn(dim, dir); }  // deliberately share decoder to reduce binary size
		GK_END_DECODER
		GK_BEGIN_ENCODER(float16)
			GK_ENCODER_APPLY { GK_ENCODER_APPLY_RANGED(as_float(x), 0, 1); }
			GK_ENCODER_FN    { return (encode_fn)fractional_encode<encoder_type>; }
		GK_END_ENCODER
		GK_BEGIN_DECODER(float16)
			GK_DECODER_APPLY { return dst_t(y); }
			GK_DECODER_FN    { return specialized_fractional_decode_fn<decoder_type>(dim, dir); }
		GK_END_DECODER
		GK_BEGIN_ENCODER(float32)
			GK_ENCODER_APPLY { GK_ENCODER_APPLY_RANGED(x, 0, 1); }
			GK_ENCODER_FN    { return (encode_fn)fractional_encode<encoder_type>; }
		GK_END_ENCODER
		GK_BEGIN_DECODER(float32)
			GK_DECODER_APPLY { return dst_t(y); }
			GK_DECODER_FN    { return specialized_fractional_decode_fn<decoder_type>(dim, dir); }
		GK_END_DECODER
	GK_END_ENCODING

#define GK_DEFINE_FRACTIONAL_UINT_ENCODING(nbits) \
	GK_BEGIN_ENCODING(u##nbits, nbits, dword_t, (0, 0, 0, 0), ((1<<nbits)-1, (1<<nbits)-1, (1<<nbits)-1, (1<<nbits)-1)) \
		GK_NULL_TRANSCODER(bool_) \
		GK_BEGIN_ENCODER(uint8) \
			GK_ENCODER_APPLY { GK_ENCODER_APPLY_RANGED(x, 0, (1<<nbits)-1); } \
			GK_ENCODER_FN    { return (encode_fn)fractional_encode<encoder_type>; } \
		GK_END_ENCODER \
		GK_BEGIN_DECODER(uint8) \
			GK_DECODER_APPLY { return (dst_t)y; } \
			GK_DECODER_FN    { return specialized_fractional_decode_fn<decoder_type>(dim, dir); } \
		GK_END_DECODER \
		GK_NULL_TRANSCODER(int8) \
		GK_BEGIN_ENCODER(float16) \
			GK_ENCODER_APPLY { GK_ENCODER_APPLY_RANGED(as_float(x), 0, (1<<nbits)-1); } \
			GK_ENCODER_FN    { return (encode_fn)fractional_encode<encoder_type>; } \
		GK_END_ENCODER  \
		GK_BEGIN_DECODER(float16) \
			GK_DECODER_APPLY { return dst_t(y); } \
			GK_DECODER_FN    { return specialized_fractional_decode_fn<decoder_type>(dim, dir); } \
		GK_END_DECODER \
		GK_BEGIN_ENCODER(float32) \
			GK_ENCODER_APPLY { GK_ENCODER_APPLY_RANGED(x, 0, (1<<nbits)-1); } \
			GK_ENCODER_FN    { return (encode_fn)fractional_encode<encoder_type>; } \
		GK_END_ENCODER \
		GK_BEGIN_DECODER(float32) \
			GK_DECODER_APPLY { return dst_t(y); } \
			GK_DECODER_FN    { return specialized_fractional_decode_fn<decoder_type>(dim, dir); } \
		GK_END_DECODER \
	GK_END_ENCODING

	GK_DEFINE_FRACTIONAL_UINT_ENCODING(2)
	GK_DEFINE_FRACTIONAL_UINT_ENCODING(3)
	GK_DEFINE_FRACTIONAL_UINT_ENCODING(4)
	GK_DEFINE_FRACTIONAL_UINT_ENCODING(5)
	GK_DEFINE_FRACTIONAL_UINT_ENCODING(6)

	GK_BEGIN_ENCODING(u8, 8, uint8_t, (0, 0, 0, 0), (255, 127, 255, 255))
		GK_NULL_TRANSCODER(bool_)
		GK_BEGIN_ENCODER(uint8)
			GK_ENCODER_APPLY { return x; }
			GK_ENCODER_FN    { return (encode_fn)generic_encode<encoder_type>; }
		GK_END_ENCODER
		GK_BEGIN_DECODER(uint8)  // TODO: use memcpy for non-reverse case
			GK_DECODER_APPLY { return y; }
			GK_DECODER_FN    { return specialized_decode_fn<decoder_type>(dim, dir); }
		GK_END_DECODER
		GK_NULL_TRANSCODER(int8)
		GK_BEGIN_ENCODER(float16)
			GK_ENCODER_APPLY { GK_ENCODER_APPLY_RANGED(as_float(x), 0, 255); }
			GK_ENCODER_FN    { return (encode_fn)generic_encode<encoder_type>; }
		GK_END_ENCODER
		GK_BEGIN_DECODER(float16)
			GK_DECODER_APPLY { return dst_t(y); }
			GK_DECODER_FN    { return specialized_decode_fn<decoder_type>(dim, dir); }
		GK_END_DECODER
		GK_BEGIN_ENCODER(float32)
			GK_ENCODER_APPLY { GK_ENCODER_APPLY_RANGED(x, 0, 255); }
			GK_ENCODER_FN    { return (encode_fn)generic_encode<encoder_type>; }
		GK_END_ENCODER
		GK_BEGIN_DECODER(float32)
			GK_DECODER_APPLY { return dst_t(y); }
			GK_DECODER_FN    { return specialized_decode_fn<decoder_type>(dim, dir); }
		GK_END_DECODER
	GK_END_ENCODING

	GK_BEGIN_ENCODING(i8, 8, int8_t, (0, -128, -128, -128), (127, 127, 127, 127))
		GK_NULL_TRANSCODER(bool_)
		GK_NULL_TRANSCODER(uint8)
		GK_BEGIN_ENCODER(int8)
			GK_ENCODER_FN    { return u8_encoding::uint8_encoder::as_encode_fn(dim); }  // deliberately share encoder to reduce binary size
		GK_END_ENCODER
		GK_BEGIN_DECODER(int8)
			GK_DECODER_FN    { return u8_encoding::uint8_decoder::as_decode_fn(dim, dir); }  // deliberately share decoder to reduce binary size
		GK_END_DECODER
		GK_BEGIN_ENCODER(float16)
			GK_ENCODER_APPLY { GK_ENCODER_APPLY_RANGED(as_float(x), -128, 127); }
			GK_ENCODER_FN    { return (encode_fn)generic_encode<encoder_type>; }
		GK_END_ENCODER
		GK_BEGIN_DECODER(float16)
			GK_DECODER_APPLY { return dst_t(y); }
			GK_DECODER_FN    { return specialized_decode_fn<decoder_type>(dim, dir); }
		GK_END_DECODER
		GK_BEGIN_ENCODER(float32)
			GK_ENCODER_APPLY { GK_ENCODER_APPLY_RANGED(x, -128, 127); }
			GK_ENCODER_FN    { return (encode_fn)generic_encode<encoder_type>; }
		GK_END_ENCODER
		GK_BEGIN_DECODER(float32)
			GK_DECODER_APPLY { return dst_t(y); }
			GK_DECODER_FN    { return specialized_decode_fn<decoder_type>(dim, dir); }
		GK_END_DECODER
	GK_END_ENCODING

#define GK_DEFINE_FRACTIONAL_DICT_ENCODING(nbits) \
	GK_BEGIN_ENCODING(f##nbits, nbits, dword_t, (0, 0, -65504, -65504), (0, 0, 65504, 65504)) /* (0, 0, ...), (0, 0, ...) because no transcoding to uint8 or int8 */ \
		GK_NULL_TRANSCODER(bool_) \
		GK_NULL_TRANSCODER(uint8) \
		GK_NULL_TRANSCODER(int8)  \
		GK_BEGIN_ENCODER(float16) \
			GK_ENCODER_APPLY { return dict.encode(x); } \
			GK_ENCODER_FN    { return (encode_fn)fractional_encode<encoder_type>; } \
		GK_END_ENCODER \
		GK_BEGIN_DECODER(float16) \
			GK_DECODER_APPLY { return dict[y]; } \
			GK_DECODER_FN    { return specialized_fractional_decode_fn<decoder_type>(dim, dir); } \
		GK_END_DECODER \
		GK_BEGIN_ENCODER(float32) \
			GK_ENCODER_APPLY { return dict.encode(x); } \
			GK_ENCODER_FN    { return (encode_fn)fractional_encode<encoder_type>; } \
		GK_END_ENCODER \
		GK_BEGIN_DECODER(float32) \
			GK_DECODER_APPLY { return dict[y]; } \
			GK_DECODER_FN    { return specialized_fractional_decode_fn<decoder_type>(dim, dir); } \
		GK_END_DECODER \
	GK_END_ENCODING

	GK_DEFINE_FRACTIONAL_DICT_ENCODING(2)
	GK_DEFINE_FRACTIONAL_DICT_ENCODING(3)
	GK_DEFINE_FRACTIONAL_DICT_ENCODING(4)
	GK_DEFINE_FRACTIONAL_DICT_ENCODING(5)
	GK_DEFINE_FRACTIONAL_DICT_ENCODING(6)

	GK_BEGIN_ENCODING(f8, 8, uint8_t, (0, 0, -65504, -65504), (0, 0, 65504, 65504)) // (0, 0, ...), (0, 0, ...) because no transcoding to uint8 or int8
		GK_NULL_TRANSCODER(bool_)
		GK_NULL_TRANSCODER(uint8)
		GK_NULL_TRANSCODER(int8)
		GK_BEGIN_ENCODER(float16)
			GK_ENCODER_APPLY { return dict.encode(x); }
			GK_ENCODER_FN    { return (encode_fn)generic_encode<encoder_type>; }
		GK_END_ENCODER
		GK_BEGIN_DECODER(float16)
			GK_DECODER_APPLY { return dict[y]; }
			GK_DECODER_FN    { return specialized_decode_fn<decoder_type>(dim, dir); }
		GK_END_DECODER
		GK_BEGIN_ENCODER(float32)
			GK_ENCODER_APPLY { return dict.encode(x); }
			GK_ENCODER_FN    { return (encode_fn)generic_encode<encoder_type>; }
		GK_END_ENCODER
		GK_BEGIN_DECODER(float32)
			GK_DECODER_APPLY { return dict[y]; }
			GK_DECODER_FN    { return specialized_decode_fn<decoder_type>(dim, dir); }
		GK_END_DECODER
	GK_END_ENCODING

	GK_BEGIN_ENCODING(f16, 16, half_t, (0, 0, -65504, -65504), (0, 0, 65504, 65504)) // (0, 0, ...), (0, 0, ...) because no transcoding to uint8 or int8
		GK_NULL_TRANSCODER(bool_)
		GK_NULL_TRANSCODER(uint8)
		GK_NULL_TRANSCODER(int8)
		GK_BEGIN_ENCODER(float16)
			GK_ENCODER_APPLY { return x; }
			GK_ENCODER_FN    { return (encode_fn)generic_encode<encoder_type>; }
		GK_END_ENCODER
		GK_BEGIN_DECODER(float16)  // TODO: use memcpy for non-reverse case
			GK_DECODER_APPLY { return y; }
			GK_DECODER_FN    { return (decode_fn)specialized_decode_fn<decoder_type>(dim, dir); }
		GK_END_DECODER
		GK_BEGIN_ENCODER(float32)
			GK_ENCODER_APPLY { return half_t(x); }
			GK_ENCODER_FN    { return (encode_fn)generic_encode<encoder_type>; }
		GK_END_ENCODER
		GK_BEGIN_DECODER(float32)
			GK_DECODER_APPLY { return as_float(y); }
			GK_DECODER_FN    { return (decode_fn)specialized_decode_fn<decoder_type>(dim, dir); }
		GK_END_DECODER
	GK_END_ENCODING

	GK_BEGIN_ENCODING(f32, 32, float, (0, 0, 0, -FLT_MAX), (0, 0, 0, FLT_MAX))
		GK_NULL_TRANSCODER(bool_)
		GK_NULL_TRANSCODER(uint8)
		GK_NULL_TRANSCODER(int8)
		GK_NULL_TRANSCODER(float16)
		GK_BEGIN_ENCODER(float32)
			GK_ENCODER_APPLY { return x; }
			GK_ENCODER_FN    { return (encode_fn)generic_encode<encoder_type>; }
		GK_END_ENCODER
		GK_BEGIN_DECODER(float32)
			GK_DECODER_APPLY { return y; }
			GK_DECODER_FN    { return (decode_fn)specialized_decode_fn<decoder_type>(dim, dir); }
		GK_END_DECODER
	GK_END_ENCODING

	mmap_file _fmap;          // Memory mapped file

	// Copies of some header fields that we access frequently, to avoid extra pointer indirection to the memory mapped file header
	int       _dim{};         // Dimensionality of the track (number of tracks bundled together).
	int       _res{};         // Resolution of the track (number of genomic positions each data value should span).
	refg_t    _refg{};   	  // The reference genome for which this is a track
							  // not optional for simplicity. open guaranteed to be called
	bool      _stranded{};    // Whether or not the track is stranded (true) or symmetric with respect to strand (false)

	// Pre-loaded information about the track data/index locations and the encoding
	chrom_map_t<strand_t, track_info_t> _tracks;  // Pointers to the track indices and data blocks in memory
	encoding     _encoding;      // Contains dictionary and callbacks to implement encoding/decoding.
	string       _sourcefile;    // Path to the GTRACK file being mapped

	// TODO: if the binary search through indices is occupying any significant fraction of execution
	//       time (as opposed to, say, decoding data), then consider implementing a a thread-local
	//       (not global) cache of data offsets for a given (query interval, index hash);
	//       of course, this would have to be at a level above any blockwise compression within a contig
	//       (i.e. block compression of is a subsequent level of structure on a contig, after the contig
	//        data offset is identified via the index)
	// TODO: mutable mutex member needed here for open_on_demand to be thread safe


	///////////////////////////////////////////////////////////////////
public:

	enum class strandedness_t : uint8_t { // how data is interpreted with respect to strand
		single_stranded, // data is applied to both strands in Interval coordinate order.
		strand_unaware, // data applied in Interval coordinate order (ignore Interval strand)
		strand_aware // data applied from 5" end to 3" end
	};


	// TODO:
	//   allow builder to be instantiated and run by loading a config file
	//   that specifies encodings, source files (e.g. bigwigs), etc.
	//
	// File format:
	//   Block index goes AFTER data, so that the data can be flushed,
	//   so that data can be streamed out even as
	//

	class builder {
	public:

		// Besides the output file and encoding type, the builder constructor
		// takes one mandatory and two optional arguments:
		//
		//    strandedness = determines the order by which the data is applied
		//    dim = dimensionality of the track (default: 1)
		//    res = resolution of data (default: 1)
		//
		// Resolution allows track data to be specified, and stored, at coarser
		// resolution than 1bp.
		//
		// For example, if a WIG file specifies track with step and span
		// of 5bp:
		//
		//    fixedStep chrom=chr1 start=1 step=5 span=5
		//    0.57
		//    0.23
		//    ...
		//
		// then load_data_from_wig will automatically call set_resolution(5)
		// so that each encoded datum gets repeated 5 times when decoded:
		//
		//    0    1    2    3    4    5    6    7    8    9      <-- coordinate
		//    0.57 0.57 0.57 0.57 0.57 0.23 0.23 0.23 0.23 0.23   <-- decoded value
		//
		// If you're calling set_data on each interval yourself, then
		// you'll want to provide the full genomic interval, but data
		// covering only interval.size() / resolution positions.
		//
		builder(string outfile, etype_t etype, strandedness_t strandedness, const genome_t&, int dim=1, int res=1);

		// Contiguous blocks of encoded data that are equal to default_value
		// can be excluded from the encoded track, thereby reducing file size.
		// For example, conservation tracks often contain large blocks that encode
		// to 0.0 or close to 0.0. By setting min_size=50, any contiguous block of
		// 0.0 longer than 50 will be encoded by omission.
		void set_sparsity(int min_run=48, float min_delta=0.0f);

		// Normally, set_data raises an error if it encounters a value that,
		// after transformation, is still outside the encodable range.
		// For example a value too large to be represented as float16 can trigger
		// an overflow error, or trying to encode "2.0" to a dict with maximum
		// value 1.0 will likewise trigger an error.
		// Calling this first will cause data to be clamped to the encodable range.
		void set_clamping();

		// Override the default [0,1] range dictionary of f2..8 etypes.
		// For f<n> the dict must contain 1<<n values respectively.
		// The last dict entry is allowed to be NaN, if so desired.
		// Must be called BEFORE and data is actually specified.
		void set_dict(const half_t* dict);
		void set_dict(const float*  dict);

		// Set the default value for regions with unspecified data.
		// The default value can be anything representable in the
		// decoded type. The value does not need to be representable
		// by the encoded type (e.g. can be an integer,
		// character, or float value represented or representable by
		// the encoded type itself). For example, if one were encoding
		// DNA as 2-bit values with 0,1,2,3 representing ACGT, one
		// could fill the un-indexed regions of the genome with N
		//       ...NNNNNNacgtacgtacgtNNNNNN...
		// by representing N with another integer, say 4, even though
		// 4 cannot be represented by 2 bits.
		//
		// Setting the default for uint8_t also sets the same 8-bit
		// binary value for int8_t, and vice versa.
		// Setting the default for half_t also sets the default for
		// float, and vice versa.
		void set_default_value(bool    value);
		void set_default_value(int     value);
		void set_default_value(half_t  value);
		void set_default_value(float   value);

		// Set the transformation that all data should undergo before
		// being encoded. Source values in range [a, b] will be affinely
		// transformed to range [c, d].
		//
		// After calling this function, set_data must only be
		// called with either half_t or float src argument,
		// since transforming integer values is not supported.
		void set_transform(float a, float b, float c, float d);

		// Restrict data to the given interval, discarding or cropping away
		// all data outside it.
		// The strand is ignored, i.e. restrictions allow data on both strands.
		// This is useful for making small version of full-sized data pipeline,
		// for the sake of creating unit tests or iterating faster during development.
		// If the restriction interval isn't aligned with the track resolution, some
		// then the restriction interval will be expanded until it's aligned.
		void set_restriction(const interval_t& restriction);

		// Set track data for a given interval.
		// The interval must not overlap any other blocks already added.
		//
		// The 'data' array must contain interval.size() * dim / res
		// values, where dim is the dimensionality passed to the constructor
		// and resolution is the resolution of the track (1 if unspecified).
		//
		// The inner-most stride of data (stride 1) must be over the values
		// for a specific position, and the outer-most stride (stride dim)
		// must be over positions in the interval, i.e. in a C array it would
		// be stored in data[position][dimension] order.
		//
		// The order of positions depends on the strandedness specified
		// in the builder constructor.
		// single_stranded: implies positive strand, and thus data is applied
		// in reference-strand order (also sense order in this case).
		// strand_unaware: ignores the Interval strand, data is applied
		// in Interval coordinate order.
		// strand_aware: data is applied from 5" end to 3" end, order depends
		// on the Interval strand.
		//
		// NOTE: If the provided data is float, it will be converted
		// to IEEE half-precision and if an overflow occurs an exception
		// will be raised. This is the default behaviour because inf
		// values are not permitted in f2..8-based types and because
		// it overloading to +/- inf is likely not the desired behaviour
		// anyway. For an etype of f16, inf is permitted, but it is
		// the user's responsibility to ensure each value is either
		// +/- inf, nan, or finite value within the range of half-precision
		// IEEE floating point.
		// If overflowing or clamping is the desired behaviour, preprocess
		// the values with overflow_for_half() or overflow_for_clamp()
		// before calling set_data.
		void set_data(const interval_t& interval, const bool*    data);
		void set_data(const interval_t& interval, const uint8_t* data);
		void set_data(const interval_t& interval, const int8_t*  data);
		void set_data(const interval_t& interval, const half_t*  data);
		void set_data(const interval_t& interval, const float*   data);
		void set_data(const interval_t& interval, const void*    data, dtype_t dtype);

		// Load all data from a WIG file.
		//
		// This function calls set_data repeatedly on the blocks of data
		// encountered in a WIG file. The resulting track is
		// stranded if two files are specified.
		//
		// The number of columns in the file must match the track dimension.
		// Furthermore, the step/span of the WIG data must match the track resolution.
		// (Exception: the last row of a WIG may have arbitrary step/span, since it
		//  usually represents a BIGWIG span being truncated at end of chromosome)
		//
		void set_data_from_wig(const string& infile); // forces single-stranded mode
		void set_data_from_wig(const string& infile_pos, const string& infile_neg);
		void set_data_from_wig(const char* infile_pos, const char* infile_neg);

		// Load all data from a BEDGRAPH file.
		//
		// Otherwise identical in behaviour to set_data_from_wig
		//
		void set_data_from_bedgraph(const string& infile); // forces single-stranded mode
		void set_data_from_bedgraph(const string& infile_pos, const string& infile_neg);
		void set_data_from_bedgraph(const char* infile_pos, const char* infile_neg);

		// Load all data from a BED file.
		//
		// This function calls set_data repeatedly on the intervals encountered
		// in a BED file. Requires `dim=1`.
		//
		// If the track is unstranded, the strand column in the BED file must either
		// not be specified or it must contain '.'. Otherwise the strand column must
		// be specified and contain either '+' or '-'.
		//
		// If categories are specified, then etype must be `u1..8`.
		//
		// The start/end of each interval must be aligned to the track resolution.
		// (Exception: the last interval may have arbitrary size, since it
		//  usually represents a BIGBED span being truncated at end of chromosome)
		//
		void set_data_from_bed(const string& infile);
		void set_data_from_bed(const string& infile, const vector<string>& categories);
		template <typename T>
		void _set_data_from_bed(const string& infile, const vector<string>& categories);

		// Can be called after, for example, a chromosome has had all its data added.
		// This causes the builder to flush all current data to the file, freeing
		// main memory in situations where the genome-wide data is too large.
		// If you call flush() and then try to set more data on a chromosome
		// that was already flushed, it will be an error.
		void flush();
		void flush(std::optional<chrom_t> chrom, strand_t strand);

		// Finalize the .gtrack file on disk.
		// Call this after all data has been added and the file is ready to be closed.
		// If this function is not called, it will be assumed that an error occurred
		// before the file was finished writing, and the .gtrack file will be deleted
		// by the genome_track_builder destructor.
		void finalize();

		// noop, kept for compatiblity
		void close() {}

		// Counts the number of data blocks (separate intervals) currently in the index.
		size_t data_size() const  { return _data_size; }
		size_t index_size() const { return _index_size; }

		INLINE int     dim()   const { return _h.dim; }
		INLINE int     res()   const { return _h.res; }
		INLINE etype_t etype() const { return _h.etype; }
		INLINE dtype_t dtype() const { return genome_track::etype_default_dtype[etype()]; }
		INLINE refg_t  refg()  const { return _h.refg; }
		INLINE const std::string refg_name()  const { return _refg_name; }
		INLINE bool    stranded() const { return _h.stranded; }
		INLINE int           sparsity_min_run()   const { return _sparsity_min_run; }
		INLINE const any_t&  sparsity_min_delta() const { return _sparsity_min_delta; }
		INLINE const any_t&  default_value() const { return _h.default_value; }
		INLINE const string& filename() const { return _outfile; }
		INLINE bool supports_dtype(dtype_t dtype) const { return _encoding.decoders[dtype][as_ordinal(neg_strand)] != nullptr; }
		INLINE bool finalized() const { return _finalized; }

		const chrom_names_t& chrom_names() const;

	private:
		// Key for the std::map that builds up all the encoded intervals, before they're flushed to disk.
		struct span_t {
			INLINE span_t(pos_t a, pos_t b): a(a), b(b) { }
			INLINE pos_t size() const { return b - a; }
			INLINE bool operator<(const span_t& s) const { return a < s.a; }
			pos_t a, b;
		};

		struct encoded_sizes_t {
			uint64_t bytes{};
			uint64_t bits{};
		};

		class track_info_t {
		private:
			using block_by_span_t = std::map<span_t, std::unique_ptr<uint8_t[]>>;

		public:
			class adder {  // optimize adding sequential subspans
			public:
				explicit adder(track_info_t* ti)
				: _ti(ti)
				, _inserter{std::inserter(_ti->_blocks, std::end(_ti->_blocks))}
				{
				}
				void validate(span_t span);
				void add(span_t span, std::unique_ptr<uint8_t[]>&& data, encoded_sizes_t encoded_sizes);

			private:
				track_info_t*                         _ti;
				std::insert_iterator<block_by_span_t> _inserter;
			};
			friend adder;

			template <std::regular_invocable<int> EncodedSizesFn>
			void flush(binary_file& out_file, EncodedSizesFn&& encoded_sizes);

			uint64_t data_size() const;
			uint64_t index_size(bool encoded_data) const;
			bool     flushed() const;

		private:
			template <class T>
			constexpr bool can_inplace(encoded_sizes_t encoded_sizes)
			{
				// data_bits != 8 * data_bytes since a datum encoding cannot cross the dword
				// boundary
				// need a bit for inplace marker
				return encoded_sizes.bits < sizeof(T) * CHAR_BIT && encoded_sizes.bytes <= sizeof(T);
			}

			bool allow_32bit() const;
			int  num_blocks() const;
			int  num_jumps() const;

			std::map<span_t, std::unique_ptr<uint8_t[]>>
										_blocks;  //  size of block is determined by builder encoding+dims
			uint64_t                    _max_offset32{};
			uint64_t                    _max_offset64{};
			pos_t                       _max_end{-1};  // Init extent to -1 so that idx.num_jumps == 0 if no blocks
			int                         _num_blocks{};
			bool                        _flushed{};
		};

		template <typename T>
		void set_dict_impl(const T* dict);

		void add_track_entry(track_info_t::adder& adder, const span_t& span, std::unique_ptr<uint8_t[]>&& data);

		template <typename T>
		void set_data_impl(const interval_t& interval, const T* data);

		void set_data_from_wig(const char* infile, strand_t strand);
		template <typename T>
		void _set_data_from_wig(const char* infile, strand_t strand);

		void set_data_from_bedgraph(const char* infile, strand_t strand);
		template <typename T>
		void _set_data_from_bedgraph(const char* infile, strand_t strand);

		void wig_bedgraph_config(const char* pos_infile, const char* neg_infile, const char* filetype);

		bool is_default(bool    x) const;
		bool is_default(uint8_t x) const;
		bool is_default(int8_t  x) const;
		bool is_default(half_t  x) const;
		bool is_default(float   x) const;

		encoded_sizes_t encoded_sizes(int length) const;

		void init_refg(refg_t refg);

		string   _outfile;  // The name of file being written to.
		header_t _h;
		any_t    _sparsity_min_delta;   // Minimum delta between default_value and datum to be considered equivalent to default
		int      _sparsity_min_run;     // Minimum run of contiguous default_values to be excluded from the encoded data track.
		strandedness_t _strandedness;   // Controls the order data is applied in.
		bool     _data_transform;       // If true, set_data should perform a linear transformation before encoding.
		bool     _data_clamp;           // If true, set_data should clamp data to range [c, d] of set_transform
		std::function<float(float)> _clamp_float;
		float    _data_scale;
		float    _data_bias;
		float    _data_min;             		// c of set_transform
		float    _data_max;             		// d of set_transform
		size_t   _data_size;            		// Number of positions spanned by encoded data
		size_t   _index_size;           		// Number of blocks that were added, i.e. the size of the index
		std::optional<interval_t> _restriction;	// Restrict all data intervals to be within this interval, cropping outside data
		genome_track::encoding _encoding;
		chrom_map_t<strand_t, track_info_t> _tracks;
		chrom_names_t _chrom_names;
		std::string _refg_name;
		bool _finalized{};
	};
};

// Useful for converting C type to dtype enum at compile type
// Can't specialize withing genome_track scope, unfortunately.
template <typename T> struct dtype_traits { };
template <> struct dtype_traits<bool>    { static const genome_track::dtype_t dtype = genome_track::bool_;   static const bool is_signed = false; };
template <> struct dtype_traits<uint8_t> { static const genome_track::dtype_t dtype = genome_track::uint8;   static const bool is_signed = false; };
template <> struct dtype_traits<int8_t > { static const genome_track::dtype_t dtype = genome_track::int8;    static const bool is_signed = true;  };
template <> struct dtype_traits<half_t > { static const genome_track::dtype_t dtype = genome_track::float16; static const bool is_signed = true;  };
template <> struct dtype_traits<float  > { static const genome_track::dtype_t dtype = genome_track::float32; static const bool is_signed = true;  };

template <typename T> INLINE const T* genome_track::float_dict::get() const { return nullptr; }
template <> INLINE const half_t* genome_track::float_dict::get<half_t>() const { return h; }
template <> INLINE const float*  genome_track::float_dict::get<float>()  const { return f; }
INLINE float  genome_track::float_dict::min() const { return f[0]; }
INLINE float  genome_track::float_dict::max()  const { return f[nval-1]; }

template <typename T> void clamp_track_data(T* data, int size, int dim, float min, float max);
template <> void clamp_track_data(float* data, int size, int dim, float min, float max);

float clamp_single(float val, float min, float max);

template <typename T>
inline void clamp_track_data(T* data, int size, int dim, float min, float max)
{
    for (size_t i = 0; i < (size_t)size*dim; ++i) {
        float x = as_float(data[i]);
        if (x < min) x = min;
        if (x > max) x = max;
        data[i] = clamp_for_half(x);
    }
}

template <> void clamp_track_data(float* data, int size, int dim, float min, float max);

END_NAMESPACE_GK

#endif // __GENOME_KIT_GENOME_TRACK__
