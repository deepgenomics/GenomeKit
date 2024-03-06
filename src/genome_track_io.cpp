/*
Copyright (C) 2016-2023 Deep Genomics Inc. All Rights Reserved.
*/
#include "genome_track.h"

#include "genome.h"
#include "strutil.h"
#include "util.h"
#include <algorithm>
#include <bit>
#include <cmath>
#include <cstring>
#include <filesystem>
#include <fstream>
#include <map>
#include <utility>
#include <variant>

#define GK_CHECK_NODATA(funcname) \
	GK_CHECK(index_size() == 0, runtime, "Cannot call " #funcname " after data has been added");

BEGIN_NAMESPACE_GK

using std::min;
using std::max;
using detail::any_t;

const unsigned short c_gtrack_sig = 0x70ac;
const unsigned short c_gtrack_ver = 0x0005;
// versions:
//   0001: initial format
//   0002: add 32-bit index mode; add inplace-data mode; move track fields into index header
//   0003: support track resolutions
//   0004: interleave index with data
//   0005: Restructured to support all chromosomes and more species
//   <--- INSERT VERSION CHANGE SUMMARIES HERE

extern int dtype_size[];
extern const char* dtype_as_cstr[];

/************************************************************************

GTRACK File Overview
--------------------

A GTRACK file contains a set of "data blocks" for each (chrom, strand) pair.
Each data block stores the values for a particular genomic interval.
Any genomic positions not represented by a data block are default_value.

For example, consider the numeric track below::

    0   1   2   3   4   5   6   7   8   9  10  11
    0   0   0  .1  .5  .3   0   0  .7  .2   0   0  ...

Using a default_value of 0, that track can be represented by two data
blocks::

    chr1:+:3-6  : [.1, .5, .3]
    chr1:+:8-10 : [.7, .2]

or, equivalently, by one data block::

    chr1:+:3-10  : [.1, .5, .3, 0, 0, .7, .2]

Data block intervals are often the same intervals that were specified by the
user when the file was created, though if genome_track::builder::set_sparsity
was configured, the final intervals may be different.

GTRACK File Format
------------------

Each GTRACK file has the following sections::

    +================================+       Header contains encoding type,
    |                                |_      dimensionality, default value,
    | header                         | \     min/max range, strandedness, etc.
    |                                |  \
    +================================+   \
    | chr1:- [data, index]           | <--|  For each (chrom, strand) pair,
    | chr1:+ [data, index]           | <--|  header contains:
    | chr2:- [data, index]           | <--|    1. file offset to its data section
    | chr2:+ [data, index]           | <--|    2. file offset to its index section
    | ...                            | <--|
    |                                | <--|  Each index contains offsets
    | chrM:- [data, index]           | <--|  into its corresponding data section,
    | chrM:+ [data, index]           | <--|  to locate individual data blocks.
    +================================+

To read about how the decoding algorithm works, see the inline comments
in genome_track::operator().

If a track is not stranded, then the header's negative-strand file offsets will
be copies of the positive-strand offsets, so that queries on both strands share
the same data.

NOTE ON ENDIAN-NESS:
    The current file format does not store the endian-ness of the architecture
    which wrote the file, so an error will be raised (signature mismatch)
    if the file was created on a different endian-ness.

HEADER
^^^^^^

The header is simply an instance of the genome_track::header_t structure.
Please see the inline documentation in genome_track.h.

DATA
^^^^

Each data section is specific to a (chrom, strand), and comprises
data blocks arranged back-to-back. Recall that a "data block" is
a chunk of contiguous data associated with a genomic interval.

For example, consider two intervals with integer data defined as::

    chr1:+:100-103 : [3,4,5]
    chr1:+:200-202 : [8,9]

Encoding these two data blocks as u8 (8-bit unsigned) results
a 5 byte chr1:+ data block sub-section::

    byte  value   bits
    [0]   3       00000011
    [1]   4       00000100
    [2]   5       00000101
    [3]   8       00001000
    [4]   9       00001001
                  ^      ^
                  MSB    LSB

Encoding these short data blocks as u4 (4-bit unsigned)
results in 8 bytes stored. This is because fractional
encodings (1/2/4-bit) are stored into 32-bit dwords::

    dword   values   bits
    [0]     3,4,5    00000000 00000000 00000101 01000011
    [1]     8,9      00000000 00000000 00000000 10011000
                     ^                                 ^
                     MSB                               LSB

In other words, the first position in a data block is always
encoded into the LSB of a new dword, not packed into the residual
bits of the previous interval's dword.

Note that on little-endian, the physical byte order within a dword
is different, but bit significance is still logically consistent
with the above.

Encoding the same data as f16 (16-bit half-precision float)
results in 10 bytes stored::

    half_t  value   bits
    [0]     3.0     01000010 00000000
    [1]     4.0     01000100 00000000
    [2]     5.0     01000101 00000000
    [3]     8.0     01001000 00000000
    [4]     9.0     01001000 10000000
                    ^               ^
                    MSB             LSB

NOTE ON RESOLUTION:
    If a track has resolution R, then a data block only encodes `dim`
    values for every R genomic positions it spans.
    For example, a track with R=5 might have a data block with only
	two values [0.57, 0.23] and it would decode as follows:

       0    1    2    3    4    5    6    7    8    9      <-- genomic position
       0.57 0.57 0.57 0.57 0.57 0.23 0.23 0.23 0.23 0.23   <-- decoded value

NOTE ON INPLACE DATA:
    Small data blocks may be packed into the index section itself, rather
    than the data section. See the INDEX documentation for details.

NOTE ON MASKS:
    For the special 'm0' encoding type, the data section is empty, since
    a 0-bit mask is entirely specified by the intervals in the index.

INDEX
^^^^^

Each index section is specific to a (chrom, strand), and immediately follows
the corresponding data section (subject to 4-byte alignment).

Each index comprises header fields and four parallel arrays::

    // INDEX HEADER
    int32_t  num_blocks;                 // Number of data blocks.
    int32_t  num_jumps;                  // Number of data block indices.
    uint8_t  is_32bit;                   // 32 or 64-bit data_offsets.
    uint8_t  reserved[3];

    // INDEX ARRAYS
    int32_t  jumps[num_jumps];           // Sorted, increasing order.
    pos_t    ends[num_blocks];           // Sorted, increasing order.
    pos_t    starts[num_blocks];         // Sorted, increasing order.
    uintNN_t data_offsets[num_blocks];   // Offset within the data sub-section
                                         // of this particular (chrom, strand).
                                         // (NN = 32 or 64 depending on is_32bit)


The number of genomic positions spanned by data block i is
`(ends[i] - starts[i])*resolution`, where the track may have
resolution >1bp.
For each position, `dim` consecutive values will be encoded
in the data block.
For the `jumps` array, jumps[j] is the index of the first element in
the `ends` array such that `ends[jumps[j]] >= j*jump_size`.

The decoding algorithm uses the index arrays as follows:

1. Use the `jumps` array to find a narrow range of block indices
   in O(1) time, merely to speed up the binary search in step 2.

2. Find the first data block with `end` beyond the query interval's
   start position, by performing binary search over `ends`.

3. Start decoding data blocks from left-to-right, using the
   `starts` and `ends` arrays to determine output positions,
   and the `data_offsets` array to find the encoded data.

4. If the track resolution is >1, then expand the coarse-grained
   data into full-resolution track.

The data for block i of (chrom, strand) is stored at file offset::

    header.tracks[chrom][strand].data_offset + data_offsets[i]

In other words, `data_offsets[i]` is relative to the start of
the data section for this particular (chrom, strand).

INPLACE DATA:
    If is_32bit is true, and data block i encodes to 31 or fewer bits,
    then the encoded data block is packed into the lower 31 bits of
    data_offset[i]. In this case, the upper-most bit is set to 1, indicating
    "data_offsets[i] stores data inplace, rather than storing an actual
    offset into the data section."
	For tracks comprising many short intervals, this scheme results
    in much smaller files.

NOTE ON MASKS:
    For the special 'm0' encoding type, the data_offsets array is empty,
    since a 0-bit mask is entirely specified by the intervals in the index.

************************************************************************/

void genome_track::open()
{
	GK_CHECK(!_fmap.is_open(), runtime, "genome_track::open() already opened");
	GK_CHECK(!_sourcefile.empty(), value, "genome_track::open() failed; no source file specified");

	// Memory map the source file
	_fmap.open(_sourcefile);

	// Copy a few key GTRACK header fields so there's less indirection when accessing data members from decoder
	const auto h = _fmap.as_ptr<header_t>(0);
	GK_CHECK(h->sig == c_gtrack_sig, file, "Expected valid GTRACK file signature {:x} but found {:x}.", c_gtrack_sig, h->sig);
	GK_CHECK(h->ver == c_gtrack_ver, file, "Expected GTRACK file version {:x} but found {:x}.", c_gtrack_ver, h->ver);
	_dim       = h->dim;
	_res       = h->res;
	_refg      = h->refg;
	_stranded  = h->stranded;

	// Initialize the encoding and then copy the float dict into the encoding
	_encoding.init(h->etype, _dim, _res, h->default_value);
	memcpy(&_encoding.dict, &h->dict, sizeof(float_dict));

	// Initialize the tracks pointers from the info. Copy the num_blocks to
	// simplify code and memory access pattern.
	_fmap.set_seek(sizeof(*h));
	std::generate_n(std::inserter(_tracks, end(_tracks)), h->num_tracks, [&]() {
		const auto info = _fmap.read<header_track_info_t>();
		return decltype(_tracks)::value_type{
			{						  info.chrom,									info.strand},
			{_fmap.as_ptr<void>(info.data_offset), _fmap.as_ptr<track_index_t>(info.index_offset)}
        };
	});

	// That's it! -------------
	// Try to avoid touching any data in the rest of the file until user requests data.
	// Do NOT try to convert index data_offsets into direct pointers to the data, since that
	// would make opening a GTRACK file much slower and page in lots of unnecessary data.
}

int genome_track::gtrack_version() { return c_gtrack_ver; }

/////////////////////////////////////////////////////////////////////

void genome_track::encoding::init_dict()
{
	if (bits_per_encoded_datum > 8)
		return;
	int size = 1 << bits_per_encoded_datum;
	auto range = (float)(size - 1);
	float f[256];
	for (int i = 0; i < size; ++i)
		f[i] = (float)i/range;
	dict.init(f, size);
}

void genome_track::encoding::init(etype_t etype, int dim, int res, any_t default_value)
{
	this->default_value = default_value;
	this->etype = etype;

	// The code below specializes the encoders[] and decoders[] arrays
	// to match the etype and the dim given.
	//
	// Looks like a pretty convoluted way to intialize what otherwise
	// seems to be a pretty standard class hierarchy, but there is some
	// small method in this madness, related to getting tighter control
	// over how polymorphism happens, and simplifying serialization
	// somewhat.
	switch (etype) {
	case m0:    ((m0_encoding*)this)->init(dim, res); break;
	case u1:    ((u1_encoding*)this)->init(dim, res); break;
	case u2:    ((u2_encoding*)this)->init(dim, res); break;
	case u3:    ((u3_encoding*)this)->init(dim, res); break;
	case u4:    ((u4_encoding*)this)->init(dim, res); break;
	case u5:    ((u5_encoding*)this)->init(dim, res); break;
	case u6:    ((u6_encoding*)this)->init(dim, res); break;
	case u8:    ((u8_encoding*)this)->init(dim, res); break;
	case i8:    ((i8_encoding*)this)->init(dim, res); break;
	case f2:    ((f2_encoding*)this)->init(dim, res); break;
	case f3:    ((f3_encoding*)this)->init(dim, res); break;
	case f4:    ((f4_encoding*)this)->init(dim, res); break;
	case f5:    ((f5_encoding*)this)->init(dim, res); break;
	case f6:    ((f6_encoding*)this)->init(dim, res); break;
	case f8:    ((f8_encoding*)this)->init(dim, res); break;
	case f16:   ((f16_encoding*)this)->init(dim, res); break;
	case f32:   ((f32_encoding*)this)->init(dim, res); break;
	default: GK_UNREACHABLE();
	}
}

genome_track::builder::builder(string outfile, etype_t etype, strandedness_t strandedness, const genome_t& genome, int dim, int res)
	: _outfile(std::move(outfile))
	, _sparsity_min_delta(0, 0, 0, 0)
	, _sparsity_min_run(0)  // 0 means no sparsity
	, _strandedness(strandedness)
	, _data_transform(false)
	, _data_clamp(false)
	, _data_scale(1.0f)
	, _data_bias(0.0f)
	, _data_min(nanf())
	, _data_max(nanf())
	, _data_size(0)
	, _index_size(0)
	, _chrom_names{genome.chrom_names()}
	, _refg_name{genome.refg_name()}
{
	GK_CHECK(dim >= 1, value, "Must have dim >= 1");

	// Initialize header. Will be written to disk at the end of finalize().
	_h.sig = c_gtrack_sig;
	_h.ver = c_gtrack_ver;
	_h.etype = etype;
	_h.dim = dim;
	_h.res = res;
	_h.refg = genome.refg();
	_h.stranded = (strandedness == strandedness_t::strand_unaware) || (strandedness == strandedness_t::strand_aware);
	_h.default_value = any_t();
	_clamp_float = [](float x) -> float { return x; };

	_encoding.init(_h.etype, _h.dim, _h.res, _h.default_value);

	_h.max_value = _encoding.range_min;  // deliberately set max_value to range_min
	_h.min_value = _encoding.range_max;  // so that update_min_max works properly
}

void genome_track::builder::set_sparsity(int min_run, float min_delta)
{
	GK_CHECK_NODATA("set_sparsity");
	GK_CHECK(min_run > 0, value, "min_run must be > 0");
	GK_CHECK(min_delta >= 0, value, "min_delta must be >= 0");
	GK_CHECK(_h.etype != m0, value, "Cannot set sparsity on m0 track");
	_sparsity_min_run = min_run;
	_sparsity_min_delta = any_t(min_delta);
}

void genome_track::builder::set_transform(float a, float b, float c, float d)
{
	GK_CHECK_NODATA("set_transform");
	GK_ASSERT(a != b);
	GK_ASSERT(c != d);
	_data_transform = true;
	_data_scale = (d - c) / (b - a);
	_data_bias  = (c*b - a*d) / (b - a);
	_data_min   = c;
	_data_max   = d;
}

template <typename T>
void reverse_track_data(T* dst, const T* src, int size, int dim)
{
	for (int i = 0; i < size; ++i)
		for (int j = 0; j < dim; ++j)
			dst[(size_t)(size-1-i)*dim + j] = src[(size_t)i*dim + j];
}

template <typename T>
void transform_track_data(T* data, int size, int dim, float scale, float bias)
{
	for (size_t i = 0; i < (size_t)size*dim; ++i)
		data[i] = as_float(data[i])*scale + bias;
}

template <> void transform_track_data(bool*    data, int size, int dim, float scale, float bias) { GK_THROW(type, "Cannot transform data of dtype bool"); }
template <> void transform_track_data(uint8_t* data, int size, int dim, float scale, float bias) { GK_THROW(type, "Cannot transform data of dtype uint8"); }
template <> void transform_track_data(int8_t*  data, int size, int dim, float scale, float bias) { GK_THROW(type, "Cannot transform data of dtype int8"); }


inline float clamp_single(float val, float min, float max) {
    if (val < min) return min;
    if (val > max) return max;
    return val;
}

template <> void clamp_track_data(float* data, int size, int dim, float min, float max) {
	for (size_t i = 0; i < (size_t)size*dim; ++i) {
        data[i] = clamp_single(data[i], min, max);
	}
}

template <> void clamp_track_data(bool*    data, int size, int dim, float min, float max) { GK_THROW(type, "Cannot clamp data of dtype bool"); }
template <> void clamp_track_data(uint8_t* data, int size, int dim, float min, float max) { GK_THROW(type, "Cannot clamp data of dtype uint8"); }
template <> void clamp_track_data(int8_t*  data, int size, int dim, float min, float max) { GK_THROW(type, "Cannot clamp data of dtype int8"); }

void update_min_max_case(bool x, any_t& min_value, any_t& max_value)
{
	min_value = any_t(0, 0, 0, 0);
	max_value = any_t(1, 1, 1, 1);
}

void update_min_max_case(uint8_t x, any_t& min_value, any_t& max_value)
{	// Don't try to detect uint8_t->int8_t overflow here since we don't know which decoders user will actually use
	if (x < min_value.u) min_value = any_t(x, (int8_t)min(x, (uint8_t)127), x, x);
	if (x > max_value.u) max_value = any_t(x, (int8_t)min(x, (uint8_t)127), x, x);
	// TODO: when user requests int8_t decoder, for example, that's the right moment to do a
	//       GK_CHECK that max.as_uint8 <= 127, thereby indicating it's safe to decode as i8
	//       without overflowing to negative values. For now, we'll just clamp the values and hope
	//       the user is aware of whether a track would over/underflow when decoding to a specific dtype.
	//       (this is somewhat pedantic, likely there will be one obvious and consistent dtype
	//        for each track, so "ideal" behaviour in these corner cases is not urgent)
}

void update_min_max_case(int8_t x, any_t& min_value, any_t& max_value)
{	// Don't try to detect int8_t->uint8_t underflow here since we don't know which decoders user will actually use
	if (x < min_value.i) min_value = any_t((uint8_t)max(x, (int8_t)0), x, x, x);
	if (x > max_value.i) max_value = any_t((uint8_t)max(x, (int8_t)0), x, x, x);
}

void update_min_max_case(half_t x, any_t& min_value, any_t& max_value)
{	// Ignore uint8 and i8 since no decoders for it
	if (is_nan(min_value.h) || x < min_value.h) min_value = any_t(0, 0, x, as_float(x));
	if (is_nan(min_value.h) || x > max_value.h) max_value = any_t(0, 0, x, as_float(x));
}

void update_min_max_case(float x, any_t& min_value, any_t& max_value)
{	// Ignore uint8 and i8 since no decoders for it
	if (is_nan(min_value.f) || x < min_value.f) min_value = any_t(0, 0, 0, x);
	if (is_nan(min_value.f) || x > max_value.f) max_value = any_t(0, 0, 0, x);
}

template <typename T>
void update_min_max(const T* data, int size, int dim, any_t& min_value, any_t& max_value)
{
	size_t n = (size_t)size*dim;
	for (size_t i = 0; i < n; ++i)
		update_min_max_case(data[i], min_value, max_value);
}

void genome_track::builder::set_default_value(bool    value) { GK_CHECK_NODATA("set_default_value"); _h.default_value = _encoding.default_value = any_t(value); }
void genome_track::builder::set_default_value(int     value) { GK_CHECK_NODATA("set_default_value"); _h.default_value = _encoding.default_value = any_t(float(value)); }
void genome_track::builder::set_default_value(half_t  value) { GK_CHECK_NODATA("set_default_value"); _h.default_value = _encoding.default_value = any_t(value); }
void genome_track::builder::set_default_value(float   value) {
	GK_CHECK_NODATA("set_default_value");
	if (_h.etype == f32)
		_h.default_value = _encoding.default_value = any_t(0, 0, 0, value);
	else
		_h.default_value = _encoding.default_value = any_t(value);
}

void genome_track::builder::set_clamping()
{
	GK_CHECK_NODATA("set_clamping");
	_data_clamp = true;
	_data_min = _encoding.dict.min();
	_data_max = _encoding.dict.max();
	_clamp_float = [&](float x) -> float {
		return as_float(clamp_single(x, _encoding.dict.min(), _encoding.dict.max()));
	};
}

template <typename T>
void genome_track::builder::set_dict_impl(const T* dict)
{
	switch (_h.etype) {
	case f2: _encoding.dict.init(dict, 4);   break;
	case f3: _encoding.dict.init(dict, 8);   break;
	case f4: _encoding.dict.init(dict, 16);  break;
	case f5: _encoding.dict.init(dict, 32);  break;
	case f6: _encoding.dict.init(dict, 64);  break;
	case f8: _encoding.dict.init(dict, 256); break;
	default: GK_THROW(value, "Cannot call set_dict on etype '{}'", etype_as_cstr[_h.etype]);
	}
	_data_min = _encoding.dict.f[0];
	_data_max = _encoding.dict.f[_encoding.dict.nval-1];
}

void genome_track::builder::set_dict(const half_t* dict) { set_dict_impl(dict); }
void genome_track::builder::set_dict(const float*  dict) { set_dict_impl(dict); }

void genome_track::builder::set_restriction(const interval_t& r)
{
	GK_CHECK_NODATA("set_restriction");
	// Align the restriction to the native resolution of the track,
	// rounding start down, and end up, to the nearest aligned position.
	_restriction = interval_t::from_dna0(r.chrom, rnddn(r.start(), _h.res),
	                                              rndup(r.end(),   _h.res), pos_strand, r.refg);
}

void genome_track::builder::set_data(const interval_t& interval, const bool*    data) { set_data_impl(interval, data); }
void genome_track::builder::set_data(const interval_t& interval, const uint8_t* data) { set_data_impl(interval, data); }
void genome_track::builder::set_data(const interval_t& interval, const int8_t*  data) { set_data_impl(interval, data); }
void genome_track::builder::set_data(const interval_t& interval, const half_t*  data) { set_data_impl(interval, data); }
void genome_track::builder::set_data(const interval_t& interval, const float*   data) { set_data_impl(interval, data); }
void genome_track::builder::set_data(const interval_t& interval, const void*    data, dtype_t dtype)
{
	switch (dtype) {
	case bool_:   set_data(interval, rcast<const bool*   >(data)); break;
	case uint8:   set_data(interval, rcast<const uint8_t*>(data)); break;
	case int8:    set_data(interval, rcast<const int8_t* >(data)); break;
	case float16: set_data(interval, rcast<const half_t* >(data)); break;
	case float32: set_data(interval, rcast<const float*  >(data)); break;
	default: GK_UNREACHABLE();
	}
}

INLINE bool genome_track::builder::is_default(bool    x) const { return x == (_h.default_value.u != 0); }
INLINE bool genome_track::builder::is_default(uint8_t x) const { return std::abs((int)x - (int)_h.default_value.u) <= (int)_sparsity_min_delta.u; }
INLINE bool genome_track::builder::is_default(int8_t  x) const { return std::abs((int)x - (int)_h.default_value.i)  <= (int)_sparsity_min_delta.i; }
INLINE bool genome_track::builder::is_default(half_t  x) const { return (is_nan(x) && is_nan(_h.default_value.f)) ||
                                                                        (fabs(_h.default_value.f - as_float(x)) <= _sparsity_min_delta.f); }
INLINE bool genome_track::builder::is_default(float   x) const { return (is_nan(x) && is_nan(_h.default_value.f)) ||
                                                                        (fabs(_h.default_value.f - x) <= _sparsity_min_delta.f); }

genome_track::builder::encoded_sizes_t genome_track::builder::encoded_sizes(int length) const
{
	return {.bytes = _encoding.num_encoded_bytes(length, _h.dim),
			.bits  = _encoding.num_required_bits(length, _h.dim)};
}

void genome_track::builder::add_track_entry(track_info_t::adder& adder, const span_t& span, std::unique_ptr<uint8_t[]>&& data)
{
	auto had_data = scast<bool>(data);
	adder.add(span, std::move(data), encoded_sizes(span.size()));
	_index_size++;
	if (had_data)
		_data_size += span.size();
}

template <typename T>
void genome_track::builder::set_data_impl(const interval_t& interval, const T* data)
{
	// Sanity checks
	GK_CHECK(_h.refg == interval.refg, value, "Mismatched reference genome");
	GK_CHECK(_h.stranded || interval.is_pos_strand(), value, "Cannot specify interval on negative strand for an unstranded track");
	GK_CHECK(!finalized(), runtime, "Cannot set data after calling finalize");

	// Figure out if we're going to need to reverse or sparsify the data before final encoding.
	bool sparsify = _sparsity_min_run && _sparsity_min_delta.f > 0;
	bool reverse  = interval.is_neg_strand() && (_strandedness == strandedness_t::strand_aware);

	// Get the interval as a span (a, b) as coarse positions.
	span_t span(interval.start(), interval.end());
	if (_h.res > 1) {
		GK_CHECK((span.a % _h.res == 0) && (span.b % _h.res == 0), value,
		         "Interval {}-{} a and b must be aligned to resolution={}", span.a, span.b, _h.res);
		span.a /= _h.res;
		span.b /= _h.res;
		// ----- From this point forward, (b-a)*dim matches number of values in 'data' -----
	}

	// Trim the span to the restriction, if applicable.
	if (_restriction) {
		//
		// Example:
		//    resolution: 2 bp
		//    interval: (0, 10)
		//    data: ABCDE
		//
		// Without restriction, this track looks as follows
		//
		//    0123456789   <-- genome positions
		//    AABBCCDDEE   <-- data
		//
		// With restriction of (2, 8) this track looks like
		//
		//    0123456789   <-- genome positions
		//      BBCCDD     <-- restricted data
		//
		GK_CHECK(_restriction->refg == _h.refg, value, "Reference genome of restriction interval does not match track.");

		// Convert restriction interval to coarse coordinates.
		span_t rspan(_restriction->start() / _h.res, _restriction->end() / _h.res);

		// If totally outside restriction zone, do nothing.
		if (interval.chrom != _restriction->chrom || span.b <= rspan.a || span.a >= rspan.b)
			return;

		// Trim left side of span.
		if (span.a < rspan.a) {
			if (!reverse)
				data += (rspan.a - span.a) * _h.dim;
			span.a = rspan.a;
		}

		// If right side of span needs to be trimmed, do so.
		if (span.b > rspan.b) {
			if (reverse)
				data += (span.b - rspan.b) * _h.dim;
			span.b = rspan.b;
		}
	}

	pos_t size = span.size();
	if (size <= 0)
		return;

	// Find the location to insert the new interval, and make sure
	// it doesn't overlap an existing interval.
	track_info_t& ti = _tracks[{interval.chrom, interval.strand}];
	GK_CHECK(!ti.flushed(), value, "Cannot set data on a chromosome ({}) that was already flushed to disk",
			 chrom_names().chrom_as_sv(interval.chrom));
	track_info_t::adder adder{&ti};
	try {
		adder.validate(span);
	}
	GK_RETHROW("for '{}'", interval);

	// Special case when there's no actual data to encode.
	if (_encoding.etype == m0) {
		GK_CHECK(!data, value, "Data for etype 'm0' must be NULL");
		add_track_entry(adder, span, nullptr);
		return;
	}
	GK_CHECK(data, value, "Data required for etype '{}'", etype_as_cstr[_encoding.etype]);

	// Allocate a new array to contain the encoded bytes
	// Encode the data into a new byte array, according to etype (via _encoding) and dtype.
	dtype_t dtype = dtype_traits<T>::dtype;
	int dim = _h.dim;
	auto dst = std::make_unique<uint8_t[]>(_encoding.num_encoded_bytes(size, dim));
	encoding::encode_fn encode = _encoding.encoders[dtype];
	GK_CHECK(encode, type, "Cannot encode as {} using decoded type {}", etype_as_cstr[_encoding.etype], dtype_as_cstr[dtype]);

	// The source pointer we'll 'encode' from, which may be original data or
	// may be redirected to transformed data.
	const T* src = data;

	// If the data needs to be transformed somehow before we encode it,
	// then create an alt_data buffer to keep a modified copy of the data.
	std::unique_ptr<T[]> alt_data;
	if (reverse || sparsify || _data_transform || _data_clamp) {
		// Negative strand and the data received in sense-strand order.
		// Reverse the order so that it's written in reference-strand order.
		// (The genome_track decoders always assumes blocks of
		//  contiguous data are stored in reference-strand order.)
		alt_data = std::make_unique<T[]>(size*dim);
		src = &alt_data[0];
		if (reverse)
			reverse_track_data(&alt_data[0], data, size, dim);
		else
			memcpy(&alt_data[0], data, size*dim*dtype_size[dtype]);

		// Transform the data if needed
		if (_data_transform)
			transform_track_data(&alt_data[0], size, dim, _data_scale, _data_bias);

		// Clamp the data if needed
		if (_data_clamp)
			clamp_track_data(&alt_data[0], size, dim, _data_min, _data_max);

		// Convert default-equivalent data to default_value so that they get encoded as such.
		// This is important for allowing _sparsity_min_delta to attract a larger range of values
		// than would normally be achieved by merely rounding to nearest encodable value.
		if (sparsify)
			for (int i = 0; i < size*dim; ++i)
				if (is_default(alt_data[i]))
					alt_data[i] = _h.default_value.as<T>();
	}

	// Finally, encode the data itself.
	encode(&dst[0], src, _encoding.dict, size, dim);

	// If encoding succeeded, update the min_value and max_value for the gtrack header
	update_min_max(src, size, dim, _h.min_value, _h.max_value);

	// If sparsify mode enabled, decode the data block and re-encode only the chunks that are not default_value
	if (_sparsity_min_run) {
		encoding::decode_fn decode = _encoding.decoders[dtype][as_ordinal(pos_strand)];  // decode forward, regardless of interval strand
		GK_ASSERT(decode);  // should always be a decoder for dtype if there was an encoder
		if (!alt_data)
			alt_data = std::make_unique<T[]>(size*dim);
		decode(&alt_data[0], &dst[0], _encoding.dict.get<T>(), size, dim, 0, 0);  // Fill `alt_data` by decoding `dst`

		// Now that we've decoded the encoded data, thereby quantizing the data as it would
		// appear in the final track, we're ready to scan the decoded copy for blocks of
		// data to insert. We skip runs of default_value so long as they're longer than
		// _sparsity_min_run. The is_default method is sensitive to _sparsity_min_delta,
		// which allows the default value to swallow nearby values, as a more aggressive
		// mechanism than mere quantization-to-default_value.

		int start = 0;  // index of next alt_data position that has a non-default value
		for (;;) {
			// Find the next position where at least ONE element IS NOT default_value
			while (start < size) {
				int j = 0;
				while (j < dim && is_default(alt_data[start*dim+j]))
					++j;
				if (j < dim)
					break;
				start++;
			}
			if (start >= size)
				break;

			// Find the next position where ALL elements are default_value.
			int end = start+1;
			while (end < size) {
				int j = 0;
				while (j < dim && is_default(alt_data[end*dim+j]))
					++j;

				// At least one non-default value? If so, keep scanning
				if (j < dim) {
					end++;
					continue;
				}

				// Now we know j==dim and so 'end' indexes a position with all-default values.
				// Scan ahead and look a spot with at least ONE element is not default_value
				int next_start = end+1;
				while (next_start < size) {
					int k = 0;
					while (k < dim && is_default(alt_data[next_start*dim+k]))
						++k;
					if (k < dim)
						break;
					next_start++;
				}
				if (next_start - end >= _sparsity_min_run)
					break;  // Encode (start, end) since the gap between end and next_start is large enough

				// Otherwise, continue looping beyond this gap.
				end = next_start;
			}

			// Add the sub-interval data after re-encoding it.
			span_t sub_span(span.a + start, span.a + end);
			auto sub_dst = std::make_unique<uint8_t[]>(_encoding.num_encoded_bytes(sub_span.size(), dim));

			// Encode the sub-interval into the smaller dst buffer
			encode(&sub_dst[0], &alt_data[start*dim], _encoding.dict, sub_span.size(), dim);

			// Insert the sub-interval into the data block map, and leave the memory be
			add_track_entry(adder, sub_span, std::move(sub_dst));

			// Move on to the next position in our search for a .
			start = end+1;
			if (start >= size)
				break;
		}

	} else {
		// Insert the new data block into the tracks data structure
		add_track_entry(adder, span, std::move(dst));
	}
}

void genome_track::builder::flush()
{
	// sort the chromosomes so files are stable between hash_map implementations
	std::vector<typename decltype(_tracks)::key_type> keys;
	keys.reserve(std::size(_tracks));
	std::transform(std::begin(_tracks), std::end(_tracks), std::back_inserter(keys), [](auto& kv) { return kv.first; });
	std::sort(std::begin(keys), std::end(keys));
	for (auto key : keys) flush(key.chrom, key.other);
}

namespace {
std::string partial_filename(std::string base, chrom_t chrom, strand_t strand)
{
	return base + "." + std::to_string(as_ordinal(chrom)) + "." + std::to_string(as_ordinal(strand));
}
}  // namespace

void genome_track::builder::flush(std::optional<chrom_t> chrom, strand_t strand)
{
	if (!chrom)
		return;

	GK_CHECK(!finalized(), runtime, "Cannot flush after calling finalize");
	auto& track_info = _tracks[{*chrom, strand}];
	if (track_info.flushed())
		return;  // avoid truncating any existing flushed files

	// flush track to disk (merge later) since don't know the total amount of chromosomes to write for the headesr
	binary_file track_segment_file(partial_filename(_outfile, *chrom, strand), "w");
	track_info.flush(track_segment_file, [&](int length) { return encoded_sizes(length); });
}

void genome_track::builder::track_info_t::adder::validate(span_t span)
{
	auto hint = _ti->_blocks.lower_bound(span);  // it = first interval equivalent-to or later-than 'span'
	if (hint != std::begin(_ti->_blocks))
		GK_CHECK(std::prev(hint)->first.b <= span.a, value, "Overlapping blocks are not allowed");
	if (hint != std::end(_ti->_blocks))
		GK_CHECK(hint->first.a >= span.b, value, "Overlapping blocks are not allowed");
	_inserter = std::inserter(_ti->_blocks, hint);
}
void genome_track::builder::track_info_t::adder::add(span_t span, std::unique_ptr<uint8_t[]>&& data,
													 genome_track::builder::encoded_sizes_t encoded_sizes)
{
	const auto outplace_mask = [&]<class T>() { return _ti->can_inplace<T>(encoded_sizes) ? 0 : 1; };
	const auto outplace64    = outplace_mask.template operator()<uint64_t>() * encoded_sizes.bytes;
	GK_CHECK(_ti->_max_offset64 + outplace64 < high_bit<uint64_t>(), runtime,
			 "Track size exceeds maximum supported (2^63-1): {}",
			 _ti->_max_offset64 + outplace64);

	_inserter = {span, std::move(data)};
	// mark outplace data sizes to determine if we exceeded the amount of inplace bits for real data offsets
	_ti->_max_offset32 += outplace_mask.template operator()<uint32_t>() * encoded_sizes.bytes;
	_ti->_max_offset64 += outplace64;
	// mark max extent for jump count
	_ti->_max_end = max(_ti->_max_end, span.b);
}
bool genome_track::builder::track_info_t::allow_32bit() const
{
	return _max_offset32 < high_bit<uint32_t>() &&
		   // if 32bit offsets are possible, choose the smaller data_size after factoring in inplace optimizations
		   _max_offset32 + num_blocks() * sizeof(uint32_t) <= _max_offset64 + num_blocks() * sizeof(uint64_t);
}
uint64_t genome_track::builder::track_info_t::data_size() const
{
	auto bytes = allow_32bit() ? _max_offset32 : _max_offset64;
	return bytes + aligned_distance<uint32_t>(bytes);  // pad to 32bit aligned (should this be 64?)
}
uint64_t genome_track::builder::track_info_t::index_size(bool encoded_data) const
{
	return sizeof(track_index_t) + num_jumps() * sizeof(int32_t) +
		   num_blocks() * (sizeof(pos_t) + sizeof(pos_t) +
						   (encoded_data ? (allow_32bit() ? sizeof(uint32_t) : sizeof(uint64_t)) : 0));
}
bool genome_track::builder::track_info_t::flushed() const { return _flushed; }
int  genome_track::builder::track_info_t::num_blocks() const
{
	return std::max(int_cast<int>(std::size(_blocks)), _num_blocks);
}
int genome_track::builder::track_info_t::num_jumps() const { return divup(_max_end + 1, (pos_t)jump_size); }

template <std::regular_invocable<int> EncodedSizesFn>
void genome_track::builder::track_info_t::flush(binary_file& out_file, EncodedSizesFn&& encoded_sizes)
{
	using std::visit;

	// If already flushed to disk, don't do it again
	if (flushed())
		return;

	// Flush the index header and the parallel arrays for this (chrom, strand).
	// The index data is written in the following format:
	//
	//   track_index_t header;               // header fields for the index
	//   int32_t  jumps[num_jumps];          // sorted, increasing order
	//   pos_t    ends[num_blocks];          // sorted, increasing order
	//   pos_t    starts[num_blocks];        // sorted, increasing order
	//   uintNN_t data_offsets[num_blocks];  // NN is 32 or 64 depending on index.is_32bit
	const track_index_t idx = {
		.num_blocks = num_blocks(),
		.num_jumps  = num_jumps(),
		.is_32bit   = allow_32bit(),
	};
	const bool encoded_data = encoded_sizes(1).bytes > 0;

	// Build dynamic arrays/write data first before writing the header
	vector<pos_t>                                    starts;
	vector<pos_t>                                    ends;
	std::variant<vector<uint32_t>, vector<uint64_t>> data_offsets;
	if (idx.is_32bit)
		data_offsets = vector<uint32_t>{};
	else
		data_offsets = vector<uint64_t>{};

	// Reserve exact capacity so that we don't end up using 2x the memory on vector padding.
	starts.reserve(idx.num_blocks);
	ends.reserve(idx.num_blocks);
	if (encoded_data)
		visit([&](auto& offsets) { offsets.reserve(idx.num_blocks); }, data_offsets);

	// Transfer the index information from the std::map into the parallel arrays/data file
	[[maybe_unused]] const auto old_data_file_size = out_file.tell();
	long long                   curr_offset        = 0;
	for (auto& [span, data] : _blocks) {
		starts.push_back(span.a);
		ends.push_back(span.b);
		if (!encoded_data)
			continue;

		visit(
			[&, span = span,
			 &data = data](  // bug in clang 12+'s C++20 support still doesn't allow capture on structured bindings
				auto& offsets) {
				using T          = typename std::decay_t<decltype(offsets)>::value_type;
				const auto sizes = encoded_sizes(span.size());

				if (can_inplace<T>(sizes)) {
					// If the bits in the offset encode the data inplace, and
					// no data was actually written to disk when this entry was
					// encountered in flush() so we don't count it towards the
					// offset.
					T x{};
					static_assert(std::endian::native == std::endian::little);
					memcpy(&x, data.get(), sizes.bytes);
					offsets.push_back(x | high_bit<T>());
					return;
				}
				// Otherwise, the data contributes to the data section, and so its
				// offset can be calculated from curr_offset.
				out_file.write(data.get(), sizes.bytes);
				offsets.push_back(int_cast<T>(curr_offset));
				curr_offset += sizes.bytes;
			},
			data_offsets);
	}
	curr_offset += out_file.write_until_align(sizeof(uint32_t));
	[[maybe_unused]] const auto old_index_file_size = out_file.tell();

	GK_DBASSERT(scast<decltype(data_size())>(curr_offset) == data_size());
	GK_DBASSERT(curr_offset == old_index_file_size - old_data_file_size);

	// Write the index arrays to the file, immediately after the header
	out_file.write(idx);
	if (idx.num_blocks > 0) {
		// Build the jump indices, i.e. jumps[i] is the index of the
		// first block that has end >= i*jump_size.
		vector<int> jumps(idx.num_jumps);
		int i = 0;
		for (int j = 0; j < idx.num_jumps; ++j) {
			while (ends[i] < j*jump_size)
				i++;
			jumps[j] = i;
		}

		out_file.write(&jumps[0], idx.num_jumps);
		out_file.write(&ends[0], idx.num_blocks);
		out_file.write(&starts[0], idx.num_blocks);

		visit(
			[&](auto& offsets) {
				if (!std::empty(offsets))
					out_file.write(&offsets[0], idx.num_blocks);
			},
			data_offsets);
	}
	GK_DBASSERT(index_size(encoded_data) == scast<decltype(index_size(true))>(out_file.tell() - old_index_file_size));

	// The std::map has been baked into the index structure, so free the memory.
	_num_blocks = int_cast<int>(std::size(_blocks));
	_blocks.clear();
	_flushed = true;
}

void genome_track::builder::wig_bedgraph_config(const char* pos_infile, const char* neg_infile, const char* filetype)
{
	GK_CHECK(!finalized(), runtime, "Cannot set data after calling finalize");
	GK_CHECK(index_size() == 0, runtime, "Cannot call set_data_from_{} after data has already been added to the track",
			 filetype);
	GK_CHECK(_strandedness != strandedness_t::strand_aware,
			 value, "Do not set 'strand_aware' when building from a {} file; such data always in reference-strand order", filetype);
	if (neg_infile)
		GK_CHECK(_strandedness != strandedness_t::single_stranded, value, "Cannot load data from two {} files in single-stranded mode", filetype);
	GK_ASSERT(pos_infile);
}

void genome_track::builder::set_data_from_wig(const string& infile)                               { set_data_from_wig(infile.c_str(), nullptr); }
void genome_track::builder::set_data_from_wig(const string& pos_infile, const string& neg_infile) { set_data_from_wig(pos_infile.c_str(), neg_infile.c_str()); }
void genome_track::builder::set_data_from_wig(const char* pos_infile, const char* neg_infile)
{
	wig_bedgraph_config(pos_infile, neg_infile, "wig");
	set_data_from_wig(pos_infile, pos_strand);
	if (neg_infile)
		set_data_from_wig(neg_infile, neg_strand);
}

void genome_track::builder::set_data_from_bedgraph(const string& infile)                               { set_data_from_bedgraph(infile.c_str(), nullptr); }
void genome_track::builder::set_data_from_bedgraph(const string& pos_infile, const string& neg_infile) { set_data_from_bedgraph(pos_infile.c_str(), neg_infile.c_str()); }
void genome_track::builder::set_data_from_bedgraph(const char* pos_infile, const char* neg_infile)
{
	wig_bedgraph_config(pos_infile, neg_infile, "bedgraph");
	GK_CHECK(_h.dim == 1, value, "Calling set_data_from_bedgraph requires dim=1");
	set_data_from_bedgraph(pos_infile, pos_strand);
	if (neg_infile)
		set_data_from_bedgraph(neg_infile, neg_strand);
}

template <typename T>
struct contig_parser { NOCOPY(contig_parser)
	contig_parser(contig_parser&&) = default;
	contig_parser& operator=(contig_parser&&) = default;

	// Parameters of the track
	genome_track::builder* bld_ptr;

	// Parameters of the interval
	std::optional<chrom_t>  chrom;
	strand_t strand;
	int      start;
	bool     is_fixed;
	bool     is_irregular_span;
	string   category;

	// Data for the above interval
	vector<T> data;

	// Temporary storage for parsing
	vector<string_view> attrs;
	vector<string_view> cols;

	contig_parser(genome_track::builder* bld, strand_t strand)
		: bld_ptr(bld)
		, strand(strand)
	{
		reset();
	}

	static float parse_float_value(string_view str)
	{
		return as_float(str);
	}

	// Parse a fixedStep line from the current line reader, making use of the temporary
	// arrays 'cols' ant 'attrs' for the parsing. Does not advance the line reader.
	template <class F>
	void parse_wig(line_reader& lr, F&& clamp_float)
	{
		if (lr.done())
			return;

		auto& bld = *bld_ptr;
		auto line = lr.line();
		is_fixed = line[0] == 'f';

		// parse "fixedStep chrom=? start=? step=? [span=?]"
		// or    "variableStep chrom=? [span=?]"
		split_view(line, ' ', attrs);
		chrom = bld.chrom_names().as_chrom(get_attr(attrs, "chrom"));
		start = as_pos(get_attr(attrs, "start", "1")) - 1;
		int span = as_int(get_attr(attrs, "span", "1"));

		// Sanity checks
		GK_CHECK(start >= 0, value, "Start {} must be non-negative", start);
		GK_CHECK(start % bld.res() == 0, value, "Start {} must be divisible by resolution={}", start, bld.res());
		if (is_fixed) {
			int step = as_int(get_attr(attrs, "step", "1"));
			GK_CHECK(span == step, value, "Expected 'span' to match 'step'");
		}

		// Take note of whether the span matches our track resolution;
		// only the last block on a chromosome is allowed to have span != resolution, which we'll
		// catch later on in try_extend.
		is_irregular_span = (span != bld.res());

		// Figure out how many columns we expect
		int num_cols_expected = is_fixed ? bld.dim() : bld.dim()+1;  // If variableStep, first column will be 'position'
		cols.resize(num_cols_expected+1);

		// Parse data rows until we hit the start of the next contig
		int num_datum = 0;
		while (!(++lr).done()) {
			// If we've encountered the next fixedStep or variableStep line, stop parsing this contig.
			line = lr.line();
			if (line[0] >= 'a' && line[0] <= 'z')
				break;

			// Split the data columns by tab separator. Own multi-dimensional extension to WIG format.
			int num_cols = split_view(lr.line(), '\t', &cols[0], num_cols_expected+1);
			GK_CHECK(num_cols == num_cols_expected, value, "Expected {} columns per line but found {}",
			                                               num_cols_expected, num_cols);

			int i = 0;
			// If variableStep, the first column is position of element to store.
			// Fill default_value until we reach that position in the values array.
			if (!is_fixed) {
				pos_t next_pos = as_pos(cols[i++])-1;
				GK_CHECK(next_pos >= 0, value, "Position {} must be non-negative", next_pos);
				if (start == 0)
					start = next_pos;
				fill_data_until(next_pos);
			}

			// Parse each float value in the current line
			for (; i < num_cols; ++i ) {
				float d = parse_float_value(cols[i]);
				T value(clamp_float(d));
				data.push_back(value);
			}

			num_datum++;
		}

		// If the fixed/variableStep data we just parsed has more than one data row, then
		// check that its span matches track resolution. (If there's only a single data row,
		// we may be at the end of a chromosome, for which the last data row is allowed to have
		// irregular span.)
		if (num_datum > 1)
			GK_CHECK(!is_irregular_span, value, "Expected span={} but found span={}", bld.res(), span);
	}

	template <class F>
	void parse_bedgraph(line_reader& lr, F&& clamp_float)
	{
		auto& bld = *bld_ptr;

		// Skip comment lines
		while (!lr.done() && lr.line()[0] == '#')
			++lr;
		if (lr.done())
			return;

		// parse "chrom start end data"
		cols.resize(5);
		int num_cols = split_view(lr.line(), '\t', &cols[0], 5);
		GK_CHECK(num_cols == 4, value, "Expected 4 columns in BEDGRAPH file but found {}", num_cols);

		// chrom, start, end
		chrom = bld.chrom_names().as_chrom(cols[0]);
		start = as_pos(cols[1]);
		pos_t end = as_pos(cols[2]);

		// Make note of the current span not being aligned with the track resolution, to detect error if we try to extend it.
		is_irregular_span = (end % bld.res() != 0);
		GK_CHECK(start >= 0, value, "Start={} must be non-negative", start);
		GK_CHECK(start < end, value, "Start={} must be less than end={}", start, end);
		GK_CHECK(start % bld.res() == 0, value, "Start {} must align to track resolution", start);

		// If no categorical / numeric value is specified, we'll default to 1.
		float d = parse_float_value(cols[3]);
		T value(clamp_float(d));

		// Fill the data array with as many values as required by the current (start, end) span.
		int num_values = udivup(end - start, bld.res());
		data.resize(data.size() + num_values, value);
		GK_DBASSERT(this->end() == rndup(end, bld.res()));

		// Advance to the next line since we've consumed the current one.
		++lr;
	}

    template <class F>
	void parse_bed(line_reader& lr, const vector<string>& categories, F&& clamp_float)
	{
		auto& bld = *bld_ptr;

		while (!lr.done()) {

			// parse "chrom start end [name] [score] [strand] [...]"
			cols.resize(7);
			int num_cols = split_view(lr.line(), '\t', &cols[0], 7);
			GK_CHECK(num_cols >= 3, value, "Expected at least 3 columns in BED file but found {}", num_cols);

			// chrom, start, end
			chrom = bld.chrom_names().as_chrom(cols[0]);
			start = as_pos(cols[1]);
			pos_t end = as_pos(cols[2]);

			// Make note of the current span not being aligned with the track resolution, to detect error if we try to extend it.
			is_irregular_span = (end % bld.res() != 0);
			GK_CHECK(start >= 0, value, "Start={} must be non-negative", start);
			GK_CHECK(start < end, value, "Start={} must be less than end={}", start, end);
			GK_CHECK(start % bld.res() == 0, value, "Start {} must align to track resolution", start);

			// If no categorical / numeric value is specified, we'll default to 1.
			T value;
			if (categories.empty()) {
				float d = num_cols >= 5 ? parse_float_value(cols[4]) : 1.0f;
				value = T(clamp_float(d));
			} else {
				// Use 'name' column value to find index of category
				GK_CHECK(num_cols >= 4, value, "Expected at least 4 columns in BED file when categories specified");
				category.assign(cols[3]);
				auto i = find(categories.begin(), categories.end(), category);
				GK_CHECK(i != categories.end(), value, "Unrecognized name \"{}\"", category);
				int category_index = (int)(i - categories.begin());
				value = T(category_index);
			}

			// strand (if any)
			strand_t s;
			switch ((num_cols >= 6) && !empty(cols[5]) ? cols[5][0] : '.') {
			case '.': s = pos_strand; GK_CHECK(!bld.stranded(), value, "Expected '+' or '-' for stranded track"); break;
			case '+': s = pos_strand; GK_CHECK(bld.stranded(), value, "Expected '.' for unstranded track"); break;
			case '-': s = neg_strand; GK_CHECK(bld.stranded(), value, "Expected '.' for unstranded track"); break;
			default: GK_THROW(value, "Unrecognized strand value '{}'", cols[5]);
			}

			// Advance to the next line before breaking or looping.
			++lr;

			// If strand matches, then we include this BED line; otherwise skip it and move on to the next.
			if (s == strand) {
				// Fill the data array with as many values as required by the current (start, end) span.
				int num_values = udivup(end - start, bld.res());
				data.resize(data.size() + num_values, value);
				GK_DBASSERT(this->end() == rndup(end, bld.res()));
				break;
			}

			// Reset state and try the next line.
			reset();
		}
	}

	// Return the genomic end position of the current contiguous block of data.
	[[nodiscard]] int end() const
	{
		auto& bld = *bld_ptr;

		GK_CHECK(data.size() % (size_t)bld.dim() == 0, value, "Expected number of values in block to be divisible by dim={}", bld.dim());
		return start + (pos_t)data.size() / bld.dim() * bld.res();
	}

	void fill_data_until(pos_t pos)
	{
		auto& bld = *bld_ptr;

		// Get starting position for this data element.
		GK_CHECK(pos % bld.res() == 0, value, "Position {} must be divisible by resolution={}", pos, bld.res());

		// Get the number of (possibly coarsened) genomic coordinates that preceed this data element.
		size_t size = (size_t)bld.dim() * (pos - start) / bld.res();
		GK_CHECK(size >= data.size(), value, "Position {} must be larger than all previous positions.", pos);

		// Fill all currently-uninitialized positions preceeding this data element with default_value.
		data.resize(size, bld.default_value().f);

		GK_DBASSERT(end() == pos);
	}

	// Sets data for the current interval on the builder object, and resets parsing state to empty.
	void set_data()
	{
		auto& bld = *bld_ptr;
		// Use is_valid_chrom to only dump if chr1..chrM, and not other loci like "chr11_gl000202_random"
		if (chrom && !data.empty()) {
			interval_t interval = interval_t::from_dna0(*chrom, start, end(), strand, bld.refg());
			const T* data_ptr   = (bld.etype() == genome_track::m0) ? nullptr : &data[0];
			bld.set_data(interval, data_ptr);
		}
	}

	void flush()
	{
		// Use is_valid_chrom to only dump if chr1..chrM, and not other loci like "chr11_gl000202_random"
		if (chrom)
			bld_ptr->flush(*chrom, strand);
	}

	void reset()
	{
		chrom.reset();
		// strand = ... deliberately don't reset strand; WIG files don't specify strand, so we need to preserve.
		start = -1;
		is_irregular_span = false;
		is_fixed = false;
		category.clear();
		data.clear();
	}

	bool try_extend(contig_parser& dst)
	{
		// If we've changed chromosome, strand, then we can't extend dst and furthermore
		// should tell the builder to flush the previous (chrom, strand) to disk to save memory.
		if (chrom != dst.chrom || strand != dst.strand)
			return false;  // Start new interval

		// Determine if the new interval should extend the old one in the temporary interval
		// structure that `bld` is building internally.
		//
		//   - If building a fixedStep track, then only extend if start == end, since that's
		//     a common occurrence in WIG files generated by bigWigToWig (i.e. it breaks tracks
		//     up into 1024 chunks, which is counter-productive and we need to merge those blocks).
		//
		//   - If building a variableStep track, then extend if start - end below some threshold.
		//     In other words, so long as the gaps in the WIG file are below some threshold, we
		//     build a big contiguous chunk and sparsify it without regard to how the WIG file
		//     originally broke up the intervals. Again, this avoids arbitrary breakpoints being
		//     carried over from the way bigWigToWig breaks up its own data.
		//
		int dst_end = dst.end();
		GK_CHECK(start >= dst_end, value, "New interval start {} is before end {} of previous interval in WIG file", start, dst_end);
		if (is_fixed) {
			if (start > dst_end)
				return false;  // Start new interval
		} else {
			if (bld_ptr->etype() == genome_track::m0)
				return false;  // For mask etype, any discontinuity must result in new interval, since no data stored
			if (start - dst_end > 8192)
				return false;  // Start new interval
		}

		// Sanity checks
		GK_CHECK(is_fixed == dst.is_fixed, value, "Cannot mix fixedStep and variableStep in single WIG file");
		GK_CHECK(!dst.is_irregular_span, value, "Only the last block on a chromosome may have irregular span");

		// Extend the `dst` data block, rather than starting a new interval. Reset our own state to empty.
		dst.fill_data_until(start);
		dst.data.insert(dst.data.end(), data.begin(), data.end());
		dst.is_irregular_span = is_irregular_span;
		reset();
		return true;  // Old interval has been extended.
	}
};

void genome_track::builder::set_data_from_wig(const char* infile, strand_t strand) {
	if (_h.etype == f32)
		_set_data_from_wig<float>(infile, strand);
	else
		_set_data_from_wig<half_t>(infile, strand);
}

template <typename T>
void genome_track::builder::_set_data_from_wig(const char* infile, strand_t strand)
{
	if (!getenv("GENOMEKIT_QUIET"))
		print("Processing {}\n", infile);

	zline_reader lr(infile);
	contig_parser<T> curr(this, strand);
	contig_parser<T> next(this, strand);
	try {
		curr.parse_wig(lr, _clamp_float);
		while (!lr.done()) {
			next.parse_wig(lr, _clamp_float);
			if (!next.try_extend(curr)) {
				curr.set_data();
				if (curr.chrom != next.chrom || curr.strand != next.strand)
					flush(curr.chrom, curr.strand);
				curr.reset();
				std::swap(curr, next);
			}
		}
	}
	GK_RETHROW("In WIG file: {}:{}", infile, lr.line_num());
	curr.set_data();
	flush(curr.chrom, curr.strand);
}

void genome_track::builder::set_data_from_bedgraph(const char* infile, strand_t strand) {
	if (_h.etype == f32)
		_set_data_from_bedgraph<float>(infile, strand);
	else
		_set_data_from_bedgraph<half_t>(infile, strand);
}

template <typename T>
void genome_track::builder::_set_data_from_bedgraph(const char* infile, strand_t strand)
{
	if (!getenv("GENOMEKIT_QUIET"))
		print("Processing {}\n", infile);

	zline_reader lr(infile);
	contig_parser<T> curr(this, strand);
	contig_parser<T> next(this, strand);
	try {
		curr.parse_bedgraph(lr, _clamp_float);
		while (!lr.done()) {
			next.parse_bedgraph(lr, _clamp_float);
			if (!next.try_extend(curr)) {
				curr.set_data();
				if (curr.chrom != next.chrom || curr.strand != next.strand)
					flush(curr.chrom, curr.strand);
				curr.reset();
				std::swap(curr, next);
			}
		}
	}
	GK_RETHROW("In BedGraph file: {}:{}", infile, lr.line_num());
	curr.set_data();
	flush(curr.chrom, curr.strand);
}

void genome_track::builder::set_data_from_bed(const string& infile) { set_data_from_bed(infile, vector<string>()); }

void genome_track::builder::set_data_from_bed(const string& infile, const vector<string>& categories)
{
	if (_h.etype == f32)
		_set_data_from_bed<float>(infile, categories);
	else
		_set_data_from_bed<half_t>(infile, categories);
}

template <typename T>
void genome_track::builder::_set_data_from_bed(const string& infile, const vector<string>& categories)
{
	GK_CHECK(!finalized(), runtime, "Cannot set data after calling finalize");
	GK_CHECK(_h.dim == 1, value, "Calling set_data_from_bed requires dim=1");
	GK_CHECK(index_size() == 0, runtime,
			 "Cannot call set_data_from_bed after data has already been added to the track");

	int num_categories = (int)categories.size();
	if (num_categories) {
		GK_CHECK(_h.etype >= u1 && _h.etype <= u8, value, "Calling set_data_from_bed with categories requires etypes u1..8");
		GK_CHECK(num_categories <= (1 << (size_t)_encoding.bits_per_encoded_datum), value, "Too many categories ({}) for encoding type", num_categories);
	}

	if (!getenv("GENOMEKIT_QUIET"))
		print("Processing {}\n", infile);

	// For BED files, we don't know strand until we've parsed the line itself.
	// The most expedient implementation is therefore to make two passes of the file,
	// one to collect positive-strand only, and one to collect negative-strand only.
	for (std::underlying_type_t<strand_t> ord_strand = _h.stranded ? as_ordinal(neg_strand) : as_ordinal(pos_strand); ord_strand < num_strand; ++ord_strand) {
		zline_reader lr(infile);
		const auto strand = strand_t { ord_strand };
		contig_parser<T> curr(this, strand);
		contig_parser<T> next(this, strand);
		try {
			curr.parse_bed(lr, categories, _clamp_float);
			while (!lr.done()) {
				next.parse_bed(lr, categories, _clamp_float);
				if (!next.try_extend(curr)) {
					curr.set_data();
					if (curr.chrom != next.chrom || curr.strand != next.strand)
						flush(curr.chrom, curr.strand);
					curr.reset();
					std::swap(curr, next);
				}
			}
		}
		GK_RETHROW("In BED file: {}:{}", infile, lr.line_num());
		curr.set_data();
		flush(curr.chrom, curr.strand);
		curr.reset();
	}
}

void genome_track::builder::finalize()
{
	constexpr auto app    = std::ios::app;
	constexpr auto binary = std::ios::binary;

	GK_CHECK(!finalized(), runtime, "Cannot finalize: already called finalize before.");

	_h.num_tracks = size(_tracks);
	// Copy the encoding dict into the header
	_h.dict = _encoding.dict;
	const bool encoded_data = encoded_sizes(1).bytes > 0;

	// Finally, go back to the start and write the header.
	binary_file file(_outfile, "w");
	file.write(_h);
	auto curr_offset = file.tell();
	curr_offset += size(_tracks) * sizeof(header_track_info_t);
	for (auto& [chrom_strand, ti] : _tracks) {
		auto [chrom, strand] = chrom_strand;
		header_track_info_t info{
			// clang 14 nonsensical complaint of long long -> unsigned long long
			.data_offset  = int_cast<uint64_t>(curr_offset),
			.index_offset = curr_offset + ti.data_size(),
			.chrom        = chrom,
			.strand       = strand,
		};
		file.write(info);
		curr_offset = info.index_offset + ti.index_size(encoded_data);
	}
	file.close();

	file.open(_outfile, "a");
	std::ofstream appended{_outfile, binary | app};
	for (auto& [chrom_strand, ti] : _tracks) {
		if (!ti.flushed()) {
			ti.flush(file, [&](int length) { return encoded_sizes(length); });
			_finalized = true;  // flush is destructive so can't repeat afterwards
			file.flush();
			continue;
		}

		const auto [chrom, strand] = chrom_strand;
		const auto track_segment_filename  = partial_filename(_outfile, chrom, strand);
		GK_DBASSERT(ti.data_size() + ti.index_size(encoded_data) == std::filesystem::file_size(track_segment_filename));

		appended << std::ifstream{track_segment_filename, binary}.rdbuf();
		GK_CHECK(appended, runtime, "Error merging track data.");
		appended.flush();
		std::filesystem::remove(track_segment_filename);
		_finalized = true;  // remove is destructive so can't repeat afterwards
	}
	GK_CHECK(appended, runtime, "Error merging track data.");
	GK_DBASSERT(curr_offset == appended.tellp() || curr_offset == file.tell());
	file.close();
}

const chrom_names_t& genome_track::builder::chrom_names() const { return _chrom_names; }

END_NAMESPACE_GK
