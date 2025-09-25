/*
Copyright (C) 2016-2023 Deep Genomics Inc. All Rights Reserved.
*/
#include "genome_dna.h"
#include "strutil.h"
#include "util.h"
#include <algorithm>
#include <filesystem>
#include <utility>

BEGIN_NAMESPACE_GK
using std::lower_bound;
using std::max;

static const unsigned c_twobit_sig_lilend = 0x1A412743;
static const unsigned c_twobit_sig_bigend = 0x4327411A;

static const char g_decode_nuc[4] = { 'T', 'C', 'A', 'G' }; // twobit nucleotide index ordering

void genome_dna::set_source(string sourcefile)
{
	// Try not to put anything here that would throw in a newly constructed object.
	// We want genome_t constructor to not throw, even if it's initialized outside
	GK_CHECK(!_fmap.is_open(), runtime, "Cannot set source when file already open.");
	_sourcefile = std::move(sourcefile);
}

void genome_dna::open()
{
	GK_CHECK(!is_open(), runtime, "genome_dna::open() already opened");
	GK_CHECK(endswith(_sourcefile, ".2bit"), value, "Unrecognized file extension for '{}'", _sourcefile);

	auto path               = std::filesystem::path{_sourcefile};
	_refg_name              = path.filename().stem().string();
	// this comes last: watch out for buffer reuse here on the path object
	const auto& data_dir = path.remove_filename().string();
	const auto& refgs       = get_refg_registry(data_dir);
	_refg                   = refgs.as_refg(_refg_name);
	const auto& chrom_names = get_chrom_names(_refg, data_dir);

	// Resolve the source file, possibly downloading it
	try {
		_sourcefile = resolve_datafile_path(_sourcefile);
	} catch (const value_error& e) {
		print("{}\n", e.what());
	}

	// Open file, both memory mapped and with regular file handle
	mmap_file fmap(_sourcefile);
	binary_file f(_sourcefile, "rb");

	// Read header
	header_t header;
	f.read(header);

	if (header.signature != c_twobit_sig_lilend) {
		if (header.signature == c_twobit_sig_bigend)
			GK_THROW(value, "Byte swapped format not supported by this implementation.");
		GK_THROW(value, "Unrecognized 2bit signature in file.");
	}
	GK_CHECK(header.version == 0, value, "2bit version unsupported (maybe 64-bit offsets): {}", header.version);

		// Read names and offsets
	for (uint32_t i = 0; i < header.seq_count; ++i) {
		// Read length of name
		uint8_t name_size;
		f.read(name_size);

		// Read name
		const string_view name{fmap.as_ptr<const char>(f.tell()), name_size};
		try {
			const auto chrom = chrom_names.as_chrom(name);
			f.set_seek(f.tell() + name_size);

			// Read offset
			uint32_t offset;
			f.read(offset);

			_seqrecs[chrom] = {offset};
		}
		GK_RETHROW("In 2bit file: '{}'.", _sourcefile);
	}

	std::swap(_fmap, fmap);
}

void genome_dna::close()
{
	_fmap.close();
	_sourcefile.clear();
	_seqrecs.clear();
}

void genome_dna::open_on_demand() const
{
	// TODO: acquire lock here and check _file for NULL again

	// Once we're inside the lock it's safe to use ncthis
	auto* ncthis = const_cast<genome_dna*>(this);
	ncthis->open();

	// TODO: release lock here (implicitly, when falls out of scope)
}

//////////////////////////////////////////

void genome_dna::operator()(const interval_t& interval, char* _dst, const bool allow_outside_chromosome) const
{
	UNPACK_INTERVAL(interval);

	GK_CHECK(interval.size() >= 0, index, "genome_dna::get({}) invalid range.", interval);
	GK_CHECK(allow_outside_chromosome || interval.start() >= 0, index, "genome_dna::get({}) index outside range.", interval);

	ensure_open();
	GK_CHECK(_refg == refg, value, "genome_dna::get({}) cannot be used on {}.", interval,
			 _refg_name);

	const auto it = _seqrecs.find(chrom);
	GK_CHECK(it != end(_seqrecs), value, "genome_dna::get({}) has no such chromosome on {}.", interval,
			 _refg_name);
	it->second.ensure_open(_fmap);
	const auto rec = it->second;

	const uint32_t a = max(interval.start(), 0);
	const uint32_t b = min((uint32_t)interval.end(), rec.num_dna_bases);

	GK_CHECK(allow_outside_chromosome || interval.end() < 0 || ((uint32_t)interval.end()) <= rec.num_dna_bases, index, "genome_dna::get({}) index larger than chromosome size {}.",
			 interval, rec.num_dna_bases);
	if (allow_outside_chromosome) {
		GK_CHECK(interval.end() > 0 && a < rec.num_dna_bases - 1, index, "genome_dna::get({}) entire range outside chromosome.", interval);
	}

	const auto size = b - a;
	if (size == 0)
		return;

	constexpr uint32_t num_per_dword = 4 * sizeof(int);
	const uint32_t     first_dword   = a / num_per_dword;                        // round down
	const uint32_t     last_dword    = (b + num_per_dword - 1) / num_per_dword;  // round up
	const uint32_t     num_dwords    = last_dword - first_dword;

	const uint32_t* RESTRICT src;

	// DNA was pre-mapped to memory, so just point to it
	src = rec.dna + first_dword;
	char* RESTRICT dst = _dst;

	auto     offset = 0;
	if (allow_outside_chromosome && interval.start() < 0) {
		offset = -interval.start();
		memset(dst, 'N', offset);
	}

	if (num_dwords <= 1) {

		const uint32_t dword = bswap32(*src);  // compiles to bswap instruction
		const uint32_t c     = a - num_per_dword * first_dword;
		for (uint32_t i = 0; i < size; ++i)
			dst[i+offset] = g_decode_nuc[(dword >> (32-2*(i+c+1))) & 3];

	} else {

		// Handle low-order bits in first dword, if any
		uint32_t i = 0;
		uint32_t c = a & (num_per_dword - 1);
		if (c > 0) {
			const uint32_t dword = bswap32(*src++);
			const uint32_t d     = num_per_dword;
			for (; i < d-c; ++i) {
				dst[i+offset] = g_decode_nuc[(dword >> (32-2*(i+c+1))) & 3];
			}
		}

		// Handle chunks that span entire dword
		uint32_t d = b & (num_per_dword - 1);
		uint32_t n = size - d;
		for (; i < n; i += num_per_dword) {
			uint32_t dword = bswap32(*src++);
			char* dst_ptr = dst + i + offset;
			dst_ptr[15] = g_decode_nuc[dword & 3]; dword >>= 2;
			dst_ptr[14] = g_decode_nuc[dword & 3]; dword >>= 2;
			dst_ptr[13] = g_decode_nuc[dword & 3]; dword >>= 2;
			dst_ptr[12] = g_decode_nuc[dword & 3]; dword >>= 2;
			dst_ptr[11] = g_decode_nuc[dword & 3]; dword >>= 2;
			dst_ptr[10] = g_decode_nuc[dword & 3]; dword >>= 2;
			dst_ptr[ 9] = g_decode_nuc[dword & 3]; dword >>= 2;
			dst_ptr[ 8] = g_decode_nuc[dword & 3]; dword >>= 2;
			dst_ptr[ 7] = g_decode_nuc[dword & 3]; dword >>= 2;
			dst_ptr[ 6] = g_decode_nuc[dword & 3]; dword >>= 2;
			dst_ptr[ 5] = g_decode_nuc[dword & 3]; dword >>= 2;
			dst_ptr[ 4] = g_decode_nuc[dword & 3]; dword >>= 2;
			dst_ptr[ 3] = g_decode_nuc[dword & 3]; dword >>= 2;
			dst_ptr[ 2] = g_decode_nuc[dword & 3]; dword >>= 2;
			dst_ptr[ 1] = g_decode_nuc[dword & 3]; dword >>= 2;
			dst_ptr[ 0] = g_decode_nuc[dword & 3];
		}

		// Handle any remaining high-order bits in final dword, if any
		if (i < size) {
			const uint32_t dword = bswap32(*src);
			do {
				dst[i+offset] = g_decode_nuc[(dword >> (32-2*(i-n+1))) & 3];
			} while (++i < size);
		}
	}

	if (allow_outside_chromosome && interval.end() > 0 && ((uint32_t)interval.end()) > rec.num_dna_bases) {
		memset(dst + rec.num_dna_bases - a + offset, 'N', interval.end() - rec.num_dna_bases);
	}

	if (rec.num_nblocks > 0) {
		const auto nblock_sizes = rec.nblock_starts + rec.num_nblocks;
		// Find the first block position (start or end) that is >= start
		auto       i_nblock     = lower_bound(rec.nblock_starts, nblock_sizes, a) - rec.nblock_starts;
		if (i_nblock != 0 && rec.nblock_starts[i_nblock - 1] + nblock_sizes[i_nblock - 1] >= a) {
			--i_nblock;
		}

		// Fill all the nblock intervals between our dna [a:b] with 'n', stopping when we've gone past 'last'
		for (; i_nblock < rec.num_nblocks; ++i_nblock) {
			auto c = rec.nblock_starts[i_nblock];
			auto d = c + nblock_sizes[i_nblock];
			if (c >= b) break;
			if (c < a) c = a;
			if (d > b) d = b;
			if (c < d)
				memset(dst + c - a, 'N', d - c);
		}
	}

	if (_want_mask && rec.num_masks > 0) {
		const auto mask_sizes = rec.mask_starts + rec.num_masks;
		auto       i_mask     = lower_bound(rec.mask_starts, mask_sizes, a) - rec.mask_starts;
		if (i_mask != 0 && rec.mask_starts[i_mask - 1] + mask_sizes[i_mask - 1] >= a) {
			--i_mask;
		}

		for (; i_mask < rec.num_masks; ++i_mask) {
			auto c = rec.mask_starts[i_mask];
			auto d = c + mask_sizes[i_mask];
			if (c >= b) break;
			if (c < a) c = a;
			if (d > b) d = b;
			if (c < d)
				for (uint32_t i = 0; i < d - c; ++i)
					dst[c - a + i] += ('a'-'A'); // conver char to lowercase
		}
	}

	if (strand == neg_strand)
		reverse_complement(dst, b-a);
}

dnastr genome_dna::operator()(const interval_t& c, const bool allow_outside_chromosome) const
{
	dnastr dna(max(c.size(), 0));
	(*this)(c, dna.data(), allow_outside_chromosome);
	return dna;
}

/////////////////////////////////////////////////

void genome_dna::seqrec_t::ensure_open(const mmap_file& fmap) const
{
	if (num_dna_bases != 0)
		return;  // already opened

	auto ncthis = as_mutable(this);
	auto curr   = offset;

	std::memcpy(&ncthis->num_dna_bases, fmap.as_ptr<void>(curr), sizeof(num_dna_bases));
	curr += sizeof(num_dna_bases);

	std::memcpy(&ncthis->num_nblocks, fmap.as_ptr<void>(curr), sizeof(num_nblocks));
	curr += sizeof(num_nblocks);
	ncthis->nblock_starts = fmap.as_ptr<std::decay_t<decltype(*nblock_starts)>>(curr);
	curr += 2 * ncthis->num_nblocks * sizeof(*nblock_starts);

	std::memcpy(&ncthis->num_masks, fmap.as_ptr<void>(curr), sizeof(num_masks));
	curr += sizeof(num_masks);
	ncthis->mask_starts = fmap.as_ptr<std::decay_t<decltype(*mask_starts)>>(curr);
	curr += 2 * ncthis->num_masks * sizeof(*mask_starts);

	curr += sizeof(uint32_t);
	ncthis->dna = fmap.as_ptr<std::decay_t<decltype(*ncthis->dna)>>(curr);
}

string default_dna_sourcefile(string_view refg_name, string_view data_dir)
{
	return prepend_dir(data_dir, std::format("{}.2bit", refg_name));
}

END_NAMESPACE_GK
