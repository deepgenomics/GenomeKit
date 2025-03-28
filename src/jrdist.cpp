/*
Copyright (C) 2016-2023 Deep Genomics Inc. All Rights Reserved.
*/
#include "jrdist.h"
#include "util.h"
#include <cstring>
#include <map>
#include <numeric>
#include <utility>

BEGIN_NAMESPACE_GK

using std::map;

const unsigned short c_rdist_sig = 0x0deb;
const unsigned short c_rdist_ver = 0x0002;
// versions:
//   0001: initial format
//   0002: Restructured to support all chromosomes and more species
//   <--- INSERT VERSION CHANGE SUMMARIES HERE

////////////////////////////////////////////////
// jrdist_t
////////////////////////////////////////////////

jrdist_t::jrdist_t(index_t              src, const jrdist_table& table) { unpack_from(table[src], table); }
jrdist_t::jrdist_t(const packed_jrdist& src, const jrdist_table& table) { unpack_from(src,        table); }

INLINE void jrdist_t::unpack_from(const packed_jrdist& src, const jrdist_table& table)
{
	this->as_interval() = src.as_interval();
	int count_size;
	_num_counts = (int)src.num_counts;
	_num_reads  = src.num_reads;
	count_size = *(const char*)(table.aux() + src.aux);
	_counts    = (const unsigned char*)(table.aux() + src.aux + 1); align_ptr(_counts, count_size);
	_ushifts   = (const unsigned char*)(_counts + count_size*_num_counts);
	_flags     = (const unsigned char*)(_ushifts + _num_counts);
}

const jrcount_t jrdist_t::operator[](int i) const
{
	GK_CHECK(i < _num_counts, index, "Invalid read index {}", i);
	unsigned char flags_byte = _flags[i / 4];
	unsigned char strand_bit = 1 << (2*(i % 4) + 0);
	unsigned char sign_bit   = 1 << (2*(i % 4) + 1);
	unsigned char count_size = _counts[-1]; // Guaranteed this byte has size of counts
	unsigned count;
	switch (count_size) {
	case 1: count = _counts[i]; break;
	case 2: count = ((unsigned short*)_counts)[i]; break;
	case 4: count = ((unsigned*)_counts)[i]; break;
	default: GK_THROW2(value, "Invalid count size byte detected -- corrupt?");
	}
	return jrcount_t(count,
	                 flags_byte & sign_bit ? -_ushifts[i] : _ushifts[i],
	                 flags_byte & strand_bit ? pos_strand : neg_strand);
}

////////////////////////////////////////////////
// read_distribution
////////////////////////////////////////////////

struct read_distributions::builder::jrdist_collector { NOCOPY(jrdist_collector)
	jrdist_collector() = default;

	vector<unsigned char>  counts8;
	vector<unsigned short> counts16;
	vector<unsigned>       counts32;
	vector<unsigned char>  ushifts;
	vector<unsigned char>  flags;

	void operator()(unsigned count, int shift, int strand)
	{
		size_t i = counts32.size();
		counts32.push_back(count);
		ushifts.push_back((unsigned char)(shift < 0 ? -shift : shift));
		if (i % 4 == 0)
			flags.push_back(0);
		if (strand)    flags[i/4] |= 1 << (2*(i % 4) + 0);
		if (shift < 0) flags[i/4] |= 1 << (2*(i % 4) + 1);
	}

	void reset()
	{
		counts8.clear();
		counts16.clear();
		counts32.clear();
		ushifts.clear();
		flags.clear();
	}

	void add_aux(interval_table<packed_jrdist>::builder& juncs, unsigned char count_size)
	{
		// Fill bytes with count_size until we're aligned with count_size
		// This is because the byte preceding the 'counts' array in aux memory
		// must contain the size of each element in that array, so that the
		// code that unpacking code knows what size was
		do { juncs.add_aux(&count_size, &count_size+1); } while (juncs.curr_aux() % count_size);

		// Then store the counts differently depending on the count_size used when
		// they were originally being built
		switch (count_size) {
			case 1: { // Store as chars
				counts8.clear();
				counts8.reserve(counts32.size());
				std::transform(std::cbegin(counts32), std::cend(counts32), back_inserter(counts8),
							   [](auto x) { return uint8_t(x); });
				juncs.add_aux(counts8);
				break;
			}
			case 2: { // Store as shorts
				counts16.clear();
				counts16.reserve(counts32.size());
				std::transform(std::cbegin(counts32), std::cend(counts32), back_inserter(counts16),
							   [](auto x) { return uint16_t(x); });
				juncs.add_aux(counts16);
				break;
			}
			case 4: {  // Store as ints
				juncs.add_aux(counts32);
				break;
			}
			default: GK_UNREACHABLE2();
		}

		// After adding the counts array, add the shifts and flags
		juncs.add_aux(ushifts);
		juncs.add_aux(flags);
	}
};

//////////////////////////////////////////////////////////////////////////////////

const jrdist_table& read_distributions::juncs() const
{
	ensure_open();
	return _juncs;
}

void read_distributions::set_source(string sourcefile)
{
	GK_CHECK2(!is_open(), runtime, "Cannot set source when file already open.");
	_sourcefile = std::move(sourcefile);
}

void read_distributions::open_on_demand() const { const_cast<read_distributions*>(this)->open(); }
void read_distributions::close() { _fmap.close(); }

void read_distributions::open()
{
	GK_CHECK2(!is_open(), runtime, "read_distributions::open() already opened");
	GK_CHECK2(!_sourcefile.empty(), value, "No file was specified");

	// Resolve the source file, possibly downloading it
	//_sourcefile = realize_datafile_path(_sourcefile); // TODO: restore this

	// Memory map the source file
	_fmap.open(_sourcefile);

	// Read the DGJUNCREADDIST file signature
	unsigned short sig, ver;
	_fmap.read(sig);
	_fmap.read(ver);
	GK_CHECK(sig == c_rdist_sig, file, "Expected valid RDIST file signature {:x} but found {:x}.", c_rdist_sig, sig);
	GK_CHECK(ver == c_rdist_ver, file, "Expected RDIST file version {:x} but found {:x}.", c_rdist_ver, ver);

	// Read number of total reads
	_fmap.read(_num_reads);

	// Read the junction table itself
	_juncs.load(_fmap);
}

int read_distributions::rdist_version() { return c_rdist_ver; }

read_distributions::builder::jrdist_entry::~jrdist_entry()
{
	clear();
}

void read_distributions::builder::jrdist_entry::clear()
{
	if (_mode != list_mode) {
		free(_counts());
		_counts() = nullptr;
	}
	_mode = list_mode;
	_size = 0;
}

template <typename T>
void read_distributions::builder::jrdist_entry::inc_count(unsigned char overhang, strand_t strand, int direction)
{
	// If the overhang is greater than the current count array size
	// (in either direction), then we must grow the count array.
	if (overhang > _size) {
		int new_size = overhang <=  75 ?  75 :
		               overhang <= 100 ? 100 : 255;
		//int new_size = overhang+1;
		int new_array_size = (2*new_size+1)*2;
		void* new_counts = malloc(sizeof(T)*new_array_size); // array is new_size in each direction, plus empty counter slot at 0
		memset(new_counts, 0, sizeof(T)*new_array_size);
		if (_counts()) {
			// Copy the old counts into the middle of this new, larger counts array
			int old_size = _size;
			int old_array_size = (2*old_size+1)*2;  // 2*size+1 to have symmetry around 0 shift, and all *2 because two strands
			memcpy((char*)new_counts + sizeof(T)*(new_size-old_size)*2, _counts(), sizeof(T)*old_array_size);
			free(_counts());
		}
		_counts() = new_counts;
		_size = new_size;
	}

	// Find the counter based on the signed shift of this overhang
	// If we haven't overflowed, then just increment and return
	int size = _size;
	int shift = direction*(int)overhang;
	T& count = ((T*)_counts())[(size+shift)*2 + (strand == pos_strand ? 1 : 0)];
	if (count < int_traits<T>::max) {
		++count;  // This is the fast path! (and by far the most frequent path)
		return;
	}

	// If we're about to overflow, then promote the counter type and transfer the
	// old counts into the new array
	using P = typename int_traits<T>::promoted_type;
	int array_size = (2*size+1)*2;
	void* new_counts = malloc(sizeof(P)*array_size); // array is new_size in each direction, plus empty counter slot at 0
	for (int i = 0; i < array_size; ++i)
		((P*)new_counts)[i] = ((T*)_counts())[i];
	free(_counts());
	_counts() = new_counts;
	_mode++;  // Update the mode to reflect an upgrade to a new counter type

	// Now increment the counter on the promoted type
	P& new_count = ((P*)new_counts)[(size+shift)*2 + (strand == pos_strand ? 1 : 0)];
	new_count++;
}

template <>
void read_distributions::builder::jrdist_entry::inc_count<unsigned long long>(unsigned char overhang, strand_t strand, int direction)
{
	GK_THROW2(not_implemented, "A read count of > 4B reads was detected; this is not currently supported.");
}

// An array of pointers to the inc_count<T> specializations that handle _mode=1,2,3
void (read_distributions::builder::jrdist_entry::*read_distributions::builder::jrdist_entry::inc_count_fn[num_modes])(
	unsigned char overhang, strand_t strand, int direction)
	= {
		  nullptr,                                                                  // [0] = list_mode
		  &read_distributions::builder::jrdist_entry::inc_count<unsigned char>,     // [1] = count8_mode
		  &read_distributions::builder::jrdist_entry::inc_count<unsigned short>,    // [2] = count16_mode
		  &read_distributions::builder::jrdist_entry::inc_count<unsigned int>,      // [3] = count32_mode
	  };

INLINE void read_distributions::builder::jrdist_entry::add(const jralign_t read)
{
	if (_mode == list_mode) {
		// If we have space in our list, just copy the value for now, rather than
		// creating an array of counts (which wastes a lot of memory for junctions
		// with few counts)
		if (_size < max_list_size) {
			_list[_size++] = read;
			return;
		}

		// Otherwise we've run out of space in the list, so convert the list to a count
		// array and fall through to include this new read in the counts
		// This will convert existing reads into counts by calling add() recursively,
		// but only after setting up count8_mode.
		ensure_counts_mode();
	}

	// Count array is now guaranteed to be large enough.
	(this->*(inc_count_fn[_mode]))(read.left,  read.strand, -1);
	(this->*(inc_count_fn[_mode]))(read.right, read.strand,  1);  // note: _mode may have upgraded since line above
}

void read_distributions::builder::jrdist_entry::ensure_counts_mode()
{
	if (_mode != list_mode)
		return;

	// Read counts are currently represented as an explicit list, rather than
	// as a count array, so convert them. Since _counts() pointer memory is
	// inside bytes 4-12 of _list, we must copy list items 1, 2, 3.
	jralign_t _list1 = _list[1];
	jralign_t _list2 = _list[2];
	jralign_t _list3 = _list[3];
	int list_size = _size;
	_mode = count8_mode;
	_size = 0;
	_counts() = nullptr;
	if (list_size > 0) add(_list[0]);  else return;
	if (list_size > 1) add(_list1);    else return;
	if (list_size > 2) add(_list2);    else return;
	if (list_size > 3) add(_list3);    else return;
	for (int i = 4; i < list_size; ++i)
		add(_list[i]);
}

void read_distributions::builder::jrdist_entry::collect(jrdist_collector& collector)
{
	ensure_counts_mode();
	collector.reset();
	switch (_mode) {
	case count8_mode:  collect_impl<unsigned char >(collector); break;
	case count16_mode: collect_impl<unsigned short>(collector); break;
	case count32_mode: collect_impl<unsigned int  >(collector); break;
	default: GK_UNREACHABLE2();
	}
}

template <typename T>
void read_distributions::builder::jrdist_entry::collect_impl(jrdist_collector& collector)
{
	if (!_counts()) // Zero reads were retained across this junction, possibly filtered out by min_overhang
		return;
	for (int shift = -_size; shift <= _size; ++shift) {
		for (int strand = 0; strand < 2; ++strand) {
			T count = ((T*)_counts())[(_size+shift)*2+strand];
			if (count > 0)
				collector(count, shift, strand); // Finally, collect the count in its final representation in the file
		}
	}
}

read_distributions::builder::builder(const char* outfile)
: _file(outfile, "w")
, _interval_filter{[&](interval_t i) {
	GK_CHECK2(!_refg || i.refg == *_refg, value, "Cannot filter {} for {}", i, get_refg_registry().refg_as_sv(*_refg));
}}
{
}

void read_distributions::builder::add(const char* infile)
{
	bool verbose = getenv("GENOMEKIT_QUIET") == nullptr;
	if (verbose)
		print("Adding {} ... ", infile);

	junction_read_alignments ralign;
	ralign.set_source(infile);
	ralign.open();

	if (_stranded == stranded::unknown) {
		_stranded = ralign.juncs().stranded() ? stranded::yes : stranded::no;
	} else {
		GK_CHECK2((_stranded == stranded::yes) == ralign.juncs().stranded(), value, "Prior jraligns were {}, but {} is {}.",
				 _stranded == stranded::yes ? "stranded" : "unstranded", infile,
				 ralign.juncs().stranded() ? "stranded" : "unstranded");
	}

	// Pull out the refg field from the first junction -- all are assumed to be identical
	if (ralign.juncs().size() > 0) {
		auto junc_refg = ralign.juncs()[0].refg;
		if (!_refg) {

			_refg = junc_refg;
			_interval_filter.validate();
		} else {
			GK_CHECK2(*_refg == junc_refg, file, "Reference genome mismatch between different input files.");
		}

		// Loop over all junctions, and all samples, and build a read distribution for each sample
		for (int j = 0; j < ralign.juncs().size(); ++j) {
			jraligns_t junc(j, ralign.juncs());

			if (!get_interval_filter().filter(junc))
				continue;

			jrdist_entry& rdist = _juncs[junc];
			for (unsigned r = 0; r < junc.num_reads(); ++r) {
				const jralign_t& align = junc[r];

				if (_filter && !_filter(ralign, j, r)) {
					continue;
				}

				// If a match on either side of the intron is less than min_overhang, drop it.
				if (align.left < _min_overhang || align.right < _min_overhang)
					continue;

				rdist.add(align);
			}
		}
	}
	if (verbose)
		print("done\n", infile);
}

void read_distributions::builder::finalize()
{
	long long num_reads_total = 0;
	jrdist_collector collector;
	jrdist_table::builder juncs(_stranded == stranded::yes);

	// Build the final table
	for (auto& [interval, rdist] : _juncs) {
		// Convert the count arrays into their final representation (3 parallel arrays: counts, ushift, flags)
		rdist.collect(collector);

		// Calculate number of reads per sample
		unsigned num_reads = std::accumulate(collector.counts32.begin(), collector.counts32.end(), 0u)/2;
		if (num_reads < _min_reads)
			continue;

		// Fill record for this junction that will be dumped to disk
		packed_jrdist junc;
		junc.as_interval() = interval; // Copy the interval fields
		junc.num_counts = int_cast<unsigned short>(collector.counts32.size());
		junc.num_reads  = num_reads;
		junc.aux = int_cast<offset_t>(juncs.curr_aux());
		juncs.add_elem(junc);  // Add it to the table

		// Fill the aux data that the record points to.
		collector.add_aux(juncs, rdist.count_size());

		num_reads_total += num_reads;
		rdist.clear();
	}

	// Dump signature
	_file.write(c_rdist_sig);
	_file.write(c_rdist_ver);

	// Dump the rest of the file and close it.
	_file.write(num_reads_total);
	juncs.dump(_file);
	_file.close();
}

END_NAMESPACE_GK
