/*
Copyright (C) 2016-2023 Deep Genomics Inc. All Rights Reserved.
*/
#pragma once
#ifndef __GENOME_KIT_GENOME_DNA__
#define __GENOME_KIT_GENOME_DNA__

#include "biostr.h"
#include "interval.h"
#include "file.h"
#include <unordered_map>

BEGIN_NAMESPACE_GK

class genome_dna {
public:
	void set_source(string sourcefile);
	void open();
	void close();
	INLINE bool is_open()     const { return _fmap.is_open(); }
	INLINE void ensure_open() const { if (!is_open()) open_on_demand(); }  // Fast check if file / indices already loaded

	dnastr operator()(const interval_t& c, bool allow_outside_chromosome=true) const;
	void   operator()(const interval_t& c, char* dst, bool allow_outside_chromosome) const;

	INLINE void enable_mask()  { _want_mask = true; }
	INLINE void disable_mask() { _want_mask = false; }  // Speeds things up if we don't bother setting masked regions as lower-case .
	INLINE const string& source() const { return _sourcefile; }
	INLINE refg_t        refg() const { ensure_open(); return _refg; }
	INLINE const string& refg_name() const { ensure_open(); return _refg_name; }

private:
	struct header_t {
		uint32_t signature;
		uint32_t version;
		uint32_t seq_count;
		uint32_t reserved;
	};

	struct seqrec_t {
		uint32_t        offset{};
		uint32_t        num_dna_bases{};
		uint32_t        num_nblocks{};
		uint32_t        num_masks{};
		const uint32_t* dna{};
		const uint32_t* nblock_starts{};
		const uint32_t* mask_starts{};

		void ensure_open(const mmap_file& fmap) const;
	};

	using seqrecs_t = std::unordered_map<chrom_t, seqrec_t>;

	void open_on_demand() const;

	// These are mutable so that they may be modified when the 2bit file
	// is opened on-demand by the ensure_open() method.
	mmap_file   _fmap;
	seqrecs_t   _seqrecs;

	string _refg_name;
	refg_t _refg{};  // not optional for simplicity. set_source guaranteed to be called
	bool   _want_mask{};
	string _sourcefile;
	// TODO mutable mutex member needed here for open_on_demand to be thread safe
};

string default_dna_sourcefile(string_view refg_name, string_view data_dir = default_data_directory);

END_NAMESPACE_GK

#endif // __GENOME_KIT_GENOME_DNA__
