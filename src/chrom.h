/*
Copyright (C) 2016-2023 Deep Genomics Inc. All Rights Reserved.
*/
#pragma once
#ifndef __GENOME_KIT_CHROM_H__
#define __GENOME_KIT_CHROM_H__

#include "defines.h"
#include "file.h"
#include "refg.h"
#include "strutil.h"

#include <cstdint>
#include <string>
#include <string_view>
#include <unordered_map>

BEGIN_NAMESPACE_GK

enum class chrom_t : std::uint32_t {};

class chrom_names_t {
public:
	using mapping_t       = std::unordered_map<chrom_t, std::string>;
	using alias_mapping_t = string_map<std::string, chrom_t>;

	explicit chrom_names_t(refg_t refg, std::string refg_name, mapping_t mapping, alias_mapping_t aliases);

	chrom_t          as_chrom(std::string_view chr) const;
	std::string_view chrom_as_sv(chrom_t chr) const;

	mapping_t::const_iterator begin() const;
	mapping_t::const_iterator end() const;

	refg_t refg() const;
	const std::string& refg_name() const;

private:
	mapping_t       _name_by_chrom;
	alias_mapping_t _aliases;
	std::string     _refg_name;
	refg_t          _refg;
};

const chrom_names_t& get_chrom_names(refg_t refg, std::string_view data_dir = default_data_directory);

END_NAMESPACE_GK

#endif
