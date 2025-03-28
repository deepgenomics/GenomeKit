/*
Copyright (C) 2016-2023 Deep Genomics Inc. All Rights Reserved.
*/
#include "chrom.h"

#include "file.h"
#include "format.h"
#include "gk_assert.h"
#include "strutil.h"
#include "util.h"

#include <algorithm>
#include <filesystem>
#include <iterator>
#include <optional>
#include <utility>
#include <vector>

BEGIN_NAMESPACE_GK

chrom_names_t::chrom_names_t(refg_t refg, std::string refg_name, mapping_t mapping, alias_mapping_t aliases)
: _name_by_chrom{std::move(mapping)}
, _aliases{std::move(aliases)}
, _refg_name{std::move(refg_name)}
, _refg{refg}
{
}

chrom_t chrom_names_t::as_chrom(std::string_view chr) const
{
	// use a consistent hash to allow users to "port" a chromosome between refgs (eg., patched references)
	chrom_t tag{fnv1a_hash32(chr)};
	if (_name_by_chrom.contains(tag)) {
		return tag;
	}

	const auto it = _aliases.find(chr);
	GK_CHECK(it != std::cend(_aliases), value, "Contig not found: '{0}' on '{1}.chrom.sizes/{1}.chromAlias.txt'.", chr,
			 _refg_name);
	return it->second;
}
std::string_view chrom_names_t::chrom_as_sv(chrom_t chr) const
{
	const auto it = _name_by_chrom.find(chr);
	GK_CHECK(it != std::cend(_name_by_chrom), value, "Chromosome tag not found: '{}'.", as_ordinal(chr));
	return it->second;
}

chrom_names_t::mapping_t::const_iterator chrom_names_t::begin() const { return std::cbegin(_name_by_chrom); }

chrom_names_t::mapping_t::const_iterator chrom_names_t::end() const { return std::cend(_name_by_chrom); }

refg_t chrom_names_t::refg() const { return _refg; }

const std::string& chrom_names_t::refg_name() const { return _refg_name; }

const chrom_names_t& get_chrom_names(refg_t refg, std::string_view data_dir)
{
	static std::unordered_map<refg_t, std::optional<chrom_names_t>> name_registry;

	auto& chrom_names = name_registry[refg];
	if (!chrom_names) {
		std::string refg_name{get_refg_registry(data_dir).refg_as_sv(refg)};
		std::string alias_path{std::format("{}.chromAlias.txt", refg_name)};

		try {
			auto resolved = resolve_datafile_path(prepend_dir(data_dir, alias_path));
			if (std::filesystem::exists(resolved)) {
				alias_path = std::move(resolved);
			}
		} catch (const value_error&) {
		}

		chrom_names_t::mapping_t       mapping;
		chrom_names_t::alias_mapping_t aliases;
		// must pre-processed to the UCSC names for matching the 2bit
		if (std::filesystem::exists(alias_path)) {
			std::vector<std::string_view> names;
			for (line_reader lr{alias_path}; !lr.done(); ++lr) {
				const auto line = lr.line();
				if (line.starts_with("#"))
					continue;

				split_view(line, '\t', names);
				GK_CHECK(std::size(names) > 1, value, "Alias file needs more than one column on names: {}.",
						 alias_path);
				const chrom_t tag{fnv1a_hash32(names[0])};
				const auto [it, inserted] = mapping.insert({tag, std::string(names[0])});
				GK_CHECK(inserted, runtime, "hash collision: '{}' on '{}'", it->second, refg_name);

				auto out = std::inserter(aliases, std::end(aliases));
				std::transform(std::next(std::cbegin(names)), std::cend(names), out,
							   [&](auto name) { return std::make_pair(std::string{name}, tag); });
			}
		} else {
			std::string file_path{std::format("{}.chrom.sizes", refg_name)};
			try {
				file_path = resolve_datafile_path(prepend_dir(data_dir, file_path));
			} catch (const value_error& e) {
				print("{}\n", e.what());
			}
			for (line_reader lr{file_path}; !lr.done(); ++lr) {
				auto name = lr.line();
				name      = name.substr(0, name.find_first_of(" \t"));
				const chrom_t tag{fnv1a_hash32(name)};
				const auto [it, inserted] = mapping.insert({tag, std::string(name)});
				GK_CHECK(inserted, runtime, "hash collision: '{}' on '{}'", it->second, refg_name);
			}
		}
		chrom_names = chrom_names_t{refg, std::move(refg_name), std::move(mapping), std::move(aliases)};
	}
	return *chrom_names;
}

END_NAMESPACE_GK
