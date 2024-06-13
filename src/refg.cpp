/*
Copyright (C) 2016-2023 Deep Genomics Inc. All Rights Reserved.
*/
#include "refg.h"

#include "file.h"
#include "format.h"
#include "gk_assert.h"

#include <filesystem>
#include <cstdio>

BEGIN_NAMESPACE_GK

refg_registry_t::refg_registry_t(std::string data_dir): _data_dir{std::move(data_dir)} {}

refg_t refg_registry_t::as_refg(std::string_view config) const
{
	auto it = _refg_by_config.find(config);
	if (it != std::end(_refg_by_config)) {
		return it->second;
	}

	std::string name;
	try {
		// try a small file to see if it's a refg
		// this is done instead of populating a dummy cfg to reduce the
		// effort required to add an assembly
		std::string path{fmt::format("{}.chrom.sizes", config)};
		if (!std::filesystem::exists(path)) {
			path = resolve_datafile_path(prepend_dir(data_dir(), path));
		}
		if (std::filesystem::exists(path)) {
			name = config;
		}
	} catch (const value_error& e) {
	}

	if (std::empty(name)) {
		// not an assembly, must be an annotation, get name from .cfg
		std::string config_path{fmt::format("{}.cfg", config)};
		try {
			config_path = resolve_datafile_path(prepend_dir(data_dir(), config_path));
		} catch (const value_error& e) {
		}
		try {
			for (line_reader lr{config_path}; !lr.done(); ++lr) {
				string_view k_v[2];
				const auto  count = split_view(lr.line(), '=', k_v, std::size(k_v));
				GK_CHECK(count == std::size(k_v), value, "Invalid line in {}:{}", config_path, lr.line_num());
				if (k_v[0] != "refg")
					continue;
				name = strip(k_v[1]);
				break;
			}
		}
		GK_RETHROW("Assembly is missing its chrom.sizes or annotation is missing its cfg file: {}", config);
	}

	const refg_t ref{fnv1a_hash64(name)};

	const auto [config_it, config_inserted] = _refg_by_config.try_emplace(std::string{config}, ref);
	GK_CHECK(config_inserted, runtime, "hash collision, try renaming one of the annotations: '{}' and '{} on '{}'",
			 config_it->first, config, name);
	const auto [name_it, name_inserted] = _names_by_refg.try_emplace(ref, name);
	GK_CHECK(name_inserted || name_it->second == name, runtime,
			 "hash collision, try renaming one of the assemblies: '{}' and '{}'", name_it->first, config);

	return ref;
}

std::string_view refg_registry_t::_try_refg_as_sv_from_file(refg_t ref) const
{
	std::string path{fmt::format("{}.hash", ref)};
	if (!std::filesystem::exists(path)) {
		path = resolve_datafile_path(prepend_dir(data_dir(), path));
	}
	if (!std::filesystem::exists(path)) {
		return {};
	}

	line_reader lr{path};
	auto refg_name = std::string(lr.line());
	refg_name.erase(refg_name.find_last_not_of(" \n\r\t") + 1);
	auto expected_ref = fnv1a_hash64(refg_name);
	GK_CHECK(expected_ref == (uint64_t)ref, runtime, "Hash mismatch in '{}' for '{}': {} != {}",
			 path, refg_name, expected_ref, ref);

	const auto [name_it, name_inserted] = _names_by_refg.try_emplace(ref, refg_name);
	GK_CHECK(name_inserted || name_it->second == refg_name, runtime,
			 "hash collision, try renaming one of the assemblies: '{}' and '{}'", name_it->first, refg_name);

	return name_it->second;
}

std::string_view refg_registry_t::refg_as_sv(refg_t ref) const
{
	const auto it = _names_by_refg.find(ref);
	if (it == std::end(_names_by_refg)) {
		const auto& refg_name = _try_refg_as_sv_from_file(ref);
		if (!refg_name.empty()) {
			return refg_name;
		}
	}
	GK_CHECK(it != std::end(_names_by_refg), value, "Could not retrieve name for {}", ref);
	return it->second;
}

const std::string& refg_registry_t::data_dir() const { return _data_dir; }

const refg_registry_t& get_refg_registry(std::string_view data_dir)
{
	// TODO: data_dir injected as a context
	static string_map<std::string, refg_registry_t> registry_by_data_dir;

	const auto it = registry_by_data_dir.find(data_dir);
	if (it != std::end(registry_by_data_dir)) {
		return it->second;
	}

	refg_registry_t registry{std::string{data_dir}};
	auto [it_insert, inserted] = registry_by_data_dir.try_emplace(registry.data_dir(), std::move(registry));
	return it_insert->second;
}

END_NAMESPACE_GK
