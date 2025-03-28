/*
Copyright (C) 2016-2023 Deep Genomics Inc. All Rights Reserved.
*/
#pragma once
#ifndef __GENOME_KIT_REFG_H__
#define __GENOME_KIT_REFG_H__

#include "file.h"
#include "strutil.h"
#include "util.h"

#include <cstdint>
#include <unordered_map>
#include <format>

BEGIN_NAMESPACE_GK

enum class refg_t : std::uint64_t {};

inline auto format_as(refg_t x) { return as_ordinal(x); }

class refg_registry_t {
	NOCOPY(refg_registry_t)  // new mappings need to be persisted
public:
	explicit refg_registry_t(std::string data_dir = default_data_directory);
	refg_registry_t(refg_registry_t&&) = default;

	refg_t           as_refg(std::string_view config) const;
	std::string_view refg_as_sv(refg_t refg) const;
	std::string _try_refg_as_sv_from_file(refg_t ref) const;

	const std::string& data_dir() const;

private:
	mutable std::unordered_map<refg_t, std::string> _names_by_refg;
	mutable string_map<std::string, refg_t>         _refg_by_config;
	std::string                                     _data_dir;

	friend const refg_registry_t& get_refg_registry(std::string_view);
};

const refg_registry_t& get_refg_registry(std::string_view data_dir = default_data_directory);

END_NAMESPACE_GK


template <>
struct std::formatter<gk::refg_t> : std::formatter<std::uint64_t> {
	template <typename FormatCtx>
	auto format(gk::refg_t x, FormatCtx& ctx) const
	{
		return std::formatter<std::uint64_t>::format(static_cast<std::uint64_t>(x), ctx);
	}
};
#endif
