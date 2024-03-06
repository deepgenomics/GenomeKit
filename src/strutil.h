/*
Copyright (C) 2016-2023 Deep Genomics Inc. All Rights Reserved.
*/
#pragma once
#ifndef __GENOME_KIT_STRUTIL_H__
#define __GENOME_KIT_STRUTIL_H__

#include "defines.h"
#include <algorithm>
#include <cctype>
#include <cstring>
#include <functional>
#include <iterator>
#include <string>
#include <string_view>
#include <unordered_map>
#include <unordered_set>
#include <vector>

BEGIN_NAMESPACE_GK
using std::string;
using std::string_view;
using std::vector;

///////////////////////////////////////////////////

inline bool startswith(string_view s, string_view start)
{
	return s.starts_with(start);
}
inline bool endswith(string_view s, string_view end)
{
	return s.ends_with(end);
}
inline bool contains(string_view s, string_view substr) { return s.find(substr) != string::npos; }
inline string_view strip(string_view s)
{
	const auto start = std::find_if(std::begin(s), std::end(s), [](auto x) { return isspace(x) == 0; });
	s = s.substr(std::distance(std::begin(s), start));  // OSX clang not supporting string_view(It, It) correctly
	const auto stop  = std::find_if(std::rbegin(s), std::rend(s), [](auto x) { return isspace(x) == 0; });
	return s.substr(0, std::distance(std::begin(s), stop.base()));
}

inline string_view skip_prefix(string_view str, string_view prefix)
{
	if (startswith(str, prefix)) {
		str.remove_prefix(size(prefix));
	}
	return str;
}

INLINE char upper(char c) { return (c >= 'a' && c <= 'z') ? c + ('A' - 'a') : c; }
INLINE char lower(char c) { return (c >= 'A' && c <= 'Z') ? c + ('a' - 'A') : c; }

///////////////////////////////////////////////////

INLINE char* find_delim(char* begin, char* end, char delim)
{
	char* p = (char*)::memchr(begin, delim, end-begin);
	return p ? p : end;
}

void split_view(string_view s, char delim, vector<string_view>& out, int max_cols=0x7fffffff);
int  split_view(string_view s, char delim, string_view* out, int max_cols = 0x7fffffff);

string_view get_attr(const vector<string_view>& attrs, string_view name, const char* default_value=nullptr, char separator='=');

// Avoid accidentally comparing cstr ptrs in templated code
INLINE bool eq_str(const string& a, const string& b) { return a == b; }
INLINE bool eq_str(const char* a, const char* b) { return strcmp(a, b) == 0; }
INLINE bool eq_str(const string& a, const char* b) { return a == b; }
INLINE bool eq_str(const char* a, const string& b) { return a == b; }

// FNV-1a hash function for byte sequences
INLINE uint64_t fnv1a_hash64(string_view s) noexcept
{
	// https://en.wikipedia.org/wiki/Fowler%E2%80%93Noll%E2%80%93Vo_hash_function#FNV_hash_parameters
	static constexpr uint64_t prime        = (uint64_t{1} << 40) + (1 << 8) + 0xb3;
	static constexpr uint64_t offset_basis = 0xcbf29ce484222325;

	uint64_t x = offset_basis;
	for (auto c : s) {
		x ^= c;
		x *= prime;
	}
	return x;
}

INLINE uint32_t fnv1a_hash32(string_view s) noexcept
{
	// https://en.wikipedia.org/wiki/Fowler%E2%80%93Noll%E2%80%93Vo_hash_function#FNV_hash_parameters
	static constexpr uint32_t prime        = (1 << 24) + (1 << 8) + 0x93;
	static constexpr uint32_t offset_basis = 0x811c9dc5;

	uint32_t x = offset_basis;
	for (auto c : s) {
		x ^= c;
		x *= prime;
	}
	return x;
}

struct string_hash {
	using hash_type      = std::hash<std::string_view>;
	using is_transparent = void;

	std::size_t operator()(const char* str) const noexcept { return hash_type{}(str); }
	std::size_t operator()(std::string_view str) const noexcept { return hash_type{}(str); }
	std::size_t operator()(const std::string& str) const noexcept { return hash_type{}(str); }
};

template <class Key, class T, class Allocator = std::allocator<std::pair<const Key, T>>>
using string_map = std::unordered_map<Key, T, string_hash, std::equal_to<>, Allocator>;

template <class Key, class Allocator = std::allocator<Key>>
using string_set = std::unordered_set<Key, string_hash, std::equal_to<>, Allocator>;

/////////////////////////////////////////////////////

int as_int(string_view s);
float as_float(string_view s);
double as_double(string_view s);

////////////////////////////////////////////////////

void reverse_complement(char* dst, int size);

END_NAMESPACE_GK

#endif // __GENOME_KIT_STRUTIL_H__
