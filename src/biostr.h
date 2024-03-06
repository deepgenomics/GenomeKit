/*
Copyright (C) 2016-2023 Deep Genomics Inc. All Rights Reserved.
*/
#pragma once
#ifndef __GENOME_KIT_BIOSTR_H__
#define __GENOME_KIT_BIOSTR_H__

#include "gk_assert.h"
#include <cstring>
#include <fmt/format.h>
#include <memory>
#include <string_view>
#include <vector>

BEGIN_NAMESPACE_GK
using std::vector;

//
// dnastr
// More compact and efficient than std::string or std::vector<char>.
// - smaller sizeof() footprint
// - optimized for fixed size, not push_back.
// - avoids redundant initialization that std::string or std::vector<char>.
// - mutable via non-const version of c_str() and data(), unlike std::string
// - size() excludes null terminator, unlike std::vector
//
class dnastr {
public:
	explicit dnastr(int size);

	INLINE int size() const { return _size; }

	INLINE       char& operator[](int i)       { GK_DBASSERT(i >= 0 && i < _size); return _data[i]; }
	INLINE const char& operator[](int i) const { GK_DBASSERT(i >= 0 && i < _size); return _data[i]; }

	INLINE       char* c_str()       { return _data.get(); }
	INLINE const char* c_str() const { return _data.get(); }
	INLINE       char* data()        { return _data.get(); }
	INLINE const char* data()  const { return _data.get(); }

	bool operator==(std::string_view s) const;
	bool operator==(const dnastr& s) const;
	bool operator!=(std::string_view s) const;
	bool operator!=(const dnastr& s) const;

private:
	std::unique_ptr<char[]> _data;
	int                     _size = 0;
};

//////////////////////////////////////////////////////////////////

inline dnastr::dnastr(int size)
	: _data{ std::make_unique<char[]>(size + 1) }
	, _size{ size }
{
	_data[size] = '\0';
}

inline bool dnastr::operator==(std::string_view s) const { return c_str() == s; }
inline bool dnastr::operator==(const dnastr& s) const
{
	return (size() == s.size()) && (0 == std::memcmp(data(), s.data(), size()));
}
inline bool dnastr::operator!=(std::string_view s) const { return c_str() != s; }
inline bool dnastr::operator!=(const dnastr& s) const
{
	return (size() != s.size()) || (0 != std::memcmp(data(), s.data(), size()));
}

END_NAMESPACE_GK

template <>
struct fmt::formatter<gk::dnastr> : fmt::formatter<const char*> {
	template <typename FormatCtx>
	auto format(const gk::dnastr& x, FormatCtx& ctx) const
	{
		return fmt::formatter<const char*>::format(x.c_str(), ctx);
	}
};

#endif // __GENOME_KIT_BIOSTR_H__
