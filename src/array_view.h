/*
Copyright (C) 2016-2023 Deep Genomics Inc. All Rights Reserved.
*/
#pragma once
#ifndef __GENOME_KIT_ARRAY_VIEW_H__
#define __GENOME_KIT_ARRAY_VIEW_H__

#include "defines.h"
#include "file.h"
#include "gk_assert.h"
#include "util.h"

BEGIN_NAMESPACE_GK

template <typename T>
class array_view {
public:
	INLINE void set_view(const T* data, size_t size) { _data = data; _size = size; }
	INLINE const T* data() const { return _data; }
	INLINE size_t data_file_offset() const { return _data_file_offset; }

	GK_DECLARE_ARRAY_TYPES(T, size_t)
	GK_DECLARE_ARRAY_CONST_METHODS(_data, _size)
	GK_DECLARE_ARRAY_CONST_ITERATOR(_data)

	// Sets the view to be of an array stored in the given memory mapped file.
	// The array in the file is assumed to have been written in the same format
	// as for binary_file::read_array, except in this case the array data itself is
	// not read into memory, but merely pointed to by the array.
	INLINE void read(mmap_file& in)
	{
		auto n = in.read<uint64_t>();
		int    s = in.read<int>();
		GK_CHECK2(s == (int)sizeof(T), file, "Expected to read array with item size '{}' but found '{}'.", sizeof(T), s);
		set_view(in.curr_ptr<T>(), n);
		_data_file_offset = in.curr_seek();
		in.move_seek(n*sizeof(T));
	}

private:
	const T* _data{};
	size_t   _size{};
	size_t   _data_file_offset{};
};

END_NAMESPACE_GK

#endif // __GENOME_KIT_ARRAY_VIEW_H__
