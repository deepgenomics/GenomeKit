/*
Copyright (C) 2016-2023 Deep Genomics Inc. All Rights Reserved.
*/
#pragma once
#ifndef __GENOME_KIT_FILE_H__
#define __GENOME_KIT_FILE_H__

#include "gk_assert.h"
#include <cstdio>
#include <cstring>
#include <functional>
#include <memory>
#include <string_view>
#include <vector>

struct gzFile_s;

BEGIN_NAMESPACE_GK
using std::memcpy;
using std::string;
using std::vector;

// If data_dir was not specified, then all files should be pulled from GENOMEKIT_DATA_DIR
// directory; the _gk_data module will replace this string with whatever that directory is,
// when we later call resolve_datafile_path() to resolve a data file path.
extern const char* default_data_directory;

bool is_file(const string& path);
string prepend_dir(std::string_view dir, std::string_view filename);

// Call this function pointer when a data file is about to be used.
// Will resolve to an absolute path, and (if Python package hooks in)
// will download the file on-demand before returning (see _gk_data.resolve_datafile_path).
extern string (*resolve_datafile_path)(string path);

template <class T>
concept serializable = std::is_trivially_copyable_v<T> && !std::is_pointer_v<T>;

///////////////////////////////////////////////////////////////////////

class binary_file {
public:
	binary_file() = default;
	binary_file(const string& path, const char* mode) { open(path, mode); }

	void open(const string& path, const char* mode);
	void close() { _fh.reset(); }
	void flush() { std::fflush(_fh.get()); }

	INLINE bool is_open() const { return scast<bool>(_fh); }

	void set_seek(size_t offset);
	long long tell() const;

	void read(       void* dst, size_t item_size, size_t num_items);
	void write(const void* src, size_t item_size, size_t num_items);

	/***
	 * Advance until next N-byte aligned position.
	 * Useful for aligning data in memory-mapped files.
	 *
	 * @return the amount of bytes written
	 */
	int write_until_align(int n);

	// Read/write several items
	template <serializable T> INLINE void read(       T* items, size_t num_items) { read( items, sizeof(T), num_items); }
	template <serializable T> INLINE void write(const T* items, size_t num_items) { write(items, sizeof(T), num_items); }

	// Read/write single item
	template <typename T> INLINE void read(       T& item) { read( &item, 1); }
	template <typename T> INLINE void write(const T& item) { write(&item, 1); }

	// Read/write a dynamic number of items.
	// The number of items is read/written.
	// The size of each item (sizeof(T::value_type)) is read/written.
	// The raw data for the items is read/written.
	template <serializable T> void write_array(const vector<T>& items);

	// Read/write a checkpoint. Useful for sanity checking file contents.
	void write_checkpoint(unsigned magic);

private:
	std::unique_ptr<std::FILE, decltype(&std::fclose)> _fh{ nullptr, [](auto) { return 0; } };
	friend class line_reader;
	friend class zline_reader;
};

/////////////////////////////////////////////////////////////////////

extern const char* stdin_path; // Special path that causes line_reader to pull from stdin

// Much faster than std::getline
// Results do not include the newline character itself.
class line_reader {
public:
	explicit line_reader(const string& path) { open(path.c_str()); }
	virtual ~line_reader() = default;

	INLINE bool  done() const { return _line == nullptr; }
	INLINE std::string_view line() const {
		auto end = _end;
		while (end > _line && end[-1] == '\0') { --end; }
		return std::string_view(_line, end - _line);
	}
	INLINE line_reader& operator++() { GK_DBASSERT(!done()); advance(); return *this; }
	INLINE long long line_num() const { return _line_num; }

protected:
	line_reader() = default;
	void open(const char* path);
	void advance();
	void refill_and_advance();
	void resize();
	void refill();
	virtual size_t fread(char* dst, unsigned bytes);

	using unique_file = std::unique_ptr<std::FILE, decltype(&std::fclose)>;

	char* _line{};                // start of current line (NULL-terminated)
	char* _end                    // end of current line (NULL terminator);
		{ rcast<char*>(-1) };     // -1 allows causes advance() to hit READ CASE 1 on first call
	char* _bufend{};              // end of valid file bytes in the buffer;
								  // not vector to shift directly without copying first on resize
	std::unique_ptr<char[]> _buf; // buffer of content currently read from file
	unique_file _fh               // file handle
		{ nullptr, [](auto) { return 0; } };
	long long _line_num{};
};

/////////////////////////////////////////////////////////////////////

class zline_reader: public line_reader {
public:
	explicit zline_reader(const string& path) { open(path.c_str()); }

private:
	void open(const char* path); // not virtual
	size_t fread(char* dst, unsigned bytes) override;

	std::unique_ptr<gzFile_s, int (*)(gzFile_s*)> _zfh // zlib file handle, used if file ends in .gz
		{ nullptr, [](auto) { return 0; } };
};

/////////////////////////////////////////////////////////////////////

// Memory-mapped read-only file handle.
//
// There are a few important advantages of using memory-mapped files
// in read-only scenarios, when compared to using malloc/fread:
//
//   + Only the parts of the file that are used will be paged in.
//     For example, if you're working with a specific set of genes,
//     or only genic intervals, you'll only page in the chunks you
//     need, while having close to in-memory performance.
//
//   + Several views of the same file can be opened simultaneously,
//     either within the same process or across multiple processes,
//     and they will ALL SHARE the same physical pages in memory.
//     That means you can have 10 processes all opening a 5GB file
//     simultaneously and the system will use at most 5GB physical
//     memory rather than 50GB.
//
//   + Since the OS knows that the memory is read-only and mapped to
//     a file that cannot change, then if it needs to swap out
//     a page of a memory-mapped file it can simply discard it,
//     rather than writing it to swap, since it knows where to find
//     the data the next time it's needed.
//
// As such, it's easy to create many instances of lightweight objects
// that use raw data and data structures stored in a file, sharing the
// memory between all instances and having that memory swapped out
// automatically when not being actively used.
//
// The idea behind adding 'read' and 'set_seek' semantics to a memory mapped
// file is that it allows the code to be very similar to a file being
// read in by fread or binary_file, where some values are 'read' from the
// file (i.e. copied from mapped memory to some other address) and some
// values are merely pointed to (e.g. read_array_view)
//
class mmap_file {
public:
	mmap_file() = default;
	explicit mmap_file(const string& path) { open(path); }

	void open(const string& path);
	void close();

	INLINE bool        is_open() const { return scast<bool>(_data); }
	INLINE const void* data()    const { return _data.get(); }
	INLINE size_t      size()    const { return _size; }

	// Returns pointer of type "const T*" at the given byte offset from start of the file
	template <typename T>
	requires serializable<T> || std::is_void_v<T> INLINE const T* as_ptr(size_t offset) const
	{
		const T* ptr = (const T*)((const char*)_data.get() + offset);
		if constexpr (!std::is_void_v<T>) {
			// data alignment doesn't matter unless it straddles a cache:
			// https://lemire.me/blog/2020/03/18/avoiding-cache-line-overlap-by-replacing-one-256-bit-store-with-two-128-bit-stores/
			// or
			// 4K aliasing on intel:
			// https://lemire.me/blog/2018/01/05/can-32-byte-alignment-alleviate-4k-aliasing/
			// GK_DBASSERT(((size_t)ptr % alignof(T)) == 0, "Unaligned value read from memory mapped file");
		}
		return ptr;
	}

	// Returns pointer of type "const T*" at the current seek offset
	template <typename T>
	requires serializable<T> || std::is_void_v<T> INLINE const T* curr_ptr() const { return as_ptr<T>(curr_seek()); }

	INLINE size_t curr_seek() const      { return _seek; }
	INLINE void set_seek(size_t offset)  { _seek =  offset; }
	INLINE void move_seek(size_t offset) { _seek += offset; }

	// Advance until next 4- or 8-byte aligned position
	// Useful for aligning data in memory-mapped files.
	INLINE void read_until_align32() { size_t residual = _seek % 4; if (residual) _seek += 4-residual; }
	INLINE void read_until_align64() { size_t residual = _seek % 8; if (residual) _seek += 8-residual; }

	// Read single item. The value is actually
	template <serializable T> INLINE T    read()      { T item; read(item); return item; }
	template <serializable T> INLINE void read(T& item) { memcpy(&item, curr_ptr<void>(), sizeof(T)); move_seek(sizeof(T)); }
	template <serializable T> INLINE void read(T* items, size_t num_items) { memcpy(items, curr_ptr<void>(), sizeof(T)*num_items); move_seek(sizeof(T)*num_items); }

	// Read a checkpoint. Useful for sanity checking file contents.
	void read_checkpoint(unsigned magic);

private:
	// clang 4.0 has an issue with extending the std::function's lifetime
	struct mmap_deleter {
		size_t size;
		void operator()(const void* ptr) const noexcept;
	};
	std::unique_ptr<const void, mmap_deleter> _data;

	size_t _size{};
	size_t _seek{};
};

///////////////////////////////////////////////////////////////////////////

template <serializable T>
void binary_file::write_array(const vector<T>& items)
{
	size_t n = items.size();
	write((uint64_t)n);
	write((int)sizeof(T));
	if (n > 0)
		write(&items[0], n);
}

///////////////////////////////////////////////////////////////////////////

END_NAMESPACE_GK

#endif // __GENOME_KIT_FILE_H__
