/*
Copyright (C) 2016-2023 Deep Genomics Inc. All Rights Reserved.
*/
#include "file.h"
#include "gk_assert.h"
#include "strutil.h"
#include <algorithm>
#include <cerrno>
#include <cstdio>
#include <memory>
#include <type_traits>
#include <utility>
#include <zlib.h>
#ifdef _WIN32
#define WIN32_LEAN_AND_MEAN
#define NOMINMAX
#include <windows.h>

#include <psapi.h>
#define __ftell64 _ftelli64
#define __fseek64 _fseeki64
#else
#include <fcntl.h>
#include <sys/mman.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <unistd.h>
// On POSIX, get 64-bit versions of ftell and fseek by setting _FILE_OFFSET_BITS=64
#if !defined(_FILE_OFFSET_BITS) || (_FILE_OFFSET_BITS != 64)
#error Must define _FILE_OFFSET_BITS=64 project-wide
#endif
#define __ftell64 ftell
#define __fseek64 fseek
#endif

BEGIN_NAMESPACE_GK

const char* default_data_directory = "{GENOMEKIT_DATA_DIR}";

string default_resolve_datafile_path(string path)
{
	// Default implementation just prepends $(GENOMEKIT_DATA_DIR) to any relative path,
	// unless that relative path is to the "tests/data" directory for unit tests
	size_t i = path.find(default_data_directory);
	if (i != string::npos) {
		const char* data_dir = getenv("GENOMEKIT_DATA_DIR");
		GK_CHECK(data_dir, runtime, "Must set GENOMEKIT_DATA_DIR environment variable.");  // destructor WILL get called; no worries
		path.replace(i, strlen(default_data_directory), data_dir);
	}
	return path;
}

string (*resolve_datafile_path)(string path) = default_resolve_datafile_path;

///////////////////////////////////////////////////////////////////////////////

bool is_file(const string& path)
{
	bool found_file;
#ifdef _WIN32
	{
		DWORD attrib = GetFileAttributesA(path.c_str());
		found_file = (attrib != INVALID_FILE_ATTRIBUTES && (attrib & FILE_ATTRIBUTE_DIRECTORY) == 0);
	}
#else
	{
		struct stat st;
		found_file = stat(path.c_str(), &st) == 0 && !S_ISDIR(st.st_mode);
	}
#endif
	return found_file;
}

string prepend_dir(std::string_view dir, std::string_view filename)
{
	if (dir.empty())
		return string{filename};
	string ret{dir};
	if (!endswith(dir, "/") && !endswith(dir, "\\"))
		ret += '/';
	ret += filename;
	return ret;
}

///////////////////////////////////////////////////////////////////////////////

void binary_file::open(const string& path, const char* mode)
{
	GK_CHECK(!is_open(), runtime, "Cannot open new file without closing old one.");
	string binmode(mode);
	if (binmode.find('b') == string::npos)
		binmode.push_back('b');
	_fh = { fopen(path.c_str(), binmode.c_str()), &fclose };
	GK_CHECK(_fh, file, "Could not open {} ({}).", path, strerror(errno));
}

void binary_file::set_seek(size_t offset)
{
	int result = __fseek64(_fh.get(), (long long)offset, SEEK_SET);
	GK_CHECK(result == 0, file, "Error seeking to position {} in file ({}).", offset, strerror(errno));
}

long long binary_file::tell() const
{
	return (long long)__ftell64(_fh.get());
}

int binary_file::write_until_align(int n)
{
	const int residual = (int)(tell() % (long long)n);
	if (!residual)
		return 0;
	static constexpr char zero = 0;
	for (int i = 0; i < n - residual; ++i)
		write(zero);
	return n - residual;
}

void binary_file::read(void* dst, size_t item_size, size_t num_items)
{
	if (num_items > 0) {
		size_t n = fread(dst, item_size, num_items, _fh.get());
		GK_CHECK(n == num_items, file, "Expected to read {} bytes, but read {} bytes ({})",
				 item_size * num_items, item_size * n, strerror(errno));
	}
}

void binary_file::write(const void* src, size_t item_size, size_t num_items)
{
	if (num_items > 0) {
		size_t n = fwrite(src, item_size, num_items, _fh.get());
		GK_CHECK(n == num_items, file, "Expected to write {} bytes, but failed ({})",
				 item_size * num_items, strerror(errno));
	}
}

void binary_file::write_checkpoint(unsigned magic) { write(magic); }

///////////////////////////////////////////////////////////

const char* stdin_path = "<stdin>";

void line_reader::open(const char* path)
{
	_fh = strcmp(path, stdin_path) != 0 ? decltype(_fh){ std::fopen(path, "rb"), &std::fclose }
										: decltype(_fh){ stdin, [](auto) { return 0; } };
	GK_CHECK(_fh, file, "Could not open {} for reading ({}).", path, strerror(errno));
	advance();
}

void line_reader::advance()
{
	// Start at the first byte after the previous line's NULL terminator.
	_line = ++_end;
	_end  = find_delim(_end, _bufend, '\n'); // During constructor, _end==_bufend=0
	if (_end == _bufend) {
		// CASE 1: initial read during constructor, and all subsequent reads
		refill_and_advance();
	} else {
		// CASE 2: all within-buffer newlines
		*_end = '\0';
		if (*_line == '\r')
			++_line;
		if (_line + 1 < _end && _end[-1] == '\r')
			_end[-1] = '\0';
	}
	_line_num++;
}

void line_reader::resize()
{
	//const size_t min_bufsize = (1 << 13)-1; // 8KB holds multiple lines clear into "I/O is the bottleneck" territory.
	const size_t min_bufsize = (1 << 17)-1; // TODO: fix bug when line longer than buffer
	const size_t max_bufsize = (1 << 23)-1; // 8MB lines is probably enough! Something's probably wrong if we blow this.

	// READ CASE 1 (RESIZE BUFFER):
	// _line == _buf means we started scanning from the start _buf
	// and STILL hit _bufend before a newline. So, grow the buffer.
	// Always allocate +1 byte so we can put EOF signal in READ CASE 2B.
	// SPECIAL CASE: first call to advance() from constructor hits this case.
	size_t old_size = _end - _line;
	size_t new_size = old_size*2 >= min_bufsize ? old_size*2 : min_bufsize;
	GK_CHECK(new_size <= max_bufsize, value, "Extremely long line encountered. Something wrong?");
	auto new_buf = std::make_unique<char[]>(new_size+1);

	// Copy the fractional line that we haven't returned yet to the
	// END of the new buffer, so that we'll fall through to READ CASE 2
	// so it will handle the read(), but with a larger buffer than we
	// originally started with.
	//
	//     contents of _buf
	//    [ABCDEFGHIJKLMNOP]
	//     ^               ^
	//     _buf,           _end,
	//     _line           _bufend
	//
	//                                   copy here
	//     contents of new_buf        vvvvvvvvvvvvvvvv
	//    [???????????????????????????ABCDEFGHIJKLMNOP]
	//     ^                          ^               ^
	//     new_buf                   (new_buf        (new_buf
	//                               +new_size       +new_size)
	//                               -old_size)
	if (_buf) {
		memmove(&new_buf[new_size-old_size], _line, old_size);
	}
	_buf = std::move(new_buf);
	_bufend = &_buf[new_size];

	// Insert non-NULL value at trailing buffer position.
	// May be overwritten with NULL if EOF lines up exactly with end of the buffer.
	*_bufend = '\1';

	// Map the _line and _end pointers to the end of the new buffer
	//
	//    [???????????????????????????ABCDEFGHIJKLMNOP]
	//     ^                          ^               ^
	//     _buf                       _line           _end,
	//                                                _bufend
	_end  = _bufend;
	_line = _bufend - old_size;
}

void line_reader::refill()
{
	// READ CASE 2 (FILL BUFFER):
	// _line > _buf means we've already returned at least one line
	// from this buffer, so a resize is not necessary.
	// So, move the fractional line to the start of _buf and read new
	// content into the remainder of the buffer.
	//
	//     contents of _buf (BEFORE, A = fractional line)
	//    [???????????????????????????AAAAAAAAAAAAAAAA]
	//     ^                          ^               ^
	//     _buf                       _line           _end,
	//                                                _bufend
	//     contents of _buf (AFTER MEMMOVE)
	//    [AAAAAAAAAAAAAAAA???????????????????????????]
	//     ^               ^                          ^
	//     _buf,           _end                       _bufend
	//     _line
	//
	size_t len = _end - _line;
	if (_end-_line > 0)
		std::memmove(&_buf[0], _line, len);
	_end  = &_buf[len];
	_line = &_buf[0];

	// Read the remainder of the buffer, if we can.
	//
	size_t nread = fread(_end, (unsigned)(_bufend - _end));
	if (nread == (size_t)(_bufend-_end)) {
		// READ CASE 2A (BUFFER IS FULL):
		// Our buffer is now full, and the new line starts at _buf.
		// case deliberately left empty
		//
		//     contents of _buf (AFTER READ, B = bytes read from disk)
		//    [AAAAAAAAAAAAAAAABBBBBBBBBBBBBBBBBBBBBBBBBBB]
		//     ^               ^                          ^
		//     _buf,           _end                       _bufend
		//     _line
	} else {
		// READ CASE 2B (BUFFER NOT FULL / EOF):
		// Our buffer is NOT full. We've hit EOF. Signify EOF with a NULL
		//
		//     contents of _buf (AFTER READ, n = newline, 0 = NULL terminator)
		//    [AAAAAAAAAAAAAAAABBBBBBBBBBBBBBBBBn0????????]
		//     ^               ^                ^
		//     _buf,           _end             _bufend
		//     _line
		//
		// Notice that the extra NULL terminator
		GK_DBASSERT(nread < (size_t)(_bufend-_end));
		_bufend = _end + nread;
		if (_bufend > &_buf[0] && *(_bufend-1) != '\n') {
			// File didn't end with newline, so add one to trigger RETURN CASE 1
			// and a NULL terminator to subsequently trigger
			*(_bufend++) = '\n';
		}
		*_bufend = '\0'; // We can do this because +1 extra bytes always allocated for buffer
	}
}

void line_reader::refill_and_advance()
{
	GK_ASSERT(_end == _bufend);

	// READ CASES occur when we've hit the end of _buf, necessitating read()
	// They fall through and continue the loop until RETURN CASES occur.

	for (;;) {
		if (_end == _bufend) {
			if (_line == &_buf[0])
				resize();
			refill();
		}

		// Try to find the next newline.
		_end = find_delim(_end, _bufend, '\n');
		if (_end < _bufend) {
			// RETURN CASES always occur when _end is within buffer extents.

			if (*_end == '\n') {
				// RETURN CASE 1: (NEWLINE FOUND)
				// We found a newline character inside the buffer. Convert it to
				// a NULL terminator and stop advancing.
				// Skip If _line currently points to a carriage-return that followed
				// the previous newline, automatically skip it (TEXT mode)
				*_end = '\0';
				if (*_line == '\r')
					++_line;
				if (_line + 1 < _end && _end[-1] == '\r')
					_end[-1] = '\0';
				return;
			}
		} else if (*_end == '\0') {
			// RETURN CASE 2: (EOF)
			//   We hit a NULL value inserted into the buffer by READ CASE 1B. Signal DONE.
			_line = _end = nullptr;
			return;
		}
	}
}

size_t line_reader::fread(char* dst, unsigned bytes)
{
	return std::fread(dst, 1, bytes, _fh.get()); // Call stdlib fread
}

///////////////////////////////////////////////////////////

void zline_reader::open(const char* path)
{
	if (endswith(path, ".gz") || endswith(path, "zip")) {
		// Open compressed file
		_zfh = { gzopen(path, "rb"), &gzclose };
		GK_CHECK(_zfh, file, "Could not open {} for reading ({}).", path, strerror(errno));
		advance();  // OK to call virtual because VMT for zline_reader will be loaded by now
	} else {
		// Open uncompressed file
		line_reader::open(path);
	}
}

size_t zline_reader::fread(char* dst, unsigned bytes)
{
	if (_zfh) {
		int result = gzread(_zfh.get(), dst, bytes);
		GK_CHECK(result >= 0, file, "I/O error reading compressed file ({}).", strerror(errno));
		return (size_t)result;
	}
	return std::fread(dst, 1, bytes, _fh.get()); // Call stdlib fread
}

///////////////////////////////////////////////////////////

void mmap_file::mmap_deleter::operator()(const void* ptr) const noexcept
{
#ifdef _WIN32
	if (ptr != nullptr)
		UnmapViewOfFile(ptr);
#else
	if (ptr != nullptr && ptr != MAP_FAILED)
		munmap(ccast<void*>(ptr), size);
#endif
}

void mmap_file::open(const string& path)
{
#ifdef GK_DEBUG_FILES
	print("mmap_file({})::open(\"{}\")\n", this, path);
#endif

#ifdef _WIN32
	using unique_handle = std::unique_ptr<std::remove_pointer_t<HANDLE>, decltype(&CloseHandle)>;
	auto fh             = unique_handle(
        CreateFileA(path.c_str(), GENERIC_READ, FILE_SHARE_READ, nullptr, OPEN_EXISTING, FILE_ATTRIBUTE_NORMAL, nullptr),
        &CloseHandle);
	GK_CHECK(fh.get() != INVALID_HANDLE_VALUE, file, "Could not open {} for reading ({}).", path, strerror(errno));

	LARGE_INTEGER file_size{};
	if (!GetFileSizeEx(fh.get(), &file_size)) {
		GK_THROW(file, "Could not determine file size for {}", path);
	}

	auto mh = unique_handle(CreateFileMappingA(fh.get(), nullptr, PAGE_READONLY, file_size.HighPart, file_size.LowPart, nullptr),
							&CloseHandle);
	if (!mh) {
		GK_THROW(file, "Could not create file mapping for {}", path);
	}

	_data.reset(MapViewOfFile(mh.get(), FILE_MAP_READ, 0, 0, file_size.QuadPart));
	if (!_data) {
		GK_THROW(file, "Could not map view of file {}", path);
	}
	_size = file_size.QuadPart;

#else
	struct unique_fd {
		NOCOPY(unique_fd)
		unique_fd(int fd) noexcept
			: fd(fd)
		{
		}
		~unique_fd() noexcept { ::close(fd); }
		int fd;
	};
	off_t file_size = 0;
	unique_fd fh { ::open(path.c_str(), O_RDONLY) };
	GK_CHECK(fh.fd >= 0, file, "Could not open {} for reading ({})", path, strerror(errno));

	file_size = lseek(fh.fd, 0, SEEK_END);
	if (file_size < 0) {
		GK_THROW(file, "Could not determined file size for {} ({})", path, strerror(errno));
	}

	_data = { mmap(nullptr, (size_t)file_size, PROT_READ, MAP_SHARED, fh.fd, 0), mmap_deleter{ (size_t)file_size } };
	if (_data.get() == MAP_FAILED) {
		GK_THROW(file, "Could not map view of file ({})", strerror(errno));
	}

	// Closing POSIX file does not cause it to be unmapped.
	_size = (size_t)file_size;
#endif
}

void mmap_file::close()
{
	_data.reset();
	_size = 0;
}

void mmap_file::read_checkpoint(unsigned magic)
{
	auto actual = read<unsigned>();
	GK_CHECK(magic == actual, file, "File I/O checkpoint expected to be '{:x}' but found '{:x}'. Binary file format does not match expectations.", magic, actual);
}

END_NAMESPACE_GK
