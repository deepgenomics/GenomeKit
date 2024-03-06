/*
Copyright (C) 2016-2023 Deep Genomics Inc. All Rights Reserved.
*/
#include "strutil.h"
#include "gk_assert.h"
#include <algorithm>
#include <charconv>
#include <iterator>
#include <string>
#include <system_error>
#include <version>

using namespace std;

BEGIN_NAMESPACE_GK

/////////////////////////////////////////////

void split_view(string_view s, char delim, vector<string_view>& out, int max_cols)
{
	out.clear();
	while (!empty(s)) {
		if ((int)size(out) + 1 < max_cols) {
			auto pos = s.find(delim);
			out.push_back(s.substr(0, pos));
			if (pos != string_view::npos) {
				s.remove_prefix(pos + 1);
			} else{
				break;
			}
		} else {
			out.push_back(s);
			break;
		}
	}
}

int split_view(string_view s, char delim, string_view* out, int max_cols)
{
	int n = 0;
	for (;!empty(s); ++n) {
		if (n + 1 < max_cols) {
			auto pos = s.find(delim);
			out[n] = s.substr(0, pos);
			if (pos != string_view::npos) {
				s.remove_prefix(pos + 1);
			} else {
				break;
			}
		} else {
			out[n] = s;
			break;
		}
	}
	return n + 1;
}

// Slow but whatever, we won't be parsing GFF3 files too often
string_view get_attr(const vector<string_view>& attrs, string_view name, const char* default_value, char separator)
{
	auto start = cbegin(attrs);
	auto stop  = cend(attrs);
	auto found = find_if(start, stop, [name, separator](auto attr) {
        return startswith(attr, name) && attr[size(name)] == separator;
	});
	if (found != stop) {
		return found->substr(size(name) + 1);
	}
	GK_CHECK(default_value, value, "Did not find expected field '{}'", name);
	return default_value;
}

template <class T>
T as_number(string_view s, const char* type_name)
{
	if (s.starts_with("+"))
		s.remove_prefix(1);

	T    val{};
	auto stop      = s.data() + s.size();
	auto [ptr, ec] = from_chars(s.data(), stop, val);
	if (ptr == stop && ec == errc{})
		return val;

	GK_CHECK(ec != errc::result_out_of_range, value, "Overflow detected when parsing \"{}\" as {}.", s,
			 type_name);
	GK_THROW(value, "Failed to parse \"{}\" as {}.", s, type_name);
	return val;
}

int as_int(string_view s) { return as_number<int>(s, "integer"); }
#if __cpp_lib_to_chars // OSX13.0 still doesn't support floating point
float  as_float(string_view s) { return as_number<float>(s, "float"); }
double as_double(string_view s) { return as_number<double>(s, "double"); }
#else
float as_float(string_view str)
{
	char buf[512];
	GK_CHECK(size(buf) > size(str), value, "buffer size too small for '{}'", str);

	auto out = copy(cbegin(str), cend(str), buf);
	*out     = '\0';
	// Check that strtof uses up the entire string (to catch "1.23abc") and that the string wasn't empty.
	char* endptr;
	auto  v = strtof(buf, &endptr);
	GK_CHECK(endptr != buf && *endptr == '\0', value, "Failed to parse '{}' as a float", buf);
	return v;
}
double as_double(string_view str)
{
	char buf[512];
	GK_CHECK(size(buf) > size(str), value, "buffer size too small for '{}'", str);

	auto out = copy(cbegin(str), cend(str), buf);
	*out     = '\0';
	// Check that strtod uses up the entire string (to catch "1.23abc") and that the string wasn't empty.
	char* endptr;
	auto v = strtod(buf, &endptr);
	GK_CHECK(endptr != buf && *endptr == '\0', value, "Failed to parse '{}' as a double", buf);
	return v;
}
#endif

static char g_dna_complement[256];

static void add_dna_complement_pair(char a, char b)
{
	g_dna_complement[(unsigned char)upper(a)] = upper(b);
	g_dna_complement[(unsigned char)lower(a)] = lower(b);
	g_dna_complement[(unsigned char)upper(b)] = upper(a);
	g_dna_complement[(unsigned char)lower(b)] = lower(a);
}

void init_dna_tables()
{
	memset(g_dna_complement, '?', 256);
	add_dna_complement_pair('A', 'T');
	add_dna_complement_pair('B', 'V');
	add_dna_complement_pair('C', 'G');
	add_dna_complement_pair('D', 'H');
	add_dna_complement_pair('K', 'M');
	add_dna_complement_pair('N', 'N');
	add_dna_complement_pair('R', 'Y');
	add_dna_complement_pair('S', 'S');
	add_dna_complement_pair('W', 'W');
	add_dna_complement_pair('X', 'X');
}

void reverse_complement(char* dst, int size)
{
	for (int i = 0; i < size/2; ++i) {
		char a = g_dna_complement[(unsigned char)dst[i]];
		char b = g_dna_complement[(unsigned char)dst[size-1-i]];
		dst[i] = b;
		dst[size-1-i] = a;
	}
	if (size & 1)
		dst[size/2] = g_dna_complement[(unsigned char)dst[size/2]];
}

END_NAMESPACE_GK
