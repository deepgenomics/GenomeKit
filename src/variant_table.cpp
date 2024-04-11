/*
Copyright (C) 2016-2023 Deep Genomics Inc. All Rights Reserved.
*/
#include "variant_table.h"

#include "genome.h"
#include "genome_dna.h"
#include "strutil.h"

#include <cctype>
#include <iterator>
#include <numeric>
#include <unordered_map>
#include <unordered_set>

using namespace std;
using namespace std::literals;

namespace {
template <class T, class P>
bool contains_if(const T& container, P pred)
{
	return find_if(cbegin(container), cend(container), pred) != cend(container);
}
}

BEGIN_NAMESPACE_GK

const unsigned short c_vcfbin_sig = 0xc0ed;
const unsigned short c_vcfbin_ver = 0x0005;
// versions:
//   0001: initial format
//   0002: remove extra variants_size field
//   0003: INFO fields
//   0004: bool dtype/sample_names
//   0005: Restructured to support all chromosomes and more species
//   <--- INSERT VERSION CHANGE SUMMARIES HERE

enum vcf_col_t {
	vcf_col_chrom,  // 0
	vcf_col_pos,    // 1
	vcf_col_id,     // 2
	vcf_col_ref,    // 3
	vcf_col_alts,   // 4
	vcf_col_qual,   // 5
	vcf_col_filter, // 6
	vcf_col_info,   // 7
	vcf_col_format  // 8
};

// For diploids, assume GT of the form 00, 01, 10, 11 only
// 		Convert gts strings to a 2-bit encoding
// 		00 -> 0     (gt_homozygous_ref)
// 		01 -> 1     (gt_heterozygous)
// 		11 -> 2     (gt_homozygous_alt)
// 		?? -> 3     (gt_unknown)
static vcf_table::gt_t as_gt_t(string_view gts, bool strict_gt = false)
{
	static const int diploid_gts_length = 3;
	static const int haploid_gts_length = 1;

	int gts_length = (int)size(gts);
	if (gts_length == diploid_gts_length) {
		char first  = gts[0];
		char delim  = gts[1];
		char second = gts[2];

		GK_CHECK((isdigit(first) || first == '.') && (isdigit(second) || second == '.'), value,
				 "GT component must an allelic index or '.': {}.", gts);
		GK_CHECK(delim == '/' || delim == '|', value, "Expected `/` or `|` as separator for GT field, Found {}", delim);

		if (strict_gt && (first == '.' || second == '.'))
			return vcf_table::gt_unknown;
		if (first != second)
			return vcf_table::gt_heterozygous;
		if (first == '0')
			return vcf_table::gt_homozygous_ref;
		if (first == '.')
			return vcf_table::gt_unknown;
		return vcf_table::gt_homozygous_alt;
	}

	GK_CHECK(gts_length == haploid_gts_length, value,
			 "Expected GT value to be defined as a haploid/diploid, but found {}\n", gts);
	if (gts[0] == '1')
		return vcf_table::gt_homozygous_alt;
	if (gts[0] == '.')
		return vcf_table::gt_unknown;
	return vcf_table::gt_homozygous_ref;
}

////////////////////////////////////////////////
// variant_t
////////////////////////////////////////////////

variant_t::variant_t(const packed_variant& src, const variant_table& table)
{
	as_interval() = src.as_interval();
	if (src.aux.hi == 255) {
		// Ref and alt were both packed inplace, into lower 4 bytes of the aux field.
		ref = rcast<const char*>(&src.aux.lo) + 0;
		alt = rcast<const char*>(&src.aux.lo) + 2;
	} else {
		// Ref and alt were stored in the aux pool.
		const char* aux = table.aux() + src.aux.as64();
		ref             = aux;
		alt             = strchr(aux, '\0') + 1;
	}
}

////////////////////////////////////////////////
// vcf_table
////////////////////////////////////////////////

void vcf_table::field_col_t::load(mmap_file& in)
{
	// Load dtype, depth
	in.read(dtype);
	in.read(depth);

	// Load ID string
	id = in.curr_ptr<char>();
	in.move_seek(strlen(id) + 1);
	in.read_until_align32();

	// Load data array
	switch (dtype) {
	case bool_: [[fallthrough]];
	case int8: {
		array_view<int8_t> view;
		view.read(in);
		data             = view.data();
		data_file_offset = view.data_file_offset();
		break;
	}
	case int16: {
		array_view<int16_t> view;
		view.read(in);
		data             = view.data();
		data_file_offset = view.data_file_offset();
		break;
	}
	case int32: {
		array_view<int32_t> view;
		view.read(in);
		data             = view.data();
		data_file_offset = view.data_file_offset();
		break;
	}
	case float16: {
		array_view<half_t> view;
		view.read(in);
		data             = view.data();
		data_file_offset = view.data_file_offset();
		break;
	}
	case float32: {
		array_view<float> view;
		view.read(in);
		data             = view.data();
		data_file_offset = view.data_file_offset();
		break;
	}
	case str: {
		in.read_checkpoint(0x31231259);
		in.read_until_align64();
		auto data_size   = in.read<uint64_t>();
		data_file_offset = in.curr_seek();
		data             = in.curr_ptr<char>();
		in.move_seek(data_size);
		in.read_checkpoint(0x31231257);
		break;
	}
	default: GK_UNREACHABLE();
	}
	in.read_until_align32();
}

const vcf_table::field_col_t* vcf_table::field_cols_t::get(const char* id) const noexcept
{
	auto it = lower_bound(cbegin(cols), cend(cols), id, [](const auto& x, auto id) { return strcmp(x.id, id) < 0; });
	if (it == cend(cols) || strcmp(id, it->id) != 0) {
		return nullptr;
	}
	return &*it;
}

void vcf_table::field_cols_t::load(mmap_file& mapped)
{
	cols.resize(mapped.read<int>());
	for (auto& x : cols) { x.load(mapped); }
}

template <typename T>
struct vcf_value_decoder {
};

template <>
struct vcf_value_decoder<int8_t> {
	int default_value;
	int8_t operator()(string_view x) { return int_cast<int8_t>(x != "." ? as_int(x) : default_value); }
};

template <>
struct vcf_value_decoder<int16_t> {
	int default_value;
	int16_t operator()(string_view x) { return int_cast<int16_t>(x != "." ? as_int(x) : default_value); }
};

template <>
struct vcf_value_decoder<int32_t> {
	int default_value;
	int32_t operator()(string_view x) { return int_cast<int32_t>(x != "." ? as_int(x) : default_value); }
};

template <>
struct vcf_value_decoder<half_t> {
	float default_value;
	half_t operator()(string_view x)
	{
		if (x != ".") {
			return half_t(as_float(x));
		}
		return half_t(default_value);
	}
};

template <>
struct vcf_value_decoder<float> {
	float default_value;
	float operator()(string_view x) { return x != "." ? as_float(x) : default_value; }
};

template <typename T, typename U>
void vcf_value_decode(vector<T>& out, const vector<string_view>& values, U default_value)
{
	vcf_value_decoder<T> decoder = { default_value };
	// don't reserve since we don't know how many more variants will follow
	transform(values.begin(), values.end(), back_inserter(out), decoder);
}

void vcf_table::builder::field_values::add_values(string_view values, vector<string_view>& split_buf,
												  long long line_num)
{
	using std::empty;

	GK_CHECK(dtype.has_value(), value, "Missing a meta-information header for ID {}.", id);

	if (id == "GT"sv) {
		// GT has a special missing code: ./.
		// VCF standard is not clear on whether a . datablock or the GT definition, ./., takes precedence. Assume GT
		// since . can mean datablock missing or haplotype missing (and we only support diploids).
		GK_ASSERT(*dtype == int8);
		GK_ASSERT(depth.value_or(0) == 1);
		data_i8.push_back(as_gt_t(values, default_int == gt_unknown));
		return;
	} else if (id == "SVTYPE"sv) {
		GK_ASSERT(*dtype == int8);
		GK_ASSERT(depth.value_or(0) == 1);
		constexpr const char* svtype_values[(size_t)svtype_t::count] = {
			".", "DEL", "INS", "DUP", "INV", "CNV", "BND", // same order as enum
		};
		auto sv_begin = cbegin(svtype_values), sv_end = cend(svtype_values),
			 found = find_if(sv_begin, sv_end, [values](auto x) { return x == values; });
		if (found != sv_end) {
			data_i8.push_back(int_cast<int8_t>(distance(sv_begin, found)));
		} else if (auto loc = _unknown_svtype.lower_bound(values); loc == cend(_unknown_svtype) || *loc != values) {
			_unknown_svtype.emplace_hint(loc, values);
			print("Unsupported SVTYPE={} value on line {}.\n", values, line_num);
		}
		return;
	}

	if (values == ".") {
		GK_CHECK(depth.has_value(), value, "Missing a meta-information header for ID {}.", id);

		switch (*dtype) {
		case bool_: data_i8.push_back(0); break;
		case int8: data_i8.resize(data_i8.size() + *depth, int_cast<int8_t>(default_int)); break;
		case int16: data_i16.resize(data_i16.size() + *depth, int_cast<int16_t>(default_int)); break;
		case int32: data_i32.resize(data_i32.size() + *depth, int_cast<int32_t>(default_int)); break;
		case float16: data_f16.resize(data_f16.size() + *depth, half_t(default_float)); break;
		case float32: data_f32.resize(data_f32.size() + *depth, default_float); break;
		case str:
			GK_ASSERT(depth.value_or(0) == 1);
			data_i8.push_back('\0');
			break;
		default: GK_UNREACHABLE();
		}
		return;
	}

	switch (*dtype) {
	case str:
		// leave the string splitting to the client to simplify encoding varying length content
		GK_ASSERT(depth.value_or(0) == 1);
		data_i8.insert(data_i8.end(), cbegin(values), cend(values));
		data_i8.push_back(0);
		return;
	case bool_:
		GK_CHECK(empty(values), value, "Expected 0 values for Flag ID \"{}\"", id);
		data_i8.push_back(1);
		return;
	default: break;
	}

	split_view(values, ',', split_buf);
	int num_values = (int)split_buf.size();
	if (!depth.has_value())
		depth = num_values;

	const char* warn_too_many
		= "This can be caused by incorrect splitting of multi-allelic variants during normalization. "
		  "For example, if a field has one number per ref + alt alleles, yet the VCF header specifies "
		  "arbitrary multiplicity ('Number=.' instead of 'Number=R'), then 'bcftools norm' will fail to split "
		  "the values correctly. If that happens, try correcting the upstream VCF header and normalizing again."
		  "See VCF 4.2 Specification sections 1.4.2 and 1.6.2.";
	GK_CHECK(*depth == num_values, value, "Expected {} values, found {}, for ID \"{}\". {}", *depth,
			 num_values, id, num_values < *depth ? "" : warn_too_many);

	switch (*dtype) {
	case int8: vcf_value_decode(data_i8, split_buf, default_int); break;
	case int16: vcf_value_decode(data_i16, split_buf, default_int); break;
	case int32: vcf_value_decode(data_i32, split_buf, default_int); break;
	case float16: vcf_value_decode(data_f16, split_buf, default_float); break;
	case float32: vcf_value_decode(data_f32, split_buf, default_float); break;
	default: GK_UNREACHABLE();
	}
}
void vcf_table::builder::field_values::dump(binary_file& out)
{
	GK_CHECK(dtype.has_value() && depth.has_value(), value, "No values or meta-information header found for ID {}.", id);

	out.write<int32_t>(*dtype);
	out.write<int32_t>(*depth);
	out.write(id.c_str(), id.size() + 1);
	out.write_until_align(4);
	switch (*dtype) {
	case str: dump_str(out); break;
	case bool_: [[fallthrough]];
	case int8: out.write_array(data_i8); break;
	case int16: out.write_array(data_i16); break;
	case int32: out.write_array(data_i32); break;
	case float16: out.write_array(data_f16); break;
	case float32: out.write_array(data_f32); break;
	default: GK_UNREACHABLE();
	}
	out.write_until_align(4);
}

void vcf_table::builder::field_values::dump_str(binary_file& out)
{
	out.write_checkpoint(0x31231259);
	size_t raw_size = data_i8.size();

	// Count strings, find maximum length string, pool strings, and count total bytes for pooled strings.
	string_set<const char*> pool;
	size_t num_str         = 0;
	size_t max_len         = 0;
	size_t pooled_raw_size = 0;
	const char* a          = raw_size > 0 ? rcast<const char*>(&data_i8[0]) : nullptr;
	const char* b          = a + raw_size;
	const char* s;
	for (s = a; s != b; ++num_str) {
		size_t len = strlen(s);

		if (len > max_len)
			max_len = len;

		if (pool.insert(s).second)
			pooled_raw_size += len + 1;

		s += len + 1;
	}

	// Compare the sizes of storing the strings three different ways:
	//  1) as raw pool of fixed-length strings;
	//  2) as raw pool of variable-length strings, with an array of offsets (one offset per variant);
	//  3) as raw pool of variable-length strings, with an array of offsets (one offset per category) and an array of
	//  category indexes (one index per variant).
	size_t num_pooled = pool.size();
	size_t pool_index_size
		= num_pooled <= 0x000000ff ? 1 : num_pooled <= 0x0000ffff ? 2 : num_pooled <= 0xffffffff ? 4 : 8;

	size_t varlen_offsets_size = (num_str + 1) * sizeof(offset40_t);
	size_t pooled_offsets_size = (num_pooled + 1) * sizeof(uint64_t);

	size_t size_as_fixlen = num_str * (max_len + 1);
	size_t size_as_varlen = varlen_offsets_size + raw_size;
	size_t size_as_pooled = pooled_offsets_size + num_str * pool_index_size + pooled_raw_size;

	// Penalize the pooling method because when used from Python it results in lots of PyString instances
	// hanging around after the user may in fact be finished with them, so heuristically require a big
	// pooling factor to justify that potential for extra memory usage.
	const size_t pool_overhead_factor = 4;

	if (size_as_fixlen <= min(size_as_varlen, size_as_pooled * pool_overhead_factor)) {
		// Mode 1: fixed-length
		vector<char> fixed_data;
		fixed_data.reserve(size_as_fixlen);
		for (s = a; s != b;) {
			size_t len = strlen(s);
			fixed_data.insert(fixed_data.end(), s, s + len);              // Fill string bytes
			fixed_data.insert(fixed_data.end(), max_len - len + 1, '\0'); // Fill with NULL
			s += len + 1;
		}

		out.write_until_align(8);
		out.write<uint64_t>(sizeof(uint32_t) * 2 + size_as_fixlen);
		out.write<uint32_t>(fixlen);
		out.write<uint32_t>((uint32_t)max_len + 1); // write stride, which includes NULL
		if (!fixed_data.empty())
			out.write(&fixed_data[0], fixed_data.size());

	} else if (size_as_varlen <= size_as_pooled * pool_overhead_factor) {
		// Mode 2: variable-length

		vector<offset40_t> offsets;
		offsets.reserve(num_str + 1);
		offsets.emplace_back(varlen_offsets_size);
		for (s = a; s != b;) {
			size_t len = strlen(s);
			offsets.emplace_back(offsets.back().as64() + len + 1);
			s += len + 1;
		}

		out.write_until_align(8);
		out.write<uint64_t>(sizeof(uint32_t) + size_as_varlen);
		out.write<uint32_t>(varlen);
		if (!offsets.empty()) {
			out.write(&offsets[0], offsets.size());
			out.write(&data_i8[0], data_i8.size());
		}
	} else {
		// Mode 3: variable-length pooled

		// Assign a unique index to each pooled string, and an offset, and compile the pooled strings into single raw
		// char array.
		unordered_map<const char*, uint64_t> indexes;
		vector<char> pooled_data;
		vector<uint64_t> offsets;
		offsets.reserve(num_pooled + 1);
		offsets.push_back(
			pooled_offsets_size
			+ num_str * pool_index_size); // raw string data follows offsets array and per-variant indices array
		for (auto k : pool) {
			auto index = (uint32_t)indexes.size();
			indexes[k] = index;
			offsets.push_back(offsets.back() + strlen(k) + 1);
			pooled_data.insert(pooled_data.end(), k, k + strlen(k) + 1);
		}

		// Build the per-variant indices into the pool offset array
		vector<uint8_t> pooled_index_1;
		vector<uint16_t> pooled_index_2;
		vector<uint32_t> pooled_index_4;
		vector<uint64_t> pooled_index_8;
		for (s = a; s != b; s += strlen(s) + 1) {
			// Find the char ptr key used by indexes, then use that to look up the index for that string
			auto k = pool.find(s);
			switch (pool_index_size) {
			case 1: pooled_index_1.push_back((uint8_t)indexes[*k]); break;
			case 2: pooled_index_2.push_back((uint16_t)indexes[*k]); break;
			case 4: pooled_index_4.push_back((uint32_t)indexes[*k]); break;
			case 8: pooled_index_8.push_back((uint64_t)indexes[*k]); break;
			}
		}

		// Dump the layout to disk
		out.write_until_align(8);
		out.write<uint64_t>(sizeof(uint32_t) * 2 + sizeof(uint64_t) * 1 + size_as_pooled);
		out.write<uint32_t>(pooled);
		out.write<uint32_t>((uint32_t)pool_index_size);
		out.write<uint64_t>((uint64_t)num_pooled);
		if (!offsets.empty()) {
			out.write(&offsets[0], offsets.size());
			switch (pool_index_size) {
			case 1: out.write(&pooled_index_1[0], pooled_index_1.size()); break;
			case 2: out.write(&pooled_index_2[0], pooled_index_2.size()); break;
			case 4: out.write(&pooled_index_4[0], pooled_index_4.size()); break;
			case 8: out.write(&pooled_index_8[0], pooled_index_8.size()); break;
			}
			out.write(&pooled_data[0], pooled_data.size());
		}
	}
	out.write_checkpoint(0x31231257);
}

vcf_table::builder::builder(const char* infile, const genome_t& genome, bool validate)
	: _infile(infile)
	, _infile_name(infile)
	, _dna{&genome.dna()}
	, _chrom_names{genome.chrom_names()}
	, _validate(validate)
{
}

void vcf_table::builder::collect_info(const char* id, optional<dtype_t> dtype, const void* default_value)
{
	GK_CHECK(!contains_if(_info_values, [id](const auto& x) { return x.id == id; }), value,
			 "Already specified INFO ID \"{}\"", id);
	GK_ASSERT(default_value == nullptr || dtype);

	if (strcmp(id, "ID") == 0) {
		_info_values.push_back({ id, str, 1 });
		return;
	}
	if (strcmp(id, "SVTYPE") == 0) {
		int val = (int)svtype_t::na;
		if (default_value) {
			if (*dtype == int32) {
				val = *(const int*)default_value;
			}
			GK_CHECK(
				*dtype == int32 && 0 <= val && val < (int)svtype_t::count, value,
				"INFO ID SVTYPE must be one of SVTYPE_NA, SVTYPE_DEL, SVTYPE_INS, SVTYPE_DUP, SVTYPE_INV, SVTYPE_CNV, "
				"or SVTYPE_BND.");
		}
		_info_values.push_back({ id, int8, 1, val });
		return;
	}

	try {
		collect_field(_info_values, id, dtype, default_value);
	}
	GK_RETHROW("FORMAT ID {}", id);
}

void vcf_table::builder::collect_fmt(const char* id, optional<dtype_t> dtype, const void* default_value)
{
	GK_CHECK(!contains_if(_fmt_values, [id](const auto& x) { return x.id == id; }), value,
			 "Already specified FORMAT ID \"{}\"", id);
	GK_ASSERT(default_value == nullptr || dtype);

	if (strcmp(id, "GT") == 0) {
		auto val = (int)gt_unknown;
		if (default_value) {
			if (*dtype == int32) {
				val = *(const int*)default_value;
			}
			GK_CHECK(*dtype == int32 && (val == (int)gt_unknown || val == (int)gt_heterozygous), value,
					 "FORMAT ID GT must be one of GT_UNKNOWN or GT_HETEROZYGOUS.");
		}
		optional<int> depth;
		if (!contains_if(_fmt_values, [](const auto& x) { return x.id == "PS"; })) {
			// unphased genotypes aggregate to single enum value
			depth = 1;
		}
		_fmt_values.push_back({ id, int8, depth, val });
		return;
	}
	if (strcmp(id, "PS") == 0) {
		using std::begin;
		using std::end;
		if (auto gt = find_if(begin(_fmt_values), end(_fmt_values), [](const auto& x) { return x.id == "GT"; });
			gt != end(_fmt_values)) { // multidimensional genotype for each allele
			gt->depth.reset();
		}
	}

	try {
		collect_field(_fmt_values, id, dtype, default_value);
	}
	GK_RETHROW("FORMAT ID {}", id);
}

void vcf_table::builder::collect_field(field_values_t& fields, const char* id, optional<dtype_t> dtype,
									   const void* default_value)
{
	if (!default_value) {
		fields.push_back({ id, dtype });
		return;
	}

	switch (*dtype) {
	case int8:
	case int16:
	case int32: fields.push_back({ id, dtype, nullopt, *(const int*)default_value }); return;
	case float16:
	case float32: fields.push_back({ id, dtype, nullopt, 0, *(const float*)default_value }); return;
	default: GK_THROW(value, "user default values are supported only for integer or floating point types.");
	}
}

static bool is_r_sized_field(const string& id) { return id == "AD"sv || id == "ADF"sv || id == "ADR"sv; }

void vcf_table::builder::parse_vcf_metainfo_line(string_view line)
{
	using std::begin;
	using std::empty;
	using std::end;
	using std::size;

	bool is_info = true;
	auto skipped = skip_prefix(line, "INFO=<");
	if (size(skipped) == size(line)) {
		skipped = skip_prefix(line, "FORMAT=<");
		if (size(skipped) == size(line)) {
			return;
		}
		is_info = false;
	}
	line = skipped;

	const char* const col_type = is_info ? "INFO" : "FORMAT";
	field_values_t& values     = is_info ? _info_values : _fmt_values;

	enum { keys_id, keys_number, keys_type, keys_description, keys_count };
	string_view keys[keys_count];

	auto get_value = [&](auto kv, string_view id) {
		GK_CHECK(startswith(kv, id) && kv[size(id)] == '=', value, "Expected {} key.", id);
		return kv.substr(size(id) + 1);
	};

	// OK to include trailing '>' since we don't use the Description key
	GK_CHECK(split_view(line, ',', keys, keys_count) >= keys_count, value, "Expected at least {} keys.",
			 as_ordinal(keys_count));

	auto value = get_value(keys[keys_id], "ID");
	auto start = begin(values);
	auto stop  = end(values);
	auto it    = find_if(start, stop, [value](const auto& x) { return x.id == value; });
	if (it == stop)
		return;

	if (!it->depth.has_value()) {
		value = get_value(keys[keys_number], "Number");
		if (!empty(value) && isdigit(value[0])) { // disallow negative numbers
			it->depth = as_int(value);
		} else if (value == "A") {
			it->depth = 1; // assert multiallelic sites normalized into multiple rows
		} else if (value == "R" || is_r_sized_field(it->id)) {
			// Number=R is defined in VCF >=4.2
			it->depth = 2; // assert multiallelic sites normalized into multiple rows
		} else if (value == "G") {
			it->depth = 2; // assert only diploids
		} else {
			// unknown/unbounded depth: defer encode size and set based on first observed row
			GK_CHECK(value == ".", value, "Unsupported Number={} encoding for {} ID {}.",
					 value, col_type, it->id);
		}
	}

	value = get_value(keys[keys_type], "Type");
	if (value == "Integer") {
		if (!it->dtype.has_value()) {
			it->dtype = int32;
		}
		switch (*it->dtype) {
		case int8: [[fallthrough]];
		case int16: [[fallthrough]];
		case int32:
			GK_CHECK(it->depth.value_or(0) > 0, value, "{} ID {} defined as an Integer type must have Number>0.",
					 col_type, it->id);
			break;
		default:
			GK_THROW(value,
					 "{} ID {} defined as an Integer type. "
					 "Default value must be an int, long, numpy.int8, numpy.int16, or numpy.int32.",
					 col_type, it->id);
		}
	} else if (value == "Float") {
		if (!it->dtype.has_value()) {
			it->dtype = float32;
		}
		switch (*it->dtype) {
		case float16: [[fallthrough]];
		case float32:
			GK_CHECK(it->depth.value_or(0) > 0, value,
					 "{} ID {} defined as a Float type must have Number>0.", col_type, it->id);
			break;
		default:
			GK_THROW(value,
					 "{} ID {} is defined as a Float type. "
					 "Default value must be an float, numpy.float16, or numpy.float32.",
					 col_type, it->id);
		}
	} else if (value == "String") {
		if (!(is_info && it->id == "SVTYPE"sv) && !(!is_info && it->id == "GT"sv)) {
			GK_CHECK(!it->dtype.has_value(), value,
					 "{} ID {} is defined as a String type. Default values are unsupported", col_type, it->id);
			GK_CHECK(is_info, value, "String type for FORMAT ID {} are unsupported.", it->id);
			it->dtype = str;
			it->depth = 1; // multi-value strings are unsupported (left to client to split)
		}
	} else if (value == "Flag") {
		GK_CHECK(!it->dtype.has_value(), value, "{} ID {} is defined as a Flag type. Default values are unsupported",
				 col_type, it->id);
		GK_CHECK(it->depth.value_or(1) == 0, value, "{} ID {} defined as a Flag type must have Number=0.", col_type,
				 it->id);
		GK_CHECK(is_info, value, "FLAG types for {} ID {} are unsupported.", col_type, it->id);
		it->dtype = bool_;
		it->depth = 1;
	} else {
		GK_THROW(value, "{} type for {} ID {} are unsupported.", value, col_type, it->id);
	}

	GK_CHECK(it->dtype.has_value(), value, "Variable number of values for {} ID {} are unsupported.",
			 col_type, it->id);
}

void vcf_table::builder::parse_vcf_header_line(string_view line)
{
	if (_fmt_values.empty()) {
		auto num_cols = count(cbegin(line), cend(line), '\t') + 1;
		GK_CHECK(num_cols >= 8, value, "Expected at least 8 tab-separated columns.");
		return;
	}

	split_view(line, '\t', _field_ids_buf);
	GK_CHECK(_field_ids_buf.size() >= 8, value, "Expected at least 8 tab-separated columns.");
	if (_field_ids_buf.size() > 9)
		_sample_names.assign(next(cbegin(_field_ids_buf), 9), cend(_field_ids_buf));
}

bool vcf_table::builder::parse_variant(const vector<string_view>& cols)
{
	using std::empty;
	using std::max;
	using std::size;

	// Parse variant-specific fields, specifically cols[0:7]
	const auto chrom = _chrom_names.as_chrom(cols[vcf_col_chrom]);
	auto ref  = cols[vcf_col_ref];
	auto alts = cols[vcf_col_alts];

	GK_CHECK(alts.find(",") == string_view::npos, value,
			 "Expected vcf file having multiallelic sites split into multiple rows. "
			 "Try setting normalize=True or manually convert with `bcftools norm -m - -f`.");
	GK_CHECK(!ref.empty() && !alts.empty(), value,
			 "Expected vcf file to have normalized variant representations, but found ref/alt of {}/{}. "
			 "Try setting normalize=True or manually convert with `bcftools norm -m - -f`.",
			 ref, alts);
	GK_CHECK(ref.find_first_not_of("actgnACTGN") == string_view::npos, value, "Invalid REF \"{}\".",
			 ref);
	GK_CHECK(ref != alts, value, "Expected '.' for no alternative alleles but found ref/alt of {}/{}.",
			 ref, alts);

	if (alts == ".") {
		// missing alt means no alternative allele (no variant)
		if (!_ancenstral_handler.notify(_infile.line_num()))
			return false;
		alts = ref;
	}

	pos_t start = as_pos(cols[vcf_col_pos]) - 1;
	pos_t end   = invalid_index;

	if (alts.front() == '<' && alts.back() == '>') { // symbolic
		fill_info_field_ids(cols[vcf_col_info]);
		if (auto end_value = get_attr(_field_ids_buf, "END", ""); !empty(end_value)) {
			ref.remove_prefix(1); // Padding base is required for symbolic
			start += 1;
			end = as_int(end_value) - 1; // END = inclusive 1-based
		} else {
			auto svlen = get_attr(_field_ids_buf, "SVLEN", "");
			GK_CHECK(!empty(svlen), value, "Unsupported SV: END/SVLEN INFO ID is required.");
			ref.remove_prefix(1); // Padding base is required for symbolic
			start += 1;
			// TODO: SVLEN in VCF 4.3 is defined as diff in len between REF and ALT (but 1000G has extensions...)
			end = start + int_cast<pos_t>(size(ref)) + as_int(svlen) - 1;
		}
	} else {
		// technically should permit * for overlapping deletions
		auto i_not_acgtn = alts.find_first_not_of("actgnACTGN");
		if (i_not_acgtn == string_view::npos) { // standard
			GK_CHECK(ref.back() != alts.back() || ref == alts, value,
					 "Expected vcf file to have normalized variant representations, but found ref/alt of {}/{}. "
					 "Try setting normalize=True or manually convert with `bcftools norm -m - -f`.",
					 ref, alts);

			// Strip left alignment padding to save storage
			auto [i_mismatch, ignore] = mismatch(cbegin(ref), cend(ref), cbegin(alts), cend(alts));
			auto offset               = int_cast<pos_t>(distance(cbegin(ref), i_mismatch));
			ref.remove_prefix(offset);
			alts.remove_prefix(offset);
			start += offset;
			end = start + int_cast<pos_t>(ref.size()) - 1;
		} else {                                                      // possibly breakend
			if (alts[i_not_acgtn] != '[' && alts[i_not_acgtn] != ']') // not a breakend
				GK_THROW(value, "Invalid ALT \"{}\".", alts);

			if (i_not_acgtn == 0) { // joined before
				auto [i_mismatch, ignore] = mismatch(crbegin(ref), crend(ref), crbegin(alts), crend(alts));
				auto offset               = int_cast<pos_t>(distance(crbegin(ref), i_mismatch));
				ref.remove_suffix(offset);
				alts.remove_suffix(offset);

				end = start + int_cast<pos_t>(ref.size()) - 1;
			} else { // joined after
				GK_CHECK(alts.back() == '[' || alts.back() == ']', value, "Invalid ALT \"{}\".",
						 alts); // not a breakend

				auto [i_mismatch, ignore] = mismatch(cbegin(ref), cend(ref), cbegin(alts), cend(alts));
				auto offset               = int_cast<pos_t>(distance(cbegin(ref), i_mismatch));
				ref.remove_prefix(offset);
				alts.remove_prefix(offset);
				start += offset;
				end = start + int_cast<pos_t>(ref.size()) - 1;
			}
		}
	}

	// Add a fixed-sized entry to the main table.
	packed_variant v;
	v.chrom  = chrom;
	v.pos5   = start;
	v.pos3   = end;
	v.strand = pos_strand;
	v.refg   = _chrom_names.refg();

	if (is_interval_in_list(v, _exclude))
		return false;
	if (!_allow.empty() && !is_interval_in_list(v, _allow))
		return false;

	if (ref.size() < 2 && alts.size() < 2) {
		// Pack ref and alt inplace, into the lower 4 bytes of aux offset. Set high bits to indicate inplace
		// ref/alt.
		v.aux.hi = 255;
		v.aux.lo = 0;
		if (!empty(ref)) {
			memcpy(&v.aux.lo, data(ref), 1);
		}
		if (!empty(alts)) {
			memcpy(rcast<char*>(&v.aux.lo) + 2, data(alts), 1);
		}
	} else {
		// Store variable-length (auxiliary) fields and the aux pool offset for this entry.
		v.aux = _variants.curr_aux();
		_variants.add_aux(ref);
		_variants.add_aux(alts); // can't pool symbolic ALTs here since the only offset is to REF
	}
	_variants.add_elem(v);

	// Validate the newly added variant's REF against the reference genome sequence
	if (_validate) {
		dnastr ref_dna = (*_dna)(v);
		GK_CHECK(ref == string_view(ref_dna.data(), ref_dna.size()), value, "Expected ref \"{}\" to match reference genome \"{}\" at position {}.",
				 ref, ref_dna, v);
	}

	return true;
}

void vcf_table::builder::fill_info_field_ids(string_view info)
{
	if (!_field_ids_buf.empty())
		return;

	split_view(info, ';', _field_ids_buf);
}

void vcf_table::builder::parse_info(const vector<string_view>& cols)
{
	using std::empty;
	using std::size;

	fill_info_field_ids(cols[vcf_col_info]);

	constexpr string_view missing{ "." };
	auto start = cbegin(_field_ids_buf);
	auto stop  = cend(_field_ids_buf);
	for (auto& dst : _info_values) {
		string_view values = missing;
		if (dst.id == "ID"sv) {
			values = cols[vcf_col_id];
		} else {
			auto found = find_if(start, stop, [&dst](auto kv) {
				if (!startswith(kv, dst.id))
					return false;
				return size(kv) == size(dst.id) || kv[dst.id.size()] == '=';
			});
			if (found != stop) {
				values = found->substr(size(dst.id));
				if (!empty(values) && values[0] == '=')
					values.remove_prefix(1);
				else
					GK_CHECK(dst.dtype == bool_, value, "Missing value for INFO ID {}.", dst.id);
			}
		}

		dst.add_values(values, _value_strs_buf, _infile.line_num());
	}
}

void vcf_table::builder::parse_format(const vector<string_view>& cols)
{
	split_view(cols[vcf_col_format], ':', _field_ids_buf);

	// Verifies that the row contains every user-specified FORMAT ID
	// Computes the index to FORMAT string of the value to be extracted from sample data cols
	for (size_t i = 0; i < _fmt_values.size(); ++i) {
		auto it_found       = std::find(_field_ids_buf.begin(), _field_ids_buf.end(), _fmt_values[i].id);
		_field_ids_order[i] = (size_t)std::distance(_field_ids_buf.begin(), it_found);
	}
	size_t field_count = _field_ids_buf.size();
	constexpr string_view missing{ "." };

	// Iterating over samples from col #9 -> \n
	for (size_t i_sample = vcf_col_format + 1; i_sample < cols.size(); ++i_sample) {
		split_view(cols[i_sample], ':', _field_ids_buf);
		GK_CHECK(_field_ids_buf.size() <= field_count, value,
				 "Expects at most {} values (from FORMAT string) but found {} for sample {}.",
				 field_count, _field_ids_buf.size(), i_sample - vcf_col_format);

		for (size_t i_values = 0; i_values < _fmt_values.size(); ++i_values) {
			field_values& dst = _fmt_values[i_values];
			size_t dst_order  = _field_ids_order[i_values];

			// GT is mandatory if specified, don't need to special case ./. for missing constant
			string_view values = dst_order < _field_ids_buf.size() ? _field_ids_buf[dst_order] : missing;
			dst.add_values(values, _value_strs_buf, _infile.line_num());
		}
	}
}

void vcf_table::builder::build(const char* outfile)
{
	using std::begin;
	using std::end;

	GK_ASSERT(_variants.size() == 0); // build can only be called once

	// For fast splitting of VCF lines at the tab-sep level
	vector<string_view> cols;

	_field_ids_order.resize(_fmt_values.size());

	for (; !_infile.done(); ++_infile) {
		auto line = _infile.line();
		if (line.size() < 2)
			continue;

		try {
			if (line[0] == '#') {
				if (line[1] == '#') {
					parse_vcf_metainfo_line(line.substr(2));
				} else {
					// Parse column header row (metadata has two ##)
					parse_vcf_header_line(line.substr(1));
				}
				continue;
			}
			GK_CHECK(_fmt_values.empty() || !_sample_names.empty(), value,
					"No FORMAT fields are specified in the supplied VCF.");

			// Split tab-separated columns into cols
			split_view(line, '\t', cols);
			GK_CHECK(cols.size() >= 8, value, "Expected at least 8 tab-separated columns but found {}",
					cols.size());

			_field_ids_buf.clear(); // mark no INFO splitted
			if (!parse_variant(cols)) {
				continue;
			}

			if (!_info_values.empty())
				parse_info(cols);

			// Parse Format/Data columns if exists
			if (!_sample_names.empty()) {
				int num_sample_cols = (int)cols.size() - (vcf_col_format + 1);
				GK_CHECK((int)_sample_names.size() == num_sample_cols, value,
						"Expects {} samples but found {}", _sample_names.size(), num_sample_cols);
				parse_format(cols);
			}
		}
		GK_RETHROW("In VCF file: {}:{}", _infile_name, _infile.line_num());
	}

	_ancenstral_handler.log();

	// Sort the INFO/FORMAT columns by ID string.
	auto less_id = [](const auto& lhs, const auto& rhs) { return lhs.id < rhs.id; };
	sort(begin(_info_values), end(_info_values), less_id);
	sort(begin(_fmt_values), end(_fmt_values), less_id);

	///////////////////////////////////////////////////////////////
	// Begin writing to file
	///////////////////////////////////////////////////////////////

	// Write the signature
	binary_file out(outfile, "w");
	out.write(c_vcfbin_sig);
	out.write(c_vcfbin_ver);

	// Write the variants table elements + aux data pool
	_variants.dump(out);
	// support easier migration of old vcfbin files
	out.write_until_align(4);

	out.write((int)_info_values.size());
	for (auto& info_value : _info_values) info_value.dump(out);

	// Write FORMAT data as separate columns
	out.write((int)_sample_names.size());
	out.write((int)_fmt_values.size());
	for (auto& fmt_value : _fmt_values) fmt_value.dump(out);

	if (!_sample_names.empty()) {
		auto total_len = int_cast<int>(accumulate(cbegin(_sample_names), cend(_sample_names), (size_t)0,
												  [](auto x, const auto& y) { return x + y.size() + 1; }));
		out.write(total_len);
		for (const auto& name : _sample_names) out.write(name.c_str(), name.size() + 1);
	}

	out.close();
}

bool vcf_table::builder::ancentral_handler::notify(long long line_number)
{
	if (action == action_t::error)
		GK_THROW(value, "Ancestral allele found: remove or build with warn/exclude.");
	if (action == action_t::exclude)
		return false;

	GK_ASSERT(action == action_t::warn);
	_lines.push_back(line_number);
	return true;
}

void vcf_table::builder::ancentral_handler::log() const {
	if (std::empty(_lines))
		return;

	print("Ancestral allele found on lines: ");
	for (auto line : _lines) { print("{}, ", line); }
	print("\n");
}

vcf_table::vcf_table(mmap_file&& mapped)
	: _fmap(std::move(mapped))
{
	auto sig = _fmap.read<unsigned short>();
	auto ver = _fmap.read<unsigned short>();
	GK_CHECK(sig == c_vcfbin_sig, file,
			 "Expected vcfbin file signature {:x} but found {:x}; not a valid vcfbin file?", c_vcfbin_sig,
			 sig);
	// Change this to last backwards compatible version
	GK_CHECK(ver <= c_vcfbin_ver, file, "Expected vcfbin file version at most {:x} but found {:x}.",
			 c_vcfbin_ver, ver);

	// Read and map each table in a specific order
	load(_fmap);
	// support easier migration of old vcfbin files
	_fmap.read_until_align32();

	_info_entries.load(_fmap);

	// Read FORMAT data, sorted by ID string.
	_num_samples = _fmap.read<int>();
	_fmt_entries.load(_fmap);

	if (_num_samples > 0) {
		_fmap.move_seek(sizeof(int)); // num bytes following
		_sample_names = _fmap.curr_ptr<char>();
	}
}

int vcf_table::file_version() noexcept { return c_vcfbin_ver; }

END_NAMESPACE_GK
