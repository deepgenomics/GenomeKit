/*
Copyright (C) 2016-2023 Deep Genomics Inc. All Rights Reserved.
*/
#pragma once
#ifndef __GENOME_KIT_VARIANT_TABLE_H__
#define __GENOME_KIT_VARIANT_TABLE_H__

#include "file.h"
#include "half.h"
#include "interval.h"
#include "table.h"
#include "util.h"
#include <cstring>
#include <optional>
#include <set>

BEGIN_NAMESPACE_GK

class genome_dna;
class genome_t;
class variant_table;

/////////////////////////////////////////////////////////////////
// packed_variant
/////////////////////////////////////////////////////////////////

#pragma pack(push, 1)
struct packed_variant : public interval_t {
	offset40_t aux;
};
#pragma pack(pop)

// packed_variant memory layout: offset name (size+padding)
//  0 pos5           (4)
//  4 pos3           (4)
//  8 refg           (8)
// 16 chrom          (4)
// 20 strand         (1)
// 21 aux            (5)
// 26 total
//
// Each packed_vcf_variant aux pool entry has a dynamic layout:
//    struct {
//        char ref[n];
//        char alt[m];
//    };

class variant_t : public interval_t {
public:
	const char* ref{};
	const char* alt{};

	variant_t(const interval_t& interval, char* ref, char* alt)
		: interval_t(interval)
		, ref(ref)
		, alt(alt)
	{
	}
	variant_t(const packed_variant& src, const variant_table& table);
};

INLINE bool operator==(const variant_t& x, const variant_t& y)
{
	return x.as_interval() == y.as_interval() && !strcmp(x.alt, y.alt) && !strcmp(x.ref, y.ref);
}
INLINE bool operator!=(const variant_t& x, const variant_t& y) { return !(x == y); }

/////////////////////////////////////////////////////////////////
// variant_table
/////////////////////////////////////////////////////////////////

class variant_table : public interval_table<packed_variant> {
};

/////////////////////////////////////////////////////////////////
// vcf_table
/////////////////////////////////////////////////////////////////

class vcf_table : public variant_table {
public:
	enum dtype_t {
		int8,
		int16,
		int32,
		float16,
		float32,
		str,
		bool_,
		num_dtype,
	};

	enum str_mode_t {
		fixlen,
		varlen,
		pooled,
	};

	enum gt_t { gt_homozygous_ref, gt_heterozygous_unphased, gt_homozygous_alt, gt_unknown, gt_heterozygous_phased_0_1, gt_heterozygous_phased_1_0 };

	enum struct svtype_t : uint8_t {
		na,
		del,
		ins,
		dup,
		inv,
		cnv,
		bnd,
		count // should be last
	};

	// INFO/FORMAT columns with offsets resolved to pointers.
	struct field_col_t {
		void load(mmap_file& in);

		const char* id;
		const void* data;
		int depth;
		dtype_t dtype;
		size_t data_file_offset;
	};
	struct field_cols_t {
		const field_col_t* get(const char* id) const noexcept;
		void load(mmap_file& mapped);

		std::vector<field_col_t> cols;
	};

	static int file_version() noexcept;

	vcf_table(mmap_file&& mapped);

	void close() { _fmap.close(); }

	const field_cols_t& info_fields() const noexcept { return _info_entries; }
	const field_cols_t& fmt_fields() const noexcept { return _fmt_entries; }

	int num_samples() const noexcept { return _num_samples; };
	std::vector<std::string_view> sample_names() const;

	class builder {
	public:
		// keep ordering consistent with string values in py_variant_table.cpp
		enum class action_t { error, warn, exclude };

		builder(const char* infile, const genome_t& genome, bool validate = true);

		void collect_info(const char* id, std::optional<dtype_t> dtype, const void* default_value = nullptr);
		void collect_fmt(const char* id, std::optional<dtype_t> dtype, const void* default_value = nullptr);

		interval_filter& get_interval_filter() { return _interval_filter; }

		void ancestral(action_t behaviour) { _ancenstral_handler._action = behaviour; }

		// Process/Parse through file given by infile
		// Accumulate fields inside private member of builder
		void build(const char* outfile);

	private:
		// Accumulate data for a particular INFO/FORMAT column
		struct field_values {
			NOCOPY(field_values)
			field_values(const char* id, std::optional<dtype_t> dtype, std::optional<int> depth = std::nullopt,
						 int default_int = 0, float default_float = 0.0)
				: id(id)
				, dtype(dtype)
				, depth(depth)
				, default_int(default_int)
				, default_float(default_float)
			{
			}
			string id;
			std::optional<dtype_t> dtype;
			std::optional<int> depth;

			// Only one of these will be used, depending on dtype
			int default_int{};
			float default_float{};
			vector<int8_t> data_i8;
			vector<int16_t> data_i16;
			vector<int32_t> data_i32;
			vector<half_t> data_f16;
			vector<float> data_f32;
			std::set<std::string, std::less<>> _unknown_svtype;

			field_values(field_values&&) = default;
			field_values& operator=(field_values&&) = default;

			void add_values(std::string_view values, vector<std::string_view>& split_buf, long long line_num);
			void dump(binary_file& out);
			void dump_str(binary_file& out);
		};
		using field_values_t = vector<field_values>;

		static void collect_field(field_values_t& fields, const char* id, std::optional<dtype_t> dtype,
								  const void* default_value);

		// Determines the INFO/FORMAT field dimensions
		void parse_vcf_metainfo_line(std::string_view line);
		// Determines sample size in vcf files
		void parse_vcf_header_line(std::string_view line);

		bool parse_variant(const vector<std::string_view>& cols);
		void fill_info_field_ids(std::string_view info);
		void parse_info(const vector<std::string_view>& cols);
		void parse_format(const vector<std::string_view>& cols);

		zline_reader _infile;
		std::string _infile_name;
		const genome_dna* _dna;
		chrom_names_t _chrom_names;

		bool _validate{};

		field_values_t _info_values;
		field_values_t _fmt_values;
		interval_filter _interval_filter;
		variant_table::builder _variants{ false };
		vector<string> _sample_names;

		struct ancentral_handler {
			action_t          _action{action_t::error};
			vector<long long> _lines;

			bool notify(long long line_number);
			void log() const;
		} _ancenstral_handler;

		// Holds the order in which FORMAT IDs were specified on a line containing FORMAT data.
		vector<size_t> _field_ids_order;
		vector<std::string_view> _field_ids_buf;
		vector<std::string_view> _value_strs_buf;
	};

private:
	mmap_file _fmap;
	string _sourcefile;

	int _num_samples{};
	field_cols_t _info_entries;
	field_cols_t _fmt_entries;
	const char* _sample_names{};
};

END_NAMESPACE_GK

#endif // __GENOME_KIT_VARIANT_TABLE_H__
