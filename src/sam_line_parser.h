/*
Copyright (C) 2016-2023 Deep Genomics Inc. All Rights Reserved.
*/
#pragma once
#ifndef GK_SAM_LINE_PARSER_H
#define GK_SAM_LINE_PARSER_H

#include "defines.h"
#include "file.h"
#include "genome_anno.h"
#include "interval.h"
#include <iterator>
#include <optional>

BEGIN_NAMESPACE_GK

class sam_line_parser {
public:
	sam_line_parser(const genome_t& genome);

	interval_filter& get_interval_filter() { return _interval_filter; }

	void detect_strand_with_library(const char* format);
	void set_include_duplicates(bool include) { _include_duplicates = include; }

	bool is_stranded() const { return _introns != nullptr || _first_read_direction != direction::unknown; }

protected:
	static constexpr int num_cols{ 12 };
	using sam_cols = std::string_view[num_cols];
	struct cigar_op {
		int length;
		char code;
	};

	void process_file(line_reader& lr);
	virtual void process_line(chrom_t chrom, pos_t pos, strand_t strand, strand_t segment_strand, const cigar_op* codes,
							  int m, int n, sam_cols cols)
		= 0;

	refg_t refg() const { return _chrom_names.refg(); }

private:
	bool infer_strand(strand_t& out, strand_t segment_strand, unsigned int flags) const;
	strand_t infer_strand(strand_t strand, bool strand_set, chrom_t chrom, pos_t pos, const cigar_op* codes,
						  int m, int n) const;

	chrom_names_t   _chrom_names;
	interval_filter _interval_filter;

	// if available, used to determine junction direction from reads
	const intr_table* _introns{};
	enum class direction { unknown, forward, reverse } _first_read_direction{};

	bool _include_duplicates{};
};

END_NAMESPACE_GK

#endif
