/*
Copyright (C) 2016-2023 Deep Genomics Inc. All Rights Reserved.
*/
#pragma once
#ifndef __GENOME_KIT_GENOME_H__
#define __GENOME_KIT_GENOME_H__

#include "defines.h"
#include "genome_anno.h"
#include "genome_dna.h"
#include "interval.h"
#include <string>

BEGIN_NAMESPACE_GK

using std::string;

class genome_t {
public:
	explicit genome_t(const string& config);

	const genome_dna&  dna()  const;
	const genome_anno& anno() const;

	// Convenient shortcuts to anno.genes, anno.trans, anno.exons, etc...
	const gene_table& genes() const;
	const tran_table& trans() const;
	const exon_table& exons() const;
	const intr_table& intrs() const;
	const cds_table&  cdss()  const;
	const utr_table&  utr5s() const;
	const utr_table&  utr3s() const;

	INLINE refg_t refg() const { return _refg; }
	INLINE const string& refg_name() const { return _refg_name; }
	INLINE const string& config() const { return _config; }
	INLINE const string& data_dir() const { return _data_dir; }
	const chrom_names_t& chrom_names() const { return get_chrom_names(_refg, _data_dir); }

	void set_config(const string& config);

	bool operator==(const genome_t& other) const;
	bool operator!=(const genome_t& other) const;

private:
	genome_dna _dna;
	genome_anno _anno;
	refg_t _refg{};
	string _refg_name;
	string _config;
	string _data_dir;
};

END_NAMESPACE_GK

#endif // __GENOME_KIT_GENOME_H__
