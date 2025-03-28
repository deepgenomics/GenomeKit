/*
Copyright (C) 2016-2023 Deep Genomics Inc. All Rights Reserved.
*/
#include "genome.h"

#include "strutil.h"
#include <cctype>

BEGIN_NAMESPACE_GK

genome_t::genome_t(const string& config)
{
	set_config(config);
}

const genome_dna&  genome_t::dna()  const { return _dna;  }
const genome_anno& genome_t::anno() const { return _anno; }

const gene_table& genome_t::genes() const { return anno().genes(); }
const tran_table& genome_t::trans() const { return anno().trans(); }
const exon_table& genome_t::exons() const { return anno().exons(); }
const intr_table& genome_t::intrs() const { return anno().intrs(); }
const cds_table&  genome_t::cdss()  const { return anno().cdss();  }
const utr_table&  genome_t::utr5s() const { return anno().utr5s(); }
const utr_table&  genome_t::utr3s() const { return anno().utr3s(); }

void genome_t::set_config(const string& config)
{
	GK_CHECK(!std::empty(config), value, "No config specified.");

	_config   = config;
	_data_dir = default_data_directory;

	const auto& refgs = get_refg_registry(_data_dir);
	_refg             = refgs.as_refg(_config);
	_refg_name        = refgs.refg_as_sv(_refg);

	_dna.set_source(default_dna_sourcefile(_refg_name, _data_dir));
	if (_config != _refg_name) {
		// Tell the annotation object which file to use, if any
		_anno.set_source(default_anno_sourcefile(_config, _data_dir));
	}
}

bool genome_t::operator==(const genome_t& other) const
{
	return _config == other._config && _data_dir == other._data_dir;
}

bool genome_t::operator!=(const genome_t& other) const { return !(*this == other); }

END_NAMESPACE_GK
