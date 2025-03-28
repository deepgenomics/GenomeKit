/*
Copyright (C) 2016-2023 Deep Genomics Inc. All Rights Reserved.
*/
#include "genome_anno.h"

#include "format.h"
#include "genome.h"
#include "strutil.h"
#include <cstring>
#include <utility>

BEGIN_NAMESPACE_GK

#define GK_GENOME_ANNO_ELEM_UNPACK_BEGIN(name) \
	name##_t::name##_t(const name##_t::packed_value& src, const name##_t::packed_table& table) \
	{ \
		unpack_from(src, table); \
	} \
	void name##_t::unpack_from(const name##_t::packed_value& src, const name##_t::packed_table& table) \
	{ \
		this->as_interval() = src.as_interval(); \
		const auto& anno    = table.anno;

#define GK_GENOME_ANNO_ELEM_UNPACK_END \
	}

////////////////////////////////////////////////
// gene_t
////////////////////////////////////////////////

GK_GENOME_ANNO_ELEM_UNPACK_BEGIN(gene)
	const char* aux = table.aux() + src.aux;
	level     = src.level;
	type      = src.type;
	num_trans = src.num_trans;
	trans.a   = &anno.trans()[src.tran0];
	trans.b   = trans.a + num_trans;
	id        = aux;                    // aux.id
	name      = strchr(id, 0) + 1;      // aux.name
GK_GENOME_ANNO_ELEM_UNPACK_END

////////////////////////////////////////////////
// tran_t
////////////////////////////////////////////////

GK_GENOME_ANNO_ELEM_UNPACK_BEGIN(tran)
	const char* aux = table.aux() + src.aux;
	level     = src.level;  // A transcripts level seems to be distinct from parent gene; not so for types
	type      = src.type;
	tsl       = src.tsl;
	num_exons = src.num_exons;
	num_cdss  = src.num_cdss;
	num_utr5s = src.num_utr5s;
	num_utr3s = src.num_utr3s;
	gene      = &anno.genes()[src.gene];
	exons.a   = &anno.exons()[src.exon0];
	exons.b   = exons.a + num_exons;
	intrs.a   = src.intr0 == invalid_index ? nullptr : &anno.intrs()[src.intr0];
	intrs.b   = src.intr0 == invalid_index ? 0 : intrs.a + num_exons-1;
	cdss.a    = src.cds0 == invalid_index  ? nullptr : &anno.cdss()[src.cds0];
	cdss.b    = src.cds0  == invalid_index ? 0 : cdss.a + num_cdss;
	utr5s.a   = src.utr50 == invalid_index ? nullptr : &anno.utr5s()[src.utr50];
	utr5s.b   = src.utr50 == invalid_index ? 0 : utr5s.a + num_utr5s;
	utr3s.a   = src.utr30 == invalid_index ? nullptr : &anno.utr3s()[src.utr30];
	utr3s.b   = src.utr30  == invalid_index ? 0 : utr3s.a + num_utr3s;
	id        = aux;                        // aux.id
	ccds_id   = strchr(id, 0) + 1;          // aux.ccds_id
	protein_id= strchr(ccds_id, 0) + 1;     // aux.protein_id
	product   = strchr(protein_id, 0) + 1;  // aux.product
GK_GENOME_ANNO_ELEM_UNPACK_END

////////////////////////////////////////////////
// exon_t
////////////////////////////////////////////////

GK_GENOME_ANNO_ELEM_UNPACK_BEGIN(exon)
	const packed_tran& parent_tran = anno.trans()[src.tran];
	index  = src.index;
	tran   = &parent_tran;
	cds    = src.cds == invalid_index ? nullptr : &anno.cdss()[src.cds];
	utr5   = src.utr5 == invalid_index ? nullptr : &anno.utr5s()[src.utr5];
	utr3   = src.utr3 == invalid_index ? nullptr : &anno.utr3s()[src.utr3];
	level  = parent_tran.level;  // inherit
	type   = parent_tran.type;   // inherit
	tsl    = parent_tran.tsl;    // inherit
	id     = table.aux() + src.id; // exon ID from string pool in aux
GK_GENOME_ANNO_ELEM_UNPACK_END

const packed_intr* get_prev_intr(const packed_exon* exon, const genome_anno& anno)
{
	const packed_tran& tran = anno.trans()[exon->tran];
	if (exon->index == 0)
		return nullptr;
	return &anno.intrs()[tran.intr0 + exon->index - 1];
}

const packed_intr* get_next_intr(const packed_exon* exon, const genome_anno& anno)
{
	const packed_tran& tran = anno.trans()[exon->tran];
	if (exon->index == (index_t)tran.num_exons-1)
		return nullptr;
	return &anno.intrs()[tran.intr0 + exon->index];
}

const packed_intr* get_prev_intr(const packed_exon* exon, const genome_t& genome) { return get_prev_intr(exon, genome.anno()); }
const packed_intr* get_next_intr(const packed_exon* exon, const genome_t& genome) { return get_prev_intr(exon, genome.anno()); }

////////////////////////////////////////////////
// intr_t
////////////////////////////////////////////////

GK_GENOME_ANNO_ELEM_UNPACK_BEGIN(intr)
	const packed_tran& parent_tran = anno.trans()[src.tran];
	index  = src.index;
	tran   = &parent_tran;
	level  = parent_tran.level;  // inherit
	type   = parent_tran.type;   // inherit
	tsl    = parent_tran.tsl;    // inherit
GK_GENOME_ANNO_ELEM_UNPACK_END

const packed_exon* get_prev_exon(const packed_intr* intr, const genome_anno& anno)
{
	return &anno.exons()[anno.trans()[intr->tran].exon0 + intr->index];
}

const packed_exon* get_next_exon(const packed_intr* intr, const genome_anno& anno)
{
	return &anno.exons()[anno.trans()[intr->tran].exon0 + intr->index + 1];
}

const packed_exon* get_prev_exon(const packed_intr* intr, const genome_t& genome) { return get_prev_exon(intr, genome.anno()); }
const packed_exon* get_next_exon(const packed_intr* intr, const genome_t& genome) { return get_prev_exon(intr, genome.anno()); }

////////////////////////////////////////////////
// cds_t
////////////////////////////////////////////////

GK_GENOME_ANNO_ELEM_UNPACK_BEGIN(cds)
	const packed_tran& parent_tran = anno.trans()[src.tran];
	const packed_exon& parent_exon = anno.exons()[parent_tran.exon0 + src.exon_index];
	phase  = src.phase;
	level  = parent_tran.level;  // inherit
	tsl    = parent_tran.tsl;    // inherit
	tran   = &parent_tran;
	exon   = &parent_exon;
GK_GENOME_ANNO_ELEM_UNPACK_END

////////////////////////////////////////////////
// utr_t
////////////////////////////////////////////////

GK_GENOME_ANNO_ELEM_UNPACK_BEGIN(utr)
	const packed_tran& parent_tran = anno.trans()[src.tran];
	const packed_exon& parent_exon = anno.exons()[parent_tran.exon0 + src.exon_index];
	tran   = &parent_tran;
	exon   = &parent_exon;
GK_GENOME_ANNO_ELEM_UNPACK_END

////////////////////////////////////////////////
// genome
////////////////////////////////////////////////

#ifdef _MSC_VER
#pragma warning (disable : 4355) // "Using 'this' in member initializer; safe in our case, so disable warning
#endif

genome_anno::genome_anno()
: _genes(*this)
, _trans(*this)
, _exons(*this)
, _intrs(*this)
, _cdss(*this)
, _utr5s(*this)
, _utr3s(*this)
{
}

const gene_table& genome_anno::genes() const { ensure_open(); return _genes; }
const tran_table& genome_anno::trans() const { ensure_open(); return _trans; }
const exon_table& genome_anno::exons() const { ensure_open(); return _exons; }
const intr_table& genome_anno::intrs() const { ensure_open(); return _intrs; }
const cds_table&  genome_anno::cdss()  const { ensure_open(); return _cdss;  }
const utr_table&  genome_anno::utr5s() const { ensure_open(); return _utr5s; }
const utr_table&  genome_anno::utr3s() const { ensure_open(); return _utr3s; }

bool genome_anno::empty() const noexcept { return std::empty(source()); }

void genome_anno::set_source(string sourcefile)
{
	GK_CHECK(!is_open(), runtime, "Cannot set source when file already open.");
	_sourcefile = std::move(sourcefile);
}

void genome_anno::open_on_demand() const
{
	// TODO: acquire lock here and check _fmap.is_open() again

	// Once we're inside the lock it's safe to use ncthis
	auto* ncthis = const_cast<genome_anno*>(this);
	ncthis->open();

	// TODO: release lock here, when falls out of scope
}

void genome_anno::close()
{
	_fmap.close();
}

string default_anno_sourcefile(string_view config, string_view data_dir)
{
	return prepend_dir(data_dir, format("{}.v{}.dganno", config, genome_anno::binary_version()));
}

END_NAMESPACE_GK
