/*
Copyright (C) 2016-2023 Deep Genomics Inc. All Rights Reserved.
*/
#pragma once
#ifndef __GENOME_KIT_PY_GENOME_ANNO_H__
#define __GENOME_KIT_PY_GENOME_ANNO_H__

#include "genome_anno.h"
#include "py_util.h"
#include "py_table.h"
#include <unordered_map>

BEGIN_NAMESPACE_GK

void Init_GenomeAnno_PyStrings();

GKPY_DECLARE_STRINGTABLE(BioType, biotype_t, biotype, num_biotype)

GKPY_TYPE_BEGIN(GenomeAnno)
	const genome_anno* anno;
	PyObject* owner;       // PyGenome instance
	PyObject* genes;       // PyGeneTable instance
	PyObject* trans;       // PyTranTable instance
	PyObject* exons;       // PyExonTable instance
	PyObject* intrs;       // PyIntrTable instance
	PyObject* cdss;        // PyCdsTable instance
	PyObject* utr5s;       // PyUtrTable instance
	PyObject* utr3s;       // PyUtrTable instance
GKPY_TYPE_END

///////////////////////////////////////////////////////////////////
// PyGenomeAnnoTable<T>
///////////////////////////////////////////////////////////////////

// T = PyGene, PyExon, etc
template <typename T>
GKPY_SUBTYPE_BEGIN(GenomeAnnoTable, IntervalTable<T>)
	PyObject* owner;
	std::unordered_map<std::string, index_t> index_by_id;
	GKPY_SUBTYPE_END

#define GKPY_DECLARE_GENOME_ANNO_TABLE(elem_name, celem_name) \
	GKPY_VALUE_SUBTYPE_BEGIN(elem_name, Interval, packed_##celem_name) \
		using Table          = PyGenomeAnnoTable<Py##elem_name>;  /* The Python type wrapping the C++ table type */ \
		using table_type     = const celem_name##_table;          /* The C++ table type */                          \
		using unpacked_value = celem_name##_t; \
		INLINE const unpacked_value unpack() { \
			return unpacked_value{value(), *((Table*)source())->table}; \
		} \
		INLINE PyGenomeAnno& anno() { \
			return *((PyGenomeAnno*)((Table*)source())->owner); \
		} \
	GKPY_VALUE_SUBTYPE_END \
\
	GKPY_SUBTYPE_BEGIN(elem_name##Table, GenomeAnnoTable< Py##elem_name >) \
	GKPY_SUBTYPE_END

	///////////////////////////////////////////////////////////////////
// PyGene, PyTran, etc
///////////////////////////////////////////////////////////////////

	GKPY_DECLARE_GENOME_ANNO_TABLE(Gene, gene)
	GKPY_DECLARE_GENOME_ANNO_TABLE(Tran, tran)
	GKPY_DECLARE_GENOME_ANNO_TABLE(Exon, exon)
	GKPY_DECLARE_GENOME_ANNO_TABLE(Intr, intr)
	GKPY_DECLARE_GENOME_ANNO_TABLE(Cds,  cds)
	GKPY_DECLARE_GENOME_ANNO_TABLE(Utr, utr)

// Aliases
using PyGenomeAnnotation = PyGenomeAnno;
using PyTranscript = PyTran;
using PyTranscriptTable = PyTranTable;
using PyIntron = PyIntr;
using PyIntronTable = PyIntrTable;

END_NAMESPACE_GK

#endif // __GENOME_KIT_PY_GENOME_ANNO_H__
