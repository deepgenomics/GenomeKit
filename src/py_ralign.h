/*
Copyright (C) 2016-2023 Deep Genomics Inc. All Rights Reserved.
*/
#pragma once
#ifndef __GENOME_KIT_PY_RALIGN_H__
#define __GENOME_KIT_PY_RALIGN_H__

#include "ralign.h"
#include "py_util.h"
#include "py_table.h"

BEGIN_NAMESPACE_GK

GKPY_TYPE_BEGIN(ReadAlignments)
	read_alignments raligns;
	bool raligns_constructed ;      // Flag indicating whether __new__ succeeded;
	PyObject* junctions;            // PyJunctionTable instance
	PyObject* alignments;           // PyAlignmentsTable instance
	PyObject* matches;              // PyAlignmentMatchTable instance
	PyObject* variants;             // PyVariantTable instance
GKPY_TYPE_END

///////////////////////////////////////////////////////////////////
// PyJunction,      PyAlignment,      PyAlignmentMatch
// PyJunctionTable, PyAlignmentTable, PyAlignmentMatchTable
///////////////////////////////////////////////////////////////////

#define GKPY_DECLARE_READ_TABLE2(elem_name, celem_name) \
	GKPY_VALUE_SUBTYPE_BEGIN(elem_name, Interval, packed_##celem_name)                                          \
		using Table          = PyIntervalTable<Py##elem_name>;   /* The Python type wrapping the C++ table type */       \
		using table_type     = const celem_name##_table;               /* The C++ table type */                                \
		using unpacked_value = celem_name##_t;                   /* The unpacked version of this element */              \
		INLINE const unpacked_value unpack() { return unpacked_value(value(), ((Table *)source())->table->raligns); } \
		INLINE const PyReadAlignments &read_alignments(); \
	GKPY_VALUE_SUBTYPE_END \
\
	GKPY_SUBTYPE_BEGIN(elem_name##Table, IntervalTable<Py##elem_name>)                                          \
		PyObject *owner; \
	GKPY_SUBTYPE_END \
\
	const PyReadAlignments &Py##elem_name::read_alignments() { return *((PyReadAlignments *)((Py##elem_name##Table *)source())->owner); }

GKPY_DECLARE_READ_TABLE2(Junction,          junction)
GKPY_DECLARE_READ_TABLE2(Alignment,         align)
GKPY_DECLARE_READ_TABLE2(AlignmentMatch,    align_match)

END_NAMESPACE_GK

#endif // __GENOME_KIT_PY_RALIGN_H__
