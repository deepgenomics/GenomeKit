/*
Copyright (C) 2016-2023 Deep Genomics Inc. All Rights Reserved.
*/
#pragma once
#ifndef __GENOME_KIT_PY_JRALIGN_H__
#define __GENOME_KIT_PY_JRALIGN_H__

#include "jralign.h"
#include "py_util.h"
#include "py_table.h"

BEGIN_NAMESPACE_GK

GKPY_TYPE_BEGIN(JReadAlignments)
	junction_read_alignments jraligns;
	bool jraligns_constructed;  // Flag indicating whether __new__ succeeded;
	PyObject *juncs;            // PyJRAlignsTable instance
	PyObject *variants;         // PyVariantTable instance
GKPY_TYPE_END

///////////////////////////////////////////////////////////////////
// PyRAlignsTable<PyJRAligns>, etc
///////////////////////////////////////////////////////////////////

#define GKPY_DECLARE_READ_TABLE(elem_name, celem_name) \
	GKPY_VALUE_SUBTYPE_BEGIN(elem_name, Interval, packed_##celem_name)                                          \
		using Table          = PyIntervalTable<Py##elem_name>; /* The Python type wrapping the C++ table type */       \
		using table_type     = const celem_name##_table;       /* The C++ table type */                                \
		using unpacked_value = celem_name##_t;                 /* The unpacked version of this element */              \
		INLINE const unpacked_value unpack() { return unpacked_value(value(), *((Table*)source())->table); }    \
		INLINE const PyJReadAlignments &junction_read_alignments(); \
	GKPY_VALUE_SUBTYPE_END \
\
	GKPY_SUBTYPE_BEGIN(elem_name##Table, IntervalTable<Py##elem_name>)                                          \
		PyObject *owner; \
	GKPY_SUBTYPE_END \
\
	const PyJReadAlignments &Py##elem_name::junction_read_alignments() { return *((PyJReadAlignments*)((Py##elem_name##Table*)source())->owner); }

GKPY_DECLARE_READ_TABLE(JRAligns, jraligns)

GKPY_VALUE_TYPE_BEGIN(JRAlign, jralign_t)
	jralign_t dummy;
	PyObject *owner;
GKPY_VALUE_TYPE_END

// Aliases
using PyJunctionReadAlignmentsTable = PyJRAlignsTable;
using PyJunctionReadAlignments = PyJRAligns;
using PyJunctionReadAlignment = PyJRAlign;

END_NAMESPACE_GK

#endif // __GENOME_KIT_PY_JRALIGN_H__
