/*
Copyright (C) 2016-2023 Deep Genomics Inc. All Rights Reserved.
*/
#pragma once
#ifndef __GENOME_KIT_PY_VARIANT_TABLE_H__
#define __GENOME_KIT_PY_VARIANT_TABLE_H__

#include "py_interval.h"
#include "py_table.h"
#include "py_util.h"
#include "variant_table.h"
#include <map>

BEGIN_NAMESPACE_GK

struct PyVariantTable;
struct PyVCFTable;

///////////////////////////////////////////////////////////////////
// PyVariant
///////////////////////////////////////////////////////////////////

struct avariant_t : public ainterval_t {
	NOCOPY(avariant_t)
	avariant_t(chrom_t chrom, pos_t start, const char* ref, size_t reflen, const char* alt, size_t altlen, refg_t refg);
	INLINE ~avariant_t()
	{
		// small-string optimization for ref/alt:
		//  - if <= 8 bytes, stored in 8-bytes of [anchor,anchor_offset]
		//  - if >  8 bytes, ref and alt both stored in single new[] pointed to by ref
		if (ref && ref != (char*)&anchor)
			delete[] ref;
	}

	char* ref;
	char* alt;
};

GKPY_VALUE_SUBTYPE_BEGIN(Variant, Interval, packed_variant)
	using Table          = PyVariantTable; /* The Python type wrapping the C++ table type */
	using table_type     = const variant_table;  /* The C++ table type */
	using unpacked_value = variant_t;      /* The unpacked version of this variant */
	unpacked_value unpack();
	// avariant_t is constructable from python,
	// Assumes that PyVariant is either as_ptr or as_value, where the value type is avariant_t instead of value_t in
	// case of as_value Rationale for the constraint: above code makes sure PyVariant still works with PyVariantTable
	// Other changes: requires setting the correct tp_basicsize such that avariant_t fits inside PyVariant inplace
	using pyconstructed_t = avariant_t;
	INLINE bool is_pyconstructed()
	{
		return !as_ptr
			   // check to see if we subclassed as a ptr type (and as_ptr is not yet initialized)
			   && ((PyTypeObject*)(PyObject_Type((PyObject*)this)))->tp_basicsize
					  >= __value_type_sizes<pyconstructed_t>::as_value;
	}
	INLINE pyconstructed_t& pyconstructed_value()
	{
		GK_DBASSERT(is_pyconstructed());
		return _get<pyconstructed_t>(0);
	}
	INLINE static pyconstructed_t& pyconstructed_value(PyObject * obj)
	{
		return ((PyVariant*)obj)->pyconstructed_value();
	}
	PyVariant() = default;
	INLINE ~PyVariant()
	{
		// FastDealloc<PyVariant> calls this as soon as the variant has refcount of zero.
		// If was constructed from Python, destruct the inplace avariant_t instance (possibly freeing ref/alt memory).
		if (is_pyconstructed())
			destruct(&pyconstructed_value());
	}
	NOCOPY(PyVariant);
	INLINE variant_t as_variant()
	{
		if (is_pyconstructed()) {
			avariant_t& v = pyconstructed_value();
			return variant_t(v, v.ref, v.alt);
		}
		return unpack();
	}
GKPY_VALUE_SUBTYPE_END

///////////////////////////////////////////////////////////////////
// PyVariantTable
///////////////////////////////////////////////////////////////////

GKPY_SUBTYPE_BEGIN(VariantTable, IntervalTable<PyVariant>)
	PyObject* owner;
GKPY_SUBTYPE_END

INLINE variant_t PyVariant::unpack() { return unpacked_value(value(), *((Table*)source())->table); }

///////////////////////////////////////////////////////////////////
// PyVCFVariant
///////////////////////////////////////////////////////////////////

GKPY_VALUE_SUBTYPE_BEGIN(VCFVariant, Variant, packed_variant)
	using Table          = PyVCFTable; /* The Python type wrapping the C++ table type */
	using table_type     = vcf_table;  /* The C++ table type */
	using unpacked_value = variant_t;  /* The unpacked version of this gene */
GKPY_VALUE_SUBTYPE_END

///////////////////////////////////////////////////////////////////
// PyVCFTable
///////////////////////////////////////////////////////////////////

GKPY_SUBTYPE_BEGIN(VCFTable, IntervalTable<PyVCFVariant>)
	static PyAutoRef numpy_memmap_function;

	std::map<const vcf_table::field_col_t*, std::vector<PyAutoRef>>* pystr_pool;
	PyObject* pystr_filename;

	PyObject* get_col_attr(const packed_variant* variant, const char* id);
	PyObject* get_col_attr_str(size_t index, const vcf_table::field_col_t* col);

GKPY_SUBTYPE_END

END_NAMESPACE_GK

#endif // __GENOME_KIT_PY_VARIANT_TABLE_H__
