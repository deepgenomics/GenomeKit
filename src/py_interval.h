/*
Copyright (C) 2016-2023 Deep Genomics Inc. All Rights Reserved.
*/
#pragma once
#ifndef __GENOME_KIT_PY_INTERVAL_H__
#define __GENOME_KIT_PY_INTERVAL_H__

#include "py_util.h"
#include "interval.h"

BEGIN_NAMESPACE_GK

class genome_t;

const pos_t invalid_pos = (pos_t)0x80000000; // largest negative int32

///////////////////////////////////////////////////////////////////
// PyInterval
///////////////////////////////////////////////////////////////////

#pragma pack(push, 1)
struct ainterval_t: public interval_t {
	pos_t anchor{invalid_pos};
	pos_t anchor_offset{};
};
#pragma pack(pop)

// ainterval_t memory layout: offset name (size+padding)
//  0 pos5          (4)
//  4 pos3          (4)
//  8 refg          (8)
// 16 chrom         (4)
// 20 strand        (1)
// 21 anchor        (4)
// 25 anchor_offset (4)
// 29 total

GKPY_VALUE_TYPE_BEGIN(Interval, ainterval_t)
	static PyObject* create(const interval_t& v, pos_t anchor=invalid_pos, pos_t anchor_offset=0);
	INLINE pos_t get_anchor()        { return as_ptr ? invalid_pos : value().anchor; }
	INLINE pos_t get_anchor_offset() { return as_ptr ? 0 : value().anchor_offset; }
GKPY_VALUE_TYPE_END

///////////////////////////////////////////////////////////////////
// global string tables
///////////////////////////////////////////////////////////////////

GKPY_DECLARE_STRINGTABLE(Strand, strand_t, strand, num_strand)

void Init_Interval_PyStrings();

INLINE interval_t PyAsInterval(PyObject* arg)
{
	if (PyInterval::check(arg)) return PyInterval::value(arg);
	else GK_THROW(type, "Expected argument of type Interval");
}

INLINE ainterval_t PyAsAnchoredInterval(PyObject* arg)
{
	if (PyInterval::check(arg)) {
		ainterval_t ret   = { (interval_t)PyInterval::value(arg) };
		ret.anchor        = ((PyInterval*)arg)->get_anchor();
		ret.anchor_offset = ((PyInterval*)arg)->get_anchor_offset();
		return ret;
	} else
		GK_THROW(type, "Expected argument of type Interval");
}

// Helper methods
pos_t as_pos(PyObject* arg, const interval_t& interval, const char* name);
refg_t as_refg(PyObject* arg);
const genome_t& as_genome(PyObject* arg);
// Python Interval.end5/end3 returns empty closed intervals, which is slightly different from C++ interval which returns
// 1-bp open intervals
INLINE static interval_t get_end5(const interval_t& i) { return i.end5().expand(0, -1); }
INLINE static interval_t get_end3(const interval_t& i) { return i.end3().expand(-1, 0); }

END_NAMESPACE_GK

#endif // __GENOME_KIT_PY_INTERVAL_H__
