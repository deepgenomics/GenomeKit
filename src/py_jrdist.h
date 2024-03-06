/*
Copyright (C) 2016-2023 Deep Genomics Inc. All Rights Reserved.
*/
#pragma once
#ifndef __GENOME_KIT_PY_JRDIST_H__
#define __GENOME_KIT_PY_JRDIST_H__

#include "jrdist.h"
#include "py_util.h"
#include "py_table.h"
#include "py_jralign.h"

BEGIN_NAMESPACE_GK

GKPY_TYPE_BEGIN(ReadDistributions)
	read_distributions rdists;
	bool rdists_constructed; // Flag indicating whether __new__ succeeded;
	PyObject* juncs;         // PyJRDist instance
GKPY_TYPE_END

///////////////////////////////////////////////////////////////////
// PyJRDistTable<PyJRDist>, etc
///////////////////////////////////////////////////////////////////

GKPY_DECLARE_READ_TABLE(JRDist, jrdist)

GKPY_VALUE_TYPE_BEGIN(JRCount, jrcount_t)
GKPY_VALUE_TYPE_END

// Aliases
using PyJunctionReadDistributionTable = PyJRDistTable;
using PyJunctionReadDistribution      = PyJRDist;
using PyJunctionReadCount             = PyJRCount;

END_NAMESPACE_GK

#endif // __GENOME_KIT_PY_JRDIST_H__
