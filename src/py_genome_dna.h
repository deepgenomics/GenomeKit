/*
Copyright (C) 2016-2023 Deep Genomics Inc. All Rights Reserved.
*/
#pragma once
#ifndef __GENOME_KIT_PY_GENOME_DNA_H__
#define __GENOME_KIT_PY_GENOME_DNA_H__

#include "py_util.h"

BEGIN_NAMESPACE_GK

GKPY_TYPE_BEGIN(GenomeDNA)
	PyObject* owner;        // The object in which we're embedded
	class genome_dna* dna;  // Pointer to the DNA object we're wrapping
GKPY_TYPE_END

END_NAMESPACE_GK

#endif // __GENOME_KIT_PY_GENOME_DNA_H__
