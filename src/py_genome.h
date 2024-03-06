/*
Copyright (C) 2016-2023 Deep Genomics Inc. All Rights Reserved.
*/
#pragma once
#ifndef __GENOME_KIT_PY_GENOME_H__
#define __GENOME_KIT_PY_GENOME_H__

#include "py_util.h"
#include "genome.h"

BEGIN_NAMESPACE_GK

GKPY_TYPE_BEGIN(Genome)
	// Lets us use PyGenericValue_RichCompare from py_util.h, for example
	using value_t = genome_t;
	INLINE static value_t& value(PyObject* obj) { return ((PyGenome*)obj)->genome; }

	// The main C++ instance that implements all the underlying
	// loading, indexing, searching, and serialization.
	genome_t genome;
	bool genome_constructed; // Flag indicating whether __new__ succeeded;
	                         // needed to distinguish between PyGenome dealloc
	                         // after an error during genome_t constructor
	                         // (in which case genome_t destructor should NOT
	                         // be called) and case where PyGenome dealloc is
	                         // occurs after genome object constructed
	                         // (in which case destructor SHOULD be called)

	// Python "sub-objects" that each refer to a different aspect
	// of the genome implementation, but are all owned by this PyGenome
	// instance and all keep references to this PyGenome instance
	// (i.e. you can do "dna = Genome(...).dna" and the dna reference
	// will keep the Genome instance alive)
	// These objects are lightweight in the sense that they merely wrap
	// and refer to a corresponding C++ object.
	// They are dynamically created the first time the user
	// accesses them by name.
	PyObject* dna;         // PyGenomeDNA instance

	PyObject* anno;        // PyGenomeAnno instance
	PyObject* genes;       // PyGeneTable instance owned and created by 'annotation'
	PyObject* trans;       // PyTranTable instance owned and created by 'annotation'
	PyObject* exons;       // PyExonTable instance owned and created by 'annotation'
	PyObject* intrs;       // PyIntrTable instance owned and created by 'annotation'
	PyObject* cdss;        // PyCdsTable instance owned and created by 'annotation'
	PyObject* utr5s;       // PyUtrTable instance owned and created by 'annotation'
	PyObject* utr3s;       // PyUtrTable instance owned and created by 'annotation'

GKPY_TYPE_END

END_NAMESPACE_GK

#endif // __GENOME_KIT_PY_GENOME_H__
