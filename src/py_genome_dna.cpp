/*
Copyright (C) 2016-2023 Deep Genomics Inc. All Rights Reserved.
*/
#include "py_util.h" // Include first to avoid POSIX warnings due to Python header inclusion order
#include "py_genome_dna.h"
#include "py_genome.h"
#include "py_interval.h"
#include "genome_dna.h"

BEGIN_NAMESPACE_GK

GKPY_INIT_BEGIN(GenomeDNA)
	// Allow GenomeDNA to be constructed one of two ways:
	//    GenomeDNA(genome : Genome)
	//    GenomeDNA(twobitfile: str)
	PyObject* arg = nullptr;
	if (!PyArg_ParseTuple(args, "O", &arg))
		return -1;
	if (PyString_Check(arg)) {
		// If argument was a string, set the source file directly
		self->dna = new genome_dna();
		self->dna->set_source(PyString_AS_STRING(arg));
	} else if (PyObject_IsInstance(arg, (PyObject*)PyGenome::Type)) {
		// If argument was a Genome object, then point to the pre-existing genome_dna instance
		Py_INCREF(arg);
		self->owner = arg;
		self->dna = &as_mutable(((PyGenome*)self->owner)->genome.dna());
	}
	self->dna->ensure_open(); // Open on-demand
GKPY_INIT_END

GKPY_DEALLOC_OWNED_BEGIN(GenomeDNA)
	// If we don't have an owner, then the genome_dna instance must be deleted
	if (self->dna && !self->owner) {
		delete self->dna;
		self->dna = nullptr;
	}
GKPY_DEALLOC_OWNED_END

GKPY_TRAVERSE_OWNED_BEGIN(GenomeDNA)
GKPY_TRAVERSE_OWNED_END

GKPY_CLEAR_OWNED_BEGIN(GenomeDNA)
GKPY_CLEAR_OWNED_END

GKPY_GETATTRO_BEGIN(GenomeDNA)
	GKPY_GETATTR_CASE("reference_genome")  return PyString_FromSV(self->dna->refg_name());
	GKPY_GETATTR_CASE("refg")              return PyString_FromSV(self->dna->refg_name());
	GKPY_GETATTR_CASE("filename")          return PyString_FromSV(self->dna->source());
	return PyObject_GenericGetAttr(selfo, attro);
GKPY_GETATTRO_END

GKPY_SETATTRO_BEGIN(GenomeDNA)
	GKPY_SETATTR_READONLY("reference_genome");
	GKPY_SETATTR_READONLY("refg");
	GKPY_SETATTR_READONLY("filename");
	return PyObject_GenericSetAttr(selfo, attro, value);
GKPY_SETATTRO_END

GKPY_OMETHOD_BEGIN(GenomeDNA, Call)
	char* arg_names[] = { "interval", "allow_outside_chromosome", nullptr };
	PyObject* interval_arg;
	int allow_outside_chromosome = 1;
	if (!PyArg_ParseTupleAndKeywords(args, kwds, "O|p", arg_names, &interval_arg, &allow_outside_chromosome))
		 return nullptr;

	// Extract DNA for an Interval arg
	const interval_t& c = PyAsInterval(interval_arg);
	return PyString_FromSV((*self->dna)(c, allow_outside_chromosome != 0));
GKPY_OMETHOD_END

GKPY_TYPEOBJ_BEGIN(GenomeDNA)
	tp_flags |= Py_TPFLAGS_BASETYPE | Py_TPFLAGS_HAVE_GC;
	tp_new = PyType_GenericNew;
	tp_init = PyGenomeDNA_Init;
	tp_dealloc = PyGenomeDNA_Dealloc;
	tp_getattro = PyGenomeDNA_GetAttro;
	tp_setattro = PyGenomeDNA_SetAttro;
	tp_traverse = PyGenomeDNA_Traverse;
	tp_clear = PyGenomeDNA_Clear;
	tp_call = PyGenomeDNA_Call;
GKPY_TYPEOBJ_END

END_NAMESPACE_GK
