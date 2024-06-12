/*
Copyright (C) 2016-2023 Deep Genomics Inc. All Rights Reserved.
*/
#include "py_genome.h"

#include "genome.h"
#include "genome_anno.h"
#include "genome_dna.h"
#include "py_genome_anno.h"
#include "py_genome_dna.h"
#include "py_util.h"
#include <structmember.h>

BEGIN_NAMESPACE_GK

GKPY_NEW_BEGIN(Genome)
	const char* config = "";
	static char* kwlist[] = { "config", nullptr };
	if (!PyArg_ParseTupleAndKeywords(args, kwds, "|s", kwlist, &config))
		return nullptr;

	// Construct the genome_t instance inplace
	GK_TENTATIVE_INPLACE_CONSTRUCT(genome, genome_t, config);
	// <--- INSERT SUB-OBJECT CONSTRUCTIONS HERE

	// Allow sub-objects to stay constructed on return
	GK_FINALIZE_CONSTRUCT(genome);
	// <--- INSERT SUB-OBJECT FINALIZATIONS HERE, IF CONSTRUCTED ABOVE

	self->genome_constructed = true;
GKPY_NEW_END

GKPY_INIT_EMPTY(Genome)

GKPY_TRAVERSE_BEGIN(Genome)
	GKPY_VISIT(dna);
	GKPY_VISIT(anno);
	GKPY_VISIT(genes);
	GKPY_VISIT(trans);
	GKPY_VISIT(exons);
	GKPY_VISIT(intrs);
	GKPY_VISIT(cdss);
	GKPY_VISIT(utr5s);
	GKPY_VISIT(utr3s);
	// <--- INSERT SUB-OBJECT VISITS HERE
GKPY_TRAVERSE_END

GKPY_CLEAR_BEGIN(Genome)
	GKPY_CLEAR(dna);
	GKPY_CLEAR(anno);
	GKPY_CLEAR(genes);
	GKPY_CLEAR(trans);
	GKPY_CLEAR(exons);
	GKPY_CLEAR(intrs);
	GKPY_CLEAR(cdss);
	GKPY_CLEAR(utr5s);
	GKPY_CLEAR(utr3s);
	// <--- INSERT SUB-OBJECT CLEARS HERE
GKPY_CLEAR_END

GKPY_DEALLOC_BEGIN(Genome)
	if (self->genome_constructed)
		destruct(&self->genome);
GKPY_DEALLOC_END

GKPY_GETATTRO_BEGIN(Genome)
	GKPY_GETATTR_MEMBER_ON_DEMAND(dna,  GenomeDNA,  "(O)", self);
	if (self->anno == nullptr && (!strcmp(attr, "anno") || !strcmp(attr, "annotation")
			|| !strcmp(attr, "genes")
			|| !strcmp(attr, "trans")
			|| !strcmp(attr, "transcripts")
			|| !strcmp(attr, "exons")
			|| !strcmp(attr, "introns")
			|| !strcmp(attr, "cdss")
		    || !strcmp(attr, "utr5s")
		    || !strcmp(attr, "utr3s")
			// <--- INSERT NAMES OF SHORTCUTS TO GenomeAnno ATTRIBUTES HERE
			)) {
		// Construct 'table' on-demand, and make shortcuts to its genes, transcripts, exons, etc
		GKPY_CONSTRUCT(anno, GenomeAnno, "(O)", self);
		self->genes = ((PyGenomeAnno*)self->anno)->genes; Py_INCREF(self->genes);
		self->trans = ((PyGenomeAnno*)self->anno)->trans; Py_INCREF(self->trans);
		self->exons = ((PyGenomeAnno*)self->anno)->exons; Py_INCREF(self->exons);
		self->intrs = ((PyGenomeAnno*)self->anno)->intrs; Py_INCREF(self->intrs);
		self->cdss  = ((PyGenomeAnno*)self->anno)->cdss;  Py_INCREF(self->cdss);
		self->utr5s = ((PyGenomeAnno*)self->anno)->utr5s; Py_INCREF(self->utr5s);
		self->utr3s = ((PyGenomeAnno*)self->anno)->utr3s; Py_INCREF(self->utr3s);
		// <--- INSERT SHORTCUTS TO GenomeAnno ATTRIBUTES HERE
	}
	GKPY_GETATTR_CASE("trans")            { GKPY_RETURN_INCREF(self->trans); }
	GKPY_GETATTR_CASE("anno")             { GKPY_RETURN_INCREF(self->anno); }
	GKPY_GETATTR_CASE("refg")             { return PyString_FromSV(self->genome.refg_name()); }
	GKPY_GETATTR_CASE("reference_genome") { return PyString_FromSV(self->genome.refg_name()); }
	GKPY_GETATTR_CASE("config")           { return PyString_FromSV(self->genome.config()); }
	GKPY_GETATTR_CASE("data_dir")         { return PyString_FromSV(self->genome.data_dir()); }
	return PyObject_GenericGetAttr(selfo, attro);
GKPY_GETATTRO_END

GKPY_SETATTRO_BEGIN(Genome)
	// Don't need to mention anything declared with GKPY_MEMBER_OBJECT, since they're already flagged as readonly
	GKPY_SETATTR_READONLY("trans");
	GKPY_SETATTR_READONLY("anno");
	GKPY_SETATTR_READONLY("refg");
	GKPY_SETATTR_READONLY("reference_genome");
	GKPY_SETATTR_READONLY("config");
	GKPY_SETATTR_READONLY("data_dir");
	return PyObject_GenericSetAttr(selfo, attro, value);
GKPY_SETATTRO_END

GKPY_MEMBERS_BEGIN(Genome)
	GKPY_MEMBER_OBJECT(Genome, dna,   "dna",         nullptr)
	GKPY_MEMBER_OBJECT(Genome, anno,  "annotation",  nullptr)
	GKPY_MEMBER_OBJECT(Genome, genes, "genes",       nullptr)
	GKPY_MEMBER_OBJECT(Genome, trans, "transcripts", nullptr)
	GKPY_MEMBER_OBJECT(Genome, exons, "exons",       nullptr)
	GKPY_MEMBER_OBJECT(Genome, intrs, "introns",     nullptr)
	GKPY_MEMBER_OBJECT(Genome, cdss,   "cdss",       nullptr)
	GKPY_MEMBER_OBJECT(Genome, utr5s,  "utr5s",      nullptr)
	GKPY_MEMBER_OBJECT(Genome, utr3s,  "utr3s",      nullptr)
GKPY_MEMBERS_END

// intentionally not documenting as a public Python api
GKPY_OMETHOD_BEGIN(Genome, _refg_hash)
		const char*  refg_name;
		static char* kwlist[] = {"refg_name", nullptr};
		if (!PyArg_ParseTupleAndKeywords(args, kwds, "s", kwlist, &refg_name))
			return nullptr;

		return PyLong_FromUnsignedLongLong(fnv1a_hash64(refg_name));
GKPY_OMETHOD_END


GKPY_METHODS_BEGIN(Genome)
		GKPY_METHOD_ENTRY(Genome, _refg_hash, METH_VARARGS | METH_KEYWORDS | METH_STATIC, nullptr)
GKPY_METHODS_END

GKPY_TYPEOBJ_BEGIN(Genome)
	tp_flags |= Py_TPFLAGS_BASETYPE | Py_TPFLAGS_HAVE_GC;
	tp_new = PyGenome_New;
	tp_init = PyGenome_Init;
	tp_dealloc = PyGenome_Dealloc;
	tp_getattro = PyGenome_GetAttro;
	tp_setattro = PyGenome_SetAttro;
	tp_traverse = PyGenome_Traverse;
	tp_clear = PyGenome_Clear;
	tp_richcompare = PyGenericValue_RichCompare<PyGenome>;
	tp_members = PyGenome_Members;
	tp_methods = PyGenome_Methods;
GKPY_TYPEOBJ_END

END_NAMESPACE_GK
