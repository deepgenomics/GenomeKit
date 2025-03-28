/*
Copyright (C) 2016-2023 Deep Genomics Inc. All Rights Reserved.
*/
#include "py_util.h" // Include first to avoid POSIX warnings due to Python header inclusion order

#include "py_genome.h"
#include "py_ralign.h"
#include "py_variant_table.h"
#include "util.h"
#include <structmember.h>

BEGIN_NAMESPACE_GK

void validate_JunctionTable(const PyAsPtrSource* self)
{
	GK_CHECK2(((PyJunctionTable*)self)->table->valid(), file,
			 "Junctions have been invalidated by ReadAlignments.close or with statement.");
}

void validate_AlignmentTable(const PyAsPtrSource* self)
{
	GK_CHECK2(((PyAlignmentTable*)self)->table->valid(), file,
			 "Alignments have been invalidated by ReadAlignments.close or with statement.");
}

void validate_AlignmentMatchTable(const PyAsPtrSource* self)
{
	GK_CHECK2(((PyAlignmentMatchTable*)self)->table->valid(), file,
			 "AlignmentMatches have been invalidated by ReadAlignments.close or with statement.");
}

#define GKPY_READ_ALIGNS_TABLE_BEGIN(name, table_name) \
	GKPY_NEW_OWNED_BEGIN(name##Table, ReadAlignments) \
		self->table = &((PyReadAlignments*)self->owner)->raligns.table_name(); \
		self->set_validator(&validate_##name##Table); \
	GKPY_NEW_OWNED_END \
\
	GKPY_DEALLOC_OWNED_BEGIN(name##Table) \
	GKPY_DEALLOC_OWNED_END \
\
	GKPY_SUBTYPEOBJ_BEGIN(name##Table) \
		tp_flags |= Py_TPFLAGS_BASETYPE; /* Allow Py_TPFLAGS_HAVE_GC to be inherited so that tp_traverse/tp_clear get inherited too */ \
		tp_new = Py##name##Table_New; \
		tp_dealloc = Py##name##Table_Dealloc; \

#define GKPY_READ_ALIGNS_TABLE_END GKPY_SUBTYPEOBJ_END

///////////////////////////////////////////////////////////////////
// PyJunction, PyJunctionTable
///////////////////////////////////////////////////////////////////

GKPY_GETATTRO_BEGIN(Junction)
	GKPY_GETATTR_CASE("interval")   { return PyInterval::create(self->value()); }
	GKPY_GETATTR_CASE("_num_alignments")  { return PyInt_FromLong(self->value().num_aligns); }
	return PyJunction::base_type::DefaultType->tp_getattro(selfo, attro);
GKPY_GETATTRO_END

GKPY_SETATTRO_BEGIN(Junction)
	GKPY_SETATTR_READONLY("interval");
	GKPY_SETATTR_READONLY("_num_alignments");
	return PyJunction::base_type::DefaultType->tp_setattro(selfo, attro, value);
GKPY_SETATTRO_END

GKPY_METHOD_BEGIN_ONEARG(Junction, _get_alignment)
	GKPY_TYPECHECK_BUILTIN(arg, PyInt);
	long index = PyInt_AS_LONG(arg);

	const junction_t& junc = self->unpack();
	GKPY_INDEXCHECK(index, junc.num_aligns);
	return PyTable_GetItem<PyAlignment>(self->read_alignments().alignments, junc.aligns[index]);
GKPY_METHOD_END

GKPY_METHODS_BEGIN(Junction)
	GKPY_METHOD_ENTRY(Junction, _get_alignment, METH_O, nullptr)
GKPY_METHODS_END

GKPY_VALUE_SUBTYPEOBJ_BEGIN(Junction, as_ptr)
	tp_flags |= Py_TPFLAGS_BASETYPE;
	tp_init = gkpy_invalid_init;
	tp_getattro = PyJunction_GetAttro;
	tp_setattro = PyJunction_SetAttro;
	tp_methods = PyJunction_Methods;
	tp_richcompare = PyGenericValue_RichCompare<PyJunction>; // for hash == comparison
	tp_hash = PyGenericValue_Hash<PyJunction>;
GKPY_VALUE_SUBTYPEOBJ_END

GKPY_READ_ALIGNS_TABLE_BEGIN(Junction, junctions)
GKPY_READ_ALIGNS_TABLE_END


///////////////////////////////////////////////////////////////////
// PyAlignment, PyAlignmentTable
///////////////////////////////////////////////////////////////////

GKPY_GETATTRO_BEGIN(Alignment)
	GKPY_GETATTR_CASE("matches")      { return PyTuple_FromSizedElemRange<PyAlignmentMatch>(self->read_alignments().matches, self->unpack().matches); }
	GKPY_GETATTR_CASE("interval")     { return PyInterval::create(self->value()); }
	return PyAlignment::base_type::DefaultType->tp_getattro(selfo, attro);
GKPY_GETATTRO_END

GKPY_SETATTRO_BEGIN(Alignment)
	GKPY_SETATTR_READONLY("matches");
	GKPY_SETATTR_READONLY("interval");
	return PyAlignment::base_type::DefaultType->tp_setattro(selfo, attro, value);
GKPY_SETATTRO_END

GKPY_VALUE_SUBTYPEOBJ_BEGIN(Alignment, as_ptr)
	tp_flags |= Py_TPFLAGS_BASETYPE;
	tp_init = gkpy_invalid_init;
	tp_getattro = PyAlignment_GetAttro;
	tp_setattro = PyAlignment_SetAttro;
	tp_richcompare = PyGenericValue_RichCompare<PyAlignment>; // for hash == comparison
	tp_hash = PyGenericValue_Hash<PyAlignment>;
GKPY_VALUE_SUBTYPEOBJ_END

GKPY_READ_ALIGNS_TABLE_BEGIN(Alignment, alignments)
GKPY_READ_ALIGNS_TABLE_END

///////////////////////////////////////////////////////////////////
// PyAlignMatch, PyAlignMatchTable
///////////////////////////////////////////////////////////////////

GKPY_GETATTRO_BEGIN(AlignmentMatch)
	GKPY_GETATTR_CASE("interval")       { return PyInterval::create(self->value()); }
	GKPY_GETATTR_CASE("variants") {
		const align_match_t match = self->unpack();
		int size = match.num_variants;
		PyObject* tuple = PyTuple_New(size);
		GKPY_TAKEREF(tuple);
		for (int i = 0; i < size; ++i) {
			PyTuple_SET_ITEM(tuple, i, PyTable_GetItem<PyVariant>(self->read_alignments().variants, match.variants[i]));
		}
		GKPY_FORGETREF(tuple); return tuple;
	}
	return PyAlignmentMatch::base_type::DefaultType->tp_getattro(selfo, attro);
GKPY_GETATTRO_END

GKPY_SETATTRO_BEGIN(AlignmentMatch)
	GKPY_SETATTR_READONLY("interval");
	GKPY_SETATTR_READONLY("variants");
	return PyAlignmentMatch::base_type::DefaultType->tp_setattro(selfo, attro, value);
GKPY_SETATTRO_END

GKPY_VALUE_SUBTYPEOBJ_BEGIN(AlignmentMatch, as_ptr)
	tp_flags |= Py_TPFLAGS_BASETYPE;
	tp_init = gkpy_invalid_init;
	tp_getattro = PyAlignmentMatch_GetAttro;
	tp_setattro = PyAlignmentMatch_SetAttro;
	tp_richcompare = PyGenericValue_RichCompare<PyAlignmentMatch>; // for hash == comparison
	tp_hash = PyGenericValue_Hash<PyAlignmentMatch>;
GKPY_VALUE_SUBTYPEOBJ_END

GKPY_READ_ALIGNS_TABLE_BEGIN(AlignmentMatch, matches)
GKPY_READ_ALIGNS_TABLE_END

///////////////////////////////////////////////////////////////////
// PyReadAlignments
///////////////////////////////////////////////////////////////////

GKPY_NEW_BEGIN(ReadAlignments)
	// We expect to be called with __new__(infile) infile is the
	// path to the .rdist file we should open
	const char* infile    = nullptr;
	static char* kwlist[] = { "infile", nullptr };
	if (!PyArg_ParseTupleAndKeywords(args, kwds, "s", kwlist, &infile))
		return nullptr;

	// Construct the genome_t instance inplace
	GK_TENTATIVE_INPLACE_CONSTRUCT(raligns, read_alignments);
	self->raligns.set_source(infile);
	self->raligns.ensure_open();
	// <--- INSERT SUB-OBJECT CONSTRUCTIONS HERE

	// Construct sub-objects
	GKPY_TENTATIVE_CONSTRUCT(junctions,     JunctionTable,       "(O)", self);
	GKPY_TENTATIVE_CONSTRUCT(alignments,    AlignmentTable,      "(O)", self);
	GKPY_TENTATIVE_CONSTRUCT(matches,       AlignmentMatchTable, "(O)", self);
	GKPY_TENTATIVE_CONSTRUCT(variants,      VariantTable,        "(On)", self, (Py_ssize_t)&self->raligns.variants());

	// Allow sub-objects to stay constructed on return
	GK_FINALIZE_CONSTRUCT(raligns);
	GK_FINALIZE_CONSTRUCT(junctions);
	GK_FINALIZE_CONSTRUCT(alignments);
	GK_FINALIZE_CONSTRUCT(matches);
	GK_FINALIZE_CONSTRUCT(variants);
	// <--- INSERT SUB-OBJECT FINALIZATIONS HERE, IF CONSTRUCTED ABOVE

	self->raligns_constructed = true;
GKPY_NEW_END

GKPY_INIT_BEGIN(ReadAlignments)
GKPY_INIT_END

GKPY_TRAVERSE_BEGIN(ReadAlignments)
	GKPY_VISIT(junctions);
	GKPY_VISIT(alignments);
	GKPY_VISIT(matches);
	GKPY_VISIT(variants);
	// <--- INSERT SUB-OBJECT VISITS HERE
GKPY_TRAVERSE_END

GKPY_CLEAR_BEGIN(ReadAlignments)
	GKPY_CLEAR(junctions);
	GKPY_CLEAR(alignments);
	GKPY_CLEAR(matches);
	GKPY_CLEAR(variants);
	// <--- INSERT SUB-OBJECT CLEARS HERE
GKPY_CLEAR_END

GKPY_DEALLOC_BEGIN(ReadAlignments)
	if (self->raligns_constructed)
		destruct(&self->raligns);
GKPY_DEALLOC_END

GKPY_GETATTRO_BEGIN(ReadAlignments)
	GKPY_GETATTR_CASE("junctions")      { GKPY_RETURN_INCREF(self->junctions); }
	GKPY_GETATTR_CASE("alignments")     { GKPY_RETURN_INCREF(self->alignments); }
	GKPY_GETATTR_CASE("matches")        { GKPY_RETURN_INCREF(self->matches); }
	GKPY_GETATTR_CASE("variants")       { GKPY_RETURN_INCREF(self->variants); }
	GKPY_GETATTR_CASE("filename")       { return PyString_FromSV(self->raligns.source()); }
	return PyObject_GenericGetAttr(selfo, attro);
GKPY_GETATTRO_END

GKPY_SETATTRO_BEGIN(ReadAlignments)
	GKPY_SETATTR_READONLY("junctions");
	GKPY_SETATTR_READONLY("alignments");
	GKPY_SETATTR_READONLY("matches");
	GKPY_SETATTR_READONLY("variants");
	GKPY_SETATTR_READONLY("filename");
	return PyObject_GenericSetAttr(selfo, attro, value);
GKPY_SETATTRO_END

static PyObject* PyReadAlignments_build_ralign(PyObject* cls, PyObject* args, PyObject* kwds)
{
	GKPY_TRY
	PyObject* infiles     = nullptr;
	PyObject* exclude   = nullptr;
	PyObject* allow   = nullptr;
	const char* outfile   = nullptr;
	auto include_duplicates = Py_False;
	const char* library_format{};
	PyObject* reference_genome{};

	static char* kwlist[] = {"outfile",        "infiles", "reference_genome", "exclude", "allow", "include_duplicates",
							 "library_format", nullptr};
	if (!PyArg_ParseTupleAndKeywords(args, kwds, "sOO|O!O!O!s", kwlist, &outfile, &infiles, &reference_genome,
									 &PyList_Type, &exclude, &PyList_Type, &allow, &PyBool_Type, &include_duplicates,
									 &library_format))
		return nullptr;

	const auto& genome = as_genome(reference_genome);
	read_alignments::builder raligns(outfile, genome);
	raligns.set_include_duplicates(include_duplicates == Py_True);

	if (exclude) {
		// Add excluded intervals
		for (Py_ssize_t i = 0; i < PyList_GET_SIZE(exclude); ++i) {
			PyObject* interval = PyList_GET_ITEM(exclude, i); // borrowed reference
			GK_CHECK2(PyInterval::check(interval), type, "Each exclude item must be an Interval");
			raligns.get_interval_filter().exclude(PyInterval::value(interval));
		}
	}

	if (allow) {
		// Add allowed intervals
		for (Py_ssize_t i = 0; i < PyList_GET_SIZE(allow); ++i) {
			PyObject* interval = PyList_GET_ITEM(allow, i); // borrowed reference
			GK_CHECK2(PyInterval::check(interval), type, "Each allow item must be an Interval");
			raligns.get_interval_filter().allow(PyInterval::value(interval));
		}
	}

	if (library_format) {
		raligns.detect_strand_with_library(library_format);
	}

	if (PyList_Check(infiles)) {
		// Add reads from individual files
		for (Py_ssize_t i = 0; i < PyList_GET_SIZE(infiles); ++i) {
			PyObject* infile = PyList_GET_ITEM(infiles, i); // borrowed reference
			GK_CHECK2(PyString_Check(infile), type, "Each item in the files list must be a string");
			raligns.add(PyString_AS_STRING(infile));
		}
	} else if (PyObject* fileno = PyObject_CallMethod(infiles, "fileno", nullptr); fileno) {
		// Add reads from standard input
		GKPY_TAKEREF(fileno);
		GK_CHECK2(PyInt_AsLong(fileno) == 0, value, "When infiles is a file, expected sys.stdin (i.e. fileno() == 0)");
		raligns.add(stdin_path);
	} else {
		GK_THROW2(type, "Expected infiles to be list of file names or sys.stdin");
	}

	// Build the final file.
	raligns.finalize();

	GKPY_RETURN_NONE;
	GKPY_CATCH_RETURN_NULL
}

static PyObject* PyReadAlignments_ralign_version(PyObject* cls, PyObject* args)
{
	if (!PyArg_UnpackTuple(args, "", 0, 0))
		return nullptr;
	return PyInt_FromLong((long)read_alignments::ralign_version());
}

GKPY_METHOD_BEGIN_NOARG(ReadAlignments, close)
	self->raligns.close();
	GKPY_RETURN_NONE;
GKPY_METHOD_END

GKPY_MEMBERS_BEGIN(ReadAlignments)
	GKPY_MEMBER_OBJECT(ReadAlignments, junctions,      "junctions",    nullptr)
	GKPY_MEMBER_OBJECT(ReadAlignments, alignments,     "alignments",   nullptr)
	GKPY_MEMBER_OBJECT(ReadAlignments, matches,        "matches",      nullptr)
	GKPY_MEMBER_OBJECT(ReadAlignments, variants,       "variants",     nullptr)
GKPY_MEMBERS_END

GKPY_METHODS_BEGIN(ReadAlignments)
	GKPY_METHOD_ENTRY(ReadAlignments, build_ralign,     METH_STATIC | METH_VARARGS | METH_KEYWORDS, nullptr)
	GKPY_METHOD_ENTRY(ReadAlignments, ralign_version,   METH_STATIC | METH_VARARGS, nullptr)
	GKPY_METHOD_ENTRY(ReadAlignments, close,            METH_NOARGS, nullptr)
GKPY_METHODS_END

GKPY_TYPEOBJ_BEGIN(ReadAlignments)
	tp_flags |= Py_TPFLAGS_BASETYPE | Py_TPFLAGS_HAVE_GC;
	tp_new = PyReadAlignments_New;
	tp_init = PyReadAlignments_Init;
	tp_dealloc = PyReadAlignments_Dealloc;
	tp_getattro = PyReadAlignments_GetAttro;
	tp_setattro = PyReadAlignments_SetAttro;
	tp_traverse = PyReadAlignments_Traverse;
	tp_clear = PyReadAlignments_Clear;
	tp_members = PyReadAlignments_Members;
	tp_methods = PyReadAlignments_Methods;
GKPY_TYPEOBJ_END

END_NAMESPACE_GK
