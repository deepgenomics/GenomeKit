/*
Copyright (C) 2016-2023 Deep Genomics Inc. All Rights Reserved.
*/
#include "py_util.h" // Include first to avoid POSIX warnings due to Python header inclusion order
#include "py_jrdist.h"
#include <structmember.h>

BEGIN_NAMESPACE_GK

void validate_JRDistTable(const PyAsPtrSource* self)
{
	GK_CHECK(((const PyJRDistTable*)self)->table->valid(), file,
			 "JRDist have been invalidated by ReadDistributions.close or with statement.");
}

#define GKPY_READ_DIST_TABLE_BEGIN(name, table_name) \
	GKPY_NEW_BEGIN(name##Table) \
		PyReadDistributions* owner = 0; \
		if (!PyArg_ParseTuple(args, "O!", PyReadDistributions::DefaultType, &owner)) \
			return NULL; \
		self->table = &owner->rdists.table_name(); \
		self->set_validator(&validate_##name##Table); \
	GKPY_NEW_END \
\
	GKPY_DEALLOC_BEGIN(name##Table) \
	GKPY_DEALLOC_END \
\
	GKPY_SUBTYPEOBJ_BEGIN(name##Table) \
		tp_flags |= Py_TPFLAGS_BASETYPE; /* Allow Py_TPFLAGS_HAVE_GC to be inherited so that tp_traverse/tp_clear get inherited too */ \
		tp_new = Py##name##Table_New; \
		tp_dealloc = Py##name##Table_Dealloc; \

#define GKPY_READ_DIST_TABLE_END GKPY_SUBTYPEOBJ_END

GKPY_NEW_BEGIN(ReadDistributions)
	// We expect to be called with __new__(infile) infile is the
	// path to the .rdist file we should open
	const char* infile    = nullptr;
	static char* kwlist[] = { "infile", nullptr };
	if (!PyArg_ParseTupleAndKeywords(args, kwds, "s", kwlist, &infile))
		return nullptr;

	// Construct the genome_t instance inplace
	GK_TENTATIVE_INPLACE_CONSTRUCT(rdists, read_distributions);
	self->rdists.set_source(infile);
	self->rdists.ensure_open();
	// <--- INSERT SUB-OBJECT CONSTRUCTIONS HERE

	// Construct sub-objects
	GKPY_TENTATIVE_CONSTRUCT(juncs, JRDistTable, "(O)", self);

	// Allow sub-objects to stay constructed on return
	GK_FINALIZE_CONSTRUCT(rdists);
	GK_FINALIZE_CONSTRUCT(juncs);
	// <--- INSERT SUB-OBJECT FINALIZATIONS HERE, IF CONSTRUCTED ABOVE

	self->rdists_constructed = true;
GKPY_NEW_END

GKPY_INIT_BEGIN(ReadDistributions)
GKPY_INIT_END

GKPY_TRAVERSE_BEGIN(ReadDistributions)
	GKPY_VISIT(juncs);
	// <--- INSERT SUB-OBJECT VISITS HERE
GKPY_TRAVERSE_END

GKPY_CLEAR_BEGIN(ReadDistributions)
	GKPY_CLEAR(juncs);
	// <--- INSERT SUB-OBJECT CLEARS HERE
GKPY_CLEAR_END

GKPY_DEALLOC_BEGIN(ReadDistributions)
	if (self->rdists_constructed)
		destruct(&self->rdists);
GKPY_DEALLOC_END

GKPY_GETATTRO_BEGIN(ReadDistributions)
	GKPY_GETATTR_CASE("juncs")    { GKPY_RETURN_INCREF(self->juncs); }
	GKPY_GETATTR_CASE("filename") { return PyString_FromSV(self->rdists.source()); }
	return PyObject_GenericGetAttr(selfo, attro);
GKPY_GETATTRO_END

GKPY_SETATTRO_BEGIN(ReadDistributions)
	GKPY_SETATTR_READONLY("juncs");
	GKPY_SETATTR_READONLY("filename");
	return PyObject_GenericSetAttr(selfo, attro, value);
GKPY_SETATTRO_END

static PyObject* PyReadDistributions_build_rdist(PyObject* cls, PyObject* args, PyObject* kwds)
{
	GKPY_TRY
	const char* outfile = nullptr;
	PyObject* infiles   = nullptr;
	int min_reads       = 0;
	int min_overhang    = 0;
	PyObject* exclude = nullptr;
	PyObject* allow = nullptr;
	PyObject* filter    = Py_None;
	static char* kwlist[]
		= { "outfile", "infiles", "min_reads", "min_overhang", "exclude", "allow", "pass_filter", nullptr };
	if (!PyArg_ParseTupleAndKeywords(args, kwds, "sO!|iiO!O!O", kwlist, &outfile, &PyList_Type, &infiles, &min_reads,
									 &min_overhang, &PyList_Type, &exclude, &PyList_Type, &allow, &filter))
		return nullptr;

	read_distributions::builder rdists(outfile);
	rdists.set_min_reads(min_reads);
	rdists.set_min_overhang(min_overhang);

	if (exclude) {
		// Add excluded intervals
		for (Py_ssize_t i = 0; i < PyList_GET_SIZE(exclude); ++i) {
			PyObject* interval = PyList_GET_ITEM(exclude, i); // borrowed reference
			GK_CHECK(PyInterval::check(interval), type, "Each exclude item must be an Interval");
			rdists.get_interval_filter().exclude(PyInterval::value(interval));
		}
	}

	if (allow) {
		// Add allowed intervals
		for (Py_ssize_t i = 0; i < PyList_GET_SIZE(allow); ++i) {
			PyObject* interval = PyList_GET_ITEM(allow, i); // borrowed reference
			GK_CHECK(PyInterval::check(interval), type, "Each allow item must be an Interval");
			rdists.get_interval_filter().allow(PyInterval::value(interval));
		}
	}

	if (filter != Py_None) {
		GK_CHECK(PyCallable_Check(filter), type, "filter must be a callable object that returns a bool");
	}

	// Add read alignments from individual files
	for (Py_ssize_t i = 0; i < PyList_GET_SIZE(infiles); ++i) {
		PyObject* infile = PyList_GET_ITEM(infiles, i); // borrowed reference
		GK_CHECK(PyString_Check(infile), type, "Each item in the files list must be a string");

		if (filter != Py_None) {
			rdists.filter([filter, infile](const junction_read_alignments& ralign, int junction, unsigned int read) {
				PyAutoRef args{ Py_BuildValue("(Oii)", infile, junction, read) };
				PyAutoRef pass{ PyObject_CallObject(filter, args.get()) };
				if (!pass.get()) {
					PyErr_Print();
					GK_THROW(value, "filter could not be called, see error above");
				}
				GK_CHECK(PyBool_Check(pass.get()), type, "filter must be a callable object that returns a bool");
				return pass.get() == Py_True;
			});
		}
		rdists.add(PyString_AS_STRING(infile));
	}

	// Build the final file.
	rdists.finalize();

	GKPY_RETURN_NONE;
	GKPY_CATCH_RETURN_NULL
}

static PyObject* PyReadDistributions_rdist_version(PyObject* cls, PyObject* args)
{
	if (!PyArg_UnpackTuple(args, "", 0, 0))
		return nullptr;
	return PyInt_FromLong((long)read_distributions::rdist_version());
}

GKPY_METHOD_BEGIN_NOARG(ReadDistributions, close)
	self->rdists.close();
	GKPY_RETURN_NONE;
GKPY_METHOD_END

GKPY_MEMBERS_BEGIN(ReadDistributions)
	GKPY_MEMBER_OBJECT(ReadDistributions, juncs, "junctions", nullptr)
GKPY_MEMBERS_END

GKPY_METHODS_BEGIN(ReadDistributions)
	GKPY_METHOD_ENTRY(ReadDistributions, build_rdist,   METH_STATIC | METH_VARARGS | METH_KEYWORDS, nullptr)
	GKPY_METHOD_ENTRY(ReadDistributions, rdist_version, METH_STATIC | METH_VARARGS, nullptr)
	GKPY_METHOD_ENTRY(ReadDistributions, close,         METH_NOARGS, nullptr)
GKPY_METHODS_END

GKPY_TYPEOBJ_BEGIN(ReadDistributions)
	tp_flags |= Py_TPFLAGS_BASETYPE | Py_TPFLAGS_HAVE_GC;
	tp_new = PyReadDistributions_New;
	tp_init = PyReadDistributions_Init;
	tp_dealloc = PyReadDistributions_Dealloc;
	tp_getattro = PyReadDistributions_GetAttro;
	tp_setattro = PyReadDistributions_SetAttro;
	tp_traverse = PyReadDistributions_Traverse;
	tp_clear = PyReadDistributions_Clear;
	tp_members = PyReadDistributions_Members;
	tp_methods = PyReadDistributions_Methods;
GKPY_TYPEOBJ_END

////////////////////////////////////////////////////////////
// JRCount, JRDist, JRDistTable
////////////////////////////////////////////////////////////

GKPY_GETATTRO_BEGIN(JRCount)
	// ******** UPDATE genome_kit.JunctionReadCount MOCK ATTRIBUTES TO MATCH ANY CHANGES HERE ********
	const jrcount_t& c = self->value();
	GKPY_GETATTR_CASE("strand")      { return PyString_FromStrand(c.strand); }
	GKPY_GETATTR_CASE("shift")       { return PyInt_FromLong(c.shift); }
	GKPY_GETATTR_CASE("count")       { return PyInt_FromLong(c.count); }
	return PyObject_GenericGetAttr(selfo, attro);
GKPY_GETATTRO_END

GKPY_SETATTRO_BEGIN(JRCount)
	GKPY_SETATTR_READONLY("strand");
	GKPY_SETATTR_READONLY("shift");
	GKPY_SETATTR_READONLY("count");
	return PyObject_GenericSetAttr(selfo, attro, value);
GKPY_SETATTRO_END

GKPY_VALUE_TYPEOBJ_BEGIN(JRCount, as_value) // as_value because we NEVER point to an entry in a jrdist_t instance, jrcount_t are created on the fly
	tp_flags |= Py_TPFLAGS_BASETYPE;
	tp_new = PyFastNew<PyJRCount>;
	tp_dealloc = PyFastDealloc<PyJRCount>;
	tp_getattro = PyJRCount_GetAttro;
	tp_setattro = PyJRCount_SetAttro;
	tp_richcompare = PyGenericValue_RichCompare<PyJRCount>; // for hash == comparison
GKPY_VALUE_TYPEOBJ_END

///////////////////////////////////////////////////////

GKPY_GETATTRO_BEGIN(JRDist)
	// ******** UPDATE genome_kit.JunctionReadDistribution MOCK ATTRIBUTES TO MATCH ANY CHANGES HERE ********
	GKPY_GETATTR_CASE("interval")    { return PyInterval::create(self->value()); }
	GKPY_GETATTR_CASE("num_reads")   { return PyInt_FromLong(self->value().num_reads); }
	GKPY_GETATTR_CASE("num_counts")  { return PyInt_FromLong(self->value().num_counts); }
	return PyJRDist::base_type::DefaultType->tp_getattro(selfo, attro);
GKPY_GETATTRO_END

GKPY_SETATTRO_BEGIN(JRDist)
	GKPY_SETATTR_READONLY("interval");
	GKPY_SETATTR_READONLY("num_reads");
	GKPY_SETATTR_READONLY("num_counts");
	return PyJRDist::base_type::DefaultType->tp_setattro(selfo, attro, value);
GKPY_SETATTRO_END

Py_ssize_t PyJRDist_Len(PyObject* self) { return (Py_ssize_t)PyJRDist::value(self).num_counts; }

PyObject* PyJRDist_GetItem(PyObject* self, Py_ssize_t index)
{
	GKPY_TRY
	const jrdist_t jrdist = ((PyJRDist*)self)->unpack();
	GKPY_INDEXCHECK(index, (Py_ssize_t)jrdist.num_counts());

	// Create a new instance of T::DefaultType and initialize it
	// to point to the element in the table.
	auto* jrcount    = (PyJRCount*)PyJRCount::DefaultType->tp_new(PyJRCount::DefaultType, nullptr, nullptr);
	jrcount->as_ptr  = nullptr;
	jrcount->value() = jrdist[(index_t)index];
	return (PyObject*)jrcount;
	GKPY_CATCH_RETURN_NULL
}

GKPY_VALUE_SUBTYPEOBJ_BEGIN(JRDist, as_ptr) // as_ptr because we always point to a packed_jrdist entry in a jrdist_table instance
	tp_flags |= Py_TPFLAGS_BASETYPE;
	tp_init = gkpy_invalid_init;
	tp_iter = PySeqIter_New;    // Standard sequence iterator calls sq_item until IndexError is raised
	sq_item = PyJRDist_GetItem;
	sq_length = PyJRDist_Len;
	tp_getattro = PyJRDist_GetAttro;
	tp_setattro = PyJRDist_SetAttro;
	tp_richcompare = PyGenericValue_RichCompare<PyJRDist>; // for hash == comparison
	tp_hash = PyGenericValue_Hash<PyJRDist>;  // uses hash(const interval_t&), which is ok since each JRDist is for a different interval
GKPY_VALUE_SUBTYPEOBJ_END

GKPY_READ_DIST_TABLE_BEGIN(JRDist, juncs)
GKPY_READ_DIST_TABLE_END

END_NAMESPACE_GK
