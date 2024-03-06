/*
Copyright (C) 2016-2023 Deep Genomics Inc. All Rights Reserved.
*/
#include "py_util.h" // Include first to avoid POSIX warnings due to Python header inclusion order

#include "py_jralign.h"
#include "py_genome.h"
#include "py_variant_table.h"
#include "util.h"
#include <structmember.h>

BEGIN_NAMESPACE_GK

void validate_JRAlignsTable(const PyAsPtrSource* self)
{
	GK_CHECK(((const PyJRAlignsTable*)self)->table->valid(), file,
			 "JRAligns have been invalidated by JReadAlignments.close or with statement.");
}

#define GKPY_JREAD_ALIGNS_TABLE_BEGIN(name, table_name) \
	GKPY_NEW_OWNED_BEGIN(name##Table, JReadAlignments) \
		self->table = &((PyJReadAlignments*)self->owner)->jraligns.table_name(); \
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

#define GKPY_JREAD_ALIGNS_TABLE_END GKPY_SUBTYPEOBJ_END

GKPY_NEW_BEGIN(JReadAlignments)
	// We expect to be called with __new__(infile) infile is the
	// path to the .jralign file we should open
	const char* infile    = nullptr;
	static char* kwlist[] = { "infile", nullptr };
	if (!PyArg_ParseTupleAndKeywords(args, kwds, "s", kwlist, &infile))
		return nullptr;

	// Construct the genome_t instance inplace
	GK_TENTATIVE_INPLACE_CONSTRUCT(jraligns, junction_read_alignments);
	self->jraligns.set_source(infile);
	self->jraligns.ensure_open();
	// <--- INSERT SUB-OBJECT CONSTRUCTIONS HERE

	// Construct sub-objects
	GKPY_TENTATIVE_CONSTRUCT(juncs, JRAlignsTable, "(O)", self);
	GKPY_TENTATIVE_CONSTRUCT(variants, VariantTable, "(On)", self, (Py_ssize_t)&self->jraligns.variants());

	// Allow sub-objects to stay constructed on return
	GK_FINALIZE_CONSTRUCT(jraligns);
	GK_FINALIZE_CONSTRUCT(juncs);
	GK_FINALIZE_CONSTRUCT(variants);
	// <--- INSERT SUB-OBJECT FINALIZATIONS HERE, IF CONSTRUCTED ABOVE

	self->jraligns_constructed = true;
GKPY_NEW_END

GKPY_INIT_BEGIN(JReadAlignments)
GKPY_INIT_END

GKPY_TRAVERSE_BEGIN(JReadAlignments)
	GKPY_VISIT(juncs);
	GKPY_VISIT(variants);
	// <--- INSERT SUB-OBJECT VISITS HERE
GKPY_TRAVERSE_END

GKPY_CLEAR_BEGIN(JReadAlignments)
	GKPY_CLEAR(juncs);
	GKPY_CLEAR(variants);
	// <--- INSERT SUB-OBJECT CLEARS HERE
GKPY_CLEAR_END

GKPY_DEALLOC_BEGIN(JReadAlignments)
	if (self->jraligns_constructed)
		destruct(&self->jraligns);
GKPY_DEALLOC_END

GKPY_GETATTRO_BEGIN(JReadAlignments)
	GKPY_GETATTR_CASE("juncs")    { GKPY_RETURN_INCREF(self->juncs); }
	GKPY_GETATTR_CASE("variants")    { GKPY_RETURN_INCREF(self->variants); }
	GKPY_GETATTR_CASE("filename") { return PyString_FromSV(self->jraligns.source()); }
	return PyObject_GenericGetAttr(selfo, attro);
GKPY_GETATTRO_END

GKPY_SETATTRO_BEGIN(JReadAlignments)
	GKPY_SETATTR_READONLY("juncs");
	GKPY_SETATTR_READONLY("variants");
	GKPY_SETATTR_READONLY("filename");
	return PyObject_GenericSetAttr(selfo, attro, value);
GKPY_SETATTRO_END

static PyObject* PyJReadAlignments_build_jralign(PyObject* cls, PyObject* args, PyObject* kwds)
{
	GKPY_TRY
	PyObject* infiles    = nullptr;
	PyObject* exclude  = nullptr;
	PyObject* allow  = nullptr;
	int min_reads = 0;
	int min_overhang = 0;
	auto include_variants = Py_False;
	auto include_duplicates = Py_False;
	const char* outfile = nullptr;
	const char* library_format{};
	PyObject* reference_genome{};
	const char* overhang_error = "error";

	static char* kwlist[] = {"outfile",        "infiles", "reference_genome", "min_reads",          "min_overhang",
							 "exclude",        "allow",   "include_variants", "include_duplicates", "library_format",
							 "overhang_error", nullptr};
	if (!PyArg_ParseTupleAndKeywords(args, kwds, "sOO|iiO!O!O!O!ss", kwlist, &outfile, &infiles, &reference_genome,
									 &min_reads, &min_overhang, &PyList_Type, &exclude, &PyList_Type, &allow,
									 &PyBool_Type, &include_variants, &PyBool_Type, &include_duplicates,
									 &library_format, &overhang_error))
		return nullptr;

	const auto&                       genome = as_genome(reference_genome);
	junction_read_alignments::builder jraligns(outfile, genome);
	jraligns.set_min_reads(min_reads);
	jraligns.set_min_overhang(min_overhang);
	jraligns.set_include_variants(include_variants == Py_True);
	jraligns.set_include_duplicates(include_duplicates == Py_True);

	if (exclude) {
		// Add excluded intervals
		for (Py_ssize_t i = 0; i < PyList_GET_SIZE(exclude); ++i) {
			PyObject* interval = PyList_GET_ITEM(exclude, i); // borrowed reference
			GK_CHECK(PyInterval::check(interval), type, "Each exclude item must be an Interval");
			jraligns.exclude(PyInterval::value(interval));
		}
	}

	if (allow) {
		// Add allowed intervals
		for (Py_ssize_t i = 0; i < PyList_GET_SIZE(allow); ++i) {
			PyObject* interval = PyList_GET_ITEM(allow, i); // borrowed reference
			GK_CHECK(PyInterval::check(interval), type, "Each allow item must be an Interval");
			jraligns.allow(PyInterval::value(interval));
		}
	}

	if (library_format) {
		jraligns.detect_strand_with_library(library_format);
	}

	if (strcmp(overhang_error, "clamp") == 0) {
		jraligns.set_overhang_error(junction_read_alignments::builder::error_handling::clamp);
	} else {
		GK_CHECK(strcmp(overhang_error, "error") == 0, value, "error_handling must be in [\"error\", \"clamp\"]");
	}

	if (PyList_Check(infiles)) {
		// Add reads from individual files
		for (Py_ssize_t i = 0; i < PyList_GET_SIZE(infiles); ++i) {
			PyObject* infile = PyList_GET_ITEM(infiles, i); // borrowed reference
			GK_CHECK(PyString_Check(infile), type, "Each item in the files list must be a string");
			jraligns.add(PyString_AS_STRING(infile));
		}
	} else if (PyObject* fileno = PyObject_CallMethod(infiles, "fileno", nullptr); fileno) {
		// Add reads from standard input
		GKPY_TAKEREF(fileno);
		GK_CHECK(PyInt_AsLong(fileno) == 0, value, "When infiles is a file, expected sys.stdin (i.e. fileno() == 0)");
		jraligns.add(stdin_path);
	} else {
		GK_THROW(type, "Expected infiles to be list of file names or sys.stdin");
	}

	// Build the final file.
	jraligns.finalize();

	GKPY_RETURN_NONE;
	GKPY_CATCH_RETURN_NULL
}

static PyObject* PyJReadAlignments_jralign_version(PyObject* cls, PyObject* args)
{
	if (!PyArg_UnpackTuple(args, "", 0, 0))
		return nullptr;
	return PyInt_FromLong((long)junction_read_alignments::jralign_version());
}

GKPY_METHOD_BEGIN_NOARG(JReadAlignments, close)
	self->jraligns.close();
	GKPY_RETURN_NONE;
GKPY_METHOD_END

GKPY_MEMBERS_BEGIN(JReadAlignments)
	GKPY_MEMBER_OBJECT(JReadAlignments, juncs, "junctions", nullptr)
	GKPY_MEMBER_OBJECT(JReadAlignments, variants, "variants", nullptr)
GKPY_MEMBERS_END

GKPY_METHODS_BEGIN(JReadAlignments)
	GKPY_METHOD_ENTRY(JReadAlignments, build_jralign,   METH_STATIC | METH_VARARGS | METH_KEYWORDS, nullptr)
	GKPY_METHOD_ENTRY(JReadAlignments, jralign_version, METH_STATIC | METH_VARARGS, nullptr)
	GKPY_METHOD_ENTRY(JReadAlignments, close,           METH_NOARGS, nullptr)
GKPY_METHODS_END

GKPY_TYPEOBJ_BEGIN(JReadAlignments)
	tp_flags |= Py_TPFLAGS_BASETYPE | Py_TPFLAGS_HAVE_GC;
	tp_new = PyJReadAlignments_New;
	tp_init = PyJReadAlignments_Init;
	tp_dealloc = PyJReadAlignments_Dealloc;
	tp_getattro = PyJReadAlignments_GetAttro;
	tp_setattro = PyJReadAlignments_SetAttro;
	tp_traverse = PyJReadAlignments_Traverse;
	tp_clear = PyJReadAlignments_Clear;
	tp_members = PyJReadAlignments_Members;
	tp_methods = PyJReadAlignments_Methods;
GKPY_TYPEOBJ_END

////////////////////////////////////////////////////////////
// JRAlign, JRAligns, JRAlignsTable
////////////////////////////////////////////////////////////

GKPY_NEW_OWNED_BEGIN(JRAlign, JRAligns)
GKPY_NEW_OWNED_END

GKPY_DEALLOC_OWNED_BEGIN(JRAlign)
GKPY_DEALLOC_OWNED_END

GKPY_GETATTRO_BEGIN(JRAlign)
	// ******** UPDATE genome_kit.JunctionReadCount MOCK ATTRIBUTES TO MATCH ANY CHANGES HERE ********
	const jralign_t& c = self->value();
	GKPY_GETATTR_CASE("left")          { return PyInt_FromLong(c.left); }
	GKPY_GETATTR_CASE("right")         { return PyInt_FromLong(c.right); }
	GKPY_GETATTR_CASE("strand")        { return PyString_FromStrand(c.strand); }
	GKPY_GETATTR_CASE("num_variants")  { return PyInt_FromLong(c.num_variants); }
	return PyObject_GenericGetAttr(selfo, attro);
GKPY_GETATTRO_END

GKPY_SETATTRO_BEGIN(JRAlign)
	GKPY_SETATTR_READONLY("left");
	GKPY_SETATTR_READONLY("right");
	GKPY_SETATTR_READONLY("strand");
	GKPY_SETATTR_READONLY("num_variants");
	return PyObject_GenericSetAttr(selfo, attro, value);
GKPY_SETATTRO_END

GKPY_OMETHOD_BEGIN(JRAlign, variants)

	const jralign_t& c = self->value();

	const PyJReadAlignments& jread_alignments = ((PyJRAligns*)self->owner)->junction_read_alignments();
	PyVariantTable::table_type& table = PyVariantTable::value(jread_alignments.variants);

	// Constructing tuple of PyVariants
	PyObject* tuple = PyTuple_New(c.num_variants);
	GKPY_TAKEREF(tuple);
	for (int i = 0; i < c.num_variants; ++i)
		PyTuple_SET_ITEM(tuple, i, PyTable_CreateItem<PyVariant>(jread_alignments.variants, &table[c.variants_indices[i]]));
	GKPY_FORGETREF(tuple); // keep list refcount at 1
	return tuple;

GKPY_OMETHOD_END

GKPY_METHODS_BEGIN(JRAlign)
	GKPY_METHOD_ENTRY(JRAlign, variants, METH_VARARGS | METH_KEYWORDS, nullptr)
GKPY_METHODS_END

GKPY_VALUE_TYPEOBJ_BEGIN(JRAlign, as_value) // as_value because we NEVER point to an entry in a jrdist_t instance, jralign_t are created on the fly
	tp_basicsize += sizeof(void*);
	tp_flags |= Py_TPFLAGS_BASETYPE;
	tp_new = PyJRAlign_New;
	tp_dealloc = PyJRAlign_Dealloc;
	tp_getattro = PyJRAlign_GetAttro;
	tp_setattro = PyJRAlign_SetAttro;
	tp_richcompare = PyGenericValue_RichCompare<PyJRAlign>; // for hash == comparison
	tp_methods = PyJRAlign_Methods;
GKPY_VALUE_TYPEOBJ_END

///////////////////////////////////////////////////////

GKPY_GETATTRO_BEGIN(JRAligns)
	// ******** UPDATE genome_kit.JunctionReadDistribution MOCK ATTRIBUTES TO MATCH ANY CHANGES HERE ********
	GKPY_GETATTR_CASE("interval")    { return PyInterval::create(self->value()); }
	GKPY_GETATTR_CASE("num_reads")   { return PyInt_FromLong(self->value().num_reads); }
	return PyJRAligns::base_type::DefaultType->tp_getattro(selfo, attro);
GKPY_GETATTRO_END

GKPY_SETATTRO_BEGIN(JRAligns)
	GKPY_SETATTR_READONLY("interval");
	GKPY_SETATTR_READONLY("num_reads");
	return PyJRAligns::base_type::DefaultType->tp_setattro(selfo, attro, value);
GKPY_SETATTRO_END

Py_ssize_t PyJRAligns_Len(PyObject* self) { return (Py_ssize_t)PyJRAligns::value(self).num_reads; }

PyObject* PyJRAligns_GetItem(PyObject* self, Py_ssize_t index)
{
	GKPY_TRY
	const jraligns_t jraligns = ((PyJRAligns*)self)->unpack();
	GKPY_INDEXCHECK(index, (Py_ssize_t)jraligns.num_reads());

	// Create a new instance of T::DefaultType and initialize it
	// to point to the element in the table.
	GKPY_BUILDVALUE_TEMP(args, "(O)", self);
	auto* jrcount    = (PyJRAlign*)PyJRAlign::DefaultType->tp_new(PyJRAlign::DefaultType, args, nullptr);
	jrcount->as_ptr  = nullptr;
	jrcount->value() = jraligns[(index_t)index];

	return (PyObject*)jrcount;
	GKPY_CATCH_RETURN_NULL
}

GKPY_VALUE_SUBTYPEOBJ_BEGIN(JRAligns, as_ptr) // as_ptr because we always point to a packed_jrdist entry in a jrdist_table instance
	tp_flags |= Py_TPFLAGS_BASETYPE;
	tp_init = gkpy_invalid_init;
	tp_iter = PySeqIter_New;    // Standard sequence iterator calls sq_item until IndexError is raised
	sq_item = PyJRAligns_GetItem;
	sq_length = PyJRAligns_Len;
	tp_getattro = PyJRAligns_GetAttro;
	tp_setattro = PyJRAligns_SetAttro;
	tp_richcompare = PyGenericValue_RichCompare<PyJRAligns>; // for hash == comparison
	tp_hash = PyGenericValue_Hash<PyJRAligns>;  // uses hash(const interval_t&), which is ok since each JRAligns is for a different interval
GKPY_VALUE_SUBTYPEOBJ_END

GKPY_JREAD_ALIGNS_TABLE_BEGIN(JRAligns, juncs)
GKPY_JREAD_ALIGNS_TABLE_END

END_NAMESPACE_GK
