/*
Copyright (C) 2016-2023 Deep Genomics Inc. All Rights Reserved.
*/
#include "py_variant_table.h"
#include "py_util.h" // Include first to avoid POSIX warnings due to Python header inclusion order

#include "genome.h"
#include <numeric>
#include <numpy/arrayobject.h>

BEGIN_NAMESPACE_GK

using namespace std;

using vt = vcf_table;

static constexpr int py_dtypes[vt::num_dtype] = {
	NPY_INT8,    // int8
	NPY_INT16,   // int16
	NPY_INT32,   // int32
	NPY_FLOAT16, // float16
	NPY_FLOAT32, // float32
	NPY_NOTYPE,  // str
	NPY_BOOL,    // bool
};

////////////////////////////////////////////////////////////
// avariant_t
////////////////////////////////////////////////////////////

bool operator==(const avariant_t& lhs, const avariant_t& rhs)
{
	return lhs.as_interval() == rhs.as_interval() && !strcmp(lhs.ref, rhs.ref) && !strcmp(lhs.alt, rhs.alt);
}

bool operator!=(const avariant_t& lhs, const avariant_t& rhs) { return !(lhs == rhs); }

////////////////////////////////////////////////////////////
// PyVariant
////////////////////////////////////////////////////////////

avariant_t::avariant_t(chrom_t chrom, pos_t start, const char* ref, size_t reflen, const char* alt, size_t altlen,
					   refg_t refg)
	: ainterval_t{ { chrom, start, start + int_cast<pos_t>(reflen) - 1, pos_strand, refg } }
{
	this->anchor        = 0; // Disabled on variants; recycled for ref/alt small-string optimization
	this->anchor_offset = 0; // Disabled on variants;

	// small-string optimization for ref/alt:
	//  - if <= 8 bytes, stored in 8-bytes of [anchor,anchor_offset]
	//  - if >  8 bytes, ref and alt both stored in single new[] pointed to by ref
	this->ref = reflen + altlen + 2 <= sizeof(anchor) + sizeof(anchor_offset)
					? (char*)&anchor
					: new char[reflen + altlen + 2]; // lump ref and alt into one malloc
	this->alt = this->ref + reflen + 1;
	memcpy(this->ref, ref, reflen + 1);
	memcpy(this->alt, alt, altlen + 1);
}

GKPY_IMETHOD_BEGIN(Variant, Init)

	char* arg_names[] = { "chromosome", "start", "ref", "alt", "reference_genome", nullptr };
	const char* chrom;
	int start;
	const char* ref;
	const char* alt;
	PyObject* pyrefg;
	if (!PyArg_ParseTupleAndKeywords(args, kwds, "sissO", arg_names, &chrom, &start, &ref, &alt, &pyrefg))
		return -1;

	const auto refg          = as_refg(pyrefg);
	PyObject*  data_dir_attr = PyObject_GetAttrString(pyrefg, "data_dir");
	PyErr_Clear();
	GKPY_TAKEREF(data_dir_attr);
	const auto& chrom_names =
		data_dir_attr ? get_chrom_names(refg, PyString_AsString(data_dir_attr)) : get_chrom_names(refg);
	// Construct the avariant_t inplace
	new (&self->pyconstructed_value())
		avariant_t(chrom_names.as_chrom(chrom), start, ref, strlen(ref), alt, strlen(alt), refg);

	return 0;
GKPY_IMETHOD_END

PyObject* PyVariant_RichCompare(PyObject* ao, PyObject* bo, int op)
{
	GKPY_TRY
	if (Py_TYPE(ao) != Py_TYPE(bo)) {
		if (op == Py_EQ)
			GKPY_RETURN_FALSE;
		if (op == Py_NE)
			GKPY_RETURN_TRUE;
		GK_THROW(type, "Incompatible arguments '{}' and '{}'.", Py_TYPE(ao)->tp_name, Py_TYPE(bo)->tp_name);
	}
	variant_t a = ((PyVariant*)ao)->as_variant();
	variant_t b = ((PyVariant*)bo)->as_variant();
	switch (op) {
	case Py_EQ: GKPY_RETURN_BOOL(a == b);
	case Py_NE: GKPY_RETURN_BOOL(a != b);
	}
	GK_THROW(type, "Unsupported comparison operator.");
	GKPY_CATCH_RETURN_NULL
}

GKPY_GETATTRO_BEGIN(Variant)
	// ******** UPDATE genome_kit.Variant MOCK ATTRIBUTES TO MATCH ANY CHANGES HERE ********
	GKPY_GETATTR_CASE("position")
	{
		return PyInt_FromLong(
			(long)(self->is_pyconstructed() ? self->pyconstructed_value().start() : self->unpack().start()));
	}
	GKPY_GETATTR_CASE("ref")
	{
		return PyString_FromString(self->is_pyconstructed() ? self->pyconstructed_value().ref : self->unpack().ref);
	}
	GKPY_GETATTR_CASE("alt")
	{
		return PyString_FromString(self->is_pyconstructed() ? self->pyconstructed_value().alt : self->unpack().alt);
	}
	// remove anchors
	GKPY_GETATTR_CASE("end5") return PyInterval::create(get_end5(self->value()));
	GKPY_GETATTR_CASE("end3") return PyInterval::create(get_end3(self->value()));

	GKPY_GETATTR_CASE("interval") { return PyInterval::create(self->value()); }
	GKPY_GETATTR_CASE("anchor") { GK_THROW(type, "Anchor attribute disabled on variants"); }
	GKPY_GETATTR_CASE("anchor_offset") { GK_THROW(type, "Anchor attribute disabled on variants"); }
	return PyVariant::base_type::DefaultType->tp_getattro(selfo, attro);
GKPY_GETATTRO_END

GKPY_SETATTRO_BEGIN(Variant)
	GKPY_SETATTR_READONLY("ref");
	GKPY_SETATTR_READONLY("alt");
	GKPY_SETATTR_READONLY("interval");
	return PyVariant::base_type::DefaultType->tp_setattro(selfo, attro, value);
GKPY_SETATTRO_END

GKPY_METHOD_BEGIN_ONEARG(Variant, shift)
	GKPY_TYPECHECK_BUILTIN(arg, PyInt);
	return PyInterval::create(self->value().shift((pos_t)PyInt_AS_LONG(arg)));
GKPY_METHOD_END
GKPY_METHOD_BEGIN_VARARGS(Variant, expand)
	pos_t up;
	pos_t dn;
	Py_ssize_t size = PyTuple_GET_SIZE(args);
	if (size == 1) {
		PyObject* arg = PyTuple_GET_ITEM(args, 0);
		GKPY_TYPECHECK_BUILTIN(arg, PyInt);
		up = dn = (pos_t)PyInt_AS_LONG(arg);
	} else if (size == 2) {
		PyObject* arg0 = PyTuple_GET_ITEM(args, 0);
		GKPY_TYPECHECK_BUILTIN(arg0, PyInt);
		PyObject* arg1 = PyTuple_GET_ITEM(args, 1);
		GKPY_TYPECHECK_BUILTIN(arg1, PyInt);
		up = (pos_t)PyInt_AS_LONG(arg0);
		dn = (pos_t)PyInt_AS_LONG(arg1);
	} else {
		GK_THROW(value, "Expected 1 or 2 arguments but got {}", size);
	}
	return PyInterval::create(self->value().expand(up, dn));
GKPY_METHOD_END
GKPY_METHOD_BEGIN_ONEARG(Variant, as_positive_strand)
	return PyInterval::create(PyInterval::value(selfo).as_pos_strand());
GKPY_METHOD_END
GKPY_METHOD_BEGIN_ONEARG(Variant, as_negative_strand)
	return PyInterval::create(PyInterval::value(selfo).as_neg_strand());
GKPY_METHOD_END
GKPY_METHOD_BEGIN_ONEARG(Variant, as_opposite_strand)
	return PyInterval::create(PyInterval::value(selfo).as_opp_strand());
GKPY_METHOD_END

GKPY_METHODS_BEGIN(Variant)
	GKPY_METHOD_ENTRY(Variant, shift, METH_O, nullptr)
	GKPY_METHOD_ENTRY(Variant, expand, METH_VARARGS, nullptr)
	GKPY_METHOD_ENTRY(Variant, as_positive_strand, METH_NOARGS, nullptr)
	GKPY_METHOD_ENTRY(Variant, as_negative_strand, METH_NOARGS, nullptr)
	GKPY_METHOD_ENTRY(Variant, as_opposite_strand, METH_NOARGS, nullptr)
GKPY_METHODS_END

GKPY_VALUE_SUBTYPEOBJ_BEGIN2(
	Variant,
	as_value_or_ptr) //  as_value_or_ptr to ensure there's enough tp_basicsize to hold an avariant_t value inplace
	tp_flags |= Py_TPFLAGS_BASETYPE;
	tp_init        = PyVariant_Init;
	tp_new         = PyFastNew<PyVariant>;
	tp_dealloc     = PyFastDealloc<PyVariant>;
	tp_getattro    = PyVariant_GetAttro;
	tp_setattro    = PyVariant_SetAttro;
	tp_richcompare = PyVariant_RichCompare;
	tp_hash        = PyGenericValue_Hash<PyVariant>;
	tp_methods     = PyVariant_Methods;
GKPY_VALUE_SUBTYPEOBJ_END

////////////////////////////////////////////////////////////
// VCFVariant
////////////////////////////////////////////////////////////

PyObject* PyVCFVariant_RichCompare(PyObject* ao, PyObject* bo, int op)
{
	GKPY_TRY
	if (Py_TYPE(ao) != Py_TYPE(bo)) {
		if (op == Py_EQ)
			GKPY_RETURN_FALSE;
		if (op == Py_NE)
			GKPY_RETURN_TRUE;
		GK_THROW(type, "Incompatible arguments '{}' and '{}'.", Py_TYPE(ao)->tp_name, Py_TYPE(bo)->tp_name);
	}
	variant_t a(PyVCFVariant::value(ao), *((PyVCFTable*)((PyVCFVariant*)ao)->source())->table);
	variant_t b(PyVCFVariant::value(bo), *((PyVCFTable*)((PyVCFVariant*)bo)->source())->table);
	switch (op) {
	case Py_EQ: GKPY_RETURN_BOOL(a == b);
	case Py_NE: GKPY_RETURN_BOOL(a != b);
	}
	GK_THROW(type, "Unsupported comparison operator.");
	GKPY_CATCH_RETURN_NULL
}

GKPY_GETATTRO_BEGIN(VCFVariant)
	PyObject* colattr = ((PyVCFTable*)self->source())->get_col_attr(&self->value(), attr);
	if (colattr)
		return colattr;
	return PyVCFVariant::base_type::DefaultType->tp_getattro(selfo, attro);
GKPY_GETATTRO_END

GKPY_SETATTRO_BEGIN(VCFVariant)
	// TODO: raise readonly error if user tries to set any recognized INFO/FORMAT column values
	return PyVCFVariant::base_type::DefaultType->tp_setattro(selfo, attro, value);
GKPY_SETATTRO_END

GKPY_VALUE_SUBTYPEOBJ_BEGIN(VCFVariant, as_ptr) // as_ptr because we always point to an entry in a vcf_table instance
	tp_flags |= Py_TPFLAGS_BASETYPE;
	tp_init        = gkpy_invalid_init;
	tp_getattro    = PyVCFVariant_GetAttro;
	tp_setattro    = PyVCFVariant_SetAttro;
	tp_richcompare = PyVCFVariant_RichCompare;
GKPY_VALUE_SUBTYPEOBJ_END

////////////////////////////////////////////////////////////
// VCFTable
////////////////////////////////////////////////////////////

PyAutoRef PyVCFTable::numpy_memmap_function;

void validate_VCFTable(const PyAsPtrSource* self)
{
	GK_CHECK(((const PyVCFTable*)self)->table->valid(), file,
			 "VCFVariants have been invalidated by VCFTable.close or with statement.");
}

GKPY_NEW_BEGIN(VCFTable)

	const char* infile    = nullptr;
	static char* kwlist[] = { "infile", nullptr };

	if (!PyArg_ParseTupleAndKeywords(args, kwds, "s", kwlist, &infile))
		return nullptr;

	if (!PyVCFTable::numpy_memmap_function) {
		PyAutoRef numpy_name{ PyString_FromString("numpy") };
		PyAutoRef numpy_module{ PyImport_Import(numpy_name.get()) };
		if (!numpy_module)
			return nullptr;
		PyVCFTable::numpy_memmap_function.reset(PyObject_GetAttrString(numpy_module.get(), "memmap"));
	}

	self->table = new vcf_table{mmap_file(infile)};
	self->set_validator(&validate_VCFTable);
	self->pystr_pool     = nullptr;
	self->pystr_filename = PyString_FromString(infile);

GKPY_NEW_END

// Required to depress
// DeprecationWarning: object.__init__() takes no parameters
GKPY_INIT_BEGIN(VCFTable)
GKPY_INIT_END

GKPY_DEALLOC_BEGIN(VCFTable)
	if (self->pystr_pool) {
		delete self->pystr_pool;
		self->pystr_pool = nullptr;
	}
	if (self->table) {
		delete self->table;
		self->table = nullptr;
	}
	Py_XDECREF(self->pystr_filename);
GKPY_DEALLOC_END

GKPY_GETATTRO_BEGIN(VCFTable)
	GKPY_GETATTR_CASE("info_ids")
	{
		const auto& fields = self->table->info_fields().cols;
		auto size          = fields.size();
		auto list          = PyList_New(size);
		for (size_t i = 0; i < size; ++i) { PyList_SET_ITEM(list, i, PyString_FromString(fields[i].id)); }
		return list;
	}
	GKPY_GETATTR_CASE("format_ids")
	{
		const auto& fields = self->table->fmt_fields().cols;
		auto size          = fields.size();
		auto list          = PyList_New(size);
		for (size_t i = 0; i < size; ++i) { PyList_SET_ITEM(list, i, PyString_FromString(fields[i].id)); }
		return list;
	}
	GKPY_GETATTR_CASE("filename")
	{
		Py_INCREF(self->pystr_filename);
		return self->pystr_filename;
	}
	GKPY_GETATTR_CASE("num_samples") { return Py_BuildValue("i", self->table->num_samples()); }
	GKPY_GETATTR_CASE("sample_names")
	{
		auto size = self->table->num_samples();
		auto list = PyList_New(size);
		auto name = self->table->sample_names();
		for (int i = 0; i < size; ++i) {
			auto len = (Py_ssize_t)strlen(name);
			PyList_SET_ITEM(list, i, PyString_FromStringAndSize(name, len));
			name += len + 1;
		}
		return list;
	}
	return PyObject_GenericGetAttr(selfo, attro);
GKPY_GETATTRO_END

GKPY_SETATTRO_BEGIN(VCFTable)
	GKPY_SETATTR_READONLY("info_ids");
	GKPY_SETATTR_READONLY("format_ids");
	GKPY_SETATTR_READONLY("filename");
	GKPY_SETATTR_READONLY("num_samples");
	GKPY_SETATTR_READONLY("sample_names");
	return PyObject_GenericSetAttr(selfo, attro, value);
GKPY_SETATTRO_END

template <typename T>
T get_numpy_scalar(PyObject* obj)
{
	T x;
	PyArray_ScalarAsCtype(obj, &x);
	return x;
}

template <class T>
void collect_ids(PyObject* ids, T collect_fn)
{
	if (!ids || ids == Py_None) {
		return;
	}

	if (PyList_Check(ids)) {
		// Extract fields specified as a list
		Py_ssize_t count = PyList_GET_SIZE(ids);
		for (Py_ssize_t i_list = 0; i_list < count; i_list++) {
			PyObject* id = PyList_GetItem(ids, i_list);

			GK_CHECK(PyString_Check(id), type, "Each ID must be a 'str'.");
			const char* id_cstr = PyString_AS_STRING(id);
			GK_CHECK(PyObject_HasAttr((PyObject*)PyVCFVariant::DefaultType, id) == 0, value,
					 "ID {} clashes with VCFVariant attribute.", id_cstr);

			collect_fn(id_cstr, nullopt, nullptr);
		}
	} else if (PyDict_Check(ids)) {
		// Extract fields specified as a dict
		PyObject* items  = PyDict_Items(ids);
		Py_ssize_t count = PyList_GET_SIZE(items);
		for (Py_ssize_t i_list = 0; i_list < count; i_list++) {
			PyObject* item          = PyList_GetItem(items, i_list);
			PyObject* id            = PyTuple_GetItem(item, 0);
			PyObject* default_value = PyTuple_GetItem(item, 1);

			GK_CHECK(PyString_Check(id), type, "Each ID must be a 'str'.");
			const char* id_cstr = PyString_AS_STRING(id);
			GK_CHECK(PyObject_HasAttr((PyObject*)PyVCFVariant::DefaultType, id) == 0, value,
					 "ID {} clashes with VCFVariant attribute.", id_cstr);

			if (default_value == Py_None) {
				collect_fn(id_cstr, nullopt, nullptr);
				continue;
			}

			if (PyInt_Check(default_value)) {
				int default_int = int_cast<int>(PyInt_AS_LONG(default_value));
				collect_fn(id_cstr, vt::int32, &default_int);
			} else if (PyLong_Check(default_value)) {
				int default_int = int_cast<int>(PyLong_AsLong(default_value));
				collect_fn(id_cstr, vt::int32, &default_int);
			} else if (PyFloat_Check(default_value)) {
				auto default_float = (float)(PyFloat_AS_DOUBLE(default_value));
				collect_fn(id_cstr, vt::float32, &default_float);
			} else if (PyObject_IsInstance(default_value, PyArray_TypeObjectFromType(NPY_INT8))) {
				int default_int = get_numpy_scalar<int8_t>(default_value);
				collect_fn(id_cstr, vt::int8, &default_int);
			} else if (PyObject_IsInstance(default_value, PyArray_TypeObjectFromType(NPY_INT16))) {
				int default_int = get_numpy_scalar<int16_t>(default_value);
				collect_fn(id_cstr, vt::int16, &default_int);
			} else if (PyObject_IsInstance(default_value, PyArray_TypeObjectFromType(NPY_INT32))) {
				int default_int = get_numpy_scalar<int32_t>(default_value);
				collect_fn(id_cstr, vt::int32, &default_int);
			} else if (PyObject_IsInstance(default_value, PyArray_TypeObjectFromType(NPY_FLOAT16))) {
				float default_float = as_float(get_numpy_scalar<half_t>(default_value));
				collect_fn(id_cstr, vt::float16, &default_float);
			} else if (PyObject_IsInstance(default_value, PyArray_TypeObjectFromType(NPY_FLOAT32))) {
				auto default_float = get_numpy_scalar<float>(default_value);
				collect_fn(id_cstr, vt::float32, &default_float);
			} else {
				GK_THROW(type, "Default field values must be one of [int, long, float, numpy.int8, "
							   "numpy.int16, numpy.int32, numpy.float16, numpy.float32].");
			}
		}
	} else {
		GK_THROW(type, "IDs must be a 'list' or 'dict'.");
	}
}

static PyObject* PyVCFTable_build_vcfbin(PyObject* cls, PyObject* args, PyObject* kwds)
{
	using namespace std;

	GKPY_TRY

	static char* kwlist[] = { "outfile", "infile", "reference_genome", "info_ids", "fmt_ids",
							  "validate", "exclude", "allow", "ancestral", nullptr };
	const char* outfile;
	PyObject* pyinfile;
	PyObject* refg;
	PyObject* info_ids  = nullptr;
	PyObject* fmt_ids   = nullptr;
	PyObject* validate  = Py_True;
	PyObject* exclude = nullptr;
	PyObject* allow = nullptr;
	const char* ancestral = "error";

	if (!PyArg_ParseTupleAndKeywords(args, kwds, "sOO|OOO!OOs", kwlist, &outfile, &pyinfile, &refg, &info_ids, &fmt_ids,
									 &PyBool_Type, &validate, &exclude, &allow, &ancestral)) {
		return nullptr;
	}

	// Convert infile to a string, even if user input sys.stdin
	const char* infile;
	if (PyObject* fileno = PyObject_CallMethod(pyinfile, "fileno", nullptr); fileno) {
		// Add reads from standard input
		GKPY_TAKEREF(fileno);
		GK_CHECK(PyInt_AsLong(fileno) == 0, value, "When infile is a file, expected sys.stdin (i.e. fileno() == 0)");
		infile = stdin_path;
	} else {
		PyErr_Clear();
		GK_CHECK(PyString_Check(pyinfile), type, "Expected infile to be a str or sys.stdin");
		infile = PyString_AS_STRING(pyinfile);
	}

	vcf_table::builder builder(infile, as_genome(refg), validate == Py_True);

	try {
		collect_ids(info_ids, [&](auto id, auto dtype, auto value) { builder.collect_info(id, dtype, value); });
	}
	GK_RETHROW("INFO field");
	try {
		collect_ids(fmt_ids, [&](auto id, auto dtype, auto value) { builder.collect_fmt(id, dtype, value); });
	}
	GK_RETHROW("FORMAT field");

	if (exclude && exclude != Py_None) {
		GK_DBASSERT(PyList_Check(exclude), "Expected exclude to be a list of intervals.");
		// Add excluded intervals
		for (Py_ssize_t i = 0; i < PyList_GET_SIZE(exclude); ++i) {
			PyObject* interval = PyList_GET_ITEM(exclude, i); // borrowed reference
			GK_CHECK(PyInterval::check(interval), type, "Each exclude item must be an Interval");
			builder.exclude(PyInterval::value(interval));
		}
	}

	if (allow && allow != Py_None) {
		GK_DBASSERT(PyList_Check(allow), "Expected allow to be a list of intervals.");
		// Add allowed intervals
		for (Py_ssize_t i = 0; i < PyList_GET_SIZE(allow); ++i) {
			PyObject* interval = PyList_GET_ITEM(allow, i); // borrowed reference
			GK_CHECK(PyInterval::check(interval), type, "Each allow item must be an Interval");
			builder.allow(PyInterval::value(interval));
		}
	}

	if (ancestral) {
		constexpr string_view actions[] = { "error"sv, "warn"sv, "exclude"sv };
		const auto start                = cbegin(actions);
		const auto stop                 = cend(actions);
		const auto found                = find(start, stop, ancestral);
		GK_CHECK(found != stop, value, "ancestral must be in [\"error\", \"warn\", \"exclude\"]");
		builder.ancestral(vcf_table::builder::action_t(distance(start, found)));
	}

	// Build the file and return.
	builder.build(outfile);
	GKPY_RETURN_NONE;

	GKPY_CATCH_RETURN_NULL
}

static PyObject* PyVCFTable_vcfbin_version(PyObject* cls, PyObject* args)
{
	if (!PyArg_UnpackTuple(args, "", 0, 0))
		return nullptr;
	return PyInt_FromLong((long)vcf_table::file_version());
}

GKPY_OMETHOD_BEGIN(VCFTable, info)
	const char* id;
	if (!PyArg_ParseTuple(args, "s", &id))
		return nullptr;

	const vcf_table::field_col_t* col = self->table->info_fields().get(id);
	GK_CHECK(col != nullptr, key, "Unrecognized INFO ID \"{}\"", id);
	GK_CHECK(col->dtype != vt::str, type, "String INFO columns must be retrieved by VCFVariant attribute.");

	// depth of 1 implies 1D, depth of > 1 implies 2D
	PyObject* shape = PyTuple_New(max(1, col->depth));
	PyTuple_SET_ITEM(shape, 0, PyInt_FromLong(self->table->size()));
	if (col->depth > 1)
		PyTuple_SET_ITEM(shape, 1, PyInt_FromLong(col->depth));

	return PyObject_CallFunction(PyVCFTable::numpy_memmap_function.get(), "OOsKO", self->pystr_filename,
								 (PyObject*)PyArray_DescrFromType(py_dtypes[col->dtype]), "r", col->data_file_offset,
								 shape);
GKPY_OMETHOD_END

GKPY_OMETHOD_BEGIN(VCFTable, format)
	const char* id;
	if (!PyArg_ParseTuple(args, "s", &id))
		return nullptr;

	const vcf_table::field_col_t* col = self->table->fmt_fields().get(id);
	GK_CHECK(col != nullptr, key, "Unrecognized FORMAT ID \"{}\"", id);

	// depth of 1 implies 2D, depth of > 1 implies 3D
	PyObject* shape = PyTuple_New(1 + max(1, col->depth));
	PyTuple_SET_ITEM(shape, 0, PyInt_FromLong(self->table->size()));
	PyTuple_SET_ITEM(shape, 1, PyInt_FromLong(self->table->num_samples()));
	if (col->depth > 1)
		PyTuple_SET_ITEM(shape, 2, PyInt_FromLong(col->depth));

	return PyObject_CallFunction(PyVCFTable::numpy_memmap_function.get(), "OOsKO", self->pystr_filename,
								 (PyObject*)PyArray_DescrFromType(py_dtypes[col->dtype]), "r", col->data_file_offset,
								 shape);
GKPY_OMETHOD_END

GKPY_METHOD_BEGIN_NOARG(VCFTable, close)
	self->table->close();
	GKPY_RETURN_NONE;
GKPY_METHOD_END

static pos_t deletion_length(variant_t variant, interval_t clip)
{
	pos_t ref_len = variant.size();
	auto alt_len  = int_cast<pos_t>(strlen(variant.alt));

	pos_t delta = clip.pos5 - variant.pos5; // clip to start
	if (delta > 0) {
		ref_len = min(ref_len - delta, 0);
		alt_len = min(alt_len - delta, 0);
	}
	delta = variant.pos3 - clip.pos3; // clip to end
	if (delta > 0) {
		ref_len = min(ref_len - delta, 0);
		alt_len = min(alt_len - delta, 0);
	}
	return max(0, ref_len - alt_len);
}

using overlap_t      = vcf_table::cursor_range_overlapping;
using length_range_t = std::pair<pos_t, overlap_t>;

static length_range_t deletion_length(const vcf_table& table, interval_t window)
{
	overlap_t range = table.find_overlapping(window);
	// find an over-estimated count in a permutation of only deletions
	pos_t length = accumulate(cbegin(range), cend(range), 0, [&table, window](auto length, auto index) {
		return length + deletion_length(variant_t(index, table), window);
	});
	return { length, range };
}

static bool append_vcf(PyObject* list, PyObject* vcftable, const packed_variant* elem)
{
	PyAutoRef item(PyTable_CreateItem<PyVCFVariant>(vcftable, elem));
	if (PyList_Append(list, item.get()) != 0)
		return false;
	item.release();
	return true;
}

GKPY_METHOD_BEGIN_ONEARG(VCFTable, sequence_variations)
	ainterval_t interval  = PyAsAnchoredInterval(arg);
	(interval_t&)interval = interval.as_pos_strand(); // slice OK, anchor retained

	const vcf_table& table = *self->table;
	overlap_t range(nullptr, nullptr, nullptr, 0);

	if (interval.anchor == invalid_pos) {
		range = table.find_overlapping(interval);
	} else {
		interval_t window = interval.with_pos(interval.anchor, interval.pos3); // window right of anchor
		while (window.size() > 0) {
			length_range_t length_range = deletion_length(table, window);
			if (range.elems == nullptr) {
				range = length_range.second;
			} else {
				range.last = length_range.second.last;
			}
			window = window.end3().expand(-1, length_range.first);
		}

		window = interval.with_pos(interval.pos5, interval.anchor - 1); // window left of anchor
		while (window.size() > 0) {
			length_range_t length_range = deletion_length(table, window);
			if (range.elems == nullptr) {
				range = length_range.second;
			} else {
				range.first = length_range.second.first;
				range.pos5  = length_range.second.pos5;
			}
			window = window.end5().expand(length_range.first, -1);
		}
	}

	// need to cover insertions at end5 (end3 is covered by find_overlapping)
	using aligned_5p_t = vcf_table::cursor_range;

	pos_t pos5                 = range.elems != nullptr ? range.pos5 : interval.pos5;
	aligned_5p_t aligned_5p    = table.find_5p_aligned(interval.with_pos(pos5, pos5));
	overlap_t::cursor_t it     = range.begin();
	overlap_t::cursor_t it_end = range.end();
	PyAutoRef variants(PyList_New(0));

	// set union
	for (aligned_5p_t::cursor_t it_5p = aligned_5p.begin(), it_5p_end = aligned_5p.end(); it_5p != it_5p_end;) {
		if (it == it_end) {
			for (; it_5p != it_5p_end; ++it_5p)
				if (it_5p->size() == 0 // only add insertions
					&& !append_vcf(variants.get(), self, &*it_5p))
					return nullptr;
			return variants.release();
		} else if (it_5p->size() != 0) { // only add insertions
			++it_5p;
			continue;
		} else if (it->pos5 < it_5p->pos5 || (it->pos5 == it_5p->pos5 && it->pos3 < it_5p->pos3)) {
			if (!append_vcf(variants.get(), self, &*it))
				return nullptr;
			++it;
		} else {
			if (!append_vcf(variants.get(), self, &*it_5p))
				return nullptr;
			if (&*it == &*it_5p)
				++it;
			++it_5p;
		}
	}
	for (; it != it_end; ++it)
		if (!append_vcf(variants.get(), self, &*it))
			return nullptr;
	return variants.release();
GKPY_METHOD_END

static vt::gt_t get_gt(const vcf_table::field_col_t& gt, size_t index)
{
	GK_ASSERT(gt.dtype == vt::int8);
	return vt::gt_t(rcast<const uint8_t*>(gt.data)[index]);
}

static float add_af(const vcf_table::field_col_t& af, index_t index_x, index_t index_y)
{
	float freq = 0;
	switch (af.dtype) {
	case vt::float16: {
		const auto* afs = rcast<const half_t*>(af.data);
		freq            = as_float(afs[index_x]) + as_float(afs[index_y]);
		break;
	}
	case vt::float32: {
		const auto* afs = rcast<const float*>(af.data);
		freq            = afs[index_x] + afs[index_y];
		break;
	}
	default: GK_THROW(value, "AF column must be defined as a float type.");
	}
	return freq;
}

GKPY_OMETHOD_BEGIN(VCFTable, variant_combinations)
	static char* kwlist[]  = { "variants", "heterozygous_alt_af_threshold", nullptr };
	PyObject* variants_seq = nullptr;
	float threshold        = 0.9f;
	if (!PyArg_ParseTupleAndKeywords(args, kwds, "O|f", kwlist, &variants_seq, &threshold))
		return nullptr;
	GK_CHECK(PySequence_Check(variants_seq), type, "variants must be a sequence.");

	Py_ssize_t variants_len = PySequence_Length(variants_seq);
	// rough upperbound: sum(nCi)<=2^n for i=0..n
	GK_CHECK(variants_len < (Py_ssize_t)sizeof(size_t) * CHAR_BIT - 1, value,
			 "Too many variants: count {} is greater than maximum count {}.", variants_len,
			 sizeof(size_t) * CHAR_BIT - 1);

	vector<PyVCFVariant*> variants(variants_len);
	for (Py_ssize_t i = 0; i < variants_len; ++i) {
		PyAutoRef item(PySequence_ITEM(variants_seq, i));
		GKPY_TYPECHECK(item.get(), PyVCFVariant::DefaultType);
		variants[i] = (PyVCFVariant*)item.get(); // ref incremented in combos below
	}

	const auto& table = *self->table;
	auto gt_cols      = table.fmt_fields().get("GT");
	auto af_cols      = table.info_fields().get("AF");

	// nCr variants, with additional filtering
	vector<bool> selection_mask;
	vector<PyObject*> combination; // vector is faster than PyTuple if combination is invalid
	PyAutoRef combos(PyList_New(0));
	for (Py_ssize_t r = 1; r <= variants_len; ++r) {
		selection_mask.clear();
		selection_mask.resize(r, true);
		selection_mask.resize(variants_len);

		do {
			combination.clear();
			for (size_t i = 0; i < variants.size(); ++i) {
				const packed_variant& variant = variants[i]->value();
				int table_index               = table.index_of(variants[i]->value());

				vt::gt_t gt = gt_cols ? get_gt(*gt_cols, table_index) : vt::gt_heterozygous;

				bool impossible_combo = false;
				// filter out independent conditions first
				if ((selection_mask[i] && gt != vt::gt_heterozygous && gt != vt::gt_homozygous_alt)
					|| (!selection_mask[i] && gt == vt::gt_homozygous_alt)) {
					impossible_combo = true;
				} else if (gt == vt::gt_heterozygous || gt == vt::gt_homozygous_alt) {
					// conditions dependent on overlapping variants
					for (size_t j = i + 1; j < variants.size(); ++j) {
						const packed_variant& next = variants[j]->value();
						if (!variant.overlaps(next))
							continue;

						if ((selection_mask[i] && selection_mask[j])
							// check for split heterozygous_alt; no representation in vcf:
							// https://github.com/samtools/hts-specs/issues/77
							|| (!selection_mask[i] && !selection_mask[j] && variant.pos5 == next.pos5 && af_cols
								&& add_af(*af_cols, table_index, table.index_of(next)) >= threshold)) {
							impossible_combo = true;
							break;
						}
					}
				}
				if (impossible_combo) {
					combination.clear();
					break;
				} else if (selection_mask[i])
					combination.push_back(variants[i]);
			}
			if (!combination.empty()) {
				PyAutoRef combo(PyTuple_New((Py_ssize_t)combination.size()));
				for (size_t i = 0; i < combination.size(); ++i) {
					PyTuple_SET_ITEM(combo.get(), (Py_ssize_t)i, combination[i]);
					Py_INCREF(combination[i]);
				}
				if (PyList_Append(combos.get(), combo.get()) != 0)
					return nullptr;
				combo.release();
			}
		} while (prev_permutation(selection_mask.begin(), selection_mask.end()));
	}
	return combos.release();
GKPY_OMETHOD_END

GKPY_METHODS_BEGIN(VCFTable)
	GKPY_METHOD_ENTRY(VCFTable, vcfbin_version, METH_VARARGS | METH_STATIC, nullptr)
	GKPY_METHOD_ENTRY(VCFTable, build_vcfbin, METH_VARARGS | METH_KEYWORDS | METH_STATIC, nullptr)
	GKPY_METHOD_ENTRY(VCFTable, info, METH_VARARGS | METH_KEYWORDS, nullptr)
	GKPY_METHOD_ENTRY(VCFTable, format, METH_VARARGS | METH_KEYWORDS, nullptr)
	GKPY_METHOD_ENTRY(VCFTable, close, METH_NOARGS, nullptr)
	GKPY_METHOD_ENTRY(VCFTable, sequence_variations, METH_O, nullptr)
	GKPY_METHOD_ENTRY(VCFTable, variant_combinations, METH_VARARGS | METH_KEYWORDS, nullptr)
GKPY_METHODS_END

GKPY_SUBTYPEOBJ_BEGIN(VCFTable)
	tp_flags |= Py_TPFLAGS_BASETYPE;
	tp_new      = PyVCFTable_New;
	tp_init     = PyVCFTable_Init;
	tp_dealloc  = PyVCFTable_Dealloc;
	tp_getattro = PyVCFTable_GetAttro;
	tp_setattro = PyVCFTable_SetAttro;
	tp_methods  = PyVCFTable_Methods;
GKPY_SUBTYPEOBJ_END

PyObject* PyVCFTable::get_col_attr(const packed_variant* variant, const char* id)
{
	auto col = table->info_fields().get(id);
	if (col == nullptr) {
		col = table->fmt_fields().get(id);
		if (col == nullptr) {
			return nullptr;
		}
	}

	GK_CHECK(col->depth == 1, value,
			 "VCF INFO/FORMAT attributes only supports scalar values (Number=1). Use VCFTable.info or VCFTable.format "
			 "for multidimensional access.");
	auto index = scast<size_t>(variant - &(*table)[0]);
	switch (col->dtype) {
	case vt::bool_: return PyReturnBool(rcast<const uint8_t*>(col->data)[index] != 0);
	case vt::int8: return PyInt_FromLong(rcast<const uint8_t*>(col->data)[index]);
	case vt::int16: return PyInt_FromLong(rcast<const uint16_t*>(col->data)[index]);
	case vt::int32: return PyInt_FromLong(rcast<const uint32_t*>(col->data)[index]);
	case vt::float16: return PyFloat_FromDouble(as_float(rcast<const half_t*>(col->data)[index]));
	case vt::float32: return PyFloat_FromDouble(rcast<const float*>(col->data)[index]);
	case vt::str: return get_col_attr_str(index, col);
	default: GK_UNREACHABLE();
	}
	return nullptr;
}

PyObject* PyVCFTable::get_col_attr_str(size_t index, const vcf_table::field_col_t* col)
{
	const auto* data = rcast<const uint8_t*>(col->data);
	auto smode       = vt::str_mode_t(*rcast<const uint32_t*>(data));
	switch (smode) {
	case vt::fixlen: {
		// Compute string offset from variant index
		uint32_t stride        = *rcast<const uint32_t*>(data + 4);
		const char* fixed_data = rcast<const char*>(data + 8);
		const char* str        = fixed_data + stride * index;
		return PyString_FromString(str); // new ref
	}
	case vt::varlen: {
		// Look up string offset by variant index
		const auto* offsets = rcast<const offset40_t*>(data + 4);
		uint64_t offset     = offsets[index].as64();
		uint64_t len        = offsets[index + 1].as64() - offset - 1;
		const char* str     = (const char*)offsets + offset;
		return PyString_FromStringAndSize(str, len); // new ref
	}
	case vt::pooled: {
		// Look up string pool index by variant index
		uint32_t pool_index_size = *rcast<const uint32_t*>(data + 4);
		uint64_t num_pooled      = *rcast<const uint64_t*>(data + 8);
		const auto* offsets      = rcast<const uint64_t*>(data + 16);
		const void* pool_indices = offsets + num_pooled + 1;

		// Get index of pooled string for this variant
		size_t pool_index;
		switch (pool_index_size) {
		case 1: pool_index = rcast<const uint8_t*>(pool_indices)[index]; break;
		case 2: pool_index = rcast<const uint16_t*>(pool_indices)[index]; break;
		case 4: pool_index = rcast<const uint32_t*>(pool_indices)[index]; break;
		case 8: pool_index = rcast<const uint64_t*>(pool_indices)[index]; break;
		default: GK_UNREACHABLE();
		}

		// Fetch the string from the PyString pool, creating a new instance if required.
		if (pystr_pool == nullptr) {
			pystr_pool = new std::remove_pointer_t<decltype(pystr_pool)>;
		}
		auto& pool = (*pystr_pool)[col];
		if (pool.empty())
			pool.resize(num_pooled);
		PyObject* pystr = pool[pool_index].get();
		if (!pystr) {
			uint64_t offset = offsets[pool_index];
			uint64_t len    = offsets[pool_index + 1] - offset - 1;
			const char* str = (const char*)offsets + offset;
			pystr           = PyString_FromStringAndSize(str, len); // new ref
			pool[pool_index].reset(pystr);
		}
		Py_INCREF(pystr);
		return pystr;
	}
	}
	return nullptr;
}

////////////////////////////////////////////////////////////
//  VariantTable
////////////////////////////////////////////////////////////

void validate_VariantTable(const PyAsPtrSource* self)
{
	GK_CHECK(((const PyVariantTable*)self)->table->valid(), file,
			 "Variants have been invalidated by close or with statement on the VariantTable source.");
}

GKPY_NEW_BEGIN(VariantTable)
	self->owner = nullptr;
	// VariantTable currently operates in only one mode:
	//    the variant_table instance that it wraps belongs to some other C++ object,
	//    the `owner`, and so this VariantTable instance merely points to that table,
	//    and does not own the table itself.
	// The above use case is different from the (yet unimplemented) use case where
	// the user creates a VariantTable in isolation, in memory, and so the underlying
	// variant_table instance is owned by the VariantTable itself.
	//
	// Currently:
	//   args[0] should be of type variant_table*, from which elements will be pulled.
	//   args[1] is optional and can be PyObject* that this instance should hold a reference to,
	//           to prevent the `owner` from being collected before this instance.
	if (!PyArg_ParseTuple(args, "On", &self->owner, &self->table))
		return nullptr;
	Py_XINCREF(self->owner);

	self->set_validator(&validate_VariantTable);

GKPY_NEW_END

GKPY_DEALLOC_OWNED_BEGIN(VariantTable)
GKPY_DEALLOC_OWNED_END

GKPY_SUBTYPEOBJ_BEGIN(VariantTable)
	tp_flags |= Py_TPFLAGS_BASETYPE;
	tp_new     = PyVariantTable_New;
	tp_dealloc = PyVariantTable_Dealloc;
GKPY_SUBTYPEOBJ_END

END_NAMESPACE_GK
