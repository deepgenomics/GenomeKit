/*
Copyright (C) 2016-2023 Deep Genomics Inc. All Rights Reserved.
*/
#include "py_util.h" // Include first to avoid POSIX warnings due to Python header inclusion order
#include "py_genome_track.h"
#include "py_genome.h"
#include "py_interval.h"
#include "genome_track.h"

#include <numpy/arrayobject.h>

BEGIN_NAMESPACE_GK

using etype_t                  = genome_track::dtype_t;
using dtype_t                  = genome_track::dtype_t;
static const dtype_t num_dtype = genome_track::num_dtype;

static int py_dtypes[num_dtype] = {
	NPY_BOOL,  // bool_
	NPY_UINT8, // uint8
	NPY_INT8,  // int8
	NPY_HALF,  // float16
	NPY_FLOAT, // float32
};

static PyTypeObject* py_dtype_type(dtype_t dtype)
{
	static PyTypeObject* _types[num_dtype] = {
		&PyBoolArrType_Type,  // bool_
		&PyUInt8ArrType_Type, // uint8
		&PyInt8ArrType_Type,  // int8
		&PyHalfArrType_Type,  // float16
		&PyFloatArrType_Type, // float32
	};
	return _types[dtype];
}

static dtype_t dtype_from_py(int py_dtype)
{
	for (int dtype = 0; dtype < num_dtype; ++dtype)
		if (py_dtypes[dtype] == py_dtype)
			return (dtype_t)dtype;
	GK_THROW(type, "data array had unrecognized dtype; try np.bool_, np.uint8, np.int8, np.float16, or np.float32");
}

/////////////////////////////////////////////////////

GKPY_NEW_BEGIN(GenomeTrackBuilder)
	// We expect to be called with __new__(*config) where each item
	// in the config tuple is a string.
	int dim = 1;
	int res = 1;
	PyObject* refg;
	const char* etype_str = nullptr;
	const char* outfile   = nullptr;
	const char*  c_strandedness = nullptr;
	static char* kwlist[] = { "outfile", "etype", "strandedness", "reference_genome", "dim", "resolution", nullptr };
	if (!PyArg_ParseTupleAndKeywords(args, kwds, "sssO|ii", kwlist, &outfile, &etype_str, &c_strandedness, &refg, &dim, &res))
		return nullptr;
	const auto is_single_stranded = strcmp("single_stranded", c_strandedness) == 0;
	const auto is_strand_unaware  = strcmp("strand_unaware", c_strandedness) == 0;
	const auto is_strand_aware    = strcmp("strand_aware", c_strandedness) == 0;
	GK_ASSERT(is_single_stranded || is_strand_unaware || is_strand_aware,
			  "strandedness must be one of 'single_stranded', 'strand_unaware', or 'strand_aware'");
	auto strandedness = is_single_stranded ? genome_track::strandedness_t::single_stranded
		: is_strand_unaware ? genome_track::strandedness_t::strand_unaware
		: genome_track::strandedness_t::strand_aware;
	self->builder = new genome_track::builder(outfile, genome_track::as_etype(etype_str), strandedness, as_genome(refg), dim, res);
	self->genome  = refg;
	Py_IncRef(self->genome);
GKPY_NEW_END

GKPY_INIT_EMPTY(GenomeTrackBuilder)

GKPY_DEALLOC_BEGIN(GenomeTrackBuilder)
	// If we don't have an owner, then the genome_track instance must be deleted
	if (self->builder)
		delete self->builder;
	Py_DecRef(self->genome);
GKPY_DEALLOC_END

GKPY_GETATTRO_BEGIN(GenomeTrackBuilder)
	GKPY_GETATTR_CASE("dim")               return PyInt_FromLong(self->builder->dim());
	GKPY_GETATTR_CASE("resolution")        return PyInt_FromLong(self->builder->res());
	GKPY_GETATTR_CASE("stranded") GKPY_RETURN_BOOL(self->builder->stranded());
	GKPY_GETATTR_CASE("etype")             return PyString_FromString(genome_track::etype_as_cstr[self->builder->etype()]);
	GKPY_GETATTR_CASE("dtype") GKPY_RETURN_INCREF(rcast<PyObject*>(py_dtype_type(self->builder->dtype())));
	GKPY_GETATTR_CASE("reference_genome")  return PyString_FromSV(self->builder->refg_name());
	GKPY_GETATTR_CASE("refg")              return PyString_FromSV(self->builder->refg_name());
	GKPY_GETATTR_CASE("filename")          return PyString_FromSV(self->builder->filename());
	GKPY_GETATTR_CASE("data_size")         return PyInt_FromSize_t(self->builder->data_size());
	GKPY_GETATTR_CASE("index_size")        return PyInt_FromSize_t(self->builder->index_size());
	return PyObject_GenericGetAttr(selfo, attro);
GKPY_GETATTRO_END

GKPY_SETATTRO_BEGIN(GenomeTrackBuilder)
	GKPY_SETATTR_READONLY("dim");
	GKPY_SETATTR_READONLY("resolution");
	GKPY_SETATTR_READONLY("stranded");
	GKPY_SETATTR_READONLY("etype");
	GKPY_SETATTR_READONLY("dtype");
	GKPY_SETATTR_READONLY("reference_genome");
	GKPY_SETATTR_READONLY("refg");
	GKPY_SETATTR_READONLY("filename");
	GKPY_SETATTR_READONLY("data_size");
	GKPY_SETATTR_READONLY("index_size");
	return PyObject_GenericSetAttr(selfo, attro, value);
GKPY_SETATTRO_END

GKPY_OMETHOD_BEGIN(GenomeTrackBuilder, set_default_value)
	float value;
	static char* kwlist[] = { "value", nullptr };
	if (!PyArg_ParseTupleAndKeywords(args, kwds, "f", kwlist, &value))
		return nullptr;
	self->builder->set_default_value(value);
	GKPY_RETURN_NONE;
GKPY_OMETHOD_END

GKPY_OMETHOD_BEGIN(GenomeTrackBuilder, set_sparsity)
	int min_run = 64;
	float min_delta = 0.0f;
	static char* kwlist[] = { "min_run", "min_delta", nullptr };
	if (!PyArg_ParseTupleAndKeywords(args, kwds, "|if", kwlist, &min_run, &min_delta))
		return nullptr;
	self->builder->set_sparsity(min_run, min_delta);
	GKPY_RETURN_NONE;
GKPY_OMETHOD_END

GKPY_METHOD_BEGIN_NOARG(GenomeTrackBuilder, set_clamping)
	self->builder->set_clamping();
	GKPY_RETURN_NONE;
GKPY_METHOD_END

GKPY_OMETHOD_BEGIN(GenomeTrackBuilder, set_dict)
	PyArrayObject* py_dict;
	static char* kwlist[] = { "dict", nullptr };
	if (!PyArg_ParseTupleAndKeywords(args, kwds, "O!", kwlist, &PyArray_Type, &py_dict))
		return nullptr;

	int dict_size = 0;
	switch (self->builder->etype()) {
	case genome_track::f2: dict_size = 4;   break;
	case genome_track::f3: dict_size = 8;   break;
	case genome_track::f4: dict_size = 16;  break;
	case genome_track::f5: dict_size = 32;  break;
	case genome_track::f6: dict_size = 64;  break;
	case genome_track::f8: dict_size = 256; break;
	default: GK_THROW(value, "Cannot call set_dict on etype '{}'", genome_track::etype_as_cstr[self->builder->etype()]);
	}

	// Check the data array to make sure it's contiguous, C-order, 1- or 2-dimensions, and the right size.
	GK_CHECK(PyArray_NDIM(py_dict) == 1, value, "dict array must be 1-dimensional");
	GK_CHECK((int)PyArray_DIMS(py_dict)[0] == dict_size, value, "dict array must be of size {} for this etype", dict_size);

	// Set the dict as either float16 or float32
	if      (PyArray_TYPE(py_dict) == NPY_HALF)  self->builder->set_dict(rcast<const half_t*>(PyArray_DATA(py_dict)));
	else if (PyArray_TYPE(py_dict) == NPY_FLOAT) self->builder->set_dict(rcast<const float* >(PyArray_DATA(py_dict)));
	else GK_THROW(type, "dict must have dtype of np.float16 or np.float32");
	GKPY_RETURN_NONE;
GKPY_OMETHOD_END

GKPY_OMETHOD_BEGIN(GenomeTrackBuilder, set_transform)
	float a, b, c, d;
	static char* kwlist[] = { "a", "b", "c", "d", nullptr };
	if (!PyArg_ParseTupleAndKeywords(args, kwds, "ffff", kwlist, &a, &b, &c, &d))
		return nullptr;
	self->builder->set_transform(a, b, c, d);
	GKPY_RETURN_NONE;
GKPY_OMETHOD_END

GKPY_OMETHOD_BEGIN(GenomeTrackBuilder, set_restriction)
	PyInterval* restriction;
	static char* kwlist[] = { "restriction", nullptr };
	if (!PyArg_ParseTupleAndKeywords(args, kwds, "O!", kwlist, PyInterval::DefaultType, &restriction))
		return nullptr;
	self->builder->set_restriction(restriction->value());
	GKPY_RETURN_NONE;
GKPY_OMETHOD_END

GKPY_OMETHOD_BEGIN(GenomeTrackBuilder, set_data)
	PyInterval* py_interval;
	PyArrayObject* py_data;
	static char* kwlist[] = { "interval", "data", nullptr };
	if (!PyArg_ParseTupleAndKeywords(args, kwds, "O!O", kwlist, PyInterval::DefaultType, &py_interval, &py_data))
		return nullptr;

	const interval_t& interval = py_interval->value();
	dtype_t dtype;
	void* data;

	if (self->builder->etype() == genome_track::m0) {

		// Mask tracks are special case, since no data.
		GK_CHECK((PyObject*)py_data == Py_None, value, "Data must be None for etype 'm0'");
		dtype = genome_track::bool_;
		data  = nullptr;

	} else {

		// Otherwise we expect a numpy array with data for interval.
		GK_CHECK(PyArray_Check(py_data), type, "Data must be numpy array");
		int ndim = PyArray_NDIM(py_data);
		int res = self->builder->res();

		// Check the data array to make sure it's contiguous, C-order, and 1- or 2-dimensions.
		GK_CHECK(ndim == 1 || ndim == 2, value, "Data must be 1- or 2-dimensional");
		GK_CHECK((int)PyArray_DIMS(py_data)[0]*res == interval.size(), value, "Data must have {} rows", interval.size());
		if (ndim > 1) {
			GK_CHECK((int)PyArray_DIMS(py_data)[1] == self->builder->dim(), value, "Data must have {} columns", self->builder->dim());
			GK_CHECK(PyArray_FLAGS(py_data) & NPY_ARRAY_CARRAY_RO, value, "Multi-dimensional data array must be C_CONTIGUOUS order");
		}

		// Set the data for the given interval using the numpy dtype, if it's supported
		dtype = dtype_from_py(PyArray_TYPE(py_data));
		// TODO: handle strided date during encoding
		size_t stride = genome_track::dtype_size[dtype];
		for (int dim = ndim; dim > 0; --dim) {
			GK_CHECK((size_t)PyArray_STRIDE(py_data, dim - 1) == stride, value,
					 "Data detected having stride!=1, possibly because it is transposed. "
					 "Consider making a copy with `np.array` to get non-strided data.");
			stride *= (size_t)PyArray_DIM(py_data, dim - 1);
		}
		data = PyArray_DATA(py_data);
	}

	self->builder->set_data(interval, data, dtype);

	GKPY_RETURN_NONE;
GKPY_OMETHOD_END

GKPY_OMETHOD_BEGIN(GenomeTrackBuilder, set_data_from_wig)
	const char* pos_strand_file;
	const char* neg_strand_file = nullptr;
	static char* kwlist[]       = { "pos_strand_file", "neg_strand_file", nullptr };
	if (!PyArg_ParseTupleAndKeywords(args, kwds, "s|s", kwlist, &pos_strand_file, &neg_strand_file))
		return nullptr;
	self->builder->set_data_from_wig(pos_strand_file, neg_strand_file);
	GKPY_RETURN_NONE;
GKPY_OMETHOD_END

GKPY_OMETHOD_BEGIN(GenomeTrackBuilder, set_data_from_bedgraph)
	const char* pos_strand_file;
	const char* neg_strand_file = nullptr;
	static char* kwlist[]       = { "pos_strand_file", "neg_strand_file", nullptr };
	if (!PyArg_ParseTupleAndKeywords(args, kwds, "s|s", kwlist, &pos_strand_file, &neg_strand_file))
		return nullptr;
	self->builder->set_data_from_bedgraph(pos_strand_file, neg_strand_file);
	GKPY_RETURN_NONE;
GKPY_OMETHOD_END

GKPY_OMETHOD_BEGIN(GenomeTrackBuilder, set_data_from_bed)
	PyObject* py_categories = nullptr;
	const char* bedfile;
	static char* kwlist[] = { "bedfile", "categories", nullptr };
	if (!PyArg_ParseTupleAndKeywords(args, kwds, "s|O!", kwlist, &bedfile,
	                                 &PyList_Type, &py_categories))
		return nullptr;

	// Convert categories to vector<string>, if specified
	vector<string> categories;
	if (py_categories) {
		for (Py_ssize_t i = 0; i < PyList_GET_SIZE(py_categories); ++i) {
			PyObject* item = PyList_GET_ITEM(py_categories, i); // borrowed reference
			GK_CHECK(PyString_Check(item), type, "Each category item must be a string");
			categories.emplace_back(PyString_AS_STRING(item));
		}
	}

	self->builder->set_data_from_bed(bedfile, categories);
	GKPY_RETURN_NONE;
GKPY_OMETHOD_END

GKPY_METHOD_BEGIN_NOARG(GenomeTrackBuilder, flush)
	self->builder->flush();
	GKPY_RETURN_NONE;
GKPY_METHOD_END

GKPY_METHOD_BEGIN_NOARG(GenomeTrackBuilder, finalize)
	self->builder->finalize();
	GKPY_RETURN_NONE;
GKPY_METHOD_END

GKPY_METHOD_BEGIN_NOARG(GenomeTrackBuilder, close)
	self->builder->close();
	GKPY_RETURN_NONE;
GKPY_METHOD_END

GKPY_METHODS_BEGIN(GenomeTrackBuilder)
	GKPY_METHOD_ENTRY(GenomeTrackBuilder, set_default_value,      METH_VARARGS | METH_KEYWORDS, nullptr)
	GKPY_METHOD_ENTRY(GenomeTrackBuilder, set_sparsity,           METH_VARARGS | METH_KEYWORDS, nullptr)
	GKPY_METHOD_ENTRY(GenomeTrackBuilder, set_clamping,           METH_VARARGS | METH_KEYWORDS, nullptr)
	GKPY_METHOD_ENTRY(GenomeTrackBuilder, set_dict,               METH_VARARGS | METH_KEYWORDS, nullptr)
	GKPY_METHOD_ENTRY(GenomeTrackBuilder, set_transform,          METH_VARARGS | METH_KEYWORDS, nullptr)
	GKPY_METHOD_ENTRY(GenomeTrackBuilder, set_restriction,        METH_VARARGS | METH_KEYWORDS, nullptr)
	GKPY_METHOD_ENTRY(GenomeTrackBuilder, set_data,               METH_VARARGS | METH_KEYWORDS, nullptr)
	GKPY_METHOD_ENTRY(GenomeTrackBuilder, set_data_from_wig,      METH_VARARGS | METH_KEYWORDS, nullptr)
	GKPY_METHOD_ENTRY(GenomeTrackBuilder, set_data_from_bedgraph, METH_VARARGS | METH_KEYWORDS, nullptr)
	GKPY_METHOD_ENTRY(GenomeTrackBuilder, set_data_from_bed,      METH_VARARGS | METH_KEYWORDS, nullptr)
	GKPY_METHOD_ENTRY(GenomeTrackBuilder, flush,                  METH_NOARGS, nullptr)
	GKPY_METHOD_ENTRY(GenomeTrackBuilder, finalize,               METH_NOARGS, nullptr)
	GKPY_METHOD_ENTRY(GenomeTrackBuilder, close,                  METH_NOARGS, nullptr)
GKPY_METHODS_END

GKPY_TYPEOBJ_BEGIN(GenomeTrackBuilder)
	tp_flags |= Py_TPFLAGS_BASETYPE | Py_TPFLAGS_HAVE_GC;
	tp_new = PyGenomeTrackBuilder_New;
	tp_init = PyGenomeTrackBuilder_Init;
	tp_dealloc = PyGenomeTrackBuilder_Dealloc;
	tp_getattro = PyGenomeTrackBuilder_GetAttro;
	tp_setattro = PyGenomeTrackBuilder_SetAttro;
	tp_methods = PyGenomeTrackBuilder_Methods;
GKPY_TYPEOBJ_END

/////////////////////////////////////////////////////////////

GKPY_INIT_BEGIN(GenomeTrack)
	// Allow GenomeTrack to be constructed one of two ways:
	//    GenomeTrack(gtrack: str)                               <-- used by users loading custom tracks
	//    GenomeTrack(genome: Genome, track_ptr: Py_ssize_t)     <-- used by Genome object for built-in tracks
	PyObject* arg           = nullptr;
	genome_track* track_ptr = nullptr;
	if (!PyArg_ParseTuple(args, "O|n", &arg, &track_ptr))
		return -1;
	if (PyString_Check(arg)) {
		// If argument was a string, set the source file directly
		self->track = new genome_track();
		self->track->set_source(PyString_AS_STRING(arg));
	} else if (PyObject_IsInstance(arg, (PyObject*)PyGenome::Type)) {
		// If argument was a Genome object, then point to the pre-existing genome_dna instance
		Py_INCREF(arg);
		self->owner = arg;
		self->track = track_ptr;
	} else
		GK_THROW(value, "GenomeTrack.__init__ could not parse arguments");
GKPY_INIT_END

GKPY_DEALLOC_OWNED_BEGIN(GenomeTrack)
	// If we don't have an owner, then the genome_track instance must be deleted
	if (self->track && !self->owner) {
		delete self->track;
		self->track = nullptr;
	}
GKPY_DEALLOC_OWNED_END

GKPY_TRAVERSE_OWNED_BEGIN(GenomeTrack)
GKPY_TRAVERSE_OWNED_END

GKPY_CLEAR_OWNED_BEGIN(GenomeTrack)
GKPY_CLEAR_OWNED_END

GKPY_GETATTRO_BEGIN(GenomeTrack)
	GKPY_GETATTR_CASE("dim")               return PyInt_FromLong(self->track->dim());
	GKPY_GETATTR_CASE("resolution")        return PyInt_FromLong(self->track->res());
	GKPY_GETATTR_CASE("stranded") GKPY_RETURN_BOOL(self->track->stranded());
	GKPY_GETATTR_CASE("etype")             return PyString_FromString(genome_track::etype_as_cstr[self->track->etype()]);
	GKPY_GETATTR_CASE("dtype") GKPY_RETURN_INCREF(rcast<PyObject*>(py_dtype_type(self->track->dtype())));
	// TODO: data_dir injected as a context
	GKPY_GETATTR_CASE("reference_genome")  return PyString_FromSV(get_refg_registry().refg_as_sv(self->track->refg()));
	GKPY_GETATTR_CASE("refg")              return PyString_FromSV(get_refg_registry().refg_as_sv(self->track->refg()));
	GKPY_GETATTR_CASE("filename")          return PyString_FromSV(self->track->source());
	GKPY_GETATTR_CASE("intervals")
	{
		auto      intervals = self->track->intervals();
		auto      count     = std::ssize(intervals);
		PyObject* list      = PyList_New(count);
		GKPY_TAKEREF(list);
		for (Py_ssize_t i = 0; i < count; ++i) {
			PyList_SET_ITEM(list, i, PyInterval::create(intervals[i]));
		}
		GKPY_FORGETREF(list);  // keep list refcount at 1
		return list;
	}
	return PyObject_GenericGetAttr(selfo, attro);
GKPY_GETATTRO_END

GKPY_SETATTRO_BEGIN(GenomeTrack)
	GKPY_SETATTR_READONLY("dim");
	GKPY_SETATTR_READONLY("resolution");
	GKPY_SETATTR_READONLY("stranded");
	GKPY_SETATTR_READONLY("etype");
	GKPY_SETATTR_READONLY("dtype");
	GKPY_SETATTR_READONLY("reference_genome");
	GKPY_SETATTR_READONLY("refg");
	GKPY_SETATTR_READONLY("filename");
	return PyObject_GenericSetAttr(selfo, attro, value);
GKPY_SETATTRO_END

// Convert from np.uint8, np.float16, etc to corresponding genome_track dtype enum.
dtype_t dtype_from_obj(PyObject* obj)
{
	GK_CHECK(PyType_Check(obj), type, "Expected numpy dtype but got '{}'", Py_TYPE(obj)->tp_name);
	dtype_t dtype;
	auto* nptype = rcast<PyTypeObject*>(obj);
	for (dtype = (dtype_t)0; dtype < num_dtype; dtype = (dtype_t)((int)dtype+1))
		if (PyType_IsSubtype(nptype, py_dtype_type(dtype)))
			break;
	GK_CHECK(dtype < num_dtype, type, "Unsupported type '{}' for track result", nptype->tp_name);
	return dtype;
}

GKPY_OMETHOD_BEGIN(GenomeTrack, Call)
	PyObject*    itv       = nullptr;
	PyObject*    dtype_arg = nullptr;
	PyObject*    out       = nullptr;
	static char* kwlist[]  = {"interval", "dtype", "out", nullptr};
	if (!PyArg_ParseTupleAndKeywords(args, kwds, "O!|OO", kwlist, PyInterval::DefaultType, &itv, &dtype_arg, &out))
		return nullptr;

	auto      c      = PyAsInterval(itv);
	int       stride = 0;
	PyObject* py_dst = nullptr;
	std::optional<dtype_t> dtype;

	if (!out || out == Py_None) {
		// Create an output array of the right type, size, and dimensionality
		npy_intp dims[2] = {c.size(), self->track->dim()};
		dtype            = dtype_arg && dtype_arg != Py_None ? dtype_from_obj(dtype_arg) : self->track->dtype();
		py_dst           = PyArray_Empty(2, dims, PyArray_DescrFromType(py_dtypes[*dtype]), 0);  // 0 => C_CONTIGUOUS
		if (!py_dst)
			return nullptr;  // Propagate the error up to interpreter immediately
	} else {
		GK_CHECK(PyArray_Check(out), type, "out must be a numpy ndarray.");

		auto out_array = rcast<PyArrayObject*>(out);

		GK_CHECK(PyArray_NDIM(out_array) == 1 || PyArray_NDIM(out_array) == 2, value,
				 "Dimension must be 1- or 2-dimensional: out is {}.", PyArray_NDIM(out_array));
		GK_CHECK(PyArray_DIM(out_array, 0) == c.size(), value, "Row mismatch: out is {} but interval is {}",
				 PyArray_DIM(out_array, 0), c.size());
		if (PyArray_NDIM(out_array) == 2) {
			GK_CHECK(PyArray_DIM(out_array, 1) == self->track->dim(), value,
					 "Column mismatch: out is {} but track is {}", PyArray_DIM(out_array, 1), self->track->dim());
		}
		GK_CHECK(PyArray_ISBEHAVED(out_array), value, "out must be writable from C.");
		dtype  = dtype_from_py(PyArray_DTYPE(out_array)->type_num);
		stride = int_cast<decltype(stride)>(PyArray_STRIDE(out_array, 0) / PyArray_ITEMSIZE(out_array));
		py_dst = out;
		Py_INCREF(py_dst);
	}
	GKPY_TAKEREF(py_dst);

	// Decode the track data into the numpy array py_dst
	(*self->track)(c, PyArray_DATA(rcast<PyArrayObject*>(py_dst)), *dtype, stride);

	GKPY_FORGETREF(py_dst);
	return py_dst;
GKPY_OMETHOD_END

static PyObject* PyGenomeTrack_gtrack_version(PyObject* cls, PyObject* args)
{
	if (!PyArg_UnpackTuple(args, "", 0, 0))
		return nullptr;
	return PyInt_FromLong((long)genome_track::gtrack_version());
}

GKPY_METHOD_BEGIN_NOARG(GenomeTrack, close)
	self->track->close();
	GKPY_RETURN_NONE;
GKPY_METHOD_END

GKPY_METHODS_BEGIN(GenomeTrack)
	GKPY_METHOD_ENTRY(GenomeTrack, gtrack_version, METH_STATIC | METH_VARARGS, nullptr)
	GKPY_METHOD_ENTRY(GenomeTrack, close,          METH_NOARGS, nullptr)
GKPY_METHODS_END

GKPY_TYPEOBJ_BEGIN(GenomeTrack)
	tp_flags |= Py_TPFLAGS_BASETYPE | Py_TPFLAGS_HAVE_GC;
	tp_new = PyType_GenericNew;
	tp_init = PyGenomeTrack_Init;
	tp_dealloc = PyGenomeTrack_Dealloc;
	tp_getattro = PyGenomeTrack_GetAttro;
	tp_setattro = PyGenomeTrack_SetAttro;
	tp_traverse = PyGenomeTrack_Traverse;
	tp_clear = PyGenomeTrack_Clear;
	tp_call = PyGenomeTrack_Call;
	tp_methods = PyGenomeTrack_Methods;
GKPY_TYPEOBJ_END

END_NAMESPACE_GK
