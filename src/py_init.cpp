/*
Copyright (C) 2016-2023 Deep Genomics Inc. All Rights Reserved.
*/
#define NO_DISABLE_IMPORT_ARRAY // ensure that the global Array_API variable for import_array() is defined in this .cpp file
#include "py_util.h"

#include "genome_kit.h"
#include "py_genome.h"
#include "py_genome_anno.h"
#include "py_genome_dna.h"
#include "py_genome_track.h"
#include "py_interval.h"
#include "py_jralign.h"
#include "py_jrdist.h"
#include "py_ralign.h"
#include "py_variant_table.h"
#ifdef _MSC_VER
#pragma warning(disable : 4190 4297)  // disable warning about C linkage, and about throwing exceptions from C functions
#endif

#include <numpy/arrayobject.h>

BEGIN_NAMESPACE_GK

/////////////////////////////////////////////////

// A reference to the _gk_data.get_file Python object
// which we'll trigger each time C++ code asks for
// a datafile path to be resolved.
static PyObject* _py_resolve_datafile_path = nullptr;
static PyObject* _gk_data_file_not_found_error = nullptr;

static string traceback_format_exc_and_clear()
{
	// format_exc doesn't work from the C API, use format_exception
	PyObject *type, *value, *traceback;
	PyErr_Fetch(&type, &value, &traceback);
	GKPY_TAKEREF(type);
	GKPY_TAKEREF(value);
	GKPY_TAKEREF(traceback);

	PyAutoRef traceback_name{PyString_FromString("traceback")};
	GK_ASSERT(traceback_name);
	PyAutoRef traceback_module{PyImport_Import(traceback_name.get())};
	if (!traceback_module)
		return {};

	PyAutoRef format_exception{PyObject_GetAttrString(traceback_module.get(), "format_exception")};
	if (!format_exception)
		return {};

	PyAutoRef result{PyObject_CallFunctionObjArgs(format_exception.get(), type, value, traceback, nullptr)};
	if (!result)
		return {};

	string ret;
	for (Py_ssize_t i = 0, len = PyList_Size(result.get()); i < len; ++i) {
		ret += PyString_AS_STRING(PyList_GetItem(result.get(), i));
	}
	return ret;
}

static string resolve_datafile_path_with_python(string path)
{
	GK_ASSERT(_py_resolve_datafile_path);
	GK_ASSERT(_gk_data_file_not_found_error);

	// new_path_obj = _gk_data.get_file(path)
	PyObject* new_path_obj = PyObject_CallFunction(_py_resolve_datafile_path, "s", path.c_str());
	GKPY_TAKEREF(new_path_obj);

	// Check for error
	if (!new_path_obj) {
		// re-raise original error in case of GKDataFileNotFoundError, otherwise always raise ValueError
		if (_gk_data_file_not_found_error && PyErr_ExceptionMatches(_gk_data_file_not_found_error)) {
			GK_THROW(gk_data_file_not_found, "{}", traceback_format_exc_and_clear());
		} else {
			GK_THROW(value, "{}", traceback_format_exc_and_clear());
		}
	}

	// Check for wrong return type
	GK_CHECK(PyString_Check(new_path_obj) != 0, type,
			 "Expected string value from resolve_datafile_path, but received '{}'", Py_TYPE(new_path_obj)->tp_name);

	return PyString_AS_STRING(new_path_obj);
}

/////////////////////////////////////////////////

static const char py_register_doc[] =
"register(cls)\n--\n\n"
"Register a Python-defined type with the C++ backend.\n"
"\n"
"    Each time the ``genome_kit`` Python package defines a new class, it should notify the\n"
"    C++ extension module by using this as a decorator. This ensures three things:\n"
"\n"
"        1. The Python type follows Python 3.5 garbage collection inheritance rules,\n"
"           even in Python 2.7 (see py_util.cpp:PyForceNewGCInheritance).\n"
"        2. The Python type is available for C++ types to inherit from (see below).\n"
"        3. All mock properties, methods, and attributes are deleted at runtime.\n"
"    \n"
"    With regard to (2) above, registering the Python-defined type allows the C++ module\n"
"    to inherit its own classes from that type. For example, C++-defined type\n"
"    ``genome_kit._cxx.Exon`` inherits from Python-defined type ``genome_kit.Interval``.\n"
"    \n"
"    In more detail, here is the type inheritance::\n"
"\n"
"       A: ``genome_kit._cxx.Interval`` # C++-defined\n"
"       B: ``genome_kit.Interval``      # Python-defined, inherits from A\n"
"       C: ``genome_kit._cxx.Exon``     # C++-defined, inherits from B\n"
"       D: ``genome_kit.Exon``          # Python-defined, inherits from C\n"
"\n"
"    These types would be defined/registered in the following steps:\n"
"\n"
"       1. ``genome_kit._cxx.Interval`` is registered immediately when ``genome_kit._cxx`` is loaded.\n"
"       2. ``genome_kit.Interval`` is defined, inheriting from ``genome_kit._cxx.Interval``\n"
"       3. the decorator then calls ``genome_kit._cxx.register(genome_kit.Interval)``\n"
"       4. ``genome_kit._cxx.Exon`` is forced to inherit from ```genome_kit.Interval``\n"
"       5. ``genome_kit.Exon`` is eventually defined, inheriting from ``genome_kit._cxx.Exon``\n"
"\n"
"    TODO: A cleaner and more future-proofed design that avoids registration would \n"
"    be to break ``genome_kit._cxx`` into different extension modules, so that the\n"
"    core C++ module provides ``Interval``, and a downstream C++ module imports the\n"
"    Python ``interval`` module so that it may inherit from the Python-defined types\n"
"    when defining ``Exon`` etc in the first place.\n"
"    However, doing so will require figuring out the multi-extension build dependencies \n"
"    and source code arrangement, which has its own complexities but may need to be\n"
"    faced some day anyway.\n"
;

#define GKPY_REGISTER_CXX_DERIVED_BEGIN(name) \
	if (!strcmp(type->tp_name, #name)) { \
		Py##name::DefaultType = type; \
		PyForceNewGCInheritance(type);  /* Adopt newer (smarter) inheritance rules for GC tracking */ \
		PyDeleteMockAttrs(type);        /* Allow mock methods to fall through to C++ implementation */

#define GKPY_REGISTER_CXX_DERIVED_END \
		Py_INCREF(arg);  /* [ch1267] - although exact mechanism is unknown, returned value of register needs an incref to prevent segv when gc randomly triggers within import and finalize */ \
		return arg; \
	}

#define GKPY_REGISTER_CXX_DERIVED(name) \
	GKPY_REGISTER_CXX_DERIVED_BEGIN(name) \
	GKPY_REGISTER_CXX_DERIVED_END

static PyObject* py_register(PyObject*, PyObject* args)
{
	GKPY_TRY
	PyObject* arg;
	if (!PyArg_ParseTuple(args, "O", &arg))
		return nullptr;

	if (PyType_Check(arg)) {
		auto* type = (PyTypeObject*)arg;

		// No need to register this type, we just need a pointer for use in error type checking
		if (!strcmp(type->tp_name, "GKDataFileNotFoundError")) {
			_gk_data_file_not_found_error = arg;
			Py_INCREF(arg);
			return arg;
		}
		GKPY_REGISTER_CXX_DERIVED_BEGIN(Interval)
			PyCheckSameBasicSize(type);        // Don't allow Python-defined Interval to grow struct size
		GKPY_REGISTER_CXX_DERIVED_END
		GKPY_REGISTER_CXX_DERIVED(Genome)
		GKPY_REGISTER_CXX_DERIVED(GenomeDNA)
		GKPY_REGISTER_CXX_DERIVED(GenomeTrack)
		GKPY_REGISTER_CXX_DERIVED(GenomeTrackBuilder)
		GKPY_REGISTER_CXX_DERIVED(GenomeAnnotation)
		GKPY_REGISTER_CXX_DERIVED(Gene)
		GKPY_REGISTER_CXX_DERIVED(GeneTable)
		GKPY_REGISTER_CXX_DERIVED(Transcript)
		GKPY_REGISTER_CXX_DERIVED(TranscriptTable)
		GKPY_REGISTER_CXX_DERIVED(Exon)
		GKPY_REGISTER_CXX_DERIVED(ExonTable)
		GKPY_REGISTER_CXX_DERIVED(Intron)
		GKPY_REGISTER_CXX_DERIVED(IntronTable)
		GKPY_REGISTER_CXX_DERIVED(Cds)
		GKPY_REGISTER_CXX_DERIVED(CdsTable)
		GKPY_REGISTER_CXX_DERIVED(Utr)
		GKPY_REGISTER_CXX_DERIVED(UtrTable)
		GKPY_REGISTER_CXX_DERIVED_BEGIN(Variant)
			PyCheckSameBasicSize(type);        // Don't allow Python-defined Interval to grow struct size
		GKPY_REGISTER_CXX_DERIVED_END
		GKPY_REGISTER_CXX_DERIVED(VariantTable)
		GKPY_REGISTER_CXX_DERIVED_BEGIN(VCFVariant)
			PyCheckSameBasicSize(type);        // Don't allow Python-defined Interval to grow struct size
		GKPY_REGISTER_CXX_DERIVED_END
		GKPY_REGISTER_CXX_DERIVED(VCFTable)
		GKPY_REGISTER_CXX_DERIVED(JReadAlignments)
		GKPY_REGISTER_CXX_DERIVED(JunctionReadAlignmentsTable)
		GKPY_REGISTER_CXX_DERIVED(JunctionReadAlignments)
		GKPY_REGISTER_CXX_DERIVED(JunctionReadAlignment)
		GKPY_REGISTER_CXX_DERIVED(ReadDistributions)
		GKPY_REGISTER_CXX_DERIVED(JunctionReadDistributionTable)
		GKPY_REGISTER_CXX_DERIVED(JunctionReadDistribution)
		GKPY_REGISTER_CXX_DERIVED(JunctionReadCount)
		GKPY_REGISTER_CXX_DERIVED(ReadAlignments)
		GKPY_REGISTER_CXX_DERIVED(Junction)
		GKPY_REGISTER_CXX_DERIVED(JunctionTable)
		GKPY_REGISTER_CXX_DERIVED(Alignment)
		GKPY_REGISTER_CXX_DERIVED(AlignmentTable)
		GKPY_REGISTER_CXX_DERIVED(AlignmentMatch)
		GKPY_REGISTER_CXX_DERIVED(AlignmentMatchTable)
		// <--- INSERT CASES FOR NEW C++-DEFINED TYPES AND SUBTYPES HERE

		if (PyType_IsSubtype(type, PyInterval::DefaultType)) {
			PyForceNewGCInheritance(type);
			Py_INCREF(arg); // [ch1267] - although exact mechanism is unknown, returned value of register needs an incref to prevent segv when gc randomly triggers within import and finalize
			return arg;
		}

		print("WARNING: genome_kit._cxx.register did not expect type '{}' to be registered.\n", type->tp_name);

	} else if (PyCallable_Check(arg)) {
		auto py_name = PyObject_GetAttrString(arg, "__name__");
		GKPY_TAKEREF(py_name);
		const char* name = PyString_AsString(py_name);
		if (!strcmp(name, "resolve_datafile_path")) {
			// Initialize datafile callback to be _gk_data.get_file
			resolve_datafile_path = resolve_datafile_path_with_python;
			Py_INCREF(arg);
			PyObject* tmp = _py_resolve_datafile_path;
			_py_resolve_datafile_path = arg;
			if (tmp)
				Py_XDECREF(tmp);
			return arg;
		}
		print("WARNING: genome_kit._cxx.register did not expect callable '{}' to be registered.\n", name);
	}

	Py_INCREF(arg); // [ch1267] - although exact mechanism is unknown, returned value of register needs an incref to prevent segv when gc randomly triggers within import and finalize
	return arg;
	GKPY_CATCH_RETURN_NULL
}

/////////////////////////////////////////////////////////////////////////

static const char py_no_gc_doc[] =
"no_gc(cls)\n--\n\n"
"Explicitly mark a type as not needing to be garbage collected.\n"
"\n"
"If a type contains a fixed set of attributes and those attributes\n"
"can never form a reference cycle under correct usage, then that type\n"
"can be safely decorated as ``@no_gc``.\n"
"\n"
"This is only important for types for which there may be enough\n"
"instances created to either slow down the garbage collector or.\n"
"waste significant per-instance memory from the GC header.\n"
;

static PyObject* py_no_gc(PyObject*, PyObject* args)
{
	GKPY_TRY
	PyTypeObject* type;
	if (!PyArg_ParseTuple(args, "O!", &PyType_Type, &type))
		return nullptr;

	PyDisableGC(type);

	return (PyObject*)type;
	GKPY_CATCH_RETURN_NULL
}

/////////////////////////////////////////////////////////////////////////

static struct PyMethodDef py_exports[] = {
	{ "register", (PyCFunction)py_register, METH_VARARGS, py_register_doc },
	{ "no_gc", (PyCFunction)py_no_gc, METH_VARARGS, py_no_gc_doc },
	{ nullptr, nullptr, 0, nullptr } // sentinel
};

bool try_relative_import(PyObject* module, char* module_name)
{
	return PyImport_ImportModuleLevel(module_name, PyModule_GetDict(module), Py_None, Py_None, 1) != nullptr;
}

END_NAMESPACE_GK

USING_NAMESPACE_GK

static int py_mod_exec(PyObject* module)
{
#ifdef NPY_1_7_API_VERSION
	auto ret = _import_array(); // Initialize numpy
	if (ret != 0)
		return ret;
#endif

	GKPY_TRY
	// Initialize genome_kit globals (string tables etc)
	init_genome_kit();

	// Interval
	Init_Interval_PyStrings();
	PyInterval::Register(module);

	// Now that we've added _cxx.Interval to the module,
	// force import of the Python counterparts so that they may do two things:
	//  1) Inherit from the C++-defined types, and
	//  2) Register the Python-defined types so that (among other things)
	//     the C++ code knows the default type to construct for Interval
	if (!try_relative_import(module, "interval")) return -1;

	// Genome, GenomeDNA
	PyGenome::Register(module);
	PyGenomeDNA::Register(module);
	PyGenomeTrack::Register(module);
	PyGenomeTrackBuilder::Register(module);

	// GenomeAnno and table types, but not row types since they subclass Python-defined Intron type
	Init_GenomeAnno_PyStrings();
	PyGenomeAnno::Register(module);

	{
		PyGene::Register(module);
		PyGeneTable::Register(module);

		PyTran::Register(module);
		PyTranTable::Register(module);

		PyExon::Register(module);
		PyExonTable::Register(module);

		PyIntr::Register(module);
		PyIntrTable::Register(module);

		PyCds::Register(module);
		PyCdsTable::Register(module);

		PyUtr::Register(module);
		PyUtrTable::Register(module);
		// <--- REGISTER NEW ANNOTATION TABLE TYPES HERE
	}

	// Variant
	PyVariant::Register(module);
	PyVariantTable::Register(module);

	// Now that we've added _cxx.Variant to the module,
	// force import of the Python counterparts, like for interval above.
	if (!try_relative_import(module, "variant")) return -1;

	// VCFVariant
	PyVCFVariant::Register(module);
	PyVCFTable::Register(module);

	// Read alignments
	PyJReadAlignments::Register(module);
	PyJRAlignsTable::Register(module);
	PyJRAligns::Register(module);
	PyJRAlign::Register(module);

	// Read distribution
	PyReadDistributions::Register(module);
	PyJRDistTable::Register(module);
	PyJRDist::Register(module);
	PyJRCount::Register(module);

	// Junction read alignments
	PyReadAlignments::Register(module);
	PyJunction::Register(module);
	PyJunctionTable::Register(module);
	PyAlignment::Register(module);
	PyAlignmentTable::Register(module);
	PyAlignmentMatch::Register(module);
	PyAlignmentMatchTable::Register(module);

	// <--- INSERT NEW C++-DEFINED TYPES (UNLESS THEY INHERIT FROM PYTHON-DEFINED TYPES -- SEE py_register)
	return 0;
	GKPY_CATCH_RETURN_VALUE(-1);
}

static int py_copy_spec_to_module(PyObject* spec, PyObject* moddict, const char* from_name, const char* to_name)
{
	auto value  = PyObject_GetAttrString(spec, from_name);
	auto result = 0;
	if (value) {
		result = PyDict_SetItemString(moddict, to_name, value);
	} else if (PyErr_ExceptionMatches(PyExc_AttributeError)) {
		PyErr_Clear();
	} else {
		result = -1;
	}
	return result;
}
static PyObject* py_mod_create(PyObject* spec, PyModuleDef*)
{
	auto modname = PyObject_GetAttrString(spec, "name");
	if (!modname)
		return nullptr;
	GKPY_TAKEREF(modname);

	auto module = PyModule_NewObject(modname);
	if (!module)
		return nullptr;
	GKPY_TAKEREF(module);

	auto moddict = PyModule_GetDict(module);
	if (py_copy_spec_to_module(spec, moddict, "loader", "__loader__") < 0)
		return nullptr;
	if (py_copy_spec_to_module(spec, moddict, "origin", "__file__") < 0)
		return nullptr;
	if (py_copy_spec_to_module(spec, moddict, "parent", "__package__") < 0)
		return nullptr;
	if (py_copy_spec_to_module(spec, moddict, "submodule_search_locations", "__path__") < 0)
		return nullptr;
	GKPY_FORGETREF(module);
	return module;
}

extern "C" {

EX_PYTHON_INIT PyObject* PyInit__cxx(void)
{
	static struct PyModuleDef_Slot py_slots[]
		= { { Py_mod_create, (void*)py_mod_create }, { Py_mod_exec, (void*)py_mod_exec }, { 0, nullptr } };
	static struct PyModuleDef moduledef
		= { PyModuleDef_HEAD_INIT, "_cxx", "GenomeKit C++ extension module", 0, py_exports, py_slots };
	return PyModuleDef_Init(&moduledef);
}

} // extern "C"
