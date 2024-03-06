/*
Copyright (C) 2016-2023 Deep Genomics Inc. All Rights Reserved.
*/
#include "py_util.h"
#include "strutil.h"
#include <vector>

BEGIN_NAMESPACE_GK

void get_nested_exception_what(std::string& out, const std::exception& e, int level)
{
	out.append(level, ' ');
	out += e.what();
	try {
		std::rethrow_if_nested(e);
	} catch (const std::exception& e) {
		out += '\n';
		get_nested_exception_what(out, e, level+2);
	}
}

void PyDisableGC(PyTypeObject* type)
{
	PyTypeObject* basetype = type->tp_base;
	type->tp_flags &= ~Py_TPFLAGS_HAVE_GC;
	type->tp_dealloc = basetype->tp_dealloc; // Change from PyType_GenericAlloc to base dealloc that matches base new/alloc
	type->tp_alloc = basetype->tp_alloc;     // Change from PyObject_GC_New to PyObject_New or whatever matches the directly-inherited tp_alloc slot
	type->tp_free = basetype->tp_free;       // Change from PyObject_GC_Del to PyObject_Del or whatever matches the directly-inherited tp_alloc slot
	PyType_Modified(type);
}

void PyForceNewGCInheritance(PyTypeObject* type)
{
	// Undo the damage that default Python 2.7/3.8+ inheritance does to the type method slots
	// https://bugs.python.org/issue41984
	// GC requires PyGC_Head header 16 byte overhead (larger than an interval itself) and incurs mark-sweep traversal costs
	PyTypeObject* basetype = type->tp_base;
	if (type->tp_basicsize == basetype->tp_basicsize && (type->tp_flags & Py_TPFLAGS_HAVE_GC) && (!(basetype->tp_flags & Py_TPFLAGS_HAVE_GC)))
		PyDisableGC(type);
}

void PyCheckSameBasicSize(PyTypeObject* type)
{
	if (type->tp_basicsize != type->tp_base->tp_basicsize)
		GK_THROW(runtime, "Subtype '{}' must have same tp_basictype as '{}'.", type->tp_name, type->tp_base->tp_name);
}

void PyDeleteMockAttrs(PyTypeObject* type)
{
	// If we're building docs, let sphinx import and see all the
	// mock methods and attributes on the original class.
	// An alternative is to use the :inherited-members: directive
	// in the .rst files, but then (a) the docstrings have to be
	// the C++, which makes the C++ ugle and is inconsistent with the way properties
	// are handled.
	static bool sphinxbuild = getenv("SPHINXBUILD") != nullptr;

	PyTypeObject* base = type->tp_base;
	if (!base || !type->tp_dict || sphinxbuild)
		return;

	// Collect the names of all attributes that have been marked as __mock__
	// by the @mock decorator. Once we have the list, we delete them all
	// (with exceptions, see below).
	vector<PyObject*> mock_names;
	PyObject* name;
	PyObject* attr;
	Py_ssize_t pos = 0;
	while (PyDict_Next(type->tp_dict, &pos, &name, &attr)) {

		if (PyObject_IsInstance(attr, (PyObject*)&PyProperty_Type)) {
			attr = PyObject_GetAttrString(attr, "fget"); // New reference
			Py_DECREF(attr); // We won't be holding on the new reference ourselves
		}

		if (PyObject_IsInstance(attr, (PyObject*)&PyStaticMethod_Type)) {
			attr = PyObject_GetAttrString(attr, "__func__"); // New reference
			Py_DECREF(attr); // We won't be holding on the new reference ourselves
		}

		if (PyObject_HasAttrString(attr, "__mock__")) {
			// Only delete the attribute if the base type ALSO has it.
			// If the base type doesn't have it, one of two situations has
			// occurred:
			//    1) The base type implements that attribute opaquely through the
			//       __getattribute__ mechanism, in which case we need to leave
			//       to mock attribute alone so that IDEs can find it through
			//       the dir() mechanism.
			//    2) Someone renamed the base type's attribute but forgot to
			//       update the mock attribute, in which case we'll still
			//       leave the mock attribute and if any code calls it there
			//       will be a RuntimeError via the mock_result mechanism.
			if (base && PyObject_HasAttr((PyObject*)base, name))
				mock_names.push_back(name);
		}
	}

	for (auto& mock_name : mock_names) PyDict_DelItem(type->tp_dict, mock_name);
}

PyObject* PyString_FromNonEmptyString(const char* str)
{
	GK_ASSERT(str);
	if (!*str)
		GKPY_RETURN_NONE;
	return PyString_FromString(str);
}

int gkpy_invalid_init(PyObject* self, PyObject*, PyObject*)
{
	auto s = PyObject_Str(PyObject_Type(self));
	PyErr_Format(PyExc_RuntimeError, "%s is an internal type and cannot be manually created.", PyString_AS_STRING(s));
	Py_DECREF(s);
	return -1;
}

END_NAMESPACE_GK
