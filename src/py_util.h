/*
Copyright (C) 2016-2023 Deep Genomics Inc. All Rights Reserved.
*/
#pragma once
#ifndef __GK_PY_UTIL_H__
#define __GK_PY_UTIL_H__

#ifdef _DEBUG
#define WAS_DEBUG
#undef _DEBUG
#endif
#include <Python.h>

#define PyInt_AsLong PyLong_AsLong
#define PyInt_AS_LONG PyLong_AS_LONG
#define PyInt_Check PyLong_Check
#define PyInt_FromLong PyLong_FromLong
#define PyInt_FromSize_t PyLong_FromSize_t
#define PyInt_Type PyLong_Type
#define PyString_AsString PyUnicode_AsUTF8
#define PyString_AS_STRING PyUnicode_AsUTF8
#define PyString_Check PyUnicode_Check
#define PyString_FromString(v) PyUnicode_DecodeUTF8(v, strlen(v), nullptr)
#define PyString_FromStringAndSize(v, l) PyUnicode_DecodeUTF8(v, l, nullptr)

template <class T>
auto PyString_FromSV(T s) {
	return PyString_FromStringAndSize(s.data(), s.size());
}

// numpy definitions for extensions that have numpy API
#define NPY_NO_DEPRECATED_API NPY_1_7_API_VERSION
#define PY_ARRAY_UNIQUE_SYMBOL GenomeKit_Array_API
#ifndef NO_DISABLE_IMPORT_ARRAY
#define NO_IMPORT_ARRAY
#endif

#ifdef WAS_DEBUG
#define _DEBUG
#endif

#include "defines.h"
#include "gk_assert.h"
#include "util.h"
#include <cstring>
#include <functional>
#include <memory>

#ifndef GKPY_LIBNAME
#error Must define GKPY_LIBNAME so that the module knows the name of the Python package it iss part of
#endif

#define _GKPY_LIBNAME_STR(name) GK_EXPAND_STR(name)
#define GKPY_LIBNAME_STR        _GKPY_LIBNAME_STR(GKPY_LIBNAME)

#ifdef GKPY_TRACE_MEM
#include <cstdio>
#define GKPY_TRACE_CTOR(self, name) print("{}: {}\n", self, #name)
#define GKPY_TRACE_DTOR(self, name) gkpy_trace_dtor __dtor_notify(#name, self)
struct gkpy_trace_dtor {
	const char* name;
	void* self;
	gkpy_trace_dtor(const char* name, void* self): name(name), self(self) { }
	~gkpy_trace_dtor() { print("{}: ~{}\n", self, name); }
};
#else
#define GKPY_TRACE_CTOR(self, name)
#define GKPY_TRACE_DTOR(self, name)
#endif

BEGIN_NAMESPACE_GK

////////////////////////////////////////////////////////////////////

void get_nested_exception_what(std::string& out, const std::exception& e, int level=0);

const char* gkpy_strip_exception_type(const char* msg);

// The GKPY_TRY mechanism converts C++ exceptions (and C structured exceptions)
// to a Python exception and then sets the global Python exception state.
// Once an exception is caught like this, it's up to the calling code to
// return control back to the interpreter without raising further exceptions
// so that the original exception is the one seen by Python.
// As such, if an internal C++ function calls a Python function (or an
// internal C++ function with its own GKPY_TRY/GKPY_CATCH) then the calling
// function must check for an error return value (NULL or -1) and return
// error state up the call stack, rather than throw a new exception.
//
#define GKPY_TRY try {
#define GKPY_CATCH_CASE(cpp_etype, py_etype, return_stmt) \
                         } catch (const cpp_etype& e) {\
                             std::string what; \
                             get_nested_exception_what(what, e); \
                             PyErr_SetString(py_etype, what.c_str()); \
                             return_stmt;

#define _GKPY_CATCH(return_stmt) \
	GKPY_CATCH_CASE(gk::assertion_error, PyExc_AssertionError, return_stmt) \
	GKPY_CATCH_CASE(gk::file_error, PyExc_OSError, return_stmt) \
	GKPY_CATCH_CASE(gk::type_error, PyExc_TypeError, return_stmt) \
	GKPY_CATCH_CASE(gk::value_error, PyExc_ValueError, return_stmt) \
	GKPY_CATCH_CASE(gk::index_error, PyExc_IndexError, return_stmt) \
	GKPY_CATCH_CASE(gk::key_error, PyExc_KeyError, return_stmt) \
	GKPY_CATCH_CASE(gk::memory_error, PyExc_MemoryError, return_stmt) \
	GKPY_CATCH_CASE(gk::not_implemented_error, PyExc_NotImplementedError, return_stmt) \
	GKPY_CATCH_CASE(gk::unreachable_code_error, PyExc_RuntimeError, return_stmt) \
	GKPY_CATCH_CASE(gk::runtime_error, PyExc_RuntimeError, return_stmt) \
	GKPY_CATCH_CASE(std::exception, PyExc_RuntimeError, return_stmt) \
    }\
	return_stmt;

#define GKPY_CATCH_RETURN_VALUE(x) _GKPY_CATCH(return x)
#define GKPY_CATCH_RETURN_NULL     _GKPY_CATCH(return NULL)
#define GKPY_CATCH                 _GKPY_CATCH(return)

////////////////////////////////////////////////////////////////////////

#define GKPY_RETURN_NONE  do { Py_INCREF(Py_None);  return Py_None; } while (0)
#define GKPY_RETURN_TRUE  do { Py_INCREF(Py_True);  return Py_True; } while (0)
#define GKPY_RETURN_FALSE do { Py_INCREF(Py_False); return Py_False; } while (0)
#define GKPY_RETURN_BOOL(cond) return PyReturnBool(cond)
#define GKPY_RETURN_INCREF(value) do { PyObject* tmp_value = value; Py_INCREF(tmp_value); return tmp_value; } while (0)

INLINE PyObject* PyReturnBool(bool cond)
{
	PyObject* r = cond ? Py_True : Py_False;
	Py_INCREF(r);
	return r;
}

// Automatically returns NULL or -1 depending on the expected
// return value of the function (PyObject* or int).
#define GKPY_RETURN_ERROR \
	do { PyImplicitError error; return error; } while (0)

struct PyImplicitError {
	INLINE operator PyTypeObject*() const { return nullptr; }
	INLINE operator PyObject*() const { return nullptr; }
	INLINE operator int() const { return -1; }
};

////////////////////////////////////////////////////////////////////

// Class to make Py_INCREF/Py_DECREF exception safe, when an error occurs and
// the stack needs to be unwound all the way to the last Python-to-C entry point.
struct PyObjectDecrementer {
	void operator()(PyObject* obj) const noexcept {
#if PY_VERSION_HEX >= 0x030B0000 // >= py311
		// see https://github.com/deepgenomics/GenomeKit/issues/102
		// and https://github.com/python/cpython/issues/126508
		if (PyThreadState_Get() != nullptr) {
			Py_XDECREF(obj);
		}
#else
		// on < py311 calling PyThreadState_Get() on exit causes:
		// Fatal Python error: PyThreadState_Get: the function must be called with the GIL held, but the GIL is released (the current Python thread state is NULL)
		Py_XDECREF(obj);
#endif
	};
};
using PyAutoRef = std::unique_ptr<PyObject, PyObjectDecrementer>;

// Use when need to call Py_INCREF on acquire
#define GKPY_MAKEREF(obj_ptr)    Py_INCREF(obj_ptr); PyAutoRef _##obj_ptr##_ref(obj_ptr);

// Use when there's no need to call Py_INCREF on acquire
#define GKPY_TAKEREF(obj_ptr)  PyAutoRef _##obj_ptr##_ref(obj_ptr);

// Use when there's no longer a need to call Py_DECREF, i.e. the responsibility is
// about to be passed on to another function, or to to Python interpreter.
#define GKPY_FORGETREF(obj_ptr)  _##obj_ptr##_ref.release();

// Use when need to call Py_INCREF on acquire
#define GKPY_BUILDVALUE_TEMP(name, fmt, ...) \
	PyObject* name = Py_BuildValue(fmt, __VA_ARGS__); \
	GKPY_TAKEREF(name);

// Use when need to call Py_INCREF on acquire
#define GKPY_BUILDVALUE_ARG(fmt, ...)    PyAutoRef(Py_BuildValue(fmt, __VA_ARGS__), false).get()

/////////////////////////////////////////////////////////////////////////////////////////////

#define GKPY_SUBTYPE_BEGIN(name, basename) \
	struct Py##name: public Py##basename { \
		using base_type = Py##basename; /* We'll inherit from base_type::DefaultType, not base_type::Type */ \
		static PyTypeObject* Type; /* C++-defined Python type */ \
		static PyTypeObject* DefaultType; /* Default type to create, may be Python-defined */ \
		static void EnsureInit(); /* Make sure Init() called at least once */ \
		static void Init(); /* Initialize 'Type' members, if not already initialized */ \
		static void Register(PyObject* module);

#define GKPY_SUBTYPE_END \
	};

#define GKPY_TYPE_BEGIN(name) GKPY_SUBTYPE_BEGIN(name, Object)
#define GKPY_TYPE_END         GKPY_SUBTYPE_END

struct PyAsPtrSource : public PyObject
{
	// use fn_ptr for dispatch since pure virtual will segfault in preprocessor macro without message
	using validate_function = void (*)(const PyAsPtrSource *); // throw if invalid
	validate_function validator{};

	void set_validator(validate_function validator) { this->validator = validator; }
	void validate() const { GK_ASSERT(validator != nullptr); (*validator)(this); }
};
#define GKPY_SOURCE_TYPE_BEGIN(name) GKPY_SUBTYPE_BEGIN(name, AsPtrSource)
#define GKPY_SOURCE_TYPE_END         GKPY_SUBTYPE_END

inline void add_fqn(PyHeapTypeObject* obj)
{
	obj->ht_qualname = obj->ht_name;
	Py_INCREF(obj->ht_qualname);
}

#define _GKPY_TYPEOBJ_BEGIN(tdecl, name, basicsize) \
	tdecl PyTypeObject* Py##name::Type = 0;\
	tdecl PyTypeObject* Py##name::DefaultType = 0;\
	tdecl void Py##name::Register(PyObject* module)\
	{\
		EnsureInit(); \
		Py_INCREF((PyObject*)Type);\
		PyModule_AddObject(module, #name, (PyObject*)Type);\
	}\
	tdecl void Py##name::EnsureInit()\
	{\
		if (!Type || !Py_TYPE(Type)) /* May have been previously initialized by a subtype's initializer */ \
			Init(); \
	}\
	tdecl void Py##name::Init()\
	{\
		bool allow_smaller_basicsize = false;\
		const char* tp_name = GKPY_LIBNAME_STR "._cxx." #name;\
		Py_ssize_t tp_basicsize = basicsize;\
		Py_ssize_t tp_itemsize = 0;\
		destructor tp_dealloc = 0;\
		getattrfunc tp_getattr = 0;\
		setattrfunc tp_setattr = 0;\
		reprfunc tp_repr = 0;\
		hashfunc tp_hash = 0;\
		ternaryfunc tp_call = 0;\
		reprfunc tp_str = 0;\
		getattrofunc tp_getattro = 0;\
		setattrofunc tp_setattro = 0;\
		PyBufferProcs* tp_as_buffer = 0;\
		long tp_flags = Py_TPFLAGS_DEFAULT;\
		const char* tp_doc = 0;\
		traverseproc tp_traverse = 0;\
		inquiry tp_clear = 0;\
		richcmpfunc tp_richcompare = 0;\
		Py_ssize_t tp_weaklistoffset = 0;\
		getiterfunc tp_iter = 0;\
		iternextfunc tp_iternext = 0;\
		PyMethodDef* tp_methods = 0;\
		PyMemberDef* tp_members = 0;\
		PyGetSetDef* tp_getset = 0;\
		PyTypeObject* tp_base = 0;\
		PyObject* tp_dict = 0;\
		descrgetfunc tp_descr_get = 0;\
		descrsetfunc tp_descr_set = 0;\
		Py_ssize_t tp_dictoffset = 0;\
		initproc tp_init = 0;\
		allocfunc tp_alloc = 0;\
		newfunc tp_new = 0;\
		freefunc tp_free = 0;\
		inquiry tp_is_gc = 0;\
		PyObject* tp_bases = 0;\
		PyObject* tp_mro = 0;\
		PyObject* tp_cache = 0;\
		PyObject* tp_subclasses = 0;\
		PyObject* tp_weaklist = 0;\
		destructor tp_del = 0;\
		\
		lenfunc sq_length = 0;\
		binaryfunc sq_concat = 0;\
		ssizeargfunc sq_repeat = 0;\
		ssizeargfunc sq_item = 0;\
		ssizeobjargproc sq_ass_item = 0;\
		objobjproc sq_contains = 0;\
		binaryfunc sq_inplace_concat = 0;\
		ssizeargfunc sq_inplace_repeat = 0;\
		\
		lenfunc mp_length = 0;\
		binaryfunc mp_subscript = 0;\
		objobjargproc mp_ass_subscript = 0;\
		{

#define GKPY_TYPEOBJ_END \
		}\
		if (tp_base && !allow_smaller_basicsize) \
			GK_CHECK(tp_base->tp_basicsize <= tp_basicsize, runtime, \
			         "subtype basicsize ({}) was less than basetype basicsize ({}); something went wrong", \
			         tp_basicsize, tp_base->tp_basicsize); \
		\
		PyHeapTypeObject* heap_type = nullptr; \
		if (!tp_base || !PyType_HasFeature(tp_base, Py_TPFLAGS_HEAPTYPE))\
		{ \
			static PyTypeObject StaticType = { PyVarObject_HEAD_INIT(nullptr, 0) };\
			Type = &StaticType;\
		} else {\
			heap_type = (PyHeapTypeObject*)PyType_GenericAlloc(&PyType_Type, 0);\
			GK_ASSERT(heap_type);\
			Type = &heap_type->ht_type;\
			tp_flags |= Py_TPFLAGS_HEAPTYPE;\
		}\
		DefaultType = Type;\
		\
		auto& _typeobj = *Type;\
		_typeobj.tp_name = tp_name;\
		_typeobj.tp_basicsize = tp_basicsize;\
		_typeobj.tp_itemsize = tp_itemsize;\
		_typeobj.tp_dealloc = tp_dealloc;\
		_typeobj.tp_getattr = tp_getattr;\
		_typeobj.tp_setattr = tp_setattr;\
		_typeobj.tp_repr = tp_repr;\
		_typeobj.tp_hash = tp_hash;\
		_typeobj.tp_call = tp_call;\
		_typeobj.tp_str = tp_str;\
		_typeobj.tp_getattro = tp_getattro;\
		_typeobj.tp_setattro = tp_setattro;\
		_typeobj.tp_as_buffer = tp_as_buffer;\
		_typeobj.tp_flags = tp_flags;\
		_typeobj.tp_doc = tp_doc;\
		_typeobj.tp_traverse = tp_traverse;\
		_typeobj.tp_clear = tp_clear;\
		_typeobj.tp_richcompare = tp_richcompare;\
		_typeobj.tp_weaklistoffset = tp_weaklistoffset;\
		_typeobj.tp_iter = tp_iter;\
		_typeobj.tp_iternext = tp_iternext;\
		_typeobj.tp_methods = tp_methods;\
		_typeobj.tp_members = tp_members;\
		_typeobj.tp_getset = tp_getset;\
		if (heap_type)\
			Py_INCREF(tp_base);\
		_typeobj.tp_base = tp_base;\
		_typeobj.tp_dict = tp_dict;\
		_typeobj.tp_descr_get = tp_descr_get;\
		_typeobj.tp_descr_set = tp_descr_set;\
		_typeobj.tp_dictoffset = tp_dictoffset;\
		_typeobj.tp_init = tp_init;\
		_typeobj.tp_alloc = tp_alloc;\
		_typeobj.tp_new = tp_new;\
		_typeobj.tp_free = tp_free;\
		_typeobj.tp_is_gc = tp_is_gc;\
		_typeobj.tp_bases = tp_bases;\
		_typeobj.tp_mro = tp_mro;\
		_typeobj.tp_cache = tp_cache;\
		_typeobj.tp_subclasses = tp_subclasses;\
		_typeobj.tp_weaklist = tp_weaklist;\
		_typeobj.tp_del = tp_del;\
		\
		static PySequenceMethods tp_as_sequence = {\
			sq_length,\
			sq_concat,\
			sq_repeat,\
			sq_item,\
			0,\
			sq_ass_item,\
			0,\
			sq_contains,\
			sq_inplace_concat,\
			sq_inplace_repeat\
		};\
		_typeobj.tp_as_sequence = (sq_length || sq_concat || sq_repeat || sq_item || sq_ass_item || sq_contains || sq_inplace_concat || sq_inplace_repeat) ? &tp_as_sequence : 0;\
		\
		static PyMappingMethods tp_as_mapping = {\
			mp_length,\
			mp_subscript,\
			mp_ass_subscript\
		};\
		_typeobj.tp_as_mapping = (mp_length || mp_subscript || mp_ass_subscript) ? &tp_as_mapping : 0;\
		\
		if (heap_type) {\
			const char* s = strrchr(tp_name, '.');\
			s = s ? s + 1 : tp_name;\
			heap_type->ht_name = PyString_FromString(s);\
			GK_ASSERT(heap_type->ht_name);\
			add_fqn(heap_type);\
			\
			/* see type_dealloc for corresponding PyObject_Free */\
			if (tp_doc) {\
				auto len = strlen(tp_doc);\
				_typeobj.tp_doc = (char*)PyObject_MALLOC(len);\
				GK_ASSERT(_typeobj.tp_doc);\
				std::memcpy((void*)_typeobj.tp_doc, tp_doc, len);\
			}\
		}\
		\
		PyType_Ready(&_typeobj);\
		if (tp_base) \
			PyForceNewGCInheritance(&_typeobj); /* Adopt newer (smarter) inheritance rules for GC tracking */ \
		if (heap_type) {\
			/* TODO: ht_cached_keys */\
			\
			const char* s = strrchr(tp_name, '.');\
			if (s) {\
				auto modname = PyString_FromStringAndSize(tp_name, s - tp_name);\
				GK_ASSERT(modname);\
				auto err = PyDict_SetItemString(_typeobj.tp_dict, "__module__", modname);\
				Py_DECREF(modname);\
				GK_ASSERT(err == 0);\
			}\
		}\
	}

// Some utility classes to compute sizes based on GKPY_VALUE_TYPE struct layouts.
struct __value_type_as_ptr: public PyObject {
	void* as_ptr; // This field will be non-NULL
	void* source; // This field will be non-NULL
};

template <typename T>
struct __value_type_as_value: public PyObject {
	void* as_ptr;  // This field will be NULL
	char value[sizeof(T)]; // This will contain the object data
};

template <typename T>
struct __value_type_sizes {
	enum { as_value = sizeof(__value_type_as_value<T>),
	       as_ptr = sizeof(__value_type_as_ptr),
	       as_value_or_ptr = as_value > as_ptr ? as_value : as_ptr };
};

#define GKPY_TEMPLATE_TYPEOBJ_BEGIN(tdecl, name)    _GKPY_TYPEOBJ_BEGIN(tdecl, name, sizeof(Py##name))
#define GKPY_TEMPLATE_SUBTYPEOBJ_BEGIN(tdecl, name) _GKPY_TYPEOBJ_BEGIN(tdecl, name, sizeof(Py##name)) \
                                                        base_type::EnsureInit(); \
                                                        tp_base = base_type::DefaultType;
#define GKPY_TYPEOBJ_BEGIN(name)                    _GKPY_TYPEOBJ_BEGIN(GK_EXPAND_EMPTY, name, sizeof(Py##name))
#define GKPY_SUBTYPEOBJ_BEGIN(name)                 _GKPY_TYPEOBJ_BEGIN(GK_EXPAND_EMPTY, name, sizeof(Py##name))  \
                                                        base_type::EnsureInit(); \
                                                        tp_base = base_type::DefaultType;
#define GKPY_VALUE_TYPEOBJ_BEGIN(name, modes)       _GKPY_TYPEOBJ_BEGIN(GK_EXPAND_EMPTY, name, __value_type_sizes<Py##name::value_t>::modes)
#define GKPY_VALUE_SUBTYPEOBJ_BEGIN(name, modes)    _GKPY_TYPEOBJ_BEGIN(GK_EXPAND_EMPTY, name, __value_type_sizes<Py##name::value_t>::modes) \
                                                        base_type::EnsureInit(); \
														tp_base = base_type::DefaultType; \
                                                        allow_smaller_basicsize = (0==strcmp(#modes, "as_ptr")); // as_ptr mode allows subtypes to be smaller than base types
#define GKPY_VALUE_SUBTYPEOBJ_BEGIN2(name, modes)   _GKPY_TYPEOBJ_BEGIN(GK_EXPAND_EMPTY, name, __value_type_sizes<Py##name::pyconstructed_t>::modes) \
                                                        base_type::EnsureInit(); \
														tp_base = base_type::DefaultType; \
                                                        allow_smaller_basicsize = (0==strcmp(#modes, "as_ptr")); // as_ptr mode allows subtypes to be smaller than base types

#define GKPY_TEMPLATE_TYPEOBJ_END    GKPY_TYPEOBJ_END
#define GKPY_TEMPLATE_SUBTYPEOBJ_END GKPY_TYPEOBJ_END
#define GKPY_SUBTYPEOBJ_END          GKPY_TYPEOBJ_END
#define GKPY_VALUE_TYPEOBJ_END       GKPY_TYPEOBJ_END
#define GKPY_VALUE_SUBTYPEOBJ_END    GKPY_TYPEOBJ_END

//////////////////////////////////////////////////////////////

#define GKPY_NEW_BEGIN(name)\
	PyObject* Py##name##_New(PyTypeObject* type, PyObject* args, PyObject* kwds)\
	{\
		GKPY_TRY\
		PyObject* selfo = type->tp_alloc(type, 0); /* If Py_TPFLAGS_HAVE_GC is set, this calls GC_Track(self) so we don't have to */ \
		Py##name* self = (Py##name*)selfo; \
		GKPY_TRACE_CTOR(self, name);\
		GKPY_TAKEREF(selfo);

#define GKPY_NEW_END \
	GKPY_FORGETREF(selfo);\
		return selfo;\
	GKPY_CATCH_RETURN_NULL \
	}

#define GKPY_DEALLOC_BEGIN(name) \
	void Py##name##_Dealloc(PyObject* selfo) \
	{ \
		GKPY_TRY \
		[[maybe_unused]] Py##name* self = *((Py##name**)&selfo); \
		GKPY_TRACE_DTOR(self, name); \
		if (PyObject_IS_GC(selfo)) \
			Py_TYPE(selfo)->tp_clear(selfo);

#define GKPY_DEALLOC_END \
		Py_TYPE(selfo)->tp_free(selfo);  /* Because Py_TPFLAGS_HAVE_GC is set, this calls GC_UnTrack(self) so we don't have to */ \
	GKPY_CATCH \
	}

/////////////////////////////////////////////////////////////////

#define GKPY_NEW_OWNED_BEGIN(name, owner_name) \
	GKPY_NEW_BEGIN(name) \
		self->owner = 0; \
		if (!PyArg_ParseTuple(args, "O!", Py##owner_name::DefaultType, &self->owner)) \
			return NULL; \
		Py_INCREF(self->owner);

#define GKPY_NEW_OWNED_END GKPY_NEW_END

#define GKPY_TRAVERSE_OWNED_BEGIN(name) \
	GKPY_TRAVERSE_BEGIN(name) \
		GKPY_VISIT(owner);

#define GKPY_TRAVERSE_OWNED_END GKPY_TRAVERSE_END

#define GKPY_CLEAR_OWNED_BEGIN(name) \
	GKPY_CLEAR_BEGIN(name) \
		/* GKPY_CLEAR(owner) -- Deliberately do not clear source until dealloc. Ensures source deallocated last. */

#define GKPY_CLEAR_OWNED_END GKPY_CLEAR_END

// Same as GKPY_DEALLOC_BEGIN but scoped to ensure GKPY_CLEAR(source) and GKPY_TRACE_DTOR
// are invoked after the owned object is both cleared and freed.
#define GKPY_DEALLOC_OWNED_BEGIN(name) \
	void Py##name##_Dealloc(PyObject* selfo)\
	{ \
		GKPY_TRY \
		[[maybe_unused]] Py##name* self = (Py##name*)selfo; \
		PyObject* __owner = self->owner; \
		{ \
			GKPY_TRACE_DTOR(self, name); \
			if (PyObject_IS_GC(selfo)) \
				Py_TYPE(selfo)->tp_clear(selfo);

#define GKPY_DEALLOC_OWNED_END \
			Py_TYPE(selfo)->tp_free(selfo);  /* Because Py_TPFLAGS_HAVE_GC is set, this calls GC_UnTrack(self) so we don't have to */ \
		} \
		Py_XDECREF(__owner); \
	GKPY_CATCH \
	}

#define GKPY_GETATTRO_BEGIN(name)\
	PyObject* Py##name##_GetAttro(PyObject* selfo, PyObject* attro)\
	{\
		GKPY_TRY\
		[[maybe_unused]] Py##name* self = *((Py##name**)&selfo);\
		[[maybe_unused]] const char* attr = PyString_AS_STRING(attro);

#define GKPY_GETATTR_CASE(name) if (!strcmp(attr, name))
#define GKPY_RETURN_ATTRO(attr) do { Py_INCREF(self->attr); return self->attr; } while (0)

#define GKPY_GETATTRO_END\
	GKPY_CATCH_RETURN_NULL\
	}

#define GKPY_SETATTRO_BEGIN(name)\
	int Py##name##_SetAttro(PyObject* selfo, PyObject* attro, PyObject* value)\
	{\
		GKPY_TRY\
		[[maybe_unused]] Py##name* self = *((Py##name**)&selfo);\
		[[maybe_unused]] const char* attr = PyString_AS_STRING(attro);

#define GKPY_SETATTR_READONLY(name) GK_CHECK(strcmp(attr, name), type, "Cannot set read-only attribute '{}' on object '{}'", name, Py_TYPE(selfo)->tp_name)

#define GKPY_SETATTRO_END\
	GKPY_CATCH_RETURN_VALUE(-1)\
	}

#define GKPY_TRAVERSE_BEGIN_COMMON(name)\
	int Py##name##_Traverse(PyObject* selfo, visitproc visit, void* arg)\
	{ \
		GKPY_TRY\
		[[maybe_unused]] Py##name* self = *((Py##name**)&selfo);

// required for Python >=3.9 (see bpo-35810, bpo-40217)
// https://docs.python.org/3/whatsnew/3.9.html#changes-in-the-c-api
#if PY_VERSION_HEX >= 0x03090000
	#define GKPY_TRAVERSE_BEGIN(name) \
	GKPY_TRAVERSE_BEGIN_COMMON(name)\
	Py_VISIT(Py_TYPE(selfo));
#else
#define GKPY_TRAVERSE_BEGIN(name) \
	GKPY_TRAVERSE_BEGIN_COMMON(name)
#endif


#define GKPY_VISIT(name)         Py_VISIT(self->name)

#define GKPY_TRAVERSE_END\
		return 0;\
	GKPY_CATCH_RETURN_VALUE(-1)\
	}

#define GKPY_CLEAR_BEGIN(name)\
	int Py##name##_Clear(PyObject* selfo)\
	{\
		GKPY_TRY\
		[[maybe_unused]] Py##name* self = *((Py##name**)&selfo);

#define GKPY_CLEAR(name)\
	do {\
		PyObject* _tmp = (PyObject*)self->name;\
		self->name = NULL;\
		Py_XDECREF(_tmp);\
	} while (0)\

#define GKPY_CLEAR_END\
		return 0;\
	GKPY_CATCH_RETURN_VALUE(-1)\
	}

#define GKPY_OMETHOD_BEGIN(name, method)\
	PyObject* Py##name##_##method(PyObject *selfo, PyObject *args, PyObject *kwds)\
	{\
		GKPY_TRY\
		[[maybe_unused]] Py##name* self = *((Py##name**)&selfo);

#define GKPY_OMETHOD_END\
	GKPY_CATCH_RETURN_NULL\
	}

#define GKPY_IMETHOD_BEGIN(name, method)\
	int Py##name##_##method(PyObject *selfo, PyObject *args, PyObject *kwds)\
	{\
		GKPY_TRY\
		[[maybe_unused]] Py##name* self = *((Py##name**)&selfo);

#define GKPY_IMETHOD_END\
		GKPY_CATCH_RETURN_VALUE(-1)\
	}

#define GKPY_RICHCOMPARE_BEGIN(name)\
	PyObject* Py##name##_RichCompare(PyObject* ao, PyObject* bo, int op) \
	{ \
		GKPY_TRY \
		if (Py_TYPE(ao) != Py_TYPE(bo)) { \
			if (op == Py_EQ) GKPY_RETURN_FALSE; \
			if (op == Py_NE) GKPY_RETURN_TRUE; \
			GK_THROW(type, "Incompatible arguments '{}' and '{}'.", Py_TYPE(ao)->tp_name, Py_TYPE(bo)->tp_name); \
		} \
		const Py##name::value_t& a = Py##name::value(ao); \
		const Py##name::value_t& b = Py##name::value(bo);

#define GKPY_RICHCOMPARE_END \
		GK_UNREACHABLE(); \
		GKPY_CATCH_RETURN_NULL \
	}

#define GKPY_INIT_BEGIN(name) \
	GKPY_IMETHOD_BEGIN(name, Init) \
		GKPY_TRACE_CTOR(self, name); \

#define GKPY_INIT_END \
		return 0;\
	GKPY_IMETHOD_END

// Used when __new__ is already defined, but we need to define __init__ too in order to avoid
// warning about hitting object.__init__. In that scenario, we don't want the debug constructor/destructor
// trace to print for both __new__ and __init__ since it's confusing.
#define GKPY_INIT_EMPTY(name) \
	GKPY_IMETHOD_BEGIN(name, Init) \
		return 0; \
	GKPY_IMETHOD_END

#define GKPY_MEMBERS_BEGIN(type_name) \
	static PyMemberDef Py##type_name##_Members[] = {

#define GKPY_MEMBER_OBJECT(type_name, member_name, attr_name, docstr) \
    {attr_name, T_OBJECT_EX, offsetof(Py##type_name, member_name), 1, docstr},

#define GKPY_MEMBERS_END \
		{NULL}  /* Sentinel */ \
	};

#define GKPY_GETATTR_MEMBER_ON_DEMAND(name, basetype, args_fmt, ...) \
	do { \
		if (self->name == 0 && !strcmp(attr, #name)) \
			GKPY_CONSTRUCT(name, basetype, args_fmt, __VA_ARGS__); \
	} while (0)

#define GKPY_METHODS_BEGIN(type_name) \
	static PyMethodDef Py##type_name##_Methods[] = {

#define GKPY_METHOD_ENTRY(type_name, method_name, method_type, docstr) \
    {#method_name, (PyCFunction)Py##type_name##_##method_name, method_type, docstr},

#define GKPY_METHODS_END \
		{NULL}  /* Sentinel */ \
	};

#define GKPY_METHOD_BEGIN_VARARGS(name, method_name) \
	PyObject* Py##name##_##method_name(PyObject* selfo, PyObject* args) \
	{ \
		GKPY_TRY \
		[[maybe_unused]] Py##name* self = (Py##name*)selfo;

#define GKPY_METHOD_BEGIN_ONEARG(name, method_name) \
	PyObject* Py##name##_##method_name(PyObject* selfo, PyObject* arg) \
	{ \
		GKPY_TRY \
		[[maybe_unused]] Py##name* self = (Py##name*)selfo;

#define GKPY_METHOD_BEGIN_NOARG(name, method_name) \
	PyObject* Py##name##_##method_name(PyObject* selfo) \
	{ \
		GKPY_TRY \
		[[maybe_unused]] Py##name* self = (Py##name*)selfo;

#define GKPY_METHOD_END \
		GKPY_CATCH_RETURN_NULL \
	}

#define GKPY_TYPECHECK(obj_ptr, type_ptr) \
	GK_CHECK(PyObject_TypeCheck(obj_ptr, type_ptr), type, "Expected subtype of '{}', not '{}'", (type_ptr)->tp_name, (obj_ptr)->ob_type->tp_name);

#define GKPY_TYPECHECK_EXACT(obj_ptr, type_ptr) \
	GK_CHECK((obj_ptr)->ob_type == (type_ptr), type, "Expected type '{}', not '{}'", (type_ptr)->tp_name, (obj_ptr)->ob_type->tp_name);

#define GKPY_INDEXCHECK(index, size) \
	GK_CHECK(index >= 0 && index < size, index, "index out of range")

#define GKPY_TYPECHECK_BUILTIN(obj_ptr, pytype) \
	GK_CHECK(pytype##_Check(obj_ptr), type, "Expected type '{}', not '{}'", pytype##_Type.tp_name, (obj_ptr)->ob_type->tp_name);

////////////////////////////////////////////////////////////////////////

#define GKPY_DECLARE_STRINGTABLE(name, ctype, cname, num) \
	extern PyObject* g_##cname##_as_pystring[num]; \
	INLINE PyObject* PyString_From##name(ctype cname) \
	{ \
		PyObject* s;                                                                     \
		if ((int)cname < (int)num)                                                       \
			s = g_##cname##_as_pystring[(int)cname]; /* If valid index, return the string */  \
		else if ((int)cname == (int)num)                                                 \
			s = Py_None;   /* If index = num, return None, by convention */              \
		else                                                                             \
			GK_THROW(index, "Invalid index {} in PyString_From" #name, cname);      \
		Py_INCREF(s);  /* Increment reference count so that we're returning a new ref */ \
		return s;      /* which can be directly returned as a result to Python */        \
	}

/////////////////////////////////////////////////////////////////////////////////////////////

#define GKPY_CONSTRUCT(name, basetype, args_fmt, ...) \
	do { \
		GKPY_BUILDVALUE_TEMP(args, args_fmt, __VA_ARGS__); \
		self->name = Py##basetype::DefaultType->tp_new(Py##basetype::DefaultType, args, NULL); \
		if (!self->name || -1 == Py##basetype::DefaultType->tp_init(self->name, args, NULL)) \
			GKPY_RETURN_ERROR; \
	} while (0)

// Tracks an object for destruction if not explicitly released.
// Useful during constructors/initialization routines where if an error occurs part-way need to
// destruct all the members that have been created so far (i.e. partial object).
#define GKPY_TENTATIVE_CONSTRUCT(name, basetype, args_fmt, ...) \
	GKPY_CONSTRUCT(name, basetype, args_fmt, __VA_ARGS__); \
	PyAutoRef name##_tag(self->name)

#define GKPY_FINALIZE_CONSTRUCT(name) \
	name##_tag.release()

////////////////////////////////////////////////////////////////////////

// Can only be used for non-base classes
template <typename T>
PyObject* PyFastNew(PyTypeObject* type, PyObject* args, PyObject* kwds)
{
	if (PyType_IS_GC(type)) {
		// Fall back to generic allocation pathway
		PyObject* obj = type->tp_alloc(type, 1);
		default_construct((T*)obj);    // default construct inplace
		return obj;
	} else {
		// Calling PyObject_Malloc directly is about 8% faster than through type->tp_alloc
		// Only works for non-garbage-collected types, obviously.
		auto* obj = (PyObject*)PyObject_Malloc(type->tp_basicsize);
		default_construct((T*)obj);    // default construct inplace
		if ((size_t)type->tp_basicsize == sizeof(T)) {
			memset(obj, 0, sizeof(T)); // allow const propagation into memset intrinsic for case sizeof(T)
		} else {
			memset(obj, 0, (size_t)type->tp_basicsize);
		}
		if (PyType_HasFeature(type, Py_TPFLAGS_HEAPTYPE))
			Py_INCREF(type); // DECREF in Python's subtype_dealloc
		return PyObject_INIT(obj, type);
	}
}

// Can only be used for objects allocated directly by PyObject_Malloc, such as for PyFastNew
template <typename T>
void PyFastDealloc(PyObject* self)
{
	destruct((T*)self);  // destruct inplace
	if (PyType_IS_GC(Py_TYPE(self)))
		Py_TYPE(self)->tp_free(self);
	else
		PyObject_Free(self);
}

/*
PyForceNewGCInheritance(type)

Disables the Py_TPFLAGS_HAVE_GC on a subclass, matching the new and
improved inheritance behaviour of Python 3.5.

MUST BE CALLED BEFORE ANY INSTANCES OF THE CLASS ARE CREATED.

Useful for types that have no internal references to other Python
objects and yet have many of instances. Disabling the flag
makes object allocation/deallocation faster and avoids having the
garbage collector track tens of millions of objects that can't
actually create cycles.

This function is not necessary in Python 3.5, where the logic for
Py_TPFLAGS_HAVE_GC inheritance was updated to avoid this problem.

    typeobject.c:2642 (Python 3.5)
    --------------------------------------------------
    // Enable GC unless this class is not adding new instance variables
    // and the base class did not use GC.
    if ((base->tp_flags & Py_TPFLAGS_HAVE_GC) ||
        type->tp_basicsize > base->tp_basicsize)
        type->tp_flags |= Py_TPFLAGS_HAVE_GC;
    --------------------------------------------------

    typeobject.c:2496 (Python 2.7)
    --------------------------------------------------
    // Enable GC unless there are really no instance variables possible
    if (!(type->tp_basicsize == sizeof(PyObject) &&
        type->tp_itemsize == 0))
        type->tp_flags |= Py_TPFLAGS_HAVE_GC;
    --------------------------------------------------
*/
void PyForceNewGCInheritance(PyTypeObject* type);
void PyDisableGC(PyTypeObject* type);

// Check that a type and its base type have matching size.
// Useful when having a subtype with larger size would either be a performance
// problem (for lightweight objects) or would cause memory corruption
// (in case of inplace objects).
void PyCheckSameBasicSize(PyTypeObject* type);

// Deletes any method or member defined on 'type' that is also and
// type->tp_base, under the assumption that the version attached
// to 'type' is merely a mock or wrapper of the base type's method.
// Useful for getting rid of overhead of wrappers without having
// to explicitly handle callables through __getattribute__.
void PyDeleteMockAttrs(PyTypeObject* type);

//////////////////////////////////////////////////////////////

// Returns a new string, or new reference to None if the argument is an empty string
PyObject* PyString_FromNonEmptyString(const char* str);

template <typename T>
bool value_equals(PyObject* ao, PyObject* bo)
{
	if constexpr (requires(T x) {
					  {
						  x.as_ptr
						  } -> std::convertible_to<void*>;
				  }) { // T = mmappable type like PyGene
		auto a = ((T*)ao)->as_ptr;
		auto b = ((T*)bo)->as_ptr;
		if (a && b) // both are from a mmapped table, which is equal if addresses are equal
			return a == b;
		// fallback to value comparisons if any not mmapped
	}
	// T = pure value type like PyGenome
	return T::value(ao) == T::value(bo);
}

template <typename T>
static PyObject* PyGenericValue_RichCompare(PyObject* ao, PyObject* bo, int op)
{
	GKPY_TRY
	if (op != Py_EQ && op != Py_NE)
		GKPY_RETURN_INCREF(Py_NotImplemented);

	if (Py_TYPE(ao) != Py_TYPE(bo))
		GKPY_RETURN_BOOL(op == Py_NE);

	const auto equals = value_equals<T>(ao, bo);
	GKPY_RETURN_BOOL(op == Py_EQ ? equals : !equals);
	GKPY_CATCH_RETURN_NULL
}

template <typename T> // T = PyGene etc
Py_hash_t PyGenericValue_Hash(PyObject *selfo) { return (Py_hash_t)hash(T::value(selfo)); }

//////////////////////////////////////////////////////////////

// Methods that need to be re-defined
#define _GKPY_VALUE_TYPE_COMMON_MEMBERS(name) \
		/* Get access to value, either via pointer or stored internally */ \
		INLINE value_t& value() { \
			if (as_ptr) { \
				source()->validate(); \
				return *(value_t*)as_ptr; \
			} \
			return _get<value_t>(0); \
		} \
		\
		INLINE static value_t&  value(PyObject* obj) { return ((Py##name*)obj)->value(); } \
		INLINE static bool      check(PyObject* obj) { return PyObject_TypeCheck(obj, DefaultType); }

#define GKPY_VALUE_TYPE_BEGIN(name, value_type) \
	GKPY_TYPE_BEGIN(name) \
		using value_t = value_type; \
		void* as_ptr; /* bytes [16:24] on x64 */ \
		\
		/* constructor/destructor */ \
		       Py##name() = default; \
		INLINE ~Py##name() { if (as_ptr) { Py_DECREF(source()); } } \
		\
		/* get field stored at an arbitrary offset relative to the first byte after the as_ptr member */ \
		template <typename T> INLINE T& _get(size_t offset) const { return *(T*)((char*)&as_ptr + sizeof(as_ptr) + offset); } \
		\
		/* Get access to the source object that owns the value pointed to, if as_ptr is not NULL */ \
		INLINE PyAsPtrSource*& source() const { GK_DBASSERT(as_ptr); return _get<PyAsPtrSource*>(0); } \
		\
		_GKPY_VALUE_TYPE_COMMON_MEMBERS(name)

#define GKPY_VALUE_TYPE_END \
	GKPY_TYPE_END

#define GKPY_VALUE_SUBTYPE_BEGIN(name, basename, value_type) \
	GKPY_SUBTYPE_BEGIN(name, basename) \
		using value_t = value_type; \
		_GKPY_VALUE_TYPE_COMMON_MEMBERS(name)

#define GKPY_VALUE_SUBTYPE_END \
	GKPY_VALUE_TYPE_END

//////////////////////////////////////////////////////////////

// Needed to override default export visibility on GCC when -fvisibility=hidden is set.
#ifdef _MSC_VER
#define EX_PYTHON_INIT
#else
#define EX_PYTHON_INIT __attribute__ ((visibility ("default")))
#endif

int gkpy_invalid_init(PyObject* self, PyObject* args, PyObject* kwargs);

END_NAMESPACE_GK

#endif // __GK_PY_UTIL_H__
