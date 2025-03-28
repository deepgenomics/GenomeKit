/*
Copyright (C) 2016-2023 Deep Genomics Inc. All Rights Reserved.
*/
#pragma once
#ifndef __GENOME_KIT_PY_TABLE_H__
#define __GENOME_KIT_PY_TABLE_H__

#include "table.h"
#include "py_util.h"
#include "py_interval.h"

#include "numpy/arrayobject.h"

BEGIN_NAMESPACE_GK

// T = PyCoord, PyInterval, or similar
template <typename T>
GKPY_SOURCE_TYPE_BEGIN(Table)
	using table_type = typename T::table_type;
	table_type* table;    // Instance of gene_table, exon_table, etc.
	INLINE static table_type& value(PyObject* self) { return *((PyTable<T>*)self)->table; }
	static PyMethodDef Methods[];
	GKPY_SOURCE_TYPE_END

	template <typename T>
PyObject* PyTable_CreateItem(PyObject* self, const typename T::value_t* as_ptr)
{
	if (!as_ptr)
		GKPY_RETURN_NONE;
	// Create a new instance of T::DefaultType and initialize it
	// to point to the element in the table.
	T* elem = (T*)T::DefaultType->tp_new(T::DefaultType, nullptr, nullptr);
	elem->as_ptr = (void*)as_ptr;
	elem->source() = (PyAsPtrSource*)self; Py_INCREF(self);
	return (PyObject*)elem;
}

template <typename T>
Py_ssize_t PyTable_Len(PyObject* self) { return (Py_ssize_t)PyTable<T>::value(self).size(); }

template <typename T>
PyObject* PyTable_GetItem(PyObject* self, Py_ssize_t index)
{
	GKPY_TRY
	typename PyTable<T>::table_type& table = PyTable<T>::value(self);
	GKPY_INDEXCHECK(index, table.size());
	return PyTable_CreateItem<T>(self, &table[(index_t)index]);
	GKPY_CATCH_RETURN_NULL
}

template <typename T>
PyObject* PyTable_index_of(PyObject* self, PyObject* arg)
{
	GKPY_TRY
	GKPY_TYPECHECK(arg, T::DefaultType);
	typename PyTable<T>::table_type& table = PyTable<T>::value(self);
	typename T::value_t& value = T::value(arg);
	return PyInt_FromLong(table.index_of(value));
	GKPY_CATCH_RETURN_NULL
}

template <typename T>
PyObject* PyTable_where(PyObject* self, PyObject* arg)
{
	GKPY_TRY
	GKPY_TYPECHECK(arg, &PyArray_Type);
	auto* mask = (PyArrayObject*)arg;
	typename PyTable<T>::table_type& table = PyTable<T>::value(self);
	GK_CHECK2(PyArray_NDIM(mask) == 1 || (PyArray_NDIM(mask) == 2 && PyArray_DIMS(mask)[1] == 1), value, "Expected 1-dimensional mask");
	GK_CHECK2(PyArray_DIMS(mask)[0] == (npy_intp)table.size(), value, "Expected mask of length {}", table.size());
	GK_CHECK2(PyArray_ISBOOL(mask), value, "Expected boolean mask");
	const bool* mask_data = (const bool*)PyArray_DATA(mask);

	// First pass to collect nonzeros
	vector<int> nz;
	for (int i = 0; i < table.size(); ++i)
		if (mask_data[i])
			nz.push_back(i);

	// Create list with a slot for each variant.
	int nnz = (int)nz.size();
	PyObject* list = PyList_New(nnz);
	GKPY_TAKEREF(list);
	for (int i = 0; i < nnz; ++i)
		PyList_SET_ITEM(list, i, PyTable_CreateItem<T>(self, &table[nz[i]]));
	GKPY_FORGETREF(list); // keep list refcount at 1
	return list;

	GKPY_CATCH_RETURN_NULL
}

template <typename T>
PyMethodDef PyTable<T>::Methods[] = {
	{"index_of",  (PyCFunction)PyTable_index_of<T>, METH_O, nullptr},
	{"where",     (PyCFunction)PyTable_where<T>,    METH_O, nullptr},
	{nullptr} // sentinel
};

GKPY_TEMPLATE_TYPEOBJ_BEGIN(template <typename T>, Table<T>)
	tp_flags |= Py_TPFLAGS_BASETYPE;
	tp_iter = PySeqIter_New;    // Standard sequence iterator calls sq_item until IndexError is raised
	sq_item = PyTable_GetItem<T>;
	sq_length = PyTable_Len<T>;
	tp_methods = PyTable<T>::Methods;
GKPY_TEMPLATE_TYPEOBJ_END

///////////////////////////////////////////////////////////////

// For range types R that have a known size(), create a pre-sized tuple
template <typename T, typename R>
PyObject* PyTuple_FromSizedElemRange(PyObject* self, const R& range)
{
	auto size = (Py_ssize_t)range.size();
	PyObject* tuple = PyTuple_New(size);
	GKPY_TAKEREF(tuple);
	typename R::cursor_t cur = range.begin();
	for (Py_ssize_t i = 0; i < size; ++i, ++cur) {
		PyTuple_SET_ITEM(tuple, i, PyTable_CreateItem<T>(self, &*cur));
	}
	GKPY_FORGETREF(tuple); // keep list refcount at 1
	return tuple;
}

// For range types R that have a known size(), create a pre-sized list
template <typename T, typename R>
PyObject* PyList_FromSizedElemRange(PyObject* self, const R& range)
{
	auto size = (Py_ssize_t)range.size();
	PyObject* list = PyList_New(size);
	GKPY_TAKEREF(list);
	typename R::cursor_t cur = range.begin();
	for (Py_ssize_t i = 0; i < size; ++i, ++cur) {
		PyList_SET_ITEM(list, i, PyTable_CreateItem<T>(self, &*cur));
	}
	GKPY_FORGETREF(list); // keep list refcount at 1
	return list;
}

// For range types R that don't a size() known ahead of time, build the list on the fly
template <typename T, typename R>
PyObject* PyList_FromElemRange(PyObject* self, const R& range)
{
	PyObject* list = PyList_New(0);
	GKPY_TAKEREF(list);
	for (const auto& cursor: range) {
		PyObject* item = PyTable_CreateItem<T>(self, &cursor);
		GKPY_TAKEREF(item);
		if (PyList_Append(list, item))  // INCREF(value)
			return nullptr;
	}                                   // DECREF(value)
	GKPY_FORGETREF(list); // keep list refcount at 1
	return list;
}

///////////////////////////////////////////////////////////////

// T = PyInterval, or similar
template <typename T>
GKPY_SUBTYPE_BEGIN(IntervalTable, Table<T>)
	static PyMethodDef Methods[];
GKPY_SUBTYPE_END

template <typename T>
PyObject* PyIntervalTable_find_5p_aligned(PyObject* selfo, PyObject* arg)
{
	GKPY_TRY
	auto* self = (PyIntervalTable<T>*)selfo;
	typename T::table_type::cursor_range r = self->table->find_5p_aligned(PyAsInterval(arg));
	return PyList_FromSizedElemRange<T>(selfo, r);
	GKPY_CATCH_RETURN_NULL
}

template <typename T>
PyObject* PyIntervalTable_find_3p_aligned(PyObject* selfo, PyObject* arg)
{
	GKPY_TRY
	auto* self = (PyIntervalTable<T>*)selfo;
	typename T::table_type::cursor_range r = self->table->find_3p_aligned(PyAsInterval(arg));
	return PyList_FromSizedElemRange<T>(selfo, r);
	GKPY_CATCH_RETURN_NULL
}

template <typename T>
PyObject* PyIntervalTable_find_5p_within(PyObject* selfo, PyObject* arg)
{
	GKPY_TRY
	auto* self = (PyIntervalTable<T>*)selfo;
	typename T::table_type::cursor_range r = self->table->find_5p_within(PyAsInterval(arg));
	return PyList_FromSizedElemRange<T>(selfo, r);
	GKPY_CATCH_RETURN_NULL
}

template <typename T>
PyObject* PyIntervalTable_find_3p_within(PyObject* selfo, PyObject* arg)
{
	GKPY_TRY
	auto* self = (PyIntervalTable<T>*)selfo;
	typename T::table_type::cursor_range r = self->table->find_3p_within(PyAsInterval(arg));
	return PyList_FromSizedElemRange<T>(selfo, r);
	GKPY_CATCH_RETURN_NULL
}

template <typename T>
PyObject* PyIntervalTable_find_within(PyObject* selfo, PyObject* arg)
{
	GKPY_TRY
	auto* self = (PyIntervalTable<T>*)selfo;
	typename T::table_type::cursor_range_within r = self->table->find_within(PyAsInterval(arg));
	return PyList_FromElemRange<T>(selfo, r);
	GKPY_CATCH_RETURN_NULL
}

template <typename T>
PyObject* PyIntervalTable_find_overlapping(PyObject* selfo, PyObject* arg)
{
	GKPY_TRY
	auto* self = (PyIntervalTable<T>*)selfo;
	typename T::table_type::cursor_range_overlapping r = self->table->find_overlapping(PyAsInterval(arg));
	return PyList_FromElemRange<T>(selfo, r);
	GKPY_CATCH_RETURN_NULL
}

template <typename T>
PyObject* PyIntervalTable_find_exact(PyObject* selfo, PyObject* arg)
{
	GKPY_TRY
	auto* self = (PyIntervalTable<T>*)selfo;
	typename T::table_type::cursor_range r = self->table->find_exact(PyAsInterval(arg));
	return PyList_FromElemRange<T>(selfo, r);
	GKPY_CATCH_RETURN_NULL
}

template <typename T>
PyObject* PyIntervalTable_getattro(PyObject* o, PyObject* attr_name)
{
	auto self = scast<PyIntervalTable<T>*>(o);
	auto attr = PyString_AS_STRING(attr_name);
	GKPY_GETATTR_CASE("stranded") GKPY_RETURN_BOOL(self->table->stranded());
	return PyObject_GenericGetAttr(o, attr_name);
}

template <typename T>
int PyIntervalTable_setattro(PyObject* o, PyObject* attr_name, PyObject* v)
{
	auto selfo = o;
	auto attr  = PyString_AS_STRING(attr_name);
	GKPY_SETATTR_READONLY("stranded");
	return PyObject_GenericSetAttr(o, attr_name, v);
}

template <typename T>
PyMethodDef PyIntervalTable<T>::Methods[] = {
	{"find_5p_aligned",  (PyCFunction)PyIntervalTable_find_5p_aligned<T>,   METH_O, nullptr},
	{"find_3p_aligned",  (PyCFunction)PyIntervalTable_find_3p_aligned<T>,   METH_O, nullptr},
	{"find_5p_within",   (PyCFunction)PyIntervalTable_find_5p_within<T>,    METH_O, nullptr},
	{"find_3p_within",   (PyCFunction)PyIntervalTable_find_3p_within<T>,    METH_O, nullptr},
	{"find_within",      (PyCFunction)PyIntervalTable_find_within<T>,       METH_O, nullptr},
	{"find_overlapping", (PyCFunction)PyIntervalTable_find_overlapping<T>,  METH_O, nullptr},
	{"find_exact",       (PyCFunction)PyIntervalTable_find_exact<T>,        METH_O, nullptr},
	{nullptr} // sentinel
};

GKPY_TEMPLATE_SUBTYPEOBJ_BEGIN(template <typename T>, IntervalTable<T>)
	tp_flags |= Py_TPFLAGS_BASETYPE;
	tp_getattro = (getattrofunc)PyIntervalTable_getattro<T>;
	tp_setattro = (setattrofunc)PyIntervalTable_setattro<T>;
	tp_methods  = PyIntervalTable<T>::Methods;
GKPY_TEMPLATE_SUBTYPEOBJ_END

END_NAMESPACE_GK

#endif // __GENOME_KIT_PY_TABLE_H__
