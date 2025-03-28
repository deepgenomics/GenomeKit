/*
Copyright (C) 2016-2023 Deep Genomics Inc. All Rights Reserved.
*/
#include "py_util.h"
#include "py_interval.h"

#include "format.h"
#include "genome.h"
#include "py_genome.h"

#define GK_CHECK_REFG(a, b) \
	GK_CHECK((a).refg == (b).refg, value, "Coordinate system mismatch, {} and {}.", (a), (b));

#define GKPY_INTERVAL_BOOL_METHOD_ONEARG(pyname, method) \
	GKPY_METHOD_BEGIN_ONEARG(pyname, method) \
		const Py##pyname::value_t& c = Py##pyname::value(selfo); \
		if (PyInterval::check(arg)) { const interval_t& a = PyInterval::value(arg); GK_CHECK_REFG(c, a); GKPY_RETURN_BOOL(c.method(a)); } \
		GK_THROW(type, "argument must be an Interval, not '{}'", Py_TYPE(arg)->tp_name); \
	GKPY_METHOD_END

#define GKPY_INTERVAL_BOOL_METHOD_NOARG(pyname, method) \
	GKPY_METHOD_BEGIN_ONEARG(pyname, method) \
		GKPY_RETURN_BOOL(Py##pyname::value(selfo).method()); \
	GKPY_METHOD_END

BEGIN_NAMESPACE_GK

PyObject* g_strand_as_pystring[num_strand];

void Init_Interval_PyStrings()
{
	for (std::underlying_type_t<strand_t> i = 0; i < num_strand; ++i) {
		char s[2] = { strand_as_char(strand_t{i}), '\0' };
		g_strand_as_pystring[i] = PyString_FromString(s); // new ref
	}
}

////////////////////////////////////////////////////////////////

refg_t as_refg(PyObject* arg)
{
	// TODO: data_dir injected as a context
	if (PyString_Check(arg))
		return get_refg_registry().as_refg(PyString_AS_STRING(arg));
	if (PyObject_TypeCheck(arg, PyGenome::DefaultType))
		return ((PyGenome*)arg)->genome.refg();
	PyObject* refg_attr = PyObject_GetAttrString(arg, "reference_genome");
	GKPY_TAKEREF(refg_attr);
	GK_CHECK2(refg_attr, value, "Expected reference_genome to either be a string or an object with a reference_genome attribute");
	if (PyString_Check(refg_attr))
		return get_refg_registry().as_refg(PyString_AS_STRING(refg_attr));
	GK_THROW2(value, "Could not determine reference genome from reference_genome argument");
}

const genome_t& as_genome(PyObject* arg)
{
	GK_CHECK(PyType_IsSubtype(Py_TYPE(arg), PyGenome::DefaultType) == 1, type, "Expected Genome, got '{}'.",
			 Py_TYPE(arg)->tp_name);
	return ((PyGenome*)arg)->genome;
}

pos_t as_pos(PyObject* arg, const interval_t& interval, const char* name)
{
	if (arg == Py_None)
		return invalid_pos;
	if (PyInt_Check(arg))
		return (pos_t)PyInt_AS_LONG(arg);
	if (PyInterval::check(arg)) {
		const interval_t& c = PyInterval::value(arg);
		GK_CHECK(c.same_strand(interval), value, "Expected {} to be in same strand / coordinate system", name);
		GK_CHECK(c.empty(), value, "Expected {} to be empty interval.", name);
		return c.start();
	}
	GK_THROW(type, "Expected {} to be of type int", name);
}

PyObject* PyInterval___getstate__(PyObject* self)
{
	GKPY_TRY
		ainterval_t v;
		*static_cast<interval_t*>(&v) = PyInterval::value(self);
		v.anchor        = ((PyInterval*)self)->get_anchor();
		v.anchor_offset = ((PyInterval*)self)->get_anchor_offset();
		return PyBytes_FromStringAndSize(reinterpret_cast<const char*>(&v), (Py_ssize_t)sizeof(v));
	GKPY_CATCH_RETURN_NULL
}

PyObject* PyInterval___setstate__(PyObject* self, PyObject* state)
{
	GKPY_TRY
		GK_CHECK2(PyBytes_Check(state), type, "Expected string type");
		GK_CHECK2(PyBytes_GET_SIZE(state) == (Py_ssize_t)sizeof(ainterval_t), value, "Expected %d bytes", (int)sizeof(ainterval_t));
		((PyInterval*)self)->as_ptr = nullptr;
		memcpy(&PyInterval::value(self), PyBytes_AsString(state), sizeof(ainterval_t));
		GKPY_RETURN_NONE;
	GKPY_CATCH_RETURN_NULL
}

GKPY_IMETHOD_BEGIN(Interval, Init)
	if (PyType_IS_GC(Py_TYPE(self))) {
		// Print to stderr so that output gets past Python unittest redirect, to ensure user sees the error msg
		print("DANGER: Must use @genome_kit.register on Interval subclass {}, expecting to crash now...", Py_TYPE(self)->tp_name);
		GK_THROW(runtime, "DANGER: Must use @genome_kit.register on Interval subclass {}, expecting to crash now...", Py_TYPE(self)->tp_name);
	}
	char* arg_names[]
		= { "chromosome", "strand", "start", "end", "reference_genome", "anchor", "anchor_offset", nullptr };
	int start, end;
	const char* strand;
	const char* chrom;
	PyObject* pyrefg;
	PyObject* anchor = Py_None;
	int anchor_offset = 0;
	if (!PyArg_ParseTupleAndKeywords(args, kwds, "ssiiO|Oi", arg_names, &chrom, &strand, &start, &end, &pyrefg, &anchor,
									 &anchor_offset))
		return -1;
	if (start > end) {
		PyErr_Format(PyExc_ValueError, "Requires start <= end but %d > %d", start, end);
		return -1;
	}
	GK_CHECK(strlen(strand) == 1, value, "Expected strand to be '+' or '-' but found '{}'.", strand);
	ainterval_t& val           = self->value();
	const auto   refg          = as_refg(pyrefg);
	PyObject*    data_dir_attr = PyObject_GetAttrString(pyrefg, "data_dir");
	PyErr_Clear();
	GKPY_TAKEREF(data_dir_attr);
	const auto& chrom_names =
		data_dir_attr ? get_chrom_names(refg, PyString_AsString(data_dir_attr)) : get_chrom_names(refg);

	val        = {interval_t::from_dna0(chrom_names.as_chrom(chrom), start, end, as_strand(strand[0]), refg), 0,
				  anchor_offset};
	val.anchor = as_pos(anchor, val, "anchor");
	return 0;
GKPY_IMETHOD_END

///////////////////////////////

// Special comparison function knows that if a PyInterval points to
// its value (as_ptr != NULL) then the actual object it's pointing to is from
// one of the C++ tables and is therefore an interval_t, not an
// ainterval_t. If either ao or bo (or both) have their interval
// stored on the Py object itself, then that means there is at least one
// actual anchor/anchor_offset to check before determining if a and b are
// equal. If both ap and bo point to their value, then their anchor/anchor_offset
// are all implicitly invalid_pos and there's no need to compare anything further.
int cmp_implicit_anchor(PyInterval* a, PyInterval* b)
{
	if (a->value() < b->value()) // First check global operator for interval_t, which ignores anchor
		return -1;
	if (a->value() > b->value())
		return 1;
	if (a->get_anchor() < b->get_anchor())
		return -1;
	if (a->get_anchor() > b->get_anchor())
		return 1;
	if (a->get_anchor_offset() < b->get_anchor_offset())
		return -1;
	if (a->get_anchor_offset() > b->get_anchor_offset())
		return 1;
	return 0;
}

GKPY_RICHCOMPARE_BEGIN(Interval)
	GK_CHECK(a.refg == b.refg,    value, "Coordinate system mismatch, {} and {}.", a, b);
	//GK_CHECK(a.sys    == b.sys,    value, "Cannot compare {} to {}: different coordinate systems.", a, b);
	//GK_CHECK(a.strand == b.strand, value, "Cannot compare {} to {}: different strands.", a, b);
	switch (op) {
	case Py_EQ: GKPY_RETURN_BOOL(cmp_implicit_anchor((PyInterval*)ao, (PyInterval*)bo) == 0);
	case Py_NE: GKPY_RETURN_BOOL(cmp_implicit_anchor((PyInterval*)ao, (PyInterval*)bo) != 0);
	case Py_LT: GKPY_RETURN_BOOL(cmp_implicit_anchor((PyInterval*)ao, (PyInterval*)bo) < 0);
	case Py_GT: GKPY_RETURN_BOOL(cmp_implicit_anchor((PyInterval*)ao, (PyInterval*)bo) > 0);
	case Py_LE: GKPY_RETURN_BOOL(cmp_implicit_anchor((PyInterval*)ao, (PyInterval*)bo) <= 0);
	case Py_GE: GKPY_RETURN_BOOL(cmp_implicit_anchor((PyInterval*)ao, (PyInterval*)bo) >= 0);
	}
GKPY_RICHCOMPARE_END

///////////////////////////////

PyObject*  PyInterval_Str(PyObject*  selfo)   { return PyString_FromSV(PyInterval::value(selfo).as_str()); }
Py_ssize_t PyInterval_Length(PyObject *selfo) { return (Py_ssize_t)PyInterval::value(selfo).size(); }

///////////////////////////////

PyObject* PyString_FromChrom(chrom_key_t<refg_t> key)
{
	// The string values are never decref'd in order to avoid a segfault on exit
	// see https://github.com/python/cpython/issues/126508
	static chrom_map_t<refg_t, PyObject*> map;

	auto& s = map[key];
	if (s == nullptr) {
		// TODO: data_dir injected as a context
		s = PyString_FromSV(get_chrom_names(key.other).chrom_as_sv(key.chrom));
	}
	Py_INCREF(s); /* Increment reference count so that we're returning a new ref */
	return s;     /* which can be directly returned as a result to Python */
}

GKPY_GETATTRO_BEGIN(Interval)
	const interval_t& c = self->value();
	switch (attr[0]) {
	case 'c':  GKPY_GETATTR_CASE("chrom")        return PyString_FromChrom({c.chrom, c.refg});
	           GKPY_GETATTR_CASE("chromosome")   return PyString_FromChrom({c.chrom, c.refg});
	           break;
	case 's':  GKPY_GETATTR_CASE("strand")       return PyString_FromStrand(c.strand);
	           GKPY_GETATTR_CASE("start")        return PyInt_FromLong((long)c.start());
			   // TODO: GEN-37 data_dir injected as a context
			   GKPY_GETATTR_CASE("sys")          return PyString_FromSV(format("{}:{}", get_refg_registry().refg_as_sv(c.refg), get_chrom_names(c.refg).chrom_as_sv(c.chrom)));
	           break;
	case 'e':  GKPY_GETATTR_CASE("end")          return PyInt_FromLong((long)c.end());
	           GKPY_GETATTR_CASE("end5")         return PyInterval::create(get_end5(c), self->get_anchor(), self->get_anchor_offset());
	           GKPY_GETATTR_CASE("end3")         return PyInterval::create(get_end3(c), self->get_anchor(), self->get_anchor_offset());
	           break;
	case 'a':  GKPY_GETATTR_CASE("anchor")        { pos_t a = self->get_anchor(); if (a == invalid_pos) GKPY_RETURN_NONE; return PyInt_FromLong((long)a); }
	           GKPY_GETATTR_CASE("anchor_offset") { return PyInt_FromLong((long)self->get_anchor_offset()); }
	           break;
			   // TODO: GEN-37 data_dir injected as a context
	case 'r':  GKPY_GETATTR_CASE("refg")             return PyString_FromSV(get_refg_registry().refg_as_sv(c.refg));
			   GKPY_GETATTR_CASE("reference_genome") return PyString_FromSV(get_refg_registry().refg_as_sv(c.refg));
			   break;
	case 'i':  GKPY_GETATTR_CASE("interval")         GKPY_RETURN_INCREF(selfo);
	           break;
	}
	return PyObject_GenericGetAttr(selfo, attro);
GKPY_GETATTRO_END

///////////////////////////////

GKPY_SETATTRO_BEGIN(Interval)
	// ******** UPDATE genome_kit.Interval MOCK PROPERTIES TO MATCH ANY CHANGES HERE ********
	GKPY_SETATTR_READONLY("chrom");
	GKPY_SETATTR_READONLY("chromosome");
	GKPY_SETATTR_READONLY("strand");
	GKPY_SETATTR_READONLY("start");
	GKPY_SETATTR_READONLY("sys");
	GKPY_SETATTR_READONLY("end");
	GKPY_SETATTR_READONLY("end5");
	GKPY_SETATTR_READONLY("end3");
	GKPY_SETATTR_READONLY("interval");
	GKPY_SETATTR_READONLY("anchor");
	GKPY_SETATTR_READONLY("anchor_offset");
	GKPY_SETATTR_READONLY("refg");
	GKPY_SETATTR_READONLY("reference_genome");
	return PyObject_GenericSetAttr(selfo, attro, value);
GKPY_SETATTRO_END

///////////////////////////////

GKPY_METHOD_BEGIN_ONEARG(Interval, shift)
	GKPY_TYPECHECK_BUILTIN(arg, PyInt);
	return PyInterval::create(self->value().shift((pos_t)PyInt_AS_LONG(arg)), self->get_anchor(), self->get_anchor_offset());
GKPY_METHOD_END
GKPY_METHOD_BEGIN_VARARGS(Interval, expand)
	pos_t up;
	pos_t dn;
	Py_ssize_t size = PyTuple_GET_SIZE(args);
	if (size == 1) {
		PyObject* arg = PyTuple_GET_ITEM(args, 0); GKPY_TYPECHECK_BUILTIN(arg, PyInt);
		up = dn = (pos_t)PyInt_AS_LONG(arg);
	} else if (size == 2) {
		PyObject* arg0 = PyTuple_GET_ITEM(args, 0); GKPY_TYPECHECK_BUILTIN(arg0, PyInt);
		PyObject* arg1 = PyTuple_GET_ITEM(args, 1); GKPY_TYPECHECK_BUILTIN(arg1, PyInt);
		up = (pos_t)PyInt_AS_LONG(arg0);
		dn = (pos_t)PyInt_AS_LONG(arg1);
	} else {
		GK_THROW(value, "Expected 1 or 2 arguments but got {}", size);
	}
	return PyInterval::create(self->value().expand(up, dn),
	                          self->get_anchor(), self->get_anchor_offset());
GKPY_METHOD_END
GKPY_METHOD_BEGIN_ONEARG(Interval, intersect)
	auto x = PyInterval::value(selfo);
	if (PyInterval::check(arg)) {
		auto y = PyInterval::value(arg);
		if (((PyInterval*)selfo)->get_anchor() != invalid_pos || ((PyInterval*)arg)->get_anchor() != invalid_pos) {
			GK_THROW(value, "anchored intersection ({}, {}) is not supported.", x, y);
		}

		if (x.end() <= y.start() || y.end() <= x.start() || x.strand != y.strand || x.refg != y.refg) {
			GKPY_RETURN_NONE;
		}

		auto start = max(x.start(), y.start());
		auto end   = min(x.end(), y.end());

		return PyInterval::create(interval_t::from_dna0(x.chrom, start, end, x.strand, x.refg));
	}
	GK_THROW(type, "argument must be Interval, not '{}'", Py_TYPE(arg)->tp_name);
GKPY_METHOD_END
GKPY_INTERVAL_BOOL_METHOD_ONEARG(Interval, upstream_of)
GKPY_INTERVAL_BOOL_METHOD_ONEARG(Interval, dnstream_of)
GKPY_INTERVAL_BOOL_METHOD_ONEARG(Interval, contains)
GKPY_INTERVAL_BOOL_METHOD_ONEARG(Interval, within)
GKPY_INTERVAL_BOOL_METHOD_ONEARG(Interval, overlaps)
GKPY_INTERVAL_BOOL_METHOD_NOARG(Interval, is_pos_strand)
GKPY_METHOD_BEGIN_ONEARG(Interval, as_positive_strand)
	return PyInterval::create(PyInterval::value(selfo).as_pos_strand(), self->get_anchor(), self->get_anchor_offset());
GKPY_METHOD_END
GKPY_METHOD_BEGIN_ONEARG(Interval, as_negative_strand)
	return PyInterval::create(PyInterval::value(selfo).as_neg_strand(), self->get_anchor(), self->get_anchor_offset());
GKPY_METHOD_END
GKPY_METHOD_BEGIN_ONEARG(Interval, as_opposite_strand)
	return PyInterval::create(PyInterval::value(selfo).as_opp_strand(), self->get_anchor(), self->get_anchor_offset());
GKPY_METHOD_END

// C++ short names to Python long names
#define PyInterval_is_positive_strand PyInterval_is_pos_strand

//////////////////////////////

GKPY_METHODS_BEGIN(Interval)
	// ******** UPDATE genome_kit.Interval MOCK METHODS TO MATCH ANY CHANGES HERE ********
	GKPY_METHOD_ENTRY(Interval, shift,  METH_O, nullptr)
	GKPY_METHOD_ENTRY(Interval, expand, METH_VARARGS, nullptr)
	GKPY_METHOD_ENTRY(Interval, intersect, METH_O, nullptr)
	GKPY_METHOD_ENTRY(Interval, upstream_of, METH_O, nullptr)
	GKPY_METHOD_ENTRY(Interval, dnstream_of, METH_O, nullptr)
	GKPY_METHOD_ENTRY(Interval, contains, METH_O, nullptr)
	GKPY_METHOD_ENTRY(Interval, within,   METH_O, nullptr)
	GKPY_METHOD_ENTRY(Interval, overlaps, METH_O, nullptr)
	GKPY_METHOD_ENTRY(Interval, is_positive_strand, METH_NOARGS, nullptr)
	GKPY_METHOD_ENTRY(Interval, as_positive_strand, METH_NOARGS, nullptr)
	GKPY_METHOD_ENTRY(Interval, as_negative_strand, METH_NOARGS, nullptr)
	GKPY_METHOD_ENTRY(Interval, as_opposite_strand, METH_NOARGS, nullptr)
	GKPY_METHOD_ENTRY(Interval, __getstate__, METH_NOARGS, nullptr)
	GKPY_METHOD_ENTRY(Interval, __setstate__, METH_O, nullptr)
GKPY_METHODS_END

//////////////////////////////

int PyInterval_Contains(PyObject* selfo, PyObject* arg)
{
	GKPY_TRY
	const interval_t& i = PyInterval::value(selfo);
	if (PyInterval::check(arg)) { const interval_t& a = PyInterval::value(arg); GK_CHECK_REFG(i, a); return i.contains(a) ? 1 : 0; }
	GK_THROW(type, "argument must be Interval, not '{}'", Py_TYPE(arg)->tp_name);
	GKPY_CATCH_RETURN_VALUE(-1)
}

///////////////////////////////

GKPY_VALUE_TYPEOBJ_BEGIN(Interval, as_value_or_ptr) // as_value_or_ptr to ensure there's enough tp_basicsize to hold an interval_t value inplace
	tp_flags |= Py_TPFLAGS_BASETYPE;
	tp_init = PyInterval_Init;
	tp_new = PyFastNew<PyInterval>;
	tp_dealloc = PyFastDealloc<PyInterval>;
	tp_getattro = PyInterval_GetAttro;
	tp_setattro = PyInterval_SetAttro;
	tp_richcompare = PyInterval_RichCompare;
	tp_str = PyInterval_Str;
	tp_hash = PyGenericValue_Hash<PyInterval>;
	tp_methods = PyInterval_Methods;
	sq_length = PyInterval_Length;
	sq_contains = PyInterval_Contains;
GKPY_TYPEOBJ_END

////////////////////////////////////////////////////////////////

PyObject* PyInterval::create(const interval_t& i, pos_t anchor, pos_t anchor_offset)
{
	// otherwise Python 3 will give `SystemError: <built-in function len> returned NULL without setting an error`
	GK_CHECK2(i.size() >= 0, value, "Intervals require a non-negative length.");

	auto* r          = (PyInterval*)PyInterval::DefaultType->tp_new(PyInterval::DefaultType, nullptr, nullptr);
	ainterval_t& val = r->value();
	val = {i, anchor, anchor_offset};
	return (PyObject*)r;
}

END_NAMESPACE_GK
