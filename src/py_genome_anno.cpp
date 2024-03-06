/*
Copyright (C) 2016-2023 Deep Genomics Inc. All Rights Reserved.
*/
#include "py_util.h" // Include first to avoid POSIX warnings due to Python header inclusion order
#include "py_genome_anno.h"

#include "genome_anno.h"
#include "py_genome.h"
#include <structmember.h>

BEGIN_NAMESPACE_GK

// anno table cannot be closed (not necessary unless a process cycle through many annos)
void validate_anno_table(const PyAsPtrSource* self) {}

#define GKPY_GENOME_ANNO_TABLE_BEGIN(name, table_name) \
	GKPY_NEW_OWNED_BEGIN(name##Table, GenomeAnno) \
		self->table = &((PyGenomeAnno*)self->owner)->anno->table_name(); \
		self->set_validator(&validate_anno_table); \
		new (&self->index_by_id) decltype(self->index_by_id); \
	GKPY_NEW_OWNED_END \
\
	GKPY_DEALLOC_OWNED_BEGIN(name##Table) \
		using Object = decltype(self->index_by_id); \
		self->index_by_id.~Object(); \
	GKPY_DEALLOC_OWNED_END \
\
	GKPY_SUBTYPEOBJ_BEGIN(name##Table) \
		tp_flags |= Py_TPFLAGS_BASETYPE; /* Allow Py_TPFLAGS_HAVE_GC to be inherited so that tp_traverse/tp_clear get inherited too */ \
		tp_new = Py##name##Table_New; \
		tp_dealloc = Py##name##Table_Dealloc; \

#define GKPY_GENOME_ANNO_TABLE_END GKPY_SUBTYPEOBJ_END

PyObject* g_biotype_as_pystring[num_biotype];

void Init_GenomeAnno_PyStrings()
{
	for (biotype_t i = 0; i < num_biotype; ++i)
		g_biotype_as_pystring[i] = PyString_FromString(biotype_as_cstr(i)); // new ref
}

PyObject* PyInt_FromElevel(elevel_t elevel)
{
	if (elevel == invalid_elevel)
		GKPY_RETURN_NONE;
	return PyInt_FromLong(1+elevel);
}

PyObject* PyInt_FromTsl(tsl_t tsl)
{
	if (tsl == invalid_tsl)
		GKPY_RETURN_NONE;
	return PyInt_FromLong(1+tsl);
}

//////////////////////////////////////////////////////////////////////

GKPY_NEW_OWNED_BEGIN(GenomeAnno, Genome)
	const auto& genome = ((PyGenome*)self->owner)->genome;
	try {
		self->anno = &genome.anno();
		self->anno->ensure_open();  // Open on-demand; our construction signals demand
	}
	GK_RETHROW("For: Genome(\"{}\")", genome.config());

	// Construct sub-objects
	GKPY_TENTATIVE_CONSTRUCT(genes, GeneTable, "(O)", self);
	GKPY_TENTATIVE_CONSTRUCT(trans, TranTable, "(O)", self);
	GKPY_TENTATIVE_CONSTRUCT(exons, ExonTable, "(O)", self);
	GKPY_TENTATIVE_CONSTRUCT(intrs, IntrTable, "(O)", self);
	GKPY_TENTATIVE_CONSTRUCT(cdss,  CdsTable,  "(O)", self);
	GKPY_TENTATIVE_CONSTRUCT(utr5s, UtrTable,  "(Os)", self, "utr5");
	GKPY_TENTATIVE_CONSTRUCT(utr3s, UtrTable,  "(Os)", self, "utr3");
	// <--- INSERT SUB-OBJECT CONSTRUCTIONS HERE

	// Allow sub-objects to stay constructed on return
	GK_FINALIZE_CONSTRUCT(genes);
	GK_FINALIZE_CONSTRUCT(trans);
	GK_FINALIZE_CONSTRUCT(exons);
	GK_FINALIZE_CONSTRUCT(intrs);
	GK_FINALIZE_CONSTRUCT(cdss);
	GK_FINALIZE_CONSTRUCT(utr5s);
	GK_FINALIZE_CONSTRUCT(utr3s);
	// <--- INSERT SUB-OBJECT FINALIZATIONS HERE, IF CONSTRUCTED ABOVE
GKPY_NEW_OWNED_END

GKPY_TRAVERSE_OWNED_BEGIN(GenomeAnno)
	GKPY_VISIT(genes);
	GKPY_VISIT(trans);
	GKPY_VISIT(exons);
	GKPY_VISIT(intrs);
	GKPY_VISIT(cdss);
	GKPY_VISIT(utr5s);
	GKPY_VISIT(utr3s);
	// <--- INSERT SUB-OBJECT VISITS HERE
GKPY_TRAVERSE_OWNED_END

GKPY_CLEAR_OWNED_BEGIN(GenomeAnno)
	GKPY_CLEAR(genes);
	GKPY_CLEAR(trans);
	GKPY_CLEAR(exons);
	GKPY_CLEAR(intrs);
	GKPY_CLEAR(cdss);
	GKPY_CLEAR(utr5s);
	GKPY_CLEAR(utr3s);
	// <--- INSERT SUB-OBJECT CLEARS HERE
GKPY_CLEAR_OWNED_END

GKPY_DEALLOC_OWNED_BEGIN(GenomeAnno)
GKPY_DEALLOC_OWNED_END

GKPY_GETATTRO_BEGIN(GenomeAnno)
	GKPY_GETATTR_CASE("trans")    { GKPY_RETURN_INCREF(self->trans); }
	GKPY_GETATTR_CASE("filename") { return PyString_FromSV(self->anno->source()); }
	return PyObject_GenericGetAttr(selfo, attro);
GKPY_GETATTRO_END

GKPY_SETATTRO_BEGIN(GenomeAnno)
	GKPY_SETATTR_READONLY("trans");
	GKPY_SETATTR_READONLY("filename");
	return PyObject_GenericSetAttr(selfo, attro, value);
GKPY_SETATTRO_END

GKPY_OMETHOD_BEGIN(GenomeAnno, build_gencode)
	const char*  infile;
	const char*  outfile;
	PyObject*    refg{};
	static char* kwlist[] = {"infile", "outfile", "reference_genome", nullptr};
	if (!PyArg_ParseTupleAndKeywords(args, kwds, "ssO", kwlist, &infile, &outfile, &refg))
		return nullptr;

	auto files = genome_anno::build_gencode(infile, outfile, as_genome(refg));
	auto list  = PyList_New(std::size(files));
	for (size_t i = 0; i < std::size(files); ++i) { PyList_SET_ITEM(list, i, PyString_FromSV(files[i])); }
	return list;
GKPY_OMETHOD_END

GKPY_OMETHOD_BEGIN(GenomeAnno, build_ucsc_refseq)
	const char* ucsc_db_dir;
	const char* outfile;
	PyObject*   refg{};
	static char* kwlist[] = {"ucsc_db_dir", "outfile", "reference_genome", nullptr};
	if (!PyArg_ParseTupleAndKeywords(args, kwds, "ssO", kwlist, &ucsc_db_dir, &outfile, &refg))
		return nullptr;

	auto files = genome_anno::build_ucsc_refseq(ucsc_db_dir, outfile, as_genome(refg));
	auto list  = PyList_New(std::size(files));
	for (size_t i = 0; i < std::size(files); ++i) { PyList_SET_ITEM(list, i, PyString_FromSV(files[i])); }
	return list;
GKPY_OMETHOD_END

GKPY_OMETHOD_BEGIN(GenomeAnno, build_ncbi_refseq)
	const char*  infile;
	const char*  outfile;
	PyObject*    refg{};
	static char* kwlist[] = {"infile", "outfile", "reference_genome", nullptr};
	if (!PyArg_ParseTupleAndKeywords(args, kwds, "ssO", kwlist, &infile, &outfile, &refg))
		return nullptr;

	auto files = genome_anno::build_ncbi_refseq(infile, outfile, as_genome(refg));
	auto list  = PyList_New(std::size(files));
	for (size_t i = 0; i < std::size(files); ++i) { PyList_SET_ITEM(list, i, PyString_FromSV(files[i])); }
	return list;
GKPY_OMETHOD_END

static PyObject* PyGenomeAnno_binary_version(PyObject* cls, PyObject* args)
{
	if (!PyArg_UnpackTuple(args, "", 0, 0))
		return nullptr;
	return PyInt_FromLong((long)genome_anno::binary_version());
}

GKPY_MEMBERS_BEGIN(GenomeAnno)
	GKPY_MEMBER_OBJECT(GenomeAnno, genes, "genes",       nullptr)
	GKPY_MEMBER_OBJECT(GenomeAnno, trans, "transcripts", nullptr)
	GKPY_MEMBER_OBJECT(GenomeAnno, exons, "exons",       nullptr)
	GKPY_MEMBER_OBJECT(GenomeAnno, intrs, "introns",     nullptr)
	GKPY_MEMBER_OBJECT(GenomeAnno, cdss, "cdss",         nullptr)
	GKPY_MEMBER_OBJECT(GenomeAnno, utr5s, "utr5s",       nullptr)
	GKPY_MEMBER_OBJECT(GenomeAnno, utr3s, "utr3s",       nullptr)
GKPY_MEMBERS_END

GKPY_METHODS_BEGIN(GenomeAnno)
	GKPY_METHOD_ENTRY(GenomeAnno, binary_version,    METH_VARARGS | METH_KEYWORDS | METH_STATIC, nullptr)
	GKPY_METHOD_ENTRY(GenomeAnno, build_gencode,     METH_VARARGS | METH_KEYWORDS | METH_STATIC, nullptr)
	GKPY_METHOD_ENTRY(GenomeAnno, build_ucsc_refseq, METH_VARARGS | METH_KEYWORDS | METH_STATIC, nullptr)
	GKPY_METHOD_ENTRY(GenomeAnno, build_ncbi_refseq, METH_VARARGS | METH_KEYWORDS | METH_STATIC, nullptr)
GKPY_METHODS_END

GKPY_TYPEOBJ_BEGIN(GenomeAnno)
	tp_flags |= Py_TPFLAGS_BASETYPE | Py_TPFLAGS_HAVE_GC;
	tp_new = PyGenomeAnno_New;
	tp_dealloc = PyGenomeAnno_Dealloc;
	tp_getattro = PyGenomeAnno_GetAttro;
	tp_setattro = PyGenomeAnno_SetAttro;
	tp_traverse = PyGenomeAnno_Traverse;
	tp_clear = PyGenomeAnno_Clear;
	tp_members = PyGenomeAnno_Members;
	tp_methods = PyGenomeAnno_Methods;
GKPY_TYPEOBJ_END

///////////////////////////////////////////////////////////////////////////

template <typename T>
int PyGenomeAnnoTable_Traverse(PyObject* selfo, visitproc visit, void* arg)
{
	GKPY_TRY
	auto* self = (PyGenomeAnnoTable<T>*)selfo;
	GKPY_VISIT(owner);
	return 0;
	GKPY_CATCH_RETURN_VALUE(-1)
}

template <typename T>
int PyGenomeAnnoTable_Clear(PyObject* selfo)
{
	/* GKPY_CLEAR(owner) -- Deliberately do not clear owner until dealloc. Ensures owner is deallocated last. */
	return 0;
}

GKPY_TEMPLATE_SUBTYPEOBJ_BEGIN(template <typename T>, GenomeAnnoTable<T>)
	tp_flags |= Py_TPFLAGS_BASETYPE | Py_TPFLAGS_HAVE_GC; /* need GC because owner and self refer to each other */
	tp_traverse = PyGenomeAnnoTable_Traverse<T>;
	tp_clear = PyGenomeAnnoTable_Clear<T>;
GKPY_TEMPLATE_SUBTYPEOBJ_END

template <typename T> // T = PyGene, PyTran -- anything that when unpacked has an 'id' member.
PyObject* PyGenomeAnnoTable_GetSubscript_ByID(PyObject* selfo, PyObject* key)
{
	auto* self = (typename T::Table*)selfo;
	if (PyInt_Check(key)) {
		auto i = (Py_ssize_t)PyInt_AS_LONG(key);
		if (i < 0)
			i = PyTable_Len<T>(self) + i;  // Python-style negative indexing
		return PyTable_GetItem<T>(selfo, i);
	} else if (PyString_Check(key)) {
		if (self->index_by_id.empty()) {
			auto& table       = *self->table;
			// return the first found
			for (index_t i = 0; i < table.size(); ++i) {  // update all versioned prefixes
				typename T::unpacked_value v{table[i], table};
				for (const char* end = v.id; end != nullptr; end = strchr(end + 1, '.')) {
					self->index_by_id.emplace(string(v.id, end), i);
				}
				self->index_by_id.emplace(v.id, i);  // full id
			}
		}

		// Allow gene/transcript lookup by ID string
		auto found = self->index_by_id.find(PyString_AS_STRING(key));
		if (found != self->index_by_id.end())
			return PyTable_GetItem<T>(selfo, found->second);
	}
	PyErr_SetObject(PyExc_KeyError, key);
	return nullptr;
}

////////////////////////////////////////////////////////////
// Gene, GeneTable
////////////////////////////////////////////////////////////

GKPY_GETATTRO_BEGIN(Gene)
	// ******** UPDATE genome_kit.Gene MOCK ATTRIBUTES TO MATCH ANY CHANGES HERE ********
	GKPY_GETATTR_CASE("type")              { return PyString_FromBioType(self->value().type); }
	GKPY_GETATTR_CASE("level")             { return PyInt_FromElevel(self->value().level); }
	GKPY_GETATTR_CASE("trans")             { return PyList_FromSizedElemRange<PyTran>(((PyGenomeAnno*)((PyGeneTable*)self->source())->owner)->trans, self->unpack().trans); }
	GKPY_GETATTR_CASE("transcripts")       { return PyList_FromSizedElemRange<PyTran>(((PyGenomeAnno*)((PyGeneTable*)self->source())->owner)->trans, self->unpack().trans); }
	GKPY_GETATTR_CASE("name")              { return PyString_FromString(self->unpack().name); }
	GKPY_GETATTR_CASE("id")                { return PyString_FromString(self->unpack().id); }
	GKPY_GETATTR_CASE("interval")          { return PyInterval::create(self->value()); }
	GKPY_GETATTR_CASE("annotation_genome") { GKPY_RETURN_INCREF(self->anno().owner); }
	return PyGene::base_type::DefaultType->tp_getattro(selfo, attro);
GKPY_GETATTRO_END

GKPY_SETATTRO_BEGIN(Gene)
	GKPY_SETATTR_READONLY("type");
	GKPY_SETATTR_READONLY("level");
	GKPY_SETATTR_READONLY("trans");
	GKPY_SETATTR_READONLY("transcripts");
	GKPY_SETATTR_READONLY("name");
	GKPY_SETATTR_READONLY("id");
	GKPY_SETATTR_READONLY("interval");
	GKPY_SETATTR_READONLY("annotation_genome");
	return PyGene::base_type::DefaultType->tp_setattro(selfo, attro, value);
GKPY_SETATTRO_END

GKPY_METHOD_BEGIN_VARARGS(Gene, _setstate)
	PyObject* table = nullptr;
	int index       = 0;
	if (!PyArg_ParseTuple(args, "Oi", &table, &index))
		return nullptr;

	self->as_ptr = (void*)(&PyGeneTable::value(table)[index]);
	Py_DecRef(self->source());
	self->source() = scast<PyAsPtrSource*>(table);
	Py_IncRef(self->source());
	GKPY_RETURN_NONE;
GKPY_METHOD_END

GKPY_METHODS_BEGIN(Gene)
	GKPY_METHOD_ENTRY(Gene, _setstate, METH_VARARGS, nullptr)
GKPY_METHODS_END

GKPY_VALUE_SUBTYPEOBJ_BEGIN(Gene, as_ptr) // as_ptr because we always point to an entry in a gene_table instance
	tp_flags |= Py_TPFLAGS_BASETYPE;
	tp_init = gkpy_invalid_init;
	tp_getattro = PyGene_GetAttro;
	tp_setattro = PyGene_SetAttro;
	tp_richcompare = PyGenericValue_RichCompare<PyGene>;
	tp_hash = PyGenericValue_Hash<PyGene>;
	tp_methods = PyGene_Methods;
GKPY_VALUE_SUBTYPEOBJ_END

GKPY_GENOME_ANNO_TABLE_BEGIN(Gene, genes)
	mp_subscript = PyGenomeAnnoTable_GetSubscript_ByID<PyGene>;
GKPY_GENOME_ANNO_TABLE_END

////////////////////////////////////////////////////////////
// Tran, TranTable
////////////////////////////////////////////////////////////

GKPY_GETATTRO_BEGIN(Tran)
	GKPY_GETATTR_CASE("tsl")               { return PyInt_FromTsl(self->value().tsl); }
	GKPY_GETATTR_CASE("type")              { return PyString_FromBioType(self->value().type); }
	GKPY_GETATTR_CASE("exons")             { return PyList_FromSizedElemRange<PyExon>(((PyGenomeAnno*)((PyTranTable*)self->source())->owner)->exons, self->unpack().exons); }
	GKPY_GETATTR_CASE("intrs")             { return PyList_FromSizedElemRange<PyIntr>(((PyGenomeAnno*)((PyTranTable*)self->source())->owner)->intrs, self->unpack().intrs); }
	GKPY_GETATTR_CASE("introns")           { return PyList_FromSizedElemRange<PyIntr>(((PyGenomeAnno*)((PyTranTable*)self->source())->owner)->intrs, self->unpack().intrs); }
	GKPY_GETATTR_CASE("cdss")              { return PyList_FromSizedElemRange<PyCds>(((PyGenomeAnno*)((PyTranTable*)self->source())->owner)->cdss, self->unpack().cdss); }
	GKPY_GETATTR_CASE("utr5s")             { return PyList_FromSizedElemRange<PyUtr>(((PyGenomeAnno*)((PyTranTable*)self->source())->owner)->utr5s, self->unpack().utr5s); }
	GKPY_GETATTR_CASE("utr3s")             { return PyList_FromSizedElemRange<PyUtr>(((PyGenomeAnno*)((PyTranTable*)self->source())->owner)->utr3s, self->unpack().utr3s); }
	GKPY_GETATTR_CASE("gene")              { return PyTable_GetItem<PyGene>(((PyGenomeAnno*)((PyTranTable*)self->source())->owner)->genes, (Py_ssize_t)self->value().gene); }
	GKPY_GETATTR_CASE("level")             { return PyInt_FromElevel(self->value().level); }
	GKPY_GETATTR_CASE("interval")          { return PyInterval::create(self->value()); }
	GKPY_GETATTR_CASE("id")                { return PyString_FromString(self->unpack().id); }
	GKPY_GETATTR_CASE("ccds_id")           { return PyString_FromNonEmptyString(self->unpack().ccds_id); }
	GKPY_GETATTR_CASE("protein_id")        { return PyString_FromNonEmptyString(self->unpack().protein_id); }
	GKPY_GETATTR_CASE("product")           { return PyString_FromNonEmptyString(self->unpack().product); }
	GKPY_GETATTR_CASE("annotation_genome") { GKPY_RETURN_INCREF(self->anno().owner); }
	return PyTran::base_type::DefaultType->tp_getattro(selfo, attro);
GKPY_GETATTRO_END

GKPY_SETATTRO_BEGIN(Tran)
	GKPY_SETATTR_READONLY("tsl");
	GKPY_SETATTR_READONLY("type");
	GKPY_SETATTR_READONLY("exons");
	GKPY_SETATTR_READONLY("intrs");
	GKPY_SETATTR_READONLY("introns");
	GKPY_SETATTR_READONLY("cdss");
	GKPY_SETATTR_READONLY("utr5s");
	GKPY_SETATTR_READONLY("utr3s");
	GKPY_SETATTR_READONLY("gene");
	GKPY_SETATTR_READONLY("level");
	GKPY_SETATTR_READONLY("interval");
	GKPY_SETATTR_READONLY("id");
	GKPY_SETATTR_READONLY("ccds_id");
	GKPY_SETATTR_READONLY("protein_id");
	GKPY_SETATTR_READONLY("product");
	GKPY_SETATTR_READONLY("annotation_genome");
	return PyTran::base_type::DefaultType->tp_setattro(selfo, attro, value);
GKPY_SETATTRO_END

GKPY_METHOD_BEGIN_VARARGS(Tran, _setstate)
	PyObject* table = nullptr;
	int index       = 0;
	if (!PyArg_ParseTuple(args, "Oi", &table, &index))
		return nullptr;

	self->as_ptr = (void*)(&PyTranTable::value(table)[index]);
	Py_DecRef(self->source());
	self->source() = scast<PyAsPtrSource*>(table);
	Py_IncRef(self->source());
	GKPY_RETURN_NONE;
GKPY_METHOD_END

GKPY_METHODS_BEGIN(Tran)
	GKPY_METHOD_ENTRY(Tran, _setstate, METH_VARARGS, nullptr)
GKPY_METHODS_END

GKPY_VALUE_SUBTYPEOBJ_BEGIN(Tran, as_ptr) // as_ptr because we always point to an entry in a gene_table instance
	tp_flags |= Py_TPFLAGS_BASETYPE;
	tp_init = gkpy_invalid_init;
	tp_getattro = PyTran_GetAttro;
	tp_setattro = PyTran_SetAttro;
	tp_richcompare = PyGenericValue_RichCompare<PyTran>;
	tp_hash = PyGenericValue_Hash<PyTran>;
	tp_methods = PyTran_Methods;
GKPY_VALUE_SUBTYPEOBJ_END

GKPY_GENOME_ANNO_TABLE_BEGIN(Tran, trans)
	mp_subscript = PyGenomeAnnoTable_GetSubscript_ByID<PyTran>;
GKPY_GENOME_ANNO_TABLE_END

////////////////////////////////////////////////////////////
// Exon, ExonTable
////////////////////////////////////////////////////////////

GKPY_GETATTRO_BEGIN(Exon)
	GKPY_GETATTR_CASE("tran")        { return PyTable_GetItem<PyTran>(self->anno().trans, (Py_ssize_t)self->value().tran); }
	GKPY_GETATTR_CASE("transcript")  { return PyTable_GetItem<PyTran>(self->anno().trans, (Py_ssize_t)self->value().tran); }
	GKPY_GETATTR_CASE("index")       { return PyInt_FromLong((long)self->value().index); }
	GKPY_GETATTR_CASE("interval")    { return PyInterval::create(self->value()); }
	GKPY_GETATTR_CASE("prev_exon")   { return PyTable_CreateItem<PyExon>(self->source(), get_prev(&self->value())); }
	GKPY_GETATTR_CASE("next_exon")   { return PyTable_CreateItem<PyExon>(self->source(), get_next(&self->value())); }
	GKPY_GETATTR_CASE("prev_intron") { return PyTable_CreateItem<PyIntr>(self->anno().intrs, get_prev_intr(&self->value(), *self->anno().anno)); }
	GKPY_GETATTR_CASE("next_intron") { return PyTable_CreateItem<PyIntr>(self->anno().intrs, get_next_intr(&self->value(), *self->anno().anno)); }
	GKPY_GETATTR_CASE("id")          { return PyString_FromNonEmptyString(self->unpack().id); }
	GKPY_GETATTR_CASE("cds")         {
		index_t cds = self->value().cds;
		if (cds == invalid_index)
			GKPY_RETURN_NONE;
		return PyTable_GetItem<PyCds>(((PyGenomeAnno*)((PyExonTable*)self->source())->owner)->cdss, (Py_ssize_t)cds);
	}
	GKPY_GETATTR_CASE("utr5")         {
		index_t utr5 = self->value().utr5;
		if (utr5 == invalid_index)
			GKPY_RETURN_NONE;
		return PyTable_GetItem<PyUtr>(((PyGenomeAnno*)((PyExonTable*)self->source())->owner)->utr5s, (Py_ssize_t)utr5);
	}
	GKPY_GETATTR_CASE("utr3")         {
		index_t utr3 = self->value().utr3;
		if (utr3 == invalid_index)
			GKPY_RETURN_NONE;
		return PyTable_GetItem<PyUtr>(((PyGenomeAnno*)((PyExonTable*)self->source())->owner)->utr3s, (Py_ssize_t)utr3);
	}
	GKPY_GETATTR_CASE("annotation_genome") { GKPY_RETURN_INCREF(self->anno().owner); }
	return PyExon::base_type::DefaultType->tp_getattro(selfo, attro);
GKPY_GETATTRO_END

GKPY_SETATTRO_BEGIN(Exon)
	GKPY_SETATTR_READONLY("tran");
	GKPY_SETATTR_READONLY("transcript");
	GKPY_SETATTR_READONLY("index");
	GKPY_SETATTR_READONLY("interval");
	GKPY_SETATTR_READONLY("prev_exon");
	GKPY_SETATTR_READONLY("next_exon");
	GKPY_SETATTR_READONLY("prev_intron");
	GKPY_SETATTR_READONLY("next_intron");
	GKPY_SETATTR_READONLY("id");
	GKPY_SETATTR_READONLY("cds");
	GKPY_SETATTR_READONLY("utr5");
	GKPY_SETATTR_READONLY("utr3");
	GKPY_SETATTR_READONLY("annotation_genome");
	return PyExon::base_type::DefaultType->tp_setattro(selfo, attro, value);
GKPY_SETATTRO_END

GKPY_METHOD_BEGIN_VARARGS(Exon, _setstate)
	PyObject* table = nullptr;
	int index       = 0;
	if (!PyArg_ParseTuple(args, "Oi", &table, &index))
		return nullptr;

	self->as_ptr = (void*)(&PyExonTable::value(table)[index]);
	Py_DecRef(self->source());
	self->source() = scast<PyAsPtrSource*>(table);
	Py_IncRef(self->source());
	GKPY_RETURN_NONE;
GKPY_METHOD_END

GKPY_METHODS_BEGIN(Exon)
	GKPY_METHOD_ENTRY(Exon, _setstate, METH_VARARGS, nullptr)
GKPY_METHODS_END

GKPY_VALUE_SUBTYPEOBJ_BEGIN(Exon, as_ptr) // as_ptr because we always point to an entry in a gene_table instance
	tp_flags |= Py_TPFLAGS_BASETYPE;
	tp_init = gkpy_invalid_init;
	tp_getattro = PyExon_GetAttro;
	tp_setattro = PyExon_SetAttro;
	tp_richcompare = PyGenericValue_RichCompare<PyExon>;
	tp_hash = PyGenericValue_Hash<PyExon>;
	tp_methods = PyExon_Methods;
GKPY_VALUE_SUBTYPEOBJ_END

GKPY_GENOME_ANNO_TABLE_BEGIN(Exon, exons)
GKPY_GENOME_ANNO_TABLE_END

////////////////////////////////////////////////////////////
// Intr, IntrTable
////////////////////////////////////////////////////////////

GKPY_GETATTRO_BEGIN(Intr)
	GKPY_GETATTR_CASE("tran")              { return PyTable_GetItem<PyTran>(((PyGenomeAnno*)((PyIntrTable*)self->source())->owner)->trans, (Py_ssize_t)self->value().tran); }
	GKPY_GETATTR_CASE("transcript")        { return PyTable_GetItem<PyTran>(((PyGenomeAnno*)((PyIntrTable*)self->source())->owner)->trans, (Py_ssize_t)self->value().tran); }
	GKPY_GETATTR_CASE("index")             { return PyInt_FromLong((long)self->value().index); }
	GKPY_GETATTR_CASE("interval")          { return PyInterval::create(self->value()); }
	GKPY_GETATTR_CASE("prev_intron")       { return PyTable_CreateItem<PyIntr>(self->source(), get_prev(&self->value())); }
	GKPY_GETATTR_CASE("next_intron")       { return PyTable_CreateItem<PyIntr>(self->source(), get_next(&self->value())); }
	GKPY_GETATTR_CASE("prev_exon")         { return PyTable_CreateItem<PyExon>(self->anno().exons, get_prev_exon(&self->value(), *self->anno().anno)); }
	GKPY_GETATTR_CASE("next_exon")         { return PyTable_CreateItem<PyExon>(self->anno().exons, get_next_exon(&self->value(), *self->anno().anno)); }
	GKPY_GETATTR_CASE("annotation_genome") { GKPY_RETURN_INCREF(self->anno().owner); }
	return PyIntr::base_type::DefaultType->tp_getattro(selfo, attro);
GKPY_GETATTRO_END

GKPY_SETATTRO_BEGIN(Intr)
	GKPY_SETATTR_READONLY("tran");
	GKPY_SETATTR_READONLY("transcript");
	GKPY_SETATTR_READONLY("index");
	GKPY_SETATTR_READONLY("interval");
	GKPY_SETATTR_READONLY("next_intron");
	GKPY_SETATTR_READONLY("prev_intron");
	GKPY_SETATTR_READONLY("next_exon");
	GKPY_SETATTR_READONLY("prev_exon");
	GKPY_SETATTR_READONLY("annotation_genome");
	return PyIntr::base_type::DefaultType->tp_setattro(selfo, attro, value);
GKPY_SETATTRO_END

GKPY_METHOD_BEGIN_VARARGS(Intr, _setstate)
	PyObject* table = nullptr;
	int index       = 0;
	if (!PyArg_ParseTuple(args, "Oi", &table, &index))
		return nullptr;

	self->as_ptr = (void*)(&PyIntrTable::value(table)[index]);
	Py_DecRef(self->source());
	self->source() = scast<PyAsPtrSource*>(table);
	Py_IncRef(self->source());
	GKPY_RETURN_NONE;
GKPY_METHOD_END

GKPY_METHODS_BEGIN(Intr)
	GKPY_METHOD_ENTRY(Intr, _setstate, METH_VARARGS, nullptr)
GKPY_METHODS_END

GKPY_VALUE_SUBTYPEOBJ_BEGIN(Intr, as_ptr) // as_ptr because we always point to an entry in a gene_table instance
	tp_flags |= Py_TPFLAGS_BASETYPE;
	tp_init = gkpy_invalid_init;
	tp_getattro = PyIntr_GetAttro;
	tp_setattro = PyIntr_SetAttro;
	tp_richcompare = PyGenericValue_RichCompare<PyIntr>;
	tp_hash = PyGenericValue_Hash<PyIntr>;
	tp_methods = PyIntr_Methods;
GKPY_VALUE_SUBTYPEOBJ_END

GKPY_GENOME_ANNO_TABLE_BEGIN(Intr, intrs)
GKPY_GENOME_ANNO_TABLE_END

////////////////////////////////////////////////////////////
// Cds, CdsTable
////////////////////////////////////////////////////////////

GKPY_GETATTRO_BEGIN(Cds)
	GKPY_GETATTR_CASE("exon")              { return PyTable_CreateItem<PyExon>(((PyGenomeAnno*)((PyCdsTable*)self->source())->owner)->exons, self->unpack().exon); }
	GKPY_GETATTR_CASE("tran")              { return PyTable_GetItem<PyTran>(((PyGenomeAnno*)((PyCdsTable*)self->source())->owner)->trans, (Py_ssize_t)self->value().tran); }
	GKPY_GETATTR_CASE("transcript")        { return PyTable_GetItem<PyTran>(((PyGenomeAnno*)((PyCdsTable*)self->source())->owner)->trans, (Py_ssize_t)self->value().tran); }
	GKPY_GETATTR_CASE("phase")             { return PyInt_FromLong((long)self->value().phase); }
	GKPY_GETATTR_CASE("interval")          { return PyInterval::create(self->value()); }
	GKPY_GETATTR_CASE("prev_cds")          { return PyTable_CreateItem<PyCds>(self->source(), get_prev(&self->value())); }
	GKPY_GETATTR_CASE("next_cds")          { return PyTable_CreateItem<PyCds>(self->source(), get_next(&self->value())); }
	GKPY_GETATTR_CASE("annotation_genome") { GKPY_RETURN_INCREF(self->anno().owner); }
	return PyCds::base_type::DefaultType->tp_getattro(selfo, attro);
GKPY_GETATTRO_END

GKPY_SETATTRO_BEGIN(Cds)
	GKPY_SETATTR_READONLY("exon");
	GKPY_SETATTR_READONLY("tran");
	GKPY_SETATTR_READONLY("transcript");
	GKPY_SETATTR_READONLY("phase");
	GKPY_SETATTR_READONLY("interval");
	GKPY_SETATTR_READONLY("prev_cds");
	GKPY_SETATTR_READONLY("next_cds");
	GKPY_SETATTR_READONLY("annotation_genome");
	return PyCds::base_type::DefaultType->tp_setattro(selfo, attro, value);
GKPY_SETATTRO_END

GKPY_METHOD_BEGIN_VARARGS(Cds, _setstate)
	PyObject* table = nullptr;
	int index       = 0;
	if (!PyArg_ParseTuple(args, "Oi", &table, &index))
		return nullptr;

	self->as_ptr = (void*)(&PyCdsTable::value(table)[index]);
	Py_DecRef(self->source());
	self->source() = scast<PyAsPtrSource*>(table);
	Py_IncRef(self->source());
	GKPY_RETURN_NONE;
GKPY_METHOD_END

GKPY_METHODS_BEGIN(Cds)
	GKPY_METHOD_ENTRY(Cds, _setstate, METH_VARARGS, nullptr)
GKPY_METHODS_END

GKPY_VALUE_SUBTYPEOBJ_BEGIN(Cds, as_ptr) // as_ptr because we always point to an entry in a gene_table instance
	tp_flags |= Py_TPFLAGS_BASETYPE;
	tp_init = gkpy_invalid_init;
	tp_getattro = PyCds_GetAttro;
	tp_setattro = PyCds_SetAttro;
	tp_richcompare = PyGenericValue_RichCompare<PyCds>;
	tp_hash = PyGenericValue_Hash<PyCds>;
	tp_methods = PyCds_Methods;
GKPY_VALUE_SUBTYPEOBJ_END

GKPY_GENOME_ANNO_TABLE_BEGIN(Cds, cdss)
GKPY_GENOME_ANNO_TABLE_END

//////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////
// Utr, UtrTable
////////////////////////////////////////////////////////////

GKPY_GETATTRO_BEGIN(Utr)
	GKPY_GETATTR_CASE("exon")              { return PyTable_CreateItem<PyExon>(((PyGenomeAnno*)((PyUtrTable*)self->source())->owner)->exons, self->unpack().exon); }
	GKPY_GETATTR_CASE("tran")              { return PyTable_GetItem<PyTran>(((PyGenomeAnno*)((PyUtrTable*)self->source())->owner)->trans, (Py_ssize_t)self->value().tran); }
	GKPY_GETATTR_CASE("transcript")        { return PyTable_GetItem<PyTran>(((PyGenomeAnno*)((PyUtrTable*)self->source())->owner)->trans, (Py_ssize_t)self->value().tran); }
	GKPY_GETATTR_CASE("interval")          { return PyInterval::create(self->value()); }
	GKPY_GETATTR_CASE("prev_utr")          { return PyTable_CreateItem<PyUtr>(self->source(), get_prev(&self->value())); }
	GKPY_GETATTR_CASE("next_utr")          { return PyTable_CreateItem<PyUtr>(self->source(), get_next(&self->value())); }
	GKPY_GETATTR_CASE("annotation_genome") { GKPY_RETURN_INCREF(self->anno().owner); }
		return PyUtr::base_type::DefaultType->tp_getattro(selfo, attro);
GKPY_GETATTRO_END

GKPY_SETATTRO_BEGIN(Utr)
		GKPY_SETATTR_READONLY("exon");
		GKPY_SETATTR_READONLY("tran");
		GKPY_SETATTR_READONLY("transcript");
		GKPY_SETATTR_READONLY("interval");
		GKPY_SETATTR_READONLY("prev_utr");
		GKPY_SETATTR_READONLY("next_utr");
		GKPY_SETATTR_READONLY("annotation_genome");
		return PyUtr::base_type::DefaultType->tp_setattro(selfo, attro, value);
GKPY_SETATTRO_END

GKPY_METHOD_BEGIN_VARARGS(Utr, _setstate)
		PyObject* table = nullptr;
		int index       = 0;
		if (!PyArg_ParseTuple(args, "Oi", &table, &index))
			return nullptr;

		self->as_ptr = (void*)(&PyUtrTable::value(table)[index]);
		Py_DecRef(self->source());
		self->source() = scast<PyAsPtrSource*>(table);
		Py_IncRef(self->source());
		GKPY_RETURN_NONE;
GKPY_METHOD_END

GKPY_METHODS_BEGIN(Utr)
		GKPY_METHOD_ENTRY(Utr, _setstate, METH_VARARGS, nullptr)
GKPY_METHODS_END

GKPY_VALUE_SUBTYPEOBJ_BEGIN(Utr, as_ptr) // as_ptr because we always point to an entry in a gene_table instance
		tp_flags |= Py_TPFLAGS_BASETYPE;
		tp_init = gkpy_invalid_init;
		tp_getattro = PyUtr_GetAttro;
		tp_setattro = PyUtr_SetAttro;
		tp_richcompare = PyGenericValue_RichCompare<PyUtr>;
		tp_hash = PyGenericValue_Hash<PyUtr>;
		tp_methods = PyUtr_Methods;
GKPY_VALUE_SUBTYPEOBJ_END

GKPY_NEW_BEGIN(UtrTable)
		self->owner         = 0;
		const char* str_arg = nullptr;
		if (!PyArg_ParseTuple(args, "O!s", PyGenomeAnno::DefaultType, &self->owner, &str_arg))
			return NULL;
		Py_INCREF(self->owner);
		if (strcmp(str_arg, "utr5") == 0) {
			self->table = &((PyGenomeAnno*)self->owner)->anno->utr5s();
		} else if (strcmp(str_arg, "utr3") == 0) {
			self->table = &((PyGenomeAnno*)self->owner)->anno->utr3s();
		} else {
			GK_THROW(value, "Invalid argument to UtrTable constructor: {}", str_arg);
		}
		self->set_validator(&validate_anno_table);
		new (&self->index_by_id) decltype(self->index_by_id);
GKPY_NEW_OWNED_END

GKPY_DEALLOC_OWNED_BEGIN(UtrTable)
		using Object = decltype(self->index_by_id);
		self->index_by_id.~Object();
GKPY_DEALLOC_OWNED_END

GKPY_SUBTYPEOBJ_BEGIN(UtrTable)
		tp_flags |= Py_TPFLAGS_BASETYPE; /* Allow Py_TPFLAGS_HAVE_GC to be inherited so that tp_traverse/tp_clear get
											inherited too */
		tp_new     = PyUtrTable_New;
		tp_dealloc = PyUtrTable_Dealloc;
GKPY_GENOME_ANNO_TABLE_END

//////////////////////////////////////////////////////////////////////

END_NAMESPACE_GK
