/*
Copyright (C) 2016-2023 Deep Genomics Inc. All Rights Reserved.
*/
// This file is only included when compiled as an executable,
// not when included as a python extension. It's intended to
// run built-in C++ unit tests.
#ifdef _WANT_MAIN

#include "genome_kit.h"
#include "gk_time.h"
#include "format.h"
#include "half.h"
#include "jralign.h"
#include "jrdist.h"
#include "strutil.h"
#include "util.h"
#include <algorithm>
#include <cstdio>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <memory>
#include <set>
#include <stdexcept>
#include <vector>

USING_NAMESPACE_GK

using std::unique_ptr;
using std::ofstream;

string custom_realize_datafile_path(string path)
{
	if (endswith(path, ".dganno")) {
		string newpath = getenv("LOCALAPPDATA");
		newpath += R"(\genome_kit\genome_kit\)";
		newpath += path.substr(0, path.size() - strlen(".dganno"));
		newpath += format(".v{}.dganno", genome_anno::binary_version());
		return newpath;
	}
	return path;
}

void walk_test()
{
	genome_t genome("gencode.v29.mini");

	for (auto& g : genome.genes()) {
		gene_t gene(g, genome.genes());
		println("{}", gene);

		for (auto& t : gene.trans) {
			tran_t tran(t, genome.trans());
			println("{}", tran);

			for (auto& e : tran.exons) {
				exon_t exon(e, genome.exons());
				println("{}", exon);

				if (exon.utr5) {
					utr_t utr5(*exon.utr5, genome.utr5s());
					println("{}", utr5);
				}

				if (exon.cds) {
					cds_t cds(*exon.cds, genome.cdss());
					println("{}", cds);
				}

				if (exon.utr3) {
					utr_t utr3(*exon.utr3, genome.utr3s());
					println("{}", utr3);
				}
			}

			for (auto& i : tran.intrs) {
				intr_t intr(i, genome.intrs());
				println("{}", intr);
			}
		}
	}
}

struct track_unittest_value_range {
	int lo; // first value to generate
	int hi; // last  value to generate (inclusive)
	int n;  // number of values to generate in range [lo, hi]
};

template <typename data_t>
data_t track_unittest_value_gen(int value_index, const track_unittest_value_range& r)
{
	double x = (r.hi-r.lo)*double(value_index % r.n)/max(1, r.n-1) + r.lo;
	return data_t(x);
}

template <> bool track_unittest_value_gen<bool>(int value_index, const track_unittest_value_range& r)
{
	if (r.n == 1)
		return true;  // Special case for m0
	return value_index % 2 != 0; // alternating true/false
}

template <typename data_t>
void test_genome_track_dtype_etype(genome_track::etype_t etype)
{
	print("track_unittest(dtype={}, etype={}) ... ", genome_track::dtype_as_cstr[as_ordinal(dtype_traits<data_t>::dtype)], genome_track::etype_as_cstr[etype]);

	track_unittest_value_range value_ranges[genome_track::num_etype] = {
		{    1,    1,    1 },  // m0   -- generate [1] always
		{    0,    1,    2 },  // u1   -- generate [0, 1] using 1 bit
		{    0,    3,    4 },  // u2   -- generate [0 ... 3] using 2 bits
		{    0,    7,    8 },  // u3   -- generate [0 ... 7] using 3 bits
		{    0,   15,   16 },  // u4   -- generate [0 ... 15] using 4 bits
		{    0,   31,   32 },  // u5   -- generate [0 ... 31] using 5 bits
		{    0,   63,   64 },  // u6   -- generate [0 ... 63] using 6 bits
		{    0,  255,  256 },  // u8   -- generate [0 ... 255] using 8 bits
		{ -128,  127,  256 },  // i8   -- generate [-128 ... 0 ... 127] using 8 bits
		{-1024,  512,    4 },  // f2   -- generate [-1024, -512, +0, +512] using 2 bits, all half_t representable
		{-1024,  768,    8 },  // f3   -- generate [-1024 ... -256, +0, +256 ... +768] using 3 bits, all half_t representable
		{-1024,  896,   16 },  // f4   -- generate [-1024 ... -128, +0, +128 ... +896] using 4 bits, all half_t representable
		{-1024,  960,   32 },  // f5   -- generate [-1024 ...  -64, +0,  +64 ... +960] using 5 bits, all half_t representable
		{-1024,  992,   64 },  // f6   -- generate [-1024 ...  -32, +0,  +32 ... +992] using 6 bits, all half_t representable
		{-1024, 1016,  256 },  // f8   -- generate [-1024 ...   -8, +0,   +8 ... +1016] using 8 bits, all half_t representable
		{-1024, 1024, 2049 },  // f16  -- generate [-1024 ...   -1, +0,   +1 ... +1024] using 10 bits, all half_t representable (10 bit mantissa)
	};

	// Special case: if encoding from data_t==int8_t and etype==u8, then don't
	// try to encode values >127 or else they'll end up being negative (un-encodable as u8)
	if (dtype_traits<data_t>::dtype == genome_track::int8) {
		value_ranges[genome_track::u8].hi = 0;
		value_ranges[genome_track::u8].n  = 128;
	}

	const int max_gap_size = 4;
	const int max_chunk_size = 17;
	const int chrom_size = 61;
	const int max_dim = 8;  // Cover specializations for dim=1..4, for dim divisible by 2,3,4,5, and for the generic code (via dim=7 case)
	const int max_res = 3;  // Cover specializations for 1bp resolution and >1bp resolution
	const int max_query_size = min(chrom_size, max_chunk_size * 3);
	data_t  src_data[2][chrom_size*max_dim];        // [neg_strand/pos_strand][position][dimension]
	data_t  dst_buf[chrom_size*max_dim*max_res+2];  // +2 so that can detect buffer overflow at ends of array, even for largest setting.
	data_t* dst_data = dst_buf+1;

	// Initialize dict for the case etype==f2..8
	float dict[256];
	if (etype >= genome_track::f2 && etype <= genome_track::f8)
		for (int i = 0; i < value_ranges[etype].n; ++i)
			dict[i] = track_unittest_value_gen<float>(i, value_ranges[etype]);

	// Several nested loops for testing different conditions exhaustively
	// (stranded/unstranded, data chunk sizes, sizes of gap between chunks, track dimension, etc)
	const genome_t genome{"hg19.mini"};
	const auto&    chrom_names = genome.chrom_names();
	auto it_chrom   = std::cbegin(chrom_names);
	bool is_mask     = (etype == genome_track::m0);
	for (int stranded = 0; stranded < 2; ++stranded) {
		for (int sparsify = 0; sparsify < 2; ++sparsify) {
			if (sparsify && is_mask)
				continue;
			for (int chunk_size = 1; chunk_size <= max_chunk_size; ++chunk_size) {
				for (int gap_size = 0; gap_size <= max_gap_size; gap_size += 3) {
					auto chrom = it_chrom->first;
					++it_chrom;
					if (it_chrom == std::cend(chrom_names))
						it_chrom = std::cbegin(chrom_names);
					int max_etype_dim = is_mask ? 1 : max_dim;
					for (int dim = 1; dim <= max_etype_dim; ++dim) {
						for (int res = 1; res <= max_res; ++res) {
							// 0 = set_data
							// 1 = set_data_from_wig
							if (is_mask)  // can't load mask from WIG file
								continue;

							const auto strandedness =
								stranded? gk::genome_track::strandedness_t::strand_unaware
										: gk::genome_track::strandedness_t::single_stranded;
							// Check that the current dtype is supported for decoding; if not skip it.
							genome_track::builder builder("track_unittest.dgtrack", etype, strandedness, genome, dim,
														  res);
							if (!builder.supports_dtype(dtype_traits<data_t>::dtype)) {
								print("(skipped; no decoder)\n");
								return; // Oops, we can't actually decode as data_t, so quit this whole thing
							}

							/////////////////////////////////////////////////////
							// 1. Fill src_dat with a mix of default_value (gap_size) and
							//    non-default value (chunk_size)
							/////////////////////////////////////////////////////

							// Choose a default_value (0) that is encodeable by all types, including our dict,
							// so that sparsity can be tested
							auto default_value = (data_t)0;
							// Fill the src_data arrays with default value
							for (int i = 0; i < chrom_size*dim; ++i) {
								src_data[as_ordinal(neg_strand)][i] = default_value;
								src_data[as_ordinal(pos_strand)][i] = default_value;
							}

							// Overwrite parts of the src_data arrays with chunks of position-specific, i.e. position modulo a prime number
							for (int start = 0; start < chrom_size; start += chunk_size+gap_size) {
								for (int i = start; i < min(start+chunk_size, chrom_size); ++i) {
									for (int j = 0; j < dim; ++j) {
										src_data[as_ordinal(neg_strand)][i*dim+j] = track_unittest_value_gen<data_t>(i*dim+j + 7, value_ranges[etype]);  // +7 ensures neg_strand differs from pos_strand
										src_data[as_ordinal(pos_strand)][i*dim+j] = track_unittest_value_gen<data_t>(i*dim+j + 0, value_ranges[etype]);
									}
								}
							}

							// (If stranded mode, then set the data for both strands. Otherwise only set positive strand.)
							vector<interval_t> intervals;
							for (std::underlying_type_t<strand_t> strand = as_ordinal(stranded ? neg_strand : pos_strand); strand < num_strand; strand += 1) {
								// Set the data in contiguous blocks of non-default value.
								// This entails finding [start, end) pairs that contain non-default value.
								for (int start = 0, count = 0;; ++count) {
									// Find the next position where at least ONE element IS NOT default_value
									while (start < chrom_size) {
										int j = 0;
										while (j < dim && src_data[strand][start*dim+j] == default_value)
											++j;
										if (j < dim)
											break;
										start++;
									}
									if (start >= chrom_size)
										break;
									// Find the next position where ALL elements are default_value
									int end = start+1;
									while (end < chrom_size) {
										int j = 0;
										while (j < dim && src_data[strand][end*dim+j] == default_value)
											++j;
										if (j == dim) {
											break;
										}
										end++;
									}

									intervals.push_back(
										interval_t::from_dna0(chrom, start, end, strand_t(strand), genome.refg()));
									start = end+1;
								}
							}

							if (sparsify) {
								// If testing sparsify, fill the middle 3rd of the track with default_value
								// so that the encoder will try to exclude these regons from the encoded track
								for (int i = 0; i < chrom_size; ++i) {
									if (i >= chrom_size/3 && i < 2*chrom_size/3) {
										for (int j = 0; j < dim; ++j) {
											src_data[as_ordinal(neg_strand)][i*dim+j] = default_value;
											src_data[as_ordinal(pos_strand)][i*dim+j] = default_value;
										}
									}
								}
							}

							/////////////////////////////////////////////////////
							// 2. Build a dgtrack from the src_data, creating intervals
							//    for all the contiguous chunks of non-default value
							/////////////////////////////////////////////////////

							// Build the track from src_data
							builder.set_default_value(default_value);
							if (sparsify)
								builder.set_sparsity(4);  // sparsify any default_value block of size >= 4
							// Set dictionary if applicable
							if (etype >= genome_track::f2 && etype <= genome_track::f8)
								builder.set_dict(dict);
							// Load the data into the track
							for (auto& interval : intervals) {
								// if res > 1, we want to still encode all of src_data, so that means
								// we'll need to scale up the intervals by a factor of 'res', and it alsp
								// means we'll be decoding more data than we encoded (since values will get
								// repeated on decoding).
								interval_t scaled =
									interval_t::from_dna0(chrom, interval.start() * res, interval.end() * res,
														  interval.strand, genome.refg());
								builder.set_data(scaled, is_mask ? nullptr : &src_data[(int)interval.strand][interval.start()*dim]);
							}
							builder.finalize();

							/////////////////////////////////////////////////////
							// 2. Load the track from disk and ensure a variety
							//    of queries decodes data matching src_data.
							/////////////////////////////////////////////////////
							genome_track track;
							track.set_source("track_unittest.dgtrack");
							track.open();
							if (track.empty())  // No data was written
								continue;
							GK_ASSERT(track.dim() == dim);
							GK_ASSERT(track.refg() == genome.refg());
							GK_ASSERT(track.stranded() == (stranded == 1));

							///////////////////////////////////////////
                            for (std::underlying_type_t<strand_t> strand = 0; strand < num_strand; strand += 1) {
								for (int query_size = 0; query_size < max_query_size; ++query_size) {
									for (int start = 0; start < chrom_size*res-query_size; ++start) {
										//////////////////////////////////////////////////
										// Finally, decode a specific query interval (start, start+query_size)
										// into  dst_data and compare it to src_data
										//////////////////////////////////////////////////

										// Start by initializing dst_data with a magic number, including the buffer
										// underrun (-1) and overrun (query_size*dim) positions
										const auto unwritten_value = (data_t)123;
										for (int i = -1; i <= query_size*dim; ++i)
											dst_data[i] = unwritten_value;  // magic number should get overwritten

										// Then build the query interval and decode it into dst_data
										interval_t interval = interval_t::from_dna0(chrom, start, start + query_size,
																					strand_t(strand), genome.refg());
										track(interval, dst_data);  // +1 to be able to detect buffer underrun via dst_data[0]

										// Check for buffer underrun/overrun
										GK_ASSERT(dst_data[-1]             == unwritten_value);
										GK_ASSERT(dst_data[query_size*dim] == unwritten_value);

										// Check that each decoded datum in dst_data is equal its corresponding value in src_data.
										// Note that this only works if the encoding format is able to represent the encoded values exactly,
										// so be careful to ensure that's the case when the track was created (e.g. don't try to encode a value
										// to a dict-based format that isn't in the dict).
										for (int i = 0; i < query_size; ++i) {
											for (int j = 0; j < dim; ++j) {
												int d = strand_t(strand) == pos_strand ? i : query_size-i-1;  // Index of dst position to compare, possibly reversed
												int s = (i+start)/res;                // Index of src position to compare
												int src_strand = stranded ? strand : as_ordinal(pos_strand); // If neg_strand and !stranded, then compare to pos_strand data
												data_t decoded = dst_data[d*dim+j];
												data_t encoded = src_data[src_strand][s*dim+j];
												GK_ASSERT(decoded == encoded);
											}
										}
									} // start
								} // query_size
							} // strand
							//////////////////////////////////////////////

						} // res
					} // dim
				} // gap_size
			} // chrom_size
		} // sparsify
	} // stranded

	print("done\n");
}

template <typename data_t>
void test_genome_track_dtype()
{
	test_genome_track_dtype_etype<data_t>(genome_track::m0);
	test_genome_track_dtype_etype<data_t>(genome_track::u1);
	test_genome_track_dtype_etype<data_t>(genome_track::u2);
	test_genome_track_dtype_etype<data_t>(genome_track::u3);
	test_genome_track_dtype_etype<data_t>(genome_track::u4);
	test_genome_track_dtype_etype<data_t>(genome_track::u5);
	test_genome_track_dtype_etype<data_t>(genome_track::u6);
	test_genome_track_dtype_etype<data_t>(genome_track::u8);
	test_genome_track_dtype_etype<data_t>(genome_track::i8);
	test_genome_track_dtype_etype<data_t>(genome_track::f2);
	test_genome_track_dtype_etype<data_t>(genome_track::f3);
	test_genome_track_dtype_etype<data_t>(genome_track::f4);
	test_genome_track_dtype_etype<data_t>(genome_track::f5);
	test_genome_track_dtype_etype<data_t>(genome_track::f6);
	test_genome_track_dtype_etype<data_t>(genome_track::f8);
	test_genome_track_dtype_etype<data_t>(genome_track::f16);
}

void test_genome_track()
{
	test_genome_track_dtype<bool   >();
	test_genome_track_dtype<uint8_t>();
	test_genome_track_dtype<int8_t >();
	test_genome_track_dtype<half_t >();
	test_genome_track_dtype<float  >();
}

void write_ralign_test(const vector<string>& samfiles)
{
	vector<interval_t> exclude;
	const genome_t     genome{"hg19.mini"};
	const auto& h37_names = genome.chrom_names();
	exclude.emplace_back(h37_names.as_chrom("chr11"), 5100000, 5682617, pos_strand,
						 genome.refg());  // whole blood-specific gene
	exclude.emplace_back(h37_names.as_chrom("chr16"), 201431, 232756, pos_strand,
						 genome.refg());  // whole blood-specific gene
	exclude.emplace_back(h37_names.as_chrom("chr7"), 141938359, 142639093, pos_strand,
						 genome.refg());  // pancreas-specific gene
	exclude.emplace_back(h37_names.as_chrom("chr1"), 22297969, 22342285, pos_strand,
						 genome.refg());  // pancreas-specific gene

	print("\n------- BUILD READ ALIGNMENTS (from BAM) -----\n");
	gk::timer_t t; t.tic();
	for (auto& samfile : samfiles) {
		string outfile = samfile + ".jralign";
		junction_read_alignments::builder ralign(outfile.c_str(), genome);
		ralign.set_min_reads(3);
		ralign.set_min_overhang(5);
		for (auto& i : exclude)
			ralign.get_interval_filter().exclude(i);
		ralign.add(samfile.c_str());
		ralign.finalize();
	}
	t.toc();
	print("Created {} junction read alignment files in {:.1f} sec ({:.1f} sec / file)\n", samfiles.size(), t.duration(), t.duration() / samfiles.size());
}

void read_ralign_test(const vector<string>& samfiles)
{
	print("\n------- TRAVERSE READ ALIGNMENTS --------\n");
	gk::timer_t t;
	for (auto& samfile : samfiles) {
		junction_read_alignments jraligns;
		jraligns.set_source(samfile + ".jralign");

		t.tic();
		jraligns.open();
		t.toc();
		print("Opened in {:.05f} sec ({} junctions)\n", t.duration(), jraligns.juncs().size());

		t.tic();
		for (int j = 0; j < jraligns.juncs().size(); ++j) {
			jraligns_t junc(j, jraligns.juncs());
			print("-----------------------------------------------\n");
			print("{} {}-{} ({} reads)\n", get_chrom_names(junc.refg).chrom_as_sv(junc.chrom),
				  junc.start(), junc.end(), junc.num_reads());
			for (unsigned r = 0; r < junc.num_reads(); ++r) {
				jralign_t read = junc[r];
				print("    {} {:2d}|{:2d}\n", strand_as_char(read.strand), read.left, read.right);
			}
		}
		t.toc();
		print("Traversed all reads in {:.04f} sec\n", t.duration());
	}
}

void write_rdist_test(const vector<string>& samfiles)
{
	print("\n------- BUILD READ DISTRIBUTIONS -----\n");
	gk::timer_t t; t.tic();
	for (auto& samfile : samfiles) {
		string infile  = samfile + ".ralign";
		string outfile = samfile + ".rdist";
		read_distributions::builder rdists(outfile.c_str());
		rdists.add(infile.c_str());
		rdists.finalize();
	}
	t.toc();
	print("Created {} junction read distribution files in {:.1f} sec {:.1f} sec / file)\n", samfiles.size(), t.duration(), t.duration() / samfiles.size());
}

void write_pooled_rdist_test(const vector<string>& samfiles)
{
	print("\n------- BUILD READ DISTRIBUTIONS -----\n");
	gk::timer_t t; t.tic();
	read_distributions::builder rdists("pooled.rdist");
	const int reps = 10;
	for (int i = 0; i < reps; ++i) {
		for (auto& samfile : samfiles) {
			string infile  = samfile + ".ralign";
			string outfile = samfile + ".rdist";
			rdists.add(infile.c_str());
		}
	}
	rdists.finalize();
	t.toc();
	print("Created {} junction read distribution files in {:.1f} sec ({:.1f} sec / file)\n", samfiles.size()*reps, t.duration(), t.duration() / samfiles.size() / reps);
}

void read_rdist_test(const vector<string>& samfiles)
{
	print("\n------- TRAVERSE READ DISTRIBUTIONS -----\n");
	gk::timer_t t;
	for (auto& samfile : samfiles) {
		read_distributions rdists;
		rdists.set_source(samfile + ".rdist");

		t.tic();
		rdists.open();
		t.toc();
		print("Opened in {:.05f} sec ({} junctions)\n", t.duration(), rdists.juncs().size());

		t.tic();
		for (int j = 0; j < rdists.juncs().size(); ++j) {
			jrdist_t rdist(j, rdists.juncs());
			print("-----------------------------------------------\n");
			print("{} {}-{} ({} counts, {} reads)\n",
				  get_chrom_names(rdist.refg).chrom_as_sv(rdist.chrom), rdist.start(), rdist.end(),
				  rdist.num_counts(), rdist.num_reads());
			for (int r = 0; r < rdist.num_counts(); ++r) {
				jrcount_t rcount = rdist[r];
				print("    {}\t{:2d}\t{:2d}\n", strand_as_char(rcount.strand), rcount.shift, rcount.count);
			}
		}
		t.toc();
		print("Traversed all read counts in {:.04f} sec\n", t.duration());
	}
}

void read_pooled_rdist_test()
{
	print("\n------- TRAVERSE READ DISTRIBUTIONS -----\n");
	gk::timer_t t;
	read_distributions rdists;
	rdists.set_source("pooled.rdist");

	t.tic();
	rdists.open();
	t.toc();
	print("Opened in {:.05f} sec ({} junctions)\n", t.duration(), rdists.juncs().size());

	t.tic();
	for (int j = 0; j < rdists.juncs().size(); ++j) {
		jrdist_t rdist(j, rdists.juncs());
		print("-----------------------------------------------\n");
		print("{} {}-{} ({} counts, {} reads)\n",
			  get_chrom_names(rdist.refg).chrom_as_sv(rdist.chrom), rdist.start(), rdist.end(),
			  rdist.num_counts(), rdist.num_reads());
		for (int r = 0; r < rdist.num_counts(); ++r) {
			jrcount_t rcount = rdist[r];
			print("    {}\t{:2d}\t{:2d}\n", strand_as_char(rcount.strand), rcount.shift, rcount.count);
		}
	}
	t.toc();
	print("Traversed all read counts in {:.04f} sec\n", t.duration());
}

void clamp_test() {
    half_t a1[] = {1.2, 2.3, 3.4, 4.5, 6.7};
    clamp_track_data(a1, 5, 1, 2.0, 5.0);
    half_t expected[] = {2.0, 2.3, 3.4, 4.5, 5.0};
    GK_ASSERT(std::equal(std::begin(a1), std::end(a1), std::begin(expected)));

    float a2[] = {1.2, 2.3, 3.4, 4.5, 6.7};
    clamp_track_data(a2, 5, 1, 2.0, 5.0);
    GK_ASSERT(std::equal(std::begin(a2), std::end(a2), std::begin(expected)));
}

genome_track::builder build_gtrack_builder(const float max_val, const genome_t& genome) {
	genome_track::builder subject("output.gtrack", genome_track::f8, gk::genome_track::strandedness_t::strand_unaware,
								  genome);
	const float delta = max_val / 255.0f;
    float dict[256];
    for (uint16_t i = 0; i < 256; ++i) {
        dict[i] = (float)i * delta;
    }
    subject.set_dict(dict);
    subject.set_default_value(0);
    subject.set_clamping();

    return subject;
}

std::string get_data_file_path(std::string relative_path) {
    return "." + relative_path;
}

void gtrack_builder_set_data_from_bed_test() {
    const float max_val = 40000.0f;
	const genome_t        genome{"hg19.mini"};
	genome_track::builder subject = build_gtrack_builder(max_val, genome);
	subject.set_data_from_bed(get_data_file_path("/tests/data/test_clamp/test.bed"));
    subject.finalize();

    genome_track track{};
    track.set_source("output.gtrack");
    track.open();
    /**
    This part of the .bed file:
chr1	565661	565662	chr1:565661..565662,+	1000	+
chr1	566953	566954	chr1:566953..566954,-	2	-
chr1	567209	567210	chr1:567209..567210,+	2000	+
chr1	567221	567222	chr1:567221..567222,+	5000	+
chr1	567226	567227	chr1:567226..567227,+	10000	+
chr1	568398	568399	chr1:568398..568399,+	20000	+
chr1	569882	569883	chr1:569882..569883,+	40000	+
chr1	569883	569884	chr1:569883..569884,+	50000	+
chr1	569886	569887	chr1:569886..569887,+	70000	+

    Dict for these values:
        1000  => 941.176514
        2000  => 2039.215698
        5000  => 5019.607910
        10000 => 10039.215820
        20000 => 20078.431641
     */
    const int start = 565661;
    half_t dst[5000] = {0};
	interval_t interval  = interval_t::from_dna0(genome.chrom_names().as_chrom("chr1"), start, start + std::size(dst),
												 strand_t::pos_strand, genome.refg());
	track(interval, dst, track.dtype());
    track.close();

    GK_ASSERT(as_float(dst[0]) == as_float(half_t(941.176514)));
    GK_ASSERT(as_float(dst[567209-start]) == as_float(half_t(2039.215698)));
    GK_ASSERT(as_float(dst[567221-start]) == as_float(half_t(5019.607910)));
    GK_ASSERT(as_float(dst[567226-start]) == as_float(half_t(10039.215820)));
    GK_ASSERT(as_float(dst[568398-start]) == as_float(half_t(20078.431641)));
    GK_ASSERT(as_float(dst[569882-start]) == as_float(half_t(max_val)));
    GK_ASSERT(as_float(dst[569883-start]) == as_float(half_t(max_val)));
    GK_ASSERT(as_float(dst[569886-start]) == as_float(half_t(max_val)));
}

void gtrack_builder_set_data_from_wig_test() {
    const float max_val = 40000.0f;
	const genome_t        genome{"hg19.mini"};
	genome_track::builder subject = build_gtrack_builder(max_val, genome);
	subject.set_data_from_wig(get_data_file_path("/tests/data/test_clamp/test.wig"));
    subject.finalize();

    genome_track track{};
    track.set_source("output.gtrack");
    track.open();
    /**
fixedStep chrom=chr1 start=1001 step=1 span=1
1000.0
fixedStep chrom=chr1 start=1002 step=1 span=1
10000.0
fixedStep chrom=chr1 start=1003 step=1 span=1
40000.0
fixedStep chrom=chr1 start=1004 step=1 span=1
70000.0

    Dict for these values:
        1000  => 941.176514
        10000 => 10039.215820
     */
    const int start = 1000;
    half_t dst[10] = {0};
	interval_t interval = interval_t::from_dna0(genome.chrom_names().as_chrom("chr1"), start, start + std::size(dst),
												strand_t::pos_strand, genome.refg());
	track(interval, dst, track.dtype());
	track.close();

    GK_ASSERT(as_float(dst[0]) == as_float(half_t(941.176514)));
    GK_ASSERT(as_float(dst[1]) == as_float(half_t(10039.215820)));
    GK_ASSERT(as_float(dst[2]) == as_float(half_t(max_val)));
    GK_ASSERT(as_float(dst[3]) == as_float(half_t(max_val)));
}

void gtrack_builder_set_data_from_bedgraph_test() {
    const float max_val = 40000.0f;
	const genome_t        genome{"hg19.mini"};
	genome_track::builder subject = build_gtrack_builder(max_val, genome);
	subject.set_data_from_bedgraph(get_data_file_path("/tests/data/test_clamp/test.bedgraph"));
    subject.finalize();

    genome_track track{};
    track.set_source("output.gtrack");
    track.open();
    /**
#bedGraph section chr1:1000-1008
chr1	1000	1002	1000
chr1	1002	1004	10000
chr1	1004	1006	40000
chr1	1006	1008	70000

    Dict for these values:
        1000  => 941.176514
        10000 => 10039.215820
     */
    const int start = 1000;
    half_t dst[10] = {0};
	interval_t interval = interval_t::from_dna0(genome.chrom_names().as_chrom("chr1"), start, start + std::size(dst),
												strand_t::pos_strand, genome.refg());
	track(interval, dst, track.dtype());
    track.close();

    GK_ASSERT(as_float(dst[0]) == as_float(half_t(941.176514)));
    GK_ASSERT(as_float(dst[2]) == as_float(half_t(10039.215820)));
    GK_ASSERT(as_float(dst[4]) == as_float(half_t(max_val)));
    GK_ASSERT(as_float(dst[6]) == as_float(half_t(max_val)));
}

void nested_exception_print(const std::exception& e, int level = 0)
{
	print("{:>{}}: {}\n", typeid(e).name(), strlen(typeid(e).name()) + level, e.what());
	try {
		std::rethrow_if_nested(e);
	} catch (const std::exception& e) {
		nested_exception_print(e, level + 2);
	}
}

int main(int argc, char* argv[])
{
	try {
		putenv("GENOMEKIT_QUIET=1");
		putenv("GENOMEKIT_DATA_DIR=tests/data/mini1");

		clamp_test();
		gtrack_builder_set_data_from_bed_test();
		gtrack_builder_set_data_from_wig_test();
		gtrack_builder_set_data_from_bedgraph_test();
		srand(18);
		init_genome_kit();
		walk_test();

		test_genome_track();

		/*
		resolve_datafile_path = custom_realize_datafile_path; // Allow C++-only mode resolve filenames similar to _gk_data_config

		//walk_test();

		vector<string> samfiles;
		//samfiles.push_back("SRR613342.splice.sorted.sam");
		//samfiles.push_back("SRR615105.splice.sorted.sam");
		//samfiles.push_back("SRR627419.splice.sorted.sam");

		//write_ralign_test(samfiles);
		//read_ralign_test(samfiles);
		//write_rdist_test(samfiles);
		//write_pooled_rdist_test(samfiles);
		//read_rdist_test(samfiles);
		//read_pooled_rdist_test();

		read_distributions::builder rdists("pooled.rdist");
		rdists.add("C:/data/rdist_bug/GTEx_demo.SRR612863.splice.sort.ralign");
		rdists.add("C:/data/rdist_bug/GTEx_demo.SRR659995.splice.sort.ralign");
		rdists.add("C:/data/rdist_bug/GTEx_demo.SRR654766.splice.sort.ralign");
		rdists.add("C:/data/rdist_bug/GTEx_demo.SRR607967.splice.sort.ralign");
		rdists.add("C:/data/rdist_bug/GTEx_demo.SRR612863.splice.sort.ralign");
		rdists.finalize();
		*/

	} catch (const std::exception& e) {
		nested_exception_print(e);
		return -1;
	}
	return 0;
}

#endif // _WANT_MAIN
