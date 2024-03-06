/*
Copyright (C) 2016-2023 Deep Genomics Inc. All Rights Reserved.
*/
#include "genome_kit.h"

BEGIN_NAMESPACE_GK

void init_genome_kit();

// Explicitly call inplace constructors of global objects if loaded by the C runtime,
// since the C runtime will not trigger construction of C++ globals. Nothing is
// done for destructors in this case, so don't make any globals that need to be
// notified of program exit for correct functionality.
static void init_globals_cruntime()
{
}

//////////////////////////////////////////////////////

void init_dna_tables(); // in dnaseq.cpp

void init_genome_kit()
{
	init_globals_cruntime();
	init_dna_tables();
}

END_NAMESPACE_GK
