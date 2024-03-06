/*
Copyright (C) 2016-2023 Deep Genomics Inc. All Rights Reserved.
*/
#pragma once
#ifndef __GENOME_KIT_API_H__
#define __GENOME_KIT_API_H__

#include "genome.h"
#include "biostr.h"
#include "interval.h"
#include "genome_dna.h"
#include "genome_anno.h"
#include "genome_track.h"
#include <vector>

BEGIN_NAMESPACE_GK

// Call one of these functions when the module is loaded.
void init_genome_kit();

END_NAMESPACE_GK

#endif // __GENOME_KIT_API_H__
