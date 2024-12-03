/*
Copyright (C) 2016-2023 Deep Genomics Inc. All Rights Reserved.
*/
#pragma once
#ifndef __GENOME_KIT_PY_GENOME_TRACK_H__
#define __GENOME_KIT_PY_GENOME_TRACK_H__

#include "py_util.h"
#include "genome_track.h"

BEGIN_NAMESPACE_GK

GKPY_TYPE_BEGIN(GenomeTrack)
	genome_track* track;        // Pointer to the track object we're wrapping
GKPY_TYPE_END

GKPY_TYPE_BEGIN(GenomeTrackBuilder)
	genome_track::builder* builder{};
	PyObject*              genome{};
GKPY_TYPE_END

END_NAMESPACE_GK

#endif // __GENOME_KIT_PY_GENOME_TRACK_H__
