/*
Copyright (C) 2016-2023 Deep Genomics Inc. All Rights Reserved.
*/
#ifndef __GENOME_KIT_TIME_H__
#define __GENOME_KIT_TIME_H__

#include "defines.h"
#ifndef _WIN32
#include <ctime>
#endif

BEGIN_NAMESPACE_GK

#ifdef _WIN32
using ticks_t = long long;
#else
using ticks_t = std::clock_t;
#endif
ticks_t ticks();
ticks_t tickrate(); // ticks per second
double  duration(ticks_t ticks);

///////////////////////////////////////////

class timer_t {
public:
	INLINE void tic()  { _tic = ticks(); }
	INLINE void toc()  { _toc = ticks(); }
	INLINE void reset(){ _tic = _toc = 0; }
	INLINE void rtic() { reset(); tic(); }
	INLINE double duration() const { return gk::duration(_toc-_tic); }

private:
	ticks_t _tic{};
	ticks_t _toc{};
};

extern timer_t g_timer;

END_NAMESPACE_GK

#endif // __GENOME_KIT_TIME_H__
