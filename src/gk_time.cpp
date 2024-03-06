/*
Copyright (C) 2016-2023 Deep Genomics Inc. All Rights Reserved.
*/
#include "gk_time.h"

#ifdef _WIN32
#define WIN32_LEAN_AND_MEAN
#define NOMINMAX
#include <windows.h>
#endif

BEGIN_NAMESPACE_GK

#ifdef _WIN32

static ticks_t g_tickrate = 0;

ticks_t ticks() // this is not thread safe
{
	ticks_t result = 0;
	QueryPerformanceCounter((LARGE_INTEGER*)&result);
	return result;
}

ticks_t tickrate()
{
	if (g_tickrate == 0)
		QueryPerformanceFrequency((LARGE_INTEGER*)&g_tickrate);
	return g_tickrate;
}

#else

ticks_t ticks()    { return ::clock(); }
ticks_t tickrate() { return CLOCKS_PER_SEC; }

#endif

double duration(ticks_t ticks)
{
	return (double)ticks / tickrate();
}

timer_t g_timer; // convenient, but obviously not safe to use concurrently from multiple threads

END_NAMESPACE_GK
