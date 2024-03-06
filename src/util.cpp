/*
Copyright (C) 2016-2023 Deep Genomics Inc. All Rights Reserved.
*/
#include "util.h"
#include "gk_assert.h"
#include <cstdarg>
#include <cstdio>
#include <ctime>

BEGIN_NAMESPACE_GK

void panic(const char* msg, ...)
{
	va_list va;
	va_start(va, msg);
	vfprintf(stderr, msg, va);
	va_end(va);
	fprintf(stderr, "\n");
	fflush(stderr);
	exit(-1);
}

// TODO: rename this function so it's clear it prints to stderr,
//       so that users know that the output printed here won't
//       be captured by standard stdout redirect.
void print(const char* msg, ...)
{
	va_list va;
	va_start(va, msg);
	vfprintf(stderr, msg, va);
	va_end(va);
	fflush(stderr);
}

void println(const char* msg, ...)
{
	va_list va;
	va_start(va, msg);
	vfprintf(stderr, msg, va);
	va_end(va);
	fprintf(stderr, "\n");
	fflush(stderr);
}

END_NAMESPACE_GK
