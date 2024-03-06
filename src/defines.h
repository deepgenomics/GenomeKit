/*
Copyright (C) 2016-2023 Deep Genomics Inc. All Rights Reserved.
*/
#pragma once
#ifndef __DEFINES_H__
#define __DEFINES_H__

#if defined(__APPLE__)
#include <machine/endian.h>
#elif defined(__GLIBC__)
#include <endian.h>
#endif
#if defined(__BYTE_ORDER) && __BYTE_ORDER == __BIG_ENDIAN
#define _BIGENDIAN
#endif

#ifndef NDEBUG
#ifndef GK_DEBUG
#define GK_DEBUG
#endif
#endif

#if (defined _MSC_VER) || (defined __CYGWIN__)
// MSVC-specific macros.
#define INLINE __forceinline
#define RESTRICT __restrict
#define LIKELY(exp) exp
#elif (defined __GNUC__)
// GCC-specific macros. Used when compiled on Linux or on Windows via MinGW (as invoked by Anaconda setup.py)
#define INLINE inline __attribute__((always_inline))
#define RESTRICT __restrict__
#define LIKELY(exp) __builtin_expect(!!(exp), 1)
#else
#error Unsupported compiler.
#endif

#include <cstdint>

#define NOCOPY(C) C(const C&) = delete; C& operator=(const C&) = delete;

#define BEGIN_NAMESPACE_GK namespace gk {
#define END_NAMESPACE_GK   }
#define USING_NAMESPACE_GK using namespace gk;

#define GK_EXPAND_EMPTY
#define GK_EXPAND(x) x
#define GK_EXPAND_STR(x) #x

#define GK_VA_HEAD(head, ...)  head
#define GK_VA_COMMA_TAIL(head, ...)  __VA_OPT__(,) __VA_ARGS__

// Shorthand so that casting is not so verbose
#define rcast reinterpret_cast
#define scast static_cast
#define ccast const_cast

#endif // __DEFINES_H__
