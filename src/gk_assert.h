/*
Copyright (C) 2016-2023 Deep Genomics Inc. All Rights Reserved.
*/
#pragma once
#ifndef __GK_ASSERT_H__
#define __GK_ASSERT_H__

#include "defines.h"
#include <format>
#include <stdexcept>
#include <string>

#define GK_ENABLE_DEBUGBREAK  // Uncomment this to enable debugbreak at outset of every GK_THROW or failed GK_ASSERT/GK_CHECK

BEGIN_NAMESPACE_GK

static const auto gk_debugbreak = getenv("GK_DEBUGBREAK") != nullptr;

class runtime_error : public std::runtime_error {
public:
	runtime_error(const char* msg, const char* file, int line): std::runtime_error(msg), file(file), line(line) {}
	runtime_error(const std::string& msg, const char* file, int line): std::runtime_error(msg), file(file), line(line) {}
	const char* what() const noexcept;

private:
	std::string buf;
	const char* file;
	int         line;
};

#define GK_DECL_ERROR_CLASS(error_class, base_class) \
	class error_class : public base_class { \
	public: \
		error_class(const char* msg, const char* file, int line): base_class(msg, file, line) { } \
		error_class(const std::string& msg, const char* file, int line): base_class(msg, file, line) { } \
	};

GK_DECL_ERROR_CLASS(assertion_error, runtime_error)
GK_DECL_ERROR_CLASS(file_error,      runtime_error)
GK_DECL_ERROR_CLASS(type_error,      runtime_error)
GK_DECL_ERROR_CLASS(value_error,     runtime_error)
GK_DECL_ERROR_CLASS(index_error,     runtime_error)
GK_DECL_ERROR_CLASS(key_error,       runtime_error)
GK_DECL_ERROR_CLASS(memory_error,    runtime_error)
GK_DECL_ERROR_CLASS(not_implemented_error,  runtime_error)
GK_DECL_ERROR_CLASS(unreachable_code_error, runtime_error)
// Special error case for when data manager fails to find a data file
GK_DECL_ERROR_CLASS(gk_data_file_not_found_error, runtime_error)

extern bool is_debugger_running();

#ifdef GK_ENABLE_DEBUGBREAK
// When a debugger is running, GK_THROW, GK_CHECK, and GK_ASSERT
	// will all trigger a debug break before throwing an exception,
	// giving the programmer a chance to inspect program state conveniently.
	#if defined(_MSC_VER)
		#define GK_DEBUGBREAK { if (gk_debugbreak && gk::is_debugger_running()) __debugbreak(); }
	#elif defined(__clang__)
		#define GK_DEBUGBREAK { if (gk_debugbreak && gk::is_debugger_running()) __builtin_debugtrap(); }
	#elif defined(__GNUC__)
		#define GK_DEBUGBREAK { if (gk_debugbreak && gk::is_debugger_running()) __builtin_trap(); }
	#else
		#error Unsupported compiler.
	#endif
#else
	#define GK_DEBUGBREAK // disabled
#endif

//! \brief Trigger a debug break on a specific condition.
//!
//! Useful for debugging optimized builds where conditional
//! breakpoints aren't reliable.
//!
#define GK_DEBUGBREAK_IF(cond) do { if (cond) { GK_DEBUGBREAK } } while (0)

//! \brief Throw an GK exception of the specified type.
//!
//! GK_THROW(etype, msg) throws an exception with a custom
//! string appended to the error message.
//!
//! GK_THROW(etype, msg, ...) throws an exception with a custom
//! formatted string appended to the error message. The format
//! for the message is the same as printf(msg, ...).
//!
#define GK_MAKE_ERROR(etype, msg, ...) \
	gk::etype##_error(std::format(msg __VA_OPT__(, ) __VA_ARGS__), __FILE__, __LINE__)
#define GK_THROW(etype, ...) throw GK_MAKE_ERROR(etype, __VA_ARGS__)
#define GK_CATCH_THROW_NESTED(etype, ...) \
	catch (const gk::etype##_error&) \
	{ \
		std::throw_with_nested(GK_MAKE_ERROR(etype, __VA_ARGS__)); \
	}
#define GK_DEBUGBREAK_THROW(etype, ...) \
	GK_DEBUGBREAK; \
	GK_THROW(etype, __VA_ARGS__)
#define GK_LIKELY_OR(cond, expr) \
	do { \
		if (LIKELY(cond)) { \
		} else { \
			expr; \
		} \
	} while (0)

// ordered from derived -> base
#define GK_RETHROW(...) \
	GK_CATCH_THROW_NESTED(assertion, __VA_ARGS__) \
	GK_CATCH_THROW_NESTED(file, __VA_ARGS__) \
	GK_CATCH_THROW_NESTED(type, __VA_ARGS__) \
	GK_CATCH_THROW_NESTED(value, __VA_ARGS__) \
	GK_CATCH_THROW_NESTED(index, __VA_ARGS__) \
	GK_CATCH_THROW_NESTED(key, __VA_ARGS__) \
	GK_CATCH_THROW_NESTED(memory, __VA_ARGS__) \
	GK_CATCH_THROW_NESTED(not_implemented, __VA_ARGS__) \
	GK_CATCH_THROW_NESTED(unreachable_code, __VA_ARGS__) \
	GK_CATCH_THROW_NESTED(gk_data_file_not_found, __VA_ARGS__) \
	GK_CATCH_THROW_NESTED(runtime, __VA_ARGS__)

//! \brief If expr is false, throw an assertion_error.
//!
//! GK_ASSERT(expr) throws an assertion_error if expr evaluates to
//! false. The failed expression is appended to the error message.
//!
//! GK_ASSERT(expr, msg) throws if expr fails, but appends a
//! custom error message that is presumably more informative than
//! just the failed expression.
//!
//! GK_ASSERT(expr, msg, ...) throws if expr fails, but appends
//! a custom formatted string appended to the error message. The
//! format for the message is the same as printf(msg, ...).
//!
#define GK_ASSERT(expr, ...) GK_LIKELY_OR(expr, GK_DEBUGBREAK_THROW(assertion, "({}): " GK_VA_HEAD(__VA_ARGS__), #expr GK_VA_COMMA_TAIL(__VA_ARGS__)))

//! \brief If expr is false, throw a specific type of exception.
//!
//! GK_CHECK(expr, etype, msg) throws if expr fails, but appends a
//! custom error message that is presumably more informative than
//! just the failed expression.
//!
//! GK_CHECK(expr, etype, msg, ...) throws if expr fails, but appends
//! a custom formatted string appended to the error message. The
//! format for the message is the same as printf(msg, ...).
//!
#define GK_CHECK(expr, etype, ...)  GK_LIKELY_OR(expr, GK_DEBUGBREAK_THROW(etype, __VA_ARGS__))

//! \brief Throw an unreachable_code_error exception.
#define GK_UNREACHABLE()            do { GK_DEBUGBREAK_THROW(unreachable_code, ""); } while (0)

//! \brief Throw a not_implemented_error exception.
#define GK_NOT_IMPLEMENTED()        do { GK_DEBUGBREAK_THROW(not_implemented, ""); } while (0)

#ifdef GK_DEBUG
#ifndef GK_ENABLE_DBASSERT
#define GK_ENABLE_DBASSERT
#endif
#endif

#ifdef GK_ENABLE_DBASSERT
#define GK_DBASSERT(...)      GK_ASSERT(__VA_ARGS__)
#else
#define GK_DBASSERT(...)      { }
#endif

END_NAMESPACE_GK

#endif  // __GK_ASSERT_H__
