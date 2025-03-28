/*
Copyright (C) 2016-2023 Deep Genomics Inc. All Rights Reserved.
*/
#pragma once
#ifndef __GENOME_KIT_UTIL_H__
#define __GENOME_KIT_UTIL_H__

#include "defines.h"
#include "gk_assert.h"
#include <algorithm>
#include <bit>
#include <cstdlib>
#include <cstring>
#include <format>
#include <iostream>
#include <functional>
#include <new>
#include <type_traits>
#include <utility>

BEGIN_NAMESPACE_GK

template<typename ...T>
void print(std::format_string<T...> fmt, T&&... args)
{
	std::cerr << std::format(fmt, std::forward<T>(args)...);
}

template <typename... T>
void println(std::format_string<T...> fmt, T&&... args)
{
	std::cerr << std::format(fmt, std::forward<T>(args)...) << '\n';
}

///////////////////////////////////////////////////

INLINE bool _int_cast_check(unsigned char     , unsigned char      x) { return true;                 }
INLINE bool _int_cast_check(unsigned char     ,          char      x) { return x >= 0;               }
INLINE bool _int_cast_check(unsigned char     , unsigned short     x) { return           x <= 0xffu; }
INLINE bool _int_cast_check(unsigned char     ,          short     x) { return x >= 0 && x <= 0xff;  }
INLINE bool _int_cast_check(unsigned char     , unsigned int       x) { return           x <= 0xffu; }
INLINE bool _int_cast_check(unsigned char     ,          int       x) { return x >= 0 && x <= 0xff;  }
INLINE bool _int_cast_check(unsigned char     , unsigned long int  x) { return           x <= 0xffu; }
INLINE bool _int_cast_check(unsigned char     ,          long int  x) { return x >= 0 && x <= 0xff;  }
INLINE bool _int_cast_check(unsigned char     , unsigned long long x) { return           x <= 0xffu; }
INLINE bool _int_cast_check(unsigned char     ,          long long x) { return x >= 0 && x <= 0xff;  }

INLINE bool _int_cast_check(         char     , unsigned char      x) { return               x <= 0x7f;    }
INLINE bool _int_cast_check(         char     ,          char      x) { return true;                       }
INLINE bool _int_cast_check(         char     , unsigned short     x) { return               x <= 0x7fu;   }
INLINE bool _int_cast_check(         char     ,          short     x) { return x >= -0x80 && x <= 0x7f;    }
INLINE bool _int_cast_check(         char     , unsigned int       x) { return               x <= 0x7fu;   }
INLINE bool _int_cast_check(         char     ,          int       x) { return x >= -0x80 && x <= 0x7f;    }
INLINE bool _int_cast_check(         char     , unsigned long int  x) { return               x <= 0x7full; }
INLINE bool _int_cast_check(         char     ,          long int  x) { return x >= -0x80 && x <= 0x7fll;  }
INLINE bool _int_cast_check(         char     , unsigned long long x) { return               x <= 0x7full; }
INLINE bool _int_cast_check(         char     ,          long long x) { return x >= -0x80 && x <= 0x7fll;  }

INLINE bool _int_cast_check(unsigned short    , unsigned char      x) { return true;                     }
INLINE bool _int_cast_check(unsigned short    ,          char      x) { return x >= 0;                   }
INLINE bool _int_cast_check(unsigned short    , unsigned short     x) { return true;                     }
INLINE bool _int_cast_check(unsigned short    ,          short     x) { return x >= 0;                   }
INLINE bool _int_cast_check(unsigned short    , unsigned int       x) { return           x <= 0xffffu;   }
INLINE bool _int_cast_check(unsigned short    ,          int       x) { return x >= 0 && x <= 0xffff;    }
INLINE bool _int_cast_check(unsigned short    , unsigned long int  x) { return           x <= 0xffffull; }
INLINE bool _int_cast_check(unsigned short    ,          long int  x) { return x >= 0 && x <= 0xffffll;  }
INLINE bool _int_cast_check(unsigned short    , unsigned long long x) { return           x <= 0xffffull; }
INLINE bool _int_cast_check(unsigned short    ,          long long x) { return x >= 0 && x <= 0xffffll;  }

INLINE bool _int_cast_check(         short    , unsigned char      x) { return true;                           }
INLINE bool _int_cast_check(         short    ,          char      x) { return true;                           }
INLINE bool _int_cast_check(         short    , unsigned short     x) { return                 x <= 0x7fffu;   }
INLINE bool _int_cast_check(         short    ,          short     x) { return true;                           }
INLINE bool _int_cast_check(         short    , unsigned int       x) { return                 x <= 0x7fffu;   }
INLINE bool _int_cast_check(         short    ,          int       x) { return x >= -0x8000 && x <= 0x7fff;    }
INLINE bool _int_cast_check(         short    , unsigned long int  x) { return                 x <= 0x7fffull; }
INLINE bool _int_cast_check(         short    ,          long int  x) { return x >= -0x8000 && x <= 0x7fffll;  }
INLINE bool _int_cast_check(         short    , unsigned long long x) { return                 x <= 0x7fffull; }
INLINE bool _int_cast_check(         short    ,          long long x) { return x >= -0x8000 && x <= 0x7fffll;  }

INLINE bool _int_cast_check(unsigned int      , unsigned char      x) { return true;                         }
INLINE bool _int_cast_check(unsigned int      ,          char      x) { return x >= 0;                       }
INLINE bool _int_cast_check(unsigned int      , unsigned short     x) { return true;                         }
INLINE bool _int_cast_check(unsigned int      ,          short     x) { return x >= 0;                       }
INLINE bool _int_cast_check(unsigned int      , unsigned int       x) { return true;                         }
INLINE bool _int_cast_check(unsigned int      ,          int       x) { return x >= 0;                       }
INLINE bool _int_cast_check(unsigned int      , unsigned long int  x) { return           x <= 0xffffffffull; }
INLINE bool _int_cast_check(unsigned int      ,          long int  x) { return x >= 0 && x <= 0xffffffffll;  }
INLINE bool _int_cast_check(unsigned int      , unsigned long long x) { return           x <= 0xffffffffull; }
INLINE bool _int_cast_check(unsigned int      ,          long long x) { return x >= 0 && x <= 0xffffffffll;  }

INLINE bool _int_cast_check(         int      , unsigned char      x) { return true;                                     }
INLINE bool _int_cast_check(         int      ,          char      x) { return true;                                     }
INLINE bool _int_cast_check(         int      , unsigned short     x) { return true;                                     }
INLINE bool _int_cast_check(         int      ,          short     x) { return true;                                     }
INLINE bool _int_cast_check(         int      , unsigned int       x) { return                       x <= 0x7fffffffu;   }
INLINE bool _int_cast_check(         int      ,          int       x) { return true;                                     }
INLINE bool _int_cast_check(         int      , unsigned long int  x) { return                       x <= 0x7fffffffull; }
INLINE bool _int_cast_check(         int      ,          long int  x) { return x >= -0x80000000ll && x <= 0x7fffffffll;  }
INLINE bool _int_cast_check(         int      , unsigned long long x) { return                       x <= 0x7fffffffull; }
INLINE bool _int_cast_check(         int      ,          long long x) { return x >= -0x80000000ll && x <= 0x7fffffffll;  }

INLINE bool _int_cast_check(unsigned long int , unsigned char      x) { return true;   }
INLINE bool _int_cast_check(unsigned long int ,          char      x) { return x >= 0; }
INLINE bool _int_cast_check(unsigned long int , unsigned short     x) { return true;   }
INLINE bool _int_cast_check(unsigned long int ,          short     x) { return x >= 0; }
INLINE bool _int_cast_check(unsigned long int , unsigned int       x) { return true;   }
INLINE bool _int_cast_check(unsigned long int ,          int       x) { return x >= 0; }
INLINE bool _int_cast_check(unsigned long int , unsigned long int  x) { return true;   }
INLINE bool _int_cast_check(unsigned long int ,          long int  x) { return x >= 0; }
INLINE bool _int_cast_check(unsigned long int , unsigned long long x) { return true;   }
INLINE bool _int_cast_check(unsigned long int ,          long long x) { return x >= 0; }

INLINE bool _int_cast_check(         long int , unsigned char      x) { return true;               }
INLINE bool _int_cast_check(         long int ,          char      x) { return true;               }
INLINE bool _int_cast_check(         long int , unsigned short     x) { return true;               }
INLINE bool _int_cast_check(         long int ,          short     x) { return true;               }
INLINE bool _int_cast_check(         long int , unsigned int       x) { return true;               }
INLINE bool _int_cast_check(         long int ,          int       x) { return true;               }
INLINE bool _int_cast_check(         long int , unsigned long int  x) { return x <= 0x7fffffffffffffffull; }
INLINE bool _int_cast_check(         long int ,          long int  x) { return true;               }
INLINE bool _int_cast_check(         long int , unsigned long long x) { return x <= 0x7fffffffffffffffull; }
INLINE bool _int_cast_check(         long int ,          long long x) { return true;               }

INLINE bool _int_cast_check(unsigned long long, unsigned char      x) { return true;   }
INLINE bool _int_cast_check(unsigned long long,          char      x) { return x >= 0; }
INLINE bool _int_cast_check(unsigned long long, unsigned short     x) { return true;   }
INLINE bool _int_cast_check(unsigned long long,          short     x) { return x >= 0; }
INLINE bool _int_cast_check(unsigned long long, unsigned int       x) { return true;   }
INLINE bool _int_cast_check(unsigned long long,          int       x) { return x >= 0; }
INLINE bool _int_cast_check(unsigned long long, unsigned long int  x) { return true;   }
INLINE bool _int_cast_check(unsigned long long,          long int  x) { return x >= 0; }
INLINE bool _int_cast_check(unsigned long long, unsigned long long x) { return true;   }
INLINE bool _int_cast_check(unsigned long long,          long long x) { return x >= 0; }

INLINE bool _int_cast_check(         long long, unsigned char      x) { return true;                       }
INLINE bool _int_cast_check(         long long,          char      x) { return true;                       }
INLINE bool _int_cast_check(         long long, unsigned short     x) { return true;                       }
INLINE bool _int_cast_check(         long long,          short     x) { return true;                       }
INLINE bool _int_cast_check(         long long, unsigned int       x) { return true;                       }
INLINE bool _int_cast_check(         long long,          int       x) { return true;                       }
INLINE bool _int_cast_check(         long long, unsigned long int  x) { return x <= 0x7fffffffffffffffull; }
INLINE bool _int_cast_check(         long long,          long int  x) { return true;                       }
INLINE bool _int_cast_check(         long long, unsigned long long x) { return x <= 0x7fffffffffffffffull; }
INLINE bool _int_cast_check(         long long,          long long x) { return true;                       }

template <typename Y, typename X>
INLINE Y int_cast(X x)
{
	Y y = Y(x);
	GK_CHECK(_int_cast_check(y, x), value, "int_cast: integer overflow when casting {}.", x);
	return y;
}

template <typename Y, typename X>
INLINE Y float_cast(X x)
{
	Y y = Y(x);
	GK_CHECK(x == X(y), value, "float_cast: loss of precision when casting {}.", x);
	return y;
}

template <typename T, typename V>
constexpr bool in_range(T x, V lo, V hi)
{
	if constexpr (std::is_unsigned_v<T>)
		return (lo == 0 || lo <= hi) && x <= hi;
	return lo <= x && x <= hi;
}

#define GK_DEFINE_INT_TRAITS(T, _min, _max, D, P, U) \
	template <> struct int_traits<T> { \
		static const T min = _min; \
		static const T max = _max; \
		using demoted_type = D; \
		using promoted_type = P; \
		using unsigned_type = U; \
	}

template <typename T> struct int_traits { };
GK_DEFINE_INT_TRAITS(unsigned char     ,                            0, 255u,                    void,           unsigned short    , unsigned char );
GK_DEFINE_INT_TRAITS(         char     ,                         -128, 127,                     void,                    short    , unsigned char );
GK_DEFINE_INT_TRAITS(unsigned short    ,                            0, 65535u,                  unsigned char,  unsigned int      , unsigned short);
GK_DEFINE_INT_TRAITS(         short    ,                       -32768, 32767,                            char,           int      , unsigned short);
GK_DEFINE_INT_TRAITS(unsigned int      ,                            0, 4294967295u,             unsigned short, unsigned long long, unsigned int  );
GK_DEFINE_INT_TRAITS(         int      ,                -2147483648ll, 2147483647,                       short,          long long, unsigned int  );
GK_DEFINE_INT_TRAITS(unsigned long     ,                            0, 18446744073709551615ull, unsigned int,   void              , unsigned long );
GK_DEFINE_INT_TRAITS(         long     ,       (long)0x800000000000ul, 9223372036854775807l,             int,   void              , unsigned long );
GK_DEFINE_INT_TRAITS(unsigned long long,                            0, 18446744073709551615ull, unsigned int,   void              , unsigned long long);
GK_DEFINE_INT_TRAITS(         long long, (long long)0x800000000000ull, 9223372036854775807ll,            int,   void              , unsigned long long);

template <class AlignedT, class T>
T aligned_distance(T x) noexcept
{
	static constexpr auto alignment = sizeof(AlignedT);
	return scast<T>((alignment - x % alignment) % alignment);
}

template <class T>
constexpr T high_bit() noexcept
{
	return std::rotr(T{1}, 1);
}

template <class Enum>
constexpr std::underlying_type_t<Enum> as_ordinal(Enum e) noexcept
{
	return static_cast<std::underlying_type_t<Enum>>(e);
}

template <class T>
T* as_mutable(const T* p) noexcept
{
	return const_cast<T*>(p);
}

template <class T>
T& as_mutable(const T& p) noexcept
{
	return const_cast<T&>(p);
}

////////////////////////////////////////////////////

template <typename T>
struct range_t {
	using cursor_t = T;

	T a;
	T b;

	INLINE       T& begin()       { return a; }
	INLINE const T& begin() const { return a; }
	INLINE       T& end()         { return b; }
	INLINE const T& end()   const { return b; }
	INLINE size_t   size()  const { return (size_t)(b - a); }
	INLINE bool     empty() const { return a == b; }
};

/////////////////////////////////////////////////////

template <typename T, typename P> void sort(T& c, const P& pred)        { std::sort(c.begin(), c.end(), pred); }
template <typename T, typename P> void sort(T& c)                       { std::sort(c.begin(), c.end(), std::less<typename T::value_type>()); }
template <typename T, typename P> void stable_sort(T& c, const P& pred) { std::stable_sort(c.begin(), c.end(), pred); }
template <typename T, typename P> void stable_sort(T& c)                { std::stable_sort(c.begin(), c.end(), std::less<typename T::value_type>()); }

template <class C, class K>
const typename C::mapped_type& find_or(const C& container, const K& key, const typename C::mapped_type& default_value)
{
	if (auto it = container.find(key); it != std::cend(container)) {
		return it->second;
	}
	return default_value;
}

////////////////////////////////////////////////////

struct uuid_t {
	uint64_t low() const noexcept {
		uint64_t ret;
		static_assert(sizeof(ret) * 2 == sizeof(bytes));
		std::memcpy(&ret, &bytes[8], sizeof(ret));
		return ret;
	}
	uint64_t high() const noexcept {
		uint64_t ret;
		static_assert(sizeof(ret) * 2 == sizeof(bytes));
		std::memcpy(&ret, &bytes[0], sizeof(ret));
		return ret;
	}

	uint8_t  bytes[16];
};
inline bool operator==(uuid_t lhs, uuid_t rhs)
{
	return std::memcmp(lhs.bytes, rhs.bytes, sizeof(lhs.bytes)) == 0;
}

////////////////////////////////////////////////////

template <typename T> void align_ptr(T*& ptr, int align)
{
	auto p = (size_t)ptr;
	size_t res = p % align;
	if (res)
		p += align - res;
	ptr = (T*)p;
}

////////////////////////////////////////////////////

template <typename T> INLINE void destruct(T* x) { x->~T(); }
template <typename T> INLINE void default_construct(T* x) { new (x) T; }

// Automatically calls destructor (but not delete!) on an object if it was not told to release the object.
template <typename T>
class constructor_tag { NOCOPY(constructor_tag)
public:
	INLINE constructor_tag(T* ptr) noexcept : _ptr(ptr) { }
	INLINE ~constructor_tag() { if (_ptr) destruct(_ptr); }
	INLINE void release() noexcept { _ptr = nullptr; }

private:
	T* _ptr;
};

#define GK_TENTATIVE_INPLACE_CONSTRUCT(member_name, type, ...) \
	new (&(self->member_name)) type(__VA_ARGS__);\
	constructor_tag<type> member_name##_tag(&self->member_name)

#define GK_FINALIZE_CONSTRUCT(tag) tag##_tag.release()

////////////////////////////////////////////////////

#ifdef _MSC_VER
#define bswap32(x) _byteswap_ulong(x)
#define bswap64(x) _byteswap_uint64(x)
#else
#define bswap32(x) __builtin_bswap32(x)
#define bswap64(x) __builtin_bswap64(x)
#endif

////////////////////////////////////////////////////

template <typename T, typename U> INLINE T divup(T x, U denom) { return (x+denom-1) / denom;     }
template <typename T, typename U> INLINE T rndup(T x, U align) { return divup(x, align) * align; }
template <typename T, typename U> INLINE T rnddn(T x, U align) { return (x/align) * align;       }

// The udiv* functions are useful for dividing signed integers with
// unsigned division. Unsigned division is faster when both 'a' and 'b'
// are known to be non-negative and at least one of them is known at compile
// time (usually the denominator).
//
// Constant propagation allows the compiler to emit simpler arithmetic
// instructions, often a single shift instruction, instead of division.
//
// One might say "why not just use unsigned types for values that
// are assumed to be negative?" There are very good reasons not to use
// unsigned types by default, even for values that "should" never be negative,
// because it's very tricky to detect underflow for those variables (input validation),
// and bugs (like accidental infinite loops) creep in too easily.
//
// See discussion https://google.github.io/styleguide/cppguide.html#Integer_Types
//
template <typename T, typename V> INLINE T udivdn(T x, V denom)
{
	using U = typename int_traits<T>::unsigned_type;
	return T((U)x / (U)denom);  // Round down
}

template <typename T, typename V> INLINE T udivup(T x, V denom)
{
	using U = typename int_traits<T>::unsigned_type;
	return T(((U)x+(U)denom-1) / (U)denom);  // Round up
}

// Unsigned modulus operator (%) applicable to signed types.
template <typename T, typename V> INLINE T umod(T x, V denom)
{
	using U = typename int_traits<T>::unsigned_type;
	return T((U)x % (U)denom);  // Round down
}

////////////////////////////////////////////////////

#define GK_DECLARE_ARRAY_TYPES(_value_type, _index_type) \
	using value_type = _value_type; \
	using index_type = _index_type;

#define GK_DECLARE_ARRAY_CONST_METHODS(data_expr, size_expr) \
	INLINE index_type size()  const { return size_expr; } \
	INLINE bool       empty() const { return size()==0; } \
	INLINE const value_type& operator[](index_type i) const { GK_DBASSERT(i < size()); return (data_expr)[i]; } \

#define GK_DECLARE_ARRAY_METHODS(data_expr, size_expr) \
	GK_DECLARE_ARRAY_CONST_METHODS(data_expr, size_expr) \
	INLINE value_type& operator[](index_type i) { GK_DBASSERT(i < size()); return (data_expr)[i]; } \

#define GK_DECLARE_ARRAY_CONST_ITERATOR(first) \
	using const_iterator = const value_type*; \
	INLINE const_iterator begin() const { return first; } \
	INLINE const_iterator end()   const { return first + size(); }

#define GK_DECLARE_ARRAY_ITERATOR(first) \
	GK_DECLARE_ARRAY_CONST_ITERATOR(first) \
	using iterator = value_type*; \
	INLINE iterator begin() { return first; } \
	INLINE iterator end()   { return first + size(); }

/////////////////////////////////////////////////////

#pragma pack(push, 1)
// A 40-bit offset type, for shaving 24 bits off of unnecessarily-large 64-bit offsets.
// Useful when trying to squeeze the most into memory mapped binary file formats.
struct offset40_t {
	uint32_t lo;
	uint8_t  hi;

	offset40_t() = default;
	INLINE offset40_t(uint64_t offset) { *this = offset; }

	INLINE offset40_t& operator=(uint64_t offset)
	{
		lo = (uint32_t)offset;
		hi = (uint8_t)(offset >> 32);
		GK_CHECK((offset >> 32) <= 255, value, "Overflow when truncating offset to 40 bits");
		return *this;
	}

	INLINE uint64_t as64() const
	{
		return ((uint64_t)hi << 32) + (uint64_t)lo;
	}
};
#pragma pack(pop)

END_NAMESPACE_GK

#endif // __GENOME_KIT_UTIL_H__
