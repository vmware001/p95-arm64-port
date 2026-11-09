#pragma once

// prefetch disabled
#define pfing(...)

#define exec(...) __VA_ARGS__
#define noexec(...)
#define nonoexec(...) __VA_ARGS__

/// __VA_ARGS__

#define clear_timers
#define start_timer(...)
#define end_timer(...)

#define loops_reset
#define loops_init(...)
#define loops_init_prefetch(...)
#define loops_undo(...)
#define loops_init_prefetch_touch(...)

#define loops_set_prefetch(A,B,C)

#define likely(x) (__builtin_expect((x)!=0, 1))
#define unlikely(x) (__builtin_expect(x, 0))

#if defined(__clang__)
#define UNROLL(N)	_Pragma("clang loop0 unroll_count(N)")
#define UNROLL_FULL	_Pragma("clang loop unroll(full)")
#define NO_UNROLL	_Pragma("clang loop unroll(disable)")

#elif defined(__GNUC__) && (__GNUC__ >= 8)
#define UNROLL(N)	_Pragma("GCC unroll N")
#define UNROLL_FULL	_Pragma("GCC unroll 64")
#define NO_UNROLL	_Pragma("GCC unroll 0")

#else
#define UNROLL(N)
#define UNROLL_FULL
#define NO_UNROLL
#endif
