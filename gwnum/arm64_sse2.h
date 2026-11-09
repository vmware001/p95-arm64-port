/* arm64_sse2.h - the tiny subset of "SSE2 double" intrinsics that polymult.c
 * uses, implemented with GCC/Clang 128-bit vector extensions on ARM64.
 *
 * Why this exists: mprime's P-1 and ECM stage 2 call polymult, and gwnum's FFT
 * data layout on ARM64 is the 128-bit (two doubles per vector) layout - exactly
 * what upstream's SSE2 polymult variant expects.  The portable scalar ("dbl")
 * polymult variant uses the 32-bit/x87 one-double layout upstream, which does
 * not match gwnum's ARM64 layout (that mismatch is what made P-1/ECM stage 2
 * crash on the ARM64 port).  Building polymult.c with -DSSE2 against this
 * header gives a correctly-laid-out 128-bit polymult using NEON code that the
 * compiler generates from the plain vector arithmetic below.
 *
 * Only the intrinsics actually referenced by polymult.c are provided:
 *   _mm_set_pd, _mm_set1_pd, _mm_add_pd, _mm_sub_pd, _mm_mul_pd, _mm_div_pd,
 *   _mm_stream_pd
 * This is NOT a general SSE/SSE2 emulation header.
 */
#ifndef GW_ARM64_SSE2_H
#define GW_ARM64_SSE2_H

/* GCC/Clang vector type: two doubles, 16-byte aligned, NEON-mapped. */
typedef double __m128d __attribute__ ((__vector_size__ (16), __aligned__ (16)));

static __inline__ __m128d _mm_set_pd (double hi, double lo)
{
	__m128d r;
	r[0] = lo;
	r[1] = hi;
	return r;
}

static __inline__ __m128d _mm_set1_pd (double a)
{
	__m128d r;
	r[0] = a;
	r[1] = a;
	return r;
}

static __inline__ __m128d _mm_add_pd (__m128d a, __m128d b) { return a + b; }
static __inline__ __m128d _mm_sub_pd (__m128d a, __m128d b) { return a - b; }
static __inline__ __m128d _mm_mul_pd (__m128d a, __m128d b) { return a * b; }
static __inline__ __m128d _mm_div_pd (__m128d a, __m128d b) { return a / b; }

/* Non-temporal store: use a plain (cacheable) 16-byte store - ARM has no
 * direct non-temporal store for this, and correctness/layout is what matters. */
static __inline__ void _mm_stream_pd (double *p, __m128d a)
{
	*(__m128d __attribute__ ((__may_alias__)) *) p = a;
}

#endif /* GW_ARM64_SSE2_H */
