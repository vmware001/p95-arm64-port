#pragma once
/* Copyright 2001-2024 - Mersenne Research, Inc.  All rights reserved */
/* Author:  George Woltman */
/* Email: woltman@alum.mit.edu */

/* These macros implement basic SSE2 building blocks that will be used by */
/* all FFT types.  This file contains the default definitions which */
/* are optimized for the Core 2 / Core i7 architecture.  Optimized versions */
/* of these macros for other architectures can be found in other files. */

/* Converted from from assembly to C code, reduced the size of this file. */


// TODO: use __builtin_prefetch(a,is_write, level);
// level 0: stream, 1: L3 cache, 2: L2 cache, 3: L1 cache

#define xtouch(A)
#define xtouch_init(A)
#define xprefetch(A)
#define xxprefetch(A)
#define xprefetchw(A)
#define L2prefetch128(A)


#define xs_complex_square xp_complex_square
#if 1
#define xs_complex_mult xp_complex_mult
#define xp_complex_mult(real1, imag1, real2, imag2, tmp1, tmp2) \
	tmp1 = real1 * real2; \
	tmp2 = real1 * imag2 + real2 * imag1; \
	real1 = tmp1 - imag1 * imag2; \
	imag1 = tmp2;
#else
#define xs_complex_mult(real1, imag1, real2m, imag2m, tmp1, tmp2) { \
	double real2 = real2m, imag2 = imag2m; \
	tmp1 = real1; \
	real1 *= real2;		/* real1 * real2 */ \
	tmp2 = imag1 * imag2; /* imag1 * imag2 */ \
	tmp1 *= imag2; /* real1 * imag2 */ \
	imag1 *= real2;		/* real2 * imag1 */ \
	real1 -= tmp2;		/* real1*real2-imag1*imag2 (new real) */ \
	imag1 += tmp1;		/* real1*imag2+real2*imag1 (new imag) */ \
}
#define xp_complex_mult(real1, imag1, real2m, imag2m, tmp1, tmp2) { \
	vec2f64 real2 = real2m, imag2 = imag2m; \
	tmp1 = real1; \
	real1 *= real2;		/* real1 * real2 */ \
	tmp2 = imag1 * imag2; /* imag1 * imag2 */ \
	tmp1 *= imag2; /* real1 * imag2 */ \
	imag1 *= real2;		/* real2 * imag1 */ \
	real1 -= tmp2;		/* real1*real2-imag1*imag2 (new real) */ \
	imag1 += tmp1;		/* real1*imag2+real2*imag1 (new imag) */ \
}
#endif
#if 1
#define xp_complex_square(real, imag, tmp) \
	tmp = imag * real; \
	real = real * real; \
	real -= imag * imag; \
	imag = tmp + tmp;
#else
#define xp_complex_square(real, imag, tmp) \
	tmp = imag * imag; \
	imag *= real; \
	real *= real; \
	real -= tmp;		/* real^2 - imag^2 (new real) */ \
	imag += imag;		/* imag * real * 2 (new imag) */
#endif


// xs4c_mulf only used in hg.h, s2cl_eight_reals_with_mult_2
#if 0
#define xs4c_mulf(r1,r2,r3,r4,r5,r6,r7,r8,m1,m2,m3,m4,m5,m6,m7,m8) { \
	vec2f64 tmp8, tmp9; \
	vec2f64 r13 = {r1, r3}; \
	vec2f64 r24 = {r2, r4}; \
	vec2f64 r57 = {r5, r7}; \
	vec2f64 r68 = {r6, r8}; \
	vec2f64 m13 = {f64ptr(m1+rbp), f64ptr(m3+rbp)}; \
	vec2f64 m24 = {f64ptr(m2+rbp), f64ptr(m4+rbp)}; \
	xp_complex_mult(r13, r24, m13, m24, tmp8, tmp9); \
	vec2f64 m57 = {f64ptr(m5+rbp), f64ptr(m7+rbp)}; \
	vec2f64 m68 = {f64ptr(m6+rbp), f64ptr(m8+rbp)}; \
	xp_complex_mult(r57, r68, m57, m68, tmp8, tmp9); \
	r1 = r13[0]; r3 = r13[1]; \
	r2 = r24[0]; r4 = r24[1]; \
	r5 = r57[0]; r7 = r57[1]; \
	r6 = r68[0]; r8 = r68[1]; \
}
#else
#define xs4c_mulf(r1,r2,r3,r4,r5,r6,r7,r8,m1,m2,m3,m4,m5,m6,m7,m8) { \
	double tmp8, tmp9; \
	xs_complex_mult(r1, r2, f64ptr(m1+rbp), f64ptr(m2+rbp), tmp8, tmp9); \
	xs_complex_mult(r3, r4, f64ptr(m3+rbp), f64ptr(m4+rbp), tmp8, tmp9); \
	xs_complex_mult(r5, r6, f64ptr(m5+rbp), f64ptr(m6+rbp), tmp8, tmp9); \
	xs_complex_mult(r7, r8, f64ptr(m7+rbp), f64ptr(m8+rbp), tmp8, tmp9); \
}
#endif
#define xp4c_mulf(r1,r2,r3,r4,r5,r6,r7,r8,m1,m2,m3,m4,m5,m6,m7,m8) { \
	vec2f64 tmp8, tmp9; \
	xp_complex_mult(r1, r2, xptr(m1+rbp), xptr(m2+rbp), tmp8, tmp9); \
	xp_complex_mult(r3, r4, xptr(m3+rbp), xptr(m4+rbp), tmp8, tmp9); \
	xp_complex_mult(r5, r6, xptr(m5+rbp), xptr(m6+rbp), tmp8, tmp9); \
	xp_complex_mult(r7, r8, xptr(m7+rbp), xptr(m8+rbp), tmp8, tmp9); \
}

#define xs4c_mulf_r(r1,r2,r3,r4,r5,r6,r7,r8,m1,m2,m3,m4,m5,m6,m7,m8) { \
	double tmp8, tmp9; \
	xs_complex_mult(r1, r2, m1, m2, tmp8, tmp9); \
	xs_complex_mult(r3, r4, m3, m4, tmp8, tmp9); \
	xs_complex_mult(r5, r6, m5, m6, tmp8, tmp9); \
	xs_complex_mult(r7, r8, m7, m8, tmp8, tmp9); \
}
#define xp4c_mulf_r(r1,r2,r3,r4,r5,r6,r7,r8,m1,m2,m3,m4,m5,m6,m7,m8) { \
	vec2f64 tmp8, tmp9; \
	xp_complex_mult(r1, r2, m1, m2, tmp8, tmp9); \
	xp_complex_mult(r3, r4, m3, m4, tmp8, tmp9); \
	xp_complex_mult(r5, r6, m5, m6, tmp8, tmp9); \
	xp_complex_mult(r7, r8, m7, m8, tmp8, tmp9); \
}


