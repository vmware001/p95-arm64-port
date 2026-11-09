/* Contains converted code in gcdhlp(...)  */
/* Copyright 1995-2026 Mersenne Research, Inc.  All rights reserved */
/* Author:  George Woltman */
/*  Email: woltman@alum.mit.edu */

/*  This file implements helper routines for the ECM code */

#include "inttypes.h"
#include <stdio.h>

#ifdef _MSC_VER
#include <intrin.h>
#define shldq(A,B,C) __shiftleft128(B,A,C)
#define shrdq(A,B,C) __shiftright128(A,B,C)
#else

#if defined(__i386__) || defined(__x86_64__) || defined(_M_X64) || defined(_M_IX86)
#include <x86intrin.h>
#endif

typedef unsigned _BitInt(128) uint128;
typedef signed _BitInt(128) int128;



#if 0 //def __x86_64__
static inline uint64_t shldq(uint64_t a, uint64_t b, uint32_t c) {
    __asm__ ( "shrdq %%cl, %1, %0" : "+r" (a) : "r" (b), "c"(c) : "cc" );
    return a;
}
static inline uint64_t shrdq(uint64_t a, uint64_t b, uint32_t c) {
    __asm__ ( "shldq %%cl, %1, %0" : "+r" (a) : "r" (b), "c"(c) : "cc" );
    return a;
}
#elif 0
static inline uint64_t shldq(uint64_t a, uint64_t b, uint32_t c) {
	uint128 r = b | ((uint128)a << 64);
	return r << (c&63) >> 64;
}
static inline uint64_t shrdq(uint64_t a, uint64_t b, uint32_t c) {
	uint128 r = a | ((uint128)b << 64);
	return r >> (c&63);
}
#else
static inline uint64_t shldq(uint64_t a, uint64_t b, uint32_t c) {
	//c &= 63;
	return a<<c | b>>1>>(63^c);
}
static inline uint64_t shrdq(uint64_t a, uint64_t b, uint32_t c) {
	//c &= 63;
	return a>>c | b<<1<<(63^c);
}
#endif
static inline uint64_t _mul128(int64_t a, int64_t b, int64_t *h) {
	int128 o = a * (int128)b;
	*h = o >> 64;
	return (uint64_t)o;
}
static inline uint64_t _umul128(uint64_t a, uint64_t b, uint64_t *h) {
	uint128 o = a * (uint128)b;
	*h = o >> 64;
	return (uint64_t)o;
}

static inline int64_t _div128(int64_t high, uint64_t low, int64_t d, int64_t *r) {
#ifdef __x86_64__
	uint64_t q;
	// the "r"(d) can be "rm", but clang is not smart
	__asm__ ("idivq %[v]": "=a" (q), "=d" (*r) : [v] "r" (d), "a" (low), "d" (high) : "cc");
	return q;
#else
	int128 n = ((uint128)(high) << 64) | low;
	*r = n % d;
	return n / d;
#endif
}
static inline uint64_t _udiv128(uint64_t high, uint64_t low, uint64_t d, uint64_t *r) {
#ifdef __x86_64__
	uint64_t q;
	__asm__ ("divq %[v]": "=a" (q), "=d" (*r) : [v] "r" (d), "a" (low), "d" (high) : "cc");
	return q;
#else
	uint128 n = ((uint128)(high) << 64) | low;
	*r = n % d;
	return n / d;
#endif
}

static inline int64_t div_u64_by_s64(uint64_t n, int64_t d) {
#ifdef __x86_64__
	uint64_t q;
	__asm__ ("xor %%edx,%%edx\n idivq %[v]": "=a" (q) : [v] "r" (d), "a" (n) : "cc", "%rdx");
	return q;
#else
	// Good for ARM64 CPU.
	int64_t r1 = n / (uint64_t)(d < 0 ? -d : d);
	if(d < 0)
		r1 = -r1;
	return r1;
#endif
}

#endif

#if 1
void addhlp (uint32_t *res, uint32_t *carry, uint32_t val) {
	uint32_t c = 0;
	*res = __builtin_addc(*res, val, c, &c);
	*carry = __builtin_addc(*carry, 0, c, &c);
}
void subhlp (uint32_t *res, uint32_t *carry, uint32_t val) {
	uint32_t c = 0;
	*res = __builtin_subc(*res, val, c, &c);
	*carry = __builtin_subc(*carry, 0, c, &c);
}
void muladdhlp (uint32_t *res, uint32_t *carryl, uint32_t *carryh, uint32_t val1, uint32_t val2) {
	uint64_t m = val1 * (uint64_t)val2;
	uint32_t c = 0;
	*res = __builtin_addc(*res, (uint32_t)m, c, &c);
	*carryl = __builtin_addc(*carryl, m >> 32, c, &c);
	*carryh = __builtin_addc(*carryh, 0, c, &c);
}
void muladd2hlp (uint32_t *res, uint32_t *carryl, uint32_t *carryh, uint32_t val1, uint32_t val2) {
	uint64_t m = val1 * (uint64_t)val2;
	uint32_t m1 = (uint32_t)m;
	uint32_t m2 = m >> 32;
	uint32_t c = 0;
	m1 = __builtin_addc(m1, m1, 0, &c);
	m2 = __builtin_addc(m2, m2, c, &c);
	*carryh = __builtin_addc(*carryh, 0, c, &c);
	*res = __builtin_addc(*res, m1, 0, &c);
	*carryl = __builtin_addc(*carryl, m2, c, &c);
	*carryh = __builtin_addc(*carryh, 0, c, &c);
}
void mulsubhlp (uint32_t *res, uint32_t *carryl, uint32_t *carryh, uint32_t val1, uint32_t val2) {
	uint64_t m = val1 * (uint64_t)val2;
	uint32_t c = 0;
	*res = __builtin_subc(*res, (uint32_t)m, c, &c);
	*carryl = __builtin_subc(*carryl, m >> 32, c, &c);
	*carryh = __builtin_subc(*carryh, 0, c, &c);
}
#else
void addhlp (uint32_t *res, uint32_t *carry, uint32_t val) {
	uint64_t a = *res + (((uint64_t) *carry) << 32);
	a += val;
	*res = (uint32_t)a;
	*carry = a >> 32;
	//puts("+");
}
void subhlp (uint32_t *res, uint32_t *carry, uint32_t val) {
	uint64_t a = *res + (((uint64_t) *carry) << 32);
	a -= val;
	*res = (uint32_t)a;
	*carry = a >> 32;
	//puts(" -");
}
void muladdhlp (uint32_t *res, uint32_t *carryl, uint32_t *carryh, uint32_t val1, uint32_t val2) {
#if 1
	uint64_t a = val1 * (uint64_t)val2;
	unsigned _BitInt(96) r =
		(((uint64_t) *carryl) << 32) | *res |
		(((unsigned _BitInt(96)) *carryh) << 64);
	r += a;
	*res = (uint32_t)r;
	*carryl = r >> 32;
	*carryh = r >> 64;
#else
	uint64_t a = val1 * (uint64_t)val2;
	uint64_t b = (((uint64_t) *carryl) << 32) + *res;
	unsigned long long c = 0;
	unsigned int c4 = 0;
	a = __builtin_addcll(a, b, 0, &c);
	*carryh = __builtin_addc(*carryh, 0, c, &c4);
	*res = (uint32_t)a;
	*carryl = a >> 32;
#endif
	//puts("   *+");
}
void muladd2hlp (uint32_t *res, uint32_t *carryl, uint32_t *carryh, uint32_t val1, uint32_t val2){
#if 1
	uint64_t a = val1 * (uint64_t)val2;
	unsigned _BitInt(96) r =
		(((uint64_t) *carryl) << 32) | *res |
		(((unsigned _BitInt(96)) *carryh) << 64);
	r += a;
	r += a;
	*res = (uint32_t)r;
	*carryl = r >> 32;
	*carryh = r >> 64;
#else
	uint64_t a = val1 * (uint64_t)val2;
	uint64_t b = (((uint64_t)*carryl) << 32) + *res;
	unsigned long long c = 0;
	unsigned int c4 = 0;
	a = __builtin_addcll(a, b, 0, &c);
	*carryh = __builtin_addc(*carryh, 0, c, &c4);
	a = __builtin_addcll(a, b, 0, &c);
	*carryh = __builtin_addc(*carryh, 0, c, &c4);
	*res = (uint32_t)a;
	*carryl = a >> 32;
#endif
	//puts("      *++");
}
void mulsubhlp (uint32_t *res, uint32_t *carryl, uint32_t *carryh, uint32_t val1, uint32_t val2){
#if 1
	uint64_t a = val1 * (uint64_t)val2;
	unsigned _BitInt(96) r =
		(((uint64_t) *carryl) << 32) | *res |
		(((unsigned _BitInt(96)) *carryh) << 64);
	r -= a;
	*res = (uint32_t)r;
	*carryl = r >> 32;
	*carryh = r >> 64;
#else
	uint64_t a = val1 * (uint64_t)val2;
	uint64_t b = (((uint64_t)*carryl) << 32) + *res;
	unsigned long long c = 0;
	unsigned int c4 = 0;
	a = __builtin_subcll(a, b, 0, &c);
	*carryh = __builtin_subc(*carryh, 0, c, &c4);
	*res = (uint32_t)a;
	*carryl = a >> 32;
#endif
	//puts("         *-");
}
#endif

const double TWOPOW32 = 4294967296.0;		 // 2^32
const double BIGVAL = 9223372036854775808.0;  // 2^63

int gcdhlp (uint32_t ulen, uint32_t *udata, uint32_t vlen, uint32_t *vdata, unsigned int *struct_ptr) {
	//puts("         gcdhlp");
	uint64_t rax,rcx,rbx,rsi,rbp,r8,r9,r10,r12,r13,r14,r15,t;
	uintptr_t c = ulen;
	r8 = vlen;

	r14 = udata[c-1];	/* U[Ulen-1] */
	r15 = 0;		/* Zero V in case U and V not same len */
	if(c == vlen)			/* Load top V word if U and V same len */
		r15 = vdata[c-1];	/* V[Ulen-1] */

	if(c == 1)			/* Are there more words to shift */
		goto simple;			/* bits from? */

	r14 <<= 32;
	r15 <<= 32;
	rax = udata[c-2];	/* U[Ulen-2] */
	rbx = vdata[c-2];	/* V[Ulen-2] */
	r14 += rax;
	r15 += rbx;

	if(c == 2)			/* Are there more words to shift */
		goto simple;			/* bits from? */

	rax = udata[c-3];	/* U[Ulen-3] */
	rbx = vdata[c-3];	/* V[Ulen-3] */

	rcx = __builtin_clzll(r14);		/* Count bits to shift U */
	/*cl ^= 63;*/			/* Turn bit # into a shift count */

	rax <<= 32;
	rbx <<= 32;
	r14 = shldq(r14, rax, rcx); //shld	r14, rax, cl;		/* Shift U */
	r15 = shldq(r15, rbx, rcx); //shld	r15, rbx, cl;		/* Shift V */

	/* Init extended GCD information */

	r8 = 1;			/* A */
	r9 = 0;			/* B */
	r10 = 0;		/* C */
	r13 = 1;			/* D */
	r12 = 0;		/* ODD */

	/* Do as many operations as we can constrained by a maximum 32-bit result */
	/*	FPU contains U, V, D, B, C, A */
	/* As Knuth suggests: */
	/*	Compute (U-B)/(V+D), the smaller quotient */
	/*	Compute (U+A)/(V-C), the larger quotient, break if not equal */
	/*	Set newB = D, newD = B+Q*D, break if newD won't fit in 32 bits */
	/*	Set newA = C, newC = A+Q*C */
	/*	Set newU = V, newV = U-Q*V */

	rsi = 1;			/* Compute 2^32 */
	rsi <<= 32;
	r14 >>= 1;			/* Lose a bit so we can use signed */
	r15 >>= 1;			/* multiply and divide instructions */

	dloop:
	rax = r14;		/* U */
	rax -= r9;			/* U-B */
	rcx = r15;		/* V */
	rcx += r13;		/* V+D */
	rax = div_u64_by_s64(rax, (intptr_t)rcx);	/* Q = (U-B)/(V+D) */

	rcx = r15;		/* V */
	rcx -= r10;		/* V-C */
	rbp = rax;		/* Copy Q */
	rbp *= rcx;		/* Q*(V-C) */
	rbp += rcx;		/* (Q+1)*(V-C) */
	rbx = (r14+r8);		/* U+A */
	if(rbx >= rbp)		/* Compare U+A to (Q+1)*(V-C) */
		goto lpdone;			/* Break if above or equal */

	rbx = r13;		/* D */
	rbx *= rax;		/* Q*D */
	rbx += r9;			/* B+Q*D (newD) */
	if(rbx >= rsi)		/* Will new D be too large? */
		goto lpdone;			/* Break if newD is 32-bits or more */
	r9 = r13;			/* newB = D */
	r13 = rbx;		/* newD = B+Q*D */

	rcx = r10;		/* C */
	rcx *= rax;		/* Q*C */
	rcx += r8;			/* A+Q*C (newC) */
	r8 = r10;			/* newA = C */
	r10 = rcx;		/* newC = A+Q*C */

	rax *= r15;		/* Q*V */
	r14 -= rax;		/* U-Q*V (newV) */
	t = r14; r14 = r15; r15 = t;	/* newU = V, newV = U-Q*V */

	r12 ^= 1;			/* Flip ODD */

	/*	Compute (U-A)/(V+C), the smaller quotient */
	/*	Compute (U+B)/(V-D), the larger quotient, break if not equal */
	/*	Set newB = D, newD = B+Q*D, break if newD won't fit in 32 bits */
	/*	Set newA = C, newC = A+Q*C */
	/*	Set newU = V, newV = U-Q*V */

	rax = r14;		/* U */
	rax -= r8;			/* U-A */
	rcx = r15;		/* V */
	rcx += r10;		/* V+C */
	rax = div_u64_by_s64(rax, (intptr_t)rcx); /* Q = (U-A)/(V+C) */

	rcx = r15;		/* V */
	rcx -= r13;		/* V-D */
	rbp = rax;		/* Copy Q */
	rbp *= rcx;		/* Q*(V-D) */
	rbp += rcx;		/* (Q+1)*(V-D) */
	rbx = (r14+r9);		/* U+B */
	if(rbx >= rbp)		/* Compare U+B to (Q+1)*(V-D) */
		goto lpdone;			/* Break if above or equal */

	rbx = r13;		/* D */
	rbx *= rax;		/* Q*D */
	rbx += r9;			/* B+Q*D (newD) */
	if(rbx >= rsi)		/* Will new D be too large? */
		goto lpdone;			/* Break if newD is 32-bits or more */
	r9 = r13;			/* newB = D */
	r13 = rbx;		/* newD = B+Q*D */

	rcx = r10;		/* C */
	rcx *= rax;		/* Q*C */
	rcx += r8;			/* A+Q*C (newC) */
	r8 = r10;			/* newA = C */
	r10 = rcx;		/* newC = A+Q*C */

	rax *= r15;		/* Q*V */
	r14 -= rax;		/* U-Q*V (newV) */
	t = r14; r14 = r15; r15 = t;	/* newU = V, newV = U-Q*V */

	r12 ^= 1;			/* Flip ODD */
	goto	dloop;

	/* The single precision case: */

	/* Init extended GCD information */

	simple:	r8 = 1;			/* A */
	r9 = 0;			/* B */
	r10 = 0;		/* C */
	r13 = 1;			/* D */
	r12 = 0;		/* ODD */

	/*	Compute U/V, the quotient */
	/*	Set newB = D, newD = B+Q*D, break if newD won't fit in 32 bits */
	/*	Set newA = C, newC = A+Q*C */
	/*	Set newU = V, newV = U-Q*V */

	rsi = 1;			/* Compute 2^32 */
	rsi <<= 32;
	sloop:
	rax = r14;		/* U */
	rax = rax / r15;			/* Q = U/V */
	rbp = rax;		/* Copy Q */

	rax *= r13;			/* Q*D */
	rax += r9;			/* B+Q*D */
	if(rax >= rsi)		/* Will new D be too large? */
		goto lpdone;		/* Break if D >= 2^32 */
	r9 = r13;			/* newB = D */
	r13 = rax;		/* newD = B+Q*D */

	rax = r10;		/* C */
	rax *= rbp;			/* Q*C */
	rax += r8;			/* A+Q*C (newC) */
	r8 = r10;			/* newA = C */
	r10 = rax;		/* newC = A+Q*C */

	rax = r15;		/* V */
	rax *= rbp;			/* Q*V */
	r14 -= rax;		/* U-Q*V (newV) */
	t = r14; r14 = r15; r15 = t;		/* newU = V, newV = U-Q*V */

	r12 ^= 1;			/* Flip ODD */

	if(r15)		/* Loop if V is not zero */
		goto sloop;

	/* Copy extended GCD info to globals */

	lpdone:
	struct_ptr[0] = r8; /* Save A */
	struct_ptr[1] = r9; /* Save B */
	struct_ptr[2] = r10; /* Save C */
	struct_ptr[3] = r13; /* Save D */
	struct_ptr[4] = r12; /* Save ODD */

	return r9 != 0 ? 1 : 0;
}
