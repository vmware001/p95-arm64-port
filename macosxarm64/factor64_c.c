/*
 * factor64_c.c -- portable C implementation of Prime95/mprime's x86-64
 * trial-factoring assembly (prime95/factor64.asm), for the Linux ARM64 port.
 *
 * Upstream ships linux64/factor64.o (prebuilt x86-64 ELF of factor64.asm).
 * This file replaces it with portable C.  It implements exactly the entry
 * points that commonb.c's 64-bit trial factoring driver calls:
 *
 *	int factor64_pass_setup (struct facasm_data *);
 *	int factor64_small      (struct facasm_data *);
 *	int factor64_sieve      (struct facasm_data *);
 *	int factor64_tf         (struct facasm_data *);
 *	int factor64            (struct facasm_data *);   (32-bit signature; unused)
 *
 * The C/asm "register file" (struct facasm_data) is private to commonb.c, so
 * it is duplicated bit-for-bit below.  See _work/factor64/CONTRACT.md for the
 * complete reverse-engineered contract, and _work/factor64/RESULTS.md for
 * validation notes.
 *
 * Copyright 1995-2022 Mersenne Research, Inc.  (original algorithm)
 * Copyright 2026 (C reimplementation for aarch64 port)
 *
 * This program is free software; you can redistribute it and/or modify it
 * under the terms of the GNU General Public License as published by the Free
 * Software Foundation; either version 2 of the License, or (at your option)
 * any later version.
 */

#include <stdint.h>
#include <stddef.h>
#include <string.h>

#ifdef X86_64

/* ----------------------------------------------------------------------- */
/* The C / assembly communication structure.                               */
/* Verbatim copy of the 64-bit definition in commonb.c (lines 4209-4235);  */
/* keep in sync with it.  It is allocated 64/128-byte aligned and memcpy'd  */
/* to auxiliary threads, so layout (3360 bytes) must match exactly.        */
/* ----------------------------------------------------------------------- */

#define SIEVE_SIZE_IN_BYTES	(12 * 1024)	/* 98304 bits per sieve area */
#define SIEVE_BITS		(SIEVE_SIZE_IN_BYTES * 8)
#define SIEVE_WORDS		(SIEVE_SIZE_IN_BYTES / 8)

struct facasm_data {
	uint64_t p;			/* Exponent of Mersenne number to factor */
	uint64_t facdists[65];		/* Multiples of 120 * p * num_siever_allocs */
	uint64_t facdist12K;		/* Factor distance covered by a 12K sieve */
	uint64_t savefac1;		/* LSW or first factor in sieve */
	uint64_t savefac0;		/* MSW or first factor in sieve */
	void	*sieve;			/* Area to sieve or already sieved area to TF */
	void	*primearray;		/* Array of sieving primes */
	void	*offsetarray;		/* Array of bits-to-clear for each sieving prime */
	void	*initsieve;		/* Array used to initialize sieve */
	double	TWO_TO_FACSIZE_PLUS_62;	/* Constant used in SSE2/AVX2/AVX512 integer TF code */

	uint32_t FACHSW;		/* High word of found factor */
	uint32_t FACMSW;		/* Middle word of found factor */
	uint32_t FACLSW;		/* Low word of found factor */
	uint32_t cpu_flags;		/* Copy of CPU_FLAGS */
	uint32_t SSE2_LOOP_COUNTER;	/* Counter used in SSE2/AVX2/FMA/AVX512 TF code */
	uint32_t FMA_INITVAL_TYPE;	/* Type of initval (more or less than 2^split).  Used in FMA code. */
	uint32_t FMA_64_MINUS_LOW_BITS;	/* 64 minus bits in a low word.  Used in FMA code. */
	uint32_t initstart;		/* First byte in initsieve to copy (see ASM code) */
	uint32_t alternate_sieve_count;	/* Number of 13-small-prime AVX2 sieve sections to process */

	uint32_t pad[3];		/* Pad to next 64-byte cache line */
	uint32_t xmm_data[80];		/* XMM/YMM/ZMM data initialized in C code */

	uint32_t other_asm_temps_and_consts[600];
};

/* Layout assertions: these offsets are what factor64.asm expects.  The 96-bit
 * factor is (savefac0 << 64) | savefac1; FACHSW is its top 32 bits. */

#define AD_ASSERT(name, cond)	typedef char ad_assert_##name[(cond) ? 1 : -1]

AD_ASSERT (size, sizeof (struct facasm_data) == 3360);
AD_ASSERT (facdists, offsetof (struct facasm_data, facdists) == 8);
AD_ASSERT (facdist12K, offsetof (struct facasm_data, facdist12K) == 66 * 8);
AD_ASSERT (savefac1, offsetof (struct facasm_data, savefac1) == 67 * 8);
AD_ASSERT (savefac0, offsetof (struct facasm_data, savefac0) == 68 * 8);
AD_ASSERT (sieve, offsetof (struct facasm_data, sieve) == 69 * 8);
AD_ASSERT (FACHSW, offsetof (struct facasm_data, FACHSW) == 74 * 8);
AD_ASSERT (initstart, offsetof (struct facasm_data, initstart) == 74 * 8 + 7 * 4);
AD_ASSERT (xmm_data, offsetof (struct facasm_data, xmm_data) == 80 * 8);
AD_ASSERT (other_asm, offsetof (struct facasm_data, other_asm_temps_and_consts) == 120 * 8);

#undef AD_ASSERT

typedef unsigned __int128 u128;

/* 2^44: upstream factor64.asm's `smfacendpt`, the upper limit of the
 * brute-force (non-sieved) small-factor search. */
#define SMALL_FACTOR_LIMIT	((uint64_t) 1 << 44)

/* ----------------------------------------------------------------------- */
/* Small helpers                                                           */
/* ----------------------------------------------------------------------- */

static inline int clz64 (uint64_t x) { return __builtin_clzll (x); }
static inline int ctz64 (uint64_t x) { return __builtin_ctzll (x); }

/* (hi * 2^64 + lo) mod s, for s < 2^32.  q is at most 96 bits. */

static inline uint32_t mod96_u32 (uint64_t hi, uint64_t lo, uint32_t s)
{
	uint32_t rh = (uint32_t) (hi % s);
	/* rh < s < 2^32, so rh << 64 < 2^96: fits in u128 */
	u128 t = ((u128) rh << 64) | lo;
	return (uint32_t) (t % s);
}

/* Modular inverse of a (mod m), both odd/positive, m < 2^32.  Extended Euclid. */

static uint32_t modinv_u32 (uint32_t a, uint32_t m)
{
	int64_t t = 0, newt = 1;
	int64_t r = m, newr = a;
	while (newr != 0) {
		int64_t q = r / newr;
		int64_t tmp = t - q * newt; t = newt; newt = tmp;
		tmp = r - q * newr; r = newr; newr = tmp;
	}
	if (t < 0) t += m;
	return (uint32_t) t;
}

/* ----------------------------------------------------------------------- */
/* Modular arithmetic used to test q | 2^p - 1.                            */
/*                                                                         */
/* Two paths:                                                              */
/*   q < 2^64      : 1-limb Montgomery multiplication (u128 products).     */
/*   2^64 <= q<2^96: 2-limb Montgomery multiplication (CIOS).              */
/* Both are branch-light and need no 128/128-bit division, which matters    */
/* on aarch64 (no hardware 128-bit divide).                                */
/* ------------------------------------------------------------------------ */

/* n^-1 mod 2^64 for odd n (Newton iteration; each step doubles the bits). */

static inline uint64_t mont_inv64 (uint64_t n)
{
	uint64_t x = n;			/* correct to 3 bits */
	x *= 2 - n * x;			/* 6 bits */
	x *= 2 - n * x;			/* 12 bits */
	x *= 2 - n * x;			/* 24 bits */
	x *= 2 - n * x;			/* 48 bits */
	x *= 2 - n * x;			/* 96 bits (capped to 64) */
	return x;
}

/* r = a * b / 2^64 mod n, with a,b < n, n odd, np = -n^-1 mod 2^64. */

static inline uint64_t mont_mul64 (uint64_t a, uint64_t b, uint64_t n, uint64_t np)
{
	u128 t = (u128) a * b;
	uint64_t t_lo = (uint64_t) t;
	uint64_t t_hi = (uint64_t) (t >> 64);
	uint64_t m = t_lo * np;
	u128 mn = (u128) m * n;
	/* t_lo + low(mn) == 0 (mod 2^64) exactly; carry into the high half. */
	uint64_t carry = (t_lo != 0);
	u128 u = (u128) t_hi + (uint64_t) (mn >> 64) + carry;
	if (u >= n) u -= n;
	return (uint64_t) u;
}

/* 2^64 mod n for odd n > 1. */

static inline uint64_t pow2_64_mod (uint64_t n)
{
	uint64_t r = (UINT64_MAX % n) + 1;
	if (r == n) r = 0;
	return r;
}

/* Return 1 if 2^p == 1 (mod n), n odd, n < 2^64.  p >= 1. */

static int pow2p_is_one_64 (uint64_t n, uint64_t p)
{
	uint64_t np = (uint64_t) 0 - mont_inv64 (n);
	uint64_t R = pow2_64_mod (n);		/* = 1 in Montgomery domain */
	uint64_t base = R + R;			/* = 2 in Montgomery domain */
	uint64_t r;
	int bit;
	if (base >= n) base -= n;
	if (p == 0) return (R == (uint64_t) (1 % n));
	bit = 63 - clz64 (p);
	r = base;				/* leading 1 bit */
	for (bit--; bit >= 0; bit--) {
		r = mont_mul64 (r, r, n, np);
		if ((p >> bit) & 1) r = mont_mul64 (r, base, n, np);
	}
	return (r == R);			/* 2^p*R == R  <=>  2^p == 1 */
}

/* Generic CIOS Montgomery multiply for up to 2 limbs (n < 2^128).  r = a*b*R^-1
 * mod n where R = 2^(64*N).  a,b,n have N limbs; n odd, np = -n[0]^-1 mod 2^64. */

static void mont_mul_n (uint64_t *r, const uint64_t *a, const uint64_t *b,
			const uint64_t *n, int N, uint64_t np)
{
	uint64_t t[4];
	int i, j;
	for (i = 0; i <= N; i++) t[i] = 0;
	for (i = 0; i < N; i++) {
		u128 acc = 0;
		uint64_t m;
		for (j = 0; j < N; j++) {
			acc += (u128) a[i] * b[j] + t[j];
			t[j] = (uint64_t) acc;
			acc >>= 64;
		}
		acc += t[N];
		t[N] = (uint64_t) acc;
		t[N + 1] = (uint64_t) (acc >> 64);
		m = t[0] * np;
		acc = (u128) m * n[0] + t[0];		/* low limb cancels to 0 */
		acc >>= 64;
		for (j = 1; j < N; j++) {
			acc += (u128) m * n[j] + t[j];
			t[j - 1] = (uint64_t) acc;
			acc >>= 64;
		}
		acc += t[N];
		t[N - 1] = (uint64_t) acc;
		t[N] = t[N + 1] + (uint64_t) (acc >> 64);
	}
	/* Result < 2n; one conditional subtraction. */
	if (t[N] || t[N - 1] > n[N - 1] ||
	    (t[N - 1] == n[N - 1] && t[N - 2] >= n[N - 2])) {
		uint64_t borrow = (t[0] < n[0]);
		r[0] = t[0] - n[0];
		r[1] = t[1] - n[1] - borrow;
	} else {
		r[0] = t[0];
		r[1] = t[1];
	}
}

/* Return 1 if 2^p == 1 (mod n), 2^64 <= n < 2^96, n odd.  p >= 1. */

static int pow2p_is_one_128 (unsigned __int128 n, uint64_t p)
{
	uint64_t nl[2], R[2], base[2], r[2], np;
	uint64_t hi = (uint64_t) (n >> 64);
	int bit, i;
	if (hi == 0) return pow2p_is_one_64 ((uint64_t) n, p);
	nl[0] = (uint64_t) n;
	nl[1] = hi;
	np = (uint64_t) 0 - mont_inv64 (nl[0]);
	/* R = 2^128 mod n, by doubling (n < 2^96 so a 2-limb register suffices) */
	R[0] = 1; R[1] = 0;
	for (i = 0; i < 128; i++) {
		uint64_t new0 = R[0] << 1;
		uint64_t new1 = (R[1] << 1) | (R[0] >> 63);
		if (new1 > nl[1] || (new1 == nl[1] && new0 >= nl[0])) {
			uint64_t borrow = (new0 < nl[0]);
			new0 -= nl[0];
			new1 -= nl[1] + borrow;
		}
		R[0] = new0; R[1] = new1;
	}
	/* base = 2 in Montgomery domain */
	{
		uint64_t new0 = R[0] << 1;
		uint64_t new1 = (R[1] << 1) | (R[0] >> 63);
		if (new1 > nl[1] || (new1 == nl[1] && new0 >= nl[0])) {
			uint64_t borrow = (new0 < nl[0]);
			new0 -= nl[0];
			new1 -= nl[1] + borrow;
		}
		base[0] = new0; base[1] = new1;
	}
	bit = 63 - clz64 (p);
	r[0] = base[0]; r[1] = base[1];
	for (bit--; bit >= 0; bit--) {
		mont_mul_n (r, r, r, nl, 2, np);
		if ((p >> bit) & 1) mont_mul_n (r, r, base, nl, 2, np);
	}
	return (r[0] == R[0] && r[1] == R[1]);
}

/* Return 1 if q divides 2^p - 1 (q odd). */

static int divides_pow2m1 (unsigned __int128 q, uint64_t p)
{
	if (p == 0) return 1;
	if ((q >> 64) == 0) return pow2p_is_one_64 ((uint64_t) q, p);
	return pow2p_is_one_128 (q, p);
}

/* ----------------------------------------------------------------------- */
/* factor64_pass_setup                                                     */
/*                                                                         */
/* Upstream asm pre-computes squaring-loop constants (initval/initshift/   */
/* shifter/initdiv/shift_count) for the TF kernels, keyed off the bit size  */
/* of savefac.  The C port recomputes everything per candidate, so there   */
/* is nothing to do here.  The caller ignores the return value.            */
/* ----------------------------------------------------------------------- */

int factor64_pass_setup (struct facasm_data *asm_data)
{
	(void) asm_data;
	return (0);
}

/* ----------------------------------------------------------------------- */
/* factor64_small                                                          */
/*                                                                         */
/* Brute-force trial dividing for factors below 2^44 without a sieve.      */
/* savefac1 holds the first factor to test (set by the caller to 2*p+1 on  */
/* the first call; incremented by 2*p between calls).  Same semantics as   */
/* factor64.asm:                                                           */
/*   - factors that are 3 or 5 mod 8 cannot divide 2^p-1 (p odd) -> skip;  */
/*   - no primality test (upstream does not do one either);                */
/*   - on a divisor: savefac1 = factor, return 1;                          */
/*   - after q >= 2^44 (or ESC check opportunity): return 2.               */
/* ----------------------------------------------------------------------- */

int factor64_small (struct facasm_data *asm_data)
{
	uint64_t p = asm_data->p;
	uint64_t step = p + p;
	uint64_t q = asm_data->savefac1;

	if (step == 0) return (2);
	/* Same loop shape as factor64.asm: the first value is always tested,
	 * then q += 2p while q < 2^44. */
	for ( ; ; ) {
		uint64_t r8 = q & 6;
		/* Only q = 1 or 7 (mod 8) can divide 2^p-1 for odd prime p. */
		if (r8 == 0 || r8 == 6) {
			if (divides_pow2m1 ((u128) q, p)) {
				asm_data->savefac1 = q;
				return (1);
			}
		}
		q += step;
		if (q >= SMALL_FACTOR_LIMIT) break;
	}
	return (2);
}

/* ----------------------------------------------------------------------- */
/* factor64_sieve                                                          */
/*                                                                         */
/* Fill the 12KB sieve for one area: bit k set  <=>  q = first + k*step may */
/* be a factor (not eliminated by a small prime).  first is a 96-bit value  */
/* (savefac0:savefac1), step = facdists[1].  Then advance savefac by        */
/* facdist12K, which is the first factor of the next area for this siever   */
/* (commonb.c reads it back as next_sieve_first_factor[]).                  */
/*                                                                         */
/* The primes 2,3,5 need not be sieved: every candidate in a pass is a      */
/* unit mod 120.  Primes dividing step make every candidate divisible by    */
/* that prime or by none (prime p: never; 7/11/13: commonb.c guarantees the */
/* pass start is not divisible).  All remaining primes come from            */
/* primearray[] (up to MaxTFSievePrime, ~1e6 by default) with offsets whose */
/* per-area recurrence is  off' = (off - 98304) mod s.                     */
/* ----------------------------------------------------------------------- */

/* The initsieve primes {7,11,13,17,19,23}: commonb.c deliberately omits the
 * prime factors of initsieve_primes from primearray because the asm's
 * initsieve pre-clears them.  The C port sieves them explicitly here. */
static const uint32_t init_sieve_primes[6] = { 7, 11, 13, 17, 19, 23 };

static inline void sieve_clear_prime (uint64_t *bits, uint32_t s, uint32_t off)
{
	uint32_t j;
	/* off is in [0,s).  Clear bits off, off+s, ... < 98304. */
	for (j = off; j < SIEVE_BITS; j += s)
		bits[j >> 6] &= ~((uint64_t) 1 << (j & 63));
}

int factor64_sieve (struct facasm_data *asm_data)
{
	uint64_t *bits = (uint64_t *) asm_data->sieve;
	uint32_t *primearray = (uint32_t *) asm_data->primearray;
	uint32_t *offsetarray = (uint32_t *) asm_data->offsetarray;
	uint64_t step = asm_data->facdists[1];
	uint64_t hi = asm_data->savefac0;
	uint64_t lo = asm_data->savefac1;
	int i;

	for (i = 0; i < SIEVE_WORDS; i++) bits[i] = ~(uint64_t) 0;

	/* Sieve the initsieve primes that are not already covered by primearray. */
	for (i = 0; i < (int) (sizeof (init_sieve_primes) / sizeof (init_sieve_primes[0])); i++) {
		uint32_t s = init_sieve_primes[i];
		uint32_t firstmod, stepmod, inv, off;
		if (s == 0) continue;
		if (step % s == 0) {
			/* Whole progression has a fixed residue mod s. */
			if (mod96_u32 (hi, lo, s) == 0) {
				memset (bits, 0, SIEVE_SIZE_IN_BYTES);
				goto advance;
			}
			continue;
		}
		firstmod = mod96_u32 (hi, lo, s);
		stepmod = (uint32_t) (step % s);
		inv = modinv_u32 (stepmod, s);
		off = s - (uint32_t) (((uint64_t) firstmod * inv) % s);
		if (off == s) off = 0;
		sieve_clear_prime (bits, s, off);
	}

	/* Sieve all remaining small primes.  offsetarray[] holds this area's
	 * starting bit for each prime; update it in place for the next area. */
	{
		uint32_t *pa, *oa;
		for (pa = primearray, oa = offsetarray; ; pa++, oa++) {
			uint32_t s = *pa, off;
			int64_t newoff;
			if (s == 0) break;
			off = *oa;
			if (off >= s) off %= s;
			sieve_clear_prime (bits, s, off);
			newoff = (int64_t) off - (int64_t) SIEVE_BITS;
			newoff %= (int64_t) s;
			if (newoff < 0) newoff += s;
			*oa = (uint32_t) newoff;
		}
	}

advance:
	/* Advance savefac (96-bit) to the first factor of the next sieve area.
	 * This is the value commonb.c reads back and stores in the siever. */
	{
		u128 q = ((u128) hi << 64) | lo;
		q += asm_data->facdist12K;
		asm_data->savefac1 = (uint64_t) q;
		asm_data->savefac0 = (uint64_t) (q >> 64);
	}
	return (0);
}

/* ----------------------------------------------------------------------- */
/* factor64_tf                                                             */
/*                                                                         */
/* Trial factor one already-sieved 12KB area: for every set bit k (in       */
/* ascending order), test q = savefac0:savefac1 + k*facdists[1].  Report    */
/* the first q dividing 2^p-1 in FACHSW/FACMSW/FACLSW and return 1;        */
/* otherwise return 2.  (commonb.c treats only 1 as "found".)              */
/*                                                                         */
/* Unlike the asm, which selects a size-class kernel from the first        */
/* factor's bit size and bails out at class boundaries (2^58, 2^64, ...),  */
/* the C port tests every candidate in the area.  That can only find more  */
/* factors, never fewer, and it never misses one the asm would find.       */
/* ----------------------------------------------------------------------- */

int factor64_tf (struct facasm_data *asm_data)
{
	const uint64_t step = asm_data->facdists[1];
	const uint64_t p = asm_data->p;
	const uint64_t *bits = (const uint64_t *) asm_data->sieve;
	u128 q = ((u128) asm_data->savefac0 << 64) | asm_data->savefac1;
	int w;

	for (w = 0; w < SIEVE_WORDS; w++) {
		uint64_t b = bits[w];
		while (b != 0) {
			int j = ctz64 (b);
			u128 cand;
			b &= b - 1;
			cand = q + (u128) (uint32_t) (w * 64 + j) * step;
			/* Upstream only has 96 bits of factor state. */
			if ((cand >> 96) != 0) continue;
			if (divides_pow2m1 (cand, p)) {
				asm_data->FACLSW = (uint32_t) cand;
				asm_data->FACMSW = (uint32_t) (cand >> 32);
				asm_data->FACHSW = (uint32_t) (cand >> 64);
				return (1);
			}
		}
	}

	/* No factor.  Advance savefac to the next area's first factor, like the
	 * asm's end-point save (the caller normally overwrites it before the
	 * next call, but keep the state faithful). */
	q += (u128) SIEVE_BITS * step;
	asm_data->savefac1 = (uint64_t) q;
	asm_data->savefac0 = (uint64_t) (q >> 64);
	return (2);
}

/* ----------------------------------------------------------------------- */
/* factor64                                                                */
/*                                                                         */
/* Only the legacy 32-bit build calls this ("do one chunk of work").  The   */
/* 64-bit build drives factor64_sieve/factor64_tf itself.  Provide the      */
/* symbol for completeness: sieve one area and TF it.                       */
/* ----------------------------------------------------------------------- */

int factor64 (struct facasm_data *asm_data)
{
	factor64_sieve (asm_data);
	return (factor64_tf (asm_data));
}

#else /* !X86_64 */

#error "factor64_c.c implements the 64-bit trial-factoring interface; compile with -DX86_64 (as the linux-arm64 Makefile does).  The 32-bit interface (setupf/factor64) is not provided by this port."

#endif /* X86_64 */
