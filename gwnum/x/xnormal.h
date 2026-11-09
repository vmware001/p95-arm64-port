#pragma once
#include "../xbasics.h"

extern void xfunc572(struct gwasm_data *g);


/* Copyright 2000-2024 - Mersenne Research, Inc.  All rights reserved. */ \
/* Author:  George Woltman */ \
/* Email: woltman@alum.mit.edu */ \
/* */ \
/* These macros efficiently implement the normalization to integers */ \
/* and multiplication by two-to-phi powers using SSE2 instructions. */ \
/* */ \
\
\
/* These macros implement the variants of the normalization routines */ \
/* in a non-pipelined way.  It is simply too much work to hand optimize */ \
/* all normalization variants. */ \
\
/* Compute the convolution error and if greater than MAXERR, set MAXERR */ \
\

#define error_check(sse4, xmmreg, tmpreg, errreg) \
	tmpreg = roundx(xmmreg); \
	tmpreg -= xmmreg;		/* This is the convolution error */ \
	tmpreg = absx(tmpreg); \
	maxpd(	errreg, tmpreg);	/* Compute maximum error */

#define error_check_interleaved(sse4, xmmreg1, tmpreg1, xmmreg2, tmpreg2, errreg) \
	tmpreg1 = roundx(xmmreg1); \
	tmpreg2 = roundx(xmmreg2); \
	tmpreg1 -= xmmreg1; \
	tmpreg2 -= xmmreg2; \
	tmpreg1 = absx(tmpreg1); \
	tmpreg2 = absx(tmpreg2); \
	maxpd(	errreg, tmpreg1); \
	maxpd(	errreg, tmpreg2);



\
/* In general, normalization routines calculate: */ \
/*		newFFTvalue = (FFTvalue * const + carry) % base */ \
/*		carry = (FFTvalue * const + carry) / base */ \
/* Since FFTvalue * const can exceed 51 bits, we instead split FFTvalue into: */ \
/*		hi = FFTvalue / base */ \
/*		lo = FFTvalue % base */ \
/* and then calculate: */ \
/*		newFFTvalue = (lo * const + carry) % base */ \
/*		carry = hi * const + (lo * const + carry) / base */ \
/* */ \
/* For the b = 2 case, we can split hi and lo using any power of 2 larger */ \
/* than the FFT base, this allows for some simpler code in this case. */ \
\
\
\
/* These routines split an FFT value and then multiplies it by the small */ \
/* constant.  If we are error-checking the input FFT value has already */ \
/* been rounded. */ \
\

#define mul_by_const2exec(echk, sse4, xmmreg, xmmreghi, xmmtmp, basereg, errreg) \
	xmmreghi = xmmreg + g->u.xmm.XMM_BIGBIGVAL;	/* Round to nearest multiple of 2^25 */ \
	xmmtmp = roundx(xmmreg); \
	xmmreghi -= g->u.xmm.XMM_BIGBIGVAL; \
echk(xmmreg -= xmmtmp);		/* This is the convolution error */ \
echk(xmmreg = absx(xmmreg)); \
echk(maxpd(errreg, xmmreg));		/* Compute maximum error */ \
	xmmtmp -= xmmreghi;	/* xmmtmp now contains low 25 bits */ \
	xmmreg = g->u.xmm.XMM_MULCONST; \
	xmmreghi *= xmmreg;	/* Multiply by the small constant */ \
	xmmreg *= xmmtmp;

#define mul_by_const2noexec(echk, sse4, xmmreg, xmmreghi, xmmtmp, basereg, errreg); \
	xmmreghi = xmmreg * xptr2(g->u.xmm.XMM_LIMIT_INVERSE, basereg); /* Compute FFTvalue / base */ \
	xmmtmp = roundx(xmmreg);	/* Round to an integer */ \
	xmmreghi = roundx(xmmreghi);	/* Round to an integer */ \
echk(xmmreg -= xmmtmp);		/* This is the convolution error */ \
echk(xmmreg = absx(xmmreg)); \
echk(maxpd(errreg, xmmreg));		/* Compute maximum error */ \
	xmmreg = xmmreghi * xptr2(g->u.xmm.XMM_LIMIT_BIGMAX, basereg); \
	xmmtmp -= xmmreg;		/* This is FFTvalue % base */ \
	xmmreg = g->u.xmm.XMM_MULCONST; \
	xmmreghi *= xmmreg;	/* Multiply by the small constant */ \
	xmmreg *= xmmtmp;

#define mul_by_const(ttp, base2, echk, sse4, xmmreg, xmmreghi, xmmtmp, biglitreg, errreg) \
	ttp(mul_by_const2##base2(echk, sse4, xmmreg, xmmreghi, xmmtmp, biglitreg, errreg)) \
	no##ttp(mul_by_const2##base2(echk, sse4, xmmreg, xmmreghi, xmmtmp, 0, errreg))

#define mul_by_const_interleaved(ttp, base2, echk, sse4, xmmreg, xmmreghi, xmmtmp, biglitreg, xmmreg2, xmmreg2hi, xmmtmp2, biglitreg2, errreg) \
	mul_by_const(ttp, base2, echk, sse4, xmmreg, xmmreghi, xmmtmp, biglitreg, errreg); \
	mul_by_const(ttp, base2, echk, sse4, xmmreg2, xmmreg2hi, xmmtmp2, biglitreg2, errreg);


/* */ \
/* These macros do the base2 and nobase2 roundings */ \
/* const - set to exec if the output of mul_by_const needs to be added in */ \
/* xmmval - input: number to round, output: value to store in the FFT */ \
/* xmmcarry - input: part of the next carry if mulbyconst set, output: the next carry */ \
/* xmmtmp - a temporary register */ \
/* */ \
\


// base2_noconst_rounding
#define rounding2execnoexec(xmmval, xmmcarry, xmmtmp, basereg) \
	xmmtmp = xptr2(g->u.xmm.XMM_LIMIT_BIGMAX, basereg); \
	xmmcarry = xmmval + xmmtmp; \
	xmmtmp = xmmcarry - xmmtmp; \
	xmmval -= xmmtmp; \
	xmmcarry *= xptr2(g->u.xmm.XMM_LIMIT_INVERSE, basereg);

// base2_const_rounding
#define rounding2execexec(xmmval, xmmcarry, xmmtmp, basereg) \
	xmmtmp = xptr2(g->u.xmm.XMM_LIMIT_BIGMAX, basereg); \
	xmmtmp += xmmval; \
	xmmcarry += xmmtmp; \
	xmmcarry *= xptr2(g->u.xmm.XMM_LIMIT_INVERSE, basereg); \
	xmmtmp -= xptr2(g->u.xmm.XMM_LIMIT_BIGMAX, basereg); \
	xmmval -= xmmtmp;

// nobase2_noconst_rounding
#define rounding2noexecnoexec(xmmval, xmmcarry, xmmtmp, basereg) \
	xmmval -= XMM_BIGVAL2; \
	xmmcarry = xptr2(g->u.xmm.XMM_LIMIT_INVERSE, basereg); \
	xmmcarry = xmmcarry * xmmval + XMM_BIGVAL2; \
	xmmtmp = xmmcarry - XMM_BIGVAL2; \
	xmmtmp *= xptr2(g->u.xmm.XMM_LIMIT_BIGMAX, basereg); \
	xmmval -= xmmtmp;

// nobase2_const_rounding
#define rounding2noexecexec(xmmval, xmmcarry, xmmtmp, basereg) \
	xmmval -= XMM_BIGVAL2; \
	xmmtmp = xptr2(g->u.xmm.XMM_LIMIT_INVERSE, basereg); \
	xmmtmp = xmmtmp * xmmval + XMM_BIGVAL2; \
	xmmcarry += xmmtmp; \
	xmmtmp -= XMM_BIGVAL2; \
	xmmtmp *= xptr2(g->u.xmm.XMM_LIMIT_BIGMAX, basereg); \
	xmmval -= xmmtmp;


#define rounding1exec(base2, const1, xmmval, xmmcarry, xmmtmp, basereg) \
	rounding2##base2##const1(xmmval, xmmcarry, xmmtmp, basereg)
#define rounding1noexec(base2, const1, xmmval, xmmcarry, xmmtmp, basereg) \
	rounding1exec(base2, const1, xmmval, xmmcarry, xmmtmp, 0)
#define rounding(ttp, base2, const1, sse4, xmmval, xmmcarry, xmmtmp, basereg) \
	rounding1##ttp(base2, const1, xmmval, xmmcarry, xmmtmp, basereg)


#define rounding_interleaved(ttp, base2, const1, sse4, xmmval, xmmcarry, xmmtmp, basereg, xmmval2, xmmcarry2, xmmtmp2, basereg2) \
	rounding(ttp, base2, const1, sse4, xmmval, xmmcarry, xmmtmp, basereg) \
	rounding(ttp, base2, const1, sse4, xmmval2, xmmcarry2, xmmtmp2, basereg2)


/* */ \
/* These macros round just one value in an XMM register.  This is done */ \
/* as part of the cleanup process where the final carry must be added */ \
/* back into the results. */ \
/* */ \
\

// base2_single_rounding
#define single_rounding2exec(xmmval, xmmcarry, xmmtmp, basereg) \
	xmmtmp = f64ptr2(g->u.xmm.XMM_LIMIT_BIGMAX, basereg); \
	xmmcarry = xmmval + xmmtmp; \
	xmmtmp = xmmcarry - xmmtmp; \
	xmmval -= xmmtmp; \
	xmmcarry *= f64ptr2(g->u.xmm.XMM_LIMIT_INVERSE, basereg);

// nobase2_single_rounding
#define single_rounding2noexec(xmmval, xmmcarry, xmmtmp, basereg) \
	xmmval -= XMM_BIGVAL1; \
	xmmcarry = f64ptr2(g->u.xmm.XMM_LIMIT_INVERSE, basereg); \
	xmmcarry = xmmcarry * xmmval + XMM_BIGVAL1; \
	xmmtmp = xmmcarry - XMM_BIGVAL1; \
	xmmtmp *= f64ptr2(g->u.xmm.XMM_LIMIT_BIGMAX, basereg); \
	xmmval -= xmmtmp;

#define single_rounding(base2, xmmval, xmmcarry, xmmtmp, basereg) \
	single_rounding2##base2(xmmval, xmmcarry, xmmtmp, basereg)

#define single_rounding_interleaved(base2, xmmval, xmmcarry, xmmtmp, basereg, xmmval2, xmmcarry2, xmmtmp2, basereg2) \
	single_rounding(base2, xmmval, xmmcarry, xmmtmp, basereg) \
	single_rounding(base2, xmmval2, xmmcarry2, xmmtmp2, basereg2)



/* */ \
/* These macros process zero-padded FFT result words.  These FFT results must */ \
/* be split into high and low parts with the high part used as a carry into */ \
/* the splitting the next FFT result word. */ \
/* */
#if 1
// base2_split_lower_zpad_word, base2_split_upper_zpad_word
#define split__zpad_wordexec(echk, sse4, xmmvalin, xmmcarry, xmmvalout, basereg, lower) \
lower(xmmvalin += xmmcarry);	/* Add in previous high FFT data */ \
	xmmcarry = g->u.xmm.XMM_BIGBIGVAL;	/* Big word rounding constant */ \
	xmmcarry += xmmvalin;		/* Round to multiple of big word */ \
	xmmcarry -= g->u.xmm.XMM_BIGBIGVAL; \
	xmmvalout = roundx(xmmvalin);	/* Round input to an integer */ \
echk(xmmvalin -= xmmvalout);	/* This is the convolution error */ \
echk(xmmvalin = absx(xmmvalin));	/* Compute absolute value */ \
echk(maxpd(xmm6, xmmvalin));		/* Compute maximum error */ \
	xmmvalout -= xmmcarry;	/* xmmvalout now contains low bigword bits */ \
lower(xmmcarry *= xptr2(g->u.xmm.XMM_LIMIT_INVERSE,basereg)); /* Saved shifted FFT hi data */ \

// nobase2_split_lower_zpad_word, nobase2_split_upper_zpad_word
#define split__zpad_wordnoexec(echk, sse4, xmmvalin, xmmcarry, xmmvalout, basereg, lower) \
lower(xmmvalin += xmmcarry);	/* Add in previous high FFT data */ \
	xmmcarry = xptr2(g->u.xmm.XMM_LIMIT_INVERSE,basereg); /* Load base */ \
	xmmcarry *= xmmvalin;	/* Compute FFTvalue / base */ \
	xmmcarry = roundx(xmmcarry);	/* Round to integer */ \
	xmmvalout = roundx(xmmvalin);	/* Round input to an integer */ \
echk(xmmvalin -= xmmvalout);	/* This is the convolution error */ \
echk(xmmvalin = absx(xmmvalin));	/* Compute absolute value */ \
echk(maxpd(xmm6, xmmvalin));		/* Compute maximum error */ \
	xmmvalin = xptr2(g->u.xmm.XMM_LIMIT_BIGMAX,basereg); \
	xmmvalin *= xmmcarry; \
	xmmvalout -= xmmvalin;	/* xmmvalout now contains FFTvalue % base */ \

#else
// base2
#define split__zpad_wordexec(echk, sse4, xmmvalin, xmmcarry, xmmvalout, basereg, lower) \
lower(xmmvalin += xmmcarry); \
	xmmcarry = xmmvalin + g->u.xmm.XMM_MINUS_C_TIMES_MULCONST; \
	xmmcarry -= g->u.xmm.XMM_MINUS_C_TIMES_MULCONST; \
	xmmvalout = roundx(xmmvalin); \
echk(xmmvalin -= xmmvalout); \
echk(xmmvalin = absx(xmmvalin)); \
echk(maxpd(xmm6, xmmvalin)); \
	xmmvalout -= xmmcarry; \
lower(xmmcarry *= xptr2(g->u.xmm.XMM_LIMIT_INVERSE, basereg));

// nobase2
#define split__zpad_wordnoexec(echk, sse4, xmmvalin, xmmcarry, xmmvalout, basereg, lower) \
lower(xmmvalin += xmmcarry;) \
	xmmcarry = xmmvalin * xptr2(g->u.xmm.XMM_LIMIT_INVERSE, basereg);\
	xmmcarry = roundx(xmmcarry); \
	xmmvalout = roundx(xmmvalin); \
echk(xmmvalin -= xmmvalout); \
echk(xmmvalin = absx(xmmvalin)); \
echk(maxpd(xmm6, xmmvalin)); \
	xmmvalin = xmmcarry * xptr2(g->u.xmm.XMM_LIMIT_BIGMAX, basereg);\
	xmmvalout -= xmmvalin;
#endif
#define split_lower_zpad_word(echk, base2, sse4, xmmvalin, xmmcarry, xmmvalout, basereg) \
	split__zpad_word##base2(echk, sse4, xmmvalin, xmmcarry, xmmvalout, basereg, exec)

\
/* Interleaved version used in 64-bit implementation of zpad normalization. */ \
/* Relies on some registers being preloaded. */ \
\

#define split_lower_zpad_word_interleaved(echk, base2, sse4, xmmvalin, xmmcarry, xmmvalout, basereg, xmmvalin2, xmmcarry2, xmmvalout2, basereg2) \
split__zpad_word##base2(echk, sse4, xmmvalin, xmmcarry, xmmvalout, basereg, exec) \
split__zpad_word##base2(echk, sse4, xmmvalin2, xmmcarry2, xmmvalout2, basereg2, exec)

\
/* Split upper is like split lower except that previous carry is not added in and */ \
/* result carry is not shifted down. */ \
\

#define split_upper_zpad_word(echk, base2, sse4, xmmvalin, xmmcarry, xmmvalout, basereg) \
split__zpad_word##base2(echk, sse4, xmmvalin, xmmcarry, xmmvalout, basereg, noexec)

\
/* Interleaved version used in 64-bit implementation of zpad normalization. */ \
/* Relies on some registers being preloaded. */ \
\

#define split_upper_zpad_word_interleaved(echk, base2, sse4, xmmvalin, xmmcarry, xmmvalout, basereg, xmmvalin2, xmmcarry2, xmmvalout2, basereg2) \
split__zpad_word##base2(echk, sse4, xmmvalin, xmmcarry, xmmvalout, basereg, noexec) \
split__zpad_word##base2(echk, sse4, xmmvalin2, xmmcarry2, xmmvalout2, basereg2, noexec)

\
/* The single word version */ \
\

#define single_split_lower_zpad_word(base2, xmmval, xmmcarry, xmmtmp, basereg) \
base2(	base2_single_split_lower_zpad_word(xmmval, xmmcarry, xmmtmp, basereg)); \
no##base2(nobase2_single_split_lower_zpad_word(xmmval, xmmcarry, xmmtmp, basereg)); \
\
\

#define base2_single_split_lower_zpad_word(xmmval, xmmcarry, xmmtmp, basereg) \
	xmmval += xmmcarry;	/* Add in previous high FFT data */ \
	xmmcarry = g->u.xmm.XMM_BIGBIGVAL[0]; /* Big word rounding constant */ \
	xmmcarry += xmmval;	/* Round to multiple of big word */ \
	xmmcarry -= g->u.xmm.XMM_BIGBIGVAL[0]; \
	xmmval = roundf64(xmmval); \
	xmmval -= xmmcarry;	/* xmmval now contains low bigword bits */ \
	xmmcarry *= f64ptr2(g->u.xmm.XMM_LIMIT_INVERSE,basereg);/* Next carry = shifted FFT hi data */ \
\
\

#define nobase2_single_split_lower_zpad_word(xmmval, xmmcarry, xmmtmp, basereg) \
	xmmval += xmmcarry;	/* Add in previous high FFT data */ \
	xmmcarry = f64ptr2(g->u.xmm.XMM_LIMIT_INVERSE,basereg); \
	xmmcarry *= xmmval;	/* Compute FFTvalue / base */ \
	xmmcarry = roundf64(xmmcarry); /* Next carry = round ( FFTvalue / base ) */ \
	xmmtmp = f64ptr2(g->u.xmm.XMM_LIMIT_BIGMAX,basereg); \
	xmmtmp *= xmmcarry; \
	xmmval = roundf64(xmmval); /* Round input to an integer */ \
	xmmval -= xmmtmp;		/* xmmval now contains FFTvalue % base */ \
\
\
/* The version for splitting the high FFT carry.  The high carry input has already */ \
/* been rounded to an integer. */ \
\

#define split_carry_zpad_word(base2, xmmcarryin, xmmcarryout, xmmtmp, basereg) \
base2(	base2_split_carry_zpad_word(xmmcarryin, xmmcarryout, xmmtmp, basereg)); \
no##base2(nobase2_split_carry_zpad_word(xmmcarryin, xmmcarryout, xmmtmp, basereg)); \
\
\

#define base2_split_carry_zpad_word(xmmcarryin, xmmcarryout, xmmtmp, basereg) \
	xmmcarryout = g->u.xmm.XMM_BIGBIGVAL; /* Big word rounding constant */ \
	xmmcarryout += xmmcarryin;	/* Round to multiple of big word */ \
	xmmcarryout -= g->u.xmm.XMM_BIGBIGVAL; \
	xmmcarryin -= xmmcarryout;	/* xmmcarryin now contains low bigword bits */ \
	xmmcarryout *= xptr2(g->u.xmm.XMM_LIMIT_INVERSE,basereg);/* Next carry = shifted carry */ \
\
\

#define nobase2_split_carry_zpad_word(xmmcarryin, xmmcarryout, xmmtmp, basereg) \
	xmmcarryout = xptr2(g->u.xmm.XMM_LIMIT_INVERSE,basereg); /* Load 1 / base */ \
	xmmcarryout *= xmmcarryin;	/* Compute carry / base */ \
	xmmcarryout = roundx(xmmcarryout); /* Next carry = round(carry / base) */ \
	xmmtmp = xptr2(g->u.xmm.XMM_LIMIT_BIGMAX,basereg); \
	xmmtmp *= xmmcarryout; \
	xmmcarryin -= xmmtmp;	/* xmmcarryin now contains carry % base */ \
\
\
/* The single float version for splitting the high FFT carry. */ \
\

#define single_split_carry_zpad_word(base2, xmmcarryin, xmmcarryout, xmmtmp, basereg) \
base2(	base2_single_split_carry_zpad_word(xmmcarryin, xmmcarryout, xmmtmp, basereg)); \
no##base2(nobase2_single_split_carry_zpad_word(xmmcarryin, xmmcarryout, xmmtmp, basereg)); \
\
\

#define base2_single_split_carry_zpad_word(xmmcarryin, xmmcarryout, xmmtmp, basereg) \
	xmmcarryout =  g->u.xmm.XMM_BIGBIGVAL[0];	/* Big word rounding constant */ \
	xmmcarryout += xmmcarryin;		/* Round to multiple of big word */ \
	xmmcarryout -= g->u.xmm.XMM_BIGBIGVAL[0]; \
	xmmcarryin -= xmmcarryout;		/* xmmcarryin now contains low bigword bits */ \
	xmmcarryout *= f64ptr2(g->u.xmm.XMM_LIMIT_INVERSE,basereg);/* Next carry = shifted carry */ \
\
\

#define nobase2_single_split_carry_zpad_word(xmmcarryin, xmmcarryout, xmmtmp, basereg) \
	xmmcarryout = f64ptr2(g->u.xmm.XMM_LIMIT_INVERSE,basereg); /* Load 1 / base */ \
	xmmcarryout *= xmmcarryin;		/* Compute carry / base */ \
	xmmcarryout = roundf64(xmmcarryout); /* Next carry = round(carry / base) */ \
	xmmtmp =  f64ptr2(g->u.xmm.XMM_LIMIT_BIGMAX,basereg); \
	xmmtmp *= xmmcarryout; \
	xmmcarryin -= xmmtmp;		/* xmmcarryin now contains carry % base */ \
\
\
/* Round the ZPAD0 - ZPAD6 values.  Simpler than other rounding macros */ \
/* in that we always round to a big word (and input value and output */ \
/* carry do not have XMM_BIGVAL added in). */ \
\

#define round_zpad7_word(base2, xmmvalin, xmmcarry, xmmvalout, basereg) \
base2(	base2_round_zpad7_word(xmmvalin, xmmcarry, xmmvalout, basereg)); \
no##base2(nobase2_round_zpad7_word(xmmvalin, xmmcarry, xmmvalout, basereg)); \
\
\

#define base2_round_zpad7_word(xmmval, xmmcarry, xmmtmp, basereg) \
	xmmval += xmmcarry;			/* Add in high part of last calculation */ \
	xmmcarry = g->u.xmm.XMM_BIGBIGVAL[0];		/* Big word rounding constant */ \
	xmmcarry += xmmval;			/* Round to multiple of big word */ \
	xmmcarry -= g->u.xmm.XMM_BIGBIGVAL[0]; \
	xmmval -= xmmcarry;			/* xmmval now contains low bigword bits */ \
	xmmcarry *= f64ptr2(g->u.xmm.XMM_LIMIT_INVERSE,basereg);	/* Shift high ZPAD data */ \
\
\

#define nobase2_round_zpad7_word(xmmval, xmmcarry, xmmtmp, basereg) \
	xmmval += xmmcarry;			/* Add in high part of last calculation */ \
	xmmcarry = f64ptr2(g->u.xmm.XMM_LIMIT_INVERSE,basereg);	/* Load base inverse */ \
	xmmcarry *= xmmval;			/* val / base */ \
	xmmcarry = roundf64(xmmcarry); /* next carry = round (val / base) */ \
	xmmtmp = f64ptr2(g->u.xmm.XMM_LIMIT_BIGMAX,basereg);	/* Load base */ \
	xmmtmp *= xmmcarry;			/* round (val / base) */ \
	xmmval -= xmmtmp;				/* new value = val - z */ \
\
\
\
/* */ \
/* Now for the actual normalization macros! */ \
/* */ \
\
\
/* For 1D macros, these registers are set on input: */ \
/* xmm6 = MAXERR */ \
/* xmm3 = carry #2 */ \
/* xmm2 = carry #1 */ \
/* rsi = pointer to the FFT data values */ \
/* rbp = pointer two-to-phi multipliers */ \
/* rdi = pointer to array of big vs. little flags */ \
/* ecx = big vs. little word flag #2 */ \
/* eax = big vs. little word flag #1 */ \
\
\
/* *************** 1D macro ****************** */ \
/* A pipelined version of this code: */ \
/*	movzx	rax, BYTE PTR [rdi]	;; Load big vs. little flags */ \
/*	xload	xmm0, [rsi+0*dist1]	;; Load values */ \
/*	mulpd	xmm0, [rbp+0]		;; Mul values1 by two-to-minus-phi */ \
/*	addpd	xmm0, xmm4		;; x = values + carry */ \
/*	xload	xmm2, XMM_LIMIT_BIGMAX[rax];; Load maximum * BIGVAL - BIGVAL */ \
/*	addpd	xmm2, xmm0		;; y = top bits of x */ \
/*	xload	xmm6, XMM_LIMIT_BIGMAX_NEG[rax];; Load -(maximum*BIGVAL-BIGVAL) */ \
/*	addpd	xmm6, xmm2		;; z = y - (maximum * BIGVAL - BIGVAL) */ \
/*	subpd	xmm0, xmm6		;; rounded value = x - z */ \
/*	mulpd	xmm2, XMM_LIMIT_INVERSE[rax];; next carry = shifted y */ \
/*	mulpd	xmm0, [rbp+16]		;; new value = val * two-to-phi */ \
/*	xstore	[rsi+0*dist1], xmm0	;; Save new value */ \
\

#define xnorm_1d(ttp, echk, const1, base2, sse4) { \
ttp(uintptr_t ax1 = u8ptr(rdi));	/* Load big vs. little flags */ \
ttp(uintptr_t cx1 = u8ptr(rdi+2)); \
	xmm0 = xptr(rsi);		/* Load values1 */ \
	unpcklo(xmm0, xptr(rsi+16)); \
	xmm0 *= xptr(rbp);		/* Mul values1 by two-to-minus-phi */ \
	xmm1 = xptr(rsi+32);		/* Load values2 */ \
	unpcklo(xmm1, xptr(rsi+48)); \
	xmm1 *= xptr(rbp+64);		/* Mul values2 by two-to-minus-phi */ \
no##const1(echk(error_check_interleaved(sse4, xmm0, xmm4, xmm1, xmm5, xmm6))); \
const1(mul_by_const(ttp, base2, echk, sse4, xmm0, xmm4, xmm5, ax1, xmm6)); \
	xmm2 += xmm0;		/* x1 = values + carry */ \
const1(mul_by_const(ttp, base2, echk, sse4, xmm1, xmm5, xmm0, cx1, xmm6)); \
	xmm3 += xmm1;		/* x2 = values + carry */ \
	rounding_interleaved(ttp, base2, const1, sse4, xmm2, xmm4, xmm0, ax1, xmm3, xmm5, xmm1, cx1); \
ttp(xmm2 *= xptr(rbp+16));		/* new value1 = val * two-to-phi */ \
ttp(xmm3 *= xptr(rbp+80));		/* new value2 = val * two-to-phi */ \
ttp(ax1 = u8ptr(rdi+1));	/* Load big vs. little flags */ \
ttp(cx1 = u8ptr(rdi+3)); \
	xmm0 = xptr(rsi);		/* Load values1 */ \
	unpckhi(xmm0, xptr(rsi+16)); \
	xmm0 *= xptr(rbp+32);		/* Mul values1 by two-to-minus-phi */ \
	xptr(rsi) = xmm2;		/* Save previous value1 */ \
	xmm1 = xptr(rsi+32);		/* Load values2 */ \
	unpckhi(xmm1, xptr(rsi+48)); \
	xmm1 *= xptr(rbp+96);		/* Mul values2 by two-to-minus-phi */ \
	xptr(rsi+32) = xmm3;		/* Save previous value2 */ \
no##const1(echk(error_check_interleaved(sse4, xmm0, xmm2, xmm1, xmm3, xmm6))); \
const1(mul_by_const(ttp, base2, echk, sse4, xmm0, xmm2, xmm3, ax1, xmm6)); \
	xmm0 += xmm4;		/* x1 = values + carry */ \
const1(mul_by_const(ttp, base2, echk, sse4, xmm1, xmm3, xmm4, cx1, xmm6)); \
	xmm1 += xmm5;		/* x2 = values + carry */ \
	rounding_interleaved(ttp, base2, const1, sse4, xmm0, xmm2, xmm4, ax1, xmm1, xmm3, xmm5, cx1); \
ttp(xmm0 *= xptr(rbp+48));		/* new value1 = val * two-to-phi */ \
ttp(xmm1 *= xptr(rbp+112));		/* new value2 = val * two-to-phi */ \
	xptr(rsi+16) = xmm0;		/* Save new value1 */ \
	xptr(rsi+48) = xmm1;		/* Save new value2 */ \
}

/* This is the normalization routine when we are computing modulo k*b^n+c */ \
/* with a zero-padded b^2n FFT.  We do this by multiplying the lower FFT */ \
/* word by k and adding in the upper word times -c.  Of course, this is made */ \
/* very tedious because we have to carefully avoid any loss of precision. */ \
/* */ \
/* xmm6 = MAXERR */ \
/* xmm3 = carry #2 (previous high FFT data - not yet mul'ed by K) */ \
/* xmm2 = carry #1 (traditional carry) */ \
/* rsi = pointer to the FFT data values */ \
/* rbp = pointer two-to-phi multipliers */ \
/* rdi = pointer to array of big vs. little flags */ \
/* eax = big vs. little word flag #1 */ \
\

#define xnorm_1d_zpad(ttp, echk, const1, base2, sse4, khi, c1, cm1) { \
no##ttp(const uintptr_t ax1=0); \
ttp(uintptr_t ax1 = u8ptr(rdi));	/* Load big vs. little flags */ \
	xmm0 = xptr(rsi);		/* Load values1 */ \
	unpcklo(xmm0, xptr(rsi+16)); \
	xmm1 = xptr(rsi+32);		/* Load values2 */ \
	unpcklo(xmm1, xptr(rsi+48)); \
	xmm0 *= xptr(rbp);		/* Mul values1 by two-to-minus-phi */ \
	xmm1 *= xptr(rbp);		/* Mul values2 by two-to-minus-phi */ \
\
	split_lower_zpad_word(echk, base2, sse4, xmm0, xmm3, xmm4, ax1); \
\
no##const1(xmm0 = g->u.xmm.XMM_K_LO); \
const1(	xmm0 = g->u.xmm.XMM_K_TIMES_MULCONST_LO); \
		xmm0 *= xmm4; \
khi(no##const1(xmm5 = g->u.xmm.XMM_K_HI)); \
khi(const1(xmm5 = g->u.xmm.XMM_K_TIMES_MULCONST_HI)); \
khi(no##base2(xmm5 *= xptr2(g->u.xmm.XMM_LIMIT_INVERSE,ax1))); /* Non-base2 rounding needs shifted carry */ \
khi(	xmm5 *= xmm4); \
\
		xmm0 += xmm2;		/* x1 = values + carry */ \
\
c1(	xmm1 *= g->u.xmm.XMM_MINUS_C);	/* Do one mul before split rather than two after split */ \
\
	split_upper_zpad_word(echk, base2, sse4, xmm1, xmm4, xmm2, ax1); \
\
no##const1(no##c1(no##cm1(xmm2 *= g->u.xmm.XMM_MINUS_C))); \
no##const1(no##c1(no##cm1(xmm4 *= g->u.xmm.XMM_MINUS_C))); \
const1(		xmm2 *= g->u.xmm.XMM_MINUS_C_TIMES_MULCONST); \
const1(		xmm4 *= g->u.xmm.XMM_MINUS_C_TIMES_MULCONST); \
\
	xmm0 += xmm2;		/* Add upper FFT word to lower FFT word */ \
khi(xmm4 += xmm5);		/* Add upper FFT word to lower FFT word */ \
\
 	rounding(ttp, base2, exec, sse4, xmm0, xmm4, xmm2, ax1); \
\
ttp(xmm0 *= xptr(rbp+16));		/* new value1 = val * two-to-phi */ \
\
ttp(ax1 = u8ptr(rdi+1));	/* Load big vs. little flags */ \
	xmm5 = xptr(rsi);		/* Load values1 */ \
	unpckhi(xmm5, xptr(rsi+16)); \
	xptr(rsi) = xmm0;		/* Save previous value1 */ \
	xmm1 = xptr(rsi+32);		/* Load values2 */ \
	unpckhi(xmm1, xptr(rsi+48)); \
	xmm5 *= xptr(rbp+32);		/* Mul values1 by two-to-minus-phi */ \
	xmm1 *= xptr(rbp+32);		/* Mul values2 by two-to-minus-phi */ \
\
	split_lower_zpad_word(echk, base2, sse4, xmm5, xmm3, xmm2, ax1); \
\
no##const1(xmm0 = g->u.xmm.XMM_K_LO); \
const1(	xmm0 = g->u.xmm.XMM_K_TIMES_MULCONST_LO); \
		xmm0 *= xmm2; \
khi(no##const1(xmm5 = g->u.xmm.XMM_K_HI)); \
khi(const1(xmm5 = g->u.xmm.XMM_K_TIMES_MULCONST_HI)); \
khi(no##base2(xmm5 *= xptr2(g->u.xmm.XMM_LIMIT_INVERSE,ax1))); /* Non-base2 rounding needs shifted carry */ \
khi(	xmm5 *= xmm2); \
\
		xmm0 += xmm4;		/* x1 = values + carry */ \
\
c1(	xmm1 *= g->u.xmm.XMM_MINUS_C);	/* Do one mul before split rather than two after split */ \
\
	split_upper_zpad_word(echk, base2, sse4, xmm1, xmm2, xmm4, ax1); \
\
no##const1(no##c1(no##cm1(xmm4 *= g->u.xmm.XMM_MINUS_C))); \
no##const1(no##c1(no##cm1(xmm2 *= g->u.xmm.XMM_MINUS_C))); \
const1(		xmm4 *= g->u.xmm.XMM_MINUS_C_TIMES_MULCONST); \
const1(		xmm2 *= g->u.xmm.XMM_MINUS_C_TIMES_MULCONST); \
\
	xmm0 += xmm4;		/* Add upper FFT word to lower FFT word */ \
khi(xmm2 += xmm5);		/* Add upper FFT word to lower FFT word */ \
\
	rounding(ttp, base2, exec, sse4, xmm0, xmm2, xmm4, ax1); \
\
ttp(xmm0 *= xptr(rbp+48));		/* new value1 = val * two-to-phi */ \
\
	xptr(rsi+16) = xmm0;		/* Save new value1 */ \
	xmm1[1] = xmm1[0] = 0.0;		/* new value2 = zero */ \
	xptr(rsi+32) = xmm1;		/* Zero previous value2 */ \
	xptr(rsi+48) = xmm1;		/* Zero current value2 */ \
}


/* *************** 1D followup macros ****************** */ \
/* This macro finishes the normalize process by adding the final */ \
/* carry from the first pass back into the lower two data values. */ \
/* xmm2,xmm3 = carries */ \
/* rsi = pointer to the FFT data values */ \
/* rbp = pointer two-to-phi multipliers */ \
/* rdi = big vs. little array pointer */ \
/* These registers are destroyed! */ \


#define xnorm012_1d_mid(ttp, base2){ \
ttp(uintptr_t ax1, cx1); \
no##ttp(const uintptr_t ax1=0, cx1=0); \
;	/* WARNING: Carry propagation into the 5th and 6th words does not work for */ \
;	/* FFT lengths 80 and 112.  gwnum.c makes sure that 5th carry will be */ \
;	/* zero for these FFT lengths.  The hack below adjusts the input pointers */ \
;	/* so that we only propagate 4 carry words */ \
if((g->FFTLEN & (64 | 16)) == (64 | 16)) { /* Length 80 and 112 have different */ \
						/* memory addresses for the 5th and 6th words */ \
	ttp(rdi += -4);			/* Modify input pointers so we */ \
		rsi += -64;		/* propagate into only 4 words */ \
	ttp(rbp += -128); \
}else{ /*not_80_or_112: */ \
	ttp(ax1 = u8ptr(rdi));	/* Load big vs. little flags */ \
	ttp(cx1 = u8ptr(rdi+2)); \
		xmm0[0] = f64ptr(rsi+8);		/* Load values1 */ \
	ttp(xmm0[0] *= f64ptr(rbp+8));		/* Mul values1 by two-to-minus-phi */ \
	ttp(xmm0[0] *= g->u.xmm.XMM_NORM012_FF[0]);	/* Mul by FFTLEN/2 */ \
		xmm1[0] = f64ptr(rsi+40);	/* Load values2 */ \
	ttp(xmm1[0] *= f64ptr(rbp+72));	/* Mul values2 by two-to-minus-phi */ \
	ttp(xmm1[0] *= g->u.xmm.XMM_NORM012_FF[0]);	/* Mul by FFTLEN/2 */ \
		xmm0[0] += xmm2[0];		/* x1 = values + carry */ \
		xmm1[0] += xmm3[0];		/* x2 = values + carry */ \
		single_rounding_interleaved(base2, xmm0[0], xmm2[0], xmm4[0], ax1+8, xmm1[0], xmm3[0], xmm5[0], cx1+8); \
	ttp(xmm0[0] *= f64ptr(rbp+24));	/* new value1 = val * two-to-phi */ \
	ttp(xmm1[0] *= f64ptr(rbp+88));	/* new value2 = val * two-to-phi */ \
		f64ptr(rsi+8) = xmm0[0];		/* Save new value1 */ \
		f64ptr(rsi+40) = xmm1[0];	/* Save new value2 */ \
	\
	ttp(ax1 = u8ptr(rdi+1));	/* Load big vs. little flags */ \
	ttp(cx1 = u8ptr(rdi+3)); \
		xmm0[0] = f64ptr(rsi+24);	/* Load values1 */ \
	ttp(xmm0[0] *= f64ptr(rbp+40));	/* Mul values1 by two-to-minus-phi */ \
	ttp(xmm0[0] *= g->u.xmm.XMM_NORM012_FF[0]);	/* Mul by FFTLEN/2 */ \
		xmm1[0] = f64ptr(rsi+56);	/* Load values2 */ \
	ttp(xmm1[0] *= f64ptr(rbp+104));	/* Mul values2 by two-to-minus-phi */ \
	ttp(xmm1[0] *= g->u.xmm.XMM_NORM012_FF[0]);	/* Mul by FFTLEN/2 */ \
		xmm0[0] += xmm2[0];		/* x1 = values + carry */ \
		xmm1[0] += xmm3[0];		/* x2 = values + carry */ \
		single_rounding_interleaved(base2, xmm0[0], xmm2[0], xmm4[0], ax1+8, xmm1[0], xmm3[0], xmm5[0], cx1+8); \
	ttp(xmm0[0] *= f64ptr(rbp+56));	/* new value1 = val * two-to-phi */ \
	ttp(xmm1[0] *= f64ptr(rbp+120));	/* new value2 = val * two-to-phi */ \
		f64ptr(rsi+24) = xmm0[0];	/* Save new value1 */ \
		f64ptr(rsi+56) = xmm1[0];	/* Save new value2 */ \
} \
/* only_do_4: */ \
ttp(ax1 = u8ptr(rdi+4));	/* Load big vs. little flags */ \
ttp(cx1 = u8ptr(rdi+6)); \
	xmm0[0] = f64ptr(rsi+64+8);	/* Load values1 */ \
ttp(xmm0[0] *= f64ptr(rbp+128+8));	/* Mul values1 by two-to-minus-phi */ \
ttp(xmm0[0] *= g->u.xmm.XMM_NORM012_FF[0]);	/* Mul by FFTLEN/2 */ \
	xmm1[0] = f64ptr(rsi+64+40);	/* Load values2 */ \
ttp(xmm1[0] *= f64ptr(rbp+128+72));	/* Mul values2 by two-to-minus-phi */ \
ttp(xmm1[0] *= g->u.xmm.XMM_NORM012_FF[0]);	/* Mul by FFTLEN/2 */ \
	xmm0[0] += xmm2[0];		/* x1 = values + carry */ \
	xmm1[0] += xmm3[0];		/* x2 = values + carry */ \
	single_rounding_interleaved(base2, xmm0[0], xmm2[0], xmm4[0], ax1+8, xmm1[0], xmm3[0], xmm5[0], cx1+8); \
ttp(xmm0[0] *= f64ptr(rbp+128+24));	/* new value1 = val * two-to-phi */ \
ttp(xmm1[0] *= f64ptr(rbp+128+88));	/* new value2 = val * two-to-phi */ \
	f64ptr(rsi+64+8) = xmm0[0];	/* Save new value1 */ \
	f64ptr(rsi+64+40) = xmm1[0];	/* Save new value2 */ \
\
ttp(ax1 = u8ptr(rdi+5));	/* Load big vs. little flags */ \
ttp(cx1 = u8ptr(rdi+7)); \
	xmm0[0] = f64ptr(rsi+64+24);	/* Load values1 */ \
ttp(xmm0[0] *= f64ptr(rbp+128+40));	/* Mul values1 by two-to-minus-phi */ \
ttp(xmm0[0] *= g->u.xmm.XMM_NORM012_FF[0]);	/* Mul by FFTLEN/2 */ \
	xmm1[0] = f64ptr(rsi+64+56);	/* Load values2 */ \
ttp(xmm1[0] *= f64ptr(rbp+128+104));	/* Mul values2 by two-to-minus-phi */ \
ttp(xmm1[0] *= g->u.xmm.XMM_NORM012_FF[0]);	/* Mul by FFTLEN/2 */ \
	xmm0[0] += xmm2[0];		/* x1 = values + carry */ \
	xmm1[0] += xmm3[0];		/* x2 = values + carry */ \
	single_rounding_interleaved(base2, xmm0[0], xmm2[0], xmm4[0], ax1+8, xmm1[0], xmm3[0], xmm5[0], cx1+8); \
ttp(xmm0[0] *= f64ptr(rbp+128+56));	/* new value1 = val * two-to-phi */ \
ttp(xmm1[0] *= f64ptr(rbp+128+120));	/* new value2 = val * two-to-phi */ \
	f64ptr(rsi+64+24) = xmm0[0];	/* Save new value1 */ \
	f64ptr(rsi+64+56) = xmm1[0];	/* Save new value2 */ \
\
ttp(ax1 = u8ptr(rdi+8));	/* Load big vs. little flags */ \
ttp(cx1 = u8ptr(rdi+10)); \
	xmm0[0] = f64ptr(rsi+128+8);	/* Load values1 */ \
ttp(xmm0[0] *= f64ptr(rbp+256+8));	/* Mul values1 by two-to-minus-phi */ \
ttp(xmm0[0] *= g->u.xmm.XMM_NORM012_FF[0]);	/* Mul by FFTLEN/2 */ \
	xmm1[0] = f64ptr(rsi+128+40);	/* Load values2 */ \
ttp(xmm1[0] *= f64ptr(rbp+256+72));	/* Mul values2 by two-to-minus-phi */ \
ttp(xmm1[0] *= g->u.xmm.XMM_NORM012_FF[0]);	/* Mul by FFTLEN/2 */ \
	xmm0[0] += xmm2[0];		/* x1 = values + carry */ \
	xmm1[0] += xmm3[0];		/* x2 = values + carry */ \
	single_rounding_interleaved(base2, xmm0[0], xmm2[0], xmm4[0], ax1+8, xmm1[0], xmm3[0], xmm5[0], cx1+8); \
ttp(xmm0[0] *= f64ptr(rbp+256+24));	/* new value1 = val * two-to-phi */ \
ttp(xmm1[0] *= f64ptr(rbp+256+88));	/* new value2 = val * two-to-phi */ \
	f64ptr(rsi+128+8) = xmm0[0];	/* Save new value1 */ \
	f64ptr(rsi+128+40) = xmm1[0];	/* Save new value2 */ \
\
	xmm2[0] -= XMM_BIGVAL1;/* Remove integer rounding constant */ \
	xmm3[0] -= XMM_BIGVAL1;/* Remove integer rounding constant */ \
ttp(xmm2[0] *= f64ptr(rbp+256+56));	/* carry *= two-to-phi */ \
ttp(xmm3[0] *= f64ptr(rbp+256+120));	/* carry *= two-to-phi */ \
	xmm2[0] += f64ptr(rsi+128+24);	/* value1 = values + carry */ \
	xmm3[0] += f64ptr(rsi+128+56);	/* value2 = values + carry */ \
	f64ptr(rsi+128+24) = xmm2[0];	/* Save new value1 */ \
	f64ptr(rsi+128+56) = xmm3[0];	/* Save new value2 */ \
\
	xmm2[0] = xmm2[1]; xmm2[1] = XMM_BIGVAL1; /* Rotate carry */ \
	xmm3[0] = xmm3[1]; xmm3[1] = XMM_BIGVAL1; /* Rotate carry */ \
}

/* This macro is similar to the above, but is for the zero padding case. */ \
/* xmm2 = carry #1 (traditional carry) */ \
/* xmm3 = carry #2 (previous high FFT data - not yet mul'ed by K) */ \
/* rbx = pointer to the FFT data values */ \
/* rdx = pointer two-to-phi multipliers */ \
/* rdi = big vs. little array pointer */ \

#define xnorm012_1d_mid_zpad(const1, base2){ \
	uintptr_t ax1; \
ax1 = u8ptr(rdi);	/* Load big vs. little flags */ \
	xmm0[0] = f64ptr(rbx+8);		/* Load values1 */ \
	xmm0[0] *= f64ptr(rdx+8);		/* Mul values1 by two-to-minus-phi */ \
	xmm0[0] *= g->u.xmm.XMM_NORM012_FF[0];	/* Mul by FFTLEN/2 */ \
	xmm0[0] += xmm2[0];		/* x1 = values1 + carry */ \
	single_split_carry_zpad_word(base2, xmm3[0], xmm1[0], xmm2[0], ax1+8); \
no##const1(xmm2[0] = g->u.xmm.XMM_K_LO[0]);	/* Calc high FFT carry times k */ \
const1(xmm2[0] = g->u.xmm.XMM_K_TIMES_MULCONST_LO[0]); \
	xmm2[0] *= xmm3[0];		/* high_FFT_carry * k_lo */ \
	xmm0[0] += xmm2[0];		/* x1 = x1 + high_FFT_carry * k_lo */ \
no##const1(xmm4[0] = g->u.xmm.XMM_K_HI[0]); \
const1(xmm4[0] = g->u.xmm.XMM_K_TIMES_MULCONST_HI[0]); \
	xmm4[0] *= f64ptr2(g->u.xmm.XMM_LIMIT_INVERSE,ax1+8); /* shift k_hi */ \
	xmm3[0] *= xmm4[0]; \
	single_rounding(base2, xmm0[0], xmm2[0], xmm4[0], ax1+8); \
	xmm2[0] += xmm3[0];		/* Carry += high_FFT_carry * shifted k_hi */ \
	xmm0[0] *= f64ptr(rdx+24);	/* new value1 = val * two-to-phi */ \
	f64ptr(rbx+8) = xmm0[0];		/* Save new value1 */ \
\
	ax1 = u8ptr(rdi+1);	/* Load big vs. little flags */ \
	xmm0[0] = f64ptr(rbx+24);	/* Load values2 */ \
	xmm0[0] *= f64ptr(rdx+40);	/* Mul values2 by two-to-minus-phi */ \
	xmm0[0] *= g->u.xmm.XMM_NORM012_FF[0];	/* Mul by FFTLEN/2 */ \
	xmm0[0] += xmm2[0];		/* x2 = values + carry */ \
	single_split_carry_zpad_word(base2, xmm1[0], xmm3[0], xmm2[0], ax1+8); \
no##const1(xmm2[0] = g->u.xmm.XMM_K_LO[0]);	/* Calc high FFT carry times k */ \
const1(xmm2[0] = g->u.xmm.XMM_K_TIMES_MULCONST_LO[0]); \
	xmm2[0] *= xmm1[0];		/* high_FFT_carry * k_lo */ \
	xmm0[0] += xmm2[0];		/* x2 = x2 + high_FFT_carry * k_lo */ \
no##const1(xmm4[0] = g->u.xmm.XMM_K_HI[0]); \
const1(xmm4[0] = g->u.xmm.XMM_K_TIMES_MULCONST_HI[0]); \
	xmm4[0] *= f64ptr2(g->u.xmm.XMM_LIMIT_INVERSE,ax1+8); /* shift k_hi */ \
	xmm1[0] *= xmm4[0]; \
	single_rounding(base2, xmm0[0], xmm2[0], xmm4[0], ax1+8); \
	xmm2[0] += xmm1[0];		/* Carry += high_FFT_carry * shifted k_hi */ \
	xmm0[0] *= f64ptr(rdx+56);	/* new value2 = val * two-to-phi */ \
	f64ptr(rbx+24) = xmm0[0];	/* Save new value2 */ \
\
	ax1 = u8ptr(rdi+4);	/* Load big vs. little flags */ \
	xmm0[0] = f64ptr(rbx+64+8);	/* Load values3 */ \
	xmm0[0] *= f64ptr(rdx+128+8);	/* Mul values3 by two-to-minus-phi */ \
	xmm0[0] *= g->u.xmm.XMM_NORM012_FF[0];	/* Mul by FFTLEN/2 */ \
	xmm0[0] += xmm2[0];		/* x3 = values + carry */ \
	single_split_carry_zpad_word(base2, xmm3[0], xmm1[0], xmm2[0], ax1+8); \
no##const1(xmm2[0] = g->u.xmm.XMM_K_LO[0]);	/* Calc high FFT carry times k */ \
const1(xmm2[0] = g->u.xmm.XMM_K_TIMES_MULCONST_LO[0]); \
	xmm2[0] *= xmm3[0];		/* high_FFT_carry * k_lo */ \
	xmm0[0] += xmm2[0];		/* x3 = x3 + high_FFT_carry * k_lo */ \
no##const1(xmm4[0] = g->u.xmm.XMM_K_HI[0]); \
const1(xmm4[0] = g->u.xmm.XMM_K_TIMES_MULCONST_HI[0]); \
	xmm4[0] *= f64ptr2(g->u.xmm.XMM_LIMIT_INVERSE,ax1+8); /* shift k_hi */ \
	xmm3[0] *= xmm4[0]; \
	single_rounding(base2, xmm0[0], xmm2[0], xmm4[0], ax1+8); \
	xmm2[0] += xmm3[0];		/* Carry += high_FFT_carry * shifted k_hi */ \
	xmm0[0] *= f64ptr(rdx+128+24);	/* new value3 = val * two-to-phi */ \
	f64ptr(rbx+64+8) = xmm0[0];	/* Save new value3 */ \
\
	ax1 = u8ptr(rdi+5);	/* Load big vs. little flags */ \
	xmm0[0] = f64ptr(rbx+64+24);	/* Load values4 */ \
	xmm0[0] *= f64ptr(rdx+128+40);	/* Mul values4 by two-to-minus-phi */ \
	xmm0[0] *= g->u.xmm.XMM_NORM012_FF[0];	/* Mul by FFTLEN/2 */ \
	xmm0[0] += xmm2[0];		/* x4 = values + carry */ \
	single_split_carry_zpad_word(base2, xmm1[0], xmm3[0], xmm2[0], ax1+8); \
no##const1(xmm2[0] = g->u.xmm.XMM_K_LO[0]);	/* Calc high FFT carry times k */ \
const1(xmm2[0] = g->u.xmm.XMM_K_TIMES_MULCONST_LO[0]); \
	xmm2[0] *= xmm1[0];		/* high_FFT_carry * k_lo */ \
	xmm0[0] += xmm2[0];		/* x4 = x4 + high_FFT_carry * k_lo */ \
no##const1(xmm4[0] = g->u.xmm.XMM_K_HI[0]); \
const1(xmm4[0] = g->u.xmm.XMM_K_TIMES_MULCONST_HI[0]); \
	xmm4[0] *= f64ptr2(g->u.xmm.XMM_LIMIT_INVERSE,ax1+8); /* shift k_hi */ \
	xmm1[0] *= xmm4[0]; \
	single_rounding(base2, xmm0[0], xmm2[0], xmm4[0], ax1+8); \
	xmm2[0] += xmm1[0];		/* Carry += high_FFT_carry * shifted k_hi */ \
	xmm0[0] *= f64ptr(rdx+128+56);	/* new value4 = val * two-to-phi */ \
	f64ptr(rbx+64+24) = xmm0[0];	/* Save new value4 */ \
\
	if(rbx != (uintptr_t)g->DESTARG || (g->FFTLEN != 80 && g->FFTLEN != 112)){ \
		/* Only the first section is funny */ \
		/* Length 80 and 112 have different */ \
		/* memory addresses for the fourth */ \
		/* and higher data elements */ \
/*notfunny:;*/ \
		ax1 = u8ptr(rdi+8);	/* Load big vs. little flags */ \
		xmm0[0] = f64ptr(rbx+128+8);	/* Load values5 */ \
		xmm0[0] *= f64ptr(rdx+256+8);	/* Mul values1 by two-to-minus-phi */ \
		xmm0[0] *= g->u.xmm.XMM_NORM012_FF[0];	/* Mul by FFTLEN/2 */ \
		xmm0[0] += xmm2[0];		/* x5 = values + carry */ \
	no##const1(xmm2[0] = g->u.xmm.XMM_K_LO[0]);	/* Calc high FFT carry times k */ \
	const1(xmm2[0] = g->u.xmm.XMM_K_TIMES_MULCONST_LO[0]); \
		xmm2[0] *= xmm3[0];		/* high_FFT_carry * k_lo */ \
		xmm0[0] += xmm2[0];		/* x5 = x5 + high_FFT_carry * k_lo */ \
	no##const1(xmm4[0] = g->u.xmm.XMM_K_HI[0]); \
	const1(xmm4[0] = g->u.xmm.XMM_K_TIMES_MULCONST_HI[0]); \
		xmm4[0] *= f64ptr2(g->u.xmm.XMM_LIMIT_INVERSE,ax1+8); /* shift k_hi */ \
		xmm3[0] *= xmm4[0]; \
		single_rounding(base2, xmm0[0], xmm2[0], xmm4[0], ax1+8); \
		xmm2[0] += xmm3[0];		/* Carry += high_FFT_carry * shifted k_hi */ \
		xmm0[0] *= f64ptr(rdx+256+24);	/* new value5 = val * two-to-phi */ \
		f64ptr(rbx+128+8) = xmm0[0];	/* Save new value5 */ \
	\
		xmm2[0] -= XMM_BIGVAL1;/* Remove integer rounding constant */ \
		xmm2[0] *= f64ptr(rdx+256+56);	/* carry *= two-to-phi */ \
		xmm2[0] += f64ptr(rbx+128+24);	/* value6 = values + carry */ \
		f64ptr(rbx+128+24) = xmm2[0];	/* Save new value6 */ \
	}else{				/* FFT length = 80 or 112 */ \
		xmm0[0] = xmm2[0];		/* Add the carry to the next section's carry */ \
		xmm0[0] -= XMM_BIGVAL1; \
		/*shufpd	xmm2, xmm2, 1;*/ \
		xmm2[1] += xmm0[0]; \
		/*shufpd	xmm2, xmm2, 1;*/ \
		xmm0[0] = xmm3[0];		/* Add the carry to the next section's carry */ \
		/*shufpd	xmm3, xmm3, 1;*/ \
		xmm3[1] += xmm0[0]; \
		/*shufpd	xmm3, xmm3, 1;*/ \
	} \
	xmm2[0] = xmm2[1]; xmm2[1] = XMM_BIGVAL1; /* Rotate carry */ \
	xmm3[0] = xmm3[1]; xmm3[1] = 0.0; \
}


/* We could take advantage of the fact that the first two-to-phi multiplier */ \
/* and the first two-to-minus-phi multiplier are one.  We also know */ \
/* the first data value is a big word (eax would be 48). */ \
/* xmm2,xmm3 = carries */ \
/* rsi = pointer to the FFT data values */ \
/* rbp = pointer two-to-phi multipliers */ \
/* rdi = big vs. little array pointer */ \
/* NOTE: If RATIONAL_FFT we could eliminate 8 multiplies. */ \
/* Input arguments are destroyed! */ \

#define xnorm012_1d(base2) { \
	double f0, f1, f4, f5; \
	uintptr_t rax, rcx; \
	xmm3[0] -= XMM_BIGVAL1; \
	xmm3[0] *= g->u.xmm.XMM_MINUS_C[0];	/* Adjust wrap-around carry */ \
	xmm3[0] += XMM_BIGVAL1; \
\
;	/* WARNING: Carry propagation into the 5th and 6th words does not work for */ \
;	/* FFT lengths 80 and 112.  gwnum.c makes sure that 5th carry will be */ \
;	/* zero for these FFT lengths.  The hack below adjusts the input pointers */ \
;	/* so that we only propagate 4 carry words */ \
if((g->FFTLEN & (64 | 16)) == (64 | 16)) { /* Length 80 and 112 have different */ \
							/* memory addresses for the 5th and 6th words */ \
	rdi += -4;			/* Modify input pointers so we */ \
	rsi += -64;		/* propagate into only 4 words */ \
	rbp += -128; \
}else{ /* not_80_or_112: */ \
	rax = u8ptr(rdi);	/* Load big vs. little flags */ \
	rcx = u8ptr(rdi+2); \
	f0 = f64ptr(rsi);		/* Load values1 */ \
	f0 *= f64ptr(rbp);		/* Mul values1 by two-to-minus-phi */ \
	f0 *= g->u.xmm.XMM_NORM012_FF[0];	/* Mul by FFTLEN/2 */ \
	f1 = f64ptr(rsi+32);	/* Load values2 */ \
	f1 *= f64ptr(rbp+64);	/* Mul values2 by two-to-minus-phi */ \
	f1 *= g->u.xmm.XMM_NORM012_FF[0];	/* Mul by FFTLEN/2 */ \
	f0 += xmm3[0];		/* x1 = values + carry */ \
	f1 += xmm2[0];		/* x2 = values + carry */ \
	single_rounding_interleaved(base2, f0, xmm3[0], f4, rax, f1, xmm2[0], f5, rcx); \
	f0 *= f64ptr(rbp+16);	/* new value1 = val * two-to-phi */ \
	f1 *= f64ptr(rbp+80);	/* new value2 = val * two-to-phi */ \
	f64ptr(rsi) = f0;		/* Save value1 */ \
	f64ptr(rsi+32) = f1;	/* Save value2 */ \
\
	rax = u8ptr(rdi+1);	/* Load big vs. little flags */ \
	rcx = u8ptr(rdi+3); \
	f0 = f64ptr(rsi+16);	/* Load values1 */ \
	f0 *= f64ptr(rbp+32);	/* Mul values1 by two-to-minus-phi */ \
	f0 *= g->u.xmm.XMM_NORM012_FF[0];	/* Mul by FFTLEN/2 */ \
	f1 = f64ptr(rsi+48);	/* Load values2 */ \
	f1 *= f64ptr(rbp+96);	/* Mul values2 by two-to-minus-phi */ \
	f1 *= g->u.xmm.XMM_NORM012_FF[0];	/* Mul by FFTLEN/2 */ \
	f0 += xmm3[0];		/* x1 = values + carry */ \
	f1 += xmm2[0];		/* x2 = values + carry */ \
	single_rounding_interleaved(base2, f0, xmm3[0], f4, rax, f1, xmm2[0], f5, rcx); \
	f0 *= f64ptr(rbp+48);	/* new value1 = val * two-to-phi */ \
	f1 *= f64ptr(rbp+112);	/* new value2 = val * two-to-phi */ \
	f64ptr(rsi+16) = f0;	/* Save value1 */ \
	f64ptr(rsi+48) = f1;	/* Save value2 */ \
\
} /*only_do_4:*/ \
	rax = u8ptr(rdi+4);	/* Load big vs. little flags */ \
	rcx = u8ptr(rdi+6); \
	f0 = f64ptr(rsi+64);	/* Load values1 */ \
	f0 *= f64ptr(rbp+128);	/* Mul values1 by two-to-minus-phi */ \
	f0 *= g->u.xmm.XMM_NORM012_FF[0];	/* Mul by FFTLEN/2 */ \
	f1 = f64ptr(rsi+64+32);	/* Load values2 */ \
	f1 *= f64ptr(rbp+128+64);	/* Mul values2 by two-to-minus-phi */ \
	f1 *= g->u.xmm.XMM_NORM012_FF[0];	/* Mul by FFTLEN/2 */ \
	f0 += xmm3[0];		/* x1 = values + carry */ \
	f1 += xmm2[0];		/* x2 = values + carry */ \
	single_rounding_interleaved(base2, f0, xmm3[0], f4, rax, f1, xmm2[0], f5, rcx); \
	f0 *= f64ptr(rbp+128+16);	/* new value1 = val * two-to-phi */ \
	f1 *= f64ptr(rbp+128+80);	/* new value2 = val * two-to-phi */ \
	f64ptr(rsi+64) = f0;	/* Save value1 */ \
	f64ptr(rsi+64+32) = f1;	/* Save value2 */ \
\
	rax = u8ptr(rdi+5);	/* Load big vs. little flags */ \
	rcx = u8ptr(rdi+7); \
	f0 = f64ptr(rsi+64+16);	/* Load values1 */ \
	f0 *= f64ptr(rbp+128+32);	/* Mul values1 by two-to-minus-phi */ \
	f0 *= g->u.xmm.XMM_NORM012_FF[0];	/* Mul by FFTLEN/2 */ \
	f1 = f64ptr(rsi+64+48);	/* Load values2 */ \
	f1 *= f64ptr(rbp+128+96);	/* Mul values2 by two-to-minus-phi */ \
	f1 *= g->u.xmm.XMM_NORM012_FF[0];	/* Mul by FFTLEN/2 */ \
	f0 += xmm3[0];		/* x1 = values + carry */ \
	f1 += xmm2[0];		/* x2 = values + carry */ \
	single_rounding_interleaved(base2, f0, xmm3[0], f4, rax, f1, xmm2[0], f5, rcx); \
	f0 *= f64ptr(rbp+128+48);	/* new value1 = val * two-to-phi */ \
	f1 *= f64ptr(rbp+128+112);	/* new value2 = val * two-to-phi */ \
	f64ptr(rsi+64+16) = f0;	/* Save value1 */ \
	f64ptr(rsi+64+48) = f1;	/* Save value2 */ \
\
	rax = u8ptr(rdi+8);	/* Load big vs. little flags */ \
	rcx = u8ptr(rdi+10); \
	f0 = f64ptr(rsi+128);	/* Load values1 */ \
	f0 *= f64ptr(rbp+256);	/* Mul values1 by two-to-minus-phi */ \
	f0 *= g->u.xmm.XMM_NORM012_FF[0];	/* Mul by FFTLEN/2 */ \
	f1 = f64ptr(rsi+128+32);	/* Load values2 */ \
	f1 *= f64ptr(rbp+256+64);	/* Mul values2 by two-to-minus-phi */ \
	f1 *= g->u.xmm.XMM_NORM012_FF[0];	/* Mul by FFTLEN/2 */ \
	f0 += xmm3[0];		/* x1 = values + carry */ \
	f1 += xmm2[0];		/* x2 = values + carry */ \
	single_rounding_interleaved(base2, f0, xmm3[0], f4, rax, f1, xmm2[0], f5, rcx); \
	f0 *= f64ptr(rbp+256+16);	/* new value1 = val * two-to-phi */ \
	f1 *= f64ptr(rbp+256+80);	/* new value2 = val * two-to-phi */ \
	f64ptr(rsi+128) = f0;	/* Save value1 */ \
	f64ptr(rsi+128+32) = f1;	/* Save value2 */ \
\
	xmm3[0] -= XMM_BIGVAL1;/* Remove integer rounding constant */ \
	xmm2[0] -= XMM_BIGVAL1;/* Remove integer rounding constant */ \
	xmm3[0] *= f64ptr(rbp+256+48);	/* carry *= two-to-phi */ \
	xmm2[0] *= f64ptr(rbp+256+112);	/* carry *= two-to-phi */ \
	xmm3[0] += f64ptr(rsi+128+16);	/* value1 = values + carry */ \
	xmm2[0] += f64ptr(rsi+128+48);	/* value2 = values + carry */ \
	f64ptr(rsi+128+16) = xmm3[0];	/* Save value1 */ \
	f64ptr(rsi+128+48) = xmm2[0];	/* Save value2 */ \
}

/* This macro is similar to the above, but is for the zero padding case. */ \
/* xmm2 = carry #1 (traditional carry) */ \
/* xmm3 = carry #2 (previous high FFT data - not yet mul'ed by K) */ \
/* rsi = pointer to the FFT data values */ \
/* rbp = pointer two-to-phi multipliers */ \
/* rdi = big vs. little array pointer */ \
/* NOTE: If RATIONAL_FFT we could eliminate 8 multiplies. */ \
\

#define xnorm012_1d_zpad(const1, base2) { \
	double f0, f1, f4, f5; \
;	/* Strip BIGVAL from the traditional carry, we'll add the traditional */ \
;	/* carry in later when we are working on the ZPAD0 - ZPAD6 values. */ \
	xmm2[0] -= XMM_BIGVAL1;/* Integerize traditional carry */ \
\
;	/* Rather than calculate high FFT carry times k and then later dividing */ \
;	/* by k, we multiply FFT high carry by const and we'll add it */ \
;	/* to the lower FFT data later (after multiplying by -c). */ \
const1(xmm3[0] *= g->u.xmm.XMM_MULCONST[0]); \
\
;	/* Work on zero-pad addin value. */ \
	f4 = g->ADDIN_VALUE;		/* Load the add in value */ \
const1(f4 *= g->u.xmm.XMM_MULCONST[0]);		/* Multiply the add in value by the small mul const */ \
	f4 += g->POSTADDIN_VALUE;		/* Add the post-mul-by-const addin value */ \
;	/* when c = -1, 1 = b^n */ \
;	/* when c = 1, -1 = b^n, 1 = -b^n */ \
;	/* when c = -3, 3 = b^n, 1 = b^n - 2 */ \
;	/* when c = 3, -3 = b^n, 3 = -b^n, 1 = -b^n - 2 */ \
;	/* The "- 2" has been precomputed in ZPAD_LSW_ADJUST.  Add ADDIN_VALUE * mul-by-const * ZPAD_LSW_ADJUST into the least significant FFT word */ \
	rsi = (uintptr_t)g->DESTARG;/* Address of squared number */ \
	f0 = g->ZPAD_LSW_ADJUST; \
	f0 *= f4; \
	f0 += f64ptr(rsi); \
	f64ptr(rsi) = f0; \
\
;	/* Multiply ZPAD0 through ZPAD6 by const * -c.  This, in essense, */ \
;	/* wraps this data from above the FFT data area to the halfway point. */ \
;	/* Later on we'll divide this by K to decide which data needs wrapping */ \
;	/* all the way down to the bottom of the FFT data. */ \
\
;	/* NOTE that ZPAD0's column multiplier is 1.0.  Also, ZPAD6 will not */ \
;	/* be bigger than a big word.  We must be careful to handle c's up */ \
;	/* to about 30 bits */ \
\
	rax = u8ptr(rdi);	/* Load big vs. little flags */ \
	f0 = g->ZPAD0_6[0];/* Load values1 */ \
	f5 = 0;		/* Create a zero high FFT data carry */ \
	single_split_lower_zpad_word(base2, f0, f5, f1, rax); \
no##const1(f0 *= g->u.xmm.XMM_MINUS_C[0]); \
const1(f0 *= g->u.xmm.XMM_MINUS_C_TIMES_MULCONST[0]); \
	f0 += f4;		/* Apply rest of ADDIN_VALUE * mul-by-const here (the FFT half way point) */ \
	round_zpad7_word(base2, f0, xmm2[0], f1, rax); \
	g->ZPAD0_6[0] = f0; \
\
	rax = u8ptr(rdi+1);	/* Load big vs. little flags */ \
	f0 = g->ZPAD0_6[1];/* Load values1 */ \
	f0 *= f64ptr(rbp+32);	/* Mul values1 by two-to-minus-phi */ \
	f0 *= g->u.xmm.XMM_NORM012_FF[0];	/* Mul by FFTLEN/2 */ \
	single_split_lower_zpad_word(base2, f0, f5, f1, rax); \
no##const1(f0 *= g->u.xmm.XMM_MINUS_C[0]); \
const1(f0 *= g->u.xmm.XMM_MINUS_C_TIMES_MULCONST[0]); \
	round_zpad7_word(base2, f0, xmm2[0], f1, rax); \
	g->ZPAD0_6[1] = f0; \
\
	rax = u8ptr(rdi+4);	/* Load big vs. little flags */ \
	f0 = g->ZPAD0_6[2];/* Load values1 */ \
	f0 *= f64ptr(rbp+128);	/* Mul values1 by two-to-minus-phi */ \
	f0 *= g->u.xmm.XMM_NORM012_FF[0];	/* Mul by FFTLEN/2 */ \
	single_split_lower_zpad_word(base2, f0, f5, f1, rax); \
no##const1(f0 *= g->u.xmm.XMM_MINUS_C[0]); \
const1(f0 *= g->u.xmm.XMM_MINUS_C_TIMES_MULCONST[0]); \
	round_zpad7_word(base2, f0, xmm2[0], f1, rax); \
	g->ZPAD0_6[2] = f0; \
\
	rax = u8ptr(rdi+5);	/* Load big vs. little flags */ \
	f0 = g->ZPAD0_6[3];/* Load values1 */ \
	f0 *= f64ptr(rbp+128+32);	/* Mul values1 by two-to-minus-phi */ \
	f0 *= g->u.xmm.XMM_NORM012_FF[0];	/* Mul by FFTLEN/2 */ \
	single_split_lower_zpad_word(base2, f0, f5, f1, rax); \
no##const1(f0 *= g->u.xmm.XMM_MINUS_C[0]); \
const1(f0 *= g->u.xmm.XMM_MINUS_C_TIMES_MULCONST[0]); \
	round_zpad7_word(base2, f0, xmm2[0], f1, rax); \
	g->ZPAD0_6[3] = f0; \
\
	if(g->FFTLEN != 80 && g->FFTLEN != 112){ \
		/* Length 80 and 112 have different */ \
		/* memory addresses for the fourth */ \
		/* and higher data elements */ \
		rax = u8ptr(rdi+8);	/* Load big vs. little flags */ \
		f0 = g->ZPAD0_6[4];/* Load values1 */ \
		f0 *= f64ptr(rbp+256);	/* Mul values1 by two-to-minus-phi */ \
		f0 *= g->u.xmm.XMM_NORM012_FF[0];	/* Mul by FFTLEN/2 */ \
		single_split_lower_zpad_word(base2, f0, f5, f1, rax); \
	no##const1(f0 *= g->u.xmm.XMM_MINUS_C[0]); \
	const1(f0 *= g->u.xmm.XMM_MINUS_C_TIMES_MULCONST[0]); \
		round_zpad7_word(base2, f0, xmm2[0], f1, rax); \
		g->ZPAD0_6[4] = f0; \
	\
		rax = u8ptr(rdi+9);	/* Load big vs. little flags */ \
		f0 = g->ZPAD0_6[5];/* Load values1 */ \
		f0 *= f64ptr(rbp+256+32);	/* Mul values1 by two-to-minus-phi */ \
		f0 *= g->u.xmm.XMM_NORM012_FF[0];	/* Mul by FFTLEN/2 */ \
		single_split_lower_zpad_word(base2, f0, f5, f1, rax); \
	no##const1(f0 *= g->u.xmm.XMM_MINUS_C[0]); \
	const1(f0 *= g->u.xmm.XMM_MINUS_C_TIMES_MULCONST[0]); \
		round_zpad7_word(base2, f0, xmm2[0], f1, rax); \
		g->ZPAD0_6[5] = f0; \
	\
		f0 = f64ptr(rbp+384);	/* Load two-to-minus-phi */ \
	}else{	/* Join common code */ \
		/* Same as the above but with different addresses required by */ \
		/* the length 80 and 112 FFT lengths */ \
		rax = u8ptr(rdi);	/* Load big vs. little flags */ \
		f0 = g->ZPAD0_6[4];/* Load values1 */ \
		f0 *= f64ptr(rbp+8);		/* Mul values1 by two-to-minus-phi */ \
		f0 *= g->u.xmm.XMM_NORM012_FF[0];	/* Mul by FFTLEN/2 */ \
		single_split_lower_zpad_word(base2, f0, f5, f1, rax+8); \
	no##const1(f0 *= g->u.xmm.XMM_MINUS_C[0]); \
	const1(f0 *= g->u.xmm.XMM_MINUS_C_TIMES_MULCONST[0]); \
		round_zpad7_word(base2, f0, xmm2[0], f1, rax+8); \
		g->ZPAD0_6[4] = f0; \
	\
		rax = u8ptr(rdi+1);	/* Load big vs. little flags */ \
		f0 = g->ZPAD0_6[5];/* Load values1 */ \
		f0 *= f64ptr(rbp+32+8);	/* Mul values1 by two-to-minus-phi */ \
		f0 *= g->u.xmm.XMM_NORM012_FF[0];	/* Mul by FFTLEN/2 */ \
		single_split_lower_zpad_word(base2, f0, f5, f1, rax+8); \
	no##const1(f0 *= g->u.xmm.XMM_MINUS_C[0]); \
	const1(f0 *= g->u.xmm.XMM_MINUS_C_TIMES_MULCONST[0]); \
		round_zpad7_word(base2, f0, xmm2[0], f1, rax+8); \
		g->ZPAD0_6[5] = f0; \
	\
		f0 = f64ptr(rbp+128+8);	/* Load two-to-minus-phi */ \
	} \
	f0 *= g->ZPAD0_6[6];/* Mul by values1 */ \
	f0 *= g->u.xmm.XMM_NORM012_FF[0];	/* Mul by FFTLEN/2 */ \
	f0 = roundf64(f0); /* Round to an integer */ \
	f0 += f5;		/* Add in shifted high ZPAD data */ \
no##const1(f0 *= g->u.xmm.XMM_MINUS_C[0]); \
const1(f0 *= g->u.xmm.XMM_MINUS_C_TIMES_MULCONST[0]); \
	f0 += xmm2[0];		/* Add in high part of last calculation */ \
	g->ZPAD0_6[6] = f0; \
\
;	/* Divide the zpad data by k.  Store the integer part in XMM_TMP */ \
;	/* and the remainder in ZPAD0.  Later we will wrap the integer part */ \
;	/* down to the bottom of the FFT data area (and multiply by -c). */ \
;	/* And we will store the remainder in the upper half of the FFT */ \
;	/* data area. */ \
\
;	/* Note there are three cases to handle.  K is smaller than a big word. */ \
;	/* K is between one and 2 big words in size.  And K is more than */ \
;	/* 2 big words in size. */ \
\
	xfunc572(g); \
\
;	/* Now normalize the data above the halfway point.  Remember that the */ \
;	/* column two-to-phi multiplier for the first value will be 1.0. */ \
\
	rax = u8ptr(rdi);	/* First word  */ \
	f0 = g->ZPAD0_6[0];/* Load remainder of divide by k */ \
	f0 += XMM_BIGVAL1; \
	single_rounding(base2, f0, xmm2[0], f4, rax); \
	f64ptr(rsi+32) = f0;	/* Save value1 */ \
\
	rax = u8ptr(rdi+1);	/* Load big vs. little flags */ \
	single_rounding(base2, xmm2[0], f0, f4, rax); \
	xmm2[0] *= f64ptr(rbp+48);	/* new value2 = val * two-to-phi */ \
	f64ptr(rsi+48) = xmm2[0];	/* Save value2 */ \
\
	rax = u8ptr(rdi+4);	/* Load big vs. little flags */ \
	single_rounding(base2, f0, xmm2[0], f4, rax); \
	f0 *= f64ptr(rbp+128+16);	/* new value3 = val * two-to-phi */ \
	f64ptr(rsi+64+32) = f0;	/* Save value3 */ \
\
	xmm2[0] -= XMM_BIGVAL1;/* Remove integer rounding constant */ \
	xmm2[0] *= f64ptr(rbp+128+48);	/* value4 = carry * two-to-phi */ \
	f64ptr(rsi+64+48) = xmm2[0];	/* Save new value4 */ \
\
;	/* Mul the integer part of (ZPAD data divided by k) by -c in */ \
;	/* preparation for adding it into the lower FFT data area. */ \
;	/* Also add in the shifted high FFT carry at this time. */ \
\
;	/* Now add in and normalize the bottom FFT data.  Remember that the */ \
;	/* column two-to-phi multiplier for the first value will be 1.0.  We  */ \
;	/* must go 6 words deep in case k is 48-50 bits and c is 32 bits. */ \
\
	vec2f64 TMP12 = g->u.xmm.XMM_TMP1; \
	vec2f64 TMP34 = g->u.xmm.XMM_TMP2; \
	vec2f64 TMP56 = g->u.xmm.XMM_TMP3; \
	TMP34 *= g->u.xmm.XMM_MINUS_C;	/* Mul by -c */ \
	TMP56 *= g->u.xmm.XMM_MINUS_C;	/* Mul by -c */ \
	rax = u8ptr(rdi);	/* First word  */ \
	f0 = TMP12[0];	/* Load integer part of divide by k */ \
	f0 += xmm3[0];		/* Add in shifted high FFT carry */ \
	f0 *= g->u.xmm.XMM_MINUS_C[0];	/* Mul by -c */ \
	f0 += XMM_BIGVAL1; \
	f0 += f64ptr(rsi);		/* Add in the FFT data */ \
	single_rounding(base2, f0, xmm2[0], f4, rax); \
	f64ptr(rsi) = f0;		/* Save value1 */ \
\
	rax = u8ptr(rdi+1);	/* Load big vs. little flags */ \
	f0 = TMP12[1];	/* Load integer part of divide by k */ \
	f0 *= g->u.xmm.XMM_MINUS_C[0];	/* Mul by -c */ \
	f1 = f64ptr(rsi+16);	/* Load FFT data */ \
	f1 *= f64ptr(rbp+32);	/* Mul values2 by two-to-minus-phi */ \
	f1 *= g->u.xmm.XMM_NORM012_FF[0];	/* Mul by FFTLEN/2 */ \
	f0 += f1;		/* Add in the FFT data */ \
	f0 += xmm2[0];		/* x2 = value + carry */ \
	single_rounding(base2, f0, xmm2[0], f4, rax); \
	f0 *= f64ptr(rbp+48);	/* new value2 = val * two-to-phi */ \
	f64ptr(rsi+16) = f0;	/* Save value2 */ \
\
	rax = u8ptr(rdi+4);	/* Load big vs. little flags */ \
	f0 = TMP34[0];	/* Load integer part of divide by k */ \
	f1 = f64ptr(rsi+64);	/* Load FFT data */ \
	f1 *= f64ptr(rbp+128);	/* Mul values3 by two-to-minus-phi */ \
	f1 *= g->u.xmm.XMM_NORM012_FF[0];	/* Mul by FFTLEN/2 */ \
	f0 += f1;		/* Add in the FFT data */ \
	f0 += xmm2[0];		/* x3 = value + carry */ \
	single_rounding(base2, f0, xmm2[0], f4, rax); \
	f0 *= f64ptr(rbp+128+16);	/* new value3 = val * two-to-phi */ \
	f64ptr(rsi+64) = f0;	/* Save value3 */ \
\
	rax = u8ptr(rdi+5);	/* Load big vs. little flags */ \
	f0 = TMP34[1];	/* Load integer part of divide by k */ \
	f1 = f64ptr(rsi+64+16);	/* Load FFT data */ \
	f1 *= f64ptr(rbp+128+32);	/* Mul values4 by two-to-minus-phi */ \
	f1 *= g->u.xmm.XMM_NORM012_FF[0];	/* Mul by FFTLEN/2 */ \
	f0 += f1;		/* Add in the FFT data */ \
	f0 += xmm2[0];		/* x4 = value + carry */ \
	single_rounding(base2, f0, xmm2[0], f4, rax); \
	f0 *= f64ptr(rbp+128+48);	/* new value4 = val * two-to-phi */ \
	f64ptr(rsi+64+16) = f0;	/* Save value4 */ \
\
	if(g->FFTLEN != 80 && g->FFTLEN != 112){ \
		/* Length 80 and 112 have different */ \
		/* memory addresses for the fourth */ \
		/* and higher data elements */ \
		rax = u8ptr(rdi+8);	/* Load big vs. little flags */ \
		f0 = TMP56[0];	/* Load integer part of divide by k */ \
		f1 = f64ptr(rsi+128);	/* Load FFT data */ \
		f1 *= f64ptr(rbp+256);	/* Mul values4 by two-to-minus-phi */ \
		f1 *= g->u.xmm.XMM_NORM012_FF[0];	/* Mul by FFTLEN/2 */ \
		f0 += f1;		/* Add in the FFT data */ \
		f0 += xmm2[0];		/* x5 = value + carry */ \
		single_rounding(base2, f0, xmm2[0], f4, rax); \
		f0 *= f64ptr(rbp+256+16);	/* new value5 = val * two-to-phi */ \
		f64ptr(rsi+128) = f0;	/* Save value5 */ \
	\
		rax = u8ptr(rdi+9);	/* Load big vs. little flags */ \
		f0 = TMP56[1];	/* Load integer part of divide by k */ \
		f1 = f64ptr(rsi+128+16);	/* Load FFT data */ \
		f1 *= f64ptr(rbp+256+32);	/* Mul values6 by two-to-minus-phi */ \
		f1 *= g->u.xmm.XMM_NORM012_FF[0];	/* Mul by FFTLEN/2 */ \
		f0 += f1;		/* Add in the FFT data */ \
		f0 += xmm2[0];		/* x6 = value + carry */ \
		single_rounding(base2, f0, xmm2[0], f4, rax); \
		f0 *= f64ptr(rbp+256+48);	/* new value6 = val * two-to-phi */ \
		f64ptr(rsi+128+16) = f0;	/* Save value6 */ \
	\
		rax = u8ptr(rdi+12);	/* Load big vs. little flags */ \
		f0 = f64ptr(rsi+192);	/* Load FFT data */ \
		f0 *= f64ptr(rbp+384);	/* Mul values7 by two-to-minus-phi */ \
		f0 *= g->u.xmm.XMM_NORM012_FF[0];	/* Mul by FFTLEN/2 */ \
		f0 += xmm2[0];		/* x7 = value + carry */ \
		single_rounding(base2, f0, xmm2[0], f4, rax); \
		f0 *= f64ptr(rbp+384+16);	/* new value7 = val * two-to-phi */ \
		f64ptr(rsi+192) = f0;	/* Save value7 */ \
	\
		xmm2[0] -= XMM_BIGVAL1;/* Remove rounding constant */ \
		xmm2[0] *= f64ptr(rbp+384+48);	/* new value8 = val * two-to-phi */ \
		xmm2[0] += f64ptr(rsi+192+16);	/* Add in FFT data */ \
		f64ptr(rsi+192+16) = xmm2[0];	/* Save value8 */ \
	}else{ \
		/* Same as the above but with different addresses required by */ \
		/* the length 80 and 112 FFT lengths */ \
		rax = u8ptr(rdi);	/* Load big vs. little flags */ \
		f0 = TMP56[0];	/* Load integer part of divide by k */ \
		f1 = f64ptr(rsi+8);		/* Load FFT data */ \
		f1 *= f64ptr(rbp+8);		/* Mul values5 by two-to-minus-phi */ \
		f1 *= g->u.xmm.XMM_NORM012_FF[0];	/* Mul by FFTLEN/2 */ \
		f0 += f1;		/* Add in the FFT data */ \
		f0 += xmm2[0];		/* x5 = value + carry */ \
		single_rounding(base2, f0, xmm2[0], f4, rax+8); \
		f0 *= f64ptr(rbp+16+8);	/* new value5 = val * two-to-phi */ \
		f64ptr(rsi+8) = f0;		/* Save value5 */ \
	\
		rax = u8ptr(rdi+1);	/* Load big vs. little flags */ \
		f0 = TMP56[1];	/* Load integer part of divide by k */ \
		f1 = f64ptr(rsi+16+8);	/* Load FFT data */ \
		f1 *= f64ptr(rbp+32+8);	/* Mul values6 by two-to-minus-phi */ \
		f1 *= g->u.xmm.XMM_NORM012_FF[0];	/* Mul by FFTLEN/2 */ \
		f0 += f1;		/* Add in the FFT data */ \
		f0 += xmm2[0];		/* x6 = value + carry */ \
		single_rounding(base2, f0, xmm2[0], f4, rax+8); \
		f0 *= f64ptr(rbp+48+8);	/* new value6 = val * two-to-phi */ \
		f64ptr(rsi+16+8) = f0;	/* Save value6 */ \
	\
		rax = u8ptr(rdi+4);	/* Load big vs. little flags */ \
		f0 = f64ptr(rsi+64+8);	/* Load FFT data */ \
		f0 *= f64ptr(rbp+128+8);	/* Mul values5 by two-to-minus-phi */ \
		f0 *= g->u.xmm.XMM_NORM012_FF[0];	/* Mul by FFTLEN/2 */ \
		f0 += xmm2[0];		/* x7 = value + carry */ \
		single_rounding(base2, f0, xmm2[0], f4, rax+8); \
		f0 *= f64ptr(rbp+128+16+8);	/* new value7 = val * two-to-phi */ \
		f64ptr(rsi+64+8) = f0;	/* Save value7 */ \
	\
		xmm2[0] -= XMM_BIGVAL1;/* Remove rounding constant */ \
		xmm2[0] *= f64ptr(rbp+128+48+8);	/* new value8 = val * two-to-phi */ \
		xmm2[0] += f64ptr(rsi+64+16+8);	/* Add in FFT data */ \
		f64ptr(rsi+64+16+8) = xmm2[0];	/* Save value8 */ \
	} \
}





/* For 2D macros, these registers are set on input: */ \
/* xmm7 = sumout */ \
/* xmm6 = maxerr */ \
/* rbp = pointer to carries */ \
/* rdi = pointer to big/little flags */ \
/* For 2D macros, these registers are set on input: */ \
/* xmm6 = maxerr */ \
/* rbp = pointer to carries */ \
/* rdi = pointer to big/little flags */ \
/* rsi = pointer to the FFT data */ \
/* rbx = pointer two-to-phi column multipliers */ \
/* rdx = pointer two-to-phi group multipliers */ \
/* ecx = big vs. little word flag #2 */ \
/* eax = big vs. little word flag #1 */ \
\
\
/* *************** 2D macro ****************** */ \
/* A pipelined version of this code: */ \
/*	movzx	rax, BYTE PTR [rdi]	;; Load big vs. little flags */ \
/*	xload	xmm0, [rsi+0*dist1]	;; Load values1 */ \
/*	xload	xmm2, [rbx]		;; col two-to-minus-phi */ \
/*	mulpd	xmm2, XMM_TTMP_FUDGE[rax];; Mul by fudge two-to-minus-phi */ \
/*	mulpd	xmm0, [rdx]		;; Mul by grp two-to-minus-phi */ \
/*	mulpd	xmm0, xmm2		;; Mul by fudged col two-to-minus-phi */ \
/*	addpd	xmm0, [rbp+0*16]	;; x1 = values + carry */ \
/*	xload	xmm2, XMM_LIMIT_BIGMAX[rax];; Load maximum * BIGVAL - BIGVAL */ \
/*	addpd	xmm2, xmm0		;; y1 = top bits of x */ \
/*	xload	xmm6, XMM_LIMIT_BIGMAX_NEG[rax];; Load -(maximum*BIGVAL-BIGVAL) */ \
/*	addpd	xmm6, xmm2		;; z1 = y1-(maximum * BIGVAL - BIGVAL) */ \
/*	subpd	xmm0, xmm6		;; rounded value = x1 - z1 */ \
/*	mulpd	xmm2, XMM_LIMIT_INVERSE[rax];; next carry = shifted y1 */ \
/*	xload	xmm4, [rbx]		;; col two-to-phi */ \
/*	mulpd	xmm4, XMM_TTP_FUDGE[rax];; mul by fudge two-to-phi */ \
/*	mulpd	xmm0, [rdx+0*32+16]	;; new value1 = val * grp two-to-phi */ \
/*	mulpd	xmm0, xmm4		;; new value1 *= fudged col two-to-phi */ \
/*	xstore	[rsi+0*dist1], xmm0	;; Save new value1 */ \
/*	xstore	[rbp+0*16], xmm2	;; Save carry */ \
/* */ \
\

#define xnorm_2d_setup(ttp)	/* Precompute FUDGE * col multipliers */ \
ttp(rbp = (rdi+256));		/* Create pointer for more one-byte offsets below */ \
	xmm0 = xptr(rbx);		/* Load col two-to-minus-phi */ \
	xptr(rdi-128) = xmm0;		/* Save ttmp * 1.0,1.0 */ \
ttp(xptr(rdi-128+16) = xmm0);	/* Save ttmp * 1.0,1.0 */ \
ttp(xptr(rdi-128+32) = xmm0);	/* Save ttmp * 1.0,1.0 */ \
ttp(xptr(rdi-128+48) = xmm0);	/* Save ttmp * 1.0,1.0 */ \
ttp(xmm0[0] *= f64ptr2(g->u.xmm.XMM_TTMP_FUDGE,64));/* Compute ttmp * 1.0,B */ \
	xmm1 = xptr(rbx+16);		/* Load col two-to-phi */ \
no##ttp(xptr(rdi+256-128) = xmm1);	/* Save ttp * 1.0,1.0 */ \
ttp(xptr(rbp-128) = xmm1);		/* Save ttp * 1.0,1.0 */ \
ttp(xptr(rbp-128+16) = xmm1);	/* Save ttp * 1.0,1.0 */ \
ttp(xptr(rbp-128+32) = xmm1);	/* Save ttp * 1.0,1.0 */ \
ttp(xptr(rbp-128+48) = xmm1);	/* Save ttp * 1.0,1.0 */ \
ttp(xmm1[0] *= f64ptr2(g->u.xmm.XMM_TTP_FUDGE,64));/* Compute ttp * 1.0,1/B */ \
ttp(xptr(rdi-128+64) = xmm0);	/* Save ttmp * 1.0,B */ \
ttp(xptr(rdi-128+64+16) = xmm0);	/* Save ttmp * 1.0,B */ \
ttp(xptr(rdi-128+64+32) = xmm0);	/* Save ttmp * 1.0,B */ \
ttp(xptr(rdi-128+64+48) = xmm0);	/* Save ttmp * 1.0,B */ \
ttp(shufpd(xmm0, xmm0, 1));		/* swizzle */ \
ttp(xptr(rbp-128+64) = xmm1);	/* Save ttp * 1.0,1/B */ \
ttp(xptr(rbp-128+64+16) = xmm1);	/* Save ttp * 1.0,1/B */ \
ttp(xptr(rbp-128+64+32) = xmm1);	/* Save ttp * 1.0,1/B */ \
ttp(xptr(rbp-128+64+48) = xmm1);	/* Save ttp * 1.0,1/B */ \
ttp(shufpd(xmm1, xmm1, 1));		/* swizzle */ \
ttp(xptr(rdi-128+128) = xmm0);	/* Save ttmp * B,1.0 */ \
ttp(xptr(rdi-128+128+16) = xmm0);	/* Save ttmp * B,1.0 */ \
ttp(xptr(rdi-128+128+32) = xmm0);	/* Save ttmp * B,1.0 */ \
ttp(xptr(rdi-128+128+48) = xmm0);	/* Save ttmp * B,1.0 */ \
ttp(xmm0[0] *= f64ptr2(g->u.xmm.XMM_TTMP_FUDGE,64)); \
ttp(xptr(rbp-128+128) = xmm1);	/* Save ttp * 1/B,1.0 */ \
ttp(xptr(rbp-128+128+16) = xmm1);	/* Save ttp * 1/B,1.0 */ \
ttp(xptr(rbp-128+128+32) = xmm1);	/* Save ttp * 1/B,1.0 */ \
ttp(xptr(rbp-128+128+48) = xmm1);	/* Save ttp * 1/B,1.0 */ \
ttp(xmm1[0] *= f64ptr2(g->u.xmm.XMM_TTP_FUDGE,64)); \
ttp(xptr(rdi-128+192) = xmm0);	/* Save ttmp * B,B */ \
ttp(xptr(rdi-128+192+16) = xmm0);	/* Save ttmp * B,B */ \
ttp(xptr(rdi-128+192+32) = xmm0);	/* Save ttmp * B,B */ \
ttp(xptr(rdi-128+192+48) = xmm0);	/* Save ttmp * B,B */ \
ttp(xptr(rbp-128+192) = xmm1);	/* Save ttp * 1/B,1/B */ \
ttp(xptr(rbp-128+192+16) = xmm1);	/* Save ttp * 1/B,1/B */ \
ttp(xptr(rbp-128+192+32) = xmm1);	/* Save ttp * 1/B,1/B */ \
ttp(xptr(rbp-128+192+48) = xmm1);	/* Save ttp * 1/B,1/B */ \


#if 0 // 32 bit
#define xnorm_2d(ttp, echk, const1, base2, sse4, CARRY0, CARRY1, CARRY2, CARRY3) \
	xmm0 = xptr(rsi+0*16);	/* Load values1 */ \
ttp(xmm2 = xptr(rdx+0*32));	/* grp two-to-minus-phi */ \
ttp(xmm2 *= xptr(rbx+rax));	/* Mul by col two-to-minus-phi */ \
no##ttp(xmm2 = xptr(rbx));		/* two-to-minus-phi */ \
	xmm0 *= xmm2;		/* Mul by fudged col two-to-minus-phi */ \
	xmm1 = xptr(rsi+1*16);	/* Load values2 */ \
ttp(xmm3 = xptr(rdx+1*32));	/* grp two-to-minus-phi */ \
ttp(xmm3 *= xptr(rbx+rcx));	/* Mul by col two-to-minus-phi */ \
no##ttp(xmm3 = xptr(rbx));		/* two-to-minus-phi */ \
	xmm1 *= xmm3;		/* Mul by fudged col two-to-minus-phi */ \
no##const1(echk(error_check_interleaved(sse4, xmm0, xmm4, xmm1, xmm5, xmm6))); \
const1(mul_by_const_interleaved(ttp, base2, echk, sse4, xmm0, xmm4, xmm2, rax, xmm1, xmm5, xmm3, rcx, xmm6)); \
	xmm0 += CARRY0;	/* x1 = values + carry */ \
	xmm1 += CARRY1;	/* x2 = values + carry */ \
\
no##base2(rounding_interleaved(ttp, base2, const1, sse4, xmm0, xmm4, xmm2, rax, xmm1, xmm5, xmm3, rcx)); \
\
base2(xmm2 = xptr2(g->u.xmm.XMM_LIMIT_BIGMAX,rax));/* Load maximum * BIGVAL - BIGVAL */ \
base2(xmm2 += xmm0);		/* y1 = top bits of x */ \
base2(xmm3 = xptr2(g->u.xmm.XMM_LIMIT_BIGMAX,rcx));/* Load maximum * BIGVAL - BIGVAL */ \
base2(xmm3 += xmm1);		/* y2 = top bits of x */ \
base2(const1(xmm4 += xmm2));		/* Add in upper mul-by-const bits */ \
base2(const1(xmm4 *= xptr2(g->u.xmm.XMM_LIMIT_INVERSE,rax)));/* next carry = shifted y1 */ \
base2(no##const1(xmm4 = xptr2(g->u.xmm.XMM_LIMIT_INVERSE,rax)));/* next carry = shifted y1 */ \
base2(no##const1(xmm4 *= xmm2));		/* next carry = shifted y1 */ \
	CARRY0 = xmm4;	/* Save carry1 */ \
base2(xmm2 -= xptr2(g->u.xmm.XMM_LIMIT_BIGMAX,rax));/* z1 = y1 - (maximum*BIGVAL-BIGVAL) */ \
ttp(xmm4 = xptr(rbx+256+rax));	/* col two-to-phi */ \
ttp(xmm4 *= xptr(rdx+0*32+16));	/* two-to-phi = col * grp */ \
ttp(rax = u8ptr(rdi+2));	/* Load next big vs. little flags */ \
base2(const1(xmm5 += xmm3));		/* Add in upper mul-by-const bits */ \
base2(const1(xmm5 *= xptr2(g->u.xmm.XMM_LIMIT_INVERSE,rcx)));/* next carry = shifted y2 */ \
base2(no##const1(xmm5 = xptr2(g->u.xmm.XMM_LIMIT_INVERSE,rcx)));/* next carry = shifted y2 */ \
base2(no##const1(xmm5 *= xmm3));		/* next carry = shifted y2 */ \
	CARRY1 = xmm5;	/* Save carry2 */ \
base2(xmm3 -= xptr2(g->u.xmm.XMM_LIMIT_BIGMAX,rcx));/* z2 = y2 - (maximum*BIGVAL-BIGVAL) */ \
ttp(xmm5 = xptr(rbx+256+rcx));	/* col two-to-phi */ \
ttp(xmm5 *= xptr(rdx+1*32+16));	/* two-to-phi = col * grp */ \
ttp(rcx = u8ptr(rdi+3));	/* Load next big vs. little flags */ \
base2(xmm0 -= xmm2);		/* rounded value = x1 - z1 */ \
base2(xmm1 -= xmm3);		/* rounded value = x2 - z2 */ \
ttp(xmm0 *= xmm4);		/* value1 = rounded value * two-to-phi */ \
	xptr(rsi+0*16) = xmm0;	/* Save new value1 */ \
ttp(xmm1 *= xmm5);		/* value2 = rounded value * two-to-phi */ \
	xptr(rsi+1*16) = xmm1;	/* Save new value2 */ \
\
	xmm0 = xptr(rsi+2*16);	/* Load values1 */ \
ttp(xmm2 = xptr(rdx+2*32));	/* grp two-to-minus-phi */ \
ttp(xmm2 *= xptr(rbx+rax));	/* Mul by col two-to-minus-phi */ \
no##ttp(xmm2 = xptr(rbx));		/* two-to-minus-phi */ \
	xmm0 *= xmm2;		/* Mul by fudged col two-to-minus-phi */ \
	xmm1 = xptr(rsi+3*16);	/* Load values2 */ \
ttp(xmm3 = xptr(rdx+3*32));	/* grp two-to-minus-phi */ \
ttp(xmm3 *= xptr(rbx+rcx));	/* Mul by col two-to-minus-phi */ \
no##ttp(xmm3 = xptr(rbx));		/* two-to-minus-phi */ \
	xmm1 *= xmm3;		/* Mul by fudged col two-to-minus-phi */ \
no##const1(echk(error_check_interleaved(sse4, xmm0, xmm4, xmm1, xmm5, xmm6))); \
const1(mul_by_const_interleaved(ttp, base2, echk, sse4, xmm0, xmm4, xmm2, rax, xmm1, xmm5, xmm3, rcx, xmm6)); \
	xmm0 += CARRY2;	/* x1 = values + carry */ \
	xmm1 += CARRY3;	/* x2 = values + carry */ \
\
no##base2(rounding_interleaved(ttp, base2, const1, sse4, xmm0, xmm4, xmm2, rax, xmm1, xmm5, xmm3, rcx)); \
\
base2(xmm2 = xptr2(g->u.xmm.XMM_LIMIT_BIGMAX,rax));/* Load maximum * BIGVAL - BIGVAL */ \
base2(xmm2 += xmm0);		/* y1 = top bits of x */ \
base2(xmm3 = xptr2(g->u.xmm.XMM_LIMIT_BIGMAX,rcx));/* Load maximum * BIGVAL - BIGVAL */ \
base2(xmm3 += xmm1);		/* y2 = top bits of x */ \
base2(const1(xmm4 += xmm2));		/* Add in upper mul-by-const bits */ \
base2(const1(xmm4 *= xptr2(g->u.xmm.XMM_LIMIT_INVERSE,rax)));/* next carry = shifted y1 */ \
base2(no##const1(xmm4 = xptr2(g->u.xmm.XMM_LIMIT_INVERSE,rax)));/* next carry = shifted y1 */ \
base2(no##const1(xmm4 *= xmm2));		/* next carry = shifted y1 */ \
	CARRY2 = xmm4;	/* Save carry1 */ \
base2(xmm2 -= xptr2(g->u.xmm.XMM_LIMIT_BIGMAX,rax));/* z1 = y1 - (maximum*BIGVAL-BIGVAL) */ \
ttp(xmm4 = xptr(rbx+256+rax));	/* col two-to-phi */ \
ttp(xmm4 *= xptr(rdx+2*32+16));	/* two-to-phi = col * grp */ \
ttp(rax = u8ptr(rdi+4));	/* Load next big vs. little flags */ \
base2(const1(xmm5 += xmm3));		/* Add in upper mul-by-const bits */ \
base2(const1(xmm5 *= xptr2(g->u.xmm.XMM_LIMIT_INVERSE,rcx)));/* next carry = shifted y2 */ \
base2(no##const1(xmm5 = xptr2(g->u.xmm.XMM_LIMIT_INVERSE,rcx)));/* next carry = shifted y2 */ \
base2(no##const1(xmm5 *= xmm3));		/* next carry = shifted y2 */ \
	CARRY3 = xmm5;	/* Save carry2 */ \
base2(xmm3 -= xptr2(g->u.xmm.XMM_LIMIT_BIGMAX,rcx));/* z2 = y2 - (maximum*BIGVAL-BIGVAL) */ \
ttp(xmm5 = xptr(rbx+256+rcx));	/* col two-to-phi */ \
ttp(xmm5 *= xptr(rdx+3*32+16));	/* two-to-phi = col * grp */ \
ttp(rcx = u8ptr(rdi+5));	/* Load next big vs. little flags */ \
base2(xmm0 -= xmm2);		/* rounded value = x1 - z1 */ \
base2(xmm1 -= xmm3);		/* rounded value = x2 - z2 */ \
ttp(xmm0 *= xmm4);		/* value1 *= rounded value * two-to-phi */ \
	xptr(rsi+2*16) = xmm0;	/* Save new value1 */ \
ttp(xmm1 *= xmm5);		/* value2 = rounded value * two-to-phi */ \
	xptr(rsi+3*16) = xmm1;	/* Save new value2 */ \

#else // 64 bit

#define xnorm_2d(ttp, echk, const1, base2, sse4, CARRY0, CARRY1, CARRY2, CARRY3) \
ttp(xmm10 = xptr(rdx+0*32));	/* grp two-to-minus-phi			;P4	;Core2 */ \
ttp(xmm10 *= xptr(rbx+rax));	/* Mul by col two-to-minus-phi		;1-6	;1-5 */ \
no##ttp(xmm10 = xptr(rbx));		/* two-to-minus-phi */ \
	xmm8 = xptr(rsi+0*16);	/* Load values1 */ \
ttp(xmm11 = xptr(rdx+1*32));	/* grp two-to-minus-phi */ \
ttp(xmm11 *= xptr(rbx+rcx));	/* Mul by col two-to-minus-phi		;3-8	;2-6 */ \
no##ttp(xmm11 = xptr(rbx));		/* two-to-minus-phi */ \
ttp(xmm2 = xptr(rdx+2*32));	/* grp two-to-minus-phi */ \
ttp(xmm2 *= xptr(rbx+r8));		/* Mul by col two-to-minus-phi		;5-10	;3-7 */ \
no##ttp(xmm2 = xptr(rbx));		/* two-to-minus-phi */ \
	xmm9 = xptr(rsi+1*16);	/* Load values2 */ \
ttp(xmm3 = xptr(rdx+3*32));	/* grp two-to-minus-phi */ \
ttp(xmm3 *= xptr(rbx+r9));		/* Mul by col two-to-minus-phi		;7-12	;4-8 */ \
no##ttp(xmm3 = xptr(rbx));		/* two-to-minus-phi */ \
	xmm8 *= xmm10;		/* Mul by fudged col two-to-minus-phi	;9-14	;6-10 */ \
	xmm1 = xptr(rsi+3*16);	/* Load values4 */ \
	xmm9 *= xmm11;		/* Mul by fudged col two-to-minus-phi	;11-16	;7-11 */ \
no##const1(echk(error_check_interleaved(sse4, xmm8, xmm12, xmm9, xmm13, xmm6))); \
const1(mul_by_const_interleaved(ttp, base2, echk, sse4, xmm8, xmm12, xmm10, rax, xmm9, xmm13, xmm11, rcx, xmm6)); \
	xmm0 = xptr(rsi+2*16);	/* Load values3 */ \
	xmm0 *= xmm2;		/* Mul by fudged col two-to-minus-phi	;13-18	;8-12 */ \
	xmm1 *= xmm3;		/* Mul by fudged col two-to-minus-phi	;15-20	;9-13 */ \
no##const1(echk(error_check_interleaved(sse4, xmm0, xmm4, xmm1, xmm5, xmm6))); \
const1(mul_by_const_interleaved(ttp, base2, echk, sse4, xmm0, xmm4, xmm2, r8, xmm1, xmm5, xmm3, r9, xmm6)); \
	xmm8 += CARRY0;	/* x1 = values + carry			;16-19	;11-13 */ \
	xmm9 += CARRY1;	/* x2 = values + carry			;18-21	;12-14 */ \
\
no##base2(rounding_interleaved(ttp, base2, const1, sse4, xmm8, xmm12, xmm10, rax, xmm9, xmm13, xmm11, rcx)); \
\
ttp(xmm14 = xptr(rbx+256+rax));	/* col two-to-phi */ \
ttp(xmm14 *= xptr(rdx+0*32+16));	/* two-to-phi = col * grp		;19-24	;12-16 */ \
	xmm0 += CARRY2;	/* x3 = values + carry			;20-23	;13-15 */ \
ttp(xmm15 = xptr(rbx+256+rcx));	/* col two-to-phi */ \
ttp(xmm15 *= xptr(rdx+1*32+16));	/* two-to-phi = col * grp		;21-26	;13-17 */ \
base2(xmm10 = xptr2(g->u.xmm.XMM_LIMIT_BIGMAX,rax));/* Load maximum * BIGVAL - BIGVAL */ \
base2(xmm10 += xmm8);		/* y1 = top bits of x			;22-25	;14-16 */ \
	xmm1 += CARRY3;	/* x4 = values + carry			;24-27	;15-17 */ \
\
no##base2(rounding_interleaved(ttp, base2, const1, sse4, xmm0, xmm4, xmm2, r8, xmm1, xmm5, xmm3, r9)); \
\
base2(xmm11 = xptr2(g->u.xmm.XMM_LIMIT_BIGMAX,rcx));/* Load maximum * BIGVAL - BIGVAL */ \
base2(xmm11 += xmm9);		/* y2 = top bits of x			;26-29	;16-18 */ \
base2(const1(xmm12 += xmm10));		/* Add in upper mul-by-const bits */ \
base2(const1(xmm12 *= xptr2(g->u.xmm.XMM_LIMIT_INVERSE,rax)));/* next carry = shifted y1 */ \
base2(no##const1(xmm12 = xptr2(g->u.xmm.XMM_LIMIT_INVERSE,rax)));/* next carry = shifted y1 */ \
base2(no##const1(xmm12 *= xmm10));	/* next carry = shifted y1		;27-32	;17-21 */ \
base2(xmm2 = xptr2(g->u.xmm.XMM_LIMIT_BIGMAX,r8));/* Load maximum * BIGVAL - BIGVAL */ \
base2(xmm2 += xmm0);		/* y3 = top bits of x			;28-31	;17-19 */ \
base2(xmm3 = xptr2(g->u.xmm.XMM_LIMIT_BIGMAX,r9));/* Load maximum * BIGVAL - BIGVAL */ \
base2(xmm3 += xmm1);		/* y4 = top bits of x			;30-33	;18-20 */ \
base2(const1(xmm13 += xmm11));		/* Add in upper mul-by-const bits */ \
base2(const1(xmm13 *= xptr2(g->u.xmm.XMM_LIMIT_INVERSE,rcx)));/* next carry = shifted y2 */ \
base2(no##const1(xmm13 = xptr2(g->u.xmm.XMM_LIMIT_INVERSE,rcx)));/* next carry = shifted y2 */ \
base2(no##const1(xmm13 *= xmm11));	/* next carry = shifted y2		;31-36	;19-23 */ \
base2(xmm10 -= xptr2(g->u.xmm.XMM_LIMIT_BIGMAX,rax));/* z1 = y1 - (maximum*BIGVAL-BIGVAL) ;32-35	;19-22 */ \
ttp(rax = u8ptr(rdi+4));	/* Load next big vs. little flags */ \
base2(const1(xmm4 += xmm2));		/* Add in upper mul-by-const bits */ \
base2(const1(xmm4 *= xptr2(g->u.xmm.XMM_LIMIT_INVERSE,r8)));/* next carry = shifted y3 */ \
base2(no##const1(xmm4 = xptr2(g->u.xmm.XMM_LIMIT_INVERSE,r8)));/* next carry = shifted y3 */ \
base2(no##const1(xmm4 *= xmm2));		/* next carry = shifted y3		;33-38	;20-24 */ \
base2(xmm11 -= xptr2(g->u.xmm.XMM_LIMIT_BIGMAX,rcx));/* z2 = y2 - (maximum*BIGVAL-BIGVAL) ;34-37	;20-22 */ \
ttp(rcx = u8ptr(rdi+5));	/* Load next big vs. little flags */ \
base2(const1(xmm5 += xmm3));		/* Add in upper mul-by-const bits */ \
base2(const1(xmm5 *= xptr2(g->u.xmm.XMM_LIMIT_INVERSE,r9)));/* next carry = shifted y4 */ \
base2(no##const1(xmm5 = xptr2(g->u.xmm.XMM_LIMIT_INVERSE,r9)));/* next carry = shifted y4 */ \
base2(no##const1(xmm5 *= xmm3));		/* next carry = shifted y4		;35-40	;21-25 */ \
	CARRY0 = xmm12;	/* Save carry1 */ \
base2(xmm2 -= xptr2(g->u.xmm.XMM_LIMIT_BIGMAX,r8));/* z3 = y3 - (maximum*BIGVAL-BIGVAL)	;36-39	;21-23 */ \
ttp(xmm12 = xptr(rbx+256+r8));	/* col two-to-phi */ \
ttp(xmm12 *= xptr(rdx+2*32+16));	/* two-to-phi = col * grp		;37-42	;22-26 */ \
ttp(r8 = u8ptr(rdi+6));	/* Load next big vs. little flags */ \
	CARRY1 = xmm13;	/* Save carry2 */ \
base2(xmm3 -= xptr2(g->u.xmm.XMM_LIMIT_BIGMAX,r9));/* z4 = y4 - (maximum*BIGVAL-BIGVAL)	;38-41	;22-24 */ \
ttp(xmm13 = xptr(rbx+256+r9));	/* col two-to-phi */ \
ttp(xmm13 *= xptr(rdx+3*32+16));	/* two-to-phi = col * grp		;39-43	;23-27 */ \
ttp(r9 = u8ptr(rdi+7));	/* Load next big vs. little flags */ \
base2(xmm8 -= xmm10);		/* rounded value = x1 - z1		;40-43	;23-25 */ \
	CARRY2 = xmm4;	/* Save carry3 */ \
base2(xmm9 -= xmm11);		/* rounded value = x2 - z2		;42-45	;24-26 */ \
	CARRY3 = xmm5;	/* Save carry4 */ \
base2(xmm0 -= xmm2);		/* rounded value = x3 - z3		;44-47	;25-27 */ \
ttp(xmm8 *= xmm14);		/* value1 = rounded value * two-to-phi	;45-50	;26-30 */ \
base2(xmm1 -= xmm3);		/* rounded value = x4 - z4		;46-49	;26-28 */ \
ttp(xmm9 *= xmm15);		/* value2 = rounded value * two-to-phi	;47-52	;27-31 */ \
ttp(xmm0 *= xmm12);		/* value3 = rounded value * two-to-phi	;49-54	;28-32 */ \
ttp(xmm1 *= xmm13);		/* value4 = rounded value * two-to-phi	;51-56	;29-33 */ \
	xptr(rsi+0*16) = xmm8;	/* Save new value1 */ \
	xptr(rsi+1*16) = xmm9;	/* Save new value2 */ \
	xptr(rsi+2*16) = xmm0;	/* Save new value3 */ \
	xptr(rsi+3*16) = xmm1;	/* Save new value4 */ \

#endif
\
/* NOTE: We'd rather store the high FFT carry without the XMM_BIGVAL added in, */ \
/* but there is too much code that expects this (like add and subtract). */ \

#if 0 // 0: 32bit   1: 64bit
#define xnorm_2d_zpad_pre_loop
#define xnorm_2d_zpad_post_loop
#define xnorm_2d_zpad(ttp, echk, const1, base2, sse4, khi, c1, cm1, CARRY0, CARRY1, CARRY2, CARRY3) \
	xmm2 = xptr(rbx+rax);	/* col two-to-minus-phi */ \
ttp(xmm2 *= xptr(rdx+0*32));	/* Mul by grp two-to-minus-phi */ \
	xmm0 = xptr(rsi);		/* Load values1 */ \
	xmm1 = xptr(rsi+2*16);	/* Load values2 */ \
	xmm0 *= xmm2;		/* Mul by fudged col two-to-minus-phi */ \
	xmm1 *= xmm2;		/* Mul by fudged col two-to-minus-phi */ \
\
	xmm3 = CARRY2;	/* Add in previous high FFT data */ \
	split_lower_zpad_word(echk, base2, sse4, xmm0, xmm3, xmm4, rax); \
	CARRY2 = xmm3; \
\
no##const1(xmm0 = g->u.xmm.XMM_K_LO); \
const1(	xmm0 = g->u.xmm.XMM_K_TIMES_MULCONST_LO); \
		xmm0 *= xmm4; \
khi(no##const1(xmm5 = g->u.xmm.XMM_K_HI)); \
khi(const1(xmm5 = g->u.xmm.XMM_K_TIMES_MULCONST_HI)); \
khi(no##base2(xmm5 *= xptr2(g->u.xmm.XMM_LIMIT_INVERSE,rax))); /* Non-base2 rounding needs shifted carry */ \
khi(	xmm5 *= xmm4); \
;  \
		xmm0 += CARRY0;	/* x1 = values + carry */ \
\
c1(	xmm1 *= g->u.xmm.XMM_MINUS_C);	/* Do one mul before split rather than two after split */ \
\
	split_upper_zpad_word(echk, base2, sse4, xmm1, xmm4, xmm2, rax); \
\
no##const1(no##c1(no##cm1(xmm2 *= g->u.xmm.XMM_MINUS_C))); \
no##const1(no##c1(no##cm1(xmm4 *= g->u.xmm.XMM_MINUS_C))); \
const1(		xmm2 *= g->u.xmm.XMM_MINUS_C_TIMES_MULCONST); \
const1(		xmm4 *= g->u.xmm.XMM_MINUS_C_TIMES_MULCONST); \
\
	xmm0 += xmm2;		/* Add upper FFT word to lower FFT word */ \
khi(xmm4 += xmm5);		/* Add upper FFT word to lower FFT word */ \
\
	rounding(ttp, base2, exec, sse4, xmm0, xmm4, xmm2, rax); \
\
ttp(xmm5 = xptr(rbx+256+rax));	/* col two-to-phi */ \
ttp(xmm5 *= xptr(rdx+0*32+16));	/* new value1 = val * grp two-to-phi */ \
ttp(rax = u8ptr(rdi+1));	/* Load next big vs. little flags */ \
ttp(xmm0 *= xmm5);		/* new value1 *= fudged col two-to-phi */ \
	CARRY0 = xmm4;	/* Save carry */ \
	xptr(rsi) = xmm0;		/* Save new value1 */ \
\
	xmm3 = xptr(rbx+rax);	/* col two-to-minus-phi */ \
ttp(xmm3 *= xptr(rdx+1*32));	/* Mul by grp two-to-minus-phi */ \
	xmm5 = xptr(rsi+16);		/* Load high values1 */ \
	xmm1 = xptr(rsi+3*16);	/* Load high values2 */ \
	xmm5 *= xmm3;		/* Mul by fudged col two-to-minus-phi */ \
	xmm1 *= xmm3;		/* Mul by fudged col two-to-minus-phi */ \
\
	xmm3 = CARRY3;	/* Add in previous high FFT data */ \
	split_lower_zpad_word(echk, base2, sse4, xmm5, xmm3, xmm2, rax); \
	CARRY3 = xmm3; \
\
no##const1(xmm0 = g->u.xmm.XMM_K_LO); \
const1(	xmm0 = g->u.xmm.XMM_K_TIMES_MULCONST_LO); \
		xmm0 *= xmm2; \
khi(no##const1(xmm5 = g->u.xmm.XMM_K_HI)); \
khi(const1(xmm5 = g->u.xmm.XMM_K_TIMES_MULCONST_HI)); \
khi(no##base2(xmm5 *= xptr2(g->u.xmm.XMM_LIMIT_INVERSE,rax))); /* Non-base2 rounding needs shifted carry */ \
khi(	xmm5 *= xmm2); \
\
		xmm0 += CARRY1;	/* x2 = values + carry */ \
\
c1(	xmm1 *= g->u.xmm.XMM_MINUS_C);	/* Do one mul before split rather than two after split */ \
\
	split_upper_zpad_word(echk, base2, sse4, xmm1, xmm2, xmm4, rax); \
\
no##const1(no##c1(no##cm1(xmm4 *= g->u.xmm.XMM_MINUS_C))); \
no##const1(no##c1(no##cm1(xmm2 *= g->u.xmm.XMM_MINUS_C))); \
const1(		xmm4 *= g->u.xmm.XMM_MINUS_C_TIMES_MULCONST); \
const1(		xmm2 *= g->u.xmm.XMM_MINUS_C_TIMES_MULCONST); \
\
	xmm0 += xmm4;		/* Add upper FFT word to lower FFT word */ \
khi(xmm2 += xmm5);		/* Add upper FFT word to lower FFT word */ \
\
	rounding(ttp, base2, exec, sse4, xmm0, xmm2, xmm4, rax); \
\
ttp(xmm3 = xptr(rbx+256+rax));	/* col two-to-phi */ \
ttp(xmm3 *= xptr(rdx+1*32+16));	/* new value2 = val * grp two-to-phi */ \
ttp(rax = u8ptr(rdi+4));	/* Load next big vs. little flags */ \
ttp(xmm0 *= xmm3);		/* new value2 *= fudged col two-to-phi */ \
	CARRY1 = xmm2;	/* Save carry */ \
	xptr(rsi+1*16) = xmm0;	/* Save new value2 */ \
\
	xmm1[1] = xmm1[0] = 0.0;		/* new high values = zero */ \
	xptr(rsi+2*16) = xmm1;	/* Zero high value1 */ \
	xptr(rsi+3*16) = xmm1;	/* Zero high value2 */ \

#else
/* 64-bit implementation using extra registers */

#define xnorm_2d_zpad(ttp, echk, const1, base2, sse4, khi, c1, cm1, CARRY0, CARRY1, CARRY2, CARRY3) \
	xmm4  = CARRY0;	/* Preload carries */ \
	xmm12 = CARRY1; \
	xmm3  = CARRY2; \
	xmm11 = CARRY3; \
	xmm2 = xptr(rbx+rax);	/* col two-to-minus-phi */ \
ttp(xmm2 *= xptr(rdx+0*32));	/* Mul by grp two-to-minus-phi */ \
	xmm0 = xptr(rsi);		/* Load values1 */ \
	xmm1 = xptr(rsi+2*16);	/* Load values2 */ \
	xmm0 *= xmm2;		/* Mul by fudged col two-to-minus-phi */ \
	xmm1 *= xmm2;		/* Mul by fudged col two-to-minus-phi */ \
\
	split_lower_zpad_word(echk, base2, sse4, xmm0, xmm3, xmm5, rax); \
\
no##const1(xmm0 = g->u.xmm.XMM_K_LO); \
const1(	xmm0 = g->u.xmm.XMM_K_TIMES_MULCONST_LO); \
		xmm0 *= xmm5; \
\
		xmm0 += xmm4;	/* x1 = values + carry */ \
\
khi(no##const1(xmm4 = g->u.xmm.XMM_K_HI)); \
khi(const1(xmm4 = g->u.xmm.XMM_K_TIMES_MULCONST_HI)); \
khi(no##base2(xmm4 *= xptr2(g->u.xmm.XMM_LIMIT_INVERSE,rax))); /* Non-base2 rounding needs shifted carry */ \
khi(	xmm5 *= xmm4); \
;  \
c1(	xmm1 *= g->u.xmm.XMM_MINUS_C);	/* Do one mul before split rather than two after split */ \
\
	split_upper_zpad_word(echk, base2, sse4, xmm1, xmm4, xmm2, rax); \
\
no##const1(no##c1(no##cm1(xmm2 *= g->u.xmm.XMM_MINUS_C))); \
no##const1(no##c1(no##cm1(xmm4 *= g->u.xmm.XMM_MINUS_C))); \
const1(		xmm2 *= g->u.xmm.XMM_MINUS_C_TIMES_MULCONST); \
const1(		xmm4 *= g->u.xmm.XMM_MINUS_C_TIMES_MULCONST); \
\
	xmm0 += xmm2;		/* Add upper FFT word to lower FFT word */ \
khi(xmm4 += xmm5);		/* Add upper FFT word to lower FFT word */ \
\
ttp(rcx = u8ptr(rdi+1));	/* Load next big vs. little flags */ \
	xmm10 = xptr(rbx+rcx);	/* col two-to-minus-phi */ \
ttp(xmm10 *= xptr(rdx+1*32));	/* Mul by grp two-to-minus-phi */ \
	xmm8 = xptr(rsi+16);		/* Load high values1 */ \
	xmm9 = xptr(rsi+3*16);	/* Load high values2 */ \
	xmm8 *= xmm10;		/* Mul by fudged col two-to-minus-phi */ \
	xmm9 *= xmm10;		/* Mul by fudged col two-to-minus-phi */ \
\
	split_lower_zpad_word(echk, base2, sse4, xmm8, xmm11, xmm13, rcx); \
\
no##const1(xmm8 = g->u.xmm.XMM_K_LO); \
const1(	xmm8 = g->u.xmm.XMM_K_TIMES_MULCONST_LO); \
		xmm8 *= xmm13; \
\
		xmm8 += xmm12;	/* x2 = values + carry */ \
\
khi(no##const1(xmm12 = g->u.xmm.XMM_K_HI)); \
khi(const1(xmm12 = g->u.xmm.XMM_K_TIMES_MULCONST_HI)); \
khi(no##base2(xmm12 *= xptr2(g->u.xmm.XMM_LIMIT_INVERSE,rcx))); /* Non-base2 rounding needs shifted carry */ \
khi(	xmm13 *= xmm12); \
\
c1(	xmm9 *= g->u.xmm.XMM_MINUS_C);	/* Do one mul before split rather than two after split */ \
\
	split_upper_zpad_word(echk, base2, sse4, xmm9, xmm12, xmm10, rcx); \
\
no##const1(no##c1(no##cm1(xmm10 *= g->u.xmm.XMM_MINUS_C))); \
no##const1(no##c1(no##cm1(xmm12 *= g->u.xmm.XMM_MINUS_C))); \
const1(		xmm10 *= g->u.xmm.XMM_MINUS_C_TIMES_MULCONST); \
const1(		xmm12 *= g->u.xmm.XMM_MINUS_C_TIMES_MULCONST); \
\
	xmm8 += xmm10;		/* Add upper FFT word to lower FFT word */ \
khi(xmm12 += xmm13);		/* Add upper FFT word to lower FFT word */ \
\
	rounding_interleaved(ttp, base2, exec, sse4, xmm0, xmm4, xmm2, rax, xmm8, xmm12, xmm10, rcx); \
\
ttp(xmm5 = xptr(rbx+256+rax));	/* col two-to-phi */ \
ttp(xmm5 *= xptr(rdx+0*32+16));	/* new value1 = val * grp two-to-phi */ \
ttp(xmm0 *= xmm5);		/* new value1 *= fudged col two-to-phi */ \
	xptr(rsi) = xmm0;		/* Save new value1 */ \
\
ttp(xmm13 = xptr(rbx+256+rcx));	/* col two-to-phi */ \
ttp(xmm13 *= xptr(rdx+1*32+16));	/* new value2 = val * grp two-to-phi */ \
ttp(rax = u8ptr(rdi+4));	/* Load next big vs. little flags */ \
ttp(xmm8 *= xmm13);		/* new value2 *= fudged col two-to-phi */ \
	xptr(rsi+1*16) = xmm8;	/* Save new value2 */ \
\
	xmm1[1] = xmm1[0] = 0.0;		/* new high values = zero */ \
	xptr(rsi+2*16) = xmm1;	/* Zero high value1 */ \
	xptr(rsi+3*16) = xmm1;	/* Zero high value2 */ \
	CARRY0 = xmm4; /* Store carries */ \
	CARRY1 = xmm12; \
	CARRY2 = xmm3; \
	CARRY3 = xmm11; \

#endif

/* *************** Top carry adjust macro ****************** */ \
/* This macro corrects the carry out of the topmost word when k is not 1. */ \
/* The problem is the top carry is from b^ceil(logb(k)+n) rather than at k*b^n. */ \
/* So we recompute the top carry by multiplying by b^ceil(logb(k)) and then */ \
/* dividing by k.  The integer part is the new carry and the remainder is */ \
/* added back to the top three words. */ \
\
/* The single-pass case, the top carry is in high word of xmm3 */ \

#define xnorm_top_carry_1d \
	xnorm_top_carry_cmn(rsi, xmm3[0], 0); \
\
\
/* The multi-pass case.  The top carry is loaded into xmm7 from the */ \
/* carries array. */ \

#define xnorm_top_carry \
	xnorm_top_carry_cmn(rsi, xmm7[0], 1); \


#define TWOPASS10(A)
#define TWOPASS11(A) A
#define TWOPASS12(A)
#define xnorm_top_carry_cmn(srcreg, xreg, twopass) \
if(g->TOP_CARRY_NEEDS_ADJUSTING == 1) { /* Does top carry need work? */ \
	double f0, f1, tmp6; \
	uintptr_t rax; \
/* Skip this code if K is 1 */ \
\
/* Two pass case - load the carry */ \
TWOPASS1##twopass(rdi = (uintptr_t)g->carries);		/* Addr of the carries */ \
TWOPASS1##twopass(rax = g->addcount1);		/* Load count of carry rows */ \
TWOPASS1##twopass(rax <<= 6);			/* Compute addr of the high carries */ \
TWOPASS1##twopass(rdi += rax); \
TWOPASS1##twopass(xreg = f64ptr(rdi-8));		/* Load very last carry */ \
\
	xreg -= XMM_BIGVAL1;/* Convert top carry from int+BIGVAL state */ \
	tmp6 = xreg; \
\
;	/* We want to calculate carry * b^ceil(logb(k)) / k and */ \
;	/* carry * b^ceil(logb(k)) % k.  This must be done very carefully as */ \
;	/* carry * b^ceil(logb(k)) may not fit in 53 bits. */ \
\
;	/* Here is a strategy that works for k values up to and including 34 bits. */ \
;	/* We do lots of modulo k operations along the way to insure all intermediate */ \
;	/* results are 51 bits or less. */ \
;	/* Calculate y = carry % k.  This will fit in 34 bits. */ \
;	/* Let z = b^ceil(logb(k)) % k.  Precalculate high_17_bits(z) and low_17_bits(z) */ \
;	/* Remainder is (high_17_bits(z) * y % k * 2^17 + low_17_bits(z) * y) % k */ \
\
	f0 = g->INVERSE_K; \
	f0 *= xreg;		/* Mul top carry by 1/k */ \
	f0 = roundf64(f0); /* Integer part */ \
	f0 *= g->K; \
	xreg -= f0;		/* y = carry % k */ \
\
	f0 = g->CARRY_ADJUST1_HI;	/* high_17_bits(z) */ \
	f0 *= xreg;		/* high_17_bits(z) * y */ \
	f1 = g->INVERSE_K; \
	f1 *= f0;		/* Mul high_17_bits(z) * y by 1/k */ \
	f1 = roundf64(f1); /* Integer part */ \
	f1 *= g->K; \
	f0 -= f1;		/* high_17_bits(z) * y % k */ \
	f0 *= 131072.0;		/* high_17_bits(z) * y % k * 2^17 */ \
	xreg *= g->CARRY_ADJUST1_LO;	/* low_17_bits(z) * y */ \
	f0 += xreg;		/* high_17_bits(z) * y % k * 2^17 + low_17_bits(z) * y */ \
\
	f1 = g->INVERSE_K; \
	f1 *= f0;		/* Mul by 1/k */ \
	f1 = roundf64(f1); /* Integer part */ \
	f1 *= g->K; \
	f0 -= f1;		/* Remainder!!! */ \
\
;	/* Finally calculate integer_part = (carry * b^ceil(logb(k)) - remainder) / k */ \
\
	xreg = tmp6;		/* Reload top carry */ \
	xreg *= g->CARRY_ADJUST1;	/* Mul by b^ceil(logb(k)) */ \
	xreg -= f0;		/* Subtract the remainder */ \
	xreg *= g->INVERSE_K;		/* Mul by 1/k */ \
	xreg = roundf64(xreg); /* Integer part of top carry over k */ \
\
;	/* Now add the remainder to the top words */ \
\
	f0 *= g->CARRY_ADJUST2;	/* Shift remainder */ \
	f1 = roundf64(f0); /* Integer part of shifted remainder */ \
	f0 -= f1;		/* Fractional part of shifted remainder */ \
	f1 *= g->CARRY_ADJUST3;	/* Weight integer part */ \
\
	TWOPASSB##twopass(srcreg) /* TWOPASSB0, TWOPASSB1, or TWOPASSB2 */ \
\
	xreg += XMM_BIGVAL1;/* Restore carry to int+BIGVAL state */ \
\
	TWOPASS1##twopass(f64ptr(rdi-8) = xreg);		/* Save very last carry */ \
}
#define TWOPASSB1(srcreg) /* Two pass scratch area case */ \
	rax = g->HIGH_SCRATCH1_OFFSET; /* Add integer part to top word */ \
	f1 += f64ptr(srcreg+rax); \
	f64ptr(srcreg+rax) = f1; \
	f0 *= g->CARRY_ADJUST4;	/* Shift fractional part */ \
	f0 = roundf64(f0); \
	f0 *= g->CARRY_ADJUST5;	/* Weight fractional part */ \
	rax = g->HIGH_SCRATCH2_OFFSET;/* Add frac part to top-1 word */ \
	f0 += f64ptr(srcreg+rax); \
	f64ptr(srcreg+rax) = f0; \

#define TWOPASSB2(srcreg) /* Two pass FFT data case */ \
	rax = g->HIGH_WORD1_OFFSET;	/* Add integer part to top word */ \
	f1 += f64ptr(srcreg+rax); \
	f64ptr(srcreg+rax) = f1; \
	f0 *= g->CARRY_ADJUST4;	/* Shift fractional part */ \
	f0 = roundf64(f0); \
	f0 *= g->CARRY_ADJUST5;	/* Weight fractional part */ \
	rax = g->HIGH_WORD2_OFFSET;	/* Add frac part to top-1 word */ \
	f0 += f64ptr(srcreg+rax); \
	f64ptr(srcreg+rax) = f0; \

#define TWOPASSB0(srcreg) /* Single pass case */ \
	rax = g->HIGH_WORD1_OFFSET;	/* Add integer part to top word */ \
	f1 += f64ptr(srcreg+rax); \
	f64ptr(srcreg+rax) = f1; \
	f0 *= g->CARRY_ADJUST4;	/* Shift fractional part */ \
	f1 = roundf64(f0); /* Integer part of shifted fractional */ \
	f0 -= f1;		/* Fractional part */ \
	f1 *= g->CARRY_ADJUST5;	/* Weight integer part */ \
	rax = g->HIGH_WORD2_OFFSET;	/* Add frac part to top-1 word */ \
	f1 += f64ptr(srcreg+rax); \
	f64ptr(srcreg+rax) = f1; \
	f0 *= g->CARRY_ADJUST6;	/* Shift fractional part */ \
	f0 = roundf64(f0); \
	f0 *= g->CARRY_ADJUST7;	/* Weight fractional part */ \
	rax = g->HIGH_WORD3_OFFSET;	/* Add frac part to top-2 word */ \
	f0 += f64ptr(srcreg+rax); \
	f64ptr(srcreg+rax) = f0; \



/* *************** 2D followup macro ****************** */ \
/* This macro finishes the normalize process by adding the final carries */ \
/* back into the appropriate FFT values. */ \
/* rsi = pointer to carries */ \
/* rbp = pointer to FFT data */ \
/* rdi = pointer to big/little flags */ \
/* rbx = pointer two-to-phi column multipliers */ \
/* rdx = pointer two-to-phi group multipliers */ \
/* rcx = big vs. little word flag #2 */ \
/* rax = big vs. little word flag #1 */ \
\

#define xnorm012_2d_part1 \
\
if(g->this_block == 0) {		/* Are we carrying into the first data block? */ \
/* If not, skip wrapping carries and negating the last carry */ \
	\
	rax = g->count3;		/* Load 3 section counts */ \
	\
	do{ \
		rdi = rsi;		/* Save section start */ \
		rbx = rax;		/* Form count for this section */ \
		rbx &= 0x7FF; \
		if(!rbx) break;			/* No cache lines to do.  We're all done! */ \
		rax >>= 11;			/* Move counts list along */ \
		do{ \
			xmm4[0] = f64ptr(rsi+0*16);	/* Load low carry word */ \
			f64ptr(rsi+0*16) = xmm0[0];	/* Save prev cache line's high carry in low word */ \
			xmm0[0] = f64ptr(rsi+0*16+8);	/* Load high carry for next cache line */ \
			f64ptr(rsi+0*16+8) = xmm4[0];	/* Move this cache line's low carry to high word */ \
			xmm4[0] = f64ptr(rsi+1*16);	/* Load low carry word */ \
			f64ptr(rsi+1*16) = xmm1[0];	/* Save prev cache line's high carry in low word */ \
			xmm1[0] = f64ptr(rsi+1*16+8);	/* Load high carry for next cache line */ \
			f64ptr(rsi+1*16+8) = xmm4[0];	/* Move this cache line's low carry to high word */ \
			xmm4[0] = f64ptr(rsi+2*16);	/* Load low carry word */ \
			f64ptr(rsi+2*16) = xmm2[0];	/* Save prev cache line's high carry in low word */ \
			xmm2[0] = f64ptr(rsi+2*16+8);	/* Load high carry for next cache line */ \
			f64ptr(rsi+2*16+8) = xmm4[0];	/* Move this cache line's low carry to high word */ \
			xmm4[0] = f64ptr(rsi+3*16);	/* Load low carry word */ \
			f64ptr(rsi+3*16) = xmm3[0];	/* Save prev cache line's high carry in low word */ \
			xmm3[0] = f64ptr(rsi+3*16+8);	/* Load high carry for next cache line */ \
			f64ptr(rsi+3*16+8) = xmm4[0];	/* Move this cache line's low carry to high word */ \
			rsi += 64;			/* Next carry cache line */ \
			rbx--;			/* Test loop counter */ \
		}while(rbx);		/* Next carry row in section */ \
	\
		f64ptr(rdi+1*16) = xmm0[0];	/* Do extra shuffling of section's last cache line */ \
		xmm0[0] = xmm1[0]; \
		f64ptr(rdi+3*16) = xmm2[0]; \
		xmm2[0] = xmm3[0]; \
	}while(1);			/* Next section */ \
	\
	xmm2[0] -= XMM_BIGVAL1; \
	xmm2[0] *= g->u.xmm.XMM_MINUS_C[0];	/* Negate the very last carry */ \
	xmm2[0] += XMM_BIGVAL1; \
	rsi = (uintptr_t)g->carries;		/* Reload carries array pointer */ \
	f64ptr(rsi+0*16) = xmm2[0];	/* Move last cache line's high carries into first cache line */ \
	f64ptr(rsi+2*16) = xmm0[0]; \
}



#define xnorm012_2d(base2) \
\
;	/* If k or c is more than one, then there will be fewer bits-per-word. */ \
;	/* This means the carry may need to be spread over 4 words instead */ \
;	/* of just 2. */ \
\
if(g->SPREAD_CARRY_OVER_EXTRA_WORDS != 1) { /* Are there few bits per word? */ \
	/* je	hard; */			/* Yes, go do it the hard way */ \
\
	rax = u8ptr(rdi+0);	/* Load big vs. little flag */ \
	rcx = u8ptr(rdi+4); \
	xmm0 = xptr(rbp+0*16);	/* FFT data */ \
	xmm6 = xptr(rbx+0*32);	/* col two-to-minus-phi */ \
	xmm6 *= g->u.xmm.XMM_NORM012_FF;	/* Mul by FFTLEN/2 */ \
	xmm0 *= xptr(rdx+0*32);	/* mul by grp two-to-minus-phi */ \
	xmm0 *= xptr2(g->u.xmm.XMM_TTMP_FUDGE,rax);/* mul by fudge two-to-minus-phi */ \
	xmm0 *= xmm6;		/* mul by col two-to-minus-phi */ \
	xmm0 += xptr(rsi+0*16);	/* Add in the carry */ \
	rounding(exec, base2, noexec, noexec, xmm0, xmm5, xmm4, rax); \
	xmm5 -= XMM_BIGVAL2; \
	xmm7 = xptr(rdx+0*32+16);	/* grp two-to-phi */ \
	xmm0 *= xmm7;		/* value *= grp two-to-phi */ \
	xmm0 *= xptr2(g->u.xmm.XMM_TTP_FUDGE,rax);/* value *= fudge two-to-phi */ \
	xmm2 = xptr(rbx+0*32+16);	/* col two-to-phi */ \
	xmm0 *= xmm2;		/* value *= col two-to-phi */ \
	xptr(rbp+0*16) = xmm0;	/* Save FFT data */ \
	xmm5 *= xmm7;		/* carry *= grp two-to-phi */ \
	xmm5 *= xptr2(g->u.xmm.XMM_TTP_FUDGE,rcx);/* carry *= fudge two-to-phi */ \
	xmm3 = xptr(rbx+1*32+16);	/* col two-to-phi */ \
	xmm5 *= xmm3;		/* carry *= col two-to-phi */ \
	xmm5 += xptr(rbp+4*16);	/* Add carry and FFT data */ \
	xptr(rbp+4*16) = xmm5;	/* Save FFT data */ \
\
	rax = u8ptr(rdi+1);	/* Load big vs. little flag */ \
	rcx = u8ptr(rdi+5); \
	xmm1 = xptr(rbp+1*16);	/* FFT data */ \
	xmm1 *= xptr(rdx+1*32);	/* mul by grp two-to-minus-phi */ \
	xmm1 *= xptr2(g->u.xmm.XMM_TTMP_FUDGE,rax);/* mul by fudge two-to-minus-phi */ \
	xmm1 *= xmm6;		/* mul by col two-to-minus-phi */ \
	xmm1 += xptr(rsi+1*16);	/* Add in the carry */ \
	rounding(exec, base2, noexec, noexec, xmm1, xmm5, xmm4, rax); \
	xmm5 -= XMM_BIGVAL2; \
	xmm7 = xptr(rdx+1*32+16);	/* grp two-to-phi */ \
	xmm1 *= xmm7;		/* value *= grp two-to-phi */ \
	xmm1 *= xptr2(g->u.xmm.XMM_TTP_FUDGE,rax);/* value *= fudge two-to-phi */ \
	xmm1 *= xmm2;		/* value *= col two-to-phi */ \
	xptr(rbp+1*16) = xmm1;	/* Save FFT data */ \
	xmm5 *= xmm7;		/* carry *= grp two-to-phi */ \
	xmm5 *= xptr2(g->u.xmm.XMM_TTP_FUDGE,rcx);/* carry *= fudge two-to-phi */ \
	xmm5 *= xmm3;		/* carry *= col two-to-phi */ \
	xmm5 += xptr(rbp+5*16);	/* Add carry and FFT data */ \
	xptr(rbp+5*16) = xmm5;	/* Save FFT data */ \
\
	rax = u8ptr(rdi+2);	/* Load big vs. little flag */ \
	rcx = u8ptr(rdi+6); \
	xmm0 = xptr(rbp+2*16);	/* FFT data */ \
	xmm0 *= xptr(rdx+2*32);	/* mul by grp two-to-minus-phi */ \
	xmm0 *= xptr2(g->u.xmm.XMM_TTMP_FUDGE,rax);/* mul by fudge two-to-minus-phi */ \
	xmm0 *= xmm6;		/* mul by col two-to-minus-phi */ \
	xmm0 += xptr(rsi+2*16);	/* Add in the carry */ \
	rounding(exec, base2, noexec, noexec, xmm0, xmm5, xmm4, rax); \
	xmm5 -= XMM_BIGVAL2; \
	xmm7 = xptr(rdx+2*32+16);	/* grp two-to-phi */ \
	xmm0 *= xmm7;		/* value *= grp two-to-phi */ \
	xmm0 *= xptr2(g->u.xmm.XMM_TTP_FUDGE,rax);/* value *= fudge two-to-phi */ \
	xmm0 *= xmm2;		/* value *= col two-to-phi */ \
	xptr(rbp+2*16) = xmm0;	/* Save FFT data */ \
	xmm5 *= xmm7;		/* carry *= grp two-to-phi */ \
	xmm5 *= xptr2(g->u.xmm.XMM_TTP_FUDGE,rcx);/* carry *= fudge two-to-phi */ \
	xmm5 *= xmm3;		/* carry *= col two-to-phi */ \
	xmm5 += xptr(rbp+6*16);	/* Add carry and FFT data */ \
	xptr(rbp+6*16) = xmm5;	/* Save FFT data */ \
\
	rax = u8ptr(rdi+3);	/* Load big vs. little flag */ \
	rcx = u8ptr(rdi+7); \
	xmm1 = xptr(rbp+3*16);	/* FFT data */ \
	xmm1 *= xptr(rdx+3*32);	/* mul by grp two-to-minus-phi */ \
	xmm1 *= xptr2(g->u.xmm.XMM_TTMP_FUDGE,rax);/* mul by fudge two-to-minus-phi */ \
	xmm1 *= xmm6;		/* mul by col two-to-minus-phi */ \
	xmm1 += xptr(rsi+3*16);	/* Add in the carry */ \
	rounding(exec, base2, noexec, noexec, xmm1, xmm5, xmm4, rax); \
	xmm5 -= XMM_BIGVAL2; \
	xmm7 = xptr(rdx+3*32+16);	/* grp two-to-phi */ \
	xmm1 *= xmm7;		/* value *= grp two-to-phi */ \
	xmm1 *= xptr2(g->u.xmm.XMM_TTP_FUDGE,rax);/* value *= fudge two-to-phi */ \
	xmm1 *= xmm2;		/* value *= col two-to-phi */ \
	xptr(rbp+3*16) = xmm1;	/* Save FFT data */ \
	xmm5 *= xmm7;		/* carry *= grp two-to-phi */ \
	xmm5 *= xptr2(g->u.xmm.XMM_TTP_FUDGE,rcx);/* carry *= fudge two-to-phi */ \
	xmm5 *= xmm3;	/* carry *= col two-to-phi */ \
	xmm5 += xptr(rbp+7*16);	/* Add carry and FFT data */ \
	xptr(rbp+7*16) = xmm5;	/* Save FFT data */ \
}else{ \
/* Same as above, but spread carry over 6 words */ \
	rax = u8ptr(rdi+0);	/* Load big vs. little flag */ \
	xmm0 = xptr(rbp+0*16);	/* FFT data */ \
	xmm6 = xptr(rdx+0*32);	/* grp two-to-minus-phi */ \
	xmm6 *= g->u.xmm.XMM_NORM012_FF;	/* Mul by FFTLEN/2 */ \
	xmm0 *= xptr2(g->u.xmm.XMM_TTMP_FUDGE,rax);/* mul by fudge two-to-minus-phi */ \
	xmm0 *= xmm6;		/* mul by grp two-to-minus-phi */ \
	xmm0 *= xptr(rbx+0*32);	/* mul by col two-to-minus-phi */ \
	xmm0 += xptr(rsi+0*16);	/* Add in the carry */ \
	rounding(exec, base2, noexec, noexec, xmm0, xmm5, xmm4, rax); \
	xmm7 = xptr(rdx+0*32+16);	/* grp two-to-phi */ \
	xmm0 *= xmm7;		/* value *= grp two-to-phi */ \
	xmm0 *= xptr2(g->u.xmm.XMM_TTP_FUDGE,rax);/* value *= fudge two-to-phi */ \
	xmm0 *= xptr(rbx+0*32+16);	/* value *= col two-to-phi */ \
	xptr(rbp+0*16) = xmm0;	/* Save FFT data */ \
\
	rax = u8ptr(rdi+4);	/* Load big vs. little flag */ \
	xmm0 = xptr(rbp+4*16);	/* FFT data */ \
	xmm0 *= xmm6;		/* mul by grp two-to-minus-phi */ \
	xmm0 *= xptr2(g->u.xmm.XMM_TTMP_FUDGE,rax);/* mul by fudge two-to-minus-phi */ \
	xmm0 *= xptr(rbx+1*32);	/* mul by col two-to-minus-phi */ \
	xmm0 += xmm5;		/* Add in the carry */ \
	rounding(exec, base2, noexec, noexec, xmm0, xmm5, xmm4, rax); \
	xmm0 *= xmm7;		/* value *= grp two-to-phi */ \
	xmm0 *= xptr2(g->u.xmm.XMM_TTP_FUDGE,rax);/* value *= fudge two-to-phi */ \
	xmm0 *= xptr(rbx+1*32+16);	/* value *= col two-to-phi */ \
	xptr(rbp+4*16) = xmm0;	/* Save FFT data */ \
\
	rcx = g->BIGLIT_INCR2;/* Different clm values step through */ \
;					/* big/lit array differently */ \
	rax = u8ptr(rdi+rcx);	/* Load big vs. little flag */ \
	xmm0 = xptr(rbp+8*16);	/* FFT data */ \
	xmm0 *= xmm6;		/* mul by grp two-to-minus-phi */ \
	xmm0 *= xptr2(g->u.xmm.XMM_TTMP_FUDGE,rax);/* mul by fudge two-to-minus-phi */ \
	xmm0 *= xptr(rbx+2*32);	/* mul by col two-to-minus-phi */ \
	xmm0 += xmm5;		/* Add in the carry */ \
	rounding(exec, base2, noexec, noexec, xmm0, xmm5, xmm4, rax); \
	xmm0 *= xmm7;		/* value *= grp two-to-phi */ \
	xmm0 *= xptr2(g->u.xmm.XMM_TTP_FUDGE,rax);/* value *= fudge two-to-phi */ \
	xmm0 *= xptr(rbx+2*32+16);	/* value *= col two-to-phi */ \
	xptr(rbp+8*16) = xmm0;	/* Save FFT data */ \
\
	rax = u8ptr(rdi+rcx+4); /* Load big vs. little flag */ \
	xmm0 = xptr(rbp+12*16);	/* FFT data */ \
	xmm0 *= xmm6;		/* mul by grp two-to-minus-phi */ \
	xmm0 *= xptr2(g->u.xmm.XMM_TTMP_FUDGE,rax);/* mul by fudge two-to-minus-phi */ \
	xmm0 *= xptr(rbx+3*32);	/* mul by col two-to-minus-phi */ \
	xmm0 += xmm5;		/* Add in the carry */ \
	rounding(exec, base2, noexec, noexec, xmm0, xmm5, xmm4, rax); \
	xmm0 *= xmm7;		/* value *= grp two-to-phi */ \
	xmm0 *= xptr2(g->u.xmm.XMM_TTP_FUDGE,rax);/* value *= fudge two-to-phi */ \
	xmm0 *= xptr(rbx+3*32+16);	/* value *= col two-to-phi */ \
	xptr(rbp+12*16) = xmm0;	/* Save FFT data */ \
\
	rcx = g->BIGLIT_INCR4;/* Different clm values step through */ \
;					/* big/lit array differently */ \
	rax = u8ptr(rdi+rcx);	/* Load big vs. little flag */ \
	xmm0 = xptr(rbp+16*16);	/* FFT data */ \
	xmm0 *= xmm6;		/* mul by grp two-to-minus-phi */ \
	xmm0 *= xptr2(g->u.xmm.XMM_TTMP_FUDGE,rax);/* mul by fudge two-to-minus-phi */ \
	xmm0 *= xptr(rbx+4*32);	/* mul by col two-to-minus-phi */ \
	xmm0 += xmm5;		/* Add in the carry */ \
	rounding(exec, base2, noexec, noexec, xmm0, xmm5, xmm4, rax); \
	xmm0 *= xmm7;		/* value *= grp two-to-phi */ \
	xmm0 *= xptr2(g->u.xmm.XMM_TTP_FUDGE,rax);/* value *= fudge two-to-phi */ \
	xmm0 *= xptr(rbx+4*32+16);	/* value *= col two-to-phi */ \
	xptr(rbp+16*16) = xmm0;	/* Save FFT data */ \
\
	rax = u8ptr(rdi+rcx+4); \
	xmm5 -= XMM_BIGVAL2; \
	xmm5 *= xmm7;		/* carry *= grp two-to-phi */ \
	xmm5 *= xptr2(g->u.xmm.XMM_TTP_FUDGE,rax);/* carry *= fudge two-to-phi */ \
	xmm5 *= xptr(rbx+5*32+16);	/* carry *= col two-to-phi */ \
	xmm5 += xptr(rbp+20*16);	/* Add carry and FFT data */ \
	xptr(rbp+20*16) = xmm5;	/* Save FFT data */ \
\
\
	rax = u8ptr(rdi+1);	/* Load big vs. little flag */ \
	xmm1 = xptr(rbp+1*16);	/* FFT data */ \
	xmm6 = xptr(rdx+1*32);	/* grp two-to-minus-phi */ \
	xmm6 *= g->u.xmm.XMM_NORM012_FF;	/* Mul by FFTLEN/2 */ \
	xmm1 *= xptr2(g->u.xmm.XMM_TTMP_FUDGE,rax);/* mul by fudge two-to-minus-phi */ \
	xmm1 *= xmm6;		/* mul by grp two-to-minus-phi */ \
	xmm1 *= xptr(rbx+0*32);	/* mul by col two-to-minus-phi */ \
	xmm1 += xptr(rsi+1*16);	/* Add in the carry */ \
	rounding(exec, base2, noexec, noexec, xmm1, xmm5, xmm4, rax); \
	xmm7 = xptr(rdx+1*32+16);	/* grp two-to-phi */ \
	xmm1 *= xmm7;		/* value *= grp two-to-phi */ \
	xmm1 *= xptr2(g->u.xmm.XMM_TTP_FUDGE,rax);/* value *= fudge two-to-phi */ \
	xmm1 *= xptr(rbx+0*32+16);	/* value *= col two-to-phi */ \
	xptr(rbp+1*16) = xmm1;	/* Save FFT data */ \
\
	rax = u8ptr(rdi+5);	/* Load big vs. little flag */ \
	xmm1 = xptr(rbp+5*16);	/* FFT data */ \
	xmm1 *= xmm6;		/* mul by grp two-to-minus-phi */ \
	xmm1 *= xptr2(g->u.xmm.XMM_TTMP_FUDGE,rax);/* mul by fudge two-to-minus-phi */ \
	xmm1 *= xptr(rbx+1*32);	/* mul by col two-to-minus-phi */ \
	xmm1 += xmm5;		/* Add in the carry */ \
	rounding(exec, base2, noexec, noexec, xmm1, xmm5, xmm4, rax); \
	xmm1 *= xmm7;		/* value *= grp two-to-phi */ \
	xmm1 *= xptr2(g->u.xmm.XMM_TTP_FUDGE,rax);/* value *= fudge two-to-phi */ \
	xmm1 *= xptr(rbx+1*32+16);	/* value *= col two-to-phi */ \
	xptr(rbp+5*16) = xmm1;	/* Save FFT data */ \
\
	rcx = g->BIGLIT_INCR2;/* Different clm values step through */ \
;					/* big/lit array differently */ \
	rax = u8ptr(rdi+rcx+1); /* Load big vs. little flag */ \
	xmm1 = xptr(rbp+9*16);	/* FFT data */ \
	xmm1 *= xmm6;		/* mul by grp two-to-minus-phi */ \
	xmm1 *= xptr2(g->u.xmm.XMM_TTMP_FUDGE,rax);/* mul by fudge two-to-minus-phi */ \
	xmm1 *= xptr(rbx+2*32);	/* mul by col two-to-minus-phi */ \
	xmm1 += xmm5;		/* Add in the carry */ \
	rounding(exec, base2, noexec, noexec, xmm1, xmm5, xmm4, rax); \
	xmm1 *= xmm7;		/* value *= grp two-to-phi */ \
	xmm1 *= xptr2(g->u.xmm.XMM_TTP_FUDGE,rax);/* value *= fudge two-to-phi */ \
	xmm1 *= xptr(rbx+2*32+16);	/* value *= col two-to-phi */ \
	xptr(rbp+9*16) = xmm1;	/* Save FFT data */ \
\
	rax = u8ptr(rdi+rcx+5); /* Load big vs. little flag */ \
	xmm1 = xptr(rbp+13*16);	/* FFT data */ \
	xmm1 *= xmm6;		/* mul by grp two-to-minus-phi */ \
	xmm1 *= xptr2(g->u.xmm.XMM_TTMP_FUDGE,rax);/* mul by fudge two-to-minus-phi */ \
	xmm1 *= xptr(rbx+3*32);	/* mul by col two-to-minus-phi */ \
	xmm1 += xmm5;		/* Add in the carry */ \
	rounding(exec, base2, noexec, noexec, xmm1, xmm5, xmm4, rax); \
	xmm1 *= xmm7;		/* value *= grp two-to-phi */ \
	xmm1 *= xptr2(g->u.xmm.XMM_TTP_FUDGE,rax);/* value *= fudge two-to-phi */ \
	xmm1 *= xptr(rbx+3*32+16);	/* value *= col two-to-phi */ \
	xptr(rbp+13*16) = xmm1;	/* Save FFT data */ \
\
	rcx = g->BIGLIT_INCR4;/* Different clm values step through */ \
;					/* big/lit array differently */ \
	rax = u8ptr(rdi+rcx+1); /* Load big vs. little flag */ \
	xmm1 = xptr(rbp+17*16);	/* FFT data */ \
	xmm1 *= xmm6;		/* mul by grp two-to-minus-phi */ \
	xmm1 *= xptr2(g->u.xmm.XMM_TTMP_FUDGE,rax);/* mul by fudge two-to-minus-phi */ \
	xmm1 *= xptr(rbx+4*32);	/* mul by col two-to-minus-phi */ \
	xmm1 += xmm5;		/* Add in the carry */ \
	rounding(exec, base2, noexec, noexec, xmm1, xmm5, xmm4, rax); \
	xmm1 *= xmm7;		/* value *= grp two-to-phi */ \
	xmm1 *= xptr2(g->u.xmm.XMM_TTP_FUDGE,rax);/* value *= fudge two-to-phi */ \
	xmm1 *= xptr(rbx+4*32+16);	/* value *= col two-to-phi */ \
	xptr(rbp+17*16) = xmm1;	/* Save FFT data */ \
\
	rax = u8ptr(rdi+rcx+5); \
	xmm5 -= XMM_BIGVAL2; \
	xmm5 *= xmm7;		/* carry *= grp two-to-phi */ \
	xmm5 *= xptr2(g->u.xmm.XMM_TTP_FUDGE,rax);/* carry *= fudge two-to-phi */ \
	xmm5 *= xptr(rbx+5*32+16);	/* carry *= col two-to-phi */ \
	xmm5 += xptr(rbp+21*16);	/* Add carry and FFT data */ \
	xptr(rbp+21*16) = xmm5;	/* Save FFT data */ \
\
	rax = u8ptr(rdi+2);	/* Load big vs. little flag */ \
	xmm2 = xptr(rbp+2*16);	/* FFT data */ \
	xmm6 = xptr(rdx+2*32);	/* grp two-to-minus-phi */ \
	xmm6 *= g->u.xmm.XMM_NORM012_FF;	/* Mul by FFTLEN/2 */ \
	xmm2 *= xptr2(g->u.xmm.XMM_TTMP_FUDGE,rax);/* mul by fudge two-to-minus-phi */ \
	xmm2 *= xmm6;		/* mul by grp two-to-minus-phi */ \
	xmm2 *= xptr(rbx+0*32);	/* mul by col two-to-minus-phi */ \
	xmm2 += xptr(rsi+2*16);	/* Add in the carry */ \
	rounding(exec, base2, noexec, noexec, xmm2, xmm5, xmm4, rax); \
	xmm7 = xptr(rdx+2*32+16);	/* grp two-to-phi */ \
	xmm2 *= xmm7;		/* value *= grp two-to-phi */ \
	xmm2 *= xptr2(g->u.xmm.XMM_TTP_FUDGE,rax);/* value *= fudge two-to-phi */ \
	xmm2 *= xptr(rbx+0*32+16);	/* value *= col two-to-phi */ \
	xptr(rbp+2*16) = xmm2;	/* Save FFT data */ \
\
	rax = u8ptr(rdi+6);	/* Load big vs. little flag */ \
	xmm2 = xptr(rbp+6*16);	/* FFT data */ \
	xmm2 *= xmm6;		/* mul by grp two-to-minus-phi */ \
	xmm2 *= xptr2(g->u.xmm.XMM_TTMP_FUDGE,rax);/* mul by fudge two-to-minus-phi */ \
	xmm2 *= xptr(rbx+1*32);	/* mul by col two-to-minus-phi */ \
	xmm2 += xmm5;		/* Add in the carry */ \
	rounding(exec, base2, noexec, noexec, xmm2, xmm5, xmm4, rax); \
	xmm2 *= xmm7;		/* value *= grp two-to-phi */ \
	xmm2 *= xptr2(g->u.xmm.XMM_TTP_FUDGE,rax);/* value *= fudge two-to-phi */ \
	xmm2 *= xptr(rbx+1*32+16);	/* value *= col two-to-phi */ \
	xptr(rbp+6*16) = xmm2;	/* Save FFT data */ \
\
	rcx = g->BIGLIT_INCR2;/* Different clm values step through */ \
;					/* big/lit array differently */ \
	rax = u8ptr(rdi+rcx+2); /* Load big vs. little flag */ \
	xmm2 = xptr(rbp+10*16);	/* FFT data */ \
	xmm2 *= xmm6;		/* mul by grp two-to-minus-phi */ \
	xmm2 *= xptr2(g->u.xmm.XMM_TTMP_FUDGE,rax);/* mul by fudge two-to-minus-phi */ \
	xmm2 *= xptr(rbx+2*32);	/* mul by col two-to-minus-phi */ \
	xmm2 += xmm5;		/* Add in the carry */ \
	rounding(exec, base2, noexec, noexec, xmm2, xmm5, xmm4, rax); \
	xmm2 *= xmm7;		/* value *= grp two-to-phi */ \
	xmm2 *= xptr2(g->u.xmm.XMM_TTP_FUDGE,rax);/* value *= fudge two-to-phi */ \
	xmm2 *= xptr(rbx+2*32+16);	/* value *= col two-to-phi */ \
	xptr(rbp+10*16) = xmm2;	/* Save FFT data */ \
\
	rax = u8ptr(rdi+rcx+6); /* Load big vs. little flag */ \
	xmm2 = xptr(rbp+14*16);	/* FFT data */ \
	xmm2 *= xmm6;		/* mul by grp two-to-minus-phi */ \
	xmm2 *= xptr2(g->u.xmm.XMM_TTMP_FUDGE,rax);/* mul by fudge two-to-minus-phi */ \
	xmm2 *= xptr(rbx+3*32);	/* mul by col two-to-minus-phi */ \
	xmm2 += xmm5;		/* Add in the carry */ \
	rounding(exec, base2, noexec, noexec, xmm2, xmm5, xmm4, rax); \
	xmm2 *= xmm7;		/* value *= grp two-to-phi */ \
	xmm2 *= xptr2(g->u.xmm.XMM_TTP_FUDGE,rax);/* value *= fudge two-to-phi */ \
	xmm2 *= xptr(rbx+3*32+16);	/* value *= col two-to-phi */ \
	xptr(rbp+14*16) = xmm2;	/* Save FFT data */ \
\
	rcx = g->BIGLIT_INCR4;/* Different clm values step through */ \
;					/* big/lit array differently */ \
	rax = u8ptr(rdi+rcx+2); /* Load big vs. little flag */ \
	xmm2 = xptr(rbp+18*16);	/* FFT data */ \
	xmm2 *= xmm6;		/* mul by grp two-to-minus-phi */ \
	xmm2 *= xptr2(g->u.xmm.XMM_TTMP_FUDGE,rax);/* mul by fudge two-to-minus-phi */ \
	xmm2 *= xptr(rbx+4*32);	/* mul by col two-to-minus-phi */ \
	xmm2 += xmm5;		/* Add in the carry */ \
	rounding(exec, base2, noexec, noexec, xmm2, xmm5, xmm4, rax); \
	xmm2 *= xmm7;		/* value *= grp two-to-phi */ \
	xmm2 *= xptr2(g->u.xmm.XMM_TTP_FUDGE,rax);/* value *= fudge two-to-phi */ \
	xmm2 *= xptr(rbx+4*32+16);	/* value *= col two-to-phi */ \
	xptr(rbp+18*16) = xmm2;	/* Save FFT data */ \
\
	rax = u8ptr(rdi+rcx+6); \
	xmm5 -= XMM_BIGVAL2; \
	xmm5 *= xmm7;		/* carry *= grp two-to-phi */ \
	xmm5 *= xptr2(g->u.xmm.XMM_TTP_FUDGE,rax);/* carry *= fudge two-to-phi */ \
	xmm5 *= xptr(rbx+5*32+16);	/* carry *= col two-to-phi */ \
	xmm5 += xptr(rbp+22*16);	/* Add carry and FFT data */ \
	xptr(rbp+22*16) = xmm5;	/* Save FFT data */ \
\
\
	rax = u8ptr(rdi+3);	/* Load big vs. little flag */ \
	xmm3 = xptr(rbp+3*16);	/* FFT data */ \
	xmm6 = xptr(rdx+3*32);	/* grp two-to-minus-phi */ \
	xmm6 *= g->u.xmm.XMM_NORM012_FF;	/* Mul by FFTLEN/2 */ \
	xmm3 *= xptr2(g->u.xmm.XMM_TTMP_FUDGE,rax);/* mul by fudge two-to-minus-phi */ \
	xmm3 *= xmm6;		/* mul by grp two-to-minus-phi */ \
	xmm3 *= xptr(rbx+0*32);	/* mul by col two-to-minus-phi */ \
	xmm3 += xptr(rsi+3*16);	/* Add in the carry */ \
	rounding(exec, base2, noexec, noexec, xmm3, xmm5, xmm4, rax); \
	xmm7 = xptr(rdx+3*32+16);	/* grp two-to-phi */ \
	xmm3 *= xmm7;		/* value *= grp two-to-phi */ \
	xmm3 *= xptr2(g->u.xmm.XMM_TTP_FUDGE,rax);/* value *= fudge two-to-phi */ \
	xmm3 *= xptr(rbx+0*32+16);	/* value *= col two-to-phi */ \
	xptr(rbp+3*16) = xmm3;	/* Save FFT data */ \
\
	rax = u8ptr(rdi+7);	/* Load big vs. little flag */ \
	xmm3 = xptr(rbp+7*16);	/* FFT data */ \
	xmm3 *= xmm6;		/* mul by grp two-to-minus-phi */ \
	xmm3 *= xptr2(g->u.xmm.XMM_TTMP_FUDGE,rax);/* mul by fudge two-to-minus-phi */ \
	xmm3 *= xptr(rbx+1*32);	/* mul by col two-to-minus-phi */ \
	xmm3 += xmm5;		/* Add in the carry */ \
	rounding(exec, base2, noexec, noexec, xmm3, xmm5, xmm4, rax); \
	xmm3 *= xmm7;		/* value *= grp two-to-phi */ \
	xmm3 *= xptr2(g->u.xmm.XMM_TTP_FUDGE,rax);/* value *= fudge two-to-phi */ \
	xmm3 *= xptr(rbx+1*32+16);	/* value *= col two-to-phi */ \
	xptr(rbp+7*16) = xmm3;	/* Save FFT data */ \
\
	rcx = g->BIGLIT_INCR2;/* Different clm values step through */ \
;					/* big/lit array differently */ \
	rax = u8ptr(rdi+rcx+3); /* Load big vs. little flag */ \
	xmm3 = xptr(rbp+11*16);	/* FFT data */ \
	xmm3 *= xmm6;		/* mul by grp two-to-minus-phi */ \
	xmm3 *= xptr2(g->u.xmm.XMM_TTMP_FUDGE,rax);/* mul by fudge two-to-minus-phi */ \
	xmm3 *= xptr(rbx+2*32);	/* mul by col two-to-minus-phi */ \
	xmm3 += xmm5;		/* Add in the carry */ \
	rounding(exec, base2, noexec, noexec, xmm3, xmm5, xmm4, rax); \
	xmm3 *= xmm7;		/* value *= grp two-to-phi */ \
	xmm3 *= xptr2(g->u.xmm.XMM_TTP_FUDGE,rax);/* value *= fudge two-to-phi */ \
	xmm3 *= xptr(rbx+2*32+16);	/* value *= col two-to-phi */ \
	xptr(rbp+11*16) = xmm3;	/* Save FFT data */ \
\
	rax = u8ptr(rdi+rcx+7); /* Load big vs. little flag */ \
	xmm3 = xptr(rbp+15*16);	/* FFT data */ \
	xmm3 *= xmm6;		/* mul by grp two-to-minus-phi */ \
	xmm3 *= xptr2(g->u.xmm.XMM_TTMP_FUDGE,rax);/* mul by fudge two-to-minus-phi */ \
	xmm3 *= xptr(rbx+3*32);	/* mul by col two-to-minus-phi */ \
	xmm3 += xmm5;		/* Add in the carry */ \
	rounding(exec, base2, noexec, noexec, xmm3, xmm5, xmm4, rax); \
	xmm3 *= xmm7;		/* value *= grp two-to-phi */ \
	xmm3 *= xptr2(g->u.xmm.XMM_TTP_FUDGE,rax);/* value *= fudge two-to-phi */ \
	xmm3 *= xptr(rbx+3*32+16);	/* value *= col two-to-phi */ \
	xptr(rbp+15*16) = xmm3;	/* Save FFT data */ \
\
	rcx = g->BIGLIT_INCR4;/* Different clm values step through */ \
;					/* big/lit array differently */ \
	rax = u8ptr(rdi+rcx+3); /* Load big vs. little flag */ \
	xmm3 = xptr(rbp+19*16);	/* FFT data */ \
	xmm3 *= xmm6;		/* mul by grp two-to-minus-phi */ \
	xmm3 *= xptr2(g->u.xmm.XMM_TTMP_FUDGE,rax);/* mul by fudge two-to-minus-phi */ \
	xmm3 *= xptr(rbx+4*32);	/* mul by col two-to-minus-phi */ \
	xmm3 += xmm5;		/* Add in the carry */ \
	rounding(exec, base2, noexec, noexec, xmm3, xmm5, xmm4, rax); \
	xmm3 *= xmm7;		/* value *= grp two-to-phi */ \
	xmm3 *= xptr2(g->u.xmm.XMM_TTP_FUDGE,rax);/* value *= fudge two-to-phi */ \
	xmm3 *= xptr(rbx+4*32+16);	/* value *= col two-to-phi */ \
	xptr(rbp+19*16) = xmm3;	/* Save FFT data */ \
\
	rax = u8ptr(rdi+rcx+7); \
	xmm5 -= XMM_BIGVAL2; \
	xmm5 *= xmm7;		/* carry *= grp two-to-phi */ \
	xmm5 *= xptr2(g->u.xmm.XMM_TTP_FUDGE,rax);/* carry *= fudge two-to-phi */ \
	xmm5 *= xptr(rbx+5*32+16);	/* carry *= col two-to-phi */ \
	xmm5 += xptr(rbp+23*16);	/* Add carry and FFT data */ \
	xptr(rbp+23*16) = xmm5;	/* Save FFT data */ \
} \
	xmm4 = XMM_BIGVAL2; \
	xptr(rsi+0*16) = xmm4;	/* Clear carry */ \
	xptr(rsi+1*16) = xmm4; \
	xptr(rsi+2*16) = xmm4; \
	xptr(rsi+3*16) = xmm4; \


/* Significantly different cleanup code for zero-padded FFTs. */ \
/* Note: The group multiplier should be 1.0 for the bottom FFT words and */ \
/* the FFT words just above the half-way point. */ \
\

#define xnorm012_2d_zpad_part1 \
\
if(g->this_block == 0) {		/* Are we carrying into the first data block? */ \
	/*jne	done;*/			/* If not, skip wrapping carries and negating the last carry */ \
	\
	rax = g->count3;		/* Load 3 section counts */ \
	\
	do { \
		rdi = rsi;		/* Save section start */ \
		rbx = rax;		/* Form count for this section */ \
		rbx &= 0x7FF; \
		if(!rbx) break;			/* No cache lines to do.  We're all done! */ \
		rax >>= 11;			/* Move counts list along */ \
		do { \
			xmm4[0] = f64ptr(rsi+0*16);	/* Load low carry word */ \
			f64ptr(rsi+0*16) = xmm0[0];	/* Save prev cache line's high carry in low word */ \
			xmm0[0] = f64ptr(rsi+0*16+8);	/* Load high carry for next cache line */ \
			f64ptr(rsi+0*16+8) = xmm4[0];	/* Move this cache line's low carry to high word */ \
			xmm4[0] = f64ptr(rsi+1*16);	/* Load low carry word */ \
			f64ptr(rsi+1*16) = xmm2[0];	/* Save prev cache line's high carry in low word */ \
			xmm2[0] = f64ptr(rsi+1*16+8);	/* Load high carry for next cache line */ \
			f64ptr(rsi+1*16+8) = xmm4[0];	/* Move this cache line's low carry to high word */ \
			xmm4[0] = f64ptr(rsi+2*16);	/* Load low carry word */ \
			f64ptr(rsi+2*16) = xmm1[0];	/* Save prev cache line's high carry in low word */ \
			xmm1[0] = f64ptr(rsi+2*16+8);	/* Load high carry for next cache line */ \
			f64ptr(rsi+2*16+8) = xmm4[0];	/* Move this cache line's low carry to high word */ \
			xmm4[0] = f64ptr(rsi+3*16);	/* Load low carry word */ \
			f64ptr(rsi+3*16) = xmm3[0];	/* Save prev cache line's high carry in low word */ \
			xmm3[0] = f64ptr(rsi+3*16+8);	/* Load high carry for next cache line */ \
			f64ptr(rsi+3*16+8) = xmm4[0];	/* Move this cache line's low carry to high word */ \
			rsi += 64;			/* Next carry cache line */ \
			rbx--;			/* Test loop counter */ \
		}while(rbx);		/* Next carry row in section */ \
	\
		f64ptr(rdi+1*16) = xmm0[0];	/* Do extra shuffling of section's last cache line */ \
		xmm0[0] = xmm2[0]; \
		f64ptr(rdi+3*16) = xmm1[0]; \
		xmm1[0] = xmm3[0]; \
	}while(1);			/* Next section */ \
\
	rsi = (uintptr_t)g->DESTARG;/* Load FFT data pointer */ \
	rdi = (uintptr_t)g->norm_biglit_array;	/* Addr of the big/little flags array */ \
	rbp = (uintptr_t)g->norm_col_mults;	/* Addr of the group multipliers */ \
	if(g->B_IS_2 == 0) {		/* Is b = 2? */ \
		xnorm012_2d_zpad_part1a(noexec); /* No, do harder rounding */ \
	} else { \
		xnorm012_2d_zpad_part1a(exec); /* Yes, do simpler rounding */ \
	} \
	rsi = (uintptr_t)g->carries;		/* Reload carries pointer */ \
	xmm7[0] = XMM_BIGVAL1;/* Clear two carries just processed by xnorm012_2d_zpad_part1_cmn */ \
	xmm6[0] = 0.0; \
	f64ptr(rsi) = xmm7[0]; \
	f64ptr(rsi+32) = xmm6[0]; \
}



#define xnorm012_2d_zpad_part1a(base2) \
	if(g->const_fft != 0) {		/* Are we also multiplying by a constant? */ \
		xnorm012_2d_zpad_part1_cmn(exec, base2); \
	}else{ \
		xnorm012_2d_zpad_part1_cmn(noexec, base2); \
	} \



#define xnorm012_2d_zpad_part1_cmn(const1, base2) { \
\
;	/* Strip BIGVAL from the traditional carry, we'll add the traditional */ \
;	/* carry in later when we are working on the ZPAD0 - ZPAD6 values. */ \
	xmm2[0] -= XMM_BIGVAL1;/* Integerize traditional carry */ \
\
;	/* Rather than calculate high FFT carry times k and then later dividing */ \
;	/* by k, we multiply FFT high carry by const and we'll add it */ \
;	/* to the lower FFT data later (after multiplying by -c). */ \
const1(xmm3[0] *= g->u.xmm.XMM_MULCONST[0]); \
\
;	/* Work on zero-pad addin value. */ \
	xmm4[0] = g->ADDIN_VALUE;		/* Load the add in value */ \
const1(xmm4[0] *= g->u.xmm.XMM_MULCONST[0]);		/* Multiply the add in value by the small mul const */ \
	xmm4[0] += g->POSTADDIN_VALUE;		/* Add the post-mul-by-const addin value */ \
;	/* when c = -1, 1 = b^n */ \
;	/* when c = 1, -1 = b^n, 1 = -b^n */ \
;	/* when c = -3, 3 = b^n, 1 = b^n - 2 */ \
;	/* when c = 3, -3 = b^n, 3 = -b^n, 1 = -b^n - 2 */ \
;	/* The "- 2" has been precomputed in ZPAD_LSW_ADJUST.  Add ADDIN_VALUE * mul-by-const * ZPAD_LSW_ADJUST into the least significant FFT word */ \
	rsi = (uintptr_t)g->DESTARG;/* Address of squared number */ \
	xmm0[0] = g->ZPAD_LSW_ADJUST; \
	xmm0[0] *= xmm4[0]; \
	xmm0[0] += f64ptr(rsi); \
	f64ptr(rsi) = xmm0[0]; \
\
;	/* Multiply ZPAD0 through ZPAD6 by const * -C.  This, in essense, */ \
;	/* wraps this data from above the FFT data area to the halfway point. */ \
;	/* Later on we'll divide this by K to decide which data needs wrapping */ \
;	/* all the way down to the bottom of the FFT data. */ \
\
;	/* NOTE that ZPAD0's column multiplier is 1.0.  Also, ZPAD6 will not */ \
;	/* be bigger than a big word.  We must be careful to handle c's up */ \
;	/* to about 30 bits */ \
\
	rax = u8ptr(rdi);	/* Load big vs. little flags */ \
	xmm0[0] = g->ZPAD0_6[0];/* Load values1 */ \
	xmm5[0] = 0.0;		/* Create a zero high FFT carry to add in */ \
	single_split_lower_zpad_word(base2, xmm0[0], xmm5[0], xmm1[0], rax); \
no##const1(xmm0[0] *= g->u.xmm.XMM_MINUS_C[0]); \
const1(xmm0[0] *= g->u.xmm.XMM_MINUS_C_TIMES_MULCONST[0]); \
	xmm0[0] += xmm4[0];		/* Apply rest of ADDIN_VALUE * mul-by-const here (the FFT half way point) */ \
	round_zpad7_word(base2, xmm0[0], xmm2[0], xmm1[0], rax); \
	g->ZPAD0_6[0] = xmm0[0]; \
\
	rax = u8ptr(rdi+4);	/* Load big vs. little flags */ \
	xmm0[0] = g->ZPAD0_6[1];/* Load values1 */ \
	xmm0[0] *= f64ptr(rbp+32);	/* Mul values1 by two-to-minus-phi */ \
	xmm0[0] *= g->u.xmm.XMM_NORM012_FF[0];	/* Mul by FFTLEN/2 */ \
	single_split_lower_zpad_word(base2, xmm0[0], xmm5[0], xmm1[0], rax); \
no##const1(xmm0[0] *= g->u.xmm.XMM_MINUS_C[0]); \
const1(xmm0[0] *= g->u.xmm.XMM_MINUS_C_TIMES_MULCONST[0]); \
	round_zpad7_word(base2, xmm0[0], xmm2[0], xmm1[0], rax); \
	g->ZPAD0_6[1] = xmm0[0]; \
\
	rcx = g->BIGLIT_INCR2;/* Different clm values step through */ \
;					/* big/lit array differently */ \
	rax = u8ptr(rdi+rcx);	/* Load big vs. little flags */ \
	xmm0[0] = g->ZPAD0_6[2];/* Load values1 */ \
	xmm0[0] *= f64ptr(rbp+2*32);	/* Mul values1 by two-to-minus-phi */ \
	xmm0[0] *= g->u.xmm.XMM_NORM012_FF[0];	/* Mul by FFTLEN/2 */ \
	single_split_lower_zpad_word(base2, xmm0[0], xmm5[0], xmm1[0], rax); \
no##const1(xmm0[0] *= g->u.xmm.XMM_MINUS_C[0]); \
const1(xmm0[0] *= g->u.xmm.XMM_MINUS_C_TIMES_MULCONST[0]); \
	round_zpad7_word(base2, xmm0[0], xmm2[0], xmm1[0], rax); \
	g->ZPAD0_6[2] = xmm0[0]; \
\
	rax = u8ptr(rdi+rcx+4); /* Load big vs. little flags */ \
	xmm0[0] = g->ZPAD0_6[3];/* Load values1 */ \
	xmm0[0] *= f64ptr(rbp+3*32);	/* Mul values1 by two-to-minus-phi */ \
	xmm0[0] *= g->u.xmm.XMM_NORM012_FF[0];	/* Mul by FFTLEN/2 */ \
	single_split_lower_zpad_word(base2, xmm0[0], xmm5[0], xmm1[0], rax); \
no##const1(xmm0[0] *= g->u.xmm.XMM_MINUS_C[0]); \
const1(xmm0[0] *= g->u.xmm.XMM_MINUS_C_TIMES_MULCONST[0]); \
	round_zpad7_word(base2, xmm0[0], xmm2[0], xmm1[0], rax); \
	g->ZPAD0_6[3] = xmm0[0]; \
\
	rcx = g->BIGLIT_INCR4;/* Different clm values step through */ \
;					/* big/lit array differently */ \
	rax = u8ptr(rdi+rcx);	/* Load big vs. little flags */ \
	xmm0[0] = g->ZPAD0_6[4];/* Load values1 */ \
	xmm0[0] *= f64ptr(rbp+4*32);	/* Mul values1 by two-to-minus-phi */ \
	xmm0[0] *= g->u.xmm.XMM_NORM012_FF[0];	/* Mul by FFTLEN/2 */ \
	single_split_lower_zpad_word(base2, xmm0[0], xmm5[0], xmm1[0], rax); \
no##const1(xmm0[0] *= g->u.xmm.XMM_MINUS_C[0]); \
const1(xmm0[0] *= g->u.xmm.XMM_MINUS_C_TIMES_MULCONST[0]); \
	round_zpad7_word(base2, xmm0[0], xmm2[0], xmm1[0], rax); \
	g->ZPAD0_6[4] = xmm0[0]; \
\
	rax = u8ptr(rdi+rcx+4); /* Load big vs. little flags */ \
	xmm0[0] = g->ZPAD0_6[5];/* Load values1 */ \
	xmm0[0] *= f64ptr(rbp+5*32);	/* Mul values1 by two-to-minus-phi */ \
	xmm0[0] *= g->u.xmm.XMM_NORM012_FF[0];	/* Mul by FFTLEN/2 */ \
	single_split_lower_zpad_word(base2, xmm0[0], xmm5[0], xmm1[0], rax); \
no##const1(xmm0[0] *= g->u.xmm.XMM_MINUS_C[0]); \
const1(xmm0[0] *= g->u.xmm.XMM_MINUS_C_TIMES_MULCONST[0]); \
	round_zpad7_word(base2, xmm0[0], xmm2[0], xmm1[0], rax); \
	g->ZPAD0_6[5] = xmm0[0]; \
\
	xmm0[0] = g->ZPAD0_6[6];/* Load values1 */ \
	xmm0[0] *= f64ptr(rbp+6*32);	/* Mul values1 by two-to-minus-phi */ \
	xmm0[0] *= g->u.xmm.XMM_NORM012_FF[0];	/* Mul by FFTLEN/2 */ \
	xmm0[0] = roundf64(xmm0[0]); /* Round to an integer */ \
	xmm0[0] += xmm5[0];		/* Add in shifted high ZPAD data */ \
no##const1(xmm0[0] *= g->u.xmm.XMM_MINUS_C[0]); \
const1(xmm0[0] *= g->u.xmm.XMM_MINUS_C_TIMES_MULCONST[0]); \
	xmm0[0] += xmm2[0];		/* Add in high part of last calculation */ \
	g->ZPAD0_6[6] = xmm0[0]; \
\
;	/* Divide the zpad data by k.  Store the integer part in XMM_TMP */ \
;	/* and the remainder in ZPAD0.  Later we will wrap the integer part */ \
;	/* down to the bottom of the FFT data area (and multiply by -c). */ \
;	/* And we will store the remainder in the upper half of the FFT */ \
;	/* data area. */ \
\
;	/* Note there are three cases to handle.  K is smaller than a big word. */ \
;	/* K is between one and 2 big words in size.  And K is more than */ \
;	/* 2 big words in size. */ \
\
	xfunc572(g); \
\
;	/* Now normalize the data above the halfway point.  Remember that the */ \
;	/* column two-to-phi multiplier for the first value will be 1.0. */ \
\
	rax = u8ptr(rdi);	/* First word  */ \
	xmm0[0] = g->ZPAD0_6[0];/* Load remainder of divide by k */ \
	xmm0[0] += XMM_BIGVAL1; \
	single_rounding(base2, xmm0[0], xmm2[0], xmm4[0], rax); \
	f64ptr(rsi+0*64+32) = xmm0[0];	/* Save value1 */ \
\
	rax = u8ptr(rdi+4);	/* Load big vs. little flags */ \
	single_rounding(base2, xmm2[0], xmm0[0], xmm4[0], rax); \
	xmm2[0] *= f64ptr(rbp+1*32+16);	/* new value2 = val * two-to-phi */ \
	f64ptr(rsi+1*64+32) = xmm2[0];	/* Save value2 */ \
\
	rcx = g->BIGLIT_INCR2;/* Different clm values step through */ \
;					/* big/lit array differently */ \
	rax = u8ptr(rdi+rcx);	/* Load big vs. little flags */ \
	single_rounding(base2, xmm0[0], xmm2[0], xmm4[0], rax); \
	xmm0[0] *= f64ptr(rbp+2*32+16);	/* new value3 = val * two-to-phi */ \
	f64ptr(rsi+2*64+32) = xmm0[0];	/* Save value3 */ \
\
	xmm2[0] -= XMM_BIGVAL1;/* Remove integer rounding constant */ \
	xmm2[0] *= f64ptr(rbp+3*32+16);	/* value4 = carry * two-to-phi */ \
	f64ptr(rsi+3*64+32) = xmm2[0];	/* Save new value4 */ \
\
;	/* Mul the integer part of (ZPAD data divided by k) by -c in */ \
;	/* preparation for adding it into the lower FFT data area. */ \
;	/* Also add in the shifted high FFT carry at this time. */ \
\
;	/* Now add in and normalize the bottom FFT data.  Remember that the */ \
;	/* column two-to-phi multiplier for the first value will be 1.0.  We  */ \
;	/* must go 6 words deep in case k is 48-50 bits and c is 32 bits. */ \
\
	vec2f64 TMP12 = g->u.xmm.XMM_TMP1; \
	vec2f64 TMP34 = g->u.xmm.XMM_TMP2; \
	vec2f64 TMP56 = g->u.xmm.XMM_TMP3; \
	TMP34 *= g->u.xmm.XMM_MINUS_C;	/* Mul by -c */ \
	TMP56 *= g->u.xmm.XMM_MINUS_C;	/* Mul by -c */ \
	rax = u8ptr(rdi);	/* First word  */ \
	xmm0[0] = TMP12[0];	/* Load integer part of divide by k */ \
	xmm0[0] += xmm3[0];		/* Add in shifted high FFT carry */ \
	xmm0[0] *= g->u.xmm.XMM_MINUS_C[0];	/* Mul by -c */ \
	xmm0[0] += XMM_BIGVAL1; \
	xmm0[0] += f64ptr(rsi+0*64);	/* Add in the FFT data */ \
	single_rounding(base2, xmm0[0], xmm2[0], xmm4[0], rax); \
	f64ptr(rsi+0*64) = xmm0[0];	/* Save value1 */ \
\
	rax = u8ptr(rdi+4);	/* Load big vs. little flags */ \
	xmm0[0] = TMP12[1];	/* Load integer part of divide by k */ \
	xmm0[0] *= g->u.xmm.XMM_MINUS_C[0];	/* Mul by -c */ \
	xmm1[0] = f64ptr(rsi+1*64);	/* Load FFT data */ \
	xmm1[0] *= f64ptr(rbp+1*32);	/* Mul values2 by two-to-minus-phi */ \
	xmm1[0] *= g->u.xmm.XMM_NORM012_FF[0];	/* Mul by FFTLEN/2 */ \
	xmm0[0] += xmm1[0];		/* Add in the FFT data */ \
	xmm0[0] += xmm2[0];		/* x2 = value + carry */ \
	single_rounding(base2, xmm0[0], xmm2[0], xmm4[0], rax); \
	xmm0[0] *= f64ptr(rbp+1*32+16);	/* new value2 = val * two-to-phi */ \
	f64ptr(rsi+1*64) = xmm0[0];	/* Save value2 */ \
\
	rcx = rdi + g->BIGLIT_INCR2;/* Different clm values step through */ \
;					/* big/lit array differently */ \
	rax = u8ptr(rcx); /* Load big vs. little flags */ \
	xmm0[0] = TMP34[0];	/* Load integer part of divide by k */ \
	xmm1[0] = f64ptr(rsi+2*64);	/* Load FFT data */ \
	xmm1[0] *= f64ptr(rbp+2*32);	/* Mul values3 by two-to-minus-phi */ \
	xmm1[0] *= g->u.xmm.XMM_NORM012_FF[0];	/* Mul by FFTLEN/2 */ \
	xmm0[0] += xmm1[0];		/* Add in the FFT data */ \
	xmm0[0] += xmm2[0];		/* x3 = value + carry */ \
	single_rounding(base2, xmm0[0], xmm2[0], xmm4[0], rax); \
	xmm0[0] *= f64ptr(rbp+2*32+16);	/* new value3 = val * two-to-phi */ \
	f64ptr(rsi+2*64) = xmm0[0];	/* Save value3 */ \
\
	rax = u8ptr(rcx+4); /* Load big vs. little flags */ \
	xmm0[0] = TMP34[1];	/* Load integer part of divide by k */ \
	xmm1[0] = f64ptr(rsi+3*64);	/* Load FFT data */ \
	xmm1[0] *= f64ptr(rbp+3*32);	/* Mul values4 by two-to-minus-phi */ \
	xmm1[0] *= g->u.xmm.XMM_NORM012_FF[0];	/* Mul by FFTLEN/2 */ \
	xmm0[0] += xmm1[0];		/* Add in the FFT data */ \
	xmm0[0] += xmm2[0];		/* x4 = value + carry */ \
	single_rounding(base2, xmm0[0], xmm2[0], xmm4[0], rax); \
	xmm0[0] *= f64ptr(rbp+3*32+16);	/* new value4 = val * two-to-phi */ \
	f64ptr(rsi+3*64) = xmm0[0];	/* Save value4 */ \
\
	rcx = rdi + g->BIGLIT_INCR4;/* Different clm values step through */ \
;					/* big/lit array differently */ \
	rax = u8ptr(rcx);	/* Load big vs. little flags */ \
	xmm0[0] = TMP56[0];	/* Load integer part of divide by k */ \
	xmm1[0] = f64ptr(rsi+4*64);	/* Load FFT data */ \
	xmm1[0] *= f64ptr(rbp+4*32);	/* Mul values4 by two-to-minus-phi */ \
	xmm1[0] *= g->u.xmm.XMM_NORM012_FF[0];	/* Mul by FFTLEN/2 */ \
	xmm0[0] += xmm1[0];		/* Add in the FFT data */ \
	xmm0[0] += xmm2[0];		/* x4 = value + carry */ \
	single_rounding(base2, xmm0[0], xmm2[0], xmm4[0], rax); \
	xmm0[0] *= f64ptr(rbp+4*32+16);	/* new value4 = val * two-to-phi */ \
	f64ptr(rsi+4*64) = xmm0[0];	/* Save value4 */ \
\
	rax = u8ptr(rcx+4); /* Load big vs. little flags */ \
	xmm0[0] = TMP56[1];	/* Load integer part of divide by k */ \
	xmm1[0] = f64ptr(rsi+5*64);	/* Load FFT data */ \
	xmm1[0] *= f64ptr(rbp+5*32);	/* Mul values5 by two-to-minus-phi */ \
	xmm1[0] *= g->u.xmm.XMM_NORM012_FF[0];	/* Mul by FFTLEN/2 */ \
	xmm0[0] += xmm1[0];		/* Add in the FFT data */ \
	xmm0[0] += xmm2[0];		/* x5 = value + carry */ \
	single_rounding(base2, xmm0[0], xmm2[0], xmm4[0], rax); \
	xmm0[0] *= f64ptr(rbp+5*32+16);	/* new value5 = val * two-to-phi */ \
	f64ptr(rsi+5*64) = xmm0[0];	/* Save value5 */ \
\
	rcx += g->BIGLIT_INCR2;/* Different clm values step through */ \
;					/* big/lit array differently */ \
	rax = u8ptr(rcx);	/* Load big vs. little flags */ \
	xmm0[0] = f64ptr(rsi+6*64);	/* Load FFT data */ \
	xmm0[0] *= f64ptr(rbp+6*32);	/* Mul values3 by two-to-minus-phi */ \
	xmm0[0] *= g->u.xmm.XMM_NORM012_FF[0];	/* Mul by FFTLEN/2 */ \
	xmm0[0] += xmm2[0];		/* x6 = value + carry */ \
	single_rounding(base2, xmm0[0], xmm2[0], xmm4[0], rax); \
	xmm0[0] *= f64ptr(rbp+6*32+16);	/* new value6 = val * two-to-phi */ \
	f64ptr(rsi+6*64) = xmm0[0];	/* Save value6 */ \
\
	xmm2[0] -= XMM_BIGVAL1;/* Remove rounding constant */ \
	xmm2[0] *= f64ptr(rbp+7*32+16);	/* new value7 = val * two-to-phi */ \
	xmm2[0] += f64ptr(rsi+7*64);	/* Add in FFT data */ \
	f64ptr(rsi+7*64) = xmm2[0];	/* Save value7 */ \
}


#define xnorm012_2d_zpad \
	if(g->B_IS_2 == 0) {		/* Is b = 2? */ \
		xnorm012_2d_zpada(noexec); /* No, do harder rounding */ \
	}else{ \
		xnorm012_2d_zpada(exec); /* Yes, do simpler rounding */ \
	}

#define xnorm012_2d_zpada(base2) \
	if(g->const_fft != 0) { /* Are we also multiplying by a constant? */ \
		xnorm012_2d_zpad_cmn(exec, base2); \
	}else{ \
		xnorm012_2d_zpad_cmn(noexec, base2); \
	}

#define xnorm012_2d_zpad_cmn(const1, base2) \
	rax = u8ptr(rdi+0);	/* Load big vs. little flag */ \
	xmm0 = xptr(rsi+0*16);	/* Load carry */ \
	xmm2 = xptr(rsi+2*16);	/* Load FFT hi data carry */ \
	xmm4 = xptr(rbp+0*64);	/* Load FFT data */ \
	xmm4 *= xptr(rdx+0*32);	/* mul by grp two-to-minus-phi */ \
	xmm4 *= xptr(rbx+0*32);	/* mul by col two-to-minus-phi */ \
	xmm4 *= xptr2(g->u.xmm.XMM_TTMP_FUDGE,rax);/* mul by fudge two-to-minus-phi */ \
	xmm4 *= g->u.xmm.XMM_NORM012_FF;	/* Mul by FFTLEN/2 */ \
	xmm0 += xmm4;		/* x1 = values1 + carry */ \
	split_carry_zpad_word(base2, xmm2, xmm6, xmm4, rax); \
no##const1(xmm4 = g->u.xmm.XMM_K_LO);		/* Calc high FFT carry times k */ \
const1(xmm4 = g->u.xmm.XMM_K_TIMES_MULCONST_LO); \
	xmm4 *= xmm2;		/* high_FFT_carry * k_lo */ \
	xmm0 += xmm4;		/* x1 = x1 + high_FFT_carry * k_lo */ \
no##const1(xmm5 = g->u.xmm.XMM_K_HI); \
const1(xmm5 = g->u.xmm.XMM_K_TIMES_MULCONST_HI); \
	xmm5 *= xptr2(g->u.xmm.XMM_LIMIT_INVERSE,rax); /* shift k_hi */ \
	xmm5 *= xmm2; \
	rounding(exec, base2, noexec, noexec, xmm0, xmm2, xmm4, rax); \
	xmm2 += xmm5;		/* Carry += high_FFT_carry * shifted k_hi */ \
	xmm0 *= xptr(rdx+0*32+16);	/* mul by grp two-to-phi */ \
	xmm0 *= xptr(rbx+0*32+16);	/* mul by col two-to-phi */ \
	xmm0 *= xptr2(g->u.xmm.XMM_TTP_FUDGE,rax);/* mul by fudge two-to-phi */ \
	xptr(rbp+0*64) = xmm0;	/* Save FFT data */ \
\
	rax = u8ptr(rdi+4);	/* Load big vs. little flag */ \
	xmm0 = xptr(rbp+1*64);	/* Load FFT data */ \
	xmm0 *= xptr(rdx+0*32);	/* mul by grp two-to-minus-phi */ \
	xmm4 = xptr(rbx+1*32);	/* col two-to-minus-phi */ \
	xmm4 *= xptr2(g->u.xmm.XMM_TTMP_FUDGE,rax);/* mul by fudge two-to-minus-phi */ \
	xmm0 *= xmm4;		/* data *= fudged col two-to-minus-phi */ \
	xmm0 *= g->u.xmm.XMM_NORM012_FF;	/* Mul by FFTLEN/2 */ \
	xmm0 += xmm2;		/* x2 = values1 + carry */ \
	split_carry_zpad_word(base2, xmm6, xmm5, xmm4, rax); \
no##const1(xmm4 = g->u.xmm.XMM_K_LO);		/* Calc high FFT carry times k */ \
const1(xmm4 = g->u.xmm.XMM_K_TIMES_MULCONST_LO); \
	xmm4 *= xmm6;		/* high_FFT_carry * k_lo */ \
	xmm0 += xmm4;		/* x2 = x2 + high_FFT_carry * k_lo */ \
no##const1(xmm4 = g->u.xmm.XMM_K_HI); \
const1(xmm4 = g->u.xmm.XMM_K_TIMES_MULCONST_HI); \
	xmm4 *= xptr2(g->u.xmm.XMM_LIMIT_INVERSE,rax); /* shift k_hi */ \
	xmm6 *= xmm4; \
	rounding(exec, base2, noexec, noexec, xmm0, xmm2, xmm4, rax); \
	xmm2 += xmm6;		/* Carry += high_FFT_carry * shifted k_hi */ \
	xmm0 *= xptr(rdx+0*32+16);	/* mul by grp two-to-phi */ \
	xmm4 = xptr(rbx+1*32+16);	/* col two-to-phi */ \
	xmm4 *= xptr2(g->u.xmm.XMM_TTP_FUDGE,rax);/* mul by fudge two-to-phi */ \
	xmm0 *= xmm4;		/* data *= fudged col two-to-phi */ \
	xptr(rbp+1*64) = xmm0;	/* Save FFT data */ \
\
	rcx = g->BIGLIT_INCR2;/* Different clm values step through */ \
;					/* big/lit array differently */ \
	rax = u8ptr(rdi+rcx);	/* Load big vs. little flag */ \
	xmm0 = xptr(rbp+2*64);	/* Load FFT data */ \
	xmm0 *= xptr(rdx+0*32);	/* mul by grp two-to-minus-phi */ \
	xmm4 = xptr(rbx+2*32);	/* col two-to-minus-phi */ \
	xmm4 *= xptr2(g->u.xmm.XMM_TTMP_FUDGE,rax);/* mul by fudge two-to-minus-phi */ \
	xmm0 *= xmm4;		/* data *= fudged col two-to-minus-phi */ \
	xmm0 *= g->u.xmm.XMM_NORM012_FF;	/* Mul by FFTLEN/2 */ \
	xmm0 += xmm2;		/* x3 = values1 + carry */ \
	split_carry_zpad_word(base2, xmm5, xmm6, xmm4, rax); \
no##const1(xmm4 = g->u.xmm.XMM_K_LO);		/* Calc high FFT carry times k */ \
const1(xmm4 = g->u.xmm.XMM_K_TIMES_MULCONST_LO); \
	xmm4 *= xmm5;		/* high_FFT_carry * k_lo */ \
	xmm0 += xmm4;		/* x3 = x3 + high_FFT_carry * k_lo */ \
no##const1(xmm4 = g->u.xmm.XMM_K_HI); \
const1(xmm4 = g->u.xmm.XMM_K_TIMES_MULCONST_HI); \
	xmm4 *= xptr2(g->u.xmm.XMM_LIMIT_INVERSE,rax); /* shift k_hi */ \
	xmm5 *= xmm4; \
	rounding(exec, base2, noexec, noexec, xmm0, xmm2, xmm4, rax); \
	xmm2 += xmm5;		/* Carry += high_FFT_carry * shifted k_hi */ \
	xmm0 *= xptr(rdx+0*32+16);	/* mul by grp two-to-phi */ \
	xmm4 = xptr(rbx+2*32+16);	/* col two-to-phi */ \
	xmm4 *= xptr2(g->u.xmm.XMM_TTP_FUDGE,rax);/* mul by fudge two-to-phi */ \
	xmm0 *= xmm4;		/* data *= fudged col two-to-phi */ \
	xptr(rbp+2*64) = xmm0;	/* Save FFT data */ \
\
	rax = u8ptr(rdi+rcx+4); /* Load big vs. little flag */ \
	xmm0 = xptr(rbp+3*64);	/* Load FFT data */ \
	xmm0 *= xptr(rdx+0*32);	/* mul by grp two-to-minus-phi */ \
	xmm4 = xptr(rbx+3*32);	/* col two-to-minus-phi */ \
	xmm4 *= xptr2(g->u.xmm.XMM_TTMP_FUDGE,rax);/* mul by fudge two-to-minus-phi */ \
	xmm0 *= xmm4;		/* data *= fudged col two-to-minus-phi */ \
	xmm0 *= g->u.xmm.XMM_NORM012_FF;	/* Mul by FFTLEN/2 */ \
	xmm0 += xmm2;		/* x4 = values1 + carry */ \
	split_carry_zpad_word(base2, xmm6, xmm5, xmm4, rax); \
no##const1(xmm4 = g->u.xmm.XMM_K_LO);		/* Calc high FFT carry times k */ \
const1(xmm4 = g->u.xmm.XMM_K_TIMES_MULCONST_LO); \
	xmm4 *= xmm6;		/* high_FFT_carry * k_lo */ \
	xmm0 += xmm4;		/* x4 = x4 + high_FFT_carry * k_lo */ \
no##const1(xmm4 = g->u.xmm.XMM_K_HI); \
const1(xmm4 = g->u.xmm.XMM_K_TIMES_MULCONST_HI); \
	xmm4 *= xptr2(g->u.xmm.XMM_LIMIT_INVERSE,rax); /* shift k_hi */ \
	xmm6 *= xmm4; \
	rounding(exec, base2, noexec, noexec, xmm0, xmm2, xmm4, rax); \
	xmm2 += xmm6;		/* Carry += high_FFT_carry * shifted k_hi */ \
	xmm0 *= xptr(rdx+0*32+16);	/* mul by grp two-to-phi */ \
	xmm4 = xptr(rbx+3*32+16);	/* col two-to-phi */ \
	xmm4 *= xptr2(g->u.xmm.XMM_TTP_FUDGE,rax);/* mul by fudge two-to-phi */ \
	xmm0 *= xmm4;		/* data *= fudged col two-to-phi */ \
	xptr(rbp+3*64) = xmm0;	/* Save FFT data */ \
\
	rcx = g->BIGLIT_INCR4;/* Different clm values step through */ \
;					/* big/lit array differently */ \
	rax = u8ptr(rdi+rcx);	/* Load big vs. little flag */ \
	xmm0 = xptr(rbp+4*64);	/* Load FFT data */ \
	xmm0 *= xptr(rdx+0*32);	/* mul by grp two-to-minus-phi */ \
	xmm4 = xptr(rbx+4*32);	/* col two-to-minus-phi */ \
	xmm4 *= xptr2(g->u.xmm.XMM_TTMP_FUDGE,rax);/* mul by fudge two-to-minus-phi */ \
	xmm0 *= xmm4;		/* data *= fudged col two-to-minus-phi */ \
	xmm0 *= g->u.xmm.XMM_NORM012_FF;	/* Mul by FFTLEN/2 */ \
	xmm0 += xmm2;		/* x5 = values1 + carry */ \
no##const1(xmm4 = g->u.xmm.XMM_K_LO);		/* Calc high FFT carry times k */ \
const1(xmm4 = g->u.xmm.XMM_K_TIMES_MULCONST_LO); \
	xmm4 *= xmm5;		/* high_FFT_carry * k_lo */ \
	xmm0 += xmm4;		/* x5 = x5 + high_FFT_carry * k_lo */ \
no##const1(xmm4 = g->u.xmm.XMM_K_HI); \
const1(xmm4 = g->u.xmm.XMM_K_TIMES_MULCONST_HI); \
	xmm4 *= xptr2(g->u.xmm.XMM_LIMIT_INVERSE,rax); /* shift k_hi */ \
	xmm5 *= xmm4; \
	rounding(exec, base2, noexec, noexec, xmm0, xmm2, xmm4, rax); \
	xmm2 += xmm5;		/* Carry += high_FFT_carry * shifted k_hi */ \
	xmm0 *= xptr(rdx+0*32+16);	/* mul by grp two-to-phi */ \
	xmm4 = xptr(rbx+4*32+16);	/* col two-to-phi */ \
	xmm4 *= xptr2(g->u.xmm.XMM_TTP_FUDGE,rax);/* mul by fudge two-to-phi */ \
	xmm0 *= xmm4;		/* data *= fudged col two-to-phi */ \
	xptr(rbp+4*64) = xmm0;	/* Save FFT data */ \
\
	rax = u8ptr(rdi+rcx+4); /* Load big vs. little flag */ \
	xmm2 -= XMM_BIGVAL2;/* Remove rounding const from carry */ \
	xmm2 *= xptr(rdx+0*32+16);	/* mul by grp two-to-phi */ \
	xmm4 = xptr(rbx+5*32+16);	/* col two-to-phi */ \
	xmm4 *= xptr2(g->u.xmm.XMM_TTP_FUDGE,rax);/* mul by fudge two-to-phi */ \
	xmm2 *= xmm4;		/* data *= fudged col two-to-phi */ \
	xmm2 += xptr(rbp+5*64);	/* Load FFT data */ \
	xptr(rbp+5*64) = xmm2;	/* Save FFT data */ \
\
	rax = u8ptr(rdi+1);	/* Load big vs. little flag */ \
	xmm1 = xptr(rsi+1*16);	/* Load carry */ \
	xmm3 = xptr(rsi+3*16);	/* Load FFT hi data carry */ \
	xmm4 = xptr(rbp+0*64+16);	/* Load FFT data */ \
	xmm4 *= xptr(rdx+1*32);	/* mul by grp two-to-minus-phi */ \
	xmm4 *= xptr(rbx+0*32);	/* mul by col two-to-minus-phi */ \
	xmm4 *= xptr2(g->u.xmm.XMM_TTMP_FUDGE,rax);/* mul by fudge two-to-minus-phi */ \
	xmm4 *= g->u.xmm.XMM_NORM012_FF;	/* Mul by FFTLEN/2 */ \
	xmm1 += xmm4;		/* x1 = values1 + carry */ \
	split_carry_zpad_word(base2, xmm3, xmm6, xmm4, rax); \
no##const1(xmm4 = g->u.xmm.XMM_K_LO);		/* Calc high FFT carry times k */ \
const1(xmm4 = g->u.xmm.XMM_K_TIMES_MULCONST_LO); \
	xmm4 *= xmm3;		/* high_FFT_carry * k_lo */ \
	xmm1 += xmm4;		/* x1 = x1 + high_FFT_carry * k_lo */ \
no##const1(xmm5 = g->u.xmm.XMM_K_HI); \
const1(xmm5 = g->u.xmm.XMM_K_TIMES_MULCONST_HI); \
	xmm5 *= xptr2(g->u.xmm.XMM_LIMIT_INVERSE,rax); /* shift k_hi */ \
	xmm5 *= xmm3; \
	rounding(exec, base2, noexec, noexec, xmm1, xmm3, xmm4, rax); \
	xmm3 += xmm5;		/* Carry += high_FFT_carry * shifted k_hi */ \
	xmm1 *= xptr(rdx+1*32+16);	/* mul by grp two-to-phi */ \
	xmm1 *= xptr(rbx+0*32+16);	/* mul by col two-to-phi */ \
	xmm1 *= xptr2(g->u.xmm.XMM_TTP_FUDGE,rax);/* mul by fudge two-to-phi */ \
	xptr(rbp+0*64+16) = xmm1;	/* Save FFT data */ \
\
	rax = u8ptr(rdi+5);	/* Load big vs. little flag */ \
	xmm1 = xptr(rbp+1*64+16);	/* Load FFT data */ \
	xmm1 *= xptr(rdx+1*32);	/* mul by grp two-to-minus-phi */ \
	xmm4 = xptr(rbx+1*32);	/* col two-to-minus-phi */ \
	xmm4 *= xptr2(g->u.xmm.XMM_TTMP_FUDGE,rax);/* mul by fudge two-to-minus-phi */ \
	xmm1 *= xmm4;		/* data *= fudged col two-to-minus-phi */ \
	xmm1 *= g->u.xmm.XMM_NORM012_FF;	/* Mul by FFTLEN/2 */ \
	xmm1 += xmm3;		/* x2 = values1 + carry */ \
	split_carry_zpad_word(base2, xmm6, xmm5, xmm4, rax); \
no##const1(xmm4 = g->u.xmm.XMM_K_LO);		/* Calc high FFT carry times k */ \
const1(xmm4 = g->u.xmm.XMM_K_TIMES_MULCONST_LO); \
	xmm4 *= xmm6;		/* high_FFT_carry * k_lo */ \
	xmm1 += xmm4;		/* x2 = x2 + high_FFT_carry * k_lo */ \
no##const1(xmm4 = g->u.xmm.XMM_K_HI); \
const1(xmm4 = g->u.xmm.XMM_K_TIMES_MULCONST_HI); \
	xmm4 *= xptr2(g->u.xmm.XMM_LIMIT_INVERSE,rax); /* shift k_hi */ \
	xmm6 *= xmm4; \
	rounding(exec, base2, noexec, noexec, xmm1, xmm3, xmm4, rax); \
	xmm3 += xmm6;		/* Carry += high_FFT_carry * shifted k_hi */ \
	xmm1 *= xptr(rdx+1*32+16);	/* mul by grp two-to-phi */ \
	xmm4 = xptr(rbx+1*32+16);	/* col two-to-phi */ \
	xmm4 *= xptr2(g->u.xmm.XMM_TTP_FUDGE,rax);/* mul by fudge two-to-phi */ \
	xmm1 *= xmm4;		/* data *= fudged col two-to-phi */ \
	xptr(rbp+1*64+16) = xmm1;	/* Save FFT data */ \
\
	rcx = g->BIGLIT_INCR2;/* Different clm values step through */ \
;					/* big/lit array differently */ \
	rax = u8ptr(rdi+rcx+1); /* Load big vs. little flag */ \
	xmm1 = xptr(rbp+2*64+16);	/* Load FFT data */ \
	xmm1 *= xptr(rdx+1*32);	/* mul by grp two-to-minus-phi */ \
	xmm4 = xptr(rbx+2*32);	/* col two-to-minus-phi */ \
	xmm4 *= xptr2(g->u.xmm.XMM_TTMP_FUDGE,rax);/* mul by fudge two-to-minus-phi */ \
	xmm1 *= xmm4;		/* data *= fudged col two-to-minus-phi */ \
	xmm1 *= g->u.xmm.XMM_NORM012_FF;	/* Mul by FFTLEN/2 */ \
	xmm1 += xmm3;		/* x3 = values1 + carry */ \
	split_carry_zpad_word(base2, xmm5, xmm6, xmm4, rax); \
no##const1(xmm4 = g->u.xmm.XMM_K_LO);		/* Calc high FFT carry times k */ \
const1(xmm4 = g->u.xmm.XMM_K_TIMES_MULCONST_LO); \
	xmm4 *= xmm5;		/* high_FFT_carry * k_lo */ \
	xmm1 += xmm4;		/* x3 = x3 + high_FFT_carry * k_lo */ \
no##const1(xmm4 = g->u.xmm.XMM_K_HI); \
const1(xmm4 = g->u.xmm.XMM_K_TIMES_MULCONST_HI); \
	xmm4 *= xptr2(g->u.xmm.XMM_LIMIT_INVERSE,rax); /* shift k_hi */ \
	xmm5 *= xmm4; \
	rounding(exec, base2, noexec, noexec, xmm1, xmm3, xmm4, rax); \
	xmm3 += xmm5;		/* Carry += high_FFT_carry * shifted k_hi */ \
	xmm1 *= xptr(rdx+1*32+16);	/* mul by grp two-to-phi */ \
	xmm4 = xptr(rbx+2*32+16);	/* col two-to-phi */ \
	xmm4 *= xptr2(g->u.xmm.XMM_TTP_FUDGE,rax);/* mul by fudge two-to-phi */ \
	xmm1 *= xmm4;		/* data *= fudged col two-to-phi */ \
	xptr(rbp+2*64+16) = xmm1;	/* Save FFT data */ \
\
	rax = u8ptr(rdi+rcx+5); /* Load big vs. little flag */ \
	xmm1 = xptr(rbp+3*64+16);	/* Load FFT data */ \
	xmm1 *= xptr(rdx+1*32);	/* mul by grp two-to-minus-phi */ \
	xmm4 = xptr(rbx+3*32);	/* col two-to-minus-phi */ \
	xmm4 *= xptr2(g->u.xmm.XMM_TTMP_FUDGE,rax);/* mul by fudge two-to-minus-phi */ \
	xmm1 *= xmm4;		/* data *= fudged col two-to-minus-phi */ \
	xmm1 *= g->u.xmm.XMM_NORM012_FF;	/* Mul by FFTLEN/2 */ \
	xmm1 += xmm3;		/* x4 = values1 + carry */ \
	split_carry_zpad_word(base2, xmm6, xmm5, xmm4, rax); \
no##const1(xmm4 = g->u.xmm.XMM_K_LO);		/* Calc high FFT carry times k */ \
const1(xmm4 = g->u.xmm.XMM_K_TIMES_MULCONST_LO); \
	xmm4 *= xmm6;		/* high_FFT_carry * k_lo */ \
	xmm1 += xmm4;		/* x4 = x4 + high_FFT_carry * k_lo */ \
no##const1(xmm4 = g->u.xmm.XMM_K_HI); \
const1(xmm4 = g->u.xmm.XMM_K_TIMES_MULCONST_HI); \
	xmm4 *= xptr2(g->u.xmm.XMM_LIMIT_INVERSE,rax); /* shift k_hi */ \
	xmm6 *= xmm4; \
	rounding(exec, base2, noexec, noexec, xmm1, xmm3, xmm4, rax); \
	xmm3 += xmm6;		/* Carry += high_FFT_carry * shifted k_hi */ \
	xmm1 *= xptr(rdx+1*32+16);	/* mul by grp two-to-phi */ \
	xmm4 = xptr(rbx+3*32+16);	/* col two-to-phi */ \
	xmm4 *= xptr2(g->u.xmm.XMM_TTP_FUDGE,rax);/* mul by fudge two-to-phi */ \
	xmm1 *= xmm4;		/* data *= fudged col two-to-phi */ \
	xptr(rbp+3*64+16) = xmm1;	/* Save FFT data */ \
\
	rcx = g->BIGLIT_INCR4;/* Different clm values step through */ \
;					/* big/lit array differently */ \
	rax = u8ptr(rdi+rcx+1); /* Load big vs. little flag */ \
	xmm1 = xptr(rbp+4*64+16);	/* Load FFT data */ \
	xmm1 *= xptr(rdx+1*32);	/* mul by grp two-to-minus-phi */ \
	xmm4 = xptr(rbx+4*32);	/* col two-to-minus-phi */ \
	xmm4 *= xptr2(g->u.xmm.XMM_TTMP_FUDGE,rax);/* mul by fudge two-to-minus-phi */ \
	xmm1 *= xmm4;		/* data *= fudged col two-to-minus-phi */ \
	xmm1 *= g->u.xmm.XMM_NORM012_FF;	/* Mul by FFTLEN/2 */ \
	xmm1 += xmm3;		/* x5 = values1 + carry */ \
no##const1(xmm4 = g->u.xmm.XMM_K_LO);		/* Calc high FFT carry times k */ \
const1(xmm4 = g->u.xmm.XMM_K_TIMES_MULCONST_LO); \
	xmm4 *= xmm5;		/* high_FFT_carry * k_lo */ \
	xmm1 += xmm4;		/* x5 = x5 + high_FFT_carry * k_lo */ \
no##const1(xmm4 = g->u.xmm.XMM_K_HI); \
const1(xmm4 = g->u.xmm.XMM_K_TIMES_MULCONST_HI); \
	xmm4 *= xptr2(g->u.xmm.XMM_LIMIT_INVERSE,rax); /* shift k_hi */ \
	xmm5 *= xmm4; \
	rounding(exec, base2, noexec, noexec, xmm1, xmm3, xmm4, rax); \
	xmm3 += xmm5;		/* Carry += high_FFT_carry * shifted k_hi */ \
	xmm1 *= xptr(rdx+1*32+16);	/* mul by grp two-to-phi */ \
	xmm4 = xptr(rbx+4*32+16);	/* col two-to-phi */ \
	xmm4 *= xptr2(g->u.xmm.XMM_TTP_FUDGE,rax);/* mul by fudge two-to-phi */ \
	xmm1 *= xmm4;		/* data *= fudged col two-to-phi */ \
	xptr(rbp+4*64+16) = xmm1;	/* Save FFT data */ \
\
	rax = u8ptr(rdi+rcx+5); /* Load big vs. little flag */ \
	xmm3 -= XMM_BIGVAL2;/* Remove rounding const from carry */ \
	xmm3 *= xptr(rdx+1*32+16);	/* mul by grp two-to-phi */ \
	xmm4 = xptr(rbx+5*32+16);	/* col two-to-phi */ \
	xmm4 *= xptr2(g->u.xmm.XMM_TTP_FUDGE,rax);/* mul by fudge two-to-phi */ \
	xmm3 *= xmm4;		/* data *= fudged col two-to-phi */ \
	xmm3 += xptr(rbp+5*64+16);	/* Load FFT data */ \
	xptr(rbp+5*64+16) = xmm3;	/* Save FFT data */ \
\
	xmm4 = XMM_BIGVAL2; \
	xptr(rsi+0*16) = xmm4;	/* Clear carry */ \
	xptr(rsi+1*16) = xmm4; \
	xptr(rsi+2*16) = xmm4; \
	xptr(rsi+3*16) = xmm4; \
\
\
\
/* For WPN macros, these registers are set on input: */ \
/* xmm7 = sumout */ \
/* xmm6 = maxerr */ \
/* rbp = pointer to carries */ \
/* rdi = pointer to big/little flags */ \
/* rsi = pointer to the FFT data */ \
/* rdx = pointer two-to-phi group multipliers */ \
/* ebx = big vs. little & fudge flags */ \
/* rcx = big vs. little word flag #2 */ \
/* rax = big vs. little word flag #1 */ \
\
/* *************** WPN macro ****************** */ \
/* A pipelined version of this code: */ \
/*	movzx	rax, BYTE PTR [rdi]	;; Load big vs. little flags */ \
/*	xload	xmm0, [rsi+0*dist1]	;; Load values1 */ \
/*	mulpd	xmm0, [rdx+0*XMM_GMD][rax] ;; Mul by fudged grp two-to-minus-phi */ \
/*	addpd	xmm0, [rbp+0*16]	;; x1 = values + carry */ \
/*	xload	xmm2, XMM_LIMIT_BIGMAX[rax];; Load maximum * BIGVAL - BIGVAL */ \
/*	addpd	xmm2, xmm0		;; y1 = top bits of x */ \
/*	xload	xmm6, XMM_LIMIT_BIGMAX_NEG[rax];; Load -(maximum*BIGVAL-BIGVAL) */ \
/*	addpd	xmm6, xmm2		;; z1 = y1-(maximum * BIGVAL - BIGVAL) */ \
/*	subpd	xmm0, xmm6		;; rounded value = x1 - z1 */ \
/*	mulpd	xmm2, XMM_LIMIT_INVERSE[rax];; next carry = shifted y1 */ \
/*	mulpd	xmm0, [rdx+0*XMM_GMD+XMM_GMD/2][rax];; new value1 = val * fudged grp two-to-phi */ \
/*	xstore	[rsi+0*dist1], xmm0	;; Save new value1 */ \
/*	xstore	[rbp+0*16], xmm2	;; Save carry */ \
/* */ \
\

#if 0
#define xnorm_wpn_preload(ttp, echk, const1, base2, sse4)


#define xnorm_wpn(ttp, echk, const1, base2, sse4, CARRY0, CARRY1, CARRY2, CARRY3) { \
ttp(	uintptr_t rcx = rbx & 0xFF);				/* Fudge flags 1-2 */ \
ttp(	rcx &= 0xF0); \
		xmm0 = xptr(rsi+0*16);		/* Load values1 */ \
ttp(	xmm0 *= xptr(rdx+0*XMM_GMD+rcx));	/* Mul by fudged grp two-to-minus-phi */ \
		xmm1 = xptr(rsi+1*16);		/* Load values2 */ \
ttp(	xmm1 *= xptr(rdx+1*XMM_GMD+rcx));	/* Mul by fudged grp two-to-minus-phi */ \
ttp(	uintptr_t rax = (rbx >> 8) & 0xFF);				/* Big/lit flags 1-4 */ \
no##ttp(base2(uintptr_t rax = 0));			/* --We should clean up base 2 rational FFT so that this isn't needed */ \
no##const1(echk(error_check_interleaved(sse4, xmm0, xmm4, xmm1, xmm5, xmm6))); \
const1(	mul_by_const_interleaved(ttp, base2, echk, sse4, xmm0, xmm4, xmm2, rax*4, xmm1, xmm5, xmm3, rax*4+16, xmm6)); \
		xmm0 += CARRY0;		/* x1 = values + carry */ \
		xmm1 += CARRY1;		/* x2 = values + carry */ \
\
no##base2(rounding_interleaved(ttp, base2, const1, sse4, xmm0, xmm4, xmm2, rax*4, xmm1, xmm5, xmm3, rax*4+16)); \
\
base2(	xmm2 = xptr2(g->u.xmm.XMM_LIMIT_BIGMAX,rax*4));	/* Load maximum * BIGVAL - BIGVAL */ \
base2(no##const1(xmm5 = xmm2));			/* Copy maximum * BIGVAL - BIGVAL */ \
base2(const1(no##echk(xmm6 = xmm2)));			/* Copy maximum * BIGVAL - BIGVAL */ \
base2(	xmm2 += xmm0);			/* y1 = top bits of x */ \
base2(	xmm3 = xptr2(g->u.xmm.XMM_LIMIT_BIGMAX,rax*4+16));/* Load maximum * BIGVAL - BIGVAL */ \
base2(	xmm7 = xmm3);			/* Copy maximum * BIGVAL - BIGVAL */ \
base2(	xmm3 += xmm1);			/* y2 = top bits of x */ \
base2(const1(xmm4 += xmm2));			/* Add in upper mul-by-const bits */ \
base2(const1(xmm4 *= xptr2(g->u.xmm.XMM_LIMIT_INVERSE,rax*4)));	/* next carry = shifted y1 */ \
base2(no##const1(xmm4 = xptr2(g->u.xmm.XMM_LIMIT_INVERSE,rax*4)));	/* next carry = shifted y1 */ \
base2(no##const1(xmm4 *= xmm2));			/* next carry = shifted y1 */ \
		CARRY0 = xmm4;		/* Save carry1 */ \
base2(const1(echk(xmm2 -= xptr2(g->u.xmm.XMM_LIMIT_BIGMAX,rax*4))));	/* z1 = y1 - (maximum*BIGVAL-BIGVAL) */ \
base2(const1(no##echk(xmm2 -= xmm6)));			/* z1 = y1 - (maximum*BIGVAL-BIGVAL) */ \
base2(no##const1(xmm2 -= xmm5));			/* z1 = y1 - (maximum*BIGVAL-BIGVAL) */ \
base2(const1(xmm5 += xmm3));			/* Add in upper mul-by-const bits */ \
base2(const1(xmm5 *= xptr2(g->u.xmm.XMM_LIMIT_INVERSE,rax*4+16)));/* next carry = shifted y2 */ \
base2(no##const1(xmm5 = xptr2(g->u.xmm.XMM_LIMIT_INVERSE,rax*4+16)));/* next carry = shifted y2 */ \
base2(no##const1(xmm5 *= xmm3));			/* next carry = shifted y2 */ \
		CARRY1 = xmm5;		/* Save carry2 */ \
base2(	xmm3 -= xmm7);			/* z2 = y2 - (maximum*BIGVAL-BIGVAL) */ \
base2(	xmm0 -= xmm2);			/* rounded value = x1 - z1 */ \
base2(	xmm1 -= xmm3);			/* rounded value = x2 - z2 */ \
ttp(	xmm0 *= xptr(rdx+0*XMM_GMD+XMM_GMD/2+rcx)); /* value1 = rounded value * fudged grp two-to-phi */ \
		xptr(rsi+0*16) = xmm0;		/* Save new value1 */ \
ttp(	xmm1 *= xptr(rdx+1*XMM_GMD+XMM_GMD/2+rcx)); /* value2 = rounded value * fudged grp two-to-phi */ \
		xptr(rsi+1*16) = xmm1;		/* Save new value2 */ \
\
ttp(	rbx &= 0x0F);/* Fudge flags 3-4 */ \
		xmm0 = xptr(rsi+2*16);		/* Load values1 */ \
ttp(	xmm0 *= xptr(rdx+2*XMM_GMD+rbx*8));	/* Mul by fudged grp two-to-minus-phi */ \
		xmm1 = xptr(rsi+3*16);		/* Load values2 */ \
ttp(	xmm1 *= xptr(rdx+3*XMM_GMD+rbx*8));	/* Mul by fudged grp two-to-minus-phi */ \
no##const1(echk(error_check_interleaved(sse4, xmm0, xmm4, xmm1, xmm5, xmm6))); \
const1(	mul_by_const_interleaved(ttp, base2, echk, sse4, xmm0, xmm4, xmm2, rax*4+32, xmm1, xmm5, xmm3, rax*4+48, xmm6)); \
		xmm0 += CARRY2;		/* x1 = values + carry */ \
		xmm1 += CARRY3;		/* x2 = values + carry */ \
\
no##base2(rounding_interleaved(ttp, base2, const1, sse4, xmm0, xmm4, xmm2, rax*4+32, xmm1, xmm5, xmm3, rax*4+48)); \
\
base2(	xmm2 = xptr2(g->u.xmm.XMM_LIMIT_BIGMAX,rax*4+32));/* Load maximum * BIGVAL - BIGVAL */ \
base2(no##const1(xmm5 = xmm2));			/* Copy maximum * BIGVAL - BIGVAL */ \
base2(const1(no##echk(xmm6 = xmm2)));			/* Copy maximum * BIGVAL - BIGVAL */ \
base2(	xmm2 += xmm0);			/* y1 = top bits of x */ \
base2(	xmm3 = xptr2(g->u.xmm.XMM_LIMIT_BIGMAX,rax*4+48));/* Load maximum * BIGVAL - BIGVAL */ \
base2(	xmm7 = xmm3);			/* Copy maximum * BIGVAL - BIGVAL */ \
base2(	xmm3 += xmm1);			/* y2 = top bits of x */ \
base2(const1(xmm4 += xmm2));			/* Add in upper mul-by-const bits */ \
base2(const1(xmm4 *= xptr2(g->u.xmm.XMM_LIMIT_INVERSE,rax*4+32)));/* next carry = shifted y1 */ \
base2(no##const1(xmm4 = xptr2(g->u.xmm.XMM_LIMIT_INVERSE,rax*4+32)));/* next carry = shifted y1 */ \
base2(no##const1(xmm4 *= xmm2));			/* next carry = shifted y1 */ \
		CARRY2 = xmm4;		/* Save carry1 */ \
base2(const1(echk(xmm2 -= xptr2(g->u.xmm.XMM_LIMIT_BIGMAX,rax*4+32))));/* z1 = y1 - (maximum*BIGVAL-BIGVAL) */ \
base2(const1(no##echk(xmm2 -= xmm6)));			/* z1 = y1 - (maximum*BIGVAL-BIGVAL) */ \
base2(no##const1(xmm2 -= xmm5));			/* z1 = y1 - (maximum*BIGVAL-BIGVAL) */ \
base2(const1(xmm5 += xmm3));			/* Add in upper mul-by-const bits */ \
base2(const1(xmm5 *= xptr2(g->u.xmm.XMM_LIMIT_INVERSE,rax*4+48)));/* next carry = shifted y2 */ \
base2(no##const1(xmm5 = xptr2(g->u.xmm.XMM_LIMIT_INVERSE,rax*4+48)));/* next carry = shifted y2 */ \
base2(no##const1(xmm5 *= xmm3));			/* next carry = shifted y2 */ \
		CARRY3 = xmm5;		/* Save carry2 */ \
base2(	xmm3 -= xmm7);			/* z2 = y2 - (maximum*BIGVAL-BIGVAL) */ \
base2(	xmm0 -= xmm2);			/* rounded value = x1 - z1 */ \
base2(	xmm1 -= xmm3);			/* rounded value = x2 - z2 */ \
ttp(	xmm0 *= xptr(rdx+2*XMM_GMD+XMM_GMD/2+rbx*8)); /* value1 *= rounded value * fudged grp two-to-phi */ \
ttp(	xmm1 *= xptr(rdx+3*XMM_GMD+XMM_GMD/2+rbx*8)); /* value2 = rounded value * fudged grp two-to-phi */ \
ttp(	rbx = *(unsigned short*)(rdi+2));		/* Load next 4 big vs. little & fudge flags */ \
		xptr(rsi+2*16) = xmm0;		/* Save new value1 */ \
		xptr(rsi+3*16) = xmm1;		/* Save new value2 */ \
}

#else

#define error_check_interleaved_64 error_check_interleaved

#define xnorm_wpn_preload(ttp, echk, const1, base2, sse4) \
/*echk(sse4(xmm15 = XMM_ABSVAL));*/ \
/*echk(no##sse4(xmm15 = XMM_BIGVAL1));*/ \

#define xnorm_wpn(ttp, echk, const1, base2, sse4, CARRY0, CARRY1, CARRY2, CARRY3) \
	no##const1(xmm12  = CARRY0);	/* Preload carries */ \
	no##const1(xmm13 = CARRY1); \
	no##const1(xmm4  = CARRY2); \
	no##const1(xmm5 = CARRY3); \
ttp(	uintptr_t rcx = rbx & 0xFF);				/* Fudge flags 1,2 */ \
ttp(	rcx &= 0xF0); \
		xmm8 = xptr(rsi+0*16);		/* Load values1				;P4	;Core2 */ \
ttp(	xmm8 *= xptr(rdx+0*XMM_GMD+rcx));	/* Mul by fudged grp two-to-minus-phi	;1-5	;1-4 */ \
		xmm9 = xptr(rsi+1*16);		/* Load values2 */ \
ttp(	xmm9 *= xptr(rdx+1*XMM_GMD+rcx));	/* Mul by fudged grp two-to-minus-phi	;3-7	;2-5 */ \
ttp(	uintptr_t rax = (rbx >> 8) & 0xFF);				/* Big/lit flags 1-4 */ \
no##ttp(base2(uintptr_t rax = 0));			/* --We should clean up base 2 rational FFT so that this isn't needed */ \
no##const1(echk(error_check_interleaved_64(sse4, xmm8, xmm10, xmm9, xmm11, xmm6))); \
const1(	mul_by_const_interleaved(ttp, base2, echk, sse4, xmm8, xmm12, xmm10, rax*4, xmm9, xmm13, xmm11, rax*4+16, xmm6)); \
ttp(	rbx &= 0x0F);/* Fudge flags 3,4 */ \
		xmm0 = xptr(rsi+2*16);		/* Load values3 */ \
ttp(	xmm0 *= xptr(rdx+2*XMM_GMD+rbx*8));	/* Mul by fudged grp two-to-minus-phi	;5-9	;3-6 */ \
		xmm1 = xptr(rsi+3*16);		/* Load values4 */ \
ttp(	xmm1 *= xptr(rdx+3*XMM_GMD+rbx*8));	/* Mul by fudged grp two-to-minus-phi	;7-11	;4-7 */ \
no##const1(echk(error_check_interleaved_64(sse4, xmm0, xmm2, xmm1, xmm3, xmm6))); \
const1(	mul_by_const_interleaved(ttp, base2, echk, sse4, xmm0, xmm4, xmm2, rax*4+32, xmm1, xmm5, xmm3, rax*4+48, xmm6)); \
const1(	xmm8 += CARRY0);		/* x1 = values + carry */ \
no##const1(xmm8 += xmm12);			/* x1 = values + carry			;6-9	;5-7 */ \
const1(	xmm9 += CARRY1);		/* x2 = values + carry */ \
no##const1(xmm9 += xmm13);			/* x2 = values + carry			;8-11	;6-8 */ \
\
no##base2(rounding_interleaved(ttp, base2, const1, sse4, xmm8, xmm12, xmm10, rax*4, xmm9, xmm13, xmm11, rax*4+16)); \
\
const1(	xmm0 += CARRY2);		/* x3 = values + carry */ \
no##const1(xmm0 += xmm4);			/* x3 = values + carry			;10-13	;7-9 */ \
base2(	xmm10 = xptr2(g->u.xmm.XMM_LIMIT_BIGMAX,rax*4));	/* Load maximum * BIGVAL - BIGVAL */ \
base2(no##const1(xmm4 = xmm10));			/* Copy maximum * BIGVAL - BIGVAL */ \
base2(const1(no##echk(xmm15 = xmm10)));			/* Copy maximum * BIGVAL - BIGVAL */ \
base2(	xmm10 += xmm8);			/* y1 = top bits of x			;12-15	;8-10 */ \
const1(	xmm1 += CARRY3);		/* x4 = values + carry */ \
no##const1(xmm1 += xmm5);			/* x4 = values + carry			;14-17	;9-11 */ \
\
no##base2(rounding_interleaved(ttp, base2, const1, sse4, xmm0, xmm4, xmm2, rax*4+32, xmm1, xmm5, xmm3, rax*4+48)); \
\
base2(	xmm11 = xptr2(g->u.xmm.XMM_LIMIT_BIGMAX,rax*4+16)); /* Load maximum * BIGVAL - BIGVAL */ \
base2(no##const1(xmm5 = xmm11));			/* Copy maximum * BIGVAL - BIGVAL */ \
base2(const1(no##echk(xmm6 = xmm11)));			/* Copy maximum * BIGVAL - BIGVAL */ \
base2(	xmm11 += xmm9);			/* y2 = top bits of x			;16-19	;10-12 */ \
base2(const1(xmm12 += xmm10));			/* Add in upper mul-by-const bits */ \
base2(const1(xmm12 *= xptr2(g->u.xmm.XMM_LIMIT_INVERSE,rax*4)));	/* next carry = shifted y1 */ \
base2(no##const1(xmm12 = xptr2(g->u.xmm.XMM_LIMIT_INVERSE,rax*4)));	/* next carry = shifted y1 */ \
base2(no##const1(xmm12 *= xmm10));			/* next carry = shifted y1		;17-22	;11-15 */ \
base2(	xmm2 = xptr2(g->u.xmm.XMM_LIMIT_BIGMAX,rax*4+32));/* Load maximum * BIGVAL - BIGVAL */ \
base2(	xmm7 = xmm2);			/* Copy maximum * BIGVAL - BIGVAL */ \
base2(	xmm2 += xmm0);			/* y3 = top bits of x			;18-21	;11-13 */ \
base2(	xmm3 = xptr2(g->u.xmm.XMM_LIMIT_BIGMAX,rax*4+48));/* Load maximum * BIGVAL - BIGVAL */ \
base2(	xmm14 = xmm3);			/* Copy maximum * BIGVAL - BIGVAL */ \
base2(	xmm3 += xmm1);			/* y4 = top bits of x			;20-23	;12-14 */ \
base2(const1(xmm13 += xmm11));			/* Add in upper mul-by-const bits */ \
base2(const1(xmm13 *= xptr2(g->u.xmm.XMM_LIMIT_INVERSE,rax*4+16)));/* next carry = shifted y2 */ \
base2(no##const1(xmm13 = xptr2(g->u.xmm.XMM_LIMIT_INVERSE,rax*4+16)));/* next carry = shifted y2 */ \
base2(no##const1(xmm13 *= xmm11));			/* next carry = shifted y2		;21-26	;13-17 */ \
base2(const1(echk(xmm10 -= xptr2(g->u.xmm.XMM_LIMIT_BIGMAX,rax*4))));	/* z1 = y1 - (maximum*BIGVAL-BIGVAL) */ \
base2(const1(no##echk(xmm10 -= xmm15)));			/* z1 = y1 - (maximum*BIGVAL-BIGVAL) */ \
base2(no##const1(xmm10 -= xmm4));			/* z1 = y1 - (maximum*BIGVAL-BIGVAL)	;22-25	;13-15 */ \
base2(const1(xmm4 += xmm2));			/* Add in upper mul-by-const bits */ \
base2(const1(xmm4 *= xptr2(g->u.xmm.XMM_LIMIT_INVERSE,rax*4+32)));/* next carry = shifted y3 */ \
base2(no##const1(xmm4 = xptr2(g->u.xmm.XMM_LIMIT_INVERSE,rax*4+32)));/* next carry = shifted y3 */ \
base2(no##const1(xmm4 *= xmm2));			/* next carry = shifted y3		;23-28	;14-18 */ \
base2(const1(echk(xmm11 -= xptr2(g->u.xmm.XMM_LIMIT_BIGMAX,rax*4+16))));/* z2 = y2 - (maximum*BIGVAL-BIGVAL) */ \
base2(const1(no##echk(xmm11 -= xmm6)));			/* z2 = y2 - (maximum*BIGVAL-BIGVAL) */ \
base2(no##const1(xmm11 -= xmm5));			/* z2 = y2 - (maximum*BIGVAL-BIGVAL)	;24-27	;14-16 */ \
base2(const1(xmm5 += xmm3));			/* Add in upper mul-by-const bits */ \
base2(const1(xmm5 *= xptr2(g->u.xmm.XMM_LIMIT_INVERSE,rax*4+48)));/* next carry = shifted y4 */ \
base2(no##const1(xmm5 = xptr2(g->u.xmm.XMM_LIMIT_INVERSE,rax*4+48)));/* next carry = shifted y4 */ \
base2(no##const1(xmm5 *= xmm3));			/* next carry = shifted y4		;25-30	;15-19 */ \
const1(	CARRY0 = xmm12);		/* Save carry1 */ \
base2(	xmm2 -= xmm7);			/* z3 = y3 - (maximum*BIGVAL-BIGVAL)	;26-29	;15-17 */ \
const1(	CARRY1 = xmm13);		/* Save carry2 */ \
base2(	xmm3 -= xmm14);			/* z4 = y4 - (maximum*BIGVAL-BIGVAL)	;28-31	;16-18 */ \
base2(	xmm8 -= xmm10);			/* rounded value = x1 - z1		;30-33	;17-19 */ \
const1(	CARRY2 = xmm4);		/* Save carry3 */ \
base2(	xmm9 -= xmm11);			/* rounded value = x2 - z2		;32-35	;18-20 */ \
const1(	CARRY3 = xmm5);		/* Save carry4 */ \
base2(	xmm0 -= xmm2);			/* rounded value = x3 - z3		;34-37	;19-21 */ \
ttp(	xmm8 *= xptr(rdx+0*XMM_GMD+XMM_GMD/2+rcx)); /* value1 *= fudged grp two-to-phi	;35-40	;20-24 */ \
base2(	xmm1 -= xmm3);			/* rounded value = x4 - z4		;36-39	;20-22 */ \
ttp(	xmm9 *= xptr(rdx+1*XMM_GMD+XMM_GMD/2+rcx)); /* value2 *= fudged grp two-to-phi	;37-42	;21-25 */ \
ttp(	xmm0 *= xptr(rdx+2*XMM_GMD+XMM_GMD/2+rbx*8)); /* value3 *= fudged grp two-to-phi;38-41	;22-26 */ \
ttp(	xmm1 *= xptr(rdx+3*XMM_GMD+XMM_GMD/2+rbx*8)); /* value4 *= fudged grp two-to-phi;40-43	;23-27 */ \
ttp(	rbx = *(unsigned short*)(rdi+2));		/* Load next 4 big vs. little & fudge flags */ \
		xptr(rsi+0*16) = xmm8;		/* Save new value1 */ \
		xptr(rsi+1*16) = xmm9;		/* Save new value2 */ \
		xptr(rsi+2*16) = xmm0;		/* Save new value3 */ \
		xptr(rsi+3*16) = xmm1;		/* Save new value4 */ \
	no##const1(CARRY0 = xmm12); /* Store carries */ \
	no##const1(CARRY1 = xmm13); \
	no##const1(CARRY2 = xmm4); \
	no##const1(CARRY3 = xmm5); \

#endif

/* NOTE: We'd rather store the high FFT carry without the XMM_BIGVAL added in, */ \
/* but there is too much code that expects this (like add and subtract). */ \
/* ALSO NOTE:  In the zero pad case, big/lit & fudge flags 1,2 and the */ \
/* same as big/lit & fudge flags 3,4. */ \

#if 0
#define xnorm_wpn_zpad_preload(ttp, echk, const1, base2, sse4, khi, c1, cm1) \

#define xnorm_wpn_zpad(ttp, echk, const1, base2, sse4, khi, c1, cm1, CARRY0, CARRY1, CARRY2, CARRY3) { \
no##ttp(uintptr_t rax = 0); \
ttp(uintptr_t rax = (rbx >> 8) & 0xFF);			/* Big/little flags 1-4 */ \
ttp(rbx &= 0xF0);/* Fudge flags 1,2 */ \
ttp(xmm2 = xptr(rdx+0*XMM_GMD+rbx)); /* Fudged grp two-to-minus-phi */ \
	xmm0 = xptr(rsi);		/* Load values1 */ \
	xmm1 = xptr(rsi+2*16);	/* Load values2 */ \
ttp(xmm0 *= xmm2);		/* Mul by fudged grp two-to-minus-phi */ \
ttp(xmm1 *= xmm2);		/* Mul by fudged grp two-to-minus-phi */ \
\
	xmm3 = CARRY2;	/* Add in previous high FFT data */ \
	split_lower_zpad_word(echk, base2, sse4, xmm0, xmm3, xmm4, rax*4); \
	CARRY2 = xmm3; \
\
no##const1(xmm0 = g->u.xmm.XMM_K_LO); \
const1(	xmm0 = g->u.xmm.XMM_K_TIMES_MULCONST_LO); \
		xmm0 *= xmm4; \
khi(no##const1(xmm5 = g->u.xmm.XMM_K_HI)); \
khi(const1(xmm5 = g->u.xmm.XMM_K_TIMES_MULCONST_HI)); \
khi(no##base2(xmm5 *= xptr2(g->u.xmm.XMM_LIMIT_INVERSE,rax*4))); /* Non-base2 rounding needs shifted carry */ \
khi(	xmm5 *= xmm4); \
;  \
		xmm0 += CARRY0;	/* x1 = values + carry */ \
\
c1(	xmm1 *= g->u.xmm.XMM_MINUS_C);	/* Do one mul before split rather than two after split */ \
\
	split_upper_zpad_word(echk, base2, sse4, xmm1, xmm4, xmm2, rax*4); \
\
no##const1(no##c1(no##cm1(xmm2 *= g->u.xmm.XMM_MINUS_C))); \
no##const1(no##c1(no##cm1(xmm4 *= g->u.xmm.XMM_MINUS_C))); \
const1(		xmm2 *= g->u.xmm.XMM_MINUS_C_TIMES_MULCONST); \
const1(		xmm4 *= g->u.xmm.XMM_MINUS_C_TIMES_MULCONST); \
\
	xmm0 += xmm2;		/* Add upper FFT word to lower FFT word */ \
khi(xmm4 += xmm5);		/* Add upper FFT word to lower FFT word */ \
\
	rounding(ttp, base2, exec, sse4, xmm0, xmm4, xmm2, rax*4); \
\
ttp(xmm0 *= xptr(rdx+0*XMM_GMD+XMM_GMD/2+rbx)); /* new value1 *= fudged grp two-to-phi */ \
	CARRY0 = xmm4;	/* Save carry */ \
	xptr(rsi) = xmm0;		/* Save new value1 */ \
\
ttp(xmm3 = xptr(rdx+1*XMM_GMD+rbx)); /* Fudged grp two-to-minus-phi */ \
	xmm5 = xptr(rsi+16);		/* Load high values1 */ \
	xmm1 = xptr(rsi+3*16);	/* Load high values2 */ \
ttp(xmm5 *= xmm3);		/* Mul by fudged grp two-to-minus-phi */ \
ttp(xmm1 *= xmm3);		/* Mul by fudged grp two-to-minus-phi */ \
\
	xmm3 = CARRY3;	/* Add in previous high FFT data */ \
	split_lower_zpad_word(echk, base2, sse4, xmm5, xmm3, xmm2, rax*4+16); \
	CARRY3 = xmm3; \
\
no##const1(xmm0 = g->u.xmm.XMM_K_LO); \
const1(	xmm0 = g->u.xmm.XMM_K_TIMES_MULCONST_LO); \
		xmm0 *= xmm2; \
khi(no##const1(xmm5 = g->u.xmm.XMM_K_HI)); \
khi(const1(xmm5 = g->u.xmm.XMM_K_TIMES_MULCONST_HI)); \
khi(no##base2(xmm5 *= xptr2(g->u.xmm.XMM_LIMIT_INVERSE,rax*4+16))); /* Non-base2 rounding needs shifted carry */ \
khi(	xmm5 *= xmm2); \
\
		xmm0 += CARRY1;	/* x2 = values + carry */ \
\
c1(	xmm1 *= g->u.xmm.XMM_MINUS_C);	/* Do one mul before split rather than two after split */ \
\
	split_upper_zpad_word(echk, base2, sse4, xmm1, xmm2, xmm4, rax*4+16); \
\
no##const1(no##c1(no##cm1(xmm4 *= g->u.xmm.XMM_MINUS_C))); \
no##const1(no##c1(no##cm1(xmm2 *= g->u.xmm.XMM_MINUS_C))); \
const1(		xmm4 *= g->u.xmm.XMM_MINUS_C_TIMES_MULCONST); \
const1(		xmm2 *= g->u.xmm.XMM_MINUS_C_TIMES_MULCONST); \
\
	xmm0 += xmm4;		/* Add upper FFT word to lower FFT word */ \
khi(xmm2 += xmm5);		/* Add upper FFT word to lower FFT word */ \
\
	rounding(ttp, base2, exec, sse4, xmm0, xmm2, xmm4, rax*4+16); \
\
ttp(xmm0 *= xptr(rdx+1*XMM_GMD+XMM_GMD/2+rbx)); /* new value2 *= fudged grp two-to-phi */ \
ttp(rbx = *(unsigned short*)(rdi+2));	/* Load next big vs. little & fudge flags */ \
	CARRY1 = xmm2;	/* Save carry */ \
	xptr(rsi+1*16) = xmm0;	/* Save new value2 */ \
\
	xmm1[1] = xmm1[0] = 0.0;		/* new high values = zero */ \
	xptr(rsi+2*16) = xmm1;	/* Zero high value1 */ \
	xptr(rsi+3*16) = xmm1;	/* Zero high value2 */ \
}

#else
/* 64-bit version using extra registers */

#define xnorm_wpn_zpad_preload(ttp, echk, const1, base2, sse4, khi, c1, cm1) \
no##const1(		xmm15 = g->u.xmm.XMM_K_LO); \
const1(			xmm15 = g->u.xmm.XMM_K_TIMES_MULCONST_LO); \
no##const1(no##c1(no##cm1(xmm14 = g->u.xmm.XMM_MINUS_C))); \
no##const1(    c1(        xmm14 = g->u.xmm.XMM_MINUS_C)); \
const1(			xmm14 = g->u.xmm.XMM_MINUS_C_TIMES_MULCONST); \
/*base2(			xmm13 = g->u.xmm.XMM_BIGBIGVAL);*/ \
/*base2(no##echk(no##sse4(	xmm6 = XMM_BIGVAL2)));*/ \
base2(no##echk(sse4(khi(no##const1(xmm6 = g->u.xmm.XMM_K_HI))))); \
base2(no##echk(sse4(khi(const1(xmm6 = g->u.xmm.XMM_K_TIMES_MULCONST_HI))))); \
/*no##base2(no##sse4(	xmm13 = XMM_BIGVAL2));*/ \
/*no##base2(sse4(echk(	xmm13 = XMM_ABSVAL)));*/ \
no##base2(no##echk(khi(no##const1(xmm6 = g->u.xmm.XMM_K_HI)))); \
no##base2(no##echk(khi(const1(xmm6 = g->u.xmm.XMM_K_TIMES_MULCONST_HI)))); \



#define xnorm_wpn_zpad(ttp, echk, const1, base2, sse4, khi, c1, cm1, CARRY0, CARRY1, CARRY2, CARRY3) \
	xmm4  = CARRY0;	/* Preload carries */ \
	xmm11 = CARRY1; \
	xmm3  = CARRY2; \
	xmm10 = CARRY3; \
no##ttp(const uintptr_t rax = 0); \
\
ttp(uintptr_t rax = (rbx >> 8) & 0xFF);			/* Big/little flags 1-4 */ \
ttp(rbx &= 0xF0);/* Fudge flags 1,2 */ \
\
ttp(xmm2 = xptr(rdx+0*XMM_GMD+rbx)); /* Fudged grp two-to-minus-phi */ \
	xmm0 = xptr(rsi);		/* Load values1 */ \
	xmm1 = xptr(rsi+2*16);	/* Load values2 */ \
ttp(xmm0 *= xmm2);		/* Mul by fudged grp two-to-minus-phi */ \
ttp(xmm1 *= xmm2);		/* Mul by fudged grp two-to-minus-phi */ \
\
ttp(xmm9 = xptr(rdx+1*XMM_GMD+rbx)); /* Fudged grp two-to-minus-phi */ \
	xmm7 = xptr(rsi+16);		/* Load high values1 */ \
	xmm8 = xptr(rsi+3*16);	/* Load high values2 */ \
ttp(xmm7 *= xmm9);		/* Mul by fudged grp two-to-minus-phi */ \
ttp(xmm8 *= xmm9);		/* Mul by fudged grp two-to-minus-phi */ \
\
;					/* Split lower word adding in previous high FFT data */ \
	split_lower_zpad_word_interleaved(echk, base2, sse4, xmm0, xmm3, xmm5, rax*4, xmm7, xmm10, xmm12, rax*4+16); \
\
no##const1(xmm0 = xmm15);	/* XMM_K_LO */ \
const1(	xmm0 = xmm15);	/* XMM_K_TIMES_MULCONST_LO */ \
		xmm0 *= xmm5; \
no##const1(xmm7 = xmm15);	/* XMM_K_LO */ \
const1(	xmm7 = xmm15);	/* XMM_K_TIMES_MULCONST_LO */ \
		xmm7 *= xmm12; \
\
		xmm0 += xmm4;	/* x1 = values + carry */ \
		xmm7 += xmm11;	/* x2 = values + carry */ \
\
khi(no##const1(no##base2(echk(	xmm4 = g->u.xmm.XMM_K_HI)))); \
khi(no##const1(no##base2(no##echk(	xmm4 = xmm6))));	/* XMM_K_HI */ \
khi(no##const1(base2(echk(		xmm4 = g->u.xmm.XMM_K_HI)))); \
khi(no##const1(base2(no##echk(no##sse4(xmm4 = g->u.xmm.XMM_K_HI))))); \
khi(no##const1(base2(no##echk(sse4(	xmm4 = xmm6)))));	/* XMM_K_HI */ \
khi(const1(no##base2(echk(		xmm4 = g->u.xmm.XMM_K_TIMES_MULCONST_HI)))); \
khi(const1(no##base2(no##echk(	xmm4 = xmm6))));	/* XMM_K_TIMES_MULCONST_HI */ \
khi(const1(base2(echk(		xmm4 = g->u.xmm.XMM_K_TIMES_MULCONST_HI)))); \
khi(const1(base2(no##echk(no##sse4(	xmm4 = g->u.xmm.XMM_K_TIMES_MULCONST_HI))))); \
khi(const1(base2(no##echk(sse4(	xmm4 = xmm6)))));	/* XMM_K_TIMES_MULCONST_HI */ \
khi(no##base2(			xmm4 *= xptr2(g->u.xmm.XMM_LIMIT_INVERSE,rax*4))); /* Non-base2 rounding needs shifted carry */ \
khi(				xmm5 *= xmm4); \
\
khi(no##const1(no##base2(echk(	xmm11 = g->u.xmm.XMM_K_HI)))); \
khi(no##const1(no##base2(no##echk(	xmm11 = xmm6))));	/* XMM_K_HI */ \
khi(no##const1(base2(echk(		xmm11 = g->u.xmm.XMM_K_HI)))); \
khi(no##const1(base2(no##echk(no##sse4(xmm11 = g->u.xmm.XMM_K_HI))))); \
khi(no##const1(base2(no##echk(sse4(	xmm11 = xmm6)))));	/* XMM_K_HI */ \
khi(const1(no##base2(echk(		xmm11 = g->u.xmm.XMM_K_TIMES_MULCONST_HI)))); \
khi(const1(no##base2(no##echk(	xmm11 = xmm6))));	/* XMM_K_TIMES_MULCONST_HI */ \
khi(const1(base2(echk(		xmm11 = g->u.xmm.XMM_K_TIMES_MULCONST_HI)))); \
khi(const1(base2(no##echk(no##sse4(	xmm11 = g->u.xmm.XMM_K_TIMES_MULCONST_HI))))); \
khi(const1(base2(no##echk(sse4(	xmm11 = xmm6)))));	/* XMM_K_TIMES_MULCONST_HI */ \
khi(no##base2(			xmm11 *= xptr2(g->u.xmm.XMM_LIMIT_INVERSE,rax*4+16))); /* Non-base2 rounding needs shifted carry */ \
khi(				xmm12 *= xmm11); \
\
c1(xmm1 *= xmm14);		/* XMM_MINUS_C	;; Do one mul before split rather than two after split */ \
c1(xmm8 *= xmm14);		/* XMM_MINUS_C	;; Do one mul before split rather than two after split */ \
\
	split_upper_zpad_word_interleaved(echk, base2, sse4, xmm1, xmm4, xmm2, rax*4, xmm8, xmm11, xmm9, rax*4+16); \
\
no##const1(no##c1(no##cm1(xmm2 *= xmm14)));	/* XMM_MINUS_C */ \
no##const1(no##c1(no##cm1(xmm4 *= xmm14)));	/* XMM_MINUS_C */ \
const1(		xmm2 *= xmm14);	/* XMM_MINUS_C_TIMES_MULCONST */ \
const1(		xmm4 *= xmm14);	/* XMM_MINUS_C_TIMES_MULCONST */ \
no##const1(no##c1(no##cm1(xmm9 *= xmm14)));	/* XMM_MINUS_C */ \
no##const1(no##c1(no##cm1(xmm11 *= xmm14)));	/* XMM_MINUS_C */ \
const1(		xmm9 *= xmm14);	/* XMM_MINUS_C_TIMES_MULCONST */ \
const1(		xmm11 *= xmm14);	/* XMM_MINUS_C_TIMES_MULCONST */ \
\
	xmm0 += xmm2;		/* Add upper FFT word to lower FFT word */ \
khi(xmm4 += xmm5);		/* Add upper FFT word to lower FFT word */ \
	xmm7 += xmm9;		/* Add upper FFT word to lower FFT word */ \
khi(xmm11 += xmm12);		/* Add upper FFT word to lower FFT word */ \
\
	rounding_interleaved(ttp, base2, exec, sse4, xmm0, xmm4, xmm2, rax*4, xmm7, xmm11, xmm9, rax*4+16); \
\
ttp(xmm0 *= xptr(rdx+0*XMM_GMD+XMM_GMD/2+rbx)); /* new value1 *= fudged grp two-to-phi */ \
ttp(xmm7 *= xptr(rdx+1*XMM_GMD+XMM_GMD/2+rbx)); /* new value2 *= fudged grp two-to-phi */ \
\
ttp(rbx = *(unsigned short*)(rdi+2));	/* Load next big vs. little & fudge flags */ \
\
	xptr(rsi) = xmm0;		/* Save new value1 */ \
	xptr(rsi+1*16) = xmm7;	/* Save new value2 */ \
\
	xmm1[1] = xmm1[0] = 0.0;		/* new high values = zero */ \
	xptr(rsi+2*16) = xmm1;	/* Zero high value1 */ \
	xptr(rsi+3*16) = xmm1;	/* Zero high value2 */ \
	CARRY0 = xmm4; /* Store carries */ \
	CARRY1 = xmm11; \
	CARRY2 = xmm3; \
	CARRY3 = xmm10; \

#endif

\
/* *************** WPN followup macro ****************** */ \
/* This macro finishes the normalize process by adding the final carries */ \
/* back into the appropriate FFT values. */ \
/* rsi = pointer to carries */ \
/* rbp = pointer to FFT data */ \
/* rdi = pointer to big/little flags */ \
/* rdx = pointer two-to-phi group multipliers */ \
/* rcx = big vs. little word flag #2 */ \
/* rax = big vs. little word flag #1 */ \
\
/* Rotate the carries array */ \

#define xnorm012_wpn_part1 \
if(g->this_block == 0) {		/* Are we carrying into the first data block? */ \
	/*jne done;*/		/* If not, skip wrapping carries and negating the last carry */ \
	int loopZ = g->addcount1;	/* Load count of cache lines in the carries array */ \
do{ \
	xmm4[0] = f64ptr(rsi+0*16);	/* Load low carry word */ \
	f64ptr(rsi+0*16) = xmm0[0];	/* Save prev cache line's high carry in low word */ \
	xmm0[0] = f64ptr(rsi+0*16+8);	/* Load high carry for next cache line */ \
	f64ptr(rsi+0*16+8) = xmm4[0];	/* Move this cache line's low carry to high word */ \
	xmm4[0] = f64ptr(rsi+1*16);	/* Load low carry word */ \
	f64ptr(rsi+1*16) = xmm1[0];	/* Save prev cache line's high carry in low word */ \
	xmm1[0] = f64ptr(rsi+1*16+8);	/* Load high carry for next cache line */ \
	f64ptr(rsi+1*16+8) = xmm4[0];	/* Move this cache line's low carry to high word */ \
	xmm4[0] = f64ptr(rsi+2*16);	/* Load low carry word */ \
	f64ptr(rsi+2*16) = xmm2[0];	/* Save prev cache line's high carry in low word */ \
	xmm2[0] = f64ptr(rsi+2*16+8);	/* Load high carry for next cache line */ \
	f64ptr(rsi+2*16+8) = xmm4[0];	/* Move this cache line's low carry to high word */ \
	xmm4[0] = f64ptr(rsi+3*16);	/* Load low carry word */ \
	f64ptr(rsi+3*16) = xmm3[0];	/* Save prev cache line's high carry in low word */ \
	xmm3[0] = f64ptr(rsi+3*16+8);	/* Load high carry for next cache line */ \
	f64ptr(rsi+3*16+8) = xmm4[0];	/* Move this cache line's low carry to high word */ \
	rsi += 64;			/* Next carry cache line */ \
	loopZ--;			/* Decrement count of cache lines */ \
}while(loopZ); \
	xmm3[0] -= XMM_BIGVAL1; \
	xmm3[0] *= g->u.xmm.XMM_MINUS_C[0];	/* Negate the very last carry */ \
	xmm3[0] += XMM_BIGVAL1; \
	rsi = (uintptr_t)g->carries;		/* Reload carries array pointer */ \
	f64ptr(rsi+0*16) = xmm3[0];	/* Move last cache line's high carries into first cache line */ \
	f64ptr(rsi+1*16) = xmm0[0]; \
	f64ptr(rsi+2*16) = xmm1[0]; \
	f64ptr(rsi+3*16) = xmm2[0]; \
} \


#define xnorm012_wpn(base2) \
\
;	/* If k or c is more than one, then there will be fewer bits-per-word. */ \
;	/* This means the carry may need to be spread over 4 words instead */ \
;	/* of just 2. */ \
\
if(g->SPREAD_CARRY_OVER_EXTRA_WORDS != 1){ /* Are there few bits per word? */ \
	/*je hard;*/			/* Yes, go do it the hard way */ \
\
	rax = u8ptr(rdi+1);	/* Load 4 big vs. little flags */ \
	rcx = u8ptr(rdi);	/* Load fudge factor flags 1,2 */ \
	rcx &= 0xF0; \
	xmm4 = xptr(rbp+0*16);	/* FFT data */ \
	xmm4 *= xptr(rdx+0*XMM_GMD+rcx); /* mul by fudged grp two-to-minus-phi */ \
	xmm0 = xptr(rsi+0*16);	/* Load carry */ \
	xmm0 += xmm4;		/* carry + FFT data */ \
	rounding(exec, base2, noexec, noexec, xmm0, xmm5, xmm4, rax*4); \
	xmm5 -= XMM_BIGVAL2; \
	xmm0 *= xptr(rdx+0*XMM_GMD+XMM_GMD/2+rcx); /* Value *= fudged grp two-to-phi */ \
	rbx = u8ptr(rdi+2);	/* Load fudge factor flags 5,6 */ \
	rbx &= 0xF0; \
	xmm5 *= xptr(rdx+0*XMM_GMD+XMM_GMD/2+rbx); /* high carry *= fudged grp two-to-phi */ \
	xmm5 += xptr(rbp+4*16);	/* Add high carry and FFT data */ \
	xptr(rbp+0*16) = xmm0;	/* Save FFT data */ \
	xptr(rbp+4*16) = xmm5;	/* Save FFT data */ \
\
	xmm4 = xptr(rbp+1*16);	/* FFT data */ \
	xmm4 *= xptr(rdx+1*XMM_GMD+rcx); /* mul by fudged grp two-to-minus-phi */ \
	xmm1 = xptr(rsi+1*16);	/* Load carry */ \
	xmm1 += xmm4;		/* carry + FFT data */ \
	rounding(exec, base2, noexec, noexec, xmm1, xmm5, xmm4, rax*4+16); \
	xmm5 -= XMM_BIGVAL2; \
	xmm1 *= xptr(rdx+1*XMM_GMD+XMM_GMD/2+rcx); /* value *= fudged grp two-to-phi */ \
	xmm5 *= xptr(rdx+1*XMM_GMD+XMM_GMD/2+rbx); /* high carry *= fudged grp two-to-phi */ \
	xmm5 += xptr(rbp+5*16);	/* Add high carry and FFT data */ \
	xptr(rbp+1*16) = xmm1;	/* Save FFT data */ \
	xptr(rbp+5*16) = xmm5;	/* Save FFT data */ \
\
	rcx = u8ptr(rdi);	/* Load fudge factor flags 3,4 */ \
	rcx &= 0x0F; \
	xmm4 = xptr(rbp+2*16);	/* FFT data */ \
	xmm4 *= xptr(rdx+2*XMM_GMD+rcx*8); /* mul by fudged grp two-to-minus-phi */ \
	xmm2 = xptr(rsi+2*16);	/* Load carry */ \
	xmm2 += xmm4;		/* carry + FFT data */ \
	rounding(exec, base2, noexec, noexec, xmm2, xmm5, xmm4, rax*4+32); \
	xmm5 -= XMM_BIGVAL2; \
	xmm2 *= xptr(rdx+2*XMM_GMD+XMM_GMD/2+rcx*8); /* value *= fudged grp two-to-phi */ \
	rbx = u8ptr(rdi+2);	/* Load fudge factor flags 7,8 */ \
	rbx &= 0x0F; \
	xmm5 *= xptr(rdx+2*XMM_GMD+XMM_GMD/2+rbx*8); /* high carry *= fudged grp two-to-phi */ \
	xmm5 += xptr(rbp+6*16);	/* Add high carry and FFT data */ \
	xptr(rbp+2*16) = xmm2;	/* Save FFT data */ \
	xptr(rbp+6*16) = xmm5;	/* Save FFT data */ \
\
	xmm4 = xptr(rbp+3*16);	/* FFT data */ \
	xmm4 *= xptr(rdx+3*XMM_GMD+rcx*8); /* mul by fudged grp two-to-minus-phi */ \
	xmm3 = xptr(rsi+3*16);	/* Load carry */ \
	xmm3 += xmm4;		/* carry + FFT data */ \
	rounding(exec, base2, noexec, noexec, xmm3, xmm5, xmm4, rax*4+48); \
	xmm5 -= XMM_BIGVAL2; \
	xmm3 *= xptr(rdx+3*XMM_GMD+XMM_GMD/2+rcx*8); /* value *= fudged grp two-to-phi */ \
	xmm5 *= xptr(rdx+3*XMM_GMD+XMM_GMD/2+rbx*8); /* high carry *= fudged grp two-to-phi */ \
	xmm5 += xptr(rbp+7*16);	/* Add high carry and FFT data */ \
	xptr(rbp+3*16) = xmm3;	/* Save FFT data */ \
	xptr(rbp+7*16) = xmm5;	/* Save FFT data */ \
\
}else{ \
	/* Same as above, but spread carry over 6 words */ \
	/*hard:*/ \
	rcx = u8ptr(rdi);	/* Load fudge factor flags 1,2 */ \
	rcx &= 0xF0; \
	rax = u8ptr(rdi+1);	/* Load 4 big vs. little flags */ \
	xmm4 = xptr(rbp+0*16);	/* FFT data */ \
	xmm4 *= xptr(rdx+0*XMM_GMD+rcx); /* mul by fudged grp two-to-minus-phi */ \
	xmm0 = xptr(rsi+0*16);	/* Load carry */ \
	xmm0 += xmm4;		/* carry + FFT data */ \
	rounding(exec, base2, noexec, noexec, xmm0, xmm5, xmm4, rax*4); \
	xmm0 *= xptr(rdx+0*XMM_GMD+XMM_GMD/2+rcx); /* mul by fudged grp two-to-phi */ \
	xptr(rbp+0*16) = xmm0;	/* Save FFT data */ \
\
	rcx = u8ptr(rdi+2);	/* Load fudge factor flags */ \
	rcx &= 0xF0; \
	rax = u8ptr(rdi+3);	/* Load big vs. little flag */ \
	xmm4 = xptr(rbp+4*16);	/* FFT data */ \
	xmm4 *= xptr(rdx+0*XMM_GMD+rcx); /* mul by fudged grp two-to-minus-phi */ \
	xmm5 += xmm4;		/* carry + FFT data */ \
	rounding(exec, base2, noexec, noexec, xmm5, xmm0, xmm4, rax*4); \
	xmm5 *= xptr(rdx+0*XMM_GMD+XMM_GMD/2+rcx); /* mul by fudged grp two-to-phi */ \
	xptr(rbp+4*16) = xmm5;	/* Save FFT data */ \
\
	rbx = g->BIGLIT_INCR2;/* Different clm values step through */ \
;					/* big/lit array differently */ \
	rcx = u8ptr(rdi+rbx);	/* Load fudge factor flags */ \
	rcx &= 0xF0; \
	rax = u8ptr(rdi+rbx+1); /* Load 4 big vs. little flags */ \
	xmm4 = xptr(rbp+8*16);	/* FFT data */ \
	xmm4 *= xptr(rdx+0*XMM_GMD+rcx); /* mul by fudged grp two-to-minus-phi */ \
	xmm0 += xmm4;		/* carry + FFT data */ \
	rounding(exec, base2, noexec, noexec, xmm0, xmm5, xmm4, rax*4); \
	xmm0 *= xptr(rdx+0*XMM_GMD+XMM_GMD/2+rcx); /* mul by fudged grp two-to-phi */ \
	xptr(rbp+8*16) = xmm0;	/* Save FFT data */ \
\
	rcx = u8ptr(rdi+rbx+2); /* Load fudge factor flags */ \
	rcx &= 0xF0; \
	rax = u8ptr(rdi+rbx+3); /* Load big vs. little flag */ \
	xmm4 = xptr(rbp+12*16);	/* FFT data */ \
	xmm4 *= xptr(rdx+0*XMM_GMD+rcx); /* mul by fudged grp two-to-minus-phi */ \
	xmm5 += xmm4;		/* carry + FFT data */ \
	rounding(exec, base2, noexec, noexec, xmm5, xmm0, xmm4, rax*4); \
	xmm5 *= xptr(rdx+0*XMM_GMD+XMM_GMD/2+rcx); /* mul by fudged grp two-to-phi */ \
	xptr(rbp+12*16) = xmm5;	/* Save FFT data */ \
\
	rbx = g->BIGLIT_INCR4;/* Different clm values step through */ \
;					/* big/lit array differently */ \
	rcx = u8ptr(rdi+rbx);	/* Load fudge factor flags */ \
	rcx &= 0xF0; \
	rax = u8ptr(rdi+rbx+1); /* Load 4 big vs. little flags */ \
	xmm4 = xptr(rbp+16*16);	/* FFT data */ \
	xmm4 *= xptr(rdx+0*XMM_GMD+rcx); /* mul by fudged grp two-to-minus-phi */ \
	xmm0 += xmm4;		/* carry + FFT data */ \
	rounding(exec, base2, noexec, noexec, xmm0, xmm5, xmm4, rax*4); \
	xmm0 *= xptr(rdx+0*XMM_GMD+XMM_GMD/2+rcx); /* mul by fudged grp two-to-phi */ \
	xptr(rbp+16*16) = xmm0;	/* Save FFT data */ \
\
	rcx = u8ptr(rdi+rbx+2); /* Load fudge factor flags */ \
	rcx &= 0xF0; \
	xmm5 -= XMM_BIGVAL2; \
	xmm5 *= xptr(rdx+0*XMM_GMD+XMM_GMD/2+rcx); /* high carry *= fudged grp two-to-phi */ \
	xmm5 += xptr(rbp+20*16);	/* Add high carry and FFT data */ \
	xptr(rbp+20*16) = xmm5;	/* Save FFT data */ \
\
\
	rcx = u8ptr(rdi);	/* Load fudge factor flags 1,2 */ \
	rcx &= 0xF0; \
	rax = u8ptr(rdi+1);	/* Load 4 big vs. little flags */ \
	xmm4 = xptr(rbp+1*16);	/* FFT data */ \
	xmm4 *= xptr(rdx+1*XMM_GMD+rcx); /* mul by fudged grp two-to-minus-phi */ \
	xmm1 = xptr(rsi+1*16);	/* Load carry */ \
	xmm1 += xmm4;		/* carry + FFT data */ \
	rounding(exec, base2, noexec, noexec, xmm1, xmm5, xmm4, rax*4+16); \
	xmm1 *= xptr(rdx+1*XMM_GMD+XMM_GMD/2+rcx); /* mul by fudged grp two-to-phi */ \
	xptr(rbp+1*16) = xmm1;	/* Save FFT data */ \
\
	rcx = u8ptr(rdi+2);	/* Load fudge factor flags */ \
	rcx &= 0xF0; \
	rax = u8ptr(rdi+3);	/* Load big vs. little flag */ \
	xmm4 = xptr(rbp+5*16);	/* FFT data */ \
	xmm4 *= xptr(rdx+1*XMM_GMD+rcx); /* mul by fudged grp two-to-minus-phi */ \
	xmm5 += xmm4;		/* carry + FFT data */ \
	rounding(exec, base2, noexec, noexec, xmm5, xmm1, xmm4, rax*4+16); \
	xmm5 *= xptr(rdx+1*XMM_GMD+XMM_GMD/2+rcx); /* mul by fudged grp two-to-phi */ \
	xptr(rbp+5*16) = xmm5;	/* Save FFT data */ \
\
	rbx = g->BIGLIT_INCR2;/* Different clm values step through */ \
;					/* big/lit array differently */ \
	rcx = u8ptr(rdi+rbx);	/* Load fudge factor flags */ \
	rcx &= 0xF0; \
	rax = u8ptr(rdi+rbx+1); /* Load 4 big vs. little flags */ \
	xmm4 = xptr(rbp+9*16);	/* FFT data */ \
	xmm4 *= xptr(rdx+1*XMM_GMD+rcx); /* mul by grp two-to-minus-phi */ \
	xmm1 += xmm4;		/* carry + FFT data */ \
	rounding(exec, base2, noexec, noexec, xmm1, xmm5, xmm4, rax*4+16); \
	xmm1 *= xptr(rdx+1*XMM_GMD+XMM_GMD/2+rcx); /* mul by fudged grp two-to-phi */ \
	xptr(rbp+9*16) = xmm1;	/* Save FFT data */ \
\
	rcx = u8ptr(rdi+rbx+2); /* Load fudge factor flags */ \
	rcx &= 0xF0; \
	rax = u8ptr(rdi+rbx+3); /* Load big vs. little flag */ \
	xmm4 = xptr(rbp+13*16);	/* FFT data */ \
	xmm4 *= xptr(rdx+1*XMM_GMD+rcx); /* mul by fudged grp two-to-minus-phi */ \
	xmm5 += xmm4;		/* carry + FFT data */ \
	rounding(exec, base2, noexec, noexec, xmm5, xmm1, xmm4, rax*4+16); \
	xmm5 *= xptr(rdx+1*XMM_GMD+XMM_GMD/2+rcx); /* mul by fudged grp two-to-phi */ \
	xptr(rbp+13*16) = xmm5;	/* Save FFT data */ \
\
	rbx = g->BIGLIT_INCR4;/* Different clm values step through */ \
;					/* big/lit array differently */ \
	rcx = u8ptr(rdi+rbx);	/* Load fudge factor flags */ \
	rcx &= 0xF0; \
	rax = u8ptr(rdi+rbx+1); /* Load 4 big vs. little flags */ \
	xmm4 = xptr(rbp+17*16);	/* FFT data */ \
	xmm4 *= xptr(rdx+1*XMM_GMD+rcx); /* mul by grp two-to-minus-phi */ \
	xmm1 += xmm4;		/* carry + FFT data */ \
	rounding(exec, base2, noexec, noexec, xmm1, xmm5, xmm4, rax*4+16); \
	xmm1 *= xptr(rdx+1*XMM_GMD+XMM_GMD/2+rcx); /* mul by fudged grp two-to-phi */ \
	xptr(rbp+17*16) = xmm1;	/* Save FFT data */ \
\
	rcx = u8ptr(rdi+rbx+2); /* Load fudge factor flags */ \
	rcx &= 0xF0; \
	xmm5 -= XMM_BIGVAL2; \
	xmm5 *= xptr(rdx+1*XMM_GMD+XMM_GMD/2+rcx); /* high carry *= fudged grp two-to-phi */ \
	xmm5 += xptr(rbp+21*16);	/* Add high carry and FFT data */ \
	xptr(rbp+21*16) = xmm5;	/* Save FFT data */ \
\
	rcx = u8ptr(rdi);	/* Load fudge factor flags 3,4 */ \
	rcx &= 0x0F; \
	rax = u8ptr(rdi+1);	/* Load 4 big vs. little flags */ \
	xmm4 = xptr(rbp+2*16);	/* FFT data */ \
	xmm4 *= xptr(rdx+2*XMM_GMD+rcx*8); /* mul by fudged grp two-to-minus-phi */ \
	xmm2 = xptr(rsi+2*16);	/* Load carry */ \
	xmm2 += xmm4;		/* carry + FFT data */ \
	rounding(exec, base2, noexec, noexec, xmm2, xmm5, xmm4, rax*4+32); \
	xmm2 *= xptr(rdx+2*XMM_GMD+XMM_GMD/2+rcx*8); /* mul by fudged grp two-to-phi */ \
	xptr(rbp+2*16) = xmm2;	/* Save FFT data */ \
\
	rcx = u8ptr(rdi+2);	/* Load fudge factor flags */ \
	rcx &= 0x0F; \
	rax = u8ptr(rdi+3);	/* Load big vs. little flag */ \
	xmm4 = xptr(rbp+6*16);	/* FFT data */ \
	xmm4 *= xptr(rdx+2*XMM_GMD+rcx*8); /* mul by fudged grp two-to-minus-phi */ \
	xmm5 += xmm4;		/* carry + FFT data */ \
	rounding(exec, base2, noexec, noexec, xmm5, xmm2, xmm4, rax*4+32); \
	xmm5 *= xptr(rdx+2*XMM_GMD+XMM_GMD/2+rcx*8); /* mul by fudged grp two-to-phi */ \
	xptr(rbp+6*16) = xmm5;	/* Save FFT data */ \
\
	rbx = g->BIGLIT_INCR2;/* Different clm values step through */ \
;					/* big/lit array differently */ \
	rcx = u8ptr(rdi+rbx);	/* Load fudge factor flags */ \
	rcx &= 0x0F; \
	rax = u8ptr(rdi+rbx+1); /* Load 4 big vs. little flags */ \
	xmm4 = xptr(rbp+10*16);	/* FFT data */ \
	xmm4 *= xptr(rdx+2*XMM_GMD+rcx*8); /* mul by fudged grp two-to-minus-phi */ \
	xmm2 += xmm4;		/* carry + FFT data */ \
	rounding(exec, base2, noexec, noexec, xmm2, xmm5, xmm4, rax*4+32); \
	xmm2 *= xptr(rdx+2*XMM_GMD+XMM_GMD/2+rcx*8); /* mul by fudged grp two-to-phi */ \
	xptr(rbp+10*16) = xmm2;	/* Save FFT data */ \
\
	rcx = u8ptr(rdi+rbx+2); /* Load fudge factor flags */ \
	rcx &= 0x0F; \
	rax = u8ptr(rdi+rbx+3); /* Load big vs. little flag */ \
	xmm4 = xptr(rbp+14*16);	/* FFT data */ \
	xmm4 *= xptr(rdx+2*XMM_GMD+rcx*8); /* mul by fudged grp two-to-minus-phi */ \
	xmm5 += xmm4;		/* carry + FFT data */ \
	rounding(exec, base2, noexec, noexec, xmm5, xmm2, xmm4, rax*4+32); \
	xmm5 *= xptr(rdx+2*XMM_GMD+XMM_GMD/2+rcx*8); /* mul by fudged grp two-to-phi */ \
	xptr(rbp+14*16) = xmm5;	/* Save FFT data */ \
\
	rbx = g->BIGLIT_INCR4;/* Different clm values step through */ \
;					/* big/lit array differently */ \
	rcx = u8ptr(rdi+rbx);	/* Load fudge factor flags */ \
	rcx &= 0x0F; \
	rax = u8ptr(rdi+rbx+1); /* Load 4 big vs. little flags */ \
	xmm4 = xptr(rbp+18*16);	/* FFT data */ \
	xmm4 *= xptr(rdx+2*XMM_GMD+rcx*8); /* mul by fudged grp two-to-minus-phi */ \
	xmm2 += xmm4;		/* carry + FFT data */ \
	rounding(exec, base2, noexec, noexec, xmm2, xmm5, xmm4, rax*4+32); \
	xmm2 *= xptr(rdx+2*XMM_GMD+XMM_GMD/2+rcx*8); /* mul by fudged grp two-to-phi */ \
	xptr(rbp+18*16) = xmm2;	/* Save FFT data */ \
\
	rcx = u8ptr(rdi+rbx+2); /* Load fudge factor flags */ \
	rcx &= 0x0F; \
	xmm5 -= XMM_BIGVAL2; \
	xmm5 *= xptr(rdx+2*XMM_GMD+XMM_GMD/2+rcx*8); /* high carry *= fudged grp two-to-phi */ \
	xmm5 += xptr(rbp+22*16);	/* Add high carry and FFT data */ \
	xptr(rbp+22*16) = xmm5;	/* Save FFT data */ \
\
\
	rcx = u8ptr(rdi);	/* Load fudge factor flags 3,4 */ \
	rcx &= 0x0F; \
	rax = u8ptr(rdi+1);	/* Load 4 big vs. little flags */ \
	xmm4 = xptr(rbp+3*16);	/* FFT data */ \
	xmm4 *= xptr(rdx+3*XMM_GMD+rcx*8); /* mul by fudged grp two-to-minus-phi */ \
	xmm3 = xptr(rsi+3*16);	/* Load carry */ \
	xmm3 += xmm4;		/* carry + FFT data */ \
	rounding(exec, base2, noexec, noexec, xmm3, xmm5, xmm4, rax*4+48); \
	xmm3 *= xptr(rdx+3*XMM_GMD+XMM_GMD/2+rcx*8); /* mul by fudged grp two-to-phi */ \
	xptr(rbp+3*16) = xmm3;	/* Save FFT data */ \
\
	rcx = u8ptr(rdi+2);	/* Load fudge factor flags */ \
	rcx &= 0x0F; \
	rax = u8ptr(rdi+3);	/* Load big vs. little flag */ \
	xmm4 = xptr(rbp+7*16);	/* FFT data */ \
	xmm4 *= xptr(rdx+3*XMM_GMD+rcx*8); /* mul by fudged grp two-to-minus-phi */ \
	xmm5 += xmm4;		/* carry + FFT data */ \
	rounding(exec, base2, noexec, noexec, xmm5, xmm3, xmm4, rax*4+48); \
	xmm5 *= xptr(rdx+3*XMM_GMD+XMM_GMD/2+rcx*8); /* mul by fudged grp two-to-phi */ \
	xptr(rbp+7*16) = xmm5;	/* Save FFT data */ \
\
	rbx = g->BIGLIT_INCR2;/* Different clm values step through */ \
;					/* big/lit array differently */ \
	rcx = u8ptr(rdi+rbx);	/* Load fudge factor flags */ \
	rcx &= 0x0F; \
	rax = u8ptr(rdi+rbx+1); /* Load 4 big vs. little flags */ \
	xmm4 = xptr(rbp+11*16);	/* FFT data */ \
	xmm4 *= xptr(rdx+3*XMM_GMD+rcx*8); /* mul by fudged grp two-to-minus-phi */ \
	xmm3 += xmm4;		/* carry + FFT data */ \
	rounding(exec, base2, noexec, noexec, xmm3, xmm5, xmm4, rax*4+48); \
	xmm3 *= xptr(rdx+3*XMM_GMD+XMM_GMD/2+rcx*8); /* mul by fudged grp two-to-phi */ \
	xptr(rbp+11*16) = xmm3;	/* Save FFT data */ \
\
	rcx = u8ptr(rdi+rbx+2); /* Load fudge factor flags */ \
	rcx &= 0x0F; \
	rax = u8ptr(rdi+rbx+3); /* Load big vs. little flag */ \
	xmm4 = xptr(rbp+15*16);	/* FFT data */ \
	xmm4 *= xptr(rdx+3*XMM_GMD+rcx*8); /* mul by fudged grp two-to-minus-phi */ \
	xmm5 += xmm4;		/* carry + FFT data */ \
	rounding(exec, base2, noexec, noexec, xmm5, xmm3, xmm4, rax*4+48); \
	xmm5 *= xptr(rdx+3*XMM_GMD+XMM_GMD/2+rcx*8); /* mul by fudged grp two-to-phi */ \
	xptr(rbp+15*16) = xmm5;	/* Save FFT data */ \
\
	rbx = g->BIGLIT_INCR4;/* Different clm values step through */ \
;					/* big/lit array differently */ \
	rcx = u8ptr(rdi+rbx);	/* Load fudge factor flags */ \
	rcx &= 0x0F; \
	rax = u8ptr(rdi+rbx+1); /* Load 4 big vs. little flags */ \
	xmm4 = xptr(rbp+19*16);	/* FFT data */ \
	xmm4 *= xptr(rdx+3*XMM_GMD+rcx*8); /* mul by fudged grp two-to-minus-phi */ \
	xmm3 += xmm4;		/* carry + FFT data */ \
	rounding(exec, base2, noexec, noexec, xmm3, xmm5, xmm4, rax*4+48); \
	xmm3 *= xptr(rdx+3*XMM_GMD+XMM_GMD/2+rcx*8); /* mul by fudged grp two-to-phi */ \
	xptr(rbp+19*16) = xmm3;	/* Save FFT data */ \
\
	rcx = u8ptr(rdi+rbx+2); /* Load fudge factor flags */ \
	rcx &= 0x0F; \
	xmm5 -= XMM_BIGVAL2; \
	xmm5 *= xptr(rdx+3*XMM_GMD+XMM_GMD/2+rcx*8); /* high carry *= fudged grp two-to-phi */ \
	xmm5 += xptr(rbp+23*16);	/* Add high carry and FFT data */ \
	xptr(rbp+23*16) = xmm5;	/* Save FFT data */ \
} \
	xmm4 = XMM_BIGVAL2; \
	xptr(rsi+0*16) = xmm4;	/* Clear carry */ \
	xptr(rsi+1*16) = xmm4; \
	xptr(rsi+2*16) = xmm4; \
	xptr(rsi+3*16) = xmm4; \



/* Significantly different cleanup code for zero-padded FFTs. */ \
/* Note: The group multiplier should be 1.0 for the bottom FFT words and */ \
/* the FFT words just above the half-way point. */ \

#define xnorm012_wpn_zpad_part1 \
if(g->this_block == 0) {		/* Are we carrying into the first data block? */ \
	/*jne done;*/			/* If not, skip wrapping carries and negating the last carry */ \
	rax = g->addcount1;		/* Load count of cache lines in the carries array */ \
	do{ \
		xmm4[0] = f64ptr(rsi+0*16);	/* Load low carry word */ \
		f64ptr(rsi+0*16) = xmm0[0];	/* Save prev cache line's high carry in low word */ \
		xmm0[0] = f64ptr(rsi+0*16+8);	/* Load high carry for next cache line */ \
		f64ptr(rsi+0*16+8) = xmm4[0];	/* Move this cache line's low carry to high word */ \
		xmm4[0] = f64ptr(rsi+1*16);	/* Load low carry word */ \
		f64ptr(rsi+1*16) = xmm2[0];	/* Save prev cache line's high carry in low word */ \
		xmm2[0] = f64ptr(rsi+1*16+8);	/* Load high carry for next cache line */ \
		f64ptr(rsi+1*16+8) = xmm4[0];	/* Move this cache line's low carry to high word */ \
		xmm4[0] = f64ptr(rsi+2*16);	/* Load low carry word */ \
		f64ptr(rsi+2*16) = xmm1[0];	/* Save prev cache line's high carry in low word */ \
		xmm1[0] = f64ptr(rsi+2*16+8);	/* Load high carry for next cache line */ \
		f64ptr(rsi+2*16+8) = xmm4[0];	/* Move this cache line's low carry to high word */ \
		xmm4[0] = f64ptr(rsi+3*16);	/* Load low carry word */ \
		f64ptr(rsi+3*16) = xmm3[0];	/* Save prev cache line's high carry in low word */ \
		xmm3[0] = f64ptr(rsi+3*16+8);	/* Load high carry for next cache line */ \
		f64ptr(rsi+3*16+8) = xmm4[0];	/* Move this cache line's low carry to high word */ \
		rsi += 64;			/* Next carry cache line */ \
		rax--;			/* Decrement count of cache lines */ \
	}while(rax); \
	rsi = (uintptr_t)g->carries;		/* Reload carries array pointer */ \
/*	movsd	Q [rsi+0*16], xmm3	;; Move last cache line's high carries into first cache line */ \
	f64ptr(rsi+1*16) = xmm0[0]; \
/*	movsd	Q [rsi+2*16], xmm2 */ \
	f64ptr(rsi+3*16) = xmm1[0]; \
\
	rsi = (uintptr_t)g->DESTARG;/* Addr of the FFT data */ \
	rdi = (uintptr_t)g->norm_biglit_array;	/* Addr of the big/little flags array */ \
\
	if(g->B_IS_2 == 0){		/* Is b = 2? */ \
		xnorm012_wpn_zpad_part1a(noexec); /* No, do harder rounding */ \
	}else{ \
		xnorm012_wpn_zpad_part1a(exec); /* Yes, do simpler rounding */ \
	} \
	rsi = (uintptr_t)g->carries;		/* Reload carries pointer */ \
	xmm0[0] = XMM_BIGVAL1;/* Clear two carries just processed by xnorm012_wpn_zpad_part1_cmn */ \
	xmm2[0] = 0.0; \
	f64ptr(rsi) = xmm0[0]; \
	f64ptr(rsi+32) = xmm2[0]; \
} \


#define xnorm012_wpn_zpad_part1a(base2) \
if(g->const_fft != 0) { /* Are we also multiplying by a constant? */ \
	xnorm012_wpn_zpad_part1_cmn(exec, base2); \
}else{/* Jump if not const */ \
	xnorm012_wpn_zpad_part1_cmn(noexec, base2); \
} \
\
\
/* On input, xmm2 and xmm3 contain high carries from the last carries array row */ \

#define xnorm012_wpn_zpad_part1_cmn(const1, base2) { \
\
;	/* Strip BIGVAL from the traditional carry, we'll add the traditional */ \
;	/* carry in later when we are working on the ZPAD0 - ZPAD6 values. */ \
	xmm2[0] -= XMM_BIGVAL1;/* Integerize traditional carry */ \
\
;	/* Rather than calculate high FFT carry times k and then later dividing */ \
;	/* by k, we multiply FFT high carry by const and we'll add it */ \
;	/* to the lower FFT data later (after multiplying by -c). */ \
const1(xmm3[0] *= g->u.xmm.XMM_MULCONST[0]); \
\
;	/* Work on zero-pad addin value. */ \
	xmm4[0] = g->ADDIN_VALUE;		/* Load the add in value */ \
const1(xmm4[0] *= g->u.xmm.XMM_MULCONST[0]);		/* Multiply the add in value by the small mul const */ \
	xmm4[0] += g->POSTADDIN_VALUE;		/* Add the post-mul-by-const addin value */ \
;	/* when c = -1, 1 = b^n */ \
;	/* when c = 1, -1 = b^n, 1 = -b^n */ \
;	/* when c = -3, 3 = b^n, 1 = b^n - 2 */ \
;	/* when c = 3, -3 = b^n, 3 = -b^n, 1 = -b^n - 2 */ \
;	/* The "- 2" has been precomputed in ZPAD_LSW_ADJUST.  Add ADDIN_VALUE * mul-by-const * ZPAD_LSW_ADJUST into the least significant FFT word */ \
	rsi = (uintptr_t)g->DESTARG;/* Address of squared number */ \
	xmm0[0] = g->ZPAD_LSW_ADJUST; \
	xmm0[0] *= xmm4[0]; \
	xmm0[0] += f64ptr(rsi); \
	f64ptr(rsi) = xmm0[0]; \
\
;	/* Multiply ZPAD0 through ZPAD6 by const * -C.  This, in essense, */ \
;	/* wraps this data from above the FFT data area to the halfway point. */ \
;	/* Later on we'll divide this by K to decide which data needs wrapping */ \
;	/* all the way down to the bottom of the FFT data. */ \
\
;	/* NOTE: ZPAD0's grp multiplier is 1.0.  Also, ZPAD6 will not */ \
;	/* be bigger than a big word.  We must be careful to handle c's up */ \
;	/* to about 30 bits */ \
\
	rax = u8ptr(rdi+1);	/* Load big vs. little flags */ \
	xmm0[0] = g->ZPAD0_6[0];/* Load values1 */ \
	xmm5[0] = 0.0;		/* Create a zero high FFT carry to add in */ \
	single_split_lower_zpad_word(base2, xmm0[0], xmm5[0], xmm1[0], rax*4); \
no##const1(xmm0[0] *= g->u.xmm.XMM_MINUS_C[0]); \
const1(xmm0[0] *= g->u.xmm.XMM_MINUS_C_TIMES_MULCONST[0]); \
	xmm0[0] += xmm4[0];		/* Apply rest of ADDIN_VALUE * mul-by-const here (the FFT half way point) */ \
	round_zpad7_word(base2, xmm0[0], xmm2[0], xmm1[0], rax*4); \
	g->ZPAD0_6[0] = xmm0[0]; \
\
	rax = u8ptr(rdi+3);	/* Load big vs. little flags */ \
	xmm0[0] = g->ZPAD0_6[1];/* Load values1 */ \
	single_split_lower_zpad_word(base2, xmm0[0], xmm5[0], xmm1[0], rax*4); \
no##const1(xmm0[0] *= g->u.xmm.XMM_MINUS_C[0]); \
const1(xmm0[0] *= g->u.xmm.XMM_MINUS_C_TIMES_MULCONST[0]); \
	round_zpad7_word(base2, xmm0[0], xmm2[0], xmm1[0], rax*4); \
	g->ZPAD0_6[1] = xmm0[0]; \
\
	rcx = g->BIGLIT_INCR2;/* Different clm values step through */ \
;					/* big/lit array differently */ \
	rax = u8ptr(rdi+rcx+1); /* Load big vs. little flags */ \
	xmm0[0] = g->ZPAD0_6[2];/* Load values1 */ \
	single_split_lower_zpad_word(base2, xmm0[0], xmm5[0], xmm1[0], rax*4); \
no##const1(xmm0[0] *= g->u.xmm.XMM_MINUS_C[0]); \
const1(xmm0[0] *= g->u.xmm.XMM_MINUS_C_TIMES_MULCONST[0]); \
	round_zpad7_word(base2, xmm0[0], xmm2[0], xmm1[0], rax*4); \
	g->ZPAD0_6[2] = xmm0[0]; \
\
	rax = u8ptr(rdi+rcx+3); /* Load big vs. little flags */ \
	xmm0[0] = g->ZPAD0_6[3];/* Load values1 */ \
	single_split_lower_zpad_word(base2, xmm0[0], xmm5[0], xmm1[0], rax*4); \
no##const1(xmm0[0] *= g->u.xmm.XMM_MINUS_C[0]); \
const1(xmm0[0] *= g->u.xmm.XMM_MINUS_C_TIMES_MULCONST[0]); \
	round_zpad7_word(base2, xmm0[0], xmm2[0], xmm1[0], rax*4); \
	g->ZPAD0_6[3] = xmm0[0]; \
\
	rcx = g->BIGLIT_INCR4;/* Different clm values step through */ \
;					/* big/lit array differently */ \
	rax = u8ptr(rdi+rcx+1); /* Load big vs. little flags */ \
	xmm0[0] = g->ZPAD0_6[4];/* Load values1 */ \
	single_split_lower_zpad_word(base2, xmm0[0], xmm5[0], xmm1[0], rax*4); \
no##const1(xmm0[0] *= g->u.xmm.XMM_MINUS_C[0]); \
const1(xmm0[0] *= g->u.xmm.XMM_MINUS_C_TIMES_MULCONST[0]); \
	round_zpad7_word(base2, xmm0[0], xmm2[0], xmm1[0], rax*4); \
	g->ZPAD0_6[4] = xmm0[0]; \
\
	rax = u8ptr(rdi+rcx+3); /* Load big vs. little flags */ \
	xmm0[0] = g->ZPAD0_6[5];/* Load values1 */ \
	single_split_lower_zpad_word(base2, xmm0[0], xmm5[0], xmm1[0], rax*4); \
no##const1(xmm0[0] *= g->u.xmm.XMM_MINUS_C[0]); \
const1(xmm0[0] *= g->u.xmm.XMM_MINUS_C_TIMES_MULCONST[0]); \
	round_zpad7_word(base2, xmm0[0], xmm2[0], xmm1[0], rax*4); \
	g->ZPAD0_6[5] = xmm0[0]; \
\
	xmm0[0] = g->ZPAD0_6[6];/* Load values1 */ \
	xmm0[0] = roundf64(xmm0[0]); /* Round to an integer */ \
	xmm0[0] += xmm5[0];		/* Add in shifted high ZPAD data */ \
no##const1(xmm0[0] *= g->u.xmm.XMM_MINUS_C[0]); \
const1(xmm0[0] *= g->u.xmm.XMM_MINUS_C_TIMES_MULCONST[0]); \
	xmm0[0] += xmm2[0];		/* Add in high part of last calculation */ \
	g->ZPAD0_6[6] = xmm0[0]; \
\
;	/* Divide the zpad data by k.  Store the integer part in XMM_TMP */ \
;	/* and the remainder in ZPAD0.  Later we will wrap the integer part */ \
;	/* down to the bottom of the FFT data area (and multiply by -c). */ \
;	/* And we will store the remainder in the upper half of the FFT */ \
;	/* data area. */ \
\
;	/* Note there are three cases to handle.  K is smaller than a big word. */ \
;	/* K is between one and 2 big words in size.  And K is more than */ \
;	/* 2 big words in size. */ \
\
	xfunc572(g); \
\
;	/* Now normalize the data above the halfway point.  Remember that the */ \
;	/* two-to-phi multipliers will be applied by the FFT. */ \
\
	rax = u8ptr(rdi+1);	/* Load big vs. little flags */ \
	xmm0[0] = g->ZPAD0_6[0];/* Load remainder of divide by k */ \
	xmm0[0] += XMM_BIGVAL1; \
	single_rounding(base2, xmm0[0], xmm2[0], xmm4[0], rax*4); \
	f64ptr(rsi+0*64+32) = xmm0[0];	/* Save value1 */ \
\
	rax = u8ptr(rdi+3);	/* Load big vs. little flags */ \
	single_rounding(base2, xmm2[0], xmm0[0], xmm4[0], rax*4); \
	f64ptr(rsi+1*64+32) = xmm2[0];	/* Save value2 */ \
\
	rcx = g->BIGLIT_INCR2;/* Different clm values step through */ \
;					/* big/lit array differently */ \
	rax = u8ptr(rdi+rcx+1); /* Load big vs. little flags */ \
	single_rounding(base2, xmm0[0], xmm2[0], xmm4[0], rax*4); \
	f64ptr(rsi+2*64+32) = xmm0[0];	/* Save value3 */ \
\
	xmm2[0] -= XMM_BIGVAL1;/* Remove integer rounding constant */ \
	f64ptr(rsi+3*64+32) = xmm2[0];	/* Save new value4 */ \
\
;	/* Mul the integer part of (ZPAD data divided by k) by -c in */ \
;	/* preparation for adding it into the lower FFT data area. */ \
;	/* Also add in the shifted high FFT carry at this time. */ \
\
;	/* Now add in and normalize the bottom FFT data.  Remember that the */ \
;	/* column two-to-phi multiplier for the first value will be 1.0.  We  */ \
;	/* must go 6 words deep in case k is 48-50 bits and c is 32 bits. */ \
\
	rax = u8ptr(rdi+1);	/* Load big vs. little flags */ \
	vec2f64 TMP12 = g->u.xmm.XMM_TMP1; \
	vec2f64 TMP34 = g->u.xmm.XMM_TMP2; \
	vec2f64 TMP56 = g->u.xmm.XMM_TMP3; \
	TMP34 *= g->u.xmm.XMM_MINUS_C;	/* Mul by -c */ \
	TMP56 *= g->u.xmm.XMM_MINUS_C;	/* Mul by -c */ \
	xmm0[0] = TMP12[0];	/* Load integer part of divide by k */ \
	xmm0[0] += xmm3[0];		/* Add in shifted high FFT carry */ \
	xmm0[0] *= g->u.xmm.XMM_MINUS_C[0];	/* Mul by -c */ \
	xmm0[0] += XMM_BIGVAL1; \
	xmm0[0] += f64ptr(rsi+0*64);	/* Add in the FFT data */ \
	single_rounding(base2, xmm0[0], xmm2[0], xmm4[0], rax*4); \
	f64ptr(rsi+0*64) = xmm0[0];	/* Save value1 */ \
\
	rax = u8ptr(rdi+3);	/* Load big vs. little flags */ \
	xmm0[0] = TMP12[1];	/* Load integer part of divide by k */ \
	xmm0[0] *= g->u.xmm.XMM_MINUS_C[0];	/* Mul by -c */ \
	xmm0[0] += f64ptr(rsi+1*64);	/* Add in the FFT data */ \
	xmm0[0] += xmm2[0];		/* x2 = value + carry */ \
	single_rounding(base2, xmm0[0], xmm2[0], xmm4[0], rax*4); \
	f64ptr(rsi+1*64) = xmm0[0];	/* Save value2 */ \
\
	rcx = rdi + g->BIGLIT_INCR2;/* Different clm values step through */ \
;					/* big/lit array differently */ \
	rax = u8ptr(rcx+1); /* Load big vs. little flags */ \
	xmm0[0] = TMP34[0];	/* Load integer part of divide by k */ \
	xmm0[0] += f64ptr(rsi+2*64);	/* Add in the FFT data */ \
	xmm0[0] += xmm2[0];		/* x3 = value + carry */ \
	single_rounding(base2, xmm0[0], xmm2[0], xmm4[0], rax*4); \
	f64ptr(rsi+2*64) = xmm0[0];	/* Save value3 */ \
\
	rax = u8ptr(rcx+3); /* Load big vs. little flags */ \
	xmm0[0] = TMP34[1];	/* Load integer part of divide by k */ \
	xmm0[0] += f64ptr(rsi+3*64);	/* Add in the FFT data */ \
	xmm0[0] += xmm2[0];		/* x4 = value + carry */ \
	single_rounding(base2, xmm0[0], xmm2[0], xmm4[0], rax*4); \
	f64ptr(rsi+3*64) = xmm0[0];	/* Save value4 */ \
\
	rcx = rdi + g->BIGLIT_INCR4;/* Different clm values step through */ \
;					/* big/lit array differently */ \
	rax = u8ptr(rcx+1); /* Load big vs. little flags */ \
	xmm0[0] = TMP56[0];	/* Load integer part of divide by k */ \
	xmm0[0] += f64ptr(rsi+4*64);	/* Add in the FFT data */ \
	xmm0[0] += xmm2[0];		/* x4 = value + carry */ \
	single_rounding(base2, xmm0[0], xmm2[0], xmm4[0], rax*4); \
	f64ptr(rsi+4*64) = xmm0[0];	/* Save value4 */ \
\
	rax = u8ptr(rcx+3); /* Load big vs. little flags */ \
	xmm0[0] = TMP56[1];	/* Load integer part of divide by k */ \
	xmm0[0] += f64ptr(rsi+5*64);	/* Add in the FFT data */ \
	xmm0[0] += xmm2[0];		/* x5 = value + carry */ \
	single_rounding(base2, xmm0[0], xmm2[0], xmm4[0], rax*4); \
	f64ptr(rsi+5*64) = xmm0[0];	/* Save value5 */ \
\
	rcx += g->BIGLIT_INCR2;/* Different clm values step through */ \
;					/* big/lit array differently */ \
	rax = u8ptr(rcx+1); /* Load big vs. little flags */ \
	xmm0[0] = f64ptr(rsi+6*64);	/* Load FFT data */ \
	xmm0[0] += xmm2[0];		/* x6 = value + carry */ \
	single_rounding(base2, xmm0[0], xmm2[0], xmm4[0], rax*4); \
	f64ptr(rsi+6*64) = xmm0[0];	/* Save value6 */ \
\
	xmm2[0] -= XMM_BIGVAL1;/* Remove rounding constant */ \
	xmm2[0] += f64ptr(rsi+7*64);	/* Add in FFT data */ \
	f64ptr(rsi+7*64) = xmm2[0];	/* Save value7 */ \
}

#define xnorm012_wpn_zpad \
	if(g->B_IS_2 == 0){		/* Is b = 2? */ \
		xnorm012_wpn_zpada(noexec);	/* No, do harder rounding */ \
	}else{ \
		xnorm012_wpn_zpada(exec); /* Yes, do simpler rounding */ \
	} \

#define xnorm012_wpn_zpada(base2) \
	if(g->const_fft != 0) {	/* Are we also multiplying by a constant? */ \
		xnorm012_wpn_zpad_cmn(exec, base2); \
	}else{ /* Jump if not const */ \
		xnorm012_wpn_zpad_cmn(noexec, base2); \
	} \

#define xnorm012_wpn_zpad_cmn(const1, base2) \
	rcx = u8ptr(rdi);	/* Load fudge factor flag */ \
	rcx &= 0xF0; \
	xmm0 = xptr(rsi+0*16);	/* Load carry */ \
	xmm2 = xptr(rsi+2*16);	/* Load hi data carry */ \
	xmm4 = xptr(rbp+0*64);	/* Load FFT data */ \
	xmm4 *= xptr(rdx+0*XMM_GMD+rcx); /* mul by fudged grp two-to-minus-phi */ \
	xmm0 += xmm4;		/* x1 = values1 + carry */ \
	rax = u8ptr(rdi+1);	/* Load big vs. little flag */ \
	split_carry_zpad_word(base2, xmm2, xmm6, xmm4, rax*4); \
no##const1(xmm4 = g->u.xmm.XMM_K_LO);		/* Calc high FFT carry times k */ \
const1(xmm4 = g->u.xmm.XMM_K_TIMES_MULCONST_LO); \
	xmm4 *= xmm2;		/* high_FFT_carry * k_lo */ \
	xmm0 += xmm4;		/* x1 = x1 + high_FFT_carry * k_lo */ \
no##const1(xmm5 = g->u.xmm.XMM_K_HI); \
const1(xmm5 = g->u.xmm.XMM_K_TIMES_MULCONST_HI); \
	xmm5 *= xptr2(g->u.xmm.XMM_LIMIT_INVERSE,rax*4); /* shift k_hi */ \
	xmm5 *= xmm2; \
	rounding(exec, base2, noexec, noexec, xmm0, xmm2, xmm4, rax*4); \
	xmm2 += xmm5;		/* Carry += high_FFT_carry * shifted k_hi */ \
	xmm0 *= xptr(rdx+0*XMM_GMD+XMM_GMD/2+rcx); /* mul by fudged grp two-to-phi */ \
	xptr(rbp+0*64) = xmm0;	/* Save FFT data */ \
\
	rcx = u8ptr(rdi+2);	/* Load fudge factor flags */ \
	rcx &= 0xF0; \
	xmm0 = xptr(rbp+1*64);	/* Load FFT data */ \
	xmm0 *= xptr(rdx+0*XMM_GMD+rcx); /* mul by fudged grp two-to-minus-phi */ \
	xmm0 += xmm2;		/* x2 = values1 + carry */ \
	rax = u8ptr(rdi+3);	/* Load big vs. little flag */ \
	split_carry_zpad_word(base2, xmm6, xmm5, xmm4, rax*4); \
no##const1(xmm4 = g->u.xmm.XMM_K_LO);		/* Calc high FFT carry times k */ \
const1(xmm4 = g->u.xmm.XMM_K_TIMES_MULCONST_LO); \
	xmm4 *= xmm6;		/* high_FFT_carry * k_lo */ \
	xmm0 += xmm4;		/* x2 = x2 + high_FFT_carry * k_lo */ \
no##const1(xmm4 = g->u.xmm.XMM_K_HI); \
const1(xmm4 = g->u.xmm.XMM_K_TIMES_MULCONST_HI); \
	xmm4 *= xptr2(g->u.xmm.XMM_LIMIT_INVERSE,rax*4); /* shift k_hi */ \
	xmm6 *= xmm4; \
	rounding(exec, base2, noexec, noexec, xmm0, xmm2, xmm4, rax*4); \
	xmm2 += xmm6;		/* Carry += high_FFT_carry * shifted k_hi */ \
	xmm0 *= xptr(rdx+0*XMM_GMD+XMM_GMD/2+rcx); /* mul by fudged grp two-to-phi */ \
	xptr(rbp+1*64) = xmm0;	/* Save FFT data */ \
\
	rcx = g->BIGLIT_INCR2;/* Different clm values step through */ \
;					/* big/lit array differently */ \
	rax = u8ptr(rdi+rcx+1); /* Load big vs. little flag */ \
	rcx = u8ptr(rdi+rcx);	/* Load fudge factor flag */ \
	rcx &= 0xF0; \
	xmm0 = xptr(rbp+2*64);	/* Load FFT data */ \
	xmm0 *= xptr(rdx+0*XMM_GMD+rcx); /* mul by fudged grp two-to-minus-phi */ \
	xmm0 += xmm2;		/* x3 = values1 + carry */ \
	split_carry_zpad_word(base2, xmm5, xmm6, xmm4, rax*4); \
no##const1(xmm4 = g->u.xmm.XMM_K_LO);		/* Calc high FFT carry times k */ \
const1(xmm4 = g->u.xmm.XMM_K_TIMES_MULCONST_LO); \
	xmm4 *= xmm5;		/* high_FFT_carry * k_lo */ \
	xmm0 += xmm4;		/* x3 = x3 + high_FFT_carry * k_lo */ \
no##const1(xmm4 = g->u.xmm.XMM_K_HI); \
const1(xmm4 = g->u.xmm.XMM_K_TIMES_MULCONST_HI); \
	xmm4 *= xptr2(g->u.xmm.XMM_LIMIT_INVERSE,rax*4); /* shift k_hi */ \
	xmm5 *= xmm4; \
	rounding(exec, base2, noexec, noexec, xmm0, xmm2, xmm4, rax*4); \
	xmm2 += xmm5;		/* Carry += high_FFT_carry * shifted k_hi */ \
	xmm0 *= xptr(rdx+0*XMM_GMD+XMM_GMD/2+rcx); /* mul by fudged grp two-to-phi */ \
	xptr(rbp+2*64) = xmm0;	/* Save FFT data */ \
\
	rcx = g->BIGLIT_INCR2;/* Different clm values step through */ \
;					/* big/lit array differently */ \
	rax = u8ptr(rdi+rcx+3); /* Load big vs. little flag */ \
	rcx = u8ptr(rdi+rcx+2); /* Load fudge factor flag */ \
	rcx &= 0xF0; \
	xmm0 = xptr(rbp+3*64);	/* Load FFT data */ \
	xmm0 *= xptr(rdx+0*XMM_GMD+rcx); /* mul by grp two-to-minus-phi */ \
	xmm0 += xmm2;		/* x4 = values1 + carry */ \
	split_carry_zpad_word(base2, xmm6, xmm5, xmm4, rax*4); \
no##const1(xmm4 = g->u.xmm.XMM_K_LO);		/* Calc high FFT carry times k */ \
const1(xmm4 = g->u.xmm.XMM_K_TIMES_MULCONST_LO); \
	xmm4 *= xmm6;		/* high_FFT_carry * k_lo */ \
	xmm0 += xmm4;		/* x4 = x4 + high_FFT_carry * k_lo */ \
no##const1(xmm4 = g->u.xmm.XMM_K_HI); \
const1(xmm4 = g->u.xmm.XMM_K_TIMES_MULCONST_HI); \
	xmm4 *= xptr2(g->u.xmm.XMM_LIMIT_INVERSE,rax*4); /* shift k_hi */ \
	xmm6 *= xmm4; \
	rounding(exec, base2, noexec, noexec, xmm0, xmm2, xmm4, rax*4); \
	xmm2 += xmm6;		/* Carry += high_FFT_carry * shifted k_hi */ \
	xmm0 *= xptr(rdx+0*XMM_GMD+XMM_GMD/2+rcx); /* mul by fudged grp two-to-phi */ \
	xptr(rbp+3*64) = xmm0;	/* Save FFT data */ \
\
	rcx = g->BIGLIT_INCR4;/* Different clm values step through */ \
;					/* big/lit array differently */ \
	rax = u8ptr(rdi+rcx+1); /* Load big vs. little flag */ \
	rcx = u8ptr(rdi+rcx);	/* Load fudge factor flag */ \
	rcx &= 0xF0; \
	xmm0 = xptr(rbp+4*64);	/* Load FFT data */ \
	xmm0 *= xptr(rdx+0*XMM_GMD+rcx); /* mul by fudged grp two-to-minus-phi */ \
	xmm0 += xmm2;		/* x5 = values1 + carry */ \
no##const1(xmm4 = g->u.xmm.XMM_K_LO);		/* Calc high FFT carry times k */ \
const1(xmm4 = g->u.xmm.XMM_K_TIMES_MULCONST_LO); \
	xmm4 *= xmm5;		/* high_FFT_carry * k_lo */ \
	xmm0 += xmm4;		/* x5 = x5 + high_FFT_carry * k_lo */ \
no##const1(xmm4 = g->u.xmm.XMM_K_HI); \
const1(xmm4 = g->u.xmm.XMM_K_TIMES_MULCONST_HI); \
	xmm4 *= xptr2(g->u.xmm.XMM_LIMIT_INVERSE,rax*4); /* shift k_hi */ \
	xmm5 *= xmm4; \
	rounding(exec, base2, noexec, noexec, xmm0, xmm2, xmm4, rax*4); \
	xmm2 += xmm5;		/* Carry += high_FFT_carry * shifted k_hi */ \
	xmm0 *= xptr(rdx+0*XMM_GMD+XMM_GMD/2+rcx); /* mul by fudged grp two-to-phi */ \
	xptr(rbp+4*64) = xmm0;	/* Save FFT data */ \
\
	rcx = g->BIGLIT_INCR4;/* Different clm values step through */ \
;					/* big/lit array differently */ \
	rcx = u8ptr(rdi+rcx+2); /* Load fudge factor flag */ \
	rcx &= 0xF0; \
	xmm2 -= XMM_BIGVAL2;/* Remove rounding const from carry */ \
	xmm2 *= xptr(rdx+0*XMM_GMD+XMM_GMD/2+rcx); /* mul by fudged grp two-to-phi */ \
	xmm2 += xptr(rbp+5*64);	/* Load FFT data */ \
	xptr(rbp+5*64) = xmm2;	/* Save FFT data */ \
\
	rax = u8ptr(rdi+1);	/* Load big vs. little flag */ \
	rcx = u8ptr(rdi);	/* Load big vs. little flag */ \
	rcx &= 0xF0; \
	xmm1 = xptr(rsi+1*16);	/* Load carry */ \
	xmm3 = xptr(rsi+3*16);	/* Load hi data carry */ \
	xmm4 = xptr(rbp+0*64+16);	/* Load FFT data */ \
	xmm4 *= xptr(rdx+1*XMM_GMD+rcx); /* mul by fudged grp two-to-minus-phi */ \
	xmm1 += xmm4;		/* x1 = values1 + carry */ \
	split_carry_zpad_word(base2, xmm3, xmm6, xmm4, rax*4+16); \
no##const1(xmm4 = g->u.xmm.XMM_K_LO);		/* Calc high FFT carry times k */ \
const1(xmm4 = g->u.xmm.XMM_K_TIMES_MULCONST_LO); \
	xmm4 *= xmm3;		/* high_FFT_carry * k_lo */ \
	xmm1 += xmm4;		/* x1 = x1 + high_FFT_carry * k_lo */ \
no##const1(xmm5 = g->u.xmm.XMM_K_HI); \
const1(xmm5 = g->u.xmm.XMM_K_TIMES_MULCONST_HI); \
	xmm5 *= xptr2(g->u.xmm.XMM_LIMIT_INVERSE,rax*4+16); /* shift k_hi */ \
	xmm5 *= xmm3; \
	rounding(exec, base2, noexec, noexec, xmm1, xmm3, xmm4, rax*4+16); \
	xmm3 += xmm5;		/* Carry += high_FFT_carry * shifted k_hi */ \
	xmm1 *= xptr(rdx+1*XMM_GMD+XMM_GMD/2+rcx); /* mul by fudged grp two-to-phi */ \
	xptr(rbp+0*64+16) = xmm1;	/* Save FFT data */ \
\
	rax = u8ptr(rdi+3);	/* Load big vs. little flag */ \
	rcx = u8ptr(rdi+2);	/* Load fudge factor flag */ \
	rcx &= 0xF0; \
	xmm1 = xptr(rbp+1*64+16);	/* Load FFT data */ \
	xmm1 *= xptr(rdx+1*XMM_GMD+rcx); /* mul by fudged grp two-to-minus-phi */ \
	xmm1 += xmm3;		/* x2 = values1 + carry */ \
	split_carry_zpad_word(base2, xmm6, xmm5, xmm4, rax*4+16); \
no##const1(xmm4 = g->u.xmm.XMM_K_LO);		/* Calc high FFT carry times k */ \
const1(xmm4 = g->u.xmm.XMM_K_TIMES_MULCONST_LO); \
	xmm4 *= xmm6;		/* high_FFT_carry * k_lo */ \
	xmm1 += xmm4;		/* x2 = x2 + high_FFT_carry * k_lo */ \
no##const1(xmm4 = g->u.xmm.XMM_K_HI); \
const1(xmm4 = g->u.xmm.XMM_K_TIMES_MULCONST_HI); \
	xmm4 *= xptr2(g->u.xmm.XMM_LIMIT_INVERSE,rax*4+16); /* shift k_hi */ \
	xmm6 *= xmm4; \
	rounding(exec, base2, noexec, noexec, xmm1, xmm3, xmm4, rax*4+16); \
	xmm3 += xmm6;		/* Carry += high_FFT_carry * shifted k_hi */ \
	xmm1 *= xptr(rdx+1*XMM_GMD+XMM_GMD/2+rcx); /* mul by fudged grp two-to-phi */ \
	xptr(rbp+1*64+16) = xmm1;	/* Save FFT data */ \
\
	rcx = g->BIGLIT_INCR2;/* Different clm values step through */ \
;					/* big/lit array differently */ \
	rax = u8ptr(rdi+rcx+1); /* Load big vs. little flag */ \
	rcx = u8ptr(rdi+rcx);	/* Load fudge factor flag */ \
	rcx &= 0xF0; \
	xmm1 = xptr(rbp+2*64+16);	/* Load FFT data */ \
	xmm1 *= xptr(rdx+1*XMM_GMD+rcx); /* mul by fudged grp two-to-minus-phi */ \
	xmm1 += xmm3;		/* x3 = values1 + carry */ \
	split_carry_zpad_word(base2, xmm5, xmm6, xmm4, rax*4+16); \
no##const1(xmm4 = g->u.xmm.XMM_K_LO);		/* Calc high FFT carry times k */ \
const1(xmm4 = g->u.xmm.XMM_K_TIMES_MULCONST_LO); \
	xmm4 *= xmm5;		/* high_FFT_carry * k_lo */ \
	xmm1 += xmm4;		/* x3 = x3 + high_FFT_carry * k_lo */ \
no##const1(xmm4 = g->u.xmm.XMM_K_HI); \
const1(xmm4 = g->u.xmm.XMM_K_TIMES_MULCONST_HI); \
	xmm4 *= xptr2(g->u.xmm.XMM_LIMIT_INVERSE,rax*4+16); /* shift k_hi */ \
	xmm5 *= xmm4; \
	rounding(exec, base2, noexec, noexec, xmm1, xmm3, xmm4, rax*4+16); \
	xmm3 += xmm5;		/* Carry += high_FFT_carry * shifted k_hi */ \
	xmm1 *= xptr(rdx+1*XMM_GMD+XMM_GMD/2+rcx); /* mul by fudged grp two-to-phi */ \
	xptr(rbp+2*64+16) = xmm1;	/* Save FFT data */ \
\
	rcx = g->BIGLIT_INCR2;/* Different clm values step through */ \
;					/* big/lit array differently */ \
	rax = u8ptr(rdi+rcx+3); /* Load big vs. little flag */ \
	rcx = u8ptr(rdi+rcx+2); /* Load fudge factor flag */ \
	rcx &= 0xF0; \
	xmm1 = xptr(rbp+3*64+16);	/* Load FFT data */ \
	xmm1 *= xptr(rdx+1*XMM_GMD+rcx); /* mul by fudged grp two-to-minus-phi */ \
	xmm1 += xmm3;		/* x4 = values1 + carry */ \
	split_carry_zpad_word(base2, xmm6, xmm5, xmm4, rax*4+16); \
no##const1(xmm4 = g->u.xmm.XMM_K_LO);		/* Calc high FFT carry times k */ \
const1(xmm4 = g->u.xmm.XMM_K_TIMES_MULCONST_LO); \
	xmm4 *= xmm6;		/* high_FFT_carry * k_lo */ \
	xmm1 += xmm4;		/* x4 = x4 + high_FFT_carry * k_lo */ \
no##const1(xmm4 = g->u.xmm.XMM_K_HI); \
const1(xmm4 = g->u.xmm.XMM_K_TIMES_MULCONST_HI); \
	xmm4 *= xptr2(g->u.xmm.XMM_LIMIT_INVERSE,rax*4+16); /* shift k_hi */ \
	xmm6 *= xmm4; \
	rounding(exec, base2, noexec, noexec, xmm1, xmm3, xmm4, rax*4+16); \
	xmm3 += xmm6;		/* Carry += high_FFT_carry * shifted k_hi */ \
	xmm1 *= xptr(rdx+1*XMM_GMD+XMM_GMD/2+rcx); /* mul by fudged grp two-to-phi */ \
	xptr(rbp+3*64+16) = xmm1;	/* Save FFT data */ \
\
	rcx = g->BIGLIT_INCR4;/* Different clm values step through */ \
;					/* big/lit array differently */ \
	rax = u8ptr(rdi+rcx+1); /* Load big vs. little flag */ \
	rcx = u8ptr(rdi+rcx);	/* Load fudge factor flag */ \
	rcx &= 0xF0; \
	xmm1 = xptr(rbp+4*64+16);	/* Load FFT data */ \
	xmm1 *= xptr(rdx+1*XMM_GMD+rcx); /* mul by fudged grp two-to-minus-phi */ \
	xmm1 += xmm3;		/* x5 = values1 + carry */ \
no##const1(xmm4 = g->u.xmm.XMM_K_LO);		/* Calc high FFT carry times k */ \
const1(xmm4 = g->u.xmm.XMM_K_TIMES_MULCONST_LO); \
	xmm4 *= xmm5;		/* high_FFT_carry * k_lo */ \
	xmm1 += xmm4;		/* x5 = x5 + high_FFT_carry * k_lo */ \
no##const1(xmm4 = g->u.xmm.XMM_K_HI); \
const1(xmm4 = g->u.xmm.XMM_K_TIMES_MULCONST_HI); \
	xmm4 *= xptr2(g->u.xmm.XMM_LIMIT_INVERSE,rax*4+16); /* shift k_hi */ \
	xmm5 *= xmm4; \
	rounding(exec, base2, noexec, noexec, xmm1, xmm3, xmm4, rax*4+16); \
	xmm3 += xmm5;		/* Carry += high_FFT_carry * shifted k_hi */ \
	xmm1 *= xptr(rdx+1*XMM_GMD+XMM_GMD/2+rcx); /* mul by fudged grp two-to-phi */ \
	xptr(rbp+4*64+16) = xmm1;	/* Save FFT data */ \
\
	rcx = g->BIGLIT_INCR4;/* Different clm values step through */ \
;					/* big/lit array differently */ \
	rcx = u8ptr(rdi+rcx+2); /* Load big vs. little flag */ \
	rcx &= 0xF0; \
	xmm3 -= XMM_BIGVAL2;/* Remove rounding const from carry */ \
	xmm3 *= xptr(rdx+1*XMM_GMD+XMM_GMD/2+rcx); /* mul by fudged grp two-to-phi */ \
	xmm3 += xptr(rbp+5*64+16);	/* Add FFT data */ \
	xptr(rbp+5*64+16) = xmm3;	/* Save FFT data */ \
\
	xmm4 = XMM_BIGVAL2; \
	xptr(rsi+0*16) = xmm4;	/* Clear carry */ \
	xptr(rsi+1*16) = xmm4; \
	xptr(rsi+2*16) = xmm4; \
	xptr(rsi+3*16) = xmm4; \
\
\
\
/* *************** 1D normalized add/sub macro ****************** */ \
/* This macro adds or subtracts, then "normalizes" eight FFT */ \
/* data values.  This involves multiplying the summed values by */ \
/* two-to-minus-phi.  Rounding the value to an integer.  Making sure */ \
/* the integer is smaller than the maximum allowable integer, generating */ \
/* a carry if necessary. Finally, the value is multiplied by two-to-phi */ \
/* and stored. */ \
/* xmm3 = carry #2 */ \
/* xmm2 = carry #1 */ \
/* rcx = pointer to the first number */ \
/* rdx = pointer to the second number */ \
/* rsi = pointer to destination */ \
/* rbp = pointer two-to-phi multipliers */ \
/* rdi = pointer to array of big vs. little flags */ \
/* rbx = big vs. little word flag #2 */ \
/* rax = big vs. little word flag #1 */ \
/* A pipelined version of this code: */ \
/*	movzx	rax, BYTE PTR [rdi]	;; Load big vs. little flags */ \
/*	xload	xmm0, [rdx+0*dist1]	;; Load second number */ \
/*	fop	xmm0, [rcx]		;; Add/sub first number */ \
/*	mulpd	xmm0, [rbp+0]		;; Mul values1 by two-to-minus-phi */ \
/*	addpd	xmm0, xmm4		;; x = values + carry */ \
/*	xload	xmm2, XMM_LIMIT_BIGMAX[rax];; Load maximum * BIGVAL - BIGVAL */ \
/*	addpd	xmm2, xmm0		;; y = top bits of x */ \
/*	xload	xmm6, XMM_LIMIT_BIGMAX_NEG[rax];; Load -(maximum*BIGVAL-BIGVAL) */ \
/*	addpd	xmm6, xmm2		;; z = y - (maximum * BIGVAL - BIGVAL) */ \
/*	subpd	xmm0, xmm6		;; rounded value = x - z */ \
/*	mulpd	xmm2, XMM_LIMIT_INVERSE[rax];; next carry = shifted y */ \
/*	mulpd	xmm0, [rbp+16]		;; new value = val * two-to-phi */ \
/*	xstore	[rsi+0*dist1], xmm0	;; Save new value */ \
\

#define xnorm_op_1d(fop, ttp, base2) \
ttp(rax = u8ptr(rdi));	/* Load big vs. little flags */ \
ttp(rbx = u8ptr(rdi+2)); \
	xmm0 = xptr(rdx);		/* Load second number */ \
	fop(xmm0, xptr(rcx));		/* Add/sub first number */ \
ttp(xmm0 *= xptr(rbp));		/* Mul values1 by two-to-minus-phi */ \
ttp(xmm0 *= g->u.xmm.XMM_NORM012_FF);	/* Mul by FFTLEN/2 */ \
	xmm1 = xptr(rdx+32);		/* Load second number */ \
	fop(xmm1, xptr(rcx+32));		/* Add/sub first number */ \
ttp(xmm1 *= xptr(rbp+64));		/* Mul values2 by two-to-minus-phi */ \
ttp(xmm1 *= g->u.xmm.XMM_NORM012_FF);	/* Mul by FFTLEN/2 */ \
	xmm0 += xmm2;		/* x1 = values + carry */ \
	xmm1 += xmm3;		/* x2 = values + carry */ \
	rounding_interleaved(ttp, base2, noexec, noexec, xmm0, xmm2, xmm4, rax, xmm1, xmm3, xmm5, rbx);  \
ttp(xmm0 *= xptr(rbp+16));		/* new value1 = val * two-to-phi */ \
ttp(xmm1 *= xptr(rbp+80));		/* new value2 = val * two-to-phi */ \
	xptr(rsi) = xmm0;		/* Save value1 */ \
	xptr(rsi+32) = xmm1;		/* Save value2 */ \
ttp(rax = u8ptr(rdi+1));	/* Load big vs. little flags */ \
ttp(rbx = u8ptr(rdi+3)); \
	xmm0 = xptr(rdx+16);		/* Load values1 */ \
	fop(xmm0, xptr(rcx+16));		/* Add/sub first number */ \
ttp(xmm0 *= xptr(rbp+32));		/* Mul values1 by two-to-minus-phi */ \
ttp(xmm0 *= g->u.xmm.XMM_NORM012_FF);	/* Mul by FFTLEN/2 */ \
	xmm1 = xptr(rdx+48);		/* Load values2 */ \
	fop(xmm1, xptr(rcx+48));		/* Add/sub first number */ \
ttp(xmm1 *= xptr(rbp+96));		/* Mul values2 by two-to-minus-phi */ \
ttp(xmm1 *= g->u.xmm.XMM_NORM012_FF);	/* Mul by FFTLEN/2 */ \
	xmm0 += xmm2;		/* x1 = values + carry */ \
	xmm1 += xmm3;		/* x2 = values + carry */ \
	rounding_interleaved(ttp, base2, noexec, noexec, xmm0, xmm2, xmm4, rax, xmm1, xmm3, xmm5, rbx);  \
ttp(xmm0 *= xptr(rbp+48));		/* new value1 = val * two-to-phi */ \
ttp(xmm1 *= xptr(rbp+112));		/* new value2 = val * two-to-phi */ \
	xptr(rsi+16) = xmm0;		/* Save new value1 */ \
	xptr(rsi+48) = xmm1;		/* Save new value2 */ \
ttp(rdi += 4);			/* Next flags ptr */ \
	rcx += 64;			/* Next src ptr */ \
	rdx += 64;			/* Next src ptr */ \
ttp(rbp += 128);		/* Next two-to-phi ptr */ \
	rsi += 64;			/* Next dest ptr */ \
\
\
\
/* This macro finishes the normalize process by adding the final */ \
/* carry from the first pass back into the lower two data values. */ \
/* xmm2,xmm3 = carries */ \
/* rax = pointer to the FFT data values */ \
/* rbx = pointer two-to-phi multipliers */ \
\

#define xnorm_op_1d_mid_cleanup \
	xmm0[0] = f64ptr(rax+8);		/* Load values1 */ \
	xmm1[0] = f64ptr(rax+40);	/* Load values2 */ \
	xmm2[0] -= XMM_BIGVAL1;/* Remove BIGVAL from carries */ \
	xmm3[0] -= XMM_BIGVAL1; \
	xmm2[0] *= f64ptr(rbx+24);	/* carry1 *= two-to-phi */ \
	xmm3[0] *= f64ptr(rbx+88);	/* carry2 *= two-to-phi */ \
	xmm0[0] += xmm2[0];		/* x1 = values + carry */ \
	xmm1[0] += xmm3[0];		/* x2 = values + carry */ \
	f64ptr(rax+8) = xmm0[0];		/* Save new value1 */ \
	f64ptr(rax+40) = xmm1[0];	/* Save new value2 */ \
	xmm2[0] = xmm2[1]; xmm2[1] = XMM_BIGVAL1; /* Rotate carry */ \
	xmm3[0] = xmm3[1]; xmm3[1] = XMM_BIGVAL1; /* Rotate carry */ \
\
\
/* *************** 1D normalized add/sub macro ****************** */ \
/* This macro adds and subtracts, then "normalizes" eight FFT */ \
/* data values.  This involves multiplying the sum values by */ \
/* two-to-minus-phi.  Adding, subtracting and rounding the value to an */ \
/* integer.  Make sure the integer is smaller than the maximum allowable */ \
/* integer, generating carries if necessary.  Finally, the values are */ \
/* multiplied by two-to-phi and stored. */ \
/* xmm7 = sub carry #2 */ \
/* xmm6 = sub carry #1 */ \
/* xmm3 = add carry #2 */ \
/* xmm2 = add carry #1 */ \
/* rcx = pointer to the first number */ \
/* rdx = pointer to the second number */ \
/* rsi = pointer to destination #1 */ \
/* rbp = pointer to destination #2 */ \
/* rbx = pointer two-to-phi multipliers */ \
/* rdi = pointer to array of big vs. little flags */ \
/* rax = big vs. little word flag #1 */ \
\

#define xnorm_addsub_1d(ttp, base2) \
ttp(rax = u8ptr(rdi));	/* Load big vs. little flags */ \
	xmm0 = xptr(rcx);		/* Load first number */ \
	xmm0 += xptr(rdx);		/* Add second number */ \
ttp(xmm5 = xptr(rbx));		/* Load fudged two-to-minus-phi */ \
ttp(xmm5 *= g->u.xmm.XMM_NORM012_FF);	/* Mul by FFTLEN/2 */ \
ttp(xmm0 *= xmm5);		/* Mul values1 by two-to-minus-phi */ \
	xmm1 = xptr(rcx);		/* Load first number */ \
	xmm1 -= xptr(rdx);		/* Sub second number */ \
ttp(xmm1 *= xmm5);		/* Mul values2 by two-to-minus-phi */ \
	xmm0 += xmm2;		/* x1 = values + carry */ \
	xmm1 += xmm6;		/* x2 = values + carry */ \
	rounding_interleaved(ttp, base2, noexec, noexec, xmm0, xmm2, xmm4, rax, xmm1, xmm6, xmm5, rax); \
ttp(xmm0 *= xptr(rbx+16));		/* new value1 = val * two-to-phi */ \
ttp(xmm1 *= xptr(rbx+16));		/* new value2 = val * two-to-phi */ \
	xptr(rsi) = xmm0;		/* Save value1 */ \
	xptr(rbp) = xmm1;		/* Save value2 */ \
\
ttp(rax = u8ptr(rdi+2));	/* Load big vs. little flags */ \
	xmm0 = xptr(rcx+32);		/* Load first number */ \
	xmm0 += xptr(rdx+32);		/* Add second number */ \
ttp(xmm5 = xptr(rbx+64));		/* Load fudged two-to-minus-phi */ \
ttp(xmm5 *= g->u.xmm.XMM_NORM012_FF);	/* Mul by FFTLEN/2 */ \
ttp(xmm0 *= xmm5);		/* Mul values1 by two-to-minus-phi */ \
	xmm1 = xptr(rcx+32);		/* Load first number */ \
	xmm1 -= xptr(rdx+32);		/* Sub second number */ \
ttp(xmm1 *= xmm5);		/* Mul values2 by two-to-minus-phi */ \
	xmm0 += xmm3;		/* x1 = values + carry */ \
	xmm1 += xmm7;		/* x2 = values + carry */ \
	rounding_interleaved(ttp, base2, noexec, noexec, xmm0, xmm3, xmm4, rax, xmm1, xmm7, xmm5, rax); \
ttp(xmm0 *= xptr(rbx+80));		/* new value1 = val * two-to-phi */ \
ttp(xmm1 *= xptr(rbx+80));		/* new value2 = val * two-to-phi */ \
	xptr(rsi+32) = xmm0;		/* Save value1 */ \
	xptr(rbp+32) = xmm1;		/* Save value2 */ \
\
ttp(rax = u8ptr(rdi+1));	/* Load big vs. little flags */ \
	xmm0 = xptr(rcx+16);		/* Load first number */ \
	xmm0 += xptr(rdx+16);		/* Add second number */ \
ttp(xmm5 = xptr(rbx+32));		/* Load fudged two-to-minus-phi */ \
ttp(xmm5 *= g->u.xmm.XMM_NORM012_FF);	/* Mul by FFTLEN/2 */ \
ttp(xmm0 *= xmm5);		/* Mul values1 by two-to-minus-phi */ \
	xmm1 = xptr(rcx+16);		/* Load first number */ \
	xmm1 -= xptr(rdx+16);		/* Sub second number */ \
ttp(xmm1 *= xmm5);		/* Mul values2 by two-to-minus-phi */ \
	xmm0 += xmm2;		/* x1 = values + carry */ \
	xmm1 += xmm6;		/* x2 = values + carry */ \
	rounding_interleaved(ttp, base2, noexec, noexec, xmm0, xmm2, xmm4, rax, xmm1, xmm6, xmm5, rax); \
ttp(xmm0 *= xptr(rbx+48));		/* new value1 = val * two-to-phi */ \
ttp(xmm1 *= xptr(rbx+48));		/* new value2 = val * two-to-phi */ \
	xptr(rsi+16) = xmm0;		/* Save value1 */ \
	xptr(rbp+16) = xmm1;		/* Save value2 */ \
\
ttp(rax = u8ptr(rdi+3));	/* Load big vs. little flags */ \
	xmm0 = xptr(rcx+48);		/* Load first number */ \
	xmm0 += xptr(rdx+48);		/* Add second number */ \
ttp(xmm5 = xptr(rbx+96));		/* Load fudged two-to-minus-phi */ \
ttp(xmm5 *= g->u.xmm.XMM_NORM012_FF);	/* Mul by FFTLEN/2 */ \
ttp(xmm0 *= xmm5);		/* Mul values1 by two-to-minus-phi */ \
	xmm1 = xptr(rcx+48);		/* Load first number */ \
	xmm1 -= xptr(rdx+48);		/* Sub second number */ \
ttp(xmm1 *= xmm5);		/* Mul values2 by two-to-minus-phi */ \
	xmm0 += xmm3;		/* x1 = values + carry */ \
	xmm1 += xmm7;		/* x2 = values + carry */ \
	rounding_interleaved(ttp, base2, noexec, noexec, xmm0, xmm3, xmm4, rax, xmm1, xmm7, xmm5, rax); \
ttp(xmm0 *= xptr(rbx+112));		/* new value1 = val * two-to-phi */ \
ttp(xmm1 *= xptr(rbx+112));		/* new value2 = val * two-to-phi */ \
	xptr(rsi+48) = xmm0;		/* Save value1 */ \
	xptr(rbp+48) = xmm1;		/* Save value2 */ \
\
ttp(rdi += 4);			/* Next flags ptr */ \
	rcx += 64;			/* Next src ptr */ \
	rdx += 64;			/* Next src ptr */ \
ttp(rbx += 128);		/* Next two-to-phi ptr */ \
	rsi += 64;			/* Next dest ptr */ \
	rbp += 64;			/* Next dest ptr */ \
\
\
/* This macro finishes the normalize process by adding the final */ \
/* carry from the first pass back into the lower two data values. */ \
/* xmm2,xmm3 = carries #1 */ \
/* xmm6,xmm7 = carries #2 */ \
/* rax = pointer to the FFT data values #1 */ \
/* top of stack = pointer to the FFT destination #1 and #2 */ \
/* rbx = pointer two-to-phi multipliers */ \
\

#define xnorm_addsub_1d_mid_cleanup(dest1, dest2) \
	rax = dest1;		/* Restore dest #1 pointer */ \
	xmm0[0] = f64ptr(rax+8);		/* Load values1 */ \
	xmm1[0] = f64ptr(rax+40);	/* Load values2 */ \
	xmm2[0] -= XMM_BIGVAL1;/* Remove BIGVAL from carries */ \
	xmm3[0] -= XMM_BIGVAL1; \
	xmm2[0] *= f64ptr(rbx+24);	/* carry1 *= two-to-phi */ \
	xmm3[0] *= f64ptr(rbx+88);	/* carry2 *= two-to-phi */ \
	xmm0[0] += xmm2[0];		/* x1 = values + carry */ \
	xmm1[0] += xmm3[0];		/* x2 = values + carry */ \
	f64ptr(rax+8) = xmm0[0];		/* Save new value1 */ \
	f64ptr(rax+40) = xmm1[0];	/* Save new value2 */ \
	xmm2[0] = xmm2[1]; xmm2[1] = XMM_BIGVAL1; /* Rotate carry */ \
	xmm3[0] = xmm3[1]; xmm3[1] = XMM_BIGVAL1; /* Rotate carry */ \
\
	rax = dest2;		/* Get FFT data pointer #2 */ \
	xmm0[0] = f64ptr(rax+8);		/* Load values1 */ \
	xmm1[0] = f64ptr(rax+40);	/* Load values2 */ \
	xmm6[0] -= XMM_BIGVAL1;/* Remove BIGVAL from carries */ \
	xmm7[0] -= XMM_BIGVAL1; \
	xmm6[0] *= f64ptr(rbx+24);	/* carry1 *= two-to-phi */ \
	xmm7[0] *= f64ptr(rbx+88);	/* carry2 *= two-to-phi */ \
	xmm0[0] += xmm6[0];		/* x1 = values + carry */ \
	xmm1[0] += xmm7[0];		/* x2 = values + carry */ \
	f64ptr(rax+8) = xmm0[0];		/* Save new value1 */ \
	f64ptr(rax+40) = xmm1[0];	/* Save new value2 */ \
	xmm6[0] = xmm6[1]; xmm6[1] = XMM_BIGVAL1; /* Rotate carry */ \
	xmm7[0] = xmm7[1]; xmm7[1] = XMM_BIGVAL1; /* Rotate carry */ \
	\
\
/* rsi = pointer to the FFT data values #1 */ \
/* rbp = pointer to the FFT data values #2 */ \
\

#define xnorm_addsub_1d_cleanup \
	xnorm_top_carry_cmn rsi, xmm3, 0;/* Adjust top carry if necessary */ \
	xnorm_top_carry_cmn rbp, xmm7, 0;/* Adjust top carry if necessary */ \
\
	xmm0[0] = f64ptr(rsi);		/* Load values1 */ \
	xmm1[0] = f64ptr(rsi+32);	/* Load values2 */ \
	xmm3[0] -= XMM_BIGVAL1;/* Remove BIGVAL from carries */ \
	xmm2[0] -= XMM_BIGVAL1; \
	xmm3[0] *= g->u.xmm.XMM_MINUS_C[0];	/* Mul wrap-around carry by -c */ \
	xmm3[0] *= f64ptr(rbx+16);	/* carry *= two-to-phi */ \
	xmm2[0] *= f64ptr(rbx+80);	/* carry *= two-to-phi */ \
	xmm0[0] += xmm3[0];		/* value1 = values + carry */ \
	xmm1[0] += xmm2[0];		/* value2 = values + carry */ \
	f64ptr(rsi) = xmm0[0];		/* Save new value1 */ \
	f64ptr(rsi+32) = xmm1[0];	/* Save new value2 */ \
\
	xmm0[0] = f64ptr(rbp);		/* Load values1 */ \
	xmm1[0] = f64ptr(rbp+32);	/* Load values2 */ \
	xmm7[0] -= XMM_BIGVAL1;/* Remove BIGVAL from carries */ \
	xmm6[0] -= XMM_BIGVAL1; \
	xmm7[0] *= g->u.xmm.XMM_MINUS_C[0];	/* Mul wrap-around carry by -c */ \
	xmm7[0] *= f64ptr(rbx+16);	/* carry *= two-to-phi */ \
	xmm6[0] *= f64ptr(rbx+80);	/* carry *= two-to-phi */ \
	xmm0[0] += xmm7[0];		/* value1 = values + carry */ \
	xmm1[0] += xmm6[0];		/* value2 = values + carry */ \
	f64ptr(rbp) = xmm0[0];		/* Save new value1 */ \
	f64ptr(rbp+32) = xmm1[0];	/* Save new value2 */ \
\
\
/* *************** 1D normalized smallmul macro ****************** */ \
/* This macro multiplies by a small constant, then "normalizes" eight FFT */ \
/* data values. */ \
/* xmm7 = small multiplier value */ \
/* xmm3 = carry #2 */ \
/* xmm2 = carry #1 */ \
/* rsi = pointer to destination */ \
/* rbp = pointer two-to-phi multipliers */ \
/* rdi = pointer to array of big vs. little flags */ \
/* rax = big vs. little word flag #1 */ \
/* rcx = big vs. little word flag #2 */ \
\

#define xnorm_smallmul_1d(ttp, base2) { \
ttp(uintptr_t ax = u8ptr(rdi));	/* Load big vs. little flags */ \
ttp(uintptr_t cx = u8ptr(rdi+2)); \
	xmm0 = xptr(rsi);		/* Load values1 */ \
	xmm0 *= xmm7;		/* Mul by small value * FFTLEN/2 */ \
ttp(xmm0 *= xptr(rbp));		/* Mul values1 by two-to-minus-phi */ \
	xmm1 = xptr(rsi+32);		/* Load values2 */ \
	xmm1 *= xmm7;		/* Mul by small value * FFTLEN/2 */ \
ttp(xmm1 *= xptr(rbp+64));		/* Mul values2 by two-to-minus-phi */ \
	xmm0 += xmm2;		/* x1 = values + carry */ \
	xmm1 += xmm3;		/* x2 = values + carry */ \
	rounding_interleaved(ttp, base2, noexec, noexec, xmm0, xmm2, xmm4, ax, xmm1, xmm3, xmm5, cx); \
ttp(xmm0 *= xptr(rbp+16));		/* new value1 = val * two-to-phi */ \
ttp(xmm1 *= xptr(rbp+80));		/* new value2 = val * two-to-phi */ \
	xptr(rsi) = xmm0;		/* Save value1 */ \
	xptr(rsi+32) = xmm1;		/* Save value2 */ \
ttp(ax = u8ptr(rdi+1));	/* Load big vs. little flags */ \
ttp(cx = u8ptr(rdi+3)); \
	xmm0 = xptr(rsi+16);		/* Load values1 */ \
	xmm0 *= xmm7;		/* Mul by small value * FFTLEN/2 */ \
ttp(xmm0 *= xptr(rbp+32));		/* Mul values1 by two-to-minus-phi */ \
	xmm1 = xptr(rsi+48);		/* Load values2 */ \
	xmm1 *= xmm7;		/* Mul by small value * FFTLEN/2 */ \
ttp(xmm1 *= xptr(rbp+96));		/* Mul values2 by two-to-minus-phi */ \
	xmm0 += xmm2;		/* x1 = values + carry */ \
	xmm1 += xmm3;		/* x2 = values + carry */ \
	rounding_interleaved(ttp, base2, noexec, noexec, xmm0, xmm2, xmm4, ax, xmm1, xmm3, xmm5, cx); \
ttp(xmm0 *= xptr(rbp+48));		/* new value1 = val * two-to-phi */ \
ttp(xmm1 *= xptr(rbp+112));		/* new value2 = val * two-to-phi */ \
	xptr(rsi+16) = xmm0;		/* Save new value1 */ \
	xptr(rsi+48) = xmm1;		/* Save new value2 */ \
ttp(rdi += 4);			/* Next flags ptr */ \
ttp(rbp += 128);		/* Next two-to-phi ptr */ \
	rsi += 64;			/* Next dest ptr */ \
}

/* This macro finishes the smallmul normalize process by adding the final */ \
/* carry from the first pass back into the lower two data values. */ \
/* xmm2,xmm3 = carries */ \
/* rsi = pointer to the FFT data values */ \
/* rbp = pointer two-to-phi multipliers */ \
/* rdi = pointer to biglit array */ \
\

#define xnorm_smallmul_1d_mid_cleanup(base2) \
	xnorm012_1d_mid(exec, base2); \
\
\

#define xnorm_smallmul_1d_cleanup(base2) \
\
	if(g->ZERO_PADDED_FFT == 0) {	/* Zero-padded FFT? */ \
		xnorm_top_carry_1d;		/* No, do a very standard carry */ \
		xnorm012_1d(base2); \
	}else{			/* Yes, do special zpad carry */ \
		xnorm_smallmul_1d_zpad_cleanup(base2);/* Do the special zpad carry */ \
	} \
\
\
/* This macro is similar to xnorm012_1d_zpad for handling zpad carries */ \
/* rsi = pointer to the FFT data values */ \
/* rbp = pointer two-to-phi multipliers */ \
/* rdi = big vs. little array pointer */ \
\

#define xnorm_smallmul_1d_zpad_cleanup(base2) { \
	double f0, f1, f2, f4; \
	uintptr_t rax; \
	f2 = xmm2[0]; \
;	/* Copy and integerize data from 7 words above halfway point to ZPAD0-ZPAD6 */ \
;	/* Clear words 5,6,7 */ \
;	/* Then we can make an exact copy of most of the xnorm012_1d_zpad code */ \
\
	rax = u8ptr(rdi);	/* Load big vs. little flags */ \
	f0 = f64ptr(rsi+32);	/* Value1 */ \
	f2 -= XMM_BIGVAL1;/* Remove XMM_BIGVAL from the carry */ \
	single_split_lower_zpad_word(base2, f0, f2, f4, rax); \
	g->ZPAD0_6[0] = f0; \
\
	rax = u8ptr(rdi+1);	/* Load big vs. little flags */ \
	f0 = f64ptr(rsi+48);	/* Value2 */ \
	f0 *= f64ptr(rbp+32);	/* Mul values2 by two-to-minus-phi */ \
	f0 *= g->u.xmm.XMM_NORM012_FF[0];	/* Mul by FFTLEN/2 */ \
	single_split_lower_zpad_word(base2, f0, f2, f4, rax); \
	g->ZPAD0_6[1] = f0; \
\
	rax = u8ptr(rdi+4);	/* Load big vs. little flags */ \
	f0 = f64ptr(rsi+64+32);	/* Value3 */ \
	f0 *= f64ptr(rbp+128);	/* Mul values3 by two-to-minus-phi */ \
	f0 *= g->u.xmm.XMM_NORM012_FF[0];	/* Mul by FFTLEN/2 */ \
	single_split_lower_zpad_word(base2, f0, f2, f4, rax); \
	g->ZPAD0_6[2] = f0; \
\
	f0 = f64ptr(rsi+64+48);	/* Value4 */ \
	f0 *= f64ptr(rbp+128+32);	/* Mul values4 by two-to-minus-phi */ \
	f0 *= g->u.xmm.XMM_NORM012_FF[0];	/* Mul by FFTLEN/2 */ \
	f0 = roundf64(f0); /* Round to an integer */ \
	f0 += f2;		/* Value4 + carry */ \
	g->ZPAD0_6[3] = f0; \
\
	rsi += g->u.xmm.ZPAD_WORD5_OFFSET; \
	rbp += g->u.xmm.ZPAD_WORD5_RBP_OFFSET; \
	f0 = f64ptr(rsi+32);	/* Value5 */ \
	f0 *= f64ptr(rbp);		/* Mul values5 by two-to-minus-phi */ \
	f0 *= g->u.xmm.XMM_NORM012_FF[0];	/* Mul by FFTLEN/2 */ \
	f0 = roundf64(f0); /* Round to an integer */ \
	g->ZPAD0_6[4] = f0; \
	f0 = f64ptr(rsi+48);	/* Value6 */ \
	f0 *= f64ptr(rbp+32);	/* Mul values6 by two-to-minus-phi */ \
	f0 *= g->u.xmm.XMM_NORM012_FF[0];	/* Mul by FFTLEN/2 */ \
	f0 = roundf64(f0); /* Round to an integer */ \
	g->ZPAD0_6[5] = f0; \
	f0 = f64ptr(rsi+64+32);	/* Value7 */ \
	f0 *= f64ptr(rbp+128);	/* Mul values7 by two-to-minus-phi */ \
	f0 *= g->u.xmm.XMM_NORM012_FF[0];	/* Mul by FFTLEN/2 */ \
	f0 = roundf64(f0); /* Round to an integer */ \
	g->ZPAD0_6[6] = f0; \
	f0 = 0;		/* Clear highest words */ \
	f64ptr(rsi+32) = f0; \
	f64ptr(rsi+48) = f0; \
	f64ptr(rsi+64+32) = f0; \
	rsi -= g->u.xmm.ZPAD_WORD5_OFFSET; \
	rbp -= g->u.xmm.ZPAD_WORD5_RBP_OFFSET; \
\
;	/* Divide the zpad data by k.  Store the integer part in XMM_TMP */ \
;	/* and the remainder in ZPAD0.  Later we will wrap the integer part */ \
;	/* down to the bottom of the FFT data area (and multiply by -c). */ \
;	/* And we will store the remainder in the upper half of the FFT */ \
;	/* data area. */ \
\
;	/* Note there are three cases to handle.  K is smaller than a big word. */ \
;	/* K is between one and 2 big words in size.  And K is more than */ \
;	/* 2 big words in size. */ \
\
	xfunc572(g); \
\
;	/* Now normalize the data above the halfway point.  Remember that the */ \
;	/* column two-to-phi multiplier for the first value will be 1.0. */ \
\
	rax = u8ptr(rdi);	/* First word  */ \
	f0 = g->ZPAD0_6[0];/* Load remainder of divide by k */ \
	f0 += XMM_BIGVAL1; \
	single_rounding(base2, f0, f2, f4, rax); \
	f64ptr(rsi+32) = f0;	/* Save value1 */ \
\
	rax = u8ptr(rdi+1);	/* Load big vs. little flags */ \
	single_rounding(base2, f2, f0, f4, rax); \
	f2 *= f64ptr(rbp+48);	/* new value2 = val * two-to-phi */ \
	f64ptr(rsi+48) = f2;	/* Save value2 */ \
\
	rax = u8ptr(rdi+4);	/* Load big vs. little flags */ \
	single_rounding(base2, f0, f2, f4, rax); \
	f0 *= f64ptr(rbp+128+16);	/* new value3 = val * two-to-phi */ \
	f64ptr(rsi+64+32) = f0;	/* Save value3 */ \
\
	f2 -= XMM_BIGVAL1;/* Remove integer rounding constant */ \
	f2 *= f64ptr(rbp+128+48);	/* value4 = carry * two-to-phi */ \
	f64ptr(rsi+64+48) = f2;	/* Save new value4 */ \
\
;	/* Mul the integer part of (ZPAD data divided by k) by -c in */ \
;	/* preparation for adding it into the lower FFT data area. */ \
;	/* Also add in the shifted high FFT carry at this time. */ \
\
;	/* Now add in and normalize the bottom FFT data.  Remember that the */ \
;	/* column two-to-phi multiplier for the first value will be 1.0.  We  */ \
;	/* must go 6 words deep in case k is 48-50 bits and c is 32 bits. */ \
\
	vec2f64 TMP12 = g->u.xmm.XMM_TMP1; \
	vec2f64 TMP34 = g->u.xmm.XMM_TMP2; \
	vec2f64 TMP56 = g->u.xmm.XMM_TMP3; \
	TMP12 *= g->u.xmm.XMM_MINUS_C;	/* Mul by -c */ \
	TMP34 *= g->u.xmm.XMM_MINUS_C;	/* Mul by -c */ \
	TMP56 *= g->u.xmm.XMM_MINUS_C;	/* Mul by -c */ \
	rax = u8ptr(rdi);	/* First word  */ \
	f0 = TMP12[0];	/* Load integer part of divide by k */ \
	f0 += XMM_BIGVAL1; \
	f0 += f64ptr(rsi);		/* Add in the FFT data */ \
	single_rounding(base2, f0, f2, f4, rax); \
	f64ptr(rsi) = f0;		/* Save value1 */ \
\
	rax = u8ptr(rdi+1);	/* Load big vs. little flags */ \
	f0 = TMP12[1];	/* Load integer part of divide by k */ \
	f1 = f64ptr(rsi+16);	/* Load FFT data */ \
	f1 *= f64ptr(rbp+32);	/* Mul values2 by two-to-minus-phi */ \
	f1 *= g->u.xmm.XMM_NORM012_FF[0];	/* Mul by FFTLEN/2 */ \
	f0 += f1;		/* Add in the FFT data */ \
	f0 += f2;		/* x2 = value + carry */ \
	single_rounding(base2, f0, f2, f4, rax); \
	f0 *= f64ptr(rbp+48);	/* new value2 = val * two-to-phi */ \
	f64ptr(rsi+16) = f0;	/* Save value2 */ \
\
	rax = u8ptr(rdi+4);	/* Load big vs. little flags */ \
	f0 = TMP34[0];	/* Load integer part of divide by k */ \
	f1 = f64ptr(rsi+64);	/* Load FFT data */ \
	f1 *= f64ptr(rbp+128);	/* Mul values3 by two-to-minus-phi */ \
	f1 *= g->u.xmm.XMM_NORM012_FF[0];	/* Mul by FFTLEN/2 */ \
	f0 += f1;		/* Add in the FFT data */ \
	f0 += f2;		/* x3 = value + carry */ \
	single_rounding(base2, f0, f2, f4, rax); \
	f0 *= f64ptr(rbp+128+16);	/* new value3 = val * two-to-phi */ \
	f64ptr(rsi+64) = f0;	/* Save value3 */ \
\
	rax = u8ptr(rdi+5);	/* Load big vs. little flags */ \
	f0 = TMP34[1];	/* Load integer part of divide by k */ \
	f1 = f64ptr(rsi+64+16);	/* Load FFT data */ \
	f1 *= f64ptr(rbp+128+32);	/* Mul values4 by two-to-minus-phi */ \
	f1 *= g->u.xmm.XMM_NORM012_FF[0];	/* Mul by FFTLEN/2 */ \
	f0 += f1;		/* Add in the FFT data */ \
	f0 += f2;		/* x4 = value + carry */ \
	single_rounding(base2, f0, f2, f4, rax); \
	f0 *= f64ptr(rbp+128+48);	/* new value4 = val * two-to-phi */ \
	f64ptr(rsi+64+16) = f0;	/* Save value4 */ \
\
	if(g->FFTLEN != 80 && g->FFTLEN != 112) {	/* Length 80 and 112 have different */ \
		/* memory addresses for the fourth */ \
		/* and higher data elements */ \
	\
		rax = u8ptr(rdi+8);	/* Load big vs. little flags */ \
		f0 = TMP56[0];	/* Load integer part of divide by k */ \
		f1 = f64ptr(rsi+128);	/* Load FFT data */ \
		f1 *= f64ptr(rbp+256);	/* Mul values4 by two-to-minus-phi */ \
		f1 *= g->u.xmm.XMM_NORM012_FF[0];	/* Mul by FFTLEN/2 */ \
		f0 += f1;		/* Add in the FFT data */ \
		f0 += f2;		/* x5 = value + carry */ \
		single_rounding(base2, f0, f2, f4, rax); \
		f0 *= f64ptr(rbp+256+16);	/* new value5 = val * two-to-phi */ \
		f64ptr(rsi+128) = f0;	/* Save value5 */ \
	\
		rax = u8ptr(rdi+9);	/* Load big vs. little flags */ \
		f0 = TMP56[1];	/* Load integer part of divide by k */ \
		f1 = f64ptr(rsi+128+16);	/* Load FFT data */ \
		f1 *= f64ptr(rbp+256+32);	/* Mul values6 by two-to-minus-phi */ \
		f1 *= g->u.xmm.XMM_NORM012_FF[0];	/* Mul by FFTLEN/2 */ \
		f0 += f1;		/* Add in the FFT data */ \
		f0 += f2;		/* x6 = value + carry */ \
		single_rounding(base2, f0, f2, f4, rax); \
		f0 *= f64ptr(rbp+256+48);	/* new value6 = val * two-to-phi */ \
		f64ptr(rsi+128+16) = f0;	/* Save value6 */ \
	\
		f2 -= XMM_BIGVAL1;/* Remove rounding constant */ \
		f2 *= f64ptr(rbp+384+16);	/* new value7 = val * two-to-phi */ \
		f2 += f64ptr(rsi+192);	/* Add in FFT data */ \
		f64ptr(rsi+192) = f2;	/* Save value7 */ \
	}else{ \
\
;	/* Same as the above but with different addresses required by */ \
;	/* the length 80 and 112 FFT lengths */ \
\
		rax = u8ptr(rdi);	/* Load big vs. little flags */ \
		f0 = TMP56[0];	/* Load integer part of divide by k */ \
		f1 = f64ptr(rsi+8);		/* Load FFT data */ \
		f1 *= f64ptr(rbp+8);		/* Mul values5 by two-to-minus-phi */ \
		f1 *= g->u.xmm.XMM_NORM012_FF[0];	/* Mul by FFTLEN/2 */ \
		f0 += f1;		/* Add in the FFT data */ \
		f0 += f2;		/* x5 = value + carry */ \
		single_rounding(base2, f0, f2, f4, rax+8); \
		f0 *= f64ptr(rbp+16+8);	/* new value5 = val * two-to-phi */ \
		f64ptr(rsi+8) = f0;		/* Save value5 */ \
	\
		rax = u8ptr(rdi+1);	/* Load big vs. little flags */ \
		f0 = TMP56[1];	/* Load integer part of divide by k */ \
		f1 = f64ptr(rsi+16+8);	/* Load FFT data */ \
		f1 *= f64ptr(rbp+32+8);	/* Mul values6 by two-to-minus-phi */ \
		f1 *= g->u.xmm.XMM_NORM012_FF[0];	/* Mul by FFTLEN/2 */ \
		f0 += f1;		/* Add in the FFT data */ \
		f0 += f2;		/* x6 = value + carry */ \
		single_rounding(base2, f0, f2, f4, rax+8); \
		f0 *= f64ptr(rbp+48+8);	/* new value6 = val * two-to-phi */ \
		f64ptr(rsi+16+8) = f0;	/* Save value6 */ \
	\
		f2 -= XMM_BIGVAL1;/* Remove rounding constant */ \
		f2 *= f64ptr(rbp+128+16+8);	/* new value7 = val * two-to-phi */ \
		f2 += f64ptr(rsi+64+8);	/* Add in FFT data */ \
		f64ptr(rsi+64+8) = f2;	/* Save value7 */ \
	} \
}




/* *************** 2D normalized add/sub macro ****************** */ \
/* This macro adds or subtracts, then "normalizes" eight FFT */ \
/* data values.  This involves multiplying the summed values by */ \
/* two-to-minus-phi.  Rounding the value to an integer.  Making sure */ \
/* the integer is smaller than the maximum allowable integer, generating */ \
/* a carry if necessary. Finally, the value is multiplied by two-to-phi */ \
/* and stored. */ \
/* rcx = pointer to the first number */ \
/* rdx = pointer to the second number */ \
/* rsi = pointer to destination */ \
/* rdi = pointer to array of big vs. little flags */ \
/* rbx = pointer to two-to-phi column multipliers */ \
/* rbp = pointer two-to-phi group multipliers */ \
/* XMM_TMP1,XMM_TMP2,XMM_TMP3,XMM_TMP4 = carries */ \
/* A pipelined version of this code: */ \
/*	xload	xmm0, [rdx]		;; Load second number */ \
/*	fop	xmm0, [rcx]		;; Add/sub first number */ \
/*	mov	mem_loc, ecx */ \
/*	movzx	rcx, BYTE PTR [rdi]	;; Load big vs. little flags */ \
/*	xload	xmm2, [rbx]		;; col two-to-minus-phi */ \
/*	mulpd	xmm2, XMM_TTMP_FUDGE[rcx];; Mul by fudge two-to-minus-phi */ \
/*	mulpd	xmm0, [rax]		;; Mul by grp two-to-minus-phi */ \
/*	mulpd	xmm0, xmm2		;; Mul by fudged col two-to-minus-phi */ \
/*	addpd	xmm0, [rbp+0*16]	;; x1 = values + carry */ \
/*	xload	xmm2, XMM_LIMIT_BIGMAX[rcx];; Load maximum * BIGVAL - BIGVAL */ \
/*	addpd	xmm2, xmm0		;; y1 = top bits of x */ \
/*	xload	xmm6, XMM_LIMIT_BIGMAX_NEG[rcx];; Load -(maximum*BIGVAL-BIGVAL) */ \
/*	addpd	xmm6, xmm2		;; z1 = y1-(maximum * BIGVAL - BIGVAL) */ \
/*	subpd	xmm0, xmm6		;; rounded value = x1 - z1 */ \
/*	mulpd	xmm2, XMM_LIMIT_INVERSE[rcx];; next carry = shifted y1 */ \
/*	xload	xmm4, [rbx]		;; col two-to-phi */ \
/*	mulpd	xmm4, XMM_TTP_FUDGE[rcx];; mul by fudge two-to-phi */ \
/*	mulpd	xmm0, [rax+0*32+16]	;; new value1 = val * grp two-to-phi */ \
/*	mulpd	xmm0, xmm4		;; new value1 *= fudged col two-to-phi */ \
/*	xstore	[rsi+0*dist1], xmm0	;; Save new value1 */ \
/*	xstore	[rbp+0*16], xmm2	;; Save carry */ \
/*	mov	ecx, mem_loc */ \
\

#define xnorm_op_2d(fop, ttp, base2, mem_loc, TMP1, TMP2, TMP3, TMP4) \
	xmm0 = xptr(rdx);		/* Load second number */ \
	fop(xmm0, xptr(rcx));		/* Add/sub first number */ \
ttp(xmm2 = xptr(rbx));		/* Col two-to-minus-phi */ \
ttp(xmm2 *= g->u.xmm.XMM_NORM012_FF);	/* Mul by FFTLEN/2 */ \
ttp(xmm0 *= xmm2);		/* Mul values1 by col two-to-minus-phi */ \
	xmm1 = xptr(rdx+16);		/* Load second number */ \
	fop(xmm1, xptr(rcx+16));		/* Add/sub first number */ \
ttp(xmm1 *= xmm2);		/* Mul values2 by col two-to-minus-phi */ \
	xmm6 = xptr(rdx+32);		/* Load second number */ \
	fop(xmm6, xptr(rcx+32));		/* Add/sub first number */ \
ttp(xmm6 *= xmm2);		/* Mul values3 by col two-to-minus-phi */ \
	xmm7 = xptr(rdx+48);		/* Load second number */ \
	fop(xmm7, xptr(rcx+48));		/* Add/sub first number */ \
ttp(xmm7 *= xmm2);		/* Mul values4 by col two-to-minus-phi */ \
	mem_loc = rdx;		/* Save second src ptr */ \
ttp(rax = u8ptr(rdi));	/* Load big vs. little flags */ \
ttp(rdx = u8ptr(rdi+1)); \
ttp(xmm0 *= xptr(rbp+0*32));	/* Mul by grp two-to-minus-phi */ \
ttp(xmm0 *= xptr2(g->u.xmm.XMM_TTMP_FUDGE,rax));/* Mul by fudge two-to-minus-phi */ \
	xmm0 += TMP1;		/* x1 = values + carry */ \
ttp(xmm1 *= xptr(rbp+1*32));	/* Mul by grp two-to-minus-phi */ \
ttp(xmm1 *= xptr2(g->u.xmm.XMM_TTMP_FUDGE,rdx));/* Mul by fudge two-to-minus-phi */ \
	xmm1 += TMP2;		/* x1 = values + carry */ \
	rounding_interleaved(ttp, base2, noexec, noexec, xmm0, xmm2, xmm4, rax, xmm1, xmm3, xmm5, rdx); \
	TMP1 = xmm2;		/* Save carry */ \
	TMP2 = xmm3;		/* Save carry */ \
ttp(xmm2 = xptr(rbx+16));		/* col two-to-phi */ \
ttp(xmm2 *= xptr2(g->u.xmm.XMM_TTP_FUDGE,rax));/* mul by fudge two-to-phi */ \
ttp(xmm0 *= xptr(rbp+0*32+16));	/* new value1 = val * grp two-to-phi */ \
ttp(xmm0 *= xmm2);		/* new value1 *= fudged col two-to-phi */ \
	xptr(rsi+0*16) = xmm0;	/* Save new value1 */ \
ttp(xmm3 = xptr(rbx+16));		/* col two-to-phi */ \
ttp(xmm3 *= xptr2(g->u.xmm.XMM_TTP_FUDGE,rdx));/* mul by fudge two-to-phi */ \
ttp(xmm1 *= xptr(rbp+1*32+16));	/* new value1 = val * grp two-to-phi */ \
ttp(xmm1 *= xmm3);		/* new value1 *= fudged col two-to-phi */ \
	xptr(rsi+1*16) = xmm1;	/* Save new value1 */ \
\
ttp(rax = u8ptr(rdi+2));	/* Load big vs. little flags */ \
ttp(rdx = u8ptr(rdi+3)); \
ttp(xmm6 *= xptr(rbp+2*32));	/* Mul by grp two-to-minus-phi */ \
ttp(xmm6 *= xptr2(g->u.xmm.XMM_TTMP_FUDGE,rax));/* Mul by fudge two-to-minus-phi */ \
	xmm6 += TMP3;		/* x1 = values + carry */ \
ttp(xmm7 *= xptr(rbp+3*32));	/* Mul by grp two-to-minus-phi */ \
ttp(xmm7 *= xptr2(g->u.xmm.XMM_TTMP_FUDGE,rdx));/* Mul by fudge two-to-minus-phi */ \
	xmm7 += TMP4;		/* x1 = values + carry */ \
	rounding_interleaved(ttp, base2, noexec, noexec, xmm6, xmm2, xmm4, rax, xmm7, xmm3, xmm5, rdx); \
	TMP3 = xmm2;		/* Save carry */ \
	TMP4 = xmm3;		/* Save carry */ \
ttp(xmm2 = xptr(rbx+16));		/* col two-to-phi */ \
ttp(xmm2 *= xptr2(g->u.xmm.XMM_TTP_FUDGE,rax));/* mul by fudge two-to-phi */ \
ttp(xmm6 *= xptr(rbp+2*32+16));	/* new value1 = val * grp two-to-phi */ \
ttp(xmm6 *= xmm2);		/* new value1 *= fudged col two-to-phi */ \
	xptr(rsi+2*16) = xmm6;	/* Save new value1 */ \
ttp(xmm3 = xptr(rbx+16));		/* col two-to-phi */ \
ttp(xmm3 *= xptr2(g->u.xmm.XMM_TTP_FUDGE,rdx));/* mul by fudge two-to-phi */ \
ttp(xmm7 *= xptr(rbp+3*32+16));	/* new value1 = val * grp two-to-phi */ \
ttp(xmm7 *= xmm3);		/* new value1 *= fudged col two-to-phi */ \
	xptr(rsi+3*16) = xmm7;	/* Save new value1 */ \
	rdx = mem_loc;		/* Restore second src ptr */ \
ttp(rdi += 4);			/* Next flags ptr */ \
ttp(rbx += 32);			/* Next column two-to-phi ptr */ \
	rsi += 64;			/* Next dest ptr */ \
	rdx += 64;			/* Next src ptr */ \
	rcx += 64;			/* Next src ptr */ \
\
\
/* *************** 2D followup macros ****************** */ \
/* This macro finishes the normalize add/sub process by adding four carries */ \
/* from the end of a block back to the start of the block.  The remaining */ \
/* four carries are rotated for starting the next block. */ \
\


#define xnorm_op_2d_blk(srcreg, screg, carry1, carry2, carry3, carry4) \
	xmm2 = XMM_BIGVAL2; \
\
	xmm0 = carry1;		/* Load carries */ \
	xmm0[0] -= XMM_BIGVAL1; \
	xmm0[0] *= f64ptr(screg+0*32+24);	/* mul by grp two-to-phi */ \
	xmm0[0] += f64ptr(srcreg+0*16+8);	/* Add in FFT data */ \
	f64ptr(srcreg+0*16+8) = xmm0[0];	/* Save FFT word */ \
	xmm0[0] = xmm0[1]; xmm0[1] = XMM_BIGVAL2[1]; /*shufpd	xmm0, xmm2, 1;*/ \
	carry1 = xmm0;		/* Save carries for start of next blk */ \
\
	xmm0 = carry2;		/* Load carries */ \
	xmm0[0] -= XMM_BIGVAL1; \
	xmm0[0] *= f64ptr(screg+1*32+24);	/* mul by grp two-to-phi */ \
	xmm0[0] += f64ptr(srcreg+1*16+8);	/* Add in FFT data */ \
	f64ptr(srcreg+1*16+8) = xmm0[0];	/* Save FFT word */ \
	xmm0[0] = xmm0[1]; xmm0[1] = XMM_BIGVAL2[1]; \
	carry2 = xmm0;		/* Save carries for start of next blk */ \
\
	xmm0 = carry3;		/* Load carries */ \
	xmm0[0] -= XMM_BIGVAL1; \
	xmm0[0] *= f64ptr(screg+2*32+24);	/* mul by grp two-to-phi */ \
	xmm0[0] += f64ptr(srcreg+2*16+8);	/* Add in FFT data */ \
	f64ptr(srcreg+2*16+8) = xmm0[0];	/* Save FFT word */ \
	xmm0[0] = xmm0[1]; xmm0[1] = XMM_BIGVAL2[1]; \
	carry3 = xmm0;		/* Save carries for start of next blk */ \
\
	xmm0 = carry4;		/* Load carries */ \
	xmm0[0] -= XMM_BIGVAL1; \
	xmm0[0] *= f64ptr(screg+3*32+24);	/* mul by grp two-to-phi */ \
	xmm0[0] += f64ptr(srcreg+3*16+8);	/* Add in FFT data */ \
	f64ptr(srcreg+3*16+8) = xmm0[0];	/* Save FFT word */ \
	xmm0[0] = xmm0[1]; xmm0[1] = XMM_BIGVAL2[1]; \
	carry4 = xmm0;		/* Save carries for start of next blk */ \


/* This macro finishes the normalize add/sub process by adding two carries */ \
/* from the end of a section back to the start of the section.  The remaining */ \
/* two carries are rotated for starting the next section. */ \
/* rax = pointer to FFT data */ \
/* rbx = pointer two-to-phi group multipliers */ \
\

#define xnorm_op_2d_sec(carry1, carry2, carry3, carry4) \
	xmm2 = XMM_BIGVAL2; \
\
	xmm0[0] = carry1[0]; /* Load carry */ \
	xmm0[0] -= XMM_BIGVAL1; \
	xmm0[0] *= f64ptr(rbx+1*32+16);	/* mul by grp two-to-phi */ \
	xmm0[0] += f64ptr(rax+1*16);	/* Add in FFT data */ \
	f64ptr(rax+1*16) = xmm0[0];	/* Save FFT word */ \
\
	xmm1 = carry2; \
	carry1 = xmm1; \
	carry2 = XMM_BIGVAL2;		/* Save carries for start of next sec */ \
\
	xmm0[0] = carry3[0];		/* Load carry */ \
	xmm0[0] -= XMM_BIGVAL1; \
	xmm0[0] *= f64ptr(rbx+3*32+16);	/* mul by grp two-to-phi */ \
	xmm0[0] += f64ptr(rax+3*16);	/* Add in FFT data */ \
	f64ptr(rax+3*16) = xmm0[0];	/* Save FFT word */ \
\
	xmm1 = carry4; \
	carry3 = xmm1; \
	carry4 = XMM_BIGVAL2;		/* Save carries for start of next sec */ \
\
\
\
/* This macro finishes the normalize add/sub process by adding the final */ \
/* two carries back into the appropriate FFT values at the start of the fft. */ \
/* rsi = pointer to FFT data */ \
/* rbp = pointer two-to-phi group multipliers */ \
/* xmm6 = non-wraparound carry */ \
/* xmm7 = wraparound carry */ \
\

#define xnorm_op_2d_fft \
	xmm2[0] = XMM_BIGVAL1; xmm2[1] = 0.0; \
\
	xmm7[0] -= xmm2[0];		/* Remove XMM_BIGVAL */ \
	xmm7[0] *= g->u.xmm.XMM_MINUS_C[0];	/* mul wrap around carry by -c */ \
	xmm7[0] *= f64ptr(rbp+0*32+16);	/* mul by grp two-to-phi */ \
	xmm7[0] += f64ptr(rsi+0*16);	/* Add in FFT data */ \
	f64ptr(rsi+0*16) = xmm7[0];	/* Save FFT word */ \
\
	xmm6[0] -= xmm2[0];		/* Remove XMM_BIGVAL */ \
	xmm6[0] *= f64ptr(rbp+2*32+16);	/* mul by grp two-to-phi */ \
	xmm6[0] += f64ptr(rsi+2*16);	/* Add in FFT data */ \
	f64ptr(rsi+2*16) = xmm6[0];	/* Save FFT word */ \
\
\
\
/* *************** 2D normalized add & sub macro ****************** */ \
/* This macro adds and subtracts, then "normalizes" eight FFT */ \
/* data values.  This involves multiplying the summed values by */ \
/* two-to-minus-phi.  Rounding the value to an integer.  Making sure */ \
/* the integer is smaller than the maximum allowable integer, generating */ \
/* a carry if necessary. Finally, the value is multiplied by two-to-phi */ \
/* and stored. */ \
/* rcx = pointer to the first number */ \
/* rdx = pointer to the second number */ \
/* rsi = pointer to destination #1 */ \
/* rbp = pointer to destination #2 */ \
/* rbx = pointer two-to-phi group multipliers */ \
/* rdi = pointer to array of big vs. little flags */ \
/* rax = pointer to two-to-phi column multipliers */ \
/* XMM_TMP1,XMM_TMP2,XMM_TMP3,XMM_TMP4 = addition carries */ \
/* XMM_TMP5,XMM_TMP6,XMM_TMP7,XMM_TMP8 = subtraction carries */ \
\

#define xnorm_addsub_2d(ttp, base2, mem_loc, TMP1, TMP2, TMP3, TMP4, TMP5, TMP6, TMP7, TMP8) \
	xmm1 = xptr(rcx+0*16);	/* Load first number */ \
	xmm0 = xptr(rdx+0*16);	/* Load second number */ \
	xmm1 -= xmm0;		/* first - second number */ \
	xmm0 += xptr(rcx+0*16);	/* first + second number */ \
ttp(xmm2 = xptr(rax));		/* Col two-to-minus-phi */ \
ttp(xmm2 *= g->u.xmm.XMM_NORM012_FF);	/* Mul by FFTLEN/2 */ \
ttp(xmm0 *= xmm2);		/* Mul values1 by col two-to-minus-phi */ \
ttp(xmm1 *= xmm2);		/* Mul values1 by col two-to-minus-phi */ \
\
	xmm7 = xptr(rcx+1*16);	/* Load first number */ \
	xmm6 = xptr(rdx+1*16);	/* Load second number */ \
	xmm7 -= xmm6;		/* first - second number */ \
	xmm6 += xptr(rcx+1*16);	/* first + second number */ \
ttp(xmm6 *= xmm2);		/* Mul values1 by col two-to-minus-phi */ \
ttp(xmm7 *= xmm2);		/* Mul values1 by col two-to-minus-phi */ \
\
	mem_loc = rcx;		/* Save first src ptr */ \
\
ttp(rcx = u8ptr(rdi+0));	/* Load big vs. little flags */ \
ttp(xmm0 *= xptr(rbx+0*32));	/* Mul by grp two-to-minus-phi */ \
ttp(xmm0 *= xptr2(g->u.xmm.XMM_TTMP_FUDGE,rcx));/* Mul by fudge two-to-minus-phi */ \
	xmm0 += TMP1;		/* x1 = values + carry */ \
ttp(xmm1 *= xptr(rbx+0*32));	/* Mul by grp two-to-minus-phi */ \
ttp(xmm1 *= xptr2(g->u.xmm.XMM_TTMP_FUDGE,rcx));/* Mul by fudge two-to-minus-phi */ \
	xmm1 += TMP5;		/* x1 = values + carry */ \
	rounding_interleaved(ttp, base2, noexec, noexec, xmm0, xmm2, xmm4, rcx, xmm1, xmm3, xmm5, rcx); \
	TMP1 = xmm2;		/* Save carry */ \
	TMP5 = xmm3;		/* Save carry */ \
ttp(xmm2 = xptr(rax+16));		/* col two-to-phi */ \
ttp(xmm2 *= xptr2(g->u.xmm.XMM_TTP_FUDGE,rcx));/* mul by fudge two-to-phi */ \
ttp(xmm0 *= xptr(rbx+0*32+16));	/* new value1 = val * grp two-to-phi */ \
ttp(xmm0 *= xmm2);		/* new value1 *= fudged col two-to-phi */ \
	xptr(rsi+0*16) = xmm0;	/* Save new value1 */ \
ttp(xmm3 = xptr(rax+16));		/* col two-to-phi */ \
ttp(xmm3 *= xptr2(g->u.xmm.XMM_TTP_FUDGE,rcx));/* mul by fudge two-to-phi */ \
ttp(xmm1 *= xptr(rbx+0*32+16));	/* new value1 = val * grp two-to-phi */ \
ttp(xmm1 *= xmm3);		/* new value1 *= fudged col two-to-phi */ \
	xptr(rbp+0*16) = xmm1;	/* Save new value1 */ \
\
ttp(rcx = u8ptr(rdi+1));	/* Load big vs. little flags */ \
ttp(xmm6 *= xptr(rbx+1*32));	/* Mul by grp two-to-minus-phi */ \
ttp(xmm6 *= xptr2(g->u.xmm.XMM_TTMP_FUDGE,rcx));/* Mul by fudge two-to-minus-phi */ \
	xmm6 += TMP2;		/* x1 = values + carry */ \
ttp(xmm7 *= xptr(rbx+1*32));	/* Mul by grp two-to-minus-phi */ \
ttp(xmm7 *= xptr2(g->u.xmm.XMM_TTMP_FUDGE,rcx));/* Mul by fudge two-to-minus-phi */ \
	xmm7 += TMP6;		/* x1 = values + carry */ \
	rounding_interleaved(ttp, base2, noexec, noexec, xmm6, xmm2, xmm4, rcx, xmm7, xmm3, xmm5, rcx); \
	TMP2 = xmm2;		/* Save carry */ \
	TMP6 = xmm3;		/* Save carry */ \
ttp(xmm2 = xptr(rax+16));		/* col two-to-phi */ \
ttp(xmm2 *= xptr2(g->u.xmm.XMM_TTP_FUDGE,rcx));/* mul by fudge two-to-phi */ \
ttp(xmm6 *= xptr(rbx+1*32+16));	/* new value1 = val * grp two-to-phi */ \
ttp(xmm6 *= xmm2);		/* new value1 *= fudged col two-to-phi */ \
	xptr(rsi+1*16) = xmm6;	/* Save new value1 */ \
ttp(xmm3 = xptr(rax+16));		/* col two-to-phi */ \
ttp(xmm3 *= xptr2(g->u.xmm.XMM_TTP_FUDGE,rcx));/* mul by fudge two-to-phi */ \
ttp(xmm7 *= xptr(rbx+1*32+16));	/* new value1 = val * grp two-to-phi */ \
ttp(xmm7 *= xmm3);		/* new value1 *= fudged col two-to-phi */ \
	xptr(rbp+1*16) = xmm7;	/* Save new value1 */ \
\
	rcx = mem_loc;		/* Restore first src ptr */ \
\
	xmm1 = xptr(rcx+2*16);	/* Load first number */ \
	xmm0 = xptr(rdx+2*16);	/* Load second number */ \
	xmm1 -= xmm0;		/* first - second number */ \
	xmm0 += xptr(rcx+2*16);	/* first + second number */ \
ttp(xmm2 = xptr(rax));		/* Col two-to-minus-phi */ \
ttp(xmm2 *= g->u.xmm.XMM_NORM012_FF);	/* Mul by FFTLEN/2 */ \
ttp(xmm0 *= xmm2);		/* Mul values1 by col two-to-minus-phi */ \
ttp(xmm1 *= xmm2);		/* Mul values1 by col two-to-minus-phi */ \
\
	xmm7 = xptr(rcx+3*16);	/* Load first number */ \
	xmm6 = xptr(rdx+3*16);	/* Load second number */ \
	xmm7 -= xmm6;		/* first - second number */ \
	xmm6 += xptr(rcx+3*16);	/* first + second number */ \
ttp(xmm6 *= xmm2);		/* Mul values1 by col two-to-minus-phi */ \
ttp(xmm7 *= xmm2);		/* Mul values1 by col two-to-minus-phi */ \
\
	mem_loc = rcx;		/* Save first src ptr */ \
\
ttp(rcx = u8ptr(rdi+2));	/* Load big vs. little flags */ \
ttp(xmm0 *= xptr(rbx+2*32));	/* Mul by grp two-to-minus-phi */ \
ttp(xmm0 *= xptr2(g->u.xmm.XMM_TTMP_FUDGE,rcx));/* Mul by fudge two-to-minus-phi */ \
	xmm0 += TMP3;		/* x1 = values + carry */ \
ttp(xmm1 *= xptr(rbx+2*32));	/* Mul by grp two-to-minus-phi */ \
ttp(xmm1 *= xptr2(g->u.xmm.XMM_TTMP_FUDGE,rcx));/* Mul by fudge two-to-minus-phi */ \
	xmm1 += TMP7;		/* x1 = values + carry */ \
	rounding_interleaved(ttp, base2, noexec, noexec, xmm0, xmm2, xmm4, rcx, xmm1, xmm3, xmm5, rcx); \
	TMP3 = xmm2;		/* Save carry */ \
	TMP7 = xmm3;		/* Save carry */ \
ttp(xmm2 = xptr(rax+16));		/* col two-to-phi */ \
ttp(xmm2 *= xptr2(g->u.xmm.XMM_TTP_FUDGE,rcx));/* mul by fudge two-to-phi */ \
ttp(xmm0 *= xptr(rbx+2*32+16));	/* new value1 = val * grp two-to-phi */ \
ttp(xmm0 *= xmm2);		/* new value1 *= fudged col two-to-phi */ \
	xptr(rsi+2*16) = xmm0;	/* Save new value1 */ \
ttp(xmm3 = xptr(rax+16));		/* col two-to-phi */ \
ttp(xmm3 *= xptr2(g->u.xmm.XMM_TTP_FUDGE,rcx));/* mul by fudge two-to-phi */ \
ttp(xmm1 *= xptr(rbx+2*32+16));	/* new value1 = val * grp two-to-phi */ \
ttp(xmm1 *= xmm3);		/* new value1 *= fudged col two-to-phi */ \
	xptr(rbp+2*16) = xmm1;	/* Save new value1 */ \
\
ttp(rcx = u8ptr(rdi+3));	/* Load big vs. little flags */ \
ttp(xmm6 *= xptr(rbx+3*32));	/* Mul by grp two-to-minus-phi */ \
ttp(xmm6 *= xptr2(g->u.xmm.XMM_TTMP_FUDGE,rcx));/* Mul by fudge two-to-minus-phi */ \
	xmm6 += TMP4;		/* x1 = values + carry */ \
ttp(xmm7 *= xptr(rbx+3*32));	/* Mul by grp two-to-minus-phi */ \
ttp(xmm7 *= xptr2(g->u.xmm.XMM_TTMP_FUDGE,rcx));/* Mul by fudge two-to-minus-phi */ \
	xmm7 += TMP8;		/* x1 = values + carry */ \
	rounding_interleaved(ttp, base2, noexec, noexec, xmm6, xmm2, xmm4, rcx, xmm7, xmm3, xmm5, rcx); \
	TMP4 = xmm2;		/* Save carry */ \
	TMP8 = xmm3;		/* Save carry */ \
ttp(xmm2 = xptr(rax+16));		/* col two-to-phi */ \
ttp(xmm2 *= xptr2(g->u.xmm.XMM_TTP_FUDGE,rcx));/* mul by fudge two-to-phi */ \
ttp(xmm6 *= xptr(rbx+3*32+16));	/* new value1 = val * grp two-to-phi */ \
ttp(xmm6 *= xmm2);		/* new value1 *= fudged col two-to-phi */ \
	xptr(rsi+3*16) = xmm6;	/* Save new value1 */ \
ttp(xmm3 = xptr(rax+16));		/* col two-to-phi */ \
ttp(xmm3 *= xptr2(g->u.xmm.XMM_TTP_FUDGE,rcx));/* mul by fudge two-to-phi */ \
ttp(xmm7 *= xptr(rbx+3*32+16));	/* new value1 = val * grp two-to-phi */ \
ttp(xmm7 *= xmm3);		/* new value1 *= fudged col two-to-phi */ \
	xptr(rbp+3*16) = xmm7;	/* Save new value1 */ \
\
	rcx = mem_loc;		/* Restore first src ptr */ \
\
	rdx += 64;			/* Next src ptr */ \
	rcx += 64;			/* Next src ptr */ \
ttp(rdi += 4);			/* Next flags ptr */ \
ttp(rax += 32);			/* Next column two-to-phi ptr */ \
	rsi += 64;			/* Next dest ptr */ \
	rbp += 64;			/* Next dest ptr */ \
\
\
/* *************** 2D normalized small mul macro ****************** */ \
/* This macro multiplies by a small value, then "normalizes" eight FFT */ \
/* data values.  */ \
/* rsi = pointer to destination */ \
/* rdi = pointer to array of big vs. little flags */ \
/* rbx = pointer to two-to-phi column multipliers */ \
/* rbp = pointer two-to-phi group multipliers */ \
/* XMM_TMP1,XMM_TMP2,XMM_TMP3,XMM_TMP4 = carries */ \
/* XMM_TMP5 = small value * optional FFTLEN/2 */ \
/* A pipelined version of this code: */ \
/*	xload	xmm0, [rsi]		;; Load second number */ \
/*	mulpd	xmm0, XMM_TMP5		;; Mul by small value */ \
/*	movzx	rcx, BYTE PTR [rdi]	;; Load big vs. little flags */ \
/*	xload	xmm2, [rbx]		;; col two-to-minus-phi */ \
/*	mulpd	xmm2, XMM_TTMP_FUDGE[rcx];; Mul by fudge two-to-minus-phi */ \
/*	mulpd	xmm0, [rax]		;; Mul by grp two-to-minus-phi */ \
/*	mulpd	xmm0, xmm2		;; Mul by fudged col two-to-minus-phi */ \
/*	addpd	xmm0, [rbp+0*16]	;; x1 = values + carry */ \
/*	xload	xmm2, XMM_LIMIT_BIGMAX[rcx];; Load maximum * BIGVAL - BIGVAL */ \
/*	addpd	xmm2, xmm0		;; y1 = top bits of x */ \
/*	xload	xmm6, XMM_LIMIT_BIGMAX_NEG[rcx];; Load -(maximum*BIGVAL-BIGVAL) */ \
/*	addpd	xmm6, xmm2		;; z1 = y1-(maximum * BIGVAL - BIGVAL) */ \
/*	subpd	xmm0, xmm6		;; rounded value = x1 - z1 */ \
/*	mulpd	xmm2, XMM_LIMIT_INVERSE[rcx];; next carry = shifted y1 */ \
/*	xload	xmm4, [rbx]		;; col two-to-phi */ \
/*	mulpd	xmm4, XMM_TTP_FUDGE[rcx];; mul by fudge two-to-phi */ \
/*	mulpd	xmm0, [rax+0*32+16]	;; new value1 = val * grp two-to-phi */ \
/*	mulpd	xmm0, xmm4		;; new value1 *= fudged col two-to-phi */ \
/*	xstore	[rsi], xmm0		;; Save new value1 */ \
/*	xstore	[rbp+0*16], xmm2	;; Save carry */ \
\

#define xnorm_smallmul_2d(ttp, base2, TMP1, TMP2, TMP3, TMP4, TMP5) \
	xmm2 = TMP5;		/* Load small value * FFTLEN/2 */ \
ttp(xmm2 *= xptr(rbx));		/* Mul by col two-to-minus-phi */ \
	xmm0 = xptr(rsi);		/* Load values1 */ \
	xmm0 *= xmm2;		/* Mul values1 by small value * col two-to-minus-phi */ \
	xmm1 = xptr(rsi+16);		/* Load values2 */ \
	xmm1 *= xmm2;		/* Mul values2 by small value * col two-to-minus-phi */ \
	xmm6 = xptr(rsi+32);		/* Load values3 */ \
	xmm6 *= xmm2;		/* Mul values3 by small value * col two-to-minus-phi */ \
	xmm7 = xptr(rsi+48);		/* Load values4 */ \
	xmm7 *= xmm2;		/* Mul values4 by small value * col two-to-minus-phi */ \
ttp(rax = u8ptr(rdi));	/* Load big vs. little flags */ \
ttp(rcx = u8ptr(rdi+1)); \
ttp(xmm0 *= xptr(rbp+0*32));	/* Mul by grp two-to-minus-phi */ \
ttp(xmm0 *= xptr2(g->u.xmm.XMM_TTMP_FUDGE,rax));/* Mul by fudge two-to-minus-phi */ \
	xmm0 += TMP1;		/* x1 = values + carry */ \
ttp(xmm1 *= xptr(rbp+1*32));	/* Mul by grp two-to-minus-phi */ \
ttp(xmm1 *= xptr2(g->u.xmm.XMM_TTMP_FUDGE,rcx));/* Mul by fudge two-to-minus-phi */ \
	xmm1 += TMP2;		/* x1 = values + carry */ \
	rounding_interleaved(ttp, base2, noexec, noexec, xmm0, xmm2, xmm4, rax, xmm1, xmm3, xmm5, rcx); \
	TMP1 = xmm2;		/* Save carry */ \
	TMP2 = xmm3;		/* Save carry */ \
ttp(xmm2 = xptr(rbx+16));		/* col two-to-phi */ \
ttp(xmm2 *= xptr2(g->u.xmm.XMM_TTP_FUDGE,rax));/* mul by fudge two-to-phi */ \
ttp(xmm0 *= xptr(rbp+0*32+16));	/* new value1 = val * grp two-to-phi */ \
ttp(xmm0 *= xmm2);		/* new value1 *= fudged col two-to-phi */ \
	xptr(rsi+0*16) = xmm0;	/* Save new value1 */ \
ttp(xmm3 = xptr(rbx+16));		/* col two-to-phi */ \
ttp(xmm3 *= xptr2(g->u.xmm.XMM_TTP_FUDGE,rcx));/* mul by fudge two-to-phi */ \
ttp(xmm1 *= xptr(rbp+1*32+16));	/* new value1 = val * grp two-to-phi */ \
ttp(xmm1 *= xmm3);		/* new value1 *= fudged col two-to-phi */ \
	xptr(rsi+1*16) = xmm1;	/* Save new value1 */ \
\
ttp(rax = u8ptr(rdi+2));	/* Load big vs. little flags */ \
ttp(rcx = u8ptr(rdi+3)); \
ttp(xmm6 *= xptr(rbp+2*32));	/* Mul by grp two-to-minus-phi */ \
ttp(xmm6 *= xptr2(g->u.xmm.XMM_TTMP_FUDGE,rax));/* Mul by fudge two-to-minus-phi */ \
	xmm6 += TMP3;		/* x1 = values + carry */ \
ttp(xmm7 *= xptr(rbp+3*32));	/* Mul by grp two-to-minus-phi */ \
ttp(xmm7 *= xptr2(g->u.xmm.XMM_TTMP_FUDGE,rcx));/* Mul by fudge two-to-minus-phi */ \
	xmm7 += TMP4;		/* x1 = values + carry */ \
	rounding_interleaved(ttp, base2, noexec, noexec, xmm6, xmm2, xmm4, rax, xmm7, xmm3, xmm5, rcx); \
	TMP3 = xmm2;		/* Save carry */ \
	TMP4 = xmm3;		/* Save carry */ \
ttp(xmm2 = xptr(rbx+16));		/* col two-to-phi */ \
ttp(xmm2 *= xptr2(g->u.xmm.XMM_TTP_FUDGE,rax));/* mul by fudge two-to-phi */ \
ttp(xmm6 *= xptr(rbp+2*32+16));	/* new value1 = val * grp two-to-phi */ \
ttp(xmm6 *= xmm2);		/* new value1 *= fudged col two-to-phi */ \
	xptr(rsi+2*16) = xmm6;	/* Save new value1 */ \
ttp(xmm3 = xptr(rbx+16));		/* col two-to-phi */ \
ttp(xmm3 *= xptr2(g->u.xmm.XMM_TTP_FUDGE,rcx));/* mul by fudge two-to-phi */ \
ttp(xmm7 *= xptr(rbp+3*32+16));	/* new value1 = val * grp two-to-phi */ \
ttp(xmm7 *= xmm3);		/* new value1 *= fudged col two-to-phi */ \
	xptr(rsi+3*16) = xmm7;	/* Save new value1 */ \
ttp(rdi += 4);			/* Next flags ptr */ \
ttp(rbx += 32);			/* Next column two-to-phi ptr */ \
	rsi += 64;			/* Next dest ptr */ \
\
\
/* This macro finishes the smallmul normalize process by adding four carries */ \
/* from the end of a block back to the start of the block.  The remaining */ \
/* four carries are rotated for starting the next block. */ \
/* rsi = Source ptr to start of block */ \
/* rbp = Group ptr of two-to-phi */ \
/* rbx = Col ptr of two-to-phi */ \
/* rdi = biglit array ptr */ \
/* rcx is destroyed */ \
\

#define xnorm_smallmul_2d_blk(base2, TMP1, TMP2, TMP3, TMP4) \
	xmm5[0] = XMM_BIGVAL1; xmm5[1] = 0.0; \
\
	xmm1[0] = TMP1[0];	/* Load carry */ \
	xnorm_smallmul_2d_prop4(base2, rsi+0*16, 8, rbp+0*32, rbx, rdi+0*1, xmm1[0]); \
	xmm1 = TMP1;		/* Load carries (prop4 will zero high word of xmm1) */ \
	xmm1[0] = xmm1[1]; xmm1[1] = XMM_BIGVAL1; /*shufpd	xmm1, xmm5, 1;*/ \
	TMP1 = xmm1;		/* Save carries for start of next blk */ \
\
	xmm1[0] = TMP2[0];	/* Load carry */ \
	xnorm_smallmul_2d_prop4(base2, rsi+1*16, 8, rbp+1*32, rbx, rdi+1*1, xmm1[0]); \
	xmm1 = TMP2;		/* Load carries */ \
	xmm1[0] = xmm1[1]; xmm1[1] = XMM_BIGVAL1; \
	TMP2 = xmm1;		/* Save carries for start of next blk */ \
\
	xmm1[0] = TMP3[0];	/* Load carry */ \
	xnorm_smallmul_2d_prop4(base2, rsi+2*16, 8, rbp+2*32, rbx, rdi+2*1, xmm1[0]); \
	xmm1 = TMP3;		/* Load carries */ \
	xmm1[0] = xmm1[1]; xmm1[1] = XMM_BIGVAL1; \
	TMP3 = xmm1;		/* Save carries for start of next blk */ \
\
	xmm1[0] = TMP4[0];	/* Load carry */ \
	xnorm_smallmul_2d_prop4(base2, rsi+3*16, 8, rbp+3*32, rbx, rdi+3*1, xmm1[0]); \
	xmm1 = TMP4;		/* Load carries */ \
	xmm1[0] = xmm1[1]; xmm1[1] = XMM_BIGVAL1; \
	TMP4 = xmm1;		/* Save carries for start of next blk */ \
\
\
/* This macro finishes the smallmul normalize process by adding two carries */ \
/* from the end of a section back to the start of the section.  The remaining */ \
/* two carries are rotated for starting the next section. */ \
/* rsi = Source ptr to start of block */ \
/* rbp = Group ptr of two-to-phi */ \
/* rbx = Col ptr of two-to-phi */ \
/* rdi = biglit array ptr */ \
/* rcx is destroyed */ \
\

#define xnorm_smallmul_2d_sec(base2, TMP1, TMP2, TMP3, TMP4) \
	xmm5 = XMM_BIGVAL2; \
\
	xmm1[0] = TMP1[0];	/* Load carry */ \
	xnorm_smallmul_2d_prop4(base2, rsi+1*16, 0, rbp+1*32, rbx, rdi+1*1, xmm1[0]); \
	xmm1 = TMP2; \
	TMP1 = xmm1; \
	TMP2 = XMM_BIGVAL2;		/* Save carries for start of next sec */ \
\
	xmm1[0] = TMP3[0];	/* Load carry */ \
	xnorm_smallmul_2d_prop4(base2, rsi+3*16, 0, rbp+3*32, rbx, rdi+3*1, xmm1[0]); \
	xmm1 = TMP4; \
	TMP3 = xmm1; \
	TMP4 = XMM_BIGVAL2;		/* Save carries for start of next sec */ \
\
\
/* This macro finishes the smallmul normalize process by adding the final */ \
/* two carries back into the appropriate FFT values at the start of the fft. */ \
/* rsi = pointer to FFT data */ \
/* rbp = pointer two-to-phi group multipliers */ \
/* rdi = pointer to big/little flags */ \
/* xmm6 = non-wraparound carry */ \
/* xmm7 = wraparound carry */ \
\

#define xnorm_smallmul_2d_fft(base2) \
	xmm7[0] -= XMM_BIGVAL1; \
	xmm7[0] *= g->u.xmm.XMM_MINUS_C[0];	/* Negate the carry */ \
	xmm7[0] += XMM_BIGVAL1; \
	xnorm_smallmul_2d_prop6(base2, rsi+0*16, 0, rbp+0*32, rbx, rdi+0*1, xmm7[0]); \
	xnorm_smallmul_2d_prop4(base2, rsi+2*16, 0, rbp+2*32, rbx, rdi+2*1, xmm6[0]); \
\
\
/* Propagate a single gwsmallmul carry across 4 words. */ \
\

#define xnorm_smallmul_2d_prop4(base2, srcptr, off8, grpptr, colptr, biglit, xcarry) { \
	double f0, xtmp, xaux; \
	rax = u8ptr(biglit);		/* First word */ \
	f0 = f64ptr(srcptr+0*64+off8);	/* Load FFT data */ \
	f0 *= f64ptr(colptr+0*32+off8);	/* Mul by col two-to-minus-phi */ \
	xaux = f64ptr(grpptr+off8);		/* Mul by grp two-to-minus-phi */ \
	xaux *= g->u.xmm.XMM_NORM012_FF[0];		/* Mul by FFTLEN/2 */ \
	f0 *= f64ptr2(g->u.xmm.XMM_TTMP_FUDGE,rax+off8);	/* Mul by fudge two-to-minus-phi */ \
	f0 *= xaux; \
	f0 += xcarry;			/* Add in the carry */ \
	single_rounding(base2, f0, xcarry, xtmp, rax+off8); \
	xtmp = f64ptr(colptr+0*32+16+off8);	/* col two-to-phi */ \
	xtmp *= f64ptr2(g->u.xmm.XMM_TTP_FUDGE,rax+off8);	/* mul by fudge two-to-phi */ \
	f0 *= f64ptr(grpptr+16+off8);	/* new value = val * grp two-to-phi */ \
	f0 *= xtmp;			/* new value *= fudged col two-to-phi */ \
	f64ptr(srcptr+0*64+off8) = f0;	/* Save value1 */ \
\
	rax = u8ptr(biglit+4);	/* Second word */ \
	f0 = f64ptr(srcptr+1*64+off8);	/* Load FFT data */ \
	f0 *= f64ptr(colptr+1*32+off8);	/* Mul by col two-to-minus-phi */ \
	xtmp = f64ptr2(g->u.xmm.XMM_TTMP_FUDGE,rax+off8);/* Mul by fudge two-to-minus-phi */ \
	xtmp *= xaux;			/* Mul by grp * FFTLEN/2 */ \
	f0 *= xtmp; \
	f0 += xcarry;			/* Add in the carry */ \
	single_rounding(base2, f0, xcarry, xtmp, rax+off8); \
	xtmp = f64ptr(colptr+1*32+16+off8);	/* col two-to-phi */ \
	xtmp *= f64ptr2(g->u.xmm.XMM_TTP_FUDGE,rax+off8);	/* mul by fudge two-to-phi */ \
	f0 *= f64ptr(grpptr+16+off8);	/* new value = val * grp two-to-phi */ \
	f0 *= xtmp;			/* new value *= fudged col two-to-phi */ \
	f64ptr(srcptr+1*64+off8) = f0;	/* Save value2 */ \
\
	rcx = g->BIGLIT_INCR2;/* Different clm values step through */ \
;						/* big/lit array differently */ \
	rax = u8ptr(biglit+rcx);	/* Third word */ \
	f0 = f64ptr(srcptr+2*64+off8);	/* Load FFT data */ \
	f0 *= f64ptr(colptr+2*32+off8);	/* Mul by col two-to-minus-phi */ \
	xtmp = f64ptr2(g->u.xmm.XMM_TTMP_FUDGE,rax+off8);/* Mul by fudge two-to-minus-phi */ \
	xtmp *= xaux;			/* Mul by grp * FFTLEN/2 */ \
	f0 *= xtmp; \
	f0 += xcarry;			/* Add in the carry */ \
	single_rounding(base2, f0, xcarry, xtmp, rax+off8); \
	xtmp = f64ptr(colptr+2*32+16+off8);	/* col two-to-phi */ \
	xtmp *= f64ptr2(g->u.xmm.XMM_TTP_FUDGE,rax+off8);	/* mul by fudge two-to-phi */ \
	f0 *= f64ptr(grpptr+16+off8);	/* new value = val * grp two-to-phi */ \
	f0 *= xtmp;			/* new value *= fudged col two-to-phi */ \
	f64ptr(srcptr+2*64+off8) = f0;	/* Save value3 */ \
\
	rax = u8ptr(biglit+rcx+4);	/* Fourth word */ \
	xcarry -= XMM_BIGVAL1;/* Make carry an integer */ \
	xtmp = f64ptr(colptr+3*32+16+off8);	/* col two-to-phi */ \
	xtmp *= f64ptr2(g->u.xmm.XMM_TTP_FUDGE,rax+off8);	/* mul by fudge two-to-phi */ \
	xcarry *= f64ptr(grpptr+16+off8);	/* carry *= grp two-to-phi */ \
	xcarry *= xtmp;			/* carry *= fudged col two-to-phi */ \
	xcarry += f64ptr(srcptr+3*64+off8);	/* Add FFT data */ \
	f64ptr(srcptr+3*64+off8) = xcarry;	/* Save value4 */ \
} \

/* Propagate a single gwsmallmul carry across 6 words. */ \


#define xnorm_smallmul_2d_prop6(base2, srcptr, off8, grpptr, colptr, biglit, xcarry) { \
	double f0, xtmp, xaux; \
	rax = u8ptr(biglit);		/* First word */ \
	f0 = f64ptr(srcptr+0*64+off8);	/* Load FFT data */ \
	f0 *= f64ptr(colptr+0*32+off8);	/* Mul by col two-to-minus-phi */ \
	xaux = f64ptr(grpptr+off8);		/* Mul by grp two-to-minus-phi */ \
	xaux *= g->u.xmm.XMM_NORM012_FF[0];		/* Mul by FFTLEN/2 */ \
	f0 *= f64ptr2(g->u.xmm.XMM_TTMP_FUDGE,rax+off8);/* Mul by fudge two-to-minus-phi */ \
	f0 *= xaux; \
	f0 += xcarry;			/* Add in the carry */ \
	single_rounding(base2, f0, xcarry, xtmp, rax+off8); \
	xtmp = f64ptr(colptr+0*32+16+off8);	/* col two-to-phi */ \
	xtmp *= f64ptr2(g->u.xmm.XMM_TTP_FUDGE,rax+off8);	/* mul by fudge two-to-phi */ \
	f0 *= f64ptr(grpptr+16+off8);	/* new value = val * grp two-to-phi */ \
	f0 *= xtmp;			/* new value *= fudged col two-to-phi */ \
	f64ptr(srcptr+0*64+off8) = f0;	/* Save value1 */ \
\
	rax = u8ptr(biglit+4);	/* Second word */ \
	f0 = f64ptr(srcptr+1*64+off8);	/* Load FFT data */ \
	f0 *= f64ptr(colptr+1*32+off8);	/* Mul by col two-to-minus-phi */ \
	xtmp = f64ptr2(g->u.xmm.XMM_TTMP_FUDGE,rax+off8);/* Mul by fudge two-to-minus-phi */ \
	xtmp *= xaux;			/* Mul by grp * FFTLEN/2 */ \
	f0 *= xtmp; \
	f0 += xcarry;			/* Add in the carry */ \
	single_rounding(base2, f0, xcarry, xtmp, rax+off8); \
	xtmp = f64ptr(colptr+1*32+16+off8);	/* col two-to-phi */ \
	xtmp *= f64ptr2(g->u.xmm.XMM_TTP_FUDGE,rax+off8);	/* mul by fudge two-to-phi */ \
	f0 *= f64ptr(grpptr+16+off8);	/* new value = val * grp two-to-phi */ \
	f0 *= xtmp;			/* new value *= fudged col two-to-phi */ \
	f64ptr(srcptr+1*64+off8) = f0;	/* Save value2 */ \
\
	rcx = g->BIGLIT_INCR2;/* Different clm values step through */ \
;						/* big/lit array differently */ \
	rax = u8ptr(biglit+rcx);	/* Third word */ \
	f0 = f64ptr(srcptr+2*64+off8);	/* Load FFT data */ \
	f0 *= f64ptr(colptr+2*32+off8);	/* Mul by col two-to-minus-phi */ \
	xtmp = f64ptr2(g->u.xmm.XMM_TTMP_FUDGE,rax+off8);/* Mul by fudge two-to-minus-phi */ \
	xtmp *= xaux;			/* Mul by grp * FFTLEN/2 */ \
	f0 *= xtmp; \
	f0 += xcarry;			/* Add in the carry */ \
	single_rounding(base2, f0, xcarry, xtmp, rax+off8); \
	xtmp = f64ptr(colptr+2*32+16+off8);	/* col two-to-phi */ \
	xtmp *= f64ptr2(g->u.xmm.XMM_TTP_FUDGE,rax+off8);	/* mul by fudge two-to-phi */ \
	f0 *= f64ptr(grpptr+16+off8);	/* new value = val * grp two-to-phi */ \
	f0 *= xtmp;			/* new value *= fudged col two-to-phi */ \
	f64ptr(srcptr+2*64+off8) = f0;	/* Save value3 */ \
\
	rax = u8ptr(biglit+rcx+4);	/* Fourth word */ \
	f0 = f64ptr(srcptr+3*64+off8);	/* Load FFT data */ \
	f0 *= f64ptr(colptr+3*32+off8);	/* Mul by col two-to-minus-phi */ \
	xtmp = f64ptr2(g->u.xmm.XMM_TTMP_FUDGE,rax+off8);/* Mul by fudge two-to-minus-phi */ \
	xtmp *= xaux;			/* Mul by grp * FFTLEN/2 */ \
	f0 *= xtmp; \
	f0 += xcarry;			/* Add in the carry */ \
	single_rounding(base2, f0, xcarry, xtmp, rax+off8); \
	xtmp = f64ptr(colptr+3*32+16+off8);	/* col two-to-phi */ \
	xtmp *= f64ptr2(g->u.xmm.XMM_TTP_FUDGE,rax+off8);	/* mul by fudge two-to-phi */ \
	f0 *= f64ptr(grpptr+16+off8);	/* new value = val * grp two-to-phi */ \
	f0 *= xtmp;			/* new value *= fudged col two-to-phi */ \
	f64ptr(srcptr+3*64+off8) = f0;	/* Save value4 */ \
\
	rcx = g->BIGLIT_INCR4;/* Different clm values step through */ \
;						/* big/lit array differently */ \
	rax = u8ptr(biglit+rcx);	/* Fifth word */ \
	f0 = f64ptr(srcptr+4*64+off8);	/* Load FFT data */ \
	f0 *= f64ptr(colptr+4*32+off8);	/* Mul by col two-to-minus-phi */ \
	xtmp = f64ptr2(g->u.xmm.XMM_TTMP_FUDGE,rax+off8);/* Mul by fudge two-to-minus-phi */ \
	xtmp *= xaux;			/* Mul by grp * FFTLEN/2 */ \
	f0 *= xtmp; \
	f0 += xcarry;			/* Add in the carry */ \
	single_rounding(base2, f0, xcarry, xtmp, rax+off8); \
	xtmp = f64ptr(colptr+4*32+16+off8);	/* col two-to-phi */ \
	xtmp *= f64ptr2(g->u.xmm.XMM_TTP_FUDGE,rax+off8);	/* mul by fudge two-to-phi */ \
	f0 *= f64ptr(grpptr+16+off8);	/* new value = val * grp two-to-phi */ \
	f0 *= xtmp;			/* new value *= fudged col two-to-phi */ \
	f64ptr(srcptr+4*64+off8) = f0;	/* Save value5 */ \
\
	rax = u8ptr(biglit+rcx+4);	/* Sixth word */ \
	xcarry -= XMM_BIGVAL1;/* Make carry an integer */ \
	xtmp = f64ptr(colptr+5*32+16+off8);	/* col two-to-phi */ \
	xtmp *= f64ptr2(g->u.xmm.XMM_TTP_FUDGE,rax+off8);	/* mul by fudge two-to-phi */ \
	xcarry *= f64ptr(grpptr+16+off8);	/* carry *= grp two-to-phi */ \
	xcarry *= xtmp;			/* carry *= fudged col two-to-phi */ \
	xcarry += f64ptr(srcptr+5*64+off8);	/* Add FFT data */ \
	f64ptr(srcptr+5*64+off8) = xcarry;	/* Save value6 */ \
}

/* This macro finishes the smallmul normalize process by adding the final */ \
/* two carries back into the appropriate FFT values at the start of the fft. */ \
/* rsi = pointer to FFT data */ \
/* rbp = pointer two-to-phi group multipliers */ \
/* rbx = pointer two-to-phi col multipliers */ \
/* rdi = pointer to big/little flags */ \
/* xmm6 = non-wraparound carry */ \
/* xmm7 = wraparound carry (this will be zero for smallmul of zero-padded number) */ \


#define xnorm_smallmul_2d_fft_zpad(base2) { \
	double f0, f1, f2, f4; \
;	/* Copy and integerize data from 7 words above halfway point to ZPAD0-ZPAD6 */ \
;	/* Clear words 5,6,7 */ \
;	/* Then we can make an exact copy of most of the xnorm012_2d_zpad code */ \
\
	rax = u8ptr(rdi);	/* Load big vs. little flags */ \
	f0 = f64ptr(rsi+0*64+32);	/* Value1 */ \
	xmm6[0] -= XMM_BIGVAL1;/* Remove XMM_BIGVAL from carry */ \
	single_split_lower_zpad_word(base2, f0, xmm6[0], f4, rax); \
	g->ZPAD0_6[0] = f0; \
\
	rax = u8ptr(rdi+4);	/* Load big vs. little flags */ \
	f0 = f64ptr(rsi+1*64+32);	/* Value2 */ \
	f0 *= f64ptr(rbx+1*32);	/* Mul values2 by col two-to-minus-phi */ \
	f0 *= g->u.xmm.XMM_NORM012_FF[0];	/* Mul by FFTLEN/2 */ \
	single_split_lower_zpad_word(base2, f0, xmm6[0], f4, rax); \
	g->ZPAD0_6[1] = f0; \
\
	rcx = g->BIGLIT_INCR2;/* Different clm values step through */ \
;					/* big/lit array differently */ \
	rax = u8ptr(rdi+rcx);	/* Load big vs. little flags */ \
	f0 = f64ptr(rsi+2*64+32);	/* Value3 */ \
	f0 *= f64ptr(rbx+2*32);	/* Mul values3 by col two-to-minus-phi */ \
	f0 *= g->u.xmm.XMM_NORM012_FF[0];	/* Mul by FFTLEN/2 */ \
	single_split_lower_zpad_word(base2, f0, xmm6[0], f4, rax); \
	g->ZPAD0_6[2] = f0; \
\
	f0 = f64ptr(rsi+3*64+32);	/* Value4 */ \
	f0 *= f64ptr(rbx+3*32);	/* Mul values4 by col two-to-minus-phi */ \
	f0 *= g->u.xmm.XMM_NORM012_FF[0];	/* Mul by FFTLEN/2 */ \
	f0 = roundf64(f0); /* Round to an integer */ \
	f0 += xmm6[0];		/* Value4 + carry */ \
	g->ZPAD0_6[3] = f0; \
\
	f0 = f64ptr(rsi+4*64+32);	/* Value5 */ \
	f0 *= f64ptr(rbx+4*32);	/* Mul values5 by col two-to-minus-phi */ \
	f0 *= g->u.xmm.XMM_NORM012_FF[0];	/* Mul by FFTLEN/2 */ \
	f0 = roundf64(f0); /* Round to an integer */ \
	g->ZPAD0_6[4] = f0; \
\
	f0 = f64ptr(rsi+5*64+32);	/* Value6 */ \
	f0 *= f64ptr(rbx+5*32);	/* Mul values6 by col two-to-minus-phi */ \
	f0 *= g->u.xmm.XMM_NORM012_FF[0];	/* Mul by FFTLEN/2 */ \
	f0 = roundf64(f0); /* Round to an integer */ \
	g->ZPAD0_6[5] = f0; \
\
	f0 = f64ptr(rsi+6*64+32);	/* Value7 */ \
	f0 *= f64ptr(rbx+6*32);	/* Mul values7 by col two-to-minus-phi */ \
	f0 *= g->u.xmm.XMM_NORM012_FF[0];	/* Mul by FFTLEN/2 */ \
	f0 = roundf64(f0); /* Round to an integer */ \
	g->ZPAD0_6[6] = f0; \
\
	f0 = 0;		/* Clear highest words */ \
	f64ptr(rsi+4*64+32) = f0; \
	f64ptr(rsi+5*64+32) = f0; \
	f64ptr(rsi+6*64+32) = f0; \
\
;	/* Divide the zpad data by k.  Store the integer part in XMM_TMP */ \
;	/* and the remainder in ZPAD0.  Later we will wrap the integer part */ \
;	/* down to the bottom of the FFT data area (and multiply by -c). */ \
;	/* And we will store the remainder in the upper half of the FFT */ \
;	/* data area. */ \
\
;	/* Note there are three cases to handle.  K is smaller than a big word. */ \
;	/* K is between one and 2 big words in size.  And K is more than */ \
;	/* 2 big words in size. */ \
\
	xfunc572(g); \
\
;	/* Now normalize the data above the halfway point.  Remember that the */ \
;	/* column two-to-phi multiplier for the first value will be 1.0. */ \
\
	rax = u8ptr(rdi);	/* First biglit flag  */ \
	f0 = g->ZPAD0_6[0];/* Load remainder of divide by k */ \
	f0 += XMM_BIGVAL1; \
	single_rounding(base2, f0, f2, f4, rax); \
	f64ptr(rsi+0*64+32) = f0;	/* Save value1 */ \
\
	rax = u8ptr(rdi+4);	/* Load big vs. little flags */ \
	single_rounding(base2, f2, f0, f4, rax); \
	f2 *= f64ptr(rbx+1*32+16);	/* new value2 = val * col two-to-phi */ \
	f64ptr(rsi+1*64+32) = f2;	/* Save value2 */ \
\
	rcx = g->BIGLIT_INCR2;/* Different clm values step through */ \
;					/* big/lit array differently */ \
	rax = u8ptr(rdi+rcx);	/* Load big vs. little flags */ \
	single_rounding(base2, f0, f2, f4, rax); \
	f0 *= f64ptr(rbx+2*32+16);	/* new value3 = val * col two-to-phi */ \
	f64ptr(rsi+2*64+32) = f0;	/* Save value3 */ \
\
	f2 -= XMM_BIGVAL1;/* Remove integer rounding constant */ \
	f2 *= f64ptr(rbx+3*32+16);	/* value4 = carry * col two-to-phi */ \
	f64ptr(rsi+3*64+32) = f2;	/* Save new value4 */ \
\
;	/* Mul the integer part of (ZPAD data divided by k) by -c in */ \
;	/* preparation for adding it into the lower FFT data area. */ \
;	/* Also add in the shifted high FFT carry at this time. */ \
\
;	/* Now add in and normalize the bottom FFT data.  Remember that the */ \
;	/* column two-to-phi multiplier for the first value will be 1.0.  We  */ \
;	/* must go 6 words deep in case k is 48-50 bits and c is 32 bits. */ \
\
	vec2f64 TMP12 = g->u.xmm.XMM_TMP1; \
	vec2f64 TMP34 = g->u.xmm.XMM_TMP2; \
	vec2f64 TMP56 = g->u.xmm.XMM_TMP3; \
	TMP12 *= g->u.xmm.XMM_MINUS_C;	/* Mul by -c */ \
	TMP34 *= g->u.xmm.XMM_MINUS_C;	/* Mul by -c */ \
	TMP56 *= g->u.xmm.XMM_MINUS_C;	/* Mul by -c */ \
	rax = u8ptr(rdi);	/* First word  */ \
	f0 = TMP12[0];	/* Load integer part of divide by k */ \
	f0 += XMM_BIGVAL1; \
	f0 += f64ptr(rsi+0*64);	/* Add in the FFT data */ \
	single_rounding(base2, f0, f2, f4, rax); \
	f64ptr(rsi+0*64) = f0;	/* Save value1 */ \
\
	rax = u8ptr(rdi+4);	/* Load big vs. little flags */ \
	f0 = TMP12[1];	/* Load integer part of divide by k */ \
	f1 = f64ptr(rsi+1*64);	/* Load FFT data */ \
	f1 *= f64ptr(rbx+1*32);	/* Mul values2 by col two-to-minus-phi */ \
	f1 *= g->u.xmm.XMM_NORM012_FF[0];	/* Mul by FFTLEN/2 */ \
	f0 += f1;		/* Add in the FFT data */ \
	f0 += f2;		/* x2 = value + carry */ \
	single_rounding(base2, f0, f2, f4, rax); \
	f0 *= f64ptr(rbx+1*32+16);	/* new value2 = val * col two-to-phi */ \
	f64ptr(rsi+1*64) = f0;	/* Save value2 */ \
\
	rcx = rdi + g->BIGLIT_INCR2;/* Different clm values step through */ \
;					/* big/lit array differently */ \
	rax = u8ptr(rcx);	/* Load big vs. little flags */ \
	f0 = TMP34[0];	/* Load integer part of divide by k */ \
	f1 = f64ptr(rsi+2*64);	/* Load FFT data */ \
	f1 *= f64ptr(rbx+2*32);	/* Mul values3 by col two-to-minus-phi */ \
	f1 *= g->u.xmm.XMM_NORM012_FF[0];	/* Mul by FFTLEN/2 */ \
	f0 += f1;		/* Add in the FFT data */ \
	f0 += f2;		/* x3 = value + carry */ \
	single_rounding(base2, f0, f2, f4, rax); \
	f0 *= f64ptr(rbx+2*32+16);	/* new value3 = val * col two-to-phi */ \
	f64ptr(rsi+2*64) = f0;	/* Save value3 */ \
\
	rax = u8ptr(rcx+4); /* Load big vs. little flags */ \
	f0 = TMP34[1];	/* Load integer part of divide by k */ \
	f1 = f64ptr(rsi+3*64);	/* Load FFT data */ \
	f1 *= f64ptr(rbx+3*32);	/* Mul values4 by col two-to-minus-phi */ \
	f1 *= g->u.xmm.XMM_NORM012_FF[0];	/* Mul by FFTLEN/2 */ \
	f0 += f1;		/* Add in the FFT data */ \
	f0 += f2;		/* x4 = value + carry */ \
	single_rounding(base2, f0, f2, f4, rax); \
	f0 *= f64ptr(rbx+3*32+16);	/* new value4 = val * col two-to-phi */ \
	f64ptr(rsi+3*64) = f0;	/* Save value4 */ \
\
	rcx = rdi + g->BIGLIT_INCR4;/* Different clm values step through */ \
;					/* big/lit array differently */ \
	rax = u8ptr(rcx);	/* Load big vs. little flags */ \
	f0 = TMP56[0];	/* Load integer part of divide by k */ \
	f1 = f64ptr(rsi+4*64);	/* Load FFT data */ \
	f1 *= f64ptr(rbx+4*32);	/* Mul values4 by col two-to-minus-phi */ \
	f1 *= g->u.xmm.XMM_NORM012_FF[0];	/* Mul by FFTLEN/2 */ \
	f0 += f1;		/* Add in the FFT data */ \
	f0 += f2;		/* x4 = value + carry */ \
	single_rounding(base2, f0, f2, f4, rax); \
	f0 *= f64ptr(rbx+4*32+16);	/* new value4 = val * col two-to-phi */ \
	f64ptr(rsi+4*64) = f0;	/* Save value4 */ \
\
	rax = u8ptr(rcx+4); /* Load big vs. little flags */ \
	f0 = TMP56[1];	/* Load integer part of divide by k */ \
	f1 = f64ptr(rsi+5*64);	/* Load FFT data */ \
	f1 *= f64ptr(rbx+5*32);	/* Mul values5 by col two-to-minus-phi */ \
	f1 *= g->u.xmm.XMM_NORM012_FF[0];	/* Mul by FFTLEN/2 */ \
	f0 += f1;		/* Add in the FFT data */ \
	f0 += f2;		/* x5 = value + carry */ \
	single_rounding(base2, f0, f2, f4, rax); \
	f0 *= f64ptr(rbx+5*32+16);	/* new value5 = val * col two-to-phi */ \
	f64ptr(rsi+5*64) = f0;	/* Save value5 */ \
\
	f2 -= XMM_BIGVAL1;/* Remove rounding constant */ \
	f2 *= f64ptr(rbx+6*32+16);	/* new value6 = val * col two-to-phi */ \
	f2 += f64ptr(rsi+6*64);	/* Add in FFT data */ \
	f64ptr(rsi+6*64) = f2;	/* Save value6 */ \
}


/* *************** WPN normalized add/sub macro ****************** */ \
/* This macro adds or subtracts, then "normalizes" eight FFT */ \
/* data values.  This involves multiplying the summed values by */ \
/* two-to-minus-phi.  Rounding the value to an integer.  Making sure */ \
/* the integer is smaller than the maximum allowable integer, generating */ \
/* a carry if necessary. Finally, the value is multiplied by two-to-phi */ \
/* and stored. */ \
/* rcx = pointer to the first number */ \
/* rdx = pointer to the second number */ \
/* rsi = pointer to destination */ \
/* rdi = pointer to array of big vs. little flags */ \
/* rbp = pointer two-to-phi group multipliers */ \
/* XMM_TMP1,XMM_TMP2,XMM_TMP3,XMM_TMP4 = carries */ \
/* A pipelined version of this code: */ \
/*	xload	xmm0, [rdx]		;; Load second number */ \
/*	fop	xmm0, [rcx]		;; Add/sub first number */ \
/*	movzx	rax, BYTE PTR [rdi]	;; Load big vs. little flags */ \
/*	mulpd	xmm0, [rbp+0*XMM_GMD][rax] ;; Mul by fudged grp two-to-minus-phi */ \
/*	addpd	xmm0, XMM_TMP1		;; x1 = values + carry */ \
/*	xload	xmm2, XMM_LIMIT_BIGMAX[rax];; Load maximum * BIGVAL - BIGVAL */ \
/*	addpd	xmm2, xmm0		;; y1 = top bits of x */ \
/*	xload	xmm6, XMM_LIMIT_BIGMAX_NEG[rax];; Load -(maximum*BIGVAL-BIGVAL) */ \
/*	addpd	xmm6, xmm2		;; z1 = y1-(maximum * BIGVAL - BIGVAL) */ \
/*	subpd	xmm0, xmm6		;; rounded value = x1 - z1 */ \
/*	mulpd	xmm2, XMM_LIMIT_INVERSE[rax];; next carry = shifted y1 */ \
/*	mulpd	xmm0, [rbp+0*XMM_GMD+XMM_GMD/2][rax] ;; new value1 = val * fudged grp two-to-phi */ \
/*	xstore	[rsi+0*dist1], xmm0	;; Save new value1 */ \
/*	xstore	XMM_TMP1, xmm2		;; Save carry */ \

#define xnorm_op_wpn(fop, ttp, base2, TMP1, TMP2, TMP3, TMP4) \
	xmm0 = xptr(rdx);		/* Load second number */ \
	fop(xmm0, xptr(rcx));		/* Add/sub first number */ \
	xmm1 = xptr(rdx+16);		/* Load second number */ \
	fop(xmm1, xptr(rcx+16));		/* Add/sub first number */ \
	xmm6 = xptr(rdx+32);		/* Load second number */ \
	fop(xmm6, xptr(rcx+32));		/* Add/sub first number */ \
	xmm7 = xptr(rdx+48);		/* Load second number */ \
	fop(xmm7, xptr(rcx+48));		/* Add/sub first number */ \
ttp(rax = u8ptr(rdi+1));	/* Load big vs. little flags 1-4 */ \
ttp(rbx = u8ptr(rdi));	/* Load fudge factor flags 1,2 */ \
ttp(rbx &= 0xF0); \
ttp(xmm0 *= xptr(rbp+0*XMM_GMD+rbx)); /* Mul by fudged grp two-to-minus-phi */ \
	xmm0 += TMP1;		/* x1 = values + carry */ \
ttp(xmm1 *= xptr(rbp+1*XMM_GMD+rbx)); /* Mul by fudged grp two-to-minus-phi */ \
	xmm1 += TMP2;		/* x1 = values + carry */ \
	rounding_interleaved(ttp, base2, noexec, noexec, xmm0, xmm2, xmm4, rax*4, xmm1, xmm3, xmm5, rax*4+16); \
	TMP1 = xmm2;		/* Save carry */ \
	TMP2 = xmm3;		/* Save carry */ \
ttp(xmm0 *= xptr(rbp+0*XMM_GMD+XMM_GMD/2+rbx)); /* value1 *= fudged grp two-to-phi */ \
	xptr(rsi+0*16) = xmm0;	/* Save new value1 */ \
ttp(xmm1 *= xptr(rbp+1*XMM_GMD+XMM_GMD/2+rbx)); /* value2 *= fudged grp two-to-phi */ \
	xptr(rsi+1*16) = xmm1;	/* Save new value1 */ \
\
ttp(rbx = u8ptr(rdi));	/* Load fudge factor flags 3,4 */ \
ttp(rbx &= 0x0F); \
ttp(xmm6 *= xptr(rbp+2*XMM_GMD+rbx*8)); /* Mul by fudged grp two-to-minus-phi */ \
	xmm6 += TMP3;		/* x1 = values + carry */ \
ttp(xmm7 *= xptr(rbp+3*XMM_GMD+rbx*8)); /* Mul by fudged grp two-to-minus-phi */ \
	xmm7 += TMP4;		/* x1 = values + carry */ \
	rounding_interleaved(ttp, base2, noexec, noexec, xmm6, xmm2, xmm4, rax*4+32, xmm7, xmm3, xmm5, rax*4+48); \
	TMP3 = xmm2;		/* Save carry */ \
	TMP4 = xmm3;		/* Save carry */ \
ttp(xmm6 *= xptr(rbp+2*XMM_GMD+XMM_GMD/2+rbx*8)); /* value1 *= fudged grp two-to-phi */ \
	xptr(rsi+2*16) = xmm6;	/* Save new value1 */ \
ttp(xmm7 *= xptr(rbp+3*XMM_GMD+XMM_GMD/2+rbx*8)); /* value2 *= fudged grp two-to-phi */ \
	xptr(rsi+3*16) = xmm7;	/* Save new value1 */ \
ttp(rdi += 2);			/* Next flags ptr */ \
	rsi += 64;			/* Next dest ptr */ \
	rdx += 64;			/* Next src ptr */ \
	rcx += 64;			/* Next src ptr */ \
\
\
/* *************** WPN followup macros ****************** */ \
/* This macro finishes the normalize add/sub process by adding four carries */ \
/* from the end of a block back to the start of the block.  The remaining */ \
/* four carries are rotated for starting the next block. */ \
\

#define xnorm_op_wpn_blk(srcreg, screg, carry1, carry2, carry3, carry4) \
	xmm2 = XMM_BIGVAL2; \
\
	rax = u8ptr(rdi);	/* Load fudge factor flags 1,2 */ \
	rax &= 0xF0; \
	xmm0 = carry1;		/* Load carries */ \
	xmm0[0] -= XMM_BIGVAL1; \
	xmm0[0] *= f64ptr(screg+0*XMM_GMD+XMM_GMD/2+8+rax); /* mul by fudged grp two-to-phi */ \
	xmm0[0] += f64ptr(srcreg+0*16+8);	/* Add in FFT data */ \
	f64ptr(srcreg+0*16+8) = xmm0[0];	/* Save FFT word */ \
	xmm0[0] = xmm0[1]; xmm0[1] = XMM_BIGVAL2[1]; /*shufpd xmm0, xmm2, 1;*/ \
	carry1 = xmm0;		/* Save carries for start of next blk */ \
\
	xmm0 = carry2;		/* Load carries */ \
	xmm0[0] -= XMM_BIGVAL1; \
	xmm0[0] *= f64ptr(screg+1*XMM_GMD+XMM_GMD/2+8+rax); /* mul by fudged grp two-to-phi */ \
	xmm0[0] += f64ptr(srcreg+1*16+8);	/* Add in FFT data */ \
	f64ptr(srcreg+1*16+8) = xmm0[0];	/* Save FFT word */ \
	xmm0[0] = xmm0[1]; xmm0[1] = XMM_BIGVAL2[1]; \
	carry2 = xmm0;		/* Save carries for start of next blk */ \
\
	rax = u8ptr(rdi);	/* Load fudge factor flags 3,4 */ \
	rax &= 0x0F; \
	xmm0 = carry3;		/* Load carries */ \
	xmm0[0] -= XMM_BIGVAL1; \
	xmm0[0] *= f64ptr(screg+2*XMM_GMD+XMM_GMD/2+8+rax*8); /* mul by fudged grp two-to-phi */ \
	xmm0[0] += f64ptr(srcreg+2*16+8);	/* Add in FFT data */ \
	f64ptr(srcreg+2*16+8) = xmm0[0];	/* Save FFT word */ \
	xmm0[0] = xmm0[1]; xmm0[1] = XMM_BIGVAL2[1]; \
	carry3 = xmm0;		/* Save carries for start of next blk */ \
\
	xmm0 = carry4;		/* Load carries */ \
	xmm0[0] -= XMM_BIGVAL1; \
	xmm0[0] *= f64ptr(screg+3*XMM_GMD+XMM_GMD/2+8+rax*8); /* mul by fudged grp two-to-phi */ \
	xmm0[0] += f64ptr(srcreg+3*16+8);	/* Add in FFT data */ \
	f64ptr(srcreg+3*16+8) = xmm0[0];	/* Save FFT word */ \
	xmm0[0] = xmm0[1]; xmm0[1] = XMM_BIGVAL2[1]; \
	carry4 = xmm0;		/* Save carries for start of next blk */ \
\
\
/* This macro finishes the normalize add/sub process by adding two carries */ \
/* from the end of a section back to the start of the section.  The remaining */ \
/* two carries are rotated for starting the next section. */ \
/* rax = pointer to FFT data */ \
/* rbx = pointer two-to-phi group multipliers */ \
\

#define xnorm_op_wpn_sec(carry1, carry2, carry3, carry4) \
	xmm2 = XMM_BIGVAL2; \
\
	rcx = u8ptr(rdi);	/* Load fudge factor flags 1,2 */ \
	rcx &= 0xF0; \
	xmm0[0] = carry1[0];		/* Load carry */ \
	xmm0[0] -= XMM_BIGVAL1; \
	xmm0[0] *= f64ptr(rbx+1*XMM_GMD+XMM_GMD/2+rcx); /* mul by fudged grp two-to-phi */ \
	xmm0[0] += f64ptr(rax+1*16);	/* Add in FFT data */ \
	f64ptr(rax+1*16) = xmm0[0];	/* Save FFT word */ \
\
	xmm1 = carry2; \
	carry1 = xmm1; \
	carry2 = xmm2;		/* Save carries for start of next sec */ \
\
	rcx = u8ptr(rdi);	/* Load fudge factor flags 3,4 */ \
	rcx &= 0x0F; \
	xmm0[0] = carry3[0];		/* Load carry */ \
	xmm0[0] -= XMM_BIGVAL1; \
	xmm0[0] *= f64ptr(rbx+3*XMM_GMD+XMM_GMD/2+rcx*8); /* mul by fudged grp two-to-phi */ \
	xmm0[0] += f64ptr(rax+3*16);	/* Add in FFT data */ \
	f64ptr(rax+3*16) = xmm0[0];	/* Save FFT word */ \
\
	xmm1 = carry4; \
	carry3 = xmm1; \
	carry4 = xmm2;		/* Save carries for start of next sec */ \
\
\
\
/* This macro finishes the normalize add/sub process by adding the final */ \
/* two carries back into the appropriate FFT values at the start of the fft. */ \
/* rsi = pointer to FFT data */ \
/* rbp = pointer two-to-phi group multipliers */ \
/* xmm6 = non-wraparound carry */ \
/* xmm7 = wraparound carry */ \
\

#define xnorm_op_wpn_fft \
	xmm2[0] = XMM_BIGVAL1; xmm2[1] = 0.0; \
\
	rax = u8ptr(rdi);	/* Load fudge factor flags 1,2 */ \
	rax &= 0xF0; \
	xmm7[0] -= XMM_BIGVAL1;		/* Remove XMM_BIGVAL */ \
	xmm7[0] *= g->u.xmm.XMM_MINUS_C[0];	/* mul wrap around carry by -c */ \
	xmm7[0] *= f64ptr(rbp+0*XMM_GMD+XMM_GMD/2+rax); /* mul by fudged grp two-to-phi */ \
	xmm7[0] += f64ptr(rsi+0*16);	/* Add in FFT data */ \
	f64ptr(rsi+0*16) = xmm7[0];	/* Save FFT word */ \
\
	rax = u8ptr(rdi);	/* Load fudge factor flags 3,4 */ \
	rax &= 0x0F; \
	xmm6[0] -= XMM_BIGVAL1;		/* Remove XMM_BIGVAL */ \
	xmm6[0] *= f64ptr(rbp+2*XMM_GMD+XMM_GMD/2+rax*8); /* mul by fudged grp two-to-phi */ \
	xmm6[0] += f64ptr(rsi+2*16);	/* Add in FFT data */ \
	f64ptr(rsi+2*16) = xmm6[0];	/* Save FFT word */ \
\
\
\
/* *************** WPN normalized add & sub macro ****************** */ \
/* This macro adds and subtracts, then "normalizes" eight FFT */ \
/* data values.  This involves multiplying the summed values by */ \
/* two-to-minus-phi.  Rounding the value to an integer.  Making sure */ \
/* the integer is smaller than the maximum allowable integer, generating */ \
/* a carry if necessary. Finally, the value is multiplied by two-to-phi */ \
/* and stored. */ \
/* rcx = pointer to the first number */ \
/* rdx = pointer to the second number */ \
/* rsi = pointer to destination #1 */ \
/* rbp = pointer to destination #2 */ \
/* rbx = pointer two-to-phi group multipliers */ \
/* rdi = pointer to array of big vs. little flags */ \
/* rax = big/lit flag */ \
/* XMM_TMP1,XMM_TMP2,XMM_TMP3,XMM_TMP4 = addition carries */ \
/* XMM_TMP5,XMM_TMP6,XMM_TMP7,XMM_TMP8 = subtraction carries */ \
\

#define xnorm_addsub_wpn(ttp, base2, TMP1, TMP2, TMP3, TMP4, TMP5, TMP6, TMP7, TMP8) \
	xmm1 = xptr(rcx+0*16);	/* Load first number */ \
	xmm0 = xptr(rdx+0*16);	/* Load second number */ \
	xmm1 -= xmm0;		/* first - second number */ \
	xmm0 += xptr(rcx+0*16);	/* first + second number */ \
\
	xmm7 = xptr(rcx+1*16);	/* Load first number */ \
	xmm6 = xptr(rdx+1*16);	/* Load second number */ \
	xmm7 -= xmm6;		/* first - second number */ \
	xmm6 += xptr(rcx+1*16);	/* first + second number */ \
\
ttp(rax = u8ptr(rdi));	/* Load fudge factor flags 1,2 */ \
ttp(rax &= 0xF0); \
ttp(xmm0 *= xptr(rbx+0*XMM_GMD+rax)); /* Mul by fudged grp two-to-minus-phi */ \
	xmm0 += TMP1;		/* x1 = values + carry */ \
ttp(xmm1 *= xptr(rbx+0*XMM_GMD+rax)); /* Mul by fudged grp two-to-minus-phi */ \
	xmm1 += TMP5;		/* x1 = values + carry */ \
ttp(xmm6 *= xptr(rbx+1*XMM_GMD+rax)); /* Mul by fudged grp two-to-minus-phi */ \
	xmm6 += TMP2;		/* x1 = values + carry */ \
ttp(xmm7 *= xptr(rbx+1*XMM_GMD+rax)); /* Mul by fudged grp two-to-minus-phi */ \
	xmm7 += TMP6;		/* x1 = values + carry */ \
\
ttp(rax = u8ptr(rdi+1));	/* Load big vs. little flags 1-4 */ \
	rounding_interleaved(ttp, base2, noexec, noexec, xmm0, xmm2, xmm4, rax*4, xmm1, xmm3, xmm5, rax*4); \
	TMP1 = xmm2;		/* Save carry */ \
	TMP5 = xmm3;		/* Save carry */ \
	rounding_interleaved(ttp, base2, noexec, noexec, xmm6, xmm2, xmm4, rax*4+16, xmm7, xmm3, xmm5, rax*4+16); \
	TMP2 = xmm2;		/* Save carry */ \
	TMP6 = xmm3;		/* Save carry */ \
\
ttp(rax = u8ptr(rdi));	/* Load fudge factor flags 1,2 */ \
ttp(rax &= 0xF0); \
ttp(xmm0 *= xptr(rbx+0*XMM_GMD+XMM_GMD/2+rax)); /* value1 *= fudged grp two-to-phi */ \
	xptr(rsi+0*16) = xmm0;	/* Save new value1 */ \
ttp(xmm1 *= xptr(rbx+0*XMM_GMD+XMM_GMD/2+rax)); /* value2 *= fudged grp two-to-phi */ \
	xptr(rbp+0*16) = xmm1;	/* Save new value1 */ \
ttp(xmm6 *= xptr(rbx+1*XMM_GMD+XMM_GMD/2+rax)); /* value1 *= fudged grp two-to-phi */ \
	xptr(rsi+1*16) = xmm6;	/* Save new value1 */ \
ttp(xmm7 *= xptr(rbx+1*XMM_GMD+XMM_GMD/2+rax)); /* value2 *= fudged grp two-to-phi */ \
	xptr(rbp+1*16) = xmm7;	/* Save new value1 */ \
\
	xmm1 = xptr(rcx+2*16);	/* Load first number */ \
	xmm0 = xptr(rdx+2*16);	/* Load second number */ \
	xmm1 -= xmm0;		/* first - second number */ \
	xmm0 += xptr(rcx+2*16);	/* first + second number */ \
\
	xmm7 = xptr(rcx+3*16);	/* Load first number */ \
	xmm6 = xptr(rdx+3*16);	/* Load second number */ \
	xmm7 -= xmm6;		/* first - second number */ \
	xmm6 += xptr(rcx+3*16);	/* first + second number */ \
\
ttp(rax = u8ptr(rdi));	/* Load fudge factor flags 3,4 */ \
ttp(rax &= 0x0F); \
ttp(xmm0 *= xptr(rbx+2*XMM_GMD+rax*8)); /* Mul by fudged grp two-to-minus-phi */ \
	xmm0 += TMP3;		/* x1 = values + carry */ \
ttp(xmm1 *= xptr(rbx+2*XMM_GMD+rax*8)); /* Mul by fudged grp two-to-minus-phi */ \
	xmm1 += TMP7;		/* x1 = values + carry */ \
ttp(xmm6 *= xptr(rbx+3*XMM_GMD+rax*8)); /* Mul by fudged grp two-to-minus-phi */ \
	xmm6 += TMP4;		/* x1 = values + carry */ \
ttp(xmm7 *= xptr(rbx+3*XMM_GMD+rax*8)); /* Mul by fudged grp two-to-minus-phi */ \
	xmm7 += TMP8;		/* x1 = values + carry */ \
\
ttp(rax = u8ptr(rdi+1));	/* Load big vs. little flags 1-4 */ \
	rounding_interleaved(ttp, base2, noexec, noexec, xmm0, xmm2, xmm4, rax*4+32, xmm1, xmm3, xmm5, rax*4+32); \
	TMP3 = xmm2;		/* Save carry */ \
	TMP7 = xmm3;		/* Save carry */ \
	rounding_interleaved(ttp, base2, noexec, noexec, xmm6, xmm2, xmm4, rax*4+48, xmm7, xmm3, xmm5, rax*4+48); \
	TMP4 = xmm2;		/* Save carry */ \
	TMP8 = xmm3;		/* Save carry */ \
\
ttp(rax = u8ptr(rdi));	/* Load fudge factor flags 3,4 */ \
ttp(rax &= 0x0F); \
ttp(xmm0 *= xptr(rbx+2*XMM_GMD+XMM_GMD/2+rax*8)); /* value1 *= fudged grp two-to-phi */ \
	xptr(rsi+2*16) = xmm0;	/* Save new value1 */ \
ttp(xmm1 *= xptr(rbx+2*XMM_GMD+XMM_GMD/2+rax*8)); /* value2 *= fudged grp two-to-phi */ \
	xptr(rbp+2*16) = xmm1;	/* Save new value1 */ \
ttp(xmm6 *= xptr(rbx+3*XMM_GMD+XMM_GMD/2+rax*8)); /* value1 *= fudged grp two-to-phi */ \
	xptr(rsi+3*16) = xmm6;	/* Save new value1 */ \
ttp(xmm7 *= xptr(rbx+3*XMM_GMD+XMM_GMD/2+rax*8)); /* value2 *= fudged grp two-to-phi */ \
	xptr(rbp+3*16) = xmm7;	/* Save new value1 */ \
\
	rdx += 64;			/* Next src ptr */ \
	rcx += 64;			/* Next src ptr */ \
ttp(rdi += 2);			/* Next flags ptr */ \
	rsi += 64;			/* Next dest ptr */ \
	rbp += 64;			/* Next dest ptr */ \
\
\
/* *************** WPN normalized small mul macro ****************** */ \
/* This macro multiplies by a small value, then "normalizes" eight FFT */ \
/* data values.  */ \
/* rsi = pointer to destination */ \
/* rdi = pointer to array of big vs. little flags */ \
/* rbp = pointer two-to-phi group multipliers */ \
/* XMM_TMP1,XMM_TMP2,XMM_TMP3,XMM_TMP4 = carries */ \
/* XMM_TMP5 = small value * optional FFTLEN/2 */ \
/* A pipelined version of this code: */ \
/*	xload	xmm0, [rsi]		;; Load second number */ \
/*	mulpd	xmm0, XMM_TMP5		;; Mul by small value */ \
/*	movzx	rcx, BYTE PTR [rdi]	;; Load big vs. little flags */ \
/*	xload	xmm2, XMM_TTMP_FUDGE[rcx];; fudge two-to-minus-phi */ \
/*	mulpd	xmm0, [rax]		;; Mul by grp two-to-minus-phi */ \
/*	mulpd	xmm0, xmm2		;; Mul by fudge two-to-minus-phi */ \
/*	addpd	xmm0, [rbp+0*16]	;; x1 = values + carry */ \
/*	xload	xmm2, XMM_LIMIT_BIGMAX[rcx];; Load maximum * BIGVAL - BIGVAL */ \
/*	addpd	xmm2, xmm0		;; y1 = top bits of x */ \
/*	xload	xmm6, XMM_LIMIT_BIGMAX_NEG[rcx];; Load -(maximum*BIGVAL-BIGVAL) */ \
/*	addpd	xmm6, xmm2		;; z1 = y1-(maximum * BIGVAL - BIGVAL) */ \
/*	subpd	xmm0, xmm6		;; rounded value = x1 - z1 */ \
/*	mulpd	xmm2, XMM_LIMIT_INVERSE[rcx];; next carry = shifted y1 */ \
/*	xload	xmm4, XMM_TTP_FUDGE[rcx];; fudge two-to-phi */ \
/*	mulpd	xmm0, [rax+0*32+16]	;; new value1 = val * grp two-to-phi */ \
/*	mulpd	xmm0, xmm4		;; new value1 *= fudge two-to-phi */ \
/*	xstore	[rsi], xmm0		;; Save new value1 */ \
/*	xstore	[rbp+0*16], xmm2	;; Save carry */ \
\

#define xnorm_smallmul_wpn(ttp, base2, TMP1, TMP2, TMP3, TMP4, TMP5) \
	xmm2 = TMP5;		/* Load small value */ \
	xmm0 = xptr(rsi);		/* Load values1 */ \
	xmm0 *= xmm2;		/* Mul values1 by small value */ \
	xmm1 = xptr(rsi+16);		/* Load values2 */ \
	xmm1 *= xmm2;		/* Mul values2 by small value */ \
	xmm6 = xptr(rsi+32);		/* Load values3 */ \
	xmm6 *= xmm2;		/* Mul values3 by small value */ \
	xmm7 = xptr(rsi+48);		/* Load values4 */ \
	xmm7 *= xmm2;		/* Mul values4 by small value */ \
ttp(rax = u8ptr(rdi+1));	/* Load big vs. little flags 1-4 */ \
ttp(rcx = u8ptr(rdi));	/* Load fudge factor flags 1,2 */ \
ttp(rcx &= 0xF0); \
ttp(xmm0 *= xptr(rbp+0*XMM_GMD+rcx)); /* Mul by fudged grp two-to-minus-phi */ \
	xmm0 += TMP1;		/* x1 = values + carry */ \
ttp(xmm1 *= xptr(rbp+1*XMM_GMD+rcx)); /* Mul by fudged grp two-to-minus-phi */ \
	xmm1 += TMP2;		/* x1 = values + carry */ \
	rounding_interleaved(ttp, base2, noexec, noexec, xmm0, xmm2, xmm4, rax*4, xmm1, xmm3, xmm5, rax*4+16); \
	TMP1 = xmm2;		/* Save carry */ \
	TMP2 = xmm3;		/* Save carry */ \
ttp(xmm0 *= xptr(rbp+0*XMM_GMD+XMM_GMD/2+rcx)); /* value1 *= fudged two-to-phi */ \
	xptr(rsi+0*16) = xmm0;	/* Save new value1 */ \
ttp(xmm1 *= xptr(rbp+1*XMM_GMD+XMM_GMD/2+rcx)); /* value2 *= fudged two-to-phi */ \
	xptr(rsi+1*16) = xmm1;	/* Save new value2 */ \
\
ttp(rcx = u8ptr(rdi));	/* Load fudge factor flags 3,4 */ \
ttp(rcx &= 0x0F); \
ttp(xmm6 *= xptr(rbp+2*XMM_GMD+rcx*8)); /* Mul by fudged grp two-to-minus-phi */ \
	xmm6 += TMP3;		/* x1 = values + carry */ \
ttp(xmm7 *= xptr(rbp+3*XMM_GMD+rcx*8)); /* Mul by fudged grp two-to-minus-phi */ \
	xmm7 += TMP4;		/* x1 = values + carry */ \
	rounding_interleaved(ttp, base2, noexec, noexec, xmm6, xmm2, xmm4, rax*4+32, xmm7, xmm3, xmm5, rax*4+48); \
	TMP3 = xmm2;		/* Save carry */ \
	TMP4 = xmm3;		/* Save carry */ \
ttp(xmm6 *= xptr(rbp+2*XMM_GMD+XMM_GMD/2+rcx*8)); /* value1 *= fudged two-to-phi */ \
	xptr(rsi+2*16) = xmm6;	/* Save new value1 */ \
ttp(xmm7 *= xptr(rbp+3*XMM_GMD+XMM_GMD/2+rcx*8)); /* value2 *= fudged two-to-phi */ \
	xptr(rsi+3*16) = xmm7;	/* Save new value2 */ \
ttp(rdi += 2);			/* Next flags ptr */ \
	rsi += 64;			/* Next dest ptr */ \
\
\
/* This macro finishes the smallmul normalize process by adding four carries */ \
/* from the end of a block back to the start of the block.  The remaining */ \
/* four carries are rotated for starting the next block. */ \
/* rsi = Source ptr to start of block */ \
/* rbp = Group ptr of two-to-phi */ \
/* rdi = biglit array ptr */ \
/* rcx is destroyed */ \
\

#define xnorm_smallmul_wpn_blk(base2, TMP1, TMP2, TMP3, TMP4) \
	xmm5[0] = XMM_BIGVAL1; xmm5[1] = 0.0; \
\
	xmm1[0] = TMP1[0];	/* Load carry */ \
	xnorm_smallmul_wpn_prop4(base2, rsi+0*16, 8, rbp+0*XMM_GMD, rdi, 0, xmm1[0]); \
	xmm1 = TMP1;		/* Load carries (prop4 will zero high word of xmm1) */ \
	xmm1[0] = xmm1[1]; xmm1[1] = XMM_BIGVAL2[1]; /*shufpd	xmm1, xmm5, 1;*/ \
	TMP1 = xmm1;		/* Save carries for start of next blk */ \
\
	xmm1[0] = TMP2[0];	/* Load carry */ \
	xnorm_smallmul_wpn_prop4(base2, rsi+1*16, 8, rbp+1*XMM_GMD, rdi, 1, xmm1[0]); \
	xmm1 = TMP2;		/* Load carries */ \
	xmm1[0] = xmm1[1]; xmm1[1] = XMM_BIGVAL2[1]; \
	TMP2 = xmm1;		/* Save carries for start of next blk */ \
\
	xmm1[0] = TMP3[0];	/* Load carry */ \
	xnorm_smallmul_wpn_prop4(base2, rsi+2*16, 8, rbp+2*XMM_GMD, rdi, 2, xmm1[0]); \
	xmm1 = TMP3;		/* Load carries */ \
	xmm1[0] = xmm1[1]; xmm1[1] = XMM_BIGVAL2[1]; \
	TMP3 = xmm1;		/* Save carries for start of next blk */ \
\
	xmm1[0] = TMP4[0];	/* Load carry */ \
	xnorm_smallmul_wpn_prop4(base2, rsi+3*16, 8, rbp+3*XMM_GMD, rdi, 3, xmm1[0]); \
	xmm1 = TMP4;		/* Load carries */ \
	xmm1[0] = xmm1[1]; xmm1[1] = XMM_BIGVAL2[1]; \
	TMP4 = xmm1;		/* Save carries for start of next blk */ \
\
\
/* This macro finishes the smallmul normalize process by adding two carries */ \
/* from the end of a section back to the start of the section.  The remaining */ \
/* two carries are rotated for starting the next section. */ \
/* rsi = Source ptr to start of block */ \
/* rbp = Group ptr of two-to-phi */ \
/* rdi = biglit array ptr */ \
/* rcx is destroyed */ \
\

#define xnorm_smallmul_wpn_sec(base2, TMP1, TMP2, TMP3, TMP4) \
	xmm5 = XMM_BIGVAL2; \
\
	xmm1[0] = TMP1[0];	/* Load carry */ \
	xnorm_smallmul_wpn_prop4(base2, rsi+1*16, 0, rbp+1*XMM_GMD, rdi, 1, xmm1[0]); \
	xmm1 = TMP2; \
	TMP1 = xmm1; \
	TMP2 = xmm5;		/* Save carries for start of next sec */ \
\
	xmm1[0] = TMP3[0];	/* Load carry */ \
	xnorm_smallmul_wpn_prop4(base2, rsi+3*16, 0, rbp+3*XMM_GMD, rdi, 3, xmm1[0]); \
	xmm1 = TMP4; \
	TMP3 = xmm1; \
	TMP4 = xmm5;		/* Save carries for start of next sec */ \
\
\
/* This macro finishes the smallmul normalize process by adding the final */ \
/* two carries back into the appropriate FFT values at the start of the fft. */ \
/* rsi = pointer to FFT data */ \
/* rbp = pointer two-to-phi group multipliers */ \
/* rdi = pointer to big/little flags */ \
/* xmm6 = non-wraparound carry */ \
/* xmm7 = wraparound carry */ \
\

#define xnorm_smallmul_wpn_fft(base2) \
	xmm7[0] -= XMM_BIGVAL1; \
	xmm7[0] *= g->u.xmm.XMM_MINUS_C[0];	/* Negate the carry */ \
	xmm7[0] += XMM_BIGVAL1; \
	xnorm_smallmul_wpn_prop6(base2, rsi+0*16, 0, rbp+0*XMM_GMD, rdi, 0, xmm7[0]); \
	xnorm_smallmul_wpn_prop4(base2, rsi+2*16, 0, rbp+2*XMM_GMD, rdi, 2, xmm6[0]); \
\
\
/* Propagate a single gwsmallmul carry across 4 words. */ \
\

#define BIGLITWORD_LE_1_0(A,B) A
#define BIGLITWORD_LE_1_1(A,B) A
#define BIGLITWORD_LE_1_2(A,B) B
#define BIGLITWORD_LE_1_3(A,B) B

#define xnorm_smallmul_wpn_prop4(base2, srcptr, off8, grpptr, biglitreg, biglitword, xcarry) { \
	double f0, xtmp; \
	rax = u8ptr(biglitreg+1);	/* Big/lit flags */ \
	rcx = u8ptr(biglitreg);	/* Fudge factor flags */ \
	rcx &= BIGLITWORD_LE_1_##biglitword(0xF0, 0x0F); \
	f0 = f64ptr(srcptr+0*64+off8);	/* Load FFT data */ \
	f0 *= f64ptr(grpptr+off8+BIGLITWORD_LE_1_##biglitword(rcx, rcx*8));	/* Mul by fudged two-to-minus-phi */ \
	f0 += xcarry;			/* Add in the carry */ \
	single_rounding(base2, f0, xcarry, xtmp, rax*4+biglitword*16+off8); \
	f0 *= f64ptr(grpptr+XMM_GMD/2+off8+BIGLITWORD_LE_1_##biglitword(rcx, rcx*8)); /* value1 *= fudged two-to-phi */ \
	f64ptr(srcptr+0*64+off8) = f0;	/* Save value1 */ \
\
	rax = u8ptr(biglitreg+3);	/* Second word big/lit flags */ \
	rcx = u8ptr(biglitreg+2);	/* Second word fudge factor flags */ \
	rcx &= BIGLITWORD_LE_1_##biglitword(0xF0, 0x0F); \
	f0 = f64ptr(srcptr+1*64+off8);	/* Load FFT data */ \
	f0 *= f64ptr(grpptr+off8+BIGLITWORD_LE_1_##biglitword(rcx, rcx*8));	/* Mul by fudged two-to-minus-phi */ \
	f0 += xcarry;			/* Add in the carry */ \
	single_rounding(base2, f0, xcarry, xtmp, rax*4+biglitword*16+off8); \
	f0 *= f64ptr(grpptr+XMM_GMD/2+off8+BIGLITWORD_LE_1_##biglitword(rcx, rcx*8)); /* value2 *= fudged two-to-phi */ \
	f64ptr(srcptr+1*64+off8) = f0;	/* Save value2 */ \
\
	rcx = g->BIGLIT_INCR2;/* Different clm values step through */ \
;						/* big/lit array differently */ \
	rax = u8ptr(biglitreg+1+rcx); /* Third word big/lit flags */ \
	rcx = u8ptr(biglitreg+rcx);	/* Third word fudge factor flags */ \
	rcx &= BIGLITWORD_LE_1_##biglitword(0xF0, 0x0F); \
	f0 = f64ptr(srcptr+2*64+off8);	/* Load FFT data */ \
	f0 *= f64ptr(grpptr+off8+BIGLITWORD_LE_1_##biglitword(rcx, rcx*8));	/* Mul by fudged two-to-minus-phi */ \
	f0 += xcarry;			/* Add in the carry */ \
	single_rounding(base2, f0, xcarry, xtmp, rax*4+biglitword*16+off8); \
	f0 *= f64ptr(grpptr+XMM_GMD/2+off8+BIGLITWORD_LE_1_##biglitword(rcx, rcx*8)); /* value3 *= fudged two-to-phi */ \
	f64ptr(srcptr+2*64+off8) = f0;	/* Save value3 */ \
\
	rcx = g->BIGLIT_INCR2;/* Different clm values step through */ \
;						/* big/lit array differently */ \
	rcx = u8ptr(biglitreg+rcx+2); /* Fourth word fudge factor flags */ \
	rcx &= BIGLITWORD_LE_1_##biglitword(0xF0, 0x0F); \
	xcarry -= XMM_BIGVAL1;/* Make carry an integer */ \
	xcarry *= f64ptr(grpptr+XMM_GMD/2+off8+BIGLITWORD_LE_1_##biglitword(rcx, rcx*8)); /* carry *= fudged two-to-phi */ \
	xcarry += f64ptr(srcptr+3*64+off8);	/* Add FFT data */ \
	f64ptr(srcptr+3*64+off8) = xcarry;	/* Save value4 */ \
}

/* Propagate a single gwsmallmul carry across 6 words. */ \


#define biglitword_only_support_00

#define xnorm_smallmul_wpn_prop6(base2, srcptr, off8, grpptr, biglitreg, biglitword, xcarry) { \
	double f0, xtmp; \
	biglitword_only_support_0##biglitword; \
\
	rax = u8ptr(biglitreg+1);	/* Big/lit flags */ \
	rcx = u8ptr(biglitreg);	/* Fudge factor flags */ \
	rcx &= 0xF0; \
	f0 = f64ptr(srcptr+0*64+off8);	/* Load FFT data */ \
	f0 *= f64ptr(grpptr+off8+rcx);	/* Mul by fudged two-to-minus-phi */ \
	f0 += xcarry;			/* Add in the carry */ \
	single_rounding(base2, f0, xcarry, xtmp, rax*4+off8); \
	f0 *= f64ptr(grpptr+XMM_GMD/2+off8+rcx); /* value1 *= fudged two-to-phi */ \
	f64ptr(srcptr+0*64+off8) = f0;	/* Save value1 */ \
\
	rax = u8ptr(biglitreg+3);	/* Big/lit flags */ \
	rcx = u8ptr(biglitreg+2);	/* Fudge factor flags */ \
	rcx &= 0xF0; \
	f0 = f64ptr(srcptr+1*64+off8);	/* Load FFT data */ \
	f0 *= f64ptr(grpptr+off8+rcx);	/* Mul by fudged two-to-minus-phi */ \
	f0 += xcarry;			/* Add in the carry */ \
	single_rounding(base2, f0, xcarry, xtmp, rax*4+off8); \
	f0 *= f64ptr(grpptr+XMM_GMD/2+off8+rcx); /* value2 *= fudged two-to-phi */ \
	f64ptr(srcptr+1*64+off8) = f0;	/* Save value2 */ \
\
	rcx = g->BIGLIT_INCR2;/* Different clm values step through */ \
;						/* big/lit array differently */ \
	rax = u8ptr(biglitreg+1+rcx); /* Big/lit flags */ \
	rcx = u8ptr(biglitreg+rcx);	/* Fudge factor flags */ \
	rcx &= 0xF0; \
	f0 = f64ptr(srcptr+2*64+off8);	/* Load FFT data */ \
	f0 *= f64ptr(grpptr+off8+rcx);	/* Mul by fudged two-to-minus-phi */ \
	f0 += xcarry;			/* Add in the carry */ \
	single_rounding(base2, f0, xcarry, xtmp, rax*4+off8); \
	f0 *= f64ptr(grpptr+XMM_GMD/2+off8+rcx); /* value3 *= fudged two-to-phi */ \
	f64ptr(srcptr+2*64+off8) = f0;	/* Save value3 */ \
\
	rcx = g->BIGLIT_INCR2;/* Different clm values step through */ \
;						/* big/lit array differently */ \
	rax = u8ptr(biglitreg+3+rcx); /* Big/lit flags */ \
	rcx = u8ptr(biglitreg+2+rcx); /* Fudge factor flags */ \
	rcx &= 0xF0; \
	f0 = f64ptr(srcptr+3*64+off8);	/* Load FFT data */ \
	f0 *= f64ptr(grpptr+off8+rcx);	/* Mul by fudged two-to-minus-phi */ \
	f0 += xcarry;			/* Add in the carry */ \
	single_rounding(base2, f0, xcarry, xtmp, rax*4+off8); \
	f0 *= f64ptr(grpptr+XMM_GMD/2+off8+rcx); /* value4 *= fudged two-to-phi */ \
	f64ptr(srcptr+3*64+off8) = f0;	/* Save value4 */ \
\
	rcx = g->BIGLIT_INCR4;/* Different clm values step through */ \
;						/* big/lit array differently */ \
	rax = u8ptr(biglitreg+1+rcx); /* Big/lit flags */ \
	rcx = u8ptr(biglitreg+rcx);	/* Fudge factor flags */ \
	rcx &= 0xF0; \
	f0 = f64ptr(srcptr+4*64+off8);	/* Load FFT data */ \
	f0 *= f64ptr(grpptr+off8+rcx);	/* Mul by fudged two-to-minus-phi */ \
	f0 += xcarry;			/* Add in the carry */ \
	single_rounding(base2, f0, xcarry, xtmp, rax*4+off8); \
	f0 *= f64ptr(grpptr+XMM_GMD/2+off8+rcx); /* value5 *= fudged two-to-phi */ \
	f64ptr(srcptr+4*64+off8) = f0;	/* Save value5 */ \
\
	rcx = g->BIGLIT_INCR4;/* Different clm values step through */ \
;						/* big/lit array differently */ \
	rcx = u8ptr(biglitreg+2+rcx); /* Fudge factor flags */ \
	rcx &= 0xF0; \
	xcarry -= XMM_BIGVAL1;/* Make carry an integer */ \
	xcarry *= f64ptr(grpptr+XMM_GMD/2+off8+rcx); /* carry *= fudged two-to-phi */ \
	xcarry += f64ptr(srcptr+5*64+off8);	/* Add FFT data */ \
	f64ptr(srcptr+5*64+off8) = xcarry;	/* Save value6 */ \
}

/* This macro finishes the smallmul normalize process by adding the final */ \
/* two carries back into the appropriate FFT values at the start of the fft. */ \
/* rsi = pointer to FFT data */ \
/* rbp = pointer two-to-phi group multipliers */ \
/* rdi = pointer to big/little flags */ \
/* xmm6 = non-wraparound carry */ \
/* xmm7 = wraparound carry (this will be zero for smallmul of zero-padded number) */ \


#define xnorm_smallmul_wpn_fft_zpad(base2) { \
	double f0, f1, f2, f4; \
\
;	/* Copy and integerize data from 7 words above halfway point to ZPAD0-ZPAD6 */ \
;	/* Clear words 5,6,7 */ \
;	/* Then we can make an exact copy of most of the xnorm012_wpn_zpad code */ \
\
	rax = u8ptr(rdi+1);	/* Load big vs. little flags */ \
	f0 = f64ptr(rsi+0*64+32);	/* Value1 */ \
	xmm6[0] -= XMM_BIGVAL1;/* Remove XMM_BIGVAL from carry */ \
	single_split_lower_zpad_word(base2, f0, xmm6[0], f4, rax*4); \
	g->ZPAD0_6[0] = f0; \
\
	rax = u8ptr(rdi+3);	/* Load big vs. little flags */ \
	f0 = f64ptr(rsi+1*64+32);	/* Value2 */ \
	single_split_lower_zpad_word(base2, f0, xmm6[0], f4, rax*4); \
	g->ZPAD0_6[1] = f0; \
\
	rcx = g->BIGLIT_INCR2;/* Different clm values step through */ \
;					/* big/lit array differently */ \
	rax = u8ptr(rdi+rcx+1); /* Load big vs. little flags */ \
	f0 = f64ptr(rsi+2*64+32);	/* Value3 */ \
	single_split_lower_zpad_word(base2, f0, xmm6[0], f4, rax*4); \
	g->ZPAD0_6[2] = f0; \
\
	f0 = f64ptr(rsi+3*64+32);	/* Value4 */ \
	f0 = roundf64(f0); /* Round to an integer */ \
	f0 += xmm6[0];		/* Value4 + carry */ \
	g->ZPAD0_6[3] = f0; \
\
	f0 = f64ptr(rsi+4*64+32);	/* Value5 */ \
	f0 = roundf64(f0); /* Round to an integer */ \
	g->ZPAD0_6[4] = f0; \
\
	f0 = f64ptr(rsi+5*64+32);	/* Value6 */ \
	f0 = roundf64(f0); /* Round to an integer */ \
	g->ZPAD0_6[5] = f0; \
\
	f0 = f64ptr(rsi+6*64+32);	/* Value7 */ \
	f0 = roundf64(f0); /* Round to an integer */ \
	g->ZPAD0_6[6] = f0; \
\
	f0 = 0;		/* Clear highest words */ \
	f64ptr(rsi+4*64+32) = f0; \
	f64ptr(rsi+5*64+32) = f0; \
	f64ptr(rsi+6*64+32) = f0; \
\
;	/* Divide the zpad data by k.  Store the integer part in XMM_TMP */ \
;	/* and the remainder in ZPAD0.  Later we will wrap the integer part */ \
;	/* down to the bottom of the FFT data area (and multiply by -c). */ \
;	/* And we will store the remainder in the upper half of the FFT */ \
;	/* data area. */ \
\
;	/* Note there are three cases to handle.  K is smaller than a big word. */ \
;	/* K is between one and 2 big words in size.  And K is more than */ \
;	/* 2 big words in size. */ \
\
	xfunc572(g); \
\
;	/* Now normalize the data above the halfway point. */ \
\
	rax = u8ptr(rdi+1);	/* First biglit flag  */ \
	f0 = g->ZPAD0_6[0];/* Load remainder of divide by k */ \
	f0 += XMM_BIGVAL1; \
	single_rounding(base2, f0, f2, f4, rax*4); \
	f64ptr(rsi+0*64+32) = f0;	/* Save value1 */ \
\
	rax = u8ptr(rdi+3);	/* Load big vs. little flags */ \
	single_rounding(base2, f2, f0, f4, rax*4); \
	f64ptr(rsi+1*64+32) = f2;	/* Save value2 */ \
\
	rcx = g->BIGLIT_INCR2;/* Different clm values step through */ \
;					/* big/lit array differently */ \
	rax = u8ptr(rdi+rcx+1); /* Load big vs. little flags */ \
	single_rounding(base2, f0, f2, f4, rax*4); \
	f64ptr(rsi+2*64+32) = f0;	/* Save value3 */ \
\
	f2 -= XMM_BIGVAL1;/* Remove integer rounding constant */ \
	f64ptr(rsi+3*64+32) = f2;	/* Save new value4 */ \
\
;	/* Mul the integer part of (ZPAD data divided by k) by -c in */ \
;	/* preparation for adding it into the lower FFT data area. */ \
;	/* Also add in the shifted high FFT carry at this time. */ \
\
;	/* Now add in and normalize the bottom FFT data.  Remember that the */ \
;	/* column two-to-phi multiplier for the first value will be 1.0.  We  */ \
;	/* must go 6 words deep in case k is 48-50 bits and c is 32 bits. */ \
\
	vec2f64 TMP12 = g->u.xmm.XMM_TMP1; \
	vec2f64 TMP34 = g->u.xmm.XMM_TMP2; \
	vec2f64 TMP56 = g->u.xmm.XMM_TMP3; \
	TMP12 *= g->u.xmm.XMM_MINUS_C;	/* Mul by -c */ \
	TMP34 *= g->u.xmm.XMM_MINUS_C;	/* Mul by -c */ \
	TMP56 *= g->u.xmm.XMM_MINUS_C;	/* Mul by -c */ \
	rax = u8ptr(rdi+1);	/* First word  */ \
	f0 = TMP12[0];	/* Load integer part of divide by k */ \
	f0 += XMM_BIGVAL1; \
	f0 += f64ptr(rsi+0*64);	/* Add in the FFT data */ \
	single_rounding(base2, f0, f2, f4, rax*4); \
	f64ptr(rsi+0*64) = f0;	/* Save value1 */ \
\
	rax = u8ptr(rdi+3);	/* Load big vs. little flags */ \
	f0 = TMP12[1];	/* Load integer part of divide by k */ \
	f1 = f64ptr(rsi+1*64);	/* Load FFT data */ \
	f0 += f1;		/* Add in the FFT data */ \
	f0 += f2;		/* x2 = value + carry */ \
	single_rounding(base2, f0, f2, f4, rax*4); \
	f64ptr(rsi+1*64) = f0;	/* Save value2 */ \
\
	rcx = rdi + g->BIGLIT_INCR2;/* Different clm values step through */ \
;					/* big/lit array differently */ \
	rax = u8ptr(rcx+1); /* Load big vs. little flags */ \
	f0 = TMP34[0];	/* Load integer part of divide by k */ \
	f1 = f64ptr(rsi+2*64);	/* Load FFT data */ \
	f0 += f1;		/* Add in the FFT data */ \
	f0 += f2;		/* x3 = value + carry */ \
	single_rounding(base2, f0, f2, f4, rax*4); \
	f64ptr(rsi+2*64) = f0;	/* Save value3 */ \
\
	rax = u8ptr(rcx+3); /* Load big vs. little flags */ \
	f0 = TMP34[1];	/* Load integer part of divide by k */ \
	f1 = f64ptr(rsi+3*64);	/* Load FFT data */ \
	f0 += f1;		/* Add in the FFT data */ \
	f0 += f2;		/* x4 = value + carry */ \
	single_rounding(base2, f0, f2, f4, rax*4); \
	f64ptr(rsi+3*64) = f0;	/* Save value4 */ \
\
	rcx = rdi + g->BIGLIT_INCR4;/* Different clm values step through */ \
;					/* big/lit array differently */ \
	rax = u8ptr(rcx+1); /* Load big vs. little flags */ \
	f0 = TMP56[0];	/* Load integer part of divide by k */ \
	f1 = f64ptr(rsi+4*64);	/* Load FFT data */ \
	f0 += f1;		/* Add in the FFT data */ \
	f0 += f2;		/* x4 = value + carry */ \
	single_rounding(base2, f0, f2, f4, rax*4); \
	f64ptr(rsi+4*64) = f0;	/* Save value4 */ \
\
	rax = u8ptr(rcx+3); /* Load big vs. little flags */ \
	f0 = TMP56[1];	/* Load integer part of divide by k */ \
	f1 = f64ptr(rsi+5*64);	/* Load FFT data */ \
	f0 += f1;		/* Add in the FFT data */ \
	f0 += f2;		/* x5 = value + carry */ \
	single_rounding(base2, f0, f2, f4, rax*4); \
	f64ptr(rsi+5*64) = f0;	/* Save value5 */ \
\
	f2 -= XMM_BIGVAL1;/* Remove rounding constant */ \
	f2 += f64ptr(rsi+6*64);	/* Add in FFT data */ \
	f64ptr(rsi+6*64) = f2;	/* Save value6 */ \
}
