#include <stdlib.h>
#include <ctype.h>
#include <stdio.h>

#include "../gwnum.h"
#include "../gwcommon.h"
#include "../unravel.h"
#include "../gwtables.h"
#include "../xarch.h"
#include "../xbasics.h"
#include "../xmult.h"
#include "../gwnum_func.h"
#include "hg.h"

/* Copyright 2001-2010 - Mersenne Research, Inc.  All rights reserved */
/* Author:  George Woltman */
/* Email: woltman@alum.mit.edu */
/* */
/* These macros take the basic XMM FFT building blocks and build even */
/* larger building blocks used in the final pass of a multi-pass FFT. */
/* */

/* *************** xmm-pass2-8-levels-real macro ****************** */
/* This macro takes 256 real values, 256 semi-real values and the */
/* first 256 complex values and performs the final 8 levels of the */
/* FFT process, squares the results, and does the first 8 levels of */
/* the inverse FFT. */

void hg_pass2_8_levels_complex(struct gwasm_data *__restrict g);
void xpass2_hg_8_levels(struct gwasm_data *__restrict g) {
	g->thread_work_routine = hg_pass2_8_levels_complex;
	pass2_wake_up_threads(g);
	if(g->NEGACYCLIC_FFT == 0) {

	uintptr_t rbx,rbp,rsi,rdi;
	vec2f64 xmm0,xmm1,xmm2,xmm3,xmm4,xmm5,xmm6,xmm7;

	/* We switch to the same format used in one-pass FFTs (that is, the */
	/* high word of xmm registers being 1 greater than the low word): */

	/* Load the TLBs for this data set and the next data set. */
	/* This data set is spread over 2 pages (8KB), the next data set also */
	/* takes 2 pages (8KB), and we touch one more page because the data */
	/* sets do not start on exact on 4KB boudaries. */

	start_timer(2);
	rsi = (uintptr_t)g->data_addr;/* Load source address */
	xtouch(xptr(rsi));			/* Read from this page (loads the TLB) */
	xtouch(xptr(rsi+4096));		/* Next page */
	xtouch(xptr(rsi+2*4096-128));	/* Next page */
	xtouch_init(rcx = (uintptr_t)g->data_prefetch);/* Prefetch pointer */
	xtouch(xptr(rcx));			/* Read from this page (loads the TLB) */
	xtouch(xptr(rcx+4096));		/* Next page */
	xtouch(xptr(rcx+2*4096-128));	/* Next page */

	/* Type 4 FFTs skip the forward FFT process */

	rbx = g->DIST_TO_FFTSRCARG;
	if likely(!(g->ffttype & 4)){

		/* Do FFT level 1 */
		/* */
		/* On input the 64-byte cache lines hold these data values: */
		/*	0	+128	256	+128	512	+128	768	+128 */
		/*	1	... */
		/*	... */
		/*	127	... */
		/* On output the 64-byte cache lines hold these data values: */
		/*	0	1	128	129	256	257	384	385 */
		/*	+512	... */
		/*	2	... */
		/*	... */
		/*	126	... */
		/*	+512	... */

		/* Do eight_reals_fft_1 on real values 0 - 255 */
		/* Do nothing on semi-real values from 256 - 511 */
		/* Do two_complex_fft on complex values from 512 - 1023 */
		/*	distance between fft data elements is 128 */

		loops_init(64);			/* Load loop counter */
		for(unsigned int loopA = 64; loopA; loopA--) {
			s2cl_eight_reals_fft_1(rsi, 2*64, 64);
			L2prefetch128(u8ptr(rsi+3*128));		/* Prefetch a few iterations ahead */
		}
		rsi += -64*2*64;		/* Restore source pointer */

		/* Do FFT levels 2,3 */
		/* Values 0-127 is real data, 128-255 is semi-real data, 256-1023 is */
		/* complex data. */
		/* */
		/* On input the 64-byte cache lines hold these data values: */
		/*	0	1	128	129	256	257	384	385 */
		/*	+512	... */
		/*	2	... */
		/*	... */
		/*	126	... */
		/*	+512	... */
		/* On output the 64-byte cache lines hold these data values: */
		/*	0	1	32	33	64	65	96	97 */
		/*	+512	... */
		/*	2	... */
		/*	... */
		/*	30	... */
		/*	+512	... */
		/*	128	... */
		/*	... */

		/* Do 16 eight_reals_fft_2 macros */
		/* Do 16 nop_two_two_complex_fft_2 macros */
		/* Do 32 four_complex_fft macros */
		/*	distance between fft data elements is 32 */
		/*	do 16 macros each processing 32 data values */

		rdi = (uintptr_t)g->u.xmm.sincos10;/* Load sin/cos pointer */
		loops_init(16);
		for(unsigned int loopA = 16; loopA; loopA--) {
			x4cl_eight_reals_fft_2(rsi, 2*64, 32*64, 64*64, rdi);
		}

		/* Do 64 four_complex_fft macros */
		/*	distance between fft data elements is 32 */
		/*	do 16 macros each processing 32 data values */

		rsi += -16*2*64+64;	/* Next source pointer */
		rdi += XMM_SCD;		/* Next sine/cosine pointer */
		loops_init_prefetch(16, 128, 1, rcx);
		for(unsigned int loopA = 16; loopA; loopA--) {
			x4cl_four_complex_fft(rsi, 2*64, 32*64, 64*64, rdi);
		}
		rsi += -16*2*64-64;	/* Restore source pointer */

		/* Do FFT levels 4 */
		/* Values 0-31 is real data, 32-63 is semi-real data, 64-1023 is */
		/* complex data. */
		/* */
		/* On input the 64-byte cache lines hold these data values: */
		/*	0	1	32	33	64	65	96	97 */
		/*	+512	... */
		/*	2	... */
		/*	... */
		/*	30	... */
		/*	+512	... */
		/*	128	... */
		/*	... */
		/* On output the 64-byte cache lines hold these data values: */
		/*	0	1	16	17	32	33	48	49 */
		/*	+512	... */
		/*	2	... */
		/*	... */
		/*	14	... */
		/*	+512	... */
		/*	64	... */
		/*	... */

		/* Do 16 eight_reals_fft_1 macros */
		/*	distance between fft data elements is 16 */
		/*	do 8 macros each processing 16 data values */

		rdi = (uintptr_t)g->u.xmm.sincos9;/* Load sin/cos pointer */
		loops_init(64);
		for(unsigned int loopA = 8; loopA; loopA--) {
			x2cl_eight_reals_fft_1(rsi, 2*64, 16*64);
		}

		/* Do 112 two_complex_fft macros */
		/*	distance between fft data elements is 64 */
		/*	do 56 macros each processing 16 data values */

		rsi += -8*2*64+32*64;	/* Load source pointer */
		rdi += XMM_SCD;		/* Next sine/cosine pointer */
		loops_reset;			/* 3 then 4 iterations of 4 */
		for(unsigned int loopA = 7; loopA; ) {
			do{
				for(unsigned int loopC = 8; loopC; loopC--) {
					x2cl_two_complex_fft(rsi, 2*64, 16*64, rdi);
				}
				rsi += -8*2*64+32*64;	/* Next source pointer */
				rdi += 2*XMM_SCD;		/* Next sine/cosine pointer */
			}while(--loopA & 3);
			rsi += -4*32*64+64;	/* Next source pointer */
		}
		rsi += -2*64;		/* Restore source pointer */

		/* Do FFT levels 5,6 */
		/* Values 0-15 is real data, 16-31 is semi-real data, 32-1023 is */
		/* complex data. */
		/* */
		/* On input the 64-byte cache lines hold these data values: */
		/*	0	1	16	17	32	33	48	49 */
		/*	+512	... */
		/*	2	... */
		/*	... */
		/*	14	... */
		/*	+512	... */
		/*	64	... */
		/*	... */
		/* On output the 64-byte cache lines hold these data values: */
		/*	0	1	4	5	8	9	12	13 */
		/*	+512	... */
		/*	2	... */
		/*	+512	... */
		/*	16	... */
		/*	... */

		/* Do 2 eight_reals_fft_2 macros */
		/* Do 2 nop_two_two_complex_fft_2 macros */
		/* Do 4 four_complex_fft macros */
		/*	distance between fft data elements is 4 */
		/*	do 2 macros each processing 32 data values */

		rdi = (uintptr_t)g->u.xmm.sincos7;/* Load sin/cos pointer */
		loops_init(32);
		for(unsigned int loopA = 2; loopA; loopA--) {
			x4cl_eight_reals_fft_2(rsi, 2*64, 4*64, 8*64, rdi);
		}

		/* Do 120 four_complex_fft macros */
		/*	distance between fft data elements is 4 */
		/*	do 30 macros each processing 32 data values */

		rsi += -2*2*64+16*64;	/* Load source pointer */
		rdi += XMM_SCD;		/* Next sine/cosine pointer */
		loops_reset;			/* 7 then 8 iterations of 2 */
		for(unsigned int loopA = 15; loopA; ) {
			do{
				for(unsigned int loopC = 2; loopC; loopC--) {
					x4cl_four_complex_fft(rsi, 2*64, 4*64, 8*64, rdi);
				}
				rsi += -2*2*64+16*64;	/* Next source pointer */
				rdi += 2*XMM_SCD;		/* Next sine/cosine pointer */
			}while(--loopA & 7);
			rsi += -8*16*64+64;	/* Next source pointer */
		}
		rsi += -2*64;		/* Restore source pointer */

		/* Do FFT levels 7,8 */
		/* Values 0-3 is real data, 4-7 is semi-real data, 8-1023 is */
		/* complex data. */
		/* */
		/* On input the 64-byte cache lines hold these data values: */
		/*	0	1	4	5	8	9	12	13 */
		/*	+512	... */
		/*	2	... */
		/*	+512	... */
		/*	16	... */
		/*	... */
		/* On output the 64-byte cache lines hold these data values: */
		/*	0	1	2	3	4	5	6	7 */
		/*	+512	... */
		/*	8	... */
		/*	... */

		/* Do 1/2 eight_reals_fft_2 macros */
		/* Do 1/2 nop_two_two_complex_fft_2 macros */
		/* Do 1 four_complex_fft macros */
		/*	distance between fft data elements is 1 */

		/* Do 126 more four_complex_fft macros */
		/*	distance between fft data elements is 1 */
		/*	do 63 macros each processing 16 data values */

		rdi = (uintptr_t)g->u.xmm.sincos6;/* Load sin/cos pointer */

		/* Execute the proper middle step */

		if likely(g->ffttype & 2) {
			if likely(!(g->ffttype & 1)) {
				s2cl_eight_reals_with_square_2(rsi, 4*64, 2*64);
				rdi += XMM_SCD;		/* Next sine/cosine pointer */
				loops_init(64, 1, 32);		/* 31 then 32 iterations */
				for(unsigned int loopA = 63; loopA; ) {
					do{
						s2cl_four_complex_with_square(rsi, 4*64, 2*64);
						rdi += XMM_SCD;		/* Next sine/cosine pointer */
					}while(--loopA & 31);
					rsi += -32*4*64+64;	/* Next source pointer */
				}
				rsi += -2*64;		/* Restore source pointer */
			}else{
				rbp = g->DIST_TO_MULSRCARG;
				s2cl_eight_reals_with_mult_2(rsi, 4*64, 2*64);
				rdi += XMM_SCD;		/* Next sine/cosine pointer */
				loops_init(64, 1, 32);		/* 31 then 32 iterations */
				for(unsigned int loopA = 63; loopA; ) {
					do{
						s2cl_four_complex_with_mult(rsi, 4*64, 2*64);
						rdi += XMM_SCD;		/* Next sine/cosine pointer */
					}while(--loopA & 31);
					rsi += -32*4*64+64;	/* Next source pointer */
				}
				rsi += -2*64;		/* Restore source pointer */
			}
		}else{
			s2cl_eight_reals_fft_2_final(rsi, 4*64, 2*64);
			rdi += XMM_SCD;		/* Next sine/cosine pointer */
			loops_init(64, 1, 32);		/* 31 then 32 iterations */
			for(unsigned int loopA = 63; loopA; ) {
				do{
					s2cl_four_complex_fft_final(rsi, 4*64, 2*64);
					rdi += XMM_SCD;		/* Next sine/cosine pointer */
				}while(--loopA & 31);
				rsi += -32*4*64+64;	/* Next source pointer */
			}
			goto	xpass2_8_real_done;
		}
	}else{
		rbp = g->DIST_TO_MULSRCARG;
		rdi = (uintptr_t)g->u.xmm.sincos6;/* Load sin/cos pointer */
		s2cl_eight_reals_with_mulf_2(rsi, 4*64, 2*64);
		rdi += XMM_SCD;		/* Next sine/cosine pointer */
		loops_init(64, 1, 32);		/* 31 then 32 iterations */
		for(unsigned int loopA = 63; loopA; ) {
			do{
				s2cl_four_complex_with_mulf(rsi, 4*64, 2*64);
				rdi += XMM_SCD;		/* Next sine/cosine pointer */
			}while(--loopA & 31);
			rsi += -32*4*64+64;	/* Next source pointer */
		}
		rsi += -2*64;		/* Restore source pointer */

		/* Do inverse FFT levels 5,6 */
		/* On input the 64-byte cache lines hold these data values: */
		/*	0	1	2	3	4	5	6	7 */
		/*	+512	... */
		/*	8	... */
		/*	... */
		/* On output the 64-byte cache lines hold these data values: */
		/*	0	1	8	9	16	17	24	25 */
		/*	+512	... */
		/*	2	... */
		/*	... */
		/*	6	... */
		/*	+512	... */
		/*	32	... */
		/*	... */

		/* Do 2 eight_reals_unfft_2 macros */
		/* Do 2 nop_two_two_complex_unfft_2 macros */
		/*	distance between fft data elements is 4 */
		/*	do 1 macro processing 32 data values */

	}
	x4cl_eight_reals_unfft_2(rsi, 8*64, 2*64, 4*64);

	/* Do 124 four_complex_unfft macros */
	/*	distance between fft data elements is 4 */
	/*	do 31 macros each processing 32 data values */

	rdi = (uintptr_t)g->u.xmm.sincos7;/* Load sin/cos pointer */
	loops_init(32, 1, 16);		/* 15 then 16 iterations */
	for(unsigned int loopA = 2*16-1; loopA; ) {
		do{
			x4cl_four_complex_unfft(rsi, 8*64, 2*64, 4*64, rdi);
			rdi += XMM_SCD;		/* Next sine/cosine pointer */
		}while(--loopA & 15);
		rsi += -16*8*64+64;	/* Next source pointer */
	}
	rsi += -2*64;		/* Restore source pointer */

	/* Do inverse FFT levels 4 */
	/* On input the 64-byte cache lines hold these data values: */
	/*	0	1	8	9	16	17	24	25 */
	/*	+512	... */
	/*	2	... */
	/*	... */
	/*	6	... */
	/*	+512	... */
	/*	32	... */
	/*	... */
	/* On output the 64-byte cache lines hold these data values: */
	/*	0	1	16	17	32	33	48	49 */
	/*	+512	... */
	/*	2	... */
	/*	... */
	/*	14	... */
	/*	+512	... */
	/*	64	... */
	/*	... */

	/* Do 8 eight_reals_unfft_1 macros */
	/*	distance between fft data elements is 16 */
	/*	do 4 macros each processing 16 data values */

	loops_init(64);
	for(unsigned int loopA = 4; loopA; loopA--) {
		x2cl_eight_reals_unfft_1(rsi, 2*64, 8*64);
	}

	/* Do 960 two_complex_unfft macros */
	/*	distance between fft data elements is 16 */
	/*	do 240 macros each processing 32 data values */

	rsi += -4*2*64+16*64;	/* Load source pointer */
	rdi = (uintptr_t)g->u.xmm.sincos9;/* Load sin/cos pointer */
	loops_reset;			/* 7 then 8 iterations of 4 */
	for(unsigned int loopA = 2*8-1; loopA; ) {
		do{
			for(unsigned int loopC = 4; loopC; loopC--) {
				x2cl_two_complex_unfft(rsi, 2*64, 8*64);
			}
			rsi += -4*2*64+16*64;	/* Next source pointer */
			rdi += XMM_SCD;		/* Next sine/cosine pointer */
		}while(--loopA & 7);
		rsi += -8*16*64+64;	/* Next source pointer */
	}
	rsi += -2*64;		/* Restore source pointer */

	/* Do inverse FFT levels 2,3 */
	/* On input the 64-byte cache lines hold these data values: */
	/*	0	1	16	17	32	33	48	49 */
	/*	+512	... */
	/*	2	... */
	/*	... */
	/*	14	... */
	/*	+512	... */
	/*	64	... */
	/*	... */
	/* On output the 64-byte cache lines hold these data values: */
	/*	0	1	64	65	128	129	192	193 */
	/*	+512	... */
	/*	2	... */
	/*	... */
	/*	62	... */
	/*	+512	... */
	/*	256	... */
	/*	... */

	/* Do 128 eight_reals_unfft_2 macros */
	/* Do 128 nop_two_two_complex_unfft_2 macros */
	/*	distance between fft data elements is 256 */
	/*	do 8 macros each processing 32 data values */

	loops_init(32);
	for(unsigned int loopA = 8; loopA; loopA--) {
		x4cl_eight_reals_unfft_2(rsi, 2*64, 16*64, 32*64);
	}

	/* Do 96 four_complex_unfft macros */
	/*	distance between fft data elements is 32 */
	/*	do 24 macros each processing 32 data values */

	rsi += -8*2*64+64*64;	/* Next source pointer */
	rdi = (uintptr_t)g->u.xmm.sincos10;/* Load sin/cos pointer */
	loops_reset;			/* 1 then 2 iterations of 64 */
	for(unsigned int loopA = 2*2-1; loopA; ) {
		do{
			for(unsigned int loopC = 8; loopC; loopC--) {
				x4cl_four_complex_unfft(rsi, 2*64, 16*64, 32*64, rdi);
			}
			rsi += -8*2*64+64*64;	/* Next source pointer */
			rdi += XMM_SCD;		/* Next sine/cosine pointer */
		}while(--loopA & 1);
		rsi += -2*64*64+64;	/* Next source pointer */
	}
	rsi += -2*64;		/* Restore source pointer */

	/* Do inverse FFT level 1 */
	/* On input the 64-byte cache lines hold these data values: */
	/*	0	1	64	65	128	129	192	193 */
	/*	+512	... */
	/*	2	... */
	/*	... */
	/*	62	... */
	/*	+512	... */
	/*	256	... */
	/*	... */
	/* On output the 64-byte cache lines hold these data values: */
	/*	0	+128	256	+128	512	+128	768	+128 */
	/*	1	... */
	/*	... */
	/*	127	... */

	/* Do 64 eight_reals_last_unfft macros */
	/* Do 64 two_two_complex_unfft macros */
	/*	distance between fft data elements is 128 */
	/*	do 32 macros each processing 16 data values */

	loops_init(32);
	for(unsigned int loopA = 32; loopA; loopA--) {
		s4cl_eight_reals_unfft_1(rsi, 2*64, 64*64, 64);
	}

	xpass2_8_real_done:
	end_timer(2);

	if(g->pass2_get_next_block(g) != 0) return;
}
hg_pass2_8_levels_complex(g);
}


/* *************** xmm-pass2-8-levels-complex macro ****************** */
/* This macro takes 2 groups of 256 complex values, premultiplies them and */
/* performs the final 8 levels of the FFT process, squares the */
/* results, and does the first 8 levels of the inverse FFT. */
/* NOTE: Rather than remembering 256 multipliers we break them up */
/* into 64 groups of 4 columns.  By remembering these 68 values we can compute */
/* all 256 multipliers using a single complex multiply.  These extra */
/* complex multiplies can probably be done faster than the computer can read the */
/* extra multipliers from memory (and are more likely stay in the L2 cache) */

void hg_pass2_8_levels_complex(struct gwasm_data *__restrict g) {
do{
	uintptr_t rbx,rbp,rsi,rdi;
	vec2f64 xmm0,xmm1,xmm2,xmm3,xmm4,xmm5,xmm6,xmm7;

	/* Load the TLBs for this data set and the next data set. */
	/* This data set is spread over 2 pages (8KB), the next data set also */
	/* takes 2 pages (8KB), and we touch one more page because the data */
	/* sets do not start on exact on 4KB boudaries.  Also, preload TLBs for */
	/* premultiplier data. */

	start_timer(3);
	xtouch_init(rcx = (uintptr_t)g->data_prefetch);/* Prefetch pointer */
	xtouch(xptr(rcx));			/* Read from this page (loads the TLB) */
	xtouch(xptr(rcx+4096));		/* Next page */
	xtouch(xptr(rcx+2*4096-128));	/* Next page */
	xtouch_init(rcx = (uintptr_t)g->premult_prefetch);/* Prefetch pointer */
	xtouch(xptr(rcx));			/* Load premultiplier TLBs */
	xtouch(xptr(rcx+4096-128));		/* Load premultiplier TLBs */
	end_timer(3);

	/* Type 4 FFTs skip the forward FFT process */

	rsi = (uintptr_t)g->data_addr;/* Load address of FFT data */
	rbx = g->DIST_TO_FFTSRCARG;
	if likely(!(g->ffttype & 4)){

		/* Do FFT level 1,2 */
		/* */
		/* We swizzle the data so that one set of data is in the high word of */
		/* an xmm register, the other set is in the low word.  This means */
		/* a distance of 4096 data values. */
		/* */
		/* On input the 64-byte cache lines hold these data values: */
		/*	0	+128	256	+128	512	+128	768	+128 */
		/*	1	... */
		/*	... */
		/*	127 */
		/* On output the 64-byte cache lines hold these data values: */
		/*	0	+512	64	+512	128	+512	192	+512 */
		/*	1	... */
		/*	... */
		/*	63	... */
		/*	256	... */
		/*	... */

		/* Do 128 four_complex_gpm_fft macros */
		/*	distance between fft data elements is 64 */
		/*	do 64 macros each processing 16 data values */

		start_timer(4);
		rdi = (uintptr_t)g->premult_addr;/* The group multipliers */
		loops_init(64);			/* 16 iterations of 4 */
		for(unsigned int loopA = 16; loopA; loopA--) {
			for(unsigned int loopB = 4; loopB; loopB--) {
				s2cl_four_complex_gpm_fft(rsi, 64, 64*64);
				L2prefetch128(u8ptr(rsi+5*64));		/* Prefetch a few iterations ahead in */
				L2prefetch128(u8ptr(rsi+64*64+5*64));	/* case they are not in the L2 cache */
			}
			rdi += XMM_PMD;		/* Next group multiplier */
			L2prefetch128(u8ptr(rdi+4*XMM_PMD));
			L2prefetch128(u8ptr(rdi+16*XMM_PMD));
		}
		rsi += -64*64;		/* Restore source pointer */
		end_timer(4);

		/* Do FFT levels 3,4 */
		/* */
		/* On input the 64-byte cache lines hold these data values: */
		/*	0	+512	64	+512	128	+512	192	+512 */
		/*	1	... */
		/*	... */
		/*	63	... */
		/*	256	... */
		/*	... */
		/* On output the 64-byte cache lines hold these data values: */
		/*	0	+512	16	+512	32	+512	48	+512 */
		/*	1	... */
		/*	... */
		/*	15	... */
		/*	64	... */
		/*	... */

		/* Do 4 groups of four_complex_fft_cpm macros */
		/* To make things run faster we precompute the column multipliers times the */
		/* 4 different sine/cosine multipliers we run into in this section.  This */
		/* costs some memory but saves us some complex multiplies. */

		/* Do 128 four_complex_cpm_fft macros */
		/*	distance between fft data elements is 16 */
		/*	do 16 macros each processing 64 data values */

		start_timer(5);
		xtouch_init(rcx = (uintptr_t)g->premult_prefetch);/* Prefetch pointer */
		loops_init_prefetch(32, 128, 1, rcx); /* 2 iters of 4 iterations of 4 */
		for(unsigned int loopA = 2; loopA; loopA--) {
			for(unsigned int loopB = 4; loopB; loopB--) {
				for(unsigned int loopC = 4; loopC; loopC--) {
					x4cl_four_complex_cpm_fft(rsi, 64, 16*64, 32*64, 4*XMM_PMD);
					rdi += XMM_PMD;		/* Next column premultiplier */
				}
				rdi += -4*XMM_PMD;		/* Restore column premultiplier */
			}
			rsi += -16*64+64*64;	/* Next source pointer */
			rdi += 8*XMM_PMD;		/* Next set of premultipliers */
		}
		rsi += -2*64*64;		/* Restore source pointer */
		end_timer(5);

		/* Do FFT levels 5,6 */
		/* */
		/* On input the 64-byte cache lines hold these data values: */
		/*	0	+512	16	+512	32	+512	48	+512 */
		/*	1	... */
		/*	... */
		/*	15	... */
		/*	64	... */
		/*	... */
		/* On output the 64-byte cache lines hold these data values: */
		/*	0	+512	4	+512	8	+512	12	+512 */
		/*	1	... */
		/*	... */
		/*	3	... */
		/*	16	... */
		/*	... */

		/* Do 128 four_complex_fft macros */
		/*	distance between fft data elements is 4 */
		/*	do 32 macros each processing 32 data values */

		start_timer(8);
		xtouch_init(rcx = (uintptr_t)g->data_prefetch);/* Prefetch pointer */
		rdi = (uintptr_t)g->xsincos_complex;	/* Load sin/cos pointer */
		loops_init_prefetch(32, 64, 1, rcx); /* 8 iterations of 4, prefetch 64 bytes every macro */
		for(unsigned int loopA = 8; loopA; loopA--) {
			for(unsigned int loopB = 4; loopB; loopB--) {
				x4cl_four_complex_fft(rsi, 64, 4*64, 8*64, rdi);
			}
			rsi += -4*64+16*64;	/* Next source pointer */
			rdi += 2*XMM_SCD;		/* Next sine/cosine pointer */
		}
		rsi += -8*16*64;		/* Restore source pointer */
		end_timer(8);

		/* Do FFT levels 7,8 as well as inverse FFT levels 7,8 */
		/* */
		/* On input the 64-byte cache lines hold these data values: */
		/*	0	+512	4	+512	8	+512	12	+512 */
		/*	1	... */
		/*	... */
		/*	3	... */
		/*	16	... */
		/*	... */
		/* On output the 64-byte cache lines hold these data values: */
		/*	0	+512	2	+512	4	+512	6	+512 */
		/*	1	... */
		/*	8	... */
		/*	... */

		/* Do 128 four_complex_fft macros */
		/*	distance between fft data elements is 1 */
		/*	do 32 macros each processing 32 data values */

		/* Execute the right middle step */

		if likely(g->ffttype & 2) {
			if(!(g->ffttype & 1)) {
				start_timer(9);
				loops_init_prefetch(32, 64, 1, rcx);
				for(unsigned int loopA = 4; loopA; loopA--) {
					rdi = (uintptr_t)g->xsincos_complex;	/* Load sin/cos pointer */
					for(unsigned int loopB = 8; loopB; loopB--) {
						x4cl_four_complex_with_square(rsi, 4*64, 64, 2*64);
						rdi += 2*XMM_SCD;		/* Next sine/cosine pointer */
					}
				}
				end_timer(9);
				rsi += -32*4*64;		/* Restore source pointer */
			}else{
				start_timer(9);
				rbp = g->DIST_TO_MULSRCARG;
				loops_init_prefetch(32, 64, 1, rcx);
				for(unsigned int loopA = 4; loopA; loopA--) {
					rdi = (uintptr_t)g->xsincos_complex;	/* Load sin/cos pointer */
					for(unsigned int loopB = 8; loopB; loopB--) {
						x4cl_four_complex_with_mult(rsi, 4*64, 64, 2*64);
						rdi += 2*XMM_SCD;		/* Next sine/cosine pointer */
					}
				}
				end_timer(9);
				rsi += -32*4*64;		/* Restore source pointer */
			}
		}else{
			start_timer(9);
			loops_init_prefetch(32, 128, 1, rcx);
			for(unsigned int loopA = 4; loopA; loopA--) {
				rdi = (uintptr_t)g->xsincos_complex;	/* Load sin/cos pointer */
				for(unsigned int loopB = 8; loopB; loopB--) {
					x4cl_four_complex_fft_final(rsi, 4*64, 64, 2*64);
					rdi += 2*XMM_SCD;		/* Next sine/cosine pointer */
				}
			}
			end_timer(9);
			goto	xpass2_8_levels_complex_done;
		}
	}else{
		start_timer(9);
		rbp = g->DIST_TO_MULSRCARG;
		xtouch_init(rcx = (uintptr_t)g->data_prefetch);/* Prefetch pointer */
		loops_init_prefetch(32, 128, 1, rcx);
		for(unsigned int loopA = 4; loopA; loopA--) {
			rdi = (uintptr_t)g->xsincos_complex;	/* Load sin/cos pointer */
			for(unsigned int loopB = 8; loopB; loopB--) {
				x4cl_four_complex_with_mulf(rsi, 4*64, 64, 2*64);
				rdi += 2*XMM_SCD;		/* Next sine/cosine pointer */
			}
		}
		end_timer(9);
		rsi += -32*4*64;		/* Restore source pointer */

		/* Do inverse FFT levels 5,6 */
		/* On input the 64-byte cache lines hold these data values: */
		/*	0	+512	2	+512	4	+512	6	+512 */
		/*	1	... */
		/*	8	... */
		/*	... */
		/* On output the 64-byte cache lines hold these data values: */
		/*	0	+512	8	+512	16	+512	24	+512 */
		/*	1	... */
		/*	... */
		/*	7	... */
		/*	32	... */
		/*	... */

		/* Do 128 four_complex_unfft macros */
		/*	distance between fft data elements is 4 */
		/*	do 32 macros each processing 32 data values */

	}
	start_timer(10);
	rdi = (uintptr_t)g->xsincos_complex;	/* Load sin/cos pointer */
	loops_init_prefetch(32, 64, 1, rcx);
	for(unsigned int loopA = 16; loopA; loopA--) {
		for(unsigned int loopB = 2; loopB; loopB--) {
			x4cl_four_complex_unfft(rsi, 64, 2*64, 4*64, rdi);
		}
		rsi += -2*64+8*64;		/* Next source pointer */
		rdi += XMM_SCD;		/* Next sine/cosine pointer */
	}
	rsi += -16*8*64;		/* Restore source pointer */
	end_timer(10);

	/* Do inverse FFT levels 3,4 */
	/* On input the 64-byte cache lines hold these data values: */
	/*	0	+512	8	+512	16	+512	24	+512 */
	/*	1	... */
	/*	... */
	/*	7	... */
	/*	32	... */
	/*	... */
	/* On output the 64-byte cache lines hold these data values: */
	/*	0	+512	32	+512	64	+512	96	+512 */
	/*	1	... */
	/*	... */
	/*	31	... */
	/*	128	... */
	/*	... */

	/* Do 4 groups of 32 four_complex_fft macros (FFT levels 3,4) */
	/* In this section we also apply the 4 column multipliers (actually three */
	/* since the first column multiplier is always one).  To make things */
	/* run faster we precompute the column multipliers times the 4 different */
	/* sine/cosine multipliers we run into in this section.  This costs some */
	/* memory but saves us some complex multiplies. */

	/* Do 128 four_complex_cpm_unfft macros */
	/*	distance between fft data elements is 16 */
	/*	do 8 macros each processing 128 data values */

	start_timer(13);
	rdi = (uintptr_t)g->premult_addr;/* Load column multipliers pointer */
	rdi += 16*XMM_PMD;		/* Load column multipliers pointer */
	loops_init_prefetch(32, 64, 1, rcx); /* 4 iters of 2 iterations of 4 */
	for(unsigned int loopA = 4; loopA; loopA--) {
		for(unsigned int loopB = 2; loopB; loopB--) {
			for(unsigned int loopC = 4; loopC; loopC--) {
				x4cl_four_complex_cpm_unfft(rsi, 64, 8*64, 16*64);
				rdi += XMM_PMD;		/* Next column premultiplier */
			}
			rdi += -4*XMM_PMD;		/* Reset premultiplier pointer */
		}
		rsi += -8*64+32*64;	/* Next source pointer */
		rdi += 4*XMM_PMD;		/* Next premultiplier pointer */
	}
	rsi += -4*32*64;		/* Restore source pointer */
	end_timer(13);

	/* Do inverse FFT levels 1,2 */
	/* On input the 64-byte cache lines hold these data values: */
	/*	0	+512	32	+512	64	+512	96	+512 */
	/*	1	... */
	/*	... */
	/*	31	... */
	/*	128	... */
	/*	... */
	/* On output the 64-byte cache lines hold these data values: */
	/*	0	+128	256	+128	512	+128	768	+128 */
	/*	1	... */
	/*	... */
	/*	255 */

	/* Do 32 four_complex_gpm_unfft macros */
	/*	distance between fft data elements is 64 */
	/*	do 32 macros each processing 32 data values */

	start_timer(14);
	xtouch_init(rcx = (uintptr_t)g->premult_prefetch);/* Prefetch pointer */
	rdi = (uintptr_t)g->premult_addr;/* Load group multiplier pointer */
	loops_init_prefetch(32, 64, 1, rcx); /* 8 iterations of 4 */
	for(unsigned int loopA = 8; loopA; loopA--) {
		for(unsigned int loopB = 4; loopB; loopB--) {
			s4cl_four_complex_gpm_unfft(rsi, 64, 32*64, 64*64, 8*XMM_PMD);
		}
		rdi += XMM_PMD;		/* Next set of multipliers */
	}
	end_timer(14);

	xpass2_8_levels_complex_done: ;
}while(g->pass2_get_next_block(g) == 0);
}
