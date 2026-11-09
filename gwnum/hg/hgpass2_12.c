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

/* *************** xmm-pass2-12-levels-real macro ****************** */
/* This macro takes 4096 real values, 4096 semi-real values and the */
/* first 4096 complex values and performs the final 12 levels of the */
/* FFT process, squares the results, and does the first 12 levels of */
/* the inverse FFT. */

void hg_pass2_12_levels_complex(struct gwasm_data *__restrict g);
void xpass2_hg_12_levels(struct gwasm_data *__restrict g) {
g->thread_work_routine = hg_pass2_12_levels_complex;
pass2_wake_up_threads(g);
if(g->NEGACYCLIC_FFT == 0) {

	uintptr_t rbx,rbp,rsi,rdi;
	vec2f64 xmm0,xmm1,xmm2,xmm3,xmm4,xmm5,xmm6,xmm7;

	/* We switch to the same format used in one-pass FFTs (that is, the */
	/* high word of xmm registers being 1 greater than the low word): */

	/* Type 4 FFTs skip the forward FFT process */

	start_timer(2);
	rsi = (uintptr_t)g->data_addr;/* Load source address */
	rbx = g->DIST_TO_FFTSRCARG;
	if likely(!(g->ffttype & 4)){

		/* Do FFT level 1 */
		/* */
		/* On input the 64-byte cache lines hold these data values: */
		/*	0	+2K	4096	+2K	8192	+2K	12288	+2K */
		/*	1	... */
		/*	... */
		/*	2047	... */
		/* On output the 64-byte cache lines hold these data values: */
		/*	0	+1	2048	+1	4096	+1	6144	+1 */
		/*	+8192	... */
		/*	2	... */
		/*	... */
		/*	2046	... */
		/*	+8192	... */

		/* Do eight_reals_fft_1 on real values 0 - 4095 */
		/* Do nothing on semi-real values from 4096 - 8191 */
		/* Do two_complex_fft on complex values from 8192 - 16383 */
		/*	distance between fft data elements is 2K */
		/*	do 1024 iterations */

		loops_init(1024);
		for(unsigned int loopA = 16; loopA; loopA--) {
			for(unsigned int loopB = 64; loopB; loopB--) {
				s2cl_eight_reals_fft_1(rsi, 2*64, 64);
				L2prefetch128(u8ptr(rsi+3*128));		/* Prefetch a few iterations ahead */
			}
			rsi += -64*2*64+dist128;	/* Next source pointer */
		}
		rsi += -16*dist128;	/* Restore source pointer */

		/* Do FFT levels 2,3 */
		/* Values 0-2047 is real data, 2048-4195 is semi-real data, 4096-16383 is */
		/* complex data. */
		/* */
		/* On input the 64-byte cache lines hold these data values: */
		/*	0	+1	2048	+1	4096	+1	6144	+1 */
		/*	+8192	... */
		/*	2	... */
		/*	... */
		/*	2046	... */
		/*	+8192	... */
		/* On output the 64-byte cache lines hold these data values: */
		/*	0	+1	512	+1	1024	+1	1536	+1 */
		/*	+8192	... */
		/*	2	... */
		/*	... */
		/*	510	... */
		/*	+8192	... */
		/*	2048	... */
		/*	... */

		/* Do 256 eight_reals_fft_2 macros */
		/* Do 256 nop_two_two_complex_fft_2 macros */
		/* Do 512 four_complex_fft macros */
		/*	distance between fft data elements is 512 */
		/*	do 256 macros each processing 32 data values */

		rdi = (uintptr_t)g->u.xmm.sincos11;/* Load sin/cos pointer */
		xtouch_init(rcx = (uintptr_t)g->data_prefetch);/* Prefetch pointer */
		loops_init(256);
		for(unsigned int loopA = 4; loopA; loopA--) {
			for(unsigned int loopB = 64; loopB; loopB--) {
				x4cl_eight_reals_fft_2(rsi, 2*64, 4*dist128, 8*dist128, rdi);
			}
			rsi += -64*2*64+dist128;	/* Next source pointer */
		}

		/* Do 1024 four_complex_fft macros */
		/*	distance between fft data elements is 512 */
		/*	do 256 macros each processing 32 data values */

		rsi += -4*dist128+64;	/* Next source pointer */
		rdi += XMM_SCD;		/* Next sine/cosine pointer */
		loops_init_prefetch(256, 128, 1, rcx);
		xtouch(xptr(rcx));			/* Preload first TLB */
		for(unsigned int loopA = 4; loopA; loopA--) {
			xtouch(xptr(rcx+4096-128));		/* Preload the TLBs */
			xtouch(xptr(rcx+2*4096-128));
			for(unsigned int loopB = 64; loopB; loopB--) {
				x4cl_four_complex_fft(rsi, 2*64, 4*dist128, 8*dist128, rdi);
			}
			rsi += -64*2*64+dist128;	/* Next source pointer */
			xtouch_init(rcx += 128);		/* Skip pad bytes */
		}
		rsi += -4*dist128-64;	/* Restore source pointer */

		/* Do FFT levels 4 */
		/* Values 0-511 is real data, 512-1023 is semi-real data, 1024-16383 is */
		/* complex data. */
		/* */
		/* On input the 64-byte cache lines hold these data values: */
		/*	0	+1	512	+1	1024	+1	1536	+1 */
		/*	+8192	... */
		/*	2	... */
		/*	... */
		/*	510	... */
		/*	+8192	... */
		/*	2048	... */
		/*	... */
		/* On output the 64-byte cache lines hold these data values: */
		/*	0	+1	256	+1	512	+1	768	+1 */
		/*	+8192	... */
		/*	2	... */
		/*	... */
		/*	254	... */
		/*	+8192	... */
		/*	1024	... */
		/*	... */

		/* Do 256 eight_reals_fft_1 macros */
		/*	distance between fft data elements is 256 */
		/*	do 128 macros each processing 16 data values */

		rdi = (uintptr_t)g->u.xmm.sincos10;/* Load sin/cos pointer */
		loops_init(1024);
		for(unsigned int loopA = 2; loopA; loopA--) {
			for(unsigned int loopB = 64; loopB; loopB--) {
				x2cl_eight_reals_fft_1(rsi, 2*64, 2*dist128);
			}
			rsi += -64*2*64+dist128;	/* Next source pointer */
		}

		/* Do 1792 two_complex_fft macros */
		/*	distance between fft data elements is 256 */
		/*	do 896 macros each processing 16 data values */

		rsi += -2*dist128+4*dist128; /* Load source pointer */
		rdi += XMM_SCD;		/* Next sine/cosine pointer */
		loops_reset;			/* 3 then 4 iterations of 128 */
		for(unsigned int loopA = 2*4-1; loopA; ) {
			do{
				for(unsigned int loopC = 2; loopC; loopC--) {
					for(unsigned int loopD = 64; loopD; loopD--) {
						x2cl_two_complex_fft(rsi, 2*64, 2*dist128, rdi);
					}
					rsi += -64*2*64+dist128;	/* Next source pointer */
				};
				rsi += -2*dist128+4*dist128; /* Next source pointer */
				rdi += 2*XMM_SCD;		/* Next sine/cosine pointer */
			}while(--loopA & 3);
			rsi += -4*4*dist128+64;	/* Next source pointer */
		}
		rsi += -2*64;		/* Restore source pointer */

		/* Do FFT levels 5,6 */
		/* Values 0-255 is real data, 256-511 is semi-real data, 512-16384 is */
		/* complex data. */
		/* */
		/* On input the 64-byte cache lines hold these data values: */
		/*	0	1	256	+1	512	+1	768	+1 */
		/*	+8192	... */
		/*	2	... */
		/*	... */
		/*	254	... */
		/*	+8192	... */
		/*	1024	... */
		/*	... */
		/* On output the 64-byte cache lines hold these data values: */
		/*	0	1	64	+1	128	+1	192	+1 */
		/*	+8192	... */
		/*	2	... */
		/*	... */
		/*	62	... */
		/*	+8192	... */
		/*	256	... */
		/*	... */

		/* Do 32 eight_reals_fft_2 macros */
		/* Do 32 nop_two_two_complex_fft_2 macros */
		/* Do 64 four_complex_fft macros */
		/*	distance between fft data elements is 64 */
		/*	do 32 macros each processing 32 data values */

		rdi = (uintptr_t)g->u.xmm.sincos9;/* Load sin/cos pointer */
		loops_init(512);
		for(unsigned int loopA = 32; loopA; loopA--) {
			x4cl_eight_reals_fft_2(rsi, 2*64, 64*64, dist128, rdi);
		}

		/* Do 1920 four_complex_fft macros */
		/*	distance between fft data elements is 64 */
		/*	do 480 macros each processing 32 data values */

		rsi += -32*2*64+2*dist128;	/* Load source pointer */
		rdi += XMM_SCD;		/* Next sine/cosine pointer */
		loops_reset;			/* 7 then 8 iterations of 32 */
		for(unsigned int loopA = 2*8-1; loopA; ) {
			do{
				for(unsigned int loopC = 32; loopC; loopC--) {
					x4cl_four_complex_fft(rsi, 2*64, 64*64, dist128, rdi);
				}
				rsi += -32*2*64+2*dist128;	/* Next source pointer */
				rdi += 2*XMM_SCD;		/* Next sine/cosine pointer */
			}while(--loopA & 7);
			rsi += -8*2*dist128+64;	/* Next source pointer */
		}
		rsi += -2*64;		/* Restore source pointer */

		/* Do FFT levels 7,8 */
		/* Values 0-63 is real data, 64-127 is semi-real data, 128-16384 is */
		/* complex data. */
		/* */
		/* On input the 64-byte cache lines hold these data values: */
		/*	0	1	64	+1	128	+1	192	+1 */
		/*	+8192	... */
		/*	2	... */
		/*	... */
		/*	62	... */
		/*	+8192	... */
		/*	256	... */
		/*	... */
		/* On output the 64-byte cache lines hold these data values: */
		/*	0	1	16	+1	32	+1	48	+1 */
		/*	+8192	... */
		/*	2	... */
		/*	... */
		/*	14	... */
		/*	+8192	... */
		/*	64	... */
		/*	... */

		/* Do 8 eight_reals_fft_2 macros */
		/* Do 8 nop_two_two_complex_fft_2 macros */
		/* Do 16 four_complex_fft macros */
		/*	distance between fft data elements is 16 */
		/*	do 8 macros each processing 32 data values */

		rdi = (uintptr_t)g->u.xmm.sincos8;/* Load sin/cos pointer */
		loops_init(512);
		for(unsigned int loopA = 8; loopA; loopA--) {
			x4cl_eight_reals_fft_2(rsi, 2*64, 16*64, 32*64, rdi);
		}

		/* Do 2016 four_complex_fft macros */
		/*	distance between fft data elements is 16 */
		/*	do 504 macros each processing 32 data values */

		rsi += -8*2*64+64*64;	/* Load source pointer */
		rdi += XMM_SCD;		/* Next sine/cosine pointer */
		loops_reset;			/* 2 iters of 16 iters of 1 or 2 of 8 */
		for(unsigned int loopA = 2*16*2-1; loopA; ) {
			do{
				do{
					for(unsigned int loopD = 8; loopD; loopD--) {
						x4cl_four_complex_fft(rsi, 2*64, 16*64, 32*64, rdi);
					}
					rsi += -8*2*64+64*64;	/* Next source pointer */
					rdi += 2*XMM_SCD;		/* Next sine/cosine pointer */
				}while(--loopA & 1);
				rsi += -2*64*64+dist128;	/* Next source pointer */
			}while(loopA & (16*2-1));
			rsi += -16*dist128+64;	/* Next source pointer */
		}
		rsi += -2*64;		/* Restore source pointer */

		/* Do FFT levels 9,10 */
		/* Values 0-15 is real data, 16-31 is semi-real data, 32-16383 is */
		/* complex data. */
		/* */
		/* On input the 64-byte cache lines hold these data values: */
		/*	0	1	16	+1	32	+1	48	+1 */
		/*	+8192	... */
		/*	2	... */
		/*	... */
		/*	14	... */
		/*	+8192	... */
		/*	64	... */
		/*	... */
		/* On output the 64-byte cache lines hold these data values: */
		/*	0	1	4	+1	8	+1	12	+1 */
		/*	+8192	... */
		/*	2	... */
		/*	+8192	... */
		/*	16	... */
		/*	... */

		/* Do 2 eight_reals_fft_2 macros */
		/* Do 2 nop_two_two_complex_fft_2 macros */
		/* Do 4 four_complex_fft macros */
		/*	distance between fft data elements is 4 */
		/*	do 2 macros each processing 32 data values */

		rdi = (uintptr_t)g->u.xmm.sincos7;/* Load sin/cos pointer */
		loops_init(512);
		for(unsigned int loopA = 2; loopA; loopA--) {
			x4cl_eight_reals_fft_2(rsi, 2*64, 4*64, 8*64, rdi);
		}

		/* Do 2040 four_complex_fft macros */
		/*	distance between fft data elements is 4 */
		/*	do 510 macros each processing 32 data values */

		rsi += -2*2*64+16*64;	/* Load source pointer */
		rdi += XMM_SCD;		/* Next sine/cosine pointer */
		loops_reset;			/* 2 iters of 16 iters of 7 or 8 of 2 */
		for(unsigned int loopA = 2*16*8-1; loopA; ) {
			do{
				do{
					for(unsigned int loopD = 2; loopD; loopD--) {
						x4cl_four_complex_fft(rsi, 2*64, 4*64, 8*64, rdi);
					}
					rsi += -2*2*64+16*64;	/* Next source pointer */
					rdi += 2*XMM_SCD;		/* Next sine/cosine pointer */
				}while(--loopA & 7);
				rsi += -8*16*64+dist128;	/* Next source pointer */
			}while(loopA & (16*8-1));
			rsi += -16*dist128+64;	/* Next source pointer */
		}
		rsi += -2*64;		/* Restore source pointer */

		/* Do FFT levels 11,12 */
		/* Values 0-3 is real data, 4-7 is semi-real data, 8-16383 is */
		/* complex data. */
		/* */
		/* On input the 64-byte cache lines hold these data values: */
		/*	0	1	4	+1	8	+1	12	+1 */
		/*	+8192	... */
		/*	2	... */
		/*	+8192	... */
		/*	16	... */
		/*	... */
		/* On output the 64-byte cache lines hold these data values: */
		/*	0	1	2	+1	4	+1	6	+1 */
		/*	8192	... */
		/*	8	... */
		/*	... */

		/* Do 1/2 eight_reals_fft_2 macros */
		/* Do 1/2 nop_two_two_complex_fft_2 macros */
		/* Do 1 four_complex_fft macros */
		/*	distance between fft data elements is 1 */

		/* Do 2046 more four_complex_fft macros */
		/*	distance between fft data elements is 1 */
		/*	do 1023 macros each processing 16 data values */

		rdi = (uintptr_t)g->u.xmm.sincos6;/* Load sin/cos pointer */

		/* Execute the proper middle step */

		if likely(g->ffttype & 2){
			if likely(!(g->ffttype & 1)){
				s2cl_eight_reals_with_square_2(rsi, 4*64, 2*64);
				rdi += XMM_SCD;		/* Next sine/cosine pointer */
				loops_init_prefetch(1024, 64, 1, rcx, 1, 32);
				for(unsigned int loopA = 2*8*2*32-1; loopA; ) {
					do{
						xtouch(xptr(rcx+4096-128));		/* Preload the TLBs */
						do{
							do{
								s2cl_four_complex_with_square(rsi, 4*64, 2*64);
								rdi += XMM_SCD;		/* Next sine/cosine pointer */
							}while(--loopA & (32-1));
							rsi += -32*4*64+dist128;	/* Next source pointer */
						}while(loopA & (2*32-1));
					}while(loopA & (8*2*32-1));
					rsi += -16*dist128+64;	/* Next source pointer */
				}
				rsi += -2*64;		/* Restore source pointer */
			}else{
				rbp = g->DIST_TO_MULSRCARG;
				s2cl_eight_reals_with_mult_2(rsi, 4*64, 2*64);
				rdi += XMM_SCD;		/* Next sine/cosine pointer */
				loops_init_prefetch(1024, 64, 1, rcx, 1, 32);
				for(unsigned int loopA = 2*8*2*32-1; loopA; ) {
					do{
						xtouch(xptr(rcx+4096-128));		/* Preload the TLBs */
						do{
							do{
								s2cl_four_complex_with_mult(rsi, 4*64, 2*64);
								rdi += XMM_SCD;		/* Next sine/cosine pointer */
							}while(--loopA & (32-1));
							rsi += -32*4*64+dist128;	/* Next source pointer */
						}while(loopA & (2*32-1));
					}while(loopA & (8*2*32-1));
					rsi += -16*dist128+64;	/* Next source pointer */
				}
				rsi += -2*64;		/* Restore source pointer */
			}
		}else{
			s2cl_eight_reals_fft_2_final(rsi, 4*64, 2*64);
			rdi += XMM_SCD;		/* Next sine/cosine pointer */
			loops_init_prefetch(1024, 64, 1, rcx, 1, 32);
			for(unsigned int loopA = 2*8*2*32-1; loopA; ) {
				do{
					xtouch(xptr(rcx+4096-128));		/* Preload the TLBs */
					do{
						do{
							s2cl_four_complex_fft_final(rsi, 4*64, 2*64);
							rdi += XMM_SCD;		/* Next sine/cosine pointer */
						}while(--loopA & (32-1));
						rsi += -32*4*64+dist128;	/* Next source pointer */
					}while(loopA & (2*32-1));
				}while(loopA & (8*2*32-1));
				rsi += -16*dist128+64;	/* Next source pointer */
			}
			goto	xpass2_12_real_done;
		}
	}else{
		rbp = g->DIST_TO_MULSRCARG;
		rdi = (uintptr_t)g->u.xmm.sincos6;/* Load sin/cos pointer */
		xtouch_init(rcx = (uintptr_t)g->data_prefetch);/* Load prefetch pointer */
		s2cl_eight_reals_with_mulf_2(rsi, 4*64, 2*64);
		rdi += XMM_SCD;		/* Next sine/cosine pointer */
		loops_init_prefetch(1024, 64, 1, rcx, 1, 32);
		xtouch(xptr(rcx));			/* Preload the TLBs */
		for(unsigned int loopA = 2*8*2*32-1; loopA; ) {
			do{
				xtouch(xptr(rcx+4096-128));		/* Preload the TLBs */
				do{
					do{
						s2cl_four_complex_with_mulf(rsi, 4*64, 2*64);
						rdi += XMM_SCD;		/* Next sine/cosine pointer */
					}while(--loopA & (32-1));
					rsi += -32*4*64+dist128;	/* Next source pointer */
				}while(loopA & (2*32-1));
			}while(loopA & (8*2*32-1));
			rsi += -16*dist128+64;	/* Next source pointer */
		}
		rsi += -2*64;		/* Restore source pointer */

		/* Do inverse FFT levels 9,10 */
		/* On input the 64-byte cache lines hold these data values: */
		/*	0	1	2	+1	4	+1	6	+1 */
		/*	8192	... */
		/*	8	... */
		/*	... */
		/* On output the 64-byte cache lines hold these data values: */
		/*	0	1	8	+1	16	+1	24	+1 */
		/*	+8192	... */
		/*	2	... */
		/*	... */
		/*	6	... */
		/*	+8192	... */
		/*	32	... */
		/*	... */

		/* Do 2 eight_reals_unfft_2 macros */
		/* Do 2 nop_two_two_complex_unfft_2 macros */
		/*	distance between fft data elements is 4 */
		/*	do 1 macro processing 32 data values */

	}
	x4cl_eight_reals_unfft_2(rsi, 8*64, 2*64, 4*64);

	/* Do 2044 four_complex_unfft macros */
	/*	distance between fft data elements is 4 */
	/*	do 511 macros each processing 32 data values */

	rdi = (uintptr_t)g->u.xmm.sincos7;/* Load sin/cos pointer */
	loops_init(512, 1, 16);		/* 2 iters of 16 iters of 15 or 16 */
	for(unsigned int loopA = 2*16*16-1; loopA; ) {
		do{
			do{
				x4cl_four_complex_unfft(rsi, 8*64, 2*64, 4*64, rdi);
				rdi += XMM_SCD;		/* Next sine/cosine pointer */
			}while(--loopA & 15);
			rsi += -16*8*64+dist128;	/* Next source pointer */
		}while(loopA & (16*16-1));
		rsi += -16*dist128+64;	/* Next source pointer */
	}
	rsi += -2*64;		/* Restore source pointer */

	/* Do inverse FFT levels 7,8 */
	/* On input the 64-byte cache lines hold these data values: */
	/*	0	1	8	+1	16	+1	24	+1 */
	/*	+8192	... */
	/*	2	... */
	/*	... */
	/*	6	... */
	/*	+8192	... */
	/*	32	... */
	/*	... */
	/* On output the 64-byte cache lines hold these data values: */
	/*	0	1	32	+1	64	+1	96	+1 */
	/*	+8192	... */
	/*	2	... */
	/*	... */
	/*	30	... */
	/*	+8192	... */
	/*	128	... */
	/*	... */

	/* Do 8 eight_reals_unfft_2 macros */
	/* Do 8 nop_two_two_complex_unfft_2 macros */
	/*	distance between fft data elements is 16 */
	/*	do 4 macros each processing 32 data values */

	loops_init(512);
	for(unsigned int loopA = 4; loopA; loopA--) {
		x4cl_eight_reals_unfft_2(rsi, 2*64, 8*64, 16*64);
	}

	/* Do 2032 four_complex_unfft macros */
	/*	distance between fft data elements is 16 */
	/*	do 508 macros each processing 32 data values */

	rdi = (uintptr_t)g->u.xmm.sincos8;/* Load sin/cos pointer */
	rsi += -4*2*64+32*64;	/* Load source pointer */
	loops_reset;			/* 2 iters of 16 of 3 or 4 of 4 */
	for(unsigned int loopA = 2*16*4-1; loopA; ) {
		do{
			do{
				for(unsigned int loopD = 4; loopD; loopD--) {
					x4cl_four_complex_unfft(rsi, 2*64, 8*64, 16*64, rdi);
				}
				rsi += -4*2*64+32*64;	/* Next source pointer */
				rdi += XMM_SCD;		/* Next sine/cosine pointer */
			}while(--loopA & 3);
			rsi += -4*32*64+dist128;	/* Next source pointer */
		}while(loopA & (16*4-1));
		rsi += -16*dist128+64;	/* Next source pointer */
	}
	rsi += -2*64;		/* Restore source pointer */

	/* Do inverse FFT levels 5,6 */
	/* On input the 64-byte cache lines hold these data values: */
	/*	0	1	32	+1	64	+1	96	+1 */
	/*	+8192	... */
	/*	2	... */
	/*	... */
	/*	30	... */
	/*	+8192	... */
	/*	128	... */
	/*	... */
	/* On output the 64-byte cache lines hold these data values: */
	/*	0	1	128	+1	256	+1	384	+1 */
	/*	+8192	... */
	/*	2	... */
	/*	... */
	/*	126	... */
	/*	+8192	... */
	/*	512	... */
	/*	... */

	/* Do 32 eight_reals_unfft_2 macros */
	/* Do 32 nop_two_two_complex_unfft_2 macros */
	/*	distance between fft data elements is 64 */
	/*	do 16 macros each processing 32 data values */

	loops_init(512);
	for(unsigned int loopA = 16; loopA; loopA--) {
		x4cl_eight_reals_unfft_2(rsi, 2*64, 32*64, 64*64);
	}

	/* Do 1984 four_complex_unfft macros */
	/*	distance between fft data elements is 16 */
	/*	do 496 macros each processing 32 data values */

	rdi = (uintptr_t)g->u.xmm.sincos9;/* Load sin/cos pointer */
	rsi += -16*2*64+dist128;	/* Load source pointer */
	loops_reset;			/* 2 iters of 15 or 16 iters of 16 */
	for(unsigned int loopA = 2*16-1; loopA; ) {
		do{
			for(unsigned int loopC = 16; loopC; loopC--) {
				x4cl_four_complex_unfft(rsi, 2*64, 32*64, 64*64, rdi);
			}
			rsi += -16*2*64+dist128;	/* Next source pointer */
			rdi += XMM_SCD;		/* Next sine/cosine pointer */
		}while(--loopA & 15);
		rsi += -16*dist128+64;	/* Next source pointer */
	}
	rsi += -2*64;		/* Restore source pointer */

	/* Do inverse FFT levels 4 */
	/* On input the 64-byte cache lines hold these data values: */
	/*	0	1	128	+1	256	+1	384	+1 */
	/*	+8192	... */
	/*	2	... */
	/*	... */
	/*	126	... */
	/*	+8192	... */
	/*	512	... */
	/*	... */
	/* On output the 64-byte cache lines hold these data values: */
	/*	0	1	256	+1	512	+1	768	+1 */
	/*	+8192	... */
	/*	2	... */
	/*	... */
	/*	254	... */
	/*	+8192	... */
	/*	1024	... */
	/*	... */

	/* Do 128 eight_reals_unfft_1 macros */
	/*	distance between fft data elements is 256 */
	/*	do 64 macros each processing 16 data values */

	loops_init(1024);
	for(unsigned int loopA = 64; loopA; loopA--) {
		x2cl_eight_reals_unfft_1(rsi, 2*64, dist128);
	}

	/* Do 1920 two_complex_unfft macros */
	/*	distance between fft data elements is 256 */
	/*	do 960 macros each processing 16 data values */

	rdi = (uintptr_t)g->u.xmm.sincos10;/* Load sin/cos pointer */
	rsi += -64*2*64+2*dist128;	/* Load source pointer */
	loops_reset;			/* 7 then 8 iterations of 64 */
	for(unsigned int loopA = 2*8-1; loopA; ) {
		do{
			for(unsigned int loopC = 64; loopC; loopC--) {
				x2cl_two_complex_unfft(rsi, 2*64, dist128);
			}
			rsi += -64*2*64+2*dist128;	/* Next source pointer */
			rdi += XMM_SCD;		/* Next sine/cosine pointer */
		}while(--loopA & 7);
		rsi += -8*2*dist128+64;	/* Next source pointer */
	}
	rsi += -2*64;		/* Restore source pointer */

	/* Do inverse FFT levels 2,3 */
	/* On input the 64-byte cache lines hold these data values: */
	/*	0	1	256	+1	512	+1	768	+1 */
	/*	+8192	... */
	/*	2	... */
	/*	... */
	/*	254	... */
	/*	+8192	... */
	/*	1024	... */
	/*	... */
	/* On output the 64-byte cache lines hold these data values: */
	/*	0	1	1024	+1	2048	+1	3072	+1 */
	/*	+8192	... */
	/*	2	... */
	/*	... */
	/*	510	... */
	/*	+8192	... */
	/*	4096	... */
	/*	... */

	/* Do 256 eight_reals_unfft_2 macros */
	/* Do 256 nop_two_two_complex_unfft_2 macros */
	/*	distance between fft data elements is 512 */
	/*	do 128 macros each processing 32 data values */

	loops_init(512);
	for(unsigned int loopA = 2; loopA; loopA--) {
		for(unsigned int loopB = 64; loopB; loopB--) {
			x4cl_eight_reals_unfft_2(rsi, 2*64, 2*dist128, 4*dist128);
		}
		rsi += -64*2*64+dist128;	/* Next source pointer */
	}

	/* Do 1536 four_complex_unfft macros */
	/*	distance between fft data elements is 512 */
	/*	do 384 macros each processing 32 data values */

	rsi += -2*dist128+8*dist128;/* Next source pointer */
	rdi = (uintptr_t)g->u.xmm.sincos11;/* Load sin/cos pointer */
	loops_reset;			/* 2 iters of 1 or 2 iters of 2 of 64 */
	for(unsigned int loopA = 2*2-1; loopA; ) {
		do{
			for(unsigned int loopC = 2; loopC; loopC--) {
				for(unsigned int loopD = 64; loopD; loopD--) {
					x4cl_four_complex_unfft(rsi, 2*64, 2*dist128, 4*dist128, rdi);
				}
				rsi += -64*2*64+dist128;	/* Next source pointer */
			}
			rsi += -2*dist128+8*dist128;/* Next source pointer */
			rdi += XMM_SCD;		/* Next sine/cosine pointer */
		}while(--loopA & 1);
		rsi += -2*8*dist128+64;	/* Next source pointer */
	}
	rsi += -2*64;		/* Restore source pointer */

	/* Do inverse FFT level 1 */
	/* On input the 64-byte cache lines hold these data values: */
	/*	0	1	1024	+1	2048	+1	3072	+1 */
	/*	+8192	... */
	/*	2	... */
	/*	... */
	/*	510	... */
	/*	+8192	... */
	/*	4096	... */
	/*	... */
	/* On output the 64-byte cache lines hold these data values: */
	/*	0	+2K	4096	+2K	8192	+2K	12288	+2K */
	/*	1	... */
	/*	... */
	/*	2047	... */

	/* Do 1024 eight_reals_last_unfft macros */
	/* Do 1024 two_two_complex_unfft macros */
	/*	distance between fft data elements is 2048 */
	/*	do 512 macros each processing 16 data values */

	loops_init(512);
	for(unsigned int loopA = 8; loopA; loopA--) {
		for(unsigned int loopB = 64; loopB; loopB--) {
			s4cl_eight_reals_unfft_1(rsi, 2*64, 8*dist128, 64);
		}
		rsi += -64*2*64+dist128;	/* Next source pointer */
	}

	xpass2_12_real_done:;
	if(g->pass2_get_next_block(g) != 0) return;
}
hg_pass2_12_levels_complex(g);
}


/* *************** xmm-pass2-12-levels-complex macro ****************** */
/* This macro takes 2 groups of 4096 complex values, premultiplies them and */
/* performs the final 12 levels of the FFT process, squares the */
/* results, and does the first 12 levels of the inverse FFT. */
/* NOTE: Rather than remembering 4096 multipliers we break them up */
/* into 1024 groups of 4 columns.  By remembering these 1028 values we can */
/* compute all 4096 multipliers using a single complex multiply.  These extra */
/* complex multiplies can probably be done faster than the computer can read the */
/* extra multipliers from memory (and are more likely stay in the L2 cache) */

void hg_pass2_12_levels_complex(struct gwasm_data *__restrict g) {
do{
	uintptr_t rbx,rbp,rsi,rdi;
	vec2f64 xmm0,xmm1,xmm2,xmm3,xmm4,xmm5,xmm6,xmm7;

	/* Type 4 FFTs skip the forward FFT process */

	rsi = (uintptr_t)g->data_addr;/* Load address of FFT data */
	rbx = g->DIST_TO_FFTSRCARG;
	if likely(!(g->ffttype & 4)){

		/* Do FFT level 1,2 */
		/* */
		/* We swizzle the data so that one set of data is in the high word of */
		/* an xmm register, the other set is in the low word.  This means */
		/* a distance of 8192 data values. */
		/* */
		/* On input the 64-byte cache lines hold these data values: */
		/*	0	+2K	4096	+2K	8192	+2K	12288	+2K */
		/*	1	... */
		/*	... */
		/*	2047 */
		/* On output the 64-byte cache lines hold these data values: */
		/*	0	+8K	1024	+8K	2048	+8K	3072	+8K */
		/*	1	... */
		/*	... */
		/*	1023	... */
		/*	4096	... */
		/*	... */

		/* Do 2048 four_complex_gpm_fft macros */
		/*	distance between fft data elements is 1024 */
		/*	do 1024 macros each processing 16 data values */

		start_timer(4);
		rdi = (uintptr_t)g->premult_addr;/* The group multipliers */
		loops_init(1024);			/* 256 iterations of 4 */
		/*	mov	rcx, sincos_complex	;; Prefetch pointer */
		for(unsigned int loopA = 8; loopA; loopA--) {
			for(unsigned int loopB = 8; loopB; loopB--) {
				for(unsigned int loopC = 16; loopC; loopC--) {
					s2cl_four_complex_gpm_fft(rsi, 64, 8*dist128);
					L2prefetch128(u8ptr(rsi+5*64));		/* Prefetch a few iterations ahead in */
					L2prefetch128(u8ptr(rsi+8*dist128+5*64));/* case they are not in the L2 cache */
				}
				rdi += XMM_PMD;		/* Next group multiplier */
				L2prefetch128(u8ptr(rdi+2*XMM_PMD));
				/*	xtouch	[rcx]			;; Touch the sin/cos data to keep it */
				/*	bump	rcx, 96			;; in the L2 data cache */
			}
			rsi += -8*16*64+dist128;	/* Next source pointer */
		}
		rsi += -8*dist128;		/* Restore source pointer */
		end_timer(4);

		/* Do FFT levels 3,4 */
		/* */
		/* On input the 64-byte cache lines hold these data values: */
		/*	0	+8K	1024	+8K	2048	+8K	3072	+8K */
		/*	1	... */
		/*	... */
		/*	1023	... */
		/*	4096	... */
		/*	... */
		/* On output the 64-byte cache lines hold these data values: */
		/*	0	+8K	256	+8K	512	+8K	768	+8K */
		/*	1	... */
		/*	... */
		/*	255	... */
		/*	1024	... */
		/*	... */

		/* Do 4 groups of four_complex_fft_cpm macros */
		/* To make things run faster we precompute the column multipliers times the */
		/* 4 different sine/cosine multipliers we run into in this section.  This */
		/* costs some memory but saves us some complex multiplies. */

		/* Do 2048 four_complex_cpm_fft macros */
		/*	distance between fft data elements is 256 */
		/*	do 256 macros each processing 64 data values */

		start_timer(5);
		xtouch_init(rcx = (uintptr_t)g->premult_prefetch);/* Prefetch pointer */
		loops_init_prefetch(512, 64, 1, rcx); /* 2 iters of 2 iters of 8 iters of 16 */
		xtouch(xptr(rcx));			/* Load first TLB */
		for(unsigned int loopA = 2; loopA; loopA--) {
			xtouch(xptr(rcx+4096-128));		/* Load prefetch TLBs */
			xtouch(xptr(rcx+2*4096-128));
			for(unsigned int loopB = 2; loopB; loopB--) {
				for(unsigned int loopC = 8; loopC; loopC--) {
					for(unsigned int loopD = 16; loopD; loopD--) {
						x4cl_four_complex_cpm_fft(rsi, 64, 2*dist128, 4*dist128, 16*XMM_PMD);
						rdi += XMM_PMD;		/* Next column premultiplier */
					}
					rdi += -16*XMM_PMD;	/* Restore column premultiplier */
				}
				rsi += -8*16*64+dist128;	/* Next source pointer */
				xtouch_init(rcx += -64*128);		/* Re-prefetch section */
			}
			rsi += -2*dist128+8*dist128;/* Next source pointer */
			rdi += 32*XMM_PMD;		/* Next set of premultipliers */
			xtouch_init(rcx += 64*128);		/* Prefetch next section */
		}
		rsi += -2*8*dist128;	/* Restore source pointer */
		end_timer(5);

		/* Do FFT levels 5,6 */
		/* */
		/* On input the 64-byte cache lines hold these data values: */
		/*	0	+8K	256	+8K	512	+8K	768	+8K */
		/*	1	... */
		/*	... */
		/*	255	... */
		/*	1024	... */
		/*	... */
		/* On output the 64-byte cache lines hold these data values: */
		/*	0	+8K	64	+8K	128	+8K	192	+8K */
		/*	1	... */
		/*	... */
		/*	63	... */
		/*	256	... */
		/*	... */
		/* */
		/* Do 2048 four_complex_fft macros */
		/*	distance between fft data elements is 64 */
		/*	do 512 macros each processing 32 data values */

		/* Do FFT levels 7,8 */
		/* */
		/* On input the 64-byte cache lines hold these data values: */
		/*	0	+8K	64	+8K	128	+8K	192	+8K */
		/*	1	... */
		/*	... */
		/*	63	... */
		/*	256	... */
		/*	... */
		/* On output the 64-byte cache lines hold these data values: */
		/*	0	+8K	16	+8K	32	+8K	48	+8K */
		/*	1	... */
		/*	... */
		/*	15	... */
		/*	64	... */
		/*	... */
		/* */
		/* Do 2048 four_complex_fft macros */
		/*	distance between fft data elements is 16 */
		/*	do 512 macros each processing 32 data values */

		/* Do FFT levels 9,10 */
		/* */
		/* On input the 64-byte cache lines hold these data values: */
		/*	0	+8K	16	+8K	32	+8K	48	+8K */
		/*	1	... */
		/*	... */
		/*	15	... */
		/*	64	... */
		/*	... */
		/* On output the 64-byte cache lines hold these data values: */
		/*	0	+8K	4	+8K	8	+8K	12	+8K */
		/*	1	... */
		/*	... */
		/*	3	... */
		/*	16	... */
		/*	... */
		/* */
		/* Do 2048 four_complex_fft macros */
		/*	distance between fft data elements is 4 */
		/*	do 512 macros each processing 32 data values */

		/* Do this 64 cache lines at a time to maximize accessing data */
		/* while it is in the L1 cache. */

		start_timer(6);
		xtouch_init(rcx = (uintptr_t)g->data_prefetch);/* Prefetch pointer */
		rdi = (uintptr_t)g->xsincos_complex;	/* Load sin/cos pointers */
		rbx = (uintptr_t)g->xsincos_complex;
		rbp = (uintptr_t)g->xsincos_complex;
		loops_init_prefetch(512, 64, 1, rcx); /* 2*2*2 iters of 4 of 64 cache lines */

		xtouch(xptr(rcx));			/* Load first TLB */
		for(unsigned int loopA = 2*2*2; loopA; ) {
			xtouch(xptr(rcx+4096-128));		/* Load six prefetch TLBs */
			xtouch(xptr(rcx+2*4096-128));
			xtouch(xptr(rcx+3*4096-128));
			xtouch(xptr(rcx+4*4096-128));
			xtouch(xptr(rcx+5*4096-128));
			xtouch(xptr(rcx+6*4096-128));

			do{
				do{
					for(unsigned int loopD = 4; loopD; loopD--) {
						for(unsigned int loopE = 16; loopE; loopE--) {
							x4cl_four_complex_fft(rsi, 4*64, 64*64, dist128, rdi);
						}
						rsi += -16*4*64;		/* Restore source pointer */
						xtouch_init(rcx += -8*128+8192+128);	/* Prefetch from next section */
						loops_undo(16);

						for(unsigned int loopE = 2; loopE; loopE--) {
							for(unsigned int loopF = 2; loopF; loopF--) {
								for(unsigned int loopG = 4; loopG; loopG--) {
									x4cl_four_complex_fft(rsi, 4*64, 16*64, 32*64, rbx);
								}
								rsi += -4*4*64+64*64;	/* Next source pointer */
								rbx += 2*XMM_SCD;		/* Next sine/cosine pointer */
							}
							rsi += -2*64*64+dist128;	/* Next source pointer */
						}
						rsi += -2*dist128;		/* Restore source pointer */
						rbx += -4*2*XMM_SCD;	/* Restore sine/cosine pointer */
						xtouch_init(rcx += -8*128+8192+128);	/* Prefetch from next section */
						loops_undo(2*2*4);

						for(unsigned int loopE = 2; loopE; loopE--) {
							for(unsigned int loopF = 8; loopF; loopF--) {
								x4cl_four_complex_fft(rsi, 16*64, 4*64, 8*64, rbp);
								rbp += 2*XMM_SCD;		/* Next sine/cosine pointer */
							}
							rsi += -8*16*64+dist128;	/* Next source pointer */
						}
						rbp += -16*2*XMM_SCD;	/* Restore sine/cosine pointer */
						xtouch_init(rcx += -2*(8192+128));	/* Back to prefetching from 1st section */

						rsi += -2*dist128+64;	/* Next source pointer */
					}
					rsi += -4*64+2*dist128;	/* Next source pointer */
					rdi += 2*XMM_SCD;		/* Next sine/cosine pointer */
					rbx += 4*2*XMM_SCD;	/* Next sine/cosine pointer */
					rbp += 16*2*XMM_SCD;	/* Next sine/cosine pointer */
				}while(--loopA & (2-1));
				xtouch_init(rcx += -64*128);		/* Re-prefetch 3 sections */
			}while(loopA & (2*2-1));
			xtouch_init(rcx += 3*(8192+128));	/* Prefetch next 3 sections */
		}
		rsi += -8*2*dist128;	/* Restore source pointer */
		end_timer(6);

		/* Do FFT levels 11,12 as well as inverse FFT levels 11,12 */
		/* */
		/* On input the 64-byte cache lines hold these data values: */
		/*	0	+8K	4	+8K	8	+8K	12	+8K */
		/*	1	... */
		/*	... */
		/*	3	... */
		/*	16	... */
		/*	... */
		/* On output the 64-byte cache lines hold these data values: */
		/*	0	+8K	2	+8K	4	+8K	6	+8K */
		/*	1	... */
		/*	8	... */
		/*	... */

		/* Do 2048 four_complex_fft macros */
		/*	distance between fft data elements is 1 */
		/*	do 512 macros each processing 32 data values */

		/* Execute the right middle step */

		if likely(g->ffttype & 2){
			if likely(!(g->ffttype & 1)){
				start_timer(9);
				loops_init_prefetch(512, 64, 1, rcx);
				for(unsigned int loopA = 2*2; loopA; ) {
					xtouch(xptr(rcx+4096-128));		/* Load prefetch TLBs */
					xtouch(xptr(rcx+2*4096-128));
					do{
						rdi = (uintptr_t)g->xsincos_complex;	/* Load sin/cos pointer */
						for(unsigned int loopC = 4; loopC; loopC--) {
							for(unsigned int loopD = 32; loopD; loopD--) {
								x4cl_four_complex_with_square(rsi, 4*64, 64, 2*64);
								rdi += 2*XMM_SCD;		/* Next sine/cosine pointer */
							}
							rsi += -16*2*4*64+dist128;	/* Next source pointer */
						}
						xtouch_init(rcx += -64*128);		/* Re-prefetch section */
					}while(--loopA & 1);
					xtouch_init(rcx += (8192+128));		/* Prefetch next section */
				}
				rsi += -16*dist128;	/* Restore source pointer */
				end_timer(9);
			}else{
				start_timer(9);
				loops_init_prefetch(512, 64, 1, rcx);
				rbp = g->DIST_TO_MULSRCARG;
				for(unsigned int loopA = 2*2; loopA; ) {
					xtouch(xptr(rcx+4096-128));		/* Load prefetch TLBs */
					xtouch(xptr(rcx+2*4096-128));
					do{
						rdi = (uintptr_t)g->xsincos_complex;	/* Load sin/cos pointer */
						for(unsigned int loopC = 4; loopC; loopC--) {
							for(unsigned int loopD = 32; loopD; loopD--) {
								x4cl_four_complex_with_mult(rsi, 4*64, 64, 2*64);
								rdi += 2*XMM_SCD;		/* Next sine/cosine pointer */
							}
							rsi += -16*2*4*64+dist128;	/* Next source pointer */
						}
						xtouch_init(rcx += -64*128);		/* Re-prefetch section */
					}while(--loopA & 1);
					xtouch_init(rcx += (8192+128));		/* Prefetch next section */
				}
				rsi += -16*dist128;	/* Restore source pointer */
				end_timer(9);
			}
		}else{
			start_timer(9);
			loops_init_prefetch(512, 64, 1, rcx);
			for(unsigned int loopA = 4; loopA; loopA--) {
				rdi = (uintptr_t)g->xsincos_complex;	/* Load sin/cos pointer */
				xtouch(xptr(rcx+4096-128));		/* Load prefetch TLBs */
				xtouch(xptr(rcx+2*4096-128));
				for(unsigned int loopB = 4; loopB; loopB--) {
					for(unsigned int loopC = 32; loopC; loopC--) {
						x4cl_four_complex_fft_final(rsi, 4*64, 64, 2*64);
						rdi += 2*XMM_SCD;		/* Next sine/cosine pointer */
					}
					rsi += -16*2*4*64+dist128;	/* Next source pointer */
				}
				xtouch_init(rcx += 128);		/* Prefetch next section */
			}
			end_timer(9);
			goto	xpass2_12_levels_complex_done;
		}
	}else{
		start_timer(9);
		xtouch_init(rcx = (uintptr_t)g->data_prefetch);/* Prefetch pointer */
		rbp = g->DIST_TO_MULSRCARG;
		loops_init_prefetch(512, 64, 1, rcx);
		xtouch(xptr(rcx));			/* Load first prefetch TLB */
		for(unsigned int loopA = 4; loopA; loopA--) {
			rdi = (uintptr_t)g->xsincos_complex;	/* Load sin/cos pointer */
			xtouch(xptr(rcx+4096-128));		/* Load prefetch TLBs */
			xtouch(xptr(rcx+2*4096-128));
			for(unsigned int loopB = 4; loopB; loopB--) {
				for(unsigned int loopC = 32; loopC; loopC--) {
					x4cl_four_complex_with_mulf(rsi, 4*64, 64, 2*64);
					rdi += 2*XMM_SCD;		/* Next sine/cosine pointer */
				}
				rsi += -16*2*4*64+dist128;	/* Next source pointer */
			}
			xtouch_init(rcx += 128);		/* Prefetch next section */
		}
		rsi += -16*dist128;	/* Restore source pointer */
		end_timer(9);

		/* Do inverse FFT levels 9,10 */
		/* On input the 64-byte cache lines hold these data values: */
		/*	0	+8K	2	+8K	4	+8K	6	+8K */
		/*	1	... */
		/*	8	... */
		/*	... */
		/* On output the 64-byte cache lines hold these data values: */
		/*	0	+8K	8	+8K	16	+8K	24	+8K */
		/*	1	... */
		/*	... */
		/*	7	... */
		/*	32	... */
		/*	... */
		/* */
		/* Do 2048 four_complex_unfft macros */
		/*	distance between fft data elements is 4 */
		/*	do 512 macros each processing 32 data values */

		/* Do inverse FFT levels 7,8 */
		/* On input the 64-byte cache lines hold these data values: */
		/*	0	+8K	8	+8K	16	+8K	24	+8K */
		/*	1	... */
		/*	... */
		/*	7	... */
		/*	32	... */
		/*	... */
		/* On output the 64-byte cache lines hold these data values: */
		/*	0	+8K	32	+8K	64	+8K	96	+8K */
		/*	1	... */
		/*	... */
		/*	31	... */
		/*	128	... */
		/*	... */
		/* */
		/* Do 2048 four_complex_unfft macros */
		/*	distance between fft data elements is 16 */
		/*	do 512 macros each processing 32 data values */

		/* Do inverse FFT levels 5,6 */
		/* On input the 64-byte cache lines hold these data values: */
		/*	0	+8K	32	+8K	64	+8K	96	+8K */
		/*	1	... */
		/*	... */
		/*	31	... */
		/*	128	... */
		/*	... */
		/* On output the 64-byte cache lines hold these data values: */
		/*	0	+8K	128	+8K	256	+8K	384	+8K */
		/*	1	... */
		/*	... */
		/*	127	... */
		/*	512	... */
		/*	... */
		/* */
		/* Do 2048 four_complex_unfft macros */
		/*	distance between fft data elements is 64 */
		/*	do 512 macros each processing 32 data values */

		/* Do this 64 cache lines at a time to maximize accessing data */
		/* while it is in the L1 cache. */

	}
	start_timer(10);
	rdi = (uintptr_t)g->xsincos_complex;	/* Load sin/cos pointers */
	rbx = (uintptr_t)g->xsincos_complex;
	rbp = (uintptr_t)g->xsincos_complex;
	loops_init_prefetch(512, 64, 1, rcx); /* 2*2*4*2 iters of 64 cache lines */

	for(unsigned int loopA = 2*2*4; loopA; ) {
		xtouch(xptr(rcx+4096-128));		/* Load six prefetch TLBs */
		xtouch(xptr(rcx+2*4096-128));
		xtouch(xptr(rcx+3*4096-128));
		xtouch(xptr(rcx+4*4096-128));
		xtouch(xptr(rcx+5*4096-128));
		xtouch(xptr(rcx+6*4096-128));

		do{
			do{
				for(unsigned int loopD = 2; loopD; loopD--) {
					for(unsigned int loopE = 16; loopE; loopE--) {
						x4cl_four_complex_unfft(rsi, 8*64, 2*64, 4*64, rbp);
						rbp += XMM_SCD;		/* Next sine/cosine pointer */
					}
					rsi += -16*8*64;		/* Restore source pointer */
					rbp += -16*XMM_SCD;	/* Next sine/cosine pointer */
					xtouch_init(rcx += -8*128+8192+128);	/* Prefetch from next section */
					loops_undo(16);

					for(unsigned int loopE = 4; loopE; loopE--) {
						for(unsigned int loopF = 4; loopF; loopF--) {
							x4cl_four_complex_unfft(rsi, 2*64, 8*64, 16*64, rbx);
						}
						rsi += -4*2*64+32*64;	/* Next source pointer */
						rbx += XMM_SCD;		/* Next sine/cosine pointer */
					}
					rsi += -4*32*64;		/* Restore source pointer */
					rbx += -4*XMM_SCD;		/* Restore sine/cosine pointer */
					xtouch_init(rcx += -8*128+8192+128);	/* Prefetch from next section */
					loops_undo(4*4);

					for(unsigned int loopE = 16; loopE; loopE--) {
						x4cl_four_complex_unfft(rsi, 2*64, 32*64, 64*64, rdi);
					}
					xtouch_init(rcx += -2*(8192+128));	/* Back to prefetching from 1st section */

					rsi += -16*2*64+64;	/* Next source pointer */
				}
				rsi += -2*64+dist128;	/* Next source pointer */
				rdi += XMM_SCD;		/* Next sine/cosine pointer */
				rbx += 4*XMM_SCD;		/* Next sine/cosine pointer */
				rbp += 16*XMM_SCD;		/* Next sine/cosine pointer */
			}while(--loopA & (4-1));
			xtouch_init(rcx += -64*128);		/* Re-prefetch 3 sections */
		}while(loopA & (2*4-1));
		xtouch_init(rcx += 3*(8192+128));	/* Prefetch next 3 sections */
	}
	rsi += -16*dist128;	/* Restore source pointer */
	end_timer(10);

	/* Do inverse FFT levels 3,4 */
	/* On input the 64-byte cache lines hold these data values: */
	/*	0	+8K	128	+8K	256	+8K	384	+8K */
	/*	1	... */
	/*	... */
	/*	127	... */
	/*	512	... */
	/*	... */
	/* On output the 64-byte cache lines hold these data values: */
	/*	0	+8K	512	+8K	1024	+8K	1536	+8K */
	/*	1	... */
	/*	... */
	/*	511	... */
	/*	2048	... */
	/*	... */

	/* Do 4 groups of 512 four_complex_fft macros (FFT levels 3,4) */
	/* In this section we also apply the 4 column multipliers (actually three */
	/* since the first column multiplier is always one).  To make things */
	/* run faster we precompute the column multipliers times the 4 different */
	/* sine/cosine multipliers we run into in this section.  This costs some */
	/* memory but saves us some complex multiplies. */

	/* Do 2048 four_complex_cpm_unfft macros */
	/*	distance between fft data elements is 256 */
	/*	do 128 macros each processing 128 data values */

	start_timer(13);
	rdi = (uintptr_t)g->premult_addr;/* Load column multipliers pointer */
	rdi += 64*XMM_PMD;		/* Load column multipliers pointer */
	loops_init_prefetch(512, 64, 1, rcx); /* 2 iters of 2 iters of 8 iters of 16 */
	for(unsigned int loopA = 2*2; loopA; ) {
		xtouch(xptr(rcx+4096-128));		/* Load prefetch TLBs */
		xtouch(xptr(rcx+2*4096-128));
		do{
			for(unsigned int loopC = 8; loopC; loopC--) {
				for(unsigned int loopD = 16; loopD; loopD--) {
					x4cl_four_complex_cpm_unfft(rsi, 64, dist128, 2*dist128);
					rdi += XMM_PMD;		/* Next column premultiplier */
				}
				rdi += -16*XMM_PMD;	/* Reset premultiplier pointer */
			}
			rsi += -128*64+4*dist128;	/* Next source pointer */
			rdi += 16*XMM_PMD;		/* Next premultiplier pointer */
			xtouch_init(rcx += -64*128);		/* Reprefetch 8KB section */
		}while(--loopA & 1);
		xtouch_init(rcx += 8192+128);		/* Next 8KB section to prefetch */
	}
	rsi += -4*4*dist128;	/* Restore source pointer */
	end_timer(13);

	/* Do inverse FFT levels 1,2 */
	/* On input the 64-byte cache lines hold these data values: */
	/*	0	+8K	512	+8K	1024	+8K	1536	+8K */
	/*	1	... */
	/*	... */
	/*	511	... */
	/*	2048	... */
	/*	... */
	/* On output the 64-byte cache lines hold these data values: */
	/*	0	+2K	4096	+2K	8192	+2K	12288	+2K */
	/*	1	... */
	/*	... */
	/*	2047 */

	/* Do 512 four_complex_gpm_unfft macros */
	/*	distance between fft data elements is 1024 */
	/*	do 512 macros each processing 32 data values */

	start_timer(14);
	xtouch_init(rcx = (uintptr_t)g->premult_prefetch);/* Prefetch pointer */
	rdi = (uintptr_t)g->premult_addr;/* Load group multiplier pointer */
	loops_init_prefetch(512, 64, 2, rcx); /* 4 iters of 8 iters of 16 */
	xtouch(xptr(rcx));			/* Load first prefetch TLB */
	for(unsigned int loopA = 4; loopA; loopA--) {
		xtouch(xptr(rcx+4096-128));		/* Load prefetch TLBs */
		for(unsigned int loopB = 8; loopB; loopB--) {
			for(unsigned int loopC = 16; loopC; loopC--) {
				s4cl_four_complex_gpm_unfft(rsi, 64, 4*dist128, 8*dist128, 32*XMM_PMD);
			}
			rdi += XMM_PMD;		/* Next set of multipliers */
		}
		rsi += -8*4*4*64+dist128;	/* Next source pointer */
	}
	end_timer(14);

	xpass2_12_levels_complex_done:;

}while(g->pass2_get_next_block(g) == 0);
}
