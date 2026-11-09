#include <stdlib.h>
#include <ctype.h>

#include "../gwnum.h"
#include "../gwcommon.h"
#include "../unravel.h"
#include "../gwtables.h"
#include "../xarch.h"
#include "../xbasics.h"
#include "../xmult.h"
#include "../gwnum_func.h"
#include "hg.h"

/* *************** xmm-pass2-11-levels-real macro ****************** */
/* This macro takes 2048 real values, 2048 semi-real values and the */
/* first 2048 complex values and performs the final 11 levels of the */
/* FFT process, squares the results, and does the first 11 levels of */
/* the inverse FFT. */

void hg_pass2_11_levels_complex(struct gwasm_data *__restrict g);
void xpass2_hg_11_levels(struct gwasm_data *__restrict g) {
g->thread_work_routine = hg_pass2_11_levels_complex;
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
		/*	0	+1K	2048	+1K	4096	+1K	6144	+1K */
		/*	1	... */
		/*	... */
		/*	1023	... */
		/* On output the 64-byte cache lines hold these data values: */
		/*	0	1	1024	+1	2048	+1	3072	+1 */
		/*	4096	... */
		/*	2	... */
		/*	... */
		/*	1022	... */
		/*	5118	... */

		/* Do eight_reals_fft_1 on real values 0 - 2047 */
		/* Do nothing on semi-real values from 2048 - 4095 */
		/* Do two_complex_fft on complex values from 4096 - 8191 */
		/*	distance between fft data elements is 1K */
		/*	do 512 iterations */

		loops_init(512);
		for(unsigned int loopA = 8; loopA; loopA--) {
			for(unsigned int loopB = 64; loopB; loopB--) {
				s2cl_eight_reals_fft_1(rsi, 2*64, 64);
				L2prefetch128(u8ptr(rsi+3*128));		/* Prefetch a few iterations ahead */
			}
			rsi += -64*2*64+dist128;	/* Next source pointer */
		}
		rsi += -8*dist128;		/* Restore source pointer */

		/* Do FFT levels 2,3 */
		/* Values 0-1023 is real data, 1024-2047 is semi-real data, 2048-8191 is */
		/* complex data. */
		/* */
		/* On input the 64-byte cache lines hold these data values: */
		/*	0	1	1024	+1	2048	+1	3072	+1 */
		/*	4096	... */
		/*	2	... */
		/*	... */
		/*	1022	... */
		/*	5118	... */
		/* On output the 64-byte cache lines hold these data values: */
		/*	0	1	256	+1	512	+1	768	+1 */
		/*	4096	... */
		/*	2	... */
		/*	... */
		/*	254	... */
		/*	4350	... */
		/*	1024	... */
		/*	... */

		/* Do 128 eight_reals_fft_2 macros */
		/* Do 128 nop_two_two_complex_fft_2 macros */
		/* Do 256 four_complex_fft macros */
		/*	distance between fft data elements is 256 */
		/*	do 128 macros each processing 32 data values */

		rdi = (uintptr_t)g->u.xmm.sincos10;/* Load sin/cos pointer */
		xtouch_init(rcx = (uintptr_t)g->data_prefetch);/* Prefetch pointer */
		loops_init(128);
		for(unsigned int loopA = 2; loopA; loopA--) {
			for(unsigned int loopB = 64; loopB; loopB--) {
				x4cl_eight_reals_fft_2(rsi, 2*64, 2*dist128, 4*dist128, rdi);
			}
			rsi += -64*2*64+dist128;	/* Next source pointer */
		}

		/* Do 512 four_complex_fft macros */
		/*	distance between fft data elements is 256 */
		/*	do 128 macros each processing 32 data values */

		rsi += -2*dist128+64;	/* Next source pointer */
		rdi += XMM_SCD;		/* Next sine/cosine pointer */
		loops_init_prefetch(128, 128, 1, rcx);
		xtouch(xptr(rcx));			/* Preload first TLB */
		for(unsigned int loopA = 2; loopA; loopA--) {
			xtouch(xptr(rcx+4096-128));		/* Preload the TLBs */
			xtouch(xptr(rcx+2*4096-128));
			for(unsigned int loopB = 64; loopB; loopB--) {
				x4cl_four_complex_fft(rsi, 2*64, 2*dist128, 4*dist128, rdi);
			}
			rsi += -64*2*64+dist128;	/* Next source pointer */
			xtouch_init(rcx += 128);		/* Skip pad bytes every 8KB */
		}
		rsi += -2*dist128-64;	/* Restore source pointer */

		/* Do FFT levels 4,5 */
		/* Values 0-255 is real data, 256-511 is semi-real data, 512-8191 is */
		/* complex data. */
		/* */
		/* On input the 64-byte cache lines hold these data values: */
		/*	0	1	256	257	512	513	768	769 */
		/*	4096	... */
		/*	2	... */
		/*	... */
		/*	254	... */
		/*	4350	... */
		/*	1024	... */
		/*	... */
		/* On output the 64-byte cache lines hold these data values: */
		/*	0	1	64	65	128	129	192	193 */
		/*	4096	... */
		/*	2	... */
		/*	... */
		/*	62	... */
		/*	4158	... */
		/*	256	... */
		/*	... */

		/* Do 32 eight_reals_fft_2 macros */
		/* Do 32 nop_two_two_complex_fft_2 macros */
		/* Do 64 four_complex_fft macros */
		/*	distance between fft data elements is 64 */
		/*	do 32 macros each processing 32 data values */

		rdi = (uintptr_t)g->u.xmm.sincos9;/* Load sin/cos pointer */
		loops_init(256);
		for(unsigned int loopA = 32; loopA; loopA--) {
			x4cl_eight_reals_fft_2(rsi, 2*64, 64*64, dist128, rdi);
		}

		/* Do 896 four_complex_fft macros */
		/*	distance between fft data elements is 64 */
		/*	do 224 macros each processing 32 data values */

		rsi += -32*2*64+2*dist128;	/* Load source pointer */
		rdi += XMM_SCD;		/* Next sine/cosine pointer */
		loops_reset;			/* 3 then 4 iterations of 32 */
		for(unsigned int loopA = 2*4-1; loopA; ) {
			do{
				for(unsigned int loopC = 32; loopC; loopC--) {
					x4cl_four_complex_fft(rsi, 2*64, 64*64, dist128, rdi);
				}
				rsi += -32*2*64+2*dist128;	/* Next source pointer */
				rdi += 2*XMM_SCD;		/* Next sine/cosine pointer */
			}while(--loopA & 3);
			rsi += -8*dist128+64;	/* Next source pointer */
		}
		rsi += -2*64;		/* Restore source pointer */

		/* Do FFT levels 6,7 */
		/* Values 0-63 is real data, 64-127 is semi-real data, 128-8191 is */
		/* complex data. */
		/* */
		/* On input the 64-byte cache lines hold these data values: */
		/*	0	1	64	65	128	129	192	193 */
		/*	4096	... */
		/*	2	... */
		/*	... */
		/*	62	... */
		/*	4158	... */
		/*	256	... */
		/*	... */
		/* On output the 64-byte cache lines hold these data values: */
		/*	0	1	16	17	32	33	48	49 */
		/*	4096	... */
		/*	2	... */
		/*	... */
		/*	14	... */
		/*	4110	... */
		/*	64	... */
		/*	... */

		/* Do 8 eight_reals_fft_2 macros */
		/* Do 8 nop_two_two_complex_fft_2 macros */
		/* Do 16 four_complex_fft macros */
		/*	distance between fft data elements is 16 */
		/*	do 8 macros each processing 32 data values */

		rdi = (uintptr_t)g->u.xmm.sincos8;/* Load sin/cos pointer */
		loops_init(256);
		for(unsigned int loopA = 8; loopA; loopA--) {
			x4cl_eight_reals_fft_2(rsi, 2*64, 16*64, 32*64, rdi);
		}

		/* Do 992 four_complex_fft macros */
		/*	distance between fft data elements is 16 */
		/*	do 248 macros each processing 32 data values */

		rsi += -8*2*64+64*64;	/* Load source pointer */
		rdi += XMM_SCD;		/* Next sine/cosine pointer */
		loops_reset;			/* 2 iters of 8 iters of 1 or 2 of 8 */
		for(unsigned int loopA = 2*8*2-1; loopA; ) {
			do{
				do{
					for(unsigned int loopD = 8; loopD; loopD--) {
						x4cl_four_complex_fft(rsi, 2*64, 16*64, 32*64, rdi);
					}
					rsi += -8*2*64+64*64;	/* Next source pointer */
					rdi += 2*XMM_SCD;		/* Next sine/cosine pointer */
				}while(--loopA & 1);
				rsi += -2*64*64+dist128;	/* Next source pointer */
			}while(loopA & (8*2-1));
			rsi += -8*dist128+64;	/* Next source pointer */
		}
		rsi += -2*64;		/* Restore source pointer */

		/* Do FFT levels 8,9 */
		/* Values 0-15 is real data, 16-31 is semi-real data, 32-8191 is */
		/* complex data. */
		/* */
		/* On input the 64-byte cache lines hold these data values: */
		/*	0	1	16	17	32	33	48	49 */
		/*	4096	... */
		/*	2	... */
		/*	... */
		/*	14	... */
		/*	4110	... */
		/*	64	... */
		/*	... */
		/* On output the 64-byte cache lines hold these data values: */
		/*	0	1	4	5	8	9	12	13 */
		/*	4096	... */
		/*	2	... */
		/*	4098	... */
		/*	16	... */
		/*	... */

		/* Do 2 eight_reals_fft_2 macros */
		/* Do 2 nop_two_two_complex_fft_2 macros */
		/* Do 4 four_complex_fft macros */
		/*	distance between fft data elements is 4 */
		/*	do 2 macros each processing 32 data values */

		rdi = (uintptr_t)g->u.xmm.sincos7;/* Load sin/cos pointer */
		loops_init(256);
		for(unsigned int loopA = 2; loopA; loopA--) {
			x4cl_eight_reals_fft_2(rsi, 2*64, 4*64, 8*64, rdi);
		}

		/* Do 1016 four_complex_fft macros */
		/*	distance between fft data elements is 4 */
		/*	do 254 macros each processing 32 data values */

		rsi += -2*2*64+16*64;	/* Load source pointer */
		rdi += XMM_SCD;		/* Next sine/cosine pointer */
		loops_reset;			/* 2 iters of 8 iters of 7 or 8 of 2 */
		for(unsigned int loopA = 2*8*8-1; loopA; ) {
			do{
				do{
					for(unsigned int loopD = 2; loopD; loopD--) {
						x4cl_four_complex_fft(rsi, 2*64, 4*64, 8*64, rdi);
					}
					rsi += -2*2*64+16*64;	/* Next source pointer */
					rdi += 2*XMM_SCD;		/* Next sine/cosine pointer */
				}while(--loopA & 7);
				rsi += -8*16*64+dist128;	/* Next source pointer */
			}while(loopA & (8*8-1));
			rsi += -8*dist128+64;	/* Next source pointer */
		}
		rsi += -2*64;		/* Restore source pointer */

		/* Do FFT levels 10,11 */
		/* Values 0-3 is real data, 4-7 is semi-real data, 8-8191 is */
		/* complex data. */
		/* */
		/* On input the 64-byte cache lines hold these data values: */
		/*	0	1	4	5	8	9	12	13 */
		/*	4096	... */
		/*	2	... */
		/*	4098	... */
		/*	16	... */
		/*	... */
		/* On output the 64-byte cache lines hold these data values: */
		/*	0	1	2	3	4	5	6	7 */
		/*	4096	... */
		/*	8	... */
		/*	... */

		/* Do 1/2 eight_reals_fft_2 macros */
		/* Do 1/2 nop_two_two_complex_fft_2 macros */
		/* Do 1 four_complex_fft macros */
		/*	distance between fft data elements is 1 */

		/* Do 1022 more four_complex_fft macros */
		/*	distance between fft data elements is 1 */
		/*	do 511 macros each processing 16 data values */

		rdi = (uintptr_t)g->u.xmm.sincos6;/* Load sin/cos pointer */

		/* Execute the proper middle step */

		if likely(g->ffttype & 2){
			if likely(!(g->ffttype & 1)){
				s2cl_eight_reals_with_square_2(rsi, 4*64, 2*64);
				rdi += XMM_SCD;		/* Next sine/cosine pointer */
				loops_init_prefetch(512, 64, 1, rcx, 1, 32);
				for(unsigned int loopA = 2*4*2*32-1; loopA; ) {
					do{
						xtouch(xptr(rcx+4096-128));		/* Preload the TLBs */
						do{
							do{
								s2cl_four_complex_with_square(rsi, 4*64, 2*64);
								rdi += XMM_SCD;		/* Next sine/cosine pointer */
							}while(--loopA & (32-1));
							rsi += -32*4*64+dist128;	/* Next source pointer */
						}while(loopA & (2*32-1));
					}while(loopA & (4*2*32-1));
					rsi += -8*dist128+64;	/* Next source pointer */
				}
				rsi += -2*64;		/* Restore source pointer */
			}else{
				rbp = g->DIST_TO_MULSRCARG;
				s2cl_eight_reals_with_mult_2(rsi, 4*64, 2*64);
				rdi += XMM_SCD;		/* Next sine/cosine pointer */
				loops_init_prefetch(512, 64, 1, rcx, 1, 32);
				for(unsigned int loopA = 2*4*2*32-1; loopA; ) {
					do{
						xtouch(xptr(rcx+4096-128));		/* Preload the TLBs */
						do{
							do{
								s2cl_four_complex_with_mult(rsi, 4*64, 2*64);
								rdi += XMM_SCD;		/* Next sine/cosine pointer */
							}while(--loopA & (32-1));
							rsi += -32*4*64+dist128;	/* Next source pointer */
						}while(loopA & (2*32-1));
					}while(loopA & (4*2*32-1));
					rsi += -8*dist128+64;	/* Next source pointer */
				}
				rsi += -2*64;		/* Restore source pointer */
			}
		}else{
			s2cl_eight_reals_fft_2_final(rsi, 4*64, 2*64);
			rdi += XMM_SCD;		/* Next sine/cosine pointer */
			loops_init_prefetch(512, 64, 1, rcx, 1, 32);
			for(unsigned int loopA = 2*4*2*32-1; loopA; ) {
				do{
					xtouch(xptr(rcx+4096-128));		/* Preload the TLBs */
					do{
						do{
							s2cl_four_complex_fft_final(rsi, 4*64, 2*64);
							rdi += XMM_SCD;		/* Next sine/cosine pointer */
						}while(--loopA & (32-1));
						rsi += -32*4*64+dist128;	/* Next source pointer */
					}while(loopA & (2*32-1));
				}while(loopA & (4*2*32-1));
				rsi += -8*dist128+64;	/* Next source pointer */
			}
			goto	xpass2_11_real_done;
		}
	}else{
		rdi = (uintptr_t)g->u.xmm.sincos6;/* Load sin/cos pointer */
		rbp = g->DIST_TO_MULSRCARG;
		xtouch_init(rcx = (uintptr_t)g->data_prefetch);/* Load prefetch pointer */
		s2cl_eight_reals_with_mulf_2(rsi, 4*64, 2*64);
		rdi += XMM_SCD;		/* Next sine/cosine pointer */
		loops_init_prefetch(512, 64, 1, rcx, 1, 32);
		xtouch(xptr(rcx));			/* Preload the TLBs */
		for(unsigned int loopA = 2*4*2*32-1; loopA; ) {
			do{
				xtouch(xptr(rcx+4096-128));		/* Preload the TLBs */
				do{
					do{
						s2cl_four_complex_with_mulf(rsi, 4*64, 2*64);
						rdi += XMM_SCD;		/* Next sine/cosine pointer */
					}while(--loopA & (32-1));
					rsi += -32*4*64+dist128;	/* Next source pointer */
				}while(loopA & (2*32-1));
			}while(loopA & (4*2*32-1));
			rsi += -8*dist128+64;	/* Next source pointer */
		}
		rsi += -2*64;		/* Restore source pointer */

		/* Do inverse FFT levels 8,9 */
		/* On input the 64-byte cache lines hold these data values: */
		/*	0	1	2	3	4	5	6	7 */
		/*	4096	... */
		/*	8	... */
		/*	... */
		/* On output the 64-byte cache lines hold these data values: */
		/*	0	1	8	9	16	17	24	25 */
		/*	4096	... */
		/*	2	... */
		/*	... */
		/*	6	... */
		/*	4102	... */
		/*	32	... */
		/*	... */

		/* Do 2 eight_reals_unfft_2 macros */
		/* Do 2 nop_two_two_complex_unfft_2 macros */
		/*	distance between fft data elements is 4 */
		/*	do 1 macro processing 32 data values */

	}
	x4cl_eight_reals_unfft_2(rsi, 8*64, 2*64, 4*64);

	/* Do 1020 four_complex_unfft macros */
	/*	distance between fft data elements is 4 */
	/*	do 255 macros each processing 32 data values */

	rdi = (uintptr_t)g->u.xmm.sincos7;/* Load sin/cos pointer */
	loops_init(256, 1, 16);		/* 2 iters of 8 iters of 15 or 16 */
	for(unsigned int loopA = 2*8*16-1; loopA; ) {
		do{
			do{
				x4cl_four_complex_unfft(rsi, 8*64, 2*64, 4*64, rdi);
				rdi += XMM_SCD;		/* Next sine/cosine pointer */
			}while(--loopA & 15);
			rsi += -16*8*64+dist128;	/* Next source pointer */
		}while(loopA & (8*16-1));
		rsi += -8*dist128+64;	/* Next source pointer */
	}
	rsi += -2*64;		/* Restore source pointer */

	/* Do inverse FFT levels 6,7 */
	/* On input the 64-byte cache lines hold these data values: */
	/*	0	1	8	9	16	17	24	25 */
	/*	4096	... */
	/*	2	... */
	/*	... */
	/*	6	... */
	/*	4102	... */
	/*	32	... */
	/*	... */
	/* On output the 64-byte cache lines hold these data values: */
	/*	0	1	32	33	64	65	96	97 */
	/*	4096	... */
	/*	2	... */
	/*	... */
	/*	30	... */
	/*	4126	... */
	/*	128	... */
	/*	... */

	/* Do 8 eight_reals_unfft_2 macros */
	/* Do 8 nop_two_two_complex_unfft_2 macros */
	/*	distance between fft data elements is 16 */
	/*	do 4 macros each processing 32 data values */

	loops_init(256);
	for(unsigned int loopA = 4; loopA; loopA--) {
		x4cl_eight_reals_unfft_2(rsi, 2*64, 8*64, 16*64);
	}

	/* Do 1008 four_complex_unfft macros */
	/*	distance between fft data elements is 16 */
	/*	do 252 macros each processing 32 data values */

	rdi = (uintptr_t)g->u.xmm.sincos8;/* Load sin/cos pointer */
	rsi += -4*2*64+32*64;	/* Load source pointer */
	loops_reset;			/* 2 iters of 8 iters of 3 or 4 of 4 */
	for(unsigned int loopA = 2*8*4-1; loopA; ) {
		do{
			do{
				for(unsigned int loopD = 4; loopD; loopD--) {
					x4cl_four_complex_unfft(rsi, 2*64, 8*64, 16*64, rdi);
				}
				rsi += -4*2*64+32*64;	/* Next source pointer */
				rdi += XMM_SCD;		/* Next sine/cosine pointer */
			}while(--loopA & 3);
			rsi += -4*32*64+dist128;	/* Next source pointer */
		}while(loopA & (8*4-1));
		rsi += -8*dist128+64;	/* Next source pointer */
	}
	rsi += -2*64;		/* Restore source pointer */

	/* Do inverse FFT levels 4,5 */
	/* On input the 64-byte cache lines hold these data values: */
	/*	0	1	32	33	64	65	96	97 */
	/*	4096	... */
	/*	2	... */
	/*	... */
	/*	30	... */
	/*	4126	... */
	/*	128	... */
	/*	... */
	/* On output the 64-byte cache lines hold these data values: */
	/*	0	1	128	129	256	257	384	385 */
	/*	4096	... */
	/*	2	... */
	/*	... */
	/*	126	... */
	/*	4222	... */
	/*	512	... */
	/*	... */

	/* Do 32 eight_reals_unfft_2 macros */
	/* Do 32 nop_two_two_complex_unfft_2 macros */
	/*	distance between fft data elements is 64 */
	/*	do 16 macros each processing 32 data values */

	loops_init(256);
	for(unsigned int loopA = 16; loopA; loopA--) {
		x4cl_eight_reals_unfft_2(rsi, 2*64, 32*64, 64*64);
	}

	/* Do 960 four_complex_unfft macros */
	/*	distance between fft data elements is 64 */
	/*	do 240 macros each processing 32 data values */

	rdi = (uintptr_t)g->u.xmm.sincos9;/* Load sin/cos pointer */
	rsi += -16*2*64+dist128;	/* Load source pointer */
	loops_reset;			/* 2 iters of 7 or 8 iterations of 16 */
	for(unsigned int loopA = 2*8-1; loopA; ) {
		do{
			for(unsigned int loopC = 16; loopC; loopC--) {
				x4cl_four_complex_unfft(rsi, 2*64, 32*64, 64*64, rdi);
			}
			rsi += -16*2*64+dist128;	/* Next source pointer */
			rdi += XMM_SCD;		/* Next sine/cosine pointer */
		}while(--loopA & 7);
		rsi += -8*dist128+64;	/* Next source pointer */
	}
	rsi += -2*64;		/* Restore source pointer */

	/* Do inverse FFT levels 2,3 */
	/* On input the 64-byte cache lines hold these data values: */
	/*	0	1	128	129	256	257	384	385 */
	/*	4096	... */
	/*	2	... */
	/*	... */
	/*	126	... */
	/*	4222	... */
	/*	512	... */
	/*	... */
	/* On output the 64-byte cache lines hold these data values: */
	/*	0	1	512	513	1024	1025	1536	1537 */
	/*	4096	... */
	/*	2	... */
	/*	... */
	/*	510	... */
	/*	4606	... */
	/*	2048	... */
	/*	... */

	/* Do 128 eight_reals_unfft_2 macros */
	/* Do 128 nop_two_two_complex_unfft_2 macros */
	/*	distance between fft data elements is 256 */
	/*	do 64 macros each processing 32 data values */

	loops_init(256);
	for(unsigned int loopA = 64; loopA; loopA--) {
		x4cl_eight_reals_unfft_2(rsi, 2*64, dist128, 2*dist128);
	}

	/* Do 768 four_complex_unfft macros */
	/*	distance between fft data elements is 256 */
	/*	do 192 macros each processing 32 data values */

	rsi += -64*2*64+4*dist128;	/* Next source pointer */
	rdi = (uintptr_t)g->u.xmm.sincos10;/* Load sin/cos pointer */
	loops_reset;			/* 1 then 2 iterations of 64 */
	for(unsigned int loopA = 3; loopA; ) {
		do{
			for(unsigned int loopC = 64; loopC; loopC--) {
				x4cl_four_complex_unfft(rsi, 2*64, dist128, 2*dist128, rdi);
			}
			rsi += -64*2*64+4*dist128;	/* Next source pointer */
			rdi += XMM_SCD;		/* Next sine/cosine pointer */
		}while(--loopA & 1);
		rsi += -2*4*dist128+64;	/* Next source pointer */
	}
	rsi += -2*64;		/* Restore source pointer */

	/* Do inverse FFT level 1 */
	/* On input the 64-byte cache lines hold these data values: */
	/*	0	1	512	513	1024	1025	1536	1537 */
	/*	4096	... */
	/*	2	... */
	/*	... */
	/*	510	... */
	/*	4606	... */
	/*	2048	... */
	/*	... */
	/* On output the 64-byte cache lines hold these data values: */
	/*	0	+1K	2048	+1K	4096	+1K	6144	+1K */
	/*	1	... */
	/*	... */
	/*	1023	... */

	/* Do 512 eight_reals_last_unfft macros */
	/* Do 512 two_two_complex_unfft macros */
	/*	distance between fft data elements is 1024 */
	/*	do 256 macros each processing 16 data values */

	loops_init(256);
	for(unsigned int loopA = 4; loopA; loopA--) {
		for(unsigned int loopB = 64; loopB; loopB--) {
			s4cl_eight_reals_unfft_1(rsi, 2*64, 4*dist128, 64);
		}
		rsi += -64*2*64+dist128;	/* Next source pointer */
	}

	xpass2_11_real_done:;
	end_timer(2);


	if(g->pass2_get_next_block(g) != 0) return;
}
hg_pass2_11_levels_complex(g);
}

/* *************** xmm-pass2-11-levels-complex macro ****************** */
/* This macro takes 2 groups of 2048 complex values, premultiplies them and */
/* performs the final 11 levels of the FFT process, squares the */
/* results, and does the first 11 levels of the inverse FFT. */
/* NOTE: Rather than remembering 2048 multipliers we break them up */
/* into 512 groups of 4 columns.  By remembering these 516 values we can compute */
/* all 2048 multipliers using a single complex multiply.  These extra */
/* complex multiplies can probably be done faster than the computer can read the */
/* extra multipliers from memory (and are more likely stay in the L2 cache) */

void hg_pass2_11_levels_complex(struct gwasm_data *__restrict g) {
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
		/* a distance of 4096 data values. */
		/* */
		/* On input the 64-byte cache lines hold these data values: */
		/*	0	+1K	2048	+1K	4096	+1K	6144	+1K */
		/*	1	... */
		/*	... */
		/*	1023 */
		/* On output the 64-byte cache lines hold these data values: */
		/*	0	+4K	512	+4K	1024	+4K	1536	+4K */
		/*	1	... */
		/*	... */
		/*	511	... */
		/*	2048	... */
		/*	... */

		/* Do 1024 four_complex_gpm_fft macros */
		/*	distance between fft data elements is 512 */
		/*	do 512 macros each processing 16 data values */

		start_timer(4);
		rdi = (uintptr_t)g->premult_addr;/* The group multipliers */
		loops_init(512);			/* 128 iterations of 4 */
		for(unsigned int loopA = 4; loopA; loopA--) {
			for(unsigned int loopB = 16; loopB; loopB--) {
				for(unsigned int loopC = 8; loopC; loopC--) {
					s2cl_four_complex_gpm_fft(rsi, 64, 4*dist128);
					L2prefetch128(u8ptr(rsi+5*64));		/* Prefetch a few iterations ahead in */
					L2prefetch128(u8ptr(rsi+4*dist128+5*64));/* case they are not in the L2 cache */
				}
				rdi += XMM_PMD;		/* Next group multiplier */
				L2prefetch128(u8ptr(rdi+2*XMM_PMD));
			}
			rsi += -16*8*64+dist128;	/* Next source pointer */
		}
		rsi += -4*dist128;		/* Restore source pointer */
		end_timer(4);

		/* Do FFT levels 3,4 */
		/* */
		/* On input the 64-byte cache lines hold these data values: */
		/*	0	+4K	512	+4K	1024	+4K	1536	+4K */
		/*	1	... */
		/*	... */
		/*	511	... */
		/*	2048	... */
		/*	... */
		/* On output the 64-byte cache lines hold these data values: */
		/*	0	+4K	128	+4K	256	+4K	384	+4K */
		/*	1	... */
		/*	... */
		/*	127	... */
		/*	512	... */
		/*	... */

		/* Do 4 groups of four_complex_fft_cpm macros */
		/* To make things run faster we precompute the column multipliers times the */
		/* 4 different sine/cosine multipliers we run into in this section.  This */
		/* costs some memory but saves us some complex multiplies. */

		/* Do 1024 four_complex_cpm_fft macros */
		/*	distance between fft data elements is 128 */
		/*	do 128 macros each processing 64 data values */

		start_timer(5);
		xtouch_init(rcx = (uintptr_t)g->premult_prefetch);/* Prefetch pointer */
		loops_init_prefetch(256, 64, 1, rcx); /* 2 iters of 2 iters of 8 iters of 8 */
		xtouch(xptr(rcx));			/* Load first TLB */
		for(unsigned int loopA = 2; loopA; loopA--) {
			for(unsigned int loopB = 2; loopB; loopB--) {
				xtouch(xptr(rcx+4096-128));		/* Load prefetch TLBs */
				for(unsigned int loopC = 8; loopC; loopC--) {
					for(unsigned int loopD = 8; loopD; loopD--) {
						x4cl_four_complex_cpm_fft(rsi, 64, dist128, 2*dist128, 8*XMM_PMD);
						rdi += XMM_PMD;		/* Next column premultiplier */
					}
					rdi += -8*XMM_PMD;		/* Restore column premultiplier */
				}
			}
			rsi += -128*64+4*dist128;	/* Next source pointer */
			rdi += 16*XMM_PMD;		/* Next set of premultipliers */
			xtouch_init(rcx += -4096);		/* Back up so we only prefetch 12KB */
		}
		rsi += -2*4*dist128;	/* Restore source pointer */
		end_timer(5);

		/* Do FFT levels 5 */
		/* */
		/* On input the 64-byte cache lines hold these data values: */
		/*	0	+4K	128	+4K	256	+4K	384	+4K */
		/*	1	... */
		/*	... */
		/*	127	... */
		/*	512	... */
		/*	... */
		/* On output the 64-byte cache lines hold these data values: */
		/*	0	+4K	64	+4K	128	+4K	192	+4K */
		/*	1	... */
		/*	... */
		/*	63	... */
		/*	256	... */
		/*	... */
		/* */
		/* Do 1024 two_two_complex_fft macros */
		/*	distance between fft data elements is 64 */
		/*	do 512 macros each processing 16 data values */

		/* Do FFT levels 6,7 */
		/* */
		/* On input the 64-byte cache lines hold these data values: */
		/*	0	+4K	64	+4K	128	+4K	192	+4K */
		/*	1	... */
		/*	... */
		/*	63	... */
		/*	256	... */
		/*	... */
		/* On output the 64-byte cache lines hold these data values: */
		/*	0	+4K	16	+4K	32	+4K	48	+4K */
		/*	1	... */
		/*	... */
		/*	15	... */
		/*	64	... */
		/*	... */
		/* */
		/* Do 1024 four_complex_fft macros */
		/*	distance between fft data elements is 16 */
		/*	do 256 macros each processing 32 data values */

		/* Do FFT levels 8,9 */
		/* */
		/* On input the 64-byte cache lines hold these data values: */
		/*	0	+4K	16	+4K	32	+4K	48	+4K */
		/*	1	... */
		/*	... */
		/*	15	... */
		/*	64	... */
		/*	... */
		/* On output the 64-byte cache lines hold these data values: */
		/*	0	+4K	4	+4K	8	+4K	12	+4K */
		/*	1	... */
		/*	... */
		/*	3	... */
		/*	16	... */
		/*	... */
		/* */
		/* Do 1024 four_complex_fft macros */
		/*	distance between fft data elements is 4 */
		/*	do 256 macros each processing 32 data values */

		/* Do this 32 cache lines at a time to maximize accessing data */
		/* while it is in the L1 cache. */

		start_timer(6);
		xtouch_init(rcx = (uintptr_t)g->data_prefetch);/* Prefetch pointer */
		rdi = (uintptr_t)g->xsincos_complex;	/* Load sin/cos pointers */
		rbx = (uintptr_t)g->xsincos_complex;
		rbp = (uintptr_t)g->xsincos_complex;
		loops_init_prefetch(512, 64, 2, rcx); /* 2 iters of 4 of 4 of 32 cache lines */

		xtouch(xptr(rcx));			/* Load first TLB */
		xtouch(xptr(rcx+4096));		/* Load six prefetch TLBs */
		xtouch(xptr(rcx+2*4096-128));
		xtouch(xptr(rcx+3*4096-128));
		xtouch(xptr(rcx+4*4096-128));
		xtouch(xptr(rcx+5*4096-128));
		xtouch(xptr(rcx+6*4096-128));

		for(unsigned int loopA = 2*4; loopA; ) {
			do{
				for(unsigned int loopC = 4; loopC; loopC--) {
					for(unsigned int loopD = 16; loopD; loopD--) {
						x2cl_two_complex_fft(rsi, 4*64, 64*64, rdi);
					}
					rsi += -16*4*64;		/* Restore source pointer */
					xtouch_init(rcx += -4*128+8192+128);	/* Prefetch from next section */
					loops_undo(16);

					for(unsigned int loopD = 2; loopD; loopD--) {
						for(unsigned int loopE = 4; loopE; loopE--) {
							x4cl_four_complex_fft(rsi, 4*64, 16*64, 32*64, rbx);
						}
						rsi += -4*4*64+64*64;
						rbx += 2*XMM_SCD;		/* Next sine/cosine pointer */
					}
					rsi += -2*64*64;		/* Restore source pointer */
					rbx += -2*2*XMM_SCD;	/* Restore sine/cosine pointer */
					xtouch_init(rcx += -4*128+8192+128);	/* Prefetch from next section */
					loops_undo(2*4*2);

					for(unsigned int loopD = 8; loopD; loopD--) {
						x4cl_four_complex_fft(rsi, 16*64, 4*64, 8*64, rbp);
						rbp += 2*XMM_SCD;		/* Next sine/cosine pointer */
					}
					rbp += -8*2*XMM_SCD;	/* Restore sine/cosine pointer */
					xtouch_init(rcx += -2*(8192+128));	/* Back to prefetching from 1st section */

					rsi += -8*16*64+64;	/* Next source pointer */
				}
				rsi += -4*64+dist128;	/* Next source pointer */
				rdi += 2*XMM_SCD;		/* Next sine/cosine pointer */
				rbx += 2*2*XMM_SCD;	/* Next sine/cosine pointer */
				rbp += 8*2*XMM_SCD;	/* Next sine/cosine pointer */
			}while(--loopA & 3);
			xtouch_init(rcx += -16*4*128);		/* Re-prefetch 3 sections */
		}
		rsi += -8*dist128;		/* Restore source pointer */
		xtouch_init(rcx += 3*(8192+128));	/* Prefetch next 3 sections */
		end_timer(6);

		/* Do FFT levels 10,11 as well as inverse FFT levels 10,11 */
		/* */
		/* On input the 64-byte cache lines hold these data values: */
		/*	0	+4K	4	+4K	8	+4K	12	+4K */
		/*	1	... */
		/*	... */
		/*	3	... */
		/*	16	... */
		/*	... */
		/* On output the 64-byte cache lines hold these data values: */
		/*	0	+4K	2	+4K	4	+4K	6	+4K */
		/*	1	... */
		/*	8	... */
		/*	... */

		/* Do 1024 four_complex_fft macros */
		/*	distance between fft data elements is 1 */
		/*	do 256 macros each processing 32 data values */

		/* Execute the right middle step */

		if likely(g->ffttype & 2){
			if likely(!(g->ffttype & 1)){
				start_timer(9);
				loops_init_prefetch(256, 64, 1, rcx);
				xtouch(xptr(rcx+4096-128));		/* Prefetch TLB */
				xtouch(xptr(rcx+2*4096-128));
				for(unsigned int loopA = 2*2; loopA; ) {
					do{
						rdi = (uintptr_t)g->xsincos_complex;	/* Load sin/cos pointer */
						for(unsigned int loopC = 2; loopC; loopC--) {
							for(unsigned int loopD = 32; loopD; loopD--) {
								x4cl_four_complex_with_square(rsi, 4*64, 64, 2*64);
								rdi += 2*XMM_SCD;		/* Next sine/cosine pointer */
							}
							rsi += -16*2*4*64+dist128;	/* Next source pointer */
						}
					}while(--loopA & 1);
					xtouch_init(rcx += -64*128);		/* Reprefetch section */
				}
				rsi += -8*dist128;		/* Restore source pointer */
				xtouch_init(rcx += 8192+128);		/* Don't prefetch pad bytes */
				end_timer(9);
			}else{
				start_timer(9);
				rbp = g->DIST_TO_MULSRCARG;
				loops_init_prefetch(256, 64, 1, rcx);
				xtouch(xptr(rcx+4096-128));		/* Prefetch TLB */
				xtouch(xptr(rcx+2*4096-128));
				for(unsigned int loopA = 2*2; loopA; ) {
					do{
						rdi = (uintptr_t)g->xsincos_complex;	/* Load sin/cos pointer */
						for(unsigned int loopC = 2; loopC; loopC--) {
							for(unsigned int loopD = 32; loopD; loopD--) {
								x4cl_four_complex_with_mult(rsi, 4*64, 64, 2*64);
								rdi += 2*XMM_SCD;		/* Next sine/cosine pointer */
							}
							rsi += -16*2*4*64+dist128;	/* Next source pointer */
						}
					}while(--loopA & 1);
					xtouch_init(rcx += -64*128);		/* Back to prefetching from 1st section */
				}
				rsi += -8*dist128;		/* Restore source pointer */
				xtouch_init(rcx += 8192+128);		/* Don't prefetch pad bytes */
				end_timer(9);
			}
		}else{
			start_timer(9);
			loops_init_prefetch(256, 64, 1, rcx);
			for(unsigned int loopA = 2*2; loopA; ) {
				xtouch(xptr(rcx+4096-128));		/* Prefetch TLB */
				xtouch(xptr(rcx+2*4096-128));
				do{
					rdi = (uintptr_t)g->xsincos_complex;	/* Load sin/cos pointer */
					for(unsigned int loopC = 2; loopC; loopC--) {
						for(unsigned int loopD = 32; loopD; loopD--) {
							x4cl_four_complex_fft_final(rsi, 4*64, 64, 2*64);
							rdi += 2*XMM_SCD;		/* Next sine/cosine pointer */
						}
						rsi += -16*2*4*64+dist128;	/* Next source pointer */
					}
				}while(--loopA & 1);
				xtouch_init(rcx += 128);		/* Don't prefetch pad bytes */
			}
			end_timer(9);
			goto	xpass2_11_levels_complex_done;
		}
	}else{
		start_timer(9);
		xtouch_init(rcx = (uintptr_t)g->data_prefetch);/* Prefetch pointer */
		rbp = g->DIST_TO_MULSRCARG;
		loops_init_prefetch(256, 64, 1, rcx);
		xtouch(xptr(rcx));			/* Load first prefetch TLB */
		for(unsigned int loopA = 2*2; loopA; ) {
			xtouch(xptr(rcx+4096-128));		/* Prefetch TLB */
			xtouch(xptr(rcx+2*4096-128));
			do{
				rdi = (uintptr_t)g->xsincos_complex;	/* Load sin/cos pointer */
				for(unsigned int loopC = 2; loopC; loopC--) {
					for(unsigned int loopD = 32; loopD; loopD--) {
						x4cl_four_complex_with_mulf(rsi, 4*64, 64, 2*64);
						rdi += 2*XMM_SCD;		/* Next sine/cosine pointer */
					}
					rsi += -16*2*4*64+dist128;	/* Next source pointer */
				}
			}while(--loopA & 1);
			xtouch_init(rcx += 128);		/* Don't prefetch pad bytes */
		}
		rsi += -8*dist128;		/* Restore source pointer */
		end_timer(9);

		/* Do inverse FFT levels 8,9 */
		/* On input the 64-byte cache lines hold these data values: */
		/*	0	+4K	2	+4K	4	+4K	6	+4K */
		/*	1	... */
		/*	8	... */
		/*	... */
		/* On output the 64-byte cache lines hold these data values: */
		/*	0	+4K	8	+4K	16	+4K	24	+4K */
		/*	1	... */
		/*	... */
		/*	7	... */
		/*	32	... */
		/*	... */
		/* */
		/* Do 1024 four_complex_unfft macros */
		/*	distance between fft data elements is 4 */
		/*	do 256 macros each processing 32 data values */

		/* Do inverse FFT levels 6,7 */
		/* On input the 64-byte cache lines hold these data values: */
		/*	0	+4K	8	+4K	16	+4K	24	+4K */
		/*	1	... */
		/*	... */
		/*	7	... */
		/*	32	... */
		/*	... */
		/* On output the 64-byte cache lines hold these data values: */
		/*	0	+4K	32	+4K	64	+4K	96	+4K */
		/*	1	... */
		/*	... */
		/*	31	... */
		/*	128	... */
		/*	... */
		/* */
		/* Do 1024 four_complex_unfft macros */
		/*	distance between fft data elements is 16 */
		/*	do 256 macros each processing 32 data values */

		/* Do inverse FFT levels 5 */
		/* On input the 64-byte cache lines hold these data values: */
		/*	0	+4K	32	+4K	64	+4K	96	+4K */
		/*	1	... */
		/*	... */
		/*	31	... */
		/*	128	... */
		/*	... */
		/* On output the 64-byte cache lines hold these data values: */
		/*	0	+4K	64	+4K	128	+4K	192	+4K */
		/*	1	... */
		/*	... */
		/*	63	... */
		/*	256	... */
		/*	... */
		/* */
		/* Do 512 two_two_complex_unfft macros */
		/*	distance between fft data elements is 64 */
		/*	do 512 macros each processing 16 data values */

		/* Do this 32 cache lines at a time to maximize accessing data */
		/* while it is in the L1 cache. */

	}
	start_timer(10);
	rdi = (uintptr_t)g->xsincos_complex;	/* Load sin/cos pointers */
	rbx = (uintptr_t)g->xsincos_complex;
	rbp = (uintptr_t)g->xsincos_complex;
	loops_init_prefetch(512, 64, 2, rcx); /* 32 iterations of 16 */

	xtouch(xptr(rcx+4096-128));		/* Load six prefetch TLBs */
	xtouch(xptr(rcx+2*4096-128));
	xtouch(xptr(rcx+3*4096-128));
	xtouch(xptr(rcx+4*4096-128));
	xtouch(xptr(rcx+5*4096-128));
	xtouch(xptr(rcx+6*4096-128));
	for(unsigned int loopA = 2*4; loopA; ) {
		do{
			for(unsigned int loopC = 2; loopC; loopC--) {
				for(unsigned int loopD = 2; loopD; loopD--) {
					for(unsigned int loopE = 8; loopE; loopE--) {
						x4cl_four_complex_unfft(rsi, 8*64, 2*64, 4*64, rbp);
						rbp += XMM_SCD;		/* Next sine/cosine pointer */
					}
					rsi += -8*8*64;		/* Restore source pointer */
					rbp += -8*XMM_SCD;		/* Next sine/cosine pointer */
					xtouch_init(rcx += -4*128+8192+128);	/* Prefetch from next section */
					loops_undo(8*2);

					for(unsigned int loopE = 2; loopE; loopE--) {
						for(unsigned int loopF = 4; loopF; loopF--) {
							x4cl_four_complex_unfft(rsi, 2*64, 8*64, 16*64, rbx);
						}
						rsi += -4*2*64+32*64;	/* Next source pointer */
						rbx += XMM_SCD;		/* Next sine/cosine pointer */
					}
					rsi += -2*32*64;		/* Restore source pointer */
					rbx += -2*XMM_SCD;		/* Restore sine/cosine pointer */
					xtouch_init(rcx += -4*128+8192+128);	/* Prefetch from next section */
					loops_undo(2*4*2);

					for(unsigned int loopE = 16; loopE; loopE--) {
						x2cl_two_complex_unfft(rsi, 2*64, 32*64);
					}
					xtouch_init(rcx += -2*(8192+128));	/* Back to prefetching from 1st section */

					rsi += -16*2*64+64;	/* Next source pointer */
				}
				rsi += -2*64+64*64;	/* Next source pointer */
				rdi += XMM_SCD;		/* Next sine/cosine pointer */
				rbx += 2*XMM_SCD;		/* Next sine/cosine pointer */
				rbp += 8*XMM_SCD;		/* Next sine/cosine pointer */
			}
			rsi += -2*64*64+dist128;	/* Next source pointer */
		}while(--loopA & 3);
		xtouch_init(rcx += -16*4*128);		/* Prefetch next 3 sections */
	}
	rsi += -8*dist128;		/* Restore source pointer */
	xtouch_init(rcx += 3*(8192+128));	/* Prefetch next section */
	end_timer(10);

	/* Do inverse FFT levels 3,4 */
	/* On input the 64-byte cache lines hold these data values: */
	/*	0	+4K	64	+4K	128	+4K	192	+4K */
	/*	1	... */
	/*	... */
	/*	63	... */
	/*	256	... */
	/*	... */
	/* On output the 64-byte cache lines hold these data values: */
	/*	0	+4K	256	+4K	512	+4K	768	+4K */
	/*	1	... */
	/*	... */
	/*	255	... */
	/*	1024	... */
	/*	... */

	/* Do 4 groups of 256 four_complex_fft macros (FFT levels 3,4) */
	/* In this section we also apply the 4 column multipliers (actually three */
	/* since the first column multiplier is always one).  To make things */
	/* run faster we precompute the column multipliers times the 4 different */
	/* sine/cosine multipliers we run into in this section.  This costs some */
	/* memory but saves us some complex multiplies. */

	/* Do 1024 four_complex_cpm_unfft macros */
	/*	distance between fft data elements is 128 */
	/*	do 64 macros each processing 128 data values */

	start_timer(13);
	rdi = (uintptr_t)g->premult_addr;/* Load column multipliers pointer */
	rdi += 64*XMM_PMD;		/* Load column multipliers pointer */
	loops_init_prefetch(256, 64, 1, rcx); /* 2 iters of 2 iters of 8 iters of 8 */
	xtouch(xptr(rcx+4096-128));		/* Load two prefetch TLBs */
	xtouch(xptr(rcx+2*4096-128));
	for(unsigned int loopA = 2*2; loopA; ) {
		do{
			for(unsigned int loopC = 8; loopC; loopC--) {
				for(unsigned int loopD = 8; loopD; loopD--) {
					x4cl_four_complex_cpm_unfft(rsi, 64, 64*64, dist128);
					rdi += XMM_PMD;		/* Next column premultiplier */
				}
				rdi += -8*XMM_PMD;		/* Reset premultiplier pointer */
			}
			rsi += -64*64+2*dist128;	/* Next source pointer */
			rdi += 8*XMM_PMD;		/* Next premultiplier pointer */
		}while(--loopA & 1);
		xtouch_init(rcx += -64*128);		/* Re-prefetch section */
	}
	rsi += -4*2*dist128;	/* Restore source pointer */
	end_timer(13);

	/* Do inverse FFT levels 1,2 */
	/* On input the 64-byte cache lines hold these data values: */
	/*	0	+4K	256	+4K	512	+4K	768	+4K */
	/*	1	... */
	/*	... */
	/*	255	... */
	/*	1024	... */
	/*	... */
	/* On output the 64-byte cache lines hold these data values: */
	/*	0	+1K	2048	+1K	4096	+1K	6144	+1K */
	/*	1	... */
	/*	... */
	/*	1023 */

	/* Do 256 four_complex_gpm_unfft macros */
	/*	distance between fft data elements is 512 */
	/*	do 256 macros each processing 32 data values */

	start_timer(14);
	xtouch_init(rcx = (uintptr_t)g->premult_prefetch);/* Prefetch pointer */
	rdi = (uintptr_t)g->premult_addr;/* Load group multiplier pointer */
	loops_init_prefetch(256, 64, 1, rcx); /* 2 iters of 16 iterations of 8 */
	xtouch(xptr(rcx));			/* Load prefetch TLB */
	xtouch(xptr(rcx+4096-128));		/* Load prefetch TLBs */
	for(unsigned int loopA = 2; loopA; loopA--) {
		xtouch(xptr(rcx+2*4096-128));
		for(unsigned int loopB = 16; loopB; loopB--) {
			for(unsigned int loopC = 8; loopC; loopC--) {
				s4cl_four_complex_gpm_unfft(rsi, 64, 2*dist128, 4*dist128, 32*XMM_PMD);
			}
			rdi += XMM_PMD;		/* Next set of multipliers */
		}
		rsi += -16*4*2*64+dist128;	/* Next source pointer */
		xtouch_init(rcx += -4096);		/* Back up so we only prefetch 12KB */
	}
	end_timer(14);

	xpass2_11_levels_complex_done:;

}while(g->pass2_get_next_block(g) == 0);
}
