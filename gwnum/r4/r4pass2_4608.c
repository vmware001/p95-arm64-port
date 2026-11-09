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

#include "r3.h"
#include "r4.h"
#include "r8.h"

/********************************************************************************************** */
/*			Macros derived from pass2-13-levels */
/********************************************************************************************** */


/* *************** radix4-pass2-4608-levels-real macro ****************** */
/* This macro takes 9216 real values and performs the final 13.17 levels of */
/* the FFT.  Independently, this macro also takes 4608 complex values and */
/* performs the final 12.17 levels of the FFT.  The results are squared, and */
/* the first 12.17/13.17 levels of the inverse FFT are done. */


void r4_pass2_4608_levels_complex(struct gwasm_data *__restrict g);
void xpass2_r4_4608_levels(struct gwasm_data *__restrict g) {
	g->thread_work_routine = r4_pass2_4608_levels_complex;
	pass2_wake_up_threads(g);
	if(g->NEGACYCLIC_FFT == 0) {
//void r4_pass2_4608_levels_real(struct gwasm_data *__restrict g)
	uintptr_t rdx,rbx,rbp,rsi,rdi,daddr;
	vec2f64 xmm0,xmm1,xmm2,xmm3,xmm4,xmm5,xmm6,xmm7;

	/* We switch to the same format used in one-pass FFTs (that is, the */
	/* high word of xmm registers being 1 greater than the low word): */

	/* Type 4 FFTs skip the forward FFT process */

	start_timer(2);
	daddr = (uintptr_t)g->data_addr;/* Load source address */
	rsi = daddr;
	rbx = g->DIST_TO_FFTSRCARG;
	if likely(!(g->ffttype & 4)) { //xpass2_4608_levels_real_4;

		/* Do FFT level 1,2 (actually 1.585 levels) */
		/* */
		/* On input the 64-byte cache lines hold these data values if pass 2 swizzles: */
		/*	0	+2304	4608	+2304	9216	+2304	13824	+2304 */
		/*	1	... */
		/*	... */
		/*	2303 */
		/* On input the 64-byte cache lines hold these data values if pass 1 swizzles: */
		/*	0	+9216	2304	+9216	4608	+9216	6912	+9216 */
		/*	1	... */
		/*	... */
		/*	2303	... */
		/* On output the 64-byte cache lines hold these data values: */
		/*	0	+9216	768	+9216	1536	+9216	2304	+9216 */
		/*	1	... */
		/*	... */
		/*	767	... */
		/*	3072	... */
		/*	... */

		/* Do 1536 six_reals_fft operations */
		/* Do 1536 three_complex_djbfft operations */
		/*	distance between fft data elements is 1536 */
		/*	do 768 iterations */

		rdi = (uintptr_t)g->xsincos_complex;	/* Ptr to the sin/cos data */
		rdx = (uintptr_t)g->sincos3;/* Ptr to real sin/cos data */
		loops_init(768);
		for(unsigned int loopA = 6; loopA; loopA--) {
			for(unsigned int loopB = 128; loopB; loopB--) {
				r3_fh3cl_six_reals_three_complex_djbfft(rsi, 64, 6*dist128, rdi, 768*XMM_SCD1, rdx, 768*XMM_HSCD1);
				rdi += XMM_SCD1;		/* Next sin/cos data */
				rdx += XMM_HSCD1;		/* Next sin/cos data */
			}
			rsi += -128*64+dist128;	/* Next source pointer */
		}
		rsi = daddr; //rsi += -6*dist128;		/* Restore source pointer */
		rdi += -768*XMM_SCD1;	/* Restore sin/cos ptr */
		rdx += -768*XMM_HSCD1;	/* Restore sin/cos ptr */

		/* Do FFT levels 3,4 (actually 1.585 levels) */
		/* */
		/* On input the 64-byte cache lines hold these data values: */
		/*	0	+9216	768	+9216	1536	+9216	2304	+9216 */
		/*	1	... */
		/*	... */
		/*	767	... */
		/*	3072	... */
		/*	... */
		/* On output the 64-byte cache lines hold these data values: */
		/*	0	+9216	256	+9216	512	+9216	768	+9216 */
		/*	1	... */
		/*	... */
		/*	255	... */
		/*	1024	... */
		/*	... */

		/* Do 512 six_reals_fft operations */
		/* Do 512 three_complex_djbfft operations */
		/*	distance between fft data elements is 512 */
		/*	do 256 macros each processing 24 data values */

		pfing(rcx = (uintptr_t)g->data_prefetch); /* Prefetch pointer */
		xtouch(xptr(rcx));			/* Preload first TLB */
		loops_init_prefetch(768, 64, 4, rcx);
		xtouch(xptr(rcx+4096-64));		/* Preload the TLBs */
		for(unsigned int loopA = 2; loopA; loopA--) {
			for(unsigned int loopB = 128; loopB; loopB--) {
				r3_h3cl_six_reals_three_complex_djbfft(rsi, 64, 2*dist128, rdi, 256*3*XMM_SCD1, rdx, 256*3*XMM_HSCD1);
				rdi += 3*XMM_SCD1;		/* Next sin/cos data */
				rdx += 3*XMM_HSCD1;	/* Next sin/cos data */
			}
			rsi += -128*64+dist128;	/* Next source pointer */
		}
		rsi += -2*dist128+6*dist128; /* Next source pointer */
		rdi += -256*3*XMM_SCD1;	/* Restore sin/cos ptr */
		rdx += 256*3*XMM_HSCD1;	/* Next sin/cos ptr */

		/* Do 2048 three_complex_fft operations */
		/*	distance between fft data elements is 512 */
		/*	do 512 macros each processing 24 data values */

		loops_reset; /* 2 more of 2 x 128 */
		xtouch(xptr(rcx+4096-64));		/* Preload the TLBs */
		xtouch(xptr(rcx+2*4096-64));
		r3_nf3cl_three_complex_djbfft_preload;
		for(unsigned int loopA = 2; loopA; loopA--) {
			for(unsigned int loopB = 2; loopB; loopB--) {
				for(unsigned int loopC = 128; loopC; loopC--) {
					r3_nf3cl_three_complex_djbfft(rsi, 64, 2*dist128, rdi, 256*3*XMM_SCD1);
					rdi += 3*XMM_SCD1;		/* Next sin/cos data */
				}
				rsi += -128*64+dist128;	/* Next source pointer */
			}
			rsi += -2*dist128+6*dist128; /* Next source pointer */
			rdi += -256*3*XMM_SCD1;	/* Restore sin/cos ptr */
		}
		rsi = daddr; //rsi += -3*6*dist128;	/* Restore source pointer */

		/* Do FFT levels 5,6 */
		/* */
		/* On input the 64-byte cache lines hold these data values: */
		/*	0	+9216	256	+9216	512	+9216	768	+9216 */
		/*	1	... */
		/*	... */
		/*	255	... */
		/*	1024	... */
		/*	... */
		/* On output the 64-byte cache lines hold these data values: */
		/*	0	+9216	128	+9216	256	+9216	384	+9216 */
		/*	1	... */
		/*	... */
		/*	127	... */
		/*	512	... */
		/*	... */

		/* Do 128 eight_reals_fft operations */
		/* Do 384 four_complex_djbfft operations */
		/*	distance between fft data elements is 128 */
		/*	do 128 macros each processing 16 data values */

		rdi += 1536*XMM_SCD1;	/* Load sin/cos pointer */
		loops_init_prefetch(1152, 64, 4, rcx);
		xtouch(xptr(rcx+4096-64));		/* Preload the TLBs */
		xtouch(xptr(rcx+2*4096-64));
		xtouch(xptr(rcx+3*4096-64));
		xtouch(xptr(rcx+4*4096-64));
		xtouch(xptr(rcx+5*4096-64));
		for(unsigned int loopA = 128; loopA; loopA--) {
			r4_h2cl_eight_reals_four_complex_djbfft(rsi, 64, dist128, rdi, rdx);
			rdi += XMM_SCD2;		/* Next sin/cos data */
			rdx += XMM_HSCD2;		/* Next sin/cos data */
		}
		rsi += -128*64+2*dist128;	/* Next source pointer */
		rdi += -128*XMM_SCD2;	/* Restore sin/cos ptr */
		rdx += -128*XMM_HSCD2;	/* Restore sin/cos ptr */

		/* Do 2048 four_complex_djbfft operations */
		/*	distance between fft data elements is 128 */
		/*	do 1024 macros each processing 16 data values */

		loops_reset;			/* 8 more iterations of 128 */
		r4_nf2cl_four_complex_djbfft_preload;
		for(unsigned int loopA = 8; loopA; loopA--) {
			for(unsigned int loopB = 128; loopB; loopB--) {
				r4_nf2cl_four_complex_djbfft(rsi, 64, dist128, rdi);
				rdi += XMM_SCD2;		/* Next sin/cos data */
			}
			rsi += -128*64+2*dist128;	/* Next source pointer */
			rdi += -128*XMM_SCD2;	/* Restore sin/cos ptr */
		}
		rsi = daddr; //rsi += -9*2*dist128;	/* Next source pointer */

		/* Do FFT levels 7,8 */
		/* */
		/* On input the 64-byte cache lines hold these data values: */
		/*	0	+9216	128	+9216	256	+9216	384	+9216 */
		/*	1	... */
		/*	... */
		/*	127	... */
		/*	512	... */
		/*	... */
		/* On output the 64-byte cache lines hold these data values: */
		/*	0	+9216	32	+9216	64	+9216	96	+9216 */
		/*	1	... */
		/*	... */
		/*	31	... */
		/*	128	... */
		/*	... */

		/* Do 32 eight_reals_fft operations */
		/* Do 96 four_complex_djbfft operations */
		/*	distance between fft data elements is 32 */
		/*	do 32 macros each processing 32 data values */

		loops_init_prefetch(576, 64, 2, rcx);
		xtouch(xptr(rcx+4096-64));		/* Preload the TLBs */
		xtouch(xptr(rcx+2*4096-64));
		xtouch(xptr(rcx+3*4096-64));
		xtouch(xptr(rcx+4*4096-64));
		xtouch(xptr(rcx+5*4096-64));
		r4_h4cl_eight_reals_four_complex_djbfft_preload;
		for(unsigned int loopA = 32; loopA; loopA--) {
			r4_h4cl_eight_reals_four_complex_djbfft(rsi, 64, 32*64, 64*64, rdi, rdx);
			rdi += 4*XMM_SCD2;		/* Next sin/cos data */
			rdx += 4*XMM_HSCD2;	/* Next sin/cos data */
		}
		rsi += -32*64+dist128;	/* Next source pointer */
		rdi += -32*4*XMM_SCD2;	/* Restore sin/cos ptr */
		rdx += -32*4*XMM_HSCD2;	/* Restore sin/cos ptr */

		/* Do 2176 four_complex_djbfft operations */
		/*	distance between fft data elements is 32 */
		/*	do 544 macros each processing 32 data values */

		loops_reset;			/* 17 more iters of 32 */
		r4_x4cl_four_complex_djbfft_preload;
		for(unsigned int loopA = 17; loopA; loopA--) {
			for(unsigned int loopB = 32; loopB; loopB--) {
				r4_x4cl_four_complex_djbfft(rsi, 64, 32*64, 64*64, rdi);
				rdi += 4*XMM_SCD2;		/* Next sin/cos data */
			}
			rsi += -32*64+dist128;	/* Next source pointer */
			rdi += -32*4*XMM_SCD2;	/* Restore sin/cos ptr */
		}
		rsi = daddr; //rsi += -18*dist128;	/* Next source pointer */

		/* Do FFT levels 9,10 */
		/* */
		/* On input the 64-byte cache lines hold these data values: */
		/*	0	+9216	32	+9216	64	+9216	96	+9216 */
		/*	1	... */
		/*	... */
		/*	31	... */
		/*	128	... */
		/*	... */
		/* On output the 64-byte cache lines hold these data values: */
		/*	0	+9216	8	+9216	16	+9216	24	+9216 */
		/*	1	... */
		/*	... */
		/*	7	... */
		/*	32	... */
		/*	... */

		/* Do 8 eight_reals_fft operations */
		/* Do 24 four_complex_djbfft operations */
		/*	distance between fft data elements is 8 */
		/*	do 8 macros each processing 32 data values */

		loops_init_prefetch(576, 64, 2, rcx);
		xtouch(xptr(rcx+4096-64));		/* Preload the TLBs */
		xtouch(xptr(rcx+2*4096-64));
		xtouch(xptr(rcx+3*4096-64));
		xtouch(xptr(rcx+4*4096-64));
		xtouch(xptr(rcx+5*4096-64));
		r4_h4cl_eight_reals_four_complex_djbfft_preload;
		for(unsigned int loopA = 8; loopA; loopA--) {
			r4_h4cl_eight_reals_four_complex_djbfft(rsi, 64, 8*64, 16*64, rdi, rdx);
			rdi += 16*XMM_SCD2;	/* Next sin/cos data */
			rdx += 16*XMM_HSCD2;	/* Next sin/cos data */
		}
		rsi += -8*64+32*64;	/* Next source pointer */
		rdi += -8*16*XMM_SCD2;	/* Restore sin/cos ptr */
		rdx += -8*16*XMM_HSCD2;	/* Restore sin/cos ptr */

		/* Do 2272 four_complex_djbfft operations */
		/*	distance between fft data elements is 8 */
		/*	do 568 macros each processing 32 data values */

		loops_reset;			/* 18 iters of 3 or 4 of 8 */
		r4_x4cl_four_complex_djbfft_preload;
		for(unsigned int loopA = 18*4-1; loopA; ) {
			do{
				for(unsigned int loopC = 8; loopC; loopC--) {
					r4_x4cl_four_complex_djbfft(rsi, 64, 8*64, 16*64, rdi);
					rdi += 16*XMM_SCD2;	/* Next sin/cos data */
				}
				rsi += -8*64+32*64;	/* Next source pointer */
				rdi += -8*16*XMM_SCD2;	/* Restore sin/cos ptr */
			}while(--loopA & 3);
			rsi += -4*32*64+dist128;	/* Next source pointer */
		}
		rsi = daddr; //rsi += -18*dist128;	/* Next source pointer */

		/* Do FFT levels 11,12,13 */
		/* Values 0-7 is real data, 8-9215 is complex data. */
		/* */
		/* On input the 64-byte cache lines hold these data values: */
		/*	0	+9216	8	+9216	16	+9216	24	+9216 */
		/*	1	... */
		/*	... */
		/*	7	... */
		/*	32	... */
		/*	... */
		/* On output the 64-byte cache lines hold these data values: */
		/*	0	+9216	4	+9216	8	+9216	12	+9216 */
		/*	1	... */
		/*	... */
		/*	3	... */
		/*	16	... */
		/*	... */

		/* Do 1 sixteen_reals_fft operations */
		/* Do 3 eight_complex_fft operations */
		/*	distance between fft data elements is 1 */

		/* Do 1148 eight_complex_fft operations */
		/*	distance between fft data elements is 1 */
		/*	do 287 macros each processing 64 data values */

		/* Execute the proper middle step */

		if likely(g->ffttype & 2) {
			if likely(!(g->ffttype & 1)) { // xpass2_4608_levels_real_2:;
				r8_h8cl_sixteen_reals_eight_complex_with_square(rsi, 8*64, 64, 2*64, 4*64);
				loops_init_prefetch(288, 64, 8, rcx, 1, 16);
				xtouch(xptr(rcx+4096-64));		/* Preload the TLBs */
				r8_x8cl_eight_complex_with_square_preload;
				for(unsigned int loopA = 18*16-1; loopA; ) {
					do{
						r8_x8cl_eight_complex_with_square(rsi, 8*64, 64, 2*64, 4*64);
					}while(--loopA & 15);
					rsi += -16*8*64+dist128;	/* Next source pointer */
				}
				//rsi += -18*dist128;	/* Next source pointer */
			}else{ // xpass2_4608_levels_real_3:;
				rbp = g->DIST_TO_MULSRCARG;
				r8_h8cl_sixteen_reals_eight_complex_with_mult(rsi, 8*64, 64, 2*64, 4*64);
				loops_init_prefetch(288, 64, 8, rcx, 1, 16);
				xtouch(xptr(rcx+4096-64));		/* Preload the TLBs */
				r8_x8cl_eight_complex_with_mult_preload;
				for(unsigned int loopA = 18*16-1; loopA; ) {
					do{
						r8_x8cl_eight_complex_with_mult(rsi, 8*64, 64, 2*64, 4*64);
					}while(--loopA & 15);
					rsi += -16*8*64+dist128;	/* Next source pointer */
				}
				//rsi += -18*dist128;	/* Next source pointer */
			}
		}else{ // xpass2_4608_levels_real_1:;
			r8_h8cl_sixteen_reals_eight_complex_fft_final(rsi, 8*64, 64, 2*64, 4*64);
			loops_init_prefetch(288, 128, 1, rcx, 1, 16);
			r8_x8cl_eight_complex_fft_final_preload;
			for(unsigned int loopA = 18*16-1; loopA; ) {
				xtouch(xptr(rcx+4096-64));		/* Preload the TLBs */
				do{
					r8_x8cl_eight_complex_fft_final(rsi, 8*64, 64, 2*64, 4*64);
				}while(--loopA & 15);
				rsi += -16*8*64+dist128;	/* Next source pointer */
			}
			goto	xpass2_4608_real_done;
		}
	}else{ // xpass2_4608_levels_real_4:;
		rbp = g->DIST_TO_MULSRCARG;
		pfing(rcx = (uintptr_t)g->data_prefetch); /* Load prefetch pointer */
		xtouch(xptr(rcx));			/* Preload first TLB */
		r8_h8cl_sixteen_reals_eight_complex_with_mulf(rsi, 8*64, 64, 2*64, 4*64);
		loops_init_prefetch(288, 64, 1, rcx, 1, 16);
		r8_x8cl_eight_complex_with_mulf_preload;
		for(unsigned int loopA = 18*16-1; loopA; ) {
			xtouch(xptr(rcx+4096-64));		/* Preload the TLBs */
			do{
				r8_x8cl_eight_complex_with_mulf(rsi, 8*64, 64, 2*64, 4*64);
			}while(--loopA & 15);
			rsi += -16*8*64+dist128;	/* Next source pointer */
		}
		//rsi += -18*dist128;	/* Next source pointer */
	}
	rsi = daddr;

	/* Do inverse FFT levels 9,10 */
	/* On input the 64-byte cache lines hold these data values: */
	/*	0	+9216	4	+9216	8	+9216	12	+9216 */
	/*	1	... */
	/*	... */
	/*	3	... */
	/*	16	... */
	/*	... */
	/* On output the 64-byte cache lines hold these data values: */
	/*	0	+9216	16	+9216	32	+4608	48	+9216 */
	/*	1	... */
	/*	... */
	/*	15	... */
	/*	64	... */
	/*	... */

	/* Do 8 eight_reals_unfft operations */
	/* Do 8 four_complex_djbunfft operations */
	/*	distance between fft data elements is 8 */
	/*	do 4 macros processing 32 data values */

	rdi = (uintptr_t)g->xsincos_complex;	/* Load sin/cos pointers */
	rdi += 1536*XMM_SCD1;
	rdx = (uintptr_t)g->sincos3;
	rdx += 1536*XMM_HSCD1;
	loops_init_prefetch(576, 64, 2, rcx);
	xtouch(xptr(rcx+4096-64));		/* Preload the TLBs */
	xtouch(xptr(rcx+2*4096-64));
	xtouch(xptr(rcx+3*4096-64));
	xtouch(xptr(rcx+4*4096-64));
	xtouch(xptr(rcx+5*4096-64));
	for(unsigned int loopA = 4; loopA; loopA--) {
		r4_h4cl_eight_reals_four_complex_djbunfft(rsi, 64, 4*64, 8*64, rdi, 4*16*XMM_SCD2, rdx, 4*16*XMM_HSCD2);
		rdi += 16*XMM_SCD2;	/* Next sin/cos data */
		rdx += 16*XMM_HSCD2;	/* Next sin/cos data */
	}
	rsi += -4*64+16*64;	/* Next source pointer */
	rdi += -4*16*XMM_SCD2;	/* Restore sin/cos ptr */
	rdx += -4*16*XMM_HSCD2;	/* Restore sin/cos ptr */

	/* Do 2288 four_complex_djbunfft operations */
	/*	distance between fft data elements is 8 */
	/*	do 572 macros each processing 32 data values */

	loops_reset;			/* 18 iters of 7 or 8 of 4 */
	r4_x4cl_four_complex_djbunfft_preload;
	for(unsigned int loopA = 18*8-1; loopA; ) {
		do{
			for(unsigned int loopC = 4; loopC; loopC--) {
				r4_x4cl_four_complex_djbunfft(rsi, 64, 4*64, 8*64, rdi, 4*16*XMM_SCD2);
				rdi += 16*XMM_SCD2;	/* Next sine/cosine pointer */
			}
			rsi += -4*64+16*64;	/* Next source pointer */
			rdi += -4*16*XMM_SCD2;	/* Restore sin/cos ptr */
		}while(--loopA & 7);
		rsi += -8*16*64+dist128;	/* Next source pointer */
	}
	rsi = daddr; //rsi += -18*dist128;	/* Next source pointer */

	/* Do inverse FFT levels 7,8 */
	/* On input the 64-byte cache lines hold these data values: */
	/*	0	+9216	16	+9216	32	+4608	48	+9216 */
	/*	1	... */
	/*	... */
	/*	15	... */
	/*	64	... */
	/*	... */
	/* On output the 64-byte cache lines hold these data values: */
	/*	0	+9216	64	+9216	128	+4608	192	+9216 */
	/*	1	... */
	/*	... */
	/*	63	... */
	/*	256	... */
	/*	... */

	/* Do 32 eight_reals_unfft operations */
	/* Do 32 four_complex_djbunfft operations */
	/*	distance between fft data elements is 32 */
	/*	do 16 macros processing 32 data values */

	loops_init_prefetch(576, 64, 2, rcx);
	xtouch(xptr(rcx+4096-64));		/* Preload the TLBs */
	xtouch(xptr(rcx+2*4096-64));
	xtouch(xptr(rcx+3*4096-64));
	xtouch(xptr(rcx+4*4096-64));
	xtouch(xptr(rcx+5*4096-64));
	for(unsigned int loopA = 16; loopA; loopA--) {
		r4_h4cl_eight_reals_four_complex_djbunfft(rsi, 64, 16*64, 32*64, rdi, 16*4*XMM_SCD2, rdx, 16*4*XMM_HSCD2);
		rdi += 4*XMM_SCD2;		/* Next sin/cos data */
		rdx += 4*XMM_HSCD2;	/* Next sin/cos data */
	}
	rsi += -16*64+64*64;	/* Next source pointer */
	rdi += -16*4*XMM_SCD2;	/* Restore sin/cos ptr */
	rdx += -16*4*XMM_HSCD2;	/* Restore sin/cos ptr */

	/* Do 2240 four_complex_djbunfft operations */
	/*	distance between fft data elements is 16 */
	/*	do 560 macros each processing 32 data values */

	loops_reset;			/* 18 iters of 1 or 2 of 16 */
	r4_x4cl_four_complex_djbunfft_preload;
	for(unsigned int loopA = 18*2-1; loopA; ) {
		do{
			for(unsigned int loopC = 16; loopC; loopC--) {
				r4_x4cl_four_complex_djbunfft(rsi, 64, 16*64, 32*64, rdi, 16*4*XMM_SCD2);
				rdi += 4*XMM_SCD2;		/* Next sin/cos data */
			}
			rsi += -16*64+64*64;	/* Next source pointer */
			rdi += -16*4*XMM_SCD2;	/* Restore sin/cos ptr */
		}while(--loopA & 1);
		rsi += -2*64*64+dist128;	/* Next source pointer */
	}
	rsi = daddr; //rsi += -18*dist128;	/* Next source pointer */

	/* Do inverse FFT levels 5,6 */
	/* On input the 64-byte cache lines hold these data values: */
	/*	0	+9216	64	+9216	128	+4608	192	+9216 */
	/*	1	... */
	/*	... */
	/*	63	... */
	/*	256	... */
	/*	... */
	/* On output the 64-byte cache lines hold these data values: */
	/*	0	+9216	256	+9216	512	+4608	768	+9216 */
	/*	1	... */
	/*	... */
	/*	255	... */
	/*	1024	... */
	/*	... */

	/* Do 128 eight_reals_unfft operations */
	/* Do 128 four_complex_djbunfft operations */
	/*	distance between fft data elements is 128 */
	/*	do 64 macros processing 32 data values */

	loops_init_prefetch(576, 64, 2, rcx);
	xtouch(xptr(rcx+4096-64));		/* Preload the TLBs */
	xtouch(xptr(rcx+2*4096-64));
	xtouch(xptr(rcx+3*4096-64));
	xtouch(xptr(rcx+4*4096-64));
	xtouch(xptr(rcx+5*4096-64));
	for(unsigned int loopA = 64; loopA; loopA--) {
		r4_h4cl_eight_reals_four_complex_djbunfft(rsi, 64, 64*64, dist128, rdi, 64*XMM_SCD2, rdx, 64*XMM_HSCD2);
		rdi += XMM_SCD2;		/* Next sin/cos data */
		rdx += XMM_HSCD2;		/* Next sin/cos data */
	}
	rsi += -64*64+2*dist128;	/* Next source pointer */
	rdi += -64*XMM_SCD2;	/* Restore sin/cos ptr */
	rdx += -64*XMM_HSCD2-1536*XMM_HSCD1; /* Restore sin/cos ptr */

	/* Do 2048 four_complex_djbunfft operations */
	/*	distance between fft data elements is 128 */
	/*	do 512 macros each processing 32 data values */

	loops_reset;			/* 8 more iters of 64 */
	r4_x4cl_four_complex_djbunfft_preload;
	for(unsigned int loopA = 8; loopA; loopA--) {
		for(unsigned int loopB = 64; loopB; loopB--) {
			r4_x4cl_four_complex_djbunfft(rsi, 64, 64*64, dist128, rdi, 64*XMM_SCD2);
			rdi += XMM_SCD2;		/* Next sin/cos data */
		}
		rsi += -64*64+2*dist128;	/* Next source pointer */
		rdi += -64*XMM_SCD2;	/* Restore sin/cos ptr */
	}
	rsi = daddr; //rsi += -9*2*dist128;	/* Next source pointer */

	/* Do inverse FFT levels 3,4 (actually 1.585 FFT levels) */
	/* On input the 64-byte cache lines hold these data values: */
	/*	0	+9216	256	+9216	512	+4608	768	+9216 */
	/*	1	... */
	/*	... */
	/*	255	... */
	/*	1024	... */
	/*	... */
	/* On output the 64-byte cache lines hold these data values: */
	/*	0	+9216	768	+9216	1536	+4608	2304	+9216 */
	/*	1	... */
	/*	... */
	/*	767	... */
	/*	3072	... */
	/*	... */

	/* Do 512 six_reals_unfft operations */
	/* Do 512 six_complex_djbunfft operations */
	/*	distance between fft data elements is 512 */
	/*	do 256 macros processing 32 data values */

	loops_init_prefetch(768, 64, 4, rcx);
	xtouch(xptr(rcx+4096-64));		/* Preload the TLBs */
	xtouch(xptr(rcx+2*4096-64));
	xtouch(xptr(rcx+3*4096-64));
	rdi = (uintptr_t)g->xsincos_complex;	/* Load sin/cos pointer */
	for(unsigned int loopA = 2; loopA; loopA--) {
		for(unsigned int loopB = 128; loopB; loopB--) {
			r3_h3cl_six_reals_three_complex_djbunfft(rsi, 64, 2*dist128, rdi, 256*3*XMM_SCD1, rdx, 256*3*XMM_HSCD1);
			rdi += 3*XMM_SCD1;		/* Next sin/cos data */
			rdx += 3*XMM_HSCD1;	/* Next sin/cos data */
		}
		rsi += -128*64+dist128;	/* Next source pointer */
	}
	rsi += -2*dist128+6*dist128; /* Next source pointer */
	rdi += -256*3*XMM_SCD1;	/* Restore sin/cos ptr */
	rdx += -256*3*XMM_HSCD1;	/* Restore sin/cos ptr */

	/* Do 2048 three_complex_djbunfft operations */
	/*	distance between fft data elements is 512 */
	/*	do 512 macros each processing 32 data values */

	loops_reset;			/* 2 more iters of 256 */
	r3_x3cl_three_complex_djbunfft_preload;
	for(unsigned int loopA = 2; loopA; loopA--) {
		for(unsigned int loopB = 2; loopB; loopB--) {
			for(unsigned int loopC = 128; loopC; loopC--) {
				r3_x3cl_three_complex_djbunfft(rsi, 64, 2*dist128, rdi, 256*3*XMM_SCD1);
				rdi += 3*XMM_SCD1;		/* Next sin/cos data */
			}
			rsi += -128*64+dist128;	/* Next source pointer */
		}
		rsi += -2*dist128+6*dist128; /* Next source pointer */
		rdi += -256*3*XMM_SCD1;	/* Restore sin/cos ptr */
	}
	rsi = daddr; //rsi += -3*6*dist128;	/* Next source pointer */

	/* Do inverse FFT level 1,2 (actually 1.585 FFT levels) */
	/* On input the 64-byte cache lines hold these data values: */
	/*	0	+9216	768	+9216	1536	+4608	2304	+9216 */
	/*	1	... */
	/*	... */
	/*	767	... */
	/*	3072	... */
	/*	... */
	/* On output the 64-byte cache lines hold these data values if pass 1 swizzles: */
	/*	0	+9216	2304	+9216	4608	+4608	6912	+9216 */
	/*	1	... */
	/*	... */
	/*	2303 */
	/* On output the 64-byte cache lines hold these data values if pass 2 swizzles: */
	/*	0	+2304	4608	+2304	9216	+2304	13824	+2304 */
	/*	1	... */
	/*	... */
	/*	2303 */

	/* Do 1536 six_reals_unfft operations */
	/* Do 1536 three_complex_djbunfft operations */
	/*	distance between fft data elements is 1536 */
	/*	do 768 macros processing 24 data values */

	pfing(rcx += -2*64);		/* Back up prefetch pointer by 2 cache lines */
	loops_init_prefetch(768, 64, 4, rcx);
	xtouch(xptr(rcx+4096-64));		/* Preload the TLBs */
	xtouch(xptr(rcx+2*4096-64));
	xtouch(xptr(rcx+3*4096-64));
	for(unsigned int loopA = 6; loopA; loopA--) {
		for(unsigned int loopB = 128; loopB; loopB--) {
			r3_h3cl_six_reals_three_complex_djbunfft(rsi, 64, 6*dist128, rdi, 768*XMM_SCD1, rdx, 768*XMM_HSCD1);
			rdi += XMM_SCD1;		/* Next sin/cos data */
			rdx += XMM_HSCD1;		/* Next sin/cos data */
		}
		rsi += -128*64+dist128;	/* Next source pointer */
	}

	xpass2_4608_real_done:;
	end_timer(2);
	if(g->pass2_get_next_block(g) != 0) return;
}
r4_pass2_4608_levels_complex(g);
}

	/* *************** radix4-pass2-4608-levels-complex macro ****************** */
	/* This macro takes 2 groups of 4608 complex values and */
	/* performs the final 12.17 levels of a standard mostly radix-4 FFT, squares the */
	/* results, and does the first 12.17 levels of the inverse FFT. */


void r4_pass2_4608_levels_complex(struct gwasm_data *__restrict g) {
do{
	uintptr_t rbx,rbp,rsi,rdi,daddr;
	vec2f64 xmm0,xmm1,xmm2,xmm3,xmm4,xmm5,xmm6,xmm7;
	/* Type 4 FFTs skip the forward FFT process */

	daddr = (uintptr_t)g->data_addr;/* Load address of FFT data */
	rsi = daddr;
	rbx = g->DIST_TO_FFTSRCARG;
	if likely(!(g->ffttype & 4)) { //xpass2_4608_levels_complex_4;

		/* Do FFT levels 1,2 (actually 1.585 levels) */
		/* */
		/* We swizzle the data so that one set of data is in the high word of */
		/* an xmm register, the other set is in the low word.  This means */
		/* a distance of 9216 data values. */
		/* */
		/* On input the 64-byte cache lines hold these data values if pass 2 swizzles: */
		/*	0	+2304	4608	+2304	9216	+2304	13824	+2304 */
		/*	1	... */
		/*	... */
		/*	2303 */
		/* On input the 64-byte cache lines hold these data values if pass 1 swizzles: */
		/*	0	+9216	2304	+9216	4608	+9216	6912	+9216 */
		/*	1	... */
		/*	... */
		/*	2303	... */
		/* On output the 64-byte cache lines hold these data values: */
		/*	0	+9216	768	+9216	1536	+9216	2304	+9216 */
		/*	1	... */
		/*	... */
		/*	767	... */
		/*	3072	... */
		/*	... */

		/* Do 3072 three_complex_fft operations */
		/*	distance between fft data elements is 1536 */
		/*	do 768 macros each processing 24 data values */

		/* Do FFT levels 3,4 (actually 1.585 levels) */
		/* */
		/* On input the 64-byte cache lines hold these data values: */
		/*	0	+9216	768	+9216	1536	+9216	2304	+9216 */
		/*	1	... */
		/*	... */
		/*	767	... */
		/*	3072	... */
		/*	... */
		/* On output the 64-byte cache lines hold these data values: */
		/*	0	+9216	256	+9216	512	+9216	768	+9216 */
		/*	1	... */
		/*	... */
		/*	255	... */
		/*	1024	... */
		/*	... */

		/* Do 2304 three_complex_fft operations */
		/*	distance between fft data elements is 512 */
		/*	do 768 macros each processing 24 data values */

#if 0 // @INSTR(,%xarch,<CORE>) NE 0; /* Core 2 faster doing separate loops */

		start_timer(4);
		rdi = (uintptr_t)g->xsincos_complex;	/* The sin/cos data */
		loops_init(768);			/* 6 iterations of 128 */
		r3_f3cl_three_complex_djbfft_preload;
		for(unsigned int loopA = 6; loopA; loopA--) {
			for(unsigned int loopB = 128; loopB; loopB--) {
				r3_f3cl_three_complex_djbfft(rsi, 64, 6*dist128, rdi, 768*XMM_SCD1);
				rdi += XMM_SCD1;		/* Next sin/cos data */
			}
			rsi += -128*64+dist128;	/* Next source pointer */
		}
		rsi = daddr; //rsi += -6*dist128;		/* Restore source pointer */
		end_timer(4);

		start_timer(5);
		pfing(rcx = (uintptr_t)g->data_prefetch); /*  Prefetch pointer (we need to prefetch  */
		;					/* 2304 cache lines + 34 pad cache lines) */
		xtouch(xptr(rcx));			/* Load first prefetch TLB */
		rdi = (uintptr_t)g->xsincos_complex;	/* Ptr to the sin/cos data */
		loops_init_prefetch(768, 64, 4, rcx); /* 3 iters of 256 */
		r3_nf3cl_three_complex_djbfft_preload;
		for(unsigned int loopA = 3; loopA; loopA--) {
			xtouch(xptr(rcx+4096-64));		/* Load prefetch TLBs */
			for(unsigned int loopB = 2; loopB; loopB--) {
				for(unsigned int loopC = 128; loopC; loopC--) {
					r3_nf3cl_three_complex_djbfft(rsi, 64, 2*dist128, rdi, 256*3*XMM_SCD1);
					rdi += 3*XMM_SCD1;		/* Next sin/cos data */
				}
				rsi += -128*64+dist128;	/* Next source pointer */
			}
			rsi += -2*dist128+6*dist128; /* Next source pointer */
			rdi += -256*3*XMM_SCD1;	/* Restore sin/cos ptr */
		}
		rsi = daddr; //rsi += -3*6*dist128;	/* Restore source pointer */
		end_timer(5);

#else /* Pentium 4 faster doing merged loops */

		start_timer(4);
		rdi = (uintptr_t)g->xsincos_complex;	/* The sin/cos data */
		rbp = rdi;
		pfing(rcx = (uintptr_t)g->data_prefetch); /*  Prefetch pointer (we need to prefetch  */
		;					/* 2304 cache lines + 34 pad cache lines) */
		xtouch(xptr(rcx));			/* Load first prefetch TLB */
		rdi = (uintptr_t)g->xsincos_complex;	/* Ptr to the sin/cos data */
		loops_init_prefetch(768, 64, 12, rcx); /* 256 iters of 3, prefetching 64 cache lines */

		r3_f3cl_three_complex_djbfft_preload;
		/*	r3_nf3cl_three_complex_djbfft_preload	;; Dangerous - assume same as r3_f3cl_three_complex_djbfft_preload */

		xtouch(xptr(rcx+4096-64));		/* Load prefetch TLBs */

		for(unsigned int loopA = 2; loopA; loopA--) {
			for(unsigned int loopB = 128; loopB; loopB--) {
				for(unsigned int loopC = 3; loopC; loopC--) {
					r3_f3cl_three_complex_djbfft(rsi, 2*dist128, 6*dist128, rdi, 768*XMM_SCD1);
					xxprefetch(u8ptr(rdi+XMM_SCD1));	/* Prefetch for next pass */
					xxprefetch(u8ptr(rdi+769*XMM_SCD1));	/* Prefetch for next pass */
					rdi += 256*XMM_SCD1;	/* Next sin/cos ptr */
				}
				rsi += -3*2*dist128;	/* Restore source pointer */
				loops_undo(3);

				for(unsigned int loopC = 3; loopC; loopC--) {
					r3_nf3cl_three_complex_djbfft(rsi, 6*dist128, 2*dist128, rbp, 256*3*XMM_SCD1);
				}

				rsi += -3*6*dist128+64;	/* Next source pointer */
				rdi += -3*256*XMM_SCD1+XMM_SCD1; /* Next sin/cos ptr */
				rbp += 3*XMM_SCD1;		/* Next sin/cos ptr */
			}
			rsi += -128*64+dist128;	/* Next source pointer */
		}
		rsi = daddr; //rsi += -2*dist128;		/* Restore source pointer */
		rdi += -256*XMM_SCD1;	/* Restore sin/cos ptr */
		end_timer(4);

#endif

		/* Do FFT levels 5,6 */
		/* */
		/* On input the 64-byte cache lines hold these data values: */
		/*	0	+9216	256	+9216	512	+9216	768	+9216 */
		/*	1	... */
		/*	... */
		/*	255	... */
		/*	1024	... */
		/*	... */
		/* On output the 64-byte cache lines hold these data values: */
		/*	0	+9216	128	+9216	256	+9216	384	+9216 */
		/*	1	... */
		/*	... */
		/*	127	... */
		/*	512	... */
		/*	... */
		/* */
		/* Do 2304 four_complex_fft operations */
		/*	distance between fft data elements is 128 */
		/*	do 1152 macros each processing 16 data values */

		/* Do FFT levels 7,8 */
		/* */
		/* On input the 64-byte cache lines hold these data values: */
		/*	0	+9216	128	+9216	256	+9216	384	+9216 */
		/*	1	... */
		/*	... */
		/*	127	... */
		/*	512	... */
		/*	... */
		/* On output the 64-byte cache lines hold these data values: */
		/*	0	+9216	32	+9216	64	+9216	96	+9216 */
		/*	1	... */
		/*	... */
		/*	31	... */
		/*	128	... */
		/*	... */
		/* */
		/* Do 2304 four_complex_fft operations */
		/*	distance between fft data elements is 32 */
		/*	do 576 macros each processing 32 data values */

		/* Do FFT levels 9,10 */
		/* */
		/* On input the 64-byte cache lines hold these data values: */
		/*	0	+9216	32	+9216	64	+9216	96	+9216 */
		/*	1	... */
		/*	... */
		/*	31	... */
		/*	128	... */
		/*	... */
		/* On output the 64-byte cache lines hold these data values: */
		/*	0	+9216	8	+9216	16	+9216	24	+9216 */
		/*	1	... */
		/*	... */
		/*	7	... */
		/*	32	... */
		/*	... */
		/* */
		/* Do 2304 four_complex_fft operations */
		/*	distance between fft data elements is 8 */
		/*	do 576 macros each processing 32 data values */

		/* Do this 32 cache lines at a time to maximize accessing data */
		/* while it is in the L1 cache. */

		start_timer(6);
		rdi += 1536*XMM_SCD1;	/* Load sin/cos pointers */
		rbx = rdi;
		rbp = rdi;
		loops_init_prefetch(1152, 64, 4, rcx); /* 9 * 8 iters of 32 cache lines */
		;	/*r4_nf2cl_four_complex_djbfft_preload -- assume same as r4_x4cl_four_complex_djbfft_preload */
		r4_x4cl_four_complex_djbfft_preload;

		for(unsigned int loopA = 9; loopA; loopA--) {
			xtouch(xptr(rcx+4096-64));		/* Load prefetch TLBs */
			xtouch(xptr(rcx+2*4096-64));

			for(unsigned int loopB = 8; loopB; loopB--) {
				for(unsigned int loopC = 16; loopC; loopC--) {
					xxprefetch(u8ptr(rdi+8*XMM_SCD2));	/* Prefetch for next pass */
					r4_nf2cl_four_complex_djbfft(rsi, 8*64, dist128, rdi);
					rdi += 8*XMM_SCD2;		/* Next sine/cosine pointer */
				}
				rsi += -16*8*64;		/* Restore source pointer */
				rdi += -16*8*XMM_SCD2;	/* Restore sine/cosine pointer */
				loops_undo(16);			/* Back up loop counter */

				for(unsigned int loopC = 2; loopC; loopC--) {
					for(unsigned int loopD = 4; loopD; loopD--) {
						r4_x4cl_four_complex_djbfft(rsi, 8*64, 32*64, 64*64, rbx);
						rbx += 4*8*XMM_SCD2;	/* Next sine/cosine pointer */
					}
					rsi += -4*8*64+dist128;	/* Next source pointer */
					rbx += -4*4*8*XMM_SCD2;	/* Restore sine/cosine pointer */
				}
				rsi += -2*dist128;		/* Restore source pointer */
				loops_undo(2*4*2);		/* Back up loop counter */

				for(unsigned int loopC = 2; loopC; loopC--) {
					for(unsigned int loopD = 4; loopD; loopD--) {
						r4_x4cl_four_complex_djbfft(rsi, 32*64, 8*64, 16*64, rbp);
					}
					rsi += -4*32*64+dist128;	/* Next source pointer */
				}

				rsi += -2*dist128+64;	/* Next source pointer */
				rdi += XMM_SCD2;		/* Next sine/cosine pointer */
				rbx += 4*XMM_SCD2;		/* Next sine/cosine pointer */
				rbp += 16*XMM_SCD2;	/* Next sine/cosine pointer */
			}
			rsi += -8*64+2*dist128;	/* Next source pointer */
			rdi += -8*XMM_SCD2;	/* Restore sine/cosine pointer */
			rbx += -8*4*XMM_SCD2;	/* Restore sine/cosine pointer */
			rbp += -8*16*XMM_SCD2;	/* Restore sine/cosine pointer */
		}
		rsi = daddr; //rsi += -9*2*dist128;	/* Restore source pointer */
		end_timer(6);

		/* Do FFT levels 11,12,13 as well as inverse FFT levels 11,12,13 */
		/* */
		/* On input the 64-byte cache lines hold these data values: */
		/*	0	+9216	8	+9216	16	+9216	24	+9216 */
		/*	1	... */
		/*	... */
		/*	7	... */
		/*	32	... */
		/*	... */
		/* On output the 64-byte cache lines hold these data values: */
		/*	0	+9216	4	+9216	8	+9216	12	+9216 */
		/*	1	... */
		/*	... */
		/*	3	... */
		/*	16	... */
		/*	... */

		/* Do 1152 eight_complex_fft operations */
		/*	distance between fft data elements is 1 */
		/*	do 288 macros each processing 64 data values */

		/* Execute the right middle step */

		if likely(g->ffttype & 2) {
			if likely(!(g->ffttype & 1)) { // xpass2_4608_levels_complex_2:;
				start_timer(9);
				loops_init_prefetch(288, 64, 8, rcx);
				xtouch(xptr(rcx+4096-64));		/* Load prefetch TLB */
				r8_x8cl_eight_complex_with_square_preload;
				for(unsigned int loopA = 18; loopA; loopA--) {
					for(unsigned int loopB = 16; loopB; loopB--) {
						r8_x8cl_eight_complex_with_square(rsi, 8*64, 64, 2*64, 4*64);
					}
					rsi += -16*8*64+dist128;	/* Next source pointer */
				}
				//rsi += -18*dist128;	/* Restore source pointer */
				end_timer(9);
			}else{ // xpass2_4608_levels_complex_3:;
				start_timer(9);
				rbp = g->DIST_TO_MULSRCARG;
				loops_init_prefetch(288, 64, 8, rcx);
				r8_x8cl_eight_complex_with_mult_preload;
				xtouch(xptr(rcx+4096-64));		/* Load prefetch TLB */
				for(unsigned int loopA = 18; loopA; loopA--) {
					for(unsigned int loopB = 16; loopB; loopB--) {
						r8_x8cl_eight_complex_with_mult(rsi, 8*64, 64, 2*64, 4*64);
					}
					rsi += -16*8*64+dist128;	/* Next source pointer */
				}
				//rsi += -18*dist128;	/* Restore source pointer */
				end_timer(9);
			}
		}else{ // xpass2_4608_levels_complex_1:;
			start_timer(9);
			loops_init_prefetch(288, 128, 1, rcx);
			r8_x8cl_eight_complex_fft_final_preload;
			for(unsigned int loopA = 18; loopA; loopA--) {
				xtouch(xptr(rcx+4096-64));		/* Load prefetch TLBs */
				for(unsigned int loopB = 16; loopB; loopB--) {
					r8_x8cl_eight_complex_fft_final(rsi, 8*64, 64, 2*64, 4*64);
				}
				rsi += -16*8*64+dist128;	/* Next source pointer */
			}
			end_timer(9);
			goto	xpass2_4608_levels_complex_done;
		}
	}else{ // xpass2_4608_levels_complex_4:;
		start_timer(9);
		pfing(rcx = (uintptr_t)g->data_prefetch); /* Prefetch pointer */
		xtouch(xptr(rcx));			/* Load first prefetch TLB */
		rbp = g->DIST_TO_MULSRCARG;
		loops_init_prefetch(288, 128, 1, rcx);
		r8_x8cl_eight_complex_with_mulf_preload;
		for(unsigned int loopA = 18; loopA; loopA--) {
			xtouch(xptr(rcx+4096-64));		/* Load prefetch TLBs */
			for(unsigned int loopB = 16; loopB; loopB--) {
				r8_x8cl_eight_complex_with_mulf(rsi, 8*64, 64, 2*64, 4*64);
			}
			rsi += -16*8*64+dist128;	/* Next source pointer */
		}
		//rsi += -18*dist128;	/* Restore source pointer */
		end_timer(9);
	}
	rsi = daddr;

	/* Do inverse FFT levels 9,10 */
	/* On input the 64-byte cache lines hold these data values: */
	/*	0	+9216	4	+9216	8	+9216	12	+9216 */
	/*	1	... */
	/*	... */
	/*	3	... */
	/*	16	... */
	/*	... */
	/* On output the 64-byte cache lines hold these data values: */
	/*	0	+9216	16	+9216	32	+4608	48	+9216 */
	/*	1	... */
	/*	... */
	/*	15	... */
	/*	64	... */
	/*	... */
	/* */
	/* Do 2304 four_complex_unfft operations */
	/*	distance between fft data elements is 8 */
	/*	do 576 macros each processing 32 data values */

	/* Do inverse FFT levels 7,8 */
	/* On input the 64-byte cache lines hold these data values: */
	/*	0	+9216	16	+9216	32	+4608	48	+9216 */
	/*	1	... */
	/*	... */
	/*	15	... */
	/*	64	... */
	/*	... */
	/* On output the 64-byte cache lines hold these data values: */
	/*	0	+9216	64	+9216	128	+4608	192	+9216 */
	/*	1	... */
	/*	... */
	/*	63	... */
	/*	256	... */
	/*	... */
	/* */
	/* Do 2304 four_complex_unfft operations */
	/*	distance between fft data elements is 32 */
	/*	do 576 macros each processing 32 data values */

	/* Do inverse FFT levels 5,6 */
	/* On input the 64-byte cache lines hold these data values: */
	/*	0	+9216	64	+9216	128	+4608	192	+9216 */
	/*	1	... */
	/*	... */
	/*	63	... */
	/*	256	... */
	/*	... */
	/* On output the 64-byte cache lines hold these data values: */
	/*	0	+9216	256	+9216	512	+4608	768	+9216 */
	/*	1	... */
	/*	... */
	/*	255	... */
	/*	1024	... */
	/*	... */
	/* */
	/* Do 2304 four_complex_unfft operations */
	/*	distance between fft data elements is 128 */
	/*	do 576 macros each processing 32 data values */

	/* Do this 64 cache lines at a time to maximize accessing data */
	/* while it is in the L1 cache. */

	start_timer(10);
	rdi = (uintptr_t)g->xsincos_complex;	/* Load sin/cos pointers */
	rdi += 1536*XMM_SCD1;
	rbx = rdi;
	rbp = rdi;
	loops_init_prefetch(576, 64, 2, rcx); /* 9*4 iters of 64 cache lines */
	r4_x4cl_four_complex_djbunfft_preload;

	for(unsigned int loopA = 9; loopA; loopA--) {
		xtouch(xptr(rcx+4096-128));		/* Load prefetch TLBs */
		xtouch(xptr(rcx+2*4096-128));
		for(unsigned int loopB = 4; loopB; loopB--) {

			for(unsigned int loopC = 2; loopC; loopC--) {
				for(unsigned int loopD = 8; loopD; loopD--) {
					r4_x4cl_four_complex_djbunfft(rsi, 16*64, 4*64, 8*64, rbp, 4*16*XMM_SCD2);
				}
				rsi += -8*16*64+dist128;	/* Next source pointer */
			}
			rsi += -2*dist128;		/* Restore source pointer */
			loops_undo(2*8);			/* Restore loop counter */

			for(unsigned int loopC = 2; loopC; loopC--) {
				for(unsigned int loopD = 2; loopD; loopD--) {
					for(unsigned int loopE = 4; loopE; loopE--) {
						r4_x4cl_four_complex_djbunfft(rsi, 4*64, 16*64, 32*64, rbx, 16*4*XMM_SCD2);
						rbx += 4*4*XMM_SCD2;	/* Next sine/cosine pointer */
					}
					rsi += -4*4*64+64*64;	/* Next source pointer */
					rbx += -4*4*4*XMM_SCD2;	/* Restore sine/cosine pointer */
				}
				rsi += -2*64*64+dist128;	/* Next source pointer */
			}
			rsi += -2*dist128;		/* Restore source pointer */
			loops_undo(2*2*4);		/* Restore loop counter */

			for(unsigned int loopC = 16; loopC; loopC--) {
				r4_x4cl_four_complex_djbunfft(rsi, 4*64, 64*64, dist128, rdi, 64*1*XMM_SCD2);
				rdi += 4*1*XMM_SCD2;	/* Next sine/cosine pointer */
			}
			rdi += -16*4*1*XMM_SCD2;	/* Restore sine/cosine pointer */

			rsi += -16*4*64+64;	/* Next source pointer */
			rdi += XMM_SCD2;		/* Next sine/cosine pointer */
			rbx += 4*XMM_SCD2;		/* Next sine/cosine pointer */
			rbp += 16*XMM_SCD2;	/* Next sine/cosine pointer */
		}
		rsi += -4*64+2*dist128;	/* Next source pointer */
		rdi += -4*XMM_SCD2;	/* Restore sine/cosine pointer */
		rbx += -4*4*XMM_SCD2;	/* Restore sine/cosine pointer */
		rbp += -4*16*XMM_SCD2;	/* Restore sine/cosine pointer */
	}
	rsi = daddr; //rsi += -9*2*dist128;	/* Restore source pointer */
	end_timer(10);

	/* Do inverse FFT levels 3,4 (actually 1.585 FFT levels) */
	/* On input the 64-byte cache lines hold these data values: */
	/*	0	+9216	256	+9216	512	+4608	768	+9216 */
	/*	1	... */
	/*	... */
	/*	255	... */
	/*	1024	... */
	/*	... */
	/* On output the 64-byte cache lines hold these data values: */
	/*	0	+9216	768	+9216	1536	+4608	2304	+9216 */
	/*	1	... */
	/*	... */
	/*	767	... */
	/*	3072	... */
	/*	... */

	/* Do 3072 three_complex_unfft operations */
	/*	distance between fft data elements is 512 */
	/*	do 768 macros each processing 32 data values */

	/* Do inverse FFT levels 1,2 (actually 1.585 FFT levels) */
	/* On input the 64-byte cache lines hold these data values: */
	/*	0	+9216	768	+9216	1536	+4608	2304	+9216 */
	/*	1	... */
	/*	... */
	/*	767	... */
	/*	3072	... */
	/*	... */
	/* On output the 64-byte cache lines hold these data values if pass 1 swizzles: */
	/*	0	+9216	2304	+9216	4608	+4608	6912	+9216 */
	/*	1	... */
	/*	... */
	/*	2303 */
	/* On output the 64-byte cache lines hold these data values if pass 2 swizzles: */
	/*	0	+2304	4608	+2304	9216	+2304	13824	+2304 */
	/*	1	... */
	/*	... */
	/*	2303 */

	/* Do 3072 three_complex_unfft operations */
	/*	distance between fft data elements is 1536 */
	/*	do 768 macros each processing 24 data values */

#if 0 // @INSTR(,%xarch,<CORE>) NE 0; /* Core 2 faster doing separate loops */

	start_timer(13);
	rdi = (uintptr_t)g->xsincos_complex;	/* Load sin/cos pointer */
	loops_init_prefetch(768, 64, 4, rcx);/* 3 iters of 128 */
	r3_x3cl_three_complex_djbunfft_preload;
	for(unsigned int loopA = 3; loopA; loopA--) {
		xtouch(xptr(rcx+4096-64));		/* Load prefetch TLBs */
		for(unsigned int loopB = 2; loopB; loopB--) {
			for(unsigned int loopC = 128; loopC; loopC--) {
				r3_x3cl_three_complex_djbunfft(rsi, 64, 2*dist128, rdi, 256*3*XMM_SCD1);
				rdi += 3*XMM_SCD1;		/* Next sine/cosine pointer */
			}
			rsi += -128*64+dist128;	/* Next source pointer */
		}
		rsi += -2*dist128+6*dist128; /* Next source pointer */
		rdi += -256*3*XMM_SCD1;	/* Reset sine/cosine pointer */
	}
	rsi = daddr; //rsi += -3*6*dist128;	/* Restore source pointer */
	end_timer(13);

	start_timer(14);
	rdi = (uintptr_t)g->xsincos_complex;	/* Load sin/cos pointer */
	pfing(rcx += -2*64);		/* Back up prefetch pointer by 2 cache lines */
	loops_init_prefetch(768, 64, 4, rcx);/* 6 iters of 128 */
	r3_x3cl_three_complex_djbunfft_preload;
	for(unsigned int loopA = 6; loopA; loopA--) {
		xtouch(xptr(rcx+4096-64));		/* Load prefetch TLBs */
		for(unsigned int loopB = 128; loopB; loopB--) {
			r3_x3cl_three_complex_djbunfft(rsi, 64, 6*dist128, rdi, 768*XMM_SCD1);
			rdi += XMM_SCD1;		/* Next set of sin/cos data */
		}
		rsi += -128*64+dist128;	/* Next source pointer */
	}
	end_timer(14);

#else /* Pentium 4 faster doing merged loops */

	start_timer(13);
	rdi = (uintptr_t)g->xsincos_complex;	/* Load sin/cos pointer */
	rbp = rdi;
	pfing(rcx += -2*64);		/* Back up prefetch pointer by 2 cache lines */
	loops_init_prefetch(768, 64, 3, rcx); /* 256 iters of 3, prefetching 256 cache lines each loop */
	r3_x3cl_three_complex_djbunfft_preload;

	for(unsigned int loopA = 2; loopA; loopA--) {
		for(unsigned int loopB = 128; loopB; loopB--) {
			xtouch(xptr(rcx+4096-64));		/* Load prefetch TLBs */

			for(unsigned int loopC = 3; loopC; loopC--) {
				r3_x3cl_three_complex_djbunfft(rsi, 6*dist128, 2*dist128, rbp, 256*3*XMM_SCD1);
			}
			rsi += -3*6*dist128;	/* Restore source pointer */
			loops_undo(3);

			for(unsigned int loopC = 3; loopC; loopC--) {
				r3_x3cl_three_complex_djbunfft(rsi, 2*dist128, 6*dist128, rdi, 768*XMM_SCD1);
				rdi += 256*XMM_SCD1;	/* Next sin/cos ptr */
			}

			rsi += -3*2*dist128+64;	/* Next source pointer */
			rdi += -3*256*XMM_SCD1+XMM_SCD1; /* Next sin/cos ptr */
			rbp += 3*XMM_SCD1;		/* Next sine/cosine pointer */
		}
		rsi += -128*64+dist128;	/* Next source pointer */
	}
	end_timer(13);

#endif

	xpass2_4608_levels_complex_done:;
}while(g->pass2_get_next_block(g) == 0);
}
