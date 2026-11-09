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

#include "r4.h"
#include "r8.h"

/* *************** radix4-pass2-11-levels-real macro ****************** */
/* This macro takes 4096 real values and performs the final 12 levels of */
/* the FFT.  Independently, this macro also takes 2048 complex values and */
/* performs the final 11 levels of the FFT.  The results are squared, and */
/* the first 11/12 levels of the inverse FFT are done. */

void r4_pass2_11_levels_complex(struct gwasm_data *__restrict g);
void xpass2_r4_11_levels(struct gwasm_data *__restrict g) {
	g->thread_work_routine = r4_pass2_11_levels_complex;
	pass2_wake_up_threads(g);
	if(g->NEGACYCLIC_FFT == 0) {

//static inline void r4_pass2_11_levels_real(struct gwasm_data *__restrict g) {
	uintptr_t rdx,rbx,rbp,rsi,rdi,daddr;
	vec2f64 xmm0,xmm1,xmm2,xmm3,xmm4,xmm5,xmm6,xmm7;

	/* We switch to the same format used in one-pass FFTs (that is, the */
	/* high word of xmm registers being 1 greater than the low word): */

	/* Type 4 FFTs skip the forward FFT process */

	start_timer(2);
	daddr = rsi = (uintptr_t)g->data_addr;/* Load source address */
	rbx = g->DIST_TO_FFTSRCARG;
	if likely(!(g->ffttype & 4)) {

		/* Do FFT level 1,2 */
		/* */
		/* On input the 64-byte cache lines hold these data values if pass 2 swizzles: */
		/*	0	+1K	2048	+1K	4096	+1K	6144	+1K */
		/*	1	... */
		/*	... */
		/*	1023 */
		/* On input the 64-byte cache lines hold these data values if pass 1 swizzles: */
		/*	0	+4K	1024	+4K	2048	+4K	3072	+4K */
		/*	1	... */
		/*	... */
		/*	1023	... */
		/* On output the 64-byte cache lines hold these data values: */
		/*	0	+4K	512	+4K	1024	+4K	1536	+4K */
		/*	1	... */
		/*	... */
		/*	511	... */
		/*	2048	... */
		/*	... */

		/* Do 512 eight_reals_fft operations */
		/* Do 512 four_complex_djbfft operations */
		/*	distance between fft data elements is 512 */
		/*	do 512 macros each processing 16 data values */

		rdi = (uintptr_t)g->xsincos_complex;	/* Ptr to the sin/cos data */
		rdx = (uintptr_t)g->sincos3;/* Ptr to real sin/cos data */
		loops_init(512);
		for(unsigned int loopA = 4; loopA; loopA--) {
			for(unsigned int loopB = 128; loopB; loopB--) {
				r4_fh2cl_eight_reals_four_complex_djbfft(rsi, 64, 4*dist128, rdi, rdx);
				rdi += XMM_SCD2;		/* Next sin/cos data */
				rdx += XMM_HSCD2;		/* Next sin/cos data */
			}
			rsi += -128*64+dist128;	/* Next source pointer */
		}
		rsi = daddr; //rsi += -4*dist128;		/* Restore source pointer */
		rdx += -512*XMM_HSCD2;	/* Restore sin/cos ptr */

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

		/* Do 128 eight_reals_fft operations */
		/* Do 384 four_complex_djbfft operations */
		/*	distance between fft data elements is 128 */
		/*	do 128 macros each processing 32 data values */

		pfing(rcx = (uintptr_t)g->data_prefetch); /* Prefetch pointer */
		xtouch(xptr(rcx));			/* Preload the first TLB */
		loops_init_prefetch(256, 64, 2, rcx);
		xtouch(xptr(rcx+4096-64));		/* Preload the TLBs */
		xtouch(xptr(rcx+2*4096-64));
		r4_h4cl_eight_reals_four_complex_djbfft_preload;
		for(unsigned int loopA = 128; loopA; loopA--) {
			r4_h4cl_eight_reals_four_complex_djbfft(rsi, 64, dist128, 2*dist128, rdi, rdx);
			rdi += XMM_SCD2;		/* Next sin/cos data */
			rdx += 4*XMM_HSCD2;	/* Next sin/cos data */
		}
		rsi += -128*64+4*dist128;	/* Next source pointer */
		rdi += -128*XMM_SCD2;	/* Restore sin/cos ptr */
		rdx += -128*4*XMM_HSCD2;	/* Restore sin/cos ptr */

		/* Do 512 four_complex_djbfft operations */
		/*	distance between fft data elements is 128 */
		/*	do 128 macros each processing 32 data values */

		loops_reset;			/* 1 more iterations of 128 */
		r4_x4cl_four_complex_djbfft_preload;
		{
			for(unsigned int loopB = 128; loopB; loopB--) {
				r4_x4cl_four_complex_djbfft(rsi, 64, dist128, 2*dist128, rdi);
				rdi += XMM_SCD2;		/* Next sin/cos data */
			}
			rsi += -128*64+4*dist128;	/* Next source pointer */
			rdi += -128*XMM_SCD2;	/* Restore sin/cos ptr */
		}
		rsi = daddr; //rsi += -2*4*dist128;	/* Next source pointer */
		pfing(rcx += 128);		/* Skip pad bytes */

		/* Do FFT levels 5,6 */
		/* */
		/* On input the 64-byte cache lines hold these data values: */
		/*	0	+4K	128	+4K	256	+4K	384	+4K */
		/*	1	... */
		/*	... */
		/*	127	... */
		/*	512	... */
		/*	... */
		/* On output the 64-byte cache lines hold these data values: */
		/*	0	+4K	32	+4K	64	+4K	96	+4K */
		/*	1	... */
		/*	... */
		/*	31	... */
		/*	128	... */
		/*	... */

		/* Do 32 eight_reals_fft operations */
		/* Do 96 four_complex_djbfft operations */
		/*	distance between fft data elements is 32 */
		/*	do 32 macros each processing 32 data values */

		loops_init_prefetch(256, 64, 2, rcx);
		xtouch(xptr(rcx+4096-64));		/* Preload the TLBs */
		xtouch(xptr(rcx+2*4096-64));
		r4_h4cl_eight_reals_four_complex_djbfft_preload;
		for(unsigned int loopA = 32; loopA; loopA--) {
			r4_h4cl_eight_reals_four_complex_djbfft(rsi, 64, 32*64, 64*64, rdi, rdx);
			rdi += 4*XMM_SCD2;		/* Next sin/cos data */
			rdx += 16*XMM_HSCD2;	/* Next sin/cos data */
		}
		rsi += -32*64+dist128;	/* Next source pointer */
		rdi += -32*4*XMM_SCD2;	/* Restore sin/cos ptr */
		rdx += -32*16*XMM_HSCD2;	/* Restore sin/cos ptr */

		/* Do 896 four_complex_djbfft operations */
		/*	distance between fft data elements is 32 */
		/*	do 224 macros each processing 32 data values */

		loops_reset;			/* 7 more iters of 32 */
		r4_x4cl_four_complex_djbfft_preload;
		for(unsigned int loopA = 7; loopA; loopA--) {
			for(unsigned int loopB = 32; loopB; loopB--) {
				r4_x4cl_four_complex_djbfft(rsi, 64, 32*64, 64*64, rdi);
				rdi += 4*XMM_SCD2;		/* Next sin/cos data */
			}
			rsi += -32*64+dist128;	/* Next source pointer */
			rdi += -32*4*XMM_SCD2;	/* Restore sin/cos ptr */
		}
		rsi = daddr; //rsi += -8*dist128;		/* Next source pointer */
		pfing(rcx += 128);		/* Skip pad bytes */

		/* Do FFT levels 7,8 */
		/* */
		/* On input the 64-byte cache lines hold these data values: */
		/*	0	+4K	32	+4K	64	+4K	96	+4K */
		/*	1	... */
		/*	... */
		/*	31	... */
		/*	128	... */
		/*	... */
		/* On output the 64-byte cache lines hold these data values: */
		/*	0	+4K	8	+4K	16	+4K	24	+4K */
		/*	1	... */
		/*	... */
		/*	7	... */
		/*	32	... */
		/*	... */

		/* Do 8 eight_reals_fft operations */
		/* Do 24 four_complex_djbfft operations */
		/*	distance between fft data elements is 8 */
		/*	do 8 macros each processing 32 data values */

		loops_init_prefetch(256, 64, 2, rcx);
		xtouch(xptr(rcx+4096-64));		/* Preload the TLBs */
		xtouch(xptr(rcx+2*4096-64));
		r4_h4cl_eight_reals_four_complex_djbfft_preload;
		for(unsigned int loopA = 8; loopA; loopA--) {
			r4_h4cl_eight_reals_four_complex_djbfft(rsi, 64, 8*64, 16*64, rdi, rdx);
			rdi += 16*XMM_SCD2;	/* Next sin/cos data */
			rdx += 64*XMM_HSCD2;	/* Next sin/cos data */
		}
		rsi += -8*64+32*64;	/* Next source pointer */
		rdi += -8*16*XMM_SCD2;	/* Restore sin/cos ptr */
		rdx += -8*64*XMM_HSCD2;	/* Restore sin/cos ptr */

		/* Do 992 four_complex_djbfft operations */
		/*	distance between fft data elements is 8 */
		/*	do 248 macros each processing 32 data values */

		loops_reset;			/* 8 iters of 3 or 4 of 8 */
		r4_x4cl_four_complex_djbfft_preload;
		for(unsigned int loopA = 8*4-1; loopA; ) {
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
		rsi = daddr; //rsi += -8*dist128;		/* Next source pointer */
		pfing(rcx += 128);		/* Skip pad bytes */

		/* Do FFT levels 9,10,11 */
		/* */
		/* On input the 64-byte cache lines hold these data values: */
		/*	0	+4K	8	+4K	16	+4K	24	+4K */
		/*	1	... */
		/*	... */
		/*	7	... */
		/*	32	... */
		/*	... */
		/* On output the 64-byte cache lines hold these data values: */
		/*	0	+4K	4	+4K	8	+4K	12	+4K */
		/*	1	... */
		/*	... */
		/*	3	... */
		/*	16	... */
		/*	... */

		/* Do 1 sixteen_reals_fft operations */
		/* Do 3 eight_complex_fft operations */
		/*	distance between fft data elements is 1 */

		/* Do 508 eight_complex_fft operations */
		/*	distance between fft data elements is 1 */
		/*	do 127 macros each processing 64 data values */

		/* Execute the proper middle step */

		if likely(g->ffttype & 2) {
			if likely(!(g->ffttype & 1)) { //xpass2_11_levels_real_2:;
				r8_h8cl_sixteen_reals_eight_complex_with_square(rsi, 8*64, 64, 2*64, 4*64);
				loops_init_prefetch(128, 64, 1, rcx, 1, 16);
				r8_x8cl_eight_complex_with_square_preload;
				xtouch(xptr(rcx+4096-128));		/* Preload the TLBs */
				xtouch(xptr(rcx+2*4096-128));	/* Preload the TLBs */
				for(unsigned int loopA = 8*16-1; loopA; ) {
					do{
						r8_x8cl_eight_complex_with_square(rsi, 8*64, 64, 2*64, 4*64);
					}while(--loopA & 15);
					rsi += -16*8*64+dist128;	/* Next source pointer */
				}
				//rsi += -8*dist128;		/* Next source pointer */
				pfing(rcx += 128+128);		/* Skip pad bytes and the one missed prefetch */
			}else{ //xpass2_11_levels_real_3:;
				rbp = g->DIST_TO_MULSRCARG;
				r8_h8cl_sixteen_reals_eight_complex_with_mult(rsi, 8*64, 64, 2*64, 4*64);
				loops_init_prefetch(128, 64, 1, rcx, 1, 16);
				r8_x8cl_eight_complex_with_mult_preload;
				xtouch(xptr(rcx+4096-64));		/* Preload the TLBs */
				xtouch(xptr(rcx+2*4096-64));
				for(unsigned int loopA = 8*16-1; loopA; ) {
					do{
						r8_x8cl_eight_complex_with_mult(rsi, 8*64, 64, 2*64, 4*64);
					}while(--loopA & 15);
					rsi += -16*8*64+dist128;	/* Next source pointer */
				}
				//rsi += -8*dist128;		/* Next source pointer */
				pfing(rcx += 128+128);		/* Skip pad bytes and the one missed prefetch */
			}
		}else{ //xpass2_11_levels_real_1:;
			r8_h8cl_sixteen_reals_eight_complex_fft_final(rsi, 8*64, 64, 2*64, 4*64);
			loops_init_prefetch(128, 128, 1, rcx, 1, 16);
			r8_x8cl_eight_complex_fft_final_preload;
			for(unsigned int loopA = 8*16-1; loopA; ) {
				xtouch(xptr(rcx+4096-64));		/* Load prefetch TLBs */
				do{
					r8_x8cl_eight_complex_fft_final(rsi, 8*64, 64, 2*64, 4*64);
				}while(--loopA & 15);
				rsi += -16*8*64+dist128;	/* Next source pointer */
			}
			goto	xpass2_11_real_done;
		}
	}else{ //xpass2_11_levels_real_4:;
		rbp = g->DIST_TO_MULSRCARG;
		pfing(rcx = (uintptr_t)g->data_prefetch);/* Load prefetch pointer */
		xtouch(xptr(rcx));			/* Preload the first TLB */
		r8_h8cl_sixteen_reals_eight_complex_with_mulf(rsi, 8*64, 64, 2*64, 4*64);
		loops_init_prefetch(128, 128, 1, rcx, 1, 16);
		r8_x8cl_eight_complex_with_mulf_preload;
		for(unsigned int loopA = 8*16-1; loopA; ) {
			xtouch(xptr(rcx+4096-64));		/* Load prefetch TLBs */
			do{
				r8_x8cl_eight_complex_with_mulf(rsi, 8*64, 64, 2*64, 4*64);
			}while(--loopA & 15);
			rsi += -16*8*64+dist128;	/* Next source pointer */
		}
		//rsi += -8*dist128;		/* Next source pointer */
		pfing(pfing(rcx += 2*128);+128);		/* Skip pad bytes and the one missed prefetch */
	}
	rsi = daddr;

	/* Do inverse FFT levels 7,8 */
	/* On input the 64-byte cache lines hold these data values: */
	/*	0	+4K	4	+4K	8	+4K	12	+4K */
	/*	1	... */
	/*	... */
	/*	3	... */
	/*	16	... */
	/*	... */
	/* On output the 64-byte cache lines hold these data values: */
	/*	0	+4K	16	+4K	32	+4K	48	+4K */
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
	rdi += 512*XMM_SCD2;
	rdx = (uintptr_t)g->sincos3;
	loops_init_prefetch(256, 64, 2, rcx);
	xtouch(xptr(rcx+4096-128));
	xtouch(xptr(rcx+2*4096-128));
	for(unsigned int loopA = 4; loopA; loopA--) {
		r4_h4cl_eight_reals_four_complex_djbunfft(rsi, 64, 4*64, 8*64, rdi, 4*16*XMM_SCD2, rdx, 4*64*XMM_HSCD2);
		rdi += 16*XMM_SCD2;	/* Next sin/cos data */
		rdx += 64*XMM_HSCD2;	/* Next sin/cos data */
	}
	rsi += -4*64+16*64;	/* Next source pointer */
	rdi += -4*16*XMM_SCD2;	/* Restore sin/cos ptr */
	rdx += -4*64*XMM_HSCD2;	/* Restore sin/cos ptr */

	/* Do 1008 four_complex_djbunfft operations */
	/*	distance between fft data elements is 8 */
	/*	do 252 macros each processing 32 data values */

	loops_reset;			/* 8 iters of 7 or 8 of 4 */
	r4_x4cl_four_complex_djbunfft_preload;
	for(unsigned int loopA = 8*8-1; loopA; ) {
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
	rsi = daddr; //rsi += -8*dist128;	/* Next source pointer */
	pfing(rcx += 128);		/* Skip pad bytes */

	/* Do inverse FFT levels 5,6 */
	/* On input the 64-byte cache lines hold these data values: */
	/*	0	+4K	16	+4K	32	+4K	48	+4K */
	/*	1	... */
	/*	... */
	/*	15	... */
	/*	64	... */
	/*	... */
	/* On output the 64-byte cache lines hold these data values: */
	/*	0	+4K	64	+4K	128	+4K	192	+4K */
	/*	1	... */
	/*	... */
	/*	63	... */
	/*	256	... */
	/*	... */

	/* Do 32 eight_reals_unfft operations */
	/* Do 32 four_complex_djbunfft operations */
	/*	distance between fft data elements is 32 */
	/*	do 16 macros processing 32 data values */

	loops_init_prefetch(256, 64, 2, rcx);
	xtouch(xptr(rcx+4096-64));		/* Preload TLBs */
	xtouch(xptr(rcx+2*4096-64));
	for(unsigned int loopA = 16; loopA; loopA--) {
		r4_h4cl_eight_reals_four_complex_djbunfft(rsi, 64, 16*64, 32*64, rdi, 16*4*XMM_SCD2, rdx, 16*16*XMM_HSCD2);
		rdi += 4*XMM_SCD2;		/* Next sin/cos data */
		rdx += 16*XMM_HSCD2;	/* Next sin/cos data */
	}
	rsi += -16*64+64*64;	/* Next source pointer */
	rdi += -16*4*XMM_SCD2;	/* Restore sin/cos ptr */
	rdx += -16*16*XMM_HSCD2;	/* Restore sin/cos ptr */

	/* Do 960 four_complex_djbunfft operations */
	/*	distance between fft data elements is 32 */
	/*	do 240 macros each processing 32 data values */

	loops_reset;			/* 8 iters of 1 or 2 of 16 */
	r4_x4cl_four_complex_djbunfft_preload;
	for(unsigned int loopA = 8*2-1; loopA; ) {
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
	rsi = daddr; //rsi += -8*dist128;		/* Next source pointer */
	pfing(rcx += 128);		/* Skip pad bytes */

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

	/* Do 128 eight_reals_unfft operations */
	/* Do 128 four_complex_djbunfft operations */
	/*	distance between fft data elements is 128 */
	/*	do 64 macros processing 32 data values */

	loops_init_prefetch(256, 64, 2, rcx);
	xtouch(xptr(rcx+4096-128));		/* Preload TLBs */
	xtouch(xptr(rcx+2*4096-128));
	for(unsigned int loopA = 64; loopA; loopA--) {
		r4_h4cl_eight_reals_four_complex_djbunfft(rsi, 64, 64*64, dist128, rdi, 64*XMM_SCD2, rdx, 64*4*XMM_HSCD2);
		rdi += XMM_SCD2;		/* Next sin/cos data */
		rdx += 4*XMM_HSCD2;	/* Next sin/cos data */
	}
	rsi += -64*64+2*dist128;	/* Next source pointer */
	rdi += -64*XMM_SCD2;	/* Restore sin/cos ptr */
	rdx += -64*4*XMM_HSCD2;	/* Restore sin/cos ptr */

	/* Do 384 four_complex_djbunfft operations */
	/*	distance between fft data elements is 128 */
	/*	do 192 macros each processing 32 data values */

	loops_reset;			/* 3 more iters of 64 */
	r4_x4cl_four_complex_djbunfft_preload;
	for(unsigned int loopA = 3; loopA; loopA--) {
		for(unsigned int loopB = 64; loopB; loopB--) {
			r4_x4cl_four_complex_djbunfft(rsi, 64, 64*64, dist128, rdi, 64*XMM_SCD2);
			rdi += XMM_SCD2;		/* Next sin/cos data */
		}
		rsi += -64*64+2*dist128;	/* Next source pointer */
		rdi += -64*XMM_SCD2;	/* Restore sin/cos ptr */
	}
	rsi = daddr; //rsi += -4*2*dist128;	/* Next source pointer */
	pfing(rcx += 128);		/* Skip pad bytes */

	/* Do inverse FFT level 1,2 */
	/* On input the 64-byte cache lines hold these data values: */
	/*	0	+4K	256	+4K	512	+4K	768	+4K */
	/*	1	... */
	/*	... */
	/*	255	... */
	/*	1024	... */
	/*	... */
	/* On output the 64-byte cache lines hold these data values if pass 1 swizzles: */
	/*	0	+4K	1024	+4K	2048	+4K	3072	+4K */
	/*	1	... */
	/*	... */
	/*	1023 */
	/* On output the 64-byte cache lines hold these data values if pass 2 swizzles: */
	/*	0	+1K	2048	+1K	4096	+1K	6144	+1K */
	/*	1	... */
	/*	... */
	/*	1023 */

	/* Do 512 eight_reals_unfft operations */
	/* Do 512 four_complex_djbunfft operations */
	/*	distance between fft data elements is 512 */
	/*	do 256 macros processing 32 data values */

	rdi = (uintptr_t)g->xsincos_complex;	/* Load sin/cos pointers */
	loops_init_prefetch(256, 64, 2, rcx);
	xtouch(xptr(rcx+4096-64));		/* Preload TLBs */
	xtouch(xptr(rcx+2*4096-64));
	for(unsigned int loopA = 2; loopA; loopA--) {
		for(unsigned int loopB = 128; loopB; loopB--) {
			r4_h4cl_eight_reals_four_complex_djbunfft(rsi, 64, 2*dist128, 4*dist128, rdi, 256*XMM_SCD2, rdx, 256*XMM_HSCD2);
			rdi += XMM_SCD2;		/* Next sin/cos data */
			rdx += XMM_HSCD2;		/* Next sin/cos data */
		}
		rsi += -128*64+dist128;	/* Next source pointer */
	}

	xpass2_11_real_done:;
	end_timer(2);
//}
		if(g->pass2_get_next_block(g) != 0) return;
	}
	r4_pass2_11_levels_complex(g);
}


/* *************** radix4-pass2-11-levels-complex macro ****************** */
/* This macro takes 2 groups of 2048 complex values and */
/* performs the final 11 levels of a standard mostly radix-4 FFT, squares the */
/* results, and does the first 11 levels of the inverse FFT. */


void r4_pass2_11_levels_complex(struct gwasm_data *__restrict g) {
do{
	uintptr_t rbx,rbp,rsi,rdi,daddr;
	vec2f64 xmm0,xmm1,xmm2,xmm3,xmm4,xmm5,xmm6,xmm7;

	/* Type 4 FFTs skip the forward FFT process */

	daddr = rsi = (uintptr_t)g->data_addr;/* Load address of FFT data */
	rbx = g->DIST_TO_FFTSRCARG;
	if likely(!(g->ffttype & 4)) {

		/* Do FFT levels 1,2 */
		/* */
		/* We swizzle the data so that one set of data is in the high word of */
		/* an xmm register, the other set is in the low word.  This means */
		/* a distance of 4096 data values. */
		/* */
		/* On input the 64-byte cache lines hold these data values if pass 2 swizzles: */
		/*	0	+1K	2048	+1K	4096	+1K	6144	+1K */
		/*	1	... */
		/*	... */
		/*	1023 */
		/* On input the 64-byte cache lines hold these data values if pass 1 swizzles: */
		/*	0	+4K	1024	+4K	2048	+4K	3072	+4K */
		/*	1	... */
		/*	... */
		/*	1023	... */
		/* On output the 64-byte cache lines hold these data values: */
		/*	0	+4K	512	+4K	1024	+4K	1536	+4K */
		/*	1	... */
		/*	... */
		/*	255	... */
		/*	1024	... */
		/*	... */

		/* Do 1024 four_complex_fft operations */
		/*	distance between fft data elements is 512 */
		/*	do 512 macros each processing 16 data values */

		start_timer(4);
		rdi = (uintptr_t)g->xsincos_complex;	/* The sin/cos data */
		loops_init(512);			/* 4 iterations of 128 */
		r4_f2cl_four_complex_djbfft_preload;
		for(unsigned int loopA = 4; loopA; loopA--) {
			for(unsigned int loopB = 128; loopB; loopB--) {
				xxprefetch(u64ptr(rdi+XMM_SCD2));	/* Prefetch for next pass */
				r4_f2cl_four_complex_djbfft(rsi, 64, 4*dist128, rdi);
				rdi += XMM_SCD2;		/* Next sin/cos data */
			}
			rsi += -128*64+dist128;	/* Next source pointer */
		}
		rsi = daddr; //rsi += -4*dist128;		/* Restore source pointer */
		end_timer(4);

		/* Do FFT levels 3,4 */
		/* */
		/* On input the 64-byte cache lines hold these data values: */
		/*	0	+4K	512	+4K	1024	+4K	1536	+4K */
		/*	1	... */
		/*	... */
		/*	255	... */
		/*	1024	... */
		/*	... */
		/* On output the 64-byte cache lines hold these data values: */
		/*	0	+4K	128	+4K	256	+4K	384	+4K */
		/*	1	... */
		/*	... */
		/*	127	... */
		/*	512	... */
		/*	... */
		/* */
		/* Do 1024 four_complex_fft operations */
		/*	distance between fft data elements is 128 */
		/*	do 256 macros each processing 32 data values */

		/* Do FFT levels 5,6 */
		/* */
		/* On input the 64-byte cache lines hold these data values: */
		/*	0	+4K	128	+4K	256	+4K	384	+4K */
		/*	1	... */
		/*	... */
		/*	127	... */
		/*	512	... */
		/*	... */
		/* On output the 64-byte cache lines hold these data values: */
		/*	0	+4K	32	+4K	64	+4K	96	+4K */
		/*	1	... */
		/*	... */
		/*	31	... */
		/*	128	... */
		/*	... */
		/* */
		/* Do 1024 four_complex_fft operations */
		/*	distance between fft data elements is 32 */
		/*	do 256 macros each processing 32 data values */

		/* Do FFT levels 7,8 */
		/* */
		/* On input the 64-byte cache lines hold these data values: */
		/*	0	+4K	32	+4K	64	+4K	96	+4K */
		/*	1	... */
		/*	... */
		/*	31	... */
		/*	128	... */
		/*	... */
		/* On output the 64-byte cache lines hold these data values: */
		/*	0	+4K	8	+4K	16	+4K	24	+4K */
		/*	1	... */
		/*	... */
		/*	7	... */
		/*	32	... */
		/*	... */
		/* */
		/* Do 1024 four_complex_fft operations */
		/*	distance between fft data elements is 8 */
		/*	do 256 macros each processing 32 data values */

		/* Do this 64 cache lines at a time to maximize accessing data */
		/* while it is in the L1 cache. */

		start_timer(6);
		pfing(rcx = (uintptr_t)g->data_prefetch);/* Prefetch pointer (we need to prefetch  */
		;					/* 1024 cache lines + 14 pad cache lines) */
		xtouch(xptr(rcx));			/* Load first prefetch TLB */
		rbx = rdi;		/* Copy sin/cos pointers */
		rbp = rdi;
		loops_init_prefetch(256, 64, 2, rcx); /* 2*8 iters of 64 cache lines */
		r4_x4cl_four_complex_djbfft_preload;

		xtouch(xptr(rcx+4096-64));		/* Load some prefetch TLBs */
		xtouch(xptr(rcx+2*4096-64));
		xtouch(xptr(rcx+3*4096-64));
		xtouch(xptr(rcx+4*4096-64));
		xtouch(xptr(rcx+5*4096-64));
		xtouch(xptr(rcx+6*4096-64));

		for(unsigned int loopA = 2; loopA; loopA--) {
			for(unsigned int loopB = 8; loopB; loopB--) {
				for(unsigned int loopC = 16; loopC; loopC--) {
					xxprefetch(u64ptr(rdi+8*XMM_SCD2));	/* Prefetch for next pass */
					r4_x4cl_four_complex_djbfft(rsi, 8*64, dist128, 2*dist128, rdi);
					rdi += 8*XMM_SCD2;		/* Next sine/cosine pointer */
				}
				rsi += -16*8*64;		/* Restore source pointer */
				rdi += -16*8*XMM_SCD2;	/* Restore sine/cosine pointer */
				pfing(rcx += -8*64+8192+128);	/* Prefetch from next 8KB page */
				loops_undo(16);			/* Back up loop counter */

				for(unsigned int loopC = 4; loopC; loopC--) {
					for(unsigned int loopD = 4; loopD; loopD--) {
						r4_x4cl_four_complex_djbfft(rsi, 8*64, 32*64, 64*64, rbx);
						rbx += 8*4*XMM_SCD2;	/* Next sine/cosine pointer */
					}
					rsi += -4*8*64+dist128;	/* Next source pointer */
					rbx += -4*8*4*XMM_SCD2;	/* Restore sine/cosine pointer */
				}
				rsi += -4*dist128;		/* Restore source pointer */
				pfing(rcx += -8*64+8192+128);	/* Prefetch from next 8KB page */
				loops_undo(4*4);			/* Back up loop counter */

				for(unsigned int loopC = 4; loopC; loopC--) {
					for(unsigned int loopD = 4; loopD; loopD--) {
						r4_x4cl_four_complex_djbfft(rsi, 32*64, 8*64, 16*64, rbp);
					}
					rsi += -4*32*64+dist128;	/* Next source pointer */
				}
				pfing(rcx += -2*(8192+128));	/* Back to prefetching from 1st 8KB page */

				rsi += -4*dist128+64;	/* Next source pointer */
				rdi += XMM_SCD2;		/* Next sine/cosine pointer */
				rbx += 4*XMM_SCD2;		/* Next sine/cosine pointer */
				rbp += 16*XMM_SCD2;	/* Next sine/cosine pointer */
			}
			rsi += -8*64+4*dist128;	/* Next source pointer */
			rdi += -8*XMM_SCD2;	/* Restore sine/cosine pointer */
			rbx += -8*4*XMM_SCD2;	/* Restore sine/cosine pointer */
			rbp += -8*16*XMM_SCD2;	/* Restore sine/cosine pointer */
		}
		rsi = daddr; //rsi += -2*4*dist128;	/* Restore source pointer */
		pfing(rcx += -2*8*8*64+3*(8192+128)); /* Move prefetch pointer along */
		end_timer(6);

		/* Do FFT levels 9,10,11 */
		/* */
		/* On input the 64-byte cache lines hold these data values: */
		/*	0	+4K	8	+4K	16	+4K	24	+4K */
		/*	1	... */
		/*	... */
		/*	7	... */
		/*	32	... */
		/*	... */
		/* On output the 64-byte cache lines hold these data values: */
		/*	0	+4K	4	+4K	8	+4K	12	+4K */
		/*	1	... */
		/*	... */
		/*	3	... */
		/*	16	... */
		/*	... */

		/* Do 512 eight_complex_fft operations */
		/*	distance between fft data elements is 1 */
		/*	do 128 macros each processing 64 data values */

		/* Execute the right middle step */

		if likely(g->ffttype & 2) {
			if likely(!(g->ffttype & 1)) { //xpass2_11_levels_complex_2:;
				start_timer(9);
				loops_init_prefetch(128, 64, 1, rcx);
				r8_x8cl_eight_complex_with_square_preload;
				xtouch(xptr(rcx+4096-64));		/* Load prefetch TLBs */
				xtouch(xptr(rcx+2*4096-64));
				for(unsigned int loopA = 8; loopA; loopA--) {
					for(unsigned int loopB = 16; loopB; loopB--) {
						r8_x8cl_eight_complex_with_square(rsi, 8*64, 64, 2*64, 4*64);
					}
					rsi += -16*8*64+dist128;	/* Next source pointer */
				}
				//rsi += -8*dist128;		/* Restore source pointer */
				pfing(rcx += 128);		/* Skip pad bytes */
				end_timer(9);
			}else{ //xpass2_11_levels_complex_3:;
				start_timer(9);
				rbp = g->DIST_TO_MULSRCARG;
				loops_init_prefetch(128, 64, 1, rcx);
				r8_x8cl_eight_complex_with_mult_preload;
				xtouch(xptr(rcx+4096-64));		/* Load prefetch TLBs */
				xtouch(xptr(rcx+2*4096-64));
				for(unsigned int loopA = 8; loopA; loopA--) {
					for(unsigned int loopB = 16; loopB; loopB--) {
						r8_x8cl_eight_complex_with_mult(rsi, 8*64, 64, 2*64, 4*64);
					}
					rsi += -16*8*64+dist128;	/* Next source pointer */
				}
				//rsi += -8*dist128;		/* Restore source pointer */
				pfing(rcx += 128);		/* Skip pad bytes */
				end_timer(9);
			}
		}else{ //xpass2_11_levels_complex_1:;
			start_timer(9);
			loops_init_prefetch(128, 128, 1, rcx);
			r8_x8cl_eight_complex_fft_final_preload;
			for(unsigned int loopA = 8; loopA; loopA--) {
				xtouch(xptr(rcx+4096-64));		/* Load prefetch TLBs */
				for(unsigned int loopB = 16; loopB; loopB--) {
					r8_x8cl_eight_complex_fft_final(rsi, 8*64, 64, 2*64, 4*64);
				}
				rsi += -16*8*64+dist128;	/* Next source pointer */
			}
			end_timer(9);
			goto	xpass2_11_levels_complex_done;
		}
	}else{ //xpass2_11_levels_complex_4:;
		start_timer(9);
		pfing(rcx = (uintptr_t)g->data_prefetch);/* Prefetch pointer */
		xtouch(xptr(rcx));			/* Load first prefetch TLB */
		rbp = g->DIST_TO_MULSRCARG;
		loops_init_prefetch(128, 128, 1, rcx);
		r8_x8cl_eight_complex_with_mulf_preload;
		for(unsigned int loopA = 8; loopA; loopA--) {
			xtouch(xptr(rcx+4096-64));		/* Load prefetch TLBs */
			for(unsigned int loopB = 16; loopB; loopB--) {
				r8_x8cl_eight_complex_with_mulf(rsi, 8*64, 64, 2*64, 4*64);
			}
			rsi += -16*8*64+dist128;	/* Next source pointer */
		}
		//rsi += -8*dist128;		/* Restore source pointer */
		pfing(rcx += 2*128);;		/* Skip pad bytes */
		end_timer(9);
	}
	rsi = daddr;

	/* Do inverse FFT levels 7,8 */
	/* On input the 64-byte cache lines hold these data values: */
	/*	0	+4K	4	+4K	8	+4K	12	+4K */
	/*	1	... */
	/*	... */
	/*	3	... */
	/*	16	... */
	/*	... */
	/* On output the 64-byte cache lines hold these data values: */
	/*	0	+4K	16	+4K	32	+4K	48	+4K */
	/*	1	... */
	/*	... */
	/*	15	... */
	/*	64	... */
	/*	... */
	/* */
	/* Do 1024 four_complex_unfft operations */
	/*	distance between fft data elements is 8 */
	/*	do 256 macros each processing 32 data values */

	/* Do inverse FFT levels 5,6 */
	/* On input the 64-byte cache lines hold these data values: */
	/*	0	+4K	16	+4K	32	+4K	48	+4K */
	/*	1	... */
	/*	... */
	/*	15	... */
	/*	64	... */
	/*	... */
	/* On output the 64-byte cache lines hold these data values: */
	/*	0	+4K	64	+4K	128	+4K	192	+4K */
	/*	1	... */
	/*	... */
	/*	63	... */
	/*	256	... */
	/*	... */
	/* */
	/* Do 1024 four_complex_unfft operations */
	/*	distance between fft data elements is 32 */
	/*	do 256 macros each processing 32 data values */

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
	/* */
	/* Do 1024 four_complex_unfft operations */
	/*	distance between fft data elements is 128 */
	/*	do 256 macros each processing 32 data values */

	/* Do this 64 cache lines at a time to maximize accessing data */
	/* while it is in the L1 cache. */

	start_timer(10);
	rdi = (uintptr_t)g->xsincos_complex;	/* Load sin/cos pointers */
	rdi += 512*XMM_SCD2;
	rbx = rdi;
	rbp = rdi;
	loops_init_prefetch(256, 64, 2, rcx); /* 4*4 iters of 64 cache lines */
	r4_x4cl_four_complex_djbunfft_preload;

	xtouch(xptr(rcx+4096-64));		/* Load some prefetch TLBs */
	xtouch(xptr(rcx+2*4096-64));
	xtouch(xptr(rcx+3*4096-64));
	xtouch(xptr(rcx+4*4096-64));
	xtouch(xptr(rcx+5*4096-64));
	xtouch(xptr(rcx+6*4096-64));

	for(unsigned int loopA = 4; loopA; loopA--) {
		for(unsigned int loopB = 4; loopB; loopB--) {
			for(unsigned int loopC = 2; loopC; loopC--) {
				for(unsigned int loopD = 8; loopD; loopD--) {
					r4_x4cl_four_complex_djbunfft(rsi, 16*64, 4*64, 8*64, rbp, 4*16*XMM_SCD2);
				}
				rsi += -8*16*64+dist128;	/* Next source pointer */
			}
			rsi += -2*dist128;		/* Restore source pointer */
			pfing(rcx += -8*64+8192+128);	/* Prefetch from next 8KB page */
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
			pfing(rcx += -8*64+8192+128);	/* Prefetch from next 8KB page */
			loops_undo(2*2*4);		/* Restore loop counter */

			for(unsigned int loopC = 16; loopC; loopC--) {
				r4_x4cl_four_complex_djbunfft(rsi, 4*64, 64*64, dist128, rdi, 64*1*XMM_SCD2);
				rdi += 4*1*XMM_SCD2;	/* Next sine/cosine pointer */
			}
			rdi += -16*4*1*XMM_SCD2;	/* Restore sine/cosine pointer */
			pfing(rcx += -2*(8192+128));	/* Back to prefetching from 1st 8KB page */

			rsi += -16*4*64+64;	/* Next source pointer */
			rdi += 1*XMM_SCD2;		/* Next sine/cosine pointer */
			rbx += 4*XMM_SCD2;		/* Next sine/cosine pointer */
			rbp += 16*XMM_SCD2;	/* Next sine/cosine pointer */
		}
		rsi += -4*64+2*dist128;	/* Next source pointer */
		rdi += -4*1*XMM_SCD2;	/* Restore sine/cosine pointer */
		rbx += -4*4*XMM_SCD2;	/* Restore sine/cosine pointer */
		rbp += -4*16*XMM_SCD2;	/* Restore sine/cosine pointer */
	}
	rsi = daddr; //rsi += -4*2*dist128;	/* Restore source pointer */
	pfing(rcx += -4*4*8*64+3*(8192+128)); /* Move prefetch pointer along */
	end_timer(10);

	/* Do inverse FFT levels 1,2 */
	/* On input the 64-byte cache lines hold these data values: */
	/*	0	+4K	256	+4K	512	+4K	768	+4K */
	/*	1	... */
	/*	... */
	/*	255	... */
	/*	1024	... */
	/*	... */
	/* On output the 64-byte cache lines hold these data values if pass 1 swizzles: */
	/*	0	+4K	1024	+4K	2048	+4K	3072	+4K */
	/*	1	... */
	/*	... */
	/*	1023 */
	/* On output the 64-byte cache lines hold these data values if pass 2 swizzles: */
	/*	0	+1K	2048	+1K	4096	+1K	6144	+1K */
	/*	1	... */
	/*	... */
	/*	1023 */

	/* Do 1024 four_complex_unfft operations */
	/*	distance between fft data elements is 512 */
	/*	do 256 macros each processing 24 data values */

	start_timer(14);
	rdi = (uintptr_t)g->xsincos_complex;	/* Load group multiplier pointer */
	loops_init_prefetch(256, 64, 2, rcx); /* 2 iters of 128 */
	xtouch(xptr(rcx+4096-64));		/* Load prefetch TLBs */
	xtouch(xptr(rcx+2*4096-64));
	r4_x4cl_four_complex_djbunfft_preload;
	for(unsigned int loopA = 2; loopA; loopA--) {
		for(unsigned int loopB = 128; loopB; loopB--) {
			r4_x4cl_four_complex_djbunfft(rsi, 64, 2*dist128, 4*dist128, rdi, 256*XMM_SCD2);
			rdi += XMM_SCD2;		/* Next set of sin/cos data */
		}
		rsi += -128*64+dist128;	/* Next source pointer */
	}
	end_timer(14);

	xpass2_11_levels_complex_done:;
}while(g->pass2_get_next_block(g) == 0);
}
