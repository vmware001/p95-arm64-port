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

/* *************** radix4-pass2-3072-levels-real macro ****************** */
/* This macro takes 6144 real values and performs the final 12.585 levels of */
/* the FFT.  Independently, this macro also takes 3072 complex values and */
/* performs the final 11.585 levels of the FFT.  The results are squared, and */
/* the first 11.585/12.585 levels of the inverse FFT are done. */

void r4_pass2_3072_levels_complex(struct gwasm_data *__restrict g);
void xpass2_r4_3072_levels(struct gwasm_data *__restrict g) {
	g->thread_work_routine = r4_pass2_3072_levels_complex;
	pass2_wake_up_threads(g);
	if(g->NEGACYCLIC_FFT == 0) {
//void r4_pass2_3072_levels_real(struct gwasm_data *__restrict g) {
	uintptr_t rdx,rbx,rbp,rsi,rdi,daddr;
	vec2f64 xmm0,xmm1,xmm2,xmm3,xmm4,xmm5,xmm6,xmm7;

	/* We switch to the same format used in one-pass FFTs (that is, the */
	/* high word of xmm registers being 1 greater than the low word): */

	/* Type 4 FFTs skip the forward FFT process */

	start_timer(2);
	daddr = (uintptr_t)g->data_addr;/* Load source address */
	rsi = daddr;
	rbx = g->DIST_TO_FFTSRCARG;
	if likely(!(g->ffttype & 4)) { //xpass2_3072_levels_real_4;

		/* Do FFT level 1,2 (actually 1.585 levels) */
		/* */
		/* On input the 64-byte cache lines hold these data values if pass 2 swizzles: */
		/*	0	+1.5K	3072	+1.5K	6144	+1.5K	9216	+1.5K */
		/*	1	... */
		/*	... */
		/*	1535 */
		/* On input the 64-byte cache lines hold these data values if pass 1 swizzles: */
		/*	0	+6K	1536	+6K	3072	+6K	4608	+6K */
		/*	1	... */
		/*	... */
		/*	1535	... */
		/* On output the 64-byte cache lines hold these data values: */
		/*	0	+6K	512	+6K	1024	+6K	1536	+6K */
		/*	1	... */
		/*	... */
		/*	511	... */
		/*	2048	... */
		/*	... */

		/* Do 1024 six_reals_fft on values 0 - 6143 */
		/* Do 1024 three_complex_djbfft on values 6144 - 12287 */
		/*	distance between fft data elements is 1K */
		/*	do 512 iterations */

		rdi = (uintptr_t)g->xsincos_complex;	/* Ptr to the sin/cos data */
		rdx = (uintptr_t)g->sincos3;/* Ptr to real sin/cos data */
		loops_init(512);
		for(unsigned int loopA = 4; loopA; loopA--) {
			for(unsigned int loopB = 128; loopB; loopB--) {
				r3_fh3cl_six_reals_three_complex_djbfft(rsi, 64, 4*dist128, rdi, 512*XMM_SCD1, rdx, 512*XMM_HSCD1);
				rdi += XMM_SCD1;		/* Next sin/cos data */
				rdx += XMM_HSCD1;		/* Next sin/cos data */
			}
			rsi += -128*64+dist128;	/* Next source pointer */
		}
		rsi = daddr; //rsi += -4*dist128;		/* Restore source pointer */
		rdi += -512*XMM_SCD1;	/* Restore sin/cos ptr */
		rdx += 512*XMM_HSCD1;	/* Next sin/cos ptr */

		/* Do FFT levels 3,4 */
		/* Values 0-2047 is real data, 2048-12287 is complex data. */
		/* */
		/* On input the 64-byte cache lines hold these data values: */
		/*	0	+6K	512	+6K	1024	+6K	1536	+6K */
		/*	1	... */
		/*	... */
		/*	511	... */
		/*	2048	... */
		/*	... */
		/* On output the 64-byte cache lines hold these data values: */
		/*	0	+6K	256	+6K	512	+6K	768	+6K */
		/*	1	... */
		/*	... */
		/*	255	... */
		/*	1024	... */
		/*	... */

		/* Do 256 eight_reals_fft operations */
		/* Do 256 four_complex_djbfft operations */
		/*	distance between fft data elements is 256 */
		/*	do 256 macros each processing 16 data values */

		pfing(rcx = (uintptr_t)g->data_prefetch); /* Prefetch pointer */
		xtouch(xptr(rcx));			/* Preload first TLB */
		rbx = rdi;
		loops_init(768);
		for(unsigned int loopA = 2; loopA; loopA--) {
			for(unsigned int loopB = 128; loopB; loopB--) {
				r4_h2cl_2sc_eight_reals_four_complex_djbfft(rsi, 64, 2*dist128, rdi, rbx, rdx);
				rdi += 3*XMM_SCD1;		/* Next sin/cos data */
				rbx += 6*XMM_SCD1;		/* Next sin/cos data */
				rdx += XMM_HSCD2;		/* Next sin/cos data */
			}
			rsi += -128*64+dist128;	/* Next source pointer */
		}
		rsi += -2*dist128+4*dist128; /* Next source pointer */
		rdi += -256*3*XMM_SCD1;	/* Restore sin/cos ptr */
		rbx += -256*6*XMM_SCD1;	/* Restore sin/cos ptr */
		rdx += -256*XMM_HSCD2;	/* Restore sin/cos ptr */

		/* Do 1024 four_complex_fft operations */
		/*	distance between fft data elements is 256 */
		/*	do 512 macros each processing 16 data values */

		loops_reset;			/* 2 more iters of 256 */
		r4_x2cl_2sc_four_complex_djbfft_preload;
		for(unsigned int loopA = 2; loopA; loopA--) {
			for(unsigned int loopB = 2; loopB; loopB--) {
				for(unsigned int loopC = 128; loopC; loopC--) {
					r4_x2cl_2sc_four_complex_djbfft(rsi, 64, 2*dist128, rdi, rbx);
					rdi += 3*XMM_SCD1;		/* Next sin/cos data */
					rbx += 6*XMM_SCD1;		/* Next sin/cos data */
				}
				rsi += -128*64+dist128;	/* Next source pointer */
			}
			rsi += -2*dist128+4*dist128; /* Next source pointer */
			rdi += -256*3*XMM_SCD1;	/* Restore sin/cos ptr */
			rbx += -256*6*XMM_SCD1;	/* Restore sin/cos ptr */
		}
		rsi = daddr; //rsi += -3*4*dist128;	/* Restore source pointer */

		/* Do FFT levels 5,6 */
		/* Values 0-511 is real data, 512-12287 is complex data. */
		/* */
		/* On input the 64-byte cache lines hold these data values: */
		/*	0	+6K	256	+6K	512	+6K	768	+6K */
		/*	1	... */
		/*	... */
		/*	255	... */
		/*	1024	... */
		/*	... */
		/* On output the 64-byte cache lines hold these data values: */
		/*	0	+6K	64	+6K	128	+6K	192	+6K */
		/*	1	... */
		/*	... */
		/*	63	... */
		/*	256	... */
		/*	... */

		/* Do 64 eight_reals_fft operations */
		/* Do 192 four_complex_djbfft operations */
		/*	distance between fft data elements is 64 */
		/*	do 64 macros each processing 32 data values */

		rdi += 1536*XMM_SCD1;	/* Load sin/cos pointer */
		loops_init_prefetch(384, 64, 2, rcx);
		xtouch(xptr(rcx+4096-64));		/* Preload the TLBs */
		xtouch(xptr(rcx+2*4096-64));
		xtouch(xptr(rcx+3*4096-64));
		r4_h4cl_eight_reals_four_complex_djbfft_preload;
		for(unsigned int loopA = 64; loopA; loopA--) {
			r4_h4cl_eight_reals_four_complex_djbfft(rsi, 64, 64*64, dist128, rdi, rdx);
			rdi += XMM_SCD2;		/* Next sin/cos data */
			rdx += 4*XMM_HSCD2;	/* Next sin/cos data */
		}
		rsi += -64*64+2*dist128;	/* Next source pointer */
		rdi += -64*XMM_SCD2;	/* Restore sin/cos ptr */
		rdx += -64*4*XMM_HSCD2;	/* Restore sin/cos ptr */

		/* Do 1280 four_complex_djbfft operations */
		/*	distance between fft data elements is 64 */
		/*	do 320 macros each processing 32 data values */

		loops_reset;			/* 5 iterations of 64 */
		r4_x4cl_four_complex_djbfft_preload;
		for(unsigned int loopA = 5; loopA; loopA--) {
			for(unsigned int loopB = 64; loopB; loopB--) {
				r4_x4cl_four_complex_djbfft(rsi, 64, 64*64, dist128, rdi);
				rdi += XMM_SCD2;		/* Next sin/cos data */
			}
			rsi += -64*64+2*dist128;	/* Next source pointer */
			rdi += -64*XMM_SCD2;	/* Restore sin/cos ptr */
		}
		rsi = daddr; //rsi += -6*2*dist128;	/* Next source pointer */

		/* Do FFT levels 7,8 */
		/* Values 0-127 is real data, 128-12687 is complex data. */
		/* */
		/* On input the 64-byte cache lines hold these data values: */
		/*	0	+6K	64	+6K	128	+6K	192	+6K */
		/*	1	... */
		/*	... */
		/*	63	... */
		/*	256	... */
		/*	... */
		/* On output the 64-byte cache lines hold these data values: */
		/*	0	+6K	16	+6K	32	+6K	48	+6K */
		/*	1	... */
		/*	... */
		/*	15	... */
		/*	64	... */
		/*	... */

		/* Do 16 eight_reals_fft operations */
		/* Do 48 four_complex_djbfft operations */
		/*	distance between fft data elements is 16 */
		/*	do 16 macros each processing 32 data values */

		loops_init_prefetch(384, 64, 2, rcx);
		xtouch(xptr(rcx+4096-64));		/* Preload the TLBs */
		xtouch(xptr(rcx+2*4096-64));
		xtouch(xptr(rcx+3*4096-64));
		r4_h4cl_eight_reals_four_complex_djbfft_preload;
		for(unsigned int loopA = 16; loopA; loopA--) {
			r4_h4cl_eight_reals_four_complex_djbfft(rsi, 64, 16*64, 32*64, rdi, rdx);
			rdi += 4*XMM_SCD2;		/* Next sin/cos data */
			rdx += 16*XMM_HSCD2;	/* Next sin/cos data */
		}
		rsi += -16*64+64*64;	/* Next source pointer */
		rdi += -16*4*XMM_SCD2;	/* Restore sin/cos ptr */
		rdx += -16*16*XMM_HSCD2;	/* Restore sin/cos ptr */

		/* Do 1472 four_complex_djbfft operations */
		/*	distance between fft data elements is 16 */
		/*	do 368 macros each processing 32 data values */

		loops_reset;			/* 12 iters of 1 or 2 of 16 */
		r4_x4cl_four_complex_djbfft_preload;
		for(unsigned int loopA = 12*2-1; loopA; ) {
			do{
				for(unsigned int loopC = 16; loopC; loopC--) {
					r4_x4cl_four_complex_djbfft(rsi, 64, 16*64, 32*64, rdi);
					rdi += 4*XMM_SCD2;		/* Next sin/cos data */
				}
				rsi += -16*64+64*64;	/* Next source pointer */
				rdi += -16*4*XMM_SCD2;	/* Restore sin/cos ptr */
			}while(--loopA & 1);
			rsi += -2*64*64+dist128;	/* Next source pointer */
		}
		rsi = daddr; //rsi += -12*dist128;	/* Next source pointer */

		/* Do FFT levels 9,10 */
		/* Values 0-31 is real data, 32-12687 is complex data. */
		/* */
		/* On input the 64-byte cache lines hold these data values: */
		/*	0	+6K	16	+6K	32	+6K	48	+6K */
		/*	1	... */
		/*	... */
		/*	15	... */
		/*	64	... */
		/*	... */
		/* On output the 64-byte cache lines hold these data values: */
		/*	0	+6K	4	+6K	8	+6K	12	+6K */
		/*	1	... */
		/*	... */
		/*	3	... */
		/*	16	... */
		/*	... */

		/* Do 4 eight_reals_fft operations */
		/* Do 12 four_complex_djbfft operations */
		/*	distance between fft data elements is 4 */
		/*	do 4 macros each processing 32 data values */

		loops_init_prefetch(384, 64, 2, rcx);
		xtouch(xptr(rcx+4096-64));		/* Preload the TLBs */
		xtouch(xptr(rcx+2*4096-64));
		xtouch(xptr(rcx+3*4096-64));
		r4_h4cl_eight_reals_four_complex_djbfft_preload;
		for(unsigned int loopA = 4; loopA; loopA--) {
			r4_h4cl_eight_reals_four_complex_djbfft(rsi, 64, 4*64, 8*64, rdi, rdx);
			rdi += 16*XMM_SCD2;	/* Next sin/cos data */
			rdx += 64*XMM_HSCD2;	/* Next sin/cos data */
		}
		rsi += -4*64+16*64;	/* Next source pointer */
		rdi += -4*16*XMM_SCD2;	/* Restore sin/cos ptr */
		rdx += -4*64*XMM_HSCD2;	/* Restore sin/cos ptr */

		/* Do 1520 four_complex_djbfft operations */
		/*	distance between fft data elements is 4 */
		/*	do 380 macros each processing 32 data values */

		loops_reset;			/* 12 iters of 7 or 8 of 4 */
		r4_x4cl_four_complex_djbfft_preload;
		for(unsigned int loopA = 12*8-1; loopA; ) {
			do{
				for(unsigned int loopC = 4; loopC; loopC--) {
					r4_x4cl_four_complex_djbfft(rsi, 64, 4*64, 8*64, rdi);
					rdi += 16*XMM_SCD2;	/* Next sin/cos data */
				}
				rsi += -4*64+16*64;	/* Next source pointer */
				rdi += -4*16*XMM_SCD2;	/* Restore sin/cos ptr */
			}while(--loopA & 7);
			rsi += -8*16*64+dist128;	/* Next source pointer */
		}
		rsi = daddr; //rsi += -12*dist128;	/* Next source pointer */

		/* Do FFT levels 11,12 */
		/* Values 0-7 is real data, 8-12687 is complex data. */
		/* */
		/* On input the 64-byte cache lines hold these data values: */
		/*	0	+6K	4	+6K	8	+6K	12	+6K */
		/*	1	... */
		/*	... */
		/*	3	... */
		/*	16	... */
		/*	... */
		/* On output the 64-byte cache lines hold these data values: */
		/*	0	+6K	2	+6K	4	+6K	6	+6K */
		/*	1	... */
		/*	8	... */
		/*	... */

		/* Do 1 eight_reals_fft operations */
		/* Do 3 four_complex_fft operations */
		/*	distance between fft data elements is 1 */

		/* Do 1532 four_complex_fft operations */
		/*	distance between fft data elements is 1 */
		/*	do 383 macros each processing 32 data values */

		/* Execute the proper middle step */

		if likely(g->ffttype & 2) {
			if likely(!(g->ffttype & 1)) { // xpass2_3072_levels_real_2:;
				r4_h4cl_eight_reals_four_complex_with_square(rsi, 4*64, 64, 2*64);
				loops_init_prefetch(384, 64, 4, rcx, 1, 32);
				xtouch(xptr(rcx+4096-64));		/* Preload the TLBs */
				xtouch(xptr(rcx+2*4096-64));
				r4_x4cl_four_complex_with_square_preload;
				for(unsigned int loopA = 12*32-1; loopA; ) {
					do{
						r4_x4cl_four_complex_with_square(rsi, 4*64, 64, 2*64);
					}while(--loopA & 31);
					rsi += -32*4*64+dist128;	/* Next source pointer */
				}
				//rsi += -12*dist128;	/* Next source pointer */
			}else{ // xpass2_3072_levels_real_3:;
				rbp = g->DIST_TO_MULSRCARG;
				r4_h4cl_eight_reals_four_complex_with_mult(rsi, 4*64, 64, 2*64);
				loops_init_prefetch(384, 64, 4, rcx, 1, 32);
				xtouch(xptr(rcx+4096-64));		/* Preload the TLBs */
				xtouch(xptr(rcx+2*4096-64));
				for(unsigned int loopA = 12*32-1; loopA; ) {
					do{
						r4_x4cl_four_complex_with_mult(rsi, 4*64, 64, 2*64);
					}while(--loopA & 31);
					rsi += -32*4*64+dist128;	/* Next source pointer */
				}
				//rsi += -12*dist128;	/* Next source pointer */
			}
		}else{ // xpass2_3072_levels_real_1:;
			r4_h4cl_eight_reals_four_complex_fft_final(rsi, 4*64, 64, 2*64);
			loops_init_prefetch(384, 128, 1, rcx, 1, 32);
			for(unsigned int loopA = 12*32-1; loopA; ) {
				xtouch(xptr(rcx+4096-64));		/* Preload the TLBs */
				do{
					r4_x4cl_four_complex_fft_final(rsi, 4*64, 64, 2*64);
				}while(--loopA & 31);
				rsi += -32*4*64+dist128;	/* Next source pointer */
			}
			goto	xpass2_3072_real_done;
		}
	}else{ // xpass2_3072_levels_real_4:;
		rbp = g->DIST_TO_MULSRCARG;
		pfing(rcx = (uintptr_t)g->data_prefetch); /* Load prefetch pointer */
		xtouch(xptr(rcx));			/* Preload first TLB */
		r4_h4cl_eight_reals_four_complex_with_mulf(rsi, 4*64, 64, 2*64);
		loops_init_prefetch(384, 64, 1, rcx, 1, 32);
		for(unsigned int loopA = 12*32-1; loopA; ) {
			xtouch(xptr(rcx+4096-64));
			do{
				r4_x4cl_four_complex_with_mulf(rsi, 4*64, 64, 2*64);
			}while(--loopA & 31);
			rsi += -32*4*64+dist128;	/* Next source pointer */
		}
		//rsi += -12*dist128;	/* Next source pointer */
		pfing(rcx += -23*64);		/* Back up prefetch pointer */
	}
	rsi = daddr;

	/* Do inverse FFT levels 9,10 */
	/* On input the 64-byte cache lines hold these data values: */
	/*	0	+6K	2	+6K	4	+6K	6	+6K */
	/*	1	... */
	/*	8	... */
	/*	... */
	/* On output the 64-byte cache lines hold these data values: */
	/*	0	+6K	8	+6K	16	+6K	24	+6K */
	/*	1	... */
	/*	... */
	/*	7	... */
	/*	32	... */
	/*	... */

	/* Do 4 eight_reals_unfft operations */
	/* Do 4 four_complex_djbunfft operations */
	/*	distance between fft data elements is 4 */
	/*	do 2 macros processing 32 data values */

	rdi = (uintptr_t)g->xsincos_complex;	/* Load sin/cos pointers */
	rdi += 1536*XMM_SCD1;
	rdx = (uintptr_t)g->sincos3;
	rdx += 1024*XMM_HSCD1;
	loops_init_prefetch(384, 64, 2, rcx);
	xtouch(xptr(rcx+4096-64));		/* Preload TLBs */
	xtouch(xptr(rcx+2*4096-64));
	xtouch(xptr(rcx+3*4096-64));
	for(unsigned int loopA = 2; loopA; loopA--) {
		r4_h4cl_eight_reals_four_complex_djbunfft(rsi, 64, 2*64, 4*64, rdi, 2*16*XMM_SCD2, rdx, 2*64*XMM_HSCD2);
		rdi += 16*XMM_SCD2;	/* Next sin/cos data */
		rdx += 64*XMM_HSCD2;	/* Next sin/cos data */
	}
	rsi += -2*64+8*64;		/* Next source pointer */
	rdi += -2*16*XMM_SCD2;	/* Restore sin/cos ptr */
	rdx += -2*64*XMM_HSCD2;	/* Restore sin/cos ptr */

	/* Do 1528 four_complex_djbunfft operations */
	/*	distance between fft data elements is 4 */
	/*	do 382 macros each processing 32 data values */

	loops_reset;			/* 12 iters of 15 or 16 of 2 */
	r4_x4cl_four_complex_djbunfft_preload;
	for(unsigned int loopA = 12*16-1; loopA; ) {
		do{
			for(unsigned int loopC = 2; loopC; loopC--) {
				r4_x4cl_four_complex_djbunfft(rsi, 64, 2*64, 4*64, rdi, 2*16*XMM_SCD2);
				rdi += 16*XMM_SCD2;	/* Next sine/cosine pointer */
			}
			rsi += -2*64+8*64;		/* Next source pointer */
			rdi += -2*16*XMM_SCD2;	/* Restore sin/cos ptr */
		}while(--loopA & 15);
		rsi += -16*8*64+dist128;	/* Next source pointer */
	}
	rsi = daddr; //rsi += -12*dist128;	/* Next source pointer */

	/* Do inverse FFT levels 7,8 */
	/* On input the 64-byte cache lines hold these data values: */
	/*	0	+6K	8	+6K	16	+6K	24	+6K */
	/*	1	... */
	/*	... */
	/*	7	... */
	/*	32	... */
	/*	... */
	/* On output the 64-byte cache lines hold these data values: */
	/*	0	+6K	32	+6K	64	+6K	96	+6K */
	/*	1	... */
	/*	... */
	/*	31	... */
	/*	128	... */
	/*	... */

	/* Do 16 eight_reals_unfft operations */
	/* Do 16 four_complex_djbunfft operations */
	/*	distance between fft data elements is 16 */
	/*	do 8 macros processing 32 data values */

	loops_init_prefetch(384, 64, 2, rcx);
	xtouch(xptr(rcx+4096-64));		/* Preload TLBs */
	xtouch(xptr(rcx+2*4096-64));
	xtouch(xptr(rcx+3*4096-64));
	for(unsigned int loopA = 8; loopA; loopA--) {
		r4_h4cl_eight_reals_four_complex_djbunfft(rsi, 64, 8*64, 16*64, rdi, 8*4*XMM_SCD2, rdx, 8*16*XMM_HSCD2);
		rdi += 4*XMM_SCD2;		/* Next sin/cos data */
		rdx += 16*XMM_HSCD2;	/* Next sin/cos data */
	}
	rsi += -8*64+32*64;	/* Next source pointer */
	rdi += -8*4*XMM_SCD2;	/* Restore sin/cos ptr */
	rdx += -8*16*XMM_HSCD2;	/* Restore sin/cos ptr */

	/* Do 1504 four_complex_djbunfft operations */
	/*	distance between fft data elements is 16 */
	/*	do 376 macros each processing 32 data values */

	loops_reset;			/* 12 iters of 3 or 4 of 4 */
	r4_x4cl_four_complex_djbunfft_preload;
	for(unsigned int loopA = 12*4-1; loopA; ) {
		do{
			for(unsigned int loopC = 8; loopC; loopC--) {
				r4_x4cl_four_complex_djbunfft(rsi, 64, 8*64, 16*64, rdi, 8*4*XMM_SCD2);
				rdi += 4*XMM_SCD2;		/* Next sin/cos data */
			}
			rsi += -8*64+32*64;	/* Next source pointer */
			rdi += -8*4*XMM_SCD2;	/* Restore sin/cos ptr */
		}while(--loopA & 3);
		rsi += -4*32*64+dist128;	/* Next source pointer */
	}
	rsi = daddr; //rsi += -12*dist128;	/* Next source pointer */

	/* Do inverse FFT levels 5,6 */
	/* On input the 64-byte cache lines hold these data values: */
	/*	0	+6K	32	+6K	64	+6K	96	+6K */
	/*	1	... */
	/*	... */
	/*	31	... */
	/*	128	... */
	/*	... */
	/* On output the 64-byte cache lines hold these data values: */
	/*	0	+6K	128	+6K	256	+6K	384	+6K */
	/*	1	... */
	/*	... */
	/*	127	... */
	/*	512	... */
	/*	... */

	/* Do 64 eight_reals_unfft operations */
	/* Do 64 four_complex_djbunfft operations */
	/*	distance between fft data elements is 64 */
	/*	do 32 macros processing 32 data values */

	loops_init_prefetch(384, 64, 2, rcx);
	xtouch(xptr(rcx+4096-64));		/* Preload TLBs */
	xtouch(xptr(rcx+2*4096-64));
	xtouch(xptr(rcx+3*4096-64));
	for(unsigned int loopA = 32; loopA; loopA--) {
		r4_h4cl_eight_reals_four_complex_djbunfft(rsi, 64, 32*64, 64*64, rdi, 32*XMM_SCD2, rdx, 32*4*XMM_HSCD2);
		rdi += XMM_SCD2;		/* Next sin/cos data */
		rdx += 4*XMM_HSCD2;	/* Next sin/cos data */
	}
	rsi += -32*64+dist128;	/* Next source pointer */
	rdi += -32*XMM_SCD2;	/* Restore sin/cos ptr */
	rdx += -32*4*XMM_HSCD2;	/* Restore sin/cos ptr */

	/* Do 1408 four_complex_djbunfft operations */
	/*	distance between fft data elements is 64 */
	/*	do 352 macros each processing 32 data values */

	loops_reset;			/* 11 iters of 32 */
	r4_x4cl_four_complex_djbunfft_preload;
	for(unsigned int loopA = 11; loopA; loopA--) {
		for(unsigned int loopB = 32; loopB; loopB--) {
			r4_x4cl_four_complex_djbunfft(rsi, 64, 32*64, 64*64, rdi, 32*XMM_SCD2);
			rdi += XMM_SCD2;		/* Next sin/cos data */
		}
		rsi += -32*64+dist128;	/* Next source pointer */
		rdi += -32*XMM_SCD2;	/* Restore sin/cos ptr */
	}
	rsi = daddr; //rsi += -12*dist128;	/* Next source pointer */

	/* Do inverse FFT levels 3,4 */
	/* On input the 64-byte cache lines hold these data values: */
	/*	0	+6K	128	+6K	256	+6K	384	+6K */
	/*	1	... */
	/*	... */
	/*	127	... */
	/*	512	... */
	/*	... */
	/* On output the 64-byte cache lines hold these data values: */
	/*	0	+6K	512	+6K	1024	+6K	1536	+6K */
	/*	1	... */
	/*	... */
	/*	511	... */
	/*	2048	... */
	/*	... */

	/* Do 256 eight_reals_unfft operations */
	/* Do 256 four_complex_djbunfft operations */
	/*	distance between fft data elements is 256 */
	/*	do 128 macros processing 32 data values */

	loops_init_prefetch(384, 64, 2, rcx);
	xtouch(xptr(rcx+4096-64));		/* Preload TLBs */
	xtouch(xptr(rcx+2*4096-64));
	xtouch(xptr(rcx+3*4096-64));
	rdi = (uintptr_t)g->xsincos_complex;	/* Load sin/cos pointer */
	rbx = rdi;
	for(unsigned int loopA = 128; loopA; loopA--) {
		r4_h4cl_2sc_eight_reals_four_complex_djbunfft(rsi, 64, dist128, 2*dist128, rdi, 128*3*XMM_SCD1, rbx, 128*6*XMM_SCD1, rdx, 128*XMM_HSCD2);
		rdi += 3*XMM_SCD1;		/* Next sin/cos data */
		rbx += 6*XMM_SCD1;		/* Next sin/cos data */
		rdx += XMM_HSCD2;		/* Next sin/cos data */
	}
	rsi += -128*64+4*dist128;	/* Next source pointer */
	rdi += -128*3*XMM_SCD1;	/* Restore sin/cos ptr */
	rbx += -128*6*XMM_SCD1;	/* Restore sin/cos ptr */
	rdx += -128*XMM_HSCD2-1024*XMM_HSCD1; /* Restore sin/cos ptr */

	/* Do 1024 four_complex_djbunfft operations */
	/*	distance between fft data elements is 256 */
	/*	do 256 macros each processing 32 data values */

	loops_reset;			/* 2 iters of 128 */
	r4_x4cl_2sc_four_complex_djbunfft_preload;
	for(unsigned int loopA = 2; loopA; loopA--) {
		for(unsigned int loopB = 128; loopB; loopB--) {
			r4_x4cl_2sc_four_complex_djbunfft(rsi, 64, dist128, 2*dist128, rdi, 128*3*XMM_SCD1, rbx, 128*6*XMM_SCD1);
			rdi += 3*XMM_SCD1;		/* Next sin/cos data */
			rbx += 6*XMM_SCD1;		/* Next sin/cos data */
		}
		rsi += -128*64+4*dist128;	/* Next source pointer */
		rdi += -128*3*XMM_SCD1;	/* Restore sin/cos ptr */
		rbx += -128*6*XMM_SCD1;	/* Restore sin/cos ptr */
	}
	rsi = daddr; //rsi += -3*4*dist128;	/* Next source pointer */

	/* Do inverse FFT level 1,2 (actually 1.585 FFT levels) */
	/* On input the 64-byte cache lines hold these data values: */
	/*	0	+6K	512	+6K	1024	+6K	1536	+6K */
	/*	1	... */
	/*	... */
	/*	511	... */
	/*	2048	... */
	/*	... */
	/* On output the 64-byte cache lines hold these data values if pass 1 swizzles: */
	/*	0	+6K	1536	+6K	3072	+6K	4608	+6K */
	/*	1	... */
	/*	... */
	/*	1535 */
	/* On output the 64-byte cache lines hold these data values if pass 2 swizzles: */
	/*	0	+1.5K	3072	+1.5K	6144	+1.5K	9216	+1.5K */
	/*	1	... */
	/*	... */
	/*	1535 */

	/* Do 1024 six_reals_unfft operations */
	/* Do 1024 three_complex_djbunfft operations */
	/*	distance between fft data elements is 1024 */
	/*	do 512 macros processing 24 data values */

	pfing(rcx += -10*64);		/* Back up prefetch pointer by 10 cache lines */
	loops_init_prefetch(512, 64, 4, rcx);
	xtouch(xptr(rcx+4096-64));		/* Preload TLBs */
	xtouch(xptr(rcx+2*4096-64));
	for(unsigned int loopA = 4; loopA; loopA--) {
		for(unsigned int loopB = 128; loopB; loopB--) {
			r3_h3cl_six_reals_three_complex_djbunfft(rsi, 64, 4*dist128, rdi, 512*XMM_SCD1, rdx, 512*XMM_HSCD1);
			rdi += XMM_SCD1;		/* Next sin/cos data */
			rdx += XMM_HSCD1;		/* Next sin/cos data */
		}
		rsi += -128*64+dist128;	/* Next source pointer */
	}

	xpass2_3072_real_done:;
	end_timer(2);
		if(g->pass2_get_next_block(g) != 0) return;
	}
	r4_pass2_3072_levels_complex(g);
}

/* *************** radix4-pass2-3072-levels-complex macro ****************** */
/* This macro takes 2 groups of 3072 complex values and */
/* performs the final 11.585 levels of a standard mostly radix-4 FFT, squares the */
/* results, and does the first 11.585 levels of the inverse FFT. */


void r4_pass2_3072_levels_complex(struct gwasm_data *__restrict g) {
do{
	uintptr_t rdx,rbx,rbp,rsi,rdi,daddr;
	vec2f64 xmm0,xmm1,xmm2,xmm3,xmm4,xmm5,xmm6,xmm7;
	/* Type 4 FFTs skip the forward FFT process */

	daddr = (uintptr_t)g->data_addr;/* Load address of FFT data */
	rsi = daddr;
	rbx = g->DIST_TO_FFTSRCARG;
	if likely(!(g->ffttype & 4)) { //xpass2_3072_levels_complex_4;

		/* Do FFT levels 1,2 (actually 1.585 levels) */
		/* */
		/* We swizzle the data so that one set of data is in the high word of */
		/* an xmm register, the other set is in the low word.  This means */
		/* a distance of 6144 data values. */
		/* */
		/* On input the 64-byte cache lines hold these data values if pass 2 swizzles: */
		/*	0	+1.5K	3072	+1.5K	6144	+1.5K	9216	+1.5K */
		/*	1	... */
		/*	... */
		/*	1535 */
		/* On input the 64-byte cache lines hold these data values if pass 1 swizzles: */
		/*	0	+6K	1536	+6K	3072	+6K	4608	+6K */
		/*	1	... */
		/*	... */
		/*	1535	... */
		/* On output the 64-byte cache lines hold these data values: */
		/*	0	+6K	512	+6K	1024	+6K	1536	+6K */
		/*	1	... */
		/*	... */
		/*	511	... */
		/*	2048	... */
		/*	... */

		/* Do 2048 three_complex_fft operations */
		/*	distance between fft data elements is 1024 */
		/*	do 512 macros each processing 24 data values */

		/* Do FFT levels 3,4 */
		/* */
		/* On input the 64-byte cache lines hold these data values: */
		/*	0	+6K	512	+6K	1024	+6K	1536	+6K */
		/*	1	... */
		/*	... */
		/*	511	... */
		/*	2048	... */
		/*	... */
		/* On output the 64-byte cache lines hold these data values: */
		/*	0	+6K	256	+6K	512	+6K	768	+6K */
		/*	1	... */
		/*	... */
		/*	255	... */
		/*	1024	... */
		/*	... */

		/* Do 1536 four_complex_fft operations */
		/*	distance between fft data elements is 256 */
		/*	do 768 macros each processing 16 data values */

#if 0 // @INSTR(,%xarch,<CORE>) NE 0; /* Core 2 faster doing separate loops */

		start_timer(4);
		rdi = (uintptr_t)g->xsincos_complex;	/* The sin/cos data */
		loops_init(512);			/* 4 iterations of 128 */
		r3_f3cl_three_complex_djbfft_preload;
		for(unsigned int loopA = 4; loopA; loopA--) {
			for(unsigned int loopB = 128; loopB; loopB--) {
				r3_f3cl_three_complex_djbfft(rsi, 64, 4*dist128, rdi, 512*XMM_SCD1);
				rdi += XMM_SCD1;		/* Next sin/cos data */
			}
			rsi += -128*64+dist128;	/* Next source pointer */
		}
		rsi = daddr; //rsi += -4*dist128;		/* Restore source pointer */
		end_timer(4);

		start_timer(5);
		rdi = (uintptr_t)g->xsincos_complex;	/* Ptr to the sin/cos data */
		rdx = rdi;
		pfing(rcx = (uintptr_t)g->data_prefetch); /* Prefetch pointer (we need to prefetch  */
		;					/* 1536 cache lines + 22 pad cache lines) */
		xtouch(xptr(rcx));			/* Load first prefetch TLB */
		loops_init(768);			/* 3 iters of 2 iters of 128 */
		r4_x2cl_2sc_four_complex_djbfft_preload;
		for(unsigned int loopA = 3; loopA; loopA--) {
			for(unsigned int loopB = 2; loopB; loopB--) {
				for(unsigned int loopC = 128; loopC; loopC--) {
					r4_x2cl_2sc_four_complex_djbfft(rsi, 64, 2*dist128, rdi, rdx);
					rdi += 3*XMM_SCD1;		/* Next sin/cos data */
					rdx += 6*XMM_SCD1;		/* Next sin/cos data */
				}
				rsi += -128*64+dist128;	/* Next source pointer */
			}
			rsi += -2*dist128+4*dist128; /* Next source pointer */
			rdi += -256*3*XMM_SCD1;	/* Restore sin/cos ptr */
			rdx += -256*6*XMM_SCD1;	/* Restore sin/cos ptr */
		}
		rsi = daddr; //rsi += -3*4*dist128;	/* Restore source pointer */
		end_timer(5);

#else /* Pentium 4 faster doing merged loops */

		start_timer(4);
		rdi = (uintptr_t)g->xsincos_complex;	/* The sin/cos data */
		rbp = rdi;
		rdx = rdi;
		pfing(rcx = (uintptr_t)g->data_prefetch); /* Prefetch pointer (we need to prefetch  */
		;					/* 1536 cache lines + 22 pad cache lines) */
		xtouch(xptr(rcx));			/* Load first prefetch TLB */
		loops_init_prefetch(1536, 64, 12, rcx); /* 256 iterations of 3 */

		for(unsigned int loopA = 2; loopA; loopA--) {
			xtouch(xptr(rcx+4096-64));		/* Load prefetch TLBs */
			for(unsigned int loopB = 128; loopB; loopB--) {
				r3_f3cl_three_complex_djbfft_preload;
				for(unsigned int loopC = 2; loopC; loopC--) {
					r3_f3cl_three_complex_djbfft(rsi, 2*dist128, 4*dist128, rdi, 512*XMM_SCD1);
					xxprefetch(u8ptr(rdi+XMM_SCD1));	/* Prefetch for next pass */
					xxprefetch(u8ptr(rdi+513*XMM_SCD1));	/* Prefetch for next pass */
					rdi += 256*XMM_SCD1;	/* Next sin/cos data */
				}
				rsi += -4*dist128;		/* Restore source pointer */
				loops_undo(6);			/* Rollback loop counter */

				r4_x2cl_2sc_four_complex_djbfft_preload;
				for(unsigned int loopC = 3; loopC; loopC--) {
					r4_x2cl_2sc_four_complex_djbfft(rsi, 4*dist128, 2*dist128, rbp, rdx);
				}

				rsi += -3*4*dist128+64;	/* Next source pointer */
				rdi += -2*256*XMM_SCD1+XMM_SCD1; /* Next sin/cos data */
				rbp += 3*XMM_SCD1;		/* Next sin/cos data */
				rdx += 6*XMM_SCD1;		/* Next sin/cos data */
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
		/*	0	+6K	256	+6K	512	+6K	768	+6K */
		/*	1	... */
		/*	... */
		/*	255	... */
		/*	1024	... */
		/*	... */
		/* On output the 64-byte cache lines hold these data values: */
		/*	0	+6K	64	+6K	128	+6K	192	+6K */
		/*	1	... */
		/*	... */
		/*	63	... */
		/*	256	... */
		/*	... */
		/* */
		/* Do 1536 four_complex_fft operations */
		/*	distance between fft data elements is 64 */
		/*	do 384 macros each processing 32 data values */

		/* Do FFT levels 7,8 */
		/* */
		/* On input the 64-byte cache lines hold these data values: */
		/*	0	+6K	64	+6K	128	+6K	192	+6K */
		/*	1	... */
		/*	... */
		/*	63	... */
		/*	256	... */
		/*	... */
		/* On output the 64-byte cache lines hold these data values: */
		/*	0	+6K	16	+6K	32	+6K	48	+6K */
		/*	1	... */
		/*	... */
		/*	15	... */
		/*	64	... */
		/*	... */
		/* */
		/* Do 1536 four_complex_fft operations */
		/*	distance between fft data elements is 16 */
		/*	do 384 macros each processing 32 data values */

		/* Do FFT levels 9,10 */
		/* */
		/* On input the 64-byte cache lines hold these data values: */
		/*	0	+6K	16	+6K	32	+6K	48	+6K */
		/*	1	... */
		/*	... */
		/*	15	... */
		/*	64	... */
		/*	... */
		/* On output the 64-byte cache lines hold these data values: */
		/*	0	+6K	4	+6K	8	+6K	12	+6K */
		/*	1	... */
		/*	... */
		/*	3	... */
		/*	16	... */
		/*	... */
		/* */
		/* Do 1536 four_complex_fft operations */
		/*	distance between fft data elements is 4 */
		/*	do 384 macros each processing 32 data values */

		/* Do this 64 cache lines at a time to maximize accessing data */
		/* while it is in the L1 cache. */

		start_timer(6);
		rdi += 1536*XMM_SCD1;	/* Load sin/cos pointers */
		rbx = rdi;
		rbp = rdi;
		loops_init_prefetch(384, 64, 2, rcx); /* 6 iters of 4 of 64 cache lines */
		r4_x4cl_four_complex_djbfft_preload;

		for(unsigned int loopA = 6; loopA; loopA--) {
			xtouch(xptr(rcx+4096-64));		/* Load prefetch TLBs */
			xtouch(xptr(rcx+2*4096-64));

			for(unsigned int loopB = 4; loopB; loopB--) {
				for(unsigned int loopC = 16; loopC; loopC--) {
					xxprefetch(u8ptr(rdi+4*XMM_SCD2));	/* Prefetch for next pass */
					r4_x4cl_four_complex_djbfft(rsi, 4*64, 64*64, dist128, rdi);
					rdi += 4*XMM_SCD2;		/* Next sine/cosine pointer */
				}
				rsi += -16*4*64;		/* Restore source pointer */
				rdi += -16*4*XMM_SCD2;	/* Restore sine/cosine pointer */
				loops_undo(16);			/* Back up loop counter */

				for(unsigned int loopC = 2; loopC; loopC--) {
					for(unsigned int loopD = 2; loopD; loopD--) {
						for(unsigned int loopE = 4; loopE; loopE--) {
							r4_x4cl_four_complex_djbfft(rsi, 4*64, 16*64, 32*64, rbx);
							rbx += 4*4*XMM_SCD2;	/* Next sine/cosine pointer */
						}
						rsi += -4*4*64+64*64;	/* Next source pointer */
						rbx += -4*4*4*XMM_SCD2;	/* Restore sine/cosine pointer */
					}
					rsi += -2*64*64+dist128;	/* Next source pointer */
				}
				rsi += -2*dist128;		/* Restore source pointer */
				loops_undo(2*2*4);		/* Back up loop counter */

				for(unsigned int loopC = 2; loopC; loopC--) {
					for(unsigned int loopD = 8; loopD; loopD--) {
						r4_x4cl_four_complex_djbfft(rsi, 16*64, 4*64, 8*64, rbp);
					}
					rsi += -8*16*64+dist128;	/* Next source pointer */
				}

				rsi += -2*dist128+64;	/* Next source pointer */
				rdi += XMM_SCD2;		/* Next sine/cosine pointer */
				rbx += 4*XMM_SCD2;		/* Next sine/cosine pointer */
				rbp += 16*XMM_SCD2;	/* Next sine/cosine pointer */
			}
			rsi += -4*64+2*dist128;	/* Next source pointer */
			rdi += -4*XMM_SCD2;	/* Restore sine/cosine pointer */
			rbx += -4*4*XMM_SCD2;	/* Restore sine/cosine pointer */
			rbp += -4*16*XMM_SCD2;	/* Restore sine/cosine pointer */
		}
		rsi = daddr; //rsi += -6*2*dist128;	/* Restore source pointer */
		end_timer(6);

		/* Do FFT levels 11,12 as well as inverse FFT levels 11,12 */
		/* */
		/* On input the 64-byte cache lines hold these data values: */
		/*	0	+6K	4	+6K	8	+6K	12	+6K */
		/*	1	... */
		/*	... */
		/*	3	... */
		/*	16	... */
		/*	... */
		/* On output the 64-byte cache lines hold these data values: */
		/*	0	+6K	2	+6K	4	+6K	6	+6K */
		/*	1	... */
		/*	8	... */
		/*	... */

		/* Do 1536 four_complex_fft operations */
		/*	distance between fft data elements is 1 */
		/*	do 384 macros each processing 32 data values */

		/* Execute the right middle step */

		if likely(g->ffttype & 2) {
			if likely(!(g->ffttype & 1)) { // xpass2_3072_levels_complex_2:;
				start_timer(9);
				loops_init_prefetch(384, 64, 4, rcx);
				xtouch(xptr(rcx+4096-64));		/* Load prefetch TLBs */
				xtouch(xptr(rcx+2*4096-64));
				r4_x4cl_four_complex_with_square_preload;
				for(unsigned int loopA = 12; loopA; loopA--) {
					for(unsigned int loopB = 32; loopB; loopB--) {
						r4_x4cl_four_complex_with_square(rsi, 4*64, 64, 2*64);
					}
					rsi += -32*4*64+dist128;	/* Next source pointer */
				}
				//rsi += -12*dist128;	/* Restore source pointer */
				end_timer(9);
			}else{ // xpass2_3072_levels_complex_3:;
				start_timer(9);
				rbp = g->DIST_TO_MULSRCARG;
				loops_init_prefetch(384, 64, 4, rcx);
				xtouch(xptr(rcx+4096-64));		/* Load prefetch TLBs */
				xtouch(xptr(rcx+2*4096-64));
				for(unsigned int loopA = 12; loopA; loopA--) {
					for(unsigned int loopB = 32; loopB; loopB--) {
						r4_x4cl_four_complex_with_mult(rsi, 4*64, 64, 2*64);
					}
					rsi += -32*4*64+dist128;	/* Next source pointer */
				}
				//rsi += -12*dist128;	/* Restore source pointer */
				end_timer(9);
			}
		}else{ // xpass2_3072_levels_complex_1:;
			start_timer(9);
			loops_init_prefetch(384, 128, 1, rcx);
			for(unsigned int loopA = 12; loopA; loopA--) {
				xtouch(xptr(rcx+4096-64));		/* Load prefetch TLBs */
				for(unsigned int loopB = 32; loopB; loopB--) {
					r4_x4cl_four_complex_fft_final(rsi, 4*64, 64, 2*64);
				}
				rsi += -32*4*64+dist128;	/* Next source pointer */
			}
			end_timer(9);
			goto	xpass2_3072_levels_complex_done;
		}
	}else{ // xpass2_3072_levels_complex_4:;
		start_timer(9);
		pfing(rcx = (uintptr_t)g->data_prefetch); /* Prefetch pointer */
		xtouch(xptr(rcx));			/* Load first prefetch TLB */
		rbp = g->DIST_TO_MULSRCARG;
		loops_init_prefetch(384, 128, 1, rcx);
		for(unsigned int loopA = 12; loopA; loopA--) {
			xtouch(xptr(rcx+4096-128));		/* Load prefetch TLBs */
			for(unsigned int loopB = 32; loopB; loopB--) {
				r4_x4cl_four_complex_with_mulf(rsi, 4*64, 64, 2*64);
			}
			rsi += -32*4*64+dist128;	/* Next source pointer */
		}
		//rsi += -12*dist128;	/* Restore source pointer */
		pfing(rcx += -96*64);		/* Back up prefetch pointer */
		end_timer(9);
	}
	rsi = daddr;

	/* Do inverse FFT levels 9,10 */
	/* On input the 64-byte cache lines hold these data values: */
	/*	0	+6K	2	+6K	4	+6K	6	+6K */
	/*	1	... */
	/*	8	... */
	/*	... */
	/* On output the 64-byte cache lines hold these data values: */
	/*	0	+6K	8	+6K	16	+6K	24	+6K */
	/*	1	... */
	/*	... */
	/*	7	... */
	/*	32	... */
	/*	... */
	/* */
	/* Do 1536 four_complex_unfft operations */
	/*	distance between fft data elements is 4 */
	/*	do 384 macros each processing 32 data values */

	/* Do inverse FFT levels 7,8 */
	/* On input the 64-byte cache lines hold these data values: */
	/*	0	+6K	8	+6K	16	+6K	24	+6K */
	/*	1	... */
	/*	... */
	/*	7	... */
	/*	32	... */
	/*	... */
	/* On output the 64-byte cache lines hold these data values: */
	/*	0	+6K	32	+6K	64	+6K	96	+6K */
	/*	1	... */
	/*	... */
	/*	31	... */
	/*	128	... */
	/*	... */
	/* */
	/* Do 1536 four_complex_unfft operations */
	/*	distance between fft data elements is 16 */
	/*	do 384 macros each processing 32 data values */

	/* Do inverse FFT levels 5,6 */
	/* On input the 64-byte cache lines hold these data values: */
	/*	0	+6K	32	+6K	64	+6K	96	+6K */
	/*	1	... */
	/*	... */
	/*	31	... */
	/*	128	... */
	/*	... */
	/* On output the 64-byte cache lines hold these data values: */
	/*	0	+6K	128	+6K	256	+6K	384	+6K */
	/*	1	... */
	/*	... */
	/*	127	... */
	/*	512	... */
	/*	... */
	/* */
	/* Do 1536 four_complex_unfft operations */
	/*	distance between fft data elements is 64 */
	/*	do 384 macros each processing 32 data values */

	/* Do this 64 cache lines at a time to maximize accessing data */
	/* while it is in the L1 cache. */

	start_timer(10);
	rdi = (uintptr_t)g->xsincos_complex;	/* Load sin/cos pointers */
	rdi += 1536*XMM_SCD1;
	rbx = rdi;
	rbp = rdi;
	loops_init_prefetch(384, 64, 2, rcx); /* 12*2 iters of 64 cache lines */
	r4_x4cl_four_complex_djbunfft_preload;

	for(unsigned int loopA = 12; loopA; loopA--) {
		xtouch(xptr(rcx+4096-64));		/* Load prefetch TLBs */

		for(unsigned int loopB = 2; loopB; loopB--) {
			for(unsigned int loopC = 16; loopC; loopC--) {
				r4_x4cl_four_complex_djbunfft(rsi, 8*64, 2*64, 4*64, rbp, 2*16*XMM_SCD2);
			}
			rsi += -16*8*64;		/* Restore source pointer */
			loops_undo(16);			/* Restore loop counter */

			for(unsigned int loopC = 4; loopC; loopC--) {
				for(unsigned int loopD = 4; loopD; loopD--) {
					r4_x4cl_four_complex_djbunfft(rsi, 2*64, 8*64, 16*64, rbx, 8*4*XMM_SCD2);
					rbx += 2*4*XMM_SCD2;	/* Next sine/cosine pointer */
				}
				rsi += -4*2*64+32*64;	/* Next source pointer */
				rbx += -4*2*4*XMM_SCD2;	/* Restore sine/cosine pointer */
			}
			rsi += -4*32*64;		/* Restore source pointer */
			loops_undo(4*4);			/* Restore loop counter */

			for(unsigned int loopC = 16; loopC; loopC--) {
				r4_x4cl_four_complex_djbunfft(rsi, 2*64, 32*64, 64*64, rdi, 32*1*XMM_SCD2);
				rdi += 2*XMM_SCD2;		/* Next sine/cosine pointer */
			}
			rdi += -16*2*XMM_SCD2;	/* Restore sine/cosine pointer */

			rsi += -16*2*64+64;	/* Next source pointer */
			rdi += XMM_SCD2;		/* Next sine/cosine pointer */
			rbx += 4*XMM_SCD2;		/* Next sine/cosine pointer */
			rbp += 16*XMM_SCD2;	/* Next sine/cosine pointer */
		}
		rsi += -2*64+dist128;	/* Next source pointer */
		rdi += -2*XMM_SCD2;	/* Restore sine/cosine pointer */
		rbx += -2*4*XMM_SCD2;	/* Restore sine/cosine pointer */
		rbp += -2*16*XMM_SCD2;	/* Restore sine/cosine pointer */
	}
	rsi = daddr; //rsi += -12*dist128;	/* Restore source pointer */
	end_timer(10);

	/* Do inverse FFT levels 3,4 */
	/* On input the 64-byte cache lines hold these data values: */
	/*	0	+6K	128	+6K	256	+6K	384	+6K */
	/*	1	... */
	/*	... */
	/*	127	... */
	/*	512	... */
	/*	... */
	/* On output the 64-byte cache lines hold these data values: */
	/*	0	+6K	512	+6K	1024	+6K	1536	+6K */
	/*	1	... */
	/*	... */
	/*	511	... */
	/*	2048	... */
	/*	... */

	/* Do 1536 four_complex_unfft operations */
	/*	distance between fft data elements is 256 */
	/*	do 384 macros each processing 32 data values */

	/* Do inverse FFT levels 1,2 (actually 1.585 FFT levels) */
	/* On input the 64-byte cache lines hold these data values: */
	/*	0	+6K	512	+6K	1024	+6K	1536	+6K */
	/*	1	... */
	/*	... */
	/*	511	... */
	/*	2048	... */
	/*	... */
	/* On output the 64-byte cache lines hold these data values if pass 1 swizzles: */
	/*	0	+6K	1536	+6K	3072	+6K	4608	+6K */
	/*	1	... */
	/*	... */
	/*	1535 */
	/* On output the 64-byte cache lines hold these data values if pass 2 swizzles: */
	/*	0	+1.5K	3072	+1.5K	6144	+1.5K	9216	+1.5K */
	/*	1	... */
	/*	... */
	/*	1535 */

	/* Do 2048 three_complex_unfft operations */
	/*	distance between fft data elements is 1024 */
	/*	do 512 macros each processing 24 data values */

#if 0 //@INSTR(,%xarch,<CORE>) NE 0; /* Core 2 faster doing separate loops */

	start_timer(13);
	rdi = (uintptr_t)g->xsincos_complex;	/* Load sin/cos pointer */
	rdx = rdi;
	loops_init_prefetch(384, 64, 2, rcx);/* 3 iters of 128 */
	xtouch(xptr(rcx+4096-64));		/* Load prefetch TLBs */
	xtouch(xptr(rcx+2*4096-64));
	xtouch(xptr(rcx+3*4096-64));
	r4_x4cl_2sc_four_complex_djbunfft_preload;
	for(unsigned int loopA = 3; loopA; loopA--) {
		for(unsigned int loopB = 128; loopB; loopB--) {
			r4_x4cl_2sc_four_complex_djbunfft(rsi, 64, dist128, 2*dist128, rdi, 128*3*XMM_SCD1, rdx, 128*6*XMM_SCD1);
			rdi += 3*XMM_SCD1;		/* Next sine/cosine pointer */
			rdx += 6*XMM_SCD1;		/* Next sine/cosine pointer */
		}
		rsi += -128*64+4*dist128;	/* Next source pointer */
		rdi += -128*3*XMM_SCD1;	/* Reset sine/cosine pointer */
		rdx += -128*6*XMM_SCD1;	/* Reset sine/cosine pointer */
	}
	rsi = daddr; //rsi += -3*4*dist128;	/* Restore source pointer */
	end_timer(13);

	start_timer(14);
	rdi = (uintptr_t)g->xsincos_complex;	/* Load sin/cos pointer */
	pfing(rcx += -10*64);		/* Back up prefetch pointer by 10 cache lines */
	loops_init_prefetch(512, 64, 4, rcx); /* 4 iters of 128 */
	xtouch(xptr(rcx+4096-64));		/* Load prefetch TLBs */
	xtouch(xptr(rcx+2*4096-64));
	r3_x3cl_three_complex_djbunfft_preload;
	for(unsigned int loopA = 4; loopA; loopA--) {
		for(unsigned int loopB = 128; loopB; loopB--) {
			r3_x3cl_three_complex_djbunfft(rsi, 64, 4*dist128, rdi, 512*XMM_SCD1);
			rdi += XMM_SCD1;		/* Next set of sin/cos data */
		}
		rsi += -128*64+dist128;	/* Next source pointer */
	}
	end_timer(14);

#else /* Pentium 4 faster doing merged loops */

	start_timer(13);
	rdi = (uintptr_t)g->xsincos_complex;	/* Load sin/cos pointer */
	rbp = rdi;
	rdx = rdi;
	pfing(rcx += -74*64);		/* Back up prefetch pointer by 74 cache lines */
	loops_init_prefetch(1536, 64, 12, rcx); /* 4 iters of 128, prefetching 128 cache lines each loop */

	for(unsigned int loopA = 128; loopA; loopA--) {
		xtouch(xptr(rcx+4096-64));		/* Load prefetch TLB */

		r4_x4cl_2sc_four_complex_djbunfft_preload;
		for(unsigned int loopB = 3; loopB; loopB--) {
			r4_x4cl_2sc_four_complex_djbunfft(rsi, 4*dist128, dist128, 2*dist128, rbp, 128*3*XMM_SCD1, rdx, 128*6*XMM_SCD1);
		}
		rsi += -3*4*dist128;	/* Restore source pointer */
		loops_undo(12);

		r3_x3cl_three_complex_djbunfft_preload;
		for(unsigned int loopB = 4; loopB; loopB--) {
			r3_x3cl_three_complex_djbunfft(rsi, dist128, 4*dist128, rdi, 512*XMM_SCD1);
			rdi += 128*XMM_SCD1;	/* Next set of sin/cos data */
		}

		rsi += -4*dist128+64;	/* Next source pointer */
		rdi += -4*128*XMM_SCD1+XMM_SCD1; /* Next sin/cos ptr */
		rbp += 3*XMM_SCD1;		/* Next sine/cosine pointer */
		rdx += 6*XMM_SCD1;		/* Next sine/cosine pointer */
	}
	end_timer(13);

#endif

	xpass2_3072_levels_complex_done:;
}while(g->pass2_get_next_block(g) == 0);
}
