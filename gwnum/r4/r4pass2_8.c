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

#include "r4.h"

/* *************** radix4-pass2-8-levels-real macro ****************** */
/* This macro takes 512 real values and performs the final 9 levels of */
/* the FFT.  Independently, this macro also takes 256 complex values and */
/* performs the final 8 levels of the FFT.  The results are squared, and */
/* the first 8/9 levels of the inverse FFT are done. */

void r4_pass2_8_levels_complex(struct gwasm_data *__restrict g);
void xpass2_r4_8_levels(struct gwasm_data *__restrict g) {
	g->thread_work_routine = r4_pass2_8_levels_complex;
	pass2_wake_up_threads(g);
	if(g->NEGACYCLIC_FFT == 0) {

//static inline void r4_pass2_8_levels_real(struct gwasm_data *__restrict g) {
	uintptr_t rdx,rbx,rbp,rsi,rdi,daddr;
	vec2f64 xmm0,xmm1,xmm2,xmm3,xmm4,xmm5,xmm6,xmm7;

	/* We switch to the same format used in one-pass FFTs (that is, the */
	/* high word of xmm registers being 1 greater than the low word): */

	/* Type 4 FFTs skip the forward FFT process */

	start_timer(2);
	daddr = (uintptr_t)g->data_addr;/* Load source address */
	rsi = daddr;
	rbx = g->DIST_TO_FFTSRCARG;
	if likely(!(g->ffttype & 4)) {

		/* Do FFT level 1,2 */
		/* */
		/* On input the 64-byte cache lines hold these data values if pass 2 swizzles: */
		/*	0	+128	256	+128	512	+128	768	+128 */
		/*	1	... */
		/*	... */
		/*	127 */
		/* On input the 64-byte cache lines hold these data values if pass 1 swizzles: */
		/*	0	+512	128	+512	256	+512	384	+512 */
		/*	1	... */
		/*	... */
		/*	127	... */
		/* On output the 64-byte cache lines hold these data values: */
		/*	0	+512	64	+512	128	+512	192	+512 */
		/*	1	... */
		/*	... */
		/*	63	... */
		/*	256	... */
		/*	... */

		/* Do 64 eight_reals_fft operations */
		/* Do 64 four_complex_djbfft operations */
		/*	distance between fft data elements is 64 */
		/*	do 64 macros each processing 16 data values */

		rdi = (uintptr_t)g->xsincos_complex;	/* Ptr to the sin/cos data */
		rdx = (uintptr_t)g->sincos3;/* Ptr to real sin/cos data */
		loops_init(64);
		for(unsigned int loopA = 64; loopA; loopA--) {
			r4_fh2cl_eight_reals_four_complex_djbfft(rsi, 64, 64*64, rdi, rdx);
			rdi += XMM_SCD2;		/* Next sin/cos data */
			rdx += XMM_HSCD2;		/* Next sin/cos data */
		}
		rsi = daddr; //rsi += -64*64;		/* Restore source pointer */
		rdi += -64*XMM_SCD2;	/* Restore sin/cos pointer */
		rdx += -64*XMM_HSCD2;	/* Restore sin/cos pointer */

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

		/* Do 16 eight_reals_fft operations */
		/* Do 48 four_complex_djbfft operations */
		/*	distance between fft data elements is 16 */
		/*	do 16 macros each processing 32 data values */

		loops_init_prefetch(32, 64, 2, rcx);
		pfing(rcx = (uintptr_t)g->data_prefetch);/* Prefetch pointer */
		xtouch(xptr(rcx));			/* Preload first TLB */
		xtouch(xptr(rcx+4096-64));		/* Preload the TLBs */
		xtouch(xptr(rcx+2*4096-64));
		r4_h4cl_eight_reals_four_complex_djbfft_preload;
		for(unsigned int loopA = 16; loopA; loopA--) {
			r4_h4cl_eight_reals_four_complex_djbfft(rsi, 64, 16*64, 32*64, rdi, rdx);
			rdi += 4*XMM_SCD2;		/* Next sin/cos data */
			rdx += 4*XMM_HSCD2;	/* Next sin/cos data */
		}
		rsi += -16*64+64*64;	/* Next source pointer */
		rdi += -16*4*XMM_SCD2;	/* Restore sin/cos ptr */
		rdx += -16*4*XMM_HSCD2;	/* Restore sin/cos ptr */

		/* Do 64 four_complex_djbfft operations */
		/*	distance between fft data elements is 16 */
		/*	do 16 macros each processing 32 data values */

		loops_reset;			/* 1 more iters of 16 */
		r4_x4cl_four_complex_djbfft_preload;
		{
			for(unsigned int loopB = 16; loopB; loopB--) {
				r4_x4cl_four_complex_djbfft(rsi, 64, 16*64, 32*64, rdi);
				rdi += 4*XMM_SCD2;		/* Next sin/cos data */
			}
			rsi += -16*64+64*64;	/* Next source pointer */
			rdi += -16*4*XMM_SCD2;	/* Restore sin/cos ptr */
		}
		rsi = daddr; //rsi += -2*64*64;		/* Next source pointer */

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

		/* Do 4 eight_reals_fft operations */
		/* Do 12 four_complex_djbfft operations */
		/*	distance between fft data elements is 4 */
		/*	do 4 macros each processing 32 data values */

		loops_init_prefetch(32, 64, 2, rcx);
		r4_h4cl_eight_reals_four_complex_djbfft_preload;
		for(unsigned int loopA = 4; loopA; loopA--) {
			r4_h4cl_eight_reals_four_complex_djbfft(rsi, 64, 4*64, 8*64, rdi, rdx);
			rdi += 16*XMM_SCD2;	/* Next sin/cos data */
			rdx += 16*XMM_HSCD2;	/* Next sin/cos data */
		}
		rsi += -4*64+16*64;	/* Next source pointer */
		rdi += -4*16*XMM_SCD2;	/* Restore sin/cos ptr */
		rdx += -4*16*XMM_HSCD2;	/* Restore sin/cos ptr */

		/* Do 112 four_complex_djbfft operations */
		/*	distance between fft data elements is 4 */
		/*	do 28 macros each processing 32 data values */

		loops_reset;			/* 7 more iters of 4 */
		r4_x4cl_four_complex_djbfft_preload;
		for(unsigned int loopA = 7; loopA; loopA--) {
			for(unsigned int loopB = 4; loopB; loopB--) {
				r4_x4cl_four_complex_djbfft(rsi, 64, 4*64, 8*64, rdi);
				rdi += 16*XMM_SCD2;	/* Next sin/cos data */
			}
			rsi += -4*64+16*64;	/* Next source pointer */
			rdi += -4*16*XMM_SCD2;	/* Restore sin/cos ptr */
		}
		rsi = daddr; //rsi += -8*16*64;		/* Next source pointer */

		/* Do FFT levels 7,8 */
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

		/* Do 1 eight_reals_fft operations */
		/* Do 3 four_complex_fft operations */
		/*	distance between fft data elements is 1 */

		/* Do 124 four_complex_fft operations */
		/*	distance between fft data elements is 1 */
		/*	do 31 macros each processing 32 data values */

		/* Execute the proper middle step */

		if likely(g->ffttype & 2) {
			if likely(!(g->ffttype & 1)) { //xpass2_8_levels_real_2:;
				r4_h4cl_eight_reals_four_complex_with_square(rsi, 4*64, 64, 2*64);
				loops_init_prefetch(32, 64, 2, rcx, 1, 32);
				r4_x4cl_four_complex_with_square_preload;
				for(unsigned int loopA = 31; loopA; loopA--) {
					r4_x4cl_four_complex_with_square(rsi, 4*64, 64, 2*64);
				}
				//rsi += -32*4*64;		/* Next source pointer */
			}else{ //xpass2_8_levels_real_3:;
				rbp = g->DIST_TO_MULSRCARG;
				r4_h4cl_eight_reals_four_complex_with_mult(rsi, 4*64, 64, 2*64);
				loops_init_prefetch(32, 64, 2, rcx, 1, 32);
				for(unsigned int loopA = 31; loopA; loopA--) {
					r4_x4cl_four_complex_with_mult(rsi, 4*64, 64, 2*64);
				}
				//rsi += -32*4*64;		/* Next source pointer */
			}
		}else{ //xpass2_8_levels_real_1:;
			r4_h4cl_eight_reals_four_complex_fft_final(rsi, 4*64, 64, 2*64);
			loops_init_prefetch(32, 64, 1, rcx, 1, 32);
			for(unsigned int loopA = 31; loopA; loopA--) {
				r4_x4cl_four_complex_fft_final(rsi, 4*64, 64, 2*64);
			}
			goto	xpass2_8_real_done;
		}
	}else{ //xpass2_8_levels_real_4:;
		rbp = g->DIST_TO_MULSRCARG;
		pfing(rcx = (uintptr_t)g->data_prefetch);/* Load prefetch pointer */
		xtouch(xptr(rcx));			/* Preload the TLBs */
		xtouch(xptr(rcx+4096-128));
		xtouch(xptr(rcx+2*4096-128));
		r4_h4cl_eight_reals_four_complex_with_mulf(rsi, 4*64, 64, 2*64);
		loops_init_prefetch(32, 64, 1, rcx, 1, 32);
		for(unsigned int loopA = 31; loopA; loopA--) {
			r4_x4cl_four_complex_with_mulf(rsi, 4*64, 64, 2*64);
		}
		//rsi += -32*4*64;		/* Next source pointer */
	}
	rsi = daddr;

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

	/* Do 4 eight_reals_unfft operations */
	/* Do 4 four_complex_djbunfft operations */
	/*	distance between fft data elements is 4 */
	/*	do 2 macros processing 32 data values */
	rdi = (uintptr_t)g->xsincos_complex;	/* Load sin/cos pointers */
	rdx = (uintptr_t)g->sincos3;
	loops_init_prefetch(32, 64, 1, rcx);
	for(unsigned int loopA = 2; loopA; loopA--) {
		r4_h4cl_eight_reals_four_complex_djbunfft(rsi, 64, 2*64, 4*64, rdi, 2*16*XMM_SCD2, rdx, 2*16*XMM_HSCD2);
		rdi += 16*XMM_SCD2;	/* Next sin/cos data */
		rdx += 16*XMM_HSCD2;	/* Next sin/cos data */
	}
	rsi += -2*64+8*64;		/* Next source pointer */
	rdi += -2*16*XMM_SCD2;	/* Restore sin/cos ptr */
	rdx += -2*16*XMM_HSCD2;	/* Restore sin/cos ptr */

	/* Do 120 four_complex_djbunfft operations */
	/*	distance between fft data elements is 4 */
	/*	do 30 macros each processing 32 data values */

	loops_reset;			/* 4 iters of 15 or 16 of 2 -- BUG in comment */
	r4_x4cl_four_complex_djbunfft_preload;
	for(unsigned int loopA = 15; loopA; loopA--) {
		for(unsigned int loopB = 2; loopB; loopB--) {
			r4_x4cl_four_complex_djbunfft(rsi, 64, 2*64, 4*64, rdi, 2*16*XMM_SCD2);
			rdi += 16*XMM_SCD2;	/* Next sine/cosine pointer */
		}
		rsi += -2*64+8*64;		/* Next source pointer */
		rdi += -2*16*XMM_SCD2;	/* Restore sin/cos ptr */
	}
	rsi = daddr; //rsi += -16*8*64;		/* Next source pointer */

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

	/* Do 16 eight_reals_unfft operations */
	/* Do 16 four_complex_djbunfft operations */
	/*	distance between fft data elements is 16 */
	/*	do 8 macros processing 32 data values */

	loops_init_prefetch(32, 64, 1, rcx);
	for(unsigned int loopA = 8; loopA; loopA--) {
		r4_h4cl_eight_reals_four_complex_djbunfft(rsi, 64, 8*64, 16*64, rdi, 8*4*XMM_SCD2, rdx, 8*4*XMM_HSCD2);
		rdi += 4*XMM_SCD2;		/* Next sin/cos data */
		rdx += 4*XMM_HSCD2;	/* Next sin/cos data */
	}
	rsi += -8*64+32*64;	/* Next source pointer */
	rdi += -8*4*XMM_SCD2;	/* Restore sin/cos ptr */
	rdx += -8*4*XMM_HSCD2;	/* Restore sin/cos ptr */

	/* Do 96 four_complex_djbunfft operations */
	/*	distance between fft data elements is 16 */
	/*	do 24 macros each processing 32 data values */

	loops_reset;			/* 3 more iters of 8 */
	r4_x4cl_four_complex_djbunfft_preload;
	for(unsigned int loopA = 3; loopA; loopA--) {
		for(unsigned int loopB = 8; loopB; loopB--) {
			r4_x4cl_four_complex_djbunfft(rsi, 64, 8*64, 16*64, rdi, 8*4*XMM_SCD2);
			rdi += 4*XMM_SCD2;		/* Next sin/cos data */
		}
		rsi += -8*64+32*64;	/* Next source pointer */
		rdi += -8*4*XMM_SCD2;	/* Restore sin/cos ptr */
	}
	rsi = daddr; //rsi += -4*32*64;		/* Next source pointer */

	/* Do inverse FFT level 1,2 */
	/* On input the 64-byte cache lines hold these data values: */
	/*	0	+512	32	+512	64	+512	96	+512 */
	/*	1	... */
	/*	... */
	/*	31	... */
	/*	128	... */
	/*	... */
	/* On output the 64-byte cache lines hold these data values if pass 1 swizzles: */
	/*	0	+512	128	+512	256	+512	384	+512 */
	/*	1	... */
	/*	... */
	/*	127 */
	/* On output the 64-byte cache lines hold these data values if pass 2 swizzles: */
	/*	0	+128	256	+128	512	+128	768	+128 */
	/*	1	... */
	/*	... */
	/*	127 */

	/* Do 64 eight_reals_unfft operations */
	/* Do 64 four_complex_djbunfft operations */
	/*	distance between fft data elements is 64 */
	/*	do 32 macros processing 32 data values */

	rdi = (uintptr_t)g->xsincos_complex;	/* Load sin/cos pointers */
	loops_init_prefetch(32, 64, 2, rcx);
	for(unsigned int loopA = 32; loopA; loopA--) {
		r4_h4cl_eight_reals_four_complex_djbunfft(rsi, 64, 32*64, 64*64, rdi, 32*XMM_SCD2, rdx, 32*XMM_HSCD2);
		rdi += XMM_SCD2;		/* Next sin/cos data */
		rdx += XMM_HSCD2;		/* Next sin/cos data */
	}

	xpass2_8_real_done:;
	end_timer(2);


		if(g->pass2_get_next_block(g) != 0) return;
	}
	r4_pass2_8_levels_complex(g);
}


/* *************** radix4-pass2-8-levels-complex macro ****************** */
/* This macro takes 2 groups of 256 complex values and */
/* performs the final 8 levels of a standard mostly radix-4 FFT, squares the */
/* results, and does the first 8 levels of the inverse FFT. */

void r4_pass2_8_levels_complex(struct gwasm_data *__restrict g) {
do{
	uintptr_t rbx,rbp,rsi,rdi,daddr;
	vec2f64 xmm0,xmm1,xmm2,xmm3,xmm4,xmm5,xmm6,xmm7;

	/* Type 4 FFTs skip the forward FFT process */

	daddr = (uintptr_t)g->data_addr;/* Load address of FFT data */
	rsi = daddr;
	rbx = g->DIST_TO_FFTSRCARG;
	if likely(!(g->ffttype & 4)) {

		/* Do FFT levels 1,2 */
		/* */
		/* We swizzle the data so that one set of data is in the high word of */
		/* an xmm register, the other set is in the low word.  This means */
		/* a distance of 512 data values. */
		/* */
		/* On input the 64-byte cache lines hold these data values if pass 2 swizzles: */
		/*	0	+128	256	+128	512	+128	768	+128 */
		/*	1	... */
		/*	... */
		/*	127 */
		/* On input the 64-byte cache lines hold these data values if pass 1 swizzles: */
		/*	0	+512	128	+512	256	+512	384	+512 */
		/*	1	... */
		/*	... */
		/*	127	... */
		/* On output the 64-byte cache lines hold these data values: */
		/*	0	+512	64	+512	128	+512	192	+512 */
		/*	1	... */
		/*	... */
		/*	63	... */
		/*	256	... */
		/*	... */

		/* Do 128 four_complex_fft operations */
		/*	distance between fft data elements is 64 */
		/*	do 64 macros each processing 16 data values */

		start_timer(4);
		rdi = (uintptr_t)g->xsincos_complex;	/* The sin/cos data */
		loops_init(64);			/* 64 iterations */
		r4_f2cl_four_complex_djbfft_preload;
		for(unsigned int loopA = 64; loopA; loopA--) {
			xxprefetch(u64ptr(rdi+XMM_SCD2));	/* Prefetch for next pass */
			r4_f2cl_four_complex_djbfft(rsi, 64, 64*64, rdi);
			rdi += XMM_SCD2;		/* Next sin/cos data */
		}
		rsi = daddr; //rsi += -64*64;		/* Restore source pointer */
		rdi += -64*XMM_SCD2;	/* Restore sin/cos data */
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
		/* */
		/* Do 128 four_complex_fft operations */
		/*	distance between fft data elements is 16 */
		/*	do 32 macros each processing 32 data values */

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
		/* */
		/* Do 128 four_complex_fft operations */
		/*	distance between fft data elements is 4 */
		/*	do 32 macros each processing 32 data values */

		/* Do this 16 cache lines at a time to maximize accessing data */
		/* while it is in the L1 cache. */

		start_timer(6);
		pfing(rcx = (uintptr_t)g->data_prefetch);/* Prefetch pointer (we need to prefetch */
		;					/* 128 cache lines + 0 pad cache lines) */
		xtouch(xptr(rcx));			/* Load all the prefetch TLBs */
		xtouch(xptr(rcx+4096-64));
		xtouch(xptr(rcx+2*4096-64));
		rbx = rdi;		/* Copy sin/cos ptr */
		loops_init_prefetch(32, 64, 2, rcx); /* 2*4 iters of 16 cache lines */
		r4_x4cl_four_complex_djbfft_preload;

		for(unsigned int loopA = 2; loopA; loopA--) {
			for(unsigned int loopB = 4; loopB; loopB--) {
				for(unsigned int loopC = 4; loopC; loopC--) {
					r4_x4cl_four_complex_djbfft(rsi, 4*64, 16*64, 32*64, rdi);
					xxprefetch(u64ptr(rdi+4*XMM_SCD2));	/* Prefetch for next pass */
					rdi += 4*4*XMM_SCD2;	/* Next sine/cosine pointer */
				}
				rsi += -4*4*64;		/* Restore source pointer */
				rdi += -4*4*4*XMM_SCD2;	/* Restore sine/cosine pointer */
				loops_undo(4);			/* Back up loop counter */

				for(unsigned int loopC = 4; loopC; loopC--) {
					r4_x4cl_four_complex_djbfft(rsi, 16*64, 4*64, 8*64, rbx);
				}

				rsi += -4*16*64+64;	/* Next source pointer */
				rdi += 4*XMM_SCD2;		/* Next sine/cosine pointer */
				rbx += 16*XMM_SCD2;	/* Next sine/cosine pointer */
			}
			rsi += -4*64+64*64;	/* Next source pointer */
			rdi += -4*4*XMM_SCD2;	/* Restore sine/cosine pointer */
			rbx += -4*16*XMM_SCD2;	/* Restore sine/cosine pointer */
		}
		rsi = daddr; //rsi += -2*64*64;		/* Restore source pointer */
		end_timer(6);

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

		/* Do 128 four_complex_fft operations */
		/*	distance between fft data elements is 1 */
		/*	do 32 macros each processing 32 data values */

		/* Execute the right middle step */

		if likely(g->ffttype & 2) {
			if likely(!(g->ffttype & 1)) { //xpass2_8_levels_complex_2:;
				start_timer(9);
				loops_init_prefetch(32, 64, 2, rcx);
				r4_x4cl_four_complex_with_square_preload;
				for(unsigned int loopA = 32; loopA; loopA--) {
					r4_x4cl_four_complex_with_square(rsi, 4*64, 64, 2*64);
				}
				//rsi += -32*4*64;		/* Restore source pointer */
				end_timer(9);
			}else{ //xpass2_8_levels_complex_3:;
				start_timer(9);
				rbp = g->DIST_TO_MULSRCARG;
				loops_init_prefetch(32, 64, 2, rcx);
				for(unsigned int loopA = 32; loopA; loopA--) {
					r4_x4cl_four_complex_with_mult(rsi, 4*64, 64, 2*64);
				}
				//rsi += -32*4*64;		/* Restore source pointer */
				end_timer(9);
			}
		}else{ //xpass2_8_levels_complex_1:;
			start_timer(9);
			loops_init_prefetch(32, 64, 1, rcx);
			for(unsigned int loopA = 32; loopA; loopA--) {
				r4_x4cl_four_complex_fft_final(rsi, 4*64, 64, 2*64);
			}
			end_timer(9);
			goto	xpass2_8_levels_complex_done;
		}
	}else{ //xpass2_8_levels_complex_4:;
		start_timer(9);
		pfing(rcx = (uintptr_t)g->data_prefetch);/* Prefetch pointer */
		xtouch(xptr(rcx));			/* Load all prefetch TLB */
		xtouch(xptr(rcx+4096-64));
		xtouch(xptr(rcx+2*4096-64));
		rbp = g->DIST_TO_MULSRCARG;
		loops_init_prefetch(32, 64, 1, rcx);
		for(unsigned int loopA = 32; loopA; loopA--) {
			r4_x4cl_four_complex_with_mulf(rsi, 4*64, 64, 2*64);
		}
		//rsi += -32*4*64;		/* Restore source pointer */
		end_timer(9);
	}
	rsi = daddr;

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
	/* */
	/* Do 128 four_complex_unfft operations */
	/*	distance between fft data elements is 4 */
	/*	do 32 macros each processing 32 data values */

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
	/* */
	/* Do 128 four_complex_unfft operations */
	/*	distance between fft data elements is 16 */
	/*	do 32 macros each processing 32 data values */

	/* Do this 16 cache lines at a time to maximize accessing data */
	/* while it is in the L1 cache. */

	start_timer(10);
	rdi = (uintptr_t)g->xsincos_complex;	/* Load sin/cos pointers */
	rbx = rdi;
	loops_init_prefetch(32, 64, 1, rcx); /* 4*2 iters of 16 cache lines */
	r4_x4cl_four_complex_djbunfft_preload;

	for(unsigned int loopA = 4; loopA; loopA--) {
		for(unsigned int loopB = 2; loopB; loopB--) {
			for(unsigned int loopC = 4; loopC; loopC--) {
				r4_x4cl_four_complex_djbunfft(rsi, 8*64, 2*64, 4*64, rbx, 2*16*XMM_SCD2);
			}
			rsi += -4*8*64;		/* Restore source pointer */
			loops_undo(4);			/* Restore loop counter */

			for(unsigned int loopC = 4; loopC; loopC--) {
				r4_x4cl_four_complex_djbunfft(rsi, 2*64, 8*64, 16*64, rdi, 8*4*XMM_SCD2);
				rdi += 2*4*XMM_SCD2;	/* Next sine/cosine pointer */
			}

			rsi += -4*2*64+64;		/* Next source pointer */
			rdi += -4*2*4*XMM_SCD2+4*XMM_SCD2; /* Next sine/cosine pointer */
			rbx += 16*XMM_SCD2;	/* Next sine/cosine pointer */
		}
		rsi += -2*64+32*64;	/* Next source pointer */
		rdi += -2*4*XMM_SCD2;	/* Restore sine/cosine pointer */
		rbx += -2*16*XMM_SCD2;	/* Restore sine/cosine pointer */
	}
	rsi = daddr; //rsi += -4*32*64;		/* Restore source pointer */
	end_timer(10);

	/* Do inverse FFT levels 1,2 */
	/* On input the 64-byte cache lines hold these data values: */
	/*	0	+512	32	+512	64	+512	96	+512 */
	/*	1	... */
	/*	... */
	/*	31	... */
	/*	128	... */
	/*	... */
	/* On output the 64-byte cache lines hold these data values if pass 1 swizzles: */
	/*	0	+512	128	+512	256	+512	384	+512 */
	/*	1	... */
	/*	... */
	/*	127 */
	/* On output the 64-byte cache lines hold these data values if pass 2 swizzles: */
	/*	0	+128	256	+128	512	+128	768	+128 */
	/*	1	... */
	/*	... */
	/*	127 */

	/* Do 128 four_complex_unfft operations */
	/*	distance between fft data elements is 64 */
	/*	do 32 macros each processing 32 data values */

	start_timer(14);
	loops_init_prefetch(32, 64, 2, rcx); /* 32 iters */
	r4_x4cl_four_complex_djbunfft_preload;
	for(unsigned int loopA = 32; loopA; loopA--) {
		r4_x4cl_four_complex_djbunfft(rsi, 64, 32*64, 64*64, rdi, 32*XMM_SCD2);
		rdi += XMM_SCD2;		/* Next set of sin/cos data */
	}
	end_timer(14);

	xpass2_8_levels_complex_done:;
}while(g->pass2_get_next_block(g) == 0);
}
