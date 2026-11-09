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
#include "r5.h"

/* *************** radix4-pass2-25600-levels-real macro ****************** */
/* This macro takes 51200 real values and performs the final 15.644 levels of */
/* the FFT.  Independently, this macro also takes 25600 complex values and */
/* performs the final 14.644 levels of the FFT.  The results are squared, and */
/* the first 14.644/15.644 levels of the inverse FFT are done. */


void r4_pass2_25600_levels_complex(struct gwasm_data *__restrict g);
void xpass2_r4_25600_levels(struct gwasm_data *__restrict g) {
	g->thread_work_routine = r4_pass2_25600_levels_complex;
	pass2_wake_up_threads(g);
	if(g->NEGACYCLIC_FFT == 0) {
//void r4_pass2_25600_levels_real(struct gwasm_data *__restrict g)
	uintptr_t rdx,rbx,rbp,rsi,rdi,daddr;
	vec2f64 xmm0,xmm1,xmm2,xmm3,xmm4,xmm5,xmm6,xmm7;

	/* We switch to the same format used in one-pass FFTs (that is, the */
	/* high word of xmm registers being 1 greater than the low word): */

	/* Type 4 FFTs skip the forward FFT process */

	start_timer(2);
	daddr = (uintptr_t)g->data_addr;/* Load source address */
	rsi = daddr;
	rbx = g->DIST_TO_FFTSRCARG;
	if likely(!(g->ffttype & 4)) { //xpass2_25600_levels_real_4;

		/* Do FFT level 1,2 (actually 2.322 levels) */
		/* */
		/* On input the 64-byte cache lines hold these data values if pass 2 swizzles: */
		/*	0	+12800	25600	+12800	51200	+12800	76800	+12800 */
		/*	1	... */
		/*	... */
		/*	12799	... */
		/* On input the 64-byte cache lines hold these data values if pass 1 swizzles: */
		/*	0	+50K	12800	+50K	25600	+50K	38400	+50K */
		/*	1	... */
		/*	... */
		/*	12799	... */
		/* On output the 64-byte cache lines hold these data values: */
		/*	0	+50K	2560	+50K	5120	+50K	7680	+50K */
		/*	1	... */
		/*	... */
		/*	2559	... */
		/*	10240	... */
		/*	... */

		/* Do 5120 ten_reals_fft operations */
		/* Do 5120 five_complex_djbfft operations */
		/*	distance between fft data elements is 5120 */
		/*	do 2560 iterations */

		rdi = (uintptr_t)g->xsincos_complex;	/* Ptr to the sin/cos data */
		rdx = (uintptr_t)g->sincos3;/* Ptr to real sin/cos data */
		loops_init(2560);
		for(unsigned int loopA = 20; loopA; loopA--) {
			for(unsigned int loopB = 128; loopB; loopB--) {
				r5_fh5cl_ten_reals_five_complex_djbfft(rsi, 64, 20*dist128, rdi, 2560*XMM_SCD2, rdx, 2560*XMM_HSCD2);
				rdi += XMM_SCD2;		/* Next sin/cos data */
				rdx += XMM_HSCD2;		/* Next sin/cos data */
			}
			rsi += -128*64+dist128;	/* Next source pointer */
		}
		rsi = daddr; //rsi += -20*dist128;	/* Restore source pointer */
		rdi += -2560*XMM_SCD2;	/* Restore sin/cos ptr */
		rdx += -2560*XMM_HSCD2;	/* Next sin/cos ptr */

		/* Do FFT levels 3,4 (actually 2.322 levels) */
		/* */
		/* On input the 64-byte cache lines hold these data values: */
		/*	0	+50K	2560	+50K	5120	+50K	7680	+50K */
		/*	1	... */
		/*	... */
		/*	2559	... */
		/*	10240	... */
		/*	... */
		/* On output the 64-byte cache lines hold these data values: */
		/*	0	+50K	512	+50K	1024	+50K	1536	+50K */
		/*	1	... */
		/*	... */
		/*	511	... */
		/*	2048	... */
		/*	... */

		/* Do 1024 ten_reals_fft operations */
		/* Do 1024 five_complex_djbfft operations */
		/*	distance between fft data elements is 1024 */
		/*	do 512 macros each processing 40 data values */

		pfing(rcx = (uintptr_t)g->data_prefetch); /* Prefetch pointer */
		xtouch(xptr(rcx));			/* Preload first TLB */
		loops_init_prefetch(2560, 64, 16, rcx);
		xtouch(xptr(rcx+4096-64));		/* Preload the TLBs */
		xtouch(xptr(rcx+2*4096-64));
		xtouch(xptr(rcx+160*64-64));
		for(unsigned int loopA = 4; loopA; loopA--) {
			for(unsigned int loopB = 128; loopB; loopB--) {
				r5_h5cl_ten_reals_five_complex_djbfft(rsi, 64, 4*dist128, rdi, 512*5*XMM_SCD2, rdx, 512*5*XMM_HSCD2);
				rdi += 5*XMM_SCD2;		/* Next sin/cos data */
				rdx += 5*XMM_HSCD2;	/* Next sin/cos data */
			}
			rsi += -128*64+dist128;	/* Next source pointer */
		}
		rsi += -4*dist128+20*dist128; /* Next source pointer */
		rdi += -512*5*XMM_SCD2;	/* Restore sin/cos ptr */
		rdx += -512*5*XMM_HSCD2+5120*XMM_HSCD2; /* Next sin/cos ptr */

		/* Do 8192 five_complex_fft operations */
		/*	distance between fft data elements is 1024 */
		/*	do 2048 macros each processing 40 data values */

		loops_reset;			/* 4 more iters of 512 */
		r5_nf5cl_five_complex_djbfft_preload;
		for(unsigned int loopA = 4; loopA; loopA--) {
			for(unsigned int loopB = 4; loopB; loopB--) {
				for(unsigned int loopC = 128; loopC; loopC--) {
					r5_nf5cl_five_complex_djbfft(rsi, 64, 4*dist128, rdi, 512*5*XMM_SCD2);
					rdi += 5*XMM_SCD2;		/* Next sin/cos data */
				}
				rsi += -128*64+dist128;	/* Next source pointer */
			}
			rsi += -4*dist128+20*dist128; /* Next source pointer */
			rdi += -512*5*XMM_SCD2;	/* Restore sin/cos ptr */
		}
		rsi = daddr; //rsi += -5*20*dist128;	/* Restore source pointer */

		/* Do FFT levels 5,6 */
		/* */
		/* On input the 64-byte cache lines hold these data values: */
		/*	0	+50K	512	+50K	1024	+50K	1536	+50K */
		/*	1	... */
		/*	... */
		/*	511	... */
		/*	2048	... */
		/*	... */
		/* On output the 64-byte cache lines hold these data values: */
		/*	0	+50K	256	+50K	512	+50K	768	+50K */
		/*	1	... */
		/*	... */
		/*	255	... */
		/*	1024	... */
		/*	... */

		/* Do 256 eight_reals_fft operations */
		/* Do 256 four_complex_djbfft operations */
		/*	distance between fft data elements is 256 */
		/*	do 256 macros each processing 16 data values */

		loops_init_prefetch(6400, 64, 2, rcx);
		for(unsigned int loopA = 2; loopA; loopA--) {
			xtouch(xptr(rcx+4096-64));		/* Preload the TLB */
			for(unsigned int loopB = 128; loopB; loopB--) {
				r4_h2cl_eight_reals_four_complex_djbfft(rsi, 64, 2*dist128, rdi, rdx);
				rdi += 25*XMM_SCD2;	/* Next sin/cos data */
				rdx += XMM_HSCD2;		/* Next sin/cos data */
			}
			rsi += -128*64+dist128;	/* Next source pointer */
		}
		rsi += -2*dist128+4*dist128; /* Next source pointer */
		rdi += -256*25*XMM_SCD2;	/* Restore sin/cos ptr */
		rdx += -256*XMM_HSCD2;	/* Restore sin/cos ptr */

		/* Do 12288 four_complex_fft operations */
		/*	distance between fft data elements is 256 */
		/*	do 6144 macros each processing 16 data values */

		loops_reset;
		r4_x2cl_four_complex_djbfft_preload;
		for(unsigned int loopA = 24; loopA; loopA--) {
			for(unsigned int loopB = 2; loopB; loopB--) {
				xtouch(xptr(rcx+4096-64));		/* Preload the TLB */
				for(unsigned int loopC = 128; loopC; loopC--) {
					r4_x2cl_four_complex_djbfft(rsi, 64, 2*dist128, rdi);
					rdi += 25*XMM_SCD2;	/* Next sin/cos data */
				}
				rsi += -128*64+dist128;	/* Next source pointer */
			}
			rsi += -2*dist128+4*dist128; /* Next source pointer */
			rdi += -256*25*XMM_SCD2;	/* Restore sin/cos ptr */
		}
		rsi = daddr; //rsi += -25*4*dist128;	/* Restore source pointer */

		/* Do FFT levels 7,8 */
		/* */
		/* On input the 64-byte cache lines hold these data values: */
		/*	0	+50K	256	+50K	512	+50K	768	+50K */
		/*	1	... */
		/*	... */
		/*	255	... */
		/*	1024	... */
		/*	... */
		/* On output the 64-byte cache lines hold these data values: */
		/*	0	+50K	64	+50K	128	+50K	192	+50K */
		/*	1	... */
		/*	... */
		/*	63	... */
		/*	256	... */
		/*	... */

		/* Do 64 eight_reals_fft operations */
		/* Do 192 four_complex_djbfft operations */
		/*	distance between fft data elements is 64 */
		/*	do 64 macros each processing 32 data values */

		rdi += 6400*XMM_SCD2;	/* Load sin/cos pointer */
		loops_init_prefetch(3200, 64, 1, rcx);
		xtouch(xptr(rcx+4096-64));		/* Preload the TLB */
		r4_h4cl_eight_reals_four_complex_djbfft_preload;
		for(unsigned int loopA = 64; loopA; loopA--) {
			r4_h4cl_eight_reals_four_complex_djbfft(rsi, 64, 64*64, dist128, rdi, rdx);
			rdi += XMM_SCD2;		/* Next sin/cos data */
			rdx += 4*XMM_HSCD2;	/* Next sin/cos data */
		}
		rsi += -64*64+2*dist128;	/* Next source pointer */
		rdi += -64*XMM_SCD2;	/* Restore sin/cos ptr */
		rdx += -64*4*XMM_HSCD2;	/* Restore sin/cos ptr */

		/* Do 12544 four_complex_djbfft operations */
		/*	distance between fft data elements is 64 */
		/*	do 3136 macros each processing 32 data values */

		loops_reset;			/* 49 more iterations of 64 */
		r4_x4cl_four_complex_djbfft_preload;
		for(unsigned int loopA = 49; loopA; loopA--) {
			xtouch(xptr(rcx+4096-64));		/* Preload the TLB */
			for(unsigned int loopB = 64; loopB; loopB--) {
				r4_x4cl_four_complex_djbfft(rsi, 64, 64*64, dist128, rdi);
				rdi += XMM_SCD2;		/* Next sin/cos data */
			}
			rsi += -64*64+2*dist128;	/* Next source pointer */
			rdi += -64*XMM_SCD2;	/* Restore sin/cos ptr */
		}
		rsi = daddr; //rsi += -50*2*dist128;	/* Restore source pointer */

		/* Do FFT levels 9,10 */
		/* */
		/* On input the 64-byte cache lines hold these data values: */
		/*	0	+50K	64	+50K	128	+50K	192	+50K */
		/*	1	... */
		/*	... */
		/*	63	... */
		/*	256	... */
		/*	... */
		/* On output the 64-byte cache lines hold these data values: */
		/*	0	+50K	16	+50K	32	+50K	48	+50K */
		/*	1	... */
		/*	... */
		/*	15	... */
		/*	64	... */
		/*	... */

		/* Do 16 eight_reals_fft operations */
		/* Do 48 four_complex_djbfft operations */
		/*	distance between fft data elements is 16 */
		/*	do 16 macros each processing 32 data values */

		loops_init(3200);
		r4_h4cl_eight_reals_four_complex_djbfft_preload;
		for(unsigned int loopA = 16; loopA; loopA--) {
			r4_h4cl_eight_reals_four_complex_djbfft(rsi, 64, 16*64, 32*64, rdi, rdx);
			rdi += 4*XMM_SCD2;		/* Next sin/cos data */
			rdx += 16*XMM_HSCD2;	/* Next sin/cos data */
		}
		rsi += -16*64+64*64;	/* Next source pointer */
		rdi += -16*4*XMM_SCD2;	/* Restore sin/cos ptr */
		rdx += -16*16*XMM_HSCD2;	/* Restore sin/cos ptr */

		/* Do 12736 four_complex_djbfft operations */
		/*	distance between fft data elements is 16 */
		/*	do 3184 macros each processing 32 data values */

		loops_reset;			/* 100 iters of 1 or 2 of 16 */
		r4_x4cl_four_complex_djbfft_preload;
		for(unsigned int loopA = 100*2-1; loopA; ) {
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
		rsi = daddr; //rsi += -100*dist128;	/* Next source pointer */

		/* Do FFT levels 11,12 */
		/* */
		/* On input the 64-byte cache lines hold these data values: */
		/*	0	+50K	16	+50K	32	+50K	48	+50K */
		/*	1	... */
		/*	... */
		/*	15	... */
		/*	64	... */
		/*	... */
		/* On output the 64-byte cache lines hold these data values: */
		/*	0	+50K	4	+50K	8	+50K	12	+50K */
		/*	1	... */
		/*	... */
		/*	3	... */
		/*	16	... */
		/*	... */

		/* Do 4 eight_reals_fft operations */
		/* Do 12 four_complex_djbfft operations */
		/*	distance between fft data elements is 4 */
		/*	do 4 macros each processing 32 data values */

		loops_init(3200);
		r4_h4cl_eight_reals_four_complex_djbfft_preload;
		for(unsigned int loopA = 4; loopA; loopA--) {
			r4_h4cl_eight_reals_four_complex_djbfft(rsi, 64, 4*64, 8*64, rdi, rdx);
			rdi += 16*XMM_SCD2;	/* Next sin/cos data */
			rdx += 64*XMM_HSCD2;	/* Next sin/cos data */
		}
		rsi += -4*64+16*64;	/* Next source pointer */
		rdi += -4*16*XMM_SCD2;	/* Restore sin/cos ptr */
		rdx += -4*64*XMM_HSCD2;	/* Restore sin/cos ptr */

		/* Do 12784 four_complex_djbfft operations */
		/*	distance between fft data elements is 4 */
		/*	do 3196 macros each processing 32 data values */

		loops_reset;			/* 100 iters of 7 or 8 of 4 */
		r4_x4cl_four_complex_djbfft_preload;
		for(unsigned int loopA = 100*8-1; loopA; ) {
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
		rsi = daddr; //rsi += -100*dist128;	/* Restore source pointer */

		/* Do FFT levels 13,14 */
		/* */
		/* On input the 64-byte cache lines hold these data values: */
		/*	0	+50K	4	+50K	8	+50K	12	+50K */
		/*	1	... */
		/*	... */
		/*	3	... */
		/*	16	... */
		/*	... */
		/* On output the 64-byte cache lines hold these data values: */
		/*	0	+50K	2	+50K	4	+50K	6	+50K */
		/*	1	... */
		/*	8	... */
		/*	... */

		/* Do 1 eight_reals_fft operations */
		/* Do 3 four_complex_fft operations */
		/*	distance between fft data elements is 1 */

		/* Do 12796 four_complex_fft operations */
		/*	distance between fft data elements is 1 */
		/*	do 3199 macros each processing 32 data values */

		/* Execute the proper middle step */

		if likely(g->ffttype & 2) {
			if likely(!(g->ffttype & 1)) { // xpass2_25600_levels_real_2:;
				r4_h4cl_eight_reals_four_complex_with_square(rsi, 4*64, 64, 2*64);
				loops_init(3200, 1, 32);
				r4_x4cl_four_complex_with_square_preload;
				for(unsigned int loopA = 100*32-1; loopA; ) {
					do{
						r4_x4cl_four_complex_with_square(rsi, 4*64, 64, 2*64);
					}while(--loopA & 31);
					rsi += -32*4*64+dist128;	/* Next source pointer */
				}
				//rsi += -100*dist128;	/* Restore source pointer */
			}else{ // xpass2_25600_levels_real_3:;
				rbp = g->DIST_TO_MULSRCARG;
				r4_h4cl_eight_reals_four_complex_with_mult(rsi, 4*64, 64, 2*64);
				loops_init(3200, 1, 32);
				for(unsigned int loopA = 100*32-1; loopA; ) {
					do{
						r4_x4cl_four_complex_with_mult(rsi, 4*64, 64, 2*64);
					}while(--loopA & 31);
					rsi += -32*4*64+dist128;	/* Next source pointer */
				}
				//rsi += -100*dist128;	/* Restore source pointer */
			}
		}else{ // xpass2_25600_levels_real_1:;
			r4_h4cl_eight_reals_four_complex_fft_final(rsi, 4*64, 64, 2*64);
			L2prefetch128(u8ptr(rcx));
			pfing(rcx += 128);
			loops_init_prefetch(3200, 128, 1, rcx, 1, 32);
			for(unsigned int loopA = 100*32-1; loopA; ) {
				xtouch(xptr(rcx+4096-64));		/* Preload the TLBs */
				do{
					r4_x4cl_four_complex_fft_final(rsi, 4*64, 64, 2*64);
				}while(--loopA & 31);
				rsi += -32*4*64+dist128;	/* Next source pointer */
			}
			goto	xpass2_25600_real_done;
		}
	}else{ // xpass2_25600_levels_real_4:;
		rbp = g->DIST_TO_MULSRCARG;
		pfing(rcx = (uintptr_t)g->data_prefetch); /* Load prefetch pointer */
		xtouch(xptr(rcx));			/* Preload the TLBs */
		r4_h4cl_eight_reals_four_complex_with_mulf(rsi, 4*64, 64, 2*64);
		L2prefetch128(u8ptr(rcx));
		pfing(rcx += 128);
		loops_init_prefetch(3200, 128, 1, rcx, 1, 32);
		for(unsigned int loopA = 100*32-1; loopA; ) {
			xtouch(xptr(rcx+4096-64));		/* Preload the TLBs */
			do{
				r4_x4cl_four_complex_with_mulf(rsi, 4*64, 64, 2*64);
			}while(--loopA & 31);
			rsi += -32*4*64+dist128;	/* Next source pointer */
		}
		//rsi += -100*dist128;	/* Restore source pointer */
	}
	rsi = daddr;

	/* Do inverse FFT levels 11,12 */
	/* On input the 64-byte cache lines hold these data values: */
	/*	0	+50K	2	+50K	4	+50K	6	+50K */
	/*	1	... */
	/*	8	... */
	/*	... */
	/* On output the 64-byte cache lines hold these data values: */
	/*	0	+50K	8	+50K	16	+50K	24	+50K */
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
	rdi += 6400*XMM_SCD2;
	rdx = (uintptr_t)g->sincos3;
	rdx += 5120*XMM_HSCD2;
	loops_init(3200);
	for(unsigned int loopA = 2; loopA; loopA--) {
		r4_h4cl_eight_reals_four_complex_djbunfft(rsi, 64, 2*64, 4*64, rdi, 2*16*XMM_SCD2, rdx, 2*64*XMM_HSCD2);
		rdi += 16*XMM_SCD2;	/* Next sin/cos data */
		rdx += 64*XMM_HSCD2;	/* Next sin/cos data */
	}
	rsi += -2*64+8*64;		/* Next source pointer */
	rdi += -2*16*XMM_SCD2;	/* Restore sin/cos ptr */
	rdx += -2*64*XMM_HSCD2;	/* Restore sin/cos ptr */

	/* Do 12792 four_complex_djbunfft operations */
	/*	distance between fft data elements is 4 */
	/*	do 3198 macros each processing 32 data values */

	loops_reset;			/* 100 iters of 15 or 16 of 2 */
	r4_x4cl_four_complex_djbunfft_preload;
	for(unsigned int loopA = 100*16-1; loopA; ) {
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
	rsi = daddr; //rsi += -100*dist128;	/* Restore source pointer */

	/* Do inverse FFT levels 9,10 */
	/* On input the 64-byte cache lines hold these data values: */
	/*	0	+50K	8	+50K	16	+50K	24	+50K */
	/*	1	... */
	/*	... */
	/*	7	... */
	/*	32	... */
	/*	... */
	/* On output the 64-byte cache lines hold these data values: */
	/*	0	+50K	32	+50K	64	+50K	96	+50K */
	/*	1	... */
	/*	... */
	/*	31	... */
	/*	128	... */
	/*	... */

	/* Do 16 eight_reals_unfft operations */
	/* Do 16 four_complex_djbunfft operations */
	/*	distance between fft data elements is 16 */
	/*	do 8 macros processing 32 data values */

	loops_init_prefetch(3200, 64, 1, rcx);
	xtouch(xptr(rcx+4096-64));		/* Preload the TLBs */
	for(unsigned int loopA = 8; loopA; loopA--) {
		r4_h4cl_eight_reals_four_complex_djbunfft(rsi, 64, 8*64, 16*64, rdi, 8*4*XMM_SCD2, rdx, 8*16*XMM_HSCD2);
		rdi += 4*XMM_SCD2;		/* Next sin/cos data */
		rdx += 16*XMM_HSCD2;	/* Next sin/cos data */
	}
	rsi += -8*64+32*64;	/* Next source pointer */
	rdi += -8*4*XMM_SCD2;	/* Restore sin/cos ptr */
	rdx += -8*16*XMM_HSCD2;	/* Restore sin/cos ptr */

	/* Do 12768 four_complex_djbunfft operations */
	/*	distance between fft data elements is 16 */
	/*	do 3192 macros each processing 32 data values */

	loops_reset;			/* 100 iters of 3 or 4 of 4 */
	r4_x4cl_four_complex_djbunfft_preload;
	for(unsigned int loopA = 100*4-1; loopA; ) {
		xtouch(xptr(rcx+4096-64));		/* Preload the TLBs */
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
	rsi = daddr; //rsi += -100*dist128;	/* Restore source pointer */

	/* Do inverse FFT levels 7,8 */
	/* On input the 64-byte cache lines hold these data values: */
	/*	0	+50K	32	+50K	64	+50K	96	+50K */
	/*	1	... */
	/*	... */
	/*	31	... */
	/*	128	... */
	/*	... */
	/* On output the 64-byte cache lines hold these data values: */
	/*	0	+50K	128	+50K	256	+50K	384	+50K */
	/*	1	... */
	/*	... */
	/*	127	... */
	/*	512	... */
	/*	... */

	/* Do 64 eight_reals_unfft operations */
	/* Do 64 four_complex_djbunfft operations */
	/*	distance between fft data elements is 64 */
	/*	do 32 macros processing 32 data values */

	loops_init_prefetch(3200, 64, 1, rcx);
	xtouch(xptr(rcx+4096-64));		/* Preload the TLB */
	for(unsigned int loopA = 32; loopA; loopA--) {
		r4_h4cl_eight_reals_four_complex_djbunfft(rsi, 64, 32*64, 64*64, rdi, 32*XMM_SCD2, rdx, 32*4*XMM_HSCD2);
		rdi += XMM_SCD2;		/* Next sin/cos data */
		rdx += 4*XMM_HSCD2;	/* Next sin/cos data */
	}
	rsi += -32*64+dist128;	/* Next source pointer */
	rdi += -32*XMM_SCD2;	/* Restore sin/cos ptr */
	rdx += -32*4*XMM_HSCD2;	/* Restore sin/cos ptr */

	/* Do 12672 four_complex_djbunfft operations */
	/*	distance between fft data elements is 64 */
	/*	do 3168 macros each processing 32 data values */

	loops_reset;			/* 99 more iters of 32 */
	r4_x4cl_four_complex_djbunfft_preload;
	for(unsigned int loopA = 99; loopA; loopA--) {
		xtouch(xptr(rcx+4096-64));		/* Preload the TLB */
		for(unsigned int loopB = 32; loopB; loopB--) {
			r4_x4cl_four_complex_djbunfft(rsi, 64, 32*64, 64*64, rdi, 32*XMM_SCD2);
			rdi += XMM_SCD2;		/* Next sin/cos data */
		}
		rsi += -32*64+dist128;	/* Next source pointer */
		rdi += -32*XMM_SCD2;	/* Restore sin/cos ptr */
	}
	rsi = daddr; //rsi += -100*dist128;	/* Next source pointer */

	/* Do inverse FFT levels 5,6 */
	/* On input the 64-byte cache lines hold these data values: */
	/*	0	+50K	128	+50K	256	+50K	384	+50K */
	/*	1	... */
	/*	... */
	/*	127	... */
	/*	512	... */
	/*	... */
	/* On output the 64-byte cache lines hold these data values: */
	/*	0	+50K	512	+50K	1024	+50K	1536	+50K */
	/*	1	... */
	/*	... */
	/*	511	... */
	/*	2048	... */
	/*	... */

	/* Do 256 eight_reals_unfft operations */
	/* Do 256 four_complex_djbunfft operations */
	/*	distance between fft data elements is 256 */
	/*	do 128 macros processing 32 data values */

	loops_init(3200);
	rdi = (uintptr_t)g->xsincos_complex;	/* Load sin/cos pointer */
	for(unsigned int loopA = 128; loopA; loopA--) {
		r4_h4cl_eight_reals_four_complex_djbunfft(rsi, 64, dist128, 2*dist128, rdi, 128*25*XMM_SCD2, rdx, 128*XMM_HSCD2);
		rdi += 25*XMM_SCD2;	/* Next sin/cos data */
		rdx += XMM_HSCD2;		/* Next sin/cos data */
	}
	rsi += -128*64+4*dist128;	/* Next source pointer */
	rdi += -128*25*XMM_SCD2;	/* Restore sin/cos ptr */
	rdx += -128*XMM_HSCD2-5120*XMM_HSCD2; /* Next sin/cos ptr */

	/* Do 12288 four_complex_djbunfft operations */
	/*	distance between fft data elements is 256 */
	/*	do 3072 macros each processing 32 data values */

	loops_reset;			/* 24 more iters of 128 */
	r4_x4cl_four_complex_djbunfft_preload;
	for(unsigned int loopA = 24; loopA; loopA--) {
		for(unsigned int loopB = 128; loopB; loopB--) {
			r4_x4cl_four_complex_djbunfft(rsi, 64, dist128, 2*dist128, rdi, 128*25*XMM_SCD2);
			rdi += 25*XMM_SCD2;	/* Next sin/cos data */
		}
		rsi += -128*64+4*dist128;	/* Next source pointer */
		rdi += -128*25*XMM_SCD2;	/* Restore sin/cos ptr */
	}
	rsi = daddr; //rsi += -25*4*dist128;	/* Next source pointer */

	/* Do inverse FFT levels 3,4 (actually 2.322 FFT levels) */
	/* On input the 64-byte cache lines hold these data values: */
	/*	0	+50K	512	+50K	1024	+50K	1536	+50K */
	/*	1	... */
	/*	... */
	/*	511	... */
	/*	2048	... */
	/*	... */
	/* On output the 64-byte cache lines hold these data values: */
	/*	0	+50K	2560	+50K	5120	+50K	7680	+50K */
	/*	1	... */
	/*	... */
	/*	2559	... */
	/*	10240	... */
	/*	... */

	/* Do 1024 ten_reals_unfft operations */
	/* Do 1024 ten_complex_djbunfft operations */
	/*	distance between fft data elements is 1024 */
	/*	do 512 macros processing 40 data values */

	pfing(rcx += -2*64);		/* Back up prefetch pointer by 2 cache lines */
	loops_init_prefetch(2560, 64, 64, rcx);
	xtouch(xptr(rcx+40*64-64));		/* Preload the TLB */
	for(unsigned int loopA = 4; loopA; loopA--) {
		for(unsigned int loopB = 128; loopB; loopB--) {
			r5_h5cl_ten_reals_five_complex_djbunfft(rsi, 64, 4*dist128, rdi, 512*5*XMM_SCD2, rdx, 512*5*XMM_HSCD2);
			rdi += 5*XMM_SCD2;		/* Next sin/cos data */
			rdx += 5*XMM_HSCD2;	/* Next sin/cos data */
		}
		rsi += -128*64+dist128;	/* Next source pointer */
	}
	rsi += -4*dist128+20*dist128; /* Next source pointer */
	rdi += -512*5*XMM_SCD2;	/* Restore sin/cos ptr */
	rdx += -512*5*XMM_HSCD2;	/* Next sin/cos ptr */

	/* Do 8192 five_complex_djbunfft operations */
	/*	distance between fft data elements is 1024 */
	/*	do 2048 macros each processing 40 data values */

	loops_reset;			/* 4 more iters of 512 */
	r5_x5cl_five_complex_djbunfft_preload;
	for(unsigned int loopA = 4; loopA; loopA--) {
		for(unsigned int loopB = 4; loopB; loopB--) {
			for(unsigned int loopC = 128; loopC; loopC--) {
				r5_x5cl_five_complex_djbunfft(rsi, 64, 4*dist128, rdi, 512*5*XMM_SCD2);
				rdi += 5*XMM_SCD2;		/* Next sin/cos data */
			}
			rsi += -128*64+dist128;	/* Next source pointer */
		}
		rsi += -4*dist128+20*dist128; /* Next source pointer */
		rdi += -512*5*XMM_SCD2;	/* Restore sin/cos ptr */
	}
	rsi = daddr; //rsi += -5*20*dist128;	/* Restore source pointer */

	/* Do inverse FFT levels 1,2 (actually 2.322 FFT levels) */
	/* On input the 64-byte cache lines hold these data values: */
	/*	0	+50K	2560	+50K	5120	+50K	7680	+50K */
	/*	1	... */
	/*	... */
	/*	2559	... */
	/*	10240	... */
	/*	... */
	/* On output the 64-byte cache lines hold these data values if pass 1 swizzles: */
	/*	0	+50K	12800	+50K	25600	+50K	38400	+50K */
	/*	1	... */
	/*	... */
	/*	12799	... */
	/* On output the 64-byte cache lines hold these data values if pass 2 swizzles: */
	/*	0	+12800	25600	+12800	51200	+12800	76800	+12800 */
	/*	1	... */
	/*	... */
	/*	12799	... */

	/* Do 5120 ten_reals_unfft operations */
	/* Do 5120 five_complex_djbunfft operations */
	/*	distance between fft data elements is 5120 */
	/*	do 2560 macros processing 40 data values */

	loops_init(2560);
	for(unsigned int loopA = 20; loopA; loopA--) {
		for(unsigned int loopB = 128; loopB; loopB--) {
			r5_h5cl_ten_reals_five_complex_djbunfft(rsi, 64, 20*dist128, rdi, 2560*XMM_SCD2, rdx, 2560*XMM_HSCD2);
			rdi += XMM_SCD2;		/* Next sin/cos data */
			rdx += XMM_HSCD2;		/* Next sin/cos data */
		}
		rsi += -128*64+dist128;	/* Next source pointer */
	}

	xpass2_25600_real_done:;
	end_timer(2);
	if(g->pass2_get_next_block(g) != 0) return;
}
r4_pass2_25600_levels_complex(g);
}


	/* *************** radix4-pass2-25600-levels-complex macro ****************** */
	/* This macro takes 2 groups of 25600 complex values and */
	/* performs the final 14.644 levels of a standard radix-4 FFT, squares the */
	/* results, and does the first 14.644 levels of the inverse FFT. */


void r4_pass2_25600_levels_complex(struct gwasm_data *__restrict g) {
do {
	uintptr_t rbx,rbp,rsi,rdi,daddr;
	vec2f64 xmm0,xmm1,xmm2,xmm3,xmm4,xmm5,xmm6,xmm7;

	/* Type 4 FFTs skip the forward FFT process */

	daddr = (uintptr_t)g->data_addr;/* Load address of FFT data */
	rsi = daddr;
	rbx = g->DIST_TO_FFTSRCARG;
	if likely(!(g->ffttype & 4)) { //xpass2_25600_levels_complex_4;

		/* Do FFT level 1,2 (actually 2.322 levels) */
		/* */
		/* We swizzle the data so that one set of data is in the high word of */
		/* an xmm register, the other set is in the low word.  This means */
		/* a distance of 12800 data values. */
		/* */
		/* On input the 64-byte cache lines hold these data values if pass 2 swizzles: */
		/*	0	+3.2K	6400	+3.2K	12800	+3.2K	19200	+3.2K */
		/*	1	... */
		/*	... */
		/*	3199 */
		/* On input the 64-byte cache lines hold these data values if pass 1 swizzles: */
		/*	0	+12.5K	3200	+12.5K	6400	+12.5K	9600	+12.5K */
		/*	1	... */
		/*	... */
		/*	3199	... */
		/* On output the 64-byte cache lines hold these data values: */
		/*	0	+12.5K	640	+12.5K	1280	+12.5K	1920	+12.5K */
		/*	1	... */
		/*	... */
		/*	639	... */
		/*	2560	... */
		/*	... */

		/* Do 10240 five_complex_fft operations */
		/*	distance between fft data elements is 5120 */
		/*	do 2560 macros each processing 40 data values */

		start_timer(3);
		rdi = (uintptr_t)g->xsincos_complex;	/* The sin/cos data */
		loops_init(2560);			/* 20 iterations of 128 */
		r5_f5cl_five_complex_djbfft_preload;
		for(unsigned int loopA = 20; loopA; loopA--) {
			for(unsigned int loopB = 128; loopB; loopB--) {
				xxprefetch(u8ptr(rdi+XMM_SCD2));	/* Prefetch for next pass */
				xxprefetch(u8ptr(rdi+2561*XMM_SCD2));	/* Prefetch for next pass */
				r5_f5cl_five_complex_djbfft(rsi, 64, 20*dist128, rdi, 2560*XMM_SCD2);
				rdi += XMM_SCD2;		/* Next sin/cos data */
			}
			rsi += -128*64+dist128;	/* Next source pointer */
		}
		rsi = daddr; //rsi += -20*dist128;	/* Restore source pointer */
		end_timer(3);

		/* Do FFT levels 3,4 (actually 2.322 levels) */
		/* */
		/* On input the 64-byte cache lines hold these data values: */
		/*	0	+50K	2560	+50K	5120	+50K	7680	+50K */
		/*	1	... */
		/*	... */
		/*	2559	... */
		/*	10240	... */
		/*	... */
		/* On output the 64-byte cache lines hold these data values: */
		/*	0	+50K	512	+50K	1024	+50K	1536	+50K */
		/*	1	... */
		/*	... */
		/*	511	... */
		/*	2048	... */
		/*	... */

		/* Do 10240 five_complex_fft operations */
		/*	distance between fft data elements is 1024 */
		/*	do 2560 macros each processing 40 data values */

		start_timer(4);
		pfing(rcx = (uintptr_t)g->data_prefetch); /*  Prefetch pointer (we need to prefetch  */
		;					/* 12800 cache lines + 198 pad cache lines) */
		xtouch(xptr(rcx));			/* Load first prefetch TLB */
		rdi = (uintptr_t)g->xsincos_complex;	/* Ptr to the sin/cos data */
		loops_init_prefetch(2560, 64, 16, rcx); /* 5 iters of 4 iters of 128 */
		xtouch(xptr(rcx+4096-64));		/* Load prefetch TLBs */
		xtouch(xptr(rcx+2*4096-64));
		xtouch(xptr(rcx+160*64-64));
		r5_nf5cl_five_complex_djbfft_preload;
		for(unsigned int loopA = 5; loopA; loopA--) {
			for(unsigned int loopB = 4; loopB; loopB--) {
				for(unsigned int loopC = 128; loopC; loopC--) {
					r5_nf5cl_five_complex_djbfft(rsi, 64, 4*dist128, rdi, 512*5*XMM_SCD2);
					rdi += 5*XMM_SCD2;		/* Next sin/cos data */
				}
				rsi += -128*64+dist128;	/* Next source pointer */
			}
			rsi += -4*dist128+20*dist128; /* Next source pointer */
			rdi += -512*5*XMM_SCD2;	/* Restore sin/cos ptr */
		}
		rsi = daddr; //rsi += -5*20*dist128;	/* Restore source pointer */
		end_timer(4);

		/* Do FFT levels 5,6 */
		/* */
		/* On input the 64-byte cache lines hold these data values: */
		/*	0	+50K	512	+50K	1024	+50K	1536	+50K */
		/*	1	... */
		/*	... */
		/*	511	... */
		/*	2048	... */
		/*	... */
		/* On output the 64-byte cache lines hold these data values: */
		/*	0	+50K	256	+50K	512	+50K	768	+50K */
		/*	1	... */
		/*	... */
		/*	255	... */
		/*	1024	... */
		/*	... */

		/* Do 12800 four_complex_fft operations */
		/*	distance between fft data elements is 256 */
		/*	do 6400 macros each processing 16 data values */

		start_timer(5);
		loops_init_prefetch(6400, 64, 4, rcx); /* 25 iters of 2 iters of 128 */
		r4_x2cl_four_complex_djbfft_preload;
		for(unsigned int loopA = 25; loopA; loopA--) {
			xtouch(xptr(rcx+4096-64));		/* Load prefetch TLB */
			for(unsigned int loopB = 2; loopB; loopB--) {
				for(unsigned int loopC = 128; loopC; loopC--) {
					r4_x2cl_four_complex_djbfft(rsi, 64, 2*dist128, rdi);
					rdi += 25*XMM_SCD2;	/* Next sin/cos data */
				}
				rsi += -128*64+dist128;	/* Next source pointer */
			}
			rsi += -2*dist128+4*dist128; /* Next source pointer */
			rdi += -256*25*XMM_SCD2;	/* Restore sin/cos ptr */
		}
		rsi = daddr; //rsi += -25*4*dist128;	/* Restore source pointer */
		end_timer(5);

		/* Do FFT levels 7,8 */
		/* */
		/* On input the 64-byte cache lines hold these data values: */
		/*	0	+50K	256	+50K	512	+50K	768	+50K */
		/*	1	... */
		/*	... */
		/*	255	... */
		/*	1024	... */
		/*	... */
		/* On output the 64-byte cache lines hold these data values: */
		/*	0	+50K	64	+50K	128	+50K	192	+50K */
		/*	1	... */
		/*	... */
		/*	63	... */
		/*	256	... */
		/*	... */
		/* */
		/* Do 12800 four_complex_fft operations */
		/*	distance between fft data elements is 64 */
		/*	do 3200 macros each processing 32 data values */

		/* Do FFT levels 9,10 */
		/* */
		/* On input the 64-byte cache lines hold these data values: */
		/*	0	+50K	64	+50K	128	+50K	192	+50K */
		/*	1	... */
		/*	... */
		/*	63	... */
		/*	256	... */
		/*	... */
		/* On output the 64-byte cache lines hold these data values: */
		/*	0	+50K	16	+50K	32	+50K	48	+50K */
		/*	1	... */
		/*	... */
		/*	15	... */
		/*	64	... */
		/*	... */
		/* */
		/* Do 12800 four_complex_fft operations */
		/*	distance between fft data elements is 16 */
		/*	do 3200 macros each processing 32 data values */

		/* Do FFT levels 11,12 */
		/* */
		/* On input the 64-byte cache lines hold these data values: */
		/*	0	+50K	16	+50K	32	+50K	48	+50K */
		/*	1	... */
		/*	... */
		/*	15	... */
		/*	64	... */
		/*	... */
		/* On output the 64-byte cache lines hold these data values: */
		/*	0	+50K	4	+50K	8	+50K	12	+50K */
		/*	1	... */
		/*	... */
		/*	3	... */
		/*	16	... */
		/*	... */
		/* */
		/* Do 12800 four_complex_fft operations */
		/*	distance between fft data elements is 4 */
		/*	do 3200 macros each processing 32 data values */

		/* Do this 64 cache lines at a time to maximize accessing data */
		/* while it is in the L1 cache. */

		start_timer(6);
		rdi += 6400*XMM_SCD2;	/* Load sin/cos pointers */
		rbx = rdi;
		rbp = rdi;
		loops_init_prefetch(3200, 64, 2, rcx); /* 2*4 iters of 4 of 64 cache lines */
		r4_x4cl_four_complex_djbfft_preload;

		for(unsigned int loopA = 50; loopA; loopA--) {
			for(unsigned int loopB = 4; loopB; loopB--) {
				xtouch(xptr(rcx+48*64-64));		/* Load prefetch TLB */

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
		rsi = daddr; //rsi += -50*2*dist128;	/* Restore source pointer */
		end_timer(6);

		/* Do FFT levels 13,14 as well as inverse FFT levels 13,14 */
		/* */
		/* On input the 64-byte cache lines hold these data values: */
		/*	0	+50K	4	+50K	8	+50K	12	+50K */
		/*	1	... */
		/*	... */
		/*	3	... */
		/*	16	... */
		/*	... */
		/* On output the 64-byte cache lines hold these data values: */
		/*	0	+50K	2	+50K	4	+50K	6	+50K */
		/*	1	... */
		/*	8	... */
		/*	... */

		/* Do 12800 four_complex_fft operations */
		/*	distance between fft data elements is 1 */
		/*	do 3200 macros each processing 32 data values */

		/* Execute the right middle step */

		if likely(g->ffttype & 2) {
			if likely(!(g->ffttype & 1)) { // xpass2_25600_levels_complex_2:;
				start_timer(9);
				loops_init(3200);
				r4_x4cl_four_complex_with_square_preload;
				for(unsigned int loopA = 100; loopA; loopA--) {
					for(unsigned int loopB = 32; loopB; loopB--) {
						r4_x4cl_four_complex_with_square(rsi, 4*64, 64, 2*64);
					}
					rsi += -32*4*64+dist128;	/* Next source pointer */
				}
				//rsi += -100*dist128;	/* Restore source pointer */
				end_timer(9);
			}else{ // xpass2_25600_levels_complex_3:;
				start_timer(9);
				rbp = g->DIST_TO_MULSRCARG;
				loops_init(3200);
				for(unsigned int loopA = 100; loopA; loopA--) {
					for(unsigned int loopB = 32; loopB; loopB--) {
						r4_x4cl_four_complex_with_mult(rsi, 4*64, 64, 2*64);
					}
					rsi += -32*4*64+dist128;	/* Next source pointer */
				}
				//rsi += -100*dist128;	/* Restore source pointer */
				end_timer(9);
			}
		}else{ // xpass2_25600_levels_complex_1:;
			start_timer(9);
			loops_init_prefetch(3200, 128, 1, rcx);
			for(unsigned int loopA = 100; loopA; loopA--) {
				xtouch(xptr(rcx+4096-128));		/* Load prefetch TLBs */
				for(unsigned int loopB = 32; loopB; loopB--) {
					r4_x4cl_four_complex_fft_final(rsi, 4*64, 64, 2*64);
				}
				rsi += -32*4*64+dist128;	/* Next source pointer */
			}
			end_timer(9);
			goto	xpass2_25600_levels_complex_done;
		}
	}else{ // xpass2_25600_levels_complex_4:;
		start_timer(9);
		pfing(rcx = (uintptr_t)g->data_prefetch); /* Prefetch pointer */
		xtouch(xptr(rcx));			/* Load first prefetch TLB */
		rbp = g->DIST_TO_MULSRCARG;
		loops_init_prefetch(3200, 128, 1, rcx);
		for(unsigned int loopA = 100; loopA; loopA--) {
			xtouch(xptr(rcx+4096-128));		/* Load prefetch TLBs */
			for(unsigned int loopB = 32; loopB; loopB--) {
				r4_x4cl_four_complex_with_mulf(rsi, 4*64, 64, 2*64);
			}
			rsi += -32*4*64+dist128;	/* Next source pointer */
		}
		//rsi += -100*dist128;	/* Restore source pointer */
		end_timer(9);
	}
	rsi = daddr;

	/* Do inverse FFT levels 11,12 */
	/* On input the 64-byte cache lines hold these data values: */
	/*	0	+50K	2	+50K	4	+50K	6	+50K */
	/*	1	... */
	/*	8	... */
	/*	... */
	/* On output the 64-byte cache lines hold these data values: */
	/*	0	+50K	8	+50K	16	+50K	24	+50K */
	/*	1	... */
	/*	... */
	/*	7	... */
	/*	32	... */
	/*	... */
	/* */
	/* Do 12800 four_complex_unfft operations */
	/*	distance between fft data elements is 4 */
	/*	do 3200 macros each processing 32 data values */

	/* Do inverse FFT levels 9,10 */
	/* On input the 64-byte cache lines hold these data values: */
	/*	0	+50K	8	+50K	16	+50K	24	+50K */
	/*	1	... */
	/*	... */
	/*	7	... */
	/*	32	... */
	/*	... */
	/* On output the 64-byte cache lines hold these data values: */
	/*	0	+50K	32	+50K	64	+50K	96	+50K */
	/*	1	... */
	/*	... */
	/*	31	... */
	/*	128	... */
	/*	... */
	/* */
	/* Do 12800 four_complex_unfft operations */
	/*	distance between fft data elements is 16 */
	/*	do 3200 macros each processing 32 data values */

	/* Do inverse FFT levels 7,8 */
	/* On input the 64-byte cache lines hold these data values: */
	/*	0	+50K	32	+50K	64	+50K	96	+50K */
	/*	1	... */
	/*	... */
	/*	31	... */
	/*	128	... */
	/*	... */
	/* On output the 64-byte cache lines hold these data values: */
	/*	0	+50K	128	+50K	256	+50K	384	+50K */
	/*	1	... */
	/*	... */
	/*	127	... */
	/*	512	... */
	/*	... */
	/* */
	/* Do 12800 four_complex_unfft operations */
	/*	distance between fft data elements is 64 */
	/*	do 3200 macros each processing 32 data values */

	/* Do this 64 cache lines at a time to maximize accessing data */
	/* while it is in the L1 cache. */

	start_timer(10);
	rdi = (uintptr_t)g->xsincos_complex;	/* Load sin/cos pointers */
	rdi += 6400*XMM_SCD2;
	rbx = rdi;
	rbp = rdi;
	loops_init_prefetch(3200, 64, 2, rcx); /* 8*2*8 iters of 64 cache lines */
	r4_x4cl_four_complex_djbunfft_preload;

	for(unsigned int loopA = 100; loopA; loopA--) {
		xtouch(xptr(rcx+48*64-64));		/* Load prefetch TLB */

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
	rsi = daddr; //rsi += -100*dist128;	/* Restore source pointer */
	end_timer(10);

	/* Do inverse FFT levels 5,6 */
	/* On input the 64-byte cache lines hold these data values: */
	/*	0	+50K	128	+50K	256	+50K	384	+50K */
	/*	1	... */
	/*	... */
	/*	127	... */
	/*	512	... */
	/*	... */
	/* On output the 64-byte cache lines hold these data values: */
	/*	0	+50K	512	+50K	1024	+50K	1536	+50K */
	/*	1	... */
	/*	... */
	/*	511	... */
	/*	2048	... */
	/*	... */

	/* Do 12800 four_complex_unfft operations */
	/*	distance between fft data elements is 256 */
	/*	do 3200 macros each processing 32 data values */

	start_timer(13);
	rdi = (uintptr_t)g->xsincos_complex;	/* Load sin/cos pointer */
	loops_init_prefetch(3200, 64, 2, rcx);/* 25 iters of 128 */
	r4_x4cl_four_complex_djbunfft_preload;
	for(unsigned int loopA = 25; loopA; loopA--) {
		xtouch(xptr(rcx+4096-64));		/* Load prefetch TLB */
		for(unsigned int loopB = 128; loopB; loopB--) {
			r4_x4cl_four_complex_djbunfft(rsi, 64, dist128, 2*dist128, rdi, 128*25*XMM_SCD2);
			rdi += 25*XMM_SCD2;	/* Next sine/cosine pointer */
		}
		rsi += -128*64+4*dist128;	/* Next source pointer */
		rdi += -128*25*XMM_SCD2;	/* Reset sine/cosine pointer */
	}
	rsi = daddr; //rsi += -25*4*dist128;	/* Restore source pointer */
	end_timer(13);

	/* Do inverse FFT levels 3,4 (actually 2.322 FFT levels) */
	/* On input the 64-byte cache lines hold these data values: */
	/*	0	+50K	512	+50K	1024	+50K	1536	+50K */
	/*	1	... */
	/*	... */
	/*	511	... */
	/*	2048	... */
	/*	... */
	/* On output the 64-byte cache lines hold these data values: */
	/*	0	+50K	2560	+50K	5120	+50K	7680	+50K */
	/*	1	... */
	/*	... */
	/*	2559	... */
	/*	10240	... */
	/*	... */

	/* Do 10240 five_complex_unfft operations */
	/*	distance between fft data elements is 1024 */
	/*	do 2560 macros each processing 40 data values */

	start_timer(14);
	pfing(rcx += -2*64);		/* Back up prefetch pointer by 2 cache lines */
	loops_init_prefetch(2560, 64, 64, rcx); /* 5 iters of 4 iters of 128 */
	xtouch(xptr(rcx+40*64-64));		/* Load prefetch TLB */
	r5_x5cl_five_complex_djbunfft_preload;
	for(unsigned int loopA = 5; loopA; loopA--) {
		for(unsigned int loopB = 4; loopB; loopB--) {
			for(unsigned int loopC = 128; loopC; loopC--) {
				r5_x5cl_five_complex_djbunfft(rsi, 64, 4*dist128, rdi, 512*5*XMM_SCD2);
				rdi += 5*XMM_SCD2;		/* Next sine/cosine pointer */
			}
			rsi += -128*64+dist128;	/* Next source pointer */
		}
		rsi += -4*dist128+20*dist128; /* Next source pointer */
		rdi += -512*5*XMM_SCD2;	/* Reset sine/cosine pointer */
	}
	rsi = daddr; //rsi += -5*20*dist128;	/* Restore source pointer */
	end_timer(14);

	/* Do inverse FFT levels 1,2 (actually 2.322 FFT levels) */
	/* On input the 64-byte cache lines hold these data values: */
	/*	0	+50K	2560	+50K	5120	+50K	7680	+50K */
	/*	1	... */
	/*	... */
	/*	2559	... */
	/*	10240	... */
	/*	... */
	/* On output the 64-byte cache lines hold these data values if pass 1 swizzles: */
	/*	0	+50K	12800	+50K	25600	+50K	38400	+50K */
	/*	1	... */
	/*	... */
	/*	12799	... */
	/* On output the 64-byte cache lines hold these data values if pass 2 swizzles: */
	/*	0	+12800	25600	+12800	51200	+12800	76800	+12800 */
	/*	1	... */
	/*	... */
	/*	12799	... */

	/* Do 10240 five_complex_unfft operations */
	/*	distance between fft data elements is 5120 */
	/*	do 2560 macros each processing 40 data values */

	start_timer(15);
	loops_init(2560);			/* 20 iters of 128 */
	r5_x5cl_five_complex_djbunfft_preload;
	for(unsigned int loopA = 20; loopA; loopA--) {
		for(unsigned int loopB = 128; loopB; loopB--) {
			r5_x5cl_five_complex_djbunfft(rsi, 64, 20*dist128, rdi, 2560*XMM_SCD2);
			rdi += XMM_SCD2;		/* Next set of sin/cos data */
		}
		rsi += -128*64+dist128;	/* Next source pointer */
	}
	end_timer(15);

	xpass2_25600_levels_complex_done:;
}while(g->pass2_get_next_block(g) == 0);
}
