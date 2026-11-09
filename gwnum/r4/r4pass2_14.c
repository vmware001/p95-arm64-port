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

/* *************** radix4-pass2-14-levels-real macro ****************** */
/* This macro takes 32768 real values and performs the final 15 levels of */
/* the FFT.  Independently, this macro also takes 16384 complex values and */
/* performs the final 14 levels of the FFT.  The results are squared, and */
/* the first 14/15 levels of the inverse FFT are done. */


void r4_pass2_14_levels_complex(struct gwasm_data *__restrict g);
void xpass2_r4_14_levels(struct gwasm_data *__restrict g) {
	g->thread_work_routine = r4_pass2_14_levels_complex;
	pass2_wake_up_threads(g);
	if(g->NEGACYCLIC_FFT == 0) {
//void r4_pass2_14_levels_real(struct gwasm_data *__restrict g)
	uintptr_t rdx,rbx,rbp,rsi,rdi,daddr;
	vec2f64 xmm0,xmm1,xmm2,xmm3,xmm4,xmm5,xmm6,xmm7;

	/* We switch to the same format used in one-pass FFTs (that is, the */
	/* high word of xmm registers being 1 greater than the low word): */

	/* Type 4 FFTs skip the forward FFT process */

	start_timer(2);
	daddr = (uintptr_t)g->data_addr;/* Load source address */
	rsi = daddr;
	rbx = g->DIST_TO_FFTSRCARG;
	if likely(!(g->ffttype & 4)) { //xpass2_14_levels_real_4;

		/* Do FFT level 1,2 */
		/* */
		/* On input the 64-byte cache lines hold these data values if pass 2 swizzles: */
		/*	0	+8K	16384	+8K	32768	+8K	49152	+8K */
		/*	1	... */
		/*	... */
		/*	8191	... */
		/* On input the 64-byte cache lines hold these data values if pass 1 swizzles: */
		/*	0	+32K	8192	+32K	16384	+32K	24576	+32K */
		/*	1	... */
		/*	... */
		/*	8191	... */
		/* On output the 64-byte cache lines hold these data values: */
		/*	0	+32K	4096	+32K	8192	+32K	12288	+32K */
		/*	1	... */
		/*	... */
		/*	4095	... */
		/*	16384	... */
		/*	... */

		/* Do 4096 eight_reals_fft operations */
		/* Do 4096 four_complex_djbfft operations */
		/*	distance between fft data elements is 4K */
		/*	do 4096 iterations */

		rdi = (uintptr_t)g->xsincos_complex;	/* Ptr to the sin/cos data */
		rdx = (uintptr_t)g->sincos3;/* Ptr to real sin/cos data */
		loops_init(4096);
		for(unsigned int loopA = 32; loopA; loopA--) {
			for(unsigned int loopB = 128; loopB; loopB--) {
				r4_fh2cl_eight_reals_four_complex_djbfft(rsi, 64, 32*dist128, rdi, rdx);
				rdi += XMM_SCD2;		/* Next sin/cos data */
				rdx += XMM_HSCD2;		/* Next sin/cos data */
			}
			rsi += -128*64+dist128;	/* Next source pointer */
		}
		rsi = daddr; //rsi += -32*dist128;	/* Restore source pointer */
		rdi += -4096*XMM_SCD2;	/* Restore sin/cos ptr */
		rdx += -4096*XMM_HSCD2;	/* Restore sin/cos ptr */

		/* Do FFT levels 3,4 */
		/* */
		/* On input the 64-byte cache lines hold these data values: */
		/*	0	+32K	4096	+32K	8192	+32K	12288	+32K */
		/*	1	... */
		/*	... */
		/*	4095	... */
		/*	16384	... */
		/*	... */
		/* On output the 64-byte cache lines hold these data values: */
		/*	0	+32K	1024	+32K	2048	+32K	3072	+32K */
		/*	1	... */
		/*	... */
		/*	1023	... */
		/*	4096	... */
		/*	... */

		/* Do 1024 eight_reals_fft operations */
		/* Do 3072 four_complex_djbfft operations */
		/*	distance between fft data elements is 1024 */
		/*	do 1024 macros each processing 32 data values */

		pfing(rcx = (uintptr_t)g->data_prefetch); /* Prefetch pointer */
		xtouch(xptr(rcx));			/* Preload first TLB */
		loops_init(1024);
		r4_h4cl_eight_reals_four_complex_djbfft_preload;
		for(unsigned int loopA = 8; loopA; loopA--) {
			for(unsigned int loopB = 128; loopB; loopB--) {
				r4_h4cl_eight_reals_four_complex_djbfft(rsi, 64, 8*dist128, 16*dist128, rdi, rdx);
				rdi += 4*XMM_SCD2;		/* Next sin/cos data */
				rdx += 4*XMM_HSCD2;	/* Next sin/cos data */
			}
			rsi += -128*64+dist128;	/* Next source pointer */
		}
		rsi += -8*dist128+32*dist128; /* Next source pointer */
		rdi += -1024*4*XMM_SCD2;	/* Restore sin/cos ptr */
		rdx += -1024*4*XMM_HSCD2;	/* Restore sin/cos ptr */

		/* Do 4096 four_complex_fft operations */
		/*	distance between fft data elements is 1024 */
		/*	do 1024 macros each processing 32 data values */

		loops_init(1024);
		r4_x4cl_four_complex_djbfft_preload;
		for(unsigned int loopA = 8; loopA; loopA--) {
			for(unsigned int loopB = 128; loopB; loopB--) {
				r4_x4cl_four_complex_djbfft(rsi, 64, 8*dist128, 16*dist128, rdi);
				rdi += 4*XMM_SCD2;		/* Next sin/cos data */
			}
			rsi += -128*64+dist128;	/* Next source pointer */
		}
		rsi = daddr; //rsi += -8*dist128-32*dist128; /* Restore source pointer */
		rdi += -1024*4*XMM_SCD2;	/* Restore sin/cos pointer */

		/* Do FFT levels 5,6 */
		/* */
		/* On input the 64-byte cache lines hold these data values: */
		/*	0	+32K	1024	+32K	2048	+32K	3072	+32K */
		/*	1	... */
		/*	... */
		/*	1023	... */
		/*	4096	... */
		/*	... */
		/* On output the 64-byte cache lines hold these data values: */
		/*	0	+32K	256	+32K	512	+32K	768	+32K */
		/*	1	... */
		/*	... */
		/*	255	... */
		/*	1024	... */
		/*	... */

		/* Do 256 eight_reals_fft operations */
		/* Do 768 four_complex_djbfft operations */
		/*	distance between fft data elements is 256 */
		/*	do 256 macros each processing 32 data values */

		loops_init_prefetch(2048, 64, 1, rcx);
		r4_h4cl_eight_reals_four_complex_djbfft_preload;
		for(unsigned int loopA = 2; loopA; loopA--) {
			for(unsigned int loopB = 128; loopB; loopB--) {
				xtouch(xptr(rcx+4096-64));		/* Preload the TLBs */
				r4_h4cl_eight_reals_four_complex_djbfft(rsi, 64, 2*dist128, 4*dist128, rdi, rdx);
				rdi += 16*XMM_SCD2;	/* Next sin/cos data */
				rdx += 16*XMM_HSCD2;	/* Next sin/cos data */
			}
			rsi += -128*64+dist128;	/* Next source pointer */
			pfing(rcx += 128);		/* Skip pad bytes */
		}
		rsi += -2*dist128+8*dist128; /* Next source pointer */
		rdi += -256*16*XMM_SCD2;	/* Restore sin/cos ptr */
		rdx += -256*16*XMM_HSCD2;	/* Restore sin/cos ptr */

		/* Do 7168 four_complex_fft operations */
		/*	distance between fft data elements is 256 */
		/*	do 1792 macros each processing 32 data values */

		loops_reset;
		r4_x4cl_four_complex_djbfft_preload;
		for(unsigned int loopA = 7; loopA; loopA--) {
			for(unsigned int loopB = 2; loopB; loopB--) {
				for(unsigned int loopC = 128; loopC; loopC--) {
					xtouch(xptr(rcx+4096-64));		/* Preload the TLBs */
					r4_x4cl_four_complex_djbfft(rsi, 64, 2*dist128, 4*dist128, rdi);
					rdi += 16*XMM_SCD2;	/* Next sin/cos data */
				}
				rsi += -128*64+dist128;	/* Next source pointer */
				pfing(rcx += 128);		/* Skip pad bytes */
			}
			rsi += -2*dist128+8*dist128; /* Next source pointer */
			rdi += -256*16*XMM_SCD2;	/* Restore sin/cos ptr */
		}
		rsi = daddr; //rsi += -8*8*dist128;	/* Restore source pointer */

		/* Do FFT levels 7,8 */
		/* */
		/* On input the 64-byte cache lines hold these data values: */
		/*	0	+32K	256	+32K	512	+32K	768	+32K */
		/*	1	... */
		/*	... */
		/*	255	... */
		/*	1024	... */
		/*	... */
		/* On output the 64-byte cache lines hold these data values: */
		/*	0	+32K	64	+32K	128	+32K	192	+32K */
		/*	1	... */
		/*	... */
		/*	63	... */
		/*	256	... */
		/*	... */

		/* Do 64 eight_reals_fft operations */
		/* Do 192 four_complex_djbfft operations */
		/*	distance between fft data elements is 64 */
		/*	do 64 macros each processing 32 data values */

		rdi += 4096*XMM_SCD2;	/* Load sin/cos pointer */
		loops_init_prefetch(2048, 64, 1, rcx);
		xtouch(xptr(rcx+4096-64));		/* Preload the TLBs */
		r4_h4cl_eight_reals_four_complex_djbfft_preload;
		for(unsigned int loopA = 64; loopA; loopA--) {
			r4_h4cl_eight_reals_four_complex_djbfft(rsi, 64, 64*64, dist128, rdi, rdx);
			rdi += XMM_SCD2;		/* Next sin/cos data */
			rdx += 64*XMM_HSCD2;	/* Next sin/cos data */
		}
		rsi += -64*64+2*dist128;	/* Next source pointer */
		rdi += -64*XMM_SCD2;	/* Restore sin/cos ptr */
		rdx += -64*64*XMM_HSCD2;	/* Restore sin/cos ptr */

		/* Do 7944 four_complex_djbfft operations */
		/*	distance between fft data elements is 64 */
		/*	do 1986 macros each processing 32 data values */

		loops_reset;			/* 31 iterations of 64 */
		r4_x4cl_four_complex_djbfft_preload;
		for(unsigned int loopA = 16*2-1; loopA; ) {
			do{
				xtouch(xptr(rcx+4096-64));		/* Preload the TLBs */
				for(unsigned int loopC = 64; loopC; loopC--) {
					r4_x4cl_four_complex_djbfft(rsi, 64, 64*64, dist128, rdi);
					rdi += XMM_SCD2;		/* Next sin/cos data */
				}
				rsi += -64*64+2*dist128;	/* Next source pointer */
				rdi += -64*XMM_SCD2;	/* Restore sin/cos ptr */
			}while(--loopA & 1);
			pfing(rcx += 128);		/* Skip pad bytes */
		}
		rsi = daddr; //rsi += -32*2*dist128;	/* Restore source pointer */

		/* Do FFT levels 9,10 */
		/* */
		/* On input the 64-byte cache lines hold these data values: */
		/*	0	+32K	64	+32K	128	+32K	192	+32K */
		/*	1	... */
		/*	... */
		/*	63	... */
		/*	256	... */
		/*	... */
		/* On output the 64-byte cache lines hold these data values: */
		/*	0	+32K	16	+32K	32	+32K	48	+32K */
		/*	1	... */
		/*	... */
		/*	15	... */
		/*	64	... */
		/*	... */

		/* Do 16 eight_reals_fft operations */
		/* Do 48 four_complex_djbfft operations */
		/*	distance between fft data elements is 16 */
		/*	do 16 macros each processing 32 data values */

		loops_init(2048);
		r4_h4cl_eight_reals_four_complex_djbfft_preload;
		for(unsigned int loopA = 16; loopA; loopA--) {
			r4_h4cl_eight_reals_four_complex_djbfft(rsi, 64, 16*64, 32*64, rdi, rdx);
			rdi += 4*XMM_SCD2;		/* Next sin/cos data */
			rdx += 256*XMM_HSCD2;	/* Next sin/cos data */
		}
		rsi += -16*64+64*64;	/* Next source pointer */
		rdi += -16*4*XMM_SCD2;	/* Restore sin/cos ptr */
		rdx += -16*256*XMM_HSCD2;	/* Restore sin/cos ptr */

		/* Do 8128 four_complex_djbfft operations */
		/*	distance between fft data elements is 16 */
		/*	do 2032 macros each processing 32 data values */

		loops_reset;			/* 64 iters of 1 or 2 of 16 */
		r4_x4cl_four_complex_djbfft_preload;
		for(unsigned int loopA = 64*2-1; loopA; ) {
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
		rsi = daddr; //rsi += -64*dist128;	/* Next source pointer */

		/* Do FFT levels 11,12 */
		/* */
		/* On input the 64-byte cache lines hold these data values: */
		/*	0	+32K	16	+32K	32	+32K	48	+32K */
		/*	1	... */
		/*	... */
		/*	15	... */
		/*	64	... */
		/*	... */
		/* On output the 64-byte cache lines hold these data values: */
		/*	0	+32K	4	+32K	8	+32K	12	+32K */
		/*	1	... */
		/*	... */
		/*	3	... */
		/*	16	... */
		/*	... */

		/* Do 4 eight_reals_fft operations */
		/* Do 12 four_complex_djbfft operations */
		/*	distance between fft data elements is 4 */
		/*	do 4 macros each processing 32 data values */

		loops_init(2048);
		r4_h4cl_eight_reals_four_complex_djbfft_preload;
		for(unsigned int loopA = 4; loopA; loopA--) {
			r4_h4cl_eight_reals_four_complex_djbfft(rsi, 64, 4*64, 8*64, rdi, rdx);
			rdi += 16*XMM_SCD2;	/* Next sin/cos data */
			rdx += 1024*XMM_HSCD2;	/* Next sin/cos data */
		}
		rsi += -4*64+16*64;	/* Next source pointer */
		rdi += -4*16*XMM_SCD2;	/* Restore sin/cos ptr */
		rdx += -4*1024*XMM_HSCD2;	/* Restore sin/cos ptr */

		/* Do 8176 four_complex_djbfft operations */
		/*	distance between fft data elements is 4 */
		/*	do 2044 macros each processing 32 data values */

		loops_reset;			/* 64 iters of 7 or 8 of 4 */
		r4_x4cl_four_complex_djbfft_preload;
		for(unsigned int loopA = 64*8-1; loopA; ) {
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
		rsi = daddr; //rsi += -64*dist128;	/* Restore source pointer */

		/* Do FFT levels 13,14 */
		/* */
		/* On input the 64-byte cache lines hold these data values: */
		/*	0	+32K	4	+32K	8	+32K	12	+32K */
		/*	1	... */
		/*	... */
		/*	3	... */
		/*	16	... */
		/*	... */
		/* On output the 64-byte cache lines hold these data values: */
		/*	0	+32K	2	+32K	4	+32K	6	+32K */
		/*	1	... */
		/*	8	... */
		/*	... */

		/* Do 1 eight_reals_fft operations */
		/* Do 3 four_complex_fft operations */
		/*	distance between fft data elements is 1 */

		/* Do 8188 four_complex_fft operations */
		/*	distance between fft data elements is 1 */
		/*	do 2047 macros each processing 32 data values */

		/* Execute the proper middle step */

		if likely(g->ffttype & 2) {
			if likely(!(g->ffttype & 1)) { // xpass2_14_levels_real_2:;
				r4_h4cl_eight_reals_four_complex_with_square(rsi, 4*64, 64, 2*64);
				loops_init(2048, 1, 32);
				r4_x4cl_four_complex_with_square_preload;
				for(unsigned int loopA = 64*32-1; loopA; ) {
					do{
						r4_x4cl_four_complex_with_square(rsi, 4*64, 64, 2*64);
					}while(--loopA & 31);
					rsi += -32*4*64+dist128;	/* Next source pointer */
				}
				//rsi += -64*dist128;	/* Restore source pointer */
			}else{ // xpass2_14_levels_real_3:;
				rbp = g->DIST_TO_MULSRCARG;
				r4_h4cl_eight_reals_four_complex_with_mult(rsi, 4*64, 64, 2*64);
				loops_init(2048, 1, 32);
				for(unsigned int loopA = 64*32-1; loopA; ) {
					do{
						r4_x4cl_four_complex_with_mult(rsi, 4*64, 64, 2*64);
					}while(--loopA & 31);
					rsi += -32*4*64+dist128;	/* Next source pointer */
				}
				//rsi += -64*dist128;	/* Restore source pointer */
			}
		}else{ // xpass2_14_levels_real_1:;
			r4_h4cl_eight_reals_four_complex_fft_final(rsi, 4*64, 64, 2*64);
			L2prefetch128(u8ptr(rcx));
			pfing(rcx += 128);
			loops_init_prefetch(2048, 128, 1, rcx, 1, 32);
			for(unsigned int loopA = 32*2*32-1; loopA; ) {
				do{
					xtouch(xptr(rcx+4096-64));		/* Preload the TLBs */
					do{
						r4_x4cl_four_complex_fft_final(rsi, 4*64, 64, 2*64);
					}while(--loopA & 31);
					rsi += -32*4*64+dist128;	/* Next source pointer */
				}while(loopA & (32*2-1));
				pfing(rcx += 128);		/* Skip pad bytes */
			}
			goto	xpass2_14_real_done;
		}
	}else{ // xpass2_14_levels_real_4:;
		rbp = g->DIST_TO_MULSRCARG;
		pfing(rcx = (uintptr_t)g->data_prefetch); /* Load prefetch pointer */
		xtouch(xptr(rcx));			/* Preload the TLBs */
		r4_h4cl_eight_reals_four_complex_with_mulf(rsi, 4*64, 64, 2*64);
		L2prefetch128(u8ptr(rcx));
		pfing(rcx += 128);
		loops_init_prefetch(2048, 128, 1, rcx, 1, 32);
		for(unsigned int loopA = 32*2*32-1; loopA; ) {
			do{
				xtouch(xptr(rcx+4096-64));		/* Preload the TLBs */
				do{
					r4_x4cl_four_complex_with_mulf(rsi, 4*64, 64, 2*64);
				}while(--loopA & 31);
				rsi += -32*4*64+dist128;	/* Next source pointer */
			}while(loopA & (32*2-1));
			pfing(rcx += 128);		/* Skip pad bytes */
		}
		//rsi += -64*dist128;	/* Restore source pointer */
	}
	rsi = daddr;

	/* Do inverse FFT levels 11,12 */
	/* On input the 64-byte cache lines hold these data values: */
	/*	0	+32K	2	+32K	4	+32K	6	+32K */
	/*	1	... */
	/*	8	... */
	/*	... */
	/* On output the 64-byte cache lines hold these data values: */
	/*	0	+32K	8	+32K	16	+32K	24	+32K */
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
	rdi += 4096*XMM_SCD2;
	rdx = (uintptr_t)g->sincos3;
	loops_init(2048);
	for(unsigned int loopA = 2; loopA; loopA--) {
		r4_h4cl_eight_reals_four_complex_djbunfft(rsi, 64, 2*64, 4*64, rdi, 2*16*XMM_SCD2, rdx, 2*1024*XMM_HSCD2);
		rdi += 16*XMM_SCD2;	/* Next sin/cos data */
		rdx += 1024*XMM_HSCD2;	/* Next sin/cos data */
	}
	rsi += -2*64+8*64;		/* Next source pointer */
	rdi += -2*16*XMM_SCD2;	/* Restore sin/cos ptr */
	rdx += -2*1024*XMM_HSCD2;	/* Restore sin/cos ptr */

	/* Do 8184 four_complex_djbunfft operations */
	/*	distance between fft data elements is 4 */
	/*	do 2046 macros each processing 32 data values */

	loops_reset;			/* 64 iters of 15 or 16 of 2 */
	r4_x4cl_four_complex_djbunfft_preload;
	for(unsigned int loopA = 64*16-1; loopA; ) {
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
	rsi = daddr; //rsi += -64*dist128;	/* Restore source pointer */

	/* Do inverse FFT levels 9,10 */
	/* On input the 64-byte cache lines hold these data values: */
	/*	0	+32K	8	+32K	16	+32K	24	+32K */
	/*	1	... */
	/*	... */
	/*	7	... */
	/*	32	... */
	/*	... */
	/* On output the 64-byte cache lines hold these data values: */
	/*	0	+32K	32	+32K	64	+32K	96	+32K */
	/*	1	... */
	/*	... */
	/*	31	... */
	/*	128	... */
	/*	... */

	/* Do 16 eight_reals_unfft operations */
	/* Do 16 four_complex_djbunfft operations */
	/*	distance between fft data elements is 16 */
	/*	do 8 macros processing 32 data values */

	loops_init_prefetch(2048, 64, 1, rcx);
	xtouch(xptr(rcx+4096-64));		/* Preload the TLBs */
	for(unsigned int loopA = 8; loopA; loopA--) {
		r4_h4cl_eight_reals_four_complex_djbunfft(rsi, 64, 8*64, 16*64, rdi, 8*4*XMM_SCD2, rdx, 8*256*XMM_HSCD2);
		rdi += 4*XMM_SCD2;		/* Next sin/cos data */
		rdx += 256*XMM_HSCD2;	/* Next sin/cos data */
	}
	rsi += -8*64+32*64;	/* Next source pointer */
	rdi += -8*4*XMM_SCD2;	/* Restore sin/cos ptr */
	rdx += -8*256*XMM_HSCD2;	/* Restore sin/cos ptr */

	/* Do 8160 four_complex_djbunfft operations */
	/*	distance between fft data elements is 16 */
	/*	do 2040 macros each processing 32 data values */

	loops_reset;			/* 64 iters of 3 or 4 of 4 */
	r4_x4cl_four_complex_djbunfft_preload;
	for(unsigned int loopA = 16*4*4-1; loopA; ) {
		do{
			xtouch(xptr(rcx+4096-64));		/* Preload the TLBs */
			do{
				for(unsigned int loopD = 8; loopD; loopD--) {
					r4_x4cl_four_complex_djbunfft(rsi, 64, 8*64, 16*64, rdi, 8*4*XMM_SCD2);
					rdi += 4*XMM_SCD2;		/* Next sin/cos data */
				}
				rsi += -8*64+32*64;	/* Next source pointer */
				rdi += -8*4*XMM_SCD2;	/* Restore sin/cos ptr */
			}while(--loopA & 3);
			rsi += -4*32*64+dist128;	/* Next source pointer */
		}while(loopA * (4*4-1));
		pfing(rcx += 128);		/* Skip pad bytes */
	}
	rsi = daddr; //rsi += -64*dist128;	/* Restore source pointer */

	/* Do inverse FFT levels 7,8 */
	/* On input the 64-byte cache lines hold these data values: */
	/*	0	+32K	32	+32K	64	+32K	96	+32K */
	/*	1	... */
	/*	... */
	/*	31	... */
	/*	128	... */
	/*	... */
	/* On output the 64-byte cache lines hold these data values: */
	/*	0	+32K	128	+32K	256	+32K	384	+32K */
	/*	1	... */
	/*	... */
	/*	127	... */
	/*	512	... */
	/*	... */

	/* Do 64 eight_reals_unfft operations */
	/* Do 64 four_complex_djbunfft operations */
	/*	distance between fft data elements is 64 */
	/*	do 32 macros processing 32 data values */

	loops_init_prefetch(2048, 64, 1, rcx);
	xtouch(xptr(rcx+4096-64));		/* Preload the TLBs */
	for(unsigned int loopA = 32; loopA; loopA--) {
		r4_h4cl_eight_reals_four_complex_djbunfft(rsi, 64, 32*64, 64*64, rdi, 32*XMM_SCD2, rdx, 32*64*XMM_HSCD2);
		rdi += XMM_SCD2;		/* Next sin/cos data */
		rdx += 64*XMM_HSCD2;	/* Next sin/cos data */
	}
	rsi += -32*64+dist128;	/* Next source pointer */
	rdi += -32*XMM_SCD2;	/* Restore sin/cos ptr */
	rdx += -32*64*XMM_HSCD2;	/* Restore sin/cos ptr */

	/* Do 8064 four_complex_djbunfft operations */
	/*	distance between fft data elements is 64 */
	/*	do 2016 macros each processing 32 data values */

	loops_reset;			/* 63 iters of 32 */
	r4_x4cl_four_complex_djbunfft_preload;
	for(unsigned int loopA = 16*4-1; loopA; ) {
		do{
			xtouch(xptr(rcx+4096-64));		/* Preload the TLBs */
			for(unsigned int loopC = 32; loopC; loopC--) {
				r4_x4cl_four_complex_djbunfft(rsi, 64, 32*64, 64*64, rdi, 32*XMM_SCD2);
				rdi += XMM_SCD2;		/* Next sin/cos data */
			}
			rsi += -32*64+dist128;	/* Next source pointer */
			rdi += -32*XMM_SCD2;	/* Restore sin/cos ptr */
		}while(--loopA & 3);
		pfing(rcx += 128);		/* Skip pad bytes */
	}
	rsi = daddr; //rsi += -64*dist128;	/* Next source pointer */

	/* Do inverse FFT levels 5,6 */
	/* On input the 64-byte cache lines hold these data values: */
	/*	0	+32K	128	+32K	256	+32K	384	+32K */
	/*	1	... */
	/*	... */
	/*	127	... */
	/*	512	... */
	/*	... */
	/* On output the 64-byte cache lines hold these data values: */
	/*	0	+32K	512	+32K	1024	+32K	1536	+32K */
	/*	1	... */
	/*	... */
	/*	511	... */
	/*	2048	... */
	/*	... */

	/* Do 256 eight_reals_unfft operations */
	/* Do 256 four_complex_djbunfft operations */
	/*	distance between fft data elements is 256 */
	/*	do 128 macros processing 32 data values */

	loops_init(2048);
	rdi = (uintptr_t)g->xsincos_complex;	/* Load sin/cos pointer */
	for(unsigned int loopA = 128; loopA; loopA--) {
		r4_h4cl_eight_reals_four_complex_djbunfft(rsi, 64, dist128, 2*dist128, rdi, 128*16*XMM_SCD2, rdx, 128*16*XMM_HSCD2);
		rdi += 16*XMM_SCD2;	/* Next sin/cos data */
		rdx += 16*XMM_HSCD2;	/* Next sin/cos data */
	}
	rsi += -128*64+4*dist128;	/* Next source pointer */
	rdi += -128*16*XMM_SCD2;	/* Restore sin/cos ptr */
	rdx += -128*16*XMM_HSCD2;	/* Restore sin/cos ptr */

	/* Do 7680 four_complex_djbunfft operations */
	/*	distance between fft data elements is 256 */
	/*	do 1920 macros each processing 32 data values */

	loops_reset;			/* 15 iters of 128 */
	r4_x4cl_four_complex_djbunfft_preload;
	for(unsigned int loopA = 15; loopA; loopA--) {
		for(unsigned int loopB = 128; loopB; loopB--) {
			r4_x4cl_four_complex_djbunfft(rsi, 64, dist128, 2*dist128, rdi, 128*16*XMM_SCD2);
			rdi += 16*XMM_SCD2;	/* Next sin/cos data */
		}
		rsi += -128*64+4*dist128;	/* Next source pointer */
		rdi += -128*16*XMM_SCD2;	/* Restore sin/cos ptr */
	}
	rsi = daddr; //rsi += -16*4*dist128;	/* Next source pointer */

	/* Do inverse FFT levels 3,4 */
	/* On input the 64-byte cache lines hold these data values: */
	/*	0	+32K	512	+32K	1024	+32K	1536	+32K */
	/*	1	... */
	/*	... */
	/*	511	... */
	/*	2048	... */
	/*	... */
	/* On output the 64-byte cache lines hold these data values: */
	/*	0	+32K	2048	+32K	4096	+32K	6144	+32K */
	/*	1	... */
	/*	... */
	/*	2047	... */
	/*	8192	... */
	/*	... */

	/* Do 1024 eight_reals_unfft operations */
	/* Do 1024 four_complex_djbunfft operations */
	/*	distance between fft data elements is 1024 */
	/*	do 512 macros processing 32 data values */

	loops_init(2048);
	rdi = (uintptr_t)g->xsincos_complex;	/* Load sin/cos pointer */
	for(unsigned int loopA = 4; loopA; loopA--) {
		for(unsigned int loopB = 128; loopB; loopB--) {
			r4_h4cl_eight_reals_four_complex_djbunfft(rsi, 64, 4*dist128, 8*dist128, rdi, 512*4*XMM_SCD2, rdx, 512*4*XMM_HSCD2);
			rdi += 4*XMM_SCD2;		/* Next sin/cos data */
			rdx += 4*XMM_HSCD2;	/* Next sin/cos data */
		}
		rsi += -128*64+dist128;	/* Next source pointer */
	}
	rsi += -4*dist128+16*dist128; /* Next source pointer */
	rdi += -512*4*XMM_SCD2;	/* Restore sin/cos ptr */
	rdx += -512*4*XMM_HSCD2;	/* Restore sin/cos ptr */

	/* Do 6144 four_complex_djbunfft operations */
	/*	distance between fft data elements is 1024 */
	/*	do 1536 macros each processing 32 data values */

	loops_reset;			/* 3 iters of 512 */
	r4_x4cl_four_complex_djbunfft_preload;
	for(unsigned int loopA = 3; loopA; loopA--) {
		for(unsigned int loopB = 4; loopB; loopB--) {
			for(unsigned int loopC = 128; loopC; loopC--) {
				r4_x4cl_four_complex_djbunfft(rsi, 64, 4*dist128, 8*dist128, rdi, 512*4*XMM_SCD2);
				rdi += 4*XMM_SCD2;		/* Next sin/cos data */
			}
			rsi += -128*64+dist128;	/* Next source pointer */
		}
		rsi += -4*dist128+16*dist128; /* Next source pointer */
		rdi += -512*4*XMM_SCD2;	/* Restore sin/cos ptr */
	}
	rsi = daddr; //rsi += -4*16*dist128;	/* Restore source pointer */

	/* Do inverse FFT levels 1,2 */
	/* On input the 64-byte cache lines hold these data values: */
	/*	0	+32K	2048	+32K	4096	+32K	6144	+32K */
	/*	1	... */
	/*	... */
	/*	2047	... */
	/*	8192	... */
	/*	... */
	/* On output the 64-byte cache lines hold these data values if pass 1 swizzles: */
	/*	0	+32K	8192	+32K	16384	+32K	24576	+32K */
	/*	1	... */
	/*	... */
	/*	8191	... */
	/* On output the 64-byte cache lines hold these data values if pass 2 swizzles: */
	/*	0	+8K	16384	+8K	32768	+8K	49152	+8K */
	/*	1	... */
	/*	... */
	/*	8191	... */

	/* Do 4096 eight_reals_unfft operations */
	/* Do 4096 four_complex_djbunfft operations */
	/*	distance between fft data elements is 4096 */
	/*	do 2048 macros processing 32 data values */

	loops_init(2048);
	for(unsigned int loopA = 16; loopA; loopA--) {
		for(unsigned int loopB = 128; loopB; loopB--) {
			r4_h4cl_eight_reals_four_complex_djbunfft(rsi, 64, 16*dist128, 32*dist128, rdi, 2048*XMM_SCD2, rdx, 2048*XMM_HSCD2);
			rdi += XMM_SCD2;		/* Next sin/cos data */
			rdx += XMM_HSCD2;		/* Next sin/cos data */
		}
		rsi += -128*64+dist128;	/* Next source pointer */
	}

	xpass2_14_real_done:;
	end_timer(2);
	if(g->pass2_get_next_block(g) != 0) return;
}
r4_pass2_14_levels_complex(g);
}



	/* *************** radix4-pass2-14-levels-complex macro ****************** */
	/* This macro takes 2 groups of 16384 complex values and */
	/* performs the final 14 levels of a standard radix-4 FFT, squares the */
	/* results, and does the first 14 levels of the inverse FFT. */


void r4_pass2_14_levels_complex(struct gwasm_data *__restrict g) {
do{
	uintptr_t rbx,rdx,rbp,rsi,rdi,daddr;
	vec2f64 xmm0,xmm1,xmm2,xmm3,xmm4,xmm5,xmm6,xmm7;

	/* Type 4 FFTs skip the forward FFT process */

	daddr = (uintptr_t)g->data_addr;/* Load address of FFT data */
	rsi = daddr;
	rbx = g->DIST_TO_FFTSRCARG;
	if likely(!(g->ffttype & 4)) { //xpass2_14_levels_complex_4;

		/* Do FFT level 1,2 */
		/* */
		/* We swizzle the data so that one set of data is in the high word of */
		/* an xmm register, the other set is in the low word.  This means */
		/* a distance of 8192 data values. */
		/* */
		/* On input the 64-byte cache lines hold these data values if pass 2 swizzles: */
		/*	0	+8K	16384	+8K	32768	+8K	49152	+8K */
		/*	1	... */
		/*	... */
		/*	8191	... */
		/* On input the 64-byte cache lines hold these data values if pass 1 swizzles: */
		/*	0	+32K	8192	+32K	16384	+32K	24576	+32K */
		/*	1	... */
		/*	... */
		/*	8191	... */
		/* On output the 64-byte cache lines hold these data values: */
		/*	0	+32K	4096	+32K	8192	+32K	12288	+32K */
		/*	1	... */
		/*	... */
		/*	4095	... */
		/*	16384	... */
		/*	... */

		/* Do 8192 four_complex_fft operations */
		/*	distance between fft data elements is 4096 */
		/*	do 4096 macros each processing 16 data values */

		/* Do FFT levels 3,4 */
		/* */
		/* On input the 64-byte cache lines hold these data values: */
		/*	0	+32K	4096	+32K	8192	+32K	12288	+32K */
		/*	1	... */
		/*	... */
		/*	4095	... */
		/*	16384	... */
		/*	... */
		/* On output the 64-byte cache lines hold these data values: */
		/*	0	+32K	1024	+32K	2048	+32K	3072	+32K */
		/*	1	... */
		/*	... */
		/*	1023	... */
		/*	4096	... */
		/*	... */

		/* Do 8192 four_complex_fft operations */
		/*	distance between fft data elements is 1024 */
		/*	do 2048 macros each processing 32 data values */

		/* Do FFT levels 5,6 */
		/* */
		/* On input the 64-byte cache lines hold these data values: */
		/*	0	+32K	1024	+32K	2048	+32K	3072	+32K */
		/*	1	... */
		/*	... */
		/*	1023	... */
		/*	4096	... */
		/*	... */
		/* On output the 64-byte cache lines hold these data values: */
		/*	0	+32K	256	+32K	512	+32K	768	+32K */
		/*	1	... */
		/*	... */
		/*	255	... */
		/*	1024	... */
		/*	... */

		/* Do 8192 four_complex_fft operations */
		/*	distance between fft data elements is 256 */
		/*	do 2048 macros each processing 32 data values */

#if 0 // @INSTR(,%xarch,<CORE>) NE 0; /* Core 2 faster doing separate loops */

		start_timer(3);
		rdi = (uintptr_t)g->xsincos_complex;	/* The sin/cos data */
		loops_init(4096);			/* 8 iterations of 128 */
		r4_f2cl_four_complex_djbfft_preload;
		for(unsigned int loopA = 32; loopA; loopA--) {
			for(unsigned int loopB = 128; loopB; loopB--) {
				r4_f2cl_four_complex_djbfft(rsi, 64, 32*dist128, rdi);
				rdi += XMM_SCD2;		/* Next sin/cos data */
			}
			rsi += -128*64+dist128;	/* Next source pointer */
		}
		rsi = daddr; //rsi += -32*dist128;	/* Restore source pointer */
		end_timer(3);

		start_timer(4);
		rdi = (uintptr_t)g->xsincos_complex;	/* Ptr to the sin/cos data */
		pfing(rcx = (uintptr_t)g->data_prefetch); /* Prefetch pointer */
		xtouch(xptr(rcx));			/* Load first TLB */
		loops_init_prefetch(2048, 64, 2, rcx); /* 2 iters of 8 iters of 128 */
		r4_x4cl_four_complex_djbfft_preload;
		for(unsigned int loopA = 2*4*2; loopA; ) {
			do{
				do{
					xtouch(xptr(rcx+4096-64));		/* Load prefetch TLB */
					for(unsigned int loopD = 128; loopD; loopD--) {
						r4_x4cl_four_complex_djbfft(rsi, 64, 8*dist128, 16*dist128, rdi);
						rdi += 4*XMM_SCD2;		/* Next sin/cos data */
					}
					rsi += -128*64+dist128;	/* Next source pointer */
				}while(--loopA & 1);
				pfing(rcx += 128);		/* Skip pad bytes */
			}while(loopA & (4*2-1));
			rsi += -8*dist128+32*dist128; /* Next source pointer */
			rdi += -1024*4*XMM_SCD2;	/* Restore sin/cos ptr */
		}
		rsi = daddr; //rsi += -2*32*dist128;	/* Restore source pointer */
		end_timer(4);

		start_timer(5);
		loops_init(2048);			/* 8 iters of 2 iters of 128 */
		r4_x4cl_four_complex_djbfft_preload;
		for(unsigned int loopA = 8; loopA; loopA--) {
			for(unsigned int loopB = 2; loopB; loopB--) {
				for(unsigned int loopC = 128; loopC; loopC--) {
					r4_x4cl_four_complex_djbfft(rsi, 64, 2*dist128, 4*dist128, rdi);
					rdi += 16*XMM_SCD2;	/* Next sin/cos data */
				}
				rsi += -128*64+dist128;	/* Next source pointer */
			}
			rsi += -2*dist128+8*dist128; /* Next source pointer */
			rdi += -256*16*XMM_SCD2;	/* Restore sin/cos ptr */
		}
		rsi = daddr; //rsi += -8*8*dist128;	/* Restore source pointer */
		end_timer(5);

#else /* Pentium 4 faster doing merged loops */

		start_timer(3);
		rdi = (uintptr_t)g->xsincos_complex;	/* The sin/cos data */
		rbp = rdi;
		rdx = rdi;
		pfing(rcx = (uintptr_t)g->data_prefetch); /* Prefetch pointer */
		xtouch(xptr(rcx));			/* Load first TLB */

		loops_init_prefetch(4096, 64, 4, rcx); /* 256 iterations of 32 */
		;	/*r4_f2cl_four_complex_djbfft_preload -- assume same as r4_x4cl_four_complex_djbfft_preload */
		r4_x4cl_four_complex_djbfft_preload;

		for(unsigned int loopA = 2*4*32; loopA; ) {
			do{
				xtouch(xptr(rcx+4096-64));		/* Load prefetch TLB */
				do{
					for(unsigned int loopD = 16; loopD; loopD--) {
						r4_f2cl_four_complex_djbfft(rsi, 2*dist128, 32*dist128, rdi);
						xxprefetch(u8ptr(rdi+XMM_SCD2));	/* Prefetch for next pass */
						rdi += 256*XMM_SCD2;	/* Next sin/cos data */
					}
					rsi += -16*2*dist128;	/* Restore source pointer */
					loops_undo(16);

					for(unsigned int loopD = 2; loopD; loopD--) {
						for(unsigned int loopE = 4; loopE; loopE--) {
							r4_x4cl_four_complex_djbfft(rsi, 2*dist128, 8*dist128, 16*dist128, rbp);
							rbp += 256*4*XMM_SCD2;	/* Next sin/cos data */
						}
						rsi += -4*2*dist128+32*dist128; /* Next source pointer */
						rbp += -4*256*4*XMM_SCD2;	/* Restore sin/cos ptr */
					}
					rsi += -2*32*dist128;	/* Restore source pointer */
					loops_undo(16);

					for(unsigned int loopD = 8; loopD; loopD--) {
						r4_x4cl_four_complex_djbfft(rsi, 8*dist128, 2*dist128, 4*dist128, rdx);
					}

					rsi += -8*8*dist128+64;	/* Next source pointer */
					rdi += -16*256*XMM_SCD2+XMM_SCD2; /* Next sin/cos data */
					rbp += 4*XMM_SCD2;		/* Next sin/cos data */
					rdx += 16*XMM_SCD2;	/* Next sin/cos data */
				}while(--loopA & (32-1));
				pfing(rcx += 128);		/* Skip pad bytes */
			}while(loopA & (4*32-1));
			rsi += -128*64+dist128;	/* Next source pointer */
		}
		rsi = daddr; //rsi += -2*dist128;		/* Restore source pointer */
		rdi += -256*XMM_SCD2;	/* Restore sin/cos data */
		end_timer(3);

#endif

		/* Do FFT levels 7,8 */
		/* */
		/* On input the 64-byte cache lines hold these data values: */
		/*	0	+32K	256	+32K	512	+32K	768	+32K */
		/*	1	... */
		/*	... */
		/*	255	... */
		/*	1024	... */
		/*	... */
		/* On output the 64-byte cache lines hold these data values: */
		/*	0	+32K	64	+32K	128	+32K	192	+32K */
		/*	1	... */
		/*	... */
		/*	63	... */
		/*	256	... */
		/*	... */
		/* */
		/* Do 8192 four_complex_fft operations */
		/*	distance between fft data elements is 64 */
		/*	do 2048 macros each processing 32 data values */

		/* Do FFT levels 9,10 */
		/* */
		/* On input the 64-byte cache lines hold these data values: */
		/*	0	+32K	64	+32K	128	+32K	192	+32K */
		/*	1	... */
		/*	... */
		/*	63	... */
		/*	256	... */
		/*	... */
		/* On output the 64-byte cache lines hold these data values: */
		/*	0	+32K	16	+32K	32	+32K	48	+32K */
		/*	1	... */
		/*	... */
		/*	15	... */
		/*	64	... */
		/*	... */
		/* */
		/* Do 8192 four_complex_fft operations */
		/*	distance between fft data elements is 16 */
		/*	do 2048 macros each processing 32 data values */

		/* Do FFT levels 11,12 */
		/* */
		/* On input the 64-byte cache lines hold these data values: */
		/*	0	+32K	16	+32K	32	+32K	48	+32K */
		/*	1	... */
		/*	... */
		/*	15	... */
		/*	64	... */
		/*	... */
		/* On output the 64-byte cache lines hold these data values: */
		/*	0	+32K	4	+32K	8	+32K	12	+32K */
		/*	1	... */
		/*	... */
		/*	3	... */
		/*	16	... */
		/*	... */
		/* */
		/* Do 8192 four_complex_fft operations */
		/*	distance between fft data elements is 4 */
		/*	do 2048 macros each processing 32 data values */

		/* Do this 64 cache lines at a time to maximize accessing data */
		/* while it is in the L1 cache. */

#if 0 // NEVER /* Never faster to do separate loops */

		start_timer(6);
		rdi += 4096*XMM_SCD2;	/* Load sin/cos pointer */
		loops_init_prefetch(2048, 64, 2, rcx); /* 8*4 iters of 64 */
		r4_x4cl_four_complex_djbfft_preload;
		for(unsigned int loopA = 8; loopA; loopA--) {
			xtouch(xptr(rcx+4096-128));		/* Load two prefetch TLBs */
			xtouch(xptr(rcx+2*4096-128));
			for(unsigned int loopB = 4; loopB; loopB--) {
				for(unsigned int loopC = 64; loopC; loopC--) {
					r4_x4cl_four_complex_djbfft(rsi, 64, 64*64, dist128, rdi);
					rdi += XMM_SCD2;		/* Next sine/cosine pointer */
				}
				rsi += -64*64+2*dist128;	/* Next source pointer */
				rdi += -64*XMM_SCD2;	/* Restore sine/cosine pointer */
			}
			rcx += -128*64+8192+128;	/* Next 8KB section to prefetch */
		}
		rsi = daddr; //rsi += -8*4*2*dist128;	/* Restore source pointer */
		end_timer(6);

		start_timer(7);
		loops_init_prefetch(2048, 64, 2, rcx); /* 8*8*2 iters of 16 */
		r4_x4cl_four_complex_djbfft_preload;
		for(unsigned int loopA = 8; loopA; loopA--) {
			xtouch(xptr(rcx+4096-128));		/* Load two prefetch TLBs */
			xtouch(xptr(rcx+2*4096-128));
			for(unsigned int loopB = 8; loopB; loopB--) {
				for(unsigned int loopC = 2; loopC; loopC--) {
					for(unsigned int loopD = 16; loopD; loopD--) {
						r4_x4cl_four_complex_djbfft(rsi, 64, 16*64, 32*64, rdi);
						rdi += 4*XMM_SCD2;		/* Next sine/cosine pointer */
					}
					rsi += -16*64+64*64;	/* Next source pointer */
					rdi += -16*4*XMM_SCD2;	/* Restore sine/cosine pointer */
				}
				rsi += -2*64*64+dist128;	/* Adjust source pointer */
			}
			rcx += -128*64+8192+128;	/* Next 8KB section to prefetch */
		}
		rsi = daddr; //rsi += -8*8*dist128;	/* Restore source pointer */
		end_timer(7);

		start_timer(8);
		loops_init_prefetch(2048, 64, 2, rcx); /* 8*8*8 iters of 4 */
		r4_x4cl_four_complex_djbfft_preload;
		for(unsigned int loopA = 8; loopA; loopA--) {
			xtouch(xptr(rcx+4096-128));		/* Load two prefetch TLBs */
			xtouch(xptr(rcx+2*4096-128));
			for(unsigned int loopB = 8; loopB; loopB--) {
				for(unsigned int loopC = 8; loopC; loopC--) {
					for(unsigned int loopD = 4; loopD; loopD--) {
						r4_x4cl_four_complex_djbfft(rsi, 64, 4*64, 8*64, rdi);
						rdi += 16*XMM_SCD2;	/* Next sine/cosine pointer */
					}
					rsi += -4*64+16*64;	/* Next source pointer */
					rdi += -4*16*XMM_SCD2;	/* Restore sine/cosine pointer */
				}
				rsi += -8*16*64+dist128;	/* Adjust source pointer */
			}
			rcx += -128*64+8192+128;	/* Next 8KB section to prefetch */
		}
		rsi = daddr; //rsi += -8*8*dist128;	/* Restore source pointer */
		end_timer(8);

#else /* Always faster to do merged loops */

		start_timer(6);
		rdi += 4096*XMM_SCD2;	/* Load sin/cos pointers */
		rbx = rdi;
		rbp = rdi;
		loops_init_prefetch(2048, 64, 2, rcx); /* 2*4 iters of 4 of 64 cache lines */
		r4_x4cl_four_complex_djbfft_preload;

		for(unsigned int loopA = 8*4*4; loopA; ) {
			xtouch(xptr(rcx+4096-128));		/* Load six prefetch TLBs */
			xtouch(xptr(rcx+2*4096-128));
			xtouch(xptr(rcx+3*4096-128));
			xtouch(xptr(rcx+4*4096-128));
			xtouch(xptr(rcx+5*4096-128));
			xtouch(xptr(rcx+6*4096-128));

			do{
				do{
					for(unsigned int loopD = 16; loopD; loopD--) {
						xxprefetch(u8ptr(rdi+4*XMM_SCD2));	/* Prefetch for next pass */
						r4_x4cl_four_complex_djbfft(rsi, 4*64, 64*64, dist128, rdi);
						rdi += 4*XMM_SCD2;		/* Next sine/cosine pointer */
					}
					rsi += -16*4*64;		/* Restore source pointer */
					pfing(rcx += -8*64+8192+128);	/* Prefetch from next section */
					loops_undo(16);			/* Back up loop counter */

					for(unsigned int loopD = 2; loopD; loopD--) {
						for(unsigned int loopE = 2; loopE; loopE--) {
							for(unsigned int loopF = 4; loopF; loopF--) {
								r4_x4cl_four_complex_djbfft(rsi, 4*64, 16*64, 32*64, rbx);
								rbx += 4*4*XMM_SCD2;	/* Next sine/cosine pointer */
							}
							rsi += -4*4*64+64*64;	/* Next source pointer */
							rbx += -4*4*4*XMM_SCD2;	/* Restore sine/cosine pointer */
						}
						rsi += -2*64*64+dist128;	/* Next source pointer */
					}
					rsi += -2*dist128;		/* Restore source pointer */
					pfing(rcx += -8*64+8192+128);	/* Prefetch from next section */
					loops_undo(2*2*4);		/* Back up loop counter */

					for(unsigned int loopD = 2; loopD; loopD--) {
						for(unsigned int loopE = 8; loopE; loopE--) {
							r4_x4cl_four_complex_djbfft(rsi, 16*64, 4*64, 8*64, rbp);
						}
						rsi += -8*16*64+dist128;	/* Next source pointer */
					}
					pfing(rcx += -2*(8192+128));	/* Back to prefetching from 1st section */

					rsi += -2*dist128+64;	/* Next source pointer */
					rdi += -16*4*XMM_SCD2+XMM_SCD2; /* Next sine/cosine pointer */
					rbx += 4*XMM_SCD2;		/* Next sine/cosine pointer */
					rbp += 16*XMM_SCD2;	/* Next sine/cosine pointer */
				}while(--loopA & (4-1));
				rsi += -4*64+2*dist128;	/* Next source pointer */
				rdi += -4*XMM_SCD2;	/* Restore sine/cosine pointer */
				rbx += -4*4*XMM_SCD2;	/* Restore sine/cosine pointer */
				rbp += -4*16*XMM_SCD2;	/* Restore sine/cosine pointer */
			}while(loopA & (4*4-1));
			pfing(rcx += -16*8*64+3*(8192+128)); /* Prefetch next 3 sections */
		}
		rsi = daddr; //rsi += -32*2*dist128;	/* Restore source pointer */
		end_timer(6);
#endif

		/* Do FFT levels 13,14 as well as inverse FFT levels 13,14 */
		/* */
		/* On input the 64-byte cache lines hold these data values: */
		/*	0	+32K	4	+32K	8	+32K	12	+32K */
		/*	1	... */
		/*	... */
		/*	3	... */
		/*	16	... */
		/*	... */
		/* On output the 64-byte cache lines hold these data values: */
		/*	0	+32K	2	+32K	4	+32K	6	+32K */
		/*	1	... */
		/*	8	... */
		/*	... */

		/* Do 8192 four_complex_fft operations */
		/*	distance between fft data elements is 1 */
		/*	do 2048 macros each processing 32 data values */

		/* Execute the right middle step */

		if likely(g->ffttype & 2) {
			if likely(!(g->ffttype & 1)) { // xpass2_14_levels_complex_2:;
				start_timer(9);
				loops_init(2048);
				r4_x4cl_four_complex_with_square_preload;
				for(unsigned int loopA = 64; loopA; loopA--) {
					for(unsigned int loopB = 32; loopB; loopB--) {
						r4_x4cl_four_complex_with_square(rsi, 4*64, 64, 2*64);
					}
					rsi += -32*4*64+dist128;	/* Next source pointer */
				}
				rsi += -64*dist128;	/* Restore source pointer */
				end_timer(9);
			}else{ // xpass2_14_levels_complex_3:;
				start_timer(9);
				rbp = g->DIST_TO_MULSRCARG;
				loops_init(2048);
				for(unsigned int loopA = 64; loopA; loopA--) {
					for(unsigned int loopB = 32; loopB; loopB--) {
						r4_x4cl_four_complex_with_mult(rsi, 4*64, 64, 2*64);
					}
					rsi += -32*4*64+dist128;	/* Next source pointer */
				}
				rsi += -64*dist128;	/* Restore source pointer */
				end_timer(9);
			}
		}else{ // xpass2_14_levels_complex_1:;
			start_timer(9);
			loops_init_prefetch(2048, 128, 1, rcx);
			for(unsigned int loopA = 64; loopA; ) {
				do{
					xtouch(xptr(rcx+4096-128));		/* Load prefetch TLB */
					for(unsigned int loopC = 32; loopC; loopC--) {
						r4_x4cl_four_complex_fft_final(rsi, 4*64, 64, 2*64);
					}
					rsi += -32*4*64+dist128;	/* Next source pointer */
				}while(--loopA & 1);
				pfing(rcx += 128);		/* Prefetch next section */
			}
			end_timer(9);
			goto	xpass2_14_levels_complex_done;
		}
	}else{ // xpass2_14_levels_complex_4:;
		start_timer(9);
		pfing(rcx = (uintptr_t)g->data_prefetch); /* Prefetch pointer */
		xtouch(xptr(rcx));			/* Load first prefetch TLB */
		rbp = g->DIST_TO_MULSRCARG;
		loops_init_prefetch(2048, 128, 1, rcx);
		for(unsigned int loopA = 64; loopA; ) {
			do{
				xtouch(xptr(rcx+4096-64));		/* Load prefetch TLB */
				for(unsigned int loopC = 32; loopC; loopC--) {
					r4_x4cl_four_complex_with_mulf(rsi, 4*64, 64, 2*64);
				}
				rsi += -32*4*64+dist128;	/* Next source pointer */
			}while(--loopA & 1);
			pfing(rcx += 128);		/* Prefetch next section */
		}
		rsi += -64*dist128;	/* Restore source pointer */
		end_timer(9);
	}
	/* Do inverse FFT levels 11,12 */
	/* On input the 64-byte cache lines hold these data values: */
	/*	0	+32K	2	+32K	4	+32K	6	+32K */
	/*	1	... */
	/*	8	... */
	/*	... */
	/* On output the 64-byte cache lines hold these data values: */
	/*	0	+32K	8	+32K	16	+32K	24	+32K */
	/*	1	... */
	/*	... */
	/*	7	... */
	/*	32	... */
	/*	... */
	/* */
	/* Do 8192 four_complex_unfft operations */
	/*	distance between fft data elements is 4 */
	/*	do 2048 macros each processing 32 data values */

	/* Do inverse FFT levels 9,10 */
	/* On input the 64-byte cache lines hold these data values: */
	/*	0	+32K	8	+32K	16	+32K	24	+32K */
	/*	1	... */
	/*	... */
	/*	7	... */
	/*	32	... */
	/*	... */
	/* On output the 64-byte cache lines hold these data values: */
	/*	0	+32K	32	+32K	64	+32K	96	+32K */
	/*	1	... */
	/*	... */
	/*	31	... */
	/*	128	... */
	/*	... */
	/* */
	/* Do 8192 four_complex_unfft operations */
	/*	distance between fft data elements is 16 */
	/*	do 2048 macros each processing 32 data values */

	/* Do inverse FFT levels 7,8 */
	/* On input the 64-byte cache lines hold these data values: */
	/*	0	+32K	32	+32K	64	+32K	96	+32K */
	/*	1	... */
	/*	... */
	/*	31	... */
	/*	128	... */
	/*	... */
	/* On output the 64-byte cache lines hold these data values: */
	/*	0	+32K	128	+32K	256	+32K	384	+32K */
	/*	1	... */
	/*	... */
	/*	127	... */
	/*	512	... */
	/*	... */
	/* */
	/* Do 8192 four_complex_unfft operations */
	/*	distance between fft data elements is 64 */
	/*	do 2048 macros each processing 32 data values */

	/* Do this 64 cache lines at a time to maximize accessing data */
	/* while it is in the L1 cache. */

#if 0 // NEVER /* Never faster to do separate loops */

	start_timer(10);
	rdi = (uintptr_t)g->xsincos_complex;	/* Load sin/cos pointers */
	rdi += 4096*XMM_SCD2;
	loops_init_prefetch(2048, 64, 2, rcx);/* 8 iters of 8 iters of 16 iters of 2 */
	r4_x4cl_four_complex_djbunfft_preload;
	for(unsigned int loopA = 8; loopA; loopA--) {
		xtouch(xptr(rcx+4096-128));		/* Load prefetch TLBs */
		xtouch(xptr(rcx+2*4096-128));
		for(unsigned int loopB = 8; loopB; loopB--) {
			for(unsigned int loopC = 16; loopC; loopC--) {
				for(unsigned int loopD = 2; loopD; loopD--) {
					r4_x4cl_four_complex_djbunfft(rsi, 64, 2*64, 4*64, rdi, 2*16*XMM_SCD2);
					rdi += 16*XMM_SCD2;	/* Next sine/cosine pointer */
				}
				rsi += -2*64+8*64;		/* Next source pointer */
				rdi += -2*16*XMM_SCD2;	/* Reset sine/cosine pointer */
			}
			rsi += -16*8*64+dist128;	/* Next source pointer */
		}
		rcx += -128*64+8192+128;	/* Next 8KB section to prefetch */
	}
	rsi = daddr; //rsi += -8*8*dist128;	/* Restore source pointer */
	end_timer(10);

	start_timer(11);
	loops_init_prefetch(2048, 64, 2, rcx);/* 8 iters of 8 iters of 4 iters of 8 */
	r4_x4cl_four_complex_djbunfft_preload;
	for(unsigned int loopA = 8; loopA; loopA--) {
		xtouch(xptr(rcx+4096-128));		/* Load prefetch TLBs */
		xtouch(xptr(rcx+2*4096-128));
		for(unsigned int loopB = 8; loopB; loopB--) {
			for(unsigned int loopC = 4; loopC; loopC--) {
				for(unsigned int loopD = 8; loopD; loopD--) {
					r4_x4cl_four_complex_djbunfft(rsi, 64, 8*64, 16*64, rdi, 8*4*XMM_SCD2);
					rdi += 4*XMM_SCD2;		/* Next sine/cosine pointer */
				}
				rsi += -8*64+32*64;	/* Next source pointer */
				rdi += -8*4*XMM_SCD2;	/* Reset sine/cosine pointer */
			}
			rsi += -4*32*64+dist128;	/* Next source pointer */
		}
		rcx += -128*64+8192+128;	/* Next 8KB section to prefetch */
	}
	rsi = daddr; //rsi += -8*8*dist128;	/* Restore source pointer */
	end_timer(11);

	start_timer(12);
	loops_init_prefetch(2048, 64, 2, rcx);/* 8 iters of 8 iters of 32 */
	r4_x4cl_four_complex_djbunfft_preload;
	for(unsigned int loopA = 8; loopA; loopA--) {
		xtouch(xptr(rcx+4096-128));		/* Load prefetch TLBs */
		xtouch(xptr(rcx+2*4096-128));
		for(unsigned int loopB = 8; loopB; loopB--) {
			for(unsigned int loopC = 32; loopC; loopC--) {
				r4_x4cl_four_complex_djbunfft(rsi, 64, 32*64, 64*64, rdi, 32*1*XMM_SCD2);
				rdi += 1*XMM_SCD2;		/* Next sine/cosine pointer */
			}
			rsi += -32*64+dist128;	/* Next source pointer */
			rdi += -32*1*XMM_SCD2;	/* Reset sine/cosine pointer */
		}
		rcx += -128*64+8192+128;	/* Next 8KB section to prefetch */
	}
	rsi = daddr; //rsi += -8*8*dist128;	/* Restore source pointer */
	end_timer(12);

#else /* Always faster to do merged loops */

	start_timer(10);
	rdi = (uintptr_t)g->xsincos_complex;	/* Load sin/cos pointers */
	rdi += 4096*XMM_SCD2;
	rbx = rdi;
	rbp = rdi;
	loops_init_prefetch(2048, 64, 2, rcx); /* 8*2*8 iters of 64 cache lines */
	r4_x4cl_four_complex_djbunfft_preload;

	for(unsigned int loopA = 8*8*2; loopA; ) {
		xtouch(xptr(rcx+4096-128));		/* Load six prefetch TLBs */
		xtouch(xptr(rcx+2*4096-128));
		xtouch(xptr(rcx+3*4096-128));
		xtouch(xptr(rcx+4*4096-128));
		xtouch(xptr(rcx+5*4096-128));
		xtouch(xptr(rcx+6*4096-128));

		do{
			do{
				for(unsigned int loopD = 16; loopD; loopD--) {
					r4_x4cl_four_complex_djbunfft(rsi, 8*64, 2*64, 4*64, rbp, 2*16*XMM_SCD2);
				}
				rsi += -16*8*64;		/* Restore source pointer */
				pfing(rcx += -8*64+8192+128);	/* Prefetch from next section */
				loops_undo(16);			/* Restore loop counter */

				for(unsigned int loopD = 4; loopD; loopD--) {
					for(unsigned int loopE = 4; loopE; loopE--) {
						r4_x4cl_four_complex_djbunfft(rsi, 2*64, 8*64, 16*64, rbx, 8*4*XMM_SCD2);
						rbx += 2*4*XMM_SCD2;	/* Next sine/cosine pointer */
					}
					rsi += -4*2*64+32*64;	/* Next source pointer */
					rbx += -4*2*4*XMM_SCD2;	/* Restore sine/cosine pointer */
				}
				rsi += -4*32*64;		/* Restore source pointer */
				pfing(rcx += -8*64+8192+128);	/* Prefetch from next section */
				loops_undo(4*4);			/* Restore loop counter */

				for(unsigned int loopD = 16; loopD; loopD--) {
					r4_x4cl_four_complex_djbunfft(rsi, 2*64, 32*64, 64*64, rdi, 32*1*XMM_SCD2);
					rdi += 2*XMM_SCD2;		/* Next sine/cosine pointer */
				}
				pfing(rcx += -2*(8192+128));	/* Back to prefetching from 1st section */

				rsi += -16*2*64+64;	/* Next source pointer */
				rdi += -16*2*XMM_SCD2+XMM_SCD2; /* Next sine/cosine pointer */
				rbx += 4*XMM_SCD2;		/* Next sine/cosine pointer */
				rbp += 16*XMM_SCD2;	/* Next sine/cosine pointer */
			}while(--loopA & (2-1));
			rsi += -2*64+dist128;	/* Next source pointer */
			rdi += -2*XMM_SCD2;	/* Restore sine/cosine pointer */
			rbx += -2*4*XMM_SCD2;	/* Restore sine/cosine pointer */
			rbp += -2*16*XMM_SCD2;	/* Restore sine/cosine pointer */
		}while(loopA & (8*2-1));
		pfing(rcx += -16*8*64+3*(8192+128)); /* Prefetch next 3 sections */
	}
	rsi = daddr; //rsi += -64*dist128;	/* Restore source pointer */
	end_timer(10);
#endif

	/* Do inverse FFT levels 5,6 */
	/* On input the 64-byte cache lines hold these data values: */
	/*	0	+32K	128	+32K	256	+32K	384	+32K */
	/*	1	... */
	/*	... */
	/*	127	... */
	/*	512	... */
	/*	... */
	/* On output the 64-byte cache lines hold these data values: */
	/*	0	+32K	512	+32K	1024	+32K	1536	+32K */
	/*	1	... */
	/*	... */
	/*	511	... */
	/*	2048	... */
	/*	... */

	/* Do 8192 four_complex_unfft operations */
	/*	distance between fft data elements is 256 */
	/*	do 2048 macros each processing 32 data values */

	/* Do inverse FFT levels 3,4 */
	/* On input the 64-byte cache lines hold these data values: */
	/*	0	+32K	512	+32K	1024	+32K	1536	+32K */
	/*	1	... */
	/*	... */
	/*	511	... */
	/*	2048	... */
	/*	... */
	/* On output the 64-byte cache lines hold these data values: */
	/*	0	+32K	2048	+32K	4096	+32K	6144	+32K */
	/*	1	... */
	/*	... */
	/*	2047	... */
	/*	8192	... */
	/*	... */

	/* Do 8192 four_complex_unfft operations */
	/*	distance between fft data elements is 1024 */
	/*	do 2048 macros each processing 32 data values */

	/* Do inverse FFT levels 1,2 */
	/* On input the 64-byte cache lines hold these data values: */
	/*	0	+32K	2048	+32K	4096	+32K	6144	+32K */
	/*	1	... */
	/*	... */
	/*	2047	... */
	/*	8192	... */
	/*	... */
	/* On output the 64-byte cache lines hold these data values if pass 1 swizzles: */
	/*	0	+32K	8192	+32K	16384	+32K	24576	+32K */
	/*	1	... */
	/*	... */
	/*	8191	... */
	/* On output the 64-byte cache lines hold these data values if pass 2 swizzles: */
	/*	0	+8K	16384	+8K	32768	+8K	49152	+8K */
	/*	1	... */
	/*	... */
	/*	8191	... */

	/* Do 8192 four_complex_unfft operations */
	/*	distance between fft data elements is 1024 */
	/*	do 2048 macros each processing 32 data values */

#if 0 // @INSTR(,%xarch,<CORE>) NE 0; /* Core 2 faster doing separate loops */

	start_timer(13);
	rdi = (uintptr_t)g->xsincos_complex;	/* Load sin/cos pointer */
	loops_init_prefetch(2048, 64, 2, rcx);/* 8 iters of 2 iters of 128 */
	r4_x4cl_four_complex_djbunfft_preload;
	for(unsigned int loopA = 8*2; loopA; ) {
		do{
			xtouch(xptr(rcx+4096-64));		/* Load prefetch TLB */
			for(unsigned int loopC = 128; loopC; loopC--) {
				r4_x4cl_four_complex_djbunfft(rsi, 64, dist128, 2*dist128, rdi, 128*16*XMM_SCD2);
				rdi += 16*XMM_SCD2;	/* Next sine/cosine pointer */
			}
			rsi += -128*64+4*dist128;	/* Next source pointer */
			rdi += -128*16*XMM_SCD2;	/* Reset sine/cosine pointer */
		}while(--loopA & 1);
		pfing(rcx += -128*64+8192+128);	/* Next 8KB section to prefetch */
	}
	rsi = daddr; //rsi += -8*2*4*dist128;	/* Restore source pointer */
	end_timer(13);

	start_timer(14);
	loops_init(2048);			/* 8 iters of 2 iters of 128 */
	r4_x4cl_four_complex_djbunfft_preload;
	for(unsigned int loopA = 4; loopA; loopA--) {
		for(unsigned int loopB = 4; loopB; loopB--) {
			for(unsigned int loopC = 128; loopC; loopC--) {
				r4_x4cl_four_complex_djbunfft(rsi, 64, 4*dist128, 8*dist128, rdi, 512*4*XMM_SCD2);
				rdi += 4*XMM_SCD2;		/* Next sine/cosine pointer */
			}
			rsi += -128*64+dist128;	/* Next source pointer */
		}
		rsi += -4*dist128+16*dist128; /* Next source pointer */
		rdi += -512*4*XMM_SCD2;	/* Reset sine/cosine pointer */
	}
	rsi = daddr; //rsi += -4*16*dist128;	/* Restore source pointer */
	end_timer(14);

	start_timer(15);
	rdi = (uintptr_t)g->xsincos_complex;	/* Load sin/cos pointer */
	loops_init(2048);			/* 16 iters of 128 */
	r4_x4cl_four_complex_djbunfft_preload;
	for(unsigned int loopA = 16; loopA; loopA--) {
		for(unsigned int loopB = 128; loopB; loopB--) {
			r4_x4cl_four_complex_djbunfft(rsi, 64, 16*dist128, 32*dist128, rdi, 2048*XMM_SCD2);
			rdi += XMM_SCD2;		/* Next set of sin/cos data */
		}
		rsi += -128*64+dist128;	/* Next source pointer */
	}
	end_timer(15);

#else /* Pentium 4 faster doing merged loops */

	start_timer(13);
	rdi = (uintptr_t)g->xsincos_complex;	/* Load sin/cos pointer */
	rbp = rdi;
	rdx = rdi;
	loops_init_prefetch(2048, 64, 2, rcx); /* 128 iters of 16 */
	r4_x4cl_four_complex_djbunfft_preload;

	for(unsigned int loopA = 8*16; loopA; ) {
		xtouch(xptr(rcx+4096-64));		/* Load prefetch TLB */
		do{
			for(unsigned int loopC = 16; loopC; loopC--) {
				r4_x4cl_four_complex_djbunfft(rsi, 4*dist128, dist128, 2*dist128, rdx, 128*16*XMM_SCD2);
			}
			rsi += -16*4*dist128;	/* Restore source pointer */
			loops_undo(16);

			for(unsigned int loopC = 4; loopC; loopC--) {
				for(unsigned int loopD = 4; loopD; loopD--) {
					r4_x4cl_four_complex_djbunfft(rsi, dist128, 4*dist128, 8*dist128, rbp, 512*4*XMM_SCD2);
					rbp += 128*4*XMM_SCD2;	/* Next sine/cosine pointer */
				}
				rsi += -4*dist128+16*dist128; /* Next source pointer */
				rbp += -4*128*4*XMM_SCD2;	/* Reset sine/cosine pointer */
			}
			rsi += -4*16*dist128;	/* Restore source pointer */
			loops_undo(16);

			for(unsigned int loopC = 16; loopC; loopC--) {
				r4_x4cl_four_complex_djbunfft(rsi, dist128, 16*dist128, 32*dist128, rdi, 2048*XMM_SCD2);
				rdi += 128*XMM_SCD2;	/* Next set of sin/cos data */
			}

			rsi += -16*dist128+64;	/* Next source pointer */
			rdi += -16*128*XMM_SCD2+XMM_SCD2; /* Next set of sin/cos data */
			rbp += 4*XMM_SCD2;		/* Next sine/cosine pointer */
			rdx += 16*XMM_SCD2;	/* Next sine/cosine pointer */
		}while(--loopA & 15);
		pfing(rcx += 128);		/* Next 8KB section to prefetch */
	}
	end_timer(13);

#endif

	xpass2_14_levels_complex_done:;
}while(g->pass2_get_next_block(g) == 0);
}
