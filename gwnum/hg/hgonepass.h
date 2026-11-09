#pragma once
void hgonepass_xmiddle_1(struct gwasm_data *__restrict g);
void hgonepass_xmiddle_1p(struct gwasm_data *__restrict g);
void hgonepass_xmiddle_2(struct gwasm_data *__restrict g);
void hgonepass_xmiddle_2p(struct gwasm_data *__restrict g);
void hgonepass_xmiddle_3(struct gwasm_data *__restrict g);
void hgonepass_xmiddle_3p(struct gwasm_data *__restrict g);
void hgonepass_xmiddle_4(struct gwasm_data *__restrict g);
void hgonepass_xmiddle_4p(struct gwasm_data *__restrict g);

/* Copyright 2001-2024 - Mersenne Research, Inc.  All rights reserved */ \
/* Author:  George Woltman */ \
/* Email: woltman@alum.mit.edu */ \
/* */ \
\
/* ******************************************************** */ \
/* ******************************************************** */ \
/* ********************  FFT MACROS  ********************** */ \
/* ******************************************************** */ \
/* ******************************************************** */ \
\
/* These ffts are all done in a single pass.  On entry the memory layout is */ \
/* quite simple.  Since all FFT data resides in the L2 cache, we can ignore */ \
/* the 128-byte L2 cache lines, we only worry about the 64-byte L1 cache lines. */ \
/* Each 64-byte cache line holds 8 doubles. */ \
\
/* These macros support 4 "types".  Type 1 if forward FFT only. */ \
/* Type 2 is forward FFT, square, inverse FFT.  Type 3 is forward */ \
/* FFT, multiply, inverse FFT.  Type 4 is multiply, and inverse FFT. */ \
\
/* rsi = destination pointer */ \
/* DIST_TO_FFTSRCARG = rbx = distance from rsi to source (zero for types 1, 2) */ \
/* DIST_TO_MULSRCARG = distance from rsi to multiplier (zero for types 2, 3) */ \
\
\
/* Perform a 32-element FFT. */ \
\

#define xfft32(type) \
\
/* Do a multiply with pre-FFTed inputs */ \
\
if(g->ffttype == 4) { \
	xmiddle_4(g); \
}else{ \
	\
	/* Do FFT levels 1,2,3 */ \
	/* Values 0-31 is real data. */ \
	/* */ \
	/* On input the 64-byte cache lines hold these data values: */ \
	/*	0	8	1	9	16	24	17	25 */ \
	/*	2	... */ \
	/*	4	... */ \
	/*	6	... */ \
	/* On output the 64-byte cache lines hold these data values: */ \
	/*	0	1	4	5	8	9	12	13 */ \
	/*	2	... */ \
	/*	16	... */ \
	/*	18	... */ \
	\
	/* Do 4 eight_reals_first_fft macros */ \
	/*	distance between fft data elements is 4 */ \
	/*	do 2 macros each processing 16 data values */ \
	\
	for(unsigned int loop=2; loop != 0; loop--){ \
		s2cl_eight_reals_first_fft(rsi, 64, 2*64); \
	} \
	rsi += -2*64;		/* Restore source pointer */ \
	\
	/* Do FFT levels 4,5 */ \
	/*	do 1 macro each processing 16 data values */ \
	/*	mov	eax, count1		;; 1 iteration */ \
	\
	xmiddle_123(g); \
} \
\
/* Do inverse FFT levels 1,2,3 */ \
/* On input the 64-byte cache lines hold these data values: */ \
/*	0	1	2	3	4	5	6	7 */ \
/*	8	9	10	11	... */ \
/*	... */ \
/* On output the 64-byte cache lines hold these data values: */ \
/*	0	1	8	9	16	17	24	25 */ \
/*	2	... */ \
/*	4	... */ \
/*	6	... */ \
\
/* Do 4 eight_reals_last_unfft macros */ \
/*	distance between fft data elements is 4 */ \
/*	do 1 macros each processing 32 data values */ \
\
x4cl_eight_reals_last_unfft(rsi, 0, 64, 2*64); \
\
xfft_3_ret; \
\
\
\
/* Perform a 32-element negacyclic FFT */ \
\

#define xfft32p(type) \
\
/* Do a multiply with pre-FFTed inputs */ \
\
if(g->ffttype == 4) { \
	xmiddle_4p(g); \
}else{ \
	\
	/* Do FFT levels 1,2 */ \
	/* */ \
	/* On input the 64-byte cache lines hold these data values: */ \
	/*	0	8	1	9	16	24	17	25 */ \
	/*	2	... */ \
	/*	4	... */ \
	/*	6	... */ \
	/* On output the 64-byte cache lines hold these data values: */ \
	/*	0	1	4	5	8	9	12	13 */ \
	/*	2	... */ \
	/*	16	... */ \
	/*	18	... */ \
	\
	/* Do 4 four_complex_first_fft macros */ \
	/*	distance between fft data elements is 4 */ \
	/*	do 2 macros each processing 16 data values */ \
	\
	rdi = (uintptr_t)g->plus1_premults;/* Address of premultiplier table */ \
	for(unsigned int loop=2; loop != 0; loop--){ \
		s2cl_four_complex_first_fft(rsi, 64, 2*64); \
		rdi += 128;		/* Next premultiplier pointer */ \
	} \
	rsi += -2*64;		/* Restore source pointer */ \
	\
	/* Do FFT levels 4,5 */ \
	/*	do 2 macro each processing 16 data values */ \
	/*	mov	eax, count1		;; 2 iterations */ \
	\
	xmiddle_123p(g); \
} \
\
/* Do inverse FFT levels 1,2,3 */ \
/* On input the 64-byte cache lines hold these data values: */ \
/*	0	1	2	3	4	5	6	7 */ \
/*	8	9	10	11	... */ \
/*	... */ \
/* On output the 64-byte cache lines hold these data values: */ \
/*	0	1	8	9	16	17	24	25 */ \
/*	2	... */ \
/*	4	... */ \
/*	6	... */ \
\
/* Do 4 four_complex_last_unfft macros */ \
/*	distance between fft data elements is 4 */ \
/*	do 1 macros each processing 32 data values */ \
\
rdi = (uintptr_t)g->plus1_premults;/* Address of premultiplier table */ \
x4cl_four_complex_last_unfft(rsi, 0, 64, 2*64, 128); \
\
xfft_3_ret; \
\
\
\
/* Perform a 48-element FFT. */ \
\

#define xfft48(type) \
\
/* Do a multiply with pre-FFTed inputs */ \
\
if(g->ffttype == 4) { \
	xmiddle_4(g); \
}else{ \
	\
	/* Do FFT levels 1,2,3 */ \
	/* Values 0-47 is real data. */ \
	/* */ \
	/* On input the 64-byte cache lines hold these data values: */ \
	/*	0	12	1	13	24	36	25	37 */ \
	/*	2 */ \
	/*	... */ \
	/*	10 */ \
	/* On output the 64-byte cache lines hold these data values (in +1 format): */ \
	/*	0	4	8	12 */ \
	/*	2	... */ \
	/*	16	24	32	40 */ \
	/*	18 */ \
	/*	... */ \
	/*	22 */ \
	\
	/* Do 8 six_reals_first_fft macros */ \
	/*	distance between fft data elements is 8 */ \
	/*	do 2 macros each processing 24 data values */ \
	\
	xsix_reals_fft_preload; \
	for(unsigned int loop=2; loop != 0; loop--){ \
		s3cl_six_reals_first_fft(rsi, 64, 2*64); \
	} \
	rsi += -2*64;		/* Restore source pointer */ \
	\
	/* Do FFT levels 4 */ \
	/* Values 0-7 is real data, 8-15 is semi-real data, 16-47 is */ \
	/* complex data. */ \
	/* */ \
	/* On input the 64-byte cache lines hold these data values: */ \
	/*	0	4	8	12 */ \
	/*	2	... */ \
	/*	16	24	32	40 */ \
	/*	18 */ \
	/*	... */ \
	/*	22 */ \
	/* On output the 64-byte cache lines hold these data values: */ \
	/*	0	4	8	12 */ \
	/*	2	... */ \
	/*	16	... */ \
	/*	... */ \
	\
	/* Do 4 four_reals_fft_1 macros */ \
	/*	distance between fft data elements is 4 */ \
	/*	do 1 macros each processing 16 data values */ \
	\
	x2cl_four_reals_fft_1(rsi, 2*64, 64); \
	\
	/* Do 4 two_two_complex_fft macros */ \
	/*	distance between fft data elements is 4 */ \
	/*	do 2 macros each processing 16 data values */ \
	\
	rdi = (uintptr_t)g->sincos2;/* Sine/cosine pointer */ \
	for(unsigned int loop=2; loop != 0; loop--){ \
		x2cl_two_complex_fft(rsi, 64, 2*64, rdi); \
	} \
	rsi += -2*64-2*64;		/* Restore source pointer */ \
	\
	/* Do middle FFT levels */ \
	\
	xmiddle_123(g); \
} \
\
/* Do inverse FFT level 4 */ \
/* On input the 64-byte cache lines hold these data values: */ \
/*	0	1	2	3	4	5	6	7 */ \
/*	8	9	10	11	... */ \
/*	... */ \
/* On output the 64-byte cache lines hold these data values: */ \
/*	0	1	4	5	8	9	12	13 */ \
/*	2	... */ \
/*	16	... */ \
/*	... */ \
\
/* Do 2 eight_reals_unfft_1 macros */ \
/*	distance between fft data elements is 4 */ \
/*	do 1 macro processing 16 data values */ \
\
x2cl_eight_reals_unfft_1(rsi, 2*64, 64); \
\
/* Do 4 two_two_complex_unfft macros */ \
/*	distance between fft data elements is 4 */ \
/*	do 2 macros each processing 16 data values */ \
\
rdi = (uintptr_t)g->sincos2;/* Load sin/cos pointer */ \
for(unsigned int loop=2; loop != 0; loop--){ \
	x2cl_two_complex_unfft(rsi, 2*64, 64); \
	rdi += XMM_SCD;		/* Next sine/cosine pointer */ \
} \
rsi += -6*64;		/* Restore source pointer */ \
\
/* Do inverse FFT levels 1,2,3 */ \
/* On input the 64-byte cache lines hold these data values: */ \
/*	0	4	8	12 */ \
/*	2	... */ \
/*	16	... */ \
/*	... */ \
/* On output the 64-byte cache lines hold these data values: */ \
/*	0	12	24	36 */ \
/*	2	... */ \
/*	... */ \
/*	10 */ \
\
/* Do 8 six_reals_last_unfft macros */ \
/*	distance between fft data elements is 16 */ \
/*	do 2 macros each processing 24 data values */ \
\
xsix_reals_unfft_preload; \
for(unsigned int loop=2; loop != 0; loop--){ \
	x3cl_six_reals_last_unfft(rsi, 64, 2*64); \
} \
\
xfft_3_ret; \
\
\
\
/* Perform a 48-element negacyclic FFT */ \
\

#define xfft48p(type) \
\
/* Do a multiply with pre-FFTed inputs */ \
\
if(g->ffttype == 4) { \
	xmiddle_4p(g); \
}else{ \
	\
	/* Do FFT levels 1,2,3 */ \
	/* */ \
	/* On input the 64-byte cache lines hold these data values: */ \
	/*	0	12	1	13	24	36	25	37 */ \
	/*	2 */ \
	/*	... */ \
	/*	10 */ \
	/* On output the 64-byte cache lines hold these data values (in +1 format): */ \
	/*	0	4	8	12 */ \
	/*	2	... */ \
	/*	16	24	32	40 */ \
	/*	18 */ \
	/*	... */ \
	/*	22 */ \
	\
	/* Do 8 three_complex_first_fft macros */ \
	/*	distance between fft data elements is 8 */ \
	/*	do 2 macros each processing 24 data values */ \
	\
	rdi = (uintptr_t)g->plus1_premults;/* Address of premultiplier table */ \
	for(unsigned int loop=2; loop != 0; loop--){ \
		s3cl_three_complex_first_fft(rsi, 64, 2*64); \
		rdi += 192;		/* Next premultiplier pointer */ \
	} \
	rsi += -2*64;		/* Restore source pointer */ \
	\
	/* Do FFT levels 4 */ \
	/* */ \
	/* On input the 64-byte cache lines hold these data values: */ \
	/*	0	4	8	12 */ \
	/*	2	... */ \
	/*	16	24	32	40 */ \
	/*	18 */ \
	/*	... */ \
	/*	22 */ \
	/* On output the 64-byte cache lines hold these data values: */ \
	/*	0	4	8	12 */ \
	/*	2	... */ \
	/*	16	... */ \
	/*	... */ \
	\
	/* Do 6 two_two_complex_fft macros */ \
	/*	distance between fft data elements is 4 */ \
	/*	do 3 macros each processing 16 data values */ \
	\
	rdi = (uintptr_t)g->sincos2;/* Sine/cosine pointer */ \
	x2cl_two_complex_fft_in_place(rsi, 2*64, 64, rdi); \
	rdi += XMM_SCD;		/* Next sin/cos pointer */ \
	for(unsigned int loop=2; loop != 0; loop--){ \
		x2cl_two_complex_fft(rsi, 64, 2*64, rdi); \
	} \
	rsi += -2*64-2*64;		/* Restore source pointer */ \
	\
	/* Do middle FFT levels */ \
	\
	xmiddle_123p(g); \
} \
\
/* Do inverse FFT level 4 */ \
/* On input the 64-byte cache lines hold these data values: */ \
/*	0	1	2	3	4	5	6	7 */ \
/*	8	9	10	11	... */ \
/*	... */ \
/* On output the 64-byte cache lines hold these data values: */ \
/*	0	1	4	5	8	9	12	13 */ \
/*	2	... */ \
/*	16	... */ \
/*	... */ \
\
/* Do 6 two_two_complex_unfft macros */ \
/*	distance between fft data elements is 4 */ \
/*	do 3 macros each processing 16 data values */ \
\
rdi = (uintptr_t)g->sincos2;/* Load sin/cos pointer */ \
for(unsigned int loop=3; loop != 0; loop--){ \
	x2cl_two_complex_unfft(rsi, 2*64, 64); \
	rdi += XMM_SCD;		/* Next sine/cosine pointer */ \
} \
rsi += -6*64;		/* Restore source pointer */ \
\
/* Do inverse FFT levels 1,2,3 */ \
/* On input the 64-byte cache lines hold these data values: */ \
/*	0	4	8	12 */ \
/*	2	... */ \
/*	16	... */ \
/*	... */ \
/* On output the 64-byte cache lines hold these data values: */ \
/*	0	12	24	36 */ \
/*	2	... */ \
/*	... */ \
/*	10 */ \
\
/* Do 8 three_complex_last_unfft macros */ \
/*	distance between fft data elements is 8 */ \
/*	do 2 macros each processing 24 data values */ \
\
rdi = (uintptr_t)g->plus1_premults;/* Address of premultiplier table */ \
for(unsigned int loop=2; loop != 0; loop--){ \
	x3cl_three_complex_last_unfft(rsi, 64, 2*64); \
	rdi += 192;		/* Next premultiplier pointer */ \
} \
\
xfft_3_ret; \
\
\
\
/* Perform a 64-element FFT. */ \
\

#define xfft64(type) \
\
/* Do a multiply with pre-FFTed inputs */ \
\
if(g->ffttype == 4) { \
	xmiddle_4(g); \
}else{ \
	\
	/* Do FFT levels 1,2,3 */ \
	/* Values 0-63 is real data. */ \
	/* */ \
	/* On input the 64-byte cache lines hold these data values: */ \
	/*	0	16	1	17	32	48	33	49 */ \
	/*	2	... */ \
	/*	... */ \
	/*	14	... */ \
	/* On output the 64-byte cache lines hold these data values: */ \
	/*	0	1	8	9	16	17	24	25 */ \
	/*	2	... */ \
	/*	4	... */ \
	/*	6	... */ \
	/*	32	... */ \
	/*	... */ \
	\
	/* Do 8 eight_reals_first_fft macros */ \
	/*	distance between fft data elements is 8 */ \
	/*	do 4 macros each processing 16 data values */ \
	\
	for(unsigned int loop=4; loop != 0; loop--){ \
		s2cl_eight_reals_first_fft(rsi, 64, 4*64); \
	} \
	rsi += -4*64;		/* Restore source pointer */ \
	\
	/* Do FFT level 4 */ \
	/* Values 0-7 is real data, 8-15 is semi-real data, 16-63 is */ \
	/* complex data. */ \
	/* */ \
	/* On input the 64-byte cache lines hold these data values: */ \
	/*	0	1	8	9	16	17	24	25 */ \
	/*	2	... */ \
	/*	4	... */ \
	/*	6	... */ \
	/*	32	... */ \
	/*	... */ \
	/* On output the 64-byte cache lines hold these data values: */ \
	/*	0	1	4	5	8	9	12	13 */ \
	/*	2	3	... */ \
	/*	16	... */ \
	/*	... */ \
	\
	/* Do 2 eight_reals_fft_1 macros */ \
	/* Do 2 two_two_complex_fft macros */ \
	/*	distance between fft data elements is 4 */ \
	/*	do 2 macros each processing 16 data values */ \
	\
	rdi = (uintptr_t)g->sincos2;/* Load sin/cos pointer */ \
	for(unsigned int loop=2; loop != 0; loop--){ \
		x2cl_eight_reals_fft_1(rsi, 64, 2*64); \
	} \
	\
	/* Do 4 two_two_complex_fft macros */ \
	/*	distance between fft data elements is 4 */ \
	/*	do 2 macros each processing 16 data values */ \
	\
	rsi += -2*64+4*64;		/* Next source pointer */ \
	rdi += XMM_SCD;		/* Next sine/cosine pointer */ \
	for(unsigned int loop=2; loop != 0; loop--){ \
		x2cl_two_complex_fft(rsi, 64, 2*64, rdi); \
	} \
	rsi += -2*64-4*64;		/* Restore source pointer */ \
	\
	/* Do middle FFT levels */ \
	\
	xmiddle_123(g); \
} \
\
/* Do FFT levels 5,6 */ \
/*	do 3 macro each processing 16 data values */ \
/*	mov	eax, count1		;; 3 iterations */ \
\
/* Do inverse FFT level 4 */ \
/* On input the 64-byte cache lines hold these data values: */ \
/*	0	1	2	3	4	5	6	7 */ \
/*	8	9	10	11	... */ \
/*	... */ \
/* On output the 64-byte cache lines hold these data values: */ \
/*	0	1	4	5	8	9	12	13 */ \
/*	2	... */ \
/*	16	... */ \
/*	... */ \
\
/* Do 2 eight_reals_unfft_1 macros */ \
/*	distance between fft data elements is 4 */ \
/*	do 1 macro processing 16 data values */ \
\
x2cl_eight_reals_unfft_1(rsi, 2*64, 64); \
\
/* Do 6 two_two_complex_unfft macros */ \
/*	distance between fft data elements is 4 */ \
/*	do 3 macros each processing 16 data values */ \
\
rdi = (uintptr_t)g->sincos2;/* Load sin/cos pointer */ \
for(unsigned int loop=3; loop != 0; loop--){ \
	x2cl_two_complex_unfft(rsi, 2*64, 64); \
	rdi += XMM_SCD;		/* Next sine/cosine pointer */ \
} \
rsi += -8*64;		/* Restore source pointer */ \
\
/* Do inverse FFT levels 1,2,3 */ \
/* On input the 64-byte cache lines hold these data values: */ \
/*	0	1	4	5	8	9	12	13 */ \
/*	2	... */ \
/*	16	... */ \
/*	... */ \
/* On output the 64-byte cache lines hold these data values: */ \
/*	0	1	16	17	32	33	48	49 */ \
/*	2	... */ \
/*	... */ \
/*	14	... */ \
\
/* Do 8 eight_reals_last_unfft macros */ \
/*	distance between fft data elements is 8 */ \
/*	do 2 macros each processing 32 data values */ \
\
for(unsigned int loop=2; loop != 0; loop--){ \
	x4cl_eight_reals_last_unfft(rsi, 64, 2*64, 4*64); \
} \
\
xfft_3_ret; \
\
\
\
/* Perform a 64-element negacyclic FFT */ \
\

#define xfft64p(type) \
\
/* Do a multiply with pre-FFTed inputs */ \
/* printf("ffttype %i\n", (int)g->ffttype); */ \
if(g->ffttype == 4) { \
	xmiddle_4p(g); \
}else{ \
	\
	/* Do FFT levels 1,2,3 */ \
	/* */ \
	/* On input the 64-byte cache lines hold these data values: */ \
	/*	0	16	1	17	32	48	33	49 */ \
	/*	2	... */ \
	/*	... */ \
	/*	14	... */ \
	/* On output the 64-byte cache lines hold these data values: */ \
	/*	0	1	8	9	16	17	24	25 */ \
	/*	2	... */ \
	/*	4	... */ \
	/*	6	... */ \
	/*	32	... */ \
	/*	... */ \
	\
	/* Do 8 four_complex_first_fft macros */ \
	/*	distance between fft data elements is 8 */ \
	/*	do 4 macros each processing 16 data values */ \
	\
	rdi = (uintptr_t)g->plus1_premults;/* Address of premultiplier table */ \
	for(unsigned int loop=4; loop != 0; loop--){ \
		s2cl_four_complex_first_fft(rsi, 64, 4*64); \
		rdi += 128;		/* Next premultiplier pointer */ \
	} \
	rsi += -4*64;		/* Restore source pointer */ \
	\
	/* Do FFT level 4 */ \
	/* */ \
	/* On input the 64-byte cache lines hold these data values: */ \
	/*	0	1	8	9	16	17	24	25 */ \
	/*	2	... */ \
	/*	4	... */ \
	/*	6	... */ \
	/*	32	... */ \
	/*	... */ \
	/* On output the 64-byte cache lines hold these data values: */ \
	/*	0	1	4	5	8	9	12	13 */ \
	/*	2	3	... */ \
	/*	16	... */ \
	/*	... */ \
	\
	/* Do 8 two_two_complex_fft macros */ \
	/*	distance between fft data elements is 4 */ \
	/*	do 4 macros each processing 16 data values */ \
	\
	rdi = (uintptr_t)g->sincos2;/* Load sin/cos pointer */ \
	for(unsigned int loop=2; loop != 0; loop--){ \
		for(unsigned int loopB=2; loopB != 0; loopB--){ \
			x2cl_two_complex_fft(rsi, 64, 2*64, rdi); \
		} \
		rsi += -2*64+4*64;		/* Restore source pointer */ \
		rdi += 2*XMM_SCD;		/* Next sine/cosine pointer */ \
	} \
	rsi += -8*64;		/* Restore source pointer */ \
	\
	/* End common FFT code */ \
	\
	xmiddle_123p(g); \
} \
\
/* Do FFT levels 5,6 */ \
/*	do 4 macro each processing 16 data values */ \
/*	mov	eax, count1		;; 4 iterations */ \
\
/* Do inverse FFT level 4 */ \
/* On input the 64-byte cache lines hold these data values: */ \
/*	0	1	2	3	4	5	6	7 */ \
/*	8	9	10	11	... */ \
/*	... */ \
/* On output the 64-byte cache lines hold these data values: */ \
/*	0	1	4	5	8	9	12	13 */ \
/*	2	... */ \
/*	16	... */ \
/*	... */ \
\
/* Do 8 two_two_complex_unfft macros */ \
/*	distance between fft data elements is 4 */ \
/*	do 4 macros each processing 16 data values */ \
\
rdi = (uintptr_t)g->sincos2;/* Load sin/cos pointer */ \
for(unsigned int loop=4; loop != 0; loop--){ \
	x2cl_two_complex_unfft(rsi, 2*64, 64); \
	rdi += XMM_SCD;		/* Next sine/cosine pointer */ \
} \
rsi += -8*64;		/* Restore source pointer */ \
\
/* Do inverse FFT levels 1,2,3 */ \
/* On input the 64-byte cache lines hold these data values: */ \
/*	0	1	4	5	8	9	12	13 */ \
/*	2	... */ \
/*	16	... */ \
/*	... */ \
/* On output the 64-byte cache lines hold these data values: */ \
/*	0	1	16	17	32	33	48	49 */ \
/*	2	... */ \
/*	... */ \
/*	14	... */ \
\
/* Do 8 four_complex_last_unfft macros */ \
/*	distance between fft data elements is 8 */ \
/*	do 2 macros each processing 32 data values */ \
\
rdi = (uintptr_t)g->plus1_premults;/* Address of premultiplier table */ \
for(unsigned int loop=2; loop != 0; loop--){ \
	x4cl_four_complex_last_unfft(rsi, 64, 2*64, 4*64, 2*128); \
	rdi += 128;		/* Next premultiplier pointer */ \
} \
\
xfft_3_ret; \
\
\
\
/* Perform a 80-element FFT. */ \
\

#define xfft80(type) \
\
/* Do a multiply with pre-FFTed inputs */ \
\
if(g->ffttype == 4) { \
	xmiddle_4(g); \
}else{ \
	\
	/* Do FFT levels 1,2,3 */ \
	/* Values 0-79 is real data. */ \
	/* */ \
	/* On input the 64-byte cache lines hold these data values: */ \
	/*	0	4	1	5	40	44	41	45 */ \
	/*	2	... */ \
	/*	8	24	9	25	48	64	49	65 */ \
	/*	10 */ \
	/*	... */ \
	/*	22 */ \
	/* On output the 64-byte cache lines hold these data values (in +1 format): */ \
	/*	0	4	8	12 */ \
	/*	2	... */ \
	/*	16	32	48	64 */ \
	/*	18 */ \
	/*	... */ \
	/*	30 */ \
	\
	/* Do 16 five_reals_first_fft macros */ \
	/*	distance between fft data elements is 16 */ \
	/*	do 2 macros each processing 40 data values */ \
	\
	xfive_reals_fft_preload;		/* Preload const1s for five real macros */ \
	for(unsigned int loop=2; loop != 0; loop--){ \
		s5cl_five_reals_first_fft(rsi, 64, 2*64); \
	} \
	rsi += -2*64;		/* Restore source pointer */ \
	\
	/* Do FFT levels 4,5 */ \
	/* Values 0-15 is real data, 16-79 is complex data. */ \
	/* */ \
	/* On input the 64-byte cache lines hold these data values: */ \
	/*	0	4	8	12 */ \
	/*	2	... */ \
	/*	16	32	48	64 */ \
	/*	18 */ \
	/*	... */ \
	/*	30 */ \
	/* On output the 64-byte cache lines hold these data values: */ \
	/*	0	4	8	12 */ \
	/*	2	... */ \
	/*	16	... */ \
	/*	... */ \
	\
	/* Do 2 eight_reals_fft_2 macros */ \
	/*	distance between fft data elements is 4 */ \
	/*	do 2 macros each processing 8 data values */ \
	\
	for(unsigned int loop=2; loop != 0; loop--){ \
		x1cl_half_eight_reals_fft_2(rsi, 64); \
	} \
	\
	/* Do 8 four_complex_fft macros */ \
	/*	distance between fft data elements is 4 */ \
	/*	do 2 macros each processing 32 data values */ \
	\
	rdi = (uintptr_t)g->sincos2; \
	for(unsigned int loop=2; loop != 0; loop--){ \
		x4cl_four_complex_fft(rsi, 64, 2*64, 4*64, rdi); \
	} \
	rsi = (uintptr_t)g->DESTARG;/* Restore source pointer */ \
	\
	/* Do middle FFT levels */ \
	\
	xmiddle_123(g); \
} \
\
/* Do inverse FFT levels 4,5 */ \
/* On input the 64-byte cache lines hold these data values: */ \
/*	0	1	2	3	4	5	6	7 */ \
/*	8	9	10	11	... */ \
/*	... */ \
/* On output the 64-byte cache lines hold these data values: */ \
/*	0	1	8	9	16	17	24	25 */ \
/*	2	... */ \
/*	4	... */ \
/*	6	... */ \
/*	32	... */ \
/*	... */ \
\
/* Do 2 eight_reals_unfft_2 macros */ \
/*	distance between fft data elements is 4 */ \
/*	do 1 macro each processing 16 data values */ \
\
x2cl_half_eight_reals_unfft_2(rsi, 0, 64); \
\
/* Do 8 four_complex_unfft macros */ \
/*	distance between fft data elements is 4 */ \
/*	do 2 macros each processing 32 data values */ \
\
rsi += 2*64;		/* Next source pointer */ \
rdi = (uintptr_t)g->sincos2;/* Load sin/cos pointer */ \
for(unsigned int loop=2; loop != 0; loop--){ \
	x4cl_four_complex_unfft(rsi, 4*64, 64, 2*64, rdi); \
	rdi += XMM_SCD;		/* Next sine/cosine pointer */ \
} \
rsi = (uintptr_t)g->DESTARG;/* Restore source pointer */ \
\
/* Do inverse FFT levels 1,2,3 */ \
/* On output the 64-byte cache lines hold these data values: */ \
/*	0	4	40	44 */ \
/*	2	... */ \
/*	8	24	48	64 */ \
/*	10 */ \
/*	... */ \
/*	22 */ \
\
/* Do 16 five_reals_last_unfft macros */ \
/*	distance between fft data elements is 16 */ \
/*	do 2 macros each processing 40 data values */ \
\
xfive_reals_unfft_preload;	/* Preload const1s for five real macros */ \
for(unsigned int loop=2; loop != 0; loop--){ \
	x5cl_five_reals_last_unfft(rsi, 64, 2*64); \
} \
xfft_3_ret; \
\
\
\
/* Perform a 96-element FFT. */ \
\

#define xfft96(type) \
\
/* Do a multiply with pre-FFTed inputs */ \
\
if(g->ffttype == 4) { \
	xmiddle_4(g); \
}else{ \
	\
	/* Do FFT levels 1,2,3 */ \
	/* Values 0-95 is real data. */ \
	/* */ \
	/* On input the 64-byte cache lines hold these data values: */ \
	/*	0	24	1	25	48	72	49	73 */ \
	/*	2 */ \
	/*	... */ \
	/*	22 */ \
	/* On output the 64-byte cache lines hold these data values (in +1 format): */ \
	/*	0	8	16	24 */ \
	/*	2	... */ \
	/*	... */ \
	/*	6 */ \
	/*	32	48	64	80 */ \
	/*	34 */ \
	/*	... */ \
	/*	46 */ \
	\
	/* Do 16 six_reals_first_fft macros */ \
	/*	distance between fft data elements is 16 */ \
	/*	do 4 macros each processing 24 data values */ \
	\
	xsix_reals_fft_preload; \
	for(unsigned int loop=4; loop != 0; loop--){ \
		s3cl_six_reals_first_fft(rsi, 64, 4*64); \
	} \
	rsi += -4*64;		/* Restore source pointer */ \
	\
	/* Do FFT levels 4,5 */ \
	/* Values 0-15 is real data, 16-31 is semi-real data, 32-95 is */ \
	/* complex data. */ \
	/* */ \
	/* On input the 64-byte cache lines hold these data values: */ \
	/*	0	8	16	24 */ \
	/*	2	... */ \
	/*	... */ \
	/*	6 */ \
	/*	32	48	64	80 */ \
	/*	34 */ \
	/*	... */ \
	/*	46 */ \
	/* On output the 64-byte cache lines hold these data values: */ \
	/*	0	4	8	12 */ \
	/*	2	... */ \
	/*	16	... */ \
	/*	... */ \
	\
	/* Do 2 eight_reals_fft_2 macros */ \
	/* Do 2 nop_two_two_complex_fft_2 macros */ \
	/* Do 4 four_complex_fft macros */ \
	/*	distance between fft data elements is 4 */ \
	/*	do 2 macros each processing 16 data values */ \
	\
	for(unsigned int loop=2; loop != 0; loop--){ \
		x2cl_eight_reals_fft_2(rsi, 64, 2*64); \
	} \
	rsi += -2*64+4*64;		/* Next source pointer */ \
	\
	/* Do 8 four_complex_fft macros */ \
	/*	distance between fft data elements is 4 */ \
	/*	do 2 macros each processing 32 data values */ \
	\
	rdi = (uintptr_t)g->sincos2; \
	for(unsigned int loop=2; loop != 0; loop--){ \
		x4cl_four_complex_fft(rsi, 64, 2*64, 4*64, rdi); \
	} \
	rsi = (uintptr_t)g->DESTARG;/* Restore source pointer */ \
	\
	/* Do middle FFT levels */ \
	\
	xmiddle_123(g); \
} \
\
/* Do inverse FFT levels 4,5 */ \
/* On input the 64-byte cache lines hold these data values: */ \
/*	0	1	2	3	4	5	6	7 */ \
/*	8	9	10	11	... */ \
/*	... */ \
/* On output the 64-byte cache lines hold these data values: */ \
/*	0	1	8	9	16	17	24	25 */ \
/*	2	... */ \
/*	4	... */ \
/*	6	... */ \
/*	32	... */ \
/*	... */ \
\
/* Do 2 eight_reals_unfft_2 macros */ \
/* Do 2 nop_two_two_complex_unfft_2 macros */ \
/*	distance between fft data elements is 4 */ \
/*	do 1 macro each processing 32 data values */ \
\
x4cl_eight_reals_unfft_2(rsi, 0, 64, 2*64); \
\
/* Do 8 four_complex_unfft macros */ \
/*	distance between fft data elements is 4 */ \
/*	do 2 macros each processing 32 data values */ \
\
rsi += 4*64;		/* Next source pointer */ \
rdi = (uintptr_t)g->sincos2;/* Load sin/cos pointer */ \
for(unsigned int loop=2; loop != 0; loop--){ \
	x4cl_four_complex_unfft(rsi, 4*64, 64, 2*64, rdi); \
	rdi += XMM_SCD;		/* Next sine/cosine pointer */ \
} \
rsi = (uintptr_t)g->DESTARG;/* Next source pointer */ \
\
/* Do inverse FFT levels 1,2,3 */ \
/* On input the 64-byte cache lines hold these data values: */ \
/*	0	8	16	24 */ \
/*	2	... */ \
/*	4	... */ \
/*	6	... */ \
/*	32	... */ \
/*	... */ \
/* On output the 64-byte cache lines hold these data values: */ \
/*	0	24	48	72 */ \
/*	2	... */ \
/*	... */ \
/*	22 */ \
\
/* Do 16 six_reals_last_unfft macros */ \
/*	distance between fft data elements is 16 */ \
/*	do 4 macros each processing 24 data values */ \
\
xsix_reals_unfft_preload; \
for(unsigned int loop=4; loop != 0; loop--){ \
	x3cl_six_reals_last_unfft(rsi, 64, 4*64); \
} \
xfft_3_ret; \
\
\
\
/* Perform a 96-element negacyclic FFT */ \
\

#define xfft96p(type) \
\
/* Do a multiply with pre-FFTed inputs */ \
/* printf("ffttype %i\n", (int)g->ffttype); */ \
\
if(g->ffttype == 4) { \
	xmiddle_4p(g); \
}else{ \
	\
	/* Do FFT levels 1,2,3 */ \
	/* */ \
	/* On input the 64-byte cache lines hold these data values: */ \
	/*	0	24	1	25	48	72	49	73 */ \
	/*	2 */ \
	/*	... */ \
	/*	22 */ \
	/* On output the 64-byte cache lines hold these data values (in +1 format): */ \
	/*	0	8	16	24 */ \
	/*	2	... */ \
	/*	... */ \
	/*	6 */ \
	/*	32	48	64	80 */ \
	/*	34 */ \
	/*	... */ \
	/*	46 */ \
	\
	/* Do 16 three_complex_first_fft macros */ \
	/*	distance between fft data elements is 16 */ \
	/*	do 4 macros each processing 24 data values */ \
	\
	rdi = (uintptr_t)g->plus1_premults;/* Address of premultiplier table */ \
	for(unsigned int loop=4; loop != 0; loop--){ \
		s3cl_three_complex_first_fft(rsi, 64, 4*64); \
		rdi += 192;		/* Next premultiplier pointer */ \
	} \
	rsi += -4*64;		/* Restore source pointer */ \
	\
	/* Do FFT levels 4,5 */ \
	/* */ \
	/* On input the 64-byte cache lines hold these data values: */ \
	/*	0	8	16	24 */ \
	/*	2	... */ \
	/*	... */ \
	/*	6 */ \
	/*	32	48	64	80 */ \
	/*	34 */ \
	/*	... */ \
	/*	46 */ \
	/* On output the 64-byte cache lines hold these data values: */ \
	/*	0	4	8	12 */ \
	/*	2	... */ \
	/*	16	... */ \
	/*	... */ \
	\
	/* Do 12 four_complex_fft macros */ \
	/*	distance between fft data elements is 4 */ \
	/*	do 2 macros each processing 16 data values */ \
	/*	do 2 macros each processing 32 data values */ \
	\
	rdi = (uintptr_t)g->sincos2; \
	for(unsigned int loop=2; loop != 0; loop--){ \
		x2cl_four_complex_fft(rsi, 64, 2*64); \
	} \
	rsi += -2*64+4*64;		/* Next source pointer */ \
	rdi += XMM_SCD;		/* Next sin/cos pointer */ \
	for(unsigned int loop=2; loop != 0; loop--){ \
		x4cl_four_complex_fft(rsi, 64, 2*64, 4*64, rdi); \
	} \
	rsi += -2*64-4*64;		/* Restore source pointer */ \
	\
	/* Do common FFT code */ \
	\
	xmiddle_123p(g); \
} \
\
/* Do FFT levels 6,7 */ \
/*	do 8 macros each processing 16 data values */ \
/*	mov	eax, count1		;; 8 iterations */ \
\
/* Do inverse FFT levels 4,5 */ \
/* On input the 64-byte cache lines hold these data values: */ \
/*	0	1	2	3	4	5	6	7 */ \
/*	8	9	10	11	... */ \
/*	... */ \
/* On output the 64-byte cache lines hold these data values: */ \
/*	0	1	8	9	16	17	24	25 */ \
/*	2	... */ \
/*	4	... */ \
/*	6	... */ \
/*	32	... */ \
/*	... */ \
\
/* Do 12 four_complex_unfft macros */ \
/*	distance between fft data elements is 4 */ \
/*	do 3 macros each processing 32 data values */ \
\
rdi = (uintptr_t)g->sincos2;/* Load sin/cos pointer */ \
for(unsigned int loop=3; loop != 0; loop--){ \
	x4cl_four_complex_unfft(rsi, 4*64, 64, 2*64, rdi); \
	rdi += XMM_SCD;		/* Next sine/cosine pointer */ \
} \
rsi = (uintptr_t)g->DESTARG;/* Next source pointer */ \
\
/* Do inverse FFT levels 1,2,3 */ \
/* On input the 64-byte cache lines hold these data values: */ \
/*	0	8	16	24 */ \
/*	2	... */ \
/*	4	... */ \
/*	6	... */ \
/*	32	... */ \
/*	... */ \
/* On output the 64-byte cache lines hold these data values: */ \
/*	0	24	48	72 */ \
/*	2	... */ \
/*	... */ \
/*	22 */ \
\
/* Do 16 three_complex_last_unfft macros */ \
/*	distance between fft data elements is 16 */ \
/*	do 4 macros each processing 24 data values */ \
\
rdi = (uintptr_t)g->plus1_premults;/* Address of premultiplier table */ \
for(unsigned int loop=4; loop != 0; loop--){ \
	x3cl_three_complex_last_unfft(rsi, 64, 4*64); \
	rdi += 192;		/* Next premultiplier pointer */ \
} \
\
xfft_3_ret; \
\
\
\
/* Perform a 112-element FFT. */ \
\

#define xfft112(type) \
\
/* Do a multiply with pre-FFTed inputs */ \
\
if(g->ffttype == 4) { \
	xmiddle_4(g); \
}else{ \
	\
	/* Do FFT levels 1,2,3 */ \
	/* Values 0-111 is real data. */ \
	/* */ \
	/* On input the 64-byte cache lines hold these data values: */ \
	/*	0	4	1	5	56	60	57	61 */ \
	/*	2	... */ \
	/*	8	16	9	17	64	72	65	73 */ \
	/*	10 */ \
	/*	... */ \
	/*	14 */ \
	/*	24	40	25	41	80	96	81	97 */ \
	/*	26 */ \
	/*	... */ \
	/*	38 */ \
	/* On output the 64-byte cache lines hold these data values (in +1 format): */ \
	/*	0	4	8	12 */ \
	/*	2	... */ \
	/*	16	24	32	40 */ \
	/*	18 */ \
	/*	... */ \
	/*	22 */ \
	/*	48	64	80	96 */ \
	/*	50 */ \
	/*	... */ \
	/*	62 */ \
	\
	/* Do 16 seven_reals_first_fft macros */ \
	/*	distance between fft data elements is 16 */ \
	/*	do 2 macros each processing 56 data values */ \
	\
	xseven_reals_fft_preload; \
	for(unsigned int loop=2; loop != 0; loop--){ \
		s7cl_seven_reals_first_fft(rsi, 64, 2*64); \
	} \
	rsi += -2*64;		/* Restore source pointer */ \
	\
	/* Do FFT levels 4,5 */ \
	/* Values 0-15 is real data, 16-111 is complex data. */ \
	\
	/* Do 2 eight_reals_fft_2 macros */ \
	/*	distance between fft data elements is 4 */ \
	/*	do 2 macros each processing 8 data values */ \
	\
	for(unsigned int loop=2; loop != 0; loop--){ \
		x1cl_half_eight_reals_fft_2(rsi, 64); \
	} \
	\
	/* Do 4 four_complex_fft macros */ \
	/*	distance between fft data elements is 4 */ \
	/*	do 2 macros each processing 16 data values */ \
	\
	rdi = (uintptr_t)g->sincos2;/* Load sin/cos pointer */ \
	for(unsigned int loop=2; loop != 0; loop--){ \
		x2cl_four_complex_fft(rsi, 64, 2*64); \
	} \
	rsi += -2*64+4*64;		/* Next source pointer */ \
	\
	/* Do 8 four_complex_fft macros */ \
	/*	distance between fft data elements is 4 */ \
	/*	do 2 macros each processing 32 data values */ \
	\
	rdi += XMM_SCD;		/* Next sine/cosine pointer */ \
	for(unsigned int loop=2; loop != 0; loop--){ \
		x4cl_four_complex_fft(rsi, 64, 2*64, 4*64, rdi); \
	} \
	rsi = (uintptr_t)g->DESTARG;/* Restore source pointer */ \
	\
	/* Do middle FFT levels */ \
	\
	xmiddle_123(g); \
} \
\
/* Do inverse FFT levels 4,5 */ \
/* On input the 64-byte cache lines hold these data values: */ \
/*	0	1	2	3	4	5	6	7 */ \
/*	8	9	10	11	... */ \
/*	... */ \
/* On output the 64-byte cache lines hold these data values: */ \
/*	0	1	8	9	16	17	24	25 */ \
/*	2	... */ \
/*	4	... */ \
/*	6	... */ \
/*	32	... */ \
/*	... */ \
\
/* Do 2 eight_reals_unfft_2 macros */ \
/*	distance between fft data elements is 4 */ \
/*	do 1 macro each processing 16 data values */ \
\
x2cl_half_eight_reals_unfft_2(rsi, 0, 64); \
\
/* Do 12 four_complex_unfft macros */ \
/*	distance between fft data elements is 4 */ \
/*	do 3 macros each processing 32 data values */ \
\
rsi += 2*64;		/* Next source pointer */ \
rdi = (uintptr_t)g->sincos2;/* Load sin/cos pointer */ \
for(unsigned int loop=3; loop != 0; loop--){ \
	x4cl_four_complex_unfft(rsi, 4*64, 64, 2*64, rdi); \
	rdi += XMM_SCD;		/* Next sine/cosine pointer */ \
} \
rsi = (uintptr_t)g->DESTARG;/* Restore source pointer */ \
\
/* Do inverse FFT levels 1,2,3 */ \
/* On input the 64-byte cache lines hold these data values: */ \
/*	0	4	8	12 */ \
/*	2	... */ \
/*	16	24	32	40 */ \
/*	... */ \
/* On output the 64-byte cache lines hold these data values: */ \
/*	0	4	56	60 */ \
/*	2	... */ \
/*	8	16	64	72 */ \
/*	10 */ \
/*	... */ \
/*	14 */ \
/*	24	40	80	96 */ \
/*	26 */ \
/*	... */ \
/*	38 */ \
\
/* Do 16 seven_reals_last_unfft macros */ \
/*	distance between fft data elements is 16 */ \
/*	do 2 macros each processing 56 data values */ \
\
xseven_reals_unfft_preload; \
for(unsigned int loop=2; loop != 0; loop--){ \
	x7cl_seven_reals_last_unfft(rsi, 64, 2*64); \
} \
\
xfft_3_ret; \
\
\
\
/* Perform a 128-element FFT. */ \
\

#define xfft128(type) \
\
/* Do a multiply with pre-FFTed inputs */ \
\
if(g->ffttype == 4) { \
	xmiddle_4(g); \
}else{ \
	\
	/* Do FFT levels 1,2,3 */ \
	/* Values 0-127 is real data. */ \
	/* */ \
	/* On input the 64-byte cache lines hold these data values: */ \
	/*	0	32	1	33	64	96	65	97 */ \
	/*	2	... */ \
	/*	... */ \
	/*	30 */ \
	/* On output the 64-byte cache lines hold these data values: */ \
	/*	0	1	16	17	32	33	48	49 */ \
	/*	2	... */ \
	/*	... */ \
	/*	14	... */ \
	/*	64	... */ \
	/*	... */ \
	\
	/* Do 16 eight_reals_first_fft macros */ \
	/*	distance between fft data elements is 16 */ \
	/*	do 8 macros each processing 16 data values */ \
	\
	for(unsigned int loop=8; loop != 0; loop--){ \
		s2cl_eight_reals_first_fft(rsi, 64, 8*64); \
	} \
	rsi += -8*64;		/* Restore source pointer */ \
	\
	/* Do FFT levels 4,5 */ \
	/* Values 0-15 is real data, 16-31 is semi-real data, 32-128 is */ \
	/* complex data. */ \
	/* */ \
	/* On input the 64-byte cache lines hold these data values: */ \
	/*	0	1	16	17	32	33	48	49 */ \
	/*	2	... */ \
	/*	... */ \
	/*	14	... */ \
	/*	64	... */ \
	/*	... */ \
	/* On output the 64-byte cache lines hold these data values: */ \
	/*	0	1	4	5	8	9	12	13 */ \
	/*	2	... */ \
	/*	16	... */ \
	/*	... */ \
	\
	/* Do 2 eight_reals_fft_2 macros */ \
	/* Do 2 nop_two_two_complex_fft_2 macros */ \
	/* Do 4 four_complex_fft macros */ \
	/*	distance between fft data elements is 4 */ \
	/*	do 2 macros each processing 32 data values */ \
	\
	rdi = (uintptr_t)g->sincos2;/* Load sin/cos pointer */ \
	for(unsigned int loop=2; loop != 0; loop--){ \
		x4cl_eight_reals_fft_2(rsi, 64, 2*64, 4*64, rdi); \
	} \
	rsi += -2*64+8*64;		/* Next source pointer */ \
	\
	/* Do 8 four_complex_fft macros */ \
	/*	distance between fft data elements is 4 */ \
	/*	do 2 macros each processing 32 data values */ \
	\
	rdi += XMM_SCD;		/* Next sine/cosine pointer */ \
	for(unsigned int loop=2; loop != 0; loop--){ \
		x4cl_four_complex_fft(rsi, 64, 2*64, 4*64, rdi); \
	} \
	rsi = (uintptr_t)g->DESTARG;/* Restore source pointer */ \
	\
	/* Do middle FFT levels */ \
	/* Do FFT levels 6,7 */ \
	/*	do 4/5/6/7 macros each processing 16 data values */ \
	/*	mov	eax, count1		;; 4/5/6/7 iterations */ \
	\
	xmiddle_123(g); \
} \
\
/* Do inverse FFT levels 4,5 */ \
/* On input the 64-byte cache lines hold these data values: */ \
/*	0	1	2	3	4	5	6	7 */ \
/*	8	9	10	11	... */ \
/*	... */ \
/* On output the 64-byte cache lines hold these data values: */ \
/*	0	1	8	9	16	17	24	25 */ \
/*	2	... */ \
/*	4	... */ \
/*	6	... */ \
/*	32	... */ \
/*	... */ \
\
/* Do 2 eight_reals_unfft_2 macros */ \
/* Do 2 nop_two_two_complex_unfft_2 macros */ \
/*	distance between fft data elements is 4 */ \
/*	do 1 macros each processing 32 data values */ \
\
x4cl_eight_reals_unfft_2(rsi, 0, 64, 2*64); \
\
/* Do 12 four_complex_unfft macros */ \
/*	distance between fft data elements is 4 */ \
/*	do 3 macros each processing 32 data values */ \
\
rsi += 4*64;		/* Next source pointer */ \
rdi = (uintptr_t)g->sincos2;/* Load sin/cos pointer */ \
for(unsigned int loop=3; loop != 0; loop--){ \
	x4cl_four_complex_unfft(rsi, 4*64, 64, 2*64, rdi); \
	rdi += XMM_SCD;		/* Next sine/cosine pointer */ \
} \
rsi = (uintptr_t)g->DESTARG;/* Restore source pointer */ \
\
/* Do inverse FFT levels 1,2,3 */ \
/* On input the 64-byte cache lines hold these data values: */ \
/*	0	1	8	9	16	17	24	25 */ \
/*	2	... */ \
/*	4	... */ \
/*	6	... */ \
/*	32	... */ \
/*	... */ \
/* On output the 64-byte cache lines hold these data values: */ \
/*	0	1	32	33	64	65	96	97 */ \
/*	2	... */ \
/*	... */ \
/*	30 */ \
\
/* Do 16 eight_reals_last_unfft macros */ \
/*	distance between fft data elements is 16 */ \
/*	do 4 macros each processing 32 data values */ \
\
for(unsigned int loop=4; loop != 0; loop--){ \
	x4cl_eight_reals_last_unfft(rsi, 64, 4*64, 8*64); \
} \
\
xfft_3_ret; \
\
\
\
/* Perform a 128-element negacyclic FFT */ \
\

#define xfft128p(type) \
\
/* Do a multiply with pre-FFTed inputs */ \
\
if(g->ffttype == 4) { \
	xmiddle_4p(g); \
}else{ \
	\
	/* Do FFT levels 1,2,3 */ \
	/* */ \
	/* On input the 64-byte cache lines hold these data values: */ \
	/*	0	32	1	33	64	96	65	97 */ \
	/*	2	... */ \
	/*	... */ \
	/*	30 */ \
	/* On output the 64-byte cache lines hold these data values: */ \
	/*	0	1	16	17	32	33	48	49 */ \
	/*	2	... */ \
	/*	... */ \
	/*	14	... */ \
	/*	64	... */ \
	/*	... */ \
	\
	/* Do 16 four_complex_first_fft macros */ \
	/*	distance between fft data elements is 16 */ \
	/*	do 8 macros each processing 16 data values */ \
	\
	rdi = (uintptr_t)g->plus1_premults;/* Address of premultiplier table */ \
	for(unsigned int loop=8; loop != 0; loop--){ \
		s2cl_four_complex_first_fft(rsi, 64, 8*64); \
		rdi += 128;		/* Next premultiplier pointer */ \
	} \
	rsi += -8*64;		/* Restore source pointer */ \
	\
	/* Do FFT levels 4,5 */ \
	/* */ \
	/* On input the 64-byte cache lines hold these data values: */ \
	/*	0	1	16	17	32	33	48	49 */ \
	/*	2	... */ \
	/*	... */ \
	/*	14	... */ \
	/*	64	... */ \
	/*	... */ \
	/* On output the 64-byte cache lines hold these data values: */ \
	/*	0	1	4	5	8	9	12	13 */ \
	/*	2	... */ \
	/*	16	... */ \
	/*	... */ \
	\
	/* Do 8 four_complex_fft macros */ \
	/*	distance between fft data elements is 4 */ \
	/*	do 2 macros each processing 32 data values */ \
	\
	/* Do FFT levels 6,7 */ \
	/*	do 8 macros each processing 16 data values */ \
	/*	mov	eax, count1		;; 8 iterations */ \
	\
	/* Join common code to perform rest of the FFT */ \
	\
	xlast_two_fft_levels_p; \
} \
\
/* Do inverse FFT levels 4,5 */ \
/* On input the 64-byte cache lines hold these data values: */ \
/*	0	1	2	3	4	5	6	7 */ \
/*	8	9	10	11	... */ \
/*	... */ \
/* On output the 64-byte cache lines hold these data values: */ \
/*	0	1	8	9	16	17	24	25 */ \
/*	2	... */ \
/*	4	... */ \
/*	6	... */ \
/*	32	... */ \
/*	... */ \
\
/* Do 16 four_complex_unfft macros */ \
/*	distance between fft data elements is 4 */ \
/*	do 4 macros each processing 32 data values */ \
\
rdi = (uintptr_t)g->sincos2;/* Load sin/cos pointer */ \
for(unsigned int loop=4; loop != 0; loop--){ \
	x4cl_four_complex_unfft(rsi, 4*64, 64, 2*64, rdi); \
	rdi += XMM_SCD;		/* Next sine/cosine pointer */ \
} \
rsi = (uintptr_t)g->DESTARG;/* Next source pointer */ \
\
/* Do inverse FFT levels 1,2,3 */ \
/* On input the 64-byte cache lines hold these data values: */ \
/*	0	1	8	9	16	17	24	25 */ \
/*	2	... */ \
/*	4	... */ \
/*	6	... */ \
/*	32	... */ \
/*	... */ \
/* On output the 64-byte cache lines hold these data values: */ \
/*	0	1	32	33	64	65	96	97 */ \
/*	2	... */ \
/*	... */ \
/*	30 */ \
\
/* Do 16 four_complex_last_unfft macros */ \
/*	distance between fft data elements is 16 */ \
/*	do 4 macros each processing 32 data values */ \
\
rdi = (uintptr_t)g->plus1_premults;/* Address of premultiplier table */ \
for(unsigned int loop=4; loop != 0; loop--){ \
	x4cl_four_complex_last_unfft(rsi, 64, 4*64, 8*64, 4*128); \
	rdi += 128;		/* Next premultiplier pointer */ \
} \
\
xfft_3_ret; \
\
\
\
/* Perform a 160-element FFT. */ \
\

#define xfft160(type) \
\
/* Do a multiply with pre-FFTed inputs */ \
\
if(g->ffttype == 4) { \
	xmiddle_4(g); \
}else{ \
	\
	/* Do FFT levels 1,2,3 */ \
	/* Values 0-159 is real data. */ \
	/* */ \
	/* On input the 64-byte cache lines hold these data values: */ \
	/*	0	8	1	9	80	88	81	89 */ \
	/*	2	... */ \
	/*	... */ \
	/*	6 */ \
	/*	16	48	17	49	96	128	97	129 */ \
	/*	18 */ \
	/*	... */ \
	/*	46 */ \
	/* On output the 64-byte cache lines hold these data values (in +1 format): */ \
	/*	0	8	16	24 */ \
	/*	2	... */ \
	/*	... */ \
	/*	6 */ \
	/*	32	64	96	128 */ \
	/*	34 */ \
	/*	... */ \
	/*	62 */ \
	\
	/* Do 32 five_reals_first_fft macros */ \
	/*	distance between fft data elements is 32 */ \
	/*	do 4 macros each processing 40 data values */ \
	\
	xfive_reals_fft_preload;		/* Preload const1s for five real macros */ \
	for(unsigned int loop=4; loop != 0; loop--){ \
		s5cl_five_reals_first_fft(rsi, 64, 4*64); \
	} \
	rsi += -4*64;		/* Restore source pointer */ \
	\
	/* Do FFT levels 4,5 */ \
	/* Values 0-31 is real data, 32-159 is complex data. */ \
	\
	/* Do 4 eight_reals_fft_2 macros */ \
	/*	distance between fft data elements is 8 */ \
	/*	do 4 macros each processing 8 data values */ \
	\
	for(unsigned int loop=4; loop != 0; loop--){ \
		x1cl_half_eight_reals_fft_2(rsi, 64); \
	} \
	\
	/* Do 8 four_complex_fft macros */ \
	/*	distance between fft data elements is 8 */ \
	/*	do 4 macros each processing 32 data values */ \
	\
	rdi = (uintptr_t)g->sincos3;		/* Load sin/cos pointer */ \
	for(unsigned int loop=4; loop != 0; loop--){ \
		x4cl_four_complex_fft(rsi, 64, 4*64, 8*64, rdi); \
	} \
	\
	/* Join common code to perform rest of the FFT */ \
	\
	cmn256_2; \
} \
/*rcx = 4*64;*/ \
cmn256(4*64); \
\
/* Finish the inverse FFT levels 4,5 */ \
\
/* Do 4 eight_reals_unfft_2 macros */ \
/*	distance between fft data elements is 8 */ \
/*	do 2 macros each processing 16 data values */ \
\
for(unsigned int loop=2; loop != 0; loop--){ \
	x2cl_half_eight_reals_unfft_2(rsi, 64, 2*64); \
} \
rsi += -2*64;		/* Restore source pointer */ \
\
/* Do inverse FFT levels 1,2,3 */ \
/* On output the 64-byte cache lines hold these data values: */ \
/*	0	8	80	88 */ \
/*	2	... */ \
/*	... */ \
/*	6 */ \
/*	16	48	96	128 */ \
/*	18 */ \
/*	... */ \
/*	46 */ \
\
/* Do 32 five_reals_last_unfft macros */ \
/*	distance between fft data elements is 32 */ \
/*	do 4 macros each processing 40 data values */ \
\
xfive_reals_unfft_preload;	/* Preload const1s for five real macros */ \
for(unsigned int loop=4; loop != 0; loop--){ \
	x5cl_five_reals_last_unfft(rsi, 64, 4*64); \
} \
\
xfft_3_ret; \
\
\
\
/* Perform a 192-element FFT. */ \
\

#define xfft192(type) \
\
/* Do a multiply with pre-FFTed inputs */ \
\
if(g->ffttype == 4) { \
	xmiddle_4(g); \
}else{ \
	\
	/* Do FFT levels 1,2,3 */ \
	/* Values 0-191 is real data. */ \
	/* */ \
	/* On input the 64-byte cache lines hold these data values: */ \
	/*	0	48	1	49	96	144	97	145 */ \
	/*	2 */ \
	/*	... */ \
	/*	46 */ \
	/* On output the 64-byte cache lines hold these data values (in +1 format): */ \
	/*	0	16	32	48 */ \
	/*	2	.. */ \
	/*	... */ \
	/*	14 */ \
	/*	64	96	128	160 */ \
	/*	66 */ \
	/*	... */ \
	/*	94 */ \
	\
	/* Do 32 six_reals_first_fft macros */ \
	/*	distance between fft data elements is 32 */ \
	/*	do 8 macros each processing 24 data values */ \
	\
	xsix_reals_fft_preload; \
	for(unsigned int loop=8; loop != 0; loop--){ \
		s3cl_six_reals_first_fft(rsi, 64, 8*64); \
	} \
	rsi += -8*64;		/* Restore source pointer */ \
	\
	/* Do FFT levels 4,5 */ \
	/* Values 0-31 is real data, 32-63 is semi-real data, 64-191 is */ \
	/* complex data. */ \
	/* */ \
	/* On input the 64-byte cache lines hold these data values: */ \
	/*	0	16	32	48 */ \
	/*	2	.. */ \
	/*	... */ \
	/*	14 */ \
	/*	64	96	128	160 */ \
	/*	66 */ \
	/*	... */ \
	/*	94 */ \
	/* On output the 64-byte cache lines hold these data values: */ \
	/*	0	8	16	24 */ \
	/*	2	... */ \
	/*	... */ \
	/*	6	... */ \
	/*	32	... */ \
	/*	... */ \
	\
	/* Do 4 eight_reals_fft_2 macros */ \
	/* Do 4 nop_two_two_complex_fft_2 macros */ \
	/* Do 8 four_complex_fft macros */ \
	/*	distance between fft data elements is 8 */ \
	/*	do 4 macros each processing 16 data values */ \
	\
	for(unsigned int loop=4; loop != 0; loop--){ \
		x2cl_eight_reals_fft_2(rsi, 64, 4*64); \
	} \
	rsi += -4*64+8*64;		/* Next source pointer */ \
	\
	/* Do 8 four_complex_fft macros */ \
	/*	distance between fft data elements is 8 */ \
	/*	do 4 macros each processing 32 data values */ \
	\
	rdi = (uintptr_t)g->sincos3;		/* Load sin/cos pointer */ \
	for(unsigned int loop=4; loop != 0; loop--){ \
		x4cl_four_complex_fft(rsi, 64, 4*64, 8*64, rdi); \
	} \
	\
	/* Join common code to perform rest of the FFT */ \
	\
	cmn256_2; \
} \
/*rcx = 8*64;*/ \
cmn256(8*64); \
\
/* Finish inverse FFT levels 4,5 */ \
\
/* Do 4 eight_reals_unfft_2 macros */ \
/* Do 4 nop_two_two_complex_unfft_2 macros */ \
/*	distance between fft data elements is 8 */ \
/*	do 2 macros each processing 32 data values */ \
\
for(unsigned int loop=2; loop != 0; loop--){ \
	x4cl_eight_reals_unfft_2(rsi, 64, 2*64, 4*64); \
} \
rsi += -2*64;		/* Restore source pointer */ \
\
/* Do inverse FFT levels 1,2,3 */ \
/* On input the 64-byte cache lines hold these data values: */ \
/*	0	1	16	17	32	33	48	49 */ \
/*	2	... */ \
/*	...	... */ \
/*	14	... */ \
/*	64	... */ \
/*	... */ \
/* On output the 64-byte cache lines hold these data values: */ \
/*	0	48	96	144 */ \
/*	2 */ \
/*	... */ \
/*	46 */ \
\
/* Do 32 six_reals_last_unfft macros */ \
/*	distance between fft data elements is 32 */ \
/*	do 8 macros each processing 24 data values */ \
\
xsix_reals_unfft_preload; \
for(unsigned int loop=8; loop != 0; loop--){ \
	x3cl_six_reals_last_unfft(rsi, 64, 8*64); \
} \
\
xfft_3_ret; \
\
\
\
/* Perform a 192-element negacyclic FFT */ \
\

#define xfft192p(type) \
\
/* Do a multiply with pre-FFTed inputs */ \
\
if(g->ffttype == 4) { \
	xmiddle_4p(g); \
}else{ \
	\
	/* Do FFT levels 1,2,3 */ \
	/* */ \
	/* On input the 64-byte cache lines hold these data values: */ \
	/*	0	48	1	49	96	144	97	145 */ \
	/*	2 */ \
	/*	... */ \
	/*	46 */ \
	/* On output the 64-byte cache lines hold these data values (in +1 format): */ \
	/*	0	16	32	48 */ \
	/*	2	.. */ \
	/*	... */ \
	/*	14 */ \
	/*	64	96	128	160 */ \
	/*	66 */ \
	/*	... */ \
	/*	94 */ \
	\
	/* Do 32 three_complex_first_fft macros */ \
	/*	distance between fft data elements is 32 */ \
	/*	do 8 macros each processing 24 data values */ \
	\
	rdi = (uintptr_t)g->plus1_premults;/* Address of premultiplier table */ \
	for(unsigned int loop=8; loop != 0; loop--){ \
		s3cl_three_complex_first_fft(rsi, 64, 8*64); \
		rdi += 192;		/* Next premultiplier pointer */ \
	} \
	rsi += -8*64;		/* Restore source pointer */ \
	\
	/* Do FFT levels 4 */ \
	/* */ \
	/* On input the 64-byte cache lines hold these data values: */ \
	/*	0	16	32	48 */ \
	/*	2	.. */ \
	/*	... */ \
	/*	14 */ \
	/*	64	96	128	160 */ \
	/*	66 */ \
	/*	... */ \
	/*	94 */ \
	/* On output the 64-byte cache lines hold these data values: */ \
	/*	0	16	32	48 */ \
	/*	2	... */ \
	/*	... */ \
	/*	14 */ \
	/*	64	... */ \
	/*	... */ \
	\
	/* Do 24 two_two_complex_fft macros */ \
	/*	distance between fft data elements is 16 */ \
	/*	do 12 macros each processing 16 data values */ \
	\
	rdi = (uintptr_t)g->sincos3;		/* Load sin/cos pointer */ \
	for(unsigned int loop=4; loop != 0; loop--){ \
		x2cl_two_complex_fft_in_place(rsi, 64, 4*64, rdi); \
	} \
	rsi += -4*64+8*64;		/* Next source pointer */ \
	rdi += XMM_SCD;		/* Next sin/cos pointer */ \
	for(unsigned int loop=8; loop != 0; loop--){ \
		x2cl_two_complex_fft(rsi, 64, 8*64, rdi); \
	} \
	rsi += -8*64-8*64;		/* Restore source pointer */ \
	\
	/* Do FFT levels 5,6 */ \
	\
	/* Join common code to perform rest of the FFT */ \
	\
	xlast_two_fft_levels_p; \
} \
cmn256p; \
\
/* Do inverse FFT levels 1,2,3 */ \
/* On input the 64-byte cache lines hold these data values: */ \
/*	0	1	16	17	32	33	48	49 */ \
/*	2	... */ \
/*	...	... */ \
/*	14	... */ \
/*	64	... */ \
/*	... */ \
/* On output the 64-byte cache lines hold these data values: */ \
/*	0	48	96	144 */ \
/*	2 */ \
/*	... */ \
/*	46 */ \
\
/* Do 32 three_complex_last_unfft macros */ \
/*	distance between fft data elements is 32 */ \
/*	do 8 macros each processing 24 data values */ \
\
rdi = (uintptr_t)g->plus1_premults;/* Address of premultiplier table */ \
for(unsigned int loop=8; loop != 0; loop--){ \
	x3cl_three_complex_last_unfft(rsi, 64, 8*64); \
	rdi += 192;		/* Next premultiplier pointer */ \
} \
\
xfft_3_ret; \
\
\
\
/* Perform a 224-element FFT. */ \
\

#define xfft224(type) \
\
/* Do a multiply with pre-FFTed inputs */ \
\
if(g->ffttype == 4) { \
	xmiddle_4(g); \
}else{ \
	\
	/* Do FFT levels 1,2,3 */ \
	/* Values 0-223 is real data. */ \
	/* */ \
	/* On input the 64-byte cache lines hold these data values: */ \
	/*	0	8	1	9	112	120	113	121 */ \
	/*	2	... */ \
	/*	... */ \
	/*	6 */ \
	/*	16	32	17	33	128	144	129	145 */ \
	/*	18 */ \
	/*	... */ \
	/*	30 */ \
	/*	48	80	49	81	160	192	161	193 */ \
	/*	50 */ \
	/*	... */ \
	/*	78 */ \
	/* On output the 64-byte cache lines hold these data values (in +1 format): */ \
	/*	0	8	16	24 */ \
	/*	2	... */ \
	/*	... */ \
	/*	6 */ \
	/*	32	48	64	80 */ \
	/*	34 */ \
	/*	... */ \
	/*	46 */ \
	/*	96	128	160	192 */ \
	/*	98 */ \
	/*	... */ \
	/*	126 */ \
	\
	/* Do 32 seven_reals_first_fft macros */ \
	/*	distance between fft data elements is 32 */ \
	/*	do 4 macros each processing 56 data values */ \
	\
	xseven_reals_fft_preload; \
	for(unsigned int loop=4; loop != 0; loop--){ \
		s7cl_seven_reals_first_fft(rsi, 64, 4*64); \
	} \
	rsi += -4*64;		/* Restore source pointer */ \
	\
	/* Do FFT levels 4,5 */ \
	/* Values 0-31 is real data, 32-223 is complex data. */ \
	\
	/* Do 4 eight_reals_fft_2 macros */ \
	/*	distance between fft data elements is 8 */ \
	/*	do 4 macros each processing 8 data values */ \
	\
	for(unsigned int loop=4; loop != 0; loop--){ \
		x1cl_half_eight_reals_fft_2(rsi, 64); \
	} \
	\
	/* Do 8 four_complex_fft macros */ \
	/*	distance between fft data elements is 8 */ \
	/*	do 4 macros each processing 16 data values */ \
	\
	rdi = (uintptr_t)g->sincos3;		/* Load sin/cos pointer */ \
	for(unsigned int loop=4; loop != 0; loop--){ \
		x2cl_four_complex_fft(rsi, 64, 4*64); \
	} \
	rsi += -4*64+8*64;		/* Next source pointer */ \
	\
	/* Do 16 four_complex_fft macros */ \
	/*	distance between fft data elements is 8 */ \
	/*	do 4 macros each processing 32 data values */ \
	\
	rdi += XMM_SCD;		/* Next sine/cosine pointer */ \
	for(unsigned int loop=4; loop != 0; loop--){ \
		x4cl_four_complex_fft(rsi, 64, 4*64, 8*64, rdi); \
	} \
	\
	/* Join common code to perform rest of the FFT */ \
	\
	cmn256_2; \
} \
/*rcx = 4*64;*/ \
cmn256(4*64); \
\
/* Finish inverse FFT levels 4,5 */ \
\
/* Do 4 eight_reals_unfft_2 macros */ \
/*	distance between fft data elements is 8 */ \
/*	do 2 macros each processing 16 data values */ \
\
for(unsigned int loop=2; loop != 0; loop--){ \
	x2cl_half_eight_reals_unfft_2(rsi, 64, 2*64); \
} \
rsi += -2*64;		/* Restore source pointer */ \
\
/* Do inverse FFT levels 1,2,3 */ \
/*	0	8	16	24 */ \
/*	2	... */ \
/*	... */ \
/*	6 */ \
/*	32	48	64	80 */ \
/*	... */ \
/* On output the 64-byte cache lines hold these data values: */ \
/*	0	8	112	120 */ \
/*	2	... */ \
/*	... */ \
/*	6 */ \
/*	16	32	128	144 */ \
/*	18 */ \
/*	... */ \
/*	30 */ \
/*	48	80	160	192 */ \
/*	50 */ \
/*	... */ \
/*	78 */ \
\
/* Do 32 seven_reals_last_unfft macros */ \
/*	distance between fft data elements is 32 */ \
/*	do 4 macros each processing 56 data values */ \
\
xseven_reals_unfft_preload; \
for(unsigned int loop=4; loop != 0; loop--){ \
	x7cl_seven_reals_last_unfft(rsi, 64, 4*64); \
} \
\
xfft_3_ret; \
\
\
\
/* Perform a 256-element FFT. */ \
\

#define xfft256(type) \
\
/* Do a multiply with pre-FFTed inputs */ \
\
if(g->ffttype == 4) { \
	xmiddle_4(g); \
}else{ \
	\
	/* Do FFT levels 1,2,3 */ \
	/* Values 0-255 is real data. */ \
	/* */ \
	/* On input the 64-byte cache lines hold these data values: */ \
	/*	0	64	1	65	128	192	129	193 */ \
	/*	2	... */ \
	/*	... */ \
	/*	62 */ \
	/* On output the 64-byte cache lines hold these data values: */ \
	/*	0	1	32	33	64	65	96	97 */ \
	/*	2	... */ \
	/*	... */ \
	/*	30	... */ \
	/*	128	... */ \
	/*	... */ \
	\
	/* Do 32 eight_reals_first_fft macros */ \
	/*	distance between fft data elements is 32 */ \
	/*	do 16 macros each processing 16 data values */ \
	\
	for(unsigned int loop=16; loop != 0; loop--){ \
		s2cl_eight_reals_first_fft(rsi, 64, 16*64); \
	} \
	rsi += -16*64;		/* Restore source pointer */ \
	\
	/* Do FFT levels 4,5 */ \
	/* Values 0-31 is real data, 32-63 is semi-real data, 64-255 is */ \
	/* complex data. */ \
	/* */ \
	/* On input the 64-byte cache lines hold these data values: */ \
	/*	0	1	32	33	64	65	96	97 */ \
	/*	2	... */ \
	/*	... */ \
	/*	30	... */ \
	/*	128	... */ \
	/*	... */ \
	/* On output the 64-byte cache lines hold these data values: */ \
	/*	0	1	8	9	16	17	24	25 */ \
	/*	2	... */ \
	/*	... */ \
	/*	6	... */ \
	/*	32	... */ \
	/*	... */ \
	\
	/* Do 4 eight_reals_fft_2 macros */ \
	/* Do 4 nop_two_two_complex_fft_2 macros */ \
	/* Do 8 four_complex_fft macros */ \
	/*	distance between fft data elements is 8 */ \
	/*	do 4 macros each processing 32 data values */ \
	\
	rdi = (uintptr_t)g->sincos3;		/* Load sin/cos pointer */ \
	for(unsigned int loop=4; loop != 0; loop--){ \
		x4cl_eight_reals_fft_2(rsi, 64, 4*64, 8*64, rdi); \
	} \
	rsi += -4*64+16*64;	/* Next source pointer */ \
	\
	/* Do 16 four_complex_fft macros */ \
	/*	distance between fft data elements is 8 */ \
	/*	do 4 macros each processing 32 data values */ \
	\
	rdi += XMM_SCD;		/* Next sine/cosine pointer */ \
	for(unsigned int loop=4; loop != 0; loop--){ \
		x4cl_four_complex_fft(rsi, 64, 4*64, 8*64, rdi); \
	} \
	\
	/* Call common code for rest of the FFT */ \
	\
	cmn256_2; \
} \
/*rcx = 8*64;*/ \
cmn256(8*64); \
\
/* Finish the inverse FFT levels 4,5 */ \
\
/* Do 4 eight_reals_unfft_2 macros */ \
/* Do 4 nop_two_two_complex_unfft_2 macros */ \
/*	distance between fft data elements is 8 */ \
/*	do 2 macros each processing 32 data values */ \
\
for(unsigned int loop=2; loop != 0; loop--){ \
	x4cl_eight_reals_unfft_2(rsi, 64, 2*64, 4*64); \
} \
rsi += -2*64;		/* Restore source pointer */ \
\
/* Do inverse FFT levels 1,2,3 */ \
/* On input the 64-byte cache lines hold these data values: */ \
/*	0	1	16	17	32	33	48	49 */ \
/*	2	... */ \
/*	...	... */ \
/*	14	... */ \
/*	64	... */ \
/*	... */ \
/* On output the 64-byte cache lines hold these data values: */ \
/*	0	1	64	65	128	129	192	193 */ \
/*	2	... */ \
/*	... */ \
/*	62 */ \
\
/* Do 32 eight_reals_last_unfft macros */ \
/*	distance between fft data elements is 32 */ \
/*	do 8 macros each processing 32 data values */ \
\
for(unsigned int loop=8; loop != 0; loop--){ \
	x4cl_eight_reals_last_unfft(rsi, 64, 8*64, 16*64); \
} \
\
xfft_3_ret; \


/* Common code for length 160, 192, 224, and 256 FFTs */

#define cmn256_2 \
/*push_amt = SZPTR;*/ \
\
/* Do FFT level 6 */ \
/* Values 0-7 is real data, 8-15 is semi-real data, 16-255 is */ \
/* complex data. */ \
/* */ \
/* On input the 64-byte cache lines hold these data values: */ \
/*	0	1	8	9	16	17	24	25 */ \
/*	2	... */ \
/*	... */ \
/*	6	... */ \
/*	32	... */ \
/*	... */ \
/* On output the 64-byte cache lines hold these data values: */ \
/*	0	1	4	5	8	9	12	13 */ \
/*	2	... */ \
/*	16	... */ \
/*	... */ \
\
/* Do 2 eight_reals_fft_1 macros */ \
/* Do 2 two_two_complex_fft macros */ \
/*	distance between fft data elements is 4 */ \
/*	do 2 macros each processing 16 data values */ \
\
rsi = (uintptr_t)g->DESTARG;/* Next source pointer */ \
rdi = (uintptr_t)g->sincos2;/* Load sin/cos pointer */ \
for(unsigned int loop=2; loop != 0; loop--){ \
	x2cl_eight_reals_fft_1(rsi, 64, 2*64); \
} \
\
/* Do 16/20/24/28 two_two_complex_fft macros */ \
/*	distance between fft data elements is 4 */ \
/*	do 8/10/12/14 macros each processing 16 data values */ \
\
rsi += -2*64+4*64;		/* Load source pointer */ \
rdi += XMM_SCD;		/* Next sine/cosine pointer */ \
for(unsigned int loop = g->count5 & 0xFF; loop != 0; loop--){ \
	for(unsigned int loop=2; loop != 0; loop--){ \
		x2cl_two_complex_fft(rsi, 64, 2*64, rdi); \
	} \
	rsi += -2*64+4*64;		/* Next source pointer */ \
	rdi += 2*XMM_SCD;		/* Next sine/cosine pointer */ \
} \
rsi = (uintptr_t)g->DESTARG;/* Restore source pointer */ \
\
/* Do FFT levels 7,8 */ \
/*	do 9/11/13/15 macros each processing 16 data values */ \
/*	mov	eax, count1		;; 1/3 then 8/12 iterations */ \
xmiddle_123(g); /* jmp xmiddle_123 */


/* Do inverse FFT level 6 */ \
/* On input the 64-byte cache lines hold these data values: */ \
/*	0	1	2	3	4	5	6	7 */ \
/*	8	9	10	11	... */ \
/*	... */ \
/* On output the 64-byte cache lines hold these data values: */ \
/*	0	1	4	5	8	9	12	13 */ \
/*	2	... */ \
/*	16	... */ \
/*	... */ \

/* Do 2 eight_reals_unfft_1 macros */ \
/*	distance between fft data elements is 4 */ \
/*	do 1 macros each processing 16 data values */ \

#define cmn256(RCX1) \
x2cl_eight_reals_unfft_1(rsi, 2*64, 64); \
/*push_amt = SZPTR;*/ \
\
/* Do 18/22/26/30 four_complex_unfft macros */ \
/*	distance between fft data elements is 4 */ \
/*	do 9/11/13/15 macros each processing 16 data values */ \
\
rdi = (uintptr_t)g->sincos2;/* Load sin/cos pointer */ \
for(unsigned int loop = g->count2 & 0xFF; loop != 0; loop--){ \
	x2cl_two_complex_unfft(rsi, 2*64, 64); \
	rdi += XMM_SCD;		/* Next sine/cosine pointer */ \
} \
rsi = (uintptr_t)g->DESTARG;/* Restore source pointer */ \
\
/* Do inverse FFT levels 4,5 */ \
/* On input the 64-byte cache lines hold these data values: */ \
/*	0	1	4	5	8	9	12	13 */ \
/*	2	... */ \
/*	16	... */ \
/*	... */ \
/* On output the 64-byte cache lines hold these data values: */ \
/*	0	1	16	17	32	33	48	49 */ \
/*	2	... */ \
/*	...	... */ \
/*	14	... */ \
/*	64	... */ \
/*	... */ \
\
/* Do 16/24 four_complex_unfft macros */ \
/*	distance between fft data elements is 8 */ \
/*	do 4/6 macros each processing 32 data values */ \
\
rsi += RCX1;		/* Next source pointer */ \
rdi = (uintptr_t)g->sincos3;		/* Load sin/cos pointer */ \
for(unsigned int loop = g->count3 & 0xFF; loop != 0; loop--){ \
	for(unsigned int loopB=2; loopB != 0; loopB--){ \
		x4cl_four_complex_unfft(rsi, 64, 2*64, 4*64, rdi); \
	} \
	rsi += -2*64+8*64;		/* Next source pointer */ \
	rdi += XMM_SCD;		/* Next sine/cosine pointer */ \
} \
rsi = (uintptr_t)g->DESTARG;/* Next source pointer */ \
\
/* This ends the common inverse FFT code */ \
\
/* ret;   push_amt = 0; */



/* Perform a 256-element negacyclic FFT */ \


#define xfft256p(type) \
\
/* Do a multiply with pre-FFTed inputs */ \
\
if(g->ffttype == 4) { \
	xmiddle_4p(g); \
}else{ \
	\
	/* Do FFT levels 1,2,3 */ \
	/* */ \
	/* On input the 64-byte cache lines hold these data values: */ \
	/*	0	64	1	65	128	192	129	193 */ \
	/*	2	... */ \
	/*	... */ \
	/*	62 */ \
	/* On output the 64-byte cache lines hold these data values: */ \
	/*	0	1	32	33	64	65	96	97 */ \
	/*	2	... */ \
	/*	... */ \
	/*	30	... */ \
	/*	128	... */ \
	/*	... */ \
	\
	/* Do 32 four_complex_first_fft macros */ \
	/*	distance between fft data elements is 32 */ \
	/*	do 16 macros each processing 16 data values */ \
	\
	rdi = (uintptr_t)g->plus1_premults;/* Address of premultiplier table */ \
	for(unsigned int loop=16; loop != 0; loop--){ \
		s2cl_four_complex_first_fft(rsi, 64, 16*64); \
		rdi += 128;		/* Next premultiplier pointer */ \
	} \
	rsi += -16*64;		/* Restore source pointer */ \
	\
	/* Do FFT level 4 */ \
	/* */ \
	/* On input the 64-byte cache lines hold these data values: */ \
	/*	0	1	32	33	64	65	96	97 */ \
	/*	2	... */ \
	/*	... */ \
	/*	30	... */ \
	/*	128	... */ \
	/*	... */ \
	/* On output the 64-byte cache lines hold these data values: */ \
	/*	0	1	16	17	32	33	48	49 */ \
	/*	2	3	... */ \
	/*	... */ \
	/*	14 */ \
	/*	64	... */ \
	/*	... */ \
	\
	/* Do 32 two_two_complex_fft macros */ \
	/*	distance between fft data elements is 16 */ \
	/*	do 16 macros each processing 16 data values */ \
	\
	rdi = (uintptr_t)g->sincos3;		/* Load sin/cos pointer */ \
	for(unsigned int loop=2; loop != 0; loop--){ \
		for(unsigned int loopB=8; loopB != 0; loopB--){ \
			x2cl_two_complex_fft(rsi, 64, 8*64, rdi); \
		} \
		rsi += -8*64+16*64;	/* Restore source pointer */ \
		rdi += 2*XMM_SCD;		/* Next sine/cosine pointer */ \
	} \
	rsi += -2*16*64;		/* Restore source pointer */ \
	\
	/* Do FFT levels 5,6 */ \
	\
	/* Join common code to perform rest of the FFT */ \
	/* Do FFT levels 7,8 */ \
	/*	do 12/16 macros each processing 16 data values */ \
	\
	xlast_two_fft_levels_p; \
} \
cmn256p; \
\
/* Do inverse FFT levels 1,2,3 */ \
/* On input the 64-byte cache lines hold these data values: */ \
/*	0	1	16	17	32	33	48	49 */ \
/*	2	... */ \
/*	...	... */ \
/*	14	... */ \
/*	64	... */ \
/*	... */ \
/* On output the 64-byte cache lines hold these data values: */ \
/*	0	1	64	65	128	129	192	193 */ \
/*	2	... */ \
/*	... */ \
/*	62 */ \
\
/* Do 32 four_complex_last_unfft macros */ \
/*	distance between fft data elements is 32 */ \
/*	do 8 macros each processing 32 data values */ \
\
rdi = (uintptr_t)g->plus1_premults;/* Address of premultiplier table */ \
for(unsigned int loop=8; loop != 0; loop--){ \
	x4cl_four_complex_last_unfft(rsi, 64, 8*64, 16*64, 8*128); \
	rdi += 128;		/* Next premultiplier pointer */ \
} \
\
xfft_3_ret; \

/* Do inverse FFT levels 5,6 */ \
/* On input the 64-byte cache lines hold these data values: */ \
/*	0	1	2	3	4	5	6	7 */ \
/*	8	9	10	11	... */ \
/*	... */ \
/* On output the 64-byte cache lines hold these data values: */ \
/*	0	1	8	9	16	17	24	25 */ \
/*	2	... */ \
/*	4	... */ \
/*	6	... */ \
/*	32	... */ \
/*	... */ \

/* Do 32 four_complex_unfft macros */ \
/*	distance between fft data elements is 4 */ \
/*	do 8 macros each processing 32 data values */ \

#define cmn256p \
/*push_amt = SZPTR;*/ \
rdi = (uintptr_t)g->sincos2;/* Load sin/cos pointer */ \
for(unsigned int loop = g->count2 & 0xFF; loop != 0; loop--){ \
	x4cl_four_complex_unfft(rsi, 4*64, 64, 2*64, rdi); \
	rdi += XMM_SCD;		/* Next sine/cosine pointer */ \
} \
rsi = (uintptr_t)g->DESTARG;/* Restore source pointer */ \
\
/* Do inverse FFT level 4 */ \
/* On input the 64-byte cache lines hold these data values: */ \
/*	0	1	8	9	16	17	24	25 */ \
/*	2	... */ \
/*	4	... */ \
/*	6	... */ \
/*	32	... */ \
/*	... */ \
/* On output the 64-byte cache lines hold these data values: */ \
/*	0	1	16	17	32	33	48	49 */ \
/*	2	... */ \
/*	...	... */ \
/*	14	... */ \
/*	64	... */ \
/*	... */ \
\
/* Do 24/32 two_two_complex_unfft macros */ \
/*	distance between fft data elements is 16 */ \
/*	do 12/16 macros each processing 16 data values */ \
\
rdi = (uintptr_t)g->sincos3;		/* Load sin/cos pointer */ \
for(unsigned int loop = g->count3 & 0xFF; loop != 0; loop--){ \
	for(unsigned int loop=4; loop != 0; loop--){ \
		x2cl_two_complex_unfft(rsi, 64, 4*64); \
	} \
	rsi += -4*64+8*64;		/* Next source pointer */ \
	rdi += XMM_SCD;		/* Next sine/cosine pointer */ \
} \
rsi = (uintptr_t)g->DESTARG;/* Restore source pointer */ \
\
/* This ends the common inverse FFT code */ \
\
/* ret; push_amt = 0; */



/* Perform a 320-element FFT. */


#define xfft320(type) \
\
/* Do a multiply with pre-FFTed inputs */ \
\
if(g->ffttype == 4) { \
	xmiddle_4(g); \
}else{ \
	\
	/* Do FFT levels 1,2,3 */ \
	/* Values 0-319 is real data. */ \
	/* */ \
	/* On input the 64-byte cache lines hold these data values: */ \
	/*	0	16	1	17	160	176	161	175 */ \
	/*	2	... */ \
	/*	... */ \
	/*	14 */ \
	/*	32	96	33	97	192	256	193	257 */ \
	/*	34 */ \
	/*	... */ \
	/*	94 */ \
	/* On output the 64-byte cache lines hold these data values (in +1 format): */ \
	/*	0	16	32	48 */ \
	/*	2	... */ \
	/*	... */ \
	/*	14 */ \
	/*	64	128	192	256 */ \
	/*	66 */ \
	/*	... */ \
	/*	126 */ \
	\
	/* Do 64 five_reals_first_fft macros */ \
	/*	distance between fft data elements is 64 */ \
	/*	do 8 macros each processing 40 data values */ \
	\
	xfive_reals_fft_preload;		/* Preload const1s for five real macros */ \
	for(unsigned int loop=8; loop != 0; loop--){ \
		s5cl_five_reals_first_fft(rsi, 64, 8*64); \
	} \
	rsi += -8*64;		/* Restore source pointer */ \
	\
	/* Do FFT levels 4,5 */ \
	/* Values 0-63 is real data, 64-319 is complex data. */ \
	\
	/* Do 8 eight_reals_fft_2 macros */ \
	/*	distance between fft data elements is 16 */ \
	/*	do 8 macros each processing 8 data values */ \
	\
	for(unsigned int loop=8; loop != 0; loop--){ \
		x1cl_half_eight_reals_fft_2(rsi, 64); \
	} \
	\
	/* Do 16 four_complex_fft macros */ \
	/*	distance between fft data elements is 16 */ \
	/*	do 8 macros each processing 32 data values */ \
	\
	rdi = (uintptr_t)g->sincos3;		/* Load sin/cos pointer */ \
	for(unsigned int loop=8; loop != 0; loop--){ \
		x4cl_four_complex_fft(rsi, 64, 8*64, 16*64, rdi); \
	} \
	\
	/* Join common code to perform rest of the FFT */ \
	\
	xlast_two_fft_levels; \
} \
/*rcx = 8*64;*/ \
cmn512(8*64); \
\
/* Finish the inverse FFT levels 4,5 */ \
\
/* Do 8 eight_reals_unfft_2 macros */ \
/*	distance between fft data elements is 16 */ \
/*	do 4 macros each processing 16 data values */ \
\
for(unsigned int loop=4; loop != 0; loop--){ \
	x2cl_half_eight_reals_unfft_2(rsi, 64, 4*64); \
} \
rsi += -4*64;		/* Restore source pointer */ \
\
/* Do inverse FFT levels 1,2,3 */ \
/* On output the 64-byte cache lines hold these data values: */ \
/*	0	16	160	176 */ \
/*	2	... */ \
/*	... */ \
/*	14 */ \
/*	32	96	192	256 */ \
/*	34 */ \
/*	... */ \
/*	92 */ \
\
/* Do 64 five_reals_last_unfft macros */ \
/*	distance between fft data elements is 64 */ \
/*	do 8 macros each processing 40 data values */ \
\
xfive_reals_unfft_preload;	/* Preload const1s for five real macros */ \
for(unsigned int loop=8; loop != 0; loop--){ \
	x5cl_five_reals_last_unfft(rsi, 64, 8*64); \
} \
\
xfft_3_ret; \
\
\
\
/* Perform a 384-element FFT. */ \
\

#define xfft384(type) \
\
/* Do a multiply with pre-FFTed inputs */ \
\
if(g->ffttype == 4) { \
	xmiddle_4(g); \
}else{ \
	\
	/* Do FFT levels 1,2,3 */ \
	/* Values 0-383 is real data. */ \
	/* */ \
	/* On input the 64-byte cache lines hold these data values: */ \
	/*	0	96	1	97	192	288	193	289 */ \
	/*	2 */ \
	/*	... */ \
	/*	94 */ \
	/* On output the 64-byte cache lines hold these data values (in +1 format): */ \
	/*	0	32	64	96 */ \
	/*	2	.. */ \
	/*	... */ \
	/*	30 */ \
	/*	128	192	256	320 */ \
	/*	130 */ \
	/*	... */ \
	/*	190 */ \
	\
	/* Do 64 six_reals_first_fft macros */ \
	/*	distance between fft data elements is 64 */ \
	/*	do 16 macros each processing 24 data values */ \
	\
	xsix_reals_fft_preload; \
	for(unsigned int loop=16; loop != 0; loop--){ \
		s3cl_six_reals_first_fft(rsi, 64, 16*64); \
	} \
	rsi += -16*64;		/* Restore source pointer */ \
	\
	/* Do FFT levels 4,5 */ \
	/* Values 0-63 is real data, 64-127 is semi-real data, 128-383 is */ \
	/* complex data. */ \
	/* */ \
	/* On input the 64-byte cache lines hold these data values: */ \
	/*	0	32	64	96 */ \
	/*	2	.. */ \
	/*	... */ \
	/*	30 */ \
	/*	128	192	256	320 */ \
	/*	130 */ \
	/*	... */ \
	/*	190 */ \
	/* On output the 64-byte cache lines hold these data values: */ \
	/*	0	16	32	48 */ \
	/*	2	... */ \
	/*	... */ \
	/*	14	... */ \
	/*	64	... */ \
	/*	... */ \
	\
	/* Do 8 eight_reals_fft_2 macros */ \
	/* Do 8 nop_two_two_complex_fft_2 macros */ \
	/* Do 16 four_complex_fft macros */ \
	/*	distance between fft data elements is 16 */ \
	/*	do 8 macros each processing 16 data values */ \
	\
	for(unsigned int loop=8; loop != 0; loop--){ \
		x2cl_eight_reals_fft_2(rsi, 64, 8*64); \
	} \
	rsi += -8*64+16*64;	/* Next source pointer */ \
	\
	/* Do 16 four_complex_fft macros */ \
	/*	distance between fft data elements is 16 */ \
	/*	do 8 macros each processing 32 data values */ \
	\
	rdi = (uintptr_t)g->sincos3;		/* Load sin/cos pointer */ \
	for(unsigned int loop=8; loop != 0; loop--){ \
		x4cl_four_complex_fft(rsi, 64, 8*64, 16*64, rdi); \
	} \
	\
	/* Join common code to perform rest of the FFT */ \
	\
	xlast_two_fft_levels; \
} \
/*rcx = 16*64*/; \
cmn512(16*64); \
\
/* Finish inverse FFT levels 4,5 */ \
\
/* Do 8 eight_reals_unfft_2 macros */ \
/* Do 8 nop_two_two_complex_unfft_2 macros */ \
/*	distance between fft data elements is 16 */ \
/*	do 4 macros each processing 32 data values */ \
\
for(unsigned int loop=4; loop != 0; loop--){ \
	x4cl_eight_reals_unfft_2(rsi, 64, 4*64, 8*64); \
} \
rsi += -4*64;		/* Restore source pointer */ \
\
/* Do inverse FFT levels 1,2,3 */ \
/* On input the 64-byte cache lines hold these data values: */ \
/*	0	1	32	33	64	65	96	97 */ \
/*	2	... */ \
/*	...	... */ \
/*	30	... */ \
/*	128	... */ \
/*	... */ \
/* On output the 64-byte cache lines hold these data values: */ \
/*	0	96	192	288 */ \
/*	2 */ \
/*	... */ \
/*	94 */ \
\
/* Do 64 six_reals_last_unfft macros */ \
/*	distance between fft data elements is 64 */ \
/*	do 16 macros each processing 24 data values */ \
\
xsix_reals_unfft_preload; \
for(unsigned int loop=16; loop != 0; loop--){ \
	x3cl_six_reals_last_unfft(rsi, 64, 16*64); \
} \
xfft_3_ret; \
\
\
\
/* Perform a 384-element negacyclic FFT */ \
\

#define xfft384p(type) \
\
/* Do a multiply with pre-FFTed inputs */ \
\
if(g->ffttype == 4) { \
	xmiddle_4p(g); \
}else{ \
	\
	/* Do FFT levels 1,2,3 */ \
	/* */ \
	/* On input the 64-byte cache lines hold these data values: */ \
	/*	0	96	1	97	192	288	193	289 */ \
	/*	2 */ \
	/*	... */ \
	/*	94 */ \
	/* On output the 64-byte cache lines hold these data values (in +1 format): */ \
	/*	0	32	64	96 */ \
	/*	2	.. */ \
	/*	... */ \
	/*	30 */ \
	/*	128	192	256	320 */ \
	/*	130 */ \
	/*	... */ \
	/*	190 */ \
	\
	/* Do 64 three_complex_first_fft macros */ \
	/*	distance between fft data elements is 64 */ \
	/*	do 16 macros each processing 24 data values */ \
	\
	rdi = (uintptr_t)g->plus1_premults;/* Address of premultiplier table */ \
	for(unsigned int loop=16; loop != 0; loop--){ \
		s3cl_three_complex_first_fft(rsi, 64, 16*64); \
		rdi += 192;		/* Next premultiplier pointer */ \
	} \
	rsi += -16*64;		/* Restore source pointer */ \
	\
	/* Do FFT levels 4,5 */ \
	/* */ \
	/* On input the 64-byte cache lines hold these data values: */ \
	/*	0	32	64	96 */ \
	/*	2	.. */ \
	/*	... */ \
	/*	30 */ \
	/*	128	192	256	320 */ \
	/*	130 */ \
	/*	... */ \
	/*	190 */ \
	/* On output the 64-byte cache lines hold these data values: */ \
	/*	0	16	32	48 */ \
	/*	2	... */ \
	/*	... */ \
	/*	14	... */ \
	/*	64	... */ \
	/*	... */ \
	\
	/* Do 24 four_complex_fft macros */ \
	/*	distance between fft data elements is 16 */ \
	/*	do 8 macros each processing 16 data values */ \
	/*	do 8 macros each processing 32 data values */ \
	\
	rdi = (uintptr_t)g->sincos3;		/* Load sin/cos pointer */ \
	for(unsigned int loop=8; loop != 0; loop--){ \
		x2cl_four_complex_fft(rsi, 64, 8*64); \
	} \
	rsi += -8*64+16*64;	/* Next source pointer */ \
	rdi += XMM_SCD;		/* Next sin/cos pointer */ \
	for(unsigned int loop=8; loop != 0; loop--){ \
		x4cl_four_complex_fft(rsi, 64, 8*64, 16*64, rdi); \
	} \
	rsi += -8*64-16*64;	/* Restore source pointer */ \
	\
	/* Join common code to perform rest of the FFT */ \
	\
	xlast_two_fft_levels_p; \
} \
cmn512p; \
\
/* Do inverse FFT levels 1,2,3 */ \
/* On input the 64-byte cache lines hold these data values: */ \
/*	0	1	32	33	64	65	96	97 */ \
/*	2	... */ \
/*	...	... */ \
/*	30	... */ \
/*	128	... */ \
/*	... */ \
/* On output the 64-byte cache lines hold these data values: */ \
/*	0	96	192	288 */ \
/*	2 */ \
/*	... */ \
/*	94 */ \
\
/* Do 64 three_complex_last_unfft macros */ \
/*	distance between fft data elements is 64 */ \
/*	do 16 macros each processing 24 data values */ \
\
rdi = (uintptr_t)g->plus1_premults;/* Address of premultiplier table */ \
for(unsigned int loop=16; loop != 0; loop--){ \
	x3cl_three_complex_last_unfft(rsi, 64, 16*64); \
	rdi += 192;		/* Next premultiplier pointer */ \
} \
\
xfft_3_ret; \
\
\
\
/* Perform a 448-element FFT. */ \
\

#define xfft448(type) \
\
/* Do a multiply with pre-FFTed inputs */ \
\
if(g->ffttype == 4) { \
	xmiddle_4(g); \
}else{ \
	\
	/* Do FFT levels 1,2,3 */ \
	/* Values 0-447 is real data. */ \
	/* */ \
	/* On input the 64-byte cache lines hold these data values: */ \
	/*	0	16	1	17	224	240	225	241 */ \
	/*	2	... */ \
	/*	... */ \
	/*	14 */ \
	/*	32	64	33	65	256	288	257	289 */ \
	/*	34 */ \
	/*	... */ \
	/*	62 */ \
	/*	96	160	97	161	320	384	321	385 */ \
	/*	98 */ \
	/*	... */ \
	/*	158 */ \
	/* On output the 64-byte cache lines hold these data values (in +1 format): */ \
	/*	0	16	32	48 */ \
	/*	2	... */ \
	/*	... */ \
	/*	14 */ \
	/*	64	96	128	160 */ \
	/*	66 */ \
	/*	... */ \
	/*	94 */ \
	/*	192	256	320	384 */ \
	/*	194 */ \
	/*	... */ \
	/*	254 */ \
	\
	/* Do 64 seven_reals_first_fft macros */ \
	/*	distance between fft data elements is 64 */ \
	/*	do 8 macros each processing 56 data values */ \
	\
	xseven_reals_fft_preload; \
	for(unsigned int loop=8; loop != 0; loop--){ \
		s7cl_seven_reals_first_fft(rsi, 64, 8*64); \
	} \
	rsi += -8*64;		/* Restore source pointer */ \
	\
	/* Do FFT levels 4,5 */ \
	/* Values 0-63 is real data, 64-447 is complex data. */ \
	\
	/* Do 8 eight_reals_fft_2 macros */ \
	/*	distance between fft data elements is 16 */ \
	/*	do 8 macros each processing 8 data values */ \
	\
	for(unsigned int loop=8; loop != 0; loop--){ \
		x1cl_half_eight_reals_fft_2(rsi, 64); \
	} \
	\
	/* Do 16 four_complex_fft macros */ \
	/*	distance between fft data elements is 16 */ \
	/*	do 8 macros each processing 16 data values */ \
	\
	rdi = (uintptr_t)g->sincos3;		/* Load sin/cos pointer */ \
	for(unsigned int loop=8; loop != 0; loop--){ \
		x2cl_four_complex_fft(rsi, 64, 8*64); \
	} \
	rsi += -8*64+16*64;	/* Next source pointer */ \
	\
	/* Do 32 four_complex_fft macros */ \
	/*	distance between fft data elements is 16 */ \
	/*	do 8 macros each processing 32 data values */ \
	\
	rdi += XMM_SCD;		/* Next sine/cosine pointer */ \
	for(unsigned int loop=8; loop != 0; loop--){ \
		x4cl_four_complex_fft(rsi, 64, 8*64, 16*64, rdi); \
	} \
	\
	/* Join common code to perform rest of the FFT */ \
	/* Do FFT levels 6,7 */ \
	/*	do 8/10/12/14 macros each processing 32 data values */ \
	/*	mov	ax, WORD PTR count4	;; 4/5/6/7 iterations of 2 */ \
	\
	xlast_two_fft_levels; \
} \
/*rcx = 8*64;*/ \
cmn512(8*64); \
\
/* Finish inverse FFT levels 4,5 */ \
\
/* Do 8 eight_reals_unfft_2 macros */ \
/*	distance between fft data elements is 16 */ \
/*	do 4 macros each processing 16 data values */ \
\
for(unsigned int loop=4; loop != 0; loop--){ \
	x2cl_half_eight_reals_unfft_2(rsi, 64, 4*64); \
} \
rsi += -4*64;		/* Restore source pointer */ \
\
/* Do inverse FFT levels 1,2,3 */ \
/* On input the 64-byte cache lines hold these data values: */ \
/*	0	16	32	48 */ \
/*	2	... */ \
/*	... */ \
/*	14 */ \
/*	64	96	128	160 */ \
/*	... */ \
/* On output the 64-byte cache lines hold these data values: */ \
/*	0	16	224	240 */ \
/*	2	... */ \
/*	... */ \
/*	14 */ \
/*	32	64	256	288 */ \
/*	34 */ \
/*	... */ \
/*	62 */ \
/*	96	160	320	384 */ \
/*	98 */ \
/*	... */ \
/*	158 */ \
\
/* Do 64 seven_reals_last_unfft macros */ \
/*	distance between fft data elements is 64 */ \
/*	do 8 macros each processing 56 data values */ \
\
xseven_reals_unfft_preload; \
for(unsigned int loop=8; loop != 0; loop--){ \
	x7cl_seven_reals_last_unfft(rsi, 64, 8*64); \
} \
\
xfft_3_ret; \
\
\
\
/* Perform a 512-element FFT. */ \
\

#define xfft512(type) \
\
/* Do a multiply with pre-FFTed inputs */ \
\
if(g->ffttype == 4) { \
	xmiddle_4(g); \
}else{ \
	\
	/* Do FFT levels 1,2,3 */ \
	/* Values 0-511 is real data. */ \
	/* */ \
	/* On input the 64-byte cache lines hold these data values: */ \
	/*	0	128	1	129	256	384	257	385 */ \
	/*	2	... */ \
	/*	... */ \
	/*	126 */ \
	/* On output the 64-byte cache lines hold these data values: */ \
	/*	0	1	64	65	128	129	192	193 */ \
	/*	2	... */ \
	/*	... */ \
	/*	62	... */ \
	/*	256	... */ \
	/*	... */ \
	\
	/* Do 64 eight_reals_first_fft macros */ \
	/*	distance between fft data elements is 64 */ \
	/*	do 32 macros each processing 16 data values */ \
	\
	for(unsigned int loop=32; loop != 0; loop--){ \
		s2cl_eight_reals_first_fft(rsi, 64, 32*64); \
	} \
	rsi += -32*64;		/* Restore source pointer */ \
	\
	/* Do FFT levels 4,5 */ \
	/* Values 0-63 is real data, 64-127 is semi-real data, 128-511 is */ \
	/* complex data. */ \
	/* */ \
	/* On input the 64-byte cache lines hold these data values: */ \
	/*	0	1	64	65	128	129	192	193 */ \
	/*	2	... */ \
	/*	... */ \
	/*	62	... */ \
	/*	256	... */ \
	/*	... */ \
	/* On output the 64-byte cache lines hold these data values: */ \
	/*	0	1	16	17	32	33	48	49 */ \
	/*	2	... */ \
	/*	... */ \
	/*	14	... */ \
	/*	64	... */ \
	/*	... */ \
	\
	/* Do 8 eight_reals_fft_2 macros */ \
	/* Do 8 nop_two_two_complex_fft_2 macros */ \
	/* Do 16 four_complex_fft macros */ \
	/*	distance between fft data elements is 16 */ \
	/*	do 8 macros each processing 32 data values */ \
	\
	rdi = (uintptr_t)g->sincos3;		/* Load sin/cos pointer */ \
	for(unsigned int loop=8; loop != 0; loop--){ \
		x4cl_eight_reals_fft_2(rsi, 64, 8*64, 16*64, rdi); \
	} \
	rsi += -8*64+32*64;	/* Next source pointer */ \
	\
	/* Do 32 four_complex_fft macros */ \
	/*	distance between fft data elements is 16 */ \
	/*	do 8 macros each processing 32 data values */ \
	\
	rdi += XMM_SCD;		/* Next sine/cosine pointer */ \
	for(unsigned int loop=8; loop != 0; loop--){ \
		x4cl_four_complex_fft(rsi, 64, 8*64, 16*64, rdi); \
	} \
	\
	/* Do FFT levels 6,7 */ \
	/*	do 8/10/12/14 macros each processing 32 data values */ \
	/*	mov	ax, WORD PTR count4	;; 4/5/6/7 iterations of 2 */ \
	/* Do FFT levels 8,9 */ \
	/*	do 19/23/27/31 macros each processing 16 data values */ \
	/*	mov	eax, count1		;; 3/7 then 16/24 iterations */ \
	\
	xlast_two_fft_levels; \
} \
/*rcx = 16*64;*/ \
cmn512(16*64); \
\
/* Finish the inverse FFT levels 4,5 */ \
\
/* Do 8 eight_reals_unfft_2 macros */ \
/* Do 8 nop_two_two_complex_unfft_2 macros */ \
/*	distance between fft data elements is 16 */ \
/*	do 4 macros each processing 32 data values */ \
\
for(unsigned int loop=4; loop != 0; loop--){ \
	x4cl_eight_reals_unfft_2(rsi, 64, 4*64, 8*64); \
} \
rsi += -4*64;		/* Restore source pointer */ \
\
/* Do inverse FFT levels 1,2,3 */ \
/* On input the 64-byte cache lines hold these data values: */ \
/*	0	1	32	33	64	65	96	97 */ \
/*	2	... */ \
/*	...	... */ \
/*	30	... */ \
/*	128	... */ \
/*	... */ \
/* On output the 64-byte cache lines hold these data values: */ \
/*	0	1	128	129	256	257	384	385 */ \
/*	2	... */ \
/*	... */ \
/*	126 */ \
\
/* Do 64 eight_reals_last_unfft macros */ \
/*	distance between fft data elements is 64 */ \
/*	do 16 macros each processing 32 data values */ \
\
for(unsigned int loop=16; loop != 0; loop--){ \
	x4cl_eight_reals_last_unfft(rsi, 64, 16*64, 32*64); \
} \
\
xfft_3_ret; \
\
\
/* Common code for length 320, 384, 448 and 512 FFTs */ \
\
/* Do inverse FFT levels 6,7 */ \
/* On input the 64-byte cache lines hold these data values: */ \
/*	0	1	2	3	4	5	6	7 */ \
/*	8	9	10	11	... */ \
/*	... */ \
/* On output the 64-byte cache lines hold these data values: */ \
/*	0	1	8	9	16	17	24	25 */ \
/*	2	... */ \
/*	4	... */ \
/*	6	... */ \
/*	32	... */ \
/*	... */ \
\
/* Do 2 eight_reals_unfft_2 macros */ \
/* Do 2 nop_two_two_complex_unfft_2 macros */ \
/*	distance between fft data elements is 4 */ \
/*	do 1 macro processing 32 data values */ \

#define cmn512(RCX1) \
/*push_amt = SZPTR;*/ \
x4cl_eight_reals_unfft_2(rsi, 4*64, 64, 2*64); \
\
/* Do 36/44/52/60 four_complex_unfft macros */ \
/*	distance between fft data elements is 4 */ \
/*	do 9/11/13/15 macros each processing 32 data values */ \
\
rdi = (uintptr_t)g->sincos2;/* Load sin/cos pointer */ \
for(unsigned int loop = g->count2 & 0xFF; loop != 0; loop--){ \
	x4cl_four_complex_unfft(rsi, 4*64, 64, 2*64, rdi); \
	rdi += XMM_SCD;		/* Next sine/cosine pointer */ \
} \
rsi = (uintptr_t)g->DESTARG;/* Restore source pointer */ \
\
/* Do inverse FFT levels 4,5 */ \
/* On input the 64-byte cache lines hold these data values: */ \
/*	0	1	8	9	16	17	24	25 */ \
/*	2	... */ \
/*	4	... */ \
/*	6	... */ \
/*	32	... */ \
/*	... */ \
/* On output the 64-byte cache lines hold these data values: */ \
/*	0	1	32	33	64	65	96	97 */ \
/*	2	... */ \
/*	...	... */ \
/*	30	... */ \
/*	128	... */ \
/*	... */ \
\
/* Do 32/48 four_complex_unfft macros */ \
/*	distance between fft data elements is 16 */ \
/*	do 8/12 macros each processing 32 data values */ \
\
rdi = (uintptr_t)g->sincos3;		/* Load sin/cos pointer */ \
rsi += RCX1;		/* Load source pointer */ \
for(unsigned int loop = g->count3 & 0xFF; loop != 0; loop--){ \
	for(unsigned int loopB=4; loopB != 0; loopB--){ \
		x4cl_four_complex_unfft(rsi, 64, 4*64, 8*64, rdi); \
	} \
	rsi += -4*64+16*64;	/* Next source pointer */ \
	rdi += XMM_SCD;		/* Next sine/cosine pointer */ \
} \
rsi = (uintptr_t)g->DESTARG;/* Restore source pointer */ \
\
/* This ends the common inverse FFT code */ \
\
/* ret; push_amt = 0; */ \
\
\
\
/* Perform a 512-element negacyclic FFT */ \
\

#define xfft512p(type) \
\
/* Do a multiply with pre-FFTed inputs */ \
\
if(g->ffttype == 4) { \
	xmiddle_4p(g); \
}else{ \
	\
	/* Do FFT levels 1,2,3 */ \
	/* */ \
	/* On input the 64-byte cache lines hold these data values: */ \
	/*	0	128	1	129	256	384	257	385 */ \
	/*	2	... */ \
	/*	... */ \
	/*	126 */ \
	/* On output the 64-byte cache lines hold these data values: */ \
	/*	0	1	64	65	128	129	192	193 */ \
	/*	2	... */ \
	/*	... */ \
	/*	62	... */ \
	/*	256	... */ \
	/*	... */ \
	\
	/* Do 64 four_complex_first_fft macros */ \
	/*	distance between fft data elements is 64 */ \
	/*	do 32 macros each processing 16 data values */ \
	\
	rdi = (uintptr_t)g->plus1_premults;/* Address of premultiplier table */ \
	for(unsigned int loop=32; loop != 0; loop--){ \
		s2cl_four_complex_first_fft(rsi, 64, 32*64); \
		rdi += 128;		/* Next premultiplier pointer */ \
	} \
	rsi += -32*64;		/* Restore source pointer */ \
	\
	/* Do FFT levels 4,5 */ \
	/* */ \
	/* On input the 64-byte cache lines hold these data values: */ \
	/*	0	1	64	65	128	129	192	193 */ \
	/*	2	... */ \
	/*	... */ \
	/*	62	... */ \
	/*	256	... */ \
	/*	... */ \
	/* On output the 64-byte cache lines hold these data values: */ \
	/*	0	1	16	17	32	33	48	49 */ \
	/*	2	... */ \
	/*	... */ \
	/*	14	... */ \
	/*	64	... */ \
	/*	... */ \
	\
	/* Do 64 four_complex_fft macros */ \
	/*	distance between fft data elements is 16 */ \
	/*	do 16 macros each processing 32 data values */ \
	\
	/* Do FFT levels 6,7 */ \
	/*	do 12/16 macros each processing 32 data values */ \
	/*	mov	ax, WORD PTR count4	;; 6/8 iterations of 2 */ \
	/* Do FFT levels 8,9 */ \
	/*	do 24/32 macros each processing 16 data values */ \
	/*	mov	eax, count1		;; 24/32 iterations */ \
	\
	xlast_four_fft_levels_p; \
} \
cmn512p; \
\
/* Do inverse FFT levels 1,2,3 */ \
/* On input the 64-byte cache lines hold these data values: */ \
/*	0	1	32	33	64	65	96	97 */ \
/*	2	... */ \
/*	...	... */ \
/*	30	... */ \
/*	128	... */ \
/*	... */ \
/* On output the 64-byte cache lines hold these data values: */ \
/*	0	1	128	129	256	257	384	385 */ \
/*	2	... */ \
/*	... */ \
/*	126 */ \
\
/* Do 64 four_complex_last_unfft macros */ \
/*	distance between fft data elements is 64 */ \
/*	do 16 macros each processing 32 data values */ \
\
rdi = (uintptr_t)g->plus1_premults;/* Address of premultiplier table */ \
for(unsigned int loop=16; loop != 0; loop--){ \
	x4cl_four_complex_last_unfft(rsi, 64, 16*64, 32*64, 16*128); \
	rdi += 128;		/* Next premultiplier pointer */ \
} \
\
xfft_3_ret;

/* Common inverse FFT code for negacyclic lengths 384 and 512 */ \

/* Do inverse FFT levels 6,7 */ \
/* On input the 64-byte cache lines hold these data values: */ \
/*	0	1	2	3	4	5	6	7 */ \
/*	8	9	10	11	... */ \
/*	... */ \
/* On output the 64-byte cache lines hold these data values: */ \
/*	0	1	8	9	16	17	24	25 */ \
/*	2	... */ \
/*	4	... */ \
/*	6	... */ \
/*	32	... */ \
/*	... */ \

/* Do 48/64 four_complex_unfft macros */ \
/*	distance between fft data elements is 4 */ \
/*	do 12/16 macros each processing 32 data values */ \

#define cmn512p \
/* push_amt = SZPTR;*/ \
rdi = (uintptr_t)g->sincos2;/* Load sin/cos pointer */ \
for(unsigned int loop = g->count2 & 0xFF; loop != 0; loop--){ \
	x4cl_four_complex_unfft(rsi, 4*64, 64, 2*64, rdi); \
	rdi += XMM_SCD;		/* Next sine/cosine pointer */ \
} \
rsi = (uintptr_t)g->DESTARG;/* Restore source pointer */ \
\
/* Do inverse FFT levels 4,5 */ \
/* On input the 64-byte cache lines hold these data values: */ \
/*	0	1	8	9	16	17	24	25 */ \
/*	2	... */ \
/*	4	... */ \
/*	6	... */ \
/*	32	... */ \
/*	... */ \
/* On output the 64-byte cache lines hold these data values: */ \
/*	0	1	32	33	64	65	96	97 */ \
/*	2	... */ \
/*	...	... */ \
/*	30	... */ \
/*	128	... */ \
/*	... */ \
\
/* Do 48/64 four_complex_unfft macros */ \
/*	distance between fft data elements is 16 */ \
/*	do 12/16 macros each processing 32 data values */ \
\
rdi = (uintptr_t)g->sincos3;		/* Load sin/cos pointer */ \
for(unsigned int loop = g->count3 & 0xFF; loop != 0; loop--){ \
	for(unsigned int loopB=4; loopB != 0; loopB--){ \
		x4cl_four_complex_unfft(rsi, 64, 4*64, 8*64, rdi); \
	} \
	rsi += -4*64+16*64;	/* Next source pointer */ \
	rdi += XMM_SCD;		/* Next sine/cosine pointer */ \
} \
rsi = (uintptr_t)g->DESTARG;/* Restore source pointer */ \
\
/* This ends the common inverse FFT code */ \
\
/* ret; push_amt = 0; */



/* Perform a 640-element FFT. */


#define xfft640(type) \
\
/* Do a multiply with pre-FFTed inputs */ \
\
if(g->ffttype == 4) { \
	xmiddle_4(g); \
}else{ \
	\
	/* Do FFT levels 1,2,3 */ \
	/* Values 0-639 is real data. */ \
	/* */ \
	/* On input the 64-byte cache lines hold these data values: */ \
	/*	0	32	1	33	320	352	321	353 */ \
	/*	2	... */ \
	/*	... */ \
	/*	30 */ \
	/*	64	192	65	193	384	512	385	513 */ \
	/*	66 */ \
	/*	... */ \
	/*	190 */ \
	/* On output the 64-byte cache lines hold these data values (in +1 format): */ \
	/*	0	32	64	96 */ \
	/*	2	... */ \
	/*	... */ \
	/*	30 */ \
	/*	128	256	384	512 */ \
	/*	130 */ \
	/*	... */ \
	/*	254 */ \
	\
	/* Do 128 five_reals_first_fft macros */ \
	/*	distance between fft data elements is 128 */ \
	/*	do 16 macros each processing 40 data values */ \
	\
	xfive_reals_fft_preload;		/* Preload const1s for five real macros */ \
	for(unsigned int loop=16; loop != 0; loop--){ \
		s5cl_five_reals_first_fft(rsi, 64, 16*64); \
	} \
	rsi += -16*64;		/* Restore source pointer */ \
	\
	/* Do FFT levels 4,5 */ \
	/* Values 0-127 is real data, 128-639 is complex data. */ \
	\
	/* Do 16 eight_reals_fft_2 macros */ \
	/*	distance between fft data elements is 32 */ \
	/*	do 16 macros each processing 8 data values */ \
	\
	for(unsigned int loop=16; loop != 0; loop--){ \
		x1cl_half_eight_reals_fft_2(rsi, 64); \
	} \
	\
	/* Do 32 four_complex_fft macros */ \
	/*	distance between fft data elements is 32 */ \
	/*	do 16 macros each processing 32 data values */ \
	\
	rdi = (uintptr_t)g->sincos4;		/* Load sin/cos pointer */ \
	for(unsigned int loop=16; loop != 0; loop--){ \
		x4cl_four_complex_fft(rsi, 64, 16*64, 32*64, rdi); \
	} \
	\
	/* Join common code to perform rest of the FFT */ \
	\
	cmn1024_2; \
} \
/*rcx = 16*64;*/ \
cmn1024(16*64); \
\
/* Finish the inverse FFT levels 4,5 */ \
\
/* Do 16 eight_reals_unfft_2 macros */ \
/*	distance between fft data elements is 32 */ \
/*	do 8 macros each processing 16 data values */ \
\
for(unsigned int loop=8; loop != 0; loop--){ \
	x2cl_half_eight_reals_unfft_2(rsi, 64, 8*64); \
} \
rsi += -8*64;		/* Restore source pointer */ \
\
/* Do inverse FFT levels 1,2,3 */ \
/* On output the 64-byte cache lines hold these data values: */ \
/*	0	32	320	352 */ \
/*	2	... */ \
/*	... */ \
/*	30 */ \
/*	64	192	384	512 */ \
/*	66 */ \
/*	... */ \
/*	190 */ \
\
/* Do 128 five_reals_last_unfft macros */ \
/*	distance between fft data elements is 128 */ \
/*	do 16 macros each processing 40 data values */ \
\
xfive_reals_unfft_preload;	/* Preload const1s for five real macros */ \
for(unsigned int loop=16; loop != 0; loop--){ \
	x5cl_five_reals_last_unfft(rsi, 64, 16*64); \
} \
\
xfft_3_ret; \
\
\
\
/* Perform a 768-element FFT. */ \
\

#define xfft768(type) \
\
/* Do a multiply with pre-FFTed inputs */ \
\
if(g->ffttype == 4) { \
	xmiddle_4(g); \
}else{ \
	\
	/* Do FFT levels 1,2,3 */ \
	/* Values 0-767 is real data. */ \
	/* */ \
	/* On input the 64-byte cache lines hold these data values: */ \
	/*	0	192	1	193	384	576	385	577 */ \
	/*	2 */ \
	/*	... */ \
	/*	190 */ \
	/* On output the 64-byte cache lines hold these data values (in +1 format): */ \
	/*	0	64	128	192 */ \
	/*	2	.. */ \
	/*	... */ \
	/*	62 */ \
	/*	256	384	512	640 */ \
	/*	258 */ \
	/*	... */ \
	/*	382 */ \
	\
	/* Do 128 six_reals_first_fft macros */ \
	/*	distance between fft data elements is 128 */ \
	/*	do 32 macros each processing 24 data values */ \
	\
	xsix_reals_fft_preload; \
	for(unsigned int loop=32; loop != 0; loop--){ \
		s3cl_six_reals_first_fft(rsi, 64, 32*64); \
	} \
	rsi += -32*64;		/* Restore source pointer */ \
	\
	/* Do FFT levels 4,5 */ \
	/* Values 0-127 is real data, 128-255 is semi-real data, 256-767 is */ \
	/* complex data. */ \
	/* */ \
	/* On input the 64-byte cache lines hold these data values: */ \
	/*	0	64	128	192 */ \
	/*	2	.. */ \
	/*	... */ \
	/*	62 */ \
	/*	256	384	512	640 */ \
	/*	258 */ \
	/*	... */ \
	/*	382 */ \
	/* On output the 64-byte cache lines hold these data values: */ \
	/*	0	1	32	33	64	65	96	97 */ \
	/*	2	... */ \
	/*	... */ \
	/*	30	... */ \
	/*	128	... */ \
	/*	... */ \
	\
	/* Do 16 eight_reals_fft_2 macros */ \
	/* Do 16 nop_two_two_complex_fft_2 macros */ \
	/* Do 32 four_complex_fft macros */ \
	/*	distance between fft data elements is 32 */ \
	/*	do 16 macros each processing 16 data values */ \
	\
	for(unsigned int loop=16; loop != 0; loop--){ \
		x2cl_eight_reals_fft_2(rsi, 64, 16*64); \
	} \
	rsi += -16*64+32*64;	/* Next source pointer */ \
	\
	/* Do 32 four_complex_fft macros */ \
	/*	distance between fft data elements is 32 */ \
	/*	do 16 macros each processing 32 data values */ \
	\
	rdi = (uintptr_t)g->sincos4;		/* Load sin/cos pointer */ \
	for(unsigned int loop=16; loop != 0; loop--){ \
		x4cl_four_complex_fft(rsi, 64, 16*64, 32*64, rdi); \
	} \
	\
	/* Join common code to perform rest of the FFT */ \
	\
	cmn1024_2; \
} \
/*rcx = 32*64;*/ \
cmn1024(32*64); \
\
/* Finish inverse FFT levels 4,5 */ \
\
/* Do 16 eight_reals_unfft_2 macros */ \
/* Do 16 nop_two_two_complex_unfft_2 macros */ \
/*	distance between fft data elements is 32 */ \
/*	do 8 macros each processing 32 data values */ \
\
for(unsigned int loop=8; loop != 0; loop--){ \
	x4cl_eight_reals_unfft_2(rsi, 64, 8*64, 16*64); \
} \
rsi += -8*64;		/* Restore source pointer */ \
\
/* Do inverse FFT levels 1,2,3 */ \
/* On input the 64-byte cache lines hold these data values: */ \
/*	0	1	64	65	128	129	192	193 */ \
/*	2	... */ \
/*	...	... */ \
/*	62	... */ \
/*	256	... */ \
/*	... */ \
/* On output the 64-byte cache lines hold these data values: */ \
/*	0	192	384	576 */ \
/*	2 */ \
/*	... */ \
/*	190 */ \
\
/* Do 128 six_reals_last_unfft macros */ \
/*	distance between fft data elements is 128 */ \
/*	do 32 macros each processing 24 data values */ \
\
xsix_reals_unfft_preload; \
for(unsigned int loop=32; loop != 0; loop--){ \
	x3cl_six_reals_last_unfft(rsi, 64, 32*64); \
} \
\
xfft_3_ret; \
\
\
\
/* Perform a 768-element negacyclic FFT */ \
\

#define xfft768p(type) \
\
/* Do a multiply with pre-FFTed inputs */ \
\
if(g->ffttype == 4) { \
	xmiddle_4p(g); \
}else{ \
	\
	/* Do FFT levels 1,2,3 */ \
	/* */ \
	/* On input the 64-byte cache lines hold these data values: */ \
	/*	0	192	1	193	384	576	385	577 */ \
	/*	2 */ \
	/*	... */ \
	/*	190 */ \
	/* On output the 64-byte cache lines hold these data values (in +1 format): */ \
	/*	0	64	128	192 */ \
	/*	2	.. */ \
	/*	... */ \
	/*	62 */ \
	/*	256	384	512	640 */ \
	/*	258 */ \
	/*	... */ \
	/*	382 */ \
	\
	/* Do 128 three_complex_first_fft macros */ \
	/*	distance between fft data elements is 128 */ \
	/*	do 32 macros each processing 24 data values */ \
	\
	rdi = (uintptr_t)g->plus1_premults;/* Address of premultiplier table */ \
	for(unsigned int loop=32; loop != 0; loop--){ \
		s3cl_three_complex_first_fft(rsi, 64, 32*64); \
		rdi += 192;		/* Next premultiplier pointer */ \
	} \
	rsi += -32*64;		/* Restore source pointer */ \
	\
	/* Do FFT levels 4 */ \
	/* */ \
	/* On input the 64-byte cache lines hold these data values: */ \
	/*	0	64	128	192 */ \
	/*	2	.. */ \
	/*	... */ \
	/*	62 */ \
	/*	256	384	512	640 */ \
	/*	258 */ \
	/*	... */ \
	/*	382 */ \
	/* On output the 64-byte cache lines hold these data values: */ \
	/*	0	64	128	192 */ \
	/*	2	... */ \
	/*	... */ \
	/*	62 */ \
	/*	256	... */ \
	/*	... */ \
	\
	/* Do 96 two_two_complex_fft macros */ \
	/*	distance between fft data elements is 64 */ \
	/*	do 48 macros each processing 16 data values */ \
	\
	rdi = (uintptr_t)g->sincos4;		/* Load sin/cos pointer */ \
	for(unsigned int loop=16; loop != 0; loop--){ \
		x2cl_two_complex_fft_in_place(rsi, 64, 16*64, rdi); \
	} \
	rsi += -16*64+32*64;	/* Next source pointer */ \
	rdi += XMM_SCD;		/* Next sin/cos pointer */ \
	for(unsigned int loop=32; loop != 0; loop--){ \
		x2cl_two_complex_fft(rsi, 64, 32*64, rdi); \
	} \
	rsi += -32*64-32*64;	/* Restore source pointer */ \
	\
	/* Join common code to perform rest of the FFT */ \
	\
	xlast_four_fft_levels_p; \
} \
cmn1024p; \
\
/* Do inverse FFT levels 1,2,3 */ \
/* On input the 64-byte cache lines hold these data values: */ \
/*	0	1	64	65	128	129	192	193 */ \
/*	2	... */ \
/*	...	... */ \
/*	62	... */ \
/*	256	... */ \
/*	... */ \
/* On output the 64-byte cache lines hold these data values: */ \
/*	0	192	384	576 */ \
/*	2 */ \
/*	... */ \
/*	190 */ \
\
/* Do 128 three_complex_last_unfft macros */ \
/*	distance between fft data elements is 128 */ \
/*	do 32 macros each processing 24 data values */ \
\
rdi = (uintptr_t)g->plus1_premults;/* Address of premultiplier table */ \
for(unsigned int loop=32; loop != 0; loop--){ \
	x3cl_three_complex_last_unfft(rsi, 64, 32*64); \
	rdi += 192;		/* Next premultiplier pointer */ \
} \
\
xfft_3_ret; \
\
\
\
/* Perform a 896-element FFT. */ \
\

#define xfft896(type) \
\
/* Do a multiply with pre-FFTed inputs */ \
\
if(g->ffttype == 4) { \
	xmiddle_4(g); \
}else{ \
	\
	/* Do FFT levels 1,2,3 */ \
	/* Values 0-895 is real data. */ \
	/* */ \
	/* On input the 64-byte cache lines hold these data values: */ \
	/*	0	32	1	33	448	480	449	481 */ \
	/*	2	... */ \
	/*	... */ \
	/*	30 */ \
	/*	64	128	65	129	512	576	513	577 */ \
	/*	66 */ \
	/*	... */ \
	/*	126 */ \
	/*	192	320	193	321	640	768	641	769 */ \
	/*	194 */ \
	/*	... */ \
	/*	318 */ \
	/* On output the 64-byte cache lines hold these data values (in +1 format): */ \
	/*	0	32	64	96 */ \
	/*	2	... */ \
	/*	... */ \
	/*	30 */ \
	/*	128	192	256	320 */ \
	/*	130 */ \
	/*	... */ \
	/*	190 */ \
	/*	384	512	640	768 */ \
	/*	386 */ \
	/*	... */ \
	/*	510 */ \
	\
	/* Do 128 seven_reals_first_fft macros */ \
	/*	distance between fft data elements is 128 */ \
	/*	do 16 macros each processing 56 data values */ \
	\
	xseven_reals_fft_preload; \
	for(unsigned int loop=16; loop != 0; loop--){ \
		s7cl_seven_reals_first_fft(rsi, 64, 16*64); \
	} \
	rsi += -16*64;		/* Restore source pointer */ \
	\
	/* Do FFT levels 4,5 */ \
	/* Values 0-127 is real data, 128-895 is complex data. */ \
	\
	/* Do 16 eight_reals_fft_2 macros */ \
	/*	distance between fft data elements is 32 */ \
	/*	do 16 macros each processing 8 data values */ \
	\
	for(unsigned int loop=16; loop != 0; loop--){ \
		x1cl_half_eight_reals_fft_2(rsi, 64); \
	} \
	\
	/* Do 32 four_complex_fft macros */ \
	/*	distance between fft data elements is 32 */ \
	/*	do 16 macros each processing 16 data values */ \
	\
	rdi = (uintptr_t)g->sincos4;		/* Load sin/cos pointer */ \
	for(unsigned int loop=16; loop != 0; loop--){ \
		x2cl_four_complex_fft(rsi, 64, 16*64); \
	} \
	rsi += -16*64+32*64;	/* Next source pointer */ \
	\
	/* Do 64 four_complex_fft macros */ \
	/*	distance between fft data elements is 32 */ \
	/*	do 16 macros each processing 32 data values */ \
	\
	rdi += XMM_SCD;		/* Next sine/cosine pointer */ \
	for(unsigned int loop=16; loop != 0; loop--){ \
		x4cl_four_complex_fft(rsi, 64, 16*64, 32*64, rdi); \
	} \
	\
	/* Join common code to perform rest of the FFT */ \
	\
	cmn1024_2; \
} \
/*rcx = 16*64;*/ \
cmn1024(16*64); \
\
/* Finish inverse FFT levels 4,5 */ \
\
/* Do 16 eight_reals_unfft_2 macros */ \
/*	distance between fft data elements is 32 */ \
/*	do 8 macros each processing 16 data values */ \
\
for(unsigned int loop=8; loop != 0; loop--){ \
	x2cl_half_eight_reals_unfft_2(rsi, 64, 8*64); \
} \
rsi += -8*64;		/* Restore source pointer */ \
\
/* Do inverse FFT levels 1,2,3 */ \
/* On input the 64-byte cache lines hold these data values: */ \
/*	0	32	64	96 */ \
/*	2	... */ \
/*	... */ \
/*	30 */ \
/*	128	192	256	320 */ \
/*	... */ \
/* On output the 64-byte cache lines hold these data values: */ \
/*	0	32	448	480 */ \
/*	2	... */ \
/*	... */ \
/*	30 */ \
/*	64	128	512	576 */ \
/*	66 */ \
/*	... */ \
/*	126 */ \
/*	192	320	640	768 */ \
/*	194 */ \
/*	... */ \
/*	318 */ \
\
/* Do 128 seven_reals_last_unfft macros */ \
/*	distance between fft data elements is 128 */ \
/*	do 16 macros each processing 56 data values */ \
\
xseven_reals_unfft_preload; \
for(unsigned int loop=16; loop != 0; loop--){ \
	x7cl_seven_reals_last_unfft(rsi, 64, 16*64); \
} \
\
xfft_3_ret; \
\
\
\
/* Perform a 1024-element FFT. */ \
\

#define xfft1024(type) \
\
/* Do a multiply with pre-FFTed inputs */ \
\
if(g->ffttype == 4) { \
	xmiddle_4(g); \
}else{ \
	\
	/* Do FFT levels 1,2,3 */ \
	/* Values 0-1023 is real data. */ \
	/* */ \
	/* On input the 64-byte cache lines hold these data values: */ \
	/*	0	256	1	257	512	768	513	769 */ \
	/*	2	... */ \
	/*	... */ \
	/*	254 */ \
	/* On output the 64-byte cache lines hold these data values: */ \
	/*	0	1	128	129	256	257	384	385 */ \
	/*	2	... */ \
	/*	... */ \
	/*	126	... */ \
	/*	512	... */ \
	/*	... */ \
	\
	/* Do 128 eight_reals_first_fft macros */ \
	/*	distance between fft data elements is 128 */ \
	/*	do 64 macros each processing 16 data values */ \
	\
	for(unsigned int loop=64; loop != 0; loop--){ \
		s2cl_eight_reals_first_fft(rsi, 64, 64*64); \
	} \
	rsi += -64*64;		/* Restore source pointer */ \
	\
	/* Do FFT levels 4,5 */ \
	/* Values 0-127 is real data, 128-255 is semi-real data, 256-1023 is */ \
	/* complex data. */ \
	/* */ \
	/* On input the 64-byte cache lines hold these data values: */ \
	/*	0	1	128	129	256	257	384	385 */ \
	/*	2	... */ \
	/*	... */ \
	/*	126	... */ \
	/*	512	... */ \
	/*	... */ \
	/* On output the 64-byte cache lines hold these data values: */ \
	/*	0	1	32	33	64	65	96	97 */ \
	/*	2	... */ \
	/*	... */ \
	/*	30	... */ \
	/*	128	... */ \
	/*	... */ \
	\
	/* Do 16 eight_reals_fft_2 macros */ \
	/* Do 16 nop_two_two_complex_fft_2 macros */ \
	/* Do 32 four_complex_fft macros */ \
	/*	distance between fft data elements is 32 */ \
	/*	do 16 macros each processing 32 data values */ \
	\
	rdi = (uintptr_t)g->sincos4;		/* Load sin/cos pointer */ \
	for(unsigned int loop=16; loop != 0; loop--){ \
		x4cl_eight_reals_fft_2(rsi, 64, 16*64, 32*64, rdi); \
	} \
	rsi += -16*64+64*64;	/* Next source pointer */ \
	\
	/* Do 64 four_complex_fft macros */ \
	/*	distance between fft data elements is 32 */ \
	/*	do 16 macros each processing 32 data values */ \
	\
	rdi += XMM_SCD;		/* Next sine/cosine pointer */ \
	for(unsigned int loop=16; loop != 0; loop--){ \
		x4cl_four_complex_fft(rsi, 64, 16*64, 32*64, rdi); \
	} \
	\
	/* Call common code for rest of FFT */ \
	/* Do FFT levels 6,7,8 */ \
	/* Do FFT levels 9,10 */ \
	/*	do 39/47/55/63 macros each processing 16 data values */ \
	/*	mov	eax, count1		;; 7/15 then 32/48 iterations */ \
	\
	cmn1024_2; \
} \
/*rcx = 32*64;*/ \
cmn1024(32*64); \
\
/* Finish the inverse FFT levels 4,5 */ \
\
/* Do 16 eight_reals_unfft_2 macros */ \
/* Do 16 nop_two_two_complex_unfft_2 macros */ \
/*	distance between fft data elements is 32 */ \
/*	do 8 macros each processing 32 data values */ \
\
for(unsigned int loop=8; loop != 0; loop--){ \
	x4cl_eight_reals_unfft_2(rsi, 64, 8*64, 16*64); \
} \
rsi += -8*64;		/* Restore source pointer */ \
\
/* Do inverse FFT levels 1,2,3 */ \
/* On input the 64-byte cache lines hold these data values: */ \
/*	0	1	64	65	128	129	192	193 */ \
/*	2	... */ \
/*	...	... */ \
/*	62	... */ \
/*	256	... */ \
/*	... */ \
/* On output the 64-byte cache lines hold these data values: */ \
/*	0	1	256	257	512	513	768	769 */ \
/*	2	... */ \
/*	... */ \
/*	254 */ \
\
/* Do 128 eight_reals_last_unfft macros */ \
/*	distance between fft data elements is 128 */ \
/*	do 32 macros each processing 32 data values */ \
\
for(unsigned int loop=32; loop != 0; loop--){ \
	x4cl_eight_reals_last_unfft(rsi, 64, 32*64, 64*64); \
} \
\
xfft_3_ret; \

/* Common code for length 640, 768, 896, and 1024 FFTs */ \

#define cmn1024_2 \
/*push_amt = SZPTR;*/ \
\
/* Do FFT level 6 */ \
/* Values 0-31 is real data, 32-63 is semi-real data, 64-1023 is */ \
/* complex data. */ \
/* */ \
/* On input the 64-byte cache lines hold these data values: */ \
/*	0	1	32	33	64	65	96	97 */ \
/*	2	... */ \
/*	... */ \
/*	30	... */ \
/*	128	... */ \
/*	... */ \
/* On output the 64-byte cache lines hold these data values: */ \
/*	0	1	16	17	32	33	48	49 */ \
/*	2	... */ \
/*	... */ \
/*	14	... */ \
/*	64	... */ \
/*	... */ \
\
/* Do 8 eight_reals_fft_2 macros */ \
/* Do 8 two_two_complex_fft macros */ \
/*	distance between fft data elements is 16 */ \
/*	do 8 macros each processing 16 data values */ \
\
rsi = (uintptr_t)g->DESTARG;/* Next source pointer */ \
rdi = (uintptr_t)g->sincos3;		/* Load sin/cos pointer */ \
for(unsigned int loop=8; loop != 0; loop--){ \
	x2cl_eight_reals_fft_1(rsi, 64, 8*64); \
} \
\
/* Do 64/80/96/112 two_two_complex_fft macros */ \
/*	distance between fft data elements is 16 */ \
/*	do 32/40/48/56 macros each processing 16 data values */ \
\
rsi += -8*64+16*64;	/* Load source pointer */ \
rdi += XMM_SCD;		/* Next sine/cosine pointer */ \
for(unsigned int loop = (g->count5 >> 8) & 0xFF; loop != 0; loop--){ \
	for(unsigned int loop=8; loop != 0; loop--){ \
		x2cl_two_complex_fft(rsi, 64, 8*64, rdi); \
	} \
	rsi += -8*64+16*64;	/* Next source pointer */ \
	rdi += 2*XMM_SCD;		/* Next sine/cosine pointer */ \
} \
\
/* Do FFT levels 7,8 */ \
/*	do 18/22/26/30 macros each processing 32 data values */ \
/*	mov	ax, WORD PTR count4	;; 1/3 then 8/12 iterations of 2 */ \
\
xlast_two_fft_levels; /* jmp */

/* Do inverse FFT levels 7,8 */ \
/* On input the 64-byte cache lines hold these data values: */ \
/*	0	1	2	3	4	5	6	7 */ \
/*	8	9	10	11	... */ \
/*	... */ \
/* On output the 64-byte cache lines hold these data values: */ \
/*	0	1	8	9	16	17	24	25 */ \
/*	2	... */ \
/*	4	... */ \
/*	6	... */ \
/*	32	... */ \
/*	... */ \

/* Do 2 eight_reals_unfft_2 macros */ \
/* Do 2 nop_two_two_complex_unfft_2 macros */ \
/*	distance between fft data elements is 4 */ \
/*	do 1 macro processing 32 data values */ \


#define cmn1024(RCX1) \
x4cl_eight_reals_unfft_2(rsi, 4*64, 64, 2*64); \
\
/* Do 76/92/108/124 four_complex_unfft macros */ \
/*	distance between fft data elements is 4 */ \
/*	do 19/23/27/31 macros each processing 32 data values */ \
\
rdi = (uintptr_t)g->sincos2;/* Load sin/cos pointer */ \
for(unsigned int loop = g->count2 & 0xFF; loop != 0; loop--){ \
	x4cl_four_complex_unfft(rsi, 4*64, 64, 2*64, rdi); \
	rdi += XMM_SCD;		/* Next sine/cosine pointer */ \
} \
rsi = (uintptr_t)g->DESTARG;/* Restore source pointer */ \
\
/* Do inverse FFT level 6 */ \
/* On input the 64-byte cache lines hold these data values: */ \
/*	0	1	8	9	16	17	24	25 */ \
/*	2	... */ \
/*	4	... */ \
/*	6	... */ \
/*	32	... */ \
/*	... */ \
/* On output the 64-byte cache lines hold these data values: */ \
/*	0	1	16	17	32	33	48	49 */ \
/*	2	... */ \
/*	...	... */ \
/*	14	... */ \
/*	64	... */ \
/*	... */ \
\
/* Do 8 eight_reals_unfft_1 macros */ \
/*	distance between fft data elements is 16 */ \
/*	do 4 macros each processing 16 data values */ \
\
for(unsigned int loop=4; loop != 0; loop--){ \
	x2cl_eight_reals_unfft_1(rsi, 64, 4*64); \
} \
\
/* Do 72/88/104/120 four_complex_unfft macros */ \
/*	distance between fft data elements is 16 */ \
/*	do 36/44/52/60 macros each processing 16 data values */ \
\
rdi = (uintptr_t)g->sincos3;		/* Load sin/cos pointer */ \
rsi += -4*64+8*64;		/* Load source pointer */ \
for(unsigned int loop = g->count3 & 0xFF; loop != 0; loop--){ \
	for(unsigned int loopB=4; loopB != 0; loopB--){ \
		x2cl_two_complex_unfft(rsi, 64, 4*64); \
	} \
	rsi += -4*64+8*64;		/* Next source pointer */ \
	rdi += XMM_SCD;		/* Next sine/cosine pointer */ \
} \
rsi = (uintptr_t)g->DESTARG;/* Restore source pointer */ \
\
/* Do inverse FFT levels 4,5 */ \
/* On input the 64-byte cache lines hold these data values: */ \
/*	0	1	16	17	32	33	48	49 */ \
/*	2	... */ \
/*	...	... */ \
/*	14	... */ \
/*	64	... */ \
/*	... */ \
/* On output the 64-byte cache lines hold these data values: */ \
/*	0	1	64	65	128	129	192	193 */ \
/*	2	... */ \
/*	...	... */ \
/*	62	... */ \
/*	256	... */ \
/*	... */ \
\
/* Do 64/96 four_complex_unfft macros */ \
/*	distance between fft data elements is 32 */ \
/*	do 16/24 macros each processing 32 data values */ \
\
rsi += RCX1;		/* Next source pointer */ \
rdi = (uintptr_t)g->sincos4;		/* Load sin/cos pointer */ \
/* 2/3 iterations of 8 */ \
for(unsigned int loop = g->count5 & 0xFF; loop != 0; loop--){ \
	for(unsigned int loop=8; loop != 0; loop--){ \
		x4cl_four_complex_unfft(rsi, 64, 8*64, 16*64, rdi); \
	} \
	rsi += -8*64+32*64;	/* Next source pointer */ \
	rdi += XMM_SCD;		/* Next sine/cosine pointer */ \
} \
rsi = (uintptr_t)g->DESTARG;/* Next source pointer */ \
\
/* This ends the common inverse FFT code */ \
\
/* ret; push_amt = 0; */ \



/* Perform a 1024-element negacyclic FFT */


#define xfft1024p(type) \
\
/* Do a multiply with pre-FFTed inputs */ \
\
if(g->ffttype == 4) { \
	xmiddle_4p(g); \
}else{ \
	\
	/* Do FFT levels 1,2,3 */ \
	/* */ \
	/* On input the 64-byte cache lines hold these data values: */ \
	/*	0	256	1	257	512	768	513	769 */ \
	/*	2	... */ \
	/*	... */ \
	/*	254 */ \
	/* On output the 64-byte cache lines hold these data values: */ \
	/*	0	1	128	129	256	257	384	385 */ \
	/*	2	... */ \
	/*	... */ \
	/*	126	... */ \
	/*	512	... */ \
	/*	... */ \
	\
	/* Do 128 four_complex_first_fft macros */ \
	/*	distance between fft data elements is 128 */ \
	/*	do 64 macros each processing 16 data values */ \
	\
	rdi = (uintptr_t)g->plus1_premults;/* Address of premultiplier table */ \
	for(unsigned int loop=64; loop != 0; loop--){ \
		s2cl_four_complex_first_fft(rsi, 64, 64*64); \
		rdi += 128;		/* Next premultiplier pointer */ \
	} \
	rsi += -64*64;		/* Restore source pointer */ \
	\
	/* Do FFT level 4 */ \
	/* */ \
	/* On input the 64-byte cache lines hold these data values: */ \
	/*	0	1	128	129	256	257	384	385 */ \
	/*	2	... */ \
	/*	... */ \
	/*	126	... */ \
	/*	128	... */ \
	/*	... */ \
	/* On output the 64-byte cache lines hold these data values: */ \
	/*	0	1	64	65	128	129	192	193 */ \
	/*	2	3	... */ \
	/*	... */ \
	/*	62 */ \
	/*	256	... */ \
	/*	... */ \
	\
	/* Do 128 two_two_complex_fft macros */ \
	/*	distance between fft data elements is 64 */ \
	/*	do 64 macros each processing 16 data values */ \
	\
	rdi = (uintptr_t)g->sincos4;		/* Load sin/cos pointer */ \
	for(unsigned int loop=2; loop != 0; loop--){ \
		for(unsigned int loopB=32; loopB != 0; loopB--){ \
			x2cl_two_complex_fft(rsi, 64, 32*64, rdi); \
		}\
		rsi += -32*64+64*64;	/* Restore source pointer */ \
		rdi += 2*XMM_SCD;		/* Next sine/cosine pointer */ \
	} \
	rsi += -2*64*64;		/* Restore source pointer */ \
	\
	/* Join common code to perform rest of the FFT */ \
	/* Do FFT levels 5,6 */ \
	/* Do FFT levels 7,8 */ \
	/* Do FFT levels 9,10 */ \
	/*	do 48/64 macros each processing 16 data values */ \
	/*	mov	eax, count1		;; 48/64 iterations */ \
	\
	xlast_four_fft_levels_p; \
} \
cmn1024p; \
\
/* Do inverse FFT levels 1,2,3 */ \
/* On input the 64-byte cache lines hold these data values: */ \
/*	0	1	64	65	128	129	192	193 */ \
/*	2	... */ \
/*	...	... */ \
/*	62	... */ \
/*	256	... */ \
/*	... */ \
/* On output the 64-byte cache lines hold these data values: */ \
/*	0	1	256	257	512	513	768	769 */ \
/*	2	... */ \
/*	... */ \
/*	254 */ \
\
/* Do 128 four_complex_last_unfft macros */ \
/*	distance between fft data elements is 128 */ \
/*	do 32 macros each processing 32 data values */ \
\
rdi = (uintptr_t)g->plus1_premults;/* Address of premultiplier table */ \
for(unsigned int loop=32; loop != 0; loop--){ \
	x4cl_four_complex_last_unfft(rsi, 64, 32*64, 64*64, 32*128); \
	rdi += 128;		/* Next premultiplier pointer */ \
} \
\
xfft_3_ret; \
\
/* Common FFT code for negacyclic length 768 and 1024 */ \
\
/* Do inverse FFT levels 7,8 */ \
/* On input the 64-byte cache lines hold these data values: */ \
/*	0	1	2	3	4	5	6	7 */ \
/*	8	9	10	11	... */ \
/*	... */ \
/* On output the 64-byte cache lines hold these data values: */ \
/*	0	1	8	9	16	17	24	25 */ \
/*	2	... */ \
/*	4	... */ \
/*	6	... */ \
/*	32	... */ \
/*	... */ \

/* Do 96/128 four_complex_unfft macros */ \
/*	distance between fft data elements is 4 */ \
/*	do 24/32 macros each processing 32 data values */ \

#define cmn1024p \
/* push_amt = SZPTR; */ \
rdi = (uintptr_t)g->sincos2;/* Load sin/cos pointer */ \
for(unsigned int loop = g->count2 & 0xFF; loop != 0; loop--){ \
	x4cl_four_complex_unfft(rsi, 4*64, 64, 2*64, rdi); \
	rdi += XMM_SCD;		/* Next sine/cosine pointer */ \
} \
rsi = (uintptr_t)g->DESTARG;/* Restore source pointer */ \
\
/* Do inverse FFT levels 5,6 */ \
/* On input the 64-byte cache lines hold these data values: */ \
/*	0	1	8	9	16	17	24	25 */ \
/*	2	... */ \
/*	4	... */ \
/*	6	... */ \
/*	32	... */ \
/*	... */ \
/* On output the 64-byte cache lines hold these data values: */ \
/*	0	1	32	33	64	65	96	97 */ \
/*	2	... */ \
/*	...	... */ \
/*	30	... */ \
/*	128	... */ \
/*	... */ \
\
/* Do 96/128 four_complex_unfft macros */ \
/*	distance between fft data elements is 16 */ \
/*	do 24/32 macros each processing 32 data values */ \
\
rdi = (uintptr_t)g->sincos3;		/* Load sin/cos pointer */ \
/* 6/8 iterations of 4 */ \
for(unsigned int loop = g->count3 & 0xFF; loop != 0; loop--){ \
	for(unsigned int loopB=4; loopB != 0; loopB--){ \
		x4cl_four_complex_unfft(rsi, 64, 4*64, 8*64, rdi); \
	} \
	rsi += -4*64+16*64;	/* Next source pointer */ \
	rdi += XMM_SCD;		/* Next sine/cosine pointer */ \
} \
rsi = (uintptr_t)g->DESTARG;/* Restore source pointer */ \
\
/* Do inverse FFT level 4 */ \
/* On input the 64-byte cache lines hold these data values: */ \
/*	0	1	32	33	64	65	96	97 */ \
/*	2	... */ \
/*	...	... */ \
/*	30	... */ \
/*	128	... */ \
/*	... */ \
/* On output the 64-byte cache lines hold these data values: */ \
/*	0	1	64	65	128	129	192	193 */ \
/*	2	... */ \
/*	...	... */ \
/*	62	... */ \
/*	256	... */ \
/*	... */ \
\
/* Do 96/128 two_two_complex_unfft macros */ \
/*	distance between fft data elements is 64 */ \
/*	do 48/64 macros each processing 16 data values */ \
\
rdi = (uintptr_t)g->sincos4;		/* Load sin/cos pointer */ \
for(unsigned int loop = (g->count3 >> 8) & 0xFF; loop != 0; loop--){ \
	for(unsigned int loopB=16; loopB != 0; loopB--){ \
		x2cl_two_complex_unfft(rsi, 64, 16*64); \
	} \
	rsi += -16*64+32*64;	/* Next source pointer */ \
	rdi += XMM_SCD;		/* Next sine/cosine pointer */ \
} \
rsi = (uintptr_t)g->DESTARG;/* Restore source pointer */ \
\
/* This ends the common inverse FFT code */ \
\
/* ret; push_amt = 0; */



/* Perform a 1280-element FFT. */


#define xfft1280(type) \
\
/* Do a multiply with pre-FFTed inputs */ \
\
if(g->ffttype == 4) { \
	xmiddle_4(g); \
}else{ \
	\
	/* Do FFT levels 1,2,3 */ \
	/* Values 0-1279 is real data. */ \
	/* */ \
	/* On input the 64-byte cache lines hold these data values: */ \
	/*	0	64	1	65	640	704	641	705 */ \
	/*	2	... */ \
	/*	... */ \
	/*	62 */ \
	/*	128	384	129	385	768	1024	769	1025 */ \
	/*	130 */ \
	/*	... */ \
	/*	382 */ \
	/* On output the 64-byte cache lines hold these data values (in +1 format): */ \
	/*	0	64	128	192 */ \
	/*	2	... */ \
	/*	... */ \
	/*	62 */ \
	/*	256	512	768	1024 */ \
	/*	258 */ \
	/*	... */ \
	/*	510 */ \
	\
	/* Do 256 five_reals_first_fft macros */ \
	/*	distance between fft data elements is 256 */ \
	/*	do 32 macros each processing 40 data values */ \
	\
	xfive_reals_fft_preload;		/* Preload const1s for five real macros */ \
	for(unsigned int loop=32; loop != 0; loop--){ \
		s5cl_five_reals_first_fft(rsi, 64, 32*64); \
	} \
	rsi += -32*64;		/* Restore source pointer */ \
	\
	/* Do FFT levels 4,5 */ \
	/* Values 0-255 is real data, 256-1279 is complex data. */ \
	\
	/* Do 32 eight_reals_fft_2 macros */ \
	/*	distance between fft data elements is 64 */ \
	/*	do 32 macros each processing 8 data values */ \
	\
	for(unsigned int loop=32; loop != 0; loop--){ \
		x1cl_half_eight_reals_fft_2(rsi, 64); \
	} \
	\
	/* Do 64 four_complex_fft macros */ \
	/*	distance between fft data elements is 64 */ \
	/*	do 32 macros each processing 32 data values */ \
	\
	rdi = (uintptr_t)g->sincos4;		/* Load sin/cos pointer */ \
	for(unsigned int loop=32; loop != 0; loop--){ \
		x4cl_four_complex_fft(rsi, 64, 32*64, 64*64, rdi); \
	} \
	\
	/* Join common code to perform rest of the FFT */ \
	\
	xlast_four_fft_levels; \
} \
/*rcx = 32*64;*/ \
cmn2048(32*64); \
\
/* Finish the inverse FFT levels 4,5 */ \
\
/* Do 32 eight_reals_unfft_2 macros */ \
/*	distance between fft data elements is 64 */ \
/*	do 16 macros each processing 16 data values */ \
\
for(unsigned int loop=16; loop != 0; loop--){ \
	x2cl_half_eight_reals_unfft_2(rsi, 64, 16*64); \
} \
rsi += -16*64;		/* Restore source pointer */ \
\
/* Do inverse FFT levels 1,2,3 */ \
/* On output the 64-byte cache lines hold these data values: */ \
/*	0	64	640	704 */ \
/*	2	... */ \
/*	... */ \
/*	62 */ \
/*	128	384	768	1024 */ \
/*	130 */ \
/*	... */ \
/*	382 */ \
\
/* Do 256 five_reals_last_unfft macros */ \
/*	distance between fft data elements is 256 */ \
/*	do 32 macros each processing 40 data values */ \
\
xfive_reals_unfft_preload;	/* Preload const1s for five real macros */ \
for(unsigned int loop=32; loop != 0; loop--){ \
	x5cl_five_reals_last_unfft(rsi, 64, 32*64); \
} \
\
xfft_3_ret; \
\
\
\
/* Perform a 1536-element FFT. */ \
\

#define xfft1536(type) \
\
/* Do a multiply with pre-FFTed inputs */ \
\
if(g->ffttype == 4) { \
	xmiddle_4(g); \
}else{ \
	\
	/* Do FFT levels 1,2,3 */ \
	/* Values 0-1535 is real data. */ \
	/* */ \
	/* On input the 64-byte cache lines hold these data values: */ \
	/*	0	384	1	385	768	1152	769	1153 */ \
	/*	2 */ \
	/*	... */ \
	/*	382 */ \
	/* On output the 64-byte cache lines hold these data values (in +1 format): */ \
	/*	0	128	256	384 */ \
	/*	2	.. */ \
	/*	... */ \
	/*	126 */ \
	/*	512	768	1024	1280 */ \
	/*	514 */ \
	/*	... */ \
	/*	766 */ \
	\
	/* Do 256 six_reals_first_fft macros */ \
	/*	distance between fft data elements is 256 */ \
	/*	do 64 macros each processing 24 data values */ \
	\
	xsix_reals_fft_preload; \
	for(unsigned int loop=64; loop != 0; loop--){ \
		s3cl_six_reals_first_fft(rsi, 64, 64*64); \
	} \
	rsi += -64*64;		/* Restore source pointer */ \
	\
	/* Do FFT levels 4,5 */ \
	/* Values 0-255 is real data, 256-511 is semi-real data, 512-1535 is */ \
	/* complex data. */ \
	/* */ \
	/* On input the 64-byte cache lines hold these data values: */ \
	/*	0	128	256	384 */ \
	/*	2	.. */ \
	/*	... */ \
	/*	126 */ \
	/*	512	768	1024	1280 */ \
	/*	514 */ \
	/*	... */ \
	/*	766 */ \
	/* On output the 64-byte cache lines hold these data values: */ \
	/*	0	1	64	65	128	129	192	193 */ \
	/*	2	... */ \
	/*	... */ \
	/*	62	... */ \
	/*	256	... */ \
	/*	... */ \
	\
	/* Do 32 eight_reals_fft_2 macros */ \
	/* Do 32 nop_two_two_complex_fft_2 macros */ \
	/* Do 64 four_complex_fft macros */ \
	/*	distance between fft data elements is 64 */ \
	/*	do 32 macros each processing 16 data values */ \
	\
	for(unsigned int loop=32; loop != 0; loop--){ \
		x2cl_eight_reals_fft_2(rsi, 64, 32*64); \
	} \
	rsi += -32*64+64*64;	/* Next source pointer */ \
	\
	/* Do 64 four_complex_fft macros */ \
	/*	distance between fft data elements is 64 */ \
	/*	do 32 macros each processing 32 data values */ \
	\
	rdi = (uintptr_t)g->sincos4;		/* Load sin/cos pointer */ \
	for(unsigned int loop=32; loop != 0; loop--){ \
		x4cl_four_complex_fft(rsi, 64, 32*64, 64*64, rdi); \
	} \
	\
	/* Join common code to perform rest of the FFT */ \
	\
	xlast_four_fft_levels; \
} \
/*rcx = 64*64;*/ \
cmn2048(64*64); \
\
/* Finish inverse FFT levels 4,5 */ \
\
/* Do 32 eight_reals_unfft_2 macros */ \
/* Do 32 nop_two_two_complex_unfft_2 macros */ \
/*	distance between fft data elements is 64 */ \
/*	do 16 macros each processing 32 data values */ \
\
for(unsigned int loop=16; loop != 0; loop--){ \
	x4cl_eight_reals_unfft_2(rsi, 64, 16*64, 32*64); \
} \
rsi += -16*64;		/* Restore source pointer */ \
\
/* Do inverse FFT levels 1,2,3 */ \
/* On input the 64-byte cache lines hold these data values: */ \
/*	0	1	128	129	256	257	384	385 */ \
/*	2	... */ \
/*	...	... */ \
/*	126	... */ \
/*	512	... */ \
/*	... */ \
/* On output the 64-byte cache lines hold these data values: */ \
/*	0	384	768	1152 */ \
/*	2 */ \
/*	... */ \
/*	382 */ \
\
/* Do 256 six_reals_last_unfft macros */ \
/*	distance between fft data elements is 256 */ \
/*	do 64 macros each processing 24 data values */ \
\
xsix_reals_unfft_preload; \
for(unsigned int loop=64; loop != 0; loop--){ \
	x3cl_six_reals_last_unfft(rsi, 64, 64*64); \
} \
\
xfft_3_ret; \
\
\
\
/* Perform a 1536-element negacyclic FFT */ \
\

#define xfft1536p(type) \
\
/* Do a multiply with pre-FFTed inputs */ \
\
if(g->ffttype == 4) { \
	xmiddle_4p(g); \
}else{ \
	\
	/* Do FFT levels 1,2,3 */ \
	/* */ \
	/* On input the 64-byte cache lines hold these data values: */ \
	/*	0	384	1	385	768	1152	769	1153 */ \
	/*	2 */ \
	/*	... */ \
	/*	382 */ \
	/* On output the 64-byte cache lines hold these data values (in +1 format): */ \
	/*	0	128	256	384 */ \
	/*	2	.. */ \
	/*	... */ \
	/*	126 */ \
	/*	512	768	1024	1280 */ \
	/*	514 */ \
	/*	... */ \
	/*	766 */ \
	\
	/* Do 256 three_complex_first_fft macros */ \
	/*	distance between fft data elements is 256 */ \
	/*	do 64 macros each processing 24 data values */ \
	\
	rdi = (uintptr_t)g->plus1_premults;/* Address of premultiplier table */ \
	for(unsigned int loop=64; loop != 0; loop--){ \
		s3cl_three_complex_first_fft(rsi, 64, 64*64); \
		rdi += 192;		/* Next premultiplier pointer */ \
	} \
	rsi += -64*64;		/* Restore source pointer */ \
	\
	/* Do FFT levels 4,5 */ \
	/* */ \
	/* On input the 64-byte cache lines hold these data values: */ \
	/*	0	128	256	384 */ \
	/*	2	.. */ \
	/*	... */ \
	/*	126 */ \
	/*	512	768	1024	1280 */ \
	/*	514 */ \
	/*	... */ \
	/*	766 */ \
	/* On output the 64-byte cache lines hold these data values: */ \
	/*	0	1	64	65	128	129	192	193 */ \
	/*	2	... */ \
	/*	... */ \
	/*	62	... */ \
	/*	256	... */ \
	/*	... */ \
	\
	/* Do 96 four_complex_fft macros */ \
	/*	distance between fft data elements is 64 */ \
	/*	do 32 macros each processing 16 data values */ \
	/*	do 32 macros each processing 32 data values */ \
	\
	rdi = (uintptr_t)g->sincos4;		/* Load sin/cos pointer */ \
	for(unsigned int loop=32; loop != 0; loop--){ \
		x2cl_four_complex_fft(rsi, 64, 32*64); \
	} \
	rsi += -32*64+64*64;	/* Next source pointer */ \
	rdi += XMM_SCD;		/* Next sin/cos pointer */ \
	for(unsigned int loop=32; loop != 0; loop--){ \
		x4cl_four_complex_fft(rsi, 64, 32*64, 64*64, rdi); \
	} \
	rsi += -32*64-64*64;	/* Restore source pointer */ \
	\
	/* Join common code to perform rest of the FFT */ \
	\
	xlast_four_fft_levels_p; \
} \
cmn2048p; \
\
/* Do inverse FFT levels 1,2,3 */ \
/* On input the 64-byte cache lines hold these data values: */ \
/*	0	1	128	129	256	257	384	385 */ \
/*	2	... */ \
/*	...	... */ \
/*	126	... */ \
/*	512	... */ \
/*	... */ \
/* On output the 64-byte cache lines hold these data values: */ \
/*	0	384	768	1152 */ \
/*	2 */ \
/*	... */ \
/*	382 */ \
\
/* Do 256 three_complex_last_unfft macros */ \
/*	distance between fft data elements is 256 */ \
/*	do 64 macros each processing 24 data values */ \
\
rdi = (uintptr_t)g->plus1_premults;/* Address of premultiplier table */ \
for(unsigned int loop=64; loop != 0; loop--){ \
	x3cl_three_complex_last_unfft(rsi, 64, 64*64); \
	rdi += 192;		/* Next premultiplier pointer */ \
} \
\
xfft_3_ret; \
\
\
\
/* Perform a 1792-element FFT. */ \
\

#define xfft1792(type) \
\
/* Do a multiply with pre-FFTed inputs */ \
\
if(g->ffttype == 4) { \
	xmiddle_4(g); \
}else{ \
	\
	/* Do FFT levels 1,2,3 */ \
	/* Values 0-1791 is real data. */ \
	/* */ \
	/* On input the 64-byte cache lines hold these data values: */ \
	/*	0	64	1	65	896	960	897	961 */ \
	/*	2	... */ \
	/*	... */ \
	/*	62 */ \
	/*	128	256	129	257	1024	1152	1025	1153 */ \
	/*	130 */ \
	/*	... */ \
	/*	254 */ \
	/*	384	640	385	641	1280	1536	1281	1537 */ \
	/*	386 */ \
	/*	... */ \
	/*	638 */ \
	/* On output the 64-byte cache lines hold these data values (in +1 format): */ \
	/*	0	64	128	192 */ \
	/*	2	... */ \
	/*	... */ \
	/*	62 */ \
	/*	256	384	512	640 */ \
	/*	258 */ \
	/*	... */ \
	/*	382 */ \
	/*	768	1024	1280	1536 */ \
	/*	770 */ \
	/*	... */ \
	/*	1022 */ \
	\
	/* Do 256 seven_reals_first_fft macros */ \
	/*	distance between fft data elements is 256 */ \
	/*	do 32 macros each processing 56 data values */ \
	\
	xseven_reals_fft_preload; \
	for(unsigned int loop=32; loop != 0; loop--){ \
		s7cl_seven_reals_first_fft(rsi, 64, 32*64); \
	} \
	rsi += -32*64;		/* Restore source pointer */ \
	\
	/* Do FFT levels 4,5 */ \
	/* Values 0-255 is real data, 256-1791 is complex data. */ \
	\
	/* Do 32 eight_reals_fft_2 macros */ \
	/*	distance between fft data elements is 64 */ \
	/*	do 32 macros each processing 8 data values */ \
	\
	for(unsigned int loop=32; loop != 0; loop--){ \
		x1cl_half_eight_reals_fft_2(rsi, 64); \
	} \
	\
	/* Do 64 four_complex_fft macros */ \
	/*	distance between fft data elements is 64 */ \
	/*	do 32 macros each processing 16 data values */ \
	\
	rdi = (uintptr_t)g->sincos4;		/* Load sin/cos pointer */ \
	for(unsigned int loop=32; loop != 0; loop--){ \
		x2cl_four_complex_fft(rsi, 64, 32*64); \
	} \
	rsi += -32*64+64*64;	/* Next source pointer */ \
	\
	/* Do 128 four_complex_fft macros */ \
	/*	distance between fft data elements is 64 */ \
	/*	do 32 macros each processing 32 data values */ \
	\
	rdi += XMM_SCD;		/* Next sine/cosine pointer */ \
	for(unsigned int loop=32; loop != 0; loop--){ \
		x4cl_four_complex_fft(rsi, 64, 32*64, 64*64, rdi); \
	} \
	\
	/* Do FFT levels 6,7 */ \
	/*	do 32/40/48/56 macros each processing 32 data values */ \
	/*	mov	ax, WORD PTR count4+2	;; 4/5/6/7 iterations of 8 */ \
	/* Do FFT levels 8,9 */ \
	/*	do 38/46/54/62 macros each processing 32 data values */ \
	/*	mov	ax, WORD PTR count4	;; 3/7 then 16/24 iterations of 2 */ \
	\
	xlast_four_fft_levels; \
} \
/*rcx = 32*64*/; \
cmn2048(32*64); \
\
/* Finish inverse FFT levels 4,5 */ \
\
/* Do 32 eight_reals_unfft_2 macros */ \
/*	distance between fft data elements is 64 */ \
/*	do 16 macros each processing 16 data values */ \
\
for(unsigned int loop=16; loop != 0; loop--){ \
	x2cl_half_eight_reals_unfft_2(rsi, 64, 16*64); \
} \
rsi += -16*64;		/* Restore source pointer */ \
\
/* Do inverse FFT levels 1,2,3 */ \
/* On input the 64-byte cache lines hold these data values: */ \
/*	0	64	128	192 */ \
/*	2	... */ \
/*	... */ \
/*	62 */ \
/*	256	384	512	640 */ \
/*	... */ \
/* On output the 64-byte cache lines hold these data values: */ \
/*	0	64	896	960 */ \
/*	2	... */ \
/*	... */ \
/*	62 */ \
/*	128	256	1024	1152 */ \
/*	130 */ \
/*	... */ \
/*	254 */ \
/*	384	640	1280	1536 */ \
/*	386 */ \
/*	... */ \
/*	638 */ \
\
/* Do 256 seven_reals_last_unfft macros */ \
/*	distance between fft data elements is 256 */ \
/*	do 32 macros each processing 56 data values */ \
\
xseven_reals_unfft_preload; \
for(unsigned int loop=32; loop != 0; loop--){ \
	x7cl_seven_reals_last_unfft(rsi, 64, 32*64); \
} \
\
xfft_3_ret; \
\
\
\
/* Perform a 2048-element FFT. */ \
\

#define xfft2048(type) \
\
/* Do a multiply with pre-FFTed inputs */ \
\
if(g->ffttype == 4) { \
	xmiddle_4(g); \
}else{ \
	\
	/* Do FFT levels 1,2,3 */ \
	/* Values 0-2047 is real data. */ \
	/* */ \
	/* On input the 64-byte cache lines hold these data values: */ \
	/*	0	512	1	513	1024	1536	1025	1537 */ \
	/*	2	... */ \
	/*	... */ \
	/*	510 */ \
	/* On output the 64-byte cache lines hold these data values: */ \
	/*	0	1	256	257	512	513	768	769 */ \
	/*	2	... */ \
	/*	... */ \
	/*	254	... */ \
	/*	1K	... */ \
	/*	... */ \
	\
	/* Do 256 eight_reals_first_fft macros */ \
	/*	distance between fft data elements is 256 */ \
	/*	do 128 macros each processing 16 data values */ \
	\
	for(unsigned int loop=128; loop != 0; loop--){ \
		s2cl_eight_reals_first_fft(rsi, 64, 128*64); \
	} \
	rsi += -128*64;		/* Restore source pointer */ \
	\
	/* Do FFT levels 4,5 */ \
	/* Values 0-255 is real data, 256-511 is semi-real data, 512-2047 is */ \
	/* complex data. */ \
	/* */ \
	/* On input the 64-byte cache lines hold these data values: */ \
	/*	0	1	256	257	512	513	768	769 */ \
	/*	2	... */ \
	/*	... */ \
	/*	254	... */ \
	/*	1K	... */ \
	/*	... */ \
	/* On output the 64-byte cache lines hold these data values: */ \
	/*	0	1	64	65	128	129	192	193 */ \
	/*	2	... */ \
	/*	... */ \
	/*	62	... */ \
	/*	256	... */ \
	/*	... */ \
	\
	/* Do 32 eight_reals_fft_2 macros */ \
	/* Do 32 nop_two_two_complex_fft_2 macros */ \
	/* Do 64 four_complex_fft macros */ \
	/*	distance between fft data elements is 64 */ \
	/*	do 32 macros each processing 32 data values */ \
	\
	rdi = (uintptr_t)g->sincos4;		/* Load sin/cos pointer */ \
	for(unsigned int loop=32; loop != 0; loop--){ \
		x4cl_eight_reals_fft_2(rsi, 64, 32*64, 64*64, rdi); \
	} \
	rsi += -32*64+128*64;	/* Next source pointer */ \
	\
	/* Do 128 four_complex_fft macros */ \
	/*	distance between fft data elements is 64 */ \
	/*	do 32 macros each processing 32 data values */ \
	\
	rdi += XMM_SCD;		/* Next sine/cosine pointer */ \
	for(unsigned int loop=32; loop != 0; loop--){ \
		x4cl_four_complex_fft(rsi, 64, 32*64, 64*64, rdi); \
	} \
	\
	/* Do FFT levels 6,7 */ \
	/*	do 32/40/48/56 macros each processing 32 data values */ \
	/*	mov	ax, WORD PTR count4+2	;; 4/5/6/7 iterations of 8 */ \
	/* Do FFT levels 8,9 */ \
	/*	do 38/46/54/62 macros each processing 32 data values */ \
	/*	mov	ax, WORD PTR count4	;; 3/7 then 16/24 iterations of 2 */ \
	/* Do FFT levels 10,11 */ \
	/*	do 79/95/111/127 macros each processing 16 data values */ \
	/*	mov	eax, count1		;; 15/31 then 64/96 iterations */ \
	\
	xlast_four_fft_levels; \
} \
/*rcx = 64*64*/; \
cmn2048(64*64); \
\
/* Finish the inverse FFT levels 4,5 */ \
\
/* Do 32 eight_reals_unfft_2 macros */ \
/* Do 32 nop_two_two_complex_unfft_2 macros */ \
/*	distance between fft data elements is 64 */ \
/*	do 16 macros each processing 32 data values */ \
\
for(unsigned int loop=16; loop != 0; loop--){ \
	x4cl_eight_reals_unfft_2(rsi, 64, 16*64, 32*64); \
} \
rsi += -16*64;		/* Restore source pointer */ \
\
/* Do inverse FFT levels 1,2,3 */ \
/* On input the 64-byte cache lines hold these data values: */ \
/*	0	1	128	129	256	257	384	385 */ \
/*	2	... */ \
/*	...	... */ \
/*	126	... */ \
/*	512	... */ \
/*	... */ \
/* On output the 64-byte cache lines hold these data values: */ \
/*	0	1	512	513	1024	1025	1536	1537 */ \
/*	2	... */ \
/*	... */ \
/*	510 */ \
\
/* Do 256 eight_reals_last_unfft macros */ \
/*	distance between fft data elements is 256 */ \
/*	do 64 macros each processing 32 data values */ \
\
for(unsigned int loop=64; loop != 0; loop--){ \
	x4cl_eight_reals_last_unfft(rsi, 64, 64*64, 128*64); \
} \
\
xfft_3_ret; \

/* Common inverse FFT code for lengths 1280, 1536, 1792, 2048 */ \

/* Do inverse FFT levels 8,9 */ \
/* On input the 64-byte cache lines hold these data values: */ \
/*	0	1	2	3	4	5	6	7 */ \
/*	8	9	10	11	... */ \
/*	... */ \
/* On output the 64-byte cache lines hold these data values: */ \
/*	0	1	8	9	16	17	24	25 */ \
/*	2	... */ \
/*	4	... */ \
/*	6	... */ \
/*	32	... */ \
/*	... */ \

/* Do 2 eight_reals_unfft_2 macros */ \
/* Do 2 nop_two_two_complex_unfft_2 macros */ \
/*	distance between fft data elements is 4 */ \
/*	do 1 macro processing 32 data values */ \

#define cmn2048(RCX1) \
/* push_amt = SZPTR; */ \
x4cl_eight_reals_unfft_2(rsi, 4*64, 64, 2*64); \
\
/* Do 156/188/220/252 four_complex_unfft macros */ \
/*	distance between fft data elements is 4 */ \
/*	do 39/47/55/63 macros each processing 32 data values */ \
\
rdi = (uintptr_t)g->sincos2;/* Load sin/cos pointer */ \
for(unsigned int loop = g->count2 & 0xFF; loop != 0; loop--){ \
	x4cl_four_complex_unfft(rsi, 4*64, 64, 2*64, rdi); \
	rdi += XMM_SCD;		/* Next sine/cosine pointer */ \
} \
rsi = (uintptr_t)g->DESTARG;/* Restore source pointer */ \
\
/* Do inverse FFT levels 6,7 */ \
/* On input the 64-byte cache lines hold these data values: */ \
/*	0	1	8	9	16	17	24	25 */ \
/*	2	... */ \
/*	4	... */ \
/*	6	... */ \
/*	32	... */ \
/*	... */ \
/* On output the 64-byte cache lines hold these data values: */ \
/*	0	1	32	33	64	65	96	97 */ \
/*	2	... */ \
/*	...	... */ \
/*	30	... */ \
/*	128	... */ \
/*	... */ \
\
/* Do 8 eight_reals_unfft_2 macros */ \
/* Do 8 nop_two_two_complex_unfft_2 macros */ \
/*	distance between fft data elements is 16 */ \
/*	do 4 macros each processing 32 data values */ \
\
for(unsigned int loop=4; loop != 0; loop--){ \
	x4cl_eight_reals_unfft_2(rsi, 64, 4*64, 8*64); \
} \
\
/* Do 144/176/208/240 four_complex_unfft macros */ \
/*	distance between fft data elements is 16 */ \
/*	do 36/44/52/60 macros each processing 32 data values */ \
\
rdi = (uintptr_t)g->sincos3;		/* Load sin/cos pointer */ \
rsi += -4*64+16*64;	/* Load source pointer */ \
for(unsigned int loop = g->count3 & 0xFF; loop != 0; loop--){ \
	for(unsigned int loopB=4; loopB != 0; loopB--){ \
		x4cl_four_complex_unfft(rsi, 64, 4*64, 8*64, rdi); \
	} \
	rsi += -4*64+16*64;	/* Next source pointer */ \
	rdi += XMM_SCD;		/* Next sine/cosine pointer */ \
} \
rsi = (uintptr_t)g->DESTARG;/* Restore source pointer */ \
\
/* Do inverse FFT levels 4,5 */ \
/* On input the 64-byte cache lines hold these data values: */ \
/*	0	1	32	33	64	65	96	97 */ \
/*	2	... */ \
/*	...	... */ \
/*	30	... */ \
/*	128	... */ \
/*	... */ \
/* On output the 64-byte cache lines hold these data values: */ \
/*	0	1	128	129	256	257	384	385 */ \
/*	2	... */ \
/*	...	... */ \
/*	126	... */ \
/*	512	... */ \
/*	... */ \
\
/* Do 128/192 four_complex_unfft macros */ \
/*	distance between fft data elements is 64 */ \
/*	do 32/48 macros each processing 32 data values */ \
\
rsi += RCX1;		/* Next source pointer */ \
rdi = (uintptr_t)g->sincos4;		/* Load sin/cos pointer */ \
/* 2/3 iterations of 16 */ \
for(unsigned int loop = g->count5 & 0xFF; loop != 0; loop--){ \
	for(unsigned int loopB=16; loopB != 0; loopB--){ \
		x4cl_four_complex_unfft(rsi, 64, 16*64, 32*64, rdi); \
	} \
	rsi += -16*64+64*64;	/* Next source pointer */ \
	rdi += XMM_SCD;		/* Next sine/cosine pointer */ \
} \
rsi = (uintptr_t)g->DESTARG;/* Next source pointer */ \
\
/* This ends the common inverse FFT code */ \
\
/* ret; push_amt = 0; */ \
\
\
\
/* Perform a 2048-element negacyclic FFT */ \
\

#define xfft2048p(type) \
\
/* Do a multiply with pre-FFTed inputs */ \
\
if(g->ffttype == 4) { \
	xmiddle_4p(g); \
}else{ \
	\
	/* Do FFT levels 1,2,3 */ \
	/* */ \
	/* On input the 64-byte cache lines hold these data values: */ \
	/*	0	512	1	513	1024	1536	1025	1537 */ \
	/*	2	... */ \
	/*	... */ \
	/*	510 */ \
	/* On output the 64-byte cache lines hold these data values: */ \
	/*	0	1	256	257	512	513	768	769 */ \
	/*	2	... */ \
	/*	... */ \
	/*	254	... */ \
	/*	1K	... */ \
	/*	... */ \
	\
	/* Do 256 four_complex_first_fft macros */ \
	/*	distance between fft data elements is 256 */ \
	/*	do 128 macros each processing 16 data values */ \
	\
	rdi = (uintptr_t)g->plus1_premults;/* Address of premultiplier table */ \
	for(unsigned int loop=128; loop != 0; loop--){ \
		s2cl_four_complex_first_fft(rsi, 64, 128*64); \
		rdi += 128;		/* Next premultiplier pointer */ \
	} \
	rsi += -128*64;		/* Restore source pointer */ \
	\
	/* Do FFT levels 4,5 */ \
	/* */ \
	/* On input the 64-byte cache lines hold these data values: */ \
	/*	0	1	256	257	512	513	768	769 */ \
	/*	2	... */ \
	/*	... */ \
	/*	254	... */ \
	/*	1K	... */ \
	/*	... */ \
	/* On output the 64-byte cache lines hold these data values: */ \
	/*	0	1	64	65	128	129	192	193 */ \
	/*	2	... */ \
	/*	... */ \
	/*	62	... */ \
	/*	256	... */ \
	/*	... */ \
	\
	/* Do 256 four_complex_fft macros */ \
	/*	distance between fft data elements is 64 */ \
	/*	do 64 macros each processing 32 data values */ \
	\
	/* Do FFT levels 6,7 */ \
	/*	do 48/64 macros each processing 32 data values */ \
	/*	mov	ax, WORD PTR count4+2	;; 6/8 iterations of 8 */ \
	/* Do FFT levels 8,9 */ \
	/*	do 48/64 macros each processing 32 data values */ \
	/*	mov	ax, WORD PTR count4	;; 24/32 iterations of 2 */ \
	/* Do FFT levels 10,11 */ \
	/*	do 96/128 macros each processing 16 data values */ \
	/*	mov	eax, count1		;; 96/128 iterations */ \
	\
	xlast_six_fft_levels_p; \
} \
cmn2048p; \
\
/* Do inverse FFT levels 1,2,3 */ \
/* On input the 64-byte cache lines hold these data values: */ \
/*	0	1	128	129	256	257	384	385 */ \
/*	2	... */ \
/*	...	... */ \
/*	126	... */ \
/*	512	... */ \
/*	... */ \
/* On output the 64-byte cache lines hold these data values: */ \
/*	0	1	512	513	1024	1025	1536	1537 */ \
/*	2	... */ \
/*	... */ \
/*	510 */ \
\
/* Do 256 four_complex_last_unfft macros */ \
/*	distance between fft data elements is 256 */ \
/*	do 64 macros each processing 32 data values */ \
\
rdi = (uintptr_t)g->plus1_premults;/* Address of premultiplier table */ \
for(unsigned int loop=64; loop != 0; loop--){ \
	x4cl_four_complex_last_unfft(rsi, 64, 64*64, 128*64, 64*128); \
	rdi += 128;		/* Next premultiplier pointer */ \
} \
\
xfft_3_ret; \
\
/* Common inverse FFT for negacyclic lengths 1536, 2048 */ \
\
/* Do inverse FFT levels 8,9 */ \
/* On input the 64-byte cache lines hold these data values: */ \
/*	0	1	2	3	4	5	6	7 */ \
/*	8	9	10	11	... */ \
/*	... */ \
/* On output the 64-byte cache lines hold these data values: */ \
/*	0	1	8	9	16	17	24	25 */ \
/*	2	... */ \
/*	4	... */ \
/*	6	... */ \
/*	32	... */ \
/*	... */ \
\
/* Do 196/256 four_complex_unfft macros */ \
/*	distance between fft data elements is 4 */ \
/*	do 48/64 macros each processing 32 data values */ \

#define cmn2048p \
/* push_amt = SZPTR; */ \
rdi = (uintptr_t)g->sincos2;/* Load sin/cos pointer */ \
for(unsigned int loop = g->count2 & 0xFF; loop != 0; loop--){ \
	x4cl_four_complex_unfft(rsi, 4*64, 64, 2*64, rdi); \
	rdi += XMM_SCD;		/* Next sine/cosine pointer */ \
} \
rsi = (uintptr_t)g->DESTARG;/* Restore source pointer */ \
\
/* Do inverse FFT levels 6,7 */ \
/* On input the 64-byte cache lines hold these data values: */ \
/*	0	1	8	9	16	17	24	25 */ \
/*	2	... */ \
/*	4	... */ \
/*	6	... */ \
/*	32	... */ \
/*	... */ \
/* On output the 64-byte cache lines hold these data values: */ \
/*	0	1	32	33	64	65	96	97 */ \
/*	2	... */ \
/*	...	... */ \
/*	30	... */ \
/*	128	... */ \
/*	... */ \
\
/* Do 192/256 four_complex_unfft macros */ \
/*	distance between fft data elements is 16 */ \
/*	do 48/64 macros each processing 32 data values */ \
\
rdi = (uintptr_t)g->sincos3;		/* Load sin/cos pointer */ \
/* 12/16 iterations of 4 */ \
for(unsigned int loop = g->count3 & 0xFF; loop != 0; loop--){ \
	for(unsigned int loopB=4; loopB != 0; loopB--){ \
		x4cl_four_complex_unfft(rsi, 64, 4*64, 8*64, rdi); \
	} \
	rsi += -4*64+16*64;	/* Next source pointer */ \
	rdi += XMM_SCD;		/* Next sine/cosine pointer */ \
} \
rsi = (uintptr_t)g->DESTARG;/* Restore source pointer */ \
\
/* Do inverse FFT levels 4,5 */ \
/* On input the 64-byte cache lines hold these data values: */ \
/*	0	1	32	33	64	65	96	97 */ \
/*	2	... */ \
/*	...	... */ \
/*	30	... */ \
/*	128	... */ \
/*	... */ \
/* On output the 64-byte cache lines hold these data values: */ \
/*	0	1	128	129	256	257	384	385 */ \
/*	2	... */ \
/*	...	... */ \
/*	126	... */ \
/*	512	... */ \
/*	... */ \
\
/* Do 192/256 four_complex_unfft macros */ \
/*	distance between fft data elements is 64 */ \
/*	do 48/64 macros each processing 32 data values */ \
\
rdi = (uintptr_t)g->sincos4;		/* Load sin/cos pointer */ \
/* 3/4 iterations of 16 */ \
for(unsigned int loop = g->count5 & 0xFF; loop != 0; loop--){ \
	for(unsigned int loopB=16; loopB != 0; loopB--){ \
		x4cl_four_complex_unfft(rsi, 64, 16*64, 32*64, rdi); \
	} \
	rsi += -16*64+64*64;	/* Next source pointer */ \
	rdi += XMM_SCD;		/* Next sine/cosine pointer */ \
} \
rsi = (uintptr_t)g->DESTARG;/* Next source pointer */ \
\
/* This ends the common inverse FFT code */ \
\
/* ret; push_amt = 0; */ \



/* Perform a 2560-element FFT. */


#define xfft2560(type) \
\
/* Do a multiply with pre-FFTed inputs */ \
\
if(g->ffttype == 4) { \
	xmiddle_4(g); \
}else{ \
	\
	/* Do FFT levels 1,2,3 */ \
	/* Values 0-2559 is real data. */ \
	/* */ \
	/* On input the 64-byte cache lines hold these data values: */ \
	/*	0	128	1	129	1280	1408	1281	1409 */ \
	/*	2	... */ \
	/*	... */ \
	/*	126 */ \
	/*	256	768	257	769	1536	2048	1537	2049 */ \
	/*	258 */ \
	/*	... */ \
	/*	766 */ \
	/* On output the 64-byte cache lines hold these data values (in +1 format): */ \
	/*	0	128	256	384 */ \
	/*	2	... */ \
	/*	... */ \
	/*	126 */ \
	/*	512	1024	1536	2048 */ \
	/*	514 */ \
	/*	... */ \
	/*	1022 */ \
	\
	/* Do 512 five_reals_first_fft macros */ \
	/*	distance between fft data elements is 512 */ \
	/*	do 64 macros each processing 40 data values */ \
	\
	xfive_reals_fft_preload;		/* Preload const1s for five real macros */ \
	for(unsigned int loop=64; loop != 0; loop--){ \
		s5cl_five_reals_first_fft(rsi, 64, 64*64); \
	} \
	rsi += -64*64;		/* Restore source pointer */ \
	\
	/* Do FFT levels 4,5 */ \
	/* Values 0-511 is real data, 512-2559 is complex data. */ \
	\
	/* Do 64 eight_reals_fft_2 macros */ \
	/*	distance between fft data elements is 128 */ \
	/*	do 64 macros each processing 8 data values */ \
	\
	for(unsigned int loop=64; loop != 0; loop--){ \
		x1cl_half_eight_reals_fft_2(rsi, 64); \
	} \
	\
	/* Do 128 four_complex_fft macros */ \
	/*	distance between fft data elements is 128 */ \
	/*	do 64 macros each processing 32 data values */ \
	\
	rdi = (uintptr_t)g->sincos5;		/* Load sin/cos pointer */ \
	for(unsigned int loop=64; loop != 0; loop--){ \
		x4cl_four_complex_fft(rsi, 64, 64*64, 128*64, rdi); \
	} \
	rsi = (uintptr_t)g->DESTARG;/* Restore source pointer */ \
	\
	/* Join common code to perform rest of the FFT */ \
	\
	cmn4096_2; \
} \
/*rcx = 64*64;*/ \
cmn4096(64*64); \
\
/* Finish the inverse FFT levels 4,5 */ \
\
/* Do 64 eight_reals_unfft_2 macros */ \
/*	distance between fft data elements is 128 */ \
/*	do 32 macros each processing 16 data values */ \
\
for(unsigned int loop=32; loop != 0; loop--){ \
	x2cl_half_eight_reals_unfft_2(rsi, 64, 32*64); \
} \
rsi += -32*64;		/* Restore source pointer */ \
\
/* Do inverse FFT levels 1,2,3 */ \
/* On output the 64-byte cache lines hold these data values: */ \
/*	0	128	1280	1408 */ \
/*	2	... */ \
/*	... */ \
/*	126 */ \
/*	256	768	1536	2048 */ \
/*	258 */ \
/*	... */ \
/*	766 */ \
\
/* Do 512 five_reals_last_unfft macros */ \
/*	distance between fft data elements is 512 */ \
/*	do 64 macros each processing 40 data values */ \
\
xfive_reals_unfft_preload;	/* Preload const1s for five real macros */ \
for(unsigned int loop=64; loop != 0; loop--){ \
	x5cl_five_reals_last_unfft(rsi, 64, 64*64); \
} \
\
xfft_3_ret; \
\
\
\
/* Perform a 3072-element FFT. */ \
\

#define xfft3072(type) \
\
/* Do a multiply with pre-FFTed inputs */ \
\
if(g->ffttype == 4) { \
	xmiddle_4(g); \
}else{ \
	\
	/* Do FFT levels 1,2,3 */ \
	/* Values 0-3071 is real data. */ \
	/* */ \
	/* On input the 64-byte cache lines hold these data values: */ \
	/*	0	768	1	769	1536	2304	1537	2305 */ \
	/*	2 */ \
	/*	... */ \
	/*	766 */ \
	/* On output the 64-byte cache lines hold these data values (in +1 format): */ \
	/*	0	256	512	768 */ \
	/*	2	.. */ \
	/*	... */ \
	/*	254 */ \
	/*	1024	1532	2048	2560 */ \
	/*	1026 */ \
	/*	... */ \
	/*	1530 */ \
	\
	/* Do 512 six_reals_first_fft macros */ \
	/*	distance between fft data elements is 512 */ \
	/*	do 128 macros each processing 24 data values */ \
	\
	xsix_reals_fft_preload; \
	for(unsigned int loop=128; loop != 0; loop--){ \
		s3cl_six_reals_first_fft(rsi, 64, 128*64); \
	} \
	rsi += -128*64;		/* Restore source pointer */ \
	\
	/* Do FFT levels 4,5 */ \
	/* Values 0-511 is real data, 512-1023 is semi-real data, 1024-3071 is */ \
	/* complex data. */ \
	/* */ \
	/* On input the 64-byte cache lines hold these data values: */ \
	/*	0	256	512	768 */ \
	/*	2	.. */ \
	/*	... */ \
	/*	254 */ \
	/*	1024	1532	2048	2560 */ \
	/*	1026 */ \
	/*	... */ \
	/*	1530 */ \
	/* On output the 64-byte cache lines hold these data values: */ \
	/*	0	1	128	129	256	257	384	385 */ \
	/*	2	... */ \
	/*	... */ \
	/*	126	... */ \
	/*	512	... */ \
	/*	... */ \
	\
	/* Do 64 eight_reals_fft_2 macros */ \
	/* Do 64 nop_two_two_complex_fft_2 macros */ \
	/* Do 128 four_complex_fft macros */ \
	/*	distance between fft data elements is 128 */ \
	/*	do 64 macros each processing 16 data values */ \
	\
	for(unsigned int loop=64; loop != 0; loop--){ \
		x2cl_eight_reals_fft_2(rsi, 64, 64*64); \
	} \
	rsi += -64*64+128*64;	/* Next source pointer */ \
	\
	/* Do 128 four_complex_fft macros */ \
	/*	distance between fft data elements is 128 */ \
	/*	do 64 macros each processing 32 data values */ \
	\
	rdi = (uintptr_t)g->sincos5;		/* Load sin/cos pointer */ \
	for(unsigned int loop=64; loop != 0; loop--){ \
		x4cl_four_complex_fft(rsi, 64, 64*64, 128*64, rdi); \
	} \
	rsi = (uintptr_t)g->DESTARG;/* Restore source pointer */ \
	\
	/* Join common code to perform rest of the FFT */ \
	\
	cmn4096_2; \
} \
/*rcx = 128*64*/; \
cmn4096(128*64); \
\
/* Finish inverse FFT levels 4,5 */ \
\
/* Do 64 eight_reals_unfft_2 macros */ \
/* Do 64 nop_two_two_complex_unfft_2 macros */ \
/*	distance between fft data elements is 128 */ \
/*	do 32 macros each processing 32 data values */ \
\
for(unsigned int loop=32; loop != 0; loop--){ \
	x4cl_eight_reals_unfft_2(rsi, 64, 32*64, 64*64); \
} \
rsi += -32*64;		/* Restore source pointer */ \
\
/* Do inverse FFT levels 1,2,3 */ \
/* On input the 64-byte cache lines hold these data values: */ \
/*	0	1	256	257	 512	513	768	769 */ \
/*	2	... */ \
/*	...	... */ \
/*	254	... */ \
/*	1024	... */ \
/*	... */ \
/* On output the 64-byte cache lines hold these data values: */ \
/*	0	768	1532	2304 */ \
/*	2 */ \
/*	... */ \
/*	766 */ \
\
/* Do 512 six_reals_last_unfft macros */ \
/*	distance between fft data elements is 512 */ \
/*	do 128 macros each processing 24 data values */ \
\
xsix_reals_unfft_preload; \
for(unsigned int loop=128; loop != 0; loop--){ \
	x3cl_six_reals_last_unfft(rsi, 64, 128*64); \
} \
\
xfft_3_ret; \
\
\
\
/* Perform a 3072-element negacyclic FFT */ \
\

#define xfft3072p(type) \
\
/* Do a multiply with pre-FFTed inputs */ \
\
if(g->ffttype == 4) { \
	xmiddle_4p(g); \
}else{ \
	\
	/* Do FFT levels 1,2,3 */ \
	/* */ \
	/* On input the 64-byte cache lines hold these data values: */ \
	/*	0	768	1	769	1536	2304	1537	2305 */ \
	/*	2 */ \
	/*	... */ \
	/*	766 */ \
	/* On output the 64-byte cache lines hold these data values (in +1 format): */ \
	/*	0	256	512	768 */ \
	/*	2	.. */ \
	/*	... */ \
	/*	254 */ \
	/*	1024	1532	2048	2560 */ \
	/*	1026 */ \
	/*	... */ \
	/*	1530 */ \
	\
	/* Do 512 three_complex_first_fft macros */ \
	/*	distance between fft data elements is 512 */ \
	/*	do 128 macros each processing 24 data values */ \
	\
	rdi = (uintptr_t)g->plus1_premults;/* Address of premultiplier table */ \
	for(unsigned int loop=128; loop != 0; loop--){ \
		s3cl_three_complex_first_fft(rsi, 64, 128*64); \
		rdi += 192;		/* Next premultiplier pointer */ \
	} \
	rsi += -128*64;		/* Restore source pointer */ \
	\
	/* Do FFT levels 4 */ \
	/* */ \
	/* On input the 64-byte cache lines hold these data values: */ \
	/*	0	256	512	768 */ \
	/*	2	.. */ \
	/*	... */ \
	/*	254 */ \
	/*	1024	1532	2048	2560 */ \
	/*	1026 */ \
	/*	... */ \
	/*	1530 */ \
	/* On output the 64-byte cache lines hold these data values: */ \
	/*	0	256	512	768 */ \
	/*	2	... */ \
	/*	... */ \
	/*	254 */ \
	/*	1024	... */ \
	/*	... */ \
	\
	/* Do 384 two_two_complex_fft macros */ \
	/*	distance between fft data elements is 256 */ \
	/*	do 192 macros each processing 16 data values */ \
	\
	rdi = (uintptr_t)g->sincos5;		/* Load sin/cos pointer */ \
	for(unsigned int loop=64; loop != 0; loop--){ \
		x2cl_two_complex_fft_in_place(rsi, 64, 64*64, rdi); \
	} \
	rsi += -64*64+128*64;	/* Next source pointer */ \
	rdi += XMM_SCD;		/* Next sin/cos pointer */ \
	for(unsigned int loop=128; loop != 0; loop--){ \
		x2cl_two_complex_fft(rsi, 64, 128*64, rdi); \
	} \
	rsi += -128*64-128*64;	/* Restore source pointer */ \
	\
	/* Join common code to perform rest of the FFT */ \
	\
	xlast_six_fft_levels_p; \
} \
cmn4096p; \
\
/* Do inverse FFT levels 1,2,3 */ \
/* On input the 64-byte cache lines hold these data values: */ \
/*	0	1	256	257	 512	513	768	769 */ \
/*	2	... */ \
/*	...	... */ \
/*	254	... */ \
/*	1024	... */ \
/*	... */ \
/* On output the 64-byte cache lines hold these data values: */ \
/*	0	768	1532	2304 */ \
/*	2 */ \
/*	... */ \
/*	766 */ \
\
/* Do 512 three_complex_last_unfft macros */ \
/*	distance between fft data elements is 512 */ \
/*	do 128 macros each processing 24 data values */ \
\
rdi = (uintptr_t)g->plus1_premults;/* Address of premultiplier table */ \
for(unsigned int loop=128; loop != 0; loop--){ \
	x3cl_three_complex_last_unfft(rsi, 64, 128*64); \
	rdi += 192;		/* Next premultiplier pointer */ \
} \
\
xfft_3_ret; \
\
\
\
/* Perform a 3584-element FFT. */ \
\

#define xfft3584(type) \
\
/* Do a multiply with pre-FFTed inputs */ \
\
if(g->ffttype == 4) { \
	xmiddle_4(g); \
}else{ \
	\
	/* Do FFT levels 1,2,3 */ \
	/* Values 0-3583 is real data. */ \
	/* */ \
	/* On input the 64-byte cache lines hold these data values: */ \
	/*	0	128	1	129	1792	1920	1793	1921 */ \
	/*	2	... */ \
	/*	... */ \
	/*	126 */ \
	/*	256	512	257	513	2048	2304	2049	2305 */ \
	/*	258 */ \
	/*	... */ \
	/*	510 */ \
	/*	768	1280	769	1281	2560	3072	2561	3073 */ \
	/*	770 */ \
	/*	... */ \
	/*	1278 */ \
	/* On output the 64-byte cache lines hold these data values (in +1 format): */ \
	/*	0	128	256	384 */ \
	/*	2	... */ \
	/*	... */ \
	/*	126 */ \
	/*	512	768	1024	1280 */ \
	/*	514 */ \
	/*	... */ \
	/*	766 */ \
	/*	1536	2048	2560	3072 */ \
	/*	1538 */ \
	/*	... */ \
	/*	2046 */ \
	\
	/* Do 512 seven_reals_first_fft macros */ \
	/*	distance between fft data elements is 512 */ \
	/*	do 64 macros each processing 56 data values */ \
	\
	xseven_reals_fft_preload; \
	for(unsigned int loop=64; loop != 0; loop--){ \
		s7cl_seven_reals_first_fft(rsi, 64, 64*64); \
	} \
	rsi += -64*64;		/* Restore source pointer */ \
	\
	/* Do FFT levels 4,5 */ \
	/* Values 0-511 is real data, 512-3583 is complex data. */ \
	\
	/* Do 64 eight_reals_fft_2 macros */ \
	/*	distance between fft data elements is 128 */ \
	/*	do 64 macros each processing 8 data values */ \
	\
	for(unsigned int loop=64; loop != 0; loop--){ \
		x1cl_half_eight_reals_fft_2(rsi, 64); \
	} \
	\
	/* Do 128 four_complex_fft macros */ \
	/*	distance between fft data elements is 128 */ \
	/*	do 64 macros each processing 16 data values */ \
	\
	rdi = (uintptr_t)g->sincos5;		/* Load sin/cos pointer */ \
	for(unsigned int loop=64; loop != 0; loop--){ \
		x2cl_four_complex_fft(rsi, 64, 64*64); \
	} \
	rsi += -64*64+128*64;	/* Next source pointer */ \
	\
	/* Do 256 four_complex_fft macros */ \
	/*	distance between fft data elements is 128 */ \
	/*	do 64 macros each processing 32 data values */ \
	\
	rdi += XMM_SCD;		/* Next sine/cosine pointer */ \
	for(unsigned int loop=64; loop != 0; loop--){ \
		x4cl_four_complex_fft(rsi, 64, 64*64, 128*64, rdi); \
	} \
	rsi = (uintptr_t)g->DESTARG;/* Restore source pointer */ \
	\
	/* Join common code to perform rest of the FFT */ \
	\
	cmn4096_2; \
} \
/*rcx = 64*64;*/ \
cmn4096(64*64); \
\
/* Finish inverse FFT levels 4,5 */ \
\
/* Do 64 eight_reals_unfft_2 macros */ \
/*	distance between fft data elements is 128 */ \
/*	do 32 macros each processing 16 data values */ \
\
for(unsigned int loop=32; loop != 0; loop--){ \
	x2cl_half_eight_reals_unfft_2(rsi, 64, 32*64); \
} \
rsi += -32*64;		/* Restore source pointer */ \
\
/* Do inverse FFT levels 1,2,3 */ \
/* On input the 64-byte cache lines hold these data values: */ \
/*	0	128	256	384 */ \
/*	2	... */ \
/*	... */ \
/*	126 */ \
/*	512	768	1024	1280 */ \
/*	... */ \
/* On output the 64-byte cache lines hold these data values: */ \
/*	0	128	1792	1920 */ \
/*	2	... */ \
/*	... */ \
/*	126 */ \
/*	256	512	2048	2304 */ \
/*	258 */ \
/*	... */ \
/*	510 */ \
/*	768	1280	2560	3072 */ \
/*	770 */ \
/*	... */ \
/*	1278 */ \
\
/* Do 512 seven_reals_last_unfft macros */ \
/*	distance between fft data elements is 512 */ \
/*	do 64 macros each processing 56 data values */ \
\
xseven_reals_unfft_preload; \
for(unsigned int loop=64; loop != 0; loop--){ \
	x7cl_seven_reals_last_unfft(rsi, 64, 64*64); \
} \
\
xfft_3_ret; \
\
\
\
/* Perform a 4096-element FFT. */ \
\

#define xfft4096(type) \
\
/* Do a multiply with pre-FFTed inputs */ \
\
if(g->ffttype == 4) { \
	xmiddle_4(g); \
}else{ \
	\
	/* Do FFT levels 1,2,3 */ \
	/* Values 0-4095 is real data. */ \
	/* */ \
	/* On input the 64-byte cache lines hold these data values: */ \
	/*	0	1024	1	1025	2048	3072	2049	3073 */ \
	/*	2	... */ \
	/*	... */ \
	/*	1022 */ \
	/* On output the 64-byte cache lines hold these data values: */ \
	/*	0	1	512	513	1024	1025	1536	1537 */ \
	/*	2	... */ \
	/*	... */ \
	/*	510	... */ \
	/*	2048	... */ \
	/*	... */ \
	\
	/* Do 512 eight_reals_first_fft macros */ \
	/*	distance between fft data elements is 512 */ \
	/*	do 256 macros each processing 16 data values */ \
	\
	for(unsigned int loop=256; loop != 0; loop--){ \
		s2cl_eight_reals_first_fft(rsi, 64, 256*64); \
	} \
	rsi += -256*64;		/* Restore source pointer */ \
	\
	/* Do FFT levels 4,5 */ \
	/* Values 0-511 is real data, 512-1023 is semi-real data, 1024-4095 is */ \
	/* complex data. */ \
	/* */ \
	/* On input the 64-byte cache lines hold these data values: */ \
	/*	0	1	512	513	1024	1025	1536	1537 */ \
	/*	2	... */ \
	/*	... */ \
	/*	510	... */ \
	/*	2048	... */ \
	/*	... */ \
	/* On output the 64-byte cache lines hold these data values: */ \
	/*	0	1	128	129	256	257	384	385 */ \
	/*	2	... */ \
	/*	... */ \
	/*	126	... */ \
	/*	512	... */ \
	/*	... */ \
	\
	/* Do 64 eight_reals_fft_2 macros */ \
	/* Do 64 nop_two_two_complex_fft_2 macros */ \
	/* Do 128 four_complex_fft macros */ \
	/*	distance between fft data elements is 128 */ \
	/*	do 64 macros each processing 32 data values */ \
	\
	rdi = (uintptr_t)g->sincos5;		/* Load sin/cos pointer */ \
	for(unsigned int loop=64; loop != 0; loop--){ \
		x4cl_eight_reals_fft_2(rsi, 64, 64*64, 128*64, rdi); \
	} \
	rsi += -64*64+256*64;	/* Next source pointer */ \
	\
	/* Do 256 four_complex_fft macros */ \
	/*	distance between fft data elements is 128 */ \
	/*	do 64 macros each processing 32 data values */ \
	\
	rdi += XMM_SCD;		/* Next sine/cosine pointer */ \
	for(unsigned int loop=64; loop != 0; loop--){ \
		x4cl_four_complex_fft(rsi, 64, 64*64, 128*64, rdi); \
	} \
	rsi = (uintptr_t)g->DESTARG;/* Restore source pointer */ \
	\
	/* Call common code for rest of FFT */ \
	\
	cmn4096_2; \
} \
/*rcx = 128*64;*/ \
cmn4096(128*64); \
\
/* Finish the inverse FFT levels 4,5 */ \
\
/* Do 32 eight_reals_unfft_2 macros */ \
/* Do 32 nop_two_two_complex_unfft_2 macros */ \
/*	distance between fft data elements is 128 */ \
/*	do 32 macros each processing 32 data values */ \
\
for(unsigned int loop=32; loop != 0; loop--){ \
	x4cl_eight_reals_unfft_2(rsi, 64, 32*64, 64*64); \
} \
rsi += -32*64;		/* Restore source pointer */ \
\
/* Do inverse FFT levels 1,2,3 */ \
/* On input the 64-byte cache lines hold these data values: */ \
/*	0	1	256	257	512	513	768	768 */ \
/*	2	... */ \
/*	...	... */ \
/*	254	... */ \
/*	1024	... */ \
/*	... */ \
/* On output the 64-byte cache lines hold these data values: */ \
/*	0	1	1024	1025	2048	2049	3072	3073 */ \
/*	2	... */ \
/*	... */ \
/*	1022 */ \
\
/* Do 512 eight_reals_last_unfft macros */ \
/*	distance between fft data elements is 512 */ \
/*	do 128 macros each processing 32 data values */ \
\
for(unsigned int loop=128; loop != 0; loop--){ \
	x4cl_eight_reals_last_unfft(rsi, 64, 128*64, 256*64); \
} \
\
xfft_3_ret; \


/* Common code for length 2560, 3072, 3584, and 4096 FFTs */ \

#define cmn4096_2 \
/*push_amt = SZPTR;*/ \
\
/* Do FFT levels 6 */ \
/* Values 0-127 is real data, 128-255 is semi-real data, 256-4095 is */ \
/* complex data. */ \
/* */ \
/* On input the 64-byte cache lines hold these data values: */ \
/*	0	1	128	129	256	257	384	385 */ \
/*	2	... */ \
/*	... */ \
/*	126	... */ \
/*	512	... */ \
/*	... */ \
/* On output the 64-byte cache lines hold these data values: */ \
/*	0	1	64	65	128	129	192	193 */ \
/*	2	... */ \
/*	... */ \
/*	62	... */ \
/*	256	... */ \
/*	... */ \
\
/* Do 32 eight_reals_fft_1 macros */ \
/* Do 32 two_two_complex_fft macros */ \
/*	distance between fft data elements is 64 */ \
/*	do 32 macros each processing 16 data values */ \
\
rdi = (uintptr_t)g->sincos4;		/* Load sin/cos pointer */ \
for(unsigned int loop=32; loop != 0; loop--){ \
	x2cl_eight_reals_fft_1(rsi, 64, 32*64); \
} \
\
/* Do 256/320/384/448 two_two_complex_fft macros */ \
/*	distance between fft data elements is 64 */ \
/*	do 128/160/192/224 macros each processing 16 data values */ \
\
rsi += -32*64+64*64;	/* Load source pointer */ \
rdi += XMM_SCD;		/* Next sine/cosine pointer */ \
/* 4/5/6/7 iterations of 32 */ \
for(unsigned int loop = (g->count5 >> 8) & 0xFF; loop != 0; loop--){ \
	for(unsigned int loopB=32; loopB != 0; loopB--){ \
		x2cl_two_complex_fft(rsi, 64, 32*64, rdi); \
	} \
	rsi += -32*64+64*64;	/* Next source pointer */ \
	rdi += 2*XMM_SCD;		/* Next sine/cosine pointer */ \
} \
\
/* Do FFT levels 7,8 */ \
/*	do 72/88/104/120 macros each processing 32 data values */ \
/*	mov	ax, WORD PTR count4+2	;; 1/3 then 8/12 iterations of 8 */ \
/* Do FFT levels 9,10 */ \
/*	do 78/94/110/126 macros each processing 32 data values */ \
/*	mov	ax, WORD PTR count4	;; 7/15 then 32/48 iterations of 2 */ \
/* Do FFT levels 11,12 */ \
/*	do 159/191/223/255 macros each processing 16 data values */ \
/*	mov	eax, count1		;; 31/63 then 128/192 iterations */ \
\
xlast_four_fft_levels; /* jmp */ \
\
/* Do inverse FFT levels 9,10 */ \
/* On input the 64-byte cache lines hold these data values: */ \
/*	0	1	2	3	4	5	6	7 */ \
/*	8	9	10	11	... */ \
/*	... */ \
/* On output the 64-byte cache lines hold these data values: */ \
/*	0	1	8	9	16	17	24	25 */ \
/*	2	... */ \
/*	4	... */ \
/*	6	... */ \
/*	32	... */ \
/*	... */ \

/* Do 2 eight_reals_unfft_2 macros */ \
/* Do 2 nop_two_two_complex_unfft_2 macros */ \
/*	distance between fft data elements is 4 */ \
/*	do 1 macro processing 32 data values */ \

#define cmn4096(RCX1) \
x4cl_eight_reals_unfft_2(rsi, 4*64, 64, 2*64); \
\
/* Do 316/380/444/508 four_complex_unfft macros */ \
/*	distance between fft data elements is 4 */ \
/*	do 79/95/111/127 macros each processing 32 data values */ \
\
rdi = (uintptr_t)g->sincos2;/* Load sin/cos pointer */ \
for(unsigned int loop = g->count2 & 0xFF; loop != 0; loop--){ \
	x4cl_four_complex_unfft(rsi, 4*64, 64, 2*64, rdi); \
	rdi += XMM_SCD;		/* Next sine/cosine pointer */ \
} \
rsi = (uintptr_t)g->DESTARG;/* Restore source pointer */ \
\
/* Do inverse FFT levels 7,8 */ \
/* On input the 64-byte cache lines hold these data values: */ \
/*	0	1	8	9	16	17	24	25 */ \
/*	2	... */ \
/*	4	... */ \
/*	6	... */ \
/*	32	... */ \
/*	... */ \
/* On output the 64-byte cache lines hold these data values: */ \
/*	0	1	32	33	64	65	96	97 */ \
/*	2	... */ \
/*	...	... */ \
/*	30	... */ \
/*	128	... */ \
/*	... */ \
\
/* Do 8 eight_reals_unfft_2 macros */ \
/* Do 8 nop_two_two_complex_unfft_2 macros */ \
/*	distance between fft data elements is 16 */ \
/*	do 4 macros each processing 32 data values */ \
\
for(unsigned int loop=4; loop != 0; loop--){ \
	x4cl_eight_reals_unfft_2(rsi, 64, 4*64, 8*64); \
} \
\
/* Do 304/368/432/496 four_complex_unfft macros */ \
/*	distance between fft data elements is 16 */ \
/*	do 76/92/108/124 macros each processing 32 data values */ \
\
rdi = (uintptr_t)g->sincos3;		/* Load sin/cos pointer */ \
rsi += -4*64+16*64;	/* Load source pointer */ \
for(unsigned int loop = g->count3 & 0xFF; loop != 0; loop--){ \
	for(unsigned int loopB=4; loopB != 0; loopB--){ \
		x4cl_four_complex_unfft(rsi, 64, 4*64, 8*64, rdi); \
	} \
	rsi += -4*64+16*64;	/* Next source pointer */ \
	rdi += XMM_SCD;		/* Next sine/cosine pointer */ \
} \
rsi = (uintptr_t)g->DESTARG;/* Restore source pointer */ \
\
/* Do inverse FFT levels 6 */ \
/* On input the 64-byte cache lines hold these data values: */ \
/*	0	1	32	33	64	65	96	97 */ \
/*	2	... */ \
/*	...	... */ \
/*	30	... */ \
/*	128	... */ \
/*	... */ \
/* On output the 64-byte cache lines hold these data values: */ \
/*	0	1	64	65	128	129	192	193 */ \
/*	2	... */ \
/*	...	... */ \
/*	62	... */ \
/*	256	... */ \
/*	... */ \
\
/* Do 32 eight_reals_unfft_1 macros */ \
/*	distance between fft data elements is 64 */ \
/*	do 16 macros each processing 16 data values */ \
\
for(unsigned int loop=16; loop != 0; loop--){ \
	x2cl_eight_reals_unfft_1(rsi, 64, 16*64); \
} \
\
/* Do 288/352/416/480 two_two_complex_unfft macros */ \
/*	distance between fft data elements is 64 */ \
/*	do 144/176/208/240 macros each processing 16 data values */ \
\
rdi = (uintptr_t)g->sincos4;		/* Load sin/cos pointer */ \
rsi += -16*64+32*64;	/* Load source pointer */ \
/* 1/3 then 8/12 iterations of 16 */ \
for(unsigned int loop = (g->count3 >> 8) & 0xFF; loop != 0; loop--){ \
	for(unsigned int loopB=16; loopB != 0; loopB--){ \
		x2cl_two_complex_unfft(rsi, 64, 16*64); \
	} \
	rsi += -16*64+32*64;	/* Next source pointer */ \
	rdi += XMM_SCD;		/* Next sine/cosine pointer */ \
} \
rsi = (uintptr_t)g->DESTARG;/* Restore source pointer */ \
\
/* Do inverse FFT levels 4,5 */ \
/* On input the 64-byte cache lines hold these data values: */ \
/*	0	1	64	65	128	129	192	193 */ \
/*	2	... */ \
/*	...	... */ \
/*	62	... */ \
/*	256	... */ \
/*	... */ \
/* On output the 64-byte cache lines hold these data values: */ \
/*	0	1	256	257	512	513	768	769 */ \
/*	2	... */ \
/*	...	... */ \
/*	254	... */ \
/*	1024	... */ \
/*	... */ \
\
/* Do 256/384 four_complex_unfft macros */ \
/*	distance between fft data elements is 128 */ \
/*	do 64/96 macros each processing 32 data values */ \
\
rsi += RCX1;		/* Next source pointer */ \
rdi = (uintptr_t)g->sincos5;		/* Load sin/cos pointer */ \
/* 2/3 iterations of 32 */ \
for(unsigned int loop = g->count5 & 0xFF; loop != 0; loop--){ \
	for(unsigned int loopB=32; loopB != 0; loopB--){ \
		x4cl_four_complex_unfft(rsi, 64, 32*64, 64*64, rdi); \
	} \
	rsi += -32*64+128*64;	/* Next source pointer */ \
	rdi += XMM_SCD;		/* Next sine/cosine pointer */ \
} \
rsi = (uintptr_t)g->DESTARG;/* Next source pointer */ \
\
/* This ends the common inverse FFT code */ \
\
/*ret; push_amt = 0; */ \



/* Perform a 4096-element negacyclic FFT */ \


#define xfft4096p(type) \
\
/* Do a multiply with pre-FFTed inputs */ \
\
if(g->ffttype == 4) { \
	xmiddle_4p(g); \
}else{ \
	\
	/* Do FFT levels 1,2,3 */ \
	/* */ \
	/* On input the 64-byte cache lines hold these data values: */ \
	/*	0	1024	1	1025	2048	3072	2049	3073 */ \
	/*	2	... */ \
	/*	... */ \
	/*	1022 */ \
	/* On output the 64-byte cache lines hold these data values: */ \
	/*	0	1	512	513	1024	1025	1536	1537 */ \
	/*	2	... */ \
	/*	... */ \
	/*	510	... */ \
	/*	2048	... */ \
	/*	... */ \
	\
	/* Do 512 four_complex_first_fft macros */ \
	/*	distance between fft data elements is 512 */ \
	/*	do 256 macros each processing 16 data values */ \
	\
	rdi = (uintptr_t)g->plus1_premults;/* Address of premultiplier table */ \
	for(unsigned int loop=256; loop != 0; loop--){ \
		s2cl_four_complex_first_fft(rsi, 64, 256*64); \
		rdi += 128;		/* Next premultiplier pointer */ \
	} \
	rsi += -256*64;		/* Restore source pointer */ \
	\
	/* Do FFT level 4 */ \
	/* */ \
	/* On input the 64-byte cache lines hold these data values: */ \
	/*	0	1	512	513	1024	1025	1536	1537 */ \
	/*	2	... */ \
	/*	... */ \
	/*	510	... */ \
	/*	2048	... */ \
	/*	... */ \
	/* On output the 64-byte cache lines hold these data values: */ \
	/*	0	1	256	257	512	513	768	769 */ \
	/*	2	3	... */ \
	/*	... */ \
	/*	254 */ \
	/*	1024	... */ \
	/*	... */ \
	\
	/* Do 512 two_two_complex_fft macros */ \
	/*	distance between fft data elements is 256 */ \
	/*	do 256 macros each processing 16 data values */ \
	\
	rdi = (uintptr_t)g->sincos5;		/* Load sin/cos pointer */ \
	for(unsigned int loop=2; loop != 0; loop--){ \
		for(unsigned int loopB=128; loopB != 0; loopB--){ \
			x2cl_two_complex_fft(rsi, 64, 128*64, rdi); \
		} \
		rsi += -128*64+256*64;	/* Restore source pointer */ \
		rdi += 2*XMM_SCD;		/* Next sine/cosine pointer */ \
	} \
	rsi += -2*256*64;		/* Restore source pointer */ \
	\
	/* Do FFT levels 5,6 */ \
	/*	do 96/128 macros each processing 32 data values */ \
	/*	mov	ax, WORD PTR count4	;; 3/4 iterations of 32 */ \
	/* Do FFT levels 7,8 */ \
	/*	do 96/128 macros each processing 32 data values */ \
	/*	mov	ax, WORD PTR count4+2	;; 12/16 iterations of 8 */ \
	/* Do FFT levels 9,10 */ \
	/*	do 96/128 macros each processing 32 data values */ \
	/*	mov	ax, WORD PTR count4	;; 48/64 iterations of 2 */ \
	/* Do FFT levels 11,12 */ \
	/*	do 192/256 macros each processing 16 data values */ \
	/*	mov	eax, count1		;; 192/256 iterations */ \
	\
	xlast_six_fft_levels_p; \
} \
cmn4096p; \
\
/* Do inverse FFT levels 1,2,3 */ \
/* On input the 64-byte cache lines hold these data values: */ \
/*	0	1	256	257	512	513	768	768 */ \
/*	2	... */ \
/*	...	... */ \
/*	254	... */ \
/*	1024	... */ \
/*	... */ \
/* On output the 64-byte cache lines hold these data values: */ \
/*	0	1	1024	1025	2048	2049	3072	3073 */ \
/*	2	... */ \
/*	... */ \
/*	1022 */ \
\
/* Do 512 four_complex_last_unfft macros */ \
/*	distance between fft data elements is 512 */ \
/*	do 128 macros each processing 32 data values */ \
\
rdi = (uintptr_t)g->plus1_premults;/* Address of premultiplier table */ \
for(unsigned int loop=128; loop != 0; loop--){ \
	x4cl_four_complex_last_unfft(rsi, 64, 128*64, 256*64, 128*128); \
	rdi += 128;		/* Next premultiplier pointer */ \
} \
\
xfft_3_ret; \
\
\
/* Common inverse FFT code for lengths 3072 and 4096 */ \
\
/* Do inverse FFT levels 9,10 */ \
/* On input the 64-byte cache lines hold these data values: */ \
/*	0	1	2	3	4	5	6	7 */ \
/*	8	9	10	11	... */ \
/*	... */ \
/* On output the 64-byte cache lines hold these data values: */ \
/*	0	1	8	9	16	17	24	25 */ \
/*	2	... */ \
/*	4	... */ \
/*	6	... */ \
/*	32	... */ \
/*	... */ \
\
/* Do 384/512 four_complex_unfft macros */ \
/*	distance between fft data elements is 4 */ \
/*	do 96/128 macros each processing 32 data values */ \

#define cmn4096p \
/* push_amt = SZPTR; */ \
rdi = (uintptr_t)g->sincos2;/* Load sin/cos pointer */ \
/* 96/128 iterations */ \
for(unsigned int loop = g->count2 & 0xFF; loop != 0; loop--){ \
	x4cl_four_complex_unfft(rsi, 4*64, 64, 2*64, rdi); \
	rdi += XMM_SCD;		/* Next sine/cosine pointer */ \
} \
rsi = (uintptr_t)g->DESTARG;/* Restore source pointer */ \
\
/* Do inverse FFT levels 7,8 */ \
/* On input the 64-byte cache lines hold these data values: */ \
/*	0	1	8	9	16	17	24	25 */ \
/*	2	... */ \
/*	4	... */ \
/*	6	... */ \
/*	32	... */ \
/*	... */ \
/* On output the 64-byte cache lines hold these data values: */ \
/*	0	1	32	33	64	65	96	97 */ \
/*	2	... */ \
/*	...	... */ \
/*	30	... */ \
/*	128	... */ \
/*	... */ \
\
/* Do 384/512 four_complex_unfft macros */ \
/*	distance between fft data elements is 16 */ \
/*	do 96/128 macros each processing 32 data values */ \
\
rdi = (uintptr_t)g->sincos3;		/* Load sin/cos pointer */ \
/* 24/32 iterations of 4 */ \
for(unsigned int loop = g->count3 & 0xFF; loop != 0; loop--){ \
	for(unsigned int loopB=4; loopB != 0; loopB--){ \
		x4cl_four_complex_unfft(rsi, 64, 4*64, 8*64, rdi); \
	} \
	rsi += -4*64+16*64;	/* Next source pointer */ \
	rdi += XMM_SCD;		/* Next sine/cosine pointer */ \
} \
rsi = (uintptr_t)g->DESTARG;/* Restore source pointer */ \
\
/* Do inverse FFT levels 5,6 */ \
/* On input the 64-byte cache lines hold these data values: */ \
/*	0	1	32	33	64	65	96	97 */ \
/*	2	... */ \
/*	...	... */ \
/*	30	... */ \
/*	128	... */ \
/*	... */ \
/* On output the 64-byte cache lines hold these data values: */ \
/*	0	1	128	129	256	257	384	385 */ \
/*	2	... */ \
/*	...	... */ \
/*	126	... */ \
/*	512	... */ \
/*	... */ \
\
/* Do 384/512 four_complex_unfft macros */ \
/*	distance between fft data elements is 16 */ \
/*	do 96/128 macros each processing 32 data values */ \
\
rdi = (uintptr_t)g->sincos4;		/* Load sin/cos pointer */ \
/* 6/8 iterations of 16 */ \
for(unsigned int loop = (g->count3 >> 8) & 0xFF; loop != 0; loop--){ \
	for(unsigned int loopA=16; loopA != 0; loopA--){ \
		x4cl_four_complex_unfft(rsi, 64, 16*64, 32*64, rdi); \
	} \
	rsi += -16*64+64*64;	/* Next source pointer */ \
	rdi += XMM_SCD;		/* Next sine/cosine pointer */ \
} \
rsi = (uintptr_t)g->DESTARG;/* Restore source pointer */ \
\
/* Do inverse FFT level 4 */ \
/* On input the 64-byte cache lines hold these data values: */ \
/*	0	1	128	129	256	257	384	385 */ \
/*	2	... */ \
/*	...	... */ \
/*	126	... */ \
/*	512	... */ \
/*	... */ \
/* On output the 64-byte cache lines hold these data values: */ \
/*	0	1	256	257	512	513	768	769 */ \
/*	2	... */ \
/*	...	... */ \
/*	254	... */ \
/*	1024	... */ \
/*	... */ \
\
/* Do 384/512 two_two_complex_unfft macros */ \
/*	distance between fft data elements is 256 */ \
/*	do 192/256 macros each processing 16 data values */ \
\
rdi = (uintptr_t)g->sincos5;		/* Load sin/cos pointer */ \
/* 3/4 iterations of 64 */ \
for(unsigned int loop = g->count5 & 0xFF; loop != 0; loop--){ \
	for(unsigned int loopB=64; loopB != 0; loopB--){ \
		x2cl_two_complex_unfft(rsi, 64, 64*64); \
	} \
	rsi += -64*64+128*64;	/* Next source pointer */ \
	rdi += XMM_SCD;		/* Next sine/cosine pointer */ \
} \
rsi = (uintptr_t)g->DESTARG;/* Restore source pointer */ \
\
/* This ends the common inverse FFT code */ \
\
/* ret; push_amt = 0; */ \



/* Perform a 5120-element FFT. */


#define xfft5120(type) \
\
/* Do a multiply with pre-FFTed inputs */ \
\
if(g->ffttype == 4) { \
	xmiddle_4(g); \
}else{ \
	\
	/* Do FFT levels 1,2,3 */ \
	/* Values 0-5119 is real data. */ \
	/* */ \
	/* On input the 64-byte cache lines hold these data values: */ \
	/*	0	256	1	257	2560	2816	2561	2817 */ \
	/*	2	... */ \
	/*	... */ \
	/*	254 */ \
	/*	512	1536	513	1537	3072	4096	3073	4097 */ \
	/*	514 */ \
	/*	... */ \
	/*	1534 */ \
	/* On output the 64-byte cache lines hold these data values (in +1 format): */ \
	/*	0	256	512	768 */ \
	/*	2	... */ \
	/*	... */ \
	/*	254 */ \
	/*	1024	2048	3072	4096 */ \
	/*	1026 */ \
	/*	... */ \
	/*	2046 */ \
	\
	/* Do 1024 five_reals_first_fft macros */ \
	/*	distance between fft data elements is 1024 */ \
	/*	do 128 macros each processing 40 data values */ \
	\
	xfive_reals_fft_preload;		/* Preload const1s for five real macros */ \
	for(unsigned int loop=128; loop != 0; loop--){ \
		s5cl_five_reals_first_fft(rsi, 64, 128*64); \
	} \
	rsi += -128*64;		/* Restore source pointer */ \
	\
	/* Do FFT levels 4,5 */ \
	/* Values 0-1023 is real data, 1024-5119 is complex data. */ \
	\
	/* Do 128 eight_reals_fft_2 macros */ \
	/*	distance between fft data elements is 256 */ \
	/*	do 128 macros each processing 8 data values */ \
	\
	for(unsigned int loop=128; loop != 0; loop--){ \
		x1cl_half_eight_reals_fft_2(rsi, 64); \
	} \
	\
	/* Do 256 four_complex_fft macros */ \
	/*	distance between fft data elements is 256 */ \
	/*	do 128 macros each processing 32 data values */ \
	\
	rdi = (uintptr_t)g->sincos5;		/* Load sin/cos pointer */ \
	for(unsigned int loop=128; loop != 0; loop--){ \
		x4cl_four_complex_fft(rsi, 64, 128*64, 256*64, rdi); \
	} \
	rsi = (uintptr_t)g->DESTARG;/* Restore source pointer */ \
	\
	/* Join common code to perform rest of the FFT */ \
	\
	cmn8192_2; \
} \
/*rcx = 128*64;*/ \
cmn8192(128*64); \
\
/* Finish the inverse FFT levels 4,5 */ \
\
/* Do 128 eight_reals_unfft_2 macros */ \
/*	distance between fft data elements is 256 */ \
/*	do 64 macros each processing 16 data values */ \
\
for(unsigned int loop=64; loop != 0; loop--){ \
	x2cl_half_eight_reals_unfft_2(rsi, 64, 64*64); \
} \
rsi += -64*64;		/* Restore source pointer */ \
\
/* Do inverse FFT levels 1,2,3 */ \
/* On output the 64-byte cache lines hold these data values: */ \
/*	0	256	2560	2816 */ \
/*	2	... */ \
/*	... */ \
/*	254 */ \
/*	512	1536	3072	4096 */ \
/*	514 */ \
/*	... */ \
/*	1534 */ \
\
/* Do 1024 five_reals_last_unfft macros */ \
/*	distance between fft data elements is 1024 */ \
/*	do 128 macros each processing 40 data values */ \
\
xfive_reals_unfft_preload;	/* Preload const1s for five real macros */ \
for(unsigned int loop=128; loop != 0; loop--){ \
	x5cl_five_reals_last_unfft(rsi, 64, 128*64); \
} \
\
xfft_3_ret; \
\
\
\
/* Perform a 6144-element FFT. */ \
\

#define xfft6144(type) \
\
/* Do a multiply with pre-FFTed inputs */ \
\
if(g->ffttype == 4) { \
	xmiddle_4(g); \
}else{ \
	\
	/* Do FFT levels 1,2,3 */ \
	/* Values 0-6143 is real data. */ \
	/* */ \
	/* On input the 64-byte cache lines hold these data values: */ \
	/*	0	1536	1	1537	3072	4608	3073	4609 */ \
	/*	2	... */ \
	/*	... */ \
	/*	1534 */ \
	/* On output the 64-byte cache lines hold these data values (in +1 format): */ \
	/*	0	512	1024	1536 */ \
	/*	2	... */ \
	/*	... */ \
	/*	510	... */ \
	/*	2048	3072	4096	5120 */ \
	/*	... */ \
	/*	3070 */ \
	\
	/* Do 1024 six_reals_first_fft macros */ \
	/*	distance between fft data elements is 1024 */ \
	/*	do 256 macros each processing 24 data values */ \
	\
	xsix_reals_fft_preload; \
	for(unsigned int loop=256; loop != 0; loop--){ \
		s3cl_six_reals_first_fft(rsi, 64, 256*64); \
	} \
	rsi += -256*64;		/* Restore source pointer */ \
	\
	/* Do FFT levels 4,5 */ \
	/* Values 0-1023 is real data, 1024-2047 is semi-real data, 2048-6143 is */ \
	/* complex data. */ \
	/* */ \
	/* On input the 64-byte cache lines hold these data values: */ \
	/*	0	512	1024	1536 */ \
	/*	2	... */ \
	/*	... */ \
	/*	510	... */ \
	/*	2048	3072	4096	5120 */ \
	/*	... */ \
	/*	3070 */ \
	/* On output the 64-byte cache lines hold these data values: */ \
	/*	0	1	256	257	512	513	768	769 */ \
	/*	2	... */ \
	/*	... */ \
	/*	254	... */ \
	/*	1024	... */ \
	/*	... */ \
	\
	/* Do 128 eight_reals_fft_2 macros */ \
	/* Do 128 nop_two_two_complex_fft_2 macros */ \
	/* Do 256 four_complex_fft macros */ \
	/*	distance between fft data elements is 256 */ \
	/*	do 128 macros each processing 16 data values */ \
	\
	for(unsigned int loop=128; loop != 0; loop--){ \
		x2cl_eight_reals_fft_2(rsi, 64, 128*64); \
	} \
	rsi += -128*64+256*64;	/* Next source pointer */ \
	\
	/* Do 256 four_complex_fft macros */ \
	/*	distance between fft data elements is 256 */ \
	/*	do 128 macros each processing 32 data values */ \
	\
	rdi = (uintptr_t)g->sincos5;		/* Load sin/cos pointer */ \
	for(unsigned int loop=128; loop != 0; loop--){ \
		x4cl_four_complex_fft(rsi, 64, 128*64, 256*64, rdi); \
	} \
	rsi = (uintptr_t)g->DESTARG;/* Restore source pointer */ \
	\
	/* Join common code to perform rest of the FFT */ \
	\
	cmn8192_2; \
} \
/*rcx = 256*64;*/ \
cmn8192(256*64); \
\
/* Finish inverse FFT levels 4,5 */ \
\
/* Do 128 eight_reals_unfft_2 macros */ \
/* Do 128 nop_two_two_complex_unfft_2 macros */ \
/*	distance between fft data elements is 256 */ \
/*	do 64 macros each processing 32 data values */ \
\
for(unsigned int loop=64; loop != 0; loop--){ \
	x4cl_eight_reals_unfft_2(rsi, 64, 64*64, 128*64); \
} \
rsi += -64*64;		/* Restore source pointer */ \
\
/* Do inverse FFT levels 1,2,3 */ \
/* On input the 64-byte cache lines hold these data values: */ \
/*	0	1	512	513	1024	1025	1536	1537 */ \
/*	2	... */ \
/*	...	... */ \
/*	510	... */ \
/*	2048	... */ \
/*	... */ \
/* On output the 64-byte cache lines hold these data values: */ \
/*	0	1532	3064	4608 */ \
/*	2 */ \
/*	... */ \
/*	1530 */ \
\
/* Do 1024 six_reals_last_unfft macros */ \
/*	distance between fft data elements is 1024 */ \
/*	do 256 macros each processing 24 data values */ \
\
xsix_reals_unfft_preload; \
for(unsigned int loop=256; loop != 0; loop--){ \
	x3cl_six_reals_last_unfft(rsi, 64, 256*64); \
} \
\
xfft_3_ret; \
\
\
\
/* Perform a 6144-element negacyclic FFT */ \
\

#define xfft6144p(type) \
\
/* Do a multiply with pre-FFTed inputs */ \
\
if(g->ffttype == 4) { \
	xmiddle_4p(g); \
}else{ \
	\
	/* Do FFT levels 1,2,3 */ \
	/* */ \
	/* On input the 64-byte cache lines hold these data values: */ \
	/*	0	1536	1	1537	3072	4608	3073	4609 */ \
	/*	2	... */ \
	/*	... */ \
	/*	1534 */ \
	/* On output the 64-byte cache lines hold these data values (in +1 format): */ \
	/*	0	512	1024	1536 */ \
	/*	2	... */ \
	/*	... */ \
	/*	510	... */ \
	/*	2048	3072	4096	5120 */ \
	/*	... */ \
	/*	3070 */ \
	\
	/* Do 1024 three_complex_first_fft macros */ \
	/*	distance between fft data elements is 1024 */ \
	/*	do 256 macros each processing 24 data values */ \
	\
	rdi = (uintptr_t)g->plus1_premults;/* Address of premultiplier table */ \
	for(unsigned int loop=256; loop != 0; loop--){ \
		s3cl_three_complex_first_fft(rsi, 64, 256*64); \
		rdi += 192;		/* Next premultiplier pointer */ \
	} \
	rsi += -256*64;		/* Restore source pointer */ \
	\
	/* Do FFT levels 4,5 */ \
	/* */ \
	/* On input the 64-byte cache lines hold these data values: */ \
	/*	0	512	1024	1536 */ \
	/*	2	... */ \
	/*	... */ \
	/*	510	... */ \
	/*	2048	3072	4096	5120 */ \
	/*	... */ \
	/*	3070 */ \
	/* On output the 64-byte cache lines hold these data values: */ \
	/*	0	1	256	257	512	513	768	769 */ \
	/*	2	... */ \
	/*	... */ \
	/*	254	... */ \
	/*	1024	... */ \
	/*	... */ \
	\
	/* Do 256 four_complex_fft macros */ \
	/*	distance between fft data elements is 256 */ \
	/*	do 128 macros each processing 32 data values */ \
	\
	rdi = (uintptr_t)g->sincos5;		/* Load sin/cos pointer */ \
	for(unsigned int loop=128; loop != 0; loop--){ \
		x2cl_four_complex_fft(rsi, 64, 128*64); \
	} \
	rsi += -128*64+256*64;	/* Next source pointer */ \
	rdi += XMM_SCD;		/* Next sin/cos pointer */ \
	for(unsigned int loop=128; loop != 0; loop--){ \
		x4cl_four_complex_fft(rsi, 64, 128*64, 256*64, rdi); \
	} \
	rsi += -128*64-256*64;	/* Restore source pointer */ \
	\
	/* Join common code to perform rest of the FFT */ \
	\
	xlast_six_fft_levels_p; \
} \
cmn8192p; \
\
/* Do inverse FFT levels 1,2,3 */ \
/* On input the 64-byte cache lines hold these data values: */ \
/*	0	1	512	513	1024	1025	1536	1537 */ \
/*	2	... */ \
/*	...	... */ \
/*	510	... */ \
/*	2048	... */ \
/*	... */ \
/* On output the 64-byte cache lines hold these data values: */ \
/*	0	1532	3064	4608 */ \
/*	2 */ \
/*	... */ \
/*	1530 */ \
\
/* Do 1024 three_complex_last_unfft macros */ \
/*	distance between fft data elements is 1024 */ \
/*	do 256 macros each processing 24 data values */ \
\
rdi = (uintptr_t)g->plus1_premults;/* Address of premultiplier table */ \
for(unsigned int loop=256; loop != 0; loop--){ \
	x3cl_three_complex_last_unfft(rsi, 64, 256*64); \
	rdi += 192;		/* Next premultiplier pointer */ \
} \
\
xfft_3_ret; \
\
\
\
/* Perform a 7168-element FFT. */ \
\

#define xfft7168(type) \
\
/* Do a multiply with pre-FFTed inputs */ \
\
if(g->ffttype == 4) { \
	xmiddle_4(g); \
}else{ \
	\
	/* Do FFT levels 1,2,3 */ \
	/* Values 0-7167 is real data. */ \
	/* */ \
	/* On input the 64-byte cache lines hold these data values: */ \
	/*	0	256	1	257	3584	3840	3585	3841 */ \
	/*	2	... */ \
	/*	... */ \
	/*	254 */ \
	/*	512	1024	513	1025	4096	4608	4097	4609 */ \
	/*	514 */ \
	/*	... */ \
	/*	1022 */ \
	/*	1536	2560	1537	2561	5120	6144	5121	6145 */ \
	/*	1538 */ \
	/*	... */ \
	/*	2558 */ \
	/* On output the 64-byte cache lines hold these data values (in +1 format): */ \
	/*	0	256	512	768 */ \
	/*	2	... */ \
	/*	... */ \
	/*	254 */ \
	/*	1024	1536	2048	2560 */ \
	/*	1026 */ \
	/*	... */ \
	/*	1534 */ \
	/*	3072	4096	5120	6144 */ \
	/*	3074 */ \
	/*	... */ \
	/*	4094 */ \
	\
	/* Do 1024 seven_reals_first_fft macros */ \
	/*	distance between fft data elements is 1024 */ \
	/*	do 128 macros each processing 56 data values */ \
	\
	xseven_reals_fft_preload; \
	for(unsigned int loop=128; loop != 0; loop--){ \
		s7cl_seven_reals_first_fft(rsi, 64, 128*64); \
	} \
	rsi += -128*64;		/* Restore source pointer */ \
	\
	/* Do FFT levels 4,5 */ \
	/* Values 0-1023 is real data, 1024-7167 is complex data. */ \
	\
	/* Do 128 eight_reals_fft_2 macros */ \
	/*	distance between fft data elements is 256 */ \
	/*	do 128 macros each processing 8 data values */ \
	\
	for(unsigned int loop=128; loop != 0; loop--){ \
		x1cl_half_eight_reals_fft_2(rsi, 64); \
	} \
	\
	/* Do 256 four_complex_fft macros */ \
	/*	distance between fft data elements is 256 */ \
	/*	do 128 macros each processing 16 data values */ \
	\
	rdi = (uintptr_t)g->sincos5;		/* Load sin/cos pointer */ \
	for(unsigned int loop=128; loop != 0; loop--){ \
		x2cl_four_complex_fft(rsi, 64, 128*64); \
	} \
	rsi += -128*64+256*64;	/* Next source pointer */ \
	\
	/* Do 512 four_complex_fft macros */ \
	/*	distance between fft data elements is 256 */ \
	/*	do 128 macros each processing 32 data values */ \
	\
	rdi += XMM_SCD;		/* Next sine/cosine pointer */ \
	for(unsigned int loop=128; loop != 0; loop--){ \
		x4cl_four_complex_fft(rsi, 64, 128*64, 256*64, rdi); \
	} \
	rsi = (uintptr_t)g->DESTARG;/* Restore source pointer */ \
	\
	/* Join common code to perform rest of the FFT */ \
	\
	cmn8192_2; \
} \
/*rcx = 128*64;*/ \
cmn8192(128*64); \
\
/* Finish inverse FFT levels 4,5 */ \
\
/* Do 128 eight_reals_unfft_2 macros */ \
/*	distance between fft data elements is 256 */ \
/*	do 64 macros each processing 16 data values */ \
\
for(unsigned int loop=64; loop != 0; loop--){ \
	x2cl_half_eight_reals_unfft_2(rsi, 64, 64*64); \
} \
rsi += -64*64;		/* Restore source pointer */ \
\
/* Do inverse FFT levels 1,2,3 */ \
/* On input the 64-byte cache lines hold these data values: */ \
/*	0	256	512	768 */ \
/*	2	... */ \
/*	... */ \
/*	254 */ \
/*	1024	1536	2048	2560 */ \
/*	... */ \
/* On output the 64-byte cache lines hold these data values: */ \
/*	0	256	3584	3840 */ \
/*	2	... */ \
/*	... */ \
/*	254 */ \
/*	512	1024	4096	4608 */ \
/*	514 */ \
/*	... */ \
/*	1022 */ \
/*	1536	2560	5120	6144 */ \
/*	1538 */ \
/*	... */ \
/*	2558 */ \
\
/* Do 1024 seven_reals_last_unfft macros */ \
/*	distance between fft data elements is 1024 */ \
/*	do 128 macros each processing 56 data values */ \
\
xseven_reals_unfft_preload; \
for(unsigned int loop=128; loop != 0; loop--){ \
	x7cl_seven_reals_last_unfft(rsi, 64, 128*64); \
} \
\
xfft_3_ret; \
\
\
\
/* Perform a 8192-element FFT. */ \
\

#define xfft8192(type) \
\
/* Do a multiply with pre-FFTed inputs */ \
\
if(g->ffttype == 4) { \
	xmiddle_4(g); \
}else{ \
	\
	/* Do FFT levels 1,2,3 */ \
	/* Values 0-8191 is real data. */ \
	/* */ \
	/* On input the 64-byte cache lines hold these data values: */ \
	/*	0	2048	1	2049	4096	6144	4097	6145 */ \
	/*	2	... */ \
	/*	... */ \
	/*	2046 */ \
	/* On output the 64-byte cache lines hold these data values: */ \
	/*	0	1	1024	1025	2048	2049	3072	3073 */ \
	/*	2	... */ \
	/*	... */ \
	/*	1022	... */ \
	/*	4096	... */ \
	/*	... */ \
	\
	/* Do 1024 eight_reals_first_fft macros */ \
	/*	distance between fft data elements is 1024 */ \
	/*	do 512 macros each processing 16 data values */ \
	\
	for(unsigned int loop=512; loop != 0; loop--){ \
		s2cl_eight_reals_first_fft(rsi, 64, 512*64); \
	} \
	rsi += -512*64;		/* Restore source pointer */ \
	\
	/* Do FFT levels 4,5 */ \
	/* Values 0-1023 is real data, 1024-2047 is semi-real data, 2048-8191 is */ \
	/* complex data. */ \
	/* */ \
	/* On input the 64-byte cache lines hold these data values: */ \
	/*	0	1	1024	1025	2048	2049	3072	3073 */ \
	/*	2	... */ \
	/*	... */ \
	/*	1022	... */ \
	/*	4096	... */ \
	/*	... */ \
	/* On output the 64-byte cache lines hold these data values: */ \
	/*	0	1	256	257	512	513	768	769 */ \
	/*	2	... */ \
	/*	... */ \
	/*	254	... */ \
	/*	1024	... */ \
	/*	... */ \
	\
	/* Do 128 eight_reals_fft_2 macros */ \
	/* Do 128 nop_two_two_complex_fft_2 macros */ \
	/* Do 256 four_complex_fft macros */ \
	/*	distance between fft data elements is 256 */ \
	/*	do 128 macros each processing 32 data values */ \
	\
	rdi = (uintptr_t)g->sincos5;		/* Load sin/cos pointer */ \
	for(unsigned int loop=128; loop != 0; loop--){ \
		x4cl_eight_reals_fft_2(rsi, 64, 128*64, 256*64, rdi); \
	} \
	rsi += -128*64+512*64;	/* Next source pointer */ \
	\
	/* Do 512 four_complex_fft macros */ \
	/*	distance between fft data elements is 256 */ \
	/*	do 128 macros each processing 32 data values */ \
	\
	rdi += XMM_SCD;		/* Next sine/cosine pointer */ \
	for(unsigned int loop=128; loop != 0; loop--){ \
		x4cl_four_complex_fft(rsi, 64, 128*64, 256*64, rdi); \
	} \
	rsi = (uintptr_t)g->DESTARG;/* Restore source pointer */ \
	\
	/* Call common code for rest of FFT */ \
	\
	cmn8192_2; \
} \
/*rcx = 256*64;*/ \
cmn8192(256*64); \
\
/* Finish the inverse FFT levels 4,5 */ \
\
/* Do 128 eight_reals_unfft_2 macros */ \
/* Do 128 nop_two_two_complex_unfft_2 macros */ \
/*	distance between fft data elements is 256 */ \
/*	do 64 macros each processing 32 data values */ \
\
for(unsigned int loop=64; loop != 0; loop--){ \
	x4cl_eight_reals_unfft_2(rsi, 64, 64*64, 128*64); \
} \
rsi += -64*64;		/* Restore source pointer */ \
\
/* Do inverse FFT levels 1,2,3 */ \
/* On input the 64-byte cache lines hold these data values: */ \
/*	0	1	512	513	1024	1025	1536	1537 */ \
/*	2	... */ \
/*	...	... */ \
/*	510	... */ \
/*	2048	... */ \
/*	... */ \
/* On output the 64-byte cache lines hold these data values: */ \
/*	0	1	2048	2049	4096	4097	6144	6145 */ \
/*	2	... */ \
/*	... */ \
/*	2046 */ \
\
/* Do 1024 eight_reals_last_unfft macros */ \
/*	distance between fft data elements is 1024 */ \
/*	do 256 macros each processing 32 data values */ \
\
for(unsigned int loop=256; loop != 0; loop--){ \
	x4cl_eight_reals_last_unfft(rsi, 64, 256*64, 512*64); \
} \
\
xfft_3_ret; \
\
\
/* Common code for length 5120, 6144, 7168, and 8192 FFTs */ \

/*push_amt = SZPTR;*/
#define cmn8192_2 \
\
/* Do FFT levels 6,7 */ \
/* Values 0-255 is real data, 256-511 is semi-real data, 512-8191 is */ \
/* complex data. */ \
/* */ \
/* On input the 64-byte cache lines hold these data values: */ \
/*	0	1	256	257	512	513	768	769 */ \
/*	2	... */ \
/*	... */ \
/*	254	... */ \
/*	1024	... */ \
/*	... */ \
/* On output the 64-byte cache lines hold these data values: */ \
/*	0	1	64	65	128	129	192	193 */ \
/*	2	... */ \
/*	... */ \
/*	62	... */ \
/*	256	... */ \
/*	... */ \
\
/* Do 32 eight_reals_fft_2 macros */ \
/* Do 32 nop_two_two_complex_fft_2 macros */ \
/* Do 64 four_complex_fft macros */ \
/*	distance between fft data elements is 64 */ \
/*	do 32 macros each processing 32 data values */ \
\
rdi = (uintptr_t)g->sincos4;		/* Load sin/cos pointer */ \
for(unsigned int loop=32; loop != 0; loop--){ \
	x4cl_eight_reals_fft_2(rsi, 64, 32*64, 64*64, rdi); \
} \
\
/* Do 512/640/768/896 four_complex_fft macros */ \
/*	distance between fft data elements is 64 */ \
/*	do 128/160/192/224 macros each processing 32 data values */ \
\
rsi += -32*64+128*64;	/* Load source pointer */ \
rdi += XMM_SCD;		/* Next sine/cosine pointer */ \
/* 4/5/6/7 iterations of 32 */ \
for(unsigned int loop = (g->count5 >> 8) & 0xFF; loop != 0; loop--){ \
	for(unsigned int loopB=32; loopB != 0; loopB--){ \
	x4cl_four_complex_fft(rsi, 64, 32*64, 64*64, rdi); \
	} \
	rsi += -32*64+128*64;	/* Next source pointer */ \
	rdi += 2*XMM_SCD;		/* Next sine/cosine pointer */ \
} \
\
/* Do FFT levels 8,9 */ \
/* Values 0-63 is real data, 64-127 is semi-real data, 128-8191 is */ \
/* complex data. */ \
/* */ \
/* On input the 64-byte cache lines hold these data values: */ \
/*	0	1	64	65	128	129	192	193 */ \
/*	2	... */ \
/*	... */ \
/*	62	... */ \
/*	256	... */ \
/*	... */ \
/* On output the 64-byte cache lines hold these data values: */ \
/*	0	1	16	17	32	33	48	49 */ \
/*	2	... */ \
/*	... */ \
/*	14	... */ \
/*	64	... */ \
/*	... */ \
\
/* Do 8 eight_reals_fft_2 macros */ \
/* Do 8 nop_two_two_complex_fft_2 macros */ \
/* Do 16 four_complex_fft macros */ \
/*	distance between fft data elements is 16 */ \
/*	do 8 macros each processing 32 data values */ \
xlast_four_fft_levels

#define xlast_four_fft_levels \
rsi = (uintptr_t)g->DESTARG;/* Next source pointer */ \
rdi = (uintptr_t)g->sincos3;		/* Load sin/cos pointer */ \
for(unsigned int loop=8; loop != 0; loop--){ \
	x4cl_eight_reals_fft_2(rsi, 64, 8*64, 16*64, rdi); \
} \
\
/* Do 608/736/864/992 four_complex_fft macros */ \
/*	distance between fft data elements is 16 */ \
/*	do 152/184/216/248 macros each processing 32 data values */ \
\
rsi += -8*64+32*64;	/* Load source pointer */ \
rdi += XMM_SCD;		/* Next sine/cosine pointer */ \
/* 19/23/27/31 iterations of 8 */ \
for(unsigned int loop = (g->count4 >> 8) & 0xFF; loop != 0; loop--){ \
	for(unsigned int loopB=8; loopB != 0; loopB--){ \
		x4cl_four_complex_fft(rsi, 64, 8*64, 16*64, rdi); \
	} \
	rsi += -8*64+32*64;	/* Next source pointer */ \
	rdi += 2*XMM_SCD;		/* Next sine/cosine pointer */ \
} \
\
/* Do FFT levels 10,11 */ \
/* Values 0-15 is real data, 16-31 is semi-real data, 32-8191 is */ \
/* complex data. */ \
/* */ \
/* On input the 64-byte cache lines hold these data values: */ \
/*	0	1	16	17	32	33	48	49 */ \
/*	2	... */ \
/*	... */ \
/*	14	... */ \
/*	64	... */ \
/*	... */ \
/* On output the 64-byte cache lines hold these data values: */ \
/*	0	1	4	5	8	9	12	13 */ \
/*	2	3	... */ \
/*	16	... */ \
/*	18	... */ \
/*	32	... */ \
/*	... */ \
\
/* Do 2 eight_reals_fft_2 macros */ \
/* Do 2 nop_two_two_complex_fft_2 macros */ \
/* Do 4 four_complex_fft macros */ \
/*	distance between fft data elements is 4 */ \
/*	do 2 macros each processing 32 data values */ \
xlast_two_fft_levels

#define xlast_two_fft_levels \
rsi = (uintptr_t)g->DESTARG;/* Restore source pointer */ \
rdi = (uintptr_t)g->sincos2;/* Load sin/cos pointer */ \
for(unsigned int loop=2; loop != 0; loop--){ \
	x4cl_eight_reals_fft_2(rsi, 64, 2*64, 4*64, rdi); \
} \
\
/* Do 632/760/888/1016 four_complex_fft macros */ \
/*	distance between fft data elements is 4 */ \
/*	do 158/190/222/254 macros each processing 32 data values */ \
\
rsi += -2*64+8*64;		/* Load source pointer */ \
rdi += XMM_SCD;		/* Next sine/cosine pointer */ \
/* 79/95/111/127 iterations of 2 */ \
for(unsigned int loop = g->count4 & 0xFF; loop != 0; loop--){ \
	for(unsigned int loopB=2; loopB != 0; loopB--){ \
		x4cl_four_complex_fft(rsi, 64, 2*64, 4*64, rdi); \
	} \
	rsi += -2*64+8*64;		/* Next source pointer */ \
	rdi += 2*XMM_SCD;		/* Next sine/cosine pointer */ \
} \
xmiddle_123(g);

/* Do FFT levels 12,13 */ \
/* Values 0-3 is real data, 4-7 is semi-real data, 8-8191 is */ \
/* complex data. */ \
/* */ \
/* On input the 64-byte cache lines hold these data values: */ \
/*	0	1	4	5	8	9	12	13 */ \
/*	2	3	... */ \
/*	16	... */ \
/*	18	... */ \
/*	32	... */ \
/*	... */ \
/* On output the 64-byte cache lines hold these data values: */ \
/*	0	1	2	3	4	5	6	7 */ \
/*	8	9	... */ \
/*	... */ \
\
/* Do 1/2 eight_reals_fft_2 macros */ \
/* Do 1/2 nop_two_two_complex_fft_2 macros */ \
/* Do 1 four_complex_fft macros */ \
/*	distance between fft data elements is 1 */ \
\
/* Do 638/766/894/1022 more four_complex_fft macros */ \
/*	distance between fft data elements is 1 */ \
/*	do 319/383/447/511 macros each processing 16 data values */ \


#define xmiddle_123(g); \
rsi = (uintptr_t)g->DESTARG;/* Restore source pointer */ \
rdi = (uintptr_t)g->sincos1;		/* Load sin/cos pointer */ \
\
if(likely(g->ffttype >= 2)){ \
	if(g->ffttype == 2){ \
		hgonepass_xmiddle_2(g); \
	}else{ \
		hgonepass_xmiddle_3(g); \
	} \
}else{ \
	hgonepass_xmiddle_1(g); \
	/*pop	rax;*/			/* Pop return address */ \
	xfft_1_ret; \
	return; \
} \
rsi = (uintptr_t)g->DESTARG;/* Restore source pointer */ \

#define xmiddle_4(g); \
	hgonepass_xmiddle_4(g); \
	rsi = (uintptr_t)g->DESTARG;/* Restore source pointer */ \
	/*ret;*/ \



/* Do inverse FFT levels 10,11 */ \
/* On input the 64-byte cache lines hold these data values: */ \
/*	0	1	2	3	4	5	6	7 */ \
/*	8	9	10	11	... */ \
/*	... */ \
/* On output the 64-byte cache lines hold these data values: */ \
/*	0	1	8	9	16	17	24	25 */ \
/*	2	... */ \
/*	4	... */ \
/*	6	... */ \
/*	32	... */ \
/*	... */ \
\
/* Do 2 eight_reals_unfft_2 macros */ \
/* Do 2 nop_two_two_complex_unfft_2 macros */ \
/*	distance between fft data elements is 4 */ \
/*	do 1 macro processing 32 data values */ \

#define cmn8192(RCX1) \
x4cl_eight_reals_unfft_2(rsi, 4*64, 64, 2*64); \
\
/* Do 636/764/892/1020 four_complex_unfft macros */ \
/*	distance between fft data elements is 4 */ \
/*	do 159/191/223/255 macros each processing 32 data values */ \
\
rdi = (uintptr_t)g->sincos2;/* Load sin/cos pointer */ \
for(unsigned int loop = g->count2 & 0xFF; loop != 0; loop--){ \
	x4cl_four_complex_unfft(rsi, 4*64, 64, 2*64, rdi); \
	rdi += XMM_SCD;		/* Next sine/cosine pointer */ \
} \
rsi = (uintptr_t)g->DESTARG;/* Restore source pointer */ \
\
/* Do inverse FFT levels 8,9 */ \
/* On input the 64-byte cache lines hold these data values: */ \
/*	0	1	8	9	16	17	24	25 */ \
/*	2	... */ \
/*	4	... */ \
/*	6	... */ \
/*	32	... */ \
/*	... */ \
/* On output the 64-byte cache lines hold these data values: */ \
/*	0	1	32	33	64	65	96	97 */ \
/*	2	... */ \
/*	...	... */ \
/*	30	... */ \
/*	128	... */ \
/*	... */ \
\
/* Do 8 eight_reals_unfft_2 macros */ \
/* Do 8 nop_two_two_complex_unfft_2 macros */ \
/*	distance between fft data elements is 16 */ \
/*	do 4 macros each processing 32 data values */ \
\
for(unsigned int loop = 4; loop != 0; loop--){ \
	x4cl_eight_reals_unfft_2(rsi, 64, 4*64, 8*64); \
} \
\
/* Do 624/752/880/1008 four_complex_unfft macros */ \
/*	distance between fft data elements is 16 */ \
/*	do 156/188/220/252 macros each processing 32 data values */ \
\
rdi = (uintptr_t)g->sincos3;		/* Load sin/cos pointer */ \
rsi += -4*64+16*64;	/* Load source pointer */ \
/* 7/15 then 32/48 iterations of 4 */ \
for(unsigned int loop = g->count3 & 0xFF; loop != 0; loop--){ \
	for(unsigned int loopB = 4; loopB != 0; loopB--){ \
		x4cl_four_complex_unfft(rsi, 64, 4*64, 8*64, rdi); \
	} \
	rsi += -4*64+16*64;	/* Next source pointer */ \
	rdi += XMM_SCD;		/* Next sine/cosine pointer */ \
} \
rsi = (uintptr_t)g->DESTARG;/* Restore source pointer */ \
\
/* Do inverse FFT levels 6,7 */ \
/* On input the 64-byte cache lines hold these data values: */ \
/*	0	1	32	33	64	65	96	97 */ \
/*	2	... */ \
/*	...	... */ \
/*	30	... */ \
/*	128	... */ \
/*	... */ \
/* On output the 64-byte cache lines hold these data values: */ \
/*	0	1	128	129	256	257	384	385 */ \
/*	2	... */ \
/*	...	... */ \
/*	126	... */ \
/*	512	... */ \
/*	... */ \
\
/* Do 32 eight_reals_unfft_2 macros */ \
/* Do 32 nop_two_two_complex_unfft_2 macros */ \
/*	distance between fft data elements is 64 */ \
/*	do 16 macros each processing 32 data values */ \
\
for(unsigned int loop = 16; loop != 0; loop--){ \
	x4cl_eight_reals_unfft_2(rsi, 64, 16*64, 32*64); \
} \
\
/* Do 584/704/832/960 four_complex_unfft macros */ \
/*	distance between fft data elements is 16 */ \
/*	do 144/176/208/240 macros each processing 32 data values */ \
\
rdi = (uintptr_t)g->sincos4;		/* Load sin/cos pointer */ \
rsi += -16*64+64*64;	/* Load source pointer */ \
/* 1/3 then 8/12 iterations of 16 */ \
for(unsigned int loop = (g->count3 >> 8) & 0xFF; loop != 0; loop--){ \
	for(unsigned int loopB = 16; loopB != 0; loopB--){ \
		x4cl_four_complex_unfft(rsi, 64, 16*64, 32*64, rdi); \
	} \
	rsi += -16*64+64*64;	/* Next source pointer */ \
	rdi += XMM_SCD;		/* Next sine/cosine pointer */ \
} \
rsi = (uintptr_t)g->DESTARG;/* Restore source pointer */ \
\
/* Do inverse FFT levels 4,5 */ \
/* On input the 64-byte cache lines hold these data values: */ \
/*	0	1	128	129	256	257	384	385 */ \
/*	2	... */ \
/*	...	... */ \
/*	126	... */ \
/*	512	... */ \
/*	... */ \
/* On output the 64-byte cache lines hold these data values: */ \
/*	0	1	512	513	1024	1025	1536	1537 */ \
/*	2	... */ \
/*	...	... */ \
/*	510	... */ \
/*	2048	... */ \
/*	... */ \
\
/* Do 512/768 four_complex_unfft macros */ \
/*	distance between fft data elements is 256 */ \
/*	do 128/192 macros each processing 32 data values */ \
\
rsi += RCX1;		/* Next source pointer */ \
rdi = (uintptr_t)g->sincos5;		/* Load sin/cos pointer */ \
/* 2/3 iterations of 64 */ \
for(unsigned int loop = g->count5 & 0xFF; loop != 0; loop--){ \
	for(unsigned int loopB = 64; loopB != 0; loopB--){ \
		x4cl_four_complex_unfft(rsi, 64, 64*64, 128*64, rdi); \
	} \
	rsi += -64*64+256*64;	/* Next source pointer */ \
	rdi += XMM_SCD;		/* Next sine/cosine pointer */ \
} \
rsi = (uintptr_t)g->DESTARG;/* Next source pointer */ \
\
/* This ends the common inverse FFT code */ \
\
/*ret; */ \
/*push_amt = 0;*/ \



/* Perform a 8192-element negacyclic FFT */ \


#define xfft8192p(type) \
\
/* Do a multiply with pre-FFTed inputs */ \
\
if(g->ffttype == 4) { \
	xmiddle_4p(g); \
}else{ \
	\
	/* Do FFT levels 1,2,3 */ \
	/* */ \
	/* On input the 64-byte cache lines hold these data values: */ \
	/*	0	2048	1	2049	4096	6144	4097	6145 */ \
	/*	2	... */ \
	/*	... */ \
	/*	2046 */ \
	/* On output the 64-byte cache lines hold these data values: */ \
	/*	0	1	1024	1025	2048	2049	3072	3073 */ \
	/*	2	... */ \
	/*	... */ \
	/*	1022	... */ \
	/*	4096	... */ \
	/*	... */ \
	\
	/* Do 1024 four_complex_first_fft macros */ \
	/*	distance between fft data elements is 1024 */ \
	/*	do 512 macros each processing 16 data values */ \
	\
	rdi = (uintptr_t)g->plus1_premults;/* Address of premultiplier table */ \
	for(unsigned int loop=512; loop != 0; loop--){ \
		s2cl_four_complex_first_fft(rsi, 64, 512*64); \
		rdi += 128;		/* Next premultiplier pointer */ \
	} \
	rsi += -512*64;		/* Restore source pointer */ \
	\
	/* Do FFT levels 4,5 */ \
	/* */ \
	/* On input the 64-byte cache lines hold these data values: */ \
	/*	0	1	1024	1025	2048	2049	3072	3073 */ \
	/*	2	... */ \
	/*	... */ \
	/*	1022	... */ \
	/*	4096	... */ \
	/*	... */ \
	/* On output the 64-byte cache lines hold these data values: */ \
	/*	0	1	256	257	512	513	768	769 */ \
	/*	2	... */ \
	/*	... */ \
	/*	254	... */ \
	/*	1024	... */ \
	/*	... */ \
	\
	/* Do 1024 four_complex_fft macros */ \
	/*	distance between fft data elements is 256 */ \
	/*	do 256 macros each processing 32 data values */ \
	\
	rdi = (uintptr_t)g->sincos5;		/* Load sine/cosine pointer */ \
	/* 2 iterations of 128 */ \
	for(unsigned int loop=2; loop != 0; loop--){ \
		for(unsigned int loopB=128; loopB != 0; loopB--){ \
			x4cl_four_complex_fft(rsi, 64, 128*64, 256*64, rdi); \
		} \
		rsi += -128*64+512*64;	/* Next source pointer */ \
		rdi += 2*XMM_SCD;		/* Next sine/cosine pointer */ \
	} \
	rsi = (uintptr_t)g->DESTARG;/* Restore source pointer */ \
	\
	/* Call common code to do last levels of FFT */ \
	\
	xlast_six_fft_levels_p; \
} \
cmn8192p; \
\
/* Do FFT levels 6,7 */ \
/* */ \
/* On input the 64-byte cache lines hold these data values: */ \
/*	0	1	256	257	512	513	768	769 */ \
/*	2	... */ \
/*	... */ \
/*	254	... */ \
/*	1024	... */ \
/*	... */ \
/* On output the 64-byte cache lines hold these data values: */ \
/*	0	1	64	65	128	129	192	193 */ \
/*	2	... */ \
/*	... */ \
/*	62	... */ \
/*	256	... */ \
/*	... */ \
\
/* Do inverse FFT levels 1,2,3 */ \
/* On input the 64-byte cache lines hold these data values: */ \
/*	0	1	512	513	1024	1025	1536	1537 */ \
/*	2	... */ \
/*	...	... */ \
/*	510	... */ \
/*	2048	... */ \
/*	... */ \
/* On output the 64-byte cache lines hold these data values: */ \
/*	0	1	2048	2049	4096	4097	6144	6145 */ \
/*	2	... */ \
/*	... */ \
/*	2046 */ \
\
/* Do 1024 four_complex_last_unfft macros */ \
/*	distance between fft data elements is 1024 */ \
/*	do 256 macros each processing 32 data values */ \
\
rdi = (uintptr_t)g->plus1_premults;/* Address of premultiplier table */ \
for(unsigned int loop=256; loop != 0; loop--){ \
	x4cl_four_complex_last_unfft(rsi, 64, 256*64, 512*64, 256*128); \
	rdi += 128;		/* Next premultiplier pointer */ \
} \
\
xfft_3_ret;

/* Common negacyclic FFT code */ \

/* Do 1024 four_complex_fft macros */ \
/*	distance between fft data elements is 64 */ \
/*	do 256 macros each processing 32 data values */ \

#define xlast_six_fft_levels_p \
/*push_amt = SZPTR;*/ \
rdi = (uintptr_t)g->sincos4;		/* Load sin/cos pointer */ \
/* 8 iterations of 32 */ \
for(unsigned int loop = (g->count5 >> 8) & 0xFF; loop != 0; loop--){ \
	for(unsigned int loopB=32; loopB != 0; loopB--){ \
		x4cl_four_complex_fft(rsi, 64, 32*64, 64*64, rdi); \
	} \
	rsi += -32*64+128*64;	/* Next source pointer */ \
	rdi += 2*XMM_SCD;		/* Next sine/cosine pointer */ \
} \
rsi = (uintptr_t)g->DESTARG;/* Next source pointer */ \
\
/* Do FFT levels 8,9 */ \
/* */ \
/* On input the 64-byte cache lines hold these data values: */ \
/*	0	1	64	65	128	129	192	193 */ \
/*	2	... */ \
/*	... */ \
/*	62	... */ \
/*	256	... */ \
/*	... */ \
/* On output the 64-byte cache lines hold these data values: */ \
/*	0	1	16	17	32	33	48	49 */ \
/*	2	... */ \
/*	... */ \
/*	14	... */ \
/*	64	... */ \
/*	... */ \
\
/* Do 1024 four_complex_fft macros */ \
/*	distance between fft data elements is 16 */ \
/*	do 256 macros each processing 32 data values */ \
\
xlast_four_fft_levels_p;

#define xlast_four_fft_levels_p \
rdi = (uintptr_t)g->sincos3;		/* Load sin/cos pointer */ \
/* 32 iterations of 8 */ \
for(unsigned int loop = (g->count4 >> 8) & 0xFF; loop != 0; loop--){ \
	for(unsigned int loopB=8; loopB != 0; loopB--){ \
		x4cl_four_complex_fft(rsi, 64, 8*64, 16*64, rdi); \
	} \
	rsi += -8*64+32*64;	/* Next source pointer */ \
	rdi += 2*XMM_SCD;		/* Next sine/cosine pointer */ \
} \
rsi = (uintptr_t)g->DESTARG;/* Restore source pointer */ \
\
/* Do FFT levels 10,11 */ \
/* */ \
/* On input the 64-byte cache lines hold these data values: */ \
/*	0	1	16	17	32	33	48	49 */ \
/*	2	... */ \
/*	... */ \
/*	14	... */ \
/*	64	... */ \
/*	... */ \
/* On output the 64-byte cache lines hold these data values: */ \
/*	0	1	4	5	8	9	12	13 */ \
/*	2	3	... */ \
/*	16	... */ \
/*	18	... */ \
/*	32	... */ \
/*	... */ \
\
/* Do 1024 four_complex_fft macros */ \
/*	distance between fft data elements is 4 */ \
/*	do 256 macros each processing 32 data values */ \
\
xlast_two_fft_levels_p;

#define xlast_two_fft_levels_p \
rdi = (uintptr_t)g->sincos2;/* Load sin/cos pointer */ \
/* 128 iterations of 2 */ \
for(unsigned int loop = g->count4 & 0xFF; loop != 0; loop--){ \
	for(unsigned int loopB=2; loopB != 0; loopB--){ \
		x4cl_four_complex_fft(rsi, 64, 2*64, 4*64, rdi); \
	} \
	rsi += -2*64+8*64;		/* Next source pointer */ \
	rdi += 2*XMM_SCD;		/* Next sine/cosine pointer */ \
} \
xmiddle_123p(g)

/* Do FFT levels 12,13 */ \
/* */ \
/* On input the 64-byte cache lines hold these data values: */ \
/*	0	1	4	5	8	9	12	13 */ \
/*	2	3	... */ \
/*	16	... */ \
/*	18	... */ \
/*	32	... */ \
/*	... */ \
/* On output the 64-byte cache lines hold these data values: */ \
/*	0	1	2	3	4	5	6	7 */ \
/*	8	9	... */ \
/*	... */ \
\
/* Do 1024 four_complex_fft macros */ \
/*	distance between fft data elements is 1 */ \
/*	do 512 macros each processing 16 data values */ \

#if 0
void xmiddle_123p_f1(struct gwasm_data *__restrict g);
void xmiddle_4p_f1(struct gwasm_data *__restrict g);
#define xmiddle_123p(g) xmiddle_123p_f1(g)
#define xmiddle_4p(g) xmiddle_4p_f1(g)
#else
#define xmiddle_123p(g); \
/*rsi = (uintptr_t)g->DESTARG;*/ /* Restore source pointer */ \
/*rdi = (uintptr_t)g->sincos1;*/ /* Load sin/cos pointer */ \
/* 512 iterations */ \
\
if(__builtin_expect(g->ffttype >= 2, 1)){ \
	if(g->ffttype == 2){ \
		hgonepass_xmiddle_2p(g); \
	}else{ \
		hgonepass_xmiddle_3p(g); \
	} /* ret */ \
}else{ \
	hgonepass_xmiddle_1p(g); \
	/*pop	rax;*/			/* Pop return address */ \
	xfft_1_ret; \
	return; \
} \
rsi = (uintptr_t)g->DESTARG;/* Restore source pointer */ \

#define xmiddle_4p(g); \
	hgonepass_xmiddle_4p(g); \
	rsi = (uintptr_t)g->DESTARG;/* Restore source pointer */ \
	/* ret */
#endif

/* Do inverse FFT levels 10,11 */ \
/* On input the 64-byte cache lines hold these data values: */ \
/*	0	1	2	3	4	5	6	7 */ \
/*	8	9	10	11	... */ \
/*	... */ \
/* On output the 64-byte cache lines hold these data values: */ \
/*	0	1	8	9	16	17	24	25 */ \
/*	2	... */ \
/*	4	... */ \
/*	6	... */ \
/*	32	... */ \
/*	... */ \
\
/* Do 1024 four_complex_unfft macros */ \
/*	distance between fft data elements is 4 */ \
/*	do 256 macros each processing 32 data values */ \

#define cmn8192p \
rdi = (uintptr_t)g->sincos2;/* Load sin/cos pointer */ \
/* 256 iterations */ \
for(unsigned int loop = g->count2; loop != 0; loop--){ \
	x4cl_four_complex_unfft(rsi, 4*64, 64, 2*64, rdi); \
	rdi += XMM_SCD;		/* Next sine/cosine pointer */ \
} \
rsi = (uintptr_t)g->DESTARG;/* Restore source pointer */ \
\
/* Do inverse FFT levels 8,9 */ \
/* On input the 64-byte cache lines hold these data values: */ \
/*	0	1	8	9	16	17	24	25 */ \
/*	2	... */ \
/*	4	... */ \
/*	6	... */ \
/*	32	... */ \
/*	... */ \
/* On output the 64-byte cache lines hold these data values: */ \
/*	0	1	32	33	64	65	96	97 */ \
/*	2	... */ \
/*	...	... */ \
/*	30	... */ \
/*	128	... */ \
/*	... */ \
\
/* Do 1024 four_complex_unfft macros */ \
/*	distance between fft data elements is 16 */ \
/*	do 256 macros each processing 32 data values */ \
\
rdi = (uintptr_t)g->sincos3;		/* Load sin/cos pointer */ \
/* 64 iterations of 4 */ \
for(unsigned int loop = g->count3 & 0xFF; loop != 0; loop--){ \
	for(unsigned int loopB=4; loopB != 0; loopB--){ \
		x4cl_four_complex_unfft(rsi, 64, 4*64, 8*64, rdi); \
	} \
	rsi += -4*64+16*64;	/* Next source pointer */ \
	rdi += XMM_SCD;		/* Next sine/cosine pointer */ \
} \
rsi = (uintptr_t)g->DESTARG;/* Restore source pointer */ \
\
/* Do inverse FFT levels 6,7 */ \
/* On input the 64-byte cache lines hold these data values: */ \
/*	0	1	32	33	64	65	96	97 */ \
/*	2	... */ \
/*	...	... */ \
/*	30	... */ \
/*	128	... */ \
/*	... */ \
/* On output the 64-byte cache lines hold these data values: */ \
/*	0	1	128	129	256	257	384	385 */ \
/*	2	... */ \
/*	...	... */ \
/*	126	... */ \
/*	512	... */ \
/*	... */ \
\
/* Do 1024 four_complex_unfft macros */ \
/*	distance between fft data elements is 16 */ \
/*	do 256 macros each processing 32 data values */ \
\
rdi = (uintptr_t)g->sincos4;		/* Load sin/cos pointer */ \
/* 16 iterations of 16 */ \
for(unsigned int loop = (g->count3 >> 8) & 0xFF; loop != 0; loop--){ \
	for(unsigned int loopB=16; loopB != 0; loopB--){ \
		x4cl_four_complex_unfft(rsi, 64, 16*64, 32*64, rdi); \
	} \
	rsi += -16*64+64*64;	/* Next source pointer */ \
	rdi += XMM_SCD;		/* Next sine/cosine pointer */ \
} \
rsi = (uintptr_t)g->DESTARG;/* Restore source pointer */ \
\
/* Do inverse FFT levels 4,5 */ \
/* On input the 64-byte cache lines hold these data values: */ \
/*	0	1	128	129	256	257	384	385 */ \
/*	2	... */ \
/*	...	... */ \
/*	126	... */ \
/*	512	... */ \
/*	... */ \
/* On output the 64-byte cache lines hold these data values: */ \
/*	0	1	512	513	1024	1025	1536	1537 */ \
/*	2	... */ \
/*	...	... */ \
/*	510	... */ \
/*	2048	... */ \
/*	... */ \
\
/* Do 1024 four_complex_unfft macros */ \
/*	distance between fft data elements is 256 */ \
/*	do 256 macros each processing 32 data values */ \
\
rdi = (uintptr_t)g->sincos5;		/* Load sin/cos pointer */ \
/* 4 iterations of 64 */ \
for(unsigned int loop = g->count5 & 0xFF; loop != 0; loop--){ \
	for(unsigned int loopB=64; loopB != 0; loopB--){ \
		x4cl_four_complex_unfft(rsi, 64, 64*64, 128*64, rdi); \
	} \
	rsi += -64*64+256*64;	/* Next source pointer */ \
	rdi += XMM_SCD;		/* Next sine/cosine pointer */ \
} \
rsi = (uintptr_t)g->DESTARG;/* Next source pointer */ \
\
/* This ends the common inverse FFT code */ \
\
/* ret; */
