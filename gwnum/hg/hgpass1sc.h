#pragma once

#define EXEC11(...) __VA_ARGS__
#define EXEC12(...)
#define EXEC14(...)
#define EXEC18(...)
#define EXEC2481(...)
#define EXEC2482(...) __VA_ARGS__
#define EXEC2484(...) __VA_ARGS__
#define EXEC2488(...) __VA_ARGS__


#define xfft_header_hg2(procname, xpass2name, pass2_levels, clm) \
void procname##_i (struct gwasm_data *__restrict, int); \
void procname##_i0 (struct gwasm_data *__restrict g) { procname##_i(g, 0);} \
void procname##_i1 (struct gwasm_data *__restrict g) { procname##_i(g, 1);} \
void procname(struct gwasm_data *__restrict g){ \
	uintptr_t rbx,rsi; \
	rsi = (uintptr_t)g->DESTARG; \
	rbx = g->DIST_TO_FFTSRCARG; \
	clear_timers; \
	if(*(int*)(rsi-28+rbx) == 0) { /* Test FFT-started flag */ \
		/* pass1: pass1_forward_fft_setup(b0b) */ \
		g->next_block = 0; \
		start_timer(1); \
		g->thread_work_routine = procname##_i1; \
		pass1_wake_up_threads(g); \
		procname##_i (g, 1); \
	}else{ \
		procname##_i (g, 4); /* Start pass2 */ \
	} \
} \
void procname##_i (struct gwasm_data *__restrict g, int branch1) \
{ \
	uintptr_t rcx,rbx,rsi,rdi; \
	vec2f64 xmm0,xmm1,xmm2,xmm3,xmm4,xmm5,xmm6,xmm7; \
	rsi = (uintptr_t)g->DESTARG; \
	rbx = g->DIST_TO_FFTSRCARG; \
	/* xmult.mac  set_FFT_constants */ \
	const uintptr_t clmblkdst = clm*128; \
	const uintptr_t clmblkdst8 = clmblkdst*8+128; \
	const uintptr_t p2cl = 1<<(pass2_levels-2); \
	const uintptr_t blkdst_tmp1 = (p2cl*128+p2cl/64*128); \
	const uintptr_t blkdst = (blkdst_tmp1 % 256 == 0) ? (blkdst_tmp1-128) : blkdst_tmp1; \
	do{ \
		if(branch1 < 1) goto c0b; /* Do another inverse FFT/norm/FFT  */ \
		if(branch1 == 1) goto b0b; /* Do another forward FFT */ \

#define xfft_footer_hg2 \
	; }while(1); \
}



/* Copyright 2001-2023 - Mersenne Research, Inc.  All rights reserved */ \
/* Author:  George Woltman */ \
/* Email: woltman@alum.mit.edu */ \
/* */ \
/* These macros do the first pass of 2 pass FFTs using SSE2 instructions. */ \
/* */ \
\
/* ******************************************************** */ \
/* ******************************************************** */ \
/* ******************  PASS 1 MACROS  ********************* */ \
/* ******************************************************** */ \
/* ******************************************************** */ \
\
/************************************************************************ */ \
/* Pass 1 macros */ \
/* */ \
/* p2cl = number of 128-byte cache lines processed from a pass 2 block */ \
/* clm = cache line multiplier (num cache lines processed each iteration) */ \
/************************************************************************ */ \
\
/* Note:  All memory layouts are for 64 byte L1 cache lines.  The high */ \
/* order XMM words are not shown (they are 1K higher than the low */ \
/* order XMM word.  Since the L2 cache reads 128 byte cache lines, two */ \
/* 64 byte cache lines are processed at a time.  The second 64 byte */ \
/* cache line contains a completely independent set of data. */ \
\
\
/* Perform a 320K-element FFT.  This is done in two passes.  Pass 1 */ \
/* does 8 levels, pass 2 uses common code to perform the last 11 levels. */ \
\

#define hg_pass1sclevels8pfa5(procname, xpass2name, pass2_levels, clm) \
 /* (fftname, pass2_levels, clmarg) */ \
\
/* set_FFT_constants(fftname, 0, 0, pass2_levels, clmarg); */ \
xfft_header_hg2(procname, xpass2name, pass2_levels, clm); \
\
/* Do the last 11 FFT levels */ \
/*	do 1 pass2_11_levels_real and 39 pass2_11_levels_complex */ \
\
/* pass2: */	start_timer(0); \
xpass2name(g); \
end_timer(0); \
if(g->ffttype == 1)		/* We're done if FFTing only */ \
	return; \
\
/* Do the inverse FFT */ \
\
pass1_inverse_fft_setup(c0b, procname); \
\
/* Do inverse FFT levels 7,8 */ \
/* On input the 64-byte cache lines hold these data values: */ \
/*	0K	2K	4K	6K */ \
/*	8K	10K	12K	14K */ \
/*	... */ \
/* On output the 64-byte cache lines hold these data values: */ \
/*	0K	4K	8K	12K */ \
/*	2K	... */ \
/*	16K	... */ \
/*	... */ \
\
/* Do 1 eight_reals_unfft_2 macros */ \
/* Do 1 nop_two_two_complex_unfft_2 macros */ \
/*	distance between fft data elements is 2K */ \
/*	do 1 macros each processing 16 data values */ \
\
c0b:	start_timer(16); \
rsi = (uintptr_t)g->data_addr;/* Get source address */ \
rcx = (uintptr_t)g->scratch_area;	/* Get address of scratch area */ \
pfing(rbp = (uintptr_t)g->data_prefetch);/* Set prefetch distance */ \
pfing(rbp -= rsi); \
for(unsigned int loopA = 2*clm; loopA; loopA--) { \
	g2cl_eight_reals_unfft_2(rsi, 64, blkdst, rcx, 64, clmblkdst); \
} \
\
/* Do 38 four_complex_unfft macros */ \
/*	distance between fft data elements is 2K */ \
/*	do 19 macros each processing 16 data values */ \
\
rsi += -2*clm*64+2*blkdst;	/* Next source pointer */ \
rcx += -2*clm*64+2*clmblkdst; /* Next dest pointer */ \
rdi = (uintptr_t)g->sincos1;/* Load sin/cos pointer */ \
/* al = 256/4+5; */ /* 19 iterations of 2*clm */ \
for(int loopA = 19-1; loopA >= 0; loopA--) { \
	for(unsigned int loopC = 2*clm; loopC; loopC--) { \
		g2cl_four_complex_unfft(rsi, 64, blkdst, rcx, 64, clmblkdst); \
	} \
	rsi += -2*clm*64+2*blkdst;	/* Next source pointer */ \
	rcx += -2*clm*64+2*clmblkdst; /* Next dest pointer */ \
	rdi += XMM_SCD;		/* Next sine/cosine pointer */ \
	if((loopA & 3) == 0) \
		rcx += -8*clmblkdst+clmblkdst8; /* Next dest pointer */ \
} \
pfing(rcx = rsi-20*2*blkdst);	/* Create prefetch pointer */ \
rsi = (uintptr_t)g->scratch_area;	/* Restore source pointer */ \
end_timer(16); \
\
/* Do inverse FFT levels 6 */ \
/* On input the 64-byte cache lines hold these data values: */ \
/*	0K	4K	8K	12K */ \
/*	2K	... */ \
/*	16K	... */ \
/*	... */ \
/* On output the 64-byte cache lines hold these data values: */ \
/*	0K	8K	16K	24K */ \
/*	2K	... */ \
/*	... */ \
/*	6K */ \
/*	32K	... */ \
/*	... */ \
/* */ \
/* Do 4 eight_reals_unfft_1 macros */ \
/*	distance between fft data elements is 8K */ \
/*	do 2 macros each processing 16 data values */ \
/* */ \
/* Do 36 two_two_complex_unfft macros */ \
/*	distance between fft data elements is 8K */ \
/*	do 18 macros each processing 16 data values */ \
\
/* Do inverse FFT levels 4,5 */ \
/* On input the 64-byte cache lines hold these data values: */ \
/*	0K	8K	16K	24K */ \
/*	2K	... */ \
/*	... */ \
/*	6K */ \
/*	32K	... */ \
/*	... */ \
/* On output the 64-byte cache lines hold these data values: */ \
/*	0K	16K	32K	48K */ \
/*	2K	... */ \
/*	... */ \
/*	14K */ \
/*	64K	96K	128K	160K */ \
/*	... */ \
/* */ \
/* Do 8 eight_reals_unfft_2 macros */ \
/*	distance between fft data elements is 16K */ \
/*	do 4 macros each processing 16 data values */ \
/* */ \
/* Do 32 four_complex_unfft macros */ \
/*	distance between fft data elements is 16K */ \
/*	do 8 macros each processing 32 data values */ \
\
/* Do inverse FFT levels 1,2,3 */ \
/* On input the 64-byte cache lines hold these data values: */ \
/*	0K	16K	32K	48K */ \
/*	2K	... */ \
/*	... */ \
/*	14K */ \
/*	64K	96K	128K	160K */ \
/*	... */ \
/* On output the 64-byte cache lines hold these data values: */ \
/*	0K	16K	160K	176K */ \
/*	2K	... */ \
/*	... */ \
/*	14K */ \
/*	32K	96K	192K	256K */ \
/*	34K */ \
/*	... */ \
/*	94K */ \
/* */ \
/* Do 64 five_reals_last_unfft macros */ \
/*	distance between fft data elements is 64K */ \
/*	do 8 macros each processing 40 data values */ \
\
start_timer(17); \
rdi = (uintptr_t)g->sincos2;/* Load sin/cos pointer */ \
rbx = (uintptr_t)g->sincos3;/* Load sin/cos pointer */ \
/* 2 iterations */ \
\
for(unsigned int loopA = 2; loopA; loopA--) { \
	xtouchlines(10);			/* Preload TLBs */ \
	\
	for(unsigned int loopB = clm; loopB; loopB--) { \
		xloop1sp(c4b, 4*clmblkdst, x2cl_eight_reals_unfft_1(rsi, 64, 2*clmblkdst)); \
		\
		/* ah = 256/2+5; */		/* 9 iterations */ \
		for(unsigned int loopC = 9; loopC; ) { \
			do{ \
				for(unsigned int loopE = 2; loopE; loopE--) { \
					x2cl_two_complex_unfft(rsi, 64, 2*clmblkdst); \
				} \
				rsi += -2*64+4*clmblkdst;	/* Next source pointer */ \
				rdi += XMM_SCD;		/* Next sine/cosine pointer */ \
			}while(--loopC & 1); \
			rsi += -2*4*clmblkdst+clmblkdst8; /* Next source pointer */ \
		} \
		rdi += -9*XMM_SCD;		/* Restore sine/cosine pointer */ \
		L2prefetch128(u8ptr(rcx+rbp)); \
		pfing(rcx += blkdst);		/* Next prefetch pointer */ \
		\
		rsi += -5*clmblkdst8;	/* Restore source pointer */ \
		for(unsigned int loopC = 2; loopC; loopC--) { \
			for(unsigned int loopD = 2; loopD; loopD--) { \
				x2cl_half_eight_reals_unfft_2(rsi, 64, 4*clmblkdst); \
			} \
			rsi += -2*64+2*clmblkdst; \
		} \
		\
		rsi += -2*2*clmblkdst+clmblkdst8; /* Next source pointer */ \
		for(unsigned int loopC = 2; loopC; loopC--) { \
			for(unsigned int loopD = 2; loopD; loopD--) { \
				xloop1sp(c3b, 2*clmblkdst, x4cl_four_complex_unfft(rsi, 64, 4*clmblkdst, clmblkdst8, rbx)); \
			} \
			rsi += -2*2*clmblkdst+2*clmblkdst8; \
			rbx += XMM_SCD;		/* Next sine/cosine pointer */ \
		} \
		rbx += -2*XMM_SCD;		/* Restore sine/cosine pointer */ \
		\
		rsi += -2*2*clmblkdst8-clmblkdst8; /* Restore source pointer */ \
		xfive_reals_unfft_preload;	/* Preload consts for five real macros */ \
		for(unsigned int loopC = 4; loopC; loopC--) { \
			; \
			xloop1sp(c1b, 2*clmblkdst, x5cl_five_reals_last_unfft(rsi, 64, clmblkdst8)); \
		} \
		\
		EXEC248##clm(rsi += -4*2*clmblkdst+2*64);/* Next source pointer */ \
		EXEC248##clm(pfing(rcx += -10*blkdst+128));	/* Next prefetch pointer */ \
	} \
	EXEC248##clm(rsi += -clm*2*64+clmblkdst);/* Next source pointer */ \
	EXEC248##clm(pfing(rcx += -clm*128+10*blkdst));	/* Next prefetch pointer */ \
	EXEC1##clm(rsi += -4*2*clmblkdst+clmblkdst);/* Next source pointer */ \
	\
} \
end_timer(17); \
\
/* Normalize these values */ \
\
pass1_normalize(1, yes_fft); \
\
/* Copy data from the scratch area back to the FFT data area */ \
\
copy_scratch_data_to_fft(40, clm); \
goto	endlp;			/* Skip FFT code if flag not set */ \
\
/* Fire up auxiliary threads to do pass 1 forward FFT */ \
\
/* pass1: */	pass1_forward_fft_setup(b0b); \
\
/* Copy data from the FFT data area back to the scratch area */ \
\
b0b:	copy_fft_data_to_scratch(40, clm); \
\
/* Do FFT levels 1,2,3 */ \
/* */ \
/* On input the 64-byte cache lines hold these data values: */ \
/*	0K	16K	160K	176K */ \
/*	2K	... */ \
/*	... */ \
/*	14K */ \
/*	32K	96K	192K	256K */ \
/*	34K */ \
/*	... */ \
/*	94K */ \
/* On output the 64-byte cache lines hold these data values: */ \
/*	0K	16K	32K	48K */ \
/*	2K	... */ \
/*	... */ \
/*	14K */ \
/*	64K	128K	192K	256K */ \
/*	66K */ \
/*	... */ \
/*	126K */ \
/* */ \
/* Do 64 five_reals_first_fft macros (FFT levels 1,2,3) */ \
/*	distance between fft data elements is 64K */ \
/*	do 8 iterations each processing 40 data values */ \
\
/* Do FFT levels 4,5 */ \
/* Values 0K-63K is real data, 64K-319K is complex data. */ \
/* */ \
/* On input the 64-byte cache lines hold these data values: */ \
/*	0K	16K	32K	48K */ \
/*	2K	... */ \
/*	... */ \
/*	14K */ \
/*	64K	128K	192K	256K */ \
/*	66K */ \
/*	... */ \
/*	126K */ \
/* On output the 64-byte cache lines hold these data values: */ \
/*	0K	16K	32K	48K */ \
/*	2K	... */ \
/*	... */ \
/*	14K	... */ \
/*	64K	... */ \
/*	... */ \
/* */ \
/* Do 8 eight_reals_fft_2 macros */ \
/*	distance between fft data elements is 16K */ \
/*	do 8 macros each processing 8 data values */ \
/* */ \
/* Do 32 four_complex_fft macros */ \
/*	distance between fft data elements is 16K */ \
/*	do 8 macros each processing 32 data values */ \
\
/* Do FFT levels 6 */ \
/* Values 0K-15K is real data, 16K-31K is semi-real data, 32K-319K is */ \
/* complex data. */ \
/* */ \
/* On input the 64-byte cache lines hold these data values: */ \
/*	0K	16K	32K	48K */ \
/*	2K	... */ \
/*	... */ \
/*	14K	... */ \
/*	64K	... */ \
/*	... */ \
/* On output the 64-byte cache lines hold these data values: */ \
/*	0K	8K	16K	24K */ \
/*	2K	... */ \
/*	... */ \
/*	6K	... */ \
/*	32K	... */ \
/*	... */ \
/* */ \
/* Do 4 eight_reals_fft_1 macros */ \
/* Do 4 two_two_complex_fft macros */ \
/*	distance between fft data elements is 8K */ \
/*	do 4 macros each processing 16 data values */ \
/* */ \
/* Do 32 two_two_complex_fft macros */ \
/*	distance between fft data elements is 8K */ \
/*	do 16 macros each processing 16 data values */ \
\
yes_fft:start_timer(21); \
rsi = (uintptr_t)g->scratch_area;	/* Get address of scratch area */ \
rcx = (uintptr_t)g->data_addr; /* Load source address for prefetching */ \
pfing(rbp = (uintptr_t)g->data_prefetch);/* Set prefetch distance */ \
pfing(rbp -= rcx); \
pfing(rcx += 20*blkdst);		/* Create prefetching pointer */ \
rdi = (uintptr_t)g->sincos3;/* Load sin/cos pointer */ \
rbx = (uintptr_t)g->sincos2;/* Load sin/cos pointer */ \
/* 4 iterations */ \
\
for(unsigned int loopA = 4; loopA; loopA--) { \
	for(unsigned int loopB = clm; loopB; loopB--) { \
		xtouchlines(5);			/* Preload TLBs */ \
		\
		xfive_reals_fft_preload;		/* Preload consts for five real macros */ \
		for(unsigned int loopC = 2; loopC; loopC--) { \
			xloop1sp(b1b, 4*clmblkdst, x5cl_five_reals_fft(rsi, 64, clmblkdst8)); \
		} \
		\
		rsi += -2*4*clmblkdst;	/* Restore source pointer */ \
		for(unsigned int loopC = 2; loopC; loopC--) { \
			for(unsigned int loopD = 2; loopD; loopD--) { \
				x1cl_half_eight_reals_fft_2(rsi, 64); \
			} \
			rsi += -2*64+4*clmblkdst;	/* Next source pointer */ \
		} \
		\
		rsi += -2*4*clmblkdst+clmblkdst8;/* Next source pointer */ \
		for(unsigned int loopC = 2; loopC; loopC--) { \
			xloop1sp(b3b, 4*clmblkdst, x4cl_four_complex_fft(rsi, 64, clmblkdst8, 2*clmblkdst8, rdi)); \
		} \
		\
		rsi += -2*4*clmblkdst-clmblkdst8;/* Restore source pointer */ \
		for(unsigned int loopC = 2; loopC; loopC--) { \
			x2cl_eight_reals_fft_1(rsi, 64, 4*clmblkdst); \
		} \
		\
		rsi += -2*64+clmblkdst8;	/* Next source pointer */ \
		rbx += XMM_SCD;		/* Next sine/cosine pointer */ \
		for(unsigned int loopC = 4; loopC; loopC--) { \
			for(unsigned int loopD = 2; loopD; loopD--) { \
				x2cl_two_complex_fft(rsi, 64, 4*clmblkdst, rbx); \
			} \
			rsi += -2*64+clmblkdst8;	/* Next source pointer */ \
			rbx += 2*XMM_SCD;		/* Next sine/cosine pointer */ \
		} \
		L2prefetch128(u8ptr(rcx+rbp)); \
		pfing(rcx += blkdst);		/* Next prefetch pointer */ \
		rbx += -9*XMM_SCD;		/* Restore sine/cosine pointer */ \
		\
		EXEC248##clm(rsi += -5*clmblkdst8+2*64);	/* Next source pointer */ \
		EXEC248##clm(rcx += -5*blkdst+128);	/* Next prefetch pointer */ \
	} \
	EXEC248##clm(rsi += -clm*2*64+clmblkdst);/* Next source pointer */ \
	EXEC248##clm(rcx += -clm*128+5*blkdst);	/* Next prefetch pointer */ \
	EXEC1##clm(rsi += -5*clmblkdst8+clmblkdst);/* Next source pointer */ \
	\
} \
rsi += -4*clmblkdst;	/* Restore source pointer */ \
pfing(rcx += -40*blkdst);		/* Restore prefetch pointer */ \
end_timer(21); \
\
/* Do FFT levels 7,8 */ \
/* Values 0K-7K is real data, 8K-15K is semi-real data, 16K-319K is */ \
/* complex data. */ \
/* */ \
/* On input the 64-byte cache lines hold these data values: */ \
/*	0K	8K	16K	24K */ \
/*	2K	... */ \
/*	... */ \
/*	6K	... */ \
/*	32K	... */ \
/*	... */ \
/* On output the 64-byte cache lines hold these data values: */ \
/*	0K	2K	4K	6K */ \
/*	8K	... */ \
/*	... */ \
\
/* Do 1 eight_reals_fft_2 macros */ \
/* Do 1 nop_two_two_complex_fft_2 macros */ \
/* Do 2 four_complex_fft macros */ \
/*	distance between fft data elements is 2K */ \
/*	do 1 macros each processing 32 data values */ \
\
start_timer(24); \
rdi = (uintptr_t)g->sincos1;/* Load sin/cos pointer */ \
for(unsigned int loopA = 2*clm; loopA; loopA--) { \
	g4cl_eight_reals_fft_2(rsi, 64, clmblkdst, 2*clmblkdst, rcx, 64, blkdst, 2*blkdst); \
} \
\
/* Do 36 four_complex_fft macros */ \
/*	distance between fft data elements is 2K */ \
/*	do 9 macros each processing 32 data values */ \
\
rsi += -2*clm*64+4*clmblkdst;/* Next source pointer */ \
rcx += -2*clm*64+4*blkdst;	/* Next dest pointer */ \
rdi += XMM_SCD;		/* Next sine/cosine pointer */ \
/*eax = 5*256+1;*/ /* 9 iterations of 2*clm */ \
for(unsigned int loopA = 9; loopA; ) { \
	do{ \
		for(unsigned int loopC = clm; loopC; loopC--) { \
			for(unsigned int loopD = 2; loopD; loopD--) { \
				g4cl_four_complex_fft(rsi, 64, clmblkdst, 2*clmblkdst, rcx, 64, blkdst, 2*blkdst); \
				L2prefetch128(u8ptr(rcx-64+rbp)); \
				L2prefetch128(u8ptr(rcx-64+blkdst+rbp)); \
				pfing(rbp += -64+2*blkdst); \
			} \
			pfing(rbp += 128-4*blkdst); \
		} \
		rsi += -2*clm*64+4*clmblkdst;/* Next source pointer */ \
		rcx += -2*clm*64+4*blkdst;	/* Next dest pointer */ \
		rdi += 2*XMM_SCD;		/* Next sine/cosine pointer */ \
	}while(--loopA & 1); \
	rsi += -8*clmblkdst+clmblkdst8; \
} \
end_timer(24); \
\
/* Work on next 4 sets of 160 values */ \
\
endlp:	xpass1_get_next_block(pass2, c0b, b0b); /* Get next block and jump */ \
/* done: */	xfft_footer_hg2; \
\
\
\
/* Perform a 384K-element FFT.  This is done in two passes.  Pass 1 */ \
/* does 8 levels, pass 2 uses common code to perform the last 11 levels. */ \
\

#define hg_pass1sclevels8pfa6(procname, xpass2name, pass2_levels, clm) \
 /* (fftname, pass2_levels, clmarg) */ \
\
/* set_FFT_constants(fftname, 0, 0, pass2_levels, clmarg); */ \
xfft_header_hg2(procname, xpass2name, pass2_levels, clm); \
\
/* Do the last 11 FFT levels */ \
/*	do 1 pass2_11_levels_real and 47 pass2_11_levels_complex */ \
\
/* pass2: */	start_timer(0); \
xpass2name(g); \
end_timer(0); \
if(g->ffttype == 1)		/* We're done if FFTing only */ \
	return; \
\
/* Do the inverse FFT */ \
\
pass1_inverse_fft_setup(c0b, procname); \
\
/* Do inverse FFT levels 7,8 */ \
/* On input the 64-byte cache lines hold these data values: */ \
/*	0K	2K	4K	6K */ \
/*	8K	10K	12K	14K */ \
/*	... */ \
/* On output the 64-byte cache lines hold these data values: */ \
/*	0K	4K	8K	12K */ \
/*	2K	... */ \
/*	16K	... */ \
/*	... */ \
\
/* Do 1 eight_reals_unfft_2 macros */ \
/* Do 1 nop_two_two_complex_unfft_2 macros */ \
/*	distance between fft data elements is 2K */ \
/*	do 1 macros each processing 16 data values */ \
\
c0b:	start_timer(16); \
rsi = (uintptr_t)g->data_addr;/* Get source address */ \
rcx = (uintptr_t)g->scratch_area;	/* Get address of scratch area */ \
pfing(rbp = (uintptr_t)g->data_prefetch);/* Set prefetch distance */ \
pfing(rbp -= rsi); \
for(unsigned int loopA = 2*clm; loopA; loopA--) { \
	g2cl_eight_reals_unfft_2(rsi, 64, blkdst, rcx, 64, clmblkdst); \
} \
\
/* Do 46 four_complex_unfft macros */ \
/*	distance between fft data elements is 2K */ \
/*	do 23 macros each processing 16 data values */ \
\
rsi += -2*clm*64+2*blkdst;	/* Next source pointer */ \
rcx += -2*clm*64+2*clmblkdst;/* Next dest pointer */ \
rdi = (uintptr_t)g->sincos1;/* Load sin/cos pointer */ \
/* al = 256/4+6; */ /* 23 iterations of 2*clm */ \
for(unsigned int loopA = 23; loopA; ) { \
	do{ \
		for(unsigned int loopC = 2*clm; loopC; loopC--) { \
			g2cl_four_complex_unfft(rsi, 64, blkdst, rcx, 64, clmblkdst); \
		} \
		rsi += -2*clm*64+2*blkdst;	/* Next source pointer */ \
		rcx += -2*clm*64+2*clmblkdst;/* Next dest pointer */ \
		rdi += XMM_SCD;		/* Next sine/cosine pointer */ \
	}while(--loopA & 3); \
	rcx += -8*clmblkdst+clmblkdst8; \
} \
pfing(rcx = rsi-24*2*blkdst);	/* Create prefetch pointer */ \
rsi = (uintptr_t)g->scratch_area;	/* Restore source pointer */ \
end_timer(16); \
\
/* Do inverse FFT levels 6 */ \
/* On input the 64-byte cache lines hold these data values: */ \
/*	0K	4K	8K	12K */ \
/*	2K	... */ \
/*	16K	... */ \
/*	... */ \
/* On output the 64-byte cache lines hold these data values: */ \
/*	0K	8K	16K	24K */ \
/*	2K	... */ \
/*	... */ \
/*	6K */ \
/*	32K	... */ \
/*	... */ \
/* */ \
/* Do 4 eight_reals_unfft_1 macros */ \
/*	distance between fft data elements is 8K */ \
/*	do 2 macros each processing 16 data values */ \
/* */ \
/* Do 44 two_two_complex_unfft macros */ \
/*	distance between fft data elements is 8K */ \
/*	do 22 macros each processing 16 data values */ \
\
/* Do inverse FFT levels 4,5 */ \
/* On input the 64-byte cache lines hold these data values: */ \
/*	0K	8K	16K	24K */ \
/*	2K	... */ \
/*	... */ \
/*	6K */ \
/*	32K	... */ \
/*	... */ \
/* On output the 64-byte cache lines hold these data values: */ \
/*	0K	32K	64K	96K */ \
/*	2K	... */ \
/*	... */ \
/*	30K */ \
/*	128K	... */ \
/*	... */ \
/* */ \
/* Do 8 eight_reals_unfft_2 macros */ \
/* Do 8 nop_two_two_complex_unfft_2 macros */ \
/*	distance between fft data elements is 16K */ \
/*	do 4 macros each processing 32 data values */ \
/* */ \
/* Do 32 four_complex_unfft macros */ \
/*	distance between fft data elements is 16K */ \
/*	do 8 macros each processing 32 data values */ \
\
/* Do inverse FFT levels 1,2,3 */ \
/* On input the 64-byte cache lines hold these data values: */ \
/*	0K	32K	64K	96K */ \
/*	2K	... */ \
/*	... */ \
/*	30K */ \
/*	128K	... */ \
/*	... */ \
/* On output the 64-byte cache lines hold these data values: */ \
/*	0K	96K	192K	288K */ \
/*	... */ \
/*	94K */ \
/* */ \
/* Do 64 six_reals_last_unfft macros */ \
/*	distance between fft data elements is 64K */ \
/*	do 16 macros each processing 24 data values */ \
\
start_timer(17); \
rdi = (uintptr_t)g->sincos2;/* Load sin/cos pointer */ \
rbx = (uintptr_t)g->sincos3;/* Load sin/cos pointer */ \
/* 2 iterations */ \
xsix_reals_unfft_preload; \
\
for(unsigned int loopA = 2; loopA; loopA--) { \
	xtouchlines(12);			/* Preload TLB */ \
	\
	for(unsigned int loopB = clm; loopB; loopB--) { \
		for(unsigned int loopC = 2; loopC; loopC--) { \
			x2cl_eight_reals_unfft_1(rsi, 64, 2*clmblkdst); \
		} \
		\
		rsi += -2*64+4*clmblkdst;	/* Next source pointer */ \
		/* ah = 256/2+6; */ /* 11 iterations */ \
		for(unsigned int loopC = 11; loopC; ) { \
			do{ \
				for(unsigned int loopE = 2; loopE; loopE--) { \
					x2cl_two_complex_unfft(rsi, 64, 2*clmblkdst); \
				} \
				rsi += -2*64+4*clmblkdst;	/* Next source pointer */ \
				rdi += XMM_SCD;		/* Next sine/cosine pointer */ \
			}while(--loopC & 1); \
			rsi += -2*4*clmblkdst+clmblkdst8;/* Next source pointer */ \
		} \
		rdi += -11*XMM_SCD;	/* Restore sine/cosine pointer */ \
		\
		rsi += -6*clmblkdst8;	/* Restore source pointer */ \
		for(unsigned int loopC = 2; loopC; loopC--) { \
			for(unsigned int loopD = 2; loopD; loopD--) { \
				x4cl_eight_reals_unfft_2(rsi, 64, 4*clmblkdst, clmblkdst8); \
			} \
			rsi += -2*64+2*clmblkdst; \
		} \
		\
		rsi += -2*2*clmblkdst+2*clmblkdst8;/* Next source pointer */ \
		for(unsigned int loopC = 2; loopC; loopC--) { \
			for(unsigned int loopD = 2; loopD; loopD--) { \
				xloop1sp(c3b, 2*clmblkdst, x4cl_four_complex_unfft(rsi, 64, 4*clmblkdst, clmblkdst8, rbx)); \
			} \
			rsi += -2*2*clmblkdst+2*clmblkdst8; \
			rbx += XMM_SCD;		/* Next sine/cosine pointer */ \
		} \
		rbx += -2*XMM_SCD;		/* Restore sine/cosine pointer */ \
		\
		rsi += -3*2*clmblkdst8;	/* Restore source pointer */ \
		for(unsigned int loopC = 2; loopC; loopC--) { \
			for(unsigned int loopD = 4; loopD; loopD--) { \
				xloop1sp(c1b, 2*clmblkdst, x3cl_six_reals_last_unfft(rsi, 64, 2*clmblkdst8)); \
			} \
			rsi += -4*2*clmblkdst+clmblkdst8;/* Next source pointer */ \
		} \
		\
		EXEC248##clm(rsi += -2*clmblkdst8+2*64);	/* Next source pointer */ \
		EXEC248##clm(pfing(rcx += -12*blkdst+128));	/* Next prefetch pointer */ \
	} \
	EXEC248##clm(rsi += -clm*2*64+clmblkdst);/* Next source pointer */ \
	EXEC248##clm(pfing(rcx += -clm*128+12*blkdst));	/* Next prefetch pointer */ \
	EXEC1##clm(rsi += -2*clmblkdst8+clmblkdst);/* Next source pointer */ \
	\
} \
end_timer(17); \
\
/* Normalize these values */ \
\
pass1_normalize(1, yes_fft); \
\
/* Copy data from the scratch area back to the FFT data area */ \
\
copy_scratch_data_to_fft(48, clm); \
goto	endlp;			/* Skip FFT code if flag not set */ \
\
/* Fire up auxiliary threads to do pass 1 forward FFT */ \
\
/* pass1: */	pass1_forward_fft_setup(b0b); \
\
/* Copy data from the FFT data area back to the scratch area */ \
\
b0b:	copy_fft_data_to_scratch(48, clm); \
\
/* Do FFT levels 1,2,3 */ \
/* */ \
/* On input the 64-byte cache lines hold these data values: */ \
/*	0K	96K	192K	288K */ \
/*	... */ \
/*	94K */ \
/* On output the 64-byte cache lines hold these data values: */ \
/*	0K	32K	64K	96K */ \
/*	2K	... */ \
/*	... */ \
/*	30K */ \
/*	128K	192K	256K	320K */ \
/*	130K */ \
/*	... */ \
/*	190K */ \
/* */ \
/* Do 64 six_reals_first_fft macros (FFT levels 1,2,3) */ \
/*	distance between fft data elements is 64K */ \
/*	do 16 iterations each processing 24 data values */ \
\
/* Do FFT levels 4,5 */ \
/* Values 0K-63K is real data, 64K-127K is semi-real data, and */ \
/* 128K-383K is complex data. */ \
/* */ \
/* On input the 64-byte cache lines hold these data values: */ \
/*	0K	32K	64K	96K */ \
/*	2K	... */ \
/*	... */ \
/*	30K */ \
/*	128K	192K	256K	320K */ \
/*	130K */ \
/*	... */ \
/*	190K */ \
/* On output the 64-byte cache lines hold these data values: */ \
/*	0K	16K	32K	48K */ \
/*	2K	... */ \
/*	... */ \
/*	14K	... */ \
/*	64K	... */ \
/*	... */ \
/* */ \
/* Do 8 eight_reals_fft_2 macros */ \
/* Do 8 nop_two_two_complex_fft_2 macros */ \
/*	distance between fft data elements is 16K */ \
/*	do 8 macros each processing 16 data values */ \
/* */ \
/* Do 32 four_complex_fft macros */ \
/*	distance between fft data elements is 16K */ \
/*	do 8 macros each processing 32 data values */ \
\
/* Do FFT levels 6 */ \
/* Values 0K-15K is real data, 16K-31K is semi-real data, 32K-383K is */ \
/* complex data. */ \
/* */ \
/* On input the 64-byte cache lines hold these data values: */ \
/*	0K	16K	32K	48K */ \
/*	2K	... */ \
/*	... */ \
/*	14K	... */ \
/*	64K	... */ \
/*	... */ \
/* On output the 64-byte cache lines hold these data values: */ \
/*	0K	8K	16K	24K */ \
/*	2K	... */ \
/*	... */ \
/*	6K	... */ \
/*	32K	... */ \
/*	... */ \
/* */ \
/* Do 4 eight_reals_fft_1 macros */ \
/* Do 4 two_two_complex_fft macros */ \
/*	distance between fft data elements is 8K */ \
/*	do 4 macros each processing 16 data values */ \
/* */ \
/* Do 40 two_two_complex_fft macros */ \
/*	distance between fft data elements is 8K */ \
/*	do 20 macros each processing 16 data values */ \
\
yes_fft:start_timer(21); \
rsi = (uintptr_t)g->scratch_area;	/* Get address of scratch area */ \
rcx = (uintptr_t)g->data_addr; /* Load source address for prefetching */ \
pfing(rbp = (uintptr_t)g->data_prefetch);/* Set prefetch distance */ \
pfing(rbp -= rcx); \
pfing(rcx += 24*blkdst);		/* Create prefetch pointer */ \
rdi = (uintptr_t)g->sincos3;/* Load sin/cos pointer */ \
rbx = (uintptr_t)g->sincos2;/* Load sin/cos pointer */ \
/* 4 iterations */ \
xsix_reals_fft_preload; \
\
for(unsigned int loopA = 4; loopA; loopA--) { \
	xtouchlines(6);			/* Preload TLBs */ \
	\
	for(unsigned int loopB = clm; loopB; loopB--) { \
		for(unsigned int loopC = 2; loopC; loopC--) { \
			for(unsigned int loopD = 2; loopD; loopD--) { \
				xloop1sp(b1b, 4*clmblkdst, x3cl_six_reals_fft(rsi, 64, 2*clmblkdst8)); \
			} \
			rsi += -2*4*clmblkdst+clmblkdst8;/* Next source pointer */ \
		} \
		\
		rsi += -2*clmblkdst8;	/* Restore source pointer */ \
		for(unsigned int loopC = 2; loopC; loopC--) { \
			for(unsigned int loopD = 2; loopD; loopD--) { \
				x2cl_eight_reals_fft_2(rsi, 64, clmblkdst8); \
			} \
			rsi += -2*64+4*clmblkdst;	/* Next source pointer */ \
		} \
		\
		rsi += -2*4*clmblkdst+2*clmblkdst8;/* Next source pointer */ \
		for(unsigned int loopC = 2; loopC; loopC--) { \
			xloop1sp(b3b, 4*clmblkdst, x4cl_four_complex_fft(rsi, 64, clmblkdst8, 2*clmblkdst8, rdi)); \
		} \
		\
		rsi += -2*4*clmblkdst-2*clmblkdst8;/* Restore source pointer */ \
		for(unsigned int loopC = 2; loopC; loopC--) { \
			x2cl_eight_reals_fft_1(rsi, 64, 4*clmblkdst); \
		} \
		\
		rsi += -2*64+clmblkdst8;	/* Load source pointer */ \
		rbx += XMM_SCD;		/* Next sine/cosine pointer */ \
		/* 5 iterations */ \
		for(unsigned int loopC = 5; loopC; loopC--) { \
			for(unsigned int loopD = 2; loopD; loopD--) { \
				x2cl_two_complex_fft(rsi, 64, 4*clmblkdst, rbx); \
			} \
			rsi += -2*64+clmblkdst8;	/* Next source pointer */ \
			rbx += 2*XMM_SCD;		/* Next sine/cosine pointer */ \
		} \
		rbx += -11*XMM_SCD;	/* Restore sine/cosine pointer */ \
		\
		EXEC248##clm(rsi += -6*clmblkdst8+2*64);	/* Next source pointer */ \
		EXEC248##clm(pfing(rcx += -6*blkdst+128));	/* Next prefetch pointer */ \
	} \
	EXEC248##clm(rsi += -clm*2*64+clmblkdst);/* Next source pointer */ \
	EXEC248##clm(pfing(rcx += -clm*128+6*blkdst));	/* Next prefetch pointer */ \
	EXEC1##clm(rsi += -6*clmblkdst8+clmblkdst);/* Next source pointer */ \
	\
} \
rsi += -4*clmblkdst;	/* Restore source pointer */ \
pfing(rcx += -48*blkdst);		/* Restore prefetch pointer */ \
end_timer(21); \
\
/* Do FFT levels 7,8 */ \
/* Values 0K-7K is real data, 8K-15K is semi-real data, 16K-383K is */ \
/* complex data. */ \
/* */ \
/* On input the 64-byte cache lines hold these data values: */ \
/*	0K	8K	16K	24K */ \
/*	2K	... */ \
/*	... */ \
/*	6K	... */ \
/*	32K	... */ \
/*	... */ \
/* On output the 64-byte cache lines hold these data values: */ \
/*	0K	2K	4K	6K */ \
/*	8K	... */ \
/*	... */ \
\
/* Do 1 eight_reals_fft_2 macros */ \
/* Do 1 nop_two_two_complex_fft_2 macros */ \
/* Do 2 four_complex_fft macros */ \
/*	distance between fft data elements is 2K */ \
/*	do 1 macros each processing 32 data values */ \
\
start_timer(24); \
rdi = (uintptr_t)g->sincos1;/* Load sin/cos pointer */ \
for(unsigned int loopA = 2*clm; loopA; loopA--) { \
	g4cl_eight_reals_fft_2(rsi, 64, clmblkdst, 2*clmblkdst, rcx, 64, blkdst, 2*blkdst); \
} \
\
/* Do 44 four_complex_fft macros */ \
/*	distance between fft data elements is 2K */ \
/*	do 11 macros each processing 32 data values */ \
\
rsi += -2*clm*64+4*clmblkdst;/* Next source pointer */ \
rcx += -2*clm*64+4*blkdst;	/* Next dest pointer */ \
rdi += XMM_SCD;		/* Next sine/cosine pointer */ \
/* eax = 6*256+1; */ /* 11 iterations of 2*clm */ \
for(unsigned int loopA = 11; loopA; ) { \
	do{ \
		for(unsigned int loopC = clm; loopC; loopC--) { \
			for(unsigned int loopD = 2; loopD; loopD--) { \
				g4cl_four_complex_fft(rsi, 64, clmblkdst, 2*clmblkdst, rcx, 64, blkdst, 2*blkdst); \
				L2prefetch128(u8ptr(rcx-64+rbp)); \
				L2prefetch128(u8ptr(rcx-64+blkdst+rbp)); \
				pfing(rbp += -64+2*blkdst); \
			} \
			pfing(rbp += 128-4*blkdst); \
		} \
		rsi += -2*clm*64+4*clmblkdst;/* Next source pointer */ \
		rcx += -2*clm*64+4*blkdst;	/* Next dest pointer */ \
		rdi += 2*XMM_SCD;	/* Next sine/cosine pointer */ \
	}while(--loopA & 1); \
	rsi += -8*clmblkdst+clmblkdst8; \
} \
end_timer(24); \
\
/* Work on next 4 sets of 192 values */ \
\
endlp:	xpass1_get_next_block(pass2, c0b, b0b); /* Get next block and jump */ \
/* done: */	xfft_footer_hg2; \
\
\
\
/* Perform a 384K-element negacyclic FFT.  This is done in two passes.  Pass 1 */ \
/* does 8 levels, pass 2 uses common code to perform the last 11 levels. */ \
\

#define hg_pass1sclevels8complex3(procname, xpass2name, pass2_levels, clm) \
 /* (fftname, pass2_levels, clmarg) */ \
\
/* set_FFT_constants(fftname, 0, 1, pass2_levels, clmarg); */ \
xfft_header_hg2(procname, xpass2name, pass2_levels, clm); \
\
/* Do the last 11 FFT levels */ \
/*	do 48 pass2_11_levels_complex */ \
\
/* pass2: */	start_timer(0); \
xpass2name(g); \
end_timer(0); \
if(g->ffttype == 1)		/* We're done if FFTing only */ \
	return; \
\
/* Do the inverse FFT */ \
\
pass1_inverse_fft_setup(c0b, procname); \
\
/* Do inverse FFT levels 7,8 */ \
/* On input the 64-byte cache lines hold these data values: */ \
/*	0K	2K	4K	6K */ \
/*	8K	10K	12K	14K */ \
/*	... */ \
/* On output the 64-byte cache lines hold these data values: */ \
/*	0K	4K	8K	12K */ \
/*	2K	... */ \
/*	16K	... */ \
/*	... */ \
\
/* Do 48 four_complex_unfft macros */ \
/*	distance between fft data elements is 2K */ \
/*	do 24 macros each processing 16 data values */ \
\
c0b:	start_timer(16); \
rsi = (uintptr_t)g->data_addr;/* Get source address */ \
rcx = (uintptr_t)g->scratch_area;	/* Get address of scratch area */ \
rdi = (uintptr_t)g->sincos1;/* Load sin/cos pointer */ \
/* 6 iterations of 4 iters of 2*clm */ \
pfing(rbp = (uintptr_t)g->data_prefetch);/* Set prefetch distance */ \
pfing(rbp -= rsi); \
for(unsigned int loopA = 6; loopA; loopA--) { \
	for(unsigned int loopB = 4; loopB; loopB--) { \
		for(unsigned int loopC = 2*clm; loopC; loopC--) { \
			g2cl_four_complex_unfft(rsi, 64, blkdst, rcx, 64, clmblkdst); \
		} \
		rsi += -2*clm*64+2*blkdst;	/* Next source pointer */ \
		rcx += -2*clm*64+2*clmblkdst;/* Next dest pointer */ \
		rdi += XMM_SCD;		/* Next sine/cosine pointer */ \
	} \
	rcx += -4*2*clmblkdst+clmblkdst8; \
} \
rcx = rsi-24*2*blkdst;	/* Create prefetch pointer */ \
rsi = (uintptr_t)g->scratch_area;	/* Restore source pointer */ \
end_timer(16); \
\
/* Do inverse FFT levels 5,6 */ \
/* On input the 64-byte cache lines hold these data values: */ \
/*	0K	4K	8K	12K */ \
/*	2K	... */ \
/*	16K	... */ \
/*	... */ \
/* On output the 64-byte cache lines hold these data values: */ \
/*	0K	16K	32K	48K */ \
/*	2K	... */ \
/*	... */ \
/*	14K */ \
/*	64K	... */ \
/*	... */ \
\
/* Do 48 four_complex_unfft macros */ \
/*	distance between fft data elements is 8K */ \
/*	do 12 macros each processing 32 data values */ \
\
start_timer(17); \
rdi = (uintptr_t)g->sincos2;/* Load sin/cos pointer */ \
/* 6 iterations of 2 iters of 2*clm */ \
for(unsigned int loopA = 6; loopA; loopA--) { \
	for(unsigned int loopB = 2; loopB; loopB--) { \
		xtouch(xptr(rcx+(clm-1)*128));	/* Preload TLB */ \
		c5b:	x4cl_four_complex_unfft(rsi, 64, 2*clmblkdst, 4*clmblkdst, rdi); \
		xloop1_clm##clm(c5b, clmblkdst); \
	} \
	rsi += -2*clmblkdst+clmblkdst8;/* Next source pointer */ \
	rdi += XMM_SCD;		/* Next sine/cosine pointer */ \
} \
rsi += -6*clmblkdst8;	/* Restore source pointer */ \
end_timer(17); \
\
/* Do inverse FFT levels 4 */ \
/* On input the 64-byte cache lines hold these data values: */ \
/*	0K	16K	32K	48K */ \
/*	2K	... */ \
/*	... */ \
/*	14K */ \
/*	64K	... */ \
/*	... */ \
/* On output the 64-byte cache lines hold these data values: */ \
/*	0K	32K	64K	96K */ \
/*	2K	... */ \
/*	... */ \
/*	30K */ \
/*	128K	... */ \
/*	... */ \
\
/* Do 48 two_two_complex_unfft macros */ \
/*	distance between fft data elements is 32K */ \
/*	do 24 macros each processing 16 data values */ \
\
start_timer(18); \
rdi = (uintptr_t)g->sincos3;/* Load sin/cos pointer */ \
/* 3 iterations of 2 iters of 4*2*clm */ \
for(unsigned int loopA = 3; loopA; loopA--) { \
	for(unsigned int loopB = 2; loopB; loopB--) { \
		xtouch(xptr(rcx+(clm-1)*128));	/* Preload TLB */ \
		c3b:	x2cl_two_complex_unfft(rsi, 64, clmblkdst8); \
		xloop4_clm##clm(c3b, clmblkdst);		/* 4*2*clm iterations */ \
	} \
	rsi += -2*4*clmblkdst+2*clmblkdst8;/* Next source pointer */ \
	rdi += XMM_SCD;		/* Next sine/cosine pointer */ \
} \
rsi += -3*2*clmblkdst8;	/* Restore source pointer */ \
end_timer(18); \
\
/* Do inverse FFT levels 1,2,3 */ \
/* On input the 64-byte cache lines hold these data values: */ \
/*	0K	32K	64K	96K */ \
/*	2K	... */ \
/*	... */ \
/*	30K */ \
/*	128K	... */ \
/*	... */ \
/* On output the 64-byte cache lines hold these data values: */ \
/*	0K	96K	192K	288K */ \
/*	... */ \
/*	94K */ \
\
/* Do 64 three_complex_last_unfft macros */ \
/*	distance between fft data elements is 64K */ \
/*	do 16 macros each processing 24 data values */ \
\
start_timer(19); \
rdi = (uintptr_t)g->plus1_premults;/* Address of premultiplier table */ \
/* 8 iterations of 2*2*clm */ \
for(unsigned int loopA = 2; loopA; loopA--) { \
	for(unsigned int loopB = 4; loopB; loopB--) { \
		xtouch(xptr(rcx+(clm-1)*128));	/* Preload TLB */ \
		xloop2edi_clm##clm(c1b, clmblkdst, 192, x3cl_three_complex_last_unfft(rsi, 64, 2*clmblkdst8));	/* 2*2*clm iterations */ \
	} \
	rsi += -4*2*clmblkdst+clmblkdst8;/* Next source pointer */ \
} \
end_timer(19); \
\
/* Normalize these values */ \
\
pass1_normalize(1, yes_fft); \
\
/* Copy data from the scratch area back to the FFT data area */ \
\
copy_scratch_data_to_fft(48, clm); \
goto	endlp;			/* Skip FFT code if flag not set */ \
\
/* Fire up auxiliary threads to do pass 1 forward FFT */ \
\
/* pass1: */	pass1_forward_fft_setup(b0b); \
\
/* Copy data from the FFT data area back to the scratch area */ \
\
b0b:	copy_fft_data_to_scratch(48, clm); \
\
/* Do FFT levels 1,2,3 */ \
/* */ \
/* On input the 64-byte cache lines hold these data values: */ \
/*	0K	96K	192K	288K */ \
/*	... */ \
/*	94K */ \
/* On output the 64-byte cache lines hold these data values: */ \
/*	0K	32K	64K	96K */ \
/*	2K	... */ \
/*	... */ \
/*	30K */ \
/*	128K	192K	256K	320K */ \
/*	130K */ \
/*	... */ \
/*	190K */ \
\
/* Do 64 three_complex_first_fft macros (FFT levels 1,2,3) */ \
/*	distance between fft data elements is 64K */ \
/*	do 16 iterations each processing 24 data values */ \
\
yes_fft:start_timer(21); \
rsi = (uintptr_t)g->scratch_area;	/* Get address of scratch area */ \
rdi = (uintptr_t)g->plus1_premults;/* Address of premultiplier table */ \
rcx = (uintptr_t)g->data_addr; /* Load source address for prefetching */ \
pfing(rbp = (uintptr_t)g->data_prefetch);/* Set prefetch distance */ \
pfing(rbp -= rcx); \
pfing(rcx += 26*blkdst);		/* Create prefetching pointer */ \
/* 2 iters of 2 iters of 4*2*clm */ \
for(unsigned int loopA = 2; loopA; loopA--) { \
	for(unsigned int loopB = 2; loopB; loopB--) { \
		xtouch(xptr(rcx+(clm-1)*128));	/* Preload TLB */ \
		b1b:	x3cl_three_complex_fft(rsi, 64, 2*clmblkdst8); \
		xloop4edi_clm##clm(b1b, clmblkdst, 192);	/* 2 * 2*clm iterations */ \
	} \
	rsi += -2*4*clmblkdst+clmblkdst8;/* Next source pointer */ \
} \
rsi += -2*clmblkdst8;	/* Restore source pointer */ \
end_timer(21); \
\
/* Do FFT levels 4 */ \
/* */ \
/* On input the 64-byte cache lines hold these data values: */ \
/*	0K	32K	64K	96K */ \
/*	2K	... */ \
/*	... */ \
/*	30K */ \
/*	128K	192K	256K	320K */ \
/*	130K */ \
/*	... */ \
/*	190K */ \
/* On output the 64-byte cache lines hold these data values: */ \
/*	0K	32K	64K	96K */ \
/*	2K	... */ \
/*	... */ \
/*	30K	... */ \
/*	128K	... */ \
/*	... */ \
\
/* Do 48 two_two_complex_fft macros */ \
/*	distance between fft data elements is 32K */ \
/*	do 24 macros each processing 16 data values */ \
\
start_timer(22); \
rdi = (uintptr_t)g->sincos3;/* Load sin/cos pointer */ \
/* 2 iterations of 4*2*clm */ \
for(unsigned int loopA = 2; loopA; loopA--) { \
	xtouch(xptr(rcx+(clm-1)*128));	/* Preload TLB */ \
	b2b:	x2cl_two_complex_fft_in_place(rsi, 64, clmblkdst8, rdi); \
	xloop4_clm##clm(b2b, clmblkdst);		/* 4 * 2*clm iterations */ \
} \
rsi += -2*4*clmblkdst+2*clmblkdst8;/* Next source pointer */ \
rdi += XMM_SCD;		/* Next sin/cos pointer */ \
/* 2 iters of 2 iters of 4*2*clm */ \
for(unsigned int loopA = 2; loopA; loopA--) { \
	for(unsigned int loopB = 2; loopB; loopB--) { \
		xtouch(xptr(rcx+(clm-1)*128));	/* Preload TLB */ \
		b3b:	x2cl_two_complex_fft(rsi, 64, 2*clmblkdst8, rdi); \
		xloop4_clm##clm(b3b, clmblkdst);		/* 4 * 2*clm iterations */ \
	} \
	rsi += -2*4*clmblkdst+clmblkdst8; \
} \
rsi += -2*clmblkdst8-2*clmblkdst8;/* Restore source pointer */ \
end_timer(22); \
\
/* Do FFT levels 5,6 */ \
/* */ \
/* On input the 64-byte cache lines hold these data values: */ \
/*	0K	32K	64K	96K */ \
/*	2K	... */ \
/*	... */ \
/*	30K	... */ \
/*	128K	... */ \
/*	... */ \
/* On output the 64-byte cache lines hold these data values: */ \
/*	0K	8K	16K	24K */ \
/*	2K	... */ \
/*	... */ \
/*	6K	... */ \
/*	32K	... */ \
/*	... */ \
\
/* Do 48 four_complex_fft macros */ \
/*	distance between fft data elements is 8K */ \
/*	do 12 macros each processing 32 data values */ \
\
start_timer(23); \
rdi = (uintptr_t)g->sincos2;/* Load sin/cos pointer */ \
/* 3 iterations of 4 iters of 2*clm */ \
for(unsigned int loopA = 3; loopA; loopA--) { \
	for(unsigned int loopB = 4; loopB; loopB--) { \
		xtouch(xptr(rcx+(clm-1)*128));	/* Preload TLB */ \
		b5b:	x4cl_four_complex_fft(rsi, 64, 4*clmblkdst, clmblkdst8, rdi); \
		xloop1_clm##clm(b5b, clmblkdst);		/* 2 * 2*clm iterations */ \
	} \
	rsi += -4*clmblkdst+2*clmblkdst8;/* Next source pointer */ \
	rdi += 2*XMM_SCD;		/* Next sine/cosine pointer */ \
} \
rsi += -3*2*clmblkdst8;	/* Restore source pointer */ \
pfing(rcx += -48*blkdst);		/* Restore prefetch pointer */ \
end_timer(23); \
\
/* Do FFT levels 7,8 */ \
/* */ \
/* On input the 64-byte cache lines hold these data values: */ \
/*	0K	8K	16K	24K */ \
/*	2K	... */ \
/*	... */ \
/*	6K	... */ \
/*	32K	... */ \
/*	... */ \
/* On output the 64-byte cache lines hold these data values: */ \
/*	0K	2K	4K	6K */ \
/*	8K	... */ \
/*	... */ \
\
/* Do 48 four_complex_fft macros */ \
/*	distance between fft data elements is 2K */ \
/*	do 12 macros each processing 32 data values */ \
\
start_timer(24); \
rdi = (uintptr_t)g->sincos1;/* Load sin/cos pointer */ \
/* 6 iterations of 2 iters of 2*clm */ \
for(unsigned int loopA = 6; loopA; loopA--) { \
	for(unsigned int loopB = 2; loopB; loopB--) { \
		for(unsigned int loopC = clm; loopC; loopC--) { \
			for(unsigned int loopD = 2; loopD; loopD--) { \
				g4cl_four_complex_fft(rsi, 64, clmblkdst, 2*clmblkdst, rcx, 64, blkdst, 2*blkdst); \
				L2prefetch128(u8ptr(rcx-64+rbp)); \
				L2prefetch128(u8ptr(rcx-64+blkdst+rbp)); \
				pfing(rbp += -64+2*blkdst); \
			} \
			pfing(rbp += 128-4*blkdst); \
		} \
		rsi += -2*clm*64+4*clmblkdst;/* Next source pointer */ \
		rcx += -2*clm*64+4*blkdst;	/* Next dest pointer */ \
		rdi += 2*XMM_SCD;		/* Next sine/cosine pointer */ \
	} \
	rsi += -8*clmblkdst+clmblkdst8; \
} \
end_timer(24); \
\
/* Work on next 4 sets of 192 values */ \
\
endlp:	xpass1_get_next_block(pass2, c0b, b0b); /* Get next block and jump */ \
/* done: */	xfft_footer_hg2; \
\
\
\
/* Perform a 448K-element FFT.  This is done in two passes.  Pass 1 */ \
/* does 8 levels, pass 2 uses common code to perform the last 11 levels. */ \
\

#define hg_pass1sclevels8pfa7(procname, xpass2name, pass2_levels, clm) \
 /* (fftname, pass2_levels, clmarg) */ \
\
/* set_FFT_constants(fftname, 0, 0, pass2_levels, clmarg); */ \
xfft_header_hg2(procname, xpass2name, pass2_levels, clm); \
\
/* Do the last 11 FFT levels */ \
/*	do 1 pass2_11_levels_real and 55 pass2_11_levels_complex */ \
\
/* pass2: */	start_timer(0); \
xpass2name(g); \
end_timer(0); \
if(g->ffttype == 1)		/* We're done if FFTing only */ \
	return; \
\
/* Do the inverse FFT */ \
\
pass1_inverse_fft_setup(c0b, procname); \
\
/* Do inverse FFT levels 7,8 */ \
/* On input the 64-byte cache lines hold these data values: */ \
/*	0K	2K	4K	6K */ \
/*	8K	10K	12K	14K */ \
/*	... */ \
/* On output the 64-byte cache lines hold these data values: */ \
/*	0K	4K	8K	12K */ \
/*	2K	... */ \
/*	16K	... */ \
/*	... */ \
\
/* Do 1 eight_reals_unfft_2 macros */ \
/* Do 1 nop_two_two_complex_unfft_2 macros */ \
/*	distance between fft data elements is 2K */ \
/*	do 1 macros each processing 16 data values */ \
\
c0b:	start_timer(16); \
rsi = (uintptr_t)g->data_addr;/* Get source address */ \
rcx = (uintptr_t)g->scratch_area;	/* Get address of scratch area */ \
pfing(rbp = (uintptr_t)g->data_prefetch);/* Set prefetch distance */ \
pfing(rbp -= rsi); \
for(unsigned int loopA = 2*clm; loopA; loopA--) { \
	g2cl_eight_reals_unfft_2(rsi, 64, blkdst, rcx, 64, clmblkdst); \
} \
\
/* Do 54 four_complex_unfft macros */ \
/*	distance between fft data elements is 2K */ \
/*	do 27 macros each processing 16 data values */ \
\
rsi += -2*clm*64+2*blkdst;	/* Next source pointer */ \
rcx += -2*clm*64+2*clmblkdst;/* Next dest pointer */ \
rdi = (uintptr_t)g->sincos1;/* Load sin/cos pointer */ \
/* al = 256/4+7; */ /* 27 iterations of 2*clm */ \
for(unsigned int loopA = 27; loopA; ) { \
	do{ \
		for(unsigned int loopC = 2*clm; loopC; loopC--) { \
			g2cl_four_complex_unfft(rsi, 64, blkdst, rcx, 64, clmblkdst); \
		} \
		rsi += -2*clm*64+2*blkdst;	/* Next source pointer */ \
		rcx += -2*clm*64+2*clmblkdst;/* Next dest pointer */ \
		rdi += XMM_SCD;		/* Next sine/cosine pointer */ \
	}while(--loopA & 3); \
	rcx += -8*clmblkdst+clmblkdst8; \
} \
rcx = rsi-28*2*blkdst;	/* Create prefetch pointer */ \
rsi = (uintptr_t)g->scratch_area;	/* Restore source pointer */ \
end_timer(16); \
\
/* Do inverse FFT levels 6 */ \
/* On input the 64-byte cache lines hold these data values: */ \
/*	0K	4K	8K	12K */ \
/*	2K	... */ \
/*	16K	... */ \
/*	... */ \
/* On output the 64-byte cache lines hold these data values: */ \
/*	0K	8K	16K	24K */ \
/*	2K	... */ \
/*	... */ \
/*	6K */ \
/*	32K	... */ \
/*	... */ \
/* */ \
/* Do 4 eight_reals_unfft_1 macros */ \
/*	distance between fft data elements is 8K */ \
/*	do 2 macros each processing 16 data values */ \
/* */ \
/* Do 52 two_two_complex_unfft macros */ \
/*	distance between fft data elements is 8K */ \
/*	do 26 macros each processing 16 data values */ \
\
/* Do inverse FFT levels 4,5 */ \
/* On input the 64-byte cache lines hold these data values: */ \
/*	0K	8K	16K	24K */ \
/*	2K	... */ \
/*	... */ \
/*	6K */ \
/*	32K	... */ \
/*	... */ \
/* On output the 64-byte cache lines hold these data values: */ \
/*	0K	16K	32K	48K */ \
/*	2K	... */ \
/*	... */ \
/*	14K */ \
/*	64K	96K	128K	160K */ \
/*	... */ \
/* */ \
/* Do 8 eight_reals_unfft_2 macros */ \
/*	distance between fft data elements is 16K */ \
/*	do 4 macros each processing 16 data values */ \
/* */ \
/* Do 48 four_complex_unfft macros */ \
/*	distance between fft data elements is 16K */ \
/*	do 12 macros each processing 32 data values */ \
\
/* Do inverse FFT levels 1,2,3 */ \
/* On input the 64-byte cache lines hold these data values: */ \
/*	0K	16K	32K	48K */ \
/*	2K	... */ \
/*	... */ \
/*	14K */ \
/*	64K	96K	128K	160K */ \
/*	... */ \
/* On output the 64-byte cache lines hold these data values: */ \
/*	0K	16K	224K	240K */ \
/*	2K	... */ \
/*	... */ \
/*	14K */ \
/*	32K	64K	256K	288K */ \
/*	34K */ \
/*	... */ \
/*	62K */ \
/*	96K	160K	320K	384K */ \
/*	98K */ \
/*	... */ \
/*	158K */ \
/* */ \
/* Do 64 seven_reals_last_unfft macros */ \
/*	distance between fft data elements is 64K */ \
/*	do 8 macros each processing 56 data values */ \
\
start_timer(17); \
rdi = (uintptr_t)g->sincos2;/* Load sin/cos pointer */ \
rbx = (uintptr_t)g->sincos3;/* Load sin/cos pointer */ \
/* 2 iterations */ \
xseven_reals_unfft_preload; \
\
for(unsigned int loopA = 2; loopA; loopA--) { \
	xtouchlines(14);			/* Preload TLBs */ \
	\
	for(unsigned int loopB = clm; loopB; loopB--) { \
		for(unsigned int loopC = 2; loopC; loopC--) { \
			x2cl_eight_reals_unfft_1(rsi, 64, 2*clmblkdst); \
		} \
		\
		rsi += -2*64+4*clmblkdst;	/* Next source pointer */ \
		/* ah = 256/2+7; */ /* 13 iterations */ \
		for(unsigned int loopC = 13; loopC; ) { \
			do{ \
				for(unsigned int loopE = 2; loopE; loopE--) { \
					x2cl_two_complex_unfft(rsi, 64, 2*clmblkdst); \
				} \
				rsi += -2*64+4*clmblkdst;	/* Next source pointer */ \
				rdi += XMM_SCD;		/* Next sine/cosine pointer */ \
			}while(--loopC & 1); \
			rsi += -2*4*clmblkdst+clmblkdst8;/* Next source pointer */ \
		} \
		rdi += -13*XMM_SCD;	/* Restore sine/cosine pointer */ \
		\
		rsi += -7*clmblkdst8;	/* Restore source pointer */ \
		for(unsigned int loopC = 2; loopC; loopC--) { \
			for(unsigned int loopD = 2; loopD; loopD--) { \
				x2cl_half_eight_reals_unfft_2(rsi, 64, 4*clmblkdst); \
			} \
			rsi += -2*64+2*clmblkdst; \
		} \
		\
		rsi += -2*2*clmblkdst+clmblkdst8;/* Next source pointer */ \
		/* 3 iterations */ \
		for(unsigned int loopC = 3; loopC; loopC--) { \
			for(unsigned int loopD = 2; loopD; loopD--) { \
				xloop1sp(c3b, 2*clmblkdst, x4cl_four_complex_unfft(rsi, 64, 4*clmblkdst, clmblkdst8, rbx)); \
			} \
			rsi += -2*2*clmblkdst+2*clmblkdst8; \
			rbx += XMM_SCD;		/* Next sine/cosine pointer */ \
		} \
		rbx += -3*XMM_SCD;		/* Restore sine/cosine pointer */ \
		\
		rsi += -3*2*clmblkdst8-clmblkdst8;/* Restore source pointer */ \
		for(unsigned int loopC = 4; loopC; loopC--) { \
			xloop1prefetch2sp(c1b, 2*clmblkdst, x7cl_seven_reals_last_unfft(rsi, 64, clmblkdst8)); \
		} \
		\
		EXEC248##clm(rsi += -4*2*clmblkdst+2*64);/* Next source pointer */ \
		EXEC248##clm(pfing(rcx += -14*blkdst+128));	/* Next prefetch pointer */ \
	} \
	EXEC248##clm(rsi += -clm*2*64+clmblkdst);/* Next source pointer */ \
	EXEC248##clm(pfing(rcx += -clm*128+14*blkdst));	/* Next prefetch pointer */ \
	EXEC1##clm(rsi += -4*2*clmblkdst+clmblkdst);/* Next source pointer */ \
	\
} \
end_timer(17); \
\
/* Normalize these values */ \
\
pass1_normalize(1, yes_fft); \
\
/* Copy data from the scratch area back to the FFT data area */ \
\
copy_scratch_data_to_fft(56, clm); \
goto	endlp;			/* Skip FFT code if flag not set */ \
\
/* Fire up auxiliary threads to do pass 1 forward FFT */ \
\
/* pass1: */	pass1_forward_fft_setup(b0b); \
\
/* Copy data from the FFT data area back to the scratch area */ \
\
b0b:	copy_fft_data_to_scratch(56, clm); \
\
/* Do FFT levels 1,2,3 */ \
/* */ \
/* On input the 64-byte cache lines hold these data values: */ \
/*	0K	16K	224K	240K */ \
/*	2K	... */ \
/*	... */ \
/*	14K */ \
/*	32K	64K	256K	288K */ \
/*	34K */ \
/*	... */ \
/*	62K */ \
/*	96K	160K	320K	384K */ \
/*	98K */ \
/*	... */ \
/*	158K */ \
/* On output the 64-byte cache lines hold these data values: */ \
/*	0K	16K	32K	48K */ \
/*	2K	... */ \
/*	... */ \
/*	14K */ \
/*	64K	96K	128K	160K */ \
/*	66K */ \
/*	... */ \
/*	94K */ \
/*	192K	256K	320K	384K */ \
/*	194K */ \
/*	... */ \
/*	254K */ \
/* */ \
/* Do 64 seven_reals_first_fft macros (FFT levels 1,2,3) */ \
/*	distance between fft data elements is 64K */ \
/*	do 8 iterations each processing 56 data values */ \
\
/* Do FFT levels 4,5 */ \
/* Values 0K-63K is real data, 64K-447K is complex data. */ \
/* */ \
/* On input the 64-byte cache lines hold these data values: */ \
/*	0K	16K	32K	48K */ \
/*	2K	... */ \
/*	... */ \
/*	14K */ \
/*	64K	96K	128K	160K */ \
/*	66K */ \
/*	... */ \
/*	94K */ \
/*	192K	256K	320K	384K */ \
/*	194K */ \
/*	... */ \
/*	254K */ \
/* On output the 64-byte cache lines hold these data values: */ \
/*	0K	16K	32K	48K */ \
/*	2K	... */ \
/*	... */ \
/*	14K	... */ \
/*	64K	... */ \
/*	... */ \
/* */ \
/* Do 8 eight_reals_fft_2 macros */ \
/*	distance between fft data elements is 16K */ \
/*	do 8 macros each processing 8 data values */ \
/* */ \
/* Do 16 four_complex_fft macros */ \
/*	distance between fft data elements is 16K */ \
/*	do 8 macros each processing 16 data values */ \
/* */ \
/* Do 32 four_complex_fft macros */ \
/*	distance between fft data elements is 16K */ \
/*	do 8 macros each processing 32 data values */ \
\
/* Do FFT levels 6 */ \
/* Values 0K-15K is real data, 16K-31K is semi-real data, 32K-447K is */ \
/* complex data. */ \
/* */ \
/* On input the 64-byte cache lines hold these data values: */ \
/*	0K	16K	32K	48K */ \
/*	2K	... */ \
/*	... */ \
/*	14K	... */ \
/*	64K	... */ \
/*	... */ \
/* On output the 64-byte cache lines hold these data values: */ \
/*	0K	8K	16K	24K */ \
/*	2K	... */ \
/*	... */ \
/*	6K	... */ \
/*	32K	... */ \
/*	... */ \
\
/* Do 4 eight_reals_fft_1 macros */ \
/* Do 4 two_two_complex_fft macros */ \
/*	distance between fft data elements is 8K */ \
/*	do 4 macros each processing 16 data values */ \
/* */ \
/* Do 48 two_two_complex_fft macros */ \
/*	distance between fft data elements is 8K */ \
/*	do 24 macros each processing 16 data values */ \
\
yes_fft:start_timer(21); \
rsi = (uintptr_t)g->scratch_area;	/* Get address of scratch area */ \
rcx = (uintptr_t)g->data_addr; /* Load source address for prefetching */ \
pfing(rbp = (uintptr_t)g->data_prefetch);/* Set prefetch distance */ \
pfing(rbp -= rcx); \
pfing(rcx += +28*blkdst);		/* Create prefetching pointer */ \
rdi = (uintptr_t)g->sincos3;/* Load sin/cos pointer */ \
rbx = (uintptr_t)g->sincos2;/* Load sin/cos pointer */ \
/* 4 iterations */ \
xseven_reals_fft_preload; \
\
for(unsigned int loopA = 4; loopA; loopA--) { \
	xtouchlines(7);			/* Preload TLBs */ \
	\
	for(unsigned int loopB = clm; loopB; loopB--) { \
		for(unsigned int loopC = 2; loopC; loopC--) { \
			xloop1prefetch2sp(b1b, 4*clmblkdst, x7cl_seven_reals_fft(rsi, 64, clmblkdst8)); \
		} \
		\
		rsi += -2*4*clmblkdst;	/* Restore source pointer */ \
		for(unsigned int loopC = 2; loopC; loopC--) { \
			for(unsigned int loopD = 2; loopD; loopD--) { \
				x1cl_half_eight_reals_fft_2(rsi, 64); \
			} \
			rsi += -2*64+4*clmblkdst;	/* Next source pointer */ \
		} \
		\
		rsi += -2*4*clmblkdst+clmblkdst8;/* Next source pointer */ \
		for(unsigned int loopC = 2; loopC; loopC--) { \
			for(unsigned int loopD = 2; loopD; loopD--) { \
				x2cl_four_complex_fft(rsi, 64, clmblkdst8); \
			} \
			rsi += -2*64+4*clmblkdst;	/* Next source pointer */ \
		} \
		L2prefetch128(u8ptr(rcx+rbp)); \
		pfing(rcx += blkdst);		/* Next prefetch pointer */ \
		\
		rsi += -2*4*clmblkdst+2*clmblkdst8;/* Next source pointer */ \
		rdi += XMM_SCD;		/* Next sine/cosine pointer */ \
		for(unsigned int loopC = 2; loopC; loopC--) { \
			xloop1sp(b3b, 4*clmblkdst, x4cl_four_complex_fft(rsi, 64, clmblkdst8, 2*clmblkdst8, rdi)); \
		} \
		rdi += -XMM_SCD;		/* Restore sine/cosine pointer */ \
		\
		rsi += -2*4*clmblkdst-3*clmblkdst8;/* Restore source pointer */ \
		for(unsigned int loopC = 2; loopC; loopC--) { \
			x2cl_eight_reals_fft_1(rsi, 64, 4*clmblkdst); \
		} \
		\
		rsi += -2*64+clmblkdst8;	/* Next source pointer */ \
		rbx += XMM_SCD;		/* Next sine/cosine pointer */ \
		/* 6 iterations */ \
		for(unsigned int loopC = 6; loopC; loopC--) { \
			for(unsigned int loopD = 2; loopD; loopD--) { \
				x2cl_two_complex_fft(rsi, 64, 4*clmblkdst, rbx); \
			} \
			rsi += -2*64+clmblkdst8;	/* Next source pointer */ \
			rbx += 2*XMM_SCD;		/* Next sine/cosine pointer */ \
		} \
		rbx += -13*XMM_SCD;	/* Restore sine/cosine pointer */ \
		\
		EXEC248##clm(rsi += -7*clmblkdst8+2*64);	/* Next source pointer */ \
		EXEC248##clm(pfing(rcx += -7*blkdst+128));	/* Next prefetch pointer */ \
	} \
	EXEC248##clm(rsi += -clm*2*64+clmblkdst);/* Next source pointer */ \
	EXEC248##clm(pfing(rcx += -clm*128+7*blkdst));	/* Next prefetch pointer */ \
	EXEC1##clm(rsi += -7*clmblkdst8+clmblkdst);/* Next source pointer */ \
	\
} \
rsi += -4*clmblkdst;	/* Restore source pointer */ \
pfing(rcx += -56*blkdst);		/* Restore prefetch pointer */ \
end_timer(21); \
\
/* Do FFT levels 7,8 */ \
/* Values 0K-7K is real data, 8K-15K is semi-real data, 16K-447K is */ \
/* complex data. */ \
/* */ \
/* On input the 64-byte cache lines hold these data values: */ \
/*	0K	8K	16K	24K */ \
/*	2K	... */ \
/*	... */ \
/*	6K	... */ \
/*	32K	... */ \
/*	... */ \
/* On output the 64-byte cache lines hold these data values: */ \
/*	0K	2K	4K	6K */ \
/*	8K	... */ \
/*	... */ \
\
/* Do 1 eight_reals_fft_2 macros */ \
/* Do 1 nop_two_two_complex_fft_2 macros */ \
/* Do 2 four_complex_fft macros */ \
/*	distance between fft data elements is 2K */ \
/*	do 1 macros each processing 32 data values */ \
\
start_timer(24); \
rdi = (uintptr_t)g->sincos1;/* Load sin/cos pointer */ \
for(unsigned int loopA = 2*clm; loopA; loopA--) { \
	g4cl_eight_reals_fft_2(rsi, 64, clmblkdst, 2*clmblkdst, rcx, 64, blkdst, 2*blkdst); \
} \
\
/* Do 52 four_complex_fft macros */ \
/*	distance between fft data elements is 2K */ \
/*	do 13 macros each processing 32 data values */ \
\
rsi += -2*clm*64+4*clmblkdst;/* Next source pointer */ \
rcx += -2*clm*64+4*blkdst;	/* Next dest pointer */ \
rdi += XMM_SCD;		/* Next sine/cosine pointer */ \
/* eax = 7*256+1; */ /* 13 iterations of 2*clm */ \
for(unsigned int loopA = 13; loopA; ) { \
	do{ \
		for(unsigned int loopC = clm; loopC; loopC--) { \
			for(unsigned int loopD = 2; loopD; loopD--) { \
				g4cl_four_complex_fft(rsi, 64, clmblkdst, 2*clmblkdst, rcx, 64, blkdst, 2*blkdst); \
				L2prefetch128(u8ptr(rcx-64+rbp)); \
				L2prefetch128(u8ptr(rcx-64+blkdst+rbp)); \
				pfing(rbp += -64+2*blkdst); \
			} \
			pfing(rbp += 128-4*blkdst); \
		} \
		rsi += -2*clm*64+4*clmblkdst;/* Next source pointer */ \
		rcx += -2*clm*64+4*blkdst;	/* Next dest pointer */ \
		rdi += 2*XMM_SCD;		/* Next sine/cosine pointer */ \
	}while(--loopA & 1); \
	rsi += -8*clmblkdst+clmblkdst8; \
} \
end_timer(24); \
\
/* Work on next 4 sets of 224 values */ \
\
endlp:	xpass1_get_next_block(pass2, c0b, b0b); /* Get next block and jump */ \
/* done: */	xfft_footer_hg2; \
\
\
\
/* Perform a 512K-element FFT.  This is done in two passes.  Pass 1 */ \
/* does 8 levels, pass 2 uses common code to perform the last 11 levels. */ \
\

#define hg_pass1sclevels8pfa8(procname, xpass2name, pass2_levels, clm) \
 /* (fftname, pass2_levels, clmarg) */ \
\
/* set_FFT_constants(fftname, 0, 0, pass2_levels, clmarg); */ \
xfft_header_hg2(procname, xpass2name, pass2_levels, clm); \
\
/* Do the last 11 FFT levels */ \
/*	do 1 pass2_11_levels_real and 63 pass2_11_levels_complex */ \
\
/* pass2: */	start_timer(0); \
xpass2name(g); \
end_timer(0); \
if(g->ffttype == 1)		/* We're done if FFTing only */ \
	return; \
\
/* Do the inverse FFT */ \
\
pass1_inverse_fft_setup(c0b, procname); \
\
/* Do inverse FFT levels 7,8 */ \
/* On input the 64-byte cache lines hold these data values: */ \
/*	0K	2K	4K	6K */ \
/*	8K	10K	12K	14K */ \
/*	... */ \
/* On output the 64-byte cache lines hold these data values: */ \
/*	0K	4K	8K	12K */ \
/*	2K	... */ \
/*	16K	... */ \
/*	... */ \
\
/* Do 1 eight_reals_unfft_2 macros */ \
/* Do 1 nop_two_two_complex_unfft_2 macros */ \
/*	distance between fft data elements is 2K */ \
/*	do 1 macros each processing 16 data values */ \
\
c0b:	start_timer(16); \
rsi = (uintptr_t)g->data_addr;/* Get source address */ \
rcx = (uintptr_t)g->scratch_area;	/* Get address of scratch area */ \
pfing(rbp = (uintptr_t)g->data_prefetch);/* Set prefetch distance */ \
pfing(rbp -= rsi); \
for(unsigned int loopA = 2*clm; loopA; loopA--) { \
	g2cl_eight_reals_unfft_2(rsi, 64, blkdst, rcx, 64, clmblkdst); \
} \
\
/* Do 62 four_complex_unfft macros */ \
/*	distance between fft data elements is 2K */ \
/*	do 31 macros each processing 16 data values */ \
\
rsi += -2*clm*64+2*blkdst;	/* Next source pointer */ \
rcx += -2*clm*64+2*clmblkdst;/* Next dest pointer */ \
rdi = (uintptr_t)g->sincos1;/* Load sin/cos pointer */ \
/* al = 256/4+8; */ /* 31 iterations of 2*clm */ \
for(unsigned int loopA = 31; loopA; ) { \
	do{ \
		for(unsigned int loopC = 2*clm; loopC; loopC--) { \
			g2cl_four_complex_unfft(rsi, 64, blkdst, rcx, 64, clmblkdst); \
		} \
		rsi += -2*clm*64+2*blkdst;	/* Next source pointer */ \
		rcx += -2*clm*64+2*clmblkdst;/* Next dest pointer */ \
		rdi += XMM_SCD;		/* Next sine/cosine pointer */ \
	}while(--loopA & 3); \
	rcx += -8*clmblkdst+clmblkdst8; \
} \
rcx = rsi-32*2*blkdst;	/* Create prefetch pointer */ \
rsi = (uintptr_t)g->scratch_area;	/* Restore source pointer */ \
end_timer(16); \
\
/* Do inverse FFT levels 6 */ \
/* On input the 64-byte cache lines hold these data values: */ \
/*	0K	4K	8K	12K */ \
/*	2K	... */ \
/*	16K	... */ \
/*	... */ \
/* On output the 64-byte cache lines hold these data values: */ \
/*	0K	8K	16K	24K */ \
/*	2K	... */ \
/*	... */ \
/*	6K */ \
/*	32K	... */ \
/*	... */ \
/* */ \
/* Do 4 eight_reals_unfft_1 macros */ \
/*	distance between fft data elements is 8K */ \
/*	do 2 macros each processing 16 data values */ \
/* */ \
/* Do 60 two_two_complex_unfft macros */ \
/*	distance between fft data elements is 8K */ \
/*	do 30 macros each processing 16 data values */ \
\
/* Do inverse FFT levels 4,5 */ \
/* On input the 64-byte cache lines hold these data values: */ \
/*	0K	8K	16K	24K */ \
/*	2K	... */ \
/*	... */ \
/*	6K */ \
/*	32K	... */ \
/*	... */ \
/* On output the 64-byte cache lines hold these data values: */ \
/*	0K	32K	64K	96K */ \
/*	2K	... */ \
/*	... */ \
/*	30K */ \
/*	128K	... */ \
/*	... */ \
/* */ \
/* Do 8 eight_reals_unfft_2 macros */ \
/* Do 8 nop_two_two_complex_unfft_2 macros */ \
/*	distance between fft data elements is 16K */ \
/*	do 4 macros each processing 32 data values */ \
/* */ \
/* Do 48 four_complex_unfft macros */ \
/*	distance between fft data elements is 16K */ \
/*	do 12 macros each processing 32 data values */ \
\
/* Do inverse FFT levels 1,2,3 */ \
/* On input the 64-byte cache lines hold these data values: */ \
/*	0K	32K	64K	96K */ \
/*	2K	... */ \
/*	... */ \
/*	30K */ \
/*	128K	... */ \
/*	... */ \
/* On output the 64-byte cache lines hold these data values: */ \
/*	0K	128K	256K	384K */ \
/*	2K	... */ \
/*	... */ \
/*	126K */ \
/* */ \
/* Do 64 eight_reals_last_unfft macros */ \
/*	distance between fft data elements is 64K */ \
/*	do 16 macros each processing 32 data values */ \
\
start_timer(17); \
rdi = (uintptr_t)g->sincos2;/* Load sin/cos pointer */ \
rbx = (uintptr_t)g->sincos3;/* Load sin/cos pointer */ \
/* 2 iterations */ \
\
for(unsigned int loopA = 2; loopA; loopA--) { \
	xtouchlines(16);			/* Preload TLBs */ \
	\
	for(unsigned int loopB = clm; loopB; loopB--) { \
		for(unsigned int loopC = 2; loopC; loopC--) { \
			x2cl_eight_reals_unfft_1(rsi, 64, 2*clmblkdst); \
		} \
		\
		rsi += -2*64+4*clmblkdst;	/* Next source pointer */ \
		/* ah = 256/2+8; */ /* 15 iterations */ \
		for(unsigned int loopC = 15; loopC; ) { \
			do{ \
				for(unsigned int loopE = 2; loopE; loopE--) { \
					x2cl_two_complex_unfft(rsi, 64, 2*clmblkdst); \
				} \
				rsi += -2*64+4*clmblkdst;	/* Next source pointer */ \
				rdi += XMM_SCD;		/* Next sine/cosine pointer */ \
			}while(--loopC & 1); \
			rsi += -2*4*clmblkdst+clmblkdst8;/* Next source pointer */ \
		} \
		rdi += -15*XMM_SCD;	/* Restore sine/cosine pointer */ \
		\
		rsi += -8*clmblkdst8;	/* Restore source pointer */ \
		for(unsigned int loopC = 2; loopC; loopC--) { \
			xloop1sp(c2b, 2*clmblkdst, x4cl_eight_reals_unfft_2(rsi, 64, 4*clmblkdst, clmblkdst8)); \
		} \
		\
		rsi += -2*2*clmblkdst+2*clmblkdst8;/* Next source pointer */ \
		/* 3 iterations */ \
		for(unsigned int loopC = 3; loopC; loopC--) { \
			for(unsigned int loopD = 2; loopD; loopD--) { \
				xloop1sp(c3b, 2*clmblkdst, x4cl_four_complex_unfft(rsi, 64, 4*clmblkdst, clmblkdst8, rbx)); \
			} \
			rsi += -2*2*clmblkdst+2*clmblkdst8; \
			rbx += XMM_SCD;		/* Next sine/cosine pointer */ \
		} \
		rbx += -3*XMM_SCD;		/* Next sine/cosine pointer */ \
		\
		rsi += -4*2*clmblkdst8;	/* Restore source pointer */ \
		for(unsigned int loopC = 2; loopC; loopC--) { \
			for(unsigned int loopD = 4; loopD; loopD--) { \
				xloop1sp(c1b, 2*clmblkdst, x4cl_eight_reals_last_unfft(rsi, 64, 2*clmblkdst8, 4*clmblkdst8)); \
			} \
			rsi += -4*2*clmblkdst+clmblkdst8;/* Next source pointer */ \
		} \
		\
		EXEC248##clm(rsi += -2*clmblkdst8+2*64);	/* Next source pointer */ \
		EXEC248##clm(pfing(rcx += -16*blkdst+128));	/* Next prefetch pointer */ \
	} \
	EXEC248##clm(rsi += -clm*2*64+clmblkdst);/* Next source pointer */ \
	EXEC248##clm(pfing(rcx += -clm*128+16*blkdst));	/* Next prefetch pointer */ \
	EXEC1##clm(rsi += -2*clmblkdst8+clmblkdst);/* Next source pointer */ \
	\
} \
end_timer(17); \
\
/* Normalize these values */ \
\
pass1_normalize(1, yes_fft); \
\
/* Copy data from the scratch area back to the FFT data area */ \
\
copy_scratch_data_to_fft(64, clm); \
goto	endlp;			/* Skip FFT code if flag not set */ \
\
/* Fire up auxiliary threads to do pass 1 forward FFT */ \
\
/* pass1: */	pass1_forward_fft_setup(b0b); \
\
/* Copy data from the FFT data area back to the scratch area */ \
\
b0b:	copy_fft_data_to_scratch(64, clm); \
\
/* Do FFT levels 1,2,3 */ \
/* */ \
/* On input the 64-byte cache lines hold these data values: */ \
/*	0K	128K	256K	384K */ \
/*	2K	... */ \
/*	... */ \
/*	126K */ \
/* On output the 64-byte cache lines hold these data values: */ \
/*	0K	64K	128K	192K */ \
/*	2K	... */ \
/*	... */ \
/*	62K	... */ \
/*	256K	... */ \
/*	... */ \
/* */ \
/* Do 64 eight_reals_first_fft macros (FFT levels 1,2,3) */ \
/*	distance between fft data elements is 64K */ \
/*	do 32 iterations each processing 16 data values */ \
\
/* Do FFT levels 4,5 */ \
/* Values 0K-63K is real data, 64K-127K is semi-real data, 128K-511K is */ \
/* complex data. */ \
/* */ \
/* On input the 64-byte cache lines hold these data values: */ \
/*	0K	64K	128K	192K */ \
/*	2K	... */ \
/*	... */ \
/*	62K	... */ \
/*	256K	... */ \
/*	... */ \
/* On output the 64-byte cache lines hold these data values: */ \
/*	0K	16K	32K	48K */ \
/*	2K	... */ \
/*	... */ \
/*	14K	... */ \
/*	64K	... */ \
/*	... */ \
/* */ \
/* Do 8 eight_reals_fft_2 macros */ \
/* Do 8 nop_two_two_complex_fft_2 macros */ \
/* Do 16 four_complex_fft macros */ \
/*	distance between fft data elements is 16K */ \
/*	do 8 macros each processing 32 data values */ \
/* */ \
/* Do 32 four_complex_fft macros */ \
/*	distance between fft data elements is 16K */ \
/*	do 8 macros each processing 32 data values */ \
\
/* Do FFT levels 6 */ \
/* Values 0K-15K is real data, 16K-31K is semi-real data, 32K-511K is */ \
/* complex data. */ \
/* */ \
/* On input the 64-byte cache lines hold these data values: */ \
/*	0K	16K	32K	48K */ \
/*	2K	... */ \
/*	... */ \
/*	14K	... */ \
/*	64K	... */ \
/*	... */ \
/* On output the 64-byte cache lines hold these data values: */ \
/*	0K	8K	16K	24K */ \
/*	2K	... */ \
/*	... */ \
/*	6K	... */ \
/*	32K	... */ \
/*	... */ \
/* */ \
/* Do 4 eight_reals_fft_1 macros */ \
/* Do 4 two_two_complex_fft macros */ \
/*	distance between fft data elements is 8K */ \
/*	do 4 macros each processing 16 data values */ \
/* */ \
/* Do 56 two_two_complex_fft macros */ \
/*	distance between fft data elements is 8K */ \
/*	do 28 macros each processing 16 data values */ \
\
yes_fft:start_timer(21); \
rsi = (uintptr_t)g->scratch_area;	/* Get address of scratch area */ \
rcx = (uintptr_t)g->data_addr; /* Load source address for prefetching */ \
pfing(rbp = (uintptr_t)g->data_prefetch);/* Set prefetch distance */ \
pfing(rbp -= rcx); \
pfing(rcx += +32*blkdst);		/* Create prefetching pointer */ \
rdi = (uintptr_t)g->sincos3;/* Load sin/cos pointer */ \
rbx = (uintptr_t)g->sincos2;/* Load sin/cos pointer */ \
/* 4 iterations */ \
\
for(unsigned int loopA = 4; loopA; loopA--) { \
	xtouchlines(8);			/* Preload TLBs */ \
	\
	for(unsigned int loopB = clm; loopB; loopB--) { \
		for(unsigned int loopC = 4; loopC; loopC--) { \
			for(unsigned int loopD = 2; loopD; loopD--) { \
				for(unsigned int loopE = 2; loopE; loopE--) { \
					x2cl_eight_reals_fft(rsi, 64, 4*clmblkdst8); \
				} \
				rsi += -2*64+4*clmblkdst;	/* Next source pointer */ \
			} \
			rsi += -2*4*clmblkdst+clmblkdst8;/* Next source pointer */ \
			L2prefetch128(u8ptr(rcx+rbp)); \
			pfing(rcx += blkdst);		/* Next prefetch pointer */ \
		} \
		\
		rsi += -4*clmblkdst8;	/* Restore source pointer */ \
		for(unsigned int loopC = 2; loopC; loopC--) { \
			xloop1sp(b2b, 4*clmblkdst, x4cl_eight_reals_fft_2(rsi, 64, clmblkdst8, 2*clmblkdst8, rdi)); \
		} \
		\
		rsi += -2*4*clmblkdst+4*clmblkdst8;/* Next source pointer */ \
		rdi += XMM_SCD;		/* Next sine/cosine pointer */ \
		for(unsigned int loopC = 2; loopC; loopC--) { \
			xloop1sp(b3b, 4*clmblkdst, x4cl_four_complex_fft(rsi, 64, clmblkdst8, 2*clmblkdst8, rdi)); \
		} \
		rdi += -XMM_SCD;		/* Restore sine/cosine pointer */ \
		\
		rsi += -2*4*clmblkdst-4*clmblkdst8;/* Restore src pointer */ \
		for(unsigned int loopC = 2; loopC; loopC--) { \
			x2cl_eight_reals_fft_1(rsi, 64, 4*clmblkdst); \
		} \
		\
		rsi += -2*64+clmblkdst8;	/* Load source pointer */ \
		rbx += XMM_SCD;		/* Next sine/cosine pointer */ \
		/* 7 iterations */ \
		for(unsigned int loopC = 7; loopC; loopC--) { \
			for(unsigned int loopD = 2; loopD; loopD--) { \
				x2cl_two_complex_fft(rsi, 64, 4*clmblkdst, rbx); \
			} \
			rsi += -2*64+clmblkdst8;	/* Next source pointer */ \
			rbx += 2*XMM_SCD;		/* Next sine/cosine pointer */ \
		} \
		rbx += -15*XMM_SCD;	/* Restore sine/cosine pointer */ \
		\
		EXEC248##clm(rsi += -8*clmblkdst8+2*64);	/* Next source pointer */ \
		EXEC248##clm(pfing(rcx += -8*blkdst+128));	/* Next prefetch pointer */ \
	} \
	EXEC248##clm(rsi += -clm*2*64+clmblkdst);/* Next source pointer */ \
	EXEC248##clm(pfing(rcx += -clm*128+8*blkdst));	/* Next prefetch pointer */ \
	EXEC1##clm(rsi += -8*clmblkdst8+clmblkdst);/* Next source pointer */ \
	\
} \
rsi += -4*clmblkdst;	/* Restore source pointer */ \
pfing(rcx += -64*blkdst);		/* Restore prefetch pointer */ \
end_timer(21); \
\
/* Do FFT levels 7,8 */ \
/* Values 0K-7K is real data, 8K-15K is semi-real data, 16K-511K is */ \
/* complex data. */ \
/* */ \
/* On input the 64-byte cache lines hold these data values: */ \
/*	0K	8K	16K	24K */ \
/*	2K	... */ \
/*	... */ \
/*	6K	... */ \
/*	32K	... */ \
/*	... */ \
/* On output the 64-byte cache lines hold these data values: */ \
/*	0K	2K	4K	6K */ \
/*	8K	... */ \
/*	... */ \
\
/* Do 1 eight_reals_fft_2 macros */ \
/* Do 1 nop_two_two_complex_fft_2 macros */ \
/* Do 2 four_complex_fft macros */ \
/*	distance between fft data elements is 2K */ \
/*	do 1 macros each processing 32 data values */ \
\
start_timer(24); \
rdi = (uintptr_t)g->sincos1;/* Load sin/cos pointer */ \
for(unsigned int loopA = 2*clm; loopA; loopA--) { \
	g4cl_eight_reals_fft_2(rsi, 64, clmblkdst, 2*clmblkdst, rcx, 64, blkdst, 2*blkdst); \
} \
\
/* Do 60 four_complex_fft macros */ \
/*	distance between fft data elements is 2K */ \
/*	do 15 macros each processing 32 data values */ \
\
rsi += -2*clm*64+4*clmblkdst;/* Next source pointer */ \
rcx += -2*clm*64+4*blkdst;	/* Next dest pointer */ \
rdi += XMM_SCD;		/* Next sine/cosine pointer */ \
/* eax = 8*256+1; */ /* 15 iterations of 2*clm */ \
for(unsigned int loopA = 15; loopA; ) { \
	do{ \
		for(unsigned int loopC = clm; loopC; loopC--) { \
			for(unsigned int loopD = 2; loopD; loopD--) { \
				g4cl_four_complex_fft(rsi, 64, clmblkdst, 2*clmblkdst, rcx, 64, blkdst, 2*blkdst); \
				L2prefetch128(u8ptr(rcx-64+rbp)); \
				L2prefetch128(u8ptr(rcx-64+blkdst+rbp)); \
				pfing(rbp += -64+2*blkdst); \
			} \
			pfing(rbp += 128-4*blkdst); \
		} \
		rsi += -2*clm*64+4*clmblkdst;/* Next source pointer */ \
		rcx += -2*clm*64+4*blkdst;	/* Next dest pointer */ \
		rdi += 2*XMM_SCD;		/* Next sine/cosine pointer */ \
	}while(--loopA & 1); \
	rsi += -8*clmblkdst+clmblkdst8; \
} \
end_timer(24); \
\
/* Work on next 4 sets of 256 values */ \
\
endlp:	xpass1_get_next_block(pass2, c0b, b0b); /* Get next block and jump */ \
/* done: */	xfft_footer_hg2; \
\
\
\
/* Perform a 512K-element negayclic FFT.  This is done in two passes.  Pass 1 */ \
/* does 8 levels, pass 2 uses common code to perform the last 11 levels. */ \
\

#define hg_pass1sclevels8complex4(procname, xpass2name, pass2_levels, clm) \
 /* (fftname, pass2_levels, clmarg) */ \
\
/* set_FFT_constants(fftname, 0, 1, pass2_levels, clmarg); */ \
xfft_header_hg2(procname, xpass2name, pass2_levels, clm); \
\
/* Do the last 11 FFT levels */ \
/*	do 64 pass2_11_levels_complex */ \
\
/* pass2: */	start_timer(0); \
xpass2name(g); \
end_timer(0); \
if(g->ffttype == 1)		/* We're done if FFTing only */ \
	return; \
\
/* Do the inverse FFT */ \
\
pass1_inverse_fft_setup(c0b, procname); \
\
/* Do inverse FFT levels 7,8 */ \
/* On input the 64-byte cache lines hold these data values: */ \
/*	0K	2K	4K	6K */ \
/*	8K	10K	12K	14K */ \
/*	... */ \
/* On output the 64-byte cache lines hold these data values: */ \
/*	0K	4K	8K	12K */ \
/*	2K	... */ \
/*	16K	... */ \
/*	... */ \
\
/* Do 64 four_complex_unfft macros */ \
/*	distance between fft data elements is 2K */ \
/*	do 32 macros each processing 16 data values */ \
\
c0b:	start_timer(16); \
rsi = (uintptr_t)g->data_addr;/* Get source address */ \
rcx = (uintptr_t)g->scratch_area;	/* Get address of scratch area */ \
rdi = (uintptr_t)g->sincos1;/* Load sin/cos pointer */ \
/* 32 iterations of 2*clm */ \
pfing(rbp = (uintptr_t)g->data_prefetch);/* Set prefetch distance */ \
pfing(rbp -= rsi); \
for(unsigned int loopA = 8; loopA; loopA--) { \
	for(unsigned int loopB = 4; loopB; loopB--) { \
		for(unsigned int loopC = 2*clm; loopC; loopC--) { \
			g2cl_four_complex_unfft(rsi, 64, blkdst, rcx, 64, clmblkdst); \
		} \
		rsi += -2*clm*64+2*blkdst;	/* Next source pointer */ \
		rcx += -2*clm*64+2*clmblkdst;/* Next dest pointer */ \
		rdi += XMM_SCD;		/* Next sine/cosine pointer */ \
	} \
	rcx += -8*clmblkdst+clmblkdst8; \
} \
rcx = rsi-32*2*blkdst;	/* Create prefetch pointer */ \
rsi = (uintptr_t)g->scratch_area;	/* Restore source pointer */ \
end_timer(16); \
\
/* Do inverse FFT levels 5,6 */ \
/* On input the 64-byte cache lines hold these data values: */ \
/*	0K	4K	8K	12K */ \
/*	2K	... */ \
/*	16K	... */ \
/*	... */ \
/* On output the 64-byte cache lines hold these data values: */ \
/*	0K	16K	32K	48K */ \
/*	2K	... */ \
/*	... */ \
/*	14K */ \
/*	64K	... */ \
/*	... */ \
\
/* Do 64 four_complex_unfft macros */ \
/*	distance between fft data elements is 8K */ \
/*	do 16 macros each processing 32 data values */ \
\
start_timer(17); \
rdi = (uintptr_t)g->sincos2;/* Load sin/cos pointer */ \
/* 8 iterations of 2 iters of 2*clm */ \
for(unsigned int loopA = 8; loopA; loopA--) { \
	for(unsigned int loopB = 2; loopB; loopB--) { \
		xtouch(xptr(rcx+(clm-1)*128));	/* Preload TLB */ \
		xloop1_clm##clm(c5b, clmblkdst, x4cl_four_complex_unfft(rsi, 64, 2*clmblkdst, 4*clmblkdst, rdi)); \
	} \
	rsi += -2*clmblkdst+clmblkdst8;/* Next source pointer */ \
	rdi += XMM_SCD;		/* Next sine/cosine pointer */ \
} \
rsi += -8*clmblkdst8;	/* Restore source pointer */ \
end_timer(17); \
\
/* Do inverse FFT levels 4 */ \
/* On input the 64-byte cache lines hold these data values: */ \
/*	0K	16K	32K	48K */ \
/*	2K	... */ \
/*	... */ \
/*	14K */ \
/*	64K	... */ \
/*	... */ \
/* On output the 64-byte cache lines hold these data values: */ \
/*	0K	32K	64K	96K */ \
/*	2K	... */ \
/*	... */ \
/*	30K */ \
/*	128K	... */ \
/*	... */ \
\
/* Do 64 two_two_complex_unfft macros */ \
/*	distance between fft data elements is 32K */ \
/*	do 32 macros each processing 16 data values */ \
\
start_timer(18); \
rdi = (uintptr_t)g->sincos3;/* Load sin/cos pointer */ \
/* 4 iterations of 2 iters of 4*2*clm */ \
for(unsigned int loopA = 4; loopA; loopA--) { \
	for(unsigned int loopB = 2; loopB; loopB--) { \
		xtouch(xptr(rcx+(clm-1)*128));	/* Preload TLB */ \
		xloop4_clm##clm(c3b, clmblkdst, x2cl_two_complex_unfft(rsi, 64, clmblkdst8));		/* 4* 2*clm iterations */ \
	} \
	rsi += -8*clmblkdst+2*clmblkdst8;/* Next source pointer */ \
	rdi += XMM_SCD;		/* Next sine/cosine pointer */ \
} \
rsi += -4*2*clmblkdst8;	/* Restore source pointer */ \
end_timer(18); \
\
/* Do inverse FFT levels 1,2,3 */ \
/* On input the 64-byte cache lines hold these data values: */ \
/*	0K	32K	64K	96K */ \
/*	2K	... */ \
/*	... */ \
/*	30K */ \
/*	128K	... */ \
/*	... */ \
/* On output the 64-byte cache lines hold these data values: */ \
/*	0K	128K	256K	384K */ \
/*	2K	... */ \
/*	... */ \
/*	126K */ \
\
/* Do 64 four_complex_last_unfft macros */ \
/*	distance between fft data elements is 64K */ \
/*	do 16 macros each processing 32 data values */ \
\
start_timer(19); \
rdi = (uintptr_t)g->plus1_premults;/* Address of premultiplier table */ \
/* 2 iters of 4 iterations of 2*2*clm */ \
for(unsigned int loopA = 2; loopA; loopA--) { \
	for(unsigned int loopB = 4; loopB; loopB--) { \
		xtouch(xptr(rcx+(clm-1)*128));	/* Preload TLB */ \
		xloop2edi_clm##clm(c1b, clmblkdst, 128, x4cl_four_complex_last_unfft(rsi, 64, 2*clmblkdst8, 4*clmblkdst8, 16*128));	/* 2 * 2*clm iterations */ \
	} \
	rsi += -4*2*clmblkdst+clmblkdst8;/* Next source pointer */ \
} \
end_timer(19); \
\
/* Normalize these values */ \
\
pass1_normalize(1, yes_fft); \
\
/* Copy data from the scratch area back to the FFT data area */ \
\
copy_scratch_data_to_fft(64, clm); \
goto	endlp;			/* Skip FFT code if flag not set */ \
\
/* Fire up auxiliary threads to do pass 1 forward FFT */ \
\
/* pass1: */	pass1_forward_fft_setup(b0b); \
\
/* Copy data from the FFT data area back to the scratch area */ \
\
b0b:	copy_fft_data_to_scratch(64, clm); \
\
/* Do FFT levels 1,2,3 */ \
/* */ \
/* On input the 64-byte cache lines hold these data values: */ \
/*	0K	128K	256K	384K */ \
/*	2K	... */ \
/*	... */ \
/*	126K */ \
/* On output the 64-byte cache lines hold these data values: */ \
/*	0K	64K	128K	192K */ \
/*	2K	... */ \
/*	... */ \
/*	62K	... */ \
/*	256K	... */ \
/*	... */ \
\
/* Do 64 four_complex_first_fft macros (FFT levels 1,2,3) */ \
/*	distance between fft data elements is 64K */ \
/*	do 32 iterations each processing 16 data values */ \
\
yes_fft:start_timer(21); \
rsi = (uintptr_t)g->scratch_area;	/* Get address of scratch area */ \
rdi = (uintptr_t)g->plus1_premults;/* Address of premultiplier table */ \
rcx = (uintptr_t)g->data_addr; /* Load source address for prefetching */ \
pfing(rbp = (uintptr_t)g->data_prefetch);/* Set prefetch distance */ \
pfing(rbp -= rcx); \
pfing(rcx += 32*blkdst);		/* Create prefetching pointer */ \
/* 4 iters of 2 iters of 4*2*clm */ \
for(unsigned int loopA = 4; loopA; loopA--) { \
	for(unsigned int loopB = 2; loopB; loopB--) { \
		xtouch(xptr(rcx+(clm-1)*128));	/* Preload TLB */ \
		xloop4edi_clm##clm(b1b, clmblkdst, 128, x2cl_four_complex_first_fft_scratch(rsi, 64, 4*clmblkdst8));	/* 4 * 2*clm iterations */ \
	} \
	rsi += -2*4*clmblkdst+clmblkdst8; \
} \
rsi += -4*clmblkdst8;	/* Restore source pointer */ \
end_timer(21); \
\
/* Do FFT levels 4 */ \
/* */ \
/* On input the 64-byte cache lines hold these data values: */ \
/*	0K	32K	64K	96K */ \
/*	2K	... */ \
/*	... */ \
/*	30K */ \
/*	128K	192K	256K	320K */ \
/*	130K */ \
/*	... */ \
/*	190K */ \
/* On output the 64-byte cache lines hold these data values: */ \
/*	0K	32K	64K	96K */ \
/*	2K	... */ \
/*	... */ \
/*	30K	... */ \
/*	128K	... */ \
/*	... */ \
\
/* Do 64 two_two_complex_fft macros */ \
/*	distance between fft data elements is 32K */ \
/*	do 32 macros each processing 16 data values */ \
\
start_timer(22); \
rdi = (uintptr_t)g->sincos3;/* Load sin/cos pointer */ \
/* 2 iters of 2 iterations of 4*2*clm */ \
for(unsigned int loopA = 2; loopA; loopA--) { \
	for(unsigned int loopB = 2; loopB; loopB--) { \
		for(unsigned int loopC = 2; loopC; loopC--) { \
			xtouch(xptr(rcx+(clm-1)*128));	/* Preload TLB */ \
			xloop4_clm##clm(b3b, clmblkdst, x2cl_two_complex_fft(rsi, 64, 2*clmblkdst8, rdi));		/* 4 * 2*clm iterations */ \
		} \
		rsi += -2*4*clmblkdst+clmblkdst8; \
	} \
	rsi += -2*clmblkdst8+4*clmblkdst8;/* Next source pointer */ \
	rdi += 2*XMM_SCD;		/* Next sin/cos pointer */ \
} \
rsi += -2*4*clmblkdst8;	/* Restore source pointer */ \
end_timer(22); \
\
/* Do FFT levels 5,6 */ \
/* */ \
/* On input the 64-byte cache lines hold these data values: */ \
/*	0K	32K	64K	96K */ \
/*	2K	... */ \
/*	... */ \
/*	30K	... */ \
/*	128K	... */ \
/*	... */ \
/* On output the 64-byte cache lines hold these data values: */ \
/*	0K	8K	16K	24K */ \
/*	2K	... */ \
/*	... */ \
/*	6K	... */ \
/*	32K	... */ \
/*	... */ \
\
/* Do 64 four_complex_fft macros */ \
/*	distance between fft data elements is 8K */ \
/*	do 16 macros each processing 32 data values */ \
\
start_timer(23); \
rdi = (uintptr_t)g->sincos2;/* Load sin/cos pointer */ \
/* 4 iterations of 4 iters of 2*2*clm */ \
for(unsigned int loopA = 4; loopA; loopA--) { \
	for(unsigned int loopB = 4; loopB; loopB--) { \
		xtouch(xptr(rcx+(clm-1)*128));	/* Preload TLB */ \
		xloop1_clm##clm(b5b, clmblkdst, x4cl_four_complex_fft(rsi, 64, 4*clmblkdst, clmblkdst8, rdi));		/* 2 * 2*clm iterations */ \
	} \
	rsi += -4*clmblkdst+2*clmblkdst8;/* Next source pointer */ \
	rdi += 2*XMM_SCD;		/* Next sine/cosine pointer */ \
} \
rsi += -4*2*clmblkdst8;	/* Restore source pointer */ \
pfing(rcx += -64*blkdst);		/* Restore prefetch pointer */ \
end_timer(23); \
\
/* Do FFT levels 7,8 */ \
/* */ \
/* On input the 64-byte cache lines hold these data values: */ \
/*	0K	8K	16K	24K */ \
/*	2K	... */ \
/*	... */ \
/*	6K	... */ \
/*	32K	... */ \
/*	... */ \
/* On output the 64-byte cache lines hold these data values: */ \
/*	0K	2K	4K	6K */ \
/*	8K	... */ \
/*	... */ \
\
/* Do 64 four_complex_fft macros */ \
/*	distance between fft data elements is 2K */ \
/*	do 16 macros each processing 32 data values */ \
\
start_timer(24); \
rdi = (uintptr_t)g->sincos1;/* Load sin/cos pointer */ \
/* 8 iters of 2 of 2*clm */ \
for(unsigned int loopA = 8; loopA; loopA--) { \
	for(unsigned int loopB = 2; loopB; loopB--) { \
		for(unsigned int loopC = clm; loopC; loopC--) { \
			for(unsigned int loopD = 2; loopD; loopD--) { \
				g4cl_four_complex_fft(rsi, 64, clmblkdst, 2*clmblkdst, rcx, 64, blkdst, 2*blkdst); \
				L2prefetch128(u8ptr(rcx-64+rbp)); \
				L2prefetch128(u8ptr(rcx-64+blkdst+rbp)); \
				pfing(rbp += -64+2*blkdst); \
			} \
			pfing(rbp += 128-4*blkdst); \
		} \
	rsi += -2*clm*64+4*clmblkdst;/* Next source pointer */ \
		rcx += -2*clm*64+4*blkdst;	/* Next dest pointer */ \
		rdi += 2*XMM_SCD;		/* Next sine/cosine pointer */ \
	} \
	rsi += -8*clmblkdst+clmblkdst8; \
} \
end_timer(24); \
\
/* Work on next 4 sets of 256 values */ \
\
endlp:	xpass1_get_next_block(pass2, c0b, b0b); /* Get next block and jump */ \
/* done: */	xfft_footer_hg2; \
\
\
\
/* Perform a 640K-element FFT.  This is done in two passes.  Pass 1 */ \
/* does 9 levels, pass 2 uses common code to perform the last 11 levels. */ \
\

#define hg_pass1sclevels9pfa5(procname, xpass2name, pass2_levels, clm) \
 /* (fftname, pass2_levels, clmarg) */ \
\
/* set_FFT_constants(fftname, 0, 0, pass2_levels, clmarg); */ \
xfft_header_hg2(procname, xpass2name, pass2_levels, clm); \
\
/* Do the last 11 FFT levels */ \
/*	do 1 pass2_11_levels_real and 79 pass2_11_levels_complex */ \
\
/* pass2: */	start_timer(0); \
xpass2name(g); \
end_timer(0); \
if(g->ffttype == 1)		/* We're done if FFTing only */ \
	return; \
\
/* Do the inverse FFT */ \
\
pass1_inverse_fft_setup(c0b, procname); \
\
/* Do inverse FFT levels 8,9 */ \
/* On input the 64-byte cache lines hold these data values: */ \
/*	0K	2K	4K	6K */ \
/*	8K	10K	12K	14K */ \
/*	... */ \
/* On output the 64-byte cache lines hold these data values: */ \
/*	0K	4K	8K	12K */ \
/*	2K	... */ \
/*	16K	... */ \
/*	... */ \
\
/* Do 1 eight_reals_unfft_2 macros */ \
/* Do 1 nop_two_two_complex_unfft_2 macros */ \
/*	distance between fft data elements is 2K */ \
/*	do 1 macros each processing 16 data values */ \
\
c0b:	start_timer(16); \
rsi = (uintptr_t)g->data_addr;/* Get source address */ \
rcx = (uintptr_t)g->scratch_area;	/* Get address of scratch area */ \
pfing(rbp = (uintptr_t)g->data_prefetch);/* Set prefetch distance */ \
pfing(rbp -= rsi); \
for(unsigned int loopA = 2*clm; loopA; loopA--) { \
	g2cl_eight_reals_unfft_2(rsi, 64, blkdst, rcx, 64, clmblkdst); \
} \
\
/* Do 78 four_complex_unfft macros */ \
/*	distance between fft data elements is 2K */ \
/*	do 39 macros each processing 16 data values */ \
\
rsi += -2*clm*64+2*blkdst;	/* Next source pointer */ \
rcx += -2*clm*64+2*clmblkdst;/* Next dest pointer */ \
rdi = (uintptr_t)g->sincos1;/* Load sin/cos pointer */ \
/* al = 256/4+10; */ /* 39 iterations of 2*clm */ \
for(unsigned int loopA = 39; loopA; ) { \
	do{ \
		for(unsigned int loopC = 2*clm; loopC; loopC--) { \
			g2cl_four_complex_unfft(rsi, 64, blkdst, rcx, 64, clmblkdst); \
		} \
		rsi += -2*clm*64+2*blkdst;	/* Next source pointer */ \
		rcx += -2*clm*64+2*clmblkdst;/* Next dest pointer */ \
		rdi += XMM_SCD;		/* Next sine/cosine pointer */ \
	}while(--loopA & 3); \
	rcx += -8*clmblkdst+clmblkdst8;/* Next dest pointer */ \
} \
rcx = rsi-40*2*blkdst;	/* Create prefetch pointer */ \
rsi = (uintptr_t)g->scratch_area;	/* Restore source pointer */ \
end_timer(16); \
\
/* Do inverse FFT levels 6,7 */ \
/* On input the 64-byte cache lines hold these data values: */ \
/*	0K	4K	8K	12K */ \
/*	2K	... */ \
/*	16K	... */ \
/*	... */ \
/* On output the 64-byte cache lines hold these data values: */ \
/*	0K	16K	32K	48K */ \
/*	2K	... */ \
/*	... */ \
/*	14K */ \
/*	64K	... */ \
/*	... */ \
\
/* Do 4 eight_reals_unfft_2 macros */ \
/* Do 4 nop_two_two_complex_unfft_2 macros */ \
/*	distance between fft data elements is 8K */ \
/*	do 2 macros each processing 32 data values */ \
\
start_timer(17); \
/* 2 iterations of 2*clm */ \
for(unsigned int loopA = 2; loopA; loopA--) { \
	for(unsigned int loopB = 2*clm; loopB; loopB--) { \
		x4cl_eight_reals_unfft_2(rsi, 64, 2*clmblkdst, 4*clmblkdst); \
	} \
	rsi += -2*clm*64+clmblkdst;/* Next source pointer */ \
} \
\
/* Do 72 four_complex_unfft macros */ \
/*	distance between fft data elements is 8K */ \
/*	do 18 macros each processing 32 data values */ \
\
rsi += -2*clmblkdst+clmblkdst8;/* Next source pointer */ \
rdi = (uintptr_t)g->sincos2;/* Load sin/cos pointer */ \
/* 9 iterations of 2 iters of 2*clm */ \
for(unsigned int loopA = 9; loopA; loopA--) { \
	for(unsigned int loopB = 2; loopB; loopB--) { \
		xtouch(xptr(rcx+(clm-1)*128));	/* Preload TLB */ \
		xloop1_clm##clm(c5b, clmblkdst, x4cl_four_complex_unfft(rsi, 64, 2*clmblkdst, 4*clmblkdst, rdi)); \
	} \
	rsi += -2*clmblkdst+clmblkdst8;/* Next source pointer */ \
	rdi += XMM_SCD;		/* Next sine/cosine pointer */ \
} \
rsi += -10*clmblkdst8;	/* Restore source pointer */ \
end_timer(17); \
\
/* Do inverse FFT levels 4,5 */ \
/* On input the 64-byte cache lines hold these data values: */ \
/*	0K	16K	32K	48K */ \
/*	2K	... */ \
/*	... */ \
/*	14K */ \
/*	64K	... */ \
/*	... */ \
/* On output the 64-byte cache lines hold these data values: */ \
/*	0K	32K	64K	96K */ \
/*	2K	... */ \
/*	... */ \
/*	30K */ \
/*	128K	192K	256K	320K */ \
/*	... */ \
/* */ \
/* Do 16 eight_reals_unfft_2 macros */ \
/*	distance between fft data elements is 32K */ \
/*	do 8 macros each processing 16 data values */ \
/* */ \
/* Do 64 four_complex_unfft macros */ \
/*	distance between fft data elements is 32K */ \
/*	do 16 macros each processing 32 data values */ \
\
/* Do inverse FFT levels 1,2,3 */ \
/* On input the 64-byte cache lines hold these data values: */ \
/*	0K	32K	64K	96K */ \
/*	2K	... */ \
/*	... */ \
/*	30K */ \
/*	128K	192K	256K	320K */ \
/*	... */ \
/* On output the 64-byte cache lines hold these data values: */ \
/*	0K	32K	320K	352K */ \
/*	2K	... */ \
/*	... */ \
/*	30K */ \
/*	64K	192K	384K	512K */ \
/*	66K */ \
/*	... */ \
/*	190K */ \
\
/* Do 128 five_reals_last_unfft macros */ \
/*	distance between fft data elements is 128K */ \
/*	do 16 macros each processing 40 data values */ \
\
start_timer(18); \
rdi = (uintptr_t)g->sincos3;/* Load sin/cos pointer */ \
/* 8 iterations */ \
\
for(unsigned int loopA = 8; loopA; loopA--) { \
	xtouchlines(3);			/* Preload TLBs */ \
	\
	for(unsigned int loopB = clm; loopB; loopB--) { \
		for(unsigned int loopC = 2; loopC; loopC--) { \
			x2cl_half_eight_reals_unfft_2(rsi, 64, clmblkdst8); \
		} \
		\
		rsi += -2*64+2*clmblkdst8;	/* Next source pointer */ \
		for(unsigned int loopC = 2; loopC; loopC--) { \
			for(unsigned int loopD = 2; loopD; loopD--) { \
				x4cl_four_complex_unfft(rsi, 64, clmblkdst8, 2*clmblkdst8, rdi); \
			} \
			rsi += -2*64+4*clmblkdst8;	/* Restore source pointer */ \
			rdi += XMM_SCD;		/* Next sine/cosine pointer */ \
		} \
		L2prefetch128(u8ptr(rcx+rbp)); \
		pfing(rcx += blkdst);		/* Next prefetch pointer */ \
		rdi += -2*XMM_SCD;		/* Restore sine/cosine pointer */ \
		\
		rsi += -2*4*clmblkdst8-2*clmblkdst8;/* Restore source pointer */ \
		xfive_reals_unfft_preload;	/* Preload consts for five real macros */ \
		for(unsigned int loopC = 2; loopC; loopC--) { \
			xloop1sp(c1b, clmblkdst8, x5cl_five_reals_last_unfft(rsi, 64, 2*clmblkdst8)); \
		} \
		\
		EXEC248##clm(rsi += -2*clmblkdst8+2*64);	/* Next source pointer */ \
		EXEC248##clm(pfing(rcx += -3*blkdst+128));	/* Next prefetch pointer */ \
	} \
	EXEC248##clm(rsi += -clm*2*64+clmblkdst);/* Next source pointer */ \
	EXEC248##clm(pfing(rcx += -clm*128+3*blkdst));	/* Next prefetch pointer */ \
	EXEC1##clm(rsi += -2*clmblkdst8+clmblkdst);/* Next source pointer */ \
	\
} \
end_timer(18); \
\
/* Normalize these values */ \
\
pass1_normalize(1, yes_fft); \
\
/* Copy data from the scratch area back to the FFT data area */ \
\
copy_scratch_data_to_fft(80, clm); \
goto	endlp;			/* Skip FFT code if flag not set */ \
\
/* Fire up auxiliary threads to do pass 1 forward FFT */ \
\
/* pass1: */	pass1_forward_fft_setup(b0b); \
\
/* Copy data from the FFT data area back to the scratch area */ \
\
b0b:	copy_fft_data_to_scratch(80, clm); \
\
/* Do FFT levels 1,2,3 */ \
/* */ \
/* On input the 64-byte cache lines hold these data values: */ \
/*	0K	32K	320K	352K */ \
/*	2K	... */ \
/*	... */ \
/*	30K */ \
/*	64K	192K	384K	512K */ \
/*	66K */ \
/*	... */ \
/*	190K */ \
/* On output the 64-byte cache lines hold these data values: */ \
/*	0K	32K	64K	96K */ \
/*	2K	... */ \
/*	... */ \
/*	30K */ \
/*	128K	256K	384K	512K */ \
/*	130K */ \
/*	... */ \
/*	254K */ \
/* */ \
/* Do 128 five_reals_first_fft macros (FFT levels 1,2,3) */ \
/*	distance between fft data elements is 128K */ \
/*	do 16 iterations each processing 40 data values */ \
\
/* Do FFT levels 4,5 */ \
/* Values 0K-127K is real data, 128K-639K is complex data. */ \
/* */ \
/* On input the 64-byte cache lines hold these data values: */ \
/*	0K	32K	64K	96K */ \
/*	2K	... */ \
/*	... */ \
/*	30K */ \
/*	128K	256K	384K	512K */ \
/*	130K */ \
/*	... */ \
/*	254K */ \
/* On output the 64-byte cache lines hold these data values: */ \
/*	0K	32K	64K	96K */ \
/*	2K	... */ \
/*	... */ \
/*	30K	... */ \
/*	128K	... */ \
/*	... */ \
/* */ \
/* Do 16 eight_reals_fft_2 macros */ \
/*	distance between fft data elements is 32K */ \
/*	do 16 macros each processing 8 data values */ \
/* */ \
/* Do 64 four_complex_fft macros */ \
/*	distance between fft data elements is 32K */ \
/*	do 16 macros each processing 32 data values */ \
\
/* Do FFT levels 6,7 */ \
/* Values 0K-31K is real data, 32K-63K is semi-real data, 64K-639K is */ \
/* complex data. */ \
/* */ \
/* On input the 64-byte cache lines hold these data values: */ \
/*	0K	32K	64K	96K */ \
/*	2K	... */ \
/*	... */ \
/*	30K	... */ \
/*	128K	... */ \
/*	... */ \
/* On output the 64-byte cache lines hold these data values: */ \
/*	0K	8K	16K	24K */ \
/*	2K	... */ \
/*	... */ \
/*	6K	... */ \
/*	32K	... */ \
/*	... */ \
/* */ \
/* Do 8 eight_reals_fft_2 macros */ \
/* Do 8 two_two_complex_fft macros */ \
/*	distance between fft data elements is 8K */ \
/*	do 4 macros each processing 32 data values */ \
/* */ \
/* Do 64 four_complex_fft macros */ \
/*	distance between fft data elements is 8K */ \
/*	do 16 macros each processing 32 data values */ \
\
yes_fft:start_timer(21); \
rsi = (uintptr_t)g->scratch_area;	/* Get address of scratch area */ \
rcx = (uintptr_t)g->data_addr; /* Load source address for prefetching */ \
pfing(rbp = (uintptr_t)g->data_prefetch);/* Set prefetch distance */ \
pfing(rbp -= rcx); \
pfing(rcx += 40*blkdst);		/* Create prefetching pointer */ \
rdi = (uintptr_t)g->sincos3;/* Load sin/cos pointer */ \
rbx = (uintptr_t)g->sincos2;/* Load sin/cos pointer */ \
/* 4 iterations */ \
\
for(unsigned int loopA = 4; loopA; loopA--) { \
	for(unsigned int loopB = clm; loopB; loopB--) { \
		xtouchlines(10);			/* Preload TLBs */ \
		\
		xfive_reals_fft_preload;		/* Preload consts for five real macros */ \
		for(unsigned int loopC = 2; loopC; loopC--) { \
			for(unsigned int loopD = 2; loopD; loopD--) { \
				xloop1sp(b1b, 4*clmblkdst, x5cl_five_reals_fft(rsi, 64, 2*clmblkdst8)); \
			} \
			rsi += -2*4*clmblkdst+clmblkdst8; \
		} \
		\
		rsi += -2*clmblkdst8;	/* Restore source pointer */ \
		for(unsigned int loopC = 2; loopC; loopC--) { \
			for(unsigned int loopD = 2; loopD; loopD--) { \
				for(unsigned int loopE = 2; loopE; loopE--) { \
					x1cl_half_eight_reals_fft_2(rsi, 64); \
				} \
				rsi += -2*64+4*clmblkdst; \
			} \
			rsi += -2*4*clmblkdst+clmblkdst8; \
		} \
		\
		/*	bump	rsi, -2*clmblkdst8+2*clmblkdst8;; Next source pointer */ \
		for(unsigned int loopC = 2; loopC; loopC--) { \
			for(unsigned int loopD = 2; loopD; loopD--) { \
				for(unsigned int loopE = 2; loopE; loopE--) { \
					x4cl_four_complex_fft(rsi, 64, 2*clmblkdst8, 4*clmblkdst8, rdi); \
				} \
				rsi += -2*64+4*clmblkdst; \
			} \
			L2prefetch128(u8ptr(rcx+rbp)); \
			pfing(rcx += blkdst);		/* Next prefetch pointer */ \
			rsi += -2*4*clmblkdst+clmblkdst8; \
		} \
		\
		rsi += -2*clmblkdst8-2*clmblkdst8;/* Restore source pointer */ \
		for(unsigned int loopC = 2; loopC; loopC--) { \
			x4cl_eight_reals_fft_2(rsi, 64, 4*clmblkdst, clmblkdst8, rbx); \
		} \
		\
		rsi += -2*64+2*clmblkdst8;	/* Next source pointer */ \
		rbx += XMM_SCD;		/* Next sine/cosine pointer */ \
		for(unsigned int loopC = 4; loopC; loopC--) { \
			xloop1sp(b5b, 2*clmblkdst8, x4cl_four_complex_fft(rsi, 64, 4*clmblkdst, clmblkdst8, rbx)); \
			rbx += 2*XMM_SCD;		/* Next sine/cosine pointer */ \
		} \
		rbx += -9*XMM_SCD;		/* Restore sine/cosine pointer */ \
		\
		EXEC248##clm(rsi += -5*2*clmblkdst8+2*64);/* Next source pointer */ \
		EXEC248##clm(pfing(rcx += -10*blkdst+128));	/* Next prefetch pointer */ \
	} \
	EXEC248##clm(rsi += -clm*2*64+clmblkdst);/* Next source pointer */ \
	EXEC248##clm(pfing(rcx += -clm*128+10*blkdst));	/* Next prefetch pointer */ \
	EXEC1##clm(rsi += -5*2*clmblkdst8+clmblkdst);/* Next source pointer */ \
	\
} \
rsi += -4*clmblkdst;	/* Restore source pointer */ \
pfing(rcx += -80*blkdst);		/* Restore prefetch pointer */ \
end_timer(21); \
\
/* Do FFT levels 8,9 */ \
/* Values 0K-7K is real data, 8K-15K is semi-real data, 16K-619K is */ \
/* complex data. */ \
/* */ \
/* On input the 64-byte cache lines hold these data values: */ \
/*	0K	8K	16K	24K */ \
/*	2K	... */ \
/*	... */ \
/*	6K	... */ \
/*	32K	... */ \
/*	... */ \
/* On output the 64-byte cache lines hold these data values: */ \
/*	0K	2K	4K	6K */ \
/*	8K	... */ \
/*	... */ \
\
/* Do 1 eight_reals_fft_2 macros */ \
/* Do 1 nop_two_two_complex_fft_2 macros */ \
/* Do 2 four_complex_fft macros */ \
/*	distance between fft data elements is 2K */ \
/*	do 1 macros each processing 32 data values */ \
\
start_timer(24); \
rdi = (uintptr_t)g->sincos1;/* Load sin/cos pointer */ \
for(unsigned int loopA = 2*clm; loopA; loopA--) { \
	g4cl_eight_reals_fft_2(rsi, 64, clmblkdst, 2*clmblkdst, rcx, 64, blkdst, 2*blkdst); \
} \
\
/* Do 76 four_complex_fft macros */ \
/*	distance between fft data elements is 2K */ \
/*	do 19 macros each processing 32 data values */ \
\
rsi += -2*clm*64+4*clmblkdst;/* Next source pointer */ \
rcx += -2*clm*64+4*blkdst;	/* Next dest pointer */ \
rdi += XMM_SCD;		/* Next sine/cosine pointer */ \
/* eax = 10*256+1; */ /* 19 iterations of 2*clm */ \
for(unsigned int loopA = 19; loopA; ) { \
	do{ \
		for(unsigned int loopC = clm; loopC; loopC--) { \
			for(unsigned int loopD = 2; loopD; loopD--) { \
				g4cl_four_complex_fft(rsi, 64, clmblkdst, 2*clmblkdst, rcx, 64, blkdst, 2*blkdst); \
				L2prefetch128(u8ptr(rcx-64+rbp)); \
				L2prefetch128(u8ptr(rcx-64+blkdst+rbp)); \
				pfing(rbp += -64+2*blkdst); \
			} \
			pfing(rbp += 128-4*blkdst);	 \
		} \
		rsi += -2*clm*64+4*clmblkdst;/* Next source pointer */ \
		rcx += -2*clm*64+4*blkdst;	/* Next dest pointer */ \
		rdi += 2*XMM_SCD;		/* Next sine/cosine pointer */ \
	}while(--loopA & 1); \
	rsi += -8*clmblkdst+clmblkdst8; \
} \
end_timer(24); \
\
/* Work on next 4 sets of 320 values */ \
\
endlp:	xpass1_get_next_block(pass2, c0b, b0b); /* Get next block and jump */ \
/* done: */	xfft_footer_hg2; \
\
\
\
/* Perform a 768K-element FFT.  This is done in two passes.  Pass 1 */ \
/* does 9 levels, pass 2 uses common code to perform the last 11 levels. */ \
\

#define hg_pass1sclevels9pfa6(procname, xpass2name, pass2_levels, clm) \
 /* (fftname, pass2_levels, clmarg) */ \
\
/* set_FFT_constants(fftname, 0, 0, pass2_levels, clmarg); */ \
xfft_header_hg2(procname, xpass2name, pass2_levels, clm); \
\
/* Do the last 11 FFT levels */ \
/*	do 1 pass2_11_levels_real and 95 pass2_11_levels_complex */ \
\
/* pass2: */	start_timer(0); \
xpass2name(g); \
end_timer(0); \
if(g->ffttype == 1)		/* We're done if FFTing only */ \
	return; \
\
/* Do the inverse FFT */ \
\
pass1_inverse_fft_setup(c0b, procname); \
\
/* Do inverse FFT levels 8,9 */ \
/* On input the 64-byte cache lines hold these data values: */ \
/*	0K	2K	4K	6K */ \
/*	8K	10K	12K	14K */ \
/*	... */ \
/* On output the 64-byte cache lines hold these data values: */ \
/*	0K	4K	8K	12K */ \
/*	2K	... */ \
/*	16K	... */ \
/*	... */ \
\
/* Do 1 eight_reals_unfft_2 macros */ \
/* Do 1 nop_two_two_complex_unfft_2 macros */ \
/*	distance between fft data elements is 2K */ \
/*	do 1 macros each processing 16 data values */ \
\
c0b:	start_timer(16); \
rsi = (uintptr_t)g->data_addr;/* Get source address */ \
rcx = (uintptr_t)g->scratch_area;	/* Get address of scratch area */ \
pfing(rbp = (uintptr_t)g->data_prefetch);/* Set prefetch distance */ \
pfing(rbp -= rsi); \
for(unsigned int loopA = 2*clm; loopA; loopA--) { \
	g2cl_eight_reals_unfft_2(rsi, 64, blkdst, rcx, 64, clmblkdst); \
} \
\
/* Do 94 four_complex_unfft macros */ \
/*	distance between fft data elements is 2K */ \
/*	do 47 macros each processing 16 data values */ \
\
rsi += -2*clm*64+2*blkdst;	/* Next source pointer */ \
rcx += -2*clm*64+2*clmblkdst;/* Next dest pointer */ \
rdi = (uintptr_t)g->sincos1;/* Load sin/cos pointer */ \
/* al = 256/4+12; */ /* 47 iterations of 2*clm */ \
for(unsigned int loopA = 47; loopA; ) { \
	do{ \
		for(unsigned int loopC = 2*clm; loopC; loopC--) { \
			g2cl_four_complex_unfft(rsi, 64, blkdst, rcx, 64, clmblkdst); \
		} \
		rsi += -2*clm*64+2*blkdst;	/* Next source pointer */ \
		rcx += -2*clm*64+2*clmblkdst;/* Next dest pointer */ \
		rdi += XMM_SCD;		/* Next sine/cosine pointer */ \
	}while(--loopA & 3); \
	rcx += -8*clmblkdst+clmblkdst8; \
} \
rcx = rsi-48*2*blkdst;	/* Create prefetch pointer */ \
rsi = (uintptr_t)g->scratch_area;	/* Restore source pointer */ \
end_timer(16); \
\
/* Do inverse FFT levels 6,7 */ \
/* On input the 64-byte cache lines hold these data values: */ \
/*	0K	4K	8K	12K */ \
/*	2K	... */ \
/*	16K	... */ \
/*	... */ \
/* On output the 64-byte cache lines hold these data values: */ \
/*	0K	16K	32K	48K */ \
/*	2K	... */ \
/*	... */ \
/*	14K */ \
/*	64K	... */ \
/*	... */ \
\
/* Do 4 eight_reals_unfft_2 macros */ \
/* Do 4 nop_two_two_complex_unfft_2 macros */ \
/*	distance between fft data elements is 8K */ \
/*	do 2 macros each processing 32 data values */ \
\
start_timer(17); \
/* 2 iterations of 2*clm */ \
for(unsigned int loopA = 2; loopA; loopA--) { \
	for(unsigned int loopB = 2*clm; loopB; loopB--) { \
		x4cl_eight_reals_unfft_2(rsi, 64, 2*clmblkdst, 4*clmblkdst); \
	} \
	rsi += -2*clm*64+clmblkdst;/* Next source pointer */ \
} \
\
/* Do 88 four_complex_unfft macros */ \
/*	distance between fft data elements is 8K */ \
/*	do 22 macros each processing 32 data values */ \
\
rsi += -2*clmblkdst+clmblkdst8;/* Next source pointer */ \
rdi = (uintptr_t)g->sincos2;/* Load sin/cos pointer */ \
/* 11 iterations of 2 iters of 2*clm */ \
for(unsigned int loopA = 11; loopA; loopA--) { \
	for(unsigned int loopB = 2; loopB; loopB--) { \
		xtouch(xptr(rcx+(clm-1)*128));	/* Preload TLB */ \
		xloop1_clm##clm(c5b, clmblkdst, x4cl_four_complex_unfft(rsi, 64, 2*clmblkdst, 4*clmblkdst, rdi)); \
	} \
	rsi += -2*clmblkdst+clmblkdst8;/* Next source pointer */ \
	rdi += XMM_SCD;		/* Next sine/cosine pointer */ \
} \
rsi += -12*clmblkdst8;	/* Restore source pointer */ \
end_timer(17); \
\
/* Do inverse FFT levels 4,5 */ \
/* On input the 64-byte cache lines hold these data values: */ \
/*	0K	16K	32K	48K */ \
/*	2K	... */ \
/*	... */ \
/*	14K */ \
/*	64K	... */ \
/*	... */ \
/* On output the 64-byte cache lines hold these data values: */ \
/*	0K	64K	128K	192K */ \
/*	2K	... */ \
/*	... */ \
/*	30K */ \
/*	256K	... */ \
/*	... */ \
/* */ \
/* Do 16 eight_reals_unfft_2 macros */ \
/* Do 16 nop_two_two_complex_unfft_2 macros */ \
/*	distance between fft data elements is 32K */ \
/*	do 8 macros each processing 32 data values */ \
/* */ \
/* Do 64 four_complex_unfft macros */ \
/*	distance between fft data elements is 32K */ \
/*	do 16 macros each processing 32 data values */ \
\
/* Do inverse FFT levels 1,2,3 */ \
/* On input the 64-byte cache lines hold these data values: */ \
/*	0K	64K	128K	192K */ \
/*	2K	... */ \
/*	... */ \
/*	30K */ \
/*	256K	... */ \
/*	... */ \
/* On output the 64-byte cache lines hold these data values: */ \
/*	0K	192K	384K	576K */ \
/*	... */ \
/*	190K */ \
/* */ \
/* Do 128 six_reals_last_unfft macros */ \
/*	distance between fft data elements is 128K */ \
/*	do 32 macros each processing 24 data values */ \
\
start_timer(18); \
rdi = (uintptr_t)g->sincos3;/* Load sin/cos pointer */ \
/* 8 iterations */ \
xsix_reals_unfft_preload; \
\
for(unsigned int loopA = 8; loopA; loopA--) { \
	xtouchlines(3);			/* Preload TLBs */ \
	\
	for(unsigned int loopB = clm; loopB; loopB--) { \
		for(unsigned int loopC = 2; loopC; loopC--) { \
			x4cl_eight_reals_unfft_2(rsi, 64, clmblkdst8, 2*clmblkdst8); \
		} \
		\
		rsi += -2*64+4*clmblkdst8;	/* Next source pointer */ \
		for(unsigned int loopC = 2; loopC; loopC--) { \
			xloop1sp(c3b, 4*clmblkdst8, x4cl_four_complex_unfft(rsi, 64, clmblkdst8, 2*clmblkdst8, rdi)); \
			rdi += XMM_SCD;		/* Next sine/cosine pointer */ \
		} \
		rdi += -2*XMM_SCD;		/* Restore sine/cosine pointer */ \
		\
		rsi += -3*4*clmblkdst8;	/* Restore source pointer */ \
		for(unsigned int loopC = 4; loopC; loopC--) { \
			for(unsigned int loopD = 2; loopD; loopD--) { \
				x3cl_six_reals_last_unfft(rsi, 64, 4*clmblkdst8); \
			} \
			rsi += -2*64+clmblkdst8; \
		} \
		L2prefetch128(u8ptr(rcx+rbp)); \
		pfing(rcx += blkdst);		/* Next prefetch pointer */ \
		\
		EXEC248##clm(rsi += -4*clmblkdst8+2*64);	/* Next source pointer */ \
		EXEC248##clm(pfing(rcx += -3*blkdst+128));	/* Next prefetch pointer */ \
	} \
	EXEC248##clm(rsi += -clm*2*64+clmblkdst);/* Next source pointer */ \
	EXEC248##clm(pfing(rcx += -clm*128+3*blkdst));	/* Next prefetch pointer */ \
	EXEC1##clm(rsi += -4*clmblkdst8+clmblkdst);/* Next source pointer */ \
	\
} \
end_timer(18); \
\
/* Normalize these values */ \
\
pass1_normalize(1, yes_fft); \
\
/* Copy data from the scratch area back to the FFT data area */ \
\
copy_scratch_data_to_fft(96, clm); \
goto	endlp;			/* Skip FFT code if flag not set */ \
\
/* Fire up auxiliary threads to do pass 1 forward FFT */ \
\
/* pass1: */	pass1_forward_fft_setup(b0b); \
\
/* Copy data from the FFT data area back to the scratch area */ \
\
b0b:	copy_fft_data_to_scratch(96, clm); \
\
/* Do FFT levels 1,2,3 */ \
/* */ \
/* On input the 64-byte cache lines hold these data values: */ \
/*	0K	192K	384K	576K */ \
/*	... */ \
/*	190K */ \
/* On output the 64-byte cache lines hold these data values: */ \
/*	0K	64K	128K	192K */ \
/*	2K	... */ \
/*	... */ \
/*	62K */ \
/*	256K	384K	512K	640K */ \
/*	258K */ \
/*	... */ \
/*	382K */ \
/* */ \
/* Do 128 six_reals_first_fft macros (FFT levels 1,2,3) */ \
/*	distance between fft data elements is 128K */ \
/*	do 32 iterations each processing 24 data values */ \
\
/* Do FFT levels 4,5 */ \
/* Values 0K-127K is real data, 128K-255K is semi-real data, and */ \
/* 256K-767K is complex data. */ \
/* */ \
/* On input the 64-byte cache lines hold these data values: */ \
/*	0K	64K	128K	192K */ \
/*	2K	... */ \
/*	... */ \
/*	62K */ \
/*	256K	384K	512K	640K */ \
/*	258K */ \
/*	... */ \
/*	382K */ \
/* On output the 64-byte cache lines hold these data values: */ \
/*	0K	16K	32K	48K */ \
/*	2K	... */ \
/*	... */ \
/*	14K	... */ \
/*	64K	... */ \
/*	... */ \
/* */ \
/* Do 16 eight_reals_fft_2 macros */ \
/* Do 16 nop_two_two_complex_fft_2 macros */ \
/*	distance between fft data elements is 32K */ \
/*	do 16 macros each processing 16 data values */ \
/* */ \
/* Do 64 four_complex_fft macros */ \
/*	distance between fft data elements is 32K */ \
/*	do 16 macros each processing 32 data values */ \
\
/* Do FFT levels 6,7 */ \
/* Values 0K-31K is real data, 32K-63K is semi-real data, 64K-767K is */ \
/* complex data. */ \
/* */ \
/* On input the 64-byte cache lines hold these data values: */ \
/*	0K	32K	64K	96K */ \
/*	2K	... */ \
/*	... */ \
/*	30K	... */ \
/*	128K	... */ \
/*	... */ \
/* On output the 64-byte cache lines hold these data values: */ \
/*	0K	8K	16K	24K */ \
/*	2K	... */ \
/*	... */ \
/*	6K	... */ \
/*	32K	... */ \
/*	... */ \
/* */ \
/* Do 8 eight_reals_fft_2 macros */ \
/* Do 8 two_two_complex_fft macros */ \
/*	distance between fft data elements is 8K */ \
/*	do 4 macros each processing 32 data values */ \
/* */ \
/* Do 80 four_complex_fft macros */ \
/*	distance between fft data elements is 8K */ \
/*	do 20 macros each processing 32 data values */ \
\
yes_fft:start_timer(21); \
rsi = (uintptr_t)g->scratch_area;	/* Get address of scratch area */ \
rcx = (uintptr_t)g->data_addr; /* Load source address for prefetching */ \
pfing(rbp = (uintptr_t)g->data_prefetch);/* Set prefetch distance */ \
pfing(rbp -= rcx); \
pfing(rcx += 44*blkdst);		/* Create prefetching pointer */ \
rdi = (uintptr_t)g->sincos3;/* Load sin/cos pointer */ \
rbx = (uintptr_t)g->sincos2;/* Load sin/cos pointer */ \
/* 4 iterations */ \
xsix_reals_fft_preload; \
\
for(unsigned int loopA = 4; loopA; loopA--) { \
	xtouchlines(13);			/* Preload TLBs */ \
	\
	for(unsigned int loopB = clm; loopB; loopB--) { \
		for(unsigned int loopC = 4; loopC; loopC--) { \
			for(unsigned int loopD = 2; loopD; loopD--) { \
				for(unsigned int loopE = 2; loopE; loopE--) { \
					x3cl_six_reals_fft(rsi, 64, 4*clmblkdst8); \
				} \
				rsi += -2*64+4*clmblkdst; \
			} \
			L2prefetch128(u8ptr(rcx+rbp)); \
			pfing(rcx += blkdst);		/* Next prefetch pointer */ \
			rsi += -2*4*clmblkdst+clmblkdst8; \
		} \
		\
		rsi += -4*clmblkdst8;	/* Restore source pointer */ \
		for(unsigned int loopC = 2; loopC; loopC--) { \
			for(unsigned int loopD = 2; loopD; loopD--) { \
				for(unsigned int loopE = 2; loopE; loopE--) { \
					x2cl_eight_reals_fft_2(rsi, 64, 2*clmblkdst8); \
				} \
				rsi += -2*64+4*clmblkdst; \
			} \
			rsi += -2*4*clmblkdst+clmblkdst8; \
		} \
		\
		rsi += -2*clmblkdst8+4*clmblkdst8;/* Next source pointer */ \
		for(unsigned int loopC = 2; loopC; loopC--) { \
			for(unsigned int loopD = 2; loopD; loopD--) { \
				xloop1sp(b3b, 4*clmblkdst, x4cl_four_complex_fft(rsi, 64, 2*clmblkdst8, 4*clmblkdst8, rdi)); \
			} \
			rsi += -2*4*clmblkdst+clmblkdst8; \
		} \
		\
		rsi += -2*clmblkdst8-4*clmblkdst8;/* Restore source pointer */ \
		for(unsigned int loopC = 2; loopC; loopC--) { \
			x4cl_eight_reals_fft_2(rsi, 64, 4*clmblkdst, clmblkdst8, rbx); \
		} \
		\
		rsi += -2*64+2*clmblkdst8;	/* Load source pointer */ \
		rbx += XMM_SCD;		/* Next sine/cosine pointer */ \
		/* 5 iterations */ \
		for(unsigned int loopC = 5; loopC; loopC--) { \
			xloop1sp(b5b, 2*clmblkdst8, x4cl_four_complex_fft(rsi, 64, 4*clmblkdst, clmblkdst8, rbx)); \
			rbx += 2*XMM_SCD;		/* Next sine/cosine pointer */ \
		} \
		rbx += -11*XMM_SCD;	/* Restore sine/cosine pointer */ \
		\
		EXEC248##clm(rsi += -6*2*clmblkdst8+2*64);/* Next source pointer */ \
		EXEC248##clm(pfing(rcx += -13*blkdst+128));	/* Next prefetch pointer */ \
	} \
	EXEC248##clm(rsi += -clm*2*64+clmblkdst);/* Next source pointer */ \
	EXEC248##clm(pfing(rcx += -clm*128+13*blkdst));	/* Next prefetch pointer */ \
	EXEC1##clm(rsi += -6*2*clmblkdst8+clmblkdst);/* Next source pointer */ \
	\
} \
rsi += -4*clmblkdst;	/* Restore source pointer */ \
pfing(rcx += -96*blkdst);		/* Restore prefetch pointer */ \
end_timer(21); \
\
/* Do FFT levels 8,9 */ \
/* Values 0K-7K is real data, 8K-15K is semi-real data, 16K-767K is */ \
/* complex data. */ \
/* */ \
/* On input the 64-byte cache lines hold these data values: */ \
/*	0K	8K	16K	24K */ \
/*	2K	... */ \
/*	... */ \
/*	6K	... */ \
/*	32K	... */ \
/*	... */ \
/* On output the 64-byte cache lines hold these data values: */ \
/*	0K	2K	4K	6K */ \
/*	8K	... */ \
/*	... */ \
\
/* Do 1 eight_reals_fft_2 macros */ \
/* Do 1 nop_two_two_complex_fft_2 macros */ \
/* Do 2 four_complex_fft macros */ \
/*	distance between fft data elements is 2K */ \
/*	do 1 macros each processing 32 data values */ \
\
start_timer(24); \
rdi = (uintptr_t)g->sincos1;/* Load sin/cos pointer */ \
for(unsigned int loopA = 2*clm; loopA; loopA--) { \
	g4cl_eight_reals_fft_2(rsi, 64, clmblkdst, 2*clmblkdst, rcx, 64, blkdst, 2*blkdst); \
} \
\
/* Do 92 four_complex_fft macros */ \
/*	distance between fft data elements is 2K */ \
/*	do 23 macros each processing 32 data values */ \
\
rsi += -2*clm*64+4*clmblkdst;/* Next source pointer */ \
rcx += -2*clm*64+4*blkdst;	/* Next dest pointer */ \
rdi += XMM_SCD;		/* Next sine/cosine pointer */ \
/* eax = 12*256+1; */ /* 23 iterations of 2*clm */ \
for(unsigned int loopA = 23; loopA; ) { \
	do{ \
		for(unsigned int loopC = clm; loopC; loopC--) { \
			for(unsigned int loopD = 2; loopD; loopD--) { \
				g4cl_four_complex_fft(rsi, 64, clmblkdst, 2*clmblkdst, rcx, 64, blkdst, 2*blkdst); \
				L2prefetch128(u8ptr(rcx-64+rbp)); \
				L2prefetch128(u8ptr(rcx-64+blkdst+rbp)); \
				pfing(rbp += -64+2*blkdst); \
			} \
			pfing(rbp += 128-4*blkdst); \
		} \
		rsi += -2*clm*64+4*clmblkdst;/* Next source pointer */ \
		rcx += -2*clm*64+4*blkdst;	/* Next dest pointer */ \
		rdi += 2*XMM_SCD;		/* Next sine/cosine pointer */ \
	}while(--loopA & 1); \
	rsi += -8*clmblkdst+clmblkdst8; \
} \
end_timer(24); \
\
/* Work on next 4 sets of 384 values */ \
\
endlp:	xpass1_get_next_block(pass2, c0b, b0b); /* Get next block and jump */ \
/* done: */	xfft_footer_hg2; \
\
\
\
/* Perform a 768K-element negacyclic FFT.  This is done in two passes.  Pass 1 */ \
/* does 9 levels, pass 2 uses common code to perform the last 11 levels. */ \
\

#define hg_pass1sclevels9complex3(procname, xpass2name, pass2_levels, clm) \
 /* (fftname, pass2_levels, clmarg) */ \
\
/* set_FFT_constants(fftname, 0, 1, pass2_levels, clmarg); */ \
xfft_header_hg2(procname, xpass2name, pass2_levels, clm); \
\
/* Do the last 11 FFT levels */ \
/*	do 96 pass2_11_levels_complex */ \
\
/* pass2: */	start_timer(0); \
xpass2name(g); \
end_timer(0); \
if(g->ffttype == 1)		/* We're done if FFTing only */ \
	return; \
\
/* Do the inverse FFT */ \
\
pass1_inverse_fft_setup(c0b, procname); \
\
/* Do inverse FFT levels 8,9 */ \
/* On input the 64-byte cache lines hold these data values: */ \
/*	0K	2K	4K	6K */ \
/*	8K	10K	12K	14K */ \
/*	... */ \
/* On output the 64-byte cache lines hold these data values: */ \
/*	0K	4K	8K	12K */ \
/*	2K	... */ \
/*	16K	... */ \
/*	... */ \
\
/* Do 96 four_complex_unfft macros */ \
/*	distance between fft data elements is 2K */ \
/*	do 48 macros each processing 16 data values */ \
\
c0b:	start_timer(16); \
rsi = (uintptr_t)g->data_addr;/* Get source address */ \
rcx = (uintptr_t)g->scratch_area;	/* Get address of scratch area */ \
rdi = (uintptr_t)g->sincos1;/* Load sin/cos pointer */ \
/* 12 iterations of 4 iters of 2*clm */ \
pfing(rbp = (uintptr_t)g->data_prefetch);/* Set prefetch distance */ \
pfing(rbp -= rsi); \
for(unsigned int loopA = 12; loopA; loopA--) { \
	for(unsigned int loopB = 4; loopB; loopB--) { \
		for(unsigned int loopC = 2*clm; loopC; loopC--) { \
			g2cl_four_complex_unfft(rsi, 64, blkdst, rcx, 64, clmblkdst); \
		} \
		rsi += -2*clm*64+2*blkdst;	/* Next source pointer */ \
		rcx += -2*clm*64+2*clmblkdst;/* Next dest pointer */ \
		rdi += XMM_SCD;		/* Next sine/cosine pointer */ \
	} \
	rcx += -8*clmblkdst+clmblkdst8; \
} \
rcx = rsi-48*2*blkdst;	/* Create prefetch pointer */ \
rsi = (uintptr_t)g->scratch_area;	/* Restore source pointer */ \
end_timer(16); \
\
/* Do inverse FFT levels 6,7 */ \
/* On input the 64-byte cache lines hold these data values: */ \
/*	0K	4K	8K	12K */ \
/*	2K	... */ \
/*	16K	... */ \
/*	... */ \
/* On output the 64-byte cache lines hold these data values: */ \
/*	0K	16K	32K	48K */ \
/*	2K	... */ \
/*	... */ \
/*	14K */ \
/*	64K	... */ \
/*	... */ \
\
/* Do 96 four_complex_unfft macros */ \
/*	distance between fft data elements is 8K */ \
/*	do 24 macros each processing 32 data values */ \
\
start_timer(17); \
rdi = (uintptr_t)g->sincos2;/* Load sin/cos pointer */ \
/* 12 iterations of 2*2*clm */ \
for(unsigned int loopA = 12; loopA; loopA--) { \
	xtouch(xptr(rcx+(clm-1)*128));	/* Preload TLB */ \
	xloop2_clm##clm(c5b, clmblkdst, x4cl_four_complex_unfft(rsi, 64, 2*clmblkdst, 4*clmblkdst, rdi));		/* 2 * 2*clm iterations */ \
	rsi += -2*clmblkdst+clmblkdst8;/* Next source pointer */ \
	rdi += XMM_SCD;		/* Next sine/cosine pointer */ \
} \
rsi += -12*clmblkdst8;	/* Restore source pointer */ \
end_timer(17); \
\
/* Do inverse FFT levels 4,5 */ \
/* On input the 64-byte cache lines hold these data values: */ \
/*	0K	16K	32K	48K */ \
/*	2K	... */ \
/*	... */ \
/*	14K */ \
/*	64K	... */ \
/*	... */ \
/* On output the 64-byte cache lines hold these data values: */ \
/*	0K	64K	128K	192K */ \
/*	2K	... */ \
/*	... */ \
/*	30K */ \
/*	256K	... */ \
/*	... */ \
\
/* Do 96 four_complex_unfft macros */ \
/*	distance between fft data elements is 32K */ \
/*	do 24 macros each processing 32 data values */ \
\
start_timer(18); \
rdi = (uintptr_t)g->sincos3;/* Load sin/cos pointer */ \
/* 3 iterations of 4 iters of 2*2*clm */ \
for(unsigned int loopA = 3; loopA; loopA--) { \
	for(unsigned int loopB = 4; loopB; loopB--) { \
		xtouch(xptr(rcx+(clm-1)*128));	/* Preload TLB */ \
		xloop2_clm##clm(c3b, clmblkdst, x4cl_four_complex_unfft(rsi, 64, clmblkdst8, 2*clmblkdst8, rdi));		/* 2* 2*clm iterations */ \
	} \
	rsi += -4*2*clmblkdst+4*clmblkdst8;/* Next source pointer */ \
	rdi += XMM_SCD;		/* Next sine/cosine pointer */ \
} \
rsi += -3*4*clmblkdst8;	/* Restore source pointer */ \
end_timer(18); \
\
/* Do inverse FFT levels 1,2,3 */ \
/* On input the 64-byte cache lines hold these data values: */ \
/*	0K	64K	128K	192K */ \
/*	2K	... */ \
/*	... */ \
/*	30K */ \
/*	256K	... */ \
/*	... */ \
/* On output the 64-byte cache lines hold these data values: */ \
/*	0K	192K	384K	576K */ \
/*	... */ \
/*	190K */ \
\
/* Do 128 three_complex_last_unfft macros */ \
/*	distance between fft data elements is 128K */ \
/*	do 32 macros each processing 24 data values */ \
\
start_timer(19); \
rdi = (uintptr_t)g->plus1_premults;/* Address of premultiplier table */ \
/* 4 iters of 8 iters of 2*clm */ \
for(unsigned int loopA = 4; loopA; loopA--) { \
	for(unsigned int loopB = 8; loopB; loopB--) { \
		xtouch(xptr(rcx+(clm-1)*128));	/* Preload TLB */ \
		xloop1_clm##clm(c1b, clmblkdst, x3cl_three_complex_last_unfft(rsi, 64, 4*clmblkdst8));		/* 2*clm iterations */ \
		rdi += 192;		/* Next premultiplier pointer */ \
	} \
	rsi += -8*clmblkdst+clmblkdst8;/* Next source pointer */ \
} \
end_timer(19); \
\
/* Normalize these values */ \
\
pass1_normalize(1, yes_fft); \
\
/* Copy data from the scratch area back to the FFT data area */ \
\
copy_scratch_data_to_fft(96, clm); \
goto	endlp;			/* Skip FFT code if flag not set */ \
\
/* Fire up auxiliary threads to do pass 1 forward FFT */ \
\
/* pass1: */	pass1_forward_fft_setup(b0b); \
\
/* Copy data from the FFT data area back to the scratch area */ \
\
b0b:	copy_fft_data_to_scratch(96, clm); \
\
/* Do FFT levels 1,2,3 */ \
/* */ \
/* On input the 64-byte cache lines hold these data values: */ \
/*	0K	192K	384K	576K */ \
/*	... */ \
/*	190K */ \
/* On output the 64-byte cache lines hold these data values: */ \
/*	0K	64K	128K	192K */ \
/*	2K	... */ \
/*	... */ \
/*	62K */ \
/*	256K	384K	512K	640K */ \
/*	258K */ \
/*	... */ \
/*	382K */ \
\
/* Do 128 three_complex_first_fft macros (FFT levels 1,2,3) */ \
/*	distance between fft data elements is 128K */ \
/*	do 32 iterations each processing 24 data values */ \
\
yes_fft:start_timer(21); \
rsi = (uintptr_t)g->scratch_area;	/* Get address of scratch area */ \
rdi = (uintptr_t)g->plus1_premults;/* Address of premultiplier table */ \
rcx = (uintptr_t)g->data_addr; /* Load source address for prefetching */ \
pfing(rbp = (uintptr_t)g->data_prefetch);/* Set prefetch distance */ \
pfing(rbp -= rcx); \
pfing(rcx += 56*blkdst);		/* Create prefetching pointer */ \
/* 4 iters of 4 iters of 2*2*clm */ \
for(unsigned int loopA = 4; loopA; loopA--) { \
	for(unsigned int loopB = 4; loopB; loopB--) { \
		xtouch(xptr(rcx+(clm-1)*128));	/* Preload TLB */ \
		xloop2edi_clm##clm(b1b, clmblkdst, 192, x3cl_three_complex_fft(rsi, 64, 4*clmblkdst8)); /* 2 * 2*clm iterations */ \
	} \
	rsi += -8*clmblkdst+clmblkdst8;/* Next source pointer */ \
} \
rsi += -4*clmblkdst8;	/* Restore source pointer */ \
end_timer(21); \
\
/* Do FFT levels 4,5 */ \
/* */ \
/* On input the 64-byte cache lines hold these data values: */ \
/*	0K	64K	128K	192K */ \
/*	2K	... */ \
/*	... */ \
/*	62K */ \
/*	256K	384K	512K	640K */ \
/*	258K */ \
/*	... */ \
/*	382K */ \
/* On output the 64-byte cache lines hold these data values: */ \
/*	0K	32K	64K	96K */ \
/*	2K	... */ \
/*	... */ \
/*	30K	... */ \
/*	128K	... */ \
/*	... */ \
\
/* Do 96 four_complex_fft macros */ \
/*	distance between fft data elements is 32K */ \
/*	do 8 macros each processing 16 data values */ \
/*	do 16 macros each processing 32 data values */ \
\
start_timer(22); \
rdi = (uintptr_t)g->sincos3;/* Load sin/cos pointer */ \
/* 2 iters of 2 of 4*2*clm */ \
for(unsigned int loopA = 2; loopA; loopA--) { \
	for(unsigned int loopB = 2; loopB; loopB--) { \
		xtouch(xptr(rcx+(clm-1)*128));	/* Preload TLB */ \
		xloop4_clm##clm(b2b, clmblkdst, x2cl_four_complex_fft(rsi, 64, 2*clmblkdst8));		/* 4 * 2*clm iterations */ \
	} \
	rsi += -2*4*clmblkdst+clmblkdst8; \
} \
rsi += -2*clmblkdst8+4*clmblkdst8;/* Next source pointer */ \
rdi += XMM_SCD;		/* Next sine/cosine pointer */ \
/* 2 iters of 4 of 2*2*clm */ \
for(unsigned int loopA = 2; loopA; loopA--) { \
	for(unsigned int loopB = 4; loopB; loopB--) { \
		xtouch(xptr(rcx+(clm-1)*128));	/* Preload TLB */ \
		xloop2_clm##clm(b3b, clmblkdst, x4cl_four_complex_fft(rsi, 64, 2*clmblkdst8, 4*clmblkdst8, rdi));		/* 2 * 2*clm iterations */ \
	} \
	rsi += -4*2*clmblkdst+clmblkdst8; \
} \
rsi += -2*clmblkdst8-4*clmblkdst8;/* Next source pointer */ \
end_timer(22); \
\
/* Do FFT levels 6,7 */ \
/* */ \
/* On input the 64-byte cache lines hold these data values: */ \
/*	0K	32K	64K	96K */ \
/*	2K	... */ \
/*	... */ \
/*	30K	... */ \
/*	128K	... */ \
/*	... */ \
/* On output the 64-byte cache lines hold these data values: */ \
/*	0K	8K	16K	24K */ \
/*	2K	... */ \
/*	... */ \
/*	6K	... */ \
/*	32K	... */ \
/*	... */ \
\
/* Do 96 four_complex_fft macros */ \
/*	distance between fft data elements is 8K */ \
/*	do 24 macros each processing 32 data values */ \
\
start_timer(23); \
rdi = (uintptr_t)g->sincos2;/* Load sin/cos pointer */ \
/* 6 iters of 2 iters of 2*2*clm */ \
for(unsigned int loopA = 6; loopA; loopA--) { \
	for(unsigned int loopB = 2; loopB; loopB--) { \
		xtouch(xptr(rcx+(clm-1)*128));	/* Preload TLB */ \
		xloop2_clm##clm(b5b, clmblkdst, x4cl_four_complex_fft(rsi, 64, 4*clmblkdst, clmblkdst8, rdi));		/* 2 * 2*clm iterations */ \
	} \
	rsi += -4*clmblkdst+2*clmblkdst8;/* Next source pointer */ \
	rdi += 2*XMM_SCD;		/* Next sine/cosine pointer */ \
} \
rsi += -6*2*clmblkdst8;	/* Restore source pointer */ \
pfing(rcx += -96*blkdst);		/* Restore prefetch pointer */ \
end_timer(23); \
\
/* Do FFT levels 8,9 */ \
/* */ \
/* On input the 64-byte cache lines hold these data values: */ \
/*	0K	8K	16K	24K */ \
/*	2K	... */ \
/*	... */ \
/*	6K	... */ \
/*	32K	... */ \
/*	... */ \
/* On output the 64-byte cache lines hold these data values: */ \
/*	0K	2K	4K	6K */ \
/*	8K	... */ \
/*	... */ \
\
/* Do 96 four_complex_fft macros */ \
/*	distance between fft data elements is 2K */ \
/*	do 24 macros each processing 32 data values */ \
\
start_timer(24); \
rdi = (uintptr_t)g->sincos1;/* Load sin/cos pointer */ \
/* 12 iterations 2 iters of 2*clm */ \
for(unsigned int loopA = 12; loopA; loopA--) { \
	for(unsigned int loopB = 2; loopB; loopB--) { \
		for(unsigned int loopC = clm; loopC; loopC--) { \
			for(unsigned int loopD = 2; loopD; loopD--) { \
				g4cl_four_complex_fft(rsi, 64, clmblkdst, 2*clmblkdst, rcx, 64, blkdst, 2*blkdst); \
				L2prefetch128(u8ptr(rcx-64+rbp)); \
				L2prefetch128(u8ptr(rcx-64+blkdst+rbp)); \
				pfing(rbp += -64+2*blkdst); \
			} \
			pfing(rbp += 128-4*blkdst); \
		} \
	rsi += -2*clm*64+4*clmblkdst;/* Next source pointer */ \
		rcx += -2*clm*64+4*blkdst;	/* Next dest pointer */ \
		rdi += 2*XMM_SCD;		/* Next sine/cosine pointer */ \
	} \
	rsi += -8*clmblkdst+clmblkdst8; \
} \
end_timer(24); \
\
/* Work on next 4 sets of 384 values */ \
\
endlp:	xpass1_get_next_block(pass2, c0b, b0b); /* Get next block and jump */ \
/* done: */	xfft_footer_hg2; \
\
\
\
/* Perform a 896K-element FFT.  This is done in two passes.  Pass 1 */ \
/* does 9 levels, pass 2 uses common code to perform the last 11 levels. */ \
\

#define hg_pass1sclevels9pfa7(procname, xpass2name, pass2_levels, clm) \
 /* (fftname, pass2_levels, clmarg) */ \
\
/* set_FFT_constants(fftname, 0, 0, pass2_levels, clmarg); */ \
xfft_header_hg2(procname, xpass2name, pass2_levels, clm); \
\
/* Do the last 11 FFT levels */ \
/*	do 1 pass2_11_levels_real and 111 pass2_11_levels_complex */ \
\
/* pass2: */	start_timer(0); \
xpass2name(g); \
end_timer(0); \
if(g->ffttype == 1)		/* We're done if FFTing only */ \
	return; \
\
/* Do the inverse FFT */ \
\
pass1_inverse_fft_setup(c0b, procname); \
\
/* Do inverse FFT levels 8,9 */ \
/* On input the 64-byte cache lines hold these data values: */ \
/*	0K	2K	4K	6K */ \
/*	8K	10K	12K	14K */ \
/*	... */ \
/* On output the 64-byte cache lines hold these data values: */ \
/*	0K	4K	8K	12K */ \
/*	2K	... */ \
/*	16K	... */ \
/*	... */ \
\
/* Do 1 eight_reals_unfft_2 macros */ \
/* Do 1 nop_two_two_complex_unfft_2 macros */ \
/*	distance between fft data elements is 2K */ \
/*	do 1 macros each processing 16 data values */ \
\
c0b:	start_timer(16); \
rsi = (uintptr_t)g->data_addr;/* Get source address */ \
rcx = (uintptr_t)g->scratch_area;	/* Get address of scratch area */ \
pfing(rbp = (uintptr_t)g->data_prefetch);/* Set prefetch distance */ \
pfing(rbp -= rsi); \
for(unsigned int loopA = 2*clm; loopA; loopA--) { \
	g2cl_eight_reals_unfft_2(rsi, 64, blkdst, rcx, 64, clmblkdst); \
} \
\
/* Do 110 four_complex_unfft macros */ \
/*	distance between fft data elements is 2K */ \
/*	do 55 macros each processing 16 data values */ \
\
rsi += -2*clm*64+2*blkdst;	/* Next source pointer */ \
rcx += -2*clm*64+2*clmblkdst;/* Next dest pointer */ \
rdi = (uintptr_t)g->sincos1;/* Load sin/cos pointer */ \
/* al = 256/4+14; */ /* 55 iterations of 2*clm */ \
for(unsigned int loopA = 55; loopA; ) { \
	do{ \
		for(unsigned int loopC = 2*clm; loopC; loopC--) { \
			g2cl_four_complex_unfft(rsi, 64, blkdst, rcx, 64, clmblkdst); \
		} \
		rsi += -2*clm*64+2*blkdst;	/* Next source pointer */ \
		rcx += -2*clm*64+2*clmblkdst;/* Next dest pointer */ \
		rdi += XMM_SCD;		/* Next sine/cosine pointer */ \
	}while(--loopA & 3); \
	rcx += -8*clmblkdst+clmblkdst8; \
} \
rcx = rsi-56*2*blkdst;	/* Create prefetch pointer */ \
rsi = (uintptr_t)g->scratch_area;	/* Restore dest pointer */ \
end_timer(16); \
\
/* Do inverse FFT levels 6,7 */ \
/* On input the 64-byte cache lines hold these data values: */ \
/*	0K	4K	8K	12K */ \
/*	2K	... */ \
/*	16K	... */ \
/*	... */ \
/* On output the 64-byte cache lines hold these data values: */ \
/*	0K	16K	32K	48K */ \
/*	2K	... */ \
/*	... */ \
/*	14K */ \
/*	64K	... */ \
/*	... */ \
\
/* Do 4 eight_reals_unfft_2 macros */ \
/* Do 4 nop_two_two_complex_unfft_2 macros */ \
/*	distance between fft data elements is 8K */ \
/*	do 2 macros each processing 32 data values */ \
\
start_timer(17); \
/* 2 iterations of 2*clm */ \
for(unsigned int loopA = 2; loopA; loopA--) { \
	for(unsigned int loopB = 2*clm; loopB; loopB--) { \
		x4cl_eight_reals_unfft_2(rsi, 64, 2*clmblkdst, 4*clmblkdst); \
	} \
	rsi += -2*clm*64+clmblkdst;/* Next source pointer */ \
} \
\
/* Do 104 four_complex_unfft macros */ \
/*	distance between fft data elements is 8K */ \
/*	do 26 macros each processing 32 data values */ \
\
rsi += -2*clmblkdst+clmblkdst8;/* Next source pointer */ \
rdi = (uintptr_t)g->sincos2;/* Load sin/cos pointer */ \
/* 13 iterations of 2 iters of 2*clm */ \
for(unsigned int loopA = 13; loopA; loopA--) { \
	for(unsigned int loopB = 2; loopB; loopB--) { \
		xtouch(xptr(rcx+(clm-1)*128));	/* Preload TLB */ \
		xloop1_clm##clm(c5b, clmblkdst, x4cl_four_complex_unfft(rsi, 64, 2*clmblkdst, 4*clmblkdst, rdi)); \
	} \
	rsi += -2*clmblkdst+clmblkdst8;/* Next source pointer */ \
	rdi += XMM_SCD;		/* Next sine/cosine pointer */ \
} \
rsi += -14*clmblkdst8;	/* Restore source pointer */ \
end_timer(17); \
\
/* Do inverse FFT levels 4,5 */ \
/* On input the 64-byte cache lines hold these data values: */ \
/*	0K	16K	32K	48K */ \
/*	2K	... */ \
/*	... */ \
/*	14K */ \
/*	64K	... */ \
/*	... */ \
/* On output the 64-byte cache lines hold these data values: */ \
/*	0K	32K	64K	96K */ \
/*	2K	... */ \
/*	... */ \
/*	30K */ \
/*	128K	192K	256K	320K */ \
/*	... */ \
/* */ \
/* Do 16 eight_reals_unfft_2 macros */ \
/*	distance between fft data elements is 32K */ \
/*	do 8 macros each processing 16 data values */ \
/* */ \
/* Do 96 four_complex_unfft macros */ \
/*	distance between fft data elements is 32K */ \
/*	do 24 macros each processing 32 data values */ \
\
/* Do inverse FFT levels 1,2,3 */ \
/* On input the 64-byte cache lines hold these data values: */ \
/*	0K	32K	64K	96K */ \
/*	2K	... */ \
/*	... */ \
/*	30K */ \
/*	128K	192K	256K	320K */ \
/*	... */ \
/* On output the 64-byte cache lines hold these data values: */ \
/*	0K	32K	448K	480K */ \
/*	2K	... */ \
/*	... */ \
/*	30K */ \
/*	64K	128K	512K	576K */ \
/*	66K */ \
/*	... */ \
/*	126K */ \
/*	192K	320K	640K	768K */ \
/*	194K */ \
/*	... */ \
/*	318K */ \
/* */ \
/* Do 128 seven_reals_last_unfft macros */ \
/*	distance between fft data elements is 128K */ \
/*	do 32 macros each processing 56 data values */ \
\
start_timer(18); \
rdi = (uintptr_t)g->sincos3;/* Load sin/cos pointer */ \
pfing(rcx += -4*blkdst);		/* Adjust prefetch pointer */ \
/* 8 iterations */ \
xseven_reals_unfft_preload; \
\
for(unsigned int loopA = 8; loopA; loopA--) { \
	xtouchlines(5);			/* Preload TLBs */ \
	\
	for(unsigned int loopB = clm; loopB; loopB--) { \
		for(unsigned int loopC = 2; loopC; loopC--) { \
			x2cl_half_eight_reals_unfft_2(rsi, 64, clmblkdst8); \
		} \
		\
		rsi += -2*64+2*clmblkdst8;	/* Next source pointer */ \
		for(unsigned int loopC = 3; loopC; loopC--) { \
			xloop1sp(c3b, 4*clmblkdst8, x4cl_four_complex_unfft(rsi, 64, clmblkdst8, 2*clmblkdst8, rdi)); \
			rdi += XMM_SCD;		/* Next sine/cosine pointer */ \
		} \
		rdi += -3*XMM_SCD;		/* Restore sine/cosine pointer */ \
		\
		rsi += -3*4*clmblkdst8-2*clmblkdst8;/* Restore source pointer */ \
		for(unsigned int loopC = 2; loopC; loopC--) { \
			xloop1sp(c1b, clmblkdst8, x7cl_seven_reals_last_unfft(rsi, 64, 2*clmblkdst8)); \
		} \
		\
		EXEC248##clm(rsi += -2*clmblkdst8+2*64);	/* Next source pointer */ \
		EXEC248##clm(pfing(rcx += -5*blkdst+128));	/* Next prefetch pointer */ \
	} \
	EXEC248##clm(rsi += -clm*2*64+clmblkdst);/* Next source pointer */ \
	EXEC248##clm(pfing(rcx += -clm*128+5*blkdst));	/* Next prefetch pointer */ \
	EXEC1##clm(rsi += -2*clmblkdst8+clmblkdst);/* Next source pointer */ \
	\
} \
end_timer(18); \
\
/* Normalize these values */ \
\
pass1_normalize(1, yes_fft); \
\
/* Copy data from the scratch area back to the FFT data area */ \
\
copy_scratch_data_to_fft(112, clm); \
goto	endlp;			/* Skip FFT code if flag not set */ \
\
/* Fire up auxiliary threads to do pass 1 forward FFT */ \
\
/* pass1: */	pass1_forward_fft_setup(b0b); \
\
/* Copy data from the FFT data area back to the scratch area */ \
\
b0b:	copy_fft_data_to_scratch(112, clm); \
\
/* Do FFT levels 1,2,3 */ \
/* */ \
/* On input the 64-byte cache lines hold these data values: */ \
/*	0K	32K	448K	480K */ \
/*	2K	... */ \
/*	... */ \
/*	30K */ \
/*	64K	128K	512K	576K */ \
/*	66K */ \
/*	... */ \
/*	126K */ \
/*	192K	320K	640K	768K */ \
/*	194K */ \
/*	... */ \
/*	318K */ \
/* On output the 64-byte cache lines hold these data values: */ \
/*	0K	32K	64K	96K */ \
/*	2K	... */ \
/*	... */ \
/*	30K */ \
/*	128K	192K	256K	320K */ \
/*	130K */ \
/*	... */ \
/*	190K */ \
/*	384K	512K	640K	768K */ \
/*	386K */ \
/*	... */ \
/*	510K */ \
/* */ \
/* Do 128 seven_reals_first_fft macros (FFT levels 1,2,3) */ \
/*	distance between fft data elements is 128K */ \
/*	do 16 iterations each processing 56 data values */ \
\
/* Do FFT levels 4,5 */ \
/* Values 0K-127K is real data, 128K-895K is complex data. */ \
/* */ \
/* On input the 64-byte cache lines hold these data values: */ \
/*	0K	32K	64K	96K */ \
/*	2K	... */ \
/*	... */ \
/*	30K */ \
/*	128K	192K	256K	320K */ \
/*	130K */ \
/*	... */ \
/*	190K */ \
/*	384K	512K	640K	768K */ \
/*	386K */ \
/*	... */ \
/*	510K */ \
/* On output the 64-byte cache lines hold these data values: */ \
/*	0K	32K	64K	96K */ \
/*	2K	... */ \
/*	... */ \
/*	30K	... */ \
/*	128K	... */ \
/*	... */ \
/* */ \
/* Do 16 eight_reals_fft_2 macros */ \
/*	distance between fft data elements is 32K */ \
/*	do 16 macros each processing 8 data values */ \
/* */ \
/* Do 64 four_complex_fft macros */ \
/*	distance between fft data elements is 32K */ \
/*	do 16 macros each processing 32 data values */ \
\
/* Do FFT levels 6,7 */ \
/* Values 0K-31K is real data, 32K-63K is semi-real data, 64K-895K is */ \
/* complex data. */ \
/* */ \
/* On input the 64-byte cache lines hold these data values: */ \
/*	0K	32K	64K	96K */ \
/*	2K	... */ \
/*	... */ \
/*	30K	... */ \
/*	128K	... */ \
/*	... */ \
/* On output the 64-byte cache lines hold these data values: */ \
/*	0K	8K	16K	24K */ \
/*	2K	... */ \
/*	... */ \
/*	6K	... */ \
/*	32K	... */ \
/*	... */ \
/* */ \
/* Do 8 eight_reals_fft_2 macros */ \
/* Do 8 two_two_complex_fft macros */ \
/*	distance between fft data elements is 8K */ \
/*	do 4 macros each processing 32 data values */ \
/* */ \
/* Do 96 four_complex_fft macros */ \
/*	distance between fft data elements is 8K */ \
/*	do 24 macros each processing 32 data values */ \
\
yes_fft:start_timer(21); \
rsi = (uintptr_t)g->scratch_area;	/* Get address of scratch area */ \
rcx = (uintptr_t)g->data_addr; /* Load source address for prefetching */ \
pfing(rbp = (uintptr_t)g->data_prefetch);/* Set prefetch distance */ \
pfing(rbp -= rcx); \
pfing(rcx += 56*blkdst);		/* Create prefetching pointer */ \
rdi = (uintptr_t)g->sincos3;/* Load sin/cos pointer */ \
rbx = (uintptr_t)g->sincos2;/* Load sin/cos pointer */ \
/* 4 iterations */ \
xseven_reals_fft_preload; \
\
for(unsigned int loopA = 4; loopA; loopA--) { \
	xtouchlines(14);			/* Preload TLBs */ \
	\
	for(unsigned int loopB = clm; loopB; loopB--) { \
		for(unsigned int loopC = 2; loopC; loopC--) { \
			for(unsigned int loopD = 2; loopD; loopD--) { \
				xloop1sp(b1b, 4*clmblkdst, x7cl_seven_reals_fft(rsi, 64, 2*clmblkdst8)); \
			} \
			rsi += -2*4*clmblkdst+clmblkdst8; \
		} \
		\
		rsi += -2*clmblkdst8;	/* Restore source pointer */ \
		for(unsigned int loopC = 2; loopC; loopC--) { \
			for(unsigned int loopD = 2; loopD; loopD--) { \
				for(unsigned int loopE = 2; loopE; loopE--) { \
					x1cl_half_eight_reals_fft_2(rsi, 64); \
				} \
				rsi += -2*64+4*clmblkdst; \
			} \
			rsi += -2*4*clmblkdst+clmblkdst8; \
		} \
		\
		/*	bump	rsi, -2*clmblkdst8+2*clmblkdst8;; Next source pointer */ \
		for(unsigned int loopC = 2; loopC; loopC--) { \
			for(unsigned int loopD = 2; loopD; loopD--) { \
				for(unsigned int loopE = 2; loopE; loopE--) { \
					x2cl_four_complex_fft(rsi, 64, 2*clmblkdst8); \
				} \
				rsi += -2*64+4*clmblkdst; \
			} \
			rsi += -2*4*clmblkdst+clmblkdst8; \
		} \
		\
		rsi += -2*clmblkdst8+4*clmblkdst8;/* Next source pointer */ \
		rdi += XMM_SCD;		/* Next sine/cosine pointer */ \
		for(unsigned int loopC = 2; loopC; loopC--) { \
			for(unsigned int loopD = 2; loopD; loopD--) { \
				xloop1sp(b3b, 4*clmblkdst, x4cl_four_complex_fft(rsi, 64, 2*clmblkdst8, 4*clmblkdst8, rdi)); \
			} \
			rsi += -2*4*clmblkdst+clmblkdst8; \
		} \
		rdi += -XMM_SCD;		/* Restore sine/cosine pointer */ \
		\
		rsi += -2*clmblkdst8-6*clmblkdst8;/* Restore source pointer */ \
		for(unsigned int loopC = 2; loopC; loopC--) { \
			x4cl_eight_reals_fft_2(rsi, 64, 4*clmblkdst, clmblkdst8, rbx); \
		} \
		\
		rsi += -2*64+2*clmblkdst8;	/* Load source pointer */ \
		rbx += XMM_SCD;		/* Next sine/cosine pointer */ \
		/* 6 iterations */ \
		for(unsigned int loopC = 6; loopC; loopC--) { \
			xloop1sp(b5b, 2*clmblkdst8, x4cl_four_complex_fft(rsi, 64, 4*clmblkdst, clmblkdst8, rbx)); \
			rbx += 2*XMM_SCD;		/* Next sine/cosine pointer */ \
		} \
		rbx += -13*XMM_SCD;	/* Restore sine/cosine pointer */ \
		\
		EXEC248##clm(rsi += -7*2*clmblkdst8+2*64);/* Next source pointer */ \
		EXEC248##clm(pfing(rcx += -14*blkdst+128));	/* Next prefetch pointer */ \
	} \
	EXEC248##clm(rsi += -clm*2*64+clmblkdst);/* Next source pointer */ \
	EXEC248##clm(pfing(rcx += -clm*128+14*blkdst));	/* Next prefetch pointer */ \
	EXEC1##clm(rsi += -7*2*clmblkdst8+clmblkdst);/* Next source pointer */ \
	\
} \
rsi += -4*clmblkdst;	/* Restore source pointer */ \
pfing(rcx += -112*blkdst);	/* Restore prefetch pointer */ \
end_timer(21); \
\
/* Do FFT levels 8,9 */ \
/* Values 0K-7K is real data, 8K-15K is semi-real data, 16K-895K is */ \
/* complex data. */ \
/* */ \
/* On input the 64-byte cache lines hold these data values: */ \
/*	0K	8K	16K	24K */ \
/*	2K	... */ \
/*	... */ \
/*	6K	... */ \
/*	32K	... */ \
/*	... */ \
/* On output the 64-byte cache lines hold these data values: */ \
/*	0K	2K	4K	6K */ \
/*	8K	... */ \
/*	... */ \
\
/* Do 1 eight_reals_fft_2 macros */ \
/* Do 1 nop_two_two_complex_fft_2 macros */ \
/* Do 2 four_complex_fft macros */ \
/*	distance between fft data elements is 2K */ \
/*	do 1 macros each processing 32 data values */ \
\
start_timer(24); \
rdi = (uintptr_t)g->sincos1;/* Load sin/cos pointer */ \
for(unsigned int loopA = 2*clm; loopA; loopA--) { \
	g4cl_eight_reals_fft_2(rsi, 64, clmblkdst, 2*clmblkdst, rcx, 64, blkdst, 2*blkdst); \
} \
\
/* Do 108 four_complex_fft macros */ \
/*	distance between fft data elements is 2K */ \
/*	do 27 macros each processing 32 data values */ \
\
rsi += -2*clm*64+4*clmblkdst;/* Next source pointer */ \
rcx += -2*clm*64+4*blkdst;	/* Next dest pointer */ \
rdi += XMM_SCD;		/* Next sine/cosine pointer */ \
/* eax = 14*256+1; */ /* 27 iterations of 2*clm */ \
for(unsigned int loopA = 27; loopA; ) { \
	do{ \
		for(unsigned int loopC = clm; loopC; loopC--) { \
			for(unsigned int loopD = 2; loopD; loopD--) { \
				g4cl_four_complex_fft(rsi, 64, clmblkdst, 2*clmblkdst, rcx, 64, blkdst, 2*blkdst); \
				L2prefetch128(u8ptr(rcx-64+rbp)); \
				L2prefetch128(u8ptr(rcx-64+blkdst+rbp)); \
				pfing(rbp += -64+2*blkdst); \
			} \
			pfing(rbp += 128-4*blkdst); \
		} \
		rsi += -2*clm*64+4*clmblkdst;/* Next source pointer */ \
		rcx += -2*clm*64+4*blkdst;	/* Next dest pointer */ \
		rdi += 2*XMM_SCD;		/* Next sine/cosine pointer */ \
	}while(--loopA & 1); \
	rsi += -8*clmblkdst+clmblkdst8; \
} \
end_timer(24); \
\
/* Work on next 4 sets of 448 values */ \
\
endlp:	xpass1_get_next_block(pass2, c0b, b0b); /* Get next block and jump */ \
/* done: */	xfft_footer_hg2; \
\
\
\
/* Perform a 1M-element FFT.  This is done in two passes.  Pass 1 */ \
/* does 9 levels, pass 2 uses common code to perform the last 11 levels. */ \
\

#define hg_pass1sclevels9pfa8(procname, xpass2name, pass2_levels, clm) \
 /* (fftname, pass2_levels, clmarg) */ \
\
/* set_FFT_constants(fftname, 0, 0, pass2_levels, clmarg); */ \
xfft_header_hg2(procname, xpass2name, pass2_levels, clm); \
\
/* Do the last 11 FFT levels */ \
/*	do 1 pass2_11_levels_real and 127 pass2_11_levels_complex */ \
\
/* pass2: */	start_timer(0); \
xpass2name(g); \
end_timer(0); \
if(g->ffttype == 1)		/* We're done if FFTing only */ \
	return; \
\
/* Do the inverse FFT */ \
\
pass1_inverse_fft_setup(c0b, procname); \
\
/* Do inverse FFT levels 8,9 */ \
/* On input the 64-byte cache lines hold these data values: */ \
/*	0K	2K	4K	6K */ \
/*	8K	10K	12K	14K */ \
/*	... */ \
/* On output the 64-byte cache lines hold these data values: */ \
/*	0K	4K	8K	12K */ \
/*	2K	... */ \
/*	16K	... */ \
/*	... */ \
\
/* Do 1 eight_reals_unfft_2 macros */ \
/* Do 1 nop_two_two_complex_unfft_2 macros */ \
/*	distance between fft data elements is 2K */ \
/*	do 1 macros each processing 16 data values */ \
\
c0b:	start_timer(16); \
rsi = (uintptr_t)g->data_addr;/* Get source address */ \
rcx = (uintptr_t)g->scratch_area;	/* Get address of scratch area */ \
pfing(rbp = (uintptr_t)g->data_prefetch);/* Set prefetch distance */ \
pfing(rbp -= rsi); \
for(unsigned int loopA = 2*clm; loopA; loopA--) { \
	g2cl_eight_reals_unfft_2(rsi, 64, blkdst, rcx, 64, clmblkdst); \
} \
\
/* Do 126 four_complex_unfft macros */ \
/*	distance between fft data elements is 2K */ \
/*	do 63 macros each processing 16 data values */ \
\
rsi += -2*clm*64+2*blkdst;	/* Next source pointer */ \
rcx += -2*clm*64+2*clmblkdst;/* Next dest pointer */ \
rdi = (uintptr_t)g->sincos1;/* Load sin/cos pointer */ \
/* al = 256/4+16; */ /* 63 iterations of 2*clm */ \
for(unsigned int loopA = 63; loopA; ) { \
	do{ \
		for(unsigned int loopC = 2*clm; loopC; loopC--) { \
			g2cl_four_complex_unfft(rsi, 64, blkdst, rcx, 64, clmblkdst); \
		} \
		rsi += -2*clm*64+2*blkdst;	/* Next source pointer */ \
		rcx += -2*clm*64+2*clmblkdst;/* Next dest pointer */ \
		rdi += XMM_SCD;		/* Next sine/cosine pointer */ \
	}while(--loopA & 3); \
	rcx += -8*clmblkdst+clmblkdst8;/* Next dest pointer */ \
} \
rcx = rsi-64*2*blkdst;	/* Create prefetch pointer */ \
rsi = (uintptr_t)g->scratch_area;	/* Restore source pointer */ \
end_timer(16); \
\
/* Do inverse FFT levels 6,7 */ \
/* On input the 64-byte cache lines hold these data values: */ \
/*	0K	4K	8K	12K */ \
/*	2K	... */ \
/*	16K	... */ \
/*	... */ \
/* On output the 64-byte cache lines hold these data values: */ \
/*	0K	16K	32K	48K */ \
/*	2K	... */ \
/*	... */ \
/*	14K */ \
/*	64K	... */ \
/*	... */ \
\
/* Do 4 eight_reals_unfft_2 macros */ \
/* Do 4 nop_two_two_complex_unfft_2 macros */ \
/*	distance between fft data elements is 8K */ \
/*	do 2 macros each processing 32 data values */ \
\
start_timer(17); \
/* 2 iterations of 2*clm */ \
for(unsigned int loopA = 2; loopA; loopA--) { \
	xtouch(xptr(rcx+(clm-1)*128));	/* Preload TLB */ \
	xloop1_clm##clm(c4b, clmblkdst, x4cl_eight_reals_unfft_2(rsi, 64, 2*clmblkdst, 4*clmblkdst)); \
} \
\
/* Do 120 four_complex_unfft macros */ \
/*	distance between fft data elements is 8K */ \
/*	do 30 macros each processing 32 data values */ \
\
rsi += -2*clmblkdst+clmblkdst8;/* Next source pointer */ \
rdi = (uintptr_t)g->sincos2;/* Load sin/cos pointer */ \
/* 15 iterations of 2 iters of 2*clm */ \
for(unsigned int loopA = 15; loopA; loopA--) { \
	for(unsigned int loopB = 2; loopB; loopB--) { \
		xtouch(xptr(rcx+(clm-1)*128));	/* Preload TLB */ \
		xloop1_clm##clm(c5b, clmblkdst, x4cl_four_complex_unfft(rsi, 64, 2*clmblkdst, 4*clmblkdst, rdi)); \
	} \
	rsi += -2*clmblkdst+clmblkdst8;/* Next source pointer */ \
	rdi += XMM_SCD;		/* Next sine/cosine pointer */ \
} \
rsi += -16*clmblkdst8;	/* Restore source pointer */ \
end_timer(17); \
\
/* Do inverse FFT levels 4,5 */ \
/* On input the 64-byte cache lines hold these data values: */ \
/*	0K	16K	32K	48K */ \
/*	2K	... */ \
/*	... */ \
/*	14K */ \
/*	64K	... */ \
/*	... */ \
/* On output the 64-byte cache lines hold these data values: */ \
/*	0K	64K	128K	192K */ \
/*	2K	... */ \
/*	... */ \
/*	62K */ \
/*	256K	... */ \
/*	... */ \
/* */ \
/* Do 16 eight_reals_unfft_2 macros */ \
/* Do 16 nop_two_two_complex_unfft_2 macros */ \
/*	distance between fft data elements is 32K */ \
/*	do 8 macros each processing 32 data values */ \
/* */ \
/* Do 96 four_complex_unfft macros */ \
/*	distance between fft data elements is 32K */ \
/*	do 24 macros each processing 32 data values */ \
\
/* Do inverse FFT levels 1,2,3 */ \
/* On input the 64-byte cache lines hold these data values: */ \
/*	0K	64K	128K	192K */ \
/*	2K	... */ \
/*	... */ \
/*	62K */ \
/*	256K	... */ \
/*	... */ \
/* On output the 64-byte cache lines hold these data values: */ \
/*	0K	256K	512K	768K */ \
/*	2K	... */ \
/*	... */ \
/*	254K */ \
/* */ \
/* Do 128 eight_reals_last_unfft macros */ \
/*	distance between fft data elements is 128K */ \
/*	do 32 macros each processing 32 data values */ \
\
start_timer(18); \
rdi = (uintptr_t)g->sincos3;/* Load sin/cos pointer */ \
/* 8 iterations */ \
\
for(unsigned int loopA = 8; loopA; loopA--) { \
	xtouchlines(4);			/* Preload TLBs */ \
	\
	for(unsigned int loopB = clm; loopB; loopB--) { \
		xloop1sp(c2b, 4*clmblkdst8, x4cl_eight_reals_unfft_2(rsi, 64, clmblkdst8, 2*clmblkdst8)); \
		\
		/* 3 iterations */ \
		for(unsigned int loopC = 3; loopC; loopC--) { \
			xloop1sp(c3b, 4*clmblkdst8, x4cl_four_complex_unfft(rsi, 64, clmblkdst8, 2*clmblkdst8, rdi)); \
			rdi += XMM_SCD;		/* Next sine/cosine pointer */ \
		} \
		rdi += -3*XMM_SCD;		/* Next sine/cosine pointer */ \
		\
		rsi += -4*4*clmblkdst8;	/* Next source pointer */ \
		for(unsigned int loopC = 4; loopC; loopC--) { \
			for(unsigned int loopD = 2; loopD; loopD--) { \
				x4cl_eight_reals_last_unfft(rsi, 64, 4*clmblkdst8, 8*clmblkdst8); \
			} \
			rsi += -2*64+clmblkdst8;	/* Restore source pointer */ \
		} \
		\
		EXEC248##clm(rsi += -4*clmblkdst8+2*64);	/* Next source pointer */ \
		EXEC248##clm(pfing(rcx += -4*blkdst+128));	/* Next prefetch pointer */ \
	} \
	EXEC248##clm(rsi += -clm*2*64+clmblkdst);/* Next source pointer */ \
	EXEC248##clm(pfing(rcx += -clm*128+4*blkdst));	/* Next prefetch pointer */ \
	EXEC1##clm(rsi += -4*clmblkdst8+clmblkdst);/* Next source pointer */ \
	\
} \
end_timer(18); \
\
/* Normalize these values */ \
\
pass1_normalize(1, yes_fft); \
\
/* Copy data from the scratch area back to the FFT data area */ \
\
copy_scratch_data_to_fft(128, clm); \
goto	endlp;			/* Skip FFT code if flag not set */ \
\
/* Fire up auxiliary threads to do pass 1 forward FFT */ \
\
/* pass1: */	pass1_forward_fft_setup(b0b); \
\
/* Copy data from the FFT data area back to the scratch area */ \
\
b0b:	copy_fft_data_to_scratch(128, clm); \
\
/* Do FFT levels 1,2,3 */ \
/* */ \
/* On input the 64-byte cache lines hold these data values: */ \
/*	0K	256K	512K	768K */ \
/*	2K	... */ \
/*	... */ \
/*	254K */ \
/* On output the 64-byte cache lines hold these data values: */ \
/*	0K	128K	256K	384K */ \
/*	2K	... */ \
/*	... */ \
/*	126K	... */ \
/*	512K	... */ \
/*	... */ \
/* */ \
/* Do 128 eight_reals_first_fft macros (FFT levels 1,2,3) */ \
/*	distance between fft data elements is 128K */ \
/*	do 64 iterations each processing 16 data values */ \
\
/* Do FFT levels 4,5 */ \
/* Values 0K-127K is real data, 128K-255K is semi-real data, 256K-1023K is */ \
/* complex data. */ \
/* */ \
/* On input the 64-byte cache lines hold these data values: */ \
/*	0K	128K	256K	384K */ \
/*	2K	... */ \
/*	... */ \
/*	126K	... */ \
/*	512K	... */ \
/*	... */ \
/* On output the 64-byte cache lines hold these data values: */ \
/*	0K	32K	64K	96K */ \
/*	2K	... */ \
/*	... */ \
/*	30K	... */ \
/*	128K	... */ \
/*	... */ \
/* */ \
/* Do 16 eight_reals_fft_2 macros */ \
/* Do 16 nop_two_two_complex_fft_2 macros */ \
/* Do 32 four_complex_fft macros */ \
/*	distance between fft data elements is 32K */ \
/*	do 16 macros each processing 32 data values */ \
/* */ \
/* Do 64 four_complex_fft macros */ \
/*	distance between fft data elements is 32K */ \
/*	do 16 macros each processing 32 data values */ \
\
/* Do FFT levels 6,7 */ \
/* Values 0K-31K is real data, 32K-63K is semi-real data, 64K-1023K is */ \
/* complex data. */ \
/* */ \
/* On input the 64-byte cache lines hold these data values: */ \
/*	0K	32K	64K	96K */ \
/*	2K	... */ \
/*	... */ \
/*	30K	... */ \
/*	128K	... */ \
/*	... */ \
/* On output the 64-byte cache lines hold these data values: */ \
/*	0K	8K	16K	24K */ \
/*	2K	... */ \
/*	... */ \
/*	6K	... */ \
/*	32K	... */ \
/*	... */ \
/* */ \
/* Do 8 eight_reals_fft_2 macros */ \
/* Do 8 two_two_complex_fft macros */ \
/*	distance between fft data elements is 8K */ \
/*	do 4 macros each processing 32 data values */ \
/* */ \
/* Do 112 four_complex_fft macros */ \
/*	distance between fft data elements is 8K */ \
/*	do 28 macros each processing 32 data values */ \
\
yes_fft:start_timer(21); \
rsi = (uintptr_t)g->scratch_area;	/* Get address of scratch area */ \
rcx = (uintptr_t)g->data_addr; /* Load source address for prefetching */ \
pfing(rbp = (uintptr_t)g->data_prefetch);/* Set prefetch distance */ \
pfing(rbp -= rcx); \
pfing(rcx += 64*blkdst);		/* Create prefetching pointer */ \
rdi = (uintptr_t)g->sincos3;/* Load sin/cos pointer */ \
rbx = (uintptr_t)g->sincos2;/* Load sin/cos pointer */ \
/* 4 iterations */ \
\
for(unsigned int loopA = 4; loopA; loopA--) { \
	xtouchlines(16);			/* Preload TLBs */ \
	\
	for(unsigned int loopB = clm; loopB; loopB--) { \
		for(unsigned int loopC = 8; loopC; loopC--) { \
			for(unsigned int loopD = 2; loopD; loopD--) { \
				for(unsigned int loopE = 2; loopE; loopE--) { \
					x2cl_eight_reals_fft(rsi, 64, 8*clmblkdst8); \
				} \
				rsi += -2*64+4*clmblkdst; \
			} \
			L2prefetch128(u8ptr(rcx+rbp)); \
			pfing(rcx += blkdst);		/* Next prefetch pointer */ \
			rsi += -2*4*clmblkdst+clmblkdst8; \
		} \
		\
		rsi += -8*clmblkdst8;	/* Restore source pointer */ \
		for(unsigned int loopC = 2; loopC; loopC--) { \
			for(unsigned int loopD = 2; loopD; loopD--) { \
				for(unsigned int loopE = 2; loopE; loopE--) { \
					x4cl_eight_reals_fft_2(rsi, 64, 2*clmblkdst8, 4*clmblkdst8, rdi); \
				} \
				rsi += -2*64+4*clmblkdst; \
			} \
			rsi += -2*4*clmblkdst+clmblkdst8; \
		} \
		\
		rsi += -2*clmblkdst8+8*clmblkdst8;/* Next source pointer */ \
		rdi += XMM_SCD;		/* Next sine/cosine pointer */ \
		for(unsigned int loopC = 2; loopC; loopC--) { \
			for(unsigned int loopD = 2; loopD; loopD--) { \
				for(unsigned int loopE = 2; loopE; loopE--) { \
					x4cl_four_complex_fft(rsi, 64, 2*clmblkdst8, 4*clmblkdst8, rdi); \
				} \
				rsi += -2*64+4*clmblkdst; \
			} \
			rsi += -2*4*clmblkdst+clmblkdst8; \
		} \
		rdi += -XMM_SCD;		/* Restore sine/cosine pointer */ \
		\
		rsi += -2*clmblkdst8-8*clmblkdst8;/* Restore source pointer */ \
		xloop1sp(b4b, 2*clmblkdst8, x4cl_eight_reals_fft_2(rsi, 64, 4*clmblkdst, clmblkdst8, rbx);); \
		\
		rbx += XMM_SCD;		/* Next sine/cosine pointer */ \
		/* 7 iterations */ \
		for(unsigned int loopC = 7; loopC; loopC--) { \
			xloop1sp(b5b, 2*clmblkdst8, x4cl_four_complex_fft(rsi, 64, 4*clmblkdst, clmblkdst8, rbx)); \
			rbx += 2*XMM_SCD;		/* Next sine/cosine pointer */ \
		} \
		rbx += -15*XMM_SCD;	/* Restore sine/cosine pointer */ \
		\
		EXEC248##clm(rsi += -8*2*clmblkdst8+2*64);/* Next source pointer */ \
		EXEC248##clm(pfing(rcx += -16*blkdst+128));	/* Next prefetch pointer */ \
	} \
	EXEC248##clm(rsi += -clm*2*64+clmblkdst);/* Next source pointer */ \
	EXEC248##clm(pfing(rcx += -clm*128+16*blkdst));	/* Next prefetch pointer */ \
	EXEC1##clm(rsi += -8*2*clmblkdst8+clmblkdst);/* Next source pointer */ \
	\
} \
rsi += -4*clmblkdst;	/* Restore source pointer */ \
pfing(rcx += -128*blkdst);	/* Restore prefetch pointer */ \
end_timer(21); \
\
/* Do FFT levels 8,9 */ \
/* Values 0K-7K is real data, 8K-15K is semi-real data, 16K-1023K is */ \
/* complex data. */ \
/* */ \
/* On input the 64-byte cache lines hold these data values: */ \
/*	0K	8K	16K	24K */ \
/*	2K	... */ \
/*	... */ \
/*	6K	... */ \
/*	32K	... */ \
/*	... */ \
/* On output the 64-byte cache lines hold these data values: */ \
/*	0K	2K	4K	6K */ \
/*	8K	... */ \
/*	... */ \
\
/* Do 1 eight_reals_fft_2 macros */ \
/* Do 1 nop_two_two_complex_fft_2 macros */ \
/* Do 2 four_complex_fft macros */ \
/*	distance between fft data elements is 2K */ \
/*	do 1 macros each processing 32 data values */ \
\
start_timer(24); \
rdi = (uintptr_t)g->sincos1;/* Load sin/cos pointer */ \
for(unsigned int loopA = 2*clm; loopA; loopA--) { \
	g4cl_eight_reals_fft_2(rsi, 64, clmblkdst, 2*clmblkdst, rcx, 64, blkdst, 2*blkdst); \
} \
\
/* Do 124 four_complex_fft macros */ \
/*	distance between fft data elements is 2K */ \
/*	do 31 macros each processing 32 data values */ \
\
rsi += -2*clm*64+4*clmblkdst;/* Next source pointer */ \
rcx += -2*clm*64+4*blkdst;	/* Next dest pointer */ \
rdi += XMM_SCD;		/* Next sine/cosine pointer */ \
/* eax = 16*256+1; */ /* 31 iterations of 2*clm */ \
for(unsigned int loopA = 31; loopA; ) { \
	do{ \
		for(unsigned int loopC = clm; loopC; loopC--) { \
			for(unsigned int loopD = 2; loopD; loopD--) { \
				g4cl_four_complex_fft(rsi, 64, clmblkdst, 2*clmblkdst, rcx, 64, blkdst, 2*blkdst); \
				L2prefetch128(u8ptr(rcx-64+rbp)); \
				L2prefetch128(u8ptr(rcx-64+blkdst+rbp)); \
				pfing(rbp += -64+2*blkdst); \
			} \
			pfing(rbp += 128-4*blkdst); \
		} \
		rsi += -2*clm*64+4*clmblkdst;/* Next source pointer */ \
		rcx += -2*clm*64+4*blkdst;	/* Next dest pointer */ \
		rdi += 2*XMM_SCD;		/* Next sine/cosine pointer */ \
	}while(--loopA & 1); \
	rsi += -8*clmblkdst+clmblkdst8; \
} \
end_timer(24); \
\
/* Work on next 4 sets of 512 values */ \
\
endlp:	xpass1_get_next_block(pass2, c0b, b0b); /* Get next block and jump */ \
/* done: */	xfft_footer_hg2; \
\
\
\
/* Perform a 1M-element negacyclic FFT.  This is done in two passes.  Pass 1 */ \
/* does 9 levels, pass 2 uses common code to perform the last 11 levels. */ \
\

#define hg_pass1sclevels9complex4(procname, xpass2name, pass2_levels, clm) \
 /* (fftname, pass2_levels, clmarg) */ \
\
/* set_FFT_constants(fftname, 0, 1, pass2_levels, clmarg); */ \
xfft_header_hg2(procname, xpass2name, pass2_levels, clm); \
\
/* Do the last 11 FFT levels */ \
/*	do 128 pass2_11_levels_complex */ \
\
/* pass2: */	start_timer(0); \
xpass2name(g); \
end_timer(0); \
if(g->ffttype == 1)		/* We're done if FFTing only */ \
	return; \
\
/* Do the inverse FFT */ \
\
pass1_inverse_fft_setup(c0b, procname); \
\
/* Do inverse FFT levels 8,9 */ \
/* On input the 64-byte cache lines hold these data values: */ \
/*	0K	2K	4K	6K */ \
/*	8K	10K	12K	14K */ \
/*	... */ \
/* On output the 64-byte cache lines hold these data values: */ \
/*	0K	4K	8K	12K */ \
/*	2K	... */ \
/*	16K	... */ \
/*	... */ \
\
/* Do 128 four_complex_unfft macros */ \
/*	distance between fft data elements is 2K */ \
/*	do 64 macros each processing 16 data values */ \
\
c0b:	start_timer(16); \
rsi = (uintptr_t)g->data_addr;/* Get source address */ \
rcx = (uintptr_t)g->scratch_area;	/* Get address of scratch area */ \
rdi = (uintptr_t)g->sincos1;/* Load sin/cos pointer */ \
/* 64 iterations */ \
pfing(rbp = (uintptr_t)g->data_prefetch);/* Set prefetch distance */ \
pfing(rbp -= rsi); \
for(unsigned int loopA = 16; loopA; loopA--) { \
	for(unsigned int loopB = 4; loopB; loopB--) { \
		for(unsigned int loopC = 2*clm; loopC; loopC--) { \
			g2cl_four_complex_unfft(rsi, 64, blkdst, rcx, 64, clmblkdst); \
		} \
		rsi += -2*clm*64+2*blkdst;	/* Next source pointer */ \
		rcx += -2*clm*64+2*clmblkdst;/* Next dest pointer */ \
		rdi += XMM_SCD;		/* Next sine/cosine pointer */ \
	} \
	rcx += -8*clmblkdst+clmblkdst8;/* Next dest pointer */ \
} \
rcx = rsi-64*2*blkdst;	/* Create prefetch pointer */ \
rsi = (uintptr_t)g->scratch_area;	/* Restore source pointer */ \
end_timer(16); \
\
/* Do inverse FFT levels 6,7 */ \
/* On input the 64-byte cache lines hold these data values: */ \
/*	0K	4K	8K	12K */ \
/*	2K	... */ \
/*	16K	... */ \
/*	... */ \
/* On output the 64-byte cache lines hold these data values: */ \
/*	0K	16K	32K	48K */ \
/*	2K	... */ \
/*	... */ \
/*	14K */ \
/*	64K	... */ \
/*	... */ \
\
/* Do 128 four_complex_unfft macros */ \
/*	distance between fft data elements is 8K */ \
/*	do 32 macros each processing 32 data values */ \
\
start_timer(17); \
rdi = (uintptr_t)g->sincos2;/* Load sin/cos pointer */ \
/* 16 iterations of 2 iters of 2*clm */ \
for(unsigned int loopA = 16; loopA; loopA--) { \
	xtouch(xptr(rcx+(clm-1)*128));	/* Preload TLB */ \
	xloop2_clm##clm(c5b, clmblkdst, x4cl_four_complex_unfft(rsi, 64, 2*clmblkdst, 4*clmblkdst, rdi));		/* 2 * 2*clm iterations */ \
	rsi += -2*clmblkdst+clmblkdst8;/* Next source pointer */ \
	rdi += XMM_SCD;		/* Next sine/cosine pointer */ \
} \
rsi += -16*clmblkdst8;	/* Restore source pointer */ \
end_timer(17); \
\
/* Do inverse FFT levels 4,5 */ \
/* On input the 64-byte cache lines hold these data values: */ \
/*	0K	16K	32K	48K */ \
/*	2K	... */ \
/*	... */ \
/*	14K */ \
/*	64K	... */ \
/*	... */ \
/* On output the 64-byte cache lines hold these data values: */ \
/*	0K	64K	128K	192K */ \
/*	2K	... */ \
/*	... */ \
/*	62K */ \
/*	256K	... */ \
/*	... */ \
\
/* Do 128 four_complex_unfft macros */ \
/*	distance between fft data elements is 32K */ \
/*	do 32 macros each processing 32 data values */ \
\
start_timer(18); \
rdi = (uintptr_t)g->sincos3;/* Load sin/cos pointer */ \
/* 4 iterations of 4 iters of 2*2*clm */ \
for(unsigned int loopA = 4; loopA; loopA--) { \
	for(unsigned int loopB = 4; loopB; loopB--) { \
		xtouch(xptr(rcx+(clm-1)*128));	/* Preload TLB */ \
		xloop2_clm##clm(c3b, clmblkdst, x4cl_four_complex_unfft(rsi, 64, clmblkdst8, 2*clmblkdst8, rdi));		/* 2* 2*clm iterations */ \
	} \
	rsi += -4*2*clmblkdst+4*clmblkdst8;/* Next source pointer */ \
	rdi += XMM_SCD;		/* Next sine/cosine pointer */ \
} \
rsi += -4*4*clmblkdst8;	/* Restore source pointer */ \
end_timer(18); \
\
/* Do inverse FFT levels 1,2,3 */ \
/* On input the 64-byte cache lines hold these data values: */ \
/*	0K	64K	128K	192K */ \
/*	2K	... */ \
/*	... */ \
/*	62K */ \
/*	256K	... */ \
/*	... */ \
/* On output the 64-byte cache lines hold these data values: */ \
/*	0K	256K	512K	768K */ \
/*	2K	... */ \
/*	... */ \
/*	254K */ \
\
/* Do 128 four_complex_last_unfft macros */ \
/*	distance between fft data elements is 128K */ \
/*	do 32 macros each processing 32 data values */ \
\
start_timer(19); \
rdi = (uintptr_t)g->plus1_premults;/* Address of premultiplier table */ \
/* 4 iters of 8 iters of 2*clm */ \
for(unsigned int loopA = 4; loopA; loopA--) { \
	for(unsigned int loopB = 8; loopB; loopB--) { \
		xtouch(xptr(rcx+(clm-1)*128));	/* Preload TLB */ \
		xloop1_clm##clm(c1b, clmblkdst, x4cl_four_complex_last_unfft(rsi, 64, 4*clmblkdst8, 8*clmblkdst8, 32*128));		/* 2*clm iterations */ \
		rdi += 128;		/* Next premultiplier pointer */ \
	} \
	rsi += -8*clmblkdst+clmblkdst8;/* Next source pointer */ \
} \
end_timer(19); \
\
/* Normalize these values */ \
\
pass1_normalize(1, yes_fft); \
\
/* Copy data from the scratch area back to the FFT data area */ \
\
copy_scratch_data_to_fft(128, clm); \
goto	endlp;			/* Skip FFT code if flag not set */ \
\
/* Fire up auxiliary threads to do pass 1 forward FFT */ \
\
/* pass1: */	pass1_forward_fft_setup(b0b); \
\
/* Copy data from the FFT data area back to the scratch area */ \
\
b0b:	copy_fft_data_to_scratch(128, clm); \
\
/* Do FFT levels 1,2,3 */ \
/* */ \
/* On input the 64-byte cache lines hold these data values: */ \
/*	0K	256K	512K	768K */ \
/*	2K	... */ \
/*	... */ \
/*	254K */ \
/* On output the 64-byte cache lines hold these data values: */ \
/*	0K	128K	256K	384K */ \
/*	2K	... */ \
/*	... */ \
/*	126K	... */ \
/*	512K	... */ \
/*	... */ \
\
/* Do 128 four_complex_first_fft macros (FFT levels 1,2,3) */ \
/*	distance between fft data elements is 128K */ \
/*	do 64 iterations each processing 16 data values */ \
\
yes_fft:start_timer(21); \
rsi = (uintptr_t)g->scratch_area;	/* Get address of scratch area */ \
rdi = (uintptr_t)g->plus1_premults;/* Address of premultiplier table */ \
rcx = (uintptr_t)g->data_addr; /* Load source address for prefetching */ \
pfing(rbp = (uintptr_t)g->data_prefetch);/* Set prefetch distance */ \
pfing(rbp -= rcx); \
pfing(rcx += 64*blkdst);		/* Create prefetching pointer */ \
/* 8 iters of 4 iters of 2*2*clm */ \
for(unsigned int loopA = 8; loopA; loopA--) { \
	for(unsigned int loopB = 4; loopB; loopB--) { \
		xtouch(xptr(rcx+(clm-1)*128));	/* Preload TLB */ \
		xloop2edi_clm##clm(b1b, clmblkdst, 128, x2cl_four_complex_first_fft_scratch(rsi, 64, 8*clmblkdst8)); /* 2 * 2*clm iterations */ \
	} \
	rsi += -8*clmblkdst+clmblkdst8;/* Next source pointer */ \
} \
rsi += -8*clmblkdst8;	/* Restore source pointer */ \
end_timer(21); \
\
/* Do FFT levels 4,5 */ \
/* */ \
/* On input the 64-byte cache lines hold these data values: */ \
/*	0K	128K	256K	384K */ \
/*	2K	... */ \
/*	... */ \
/*	126K	... */ \
/*	512K	... */ \
/*	... */ \
/* On output the 64-byte cache lines hold these data values: */ \
/*	0K	32K	64K	96K */ \
/*	2K	... */ \
/*	... */ \
/*	30K	... */ \
/*	128K	... */ \
/*	... */ \
\
/* Do 128 four_complex_fft macros */ \
/*	distance between fft data elements is 32K */ \
/*	do 32 macros each processing 32 data values */ \
\
start_timer(22); \
rdi = (uintptr_t)g->sincos3;/* Load sin/cos pointer */ \
/* 2 iters of 2 iters of 4 of 2*2*clm */ \
for(unsigned int loopA = 2; loopA; loopA--) { \
	for(unsigned int loopB = 2; loopB; loopB--) { \
		for(unsigned int loopC = 4; loopC; loopC--) { \
			xtouch(xptr(rcx+(clm-1)*128));	/* Preload TLB */ \
			xloop2_clm##clm(b3b, clmblkdst, x4cl_four_complex_fft(rsi, 64, 2*clmblkdst8, 4*clmblkdst8, rdi));		/* 2 * 2*clm iterations */ \
		} \
		rsi += -8*clmblkdst+clmblkdst8; \
	} \
	rsi += -2*clmblkdst8+8*clmblkdst8;/* Next source pointer */ \
	rdi += 2*XMM_SCD;		/* Next sine/cosine pointer */ \
} \
rsi += -2*8*clmblkdst8;	/* Restore source pointer */ \
end_timer(22); \
\
/* Do FFT levels 6,7 */ \
/* */ \
/* On input the 64-byte cache lines hold these data values: */ \
/*	0K	32K	64K	96K */ \
/*	2K	... */ \
/*	... */ \
/*	30K	... */ \
/*	128K	... */ \
/*	... */ \
/* On output the 64-byte cache lines hold these data values: */ \
/*	0K	8K	16K	24K */ \
/*	2K	... */ \
/*	... */ \
/*	6K	... */ \
/*	32K	... */ \
/*	... */ \
\
/* Do 128 four_complex_fft macros */ \
/*	distance between fft data elements is 8K */ \
/*	do 32 macros each processing 32 data values */ \
\
start_timer(23); \
rdi = (uintptr_t)g->sincos2;/* Load sin/cos pointer */ \
/* 8 iters of 2 iters of 2*2*clm */ \
for(unsigned int loopA = 8; loopA; loopA--) { \
	for(unsigned int loopB = 2; loopB; loopB--) { \
		xtouch(xptr(rcx+(clm-1)*128));	/* Preload TLB */ \
		xloop2_clm##clm(b5b, clmblkdst, x4cl_four_complex_fft(rsi, 64, 4*clmblkdst, clmblkdst8, rdi));		/* 2 * 2*clm iterations */ \
	} \
	rsi += -4*clmblkdst+2*clmblkdst8;/* Next source pointer */ \
	rdi += 2*XMM_SCD;		/* Next sine/cosine pointer */ \
} \
rsi += -8*2*clmblkdst8;	/* Restore source pointer */ \
pfing(rcx += -128*blkdst);	/* Restore prefetch pointer */ \
end_timer(23); \
\
/* Do FFT levels 8,9 */ \
/* */ \
/* On input the 64-byte cache lines hold these data values: */ \
/*	0K	8K	16K	24K */ \
/*	2K	... */ \
/*	... */ \
/*	6K	... */ \
/*	32K	... */ \
/*	... */ \
/* On output the 64-byte cache lines hold these data values: */ \
/*	0K	2K	4K	6K */ \
/*	8K	... */ \
/*	... */ \
\
/* Do 128 four_complex_fft macros */ \
/*	distance between fft data elements is 2K */ \
/*	do 32 macros each processing 32 data values */ \
\
start_timer(24); \
rdi = (uintptr_t)g->sincos1;/* Load sin/cos pointer */ \
/* 32 iterations of 2*clm */ \
for(unsigned int loopA = 16; loopA; loopA--) { \
	for(unsigned int loopB = 2; loopB; loopB--) { \
		for(unsigned int loopC = clm; loopC; loopC--) { \
			for(unsigned int loopD = 2; loopD; loopD--) { \
				g4cl_four_complex_fft(rsi, 64, clmblkdst, 2*clmblkdst, rcx, 64, blkdst, 2*blkdst); \
				L2prefetch128(u8ptr(rcx-64+rbp)); \
				L2prefetch128(u8ptr(rcx-64+blkdst+rbp)); \
				pfing(rbp += -64+2*blkdst); \
			} \
			pfing(rbp += 128-4*blkdst); \
		} \
	rsi += -2*clm*64+4*clmblkdst;/* Next source pointer */ \
		rcx += -2*clm*64+4*blkdst;	/* Next dest pointer */ \
		rdi += 2*XMM_SCD;		/* Next sine/cosine pointer */ \
	} \
	rsi += -8*clmblkdst+clmblkdst8; \
} \
end_timer(24); \
\
/* Work on next 4 sets of 512 values */ \
\
endlp:	xpass1_get_next_block(pass2, c0b, b0b); /* Get next block and jump */ \
/* done: */	xfft_footer_hg2; \
\
\
\
/* Perform a 1280K-element FFT.  This is done in two passes.  Pass 1 */ \
/* does 10 levels, pass 2 uses common code to perform the last 11 levels. */ \
\

#define hg_pass1sclevels10pfa5(procname, xpass2name, pass2_levels, clm) \
 /* (fftname, pass2_levels, clmarg) */ \
\
/* set_FFT_constants(fftname, 0, 0, pass2_levels, clmarg); */ \
xfft_header_hg2(procname, xpass2name, pass2_levels, clm); \
\
/* Do the last 11 FFT levels */ \
/*	do 1 pass2_11_levels_real and 159 pass2_11_levels_complex */ \
\
/* pass2: */	start_timer(0); \
xpass2name(g); \
end_timer(0); \
if(g->ffttype == 1)		/* We're done if FFTing only */ \
	return; \
\
/* Do the inverse FFT */ \
\
pass1_inverse_fft_setup(c0b, procname); \
\
/* Do inverse FFT levels 9,10 */ \
/* On input the 64-byte cache lines hold these data values: */ \
/*	0K	2K	4K	6K */ \
/*	8K	10K	12K	14K */ \
/*	... */ \
/* On output the 64-byte cache lines hold these data values: */ \
/*	0K	4K	8K	12K */ \
/*	2K	... */ \
/*	16K	... */ \
/*	... */ \
\
/* Do 1 eight_reals_unfft_2 macros */ \
/* Do 1 nop_two_two_complex_unfft_2 macros */ \
/*	distance between fft data elements is 2K */ \
/*	do 1 macros each processing 16 data values */ \
\
c0b:	start_timer(16); \
rsi = (uintptr_t)g->data_addr;/* Get source address */ \
rcx = (uintptr_t)g->scratch_area;	/* Get address of scratch area */ \
pfing(rbp = (uintptr_t)g->data_prefetch);/* Set prefetch distance */ \
pfing(rbp -= rsi); \
for(unsigned int loopA = 2*clm; loopA; loopA--) { \
	g2cl_eight_reals_unfft_2(rsi, 64, blkdst, rcx, 64, clmblkdst); \
} \
\
/* Do 158 four_complex_unfft macros */ \
/*	distance between fft data elements is 2K */ \
/*	do 79 macros each processing 16 data values */ \
\
rsi += -2*clm*64+2*blkdst;	/* Next source pointer */ \
rcx += -2*clm*64+2*clmblkdst;/* Next dest pointer */ \
rdi = (uintptr_t)g->sincos1;/* Load sin/cos pointer */ \
/* al = 256/4+20; */ /* 79 iterations of 2*clm */ \
for(unsigned int loopA = 79; loopA; ) { \
	do{ \
		for(unsigned int loopC = 2*clm; loopC; loopC--) { \
			g2cl_four_complex_unfft(rsi, 64, blkdst, rcx, 64, clmblkdst); \
		} \
		rsi += -2*clm*64+2*blkdst;	/* Next source pointer */ \
		rcx += -2*clm*64+2*clmblkdst;/* Next dest pointer */ \
		rdi += XMM_SCD;		/* Next sine/cosine pointer */ \
	}while(--loopA & 3); \
	rcx += -8*clmblkdst+clmblkdst8; \
} \
rcx = rsi-80*2*blkdst;	/* Create prefetch pointer */ \
rsi = (uintptr_t)g->scratch_area;	/* Restore source pointer */ \
end_timer(16); \
\
/* Do inverse FFT levels 7,8 */ \
/* On input the 64-byte cache lines hold these data values: */ \
/*	0K	4K	8K	12K */ \
/*	2K	... */ \
/*	16K	... */ \
/*	... */ \
/* On output the 64-byte cache lines hold these data values: */ \
/*	0K	16K	32K	48K */ \
/*	2K	... */ \
/*	... */ \
/*	14K */ \
/*	64K	... */ \
/*	... */ \
\
/* Do 4 eight_reals_unfft_2 macros */ \
/* Do 4 nop_two_two_complex_unfft_2 macros */ \
/*	distance between fft data elements is 8K */ \
/*	do 2 macros each processing 32 data values */ \
\
start_timer(17); \
/* 2 iterations of 2*clm */ \
for(unsigned int loopA = 2; loopA; loopA--) { \
	for(unsigned int loopB = 2*clm; loopB; loopB--) { \
		x4cl_eight_reals_unfft_2(rsi, 64, 2*clmblkdst, 4*clmblkdst); \
	} \
	rsi += -2*clm*64+clmblkdst;/* Next source pointer */ \
} \
\
/* Do 152 four_complex_unfft macros */ \
/*	distance between fft data elements is 8K */ \
/*	do 38 macros each processing 32 data values */ \
\
rsi += -2*clmblkdst+clmblkdst8;/* Next source pointer */ \
rdi = (uintptr_t)g->sincos2;/* Load sin/cos pointer */ \
/* 19 iterations of 2 iters of 2*clm */ \
for(unsigned int loopA = 19; loopA; loopA--) { \
	for(unsigned int loopB = 2; loopB; loopB--) { \
		xtouch(xptr(rcx+(clm-1)*128));	/* Preload TLB */ \
		xloop1_clm##clm(c7b, clmblkdst, x4cl_four_complex_unfft(rsi, 64, 2*clmblkdst, 4*clmblkdst, rdi)); \
	} \
	rsi += -2*clmblkdst+clmblkdst8;/* Next source pointer */ \
	rdi += XMM_SCD;		/* Next sine/cosine pointer */ \
} \
rsi += -20*clmblkdst8;	/* Restore source pointer */ \
end_timer(17); \
\
/* Do inverse FFT levels 6 */ \
/* On input the 64-byte cache lines hold these data values: */ \
/*	0K	16K	32K	48K */ \
/*	2K	... */ \
/*	... */ \
/*	14K */ \
/*	64K	... */ \
/*	... */ \
/* On output the 64-byte cache lines hold these data values: */ \
/*	0K	32K	64K	96K */ \
/*	2K	... */ \
/*	... */ \
/*	30K */ \
/*	128K	... */ \
/*	... */ \
/* */ \
/* Do 16 eight_reals_unfft_1 macros */ \
/*	distance between fft data elements is 32K */ \
/*	do 8 macros each processing 16 data values */ \
/* */ \
/* Do 152 four_complex_unfft macros */ \
/*	distance between fft data elements is 32K */ \
/*	do 76 macros each processing 16 data values */ \
\
/* Do inverse FFT levels 4,5 */ \
/* On input the 64-byte cache lines hold these data values: */ \
/*	0K	32K	64K	96K */ \
/*	2K	... */ \
/*	... */ \
/*	30K */ \
/*	128K	... */ \
/*	... */ \
/* On output the 64-byte cache lines hold these data values: */ \
/*	0K	64K	128K	192K */ \
/*	2K	... */ \
/*	... */ \
/*	62K */ \
/*	256K	384K	512K	768K */ \
/*	... */ \
/* */ \
/* Do 32 eight_reals_unfft_2 macros */ \
/*	distance between fft data elements is 64K */ \
/*	do 16 macros each processing 16 data values */ \
/* */ \
/* Do 128 four_complex_unfft macros */ \
/*	distance between fft data elements is 64K */ \
/*	do 32 macros each processing 32 data values */ \
\
/* Do inverse FFT levels 1,2,3 */ \
/* On input the 64-byte cache lines hold these data values: */ \
/*	0K	64K	128K	192K */ \
/*	2K	... */ \
/*	... */ \
/*	62K */ \
/*	256K	384K	512K	768K */ \
/*	... */ \
/* On output the 64-byte cache lines hold these data values: */ \
/*	0K	64K	640K	704K */ \
/*	2K	... */ \
/*	... */ \
/*	62K */ \
/*	128K	384K	768K	1024K */ \
/*	130K */ \
/*	... */ \
/*	382K */ \
/* */ \
/* Do 256 five_reals_last_unfft macros */ \
/*	distance between fft data elements is 256K */ \
/*	do 32 macros each processing 40 data values */ \
\
start_timer(18); \
rdi = (uintptr_t)g->sincos3;/* Load sin/cos pointer */ \
rbx = (uintptr_t)g->sincos4;/* Load sin/cos pointer */ \
/* 8 iterations  */ \
\
for(unsigned int loopA = 8; loopA; loopA--) { \
	xtouchlines(8);			/* Preload TLBs */ \
	\
	for(unsigned int loopB = clm; loopB; loopB--) { \
		for(unsigned int loopC = 2; loopC; loopC--) { \
			x2cl_eight_reals_unfft_1(rsi, 64, clmblkdst8); \
		} \
		\
		rsi += -2*64+2*clmblkdst8;	/* Next source pointer */ \
		/* 9 iterations */ \
		for(unsigned int loopC = 9; loopC; loopC--) { \
			for(unsigned int loopD = 2; loopD; loopD--) { \
				x2cl_two_complex_unfft(rsi, 64, clmblkdst8); \
			} \
			rsi += -2*64+2*clmblkdst8;	/* Next source pointer */ \
			rdi += XMM_SCD;		/* Next sine/cosine pointer */ \
		} \
		rdi += -9*XMM_SCD;		/* Restore sine/cosine pointer */ \
		\
		rsi += -10*2*clmblkdst8;	/* Restore source pointer */ \
		for(unsigned int loopC = 2; loopC; loopC--) { \
			for(unsigned int loopD = 2; loopD; loopD--) { \
				x2cl_half_eight_reals_unfft_2(rsi, 64, 2*clmblkdst8); \
			} \
			rsi += -2*64+clmblkdst8; \
		} \
		\
		rsi += -2*clmblkdst8+4*clmblkdst8;/* Next source pointer */ \
		for(unsigned int loopC = 2; loopC; loopC--) { \
			for(unsigned int loopD = 2; loopD; loopD--) { \
				xloop1sp(c3b, clmblkdst8, x4cl_four_complex_unfft(rsi, 64, 2*clmblkdst8, 4*clmblkdst8, rbx)); \
			} \
			rsi += -2*clmblkdst8+8*clmblkdst8; \
			rbx += XMM_SCD;		/* Next sine/cosine pointer */ \
		} \
		rbx += -2*XMM_SCD;		/* Restore sine/cosine pointer */ \
		\
		rsi += -2*8*clmblkdst8-4*clmblkdst8;/* Restore source pointer */ \
		xfive_reals_unfft_preload;	/* Preload consts for five real macros */ \
		for(unsigned int loopC = 4; loopC; loopC--) { \
			xloop1sp(c1b, clmblkdst8, x5cl_five_reals_last_unfft(rsi, 64, 4*clmblkdst8)); \
		} \
		\
		EXEC248##clm(rsi += -4*clmblkdst8+2*64);	/* Next source pointer */ \
		EXEC248##clm(pfing(rcx += -8*blkdst+128));	/* Next prefetch pointer */ \
	} \
	EXEC248##clm(rsi += -clm*2*64+clmblkdst);/* Next source pointer */ \
	EXEC248##clm(pfing(rcx += -clm*128+8*blkdst));	/* Next prefetch pointer */ \
	EXEC1##clm(rsi += -4*clmblkdst8+clmblkdst);/* Next source pointer */ \
	\
} \
end_timer(18); \
\
/* Normalize these values */ \
\
pass1_normalize(1, yes_fft); \
\
/* Copy data from the scratch area back to the FFT data area */ \
\
copy_scratch_data_to_fft(160, clm); \
goto	endlp;			/* Skip FFT code if flag not set */ \
\
/* Fire up auxiliary threads to do pass 1 forward FFT */ \
\
/* pass1: */	pass1_forward_fft_setup(b0b); \
\
/* Copy data from the FFT data area back to the scratch area */ \
\
b0b:	copy_fft_data_to_scratch(160, clm); \
\
/* Do FFT levels 1,2,3 */ \
/* */ \
/* On input the 64-byte cache lines hold these data values: */ \
/*	0K	64K	640K	704K */ \
/*	2K	... */ \
/*	... */ \
/*	62K */ \
/*	128K	384K	768K	1024K */ \
/*	130K */ \
/*	... */ \
/*	382K */ \
/* On output the 64-byte cache lines hold these data values: */ \
/*	0K	64K	128K	192K */ \
/*	2K	... */ \
/*	... */ \
/*	62K */ \
/*	256K	512K	768K	1024K */ \
/*	258K */ \
/*	... */ \
/*	510K */ \
/* */ \
/* Do 256 five_reals_first_fft macros (FFT levels 1,2,3) */ \
/*	distance between fft data elements is 256K */ \
/*	do 32 iterations each processing 40 data values */ \
\
/* Do FFT levels 4,5 */ \
/* Values 0K-255K is real data, 256K-1279K is complex data. */ \
/* */ \
/* On input the 64-byte cache lines hold these data values: */ \
/*	0K	64K	128K	192K */ \
/*	2K	... */ \
/*	... */ \
/*	62K */ \
/*	256K	512K	768K	1024K */ \
/*	258K */ \
/*	... */ \
/*	510K */ \
/* On output the 64-byte cache lines hold these data values: */ \
/*	0K	64K	128K	192K */ \
/*	2K	... */ \
/*	... */ \
/*	62K	... */ \
/*	256K	... */ \
/*	... */ \
/* */ \
/* Do 32 eight_reals_fft_2 macros */ \
/*	distance between fft data elements is 64K */ \
/*	do 32 macros each processing 8 data values */ \
/* */ \
/* Do 128 four_complex_fft macros */ \
/*	distance between fft data elements is 64K */ \
/*	do 32 macros each processing 32 data values */ \
\
/* Do FFT levels 6 */ \
/* Values 0K-63K is real data, 64K-127K is semi-real data, 128K-1279K is */ \
/* complex data. */ \
/* */ \
/* On input the 64-byte cache lines hold these data values: */ \
/*	0K	64K	128K	192K */ \
/*	2K	... */ \
/*	... */ \
/*	62K	... */ \
/*	256K	... */ \
/*	... */ \
/* On output the 64-byte cache lines hold these data values: */ \
/*	0K	32K	64K	96K */ \
/*	2K	... */ \
/*	... */ \
/*	30K	... */ \
/*	128K	... */ \
/*	... */ \
/* */ \
/* Do 16 eight_reals_fft_1 macros */ \
/* Do 16 two_two_complex_fft macros */ \
/*	distance between fft data elements is 32K */ \
/*	do 16 macros each processing 16 data values */ \
/* */ \
/* Do 128 two_two_complex_fft macros */ \
/*	distance between fft data elements is 32K */ \
/*	do 64 macros each processing 16 data values */ \
\
yes_fft:start_timer(22); \
rsi = (uintptr_t)g->scratch_area;	/* Get address of scratch area */ \
rcx = (uintptr_t)g->data_addr; /* Load source address for prefetching */ \
pfing(rbp = (uintptr_t)g->data_prefetch);/* Set prefetch distance */ \
pfing(rbp -= rcx); \
pfing(rcx += 60*blkdst);		/* Create prefetching pointer */ \
rdi = (uintptr_t)g->sincos4;/* Load sin/cos pointer */ \
rbx = (uintptr_t)g->sincos3;/* Load sin/cos pointer */ \
/* 16 iterations */ \
\
for(unsigned int loopA = 16; loopA; loopA--) { \
	for(unsigned int loopB = 8; loopB; loopB--) { \
		for(unsigned int loopC = clm; loopC; loopC--) { \
			xtouchlines(4);			/* Preload TLBs */ \
			\
			xfive_reals_fft_preload;		/* Preload consts for five real macros */ \
			for(unsigned int loopD = 2; loopD; loopD--) { \
				xloop1sp(b1b, 2*clmblkdst8, x5cl_five_reals_fft(rsi, 64, 4*clmblkdst8)); \
			} \
			\
			rsi += -2*2*clmblkdst8;	/* Restore source pointer */ \
			for(unsigned int loopD = 2; loopD; loopD--) { \
				for(unsigned int loopE = 2; loopE; loopE--) { \
					x1cl_half_eight_reals_fft_2(rsi, 64); \
				} \
				rsi += -2*64+2*clmblkdst8;	/* Next source pointer */ \
			} \
			\
			/*	bump	rsi, -2*2*clmblkdst8+4*clmblkdst8;; Next source pointer */ \
			for(unsigned int loopD = 2; loopD; loopD--) { \
				xloop1sp(b3b, 2*clmblkdst8, x4cl_four_complex_fft(rsi, 64, 4*clmblkdst8, 8*clmblkdst8, rdi)); \
			} \
			\
			rsi += -4*clmblkdst8-4*clmblkdst8;/* Restore source pointer */ \
			for(unsigned int loopD = 2; loopD; loopD--) { \
				x2cl_eight_reals_fft_1(rsi, 64, 2*clmblkdst8); \
			} \
			\
			rsi += -2*64+4*clmblkdst8;	/* Next source pointer */ \
			rbx += XMM_SCD;		/* Next sine/cosine pointer */ \
			for(unsigned int loopD = 4; loopD; loopD--) { \
				for(unsigned int loopE = 2; loopE; loopE--) { \
					x2cl_two_complex_fft(rsi, 64, 2*clmblkdst8, rbx); \
				} \
				rsi += -2*64+4*clmblkdst8;	/* Next source pointer */ \
				rbx += 2*XMM_SCD;		/* Next sine/cosine pointer */ \
			} \
			rbx += -9*XMM_SCD;		/* Restore sine/cosine pointer */ \
			\
				EXEC248##clm(rsi += -5*4*clmblkdst8+2*64);/* Next source pointer */ \
			EXEC248##clm(pfing(rcx += -4*blkdst+128));	/* Next prefetch pointer */ \
		} \
		EXEC248##clm(rsi += -clm*2*64+clmblkdst);/* Next source pointer */ \
		EXEC248##clm(pfing(rcx += -clm*128+4*blkdst));	/* Next prefetch pointer */ \
		EXEC1##clm(rsi += -5*4*clmblkdst8+clmblkdst);/* Next source pointer */ \
		\
	} \
	rsi += -8*clmblkdst+clmblkdst8;/* Next source pointer */ \
} \
rsi += -2*clmblkdst8;	/* Restore source pointer */ \
end_timer(22); \
\
/* Do FFT levels 7,8 */ \
/* Values 0K-31K is real data, 32K-63K is semi-real data, 64K-1279K is */ \
/* complex data. */ \
/* */ \
/* On input the 64-byte cache lines hold these data values: */ \
/*	0K	32K	64K	96K */ \
/*	2K	... */ \
/*	... */ \
/*	30K	... */ \
/*	128K	... */ \
/*	... */ \
/* On output the 64-byte cache lines hold these data values: */ \
/*	0K	8K	16K	24K */ \
/*	2K	... */ \
/*	... */ \
/*	6K	... */ \
/*	32K	... */ \
/*	... */ \
\
/* Do 8 eight_reals_fft_2 macros */ \
/* Do 8 two_two_complex_fft macros */ \
/*	distance between fft data elements is 8K */ \
/*	do 4 macros each processing 32 data values */ \
\
start_timer(25); \
rdi = (uintptr_t)g->sincos2;/* Load sin/cos pointer */ \
/* 4 iterations of 2*clm */ \
for(unsigned int loopA = 4; loopA; loopA--) { \
	for(unsigned int loopB = 2*clm; loopB; loopB--) { \
		x4cl_eight_reals_fft_2(rsi, 64, 4*clmblkdst, clmblkdst8, rdi); \
	} \
	rsi += -2*clm*64+clmblkdst;/* Next source pointer */ \
} \
\
/* Do 144 four_complex_fft macros */ \
/*	distance between fft data elements is 8K */ \
/*	do 36 macros each processing 32 data values */ \
\
rsi += -4*clmblkdst+2*clmblkdst8;/* Load source pointer */ \
rdi += XMM_SCD;		/* Next sine/cosine pointer */ \
/* 9 iterations of 4 iters of 2*clm */ \
for(unsigned int loopA = 9; loopA; loopA--) { \
	for(unsigned int loopB = 4; loopB; loopB--) { \
		xtouch(xptr(rcx+(clm-1)*128));	/* Preload TLB */ \
		xloop1_clm##clm(b7b, clmblkdst, x4cl_four_complex_fft(rsi, 64, 4*clmblkdst, clmblkdst8, rdi)); \
	} \
	rsi += -4*clmblkdst+2*clmblkdst8;/* Next source pointer */ \
	rdi += 2*XMM_SCD;		/* Next sine/cosine pointer */ \
} \
rsi += -10*2*clmblkdst8;	/* Restore source pointer */ \
pfing(rcx += -160*blkdst);	/* Restore prefetch pointer */ \
end_timer(25); \
\
/* Do FFT levels 9,10 */ \
/* Values 0K-7K is real data, 8K-15K is semi-real data, 16K-1279K is */ \
/* complex data. */ \
/* */ \
/* On input the 64-byte cache lines hold these data values: */ \
/*	0K	8K	16K	24K */ \
/*	2K	... */ \
/*	... */ \
/*	6K	... */ \
/*	32K	... */ \
/*	... */ \
/* On output the 64-byte cache lines hold these data values: */ \
/*	0K	2K	4K	6K */ \
/*	8K	... */ \
/*	... */ \
\
/* Do 1 eight_reals_fft_2 macros */ \
/* Do 1 nop_two_two_complex_fft_2 macros */ \
/* Do 2 four_complex_fft macros */ \
/*	distance between fft data elements is 2K */ \
/*	do 1 macros each processing 32 data values */ \
\
start_timer(26); \
rdi = (uintptr_t)g->sincos1;/* Load sin/cos pointer */ \
for(unsigned int loopA = 2*clm; loopA; loopA--) { \
	g4cl_eight_reals_fft_2(rsi, 64, clmblkdst, 2*clmblkdst, rcx, 64, blkdst, 2*blkdst); \
} \
\
/* Do 156 four_complex_fft macros */ \
/*	distance between fft data elements is 2K */ \
/*	do 39 macros each processing 32 data values */ \
\
rsi += -2*clm*64+4*clmblkdst;/* Next source pointer */ \
rcx += -2*clm*64+4*blkdst;	/* Next dest pointer */ \
rdi += XMM_SCD;		/* Next sine/cosine pointer */ \
/* eax = 20*256+1; */ /* 39 iterations of 2*clm */ \
for(unsigned int loopA = 39; loopA; ) { \
	do{ \
		for(unsigned int loopC = clm; loopC; loopC--) { \
			for(unsigned int loopD = 2; loopD; loopD--) { \
				g4cl_four_complex_fft(rsi, 64, clmblkdst, 2*clmblkdst, rcx, 64, blkdst, 2*blkdst); \
				L2prefetch128(u8ptr(rcx-64+rbp)); \
				L2prefetch128(u8ptr(rcx-64+blkdst+rbp)); \
				pfing(rbp += -64+2*blkdst); \
			} \
			pfing(rbp += 128-4*blkdst); \
		} \
		rsi += -2*clm*64+4*clmblkdst;/* Next source pointer */ \
		rcx += -2*clm*64+4*blkdst;	/* Next dest pointer */ \
		rdi += 2*XMM_SCD;		/* Next sine/cosine pointer */ \
	}while(--loopA & 1); \
	rsi += -8*clmblkdst+clmblkdst8; \
} \
end_timer(26); \
\
/* Work on next 4 sets of 640 values */ \
\
endlp:	xpass1_get_next_block(pass2, c0b, b0b); /* Get next block and jump */ \
/* done: */	xfft_footer_hg2; \
\
\
\
/* Perform a 1536K-element FFT.  This is done in two passes.  Pass 1 */ \
/* does 10 levels, pass 2 uses common code to perform the last 11 levels. */ \
\

#define hg_pass1sclevels10pfa6(procname, xpass2name, pass2_levels, clm) \
 /* (fftname, pass2_levels, clmarg) */ \
\
/* set_FFT_constants(fftname, 0, 0, pass2_levels, clmarg); */ \
xfft_header_hg2(procname, xpass2name, pass2_levels, clm); \
\
/* Do the last 11 FFT levels */ \
/*	do 1 pass2_11_levels_real and 191 pass2_11_levels_complex */ \
\
/* pass2: */	start_timer(0); \
xpass2name(g); \
end_timer(0); \
if(g->ffttype == 1)		/* We're done if FFTing only */ \
	return; \
\
/* Do the inverse FFT */ \
\
pass1_inverse_fft_setup(c0b, procname); \
\
/* Do inverse FFT levels 9,10 */ \
/* On input the 64-byte cache lines hold these data values: */ \
/*	0K	2K	4K	6K */ \
/*	8K	10K	12K	14K */ \
/*	... */ \
/* On output the 64-byte cache lines hold these data values: */ \
/*	0K	4K	8K	12K */ \
/*	2K	... */ \
/*	16K	... */ \
/*	... */ \
\
/* Do 1 eight_reals_unfft_2 macros */ \
/* Do 1 nop_two_two_complex_unfft_2 macros */ \
/*	distance between fft data elements is 2K */ \
/*	do 1 macros each processing 16 data values */ \
\
c0b:	start_timer(16); \
rsi = (uintptr_t)g->data_addr;/* Get source address */ \
rcx = (uintptr_t)g->scratch_area;	/* Get address of scratch area */ \
pfing(rbp = (uintptr_t)g->data_prefetch);/* Set prefetch distance */ \
pfing(rbp -= rsi); \
for(unsigned int loopA = 2*clm; loopA; loopA--) { \
	g2cl_eight_reals_unfft_2(rsi, 64, blkdst, rcx, 64, clmblkdst); \
} \
\
/* Do 190 four_complex_unfft macros */ \
/*	distance between fft data elements is 2K */ \
/*	do 95 macros each processing 16 data values */ \
\
rsi += -2*clm*64+2*blkdst;	/* Next source pointer */ \
rcx += -2*clm*64+2*clmblkdst;/* Next dest pointer */ \
rdi = (uintptr_t)g->sincos1;/* Load sin/cos pointer */ \
/* al = 256/4+24; */ /* 95 iterations of 2*clm */ \
for(unsigned int loopA = 95; loopA; ) { \
	do{ \
		for(unsigned int loopC = 2*clm; loopC; loopC--) { \
			g2cl_four_complex_unfft(rsi, 64, blkdst, rcx, 64, clmblkdst); \
		} \
		rsi += -2*clm*64+2*blkdst;	/* Next source pointer */ \
		rcx += -2*clm*64+2*clmblkdst;/* Next dest pointer */ \
		rdi += XMM_SCD;		/* Next sine/cosine pointer */ \
	}while(--loopA & 3); \
	rcx += -8*clmblkdst+clmblkdst8; \
} \
rcx = rsi-96*2*blkdst;	/* Create prefetch pointer */ \
rsi = (uintptr_t)g->scratch_area;	/* Restore source pointer */ \
end_timer(16); \
\
/* Do inverse FFT levels 7,8 */ \
/* On input the 64-byte cache lines hold these data values: */ \
/*	0K	4K	8K	12K */ \
/*	2K	... */ \
/*	16K	... */ \
/*	... */ \
/* On output the 64-byte cache lines hold these data values: */ \
/*	0K	16K	32K	48K */ \
/*	2K	... */ \
/*	... */ \
/*	14K */ \
/*	64K	... */ \
/*	... */ \
\
/* Do 4 eight_reals_unfft_2 macros */ \
/* Do 4 nop_two_two_complex_unfft_2 macros */ \
/*	distance between fft data elements is 8K */ \
/*	do 2 macros each processing 32 data values */ \
\
start_timer(17); \
/* 2 iterations of 2*clm */ \
for(unsigned int loopA = 2; loopA; loopA--) { \
	for(unsigned int loopB = 2*clm; loopB; loopB--) { \
		x4cl_eight_reals_unfft_2(rsi, 64, 2*clmblkdst, 4*clmblkdst); \
	} \
	rsi += -2*clm*64+clmblkdst;/* Next source pointer */ \
} \
\
/* Do 184 four_complex_unfft macros */ \
/*	distance between fft data elements is 8K */ \
/*	do 46 macros each processing 32 data values */ \
\
rsi += -2*clmblkdst+clmblkdst8;/* Next source pointer */ \
rdi = (uintptr_t)g->sincos2;/* Load sin/cos pointer */ \
/* 23 iterations of 2 iters of 2*clm */ \
for(unsigned int loopA = 23; loopA; loopA--) { \
	for(unsigned int loopB = 2; loopB; loopB--) { \
		xtouch(xptr(rcx+(clm-1)*128));	/* Preload TLB */ \
		c7b:	x4cl_four_complex_unfft(rsi, 64, 2*clmblkdst, 4*clmblkdst, rdi); \
		xloop1_clm##clm(c7b, clmblkdst); \
	} \
	rsi += -2*clmblkdst+clmblkdst8;/* Next source pointer */ \
	rdi += XMM_SCD;		/* Next sine/cosine pointer */ \
} \
rsi += -24*clmblkdst8;	/* Restore source pointer */ \
end_timer(17); \
\
/* Do inverse FFT levels 6 */ \
/* On input the 64-byte cache lines hold these data values: */ \
/*	0K	16K	32K	48K */ \
/*	2K	... */ \
/*	... */ \
/*	14K */ \
/*	64K	... */ \
/*	... */ \
/* On output the 64-byte cache lines hold these data values: */ \
/*	0K	32K	64K	96K */ \
/*	2K	... */ \
/*	... */ \
/*	30K */ \
/*	128K	... */ \
/*	... */ \
/* */ \
/* Do 16 eight_reals_unfft_1 macros */ \
/*	distance between fft data elements is 32K */ \
/*	do 8 macros each processing 16 data values */ \
/* */ \
/* Do 176 four_complex_unfft macros */ \
/*	distance between fft data elements is 32K */ \
/*	do 88 macros each processing 16 data values */ \
\
/* Do inverse FFT levels 4,5 */ \
/* On input the 64-byte cache lines hold these data values: */ \
/*	0K	32K	64K	96K */ \
/*	2K	... */ \
/*	... */ \
/*	30K */ \
/*	128K	... */ \
/*	... */ \
/* On output the 64-byte cache lines hold these data values: */ \
/*	0K	128K	256K	384K */ \
/*	2K	... */ \
/*	... */ \
/*	126K */ \
/*	512K	... */ \
/*	... */ \
/* */ \
/* Do 32 eight_reals_unfft_2 macros */ \
/* Do 32 nop_two_two_complex_unfft_2 macros */ \
/*	distance between fft data elements is 64K */ \
/*	do 16 macros each processing 32 data values */ \
/* */ \
/* Do 128 four_complex_unfft macros */ \
/*	distance between fft data elements is 64K */ \
/*	do 32 macros each processing 32 data values */ \
\
/* Do inverse FFT levels 1,2,3 */ \
/* On input the 64-byte cache lines hold these data values: */ \
/*	0K	128K	256K	384K */ \
/*	2K	... */ \
/*	... */ \
/*	126K */ \
/*	512K	768K	1024K	1280K */ \
/*	514K */ \
/*	... */ \
/*	766K */ \
/* On output the 64-byte cache lines hold these data values: */ \
/*	0K	384K	768K	1152K */ \
/*	... */ \
/*	382K */ \
/* */ \
/* Do 256 six_reals_last_unfft macros */ \
/*	distance between fft data elements is 256K */ \
/*	do 64 macros each processing 24 data values */ \
\
start_timer(18); \
rdi = (uintptr_t)g->sincos3;/* Load sin/cos pointer */ \
rbx = (uintptr_t)g->sincos4;/* Load sin/cos pointer */ \
/* 8 iterations  */ \
xsix_reals_unfft_preload; \
\
for(unsigned int loopA = 8; loopA; loopA--) { \
	xtouchlines(12);			/* Preload TLB */ \
	\
	for(unsigned int loopB = clm; loopB; loopB--) { \
		for(unsigned int loopC = 2; loopC; loopC--) { \
			x2cl_eight_reals_unfft_1(rsi, 64, clmblkdst8); \
		} \
		\
		rsi += -2*64+2*clmblkdst8;	/* Next source pointer */ \
		/* 11 iterations */ \
		for(unsigned int loopC = 11; loopC; loopC--) { \
			for(unsigned int loopD = 2; loopD; loopD--) { \
				x2cl_two_complex_unfft(rsi, 64, clmblkdst8); \
			} \
			rsi += -2*64+2*clmblkdst8;	/* Next source pointer */ \
			rdi += XMM_SCD;		/* Next sine/cosine pointer */ \
		} \
		rdi += -11*XMM_SCD;	/* Restore sine/cosine pointer */ \
		\
		rsi += -12*2*clmblkdst8;	/* Restore source pointer */ \
		for(unsigned int loopC = 2; loopC; loopC--) { \
			for(unsigned int loopD = 2; loopD; loopD--) { \
				x4cl_eight_reals_unfft_2(rsi, 64, 2*clmblkdst8, 4*clmblkdst8); \
			} \
			rsi += -2*64+clmblkdst8; \
		} \
		\
		rsi += -2*clmblkdst8+8*clmblkdst8;/* Next source pointer */ \
		for(unsigned int loopC = 2; loopC; loopC--) { \
			for(unsigned int loopD = 2; loopD; loopD--) { \
				xloop1sp(c3b, clmblkdst8, x4cl_four_complex_unfft(rsi, 64, 2*clmblkdst8, 4*clmblkdst8, rbx)); \
			} \
			rsi += -2*clmblkdst8+8*clmblkdst8; \
			rbx += XMM_SCD;		/* Next sine/cosine pointer */ \
		} \
		rbx += -2*XMM_SCD;		/* Restore sine/cosine pointer */ \
		\
		rsi += -3*8*clmblkdst8;	/* Restore source pointer */ \
		for(unsigned int loopC = 8; loopC; loopC--) { \
			xloop1sp(c1b, clmblkdst8, x3cl_six_reals_last_unfft(rsi, 64, 8*clmblkdst8)); \
		} \
		\
		EXEC248##clm(rsi += -8*clmblkdst8+2*64);	/* Next source pointer */ \
		EXEC248##clm(pfing(rcx += -12*blkdst+128));	/* Next prefetch pointer */ \
	} \
	EXEC248##clm(rsi += -clm*2*64+clmblkdst);/* Next source pointer */ \
	EXEC248##clm(pfing(rcx += -clm*128+12*blkdst));	/* Next prefetch pointer */ \
	EXEC1##clm(rsi += -8*clmblkdst8+clmblkdst);/* Next source pointer */ \
	\
} \
end_timer(18); \
\
/* Normalize these values */ \
\
pass1_normalize(1, yes_fft); \
\
/* Copy data from the scratch area back to the FFT data area */ \
\
copy_scratch_data_to_fft(192, clm); \
goto	endlp;			/* Skip FFT code if flag not set */ \
\
/* Fire up auxiliary threads to do pass 1 forward FFT */ \
\
/* pass1: */	pass1_forward_fft_setup(b0b); \
\
/* Copy data from the FFT data area back to the scratch area */ \
\
b0b:	copy_fft_data_to_scratch(192, clm); \
\
/* Do FFT levels 1,2,3 */ \
/* */ \
/* On input the 64-byte cache lines hold these data values: */ \
/*	0K	384K	768K	1152K */ \
/*	... */ \
/*	382K */ \
/* On output the 64-byte cache lines hold these data values: */ \
/*	0K	128K	256K	384K */ \
/*	2K	... */ \
/*	... */ \
/*	126K */ \
/*	512K	768K	1024K	1280K */ \
/*	514K */ \
/*	... */ \
/*	766K */ \
/* */ \
/* Do 256 six_reals_first_fft macros (FFT levels 1,2,3) */ \
/*	distance between fft data elements is 256K */ \
/*	do 64 iterations each processing 24 data values */ \
\
/* Do FFT levels 4,5 */ \
/* Values 0K-255K is real data, 256K-511K is semi-real data, 512K-1535K is */ \
/* complex data. */ \
/* */ \
/* On input the 64-byte cache lines hold these data values: */ \
/*	0K	128K	256K	384K */ \
/*	2K	... */ \
/*	... */ \
/*	126K */ \
/*	512K	768K	1024K	1280K */ \
/*	514K */ \
/*	... */ \
/*	766K */ \
/* On output the 64-byte cache lines hold these data values: */ \
/*	0K	64K	128K	192K */ \
/*	2K	... */ \
/*	... */ \
/*	62K	... */ \
/*	256K	... */ \
/*	... */ \
/* */ \
/* Do 32 eight_reals_fft_2 macros */ \
/* Do 32 nop_two_two_complex_fft_2 macros */ \
/*	distance between fft data elements is 64K */ \
/*	do 32 macros each processing 16 data values */ \
/* */ \
/* Do 128 four_complex_fft macros */ \
/*	distance between fft data elements is 64K */ \
/*	do 32 macros each processing 32 data values */ \
\
/* Do FFT levels 6 */ \
/* Values 0K-63K is real data, 64K-127K is semi-real data, 128K-1535K is */ \
/* complex data. */ \
/* */ \
/* On input the 64-byte cache lines hold these data values: */ \
/*	0K	64K	128K	192K */ \
/*	2K	... */ \
/*	... */ \
/*	62K	... */ \
/*	256K	... */ \
/*	... */ \
/* On output the 64-byte cache lines hold these data values: */ \
/*	0K	32K	64K	96K */ \
/*	2K	... */ \
/*	... */ \
/*	30K	... */ \
/*	128K	... */ \
/*	... */ \
/* */ \
/* Do 16 eight_reals_fft_1 macros */ \
/* Do 16 two_two_complex_fft macros */ \
/*	distance between fft data elements is 32K */ \
/*	do 16 macros each processing 16 data values */ \
/* */ \
/* Do 160 two_two_complex_fft macros */ \
/*	distance between fft data elements is 32K */ \
/*	do 80 macros each processing 16 data values */ \
\
yes_fft:start_timer(22); \
rsi = (uintptr_t)g->scratch_area;	/* Get address of scratch area */ \
rcx = (uintptr_t)g->data_addr; /* Load source address for prefetching */ \
pfing(rbp = (uintptr_t)g->data_prefetch);/* Set prefetch distance */ \
pfing(rbp -= rcx); \
pfing(rcx += 116*blkdst);		/* Create prefetch pointer */ \
rdi = (uintptr_t)g->sincos4;/* Load sin/cos pointer */ \
rbx = (uintptr_t)g->sincos3;/* Load sin/cos pointer */ \
/* 16 iterations */ \
xsix_reals_fft_preload; \
\
for(unsigned int loopA = 2; loopA; loopA--) { \
	for(unsigned int loopB = 8; loopB; loopB--) { \
		xtouchlines(2);			/* Preload TLBs */ \
		\
		for(unsigned int loopC = clm; loopC; loopC--) { \
			for(unsigned int loopD = 4; loopD; loopD--) { \
				for(unsigned int loopE = 2; loopE; loopE--) { \
					x3cl_six_reals_fft(rsi, 64, 8*clmblkdst8); \
				} \
				rsi += -2*64+2*clmblkdst8;	/* Next source pointer */ \
			} \
			\
			rsi += -4*2*clmblkdst8;	/* Restore source pointer */ \
			for(unsigned int loopD = 2; loopD; loopD--) { \
				for(unsigned int loopE = 2; loopE; loopE--) { \
					x2cl_eight_reals_fft_2(rsi, 64, 4*clmblkdst8); \
				} \
				rsi += -2*64+2*clmblkdst8;	/* Next source pointer */ \
			} \
			\
			rsi += -2*2*clmblkdst8+8*clmblkdst8;/* Next source pointer */ \
			for(unsigned int loopD = 2; loopD; loopD--) { \
				; \
				xloop1sp(b3b, 2*clmblkdst8, x4cl_four_complex_fft(rsi, 64, 4*clmblkdst8, 8*clmblkdst8, rdi)); \
			} \
			\
			rsi += -2*2*clmblkdst8-8*clmblkdst8;/* Restore source pointer */ \
			for(unsigned int loopD = 2; loopD; loopD--) { \
				x2cl_eight_reals_fft_1(rsi, 64, 2*clmblkdst8); \
			} \
			\
			rsi += -2*64+4*clmblkdst8;	/* Load source pointer */ \
			rbx += XMM_SCD;		/* Next sine/cosine pointer */ \
			/* 5 iterations */ \
			for(unsigned int loopD = 5; loopD; loopD--) { \
				for(unsigned int loopE = 2; loopE; loopE--) { \
					x2cl_two_complex_fft(rsi, 64, 2*clmblkdst8, rbx); \
				} \
				rsi += -2*64+4*clmblkdst8;	/* Next source pointer */ \
				rbx += 2*XMM_SCD;		/* Next sine/cosine pointer */ \
			} \
			rbx += -11*XMM_SCD;	/* Restore sine/cosine pointer */ \
			\
				EXEC248##clm(rsi += -6*4*clmblkdst8+2*64);/* Next source pointer */ \
			EXEC248##clm(pfing(rcx += -2*blkdst+128));	/* Next prefetch pointer */ \
		} \
		EXEC248##clm(rsi += -clm*2*64+clmblkdst);/* Next source pointer */ \
		EXEC248##clm(pfing(rcx += -clm*128+2*blkdst));	/* Next prefetch pointer */ \
		EXEC1##clm(rsi += -6*4*clmblkdst8+clmblkdst);/* Next source pointer */ \
		\
	} \
	rsi += -8*clmblkdst+clmblkdst8;/* Next source pointer */ \
} \
rsi += -2*clmblkdst8;	/* Restore source pointer */ \
end_timer(22); \
\
/* Do FFT levels 7,8 */ \
/* Values 0K-31K is real data, 32K-63K is semi-real data, 64K-1535K is */ \
/* complex data. */ \
/* */ \
/* On input the 64-byte cache lines hold these data values: */ \
/*	0K	32K	64K	96K */ \
/*	2K	... */ \
/*	... */ \
/*	30K	... */ \
/*	128K	... */ \
/*	... */ \
/* On output the 64-byte cache lines hold these data values: */ \
/*	0K	8K	16K	24K */ \
/*	2K	... */ \
/*	... */ \
/*	6K	... */ \
/*	32K	... */ \
/*	... */ \
\
/* Do 8 eight_reals_fft_2 macros */ \
/* Do 8 two_two_complex_fft macros */ \
/*	distance between fft data elements is 8K */ \
/*	do 4 macros each processing 32 data values */ \
\
start_timer(25); \
rdi = (uintptr_t)g->sincos2;/* Load sin/cos pointer */ \
/* 4 iterations of 2*clm */ \
for(unsigned int loopA = 4; loopA; loopA--) { \
	for(unsigned int loopB = 2*clm; loopB; loopB--) { \
		x4cl_eight_reals_fft_2(rsi, 64, 4*clmblkdst, clmblkdst8, rdi); \
	} \
	rsi += -2*clm*64+clmblkdst;/* Next source pointer */ \
} \
\
/* Do 176 four_complex_fft macros */ \
/*	distance between fft data elements is 8K */ \
/*	do 44 macros each processing 32 data values */ \
\
rsi += -4*clmblkdst+2*clmblkdst8;/* Load source pointer */ \
rdi += XMM_SCD;		/* Next sine/cosine pointer */ \
/* 11 iterations of 4 iters of 2*clm */ \
for(unsigned int loopA = 11; loopA; loopA--) { \
	for(unsigned int loopB = 4; loopB; loopB--) { \
		xtouch(xptr(rcx+(clm-1)*128));	/* Preload TLB */ \
		xloop1_clm##clm(b7b, clmblkdst, x4cl_four_complex_fft(rsi, 64, 4*clmblkdst, clmblkdst8, rdi);); \
	} \
	rsi += -4*clmblkdst+2*clmblkdst8;/* Next source pointer */ \
	rdi += 2*XMM_SCD;		/* Next sine/cosine pointer */ \
} \
rsi += -12*2*clmblkdst8;	/* Restore source pointer */ \
pfing(rcx += -192*blkdst);	/* Restore prefetch pointer */ \
end_timer(25); \
\
/* Do FFT levels 9,10 */ \
/* Values 0K-7K is real data, 8K-15K is semi-real data, 16K-1535K is */ \
/* complex data. */ \
/* */ \
/* On input the 64-byte cache lines hold these data values: */ \
/*	0K	8K	16K	24K */ \
/*	2K	... */ \
/*	... */ \
/*	6K	... */ \
/*	32K	... */ \
/*	... */ \
/* On output the 64-byte cache lines hold these data values: */ \
/*	0K	2K	4K	6K */ \
/*	8K	... */ \
/*	... */ \
\
/* Do 1 eight_reals_fft_2 macros */ \
/* Do 1 nop_two_two_complex_fft_2 macros */ \
/* Do 2 four_complex_fft macros */ \
/*	distance between fft data elements is 2K */ \
/*	do 1 macros each processing 32 data values */ \
\
start_timer(26); \
rdi = (uintptr_t)g->sincos1;/* Load sin/cos pointer */ \
for(unsigned int loopA = 2*clm; loopA; loopA--) { \
	g4cl_eight_reals_fft_2(rsi, 64, clmblkdst, 2*clmblkdst, rcx, 64, blkdst, 2*blkdst); \
} \
\
/* Do 188 four_complex_fft macros */ \
/*	distance between fft data elements is 2K */ \
/*	do 47 macros each processing 32 data values */ \
\
rsi += -2*clm*64+4*clmblkdst;/* Next source pointer */ \
rcx += -2*clm*64+4*blkdst;	/* Next dest pointer */ \
rdi += XMM_SCD;		/* Next sine/cosine pointer */ \
/* eax = 24*256+1; */ /* 47 iterations of 2*clm */ \
for(unsigned int loopA = 47; loopA; ) { \
	do{ \
		for(unsigned int loopC = clm; loopC; loopC--) { \
			for(unsigned int loopD = 2; loopD; loopD--) { \
				g4cl_four_complex_fft(rsi, 64, clmblkdst, 2*clmblkdst, rcx, 64, blkdst, 2*blkdst); \
				L2prefetch128(u8ptr(rcx-64+rbp)); \
				L2prefetch128(u8ptr(rcx-64+blkdst+rbp)); \
				pfing(rbp += -64+2*blkdst); \
			} \
			pfing(rbp += 128-4*blkdst); \
		} \
		rsi += -2*clm*64+4*clmblkdst;/* Next source pointer */ \
		rcx += -2*clm*64+4*blkdst;	/* Next dest pointer */ \
		rdi += 2*XMM_SCD;		/* Next sine/cosine pointer */ \
	}while(--loopA & 1); \
	rsi += -8*clmblkdst+clmblkdst8; \
} \
end_timer(26); \
\
/* Work on next 4 sets of 768 values */ \
\
endlp:	xpass1_get_next_block(pass2, c0b, b0b); /* Get next block and jump */ \
/* done: */	xfft_footer_hg2; \
\
\
\
/* Perform a 1536K-element negacyclic FFT.  This is done in two passes.  Pass */ \
/* 1 does 10 levels, pass 2 uses common code to perform the last 11 levels. */ \
\

#define hg_pass1sclevels10complex3(procname, xpass2name, pass2_levels, clm) \
 /* (fftname, pass2_levels, clmarg) */ \
\
/* set_FFT_constants(fftname, 0, 1, pass2_levels, clmarg); */ \
xfft_header_hg2(procname, xpass2name, pass2_levels, clm); \
\
/* Do the last 11 FFT levels */ \
/*	do 192 pass2_11_levels_complex */ \
\
/* pass2: */	start_timer(0); \
xpass2name(g); \
end_timer(0); \
if(g->ffttype == 1)		/* We're done if FFTing only */ \
	return; \
\
/* Do the inverse FFT */ \
\
pass1_inverse_fft_setup(c0b, procname); \
\
/* Do inverse FFT levels 9,10 */ \
/* On input the 64-byte cache lines hold these data values: */ \
/*	0K	2K	4K	6K */ \
/*	8K	10K	12K	14K */ \
/*	... */ \
/* On output the 64-byte cache lines hold these data values: */ \
/*	0K	4K	8K	12K */ \
/*	2K	... */ \
/*	16K	... */ \
/*	... */ \
\
/* Do 192 four_complex_unfft macros */ \
/*	distance between fft data elements is 2K */ \
/*	do 96 macros each processing 16 data values */ \
\
c0b:	start_timer(16); \
rsi = (uintptr_t)g->data_addr;/* Get source address */ \
rcx = (uintptr_t)g->scratch_area;	/* Get address of scratch area */ \
rdi = (uintptr_t)g->sincos1;/* Load sin/cos pointer */ \
/* 24 iterations of 4 iters of 2*clm */ \
pfing(rbp = (uintptr_t)g->data_prefetch);/* Set prefetch distance */ \
pfing(rbp -= rsi); \
for(unsigned int loopA = 24; loopA; loopA--) { \
	for(unsigned int loopB = 4; loopB; loopB--) { \
		for(unsigned int loopC = 2*clm; loopC; loopC--) { \
			g2cl_four_complex_unfft(rsi, 64, blkdst, rcx, 64, clmblkdst); \
		} \
		rsi += -2*clm*64+2*blkdst;	/* Next source pointer */ \
		rcx += -2*clm*64+2*clmblkdst;/* Next dest pointer */ \
		rdi += XMM_SCD;		/* Next sine/cosine pointer */ \
	} \
	rcx += -8*clmblkdst+clmblkdst8; \
} \
rcx = rsi-96*2*blkdst;	/* Create prefetch pointer */ \
rsi = (uintptr_t)g->scratch_area;	/* Restore source pointer */ \
end_timer(16); \
\
/* Do inverse FFT levels 7,8 */ \
/* On input the 64-byte cache lines hold these data values: */ \
/*	0K	4K	8K	12K */ \
/*	2K	... */ \
/*	16K	... */ \
/*	... */ \
/* On output the 64-byte cache lines hold these data values: */ \
/*	0K	16K	32K	48K */ \
/*	2K	... */ \
/*	... */ \
/*	14K */ \
/*	64K	... */ \
/*	... */ \
\
/* Do 192 four_complex_unfft macros */ \
/*	distance between fft data elements is 8K */ \
/*	do 48 macros each processing 32 data values */ \
\
start_timer(17); \
rdi = (uintptr_t)g->sincos2;/* Load sin/cos pointer */ \
/* 24 iterations of 2*2*clm */ \
for(unsigned int loopA = 24; loopA; loopA--) { \
	xtouch(xptr(rcx+(clm-1)*128));	/* Preload TLB */ \
	xloop2_clm##clm(c7b, clmblkdst, x4cl_four_complex_unfft(rsi, 64, 2*clmblkdst, 4*clmblkdst, rdi));		/* 2*2*clm iterations */ \
	rsi += -2*clmblkdst+clmblkdst8;/* Next source pointer */ \
	rdi += XMM_SCD;		/* Next sine/cosine pointer */ \
} \
rsi += -24*clmblkdst8;	/* Restore source pointer */ \
end_timer(17); \
\
/* Do inverse FFT levels 5,6 */ \
/* On input the 64-byte cache lines hold these data values: */ \
/*	0K	16K	32K	48K */ \
/*	2K	... */ \
/*	... */ \
/*	14K */ \
/*	64K	... */ \
/*	... */ \
/* On output the 64-byte cache lines hold these data values: */ \
/*	0K	64K	128K	192K */ \
/*	2K	... */ \
/*	... */ \
/*	62K */ \
/*	256K	... */ \
/*	... */ \
\
/* Do 192 four_complex_unfft macros */ \
/*	distance between fft data elements is 32K */ \
/*	do 48 macros each processing 32 data values */ \
\
start_timer(18); \
rdi = (uintptr_t)g->sincos3;/* Load sin/cos pointer */ \
/* 6 iterations of 4 iters of 2*2*clm */ \
for(unsigned int loopA = 6; loopA; loopA--) { \
	for(unsigned int loopB = 4; loopB; loopB--) { \
		xtouch(xptr(rcx+(clm-1)*128));	/* Preload TLB */ \
		xloop2_clm##clm(c5b, clmblkdst, x4cl_four_complex_unfft(rsi, 64, clmblkdst8, 2*clmblkdst8, rdi));		/* 2 * 2*clm iterations */ \
	} \
	rsi += -4*2*clmblkdst+4*clmblkdst8;/* Next source pointer */ \
	rdi += XMM_SCD;		/* Next sine/cosine pointer */ \
} \
rsi += -6*4*clmblkdst8;	/* Restore source pointer */ \
end_timer(18); \
\
/* Do inverse FFT levels 4 */ \
/* On input the 64-byte cache lines hold these data values: */ \
/*	0K	64K	128K	192K */ \
/*	2K	... */ \
/*	... */ \
/*	62K */ \
/*	256K	... */ \
/*	... */ \
/* On output the 64-byte cache lines hold these data values: */ \
/*	0K	128K	256K	384K */ \
/*	2K	... */ \
/*	... */ \
/*	126K */ \
/*	512K	... */ \
/*	... */ \
\
/* Do 192 two_two_complex_unfft macros */ \
/*	distance between fft data elements is 64K */ \
/*	do 96 macros each processing 16 data values */ \
\
start_timer(19); \
rdi = (uintptr_t)g->sincos4;/* Load sin/cos pointer */ \
/* 3 iters of 4 iters 2 of 4*2*clm */ \
for(unsigned int loopA = 3; loopA; loopA--) { \
	for(unsigned int loopB = 4; loopB; loopB--) { \
		for(unsigned int loopC = 2; loopC; loopC--) { \
			xtouch(xptr(rcx+(clm-1)*128));	/* Preload TLB */ \
			xloop4_clm##clm(c3b, clmblkdst, x2cl_two_complex_unfft(rsi, 64, 4*clmblkdst8));		/* 4* 2*clm iterations */ \
		} \
		rsi += -2*4*clmblkdst+clmblkdst8;/* Next source pointer */ \
	} \
	rsi += -4*clmblkdst8+8*clmblkdst8;/* Next source pointer */ \
	rdi += XMM_SCD;		/* Next sine/cosine pointer */ \
} \
rsi += -3*8*clmblkdst8;	/* Restore source pointer */ \
end_timer(19); \
\
/* Do inverse FFT levels 1,2,3 */ \
/* On input the 64-byte cache lines hold these data values: */ \
/*	0K	128K	256K	384K */ \
/*	2K	... */ \
/*	... */ \
/*	126K */ \
/*	512K	768K	1024K	1280K */ \
/*	514K */ \
/*	... */ \
/*	766K */ \
/* On output the 64-byte cache lines hold these data values: */ \
/*	0K	384K	768K	1152K */ \
/*	... */ \
/*	382K */ \
\
/* Do 256 three_complex_last_unfft macros */ \
/*	distance between fft data elements is 256K */ \
/*	do 64 macros each processing 24 data values */ \
\
start_timer(20); \
rdi = (uintptr_t)g->plus1_premults;/* Address of premultiplier table */ \
/* 8 iters of 4 iters of 2*2*clm */ \
for(unsigned int loopA = 8; loopA; loopA--) { \
	for(unsigned int loopB = 4; loopB; loopB--) { \
		xtouch(xptr(rcx+(clm-1)*128));	/* Preload TLB */ \
		xloop2edi_clm##clm(c1b, clmblkdst, 192, x3cl_three_complex_last_unfft(rsi, 64, 8*clmblkdst8)); /* 2 * 2*clm iterations */ \
	} \
	rsi += -4*2*clmblkdst+clmblkdst8;/* Next source pointer */ \
} \
end_timer(20); \
\
/* Normalize these values */ \
\
pass1_normalize(1, yes_fft); \
\
/* Copy data from the scratch area back to the FFT data area */ \
\
copy_scratch_data_to_fft(192, clm); \
goto	endlp;			/* Skip FFT code if flag not set */ \
\
/* Fire up auxiliary threads to do pass 1 forward FFT */ \
\
/* pass1: */	pass1_forward_fft_setup(b0b); \
\
/* Copy data from the FFT data area back to the scratch area */ \
\
b0b:	copy_fft_data_to_scratch(192, clm); \
\
/* Do FFT levels 1,2,3 */ \
/* */ \
/* On input the 64-byte cache lines hold these data values: */ \
/*	0K	384K	768K	1152K */ \
/*	... */ \
/*	382K */ \
/* On output the 64-byte cache lines hold these data values: */ \
/*	0K	128K	256K	384K */ \
/*	2K	... */ \
/*	... */ \
/*	126K */ \
/*	512K	768K	1024K	1280K */ \
/*	514K */ \
/*	... */ \
/*	766K */ \
\
/* Do 256 three_complex_first_fft macros (FFT levels 1,2,3) */ \
/*	distance between fft data elements is 256K */ \
/*	do 64 iterations each processing 24 data values */ \
\
yes_fft:start_timer(22); \
rsi = (uintptr_t)g->scratch_area;	/* Get address of scratch area */ \
rdi = (uintptr_t)g->plus1_premults;/* Address of premultiplier table */ \
rcx = (uintptr_t)g->data_addr; /* Load source address for prefetching */ \
pfing(rbp = (uintptr_t)g->data_prefetch);/* Set prefetch distance */ \
pfing(rbp -= rcx); \
pfing(rcx += 104*blkdst);		/* Create prefetching pointer */ \
/* 8 iters of 2 iters of 4*2*clm */ \
for(unsigned int loopA = 8; loopA; loopA--) { \
	for(unsigned int loopB = 2; loopB; loopB--) { \
		xtouch(xptr(rcx+(clm-1)*128));	/* Preload TLB */ \
		xloop4edi_clm##clm(b1b, clmblkdst, 192, x3cl_three_complex_fft(rsi, 64, 8*clmblkdst8)); /* 4 * 2*clm iterations */ \
	} \
	rsi += -2*4*clmblkdst+clmblkdst8; \
} \
rsi += -8*clmblkdst8;	/* Restore source pointer */ \
end_timer(22); \
\
/* Do FFT levels 4 */ \
/* */ \
/* On input the 64-byte cache lines hold these data values: */ \
/*	0K	128K	256K	384K */ \
/*	2K	... */ \
/*	... */ \
/*	126K */ \
/*	512K	768K	1024K	1280K */ \
/*	514K */ \
/*	... */ \
/*	766K */ \
/* On output the 64-byte cache lines hold these data values: */ \
/*	0K	128K	256K	384K */ \
/*	2K	... */ \
/*	... */ \
/*	126K	... */ \
/*	512K	... */ \
/*	... */ \
\
/* Do 192 two_two_complex_fft macros */ \
/*	distance between fft data elements is 64K */ \
/*	do 96 macros each processing 16 data values */ \
\
start_timer(23); \
rdi = (uintptr_t)g->sincos4;/* Load sin/cos pointer */ \
/* 4 iters of 2 of 4*2*clm */ \
for(unsigned int loopA = 4; loopA; loopA--) { \
	for(unsigned int loopB = 2; loopB; loopB--) { \
		xtouch(xptr(rcx+(clm-1)*128));	/* Preload TLB */ \
		xloop4_clm##clm(b2b, clmblkdst, x2cl_two_complex_fft_in_place(rsi, 64, 4*clmblkdst8, rdi));		/* 4 * 2*clm iterations */ \
	} \
	rsi += -2*4*clmblkdst+clmblkdst8; \
} \
rsi += -4*clmblkdst8+8*clmblkdst8;/* Next source pointer */ \
rdi += XMM_SCD;		/* Next sine/cosine pointer */ \
/* 8 iters of 2 of 4*2*clm */ \
for(unsigned int loopA = 8; loopA; loopA--) { \
	for(unsigned int loopB = 2; loopB; loopB--) { \
		xtouch(xptr(rcx+(clm-1)*128));	/* Preload TLB */ \
		xloop4_clm##clm(b3b, clmblkdst, x2cl_two_complex_fft(rsi, 64, 8*clmblkdst8, rdi));		/* 4 * 2*clm iterations */ \
	} \
	rsi += -2*4*clmblkdst+clmblkdst8; \
} \
rsi += -8*clmblkdst8-8*clmblkdst8;/* Restore source pointer */ \
end_timer(23); \
\
/* Do FFT levels 5,6 */ \
/* */ \
/* On input the 64-byte cache lines hold these data values: */ \
/*	0K	128K	256K	384K */ \
/*	2K	... */ \
/*	... */ \
/*	126K	... */ \
/*	512K	... */ \
/*	... */ \
/* On output the 64-byte cache lines hold these data values: */ \
/*	0K	32K	64K	96K */ \
/*	2K	... */ \
/*	... */ \
/*	30K	... */ \
/*	128K	... */ \
/*	... */ \
\
/* Do 192 four_complex_fft macros */ \
/*	distance between fft data elements is 32K */ \
/*	do 48 macros each processing 32 data values */ \
\
start_timer(24); \
rdi = (uintptr_t)g->sincos3;/* Load sin/cos pointer */ \
/* 3 iters of 2 iters of 4 of 2*2*clm */ \
for(unsigned int loopA = 3; loopA; loopA--) { \
	for(unsigned int loopB = 2; loopB; loopB--) { \
		for(unsigned int loopC = 4; loopC; loopC--) { \
			xtouch(xptr(rcx+(clm-1)*128));	/* Preload TLB */ \
			xloop2_clm##clm(b5b, clmblkdst, x4cl_four_complex_fft(rsi, 64, 2*clmblkdst8, 4*clmblkdst8, rdi));		/* 2 * 2*clm iterations */ \
		} \
		rsi += -4*2*clmblkdst+clmblkdst8; \
	} \
	rsi += -2*clmblkdst8+8*clmblkdst8;/* Next source pointer */ \
	rdi += 2*XMM_SCD;		/* Next sine/cosine pointer */ \
} \
rsi += -3*8*clmblkdst8;	/* Restore source pointer */ \
end_timer(24); \
\
/* Do FFT levels 7,8 */ \
/* */ \
/* On input the 64-byte cache lines hold these data values: */ \
/*	0K	32K	64K	96K */ \
/*	2K	... */ \
/*	... */ \
/*	30K	... */ \
/*	128K	... */ \
/*	... */ \
/* On output the 64-byte cache lines hold these data values: */ \
/*	0K	8K	16K	24K */ \
/*	2K	... */ \
/*	... */ \
/*	6K	... */ \
/*	32K	... */ \
/*	... */ \
\
/* Do 192 four_complex_fft macros */ \
/*	distance between fft data elements is 8K */ \
/*	do 48 macros each processing 32 data values */ \
\
start_timer(25); \
rdi = (uintptr_t)g->sincos2;/* Load sin/cos pointer */ \
/* 12 iterations of 2 iters of 2*2*clm */ \
for(unsigned int loopA = 12; loopA; loopA--) { \
	for(unsigned int loopB = 2; loopB; loopB--) { \
		xtouch(xptr(rcx+(clm-1)*128));	/* Preload TLB */ \
		xloop2_clm##clm(b7b, clmblkdst, x4cl_four_complex_fft(rsi, 64, 4*clmblkdst, clmblkdst8, rdi));		/* 2*2*clm iterations */ \
	} \
	rsi += -2*2*clmblkdst+2*clmblkdst8;/* Next source pointer */ \
	rdi += 2*XMM_SCD;		/* Next sine/cosine pointer */ \
} \
rsi += -12*2*clmblkdst8;	/* Restore source pointer */ \
pfing(rcx += -192*blkdst);	/* Restore prefetch pointer */ \
end_timer(25); \
\
/* Do FFT levels 9,10 */ \
/* */ \
/* On input the 64-byte cache lines hold these data values: */ \
/*	0K	8K	16K	24K */ \
/*	2K	... */ \
/*	... */ \
/*	6K	... */ \
/*	32K	... */ \
/*	... */ \
/* On output the 64-byte cache lines hold these data values: */ \
/*	0K	2K	4K	6K */ \
/*	8K	... */ \
/*	... */ \
\
/* Do 192 four_complex_fft macros */ \
/*	distance between fft data elements is 2K */ \
/*	do 48 macros each processing 32 data values */ \
\
start_timer(26); \
rdi = (uintptr_t)g->sincos1;/* Load sin/cos pointer */ \
/* 24 iterations of 2 iters of 2*clm */ \
for(unsigned int loopA = 24; loopA; loopA--) { \
	for(unsigned int loopB = 2; loopB; loopB--) { \
		for(unsigned int loopC = clm; loopC; loopC--) { \
			for(unsigned int loopD = 2; loopD; loopD--) { \
				g4cl_four_complex_fft(rsi, 64, clmblkdst, 2*clmblkdst, rcx, 64, blkdst, 2*blkdst); \
				L2prefetch128(u8ptr(rcx-64+rbp)); \
				L2prefetch128(u8ptr(rcx-64+blkdst+rbp)); \
				pfing(rbp += -64+2*blkdst); \
			} \
			pfing(rbp += 128-4*blkdst); \
		} \
	rsi += -2*clm*64+4*clmblkdst;/* Next source pointer */ \
		rcx += -2*clm*64+4*blkdst;	/* Next dest pointer */ \
		rdi += 2*XMM_SCD;		/* Next sine/cosine pointer */ \
	} \
	rsi += -8*clmblkdst+clmblkdst8; \
} \
end_timer(26); \
\
/* Work on next 4 sets of 768 values */ \
\
endlp:	xpass1_get_next_block(pass2, c0b, b0b); /* Get next block and jump */ \
/* done: */	xfft_footer_hg2; \
\
\
\
/* Perform a 1792K-element FFT.  This is done in two passes.  Pass 1 */ \
/* does 10 levels, pass 2 uses common code to perform the last 11 levels. */ \
\

#define hg_pass1sclevels10pfa7(procname, xpass2name, pass2_levels, clm) \
 /* (fftname, pass2_levels, clmarg) */ \
\
/* set_FFT_constants(fftname, 0, 0, pass2_levels, clmarg); */ \
xfft_header_hg2(procname, xpass2name, pass2_levels, clm); \
\
/* Do the last 11 FFT levels */ \
/*	do 1 pass2_11_levels_real and 223 pass2_11_levels_complex */ \
\
/* pass2: */	start_timer(0); \
xpass2name(g); \
end_timer(0); \
if(g->ffttype == 1)		/* We're done if FFTing only */ \
	return; \
\
/* Do the inverse FFT */ \
\
pass1_inverse_fft_setup(c0b, procname); \
\
/* Do inverse FFT levels 9,10 */ \
/* On input the 64-byte cache lines hold these data values: */ \
/*	0K	2K	4K	6K */ \
/*	8K	10K	12K	14K */ \
/*	... */ \
/* On output the 64-byte cache lines hold these data values: */ \
/*	0K	4K	8K	12K */ \
/*	2K	... */ \
/*	16K	... */ \
/*	... */ \
\
/* Do 1 eight_reals_unfft_2 macros */ \
/* Do 1 nop_two_two_complex_unfft_2 macros */ \
/*	distance between fft data elements is 2K */ \
/*	do 1 macros each processing 16 data values */ \
\
c0b:	start_timer(16); \
rsi = (uintptr_t)g->data_addr;/* Get source address */ \
rcx = (uintptr_t)g->scratch_area;	/* Get address of scratch area */ \
pfing(rbp = (uintptr_t)g->data_prefetch);/* Set prefetch distance */ \
pfing(rbp -= rsi); \
for(unsigned int loopA = 2*clm; loopA; loopA--) { \
	g2cl_eight_reals_unfft_2(rsi, 64, blkdst, rcx, 64, clmblkdst); \
} \
\
/* Do 222 four_complex_unfft macros */ \
/*	distance between fft data elements is 2K */ \
/*	do 111 macros each processing 16 data values */ \
\
rsi += -2*clm*64+2*blkdst;	/* Next source pointer */ \
rcx += -2*clm*64+2*clmblkdst;/* Next dest pointer */ \
rdi = (uintptr_t)g->sincos1;/* Load sin/cos pointer */ \
/* al = 256/4+28; */ /* 111 iterations of 2*clm */ \
for(unsigned int loopA = 111; loopA; ) { \
	for(unsigned int loopB = 4; loopB; loopB--) { \
		for(unsigned int loopC = 2*clm; loopC; loopC--) { \
			g2cl_four_complex_unfft(rsi, 64, blkdst, rcx, 64, clmblkdst); \
		} \
		rsi += -2*clm*64+2*blkdst;	/* Next source pointer */ \
		rcx += -2*clm*64+2*clmblkdst;/* Next dest pointer */ \
		rdi += XMM_SCD;		/* Next sine/cosine pointer */ \
	}while(--loopA & 3); \
	rcx += -8*clmblkdst+clmblkdst8; \
} \
rcx = rsi-112*2*blkdst;	/* Create prefetch pointer */ \
rsi = (uintptr_t)g->scratch_area;	/* Restore source pointer */ \
end_timer(16); \
\
/* Do inverse FFT levels 7,8 */ \
/* On input the 64-byte cache lines hold these data values: */ \
/*	0K	4K	8K	12K */ \
/*	2K	... */ \
/*	16K	... */ \
/*	... */ \
/* On output the 64-byte cache lines hold these data values: */ \
/*	0K	16K	32K	48K */ \
/*	2K	... */ \
/*	... */ \
/*	14K */ \
/*	64K	... */ \
/*	... */ \
\
/* Do 4 eight_reals_unfft_2 macros */ \
/* Do 4 nop_two_two_complex_unfft_2 macros */ \
/*	distance between fft data elements is 8K */ \
/*	do 2 macros each processing 32 data values */ \
\
start_timer(17); \
/* 2 iterations of 2*clm */ \
for(unsigned int loopA = 2; loopA; loopA--) { \
	xtouch(xptr(rcx+(clm-1)*128));	/* Preload TLB */ \
	xloop1_clm##clm(c6b, clmblkdst, x4cl_eight_reals_unfft_2(rsi, 64, 2*clmblkdst, 4*clmblkdst)); \
} \
\
/* Do 216 four_complex_unfft macros */ \
/*	distance between fft data elements is 8K */ \
/*	do 54 macros each processing 32 data values */ \
\
rsi += -2*clmblkdst+clmblkdst8;/* Next source pointer */ \
rdi = (uintptr_t)g->sincos2;/* Load sin/cos pointer */ \
/* 27 iterations of 2 iters of 2*clm */ \
for(unsigned int loopA = 27; loopA; loopA--) { \
	for(unsigned int loopB = 2; loopB; loopB--) { \
		xtouch(xptr(rcx+(clm-1)*128));	/* Preload TLB */ \
		xloop1_clm##clm(c7b, clmblkdst, x4cl_four_complex_unfft(rsi, 64, 2*clmblkdst, 4*clmblkdst, rdi)); \
	} \
	rsi += -2*clmblkdst+clmblkdst8;/* Next source pointer */ \
	rdi += XMM_SCD;		/* Next sine/cosine pointer */ \
} \
rsi += -28*clmblkdst8;	/* Restore source pointer */ \
end_timer(17); \
\
/* Do inverse FFT levels 6 */ \
/* On input the 64-byte cache lines hold these data values: */ \
/*	0K	16K	32K	48K */ \
/*	2K	... */ \
/*	... */ \
/*	14K */ \
/*	64K	... */ \
/*	... */ \
/* On output the 64-byte cache lines hold these data values: */ \
/*	0K	32K	64K	96K */ \
/*	2K	... */ \
/*	... */ \
/*	30K */ \
/*	128K	... */ \
/*	... */ \
\
/* Do 16 eight_reals_unfft_1 macros */ \
/*	distance between fft data elements is 32K */ \
/*	do 8 macros each processing 16 data values */ \
/* */ \
/* Do 208 four_complex_unfft macros */ \
/*	distance between fft data elements is 32K */ \
/*	do 104 macros each processing 16 data values */ \
\
/* Do inverse FFT levels 4,5 */ \
/* On input the 64-byte cache lines hold these data values: */ \
/*	0K	32K	64K	96K */ \
/*	2K	... */ \
/*	... */ \
/*	30K */ \
/*	128K	... */ \
/*	... */ \
/* On output the 64-byte cache lines hold these data values: */ \
/*	0K	64K	128K	192K */ \
/*	2K	... */ \
/*	... */ \
/*	62K */ \
/*	256K	384K	512K	640K */ \
/*	... */ \
/* */ \
/* Do 32 eight_reals_unfft_2 macros */ \
/*	distance between fft data elements is 64K */ \
/*	do 16 macros each processing 16 data values */ \
/* */ \
/* Do 192 four_complex_unfft macros */ \
/*	distance between fft data elements is 64K */ \
/*	do 48 macros each processing 32 data values */ \
/* */ \
/* Do inverse FFT levels 1,2,3 */ \
/* On input the 64-byte cache lines hold these data values: */ \
/*	0K	64K	128K	192K */ \
/*	2K	... */ \
/*	... */ \
/*	62K */ \
/*	256K	384K	512K	640K */ \
/*	... */ \
/* On output the 64-byte cache lines hold these data values: */ \
/*	0K	64K	896K	960K */ \
/*	2K	... */ \
/*	... */ \
/*	62K */ \
/*	128K	256K	1024K	1152K */ \
/*	130K */ \
/*	... */ \
/*	254K */ \
/*	384K	640K	1280K	1536K */ \
/*	386K */ \
/*	... */ \
/*	638K */ \
/* */ \
/* Do 256 seven_reals_last_unfft macros */ \
/*	distance between fft data elements is 256K */ \
/*	do 32 macros each processing 56 data values */ \
\
start_timer(18); \
rdi = (uintptr_t)g->sincos3;/* Load sin/cos pointer */ \
rbx = (uintptr_t)g->sincos4;/* Load sin/cos pointer */ \
/* 8 iterations */ \
xseven_reals_unfft_preload; \
\
for(unsigned int loopA = 8; loopA; loopA--) { \
	xtouchlines(14);			/* Preload TLBs */ \
	\
	for(unsigned int loopB = clm; loopB; loopB--) { \
		for(unsigned int loopC = 2; loopC; loopC--) { \
			x2cl_eight_reals_unfft_1(rsi, 64, clmblkdst8); \
		} \
		\
		rsi += -2*64+2*clmblkdst8;	/* Next source pointer */ \
		/* 13 iterations */ \
		for(unsigned int loopC = 13; loopC; loopC--) { \
			for(unsigned int loopD = 2; loopD; loopD--) { \
				x2cl_two_complex_unfft(rsi, 64, clmblkdst8); \
			} \
			rsi += -2*64+2*clmblkdst8;	/* Next source pointer */ \
			rdi += XMM_SCD;		/* Next sine/cosine pointer */ \
		} \
		rdi += -13*XMM_SCD;	/* Restore sine/cosine pointer */ \
		\
		rsi += -14*2*clmblkdst8;	/* Restore source pointer */ \
		for(unsigned int loopC = 2; loopC; loopC--) { \
			for(unsigned int loopD = 2; loopD; loopD--) { \
				x2cl_half_eight_reals_unfft_2(rsi, 64, 2*clmblkdst8); \
			} \
			rsi += -2*64+clmblkdst8; \
		} \
		\
		rsi += -2*clmblkdst8+4*clmblkdst8;/* Next source pointer */ \
		/* 3 iterations */ \
		for(unsigned int loopC = 3; loopC; loopC--) { \
			for(unsigned int loopD = 2; loopD; loopD--) { \
				xloop1sp(c3b, clmblkdst8, x4cl_four_complex_unfft(rsi, 64, 2*clmblkdst8, 4*clmblkdst8, rbx)); \
			} \
			rsi += -2*clmblkdst8+8*clmblkdst8; \
			rbx += XMM_SCD;		/* Next sine/cosine pointer */ \
		} \
		rbx += -3*XMM_SCD;		/* Restore sine/cosine pointer */ \
		\
		rsi += -3*8*clmblkdst8-4*clmblkdst8;/* Restore source pointer */ \
		for(unsigned int loopC = 4; loopC; loopC--) { \
			xloop1prefetch2sp(c1b, clmblkdst8, x7cl_seven_reals_last_unfft(rsi, 64, 4*clmblkdst8)); \
		} \
		\
		EXEC248##clm(rsi += -4*clmblkdst8+2*64);	/* Next source pointer */ \
		EXEC248##clm(pfing(rcx += -14*blkdst+128));	/* Next prefetch pointer */ \
	} \
	EXEC248##clm(rsi += -clm*2*64+clmblkdst);/* Next source pointer */ \
	EXEC248##clm(pfing(rcx += -clm*128+14*blkdst));	/* Next prefetch pointer */ \
	EXEC1##clm(rsi += -4*clmblkdst8+clmblkdst);/* Next source pointer */ \
	\
} \
end_timer(18); \
\
/* Normalize these values */ \
\
pass1_normalize(1, yes_fft); \
\
/* Copy data from the scratch area back to the FFT data area */ \
\
copy_scratch_data_to_fft(224, clm); \
goto	endlp;			/* Skip FFT code if flag not set */ \
\
/* Fire up auxiliary threads to do pass 1 forward FFT */ \
\
/* pass1: */	pass1_forward_fft_setup(b0b); \
\
/* Copy data from the FFT data area back to the scratch area */ \
\
b0b:	copy_fft_data_to_scratch(224, clm); \
\
/* Do FFT levels 1,2,3 */ \
/* */ \
/* On input the 64-byte cache lines hold these data values: */ \
/*	0K	64K	896K	960K */ \
/*	2K	... */ \
/*	... */ \
/*	62K */ \
/*	128K	256K	1024K	1152K */ \
/*	130K */ \
/*	... */ \
/*	254K */ \
/*	384K	640K	1280K	1536K */ \
/*	386K */ \
/*	... */ \
/*	638K */ \
/* On output the 64-byte cache lines hold these data values: */ \
/*	0K	64K	128K	192K */ \
/*	2K	... */ \
/*	... */ \
/*	62K */ \
/*	256K	384K	512K	640K */ \
/*	258K */ \
/*	... */ \
/*	382K */ \
/*	768K	1024K	1280K	1536K */ \
/*	770K */ \
/*	... */ \
/*	1022K */ \
/* */ \
/* Do 256 seven_reals_first_fft macros (FFT levels 1,2,3) */ \
/*	distance between fft data elements is 256K */ \
/*	do 32 iterations each processing 56 data values */ \
\
/* Do FFT levels 4,5 */ \
/* Values 0K-255K is real data, 256K-511K is semi-real data, 512K-1791K is */ \
/* complex data. */ \
/* */ \
/* On input the 64-byte cache lines hold these data values: */ \
/*	0K	64K	128K	192K */ \
/*	2K	... */ \
/*	... */ \
/*	62K */ \
/*	256K	384K	512K	640K */ \
/*	258K */ \
/*	... */ \
/*	382K */ \
/*	768K	1024K	1280K	1536K */ \
/*	770K */ \
/*	... */ \
/*	1022K */ \
/* On output the 64-byte cache lines hold these data values: */ \
/*	0K	64K	128K	192K */ \
/*	2K	... */ \
/*	... */ \
/*	62K	... */ \
/*	256K	... */ \
/*	... */ \
/* */ \
/* Do 32 eight_reals_fft_2 macros */ \
/*	distance between fft data elements is 64K */ \
/*	do 32 macros each processing 8 data values */ \
/* */ \
/* Do 192 four_complex_fft macros */ \
/*	distance between fft data elements is 64K */ \
/*	do 48 macros each processing 32 data values */ \
\
/* Do FFT levels 6 */ \
/* Values 0K-63K is real data, 64K-127K is semi-real data, 128K-1791K is */ \
/* complex data. */ \
/* */ \
/* On input the 64-byte cache lines hold these data values: */ \
/*	0K	64K	128K	192K */ \
/*	2K	... */ \
/*	... */ \
/*	62K	... */ \
/*	256K	... */ \
/*	... */ \
/* On output the 64-byte cache lines hold these data values: */ \
/*	0K	32K	64K	96K */ \
/*	2K	... */ \
/*	... */ \
/*	30K	... */ \
/*	128K	... */ \
/*	... */ \
/* */ \
/* Do 16 eight_reals_fft_1 macros */ \
/* Do 16 two_two_complex_fft macros */ \
/*	distance between fft data elements is 32K */ \
/*	do 16 macros each processing 16 data values */ \
/* */ \
/* Do 192 two_two_complex_fft macros */ \
/*	distance between fft data elements is 32K */ \
/*	do 96 macros each processing 16 data values */ \
\
yes_fft:start_timer(22); \
rsi = (uintptr_t)g->scratch_area;	/* Get address of scratch area */ \
rcx = (uintptr_t)g->data_addr; /* Load source address for prefetching */ \
pfing(rbp = (uintptr_t)g->data_prefetch);/* Set prefetch distance */ \
pfing(rbp -= rcx); \
pfing(rcx += 72*blkdst);		/* Create prefetching pointer */ \
rdi = (uintptr_t)g->sincos4;/* Load sin/cos pointer */ \
rbx = (uintptr_t)g->sincos3;/* Load sin/cos pointer */ \
/* 16 iterations */ \
xseven_reals_fft_preload; \
\
for(unsigned int loopA = 2; loopA; loopA--) { \
	for(unsigned int loopB = 8; loopB; loopB--) { \
		xtouchlines(6);			/* Preload TLBs */ \
		\
		for(unsigned int loopC = clm; loopC; loopC--) { \
			for(unsigned int loopD = 2; loopD; loopD--) { \
				xloop1prefetch2sp(b1b, 2*clmblkdst8, x7cl_seven_reals_fft(rsi, 64, 4*clmblkdst8)); \
			} \
			\
			rsi += -2*2*clmblkdst8;	/* Restore source pointer */ \
			for(unsigned int loopD = 2; loopD; loopD--) { \
				for(unsigned int loopE = 2; loopE; loopE--) { \
					x1cl_half_eight_reals_fft_2(rsi, 64); \
				} \
				rsi += -2*64+2*clmblkdst8;	/* Next source pointer */ \
			} \
			\
			/*	bump	rsi, -2*2*clmblkdst8+4*clmblkdst8;; Next source pointer */ \
			for(unsigned int loopD = 2; loopD; loopD--) { \
				for(unsigned int loopE = 2; loopE; loopE--) { \
					x2cl_four_complex_fft(rsi, 64, 4*clmblkdst8); \
				} \
				rsi += -2*64+2*clmblkdst8;	/* Next source pointer */ \
			} \
			\
			rsi += -2*2*clmblkdst8+8*clmblkdst8;/* Next source pointer */ \
			rdi += XMM_SCD;		/* Next sine/cosine pointer */ \
			for(unsigned int loopD = 2; loopD; loopD--) { \
				xloop1sp(b3b, 2*clmblkdst8, x4cl_four_complex_fft(rsi, 64, 4*clmblkdst8, 8*clmblkdst8, rdi)); \
			} \
			rdi += -XMM_SCD;		/* Restore sine/cosine pointer */ \
			\
			rsi += -2*2*clmblkdst8-12*clmblkdst8;/* Restore src pointer */ \
			for(unsigned int loopD = 2; loopD; loopD--) { \
				x2cl_eight_reals_fft_1(rsi, 64, 2*clmblkdst8); \
			} \
			\
			rsi += -2*64+4*clmblkdst8;	/* Next source pointer */ \
			rbx += XMM_SCD;		/* Next sine/cosine pointer */ \
			/* 6 iterations */ \
			for(unsigned int loopD = 6; loopD; loopD--) { \
				for(unsigned int loopE = 2; loopE; loopE--) { \
					x2cl_two_complex_fft(rsi, 64, 2*clmblkdst8, rbx); \
				} \
				rsi += -2*64+4*clmblkdst8;	/* Next source pointer */ \
				rbx += 2*XMM_SCD;		/* Next sine/cosine pointer */ \
			} \
			rbx += -13*XMM_SCD;	/* Restore sine/cosine pointer */ \
			\
				EXEC248##clm(rsi += -7*4*clmblkdst8+2*64);/* Next source pointer */ \
			EXEC248##clm(pfing(rcx += -6*blkdst+128));	/* Next prefetch pointer */ \
		} \
		EXEC248##clm(rsi += -clm*2*64+clmblkdst);/* Next source pointer */ \
		EXEC248##clm(pfing(rcx += -clm*128+6*blkdst));	/* Next prefetch pointer */ \
		EXEC1##clm(rsi += -7*4*clmblkdst8+clmblkdst);/* Next source pointer */ \
		\
	} \
	rsi += -8*clmblkdst+clmblkdst8;/* Next source pointer */ \
} \
rsi += -2*clmblkdst8;	/* Restore source pointer */ \
end_timer(22); \
\
/* Do FFT levels 7,8 */ \
/* Values 0K-31K is real data, 32K-63K is semi-real data, 64K-1791K is */ \
/* complex data. */ \
/* */ \
/* On input the 64-byte cache lines hold these data values: */ \
/*	0K	32K	64K	96K */ \
/*	2K	... */ \
/*	... */ \
/*	30K	... */ \
/*	128K	... */ \
/*	... */ \
/* On output the 64-byte cache lines hold these data values: */ \
/*	0K	8K	16K	24K */ \
/*	2K	... */ \
/*	... */ \
/*	6K	... */ \
/*	32K	... */ \
/*	... */ \
\
/* Do 8 eight_reals_fft_2 macros */ \
/* Do 8 two_two_complex_fft macros */ \
/*	distance between fft data elements is 8K */ \
/*	do 4 macros each processing 32 data values */ \
\
start_timer(25); \
rdi = (uintptr_t)g->sincos2;/* Load sin/cos pointer */ \
/* 4 iterations of 2*clm */ \
for(unsigned int loopA = 4; loopA; loopA--) { \
	xtouch(xptr(rcx+(clm-1)*128));	/* Preload TLB */ \
	xloop1_clm##clm(b6b, clmblkdst, x4cl_eight_reals_fft_2(rsi, 64, 4*clmblkdst, clmblkdst8, rdi)); \
} \
\
/* Do 208 four_complex_fft macros */ \
/*	distance between fft data elements is 8K */ \
/*	do 52 macros each processing 32 data values */ \
\
rsi += -4*clmblkdst+2*clmblkdst8;/* Load source pointer */ \
rdi += XMM_SCD;		/* Next sine/cosine pointer */ \
/* 13 iterations of 4 iters of 2*clm */ \
for(unsigned int loopA = 13; loopA; loopA--) { \
	for(unsigned int loopB = 4; loopB; loopB--) { \
		xtouch(xptr(rcx+(clm-1)*128));	/* Preload TLB */ \
		xloop1_clm##clm(b7b, clmblkdst, x4cl_four_complex_fft(rsi, 64, 4*clmblkdst, clmblkdst8, rdi)); \
	} \
	rsi += -4*clmblkdst+2*clmblkdst8;/* Next source pointer */ \
	rdi += 2*XMM_SCD;		/* Next sine/cosine pointer */ \
} \
rsi += -14*2*clmblkdst8;	/* Restore source pointer */ \
pfing(rcx += -224*blkdst);	/* Restore prefetch pointer */ \
end_timer(25); \
\
/* Do FFT levels 9,10 */ \
/* Values 0K-7K is real data, 8K-15K is semi-real data, 16K-1791K is */ \
/* complex data. */ \
/* */ \
/* On input the 64-byte cache lines hold these data values: */ \
/*	0K	8K	16K	24K */ \
/*	2K	... */ \
/*	... */ \
/*	6K	... */ \
/*	32K	... */ \
/*	... */ \
/* On output the 64-byte cache lines hold these data values: */ \
/*	0K	2K	4K	6K */ \
/*	8K	... */ \
/*	... */ \
\
/* Do 1 eight_reals_fft_2 macros */ \
/* Do 1 nop_two_two_complex_fft_2 macros */ \
/* Do 2 four_complex_fft macros */ \
/*	distance between fft data elements is 2K */ \
/*	do 1 macros each processing 32 data values */ \
\
start_timer(26); \
rdi = (uintptr_t)g->sincos1;/* Load sin/cos pointer */ \
for(unsigned int loopA = 2*clm; loopA; loopA--) { \
	g4cl_eight_reals_fft_2(rsi, 64, clmblkdst, 2*clmblkdst, rcx, 64, blkdst, 2*blkdst); \
} \
\
/* Do 220 four_complex_fft macros */ \
/*	distance between fft data elements is 2K */ \
/*	do 55 macros each processing 32 data values */ \
\
rsi += -2*clm*64+4*clmblkdst;/* Next source pointer */ \
rcx += -2*clm*64+4*blkdst;	/* Next dest pointer */ \
rdi += XMM_SCD;		/* Next sine/cosine pointer */ \
/* eax = 28*256+1; */		/* 55 iterations of 2*clm */ \
for(unsigned int loopA = 55; loopA; ) { \
	do{ \
		for(unsigned int loopC = clm; loopC; loopC--) { \
			for(unsigned int loopD = 2; loopD; loopD--) { \
				g4cl_four_complex_fft(rsi, 64, clmblkdst, 2*clmblkdst, rcx, 64, blkdst, 2*blkdst); \
				L2prefetch128(u8ptr(rcx-64+rbp)); \
				L2prefetch128(u8ptr(rcx-64+blkdst+rbp)); \
				pfing(rbp += -64+2*blkdst); \
			} \
			pfing(rbp += 128-4*blkdst); \
		} \
		rsi += -2*clm*64+4*clmblkdst;/* Next source pointer */ \
		rcx += -2*clm*64+4*blkdst;	/* Next dest pointer */ \
		rdi += 2*XMM_SCD;		/* Next sine/cosine pointer */ \
	}while(--loopA & 1); \
	rsi += -8*clmblkdst+clmblkdst8; \
} \
end_timer(26); \
\
/* Work on next 4 sets of 896 values */ \
\
endlp:	xpass1_get_next_block(pass2, c0b, b0b); /* Get next block and jump */ \
/* done: */	xfft_footer_hg2; \
\
\
\
/* Perform a 2M-element FFT.  This is done in two passes.  Pass 1 */ \
/* does 10 levels, pass 2 uses common code to perform the last 11 levels. */ \
\

#define hg_pass1sclevels10pfa8(procname, xpass2name, pass2_levels, clm) \
 /* (fftname, pass2_levels, clmarg) */ \
\
/* set_FFT_constants(fftname, 0, 0, pass2_levels, clmarg); */ \
xfft_header_hg2(procname, xpass2name, pass2_levels, clm); \
\
/* Do the last 11 FFT levels */ \
/*	do 1 pass2_11_levels_real and 255 pass2_11_levels_complex */ \
\
/* pass2: */	start_timer(0); \
xpass2name(g); \
end_timer(0); \
if(g->ffttype == 1)		/* We're done if FFTing only */ \
	return; \
\
/* Do the inverse FFT */ \
\
pass1_inverse_fft_setup(c0b, procname); \
\
/* Do inverse FFT levels 9,10 */ \
/* On input the 64-byte cache lines hold these data values: */ \
/*	0K	2K	4K	6K */ \
/*	8K	10K	12K	14K */ \
/*	... */ \
/* On output the 64-byte cache lines hold these data values: */ \
/*	0K	4K	8K	12K */ \
/*	2K	... */ \
/*	16K	... */ \
/*	... */ \
\
/* Do 1 eight_reals_unfft_2 macros */ \
/* Do 1 nop_two_two_complex_unfft_2 macros */ \
/*	distance between fft data elements is 2K */ \
/*	do 1 macros each processing 16 data values */ \
\
c0b:	start_timer(16); \
rsi = (uintptr_t)g->data_addr;/* Get source address */ \
rcx = (uintptr_t)g->scratch_area;	/* Get address of scratch area */ \
pfing(rbp = (uintptr_t)g->data_prefetch);/* Set prefetch distance */ \
pfing(rbp -= rsi); \
for(unsigned int loopA = 2*clm; loopA; loopA--) { \
	g2cl_eight_reals_unfft_2(rsi, 64, blkdst, rcx, 64, clmblkdst); \
} \
\
/* Do 254 four_complex_unfft macros */ \
/*	distance between fft data elements is 2K */ \
/*	do 127 macros each processing 16 data values */ \
\
rsi += -2*clm*64+2*blkdst;	/* Next source pointer */ \
rcx += -2*clm*64+2*clmblkdst;/* Next dest pointer */ \
rdi = (uintptr_t)g->sincos1;/* Load sin/cos pointer */ \
/* al = 256/4+32; */ /* 127 iterations of 2*clm */ \
for(unsigned int loopA = 127; loopA; ) { \
	do{ \
		for(unsigned int loopC = 2*clm; loopC; loopC--) { \
			g2cl_four_complex_unfft(rsi, 64, blkdst, rcx, 64, clmblkdst); \
		} \
		rsi += -2*clm*64+2*blkdst;	/* Next source pointer */ \
		rcx += -2*clm*64+2*clmblkdst;/* Next dest pointer */ \
		rdi += XMM_SCD;		/* Next sine/cosine pointer */ \
	}while(--loopA & 3); \
	rcx += -8*clmblkdst+clmblkdst8; \
} \
rcx = rsi-128*2*blkdst;	/* Create prefetch pointer */ \
rsi = (uintptr_t)g->scratch_area;	/* Restore source pointer */ \
end_timer(16); \
\
/* Do inverse FFT levels 7,8 */ \
/* On input the 64-byte cache lines hold these data values: */ \
/*	0K	4K	8K	12K */ \
/*	2K	... */ \
/*	16K	... */ \
/*	... */ \
/* On output the 64-byte cache lines hold these data values: */ \
/*	0K	16K	32K	48K */ \
/*	2K	... */ \
/*	... */ \
/*	14K */ \
/*	64K	... */ \
/*	... */ \
\
/* Do 4 eight_reals_unfft_2 macros */ \
/* Do 4 nop_two_two_complex_unfft_2 macros */ \
/*	distance between fft data elements is 8K */ \
/*	do 2 macros each processing 32 data values */ \
\
start_timer(17); \
/* 2 iterations of 2*clm */ \
for(unsigned int loopA = 2; loopA; loopA--) { \
	xtouch(xptr(rcx+(clm-1)*128));	/* Preload TLB */ \
	xloop1_clm##clm(c6b, clmblkdst, x4cl_eight_reals_unfft_2(rsi, 64, 2*clmblkdst, 4*clmblkdst)); \
} \
\
/* Do 248 four_complex_unfft macros */ \
/*	distance between fft data elements is 8K */ \
/*	do 62 macros each processing 32 data values */ \
\
rsi += -2*clmblkdst+clmblkdst8;/* Next source pointer */ \
rdi = (uintptr_t)g->sincos2;/* Load sin/cos pointer */ \
/* 31 iterations of 2 iters of 2*clm */ \
for(unsigned int loopA = 31; loopA; loopA--) { \
	for(unsigned int loopB = 2; loopB; loopB--) { \
		xtouch(xptr(rcx+(clm-1)*128));	/* Preload TLB */ \
		xloop1_clm##clm(c7b, clmblkdst, x4cl_four_complex_unfft(rsi, 64, 2*clmblkdst, 4*clmblkdst, rdi)); \
	} \
	rsi += -2*clmblkdst+clmblkdst8;/* Next source pointer */ \
	rdi += XMM_SCD;		/* Next sine/cosine pointer */ \
} \
rsi += -32*clmblkdst8;	/* Restore source pointer */ \
end_timer(17); \
\
/* Do inverse FFT levels 6 */ \
/* On input the 64-byte cache lines hold these data values: */ \
/*	0K	16K	32K	48K */ \
/*	2K	... */ \
/*	... */ \
/*	14K */ \
/*	64K	... */ \
/*	... */ \
/* On output the 64-byte cache lines hold these data values: */ \
/*	0K	32K	64K	96K */ \
/*	2K	... */ \
/*	... */ \
/*	30K */ \
/*	128K	... */ \
/*	... */ \
/* */ \
/* Do 16 eight_reals_unfft_1 macros */ \
/*	distance between fft data elements is 32K */ \
/*	do 8 macros each processing 16 data values */ \
/* */ \
/* Do 240 four_complex_unfft macros */ \
/*	distance between fft data elements is 32K */ \
/*	do 120 macros each processing 16 data values */ \
\
/* Do inverse FFT levels 4,5 */ \
/* On input the 64-byte cache lines hold these data values: */ \
/*	0K	32K	64K	96K */ \
/*	2K	... */ \
/*	... */ \
/*	30K */ \
/*	128K	... */ \
/*	... */ \
/* On output the 64-byte cache lines hold these data values: */ \
/*	0K	128K	256K	384K */ \
/*	2K	... */ \
/*	... */ \
/*	126K */ \
/*	512K	... */ \
/*	... */ \
/* */ \
/* Do 32 eight_reals_unfft_2 macros */ \
/* Do 32 nop_two_two_complex_unfft_2 macros */ \
/*	distance between fft data elements is 64K */ \
/*	do 16 macros each processing 32 data values */ \
/* */ \
/* Do 192 four_complex_unfft macros */ \
/*	distance between fft data elements is 64K */ \
/*	do 48 macros each processing 32 data values */ \
\
/* Do inverse FFT levels 1,2,3 */ \
/* On input the 64-byte cache lines hold these data values: */ \
/*	0K	128K	256K	384K */ \
/*	2K	... */ \
/*	... */ \
/*	126K */ \
/*	512K	... */ \
/*	... */ \
/* On output the 64-byte cache lines hold these data values: */ \
/*	0K	512K	1024K	1536K */ \
/*	2K	... */ \
/*	... */ \
/*	510K */ \
/* */ \
/* Do 256 eight_reals_last_unfft macros */ \
/*	distance between fft data elements is 256K */ \
/*	do 64 macros each processing 32 data values */ \
\
start_timer(18); \
rdi = (uintptr_t)g->sincos3;/* Load sin/cos pointer */ \
rbx = (uintptr_t)g->sincos4;/* Load sin/cos pointer */ \
/* 8 iterations */ \
\
for(unsigned int loopA = 8; loopA; loopA--) { \
	xtouchlines(16);			/* Preload TLBs */ \
	\
	for(unsigned int loopB = clm; loopB; loopB--) { \
		for(unsigned int loopC = 2; loopC; loopC--) { \
			x2cl_eight_reals_unfft_1(rsi, 64, clmblkdst8); \
		} \
		\
		rsi += -2*64+2*clmblkdst8;	/* Next source pointer */ \
		/* 15 iterations */ \
		for(unsigned int loopC = 15; loopC; loopC--) { \
			for(unsigned int loopD = 2; loopD; loopD--) { \
				x2cl_two_complex_unfft(rsi, 64, clmblkdst8); \
			} \
			rsi += -2*64+2*clmblkdst8;	/* Next source pointer */ \
			rdi += XMM_SCD;		/* Next sine/cosine pointer */ \
		} \
		rdi += -15*XMM_SCD;	/* Restore sine/cosine pointer */ \
		\
		rsi += -16*2*clmblkdst8;	/* Restore source pointer */ \
		for(unsigned int loopC = 2; loopC; loopC--) { \
			xloop1sp(c2b, clmblkdst8, x4cl_eight_reals_unfft_2(rsi, 64, 2*clmblkdst8, 4*clmblkdst8)); \
		} \
		\
		rsi += -2*clmblkdst8+8*clmblkdst8;/* Next source pointer */ \
		/* 3 iterations */ \
		for(unsigned int loopC = 3; loopC; loopC--) { \
			for(unsigned int loopD = 2; loopD; loopD--) { \
				; \
				xloop1sp(c3b, clmblkdst8, x4cl_four_complex_unfft(rsi, 64, 2*clmblkdst8, 4*clmblkdst8, rbx)); \
			} \
			rsi += -2*clmblkdst8+8*clmblkdst8; \
			rbx += XMM_SCD;		/* Next sine/cosine pointer */ \
		} \
		rbx += -3*XMM_SCD;		/* Next sine/cosine pointer */ \
		\
		rsi += -4*8*clmblkdst8;	/* Restore source pointer */ \
		for(unsigned int loopC = 8; loopC; loopC--) { \
			xloop1sp(c1b, clmblkdst8, x4cl_eight_reals_last_unfft(rsi, 64, 8*clmblkdst8, 16*clmblkdst8)); \
		} \
		\
		EXEC248##clm(rsi += -8*clmblkdst8+2*64);	/* Next source pointer */ \
		EXEC248##clm(pfing(rcx += -16*blkdst+128));	/* Next prefetch pointer */ \
	} \
	EXEC248##clm(rsi += -clm*2*64+clmblkdst);/* Next source pointer */ \
	EXEC248##clm(pfing(rcx += -clm*128+16*blkdst));	/* Next prefetch pointer */ \
	EXEC1##clm(rsi += -8*clmblkdst8+clmblkdst);/* Next source pointer */ \
	\
} \
end_timer(18); \
\
/* Normalize these values */ \
\
pass1_normalize(1, yes_fft); \
\
/* Copy data from the scratch area back to the FFT data area */ \
\
copy_scratch_data_to_fft(256, clm); \
goto	endlp;			/* Skip FFT code if flag not set */ \
\
/* Fire up auxiliary threads to do pass 1 forward FFT */ \
\
/* pass1: */	pass1_forward_fft_setup(b0b); \
\
/* Copy data from the FFT data area back to the scratch area */ \
\
b0b:	copy_fft_data_to_scratch(256, clm); \
\
/* Do FFT levels 1,2,3 */ \
/* */ \
/* On input the 64-byte cache lines hold these data values: */ \
/*	0K	512K	1024K	1536K */ \
/*	2K	... */ \
/*	... */ \
/*	510K */ \
/* On output the 64-byte cache lines hold these data values: */ \
/*	0K	256K	512K	768K */ \
/*	2K	... */ \
/*	... */ \
/*	254K	... */ \
/*	1024K	... */ \
/*	... */ \
/* */ \
/* Do 256 eight_reals_first_fft macros (FFT levels 1,2,3) */ \
/*	distance between fft data elements is 256K */ \
/*	do 64 iterations each processing 16 data values */ \
\
/* Do FFT levels 4,5 */ \
/* Values 0K-255K is real data, 256K-511K is semi-real data, 512K-2047K is */ \
/* complex data. */ \
/* */ \
/* On input the 64-byte cache lines hold these data values: */ \
/*	0K	256K	512K	768K */ \
/*	2K	... */ \
/*	... */ \
/*	254K	... */ \
/*	1024K	... */ \
/*	... */ \
/* On output the 64-byte cache lines hold these data values: */ \
/*	0K	64K	128K	192K */ \
/*	2K	... */ \
/*	... */ \
/*	62K	... */ \
/*	256K	... */ \
/*	... */ \
/* */ \
/* Do 32 eight_reals_fft_2 macros */ \
/* Do 32 nop_two_two_complex_fft_2 macros */ \
/* Do 64 four_complex_fft macros */ \
/*	distance between fft data elements is 64K */ \
/*	do 32 macros each processing 32 data values */ \
/* */ \
/* Do 128 four_complex_fft macros */ \
/*	distance between fft data elements is 64K */ \
/*	do 32 macros each processing 32 data values */ \
\
/* Do FFT levels 6 */ \
/* Values 0K-63K is real data, 64K-127K is semi-real data, 128K-2047K is */ \
/* complex data. */ \
/* */ \
/* On input the 64-byte cache lines hold these data values: */ \
/*	0K	64K	128K	192K */ \
/*	2K	... */ \
/*	... */ \
/*	62K	... */ \
/*	256K	... */ \
/*	... */ \
/* On output the 64-byte cache lines hold these data values: */ \
/*	0K	32K	64K	96K */ \
/*	2K	... */ \
/*	... */ \
/*	30K	... */ \
/*	128K	... */ \
/*	... */ \
/* */ \
/* Do 16 eight_reals_fft_1 macros */ \
/* Do 16 two_two_complex_fft macros */ \
/*	distance between fft data elements is 32K */ \
/*	do 16 macros each processing 16 data values */ \
/* */ \
/* Do 224 two_two_complex_fft macros */ \
/*	distance between fft data elements is 32K */ \
/*	do 112 macros each processing 16 data values */ \
\
yes_fft:start_timer(22); \
rsi = (uintptr_t)g->scratch_area;	/* Get address of scratch area */ \
rcx = (uintptr_t)g->data_addr; /* Load source address for prefetching */ \
pfing(rbp = (uintptr_t)g->data_prefetch);/* Set prefetch distance */ \
pfing(rbp -= rcx); \
pfing(rcx += 128*blkdst);		/* Create prefetching pointer */ \
rdi = (uintptr_t)g->sincos4;/* Load sin/cos pointer */ \
rbx = (uintptr_t)g->sincos3;/* Load sin/cos pointer */ \
/* 16 iterations */ \
\
for(unsigned int loopA = 2; loopA; loopA--) { \
	for(unsigned int loopB = 8; loopB; loopB--) { \
		xtouchlines(4);			/* Preload TLBs */ \
		\
		for(unsigned int loopC = clm; loopC; loopC--) { \
			for(unsigned int loopD = 8; loopD; loopD--) { \
				for(unsigned int loopE = 2; loopE; loopE--) { \
					x2cl_eight_reals_fft(rsi, 64, 16*clmblkdst8); \
				} \
				rsi += -2*64+2*clmblkdst8;	/* Next source pointer */ \
			} \
			\
			rsi += -8*2*clmblkdst8;	/* Restore source pointer */ \
			for(unsigned int loopD = 2; loopD; loopD--) { \
				xloop1sp(b2b, 2*clmblkdst8, x4cl_eight_reals_fft_2(rsi, 64, 4*clmblkdst8, 8*clmblkdst8, rdi)); \
			} \
			\
			rsi += -2*2*clmblkdst8+16*clmblkdst8;/* Next source pointer */ \
			rdi += XMM_SCD;		/* Next sine/cosine pointer */ \
			for(unsigned int loopD = 2; loopD; loopD--) { \
				xloop1sp(b3b, 2*clmblkdst8, x4cl_four_complex_fft(rsi, 64, 4*clmblkdst8, 8*clmblkdst8, rdi)); \
			} \
			rdi += -XMM_SCD;		/* Restore sine/cosine pointer */ \
			\
			rsi += -2*2*clmblkdst8-16*clmblkdst8;/* Restore src pointer */ \
			for(unsigned int loopD = 2; loopD; loopD--) { \
				x2cl_eight_reals_fft_1(rsi, 64, 2*clmblkdst8); \
			} \
			\
			rsi += -2*64+4*clmblkdst8;	/* Load source pointer */ \
			rbx += XMM_SCD;		/* Next sine/cosine pointer */ \
			/* 7 iterations */ \
			for(unsigned int loopD = 7; loopD; loopD--) { \
				for(unsigned int loopE = 2; loopE; loopE--) { \
					x2cl_two_complex_fft(rsi, 64, 2*clmblkdst8, rbx); \
				} \
				rsi += -2*64+4*clmblkdst8;	/* Next source pointer */ \
				rbx += 2*XMM_SCD;		/* Next sine/cosine pointer */ \
			} \
			rbx += -15*XMM_SCD;	/* Restore sine/cosine pointer */ \
			\
				EXEC248##clm(rsi += -8*4*clmblkdst8+2*64);/* Next source pointer */ \
			EXEC248##clm(pfing(rcx += -4*blkdst+128));	/* Next prefetch pointer */ \
		} \
		EXEC248##clm(rsi += -clm*2*64+clmblkdst);/* Next source pointer */ \
		EXEC248##clm(pfing(rcx += -clm*128+4*blkdst));	/* Next prefetch pointer */ \
		EXEC1##clm(rsi += -8*4*clmblkdst8+clmblkdst);/* Next source pointer */ \
		\
	} \
	rsi += -8*clmblkdst+clmblkdst8;/* Next source pointer */ \
} \
rsi += -2*clmblkdst8;	/* Restore source pointer */ \
end_timer(22); \
\
/* Do FFT levels 7,8 */ \
/* Values 0K-31K is real data, 32K-63K is semi-real data, 64K-2047K is */ \
/* complex data. */ \
/* */ \
/* On input the 64-byte cache lines hold these data values: */ \
/*	0K	32K	64K	96K */ \
/*	2K	... */ \
/*	... */ \
/*	30K	... */ \
/*	128K	... */ \
/*	... */ \
/* On output the 64-byte cache lines hold these data values: */ \
/*	0K	8K	16K	24K */ \
/*	2K	... */ \
/*	... */ \
/*	6K	... */ \
/*	32K	... */ \
/*	... */ \
\
/* Do 8 eight_reals_fft_2 macros */ \
/* Do 8 two_two_complex_fft macros */ \
/*	distance between fft data elements is 8K */ \
/*	do 4 macros each processing 32 data values */ \
\
start_timer(25); \
rdi = (uintptr_t)g->sincos2;/* Load sin/cos pointer */ \
/* 4 iterations of 2*clm */ \
for(unsigned int loopA = 4; loopA; loopA--) { \
	xtouch(xptr(rcx+(clm-1)*128));	/* Preload TLB */ \
	xloop1_clm##clm(b6b, clmblkdst, x4cl_eight_reals_fft_2(rsi, 64, 4*clmblkdst, clmblkdst8, rdi)); \
} \
\
/* Do 240 four_complex_fft macros */ \
/*	distance between fft data elements is 8K */ \
/*	do 60 macros each processing 32 data values */ \
\
rsi += -4*clmblkdst+2*clmblkdst8;/* Load source pointer */ \
rdi += XMM_SCD;		/* Next sine/cosine pointer */ \
/* 15 iterations of 4 iters of 2*clm */ \
for(unsigned int loopA = 15; loopA; loopA--) { \
	for(unsigned int loopB = 4; loopB; loopB--) { \
		xtouch(xptr(rcx+(clm-1)*128));	/* Preload TLB */ \
		xloop1_clm##clm(b7b, clmblkdst, x4cl_four_complex_fft(rsi, 64, 4*clmblkdst, clmblkdst8, rdi)); \
	} \
	rsi += -4*clmblkdst+2*clmblkdst8;/* Next source pointer */ \
	rdi += 2*XMM_SCD;		/* Next sine/cosine pointer */ \
} \
rsi += -16*2*clmblkdst8;	/* Restore source pointer */ \
pfing(rcx += -256*blkdst);	/* Restore prefetch pointer */ \
end_timer(25); \
\
/* Do FFT levels 9,10 */ \
/* Values 0K-7K is real data, 8K-15K is semi-real data, 16K-2047K is */ \
/* complex data. */ \
/* */ \
/* On input the 64-byte cache lines hold these data values: */ \
/*	0K	8K	16K	24K */ \
/*	2K	... */ \
/*	... */ \
/*	6K	... */ \
/*	32K	... */ \
/*	... */ \
/* On output the 64-byte cache lines hold these data values: */ \
/*	0K	2K	4K	6K */ \
/*	8K	... */ \
/*	... */ \
\
/* Do 1 eight_reals_fft_2 macros */ \
/* Do 1 nop_two_two_complex_fft_2 macros */ \
/* Do 2 four_complex_fft macros */ \
/*	distance between fft data elements is 2K */ \
/*	do 1 macros each processing 32 data values */ \
\
start_timer(26); \
rdi = (uintptr_t)g->sincos1;/* Load sin/cos pointer */ \
for(unsigned int loopA = 2*clm; loopA; loopA--) { \
	g4cl_eight_reals_fft_2(rsi, 64, clmblkdst, 2*clmblkdst, rcx, 64, blkdst, 2*blkdst); \
} \
\
/* Do 252 four_complex_fft macros */ \
/*	distance between fft data elements is 2K */ \
/*	do 63 macros each processing 32 data values */ \
\
rsi += -2*clm*64+4*clmblkdst;/* Next source pointer */ \
rcx += -2*clm*64+4*blkdst;	/* Next dest pointer */ \
rdi += XMM_SCD;		/* Next sine/cosine pointer */ \
/* eax = 32*256+1; */ /* 63 iterations of 2*clm */ \
for(unsigned int loopA = 63; loopA; ) { \
	do{ \
		for(unsigned int loopC = clm; loopC; loopC--) { \
			for(unsigned int loopD = 2; loopD; loopD--) { \
				g4cl_four_complex_fft(rsi, 64, clmblkdst, 2*clmblkdst, rcx, 64, blkdst, 2*blkdst); \
				L2prefetch128(u8ptr(rcx-64+rbp)); \
				L2prefetch128(u8ptr(rcx-64+blkdst+rbp)); \
				pfing(rbp += -64+2*blkdst); \
			} \
			pfing(rbp += 128-4*blkdst); \
		} \
		rsi += -2*clm*64+4*clmblkdst;/* Next source pointer */ \
		rcx += -2*clm*64+4*blkdst;	/* Next dest pointer */ \
		rdi += 2*XMM_SCD;		/* Next sine/cosine pointer */ \
	}while(--loopA & 1); \
	rsi += -8*clmblkdst+clmblkdst8; \
} \
end_timer(26); \
\
/* Work on next 4 sets of 1024 values */ \
\
endlp:	xpass1_get_next_block(pass2, c0b, b0b); /* Get next block and jump */ \
/* done: */	xfft_footer_hg2; \
\
\
\
/* Perform a 2M-element negacyclic FFT.  This is done in two passes.  Pass 1 */ \
/* does 10 levels, pass 2 uses common code to perform the last 11 levels. */ \
\

#define hg_pass1sclevels10complex4(procname, xpass2name, pass2_levels, clm) \
 /* (fftname, pass2_levels, clmarg) */ \
\
/* set_FFT_constants(fftname, 0, 1, pass2_levels, clmarg); */ \
xfft_header_hg2(procname, xpass2name, pass2_levels, clm); \
\
/* Do the last 11 FFT levels */ \
/*	do 256 pass2_11_levels_complex */ \
\
/* pass2: */	start_timer(0); \
xpass2name(g); \
end_timer(0); \
if(g->ffttype == 1)		/* We're done if FFTing only */ \
	return; \
\
/* Do the inverse FFT */ \
\
pass1_inverse_fft_setup(c0b, procname); \
\
/* Do inverse FFT levels 9,10 */ \
/* On input the 64-byte cache lines hold these data values: */ \
/*	0K	2K	4K	6K */ \
/*	8K	10K	12K	14K */ \
/*	... */ \
/* On output the 64-byte cache lines hold these data values: */ \
/*	0K	4K	8K	12K */ \
/*	2K	... */ \
/*	16K	... */ \
/*	... */ \
\
/* Do 256 four_complex_unfft macros */ \
/*	distance between fft data elements is 2K */ \
/*	do 128 macros each processing 16 data values */ \
\
c0b:	start_timer(16); \
rsi = (uintptr_t)g->data_addr;/* Get source address */ \
rcx = (uintptr_t)g->scratch_area;	/* Get address of scratch area */ \
rdi = (uintptr_t)g->sincos1;/* Load sin/cos pointer */ \
/* 128 iterations */ \
pfing(rbp = (uintptr_t)g->data_prefetch);/* Set prefetch distance */ \
pfing(rbp -= rsi); \
for(unsigned int loopA = 32; loopA; loopA--) { \
	for(unsigned int loopB = 4; loopB; loopB--) { \
		for(unsigned int loopC = 2*clm; loopC; loopC--) { \
			g2cl_four_complex_unfft(rsi, 64, blkdst, rcx, 64, clmblkdst); \
		} \
		rsi += -2*clm*64+2*blkdst;	/* Next source pointer */ \
		rcx += -2*clm*64+2*clmblkdst;/* Next dest pointer */ \
		rdi += XMM_SCD;		/* Next sine/cosine pointer */ \
	} \
	rcx += -8*clmblkdst+clmblkdst8; \
} \
rcx = rsi-128*2*blkdst;	/* Create prefetch pointer */ \
rsi = (uintptr_t)g->scratch_area;	/* Restore source pointer */ \
end_timer(16); \
\
/* Do inverse FFT levels 7,8 */ \
/* On input the 64-byte cache lines hold these data values: */ \
/*	0K	4K	8K	12K */ \
/*	2K	... */ \
/*	16K	... */ \
/*	... */ \
/* On output the 64-byte cache lines hold these data values: */ \
/*	0K	16K	32K	48K */ \
/*	2K	... */ \
/*	... */ \
/*	14K */ \
/*	64K	... */ \
/*	... */ \
\
/* Do 256 four_complex_unfft macros */ \
/*	distance between fft data elements is 8K */ \
/*	do 64 macros each processing 32 data values */ \
\
start_timer(17); \
rdi = (uintptr_t)g->sincos2;/* Load sin/cos pointer */ \
/* 32 iterations of 2*2*clm */ \
for(unsigned int loopA = 32; loopA; loopA--) { \
	xtouch(xptr(rcx+(clm-1)*128));	/* Preload TLB */ \
	xloop2_clm##clm(c7b, clmblkdst, x4cl_four_complex_unfft(rsi, 64, 2*clmblkdst, 4*clmblkdst, rdi));		/* 2 * 2*clm iterations */ \
	rsi += -2*clmblkdst+clmblkdst8;/* Next source pointer */ \
	rdi += XMM_SCD;		/* Next sine/cosine pointer */ \
} \
rsi += -32*clmblkdst8;	/* Restore source pointer */ \
end_timer(17); \
\
/* Do inverse FFT levels 5,6 */ \
/* On input the 64-byte cache lines hold these data values: */ \
/*	0K	16K	32K	48K */ \
/*	2K	... */ \
/*	... */ \
/*	14K */ \
/*	64K	... */ \
/*	... */ \
/* On output the 64-byte cache lines hold these data values: */ \
/*	0K	64K	128K	192K */ \
/*	2K	... */ \
/*	... */ \
/*	62K */ \
/*	256K	... */ \
/*	... */ \
\
/* Do 256 four_complex_unfft macros */ \
/*	distance between fft data elements is 32K */ \
/*	do 64 macros each processing 32 data values */ \
\
start_timer(18); \
rdi = (uintptr_t)g->sincos3;/* Load sin/cos pointer */ \
/* 8 iterations of 4 iters of 2*2*clm */ \
for(unsigned int loopA = 8; loopA; loopA--) { \
	for(unsigned int loopB = 4; loopB; loopB--) { \
		xtouch(xptr(rcx+(clm-1)*128));	/* Preload TLB */ \
		xloop2_clm##clm(c5b, clmblkdst, x4cl_four_complex_unfft(rsi, 64, clmblkdst8, 2*clmblkdst8, rdi));		/* 2 * 2*clm iterations */ \
	} \
	rsi += -4*2*clmblkdst+4*clmblkdst8;/* Next source pointer */ \
	rdi += XMM_SCD;		/* Next sine/cosine pointer */ \
} \
rsi += -8*4*clmblkdst8;	/* Restore source pointer */ \
end_timer(18); \
\
/* Do inverse FFT levels 4 */ \
/* On input the 64-byte cache lines hold these data values: */ \
/*	0K	64K	128K	192K */ \
/*	2K	... */ \
/*	... */ \
/*	62K */ \
/*	256K	... */ \
/*	... */ \
/* On output the 64-byte cache lines hold these data values: */ \
/*	0K	128K	256K	384K */ \
/*	2K	... */ \
/*	... */ \
/*	126K */ \
/*	512K	... */ \
/*	... */ \
\
/* Do 256 two_two_complex_unfft macros */ \
/*	distance between fft data elements is 64K */ \
/*	do 128 macros each processing 16 data values */ \
\
start_timer(19); \
rdi = (uintptr_t)g->sincos4;/* Load sin/cos pointer */ \
/* 4 iters of 4 iters 2 of 4*2*clm */ \
for(unsigned int loopA = 4; loopA; loopA--) { \
	for(unsigned int loopB = 4; loopB; loopB--) { \
		for(unsigned int loopC = 2; loopC; loopC--) { \
			xtouch(xptr(rcx+(clm-1)*128));	/* Preload TLB */ \
			xloop4_clm##clm(c3b, clmblkdst, x2cl_two_complex_unfft(rsi, 64, 4*clmblkdst8));		/* 4* 2*clm iterations */ \
		} \
		rsi += -2*4*clmblkdst+clmblkdst8;/* Next source pointer */ \
	} \
	rsi += -4*clmblkdst8+8*clmblkdst8;/* Next source pointer */ \
	rdi += XMM_SCD;		/* Next sine/cosine pointer */ \
} \
rsi += -4*8*clmblkdst8;	/* Restore source pointer */ \
end_timer(19); \
\
/* Do inverse FFT levels 1,2,3 */ \
/* On input the 64-byte cache lines hold these data values: */ \
/*	0K	128K	256K	384K */ \
/*	2K	... */ \
/*	... */ \
/*	126K */ \
/*	512K	... */ \
/*	... */ \
/* On output the 64-byte cache lines hold these data values: */ \
/*	0K	512K	1024K	1536K */ \
/*	2K	... */ \
/*	... */ \
/*	510K */ \
\
/* Do 256 four_complex_unfft macros */ \
/*	distance between fft data elements is 256K */ \
/*	do 64 macros each processing 32 data values */ \
\
start_timer(20); \
rdi = (uintptr_t)g->plus1_premults;/* Address of premultiplier table */ \
/* 8 iters of 4 iters of 2*2*clm */ \
for(unsigned int loopA = 8; loopA; loopA--) { \
	for(unsigned int loopB = 4; loopB; loopB--) { \
		xtouch(xptr(rcx+(clm-1)*128));	/* Preload TLB */ \
		xloop2edi_clm##clm(c1b, clmblkdst, 128, x4cl_four_complex_last_unfft(rsi, 64, 8*clmblkdst8, 16*clmblkdst8, 64*128));	/* 2 * 2*clm iterations */ \
	} \
	rsi += -4*2*clmblkdst+clmblkdst8;/* Next source pointer */ \
} \
end_timer(20); \
\
/* Normalize these values */ \
\
pass1_normalize(1, yes_fft); \
\
/* Copy data from the scratch area back to the FFT data area */ \
\
copy_scratch_data_to_fft(256, clm); \
goto	endlp;			/* Skip FFT code if flag not set */ \
\
/* Fire up auxiliary threads to do pass 1 forward FFT */ \
\
/* pass1: */	pass1_forward_fft_setup(b0b); \
\
/* Copy data from the FFT data area back to the scratch area */ \
\
b0b:	copy_fft_data_to_scratch(256, clm); \
\
/* Do FFT levels 1,2,3 */ \
/* */ \
/* On input the 64-byte cache lines hold these data values: */ \
/*	0K	512K	1024K	1536K */ \
/*	2K	... */ \
/*	... */ \
/*	510K */ \
/* On output the 64-byte cache lines hold these data values: */ \
/*	0K	256K	512K	768K */ \
/*	2K	... */ \
/*	... */ \
/*	254K	... */ \
/*	1024K	... */ \
/*	... */ \
\
/* Do 256 four_complex_first_fft macros (FFT levels 1,2,3) */ \
/*	distance between fft data elements is 256K */ \
/*	do 128 iterations each processing 16 data values */ \
\
yes_fft:start_timer(22); \
rsi = (uintptr_t)g->scratch_area;	/* Get address of scratch area */ \
rdi = (uintptr_t)g->plus1_premults;/* Address of premultiplier table */ \
rcx = (uintptr_t)g->data_addr; /* Load source address for prefetching */ \
pfing(rbp = (uintptr_t)g->data_prefetch);/* Set prefetch distance */ \
pfing(rbp -= rcx); \
pfing(rcx += 128*blkdst);		/* Create prefetching pointer */ \
/* 16 iters of 2 iters of 4*2*clm */ \
for(unsigned int loopA = 16; loopA; loopA--) { \
	for(unsigned int loopB = 2; loopB; loopB--) { \
		xtouch(xptr(rcx+(clm-1)*128));	/* Preload TLB */ \
		xloop4edi_clm##clm(b1b, clmblkdst, 128, x2cl_four_complex_first_fft_scratch(rsi, 64, 16*clmblkdst8));	/* 4 * 2*clm iterations */ \
	} \
	rsi += -2*4*clmblkdst+clmblkdst8; \
} \
rsi += -16*clmblkdst8;	/* Restore source pointer */ \
end_timer(22); \
\
/* Do FFT level 4 */ \
/* */ \
/* On input the 64-byte cache lines hold these data values: */ \
/*	0K	256K	512K	768K */ \
/*	2K	... */ \
/*	... */ \
/*	254K	... */ \
/*	1024K	... */ \
/*	... */ \
/* On output the 64-byte cache lines hold these data values: */ \
/*	0K	128K	256K	384K */ \
/*	2K	... */ \
/*	... */ \
/*	126K	... */ \
/*	512K	... */ \
/*	... */ \
\
/* Do 256 two_two_complex_fft macros */ \
/*	distance between fft data elements is 64K */ \
/*	do 128 macros each processing 16 data values */ \
\
start_timer(23); \
rdi = (uintptr_t)g->sincos4;/* Load sin/cos pointer */ \
/* 2 iters of 8 iters of 2 of 4*2*clm */ \
for(unsigned int loopA = 2; loopA; loopA--) { \
	for(unsigned int loopB = 8; loopB; loopB--) { \
		for(unsigned int loopC = 2; loopC; loopC--) { \
			xtouch(xptr(rcx+(clm-1)*128));	/* Preload TLB */ \
			b3b:	; \
			xloop4_clm##clm(b3b, clmblkdst, x2cl_two_complex_fft(rsi, 64, 8*clmblkdst8, rdi));		/* 4 * 2*clm iterations */ \
		} \
		rsi += -2*4*clmblkdst+clmblkdst8; \
	} \
	rsi += -8*clmblkdst8+16*clmblkdst8;/* Next source pointer */ \
	rdi += 2*XMM_SCD;		/* Next sine/cosine pointer */ \
} \
rsi += -2*16*clmblkdst8;	/* Restore source pointer */ \
end_timer(23); \
\
/* Do FFT levels 5,6 */ \
/* */ \
/* On input the 64-byte cache lines hold these data values: */ \
/*	0K	128K	256K	384K */ \
/*	2K	... */ \
/*	... */ \
/*	126K	... */ \
/*	512K	... */ \
/*	... */ \
/* On output the 64-byte cache lines hold these data values: */ \
/*	0K	32K	64K	96K */ \
/*	2K	... */ \
/*	... */ \
/*	30K	... */ \
/*	128K	... */ \
/*	... */ \
\
/* Do 256 four_complex_fft macros */ \
/*	distance between fft data elements is 32K */ \
/*	do 64 macros each processing 32 data values */ \
\
start_timer(24); \
rdi = (uintptr_t)g->sincos3;/* Load sin/cos pointer */ \
/* 4 iters of 2 iters of 4 of 2*2*clm */ \
for(unsigned int loopA = 4; loopA; loopA--) { \
	for(unsigned int loopB = 2; loopB; loopB--) { \
		for(unsigned int loopC = 4; loopC; loopC--) { \
			xtouch(xptr(rcx+(clm-1)*128));	/* Preload TLB */ \
			xloop2_clm##clm(b5b, clmblkdst, x4cl_four_complex_fft(rsi, 64, 2*clmblkdst8, 4*clmblkdst8, rdi));		/* 2 * 2*clm iterations */ \
		} \
		rsi += -4*2*clmblkdst+clmblkdst8; \
	} \
	rsi += -2*clmblkdst8+8*clmblkdst8;/* Next source pointer */ \
	rdi += 2*XMM_SCD;		/* Next sine/cosine pointer */ \
} \
rsi += -4*8*clmblkdst8;	/* Restore source pointer */ \
end_timer(24); \
\
/* Do FFT levels 7,8 */ \
/* */ \
/* On input the 64-byte cache lines hold these data values: */ \
/*	0K	32K	64K	96K */ \
/*	2K	... */ \
/*	... */ \
/*	30K	... */ \
/*	128K	... */ \
/*	... */ \
/* On output the 64-byte cache lines hold these data values: */ \
/*	0K	8K	16K	24K */ \
/*	2K	... */ \
/*	... */ \
/*	6K	... */ \
/*	32K	... */ \
/*	... */ \
\
/* Do 256 four_complex_fft macros */ \
/*	distance between fft data elements is 8K */ \
/*	do 64 macros each processing 32 data values */ \
\
start_timer(25); \
rdi = (uintptr_t)g->sincos2;/* Load sin/cos pointer */ \
/* 16 iterations of 2 iters of 2*2*clm */ \
for(unsigned int loopA = 16; loopA; loopA--) { \
	for(unsigned int loopB = 2; loopB; loopB--) { \
		xtouch(xptr(rcx+(clm-1)*128));	/* Preload TLB */ \
		xloop2_clm##clm(b7b, clmblkdst, x4cl_four_complex_fft(rsi, 64, 4*clmblkdst, clmblkdst8, rdi));		/* 2 * 2*clm iterations */ \
	} \
	rsi += -2*2*clmblkdst+2*clmblkdst8;/* Next source pointer */ \
	rdi += 2*XMM_SCD;		/* Next sine/cosine pointer */ \
} \
rsi += -16*2*clmblkdst8;	/* Restore source pointer */ \
pfing(rcx += -256*blkdst);	/* Restore prefetch pointer */ \
end_timer(25); \
\
/* Do FFT levels 9,10 */ \
/* */ \
/* On input the 64-byte cache lines hold these data values: */ \
/*	0K	8K	16K	24K */ \
/*	2K	... */ \
/*	... */ \
/*	6K	... */ \
/*	32K	... */ \
/*	... */ \
/* On output the 64-byte cache lines hold these data values: */ \
/*	0K	2K	4K	6K */ \
/*	8K	... */ \
/*	... */ \
\
/* Do 256 four_complex_fft macros */ \
/*	distance between fft data elements is 2K */ \
/*	do 64 macros each processing 32 data values */ \
\
start_timer(26); \
rdi = (uintptr_t)g->sincos1;/* Load sin/cos pointer */ \
/* 32 iters of 2 of 2*clm */ \
for(unsigned int loopA = 32; loopA; loopA--) { \
	for(unsigned int loopB = 2; loopB; loopB--) { \
		for(unsigned int loopC = clm; loopC; loopC--) { \
			for(unsigned int loopD = 2; loopD; loopD--) { \
				g4cl_four_complex_fft(rsi, 64, clmblkdst, 2*clmblkdst, rcx, 64, blkdst, 2*blkdst); \
				L2prefetch128(u8ptr(rcx-64+rbp)); \
				L2prefetch128(u8ptr(rcx-64+blkdst+rbp)); \
				pfing(rbp += -64+2*blkdst); \
			} \
			pfing(rbp += +128-4*blkdst); \
		} \
	rsi += -2*clm*64+4*clmblkdst;/* Next source pointer */ \
		rcx += -2*clm*64+4*blkdst;	/* Next dest pointer */ \
		rdi += 2*XMM_SCD;		/* Next sine/cosine pointer */ \
	} \
	rsi += -8*clmblkdst+clmblkdst8; \
} \
end_timer(26); \
\
/* Work on next 4 sets of 1024 values */ \
\
endlp:	xpass1_get_next_block(pass2, c0b, b0b); /* Get next block and jump */ \
/* done: */	xfft_footer_hg2; \
\
\
\
/* Perform a 2560K-element FFT.  This is done in two passes.  Pass 1 */ \
/* does 11 levels, pass 2 uses common code to perform the last 11 levels. */ \
\

#define hg_pass1sclevels11pfa5(procname, xpass2name, pass2_levels, clm) \
 /* (fftname, pass2_levels, clmarg) */ \
\
/* set_FFT_constants(fftname, 0, 0, pass2_levels, clmarg); */ \
xfft_header_hg2(procname, xpass2name, pass2_levels, clm); \
\
/* Do the last 11 FFT levels */ \
/*	do 1 pass2_11_levels_real and 319 pass2_11_levels_complex */ \
\
/* pass2: */	start_timer(0); \
xpass2name(g); \
end_timer(0); \
if(g->ffttype == 1)		/* We're done if FFTing only */ \
	return; \
\
/* Do the inverse FFT */ \
\
pass1_inverse_fft_setup(c0b, procname); \
\
/* Do inverse FFT levels 10,11 */ \
/* On input the 64-byte cache lines hold these data values: */ \
/*	0K	2K	4K	6K */ \
/*	8K	10K	12K	14K */ \
/*	... */ \
/* On output the 64-byte cache lines hold these data values: */ \
/*	0K	4K	8K	12K */ \
/*	2K	... */ \
/*	16K	... */ \
/*	... */ \
\
/* Do 1 eight_reals_unfft_2 macros */ \
/* Do 1 nop_two_two_complex_unfft_2 macros */ \
/*	distance between fft data elements is 2K */ \
/*	do 1 macros each processing 16 data values */ \
\
c0b:	start_timer(16); \
rsi = (uintptr_t)g->data_addr;/* Get source address */ \
rcx = (uintptr_t)g->scratch_area;	/* Get address of scratch area */ \
pfing(rbp = (uintptr_t)g->data_prefetch);/* Set prefetch distance */ \
pfing(rbp -= rsi); \
for(unsigned int loopA = 2*clm; loopA; loopA--) { \
	g2cl_eight_reals_unfft_2(rsi, 64, blkdst, rcx, 64, clmblkdst); \
} \
\
/* Do 318 four_complex_unfft macros */ \
/*	distance between fft data elements is 2K */ \
/*	do 159 macros each processing 16 data values */ \
\
rsi += -2*clm*64+2*blkdst;	/* Next source pointer */ \
rcx += -2*clm*64+2*clmblkdst;/* Next dest pointer */ \
rdi = (uintptr_t)g->sincos1;/* Load sin/cos pointer */ \
/* al = 256/4+40; */ /* 159 iterations of 2*clm */ \
for(unsigned int loopA = 159; loopA; ) { \
	do{ \
		for(unsigned int loopC = 2*clm; loopC; loopC--) { \
			g2cl_four_complex_unfft(rsi, 64, blkdst, rcx, 64, clmblkdst); \
		} \
		rsi += -2*clm*64+2*blkdst;	/* Next source pointer */ \
		rcx += -2*clm*64+2*clmblkdst;/* Next dest pointer */ \
		rdi += XMM_SCD;		/* Next sine/cosine pointer */ \
	}while(--loopA & 3); \
	rcx += -8*clmblkdst+clmblkdst8; \
} \
rcx = rsi-160*2*blkdst;	/* Create prefetch pointer */ \
rsi = (uintptr_t)g->scratch_area;	/* Restore source pointer */ \
end_timer(16); \
\
/* Do inverse FFT levels 8,9 */ \
/* On input the 64-byte cache lines hold these data values: */ \
/*	0K	4K	8K	12K */ \
/*	2K	... */ \
/*	16K	... */ \
/*	... */ \
/* On output the 64-byte cache lines hold these data values: */ \
/*	0K	16K	32K	48K */ \
/*	2K	... */ \
/*	... */ \
/*	14K */ \
/*	64K	... */ \
/*	... */ \
\
/* Do 4 eight_reals_unfft_2 macros */ \
/* Do 4 nop_two_two_complex_unfft_2 macros */ \
/*	distance between fft data elements is 8K */ \
/*	do 2 macros each processing 32 data values */ \
\
start_timer(17); \
/* 2 iterations of 2*clm */ \
for(unsigned int loopA = 2; loopA; loopA--) { \
	for(unsigned int loopB = 2*clm; loopB; loopB--) { \
		x4cl_eight_reals_unfft_2(rsi, 64, 2*clmblkdst, 4*clmblkdst); \
	} \
	rsi += -2*clm*64+clmblkdst;/* Next source pointer */ \
} \
\
/* Do 312 four_complex_unfft macros */ \
/*	distance between fft data elements is 8K */ \
/*	do 78 macros each processing 32 data values */ \
\
rsi += -2*clmblkdst+clmblkdst8;/* Next source pointer */ \
rdi = (uintptr_t)g->sincos2;/* Load sin/cos pointer */ \
/* 39 iterations of 2*2*clm */ \
for(unsigned int loopA = 39; loopA; loopA--) { \
	xtouch(xptr(rcx+(clm-1)*128));	/* Preload TLB */ \
	xloop2_clm##clm(c7b, clmblkdst, x4cl_four_complex_unfft(rsi, 64, 2*clmblkdst, 4*clmblkdst, rdi)); \
	rsi += -2*clmblkdst+clmblkdst8;/* Next source pointer */ \
	rdi += XMM_SCD;		/* Next sine/cosine pointer */ \
} \
rsi += -40*clmblkdst8;	/* Restore source pointer */ \
end_timer(17); \
\
/* Do inverse FFT levels 6,7 */ \
/* On input the 64-byte cache lines hold these data values: */ \
/*	0K	16K	32K	48K */ \
/*	2K	... */ \
/*	... */ \
/*	14K */ \
/*	64K	... */ \
/*	... */ \
/* On output the 64-byte cache lines hold these data values: */ \
/*	0K	64K	128K	192K */ \
/*	2K	... */ \
/*	... */ \
/*	62K */ \
/*	256K	... */ \
/*	... */ \
\
/* Do 16 eight_reals_unfft_2 macros */ \
/* Do 16 nop_two_two_complex_unfft_2 macros */ \
/*	distance between fft data elements is 32K */ \
/*	do 8 macros each processing 32 data values */ \
\
start_timer(18); \
/* 8 iterations of 2*clm */ \
for(unsigned int loopA = 8; loopA; loopA--) { \
	for(unsigned int loopB = 2*clm; loopB; loopB--) { \
		x4cl_eight_reals_unfft_2(rsi, 64, clmblkdst8, 2*clmblkdst8); \
	} \
	rsi += -2*clm*64+clmblkdst;/* Next source pointer */ \
} \
\
/* Do 288 four_complex_unfft macros */ \
/*	distance between fft data elements is 8K */ \
/*	do 72 macros each processing 32 data values */ \
\
rsi += -8*clmblkdst+4*clmblkdst8;/* Next source pointer */ \
rdi = (uintptr_t)g->sincos3;/* Load sin/cos pointer */ \
/* 9 iterations of 4 iters of 2*2*clm */ \
for(unsigned int loopA = 9; loopA; loopA--) { \
	for(unsigned int loopB = 4; loopB; loopB--) { \
		xtouch(xptr(rcx+(clm-1)*128));	/* Preload TLB */ \
		xloop2_clm##clm(c5b, clmblkdst, x4cl_four_complex_unfft(rsi, 64, clmblkdst8, 2*clmblkdst8, rdi)); \
	} \
	rsi += -4*2*clmblkdst+4*clmblkdst8;/* Next source pointer */ \
	rdi += XMM_SCD;		/* Next sine/cosine pointer */ \
} \
rsi += -10*4*clmblkdst8;	/* Restore source pointer */ \
end_timer(18); \
\
/* Do inverse FFT levels 4,5 */ \
/* On input the 64-byte cache lines hold these data values: */ \
/*	0K	64K	128K	192K */ \
/*	2K	... */ \
/*	... */ \
/*	62K */ \
/*	256K	... */ \
/*	... */ \
/* On output the 64-byte cache lines hold these data values: */ \
/*	0K	128K	256K	384K */ \
/*	2K	... */ \
/*	... */ \
/*	126K */ \
/*	512K	768K	1024K	1280K */ \
/*	... */ \
/* */ \
/* Do 64 eight_reals_unfft_2 macros */ \
/*	distance between fft data elements is 128K */ \
/*	do 32 macros each processing 16 data values */ \
/* */ \
/* Do 256 four_complex_unfft macros */ \
/*	distance between fft data elements is 128K */ \
/*	do 64 macros each processing 32 data values */ \
\
/* Do inverse FFT levels 1,2,3 */ \
/* On input the 64-byte cache lines hold these data values: */ \
/*	0K	128K	256K	384K */ \
/*	2K	... */ \
/*	... */ \
/*	126K */ \
/*	512K	768K	1024K	1280K */ \
/*	... */ \
/* On output the 64-byte cache lines hold these data values: */ \
/*	0K	128K	1280K	1408K */ \
/*	2K	... */ \
/*	... */ \
/*	126K */ \
/*	256K	768K	1536K	2048K */ \
/*	258K */ \
/*	... */ \
/*	766K */ \
/* */ \
/* Do 512 five_reals_last_unfft macros */ \
/*	distance between fft data elements is 512K */ \
/*	do 64 macros each processing 40 data values */ \
\
start_timer(19); \
rdi = (uintptr_t)g->sincos4;/* Load sin/cos pointer */ \
/* 32 iterations */ \
\
for(unsigned int loopA = 4; loopA; loopA--) { \
	for(unsigned int loopB = 8; loopB; loopB--) { \
		xtouchlines(3);			/* Preload TLBs */ \
		\
		for(unsigned int loopC = clm; loopC; loopC--) { \
			for(unsigned int loopD = 2; loopD; loopD--) { \
				x2cl_half_eight_reals_unfft_2(rsi, 64, 4*clmblkdst8); \
			} \
			\
			rsi += -2*64+8*clmblkdst8;	/* Next source pointer */ \
			for(unsigned int loopD = 2; loopD; loopD--) { \
				for(unsigned int loopE = 2; loopE; loopE--) { \
					x4cl_four_complex_unfft(rsi, 64, 4*clmblkdst8, 8*clmblkdst8, rdi); \
				} \
				rsi += -2*64+16*clmblkdst8;/* Next source pointer */ \
				rdi += XMM_SCD;		/* Next sine/cosine pointer */ \
			} \
			L2prefetch128(u8ptr(rcx+rbp)); \
			pfing(rcx += blkdst);		/* Next prefetch pointer */ \
			rdi += -2*XMM_SCD;		/* Restore sine/cosine pointer */ \
			\
			rsi += -2*16*clmblkdst8-8*clmblkdst8;/* Restore src pointer */ \
			xfive_reals_unfft_preload;	/* Preload consts for five real macros */ \
			for(unsigned int loopD = 2; loopD; loopD--) { \
				xloop1sp(c1b, 4*clmblkdst8, x5cl_five_reals_last_unfft(rsi, 64, 8*clmblkdst8)); \
			} \
			\
				EXEC248##clm(rsi += -2*4*clmblkdst8+2*64);/* Next source pointer */ \
			EXEC248##clm(pfing(rcx += -3*blkdst+128));	/* Next prefetch pointer */ \
		} \
		EXEC248##clm(rsi += -clm*2*64+clmblkdst);/* Next source pointer */ \
		EXEC248##clm(pfing(rcx += -clm*128+3*blkdst));	/* Next prefetch pointer */ \
		EXEC1##clm(rsi += -2*4*clmblkdst8+clmblkdst;)/* Next source pointer */ \
		\
	} \
	rsi += -8*clmblkdst+clmblkdst8;/* Next source pointer */ \
} \
end_timer(19); \
\
/* Normalize these values */ \
\
pass1_normalize(1, yes_fft); \
\
/* Copy data from the scratch area back to the FFT data area */ \
\
copy_scratch_data_to_fft(320, clm); \
goto	endlp;			/* Skip FFT code if flag not set */ \
\
/* Fire up auxiliary threads to do pass 1 forward FFT */ \
\
/* pass1: */	pass1_forward_fft_setup(b0b); \
\
/* Copy data from the FFT data area back to the scratch area */ \
\
b0b:	copy_fft_data_to_scratch(320, clm); \
\
/* Do FFT levels 1,2,3 */ \
/* */ \
/* On input the 64-byte cache lines hold these data values: */ \
/*	0K	128K	1280K	1408K */ \
/*	2K	... */ \
/*	... */ \
/*	126K */ \
/*	256K	768K	1536K	2048K */ \
/*	258K */ \
/*	... */ \
/*	766K */ \
/* On output the 64-byte cache lines hold these data values: */ \
/*	0K	128K	256K	384K */ \
/*	2K	... */ \
/*	... */ \
/*	126K */ \
/*	512K	1024K	1536K	2048K */ \
/*	514K */ \
/*	... */ \
/*	1022K */ \
/* */ \
/* Do 512 five_reals_first_fft macros (FFT levels 1,2,3) */ \
/*	distance between fft data elements is 512K */ \
/*	do 64 iterations each processing 40 data values */ \
\
/* Do FFT levels 4,5 */ \
/* Values 0K-511K is real data, 512K-2559K is complex data. */ \
/* */ \
/* On input the 64-byte cache lines hold these data values: */ \
/*	0K	128K	256K	384K */ \
/*	2K	... */ \
/*	... */ \
/*	126K */ \
/*	512K	1024K	1536K	2048K */ \
/*	514K */ \
/*	... */ \
/*	1022K */ \
/* On output the 64-byte cache lines hold these data values: */ \
/*	0K	128K	256K	384K */ \
/*	2K	... */ \
/*	... */ \
/*	126K	... */ \
/*	512K	... */ \
/*	... */ \
/* */ \
/* Do 64 eight_reals_fft_2 macros */ \
/*	distance between fft data elements is 128K */ \
/*	do 64 macros each processing 8 data values */ \
/* */ \
/* Do 256 four_complex_fft macros */ \
/*	distance between fft data elements is 128K */ \
/*	do 64 macros each processing 32 data values */ \
\
/* Do FFT levels 6,7 */ \
/* Values 0K-127K is real data, 128K-255K is semi-real data, 256K-2559K is */ \
/* complex data. */ \
/* */ \
/* On input the 64-byte cache lines hold these data values: */ \
/*	0K	128K	256K	384K */ \
/*	2K	... */ \
/*	... */ \
/*	126K	... */ \
/*	512K	... */ \
/*	... */ \
/* On output the 64-byte cache lines hold these data values: */ \
/*	0K	32K	64K	96K */ \
/*	2K	... */ \
/*	... */ \
/*	30K	... */ \
/*	128K	... */ \
/*	... */ \
/* */ \
/* Do 32 eight_reals_fft_2 macros */ \
/* Do 32 two_two_complex_fft macros */ \
/*	distance between fft data elements is 32K */ \
/*	do 16 macros each processing 32 data values */ \
/* */ \
/* Do 256 four_complex_fft macros */ \
/*	distance between fft data elements is 8K */ \
/*	do 64 macros each processing 32 data values */ \
\
yes_fft:start_timer(22); \
rsi = (uintptr_t)g->scratch_area;	/* Get address of scratch area */ \
rcx = (uintptr_t)g->data_addr; /* Load source address for prefetching */ \
pfing(rbp = (uintptr_t)g->data_prefetch);/* Set prefetch distance */ \
pfing(rbp -= rcx); \
pfing(rcx += 170*blkdst);		/* Create prefetching pointer */ \
rdi = (uintptr_t)g->sincos4;/* Load sin/cos pointer */ \
rbx = (uintptr_t)g->sincos3;/* Load sin/cos pointer */ \
/* 16 iterations */ \
\
for(unsigned int loopA = 2; loopA; loopA--) { \
	for(unsigned int loopB = 8; loopB; loopB--) { \
		for(unsigned int loopC = clm; loopC; loopC--) { \
			xtouchlines(7);			/* Preload TLBs */ \
			\
			xfive_reals_fft_preload;		/* Preload consts for five real macros */ \
			for(unsigned int loopD = 4; loopD; loopD--) { \
				xloop1sp(b1b, 2*clmblkdst8, x5cl_five_reals_fft(rsi, 64, 8*clmblkdst8)); \
			} \
			\
			rsi += -4*2*clmblkdst8;	/* Restore source pointer */ \
			for(unsigned int loopD = 4; loopD; loopD--) { \
				for(unsigned int loopE = 2; loopE; loopE--) { \
					x1cl_half_eight_reals_fft_2(rsi, 64); \
				} \
				rsi += -2*64+2*clmblkdst8; \
			} \
			L2prefetch128(u8ptr(rcx+rbp)); \
			pfing(rcx += blkdst);		/* Next prefetch pointer */ \
			\
			/*	bump	rsi, -4*2*clmblkdst8+8*clmblkdst8 */ \
			for(unsigned int loopD = 4; loopD; loopD--) { \
				for(unsigned int loopE = 2; loopE; loopE--) { \
					x4cl_four_complex_fft(rsi, 64, 8*clmblkdst8, 16*clmblkdst8, rdi); \
				} \
				rsi += -2*64+2*clmblkdst8; \
			} \
			L2prefetch128(u8ptr(rcx+rbp)); \
			pfing(rcx += blkdst);		/* Next prefetch pointer */ \
			\
			rsi += -4*2*clmblkdst8-8*clmblkdst8;/*Restore src pointer */ \
			for(unsigned int loopD = 2; loopD; loopD--) { \
				x4cl_eight_reals_fft_2(rsi, 64, 2*clmblkdst8, 4*clmblkdst8, rbx); \
			} \
			L2prefetch128(u8ptr(rcx+rbp)); \
			pfing(rcx += blkdst);		/* Next prefetch pointer */ \
			\
			rsi += -2*64+8*clmblkdst8;	/* Next source pointer */ \
			rbx += XMM_SCD;		/* Next sine/cosine pointer */ \
			for(unsigned int loopD = 4; loopD; loopD--) { \
				for(unsigned int loopE = 2; loopE; loopE--) { \
					x4cl_four_complex_fft(rsi, 64, 2*clmblkdst8, 4*clmblkdst8, rbx); \
				} \
				rsi += -2*64+8*clmblkdst8; \
				rbx += 2*XMM_SCD;		/* Next sine/cosine pointer */ \
			} \
			rbx += -9*XMM_SCD;		/* Restore sine/cosine pointer */ \
			\
				EXEC248##clm(rsi += -5*8*clmblkdst8+2*64);/* Next source pointer */ \
			EXEC248##clm(pfing(rcx += -7*blkdst+128));	/* Next prefetch pointer */ \
		} \
		EXEC248##clm(rsi += -clm*2*64+clmblkdst);/* Next source pointer */ \
		EXEC248##clm(pfing(rcx += -clm*128+7*blkdst));	/* Next prefetch pointer */ \
		EXEC1##clm(rsi += -5*8*clmblkdst8+clmblkdst);/* Next source pointer */ \
		\
	} \
	rsi += -8*clmblkdst+clmblkdst8;/* Next source pointer */ \
} \
rsi += -2*clmblkdst8;	/* Restore source pointer */ \
end_timer(22); \
\
/* Do FFT levels 8,9 */ \
/* Values 0K-31K is real data, 32K-63K is semi-real data, 64K-2559K is */ \
/* complex data. */ \
/* */ \
/* On input the 64-byte cache lines hold these data values: */ \
/*	0K	32K	64K	96K */ \
/*	2K	... */ \
/*	... */ \
/*	30K	... */ \
/*	128K	... */ \
/*	... */ \
/* On output the 64-byte cache lines hold these data values: */ \
/*	0K	8K	16K	24K */ \
/*	2K	... */ \
/*	... */ \
/*	6K	... */ \
/*	32K	... */ \
/*	... */ \
\
/* Do 8 eight_reals_fft_2 macros */ \
/* Do 8 two_two_complex_fft macros */ \
/*	distance between fft data elements is 8K */ \
/*	do 4 macros each processing 32 data values */ \
\
start_timer(25); \
rdi = (uintptr_t)g->sincos2;/* Load sin/cos pointer */ \
/* 4 iterations of 2*clm */ \
for(unsigned int loopA = 4; loopA; loopA--) { \
	for(unsigned int loopB = 2*clm; loopB; loopB--) { \
		x4cl_eight_reals_fft_2(rsi, 64, 4*clmblkdst, clmblkdst8, rdi); \
	} \
	rsi += -2*clm*64+clmblkdst;/* Next source pointer */ \
} \
\
/* Do 304 four_complex_fft macros */ \
/*	distance between fft data elements is 8K */ \
/*	do 76 macros each processing 32 data values */ \
\
rsi += -4*clmblkdst+2*clmblkdst8;/* Load source pointer */ \
rdi += XMM_SCD;		/* Next sine/cosine pointer */ \
/* 19 iterations of 2 iters of 2*2*clm */ \
for(unsigned int loopA = 19; loopA; loopA--) { \
	for(unsigned int loopB = 2; loopB; loopB--) { \
		xtouch(xptr(rcx+(clm-1)*128));	/* Preload TLB */ \
		xloop2_clm##clm(b7b, clmblkdst, x4cl_four_complex_fft(rsi, 64, 4*clmblkdst, clmblkdst8, rdi)); \
	} \
	rsi += -2*2*clmblkdst+2*clmblkdst8;/* Next source pointer */ \
	rdi += 2*XMM_SCD;		/* Next sine/cosine pointer */ \
} \
rsi += -20*2*clmblkdst8;	/* Restore source pointer */ \
pfing(rcx += -320*blkdst);	/* Restore prefetch pointer */ \
end_timer(25); \
\
/* Do FFT levels 10,11 */ \
/* Values 0K-7K is real data, 8K-15K is semi-real data, 16K-2559K is */ \
/* complex data. */ \
/* */ \
/* On input the 64-byte cache lines hold these data values: */ \
/*	0K	8K	16K	24K */ \
/*	2K	... */ \
/*	... */ \
/*	6K	... */ \
/*	32K	... */ \
/*	... */ \
/* On output the 64-byte cache lines hold these data values: */ \
/*	0K	2K	4K	6K */ \
/*	8K	... */ \
/*	... */ \
\
/* Do 1 eight_reals_fft_2 macros */ \
/* Do 1 nop_two_two_complex_fft_2 macros */ \
/* Do 2 four_complex_fft macros */ \
/*	distance between fft data elements is 2K */ \
/*	do 1 macros each processing 32 data values */ \
\
start_timer(26); \
rdi = (uintptr_t)g->sincos1;/* Load sin/cos pointer */ \
for(unsigned int loopA = 2*clm; loopA; loopA--) { \
	g4cl_eight_reals_fft_2(rsi, 64, clmblkdst, 2*clmblkdst, rcx, 64, blkdst, 2*blkdst); \
} \
\
/* Do 316 four_complex_fft macros */ \
/*	distance between fft data elements is 2K */ \
/*	do 79 macros each processing 32 data values */ \
\
rsi += -2*clm*64+4*clmblkdst;/* Next source pointer */ \
rcx += -2*clm*64+4*blkdst;	/* Next dest pointer */ \
rdi += XMM_SCD;		/* Next sine/cosine pointer */ \
/* eax = 40*256+1; */ /* 79 iterations of 2*clm */ \
for(unsigned int loopA = 79; loopA; ) { \
	do{ \
		for(unsigned int loopC = clm; loopC; loopC--) { \
			for(unsigned int loopD = 2; loopD; loopD--) { \
				g4cl_four_complex_fft(rsi, 64, clmblkdst, 2*clmblkdst, rcx, 64, blkdst, 2*blkdst); \
				L2prefetch128(u8ptr(rcx-64+rbp)); \
				L2prefetch128(u8ptr(rcx-64+blkdst+rbp)); \
				pfing(rbp += -64+2*blkdst); \
			} \
			pfing(rbp += 128-4*blkdst); \
		} \
		rsi += -2*clm*64+4*clmblkdst;/* Next source pointer */ \
		rcx += -2*clm*64+4*blkdst;	/* Next source pointer */ \
		rdi += 2*XMM_SCD;		/* Next sine/cosine pointer */ \
	}while(--loopA & 1); \
	rsi += -8*clmblkdst+clmblkdst8; \
} \
end_timer(26); \
\
/* Work on next 4 sets of 1280 values */ \
\
endlp:	xpass1_get_next_block(pass2, c0b, b0b); /* Get next block and jump */ \
/* done: */	xfft_footer_hg2; \
\
\
\
/* Perform a 3M-element FFT.  This is done in two passes.  Pass 1 */ \
/* does 11 levels, pass 2 uses common code to perform the last 11 levels. */ \
\

#define hg_pass1sclevels11pfa6(procname, xpass2name, pass2_levels, clm) \
 /* (fftname, pass2_levels, clmarg) */ \
\
/* set_FFT_constants(fftname, 0, 0, pass2_levels, clmarg); */ \
xfft_header_hg2(procname, xpass2name, pass2_levels, clm); \
\
/* Do the last 11 FFT levels */ \
/*	do 1 pass2_11_levels_real and 383 pass2_11_levels_complex */ \
\
/* pass2: */	start_timer(0); \
xpass2name(g); \
end_timer(0); \
if(g->ffttype == 1)		/* We're done if FFTing only */ \
	return; \
\
/* Do the inverse FFT */ \
\
pass1_inverse_fft_setup(c0b, procname); \
\
/* Do inverse FFT levels 10,11 */ \
/* On input the 64-byte cache lines hold these data values: */ \
/*	0K	2K	4K	6K */ \
/*	8K	10K	12K	14K */ \
/*	... */ \
/* On output the 64-byte cache lines hold these data values: */ \
/*	0K	4K	8K	12K */ \
/*	2K	... */ \
/*	16K	... */ \
/*	... */ \
\
/* Do 1 eight_reals_unfft_2 macros */ \
/* Do 1 nop_two_two_complex_unfft_2 macros */ \
/*	distance between fft data elements is 2K */ \
/*	do 1 macros each processing 16 data values */ \
\
c0b:	start_timer(16); \
rsi = (uintptr_t)g->data_addr;/* Get source address */ \
rcx = (uintptr_t)g->scratch_area;	/* Get address of scratch area */ \
pfing(rbp = (uintptr_t)g->data_prefetch);/* Set prefetch distance */ \
pfing(rbp -= rsi); \
for(unsigned int loopA = 2*clm; loopA; loopA--) { \
	g2cl_eight_reals_unfft_2(rsi, 64, blkdst, rcx, 64, clmblkdst); \
} \
\
/* Do 382 four_complex_unfft macros */ \
/*	distance between fft data elements is 2K */ \
/*	do 191 macros each processing 16 data values */ \
\
rsi += -2*clm*64+2*blkdst;	/* Next source pointer */ \
rcx += -2*clm*64+2*clmblkdst;/* Next dest pointer */ \
rdi = (uintptr_t)g->sincos1;/* Load sin/cos pointer */ \
/* al = 256/4+48; */ /* 191 iterations of 2*clm */ \
for(unsigned int loopA = 256; loopA; ) { \
	do{ \
		for(unsigned int loopC = 2*clm; loopC; loopC--) { \
			g2cl_four_complex_unfft(rsi, 64, blkdst, rcx, 64, clmblkdst); \
		} \
		rsi += -2*clm*64+2*blkdst;	/* Next source pointer */ \
		rcx += -2*clm*64+2*clmblkdst;/* Next dest pointer */ \
		rdi += XMM_SCD;		/* Next sine/cosine pointer */ \
	}while(--loopA & 3); \
	rcx += -8*clmblkdst+clmblkdst8; \
} \
rcx = rsi-192*2*blkdst;	/* Create prefetch pointer */ \
rsi = (uintptr_t)g->scratch_area;	/* Restore source pointer */ \
end_timer(16); \
\
/* Do inverse FFT levels 8,9 */ \
/* On input the 64-byte cache lines hold these data values: */ \
/*	0K	4K	8K	12K */ \
/*	2K	... */ \
/*	16K	... */ \
/*	... */ \
/* On output the 64-byte cache lines hold these data values: */ \
/*	0K	16K	32K	48K */ \
/*	2K	... */ \
/*	... */ \
/*	14K */ \
/*	64K	... */ \
/*	... */ \
\
/* Do 4 eight_reals_unfft_2 macros */ \
/* Do 4 nop_two_two_complex_unfft_2 macros */ \
/*	distance between fft data elements is 8K */ \
/*	do 2 macros each processing 32 data values */ \
\
start_timer(17); \
/* 2 iterations  */ \
for(unsigned int loopA = 2; loopA; loopA--) { \
	for(unsigned int loopB = 2*clm; loopB; loopB--) { \
		x4cl_eight_reals_unfft_2(rsi, 64, 2*clmblkdst, 4*clmblkdst); \
	} \
	rsi += -2*clm*64+clmblkdst;/* Next source pointer */ \
} \
\
/* Do 376 four_complex_unfft macros */ \
/*	distance between fft data elements is 8K */ \
/*	do 94 macros each processing 32 data values */ \
\
rsi += -2*clmblkdst+clmblkdst8;/* Next source pointer */ \
rdi = (uintptr_t)g->sincos2;/* Load sin/cos pointer */ \
/* 47 iterations of 2*2*clm */ \
for(unsigned int loopA = 47; loopA; loopA--) { \
	xtouch(xptr(rcx+(clm-1)*128));	/* Preload TLB */ \
	xloop2_clm##clm(c7b, clmblkdst, x4cl_four_complex_unfft(rsi, 64, 2*clmblkdst, 4*clmblkdst, rdi)); \
	rsi += -2*clmblkdst+clmblkdst8;/* Next source pointer */ \
	rdi += XMM_SCD;		/* Next sine/cosine pointer */ \
} \
rsi += -48*clmblkdst8;	/* Restore source pointer */ \
end_timer(17); \
\
/* Do inverse FFT levels 6,7 */ \
/* On input the 64-byte cache lines hold these data values: */ \
/*	0K	16K	32K	48K */ \
/*	2K	... */ \
/*	... */ \
/*	14K */ \
/*	64K	... */ \
/*	... */ \
/* On output the 64-byte cache lines hold these data values: */ \
/*	0K	64K	128K	192K */ \
/*	2K	... */ \
/*	... */ \
/*	62K */ \
/*	256K	... */ \
/*	... */ \
\
/* Do 16 eight_reals_unfft_2 macros */ \
/* Do 16 nop_two_two_complex_unfft_2 macros */ \
/*	distance between fft data elements is 32K */ \
/*	do 8 macros each processing 32 data values */ \
\
start_timer(18); \
/* 8 iterations of 2*clm */ \
for(unsigned int loopA = 8; loopA; loopA--) { \
	for(unsigned int loopB = 2*clm; loopB; loopB--) { \
		x4cl_eight_reals_unfft_2(rsi, 64, clmblkdst8, 2*clmblkdst8); \
	} \
	rsi += -2*clm*64+clmblkdst;/* Next source pointer */ \
} \
\
/* Do 352 four_complex_unfft macros */ \
/*	distance between fft data elements is 32K */ \
/*	do 88 macros each processing 32 data values */ \
\
rsi += -8*clmblkdst+4*clmblkdst8;/* Next source pointer */ \
rdi = (uintptr_t)g->sincos3;/* Load sin/cos pointer */ \
/* 11 iterations of 4 iters of 2*2*clm */ \
for(unsigned int loopA = 11; loopA; loopA--) { \
	for(unsigned int loopB = 4; loopB; loopB--) { \
		xtouch(xptr(rcx+(clm-1)*128));	/* Preload TLB */ \
		xloop2_clm##clm(c5b, clmblkdst, x4cl_four_complex_unfft(rsi, 64, clmblkdst8, 2*clmblkdst8, rdi)); \
	} \
	rsi += -4*2*clmblkdst+4*clmblkdst8;/* Next source pointer */ \
	rdi += XMM_SCD;		/* Next sine/cosine pointer */ \
} \
rsi += -12*4*clmblkdst8;	/* Restore source pointer */ \
end_timer(18); \
\
/* Do inverse FFT levels 4,5 */ \
/* On input the 64-byte cache lines hold these data values: */ \
/*	0K	64K	128K	192K */ \
/*	2K	... */ \
/*	... */ \
/*	62K */ \
/*	256K	... */ \
/*	... */ \
/* On output the 64-byte cache lines hold these data values: */ \
/*	0K	256K	512K	768K */ \
/*	2K	... */ \
/*	... */ \
/*	254K */ \
/*	1024K	... */ \
/*	... */ \
/* */ \
/* Do 64 eight_reals_unfft_2 macros */ \
/* Do 64 nop_two_two_complex_unfft_2 macros */ \
/*	distance between fft data elements is 128K */ \
/*	do 32 macros each processing 32 data values */ \
/* */ \
/* Do 256 four_complex_unfft macros */ \
/*	distance between fft data elements is 128K */ \
/*	do 64 macros each processing 32 data values */ \
\
/* Do inverse FFT levels 1,2,3 */ \
/* On input the 64-byte cache lines hold these data values: */ \
/*	0K	256K	512K	768K */ \
/*	2K	... */ \
/*	... */ \
/*	254K */ \
/*	1024K	1536K	2048K	2560K */ \
/*	1026K */ \
/*	... */ \
/*	1534K */ \
/* On output the 64-byte cache lines hold these data values: */ \
/*	0K	768K	1536K	2304K */ \
/*	... */ \
/*	766K */ \
/* */ \
/* Do 512 six_reals_last_unfft macros */ \
/*	distance between fft data elements is 512K */ \
/*	do 128 macros each processing 24 data values */ \
\
start_timer(19); \
rdi = (uintptr_t)g->sincos4;/* Load sin/cos pointer */ \
/* 32 iterations */ \
xsix_reals_unfft_preload; \
\
for(unsigned int loopA = 4; loopA; loopA--) { \
	for(unsigned int loopB = 8; loopB; loopB--) { \
		xtouchlines(3);			/* Preload TLBs */ \
		\
		for(unsigned int loopC = clm; loopC; loopC--) { \
			for(unsigned int loopD = 2; loopD; loopD--) { \
				x4cl_eight_reals_unfft_2(rsi, 64, 4*clmblkdst8, 8*clmblkdst8); \
			} \
			\
			rsi += -2*64+16*clmblkdst8;/* Next source pointer */ \
			for(unsigned int loopD = 2; loopD; loopD--) { \
				for(unsigned int loopE = 2; loopE; loopE--) { \
					x4cl_four_complex_unfft(rsi, 64, 4*clmblkdst8, 8*clmblkdst8, rdi); \
				} \
				rsi += -2*64+16*clmblkdst8;/* Next source pointer */ \
				rdi += XMM_SCD;		/* Next sine/cosine pointer */ \
			} \
			L2prefetch128(u8ptr(rcx+rbp)); \
			pfing(rcx += blkdst);		/* Next prefetch pointer */ \
			rdi += -2*XMM_SCD;		/* Restore sine/cosine pointer */ \
			\
			rsi += -3*16*clmblkdst8;	/* Restore source pointer */ \
			for(unsigned int loopD = 2; loopD; loopD--) { \
				for(unsigned int loopE = 2; loopE; loopE--) { \
					for(unsigned int loopF = 2; loopF; loopF--) { \
						x3cl_six_reals_last_unfft(rsi, 64, 16*clmblkdst8); \
					} \
					rsi += -2*64+4*clmblkdst8;	/* Next source pointer */ \
				} \
				L2prefetch128(u8ptr(rcx+rbp)); \
				pfing(rcx += blkdst);		/* Next prefetch pointer */ \
			} \
			\
				EXEC248##clm(rsi += -4*4*clmblkdst8+2*64);/* Next source pointer */ \
			EXEC248##clm(pfing(rcx += -3*blkdst+128));	/* Next prefetch pointer */ \
		} \
		EXEC248##clm(rsi += -clm*2*64+clmblkdst);/* Next source pointer */ \
		EXEC248##clm(pfing(rcx += -clm*128+3*blkdst));	/* Next prefetch pointer */ \
		EXEC1##clm(rsi += -4*4*clmblkdst8+clmblkdst);/* Next source pointer */ \
		\
	} \
	rsi += -8*clmblkdst+clmblkdst8;/* Next source pointer */ \
} \
end_timer(19); \
\
/* Normalize these values */ \
\
pass1_normalize(1, yes_fft); \
\
/* Copy data from the scratch area back to the FFT data area */ \
\
copy_scratch_data_to_fft(384, clm); \
goto	endlp;			/* Skip FFT code if flag not set */ \
\
/* Fire up auxiliary threads to do pass 1 forward FFT */ \
\
/* pass1: */	pass1_forward_fft_setup(b0b); \
\
/* Copy data from the FFT data area back to the scratch area */ \
\
b0b:	copy_fft_data_to_scratch(384, clm); \
\
/* Do FFT levels 1,2,3 */ \
/* */ \
/* On input the 64-byte cache lines hold these data values: */ \
/*	0K	768K	1536K	2304K */ \
/*	... */ \
/*	766K */ \
/* On output the 64-byte cache lines hold these data values: */ \
/*	0K	256K	512K	768K */ \
/*	2K	... */ \
/*	... */ \
/*	254K */ \
/*	1024K	1536K	2048K	2560K */ \
/*	1026K */ \
/*	... */ \
/*	1534K */ \
\
/* Do 512 six_reals_first_fft macros (FFT levels 1,2,3) */ \
/*	distance between fft data elements is 512K */ \
/*	do 128 iterations each processing 24 data values */ \
\
yes_fft:start_timer(22); \
rsi = (uintptr_t)g->scratch_area;	/* Get address of scratch area */ \
rcx = (uintptr_t)g->data_addr; /* Load source address for prefetching */ \
pfing(rbp = (uintptr_t)g->data_prefetch);/* Set prefetch distance */ \
pfing(rbp -= rcx); \
pfing(rcx += 186*blkdst);		/* Create prefetching pointer */ \
/* 16 iterations of 4 iters of 2*2*clm */ \
xsix_reals_fft_preload; \
for(unsigned int loopA = 16; loopA; loopA--) { \
	for(unsigned int loopB = 4; loopB; loopB--) { \
		xtouch(xptr(rcx+(clm-1)*128));	/* Preload TLB */ \
		xloop2_clm##clm(b1b, clmblkdst, x3cl_six_reals_fft(rsi, 64, 16*clmblkdst8)); \
	} \
	rsi += -8*clmblkdst+clmblkdst8; \
} \
rsi += -16*clmblkdst8;	/* Restore source pointer */ \
end_timer(22); \
\
/* Do FFT levels 4,5 */ \
/* Values 0K-511K is real data, 512K-1023K is semi-real data, 1024K-3071K is */ \
/* complex data. */ \
/* */ \
/* On input the 64-byte cache lines hold these data values: */ \
/*	0K	256K	512K	768K */ \
/*	2K	... */ \
/*	... */ \
/*	254K */ \
/*	1024K	1536K	2048K	2560K */ \
/*	1026K */ \
/*	... */ \
/*	1534K */ \
/* On output the 64-byte cache lines hold these data values: */ \
/*	0K	128K	256K	384K */ \
/*	2K	... */ \
/*	... */ \
/*	126K	... */ \
/*	512K	... */ \
/*	... */ \
\
/* Do 32 eight_reals_fft_2 macros */ \
/* Do 32 nop_two_two_complex_fft_2 macros */ \
/* Do 64 four_complex_fft macros */ \
/*	distance between fft data elements is 128K */ \
/*	do 64 macros each processing 16 data values */ \
\
start_timer(23); \
/* 8 iterations of 2 iters of 4*2*clm */ \
for(unsigned int loopA = 8; loopA; loopA--) { \
	for(unsigned int loopB = 2; loopB; loopB--) { \
		xtouch(xptr(rcx+(clm-1)*128));	/* Preload TLB */ \
		xloop4_clm##clm(b2b, clmblkdst, x2cl_eight_reals_fft_2(rsi, 64, 8*clmblkdst8)); \
	} \
	rsi += -8*clmblkdst+clmblkdst8; \
} \
\
/* Do 256 four_complex_fft macros */ \
/*	distance between fft data elements is 128K */ \
/*	do 64 macros each processing 32 data values */ \
\
rsi += -8*clmblkdst8+16*clmblkdst8;/* Next source pointer */ \
rdi = (uintptr_t)g->sincos4;/* Load sin/cos pointer */ \
/* 8 iters of 4 iters of 2*2*clm */ \
for(unsigned int loopA = 8; loopA; loopA--) { \
	for(unsigned int loopB = 4; loopB; loopB--) { \
		xtouch(xptr(rcx+(clm-1)*128));	/* Preload TLB */ \
		xloop2_clm##clm(b3b, clmblkdst, x4cl_four_complex_fft(rsi, 64, 8*clmblkdst8, 16*clmblkdst8, rdi)); \
	} \
	rsi += -8*clmblkdst+clmblkdst8; \
} \
rsi += -8*clmblkdst8-16*clmblkdst8;/* Restore src pointer */ \
end_timer(23); \
\
/* Do FFT levels 6,7 */ \
/* Values 0K-127K is real data, 128K-255K is semi-real data, 256K-4095K is */ \
/* complex data. */ \
/* */ \
/* On input the 64-byte cache lines hold these data values: */ \
/*	0K	128K	256K	384K */ \
/*	2K	... */ \
/*	... */ \
/*	126K	... */ \
/*	512K	... */ \
/*	... */ \
/* On output the 64-byte cache lines hold these data values: */ \
/*	0K	32K	64K	96K */ \
/*	2K	... */ \
/*	... */ \
/*	30K	... */ \
/*	128K	... */ \
/*	... */ \
\
/* Do 32 eight_reals_fft_2 macros */ \
/* Do 32 two_two_complex_fft macros */ \
/*	distance between fft data elements is 32K */ \
/*	do 16 macros each processing 32 data values */ \
\
start_timer(24); \
rdi = (uintptr_t)g->sincos3;/* Load sin/cos pointer */ \
/* 2 iters of 8 iters of 2*clm */ \
for(unsigned int loopA = 2; loopA; loopA--) { \
	for(unsigned int loopB = 8; loopB; loopB--) { \
		for(unsigned int loopC = 2*clm; loopC; loopC--) { \
			x4cl_eight_reals_fft_2(rsi, 64, 2*clmblkdst8, 4*clmblkdst8, rdi); \
		} \
		rsi += -2*clm*64+clmblkdst;/* Next source pointer */ \
	} \
	rsi += -8*clmblkdst+clmblkdst8; \
} \
\
/* Do 320 four_complex_fft macros */ \
/*	distance between fft data elements is 8K */ \
/*	do 80 macros each processing 32 data values */ \
\
rsi += -2*clmblkdst8+8*clmblkdst8;/* Load source pointer */ \
rdi += XMM_SCD;		/* Next sine/cosine pointer */ \
/* 5 iterations of 8 iters of 2*2*clm */ \
for(unsigned int loopA = 5; loopA; loopA--) { \
	for(unsigned int loopB = 2; loopB; loopB--) { \
		for(unsigned int loopC = 4; loopC; loopC--) { \
			xtouch(xptr(rcx+(clm-1)*128));	/* Preload TLB */ \
			xloop2_clm##clm(b5b, clmblkdst, x4cl_four_complex_fft(rsi, 64, 2*clmblkdst8, 4*clmblkdst8, rdi)); \
		} \
		rsi += -8*clmblkdst+clmblkdst8; \
	} \
	rsi += -2*clmblkdst8+8*clmblkdst8;/* Next source pointer */ \
	rdi += 2*XMM_SCD;		/* Next sine/cosine pointer */ \
} \
rsi += -6*8*clmblkdst8;	/* Restore source pointer */ \
end_timer(24); \
\
/* Do FFT levels 8,9 */ \
/* Values 0K-31K is real data, 32K-63K is semi-real data, 64K-3071K is */ \
/* complex data. */ \
/* */ \
/* On input the 64-byte cache lines hold these data values: */ \
/*	0K	32K	64K	96K */ \
/*	2K	... */ \
/*	... */ \
/*	30K	... */ \
/*	128K	... */ \
/*	... */ \
/* On output the 64-byte cache lines hold these data values: */ \
/*	0K	8K	16K	24K */ \
/*	2K	... */ \
/*	... */ \
/*	6K	... */ \
/*	32K	... */ \
/*	... */ \
\
/* Do 8 eight_reals_fft_2 macros */ \
/* Do 8 two_two_complex_fft macros */ \
/*	distance between fft data elements is 8K */ \
/*	do 4 macros each processing 32 data values */ \
\
start_timer(25); \
rdi = (uintptr_t)g->sincos2;/* Load sin/cos pointer */ \
/* 4 iterations of 2*clm */ \
for(unsigned int loopA = 4; loopA; loopA--) { \
	for(unsigned int loopB = 2*clm; loopB; loopB--) { \
		x4cl_eight_reals_fft_2(rsi, 64, 4*clmblkdst, clmblkdst8, rdi); \
	} \
	rsi += -2*clm*64+clmblkdst;/* Next source pointer */ \
} \
\
/* Do 368 two_two_complex_fft macros */ \
/*	distance between fft data elements is 8K */ \
/*	do 92 macros each processing 32 data values */ \
\
rsi += -4*clmblkdst+2*clmblkdst8;/* Load source pointer */ \
rdi += XMM_SCD;		/* Next sine/cosine pointer */ \
/* 23 iterations of 2 iters of 2*2*clm */ \
for(unsigned int loopA = 23; loopA; loopA--) { \
	for(unsigned int loopB = 2; loopB; loopB--) { \
		xtouch(xptr(rcx+(clm-1)*128));	/* Preload TLB */ \
		xloop2_clm##clm(b7b, clmblkdst, x4cl_four_complex_fft(rsi, 64, 4*clmblkdst, clmblkdst8, rdi)); \
	} \
	rsi += -2*2*clmblkdst+2*clmblkdst8;/* Next source pointer */ \
	rdi += 2*XMM_SCD;		/* Next sine/cosine pointer */ \
} \
rsi += -24*2*clmblkdst8;	/* Restore source pointer */ \
pfing(rcx += -384*blkdst);	/* Restore prefetch pointer */ \
end_timer(25); \
\
/* Do FFT levels 10,11 */ \
/* Values 0K-7K is real data, 8K-15K is semi-real data, 16K-3071K is */ \
/* complex data. */ \
/* */ \
/* On input the 64-byte cache lines hold these data values: */ \
/*	0K	8K	16K	24K */ \
/*	2K	... */ \
/*	... */ \
/*	6K	... */ \
/*	32K	... */ \
/*	... */ \
/* On output the 64-byte cache lines hold these data values: */ \
/*	0K	2K	4K	6K */ \
/*	8K	... */ \
/*	... */ \
\
/* Do 1 eight_reals_fft_2 macros */ \
/* Do 1 nop_two_two_complex_fft_2 macros */ \
/* Do 2 four_complex_fft macros */ \
/*	distance between fft data elements is 2K */ \
/*	do 1 macros each processing 32 data values */ \
\
start_timer(26); \
rdi = (uintptr_t)g->sincos1;/* Load sin/cos pointer */ \
for(unsigned int loopA = 2*clm; loopA; loopA--) { \
	g4cl_eight_reals_fft_2(rsi, 64, clmblkdst, 2*clmblkdst, rcx, 64, blkdst, 2*blkdst); \
} \
\
/* Do 380 four_complex_fft macros */ \
/*	distance between fft data elements is 2K */ \
/*	do 95 macros each processing 32 data values */ \
\
rsi += -2*clm*64+4*clmblkdst;/* Next source pointer */ \
rcx += -2*clm*64+4*blkdst;	/* Next dest pointer */ \
rdi += XMM_SCD;		/* Next sine/cosine pointer */ \
/* eax = 48*256+1; */ /* 95 iterations of 2*clm */ \
for(unsigned int loopA = 95; loopA; ) { \
	do{ \
		for(unsigned int loopC = clm; loopC; loopC--) { \
			for(unsigned int loopD = 2; loopD; loopD--) { \
				g4cl_four_complex_fft(rsi, 64, clmblkdst, 2*clmblkdst, rcx, 64, blkdst, 2*blkdst); \
				L2prefetch128(u8ptr(rcx-64+rbp)); \
				L2prefetch128(u8ptr(rcx-64+blkdst+rbp)); \
				pfing(rbp += -64+2*blkdst); \
			} \
			pfing(rbp += 128-4*blkdst); \
		} \
		rsi += -2*clm*64+4*clmblkdst;/* Next source pointer */ \
		rcx += -2*clm*64+4*blkdst;	/* Next dest pointer */ \
		rdi += 2*XMM_SCD;		/* Next sine/cosine pointer */ \
	}while(--loopA & 1); \
	rsi += -8*clmblkdst+clmblkdst8; \
} \
end_timer(26); \
\
/* Work on next 4 sets of 1536 values */ \
\
endlp:	xpass1_get_next_block(pass2, c0b, b0b); /* Get next block and jump */ \
/* done: */	xfft_footer_hg2; \
\
\
\
/* Perform a 3M-element negacyclic FFT.  This is done in two passes.  Pass 1 */ \
/* does 11 levels, pass 2 uses common code to perform the last 11 levels. */ \
\

#define hg_pass1sclevels11complex3(procname, xpass2name, pass2_levels, clm) \
 /* (fftname, pass2_levels, clmarg) */ \
\
/* set_FFT_constants(fftname, 0, 1, pass2_levels, clmarg); */ \
xfft_header_hg2(procname, xpass2name, pass2_levels, clm); \
\
/* Do the last 11 FFT levels */ \
/*	do 384 pass2_11_levels_complex */ \
\
/* pass2: */	start_timer(0); \
xpass2name(g); \
end_timer(0); \
if(g->ffttype == 1)		/* We're done if FFTing only */ \
	return; \
\
/* Do the inverse FFT */ \
\
pass1_inverse_fft_setup(c0b, procname); \
\
/* Do inverse FFT levels 10,11 */ \
/* On input the 64-byte cache lines hold these data values: */ \
/*	0K	2K	4K	6K */ \
/*	8K	10K	12K	14K */ \
/*	... */ \
/* On output the 64-byte cache lines hold these data values: */ \
/*	0K	4K	8K	12K */ \
/*	2K	... */ \
/*	16K	... */ \
/*	... */ \
\
/* Do 384 four_complex_unfft macros */ \
/*	distance between fft data elements is 2K */ \
/*	do 192 macros each processing 16 data values */ \
\
c0b:	start_timer(16); \
rsi = (uintptr_t)g->data_addr;/* Get source address */ \
rcx = (uintptr_t)g->scratch_area;	/* Get address of scratch area */ \
rdi = (uintptr_t)g->sincos1;/* Load sin/cos pointer */ \
/* 48 iterations of 4 of 2*clm */ \
pfing(rbp = (uintptr_t)g->data_prefetch);/* Set prefetch distance */ \
pfing(rbp -= rsi); \
for(unsigned int loopA = 48; loopA; loopA--) { \
	for(unsigned int loopB = 4; loopB; loopB--) { \
		for(unsigned int loopC = 2*clm; loopC; loopC--) { \
			g2cl_four_complex_unfft(rsi, 64, blkdst, rcx, 64, clmblkdst); \
		} \
		rsi += -2*clm*64+2*blkdst;	/* Next source pointer */ \
		rcx += -2*clm*64+2*clmblkdst;/* Next dest pointer */ \
		rdi += XMM_SCD;		/* Next sine/cosine pointer */ \
	} \
	rcx += -8*clmblkdst+clmblkdst8; \
} \
rcx = rsi-192*2*blkdst;	/* Create prefetch pointer */ \
rsi = (uintptr_t)g->scratch_area;	/* Restore source pointer */ \
end_timer(16); \
\
/* Do inverse FFT levels 8,9 */ \
/* On input the 64-byte cache lines hold these data values: */ \
/*	0K	4K	8K	12K */ \
/*	2K	... */ \
/*	16K	... */ \
/*	... */ \
/* On output the 64-byte cache lines hold these data values: */ \
/*	0K	16K	32K	48K */ \
/*	2K	... */ \
/*	... */ \
/*	14K */ \
/*	64K	... */ \
/*	... */ \
\
/* Do 384 four_complex_unfft macros */ \
/*	distance between fft data elements is 8K */ \
/*	do 96 macros each processing 32 data values */ \
\
start_timer(17); \
rdi = (uintptr_t)g->sincos2;/* Load sin/cos pointer */ \
/* 48 iterations of 2*2*clm */ \
for(unsigned int loopA = 48; loopA; loopA--) { \
	xtouch(xptr(rcx+(clm-1)*128));	/* Preload TLB */ \
	xloop2_clm##clm(c7b, clmblkdst, x4cl_four_complex_unfft(rsi, 64, 2*clmblkdst, 4*clmblkdst, rdi)); \
	rsi += -2*clmblkdst+clmblkdst8;/* Next source pointer */ \
	rdi += XMM_SCD;		/* Next sine/cosine pointer */ \
} \
rsi += -48*clmblkdst8;	/* Restore source pointer */ \
end_timer(17); \
\
/* Do inverse FFT levels 6,7 */ \
/* On input the 64-byte cache lines hold these data values: */ \
/*	0K	16K	32K	48K */ \
/*	2K	... */ \
/*	... */ \
/*	14K */ \
/*	64K	... */ \
/*	... */ \
/* On output the 64-byte cache lines hold these data values: */ \
/*	0K	64K	128K	192K */ \
/*	2K	... */ \
/*	... */ \
/*	62K */ \
/*	256K	... */ \
/*	... */ \
\
/* Do 384 four_complex_unfft macros */ \
/*	distance between fft data elements is 8K */ \
/*	do 96 macros each processing 32 data values */ \
\
start_timer(18); \
rdi = (uintptr_t)g->sincos3;/* Load sin/cos pointer */ \
/* 12 iterations of 4 iters of 2*2*clm */ \
for(unsigned int loopA = 12; loopA; loopA--) { \
	for(unsigned int loopB = 4; loopB; loopB--) { \
		xtouch(xptr(rcx+(clm-1)*128));	/* Preload TLB */ \
		xloop2_clm##clm(c5b, clmblkdst, x4cl_four_complex_unfft(rsi, 64, clmblkdst8, 2*clmblkdst8, rdi));		/* 2*2*clm iterations */ \
	} \
	rsi += -4*2*clmblkdst+4*clmblkdst8;/* Next source pointer */ \
	rdi += XMM_SCD;		/* Next sine/cosine pointer */ \
} \
rsi += -12*4*clmblkdst8;	/* Restore source pointer */ \
end_timer(18); \
\
/* Do inverse FFT levels 4,5 */ \
/* On input the 64-byte cache lines hold these data values: */ \
/*	0K	64K	128K	192K */ \
/*	2K	... */ \
/*	... */ \
/*	62K */ \
/*	256K	... */ \
/*	... */ \
/* On output the 64-byte cache lines hold these data values: */ \
/*	0K	256K	512K	768K */ \
/*	2K	... */ \
/*	... */ \
/*	254K */ \
/*	1024K	... */ \
/*	... */ \
\
/* Do 384 four_complex_unfft macros */ \
/*	distance between fft data elements is 128K */ \
/*	do 96 macros each processing 32 data values */ \
\
start_timer(19); \
rdi = (uintptr_t)g->sincos4;/* Load sin/cos pointer */ \
/* 3 iters of 4 iters of 4 of 2*2*clm */ \
for(unsigned int loopA = 3; loopA; loopA--) { \
	for(unsigned int loopB = 4; loopB; loopB--) { \
		xtouch(xptr(rcx+(clm-1)*128));	/* Preload TLB */ \
		for(unsigned int loopC = 4; loopC; loopC--) { \
			xloop2_clm##clm(c3b, clmblkdst, x4cl_four_complex_unfft(rsi, 64, 4*clmblkdst8, 8*clmblkdst8, rdi));		/* 2* 2*clm iterations */ \
		} \
		rsi += -4*2*clmblkdst+clmblkdst8;/* Next source pointer */ \
	} \
	rsi += -4*clmblkdst8+16*clmblkdst8;/* Next source pointer */ \
	rdi += XMM_SCD;		/* Next sine/cosine pointer */ \
} \
rsi += -3*16*clmblkdst8;	/* Restore source pointer */ \
end_timer(19); \
\
/* Do inverse FFT levels 1,2,3 */ \
/* On input the 64-byte cache lines hold these data values: */ \
/*	0K	256K	512K	768K */ \
/*	2K	... */ \
/*	... */ \
/*	254K */ \
/*	1024K	1536K	2048K	2560K */ \
/*	1026K */ \
/*	... */ \
/*	1534K */ \
/* On output the 64-byte cache lines hold these data values: */ \
/*	0K	768K	1536K	2304K */ \
/*	... */ \
/*	766K */ \
/* */ \
/* Do 512 three_complex_last_unfft macros */ \
/*	distance between fft data elements is 512K */ \
/*	do 128 macros each processing 24 data values */ \
\
start_timer(20); \
rdi = (uintptr_t)g->plus1_premults;/* Address of premultiplier table */ \
/* 16 iters of 4 iters of 2*2*clm */ \
for(unsigned int loopA = 16; loopA; loopA--) { \
	for(unsigned int loopB = 4; loopB; loopB--) { \
		xtouch(xptr(rcx+(clm-1)*128));	/* Preload TLB */ \
		xloop2edi_clm##clm(c1b, clmblkdst, 192, x3cl_three_complex_last_unfft(rsi, 64, 16*clmblkdst8)); /* 2*2*clm iterations */ \
	} \
	rsi += -8*clmblkdst+clmblkdst8;/* Next source pointer */ \
} \
end_timer(20); \
\
/* Normalize these values */ \
\
pass1_normalize(1, yes_fft); \
\
/* Copy data from the scratch area back to the FFT data area */ \
\
copy_scratch_data_to_fft(384, clm); \
goto	endlp;			/* Skip FFT code if flag not set */ \
\
/* Fire up auxiliary threads to do pass 1 forward FFT */ \
\
/* pass1: */	pass1_forward_fft_setup(b0b); \
\
/* Copy data from the FFT data area back to the scratch area */ \
\
b0b:	copy_fft_data_to_scratch(384, clm); \
\
/* Do FFT levels 1,2,3 */ \
/* */ \
/* On input the 64-byte cache lines hold these data values: */ \
/*	0K	768K	1536K	2304K */ \
/*	... */ \
/*	766K */ \
/* On output the 64-byte cache lines hold these data values: */ \
/*	0K	256K	512K	768K */ \
/*	2K	... */ \
/*	... */ \
/*	254K */ \
/*	1024K	1536K	2048K	2560K */ \
/*	1026K */ \
/*	... */ \
/*	1534K */ \
\
/* Do 512 three_complex_first_fft macros (FFT levels 1,2,3) */ \
/*	distance between fft data elements is 512K */ \
/*	do 128 iterations each processing 24 data values */ \
\
yes_fft:start_timer(22); \
rsi = (uintptr_t)g->scratch_area;	/* Get address of scratch area */ \
rdi = (uintptr_t)g->plus1_premults;/* Address of premultiplier table */ \
rcx = (uintptr_t)g->data_addr; /* Load source address for prefetching */ \
pfing(rbp = (uintptr_t)g->data_prefetch);/* Set prefetch distance */ \
pfing(rbp -= rcx); \
pfing(rcx += 208*blkdst);		/* Create prefetching pointer */ \
/* 16 iterations of 2 iters of 4*2*clm */ \
for(unsigned int loopA = 16; loopA; loopA--) { \
	for(unsigned int loopB = 2; loopB; loopB--) { \
		xtouch(xptr(rcx+(clm-1)*128));	/* Preload TLB */ \
		xloop4edi_clm##clm(b1b, clmblkdst, 192, x3cl_three_complex_fft(rsi, 64, 16*clmblkdst8)); /* 4*2*clm iterations */ \
	} \
	rsi += -2*4*clmblkdst+clmblkdst8; \
} \
rsi += -16*clmblkdst8;	/* Restore source pointer */ \
end_timer(22); \
\
/* Do FFT levels 4,5 */ \
/* */ \
/* On input the 64-byte cache lines hold these data values: */ \
/*	0K	256K	512K	768K */ \
/*	2K	... */ \
/*	... */ \
/*	254K */ \
/*	1024K	1536K	2048K	2560K */ \
/*	1026K */ \
/*	... */ \
/*	1534K */ \
/* On output the 64-byte cache lines hold these data values: */ \
/*	0K	128K	256K	384K */ \
/*	2K	... */ \
/*	... */ \
/*	126K	... */ \
/*	512K	... */ \
/*	... */ \
\
/* Do 384 four_complex_fft macros */ \
/*	distance between fft data elements is 128K */ \
/*	do 32 macros each processing 16 data values */ \
/*	do 64 macros each processing 32 data values */ \
\
start_timer(23); \
rdi = (uintptr_t)g->sincos4;/* Load sin/cos pointer */ \
/* 8 iterations of 2 iters of 4*2*clm */ \
for(unsigned int loopA = 8; loopA; loopA--) { \
	for(unsigned int loopB = 2; loopB; loopB--) { \
		xtouch(xptr(rcx+(clm-1)*128));	/* Preload TLB */ \
		xloop4_clm##clm(b2b, clmblkdst, x2cl_four_complex_fft(rsi, 64, 8*clmblkdst8));		/* 4*2*clm iterations */ \
	} \
	rsi += -8*clmblkdst+clmblkdst8; \
} \
rsi += -8*clmblkdst8+16*clmblkdst8;/* Next source pointer */ \
rdi += XMM_SCD;		/* Next sine/cosine pointer */ \
/* 8 iters of 4 iters of 2*2*clm */ \
for(unsigned int loopA = 8; loopA; loopA--) { \
	for(unsigned int loopB = 4; loopB; loopB--) { \
		xtouch(xptr(rcx+(clm-1)*128));	/* Preload TLB */ \
		xloop2_clm##clm(b3b, clmblkdst, x4cl_four_complex_fft(rsi, 64, 8*clmblkdst8, 16*clmblkdst8, rdi));		/* 2*2*clm iterations */ \
	} \
	rsi += -8*clmblkdst+clmblkdst8; \
} \
rsi += -8*clmblkdst8-16*clmblkdst8;/* Restore src pointer */ \
end_timer(23); \
\
/* Do FFT levels 6,7 */ \
/* */ \
/* On input the 64-byte cache lines hold these data values: */ \
/*	0K	128K	256K	384K */ \
/*	2K	... */ \
/*	... */ \
/*	126K	... */ \
/*	512K	... */ \
/*	... */ \
/* On output the 64-byte cache lines hold these data values: */ \
/*	0K	32K	64K	96K */ \
/*	2K	... */ \
/*	... */ \
/*	30K	... */ \
/*	128K	... */ \
/*	... */ \
\
/* Do 384 four_complex_fft macros */ \
/*	distance between fft data elements is 8K */ \
/*	do 96 macros each processing 32 data values */ \
\
start_timer(24); \
rdi = (uintptr_t)g->sincos3;/* Load sin/cos pointer */ \
/* 6 iters of 2 iters of 4 of 2*2*clm */ \
for(unsigned int loopA = 6; loopA; loopA--) { \
	for(unsigned int loopB = 2; loopB; loopB--) { \
		for(unsigned int loopC = 4; loopC; loopC--) { \
			xtouch(xptr(rcx+(clm-1)*128));	/* Preload TLB */ \
			xloop2_clm##clm(b5b, clmblkdst, x4cl_four_complex_fft(rsi, 64, 2*clmblkdst8, 4*clmblkdst8, rdi)); \
		} \
		rsi += -8*clmblkdst+clmblkdst8; \
	} \
	rsi += -2*clmblkdst8+8*clmblkdst8;/* Next source pointer */ \
	rdi += 2*XMM_SCD;		/* Next sine/cosine pointer */ \
} \
rsi += -6*8*clmblkdst8;	/* Restore source pointer */ \
end_timer(24); \
\
/* Do FFT levels 8,9 */ \
/* */ \
/* On input the 64-byte cache lines hold these data values: */ \
/*	0K	32K	64K	96K */ \
/*	2K	... */ \
/*	... */ \
/*	30K	... */ \
/*	128K	... */ \
/*	... */ \
/* On output the 64-byte cache lines hold these data values: */ \
/*	0K	8K	16K	24K */ \
/*	2K	... */ \
/*	... */ \
/*	6K	... */ \
/*	32K	... */ \
/*	... */ \
\
/* Do 384 two_two_complex_fft macros */ \
/*	distance between fft data elements is 8K */ \
/*	do 96 macros each processing 32 data values */ \
\
start_timer(25); \
rdi = (uintptr_t)g->sincos2;/* Load sin/cos pointer */ \
/* 24 iterations of 2 iters of 2*2*clm */ \
for(unsigned int loopA = 24; loopA; loopA--) { \
	for(unsigned int loopB = 2; loopB; loopB--) { \
		xtouch(xptr(rcx+(clm-1)*128));	/* Preload TLB */ \
		xloop2_clm##clm(b7b, clmblkdst, x4cl_four_complex_fft(rsi, 64, 4*clmblkdst, clmblkdst8, rdi));		/* 2*2*clm iterations */ \
	} \
	rsi += -2*2*clmblkdst+2*clmblkdst8;/* Next source pointer */ \
	rdi += 2*XMM_SCD;		/* Next sine/cosine pointer */ \
} \
rsi += -24*2*clmblkdst8;	/* Restore source pointer */ \
pfing(rcx += -384*blkdst);	/* Restore prefetch pointer */ \
end_timer(25); \
\
/* Do FFT levels 10,11 */ \
/* */ \
/* On input the 64-byte cache lines hold these data values: */ \
/*	0K	8K	16K	24K */ \
/*	2K	... */ \
/*	... */ \
/*	6K	... */ \
/*	32K	... */ \
/*	... */ \
/* On output the 64-byte cache lines hold these data values: */ \
/*	0K	2K	4K	6K */ \
/*	8K	... */ \
/*	... */ \
\
/* Do 384 four_complex_fft macros */ \
/*	distance between fft data elements is 2K */ \
/*	do 96 macros each processing 32 data values */ \
\
start_timer(26); \
rdi = (uintptr_t)g->sincos1;/* Load sin/cos pointer */ \
/* 48 iterations of 2 iters of 2*clm */ \
for(unsigned int loopA = 48; loopA; loopA--) { \
	for(unsigned int loopB = 2; loopB; loopB--) { \
		for(unsigned int loopC = clm; loopC; loopC--) { \
			for(unsigned int loopD = 2; loopD; loopD--) { \
				g4cl_four_complex_fft(rsi, 64, clmblkdst, 2*clmblkdst, rcx, 64, blkdst, 2*blkdst); \
				L2prefetch128(u8ptr(rcx-64+rbp)); \
				L2prefetch128(u8ptr(rcx-64+blkdst+rbp)); \
				pfing(rbp += -64+2*blkdst); \
			} \
			pfing(rbp += 128-4*blkdst); \
		} \
	rsi += -2*clm*64+4*clmblkdst;/* Next source pointer */ \
		rcx += -2*clm*64+4*blkdst;	/* Next dest pointer */ \
		rdi += 2*XMM_SCD;		/* Next sine/cosine pointer */ \
	} \
	rsi += -8*clmblkdst+clmblkdst8; \
} \
end_timer(26); \
\
/* Work on next 4 sets of 1536 values */ \
\
endlp:	xpass1_get_next_block(pass2, c0b, b0b); /* Get next block and jump */ \
/* done: */	xfft_footer_hg2; \
\
\
\
/* Perform a 3584K-element FFT.  This is done in two passes.  Pass 1 */ \
/* does 11 levels, pass 2 uses common code to perform the last 11 levels. */ \
\

#define hg_pass1sclevels11pfa7(procname, xpass2name, pass2_levels, clm) \
 /* (fftname, pass2_levels, clmarg) */ \
\
/* set_FFT_constants(fftname, 0, 0, pass2_levels, clmarg); */ \
xfft_header_hg2(procname, xpass2name, pass2_levels, clm); \
\
/* Do the last 11 FFT levels */ \
/*	do 1 pass2_11_levels_real and 447 pass2_11_levels_complex */ \
\
/* pass2: */	start_timer(0); \
xpass2name(g); \
end_timer(0); \
if(g->ffttype == 1)		/* We're done if FFTing only */ \
	return; \
\
/* Do the inverse FFT */ \
\
pass1_inverse_fft_setup(c0b, procname); \
\
/* Do inverse FFT levels 10,11 */ \
/* On input the 64-byte cache lines hold these data values: */ \
/*	0K	2K	4K	6K */ \
/*	8K	10K	12K	14K */ \
/*	... */ \
/* On output the 64-byte cache lines hold these data values: */ \
/*	0K	4K	8K	12K */ \
/*	2K	... */ \
/*	16K	... */ \
/*	... */ \
\
/* Do 1 eight_reals_unfft_2 macros */ \
/* Do 1 nop_two_two_complex_unfft_2 macros */ \
/*	distance between fft data elements is 2K */ \
/*	do 1 macros each processing 16 data values */ \
\
c0b:	start_timer(16); \
rsi = (uintptr_t)g->data_addr;/* Get source address */ \
rcx = (uintptr_t)g->scratch_area;	/* Get address of scratch area */ \
pfing(rbp = (uintptr_t)g->data_prefetch);/* Set prefetch distance */ \
pfing(rbp -= rsi); \
for(unsigned int loopA = 2*clm; loopA; loopA--) { \
	g2cl_eight_reals_unfft_2(rsi, 64, blkdst, rcx, 64, clmblkdst); \
} \
\
/* Do 446 four_complex_unfft macros */ \
/*	distance between fft data elements is 2K */ \
/*	do 223 macros each processing 16 data values */ \
\
rsi += -2*clm*64+2*blkdst;	/* Next source pointer */ \
rcx += -2*clm*64+2*clmblkdst;/* Next dest pointer */ \
rdi = (uintptr_t)g->sincos1;/* Load sin/cos pointer */ \
/* al = 256/4+56; */ /* 223 iterations of 2*clm */ \
for(unsigned int loopA = 223; loopA; ) { \
	do{ \
		for(unsigned int loopC = 2*clm; loopC; loopC--) { \
			g2cl_four_complex_unfft(rsi, 64, blkdst, rcx, 64, clmblkdst); \
		} \
		rsi += -2*clm*64+2*blkdst;	/* Next source pointer */ \
		rcx += -2*clm*64+2*clmblkdst;/* Next dest pointer */ \
		rdi += XMM_SCD;		/* Next sine/cosine pointer */ \
	}while(--loopA & 3); \
	rcx += -8*clmblkdst+clmblkdst8; \
} \
rcx = rsi-224*2*blkdst;	/* Create prefetch pointer */ \
rsi = (uintptr_t)g->scratch_area;	/* Restore source pointer */ \
end_timer(16); \
\
/* Do inverse FFT levels 8,9 */ \
/* On input the 64-byte cache lines hold these data values: */ \
/*	0K	4K	8K	12K */ \
/*	2K	... */ \
/*	16K	... */ \
/*	... */ \
/* On output the 64-byte cache lines hold these data values: */ \
/*	0K	16K	32K	48K */ \
/*	2K	... */ \
/*	... */ \
/*	14K */ \
/*	64K	... */ \
/*	... */ \
\
/* Do 4 eight_reals_unfft_2 macros */ \
/* Do 4 nop_two_two_complex_unfft_2 macros */ \
/*	distance between fft data elements is 8K */ \
/*	do 2 macros each processing 32 data values */ \
\
start_timer(17); \
/* 2 iterations of 2*clm */ \
for(unsigned int loopA = 2; loopA; loopA--) { \
	xtouch(xptr(rcx+(clm-1)*128));	/* Preload TLB */ \
	xloop1_clm##clm(c6b, clmblkdst, x4cl_eight_reals_unfft_2(rsi, 64, 2*clmblkdst, 4*clmblkdst)); \
} \
\
/* Do 440 four_complex_unfft macros */ \
/*	distance between fft data elements is 8K */ \
/*	do 110 macros each processing 32 data values */ \
\
rsi += -2*clmblkdst+clmblkdst8;/* Next source pointer */ \
rdi = (uintptr_t)g->sincos2;/* Load sin/cos pointer */ \
/* 55 iterations of 2*2*clm */ \
for(unsigned int loopA = 55; loopA; loopA--) { \
	xtouch(xptr(rcx+(clm-1)*128));	/* Preload TLB */ \
	xloop2_clm##clm(c7b, clmblkdst, x4cl_four_complex_unfft(rsi, 64, 2*clmblkdst, 4*clmblkdst, rdi)); \
	rsi += -2*clmblkdst+clmblkdst8;/* Next source pointer */ \
	rdi += XMM_SCD;		/* Next sine/cosine pointer */ \
} \
rsi += -56*clmblkdst8;	/* Restore source pointer */ \
end_timer(17); \
\
/* Do inverse FFT levels 6,7 */ \
/* On input the 64-byte cache lines hold these data values: */ \
/*	0K	16K	32K	48K */ \
/*	2K	... */ \
/*	... */ \
/*	14K */ \
/*	64K	... */ \
/*	... */ \
/* On output the 64-byte cache lines hold these data values: */ \
/*	0K	64K	128K	192K */ \
/*	2K	... */ \
/*	... */ \
/*	62K */ \
/*	256K	... */ \
/*	... */ \
\
/* Do 16 eight_reals_unfft_2 macros */ \
/* Do 16 nop_two_two_complex_unfft_2 macros */ \
/*	distance between fft data elements is 32K */ \
/*	do 8 macros each processing 32 data values */ \
\
start_timer(18); \
/* 8 iterations of 2*clm */ \
for(unsigned int loopA = 8; loopA; loopA--) { \
	xtouch(xptr(rcx+(clm-1)*128));	/* Preload TLB */ \
	xloop1_clm##clm(c4b, clmblkdst, x4cl_eight_reals_unfft_2(rsi, 64, clmblkdst8, 2*clmblkdst8)); \
} \
\
/* Do 416 four_complex_unfft macros */ \
/*	distance between fft data elements is 32K */ \
/*	do 104 macros each processing 32 data values */ \
\
rsi += -8*clmblkdst+4*clmblkdst8;/* Next source pointer */ \
rdi = (uintptr_t)g->sincos3;/* Load sin/cos pointer */ \
/* 13 iterations of 4 iters of 2*2*clm */ \
for(unsigned int loopA = 13; loopA; loopA--) { \
	for(unsigned int loopB = 4; loopB; loopB--) { \
		xtouch(xptr(rcx+(clm-1)*128));	/* Preload TLB */ \
		xloop2_clm##clm(c5b, clmblkdst, x4cl_four_complex_unfft(rsi, 64, clmblkdst8, 2*clmblkdst8, rdi)); \
	} \
	rsi += -4*2*clmblkdst+4*clmblkdst8;/* Next source pointer */ \
	rdi += XMM_SCD;		/* Next sine/cosine pointer */ \
} \
rsi += -14*4*clmblkdst8;	/* Restore source pointer */ \
end_timer(18); \
\
/* Do inverse FFT levels 4,5 */ \
/* On input the 64-byte cache lines hold these data values: */ \
/*	0K	64K	128K	192K */ \
/*	2K	... */ \
/*	... */ \
/*	62K */ \
/*	256K	... */ \
/*	... */ \
/* On output the 64-byte cache lines hold these data values: */ \
/*	0K	128K	256K	384K */ \
/*	2K	... */ \
/*	... */ \
/*	126K */ \
/*	512K	768K	1024K	1280K */ \
/*	... */ \
/* */ \
/* Do 64 eight_reals_unfft_2 macros */ \
/*	distance between fft data elements is 128K */ \
/*	do 32 macros each processing 16 data values */ \
/* */ \
/* Do 384 four_complex_unfft macros */ \
/*	distance between fft data elements is 128K */ \
/*	do 96 macros each processing 32 data values */ \
\
/* Do inverse FFT levels 1,2,3 */ \
/* On input the 64-byte cache lines hold these data values: */ \
/*	0K	128K	256K	384K */ \
/*	2K	... */ \
/*	... */ \
/*	126K */ \
/*	512K	768K	1024K	1280K */ \
/*	... */ \
/* On output the 64-byte cache lines hold these data values: */ \
/*	0K	128K	1792K	1920K */ \
/*	2K	... */ \
/*	... */ \
/*	126K */ \
/*	256K	512K	2048K	2304K */ \
/*	258K */ \
/*	... */ \
/*	510K */ \
/*	768K	1280K	2560K	3072K */ \
/*	770K */ \
/*	... */ \
/*	1278K */ \
/* */ \
/* Do 512 seven_reals_last_unfft macros */ \
/*	distance between fft data elements is 512K */ \
/*	do 64 macros each processing 56 data values */ \
\
start_timer(19); \
rdi = (uintptr_t)g->sincos4;/* Load sin/cos pointer */ \
/* 32 iterations */ \
xseven_reals_unfft_preload; \
\
for(unsigned int loopA = 4; loopA; loopA--) { \
	for(unsigned int loopB = 8; loopB; loopB--) { \
		xtouchlines(5);			/* Preload TLBs */ \
		\
		for(unsigned int loopC = clm; loopC; loopC--) { \
			for(unsigned int loopD = 2; loopD; loopD--) { \
				x2cl_half_eight_reals_unfft_2(rsi, 64, 4*clmblkdst8); \
			} \
			\
			rsi += -2*64+8*clmblkdst8;	/* Next source pointer */ \
			/* 3 iterations */ \
			for(unsigned int loopD = 3; loopD; loopD--) { \
				for(unsigned int loopE = 2; loopE; loopE--) { \
					x4cl_four_complex_unfft(rsi, 64, 4*clmblkdst8, 8*clmblkdst8, rdi); \
				} \
				rsi += -2*64+16*clmblkdst8;/* Next source pointer */ \
				rdi += XMM_SCD;		/* Next sine/cosine pointer */ \
			} \
			L2prefetch128(u8ptr(rcx+rbp)); \
			pfing(rcx += blkdst);		/* Next prefetch pointer */ \
			rdi += -3*XMM_SCD;		/* Restore sine/cosine pointer */ \
			\
			rsi += -3*16*clmblkdst8-8*clmblkdst8;/* Restore src pointer */ \
			for(unsigned int loopD = 2; loopD; loopD--) { \
				; \
				xloop1prefetch2sp(c1b, 4*clmblkdst8, x7cl_seven_reals_last_unfft(rsi, 64, 8*clmblkdst8)); \
			} \
			\
				EXEC248##clm(rsi += -2*4*clmblkdst8+2*64);/* Next source pointer */ \
			EXEC248##clm(pfing(rcx += -5*blkdst+128));	/* Next prefetch pointer */ \
		} \
		EXEC248##clm(rsi += -clm*2*64+clmblkdst);/* Next source pointer */ \
		EXEC248##clm(pfing(rcx += -clm*128+5*blkdst));	/* Next prefetch pointer */ \
		EXEC1##clm(rsi += -2*4*clmblkdst8+clmblkdst);/* Next source pointer */ \
		\
	} \
	rsi += -8*clmblkdst+clmblkdst8;/* Next source pointer */ \
} \
end_timer(19); \
\
/* Normalize these values */ \
\
pass1_normalize(1, yes_fft); \
\
/* Copy data from the scratch area back to the FFT data area */ \
\
copy_scratch_data_to_fft(448, clm); \
goto	endlp;			/* Skip FFT code if flag not set */ \
\
/* Fire up auxiliary threads to do pass 1 forward FFT */ \
\
/* pass1: */	pass1_forward_fft_setup(b0b); \
\
/* Copy data from the FFT data area back to the scratch area */ \
\
b0b:	copy_fft_data_to_scratch(448, clm); \
\
/* Do FFT levels 1,2,3 */ \
/* */ \
/* On input the 64-byte cache lines hold these data values: */ \
/*	0K	128K	1792K	1920K */ \
/*	2K	... */ \
/*	... */ \
/*	126K */ \
/*	256K	512K	2048K	2304K */ \
/*	258K */ \
/*	... */ \
/*	510K */ \
/*	768K	1280K	2560K	3072K */ \
/*	770K */ \
/*	... */ \
/*	1278K */ \
/* On output the 64-byte cache lines hold these data values: */ \
/*	0K	128K	256K	384K */ \
/*	2K	... */ \
/*	... */ \
/*	126K */ \
/*	512K	768K	1024K	1280K */ \
/*	514K */ \
/*	... */ \
/*	766K */ \
/*	1536K	2048K	2560K	3072K */ \
/*	1538K */ \
/*	... */ \
/*	2046K */ \
/* */ \
/* Do 512 seven_reals_first_fft macros (FFT levels 1,2,3) */ \
/*	distance between fft data elements is 512K */ \
/*	do 64 iterations each processing 56 data values */ \
\
yes_fft:start_timer(22); \
rsi = (uintptr_t)g->scratch_area;	/* Get address of scratch area */ \
rcx = (uintptr_t)g->data_addr; /* Load source address for prefetching */ \
pfing(rbp = (uintptr_t)g->data_prefetch);/* Set prefetch distance */ \
pfing(rbp -= rcx); \
pfing(rcx += 162*blkdst);		/* Create prefetching pointer */ \
/* 8 iterations 8 iters of 2*clm */ \
xseven_reals_fft_preload; \
for(unsigned int loopA = 8; loopA; loopA--) { \
	for(unsigned int loopB = 8; loopB; loopB--) { \
		xtouch(xptr(rcx+(clm-1)*128));	/* Preload TLB */ \
		xtouch(xptr(rcx+(clm-1)*128+blkdst)); \
		xloop1prefetch2_clm##clm(b1b, clmblkdst, x7cl_seven_reals_fft(rsi, 64, 8*clmblkdst8)); \
	} \
	rsi += -8*clmblkdst+clmblkdst8; \
} \
rsi += -8*clmblkdst8;	/* Restore source pointer */ \
end_timer(22); \
\
/* Do FFT levels 4,5 */ \
/* Values 0K-511K is real data, 512K-1023K is semi-real data, 1024K-3583K is */ \
/* complex data. */ \
/* */ \
/* On input the 64-byte cache lines hold these data values: */ \
/*	0K	128K	256K	384K */ \
/*	2K	... */ \
/*	... */ \
/*	126K */ \
/*	512K	768K	1024K	1280K */ \
/*	514K */ \
/*	... */ \
/*	766K */ \
/*	1536K	2048K	2560K	3072K */ \
/*	1538K */ \
/*	... */ \
/*	2046K */ \
/* On output the 64-byte cache lines hold these data values: */ \
/*	0K	128K	256K	384K */ \
/*	2K	... */ \
/*	... */ \
/*	126K	... */ \
/*	512K	... */ \
/*	... */ \
\
/* Do 64 eight_reals_fft_2 macros */ \
/*	distance between fft data elements is 128K */ \
/*	do 64 macros each processing 8 data values */ \
\
start_timer(23); \
/* 8 iterations of 8 iters of 2*clm */ \
for(unsigned int loopA = 8; loopA; loopA--) { \
	for(unsigned int loopB = 8; loopB; loopB--) { \
		for(unsigned int loopC = 2*clm; loopC; loopC--) { \
			x1cl_half_eight_reals_fft_2(rsi, 64); \
		} \
		rsi += -2*clm*64+clmblkdst;/* Next source pointer */ \
	} \
	rsi += -8*clmblkdst+clmblkdst8; \
} \
\
/* Do 128 four_complex_fft macros */ \
/*	distance between fft data elements is 128K */ \
/*	do 64 macros each processing 16 data values */ \
\
rdi = (uintptr_t)g->sincos4;/* Load sin/cos pointer */ \
/* 8 iterations of 2 iters of 4*2*clm */ \
for(unsigned int loopA = 8; loopA; loopA--) { \
	for(unsigned int loopB = 2; loopB; loopB--) { \
		xtouch(xptr(rcx+(clm-1)*128));	/* Preload TLB */ \
		xloop4_clm##clm(b3a, clmblkdst, x2cl_four_complex_fft(rsi, 64, 8*clmblkdst8)); \
	} \
	rsi += -8*clmblkdst+clmblkdst8; \
} \
\
/* Do 256 four_complex_fft macros */ \
/*	distance between fft data elements is 128K */ \
/*	do 64 macros each processing 32 data values */ \
\
rsi += -8*clmblkdst8+16*clmblkdst8;/* Next source pointer */ \
rdi += XMM_SCD;		/* Next sine/cosine pointer */ \
/* 8 iterations of 4 of 2*2*clm */ \
for(unsigned int loopA = 8; loopA; loopA--) { \
	for(unsigned int loopB = 4; loopB; loopB--) { \
		xtouch(xptr(rcx+(clm-1)*128));	/* Preload TLB */ \
		xloop2_clm##clm(b3b, clmblkdst, x4cl_four_complex_fft(rsi, 64, 8*clmblkdst8, 16*clmblkdst8, rdi)); \
	} \
	rsi += -8*clmblkdst+clmblkdst8; \
} \
rsi += -8*clmblkdst8-24*clmblkdst8;/* Restore src pointer */ \
end_timer(23); \
\
/* Do FFT levels 6,7 */ \
/* Values 0K-127K is real data, 128K-255K is semi-real data, 256K-3583K is */ \
/* complex data. */ \
/* */ \
/* On input the 64-byte cache lines hold these data values: */ \
/*	0K	128K	256K	384K */ \
/*	2K	... */ \
/*	... */ \
/*	126K	... */ \
/*	512K	... */ \
/*	... */ \
/* On output the 64-byte cache lines hold these data values: */ \
/*	0K	32K	64K	96K */ \
/*	2K	... */ \
/*	... */ \
/*	30K	... */ \
/*	128K	... */ \
/*	... */ \
\
/* Do 32 eight_reals_fft_2 macros */ \
/* Do 32 two_two_complex_fft macros */ \
/*	distance between fft data elements is 32K */ \
/*	do 16 macros each processing 32 data values */ \
\
start_timer(24); \
rdi = (uintptr_t)g->sincos3;/* Load sin/cos pointer */ \
/* 2 iterations of 4 iters of 2*2*clm */ \
for(unsigned int loopA = 2; loopA; loopA--) { \
	for(unsigned int loopB = 4; loopB; loopB--) { \
		xtouch(xptr(rcx+(clm-1)*128));	/* Preload TLB */ \
		xloop2_clm##clm(b4b, clmblkdst, x4cl_eight_reals_fft_2(rsi, 64, 2*clmblkdst8, 4*clmblkdst8, rdi)); \
	} \
	rsi += -8*clmblkdst+clmblkdst8; \
} \
\
/* Do 384 four_complex_fft macros */ \
/*	distance between fft data elements is 8K */ \
/*	do 96 macros each processing 32 data values */ \
\
rsi += -2*clmblkdst8+8*clmblkdst8;/* Load source pointer */ \
rdi += XMM_SCD;		/* Next sine/cosine pointer */ \
/* 6 iters of 2 iters of 4 of 2*2*clm */ \
for(unsigned int loopA = 6; loopA; loopA--) { \
	for(unsigned int loopB = 2; loopB; loopB--) { \
		for(unsigned int loopC = 4; loopC; loopC--) { \
			xtouch(xptr(rcx+(clm-1)*128));	/* Preload TLB */ \
			xloop2_clm##clm(b5b, clmblkdst, x4cl_four_complex_fft(rsi, 64, 2*clmblkdst8, 4*clmblkdst8, rdi)); \
		} \
		rsi += -8*clmblkdst+clmblkdst8; \
	} \
	rsi += -2*clmblkdst8+8*clmblkdst8;/* Next source pointer */ \
	rdi += 2*XMM_SCD;		/* Next sine/cosine pointer */ \
} \
rsi += -7*8*clmblkdst8;	/* Restore source pointer */ \
end_timer(24); \
\
/* Do FFT levels 8,9 */ \
/* Values 0K-31K is real data, 32K-63K is semi-real data, 64K-3583K is */ \
/* complex data. */ \
/* */ \
/* On input the 64-byte cache lines hold these data values: */ \
/*	0K	32K	64K	96K */ \
/*	2K	... */ \
/*	... */ \
/*	30K	... */ \
/*	128K	... */ \
/*	... */ \
/* On output the 64-byte cache lines hold these data values: */ \
/*	0K	8K	16K	24K */ \
/*	2K	... */ \
/*	... */ \
/*	6K	... */ \
/*	32K	... */ \
/*	... */ \
\
/* Do 8 eight_reals_fft_2 macros */ \
/* Do 8 two_two_complex_fft macros */ \
/*	distance between fft data elements is 8K */ \
/*	do 4 macros each processing 32 data values */ \
\
start_timer(25); \
rdi = (uintptr_t)g->sincos2;/* Load sin/cos pointer */ \
/* 4 iterations of 2*clm */ \
for(unsigned int loopA = 4; loopA; loopA--) { \
	for(unsigned int loopB = 2*clm; loopB; loopB--) { \
		x4cl_eight_reals_fft_2(rsi, 64, 4*clmblkdst, clmblkdst8, rdi); \
	} \
	rsi += -2*clm*64+clmblkdst;/* Next source pointer */ \
} \
\
/* Do 432 four_complex_fft macros */ \
/*	distance between fft data elements is 8K */ \
/*	do 108 macros each processing 32 data values */ \
\
rsi += -4*clmblkdst+2*clmblkdst8;/* Load source pointer */ \
rdi += XMM_SCD;		/* Next sine/cosine pointer */ \
/* 27 iterations of 2 iters of 2*2*clm */ \
for(unsigned int loopA = 27; loopA; loopA--) { \
	for(unsigned int loopB = 2; loopB; loopB--) { \
		xtouch(xptr(rcx+(clm-1)*128));	/* Preload TLB */ \
		xloop2_clm##clm(b7b, clmblkdst, x4cl_four_complex_fft(rsi, 64, 4*clmblkdst, clmblkdst8, rdi)); \
	} \
	rsi += -2*2*clmblkdst+2*clmblkdst8;/* Next source pointer */ \
	rdi += 2*XMM_SCD;		/* Next sine/cosine pointer */ \
} \
rsi += -28*2*clmblkdst8;	/* Restore source pointer */ \
pfing(rcx += -448*blkdst);	/* Restore prefetching pointer */ \
end_timer(25); \
\
/* Do FFT levels 10,11 */ \
/* Values 0K-7K is real data, 8K-15K is semi-real data, 16K-3583K is */ \
/* complex data. */ \
/* */ \
/* On input the 64-byte cache lines hold these data values: */ \
/*	0K	8K	16K	24K */ \
/*	2K	... */ \
/*	... */ \
/*	6K	... */ \
/*	32K	... */ \
/*	... */ \
/* On output the 64-byte cache lines hold these data values: */ \
/*	0K	2K	4K	6K */ \
/*	8K	... */ \
/*	... */ \
\
/* Do 1 eight_reals_fft_2 macros */ \
/* Do 1 nop_two_two_complex_fft_2 macros */ \
/* Do 2 four_complex_fft macros */ \
/*	distance between fft data elements is 2K */ \
/*	do 1 macros each processing 32 data values */ \
\
start_timer(26); \
rdi = (uintptr_t)g->sincos1;/* Load sin/cos pointer */ \
for(unsigned int loopA = 2*clm; loopA; loopA--) { \
	g4cl_eight_reals_fft_2(rsi, 64, clmblkdst, 2*clmblkdst, rcx, 64, blkdst, 2*blkdst); \
} \
\
/* Do 444 four_complex_fft macros */ \
/*	distance between fft data elements is 2K */ \
/*	do 111 macros each processing 32 data values */ \
\
rsi += -2*clm*64+4*clmblkdst;/* Next source pointer */ \
rcx += -2*clm*64+4*blkdst;	/* Next dest pointer */ \
rdi += XMM_SCD;		/* Next sine/cosine pointer */ \
/* eax = 56*256+1; */ /* 111 iterations of 2*clm */ \
for(unsigned int loopA = 111; loopA; ) { \
	do{ \
		for(unsigned int loopC = clm; loopC; loopC--) { \
			for(unsigned int loopD = 2; loopD; loopD--) { \
				g4cl_four_complex_fft(rsi, 64, clmblkdst, 2*clmblkdst, rcx, 64, blkdst, 2*blkdst); \
				L2prefetch128(u8ptr(rcx-64+rbp)); \
				L2prefetch128(u8ptr(rcx-64+blkdst+rbp)); \
				pfing(rbp += -64+2*blkdst); \
			} \
			pfing(rbp += 128-4*blkdst); \
		} \
		rsi += -2*clm*64+4*clmblkdst;/* Next source pointer */ \
		rcx += -2*clm*64+4*blkdst;	/* Next dest pointer */ \
		rdi += 2*XMM_SCD;		/* Next sine/cosine pointer */ \
	}while(--loopA & 1); \
	rsi += -8*clmblkdst+clmblkdst8; \
} \
end_timer(26); \
\
/* Work on next 4 sets of 1792 values */ \
\
endlp:	xpass1_get_next_block(pass2, c0b, b0b); /* Get next block and jump */ \
/* done: */	xfft_footer_hg2; \
\
\
\
/* Perform a 4M-element FFT.  This is done in two passes.  Pass 1 */ \
/* does 11 levels, pass 2 uses common code to perform the last 11 levels. */ \
\

#define hg_pass1sclevels11pfa8(procname, xpass2name, pass2_levels, clm) \
 /* (fftname, pass2_levels, clmarg) */ \
\
/* set_FFT_constants(fftname, 0, 0, pass2_levels, clmarg); */ \
xfft_header_hg2(procname, xpass2name, pass2_levels, clm); \
\
/* Do the last 11 FFT levels */ \
/*	do 1 pass2_11_levels_real and 511 pass2_11_levels_complex */ \
\
/* pass2: */	start_timer(0); \
xpass2name(g); \
end_timer(0); \
if(g->ffttype == 1)		/* We're done if FFTing only */ \
	return; \
\
/* Do the inverse FFT */ \
\
pass1_inverse_fft_setup(c0b, procname); \
\
/* Do inverse FFT levels 10,11 */ \
/* On input the 64-byte cache lines hold these data values: */ \
/*	0K	2K	4K	6K */ \
/*	8K	10K	12K	14K */ \
/*	... */ \
/* On output the 64-byte cache lines hold these data values: */ \
/*	0K	4K	8K	12K */ \
/*	2K	... */ \
/*	16K	... */ \
/*	... */ \
\
/* Do 1 eight_reals_unfft_2 macros */ \
/* Do 1 nop_two_two_complex_unfft_2 macros */ \
/*	distance between fft data elements is 2K */ \
/*	do 1 macros each processing 16 data values */ \
\
c0b:	start_timer(16); \
rsi = (uintptr_t)g->data_addr;/* Get source address */ \
rcx = (uintptr_t)g->scratch_area;	/* Get address of scratch area */ \
pfing(rbp = (uintptr_t)g->data_prefetch);/* Set prefetch distance */ \
pfing(rbp -= rsi); \
for(unsigned int loopA = 2*clm; loopA; loopA--) { \
	g2cl_eight_reals_unfft_2(rsi, 64, blkdst, rcx, 64, clmblkdst); \
} \
\
/* Do 510 four_complex_unfft macros */ \
/*	distance between fft data elements is 2K */ \
/*	do 255 macros each processing 16 data values */ \
\
rsi += -2*clm*64+2*blkdst;	/* Next source pointer */ \
rcx += -2*clm*64+2*clmblkdst;/* Next dest pointer */ \
rdi = (uintptr_t)g->sincos1;/* Load sin/cos pointer */ \
/* eax = 64*256+256/4; */ /* 255 iterations of 2*clm */ \
for(unsigned int loopA = 255; loopA; ) { \
	do{ \
		for(unsigned int loopC = 2*clm; loopC; loopC--) { \
			g2cl_four_complex_unfft(rsi, 64, blkdst, rcx, 64, clmblkdst); \
		} \
		rsi += -2*clm*64+2*blkdst;	/* Next source pointer */ \
		rcx += -2*clm*64+2*clmblkdst;/* Next dest pointer */ \
		rdi += XMM_SCD;		/* Next sine/cosine pointer */ \
	}while(--loopA & 3); \
	rcx += -8*clmblkdst+clmblkdst8; \
} \
rcx = rsi-256*2*blkdst;	/* Create prefetch pointer */ \
rsi = (uintptr_t)g->scratch_area;	/* Restore source pointer */ \
end_timer(16); \
\
/* Do inverse FFT levels 8,9 */ \
/* On input the 64-byte cache lines hold these data values: */ \
/*	0K	4K	8K	12K */ \
/*	2K	... */ \
/*	16K	... */ \
/*	... */ \
/* On output the 64-byte cache lines hold these data values: */ \
/*	0K	16K	32K	48K */ \
/*	2K	... */ \
/*	... */ \
/*	14K */ \
/*	64K	... */ \
/*	... */ \
\
/* Do 4 eight_reals_unfft_2 macros */ \
/* Do 4 nop_two_two_complex_unfft_2 macros */ \
/*	distance between fft data elements is 8K */ \
/*	do 2 macros each processing 32 data values */ \
\
start_timer(17); \
xtouch(xptr(rcx+(clm-1)*128));	/* Preload TLB */ \
xloop2_clm##clm(c6b, clmblkdst, x4cl_eight_reals_unfft_2(rsi, 64, 2*clmblkdst, 4*clmblkdst));		/* 2*2*clm iterations */ \
\
/* Do 504 four_complex_unfft macros */ \
/*	distance between fft data elements is 8K */ \
/*	do 126 macros each processing 32 data values */ \
\
rsi += -2*clmblkdst+clmblkdst8;/* Next source pointer */ \
rdi = (uintptr_t)g->sincos2;/* Load sin/cos pointer */ \
/* 63 iterations of 2*2*clm */ \
for(unsigned int loopA = 63; loopA; loopA--) { \
	xtouch(xptr(rcx+(clm-1)*128));	/* Preload TLB */ \
	xloop2_clm##clm(c7b, clmblkdst, x4cl_four_complex_unfft(rsi, 64, 2*clmblkdst, 4*clmblkdst, rdi)); \
	rsi += -2*clmblkdst+clmblkdst8;/* Next source pointer */ \
	rdi += XMM_SCD;		/* Next sine/cosine pointer */ \
} \
rsi += -64*clmblkdst8;	/* Restore source pointer */ \
end_timer(17); \
\
/* Do inverse FFT levels 6,7 */ \
/* On input the 64-byte cache lines hold these data values: */ \
/*	0K	16K	32K	48K */ \
/*	2K	... */ \
/*	... */ \
/*	14K */ \
/*	64K	... */ \
/*	... */ \
/* On output the 64-byte cache lines hold these data values: */ \
/*	0K	64K	128K	192K */ \
/*	2K	... */ \
/*	... */ \
/*	62K */ \
/*	256K	... */ \
/*	... */ \
\
/* Do 16 eight_reals_unfft_2 macros */ \
/* Do 16 nop_two_two_complex_unfft_2 macros */ \
/*	distance between fft data elements is 32K */ \
/*	do 8 macros each processing 32 data values */ \
\
start_timer(18); \
/* 4 iterations of 2*2*clm */ \
for(unsigned int loopA = 4; loopA; loopA--) { \
	xtouch(xptr(rcx+(clm-1)*128));	/* Preload TLB */ \
	xloop2_clm##clm(c4b, clmblkdst, x4cl_eight_reals_unfft_2(rsi, 64, clmblkdst8, 2*clmblkdst8)); \
} \
\
/* Do 480 four_complex_unfft macros */ \
/*	distance between fft data elements is 32K */ \
/*	do 120 macros each processing 32 data values */ \
\
rsi += -8*clmblkdst+4*clmblkdst8;/* Next source pointer */ \
rdi = (uintptr_t)g->sincos3;/* Load sin/cos pointer */ \
/* 15 iterations of 4 iters of 2*2*clm */ \
for(unsigned int loopA = 15; loopA; loopA--) { \
	for(unsigned int loopB = 4; loopB; loopB--) { \
		xtouch(xptr(rcx+(clm-1)*128));	/* Preload TLB */ \
		xloop2_clm##clm(c5b, clmblkdst, x4cl_four_complex_unfft(rsi, 64, clmblkdst8, 2*clmblkdst8, rdi)); \
	} \
	rsi += -4*2*clmblkdst+4*clmblkdst8;/* Next source pointer */ \
	rdi += XMM_SCD;		/* Next sine/cosine pointer */ \
} \
rsi += -16*4*clmblkdst8;	/* Restore source pointer */ \
end_timer(18); \
\
/* Do inverse FFT levels 4,5 */ \
/* On input the 64-byte cache lines hold these data values: */ \
/*	0K	64K	128K	192K */ \
/*	2K	... */ \
/*	... */ \
/*	62K */ \
/*	256K	... */ \
/*	... */ \
/* On output the 64-byte cache lines hold these data values: */ \
/*	0K	256K	512K	768K */ \
/*	2K	... */ \
/*	... */ \
/*	254K */ \
/*	1024K	... */ \
/*	... */ \
/* */ \
/* Do 64 eight_reals_unfft_2 macros */ \
/* Do 64 nop_two_two_complex_unfft_2 macros */ \
/*	distance between fft data elements is 128K */ \
/*	do 32 macros each processing 32 data values */ \
/* */ \
/* Do 384 four_complex_unfft macros */ \
/*	distance between fft data elements is 128K */ \
/*	do 96 macros each processing 32 data values */ \
\
/* Do inverse FFT levels 1,2,3 */ \
/* On input the 64-byte cache lines hold these data values: */ \
/*	0K	256K	512K	768K */ \
/*	2K	... */ \
/*	... */ \
/*	254K */ \
/*	1024K	... */ \
/*	... */ \
/* On output the 64-byte cache lines hold these data values: */ \
/*	0K	1024K	2048K	3072K */ \
/*	2K	... */ \
/*	... */ \
/*	1022K */ \
/* */ \
/* Do 512 eight_reals_last_unfft macros */ \
/*	distance between fft data elements is 512K */ \
/*	do 128 macros each processing 32 data values */ \
\
start_timer(19); \
rdi = (uintptr_t)g->sincos4;/* Load sin/cos pointer */ \
/* 32 iterations */ \
\
for(unsigned int loopA = 4; loopA; loopA--) { \
	for(unsigned int loopB = 8; loopB; loopB--) { \
		xtouchlines(4);			/* Preload TLBs */ \
		\
		for(unsigned int loopC = clm; loopC; loopC--) { \
			for(unsigned int loopD = 2; loopD; loopD--) { \
				x4cl_eight_reals_unfft_2(rsi, 64, 4*clmblkdst8, 8*clmblkdst8); \
			} \
			L2prefetch128(u8ptr(rcx+rbp)); \
			pfing(rcx += blkdst);		/* Next prefetch pointer */ \
			\
			rsi += -2*64+16*clmblkdst8;/* Next source pointer */ \
			/* 3 iterations */ \
			for(unsigned int loopD = 3; loopD; loopD--) { \
				for(unsigned int loopE = 2; loopE; loopE--) { \
					x4cl_four_complex_unfft(rsi, 64, 4*clmblkdst8, 8*clmblkdst8, rdi); \
				} \
				rsi += -2*64+16*clmblkdst8;/* Next source pointer */ \
				rdi += XMM_SCD;		/* Next sine/cosine pointer */ \
			} \
			L2prefetch128(u8ptr(rcx+rbp)); \
			pfing(rcx += blkdst);		/* Next prefetch pointer */ \
			rdi += -3*XMM_SCD;		/* Restore sine/cosine pointer */ \
			\
			rsi += -4*16*clmblkdst8;	/* Restore source pointer */ \
			for(unsigned int loopD = 2; loopD; loopD--) { \
				for(unsigned int loopE = 2; loopE; loopE--) { \
					for(unsigned int loopF = 2; loopF; loopF--) { \
						x4cl_eight_reals_last_unfft(rsi, 64, 16*clmblkdst8, 32*clmblkdst8); \
					} \
					rsi += -2*64+4*clmblkdst8; \
				} \
				L2prefetch128(u8ptr(rcx+rbp)); \
				pfing(rcx += blkdst);		/* Next prefetch pointer */ \
			} \
			\
				EXEC248##clm(rsi += -4*4*clmblkdst8+2*64);/* Next source pointer */ \
			EXEC248##clm(pfing(rcx += -4*blkdst+128));	/* Next prefetch pointer */ \
		} \
		EXEC248##clm(rsi += -clm*2*64+clmblkdst);/* Next source pointer */ \
		EXEC248##clm(pfing(rcx += -clm*128+4*blkdst));	/* Next prefetch pointer */ \
		EXEC1##clm(rsi += -4*4*clmblkdst8+clmblkdst);/* Next source pointer */ \
		\
	} \
	rsi += -8*clmblkdst+clmblkdst8;/* Next source pointer */ \
} \
end_timer(19); \
\
/* Normalize these values */ \
\
pass1_normalize(1, yes_fft); \
\
/* Copy data from the scratch area back to the FFT data area */ \
\
copy_scratch_data_to_fft(512, clm); \
goto	endlp;			/* Skip FFT code if flag not set */ \
\
/* Fire up auxiliary threads to do pass 1 forward FFT */ \
\
/* pass1: */	pass1_forward_fft_setup(b0b); \
\
/* Copy data from the FFT data area back to the scratch area */ \
\
b0b:	copy_fft_data_to_scratch(512, clm); \
\
/* Do FFT levels 1,2,3 */ \
/* */ \
/* On input the 64-byte cache lines hold these data values: */ \
/*	0K	1024K	2048K	3072K */ \
/*	2K	... */ \
/*	... */ \
/*	1022K */ \
/* On output the 64-byte cache lines hold these data values: */ \
/*	0K	512K	1024K	1536K */ \
/*	2K	... */ \
/*	... */ \
/*	510K	... */ \
/*	2048K	... */ \
/*	... */ \
/* */ \
/* Do 512 eight_reals_first_fft macros (FFT levels 1,2,3) */ \
/*	distance between fft data elements is 512K */ \
/*	do 256 iterations each processing 16 data values */ \
\
yes_fft:start_timer(22); \
rsi = (uintptr_t)g->scratch_area;	/* Get address of scratch area */ \
rcx = (uintptr_t)g->data_addr; /* Load source address for prefetching */ \
pfing(rbp = (uintptr_t)g->data_prefetch);/* Set prefetch distance */ \
pfing(rbp -= rcx); \
pfing(rcx += 256*blkdst);		/* Create prefetching pointer */ \
/* 32 iterations of 2 iters of 4*2*clm */ \
for(unsigned int loopA = 32; loopA; loopA--) { \
	for(unsigned int loopB = 2; loopB; loopB--) { \
		xtouch(xptr(rcx+(clm-1)*128));	/* Preload TLB */ \
		xloop4_clm##clm(b1b, clmblkdst, x2cl_eight_reals_fft(rsi, 64, 32*clmblkdst8)); \
	} \
	rsi += -8*clmblkdst+clmblkdst8; \
} \
rsi += -32*clmblkdst8;	/* Restore source pointer */ \
end_timer(22); \
\
/* Do FFT levels 4,5 */ \
/* Values 0K-511K is real data, 512K-1023K is semi-real data, 1024K-4095K is */ \
/* complex data. */ \
/* */ \
/* On input the 64-byte cache lines hold these data values: */ \
/*	0K	512K	1024K	1536K */ \
/*	2K	... */ \
/*	... */ \
/*	510K	... */ \
/*	2048K	... */ \
/*	... */ \
/* On output the 64-byte cache lines hold these data values: */ \
/*	0K	128K	256K	384K */ \
/*	2K	... */ \
/*	... */ \
/*	126K	... */ \
/*	512K	... */ \
/*	... */ \
\
/* Do 64 eight_reals_fft_2 macros */ \
/* Do 64 nop_two_two_complex_fft_2 macros */ \
/* Do 128 four_complex_fft macros */ \
/*	distance between fft data elements is 128K */ \
/*	do 64 macros each processing 32 data values */ \
\
start_timer(23); \
rdi = (uintptr_t)g->sincos4;/* Load sin/cos pointer */ \
/* 8 iterations of 4 iters of 2*2*clm */ \
for(unsigned int loopA = 8; loopA; loopA--) { \
	for(unsigned int loopB = 4; loopB; loopB--) { \
		xtouch(xptr(rcx+(clm-1)*128));	/* Preload TLB */ \
		xloop2_clm##clm(b2b, clmblkdst, x4cl_eight_reals_fft_2(rsi, 64, 8*clmblkdst8, 16*clmblkdst8, rdi)); \
	} \
	rsi += -8*clmblkdst+clmblkdst8; \
} \
\
/* Do 256 four_complex_fft macros */ \
/*	distance between fft data elements is 128K */ \
/*	do 64 macros each processing 32 data values */ \
\
rsi += -8*clmblkdst8+32*clmblkdst8;/* Next source pointer */ \
rdi += XMM_SCD;		/* Next sine/cosine pointer */ \
/* 8 iterations of 4 iters of 2*2*clm */ \
for(unsigned int loopA = 8; loopA; loopA--) { \
	for(unsigned int loopB = 4; loopB; loopB--) { \
		xtouch(xptr(rcx+(clm-1)*128));	/* Preload TLB */ \
		xloop2_clm##clm(b3b, clmblkdst, x4cl_four_complex_fft(rsi, 64, 8*clmblkdst8, 16*clmblkdst8, rdi)); \
	} \
	rsi += -8*clmblkdst+clmblkdst8; \
} \
rsi += -8*clmblkdst8-32*clmblkdst8;/* Restore source pointer */ \
end_timer(23); \
\
/* Do FFT levels 6,7 */ \
/* Values 0K-127K is real data, 128K-255K is semi-real data, 256K-4095K is */ \
/* complex data. */ \
/* */ \
/* On input the 64-byte cache lines hold these data values: */ \
/*	0K	128K	256K	384K */ \
/*	2K	... */ \
/*	... */ \
/*	126K	... */ \
/*	512K	... */ \
/*	... */ \
/* On output the 64-byte cache lines hold these data values: */ \
/*	0K	32K	64K	96K */ \
/*	2K	... */ \
/*	... */ \
/*	30K	... */ \
/*	128K	... */ \
/*	... */ \
\
/* Do 32 eight_reals_fft_2 macros */ \
/* Do 32 two_two_complex_fft macros */ \
/*	distance between fft data elements is 32K */ \
/*	do 16 macros each processing 32 data values */ \
\
start_timer(24); \
rdi = (uintptr_t)g->sincos3;/* Load sin/cos pointer */ \
/* 2 iterations of 4 iters of 2*2*clm */ \
for(unsigned int loopA = 2; loopA; loopA--) { \
	for(unsigned int loopB = 4; loopB; loopB--) { \
		xtouch(xptr(rcx+(clm-1)*128));	/* Preload TLB */ \
		xloop2_clm##clm(b4b, clmblkdst, x4cl_eight_reals_fft_2(rsi, 64, 2*clmblkdst8, 4*clmblkdst8, rdi)); \
	} \
	rsi += -8*clmblkdst+clmblkdst8; \
} \
\
/* Do 448 four_complex_fft macros */ \
/*	distance between fft data elements is 8K */ \
/*	do 112 macros each processing 32 data values */ \
\
rsi += -2*clmblkdst8+8*clmblkdst8;/* Load source pointer */ \
rdi += XMM_SCD;		/* Next sine/cosine pointer */ \
/* 7 iters of 2 iters of 4 of 2*2*clm */ \
for(unsigned int loopA = 7; loopA; loopA--) { \
	for(unsigned int loopB = 2; loopB; loopB--) { \
		for(unsigned int loopC = 4; loopC; loopC--) { \
			xtouch(xptr(rcx+(clm-1)*128));	/* Preload TLB */ \
			xloop2_clm##clm(b5b, clmblkdst, x4cl_four_complex_fft(rsi, 64, 2*clmblkdst8, 4*clmblkdst8, rdi)); \
		} \
		rsi += -8*clmblkdst+clmblkdst8; \
	} \
	rsi += -2*clmblkdst8+8*clmblkdst8;/* Next source pointer */ \
	rdi += 2*XMM_SCD;		/* Next sine/cosine pointer */ \
} \
rsi += -8*8*clmblkdst8;	/* Restore source pointer */ \
end_timer(24); \
\
/* Do FFT levels 8,9 */ \
/* Values 0K-31K is real data, 32K-63K is semi-real data, 64K-4095K is */ \
/* complex data. */ \
/* */ \
/* On input the 64-byte cache lines hold these data values: */ \
/*	0K	32K	64K	96K */ \
/*	2K	... */ \
/*	... */ \
/*	30K	... */ \
/*	128K	... */ \
/*	... */ \
/* On output the 64-byte cache lines hold these data values: */ \
/*	0K	8K	16K	24K */ \
/*	2K	... */ \
/*	... */ \
/*	6K	... */ \
/*	32K	... */ \
/*	... */ \
\
/* Do 8 eight_reals_fft_2 macros */ \
/* Do 8 two_two_complex_fft macros */ \
/*	distance between fft data elements is 8K */ \
/*	do 4 macros each processing 32 data values */ \
\
start_timer(25); \
rdi = (uintptr_t)g->sincos2;/* Load sin/cos pointer */ \
/* 2 iterations of 2*2*clm */ \
for(unsigned int loopA = 2; loopA; loopA--) { \
	xtouch(xptr(rcx+(clm-1)*128));	/* Preload TLB */ \
	xloop2_clm##clm(b6b, clmblkdst, x4cl_eight_reals_fft_2(rsi, 64, 4*clmblkdst, clmblkdst8, rdi)); \
} \
\
/* Do 496 four_complex_fft macros */ \
/*	distance between fft data elements is 8K */ \
/*	do 124 macros each processing 32 data values */ \
\
rsi += -2*2*clmblkdst+2*clmblkdst8;/* Load source pointer */ \
rdi += XMM_SCD;		/* Next sine/cosine pointer */ \
/* 31 iterations of 2 iters of 2*2*clm */ \
for(unsigned int loopA = 31; loopA; loopA--) { \
	for(unsigned int loopB = 2; loopB; loopB--) { \
		xtouch(xptr(rcx+(clm-1)*128));	/* Preload TLB */ \
		xloop2_clm##clm(b7b, clmblkdst, x4cl_four_complex_fft(rsi, 64, 4*clmblkdst, clmblkdst8, rdi)); \
	} \
	rsi += -2*2*clmblkdst+2*clmblkdst8;/* Next source pointer */ \
	rdi += 2*XMM_SCD;		/* Next sine/cosine pointer */ \
} \
rsi += -32*2*clmblkdst8;	/* Restore source pointer */ \
pfing(rcx += -512*blkdst);	/* Restore prefetch pointer */ \
end_timer(25); \
\
/* Do FFT levels 10,11 */ \
/* Values 0K-7K is real data, 8K-15K is semi-real data, 16K-4095K is */ \
/* complex data. */ \
/* */ \
/* On input the 64-byte cache lines hold these data values: */ \
/*	0K	8K	16K	24K */ \
/*	2K	... */ \
/*	... */ \
/*	6K	... */ \
/*	32K	... */ \
/*	... */ \
/* On output the 64-byte cache lines hold these data values: */ \
/*	0K	2K	4K	6K */ \
/*	8K	... */ \
/*	... */ \
\
/* Do 1 eight_reals_fft_2 macros */ \
/* Do 1 nop_two_two_complex_fft_2 macros */ \
/* Do 2 four_complex_fft macros */ \
/*	distance between fft data elements is 2K */ \
/*	do 1 macros each processing 32 data values */ \
\
start_timer(26); \
rdi = (uintptr_t)g->sincos1;/* Load sin/cos pointer */ \
for(unsigned int loopA = 2*clm; loopA; loopA--) { \
	g4cl_eight_reals_fft_2(rsi, 64, clmblkdst, 2*clmblkdst, rcx, 64, blkdst, 2*blkdst); \
} \
\
/* Do 508 four_complex_fft macros */ \
/*	distance between fft data elements is 2K */ \
/*	do 127 macros each processing 32 data values */ \
\
rsi += -2*clm*64+4*clmblkdst;/* Next source pointer */ \
rcx += -2*clm*64+4*blkdst;	/* Next dest pointer */ \
rdi += XMM_SCD;		/* Next sine/cosine pointer */ \
/* eax = 64*256+1; */ /* 127 iterations of 2*clm */ \
for(unsigned int loopA = 127; loopA; ) { \
	do{ \
		for(unsigned int loopC = clm; loopC; loopC--) { \
			for(unsigned int loopD = 2; loopD; loopD--) { \
				g4cl_four_complex_fft(rsi, 64, clmblkdst, 2*clmblkdst, rcx, 64, blkdst, 2*blkdst); \
				L2prefetch128(u8ptr(rcx-64+rbp)); \
				L2prefetch128(u8ptr(rcx-64+blkdst+rbp)); \
				pfing(rbp += -64+2*blkdst); \
			} \
			pfing(rbp += 128-4*blkdst); \
		} \
		rsi += -2*clm*64+4*clmblkdst;/* Next source pointer */ \
		rcx += -2*clm*64+4*blkdst;	/* Next dest pointer */ \
		rdi += 2*XMM_SCD;		/* Next sine/cosine pointer */ \
	}while(--loopA & 1); \
	rsi += -8*clmblkdst+clmblkdst8; \
} \
end_timer(26); \
\
/* Work on next 4 sets of 2048 values */ \
\
endlp:	xpass1_get_next_block(pass2, c0b, b0b); /* Get next block and jump */ \
/* done: */	xfft_footer_hg2; \
\
\
\
/* Perform a 4M-element negacyclic FFT.  This is done in two passes.  Pass 1 */ \
/* does 11 levels, pass 2 uses common code to perform the last 11 levels. */ \
\

#define hg_pass1sclevels11complex4(procname, xpass2name, pass2_levels, clm) \
 /* (fftname, pass2_levels, clmarg) */ \
\
/* set_FFT_constants(fftname, 0, 1, pass2_levels, clmarg); */ \
xfft_header_hg2(procname, xpass2name, pass2_levels, clm); \
\
/* Do the last 11 FFT levels */ \
/*	do 512 pass2_11_levels_complex */ \
\
/* pass2: */	start_timer(0); \
xpass2name(g); \
end_timer(0); \
if(g->ffttype == 1)		/* We're done if FFTing only */ \
	return; \
\
/* Do the inverse FFT */ \
\
pass1_inverse_fft_setup(c0b, procname); \
\
/* Do inverse FFT levels 10,11 */ \
/* On input the 64-byte cache lines hold these data values: */ \
/*	0K	2K	4K	6K */ \
/*	8K	10K	12K	14K */ \
/*	... */ \
/* On output the 64-byte cache lines hold these data values: */ \
/*	0K	4K	8K	12K */ \
/*	2K	... */ \
/*	16K	... */ \
/*	... */ \
\
/* Do 512 four_complex_unfft macros */ \
/*	distance between fft data elements is 2K */ \
/*	do 256 macros each processing 16 data values */ \
\
c0b:	start_timer(16); \
rsi = (uintptr_t)g->data_addr;/* Get source address */ \
rcx = (uintptr_t)g->scratch_area;	/* Get address of scratch area */ \
rdi = (uintptr_t)g->sincos1;/* Load sin/cos pointer */ \
/* 256 iterations of 2*clm */ \
pfing(rbp = (uintptr_t)g->data_prefetch);/* Set prefetch distance */ \
pfing(rbp -= rsi); \
for(unsigned int loopA = 64; loopA; loopA--) { \
	for(unsigned int loopB = 4; loopB; loopB--) { \
		for(unsigned int loopC = 2*clm; loopC; loopC--) { \
			g2cl_four_complex_unfft(rsi, 64, blkdst, rcx, 64, clmblkdst); \
		} \
		rsi += -2*clm*64+2*blkdst;	/* Next source pointer */ \
		rcx += -2*clm*64+2*clmblkdst;/* Next dest pointer */ \
		rdi += XMM_SCD;		/* Next sine/cosine pointer */ \
	} \
	rcx += -8*clmblkdst+clmblkdst8; \
} \
rcx = rsi-256*2*blkdst;	/* Create prefetch pointer */ \
rsi = (uintptr_t)g->scratch_area;	/* Restore source pointer */ \
end_timer(16); \
\
/* Do inverse FFT levels 8,9 */ \
/* On input the 64-byte cache lines hold these data values: */ \
/*	0K	4K	8K	12K */ \
/*	2K	... */ \
/*	16K	... */ \
/*	... */ \
/* On output the 64-byte cache lines hold these data values: */ \
/*	0K	16K	32K	48K */ \
/*	2K	... */ \
/*	... */ \
/*	14K */ \
/*	64K	... */ \
/*	... */ \
\
/* Do 512 four_complex_unfft macros */ \
/*	distance between fft data elements is 8K */ \
/*	do 128 macros each processing 32 data values */ \
\
start_timer(17); \
rdi = (uintptr_t)g->sincos2;/* Load sin/cos pointer */ \
/* 64 iterations of 2*2*clm */ \
for(unsigned int loopA = 64; loopA; loopA--) { \
	xtouch(xptr(rcx+(clm-1)*128));	/* Preload TLB */ \
	xloop2_clm##clm(c7b, clmblkdst, x4cl_four_complex_unfft(rsi, 64, 2*clmblkdst, 4*clmblkdst, rdi));		/* 2 * 2*clm iterations */ \
	rsi += -2*clmblkdst+clmblkdst8;/* Next source pointer */ \
	rdi += XMM_SCD;		/* Next sine/cosine pointer */ \
} \
rsi += -64*clmblkdst8;	/* Restore source pointer */ \
end_timer(17); \
\
/* Do inverse FFT levels 6,7 */ \
/* On input the 64-byte cache lines hold these data values: */ \
/*	0K	16K	32K	48K */ \
/*	2K	... */ \
/*	... */ \
/*	14K */ \
/*	64K	... */ \
/*	... */ \
/* On output the 64-byte cache lines hold these data values: */ \
/*	0K	64K	128K	192K */ \
/*	2K	... */ \
/*	... */ \
/*	62K */ \
/*	256K	... */ \
/*	... */ \
\
/* Do 512 four_complex_unfft macros */ \
/*	distance between fft data elements is 8K */ \
/*	do 128 macros each processing 32 data values */ \
\
start_timer(18); \
rdi = (uintptr_t)g->sincos3;/* Load sin/cos pointer */ \
/* 16 iterations of 4 iters of 2*2*clm */ \
for(unsigned int loopA = 16; loopA; loopA--) { \
	for(unsigned int loopB = 4; loopB; loopB--) { \
		xtouch(xptr(rcx+(clm-1)*128));	/* Preload TLB */ \
		xloop2_clm##clm(c5b, clmblkdst, x4cl_four_complex_unfft(rsi, 64, clmblkdst8, 2*clmblkdst8, rdi));		/* 2 * 2*clm iterations */ \
	} \
	rsi += -4*2*clmblkdst+4*clmblkdst8;/* Next source pointer */ \
	rdi += XMM_SCD;		/* Next sine/cosine pointer */ \
} \
rsi += -16*4*clmblkdst8;	/* Restore source pointer */ \
end_timer(18); \
\
/* Do inverse FFT levels 4,5 */ \
/* On input the 64-byte cache lines hold these data values: */ \
/*	0K	64K	128K	192K */ \
/*	2K	... */ \
/*	... */ \
/*	62K */ \
/*	256K	... */ \
/*	... */ \
/* On output the 64-byte cache lines hold these data values: */ \
/*	0K	256K	512K	768K */ \
/*	2K	... */ \
/*	... */ \
/*	254K */ \
/*	1024K	... */ \
/*	... */ \
\
/* Do 512 four_complex_unfft macros */ \
/*	distance between fft data elements is 128K */ \
/*	do 128 macros each processing 32 data values */ \
\
start_timer(19); \
rdi = (uintptr_t)g->sincos4;/* Load sin/cos pointer */ \
/* 4 iters of 4 iters 4 of 2*2*clm */ \
for(unsigned int loopA = 4; loopA; loopA--) { \
	for(unsigned int loopB = 4; loopB; loopB--) { \
		for(unsigned int loopC = 4; loopC; loopC--) { \
			xtouch(xptr(rcx+(clm-1)*128));	/* Preload TLB */ \
			xloop2_clm##clm(c3b, clmblkdst, x4cl_four_complex_unfft(rsi, 64, 4*clmblkdst8, 8*clmblkdst8, rdi));		/* 2* 2*clm iterations */ \
		} \
		rsi += -4*2*clmblkdst+clmblkdst8;/* Next source pointer */ \
	} \
	rsi += -4*clmblkdst8+16*clmblkdst8;/* Next source pointer */ \
	rdi += XMM_SCD;		/* Next sine/cosine pointer */ \
} \
rsi += -4*16*clmblkdst8;	/* Restore source pointer */ \
end_timer(19); \
\
/* Do inverse FFT levels 1,2,3 */ \
/* On input the 64-byte cache lines hold these data values: */ \
/*	0K	256K	512K	768K */ \
/*	2K	... */ \
/*	... */ \
/*	254K */ \
/*	1024K	... */ \
/*	... */ \
/* On output the 64-byte cache lines hold these data values: */ \
/*	0K	1024K	2048K	3072K */ \
/*	2K	... */ \
/*	... */ \
/*	1022K */ \
\
/* Do 512 four_complex_last_unfft macros */ \
/*	distance between fft data elements is 512K */ \
/*	do 128 macros each processing 32 data values */ \
\
start_timer(20); \
rdi = (uintptr_t)g->plus1_premults;/* Address of premultiplier table */ \
/* 16 iters of 4 iters of 2*2*clm */ \
for(unsigned int loopA = 16; loopA; loopA--) { \
	for(unsigned int loopB = 4; loopB; loopB--) { \
		xtouch(xptr(rcx+(clm-1)*128));	/* Preload TLB */ \
		xloop2edi_clm##clm(c1b, clmblkdst, 128, x4cl_four_complex_last_unfft(rsi, 64, 16*clmblkdst8, 32*clmblkdst8, 128*128)); /* 2 * 2*clm iterations */ \
	} \
	rsi += -4*2*clmblkdst+clmblkdst8;/* Next source pointer */ \
} \
end_timer(20); \
\
/* Normalize these values */ \
\
pass1_normalize(1, yes_fft); \
\
/* Copy data from the scratch area back to the FFT data area */ \
\
copy_scratch_data_to_fft(512, clm); \
goto	endlp;			/* Skip FFT code if flag not set */ \
\
/* Fire up auxiliary threads to do pass 1 forward FFT */ \
\
/* pass1: */	pass1_forward_fft_setup(b0b); \
\
/* Copy data from the FFT data area back to the scratch area */ \
\
b0b:	copy_fft_data_to_scratch(512, clm); \
\
/* Do FFT levels 1,2,3 */ \
/* */ \
/* On input the 64-byte cache lines hold these data values: */ \
/*	0K	1024K	2048K	3072K */ \
/*	2K	... */ \
/*	... */ \
/*	1022K */ \
/* On output the 64-byte cache lines hold these data values: */ \
/*	0K	512K	1024K	1536K */ \
/*	2K	... */ \
/*	... */ \
/*	510K	... */ \
/*	2048K	... */ \
/*	... */ \
\
/* Do 512 four_complex_first_fft macros (FFT levels 1,2,3) */ \
/*	distance between fft data elements is 512K */ \
/*	do 256 iterations each processing 16 data values */ \
\
yes_fft:start_timer(22); \
rsi = (uintptr_t)g->scratch_area;	/* Get address of scratch area */ \
rdi = (uintptr_t)g->plus1_premults;/* Address of premultiplier table */ \
rcx = (uintptr_t)g->data_addr; /* Load source address for prefetching */ \
pfing(rbp = (uintptr_t)g->data_prefetch);/* Set prefetch distance */ \
pfing(rbp -= rcx); \
pfing(rcx += 256*blkdst);		/* Create prefetching pointer */ \
/* 32 iters of 2 iters of 4*2*clm */ \
for(unsigned int loopA = 32; loopA; loopA--) { \
	for(unsigned int loopB = 2; loopB; loopB--) { \
		xtouch(xptr(rcx+(clm-1)*128));	/* Preload TLB */ \
		xloop4edi_clm##clm(b1b, clmblkdst, 128, x2cl_four_complex_first_fft_scratch(rsi, 64, 32*clmblkdst8)); /* 4 * 2*clm iterations */ \
	} \
	rsi += -2*4*clmblkdst+clmblkdst8; \
} \
rsi += -32*clmblkdst8;	/* Restore source pointer */ \
end_timer(22); \
\
/* Do FFT levels 4,5 */ \
/* */ \
/* On input the 64-byte cache lines hold these data values: */ \
/*	0K	512K	1024K	1536K */ \
/*	2K	... */ \
/*	... */ \
/*	510K	... */ \
/*	2048K	... */ \
/*	... */ \
/* On output the 64-byte cache lines hold these data values: */ \
/*	0K	128K	256K	384K */ \
/*	2K	... */ \
/*	... */ \
/*	126K	... */ \
/*	512K	... */ \
/*	... */ \
\
/* Do 512 four_complex_fft macros */ \
/*	distance between fft data elements is 128K */ \
/*	do 128 macros each processing 32 data values */ \
\
start_timer(23); \
rdi = (uintptr_t)g->sincos4;/* Load sin/cos pointer */ \
/* 2 iters of 8 iters of 4 of 2*2*clm */ \
for(unsigned int loopA = 2; loopA; loopA--) { \
	for(unsigned int loopB = 8; loopB; loopB--) { \
		for(unsigned int loopC = 4; loopC; loopC--) { \
			xtouch(xptr(rcx+(clm-1)*128));	/* Preload TLB */ \
			xloop2_clm##clm(b3b, clmblkdst, x4cl_four_complex_fft(rsi, 64, 8*clmblkdst8, 16*clmblkdst8, rdi));		/* 2 * 2*clm iterations */ \
		} \
		rsi += -4*2*clmblkdst+clmblkdst8; \
	} \
	rsi += -8*clmblkdst8+32*clmblkdst8;/* Next source pointer */ \
	rdi += 2*XMM_SCD;		/* Next sine/cosine pointer */ \
} \
rsi += -2*32*clmblkdst8;	/* Restore source pointer */ \
end_timer(23); \
\
/* Do FFT levels 6,7 */ \
/* */ \
/* On input the 64-byte cache lines hold these data values: */ \
/*	0K	128K	256K	384K */ \
/*	2K	... */ \
/*	... */ \
/*	126K	... */ \
/*	512K	... */ \
/*	... */ \
/* On output the 64-byte cache lines hold these data values: */ \
/*	0K	32K	64K	96K */ \
/*	2K	... */ \
/*	... */ \
/*	30K	... */ \
/*	128K	... */ \
/*	... */ \
\
/* Do 4512 four_complex_fft macros */ \
/*	distance between fft data elements is 8K */ \
/*	do 128 macros each processing 32 data values */ \
\
start_timer(24); \
rdi = (uintptr_t)g->sincos3;/* Load sin/cos pointer */ \
/* 8 iters of 2 iters of 4 of 2*2*clm */ \
for(unsigned int loopA = 8; loopA; loopA--) { \
	for(unsigned int loopB = 2; loopB; loopB--) { \
		for(unsigned int loopC = 4; loopC; loopC--) { \
			xtouch(xptr(rcx+(clm-1)*128));	/* Preload TLB */ \
			xloop2_clm##clm(b5b, clmblkdst, x4cl_four_complex_fft(rsi, 64, 2*clmblkdst8, 4*clmblkdst8, rdi));		/* 2 * 2*clm iterations */ \
		} \
		rsi += -4*2*clmblkdst+clmblkdst8; \
	} \
	rsi += -2*clmblkdst8+8*clmblkdst8;/* Next source pointer */ \
	rdi += 2*XMM_SCD;		/* Next sine/cosine pointer */ \
} \
rsi += -8*8*clmblkdst8;	/* Restore source pointer */ \
end_timer(24); \
\
/* Do FFT levels 8,9 */ \
/* */ \
/* On input the 64-byte cache lines hold these data values: */ \
/*	0K	32K	64K	96K */ \
/*	2K	... */ \
/*	... */ \
/*	30K	... */ \
/*	128K	... */ \
/*	... */ \
/* On output the 64-byte cache lines hold these data values: */ \
/*	0K	8K	16K	24K */ \
/*	2K	... */ \
/*	... */ \
/*	6K	... */ \
/*	32K	... */ \
/*	... */ \
\
/* Do 512 four_complex_fft macros */ \
/*	distance between fft data elements is 8K */ \
/*	do 128 macros each processing 32 data values */ \
\
start_timer(25); \
rdi = (uintptr_t)g->sincos2;/* Load sin/cos pointer */ \
/* 32 iterations of 2 iters of 2*2*clm */ \
for(unsigned int loopA = 32; loopA; loopA--) { \
	for(unsigned int loopB = 2; loopB; loopB--) { \
		xtouch(xptr(rcx+(clm-1)*128));	/* Preload TLB */ \
		xloop2_clm##clm(b7b, clmblkdst, x4cl_four_complex_fft(rsi, 64, 4*clmblkdst, clmblkdst8, rdi));		/* 2 * 2*clm iterations */ \
	} \
	rsi += -2*2*clmblkdst+2*clmblkdst8;/* Next source pointer */ \
	rdi += 2*XMM_SCD;		/* Next sine/cosine pointer */ \
} \
rsi += -32*2*clmblkdst8;	/* Restore source pointer */ \
pfing(rcx += -512*blkdst);	/* Restore prefetch pointer */ \
end_timer(25); \
\
/* Do FFT levels 10,11 */ \
/* */ \
/* On input the 64-byte cache lines hold these data values: */ \
/*	0K	8K	16K	24K */ \
/*	2K	... */ \
/*	... */ \
/*	6K	... */ \
/*	32K	... */ \
/*	... */ \
/* On output the 64-byte cache lines hold these data values: */ \
/*	0K	2K	4K	6K */ \
/*	8K	... */ \
/*	... */ \
\
/* Do 512 four_complex_fft macros */ \
/*	distance between fft data elements is 2K */ \
/*	do 128 macros each processing 32 data values */ \
\
start_timer(26); \
rdi = (uintptr_t)g->sincos1;/* Load sin/cos pointer */ \
/* 64 iters of 2 of 2*clm */ \
for(unsigned int loopA = 64; loopA; loopA--) { \
	for(unsigned int loopB = 2; loopB; loopB--) { \
		for(unsigned int loopC = clm; loopC; loopC--) { \
			for(unsigned int loopD = 2; loopD; loopD--) { \
				g4cl_four_complex_fft(rsi, 64, clmblkdst, 2*clmblkdst, rcx, 64, blkdst, 2*blkdst); \
				L2prefetch128(u8ptr(rcx-64+rbp)); \
				L2prefetch128(u8ptr(rcx-64+blkdst+rbp)); \
				pfing(rbp += -64+2*blkdst); \
			} \
			pfing(rbp += 128-4*blkdst); \
		} \
	rsi += -2*clm*64+4*clmblkdst;/* Next source pointer */ \
		rcx += -2*clm*64+4*blkdst;	/* Next dest pointer */ \
		rdi += 2*XMM_SCD;		/* Next sine/cosine pointer */ \
	} \
	rsi += -8*clmblkdst+clmblkdst8; \
} \
end_timer(26); \
\
/* Work on next 4 sets of 2048 values */ \
\
endlp:	xpass1_get_next_block(pass2, c0b, b0b); /* Get next block and jump */ \
/* done: */	xfft_footer_hg2; \
\
\
\
/* Perform a 5120K-element FFT.  This is done in two passes.  Pass 1 */ \
/* does 12 levels, pass 2 uses common code to perform the last 11 levels. */ \
\

#define hg_pass1sclevels12pfa5(procname, xpass2name, pass2_levels, clm) \
 /* (fftname, pass2_levels, clmarg) */ \
\
/* set_FFT_constants(fftname, 0, 0, pass2_levels, clmarg); */ \
xfft_header_hg2(procname, xpass2name, pass2_levels, clm); \
\
/* Do the last 11 FFT levels */ \
/*	do 1 pass2_11_levels_real and 159 pass2_11_levels_complex */ \
\
/* pass2: */	start_timer(0); \
xpass2name(g); \
end_timer(0); \
if(g->ffttype == 1)		/* We're done if FFTing only */ \
	return; \
\
/* Do the inverse FFT */ \
\
pass1_inverse_fft_setup(c0b, procname); \
\
/* Do inverse FFT levels 11,12 */ \
/* On input the 64-byte cache lines hold these data values: */ \
/*	0K	2K	4K	6K */ \
/*	8K	10K	12K	14K */ \
/*	... */ \
/* On output the 64-byte cache lines hold these data values: */ \
/*	0K	4K	8K	12K */ \
/*	2K	... */ \
/*	16K	... */ \
/*	... */ \
\
/* Do 1 eight_reals_unfft_2 macros */ \
/* Do 1 nop_two_two_complex_unfft_2 macros */ \
/*	distance between fft data elements is 2K */ \
/*	do 1 macros each processing 16 data values */ \
\
c0b:	start_timer(16); \
rsi = (uintptr_t)g->data_addr;/* Get source address */ \
rcx = (uintptr_t)g->scratch_area;	/* Get address of scratch area */ \
pfing(rbp = (uintptr_t)g->data_prefetch);/* Set prefetch distance */ \
pfing(rbp -= rsi); \
for(unsigned int loopA = 2*clm; loopA; loopA--) { \
	g2cl_eight_reals_unfft_2(rsi, 64, blkdst, rcx, 64, clmblkdst); \
} \
\
/* Do 638 four_complex_unfft macros */ \
/*	distance between fft data elements is 2K */ \
/*	do 319 macros each processing 16 data values */ \
\
rsi += -2*clm*64+2*blkdst;	/* Next source pointer */ \
rcx += -2*clm*64+2*clmblkdst;/* Next dest pointer */ \
rdi = (uintptr_t)g->sincos1;/* Load sin/cos pointer */ \
/* eax = 80*256+256/4; */ /* 319 iterations of 2*clm */ \
for(unsigned int loopA = 319; loopA; ) { \
	do{ \
		for(unsigned int loopC = 2*clm; loopC; loopC--) { \
			g2cl_four_complex_unfft(rsi, 64, blkdst, rcx, 64, clmblkdst); \
		} \
		rsi += -2*clm*64+2*blkdst;	/* Next source pointer */ \
		rcx += -2*clm*64+2*clmblkdst;/* Next dest pointer */ \
		rdi += XMM_SCD;		/* Next sine/cosine pointer */ \
	}while(--loopA & 3); \
	rcx += -4*2*clmblkdst+clmblkdst8; \
} \
rcx = rsi-320*2*blkdst;	/* Create prefetch pointer */ \
rsi = (uintptr_t)g->scratch_area;	/* Restore source pointer */ \
end_timer(16); \
\
/* Do inverse FFT levels 9,10 */ \
/* On input the 64-byte cache lines hold these data values: */ \
/*	0K	4K	8K	12K */ \
/*	2K	... */ \
/*	16K	... */ \
/*	... */ \
/* On output the 64-byte cache lines hold these data values: */ \
/*	0K	16K	32K	48K */ \
/*	2K	... */ \
/*	... */ \
/*	14K */ \
/*	64K	... */ \
/*	... */ \
\
/* Do 4 eight_reals_unfft_2 macros */ \
/* Do 4 nop_two_two_complex_unfft_2 macros */ \
/*	distance between fft data elements is 8K */ \
/*	do 2 macros each processing 32 data values */ \
\
start_timer(17); \
/* 2 iterations of 2*clm */ \
for(unsigned int loopA = 2; loopA; loopA--) { \
	for(unsigned int loopB = 2*clm; loopB; loopB--) { \
		x4cl_eight_reals_unfft_2(rsi, 64, 2*clmblkdst, 4*clmblkdst); \
	} \
	rsi += -2*clm*64+clmblkdst;/* Next source pointer */ \
} \
\
/* Do 632 four_complex_unfft macros */ \
/*	distance between fft data elements is 8K */ \
/*	do 158 macros each processing 32 data values */ \
\
rsi += -2*clmblkdst+clmblkdst8;/* Next source pointer */ \
rdi = (uintptr_t)g->sincos2;/* Load sin/cos pointer */ \
/* 79 iterations of 2*2*clm */ \
for(unsigned int loopA = 79; loopA; loopA--) { \
	xtouch(xptr(rcx+(clm-1)*128));	/* Preload TLB */ \
	xloop2_clm##clm(c7b, clmblkdst, x4cl_four_complex_unfft(rsi, 64, 2*clmblkdst, 4*clmblkdst, rdi)); \
	rsi += -2*clmblkdst+clmblkdst8;/* Next source pointer */ \
	rdi += XMM_SCD;		/* Next sine/cosine pointer */ \
} \
rsi += -80*clmblkdst8;	/* Restore source pointer */ \
end_timer(17); \
\
/* Do inverse FFT levels 7,8 */ \
/* On input the 64-byte cache lines hold these data values: */ \
/*	0K	16K	32K	48K */ \
/*	2K	... */ \
/*	... */ \
/*	14K */ \
/*	64K	... */ \
/*	... */ \
/* On output the 64-byte cache lines hold these data values: */ \
/*	0K	64K	128K	192K */ \
/*	2K	... */ \
/*	... */ \
/*	62K */ \
/*	256K	... */ \
/*	... */ \
\
/* Do 16 eight_reals_unfft_2 macros */ \
/* Do 16 nop_two_two_complex_unfft_2 macros */ \
/*	distance between fft data elements is 32K */ \
/*	do 8 macros each processing 32 data values */ \
\
start_timer(18); \
/* 8 iterations of 2*clm */ \
for(unsigned int loopA = 8; loopA; loopA--) { \
	for(unsigned int loopB = 2*clm; loopB; loopB--) { \
		x4cl_eight_reals_unfft_2(rsi, 64, clmblkdst8, 2*clmblkdst8); \
	} \
	rsi += -2*clm*64+clmblkdst;/* Next source pointer */ \
} \
\
/* Do 608 four_complex_unfft macros */ \
/*	distance between fft data elements is 8K */ \
/*	do 152 macros each processing 32 data values */ \
\
rsi += -8*clmblkdst+4*clmblkdst8;/* Next source pointer */ \
rdi = (uintptr_t)g->sincos3;/* Load sin/cos pointer */ \
/* 19 iterations of 4 iters of 2*2*clm */ \
for(unsigned int loopA = 19; loopA; loopA--) { \
	for(unsigned int loopB = 4; loopB; loopB--) { \
		xtouch(xptr(rcx+(clm-1)*128));	/* Preload TLB */ \
		xloop2_clm##clm(c5d, clmblkdst, x4cl_four_complex_unfft(rsi, 64, clmblkdst8, 2*clmblkdst8, rdi)); \
	} \
	rsi += -4*2*clmblkdst+4*clmblkdst8;/* Next source pointer */ \
	rdi += XMM_SCD;		/* Next sine/cosine pointer */ \
} \
rsi += -20*4*clmblkdst8;	/* Restore source pointer */ \
end_timer(18); \
\
/* Do inverse FFT levels 6 */ \
/* On input the 64-byte cache lines hold these data values: */ \
/*	0K	64K	128K	192K */ \
/*	2K	... */ \
/*	... */ \
/*	62K */ \
/*	256K	... */ \
/*	... */ \
/* On output the 64-byte cache lines hold these data values: */ \
/*	0K	128K	256K	384K */ \
/*	2K	... */ \
/*	... */ \
/*	126K */ \
/*	512K	... */ \
/*	... */ \
/* */ \
/* Do 64 eight_reals_unfft_1 macros */ \
/*	distance between fft data elements is 128K */ \
/*	do 32 macros each processing 16 data values */ \
/* */ \
/* Do 576 four_complex_unfft macros */ \
/*	distance between fft data elements is 128K */ \
/*	do 288 macros each processing 16 data values */ \
\
/* Do inverse FFT levels 4,5 */ \
/* On input the 64-byte cache lines hold these data values: */ \
/*	0K	128K	256K	384K */ \
/*	2K	... */ \
/*	... */ \
/*	126K */ \
/*	512K	... */ \
/*	... */ \
/* On output the 64-byte cache lines hold these data values: */ \
/*	0K	256K	512K	768K */ \
/*	2K	... */ \
/*	... */ \
/*	254K */ \
/*	1024K	1536K	2048K	2304K */ \
/*	... */ \
/* */ \
/* Do 128 eight_reals_unfft_2 macros */ \
/*	distance between fft data elements is 256K */ \
/*	do 64 macros each processing 16 data values */ \
/* */ \
/* Do 512 four_complex_unfft macros */ \
/*	distance between fft data elements is 256K */ \
/*	do 128 macros each processing 32 data values */ \
\
/* Do inverse FFT levels 1,2,3 */ \
/* On input the 64-byte cache lines hold these data values: */ \
/*	0K	256K	512K	768K */ \
/*	2K	... */ \
/*	... */ \
/*	254K */ \
/*	1024K	1536K	2048K	2560K */ \
/*	... */ \
/* On output the 64-byte cache lines hold these data values: */ \
/*	0K	256K	2560K	2816K */ \
/*	2K	... */ \
/*	... */ \
/*	254K */ \
/*	512K	1536K	3072K	4096K */ \
/*	514K */ \
/*	... */ \
/*	1534K */ \
/* */ \
/* Do 1024 five_reals_last_unfft macros */ \
/*	distance between fft data elements is 1024K */ \
/*	do 128 macros each processing 40 data values */ \
\
start_timer(19); \
rdi = (uintptr_t)g->sincos4;/* Load sin/cos pointer */ \
rbx = (uintptr_t)g->sincos5;/* Load sin/cos pointer */ \
/* 4 iterations  */ \
\
for(unsigned int loopA = 4; loopA; loopA--) { \
	for(unsigned int loopB = 8; loopB; loopB--) { \
		xtouchlines(6);			/* Preload TLBs */ \
		\
		for(unsigned int loopC = clm; loopC; loopC--) { \
			for(unsigned int loopD = 2; loopD; loopD--) { \
				x2cl_eight_reals_unfft_1(rsi, 64, 4*clmblkdst8); \
			} \
			\
			rsi += -2*64+8*clmblkdst8;	/* Next source pointer */ \
			/* 9 iterations */ \
			for(unsigned int loopD = 9; loopD; loopD--) { \
				for(unsigned int loopE = 2; loopE; loopE--) { \
					x2cl_two_complex_unfft(rsi, 64, 4*clmblkdst8); \
				} \
				rsi += -2*64+8*clmblkdst8;	/* Next source pointer */ \
				rdi += XMM_SCD;		/* Next sine/cosine pointer */ \
			} \
			rdi += -9*XMM_SCD;		/* Restore sine/cosine pointer */ \
			\
			rsi += -10*8*clmblkdst8;	/* Restore source pointer */ \
			for(unsigned int loopD = 2; loopD; loopD--) { \
				for(unsigned int loopE = 2; loopE; loopE--) { \
					x2cl_half_eight_reals_unfft_2(rsi, 64, 8*clmblkdst8); \
				} \
				rsi += -2*64+4*clmblkdst8; \
			} \
			\
			rsi += -2*4*clmblkdst8+16*clmblkdst8;/* Next source pointer */ \
			for(unsigned int loopD = 2; loopD; loopD--) { \
				xloop2sp(c3b, 4*clmblkdst8, x4cl_four_complex_unfft(rsi, 64, 8*clmblkdst8, 16*clmblkdst8, rbx)); \
				rsi += -2*4*clmblkdst8+32*clmblkdst8; \
				rbx += XMM_SCD;		/* Next sine/cosine pointer */ \
			} \
			rbx += -2*XMM_SCD;		/* Restore sine/cosine pointer */ \
			\
			rsi += -2*32*clmblkdst8-16*clmblkdst8;/* Restore src pointer */ \
			xfive_reals_unfft_preload;	/* Preload consts for five real macros */ \
			for(unsigned int loopD = 4; loopD; loopD--) { \
				xloop1sp(c1b, 4*clmblkdst8, x5cl_five_reals_last_unfft(rsi, 64, 16*clmblkdst8)); \
			} \
			\
				EXEC248##clm(rsi += -4*4*clmblkdst8+2*64);/* Next source pointer */ \
			EXEC248##clm(pfing(rcx += -6*blkdst+128));	/* Next prefetch pointer */ \
		} \
		EXEC248##clm(rsi += -clm*2*64+clmblkdst);/* Next source pointer */ \
		EXEC248##clm(pfing(rcx += -clm*128+6*blkdst));	/* Next prefetch pointer */ \
		EXEC1##clm(rsi += -4*4*clmblkdst8+clmblkdst);/* Next source pointer */ \
		\
	} \
	rsi += -8*clmblkdst+clmblkdst8;/* Next source pointer */ \
} \
end_timer(19); \
\
/* Normalize these values */ \
\
pass1_normalize(1, yes_fft); \
\
/* Copy data from the scratch area back to the FFT data area */ \
\
copy_scratch_data_to_fft(640, clm); \
goto	endlp;			/* Skip FFT code if flag not set */ \
\
/* Fire up auxiliary threads to do pass 1 forward FFT */ \
\
/* pass1: */	pass1_forward_fft_setup(b0b); \
\
/* Copy data from the FFT data area back to the scratch area */ \
\
b0b:	copy_fft_data_to_scratch(640, clm); \
\
/* Do FFT levels 1,2,3 */ \
/* */ \
/* On input the 64-byte cache lines hold these data values: */ \
/*	0K	256K	2560K	2816K */ \
/*	2K	... */ \
/*	... */ \
/*	254K */ \
/*	512K	1536K	3072K	4096K */ \
/*	514K */ \
/*	... */ \
/*	1534K */ \
/* On output the 64-byte cache lines hold these data values: */ \
/*	0K	256K	512K	768K */ \
/*	2K	... */ \
/*	... */ \
/*	254K */ \
/*	1024K	2048K	3072K	4096K */ \
/*	1026K */ \
/*	... */ \
/*	2046K */ \
/* */ \
/* Do 1024 five_reals_first_fft macros (FFT levels 1,2,3) */ \
/*	distance between fft data elements is 1024K */ \
/*	do 128 iterations each processing 40 data values */ \
\
/* Do FFT levels 4,5 */ \
/* Values 0K-1023K is real data, 1024K-5119K is complex data. */ \
/* */ \
/* On input the 64-byte cache lines hold these data values: */ \
/*	0K	256K	512K	768K */ \
/*	2K	... */ \
/*	... */ \
/*	254K */ \
/*	1024K	2048K	3072K	4096K */ \
/*	1026K */ \
/*	... */ \
/*	2046K */ \
/* On output the 64-byte cache lines hold these data values: */ \
/*	0K	256K	512K	768K */ \
/*	2K	... */ \
/*	... */ \
/*	254K	... */ \
/*	1024K	... */ \
/*	... */ \
/* */ \
/* Do 128 eight_reals_fft_2 macros */ \
/*	distance between fft data elements is 256K */ \
/*	do 128 macros each processing 8 data values */ \
/* */ \
/* Do 512 four_complex_fft macros */ \
/*	distance between fft data elements is 256K */ \
/*	do 128 macros each processing 32 data values */ \
\
/* Do FFT levels 6 */ \
/* Values 0K-255K is real data, 256K-511K is semi-real data, 512K-5119K is */ \
/* complex data. */ \
/* */ \
/* On input the 64-byte cache lines hold these data values: */ \
/*	0K	256K	512K	768K */ \
/*	2K	... */ \
/*	... */ \
/*	254K	... */ \
/*	1024K	... */ \
/*	... */ \
/* On output the 64-byte cache lines hold these data values: */ \
/*	0K	128K	256K	384K */ \
/*	2K	... */ \
/*	... */ \
/*	126K	... */ \
/*	512K	... */ \
/*	... */ \
/* */ \
/* Do 64 eight_reals_fft_1 macros */ \
/* Do 64 two_two_complex_fft macros */ \
/*	distance between fft data elements is 128K */ \
/*	do 64 macros each processing 16 data values */ \
/* */ \
/* Do 512 two_two_complex_fft macros */ \
/*	distance between fft data elements is 128K */ \
/*	do 256 macros each processing 16 data values */ \
\
yes_fft:start_timer(22); \
rsi = (uintptr_t)g->scratch_area;	/* Get address of scratch area */ \
rcx = (uintptr_t)g->data_addr; /* Load source address for prefetching */ \
pfing(rbp = (uintptr_t)g->data_prefetch);/* Set prefetch distance */ \
pfing(rbp -= rcx); \
pfing(rcx += 337*blkdst);		/* Create prefetching pointer */ \
rdi = (uintptr_t)g->sincos5;/* Load sin/cos pointer */ \
rbx = (uintptr_t)g->sincos4;/* Load sin/cos pointer */ \
/* 8 iterations */ \
\
for(unsigned int loopA = 8; loopA; loopA--) { \
	for(unsigned int loopB = 2*4; loopB; loopB--) { \
		for(unsigned int loopC = clm; loopC; loopC--) { \
			xtouchlines(3);			/* Preload TLBs */ \
			\
			xfive_reals_fft_preload;		/* Preload consts for five real macros */ \
			for(unsigned int loopD = 2; loopD; loopD--) { \
				xloop1sp(b1b, 8*clmblkdst8, x5cl_five_reals_fft(rsi, 64, 16*clmblkdst8)); \
			} \
			\
			rsi += -2*8*clmblkdst8;	/* Restore source pointer */ \
			for(unsigned int loopD = 2; loopD; loopD--) { \
				for(unsigned int loopE = 2; loopE; loopE--) { \
					x1cl_half_eight_reals_fft_2(rsi, 64); \
				} \
				rsi += -2*64+8*clmblkdst8;	/* Next source pointer */ \
			} \
			\
			rsi += -2*8*clmblkdst8+16*clmblkdst8; /* Next source pointer */ \
			xloop2sp(b3b, 8*clmblkdst8, x4cl_four_complex_fft(rsi, 64, 16*clmblkdst8, 32*clmblkdst8, rdi)); \
			\
			rsi += -2*8*clmblkdst8-16*clmblkdst8;/* Restore src pointer */ \
			for(unsigned int loopD = 2; loopD; loopD--) { \
				x2cl_eight_reals_fft_1(rsi, 64, 8*clmblkdst8); \
			} \
			\
			rsi += -2*64+16*clmblkdst8;/* Next source pointer */ \
			rbx += XMM_SCD;		/* Next sine/cosine pointer */ \
			for(unsigned int loopD = 4; loopD; loopD--) { \
				for(unsigned int loopE = 2; loopE; loopE--) { \
					x2cl_two_complex_fft(rsi, 64, 8*clmblkdst8, rbx); \
				} \
				rsi += -2*64+16*clmblkdst8;/* Next source pointer */ \
				rbx += +2*XMM_SCD;		/* Next sine/cosine pointer */ \
			} \
			rbx += -9*XMM_SCD;		/* Restore sine/cosine pointer */ \
			\
				EXEC248##clm(rsi += -5*16*clmblkdst8+2*64);/* Next source pointer */ \
			EXEC248##clm(pfing(rcx += -3*blkdst+128));	/* Next prefetch pointer */ \
		} \
		EXEC248##clm(rsi += -clm*2*64+clmblkdst);/* Next source pointer */ \
		EXEC248##clm(pfing(rcx += -clm*128+3*blkdst));	/* Next prefetch pointer */ \
		EXEC1##clm(rsi += -5*16*clmblkdst8+clmblkdst);/* Next source pointer */ \
		\
	} \
	rsi += -8*clmblkdst+clmblkdst8;/* Next source pointer */ \
} \
rsi += -8*clmblkdst8;	/* Restore source pointer */ \
end_timer(22); \
\
/* Do FFT levels 7,8 */ \
/* Values 0K-127K is real data, 128K-255K is semi-real data, 256K-5119K is */ \
/* complex data. */ \
/* */ \
/* On input the 64-byte cache lines hold these data values: */ \
/*	0K	128K	256K	384K */ \
/*	2K	... */ \
/*	... */ \
/*	126K	... */ \
/*	512K	... */ \
/*	... */ \
/* On output the 64-byte cache lines hold these data values: */ \
/*	0K	32K	64K	96K */ \
/*	2K	... */ \
/*	... */ \
/*	30K	... */ \
/*	128K	... */ \
/*	... */ \
\
/* Do 32 eight_reals_fft_2 macros */ \
/* Do 32 two_two_complex_fft macros */ \
/*	distance between fft data elements is 32K */ \
/*	do 16 macros each processing 32 data values */ \
\
start_timer(24); \
rdi = (uintptr_t)g->sincos3;/* Load sin/cos pointer */ \
/* 2 iterations of 8 iters of 2*clm */ \
for(unsigned int loopA = 2; loopA; loopA--) { \
	for(unsigned int loopB = 8; loopB; loopB--) { \
		for(unsigned int loopC = 2*clm; loopC; loopC--) { \
			x4cl_eight_reals_fft_2(rsi, 64, 2*clmblkdst8, 4*clmblkdst8, rdi); \
		} \
		rsi += -2*clm*64+clmblkdst;/* Next source pointer */ \
	} \
	rsi += -8*clmblkdst+clmblkdst8;/* Next source pointer */ \
} \
\
/* Do 576 four_complex_fft macros */ \
/*	distance between fft data elements is 32K */ \
/*	do 144 macros each processing 32 data values */ \
\
rsi += -2*clmblkdst8+8*clmblkdst8;/* Load source pointer */ \
rdi += XMM_SCD;		/* Next sine/cosine pointer */ \
/* 9 iters of 2*4 iters of 2*2*clm */ \
for(unsigned int loopA = 9; loopA; loopA--) { \
	for(unsigned int loopB = 2; loopB; loopB--) { \
		for(unsigned int loopC = 4; loopC; loopC--) { \
			xtouch(xptr(rcx+(clm-1)*128));	/* Preload TLB */ \
			xloop2_clm##clm(b5d, clmblkdst, x4cl_four_complex_fft(rsi, 64, 2*clmblkdst8, 4*clmblkdst8, rdi)); \
		} \
		rsi += -8*clmblkdst+clmblkdst8;/* Next source pointer */ \
	} \
	rsi += -2*clmblkdst8+8*clmblkdst8;/* Next source pointer */ \
	rdi += 2*XMM_SCD;		/* Next sine/cosine pointer */ \
} \
rsi += -10*8*clmblkdst8;	/* Restore source pointer */ \
end_timer(24); \
\
/* Do FFT levels 9,10 */ \
/* Values 0K-31K is real data, 32K-63K is semi-real data, 64K-5119K is */ \
/* complex data. */ \
/* */ \
/* On input the 64-byte cache lines hold these data values: */ \
/*	0K	32K	64K	96K */ \
/*	2K	... */ \
/*	... */ \
/*	30K	... */ \
/*	128K	... */ \
/*	... */ \
/* On output the 64-byte cache lines hold these data values: */ \
/*	0K	8K	16K	24K */ \
/*	2K	... */ \
/*	... */ \
/*	6K	... */ \
/*	32K	... */ \
/*	... */ \
\
/* Do 8 eight_reals_fft_2 macros */ \
/* Do 8 two_two_complex_fft macros */ \
/*	distance between fft data elements is 8K */ \
/*	do 4 macros each processing 32 data values */ \
\
start_timer(25); \
rdi = (uintptr_t)g->sincos2;/* Load sin/cos pointer */ \
/* 4 iterations of 2*clm */ \
for(unsigned int loopA = 4; loopA; loopA--) { \
	for(unsigned int loopB = 2*clm; loopB; loopB--) { \
		x4cl_eight_reals_fft_2(rsi, 64, 4*clmblkdst, clmblkdst8, rdi); \
	} \
	rsi += -2*clm*64+clmblkdst;/* Next source pointer */ \
} \
\
/* Do 624 four_complex_fft macros */ \
/*	distance between fft data elements is 8K */ \
/*	do 156 macros each processing 32 data values */ \
\
rsi += -4*clmblkdst+2*clmblkdst8;/* Load source pointer */ \
rdi += XMM_SCD;		/* Next sine/cosine pointer */ \
/* 39 iterations of 2 iters of 2*2*clm */ \
for(unsigned int loopA = 39; loopA; loopA--) { \
	xtouch(xptr(rcx+(clm-1)*128));	/* Preload TLB */ \
	xloop4_clm##clm(b7b, clmblkdst, x4cl_four_complex_fft(rsi, 64, 4*clmblkdst, clmblkdst8, rdi)); \
	rsi += -4*clmblkdst+2*clmblkdst8;/* Next source pointer */ \
	rdi += 2*XMM_SCD;		/* Next sine/cosine pointer */ \
} \
rsi += -40*2*clmblkdst8;	/* Restore source pointer */ \
pfing(rcx += -640*blkdst);	/* Restore prefetch pointer */ \
end_timer(25); \
\
/* Do FFT levels 11,12 */ \
/* Values 0K-7K is real data, 8K-15K is semi-real data, 16K-5119K is */ \
/* complex data. */ \
/* */ \
/* On input the 64-byte cache lines hold these data values: */ \
/*	0K	8K	16K	24K */ \
/*	2K	... */ \
/*	... */ \
/*	6K	... */ \
/*	32K	... */ \
/*	... */ \
/* On output the 64-byte cache lines hold these data values: */ \
/*	0K	2K	4K	6K */ \
/*	8K	... */ \
/*	... */ \
\
/* Do 1 eight_reals_fft_2 macros */ \
/* Do 1 nop_two_two_complex_fft_2 macros */ \
/* Do 2 four_complex_fft macros */ \
/*	distance between fft data elements is 2K */ \
/*	do 1 macros each processing 32 data values */ \
\
start_timer(26); \
rdi = (uintptr_t)g->sincos1;/* Load sin/cos pointer */ \
for(unsigned int loopA = 2*clm; loopA; loopA--) { \
	g4cl_eight_reals_fft_2(rsi, 64, clmblkdst, 2*clmblkdst, rcx, 64, blkdst, 2*blkdst); \
} \
\
/* Do 636 four_complex_fft macros */ \
/*	distance between fft data elements is 2K */ \
/*	do 159 macros each processing 32 data values */ \
\
rsi += -2*clm*64+4*clmblkdst;/* Next source pointer */ \
rcx += -2*clm*64+4*blkdst;	/* Next dest pointer */ \
rdi += XMM_SCD;		/* Next sine/cosine pointer */ \
/* eax = 80*256+1; */ /* 159 iterations of 2*clm */ \
for(unsigned int loopA = 159; loopA; ) { \
	do{ \
		for(unsigned int loopC = clm; loopC; loopC--) { \
			for(unsigned int loopD = 2; loopD; loopD--) { \
				g4cl_four_complex_fft(rsi, 64, clmblkdst, 2*clmblkdst, rcx, 64, blkdst, 2*blkdst); \
				L2prefetch128(u8ptr(rcx-64+rbp)); \
				L2prefetch128(u8ptr(rcx-64+blkdst+rbp)); \
				pfing(rbp += -64+2*blkdst); \
			} \
			pfing(rbp += 128-4*blkdst); \
		} \
		rsi += -2*clm*64+4*clmblkdst;/* Next source pointer */ \
		rcx += -2*clm*64+4*blkdst;	/* Next dest pointer */ \
		rdi += 2*XMM_SCD;		/* Next sine/cosine pointer */ \
	}while(--loopA & 1); \
	rsi += -8*clmblkdst+clmblkdst8; \
} \
end_timer(26); \
\
/* Work on next 4 sets of 2560 values */ \
\
endlp:	xpass1_get_next_block(pass2, c0b, b0b); /* Get next block and jump */ \
/* done: */	xfft_footer_hg2; \
\
\
\
/* Perform a 6144K-element FFT.  This is done in two passes.  Pass 1 */ \
/* does 12 levels, pass 2 uses common code to perform the last 11 levels. */ \
\

#define hg_pass1sclevels12pfa6(procname, xpass2name, pass2_levels, clm) \
 /* (fftname, pass2_levels, clmarg) */ \
\
/* set_FFT_constants(fftname, 0, 0, pass2_levels, clmarg); */ \
xfft_header_hg2(procname, xpass2name, pass2_levels, clm); \
\
/* Do the last 11 FFT levels */ \
/*	do 1 pass2_11_levels_real and 191 pass2_11_levels_complex */ \
\
/* pass2: */	start_timer(0); \
xpass2name(g); \
end_timer(0); \
if(g->ffttype == 1)		/* We're done if FFTing only */ \
	return; \
\
/* Do the inverse FFT */ \
\
pass1_inverse_fft_setup(c0b, procname); \
\
/* Do inverse FFT levels 11,12 */ \
/* On input the 64-byte cache lines hold these data values: */ \
/*	0K	2K	4K	6K */ \
/*	8K	10K	12K	14K */ \
/*	... */ \
/* On output the 64-byte cache lines hold these data values: */ \
/*	0K	4K	8K	12K */ \
/*	2K	... */ \
/*	16K	... */ \
/*	... */ \
\
/* Do 1 eight_reals_unfft_2 macros */ \
/* Do 1 nop_two_two_complex_unfft_2 macros */ \
/*	distance between fft data elements is 2K */ \
/*	do 1 macros each processing 16 data values */ \
\
c0b:	start_timer(16); \
rsi = (uintptr_t)g->data_addr;/* Get source address */ \
rcx = (uintptr_t)g->scratch_area;	/* Get address of scratch area */ \
pfing(rbp = (uintptr_t)g->data_prefetch);/* Set prefetch distance */ \
pfing(rbp -= rsi); \
for(unsigned int loopA = 2*clm; loopA; loopA--) { \
	g2cl_eight_reals_unfft_2(rsi, 64, blkdst, rcx, 64, clmblkdst); \
} \
\
/* Do 766 four_complex_unfft macros */ \
/*	distance between fft data elements is 2K */ \
/*	do 383 macros each processing 16 data values */ \
\
rsi += -2*clm*64+2*blkdst;	/* Next source pointer */ \
rcx += -2*clm*64+2*clmblkdst;/* Next dest pointer */ \
rdi = (uintptr_t)g->sincos1;/* Load sin/cos pointer */ \
/* eax = 96*256+256/4; */ /* 383 iterations of 2*clm */ \
for(unsigned int loopA = 383; loopA; ) { \
	do{ \
		for(unsigned int loopC = 2*clm; loopC; loopC--) { \
			g2cl_four_complex_unfft(rsi, 64, blkdst, rcx, 64, clmblkdst); \
		} \
		rsi += -2*clm*64+2*blkdst;	/* Next source pointer */ \
		rcx += -2*clm*64+2*clmblkdst;/* Next dest pointer */ \
		rdi += XMM_SCD;		/* Next sine/cosine pointer */ \
	}while(--loopA & 3); \
	rcx += -8*clmblkdst+clmblkdst8; \
} \
rcx = rsi-384*2*blkdst;	/* Create prefetch pointer */ \
rsi = (uintptr_t)g->scratch_area;	/* Restore source pointer */ \
end_timer(16); \
\
/* Do inverse FFT levels 9,10 */ \
/* On input the 64-byte cache lines hold these data values: */ \
/*	0K	4K	8K	12K */ \
/*	2K	... */ \
/*	16K	... */ \
/*	... */ \
/* On output the 64-byte cache lines hold these data values: */ \
/*	0K	16K	32K	48K */ \
/*	2K	... */ \
/*	... */ \
/*	14K */ \
/*	64K	... */ \
/*	... */ \
\
/* Do 4 eight_reals_unfft_2 macros */ \
/* Do 4 nop_two_two_complex_unfft_2 macros */ \
/*	distance between fft data elements is 8K */ \
/*	do 2 macros each processing 32 data values */ \
\
start_timer(17); \
/* 2 iterations of 2*clm */ \
for(unsigned int loopA = 2; loopA; loopA--) { \
	for(unsigned int loopB = 2*clm; loopB; loopB--) { \
		x4cl_eight_reals_unfft_2(rsi, 64, 2*clmblkdst, 4*clmblkdst); \
	} \
	rsi += -2*clm*64+clmblkdst;/* Next source pointer */ \
} \
\
/* Do 760 four_complex_unfft macros */ \
/*	distance between fft data elements is 8K */ \
/*	do 190 macros each processing 32 data values */ \
\
rsi += -2*clmblkdst+clmblkdst8;/* Next source pointer */ \
rdi = (uintptr_t)g->sincos2;/* Load sin/cos pointer */ \
/* 95 iterations of 2*2*clm */ \
for(unsigned int loopA = 95; loopA; loopA--) { \
	xtouch(xptr(rcx+(clm-1)*128));	/* Preload TLB */ \
	xloop2_clm##clm(c7b, clmblkdst, x4cl_four_complex_unfft(rsi, 64, 2*clmblkdst, 4*clmblkdst, rdi)); \
	rsi += -2*clmblkdst+clmblkdst8;/* Next source pointer */ \
	rdi += XMM_SCD;		/* Next sine/cosine pointer */ \
} \
rsi += -96*clmblkdst8;	/* Restore source pointer */ \
end_timer(17); \
\
/* Do inverse FFT levels 7,8 */ \
/* On input the 64-byte cache lines hold these data values: */ \
/*	0K	16K	32K	48K */ \
/*	2K	... */ \
/*	... */ \
/*	14K */ \
/*	64K	... */ \
/*	... */ \
/* On output the 64-byte cache lines hold these data values: */ \
/*	0K	64K	128K	192K */ \
/*	2K	... */ \
/*	... */ \
/*	62K */ \
/*	256K	... */ \
/*	... */ \
\
/* Do 16 eight_reals_unfft_2 macros */ \
/* Do 16 nop_two_two_complex_unfft_2 macros */ \
/*	distance between fft data elements is 32K */ \
/*	do 8 macros each processing 32 data values */ \
\
start_timer(18); \
/* 8 iterations of 2*clm */ \
for(unsigned int loopA = 8; loopA; loopA--) { \
	for(unsigned int loopB = 2*clm; loopB; loopB--) { \
		x4cl_eight_reals_unfft_2(rsi, 64, clmblkdst8, 2*clmblkdst8); \
	} \
	rsi += -2*clm*64+clmblkdst;/* Next source pointer */ \
} \
\
/* Do 736 four_complex_unfft macros */ \
/*	distance between fft data elements is 32K */ \
/*	do 184 macros each processing 32 data values */ \
\
rsi += -8*clmblkdst+4*clmblkdst8;/* Next source pointer */ \
rdi = (uintptr_t)g->sincos3;/* Load sin/cos pointer */ \
/* 23 iterations of 4 iters of 2*2*clm */ \
for(unsigned int loopA = 23; loopA; loopA--) { \
	for(unsigned int loopB = 4; loopB; loopB--) { \
		xtouch(xptr(rcx+(clm-1)*128));	/* Preload TLB */ \
		xloop2_clm##clm(c5d, clmblkdst, x4cl_four_complex_unfft(rsi, 64, clmblkdst8, 2*clmblkdst8, rdi)); \
	} \
	rsi += -4*2*clmblkdst+4*clmblkdst8;/* Next source pointer */ \
	rdi += XMM_SCD;		/* Next sine/cosine pointer */ \
} \
rsi += -24*4*clmblkdst8;	/* Restore source pointer */ \
end_timer(18); \
\
/* Do inverse FFT levels 6 */ \
/* On input the 64-byte cache lines hold these data values: */ \
/*	0K	64K	128K	192K */ \
/*	2K	... */ \
/*	... */ \
/*	62K */ \
/*	256K	... */ \
/*	... */ \
/* On output the 64-byte cache lines hold these data values: */ \
/*	0K	128K	256K	384K */ \
/*	2K	... */ \
/*	... */ \
/*	126K */ \
/*	512K	... */ \
/*	... */ \
/* */ \
/* Do 64 eight_reals_unfft_1 macros */ \
/*	distance between fft data elements is 128K */ \
/*	do 32 macros each processing 16 data values */ \
/* */ \
/* Do 704 four_complex_unfft macros */ \
/*	distance between fft data elements is 128K */ \
/*	do 352 macros each processing 16 data values */ \
\
/* Do inverse FFT levels 4,5 */ \
/* On input the 64-byte cache lines hold these data values: */ \
/*	0K	128K	256K	384K */ \
/*	2K	... */ \
/*	... */ \
/*	126K */ \
/*	512K	... */ \
/*	... */ \
/* On output the 64-byte cache lines hold these data values: */ \
/*	0K	512K	1024K	1536K */ \
/*	2K	... */ \
/*	... */ \
/*	510K */ \
/*	2048K	3072K	4096K	5120K */ \
/*	2050K */ \
/*	... */ \
/*	3070K */ \
/* */ \
/* Do 128 eight_reals_unfft_2 macros */ \
/* Do 128 nop_two_two_complex_unfft_2 macros */ \
/*	distance between fft data elements is 256K */ \
/*	do 64 macros each processing 32 data values */ \
/* */ \
/* Do 512 four_complex_unfft macros */ \
/*	distance between fft data elements is 256K */ \
/*	do 128 macros each processing 32 data values */ \
\
/* Do inverse FFT levels 1,2,3 */ \
/* On input the 64-byte cache lines hold these data values: */ \
/*	0K	512K	1024K	1536K */ \
/*	2K	... */ \
/*	... */ \
/*	510K */ \
/*	2048K	3072K	4096K	5120K */ \
/*	2050K */ \
/*	... */ \
/*	3070K */ \
/* On output the 64-byte cache lines hold these data values: */ \
/*	0K	1536K	3072K	4608K */ \
/*	... */ \
/*	1534K */ \
/* */ \
/* Do 1024 six_reals_last_unfft macros */ \
/*	distance between fft data elements is 1024K */ \
/*	do 256 macros each processing 24 data values */ \
\
start_timer(19); \
rdi = (uintptr_t)g->sincos4;/* Load sin/cos pointer */ \
rbx = (uintptr_t)g->sincos5;/* Load sin/cos pointer */ \
/* 4 iterations  */ \
xsix_reals_unfft_preload; \
\
for(unsigned int loopA = 4; loopA; loopA--) { \
	for(unsigned int loopB = 2*4; loopB; loopB--) { \
		xtouchlines(7);			/* Preload TLBs */ \
		\
		for(unsigned int loopC = clm; loopC; loopC--) { \
			xloop1sp(c4b, 8*clmblkdst8, x2cl_eight_reals_unfft_1(rsi, 64, 4*clmblkdst8)); \
			\
			/* 11 iterations */ \
			for(unsigned int loopD = 11; loopD; loopD--) { \
				for(unsigned int loopE = 2; loopE; loopE--) { \
					x2cl_two_complex_unfft(rsi, 64, 4*clmblkdst8); \
				} \
				rsi += -2*64+8*clmblkdst8;	/* Next source pointer */ \
				rdi += XMM_SCD;		/* Next sine/cosine pointer */ \
			} \
			rdi += -11*XMM_SCD;	/* Restore sine/cosine pointer */ \
			\
			rsi += -12*8*clmblkdst8;	/* Restore source pointer */ \
			for(unsigned int loopD = 2; loopD; loopD--) { \
				for(unsigned int loopE = 2; loopE; loopE--) { \
					x4cl_eight_reals_unfft_2(rsi, 64, 8*clmblkdst8, 16*clmblkdst8); \
				} \
				rsi += -2*64+4*clmblkdst8; \
			} \
			\
			rsi += -8*clmblkdst8+32*clmblkdst8;/* Next source pointer */ \
			for(unsigned int loopD = 2; loopD; loopD--) { \
				xloop2sp(c3b, 4*clmblkdst8, x4cl_four_complex_unfft(rsi, 64, 8*clmblkdst8, 16*clmblkdst8, rbx)); \
				rsi += -2*4*clmblkdst8+32*clmblkdst8; \
				rbx += XMM_SCD;		/* Next sine/cosine pointer */ \
			} \
			rbx += -2*XMM_SCD;		/* Restore sine/cosine pointer */ \
			\
			rsi += -3*32*clmblkdst8;	/* Restore source pointer */ \
			for(unsigned int loopD = 4; loopD; loopD--) { \
				xloop2sp(c1b, 4*clmblkdst8, x3cl_six_reals_last_unfft(rsi, 64, 32*clmblkdst8)); \
			} \
			\
				EXEC248##clm(rsi += -8*4*clmblkdst8+2*64);/* Next source pointer */ \
			EXEC248##clm(pfing(rcx += -7*blkdst+128));	/* Next prefetch pointer */ \
		} \
		EXEC248##clm(rsi += -clm*2*64+clmblkdst);/* Next source pointer */ \
		EXEC248##clm(pfing(rcx += -clm*128+7*blkdst));	/* Next prefetch pointer */ \
		EXEC1##clm(rsi += -8*4*clmblkdst8+clmblkdst);/* Next source pointer */ \
		\
	} \
	rsi += -8*clmblkdst+clmblkdst8;/* Next source pointer */ \
} \
end_timer(19); \
\
/* Normalize these values */ \
\
pass1_normalize(1, yes_fft); \
\
/* Copy data from the scratch area back to the FFT data area */ \
\
copy_scratch_data_to_fft(768, clm); \
goto	endlp;			/* Skip FFT code if flag not set */ \
\
/* Fire up auxiliary threads to do pass 1 forward FFT */ \
\
/* pass1: */	pass1_forward_fft_setup(b0b); \
\
/* Copy data from the FFT data area back to the scratch area */ \
\
b0b:	copy_fft_data_to_scratch(768, clm); \
\
/* Do FFT levels 1,2,3 */ \
/* */ \
/* On input the 64-byte cache lines hold these data values: */ \
/*	0K	1536K	3072K	4608K */ \
/*	... */ \
/*	1534K */ \
/* On output the 64-byte cache lines hold these data values: */ \
/*	0K	512K	1024K	1536K */ \
/*	2K	... */ \
/*	... */ \
/*	510K */ \
/*	2048K	3072K	4096K	5120K */ \
/*	2050K */ \
/*	... */ \
/*	3070K */ \
/* */ \
/* Do 1024 six_reals_first_fft macros (FFT levels 1,2,3) */ \
/*	distance between fft data elements is 1024K */ \
/*	do 256 iterations each processing 24 data values */ \
\
/* Do FFT levels 4,5 */ \
/* Values 0K-1023K is real data, 1024K-2047K is semi-real data, 2048K-6143K is */ \
/* complex data. */ \
/* */ \
/* On input the 64-byte cache lines hold these data values: */ \
/*	0K	512K	1024K	1536K */ \
/*	2K	... */ \
/*	... */ \
/*	510K */ \
/*	2048K	3072K	4096K	5120K */ \
/*	2050K */ \
/*	... */ \
/*	3070K */ \
/* On output the 64-byte cache lines hold these data values: */ \
/*	0K	256K	512K	768K */ \
/*	2K	... */ \
/*	... */ \
/*	254K	... */ \
/*	1024K	... */ \
/*	... */ \
/* */ \
/* Do 128 eight_reals_fft_2 macros */ \
/* Do 128 nop_two_two_complex_fft_2 macros */ \
/*	distance between fft data elements is 256K */ \
/*	do 128 macros each processing 16 data values */ \
/* */ \
/* Do 512 four_complex_fft macros */ \
/*	distance between fft data elements is 256K */ \
/*	do 128 macros each processing 32 data values */ \
\
/* Do FFT levels 6 */ \
/* Values 0K-255K is real data, 256K-511K is semi-real data, 512K-6143K is */ \
/* complex data. */ \
/* */ \
/* On input the 64-byte cache lines hold these data values: */ \
/*	0K	256K	512K	768K */ \
/*	2K	... */ \
/*	... */ \
/*	254K	... */ \
/*	1024K	... */ \
/*	... */ \
/* On output the 64-byte cache lines hold these data values: */ \
/*	0K	128K	256K	384K */ \
/*	2K	... */ \
/*	... */ \
/*	126K	... */ \
/*	512K	... */ \
/*	... */ \
/* */ \
/* Do 64 eight_reals_fft_1 macros */ \
/* Do 64 two_two_complex_fft macros */ \
/*	distance between fft data elements is 128K */ \
/*	do 64 macros each processing 16 data values */ \
/* */ \
/* Do 640 two_two_complex_fft macros */ \
/*	distance between fft data elements is 128K */ \
/*	do 320 macros each processing 16 data values */ \
\
yes_fft:start_timer(22); \
rsi = (uintptr_t)g->scratch_area;	/* Get address of scratch area */ \
rcx = (uintptr_t)g->data_addr; /* Load source address for prefetching */ \
pfing(rbp = (uintptr_t)g->data_prefetch);/* Set prefetch distance */ \
pfing(rbp -= rcx); \
pfing(rcx += 394*blkdst);		/* Create prefetch pointer */ \
rdi = (uintptr_t)g->sincos5;/* Load sin/cos pointer */ \
rbx = (uintptr_t)g->sincos4;/* Load sin/cos pointer */ \
/* 8 iterations */ \
xsix_reals_fft_preload; \
\
for(unsigned int loopA = 8; loopA; loopA--) { \
	for(unsigned int loopB = 2*4; loopB; loopB--) { \
		xtouchlines(3);			/* Preload TLBs */ \
		\
		for(unsigned int loopC = clm; loopC; loopC--) { \
			xloop4sp(b1b, 8*clmblkdst8, x3cl_six_reals_fft(rsi, 64, 32*clmblkdst8)); \
			\
			rsi += -4*8*clmblkdst8;	/* Restore source pointer */ \
			for(unsigned int loopD = 2; loopD; loopD--) { \
				for(unsigned int loopE = 2; loopE; loopE--) { \
					x2cl_eight_reals_fft_2(rsi, 64, 16*clmblkdst8); \
				} \
				rsi += -2*64+8*clmblkdst8;	/* Next source pointer */ \
			} \
			\
			rsi += -2*8*clmblkdst8+32*clmblkdst8;/* Next source pointer */ \
			for(unsigned int loopD = 2; loopD; loopD--) { \
				xloop1sp(b3b, 8*clmblkdst8, x4cl_four_complex_fft(rsi, 64, 16*clmblkdst8, 32*clmblkdst8, rdi)); \
			} \
			\
			rsi += -2*8*clmblkdst8-32*clmblkdst8;/* Restore src pointer */ \
			for(unsigned int loopD = 2; loopD; loopD--) { \
				x2cl_eight_reals_fft_1(rsi, 64, 8*clmblkdst8); \
			} \
			\
			rsi += -2*64+16*clmblkdst8;/* Load source pointer */ \
			rbx += XMM_SCD;		/* Next sine/cosine pointer */ \
			/* 5 iterations */ \
			for(unsigned int loopD = 5; loopD; loopD--) { \
				for(unsigned int loopE = 2; loopE; loopE--) { \
					x2cl_two_complex_fft(rsi, 64, 8*clmblkdst8, rbx); \
				} \
				rsi += -2*64+16*clmblkdst8;/* Next source pointer */ \
				rbx += 2*XMM_SCD;		/* Next sine/cosine pointer */ \
			} \
			rbx += -11*XMM_SCD;	/* Restore sine/cosine pointer */ \
			\
				EXEC248##clm(rsi += -6*16*clmblkdst8+2*64);/* Next source pointer */ \
			EXEC248##clm(pfing(rcx += -3*blkdst+128));	/* Next prefetch pointer */ \
		} \
		EXEC248##clm(rsi += -clm*2*64+clmblkdst);/* Next source pointer */ \
		EXEC248##clm(pfing(rcx += -clm*128+3*blkdst));	/* Next prefetch pointer */ \
		EXEC1##clm(rsi += -6*16*clmblkdst8+clmblkdst);/* Next source pointer */ \
		\
	} \
	rsi += -8*clmblkdst+clmblkdst8;/* Next source pointer */ \
} \
rsi += -8*clmblkdst8;	/* Restore source pointer */ \
end_timer(22); \
\
/* Do FFT levels 7,8 */ \
/* Values 0K-127K is real data, 128K-255K is semi-real data, 256K-6143K is */ \
/* complex data. */ \
/* */ \
/* On input the 64-byte cache lines hold these data values: */ \
/*	0K	128K	256K	384K */ \
/*	2K	... */ \
/*	... */ \
/*	126K	... */ \
/*	512K	... */ \
/*	... */ \
/* On output the 64-byte cache lines hold these data values: */ \
/*	0K	32K	64K	96K */ \
/*	2K	... */ \
/*	... */ \
/*	30K	... */ \
/*	128K	... */ \
/*	... */ \
\
/* Do 32 eight_reals_fft_2 macros */ \
/* Do 32 two_two_complex_fft macros */ \
/*	distance between fft data elements is 32K */ \
/*	do 16 macros each processing 32 data values */ \
\
start_timer(24); \
rdi = (uintptr_t)g->sincos3;/* Load sin/cos pointer */ \
/* 2 iterations of 8 iters of 2*clm */ \
for(unsigned int loopA = 2; loopA; loopA--) { \
	for(unsigned int loopB = 8; loopB; loopB--) { \
		for(unsigned int loopC = 2*clm; loopC; loopC--) { \
			x4cl_eight_reals_fft_2(rsi, 64, 2*clmblkdst8, 4*clmblkdst8, rdi); \
		} \
		rsi += -2*clm*64+clmblkdst;/* Next source pointer */ \
	} \
	rsi += -8*clmblkdst+clmblkdst8;/* Next source pointer */ \
} \
\
/* Do 704 four_complex_fft macros */ \
/*	distance between fft data elements is 32K */ \
/*	do 176 macros each processing 32 data values */ \
\
rsi += -2*clmblkdst8+8*clmblkdst8;/* Load source pointer */ \
rdi += XMM_SCD;		/* Next sine/cosine pointer */ \
/* 11 iters of 2*4 iters of 2*2*clm */ \
for(unsigned int loopA = 11; loopA; loopA--) { \
	for(unsigned int loopB = 2; loopB; loopB--) { \
		for(unsigned int loopC = 4; loopC; loopC--) { \
			xtouch(xptr(rcx+(clm-1)*128));	/* Preload TLB */ \
			xloop2_clm##clm(b5d, clmblkdst, x4cl_four_complex_fft(rsi, 64, 2*clmblkdst8, 4*clmblkdst8, rdi)); \
		} \
		rsi += -8*clmblkdst+clmblkdst8;/* Next source pointer */ \
	} \
	rsi += -2*clmblkdst8+8*clmblkdst8;/* Next source pointer */ \
	rdi += 2*XMM_SCD;		/* Next sine/cosine pointer */ \
} \
rsi += -12*8*clmblkdst8;	/* Restore source pointer */ \
end_timer(24); \
\
/* Do FFT levels 9,10 */ \
/* Values 0K-31K is real data, 32K-63K is semi-real data, 64K-6143K is */ \
/* complex data. */ \
/* */ \
/* On input the 64-byte cache lines hold these data values: */ \
/*	0K	32K	64K	96K */ \
/*	2K	... */ \
/*	... */ \
/*	30K	... */ \
/*	128K	... */ \
/*	... */ \
/* On output the 64-byte cache lines hold these data values: */ \
/*	0K	8K	16K	24K */ \
/*	2K	... */ \
/*	... */ \
/*	6K	... */ \
/*	32K	... */ \
/*	... */ \
\
/* Do 8 eight_reals_fft_2 macros */ \
/* Do 8 two_two_complex_fft macros */ \
/*	distance between fft data elements is 8K */ \
/*	do 4 macros each processing 32 data values */ \
\
start_timer(25); \
rdi = (uintptr_t)g->sincos2;/* Load sin/cos pointer */ \
/* 4 iterations of 2*clm */ \
for(unsigned int loopA = 4; loopA; loopA--) { \
	for(unsigned int loopB = 2*clm; loopB; loopB--) { \
		x4cl_eight_reals_fft_2(rsi, 64, 4*clmblkdst, clmblkdst8, rdi); \
	} \
	rsi += -2*clm*64+clmblkdst;/* Next source pointer */ \
} \
\
/* Do 752 four_complex_fft macros */ \
/*	distance between fft data elements is 8K */ \
/*	do 188 macros each processing 32 data values */ \
\
rsi += -4*clmblkdst+2*clmblkdst8;/* Load source pointer */ \
rdi += XMM_SCD;		/* Next sine/cosine pointer */ \
/* 47 iterations of 4 iters of 2*clm */ \
for(unsigned int loopA = 47; loopA; loopA--) { \
	for(unsigned int loopB = 2; loopB; loopB--) { \
		xtouch(xptr(rcx+(clm-1)*128));	/* Preload TLB */ \
		xloop2_clm##clm(b7b, clmblkdst, x4cl_four_complex_fft(rsi, 64, 4*clmblkdst, clmblkdst8, rdi)); \
	} \
	rsi += -4*clmblkdst+2*clmblkdst8;/* Next source pointer */ \
	rdi += 2*XMM_SCD;		/* Next sine/cosine pointer */ \
} \
rsi += -48*2*clmblkdst8;	/* Restore source pointer */ \
pfing(rcx += -768*blkdst);	/* Restore prefetch pointer */ \
end_timer(25); \
\
/* Do FFT levels 11,12 */ \
/* Values 0K-7K is real data, 8K-15K is semi-real data, 16K-6143K is */ \
/* complex data. */ \
/* */ \
/* On input the 64-byte cache lines hold these data values: */ \
/*	0K	8K	16K	24K */ \
/*	2K	... */ \
/*	... */ \
/*	6K	... */ \
/*	32K	... */ \
/*	... */ \
/* On output the 64-byte cache lines hold these data values: */ \
/*	0K	2K	4K	6K */ \
/*	8K	... */ \
/*	... */ \
\
/* Do 1 eight_reals_fft_2 macros */ \
/* Do 1 nop_two_two_complex_fft_2 macros */ \
/* Do 2 four_complex_fft macros */ \
/*	distance between fft data elements is 2K */ \
/*	do 1 macros each processing 32 data values */ \
\
start_timer(26); \
rdi = (uintptr_t)g->sincos1;/* Load sin/cos pointer */ \
for(unsigned int loopA = 2*clm; loopA; loopA--) { \
	g4cl_eight_reals_fft_2(rsi, 64, clmblkdst, 2*clmblkdst, rcx, 64, blkdst, 2*blkdst); \
} \
\
/* Do 764 four_complex_fft macros */ \
/*	distance between fft data elements is 2K */ \
/*	do 191 macros each processing 32 data values */ \
\
rsi += -2*clm*64+4*clmblkdst;/* Next source pointer */ \
rcx += -2*clm*64+4*blkdst;	/* Next dest pointer */ \
rdi += XMM_SCD;		/* Next sine/cosine pointer */ \
/* eax = 96*256+1; */ /* 95 iterations of 2*clm     WARNING in this comment */ \
for(unsigned int loopA = 191; loopA; ) { \
	do{ \
		for(unsigned int loopC = clm; loopC; loopC--) { \
			for(unsigned int loopD = 2; loopD; loopD--) { \
				g4cl_four_complex_fft(rsi, 64, clmblkdst, 2*clmblkdst, rcx, 64, blkdst, 2*blkdst); \
				L2prefetch128(u8ptr(rcx-64+rbp)); \
				L2prefetch128(u8ptr(rcx-64+blkdst+rbp)); \
				pfing(rbp += -64+2*blkdst); \
			} \
			pfing(rbp += 128-4*blkdst); \
		} \
		rsi += -2*clm*64+4*clmblkdst;/* Next source pointer */ \
		rcx += -2*clm*64+4*blkdst;	/* Next dest pointer */ \
		rdi += 2*XMM_SCD;		/* Next sine/cosine pointer */ \
	}while(--loopA & 1); \
	rsi += -8*clmblkdst+clmblkdst8; \
} \
end_timer(26); \
\
/* Work on next 4 sets of 768 values */ \
\
endlp:	xpass1_get_next_block(pass2, c0b, b0b); /* Get next block and jump */ \
/* done: */	xfft_footer_hg2; \
\
\
\
/* Perform a 6144K-element negacyclic FFT.  This is done in two passes.  Pass */ \
/* 1 does 12 levels, pass 2 uses common code to perform the last 11 levels. */ \
\

#define hg_pass1sclevels12complex3(procname, xpass2name, pass2_levels, clm) \
 /* (fftname, pass2_levels, clmarg) */ \
\
/* set_FFT_constants(fftname, 0, 1, pass2_levels, clmarg); */ \
xfft_header_hg2(procname, xpass2name, pass2_levels, clm); \
\
/* Do the last 11 FFT levels */ \
/*	do 192 pass2_11_levels_complex */ \
\
/* pass2: */	start_timer(0); \
xpass2name(g); \
end_timer(0); \
if(g->ffttype == 1)		/* We're done if FFTing only */ \
	return; \
\
/* Do the inverse FFT */ \
\
pass1_inverse_fft_setup(c0b, procname); \
\
/* Do inverse FFT levels 11,12 */ \
/* On input the 64-byte cache lines hold these data values: */ \
/*	0K	2K	4K	6K */ \
/*	8K	10K	12K	14K */ \
/*	... */ \
/* On output the 64-byte cache lines hold these data values: */ \
/*	0K	4K	8K	12K */ \
/*	2K	... */ \
/*	16K	... */ \
/*	... */ \
\
/* Do 768 four_complex_unfft macros */ \
/*	distance between fft data elements is 2K */ \
/*	do 384 macros each processing 16 data values */ \
\
c0b:	start_timer(15); \
rsi = (uintptr_t)g->data_addr;/* Get source address */ \
rcx = (uintptr_t)g->scratch_area;	/* Get address of scratch area */ \
rdi = (uintptr_t)g->sincos1;/* Load sin/cos pointer */ \
/* 96 iterations of 4 iters of 2*clm */ \
pfing(rbp = (uintptr_t)g->data_prefetch);/* Set prefetch distance */ \
pfing(rbp -= rsi); \
for(unsigned int loopA = 96; loopA; loopA--) { \
	for(unsigned int loopB = 4; loopB; loopB--) { \
		for(unsigned int loopC = 2*clm; loopC; loopC--) { \
			g2cl_four_complex_unfft(rsi, 64, blkdst, rcx, 64, clmblkdst); \
		} \
		rsi += -2*clm*64+2*blkdst;	/* Next source pointer */ \
		rcx += -2*clm*64+2*clmblkdst;/* Next dest pointer */ \
		rdi += XMM_SCD;		/* Next sine/cosine pointer */ \
	} \
	rcx += -8*clmblkdst+clmblkdst8; \
} \
rcx = rsi-384*2*blkdst;	/* Create prefetch pointer */ \
rsi = (uintptr_t)g->scratch_area;	/* Restore source pointer */ \
end_timer(15); \
\
/* Do inverse FFT levels 9,10 */ \
/* On input the 64-byte cache lines hold these data values: */ \
/*	0K	4K	8K	12K */ \
/*	2K	... */ \
/*	16K	... */ \
/*	... */ \
/* On output the 64-byte cache lines hold these data values: */ \
/*	0K	16K	32K	48K */ \
/*	2K	... */ \
/*	... */ \
/*	14K */ \
/*	64K	... */ \
/*	... */ \
\
/* Do 768 four_complex_unfft macros */ \
/*	distance between fft data elements is 8K */ \
/*	do 192 macros each processing 32 data values */ \
\
start_timer(16); \
rdi = (uintptr_t)g->sincos2;/* Load sin/cos pointer */ \
/* 96 iterations of 2*2*clm */ \
for(unsigned int loopA = 96; loopA; loopA--) { \
	xtouch(xptr(rcx+(clm-1)*128));	/* Preload TLB */ \
	xloop2_clm##clm(c7b, clmblkdst, x4cl_four_complex_unfft(rsi, 64, 2*clmblkdst, 4*clmblkdst, rdi));		/* 2*2*clm iterations */ \
	rsi += -2*clmblkdst+clmblkdst8;/* Next source pointer */ \
	rdi += XMM_SCD;		/* Next sine/cosine pointer */ \
} \
rsi += -96*clmblkdst8;	/* Restore source pointer */ \
end_timer(16); \
\
/* Do inverse FFT levels 7,8 */ \
/* On input the 64-byte cache lines hold these data values: */ \
/*	0K	16K	32K	48K */ \
/*	2K	... */ \
/*	... */ \
/*	14K */ \
/*	64K	... */ \
/*	... */ \
/* On output the 64-byte cache lines hold these data values: */ \
/*	0K	64K	128K	192K */ \
/*	2K	... */ \
/*	... */ \
/*	62K */ \
/*	256K	... */ \
/*	... */ \
\
/* Do 768 four_complex_unfft macros */ \
/*	distance between fft data elements is 32K */ \
/*	do 192 macros each processing 32 data values */ \
\
start_timer(17); \
rdi = (uintptr_t)g->sincos3;/* Load sin/cos pointer */ \
/* 24 iterations of 4 iters of 2*2*clm */ \
for(unsigned int loopA = 24; loopA; loopA--) { \
	for(unsigned int loopB = 4; loopB; loopB--) { \
		xtouch(xptr(rcx+(clm-1)*128));	/* Preload TLB */ \
		xloop2_clm##clm(c5b, clmblkdst, x4cl_four_complex_unfft(rsi, 64, clmblkdst8, 2*clmblkdst8, rdi));		/* 2*2*clm iterations */ \
	} \
	rsi += -4*2*clmblkdst+4*clmblkdst8;/* Next source pointer */ \
	rdi += XMM_SCD;		/* Next sine/cosine pointer */ \
} \
rsi += -24*4*clmblkdst8;	/* Restore source pointer */ \
end_timer(17); \
\
/* Do inverse FFT levels 5,6 */ \
/* On input the 64-byte cache lines hold these data values: */ \
/*	0K	64K	128K	192K */ \
/*	2K	... */ \
/*	... */ \
/*	62K */ \
/*	256K	... */ \
/*	... */ \
/* On output the 64-byte cache lines hold these data values: */ \
/*	0K	256K	512K	768K */ \
/*	2K	... */ \
/*	... */ \
/*	254K */ \
/*	1024K	... */ \
/*	... */ \
\
/* Do 768 four_complex_unfft macros */ \
/*	distance between fft data elements is 128K */ \
/*	do 192 macros each processing 32 data values */ \
\
start_timer(18); \
rdi = (uintptr_t)g->sincos4;/* Load sin/cos pointer */ \
/* 6 iters of 4 iters of 4 of 2*2*clm */ \
for(unsigned int loopA = 6; loopA; loopA--) { \
	for(unsigned int loopB = 4; loopB; loopB--) { \
		for(unsigned int loopC = 4; loopC; loopC--) { \
			xtouch(xptr(rcx+(clm-1)*128));	/* Preload TLB */ \
			xloop2_clm##clm(c4b, clmblkdst, x4cl_four_complex_unfft(rsi, 64, 4*clmblkdst8, 8*clmblkdst8, rdi));		/* 2 * 2*clm iterations */ \
		} \
		rsi += -4*2*clmblkdst+clmblkdst8;/* Next source pointer */ \
	} \
	rsi += -4*clmblkdst8+16*clmblkdst8;/* Next source pointer */ \
	rdi += XMM_SCD;		/* Next sine/cosine pointer */ \
} \
rsi += -6*16*clmblkdst8;	/* Restore source pointer */ \
end_timer(19); \
\
/* Do inverse FFT levels 4 */ \
/* On input the 64-byte cache lines hold these data values: */ \
/*	0K	256K	512K	768K */ \
/*	2K	... */ \
/*	... */ \
/*	254K */ \
/*	1024K	... */ \
/*	... */ \
/* On output the 64-byte cache lines hold these data values: */ \
/*	0K	512K	1024K	1536K */ \
/*	2K	... */ \
/*	... */ \
/*	510K */ \
/*	2048K	3072K	4096K	5120K */ \
/*	2050K */ \
/*	... */ \
/*	3070K */ \
\
/* Do 768 two_two_complex_unfft macros */ \
/*	distance between fft data elements is 256K */ \
/*	do 384 macros each processing 16 data values */ \
\
start_timer(19); \
rdi = (uintptr_t)g->sincos5;/* Load sin/cos pointer */ \
/* 3 iters of 16 iters of 8*2*clm */ \
for(unsigned int loopA = 3; loopA; loopA--) { \
	for(unsigned int loopB = 16; loopB; loopB--) { \
		xtouch(xptr(rcx+(clm-1)*128));	/* Preload TLB */ \
		xloop8_clm##clm(c3b, clmblkdst, x2cl_two_complex_unfft(rsi, 64, 16*clmblkdst8));		/* 8*2*clm iterations */ \
		rsi += -8*clmblkdst+clmblkdst8;/* Next source pointer */ \
	} \
	rsi += -16*clmblkdst8+32*clmblkdst8;/* Next source pointer */ \
	rdi += XMM_SCD;		/* Next sine/cosine pointer */ \
} \
rsi += -3*32*clmblkdst8;	/* Restore source pointer */ \
end_timer(19); \
\
/* Do inverse FFT levels 1,2,3 */ \
/* On input the 64-byte cache lines hold these data values: */ \
/*	0K	512K	1024K	1536K */ \
/*	2K	... */ \
/*	... */ \
/*	510K */ \
/*	2048K	3072K	4096K	5120K */ \
/*	2050K */ \
/*	... */ \
/*	3070K */ \
/* On output the 64-byte cache lines hold these data values: */ \
/*	0K	1536K	3072K	4608K */ \
/*	... */ \
/*	1534K */ \
\
/* Do 1024 three_complex_last_unfft macros */ \
/*	distance between fft data elements is 1024K */ \
/*	do 256 macros each processing 24 data values */ \
\
start_timer(20); \
rdi = (uintptr_t)g->plus1_premults;/* Address of premultiplier table */ \
/* 32 iters of 2 iters of 4*2*clm */ \
for(unsigned int loopA = 32; loopA; loopA--) { \
	for(unsigned int loopB = 2; loopB; loopB--) { \
		xtouch(xptr(rcx+(clm-1)*128));	/* Preload TLB */ \
		xloop4edi_clm##clm(c1b, clmblkdst, 192, x3cl_three_complex_last_unfft(rsi, 64, 32*clmblkdst8)); /* 2 * 2*clm iterations */ \
	} \
	rsi += -4*2*clmblkdst+clmblkdst8;/* Next source pointer */ \
} \
end_timer(20); \
\
/* Normalize these values */ \
\
pass1_normalize(1, yes_fft); \
\
/* Copy data from the scratch area back to the FFT data area */ \
\
copy_scratch_data_to_fft(768, clm); \
goto	endlp;			/* Skip FFT code if flag not set */ \
\
/* Fire up auxiliary threads to do pass 1 forward FFT */ \
\
/* pass1: */	pass1_forward_fft_setup(b0b); \
\
/* Copy data from the FFT data area back to the scratch area */ \
\
b0b:	copy_fft_data_to_scratch(768, clm); \
\
/* Do FFT levels 1,2,3 */ \
/* */ \
/* On input the 64-byte cache lines hold these data values: */ \
/*	0K	1536K	3072K	4608K */ \
/*	... */ \
/*	1534K */ \
/* On output the 64-byte cache lines hold these data values: */ \
/*	0K	512K	1024K	1536K */ \
/*	2K	... */ \
/*	... */ \
/*	510K */ \
/*	2048K	3072K	4096K	5120K */ \
/*	2050K */ \
/*	... */ \
/*	3070K */ \
\
/* Do 1024 three_complex_first_fft macros (FFT levels 1,2,3) */ \
/*	distance between fft data elements is 1024K */ \
/*	do 256 iterations each processing 24 data values */ \
\
yes_fft:start_timer(22); \
rsi = (uintptr_t)g->scratch_area;	/* Get address of scratch area */ \
rdi = (uintptr_t)g->plus1_premults;/* Address of premultiplier table */ \
rcx = (uintptr_t)g->data_addr; /* Load source address for prefetching */ \
pfing(rbp = (uintptr_t)g->data_prefetch);/* Set prefetch distance */ \
pfing(rbp -= rcx); \
pfing(rcx += 400*blkdst);		/* Create prefetching pointer */ \
/* 32 iters of 2 iters of 4*2*clm */ \
for(unsigned int loopA = 32; loopA; loopA--) { \
	xtouch(xptr(rcx+(clm-1)*128));	/* Preload TLB */ \
	xloop8edi_clm##clm(b1b, clmblkdst, 192, x3cl_three_complex_fft(rsi, 64, 32*clmblkdst8));	/* 4 * 2*clm iterations */ \
	rsi += -8*clmblkdst+clmblkdst8; \
} \
rsi += -32*clmblkdst8;	/* Restore source pointer */ \
end_timer(22); \
\
/* Do FFT levels 4 */ \
/* */ \
/* On input the 64-byte cache lines hold these data values: */ \
/*	0K	512K	1024K	1536K */ \
/*	2K	... */ \
/*	... */ \
/*	510K */ \
/*	2048K	3072K	4096K	5120K */ \
/*	2050K */ \
/*	... */ \
/*	3070K */ \
/* On output the 64-byte cache lines hold these data values: */ \
/*	0K	512K	1024K	1536K */ \
/*	2K	... */ \
/*	... */ \
/*	510K	... */ \
/*	2048K	... */ \
/*	... */ \
\
/* Do 768 two_two_complex_fft macros */ \
/*	distance between fft data elements is 256K */ \
/*	do 384 macros each processing 16 data values */ \
\
start_timer(23); \
rdi = (uintptr_t)g->sincos5;/* Load sin/cos pointer */ \
/* 16 iters of 2 of 4*2*clm */ \
for(unsigned int loopA = 16; loopA; loopA--) { \
	xtouch(xptr(rcx+(clm-1)*128));	/* Preload TLB */ \
	xloop8_clm##clm(b2b, clmblkdst, x2cl_two_complex_fft_in_place(rsi, 64, 16*clmblkdst8, rdi));		/* 4 * 2*clm iterations */ \
	rsi += -8*clmblkdst+clmblkdst8; \
} \
rsi += -16*clmblkdst8+32*clmblkdst8;/* Next source pointer */ \
rdi += XMM_SCD;		/* Next sine/cosine pointer */ \
/* 32 iters of 2 of 4*2*clm */ \
for(unsigned int loopA = 32; loopA; loopA--) { \
	xtouch(xptr(rcx+(clm-1)*128));	/* Preload TLB */ \
	xloop8_clm##clm(b3b, clmblkdst, x2cl_two_complex_fft(rsi, 64, 32*clmblkdst8, rdi));		/* 4 * 2*clm iterations */ \
	rsi += -8*clmblkdst+clmblkdst8; \
} \
rsi += -32*clmblkdst8-32*clmblkdst8;/* Restore source pointer */ \
end_timer(23); \
\
/* Do FFT levels 5,6 */ \
/* */ \
/* On input the 64-byte cache lines hold these data values: */ \
/*	0K	512K	1024K	1536K */ \
/*	2K	... */ \
/*	... */ \
/*	510K	... */ \
/*	2048K	... */ \
/*	... */ \
/* On output the 64-byte cache lines hold these data values: */ \
/*	0K	128K	256K	384K */ \
/*	2K	... */ \
/*	... */ \
/*	126K	... */ \
/*	512K	... */ \
/*	... */ \
\
/* Do 768 four_complex_fft macros */ \
/*	distance between fft data elements is 128K */ \
/*	do 192 macros each processing 32 data values */ \
\
start_timer(24); \
rdi = (uintptr_t)g->sincos4;/* Load sin/cos pointer */ \
/* 3 iters of 8 iters of 4 of 2*2*clm */ \
for(unsigned int loopA = 3; loopA; loopA--) { \
	for(unsigned int loopB = 8; loopB; loopB--) { \
		for(unsigned int loopC = 4; loopC; loopC--) { \
			xtouch(xptr(rcx+(clm-1)*128));	/* Preload TLB */ \
			xloop2_clm##clm(b4b, clmblkdst, x4cl_four_complex_fft(rsi, 64, 8*clmblkdst8, 16*clmblkdst8, rdi));		/* 2 * 2*clm iterations */ \
		} \
		rsi += -4*2*clmblkdst+clmblkdst8; \
	} \
	rsi += -8*clmblkdst8+32*clmblkdst8;/* Next source pointer */ \
	rdi += 2*XMM_SCD;		/* Next sine/cosine pointer */ \
} \
rsi += -3*32*clmblkdst8;	/* Restore source pointer */ \
end_timer(24); \
\
/* Do FFT levels 7,8 */ \
/* */ \
/* On input the 64-byte cache lines hold these data values: */ \
/*	0K	128K	256K	384K */ \
/*	2K	... */ \
/*	... */ \
/*	126K	... */ \
/*	512K	... */ \
/*	... */ \
/* On output the 64-byte cache lines hold these data values: */ \
/*	0K	32K	64K	96K */ \
/*	2K	... */ \
/*	... */ \
/*	30K	... */ \
/*	128K	... */ \
/*	... */ \
\
/* Do 768 four_complex_fft macros */ \
/*	distance between fft data elements is 32K */ \
/*	do 192 macros each processing 32 data values */ \
\
start_timer(24); \
rdi = (uintptr_t)g->sincos3;/* Load sin/cos pointer */ \
/* 12 iters of 2 iters of 4 of 2*2*clm */ \
for(unsigned int loopA = 12; loopA; loopA--) { \
	for(unsigned int loopB = 2; loopB; loopB--) { \
		for(unsigned int loopC = 4; loopC; loopC--) { \
			xtouch(xptr(rcx+(clm-1)*128));	/* Preload TLB */ \
			xloop2_clm##clm(b5b, clmblkdst, x4cl_four_complex_fft(rsi, 64, 2*clmblkdst8, 4*clmblkdst8, rdi));		/* 2 * 2*clm iterations */ \
		} \
		rsi += -4*2*clmblkdst+clmblkdst8; \
	} \
	rsi += -2*clmblkdst8+8*clmblkdst8;/* Next source pointer */ \
	rdi += 2*XMM_SCD;		/* Next sine/cosine pointer */ \
} \
rsi += -12*8*clmblkdst8;	/* Restore source pointer */ \
end_timer(24); \
\
/* Do FFT levels 9,10 */ \
/* */ \
/* On input the 64-byte cache lines hold these data values: */ \
/*	0K	32K	64K	96K */ \
/*	2K	... */ \
/*	... */ \
/*	30K	... */ \
/*	128K	... */ \
/*	... */ \
/* On output the 64-byte cache lines hold these data values: */ \
/*	0K	8K	16K	24K */ \
/*	2K	... */ \
/*	... */ \
/*	6K	... */ \
/*	32K	... */ \
/*	... */ \
\
/* Do 768 four_complex_fft macros */ \
/*	distance between fft data elements is 8K */ \
/*	do 192 macros each processing 32 data values */ \
\
start_timer(25); \
rdi = (uintptr_t)g->sincos2;/* Load sin/cos pointer */ \
/* 48 iterations of 2 iters of 2*2*clm */ \
for(unsigned int loopA = 48; loopA; loopA--) { \
	for(unsigned int loopB = 2; loopB; loopB--) { \
		xtouch(xptr(rcx+(clm-1)*128));	/* Preload TLB */ \
		xloop2_clm##clm(b7b, clmblkdst, x4cl_four_complex_fft(rsi, 64, 4*clmblkdst, clmblkdst8, rdi));		/* 2*2*clm iterations */ \
	} \
	rsi += -2*2*clmblkdst+2*clmblkdst8;/* Next source pointer */ \
	rdi += 2*XMM_SCD;		/* Next sine/cosine pointer */ \
} \
rsi += -48*2*clmblkdst8;	/* Restore source pointer */ \
pfing(rcx += -768*blkdst);	/* Restore prefetch pointer */ \
end_timer(25); \
\
/* Do FFT levels 11,12 */ \
/* */ \
/* On input the 64-byte cache lines hold these data values: */ \
/*	0K	8K	16K	24K */ \
/*	2K	... */ \
/*	... */ \
/*	6K	... */ \
/*	32K	... */ \
/*	... */ \
/* On output the 64-byte cache lines hold these data values: */ \
/*	0K	2K	4K	6K */ \
/*	8K	... */ \
/*	... */ \
\
/* Do 768 four_complex_fft macros */ \
/*	distance between fft data elements is 2K */ \
/*	do 192 macros each processing 32 data values */ \
\
start_timer(26); \
rdi = (uintptr_t)g->sincos1;/* Load sin/cos pointer */ \
/* 96 iterations of 2 iters of 2*clm */ \
for(unsigned int loopA = 96; loopA; loopA--) { \
	for(unsigned int loopB = 2; loopB; loopB--) { \
		for(unsigned int loopC = clm; loopC; loopC--) { \
			for(unsigned int loopD = 2; loopD; loopD--) { \
				g4cl_four_complex_fft(rsi, 64, clmblkdst, 2*clmblkdst, rcx, 64, blkdst, 2*blkdst); \
				L2prefetch128(u8ptr(rcx-64+rbp)); \
				L2prefetch128(u8ptr(rcx-64+blkdst+rbp)); \
				pfing(rbp += -64+2*blkdst); \
			} \
			pfing(rbp += 128-4*blkdst); \
		} \
	rsi += -2*clm*64+4*clmblkdst;/* Next source pointer */ \
		rcx += -2*clm*64+4*blkdst;	/* Next dest pointer */ \
		rdi += 2*XMM_SCD;		/* Next sine/cosine pointer */ \
	} \
	rsi += -8*clmblkdst+clmblkdst8; \
} \
end_timer(26); \
\
/* Work on next 4 sets of 768 values */ \
\
endlp:	xpass1_get_next_block(pass2, c0b, b0b); /* Get next block and jump */ \
/* done: */	xfft_footer_hg2; \
\
\
\
/* Perform a 7168K-element FFT.  This is done in two passes.  Pass 1 */ \
/* does 12 levels, pass 2 uses common code to perform the last 11 levels. */ \
\

#define hg_pass1sclevels12pfa7(procname, xpass2name, pass2_levels, clm) \
 /* (fftname, pass2_levels, clmarg) */ \
\
/* set_FFT_constants(fftname, 0, 0, pass2_levels, clmarg); */ \
xfft_header_hg2(procname, xpass2name, pass2_levels, clm); \
\
/* Do the last 11 FFT levels */ \
/*	do 1 pass2_11_levels_real and 223 pass2_11_levels_complex */ \
\
/* pass2: */	start_timer(0); \
xpass2name(g); \
end_timer(0); \
if(g->ffttype == 1)		/* We're done if FFTing only */ \
	return; \
\
/* Do the inverse FFT */ \
\
pass1_inverse_fft_setup(c0b, procname); \
\
/* Do inverse FFT levels 11,12 */ \
/* On input the 64-byte cache lines hold these data values: */ \
/*	0K	2K	4K	6K */ \
/*	8K	10K	12K	14K */ \
/*	... */ \
/* On output the 64-byte cache lines hold these data values: */ \
/*	0K	4K	8K	12K */ \
/*	2K	... */ \
/*	16K	... */ \
/*	... */ \
\
/* Do 1 eight_reals_unfft_2 macros */ \
/* Do 1 nop_two_two_complex_unfft_2 macros */ \
/*	distance between fft data elements is 2K */ \
/*	do 1 macros each processing 16 data values */ \
\
c0b:	start_timer(16); \
rsi = (uintptr_t)g->data_addr;/* Get source address */ \
rcx = (uintptr_t)g->scratch_area;	/* Get address of scratch area */ \
pfing(rbp = (uintptr_t)g->data_prefetch);/* Set prefetch distance */ \
pfing(rbp -= rsi); \
for(unsigned int loopA = 2*clm; loopA; loopA--) { \
	g2cl_eight_reals_unfft_2(rsi, 64, blkdst, rcx, 64, clmblkdst); \
} \
\
/* Do 894 four_complex_unfft macros */ \
/*	distance between fft data elements is 2K */ \
/*	do 447 macros each processing 16 data values */ \
\
rsi += -2*clm*64+2*blkdst;	/* Next source pointer */ \
rcx += -2*clm*64+2*clmblkdst;/* Next dest pointer */ \
rdi = (uintptr_t)g->sincos1;/* Load sin/cos pointer */ \
/* eax = 112*65536+256/4; */ /* 447 iterations of 2*clm */ \
for(unsigned int loopA = 447; loopA; ) { \
	do{ \
		for(unsigned int loopC = 2*clm; loopC; loopC--) { \
			g2cl_four_complex_unfft(rsi, 64, blkdst, rcx, 64, clmblkdst); \
		} \
		rsi += -2*clm*64+2*blkdst;	/* Next source pointer */ \
		rcx += -2*clm*64+2*clmblkdst;/* Next dest pointer */ \
		rdi += XMM_SCD;		/* Next sine/cosine pointer */ \
	}while(--loopA & 3); \
	rcx += -4*2*clmblkdst+clmblkdst8; \
} \
rcx = rsi-448*2*blkdst;	/* Create prefetch pointer */ \
rsi = (uintptr_t)g->scratch_area;	/* Restore source pointer */ \
end_timer(16); \
\
/* Do inverse FFT levels 9,10 */ \
/* On input the 64-byte cache lines hold these data values: */ \
/*	0K	4K	8K	12K */ \
/*	2K	... */ \
/*	16K	... */ \
/*	... */ \
/* On output the 64-byte cache lines hold these data values: */ \
/*	0K	16K	32K	48K */ \
/*	2K	... */ \
/*	... */ \
/*	14K */ \
/*	64K	... */ \
/*	... */ \
\
/* Do 4 eight_reals_unfft_2 macros */ \
/* Do 4 nop_two_two_complex_unfft_2 macros */ \
/*	distance between fft data elements is 8K */ \
/*	do 2 macros each processing 32 data values */ \
\
start_timer(17); \
/* 2 iterations of 2*clm */ \
for(unsigned int loopA = 2; loopA; loopA--) { \
	for(unsigned int loopB = 2*clm; loopB; loopB--) { \
		x4cl_eight_reals_unfft_2(rsi, 64, 2*clmblkdst, 4*clmblkdst); \
	} \
	rsi += -clm*2*64+clmblkdst;/* Next source pointer */ \
} \
\
/* Do 888 four_complex_unfft macros */ \
/*	distance between fft data elements is 8K */ \
/*	do 222 macros each processing 32 data values */ \
\
rsi += -2*clmblkdst+clmblkdst8;/* Next source pointer */ \
rdi = (uintptr_t)g->sincos2;/* Load sin/cos pointer */ \
/* 111 iterations of 2*2*clm */ \
for(unsigned int loopA = 111; loopA; loopA--) { \
	xtouch(xptr(rcx+(clm-1)*128));	/* Preload TLB */ \
	xloop2_clm##clm(c7b, clmblkdst, x4cl_four_complex_unfft(rsi, 64, 2*clmblkdst, 4*clmblkdst, rdi)); \
	rsi += -2*clmblkdst+clmblkdst8;/* Next source pointer */ \
	rdi += XMM_SCD;		/* Next sine/cosine pointer */ \
} \
rsi += -112*clmblkdst8;	/* Restore source pointer */ \
end_timer(17); \
\
/* Do inverse FFT levels 7,8 */ \
/* On input the 64-byte cache lines hold these data values: */ \
/*	0K	16K	32K	48K */ \
/*	2K	... */ \
/*	... */ \
/*	14K */ \
/*	64K	... */ \
/*	... */ \
/* On output the 64-byte cache lines hold these data values: */ \
/*	0K	64K	128K	192K */ \
/*	2K	... */ \
/*	... */ \
/*	62K */ \
/*	256K	... */ \
/*	... */ \
\
/* Do 16 eight_reals_unfft_2 macros */ \
/* Do 16 nop_two_two_complex_unfft_2 macros */ \
/*	distance between fft data elements is 32K */ \
/*	do 8 macros each processing 32 data values */ \
\
start_timer(18); \
/* 8 iterations of 2*clm */ \
for(unsigned int loopA = 8; loopA; loopA--) { \
	for(unsigned int loopB = 2*clm; loopB; loopB--) { \
		x4cl_eight_reals_unfft_2(rsi, 64, clmblkdst8, 2*clmblkdst8); \
	} \
	rsi += -clm*2*64+clmblkdst;/* Next source pointer */ \
} \
\
/* Do 864 four_complex_unfft macros */ \
/*	distance between fft data elements is 32K */ \
/*	do 216 macros each processing 32 data values */ \
\
rsi += -8*clmblkdst+4*clmblkdst8;/* Next source pointer */ \
rdi = (uintptr_t)g->sincos3;/* Load sin/cos pointer */ \
/* 27 iterations of 4 iters of 2*2*clm */ \
for(unsigned int loopA = 27; loopA; loopA--) { \
	for(unsigned int loopB = 4; loopB; loopB--) { \
		xtouch(xptr(rcx+(clm-1)*128));	/* Preload TLB */ \
		xloop2_clm##clm(c5d, clmblkdst, x4cl_four_complex_unfft(rsi, 64, clmblkdst8, 2*clmblkdst8, rdi)); \
	} \
	rsi += -4*2*clmblkdst+4*clmblkdst8;/* Next source pointer */ \
	rdi += XMM_SCD;		/* Next sine/cosine pointer */ \
} \
rsi += -28*4*clmblkdst8;	/* Restore source pointer */ \
end_timer(18); \
\
/* Do inverse FFT levels 6 */ \
/* On input the 64-byte cache lines hold these data values: */ \
/*	0K	64K	128K	192K */ \
/*	2K	... */ \
/*	... */ \
/*	62K */ \
/*	256K	... */ \
/*	... */ \
/* On output the 64-byte cache lines hold these data values: */ \
/*	0K	128K	256K	384K */ \
/*	2K	... */ \
/*	... */ \
/*	126K */ \
/*	512K	... */ \
/*	... */ \
\
/* Do 64 eight_reals_unfft_1 macros */ \
/*	distance between fft data elements is 128K */ \
/*	do 32 macros each processing 16 data values */ \
/* */ \
/* Do 832 four_complex_unfft macros */ \
/*	distance between fft data elements is 128K */ \
/*	do 416 macros each processing 16 data values */ \
\
/* Do inverse FFT levels 4,5 */ \
/* On input the 64-byte cache lines hold these data values: */ \
/*	0K	128K	256K	384K */ \
/*	2K	... */ \
/*	... */ \
/*	126K */ \
/*	512K	... */ \
/*	... */ \
/* On output the 64-byte cache lines hold these data values: */ \
/*	0K	256K	512K	768K */ \
/*	2K	... */ \
/*	... */ \
/*	254K */ \
/*	1024K	1536K	2048K	2560K */ \
/*	... */ \
/* */ \
/* Do 128 eight_reals_unfft_2 macros */ \
/*	distance between fft data elements is 256K */ \
/*	do 64 macros each processing 16 data values */ \
/* */ \
/* Do 768 four_complex_unfft macros */ \
/*	distance between fft data elements is 256K */ \
/*	do 192 macros each processing 32 data values */ \
/* */ \
/* Do inverse FFT levels 1,2,3 */ \
/* On input the 64-byte cache lines hold these data values: */ \
/*	0K	256K	512K	768K */ \
/*	2K	... */ \
/*	... */ \
/*	254K */ \
/*	1024K	1536K	2048K	2560K */ \
/*	... */ \
/* On output the 64-byte cache lines hold these data values: */ \
/*	0K	256K	3584K	3840K */ \
/*	2K	... */ \
/*	... */ \
/*	254K */ \
/*	512K	1024K	4096K	4608K */ \
/*	514K */ \
/*	... */ \
/*	1022K */ \
/*	1536K	2560K	5120K	6144K */ \
/*	1538K */ \
/*	... */ \
/*	2558K */ \
/* */ \
/* Do 1024 seven_reals_last_unfft macros */ \
/*	distance between fft data elements is 1024K */ \
/*	do 128 macros each processing 56 data values */ \
\
start_timer(19); \
rdi = (uintptr_t)g->sincos4;/* Load sin/cos pointer */ \
rbx = (uintptr_t)g->sincos5;/* Load sin/cos pointer */ \
/* 4 iterations */ \
xseven_reals_unfft_preload; \
\
for(unsigned int loopA = 4; loopA; loopA--) { \
	for(unsigned int loopB = 2*4; loopB; loopB--) { \
		xtouchlines(7);			/* Preload TLBs */ \
		\
		for(unsigned int loopC = clm; loopC; loopC--) { \
			for(unsigned int loopD = 2; loopD; loopD--) { \
				x2cl_eight_reals_unfft_1(rsi, 64, 4*clmblkdst8); \
			} \
			\
			rsi += -2*64+8*clmblkdst8;	/* Next source pointer */ \
			/* 13 iterations */ \
			for(unsigned int loopD = 13; loopD; loopD--) { \
				for(unsigned int loopE = 2; loopE; loopE--) { \
					x2cl_two_complex_unfft(rsi, 64, 4*clmblkdst8); \
				} \
				rsi += -2*64+8*clmblkdst8;	/* Next source pointer */ \
				rdi += XMM_SCD;		/* Next sine/cosine pointer */ \
			} \
			rdi += -13*XMM_SCD;	/* Restore sine/cosine pointer */ \
			\
			rsi += -14*8*clmblkdst8;	/* Restore source pointer */ \
			for(unsigned int loopD = 2; loopD; loopD--) { \
				for(unsigned int loopE = 2; loopE; loopE--) { \
					x2cl_half_eight_reals_unfft_2(rsi, 64, 8*clmblkdst8); \
				} \
				rsi += -2*64+4*clmblkdst8; \
			} \
			\
			rsi += -2*4*clmblkdst8+16*clmblkdst8;/* Next source pointer */ \
			/* 3 iterations */ \
			for(unsigned int loopD = 3; loopD; loopD--) { \
				xloop2sp(c3b, 4*clmblkdst8, x4cl_four_complex_unfft(rsi, 64, 8*clmblkdst8, 16*clmblkdst8, rbx)); \
				rsi += -2*4*clmblkdst8+32*clmblkdst8; \
				rbx += XMM_SCD;		/* Next sine/cosine pointer */ \
			} \
			rbx += -3*XMM_SCD;		/* Restore sine/cosine pointer */ \
			\
			rsi += -3*32*clmblkdst8-16*clmblkdst8;/* Restore src pointer */ \
			for(unsigned int loopD = 4; loopD; loopD--) { \
				xloop1sp(c1b, 4*clmblkdst8, x7cl_seven_reals_last_unfft(rsi, 64, 16*clmblkdst8)); \
			} \
			\
				EXEC248##clm(rsi += -4*4*clmblkdst8+2*64);/* Next source pointer */ \
			EXEC248##clm(pfing(rcx += -7*blkdst+128));	/* Next prefetch pointer */ \
		} \
		EXEC248##clm(rsi += -clm*2*64+clmblkdst);/* Next source pointer */ \
		EXEC248##clm(pfing(rcx += -clm*128+7*blkdst));	/* Next prefetch pointer */ \
		EXEC1##clm(rsi += -4*4*clmblkdst8+clmblkdst);/* Next source pointer */ \
		\
	} \
	rsi += -8*clmblkdst+clmblkdst8;/* Next source pointer */ \
} \
end_timer(19); \
\
/* Normalize these values */ \
\
pass1_normalize(1, yes_fft); \
\
/* Copy data from the scratch area back to the FFT data area */ \
\
copy_scratch_data_to_fft(896, clm); \
goto	endlp;			/* Skip FFT code if flag not set */ \
\
/* Fire up auxiliary threads to do pass 1 forward FFT */ \
\
/* pass1: */	pass1_forward_fft_setup(b0b); \
\
/* Copy data from the FFT data area back to the scratch area */ \
\
b0b:	copy_fft_data_to_scratch(896, clm); \
\
/* Do FFT levels 1,2,3 */ \
/* */ \
/* On input the 64-byte cache lines hold these data values: */ \
/*	0K	256K	3584K	3840K */ \
/*	2K	... */ \
/*	... */ \
/*	254K */ \
/*	512K	1024K	4096K	4608K */ \
/*	514K */ \
/*	... */ \
/*	1022K */ \
/*	1536K	2560K	5120K	6144K */ \
/*	1538K */ \
/*	... */ \
/*	2558K */ \
/* On output the 64-byte cache lines hold these data values: */ \
/*	0K	256K	512K	768K */ \
/*	2K	... */ \
/*	... */ \
/*	254K */ \
/*	1024K	1536K	2048K	2560K */ \
/*	1026K */ \
/*	... */ \
/*	1534K */ \
/*	3072K	4096K	5120K	6144K */ \
/*	3074K */ \
/*	... */ \
/*	4094K */ \
/* */ \
/* Do 1024 seven_reals_first_fft macros (FFT levels 1,2,3) */ \
/*	distance between fft data elements is 1024K */ \
/*	do 128 iterations each processing 56 data values */ \
\
/* Do FFT levels 4,5 */ \
/* Values 0K-1023K is real data, 1024K-2047K is semi-real data, 2048K-7167K is */ \
/* complex data. */ \
/* */ \
/* On input the 64-byte cache lines hold these data values: */ \
/*	0K	256K	512K	768K */ \
/*	2K	... */ \
/*	... */ \
/*	254K */ \
/*	1024K	1536K	2048K	2560K */ \
/*	1026K */ \
/*	... */ \
/*	1534K */ \
/*	3072K	4096K	5120K	6144K */ \
/*	3074K */ \
/*	... */ \
/*	4094K */ \
/* On output the 64-byte cache lines hold these data values: */ \
/*	0K	256K	512K	768K */ \
/*	2K	... */ \
/*	... */ \
/*	254K	... */ \
/*	1024K	... */ \
/*	... */ \
/* */ \
/* Do 128 eight_reals_fft_2 macros */ \
/*	distance between fft data elements is 256K */ \
/*	do 128 macros each processing 8 data values */ \
/* */ \
/* Do 768 four_complex_fft macros */ \
/*	distance between fft data elements is 256K */ \
/*	do 192 macros each processing 32 data values */ \
\
/* Do FFT levels 6 */ \
/* Values 0K-255K is real data, 256K-511K is semi-real data, 512K-7167K is */ \
/* complex data. */ \
/* */ \
/* On input the 64-byte cache lines hold these data values: */ \
/*	0K	256K	512K	768K */ \
/*	2K	... */ \
/*	... */ \
/*	254K	... */ \
/*	1024K	... */ \
/*	... */ \
/* On output the 64-byte cache lines hold these data values: */ \
/*	0K	128K	256K	384K */ \
/*	2K	... */ \
/*	... */ \
/*	126K	... */ \
/*	512K	... */ \
/*	... */ \
/* */ \
/* Do 64 eight_reals_fft_1 macros */ \
/* Do 64 two_two_complex_fft macros */ \
/*	distance between fft data elements is 128K */ \
/*	do 64 macros each processing 16 data values */ \
/* */ \
/* Do 768 two_two_complex_fft macros */ \
/*	distance between fft data elements is 128K */ \
/*	do 384 macros each processing 16 data values */ \
\
yes_fft:start_timer(22); \
rsi = (uintptr_t)g->scratch_area;	/* Get address of scratch area */ \
rcx = (uintptr_t)g->data_addr; /* Load source address for prefetching */ \
pfing(rbp = (uintptr_t)g->data_prefetch);/* Set prefetch distance */ \
pfing(rbp -= rcx); \
pfing(rcx += 426*blkdst);		/* Create prefetching pointer */ \
rdi = (uintptr_t)g->sincos5;/* Load sin/cos pointer */ \
rbx = (uintptr_t)g->sincos4;/* Load sin/cos pointer */ \
/* 8 iterations */ \
xseven_reals_fft_preload; \
\
for(unsigned int loopA = 8; loopA; loopA--) { \
	for(unsigned int loopB = 2*4; loopB; loopB--) { \
		xtouchlines(4);			/* Preload TLBs */ \
		\
		for(unsigned int loopC = clm; loopC; loopC--) { \
			for(unsigned int loopD = 2; loopD; loopD--) { \
				xloop1sp(b1b, 8*clmblkdst8, x7cl_seven_reals_fft(rsi, 64, 16*clmblkdst8)); \
			} \
			\
			rsi += -2*8*clmblkdst8;	/* Restore source pointer */ \
			for(unsigned int loopD = 2; loopD; loopD--) { \
				for(unsigned int loopE = 2; loopE; loopE--) { \
					x1cl_half_eight_reals_fft_2(rsi, 64); \
				} \
				rsi += -2*64+8*clmblkdst8;	/* Next source pointer */ \
			} \
			\
			rsi += -2*8*clmblkdst8+16*clmblkdst8;/* Next source pointer */ \
			for(unsigned int loopD = 2; loopD; loopD--) { \
				for(unsigned int loopE = 2; loopE; loopE--) { \
					x2cl_four_complex_fft(rsi, 64, 16*clmblkdst8); \
				} \
				rsi += -2*64+8*clmblkdst8;	/* Next source pointer */ \
			} \
			\
			rsi += -2*8*clmblkdst8+32*clmblkdst8;/* Next source pointer */ \
			rdi += XMM_SCD;		/* Next sine/cosine pointer */ \
			for(unsigned int loopD = 2; loopD; loopD--) { \
				xloop1sp(b3b, 8*clmblkdst8, x4cl_four_complex_fft(rsi, 64, 16*clmblkdst8, 32*clmblkdst8, rdi)); \
			} \
			rdi += -XMM_SCD;		/* Restore sine/cosine pointer */ \
			\
			rsi += -2*8*clmblkdst8-48*clmblkdst8;/* Restore src pointer */ \
			for(unsigned int loopD = 2; loopD; loopD--) { \
				x2cl_eight_reals_fft_1(rsi, 64, 8*clmblkdst8); \
			} \
			\
			rsi += -2*64+16*clmblkdst8;/* Next source pointer */ \
			rbx += XMM_SCD;		/* Next sine/cosine pointer */ \
			/* 6 iterations */ \
			for(unsigned int loopD = 6; loopD; loopD--) { \
				for(unsigned int loopE = 2; loopE; loopE--) { \
					x2cl_two_complex_fft(rsi, 64, 8*clmblkdst8, rbx); \
				} \
				rsi += -2*64+16*clmblkdst8;/* Next source pointer */ \
				rbx += 2*XMM_SCD;		/* Next sine/cosine pointer */ \
			} \
			rbx += -13*XMM_SCD;	/* Restore sine/cosine pointer */ \
			\
				EXEC248##clm(rsi += -7*16*clmblkdst8+2*64);/* Next source pointer */ \
			EXEC248##clm(pfing(rcx += -4*blkdst+128));	/* Next prefetch pointer */ \
		} \
		EXEC248##clm(rsi += -clm*2*64+clmblkdst);/* Next source pointer */ \
		EXEC248##clm(pfing(rcx += -clm*128+4*blkdst));	/* Next prefetch pointer */ \
		EXEC1##clm(rsi += -7*16*clmblkdst8+clmblkdst);/* Next source pointer */ \
		\
	} \
	rsi += -8*clmblkdst+clmblkdst8;/* Next source pointer */ \
} \
rsi += -8*clmblkdst8;	/* Restore source pointer */ \
end_timer(22); \
\
/* Do FFT levels 7,8 */ \
/* Values 0K-127K is real data, 128K-255K is semi-real data, 256K-7167K is */ \
/* complex data. */ \
/* */ \
/* On input the 64-byte cache lines hold these data values: */ \
/*	0K	128K	256K	384K */ \
/*	2K	... */ \
/*	... */ \
/*	126K	... */ \
/*	512K	... */ \
/*	... */ \
/* On output the 64-byte cache lines hold these data values: */ \
/*	0K	32K	64K	96K */ \
/*	2K	... */ \
/*	... */ \
/*	30K	... */ \
/*	128K	... */ \
/*	... */ \
\
/* Do 32 eight_reals_fft_2 macros */ \
/* Do 32 two_two_complex_fft macros */ \
/*	distance between fft data elements is 32K */ \
/*	do 16 macros each processing 32 data values */ \
\
start_timer(24); \
rdi = (uintptr_t)g->sincos3;/* Load sin/cos pointer */ \
/* 2 iterations of 8 iters of 2*clm */ \
for(unsigned int loopA = 2; loopA; loopA--) { \
	for(unsigned int loopB = 8; loopB; loopB--) { \
		for(unsigned int loopC = 2*clm; loopC; loopC--) { \
			x4cl_eight_reals_fft_2(rsi, 64, 2*clmblkdst8, 4*clmblkdst8, rdi); \
		} \
		rsi += -2*clm*64+clmblkdst;/* Next source pointer */ \
	} \
	rsi += -8*clmblkdst+clmblkdst8;/* Next source pointer */ \
} \
\
/* Do 832 four_complex_fft macros */ \
/*	distance between fft data elements is 32K */ \
/*	do 208 macros each processing 32 data values */ \
\
rsi += -2*clmblkdst8+8*clmblkdst8;/* Load source pointer */ \
rdi += XMM_SCD;		/* Next sine/cosine pointer */ \
/* 13 iters of 2*4 iters of 2*2*clm */ \
for(unsigned int loopA = 13; loopA; loopA--) { \
	for(unsigned int loopB = 2; loopB; loopB--) { \
		for(unsigned int loopC = 4; loopC; loopC--) { \
			xtouch(xptr(rcx+(clm-1)*128));	/* Preload TLB */ \
			xloop2_clm##clm(b5d, clmblkdst, x4cl_four_complex_fft(rsi, 64, 2*clmblkdst8, 4*clmblkdst8, rdi)); \
		} \
		rsi += -8*clmblkdst+clmblkdst8;/* Next source pointer */ \
	} \
	rsi += -2*clmblkdst8+8*clmblkdst8;/* Next source pointer */ \
	rdi += 2*XMM_SCD;		/* Next sine/cosine pointer */ \
} \
rsi += -14*8*clmblkdst8;	/* Restore source pointer */ \
end_timer(24); \
\
/* Do FFT levels 9,10 */ \
/* Values 0K-31K is real data, 32K-63K is semi-real data, 64K-7167K is */ \
/* complex data. */ \
/* */ \
/* On input the 64-byte cache lines hold these data values: */ \
/*	0K	32K	64K	96K */ \
/*	2K	... */ \
/*	... */ \
/*	30K	... */ \
/*	128K	... */ \
/*	... */ \
/* On output the 64-byte cache lines hold these data values: */ \
/*	0K	8K	16K	24K */ \
/*	2K	... */ \
/*	... */ \
/*	6K	... */ \
/*	32K	... */ \
/*	... */ \
\
/* Do 8 eight_reals_fft_2 macros */ \
/* Do 8 two_two_complex_fft macros */ \
/*	distance between fft data elements is 8K */ \
/*	do 4 macros each processing 32 data values */ \
\
start_timer(25); \
rdi = (uintptr_t)g->sincos2;/* Load sin/cos pointer */ \
/* 4 iterations of 2*clm */ \
for(unsigned int loopA = 4; loopA; loopA--) { \
	for(unsigned int loopB = 2*clm; loopB; loopB--) { \
		x4cl_eight_reals_fft_2(rsi, 64, 4*clmblkdst, clmblkdst8, rdi); \
	} \
	rsi += -2*clm*64+clmblkdst;/* Next source pointer */ \
} \
\
/* Do 880 four_complex_fft macros */ \
/*	distance between fft data elements is 8K */ \
/*	do 220 macros each processing 32 data values */ \
\
rsi += -4*clmblkdst+2*clmblkdst8;/* Load source pointer */ \
rdi += XMM_SCD;		/* Next sine/cosine pointer */ \
/* 55 iterations of 4 iters of 2*clm */ \
for(unsigned int loopA = 55; loopA; loopA--) { \
	for(unsigned int loopB = 2; loopB; loopB--) { \
		xtouch(xptr(rcx+(clm-1)*128));	/* Preload TLB */ \
		xloop2_clm##clm(b7b, clmblkdst, x4cl_four_complex_fft(rsi, 64, 4*clmblkdst, clmblkdst8, rdi)); \
	} \
	rsi += -4*clmblkdst+2*clmblkdst8;/* Next source pointer */ \
	rdi += 2*XMM_SCD;		/* Next sine/cosine pointer */ \
} \
rsi += -56*2*clmblkdst8;	/* Restore source pointer */ \
pfing(rcx += -896*blkdst);	/* Restore prefetch pointer */ \
end_timer(25); \
\
/* Do FFT levels 11,12 */ \
/* Values 0K-7K is real data, 8K-15K is semi-real data, 16K-7167K is */ \
/* complex data. */ \
/* */ \
/* On input the 64-byte cache lines hold these data values: */ \
/*	0K	8K	16K	24K */ \
/*	2K	... */ \
/*	... */ \
/*	6K	... */ \
/*	32K	... */ \
/*	... */ \
/* On output the 64-byte cache lines hold these data values: */ \
/*	0K	2K	4K	6K */ \
/*	8K	... */ \
/*	... */ \
\
/* Do 1 eight_reals_fft_2 macros */ \
/* Do 1 nop_two_two_complex_fft_2 macros */ \
/* Do 2 four_complex_fft macros */ \
/*	distance between fft data elements is 2K */ \
/*	do 1 macros each processing 32 data values */ \
\
start_timer(26); \
rdi = (uintptr_t)g->sincos1;/* Load sin/cos pointer */ \
for(unsigned int loopA = 2*clm; loopA; loopA--) { \
	g4cl_eight_reals_fft_2(rsi, 64, clmblkdst, 2*clmblkdst, rcx, 64, blkdst, 2*blkdst); \
} \
\
/* Do 892 four_complex_fft macros */ \
/*	distance between fft data elements is 2K */ \
/*	do 223 macros each processing 32 data values */ \
\
rsi += -2*clm*64+4*clmblkdst;/* Next source pointer */ \
rcx += -2*clm*64+4*blkdst;	/* Next dest pointer */ \
rdi += XMM_SCD;		/* Next sine/cosine pointer */ \
/* eax = 112*65536+1; */ /* 223 iterations of 2*clm */ \
for(unsigned int loopA = 112; loopA; ) { \
	do{ \
		for(unsigned int loopC = clm; loopC; loopC--) { \
			for(unsigned int loopD = 2; loopD; loopD--) { \
				g4cl_four_complex_fft(rsi, 64, clmblkdst, 2*clmblkdst, rcx, 64, blkdst, 2*blkdst); \
				L2prefetch128(u8ptr(rcx-64+rbp)); \
				L2prefetch128(u8ptr(rcx-64+blkdst+rbp)); \
				pfing(rbp += -64+2*blkdst); \
			} \
			pfing(rbp += 128-4*blkdst); \
		} \
		rsi += -2*clm*64+4*clmblkdst;/* Next source pointer */ \
		rcx += -2*clm*64+4*blkdst;	/* Next dest pointer */ \
		rdi += 2*XMM_SCD;		/* Next sine/cosine pointer */ \
	}while(--loopA & 1); \
	rsi += -8*clmblkdst+clmblkdst8; \
} \
end_timer(26); \
\
/* Work on next 4 sets of 896 values */ \
\
endlp:	xpass1_get_next_block(pass2, c0b, b0b); /* Get next block and jump */ \
/* done: */	xfft_footer_hg2; \
\
\
\
/* Perform a 8M-element FFT.  This is done in two passes.  Pass 1 */ \
/* does 12 levels, pass 2 uses common code to perform the last 11 levels. */ \
\

#define hg_pass1sclevels12pfa8(procname, xpass2name, pass2_levels, clm) \
 /* (fftname, pass2_levels, clmarg) */ \
\
/* set_FFT_constants(fftname, 0, 0, pass2_levels, clmarg); */ \
xfft_header_hg2(procname, xpass2name, pass2_levels, clm); \
\
/* Do the last 11 FFT levels */ \
/*	do 1 pass2_11_levels_real and 255 pass2_11_levels_complex */ \
\
/* pass2: */	start_timer(0); \
xpass2name(g); \
end_timer(0); \
if(g->ffttype == 1)		/* We're done if FFTing only */ \
	return; \
\
/* Do the inverse FFT */ \
\
pass1_inverse_fft_setup(c0b, procname); \
\
/* Do inverse FFT levels 10,11 */ \
/* On input the 64-byte cache lines hold these data values: */ \
/*	0K	2K	4K	6K */ \
/*	8K	10K	12K	14K */ \
/*	... */ \
/* On output the 64-byte cache lines hold these data values: */ \
/*	0K	4K	8K	12K */ \
/*	2K	... */ \
/*	16K	... */ \
/*	... */ \
\
/* Do 1 eight_reals_unfft_2 macros */ \
/* Do 1 nop_two_two_complex_unfft_2 macros */ \
/*	distance between fft data elements is 2K */ \
/*	do 1 macros each processing 16 data values */ \
\
c0b:	start_timer(16); \
rsi = (uintptr_t)g->data_addr;/* Get source address */ \
rcx = (uintptr_t)g->scratch_area;	/* Get address of scratch area */ \
pfing(rbp = (uintptr_t)g->data_prefetch);/* Set prefetch distance */ \
pfing(rbp -= rsi); \
for(unsigned int loopA = 2*clm; loopA; loopA--) { \
	g2cl_eight_reals_unfft_2(rsi, 64, blkdst, rcx, 64, clmblkdst); \
} \
\
/* Do 1022 four_complex_unfft macros */ \
/*	distance between fft data elements is 2K */ \
/*	do 511 macros each processing 16 data values */ \
\
rsi += -2*clm*64+2*blkdst;	/* Next source pointer */ \
rcx += -2*clm*64+2*clmblkdst;/* Next dest pointer */ \
rdi = (uintptr_t)g->sincos1;/* Load sin/cos pointer */ \
/* eax = 128*65536+256/4; */ /* 511 iterations of 2*clm */ \
for(unsigned int loopA = 511; loopA; ) { \
	do{ \
		for(unsigned int loopC = 2*clm; loopC; loopC--) { \
			g2cl_four_complex_unfft(rsi, 64, blkdst, rcx, 64, clmblkdst); \
		} \
		rsi += -2*clm*64+2*blkdst;	/* Next source pointer */ \
		rcx += -2*clm*64+2*clmblkdst;/* Next dest pointer */ \
		rdi += XMM_SCD;		/* Next sine/cosine pointer */ \
	}while(--loopA & 3); \
	rcx += -8*clmblkdst+clmblkdst8; \
} \
rcx = rsi-512*2*blkdst;	/* Create prefetch pointer */ \
rsi = (uintptr_t)g->scratch_area;	/* Restore source pointer */ \
end_timer(16); \
\
/* Do inverse FFT levels 9,10 */ \
/* On input the 64-byte cache lines hold these data values: */ \
/*	0K	4K	8K	12K */ \
/*	2K	... */ \
/*	16K	... */ \
/*	... */ \
/* On output the 64-byte cache lines hold these data values: */ \
/*	0K	16K	32K	48K */ \
/*	2K	... */ \
/*	... */ \
/*	14K */ \
/*	64K	... */ \
/*	... */ \
\
/* Do 4 eight_reals_unfft_2 macros */ \
/* Do 4 nop_two_two_complex_unfft_2 macros */ \
/*	distance between fft data elements is 8K */ \
/*	do 2 macros each processing 32 data values */ \
\
start_timer(17); \
/* 2 iterations of 2*clm */ \
for(unsigned int loopA = 2; loopA; loopA--) { \
	for(unsigned int loopB = 2*clm; loopB; loopB--) { \
		x4cl_eight_reals_unfft_2(rsi, 64, 2*clmblkdst, 4*clmblkdst); \
	} \
	rsi += -2*clm*64+clmblkdst;/* Next source pointer */ \
} \
\
/* Do 1016 four_complex_unfft macros */ \
/*	distance between fft data elements is 8K */ \
/*	do 254 macros each processing 32 data values */ \
\
rsi += -2*clmblkdst+clmblkdst8;/* Next source pointer */ \
rdi = (uintptr_t)g->sincos2;/* Load sin/cos pointer */ \
/* 127 iterations of 2*2*clm */ \
for(unsigned int loopA = 127; loopA; loopA--) { \
	xtouch(xptr(rcx+(clm-1)*128));	/* Preload TLB */ \
	xloop2_clm##clm(c7b, clmblkdst, x4cl_four_complex_unfft(rsi, 64, 2*clmblkdst, 4*clmblkdst, rdi)); \
	rsi += -2*clmblkdst+clmblkdst8;/* Next source pointer */ \
	rdi += XMM_SCD;		/* Next sine/cosine pointer */ \
} \
rsi += -128*clmblkdst8;	/* Restore source pointer */ \
end_timer(17); \
\
/* Do inverse FFT levels 7,8 */ \
/* On input the 64-byte cache lines hold these data values: */ \
/*	0K	16K	32K	48K */ \
/*	2K	... */ \
/*	... */ \
/*	14K */ \
/*	64K	... */ \
/*	... */ \
/* On output the 64-byte cache lines hold these data values: */ \
/*	0K	64K	128K	192K */ \
/*	2K	... */ \
/*	... */ \
/*	62K */ \
/*	256K	... */ \
/*	... */ \
\
/* Do 16 eight_reals_unfft_2 macros */ \
/* Do 16 nop_two_two_complex_unfft_2 macros */ \
/*	distance between fft data elements is 32K */ \
/*	do 8 macros each processing 32 data values */ \
\
start_timer(18); \
/* 8 iterations of 2*clm */ \
for(unsigned int loopA = 8; loopA; loopA--) { \
	for(unsigned int loopB = 2*clm; loopB; loopB--) { \
		x4cl_eight_reals_unfft_2(rsi, 64, clmblkdst8, 2*clmblkdst8); \
	} \
	rsi += -2*clm*64+clmblkdst;/* Next source pointer */ \
} \
\
/* Do 992 four_complex_unfft macros */ \
/*	distance between fft data elements is 32K */ \
/*	do 248 macros each processing 32 data values */ \
\
rsi += -8*clmblkdst+4*clmblkdst8;/* Next source pointer */ \
rdi = (uintptr_t)g->sincos3;/* Load sin/cos pointer */ \
/* 31 iterations of 4 iters of 2*2*clm */ \
for(unsigned int loopA = 31; loopA; loopA--) { \
	for(unsigned int loopB = 4; loopB; loopB--) { \
		xtouch(xptr(rcx+(clm-1)*128));	/* Preload TLB */ \
		xloop2_clm##clm(c5d, clmblkdst, x4cl_four_complex_unfft(rsi, 64, clmblkdst8, 2*clmblkdst8, rdi)); \
	} \
	rsi += -4*2*clmblkdst+4*clmblkdst8;/* Next source pointer */ \
	rdi += XMM_SCD;		/* Next sine/cosine pointer */ \
} \
rsi += -32*4*clmblkdst8;	/* Restore source pointer */ \
end_timer(18); \
\
/* Do inverse FFT levels 6 */ \
/* On input the 64-byte cache lines hold these data values: */ \
/*	0K	64K	128K	192K */ \
/*	2K	... */ \
/*	... */ \
/*	62K */ \
/*	256K	... */ \
/*	... */ \
/* On output the 64-byte cache lines hold these data values: */ \
/*	0K	128K	256K	384K */ \
/*	2K	... */ \
/*	... */ \
/*	126K */ \
/*	512K	... */ \
/*	... */ \
/* */ \
/* Do 64 eight_reals_unfft_1 macros */ \
/*	distance between fft data elements is 128K */ \
/*	do 32 macros each processing 16 data values */ \
/* */ \
/* Do 960 four_complex_unfft macros */ \
/*	distance between fft data elements is 128K */ \
/*	do 480 macros each processing 16 data values */ \
\
/* Do inverse FFT levels 4,5 */ \
/* On input the 64-byte cache lines hold these data values: */ \
/*	0K	128K	256K	384K */ \
/*	2K	... */ \
/*	... */ \
/*	126K */ \
/*	512K	... */ \
/*	... */ \
/* On output the 64-byte cache lines hold these data values: */ \
/*	0K	512K	1024K	1536K */ \
/*	2K	... */ \
/*	... */ \
/*	510K */ \
/*	2048K	... */ \
/*	... */ \
/* */ \
/* Do 128 eight_reals_unfft_2 macros */ \
/* Do 128 nop_two_two_complex_unfft_2 macros */ \
/*	distance between fft data elements is 256K */ \
/*	do 64 macros each processing 32 data values */ \
/* */ \
/* Do 768 four_complex_unfft macros */ \
/*	distance between fft data elements is 256K */ \
/*	do 192 macros each processing 32 data values */ \
\
/* Do inverse FFT levels 1,2,3 */ \
/* On input the 64-byte cache lines hold these data values: */ \
/*	0K	512K	1024K	1536K */ \
/*	2K	... */ \
/*	... */ \
/*	510K */ \
/*	2048K	... */ \
/*	... */ \
/* On output the 64-byte cache lines hold these data values: */ \
/*	0K	2048K	4096K	6144K */ \
/*	2K	... */ \
/*	... */ \
/*	2046K */ \
/* */ \
/* Do 1024 eight_reals_last_unfft macros */ \
/*	distance between fft data elements is 1024K */ \
/*	do 256 macros each processing 32 data values */ \
\
start_timer(19); \
rdi = (uintptr_t)g->sincos4;/* Load sin/cos pointer */ \
rbx = (uintptr_t)g->sincos5;/* Load sin/cos pointer */ \
/* 4 iterations */ \
\
for(unsigned int loopA = 4; loopA; loopA--) { \
	for(unsigned int loopB = 2*4; loopB; loopB--) { \
		xtouchlines(9);			/* Preload TLBs */ \
		\
		for(unsigned int loopC = clm; loopC; loopC--) { \
			for(unsigned int loopD = 2; loopD; loopD--) { \
				x2cl_eight_reals_unfft_1(rsi, 64, 4*clmblkdst8); \
			} \
			\
			rsi += -2*64+8*clmblkdst8;	/* Next source pointer */ \
			/* 15 iterations */ \
			for(unsigned int loopD = 15; loopD; loopD--) { \
				for(unsigned int loopE = 2; loopE; loopE--) { \
					x2cl_two_complex_unfft(rsi, 64, 4*clmblkdst8); \
				} \
				rsi += -2*64+8*clmblkdst8;	/* Next source pointer */ \
				rdi += XMM_SCD;		/* Next sine/cosine pointer */ \
			} \
			rdi += -15*XMM_SCD;	/* Restore sine/cosine pointer */ \
			\
			rsi += -16*8*clmblkdst8;	/* Restore source pointer */ \
			for(unsigned int loopD = 2; loopD; loopD--) { \
				xloop1sp(c2b, 4*clmblkdst8, x4cl_eight_reals_unfft_2(rsi, 64, 8*clmblkdst8, 16*clmblkdst8)); \
			} \
			\
			rsi += -2*4*clmblkdst8+32*clmblkdst8;/* Next source pointer */ \
			/* 3 iterations */ \
			for(unsigned int loopD = 3; loopD; loopD--) { \
				xloop2sp(c3b, 4*clmblkdst8, x4cl_four_complex_unfft(rsi, 64, 8*clmblkdst8, 16*clmblkdst8, rbx)); \
				rsi += -2*4*clmblkdst8+32*clmblkdst8; \
				rbx += XMM_SCD;		/* Next sine/cosine pointer */ \
			} \
			rbx += -3*XMM_SCD;		/* Next sine/cosine pointer */ \
			\
			rsi += -4*32*clmblkdst8;	/* Restore source pointer */ \
			for(unsigned int loopD = 4; loopD; loopD--) { \
				xloop2sp(c1b, 4*clmblkdst8, x4cl_eight_reals_last_unfft(rsi, 64, 32*clmblkdst8, 64*clmblkdst8)); \
			} \
			\
				EXEC248##clm(rsi += -8*4*clmblkdst8+2*64);/* Next source pointer */ \
			EXEC248##clm(pfing(rcx += -9*blkdst+128));	/* Next prefetch pointer */ \
		} \
		EXEC248##clm(rsi += -clm*2*64+clmblkdst);/* Next source pointer */ \
		EXEC248##clm(pfing(rcx += -clm*128+9*blkdst));	/* Next prefetch pointer */ \
		EXEC1##clm(rsi += -8*4*clmblkdst8+clmblkdst);/* Next source pointer */ \
		\
	} \
	rsi += -8*clmblkdst+clmblkdst8;/* Next source pointer */ \
} \
end_timer(19); \
\
/* Normalize these values */ \
\
pass1_normalize(1, yes_fft); \
\
/* Copy data from the scratch area back to the FFT data area */ \
\
copy_scratch_data_to_fft(1024, clm); \
goto	endlp;			/* Skip FFT code if flag not set */ \
\
/* Fire up auxiliary threads to do pass 1 forward FFT */ \
\
/* pass1: */	pass1_forward_fft_setup(b0b); \
\
/* Copy data from the FFT data area back to the scratch area */ \
\
b0b:	copy_fft_data_to_scratch(1024, clm); \
\
/* Do FFT levels 1,2,3 */ \
/* */ \
/* On input the 64-byte cache lines hold these data values: */ \
/*	0K	2048K	4096K	6144K */ \
/*	2K	... */ \
/*	... */ \
/*	2046K */ \
/* On output the 64-byte cache lines hold these data values: */ \
/*	0K	1024K	2048K	3072K */ \
/*	2K	... */ \
/*	... */ \
/*	1022K	... */ \
/*	4096K	... */ \
/*	... */ \
/* */ \
/* Do 1024 eight_reals_first_fft macros (FFT levels 1,2,3) */ \
/*	distance between fft data elements is 1024K */ \
/*	do 256 iterations each processing 16 data values */ \
\
/* Do FFT levels 4,5 */ \
/* Values 0K-1023K is real data, 1024K-2047K is semi-real data, 2048K-8191K is */ \
/* complex data. */ \
/* */ \
/* On input the 64-byte cache lines hold these data values: */ \
/*	0K	1024K	2048K	3072K */ \
/*	2K	... */ \
/*	... */ \
/*	1022K	... */ \
/*	4096K	... */ \
/*	... */ \
/* On output the 64-byte cache lines hold these data values: */ \
/*	0K	256K	512K	768K */ \
/*	2K	... */ \
/*	... */ \
/*	254K	... */ \
/*	1024K	... */ \
/*	... */ \
/* */ \
/* Do 128 eight_reals_fft_2 macros */ \
/* Do 128 nop_two_two_complex_fft_2 macros */ \
/* Do 256 four_complex_fft macros */ \
/*	distance between fft data elements is 256K */ \
/*	do 128 macros each processing 32 data values */ \
/* */ \
/* Do 512 four_complex_fft macros */ \
/*	distance between fft data elements is 256K */ \
/*	do 128 macros each processing 32 data values */ \
\
/* Do FFT levels 6 */ \
/* Values 0K-255K is real data, 256K-511K is semi-real data, 512K-8191K is */ \
/* complex data. */ \
/* */ \
/* On input the 64-byte cache lines hold these data values: */ \
/*	0K	256K	512K	768K */ \
/*	2K	... */ \
/*	... */ \
/*	254K	... */ \
/*	1024K	... */ \
/*	... */ \
/* On output the 64-byte cache lines hold these data values: */ \
/*	0K	128K	256K	384K */ \
/*	2K	... */ \
/*	... */ \
/*	126K	... */ \
/*	512K	... */ \
/*	... */ \
/* */ \
/* Do 64 eight_reals_fft_1 macros */ \
/* Do 64 two_two_complex_fft macros */ \
/*	distance between fft data elements is 128K */ \
/*	do 64 macros each processing 16 data values */ \
/* */ \
/* Do 896 two_two_complex_fft macros */ \
/*	distance between fft data elements is 128K */ \
/*	do 448 macros each processing 16 data values */ \
\
yes_fft:start_timer(22); \
rsi = (uintptr_t)g->scratch_area;	/* Get address of scratch area */ \
rcx = (uintptr_t)g->data_addr; /* Load source address for prefetching */ \
pfing(rbp = (uintptr_t)g->data_prefetch);/* Set prefetch distance */ \
pfing(rbp -= rcx); \
pfing(rcx += 522*blkdst);		/* Create prefetching pointer */ \
rdi = (uintptr_t)g->sincos5;/* Load sin/cos pointer */ \
rbx = (uintptr_t)g->sincos4;/* Load sin/cos pointer */ \
/* 8 iterations */ \
\
for(unsigned int loopA = 8; loopA; loopA--) { \
	for(unsigned int loopB = 2*4; loopB; loopB--) { \
		xtouchlines(4);			/* Preload TLBs */ \
		\
		for(unsigned int loopC = clm; loopC; loopC--) { \
			for(unsigned int loopD = 8; loopD; loopD--) { \
				for(unsigned int loopE = 2; loopE; loopE--) { \
					x2cl_eight_reals_fft(rsi, 64, 64*clmblkdst8); \
				} \
				rsi += -2*64+8*clmblkdst8;	/* Next source pointer */ \
			} \
			\
			rsi += -8*8*clmblkdst8;	/* Restore source pointer */ \
			for(unsigned int loopD = 2; loopD; loopD--) { \
				xloop1sp(b2b, 8*clmblkdst8, x4cl_eight_reals_fft_2(rsi, 64, 16*clmblkdst8, 32*clmblkdst8, rdi)); \
			} \
			\
			rsi += -2*8*clmblkdst8+64*clmblkdst8;/* Next source pointer */ \
			rdi += XMM_SCD;		/* Next sine/cosine pointer */ \
			for(unsigned int loopD = 2; loopD; loopD--) { \
				xloop1sp(b3b, 8*clmblkdst8, x4cl_four_complex_fft(rsi, 64, 16*clmblkdst8, 32*clmblkdst8, rdi)); \
			} \
			rdi += -XMM_SCD;		/* Restore sine/cosine pointer */ \
			\
			rsi += -2*8*clmblkdst8-64*clmblkdst8;/* Restore src pointer */ \
			for(unsigned int loopD = 2; loopD; loopD--) { \
				x2cl_eight_reals_fft_1(rsi, 64, 8*clmblkdst8); \
			} \
			\
			rsi += -2*64+16*clmblkdst8;/* Load source pointer */ \
			rbx += XMM_SCD;		/* Next sine/cosine pointer */ \
			/* 7 iterations */ \
			for(unsigned int loopD = 7; loopD; loopD--) { \
				for(unsigned int loopE = 2; loopE; loopE--) { \
					x2cl_two_complex_fft(rsi, 64, 8*clmblkdst8, rbx); \
				} \
				rsi += -2*64+16*clmblkdst8;/* Next source pointer */ \
				rbx += 2*XMM_SCD;		/* Next sine/cosine pointer */ \
			} \
			rbx += -15*XMM_SCD;	/* Restore sine/cosine pointer */ \
			\
				EXEC248##clm(rsi += -8*16*clmblkdst8+2*64);/* Next source pointer */ \
			EXEC248##clm(pfing(rcx += -4*blkdst+128));	/* Next prefetch pointer */ \
		} \
		EXEC248##clm(rsi += -clm*2*64+clmblkdst);/* Next source pointer */ \
		EXEC248##clm(pfing(rcx += -clm*128+4*blkdst));	/* Next prefetch pointer */ \
		EXEC1##clm(rsi += -8*16*clmblkdst8+clmblkdst);/* Next source pointer */ \
		\
	} \
	rsi += -8*clmblkdst+clmblkdst8;/* Next source pointer */ \
} \
rsi += -8*clmblkdst8;	/* Restore source pointer */ \
end_timer(22); \
\
/* Do FFT levels 7,8 */ \
/* Values 0K-127K is real data, 128K-255K is semi-real data, 256K-8191K is */ \
/* complex data. */ \
/* */ \
/* On input the 64-byte cache lines hold these data values: */ \
/*	0K	128K	256K	384K */ \
/*	2K	... */ \
/*	... */ \
/*	126K	... */ \
/*	512K	... */ \
/*	... */ \
/* On output the 64-byte cache lines hold these data values: */ \
/*	0K	32K	64K	96K */ \
/*	2K	... */ \
/*	... */ \
/*	30K	... */ \
/*	128K	... */ \
/*	... */ \
\
/* Do 32 eight_reals_fft_2 macros */ \
/* Do 32 two_two_complex_fft macros */ \
/*	distance between fft data elements is 32K */ \
/*	do 16 macros each processing 32 data values */ \
\
start_timer(24); \
rdi = (uintptr_t)g->sincos3;/* Load sin/cos pointer */ \
/* 2 iterations of 8 iters of 2*clm */ \
for(unsigned int loopA = 2; loopA; loopA--) { \
	for(unsigned int loopB = 8; loopB; loopB--) { \
		for(unsigned int loopC = 2*clm; loopC; loopC--) { \
			x4cl_eight_reals_fft_2(rsi, 64, 2*clmblkdst8, 4*clmblkdst8, rdi); \
		} \
		rsi += -2*clm*64+clmblkdst;/* Next source pointer */ \
	} \
	rsi += -8*clmblkdst+clmblkdst8;/* Next source pointer */ \
} \
\
/* Do 960 four_complex_fft macros */ \
/*	distance between fft data elements is 32K */ \
/*	do 240 macros each processing 32 data values */ \
\
rsi += -2*clmblkdst8+8*clmblkdst8;/* Load source pointer */ \
rdi += XMM_SCD;		/* Next sine/cosine pointer */ \
for(unsigned int loopA = 15; loopA; loopA--) { \
	for(unsigned int loopB = 2; loopB; loopB--) { \
		for(unsigned int loopC = 4; loopC; loopC--) { \
			xtouch(xptr(rcx+(clm-1)*128));	/* Preload TLB */ \
			xloop2_clm##clm(b5d, clmblkdst, x4cl_four_complex_fft(rsi, 64, 2*clmblkdst8, 4*clmblkdst8, rdi)); \
		} \
		rsi += -8*clmblkdst+clmblkdst8;/* Next source pointer */ \
	} \
	rsi += -2*clmblkdst8+8*clmblkdst8;/* Next source pointer */ \
	rdi += 2*XMM_SCD;		/* Next sine/cosine pointer */ \
} \
rsi += -16*8*clmblkdst8;	/* Restore source pointer */ \
end_timer(24); \
\
/* Do FFT levels 9,10 */ \
/* Values 0K-31K is real data, 32K-63K is semi-real data, 64K-8191K is */ \
/* complex data. */ \
/* */ \
/* On input the 64-byte cache lines hold these data values: */ \
/*	0K	32K	64K	96K */ \
/*	2K	... */ \
/*	... */ \
/*	30K	... */ \
/*	128K	... */ \
/*	... */ \
/* On output the 64-byte cache lines hold these data values: */ \
/*	0K	8K	16K	24K */ \
/*	2K	... */ \
/*	... */ \
/*	6K	... */ \
/*	32K	... */ \
/*	... */ \
\
/* Do 8 eight_reals_fft_2 macros */ \
/* Do 8 two_two_complex_fft macros */ \
/*	distance between fft data elements is 8K */ \
/*	do 4 macros each processing 32 data values */ \
\
start_timer(25); \
rdi = (uintptr_t)g->sincos2;/* Load sin/cos pointer */ \
/* 4 iterations of 2*clm */ \
for(unsigned int loopA = 4; loopA; loopA--) { \
	for(unsigned int loopB = 2*clm; loopB; loopB--) { \
		x4cl_eight_reals_fft_2(rsi, 64, 4*clmblkdst, clmblkdst8, rdi); \
	} \
	rsi += -2*clm*64+clmblkdst;/* Next source pointer */ \
} \
\
/* Do 1008 four_complex_fft macros */ \
/*	distance between fft data elements is 8K */ \
/*	do 252 macros each processing 32 data values */ \
\
rsi += -4*clmblkdst+2*clmblkdst8;/* Load source pointer */ \
rdi += XMM_SCD;		/* Next sine/cosine pointer */ \
/* 63 iterations of 2 iters of 2*2*clm */ \
for(unsigned int loopA = 63; loopA; loopA--) { \
	for(unsigned int loopB = 2; loopB; loopB--) { \
		xtouch(xptr(rcx+(clm-1)*128));	/* Preload TLB */ \
		xloop2_clm##clm(b7b, clmblkdst, x4cl_four_complex_fft(rsi, 64, 4*clmblkdst, clmblkdst8, rdi)); \
	} \
	rsi += -2*2*clmblkdst+2*clmblkdst8;/* Next source pointer */ \
	rdi += 2*XMM_SCD;		/* Next sine/cosine pointer */ \
} \
rsi += -64*2*clmblkdst8;	/* Restore source pointer */ \
pfing(rcx += -1024*blkdst);	/* Restore prefetch pointer */ \
end_timer(25); \
\
/* Do FFT levels 9,10 */ \
/* Values 0K-7K is real data, 8K-15K is semi-real data, 16K-8191K is */ \
/* complex data. */ \
/* */ \
/* On input the 64-byte cache lines hold these data values: */ \
/*	0K	8K	16K	24K */ \
/*	2K	... */ \
/*	... */ \
/*	6K	... */ \
/*	32K	... */ \
/*	... */ \
/* On output the 64-byte cache lines hold these data values: */ \
/*	0K	2K	4K	6K */ \
/*	8K	... */ \
/*	... */ \
\
/* Do 1 eight_reals_fft_2 macros */ \
/* Do 1 nop_two_two_complex_fft_2 macros */ \
/* Do 2 four_complex_fft macros */ \
/*	distance between fft data elements is 2K */ \
/*	do 1 macros each processing 32 data values */ \
\
start_timer(26); \
rdi = (uintptr_t)g->sincos1;/* Load sin/cos pointer */ \
for(unsigned int loopA = 2*clm; loopA; loopA--) { \
	g4cl_eight_reals_fft_2(rsi, 64, clmblkdst, 2*clmblkdst, rcx, 64, blkdst, 2*blkdst); \
} \
\
/* Do 1020 four_complex_fft macros */ \
/*	distance between fft data elements is 2K */ \
/*	do 255 macros each processing 32 data values */ \
\
rsi += -2*clm*64+4*clmblkdst;/* Next source pointer */ \
rcx += -2*clm*64+4*blkdst;	/* Next dest pointer */ \
rdi += XMM_SCD;		/* Next sine/cosine pointer */ \
/* eax = 128*65536+1; */ /* 255 iterations of 2*clm */ \
for(unsigned int loopA = 255; loopA; ) { \
	do{ \
		for(unsigned int loopC = clm; loopC; loopC--) { \
			for(unsigned int loopD = 2; loopD; loopD--) { \
				g4cl_four_complex_fft(rsi, 64, clmblkdst, 2*clmblkdst, rcx, 64, blkdst, 2*blkdst); \
				L2prefetch128(u8ptr(rcx-64+rbp)); \
				L2prefetch128(u8ptr(rcx-64+blkdst+rbp)); \
				pfing(rbp += -64+2*blkdst); \
			} \
			pfing(rbp += 128-4*blkdst); \
		} \
		rsi += -2*clm*64+4*clmblkdst;/* Next source pointer */ \
		rcx += -2*clm*64+4*blkdst;	/* Next dest pointer */ \
		rdi += 2*XMM_SCD;		/* Next sine/cosine pointer */ \
	}while(--loopA & 1); \
	rsi += -8*clmblkdst+clmblkdst8; \
} \
end_timer(26); \
\
/* Work on next 4 sets of 1024 values */ \
\
endlp:	xpass1_get_next_block(pass2, c0b, b0b); /* Get next block and jump */ \
/* done: */	xfft_footer_hg2; \
\
\
\
/* Perform a 8M-element negacyclic FFT.  This is done in two passes.  Pass 1 */ \
/* does 12 levels, pass 2 uses common code to perform the last 11 levels. */ \
\

#define hg_pass1sclevels12complex4(procname, xpass2name, pass2_levels, clm) \
 /* (fftname, pass2_levels, clmarg) */ \
\
/* set_FFT_constants(fftname, 0, 1, pass2_levels, clmarg); */ \
xfft_header_hg2(procname, xpass2name, pass2_levels, clm); \
\
/* Do the last 11 FFT levels */ \
/*	do 256 pass2_11_levels_complex */ \
\
/* pass2: */	start_timer(0); \
xpass2name(g); \
end_timer(0); \
if(g->ffttype == 1)		/* We're done if FFTing only */ \
	return; \
\
/* Do the inverse FFT */ \
\
pass1_inverse_fft_setup(c0b, procname); \
\
/* Do inverse FFT levels 11,12 */ \
/* On input the 64-byte cache lines hold these data values: */ \
/*	0K	2K	4K	6K */ \
/*	8K	10K	12K	14K */ \
/*	... */ \
/* On output the 64-byte cache lines hold these data values: */ \
/*	0K	4K	8K	12K */ \
/*	2K	... */ \
/*	16K	... */ \
/*	... */ \
\
/* Do 1024 four_complex_unfft macros */ \
/*	distance between fft data elements is 2K */ \
/*	do 512 macros each processing 16 data values */ \
\
c0b:	start_timer(15); \
rsi = (uintptr_t)g->data_addr;/* Get source address */ \
rcx = (uintptr_t)g->scratch_area;	/* Get address of scratch area */ \
rdi = (uintptr_t)g->sincos1;/* Load sin/cos pointer */ \
/* 128 iterations */ \
pfing(rbp = (uintptr_t)g->data_prefetch);/* Set prefetch distance */ \
pfing(rbp -= rsi); \
for(unsigned int loopA = 128; loopA; loopA--) { \
	for(unsigned int loopB = 4; loopB; loopB--) { \
		for(unsigned int loopC = 2*clm; loopC; loopC--) { \
			g2cl_four_complex_unfft(rsi, 64, blkdst, rcx, 64, clmblkdst); \
		} \
		rsi += -2*clm*64+2*blkdst;	/* Next source pointer */ \
		rcx += -2*clm*64+2*clmblkdst;/* Next dest pointer */ \
		rdi += XMM_SCD;		/* Next sine/cosine pointer */ \
	} \
	rcx += -8*clmblkdst+clmblkdst8; \
} \
rcx = rsi-512*2*blkdst;	/* Create prefetch pointer */ \
rsi = (uintptr_t)g->scratch_area;	/* Restore source pointer */ \
end_timer(15); \
\
/* Do inverse FFT levels 9,10 */ \
/* On input the 64-byte cache lines hold these data values: */ \
/*	0K	4K	8K	12K */ \
/*	2K	... */ \
/*	16K	... */ \
/*	... */ \
/* On output the 64-byte cache lines hold these data values: */ \
/*	0K	16K	32K	48K */ \
/*	2K	... */ \
/*	... */ \
/*	14K */ \
/*	64K	... */ \
/*	... */ \
\
/* Do 1024 four_complex_unfft macros */ \
/*	distance between fft data elements is 8K */ \
/*	do 256 macros each processing 32 data values */ \
\
start_timer(16); \
rdi = (uintptr_t)g->sincos2;/* Load sin/cos pointer */ \
/* 128 iterations of 2*2*clm */ \
for(unsigned int loopA = 128; loopA; loopA--) { \
	xtouch(xptr(rcx+(clm-1)*128));	/* Preload TLB */ \
	xloop2_clm##clm(c7b, clmblkdst, x4cl_four_complex_unfft(rsi, 64, 2*clmblkdst, 4*clmblkdst, rdi));		/* 2 * 2*clm iterations */ \
	rsi += -2*clmblkdst+clmblkdst8;/* Next source pointer */ \
	rdi += XMM_SCD;		/* Next sine/cosine pointer */ \
} \
rsi += -128*clmblkdst8;	/* Restore source pointer */ \
end_timer(16); \
\
/* Do inverse FFT levels 7,8 */ \
/* On input the 64-byte cache lines hold these data values: */ \
/*	0K	16K	32K	48K */ \
/*	2K	... */ \
/*	... */ \
/*	14K */ \
/*	64K	... */ \
/*	... */ \
/* On output the 64-byte cache lines hold these data values: */ \
/*	0K	64K	128K	192K */ \
/*	2K	... */ \
/*	... */ \
/*	62K */ \
/*	256K	... */ \
/*	... */ \
\
/* Do 1024 four_complex_unfft macros */ \
/*	distance between fft data elements is 32K */ \
/*	do 256 macros each processing 32 data values */ \
\
start_timer(17); \
rdi = (uintptr_t)g->sincos3;/* Load sin/cos pointer */ \
/* eax = 32*256+4; */ /* 24 iterations of 4 iters of 2*2*clm   WARNING comment */ \
for(unsigned int loopA = 32; loopA; loopA--) { \
	for(unsigned int loopB = 4; loopB; loopB--) { \
		xtouch(xptr(rcx+(clm-1)*128));	/* Preload TLB */ \
		xloop2_clm##clm(c6b, clmblkdst, x4cl_four_complex_unfft(rsi, 64, clmblkdst8, 2*clmblkdst8, rdi));		/* 2*2*clm iterations */ \
	} \
	rsi += -4*2*clmblkdst+4*clmblkdst8;/* Next source pointer */ \
	rdi += XMM_SCD;		/* Next sine/cosine pointer */ \
} \
rsi += -32*4*clmblkdst8;	/* Restore source pointer */ \
end_timer(17); \
\
/* Do inverse FFT levels 5,6 */ \
/* On input the 64-byte cache lines hold these data values: */ \
/*	0K	64K	128K	192K */ \
/*	2K	... */ \
/*	... */ \
/*	62K */ \
/*	256K	... */ \
/*	... */ \
/* On output the 64-byte cache lines hold these data values: */ \
/*	0K	256K	512K	768K */ \
/*	2K	... */ \
/*	... */ \
/*	254K */ \
/*	1024K	... */ \
/*	... */ \
\
/* Do 1024 four_complex_unfft macros */ \
/*	distance between fft data elements is 128K */ \
/*	do 256 macros each processing 32 data values */ \
\
start_timer(18); \
rdi = (uintptr_t)g->sincos4;/* Load sin/cos pointer */ \
/* 8 iterations of 4 iters of 2*2*clm */ \
for(unsigned int loopA = 8; loopA; loopA--) { \
	for(unsigned int loopB = 4; loopB; loopB--) { \
		for(unsigned int loopC = 2; loopC; loopC--) { \
			xtouch(xptr(rcx+(clm-1)*128));	/* Preload TLB */ \
			xloop4_clm##clm(c5b, clmblkdst, x4cl_four_complex_unfft(rsi, 64, 4*clmblkdst8, 8*clmblkdst8, rdi));		/* 4 * 2*clm iterations */ \
		} \
		rsi += -2*4*clmblkdst+clmblkdst8;/* Next source pointer */ \
	} \
	rsi += -4*clmblkdst8+16*clmblkdst8;/* Next source pointer */ \
	rdi += XMM_SCD;		/* Next sine/cosine pointer */ \
} \
rsi += -8*16*clmblkdst8;	/* Restore source pointer */ \
end_timer(18); \
\
/* Do inverse FFT levels 4 */ \
/* On input the 64-byte cache lines hold these data values: */ \
/*	0K	256K	512K	768K */ \
/*	2K	... */ \
/*	... */ \
/*	254K */ \
/*	1024K	... */ \
/*	... */ \
/* On output the 64-byte cache lines hold these data values: */ \
/*	0K	512K	1024K	1536K */ \
/*	2K	... */ \
/*	... */ \
/*	510K */ \
/*	2048K	... */ \
/*	... */ \
\
/* Do 1024 two_two_complex_unfft macros */ \
/*	distance between fft data elements is 256K */ \
/*	do 512 macros each processing 16 data values */ \
\
start_timer(19); \
rdi = (uintptr_t)g->sincos5;/* Load sin/cos pointer */ \
/* 4 iters of 16 iters 2 of 4*2*clm */ \
for(unsigned int loopA = 4; loopA; loopA--) { \
	for(unsigned int loopB = 16; loopB; loopB--) { \
		xtouch(xptr(rcx+(clm-1)*128));	/* Preload TLB */ \
		xloop8_clm##clm(c3b, clmblkdst, x2cl_two_complex_unfft(rsi, 64, 16*clmblkdst8));		/* 8 * 2*clm iterations */ \
		rsi += -8*clmblkdst+clmblkdst8;/* Next source pointer */ \
	} \
	rsi += -16*clmblkdst8+32*clmblkdst8;/* Next source pointer */ \
	rdi += XMM_SCD;		/* Next sine/cosine pointer */ \
} \
rsi += -4*32*clmblkdst8;	/* Restore source pointer */ \
end_timer(19); \
\
/* Do inverse FFT levels 1,2,3 */ \
/* On input the 64-byte cache lines hold these data values: */ \
/*	0K	512K	1024K	1536K */ \
/*	2K	... */ \
/*	... */ \
/*	510K */ \
/*	2048K	... */ \
/*	... */ \
/* On output the 64-byte cache lines hold these data values: */ \
/*	0K	2048K	4096K	6144K */ \
/*	2K	... */ \
/*	... */ \
/*	2046K */ \
\
/* Do 1024 four_complex_unfft macros */ \
/*	distance between fft data elements is 1024K */ \
/*	do 256 macros each processing 32 data values */ \
\
start_timer(20); \
rdi = (uintptr_t)g->plus1_premults;/* Address of premultiplier table */ \
/* eax = 32*65536+4; */ /* 8 iters of 4 iters of 2*2*clm   WARNING comment */ \
for(unsigned int loopA = 32; loopA; loopA--) { \
	for(unsigned int loopB = 4; loopB; loopB--) { \
		xtouch(xptr(rcx+(clm-1)*128));	/* Preload TLB */ \
		xloop2edi_clm##clm(c1b, clmblkdst, 128, x4cl_four_complex_last_unfft(rsi, 64, 32*clmblkdst8, 64*clmblkdst8, 256*128));	/* 2 * 2*clm iterations */ \
	} \
	rsi += -4*2*clmblkdst+clmblkdst8;/* Next source pointer */ \
} \
end_timer(20); \
\
/* Normalize these values */ \
\
pass1_normalize(1, yes_fft); \
\
/* Copy data from the scratch area back to the FFT data area */ \
\
copy_scratch_data_to_fft(1024, clm); \
goto	endlp;			/* Skip FFT code if flag not set */ \
\
/* Fire up auxiliary threads to do pass 1 forward FFT */ \
\
/* pass1: */	pass1_forward_fft_setup(b0b); \
\
/* Copy data from the FFT data area back to the scratch area */ \
\
b0b:	copy_fft_data_to_scratch(1024, clm); \
\
/* Do FFT levels 1,2,3 */ \
/* */ \
/* On input the 64-byte cache lines hold these data values: */ \
/*	0K	2048K	4096K	6144K */ \
/*	2K	... */ \
/*	... */ \
/*	2046K */ \
/* On output the 64-byte cache lines hold these data values: */ \
/*	0K	1024K	2048K	3072K */ \
/*	2K	... */ \
/*	... */ \
/*	1022K	... */ \
/*	4096K	... */ \
/*	... */ \
\
/* Do 1024 four_complex_first_fft macros (FFT levels 1,2,3) */ \
/*	distance between fft data elements is 1024K */ \
/*	do 512 iterations each processing 16 data values */ \
\
yes_fft:start_timer(22); \
rsi = (uintptr_t)g->scratch_area;	/* Get address of scratch area */ \
rdi = (uintptr_t)g->plus1_premults;/* Address of premultiplier table */ \
rcx = (uintptr_t)g->data_addr; /* Load source address for prefetching */ \
pfing(rbp = (uintptr_t)g->data_prefetch);/* Set prefetch distance */ \
pfing(rbp -= rcx); \
pfing(rcx += 512*blkdst);		/* Create prefetching pointer */ \
/* 64 iters of 2 iters of 4*2*clm */ \
for(unsigned int loopA = 64; loopA; loopA--) { \
	for(unsigned int loopB = 2; loopB; loopB--) { \
		xtouch(xptr(rcx+(clm-1)*128));	/* Preload TLB */ \
		xloop4edi_clm##clm(b1b, clmblkdst, 128, x2cl_four_complex_first_fft_scratch(rsi, 64, 64*clmblkdst8));	/* 4 * 2*clm iterations */ \
	} \
	rsi += -2*4*clmblkdst+clmblkdst8; \
} \
rsi += -64*clmblkdst8;	/* Restore source pointer */ \
end_timer(22); \
\
/* Do FFT levels 4 */ \
/* */ \
/* On input the 64-byte cache lines hold these data values: */ \
/*	0K	1024K	2048K	3072K */ \
/*	2K	... */ \
/*	... */ \
/*	1022K	... */ \
/*	4096K	... */ \
/*	... */ \
/* On output the 64-byte cache lines hold these data values: */ \
/*	0K	512K	1024K	1536K */ \
/*	2K	... */ \
/*	... */ \
/*	510K	... */ \
/*	2048K	... */ \
/*	... */ \
\
/* Do 1024 two_two_complex_fft macros */ \
/*	distance between fft data elements is 256K */ \
/*	do 512 macros each processing 16 data values */ \
\
start_timer(23); \
rdi = (uintptr_t)g->sincos5;/* Load sin/cos pointer */ \
/* 2 iters of 32 iters of 2 of 4*2*clm */ \
for(unsigned int loopA = 2; loopA; loopA--) { \
	for(unsigned int loopB = 32; loopB; loopB--) { \
		xtouch(xptr(rcx+(clm-1)*128));	/* Preload TLB */ \
		xloop8_clm##clm(b3b, clmblkdst, x2cl_two_complex_fft(rsi, 64, 32*clmblkdst8, rdi));		/* 8 * 2*clm iterations */ \
		rsi += -8*clmblkdst+clmblkdst8; \
	} \
	rsi += -32*clmblkdst8+64*clmblkdst8;/* Next source pointer */ \
	rdi += 2*XMM_SCD;		/* Next sine/cosine pointer */ \
} \
rsi += -2*64*clmblkdst8;	/* Restore source pointer */ \
end_timer(23); \
\
/* Do FFT levels 5,6 */ \
/* */ \
/* On input the 64-byte cache lines hold these data values: */ \
/*	0K	512K	1024K	1536K */ \
/*	2K	... */ \
/*	... */ \
/*	510K	... */ \
/*	2048K	... */ \
/*	... */ \
/* On output the 64-byte cache lines hold these data values: */ \
/*	0K	128K	256K	384K */ \
/*	2K	... */ \
/*	... */ \
/*	126K	... */ \
/*	512K	... */ \
/*	... */ \
\
/* Do 1024 four_complex_fft macros */ \
/*	distance between fft data elements is 128K */ \
/*	do 256 macros each processing 32 data values */ \
\
start_timer(24); \
rdi = (uintptr_t)g->sincos4;/* Load sin/cos pointer */ \
/* 4 iters of 8 iters of 4 of 2*2*clm */ \
for(unsigned int loopA = 4; loopA; loopA--) { \
	for(unsigned int loopB = 8; loopB; loopB--) { \
		for(unsigned int loopC = 4; loopC; loopC--) { \
			xtouch(xptr(rcx+(clm-1)*128));	/* Preload TLB */ \
			xloop2_clm##clm(b4b, clmblkdst, x4cl_four_complex_fft(rsi, 64, 8*clmblkdst8, 16*clmblkdst8, rdi));		/* 2 * 2*clm iterations */ \
		} \
		rsi += -4*2*clmblkdst+clmblkdst8; \
	} \
	rsi += -8*clmblkdst8+32*clmblkdst8;/* Next source pointer */ \
	rdi += 2*XMM_SCD;		/* Next sine/cosine pointer */ \
} \
rsi += -4*32*clmblkdst8;	/* Restore source pointer */ \
end_timer(24); \
\
/* Do FFT levels 7,8 */ \
/* */ \
/* On input the 64-byte cache lines hold these data values: */ \
/*	0K	128K	256K	384K */ \
/*	2K	... */ \
/*	... */ \
/*	126K	... */ \
/*	512K	... */ \
/*	... */ \
/* On output the 64-byte cache lines hold these data values: */ \
/*	0K	32K	64K	96K */ \
/*	2K	... */ \
/*	... */ \
/*	30K	... */ \
/*	128K	... */ \
/*	... */ \
\
/* Do 1024 four_complex_fft macros */ \
/*	distance between fft data elements is 32K */ \
/*	do 256 macros each processing 32 data values */ \
\
start_timer(24); \
rdi = (uintptr_t)g->sincos3;/* Load sin/cos pointer */ \
/* 4 iters of 2 iters of 4 of 2*2*clm */ \
for(unsigned int loopA = 16; loopA; loopA--) { \
	for(unsigned int loopB = 2; loopB; loopB--) { \
		for(unsigned int loopC = 4; loopC; loopC--) { \
			xtouch(xptr(rcx+(clm-1)*128));	/* Preload TLB */ \
			xloop2_clm##clm(b5b, clmblkdst, x4cl_four_complex_fft(rsi, 64, 2*clmblkdst8, 4*clmblkdst8, rdi));		/* 2 * 2*clm iterations */ \
		} \
		rsi += -4*2*clmblkdst+clmblkdst8; \
	} \
	rsi += -2*clmblkdst8+8*clmblkdst8;/* Next source pointer */ \
	rdi += 2*XMM_SCD;		/* Next sine/cosine pointer */ \
} \
rsi += -16*8*clmblkdst8;	/* Restore source pointer */ \
end_timer(24); \
\
/* Do FFT levels 9,10 */ \
/* */ \
/* On input the 64-byte cache lines hold these data values: */ \
/*	0K	32K	64K	96K */ \
/*	2K	... */ \
/*	... */ \
/*	30K	... */ \
/*	128K	... */ \
/*	... */ \
/* On output the 64-byte cache lines hold these data values: */ \
/*	0K	8K	16K	24K */ \
/*	2K	... */ \
/*	... */ \
/*	6K	... */ \
/*	32K	... */ \
/*	... */ \
\
/* Do 1024 four_complex_fft macros */ \
/*	distance between fft data elements is 8K */ \
/*	do 256 macros each processing 32 data values */ \
\
start_timer(25); \
rdi = (uintptr_t)g->sincos2;/* Load sin/cos pointer */ \
/* 64 iterations of 4*2*clm */ \
for(unsigned int loopA = 64; loopA; loopA--) { \
	xtouch(xptr(rcx+(clm-1)*128));	/* Preload TLB */ \
	xloop4_clm##clm(b7b, clmblkdst, x4cl_four_complex_fft(rsi, 64, 4*clmblkdst, clmblkdst8, rdi));		/* 2 * 2*clm iterations */ \
	rsi += -4*clmblkdst+2*clmblkdst8;/* Next source pointer */ \
	rdi += 2*XMM_SCD;		/* Next sine/cosine pointer */ \
} \
rsi += -64*2*clmblkdst8;	/* Restore source pointer */ \
pfing(rcx += -1024*blkdst);	/* Restore prefetch pointer */ \
end_timer(25); \
\
/* Do FFT levels 11,12 */ \
/* */ \
/* On input the 64-byte cache lines hold these data values: */ \
/*	0K	8K	16K	24K */ \
/*	2K	... */ \
/*	... */ \
/*	6K	... */ \
/*	32K	... */ \
/*	... */ \
/* On output the 64-byte cache lines hold these data values: */ \
/*	0K	2K	4K	6K */ \
/*	8K	... */ \
/*	... */ \
\
/* Do 1024 four_complex_fft macros */ \
/*	distance between fft data elements is 2K */ \
/*	do 256 macros each processing 32 data values */ \
\
start_timer(26); \
rdi = (uintptr_t)g->sincos1;/* Load sin/cos pointer */ \
/* 128 iters of 2 of 2*clm */ \
for(unsigned int loopA = 128; loopA; loopA--) { \
	for(unsigned int loopB = 2; loopB; loopB--) { \
		for(unsigned int loopC = clm; loopC; loopC--) { \
			for(unsigned int loopD = 2; loopD; loopD--) { \
				g4cl_four_complex_fft(rsi, 64, clmblkdst, 2*clmblkdst, rcx, 64, blkdst, 2*blkdst); \
				L2prefetch128(u8ptr(rcx-64+rbp)); \
				L2prefetch128(u8ptr(rcx-64+blkdst+rbp)); \
				pfing(rbp += -64+2*blkdst); \
			} \
			pfing(rbp += 128-4*blkdst); \
		} \
	rsi += -2*clm*64+4*clmblkdst;/* Next source pointer */ \
		rcx += -2*clm*64+4*blkdst;	/* Next dest pointer */ \
		rdi += 2*XMM_SCD;		/* Next sine/cosine pointer */ \
	} \
	rsi += -8*clmblkdst+clmblkdst8; \
} \
end_timer(26); \
\
/* Work on next 4 sets of 1024 values */ \
\
endlp:	xpass1_get_next_block(pass2, c0b, b0b); /* Get next block and jump */ \
/* done: */	xfft_footer_hg2; \
\
\
/* */ \
/* Utility macros */ \
/* */ \
\
/* Macros that help in looping.  They are called right after a building block */ \
/* macro. */ \
\
/* This macro then loops 2*clm times. If prefetching, a prefetch is done */ \
/* every 2 times the building block macro is called. */ \

#define xloop1_clm1(L486, incr, ...) \
for(unsigned int loopZ = 2; loopZ; loopZ--) { \
	__VA_ARGS__; \
} \
rsi += -2*1*64+incr;	/* Next source pointer */ \
L2prefetch128(u8ptr(rcx+rbp)); \
pfing(rcx += blkdst);		/* Next prefetch pointer */ \
\
\

#define xloop1_clmGE2(clm, incr, ...) \
for(unsigned int loopZ = clm; loopZ; loopZ--) { \
	for(unsigned int loopY = 2; loopY; loopY--) { \
		__VA_ARGS__; \
	} \
	L2prefetch128(u8ptr(rcx+rbp)); \
	rcx += 128; \
} \
rsi += -2*clm*64+incr;	/* Next source pointer */ \
rcx += -clm*128+blkdst;	/* Next prefetch pointer */ \

#define xloop1_clm2(L486, ...) xloop1_clmGE2(2, __VA_ARGS__)
#define xloop1_clm4(L486, ...) xloop1_clmGE2(4, __VA_ARGS__)
#define xloop1_clm8(L486, ...) xloop1_clmGE2(8, __VA_ARGS__)


/* This macro then loops 2*clm times. If prefetching, two prefetches are done */ \
/* every 2 times the building block macro is called. */ \

#define xloop1prefetch2_clm1(L486, incr, ...) \
for(unsigned int loopZ = 2; loopZ; loopZ--) { \
	__VA_ARGS__; \
	L2prefetch128(u8ptr(rcx+rbp)); \
} \
pfing(rcx += blkdst);		/* Next prefetch pointer */ \
rsi += -2*1*64+incr;	/* Next source pointer */ \
\
\

#define xloop1prefetch2_clmGE2(clm, incr, ...) \
for(unsigned int loopZ = clm; loopZ; loopZ--) { \
	for(unsigned int loopY = 2; loopY; loopY--) { \
		__VA_ARGS__; \
		L2prefetch128(u8ptr(rcx+rbp)); \
		rcx += blkdst;		/* Next prefetch pointer */ \
	} \
	rcx += -2*blkdst+128;	/* Next prefetch pointer */ \
} \
rsi += -2*clm*64+incr;	/* Next source pointer */ \
rcx += -clm*128+2*blkdst;	/* Next prefetch pointer */ \

#define xloop1prefetch2_clm2(L486, incr, ...) xloop1prefetch2_clmGE2(2, incr, __VA_ARGS__)
#define xloop1prefetch2_clm4(L486, incr, ...) xloop1prefetch2_clmGE2(4, incr, __VA_ARGS__)
#define xloop1prefetch2_clm8(L486, incr, ...) xloop1prefetch2_clmGE2(8, incr, __VA_ARGS__)

/* This macro then loops 4*clm times. If prefetching, a prefetch is done */ \
/* every 4 times the building block macro is called. */ \

#define xloop2_clm1(L486, incr, ...) xloop2edi_clm1(L486, incr, 0, __VA_ARGS__)
#define xloop2_clm2(L486, incr, ...) xloop2edi_clm2(L486, incr, 0, __VA_ARGS__)
#define xloop2_clm4(L486, incr, ...) xloop2edi_clmGE4(4, incr, 0, __VA_ARGS__)
#define xloop2_clm8(L486, incr, ...) xloop2edi_clmGE4(8, incr, 0, __VA_ARGS__)

#define xloop2edi_clm1(L486, incr, incr_edi, ...) \
for(unsigned int loopZ = 2; loopZ; loopZ--) { \
	for(unsigned int loopY = 2; loopY; loopY--) { \
		__VA_ARGS__; \
	} \
	rsi += -2*1*64+incr;	/* Next source pointer */ \
	rdi += incr_edi;		/* Next premultiplier pointer */ \
} \
L2prefetch128(u8ptr(rcx+rbp)); \
pfing(rcx += blkdst);		/* Next prefetch pointer */ \
\

#define xloop2edi_clm2(L486, incr, incr_edi, ...) \
for(unsigned int loopZ = 2; loopZ; loopZ--) { \
	for(unsigned int loopY = 2*2; loopY; loopY--) { \
		__VA_ARGS__; \
	} \
	L2prefetch128(u8ptr(rcx+rbp)); \
	rcx += 128; \
	rsi += -2*2*64+incr;	/* Next source pointer */ \
	rdi += incr_edi;		/* Next premultiplier pointer */ \
} \
rcx += -2*128+blkdst;	/* Next prefetch pointer */ \
\

#define xloop2edi_clmGE4(clm, incr, incr_edi, ...) \
for(unsigned int loopZ = 2; loopZ; loopZ--) { \
	for(unsigned int loopY = clm*2; loopY; loopY--) { \
		for(unsigned int loopX = 4; loopX; loopX--) { \
			__VA_ARGS__; \
		} \
		L2prefetch128(u8ptr(rcx+rbp)); \
		rcx += 128; \
	} \
	rsi += -2*clm*64+incr;	/* Next source pointer */ \
	rdi += incr_edi;		/* Next premultiplier pointer */ \
} \
rcx += -clm*128+blkdst;	/* Next prefetch pointer */ \

#define xloop2edi_clm4(L486, ...) xloop2edi_clmGE4(4, __VA_ARGS__)
#define xloop2edi_clm8(L486, ...) xloop2edi_clmGE4(8, __VA_ARGS__)

/* This macro then loops 8*clm times. If prefetching, a prefetch is done */ \
/* every 8 times the building block macro is called. */ \

#define xloop4_clm1(L486, incr, ...) xloop4edi_clm1(L486, incr, 0, __VA_ARGS__)
#define xloop4_clm2(L486, incr, ...) xloop4edi_clm2(L486, incr, 0, __VA_ARGS__)
#define xloop4_clm4(L486, incr, ...) xloop4edi_clm4(L486, incr, 0, __VA_ARGS__)
#define xloop4_clm8(L486, incr, ...) xloop4edi_clm8(L486, incr, 0, __VA_ARGS__)

#define xloop4edi_clm1(L486, incr, incr_edi, ...) \
for(unsigned int loopZ = 4; loopZ; loopZ--) { \
	for(unsigned int loopY = 2*1; loopY; loopY--) { \
		__VA_ARGS__; \
	} \
	rsi += -2*1*64+incr;	/* Next source pointer */ \
	rdi += incr_edi;		/* Next premultiplier pointer */ \
} \
L2prefetch128(u8ptr(rcx+rbp)); \
pfing(rcx += blkdst);		/* Next prefetch pointer */ \
\
\

#define xloop4edi_clm2(L486, incr, incr_edi, ...) \
for(unsigned int loopZ = 2; loopZ; loopZ--) { \
	for(unsigned int loopY = 2; loopY; loopY--) { \
		for(unsigned int loopX = 2*2; loopX; loopX--) { \
			__VA_ARGS__; \
		} \
		rsi += -2*2*64+incr;	/* Next source pointer */ \
		rdi += incr_edi;		/* Next premultiplier pointer */ \
	} \
	L2prefetch128(u8ptr(rcx+rbp)); \
	rcx += 128; \
} \
rcx += -2*128+blkdst;	/* Next prefetch pointer */ \
\
\

#define xloop4edi_clm4(L486, incr, incr_edi, ...) \
for(unsigned int loopZ = 4; loopZ; loopZ--) { \
	for(unsigned int loopY = 2*4; loopY; loopY--) { \
		__VA_ARGS__; \
	} \
	rsi += -2*4*64+incr;	/* Next source pointer */ \
	rdi += incr_edi;		/* Next premultiplier pointer */ \
	L2prefetch128(u8ptr(rcx+rbp)); \
	rcx += 128; \
} \
rcx += -4*128+blkdst;	/* Next prefetch pointer */ \
\
\

#define xloop4edi_clm8(L486, incr, incr_edi, ...) \
for(unsigned int loopZ = 4; loopZ; loopZ--) { \
	for(unsigned int loopY = (8*4); loopY; loopY--) { \
		for(unsigned int loopX = 2*4; loopX; loopX--) { \
			__VA_ARGS__; \
		} \
		L2prefetch128(u8ptr(rcx+rbp)); \
		rcx += 128; \
	} \
	rsi += -2*8*64+incr;	/* Next source pointer */ \
	rdi += incr_edi;		/* Next premultiplier pointer */ \
} \
rcx += -8*128+blkdst;	/* Next prefetch pointer */ \
\
\
/* This macro then loops 8*2*clm times. If prefetching, a prefetch is done */ \
/* every 16 times the building block macro is called. */ \

#define xloop8_clm1(L486, incr, ...) xloop8edi_clm1(L486, incr, 0, __VA_ARGS__)
#define xloop8_clm2(L486, incr, ...) xloop8edi_clm2(L486, incr, 0, __VA_ARGS__)
#define xloop8_clm4(L486, incr, ...) xloop8edi_clm4(L486, incr, 0, __VA_ARGS__)
#define xloop8_clm8(L486, incr, ...) xloop8edi_clm8(L486, incr, 0, __VA_ARGS__)

#define xloop8edi_clm1(L486, incr, incr_edi, ...) \
for(unsigned int loopZ = 8; loopZ; loopZ--) { \
	for(unsigned int loopY = 2*1; loopY; loopY--) { \
		__VA_ARGS__; \
	} \
	rsi += -2*1*64+incr;	/* Next source pointer */ \
	rdi += incr_edi;		/* Next premultiplier pointer */ \
} \
L2prefetch128(u8ptr(rcx+rbp)); \
pfing(rcx += blkdst);		/* Next prefetch pointer */ \
\
\

#define xloop8edi_clm2(L486, incr, incr_edi, ...) \
for(unsigned int loopZ = 2; loopZ; loopZ--) { \
	for(unsigned int loopY = 4; loopY; loopY--) { \
		for(unsigned int loopX = 4; loopX; loopX--) { \
			__VA_ARGS__; \
		} \
		rsi += -2*2*64+incr;	/* Next source pointer */ \
		rdi += incr_edi;		/* Next premultiplier pointer */ \
	} \
	L2prefetch128(u8ptr(rcx+rbp)); \
	rcx += 128; \
} \
rcx += -2*128+blkdst;	/* Next prefetch pointer */ \
\
\

#define xloop8edi_clm4(L486, incr, incr_edi, ...) \
for(unsigned int loopZ = 4; loopZ; loopZ--) { \
	for(unsigned int loopY = 2; loopY; loopY--) { \
		for(unsigned int loopX = 2*4; loopX; loopX--) { \
			__VA_ARGS__; \
		} \
		rsi += -2*4*64+incr;	/* Next source pointer */ \
		rdi += incr_edi;		/* Next premultiplier pointer */ \
	} \
	L2prefetch128(u8ptr(rcx+rbp)); \
	rcx += 128; \
} \
rcx += -4*128+blkdst;	/* Next prefetch pointer */ \
\
\

#define xloop8edi_clm8(L486, incr, incr_edi, ...) \
for(unsigned int loopZ = 8; loopZ; loopZ--) { \
	for(unsigned int loopY = 2*8; loopY; loopY--) { \
		__VA_ARGS__; \
	} \
	rsi += -2*8*64+incr;	/* Next source pointer */ \
	rdi += incr_edi;		/* Next premultiplier pointer */ \
	L2prefetch128(u8ptr(rcx+rbp)); \
	rcx += 128; \
} \
rcx += -8*128+blkdst;	/* Next prefetch pointer */ \
\
\
\

#define xtouchlines(num) \
\
\

#define xloop1sp(L486, incr, func1) \
for(unsigned int loopZ = 2; loopZ; loopZ--) { \
	func1; \
} \
rsi += -2*64+incr;		/* Next source pointer */ \
L2prefetch128(u8ptr(rcx+rbp)); \
pfing(rcx += blkdst);		/* Next prefetch pointer */ \

#define xloop1prefetch2sp(L486, incr, func1) \
for(unsigned int loopZ = 2; loopZ; loopZ--) { \
	func1; \
	L2prefetch128(u8ptr(rcx+rbp)); \
	pfing(rcx += blkdst);		/* Next prefetch pointer */ \
} \
rsi += -2*64+incr;		/* Next source pointer */ \

#define xloop2sp(L486, incr, func1) \
for(unsigned int loopZ = 2; loopZ; loopZ--) { \
	for(unsigned int loopY = 2; loopY; loopY--) { \
		func1; \
	} \
	rsi = rsi-2*64+incr;	/* Next source pointer */ \
} \
L2prefetch128(u8ptr(rcx+rbp, ...)); \
pfing(rcx += blkdst);		/* Next prefetch pointer */ \

#define xloop4sp(L486, incr, func1) \
for(unsigned int loopZ = 4; loopZ; loopZ--) { \
	for(unsigned int loopY = 2; loopY; loopY--) { \
		func1; \
	} \
	rsi = rsi-2*64+incr;	/* Next source pointer */ \
} \
L2prefetch128(u8ptr(rcx+rbp)); \
pfing(rcx += blkdst);		/* Next prefetch pointer */ \
\
