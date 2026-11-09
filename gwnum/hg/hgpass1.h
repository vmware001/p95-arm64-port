#pragma once

#define xfft_header_hg1(procname, xpass2name, pass2_levels, clm) \
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
		/*puts("fft-started flag is 0, pass1 forward");*/ \
		g->next_block = 0; \
		/*start_timer(1);*/ \
		g->thread_work_routine = procname##_i1; \
		pass1_wake_up_threads(g); \
		procname##_i (g, 1); \
	}else{ \
		/*puts("start pass2");*/ \
		procname##_i (g, 4); /* Start pass2 */ \
	} \
} \
void procname##_i (struct gwasm_data *__restrict g, int branch1) \
{ \
	uintptr_t rbx,rsi,rdi; \
	vec2f64 xmm0,xmm1,xmm2,xmm3,xmm4,xmm5,xmm6,xmm7; \
	rsi = (uintptr_t)g->DESTARG; \
	/*rbx = g->DIST_TO_FFTSRCARG;*/ \
	/* xmult.mac  set_FFT_constants */ \
	/*const uintptr_t clmblkdst = clm*128;*/ \
	/*const uintptr_t clmblkdst8 = clmblkdst*8+128;*/ \
	const uintptr_t p2cl = (pass2_levels >= 64) ? (pass2_levels / 4) : (1<<(pass2_levels-2)); \
	const uintptr_t blkdst_tmp1 = (p2cl*128+p2cl/64*128); \
	const uintptr_t blkdst = (blkdst_tmp1 % 256 == 0) ? (blkdst_tmp1-128) : blkdst_tmp1; \
	P1(procname); /*INITPASS_hg1(pass2_levels, clm);*/ \
	do{ \
		if(branch1 == 0) goto c0b; /* Do another inverse FFT/norm/FFT  */ \
		if(branch1 & 1) goto b0b; /* Do another forward FFT */ \

#define xfft_footer_hg1 \
	; }while(1); \
}

#ifdef DEBUG524
#define N452(A) #A
#define N_TO_STRING(A) N452(A)
#define P1(N) static char Z1=0; if(Z1 == 0) \
	{Z1=1; puts(#N " " N_TO_STRING(__FILE__)  " " N_TO_STRING(__LINE__));}
#else
#define P1(N)
#endif

#if 1
// branch1: 0=inverse 1=forward 2=done 3=exit 4=pass2 5=carries
#define xpass1_get_next_block(pass2, c0b, b0b) \
do{ \
	branch1 = g->pass1_get_next_block(g); \
	if(branch1 & 2) /* 2 or 3 */ \
		{return;} \
	if(branch1 != 5) \
		{break;} \
	start_timer(29); \
	xgw_carries(g); \
	end_timer(29); \
}while(1);
#else
#define xpass1_get_next_block(pass2, c0b, b0b) \
do{ /*again1:*/ \
	branch1 = g->pass1_get_next_block(g); \
	if(branch1 == 2 || branch1 == 3) \
		return; \
	if(branch1 > 4) { \
		start_timer(29); \
		xgw_carries(g); \
		end_timer(29); \
		/*goto again1;*/ \
	}else \
		break; \
}while(1);
#endif

#define pass1_forward_fft_setup(LABEL1)

#define pass1_inverse_fft_setup(LABEL1, procname) \
	g->next_block = 1; \
	/*start_timer(1);*/ \
	g->thread_work_routine = procname##_i0 ; \
	pass1_wake_up_threads(g); \


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
/* Perform a 40K-element FFT.  This is done in two passes.  Pass 1 */ \
/* does 5 levels, pass 2 uses common code to perform the last 11 levels. */ \
\

#define hg_pass1levels5pfa5(procname, xpass2name, pass2_levels, clm) \
/* (fftname, pass2_levels, clmarg) */ \
\
/* set_FFT_constants(fftname, 1, 0, pass2_levels, clmarg); */ \
xfft_header_hg1(procname, xpass2name, pass2_levels, clm); \
\
/* Do the last 11 FFT levels */ \
/*	do 1 pass2_11_levels_real and 4 pass2_11_levels_complex */ \
\
pass2:	start_timer(0); /* pass1_inverse_fft_setup */ \
xpass2name(g); \
end_timer(0); \
if(g->ffttype == 1)		/* We're done if FFTing only */ \
	return; \
	\
\
\
/* Do the inverse FFT */ \
\
pass1_inverse_fft_setup(c0b, procname); \
\
\
\
/* Do inverse FFT levels 4,5 */ \
/* On input the 64-byte cache lines hold these data values: */ \
/*	0K	2K	4K	6K */ \
/*	8K	10K	12K	14K */ \
/*	... */ \
/* On output the 64-byte cache lines hold these data values: */ \
/*	0K	2K	4K	6K */ \
/*	8K	12K	16K	20K */ \
/*	... */ \
\
/* Do 1 eight_reals_unfft_2 macros */ \
/*	distance between fft data elements is 2K */ \
/*	do 1 macros each processing 8 data values */ \
\
c0b:	start_timer(18); \
rsi = (uintptr_t)g->data_addr;/* Get source address */ \
/*rax = 0;*/ \
/*rdx = 0;*/ \
pfing(rbp = (uintptr_t)g->data_prefetch);/* Set prefetch distance */ \
pfing(rbp -= rsi); \
for(unsigned int loopA = 2*clm; loopA; loopA--) { \
	x1cl_half_eight_reals_unfft_2(rsi, 64); \
} \
\
/* Do 4 four_complex_unfft macros */ \
/*	distance between fft data elements is 2K */ \
/*	do 2 macros each processing 16 data values */ \
\
rsi += -2*clm*64+blkdst;	/* Next source pointer */ \
rdi = (uintptr_t)g->sincos1;/* Load sin/cos pointer */ \
/* 2 iterations of 2*clm */ \
for(unsigned int loopA = 2; loopA; loopA--) { \
	for(unsigned int loopB = 2*clm; loopB; loopB--) { \
		x2cl_four_complex_unfft(rsi, 64, blkdst); \
	} \
	rsi += -2*clm*64+2*blkdst;	/* Next source pointer */ \
	rdi += XMM_SCD;		/* Next sine/cosine pointer */ \
} \
rsi += -2*2*blkdst-blkdst;	/* Restore source pointer */ \
end_timer(18); \
\
/* Do inverse FFT levels 1,2,3 */ \
/* On input the 64-byte cache lines hold these data values: */ \
/*	0K	2K	4K	6K */ \
/*	8K	12K	16K	20K */ \
/*	... */ \
/* On output the 64-byte cache lines hold these data values: */ \
/*	0K	2K	20K	22K */ \
/*	4K	12K	24K	32K */ \
/*	6K */ \
/*	... */ \
/*	10K */ \
\
/* Do 8 five_reals_last_unfft macros */ \
/*	distance between fft data elements is 8K */ \
/*	do 1 macros each processing 40 data values */ \
\
start_timer(19); \
xfive_reals_unfft_preload;	/* Preload consts for five real macros */ \
for(unsigned int loopA = clm*2; loopA; ) { \
	do{ \
		x5cl_five_reals_last_unfft(rsi, 64, blkdst); \
		L2prefetch128(u8ptr(rsi-64+rbp)); \
		pfing(rbp += -64+blkdst); \
	}while(--loopA & 1); \
	pfing(rbp += -2*blkdst+2*64);	 /* Restore prefetching offset */ \
} \
end_timer(19); \
\
/* Normalize these values */ \
\
pass1_normalize(0, b0b); \
goto	endlp; \
\
/* Do FFT levels 1,2,3 */ \
/* */ \
/* On input the 64-byte cache lines hold these data values: */ \
/*	0K	2K	20K	22K */ \
/*	4K	12K	24K	32K */ \
/*	6K */ \
/*	... */ \
/*	10K */ \
/* On output the 64-byte cache lines hold these data values: */ \
/*	0K	2K	4K	6K */ \
/*	8K	16K	24K	32K */ \
/*	10K */ \
/*	... */ \
/*	14K */ \
\
/* Fire up auxiliary threads to do pass 1 forward FFT */ \
\
pass1:	/*pass1_forward_fft_setup(b0b);*/ \
\
/* Do 8 five_reals_first_fft macros (FFT levels 1,2,3) */ \
/*	distance between fft data elements is 8K */ \
/*	do 1 iterations each processing 40 data values */ \
\
b0b:	start_timer(21); \
rsi = (uintptr_t)g->data_addr;/* Load source address for prefetching */ \
rbx = g->DIST_TO_FFTSRCARG; \
pfing(rbp = (uintptr_t)g->data_prefetch);/* Set prefetch distance */ \
pfing(rbp -= rsi); \
/*rax = 0;*/		/* 1 iteration of 2*clm */ \
/*rdx = 0;*/ \
xfive_reals_fft_preload;		/* Preload consts for five real macros */ \
for(unsigned int loopA = clm*2; loopA; ) { \
	do{ \
		x5cl_five_reals_first_fft(rsi, 64, blkdst); \
		L2prefetch128(u8ptr(rsi-64+2*blkdst+rbp)); \
		pfing(rbp += -64+blkdst); \
	}while(--loopA & 1); \
	pfing(rbp += -2*blkdst+2*64);	 /* Restore prefetching offset */ \
} \
rsi += -2*clm*64;		/* Restore source pointer */ \
end_timer(21); \
\
/* Do FFT levels 4,5 */ \
/* Values 0K-7K is real data, 16K-39K is complex data. */ \
/* */ \
/* On input the 64-byte cache lines hold these data values: */ \
/*	0K	8K	16K	24K */ \
/*	2K	... */ \
/*	4K	... */ \
/*	6K	... */ \
/*	32K	... */ \
/*	... */ \
/* On output the 64-byte cache lines hold these data values: */ \
/*	0K	2K	4K	6K */ \
/*	8K	... */ \
/*	... */ \
\
/* Do 1 eight_reals_fft_2 macros */ \
/*	distance between fft data elements is 2K */ \
/*	do 1 macros each processing 8 data values */ \
\
start_timer(22); \
for(unsigned int loopA = clm*2; loopA; ) { \
	do{ \
		x1cl_half_eight_reals_fft_2(rsi, 64); \
	}while(--loopA & 1); \
	L2prefetch128(u8ptr(rsi-2*64+rbp)); \
} \
\
/* Do 4 four_complex_fft macros */ \
/*	distance between fft data elements is 2K */ \
/*	do 1 macros each processing 32 data values */ \
\
rsi += -2*clm*64+blkdst;	/* Next source pointer */ \
rdi = (uintptr_t)g->sincos1;/* Load sin/cos pointer */ \
for(unsigned int loopA = clm*2; loopA; ) { \
	do{ \
		x4cl_four_complex_fft(rsi, 64, blkdst, 2*blkdst, rdi); \
		L2prefetch128(u8ptr(rsi-64+rbp)); \
		L2prefetch128(u8ptr(rsi-64+blkdst+rbp)); \
		pfing(rbp += -64+2*blkdst); \
	}while(--loopA & 1); \
	pfing(rbp += -4*blkdst+2*64);	 /* Restore prefetching offset */ \
} \
end_timer(22); \
\
/* Work on next 4 sets of 20 values */ \
\
endlp: 	xpass1_get_next_block(pass2, c0b, b0b); /* Get next block and jump */ \
done:	xfft_footer_hg1;  \

\
/* Perform a 48K-element FFT.  This is done in two passes.  Pass 1 */ \
/* does 5 levels, pass 2 uses common code to perform the last 11 levels. */ \
\

#define hg_pass1levels5pfa6(procname, xpass2name, pass2_levels, clm) \
/* (fftname, pass2_levels, clmarg) */ \
\
/* set_FFT_constants(fftname, 1, 0, pass2_levels, clmarg); */ \
xfft_header_hg1(procname, xpass2name, pass2_levels, clm); \
\
/* Do the last 11 FFT levels */ \
/*	do 1 pass2_11_levels_real and 5 pass2_11_levels_complex */ \
\
pass2:	start_timer(0); \
xpass2name(g); \
end_timer(0); \
if(g->ffttype == 1)		/* We're done if FFTing only */ \
	return; /* goto done; */ \
\
/* Do the inverse FFT */ \
\
pass1_inverse_fft_setup(c0b, procname); \
\
/* Do inverse FFT levels 4,5 */ \
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
c0b:	start_timer(18); \
rsi = (uintptr_t)g->data_addr;/* Get source address */ \
/*rax = 0;*/ \
/*rdx = 0;*/ \
pfing(rbp = (uintptr_t)g->data_prefetch);/* Set prefetch distance */ \
pfing(rbp -= rsi); \
for(unsigned int loopA = 2*clm; loopA; loopA--) { \
	x2cl_eight_reals_unfft_2(rsi, 64, blkdst); \
} \
\
/* Do 4 four_complex_unfft macros */ \
/*	distance between fft data elements is 2K */ \
/*	do 2 macros each processing 16 data values */ \
\
rsi += -2*clm*64+2*blkdst;	/* Next source pointer */ \
rdi = (uintptr_t)g->sincos1;/* Load sin/cos pointer */ \
/* 2 iterations of 2*clm */ \
for(unsigned int loopA = 2; loopA; loopA--) { \
	for(unsigned int loopB = 2*clm; loopB; loopB--) { \
		x2cl_four_complex_unfft(rsi, 64, blkdst); \
	} \
	rsi += -2*clm*64+2*blkdst;	/* Next source pointer */ \
	rdi += XMM_SCD;		/* Next sine/cosine pointer */ \
} \
rsi += -3*2*blkdst;	/* Restore source pointer */ \
end_timer(18); \
\
/* Do inverse FFT levels 1,2,3 */ \
/* On input the 64-byte cache lines hold these data values: */ \
/*	0K	4K	8K	12K */ \
/*	2K	... */ \
/*	16K	... */ \
/*	... */ \
/* On output the 64-byte cache lines hold these data values: */ \
/*	0K	12K	24K	36K */ \
/*	... */ \
/*	10K */ \
\
/* Do 8 six_reals_last_unfft macros */ \
/*	distance between fft data elements is 8K */ \
/*	do 2 macros each processing 24 data values */ \
\
start_timer(19); \
/* 2 iterations of 2*clm */ \
xsix_reals_unfft_preload; \
for(unsigned int loopA = 2; loopA; loopA--) { \
	for(unsigned int loopB = clm*2; loopB; ) { \
		do{ \
			x3cl_six_reals_last_unfft(rsi, 64, 2*blkdst); \
			L2prefetch128(u8ptr(rsi-64+rbp)); \
			pfing(rbp += -64+2*blkdst); \
		}while(--loopB & 1); \
		L2prefetch128(u8ptr(rsi+rbp)); \
		pfing(rbp += -4*blkdst+2*64);	 /* Restore prefetching offset */ \
	} \
	rsi += -2*clm*64+blkdst;	/* Next source pointer */ \
} \
end_timer(19); \
\
/* Normalize these values */ \
\
pass1_normalize(0, b0b); \
goto	endlp; \
\
/* Do FFT levels 1,2,3 */ \
/* */ \
/* On input the 64-byte cache lines hold these data values: */ \
/*	0K	12K	24K	36K */ \
/*	... */ \
/*	10K */ \
/* On output the 64-byte cache lines hold these data values: */ \
/*	0K	4K	8K	12K */ \
/*	2K	... */ \
/*	16K	... */ \
/*	... */ \
\
/* Fire up auxiliary threads to do pass 1 forward FFT */ \
\
pass1:	pass1_forward_fft_setup(b0b); \
\
/* Do 8 six_reals_first_fft macros (FFT levels 1,2,3) */ \
/*	distance between fft data elements is 8K */ \
/*	do 2 iterations each processing 24 data values */ \
\
b0b:	start_timer(21); \
rsi = (uintptr_t)g->data_addr;/* Load source address for prefetching */ \
rbx = g->DIST_TO_FFTSRCARG; \
pfing(rbp = (uintptr_t)g->data_prefetch);/* Set prefetch distance */ \
pfing(rbp -= rsi); \
/* 2 iterations of 2*clm */ \
/*rdx = 0;*/ \
xsix_reals_fft_preload; \
for(unsigned int loopA = 2; loopA; loopA--) { \
	for(unsigned int loopB = clm*2; loopB; ) { \
		do{ \
			x3cl_six_reals_first_fft(rsi, 64, 2*blkdst); \
			L2prefetch128(u8ptr(rsi-64+rbp)); \
			pfing(rbp += -64+2*blkdst); \
		}while(--loopB & 1); \
		L2prefetch128(u8ptr(rsi+rbp)); \
		pfing(rbp += -4*blkdst+2*64);	/* Restore prefetching offset */ \
	} \
	rsi += -2*clm*64+blkdst;	/* Next source pointer */ \
} \
rsi += -2*blkdst;		/* Restore source pointer */ \
end_timer(21); \
\
/* Do FFT levels 4,5 */ \
/* Values 0K-7K is real data, 8K-15K is semi-real data, 16K-47K is */ \
/* complex data. */ \
/* */ \
/* On input the 64-byte cache lines hold these data values: */ \
/*	0K	4K	8K	12K */ \
/*	2K	... */ \
/*	16K	24K	32K	40K */ \
/*	18K */ \
/*	... */ \
/*	22K */ \
/* On output the 64-byte cache lines hold these data values: */ \
/*	0K	2K	4K	6K */ \
/*	8K	... */ \
/*	... */ \
\
/* Do 1 eight_reals_fft_2 macros */ \
/* Do 1 nop_two_two_complex_fft_2 macros */ \
/*	distance between fft data elements is 2K */ \
/*	do 1 macros each processing 16 data values */ \
\
start_timer(22); \
for(unsigned int loopA = clm*2; loopA; ) { \
	do{ \
		x2cl_eight_reals_fft_2(rsi, 64, blkdst); \
		L2prefetch128(u8ptr(rsi-64+rbp)); \
		pfing(rbp += -64+blkdst); \
	}while(--loopA & 1); \
	pfing(rbp += -2*blkdst+2*64);	/* Restore prefetching offset */ \
} \
\
/* Do 4 four_complex_fft macros */ \
/*	distance between fft data elements is 2K */ \
/*	do 1 macros each processing 32 data values */ \
\
rsi += -2*clm*64+2*blkdst;	/* Next source pointer */ \
rdi = (uintptr_t)g->sincos1;/* Load sin/cos pointer */ \
for(unsigned int loopA = clm*2; loopA; ) { \
	do{ \
		x4cl_four_complex_fft(rsi, 64, blkdst, 2*blkdst, rdi); \
		L2prefetch128(u8ptr(rsi-64+rbp)); \
		L2prefetch128(u8ptr(rsi-64+blkdst+rbp)); \
		pfing(rbp += -64+2*blkdst); \
	}while(--loopA & 1); \
	pfing(rbp += -4*blkdst+2*64);	/* Restore prefetching offset */ \
} \
end_timer(22); \
\
/* Work on next 4 sets of 24 values */ \
\
endlp:	xpass1_get_next_block(pass2, c0b, b0b); /* Get next block and jump */ \
done:	xfft_footer_hg1; \
\
\
\
/* Perform a 48K-element negacyclic FFT.  This is done in two passes.  Pass 1 */ \
/* does 5 levels, pass 2 uses common code to perform the last 11 levels. */ \
\

#define hg_pass1levels5complex3(procname, xpass2name, pass2_levels, clm) \
/* (fftname, pass2_levels, clmarg) */ \
\
/* set_FFT_constants(fftname, 1, 1, pass2_levels, clmarg); */ \
xfft_header_hg1(procname, xpass2name, pass2_levels, clm); \
\
/* Do the last 11 FFT levels */ \
/*	do 6 pass2_11_levels_complex */ \
\
pass2:	start_timer(0); \
xpass2name(g); \
end_timer(0); \
if(g->ffttype == 1)		/* We're done if FFTing only */ \
	return; /* goto done; */ \
\
/* Do the inverse FFT */ \
\
pass1_inverse_fft_setup(c0b, procname); \
\
/* Do inverse FFT levels 4,5 */ \
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
/* Do 6 four_complex_unfft macros */ \
/*	distance between fft data elements is 2K */ \
/*	do 3 macros each processing 16 data values */ \
\
c0b:	start_timer(18); \
rsi = (uintptr_t)g->data_addr;/* Get source address */ \
rdi = (uintptr_t)g->sincos1;/* Load sin/cos pointer */ \
/* 3 iterations of 2*clm */ \
/* rdx = 0 */ \
pfing(rbp = (uintptr_t)g->data_prefetch);/* Set prefetch distance */ \
pfing(rbp -= rsi); \
for(unsigned int loopA = 3; loopA; loopA--) { \
	for(unsigned int loopB = 2*clm; loopB; loopB--) { \
		x2cl_four_complex_unfft(rsi, 64, blkdst); \
	} \
	rsi += -2*clm*64+2*blkdst;	/* Next source pointer */ \
	rdi += XMM_SCD;		/* Next sine/cosine pointer */ \
} \
rsi += -3*2*blkdst;	/* Restore source pointer */ \
end_timer(18); \
\
/* Do inverse FFT levels 1,2,3 */ \
/* On input the 64-byte cache lines hold these data values: */ \
/*	0K	4K	8K	12K */ \
/*	2K	... */ \
/*	16K	... */ \
/*	... */ \
/* On output the 64-byte cache lines hold these data values: */ \
/*	0K	12K	24K	36K */ \
/*	... */ \
/*	10K */ \
\
/* Do 8 three_complex_last_unfft macros */ \
/*	distance between fft data elements is 8K */ \
/*	do 2 macros each processing 24 data values */ \
\
start_timer(19); \
rdi = (uintptr_t)g->plus1_premults;/* Address of premultiplier table */ \
/* 2 iterations of 2*clm */ \
for(unsigned int loopA = 2; loopA; loopA--) { \
	for(unsigned int loopB = clm*2; loopB; ) { \
		do{ \
			x3cl_three_complex_last_unfft(rsi, 64, 2*blkdst); \
			L2prefetch128(u8ptr(rsi-64+rbp)); \
			pfing(rbp += -64+2*blkdst); \
		}while(--loopB & 1); \
		L2prefetch128(u8ptr(rsi+rbp)); \
		pfing(rbp += -4*blkdst+2*64);	/* Restore prefetching offset */ \
	} \
	rsi += -2*clm*64+blkdst;	/* Next source pointer */ \
	rdi += 192;		/* Next premultiplier pointer */ \
} \
end_timer(19); \
\
/* Normalize these values */ \
\
pass1_normalize(0, b0b); \
goto	endlp; \
\
/* Do FFT levels 1,2,3 */ \
/* */ \
/* On input the 64-byte cache lines hold these data values: */ \
/*	0K	12K	24K	36K */ \
/*	... */ \
/*	10K */ \
/* On output the 64-byte cache lines hold these data values: */ \
/*	0K	4K	8K	12K */ \
/*	2K	... */ \
/*	16K	... */ \
/*	... */ \
\
/* Fire up auxiliary threads to do pass 1 forward FFT */ \
\
pass1:	pass1_forward_fft_setup(b0b); \
\
/* Do 8 three_complex_first_fft macros (FFT levels 1,2,3) */ \
/*	distance between fft data elements is 8K */ \
/*	do 2 iterations each processing 24 data values */ \
\
b0b:	start_timer(21); \
rsi = (uintptr_t)g->data_addr;/* Load source address for prefetching */ \
rbx = g->DIST_TO_FFTSRCARG; \
pfing(rbp = (uintptr_t)g->data_prefetch);/* Set prefetch distance */ \
pfing(rbp -= rsi); \
rdi = (uintptr_t)g->plus1_premults;/* Address of premultiplier table */ \
/* 2 iterations of 2*clm */ \
/* rdx = 0 */ \
for(unsigned int loopA = 2; loopA; loopA--) { \
	for(unsigned int loopB = clm*2; loopB; ) { \
		do{ \
			x3cl_three_complex_first_fft(rsi, 64, 2*blkdst); \
			L2prefetch128(u8ptr(rsi-64+rbp)); \
			pfing(rbp += -64+2*blkdst); \
		}while(--loopB & 1); \
		L2prefetch128(u8ptr(rsi+rbp)); \
		pfing(rbp += -4*blkdst+2*64);	/* Restore prefetching offset */ \
	} \
	rsi += -2*clm*64+blkdst;	/* Next source pointer */ \
	rdi += 192;		/* Next premultiplier pointer */ \
} \
rsi += -2*blkdst;		/* Restore source pointer */ \
end_timer(21); \
\
/* Do FFT levels 4,5 */ \
/* */ \
/* On input the 64-byte cache lines hold these data values: */ \
/*	0K	4K	8K	12K */ \
/*	2K	... */ \
/*	16K	24K	32K	40K */ \
/*	18K */ \
/*	... */ \
/*	22K */ \
/* On output the 64-byte cache lines hold these data values: */ \
/*	0K	2K	4K	6K */ \
/*	8K	... */ \
/*	... */ \
\
/* Do 6 four_complex_fft macros */ \
/*	distance between fft data elements is 2K */ \
/*	do 1 macros each processing 16 data values */ \
/*	do 1 macros each processing 32 data values */ \
\
start_timer(22); \
rdi = (uintptr_t)g->sincos1;/* Load sin/cos pointer */ \
for(unsigned int loopA = clm*2; loopA; ) { \
	do{ \
		x2cl_four_complex_fft(rsi, 64, blkdst); \
		L2prefetch128(u8ptr(rsi-64+rbp)); \
		pfing(rbp += -64+blkdst); \
	}while(--loopA & 1); \
	pfing(rbp += -2*blkdst+2*64);	/* Restore prefetching offset */ \
} \
rsi += -2*clm*64+2*blkdst;	/* Next source pointer */ \
rdi += XMM_SCD;		/* Next sin/cos pointer */ \
for(unsigned int loopA = clm*2; loopA; ) { \
	do{ \
		x4cl_four_complex_fft(rsi, 64, blkdst, 2*blkdst, rdi); \
		L2prefetch128(u8ptr(rsi-64+rbp)); \
		L2prefetch128(u8ptr(rsi-64+blkdst+rbp)); \
		pfing(rbp += -64+2*blkdst); \
	}while(--loopA & 1); \
	pfing(rbp += -4*blkdst+2*64);	/* Restore prefetching offset */ \
} \
end_timer(22); \
\
/* Work on next 4 sets of 24 values */ \
\
endlp:	xpass1_get_next_block(pass2, c0b, b0b); /* Get next block and jump */ \
done:	xfft_footer_hg1; \
\
\
\
/* Perform a 56K-element FFT.  This is done in two passes.  Pass 1 */ \
/* does 5 levels, pass 2 uses common code to perform the last 11 levels. */ \
\

#define hg_pass1levels5pfa7(procname, xpass2name, pass2_levels, clm) \
/* (fftname, pass2_levels, clmarg) */ \
\
/* set_FFT_constants(fftname, 1, 0, pass2_levels, clmarg); */ \
xfft_header_hg1(procname, xpass2name, pass2_levels, clm); \
\
/* Do the last 11 FFT levels */ \
/*	do 1 pass2_11_levels_real and 6 pass2_11_levels_complex */ \
\
pass2:	start_timer(0); \
xpass2name(g); \
end_timer(0); \
if(g->ffttype == 1)		/* We're done if FFTing only */ \
	return; /* goto done; */ \
\
/* Do the inverse FFT */ \
\
pass1_inverse_fft_setup(c0b, procname); \
\
/* Do inverse FFT levels 4,5 */ \
/* On input the 64-byte cache lines hold these data values: */ \
/*	0K	2K	4K	6K */ \
/*	8K	10K	12K	14K */ \
/*	... */ \
/* On output the 64-byte cache lines hold these data values: */ \
/*	0K	2K	4K	6K */ \
/*	8K	12K	16K	20K */ \
/*	10K */ \
/*	24K */ \
/*	... */ \
\
/* Do 1 eight_reals_unfft_2 macros */ \
/*	distance between fft data elements is 2K */ \
/*	do 1 macros each processing 16 data values */ \
\
c0b:	start_timer(18); \
rsi = (uintptr_t)g->data_addr;/* Get source address */ \
/*rax = 0;*/ \
/*rdx = 0;*/ \
pfing(rbp = (uintptr_t)g->data_prefetch);/* Set prefetch distance */ \
pfing(rbp -= rsi); \
for(unsigned int loopA = 2*clm; loopA; loopA--) { \
	x1cl_half_eight_reals_unfft_2(rsi, 64); \
} \
\
/* Do 6 four_complex_unfft macros */ \
/*	distance between fft data elements is 2K */ \
/*	do 3 macros each processing 16 data values */ \
\
rsi += -2*clm*64+blkdst;	/* Next source pointer */ \
rdi = (uintptr_t)g->sincos1;/* Load sin/cos pointer */ \
/* 3 iterations of 2*clm */ \
for(unsigned int loopA = 3; loopA; loopA--) { \
	for(unsigned int loopB = 2*clm; loopB; loopB--) { \
		x2cl_four_complex_unfft(rsi, 64, blkdst); \
	} \
	rsi += -2*clm*64+2*blkdst;	/* Next source pointer */ \
	rdi += XMM_SCD;		/* Next sine/cosine pointer */ \
} \
rsi += -3*2*blkdst-blkdst;	/* Restore source pointer */ \
end_timer(18); \
\
/* Do inverse FFT levels 1,2,3 */ \
/* On input the 64-byte cache lines hold these data values: */ \
/*	0K	2K	4K	6K */ \
/*	8K	12K	16K	20K */ \
/*	10K */ \
/*	24K */ \
/*	... */ \
/* On output the 64-byte cache lines hold these data values: */ \
/*	0K	2K	28K	30K */ \
/*	4K	8K	32K	36K */ \
/*	6K */ \
/*	12K	20K	40K	48K */ \
/*	14K */ \
/*	... */ \
/*	18K */ \
\
/* Do 8 seven_reals_last_unfft macros */ \
/*	distance between fft data elements is 8K */ \
/*	do 1 macros each processing 56 data values */ \
\
start_timer(19); \
xseven_reals_unfft_preload; \
for(unsigned int loopA = clm*2; loopA; ) { \
	do{ \
		x7cl_seven_reals_last_unfft(rsi, 64, blkdst); \
		L2prefetch128(u8ptr(rsi-64+rbp)); \
		L2prefetch128(u8ptr(rsi-64+blkdst+rbp)); \
		pfing(rbp += -64+2*blkdst); \
	}while(--loopA & 1); \
	pfing(rbp += -4*blkdst+2*64);	/* Restore prefetching offset */ \
} \
end_timer(19); \
\
/* Normalize these values */ \
\
pass1_normalize(0, b0b); \
goto	endlp; \
\
/* Do FFT levels 1,2,3 */ \
/* */ \
/* On input the 64-byte cache lines hold these data values: */ \
/*	0K	2K	28K	30K */ \
/*	4K	8K	32K	36K */ \
/*	6K */ \
/*	12K	20K	40K	48K */ \
/*	14K */ \
/*	... */ \
/*	18K */ \
/* On output the 64-byte cache lines hold these data values: */ \
/*	0K	2K	4K	6K */ \
/*	8K	12K	16K	20K */ \
/*	10K */ \
/*	24K	32K	40K	48K */ \
/*	26K */ \
/*	... */ \
/*	30K */ \
\
/* Fire up auxiliary threads to do pass 1 forward FFT */ \
\
pass1:	pass1_forward_fft_setup(b0b); \
\
/* Do 8 seven_reals_first_fft macros (FFT levels 1,2,3) */ \
/*	distance between fft data elements is 8K */ \
/*	do 1 iterations each processing 56 data values */ \
\
b0b:	start_timer(21); \
rsi = (uintptr_t)g->data_addr;/* Load source address for prefetching */ \
rbx = g->DIST_TO_FFTSRCARG; \
pfing(rbp = (uintptr_t)g->data_prefetch);/* Set prefetch distance */ \
pfing(rbp -= rsi); \
/* rax = 0 */	/* 1 iteration of 2*clm */ \
/* rdx = 0 */ \
xseven_reals_fft_preload; \
for(unsigned int loopA = clm*2; loopA; ) { \
	do{ \
		x7cl_seven_reals_first_fft(rsi, 64, blkdst); \
		L2prefetch128(u8ptr(rsi-64+3*blkdst+rbp)); \
		L2prefetch128(u8ptr(rsi-64+4*blkdst+rbp)); \
		pfing(rbp += -64+2*blkdst); \
	}while(--loopA & 1); \
	pfing(rbp += -4*blkdst+2*64);	/* Restore prefetching offset */ \
} \
rsi += -2*clm*64;		/* Restore source pointer */ \
end_timer(21); \
\
/* Do FFT levels 4,5 */ \
/* Values 0K-7K is real data, 8K-15K is semi-real data, 16K-55K is */ \
/* complex data. */ \
/* */ \
/* On input the 64-byte cache lines hold these data values: */ \
/*	0K	2K	4K	6K */ \
/*	8K	12K	16K	20K */ \
/*	10K */ \
/*	24K	32K	40K	48K */ \
/*	26K */ \
/*	... */ \
/*	30K */ \
/* On output the 64-byte cache lines hold these data values: */ \
/*	0K	2K	4K	6K */ \
/*	8K	... */ \
/*	... */ \
\
/* Do 1 eight_reals_fft_2 macros */ \
/*	distance between fft data elements is 2K */ \
/*	do 1 macros each processing 8 data values */ \
\
start_timer(22); \
for(unsigned int loopA = clm*2; loopA; ) { \
	do{ \
		x1cl_half_eight_reals_fft_2(rsi, 64); \
	}while(--loopA & 1); \
	L2prefetch128(u8ptr(rsi-128+rbp)); \
} \
\
/* Do 2 four_complex_fft macros */ \
/*	distance between fft data elements is 2K */ \
/*	do 1 macros each processing 16 data values */ \
\
rsi += -2*clm*64+blkdst;	/* Next source pointer */ \
rdi = (uintptr_t)g->sincos1;/* Load sin/cos pointer */ \
for(unsigned int loopA = clm*2; loopA; ) { \
	do{ \
		x2cl_four_complex_fft(rsi, 64, blkdst); \
		L2prefetch128(u8ptr(rsi-64+rbp)); \
		pfing(rbp += -64+blkdst); \
	}while(--loopA & 1); \
	pfing(rbp += -2*blkdst+2*64);	/* Restore prefetching offset */ \
} \
\
/* Do 4 four_complex_fft macros */ \
/*	distance between fft data elements is 2K */ \
/*	do 1 macros each processing 32 data values */ \
\
rsi += -2*clm*64+2*blkdst;	/* Next source pointer */ \
rdi += XMM_SCD;		/* Next sine/cosine pointer */ \
for(unsigned int loopA = clm*2; loopA; ) { \
	do{ \
		x4cl_four_complex_fft(rsi, 64, blkdst, 2*blkdst, rdi); \
		L2prefetch128(u8ptr(rsi-64+rbp)); \
		L2prefetch128(u8ptr(rsi-64+blkdst+rbp)); \
		pfing(rbp += -64+2*blkdst); \
	}while(--loopA & 1); \
	pfing(rbp += -4*blkdst+2*64);	/* Restore prefetching offset */ \
} \
end_timer(22); \
\
/* Work on next 4 sets of 28 values */ \
\
endlp:	xpass1_get_next_block(pass2, c0b, b0b); /* Get next block and jump */ \
done:	xfft_footer_hg1; \
\
\
\
/* Perform a 64K-element FFT.  This is done in two passes.  Pass 1 */ \
/* does 5 levels, pass 2 uses common code to perform the last 11 levels. */ \
\

#define hg_pass1levels5pfa8(procname, xpass2name, pass2_levels, clm) \
/* (fftname, pass2_levels, clmarg) */ \
\
/* set_FFT_constants(fftname, 1, 0, pass2_levels, clmarg); */ \
xfft_header_hg1(procname, xpass2name, pass2_levels, clm); \
\
/* Do the last 11 FFT levels */ \
/*	do 1 pass2_11_levels_real and 7 pass2_11_levels_complex */ \
\
pass2:	start_timer(0); \
xpass2name(g); \
end_timer(0); \
if(g->ffttype == 1)		/* We're done if FFTing only */ \
	return; /* goto done; */ \
\
/* Do the inverse FFT */ \
\
pass1_inverse_fft_setup(c0b, procname); \
\
/* Do inverse FFT levels 4,5 */ \
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
c0b:	start_timer(18); \
rsi = (uintptr_t)g->data_addr;/* Get source address */ \
/* rax = 0 */ \
/* rdx = 0 */ \
pfing(rbp = (uintptr_t)g->data_prefetch);/* Set prefetch distance */ \
pfing(rbp -= rsi); \
for(unsigned int loopA = 2*clm; loopA; loopA--) { \
	x2cl_eight_reals_unfft_2(rsi, 64, blkdst); \
} \
\
/* Do 6 four_complex_unfft macros */ \
/*	distance between fft data elements is 2K */ \
/*	do 3 macros each processing 16 data values */ \
\
rsi += -2*clm*64+2*blkdst;	/* Next source pointer */ \
rdi = (uintptr_t)g->sincos1;/* Load sin/cos pointer */ \
/* 3 iterations of 2*clm */ \
for(unsigned int loopA = 3; loopA; loopA--) { \
	for(unsigned int loopB = 2*clm; loopB; loopB--) { \
		x2cl_four_complex_unfft(rsi, 64, blkdst); \
	} \
	rsi += -2*clm*64+2*blkdst;	/* Next source pointer */ \
	rdi += XMM_SCD;		/* Next sine/cosine pointer */ \
} \
rsi += -4*2*blkdst;	/* Restore source pointer */ \
end_timer(18); \
\
/* Do inverse FFT levels 1,2,3 */ \
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
\
/* Do 8 eight_reals_last_unfft macros */ \
/*	distance between fft data elements is 8K */ \
/*	do 2 macros each processing 32 data values */ \
\
start_timer(19); \
/* 2 iterations of 2*clm */ \
for(unsigned int loopA = 2; loopA; loopA--) { \
	for(unsigned int loopB = clm*2; loopB; ) { \
		do{ \
			x4cl_eight_reals_last_unfft(rsi, 64, 2*blkdst, 4*blkdst); \
			L2prefetch128(u8ptr(rsi-64+rbp)); \
			L2prefetch128(u8ptr(rsi-64+2*blkdst+rbp)); \
			pfing(rbp += -64+4*blkdst); \
		}while(--loopB & 1); \
		pfing(rbp += -8*blkdst+2*64);	/* Restore prefetching offset */ \
	} \
	rsi += -2*clm*64+blkdst;	/* Next source pointer */ \
} \
end_timer(19); \
\
/* Normalize these values */ \
\
pass1_normalize(0, b0b); \
goto	endlp; \
\
/* Do FFT levels 1,2,3 */ \
/* */ \
/* On input the 64-byte cache lines hold these data values: */ \
/*	0K	16K	32K	48K */ \
/*	2K	... */ \
/*	... */ \
/*	14K */ \
/* On output the 64-byte cache lines hold these data values: */ \
/*	0K	8K	16K	24K */ \
/*	2K	... */ \
/*	4K	... */ \
/*	6K	... */ \
/*	32K	... */ \
/*	... */ \
\
/* Fire up auxiliary threads to do pass 1 forward FFT */ \
\
pass1:	pass1_forward_fft_setup(b0b); \
\
/* Do 8 eight_reals_first_fft macros (FFT levels 1,2,3) */ \
/*	distance between fft data elements is 8K */ \
/*	do 4 iterations each processing 16 data values */ \
\
b0b:	start_timer(21); \
rsi = (uintptr_t)g->data_addr;/* Load source address for prefetching */ \
rbx = g->DIST_TO_FFTSRCARG; \
pfing(rbp = (uintptr_t)g->data_prefetch);/* Set prefetch distance */ \
pfing(rbp -= rsi); \
/* eax = 4; */			/* 4 iterations of 2*clm */ \
/* rdx = 0 */ \
for(unsigned int loopA = 4; loopA; loopA--) { \
	for(unsigned int loopB = clm*2; loopB; ) { \
		do{ \
			x2cl_eight_reals_first_fft(rsi, 64, 4*blkdst); \
			L2prefetch128(u8ptr(rsi-64+rbp)); \
			pfing(rbp += -64+4*blkdst); \
		}while(--loopB & 1); \
		pfing(rbp += -8*blkdst+2*64);	/* Restore prefetching offset */ \
	} \
	rsi += -2*clm*64+blkdst;	/* Next source pointer */ \
} \
rsi += -4*blkdst;		/* Restore source pointer */ \
end_timer(21); \
\
/* Do FFT levels 4,5 */ \
/* Values 0K-7K is real data, 8K-15K is semi-real data, 16K-63K is */ \
/* complex data. */ \
/* */ \
/* On input the 64-byte cache lines hold these data values: */ \
/*	0K	8K	16K	24K */ \
/*	2K	... */ \
/*	4K	... */ \
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
start_timer(22); \
rdi = (uintptr_t)g->sincos1;/* Load sin/cos pointer */ \
for(unsigned int loopA = clm*2; loopA; ) { \
	do{ \
		x4cl_eight_reals_fft_2(rsi, 64, blkdst, 2*blkdst, rdi); \
		L2prefetch128(u8ptr(rsi-64+rbp)); \
		L2prefetch128(u8ptr(rsi-64+blkdst+rbp)); \
		pfing(rbp += -64+2*blkdst); \
	}while(--loopA & 1); \
	pfing(rbp += -4*blkdst+2*64);	/* Restore prefetching offset */ \
} \
\
/* Do 2 four_complex_fft macros */ \
/*	distance between fft data elements is 2K */ \
/*	do 1 macros each processing 32 data values */ \
\
rsi += -2*clm*64+4*blkdst;	/* Next source pointer */ \
rdi += XMM_SCD;		/* Next sine/cosine pointer */ \
for(unsigned int loopA = clm*2; loopA; ) { \
	do{ \
		x4cl_four_complex_fft(rsi, 64, blkdst, 2*blkdst, rdi); \
		L2prefetch128(u8ptr(rsi-64+rbp)); \
		L2prefetch128(u8ptr(rsi-64+blkdst+rbp)); \
		pfing(rbp += -64+2*blkdst); \
	}while(--loopA & 1); \
	pfing(rbp += -4*blkdst+2*64);	/* Restore prefetching offset */ \
} \
end_timer(22); \
\
/* Work on next 4 sets of 32 values */ \
\
endlp:	xpass1_get_next_block(pass2, c0b, b0b); /* Get next block and jump */ \
done:	xfft_footer_hg1; \
\
\
\
/* Perform a 64K-element negacyclic FFT.  This is done in two passes.  Pass 1 */ \
/* does 5 levels, pass 2 uses common code to perform the last 11 levels. */ \
\

#define hg_pass1levels5complex4(procname, xpass2name, pass2_levels, clm) \
/* (fftname, pass2_levels, clmarg) */ \
\
/* set_FFT_constants(fftname, 1, 1, pass2_levels, clmarg); */ \
xfft_header_hg1(procname, xpass2name, pass2_levels, clm); \
\
/* Do the last 11 FFT levels */ \
/*	do 8 pass2_11_levels_complex */ \
\
pass2:	start_timer(0); \
xpass2name(g); \
end_timer(0); \
if(g->ffttype == 1)		/* We're done if FFTing only */ \
	return; /* goto done; */ \
\
/* Do the inverse FFT */ \
\
pass1_inverse_fft_setup(c0b, procname); \
\
/* Do inverse FFT levels 4,5 */ \
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
/* Do 8 four_complex_unfft macros */ \
/*	distance between fft data elements is 2K */ \
/*	do 4 macros each processing 16 data values */ \
\
c0b:	start_timer(18); \
rsi = (uintptr_t)g->data_addr;/* Get source address */ \
rdi = (uintptr_t)g->sincos1;/* Load sin/cos pointer */ \
/* 4 iterations of 2*clm */ \
/* rdx = 0 */ \
pfing(rbp = (uintptr_t)g->data_prefetch);/* Set prefetch distance */ \
pfing(rbp -= rsi); \
for(unsigned int loopA = 4; loopA; loopA--) { \
	for(unsigned int loopB = 2*clm; loopB; loopB--) { \
		x2cl_four_complex_unfft(rsi, 64, blkdst); \
	} \
	rsi += -2*clm*64+2*blkdst;	/* Next source pointer */ \
	rdi += XMM_SCD;		/* Next sine/cosine pointer */ \
} \
rsi += -4*2*blkdst;	/* Restore source pointer */ \
end_timer(18); \
\
/* Do inverse FFT levels 1,2,3 */ \
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
\
/* Do 8 four_complex_last_unfft macros */ \
/*	distance between fft data elements is 8K */ \
/*	do 2 macros each processing 32 data values */ \
\
start_timer(19); \
rdi = (uintptr_t)g->plus1_premults;/* Address of premultiplier table */ \
/* 2 iterations of 2*clm */ \
for(unsigned int loopA = 2; loopA; loopA--) { \
	for(unsigned int loopB = clm*2; loopB; ) { \
		do{ \
			x4cl_four_complex_last_unfft(rsi, 64, 2*blkdst, 4*blkdst, 2*128); \
			L2prefetch128(u8ptr(rsi-64+rbp)); \
			L2prefetch128(u8ptr(rsi-64+2*blkdst+rbp)); \
			pfing(rbp += -64+4*blkdst); \
		}while(--loopB & 1); \
		pfing(rbp += -8*blkdst+2*64);	/* Restore prefetching offset */ \
	} \
	rsi += -2*clm*64+blkdst;	/* Next source pointer */ \
	rdi += 128;		/* Next premultiplier pointer */ \
} \
end_timer(19); \
\
/* Normalize these values */ \
\
pass1_normalize(0, b0b); \
goto	endlp; \
\
/* Do FFT levels 1,2,3 */ \
/* */ \
/* On input the 64-byte cache lines hold these data values: */ \
/*	0K	16K	32K	48K */ \
/*	2K	... */ \
/*	... */ \
/*	14K */ \
/* On output the 64-byte cache lines hold these data values: */ \
/*	0K	8K	16K	24K */ \
/*	2K	... */ \
/*	4K	... */ \
/*	6K	... */ \
/*	32K	... */ \
/*	... */ \
\
/* Fire up auxiliary threads to do pass 1 forward FFT */ \
\
pass1:	pass1_forward_fft_setup(b0b); \
\
/* Do 8 four_complex_first_fft macros (FFT levels 1,2,3) */ \
/*	distance between fft data elements is 8K */ \
/*	do 4 iterations each processing 16 data values */ \
\
b0b:	start_timer(21); \
rsi = (uintptr_t)g->data_addr;/* Load source address for prefetching */ \
rbx = g->DIST_TO_FFTSRCARG; \
pfing(rbp = (uintptr_t)g->data_prefetch);/* Set prefetch distance */ \
pfing(rbp -= rsi); \
rdi = (uintptr_t)g->plus1_premults;/* Address of premultiplier table */ \
/* 4 iterations of 2*clm */ \
/* rdx = 0 */ \
for(unsigned int loopA = 4; loopA; loopA--) { \
	for(unsigned int loopB = clm*2; loopB; ) { \
		do{ \
			x2cl_four_complex_first_fft(rsi, 64, 4*blkdst); \
			L2prefetch128(u8ptr(rsi-64+rbp)); \
			pfing(rbp += -64+4*blkdst); \
		}while(--loopB & 1); \
		pfing(rbp += -8*blkdst+2*64);	/* Restore prefetching offset */ \
	} \
	rsi += -2*clm*64+blkdst;	/* Next source pointer */ \
	rdi += 128;		/* Next premultiplier pointer */ \
} \
rsi += -4*blkdst;		/* Restore source pointer */ \
end_timer(21); \
\
/* Do FFT levels 4,5 */ \
/* */ \
/* On input the 64-byte cache lines hold these data values: */ \
/*	0K	8K	16K	24K */ \
/*	2K	... */ \
/*	4K	... */ \
/*	6K	... */ \
/*	32K	... */ \
/*	... */ \
/* On output the 64-byte cache lines hold these data values: */ \
/*	0K	2K	4K	6K */ \
/*	8K	... */ \
/*	... */ \
\
/* Do 8 four_complex_fft macros */ \
/*	distance between fft data elements is 2K */ \
/*	do 2 macros each processing 32 data values */ \
\
start_timer(22); \
rdi = (uintptr_t)g->sincos1;/* Load sin/cos pointer */ \
/* 2 iterations of 2*clm */ \
for(unsigned int loopA = 2; loopA; loopA--) { \
	for(unsigned int loopB = clm*2; loopB; ) { \
		do{ \
			x4cl_four_complex_fft(rsi, 64, blkdst, 2*blkdst, rdi); \
			L2prefetch128(u8ptr(rsi-64+rbp)); \
			L2prefetch128(u8ptr(rsi-64+blkdst+rbp)); \
			pfing(rbp += -64+2*blkdst); \
		}while(--loopB & 1); \
		pfing(rbp += -4*blkdst+2*64);	/* Restore prefetching offset */ \
	} \
	rsi += -2*clm*64+4*blkdst;	/* Next source pointer */ \
	rdi += 2*XMM_SCD;		/* Next sine/cosine pointer */ \
} \
end_timer(22); \
\
/* Work on next 4 sets of 32 values */ \
\
endlp:	xpass1_get_next_block(pass2, c0b, b0b); /* Get next block and jump */ \
done:	xfft_footer_hg1; \
\
\
\
/* Perform a 80K-element FFT.  This is done in two passes.  Pass 1 */ \
/* does 6 levels, pass 2 uses common code to perform the last 11 levels. */ \
\

#define hg_pass1levels6pfa5(procname, xpass2name, pass2_levels, clm) \
/* (fftname, pass2_levels, clmarg) */ \
\
/* set_FFT_constants(fftname, 1, 0, pass2_levels, clmarg); */ \
xfft_header_hg1(procname, xpass2name, pass2_levels, clm); \
\
/* Do the last 11 FFT levels */ \
/*	do 1 pass2_11_levels_real and 9 pass2_11_levels_complex */ \
\
pass2:	start_timer(0); \
xpass2name(g); \
end_timer(0); \
if(g->ffttype == 1)		/* We're done if FFTing only */ \
	return; /* goto done; */ \
\
/* Do the inverse FFT */ \
\
pass1_inverse_fft_setup(c0b, procname); \
\
/* Do inverse FFT levels 6 */ \
/* On input and output the 64-byte cache lines hold these data values: */ \
/*	0K	2K	4K	6K */ \
/*	8K	10K	12K	14K */ \
/*	... */ \
\
/* Do 1 eight_reals_unfft_1 macros */ \
/*	distance between fft data elements is 2K */ \
/*	do 1 macros each processing 8 data values */ \
\
c0b:	start_timer(16); /*puts("cob");*/ \
rsi = (uintptr_t)g->data_addr;/* Get source address */ \
/* rax = 0 */ \
/* rdx = 0 */ \
pfing(rbp = (uintptr_t)g->data_prefetch);/* Set prefetch distance */ \
pfing(rbp -= rsi); \
for(unsigned int loopA = 2*clm; loopA; loopA--) { \
	x1cl_eight_reals_unfft_1(rsi, 64); \
} \
\
/* Do 9 two_two_complex_unfft macros */ \
/*	distance between fft data elements is 2K */ \
/*	do 9 macros each processing 8 data values */ \
\
rsi += -2*clm*64+blkdst;	/* Next source pointer */ \
rdi = (uintptr_t)g->sincos1;/* Load sin/cos pointer */ \
/* eax = 9; */			/* 9 iterations of 2*clm */ \
for(unsigned int loopA = 9; loopA; loopA--) { \
	for(unsigned int loopB = 2*clm; loopB; loopB--) { \
		x1cl_two_complex_unfft(rsi, 64); \
	} \
	rsi += -2*clm*64+blkdst;	/* Next source pointer */ \
	rdi += XMM_SCD;		/* Next sine/cosine pointer */ \
} \
rsi += -10*blkdst;		/* Restore source pointer */ \
end_timer(16); \
\
/* Do inverse FFT levels 4,5 */ \
/* On input the 64-byte cache lines hold these data values: */ \
/*	0K	2K	4K	6K */ \
/*	8K	10K	12K	14K */ \
/*	... */ \
/* On output the 64-byte cache lines hold these data values: */ \
/*	0K	4K	8K	12K */ \
/*	2K	... */ \
/*	16K	24K	32K	40K */ \
/*	... */ \
\
/* Do 2 eight_reals_unfft_2 macros */ \
/*	distance between fft data elements is 4K */ \
/*	do 1 macros each processing 16 data values */ \
\
start_timer(18); \
for(unsigned int loopA = clm*2; loopA; ) { \
	do{ \
		x2cl_half_eight_reals_unfft_2(rsi, 64, blkdst); \
		L2prefetch128(u8ptr(rsi-64+rbp)); \
		pfing(rbp += -64+blkdst); \
	}while(--loopA & 1); \
	pfing(rbp += -2*blkdst+2*64);	/* Restore prefetching offset */ \
} \
\
/* Do 8 four_complex_unfft macros */ \
/*	distance between fft data elements is 4K */ \
/*	do 2 macros each processing 32 data values */ \
\
rsi += -2*clm*64+2*blkdst;	/* Next source pointer */ \
rdi = (uintptr_t)g->sincos2;/* Load sin/cos pointer */ \
/* 2 iterations of 2*clm */ \
for(unsigned int loopA = 2; loopA; loopA--) { \
	for(unsigned int loopB = clm*2; loopB; ) { \
		do{ \
			x4cl_four_complex_unfft(rsi, 64, blkdst, 2*blkdst, rdi); \
			L2prefetch128(u8ptr(rsi-64+rbp)); \
			L2prefetch128(u8ptr(rsi-64+blkdst+rbp)); \
			pfing(rbp += -64+2*blkdst); \
		}while(--loopB & 1); \
		pfing(rbp += -4*blkdst+2*64);	/* Restore prefetching offset */ \
	} \
	rsi += -2*clm*64+4*blkdst;	/* Next source pointer */ \
	rdi += XMM_SCD;		/* Next sine/cosine pointer */ \
} \
rsi += -2*4*blkdst-2*blkdst;/* Restore source pointer */ \
end_timer(18); \
\
/* Do inverse FFT levels 1,2,3 */ \
/* On input the 64-byte cache lines hold these data values: */ \
/*	0K	4K	8K	12K */ \
/*	2K	... */ \
/*	16K	24K	32K	40K */ \
/*	... */ \
/* On output the 64-byte cache lines hold these data values: */ \
/*	0K	4K	40K	44K */ \
/*	2K	... */ \
/*	8K	24K	48K	64K */ \
/*	10K */ \
/*	... */ \
/*	22K */ \
\
/* Do 16 five_reals_last_unfft macros */ \
/*	distance between fft data elements is 16K */ \
/*	do 2 macros each processing 40 data values */ \
\
start_timer(19); \
/* 2 iterations of 2*clm */ \
xfive_reals_unfft_preload;	/* Preload consts for five real macros */ \
for(unsigned int loopA = 2; loopA; loopA--) { \
	for(unsigned int loopB = clm*2; loopB; ) { \
		do{ \
			x5cl_five_reals_last_unfft(rsi, 64, 2*blkdst); \
			L2prefetch128(u8ptr(rsi-64+rbp)); \
			L2prefetch128(u8ptr(rsi-64+2*blkdst+rbp)); \
			pfing(rbp += -64+4*blkdst); \
		}while(--loopB & 1); \
		pfing(rbp += -8*blkdst+2*64);	/* Restore prefetching offset */ \
	} \
	rsi += -2*clm*64+blkdst;	/* Next source pointer */ \
} \
end_timer(19); \
\
/* Normalize these values */ \
\
pass1_normalize(0, b0b); \
goto	endlp; \
\
/* Do FFT levels 1,2,3 */ \
/* */ \
/* On input the 64-byte cache lines hold these data values: */ \
/*	0K	4K	40K	44K */ \
/*	2K	... */ \
/*	8K	24K	48K	64K */ \
/*	10K */ \
/*	... */ \
/*	22K */ \
/* On output the 64-byte cache lines hold these data values: */ \
/*	0K	4K	8K	12K */ \
/*	2K	... */ \
/*	16K	32K	48K	64K */ \
/*	18K */ \
/*	... */ \
/*	30K */ \
\
/* Fire up auxiliary threads to do pass 1 forward FFT */ \
\
pass1:	pass1_forward_fft_setup(b0b); \
\
/* Do 16 five_reals_first_fft macros (FFT levels 1,2,3) */ \
/*	distance between fft data elements is 16K */ \
/*	do 2 iterations each processing 40 data values */ \
\
b0b:	start_timer(21); /*puts("bob");*/ \
rsi = (uintptr_t)g->data_addr;/* Load source address for prefetching */ \
rbx = g->DIST_TO_FFTSRCARG; \
pfing(rbp = (uintptr_t)g->data_prefetch);/* Set prefetch distance */ \
pfing(rbp -= rsi); \
/* 2 iterations of 2*clm */ \
/* rdx = 0 */ \
xfive_reals_fft_preload;		/* Preload consts for five real macros */ \
for(unsigned int loopA = 2; loopA; loopA--) { \
	for(unsigned int loopB = clm*2; loopB; ) { \
		do{ \
			x5cl_five_reals_first_fft(rsi, 64, 2*blkdst); \
			L2prefetch128(u8ptr(rsi-64+2*blkdst+rbp)); \
			L2prefetch128(u8ptr(rsi-64+4*blkdst+rbp)); \
			pfing(rbp += -64+4*blkdst); \
		}while(--loopB & 1); \
		pfing(rbp += -8*blkdst+2*64);	/* Restore prefetching offset */ \
	} \
	rsi += -2*clm*64+blkdst;	/* Next source pointer */ \
} \
rsi += -2*blkdst;		/* Restore source pointer */ \
end_timer(21); \
\
/* Do FFT levels 4,5 */ \
/* Values 0K-15K is real data, 16K-31K is semi-real data, 32K-127K is */ \
/* complex data. */ \
/* */ \
/* On input the 64-byte cache lines hold these data values: */ \
/*	0K	4K	8K	12K */ \
/*	2K	... */ \
/*	16K	32K	48K	64K */ \
/*	18K */ \
/*	... */ \
/*	30K */ \
/* On output the 64-byte cache lines hold these data values: */ \
/*	0K	4K	8K	12K */ \
/*	2K	... */ \
/*	16K	... */ \
/*	... */ \
\
/* Do 2 eight_reals_fft_2 macros */ \
/*	distance between fft data elements is 4K */ \
/*	do 2 macros each processing 8 data values */ \
\
start_timer(22); \
/* 2 iterations of 2*clm */ \
for(unsigned int loopA = 2; loopA; loopA--) { \
	for(unsigned int loopB = clm*2; loopB; ) { \
		do{ \
			x1cl_half_eight_reals_fft_2(rsi, 64); \
		}while(--loopB & 1); \
		L2prefetch128(u8ptr(rsi-128+rbp)); \
	} \
	rsi += -2*clm*64+blkdst;	/* Next source pointer */ \
} \
\
/* Do 8 four_complex_fft macros */ \
/*	distance between fft data elements is 4K */ \
/*	do 2 macros each processing 32 data values */ \
\
rdi = (uintptr_t)g->sincos2;/* Load sin/cos pointer */ \
/* 2 iterations of 2*clm */ \
for(unsigned int loopA = 2; loopA; loopA--) { \
	for(unsigned int loopB = clm*2; loopB; ) { \
		do{ \
			x4cl_four_complex_fft(rsi, 64, 2*blkdst, 4*blkdst, rdi); \
			L2prefetch128(u8ptr(rsi-64+rbp)); \
			L2prefetch128(u8ptr(rsi-64+2*blkdst+rbp)); \
			pfing(rbp += -64+4*blkdst); \
		}while(--loopB & 1); \
		pfing(rbp += -8*blkdst+2*64);	/* Restore prefetching offset */ \
	} \
	rsi += -2*clm*64+blkdst;	/* Next source pointer */ \
} \
rsi += -2*blkdst-2*blkdst;	/* Restore source pointer */ \
end_timer(22); \
\
/* Do FFT levels 6 */ \
/* Values 0K-3K is real data, 4K-7K is semi-real data, 8K-127K is */ \
/* complex data. */ \
/* */ \
/* On input the 64-byte cache lines hold these data values: */ \
/*	0K	4K	8K	12K */ \
/*	2K	... */ \
/*	16K	... */ \
/*	... */ \
/* On output the 64-byte cache lines hold these data values: */ \
/*	0K	2K	4K	6K */ \
/*	8K	... */ \
/*	... */ \
\
/* Do 2 eight_reals_fft_1 macros */ \
/*	distance between fft data elements is 2K */ \
/*	do 1 macros each processing 16 data values */ \
\
start_timer(24); \
rdi = (uintptr_t)g->sincos1;/* Load sin/cos pointer */ \
for(unsigned int loopA = clm*2; loopA; ) { \
	do{ \
		x2cl_eight_reals_fft_1(rsi, 64, blkdst); \
	}while(--loopA & 1); \
} \
\
/* Do 8 two_two_complex_fft macros */ \
/*	distance between fft data elements is 2K */ \
/*	do 4 macros each processing 16 data values */ \
\
rsi += -2*clm*64+2*blkdst;	/* Next source pointer */ \
rdi += XMM_SCD;		/* Next sine/cosine pointer */ \
/* 4 iterations of 2*clm */ \
for(unsigned int loopA = 4; loopA; loopA--) { \
	for(unsigned int loopB = clm*2; loopB; ) { \
		do{ \
			x2cl_two_complex_fft(rsi, 64, blkdst, rdi); \
			L2prefetch128(u8ptr(rsi-64+rbp)); \
			pfing(rbp += -64+blkdst); \
		}while(--loopB & 1); \
		pfing(rbp += -2*blkdst+2*64);	/* Restore prefetching offset */ \
	} \
	rsi += -2*clm*64+2*blkdst;	/* Next source pointer */ \
	rdi += 2*XMM_SCD;		/* Next sine/cosine pointer */ \
} \
end_timer(24); \
\
/* Work on next 4 sets of 40 values */ \
\
endlp:	xpass1_get_next_block(pass2, c0b, b0b); /* Get next block and jump */ \
done:	xfft_footer_hg1; \
\
\
\
/* Perform a 96K-element FFT.  This is done in two passes.  Pass 1 */ \
/* does 6 levels, pass 2 uses common code to perform the last 11 levels. */ \
\

#define hg_pass1levels6pfa6(procname, xpass2name, pass2_levels, clm) \
/* (fftname, pass2_levels, clmarg) */ \
\
/* set_FFT_constants(fftname, 1, 0, pass2_levels, clmarg); */ \
xfft_header_hg1(procname, xpass2name, pass2_levels, clm); \
\
/* Do the last 11 FFT levels */ \
/*	do 1 pass2_11_levels_real and 11 pass2_11_levels_complex */ \
\
pass2:	start_timer(0); \
xpass2name(g); \
end_timer(0); \
if(g->ffttype == 1)		/* We're done if FFTing only */ \
	return; /* goto done; */ \
\
/* Do the inverse FFT */ \
\
pass1_inverse_fft_setup(c0b, procname); \
\
/* Do inverse FFT levels 6 */ \
/* On input and output the 64-byte cache lines hold these data values: */ \
/*	0K	2K	4K	6K */ \
/*	8K	10K	12K	14K */ \
/*	... */ \
\
/* Do 1 eight_reals_unfft_1 macros */ \
/*	distance between fft data elements is 2K */ \
/*	do 1 macros each processing 8 data values */ \
\
c0b:	start_timer(16); \
rsi = (uintptr_t)g->data_addr;/* Get source address */ \
/* rax = 0 */ \
/* rdx = 0 */ \
pfing(rbp = (uintptr_t)g->data_prefetch);/* Set prefetch distance */ \
pfing(rbp -= rsi); \
for(unsigned int loopA = 2*clm; loopA; loopA--) { \
	x1cl_eight_reals_unfft_1(rsi, 64); \
} \
\
/* Do 11 two_two_complex_unfft macros */ \
/*	distance between fft data elements is 2K */ \
/*	do 11 macros each processing 8 data values */ \
\
rsi += -2*clm*64+blkdst;	/* Next source pointer */ \
rdi = (uintptr_t)g->sincos1;/* Load sin/cos pointer */ \
/* 11 iterations of 2*clm */ \
for(unsigned int loopA = 11; loopA; loopA--) { \
	for(unsigned int loopB = 2*clm; loopB; loopB--) { \
		x1cl_two_complex_unfft(rsi, 64); \
	} \
	rsi += -2*clm*64+blkdst;	/* Next source pointer */ \
	rdi += XMM_SCD;		/* Next sine/cosine pointer */ \
} \
rsi += -12*blkdst;		/* Restore source pointer */ \
end_timer(16); \
\
/* Do inverse FFT levels 4,5 */ \
/* On input the 64-byte cache lines hold these data values: */ \
/*	0K	2K	4K	6K */ \
/*	8K	10K	12K	14K */ \
/*	... */ \
/* On output the 64-byte cache lines hold these data values: */ \
/*	0K	8K	16K	24K */ \
/*	2K	... */ \
/*	... */ \
/*	6K */ \
/*	32K	... */ \
/*	... */ \
\
/* Do 2 eight_reals_unfft_2 macros */ \
/* Do 2 nop_two_two_complex_unfft_2 macros */ \
/*	distance between fft data elements is 4K */ \
/*	do 1 macros each processing 32 data values */ \
\
start_timer(18); \
for(unsigned int loopA = clm*2; loopA; ) { \
	do{ \
		x4cl_eight_reals_unfft_2(rsi, 64, blkdst, 2*blkdst); \
		L2prefetch128(u8ptr(rsi-64+rbp)); \
		L2prefetch128(u8ptr(rsi-64+blkdst+rbp)); \
		pfing(rbp += -64+2*blkdst); \
	}while(--loopA & 1); \
	pfing(rbp += -4*blkdst+2*64);	/* Restore prefetching offset */ \
} \
\
/* Do 8 four_complex_unfft macros */ \
/*	distance between fft data elements is 4K */ \
/*	do 2 macros each processing 32 data values */ \
\
rsi += -2*clm*64+4*blkdst;	/* Next source pointer */ \
rdi = (uintptr_t)g->sincos2;/* Load sin/cos pointer */ \
/* 2 iterations of 2*clm */ \
for(unsigned int loopA = 2; loopA; loopA--) { \
	for(unsigned int loopB = clm*2; loopB; ) { \
		do{ \
			x4cl_four_complex_unfft(rsi, 64, blkdst, 2*blkdst, rdi); \
			L2prefetch128(u8ptr(rsi-64+rbp)); \
			L2prefetch128(u8ptr(rsi-64+blkdst+rbp)); \
			pfing(rbp += -64+2*blkdst); \
		}while(--loopB & 1); \
		pfing(rbp += -4*blkdst+2*64);	/* Restore prefetching offset */ \
	} \
	rsi += -2*clm*64+4*blkdst;	/* Next source pointer */ \
	rdi += XMM_SCD;		/* Next sine/cosine pointer */ \
} \
rsi += -3*4*blkdst;	/* Restore source pointer */ \
end_timer(18); \
\
/* Do inverse FFT levels 1,2,3 */ \
/* On input the 64-byte cache lines hold these data values: */ \
/*	0K	8K	16K	24K */ \
/*	2K	... */ \
/*	... */ \
/*	6K */ \
/*	32K	... */ \
/*	... */ \
/* On output the 64-byte cache lines hold these data values: */ \
/*	0K	24K	48K	72K */ \
/*	... */ \
/*	22K */ \
\
/* Do 16 six_reals_last_unfft macros */ \
/*	distance between fft data elements is 16K */ \
/*	do 4 macros each processing 24 data values */ \
\
start_timer(19); \
/* 4 iterations of 2*clm */ \
xsix_reals_unfft_preload; \
for(unsigned int loopA = 4; loopA; loopA--) { \
	for(unsigned int loopB = clm*2; loopB; ) { \
		do{ \
			x3cl_six_reals_last_unfft(rsi, 64, 4*blkdst); \
			L2prefetch128(u8ptr(rsi-64+rbp)); \
			pfing(rbp += -64+4*blkdst); \
		}while(--loopB & 1); \
		L2prefetch128(u8ptr(rsi+rbp)); \
		pfing(rbp += -8*blkdst+2*64);	/* Restore prefetching offset */ \
	} \
	rsi += -2*clm*64+blkdst;	/* Next source pointer */ \
} \
end_timer(19); \
\
/* Normalize these values */ \
\
pass1_normalize(0, b0b); \
goto	endlp; \
\
/* Do FFT levels 1,2,3 */ \
/* */ \
/* On input the 64-byte cache lines hold these data values: */ \
/*	0K	24K	48K	72K */ \
/*	... */ \
/*	22K */ \
/* On output the 64-byte cache lines hold these data values: */ \
/*	0K	8K	16K	24K */ \
/*	2K	... */ \
/*	... */ \
/*	6K */ \
/*	32K	48K	64K	80K */ \
/*	34K */ \
/*	... */ \
/*	46K */ \
\
/* Fire up auxiliary threads to do pass 1 forward FFT */ \
\
pass1:	pass1_forward_fft_setup(b0b); \
\
/* Do 16 six_reals_first_fft macros (FFT levels 1,2,3) */ \
/*	distance between fft data elements is 16K */ \
/*	do 4 iterations each processing 24 data values */ \
\
b0b:	start_timer(21); \
rsi = (uintptr_t)g->data_addr;/* Load source address for prefetching */ \
rbx = g->DIST_TO_FFTSRCARG; \
pfing(rbp = (uintptr_t)g->data_prefetch);/* Set prefetch distance */ \
pfing(rbp -= rsi); \
/* 4 iterations of 2*clm */ \
/* rdx = 0 */ \
xsix_reals_fft_preload; \
for(unsigned int loopA = 4; loopA; loopA--) { \
	for(unsigned int loopB = clm*2; loopB; ) { \
		do{ \
			x3cl_six_reals_first_fft(rsi, 64, 4*blkdst); \
			L2prefetch128(u8ptr(rsi-64+rbp)); \
			pfing(rbp += -64+4*blkdst); \
		}while(--loopB & 1); \
		L2prefetch128(u8ptr(rsi+rbp)); \
		pfing(rbp += -8*blkdst+2*64);	/* Restore prefetching offset */ \
	} \
	rsi += -2*clm*64+blkdst;	/* Next source pointer */ \
} \
rsi += -4*blkdst;		/* Restore source pointer */ \
end_timer(21); \
\
/* Do FFT levels 4,5 */ \
/* Values 0K-15K is real data, 16K-31K is semi-real data, 32K-95K is */ \
/* complex data. */ \
/* */ \
/* On input the 64-byte cache lines hold these data values: */ \
/*	0K	8K	16K	24K */ \
/*	2K	... */ \
/*	... */ \
/*	6K */ \
/*	32K	48K	64K	80K */ \
/*	34K */ \
/*	... */ \
/*	46K */ \
/* On output the 64-byte cache lines hold these data values: */ \
/*	0K	4K	8K	12K */ \
/*	2K	... */ \
/*	16K	... */ \
/*	... */ \
\
/* Do 2 eight_reals_fft_2 macros */ \
/* Do 2 nop_two_two_complex_fft_2 macros */ \
/*	distance between fft data elements is 4K */ \
/*	do 2 macros each processing 16 data values */ \
\
start_timer(22); \
/* 2 iterations of 2*clm */ \
for(unsigned int loopA = 2; loopA; loopA--) { \
	for(unsigned int loopB = clm*2; loopB; ) { \
		do{ \
			x2cl_eight_reals_fft_2(rsi, 64, 2*blkdst); \
			L2prefetch128(u8ptr(rsi-64+rbp)); \
			pfing(rbp += -64+2*blkdst); \
		}while(--loopB & 1); \
		pfing(rbp += -4*blkdst+2*64);	/* Restore prefetching offset */ \
	} \
	rsi += -2*clm*64+blkdst;	/* Next source pointer */ \
} \
\
/* Do 8 four_complex_fft macros */ \
/*	distance between fft data elements is 4K */ \
/*	do 2 macros each processing 32 data values */ \
\
rsi += -2*blkdst+4*blkdst;	/* Next source pointer */ \
rdi = (uintptr_t)g->sincos2;/* Load sin/cos pointer */ \
/* 2 iterations of 2*clm */ \
for(unsigned int loopA = 2; loopA; loopA--) { \
	for(unsigned int loopB = clm*2; loopB; ) { \
		do{ \
			x4cl_four_complex_fft(rsi, 64, 2*blkdst, 4*blkdst, rdi); \
			L2prefetch128(u8ptr(rsi-64+rbp)); \
			L2prefetch128(u8ptr(rsi-64+2*blkdst+rbp)); \
			pfing(rbp += -64+4*blkdst); \
		}while(--loopB & 1); \
		pfing(rbp += -8*blkdst+2*64);	/* Restore prefetching offset */ \
	} \
	rsi += -2*clm*64+blkdst;	/* Next source pointer */ \
} \
rsi += -2*blkdst-4*blkdst;	/* Restore source pointer */ \
end_timer(22); \
\
/* Do FFT levels 6 */ \
/* Values 0K-3K is real data, 4K-7K is semi-real data, 8K-95K is */ \
/* complex data. */ \
/* */ \
/* On input the 64-byte cache lines hold these data values: */ \
/*	0K	4K	8K	12K */ \
/*	2K	... */ \
/*	16K	... */ \
/*	... */ \
/* On output the 64-byte cache lines hold these data values: */ \
/*	0K	2K	4K	6K */ \
/*	8K	... */ \
/*	... */ \
\
/* Do 2 eight_reals_fft_1 macros */ \
/*	distance between fft data elements is 2K */ \
/*	do 1 macros each processing 16 data values */ \
\
start_timer(24); \
rdi = (uintptr_t)g->sincos1;/* Load sin/cos pointer */ \
for(unsigned int loopA = clm*2; loopA; ) { \
	do{ \
		x2cl_eight_reals_fft_1(rsi, 64, blkdst); \
	}while(--loopA & 1); \
} \
\
/* Do 10 two_two_complex_fft macros */ \
/*	distance between fft data elements is 2K */ \
/*	do 5 macros each processing 16 data values */ \
\
rsi += -2*clm*64+2*blkdst;	/* Next source pointer */ \
rdi += XMM_SCD;		/* Next sine/cosine pointer */ \
/* 5 iterations of 2*clm */ \
for(unsigned int loopA = 5; loopA; loopA--) { \
	for(unsigned int loopB = clm*2; loopB; ) { \
		do{ \
			x2cl_two_complex_fft(rsi, 64, blkdst, rdi); \
			L2prefetch128(u8ptr(rsi-64+rbp)); \
			pfing(rbp += -64+blkdst); \
		}while(--loopB & 1); \
		pfing(rbp += -2*blkdst+2*64);	/* Restore prefetching offset */ \
	} \
	rsi += -2*clm*64+2*blkdst;	/* Next source pointer */ \
	rdi += 2*XMM_SCD;		/* Next sine/cosine pointer */ \
} \
end_timer(24); \
\
/* Work on next 4 sets of 48 values */ \
\
endlp:	xpass1_get_next_block(pass2, c0b, b0b); /* Get next block and jump */ \
done:	xfft_footer_hg1; \
\
\
\
/* Perform a 96K-element negacyclic FFT.  This is done in two passes.  Pass 1 */ \
/* does 6 levels, pass 2 uses common code to perform the last 11 levels. */ \
\

#define hg_pass1levels6complex3(procname, xpass2name, pass2_levels, clm) \
/* (fftname, pass2_levels, clmarg) */ \
\
/* set_FFT_constants(fftname, 1, 1, pass2_levels, clmarg); */ \
xfft_header_hg1(procname, xpass2name, pass2_levels, clm); \
\
/* Do the last 11 FFT levels */ \
/*	do 12 pass2_11_levels_complex */ \
\
pass2:	start_timer(0); \
xpass2name(g); \
end_timer(0); \
if(g->ffttype == 1)		/* We're done if FFTing only */ \
	return; /* goto done; */ \
\
/* Do the inverse FFT */ \
\
pass1_inverse_fft_setup(c0b, procname); \
\
/* Do inverse FFT levels 5,6 */ \
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
/* Do 12 four_complex_unfft macros */ \
/*	distance between fft data elements is 2K */ \
/*	do 6 macros each processing 16 data values */ \
\
c0b:	start_timer(16); \
rsi = (uintptr_t)g->data_addr;/* Get source address */ \
rdi = (uintptr_t)g->sincos1;/* Load sin/cos pointer */ \
/* 6 iterations of 2*clm */ \
/* rdx = 0 */ \
pfing(rbp = (uintptr_t)g->data_prefetch);/* Set prefetch distance */ \
pfing(rbp -= rsi); \
for(unsigned int loopA = 6; loopA; loopA--) { \
	for(unsigned int loopB = 2*clm; loopB; loopB--) { \
		x2cl_four_complex_unfft(rsi, 64, blkdst); \
	} \
	rsi += -2*clm*64+2*blkdst;	/* Next source pointer */ \
	rdi += XMM_SCD;		/* Next sine/cosine pointer */ \
} \
rsi += -6*2*blkdst;	/* Restore source pointer */ \
end_timer(16); \
\
/* Do inverse FFT levels 4 */ \
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
\
/* Do 12 two_two_complex_unfft macros */ \
/*	distance between fft data elements is 8K */ \
/*	do 6 macros each processing 16 data values */ \
\
start_timer(18); \
rdi = (uintptr_t)g->sincos2;/* Load sin/cos pointer */ \
/* eax = 3*256+2; */ /* 3 iters of 2 iters of 2*clm */ \
for(unsigned int loopA = 3; loopA; loopA--) { \
	for(unsigned int loopB = 2; loopB; loopB--) { \
		for(unsigned int loopC = clm*2; loopC; ) { \
			do{ \
				x2cl_two_complex_unfft(rsi, 64, 2*blkdst); \
				L2prefetch128(u8ptr(rsi-64+rbp)); \
				pfing(rbp += -64+2*blkdst); \
			}while(--loopC & 1); \
			pfing(rbp += -4*blkdst+2*64);	/* Restore prefetching offset */ \
		} \
		rsi += -2*clm*64+blkdst;	/* Next source pointer */ \
	} \
	rsi += -2*blkdst+4*blkdst;	/* Next source pointer */ \
	rdi += XMM_SCD;		/* Next sine/cosine pointer */ \
} \
rsi += -3*4*blkdst;	/* Restore source pointer */ \
end_timer(18); \
\
/* Do inverse FFT levels 1,2,3 */ \
/* On input the 64-byte cache lines hold these data values: */ \
/*	0K	8K	16K	24K */ \
/*	2K	... */ \
/*	... */ \
/*	6K */ \
/*	32K	... */ \
/*	... */ \
/* On output the 64-byte cache lines hold these data values: */ \
/*	0K	24K	48K	72K */ \
/*	... */ \
/*	22K */ \
\
/* Do 8 three_complex_last_unfft macros */ \
/*	distance between fft data elements is 16K */ \
/*	do 4 macros each processing 24 data values */ \
\
start_timer(19); \
rdi = (uintptr_t)g->plus1_premults;/* Address of premultiplier table */ \
/* 4 iterations of 2*clm */ \
for(unsigned int loopA = 4; loopA; loopA--) { \
	for(unsigned int loopB = clm*2; loopB; ) { \
		do{ \
			x3cl_three_complex_last_unfft(rsi, 64, 4*blkdst); \
			L2prefetch128(u8ptr(rsi-64+rbp)); \
			pfing(rbp += -64+4*blkdst); \
		}while(--loopB & 1); \
		L2prefetch128(u8ptr(rsi+rbp)); \
		pfing(rbp += -8*blkdst+2*64);	/* Restore prefetching offset */ \
	} \
	rsi += -2*clm*64+blkdst;	/* Next source pointer */ \
	rdi += 192;		/* Next premultiplier pointer */ \
} \
end_timer(19); \
\
/* Normalize these values */ \
\
pass1_normalize(0, b0b); \
goto	endlp; \
\
/* Do FFT levels 1,2,3 */ \
/* */ \
/* On input the 64-byte cache lines hold these data values: */ \
/*	0K	24K	48K	72K */ \
/*	... */ \
/*	22K */ \
/* On output the 64-byte cache lines hold these data values: */ \
/*	0K	8K	16K	24K */ \
/*	2K	... */ \
/*	... */ \
/*	6K */ \
/*	32K	48K	64K	80K */ \
/*	34K */ \
/*	... */ \
/*	46K */ \
\
/* Fire up auxiliary threads to do pass 1 forward FFT */ \
\
pass1:	pass1_forward_fft_setup(b0b); \
\
/* Do 16 three_complex_first_fft macros (FFT levels 1,2,3) */ \
/*	distance between fft data elements is 16K */ \
/*	do 4 iterations each processing 24 data values */ \
\
b0b:	start_timer(21); \
rsi = (uintptr_t)g->data_addr;/* Load source address for prefetching */ \
rbx = g->DIST_TO_FFTSRCARG; \
pfing(rbp = (uintptr_t)g->data_prefetch);/* Set prefetch distance */ \
pfing(rbp -= rsi); \
rdi = (uintptr_t)g->plus1_premults;/* Address of premultiplier table */ \
/* 4 iterations of 2*clm */ \
/* rdx = 0 */ \
for(unsigned int loopA = 4; loopA; loopA--) { \
	for(unsigned int loopB = clm*2; loopB; ) { \
		do{ \
			x3cl_three_complex_first_fft(rsi, 64, 4*blkdst); \
			L2prefetch128(u8ptr(rsi-64+rbp)); \
			pfing(rbp += -64+4*blkdst); \
		}while(--loopB & 1); \
		L2prefetch128(u8ptr(rsi+rbp)); \
		pfing(rbp += -8*blkdst+2*64);	/* Restore prefetching offset */ \
	} \
	rsi += -2*clm*64+blkdst;	/* Next source pointer */ \
	rdi += 192;		/* Next premultiplier pointer */ \
} \
rsi += -4*blkdst;		/* Restore source pointer */ \
end_timer(21); \
\
/* Do FFT levels 4 */ \
/* */ \
/* On input the 64-byte cache lines hold these data values: */ \
/*	0K	8K	16K	24K */ \
/*	2K	... */ \
/*	... */ \
/*	6K */ \
/*	32K	48K	64K	80K */ \
/*	34K */ \
/*	... */ \
/*	46K */ \
/* On output the 64-byte cache lines hold these data values: */ \
/*	0K	8K	16K	24K */ \
/*	2K	... */ \
/*	... */ \
/*	6K */ \
/*	32K	... */ \
/*	... */ \
\
/* Do 12 two_two_complex_fft macros */ \
/*	distance between fft data elements is 8K */ \
/*	do 6 macros each processing 16 data values */ \
\
start_timer(22); \
rdi = (uintptr_t)g->sincos2;/* Load sin/cos pointer */ \
/* 2 iterations of 2*clm */ \
for(unsigned int loopA = 2; loopA; loopA--) { \
	for(unsigned int loopB = clm*2; loopB; ) { \
		do{ \
			x2cl_two_complex_fft_in_place(rsi, 64, 2*blkdst, rdi); \
			L2prefetch128(u8ptr(rsi-64+rbp)); \
			pfing(rbp += -64+2*blkdst); \
		}while(--loopB & 1); \
		pfing(rbp += -4*blkdst+2*64);	/* Restore prefetching offset */ \
	} \
	rsi += -2*clm*64+blkdst;	/* Next source pointer */ \
} \
rsi += -2*blkdst+4*blkdst;	/* Next source pointer */ \
rdi += XMM_SCD;		/* Next sin/cos pointer */ \
/* 4 iterations of 2*clm */ \
for(unsigned int loopA = 4; loopA; loopA--) { \
	for(unsigned int loopB = clm*2; loopB; ) { \
		do{ \
			x2cl_two_complex_fft(rsi, 64, 4*blkdst, rdi); \
			L2prefetch128(u8ptr(rsi-64+rbp)); \
			pfing(rbp += -64+4*blkdst); \
		}while(--loopB & 1); \
		pfing(rbp += -8*blkdst+2*64);	/* Restore prefetching offset */ \
	} \
	rsi += -2*clm*64+blkdst;	/* Next source pointer */ \
} \
rsi += -4*blkdst-4*blkdst;	/* Restore source pointer */ \
end_timer(22); \
\
/* Do FFT levels 5,6 */ \
/* */ \
/* On input the 64-byte cache lines hold these data values: */ \
/*	0K	8K	16K	24K */ \
/*	2K	... */ \
/*	... */ \
/*	6K */ \
/*	32K	... */ \
/*	... */ \
/* On output the 64-byte cache lines hold these data values: */ \
/*	0K	2K	4K	6K */ \
/*	8K	... */ \
/*	... */ \
\
/* Do 12 four_complex_fft macros */ \
/*	distance between fft data elements is 2K */ \
/*	do 3 macros each processing 32 data values */ \
\
start_timer(24); \
rdi = (uintptr_t)g->sincos1;/* Load sin/cos pointer */ \
/* 3 iterations of 2*clm */ \
for(unsigned int loopA = 3; loopA; loopA--) { \
	for(unsigned int loopB = clm*2; loopB; ) { \
		do{ \
			x4cl_four_complex_fft(rsi, 64, blkdst, 2*blkdst, rdi); \
			L2prefetch128(u8ptr(rsi-64+rbp)); \
			L2prefetch128(u8ptr(rsi-64+blkdst+rbp)); \
			pfing(rbp += -64+2*blkdst); \
		}while(--loopB & 1); \
		pfing(rbp += -4*blkdst+2*64);	/* Restore prefetching offset */ \
	} \
	rsi += -2*clm*64+4*blkdst;	/* Next source pointer */ \
	rdi += 2*XMM_SCD;		/* Next sine/cosine pointer */ \
} \
end_timer(24); \
\
/* Work on next 4 sets of 48 values */ \
\
endlp:	xpass1_get_next_block(pass2, c0b, b0b); /* Get next block and jump */ \
done:	xfft_footer_hg1; \
\
\
\
/* Perform a 112K-element FFT.  This is done in two passes.  Pass 1 */ \
/* does 6 levels, pass 2 uses common code to perform the last 11 levels. */ \
\

#define hg_pass1levels6pfa7(procname, xpass2name, pass2_levels, clm) \
/* (fftname, pass2_levels, clmarg) */ \
\
/* set_FFT_constants(fftname, 1, 0, pass2_levels, clmarg); */ \
xfft_header_hg1(procname, xpass2name, pass2_levels, clm); \
\
/* Do the last 11 FFT levels */ \
/*	do 1 pass2_11_levels_real and 13 pass2_11_levels_complex */ \
\
pass2:	start_timer(0); \
xpass2name(g); \
end_timer(0); \
if(g->ffttype == 1)		/* We're done if FFTing only */ \
	return; /* goto done; */ \
\
/* Do the inverse FFT */ \
\
pass1_inverse_fft_setup(c0b, procname); \
\
/* Do inverse FFT levels 6 */ \
/* On input and output the 64-byte cache lines hold these data values: */ \
/*	0K	2K	4K	6K */ \
/*	8K	10K	12K	14K */ \
/*	... */ \
\
/* Do 1 eight_reals_unfft_1 macros */ \
/*	distance between fft data elements is 2K */ \
/*	do 1 macros each processing 8 data values */ \
\
c0b:	start_timer(16); \
rsi = (uintptr_t)g->data_addr;/* Get source address */ \
/* rax = 0 */ \
/* rdx = 0 */ \
pfing(rbp = (uintptr_t)g->data_prefetch);/* Set prefetch distance */ \
pfing(rbp -= rsi); \
for(unsigned int loopA = 2*clm; loopA; loopA--) { \
	x1cl_eight_reals_unfft_1(rsi, 64); \
} \
\
/* Do 13 two_two_complex_unfft macros */ \
/*	distance between fft data elements is 2K */ \
/*	do 13 macros each processing 8 data values */ \
\
rsi += -2*clm*64+blkdst;	/* Next source pointer */ \
rdi = (uintptr_t)g->sincos1;/* Load sin/cos pointer */ \
/* 13 iterations of 2*clm */ \
for(unsigned int loopA = 13; loopA; loopA--) { \
	for(unsigned int loopB = 2*clm; loopB; loopB--) { \
		x1cl_two_complex_unfft(rsi, 64); \
	} \
	rsi += -2*clm*64+blkdst;	/* Next source pointer */ \
	rdi += XMM_SCD;		/* Next sine/cosine pointer */ \
} \
rsi += -14*blkdst;		/* Restore source pointer */ \
end_timer(16); \
\
/* Do inverse FFT levels 4,5 */ \
/* On input the 64-byte cache lines hold these data values: */ \
/*	0K	2K	4K	6K */ \
/*	8K	10K	12K	14K */ \
/*	... */ \
/* On output the 64-byte cache lines hold these data values: */ \
/*	0K	4K	8K	12K */ \
/*	2K	... */ \
/*	16K	24K	32K	40K */ \
/*	... */ \
\
/* Do 2 eight_reals_unfft_2 macros */ \
/*	distance between fft data elements is 4K */ \
/*	do 1 macros each processing 16 data values */ \
\
start_timer(18); \
for(unsigned int loopA = clm*2; loopA; ) { \
	do{ \
		x2cl_half_eight_reals_unfft_2(rsi, 64, blkdst); \
		L2prefetch128(u8ptr(rsi-64+rbp)); \
		pfing(rbp += -64+blkdst); \
	}while(--loopA & 1); \
	pfing(rbp += -2*blkdst+2*64);	/* Restore prefetching offset */ \
} \
\
/* Do 12 four_complex_unfft macros */ \
/*	distance between fft data elements is 4K */ \
/*	do 3 macros each processing 32 data values */ \
\
rsi += -2*clm*64+2*blkdst;	/* Next source pointer */ \
rdi = (uintptr_t)g->sincos2;/* Load sin/cos pointer */ \
/* 3 iterations of 2*clm */ \
for(unsigned int loopA = 3; loopA; loopA--) { \
	for(unsigned int loopB = clm*2; loopB; ) { \
		do{ \
			x4cl_four_complex_unfft(rsi, 64, blkdst, 2*blkdst, rdi); \
			L2prefetch128(u8ptr(rsi-64+rbp)); \
			L2prefetch128(u8ptr(rsi-64+blkdst+rbp)); \
			pfing(rbp += -64+2*blkdst); \
		}while(--loopB & 1); \
		pfing(rbp += -4*blkdst+2*64);	/* Restore prefetching offset */ \
	} \
	rsi += -2*clm*64+4*blkdst;	/* Next source pointer */ \
	rdi += XMM_SCD;		/* Next sine/cosine pointer */ \
} \
rsi += -3*4*blkdst-2*blkdst;/* Restore source pointer */ \
end_timer(18); \
\
/* Do inverse FFT levels 1,2,3 */ \
/* On input the 64-byte cache lines hold these data values: */ \
/*	0K	4K	8K	12K */ \
/*	2K	... */ \
/*	16K	24K	32K	40K */ \
/*	... */ \
/* On output the 64-byte cache lines hold these data values: */ \
/*	0K	4K	56K	60K */ \
/*	2K	... */ \
/*	8K	16K	64K	72K */ \
/*	10K */ \
/*	... */ \
/*	14K */ \
/*	24K	40K	80K	96K */ \
/*	26K */ \
/*	... */ \
/*	38K */ \
\
/* Do 16 seven_reals_last_unfft macros */ \
/*	distance between fft data elements is 16K */ \
/*	do 2 macros each processing 56 data values */ \
\
start_timer(19); \
/* 2 iterations of 2*clm */ \
xseven_reals_unfft_preload; \
for(unsigned int loopA = 2; loopA; loopA--) { \
	for(unsigned int loopB = clm*2; loopB; ) { \
		do{ \
			x7cl_seven_reals_last_unfft(rsi, 64, 2*blkdst); \
			L2prefetch128(u8ptr(rsi-64+rbp)); \
			L2prefetch128(u8ptr(rsi-64+2*blkdst+rbp)); \
			pfing(rbp += -64+4*blkdst); \
		}while(--loopB & 1); \
		pfing(rbp += -8*blkdst+2*64);	/* Restore prefetching offset */ \
	} \
	rsi += -2*clm*64+blkdst;	/* Next source pointer */ \
} \
end_timer(19); \
\
/* Normalize these values */ \
\
pass1_normalize(0, b0b); \
goto	endlp; \
\
/* Do FFT levels 1,2,3 */ \
/* */ \
/* On input the 64-byte cache lines hold these data values: */ \
/*	0K	4K	56K	60K */ \
/*	2K	... */ \
/*	8K	16K	64K	72K */ \
/*	10K */ \
/*	... */ \
/*	14K */ \
/*	24K	40K	80K	96K */ \
/*	26K */ \
/*	... */ \
/*	38K */ \
/* On output the 64-byte cache lines hold these data values: */ \
/*	0K	4K	8K	12K */ \
/*	2K	... */ \
/*	16K	24K	32K	40K */ \
/*	18K */ \
/*	... */ \
/*	22K */ \
/*	48K	64K	80K	96K */ \
/*	50K */ \
/*	... */ \
/*	62K */ \
\
/* Fire up auxiliary threads to do pass 1 forward FFT */ \
\
pass1:	pass1_forward_fft_setup(b0b); \
\
/* Do 16 seven_reals_first_fft macros (FFT levels 1,2,3) */ \
/*	distance between fft data elements is 16K */ \
/*	do 2 iterations each processing 56 data values */ \
\
b0b:	start_timer(21); \
rsi = (uintptr_t)g->data_addr;/* Load source address for prefetching */ \
rbx = g->DIST_TO_FFTSRCARG; \
pfing(rbp = (uintptr_t)g->data_prefetch);/* Set prefetch distance */ \
pfing(rbp -= rsi); \
/* 2 iterations of 2*clm */ \
/* rdx = 0 */ \
xseven_reals_fft_preload; \
for(unsigned int loopA = 2; loopA; loopA--) { \
	for(unsigned int loopB = clm*2; loopB; ) { \
		do{ \
			x7cl_seven_reals_first_fft(rsi, 64, 2*blkdst); \
			L2prefetch128(u8ptr(rsi-64+6*blkdst+rbp)); \
			L2prefetch128(u8ptr(rsi-64+8*blkdst+rbp)); \
			pfing(rbp += -64+4*blkdst); \
		}while(--loopB & 1); \
		pfing(rbp += -8*blkdst+2*64);	/* Restore prefetching offset */ \
	} \
	rsi += -2*clm*64+blkdst;	/* Next source pointer */ \
} \
rsi += -2*blkdst;		/* Restore source pointer */ \
end_timer(21); \
\
/* Do FFT levels 4,5 */ \
/* Values 0K-15K is real data, 16K-31K is semi-real data, 32K-111K is */ \
/* complex data. */ \
/* */ \
/* On input the 64-byte cache lines hold these data values: */ \
/*	0K	4K	8K	12K */ \
/*	2K	... */ \
/*	16K	24K	32K	40K */ \
/*	18K */ \
/*	... */ \
/*	22K */ \
/*	48K	64K	80K	96K */ \
/*	50K */ \
/*	... */ \
/*	62K */ \
/* On output the 64-byte cache lines hold these data values: */ \
/*	0K	4K	8K	12K */ \
/*	2K	... */ \
/*	16K	... */ \
/*	... */ \
\
/* Do 2 eight_reals_fft_2 macros */ \
/*	distance between fft data elements is 4K */ \
/*	do 2 macros each processing 8 data values */ \
\
start_timer(22); \
/* 2 iterations of 2*clm */ \
for(unsigned int loopA = 2; loopA; loopA--) { \
	for(unsigned int loopB = clm*2; loopB; ) { \
		do{ \
			x1cl_half_eight_reals_fft_2(rsi, 64); \
		}while(--loopB & 1); \
		L2prefetch128(u8ptr(rsi-128+rbp)); \
	} \
	rsi += -2*clm*64+blkdst;	/* Next source pointer */ \
} \
\
/* Do 4 four_complex_fft macros */ \
/*	distance between fft data elements is 4K */ \
/*	do 2 macros each processing 16 data values */ \
\
rdi = (uintptr_t)g->sincos2;/* Load sin/cos pointer */ \
/* 2 iterations of 2*clm */ \
for(unsigned int loopA = 2; loopA; loopA--) { \
	for(unsigned int loopB = clm*2; loopB; ) { \
		do{ \
			x2cl_four_complex_fft(rsi, 64, 2*blkdst); \
			L2prefetch128(u8ptr(rsi-64+rbp)); \
			pfing(rbp += -64+2*blkdst); \
		}while(--loopB & 1); \
		pfing(rbp += -4*blkdst+2*64);	/* Restore prefetching offset */ \
	} \
	rsi += -2*clm*64+blkdst;	/* Next source pointer */ \
} \
\
/* Do 8 four_complex_fft macros */ \
/*	distance between fft data elements is 4K */ \
/*	do 2 macros each processing 32 data values */ \
\
rsi += -2*blkdst+4*blkdst;	/* Next source pointer */ \
rdi += XMM_SCD;		/* Next sine/cosine pointer */ \
/* 2 iterations of 2*clm */ \
for(unsigned int loopA = 2; loopA; loopA--) { \
	for(unsigned int loopB = clm*2; loopB; ) { \
		do{ \
			x4cl_four_complex_fft(rsi, 64, 2*blkdst, 4*blkdst, rdi); \
			L2prefetch128(u8ptr(rsi-64+rbp)); \
			L2prefetch128(u8ptr(rsi-64+2*blkdst+rbp)); \
			pfing(rbp += -64+4*blkdst); \
		}while(--loopB & 1); \
		pfing(rbp += -8*blkdst+2*64);	/* Restore prefetching offset */ \
	} \
	rsi += -2*clm*64+blkdst;	/* Next source pointer */ \
} \
rsi += -2*blkdst-6*blkdst;	/* Restore source pointer */ \
end_timer(22); \
\
/* Do FFT levels 6 */ \
/* Values 0K-3K is real data, 4K-7K is semi-real data, 8K-111K is */ \
/* complex data. */ \
/* */ \
/* On input the 64-byte cache lines hold these data values: */ \
/*	0K	4K	8K	12K */ \
/*	2K	... */ \
/*	16K	... */ \
/*	... */ \
/* On output the 64-byte cache lines hold these data values: */ \
/*	0K	2K	4K	6K */ \
/*	8K	... */ \
/*	... */ \
\
/* Do 2 eight_reals_fft_1 macros */ \
/*	distance between fft data elements is 2K */ \
/*	do 1 macros each processing 16 data values */ \
\
start_timer(24); \
rdi = (uintptr_t)g->sincos1;/* Load sin/cos pointer */ \
for(unsigned int loopA = clm*2; loopA; ) { \
	do{ \
		x2cl_eight_reals_fft_1(rsi, 64, blkdst); \
	}while(--loopA & 1); \
} \
\
/* Do 12 two_two_complex_fft macros */ \
/*	distance between fft data elements is 2K */ \
/*	do 6 macros each processing 16 data values */ \
\
rsi += -2*clm*64+2*blkdst;	/* Next source pointer */ \
rdi += XMM_SCD;		/* Next sine/cosine pointer */ \
/* 6 iterations of 2*clm */ \
for(unsigned int loopA = 6; loopA; loopA--) { \
	for(unsigned int loopB = clm*2; loopB; ) { \
		do{ \
			x2cl_two_complex_fft(rsi, 64, blkdst, rdi); \
			L2prefetch128(u8ptr(rsi-64+rbp)); \
			pfing(rbp += -64+blkdst); \
		}while(--loopB & 1); \
		pfing(rbp += -2*blkdst+2*64);	/* Restore prefetching offset */ \
	} \
	rsi += -2*clm*64+2*blkdst;	/* Next source pointer */ \
	rdi += 2*XMM_SCD;		/* Next sine/cosine pointer */ \
} \
end_timer(24); \
\
/* Work on next 4 sets of 56 values */ \
\
endlp:	xpass1_get_next_block(pass2, c0b, b0b); /* Get next block and jump */ \
done:	xfft_footer_hg1; \
\
\
\
/* Perform a 128K-element FFT.  This is done in two passes.  Pass 1 */ \
/* does 6 levels, pass 2 uses common code to perform the last 11 levels. */ \
\

#define hg_pass1levels6pfa8(procname, xpass2name, pass2_levels, clm) \
/* (fftname, pass2_levels, clmarg) */ \
\
/* set_FFT_constants(fftname, 1, 0, pass2_levels, clmarg); */ \
xfft_header_hg1(procname, xpass2name, pass2_levels, clm); \
\
/* Do the last 11 FFT levels */ \
/*	do 1 pass2_11_levels_real and 15 pass2_11_levels_complex */ \
\
pass2:	start_timer(0); \
xpass2name(g); \
end_timer(0); \
if(g->ffttype == 1)		/* We're done if FFTing only */ \
	return; /* goto done; */ \
\
/* Do the inverse FFT */ \
\
pass1_inverse_fft_setup(c0b, procname); \
\
/* Do inverse FFT levels 6 */ \
/* On input and output the 64-byte cache lines hold these data values: */ \
/*	0K	2K	4K	6K */ \
/*	8K	10K	12K	14K */ \
/*	... */ \
\
/* Do 1 eight_reals_unfft_1 macros */ \
/*	distance between fft data elements is 2K */ \
/*	do 1 macros each processing 8 data values */ \
\
c0b:	start_timer(16); \
rsi = (uintptr_t)g->data_addr;/* Get source address */ \
/* rax = 0 */ \
/* rdx = 0 */ \
pfing(rbp = (uintptr_t)g->data_prefetch);/* Set prefetch distance */ \
pfing(rbp -= rsi); \
for(unsigned int loopA = 2*clm; loopA; loopA--) { \
	x1cl_eight_reals_unfft_1(rsi, 64); \
} \
\
/* Do 15 two_two_complex_unfft macros */ \
/*	distance between fft data elements is 2K */ \
/*	do 15 macros each processing 8 data values */ \
\
rsi += -2*clm*64+blkdst;	/* Next source pointer */ \
rdi = (uintptr_t)g->sincos1;/* Load sin/cos pointer */ \
/* 15 iterations of 2*clm */ \
for(unsigned int loopA = 15; loopA; loopA--) { \
	for(unsigned int loopB = 2*clm; loopB; loopB--) { \
		x1cl_two_complex_unfft(rsi, 64); \
	} \
	rsi += -2*clm*64+blkdst;	/* Next source pointer */ \
	rdi += XMM_SCD;		/* Next sine/cosine pointer */ \
} \
rsi += -16*blkdst;		/* Restore source pointer */ \
end_timer(16); \
\
/* Do inverse FFT levels 4,5 */ \
/* On input the 64-byte cache lines hold these data values: */ \
/*	0K	2K	4K	6K */ \
/*	8K	10K	12K	14K */ \
/*	... */ \
/* On output the 64-byte cache lines hold these data values: */ \
/*	0K	8K	16K	24K */ \
/*	2K	... */ \
/*	... */ \
/*	6K */ \
/*	32K	... */ \
/*	... */ \
\
/* Do 2 eight_reals_unfft_2 macros */ \
/* Do 2 nop_two_two_complex_unfft_2 macros */ \
/*	distance between fft data elements is 4K */ \
/*	do 1 macros each processing 32 data values */ \
\
start_timer(18); \
for(unsigned int loopA = clm*2; loopA; ) { \
	do{ \
		x4cl_eight_reals_unfft_2(rsi, 64, blkdst, 2*blkdst); \
		L2prefetch128(u8ptr(rsi-64+rbp)); \
		L2prefetch128(u8ptr(rsi-64+blkdst+rbp)); \
		pfing(rbp += -64+2*blkdst); \
	}while(--loopA & 1); \
	pfing(rbp += -4*blkdst+2*64);	/* Restore prefetching offset */ \
} \
\
/* Do 12 four_complex_unfft macros */ \
/*	distance between fft data elements is 4K */ \
/*	do 3 macros each processing 32 data values */ \
\
rsi += -2*clm*64+4*blkdst;	/* Next source pointer */ \
rdi = (uintptr_t)g->sincos2;/* Load sin/cos pointer */ \
/* 3 iterations of 2*clm */ \
for(unsigned int loopA = 3; loopA; loopA--) { \
	for(unsigned int loopB = clm*2; loopB; ) { \
		do{ \
			x4cl_four_complex_unfft(rsi, 64, blkdst, 2*blkdst, rdi); \
			L2prefetch128(u8ptr(rsi-64+rbp)); \
			L2prefetch128(u8ptr(rsi-64+blkdst+rbp)); \
			pfing(rbp += -64+2*blkdst); \
		}while(--loopB & 1); \
		pfing(rbp += -4*blkdst+2*64);	/* Restore prefetching offset */ \
	} \
	rsi += -2*clm*64+4*blkdst;	/* Next source pointer */ \
	rdi += XMM_SCD;		/* Next sine/cosine pointer */ \
} \
rsi += -4*4*blkdst;	/* Restore source pointer */ \
end_timer(18); \
\
/* Do inverse FFT levels 1,2,3 */ \
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
\
/* Do 16 eight_reals_last_unfft macros */ \
/*	distance between fft data elements is 16K */ \
/*	do 4 macros each processing 32 data values */ \
\
start_timer(19); \
/* 4 iterations of 2*clm */ \
for(unsigned int loopA = 4; loopA; loopA--) { \
	for(unsigned int loopB = clm*2; loopB; ) { \
		do{ \
			x4cl_eight_reals_last_unfft(rsi, 64, 4*blkdst, 8*blkdst); \
			L2prefetch128(u8ptr(rsi-64+rbp)); \
			L2prefetch128(u8ptr(rsi-64+4*blkdst+rbp)); \
			pfing(rbp += -64+8*blkdst); \
		}while(--loopB & 1); \
		pfing(rbp += -16*blkdst+2*64);	/* Restore prefetching offset */ \
	} \
	rsi += -2*clm*64+blkdst;	/* Next source pointer */ \
} \
end_timer(19); \
\
/* Normalize these values */ \
\
pass1_normalize(0, b0b); \
goto	endlp; \
\
/* Do FFT levels 1,2,3 */ \
/* */ \
/* On input the 64-byte cache lines hold these data values: */ \
/*	0K	32K	64K	96K */ \
/*	2K	... */ \
/*	... */ \
/*	30K */ \
/* On output the 64-byte cache lines hold these data values: */ \
/*	0K	16K	32K	48K */ \
/*	2K	... */ \
/*	... */ \
/*	14K	... */ \
/*	64K	... */ \
/*	... */ \
\
/* Fire up auxiliary threads to do pass 1 forward FFT */ \
\
pass1:	pass1_forward_fft_setup(b0b); \
\
/* Do 16 eight_reals_first_fft macros (FFT levels 1,2,3) */ \
/*	distance between fft data elements is 16K */ \
/*	do 8 iterations each processing 16 data values */ \
\
b0b:	start_timer(21); \
rsi = (uintptr_t)g->data_addr;/* Load source address for prefetching */ \
rbx = g->DIST_TO_FFTSRCARG; \
pfing(rbp = (uintptr_t)g->data_prefetch);/* Set prefetch distance */ \
pfing(rbp -= rsi); \
/* 8 iterations of 2*clm */ \
/* rdx = 0 */ \
for(unsigned int loopA = 8; loopA; loopA--) { \
	for(unsigned int loopB = clm*2; loopB; ) { \
		do{ \
			x2cl_eight_reals_first_fft(rsi, 64, 8*blkdst); \
			L2prefetch128(u8ptr(rsi-64+rbp)); \
			pfing(rbp += -64+8*blkdst); \
		}while(--loopB & 1); \
		pfing(rbp += -16*blkdst+2*64);	/* Restore prefetching offset */ \
	} \
	rsi += -2*clm*64+blkdst;	/* Next source pointer */ \
} \
rsi += -8*blkdst;		/* Restore source pointer */ \
end_timer(21); \
\
/* Do FFT levels 4,5 */ \
/* Values 0K-15K is real data, 16K-31K is semi-real data, 32K-127K is */ \
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
/*	0K	4K	8K	12K */ \
/*	2K	... */ \
/*	16K	... */ \
/*	... */ \
\
/* Do 2 eight_reals_fft_2 macros */ \
/* Do 2 nop_two_two_complex_fft_2 macros */ \
/* Do 4 four_complex_fft macros */ \
/*	distance between fft data elements is 4K */ \
/*	do 2 macros each processing 32 data values */ \
\
start_timer(22); \
rdi = (uintptr_t)g->sincos2;/* Load sin/cos pointer */ \
/* 2 iterations of 2*clm */ \
for(unsigned int loopA = 2; loopA; loopA--) { \
	for(unsigned int loopB = clm*2; loopB; ) { \
		do{ \
			x4cl_eight_reals_fft_2(rsi, 64, 2*blkdst, 4*blkdst, rdi); \
			L2prefetch128(u8ptr(rsi-64+rbp)); \
			L2prefetch128(u8ptr(rsi-64+2*blkdst+rbp)); \
			pfing(rbp += -64+4*blkdst); \
		}while(--loopB & 1); \
		pfing(rbp += -8*blkdst+2*64);	/* Restore prefetching offset */ \
	} \
	rsi += -2*clm*64+blkdst;	/* Next source pointer */ \
} \
\
/* Do 8 four_complex_fft macros */ \
/*	distance between fft data elements is 4K */ \
/*	do 2 macros each processing 32 data values */ \
\
rsi += -2*blkdst+8*blkdst;	/* Next source pointer */ \
rdi += XMM_SCD;		/* Next sine/cosine pointer */ \
/* 2 iterations of 2*clm */ \
for(unsigned int loopA = 2; loopA; loopA--) { \
	for(unsigned int loopB = clm*2; loopB; ) { \
		do{ \
			x4cl_four_complex_fft(rsi, 64, 2*blkdst, 4*blkdst, rdi); \
			L2prefetch128(u8ptr(rsi-64+rbp)); \
			L2prefetch128(u8ptr(rsi-64+2*blkdst+rbp)); \
			pfing(rbp += -64+4*blkdst); \
		}while(--loopB & 1); \
		pfing(rbp += -8*blkdst+2*64);	/* Restore prefetching offset */ \
	} \
	rsi += -2*clm*64+blkdst;	/* Next source pointer */ \
} \
rsi += -2*blkdst-8*blkdst;	/* Restore source pointer */ \
end_timer(22); \
\
/* Do FFT levels 6 */ \
/* Values 0K-3K is real data, 4K-7K is semi-real data, 8K-127K is */ \
/* complex data. */ \
/* */ \
/* On input the 64-byte cache lines hold these data values: */ \
/*	0K	4K	8K	12K */ \
/*	2K	... */ \
/*	16K	... */ \
/*	... */ \
/* On output the 64-byte cache lines hold these data values: */ \
/*	0K	2K	4K	6K */ \
/*	8K	... */ \
/*	... */ \
\
/* Do 2 eight_reals_fft_1 macros */ \
/*	distance between fft data elements is 2K */ \
/*	do 1 macros each processing 16 data values */ \
\
start_timer(24); \
rdi = (uintptr_t)g->sincos1;/* Load sin/cos pointer */ \
for(unsigned int loopA = 2*clm; loopA; loopA--) { \
	x2cl_eight_reals_fft_1(rsi, 64, blkdst); \
} \
\
/* Do 14 two_two_complex_fft macros */ \
/*	distance between fft data elements is 2K */ \
/*	do 7 macros each processing 16 data values */ \
\
rsi += -2*clm*64+2*blkdst;	/* Next source pointer */ \
rdi += XMM_SCD;		/* Next sine/cosine pointer */ \
/* 7 iterations of 2*clm */ \
for(unsigned int loopA = 7; loopA; loopA--) { \
	for(unsigned int loopB = clm*2; loopB; ) { \
		do{ \
			x2cl_two_complex_fft(rsi, 64, blkdst, rdi); \
			L2prefetch128(u8ptr(rsi-64+rbp)); \
			pfing(rbp += -64+blkdst); \
		}while(--loopB & 1); \
		pfing(rbp += -2*blkdst+2*64);	/* Restore prefetching offset */ \
	} \
	rsi += -2*clm*64+2*blkdst;	/* Next source pointer */ \
	rdi += 2*XMM_SCD;		/* Next sine/cosine pointer */ \
} \
end_timer(24); \
\
/* Work on next 4 sets of 64 values */ \
\
endlp:	xpass1_get_next_block(pass2, c0b, b0b); /* Get next block and jump */ \
done:	xfft_footer_hg1; \
\
\
\
/* Perform a 128K-element negacyclic FFT.  This is done in two passes.  Pass 1 */ \
/* does 6 levels, pass 2 uses common code to perform the last 11 levels. */ \
\

#define hg_pass1levels6complex4(procname, xpass2name, pass2_levels, clm) \
/* (fftname, pass2_levels, clmarg) */ \
\
/* set_FFT_constants(fftname, 1, 1, pass2_levels, clmarg); */ \
xfft_header_hg1(procname, xpass2name, pass2_levels, clm); \
\
/* Do the last 11 FFT levels */ \
/*	do 16 pass2_11_levels_complex */ \
\
pass2:	start_timer(0); \
xpass2name(g); \
end_timer(0); \
if(g->ffttype == 1)		/* We're done if FFTing only */ \
	return; /* goto done; */ \
\
/* Do the inverse FFT */ \
\
pass1_inverse_fft_setup(c0b, procname); \
\
/* Do inverse FFT levels 5,6 */ \
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
/* Do 16 four_complex_unfft macros */ \
/*	distance between fft data elements is 2K */ \
/*	do 8 macros each processing 16 data values */ \
\
c0b:	start_timer(16); \
rsi = (uintptr_t)g->data_addr;/* Get source address */ \
rdi = (uintptr_t)g->sincos1;/* Load sin/cos pointer */ \
/* 8 iterations of 2*clm */ \
/* rdx = 0 */ \
pfing(rbp = (uintptr_t)g->data_prefetch);/* Set prefetch distance */ \
pfing(rbp -= rsi); \
for(unsigned int loopA = 8; loopA; loopA--) { \
	for(unsigned int loopB = 2*clm; loopB; loopB--) { \
		x2cl_four_complex_unfft(rsi, 64, blkdst); \
	} \
	rsi += -2*clm*64+2*blkdst;	/* Next source pointer */ \
	rdi += XMM_SCD;		/* Next sine/cosine pointer */ \
} \
rsi += -8*2*blkdst;	/* Restore source pointer */ \
end_timer(16); \
\
/* Do inverse FFT levels 4 */ \
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
\
/* Do 16 two_two_complex_unfft macros */ \
/*	distance between fft data elements is 4K */ \
/*	do 8 macros each processing 16 data values */ \
\
start_timer(18); \
rdi = (uintptr_t)g->sincos2;/* Load sin/cos pointer */ \
/* eax = 4*256+2; */		/* 4 iters of 2 iters of 2*clm */ \
for(unsigned int loopA = 4; loopA; loopA--) { \
	for(unsigned int loopB = 2; loopB; loopB--) { \
		for(unsigned int loopC = clm*2; loopC; ) { \
			do{ \
				x2cl_two_complex_unfft(rsi, 64, 2*blkdst); \
				L2prefetch128(u8ptr(rsi-64+rbp)); \
				pfing(rbp += -64+2*blkdst); \
			}while(--loopC & 1); \
			pfing(rbp += -4*blkdst+2*64);	/* Restore prefetching offset */ \
		} \
		rsi += -2*clm*64+blkdst;	/* Next source pointer */ \
	} \
	rsi += -2*blkdst+4*blkdst;	/* Next source pointer */ \
	rdi += XMM_SCD;		/* Next sine/cosine pointer */ \
} \
rsi += -4*4*blkdst;	/* Restore source pointer */ \
end_timer(18); \
\
/* Do inverse FFT levels 1,2,3 */ \
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
\
/* Do 16 four_complex_last_unfft macros */ \
/*	distance between fft data elements is 16K */ \
/*	do 4 macros each processing 32 data values */ \
\
start_timer(19); \
rdi = (uintptr_t)g->plus1_premults;/* Address of premultiplier table */ \
/* 4 iterations of 2*clm */ \
for(unsigned int loopA = 4; loopA; loopA--) { \
	for(unsigned int loopB = clm*2; loopB; ) { \
		do{ \
			x4cl_four_complex_last_unfft(rsi, 64, 4*blkdst, 8*blkdst, 4*128); \
			L2prefetch128(u8ptr(rsi-64+rbp)); \
			L2prefetch128(u8ptr(rsi-64+4*blkdst+rbp)); \
			pfing(rbp += -64+8*blkdst); \
		}while(--loopB & 1); \
		pfing(rbp += -16*blkdst+2*64);	/* Restore prefetching offset */ \
	} \
	rsi += -2*clm*64+blkdst;	/* Next source pointer */ \
	rdi += 128;		/* Next premultiplier pointer */ \
} \
end_timer(19); \
\
/* Normalize these values */ \
\
pass1_normalize(0, b0b); \
goto	endlp; \
\
/* Do FFT levels 1,2,3 */ \
/* */ \
/* On input the 64-byte cache lines hold these data values: */ \
/*	0K	32K	64K	96K */ \
/*	2K	... */ \
/*	... */ \
/*	30K */ \
/* On output the 64-byte cache lines hold these data values: */ \
/*	0K	16K	32K	48K */ \
/*	2K	... */ \
/*	... */ \
/*	14K	... */ \
/*	64K	... */ \
/*	... */ \
\
/* Fire up auxiliary threads to do pass 1 forward FFT */ \
\
pass1:	pass1_forward_fft_setup(b0b); \
\
/* Do 16 four_complex_first_fft macros (FFT levels 1,2,3) */ \
/*	distance between fft data elements is 16K */ \
/*	do 8 iterations each processing 16 data values */ \
\
b0b:	start_timer(21); \
rsi = (uintptr_t)g->data_addr;/* Load source address for prefetching */ \
rbx = g->DIST_TO_FFTSRCARG; \
pfing(rbp = (uintptr_t)g->data_prefetch);/* Set prefetch distance */ \
pfing(rbp -= rsi); \
rdi = (uintptr_t)g->plus1_premults;/* Address of premultiplier table */ \
/* 8 iterations of 2*clm */ \
/* rdx = 0 */ \
for(unsigned int loopA = 8; loopA; loopA--) { \
	for(unsigned int loopB = clm*2; loopB; ) { \
		do{ \
			x2cl_four_complex_first_fft(rsi, 64, 8*blkdst); \
			L2prefetch128(u8ptr(rsi-64+rbp)); \
			pfing(rbp += -64+8*blkdst); \
		}while(--loopB & 1); \
		pfing(rbp += -16*blkdst+2*64);	/* Restore prefetching offset */ \
	} \
	rsi += -2*clm*64+blkdst;	/* Next source pointer */ \
	rdi += 128;		/* Next premultiplier pointer */ \
} \
rsi += -8*blkdst;		/* Restore source pointer */ \
end_timer(21); \
\
/* Do FFT levels 4 */ \
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
/*	6K */ \
/*	32K	... */ \
/*	... */ \
\
/* Do 16 two_two_complex_fft macros */ \
/*	distance between fft data elements is 8K */ \
/*	do 8 macros each processing 16 data values */ \
\
start_timer(22); \
rdi = (uintptr_t)g->sincos2;/* Load sin/cos pointer */ \
/* eax = 2*256+4; */ /* 2 iters of 4 iters of 2*clm */ \
for(unsigned int loopA = 2; loopA; loopA--) { \
	for(unsigned int loopB = 4; loopB; loopB--) { \
		for(unsigned int loopC = clm*2; loopC; ) { \
			do{ \
				x2cl_two_complex_fft(rsi, 64, 4*blkdst, rdi); \
				L2prefetch128(u8ptr(rsi-64+rbp)); \
				pfing(rbp += -64+4*blkdst); \
			}while(--loopC & 1); \
			pfing(rbp += -8*blkdst+2*64);	/* Restore prefetching offset */ \
		} \
		rsi += -2*clm*64+blkdst;	/* Next source pointer */ \
	} \
	rsi += -4*blkdst+8*blkdst;	/* Next source pointer */ \
	rdi += 2*XMM_SCD;		/* Next sin/cos pointer */ \
} \
rsi += -2*8*blkdst;	/* Restore source pointer */ \
end_timer(22); \
\
/* Do FFT levels 5,6 */ \
/* */ \
/* On input the 64-byte cache lines hold these data values: */ \
/*	0K	8K	16K	24K */ \
/*	2K	... */ \
/*	... */ \
/*	6K */ \
/*	32K	... */ \
/*	... */ \
/* On output the 64-byte cache lines hold these data values: */ \
/*	0K	2K	4K	6K */ \
/*	8K	... */ \
/*	... */ \
\
/* Do 16 four_complex_fft macros */ \
/*	distance between fft data elements is 2K */ \
/*	do 4 macros each processing 32 data values */ \
\
start_timer(24); \
rdi = (uintptr_t)g->sincos1;/* Load sin/cos pointer */ \
/* 4 iterations of 2*clm */ \
for(unsigned int loopA = 4; loopA; loopA--) { \
	for(unsigned int loopB = clm*2; loopB; ) { \
		do{ \
			x4cl_four_complex_fft(rsi, 64, blkdst, 2*blkdst, rdi); \
			L2prefetch128(u8ptr(rsi-64+rbp)); \
			L2prefetch128(u8ptr(rsi-64+blkdst+rbp)); \
			pfing(rbp += -64+2*blkdst); \
		}while(--loopB & 1); \
		pfing(rbp += -4*blkdst+2*64);	/* Restore prefetching offset */ \
	} \
	rsi += -2*clm*64+4*blkdst;	/* Next source pointer */ \
	rdi += 2*XMM_SCD;		/* Next sine/cosine pointer */ \
} \
end_timer(24); \
\
/* Work on next 4 sets of 64 values */ \
\
endlp:	xpass1_get_next_block(pass2, c0b, b0b); /* Get next block and jump */ \
done:	xfft_footer_hg1; \
\
\
\
/* Perform a 160K-element FFT.  This is done in two passes.  Pass 1 */ \
/* does 7 levels, pass 2 uses common code to perform the last 11 levels. */ \
\

#define hg_pass1levels7pfa5(procname, xpass2name, pass2_levels, clm) \
/* (fftname, pass2_levels, clmarg) */ \
\
/* set_FFT_constants(fftname, 1, 0, pass2_levels, clmarg); */ \
xfft_header_hg1(procname, xpass2name, pass2_levels, clm); \
\
/* Do the last 11 FFT levels */ \
/*	do 1 pass2_11_levels_real and 19 pass2_11_levels_complex */ \
\
pass2:	start_timer(0); \
xpass2name(g); \
end_timer(0); \
if(g->ffttype == 1)		/* We're done if FFTing only */ \
	return; /* goto done; */ \
\
/* Do the inverse FFT */ \
\
pass1_inverse_fft_setup(c0b, procname); \
\
/* Do inverse FFT levels 6,7 */ \
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
/* rax = 0 */ \
/* rdx = 0 */ \
pfing(rbp = (uintptr_t)g->data_prefetch);/* Set prefetch distance */ \
pfing(rbp -= rsi); \
for(unsigned int loopA = 2*clm; loopA; loopA--) { \
	x2cl_eight_reals_unfft_2(rsi, 64, blkdst); \
} \
\
/* Do 18 four_complex_unfft macros */ \
/*	distance between fft data elements is 2K */ \
/*	do 9 macros each processing 16 data values */ \
\
rsi += -2*clm*64+2*blkdst;	/* Next source pointer */ \
rdi = (uintptr_t)g->sincos1;/* Load sin/cos pointer */ \
/* 9 iterations of 2*clm */ \
for(unsigned int loopA = 9; loopA; loopA--) { \
	for(unsigned int loopB = 2*clm; loopB; loopB--) { \
		x2cl_four_complex_unfft(rsi, 64, blkdst); \
	} \
	rsi += -2*clm*64+2*blkdst;	/* Next source pointer */ \
	rdi += XMM_SCD;		/* Next sine/cosine pointer */ \
} \
rsi += -10*2*blkdst;	/* Restore source pointer */ \
end_timer(16); \
\
/* Do inverse FFT levels 4,5 */ \
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
/*	32K	48K	64K	80K */ \
/*	... */ \
\
/* Do 4 eight_reals_unfft_2 macros */ \
/*	distance between fft data elements is 8K */ \
/*	do 2 macros each processing 16 data values */ \
\
start_timer(18); \
/* 2 iterations of 2*clm */ \
for(unsigned int loopA = 2; loopA; loopA--) { \
	for(unsigned int loopB = clm*2; loopB; ) { \
		do{ \
			x2cl_half_eight_reals_unfft_2(rsi, 64, 2*blkdst); \
			L2prefetch128(u8ptr(rsi-64+rbp)); \
			pfing(rbp += -64+2*blkdst); \
		}while(--loopB & 1); \
		pfing(rbp += -4*blkdst+2*64);	/* Restore prefetching offset */ \
	} \
	rsi += -2*clm*64+blkdst;	/* Next source pointer */ \
} \
\
/* Do 16 four_complex_unfft macros */ \
/*	distance between fft data elements is 8K */ \
/*	do 4 macros each processing 32 data values */ \
\
rsi += -2*blkdst+4*blkdst;	/* Next source pointer */ \
rdi = (uintptr_t)g->sincos2;/* Load sin/cos pointer */ \
/* eax = 2*256+2; */ /* 2 iterations of 2 iters of 2*clm */ \
for(unsigned int loopA = 2; loopA; loopA--) { \
	for(unsigned int loopB = 2; loopB; loopB--) { \
		for(unsigned int loopC = clm*2; loopC; ) { \
			do{ \
				x4cl_four_complex_unfft(rsi, 64, 2*blkdst, 4*blkdst, rdi); \
				L2prefetch128(u8ptr(rsi-64+rbp)); \
				L2prefetch128(u8ptr(rsi-64+2*blkdst+rbp)); \
				pfing(rbp += -64+4*blkdst); \
			}while(--loopC & 1); \
			pfing(rbp += -8*blkdst+2*64);	/* Restore prefetching offset */ \
		} \
		rsi += -2*clm*64+blkdst;	/* Next source pointer */ \
	} \
	rsi += -2*blkdst+8*blkdst;	/* Next source pointer */ \
	rdi += XMM_SCD;		/* Next sine/cosine pointer */ \
} \
rsi += -2*8*blkdst-4*blkdst;/* Restore source pointer */ \
end_timer(18); \
\
/* Do inverse FFT levels 1,2,3 */ \
/* On input the 64-byte cache lines hold these data values: */ \
/*	0K	8K	16K	24K */ \
/*	2K	... */ \
/*	... */ \
/*	6K */ \
/*	32K	48K	64K	80K */ \
/*	... */ \
/* On output the 64-byte cache lines hold these data values: */ \
/*	0K	8K	80K	88K */ \
/*	2K	... */ \
/*	... */ \
/*	6K */ \
/*	16K	48K	96K	128K */ \
/*	18K */ \
/*	... */ \
/*	46K */ \
\
/* Do 32 five_reals_last_unfft macros */ \
/*	distance between fft data elements is 32K */ \
/*	do 4 macros each processing 40 data values */ \
\
start_timer(19); \
/* 4 iterations of 2*clm */ \
xfive_reals_unfft_preload;	/* Preload consts for five real macros */ \
for(unsigned int loopA = 4; loopA; loopA--) { \
	for(unsigned int loopB = clm*2; loopB; ) { \
		do{ \
			x5cl_five_reals_last_unfft(rsi, 64, 4*blkdst); \
			L2prefetch128(u8ptr(rsi-64+rbp)); \
			L2prefetch128(u8ptr(rsi-64+4*blkdst+rbp)); \
			pfing(rbp += -64+8*blkdst); \
		}while(--loopB & 1); \
		pfing(rbp += -16*blkdst+2*64);	/* Restore prefetching offset */ \
	} \
	rsi += -2*clm*64+blkdst;	/* Next source pointer */ \
} \
end_timer(19); \
\
/* Normalize these values */ \
\
pass1_normalize(0, b0b); \
goto	endlp; \
\
/* Do FFT levels 1,2,3 */ \
/* */ \
/* On input the 64-byte cache lines hold these data values: */ \
/*	0K	8K	80K	88K */ \
/*	2K	... */ \
/*	... */ \
/*	6K */ \
/*	16K	48K	96K	128K */ \
/*	18K */ \
/*	... */ \
/*	46K */ \
/* On output the 64-byte cache lines hold these data values: */ \
/*	0K	8K	16K	24K */ \
/*	2K	... */ \
/*	... */ \
/*	6K */ \
/*	32K	64K	96K	128K */ \
/*	34K */ \
/*	... */ \
/*	62K */ \
\
/* Fire up auxiliary threads to do pass 1 forward FFT */ \
\
pass1:	pass1_forward_fft_setup(b0b); \
\
/* Do 32 five_reals_first_fft macros (FFT levels 1,2,3) */ \
/*	distance between fft data elements is 32K */ \
/*	do 4 iterations each processing 40 data values */ \
\
b0b:	start_timer(21); \
rsi = (uintptr_t)g->data_addr;/* Load source address for prefetching */ \
rbx = g->DIST_TO_FFTSRCARG; \
pfing(rbp = (uintptr_t)g->data_prefetch);/* Set prefetch distance */ \
pfing(rbp -= rsi); \
/* 4 iterations of 2*clm */ \
/* rdx = 0 */ \
xfive_reals_fft_preload;		/* Preload consts for five real macros */ \
for(unsigned int loopA = 4; loopA; loopA--) { \
	for(unsigned int loopB = clm*2; loopB; ) { \
		do{ \
			x5cl_five_reals_first_fft(rsi, 64, 4*blkdst); \
			L2prefetch128(u8ptr(rsi-64+4*blkdst+rbp)); \
			L2prefetch128(u8ptr(rsi-64+8*blkdst+rbp)); \
			pfing(rbp += -64+8*blkdst); \
		}while(--loopB & 1); \
		pfing(rbp += -16*blkdst+2*64);	/* Restore prefetching offset */ \
	} \
	rsi += -2*clm*64+blkdst;	/* Next source pointer */ \
} \
rsi += -4*blkdst;		/* Restore source pointer */ \
end_timer(21); \
\
/* Do FFT levels 4,5 */ \
/* Values 0K-31K is real data, 32K-159K is complex data. */ \
/* */ \
/* On input the 64-byte cache lines hold these data values: */ \
/*	0K	8K	16K	24K */ \
/*	2K	... */ \
/*	... */ \
/*	6K */ \
/*	32K	64K	96K	128K */ \
/*	34K */ \
/*	... */ \
/*	62K */ \
/* On output the 64-byte cache lines hold these data values: */ \
/*	0K	8K	16K	24K */ \
/*	2K	... */ \
/*	... */ \
/*	6K	... */ \
/*	32K	... */ \
/*	... */ \
\
/* Do 4 eight_reals_fft_2 macros */ \
/*	distance between fft data elements is 8K */ \
/*	do 4 macros each processing 8 data values */ \
\
start_timer(22); \
/* 4 iterations of 2*clm */ \
for(unsigned int loopA = 4; loopA; loopA--) { \
	for(unsigned int loopB = clm*2; loopB; ) { \
		do{ \
			x1cl_half_eight_reals_fft_2(rsi, 64); \
		}while(--loopB & 1); \
		L2prefetch128(u8ptr(rsi-128+rbp)); \
	} \
	rsi += -2*clm*64+blkdst;	/* Next source pointer */ \
} \
\
/* Do 16 four_complex_fft macros */ \
/*	distance between fft data elements is 8K */ \
/*	do 4 macros each processing 32 data values */ \
\
rdi = (uintptr_t)g->sincos2;/* Load sin/cos pointer */ \
/* 4 iterations of 2*clm */ \
for(unsigned int loopA = 4; loopA; loopA--) { \
	for(unsigned int loopB = clm*2; loopB; ) { \
		do{ \
			x4cl_four_complex_fft(rsi, 64, 4*blkdst, 8*blkdst, rdi); \
			L2prefetch128(u8ptr(rsi-64+rbp)); \
			L2prefetch128(u8ptr(rsi-64+4*blkdst+rbp)); \
			pfing(rbp += -64+8*blkdst); \
		}while(--loopB & 1); \
		pfing(rbp += -16*blkdst+2*64);	/* Restore prefetching offset */ \
	} \
	rsi += -2*clm*64+blkdst;	/* Next source pointer */ \
} \
rsi += -4*blkdst-4*blkdst;	/* Restore source pointer */ \
end_timer(22); \
\
/* Do FFT levels 6,7 */ \
/* Values 0K-7K is real data, 8K-15K is semi-real data, 16K-159K is */ \
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
	x4cl_eight_reals_fft_2(rsi, 64, blkdst, 2*blkdst, rdi); \
} \
\
/* Do 16 four_complex_fft macros */ \
/*	distance between fft data elements is 2K */ \
/*	do 4 macros each processing 32 data values */ \
\
rsi += -2*clm*64+4*blkdst;	/* Next source pointer */ \
rdi += XMM_SCD;		/* Next sine/cosine pointer */ \
/* 4 iterations of 2*clm */ \
for(unsigned int loopA = 4; loopA; loopA--) { \
	for(unsigned int loopB = clm*2; loopB; ) { \
		do{ \
			x4cl_four_complex_fft(rsi, 64, blkdst, 2*blkdst, rdi); \
			L2prefetch128(u8ptr(rsi-64+rbp)); \
			L2prefetch128(u8ptr(rsi-64+blkdst+rbp)); \
			pfing(rbp += -64+2*blkdst); \
		}while(--loopB & 1); \
		pfing(rbp += -4*blkdst+2*64);	/* Restore prefetching offset */ \
	} \
	rsi += -2*clm*64+4*blkdst;	/* Next source pointer */ \
	rdi += 2*XMM_SCD;		/* Next sine/cosine pointer */ \
} \
end_timer(24); \
\
/* Work on next 4 sets of 80 values */ \
\
endlp:	xpass1_get_next_block(pass2, c0b, b0b); /* Get next block and jump */ \
done:	xfft_footer_hg1; \
\
\
\
/* Perform a 192K-element FFT.  This is done in two passes.  Pass 1 */ \
/* does 7 levels, pass 2 uses common code to perform the last 11 levels. */ \
\

#define hg_pass1levels7pfa6(procname, xpass2name, pass2_levels, clm) \
/* (fftname, pass2_levels, clmarg) */ \
\
/* set_FFT_constants(fftname, 1, 0, pass2_levels, clmarg); */ \
xfft_header_hg1(procname, xpass2name, pass2_levels, clm); \
\
/* Do the last 11 FFT levels */ \
/*	do 1 pass2_11_levels_real and 23 pass2_11_levels_complex */ \
\
pass2:	start_timer(0); \
xpass2name(g); \
end_timer(0); \
if(g->ffttype == 1)		/* We're done if FFTing only */ \
	return; /* goto done; */ \
\
/* Do the inverse FFT */ \
\
pass1_inverse_fft_setup(c0b, procname); \
\
/* Do inverse FFT levels 6,7 */ \
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
/* rax = 0 */ \
/* rdx = 0 */ \
pfing(rbp = (uintptr_t)g->data_prefetch);/* Set prefetch distance */ \
pfing(rbp -= rsi); \
for(unsigned int loopA = 2*clm; loopA; loopA--) { \
	x2cl_eight_reals_unfft_2(rsi, 64, blkdst); \
} \
\
/* Do 22 four_complex_unfft macros */ \
/*	distance between fft data elements is 2K */ \
/*	do 11 macros each processing 16 data values */ \
\
rsi += -2*clm*64+2*blkdst;	/* Next source pointer */ \
rdi = (uintptr_t)g->sincos1;/* Load sin/cos pointer */ \
/* 11 iterations of 2*clm */ \
for(unsigned int loopA = 11; loopA; loopA--) { \
	for(unsigned int loopB = 2*clm; loopB; loopB--) { \
		x2cl_four_complex_unfft(rsi, 64, blkdst); \
	} \
	rsi += -2*clm*64+2*blkdst;	/* Next source pointer */ \
	rdi += XMM_SCD;		/* Next sine/cosine pointer */ \
} \
rsi += -12*2*blkdst;	/* Restore source pointer */ \
end_timer(16); \
\
/* Do inverse FFT levels 4,5 */ \
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
start_timer(18); \
/* 2 iterations of 2*clm */ \
for(unsigned int loopA = 2; loopA; loopA--) { \
	for(unsigned int loopB = clm*2; loopB; ) { \
		do{ \
			x4cl_eight_reals_unfft_2(rsi, 64, 2*blkdst, 4*blkdst); \
			L2prefetch128(u8ptr(rsi-64+rbp)); \
			L2prefetch128(u8ptr(rsi-64+2*blkdst+rbp)); \
			pfing(rbp += -64+4*blkdst); \
		}while(--loopB & 1); \
		pfing(rbp += -8*blkdst+2*64);	/* Restore prefetching offset */ \
	} \
	rsi += -2*clm*64+blkdst;	/* Next source pointer */ \
} \
\
/* Do 16 four_complex_unfft macros */ \
/*	distance between fft data elements is 8K */ \
/*	do 4 macros each processing 32 data values */ \
\
rsi += -2*blkdst+8*blkdst;	/* Next source pointer */ \
rdi = (uintptr_t)g->sincos2;/* Load sin/cos pointer */ \
/* eax = 2*256+2; */ /* 2 iterations of 2 iters of 2*clm */ \
for(unsigned int loopA = 2; loopA; loopA--) { \
	for(unsigned int loopB = 2; loopB; loopB--) { \
		for(unsigned int loopC = clm*2; loopC; ) { \
			do{ \
				x4cl_four_complex_unfft(rsi, 64, 2*blkdst, 4*blkdst, rdi); \
				L2prefetch128(u8ptr(rsi-64+rbp)); \
				L2prefetch128(u8ptr(rsi-64+2*blkdst+rbp)); \
				pfing(rbp += -64+4*blkdst); \
			}while(--loopC & 1); \
			pfing(rbp += -8*blkdst+2*64);	/* Restore prefetching offset */ \
		} \
		rsi += -2*clm*64+blkdst;	/* Next source pointer */ \
	} \
	rsi += -2*blkdst+8*blkdst;	/* Next source pointer */ \
	rdi += XMM_SCD;		/* Next sine/cosine pointer */ \
} \
rsi += -3*8*blkdst;	/* Restore source pointer */ \
end_timer(18); \
\
/* Do inverse FFT levels 1,2,3 */ \
/* On input the 64-byte cache lines hold these data values: */ \
/*	0K	16K	32K	48K */ \
/*	2K	... */ \
/*	... */ \
/*	14K */ \
/*	64K	... */ \
/*	... */ \
/* On output the 64-byte cache lines hold these data values: */ \
/*	0K	48K	96K	144K */ \
/*	... */ \
/*	46K */ \
\
/* Do 32 six_reals_last_unfft macros */ \
/*	distance between fft data elements is 32K */ \
/*	do 8 macros each processing 24 data values */ \
\
start_timer(19); \
/* 8 iterations of 2*clm */ \
xsix_reals_unfft_preload; \
for(unsigned int loopA = 8; loopA; loopA--) { \
	for(unsigned int loopB = clm*2; loopB; ) { \
		do{ \
			x3cl_six_reals_last_unfft(rsi, 64, 8*blkdst); \
			L2prefetch128(u8ptr(rsi-64+rbp)); \
			pfing(rbp += -64+8*blkdst); \
		}while(--loopB & 1); \
		L2prefetch128(u8ptr(rsi+rbp)); \
		pfing(rbp += -16*blkdst+2*64);	/* Restore prefetching offset */ \
	} \
	rsi += -2*clm*64+blkdst;	/* Next source pointer */ \
} \
end_timer(19); \
\
/* Normalize these values */ \
\
pass1_normalize(0, b0b); \
goto	endlp; \
\
/* Do FFT levels 1,2,3 */ \
/* */ \
/* On input the 64-byte cache lines hold these data values: */ \
/*	0K	48K	96K	144K */ \
/*	... */ \
/*	46K */ \
/* On output the 64-byte cache lines hold these data values: */ \
/*	0K	16K	32K	48K */ \
/*	2K	... */ \
/*	... */ \
/*	14K */ \
/*	64K	96K	128K	160K */ \
/*	66K */ \
/*	... */ \
/*	94K */ \
\
/* Fire up auxiliary threads to do pass 1 forward FFT */ \
\
pass1:	pass1_forward_fft_setup(b0b); \
\
/* Do 32 six_reals_first_fft macros (FFT levels 1,2,3) */ \
/*	distance between fft data elements is 32K */ \
/*	do 8 iterations each processing 24 data values */ \
\
b0b:	start_timer(21); \
rsi = (uintptr_t)g->data_addr;/* Load source address for prefetching */ \
rbx = g->DIST_TO_FFTSRCARG; \
pfing(rbp = (uintptr_t)g->data_prefetch);/* Set prefetch distance */ \
pfing(rbp -= rsi); \
/* 8 iterations of 2*clm */ \
/* rdx = 0 */ \
xsix_reals_fft_preload; \
for(unsigned int loopA = 8; loopA; loopA--) { \
	for(unsigned int loopB = clm*2; loopB; ) { \
		do{ \
			x3cl_six_reals_first_fft(rsi, 64, 8*blkdst); \
			L2prefetch128(u8ptr(rsi-64+rbp)); \
			pfing(rbp += -64+8*blkdst); \
		}while(--loopB & 1); \
		L2prefetch128(u8ptr(rsi+rbp)); \
		pfing(rbp += -16*blkdst+2*64);	/* Restore prefetching offset */ \
	} \
	rsi += -2*clm*64+blkdst;	/* Next source pointer */ \
} \
rsi += -8*blkdst;		/* Restore source pointer */ \
end_timer(21); \
\
/* Do FFT levels 4,5 */ \
/* Values 0K-31K is real data, 32K-63K is semi-real data, 64K-191K is */ \
/* complex data. */ \
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
/* On output the 64-byte cache lines hold these data values: */ \
/*	0K	8K	16K	24K */ \
/*	2K	... */ \
/*	... */ \
/*	6K	... */ \
/*	32K	... */ \
/*	... */ \
\
/* Do 4 eight_reals_fft_2 macros */ \
/* Do 4 nop_two_two_complex_fft_2 macros */ \
/*	distance between fft data elements is 8K */ \
/*	do 4 macros each processing 16 data values */ \
\
start_timer(22); \
/* 4 iterations of 2*clm */ \
for(unsigned int loopA = 4; loopA; loopA--) { \
	for(unsigned int loopB = clm*2; loopB; ) { \
		do{ \
			x2cl_eight_reals_fft_2(rsi, 64, 4*blkdst); \
			L2prefetch128(u8ptr(rsi-64+rbp)); \
			pfing(rbp += -64+4*blkdst); \
		}while(--loopB & 1); \
		pfing(rbp += -8*blkdst+2*64);	/* Restore prefetching offset */ \
	} \
	rsi += -2*clm*64+blkdst;	/* Next source pointer */ \
} \
\
/* Do 16 four_complex_fft macros */ \
/*	distance between fft data elements is 8K */ \
/*	do 4 macros each processing 32 data values */ \
\
rsi += -4*blkdst+8*blkdst;	/* Next source pointer */ \
rdi = (uintptr_t)g->sincos2;/* Load sin/cos pointer */ \
/* 4 iterations of 2*clm */ \
for(unsigned int loopA = 4; loopA; loopA--) { \
	for(unsigned int loopB = clm*2; loopB; ) { \
		do{ \
			x4cl_four_complex_fft(rsi, 64, 4*blkdst, 8*blkdst, rdi); \
			L2prefetch128(u8ptr(rsi-64+rbp)); \
			L2prefetch128(u8ptr(rsi-64+4*blkdst+rbp)); \
			pfing(rbp += -64+8*blkdst); \
		}while(--loopB & 1); \
		pfing(rbp += -16*blkdst+2*64);	/* Restore prefetching offset */ \
	} \
	rsi += -2*clm*64+blkdst;	/* Next source pointer */ \
} \
rsi += -4*blkdst-8*blkdst;	/* Restore source pointer */ \
end_timer(22); \
\
/* Do FFT levels 6,7 */ \
/* Values 0K-7K is real data, 8K-15K is semi-real data, 16K-191K is */ \
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
	x4cl_eight_reals_fft_2(rsi, 64, blkdst, 2*blkdst, rdi); \
} \
\
/* Do 20 four_complex_fft macros */ \
/*	distance between fft data elements is 2K */ \
/*	do 5 macros each processing 32 data values */ \
\
rsi += -2*clm*64+4*blkdst;	/* Next source pointer */ \
rdi += XMM_SCD;		/* Next sine/cosine pointer */ \
/* 5 iterations of 2*clm */ \
for(unsigned int loopA = 5; loopA; loopA--) { \
	for(unsigned int loopB = clm*2; loopB; ) { \
		do{ \
			x4cl_four_complex_fft(rsi, 64, blkdst, 2*blkdst, rdi); \
			L2prefetch128(u8ptr(rsi-64+rbp)); \
			L2prefetch128(u8ptr(rsi-64+blkdst+rbp)); \
			pfing(rbp += -64+2*blkdst); \
		}while(--loopB & 1); \
		pfing(rbp += -4*blkdst+2*64);	/* Restore prefetching offset */ \
	} \
	rsi += -2*clm*64+4*blkdst;	/* Next source pointer */ \
	rdi += 2*XMM_SCD;		/* Next sine/cosine pointer */ \
} \
end_timer(24); \
\
/* Work on next 4 sets of 96 values */ \
\
endlp:	xpass1_get_next_block(pass2, c0b, b0b); /* Get next block and jump */ \
done:	xfft_footer_hg1; \
\
\
\
/* Perform a 192K-element negacyclic FFT.  This is done in two passes.  Pass 1 */ \
/* does 7 levels, pass 2 uses common code to perform the last 11 levels. */ \
\

#define hg_pass1levels7complex3(procname, xpass2name, pass2_levels, clm) \
/* (fftname, pass2_levels, clmarg) */ \
\
/* set_FFT_constants(fftname, 1, 1, pass2_levels, clmarg); */ \
xfft_header_hg1(procname, xpass2name, pass2_levels, clm); \
\
/* Do the last 11 FFT levels */ \
/*	do 24 pass2_11_levels_complex */ \
\
pass2:	start_timer(0); \
xpass2name(g); \
end_timer(0); \
if(g->ffttype == 1)		/* We're done if FFTing only */ \
	return; /* goto done; */ \
\
/* Do the inverse FFT */ \
\
pass1_inverse_fft_setup(c0b, procname); \
\
/* Do inverse FFT levels 6,7 */ \
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
/* Do 24 four_complex_unfft macros */ \
/*	distance between fft data elements is 2K */ \
/*	do 12 macros each processing 16 data values */ \
\
c0b:	start_timer(16); \
rsi = (uintptr_t)g->data_addr;/* Get source address */ \
rdi = (uintptr_t)g->sincos1;/* Load sin/cos pointer */ \
/* 12 iterations of 2*clm */ \
/* rdx = 0 */ \
pfing(rbp = (uintptr_t)g->data_prefetch);/* Set prefetch distance */ \
pfing(rbp -= rsi); \
for(unsigned int loopA = 12; loopA; loopA--) { \
	for(unsigned int loopB = 2*clm; loopB; loopB--) { \
		x2cl_four_complex_unfft(rsi, 64, blkdst); \
	} \
	rsi += -2*clm*64+2*blkdst;	/* Next source pointer */ \
	rdi += XMM_SCD;		/* Next sine/cosine pointer */ \
} \
rsi += -12*2*blkdst;	/* Restore source pointer */ \
end_timer(16); \
\
/* Do inverse FFT levels 4,5 */ \
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
/* Do 24 four_complex_unfft macros */ \
/*	distance between fft data elements is 8K */ \
/*	do 6 macros each processing 32 data values */ \
\
start_timer(18); \
rdi = (uintptr_t)g->sincos2;/* Load sin/cos pointer */ \
/* eax = 3*256+2; */ /* 3 iterations of 2 iters of 2*clm */ \
for(unsigned int loopA = 3; loopA; loopA--) { \
	for(unsigned int loopB = 2; loopB; loopB--) { \
		for(unsigned int loopC = clm*2; loopC; ) { \
			do{ \
				x4cl_four_complex_unfft(rsi, 64, 2*blkdst, 4*blkdst, rdi); \
				L2prefetch128(u8ptr(rsi-64+rbp)); \
				L2prefetch128(u8ptr(rsi-64+2*blkdst+rbp)); \
				pfing(rbp += -64+4*blkdst); \
			}while(--loopC & 1); \
			pfing(rbp += -8*blkdst+2*64);	/* Restore prefetching offset */ \
		} \
		rsi += -2*clm*64+blkdst;	/* Next source pointer */ \
	} \
	rsi += -2*blkdst+8*blkdst;	/* Next source pointer */ \
	rdi += XMM_SCD;		/* Next sine/cosine pointer */ \
} \
rsi += -3*8*blkdst;	/* Restore source pointer */ \
end_timer(18); \
\
/* Do inverse FFT levels 1,2,3 */ \
/* On input the 64-byte cache lines hold these data values: */ \
/*	0K	16K	32K	48K */ \
/*	2K	... */ \
/*	... */ \
/*	14K */ \
/*	64K	... */ \
/*	... */ \
/* On output the 64-byte cache lines hold these data values: */ \
/*	0K	48K	96K	144K */ \
/*	... */ \
/*	46K */ \
\
/* Do 32 three_complex_last_unfft macros */ \
/*	distance between fft data elements is 32K */ \
/*	do 8 macros each processing 24 data values */ \
\
start_timer(19); \
rdi = (uintptr_t)g->plus1_premults;/* Address of premultiplier table */ \
/* 8 iterations of 2*clm */ \
for(unsigned int loopA = 8; loopA; loopA--) { \
	for(unsigned int loopB = clm*2; loopB; ) { \
		do{ \
			x3cl_three_complex_last_unfft(rsi, 64, 8*blkdst); \
			L2prefetch128(u8ptr(rsi-64+rbp)); \
			pfing(rbp += -64+8*blkdst); \
		}while(--loopB & 1); \
		L2prefetch128(u8ptr(rsi+rbp)); \
		pfing(rbp += -16*blkdst+2*64);	/* Restore prefetching offset */ \
	} \
	rsi += -2*clm*64+blkdst;	/* Next source pointer */ \
	rdi += 192;		/* Next premultiplier pointer */ \
} \
end_timer(19); \
\
/* Normalize these values */ \
\
pass1_normalize(0, b0b); \
goto	endlp; \
\
/* Do FFT levels 1,2,3 */ \
/* */ \
/* On input the 64-byte cache lines hold these data values: */ \
/*	0K	48K	96K	144K */ \
/*	... */ \
/*	46K */ \
/* On output the 64-byte cache lines hold these data values: */ \
/*	0K	16K	32K	48K */ \
/*	2K	... */ \
/*	... */ \
/*	14K */ \
/*	64K	96K	128K	160K */ \
/*	66K */ \
/*	... */ \
/*	94K */ \
\
/* Fire up auxiliary threads to do pass 1 forward FFT */ \
\
pass1:	pass1_forward_fft_setup(b0b); \
\
/* Do 32 three_complex_first_fft macros (FFT levels 1,2,3) */ \
/*	distance between fft data elements is 32K */ \
/*	do 8 iterations each processing 24 data values */ \
\
b0b:	start_timer(21); \
rsi = (uintptr_t)g->data_addr;/* Load source address for prefetching */ \
rbx = g->DIST_TO_FFTSRCARG; \
pfing(rbp = (uintptr_t)g->data_prefetch);/* Set prefetch distance */ \
pfing(rbp -= rsi); \
rdi = (uintptr_t)g->plus1_premults;/* Address of premultiplier table */ \
/* 8 iterations of 2*clm */ \
/* rdx = 0 */ \
for(unsigned int loopA = 8; loopA; loopA--) { \
	for(unsigned int loopB = clm*2; loopB; ) { \
		do{ \
			x3cl_three_complex_first_fft(rsi, 64, 8*blkdst); \
			L2prefetch128(u8ptr(rsi-64+rbp)); \
			pfing(rbp += -64+8*blkdst); \
		}while(--loopB & 1); \
		L2prefetch128(u8ptr(rsi+rbp)); \
		pfing(rbp += -16*blkdst+2*64);	/* Restore prefetching offset */ \
	} \
	rsi += -2*clm*64+blkdst;	/* Next source pointer */ \
	rdi += 192;		/* Next premultiplier pointer */ \
} \
rsi += -8*blkdst;		/* Restore source pointer */ \
end_timer(21); \
\
/* Do FFT levels 4,5 */ \
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
/* On output the 64-byte cache lines hold these data values: */ \
/*	0K	8K	16K	24K */ \
/*	2K	... */ \
/*	... */ \
/*	6K	... */ \
/*	32K	... */ \
/*	... */ \
\
/* Do 24 four_complex_fft macros */ \
/*	distance between fft data elements is 8K */ \
/*	do 4 macros each processing 16 data values */ \
/*	do 4 macros each processing 32 data values */ \
\
start_timer(22); \
rdi = (uintptr_t)g->sincos2;/* Load sin/cos pointer */ \
/* 4 iterations of 2*clm */ \
for(unsigned int loopA = 4; loopA; loopA--) { \
	for(unsigned int loopB = clm*2; loopB; ) { \
		do{ \
			x2cl_four_complex_fft(rsi, 64, 4*blkdst); \
			L2prefetch128(u8ptr(rsi-64+rbp)); \
			pfing(rbp += -64+4*blkdst); \
		}while(--loopB & 1); \
		pfing(rbp += -8*blkdst+2*64);	/* Restore prefetching offset */ \
	} \
	rsi += -2*clm*64+blkdst;	/* Next source pointer */ \
} \
rsi += -4*blkdst+8*blkdst;	/* Next source pointer */ \
rdi += XMM_SCD;		/* Next sin/cos pointer */ \
/* 4 iterations of 2*clm */ \
for(unsigned int loopA = 4; loopA; loopA--) { \
	for(unsigned int loopB = clm*2; loopB; ) { \
		do{ \
			x4cl_four_complex_fft(rsi, 64, 4*blkdst, 8*blkdst, rdi); \
			L2prefetch128(u8ptr(rsi-64+rbp)); \
			L2prefetch128(u8ptr(rsi-64+4*blkdst+rbp)); \
			pfing(rbp += -64+8*blkdst); \
		}while(--loopB & 1); \
		pfing(rbp += -16*blkdst+2*64);	/* Restore prefetching offset */ \
	} \
	rsi += -2*clm*64+blkdst;	/* Next source pointer */ \
} \
rsi += -4*blkdst-8*blkdst;	/* Restore source pointer */ \
end_timer(22); \
\
/* Do FFT levels 6,7 */ \
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
/* Do 24 four_complex_fft macros */ \
/*	distance between fft data elements is 2K */ \
/*	do 6 macros each processing 32 data values */ \
\
start_timer(24); \
rdi = (uintptr_t)g->sincos1;/* Load sin/cos pointer */ \
/* 6 iterations of 2*clm */ \
for(unsigned int loopA = 6; loopA; loopA--) { \
	for(unsigned int loopB = clm*2; loopB; ) { \
		do{ \
			x4cl_four_complex_fft(rsi, 64, blkdst, 2*blkdst, rdi); \
			L2prefetch128(u8ptr(rsi-64+rbp)); \
			L2prefetch128(u8ptr(rsi-64+blkdst+rbp)); \
			pfing(rbp += -64+2*blkdst); \
		}while(--loopB & 1); \
		pfing(rbp += -4*blkdst+2*64);	/* Restore prefetching offset */ \
	} \
	rsi += -2*clm*64+4*blkdst;	/* Next source pointer */ \
	rdi += 2*XMM_SCD;		/* Next sine/cosine pointer */ \
} \
end_timer(24); \
\
/* Work on next 4 sets of 96 values */ \
\
endlp:	xpass1_get_next_block(pass2, c0b, b0b); /* Get next block and jump */ \
done:	xfft_footer_hg1; \
\
\
\
/* Perform a 224K-element FFT.  This is done in two passes.  Pass 1 */ \
/* does 7 levels, pass 2 uses common code to perform the last 11 levels. */ \
\

#define hg_pass1levels7pfa7(procname, xpass2name, pass2_levels, clm) \
/* (fftname, pass2_levels, clmarg) */ \
\
/* set_FFT_constants(fftname, 1, 0, pass2_levels, clmarg); */ \
xfft_header_hg1(procname, xpass2name, pass2_levels, clm); \
\
/* Do the last 11 FFT levels */ \
/*	do 1 pass2_11_levels_real and 27 pass2_11_levels_complex */ \
\
pass2:	start_timer(0); \
xpass2name(g); \
end_timer(0); \
if(g->ffttype == 1)		/* We're done if FFTing only */ \
	return; /* goto done; */ \
\
/* Do the inverse FFT */ \
\
pass1_inverse_fft_setup(c0b, procname); \
\
/* Do inverse FFT levels 6,7 */ \
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
/* rax = 0 */ \
/* rdx = 0 */ \
pfing(rbp = (uintptr_t)g->data_prefetch);/* Set prefetch distance */ \
pfing(rbp -= rsi); \
for(unsigned int loopA = 2*clm; loopA; loopA--) { \
	x2cl_eight_reals_unfft_2(rsi, 64, blkdst); \
} \
\
/* Do 26 four_complex_unfft macros */ \
/*	distance between fft data elements is 2K */ \
/*	do 13 macros each processing 16 data values */ \
\
rsi += -2*clm*64+2*blkdst;	/* Next source pointer */ \
rdi = (uintptr_t)g->sincos1;/* Load sin/cos pointer */ \
/* 13 iterations of 2*clm */ \
for(unsigned int loopA = 13; loopA; loopA--) { \
	for(unsigned int loopB = 2*clm; loopB; loopB--) { \
		x2cl_four_complex_unfft(rsi, 64, blkdst); \
	} \
	rsi += -2*clm*64+2*blkdst;	/* Next source pointer */ \
	rdi += XMM_SCD;		/* Next sine/cosine pointer */ \
} \
rsi += -14*2*blkdst;	/* Restore source pointer */ \
end_timer(16); \
\
/* Do inverse FFT levels 4,5 */ \
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
/*	32K	48K	64K	80K */ \
/*	... */ \
\
/* Do 4 eight_reals_unfft_2 macros */ \
/*	distance between fft data elements is 8K */ \
/*	do 2 macros each processing 16 data values */ \
\
start_timer(18); \
/* 2 iterations of 2*clm */ \
for(unsigned int loopA = 2; loopA; loopA--) { \
	for(unsigned int loopB = clm*2; loopB; ) { \
		do{ \
			x2cl_half_eight_reals_unfft_2(rsi, 64, 2*blkdst); \
			L2prefetch128(u8ptr(rsi-64+rbp)); \
			pfing(rbp += -64+2*blkdst); \
		}while(--loopB & 1); \
		pfing(rbp += -4*blkdst+2*64);	/* Restore prefetching offset */ \
	} \
	rsi += -2*clm*64+blkdst;	/* Next source pointer */ \
} \
\
/* Do 24 four_complex_unfft macros */ \
/*	distance between fft data elements is 8K */ \
/*	do 6 macros each processing 32 data values */ \
\
rsi += -2*blkdst+4*blkdst;	/* Next source pointer */ \
rdi = (uintptr_t)g->sincos2;/* Load sin/cos pointer */ \
/* eax = 3*256+2; */ /* 3 iterations of 2 iters of 2*clm */ \
for(unsigned int loopA = 3; loopA; loopA--) { \
	for(unsigned int loopB = 2; loopB; loopB--) { \
		for(unsigned int loopC = clm*2; loopC; ) { \
			do{ \
				x4cl_four_complex_unfft(rsi, 64, 2*blkdst, 4*blkdst, rdi); \
				L2prefetch128(u8ptr(rsi-64+rbp)); \
				L2prefetch128(u8ptr(rsi-64+2*blkdst+rbp)); \
				pfing(rbp += -64+4*blkdst); \
			}while(--loopC & 1); \
			pfing(rbp += -8*blkdst+2*64);	/* Restore prefetching offset */ \
		} \
		rsi += -2*clm*64+blkdst;	/* Next source pointer */ \
	} \
	rsi += -2*blkdst+8*blkdst;	/* Next source pointer */ \
	rdi += XMM_SCD;		/* Next sine/cosine pointer */ \
} \
rsi += -3*8*blkdst-4*blkdst;/* Restore source pointer */ \
end_timer(18); \
\
/* Do inverse FFT levels 1,2,3 */ \
/* On input the 64-byte cache lines hold these data values: */ \
/*	0K	8K	16K	24K */ \
/*	2K	... */ \
/*	... */ \
/*	6K */ \
/*	32K	48K	64K	80K */ \
/*	34K */ \
/*	... */ \
/*	46K */ \
/*	96K	112K	128K	144K */ \
/*	98K */ \
/*	... */ \
/*	110K */ \
/*	160K	176K	192K	208K */ \
/*	162K */ \
/*	... */ \
/*	174K */ \
/* On output the 64-byte cache lines hold these data values: */ \
/*	0K	8K	112K	120K */ \
/*	2K	... */ \
/*	... */ \
/*	6K */ \
/*	16K	32K	128K	144K */ \
/*	18K */ \
/*	... */ \
/*	30K */ \
/*	48K	80K	160K	192K */ \
/*	50K */ \
/*	... */ \
/*	78K */ \
\
/* Do 32 seven_reals_last_unfft macros */ \
/*	distance between fft data elements is 32K */ \
/*	do 4 macros each processing 56 data values */ \
\
start_timer(19); \
/* 4 iterations of 2*clm */ \
xseven_reals_unfft_preload; \
for(unsigned int loopA = 4; loopA; loopA--) { \
	for(unsigned int loopB = clm*2; loopB; ) { \
		do{ \
			x7cl_seven_reals_last_unfft(rsi, 64, 4*blkdst); \
			L2prefetch128(u8ptr(rsi-64+rbp)); \
			L2prefetch128(u8ptr(rsi-64+4*blkdst+rbp)); \
			pfing(rbp += -64+8*blkdst); \
		}while(--loopB & 1); \
		pfing(rbp += -16*blkdst+2*64);	/* Restore prefetching offset */ \
	} \
	rsi += -2*clm*64+blkdst;	/* Next source pointer */ \
} \
end_timer(19); \
\
/* Normalize these values */ \
\
pass1_normalize(0, b0b); \
goto	endlp; \
\
/* Do FFT levels 1,2,3 */ \
/* */ \
/* On input the 64-byte cache lines hold these data values: */ \
/*	0K	8K	112K	120K */ \
/*	2K	... */ \
/*	... */ \
/*	6K */ \
/*	16K	32K	128K	144K */ \
/*	18K */ \
/*	... */ \
/*	30K */ \
/*	48K	80K	160K	192K */ \
/*	50K */ \
/*	... */ \
/*	78K */ \
/* On output the 64-byte cache lines hold these data values: */ \
/*	0K	8K	16K	24K */ \
/*	2K	... */ \
/*	... */ \
/*	6K */ \
/*	32K	48K	64K	80K */ \
/*	34K */ \
/*	... */ \
/*	46K */ \
/*	96K	128K	160K	192K */ \
/*	98K */ \
/*	... */ \
/*	126K */ \
\
/* Fire up auxiliary threads to do pass 1 forward FFT */ \
\
pass1:	pass1_forward_fft_setup(b0b); \
\
/* Do 32 seven_reals_first_fft macros (FFT levels 1,2,3) */ \
/*	distance between fft data elements is 32K */ \
/*	do 4 iterations each processing 56 data values */ \
\
b0b:	start_timer(21); \
rsi = (uintptr_t)g->data_addr;/* Load source address for prefetching */ \
rbx = g->DIST_TO_FFTSRCARG; \
pfing(rbp = (uintptr_t)g->data_prefetch);/* Set prefetch distance */ \
pfing(rbp -= rsi); \
/* 4 iterations of 2*clm */ \
/* rdx = 0 */ \
xseven_reals_fft_preload; \
for(unsigned int loopA = 4; loopA; loopA--) { \
	for(unsigned int loopB = clm*2; loopB; ) { \
		do{ \
			x7cl_seven_reals_first_fft(rsi, 64, 4*blkdst); \
			L2prefetch128(u8ptr(rsi-64+12*blkdst+rbp)); \
			L2prefetch128(u8ptr(rsi-64+16*blkdst+rbp)); \
			pfing(rbp += -64+8*blkdst); \
		}while(--loopB & 1); \
		pfing(rbp += -16*blkdst+2*64);	/* Restore prefetching offset */ \
	} \
	rsi += -2*clm*64+blkdst;	/* Next source pointer */ \
} \
rsi += -4*blkdst;		/* Restore source pointer */ \
end_timer(21); \
\
/* Do FFT levels 4,5 */ \
/* Values 0K-31K is real data, 32K-63K is semi-real data, 64K-223K is */ \
/* complex data. */ \
/* */ \
/* On input the 64-byte cache lines hold these data values: */ \
/*	0K	8K	16K	24K */ \
/*	2K	... */ \
/*	... */ \
/*	6K */ \
/*	32K	48K	64K	80K */ \
/*	34K */ \
/*	... */ \
/*	46K */ \
/*	96K	128K	160K	192K */ \
/*	98K */ \
/*	... */ \
/*	126K */ \
/* On output the 64-byte cache lines hold these data values: */ \
/*	0K	8K	16K	24K */ \
/*	2K	... */ \
/*	... */ \
/*	6K	... */ \
/*	32K	... */ \
/*	... */ \
\
/* Do 4 eight_reals_fft_2 macros */ \
/*	distance between fft data elements is 8K */ \
/*	do 4 macros each processing 8 data values */ \
\
start_timer(22); \
/* 4 iterations of 2*clm */ \
for(unsigned int loopA = 4; loopA; loopA--) { \
	for(unsigned int loopB = clm*2; loopB; ) { \
		do{ \
			x1cl_half_eight_reals_fft_2(rsi, 64); \
		}while(--loopB & 1); \
		L2prefetch128(u8ptr(rsi-128+rbp)); \
	} \
	rsi += -2*clm*64+blkdst;	/* Next source pointer */ \
} \
\
/* Do 8 four_complex_fft macros */ \
/*	distance between fft data elements is 8K */ \
/*	do 4 macros each processing 16 data values */ \
\
rdi = (uintptr_t)g->sincos2;/* Load sin/cos pointer */ \
/* 4 iterations of 2*clm */ \
for(unsigned int loopA = 4; loopA; loopA--) { \
	for(unsigned int loopB = clm*2; loopB; ) { \
		do{ \
			x2cl_four_complex_fft(rsi, 64, 4*blkdst); \
			L2prefetch128(u8ptr(rsi-64+rbp)); \
			pfing(rbp += -64+4*blkdst); \
		}while(--loopB & 1); \
		pfing(rbp += -8*blkdst+2*64);	/* Restore prefetching offset */ \
	} \
	rsi += -2*clm*64+blkdst;	/* Next source pointer */ \
} \
\
/* Do 16 four_complex_fft macros */ \
/*	distance between fft data elements is 8K */ \
/*	do 4 macros each processing 32 data values */ \
\
rsi += -4*blkdst+8*blkdst;	/* Next source pointer */ \
rdi += XMM_SCD;		/* Next sine/cosine pointer */ \
/* 4 iterations of 2*clm */ \
for(unsigned int loopA = 4; loopA; loopA--) { \
	for(unsigned int loopB = clm*2; loopB; ) { \
		do{ \
			x4cl_four_complex_fft(rsi, 64, 4*blkdst, 8*blkdst, rdi); \
			L2prefetch128(u8ptr(rsi-64+rbp)); \
			L2prefetch128(u8ptr(rsi-64+4*blkdst+rbp)); \
			pfing(rbp += -64+8*blkdst); \
		}while(--loopB & 1); \
		pfing(rbp += -16*blkdst+2*64);	/* Restore prefetching offset */ \
	} \
	rsi += -2*clm*64+blkdst;	/* Next source pointer */ \
} \
rsi += -4*blkdst-12*blkdst;/* Restore source pointer */ \
end_timer(22); \
\
/* Do FFT levels 6,7 */ \
/* Values 0K-7K is real data, 8K-15K is semi-real data, 16K-223K is */ \
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
	x4cl_eight_reals_fft_2(rsi, 64, blkdst, 2*blkdst, rdi); \
} \
\
/* Do 24 four_complex_fft macros */ \
/*	distance between fft data elements is 2K */ \
/*	do 6 macros each processing 32 data values */ \
\
rsi += -2*clm*64+4*blkdst;	/* Next source pointer */ \
rdi += XMM_SCD;		/* Next sine/cosine pointer */ \
/* 6 iterations of 2*clm */ \
for(unsigned int loopA = 6; loopA; loopA--) { \
	for(unsigned int loopB = clm*2; loopB; ) { \
		do{ \
			x4cl_four_complex_fft(rsi, 64, blkdst, 2*blkdst, rdi); \
			L2prefetch128(u8ptr(rsi-64+rbp)); \
			L2prefetch128(u8ptr(rsi-64+blkdst+rbp)); \
			pfing(rbp += -64+2*blkdst); \
		}while(--loopB & 1); \
		pfing(rbp += -4*blkdst+2*64);	/* Restore prefetching offset */ \
	} \
	rsi += -2*clm*64+4*blkdst;	/* Next source pointer */ \
	rdi += 2*XMM_SCD;		/* Next sine/cosine pointer */ \
} \
end_timer(24); \
\
/* Work on next 4 sets of 112 values */ \
\
endlp:	xpass1_get_next_block(pass2, c0b, b0b); /* Get next block and jump */ \
done:	xfft_footer_hg1; \
\
\
\
/* Perform a 256K-element FFT.  This is done in two passes.  Pass 1 */ \
/* does 7 levels, pass 2 uses common code to perform the last 11 levels. */ \
\

#define hg_pass1levels7pfa8(procname, xpass2name, pass2_levels, clm) \
/* (fftname, pass2_levels, clmarg) */ \
\
/* set_FFT_constants(fftname, 1, 0, pass2_levels, clmarg); */ \
xfft_header_hg1(procname, xpass2name, pass2_levels, clm); \
\
/* Do the last 11 FFT levels */ \
/*	do 1 pass2_11_levels_real and 31 pass2_11_levels_complex */ \
\
pass2:	start_timer(0); \
xpass2name(g); \
end_timer(0); \
if(g->ffttype == 1)		/* We're done if FFTing only */ \
	return; /* goto done; */ \
\
/* Do the inverse FFT */ \
\
pass1_inverse_fft_setup(c0b, procname); \
\
/* Do inverse FFT levels 6,7 */ \
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
/* rax = 0 */ \
/* rdx = 0 */ \
pfing(rbp = (uintptr_t)g->data_prefetch);/* Set prefetch distance */ \
pfing(rbp -= rsi); \
for(unsigned int loopA = 2*clm; loopA; loopA--) { \
	x2cl_eight_reals_unfft_2(rsi, 64, blkdst); \
} \
\
/* Do 30 four_complex_unfft macros */ \
/*	distance between fft data elements is 2K */ \
/*	do 15 macros each processing 16 data values */ \
\
rsi += -2*clm*64+2*blkdst;	/* Next source pointer */ \
rdi = (uintptr_t)g->sincos1;/* Load sin/cos pointer */ \
/* 15 iterations of 2*clm */ \
for(unsigned int loopA = 15; loopA; loopA--) { \
	for(unsigned int loopB = 2*clm; loopB; loopB--) { \
		x2cl_four_complex_unfft(rsi, 64, blkdst); \
	} \
	rsi += -2*clm*64+2*blkdst;	/* Next source pointer */ \
	rdi += XMM_SCD;		/* Next sine/cosine pointer */ \
} \
rsi += -16*2*blkdst;	/* Restore source pointer */ \
end_timer(16); \
\
/* Do inverse FFT levels 4,5 */ \
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
start_timer(18); \
/* 2 iterations of 2*clm */ \
for(unsigned int loopA = 2; loopA; loopA--) { \
	for(unsigned int loopB = clm*2; loopB; ) { \
		do{ \
			x4cl_eight_reals_unfft_2(rsi, 64, 2*blkdst, 4*blkdst); \
			L2prefetch128(u8ptr(rsi-64+rbp)); \
			L2prefetch128(u8ptr(rsi-64+2*blkdst+rbp)); \
			pfing(rbp += -64+4*blkdst); \
		}while(--loopB & 1); \
		pfing(rbp += -8*blkdst+2*64);	/* Restore prefetching offset */ \
	} \
	rsi += -2*clm*64+blkdst;	/* Next source pointer */ \
} \
\
/* Do 24 four_complex_unfft macros */ \
/*	distance between fft data elements is 8K */ \
/*	do 6 macros each processing 32 data values */ \
\
rsi += -2*blkdst+8*blkdst;	/* Next source pointer */ \
rdi = (uintptr_t)g->sincos2;/* Load sin/cos pointer */ \
/* eax = 3*256+2; */ /* 3 iterations of 2 iters of 2*clm */ \
for(unsigned int loopA = 3; loopA; loopA--) { \
	for(unsigned int loopB = 2; loopB; loopB--) { \
		for(unsigned int loopC = clm*2; loopC; ) { \
			do{ \
				x4cl_four_complex_unfft(rsi, 64, 2*blkdst, 4*blkdst, rdi); \
				L2prefetch128(u8ptr(rsi-64+rbp)); \
				L2prefetch128(u8ptr(rsi-64+2*blkdst+rbp)); \
				pfing(rbp += -64+4*blkdst); \
			}while(--loopC & 1); \
			pfing(rbp += -8*blkdst+2*64);	/* Restore prefetching offset */ \
		} \
		rsi += -2*clm*64+blkdst;	/* Next source pointer */ \
	} \
	rsi += -2*blkdst+8*blkdst;	/* Next source pointer */ \
	rdi += XMM_SCD;		/* Next sine/cosine pointer */ \
} \
rsi += -4*8*blkdst;	/* Restore source pointer */ \
end_timer(18); \
\
/* Do inverse FFT levels 1,2,3 */ \
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
\
/* Do 32 eight_reals_last_unfft macros */ \
/*	distance between fft data elements is 32K */ \
/*	do 8 macros each processing 32 data values */ \
\
start_timer(19); \
/* 8 iterations of 2*clm */ \
for(unsigned int loopA = 8; loopA; loopA--) { \
	for(unsigned int loopB = clm*2; loopB; ) { \
		do{ \
			x4cl_eight_reals_last_unfft(rsi, 64, 8*blkdst, 16*blkdst); \
			L2prefetch128(u8ptr(rsi-64+rbp)); \
			L2prefetch128(u8ptr(rsi-64+8*blkdst+rbp)); \
			pfing(rbp += -64+16*blkdst); \
		}while(--loopB & 1); \
		pfing(rbp += -32*blkdst+2*64);	/* Restore prefetching offset */ \
	} \
	rsi += -2*clm*64+blkdst;	/* Next source pointer */ \
} \
end_timer(19); \
\
/* Normalize these values */ \
\
pass1_normalize(0, b0b); \
goto	endlp; \
\
/* Do FFT levels 1,2,3 */ \
/* */ \
/* On input the 64-byte cache lines hold these data values: */ \
/*	0K	64K	128K	192K */ \
/*	2K	... */ \
/*	... */ \
/*	62K */ \
/* On output the 64-byte cache lines hold these data values: */ \
/*	0K	32K	64K	96K */ \
/*	2K	... */ \
/*	... */ \
/*	30K	... */ \
/*	128K	... */ \
/*	... */ \
\
/* Fire up auxiliary threads to do pass 1 forward FFT */ \
\
pass1:	pass1_forward_fft_setup(b0b); \
\
/* Do 32 eight_reals_first_fft macros (FFT levels 1,2,3) */ \
/*	distance between fft data elements is 32K */ \
/*	do 16 iterations each processing 16 data values */ \
\
b0b:	start_timer(21); \
rsi = (uintptr_t)g->data_addr;/* Load source address for prefetching */ \
rbx = g->DIST_TO_FFTSRCARG; \
pfing(rbp = (uintptr_t)g->data_prefetch);/* Set prefetch distance */ \
pfing(rbp -= rsi); \
/* 16 iterations of 2*clm */ \
/* rdx = 0 */ \
for(unsigned int loopA = 16; loopA; loopA--) { \
	for(unsigned int loopB = clm*2; loopB; ) { \
		do{ \
			x2cl_eight_reals_first_fft(rsi, 64, 16*blkdst); \
			L2prefetch128(u8ptr(rsi-64+rbp)); \
			pfing(rbp += -64+16*blkdst); \
		}while(--loopB & 1); \
		pfing(rbp += -32*blkdst+2*64);	/* Restore prefetching offset */ \
	} \
	rsi += -2*clm*64+blkdst;	/* Next source pointer */ \
} \
rsi += -16*blkdst;		/* Restore source pointer */ \
end_timer(21); \
\
/* Do FFT levels 4,5 */ \
/* Values 0K-31K is real data, 32K-63K is semi-real data, 64K-255K is */ \
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
/* Do 4 eight_reals_fft_2 macros */ \
/* Do 4 nop_two_two_complex_fft_2 macros */ \
/* Do 8 four_complex_fft macros */ \
/*	distance between fft data elements is 8K */ \
/*	do 4 macros each processing 32 data values */ \
\
start_timer(22); \
rdi = (uintptr_t)g->sincos2;/* Load sin/cos pointer */ \
/* 4 iterations of 2*clm */ \
for(unsigned int loopA = 4; loopA; loopA--) { \
	for(unsigned int loopB = clm*2; loopB; ) { \
		do{ \
			x4cl_eight_reals_fft_2(rsi, 64, 4*blkdst, 8*blkdst, rdi); \
			L2prefetch128(u8ptr(rsi-64+rbp)); \
			L2prefetch128(u8ptr(rsi-64+4*blkdst+rbp)); \
			pfing(rbp += -64+8*blkdst); \
		}while(--loopB & 1); \
		pfing(rbp += -16*blkdst+2*64);	/* Restore prefetching offset */ \
	} \
	rsi += -2*clm*64+blkdst;	/* Next source pointer */ \
} \
\
/* Do 16 four_complex_fft macros */ \
/*	distance between fft data elements is 8K */ \
/*	do 4 macros each processing 32 data values */ \
\
rsi += -4*blkdst+16*blkdst;/* Next source pointer */ \
rdi += XMM_SCD;		/* Next sine/cosine pointer */ \
/* 4 iterations of 2*clm */ \
for(unsigned int loopA = 4; loopA; loopA--) { \
	for(unsigned int loopB = clm*2; loopB; ) { \
		do{ \
			x4cl_four_complex_fft(rsi, 64, 4*blkdst, 8*blkdst, rdi); \
			L2prefetch128(u8ptr(rsi-64+rbp)); \
			L2prefetch128(u8ptr(rsi-64+4*blkdst+rbp)); \
			pfing(rbp += -64+8*blkdst); \
		}while(--loopB & 1); \
		pfing(rbp += -16*blkdst+2*64);	/* Restore prefetching offset */ \
	} \
	rsi += -2*clm*64+blkdst;	/* Next source pointer */ \
} \
rsi += -4*blkdst-16*blkdst;/* Restore source pointer */ \
end_timer(22); \
\
/* Do FFT levels 6,7 */ \
/* Values 0K-7K is real data, 8K-15K is semi-real data, 16K-255K is */ \
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
	x4cl_eight_reals_fft_2(rsi, 64, blkdst, 2*blkdst, rdi); \
} \
\
/* Do 28 four_complex_fft macros */ \
/*	distance between fft data elements is 2K */ \
/*	do 7 macros each processing 32 data values */ \
\
rsi += -2*clm*64+4*blkdst;	/* Next source pointer */ \
rdi += XMM_SCD;		/* Next sine/cosine pointer */ \
/* 7 iterations of 2*clm */ \
for(unsigned int loopA = 7; loopA; loopA--) { \
	for(unsigned int loopB = clm*2; loopB; ) { \
		do{ \
			x4cl_four_complex_fft(rsi, 64, blkdst, 2*blkdst, rdi); \
			L2prefetch128(u8ptr(rsi-64+rbp)); \
			L2prefetch128(u8ptr(rsi-64+blkdst+rbp)); \
			pfing(rbp += -64+2*blkdst); \
		}while(--loopB & 1); \
		pfing(rbp += -4*blkdst+2*64);	/* Restore prefetching offset */ \
	} \
	rsi += -2*clm*64+4*blkdst;	/* Next source pointer */ \
	rdi += 2*XMM_SCD;		/* Next sine/cosine pointer */ \
} \
end_timer(24); \
\
/* Work on next 4 sets of 128 values */ \
\
endlp:	xpass1_get_next_block(pass2, c0b, b0b); /* Get next block and jump */ \
done:	xfft_footer_hg1; \
\
\
\
/* Perform a 256K-element negacyclic FFT.  This is done in two passes.  Pass 1 */ \
/* does 7 levels, pass 2 uses common code to perform the last 11 levels. */ \
\

#define hg_pass1levels7complex4(procname, xpass2name, pass2_levels, clm) \
/* (fftname, pass2_levels, clmarg) */ \
\
/* set_FFT_constants(fftname, 1, 1, pass2_levels, clmarg); */ \
xfft_header_hg1(procname, xpass2name, pass2_levels, clm); \
\
/* Do the last 11 FFT levels */ \
/*	do 32 pass2_11_levels_complex */ \
\
pass2:	start_timer(0); \
xpass2name(g); \
end_timer(0); \
if(g->ffttype == 1)		/* We're done if FFTing only */ \
	return; /* goto done; */ \
\
/* Do the inverse FFT */ \
\
pass1_inverse_fft_setup(c0b, procname); \
\
/* Do inverse FFT levels 6,7 */ \
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
/* Do 32 four_complex_unfft macros */ \
/*	distance between fft data elements is 2K */ \
/*	do 16 macros each processing 16 data values */ \
\
c0b:	start_timer(16); \
rsi = (uintptr_t)g->data_addr;/* Get source address */ \
rdi = (uintptr_t)g->sincos1;/* Load sin/cos pointer */ \
/* eax = 16*256; */ /* 16 iterations of 2*clm */ \
/* rdx = 0 */ \
pfing(rbp = (uintptr_t)g->data_prefetch);/* Set prefetch distance */ \
pfing(rbp -= rsi); \
for(unsigned int loopA = 16; loopA; loopA--) { \
	for(unsigned int loopB = 2*clm; loopB; loopB--) { \
		x2cl_four_complex_unfft(rsi, 64, blkdst); \
	} \
	rsi += -2*clm*64+2*blkdst;	/* Next source pointer */ \
	rdi += XMM_SCD;		/* Next sine/cosine pointer */ \
} \
rsi += -16*2*blkdst;	/* Restore source pointer */ \
end_timer(16); \
\
/* Do inverse FFT levels 4,5 */ \
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
/* Do 32 four_complex_unfft macros */ \
/*	distance between fft data elements is 8K */ \
/*	do 8 macros each processing 32 data values */ \
\
start_timer(18); \
rdi = (uintptr_t)g->sincos2;/* Load sin/cos pointer */ \
/* eax = 4*256+2; */ /* 4 iterations of 2 iters of 2*clm */ \
for(unsigned int loopA = 4; loopA; loopA--) { \
	for(unsigned int loopB = 2; loopB; loopB--) { \
		for(unsigned int loopC = clm*2; loopC; ) { \
			do{ \
				x4cl_four_complex_unfft(rsi, 64, 2*blkdst, 4*blkdst, rdi); \
				L2prefetch128(u8ptr(rsi-64+rbp)); \
				L2prefetch128(u8ptr(rsi-64+2*blkdst+rbp)); \
				pfing(rbp += -64+4*blkdst); \
			}while(--loopC & 1); \
			pfing(rbp += -8*blkdst+2*64);	/* Restore prefetching offset */ \
		} \
		rsi += -2*clm*64+blkdst;	/* Next source pointer */ \
	} \
	rsi += -2*blkdst+8*blkdst;	/* Next source pointer */ \
	rdi += XMM_SCD;		/* Next sine/cosine pointer */ \
} \
rsi += -4*8*blkdst;	/* Restore source pointer */ \
end_timer(18); \
\
/* Do inverse FFT levels 1,2,3 */ \
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
\
/* Do 32 four_complex_last_unfft macros */ \
/*	distance between fft data elements is 32K */ \
/*	do 8 macros each processing 32 data values */ \
\
start_timer(19); \
rdi = (uintptr_t)g->plus1_premults;/* Address of premultiplier table */ \
/* 8 iterations of 2*clm */ \
for(unsigned int loopA = 8; loopA; loopA--) { \
	for(unsigned int loopB = clm*2; loopB; ) { \
		do{ \
			x4cl_four_complex_last_unfft(rsi, 64, 8*blkdst, 16*blkdst, 8*128); \
			L2prefetch128(u8ptr(rsi-64+rbp)); \
			L2prefetch128(u8ptr(rsi-64+8*blkdst+rbp)); \
			pfing(rbp += -64+16*blkdst); \
		}while(--loopB & 1); \
		pfing(rbp += -32*blkdst+2*64);	/* Restore prefetching offset */ \
	} \
	rsi += -2*clm*64+blkdst;	/* Next source pointer */ \
	rdi += 128;		/* Next premultiplier pointer */ \
} \
end_timer(19); \
\
/* Normalize these values */ \
\
pass1_normalize(0, b0b); \
goto	endlp; \
\
/* Do FFT levels 1,2,3 */ \
/* */ \
/* On input the 64-byte cache lines hold these data values: */ \
/*	0K	64K	128K	192K */ \
/*	2K	... */ \
/*	... */ \
/*	62K */ \
/* On output the 64-byte cache lines hold these data values: */ \
/*	0K	32K	64K	96K */ \
/*	2K	... */ \
/*	... */ \
/*	30K	... */ \
/*	128K	... */ \
/*	... */ \
\
/* Fire up auxiliary threads to do pass 1 forward FFT */ \
\
pass1:	pass1_forward_fft_setup(b0b); \
\
/* Do 32 four_complex_first_fft macros (FFT levels 1,2,3) */ \
/*	distance between fft data elements is 32K */ \
/*	do 16 iterations each processing 16 data values */ \
\
b0b:	start_timer(21); \
rsi = (uintptr_t)g->data_addr;/* Load source address for prefetching */ \
rbx = g->DIST_TO_FFTSRCARG; \
pfing(rbp = (uintptr_t)g->data_prefetch);/* Set prefetch distance */ \
pfing(rbp -= rsi); \
rdi = (uintptr_t)g->plus1_premults;/* Address of premultiplier table */ \
/* 16 iterations of 2*clm */ \
/* rdx = 0 */ \
for(unsigned int loopA = 16; loopA; loopA--) { \
	for(unsigned int loopB = clm*2; loopB; ) { \
		do{ \
			x2cl_four_complex_first_fft(rsi, 64, 16*blkdst); \
			L2prefetch128(u8ptr(rsi-64+rbp)); \
			pfing(rbp += -64+16*blkdst); \
		}while(--loopB & 1); \
		pfing(rbp += -32*blkdst+2*64);	/* Restore prefetching offset */ \
	} \
	rsi += -2*clm*64+blkdst;	/* Next source pointer */ \
	rdi += 128;		/* Next premultiplier pointer */ \
} \
rsi += -16*blkdst;		/* Restore source pointer */ \
end_timer(21); \
\
/* Do FFT levels 4,5 */ \
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
/* Do 32 four_complex_fft macros */ \
/*	distance between fft data elements is 8K */ \
/*	do 8 macros each processing 32 data values */ \
\
start_timer(22); \
rdi = (uintptr_t)g->sincos2;/* Load sin/cos pointer */ \
/* eax = 2*256+4; */ /* 2 iterations 4 iters of 2*clm */ \
for(unsigned int loopA = 4; loopA; loopA--) { \
	for(unsigned int loopB = 2; loopB; loopB--) { \
		for(unsigned int loopC = clm*2; loopC; ) { \
			do{ \
				x4cl_four_complex_fft(rsi, 64, 4*blkdst, 8*blkdst, rdi); \
				L2prefetch128(u8ptr(rsi-64+rbp)); \
				L2prefetch128(u8ptr(rsi-64+4*blkdst+rbp)); \
				pfing(rbp += -64+8*blkdst); \
			}while(--loopC & 1); \
			pfing(rbp += -16*blkdst+2*64);	/* Restore prefetching offset */ \
		} \
		rsi += -2*clm*64+blkdst;	/* Next source pointer */ \
	} \
	rsi += -4*blkdst+16*blkdst;/* Next source pointer */ \
	rdi += 2*XMM_SCD;		/* Next sine/cosine pointer */ \
} \
rsi += -2*16*blkdst;	/* Restore source pointer */ \
end_timer(22); \
\
/* Do FFT levels 6,7 */ \
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
/* Do 32 four_complex_fft macros */ \
/*	distance between fft data elements is 2K */ \
/*	do 8 macros each processing 32 data values */ \
\
start_timer(24); \
rdi = (uintptr_t)g->sincos1;/* Load sin/cos pointer */ \
/* 8 iterations of 2*clm */ \
for(unsigned int loopA = 8; loopA; loopA--) { \
	for(unsigned int loopB = clm*2; loopB; ) { \
		do{ \
			x4cl_four_complex_fft(rsi, 64, blkdst, 2*blkdst, rdi); \
			L2prefetch128(u8ptr(rsi-64+rbp)); \
			L2prefetch128(u8ptr(rsi-64+blkdst+rbp)); \
			pfing(rbp += -64+2*blkdst); \
		}while(--loopB & 1); \
		pfing(rbp += -4*blkdst+2*64);	/* Restore prefetching offset */ \
	} \
	rsi += -2*clm*64+4*blkdst;	/* Next source pointer */ \
	rdi += 2*XMM_SCD;		/* Next sine/cosine pointer */ \
} \
end_timer(24); \
\
/* Work on next 4 sets of 128 values */ \
\
endlp:	xpass1_get_next_block(pass2, c0b, b0b); /* Get next block and jump */ \
done:	xfft_footer_hg1; \
\
\
\
/* Perform a 320K-element FFT.  This is done in two passes.  Pass 1 */ \
/* does 8 levels, pass 2 uses common code to perform the last 11 levels. */ \
\

#define hg_pass1levels8pfa5(procname, xpass2name, pass2_levels, clm) \
/* (fftname, pass2_levels, clmarg) */ \
\
/* set_FFT_constants(fftname, 1, 0, pass2_levels, clmarg); */ \
xfft_header_hg1(procname, xpass2name, pass2_levels, clm); \
\
/* Do the last 11 FFT levels */ \
/*	do 1 pass2_11_levels_real and 39 pass2_11_levels_complex */ \
\
pass2:	start_timer(0); \
xpass2name(g); \
end_timer(0); \
if(g->ffttype == 1)		/* We're done if FFTing only */ \
	return; /* goto done; */ \
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
/* rax = 0 */ \
/* rdx = 0 */ \
pfing(rbp = (uintptr_t)g->data_prefetch);/* Set prefetch distance */ \
pfing(rbp -= rsi); \
for(unsigned int loopA = 2*clm; loopA; loopA--) { \
	x2cl_eight_reals_unfft_2(rsi, 64, blkdst); \
} \
\
/* Do 38 four_complex_unfft macros */ \
/*	distance between fft data elements is 2K */ \
/*	do 19 macros each processing 16 data values */ \
\
rsi += -2*clm*64+2*blkdst;	/* Next source pointer */ \
rdi = (uintptr_t)g->sincos1;/* Load sin/cos pointer */ \
/* 19 iterations of 2*clm */ \
for(unsigned int loopA = 19; loopA; loopA--) { \
	for(unsigned int loopB = 2*clm; loopB; loopB--) { \
		x2cl_four_complex_unfft(rsi, 64, blkdst); \
	} \
	rsi += -2*clm*64+2*blkdst;	/* Next source pointer */ \
	rdi += XMM_SCD;		/* Next sine/cosine pointer */ \
} \
rsi += -20*2*blkdst;	/* Restore source pointer */ \
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
\
/* Do 4 eight_reals_unfft_1 macros */ \
/*	distance between fft data elements is 8K */ \
/*	do 2 macros each processing 16 data values */ \
\
start_timer(17); \
/* 2 iterations of 2*clm */ \
for(unsigned int loopA = 2; loopA; loopA--) { \
	for(unsigned int loopB = clm*2; loopB; ) { \
		do{ \
			x2cl_eight_reals_unfft_1(rsi, 64, 2*blkdst); \
		}while(--loopB & 1); \
		L2prefetch128(u8ptr(rsi-128+rbp)); \
	} \
	rsi += -2*clm*64+blkdst;	/* Next source pointer */ \
} \
\
/* Do 36 two_two_complex_unfft macros */ \
/*	distance between fft data elements is 8K */ \
/*	do 18 macros each processing 16 data values */ \
\
rsi += -2*blkdst+4*blkdst;	/* Next source pointer */ \
rdi = (uintptr_t)g->sincos2;/* Load sin/cos pointer */ \
/* eax = 9*256+2; */ /* 9 iterations of 2 iters of 2*clm */ \
for(unsigned int loopA = 9; loopA; loopA--) { \
	for(unsigned int loopB = 2; loopB; loopB--) { \
		for(unsigned int loopC = clm*2; loopC; ) { \
			do{ \
				x2cl_two_complex_unfft(rsi, 64, 2*blkdst); \
			}while(--loopC & 1); \
			L2prefetch128(u8ptr(rsi-128+rbp)); \
		} \
		rsi += -2*clm*64+blkdst;	/* Next source pointer */ \
	} \
	rsi += -2*blkdst+4*blkdst;	/* Next source pointer */ \
	rdi += XMM_SCD;		/* Next sine/cosine pointer */ \
} \
rsi += -10*4*blkdst;	/* Restore source pointer */ \
end_timer(17); \
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
\
/* Do 8 eight_reals_unfft_2 macros */ \
/*	distance between fft data elements is 16K */ \
/*	do 4 macros each processing 16 data values */ \
\
start_timer(18); \
/* 4 iterations of 2*clm */ \
for(unsigned int loopA = 4; loopA; loopA--) { \
	for(unsigned int loopB = clm*2; loopB; ) { \
		do{ \
			x2cl_half_eight_reals_unfft_2(rsi, 64, 4*blkdst); \
			L2prefetch128(u8ptr(rsi-64+rbp)); \
			pfing(rbp += -64+4*blkdst); \
		}while(--loopB & 1); \
		pfing(rbp += -8*blkdst+2*64);	/* Restore prefetching offset */ \
	} \
	rsi += -2*clm*64+blkdst;	/* Next source pointer */ \
} \
\
/* Do 32 four_complex_unfft macros */ \
/*	distance between fft data elements is 16K */ \
/*	do 8 macros each processing 32 data values */ \
\
rsi += -4*blkdst+8*blkdst;	/* Next source pointer */ \
rdi = (uintptr_t)g->sincos3;/* Load sin/cos pointer */ \
/* eax = 2*256+4; */ /* 2 iterations of 4 iters of 2*clm */ \
for(unsigned int loopA = 2; loopA; loopA--) { \
	for(unsigned int loopB = 4; loopB; loopB--) { \
		for(unsigned int loopC = clm*2; loopC; ) { \
			do{ \
				x4cl_four_complex_unfft(rsi, 64, 4*blkdst, 8*blkdst, rdi); \
				L2prefetch128(u8ptr(rsi-64+rbp)); \
				L2prefetch128(u8ptr(rsi-64+4*blkdst+rbp)); \
				pfing(rbp += -64+8*blkdst); \
			}while(--loopC & 1); \
			pfing(rbp += -16*blkdst+2*64);	/* Restore prefetching offset */ \
		} \
		rsi += -2*clm*64+blkdst;	/* Next source pointer */ \
	} \
	rsi += -4*blkdst+16*blkdst;/* Next source pointer */ \
	rdi += XMM_SCD;		/* Next sine/cosine pointer */ \
} \
rsi += -2*16*blkdst-8*blkdst;/* Restore source pointer */ \
end_timer(18); \
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
\
/* Do 64 five_reals_last_unfft macros */ \
/*	distance between fft data elements is 64K */ \
/*	do 8 macros each processing 40 data values */ \
\
start_timer(19); \
/* 8 iterations of 2*clm */ \
xfive_reals_unfft_preload;	/* Preload consts for five real macros */ \
for(unsigned int loopA = 8; loopA; loopA--) { \
	for(unsigned int loopB = clm*2; loopB; ) { \
		do{ \
			x5cl_five_reals_last_unfft(rsi, 64, 8*blkdst); \
			L2prefetch128(u8ptr(rsi-64+rbp)); \
			L2prefetch128(u8ptr(rsi-64+8*blkdst+rbp)); \
			pfing(rbp += -64+16*blkdst); \
		}while(--loopB & 1); \
		pfing(rbp += -32*blkdst+2*64);	/* Restore prefetching offset */ \
	} \
	rsi += -2*clm*64+blkdst;	/* Next source pointer */ \
} \
end_timer(19); \
\
/* Normalize these values */ \
\
pass1_normalize(0, b0b); \
goto	endlp; \
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
\
/* Fire up auxiliary threads to do pass 1 forward FFT */ \
\
pass1:	pass1_forward_fft_setup(b0b); \
\
/* Do 64 five_reals_first_fft macros (FFT levels 1,2,3) */ \
/*	distance between fft data elements is 64K */ \
/*	do 8 iterations each processing 40 data values */ \
\
b0b:	start_timer(21); \
rsi = (uintptr_t)g->data_addr;/* Load source address for prefetching */ \
rbx = g->DIST_TO_FFTSRCARG; \
pfing(rbp = (uintptr_t)g->data_prefetch);/* Set prefetch distance */ \
pfing(rbp -= rsi); \
/* 8 iterations of 2*clm */ \
/* rdx = 0 */ \
xfive_reals_fft_preload;		/* Preload consts for five real macros */ \
for(unsigned int loopA = 8; loopA; loopA--) { \
	for(unsigned int loopB = clm*2; loopB; ) { \
		do{ \
			x5cl_five_reals_first_fft(rsi, 64, 8*blkdst); \
			L2prefetch128(u8ptr(rsi-64+8*blkdst+rbp)); \
			L2prefetch128(u8ptr(rsi-64+16*blkdst+rbp)); \
			pfing(rbp += -64+16*blkdst); \
		}while(--loopB & 1); \
		pfing(rbp += -32*blkdst+2*64);	/* Restore prefetching offset */ \
	} \
	rsi += -2*clm*64+blkdst;	/* Next source pointer */ \
} \
rsi += -8*blkdst;		/* Restore source pointer */ \
end_timer(21); \
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
\
/* Do 8 eight_reals_fft_2 macros */ \
/*	distance between fft data elements is 16K */ \
/*	do 8 macros each processing 8 data values */ \
\
start_timer(22); \
/* 8 iterations of 2*clm */ \
for(unsigned int loopA = 8; loopA; loopA--) { \
	for(unsigned int loopB = clm*2; loopB; ) { \
		do{ \
			x1cl_half_eight_reals_fft_2(rsi, 64); \
		}while(--loopB & 1); \
		L2prefetch128(u8ptr(rsi-128+rbp)); \
	} \
	rsi += -2*clm*64+blkdst;	/* Next source pointer */ \
} \
\
/* Do 32 four_complex_fft macros */ \
/*	distance between fft data elements is 16K */ \
/*	do 8 macros each processing 32 data values */ \
\
rdi = (uintptr_t)g->sincos3;/* Load sin/cos pointer */ \
/* 8 iterations of 2*clm */ \
for(unsigned int loopA = 8; loopA; loopA--) { \
	for(unsigned int loopB = clm*2; loopB; ) { \
		do{ \
			x4cl_four_complex_fft(rsi, 64, 8*blkdst, 16*blkdst, rdi); \
			L2prefetch128(u8ptr(rsi-64+rbp)); \
			L2prefetch128(u8ptr(rsi-64+8*blkdst+rbp)); \
			pfing(rbp += -64+16*blkdst); \
		}while(--loopB & 1); \
		pfing(rbp += -32*blkdst+2*64);	/* Restore prefetching offset */ \
	} \
	rsi += -2*clm*64+blkdst;	/* Next source pointer */ \
} \
rsi += -8*blkdst-8*blkdst;	/* Restore source pointer */ \
end_timer(22); \
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
\
/* Do 4 eight_reals_fft_1 macros */ \
/* Do 4 two_two_complex_fft macros */ \
/*	distance between fft data elements is 8K */ \
/*	do 4 macros each processing 16 data values */ \
\
start_timer(23); \
rdi = (uintptr_t)g->sincos2;/* Load sin/cos pointer */ \
/* 4 iterations of 2*clm */ \
for(unsigned int loopA = 4; loopA; loopA--) { \
	for(unsigned int loopB = clm*2; loopB; ) { \
		do{ \
			x2cl_eight_reals_fft_1(rsi, 64, 4*blkdst); \
		}while(--loopB & 1); \
		L2prefetch128(u8ptr(rsi-128+rbp)); \
	} \
	rsi += -2*clm*64+blkdst;	/* Next source pointer */ \
} \
\
/* Do 32 two_two_complex_fft macros */ \
/*	distance between fft data elements is 8K */ \
/*	do 16 macros each processing 16 data values */ \
\
rsi += -4*blkdst+8*blkdst;	/* Next source pointer */ \
rdi += XMM_SCD;		/* Next sine/cosine pointer */ \
/* eax = 4*256+4; */ /* 4 iterations of 4 iters of 2*clm */ \
for(unsigned int loopA = 4; loopA; loopA--) { \
	for(unsigned int loopB = 4; loopB; loopB--) { \
		for(unsigned int loopC = 2; loopC; loopC--) { \
			for(unsigned int loopD = 2; loopD; loopD--) { \
				x2cl_two_complex_fft(rsi, 64, 4*blkdst, rdi); \
			} \
			L2prefetch128(u8ptr(rsi-128+rbp)); \
		} \
		rsi += -2*clm*64+blkdst;	/* Next source pointer */ \
	} \
	rsi += -4*blkdst+8*blkdst;	/* Next source pointer */ \
	rdi += 2*XMM_SCD;		/* Next sine/cosine pointer */ \
} \
rsi += -5*8*blkdst;	/* Restore source pointer */ \
end_timer(23); \
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
	x4cl_eight_reals_fft_2(rsi, 64, blkdst, 2*blkdst, rdi); \
} \
\
/* Do 36 four_complex_fft macros */ \
/*	distance between fft data elements is 2K */ \
/*	do 9 macros each processing 32 data values */ \
\
rsi += -2*clm*64+4*blkdst;	/* Next source pointer */ \
rdi += XMM_SCD;		/* Next sine/cosine pointer */ \
/* 9 iterations of 2*clm */ \
for(unsigned int loopA = 9; loopA; loopA--) { \
	for(unsigned int loopB = clm*2; loopB; ) { \
		do{ \
			x4cl_four_complex_fft(rsi, 64, blkdst, 2*blkdst, rdi); \
			L2prefetch128(u8ptr(rsi-64+rbp)); \
			L2prefetch128(u8ptr(rsi-64+blkdst+rbp)); \
			pfing(rbp += -64+2*blkdst); \
		}while(--loopB & 1); \
		pfing(rbp += -4*blkdst+2*64);	/* Restore prefetching offset */ \
	} \
	rsi += -2*clm*64+4*blkdst;	/* Next source pointer */ \
	rdi += 2*XMM_SCD;		/* Next sine/cosine pointer */ \
} \
end_timer(24); \
\
/* Work on next 4 sets of 160 values */ \
\
endlp:	xpass1_get_next_block(pass2, c0b, b0b); /* Get next block and jump */ \
done:	xfft_footer_hg1; \
\
\
\
/* Perform a 384K-element FFT.  This is done in two passes.  Pass 1 */ \
/* does 8 levels, pass 2 uses common code to perform the last 11 levels. */ \
\

#define hg_pass1levels8pfa6(procname, xpass2name, pass2_levels, clm) \
/* (fftname, pass2_levels, clmarg) */ \
\
/* set_FFT_constants(fftname, 1, 0, pass2_levels, clmarg); */ \
xfft_header_hg1(procname, xpass2name, pass2_levels, clm); \
\
/* Do the last 11 FFT levels */ \
/*	do 1 pass2_11_levels_real and 47 pass2_11_levels_complex */ \
\
pass2:	start_timer(0); \
xpass2name(g); \
end_timer(0); \
if(g->ffttype == 1)		/* We're done if FFTing only */ \
	return; /* goto done; */ \
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
/* rax = 0 */ \
/* rdx = 0 */ \
pfing(rbp = (uintptr_t)g->data_prefetch);/* Set prefetch distance */ \
pfing(rbp -= rsi); \
for(unsigned int loopA = 2*clm; loopA; loopA--) { \
	x2cl_eight_reals_unfft_2(rsi, 64, blkdst); \
} \
\
/* Do 46 four_complex_unfft macros */ \
/*	distance between fft data elements is 2K */ \
/*	do 23 macros each processing 16 data values */ \
\
rsi += -2*clm*64+2*blkdst;	/* Next source pointer */ \
rdi = (uintptr_t)g->sincos1;/* Load sin/cos pointer */ \
/* 23 iterations of 2*clm */ \
for(unsigned int loopA = 23; loopA; loopA--) { \
	for(unsigned int loopB = 2*clm; loopB; loopB--) { \
		x2cl_four_complex_unfft(rsi, 64, blkdst); \
	} \
	rsi += -2*clm*64+2*blkdst;	/* Next source pointer */ \
	rdi += XMM_SCD;		/* Next sine/cosine pointer */ \
} \
rsi += -24*2*blkdst;	/* Restore source pointer */ \
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
\
/* Do 4 eight_reals_unfft_1 macros */ \
/*	distance between fft data elements is 8K */ \
/*	do 2 macros each processing 16 data values */ \
\
start_timer(17); \
/* 2 iterations of 2*clm */ \
for(unsigned int loopA = 2; loopA; loopA--) { \
	for(unsigned int loopB = clm*2; loopB; ) { \
		do{ \
			x2cl_eight_reals_unfft_1(rsi, 64, 2*blkdst); \
		}while(--loopB & 1); \
		L2prefetch128(u8ptr(rsi-128+rbp)); \
	} \
	rsi += -2*clm*64+blkdst;	/* Next source pointer */ \
} \
\
/* Do 44 two_two_complex_unfft macros */ \
/*	distance between fft data elements is 8K */ \
/*	do 22 macros each processing 16 data values */ \
\
rsi += -2*blkdst+4*blkdst;	/* Next source pointer */ \
rdi = (uintptr_t)g->sincos2;/* Load sin/cos pointer */ \
/* eax = 11*256+2; */ /* 11 iterations of 2 iters of 2*clm */ \
for(unsigned int loopA = 11; loopA; loopA--) { \
	for(unsigned int loopB = 2; loopB; loopB--) { \
		for(unsigned int loopC = clm*2; loopC; ) { \
			do{ \
				x2cl_two_complex_unfft(rsi, 64, 2*blkdst); \
			}while(--loopC & 1); \
			L2prefetch128(u8ptr(rsi-128+rbp)); \
		} \
		rsi += -2*clm*64+blkdst;	/* Next source pointer */ \
	} \
	rsi += -2*blkdst+4*blkdst;	/* Next source pointer */ \
	rdi += XMM_SCD;		/* Next sine/cosine pointer */ \
} \
rsi += -12*4*blkdst;	/* Restore source pointer */ \
end_timer(17); \
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
\
/* Do 8 eight_reals_unfft_2 macros */ \
/* Do 8 nop_two_two_complex_unfft_2 macros */ \
/*	distance between fft data elements is 16K */ \
/*	do 4 macros each processing 32 data values */ \
\
start_timer(18); \
/* 4 iterations of 2*clm */ \
for(unsigned int loopA = 4; loopA; loopA--) { \
	for(unsigned int loopB = clm*2; loopB; ) { \
		do{ \
			x4cl_eight_reals_unfft_2(rsi, 64, 4*blkdst, 8*blkdst); \
			L2prefetch128(u8ptr(rsi-64+rbp)); \
			L2prefetch128(u8ptr(rsi-64+4*blkdst+rbp)); \
			pfing(rbp += -64+8*blkdst); \
		}while(--loopB & 1); \
		pfing(rbp += -16*blkdst+2*64);	/* Restore prefetching offset */ \
	} \
	rsi += -2*clm*64+blkdst;	/* Next source pointer */ \
} \
\
/* Do 32 four_complex_unfft macros */ \
/*	distance between fft data elements is 16K */ \
/*	do 8 macros each processing 32 data values */ \
\
rsi += -4*blkdst+16*blkdst;/* Next source pointer */ \
rdi = (uintptr_t)g->sincos3;/* Load sin/cos pointer */ \
/* eax = 2*256+4; */ /* 2 iterations of 4 iters of 2*clm */ \
for(unsigned int loopA = 2; loopA; loopA--) { \
	for(unsigned int loopB = 4; loopB; loopB--) { \
		for(unsigned int loopC = clm*2; loopC; ) { \
			do{ \
				x4cl_four_complex_unfft(rsi, 64, 4*blkdst, 8*blkdst, rdi); \
				L2prefetch128(u8ptr(rsi-64+rbp)); \
				L2prefetch128(u8ptr(rsi-64+4*blkdst+rbp)); \
				pfing(rbp += -64+8*blkdst); \
			}while(--loopC & 1); \
			pfing(rbp += -16*blkdst+2*64);	/* Restore prefetching offset */ \
		} \
		rsi += -2*clm*64+blkdst;	/* Next source pointer */ \
	} \
	rsi += -4*blkdst+16*blkdst;/* Next source pointer */ \
	rdi += XMM_SCD;		/* Next sine/cosine pointer */ \
} \
rsi += -3*16*blkdst;	/* Restore source pointer */ \
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
/* Do 64 six_reals_last_unfft macros */ \
/*	distance between fft data elements is 64K */ \
/*	do 16 macros each processing 24 data values */ \
\
start_timer(19); \
/* 16 iterations of 2*clm */ \
xsix_reals_unfft_preload; \
for(unsigned int loopA = 16; loopA; loopA--) { \
	for(unsigned int loopB = clm*2; loopB; ) { \
		do{ \
			x3cl_six_reals_last_unfft(rsi, 64, 16*blkdst); \
			L2prefetch128(u8ptr(rsi-64+rbp)); \
			pfing(rbp += -64+16*blkdst); \
		}while(--loopB & 1); \
		L2prefetch128(u8ptr(rsi+rbp)); \
		pfing(rbp += -32*blkdst+2*64);	/* Restore prefetching offset */ \
	} \
	rsi += -2*clm*64+blkdst;	/* Next source pointer */ \
} \
end_timer(19); \
\
/* Normalize these values */ \
\
pass1_normalize(0, b0b); \
goto	endlp; \
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
/* Fire up auxiliary threads to do pass 1 forward FFT */ \
\
pass1:	pass1_forward_fft_setup(b0b); \
\
/* Do 64 six_reals_first_fft macros (FFT levels 1,2,3) */ \
/*	distance between fft data elements is 64K */ \
/*	do 16 iterations each processing 24 data values */ \
\
b0b:	start_timer(21); \
rsi = (uintptr_t)g->data_addr;/* Load source address for prefetching */ \
rbx = g->DIST_TO_FFTSRCARG; \
pfing(rbp = (uintptr_t)g->data_prefetch);/* Set prefetch distance */ \
pfing(rbp -= rsi); \
/* 16 iterations of 2*clm */ \
/* rdx = 0 */ \
xsix_reals_fft_preload; \
for(unsigned int loopA = 16; loopA; loopA--) { \
	for(unsigned int loopB = clm*2; loopB; ) { \
		do{ \
			x3cl_six_reals_first_fft(rsi, 64, 16*blkdst); \
			L2prefetch128(u8ptr(rsi-64+rbp)); \
			pfing(rbp += -64+16*blkdst); \
		}while(--loopB & 1); \
		L2prefetch128(u8ptr(rsi+rbp)); \
		pfing(rbp += -32*blkdst+2*64);	/* Restore prefetching offset */ \
	} \
	rsi += -2*clm*64+blkdst;	/* Next source pointer */ \
} \
rsi += -16*blkdst;		/* Restore source pointer */ \
end_timer(21); \
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
\
/* Do 8 eight_reals_fft_2 macros */ \
/* Do 8 nop_two_two_complex_fft_2 macros */ \
/*	distance between fft data elements is 16K */ \
/*	do 8 macros each processing 16 data values */ \
\
start_timer(22); \
/* 8 iterations of 2*clm */ \
for(unsigned int loopA = 8; loopA; loopA--) { \
	for(unsigned int loopB = clm*2; loopB; ) { \
		do{ \
			x2cl_eight_reals_fft_2(rsi, 64, 8*blkdst); \
			L2prefetch128(u8ptr(rsi-64+rbp)); \
			pfing(rbp += -64+8*blkdst); \
		}while(--loopB & 1); \
		pfing(rbp += -16*blkdst+2*64);	/* Restore prefetching offset */ \
	} \
	rsi += -2*clm*64+blkdst;	/* Next source pointer */ \
} \
\
/* Do 32 four_complex_fft macros */ \
/*	distance between fft data elements is 16K */ \
/*	do 8 macros each processing 32 data values */ \
\
rsi += -8*blkdst+16*blkdst;/* Next source pointer */ \
rdi = (uintptr_t)g->sincos3;/* Load sin/cos pointer */ \
/* 8 iterations of 2*clm */ \
for(unsigned int loopA = 8; loopA; loopA--) { \
	for(unsigned int loopB = clm*2; loopB; ) { \
		do{ \
			x4cl_four_complex_fft(rsi, 64, 8*blkdst, 16*blkdst, rdi); \
			L2prefetch128(u8ptr(rsi-64+rbp)); \
			L2prefetch128(u8ptr(rsi-64+8*blkdst+rbp)); \
			pfing(rbp += -64+16*blkdst); \
		}while(--loopB & 1); \
		pfing(rbp += -32*blkdst+2*64);	/* Restore prefetching offset */ \
	} \
	rsi += -2*clm*64+blkdst;	/* Next source pointer */ \
} \
rsi += -8*blkdst-16*blkdst;/* Restore source pointer */ \
end_timer(22); \
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
\
/* Do 4 eight_reals_fft_1 macros */ \
/* Do 4 two_two_complex_fft macros */ \
/*	distance between fft data elements is 8K */ \
/*	do 4 macros each processing 16 data values */ \
\
start_timer(23); \
rdi = (uintptr_t)g->sincos2;/* Load sin/cos pointer */ \
/* 4 iterations of 2*clm */ \
for(unsigned int loopA = 4; loopA; loopA--) { \
	for(unsigned int loopB = clm*2; loopB; ) { \
		do{ \
			x2cl_eight_reals_fft_1(rsi, 64, 4*blkdst); \
		}while(--loopB & 1); \
		L2prefetch128(u8ptr(rsi-128+rbp)); \
	} \
	rsi += -2*clm*64+blkdst;	/* Next source pointer */ \
} \
\
/* Do 40 two_two_complex_fft macros */ \
/*	distance between fft data elements is 8K */ \
/*	do 20 macros each processing 16 data values */ \
\
rsi += -4*blkdst+8*blkdst;	/* Next source pointer */ \
rdi += XMM_SCD;		/* Next sine/cosine pointer */ \
/* eax = 5*256+4; */ /* 5 iterations of 4 iters of 2*clm */ \
for(unsigned int loopA = 5; loopA; loopA--) { \
	for(unsigned int loopB = 4; loopB; loopB--) { \
		for(unsigned int loopC = clm*2; loopC; ) { \
			do{ \
				x2cl_two_complex_fft(rsi, 64, 4*blkdst, rdi); \
			}while(--loopC & 1); \
			L2prefetch128(u8ptr(rsi-128+rbp)); \
		} \
		rsi += -2*clm*64+blkdst;	/* Next source pointer */ \
	} \
	rsi += -4*blkdst+8*blkdst;	/* Next source pointer */ \
	rdi += 2*XMM_SCD;		/* Next sine/cosine pointer */ \
} \
rsi += -6*8*blkdst;	/* Restore source pointer */ \
end_timer(23); \
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
	x4cl_eight_reals_fft_2(rsi, 64, blkdst, 2*blkdst, rdi); \
} \
\
/* Do 44 four_complex_fft macros */ \
/*	distance between fft data elements is 2K */ \
/*	do 11 macros each processing 32 data values */ \
\
rsi += -2*clm*64+4*blkdst;	/* Next source pointer */ \
rdi += XMM_SCD;		/* Next sine/cosine pointer */ \
/* 11 iterations of 2*clm */ \
for(unsigned int loopA = 11; loopA; loopA--) { \
	for(unsigned int loopB = clm*2; loopB; ) { \
		do{ \
			x4cl_four_complex_fft(rsi, 64, blkdst, 2*blkdst, rdi); \
			L2prefetch128(u8ptr(rsi-64+rbp)); \
			L2prefetch128(u8ptr(rsi-64+blkdst+rbp)); \
			pfing(rbp += -64+2*blkdst); \
		}while(--loopB & 1); \
		pfing(rbp += -4*blkdst+2*64);	/* Restore prefetching offset */ \
	} \
	rsi += -2*clm*64+4*blkdst;	/* Next source pointer */ \
	rdi += 2*XMM_SCD;		/* Next sine/cosine pointer */ \
} \
end_timer(24); \
\
/* Work on next 4 sets of 192 values */ \
\
endlp:	xpass1_get_next_block(pass2, c0b, b0b); /* Get next block and jump */ \
done:	xfft_footer_hg1; \
\
\
\
/* Perform a 384K-element negacyclic FFT.  This is done in two passes.  Pass 1 */ \
/* does 8 levels, pass 2 uses common code to perform the last 11 levels. */ \
\

#define hg_pass1levels8complex3(procname, xpass2name, pass2_levels, clm) \
/* (fftname, pass2_levels, clmarg) */ \
\
/* set_FFT_constants(fftname, 1, 1, pass2_levels, clmarg); */ \
xfft_header_hg1(procname, xpass2name, pass2_levels, clm); \
\
/* Do the last 11 FFT levels */ \
/*	do 48 pass2_11_levels_complex */ \
\
pass2:	start_timer(0); \
xpass2name(g); \
end_timer(0); \
if(g->ffttype == 1)		/* We're done if FFTing only */ \
	return; /* goto done; */ \
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
rdi = (uintptr_t)g->sincos1;/* Load sin/cos pointer */ \
/* 24 iterations of 2*clm */ \
/* rdx = 0 */ \
pfing(rbp = (uintptr_t)g->data_prefetch);/* Set prefetch distance */ \
pfing(rbp -= rsi); \
for(unsigned int loopA = 24; loopA; loopA--) { \
	for(unsigned int loopB = 2*clm; loopB; loopB--) { \
		x2cl_four_complex_unfft(rsi, 64, blkdst); \
	} \
	rsi += -2*clm*64+2*blkdst;	/* Next source pointer */ \
	rdi += XMM_SCD;		/* Next sine/cosine pointer */ \
} \
rsi += -24*2*blkdst;	/* Restore source pointer */ \
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
/* eax = 6*256+2; */ /* 6 iterations of 2 iters of 2*clm */ \
for(unsigned int loopA = 6; loopA; loopA--) { \
	for(unsigned int loopB = 2; loopB; loopB--) { \
		for(unsigned int loopC = clm*2; loopC; ) { \
			do{ \
				x4cl_four_complex_unfft(rsi, 64, 2*blkdst, 4*blkdst, rdi); \
				L2prefetch128(u8ptr(rsi-64+rbp)); \
				L2prefetch128(u8ptr(rsi-64+2*blkdst+rbp)); \
				pfing(rbp += -64+4*blkdst); \
			}while(--loopC & 1); \
			pfing(rbp += -8*blkdst+2*64);	/* Restore prefetching offset */ \
		} \
		rsi += -2*clm*64+blkdst;	/* Next source pointer */ \
	} \
	rsi += -2*blkdst+8*blkdst;	/* Next source pointer */ \
	rdi += XMM_SCD;		/* Next sine/cosine pointer */ \
} \
rsi += -6*8*blkdst;	/* Restore source pointer */ \
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
/* eax = 3*256+8; */ /* 3 iterations of 8 iters of 2*clm */ \
for(unsigned int loopA = 3; loopA; loopA--) { \
	for(unsigned int loopB = 8; loopB; loopB--) { \
		for(unsigned int loopC = clm*2; loopC; ) { \
			do{ \
				x2cl_two_complex_unfft(rsi, 64, 8*blkdst); \
			}while(--loopC & 1); \
			L2prefetch128(u8ptr(rsi-128+rbp)); \
		} \
		rsi += -2*clm*64+blkdst;	/* Next source pointer */ \
	} \
	rsi += -8*blkdst+16*blkdst;/* Next source pointer */ \
	rdi += XMM_SCD;		/* Next sine/cosine pointer */ \
} \
rsi += -3*16*blkdst;	/* Restore source pointer */ \
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
/* 16 iterations of 2*clm */ \
for(unsigned int loopA = 16; loopA; loopA--) { \
	for(unsigned int loopB = clm*2; loopB; ) { \
		do{ \
			x3cl_three_complex_last_unfft(rsi, 64, 16*blkdst); \
			L2prefetch128(u8ptr(rsi-64+rbp)); \
			pfing(rbp += -64+16*blkdst); \
		}while(--loopB & 1); \
		L2prefetch128(u8ptr(rsi+rbp)); \
		pfing(rbp += -32*blkdst+2*64);	/* Restore prefetching offset */ \
	} \
	rsi += -2*clm*64+blkdst;	/* Next source pointer */ \
	rdi += 192;		/* Next premultiplier pointer */ \
} \
end_timer(19); \
\
/* Normalize these values */ \
\
pass1_normalize(0, b0b); \
goto	endlp; \
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
/* Fire up auxiliary threads to do pass 1 forward FFT */ \
\
pass1:	pass1_forward_fft_setup(b0b); \
\
/* Do 64 three_complex_first_fft macros (FFT levels 1,2,3) */ \
/*	distance between fft data elements is 64K */ \
/*	do 16 iterations each processing 24 data values */ \
\
b0b:	start_timer(21); \
rsi = (uintptr_t)g->data_addr;/* Load source address for prefetching */ \
rbx = g->DIST_TO_FFTSRCARG; \
pfing(rbp = (uintptr_t)g->data_prefetch);/* Set prefetch distance */ \
pfing(rbp -= rsi); \
rdi = (uintptr_t)g->plus1_premults;/* Address of premultiplier table */ \
/* 16 iterations of 2*clm */ \
/* rdx = 0 */ \
for(unsigned int loopA = 16; loopA; loopA--) { \
	for(unsigned int loopB = clm*2; loopB; ) { \
		do{ \
			x3cl_three_complex_first_fft(rsi, 64, 16*blkdst); \
			L2prefetch128(u8ptr(rsi-64+rbp)); \
			pfing(rbp += -64+16*blkdst); \
		}while(--loopB & 1); \
		L2prefetch128(u8ptr(rsi+rbp)); \
		pfing(rbp += -32*blkdst+2*64);	/* Restore prefetching offset */ \
	} \
	rsi += -2*clm*64+blkdst;	/* Next source pointer */ \
	rdi += 192;		/* Next premultiplier pointer */ \
} \
rsi += -16*blkdst;		/* Restore source pointer */ \
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
/* 8 iterations of 2*clm */ \
for(unsigned int loopA = 8; loopA; loopA--) { \
	for(unsigned int loopB = clm*2; loopB; ) { \
		do{ \
			x2cl_two_complex_fft_in_place(rsi, 64, 8*blkdst, rdi); \
		}while(--loopB & 1); \
		L2prefetch128(u8ptr(rsi-128+rbp)); \
	} \
	rsi += -2*clm*64+blkdst;	/* Next source pointer */ \
} \
rsi += -8*blkdst+16*blkdst;/* Next source pointer */ \
rdi += XMM_SCD;		/* Next sin/cos pointer */ \
/* 16 iterations of 2*clm */ \
for(unsigned int loopA = 16; loopA; loopA--) { \
	for(unsigned int loopB = clm*2; loopB; ) { \
		do{ \
			x2cl_two_complex_fft(rsi, 64, 16*blkdst, rdi); \
		}while(--loopB & 1); \
		L2prefetch128(u8ptr(rsi-128+rbp)); \
	} \
	rsi += -2*clm*64+blkdst;	/* Next source pointer */ \
} \
rsi += -16*blkdst-16*blkdst;/* Restore source pointer */ \
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
/* eax = 3*256+4; */ /* 3 iterations of 4 iters of 2*clm */ \
for(unsigned int loopA = 3; loopA; loopA--) { \
	for(unsigned int loopB = 4; loopB; loopB--) { \
		for(unsigned int loopC = clm*2; loopC; ) { \
			do{ \
				x4cl_four_complex_fft(rsi, 64, 4*blkdst, 8*blkdst, rdi); \
				L2prefetch128(u8ptr(rsi-64+rbp)); \
				L2prefetch128(u8ptr(rsi-64+4*blkdst+rbp)); \
				pfing(rbp += -64+8*blkdst); \
			}while(--loopC & 1); \
			pfing(rbp += -16*blkdst+2*64);	/* Restore prefetching offset */ \
		} \
		rsi += -2*clm*64+blkdst;	/* Next source pointer */ \
	} \
	rsi += -4*blkdst+16*blkdst;/* Next source pointer */ \
	rdi += 2*XMM_SCD;		/* Next sine/cosine pointer */ \
} \
rsi += -3*16*blkdst;	/* Restore source pointer */ \
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
/* 12 iterations of 2*clm */ \
for(unsigned int loopA = 12; loopA; loopA--) { \
	for(unsigned int loopB = clm*2; loopB; ) { \
		do{ \
			x4cl_four_complex_fft(rsi, 64, blkdst, 2*blkdst, rdi); \
			L2prefetch128(u8ptr(rsi-64+rbp)); \
			L2prefetch128(u8ptr(rsi-64+blkdst+rbp)); \
			pfing(rbp += -64+2*blkdst); \
		}while(--loopB & 1); \
		pfing(rbp += -4*blkdst+2*64);	/* Restore prefetching offset */ \
	} \
	rsi += -2*clm*64+4*blkdst;	/* Next source pointer */ \
	rdi += 2*XMM_SCD;		/* Next sine/cosine pointer */ \
} \
end_timer(24); \
\
/* Work on next 4 sets of 192 values */ \
\
endlp:	xpass1_get_next_block(pass2, c0b, b0b); /* Get next block and jump */ \
done:	xfft_footer_hg1; \
\
\
\
/* Perform a 448K-element FFT.  This is done in two passes.  Pass 1 */ \
/* does 8 levels, pass 2 uses common code to perform the last 11 levels. */ \
\

#define hg_pass1levels8pfa7(procname, xpass2name, pass2_levels, clm) \
/* (fftname, pass2_levels, clmarg) */ \
\
/* set_FFT_constants(fftname, 1, 0, pass2_levels, clmarg); */ \
xfft_header_hg1(procname, xpass2name, pass2_levels, clm); \
\
/* Do the last 11 FFT levels */ \
/*	do 1 pass2_11_levels_real and 55 pass2_11_levels_complex */ \
\
pass2:	start_timer(0); \
xpass2name(g); \
end_timer(0); \
if(g->ffttype == 1)		/* We're done if FFTing only */ \
	return; /* goto done; */ \
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
/* rax = 0 */ \
/* rdx = 0 */ \
pfing(rbp = (uintptr_t)g->data_prefetch);/* Set prefetch distance */ \
pfing(rbp -= rsi); \
for(unsigned int loopA = 2*clm; loopA; loopA--) { \
	x2cl_eight_reals_unfft_2(rsi, 64, blkdst); \
} \
\
/* Do 54 four_complex_unfft macros */ \
/*	distance between fft data elements is 2K */ \
/*	do 27 macros each processing 16 data values */ \
\
rsi += -2*clm*64+2*blkdst;	/* Next source pointer */ \
rdi = (uintptr_t)g->sincos1;/* Load sin/cos pointer */ \
/* eax = 27*256; */ /* 27 iterations of 2*clm */ \
for(unsigned int loopA = 27; loopA; loopA--) { \
	for(unsigned int loopB = 2*clm; loopB; loopB--) { \
		x2cl_four_complex_unfft(rsi, 64, blkdst); \
	} \
	rsi += -2*clm*64+2*blkdst;	/* Next source pointer */ \
	rdi += XMM_SCD;		/* Next sine/cosine pointer */ \
} \
rsi += -28*2*blkdst;	/* Restore source pointer */ \
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
\
/* Do 4 eight_reals_unfft_1 macros */ \
/*	distance between fft data elements is 8K */ \
/*	do 2 macros each processing 16 data values */ \
\
start_timer(17); \
/* 2 iterations of 2*clm */ \
for(unsigned int loopA = 2; loopA; loopA--) { \
	for(unsigned int loopB = clm*2; loopB; ) { \
		do{ \
			x2cl_eight_reals_unfft_1(rsi, 64, 2*blkdst); \
		}while(--loopB & 1); \
		L2prefetch128(u8ptr(rsi-128+rbp)); \
	} \
	rsi += -2*clm*64+blkdst;	/* Next source pointer */ \
} \
\
/* Do 52 two_two_complex_unfft macros */ \
/*	distance between fft data elements is 8K */ \
/*	do 26 macros each processing 16 data values */ \
\
rsi += -2*blkdst+4*blkdst;	/* Next source pointer */ \
rdi = (uintptr_t)g->sincos2;/* Load sin/cos pointer */ \
/* eax = 13*256+2; */ /* 13 iterations of 2 iters of 2*clm */ \
for(unsigned int loopA = 13; loopA; loopA--) { \
	for(unsigned int loopB = 2; loopB; loopB--) { \
		for(unsigned int loopC = clm*2; loopC; ) { \
			do{ \
				x2cl_two_complex_unfft(rsi, 64, 2*blkdst); \
			}while(--loopC & 1); \
			L2prefetch128(u8ptr(rsi-128+rbp)); \
		} \
		rsi += -2*clm*64+blkdst;	/* Next source pointer */ \
	} \
	rsi += -2*blkdst+4*blkdst;	/* Next source pointer */ \
	rdi += XMM_SCD;		/* Next sine/cosine pointer */ \
} \
rsi += -14*4*blkdst;	/* Restore source pointer */ \
end_timer(17); \
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
\
/* Do 8 eight_reals_unfft_2 macros */ \
/*	distance between fft data elements is 16K */ \
/*	do 4 macros each processing 16 data values */ \
\
start_timer(18); \
/* 4 iterations of 2*clm */ \
for(unsigned int loopA = 4; loopA; loopA--) { \
	for(unsigned int loopB = clm*2; loopB; ) { \
		do{ \
			x2cl_half_eight_reals_unfft_2(rsi, 64, 4*blkdst); \
			L2prefetch128(u8ptr(rsi-64+rbp)); \
			pfing(rbp += -64+4*blkdst); \
		}while(--loopB & 1); \
		pfing(rbp += -8*blkdst+2*64);	/* Restore prefetching offset */ \
	} \
	rsi += -2*clm*64+blkdst;	/* Next source pointer */ \
} \
\
/* Do 48 four_complex_unfft macros */ \
/*	distance between fft data elements is 16K */ \
/*	do 12 macros each processing 32 data values */ \
\
rsi += -4*blkdst+8*blkdst;	/* Next source pointer */ \
rdi = (uintptr_t)g->sincos3;/* Load sin/cos pointer */ \
/* eax = 3*256+4; */ /* 3 iterations of 4 iters of 2*clm */ \
for(unsigned int loopA = 3; loopA; loopA--) { \
	for(unsigned int loopB = 4; loopB; loopB--) { \
		for(unsigned int loopC = clm*2; loopC; ) { \
			do{ \
				x4cl_four_complex_unfft(rsi, 64, 4*blkdst, 8*blkdst, rdi); \
				L2prefetch128(u8ptr(rsi-64+rbp)); \
				L2prefetch128(u8ptr(rsi-64+4*blkdst+rbp)); \
				pfing(rbp += -64+8*blkdst); \
			}while(--loopC & 1); \
			pfing(rbp += -16*blkdst+2*64);	/* Restore prefetching offset */ \
		} \
		rsi += -2*clm*64+blkdst;	/* Next source pointer */ \
	} \
	rsi += -4*blkdst+16*blkdst;/* Next source pointer */ \
	rdi += XMM_SCD;		/* Next sine/cosine pointer */ \
} \
rsi += -3*16*blkdst-8*blkdst;/* Restore source pointer */ \
end_timer(18); \
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
\
/* Do 64 seven_reals_last_unfft macros */ \
/*	distance between fft data elements is 64K */ \
/*	do 8 macros each processing 56 data values */ \
\
start_timer(19); \
/* 8 iterations of 2*clm */ \
xseven_reals_unfft_preload; \
for(unsigned int loopA = 8; loopA; loopA--) { \
	for(unsigned int loopB = clm*2; loopB; ) { \
		do{ \
			x7cl_seven_reals_last_unfft(rsi, 64, 8*blkdst); \
			L2prefetch128(u8ptr(rsi-64+rbp)); \
			L2prefetch128(u8ptr(rsi-64+8*blkdst+rbp)); \
			pfing(rbp += -64+16*blkdst); \
		}while(--loopB & 1); \
		pfing(rbp += -32*blkdst+2*64);	/* Restore prefetching offset */ \
	} \
	rsi += -2*clm*64+blkdst;	/* Next source pointer */ \
} \
end_timer(19); \
\
/* Normalize these values */ \
\
pass1_normalize(0, b0b); \
goto	endlp; \
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
\
/* Fire up auxiliary threads to do pass 1 forward FFT */ \
\
pass1:	pass1_forward_fft_setup(b0b); \
\
/* Do 64 seven_reals_first_fft macros (FFT levels 1,2,3) */ \
/*	distance between fft data elements is 64K */ \
/*	do 8 iterations each processing 56 data values */ \
\
b0b:	start_timer(21); \
rsi = (uintptr_t)g->data_addr;/* Load source address for prefetching */ \
rbx = g->DIST_TO_FFTSRCARG; \
pfing(rbp = (uintptr_t)g->data_prefetch);/* Set prefetch distance */ \
pfing(rbp -= rsi); \
/* 8 iterations of 2*clm */ \
/* rdx = 0 */ \
xseven_reals_fft_preload; \
for(unsigned int loopA = 8; loopA; loopA--) { \
	for(unsigned int loopB = clm*2; loopB; ) { \
		do{ \
			x7cl_seven_reals_first_fft(rsi, 64, 8*blkdst); \
			L2prefetch128(u8ptr(rsi-64+24*blkdst+rbp)); \
			L2prefetch128(u8ptr(rsi-64+32*blkdst+rbp)); \
			pfing(rbp += -64+16*blkdst); \
		}while(--loopB & 1); \
		pfing(rbp += -32*blkdst+2*64);	/* Restore prefetching offset */ \
	} \
	rsi += -2*clm*64+blkdst;	/* Next source pointer */ \
} \
rsi += -8*blkdst;		/* Restore source pointer */ \
end_timer(21); \
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
\
/* Do 8 eight_reals_fft_2 macros */ \
/*	distance between fft data elements is 16K */ \
/*	do 8 macros each processing 8 data values */ \
\
start_timer(22); \
/* 8 iterations of 2*clm */ \
for(unsigned int loopA = 8; loopA; loopA--) { \
	for(unsigned int loopB = clm*2; loopB; ) { \
		do{ \
			x1cl_half_eight_reals_fft_2(rsi, 64); \
		}while(--loopB & 1); \
		L2prefetch128(u8ptr(rsi-128+rbp)); \
	} \
	rsi += -2*clm*64+blkdst;	/* Next source pointer */ \
} \
\
/* Do 16 four_complex_fft macros */ \
/*	distance between fft data elements is 16K */ \
/*	do 8 macros each processing 16 data values */ \
\
rdi = (uintptr_t)g->sincos3;/* Load sin/cos pointer */ \
/* 8 iterations of 2*clm */ \
for(unsigned int loopA = 8; loopA; loopA--) { \
	for(unsigned int loopB = clm*2; loopB; ) { \
		do{ \
			x2cl_four_complex_fft(rsi, 64, 8*blkdst); \
			L2prefetch128(u8ptr(rsi-64+rbp)); \
			pfing(rbp += -64+8*blkdst); \
		}while(--loopB & 1); \
		pfing(rbp += -16*blkdst+2*64);	/* Restore prefetching offset */ \
	} \
	rsi += -2*clm*64+blkdst;	/* Next source pointer */ \
} \
\
/* Do 32 four_complex_fft macros */ \
/*	distance between fft data elements is 16K */ \
/*	do 8 macros each processing 32 data values */ \
\
rsi += -8*blkdst+16*blkdst;/* Next source pointer */ \
rdi += XMM_SCD;		/* Next sine/cosine pointer */ \
/* 8 iterations of 2*clm */ \
for(unsigned int loopA = 8; loopA; loopA--) { \
	for(unsigned int loopB = clm*2; loopB; ) { \
		do{ \
			x4cl_four_complex_fft(rsi, 64, 8*blkdst, 16*blkdst, rdi); \
			L2prefetch128(u8ptr(rsi-64+rbp)); \
			L2prefetch128(u8ptr(rsi-64+8*blkdst+rbp)); \
			pfing(rbp += -64+16*blkdst); \
		}while(--loopB & 1); \
		pfing(rbp += -32*blkdst+2*64);	/* Restore prefetching offset */ \
	} \
	rsi += -2*clm*64+blkdst;	/* Next source pointer */ \
} \
rsi += -8*blkdst-24*blkdst;/* Restore source pointer */ \
end_timer(22); \
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
\
start_timer(23); \
rdi = (uintptr_t)g->sincos2;/* Load sin/cos pointer */ \
/* 4 iterations of 2*clm */ \
for(unsigned int loopA = 4; loopA; loopA--) { \
	for(unsigned int loopB = clm*2; loopB; ) { \
		do{ \
			x2cl_eight_reals_fft_1(rsi, 64, 4*blkdst); \
		}while(--loopB & 1); \
		L2prefetch128(u8ptr(rsi-128+rbp)); \
	} \
	rsi += -2*clm*64+blkdst;	/* Next source pointer */ \
} \
\
/* Do 48 two_two_complex_fft macros */ \
/*	distance between fft data elements is 8K */ \
/*	do 24 macros each processing 16 data values */ \
\
rsi += -4*blkdst+8*blkdst;	/* Next source pointer */ \
rdi += XMM_SCD;		/* Next sine/cosine pointer */ \
/* eax = 6*256+4; */ /* 6 iterations of 4 iters of 2*clm */ \
for(unsigned int loopA = 6; loopA; loopA--) { \
	for(unsigned int loopB = 4; loopB; loopB--) { \
		for(unsigned int loopC = clm*2; loopC; ) { \
			do{ \
				x2cl_two_complex_fft(rsi, 64, 4*blkdst, rdi); \
			}while(--loopC & 1); \
			L2prefetch128(u8ptr(rsi-128+rbp)); \
		} \
		rsi += -2*clm*64+blkdst;	/* Next source pointer */ \
	} \
	rsi += -4*blkdst+8*blkdst;	/* Next source pointer */ \
	rdi += 2*XMM_SCD;		/* Next sine/cosine pointer */ \
} \
rsi += -7*8*blkdst;	/* Restore source pointer */ \
end_timer(23); \
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
	x4cl_eight_reals_fft_2(rsi, 64, blkdst, 2*blkdst, rdi); \
} \
\
/* Do 52 four_complex_fft macros */ \
/*	distance between fft data elements is 2K */ \
/*	do 13 macros each processing 32 data values */ \
\
rsi += -2*clm*64+4*blkdst;	/* Next source pointer */ \
rdi += XMM_SCD;		/* Next sine/cosine pointer */ \
/* 13 iterations of 2*clm */ \
for(unsigned int loopA = 13; loopA; loopA--) { \
	for(unsigned int loopB = clm*2; loopB; ) { \
		do{ \
			x4cl_four_complex_fft(rsi, 64, blkdst, 2*blkdst, rdi); \
			L2prefetch128(u8ptr(rsi-64+rbp)); \
			L2prefetch128(u8ptr(rsi-64+blkdst+rbp)); \
			pfing(rbp += -64+2*blkdst); \
		}while(--loopB & 1); \
		pfing(rbp += -4*blkdst+2*64);	/* Restore prefetching offset */ \
	} \
	rsi += -2*clm*64+4*blkdst;	/* Next source pointer */ \
	rdi += 2*XMM_SCD;		/* Next sine/cosine pointer */ \
} \
end_timer(24); \
\
/* Work on next 4 sets of 224 values */ \
\
endlp:	xpass1_get_next_block(pass2, c0b, b0b); /* Get next block and jump */ \
done:	xfft_footer_hg1; \
\
\
\
/* Perform a 512K-element FFT.  This is done in two passes.  Pass 1 */ \
/* does 8 levels, pass 2 uses common code to perform the last 11 levels. */ \
\

#define hg_pass1levels8pfa8(procname, xpass2name, pass2_levels, clm) \
/* (fftname, pass2_levels, clmarg) */ \
\
/* set_FFT_constants(fftname, 1, 0, pass2_levels, clmarg); */ \
xfft_header_hg1(procname, xpass2name, pass2_levels, clm); \
\
/* Do the last 11 FFT levels */ \
/*	do 1 pass2_11_levels_real and 63 pass2_11_levels_complex */ \
\
pass2:	start_timer(0); \
xpass2name(g); \
end_timer(0); \
if(g->ffttype == 1)		/* We're done if FFTing only */ \
	return; /* goto done; */ \
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
/* rax = 0 */ \
/* rdx = 0 */ \
pfing(rbp = (uintptr_t)g->data_prefetch);/* Set prefetch distance */ \
pfing(rbp -= rsi); \
for(unsigned int loopA = 2*clm; loopA; loopA--) { \
	x2cl_eight_reals_unfft_2(rsi, 64, blkdst); \
} \
\
/* Do 62 four_complex_unfft macros */ \
/*	distance between fft data elements is 2K */ \
/*	do 31 macros each processing 16 data values */ \
\
rsi += -2*clm*64+2*blkdst;	/* Next source pointer */ \
rdi = (uintptr_t)g->sincos1;/* Load sin/cos pointer */ \
/* 31 iterations of 2*clm */ \
for(unsigned int loopA = 31; loopA; loopA--) { \
	for(unsigned int loopB = 2*clm; loopB; loopB--) { \
		x2cl_four_complex_unfft(rsi, 64, blkdst); \
	} \
	rsi += -2*clm*64+2*blkdst;	/* Next source pointer */ \
	rdi += XMM_SCD;		/* Next sine/cosine pointer */ \
} \
rsi += -32*2*blkdst;	/* Restore source pointer */ \
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
\
/* Do 4 eight_reals_unfft_1 macros */ \
/*	distance between fft data elements is 8K */ \
/*	do 2 macros each processing 16 data values */ \
\
start_timer(17); \
/* 2 iterations of 2*clm */ \
for(unsigned int loopA = 2; loopA; loopA--) { \
	for(unsigned int loopB = clm*2; loopB; ) { \
		do{ \
			x2cl_eight_reals_unfft_1(rsi, 64, 2*blkdst); \
		}while(--loopB & 1); \
		L2prefetch128(u8ptr(rsi-128+rbp)); \
	} \
	rsi += -2*clm*64+blkdst;	/* Next source pointer */ \
} \
\
/* Do 60 two_two_complex_unfft macros */ \
/*	distance between fft data elements is 8K */ \
/*	do 30 macros each processing 16 data values */ \
\
rsi += -2*blkdst+4*blkdst;	/* Next source pointer */ \
rdi = (uintptr_t)g->sincos2;/* Load sin/cos pointer */ \
/* eax = 15*256+2; */ /* 15 iterations of 2 iters of 2*clm */ \
for(unsigned int loopA = 15; loopA; loopA--) { \
	for(unsigned int loopB = 2; loopB; loopB--) { \
		for(unsigned int loopC = clm*2; loopC; ) { \
			do{ \
				x2cl_two_complex_unfft(rsi, 64, 2*blkdst); \
			}while(--loopC & 1); \
			L2prefetch128(u8ptr(rsi-128+rbp)); \
		} \
		rsi += -2*clm*64+blkdst;	/* Next source pointer */ \
	} \
	rsi += -2*blkdst+4*blkdst;	/* Next source pointer */ \
	rdi += XMM_SCD;		/* Next sine/cosine pointer */ \
} \
rsi += -16*4*blkdst;	/* Restore source pointer */ \
end_timer(17); \
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
\
/* Do 8 eight_reals_unfft_2 macros */ \
/* Do 8 nop_two_two_complex_unfft_2 macros */ \
/*	distance between fft data elements is 16K */ \
/*	do 4 macros each processing 32 data values */ \
\
start_timer(18); \
/* 4 iterations of 2*clm */ \
for(unsigned int loopA = 4; loopA; loopA--) { \
	for(unsigned int loopB = clm*2; loopB; ) { \
		do{ \
			x4cl_eight_reals_unfft_2(rsi, 64, 4*blkdst, 8*blkdst); \
			L2prefetch128(u8ptr(rsi-64+rbp)); \
			L2prefetch128(u8ptr(rsi-64+4*blkdst+rbp)); \
			pfing(rbp += -64+8*blkdst); \
		}while(--loopB & 1); \
		pfing(rbp += -16*blkdst+2*64);	/* Restore prefetching offset */ \
	} \
	rsi += -2*clm*64+blkdst;	/* Next source pointer */ \
} \
\
/* Do 48 four_complex_unfft macros */ \
/*	distance between fft data elements is 16K */ \
/*	do 12 macros each processing 32 data values */ \
\
rsi += -4*blkdst+16*blkdst;/* Next source pointer */ \
rdi = (uintptr_t)g->sincos3;/* Load sin/cos pointer */ \
/* eax = 3*256+4; */ /* 3 iterations of 4 iters of 2*clm */ \
for(unsigned int loopA = 3; loopA; loopA--) { \
	for(unsigned int loopB = 4; loopB; loopB--) { \
		for(unsigned int loopC = clm*2; loopC; ) { \
			do{ \
				x4cl_four_complex_unfft(rsi, 64, 4*blkdst, 8*blkdst, rdi); \
				L2prefetch128(u8ptr(rsi-64+rbp)); \
				L2prefetch128(u8ptr(rsi-64+4*blkdst+rbp)); \
				pfing(rbp += -64+8*blkdst); \
			}while(--loopC & 1); \
			pfing(rbp += -16*blkdst+2*64);	/* Restore prefetching offset */ \
		} \
		rsi += -2*clm*64+blkdst;	/* Next source pointer */ \
	} \
	rsi += -4*blkdst+16*blkdst;/* Next source pointer */ \
	rdi += XMM_SCD;		/* Next sine/cosine pointer */ \
} \
rsi += -4*16*blkdst;	/* Restore source pointer */ \
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
/* Do 64 eight_reals_last_unfft macros */ \
/*	distance between fft data elements is 64K */ \
/*	do 16 macros each processing 32 data values */ \
\
start_timer(19); \
/* 16 iterations of 2*clm */ \
for(unsigned int loopA = 16; loopA; loopA--) { \
	for(unsigned int loopB = clm*2; loopB; ) { \
		do{ \
			x4cl_eight_reals_last_unfft(rsi, 64, 16*blkdst, 32*blkdst); \
			/*	L2prefetch128 [rsi-64][rbp] */ \
			L2prefetch128(u8ptr(rsi-64+16*blkdst+rbp)); \
			pfing(rbp += -64+32*blkdst); \
		}while(--loopB & 1); \
		pfing(rbp += -64*blkdst+2*64);	/* Restore prefetching offset */ \
	} \
	rsi += -2*clm*64+blkdst;	/* Next source pointer */ \
} \
end_timer(19); \
\
/* Normalize these values */ \
\
pass1_normalize(0, b0b); \
goto	endlp; \
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
/* Fire up auxiliary threads to do pass 1 forward FFT */ \
\
pass1:	pass1_forward_fft_setup(b0b); \
\
/* Do 64 eight_reals_first_fft macros (FFT levels 1,2,3) */ \
/*	distance between fft data elements is 64K */ \
/*	do 32 iterations each processing 16 data values */ \
\
b0b:	start_timer(21); \
rsi = (uintptr_t)g->data_addr;/* Load source address for prefetching */ \
rbx = g->DIST_TO_FFTSRCARG; \
pfing(rbp = (uintptr_t)g->data_prefetch);/* Set prefetch distance */ \
pfing(rbp -= rsi); \
/* 32 iterations of 2*clm */ \
/* rdx = 0 */ \
for(unsigned int loopA = 32; loopA; loopA--) { \
	for(unsigned int loopB = clm*2; loopB; ) { \
		do{ \
			x2cl_eight_reals_first_fft(rsi, 64, 32*blkdst); \
			L2prefetch128(u8ptr(rsi-64+rbp)); \
			pfing(rbp += -64+32*blkdst); \
		}while(--loopB & 1); \
		pfing(rbp += -64*blkdst+2*64);	/* Restore prefetching offset */ \
	} \
	rsi += -2*clm*64+blkdst;	/* Next source pointer */ \
} \
rsi += -32*blkdst;		/* Restore source pointer */ \
end_timer(21); \
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
\
/* Do 8 eight_reals_fft_2 macros */ \
/* Do 8 nop_two_two_complex_fft_2 macros */ \
/* Do 16 four_complex_fft macros */ \
/*	distance between fft data elements is 16K */ \
/*	do 8 macros each processing 32 data values */ \
\
start_timer(22); \
rdi = (uintptr_t)g->sincos3;/* Load sin/cos pointer */ \
/* 8 iterations of 2*clm */ \
for(unsigned int loopA = 8; loopA; loopA--) { \
	for(unsigned int loopB = clm*2; loopB; ) { \
		do{ \
			x4cl_eight_reals_fft_2(rsi, 64, 8*blkdst, 16*blkdst, rdi); \
			L2prefetch128(u8ptr(rsi-64+rbp)); \
			/*	L2prefetch128 [rsi-64+8*blkdst][rbp] */ \
			pfing(rbp += -64+16*blkdst); \
		}while(--loopB & 1); \
		pfing(rbp += -32*blkdst+2*64);	/* Restore prefetching offset */ \
	} \
	rsi += -2*clm*64+blkdst;	/* Next source pointer */ \
} \
\
/* Do 32 four_complex_fft macros */ \
/*	distance between fft data elements is 16K */ \
/*	do 8 macros each processing 32 data values */ \
\
rsi += -8*blkdst+32*blkdst;/* Next source pointer */ \
rdi += XMM_SCD;		/* Next sine/cosine pointer */ \
/* 8 iterations of 2*clm */ \
for(unsigned int loopA = 8; loopA; loopA--) { \
	for(unsigned int loopB = clm*2; loopB; ) { \
		do{ \
			x4cl_four_complex_fft(rsi, 64, 8*blkdst, 16*blkdst, rdi); \
			L2prefetch128(u8ptr(rsi-64+rbp)); \
			/*	L2prefetch128 [rsi-64+8*blkdst][rbp] */ \
			pfing(rbp += -64+16*blkdst); \
		}while(--loopB & 1); \
		pfing(rbp += -32*blkdst+2*64);	/* Restore prefetching offset */ \
	} \
	rsi += -2*clm*64+blkdst;	/* Next source pointer */ \
} \
rsi += -8*blkdst-32*blkdst;/* Restore source pointer */ \
end_timer(22); \
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
\
/* Do 4 eight_reals_fft_1 macros */ \
/* Do 4 two_two_complex_fft macros */ \
/*	distance between fft data elements is 8K */ \
/*	do 4 macros each processing 16 data values */ \
\
start_timer(23); \
rdi = (uintptr_t)g->sincos2;/* Load sin/cos pointer */ \
/* 4 iterations of 2*clm */ \
for(unsigned int loopA = 4; loopA; loopA--) { \
	for(unsigned int loopB = clm*2; loopB; ) { \
		do{ \
			x2cl_eight_reals_fft_1(rsi, 64, 4*blkdst); \
		}while(--loopB & 1); \
		L2prefetch128(u8ptr(rsi-128+rbp)); \
	} \
	rsi += -2*clm*64+blkdst;	/* Next source pointer */ \
} \
\
/* Do 56 two_two_complex_fft macros */ \
/*	distance between fft data elements is 8K */ \
/*	do 28 macros each processing 16 data values */ \
\
rsi += -4*blkdst+8*blkdst;	/* Load source pointer */ \
rdi += XMM_SCD;		/* Next sine/cosine pointer */ \
/* eax = 7*256+4; */ /* 7 iterations of 4 iters of 2*clm */ \
for(unsigned int loopA = 7; loopA; loopA--) { \
	for(unsigned int loopB = 4; loopB; loopB--) { \
		for(unsigned int loopC = clm*2; loopC; ) { \
			do{ \
				x2cl_two_complex_fft(rsi, 64, 4*blkdst, rdi); \
			}while(--loopC & 1); \
			L2prefetch128(u8ptr(rsi-128+rbp)); \
		} \
		rsi += -2*clm*64+blkdst;	/* Next source pointer */ \
	} \
	rsi += -4*blkdst+8*blkdst;	/* Next source pointer */ \
	rdi += 2*XMM_SCD;		/* Next sine/cosine pointer */ \
} \
rsi += -8*8*blkdst;	/* Restore source pointer */ \
end_timer(23); \
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
	x4cl_eight_reals_fft_2(rsi, 64, blkdst, 2*blkdst, rdi); \
} \
\
/* Do 60 four_complex_fft macros */ \
/*	distance between fft data elements is 2K */ \
/*	do 15 macros each processing 32 data values */ \
\
rsi += -2*clm*64+4*blkdst;	/* Next source pointer */ \
rdi += XMM_SCD;		/* Next sine/cosine pointer */ \
/* 15 iterations of 2*clm */ \
for(unsigned int loopA = 15; loopA; loopA--) { \
	for(unsigned int loopB = clm*2; loopB; ) { \
		do{ \
			x4cl_four_complex_fft(rsi, 64, blkdst, 2*blkdst, rdi); \
			L2prefetch128(u8ptr(rsi-64+rbp)); \
			L2prefetch128(u8ptr(rsi-64+blkdst+rbp)); \
			pfing(rbp += -64+2*blkdst); \
		}while(--loopB & 1); \
		pfing(rbp += -4*blkdst+2*64);	/* Restore prefetching offset */ \
	} \
	rsi += -2*clm*64+4*blkdst;	/* Next source pointer */ \
	rdi += 2*XMM_SCD;		/* Next sine/cosine pointer */ \
} \
end_timer(24); \
\
/* Work on next 4 sets of 256 values */ \
\
endlp:	xpass1_get_next_block(pass2, c0b, b0b); /* Get next block and jump */ \
done:	xfft_footer_hg1; \
\
\
\
/* Perform a 512K-element negacyclic FFT.  This is done in two passes.  Pass 1 */ \
/* does 8 levels, pass 2 uses common code to perform the last 11 levels. */ \
\

#define hg_pass1levels8complex4(procname, xpass2name, pass2_levels, clm) \
/* (fftname, pass2_levels, clmarg) */ \
\
/* set_FFT_constants(fftname, 1, 1, pass2_levels, clmarg); */ \
xfft_header_hg1(procname, xpass2name, pass2_levels, clm); \
\
/* Do the last 11 FFT levels */ \
/*	do 64 pass2_11_levels_complex */ \
\
pass2:	start_timer(0); \
xpass2name(g); \
end_timer(0); \
if(g->ffttype == 1)		/* We're done if FFTing only */ \
	return; /* goto done; */ \
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
rdi = (uintptr_t)g->sincos1;/* Load sin/cos pointer */ \
/* 32 iterations of 2*clm */ \
/* rdx = 0 */ \
pfing(rbp = (uintptr_t)g->data_prefetch);/* Set prefetch distance */ \
pfing(rbp -= rsi); \
for(unsigned int loopA = 32; loopA; loopA--) { \
	for(unsigned int loopB = 2*clm; loopB; loopB--) { \
		x2cl_four_complex_unfft(rsi, 64, blkdst); \
	} \
	rsi += -2*clm*64+2*blkdst;	/* Next source pointer */ \
	rdi += XMM_SCD;		/* Next sine/cosine pointer */ \
} \
rsi += -32*2*blkdst;	/* Restore source pointer */ \
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
/* eax = 8*256+2; */ /* 8 iterations of 2 iters of 2*clm */ \
for(unsigned int loopA = 8; loopA; loopA--) { \
	for(unsigned int loopB = 2; loopB; loopB--) { \
		for(unsigned int loopC = clm*2; loopC; ) { \
			do{ \
				x4cl_four_complex_unfft(rsi, 64, 2*blkdst, 4*blkdst, rdi); \
				L2prefetch128(u8ptr(rsi-64+rbp)); \
				/*	L2prefetch128 [rsi-64+2*blkdst][rbp] */ \
				pfing(rbp += -64+4*blkdst); \
			}while(--loopC & 1); \
			pfing(rbp += -8*blkdst+2*64);	/* Restore prefetching offset */ \
		} \
		rsi += -2*clm*64+blkdst;	/* Next source pointer */ \
	} \
	rsi += -2*blkdst+8*blkdst;	/* Next source pointer */ \
	rdi += XMM_SCD;		/* Next sine/cosine pointer */ \
} \
rsi += -8*8*blkdst;	/* Restore source pointer */ \
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
/* eax = 4*256+8; */ /* 4 iterations of 8 iters of 2*clm */ \
for(unsigned int loopA = 4; loopA; loopA--) { \
	for(unsigned int loopB = 8; loopB; loopB--) { \
		for(unsigned int loopC = clm*2; loopC; ) { \
			do{ \
				x2cl_two_complex_unfft(rsi, 64, 8*blkdst); \
			}while(--loopC & 1); \
			L2prefetch128(u8ptr(rsi-128+rbp)); \
		} \
		rsi += -2*clm*64+blkdst;	/* Next source pointer */ \
	} \
	rsi += -8*blkdst+16*blkdst;/* Next source pointer */ \
	rdi += XMM_SCD;		/* Next sine/cosine pointer */ \
} \
rsi += -4*16*blkdst;	/* Restore source pointer */ \
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
/* 16 iterations of 2*clm */ \
for(unsigned int loopA = 16; loopA; loopA--) { \
	for(unsigned int loopB = clm*2; loopB; ) { \
		do{ \
			x4cl_four_complex_last_unfft(rsi, 64, 16*blkdst, 32*blkdst, 16*128); \
			/*	L2prefetch128 [rsi-64][rbp] */ \
			L2prefetch128(u8ptr(rsi-64+16*blkdst+rbp)); \
			pfing(rbp += -64+32*blkdst); \
		}while(--loopB & 1); \
		pfing(rbp += -64*blkdst+2*64);	/* Restore prefetching offset */ \
	} \
	rsi += -2*clm*64+blkdst;	/* Next source pointer */ \
	rdi += 128;		/* Next premultiplier pointer */ \
} \
end_timer(19); \
\
/* Normalize these values */ \
\
pass1_normalize(0, b0b); \
goto	endlp; \
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
/* Fire up auxiliary threads to do pass 1 forward FFT */ \
\
pass1:	pass1_forward_fft_setup(b0b); \
\
/* Do 64 four_complex_first_fft macros (FFT levels 1,2,3) */ \
/*	distance between fft data elements is 64K */ \
/*	do 32 iterations each processing 16 data values */ \
\
b0b:	start_timer(21); \
rsi = (uintptr_t)g->data_addr;/* Load source address for prefetching */ \
rbx = g->DIST_TO_FFTSRCARG; \
pfing(rbp = (uintptr_t)g->data_prefetch);/* Set prefetch distance */ \
pfing(rbp -= rsi); \
rdi = (uintptr_t)g->plus1_premults;/* Address of premultiplier table */ \
/* 32 iterations of 2*clm */ \
/* rdx = 0 */ \
for(unsigned int loopA = 32; loopA; loopA--) { \
	for(unsigned int loopB = clm*2; loopB; ) { \
		do{ \
			x2cl_four_complex_first_fft(rsi, 64, 32*blkdst); \
			L2prefetch128(u8ptr(rsi-64+rbp)); \
			pfing(rbp += -64+32*blkdst); \
		}while(--loopB & 1); \
		pfing(rbp += -64*blkdst+2*64);	/* Restore prefetching offset */ \
	} \
	rsi += -2*clm*64+blkdst;	/* Next source pointer */ \
	rdi += 128;		/* Next premultiplier pointer */ \
} \
rsi += -32*blkdst;		/* Restore source pointer */ \
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
/* eax = 2*256+16; */ /* 2 iters of 16 iterations of 2*clm */ \
for(unsigned int loopA = 2; loopA; loopA--) { \
	for(unsigned int loopB = 16; loopB; loopB--) { \
		for(unsigned int loopC = clm*2; loopC; ) { \
			do{ \
				x2cl_two_complex_fft(rsi, 64, 16*blkdst, rdi); \
			}while(--loopC & 1); \
			L2prefetch128(u8ptr(rsi-128+rbp)); \
		} \
		rsi += -2*clm*64+blkdst;	/* Next source pointer */ \
	} \
	rsi += -16*blkdst+32*blkdst;/* Next source pointer */ \
	rdi += 2*XMM_SCD;		/* Next sin/cos pointer */ \
} \
rsi += -2*32*blkdst;	/* Restore source pointer */ \
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
/* eax = 4*256+4; */ /* 4 iterations of 4 iters of 2*clm */ \
for(unsigned int loopA = 4; loopA; loopA--) { \
	for(unsigned int loopB = 4; loopB; loopB--) { \
		for(unsigned int loopC = clm*2; loopC; ) { \
			do{ \
				x4cl_four_complex_fft(rsi, 64, 4*blkdst, 8*blkdst, rdi); \
				L2prefetch128(u8ptr(rsi-64+rbp)); \
				/*	L2prefetch128 [rsi-64+4*blkdst][rbp] */ \
				pfing(rbp += -64+8*blkdst); \
			}while(--loopC & 1); \
			pfing(rbp += -16*blkdst+2*64);	/* Restore prefetching offset */ \
		} \
		rsi += -2*clm*64+blkdst;	/* Next source pointer */ \
	} \
	rsi += -4*blkdst+16*blkdst;/* Next source pointer */ \
	rdi += 2*XMM_SCD;		/* Next sine/cosine pointer */ \
} \
rsi += -4*16*blkdst;	/* Restore source pointer */ \
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
/* 16 iterations of 2*clm */ \
for(unsigned int loopA = 16; loopA; loopA--) { \
	for(unsigned int loopB = clm*2; loopB; ) { \
		do{ \
			x4cl_four_complex_fft(rsi, 64, blkdst, 2*blkdst, rdi); \
			L2prefetch128(u8ptr(rsi-64+rbp)); \
			L2prefetch128(u8ptr(rsi-64+blkdst+rbp)); \
			pfing(rbp += -64+2*blkdst); \
		}while(--loopB & 1); \
		pfing(rbp += -4*blkdst+2*64);	/* Restore prefetching offset */ \
	} \
	rsi += -2*clm*64+4*blkdst;	/* Next source pointer */ \
	rdi += 2*XMM_SCD;		/* Next sine/cosine pointer */ \
} \
end_timer(24); \
\
/* Work on next 4 sets of 256 values */ \
\
endlp:	xpass1_get_next_block(pass2, c0b, b0b); /* Get next block and jump */ \
done:	xfft_footer_hg1; \
\
