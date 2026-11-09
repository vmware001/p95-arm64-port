#pragma once

#define xfft_header_r(procname, xpass2name, pass2_levels, clm) \
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
		/*start_timer(1);*/ \
		g->thread_work_routine = procname##_i1; \
		pass1_wake_up_threads(g); \
		procname##_i (g, 1); \
	}else{ \
		procname##_i (g, 4); /* Start pass2 */ \
	} \
} \
void procname##_i (struct gwasm_data *__restrict g, int branch1) \
{ \
	uintptr_t rcx,rdx,rsi,rdi; \
	vec2f64 xmm0,xmm1,xmm2,xmm3,xmm4,xmm5,xmm6,xmm7; \
	rsi = (uintptr_t)g->DESTARG; \
	/*rbx = g->DIST_TO_FFTSRCARG;*/ \
	/* xmult.mac  set_FFT_constants */ \
	const uintptr_t clmblkdst = clm*128; \
	const uintptr_t clmblkdst8 = clmblkdst*8+128; \
	const uintptr_t p2cl = (pass2_levels >= 64) ? (pass2_levels / 4) : (1<<(pass2_levels-2)); \
	const uintptr_t blkdst_tmp1 = (p2cl*128+p2cl/64*128); \
	const uintptr_t blkdst = (blkdst_tmp1 % 256 == 0) ? (blkdst_tmp1-128) : blkdst_tmp1; \
	P1(procname); \
	do{ \
		if(branch1 == 0) goto c0b; /* Do another inverse FFT/norm/FFT  */ \
		if(branch1 & 1) goto b0b; /* Do another forward FFT */ \

#define xfft_footer_r \
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




#define pass1_forward_fft_setup(LABEL1)

#define pass1_inverse_fft_setup(LABEL1, procname) \
	g->next_block = 1; \
	/*start_timer(1);*/ \
	g->thread_work_routine = procname##_i0 ; \
	pass1_wake_up_threads(g); \


#if 0  // looks like extra parameters gone unused
#define newxloop1nopre(clm, label2, incr_rsi, pm_incr_reg, pm_incr_amt, sc_incr_reg, sc_incr_amt) \
	rsi += -2*clm*64+incr_rsi; /* Next source pointer */ \
	pm_incr_reg += pm_incr_amt; /* Next premultiplier pointer */ \
	sc_incr_reg += sc_incr_amt; /* Next sin/cos pointer */

#define newxloop1(clm, label2, incr_rsi, pm_incr_reg, pm_incr_amt, sc_incr_reg, sc_incr_amt) \
newxloop1nopre(clm, label2, incr_rsi, pm_incr_reg, pm_incr_amt, sc_incr_reg, sc_incr_amt)
#define newxloop2(clm, label2, incr_rsi, pm_incr_reg, pm_incr_amt, sc_incr_reg, sc_incr_amt) \
newxloop1nopre(clm, label2, incr_rsi, pm_incr_reg, pm_incr_amt, sc_incr_reg, sc_incr_amt)
#define newxloop4(clm, label2, incr_rsi, pm_incr_reg, pm_incr_amt, sc_incr_reg, sc_incr_amt) \
newxloop1nopre(clm, label2, incr_rsi, pm_incr_reg, pm_incr_amt, sc_incr_reg, sc_incr_amt)
#define newxloop1doublepre(clm, label2, label1, incr_rsi, pm_incr_reg, pm_incr_amt, sc_incr_reg, sc_incr_amt) \
newxloop1nopre(clm, label2, incr_rsi, pm_incr_reg, pm_incr_amt, sc_incr_reg, sc_incr_amt)
#define newxloop1quadpre(clm, label2, label1, incr_rsi, pm_incr_reg, pm_incr_amt, sc_incr_reg, sc_incr_amt) \
newxloop1nopre(clm, label2, incr_rsi, pm_incr_reg, pm_incr_amt, sc_incr_reg, sc_incr_amt)

#else
#define newxloop1nopre(clm, label2, incr_rsi) \
rsi += -2*clm*64+incr_rsi;
#define newxloop1(clm, label2, incr_rsi) \
newxloop1nopre(clm, label2, incr_rsi)
#define newxloop2(clm, label2, incr_rsi) \
newxloop1nopre(clm, label2, incr_rsi)
#define newxloop4(clm, label2, incr_rsi) \
newxloop1nopre(clm, label2, incr_rsi)
#define newxloop1doublepre(clm, label2, label1, incr_rsi) \
newxloop1nopre(clm, label2, incr_rsi)
#define newxloop1quadpre(clm, label2, label1, incr_rsi) \
newxloop1nopre(clm, label2, incr_rsi)
#endif

// clang register thrashing slowdown workaround in x86_64
//#define REGZ(A) __asm__ __volatile__("" : "=r"(A) : "0"(A));
#define REGZ(A)

/* Copyright 2009-2023 - Mersenne Research, Inc.  All rights reserved */ \
/* Author:  George Woltman */ \
/* Email: woltman@alum.mit.edu */ \
/* */ \
/* These macros do the first pass of standard radix-4 2 pass DJB FFTs using SSE2 instructions. */ \
/* */ \
\
/* ******************************************************** */ \
/* ******************************************************** */ \
/* ******************  PASS 1 MACROS  ********************* */ \
/* ******************************************************** */ \
/* ******************************************************** */ \
\
\
/************************************************************************ */ \
/* Pass 1 macros */ \
/* */ \
/* p2cl = number of 128-byte cache lines processed from a pass 2 block */ \
/* clm = cache line multiplier (num cache lines processed each iteration) */ \
/************************************************************************ */ \
\
/* Note:  All memory layouts are for 64 byte L1 cache lines.  The high */ \
/* order XMM words are not shown (they are 512 higher than the low */ \
/* order XMM word assuming 10 levels done in pass 2). */ \
\
\
/* Perform a pass 1 FFT of length 32 doubles. */ \
/* Memory examples are for a 32K FFT with 10 levels done in pass 2. */ \
/* Sin/cos data needed: 4*2*clm*XMM_SCD3 + 2*clm*XMM_SCD2 + 2*clm*XMM_SCD2 = 8*2*clm 64-byte cache lines */ \
\

#define r4_pass1sc32(procname, xpass2name, pass2_levels, clm) /* (fftname, pass2_levels, clmarg) */ \
\
/* set_FFT_constants(fftname, 0, 0, pass2_levels, clmarg); */ \
xfft_header_r(procname, xpass2name, pass2_levels, clm); \
\
/* Do pass 2 (the last 10 FFT levels if this is a 32K FFT) */ \
\
/* pass2: */	start_timer(0); \
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
/* On input the 64-byte cache lines hold these data values if pass 2 swizzles: */ \
/*	0K	+512	1K	+512	2K	+512	3K	+512 */ \
/*	4K */ \
/*	... */ \
/* On input the 64-byte cache lines hold these data values if pass 1 swizzles: */ \
/*	0K	+2K	512	+2K	1K	+2K	1536	+2K */ \
/*	4K */ \
/*	... */ \
/* On output the 64-byte cache lines hold these data values (+512 in high word): */ \
/*	0K	2K	4K	6K */ \
/*	1K	... */ \
/*	8K	... */ \
/*	... */ \
\
/* Do 1 eight_reals_unfft operations */ \
/*	distance between fft data elements is 1K */ \
/*	do 1 macros each processing 16 data values */ \
/* Do 3 four_complex_unfft4 operations */ \
/*	distance between fft data elements is 1K */ \
/*	do 3 macros each processing 16 data values */ \
\
c0b:	start_timer(16); \
rsi = (uintptr_t)g->data_addr;/* Get source address */ \
rcx = (uintptr_t)g->scratch_area;	/* Get address of scratch area */ \
rdx = (uintptr_t)g->premult_addr;/* Real data sin/cos pointer */ \
rdi = (rdx+2*clm*XMM_SCD2); /* Complex data sin/cos ptr */ \
pfing(rbp = (rdx+(8-4)*2*clm*64)); /* Finish prefetch the sin/cos/premult data */ \
loops_init_prefetch_touch(4*2*clm, 64, 1, rbp, 4096, c9a); /* 4*2*clm iterations */ \
xtouch(xptr(rbp));			/* Preload TLB */ \
for(unsigned int loopA = 2*clm; loopA; loopA--) { \
	r4_sg2cl_eight_reals_djbunfft(rsi, 64, blkdst, rcx, 64, clmblkdst, rdi, rdx); \
	rdi += XMM_SCD2;		/* Next premult/sin/cos data */ \
	rdx += XMM_SCD2;		/* Next sine/cosine pointer */ \
} \
rsi += -2*clm*64+2*blkdst;	/* Next source pointer */ \
rcx += -2*clm*64+2*clmblkdst; /* Next dest pointer */ \
rdi += -2*clm*XMM_SCD2;	/* Restore premult/sine/cosine pointer */ \
\
loops_reset; \
xtouch(xptr(rbp+4096));		/* Preload TLB */ \
for(unsigned int loopA = 3; loopA; loopA--) { \
	for(unsigned int loopB = 2*clm; loopB; loopB--) { \
		r4_sg2cl_four_complex_djbunfft(rsi, 64, blkdst, rcx, 64, clmblkdst, rdi); \
		rdi += XMM_SCD2;		/* Next premult/sin/cos data */ \
	} \
	rsi += -2*clm*64+2*blkdst;	/* Next source pointer */ \
	rcx += -2*clm*64+2*clmblkdst; /* Next dest pointer */ \
	rdi += -2*clm*XMM_SCD2;	/* Restore premult/sin/cos data pointer */ \
} \
rsi = (uintptr_t)g->scratch_area;	/* Restore source pointer */ \
rdi += 2*clm*XMM_SCD2;	/* Next sin/cos pointer */ \
end_timer(16); \
\
/* Do inverse FFT levels 1,2,3 */ \
/* On input the 64-byte cache lines hold these data values: */ \
/*	0K	2K	4K	6K */ \
/*	1K	... */ \
/*	8K	... */ \
/*	... */ \
/* On output the 64-byte cache lines hold these data values: */ \
/*	0K	8K	16K	24K */ \
/*	1K	... */ \
/*	... */ \
/*	7K */ \
\
/* Do 4 eight_reals_last_unfft operations */ \
/*	distance between fft data elements is 4K */ \
/*	do 2 macros each processing 32 data values */ \
\
start_timer(20); \
set_data_prefetch_ptrs;		/* Init prefetching pointers */ \
loops_init(2*2*clm);		/* 2 iterations of 2*clm */ \
r4_x4cl_eight_reals_last_unfft_preload; \
for(unsigned int loopA = 2*2*clm; loopA; ) { \
	xtouch(xptr(rcx+rbp+(clm-1)*128));	/* Preload TLB */ \
	do{ \
		r4_x4cl_eight_reals_last_unfft(rsi, 64, 2*clmblkdst, 4*clmblkdst, rdi, 2*2*clm*XMM_SCD3); \
		rdi += XMM_SCD3;		/* Next sin/cos pointer */ \
	}while(--loopA & (2*clm-1)); \
	newxloop1doublepre(clm, c1b, c1a, clmblkdst); /* 2*clm iterations */ \
} \
end_timer(20); \
\
/* Normalize these values */ \
\
pass1_normalize(1, yes_fft); \
\
/* Copy data from the scratch area back to the FFT data area */ \
\
copy_scratch_data_to_fft(8, clm); \
goto	endlp;			/* Skip FFT code if flag not set */ \
\
/* Fire up auxiliary threads to do pass 1 forward FFT */ \
\
/* pass1: */	pass1_forward_fft_setup(b0b); \
\
/* Copy data from the FFT data area back to the scratch area */ \
\
b0b:	copy_fft_data_to_scratch(8, clm); \
\
/* Do FFT levels 1,2,3 */ \
/* */ \
/* On input the 64-byte cache lines hold these data values: */ \
/*	0K	8K	16K	24K */ \
/*	1K	... */ \
/*	... */ \
/*	7K */ \
/* On output the 64-byte cache lines hold these data values: */ \
/*	0K	4K	8K	12K */ \
/*	1K	... */ \
/*	... */ \
/*	3K	... */ \
/*	16K	... */ \
/*	... */ \
\
/* Do 4 eight_reals_first_fft operations */ \
/*	distance between fft data elements is 4K */ \
/*	do 4 iterations each processing 16 data values */ \
\
yes_fft:start_timer(22); \
rsi = (uintptr_t)g->scratch_area;	/* Get address of scratch area */ \
rdi = (uintptr_t)g->premult_addr;/* Address of premultiplier/sin/cos table */ \
rdi += 2*clm*XMM_SCD2 + 2*clm*XMM_SCD2; \
set_data_prefetch_ptrs;		/* Init prefetching pointers */ \
pfing(rcx += 4*blkdst);		/* We've already prefetched 4 blocks */ \
loops_init(4*2*clm);		/* 4 iters of 2*clm */ \
r4_x2cl_eight_reals_first_fft_scratch_preload; \
for(unsigned int loopA = 4*2*clm; loopA; ) { \
	xtouch(xptr(rcx+rbp+(clm-1)*128));	/* Preload TLB */ \
	do{ \
		r4_x2cl_eight_reals_first_fft_scratch(rsi, 64, 4*clmblkdst, rdi); \
		rdi += XMM_SCD3;		/* Next sin/cos pointer */ \
	}while(--loopA & (2*clm-1)); \
	newxloop1(clm, b1b, clmblkdst);	/* 2*clm iters */ \
} \
rsi += -4*clmblkdst;	/* Restore source pointer */ \
end_timer(22); \
\
/* Do FFT levels 4,5 */ \
/* */ \
/* On input the 64-byte cache lines hold these data values (+512 in high word): */ \
/*	0K	4K	8K	12K */ \
/*	1K	... */ \
/*	... */ \
/*	3K	... */ \
/*	16K	... */ \
/*	... */ \
/* On output the 64-byte cache lines hold these data values if pass 2 swizzles: */ \
/*	0K	+512	1K	+512	2K	+512	3K	+512 */ \
/*	4K	... */ \
/*	... */ \
/* On output the 64-byte cache lines hold these data values if pass 1 swizzles: */ \
/*	0K	+2K	512	+2K	1K	+2K	1536	+2K */ \
/*	4K	... */ \
/*	... */ \
\
/* Do 1 eight_reals_fft operations */ \
/* Do 1 four_complex_djbfft operations */ \
/*	distance between fft data elements is 1K */ \
/*	do 1 macros each processing 32 data values */ \
/* Do 2 four_complex_djbfft operations */ \
/*	distance between fft data elements is 1K */ \
/*	do 1 macros each processing 32 data values */ \
\
start_timer(26); \
rdi += -4*2*clm*XMM_SCD3-2*clm*XMM_SCD2; /* Complex data sin/cos ptr */ \
rdx = (rdi-2*clm*XMM_SCD2); /* Real data sin/cos ptr */ \
rcx = (uintptr_t)g->data_addr;/* Destination pointer */ \
loops_init_prefetch(2*2*clm, 128, 1, rbp); /* 2 iters of 2*clm */ \
pfing(rbp = (uintptr_t)g->premult_prefetch);/* Load premult/sin/cos prefetch pointer */ \
r4_sg4cl_eight_reals_four_complex_djbfft_preload; \
for(unsigned int loopA = 2*clm; loopA; loopA--) { \
	r4_sg4cl_eight_reals_four_complex_djbfft(rsi, 64, clmblkdst, 2*clmblkdst, rcx, 64, blkdst, 2*blkdst, rdi, rdx); \
	rdi += XMM_SCD2;		/* Next sine/cosine pointer */ \
	rdx += XMM_SCD2;		/* Next sine/cosine pointer */ \
} \
rsi += -2*clm*64+4*clmblkdst; /* Next source pointer */ \
rcx += -2*clm*64+4*blkdst;	/* Next dest pointer */ \
rdi += -2*clm*XMM_SCD2;	/* Restore premult/sine/cosine pointer */ \
\
loops_reset; \
r4_sg4cl_four_complex_djbfft_preload; \
{ \
	for(unsigned int loopB = 2*clm; loopB; loopB--) { \
		r4_sg4cl_four_complex_djbfft(rsi, 64, clmblkdst, 2*clmblkdst, rcx, 64, blkdst, 2*blkdst, rdi); \
		rdi += XMM_SCD2;		/* Next premult/sin/cos ptr */ \
	} \
	rsi += -2*clm*64+4*clmblkdst; /* Next source pointer */ \
	rcx += -2*clm*64+4*blkdst;	/* Next dest pointer */ \
	rdi += -2*clm*XMM_SCD2;	/* Restore premult/sine/cosine pointer */ \
} \
end_timer(26); \
\
/* Work on next 4*clm sets of 32 values */ \
\
endlp:	xpass1_get_next_block(pass2, c0b, b0b); /* Get next block and jump */ \
/* done: */ xfft_footer_r; \
\
\
\
/* Perform a pass 1 FFT of length 64 doubles. */ \
/* Memory examples are for a 64K FFT with 10 levels done in pass 2. */ \
/* Sin/cos data needed: 8*2*clm*XMM_SCD3 + 2*clm*XMM_SCD4 + 2*clm*XMM_SCD8 = 18*2*clm 64-byte cache lines */ \
/* BUG:	we under-prefetch the sin/cos data by 6*2*clm*64 bytes */ \
\

#define r4_pass1sc64(procname, xpass2name, pass2_levels, clm) /* (fftname, pass2_levels, clmarg) */ \
\
/* set_FFT_constants(fftname, 0, 0, pass2_levels, clmarg); */ \
xfft_header_r(procname, xpass2name, pass2_levels, clm); \
\
/* Do pass 2 (the last 10 FFT levels if this is a 64K FFT) */ \
\
/* pass2: */	start_timer(0); \
xpass2name(g); \
end_timer(0); \
if(g->ffttype == 1)		/* We're done if FFTing only */ \
	return; /* goto done; */ \
\
/* Do the inverse FFT */ \
\
pass1_inverse_fft_setup(c0b, procname); \
\
/* Do inverse FFT levels 4,5,6 */ \
/* On input the 64-byte cache lines hold these data values if pass 2 swizzles: */ \
/*	0K	+512	1K	+512	2K	+512	3K	+512 */ \
/*	4K */ \
/*	... */ \
/* On input the 64-byte cache lines hold these data values if pass 1 swizzles: */ \
/*	0K	+2K	512	+2K	1K	+2K	1536	+2K */ \
/*	4K */ \
/*	... */ \
/* On output the 64-byte cache lines hold these data values (+512 in high word): */ \
/*	0K	4K	8K	12K */ \
/*	1K	... */ \
/*	... */ \
/*	3K */ \
/*	16K	... */ \
/*	... */ \
\
/* Do 1 sixteen_reals_unfft operations */ \
/*	distance between fft data elements is 1K */ \
/*	do 1 macros each processing 32 data values */ \
/* Do 3 eight_complex_unfft8 operations */ \
/*	distance between fft data elements is 1K */ \
/*	do 3 macros each processing 32 data values */ \
\
c0b:	start_timer(16); \
rsi = (uintptr_t)g->data_addr;/* Get source address */ \
rcx = (uintptr_t)g->scratch_area;	/* Get address of scratch area */ \
rdx = (uintptr_t)g->premult_addr;/* Real data sin/cos ptr */ \
rdi = (rdx+2*clm*XMM_SCD4+128); /* Complex data sin/cos pointer and use -128 */ \
;					/* addressing to generate 1 byte offsets in radix-8 macro */ \
pfing(rbp = (rdx+(18-8)*2*clm*64)); /* Finish prefetching the sin/cos/premult data */ \
loops_init_prefetch(4*2*clm, 128, 1, rbp); /* 4*2*clm iterations */ \
/*bug	xtouchpages rbp, 4*2*clm * 128 / 1 ;; Prefetch TLBs */ \
for(unsigned int loopA = 2*clm; loopA; loopA--) { \
	r8_sg4cl_sixteen_reals_unfft8(rsi, 64, blkdst, 2*blkdst, rcx, 64, clmblkdst, 2*clmblkdst, rdi-128, rdx); \
	rdi += XMM_SCD8;		/* Next premult/sin/cos data */ \
	rdx += XMM_SCD4;		/* Next sine/cosine pointer */ \
} \
rsi += -2*clm*64+4*blkdst;	/* Next source pointer */ \
rcx += -2*clm*64+4*clmblkdst; /* Next dest pointer */ \
rdi += -2*clm*XMM_SCD8;	/* Restore sin/cos pointer */ \
\
loops_reset; \
c9a:	xtouch(xptr(rbp+4096));		/* Prefetch TLB */ \
for(unsigned int loopA = 2*2-1; loopA; ) { \
	do{ \
		for(unsigned int loopC = 2*clm; loopC; loopC--) { \
			r8_sg4cl_eight_complex_unfft8(rsi, 64, blkdst, 2*blkdst, rcx, 64, clmblkdst, 2*clmblkdst, rdi-128); \
			rdi += XMM_SCD8;		/* Next premult/sin/cos data */ \
		} \
		rsi += -2*clm*64+4*blkdst;	/* Next source pointer */ \
		rcx += -2*clm*64+4*clmblkdst; /* Next dest pointer */ \
		rdi += -2*clm*XMM_SCD8;	/* Restore premult/sin/cos data pointer */ \
	}while(--loopA & 1); \
	rcx += -2*4*clmblkdst+clmblkdst8; \
} \
rsi = (uintptr_t)g->scratch_area;	/* Restore source pointer */ \
rdi += -128+2*clm*XMM_SCD8; /* Next sin/cos pointer */ \
end_timer(16); \
\
/* Do inverse FFT levels 1,2,3 */ \
/* On input the 64-byte cache lines hold these data values: */ \
/*	0K	4K	8K	12K */ \
/*	1K	... */ \
/*	... */ \
/*	3K */ \
/*	16K	... */ \
/*	... */ \
/* On output the 64-byte cache lines hold these data values: */ \
/*	0K	16K	32K	48K */ \
/*	1K	... */ \
/*	... */ \
/*	15K */ \
\
/* Do 8 eight_reals_last_unfft operations */ \
/*	distance between fft data elements is 8K */ \
/*	do 4 macros each processing 32 data values */ \
\
start_timer(20); \
set_data_prefetch_ptrs;		/* Init prefetching pointers */ \
loops_init(4*2*clm);		/* 4 iterations of 2*clm */ \
r4_x4cl_eight_reals_last_unfft_preload; \
for(unsigned int loopA = 4*2*clm; loopA; ) { \
	xtouch(xptr(rcx+rbp+(clm-1)*128));	/* Preload TLB */ \
	do{ \
		r4_x4cl_eight_reals_last_unfft(rsi, 64, 4*clmblkdst, clmblkdst8, rdi, 4*2*clm*XMM_SCD3); \
		rdi += XMM_SCD3;		/* Next sin/cos ptr */ \
	}while(--loopA & (2*clm-1)); \
	newxloop1doublepre(clm, c1b, c1a, clmblkdst); /* 2*clm iterations */ \
} \
end_timer(20); \
\
/* Normalize these values */ \
\
pass1_normalize(1, yes_fft); \
\
/* Copy data from the scratch area back to the FFT data area */ \
\
copy_scratch_data_to_fft(16, clm); \
goto	endlp;			/* Skip FFT code if flag not set */ \
\
/* Fire up auxiliary threads to do pass 1 forward FFT */ \
\
/* pass1: */	pass1_forward_fft_setup(b0b); \
\
/* Copy data from the FFT data area back to the scratch area */ \
\
b0b:	copy_fft_data_to_scratch(16, clm); \
\
/* Do FFT levels 1,2,3 */ \
/* */ \
/* On input the 64-byte cache lines hold these data values: */ \
/*	0K	16K	32K	48K */ \
/*	1K	... */ \
/*	... */ \
/*	15K */ \
/* On output the 64-byte cache lines hold these data values: */ \
/*	0K	8K	16K	24K */ \
/*	1K	... */ \
/*	... */ \
/*	7K	... */ \
/*	32K	... */ \
/*	... */ \
\
/* Do 8 eight_reals_first_fft operations */ \
/*	distance between fft data elements is 8K */ \
/*	do 8 iterations each processing 16 data values */ \
\
yes_fft:start_timer(22); \
rsi = (uintptr_t)g->scratch_area;	/* Get address of scratch area */ \
rdi = (uintptr_t)g->premult_addr;/* Address of premultiplier/sin/cos table */ \
rdi += 2*clm*XMM_SCD4 + 2*clm*XMM_SCD8; \
set_data_prefetch_ptrs;		/* Init prefetching pointers */ \
pfing(rcx += 8*blkdst);		/* We've already prefetched 8 blocks */ \
loops_init(8*2*clm);		/* 8 iters of 2*clm */ \
r4_x2cl_eight_reals_first_fft_scratch_preload; \
for(unsigned int loopA = 8*2*clm; loopA; ) { \
	xtouch(xptr(rcx+rbp+(clm-1)*128));	/* Preload TLB */ \
	do{ \
		r4_x2cl_eight_reals_first_fft_scratch(rsi, 64, clmblkdst8, rdi); \
		rdi += XMM_SCD3;		/* Next sin/cos ptr */ \
	}while(--loopA & (2*clm-1)); \
	newxloop1(clm, b1b, clmblkdst);	/* 2*clm iters */ \
} \
rsi += -8*clmblkdst;	/* Restore source pointer */ \
end_timer(22); \
\
/* Do FFT levels 4,5,6 */ \
/* */ \
/* On input the 64-byte cache lines hold these data values (+512 in high word): */ \
/*	0K	8K	16K	24K */ \
/*	1K	... */ \
/*	... */ \
/*	7K	... */ \
/*	32K	... */ \
/*	... */ \
/* On output the 64-byte cache lines hold these data values if pass 2 swizzles: */ \
/*	0K	+512	1K	+512	2K	+512	3K	+512 */ \
/*	4K	... */ \
/*	... */ \
/* On output the 64-byte cache lines hold these data values if pass 1 swizzles: */ \
/*	0K	+2K	512	+2K	1K	+2K	1536	+2K */ \
/*	4K	... */ \
/*	... */ \
\
/* Do 1 sixteen_reals_fft operations */ \
/* Do 1 eight_complex_fft8 operations */ \
/*	distance between fft data elements is 1K */ \
/*	do 1 macros each processing 64 data values */ \
/* Do 2 eight_complex_fft8 operations */ \
/*	distance between fft data elements is 1K */ \
/*	do 1 macros each processing 64 data values */ \
\
start_timer(26); \
rdi += -8*2*clm*XMM_SCD3-2*clm*XMM_SCD8+128; /* Complex data sin/cos ptr and use -128 */ \
;					/* addressing to generate 1 byte offsets in radix-8 macro */ \
rdx = (rdi-128-2*clm*XMM_SCD4); /* Real data sin/cos ptr */ \
rcx = (uintptr_t)g->data_addr;/* Destination pointer */ \
loops_init_prefetch(2*2*clm, 128, 1, rbp); /* 2 iters of 2*clm */ \
pfing(rbp = (uintptr_t)g->premult_prefetch);/* Load premult/sin/cos prefetch pointer */ \
r8_sg8cl_sixteen_reals_eight_complex_fft8_preload; \
for(unsigned int loopA = 2*clm; loopA; loopA--) { \
	r8_sg8cl_sixteen_reals_eight_complex_fft8(rsi, 64, clmblkdst, 2*clmblkdst, 4*clmblkdst, rcx, 64, blkdst, 2*blkdst, 4*blkdst, rdi-128, rdx); \
	rdi += XMM_SCD8;		/* Next sine/cosine pointer */ \
	rdx += XMM_SCD4;		/* Next sine/cosine pointer */ \
	REGZ(rdi); \
} \
rsi += -2*clm*64+clmblkdst8; /* Next source pointer */ \
rcx += -2*clm*64+8*blkdst;	/* Next dest pointer */ \
rdi += -2*clm*XMM_SCD8;	/* Restore premult/sine/cosine pointer */ \
\
loops_reset; \
{ \
	for(unsigned int loopB = 2*clm; loopB; loopB--) { \
		r8_sg8cl_eight_complex_fft8(rsi, 64, clmblkdst, 2*clmblkdst, 4*clmblkdst, rcx, 64, blkdst, 2*blkdst, 4*blkdst, rdi-128); \
		rdi += XMM_SCD8;		/* Next premult/sin/cos ptr */ \
	} \
	rsi += -2*clm*64+clmblkdst8; /* Next source pointer */ \
	rcx += -2*clm*64+8*blkdst;	/* Next dest pointer */ \
	rdi += -2*clm*XMM_SCD8;	/* Restore premult/sine/cosine pointer */ \
} \
end_timer(26); \
\
/* Work on next 4*clm sets of 64 values */ \
\
endlp:	xpass1_get_next_block(pass2, c0b, b0b); /* Get next block and jump */ \
/* done: */ xfft_footer_r; \
\
\
\
/* Perform a pass 1 FFT of length 80 doubles. */ \
/* Memory examples are for a 80K FFT with 10 levels done in pass 2. */ \
/* Sin/cos data needed: 4*2*clm*XMM_SCD9 + 2*clm*XMM_SCD2 + 2*clm*XMM_SCD2 = 20*2*clm 64-byte cache lines */ \
\

#define r4_pass1sc80(procname, xpass2name, pass2_levels, clm) /* (fftname, pass2_levels, clmarg) */ \
\
/* set_FFT_constants(fftname, 0, 0, pass2_levels, clmarg); */ \
xfft_header_r(procname, xpass2name, pass2_levels, clm); \
\
/* Do pass 2 (the last 10 FFT levels if this is a 80K FFT) */ \
\
/* pass2: */	start_timer(0); \
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
/* On input the 64-byte cache lines hold these data values if pass 2 swizzles: */ \
/*	0K	+512	1K	+512	2K	+512	3K	+512 */ \
/*	4K */ \
/*	... */ \
/* On input the 64-byte cache lines hold these data values if pass 1 swizzles: */ \
/*	0K	+2K	512	+2K	1K	+2K	1536	+2K */ \
/*	4K */ \
/*	... */ \
/* On output the 64-byte cache lines hold these data values (+512 in high word): */ \
/*	0K	2K	4K	6K */ \
/*	1K	... */ \
/*	8K	... */ \
/*	... */ \
\
/* Do 1 eight_reals_unfft operations */ \
/*	distance between fft data elements is 1K */ \
/*	do 1 macros each processing 16 data values */ \
/* Do 9 four_complex_djbunfft operations */ \
/*	distance between fft data elements is 1K */ \
/*	do 9 macros each processing 16 data values */ \
\
c0b:	start_timer(16); \
rsi = (uintptr_t)g->data_addr;/* Get source address */ \
rcx = (uintptr_t)g->scratch_area;	/* Get address of scratch area */ \
rdx = (uintptr_t)g->premult_addr;/* Real data sin/cos pointer */ \
rdi = (rdx+2*clm*XMM_SCD2); /* Complex data sin/cos ptr */ \
pfing(rbp = (rdx+(20-10)*2*clm*64)); /* Finish prefetch the sin/cos/premult data */ \
loops_init_prefetch(10*2*clm, 64, 1, rbp); /* 10*2*clm iterations */ \
/*	xtouchpages [rbp]		;; Prefetch TLBs */ \
for(unsigned int loopA = 2*clm; loopA; loopA--) { \
	r4_sg2cl_eight_reals_djbunfft(rsi, 64, blkdst, rcx, 64, clmblkdst, rdi, rdx); \
	rdi += XMM_SCD2;		/* Next premult/sin/cos data */ \
	rdx += XMM_SCD2;		/* Next sine/cosine pointer */ \
} \
rsi += -2*clm*64+2*blkdst;	/* Next source pointer */ \
rcx += -2*clm*64+2*clmblkdst; /* Next dest pointer */ \
rdi += -2*clm*XMM_SCD2;	/* Restore premult/sine/cosine pointer */ \
\
loops_reset; \
for(unsigned int loopA = 9; loopA; loopA--) { \
	for(unsigned int loopB = 2*clm; loopB; loopB--) { \
		r4_sg2cl_four_complex_djbunfft(rsi, 64, blkdst, rcx, 64, clmblkdst, rdi); \
		rdi += XMM_SCD2;		/* Next premult/sin/cos data */ \
	} \
	rsi += -2*clm*64+2*blkdst;	/* Next source pointer */ \
	rcx += -2*clm*64+2*clmblkdst; /* Next dest pointer */ \
	rdi += -2*clm*XMM_SCD2;	/* Restore premult/sin/cos data pointer */ \
} \
rsi = (uintptr_t)g->scratch_area;	/* Restore source pointer */ \
rdi += 2*clm*XMM_SCD2;	/* Next sin/cos pointer */ \
end_timer(16); \
\
/* Do inverse FFT levels 1,2,3,4,5 (actually 4.32 levels) */ \
/* On input the 64-byte cache lines hold these data values: */ \
/*	0K	2K	4K	6K */ \
/*	1K	... */ \
/*	8K	... */ \
/*	... */ \
/* On output the 64-byte cache lines hold these data values: */ \
/*	0K	20K	40K	60K */ \
/*	1K	... */ \
/*	... */ \
/*	19K */ \
\
/* Do 16 20_reals_last_unfft operations */ \
/*	distance between fft data elements is 4K */ \
/*	do 2 macros each processing 80 data values */ \
\
start_timer(20); \
set_data_prefetch_ptrs;		/* Init prefetching pointers */ \
loops_init(2*2*clm);		/* 2 iterations of 2*clm */ \
r5_x10cl_20_reals_last_unfft_preload; \
for(unsigned int loopA = 2*2*clm; loopA; ) { \
	xtouch(xptr(rcx+rbp+(clm-1)*128));	/* Preload TLB */ \
	do{ \
		r5_x10cl_20_reals_last_unfft(rsi, 64, 2*clmblkdst, rdi, 2*2*clm*XMM_SCD9); \
		rdi += XMM_SCD9;		/* Next sin/cos pointer */ \
	}while(--loopA & (2*clm-1)); \
	newxloop1quadpre(clm, c1b, c1a, clmblkdst); /* 2*clm iterations */ \
} \
end_timer(20); \
\
/* Normalize these values */ \
\
pass1_normalize(1, yes_fft); \
\
/* Copy data from the scratch area back to the FFT data area */ \
\
copy_scratch_no_pad_data_to_fft(20, clm); \
goto	endlp;			/* Skip FFT code if flag not set */ \
\
/* Fire up auxiliary threads to do pass 1 forward FFT */ \
\
/* pass1: */	pass1_forward_fft_setup(b0b); \
\
/* Copy data from the FFT data area back to the scratch area */ \
\
b0b:	copy_fft_data_to_scratch_no_pad(20, clm); \
\
/* Do FFT levels 1,2,3,4,5 (actually 4.32 levels) */ \
/* */ \
/* On input the 64-byte cache lines hold these data values: */ \
/*	0K	20K	40K	60K */ \
/*	1K	... */ \
/*	... */ \
/*	19K */ \
/* On output the 64-byte cache lines hold these data values: */ \
/*	0K	4K	8K	12K */ \
/*	1K	... */ \
/*	... */ \
/*	3K	... */ \
/*	16K	... */ \
/*	... */ \
\
/* Do 16 20_reals_first_fft operations */ \
/*	distance between fft data elements is 4K */ \
/*	do 4 iterations each processing 40 data values */ \
\
yes_fft:start_timer(22); \
rsi = (uintptr_t)g->scratch_area;	/* Get address of scratch area */ \
rdi = (uintptr_t)g->premult_addr;/* Address of premultiplier/sin/cos table */ \
rdi += 2*clm*XMM_SCD2 + 2*clm*XMM_SCD2; \
set_data_prefetch_ptrs;		/* Init prefetching pointers */ \
pfing(rcx += (8-4)*blkdst);	/* We've already prefetched 8 blocks, but */ \
;					/* we'll prefetch 4 blocks too many below */ \
loops_init(4*2*clm);		/* 4 iters of 2*clm */ \
r5_x5cl_20_reals_first_fft_scratch_preload; \
for(unsigned int loopA = 4*2*clm; loopA; ) { \
	xtouch(xptr(rcx+rbp+(clm-1)*128));	/* Preload TLB */ \
	do{ \
		r5_x5cl_20_reals_first_fft_scratch(rsi, 64, 4*clmblkdst, rdi); \
		rdi += XMM_SCD9;		/* Next sin/cos pointer */ \
	}while(--loopA & (2*clm-1)); \
	newxloop1quadpre(clm, b1b, b1a, clmblkdst); /* 2*clm iters prefetching 8*clm cache lines */ \
} \
rsi += -4*clmblkdst;	/* Restore source pointer */ \
end_timer(22); \
\
/* Do FFT levels 6,7 */ \
/* */ \
/* On input the 64-byte cache lines hold these data values (+512 in high word): */ \
/*	0K	4K	8K	12K */ \
/*	1K	... */ \
/*	... */ \
/*	3K	... */ \
/*	16K	... */ \
/*	... */ \
/* On output the 64-byte cache lines hold these data values if pass 2 swizzles: */ \
/*	0K	+512	1K	+512	2K	+512	3K	+512 */ \
/*	4K	... */ \
/*	... */ \
/* On output the 64-byte cache lines hold these data values if pass 1 swizzles: */ \
/*	0K	+2K	512	+2K	1K	+2K	1536	+2K */ \
/*	4K	... */ \
/*	... */ \
\
/* Do 1 eight_reals_fft operations */ \
/* Do 1 four_complex_djbfft operations */ \
/*	distance between fft data elements is 1K */ \
/*	do 1 macros each processing 32 data values */ \
/* Do 8 four_complex_djbfft operations */ \
/*	distance between fft data elements is 1K */ \
/*	do 4 macros each processing 32 data values */ \
\
start_timer(26); \
rdi += -4*2*clm*XMM_SCD9-2*clm*XMM_SCD2; /* Complex data sin/cos ptr */ \
rdx = (rdi-2*clm*XMM_SCD2); /* Real data sin/cos ptr */ \
rcx = (uintptr_t)g->data_addr;/* Destination pointer */ \
loops_init_prefetch(5*2*clm, 128, 1, rbp); /* 5 iters of 2*clm */ \
pfing(rbp = (uintptr_t)g->premult_prefetch);/* Load premult/sin/cos prefetch pointer */ \
r4_sg4cl_eight_reals_four_complex_djbfft_preload; \
for(unsigned int loopA = 2*clm; loopA; loopA--) { \
	r4_sg4cl_eight_reals_four_complex_djbfft(rsi, 64, clmblkdst, 2*clmblkdst, rcx, 64, blkdst, 2*blkdst, rdi, rdx); \
	rdi += XMM_SCD2;		/* Next sine/cosine pointer */ \
	rdx += XMM_SCD2;		/* Next sine/cosine pointer */ \
} \
rsi += -2*clm*64+4*clmblkdst; /* Next source pointer */ \
rcx += -2*clm*64+4*blkdst;	/* Next dest pointer */ \
rdi += -2*clm*XMM_SCD2;	/* Restore premult/sine/cosine pointer */ \
\
loops_reset; \
r4_sg4cl_four_complex_djbfft_preload; \
for(unsigned int loopA = 4; loopA; loopA--) { \
	for(unsigned int loopB = 2*clm; loopB; loopB--) { \
		r4_sg4cl_four_complex_djbfft(rsi, 64, clmblkdst, 2*clmblkdst, rcx, 64, blkdst, 2*blkdst, rdi); \
		rdi += XMM_SCD2;		/* Next premult/sin/cos ptr */ \
	} \
	rsi += -2*clm*64+4*clmblkdst; /* Next source pointer */ \
	rcx += -2*clm*64+4*blkdst;	/* Next dest pointer */ \
	rdi += -2*clm*XMM_SCD2;	/* Restore premult/sine/cosine pointer */ \
} \
end_timer(26); \
\
/* Work on next 4*clm sets of 80 values */ \
\
endlp:	xpass1_get_next_block(pass2, c0b, b0b); /* Get next block and jump */ \
/* done: */ xfft_footer_r; \
\
\
\
/* Perform a pass 1 FFT of length 96 doubles. */ \
/* Memory examples are for a 96K FFT with 10 levels done in pass 2. */ \
/* Sin/cos data needed: 12*2*clm*XMM_SCD3 + 4*2*clm*XMM_SCD1 + 4*2*clm*XMM_SCD1 + 2*clm*XMM_SCD2 + 2*clm*XMM_SCD2 = 24*2*clm 64-byte cache lines */ \
\

#define r4_pass1sc96(procname, xpass2name, pass2_levels, clm) /* (fftname, pass2_levels, clmarg) */ \
\
/* set_FFT_constants(fftname, 0, 0, pass2_levels, clmarg); */ \
xfft_header_r(procname, xpass2name, pass2_levels, clm); \
\
/* Do pass 2 (the last 10 FFT levels if this is a 128K FFT) */ \
\
/* pass2: */	start_timer(0); \
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
/* On input the 64-byte cache lines hold these data values if pass 2 swizzles: */ \
/*	0K	+512	1K	+512	2K	+512	3K	+512 */ \
/*	4K */ \
/*	... */ \
/* On input the 64-byte cache lines hold these data values if pass 1 swizzles: */ \
/*	0K	+2K	512	+2K	1K	+2K	1536	+2K */ \
/*	4K */ \
/*	... */ \
/* On output the 64-byte cache lines hold these data values (+512 in high word): */ \
/*	0K	2K	4K	6K */ \
/*	1K	... */ \
/*	8K	... */ \
/*	... */ \
\
/* Do 1 eight_reals_unfft operations */ \
/*	distance between fft data elements is 1K */ \
/*	do 1 macros each processing 16 data values */ \
/* Do 11 four_complex_djbunfft operations */ \
/*	distance between fft data elements is 1K */ \
/*	do 11 macros each processing 16 data values */ \
\
c0b:	start_timer(16); \
rsi = (uintptr_t)g->data_addr;/* Get source address */ \
rcx = (uintptr_t)g->scratch_area;	/* Get address of scratch area */ \
rdx = (uintptr_t)g->premult_addr;/* Real data sin/cos pointer */ \
rdi = (rdx+2*clm*XMM_SCD2); /* Complex data sin/cos ptr */ \
pfing(rbp = (rdx+(24-12)*2*clm*64)); /* Finish prefetch the sin/cos/premult data */ \
loops_init_prefetch(12*2*clm, 64, 1, rbp); /* 12*2*clm iterations */ \
/*	xtouchpages [rbp]		;; Prefetch TLBs */ \
for(unsigned int loopA = 2*clm; loopA; loopA--) { \
	r4_sg2cl_eight_reals_djbunfft(rsi, 64, blkdst, rcx, 64, clmblkdst, rdi, rdx); \
	rdi += XMM_SCD2;		/* Next premult/sin/cos data */ \
	rdx += XMM_SCD2;		/* Next sine/cosine pointer */ \
} \
rsi += -2*clm*64+2*blkdst;	/* Next source pointer */ \
rcx += -2*clm*64+2*clmblkdst; /* Next dest pointer */ \
rdi += -2*clm*XMM_SCD2;	/* Restore premult/sine/cosine pointer */ \
\
loops_reset; \
for(unsigned int loopA = 11; loopA; loopA--) { \
	for(unsigned int loopB = 2*clm; loopB; loopB--) { \
		r4_sg2cl_four_complex_djbunfft(rsi, 64, blkdst, rcx, 64, clmblkdst, rdi); \
		rdi += XMM_SCD2;		/* Next premult/sin/cos data */ \
	} \
	rsi += -2*clm*64+2*blkdst;	/* Next source pointer */ \
	rcx += -2*clm*64+2*clmblkdst; /* Next dest pointer */ \
	rdi += -2*clm*XMM_SCD2;	/* Restore premult/sin/cos data pointer */ \
} \
rsi = (uintptr_t)g->scratch_area;	/* Restore source pointer */ \
rdx = (rdi+2*clm*XMM_SCD2); /* Next real data sin/cos pointer */ \
end_timer(16); \
\
/* Do inverse FFT levels 4,5 (actually 1.585 levels) */ \
/* On input the 64-byte cache lines hold these data values: */ \
/*	0K	2K	4K	6K */ \
/*	1K	... */ \
/*	8K	... */ \
/*	... */ \
/* On output the 64-byte cache lines hold these data values: */ \
/*	0K	6K	12K	18K */ \
/*	1K	... */ \
/*	... */ \
/*	5K */ \
/*	24K	... */ \
/*	... */ \
\
/* Do 4 six_reals_unfft operations */ \
/*	distance between fft data elements is 4K */ \
/*	do 2 macros each processing 24 data values */ \
/* Do 12 three_complex_unfft operations */ \
/*	distance between fft data elements is 4K */ \
/*	do 6 macros each processing 24 data values */ \
\
start_timer(17); \
rdi = (rdx+4*2*clm*XMM_SCD1); /* Next complex data sin/cos ptr */ \
set_data_prefetch_ptrs;		/* Init prefetching pointers */ \
loops_init(8*2*clm);		/* 8*2*clm iterations */ \
\
xtouch(xptr(rcx+rbp+(clm-1)*128));	/* Preload TLB */ \
for(unsigned int loopA = 2*2*clm; loopA; ) { \
	do{ \
		r3_x3cl_six_reals_unfft(rsi, 64, 2*clmblkdst, rdi, 2*2*clm*XMM_SCD1, rdx, 2*2*clm*XMM_SCD1); \
		rdi += XMM_SCD1;		/* Next sine/cosine pointer */ \
		rdx += XMM_SCD1;		/* Next sine/cosine pointer */ \
	}while(--loopA & (2*clm-1)); \
	newxloop2(clm, c6b, clmblkdst);	/* 2*2*clm iterations */ \
} \
rsi += -2*clmblkdst+6*clmblkdst; /* Next source pointer */ \
rdi += -2*2*clm*XMM_SCD1;	/* Restore sine/cosine pointer */ \
\
loops_reset; \
r3_x3cl_three_complex_djbunfft_preload; \
for(unsigned int loopA = 3*2*2*clm; loopA; ) { \
	xtouch(xptr(rcx+rbp+(clm-1)*128));	/* Preload TLB */ \
	do{ \
		do{ \
			r3_x3cl_three_complex_djbunfft(rsi, 64, 2*clmblkdst, rdi, 2*2*clm*XMM_SCD1); \
			rdi += XMM_SCD1;		/* Next clm's sin/cos data */ \
		}while(--loopA & (2*clm-1)); \
		newxloop2(clm, c7b, clmblkdst);	/* 2*2*clm iterations */ \
	}while(loopA & (2*2*clm-1)); \
	rsi += -2*clmblkdst+6*clmblkdst; /* Next source pointer */ \
	rdi += -2*2*clm*XMM_SCD1;	/* Restore sine/cosine pointer */ \
} \
rsi += -4*6*clmblkdst;	/* Restore source pointer */ \
rdi += 4*2*clm*XMM_SCD1;	/* Next sin/cos pointer */ \
end_timer(17); \
\
/* Do inverse FFT levels 1,2,3 */ \
/* On input the 64-byte cache lines hold these data values: */ \
/*	0K	6K	12K	18K */ \
/*	1K	... */ \
/*	... */ \
/*	5K */ \
/*	24K	... */ \
/*	... */ \
/* On output the 64-byte cache lines hold these data values: */ \
/*	0K	24K	48K	72K */ \
/*	1K	... */ \
/*	... */ \
/*	23K */ \
\
/* Do 12 eight_reals_last_unfft operations */ \
/*	distance between fft data elements is 12K */ \
/*	do 6 macros each processing 32 data values */ \
\
start_timer(20); \
loops_init(6*2*clm);		/* 6 iterations of 2*clm */ \
r4_x4cl_eight_reals_last_unfft_preload; \
for(unsigned int loopA = 6*2*clm; loopA; ) { \
	xtouch(xptr(rcx+rbp+(clm-1)*128));	/* Preload TLB */ \
	do{ \
		r4_x4cl_eight_reals_last_unfft(rsi, 64, 6*clmblkdst, 12*clmblkdst, rdi, 6*2*clm*XMM_SCD3); \
		rdi += XMM_SCD3;		/* Next sin/cos pointer */ \
	}while(--loopA & (2*clm-1)); \
	newxloop1(clm, c1b, clmblkdst);	/* 2*clm iterations */ \
} \
end_timer(20); \
\
/* Normalize these values */ \
\
pass1_normalize(1, yes_fft); \
\
/* Copy data from the scratch area back to the FFT data area */ \
\
copy_scratch_no_pad_data_to_fft(24, clm); \
goto	endlp;			/* Skip FFT code if flag not set */ \
\
/* Fire up auxiliary threads to do pass 1 forward FFT */ \
\
/* pass1: */	pass1_forward_fft_setup(b0b); \
\
/* Copy data from the FFT data area back to the scratch area */ \
\
b0b:	copy_fft_data_to_scratch_no_pad(24, clm); \
\
/* Do FFT levels 1,2,3 */ \
/* */ \
/* On input the 64-byte cache lines hold these data values: */ \
/*	0K	24K	48K	72K */ \
/*	1K	... */ \
/*	... */ \
/*	23K */ \
/* On output the 64-byte cache lines hold these data values: */ \
/*	0K	12K	24K	36K */ \
/*	1K	... */ \
/*	... */ \
/*	11K	... */ \
/*	48K	... */ \
/*	... */ \
\
/* Do 12 eight_reals_first_fft operations */ \
/*	distance between fft data elements is 12K */ \
/*	do 12 iterations each processing 16 data values */ \
\
yes_fft:start_timer(22); \
rsi = (uintptr_t)g->scratch_area;	/* Get address of scratch area */ \
rdi = (uintptr_t)g->premult_addr;/* Address of premultiplier/sin/cos table */ \
rdi += 4*2*clm*XMM_SCD1 + 4*2*clm*XMM_SCD1 + 2*clm*XMM_SCD2 + 2*clm*XMM_SCD2; \
set_data_prefetch_ptrs;		/* Init prefetching pointers */ \
pfing(rcx += 10*blkdst);		/* We've already prefetched 10 blocks */ \
loops_init(12*2*clm);		/* 6 iters of 2*2*clm */ \
r4_x2cl_eight_reals_first_fft_scratch_preload; \
for(unsigned int loopA = 6*2*2*clm; loopA; ) { \
	xtouch(xptr(rcx+rbp+(clm-1)*128));	/* Preload TLB */ \
	do{ \
		do{ \
			r4_x2cl_eight_reals_first_fft_scratch(rsi, 64, 12*clmblkdst, rdi); \
			rdi += XMM_SCD3;		/* Next sin/cos pointer */ \
		}while(--loopA & (2*clm-1)); \
		newxloop2(clm, b1b, clmblkdst);	/* 2*2*clm iters prefetching 2*clm cache lines */ \
	}while(loopA & (2*2*clm-1)); \
} \
rsi += -6*2*clmblkdst;	/* Restore source pointer */ \
end_timer(22); \
\
/* Do FFT levels 4,5 (actually 1.585 levels) */ \
/* */ \
/* On input the 64-byte cache lines hold these data values: */ \
/*	0K	12K	24K	36K */ \
/*	1K	... */ \
/*	... */ \
/*	11K	... */ \
/*	48K	... */ \
/*	... */ \
/* On output the 64-byte cache lines hold these data values: */ \
/*	0K	4K	8K	12K */ \
/*	1K	... */ \
/*	... */ \
/*	3K	... */ \
/*	16K	... */ \
/*	... */ \
\
/* Do 4 six_reals_fft operations */ \
/* Do 4 three_complex_fft operations */ \
/*	distance between fft data elements is 4K */ \
/*	do 4 macros each processing 24 data values */ \
/* Do 8 three_complex_fft operations */ \
/*	distance between fft data elements is 4K */ \
/*	do 4 macros each processing 24 data values */ \
\
start_timer(25); \
rdi += -12*2*clm*XMM_SCD3-4*2*clm*XMM_SCD1; /* Complex data sin/cos ptr */ \
rdx = (rdi-4*2*clm*XMM_SCD1); /* Real data sin/cos ptr */ \
loops_init(8*2*clm);		/* 2 iterations of 4 iters of 2*clm */ \
r3_x3cl_six_reals_three_complex_djbfft_preload; \
for(unsigned int loopA = 4*2*clm; loopA; ) { \
	xtouch(xptr(rcx+rbp+(clm-1)*128));	/* Preload TLB */ \
	do{ \
		r3_x3cl_six_reals_three_complex_djbfft(rsi, 64, 4*clmblkdst, rdi, rdx); \
		rdi += XMM_SCD1;		/* Next sine/cosine pointer */ \
		rdx += XMM_SCD1;		/* Next sine/cosine pointer */ \
	}while(--loopA & (2*clm-1)); \
	newxloop1(clm, b6b, clmblkdst);	/* 2*clm iterations */ \
} \
rsi += -4*clmblkdst+12*clmblkdst; \
rdi += -4*2*clm*XMM_SCD1;	/* Restore sin/cos ptr */ \
\
loops_reset; \
r3_x3cl_three_complex_djbfft_preload; \
{ \
	for(unsigned int loopB = 4*2*clm; loopB; ) { \
		xtouch(xptr(rcx+rbp+(clm-1)*128));	/* Preload TLB */ \
		do{ \
			r3_x3cl_three_complex_djbfft(rsi, 64, 4*clmblkdst, rdi); \
			rdi += XMM_SCD1;		/* Next sine/cosine pointer */ \
		}while(--loopB & (2*clm-1)); \
		newxloop1(clm, b7b, clmblkdst);	/* 2*clm iterations */ \
	} \
	rsi += -4*clmblkdst+12*clmblkdst; /* Next source pointer */ \
	rdi += -4*2*clm*XMM_SCD1;	/* Restore sin/cos ptr */ \
} \
rsi += -2*12*clmblkdst;	/* Restore source pointer */ \
end_timer(25); \
\
/* Do FFT levels 6,7 */ \
/* */ \
/* On input the 64-byte cache lines hold these data values (+512 in high word): */ \
/*	0K	4K	8K	12K */ \
/*	1K	... */ \
/*	... */ \
/*	3K	... */ \
/*	16K	... */ \
/*	... */ \
/* On output the 64-byte cache lines hold these data values if pass 2 swizzles: */ \
/*	0K	+512	1K	+512	2K	+512	3K	+512 */ \
/*	4K	... */ \
/*	... */ \
/* On output the 64-byte cache lines hold these data values if pass 1 swizzles: */ \
/*	0K	+2K	512	+2K	1K	+2K	1536	+2K */ \
/*	4K	... */ \
/*	... */ \
\
/* Do 1 eight_reals_fft operations */ \
/* Do 1 four_complex_djbfft operations */ \
/*	distance between fft data elements is 1K */ \
/*	do 1 macros each processing 32 data values */ \
/* Do 10 four_complex_djbfft operations */ \
/*	distance between fft data elements is 1K */ \
/*	do 5 macros each processing 32 data values */ \
\
start_timer(26); \
rdi += -4*2*clm*XMM_SCD1-2*clm*XMM_SCD2; /* Complex data sin/cos ptr */ \
rdx = (rdi-2*clm*XMM_SCD2); /* Real data sin/cos ptr */ \
rcx = (uintptr_t)g->data_addr;/* Destination pointer */ \
loops_init_prefetch(6*2*clm, 128, 1, rbp); /* 6 iters of 2*clm */ \
pfing(rbp = (uintptr_t)g->premult_prefetch);/* Load premult/sin/cos prefetch pointer */ \
r4_sg4cl_eight_reals_four_complex_djbfft_preload; \
for(unsigned int loopA = 2*clm; loopA; loopA--) { \
	r4_sg4cl_eight_reals_four_complex_djbfft(rsi, 64, clmblkdst, 2*clmblkdst, rcx, 64, blkdst, 2*blkdst, rdi, rdx); \
	rdi += XMM_SCD2;		/* Next sine/cosine pointer */ \
	rdx += XMM_SCD2;		/* Next sine/cosine pointer */ \
} \
rsi += -2*clm*64+4*clmblkdst; /* Next source pointer */ \
rcx += -2*clm*64+4*blkdst;	/* Next dest pointer */ \
rdi += -2*clm*XMM_SCD2;	/* Restore premult/sine/cosine pointer */ \
\
loops_reset; \
r4_sg4cl_four_complex_djbfft_preload; \
for(unsigned int loopA = 5; loopA; loopA--) { \
	for(unsigned int loopB = 2*clm; loopB; loopB--) { \
		r4_sg4cl_four_complex_djbfft(rsi, 64, clmblkdst, 2*clmblkdst, rcx, 64, blkdst, 2*blkdst, rdi); \
		rdi += XMM_SCD2;		/* Next premult/sin/cos ptr */ \
	} \
	rsi += -2*clm*64+4*clmblkdst; /* Next source pointer */ \
	rcx += -2*clm*64+4*blkdst;	/* Next dest pointer */ \
	rdi += -2*clm*XMM_SCD2;	/* Restore premult/sine/cosine pointer */ \
} \
end_timer(26); \
\
/* Work on next 4*clm sets of 96 values */ \
\
endlp:	xpass1_get_next_block(pass2, c0b, b0b); /* Get next block and jump */ \
/* done: */ xfft_footer_r; \
\
\
\
/* Perform a pass 1 FFT of length 112 doubles. */ \
/* Memory examples are for a 112K FFT with 10 levels done in pass 2. */ \
/* Sin/cos data needed: 4*2*clm*XMM_SCD13 + 2*clm*XMM_SCD2 + 2*clm*XMM_SCD2 = 28*2*clm 64-byte cache lines */ \
/* NOTE:  We under-prefetch FFT data by 4 blocks */ \
\

#define r4_pass1sc112(procname, xpass2name, pass2_levels, clm) /* (fftname, pass2_levels, clmarg) */ \
\
/* set_FFT_constants(fftname, 0, 0, pass2_levels, clmarg); */ \
xfft_header_r(procname, xpass2name, pass2_levels, clm); \
\
/* Do pass 2 (the last 10 FFT levels if this is a 112K FFT) */ \
\
/* pass2: */	start_timer(0); \
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
/* On input the 64-byte cache lines hold these data values if pass 2 swizzles: */ \
/*	0K	+512	1K	+512	2K	+512	3K	+512 */ \
/*	4K */ \
/*	... */ \
/* On input the 64-byte cache lines hold these data values if pass 1 swizzles: */ \
/*	0K	+2K	512	+2K	1K	+2K	1536	+2K */ \
/*	4K */ \
/*	... */ \
/* On output the 64-byte cache lines hold these data values (+512 in high word): */ \
/*	0K	2K	4K	6K */ \
/*	1K	... */ \
/*	8K	... */ \
/*	... */ \
\
/* Do 1 eight_reals_unfft operations */ \
/*	distance between fft data elements is 1K */ \
/*	do 1 macros each processing 16 data values */ \
/* Do 13 four_complex_djbunfft operations */ \
/*	distance between fft data elements is 1K */ \
/*	do 13 macros each processing 16 data values */ \
\
c0b:	start_timer(16); \
rsi = (uintptr_t)g->data_addr;/* Get source address */ \
rcx = (uintptr_t)g->scratch_area;	/* Get address of scratch area */ \
rdx = (uintptr_t)g->premult_addr;/* Real data sin/cos pointer */ \
rdi = (rdx+2*clm*XMM_SCD2); /* Complex data sin/cos ptr */ \
pfing(rbp = (rdx+(28-14)*2*clm*64)); /* Finish prefetch the sin/cos/premult data */ \
loops_init_prefetch(14*2*clm, 64, 1, rbp); /* 14*2*clm iterations */ \
/*	xtouchpages [rbp]		;; Prefetch TLBs */ \
for(unsigned int loopA = 2*clm; loopA; loopA--) { \
	r4_sg2cl_eight_reals_djbunfft(rsi, 64, blkdst, rcx, 64, clmblkdst, rdi, rdx); \
	rdi += XMM_SCD2;		/* Next premult/sin/cos data */ \
	rdx += XMM_SCD2;		/* Next sine/cosine pointer */ \
} \
rsi += -2*clm*64+2*blkdst;	/* Next source pointer */ \
rcx += -2*clm*64+2*clmblkdst; /* Next dest pointer */ \
rdi += -2*clm*XMM_SCD2;	/* Restore premult/sine/cosine pointer */ \
\
loops_reset; \
for(unsigned int loopA = 13; loopA; loopA--) { \
	for(unsigned int loopB = 2*clm; loopB; loopB--) { \
		r4_sg2cl_four_complex_djbunfft(rsi, 64, blkdst, rcx, 64, clmblkdst, rdi); \
		rdi += XMM_SCD2;		/* Next premult/sin/cos data */ \
	} \
	rsi += -2*clm*64+2*blkdst;	/* Next source pointer */ \
	rcx += -2*clm*64+2*clmblkdst; /* Next dest pointer */ \
	rdi += -2*clm*XMM_SCD2;	/* Restore premult/sin/cos data pointer */ \
} \
rsi = (uintptr_t)g->scratch_area;	/* Restore source pointer */ \
rdi += 2*clm*XMM_SCD2;	/* Next sin/cos pointer */ \
end_timer(16); \
\
/* Do inverse FFT levels 1,2,3,4,5 (actually 4.807 levels) */ \
/* On input the 64-byte cache lines hold these data values: */ \
/*	0K	2K	4K	6K */ \
/*	1K	... */ \
/*	8K	... */ \
/*	... */ \
/* On output the 64-byte cache lines hold these data values: */ \
/*	0K	28K	56K	84K */ \
/*	1K	... */ \
/*	... */ \
/*	27K */ \
\
/* Do 16 28_reals_last_unfft operations */ \
/*	distance between fft data elements is 4K */ \
/*	do 2 macros each processing 112 data values */ \
\
start_timer(20); \
set_data_prefetch_ptrs;		/* Init prefetching pointers */ \
loops_init(2*2*clm);		/* 2 iterations of 2*clm */ \
r7_x14cl_28_reals_last_unfft_preload; \
for(unsigned int loopA = 2*2*clm; loopA; ) { \
	xtouch(xptr(rcx+rbp+(clm-1)*128));	/* Preload TLB */ \
	do{ \
		r7_x14cl_28_reals_last_unfft(rsi, 64, 2*clmblkdst, rdi, 2*2*clm*XMM_SCD13); \
		rdi += XMM_SCD13;		/* Next sin/cos pointer */ \
	}while(--loopA & (2*clm-1)); \
	newxloop1quadpre(clm, c1b, c1a, clmblkdst); /* 2*clm iterations */ \
} \
end_timer(20); \
\
/* Normalize these values */ \
\
pass1_normalize(1, yes_fft); \
\
/* Copy data from the scratch area back to the FFT data area */ \
\
copy_scratch_no_pad_data_to_fft(28, clm); \
goto	endlp;			/* Skip FFT code if flag not set */ \
\
/* Fire up auxiliary threads to do pass 1 forward FFT */ \
\
/* pass1: */	pass1_forward_fft_setup(b0b); \
\
/* Copy data from the FFT data area back to the scratch area */ \
\
b0b:	copy_fft_data_to_scratch_no_pad(28, clm); \
\
/* Do FFT levels 1,2,3,4,5 (actually 4.807 levels) */ \
/* */ \
/* On input the 64-byte cache lines hold these data values: */ \
/*	0K	28K	56K	84K */ \
/*	1K	... */ \
/*	... */ \
/*	27K */ \
/* On output the 64-byte cache lines hold these data values: */ \
/*	0K	4K	8K	12K */ \
/*	1K	... */ \
/*	... */ \
/*	3K	... */ \
/*	16K	... */ \
/*	... */ \
\
/* Do 16 28_reals_first_fft operations */ \
/*	distance between fft data elements is 4K */ \
/*	do 4 iterations each processing 56 data values */ \
\
yes_fft:start_timer(22); \
rsi = (uintptr_t)g->scratch_area;	/* Get address of scratch area */ \
rdi = (uintptr_t)g->premult_addr;/* Address of premultiplier/sin/cos table */ \
rdi += 2*clm*XMM_SCD2 + 2*clm*XMM_SCD2; \
set_data_prefetch_ptrs;		/* Init prefetching pointers */ \
pfing(rcx += 8*blkdst);		/* We've already prefetched 8 blocks */ \
loops_init(4*2*clm);		/* 4 iters of 2*clm */ \
r7_x7cl_28_reals_first_fft_scratch_preload; \
for(unsigned int loopA = 4*2*clm; loopA; ) { \
	xtouch(xptr(rcx+rbp+(clm-1)*128));	/* Preload TLB */ \
	do{ \
		r7_x7cl_28_reals_first_fft_scratch(rsi, 64, 4*clmblkdst, rdi); \
		rdi += XMM_SCD13;		/* Next sin/cos pointer */ \
	}while(--loopA & (2*clm-1)); \
	newxloop1quadpre(clm, b1b, b1a, clmblkdst); /* 2*clm iters prefetching 8*clm cache lines */ \
} \
rsi += -4*clmblkdst;	/* Restore source pointer */ \
end_timer(22); \
\
/* Do FFT levels 6,7 */ \
/* */ \
/* On input the 64-byte cache lines hold these data values (+512 in high word): */ \
/*	0K	4K	8K	12K */ \
/*	1K	... */ \
/*	... */ \
/*	3K	... */ \
/*	16K	... */ \
/*	... */ \
/* On output the 64-byte cache lines hold these data values if pass 2 swizzles: */ \
/*	0K	+512	1K	+512	2K	+512	3K	+512 */ \
/*	4K	... */ \
/*	... */ \
/* On output the 64-byte cache lines hold these data values if pass 1 swizzles: */ \
/*	0K	+2K	512	+2K	1K	+2K	1536	+2K */ \
/*	4K	... */ \
/*	... */ \
\
/* Do 1 eight_reals_fft operations */ \
/* Do 1 four_complex_djbfft operations */ \
/*	distance between fft data elements is 1K */ \
/*	do 1 macros each processing 32 data values */ \
/* Do 12 four_complex_djbfft operations */ \
/*	distance between fft data elements is 1K */ \
/*	do 6 macros each processing 32 data values */ \
\
start_timer(26); \
rdi += -4*2*clm*XMM_SCD13-2*clm*XMM_SCD2; /* Complex data sin/cos ptr */ \
rdx = (rdi-2*clm*XMM_SCD2); /* Real data sin/cos ptr */ \
rcx = (uintptr_t)g->data_addr;/* Destination pointer */ \
loops_init_prefetch(7*2*clm, 128, 1, rbp); /* 7 iters of 2*clm */ \
pfing(rbp = (uintptr_t)g->premult_prefetch);/* Load premult/sin/cos prefetch pointer */ \
r4_sg4cl_eight_reals_four_complex_djbfft_preload; \
for(unsigned int loopA = 2*clm; loopA; loopA--) { \
	r4_sg4cl_eight_reals_four_complex_djbfft(rsi, 64, clmblkdst, 2*clmblkdst, rcx, 64, blkdst, 2*blkdst, rdi, rdx); \
	rdi += XMM_SCD2;		/* Next sine/cosine pointer */ \
	rdx += XMM_SCD2;		/* Next sine/cosine pointer */ \
} \
rsi += -2*clm*64+4*clmblkdst; /* Next source pointer */ \
rcx += -2*clm*64+4*blkdst;	/* Next dest pointer */ \
rdi += -2*clm*XMM_SCD2;	/* Restore premult/sine/cosine pointer */ \
\
loops_reset; \
r4_sg4cl_four_complex_djbfft_preload; \
for(unsigned int loopA = 6; loopA; loopA--) { \
	for(unsigned int loopB = 2*clm; loopB; loopB--) { \
		r4_sg4cl_four_complex_djbfft(rsi, 64, clmblkdst, 2*clmblkdst, rcx, 64, blkdst, 2*blkdst, rdi); \
		rdi += XMM_SCD2;		/* Next premult/sin/cos ptr */ \
	} \
	rsi += -2*clm*64+4*clmblkdst; /* Next source pointer */ \
	rcx += -2*clm*64+4*blkdst;	/* Next dest pointer */ \
	rdi += -2*clm*XMM_SCD2;	/* Restore premult/sine/cosine pointer */ \
} \
end_timer(26); \
\
/* Work on next 4*clm sets of 112 values */ \
\
endlp:	xpass1_get_next_block(pass2, c0b, b0b); /* Get next block and jump */ \
/* done: */ xfft_footer_r; \
\
\
\
/* Perform a pass 1 FFT of length 128 doubles. */ \
/* Memory examples are for a 128K FFT with 10 levels done in pass 2. */ \
/* Sin/cos data needed: 16*2*clm*XMM_SCD3 + 4*2*clm*XMM_SCD2 + 4*2*clm*XMM_SCD2 + 2*clm*XMM_SCD2 + 2*clm*XMM_SCD2 = 34*2*clm 64-byte cache lines */ \
/* BUG:	we under-prefetch the sin/cos data by 2*2*clm*64 bytes */ \
\

#define r4_pass1sc128(procname, xpass2name, pass2_levels, clm) /* (fftname, pass2_levels, clmarg) */ \
\
/* set_FFT_constants(fftname, 0, 0, pass2_levels, clmarg); */ \
xfft_header_r(procname, xpass2name, pass2_levels, clm); \
\
/* Do pass 2 (the last 10 FFT levels if this is a 128K FFT) */ \
\
/* pass2: */	start_timer(0); \
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
/* On input the 64-byte cache lines hold these data values if pass 2 swizzles: */ \
/*	0K	+512	1K	+512	2K	+512	3K	+512 */ \
/*	4K */ \
/*	... */ \
/* On input the 64-byte cache lines hold these data values if pass 1 swizzles: */ \
/*	0K	+2K	512	+2K	1K	+2K	1536	+2K */ \
/*	4K */ \
/*	... */ \
/* On output the 64-byte cache lines hold these data values (+512 in high word): */ \
/*	0K	2K	4K	6K */ \
/*	1K	... */ \
/*	8K	... */ \
/*	... */ \
\
/* Do 1 eight_reals_unfft operations */ \
/*	distance between fft data elements is 1K */ \
/*	do 1 macros each processing 16 data values */ \
/* Do 15 four_complex_djbunfft operations */ \
/*	distance between fft data elements is 1K */ \
/*	do 15 macros each processing 16 data values */ \
\
c0b:	start_timer(16); \
rsi = (uintptr_t)g->data_addr;/* Get source address */ \
rcx = (uintptr_t)g->scratch_area;	/* Get address of scratch area */ \
rdx = (uintptr_t)g->premult_addr;/* Real data sin/cos pointer */ \
rdi = (rdx+2*clm*XMM_SCD2); /* Complex data sin/cos ptr */ \
pfing(rbp = (rdx+(34-16)*2*clm*64)); /* Finish prefetch the sin/cos/premult data */ \
loops_init_prefetch_touch(16*2*clm, 64, 1, rbp, 4096, c9a); /* 16*2*clm iterations */ \
xtouch(xptr(rbp));			/* Prefetch TLB */ \
for(unsigned int loopA = 2*clm; loopA; loopA--) { \
	r4_sg2cl_eight_reals_djbunfft(rsi, 64, blkdst, rcx, 64, clmblkdst, rdi, rdx); \
	rdi += XMM_SCD2;		/* Next premult/sin/cos data */ \
	rdx += XMM_SCD2;		/* Next sine/cosine pointer */ \
} \
rsi += -2*clm*64+2*blkdst;	/* Next source pointer */ \
rcx += -2*clm*64+2*clmblkdst; /* Next dest pointer */ \
rdi += -2*clm*XMM_SCD2;	/* Restore premult/sine/cosine pointer */ \
\
loops_reset; \
c9a:	xtouch(xptr(rbp+4096));		/* Prefetch TLB */ \
for(unsigned int loopA = 4*4-1; loopA; ) { \
	do{ \
		for(unsigned int loopC = 2*clm; loopC; loopC--) { \
			r4_sg2cl_four_complex_djbunfft(rsi, 64, blkdst, rcx, 64, clmblkdst, rdi); \
			rdi += XMM_SCD2;		/* Next premult/sin/cos data */ \
		} \
		rsi += -2*clm*64+2*blkdst;	/* Next source pointer */ \
		rcx += -2*clm*64+2*clmblkdst; /* Next dest pointer */ \
		rdi += -2*clm*XMM_SCD2;	/* Restore premult/sin/cos data pointer */ \
	}while(--loopA & 3); \
	rcx += -4*2*clmblkdst+clmblkdst8; \
} \
rsi = (uintptr_t)g->scratch_area;	/* Restore source pointer */ \
rdx = (rdi+2*clm*XMM_SCD2); /* Next real data sin/cos pointer */ \
end_timer(16); \
\
/* Do inverse FFT levels 4,5 */ \
/* On input the 64-byte cache lines hold these data values: */ \
/*	0K	2K	4K	6K */ \
/*	1K	... */ \
/*	8K	... */ \
/*	... */ \
/* On output the 64-byte cache lines hold these data values: */ \
/*	0K	8K	16K	24K */ \
/*	1K	... */ \
/*	... */ \
/*	7K */ \
/*	32K	... */ \
/*	... */ \
\
/* Do 4 eight_reals_unfft operations */ \
/*	distance between fft data elements is 4K */ \
/*	do 2 macros each processing 32 data values */ \
/* Do 12 four_complex_unfft operations */ \
/*	distance between fft data elements is 4K */ \
/*	do 6 macros each processing 32 data values */ \
\
start_timer(17); \
rdi = (rdx+4*2*clm*XMM_SCD2); /* Next complex data sin/cos ptr */ \
set_data_prefetch_ptrs;		/* Init prefetching pointers */ \
loops_init(8*2*clm);		/* 8*2*clm iterations */ \
\
r4_x4cl_eight_reals_unfft_preload; \
for(unsigned int loopA = 2*2*clm; loopA; ) { \
	xtouch(xptr(rcx+rbp+(clm-1)*128));	/* Preload TLB */ \
	do{ \
		r4_x4cl_eight_reals_unfft(rsi, 64, 2*clmblkdst, 4*clmblkdst, rdi, 2*2*clm*XMM_SCD2, rdx, 2*2*clm*XMM_SCD2); \
		rdi += XMM_SCD2;		/* Next sine/cosine pointer */ \
		rdx += XMM_SCD2;		/* Next sine/cosine pointer */ \
	}while(--loopA & (2*clm-1)); \
	newxloop1(clm, c6b, clmblkdst);	/* 2*clm iterations */ \
} \
rsi += -2*clmblkdst+clmblkdst8; /* Next source pointer */ \
rdi += -2*2*clm*XMM_SCD2;	/* Restore sine/cosine pointer */ \
\
loops_reset; \
r4_x4cl_four_complex_djbunfft_preload; \
for(unsigned int loopA = 3; loopA; loopA--) { \
	for(unsigned int loopB = 2*2*clm; loopB; ) { \
		xtouch(xptr(rcx+rbp+(clm-1)*128));	/* Preload TLB */ \
		do{ \
			r4_x4cl_four_complex_djbunfft(rsi, 64, 2*clmblkdst, 4*clmblkdst, rdi, 2*2*clm*XMM_SCD2); \
			rdi += XMM_SCD2;		/* Next clm's sin/cos data */ \
		}while(--loopB & (2*clm-1)); \
		newxloop1(clm, c7b, clmblkdst);	/* 2*clm iterations */ \
	} \
	rsi += -2*clmblkdst+clmblkdst8; /* Next source pointer */ \
	rdi += -2*2*clm*XMM_SCD2;	/* Restore sine/cosine pointer */ \
} \
rsi += -4*clmblkdst8;	/* Restore source pointer */ \
rdi += 4*2*clm*XMM_SCD2;	/* Next sin/cos pointer */ \
end_timer(17); \
\
/* Do inverse FFT levels 1,2,3 */ \
/* On input the 64-byte cache lines hold these data values: */ \
/*	0K	8K	16K	24K */ \
/*	1K	... */ \
/*	... */ \
/*	7K */ \
/*	32K	... */ \
/*	... */ \
/* On output the 64-byte cache lines hold these data values: */ \
/*	0K	32K	64K	96K */ \
/*	1K	... */ \
/*	... */ \
/*	31K */ \
\
/* Do 16 eight_reals_last_unfft operations */ \
/*	distance between fft data elements is 16K */ \
/*	do 8 macros each processing 32 data values */ \
\
start_timer(20); \
loops_init(8*2*clm);		/* 8 iterations of 2*clm */ \
r4_x4cl_eight_reals_last_unfft_preload; \
for(unsigned int loopA = 8*2*clm; loopA; ) { \
	xtouch(xptr(rcx+rbp+(clm-1)*128));	/* Preload TLB */ \
	do{ \
		r4_x4cl_eight_reals_last_unfft(rsi, 64, clmblkdst8, 2*clmblkdst8, rdi, 8*2*clm*XMM_SCD3); \
		rdi += XMM_SCD3;		/* Next sin/cos pointer */ \
	}while(--loopA & (2*clm-1)); \
	newxloop1(clm, c1b, clmblkdst);	/* 2*clm iterations */ \
} \
end_timer(20); \
\
/* Normalize these values */ \
\
pass1_normalize(1, yes_fft); \
\
/* Copy data from the scratch area back to the FFT data area */ \
\
copy_scratch_data_to_fft(32, clm); \
goto	endlp;			/* Skip FFT code if flag not set */ \
\
/* Fire up auxiliary threads to do pass 1 forward FFT */ \
\
/* pass1: */	pass1_forward_fft_setup(b0b); \
\
/* Copy data from the FFT data area back to the scratch area */ \
\
b0b:	copy_fft_data_to_scratch(32, clm); \
\
/* Do FFT levels 1,2,3 */ \
/* */ \
/* On input the 64-byte cache lines hold these data values: */ \
/*	0K	32K	64K	96K */ \
/*	1K	... */ \
/*	... */ \
/*	31K */ \
/* On output the 64-byte cache lines hold these data values: */ \
/*	0K	16K	32K	48K */ \
/*	1K	... */ \
/*	... */ \
/*	15K	... */ \
/*	64K	... */ \
/*	... */ \
\
/* Do 16 eight_reals_first_fft operations */ \
/*	distance between fft data elements is 16K */ \
/*	do 16 iterations each processing 16 data values */ \
\
yes_fft:start_timer(22); \
rsi = (uintptr_t)g->scratch_area;	/* Get address of scratch area */ \
rdi = (uintptr_t)g->premult_addr;/* Address of premultiplier/sin/cos table */ \
rdi += 4*2*clm*XMM_SCD2 + 4*2*clm*XMM_SCD2 + 2*clm*XMM_SCD2 + 2*clm*XMM_SCD2; \
set_data_prefetch_ptrs;		/* Init prefetching pointers */ \
pfing(rcx += 16*blkdst);		/* We've already prefetched 16 blocks */ \
loops_init(16*2*clm);		/* 2 iters of 4 iters of 2*2*clm */ \
r4_x2cl_eight_reals_first_fft_scratch_preload; \
for(unsigned int loopA = 2; loopA; loopA--) { \
	for(unsigned int loopB = 4*2*2*clm; loopB; ) { \
		xtouch(xptr(rcx+rbp+(clm-1)*128));	/* Preload TLB */ \
		do{ \
			do{ \
				r4_x2cl_eight_reals_first_fft_scratch(rsi, 64, 2*clmblkdst8, rdi); \
				rdi += XMM_SCD3;		/* Next sin/cos pointer */ \
			}while(--loopB & (2*clm-1)); \
			newxloop2(clm, b1b, clmblkdst);	/* 2*2*clm iters prefetching 2*clm cache lines */ \
		}while(loopB & (2*2*clm-1)); \
	} \
	rsi += -4*2*clmblkdst+clmblkdst8; \
} \
rsi += -2*clmblkdst8;	/* Restore source pointer */ \
end_timer(22); \
\
/* Do FFT levels 4,5 */ \
/* */ \
/* On input the 64-byte cache lines hold these data values: */ \
/*	0K	16K	32K	48K */ \
/*	1K	... */ \
/*	... */ \
/*	15K	... */ \
/*	64K	... */ \
/*	... */ \
/* On output the 64-byte cache lines hold these data values: */ \
/*	0K	4K	8K	12K */ \
/*	1K	... */ \
/*	... */ \
/*	3K	... */ \
/*	16K	... */ \
/*	... */ \
\
/* Do 4 eight_reals_fft operations */ \
/* Do 4 four_complex_fft operations */ \
/*	distance between fft data elements is 4K */ \
/*	do 4 macros each processing 32 data values */ \
/* Do 8 four_complex_fft operations */ \
/*	distance between fft data elements is 4K */ \
/*	do 4 macros each processing 32 data values */ \
\
start_timer(25); \
rdi += -16*2*clm*XMM_SCD3-4*2*clm*XMM_SCD2; /* Complex data sin/cos ptr */ \
rdx = (rdi-4*2*clm*XMM_SCD2); /* Real data sin/cos ptr */ \
loops_init(8*2*clm);		/* 2 iterations of 4 iters of 2*clm */ \
r4_x4cl_eight_reals_four_complex_djbfft_preload; \
for(unsigned int loopA = 4*2*clm; loopA; ) { \
	xtouch(xptr(rcx+rbp+(clm-1)*128));	/* Preload TLB */ \
	do{ \
		r4_x4cl_eight_reals_four_complex_djbfft(rsi, 64, 4*clmblkdst, clmblkdst8, rdi, rdx); \
		rdi += XMM_SCD2;		/* Next sine/cosine pointer */ \
		rdx += XMM_SCD2;		/* Next sine/cosine pointer */ \
	}while(--loopA & (2*clm-1)); \
	newxloop1(clm, b6b, clmblkdst);	/* 2*clm iterations */ \
} \
rsi += -4*clmblkdst+2*clmblkdst8; \
rdi += -4*2*clm*XMM_SCD2;	/* Restore sin/cos ptr */ \
\
loops_reset; \
r4_x4cl_four_complex_djbfft_preload; \
{ \
	for(unsigned int loopB = 4*2*clm; loopB; ) { \
		xtouch(xptr(rcx+rbp+(clm-1)*128));	/* Preload TLB */ \
		do{ \
			r4_x4cl_four_complex_djbfft(rsi, 64, 4*clmblkdst, clmblkdst8, rdi); \
			rdi += XMM_SCD2;		/* Next sine/cosine pointer */ \
		}while(--loopB & (2*clm-1)); \
		newxloop1(clm, b7b, clmblkdst);	/* 2*clm iterations */ \
	} \
	rsi += -4*clmblkdst+2*clmblkdst8; /* Next source pointer */ \
	rdi += -4*2*clm*XMM_SCD2;	/* Restore sin/cos ptr */ \
} \
rsi += -2*2*clmblkdst8;	/* Restore source pointer */ \
end_timer(25); \
\
/* Do FFT levels 6,7 */ \
/* */ \
/* On input the 64-byte cache lines hold these data values (+512 in high word): */ \
/*	0K	4K	8K	12K */ \
/*	1K	... */ \
/*	... */ \
/*	3K	... */ \
/*	16K	... */ \
/*	... */ \
/* On output the 64-byte cache lines hold these data values if pass 2 swizzles: */ \
/*	0K	+512	1K	+512	2K	+512	3K	+512 */ \
/*	4K	... */ \
/*	... */ \
/* On output the 64-byte cache lines hold these data values if pass 1 swizzles: */ \
/*	0K	+2K	512	+2K	1K	+2K	1536	+2K */ \
/*	4K	... */ \
/*	... */ \
\
/* Do 1 eight_reals_fft operations */ \
/* Do 1 four_complex_djbfft operations */ \
/*	distance between fft data elements is 1K */ \
/*	do 1 macros each processing 32 data values */ \
/* Do 14 four_complex_djbfft operations */ \
/*	distance between fft data elements is 1K */ \
/*	do 7 macros each processing 32 data values */ \
\
start_timer(26); \
rdi += -4*2*clm*XMM_SCD2-2*clm*XMM_SCD2; /* Complex data sin/cos ptr */ \
rdx = (rdi-2*clm*XMM_SCD2); /* Real data sin/cos ptr */ \
rcx = (uintptr_t)g->data_addr;/* Destination pointer */ \
loops_init_prefetch(8*2*clm, 128, 1, rbp); /* 8 iters of 2*clm */ \
pfing(rbp = (uintptr_t)g->premult_prefetch);/* Load premult/sin/cos prefetch pointer */ \
r4_sg4cl_eight_reals_four_complex_djbfft_preload; \
for(unsigned int loopA = 2*clm; loopA; loopA--) { \
	r4_sg4cl_eight_reals_four_complex_djbfft(rsi, 64, clmblkdst, 2*clmblkdst, rcx, 64, blkdst, 2*blkdst, rdi, rdx); \
	rdi += XMM_SCD2;		/* Next sine/cosine pointer */ \
	rdx += XMM_SCD2;		/* Next sine/cosine pointer */ \
} \
rsi += -2*clm*64+4*clmblkdst; /* Next source pointer */ \
rcx += -2*clm*64+4*blkdst;	/* Next dest pointer */ \
rdi += -2*clm*XMM_SCD2;	/* Restore premult/sine/cosine pointer */ \
\
loops_reset; \
r4_sg4cl_four_complex_djbfft_preload; \
for(unsigned int loopA = 4*2-1; loopA; ) { \
	do{ \
		for(unsigned int loopC = 2*clm; loopC; loopC--) { \
			r4_sg4cl_four_complex_djbfft(rsi, 64, clmblkdst, 2*clmblkdst, rcx, 64, blkdst, 2*blkdst, rdi); \
			rdi += XMM_SCD2;		/* Next premult/sin/cos ptr */ \
		} \
		rsi += -2*clm*64+4*clmblkdst; /* Next source pointer */ \
		rcx += -2*clm*64+4*blkdst;	/* Next dest pointer */ \
		rdi += -2*clm*XMM_SCD2;	/* Restore premult/sine/cosine pointer */ \
	}while(--loopA & 1); \
	rsi += -2*4*clmblkdst+clmblkdst8; /* Next source pointer */ \
} \
end_timer(26); \
\
/* Work on next 4*clm sets of 128 values */ \
\
endlp:	xpass1_get_next_block(pass2, c0b, b0b); /* Get next block and jump */ \
/* done: */ xfft_footer_r; \
\
\
\
/* Perform a pass 1 FFT of length 224 doubles. */ \
/* Memory examples are for a 224K FFT with 10 levels done in pass 2. */ \
/* Sin/cos data needed: 8*2*clm*XMM_SCD13 + 2*clm*XMM_SCD4 + 2*clm*XMM_SCD8 = 58*2*clm 64-byte cache lines */ \
/* NOTE:  We under-prefetch sin/cos data by several cache lines */ \
/* NOTE:  We under-prefetch FFT data by 8 blocks */ \
\

#define r4_pass1sc224(procname, xpass2name, pass2_levels, clm) /* (fftname, pass2_levels, clmarg) */ \
\
/* set_FFT_constants(fftname, 0, 0, pass2_levels, clmarg); */ \
xfft_header_r(procname, xpass2name, pass2_levels, clm); \
\
/* Do pass 2 (the last 10 FFT levels if this is a 224K FFT) */ \
\
/* pass2: */	start_timer(0); \
xpass2name(g); \
end_timer(0); \
if(g->ffttype == 1)		/* We're done if FFTing only */ \
	return; /* goto done; */ \
\
/* Do the inverse FFT */ \
\
pass1_inverse_fft_setup(c0b, procname); \
\
/* Do inverse FFT levels 6,7,8 */ \
/* On input the 64-byte cache lines hold these data values if pass 2 swizzles: */ \
/*	0K	+512	1K	+512	2K	+512	3K	+512 */ \
/*	4K */ \
/*	... */ \
/* On input the 64-byte cache lines hold these data values if pass 1 swizzles: */ \
/*	0K	+2K	512	+2K	1K	+2K	1536	+2K */ \
/*	4K */ \
/*	... */ \
/* On output the 64-byte cache lines hold these data values (+512 in high word): */ \
/*	0K	4K	8K	12K */ \
/*	1K	... */ \
/*	... */ \
/*	3K */ \
/*	16K	... */ \
/*	... */ \
\
/* Do 1 sixteen_reals_unfft operations */ \
/*	distance between fft data elements is 1K */ \
/*	do 1 macros each processing 32 data values */ \
/* Do 13 eight_complex_unfft8 operations */ \
/*	distance between fft data elements is 1K */ \
/*	do 13 macros each processing 32 data values */ \
\
c0b:	start_timer(16); \
rsi = (uintptr_t)g->data_addr;/* Get source address */ \
rcx = (uintptr_t)g->scratch_area;	/* Get address of scratch area */ \
rdx = (uintptr_t)g->premult_addr;/* Real data sin/cos ptr */ \
rdi = (rdx+2*clm*XMM_SCD4+128); /* Complex data sin/cos pointer and use -128 */ \
;					/* addressing to generate 1 byte offsets in radix-8 macro */ \
pfing(rbp = (rdx+(58-28)*2*clm*64)); /* Finish prefetching the sin/cos/premult data */ \
loops_init_prefetch(14*2*clm, 128, 1, rbp); /* 14*2*clm iterations */ \
/*bug	xtouchpages rbp, 14*2*clm * 128 / 1 ;; Prefetch TLBs */ \
for(unsigned int loopA = 2*clm; loopA; loopA--) { \
	r8_sg4cl_sixteen_reals_unfft8(rsi, 64, blkdst, 2*blkdst, rcx, 64, clmblkdst, 2*clmblkdst, rdi-128, rdx); \
	rdi += XMM_SCD8;		/* Next premult/sin/cos data */ \
	rdx += XMM_SCD4;		/* Next sine/cosine pointer */ \
} \
rsi += -2*clm*64+4*blkdst;	/* Next source pointer */ \
rcx += -2*clm*64+4*clmblkdst; /* Next dest pointer */ \
rdi += -2*clm*XMM_SCD8;	/* Restore sin/cos pointer */ \
\
loops_reset; \
c9a:	xtouch(xptr(rbp+4096));		/* Prefetch TLB */ \
for(unsigned int loopA = 13; loopA; loopA--) { \
	for(unsigned int loopB = 2*clm; loopB; loopB--) { \
		r8_sg4cl_eight_complex_unfft8(rsi, 64, blkdst, 2*blkdst, rcx, 64, clmblkdst, 2*clmblkdst, rdi-128); \
		rdi += XMM_SCD8;		/* Next premult/sin/cos data */ \
	} \
	rsi += -2*clm*64+4*blkdst;	/* Next source pointer */ \
	rcx += -2*clm*64+4*clmblkdst; /* Next dest pointer */ \
	rdi += -2*clm*XMM_SCD8;	/* Restore premult/sin/cos data pointer */ \
} \
rsi = (uintptr_t)g->scratch_area;	/* Restore source pointer */ \
rdi += -128+2*clm*XMM_SCD8; /* Next sin/cos pointer */ \
end_timer(16); \
\
/* Do inverse FFT levels 1,2,3,4,5 (actually 4.807 levels) */ \
/* On input the 64-byte cache lines hold these data values: */ \
/*	0K	4K	8K	12K */ \
/*	1K	... */ \
/*	... */ \
/*	3K */ \
/*	16K	... */ \
/*	... */ \
/* On output the 64-byte cache lines hold these data values: */ \
/*	0K	56K	112K	168K */ \
/*	1K	... */ \
/*	... */ \
/*	55K */ \
\
/* Do 32 28_reals_last_unfft operations */ \
/*	distance between fft data elements is 8K */ \
/*	do 4 macros each processing 112 data values */ \
\
start_timer(20); \
set_data_prefetch_ptrs;		/* Init prefetching pointers */ \
loops_init(4*2*clm);		/* 4 iterations of 2*clm */ \
r7_x14cl_28_reals_last_unfft_preload; \
for(unsigned int loopA = 4*2*clm; loopA; ) { \
	xtouch(xptr(rcx+rbp+(clm-1)*128));	/* Preload TLB */ \
	do{ \
		r7_x14cl_28_reals_last_unfft(rsi, 64, 4*clmblkdst, rdi, 4*2*clm*XMM_SCD13); \
		rdi += XMM_SCD13;		/* Next sin/cos pointer */ \
	}while(--loopA & (2*clm-1)); \
	newxloop1quadpre(clm, c1b, c1a, clmblkdst); /* 2*clm iterations */ \
} \
end_timer(20); \
\
/* Normalize these values */ \
\
pass1_normalize(1, yes_fft); \
\
/* Copy data from the scratch area back to the FFT data area */ \
\
copy_scratch_no_pad_data_to_fft(56, clm); \
goto	endlp;			/* Skip FFT code if flag not set */ \
\
/* Fire up auxiliary threads to do pass 1 forward FFT */ \
\
/* pass1: */	pass1_forward_fft_setup(b0b); \
\
/* Copy data from the FFT data area back to the scratch area */ \
\
b0b:	copy_fft_data_to_scratch_no_pad(56, clm); \
\
/* Do FFT levels 1,2,3,4,5 (actually 4.807 levels) */ \
/* */ \
/* On input the 64-byte cache lines hold these data values: */ \
/*	0K	56K	112K	168K */ \
/*	1K	... */ \
/*	... */ \
/*	55K */ \
/* On output the 64-byte cache lines hold these data values: */ \
/*	0K	8K	16K	24K */ \
/*	1K	... */ \
/*	... */ \
/*	7K	... */ \
/*	32K	... */ \
/*	... */ \
\
/* Do 32 28_reals_first_fft operations */ \
/*	distance between fft data elements is 8K */ \
/*	do 8 iterations each processing 56 data values */ \
\
yes_fft:start_timer(22); \
rsi = (uintptr_t)g->scratch_area;	/* Get address of scratch area */ \
rdi = (uintptr_t)g->premult_addr;/* Address of premultiplier/sin/cos table */ \
rdi += 2*clm*XMM_SCD4 + 2*clm*XMM_SCD8; \
set_data_prefetch_ptrs;		/* Init prefetching pointers */ \
pfing(rcx += 16*blkdst);		/* We've already prefetched 16 blocks */ \
loops_init(8*2*clm);		/* 8 iters of 2*clm */ \
r7_x7cl_28_reals_first_fft_scratch_preload; \
for(unsigned int loopA = 8*2*clm; loopA; ) { \
	xtouch(xptr(rcx+rbp+(clm-1)*128));	/* Preload TLB */ \
	do{ \
		r7_x7cl_28_reals_first_fft_scratch(rsi, 64, 8*clmblkdst, rdi); \
		rdi += XMM_SCD13;		/* Next sin/cos pointer */ \
	}while(--loopA & (2*clm-1)); \
	newxloop1quadpre(clm, b1b, b1a, clmblkdst); /* 2*clm iters prefetching 8*clm cache lines */ \
} \
rsi += -8*clmblkdst;	/* Restore source pointer */ \
end_timer(22); \
\
/* Do FFT levels 6,7,8 */ \
/* */ \
/* On input the 64-byte cache lines hold these data values (+512 in high word): */ \
/*	0K	8K	16K	24K */ \
/*	1K	... */ \
/*	... */ \
/*	7K	... */ \
/*	32K	... */ \
/*	... */ \
/* On output the 64-byte cache lines hold these data values if pass 2 swizzles: */ \
/*	0K	+512	1K	+512	2K	+512	3K	+512 */ \
/*	4K	... */ \
/*	... */ \
/* On output the 64-byte cache lines hold these data values if pass 1 swizzles: */ \
/*	0K	+2K	512	+2K	1K	+2K	1536	+2K */ \
/*	4K	... */ \
/*	... */ \
\
/* Do 1 sixteen_reals_fft operations */ \
/* Do 1 eight_complex_fft8 operations */ \
/*	distance between fft data elements is 1K */ \
/*	do 1 macros each processing 64 data values */ \
/* Do 12 eight_complex_fft8 operations */ \
/*	distance between fft data elements is 1K */ \
/*	do 6 macros each processing 64 data values */ \
\
start_timer(26); \
rdi += -8*2*clm*XMM_SCD13-2*clm*XMM_SCD8+128; /* Complex data sin/cos ptr and use -128 */ \
;					/* addressing to generate 1 byte offsets in radix-8 macro */ \
rdx = (rdi-128-2*clm*XMM_SCD4); /* Real data sin/cos ptr */ \
rcx = (uintptr_t)g->data_addr;/* Destination pointer */ \
loops_init_prefetch(7*2*clm, 128, 1, rbp); /* 7 iters of 2*clm */ \
pfing(rbp = (uintptr_t)g->premult_prefetch);/* Load premult/sin/cos prefetch pointer */ \
r8_sg8cl_sixteen_reals_eight_complex_fft8_preload; \
for(unsigned int loopA = 2*clm; loopA; loopA--) { \
	r8_sg8cl_sixteen_reals_eight_complex_fft8(rsi, 64, clmblkdst, 2*clmblkdst, 4*clmblkdst, rcx, 64, blkdst, 2*blkdst, 4*blkdst, rdi-128, rdx); \
	rdi += XMM_SCD8;		/* Next sine/cosine pointer */ \
	rdx += XMM_SCD4;		/* Next sine/cosine pointer */ \
} \
rsi += -2*clm*64+8*clmblkdst; /* Next source pointer */ \
rcx += -2*clm*64+8*blkdst;	/* Next dest pointer */ \
rdi += -2*clm*XMM_SCD8;	/* Restore premult/sine/cosine pointer */ \
\
loops_reset; \
for(unsigned int loopA = 6; loopA; loopA--) { \
	for(unsigned int loopB = 2*clm; loopB; loopB--) { \
		r8_sg8cl_eight_complex_fft8(rsi, 64, clmblkdst, 2*clmblkdst, 4*clmblkdst, rcx, 64, blkdst, 2*blkdst, 4*blkdst, rdi-128); \
		rdi += XMM_SCD8;		/* Next premult/sin/cos ptr */ \
	} \
	rsi += -2*clm*64+8*clmblkdst; /* Next source pointer */ \
	rcx += -2*clm*64+8*blkdst;	/* Next dest pointer */ \
	rdi += -2*clm*XMM_SCD8;	/* Restore premult/sine/cosine pointer */ \
} \
end_timer(26); \
\
/* Work on next 4*clm sets of 224 values */ \
\
endlp:	xpass1_get_next_block(pass2, c0b, b0b); /* Get next block and jump */ \
/* done: */ xfft_footer_r; \
\
\
\
/* Perform a pass 1 FFT of length 256 doubles. */ \
/* Memory examples are for a 256K FFT with 10 levels done in pass 2. */ \
/* Sin/cos data needed: 32*2*clm*XMM_SCD3 + 8*2*clm*XMM_SCD2 + 8*2*clm*XMM_SCD2 + 2*clm*XMM_SCD4 + 2*clm*XMM_SCD8 = 68*2*clm 64-byte cache lines */ \
/* BUG:	we under-prefetch the sin/cos data by 20*2*clm*64 bytes */ \
\

#define r4_pass1sc256(procname, xpass2name, pass2_levels, clm) /* (fftname, pass2_levels, clmarg) */ \
\
/* set_FFT_constants(fftname, 0, 0, pass2_levels, clmarg); */ \
xfft_header_r(procname, xpass2name, pass2_levels, clm); \
\
/* Do pass 2 (the last 10 FFT levels if this is a 256K FFT) */ \
\
/* pass2: */	start_timer(0); \
xpass2name(g); \
end_timer(0); \
if(g->ffttype == 1)		/* We're done if FFTing only */ \
	return; /* goto done; */ \
\
/* Do the inverse FFT */ \
\
pass1_inverse_fft_setup(c0b, procname); \
\
/* Do inverse FFT levels 6,7,8 */ \
/* On input the 64-byte cache lines hold these data values if pass 2 swizzles: */ \
/*	0K	+512	1K	+512	2K	+512	3K	+512 */ \
/*	4K */ \
/*	... */ \
/* On input the 64-byte cache lines hold these data values if pass 1 swizzles: */ \
/*	0K	+2K	512	+2K	1K	+2K	1536	+2K */ \
/*	4K */ \
/*	... */ \
/* On output the 64-byte cache lines hold these data values (+512 in high word): */ \
/*	0K	4K	8K	12K */ \
/*	1K	... */ \
/*	... */ \
/*	3K */ \
/*	16K	... */ \
/*	... */ \
\
/* Do 1 sixteen_reals_unfft operations */ \
/*	distance between fft data elements is 1K */ \
/*	do 1 macros each processing 32 data values */ \
/* Do 15 eight_complex_unfft8 operations */ \
/*	distance between fft data elements is 1K */ \
/*	do 15 macros each processing 32 data values */ \
\
c0b:	start_timer(16); \
rsi = (uintptr_t)g->data_addr;/* Get source address */ \
rcx = (uintptr_t)g->scratch_area;	/* Get address of scratch area */ \
rdx = (uintptr_t)g->premult_addr;/* Real data sin/cos ptr */ \
rdi = (rdx+2*clm*XMM_SCD4+128); /* Complex data sin/cos pointer and use -128 */ \
;					/* addressing to generate 1 byte offsets in radix-8 macro */ \
pfing(rbp = (rdx+(68-32)*2*clm*64)); /* Finish prefetching the sin/cos/premult data */ \
loops_init_prefetch(16*2*clm, 128, 1, rbp); /* 16*2*clm iterations */ \
/*bug	xtouchpages rbp, 16*2*clm * 128 / 1 ;; Prefetch TLBs */ \
for(unsigned int loopA = 2*clm; loopA; loopA--) { \
	r8_sg4cl_sixteen_reals_unfft8(rsi, 64, blkdst, 2*blkdst, rcx, 64, clmblkdst, 2*clmblkdst, rdi-128, rdx); \
	rdi += XMM_SCD8;		/* Next premult/sin/cos data */ \
	rdx += XMM_SCD4;		/* Next sine/cosine pointer */ \
} \
rsi += -2*clm*64+4*blkdst;	/* Next source pointer */ \
rcx += -2*clm*64+4*clmblkdst; /* Next dest pointer */ \
rdi += -2*clm*XMM_SCD8;	/* Restore sin/cos pointer */ \
\
loops_reset; \
c9a:	xtouch(xptr(rbp+4096));		/* Prefetch TLB */ \
for(unsigned int loopA = 8*2-1; loopA; ) { \
	do{ \
		for(unsigned int loopC = 2*clm; loopC; loopC--) { \
			r8_sg4cl_eight_complex_unfft8(rsi, 64, blkdst, 2*blkdst, rcx, 64, clmblkdst, 2*clmblkdst, rdi-128); \
			rdi += XMM_SCD8;		/* Next premult/sin/cos data */ \
		} \
		rsi += -2*clm*64+4*blkdst;	/* Next source pointer */ \
		rcx += -2*clm*64+4*clmblkdst; /* Next dest pointer */ \
		rdi += -2*clm*XMM_SCD8;	/* Restore premult/sin/cos data pointer */ \
	}while(--loopA & 1); \
	rcx += -2*4*clmblkdst+clmblkdst8; \
} \
rsi = (uintptr_t)g->scratch_area;	/* Restore source pointer */ \
rdx = (rdi-128+2*clm*XMM_SCD8); /* Next real data sin/cos pointer */ \
end_timer(16); \
\
/* Do inverse FFT levels 4,5 */ \
/* On input the 64-byte cache lines hold these data values: */ \
/*	0K	4K	8K	12K */ \
/*	1K	... */ \
/*	... */ \
/*	3K */ \
/*	16K	... */ \
/*	... */ \
/* On output the 64-byte cache lines hold these data values: */ \
/*	0K	16K	32K	48K */ \
/*	1K	... */ \
/*	... */ \
/*	15K */ \
/*	64K	... */ \
/*	... */ \
\
/* Do 8 eight_reals_unfft operations */ \
/*	distance between fft data elements is 8K */ \
/*	do 4 macros each processing 32 data values */ \
/* Do 24 four_complex_unfft operations */ \
/*	distance between fft data elements is 8K */ \
/*	do 12 macros each processing 32 data values */ \
\
start_timer(17); \
rdi = (rdx+8*2*clm*XMM_SCD2); /* Complex data sin/cos ptr */ \
set_data_prefetch_ptrs;		/* Init prefetching pointers */ \
loops_init(16*2*clm);		/* 16 iterations of 2*clm */ \
r4_x4cl_eight_reals_unfft_preload; \
for(unsigned int loopA = 4*2*clm; loopA; ) { \
	xtouch(xptr(rcx+rbp+(clm-1)*128));	/* Preload TLB */ \
	do{ \
		r4_x4cl_eight_reals_unfft(rsi, 64, 4*clmblkdst, clmblkdst8, rdi, 4*2*clm*XMM_SCD2, rdx, 4*2*clm*XMM_SCD2); \
		rdi += XMM_SCD2;		/* Next sine/cosine pointer */ \
		rdx += XMM_SCD2;		/* Next sine/cosine pointer */ \
	}while(--loopA & (2*clm-1)); \
	newxloop1(clm, c6b, clmblkdst);	/* 2*clm iterations */ \
} \
rsi += -4*clmblkdst+2*clmblkdst8; /* Next source pointer */ \
rdi += -4*2*clm*XMM_SCD2;	/* Restore sine/cosine pointer */ \
\
loops_reset; \
r4_x4cl_four_complex_djbunfft_preload; \
for(unsigned int loopA = 3; loopA; loopA--) { \
	for(unsigned int loopB = 4*2*clm; loopB; ) { \
		xtouch(xptr(rcx+rbp+(clm-1)*128));	/* Preload TLB */ \
		do{ \
			r4_x4cl_four_complex_djbunfft(rsi, 64, 4*clmblkdst, clmblkdst8, rdi, 4*2*clm*XMM_SCD2); \
			rdi += XMM_SCD2;		/* Next clm's sin/cos data */ \
		}while(--loopB & (2*clm-1)); \
		newxloop1(clm, c7b, clmblkdst);	/* 2*clm iterations */ \
	} \
	rsi += -4*clmblkdst+2*clmblkdst8; /* Next source pointer */ \
	rdi += -4*2*clm*XMM_SCD2;	/* Restore sine/cosine pointer */ \
} \
rsi += -4*2*clmblkdst8;	/* Restore source pointer */ \
rdi += 8*2*clm*XMM_SCD2;	/* Next sin/cos ptr */ \
end_timer(17); \
\
/* Do inverse FFT levels 1,2,3 */ \
/* On input the 64-byte cache lines hold these data values: */ \
/*	0K	16K	32K	48K */ \
/*	1K	... */ \
/*	... */ \
/*	15K */ \
/*	64K	... */ \
/*	... */ \
/* On output the 64-byte cache lines hold these data values: */ \
/*	0K	64K	128K	192K */ \
/*	1K	... */ \
/*	... */ \
/*	63K */ \
\
/* Do 32 eight_reals_last_unfft operations */ \
/*	distance between fft data elements is 32K */ \
/*	do 16 macros each processing 32 data values */ \
\
start_timer(20); \
loops_init(16*2*clm);		/* 16 iterations of 2*clm */ \
r4_x4cl_eight_reals_last_unfft_preload; \
for(unsigned int loopA = 2; loopA; loopA--) { \
	for(unsigned int loopB = 8*2*clm; loopB; ) { \
		xtouch(xptr(rcx+rbp+(clm-1)*128));	/* Preload TLB */ \
		do{ \
			r4_x4cl_eight_reals_last_unfft(rsi, 64, 2*clmblkdst8, 4*clmblkdst8, rdi, 16*2*clm*XMM_SCD3); \
			rdi += XMM_SCD3;		/* Next sin/cos ptr */ \
		}while(--loopB & (2*clm-1)); \
		newxloop1(clm, c1b, clmblkdst);	/* 2*clm iterations */ \
	} \
	rsi += -8*clmblkdst+clmblkdst8; /* Next source pointer */ \
} \
end_timer(20); \
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
/*	0K	64K	128K	192K */ \
/*	1K	... */ \
/*	... */ \
/*	63K */ \
/* On output the 64-byte cache lines hold these data values: */ \
/*	0K	32K	64K	96K */ \
/*	1K	... */ \
/*	... */ \
/*	31K	... */ \
/*	128K	... */ \
/*	... */ \
\
/* Do 32 eight_reals_first_fft operations */ \
/*	distance between fft data elements is 32K */ \
/*	do 32 iterations each processing 16 data values */ \
\
yes_fft:start_timer(22); \
rsi = (uintptr_t)g->scratch_area;	/* Get address of scratch area */ \
rdi = (uintptr_t)g->premult_addr;/* Address of premultiplier/sin/cos table */ \
rdi += 8*2*clm*XMM_SCD2 + 8*2*clm*XMM_SCD2 + 2*clm*XMM_SCD4 + 2*clm*XMM_SCD8; \
set_data_prefetch_ptrs;		/* Init prefetching pointers */ \
pfing(rcx += 32*blkdst);		/* We've already prefetched 32 blocks */ \
loops_init(32*2*clm);		/* 4 iters of 4 iters of 2*2*clm */ \
r4_x2cl_eight_reals_first_fft_scratch_preload; \
for(unsigned int loopA = 4; loopA; loopA--) { \
	for(unsigned int loopB = 4*2*2*clm; loopB; ) { \
		xtouch(xptr(rcx+rbp+(clm-1)*128));	/* Preload TLB */ \
		do{ \
			do{ \
				r4_x2cl_eight_reals_first_fft_scratch(rsi, 64, 4*clmblkdst8, rdi); \
				rdi += XMM_SCD3;		/* Next sin/cos ptr */ \
			}while(--loopB & (2*clm-1)); \
			newxloop2(clm, b1b, clmblkdst);	/* 2*2*clm iters prefetching 2*clm cache lines */ \
		}while(loopB & (2*2*clm-1)); \
	} \
	rsi += -4*2*clmblkdst+clmblkdst8; \
} \
rsi += -4*clmblkdst8;	/* Restore source pointer */ \
end_timer(22); \
\
/* Do FFT levels 4,5 */ \
/* */ \
/* On input the 64-byte cache lines hold these data values: */ \
/*	0K	32K	64K	96K */ \
/*	1K	... */ \
/*	... */ \
/*	31K	... */ \
/*	128K	... */ \
/*	... */ \
/* On output the 64-byte cache lines hold these data values: */ \
/*	0K	8K	16K	24K */ \
/*	1K	... */ \
/*	... */ \
/*	7K	... */ \
/*	32K	... */ \
/*	... */ \
\
/* Do 2 eight_reals_fft operations */ \
/* Do 2 four_complex_fft operations */ \
/*	distance between fft data elements is 8K */ \
/*	do 2 macros each processing 32 data values */ \
/* Do 28 four_complex_fft operations */ \
/*	distance between fft data elements is 8K */ \
/*	do 14 macros each processing 32 data values */ \
/*  BUG in comment, numbers in comment may be wrong */ \
start_timer(25); \
rdi += -32*2*clm*XMM_SCD3 - 8*2*clm*XMM_SCD2; /* Complex data sin/cos ptr */ \
rdx = (rdi-8*2*clm*XMM_SCD2); /* Real data sin/cos ptr */ \
loops_init(16*2*clm);		/* 2 iterations of 8 iters of 2*clm */ \
r4_x4cl_eight_reals_four_complex_djbfft_preload; \
for(unsigned int loopA = 8*2*clm; loopA; ) { \
	xtouch(xptr(rcx+rbp+(clm-1)*128));	/* Preload TLB */ \
	do{ \
		r4_x4cl_eight_reals_four_complex_djbfft(rsi, 64, clmblkdst8, 2*clmblkdst8, rdi, rdx); \
		rdi += XMM_SCD2;		/* Next sine/cosine pointer */ \
		rdx += XMM_SCD2;		/* Next sine/cosine pointer */ \
	}while(--loopA & (2*clm-1)); \
	newxloop1(clm, b6b, clmblkdst);	/* 2*clm iterations */ \
} \
rsi += -8*clmblkdst+4*clmblkdst8; \
rdi += -8*2*clm*XMM_SCD2;	/* Restore sin/cos ptr */ \
\
loops_reset; \
r4_x4cl_four_complex_djbfft_preload; \
{ \
	for(unsigned int loopB = 8*2*clm; loopB; ) { \
		xtouch(xptr(rcx+rbp+(clm-1)*128));	/* Preload TLB */ \
		do{ \
			r4_x4cl_four_complex_djbfft(rsi, 64, clmblkdst8, 2*clmblkdst8, rdi); \
			rdi += XMM_SCD2;		/* Next sine/cosine pointer */ \
		}while(--loopB & (2*clm-1)); \
		newxloop1(clm, b7b, clmblkdst);	/* 2*clm iterations */ \
	} \
	rsi += -8*clmblkdst+4*clmblkdst8; /* Next source pointer */ \
	rdi += -8*2*clm*XMM_SCD2;	/* Restore sin/cos ptr */ \
} \
rsi += -2*4*clmblkdst8;	/* Restore source pointer */ \
end_timer(25); \
\
/* Do FFT levels 6,7,8 */ \
/* */ \
/* On input the 64-byte cache lines hold these data values (+512 in high word): */ \
/*	0K	8K	16K	24K */ \
/*	1K	... */ \
/*	... */ \
/*	7K	... */ \
/*	32K	... */ \
/*	... */ \
/* On output the 64-byte cache lines hold these data values if pass 2 swizzles: */ \
/*	0K	+512	1K	+512	2K	+512	3K	+512 */ \
/*	4K	... */ \
/*	... */ \
/* On output the 64-byte cache lines hold these data values if pass 1 swizzles: */ \
/*	0K	+2K	512	+2K	1K	+2K	1536	+2K */ \
/*	4K	... */ \
/*	... */ \
\
/* Do 1 sixteen_reals_fft operations */ \
/* Do 1 eight_complex_fft8 operations */ \
/*	distance between fft data elements is 1K */ \
/*	do 1 macros each processing 64 data values */ \
/* Do 14 eight_complex_fft8 operations */ \
/*	distance between fft data elements is 1K */ \
/*	do 7 macros each processing 64 data values */ \
\
start_timer(26); \
rdi += -8*2*clm*XMM_SCD2-2*clm*XMM_SCD8+128; /* Complex data sin/cos ptr and use -128 */ \
;					/* addressing to generate 1 byte offsets in radix-8 macro */ \
rdx = (rdi-128-2*clm*XMM_SCD4); /* Real data sin/cos ptr */ \
rcx = (uintptr_t)g->data_addr;/* Destination pointer */ \
loops_init_prefetch(8*2*clm, 128, 1, rbp); /* 8 iters of 2*clm */ \
pfing(rbp = (uintptr_t)g->premult_prefetch);/* Load premult/sin/cos prefetch pointer */ \
r8_sg8cl_sixteen_reals_eight_complex_fft8_preload; \
for(unsigned int loopA = 2*clm; loopA; loopA--) { \
	r8_sg8cl_sixteen_reals_eight_complex_fft8(rsi, 64, clmblkdst, 2*clmblkdst, 4*clmblkdst, rcx, 64, blkdst, 2*blkdst, 4*blkdst, rdi-128, rdx); \
	rdi += XMM_SCD8;		/* Next sine/cosine pointer */ \
	rdx += XMM_SCD4;		/* Next sine/cosine pointer */ \
} \
rsi += -2*clm*64+clmblkdst8; /* Next source pointer */ \
rcx += -2*clm*64+8*blkdst;	/* Next dest pointer */ \
rdi += -2*clm*XMM_SCD8;	/* Restore premult/sine/cosine pointer */ \
\
loops_reset; \
for(unsigned int loopA = 7; loopA; loopA--) { \
	for(unsigned int loopB = 2*clm; loopB; loopB--) { \
		r8_sg8cl_eight_complex_fft8(rsi, 64, clmblkdst, 2*clmblkdst, 4*clmblkdst, rcx, 64, blkdst, 2*blkdst, 4*blkdst, rdi-128); \
		rdi += XMM_SCD8;		/* Next premult/sin/cos ptr */ \
	} \
	rsi += -2*clm*64+clmblkdst8; /* Next source pointer */ \
	rcx += -2*clm*64+8*blkdst;	/* Next dest pointer */ \
	rdi += -2*clm*XMM_SCD8;	/* Restore premult/sine/cosine pointer */ \
} \
end_timer(26); \
\
/* Work on next 4*clm sets of 256 values */ \
\
endlp:	xpass1_get_next_block(pass2, c0b, b0b); /* Get next block and jump */ \
/* done: */ xfft_footer_r; \
\
\
\
/* Perform a pass 1 negacyclic FFT of length 32 doubles. */ \
/* Memory examples are for a 32K FFT with 10 levels done in pass 2. */ \
/* Sin/cos data needed: 4*2*clm*XMM_SCD4 + 2*clm*XMM_SCD2 = 9*2*clm 64-byte cache lines */ \
/* BUG:	we under-prefetch the sin/cos data by 1*2*clm 64-byte cache lines */ \
\

#define r4_pass1sc32ac(procname, xpass2name, pass2_levels, clm) /* (fftname, pass2_levels, clmarg) */ \
\
/* set_FFT_constants(fftname, 0, 1, pass2_levels, clmarg); */ \
xfft_header_r(procname, xpass2name, pass2_levels, clm); \
\
/* Do pass 2 (the last 10 FFT levels if this is a 32K FFT) */ \
\
/* pass2: */	start_timer(0); \
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
/* On input the 64-byte cache lines hold these data values if pass 2 swizzles: */ \
/*	0K	+512	1K	+512	2K	+512	3K	+512 */ \
/*	4K */ \
/*	... */ \
/* On input the 64-byte cache lines hold these data values if pass 1 swizzles: */ \
/*	0K	+2K	512	+2K	1K	+2K	1536	+2K */ \
/*	4K */ \
/*	... */ \
/* On output the 64-byte cache lines hold these data values (+512 in high word): */ \
/*	0K	2K	4K	6K */ \
/*	1K	... */ \
/*	8K	... */ \
/*	... */ \
\
/* Do 4 four_complex_djbunfft operations */ \
/*	distance between fft data elements is 1K */ \
/*	do 4 macros each processing 16 data values */ \
\
c0b:	start_timer(16); \
rsi = (uintptr_t)g->data_addr;/* Get source address */ \
rcx = (uintptr_t)g->scratch_area;	/* Get address of scratch area */ \
rdi = (uintptr_t)g->premult_addr;/* Load premult/sin/cos pointer */ \
pfing(rbp = (rdi+(9-4)*2*clm*64)); /* Finish prefetching the sin/cos/premult data */ \
loops_init_prefetch(4*2*clm, 64, 1, rbp); /* 4*2*clm iterations */ \
/*bug	xtouchpages rbp, 4*2*clm * 64 / 1 ;; Prefetch TLBs */ \
for(unsigned int loopA = 4; loopA; loopA--) { \
	for(unsigned int loopB = 2*clm; loopB; loopB--) { \
		r4_sg2cl_four_complex_djbunfft(rsi, 64, blkdst, rcx, 64, clmblkdst, rdi); \
		rdi += XMM_SCD2;		/* Next premult/sin/cos data */ \
	} \
	rsi += -2*clm*64+2*blkdst;	/* Next source pointer */ \
	rcx += -2*clm*64+2*clmblkdst; /* Next dest pointer */ \
	rdi += -2*clm*XMM_SCD2;	/* Restore premult/sin/cos data pointer */ \
} \
rsi = (uintptr_t)g->scratch_area;	/* Restore source pointer */ \
rdi += 2*clm*XMM_SCD2;	/* Next premult/sin/cos data pointer */ \
end_timer(16); \
\
/* Do inverse FFT levels 1,2,3 */ \
/* On input the 64-byte cache lines hold these data values: */ \
/*	0K	2K	4K	6K */ \
/*	1K	... */ \
/*	8K	... */ \
/*	... */ \
/* On output the 64-byte cache lines hold these data values: */ \
/*	0K	8K	16K	24K */ \
/*	1K	... */ \
/*	... */ \
/*	7K */ \
\
/* Do 4 four_complex_unfft operations */ \
/*	distance between fft data elements is 4K */ \
/*	do 2 macros each processing 32 data values */ \
\
start_timer(20); \
set_data_prefetch_ptrs;		/* Init prefetching pointers */ \
loops_init(2*2*clm);		/* 2 iterations of 2*clm */ \
for(unsigned int loopA = 2*2*clm; loopA; ) { \
	xtouch(xptr(rcx+rbp+(clm-1)*128));	/* Preload TLB */ \
	do{ \
		r4_x4cl_four_complex_last_unfft4(rsi, 64, 2*clmblkdst, 4*clmblkdst, rdi, 2*2*clm*XMM_SCD4); \
		rdi += XMM_SCD4;		/* Next sin/cos/premult ptr */ \
	}while(--loopA & (2*clm-1)); \
	newxloop1doublepre(clm, c1b, c1a, clmblkdst);	/* 2*clm iterations */ \
} \
end_timer(20); \
\
/* Normalize these values */ \
\
pass1_normalize(1, yes_fft); \
\
/* Copy data from the scratch area back to the FFT data area */ \
\
copy_scratch_data_to_fft(8, clm); \
goto	endlp;			/* Skip FFT code if flag not set */ \
\
/* Fire up auxiliary threads to do pass 1 forward FFT */ \
\
/* pass1: */	pass1_forward_fft_setup(b0b); \
\
/* Copy data from the FFT data area back to the scratch area */ \
\
b0b:	copy_fft_data_to_scratch(8, clm); \
\
/* Do FFT levels 1,2,3 */ \
/* */ \
/* On input the 64-byte cache lines hold these data values: */ \
/*	0K	8K	16K	24K */ \
/*	1K	... */ \
/*	... */ \
/*	7K */ \
/* On output the 64-byte cache lines hold these data values: */ \
/*	0K	4K	8K	12K */ \
/*	1K	... */ \
/*	... */ \
/*	3K	... */ \
/*	16K	... */ \
/*	... */ \
\
/* Do 4 four_complex_first_fft operations */ \
/*	distance between fft data elements is 4K */ \
/*	do 4 iterations each processing 16 data values */ \
\
yes_fft:start_timer(22); \
rsi = (uintptr_t)g->scratch_area;	/* Get address of scratch area */ \
rdi = (uintptr_t)g->premult_addr;/* Address of premultiplier/sin/cos table */ \
rdi += 2*clm*XMM_SCD2; \
set_data_prefetch_ptrs;		/* Init prefetching pointers */ \
pfing(rcx += 4*blkdst);		/* We've already prefetched 16 blocks */ \
loops_init(4*2*clm);		/* 4 iters of 2*clm */ \
for(unsigned int loopA = 4*2*clm; loopA; ) { \
	xtouch(xptr(rcx+rbp+(clm-1)*128));	/* Preload TLB */ \
	do{ \
		r4_x2cl_four_complex_first_fft4_scratch(rsi, 64, 4*clmblkdst, rdi); \
		rdi += XMM_SCD4;		/* Next sin/cos/premult ptr */ \
	}while(--loopA & (2*clm-1)); \
	newxloop1(clm, b1b, clmblkdst);	/* 2*clm iters */ \
} \
rsi += -4*clmblkdst;	/* Restore source pointer */ \
end_timer(22); \
\
/* Do FFT levels 4,5 */ \
/* */ \
/* On input the 64-byte cache lines hold these data values (+512 in high word): */ \
/*	0K	4K	8K	12K */ \
/*	1K	... */ \
/*	... */ \
/*	3K	... */ \
/*	16K	... */ \
/*	... */ \
/* On output the 64-byte cache lines hold these data values if pass 2 swizzles: */ \
/*	0K	+512	1K	+512	2K	+512	3K	+512 */ \
/*	4K	... */ \
/*	... */ \
/* On output the 64-byte cache lines hold these data values if pass 1 swizzles: */ \
/*	0K	+2K	512	+2K	1K	+2K	1536	+2K */ \
/*	4K	... */ \
/*	... */ \
\
/* Do 4 four_complex_fft operations */ \
/*	distance between fft data elements is 1K */ \
/*	do 2 macros each processing 32 data values */ \
\
start_timer(26); \
rdi += -4*2*clm*XMM_SCD4-2*clm*XMM_SCD2; /* Next premult/sin/cos ptr */ \
rcx = (uintptr_t)g->data_addr;/* Destination pointer */ \
loops_init_prefetch(2*2*clm, 128, 1, rbp); /* 2 iters of 2*clm */ \
pfing(rbp = (uintptr_t)g->premult_prefetch);/* Load premult/sin/cos prefetch pointer */ \
r4_sg4cl_four_complex_djbfft_preload; \
for(unsigned int loopA = 2; loopA; loopA--) { \
	for(unsigned int loopB = 2*clm; loopB; loopB--) { \
		r4_sg4cl_four_complex_djbfft(rsi, 64, clmblkdst, 2*clmblkdst, rcx, 64, blkdst, 2*blkdst, rdi); \
		rdi += XMM_SCD2;		/* Next premult/sin/cos ptr */ \
	} \
	rsi += -2*clm*64+4*clmblkdst; /* Next source pointer */ \
	rcx += -2*clm*64+4*blkdst;	/* Next dest pointer */ \
	rdi += -2*clm*XMM_SCD2;	/* Restore premult/sine/cosine pointer */ \
} \
end_timer(26); \
\
/* Work on next 4*clm sets of 32 values */ \
\
endlp:	xpass1_get_next_block(pass2, c0b, b0b); /* Get next block and jump */ \
/* done: */ xfft_footer_r; \
\
\
\
/* Perform a pass 1 negacyclic FFT of length 64 doubles. */ \
/* Memory examples are for a 64K FFT with 10 levels done in pass 2. */ \
/* Sin/cos data needed: 8*2*clm*XMM_SCD4 + 2*clm*XMM_SCD8 = 20*2*clm 64-byte cache lines */ \
/* BUG - we substantially underprefetch the sin/cos data */ \
\

#define r4_pass1sc64ac(procname, xpass2name, pass2_levels, clm) /* (fftname, pass2_levels, clmarg) */ \
\
/* set_FFT_constants(fftname, 0, 1, pass2_levels, clmarg); */ \
xfft_header_r(procname, xpass2name, pass2_levels, clm); \
\
/* Do pass 2 (the last 10 FFT levels if this is a 64K FFT) */ \
\
/* pass2: */	start_timer(0); \
xpass2name(g); \
end_timer(0); \
if(g->ffttype == 1)		/* We're done if FFTing only */ \
	return; /* goto done; */ \
\
/* Do the inverse FFT */ \
\
pass1_inverse_fft_setup(c0b, procname); \
\
/* Do inverse FFT levels 4,5,6 */ \
/* On input the 64-byte cache lines hold these data values if pass 2 swizzles: */ \
/*	0K	+512	1K	+512	2K	+512	3K	+512 */ \
/*	4K */ \
/*	... */ \
/* On input the 64-byte cache lines hold these data values if pass 1 swizzles: */ \
/*	0K	+2K	512	+2K	1K	+2K	1536	+2K */ \
/*	4K */ \
/*	... */ \
/* On output the 64-byte cache lines hold these data values (+512 in high word): */ \
/*	0K	4K	8K	12K */ \
/*	1K	... */ \
/*	... */ \
/*	3K */ \
/*	16K	... */ \
/*	... */ \
\
/* Do 4 eight_complex_unfft8 operations */ \
/*	distance between fft data elements is 1K */ \
/*	do 4 macros each processing 32 data values */ \
\
c0b:	start_timer(16); \
rsi = (uintptr_t)g->data_addr;/* Get source address */ \
rcx = (uintptr_t)g->scratch_area;	/* Get address of scratch area */ \
rdi = (uintptr_t)g->premult_addr;/* Load premult/sin/cos pointer */ \
pfing(rbp = (rdi+(20-8)*2*clm*64)); /* Finish prefetching the sin/cos/premult data */ \
loops_init_prefetch(4*2*clm, 128, 1, rbp); /* 4*2*clm iterations */ \
/*bug	xtouchpages rbp, 4*2*clm * 128 / 1 ;; Prefetch TLBs */ \
rdi += 128;		/* Generates 1 byte offsets in radix-8 macro */ \
for(unsigned int loopA = 2; loopA; loopA--) { \
	for(unsigned int loopB = 2; loopB; loopB--) { \
		for(unsigned int loopC = 2*clm; loopC; loopC--) { \
			r8_sg4cl_eight_complex_unfft8(rsi, 64, blkdst, 2*blkdst, rcx, 64, clmblkdst, 2*clmblkdst, rdi-128); \
			rdi += XMM_SCD8;		/* Next premult/sin/cos data */ \
		} \
		rsi += -2*clm*64+4*blkdst;	/* Next source pointer */ \
		rcx += -2*clm*64+4*clmblkdst; /* Next dest pointer */ \
		rdi += -2*clm*XMM_SCD8;	/* Restore premult/sin/cos data pointer */ \
	} \
	rcx += -2*4*clmblkdst+clmblkdst8; \
} \
rsi = (uintptr_t)g->scratch_area;	/* Restore source pointer */ \
rdi += -128+2*clm*XMM_SCD8; /* Next sin/cos data pointer */ \
end_timer(16); \
\
/* Do inverse FFT levels 1,2,3 */ \
/* On input the 64-byte cache lines hold these data values: */ \
/*	0K	4K	8K	12K */ \
/*	1K	... */ \
/*	... */ \
/*	3K */ \
/*	16K	... */ \
/*	... */ \
/* On output the 64-byte cache lines hold these data values: */ \
/*	0K	16K	32K	48K */ \
/*	1K	... */ \
/*	... */ \
/*	15K */ \
\
/* Do 8 four_complex_unfft operations */ \
/*	distance between fft data elements is 8K */ \
/*	do 4 macros each processing 32 data values */ \
\
start_timer(20); \
set_data_prefetch_ptrs;		/* Init prefetching pointers */ \
loops_init(4*2*clm);		/* 4 iterations of 2*clm */ \
for(unsigned int loopA = 4*2*clm; loopA; ) { \
	xtouch(xptr(rcx+rbp+(clm-1)*128));	/* Preload TLB */ \
	do{ \
		r4_x4cl_four_complex_last_unfft4(rsi, 64, 4*clmblkdst, clmblkdst8, rdi, 4*2*clm*XMM_SCD4); \
		rdi += XMM_SCD4;		/* Next sin/cos/premult ptr */ \
	}while(--loopA & (2*clm-1)); \
	newxloop1doublepre(clm, c1b, c1a, clmblkdst); /* 2*clm iterations */ \
} \
end_timer(20); \
\
/* Normalize these values */ \
\
pass1_normalize(1, yes_fft); \
\
/* Copy data from the scratch area back to the FFT data area */ \
\
copy_scratch_data_to_fft(16, clm); \
goto	endlp;			/* Skip FFT code if flag not set */ \
\
/* Fire up auxiliary threads to do pass 1 forward FFT */ \
\
/* pass1: */	pass1_forward_fft_setup(b0b); \
\
/* Copy data from the FFT data area back to the scratch area */ \
\
b0b:	copy_fft_data_to_scratch(16, clm); \
\
/* Do FFT levels 1,2,3 */ \
/* */ \
/* On input the 64-byte cache lines hold these data values: */ \
/*	0K	16K	32K	48K */ \
/*	1K	... */ \
/*	... */ \
/*	15K */ \
/* On output the 64-byte cache lines hold these data values: */ \
/*	0K	8K	16K	24K */ \
/*	1K	... */ \
/*	... */ \
/*	7K	... */ \
/*	32K	... */ \
/*	... */ \
\
/* Do 8 four_complex_first_fft operations */ \
/*	distance between fft data elements is 8K */ \
/*	do 8 iterations each processing 16 data values */ \
\
yes_fft:start_timer(22); \
rsi = (uintptr_t)g->scratch_area;	/* Get address of scratch area */ \
rdi = (uintptr_t)g->premult_addr;/* Address of sin/cos table */ \
rdi += 2*clm*XMM_SCD8; \
set_data_prefetch_ptrs;		/* Init prefetching pointers */ \
pfing(rcx += 8*blkdst);		/* We've already prefetched 8 blocks */ \
loops_init(8*2*clm);		/* 8 iters of 2*clm */ \
for(unsigned int loopA = 8*2*clm; loopA; ) { \
	xtouch(xptr(rcx+rbp+(clm-1)*128));	/* Preload TLB */ \
	do{ \
		r4_x2cl_four_complex_first_fft4_scratch(rsi, 64, clmblkdst8, rdi); \
		rdi += XMM_SCD4;		/* Next sin/cos/premult ptr */ \
	}while(--loopA & (2*clm-1)); \
	newxloop1(clm, b1b, clmblkdst);	/* 2*clm iters */ \
} \
rsi += -8*clmblkdst;	/* Restore source pointer */ \
end_timer(22); \
\
/* Do FFT levels 4,5,6 */ \
/* */ \
/* On input the 64-byte cache lines hold these data values (+512 in high word): */ \
/*	0K	8K	16K	24K */ \
/*	1K	... */ \
/*	... */ \
/*	7K	... */ \
/*	32K	... */ \
/*	... */ \
/* On output the 64-byte cache lines hold these data values if pass 2 swizzles: */ \
/*	0K	+512	1K	+512	2K	+512	3K	+512 */ \
/*	4K	... */ \
/*	... */ \
/* On output the 64-byte cache lines hold these data values if pass 1 swizzles: */ \
/*	0K	+2K	512	+2K	1K	+2K	1536	+2K */ \
/*	4K	... */ \
/*	... */ \
\
/* Do 16 eight_complex_fft operations */ \
/*	distance between fft data elements is 1K */ \
/*	do 2 macros each processing 64 data values */ \
\
start_timer(26); \
rdi += -8*2*clm*XMM_SCD4-2*clm*XMM_SCD8+128; /* Next premult/sin/cos ptr (gen 1 byte offsets in radix-8 macro) */ \
rcx = (uintptr_t)g->data_addr;/* Destination pointer */ \
loops_init_prefetch(2*2*clm, 128, 1, rbp); /* 2 iters of 2*clm */ \
pfing(rbp = (uintptr_t)g->premult_prefetch);/* Load premult/sin/cos prefetch pointer */ \
for(unsigned int loopA = 2; loopA; loopA--) { \
	for(unsigned int loopB = 2*clm; loopB; loopB--) { \
		r8_sg8cl_eight_complex_fft8(rsi, 64, clmblkdst, 2*clmblkdst, 4*clmblkdst, rcx, 64, blkdst, 2*blkdst, 4*blkdst, rdi-128); \
		rdi += XMM_SCD8;		/* Next premult/sin/cos ptr */ \
	} \
	rsi += -2*clm*64+clmblkdst8; /* Next source pointer */ \
	rcx += -2*clm*64+8*blkdst;	/* Next dest pointer */ \
	rdi += -2*clm*XMM_SCD8;	/* Restore premult/sine/cosine pointer */ \
} \
end_timer(26); \
\
/* Work on next 4*clm sets of 64 values */ \
\
endlp:	xpass1_get_next_block(pass2, c0b, b0b); /* Get next block and jump */ \
/* done: */ xfft_footer_r; \
\
\
\
/* Perform a pass 1 negacyclic FFT of length 96 doubles. */ \
/* Memory examples are for a 96K FFT with 10 levels done in pass 2. */ \
/* Sin/cos data needed: 12*2*clm*XMM_SCD4 + 4*2*clm*XMM_SCD1 + 2*clm*XMM_SCD2 = 27*2*clm 64-byte cache lines */ \
/* BUG:	we under-prefetch the sin/cos data by 3*2*clm 64-byte cache lines */ \
\

#define r4_pass1sc96ac(procname, xpass2name, pass2_levels, clm) /* (fftname, pass2_levels, clmarg) */ \
\
/* set_FFT_constants(fftname, 0, 1, pass2_levels, clmarg); */ \
xfft_header_r(procname, xpass2name, pass2_levels, clm); \
\
/* Do pass 2 (the last 10 FFT levels if this is a 96K FFT) */ \
\
/* pass2: */	start_timer(0); \
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
/* On input the 64-byte cache lines hold these data values if pass 2 swizzles: */ \
/*	0K	+512	1K	+512	2K	+512	3K	+512 */ \
/*	4K */ \
/*	... */ \
/* On input the 64-byte cache lines hold these data values if pass 1 swizzles: */ \
/*	0K	+2K	512	+2K	1K	+2K	1536	+2K */ \
/*	4K */ \
/*	... */ \
/* On output the 64-byte cache lines hold these data values (+512 in high word): */ \
/*	0K	2K	4K	6K */ \
/*	1K	... */ \
/*	8K	... */ \
/*	... */ \
\
/* Do 12 four_complex_djbunfft operations */ \
/*	distance between fft data elements is 1K */ \
/*	do 12 macros each processing 16 data values */ \
\
c0b:	start_timer(16); \
rsi = (uintptr_t)g->data_addr;/* Get source address */ \
rcx = (uintptr_t)g->scratch_area;	/* Get address of scratch area */ \
rdi = (uintptr_t)g->premult_addr;/* Load premult/sin/cos pointer */ \
pfing(rbp = (rdi+(27-12)*2*clm*64)); /* Finish prefetching the sin/cos/premult data */ \
loops_init_prefetch(12*2*clm, 64, 1, rbp); /* 12*2*clm iterations */ \
/*bug	xtouchpages rbp, 12*2*clm * 64 / 1 ;; Prefetch TLBs */ \
for(unsigned int loopA = 12; loopA; loopA--) { \
	for(unsigned int loopB = 2*clm; loopB; loopB--) { \
		r4_sg2cl_four_complex_djbunfft(rsi, 64, blkdst, rcx, 64, clmblkdst, rdi); \
		rdi += XMM_SCD2;		/* Next premult/sin/cos data */ \
	} \
	rsi += -2*clm*64+2*blkdst;	/* Next source pointer */ \
	rcx += -2*clm*64+2*clmblkdst; /* Next dest pointer */ \
	rdi += -2*clm*XMM_SCD2;	/* Restore premult/sin/cos data pointer */ \
} \
rsi = (uintptr_t)g->scratch_area;	/* Restore source pointer */ \
rdi += 2*clm*XMM_SCD2;	/* Next premult/sin/cos data pointer */ \
end_timer(16); \
\
/* Do inverse FFT levels 4,5 (actually 1.585 levels) */ \
/* On input the 64-byte cache lines hold these data values: */ \
/*	0K	2K	4K	6K */ \
/*	1K	... */ \
/*	8K	... */ \
/*	... */ \
/* On output the 64-byte cache lines hold these data values: */ \
/*	0K	6K	12K	18K */ \
/*	1K	... */ \
/*	... */ \
/*	5K */ \
/*	24K	... */ \
/*	... */ \
\
/* Do 16 three_complex_unfft operations */ \
/*	distance between fft data elements is 4K */ \
/*	do 8 macros each processing 24 data values */ \
\
start_timer(17); \
set_data_prefetch_ptrs;		/* Init prefetching pointers */ \
loops_init(8*2*clm);		/* 8*2*clm iterations */ \
r3_x3cl_three_complex_djbunfft_preload; \
for(unsigned int loopA = 4*2*2*clm; loopA; ) { \
	xtouch(xptr(rcx+rbp+(clm-1)*128));	/* Preload TLB */ \
	do{ \
		do{ \
			r3_x3cl_three_complex_djbunfft(rsi, 64, 2*clmblkdst, rdi, 2*2*clm*XMM_SCD1); \
			rdi += XMM_SCD1;		/* Next clm's sin/cos data */ \
		}while(--loopA & (2*clm-1)); \
		newxloop2(clm, c7b, clmblkdst);	/* 2*2*clm iterations */ \
	}while(loopA & (2*2*clm-1)); \
	rsi += -2*clmblkdst+6*clmblkdst; /* Next source pointer */ \
	rdi += -2*2*clm*XMM_SCD1;	/* Restore sine/cosine pointer */ \
} \
rsi += -4*6*clmblkdst;	/* Restore source pointer */ \
rdi += 4*2*clm*XMM_SCD1;	/* Next sine/cosine pointer */ \
end_timer(17); \
\
/* Do inverse FFT levels 1,2,3 */ \
/* On input the 64-byte cache lines hold these data values: */ \
/*	0K	6K	12K	18K */ \
/*	1K	... */ \
/*	... */ \
/*	5K */ \
/*	24K	... */ \
/*	... */ \
/* On output the 64-byte cache lines hold these data values: */ \
/*	0K	24K	48K	72K */ \
/*	1K	... */ \
/*	... */ \
/*	23K */ \
\
/* Do 12 four_complex_unfft operations */ \
/*	distance between fft data elements is 12K */ \
/*	do 6 macros each processing 32 data values */ \
\
start_timer(20); \
loops_init(6*2*clm);		/* 6 iterations of 2*clm */ \
for(unsigned int loopA = 6*2*clm; loopA; ) { \
	xtouch(xptr(rcx+rbp+(clm-1)*128));	/* Preload TLB */ \
	do{ \
		r4_x4cl_four_complex_last_unfft4(rsi, 64, 6*clmblkdst, 12*clmblkdst, rdi, 6*2*clm*XMM_SCD4); \
		rdi += XMM_SCD4;		/* Next sin/cos/premult ptr */ \
	}while(--loopA & (2*clm-1)); \
	newxloop1(clm, c1b, clmblkdst);	/* 2*clm iterations */ \
} \
end_timer(20); \
\
/* Normalize these values */ \
\
pass1_normalize(1, yes_fft); \
\
/* Copy data from the scratch area back to the FFT data area */ \
\
copy_scratch_no_pad_data_to_fft(24, clm); \
goto	endlp;			/* Skip FFT code if flag not set */ \
\
/* Fire up auxiliary threads to do pass 1 forward FFT */ \
\
/* pass1: */	pass1_forward_fft_setup(b0b); \
\
/* Copy data from the FFT data area back to the scratch area */ \
\
b0b:	copy_fft_data_to_scratch_no_pad(24, clm); \
\
/* Do FFT levels 1,2,3 */ \
/* */ \
/* On input the 64-byte cache lines hold these data values: */ \
/*	0K	24K	48K	72K */ \
/*	1K	... */ \
/*	... */ \
/*	23K */ \
/* On output the 64-byte cache lines hold these data values: */ \
/*	0K	12K	24K	36K */ \
/*	1K	... */ \
/*	... */ \
/*	11K	... */ \
/*	48K	... */ \
/*	... */ \
\
/* Do 12 four_complex_first_fft operations */ \
/*	distance between fft data elements is 12K */ \
/*	do 12 iterations each processing 16 data values */ \
\
yes_fft:start_timer(22); \
rsi = (uintptr_t)g->scratch_area;	/* Get address of scratch area */ \
rdi = (uintptr_t)g->premult_addr;/* Address of premultiplier/sin/cos table */ \
rdi += 4*2*clm*XMM_SCD1 + 2*clm*XMM_SCD2; \
set_data_prefetch_ptrs;		/* Init prefetching pointers */ \
pfing(rcx += 10*blkdst);		/* We've already prefetched 10 blocks */ \
loops_init(12*2*clm);		/* 6 iters of 2*2*clm */ \
for(unsigned int loopA = 6*2*2*clm; loopA; ) { \
	xtouch(xptr(rcx+rbp+(clm-1)*128));	/* Preload TLB */ \
	do{ \
		do{ \
			r4_x2cl_four_complex_first_fft4_scratch(rsi, 64, 12*clmblkdst, rdi); \
			rdi += XMM_SCD4;		/* Next sin/cos/premult ptr */ \
		}while(--loopA & (2*clm-1)); \
		newxloop2(clm, b1b, clmblkdst);	/* 2*2*clm iters prefetching 2*clm cache lines */ \
	}while(loopA & (2*2*clm-1)); \
} \
rsi += -6*2*clmblkdst;	/* Restore source pointer */ \
end_timer(22); \
\
/* Do FFT levels 4,5 (actually 1.585 levels) */ \
/* */ \
/* On input the 64-byte cache lines hold these data values: */ \
/*	0K	12K	24K	36K */ \
/*	1K	... */ \
/*	... */ \
/*	11K	... */ \
/*	48K	... */ \
/*	... */ \
/* On output the 64-byte cache lines hold these data values: */ \
/*	0K	4K	8K	12K */ \
/*	1K	... */ \
/*	... */ \
/*	3K	... */ \
/*	16K	... */ \
/*	... */ \
\
/* Do 16 three_complex_fft operations */ \
/*	distance between fft data elements is 4K */ \
/*	do 8 macros each processing 24 data values */ \
\
start_timer(25); \
rdi += -12*2*clm*XMM_SCD4 - 4*2*clm*XMM_SCD1; \
loops_init(8*2*clm);		/* 2 iterations of 4 iters of 2*clm */ \
r3_x3cl_three_complex_djbfft_preload; \
for(unsigned int loopA = 2; loopA; loopA--) { \
	for(unsigned int loopB = 4*2*clm; loopB; ) { \
		xtouch(xptr(rcx+rbp+(clm-1)*128));	/* Preload TLB */ \
		do{ \
			r3_x3cl_three_complex_djbfft(rsi, 64, 4*clmblkdst, rdi); \
			rdi += XMM_SCD1;		/* Next sine/cosine pointer */ \
		}while(--loopB & (2*clm-1)); \
		newxloop1(clm, b7b, clmblkdst);	/* 2*clm iterations */ \
	} \
	rsi += -4*clmblkdst+12*clmblkdst; /* Next source pointer */ \
	rdi += -4*2*clm*XMM_SCD1;	/* Restore sin/cos ptr */ \
} \
rsi += -2*12*clmblkdst;	/* Restore source pointer */ \
end_timer(25); \
\
/* Do FFT levels 6,7 */ \
/* */ \
/* On input the 64-byte cache lines hold these data values (+512 in high word): */ \
/*	0K	4K	8K	12K */ \
/*	1K	... */ \
/*	... */ \
/*	3K	... */ \
/*	16K	... */ \
/*	... */ \
/* On output the 64-byte cache lines hold these data values if pass 2 swizzles: */ \
/*	0K	+512	1K	+512	2K	+512	3K	+512 */ \
/*	4K	... */ \
/*	... */ \
/* On output the 64-byte cache lines hold these data values if pass 1 swizzles: */ \
/*	0K	+2K	512	+2K	1K	+2K	1536	+2K */ \
/*	4K	... */ \
/*	... */ \
\
/* Do 12 four_complex_fft operations */ \
/*	distance between fft data elements is 1K */ \
/*	do 6 macros each processing 32 data values */ \
\
start_timer(26); \
rdi += -2*clm*XMM_SCD2;	/* Next premult/sin/cos ptr */ \
rcx = (uintptr_t)g->data_addr;/* Destination pointer */ \
loops_init_prefetch(6*2*clm, 128, 1, rbp); /* 6 iters of 2*clm */ \
pfing(rbp = (uintptr_t)g->premult_prefetch);/* Load premult/sin/cos prefetch pointer */ \
r4_sg4cl_four_complex_djbfft_preload; \
for(unsigned int loopA = 6; loopA; loopA--) { \
	for(unsigned int loopB = 2*clm; loopB; loopB--) { \
		r4_sg4cl_four_complex_djbfft(rsi, 64, clmblkdst, 2*clmblkdst, rcx, 64, blkdst, 2*blkdst, rdi); \
		rdi += XMM_SCD2;		/* Next premult/sin/cos ptr */ \
	} \
	rsi += -2*clm*64+4*clmblkdst; /* Next source pointer */ \
	rcx += -2*clm*64+4*blkdst;	/* Next dest pointer */ \
	rdi += -2*clm*XMM_SCD2;	/* Restore premult/sine/cosine pointer */ \
} \
end_timer(26); \
\
/* Work on next 4*clm sets of 96 values */ \
\
endlp:	xpass1_get_next_block(pass2, c0b, b0b); /* Get next block and jump */ \
/* done: */ xfft_footer_r; \
\
\
\
/* Perform a pass 1 negacyclic FFT of length 128 doubles. */ \
/* Memory examples are for a 128K FFT with 10 levels done in pass 2. */ \
/* Sin/cos data needed: 16*2*clm*XMM_SCD4 + 4*2*clm*XMM_SCD2 + 2*clm*XMM_SCD2 = 37*2*clm 64-byte cache lines */ \
/* BUG:	we under-prefetch the sin/cos data by 5*2*clm 64-byte cache lines */ \
\

#define r4_pass1sc128ac(procname, xpass2name, pass2_levels, clm) /* (fftname, pass2_levels, clmarg) */ \
\
/* set_FFT_constants(fftname, 0, 1, pass2_levels, clmarg); */ \
xfft_header_r(procname, xpass2name, pass2_levels, clm); \
\
/* Do pass 2 (the last 10 FFT levels if this is a 128K FFT) */ \
\
/* pass2: */	start_timer(0); \
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
/* On input the 64-byte cache lines hold these data values if pass 2 swizzles: */ \
/*	0K	+512	1K	+512	2K	+512	3K	+512 */ \
/*	4K */ \
/*	... */ \
/* On input the 64-byte cache lines hold these data values if pass 1 swizzles: */ \
/*	0K	+2K	512	+2K	1K	+2K	1536	+2K */ \
/*	4K */ \
/*	... */ \
/* On output the 64-byte cache lines hold these data values (+512 in high word): */ \
/*	0K	2K	4K	6K */ \
/*	1K	... */ \
/*	8K	... */ \
/*	... */ \
\
/* Do 16 four_complex_djbunfft operations */ \
/*	distance between fft data elements is 1K */ \
/*	do 16 macros each processing 16 data values */ \
\
c0b:	start_timer(16); \
rsi = (uintptr_t)g->data_addr;/* Get source address */ \
rcx = (uintptr_t)g->scratch_area;	/* Get address of scratch area */ \
rdi = (uintptr_t)g->premult_addr;/* Load premult/sin/cos pointer */ \
pfing(rbp = (rdi+(37-16)*2*clm*64)); /* Finish prefetching the sin/cos/premult data */ \
loops_init_prefetch(16*2*clm, 64, 1, rbp); /* 16*2*clm iterations */ \
/*bug	xtouchpages rbp, 16*2*clm * 64 / 1 ;; Prefetch TLBs */ \
for(unsigned int loopA = 4; loopA; loopA--) { \
	for(unsigned int loopB = 4; loopB; loopB--) { \
		for(unsigned int loopC = 2*clm; loopC; loopC--) { \
			r4_sg2cl_four_complex_djbunfft(rsi, 64, blkdst, rcx, 64, clmblkdst, rdi); \
			rdi += XMM_SCD2;		/* Next premult/sin/cos data */ \
		} \
		rsi += -2*clm*64+2*blkdst;	/* Next source pointer */ \
		rcx += -2*clm*64+2*clmblkdst; /* Next dest pointer */ \
		rdi += -2*clm*XMM_SCD2;	/* Restore premult/sin/cos data pointer */ \
	} \
	rcx += -4*2*clmblkdst+clmblkdst8; \
} \
rsi = (uintptr_t)g->scratch_area;	/* Restore source pointer */ \
rdi += 2*clm*XMM_SCD2;	/* Next premult/sin/cos data pointer */ \
end_timer(16); \
\
/* Do inverse FFT levels 4,5 */ \
/* On input the 64-byte cache lines hold these data values: */ \
/*	0K	2K	4K	6K */ \
/*	1K	... */ \
/*	8K	... */ \
/*	... */ \
/* On output the 64-byte cache lines hold these data values: */ \
/*	0K	8K	16K	24K */ \
/*	1K	... */ \
/*	... */ \
/*	7K */ \
/*	32K	... */ \
/*	... */ \
\
/* Do 16 four_complex_unfft operations */ \
/*	distance between fft data elements is 4K */ \
/*	do 8 macros each processing 32 data values */ \
\
start_timer(17); \
set_data_prefetch_ptrs;		/* Init prefetching pointers */ \
loops_init(8*2*clm);		/* 8*2*clm iterations */ \
r4_x4cl_four_complex_djbunfft_preload; \
for(unsigned int loopA = 4; loopA; loopA--) { \
	for(unsigned int loopB = 2*2*clm; loopB; ) { \
		xtouch(xptr(rcx+rbp+(clm-1)*128));	/* Preload TLB */ \
		do{ \
			r4_x4cl_four_complex_djbunfft(rsi, 64, 2*clmblkdst, 4*clmblkdst, rdi, 2*2*clm*XMM_SCD2); \
			rdi += XMM_SCD2;		/* Next clm's sin/cos data */ \
		}while(--loopB & (2*clm-1)); \
		newxloop1(clm, c7b, clmblkdst);	/* 2*clm iterations */ \
	} \
	rsi += -2*clmblkdst+clmblkdst8; /* Next source pointer */ \
	rdi += -2*2*clm*XMM_SCD2;	/* Restore sine/cosine pointer */ \
} \
rsi += -4*clmblkdst8;	/* Restore source pointer */ \
rdi += 4*2*clm*XMM_SCD2;	/* Next sine/cosine pointer */ \
end_timer(17); \
\
/* Do inverse FFT levels 1,2,3 */ \
/* On input the 64-byte cache lines hold these data values: */ \
/*	0K	8K	16K	24K */ \
/*	1K	... */ \
/*	... */ \
/*	7K */ \
/*	32K	... */ \
/*	... */ \
/* On output the 64-byte cache lines hold these data values: */ \
/*	0K	32K	64K	96K */ \
/*	1K	... */ \
/*	... */ \
/*	31K */ \
\
/* Do 16 four_complex_unfft operations */ \
/*	distance between fft data elements is 16K */ \
/*	do 8 macros each processing 32 data values */ \
\
start_timer(20); \
loops_init(8*2*clm);		/* 8 iterations of 2*clm */ \
for(unsigned int loopA = 8*2*clm; loopA; ) { \
	xtouch(xptr(rcx+rbp+(clm-1)*128));	/* Preload TLB */ \
	do{ \
		r4_x4cl_four_complex_last_unfft4(rsi, 64, clmblkdst8, 2*clmblkdst8, rdi, 8*2*clm*XMM_SCD4); \
		rdi += XMM_SCD4;		/* Next sin/cos/premult ptr */ \
	}while(--loopA & (2*clm-1)); \
	newxloop1(clm, c1b, clmblkdst);	/* 2*clm iterations */ \
} \
end_timer(20); \
\
/* Normalize these values */ \
\
pass1_normalize(1, yes_fft); \
\
/* Copy data from the scratch area back to the FFT data area */ \
\
copy_scratch_data_to_fft(32, clm); \
goto	endlp;			/* Skip FFT code if flag not set */ \
\
/* Fire up auxiliary threads to do pass 1 forward FFT */ \
\
/* pass1: */	pass1_forward_fft_setup(b0b); \
\
/* Copy data from the FFT data area back to the scratch area */ \
\
b0b:	copy_fft_data_to_scratch(32, clm); \
\
/* Do FFT levels 1,2,3 */ \
/* */ \
/* On input the 64-byte cache lines hold these data values: */ \
/*	0K	32K	64K	96K */ \
/*	1K	... */ \
/*	... */ \
/*	31K */ \
/* On output the 64-byte cache lines hold these data values: */ \
/*	0K	16K	32K	48K */ \
/*	1K	... */ \
/*	... */ \
/*	15K	... */ \
/*	64K	... */ \
/*	... */ \
\
/* Do 16 four_complex_first_fft operations */ \
/*	distance between fft data elements is 16K */ \
/*	do 16 iterations each processing 16 data values */ \
\
yes_fft:start_timer(22); \
rsi = (uintptr_t)g->scratch_area;	/* Get address of scratch area */ \
rdi = (uintptr_t)g->premult_addr;/* Address of premultiplier/sin/cos table */ \
rdi += 4*2*clm*XMM_SCD2 + 2*clm*XMM_SCD2; \
set_data_prefetch_ptrs;		/* Init prefetching pointers */ \
pfing(rcx += 16*blkdst);		/* We've already prefetched 16 blocks */ \
loops_init(16*2*clm);		/* 2 iters of 4 iters of 2*2*clm */ \
for(unsigned int loopA = 2; loopA; loopA--) { \
	for(unsigned int loopB = 4*2*2*clm; loopB; ) { \
		xtouch(xptr(rcx+rbp+(clm-1)*128));	/* Preload TLB */ \
		do{ \
			do{ \
				r4_x2cl_four_complex_first_fft4_scratch(rsi, 64, 2*clmblkdst8, rdi); \
				rdi += XMM_SCD4;		/* Next sin/cos/premult ptr */ \
			}while(--loopB & (2*clm-1)); \
			newxloop2(clm, b1b, clmblkdst);	/* 2*2*clm iters prefetching 2*clm cache lines */ \
		}while(loopB & (2*2*clm-1)); \
	} \
	rsi += -4*2*clmblkdst+clmblkdst8; \
} \
rsi += -2*clmblkdst8;	/* Restore source pointer */ \
end_timer(22); \
\
/* Do FFT levels 4,5 */ \
/* */ \
/* On input the 64-byte cache lines hold these data values: */ \
/*	0K	16K	32K	48K */ \
/*	1K	... */ \
/*	... */ \
/*	15K	... */ \
/*	64K	... */ \
/*	... */ \
/* On output the 64-byte cache lines hold these data values: */ \
/*	0K	4K	8K	12K */ \
/*	1K	... */ \
/*	... */ \
/*	3K	... */ \
/*	16K	... */ \
/*	... */ \
\
/* Do 16 four_complex_fft operations */ \
/*	distance between fft data elements is 4K */ \
/*	do 8 macros each processing 32 data values */ \
\
start_timer(25); \
rdi += -16*2*clm*XMM_SCD4 - 4*2*clm*XMM_SCD2; \
loops_init(8*2*clm);		/* 2 iterations of 4 iters of 2*clm */ \
r4_x4cl_four_complex_djbfft_preload; \
for(unsigned int loopA = 2; loopA; loopA--) { \
	for(unsigned int loopB = 4*2*clm; loopB; ) { \
		xtouch(xptr(rcx+rbp+(clm-1)*128));	/* Preload TLB */ \
		do{ \
			r4_x4cl_four_complex_djbfft(rsi, 64, 4*clmblkdst, clmblkdst8, rdi); \
			rdi += XMM_SCD2;		/* Next sine/cosine pointer */ \
		}while(--loopB & (2*clm-1)); \
		newxloop1(clm, b7b, clmblkdst);	/* 2*clm iterations */ \
	} \
	rsi += -4*clmblkdst+2*clmblkdst8; /* Next source pointer */ \
	rdi += -4*2*clm*XMM_SCD2;	/* Restore sin/cos ptr */ \
} \
rsi += -2*2*clmblkdst8;	/* Restore source pointer */ \
end_timer(25); \
\
/* Do FFT levels 6,7 */ \
/* */ \
/* On input the 64-byte cache lines hold these data values (+512 in high word): */ \
/*	0K	4K	8K	12K */ \
/*	1K	... */ \
/*	... */ \
/*	3K	... */ \
/*	16K	... */ \
/*	... */ \
/* On output the 64-byte cache lines hold these data values if pass 2 swizzles: */ \
/*	0K	+512	1K	+512	2K	+512	3K	+512 */ \
/*	4K	... */ \
/*	... */ \
/* On output the 64-byte cache lines hold these data values if pass 1 swizzles: */ \
/*	0K	+2K	512	+2K	1K	+2K	1536	+2K */ \
/*	4K	... */ \
/*	... */ \
\
/* Do 16 four_complex_fft operations */ \
/*	distance between fft data elements is 1K */ \
/*	do 8 macros each processing 32 data values */ \
\
start_timer(26); \
rdi += -2*clm*XMM_SCD2;	/* Next premult/sin/cos ptr */ \
rcx = (uintptr_t)g->data_addr;/* Destination pointer */ \
loops_init_prefetch(8*2*clm, 128, 1, rbp); /* 8 iters of 2*clm */ \
pfing(rbp = (uintptr_t)g->premult_prefetch);/* Load premult/sin/cos prefetch pointer */ \
r4_sg4cl_four_complex_djbfft_preload; \
for(unsigned int loopA = 4; loopA; loopA--) { \
	for(unsigned int loopB = 2; loopB; loopB--) { \
		for(unsigned int loopC = 2*clm; loopC; loopC--) { \
			r4_sg4cl_four_complex_djbfft(rsi, 64, clmblkdst, 2*clmblkdst, rcx, 64, blkdst, 2*blkdst, rdi); \
			rdi += XMM_SCD2;		/* Next premult/sin/cos ptr */ \
		} \
		rsi += -2*clm*64+4*clmblkdst; /* Next source pointer */ \
		rcx += -2*clm*64+4*blkdst;	/* Next dest pointer */ \
		rdi += -2*clm*XMM_SCD2;	/* Restore premult/sine/cosine pointer */ \
	} \
	rsi += -2*4*clmblkdst+clmblkdst8; /* Next source pointer */ \
} \
end_timer(26); \
\
/* Work on next 4*clm sets of 128 values */ \
\
endlp:	xpass1_get_next_block(pass2, c0b, b0b); /* Get next block and jump */ \
/* done: */ xfft_footer_r; \
\
\
\
/* Perform a pass 1 negacyclic FFT of length 256 doubles. */ \
/* Memory examples are for a 256K FFT with 10 levels done in pass 2. */ \
/* Sin/cos data needed: 32*2*clm*XMM_SCD4 + 8*2*clm*XMM_SCD2 + 2*clm*XMM_SCD8 = 76*2*clm 64-byte cache lines */ \
/* BUG - we substantially underprefetch the sin/cos data */ \
\

#define r4_pass1sc256ac(procname, xpass2name, pass2_levels, clm) /* (fftname, pass2_levels, clmarg) */ \
\
/* set_FFT_constants(fftname, 0, 1, pass2_levels, clmarg); */ \
xfft_header_r(procname, xpass2name, pass2_levels, clm); \
\
/* Do pass 2 (the last 10 FFT levels if this is a 256K FFT) */ \
\
/* pass2: */	start_timer(0); \
xpass2name(g); \
end_timer(0); \
if(g->ffttype == 1)		/* We're done if FFTing only */ \
	return; /* goto done; */ \
\
/* Do the inverse FFT */ \
\
pass1_inverse_fft_setup(c0b, procname); \
\
/* Do inverse FFT levels 6,7,8 */ \
/* On input the 64-byte cache lines hold these data values if pass 2 swizzles: */ \
/*	0K	+512	1K	+512	2K	+512	3K	+512 */ \
/*	4K */ \
/*	... */ \
/* On input the 64-byte cache lines hold these data values if pass 1 swizzles: */ \
/*	0K	+2K	512	+2K	1K	+2K	1536	+2K */ \
/*	4K */ \
/*	... */ \
/* On output the 64-byte cache lines hold these data values (+512 in high word): */ \
/*	0K	4K	8K	12K */ \
/*	1K	... */ \
/*	... */ \
/*	3K */ \
/*	16K	... */ \
/*	... */ \
\
/* Do 16 eight_complex_unfft8 operations */ \
/*	distance between fft data elements is 1K */ \
/*	do 16 macros each processing 32 data values */ \
\
c0b:	start_timer(16); \
rsi = (uintptr_t)g->data_addr;/* Get source address */ \
rcx = (uintptr_t)g->scratch_area;	/* Get address of scratch area */ \
rdi = (uintptr_t)g->premult_addr;/* Load premult/sin/cos pointer */ \
pfing(rbp = (rdi+(76-32)*2*clm*64)); /* Finish prefetching the sin/cos/premult data */ \
loops_init_prefetch(16*2*clm, 128, 1, rbp); /* 16*2*clm iterations */ \
/*bug	xtouchpages rbp, 16*2*clm * 128 / 1 ;; Prefetch TLBs */ \
rdi += 128;		/* Generates 1 byte offsets in radix-8 macro */ \
for(unsigned int loopA = 8; loopA; loopA--) { \
	for(unsigned int loopB = 2; loopB; loopB--) { \
		for(unsigned int loopC = 2*clm; loopC; loopC--) { \
			r8_sg4cl_eight_complex_unfft8(rsi, 64, blkdst, 2*blkdst, rcx, 64, clmblkdst, 2*clmblkdst, rdi-128); \
			rdi += XMM_SCD8;		/* Next premult/sin/cos data */ \
		} \
		rsi += -2*clm*64+4*blkdst;	/* Next source pointer */ \
		rcx += -2*clm*64+4*clmblkdst; /* Next dest pointer */ \
		rdi += -2*clm*XMM_SCD8;	/* Restore premult/sin/cos data pointer */ \
	} \
	rcx += -2*4*clmblkdst+clmblkdst8; \
} \
rsi = (uintptr_t)g->scratch_area;	/* Restore source pointer */ \
rdi += -128+2*clm*XMM_SCD8; /* Next sin/cos data pointer */ \
end_timer(16); \
\
/* Do inverse FFT levels 4,5 */ \
/* On input the 64-byte cache lines hold these data values: */ \
/*	0K	4K	8K	12K */ \
/*	1K	... */ \
/*	... */ \
/*	3K */ \
/*	16K	... */ \
/*	... */ \
/* On output the 64-byte cache lines hold these data values: */ \
/*	0K	16K	32K	48K */ \
/*	1K	... */ \
/*	... */ \
/*	15K */ \
/*	64K	... */ \
/*	... */ \
\
/* Do 32 four_complex_unfft operations */ \
/*	distance between fft data elements is 8K */ \
/*	do 16 macros each processing 32 data values */ \
\
start_timer(17); \
set_data_prefetch_ptrs;		/* Init prefetching pointers */ \
loops_init(16*2*clm);		/* 16 iterations of 2*clm */ \
r4_x4cl_four_complex_djbunfft_preload; \
for(unsigned int loopA = 4; loopA; loopA--) { \
	for(unsigned int loopB = 4*2*clm; loopB; ) { \
		xtouch(xptr(rcx+rbp+(clm-1)*128));	/* Preload TLB */ \
		do{ \
			r4_x4cl_four_complex_djbunfft(rsi, 64, 4*clmblkdst, clmblkdst8, rdi, 4*2*clm*XMM_SCD2); \
			rdi += XMM_SCD2;		/* Next clm's sin/cos data */ \
		}while(--loopB & (2*clm-1)); \
		newxloop1(clm, c7b, clmblkdst);	/* Loop 2*clm times */ \
	} \
	rsi += -4*clmblkdst+2*clmblkdst8; /* Next source pointer */ \
	rdi += -4*2*clm*XMM_SCD2;	/* Restore sine/cosine pointer */ \
} \
rsi += -4*2*clmblkdst8;	/* Restore source pointer */ \
rdi += 8*2*clm*XMM_SCD2;	/* Next sin/cos data pointer */ \
end_timer(17); \
\
/* Do inverse FFT levels 1,2,3 */ \
/* On input the 64-byte cache lines hold these data values: */ \
/*	0K	16K	32K	48K */ \
/*	1K	... */ \
/*	... */ \
/*	15K */ \
/*	64K	... */ \
/*	... */ \
/* On output the 64-byte cache lines hold these data values: */ \
/*	0K	64K	128K	192K */ \
/*	1K	... */ \
/*	... */ \
/*	63K */ \
\
/* Do 32 four_complex_unfft operations */ \
/*	distance between fft data elements is 32K */ \
/*	do 16 macros each processing 32 data values */ \
\
start_timer(20); \
loops_init(16*2*clm);		/* 16 iterations of 2*clm */ \
for(unsigned int loopA = 2; loopA; loopA--) { \
	for(unsigned int loopB = 8*2*clm; loopB; ) { \
		xtouch(xptr(rcx+rbp+(clm-1)*128));	/* Preload TLB */ \
		do{ \
			r4_x4cl_four_complex_last_unfft4(rsi, 64, 2*clmblkdst8, 4*clmblkdst8, rdi, 16*2*clm*XMM_SCD4); \
			rdi += XMM_SCD4;		/* Next sin/cos/premult ptr */ \
		}while(--loopB & (2*clm-1)); \
		newxloop1(clm, c1b, clmblkdst);	/* 2*clm iterations */ \
	} \
	rsi += -8*clmblkdst+clmblkdst8; /* Next source pointer */ \
} \
end_timer(20); \
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
/*	0K	64K	128K	192K */ \
/*	1K	... */ \
/*	... */ \
/*	63K */ \
/* On output the 64-byte cache lines hold these data values: */ \
/*	0K	32K	64K	96K */ \
/*	1K	... */ \
/*	... */ \
/*	31K	... */ \
/*	128K	... */ \
/*	... */ \
\
/* Do 32 four_complex_first_fft operations */ \
/*	distance between fft data elements is 32K */ \
/*	do 32 iterations each processing 16 data values */ \
\
yes_fft:start_timer(22); \
rsi = (uintptr_t)g->scratch_area;	/* Get address of scratch area */ \
rdi = (uintptr_t)g->premult_addr;/* Address of sin/cos table */ \
rdi += 8*2*clm*XMM_SCD2 + 2*clm*XMM_SCD8; \
set_data_prefetch_ptrs;		/* Init prefetching pointers */ \
pfing(rcx += 32*blkdst);		/* We've already prefetched 32 blocks */ \
loops_init(32*2*clm);		/* 4 iters of 4 iters of 2*2*clm */ \
for(unsigned int loopA = 4; loopA; loopA--) { \
	for(unsigned int loopB = 4*2*2*clm; loopB; ) { \
		xtouch(xptr(rcx+rbp+(clm-1)*128));	/* Preload TLB */ \
		do{ \
			do{ \
				r4_x2cl_four_complex_first_fft4_scratch(rsi, 64, 4*clmblkdst8, rdi); \
				rdi += XMM_SCD4;		/* Next sin/cos/premult ptr */ \
			}while(--loopB & (2*clm-1)); \
			newxloop2(clm, b1b, clmblkdst);	/* 2*2*clm iters prefetching 2*clm cache lines */ \
		}while(--loopB & (2*2*clm-1)); \
	} \
	rsi += -4*2*clmblkdst+clmblkdst8; \
} \
rsi += -4*clmblkdst8;	/* Restore source pointer */ \
end_timer(22); \
\
/* Do FFT levels 4,5 */ \
/* */ \
/* On input the 64-byte cache lines hold these data values: */ \
/*	0K	32K	64K	96K */ \
/*	1K	... */ \
/*	... */ \
/*	31K	... */ \
/*	128K	... */ \
/*	... */ \
/* On output the 64-byte cache lines hold these data values: */ \
/*	0K	8K	16K	24K */ \
/*	1K	... */ \
/*	... */ \
/*	7K	... */ \
/*	32K	... */ \
/*	... */ \
\
/* Do 32 four_complex_fft operations */ \
/*	distance between fft data elements is 8K */ \
/*	do 16 macros each processing 32 data values */ \
\
start_timer(25); \
rdi += -32*2*clm*XMM_SCD4 - 8*2*clm*XMM_SCD2; /* Next sin/cos ptr */ \
loops_init(16*2*clm);		/* 2 iterations of 8 iters of 2*clm */ \
r4_x4cl_four_complex_djbfft_preload; \
for(unsigned int loopA = 2; loopA; loopA--) { \
	for(unsigned int loopB = 8*2*clm; loopB; ) { \
		xtouch(xptr(rcx+rbp+(clm-1)*128));	/* Preload TLB */ \
		do{ \
			r4_x4cl_four_complex_djbfft(rsi, 64, clmblkdst8, 2*clmblkdst8, rdi); \
			rdi += XMM_SCD2;		/* Next sine/cosine pointer */ \
		}while(--loopB & (2*clm-1)); \
		newxloop1(clm, b7b, clmblkdst);	/* 2*clm iterations */ \
	} \
	rsi += -8*clmblkdst+4*clmblkdst8; /* Next source pointer */ \
	rdi += -8*2*clm*XMM_SCD2;	/* Restore sin/cos ptr */ \
} \
rsi += -2*4*clmblkdst8;	/* Restore source pointer */ \
end_timer(25); \
\
/* Do FFT levels 6,7,8 */ \
/* */ \
/* On input the 64-byte cache lines hold these data values (+512 in high word): */ \
/*	0K	8K	16K	24K */ \
/*	1K	... */ \
/*	... */ \
/*	7K	... */ \
/*	32K	... */ \
/*	... */ \
/* On output the 64-byte cache lines hold these data values if pass 2 swizzles: */ \
/*	0K	+512	1K	+512	2K	+512	3K	+512 */ \
/*	4K	... */ \
/*	... */ \
/* On output the 64-byte cache lines hold these data values if pass 1 swizzles: */ \
/*	0K	+2K	512	+2K	1K	+2K	1536	+2K */ \
/*	4K	... */ \
/*	... */ \
\
/* Do 64 eight_complex_fft operations */ \
/*	distance between fft data elements is 1K */ \
/*	do 32 macros each processing 64 data values */ \
\
start_timer(26); \
rdi += -2*clm*XMM_SCD8+128; /* Next premult/sin/cos ptr (gen 1 byte offsets in radix-8 macro) */ \
rcx = (uintptr_t)g->data_addr;/* Destination pointer */ \
loops_init_prefetch(8*2*clm, 128, 1, rbp); /* 8 iters of 2*clm */ \
pfing(rbp = (uintptr_t)g->premult_prefetch);/* Load premult/sin/cos prefetch pointer */ \
for(unsigned int loopA = 8; loopA; loopA--) { \
	for(unsigned int loopB = 2*clm; loopB; loopB--) { \
		r8_sg8cl_eight_complex_fft8(rsi, 64, clmblkdst, 2*clmblkdst, 4*clmblkdst, rcx, 64, blkdst, 2*blkdst, 4*blkdst, rdi-128); \
		rdi += XMM_SCD8;		/* Next premult/sin/cos ptr */ \
	} \
	rsi += -2*clm*64+clmblkdst8; /* Next source pointer */ \
	rcx += -2*clm*64+8*blkdst;	/* Next dest pointer */ \
	rdi += -2*clm*XMM_SCD8;	/* Restore premult/sine/cosine pointer */ \
} \
end_timer(26); \
\
/* Work on next 4*clm sets of 256 values */ \
\
endlp:	xpass1_get_next_block(pass2, c0b, b0b); /* Get next block and jump */ \
/* done: */ xfft_footer_r; \
\
\
