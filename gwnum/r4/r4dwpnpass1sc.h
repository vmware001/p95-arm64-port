#pragma once

#define xfft_header_rn(procname, xpass2name, pass2_levels, clm) \
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

#define xfft_footer_rn \
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
#define xpass1_get_next_block(pass2, c0b, b0b, wpn) \
do{ \
	branch1 = g->pass1_get_next_block(g); \
	if(branch1 & 2) /* 2 or 3 */ \
		{return;} \
	if(branch1 != 5) \
		{break;} \
	start_timer(29); \
	xgw_carries_wpn(g); \
	end_timer(29); \
}while(1);




#define pass1_forward_fft_setup(LABEL1)

#define pass1_inverse_fft_setup(LABEL1, procname) \
	g->next_block = 1; \
	start_timer(1); \
	g->thread_work_routine = procname##_i0 ; \
	pass1_wake_up_threads(g); \



#define newxloop_7arg(clm, label2, incr_rsi, pm_incr_reg, pm_incr_amt, sc_incr_reg, sc_incr_amt) \
	rsi += -2*clm*64+incr_rsi; /* Next source pointer */ \
	pm_incr_reg += pm_incr_amt; /* Next premultiplier pointer */ \
	sc_incr_reg += sc_incr_amt; /* Next sin/cos pointer */
#define newxloop_5arg(clm, label2, incr_rsi, pm_incr_reg, pm_incr_amt) \
	rsi += -2*clm*64+incr_rsi; \
	pm_incr_reg += pm_incr_amt;
#define newxloop_3arg(clm, label2, incr_rsi) \
	rsi += -2*clm*64+incr_rsi;
#define newxloop_2arg(clm, label2)


#define GET_MACRO(A,B,C,D,E,F,G,H, ...) H
#define newxloop_N(...)    GET_MACRO(__VA_ARGS__, newxloop_7arg, newxloop_6arg, newxloop_5arg, newxloop_4arg, newxloop_3arg, newxloop_2arg)(__VA_ARGS__)



#define newxloop1nopre(...) \
	newxloop_N(__VA_ARGS__)
#define newxloop1(...) \
	newxloop_N(__VA_ARGS__)
#define newxloop2(...) \
	newxloop_N(__VA_ARGS__)
#define newxloop4(...) \
	newxloop_N(__VA_ARGS__)
#define newxloop1doublepre(clm, label2, label1, ...) \
	newxloop1nopre(clm, label2, __VA_ARGS__)
#define newxloop1quadpre(clm, label2, label1, ...) \
	newxloop1nopre(clm, label2, __VA_ARGS__)




/* Copyright 2009-2023 - Mersenne Research, Inc.  All rights reserved */ \
/* Author:  George Woltman */ \
/* Email: woltman@alum.mit.edu */ \
/* */ \
/* These macros do the first pass of a primarily radix-4 DJB "r4dwpn" FFTs. */ \
/* This is like the radix-4 "r4delay" FFTs except that the group ttp multiplier */ \
/* is merged with the sin/cos twiddles in one of the FFT levels.  This saves */ \
/* us two multiply operations during normalization. */ \
/* */ \
\
/* ******************************************************** */ \
/* ******************************************************** */ \
/* ******************  PASS 1 MACROS  ********************* */ \
/* ******************************************************** */ \
/* ******************************************************** */ \
\
/* */ \
/* Estimated memory requirements for each real pass 1 macro: */ \
/* */ \
/* Pass 1 size:			128	256	320	384	448	512	640	768	896	1024	1280	1536	1792	2048	2560	3072	3584	4096 */ \
/* */ \
/* FFT data:			16K	32K	40K	48K	56K	64K	80K	96K	112K	128K	160K	192K	224K	256K	320K	384K	448K	512K */ \
/* Scratch size:		16K	32K	40K	48K	56K	64K	80K	96K	112K	128K	160K	192K	224K	256K	320K	384K	448K	512K */ \
/* Variable sin/cos data:	11.5K	23K	17.5K	19.5K	21.5K	23.5K	35K	39K	43K	47K	33.5K	31.5K	46.5K	48.5K	73K	63K	75K	79K */ \
/* Fixed sin/cos data:		1.5K	1.5K	4.5K	4.5K	6.5K	6K	9K	9K	13K	16K	18K	18K	26K	32K	36K	48K	52K	64K */ \
/* Normalization data:		5K	5K	12.5K	15K	17.5K	20K	12.5K	15K	17.5K	20K	25K	60K	17.5K	20K	50K	60K	70K	80K */ \
/* Big/lit flags:		0.5K	1K	1.25K	1.5K	1.75K	2K	2.5K	3K	3.5K	4K	5K	6K	7K	8K	10K	12K	14K	16K */ \
/* Data prefetched early:	8K	16K	16K	12K	16K	16K	32K	24K	32K	32K	16K	0K	64K	0K	32K	0K	128K	0K */ \
/* */ \
/* Min L2 cache size (clm=4):	86K	165.5K	189.25K	216K	252.75K	283K	366K	417K	488K	550K	611K	723K	879.5K	928K	1234K	1398K	1758K	1858K */ \
/*     L2 cache size (clm=2):		86K	103.13K	117.75K	138.38K	154.5K	193.75K	220.5K	259.25K	293K	327K	400.5K	461.5K	488.5K	660K	753K	940K	999K */ \
/*     L2 cache size (clm=1):		46.25K	60K	68.63K	81.25K	90.25K	107.63K	122.25K	144.88K	164.5K	185K	239.25K	252.5K	270.25K	373K	430.5K	531K	571.5K */ \
/* */ \
/* NOTE:  To calculate the minimum efficient L2 cache size, we double the */ \
/* data and variable sin/cos data sizes as both the current data */ \
/* and next block data will be in the cache at the same time.  Any data */ \
/* prefetched early (during the inverse FFT process) runs the danger of */ \
/* pushing out the fixed sin/cos data and normalization data because of the */ \
/* LRU cache-eviction algorithm.  The L2 cache size formula is: */ \
/*	(FFT data size + Variable sin/cos data) * 2 + */ \
/*	(Scratch size + Fixed sin/cos data + Normalization data + big/lit data) + */ \
/*	(amount of data prefetched during inverse FFT process) */ \
/* ALSO NOTE:  Numbers above are for clm = 4.  Clm affects the FFT data size, */ \
/* scratch size, variable sin/cos size, big/lit flags, and data prefetched early. */ \
\
/* Estimated memory requirements for each negacyclic pass 1 macro: */ \
/* */ \
/* Pass 1 size:			128	256	384	512	640	768	1024	1280	1536	2048	2560	3072	4096	5120 */ \
/* */ \
/* FFT data:			16K	32K	48K	64K	80K	96K	128K	160K	192K	256K	320K	384K	512K	640K */ \
/* Scratch size:		16K	32K	48K	64K	80K	96K	128K	160K	192K	256K	320K	384K	512K	640K */ \
/* Variable sin/cos data:	8K	12K	12K	16K	16K	24K	32K	32K	24K	34K	28K	48K	56K	56K */ \
/* Fixed sin/cos data:		5K	5K	9K	12K	15K	18K	24K	30K	39K	52K	65K	78K	104K	130K */ \
/* Normalization data:		5K	5K	15K	20K	25K	15K	20K	25K	60K	20K	100K	60K	80K	100K */ \
/* Big/lit flags:		0.5K	1K	1.5K	2K	2.5K	3K	4K	5K	6K	8K	10K	12K	16K	20K */ \
/* Data prefetched early:	8K	16K	12K	16K	10K	24K	32K	20K	0K	0K	0K	0K	0K	0K */ \
/* */ \
/* Min L2 cache size (clm=4):	82.5K	147K	205.5K	274K	324.5K	396K	528K	624K	729K	916K	1191K	1398K	1848K	2282K */ \
/*     L2 cache size (clm=2):	46.25K	78.5K	114.75K	153K	182.25K	214.5K	286K	339.5K	414K	494K	678K	768K	1016K	1060K */ \
/*     L2 cache size (clm=1):		44.25K	69.38K	92.5K	111.13K	123.75K	165K	197.25K	256.5K	283K	421.5K	453K	600K	645K */ \
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
/* Perform a pass 1 FFT of length 128 doubles. */ \
/* Memory examples are for a 128K FFT with 10 levels done in pass 2. */ \
/* Sin/cos data needed: 4*2*clm*XMM_SCND2 + 4*2*clm*XMM_SCND2R + 2*clm*XMM_SCD2 + 4*2*clm*XMM_SCD4 = 23*2*clm 64-byte cache lines */ \
\

#define r4dwpn_pass1sc128(procname, xpass2name, pass2_levels, clm) /* (fftname, pass2_levels, clmarg) */ \
\
/* set_FFT_constants(fftname, 0, 0, pass2_levels, clmarg); */ \
xfft_header_rn(procname, xpass2name, pass2_levels, clm); \
\
/* Do pass 2 (the last 10 FFT levels if this is a 128K FFT) */ \
\
/* pass2: */start_timer(0); \
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
/* Do 15 four_complex_unfft4 operations */ \
/*	distance between fft data elements is 1K */ \
/*	do 15 macros each processing 16 data values */ \
\
c0b:	start_timer(16); \
rsi = (uintptr_t)g->data_addr;/* Get source address */ \
rcx = (uintptr_t)g->scratch_area;	/* Get address of scratch area */ \
rdx = (uintptr_t)g->premult_addr;/* Real data sin/cos pointer */ \
rdi = (rdx+2*clm*XMM_SCD2); /* Complex data sin/cos ptr */ \
pfing(rbp = (rdx+(23-8)*2*clm*64)); /* Finish prefetch the sin/cos/premult data */ \
loops_init_prefetch_touch(16*2*clm, 64, 2, rbp, 4096, c9a); /* 16*2*clm iterations */ \
xtouch(xptr(rbp));			/* Prefetch TLB */ \
for(unsigned int loopA = 2*clm; loopA; loopA--) { \
	r4_sg2cl_eight_reals_unfft4(rsi, 64, blkdst, rcx, 64, clmblkdst, rdi, rdx); \
	rdi += XMM_SCD4;		/* Next premult/sin/cos data */ \
	rdx += XMM_SCD2;		/* Next sine/cosine pointer */ \
} \
rsi += -2*clm*64+2*blkdst;	/* Next source pointer */ \
rcx += -2*clm*64+2*clmblkdst; /* Next dest pointer */ \
rdi += -2*clm*XMM_SCD4;	/* Restore premult/sine/cosine pointer */ \
\
loops_reset; \
/* c9a: */ xtouch(xptr(rbp+4096));		/* Prefetch TLB */ \
for(unsigned int loopA = 4*4-1; loopA; ) { \
	do{ \
		for(unsigned int loopC = 2*clm; loopC; loopC--) { \
			r4_sg2cl_four_complex_unfft4(rsi, 64, blkdst, rcx, 64, clmblkdst, rdi); \
			rdi += XMM_SCD4;		/* Next premult/sin/cos data */ \
		} \
		rsi += -2*clm*64+2*blkdst;	/* Next source pointer */ \
		rcx += -2*clm*64+2*clmblkdst; /* Next dest pointer */ \
		rdi += -2*clm*XMM_SCD4;	/* Restore premult/sin/cos data pointer */ \
	}while(--loopA & 3); \
	rcx += -4*2*clmblkdst+clmblkdst8; \
	rdi += 2*clm*XMM_SCD4;	/* Next premult/sin/cos pointer */ \
} \
rsi = (uintptr_t)g->scratch_area;	/* Restore source pointer */ \
rdx = rdi;		/* Next real data sin/cos pointer */ \
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
rdi = (rdx+4*2*clm*XMM_SCND2R); /* Next complex data sin/cos ptr */ \
set_data_prefetch_ptrs;		/* Init prefetching pointers */ \
loops_init(8*2*clm);		/* 8*2*clm iterations */ \
\
r4_x4cl_wpn_eight_reals_unfft_preload; \
for(unsigned int loopA = 2*2*clm; loopA; ) { \
	xtouch(xptr(rcx+rbp+(clm-1)*128));	/* Preload TLB */ \
	do{ \
		r4_x4cl_wpn_eight_reals_unfft(rsi, 64, 2*clmblkdst, 4*clmblkdst, rdi, 2*2*clm*XMM_SCND2, rdx, 2*2*clm*XMM_SCND2R); \
		rdi += XMM_SCND2;		/* Next sine/cosine pointer */ \
		rdx += XMM_SCND2R;		/* Next sine/cosine pointer */ \
	}while(--loopA & (2*clm-1)); \
	newxloop1(clm, c6b, clmblkdst);	/* 2*clm iterations */ \
} \
rsi += -2*clmblkdst+clmblkdst8; /* Next source pointer */ \
rdi += -2*2*clm*XMM_SCND2;	/* Restore sine/cosine pointer */ \
\
loops_reset; \
r4_x4cl_wpn_four_complex_djbunfft_preload; \
for(unsigned int loopA = 3; loopA; loopA--) { \
	for(unsigned int loopB = 2*2*clm; loopB; ) { \
		xtouch(xptr(rcx+rbp+(clm-1)*128));	/* Preload TLB */ \
		do{ \
			r4_x4cl_wpn_four_complex_djbunfft(rsi, 64, 2*clmblkdst, 4*clmblkdst, rdi, 2*2*clm*XMM_SCND2); \
			rdi += XMM_SCND2;		/* Next clm's sin/cos data */ \
		}while(--loopB & (2*clm-1)); \
		newxloop1(clm, c7b, clmblkdst);	/* 2*clm iterations */ \
	} \
	rsi += -2*clmblkdst+clmblkdst8; /* Next source pointer */ \
	rdi += -2*2*clm*XMM_SCND2;	/* Restore sine/cosine pointer */ \
} \
rsi += -4*clmblkdst8;	/* Restore source pointer */ \
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
rdi = (uintptr_t)g->sincos2;/* Load fixed sin/cos data ptr */ \
loops_init(8*2*clm);		/* 8 iterations of 2*clm */ \
r4_x4cl_eight_reals_last_unfft_preload; \
for(unsigned int loopA = 8*2*clm; loopA; ) { \
	xtouch(xptr(rcx+rbp+(clm-1)*128));	/* Preload TLB */ \
	do{ \
		r4_x4cl_eight_reals_last_unfft(rsi, 64, clmblkdst8, 2*clmblkdst8, rdi, 8*XMM_SCD3); \
	}while(--loopA & (2*clm-1)); \
	newxloop1(clm, c1b, clmblkdst, rdi, XMM_SCD3); /* 2*clm iterations */ \
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
/* pass1: */pass1_forward_fft_setup(b0b); \
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
rdi = (uintptr_t)g->sincos2;/* Load fixed sin/cos postmultipliers */ \
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
			}while(--loopB & (2*clm-1)); \
			newxloop2(clm, b1b, clmblkdst, rdi, XMM_SCD3); /* 2*2*clm iters prefetching 2*clm cache lines */ \
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
rdx = (uintptr_t)g->premult_addr;/* Address of premultiplier/sin/cos table */ \
rdx += 2*clm*XMM_SCD2 + 4*2*clm*XMM_SCD4; /* Real data sin/cos ptr */ \
rdi = (rdx+4*2*clm*XMM_SCND2R); /* Complex data sin/cos ptr */ \
loops_init(8*2*clm);		/* 2 iterations of 4 iters of 2*clm */ \
r4_x4cl_wpn_eight_reals_four_complex_djbfft_preload; \
for(unsigned int loopA = 4*2*clm; loopA; ) { \
	xtouch(xptr(rcx+rbp+(clm-1)*128));	/* Preload TLB */ \
	do{ \
		r4_x4cl_wpn_eight_reals_four_complex_djbfft(rsi, 64, 4*clmblkdst, clmblkdst8, rdi, rdx); \
		rdi += XMM_SCND2;		/* Next sine/cosine pointer */ \
		rdx += XMM_SCND2R;		/* Next sine/cosine pointer */ \
	}while(--loopA & (2*clm-1)); \
	newxloop1(clm, b6b, clmblkdst);	/* 2*clm iterations */ \
} \
rsi += -4*clmblkdst+2*clmblkdst8; \
rdi += -4*2*clm*XMM_SCND2;	/* Restore sin/cos ptr */ \
\
loops_reset; \
r4_x4cl_wpn_four_complex_djbfft_preload; \
{ \
	for(unsigned int loopB = 4*2*clm; loopB; ) { \
		xtouch(xptr(rcx+rbp+(clm-1)*128));	/* Preload TLB */ \
		do{ \
			r4_x4cl_wpn_four_complex_djbfft(rsi, 64, 4*clmblkdst, clmblkdst8, rdi); \
			rdi += XMM_SCND2;		/* Next sine/cosine pointer */ \
		}while(--loopB & (2*clm-1)); \
		newxloop1(clm, b7b, clmblkdst);	/* 2*clm iterations */ \
	} \
	rsi += -4*clmblkdst+2*clmblkdst8; /* Next source pointer */ \
	rdi += -4*2*clm*XMM_SCND2;	/* Restore sin/cos ptr */ \
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
/* Do 1 four_complex_fft4 operations */ \
/*	distance between fft data elements is 1K */ \
/*	do 1 macros each processing 32 data values */ \
/* Do 14 four_complex_fft4 operations */ \
/*	distance between fft data elements is 1K */ \
/*	do 7 macros each processing 32 data values */ \
\
start_timer(26); \
rdi += -4*2*clm*XMM_SCND2R-4*2*clm*XMM_SCD4; /* Complex data sin/cos ptr */ \
rdx = (rdi-2*clm*XMM_SCD2); /* Real data sin/cos ptr */ \
rcx = (uintptr_t)g->data_addr;/* Destination pointer */ \
loops_init_prefetch(8*2*clm, 128, 1, rbp); /* 8 iters of 2*clm */ \
pfing(rbp = (uintptr_t)g->premult_prefetch);/* Load premult/sin/cos prefetch pointer */ \
r4_sg4cl_eight_reals_four_complex_fft4_preload; \
for(unsigned int loopA = 2*clm; loopA; loopA--) { \
	r4_sg4cl_eight_reals_four_complex_fft4(rsi, 64, clmblkdst, 2*clmblkdst, rcx, 64, blkdst, 2*blkdst, rdi, rdx); \
	rdi += XMM_SCD4;		/* Next sine/cosine pointer */ \
	rdx += XMM_SCD2;		/* Next sine/cosine pointer */ \
} \
rsi += -2*clm*64+4*clmblkdst; /* Next source pointer */ \
rcx += -2*clm*64+4*blkdst;	/* Next dest pointer */ \
rdi += -2*clm*XMM_SCD4;	/* Restore premult/sine/cosine pointer */ \
\
loops_reset; \
for(unsigned int loopA = 4*2-1; loopA; ) { \
	do{ \
		for(unsigned int loopC = 2*clm; loopC; loopC--) { \
			r4_sg4cl_four_complex_fft4(rsi, 64, clmblkdst, 2*clmblkdst, rcx, 64, blkdst, 2*blkdst, rdi); \
			rdi += XMM_SCD4;		/* Next premult/sin/cos ptr */ \
		} \
		rsi += -2*clm*64+4*clmblkdst; /* Next source pointer */ \
		rcx += -2*clm*64+4*blkdst;	/* Next dest pointer */ \
		rdi += -2*clm*XMM_SCD4;	/* Restore premult/sine/cosine pointer */ \
	}while(--loopA & 1); \
	rsi += -2*4*clmblkdst+clmblkdst8; /* Next source pointer */ \
	rdi += 2*clm*XMM_SCD4;	/* Next premult/sine/cosine pointer */ \
} \
end_timer(26); \
\
/* Work on next 4*clm sets of 128 values */ \
\
endlp:	xpass1_get_next_block(pass2, c0b, b0b, _wpn); /* Get next block and jump */ \
/* done: */ xfft_footer_rn; \



/* Perform a pass 1 FFT of length 256 doubles. */ \
/* Memory examples are for a 256K FFT with 10 levels done in pass 2. */ \
/* Sin/cos data needed: 8*2*clm*XMM_SCND2 + 8*2*clm*XMM_SCND2R + 2*clm*XMM_SCD4 + 4*2*clm*XMM_SCD8 = 46*2*clm 64-byte cache lines */ \


#define r4dwpn_pass1sc256(procname, xpass2name, pass2_levels, clm) /* (fftname, pass2_levels, clmarg) */ \
\
/* set_FFT_constants(fftname, 0, 0, pass2_levels, clmarg); */ \
xfft_header_rn(procname, xpass2name, pass2_levels, clm); \
\
/* Do pass 2 (the last 10 FFT levels if this is a 256K FFT) */ \
\
/* pass2: */start_timer(0); \
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
pfing(rbp = (rdx+(46-32)*2*clm*64)); /* Finish prefetching the sin/cos/premult data */ \
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
/* c9a: */ xtouch(xptr(rbp+4096));		/* Prefetch TLB */ \
for(unsigned int loopA = 4*2*2-1; loopA; ) { \
	do{ \
		do{ \
			for(unsigned int loopD = 2*clm; loopD; loopD--) { \
				r8_sg4cl_eight_complex_unfft8(rsi, 64, blkdst, 2*blkdst, rcx, 64, clmblkdst, 2*clmblkdst, rdi-128); \
				rdi += XMM_SCD8;		/* Next premult/sin/cos data */ \
			} \
			rsi += -2*clm*64+4*blkdst;	/* Next source pointer */ \
			rcx += -2*clm*64+4*clmblkdst; /* Next dest pointer */ \
			rdi += -2*clm*XMM_SCD8;	/* Restore premult/sin/cos data pointer */ \
		}while(--loopA & 1); \
		rcx += -2*4*clmblkdst+clmblkdst8; \
	}while(loopA & (2*2-1)); \
	rdi += 2*clm*XMM_SCD8;	/* Next premult/sin/cos pointer */ \
} \
rsi = (uintptr_t)g->scratch_area;	/* Restore source pointer */ \
rdx = (rdi-128);		/* Next real data sin/cos pointer */ \
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
rdi = (rdx+8*2*clm*XMM_SCND2R); /* Complex data sin/cos ptr */ \
set_data_prefetch_ptrs;		/* Init prefetching pointers */ \
loops_init(16*2*clm);		/* 16 iterations of 2*clm */ \
r4_x4cl_wpn_eight_reals_unfft_preload; \
for(unsigned int loopA = 4*2*clm; loopA; ) { \
	xtouch(xptr(rcx+rbp+(clm-1)*128));	/* Preload TLB */ \
	do{ \
		r4_x4cl_wpn_eight_reals_unfft(rsi, 64, 4*clmblkdst, clmblkdst8, rdi, 4*2*clm*XMM_SCND2, rdx, 4*2*clm*XMM_SCND2R); \
		rdi += XMM_SCND2;		/* Next sine/cosine pointer */ \
		rdx += XMM_SCND2R;		/* Next sine/cosine pointer */ \
	}while(--loopA & (2*clm-1)); \
	newxloop1(clm, c6b, clmblkdst);	/* 2*clm iterations */ \
} \
rsi += -4*clmblkdst+2*clmblkdst8; /* Next source pointer */ \
rdi += -4*2*clm*XMM_SCND2;	/* Restore sine/cosine pointer */ \
\
loops_reset; \
r4_x4cl_wpn_four_complex_djbunfft_preload; \
for(unsigned int loopA = 3; loopA; loopA--) { \
	for(unsigned int loopB = 4*2*clm; loopB; ) { \
		xtouch(xptr(rcx+rbp+(clm-1)*128));	/* Preload TLB */ \
		do{ \
			r4_x4cl_wpn_four_complex_djbunfft(rsi, 64, 4*clmblkdst, clmblkdst8, rdi, 4*2*clm*XMM_SCND2); \
			rdi += XMM_SCND2;		/* Next clm's sin/cos data */ \
		}while(--loopB & (2*clm-1)); \
		newxloop1(clm, c7b, clmblkdst);	/* 2*clm iterations */ \
	} \
	rsi += -4*clmblkdst+2*clmblkdst8; /* Next source pointer */ \
	rdi += -4*2*clm*XMM_SCND2;	/* Restore sine/cosine pointer */ \
} \
rsi += -4*2*clmblkdst8;	/* Restore source pointer */ \
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
rdi = (uintptr_t)g->sincos2;/* Load fixed sin/cos postmultipliers */ \
loops_init(16*2*clm);		/* 16 iterations of 2*clm */ \
r4_x4cl_eight_reals_last_unfft_preload; \
for(unsigned int loopA = 2; loopA; loopA--) { \
	for(unsigned int loopB = 8*2*clm; loopB; ) { \
		xtouch(xptr(rcx+rbp+(clm-1)*128));	/* Preload TLB */ \
		do{ \
			r4_x4cl_eight_reals_last_unfft(rsi, 64, 2*clmblkdst8, 4*clmblkdst8, rdi, 16*XMM_SCD3); \
		}while(--loopB & (2*clm-1)); \
		newxloop1(clm, c1b, clmblkdst, rdi, XMM_SCD3); /* 2*clm iterations */ \
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
/* pass1: */pass1_forward_fft_setup(b0b); \
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
rdi = (uintptr_t)g->sincos2;/* Load fixed sin/cos postmultipliers */ \
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
			}while(--loopB & (2*clm-1)); \
			newxloop2(clm, b1b, clmblkdst, rdi, XMM_SCD3); /* 2*2*clm iters prefetching 2*clm cache lines */ \
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
/*  BUG numbers in comment may be wrong */ \
start_timer(25); \
rdx = (uintptr_t)g->premult_addr;/* Address of premultiplier/sin/cos table */ \
rdx += 2*clm*XMM_SCD4 + 4*2*clm*XMM_SCD8; /* Real data sin/cos ptr */ \
rdi = (rdx+8*2*clm*XMM_SCND2R); /* Complex data sin/cos ptr */ \
loops_init(16*2*clm);		/* 2 iterations of 8 iters of 2*clm */ \
r4_x4cl_wpn_eight_reals_four_complex_djbfft_preload; \
for(unsigned int loopA = 8*2*clm; loopA; ) { \
	xtouch(xptr(rcx+rbp+(clm-1)*128));	/* Preload TLB */ \
	do{ \
		r4_x4cl_wpn_eight_reals_four_complex_djbfft(rsi, 64, clmblkdst8, 2*clmblkdst8, rdi, rdx); \
		rdi += XMM_SCND2;		/* Next sine/cosine pointer */ \
		rdx += XMM_SCND2R;		/* Next sine/cosine pointer */ \
	}while(--loopA & (2*clm-1)); \
	newxloop1(clm, b6b, clmblkdst);	/* 2*clm iterations */ \
} \
rsi += -8*clmblkdst+4*clmblkdst8; \
rdi += -8*2*clm*XMM_SCND2;	/* Restore sin/cos ptr */ \
\
loops_reset; \
r4_x4cl_wpn_four_complex_djbfft_preload; \
{ \
	for(unsigned int loopB = 8*2*clm; loopB; ) { \
		xtouch(xptr(rcx+rbp+(clm-1)*128));	/* Preload TLB */ \
		do{ \
			r4_x4cl_wpn_four_complex_djbfft(rsi, 64, clmblkdst8, 2*clmblkdst8, rdi); \
			rdi += XMM_SCND2;		/* Next sine/cosine pointer */ \
		}while(--loopB & (2*clm-1)); \
		newxloop1(clm, b7b, clmblkdst);	/* 2*clm iterations */ \
	} \
	rsi += -8*clmblkdst+4*clmblkdst8; /* Next source pointer */ \
	rdi += -8*2*clm*XMM_SCND2;	/* Restore sin/cos ptr */ \
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
/* Do 30 eight_complex_fft8 operations */ \
/*	distance between fft data elements is 1K */ \
/*	do 15 macros each processing 64 data values */ \
\
start_timer(26); \
rdi += -8*2*clm*XMM_SCND2R-4*2*clm*XMM_SCD8+128; /* Complex data sin/cos ptr and use -128 */ \
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
for(unsigned int loopA = 4*2-1; loopA; ) { \
	do{ \
		for(unsigned int loopC = 2*clm; loopC; loopC--) { \
			r8_sg8cl_eight_complex_fft8(rsi, 64, clmblkdst, 2*clmblkdst, 4*clmblkdst, rcx, 64, blkdst, 2*blkdst, 4*blkdst, rdi-128); \
			rdi += XMM_SCD8;		/* Next premult/sin/cos ptr */ \
		} \
		rsi += -2*clm*64+clmblkdst8; /* Next source pointer */ \
		rcx += -2*clm*64+8*blkdst;	/* Next dest pointer */ \
		rdi += -2*clm*XMM_SCD8;	/* Restore premult/sine/cosine pointer */ \
	}while(--loopA & 1); \
	rdi += 2*clm*XMM_SCD8;	/* Next premult/sine/cosine pointer */ \
} \
end_timer(26); \
\
/* Work on next 4*clm sets of 256 values */ \
\
endlp:	xpass1_get_next_block(pass2, c0b, b0b, _wpn); /* Get next block and jump */ \
/* done: */ xfft_footer_rn; \
\
\
\
/* Perform a pass 1 FFT of length 320 doubles. */ \
/* Memory examples are for a 320K FFT with 10 levels done in pass 2. */ \
/* Sin/cos data needed: 4*2*clm*XMM_SCND2 + 4*2*clm*XMM_SCND2R + 2*clm*XMM_SCD2 + 10*2*clm*XMM_SCD4 = 35*2*clm 64-byte cache lines */ \
\

#define r4dwpn_pass1sc320(procname, xpass2name, pass2_levels, clm) /* (fftname, pass2_levels, clmarg) */ \
\
/* set_FFT_constants(fftname, 0, 0, pass2_levels, clmarg); */ \
xfft_header_rn(procname, xpass2name, pass2_levels, clm); \
\
/* Do pass 2 (the last 10 FFT levels if this is a 320K FFT) */ \
\
/* pass2: */start_timer(0); \
xpass2name(g); \
end_timer(0); \
if(g->ffttype == 1)		/* We're done if FFTing only */ \
	return; /* goto done; */ \
\
/* Do the inverse FFT */ \
\
pass1_inverse_fft_setup(c0b, procname); \
\
/* Do inverse FFT levels 8,9 */ \
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
/* Do 39 four_complex_unfft4 operations */ \
/*	distance between fft data elements is 1K */ \
/*	do 39 macros each processing 16 data values */ \
\
c0b:	start_timer(16); \
rsi = (uintptr_t)g->data_addr;/* Get source address */ \
rcx = (uintptr_t)g->scratch_area;	/* Get address of scratch area */ \
rdx = (uintptr_t)g->premult_addr;/* Real data sin/cos pointer */ \
rdi = (rdx+2*clm*XMM_SCD2); /* Complex data sin/cos ptr */ \
pfing(rbp = (rdx+(35-20)*2*clm*64)); /* Finish prefetching the sin/cos data */ \
loops_init_prefetch(40*2*clm, 64, 2, rbp); /* 40*2*clm iterations */ \
/*bug	xtouchpages rbp, 40*2*clm * 64 / 2 ;; Prefetch TLBs */ \
for(unsigned int loopA = 2*clm; loopA; loopA--) { \
	r4_sg2cl_eight_reals_unfft4(rsi, 64, blkdst, rcx, 64, clmblkdst, rdi, rdx); \
	rdi += XMM_SCD4;		/* Next premult/sin/cos data */ \
	rdx += XMM_SCD2;		/* Next sine/cosine pointer */ \
} \
rsi += -2*clm*64+2*blkdst;	/* Next source pointer */ \
rcx += -2*clm*64+2*clmblkdst; /* Next dest pointer */ \
rdi += -2*clm*XMM_SCD4;	/* Restore premult/sine/cosine pointer */ \
\
loops_reset; \
for(unsigned int loopA = 10*4-1; loopA; ) { \
	do{ \
		for(unsigned int loopC = 2*clm; loopC; loopC--) { \
			r4_sg2cl_four_complex_unfft4(rsi, 64, blkdst, rcx, 64, clmblkdst, rdi); \
			rdi += XMM_SCD4;		/* Next premult/sin/cos data */ \
		} \
		rsi += -2*clm*64+2*blkdst;	/* Next source pointer */ \
		rcx += -2*clm*64+2*clmblkdst; /* Next dest pointer */ \
		rdi += -2*clm*XMM_SCD4;	/* Restore premult/sin/cos data pointer */ \
	}while(--loopA & 3); \
	rcx += -4*2*clmblkdst+clmblkdst8; \
	rdi += 2*clm*XMM_SCD4;	/* Next premult/sin/cos pointer */ \
} \
rsi = (uintptr_t)g->scratch_area;	/* Restore source pointer */ \
rdx = rdi;		/* Next real data sin/cos pointer */ \
end_timer(16); \
\
/* Do inverse FFT levels 6,7 */ \
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
/* Do 36 four_complex_unfft operations */ \
/*	distance between fft data elements is 4K */ \
/*	do 18 macros each processing 32 data values */ \
\
start_timer(17); \
rdi = (rdx+4*2*clm*XMM_SCND2R); /* Next sin/cos ptr */ \
loops_init(20*2*clm);		/* 20*2*clm iterations */ \
\
r4_x4cl_wpn_eight_reals_unfft_preload; \
for(unsigned int loopA = 2; loopA; loopA--) { \
	for(unsigned int loopB = 2*clm; loopB; loopB--) { \
		r4_x4cl_wpn_eight_reals_unfft(rsi, 64, 2*clmblkdst, 4*clmblkdst, rdi, 2*2*clm*XMM_SCND2, rdx, 2*2*clm*XMM_SCND2R); \
		rdi += XMM_SCND2;		/* Next sine/cosine pointer */ \
		rdx += XMM_SCND2R;		/* Next sine/cosine pointer */ \
	} \
	rsi += -2*clm*64+clmblkdst; /* Next source pointer */ \
} \
rsi += -2*clmblkdst+clmblkdst8; /* Next source pointer */ \
rdi += -2*2*clm*XMM_SCND2;	/* Restore sine/cosine pointer */ \
\
loops_reset; \
r4_x4cl_wpn_four_complex_djbunfft_preload; \
for(unsigned int loopA = 9; loopA; loopA--) { \
	for(unsigned int loopB = 2; loopB; loopB--) { \
		for(unsigned int loopC = 2*clm; loopC; loopC--) { \
			r4_x4cl_wpn_four_complex_djbunfft(rsi, 64, 2*clmblkdst, 4*clmblkdst, rdi, 2*2*clm*XMM_SCND2); \
			rdi += XMM_SCND2;		/* Next clm's sin/cos data */ \
		} \
		rsi += -2*clm*64+clmblkdst; /* Next source pointer */ \
	} \
	rsi += -2*clmblkdst+clmblkdst8; /* Next source pointer */ \
	rdi += -2*2*clm*XMM_SCND2;	/* Restore sine/cosine pointer */ \
} \
rsi += -10*clmblkdst8;	/* Restore source pointer */ \
end_timer(17); \
\
/* Do inverse FFT levels 1,2,3,4,5 (actually 4.32 levels) */ \
/* On input the 64-byte cache lines hold these data values: */ \
/*	0K	8K	16K	24K */ \
/*	1K	... */ \
/*	... */ \
/*	7K */ \
/*	32K	... */ \
/*	... */ \
/* On output the 64-byte cache lines hold these data values: */ \
/*	0K	80K	160K	240K */ \
/*	1K	... */ \
/*	... */ \
/*	79K */ \
\
/* Do 64 20_reals_last_unfft operations */ \
/*	distance between fft data elements is 16K */ \
/*	do 8 macros each processing 80 data values */ \
\
start_timer(20); \
rdi = (uintptr_t)g->sincos2;/* Load fixed sin/cos postmultipliers */ \
set_data_prefetch_ptrs;		/* Init prefetching pointers */ \
loops_init(8*2*clm);		/* 8 iterations of 2*clm */ \
r5_x10cl_20_reals_last_unfft_preload; \
for(unsigned int loopA = 8*2*clm; loopA; ) { \
	xtouch(xptr(rcx+rbp+(clm-1)*128));	/* Preload TLB */ \
	do{ \
		r5_x10cl_20_reals_last_unfft(rsi, 64, clmblkdst8, rdi, 8*XMM_SCD9); \
	}while(--loopA & (2*clm-1)); \
	newxloop1quadpre(clm, c1b, c1a, clmblkdst, rdi, XMM_SCD9); /* 2*clm iterations */ \
} \
end_timer(20); \
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
/* pass1: */pass1_forward_fft_setup(b0b); \
\
/* Copy data from the FFT data area back to the scratch area */ \
\
b0b:	copy_fft_data_to_scratch(80, clm); \
\
/* Do FFT levels 1,2,3,4,5 (actually 4.32 levels) */ \
/* */ \
/* On input the 64-byte cache lines hold these data values: */ \
/*	0K	80K	160K	240K */ \
/*	1K	... */ \
/*	... */ \
/*	79K */ \
/* On output the 64-byte cache lines hold these data values: */ \
/*	0K	16K	32K	48K */ \
/*	1K	... */ \
/*	... */ \
/*	15K	... */ \
/*	64K	... */ \
/*	... */ \
\
/* Do 64 20_reals_first_fft operations */ \
/*	distance between fft data elements is 16K */ \
/*	do 16 iterations each processing 40 data values */ \
\
yes_fft:start_timer(22); \
rsi = (uintptr_t)g->scratch_area;	/* Get address of scratch area */ \
rdi = (uintptr_t)g->sincos2;/* Load fixed sin/cos postmultipliers */ \
set_data_prefetch_ptrs;		/* Init prefetching pointers */ \
pfing(rcx += (32-4)*blkdst);	/* We've already prefetched 32 blocks, but */ \
;					/* code below will fetch 4 blocks too many */ \
loops_init(16*2*clm);		/* 16 iters of 2*clm */ \
r5_x5cl_20_reals_first_fft_scratch_preload; \
for(unsigned int loopA = 2; loopA; loopA--) { \
	for(unsigned int loopB = 8*2*clm; loopB; ) { \
		xtouch(xptr(rcx+rbp+(clm-1)*128));	/* Preload TLB */ \
		do{ \
			r5_x5cl_20_reals_first_fft_scratch(rsi, 64, 2*clmblkdst8, rdi); \
		}while(--loopB & (2*clm-1)); \
		newxloop1doublepre(clm, b1b, b1a, clmblkdst, rdi, XMM_SCD9); /* 2*clm iters prefetching 4*clm cache lines */ \
	} \
	rsi += -8*clmblkdst+clmblkdst8; \
} \
rsi += -2*clmblkdst8;	/* Restore source pointer */ \
end_timer(22); \
\
/* Do FFT levels 6,7 */ \
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
/* Do 32 four_complex_fft operations */ \
/*	distance between fft data elements is 4K */ \
/*	do 16 macros each processing 32 data values */ \
\
start_timer(25); \
rdx = (uintptr_t)g->premult_addr;/* Address of premultiplier/sin/cos table */ \
rdx += 2*clm*XMM_SCD2 + 10*2*clm*XMM_SCD4; /* Real data sin/cos ptr */ \
rdi = (rdx+4*2*clm*XMM_SCND2R); /* Complex data sin/cos ptr */ \
loops_init(20*2*clm);		/* 5 iterations of 4 iters of 2*clm */ \
r4_x4cl_wpn_eight_reals_four_complex_djbfft_preload; \
for(unsigned int loopA = 4*2*clm; loopA; ) { \
	xtouch(xptr(rcx+rbp+(clm-1)*128));	/* Preload TLB */ \
	do{ \
		r4_x4cl_wpn_eight_reals_four_complex_djbfft(rsi, 64, 4*clmblkdst, clmblkdst8, rdi, rdx); \
		rdi += XMM_SCND2;		/* Next sine/cosine pointer */ \
		rdx += XMM_SCND2R;		/* Next sine/cosine pointer */ \
	}while(--loopA & (2*clm-1)); \
	newxloop1(clm, b6b, clmblkdst);	/* 2*clm iterations */ \
} \
rsi += -4*clmblkdst+2*clmblkdst8; \
rdi += -4*2*clm*XMM_SCND2;	/* Restore sin/cos ptr */ \
\
loops_reset; \
r4_x4cl_wpn_four_complex_djbfft_preload; \
for(unsigned int loopA = 4; loopA; loopA--) { \
	for(unsigned int loopB = 4*2*clm; loopB; ) { \
		xtouch(xptr(rcx+rbp+(clm-1)*128));	/* Preload TLB */ \
		do{ \
			r4_x4cl_wpn_four_complex_djbfft(rsi, 64, 4*clmblkdst, clmblkdst8, rdi); \
			rdi += XMM_SCND2;		/* Next sine/cosine pointer */ \
		}while(--loopB & (2*clm-1)); \
		newxloop1(clm, b7b, clmblkdst);	/* 2*clm iterations */ \
	} \
	rsi += -4*clmblkdst+2*clmblkdst8; /* Next source pointer */ \
	rdi += -4*2*clm*XMM_SCND2;	/* Restore sin/cos ptr */ \
} \
rsi += -5*2*clmblkdst8;	/* Restore source pointer */ \
end_timer(25); \
\
/* Do FFT levels 8,9 */ \
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
/* Do 1 four_complex_fft4 operations */ \
/*	distance between fft data elements is 1K */ \
/*	do 1 macros each processing 32 data values */ \
/* Do 38 four_complex_fft4 operations */ \
/*	distance between fft data elements is 1K */ \
/*	do 19 macros each processing 32 data values */ \
\
start_timer(26); \
rdi += -4*2*clm*XMM_SCND2R-10*2*clm*XMM_SCD4; /* Complex data sin/cos ptr */ \
rdx = (rdi-2*clm*XMM_SCD2); /* Real data sin/cos ptr */ \
rcx = (uintptr_t)g->data_addr;/* Destination pointer */ \
loops_init_prefetch(20*2*clm, 64, 1, rbp); /* 20 iters of 2*clm */ \
pfing(rbp = (uintptr_t)g->premult_prefetch);/* Load premult/sin/cos prefetch pointer */ \
r4_sg4cl_eight_reals_four_complex_fft4_preload; \
for(unsigned int loopA = 2*clm; loopA; loopA--) { \
	r4_sg4cl_eight_reals_four_complex_fft4(rsi, 64, clmblkdst, 2*clmblkdst, rcx, 64, blkdst, 2*blkdst, rdi, rdx); \
	rdi += XMM_SCD4;		/* Next sine/cosine pointer */ \
	rdx += XMM_SCD2;		/* Next sine/cosine pointer */ \
} \
rsi += -2*clm*64+4*clmblkdst; /* Next source pointer */ \
rcx += -2*clm*64+4*blkdst;	/* Next dest pointer */ \
rdi += -2*clm*XMM_SCD4;	/* Restore premult/sine/cosine pointer */ \
\
loops_reset; \
for(unsigned int loopA = 10*2-1; loopA; ) { \
	do{ \
		for(unsigned int loopC = 2*clm; loopC; loopC--) { \
			r4_sg4cl_four_complex_fft4(rsi, 64, clmblkdst, 2*clmblkdst, rcx, 64, blkdst, 2*blkdst, rdi); \
			rdi += XMM_SCD4;		/* Next premult/sin/cos ptr */ \
		} \
		rsi += -2*clm*64+4*clmblkdst; /* Next source pointer */ \
		rcx += -2*clm*64+4*blkdst;	/* Next dest pointer */ \
		rdi += -2*clm*XMM_SCD4;	/* Restore premult/sine/cosine pointer */ \
	}while(--loopA & 1); \
	rsi += -2*4*clmblkdst+clmblkdst8; /* Next source pointer */ \
	rdi += 2*clm*XMM_SCD4;	/* Next premult/sine/cosine pointer */ \
} \
end_timer(26); \
\
/* Work on next 4*clm sets of 320 values */ \
\
endlp:	xpass1_get_next_block(pass2, c0b, b0b, _wpn); /* Get next block and jump */ \
/* done: */ xfft_footer_rn; \
\
\
\
/* Perform a pass 1 FFT of length 384 doubles. */ \
/* Memory examples are for a 384K FFT with 10 levels done in pass 2. */ \
/* Sin/cos data needed: 16*2*clm*XMM_SCD1 + 16*2*clm*XMM_SCD1 + 4*2*clm*XMM_SCND2 + 4*2*clm*XMM_SCND2R + */ \
/*			2*clm*XMM_SCD2 + 4*2*clm*XMM_SCD4 = 39*2*clm 64-byte cache lines */ \
\

#define r4dwpn_pass1sc384(procname, xpass2name, pass2_levels, clm) /* (fftname, pass2_levels, clmarg) */ \
\
/* set_FFT_constants(fftname, 0, 0, pass2_levels, clmarg); */ \
xfft_header_rn(procname, xpass2name, pass2_levels, clm); \
\
/* Do pass 2 (the last 10 FFT levels if this is a 384K FFT) */ \
\
/* pass2: */start_timer(0); \
xpass2name(g); \
end_timer(0); \
if(g->ffttype == 1)		/* We're done if FFTing only */ \
	return; /* goto done; */ \
\
/* Do the inverse FFT */ \
\
pass1_inverse_fft_setup(c0b, procname); \
\
/* Do inverse FFT levels 8,9 */ \
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
/* Do 47 four_complex_unfft4 operations */ \
/*	distance between fft data elements is 1K */ \
/*	do 47 macros each processing 16 data values */ \
\
c0b:	start_timer(16); \
rsi = (uintptr_t)g->data_addr;/* Get source address */ \
rcx = (uintptr_t)g->scratch_area;	/* Get address of scratch area */ \
rdx = (uintptr_t)g->premult_addr;/* Real data sin/cos pointer */ \
rdi = (rdx+2*clm*XMM_SCD2); /* Complex data sin/cos ptr */ \
pfing(rbp = (rdx+(39-24)*2*clm*64)); /* Finish prefetching the sin/cos data */ \
loops_init_prefetch(48*2*clm, 64, 2, rbp); /* 48*2*clm iterations */ \
/*bug	xtouchpages rbp, 48*2*clm * 64 / 2 ;; Prefetch TLBs */ \
for(unsigned int loopA = 2*clm; loopA; loopA--) { \
	r4_sg2cl_eight_reals_unfft4(rsi, 64, blkdst, rcx, 64, clmblkdst, rdi, rdx); \
	rdi += XMM_SCD4;		/* Next premult/sin/cos data */ \
	rdx += XMM_SCD2;		/* Next sine/cosine pointer */ \
} \
rsi += -2*clm*64+2*blkdst;	/* Next source pointer */ \
rcx += -2*clm*64+2*clmblkdst; /* Next dest pointer */ \
rdi += -2*clm*XMM_SCD4;	/* Restore premult/sine/cosine pointer */ \
\
loops_reset; \
for(unsigned int loopA = 4*3*4-1; loopA; ) { \
	for(unsigned int loopB = 3; loopB; loopB--) { \
		do{ \
			for(unsigned int loopD = 2*clm; loopD; loopD--) { \
				r4_sg2cl_four_complex_unfft4(rsi, 64, blkdst, rcx, 64, clmblkdst, rdi); \
				rdi += XMM_SCD4;		/* Next premult/sin/cos data */ \
			} \
			rsi += -2*clm*64+2*blkdst;	/* Next source pointer */ \
			rcx += -2*clm*64+2*clmblkdst; /* Next dest pointer */ \
			rdi += -2*clm*XMM_SCD4;	/* Restore premult/sin/cos data pointer */ \
		}while(--loopA & 3); \
		rcx += -4*2*clmblkdst+clmblkdst8; \
	} \
	rdi += 2*clm*XMM_SCD4;	/* Next premult/sin/cos pointer */ \
} \
rsi = (uintptr_t)g->scratch_area;	/* Restore source pointer */ \
rdx = rdi;		/* Next real data sin/cos pointer */ \
end_timer(16); \
\
/* Do inverse FFT levels 6,7 */ \
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
/* Do 44 four_complex_unfft operations */ \
/*	distance between fft data elements is 4K */ \
/*	do 22 macros each processing 32 data values */ \
\
start_timer(17); \
rdi = (rdx+4*2*clm*XMM_SCND2R); /* Complex data sin/cos ptr */ \
loops_init(24*2*clm);		/* 24*2*clm iterations */ \
\
r4_x4cl_wpn_eight_reals_unfft_preload; \
for(unsigned int loopA = 2; loopA; loopA--) { \
	for(unsigned int loopB = 2*clm; loopB; loopB--) { \
		r4_x4cl_wpn_eight_reals_unfft(rsi, 64, 2*clmblkdst, 4*clmblkdst, rdi, 2*2*clm*XMM_SCND2, rdx, 2*2*clm*XMM_SCND2R); \
		rdi += XMM_SCND2;		/* Next sine/cosine pointer */ \
		rdx += XMM_SCND2R;		/* Next sine/cosine pointer */ \
	} \
	rsi += -2*clm*64+clmblkdst; /* Next source pointer */ \
} \
rsi += -2*clmblkdst+clmblkdst8; /* Next source pointer */ \
rdi += -2*2*clm*XMM_SCND2;	/* Restore sine/cosine pointer */ \
\
loops_reset; \
r4_x4cl_wpn_four_complex_djbunfft_preload; \
for(unsigned int loopA = 11; loopA; loopA--) { \
	for(unsigned int loopB = 2; loopB; loopB--) { \
		for(unsigned int loopC = 2*clm; loopC; loopC--) { \
			r4_x4cl_wpn_four_complex_djbunfft(rsi, 64, 2*clmblkdst, 4*clmblkdst, rdi, 2*2*clm*XMM_SCND2); \
			rdi += XMM_SCND2;		/* Next clm's sin/cos data */ \
		} \
		rsi += -2*clm*64+clmblkdst; /* Next source pointer */ \
	} \
	rsi += -2*clmblkdst+clmblkdst8; /* Next source pointer */ \
	rdi += -2*2*clm*XMM_SCND2;	/* Restore sine/cosine pointer */ \
} \
rsi += -12*clmblkdst8;	/* Restore source pointer */ \
rdx = (rdi+4*2*clm*XMM_SCND2); /* Next real data sin/cos pointer */ \
end_timer(17); \
\
/* Do inverse FFT levels 4,5 (actually 1.585 levels) */ \
/* On input the 64-byte cache lines hold these data values: */ \
/*	0K	8K	16K	24K */ \
/*	1K	... */ \
/*	... */ \
/*	7K */ \
/*	32K	... */ \
/*	... */ \
/* On output the 64-byte cache lines hold these data values: */ \
/*	0K	24K	48K	72K */ \
/*	1K	... */ \
/*	... */ \
/*	23K */ \
/*	96K	... */ \
/*	... */ \
\
/* Do 16 six_reals_unfft operations */ \
/*	distance between fft data elements is 16K */ \
/*	do 8 macros each processing 24 data values */ \
/* Do 48 three_complex_unfft operations */ \
/*	distance between fft data elements is 16K */ \
/*	do 24 macros each processing 24 data values */ \
\
start_timer(18); \
rdi = (rdx+16*2*clm*XMM_SCD1); /* Complex data sin/cos ptr */ \
loops_init(32*2*clm);		/* 32 iterations of 2*clm */ \
for(unsigned int loopA = 8; loopA; loopA--) { \
	for(unsigned int loopB = 2*clm; loopB; loopB--) { \
		r3_x3cl_six_reals_unfft(rsi, 64, clmblkdst8, rdi, 8*2*clm*XMM_SCD1, rdx, 8*2*clm*XMM_SCD1); \
		rdi += XMM_SCD1;		/* Next sine/cosine pointer */ \
		rdx += XMM_SCD1;		/* Next sine/cosine pointer */ \
	} \
	rsi += -2*clm*64+clmblkdst; /* Next source pointer */ \
} \
rsi += -8*clmblkdst+3*clmblkdst8; /* Next source pointer */ \
rdi += -8*2*clm*XMM_SCD1;	/* Restore sine/cosine pointer */ \
\
loops_reset; \
r3_x3cl_three_complex_djbunfft_preload; \
for(unsigned int loopA = 3; loopA; loopA--) { \
	for(unsigned int loopB = 8; loopB; loopB--) { \
		for(unsigned int loopC = 2*clm; loopC; loopC--) { \
			r3_x3cl_three_complex_djbunfft(rsi, 64, clmblkdst8, rdi, 8*2*clm*XMM_SCD1); \
			rdi += XMM_SCD1;		/* Next clm's sin/cos data */ \
		} \
		rsi += -2*clm*64+clmblkdst; /* Next source pointer */ \
	} \
	rsi += -8*clmblkdst+3*clmblkdst8; /* Next source pointer */ \
	rdi += -8*2*clm*XMM_SCD1;	/* Restore sine/cosine pointer */ \
} \
rsi += -4*3*clmblkdst8;	/* Restore source pointer */ \
end_timer(18); \
\
/* Do inverse FFT levels 1,2,3 */ \
/* On input the 64-byte cache lines hold these data values: */ \
/*	0K	24K	48K	72K */ \
/*	1K	... */ \
/*	... */ \
/*	23K */ \
/*	96K	... */ \
/*	... */ \
/* On output the 64-byte cache lines hold these data values: */ \
/*	0K	96K	192K	288K */ \
/*	1K	... */ \
/*	... */ \
/*	95K */ \
\
/* Do 48 eight_reals_last_unfft operations */ \
/*	distance between fft data elements is 48K */ \
/*	do 24 macros each processing 32 data values */ \
\
start_timer(20); \
rdi = (uintptr_t)g->sincos2;/* Load fixed sin/cos postmultipliers */ \
set_data_prefetch_ptrs;		/* Init prefetching pointers */ \
loops_init(24*2*clm);		/* 24 iterations of 2*clm */ \
r4_x4cl_eight_reals_last_unfft_preload; \
for(unsigned int loopA = 3; loopA; loopA--) { \
	for(unsigned int loopB = 8*2*clm; loopB; ) { \
		xtouch(xptr(rcx+rbp+(clm-1)*128));	/* Preload TLB */ \
		do{ \
			r4_x4cl_eight_reals_last_unfft(rsi, 64, 3*clmblkdst8, 6*clmblkdst8, rdi, 24*XMM_SCD3); \
		}while(--loopB & (2*clm-1)); \
		newxloop1(clm, c1b, clmblkdst, rdi, XMM_SCD3); /* 2*clm iterations */ \
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
copy_scratch_data_to_fft(96, clm); \
goto	endlp;			/* Skip FFT code if flag not set */ \
\
/* Fire up auxiliary threads to do pass 1 forward FFT */ \
\
/* pass1: */pass1_forward_fft_setup(b0b); \
\
/* Copy data from the FFT data area back to the scratch area */ \
\
b0b:	copy_fft_data_to_scratch(96, clm); \
\
/* Do FFT levels 1,2,3 */ \
/* */ \
/* On input the 64-byte cache lines hold these data values: */ \
/*	0K	96K	192K	288K */ \
/*	1K	... */ \
/*	... */ \
/*	95K */ \
/* On output the 64-byte cache lines hold these data values: */ \
/*	0K	48K	96K	144K */ \
/*	1K	... */ \
/*	... */ \
/*	47K	... */ \
/*	192K	... */ \
/*	... */ \
\
/* Do 48 eight_reals_first_fft operations */ \
/*	distance between fft data elements is 48K */ \
/*	do 48 iterations each processing 16 data values */ \
\
yes_fft:start_timer(22); \
rsi = (uintptr_t)g->scratch_area;	/* Get address of scratch area */ \
rdi = (uintptr_t)g->sincos2;/* Load fixed sin/cos postmultipliers */ \
set_data_prefetch_ptrs;		/* Init prefetching pointers */ \
pfing(rcx += (24-4)*blkdst);	/* We've already prefetched 24 blocks, but */ \
;					/* code below prefetches 4 blocks too many */ \
loops_init(48*2*clm);		/* 6 iters of 4 iters of 2*2*clm */ \
r4_x2cl_eight_reals_first_fft_scratch_preload; \
for(unsigned int loopA = 6; loopA; loopA--) { \
	for(unsigned int loopB = 4*2*2*clm; loopB; ) { \
		xtouch(xptr(rcx+rbp+(clm-1)*128));	/* Preload TLB */ \
		do{ \
			do{ \
				r4_x2cl_eight_reals_first_fft_scratch(rsi, 64, 6*clmblkdst8, rdi); \
			}while(--loopB & (2*clm-1)); \
			newxloop2(clm, b1b, clmblkdst, rdi, XMM_SCD3); /* 2*2*clm iters prefetching 2*clm cache lines */ \
		}while(loopB & (2*2*clm-1)); \
	} \
	rsi += -4*2*clmblkdst+clmblkdst8; \
} \
rsi += -6*clmblkdst8;	/* Restore source pointer */ \
end_timer(22); \
\
/* Do FFT levels 4,5 (actually 1.585 levels) */ \
/* */ \
/* On input the 64-byte cache lines hold these data values: */ \
/*	0K	48K	96K	144K */ \
/*	1K	... */ \
/*	... */ \
/*	47K	... */ \
/*	192K	... */ \
/*	... */ \
/* On output the 64-byte cache lines hold these data values: */ \
/*	0K	16K	32K	48K */ \
/*	1K	... */ \
/*	... */ \
/*	15K	... */ \
/*	64K	... */ \
/*	... */ \
\
/* Do 16 six_reals_fft operations */ \
/* Do 16 three_complex_fft operations */ \
/*	distance between fft data elements is 16K */ \
/*	do 16 macros each processing 24 data values */ \
/* Do 32 three_complex_fft operations */ \
/*	distance between fft data elements is 16K */ \
/*	do 16 macros each processing 24 data values */ \
\
start_timer(24); \
rdx = (uintptr_t)g->premult_addr;/* Address of premultiplier/sin/cos table */ \
rdx += 4*2*clm*XMM_SCND2 + 4*2*clm*XMM_SCND2R + 2*clm*XMM_SCD2 + 4*2*clm*XMM_SCD4; /* Real data sin/cos ptr */ \
rdi = (rdx+16*2*clm*XMM_SCD1); /* Complex data sin/cos ptr */ \
loops_init(32*2*clm);		/* 2 sets of 2 iters of 8 of 2*clm */ \
r3_x3cl_six_reals_three_complex_djbfft_preload; \
for(unsigned int loopA = 2; loopA; loopA--) { \
	for(unsigned int loopB = 8*2*clm; loopB; ) { \
		xtouch(xptr(rcx+rbp+(clm-1)*128));	/* Preload TLB */ \
		do{ \
			r3_x3cl_six_reals_three_complex_djbfft(rsi, 64, 2*clmblkdst8, rdi, rdx); \
			rdi += XMM_SCD1;		/* Next sine/cosine pointer */ \
			rdx += XMM_SCD1;		/* Next sine/cosine pointer */ \
		}while(--loopB & (2*clm-1)); \
		newxloop1(clm, b4b, clmblkdst);	/* 2*clm iterations */ \
	} \
	rsi += -8*clmblkdst+clmblkdst8; \
} \
rsi += -2*clmblkdst8+6*clmblkdst8; /* Next source pointer */ \
rdi += -16*2*clm*XMM_SCD1;	/* Restore sin/cos ptr */ \
\
loops_reset; /* TODO: check if loops are correct. */ \
r3_x3cl_three_complex_djbfft_preload; \
for(unsigned int loopA = 2; loopA; loopA--) { \
	for(unsigned int loopB = 8*2*clm; loopB; ) { \
		xtouch(xptr(rcx+rbp+(clm-1)*128));	/* Preload TLB */ \
		do{ \
			r3_x3cl_three_complex_djbfft(rsi, 64, 2*clmblkdst8, rdi); \
			rdi += XMM_SCD1;		/* Next sine/cosine pointer */ \
		}while(--loopB & (2*clm-1)); \
		newxloop1(clm, b5b, clmblkdst);	/* 2*clm iterations */ \
	} \
	rsi += -8*clmblkdst+clmblkdst8; \
} \
rsi += -2*clmblkdst8-6*clmblkdst8; /* Restore source pointer */ \
end_timer(24); \
\
/* Do FFT levels 6,7 */ \
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
/* Do 40 four_complex_fft operations */ \
/*	distance between fft data elements is 4K */ \
/*	do 20 macros each processing 32 data values */ \
\
start_timer(25); \
rdi += -16*2*clm*XMM_SCD1-16*2*clm*XMM_SCD1-4*2*clm*XMM_SCND2; /* Complex data sin/cos ptr */ \
rdx = (rdi-4*2*clm*XMM_SCND2R); /* Real data sin/cos ptr */ \
loops_init(24*2*clm);		/* 6 iterations of 4 iters of 2*clm */ \
r4_x4cl_wpn_eight_reals_four_complex_djbfft_preload; \
for(unsigned int loopA = 4*2*clm; loopA; ) { \
	xtouch(xptr(rcx+rbp+(clm-1)*128));	/* Preload TLB */ \
	do{ \
		r4_x4cl_wpn_eight_reals_four_complex_djbfft(rsi, 64, 4*clmblkdst, clmblkdst8, rdi, rdx); \
		rdi += XMM_SCND2;		/* Next sine/cosine pointer */ \
		rdx += XMM_SCND2R;		/* Next sine/cosine pointer */ \
	}while(--loopA & (2*clm-1)); \
	newxloop1nopre(clm, b6b, clmblkdst);	/* 2*clm iterations */ \
} \
rsi += -4*clmblkdst+2*clmblkdst8; \
rdi += -4*2*clm*XMM_SCND2;	/* Restore sin/cos ptr */ \
\
loops_reset; \
r4_x4cl_wpn_four_complex_djbfft_preload; \
for(unsigned int loopA = 5; loopA; loopA--) { \
	for(unsigned int loopB = 4*2*clm; loopB; ) { \
		xtouch(xptr(rcx+rbp+(clm-1)*128));	/* Preload TLB */ \
		do{ \
			r4_x4cl_wpn_four_complex_djbfft(rsi, 64, 4*clmblkdst, clmblkdst8, rdi); \
			rdi += XMM_SCND2;		/* Next sine/cosine pointer */ \
		}while(--loopB & (2*clm-1)); \
		newxloop1(clm, b7b, clmblkdst);	/* 2*clm iterations */ \
	} \
	rsi += -4*clmblkdst+2*clmblkdst8; /* Next source pointer */ \
	rdi += -4*2*clm*XMM_SCND2;	/* Restore sin/cos ptr */ \
} \
rsi += -6*2*clmblkdst8;	/* Restore source pointer */ \
end_timer(25); \
\
/* Do FFT levels 8,9 */ \
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
/* Do 1 four_complex_fft4 operations */ \
/*	distance between fft data elements is 1K */ \
/*	do 1 macros each processing 32 data values */ \
/* Do 46 four_complex_fft4 operations */ \
/*	distance between fft data elements is 1K */ \
/*	do 23 macros each processing 32 data values */ \
\
start_timer(26); \
rdi += -4*2*clm*XMM_SCND2R-4*2*clm*XMM_SCD4; /* Complex data sin/cos ptr */ \
rdx = (rdi-2*clm*XMM_SCD2); /* Real data sin/cos ptr */ \
rcx = (uintptr_t)g->data_addr;/* Destination pointer */ \
loops_init_prefetch(24*2*clm, 64, 1, rbp); /* 24 iters of 2*clm */ \
pfing(rbp = (uintptr_t)g->premult_prefetch);/* Load premult/sin/cos prefetch pointer */ \
r4_sg4cl_eight_reals_four_complex_fft4_preload; \
for(unsigned int loopA = 2*clm; loopA; loopA--) { \
	r4_sg4cl_eight_reals_four_complex_fft4(rsi, 64, clmblkdst, 2*clmblkdst, rcx, 64, blkdst, 2*blkdst, rdi, rdx); \
	rdi += XMM_SCD4;		/* Next sine/cosine pointer */ \
	rdx += XMM_SCD2;		/* Next sine/cosine pointer */ \
} \
rsi += -2*clm*64+4*clmblkdst; /* Next source pointer */ \
rcx += -2*clm*64+4*blkdst;	/* Next dest pointer */ \
rdi += -2*clm*XMM_SCD4;	/* Restore premult/sine/cosine pointer */ \
\
loops_reset; \
for(unsigned int loopA = 4*3*2-1; loopA; ) { \
	for(unsigned int loopB = 3; loopB; loopB--) { \
		do{ \
			for(unsigned int loopD = 2*clm; loopD; loopD--) { \
				r4_sg4cl_four_complex_fft4(rsi, 64, clmblkdst, 2*clmblkdst, rcx, 64, blkdst, 2*blkdst, rdi); \
				rdi += XMM_SCD4;		/* Next premult/sin/cos ptr */ \
			} \
			rsi += -2*clm*64+4*clmblkdst; /* Next source pointer */ \
			rcx += -2*clm*64+4*blkdst;	/* Next dest pointer */ \
			rdi += -2*clm*XMM_SCD4;	/* Restore premult/sine/cosine pointer */ \
		}while(--loopA & 1); \
		rsi += -2*4*clmblkdst+clmblkdst8; /* Next source pointer */ \
	} \
	rdi += 2*clm*XMM_SCD4;	/* Next premult/sine/cosine pointer */ \
} \
end_timer(26); \
\
/* Work on next 4*clm sets of 384 values */ \
\
endlp:	xpass1_get_next_block(pass2, c0b, b0b, _wpn); /* Get next block and jump */ \
/* done: */ xfft_footer_rn; \
\
\
\
/* Perform a pass 1 FFT of length 448 doubles. */ \
/* Memory examples are for a 448K FFT with 10 levels done in pass 2. */ \
/* Sin/cos data needed: 4*2*clm*XMM_SCND2 + 4*2*clm*XMM_SCND2R + 2*clm*XMM_SCD2 + 14*2*clm*XMM_SCD4 = 43*2*clm 64-byte cache lines */ \
\

#define r4dwpn_pass1sc448(procname, xpass2name, pass2_levels, clm) /* (fftname, pass2_levels, clmarg) */ \
\
/* set_FFT_constants(fftname, 0, 0, pass2_levels, clmarg); */ \
xfft_header_rn(procname, xpass2name, pass2_levels, clm); \
\
/* Do pass 2 (the last 10 FFT levels if this is a 448K FFT) */ \
\
/* pass2: */start_timer(0); \
xpass2name(g); \
end_timer(0); \
if(g->ffttype == 1)		/* We're done if FFTing only */ \
	return; /* goto done; */ \
\
/* Do the inverse FFT */ \
\
pass1_inverse_fft_setup(c0b, procname); \
\
/* Do inverse FFT levels 8,9 */ \
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
/* Do 55 four_complex_unfft4 operations */ \
/*	distance between fft data elements is 1K */ \
/*	do 55 macros each processing 16 data values */ \
\
c0b:	start_timer(16); \
rsi = (uintptr_t)g->data_addr;/* Get source address */ \
rcx = (uintptr_t)g->scratch_area;	/* Get address of scratch area */ \
rdx = (uintptr_t)g->premult_addr;/* Real data sin/cos pointer */ \
rdi = (rdx+2*clm*XMM_SCD2); /* Complex data sin/cos ptr */ \
pfing(rbp = (rdx+(43-28)*2*clm*64)); /* Finish prefetching the sin/cos data */ \
loops_init_prefetch(56*2*clm, 64, 2, rbp); /* 56*2*clm iterations */ \
/*bug	xtouchpages rbp, 56*2*clm * 64 / 4 ;; Prefetch TLBs */ \
for(unsigned int loopA = 2*clm; loopA; loopA--) { \
	r4_sg2cl_eight_reals_unfft4(rsi, 64, blkdst, rcx, 64, clmblkdst, rdi, rdx); \
	rdi += XMM_SCD4;		/* Next premult/sin/cos data */ \
	rdx += XMM_SCD2;		/* Next sine/cosine pointer */ \
} \
rsi += -2*clm*64+2*blkdst;	/* Next source pointer */ \
rcx += -2*clm*64+2*clmblkdst; /* Next dest pointer */ \
rdi += -2*clm*XMM_SCD4;	/* Restore premult/sine/cosine pointer */ \
\
loops_reset; \
for(unsigned int loopA = 14*4-1; loopA; ) { \
	do{ \
		for(unsigned int loopC = 2*clm; loopC; loopC--) { \
			r4_sg2cl_four_complex_unfft4(rsi, 64, blkdst, rcx, 64, clmblkdst, rdi); \
			rdi += XMM_SCD4;		/* Next premult/sin/cos data */ \
		} \
		rsi += -2*clm*64+2*blkdst;	/* Next source pointer */ \
		rcx += -2*clm*64+2*clmblkdst; /* Next dest pointer */ \
		rdi += -2*clm*XMM_SCD4;	/* Restore premult/sin/cos data pointer */ \
	}while(--loopA & 3); \
	rcx += -4*2*clmblkdst+clmblkdst8; \
	rdi += 2*clm*XMM_SCD4;	/* Next premult/sin/cos pointer */ \
} \
rsi = (uintptr_t)g->scratch_area;	/* Restore source pointer */ \
rdx = rdi;		/* Next real data sin/cos pointer */ \
end_timer(16); \
\
/* Do inverse FFT levels 6,7 */ \
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
/* Do 52 four_complex_unfft operations */ \
/*	distance between fft data elements is 4K */ \
/*	do 26 macros each processing 32 data values */ \
\
start_timer(17); \
rdi = (rdx+4*2*clm*XMM_SCND2R); /* Complex data sin/cos ptr */ \
loops_init(28*2*clm);		/* 28*2*clm iterations */ \
\
r4_x4cl_wpn_eight_reals_unfft_preload; \
for(unsigned int loopA = 2; loopA; loopA--) { \
	for(unsigned int loopB = 2*clm; loopB; loopB--) { \
		r4_x4cl_wpn_eight_reals_unfft(rsi, 64, 2*clmblkdst, 4*clmblkdst, rdi, 2*2*clm*XMM_SCND2, rdx, 2*2*clm*XMM_SCND2R); \
		rdi += XMM_SCND2;		/* Next sine/cosine pointer */ \
		rdx += XMM_SCND2R;		/* Next sine/cosine pointer */ \
	} \
	rsi += -2*clm*64+clmblkdst; /* Next source pointer */ \
} \
rsi += -2*clmblkdst+clmblkdst8; /* Next source pointer */ \
rdi += -2*2*clm*XMM_SCND2;	/* Restore sine/cosine pointer */ \
\
loops_reset; \
r4_x4cl_wpn_four_complex_djbunfft_preload; \
for(unsigned int loopA = 13; loopA; loopA--) { \
	for(unsigned int loopB = 2; loopB; loopB--) { \
		for(unsigned int loopC = 2*clm; loopC; loopC--) { \
			r4_x4cl_wpn_four_complex_djbunfft(rsi, 64, 2*clmblkdst, 4*clmblkdst, rdi, 2*2*clm*XMM_SCND2); \
			rdi += XMM_SCND2;		/* Next clm's sin/cos data */ \
		} \
		rsi += -2*clm*64+clmblkdst; /* Next source pointer */ \
	} \
	rsi += -2*clmblkdst+clmblkdst8; /* Next source pointer */ \
	rdi += -2*2*clm*XMM_SCND2;	/* Restore sine/cosine pointer */ \
} \
rsi += -14*clmblkdst8;	/* Restore source pointer */ \
end_timer(17); \
\
/* Do inverse FFT levels 1,2,3,4,5 (actually 4.807 levels) */ \
/* On input the 64-byte cache lines hold these data values: */ \
/*	0K	8K	16K	24K */ \
/*	1K	... */ \
/*	... */ \
/*	7K */ \
/*	32K	... */ \
/*	... */ \
/* On output the 64-byte cache lines hold these data values: */ \
/*	0K	112K	224K	336K */ \
/*	1K	... */ \
/*	... */ \
/*	111K */ \
\
/* Do 64 28_reals_last_unfft operations */ \
/*	distance between fft data elements is 16K */ \
/*	do 8 macros each processing 112 data values */ \
\
start_timer(20); \
rdi = (uintptr_t)g->sincos2;/* Load fixed sin/cos postmultipliers */ \
set_data_prefetch_ptrs;		/* Init prefetching pointers */ \
loops_init(8*2*clm);		/* 8 iterations of 2*clm */ \
r7_x14cl_28_reals_last_unfft_preload; \
for(unsigned int loopA = 8*2*clm; loopA; ) { \
	xtouch(xptr(rcx+rbp+(clm-1)*128));	/* Preload TLB */ \
	do{ \
		r7_x14cl_28_reals_last_unfft(rsi, 64, clmblkdst8, rdi, 8*XMM_SCD13); \
	}while(--loopA & (2*clm-1)); \
	newxloop1quadpre(clm, c1b, c1a, clmblkdst, rdi, XMM_SCD13); /* 2*clm iterations */ \
} \
end_timer(20); \
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
/* pass1: */pass1_forward_fft_setup(b0b); \
\
/* Copy data from the FFT data area back to the scratch area */ \
\
b0b:	copy_fft_data_to_scratch(112, clm); \
\
/* Do FFT levels 1,2,3,4,5 (actually 4.807 levels) */ \
/* */ \
/* On input the 64-byte cache lines hold these data values: */ \
/*	0K	112K	224K	336K */ \
/*	1K	... */ \
/*	... */ \
/*	111K */ \
/* On output the 64-byte cache lines hold these data values: */ \
/*	0K	16K	32K	48K */ \
/*	1K	... */ \
/*	... */ \
/*	15K	... */ \
/*	64K	... */ \
/*	... */ \
\
/* Do 64 28_reals_first_fft operations */ \
/*	distance between fft data elements is 16K */ \
/*	do 16 iterations each processing 56 data values */ \
\
yes_fft:start_timer(22); \
rsi = (uintptr_t)g->scratch_area;	/* Get address of scratch area */ \
rdi = (uintptr_t)g->sincos2;/* Load fixed sin/cos postmultipliers */ \
set_data_prefetch_ptrs;		/* Init prefetching pointers */ \
pfing(rcx += (32-4)*blkdst);	/* We've already prefetched 32 blocks, but */ \
;					/* code below will fetch 4 blocks too many */ \
loops_init(16*2*clm);		/* 16 iters of 2*clm */ \
r7_x7cl_28_reals_first_fft_scratch_preload; \
for(unsigned int loopA = 2; loopA; loopA--) { \
	for(unsigned int loopB = 8*2*clm; loopB; ) { \
		xtouch(xptr(rcx+rbp+(clm-1)*128));	/* Preload TLB */ \
		do{ \
			r7_x7cl_28_reals_first_fft_scratch(rsi, 64, 2*clmblkdst8, rdi); \
		}while(--loopB & (2*clm-1)); \
		newxloop1doublepre(clm, b1b, b1a, clmblkdst, rdi, XMM_SCD13); /* 2*clm iters prefetching 4*clm cache lines */ \
	} \
	rsi += -8*clmblkdst+clmblkdst8; \
} \
rsi += -2*clmblkdst8;	/* Restore source pointer */ \
end_timer(22); \
\
/* Do FFT levels 6,7 */ \
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
/* Do 48 four_complex_fft operations */ \
/*	distance between fft data elements is 4K */ \
/*	do 24 macros each processing 32 data values */ \
\
start_timer(25); \
rdx = (uintptr_t)g->premult_addr;/* Address of premultiplier/sin/cos table */ \
rdx += 2*clm*XMM_SCD2 + 14*2*clm*XMM_SCD4; /* Real data sin/cos ptr */ \
rdi = (rdx+4*2*clm*XMM_SCND2R); /* Complex data sin/cos ptr */ \
loops_init(28*2*clm);		/* 7 iterations of 4 iters of 2*clm */ \
r4_x4cl_wpn_eight_reals_four_complex_djbfft_preload; \
for(unsigned int loopA = 4*2*clm; loopA; ) { \
	xtouch(xptr(rcx+rbp+(clm-1)*128));	/* Preload TLB */ \
	do{ \
		r4_x4cl_wpn_eight_reals_four_complex_djbfft(rsi, 64, 4*clmblkdst, clmblkdst8, rdi, rdx); \
		rdi += XMM_SCND2;		/* Next sine/cosine pointer */ \
		rdx += XMM_SCND2R;		/* Next sine/cosine pointer */ \
	}while(--loopA & (2*clm-1)); \
	newxloop1(clm, b6b, clmblkdst); /* 2*clm iterations */ \
} \
rsi += -4*clmblkdst+2*clmblkdst8; \
rdi += -4*2*clm*XMM_SCND2;	/* Restore sin/cos ptr */ \
\
loops_reset; \
r4_x4cl_wpn_four_complex_djbfft_preload; \
for(unsigned int loopA = 6; loopA; loopA--) { \
	for(unsigned int loopB = 4*2*clm; loopB; ) { \
		xtouch(xptr(rcx+rbp+(clm-1)*128));	/* Preload TLB */ \
		do{ \
			r4_x4cl_wpn_four_complex_djbfft(rsi, 64, 4*clmblkdst, clmblkdst8, rdi); \
			rdi += XMM_SCND2;		/* Next sine/cosine pointer */ \
		}while(--loopB & (2*clm-1)); \
		newxloop1doublepre(clm, b7b, b7a, clmblkdst); /* 2*clm iterations */ \
	} \
	rsi += -4*clmblkdst+2*clmblkdst8; /* Next source pointer */ \
	rdi += -4*2*clm*XMM_SCND2;	/* Restore sin/cos ptr */ \
} \
rsi += -7*2*clmblkdst8;	/* Restore source pointer */ \
end_timer(25); \
\
/* Do FFT levels 8,9 */ \
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
/* Do 1 four_complex_fft4 operations */ \
/*	distance between fft data elements is 1K */ \
/*	do 1 macros each processing 32 data values */ \
/* Do 54 four_complex_fft4 operations */ \
/*	distance between fft data elements is 1K */ \
/*	do 27 macros each processing 32 data values */ \
\
start_timer(26); \
rdi += -4*2*clm*XMM_SCND2R-14*2*clm*XMM_SCD4; /* Complex data sin/cos ptr */ \
rdx = (rdi-2*clm*XMM_SCD2); /* Real data sin/cos ptr */ \
rcx = (uintptr_t)g->data_addr;/* Destination pointer */ \
loops_init_prefetch(28*2*clm, 64, 1, rbp); /* 28 iters of 2*clm */ \
pfing(rbp = (uintptr_t)g->premult_prefetch);/* Load premult/sin/cos prefetch pointer */ \
r4_sg4cl_eight_reals_four_complex_fft4_preload; \
for(unsigned int loopA = 2*clm; loopA; loopA--) { \
	r4_sg4cl_eight_reals_four_complex_fft4(rsi, 64, clmblkdst, 2*clmblkdst, rcx, 64, blkdst, 2*blkdst, rdi, rdx); \
	rdi += XMM_SCD4;		/* Next sine/cosine pointer */ \
	rdx += XMM_SCD2;		/* Next sine/cosine pointer */ \
} \
rsi += -2*clm*64+4*clmblkdst; /* Next source pointer */ \
rcx += -2*clm*64+4*blkdst;	/* Next dest pointer */ \
rdi += -2*clm*XMM_SCD4;	/* Restore premult/sine/cosine pointer */ \
\
loops_reset; \
for(unsigned int loopA = 14*2-1; loopA; ) { \
	do{ \
		for(unsigned int loopC = 2*clm; loopC; loopC--) { \
			r4_sg4cl_four_complex_fft4(rsi, 64, clmblkdst, 2*clmblkdst, rcx, 64, blkdst, 2*blkdst, rdi); \
			rdi += XMM_SCD4;		/* Next premult/sin/cos ptr */ \
		} \
		rsi += -2*clm*64+4*clmblkdst; /* Next source pointer */ \
		rcx += -2*clm*64+4*blkdst;	/* Next dest pointer */ \
		rdi += -2*clm*XMM_SCD4;	/* Restore premult/sine/cosine pointer */ \
	}while(--loopA & 1); \
	rsi += -2*4*clmblkdst+clmblkdst8; /* Next source pointer */ \
	rdi += 2*clm*XMM_SCD4;	/* Next premult/sine/cosine pointer */ \
} \
end_timer(26); \
\
/* Work on next 4*clm sets of 448 values */ \
\
endlp:	xpass1_get_next_block(pass2, c0b, b0b, _wpn); /* Get next block and jump */ \
/* done: */ xfft_footer_rn; \
\
\
\
/* Perform a pass 1 FFT of length 512 doubles. */ \
/* Memory examples are for a 512K FFT with 10 levels done in pass 2. */ \
/* Sin/cos data needed: 4*2*clm*XMM_SCND2 + 4*2*clm*XMM_SCND2R + 2*clm*XMM_SCD2 + 16*2*clm*XMM_SCD4 = 47*2*clm 64-byte cache lines */ \
\

#define r4dwpn_pass1sc512(procname, xpass2name, pass2_levels, clm) /* (fftname, pass2_levels, clmarg) */ \
\
/* set_FFT_constants(fftname, 0, 0, pass2_levels, clmarg); */ \
xfft_header_rn(procname, xpass2name, pass2_levels, clm); \
\
/* Do pass 2 (the last 10 FFT levels if this is a 512K FFT) */ \
\
/* pass2: */start_timer(0); \
xpass2name(g); \
end_timer(0); \
if(g->ffttype == 1)		/* We're done if FFTing only */ \
	return; /* goto done; */ \
\
/* Do the inverse FFT */ \
\
pass1_inverse_fft_setup(c0b, procname); \
\
/* Do inverse FFT levels 8,9 */ \
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
/* Do 63 four_complex_unfft4 operations */ \
/*	distance between fft data elements is 1K */ \
/*	do 63 macros each processing 16 data values */ \
\
c0b:	start_timer(16); \
rsi = (uintptr_t)g->data_addr;/* Get source address */ \
rcx = (uintptr_t)g->scratch_area;	/* Get address of scratch area */ \
rdx = (uintptr_t)g->premult_addr;/* Real data sin/cos pointer */ \
rdi = (rdx+2*clm*XMM_SCD2); /* Complex data sin/cos ptr */ \
pfing(rbp = (rdx+(47-16)*2*clm*64)); /* Finish prefetch the sin/cos/premult data */ \
loops_init_prefetch_touch(64*2*clm, 64, 4, rbp, 4096, c9a); /* 64*2*clm iterations */ \
xtouch(xptr(rbp));			/* Prefetch TLB */ \
for(unsigned int loopA = 2*clm; loopA; loopA--) { \
	r4_sg2cl_eight_reals_unfft4(rsi, 64, blkdst, rcx, 64, clmblkdst, rdi, rdx); \
	rdi += XMM_SCD4;		/* Next premult/sin/cos data */ \
	rdx += XMM_SCD2;		/* Next sine/cosine pointer */ \
} \
rsi += -2*clm*64+2*blkdst;	/* Next source pointer */ \
rcx += -2*clm*64+2*clmblkdst; /* Next dest pointer */ \
rdi += -2*clm*XMM_SCD4;	/* Restore premult/sine/cosine pointer */ \
\
loops_reset; \
/* c9a: */ xtouch(xptr(rbp+4096));		/* Prefetch TLB */ \
for(unsigned int loopA = 16*4-1; loopA; ) { \
	do{ \
		for(unsigned int loopC = 2*clm; loopC; loopC--) { \
			r4_sg2cl_four_complex_unfft4(rsi, 64, blkdst, rcx, 64, clmblkdst, rdi); \
			rdi += XMM_SCD4;		/* Next premult/sin/cos data */ \
		} \
		rsi += -2*clm*64+2*blkdst;	/* Next source pointer */ \
		rcx += -2*clm*64+2*clmblkdst; /* Next dest pointer */ \
		rdi += -2*clm*XMM_SCD4;	/* Restore premult/sin/cos data pointer */ \
	}while(--loopA & 3); \
	rcx += -4*2*clmblkdst+clmblkdst8; \
	rdi += 2*clm*XMM_SCD4;	/* Next premult/sin/cos pointer */ \
} \
rsi = (uintptr_t)g->scratch_area;	/* Restore source pointer */ \
rdx = rdi;		/* Next real data sin/cos pointer */ \
end_timer(16); \
\
/* Do inverse FFT levels 6,7 */ \
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
/* Do 60 four_complex_unfft operations */ \
/*	distance between fft data elements is 4K */ \
/*	do 30 macros each processing 32 data values */ \
\
start_timer(17); \
loops_init(32*2*clm);		/* 32*2*clm iterations */ \
rdi = (rdx+4*2*clm*XMM_SCND2R); /* Next complex data sin/cos ptr */ \
\
r4_x4cl_wpn_eight_reals_unfft_preload; \
for(unsigned int loopA = 2; loopA; loopA--) { \
	for(unsigned int loopB = 2*clm; loopB; loopB--) { \
		r4_x4cl_wpn_eight_reals_unfft(rsi, 64, 2*clmblkdst, 4*clmblkdst, rdi, 2*2*clm*XMM_SCND2, rdx, 2*2*clm*XMM_SCND2R); \
		rdi += XMM_SCND2;		/* Next sine/cosine pointer */ \
		rdx += XMM_SCND2R;		/* Next sine/cosine pointer */ \
	} \
	rsi += -2*clm*64+clmblkdst; /* Next source pointer */ \
} \
rsi += -2*clmblkdst+clmblkdst8; /* Next source pointer */ \
rdi += -2*2*clm*XMM_SCND2;	/* Restore sine/cosine pointer */ \
\
loops_reset; \
r4_x4cl_wpn_four_complex_djbunfft_preload; \
for(unsigned int loopA = 15; loopA; loopA--) { \
	for(unsigned int loopB = 2; loopB; loopB--) { \
		for(unsigned int loopC = 2*clm; loopC; loopC--) { \
			r4_x4cl_wpn_four_complex_djbunfft(rsi, 64, 2*clmblkdst, 4*clmblkdst, rdi, 2*2*clm*XMM_SCND2); \
			rdi += XMM_SCND2;		/* Next clm's sin/cos data */ \
		} \
		rsi += -2*clm*64+clmblkdst; /* Next source pointer */ \
	} \
	rsi += -2*clmblkdst+clmblkdst8; /* Next source pointer */ \
	rdi += -2*2*clm*XMM_SCND2;	/* Restore sine/cosine pointer */ \
} \
rsi += -16*clmblkdst8;	/* Restore source pointer */ \
end_timer(17); \
\
/* Do inverse FFT levels 4,5 */ \
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
/*	128K	... */ \
/*	... */ \
\
/* Do 16 eight_reals_unfft operations */ \
/*	distance between fft data elements is 16K */ \
/*	do 8 macros each processing 32 data values */ \
/* Do 48 four_complex_unfft operations */ \
/*	distance between fft data elements is 16K */ \
/*	do 24 macros each processing 32 data values */ \
\
start_timer(18); \
rdx = (uintptr_t)g->sincos2;/* Load fixed sin/cos data ptr */ \
rdx += 64*XMM_SCD3;	/* Real data sin/cos ptr */ \
rdi = (rdx+16*XMM_SCD2);	/* Complex data sin/cos ptr */ \
loops_init(32*2*clm);		/* 32 iterations of 2*clm */ \
r4_x4cl_eight_reals_unfft_preload; \
for(unsigned int loopA = 8; loopA; loopA--) { \
	for(unsigned int loopB = 2*clm; loopB; loopB--) { \
		r4_x4cl_eight_reals_unfft(rsi, 64, clmblkdst8, 2*clmblkdst8, rdi, 8*XMM_SCD2, rdx, 8*XMM_SCD2); \
	} \
	rsi += -2*clm*64+clmblkdst; /* Next source pointer */ \
	rdi += XMM_SCD2;		/* Next sine/cosine pointer */ \
	rdx += XMM_SCD2;		/* Next sine/cosine pointer */ \
} \
rsi += -8*clmblkdst+4*clmblkdst8; /* Next source pointer */ \
rdi += -8*XMM_SCD2;	/* Restore sine/cosine pointer */ \
\
loops_reset; \
r4_x4cl_four_complex_djbunfft_preload; \
for(unsigned int loopA = 3; loopA; loopA--) { \
	for(unsigned int loopB = 8; loopB; loopB--) { \
		for(unsigned int loopC = 2*clm; loopC; loopC--) { \
			r4_x4cl_four_complex_djbunfft(rsi, 64, clmblkdst8, 2*clmblkdst8, rdi, 8*XMM_SCD2); \
		} \
		rsi += -2*clm*64+clmblkdst; /* Next source pointer */ \
		rdi += XMM_SCD2;		/* Next clm's sin/cos data */ \
	} \
	rsi += -8*clmblkdst+4*clmblkdst8; /* Next source pointer */ \
	rdi += -8*XMM_SCD2;	/* Restore sine/cosine pointer */ \
} \
rsi += -4*4*clmblkdst8;	/* Restore source pointer */ \
end_timer(18); \
\
/* Do inverse FFT levels 1,2,3 */ \
/* On input the 64-byte cache lines hold these data values: */ \
/*	0K	32K	64K	96K */ \
/*	1K	... */ \
/*	... */ \
/*	31K */ \
/*	128K	... */ \
/*	... */ \
/* On output the 64-byte cache lines hold these data values: */ \
/*	0K	128K	256K	384K */ \
/*	1K	... */ \
/*	... */ \
/*	127K */ \
\
/* Do 64 eight_reals_last_unfft operations */ \
/*	distance between fft data elements is 64K */ \
/*	do 32 macros each processing 32 data values */ \
\
start_timer(20); \
rdi = (uintptr_t)g->sincos2;/* Load fixed sin/cos data ptr */ \
set_data_prefetch_ptrs;		/* Init prefetching pointers */ \
loops_init(32*2*clm);		/* 32 iterations of 2*clm */ \
r4_x4cl_eight_reals_last_unfft_preload; \
for(unsigned int loopA = 4; loopA; loopA--) { \
	for(unsigned int loopB = 8*2*clm; loopB; ) { \
		xtouch(xptr(rcx+rbp+(clm-1)*128));	/* Preload TLB */ \
		do{ \
			r4_x4cl_eight_reals_last_unfft(rsi, 64, 4*clmblkdst8, 8*clmblkdst8, rdi, 32*XMM_SCD3); \
		}while(--loopB & (2*clm-1)); \
		newxloop1(clm, c1b, clmblkdst, rdi, XMM_SCD3); /* 2*clm iterations */ \
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
copy_scratch_data_to_fft(128, clm); \
goto	endlp;			/* Skip FFT code if flag not set */ \
\
/* Fire up auxiliary threads to do pass 1 forward FFT */ \
\
/* pass1: */pass1_forward_fft_setup(b0b); \
\
/* Copy data from the FFT data area back to the scratch area */ \
\
b0b:	copy_fft_data_to_scratch(128, clm); \
\
/* Do FFT levels 1,2,3 */ \
/* */ \
/* On input the 64-byte cache lines hold these data values: */ \
/*	0K	128K	256K	384K */ \
/*	1K	... */ \
/*	... */ \
/*	127K */ \
/* On output the 64-byte cache lines hold these data values: */ \
/*	0K	64K	128K	192K */ \
/*	1K	... */ \
/*	... */ \
/*	63K	... */ \
/*	256K	... */ \
/*	... */ \
\
/* Do 64 eight_reals_first_fft operations */ \
/*	distance between fft data elements is 64K */ \
/*	do 64 iterations each processing 16 data values */ \
\
yes_fft:start_timer(22); \
rsi = (uintptr_t)g->scratch_area;	/* Get address of scratch area */ \
rdi = (uintptr_t)g->sincos2;/* Load fixed sin/cos postmultipliers */ \
set_data_prefetch_ptrs;		/* Init prefetching pointers */ \
pfing(rcx += 32*blkdst);		/* We've already prefetched 32 blocks */ \
loops_init(64*2*clm);		/* 8 iters of 4 iters of 2*2*clm */ \
r4_x2cl_eight_reals_first_fft_scratch_preload; \
for(unsigned int loopA = 8; loopA; loopA--) { \
	for(unsigned int loopB = 4*2*2*clm; loopB; ) { \
		xtouch(xptr(rcx+rbp+(clm-1)*128));	/* Preload TLB */ \
		do{ \
			do{ \
				r4_x2cl_eight_reals_first_fft_scratch(rsi, 64, 8*clmblkdst8, rdi); \
			}while(--loopB & (2*clm-1)); \
			newxloop2(clm, b1b, clmblkdst, rdi, XMM_SCD3); /* 2*2*clm iters prefetching 2*clm cache lines */ \
		}while(loopB & (2*2*clm-1)); \
	} \
	rsi += -4*2*clmblkdst+clmblkdst8; \
} \
rsi += -8*clmblkdst8;	/* Restore source pointer */ \
rdx = rdi;		/* Next real data sin/cos ptr */ \
end_timer(22); \
\
/* Do FFT levels 4,5 */ \
/* */ \
/* On input the 64-byte cache lines hold these data values: */ \
/*	0K	64K	128K	192K */ \
/*	1K	... */ \
/*	... */ \
/*	63K	... */ \
/*	256K	... */ \
/*	... */ \
/* On output the 64-byte cache lines hold these data values: */ \
/*	0K	16K	32K	48K */ \
/*	1K	... */ \
/*	... */ \
/*	15K	... */ \
/*	64K	... */ \
/*	... */ \
\
/* Do 16 eight_reals_fft operations */ \
/* Do 16 four_complex_fft operations */ \
/*	distance between fft data elements is 16K */ \
/*	do 16 macros each processing 32 data values */ \
/* Do 32 four_complex_fft operations */ \
/*	distance between fft data elements is 16K */ \
/*	do 16 macros each processing 32 data values */ \
\
start_timer(24); \
rdi = (rdx+16*XMM_SCD2);	/* Complex data sin/cos ptr */ \
loops_init(32*2*clm);		/* 4 iters of 8 of 2*clm */ \
r4_x4cl_eight_reals_four_complex_djbfft_preload; \
for(unsigned int loopA = 2; loopA; loopA--) { \
	for(unsigned int loopB = 8*2*clm; loopB; ) { \
		xtouch(xptr(rcx+rbp+(clm-1)*128));	/* Preload TLB */ \
		do{ \
			r4_x4cl_eight_reals_four_complex_djbfft(rsi, 64, 2*clmblkdst8, 4*clmblkdst8, rdi, rdx); \
		}while(--loopB & (2*clm-1)); \
		newxloop1(clm, b4b, clmblkdst, rdi, XMM_SCD2, rdx, XMM_SCD2); /* 2*clm iterations */ \
	} \
	rsi += -8*clmblkdst+clmblkdst8; \
} \
rsi += -2*clmblkdst8+8*clmblkdst8; /* Next source pointer */ \
rdi += -16*XMM_SCD2;	/* Restore sin/cos ptr */ \
\
loops_reset; \
r4_x4cl_four_complex_djbfft_preload; \
for(unsigned int loopA = 2; loopA; loopA--) { \
	for(unsigned int loopB = 8*2*clm; loopB; ) { \
		xtouch(xptr(rcx+rbp+(clm-1)*128));	/* Preload TLB */ \
		do{ \
			r4_x4cl_four_complex_djbfft(rsi, 64, 2*clmblkdst8, 4*clmblkdst8, rdi); \
		}while(--loopB & (2*clm-1)); \
		newxloop1(clm, b5b, clmblkdst, rdi, XMM_SCD2); /* 2*clm iterations */ \
	} \
	rsi += -8*clmblkdst+clmblkdst8; \
} \
rsi += -2*clmblkdst8-8*clmblkdst8; /* Restore source pointer */ \
end_timer(24); \
\
/* Do FFT levels 6,7 */ \
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
/* Do 56 four_complex_fft operations */ \
/*	distance between fft data elements is 4K */ \
/*	do 28 macros each processing 32 data values */ \
\
start_timer(25); \
rdx = (uintptr_t)g->premult_addr;/* Address of premultiplier/sin/cos table */ \
rdx += 2*clm*XMM_SCD2 + 16*2*clm*XMM_SCD4; /* Real data sin/cos ptr */ \
rdi = (rdx+4*2*clm*XMM_SCND2R); /* Complex data sin/cos ptr */ \
loops_init(32*2*clm);		/* 8 iterations of 4 iters of 2*clm */ \
r4_x4cl_wpn_eight_reals_four_complex_djbfft_preload; \
for(unsigned int loopA = 4*2*clm; loopA; ) { \
	xtouch(xptr(rcx+rbp+(clm-1)*128));	/* Preload TLB */ \
	do{ \
		r4_x4cl_wpn_eight_reals_four_complex_djbfft(rsi, 64, 4*clmblkdst, clmblkdst8, rdi, rdx); \
		rdi += XMM_SCND2;		/* Next sine/cosine pointer */ \
		rdx += XMM_SCND2R;		/* Next sine/cosine pointer */ \
	}while(--loopA & (2*clm-1)); \
	newxloop1(clm, b6b, clmblkdst);	/* 2*clm iterations */ \
} \
rsi += -4*clmblkdst+2*clmblkdst8; \
rdi += -4*2*clm*XMM_SCND2;	/* Restore sin/cos ptr */ \
\
loops_reset; \
r4_x4cl_wpn_four_complex_djbfft_preload; \
for(unsigned int loopA = 7; loopA; loopA--) { \
	for(unsigned int loopB = 4*2*clm; loopB; ) { \
		xtouch(xptr(rcx+rbp+(clm-1)*128));	/* Preload TLB */ \
		do{ \
			r4_x4cl_wpn_four_complex_djbfft(rsi, 64, 4*clmblkdst, clmblkdst8, rdi); \
			rdi += XMM_SCND2;		/* Next sine/cosine pointer */ \
		}while(--loopB & (2*clm-1)); \
		newxloop1(clm, b7b, clmblkdst);	/* 2*clm iterations */ \
	} \
	rsi += -4*clmblkdst+2*clmblkdst8; /* Next source pointer */ \
	rdi += -4*2*clm*XMM_SCND2;	/* Restore sin/cos ptr */ \
} \
rsi += -8*2*clmblkdst8;	/* Restore source pointer */ \
end_timer(25); \
\
/* Do FFT levels 8,9 */ \
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
/* Do 1 four_complex_fft4 operations */ \
/*	distance between fft data elements is 1K */ \
/*	do 1 macros each processing 32 data values */ \
/* Do 62 four_complex_fft4 operations */ \
/*	distance between fft data elements is 1K */ \
/*	do 31 macros each processing 32 data values */ \
\
start_timer(26); \
rdi += -4*2*clm*XMM_SCND2R-16*2*clm*XMM_SCD4; /* Complex data sin/cos ptr */ \
rdx = (rdi-2*clm*XMM_SCD2); /* Real data sin/cos ptr */ \
rcx = (uintptr_t)g->data_addr;/* Destination pointer */ \
loops_init_prefetch(32*2*clm, 64, 1, rbp); /* 32 iters of 2*clm */ \
pfing(rbp = (uintptr_t)g->premult_prefetch);/* Load premult/sin/cos prefetch pointer */ \
r4_sg4cl_eight_reals_four_complex_fft4_preload; \
for(unsigned int loopA = 2*clm; loopA; loopA--) { \
	r4_sg4cl_eight_reals_four_complex_fft4(rsi, 64, clmblkdst, 2*clmblkdst, rcx, 64, blkdst, 2*blkdst, rdi, rdx); \
	rdi += XMM_SCD4;		/* Next sine/cosine pointer */ \
	rdx += XMM_SCD2;		/* Next sine/cosine pointer */ \
} \
rsi += -2*clm*64+4*clmblkdst; /* Next source pointer */ \
rcx += -2*clm*64+4*blkdst;	/* Next dest pointer */ \
rdi += -2*clm*XMM_SCD4;	/* Restore premult/sine/cosine pointer */ \
\
loops_reset; \
for(unsigned int loopA = 16*2-1; loopA; ) { \
	do{ \
		for(unsigned int loopC = 2*clm; loopC; loopC--) { \
			r4_sg4cl_four_complex_fft4(rsi, 64, clmblkdst, 2*clmblkdst, rcx, 64, blkdst, 2*blkdst, rdi); \
			rdi += XMM_SCD4;		/* Next premult/sin/cos ptr */ \
		} \
		rsi += -2*clm*64+4*clmblkdst; /* Next source pointer */ \
		rcx += -2*clm*64+4*blkdst;	/* Next dest pointer */ \
		rdi += -2*clm*XMM_SCD4;	/* Restore premult/sine/cosine pointer */ \
	}while(--loopA & 1); \
	rsi += -2*4*clmblkdst+clmblkdst8; /* Next source pointer */ \
	rdi += 2*clm*XMM_SCD4;	/* Next premult/sine/cosine pointer */ \
} \
end_timer(26); \
\
/* Work on next 4*clm sets of 512 values */ \
\
endlp:	xpass1_get_next_block(pass2, c0b, b0b, _wpn); /* Get next block and jump */ \
/* done: */ xfft_footer_rn; \
\
\
\
/* Perform a pass 1 FFT of length 640 doubles. */ \
/* Memory examples are for a 640K FFT with 10 levels done in pass 2. */ \
/* Sin/cos data needed: 8*2*clm*XMM_SCND2 + 8*2*clm*XMM_SCND2R + 2*clm*XMM_SCD4 + 10*2*clm*XMM_SCD8 = 70*2*clm 64-byte cache lines */ \
\

#define r4dwpn_pass1sc640(procname, xpass2name, pass2_levels, clm) /* (fftname, pass2_levels, clmarg) */ \
\
/* set_FFT_constants(fftname, 0, 0, pass2_levels, clmarg); */ \
xfft_header_rn(procname, xpass2name, pass2_levels, clm); \
\
/* Do pass 2 (the last 10 FFT levels if this is a 640K FFT) */ \
\
/* pass2: */start_timer(0); \
xpass2name(g); \
end_timer(0); \
if(g->ffttype == 1)		/* We're done if FFTing only */ \
	return; /* goto done; */ \
\
/* Do the inverse FFT */ \
\
pass1_inverse_fft_setup(c0b, procname); \
\
/* Do inverse FFT levels 8,9,10 */ \
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
/* Do 39 eight_complex_unfft8 operations */ \
/*	distance between fft data elements is 1K */ \
/*	do 39 macros each processing 32 data values */ \
\
c0b:	start_timer(16); \
rsi = (uintptr_t)g->data_addr;/* Get source address */ \
rcx = (uintptr_t)g->scratch_area;	/* Get address of scratch area */ \
rdx = (uintptr_t)g->premult_addr;/* Real data sin/cos pointer */ \
rdi = (rdx+2*clm*XMM_SCD4+128); /* Complex data sin/cos ptr and use -128 */ \
;					/* addressing to generate 1 byte offsets in radix-8 macro */ \
pfing(rbp = (rdx+(70-40)*2*clm*64)); /* Finish prefetching the sin/cos data */ \
loops_init_prefetch(40*2*clm, 64, 1, rbp); /* 40*2*clm iterations */ \
/*bug	xtouchpages rbp, 40*2*clm * 64 / 1 ;; Prefetch TLBs */ \
for(unsigned int loopA = 2*clm; loopA; loopA--) { \
	r8_sg4cl_sixteen_reals_unfft8(rsi, 64, blkdst, 2*blkdst, rcx, 64, clmblkdst, 2*clmblkdst, rdi-128, rdx); \
	rdi += XMM_SCD8;		/* Next premult/sin/cos data */ \
	rdx += XMM_SCD4;		/* Next sine/cosine pointer */ \
} \
rsi += -2*clm*64+4*blkdst;	/* Next source pointer */ \
rcx += -2*clm*64+4*clmblkdst; /* Next dest pointer */ \
rdi += -2*clm*XMM_SCD8;	/* Restore premult/sine/cosine pointer */ \
\
loops_reset; \
for(unsigned int loopA = 10*2*2-1; loopA; ) { \
	do{ \
		do{ \
			for(unsigned int loopD = 2*clm; loopD; loopD--) { \
				r8_sg4cl_eight_complex_unfft8(rsi, 64, blkdst, 2*blkdst, rcx, 64, clmblkdst, 2*clmblkdst, rdi-128); \
				rdi += XMM_SCD8;		/* Next premult/sin/cos data */ \
			} \
			rsi += -2*clm*64+4*blkdst;	/* Next source pointer */ \
			rcx += -2*clm*64+4*clmblkdst; /* Next dest pointer */ \
			rdi += -2*clm*XMM_SCD8;	/* Restore premult/sin/cos data pointer */ \
		}while(--loopA & 1); \
		rcx += -2*4*clmblkdst+clmblkdst8; \
	}while(loopA & (2*2-1)); \
	rdi += 2*clm*XMM_SCD8;	/* Next premult/sin/cos pointer */ \
} \
rsi = (uintptr_t)g->scratch_area;	/* Restore source pointer */ \
rdx = (rdi-128);		/* Next real data sin/cos pointer */ \
end_timer(16); \
\
/* Do inverse FFT levels 6,7 */ \
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
/* Do 72 four_complex_unfft operations */ \
/*	distance between fft data elements is 8K */ \
/*	do 36 macros each processing 32 data values */ \
\
start_timer(17); \
rdi = (rdx+8*2*clm*XMM_SCND2R); /* Complex data sin/cos ptr */ \
loops_init(40*2*clm);		/* 40 iterations of 2*clm */ \
\
r4_x4cl_wpn_eight_reals_unfft_preload; \
for(unsigned int loopA = 4; loopA; loopA--) { \
	for(unsigned int loopB = 2*clm; loopB; loopB--) { \
		r4_x4cl_wpn_eight_reals_unfft(rsi, 64, 4*clmblkdst, clmblkdst8, rdi, 4*2*clm*XMM_SCND2, rdx, 4*2*clm*XMM_SCND2R); \
		rdi += XMM_SCND2;		/* Next sine/cosine pointer */ \
		rdx += XMM_SCND2R;		/* Next sine/cosine pointer */ \
	} \
	rsi += -2*clm*64+clmblkdst; /* Next source pointer */ \
} \
rsi += -4*clmblkdst+2*clmblkdst8; /* Next source pointer */ \
rdi += -4*2*clm*XMM_SCND2;	/* Restore sine/cosine pointer */ \
\
loops_reset; \
r4_x4cl_wpn_four_complex_djbunfft_preload; \
for(unsigned int loopA = 9; loopA; loopA--) { \
	for(unsigned int loopB = 4; loopB; loopB--) { \
		for(unsigned int loopC = 2*clm; loopC; loopC--) { \
			r4_x4cl_wpn_four_complex_djbunfft(rsi, 64, 4*clmblkdst, clmblkdst8, rdi, 4*2*clm*XMM_SCND2); \
			rdi += XMM_SCND2;		/* Next clm's sin/cos data */ \
		} \
		rsi += -2*clm*64+clmblkdst; /* Next source pointer */ \
	} \
	rsi += -4*clmblkdst+2*clmblkdst8; /* Next source pointer */ \
	rdi += -4*2*clm*XMM_SCND2;	/* Restore sine/cosine pointer */ \
} \
rsi += -10*2*clmblkdst8;	/* Restore source pointer */ \
end_timer(17); \
\
/* Do inverse FFT levels 1,2,3,4,5 (actually 4.32 levels) */ \
/* On input the 64-byte cache lines hold these data values: */ \
/*	0K	16K	32K	48K */ \
/*	1K	... */ \
/*	... */ \
/*	15K */ \
/*	64K	... */ \
/*	... */ \
/* On output the 64-byte cache lines hold these data values: */ \
/*	0K	160K	320K	480K */ \
/*	1K	... */ \
/*	... */ \
/*	159K */ \
\
/* Do 128 20_reals_last_unfft operations */ \
/*	distance between fft data elements is 32K */ \
/*	do 16 macros each processing 80 data values */ \
\
start_timer(20); \
rdi = (uintptr_t)g->sincos2;/* Load fixed sin/cos postmultipliers */ \
set_data_prefetch_ptrs;		/* Init prefetching pointers */ \
loops_init(16*2*clm);		/* 16 iterations of 2*clm */ \
r5_x10cl_20_reals_last_unfft_preload; \
for(unsigned int loopA = 2; loopA; loopA--) { \
	for(unsigned int loopB = 8*2*clm; loopB; ) { \
		xtouch(xptr(rcx+rbp+(clm-1)*128));	/* Preload TLB */ \
		do{ \
			r5_x10cl_20_reals_last_unfft(rsi, 64, 2*clmblkdst8, rdi, 16*XMM_SCD9); \
		}while(--loopB & (2*clm-1)); \
		newxloop1quadpre(clm, c1b, c1a, clmblkdst, rdi, XMM_SCD9); /* 2*clm iterations */ \
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
copy_scratch_data_to_fft(160, clm); \
goto	endlp;			/* Skip FFT code if flag not set */ \
\
/* Fire up auxiliary threads to do pass 1 forward FFT */ \
\
/* pass1: */pass1_forward_fft_setup(b0b); \
\
/* Copy data from the FFT data area back to the scratch area */ \
\
b0b:	copy_fft_data_to_scratch(160, clm); \
\
/* Do FFT levels 1,2,3,4,5 (actually 4.32 levels) */ \
/* */ \
/* On input the 64-byte cache lines hold these data values: */ \
/*	0K	160K	320K	480K */ \
/*	1K	... */ \
/*	... */ \
/*	159K */ \
/* On output the 64-byte cache lines hold these data values: */ \
/*	0K	32K	64K	96K */ \
/*	1K	... */ \
/*	... */ \
/*	31K	... */ \
/*	128K	... */ \
/*	... */ \
\
/* Do 128 20_reals_first_fft operations */ \
/*	distance between fft data elements is 32K */ \
/*	do 32 iterations each processing 40 data values */ \
\
yes_fft:start_timer(22); \
rsi = (uintptr_t)g->scratch_area;	/* Get address of scratch area */ \
rdi = (uintptr_t)g->sincos2;/* Load fixed sin/cos postmultipliers */ \
set_data_prefetch_ptrs;		/* Init prefetching pointers */ \
pfing(rcx += (64-8)*blkdst);	/* We've already prefetched 64 blocks, but */ \
;					/* code below will fetch 8 blocks too many */ \
loops_init(32*2*clm);		/* 32 iters of 2*clm */ \
r5_x5cl_20_reals_first_fft_scratch_preload; \
for(unsigned int loopA = 4; loopA; loopA--) { \
	for(unsigned int loopB = 8*2*clm; loopB; ) { \
		xtouch(xptr(rcx+rbp+(clm-1)*128));	/* Preload TLB */ \
		do{ \
			r5_x5cl_20_reals_first_fft_scratch(rsi, 64, 4*clmblkdst8, rdi); \
		}while(--loopB & (2*clm-1)); \
		newxloop1doublepre(clm, b1b, b1a, clmblkdst, rdi, XMM_SCD9); /* 2*clm iters prefetching 4*clm cache lines */ \
	} \
	rsi += -8*clmblkdst+clmblkdst8; \
} \
rsi += -4*clmblkdst8;	/* Restore source pointer */ \
end_timer(22); \
\
/* Do FFT levels 6,7 */ \
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
/* Do 8 eight_reals_fft operations */ \
/* Do 8 four_complex_fft operations */ \
/*	distance between fft data elements is 8K */ \
/*	do 8 macros each processing 32 data values */ \
/* Do 64 four_complex_fft operations */ \
/*	distance between fft data elements is 8K */ \
/*	do 32 macros each processing 32 data values */ \
\
start_timer(25); \
rdx = (uintptr_t)g->premult_addr;/* Address of premultiplier/sin/cos table */ \
rdx += 2*clm*XMM_SCD4 + 10*2*clm*XMM_SCD8; /* Real data sin/cos ptr */ \
rdi = (rdx+8*2*clm*XMM_SCND2R); /* Complex data sin/cos ptr */ \
loops_init(40*2*clm);		/* 5 iterations of 8 iters of 2*clm */ \
r4_x4cl_wpn_eight_reals_four_complex_djbfft_preload; \
for(unsigned int loopA = 8*2*clm; loopA; ) { \
	xtouch(xptr(rcx+rbp+(clm-1)*128));	/* Preload TLB */ \
	do{ \
		r4_x4cl_wpn_eight_reals_four_complex_djbfft(rsi, 64, clmblkdst8, 2*clmblkdst8, rdi, rdx); \
		rdi += XMM_SCND2;		/* Next sine/cosine pointer */ \
		rdx += XMM_SCND2R;		/* Next sine/cosine pointer */ \
	}while(--loopA & (2*clm-1)); \
	newxloop1(clm, b6b, clmblkdst);	/* 2*clm iterations */ \
} \
rsi += -8*clmblkdst+4*clmblkdst8; \
rdi += -8*2*clm*XMM_SCND2;	/* Restore sin/cos ptr */ \
\
loops_reset; \
r4_x4cl_wpn_four_complex_djbfft_preload; \
for(unsigned int loopA = 4; loopA; loopA--) { \
	for(unsigned int loopB = 8*2*clm; loopB; ) { \
		xtouch(xptr(rcx+rbp+(clm-1)*128));	/* Preload TLB */ \
		do{ \
			r4_x4cl_wpn_four_complex_djbfft(rsi, 64, clmblkdst8, 2*clmblkdst8, rdi); \
			rdi += XMM_SCND2;		/* Next sine/cosine pointer */ \
		}while(--loopB & (2*clm-1)); \
		newxloop1(clm, b7b, clmblkdst);	/* 2*clm iterations */ \
	} \
	rsi += -8*clmblkdst+4*clmblkdst8; /* Next source pointer */ \
	rdi += -8*2*clm*XMM_SCND2;	/* Restore sin/cos ptr */ \
} \
rsi += -5*4*clmblkdst8;	/* Restore source pointer */ \
end_timer(25); \
\
/* Do FFT levels 8,9,10 */ \
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
/* Do 38 eight_complex_fft8 operations */ \
/*	distance between fft data elements is 1K */ \
/*	do 19 macros each processing 64 data values */ \
\
start_timer(26); \
rdi += -8*2*clm*XMM_SCND2R-10*2*clm*XMM_SCD8+128; /* Complex data premult/sin/cos ptr and use -128 */ \
;					/* addressing to generate 1 byte offsets in radix-8 macro */ \
rdx = (rdi-128-2*clm*XMM_SCD4); /* Real data sin/cos ptr */ \
rcx = (uintptr_t)g->data_addr;/* Destination pointer */ \
loops_init_prefetch(20*2*clm, 128, 1, rbp); /* 20 iters of 2*clm */ \
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
for(unsigned int loopA = 10*2-1; loopA; ) { \
	do{ \
		for(unsigned int loopC = 2*clm; loopC; loopC--) { \
			r8_sg8cl_eight_complex_fft8(rsi, 64, clmblkdst, 2*clmblkdst, 4*clmblkdst, rcx, 64, blkdst, 2*blkdst, 4*blkdst, rdi-128); \
			rdi += XMM_SCD8;		/* Next premult/sin/cos ptr */ \
		} \
		rsi += -2*clm*64+clmblkdst8; /* Next source pointer */ \
		rcx += -2*clm*64+8*blkdst;	/* Next dest pointer */ \
		rdi += -2*clm*XMM_SCD8;	/* Restore premult/sine/cosine pointer */ \
	}while(--loopA & 1); \
	rdi += 2*clm*XMM_SCD8;	/* Next premult/sine/cosine pointer */ \
} \
end_timer(26); \
\
/* Work on next 4*clm sets of 640 values */ \
\
endlp:	xpass1_get_next_block(pass2, c0b, b0b, _wpn); /* Get next block and jump */ \
/* done: */ xfft_footer_rn; \
\
\
\
/* Perform a pass 1 FFT of length 768 doubles. */ \
/* Memory examples are for a 768K FFT with 10 levels done in pass 2. */ \
/* Sin/cos data needed: 32*2*clm*XMM_SCD1 + 32*2*clm*XMM_SCD1 + 8*2*clm*XMM_SCND2 + 8*2*clm*XMM_SCND2R + */ \
/*			2*clm*XMM_SCD4 + 4*2*clm*XMM_SCD8 = 78*2*clm 64-byte cache lines */ \
\

#define r4dwpn_pass1sc768(procname, xpass2name, pass2_levels, clm) /* (fftname, pass2_levels, clmarg) */ \
\
/* set_FFT_constants(fftname, 0, 0, pass2_levels, clmarg); */ \
xfft_header_rn(procname, xpass2name, pass2_levels, clm); \
\
/* Do pass 2 (the last 10 FFT levels if this is a 768K FFT) */ \
\
/* pass2: */start_timer(0); \
xpass2name(g); \
end_timer(0); \
if(g->ffttype == 1)		/* We're done if FFTing only */ \
	return; /* goto done; */ \
\
/* Do the inverse FFT */ \
\
pass1_inverse_fft_setup(c0b, procname); \
\
/* Do inverse FFT levels 8,9,10 */ \
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
/* Do 47 eight_complex_unfft8 operations */ \
/*	distance between fft data elements is 1K */ \
/*	do 47 macros each processing 32 data values */ \
\
c0b:	start_timer(16); \
rsi = (uintptr_t)g->data_addr;/* Get source address */ \
rcx = (uintptr_t)g->scratch_area;	/* Get address of scratch area */ \
rdx = (uintptr_t)g->premult_addr;/* Real data sin/cos pointer */ \
rdi = (rdx+2*clm*XMM_SCD4+128); /* Complex data sin/cos ptr and use -128 */ \
;					/* addressing to generate 1 byte offsets in radix-8 macro */ \
pfing(rbp = (rdx+(78-48)*2*clm*64)); /* Finish prefetching the sin/cos data */ \
loops_init_prefetch(48*2*clm, 64, 1, rbp); /* 48*2*clm iterations */ \
/*bug	xtouchpages rbp, 48*2*clm * 64 / 1 ;; Prefetch TLBs */ \
for(unsigned int loopA = 2*clm; loopA; loopA--) { \
	r8_sg4cl_sixteen_reals_unfft8(rsi, 64, blkdst, 2*blkdst, rcx, 64, clmblkdst, 2*clmblkdst, rdi-128, rdx); \
	rdi += XMM_SCD8;		/* Next premult/sin/cos data */ \
	rdx += XMM_SCD4;		/* Next sine/cosine pointer */ \
} \
rsi += -2*clm*64+4*blkdst;	/* Next source pointer */ \
rcx += -2*clm*64+4*clmblkdst; /* Next dest pointer */ \
rdi += -2*clm*XMM_SCD8;	/* Restore premult/sine/cosine pointer */ \
\
loops_reset; \
for(unsigned int loopA = 4*6*2-1; loopA; ) { \
	for(unsigned int loopB = 6; loopB; loopB--) { \
		do{ \
			for(unsigned int loopD = 2*clm; loopD; loopD--) { \
				r8_sg4cl_eight_complex_unfft8(rsi, 64, blkdst, 2*blkdst, rcx, 64, clmblkdst, 2*clmblkdst, rdi-128); \
				rdi += XMM_SCD8;		/* Next premult/sin/cos data */ \
			} \
			rsi += -2*clm*64+4*blkdst;	/* Next source pointer */ \
			rcx += -2*clm*64+4*clmblkdst; /* Next dest pointer */ \
			rdi += -2*clm*XMM_SCD8;	/* Restore premult/sin/cos data pointer */ \
		}while(--loopA & 1); \
		rcx += -2*4*clmblkdst+clmblkdst8; \
	} \
	rdi += 2*clm*XMM_SCD8;	/* Next premult/sin/cos pointer */ \
} \
rsi = (uintptr_t)g->scratch_area;	/* Restore source pointer */ \
rdx = (rdi-128);		/* Next real data sin/cos pointer */ \
end_timer(16); \
\
/* Do inverse FFT levels 6,7 */ \
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
/* Do 88 four_complex_unfft operations */ \
/*	distance between fft data elements is 8K */ \
/*	do 44 macros each processing 32 data values */ \
\
start_timer(17); \
rdi = (rdx+8*2*clm*XMM_SCND2R); /* Complex data sin/cos ptr */ \
loops_init(48*2*clm);		/* 48 iterations of 2*clm */ \
\
r4_x4cl_wpn_eight_reals_unfft_preload; \
for(unsigned int loopA = 4; loopA; loopA--) { \
	for(unsigned int loopB = 2*clm; loopB; loopB--) { \
		r4_x4cl_wpn_eight_reals_unfft(rsi, 64, 4*clmblkdst, clmblkdst8, rdi, 4*2*clm*XMM_SCND2, rdx, 4*2*clm*XMM_SCND2R); \
		rdi += XMM_SCND2;		/* Next sine/cosine pointer */ \
		rdx += XMM_SCND2R;		/* Next sine/cosine pointer */ \
	} \
	rsi += -2*clm*64+clmblkdst; /* Next source pointer */ \
} \
rsi += -4*clmblkdst+2*clmblkdst8; /* Next source pointer */ \
rdi += -4*2*clm*XMM_SCND2;	/* Restore sine/cosine pointer */ \
\
loops_reset; \
r4_x4cl_wpn_four_complex_djbunfft_preload; \
for(unsigned int loopA = 11; loopA; loopA--) { \
	for(unsigned int loopB = 4; loopB; loopB--) { \
		for(unsigned int loopC = 2*clm; loopC; loopC--) { \
			r4_x4cl_wpn_four_complex_djbunfft(rsi, 64, 4*clmblkdst, clmblkdst8, rdi, 4*2*clm*XMM_SCND2); \
			rdi += XMM_SCND2;		/* Next clm's sin/cos data */ \
		} \
		rsi += -2*clm*64+clmblkdst; /* Next source pointer */ \
	} \
	rsi += -4*clmblkdst+2*clmblkdst8; /* Next source pointer */ \
	rdi += -4*2*clm*XMM_SCND2;	/* Restore sine/cosine pointer */ \
} \
rsi += -12*2*clmblkdst8;	/* Restore source pointer */ \
rdx = (rdi+8*2*clm*XMM_SCND2); /* Next real data sin/cos pointer */ \
end_timer(17); \
\
/* Do inverse FFT levels 4,5 (actually 1.585 levels) */ \
/* On input the 64-byte cache lines hold these data values: */ \
/*	0K	16K	32K	48K */ \
/*	1K	... */ \
/*	... */ \
/*	15K */ \
/*	64K	... */ \
/*	... */ \
/* On output the 64-byte cache lines hold these data values: */ \
/*	0K	48K	96K	144K */ \
/*	1K	... */ \
/*	... */ \
/*	47K */ \
/*	192K	... */ \
/*	... */ \
\
/* Do 32 six_reals_unfft operations */ \
/*	distance between fft data elements is 32K */ \
/*	do 16 macros each processing 32 data values */ \
/* Do 96 three_complex_unfft operations */ \
/*	distance between fft data elements is 32K */ \
/*	do 48 macros each processing 32 data values */ \
\
start_timer(18); \
rdi = (rdx+32*2*clm*XMM_SCD1); /* Complex data sin/cos ptr */ \
loops_init(64*2*clm);		/* 64 iterations of 2*clm */ \
\
for(unsigned int loopA = 2; loopA; loopA--) { \
	for(unsigned int loopB = 8; loopB; loopB--) { \
		for(unsigned int loopC = 2*clm; loopC; loopC--) { \
			r3_x3cl_six_reals_unfft(rsi, 64, 2*clmblkdst8, rdi, 16*2*clm*XMM_SCD1, rdx, 16*2*clm*XMM_SCD1); \
			rdi += XMM_SCD1;		/* Next sine/cosine pointer */ \
			rdx += XMM_SCD1;		/* Next sine/cosine pointer */ \
		} \
		rsi += -2*clm*64+clmblkdst; /* Next source pointer */ \
	} \
	rsi += -8*clmblkdst+clmblkdst8; /* Next source pointer */ \
} \
rsi += -2*clmblkdst8+6*clmblkdst8; /* Next source pointer */ \
rdi += -2*8*2*clm*XMM_SCD1; /* Restore sine/cosine pointer */ \
\
loops_reset; \
r3_x3cl_three_complex_djbunfft_preload; \
for(unsigned int loopA = 3; loopA; loopA--) { \
	for(unsigned int loopB = 2; loopB; loopB--) { \
		for(unsigned int loopC = 8; loopC; loopC--) { \
			for(unsigned int loopD = 2*clm; loopD; loopD--) { \
				r3_x3cl_three_complex_djbunfft(rsi, 64, 2*clmblkdst8, rdi, 16*2*clm*XMM_SCD1); \
				rdi += XMM_SCD1;		/* Next clm's sin/cos data */ \
			} \
			rsi += -2*clm*64+clmblkdst; /* Next source pointer */ \
		} \
		rsi += -8*clmblkdst+clmblkdst8; /* Next source pointer */ \
	} \
	rsi += -2*clmblkdst8+6*clmblkdst8; /* Next source pointer */ \
	rdi += -2*8*2*clm*XMM_SCD1; /* Restore sine/cosine pointer */ \
} \
rsi += -4*6*clmblkdst8;	/* Restore source pointer */ \
end_timer(18); \
\
/* Do inverse FFT levels 1,2,3 */ \
/* On input the 64-byte cache lines hold these data values: */ \
/*	0K	48K	96K	144K */ \
/*	1K	... */ \
/*	... */ \
/*	47K */ \
/*	192K	... */ \
/*	... */ \
/* On output the 64-byte cache lines hold these data values: */ \
/*	0K	192K	384K	576K */ \
/*	1K	... */ \
/*	... */ \
/*	191K */ \
\
/* Do 96 eight_reals_last_unfft operations */ \
/*	distance between fft data elements is 96K */ \
/*	do 48 macros each processing 32 data values */ \
\
start_timer(20); \
rdi = (uintptr_t)g->sincos2;/* Load fixed sin/cos postmultipliers */ \
set_data_prefetch_ptrs;		/* Init prefetching pointers */ \
loops_init(48*2*clm);		/* 48 iterations of 2*clm */ \
r4_x4cl_eight_reals_last_unfft_preload; \
for(unsigned int loopA = 6; loopA; loopA--) { \
	for(unsigned int loopB = 8*2*clm; loopB; ) { \
		xtouch(xptr(rcx+rbp+(clm-1)*128));	/* Preload TLB */ \
		do{ \
			r4_x4cl_eight_reals_last_unfft(rsi, 64, 6*clmblkdst8, 12*clmblkdst8, rdi, 48*XMM_SCD3); \
		}while(--loopB & (2*clm-1)); \
		newxloop1(clm, c1b, clmblkdst, rdi, XMM_SCD3); /* 2*clm iterations */ \
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
copy_scratch_data_to_fft(192, clm); \
goto	endlp;			/* Skip FFT code if flag not set */ \
\
/* Fire up auxiliary threads to do pass 1 forward FFT */ \
\
/* pass1: */pass1_forward_fft_setup(b0b); \
\
/* Copy data from the FFT data area back to the scratch area */ \
\
b0b:	copy_fft_data_to_scratch(192, clm); \
\
/* Do FFT levels 1,2,3 */ \
/* */ \
/* On input the 64-byte cache lines hold these data values: */ \
/*	0K	192K	384K	576K */ \
/*	1K	... */ \
/*	... */ \
/*	191K */ \
/* On output the 64-byte cache lines hold these data values: */ \
/*	0K	96K	192K	288K */ \
/*	1K	... */ \
/*	... */ \
/*	95K	... */ \
/*	384K	... */ \
/*	... */ \
\
/* Do 96 eight_reals_first_fft operations */ \
/*	distance between fft data elements is 96K */ \
/*	do 96 iterations each processing 16 data values */ \
\
yes_fft:start_timer(22); \
rsi = (uintptr_t)g->scratch_area;	/* Get address of scratch area */ \
rdi = (uintptr_t)g->sincos2;/* Load fixed sin/cos postmultipliers */ \
set_data_prefetch_ptrs;		/* Init prefetching pointers */ \
pfing(rcx += (48-8)*blkdst);	/* We've already prefetched 48 blocks, but */ \
;					/* code below prefetches 8 blocks too many */ \
loops_init(96*2*clm);		/* 12 iters of 4 iters of 2*2*clm */ \
r4_x2cl_eight_reals_first_fft_scratch_preload; \
for(unsigned int loopA = 12; loopA; loopA--) { \
	for(unsigned int loopB = 4*2*2*clm; loopB; ) { \
		xtouch(xptr(rcx+rbp+(clm-1)*128));	/* Preload TLB */ \
		do{ \
			do{ \
				r4_x2cl_eight_reals_first_fft_scratch(rsi, 64, 12*clmblkdst8, rdi); \
			}while(--loopB & (2*clm-1)); \
			newxloop2(clm, b1b, clmblkdst, rdi, XMM_SCD3); /* 2*2*clm iters prefetching 2*clm cache lines */ \
		}while(loopB & (2*2*clm-1)); \
	} \
	rsi += -4*2*clmblkdst+clmblkdst8; \
} \
rsi += -12*clmblkdst8;	/* Restore source pointer */ \
end_timer(22); \
\
/* Do FFT levels 4,5 (actually 1.585 levels) */ \
/* */ \
/* On input the 64-byte cache lines hold these data values: */ \
/*	0K	96K	192K	288K */ \
/*	1K	... */ \
/*	... */ \
/*	95K	... */ \
/*	384K	... */ \
/*	... */ \
/* On output the 64-byte cache lines hold these data values: */ \
/*	0K	32K	64K	96K */ \
/*	1K	... */ \
/*	... */ \
/*	31K	... */ \
/*	128K	... */ \
/*	... */ \
\
/* Do 32 six_reals_fft operations */ \
/* Do 32 three_complex_fft operations */ \
/*	distance between fft data elements is 32K */ \
/*	do 32 macros each processing 32 data values */ \
/* Do 64 three_complex_fft operations */ \
/*	distance between fft data elements is 32K */ \
/*	do 32 macros each processing 32 data values */ \
\
start_timer(24); \
rdx = (uintptr_t)g->premult_addr;/* Address of premultiplier/sin/cos table */ \
rdx += 8*2*clm*XMM_SCND2 + 8*2*clm*XMM_SCND2R + 2*clm*XMM_SCD4 + 4*2*clm*XMM_SCD8; /* Real data sin/cos ptr */ \
rdi = (rdx+32*2*clm*XMM_SCD1); /* Complex data sin/cos ptr */ \
loops_init(64*2*clm);		/* 8 iters of 8 of 2*clm */ \
r3_x3cl_six_reals_three_complex_djbfft_preload; \
for(unsigned int loopA = 4; loopA; loopA--) { \
	for(unsigned int loopB = 8*2*clm; loopB; ) { \
		xtouch(xptr(rcx+rbp+(clm-1)*128));	/* Preload TLB */ \
		do{ \
			r3_x3cl_six_reals_three_complex_djbfft(rsi, 64, 4*clmblkdst8, rdi, rdx); \
			rdi += XMM_SCD1;		/* Next sine/cosine pointer */ \
			rdx += XMM_SCD1;		/* Next sine/cosine pointer */ \
		}while(--loopB & (2*clm-1)); \
		newxloop1(clm, b4b, clmblkdst);	/* 2*clm iterations */ \
	} \
	rsi += -8*clmblkdst+clmblkdst8; \
} \
rsi += -4*clmblkdst8+12*clmblkdst8; /* Next source pointer */ \
rdi += -32*2*clm*XMM_SCD1;	/* Restore sin/cos ptr */ \
\
loops_reset; \
r3_x3cl_three_complex_djbfft_preload; \
for(unsigned int loopA = 4; loopA; loopA--) { \
	for(unsigned int loopB = 8*2*clm; loopB; ) { \
		xtouch(xptr(rcx+rbp+(clm-1)*128));	/* Preload TLB */ \
		do{ \
			r3_x3cl_three_complex_djbfft(rsi, 64, 4*clmblkdst8, rdi); \
			rdi += XMM_SCD1;		/* Next sine/cosine pointer */ \
		}while(--loopB & (2*clm-1)); \
		newxloop1(clm, b5b, clmblkdst);	/* 2*clm iterations */ \
	} \
	rsi += -8*clmblkdst+clmblkdst8; \
} \
rsi += -4*clmblkdst8-12*clmblkdst8; /* Restore source pointer */ \
end_timer(24); \
\
/* Do FFT levels 6,7 */ \
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
/* Do 8 eight_reals_fft operations */ \
/* Do 8 four_complex_fft operations */ \
/*	distance between fft data elements is 8K */ \
/*	do 8 macros each processing 32 data values */ \
/* Do 80 four_complex_fft operations */ \
/*	distance between fft data elements is 8K */ \
/*	do 40 macros each processing 32 data values */ \
\
start_timer(25); \
rdi += -32*2*clm*XMM_SCD1-32*2*clm*XMM_SCD1-8*2*clm*XMM_SCND2; /* Complex data sin/cos ptr */ \
rdx = (rdi-8*2*clm*XMM_SCND2R); /* Real data sin/cos ptr */ \
loops_init(48*2*clm);		/* 6 iterations of 8 iters of 2*clm */ \
r4_x4cl_eight_reals_four_complex_djbfft_preload; \
for(unsigned int loopA = 8*2*clm; loopA; ) { \
	xtouch(xptr(rcx+rbp+(clm-1)*128));	/* Preload TLB */ \
	do{ \
		r4_x4cl_wpn_eight_reals_four_complex_djbfft(rsi, 64, clmblkdst8, 2*clmblkdst8, rdi, rdx); \
		rdi += XMM_SCND2;		/* Next sine/cosine pointer */ \
		rdx += XMM_SCND2R;		/* Next sine/cosine pointer */ \
	}while(--loopA & (2*clm-1)); \
	newxloop1nopre(clm, b6b, clmblkdst);	/* 2*clm iterations */ \
} \
rsi += -8*clmblkdst+4*clmblkdst8; \
rdi += -8*2*clm*XMM_SCND2;	/* Restore sin/cos ptr */ \
\
loops_reset; \
r4_x4cl_wpn_four_complex_djbfft_preload; \
for(unsigned int loopA = 5; loopA; loopA--) { \
	for(unsigned int loopB = 8*2*clm; loopB; ) { \
		xtouch(xptr(rcx+rbp+(clm-1)*128));	/* Preload TLB */ \
		do{ \
			r4_x4cl_wpn_four_complex_djbfft(rsi, 64, clmblkdst8, 2*clmblkdst8, rdi); \
			rdi += XMM_SCND2;		/* Next sine/cosine pointer */ \
		}while(--loopB & (2*clm-1)); \
		newxloop1(clm, b7b, clmblkdst);	/* 2*clm iterations */ \
	} \
	rsi += -8*clmblkdst+4*clmblkdst8; /* Next source pointer */ \
	rdi += -8*2*clm*XMM_SCND2;	/* Restore sin/cos ptr */ \
} \
rsi += -6*4*clmblkdst8;	/* Restore source pointer */ \
end_timer(25); \
\
/* Do FFT levels 8,9,10 */ \
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
/* Do 46 eight_complex_fft8 operations */ \
/*	distance between fft data elements is 1K */ \
/*	do 23 macros each processing 64 data values */ \
\
start_timer(26); \
rdi += -8*2*clm*XMM_SCND2R-4*2*clm*XMM_SCD8+128; /* Complex data sin/cos ptr and use -128 */ \
;					/* addressing to generate 1 byte offsets in radix-8 macro */ \
rdx = (rdi-128-2*clm*XMM_SCD4); /* Real data sin/cos ptr */ \
rcx = (uintptr_t)g->data_addr;/* Destination pointer */ \
loops_init_prefetch(24*2*clm, 128, 1, rbp); /* 24 iters of 2*clm */ \
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
for(unsigned int loopA = 4*(6+2)-1; loopA; ) { \
	loopA -= 2; \
	do{ /* 5 or 6 loops */ \
		for(unsigned int loopC = 2*clm; loopC; loopC--) { \
			r8_sg8cl_eight_complex_fft8(rsi, 64, clmblkdst, 2*clmblkdst, 4*clmblkdst, rcx, 64, blkdst, 2*blkdst, 4*blkdst, rdi-128); \
			rdi += XMM_SCD8;		/* Next premult/sin/cos ptr */ \
		} \
		rsi += -2*clm*64+clmblkdst8; /* Next source pointer */ \
		rcx += -2*clm*64+8*blkdst;	/* Next dest pointer */ \
		rdi += -2*clm*XMM_SCD8;	/* Restore premult/sine/cosine pointer */ \
	}while(--loopA & 7); \
	rdi += 2*clm*XMM_SCD8;	/* Next premult/sine/cosine pointer */ \
} \
end_timer(26); \
\
/* Work on next 4*clm sets of 768 values */ \
\
endlp:	xpass1_get_next_block(pass2, c0b, b0b, _wpn); /* Get next block and jump */ \
/* done: */ xfft_footer_rn; \
\
\
\
/* Perform a pass 1 FFT of length 896 doubles. */ \
/* Memory examples are for a 896K FFT with 10 levels done in pass 2. */ \
/* Sin/cos data needed: 8*2*clm*XMM_SCND2 + 8*2*clm*XMM_SCND2R + 2*clm*XMM_SCD4 + 14*2*clm*XMM_SCD8 = 86*2*clm 64-byte cache lines */ \
\

#define r4dwpn_pass1sc896(procname, xpass2name, pass2_levels, clm) /* (fftname, pass2_levels, clmarg) */ \
\
/* set_FFT_constants(fftname, 0, 0, pass2_levels, clmarg); */ \
xfft_header_rn(procname, xpass2name, pass2_levels, clm); \
\
/* Do pass 2 (the last 10 FFT levels if this is a 448K FFT) */ \
\
/* pass2: */start_timer(0); \
xpass2name(g); \
end_timer(0); \
if(g->ffttype == 1)		/* We're done if FFTing only */ \
	return; /* goto done; */ \
\
/* Do the inverse FFT */ \
\
pass1_inverse_fft_setup(c0b, procname); \
\
/* Do inverse FFT levels 8,9,10 */ \
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
/* Do 55 eight_complex_unfft8 operations */ \
/*	distance between fft data elements is 1K */ \
/*	do 55 macros each processing 32 data values */ \
\
c0b:	start_timer(16); \
rsi = (uintptr_t)g->data_addr;/* Get source address */ \
rcx = (uintptr_t)g->scratch_area;	/* Get address of scratch area */ \
rdx = (uintptr_t)g->premult_addr;/* Real data sin/cos pointer */ \
rdi = (rdx+2*clm*XMM_SCD4+128); /* Complex data sin/cos ptr and use -128 */ \
;					/* addressing to generate 1 byte offsets in radix-8 macro */ \
pfing(rbp = (rdx+(86-56)*2*clm*64)); /* Finish prefetching the sin/cos data */ \
loops_init_prefetch(56*2*clm, 64, 1, rbp); /* 56*2*clm iterations */ \
/*bug	xtouchpages rbp, 56*2*clm * 64 / 1 ;; Prefetch TLBs */ \
for(unsigned int loopA = 2*clm; loopA; loopA--) { \
	r8_sg4cl_sixteen_reals_unfft8(rsi, 64, blkdst, 2*blkdst, rcx, 64, clmblkdst, 2*clmblkdst, rdi-128, rdx); \
	rdi += XMM_SCD8;		/* Next premult/sin/cos data */ \
	rdx += XMM_SCD4;		/* Next sine/cosine pointer */ \
} \
rsi += -2*clm*64+4*blkdst;	/* Next source pointer */ \
rcx += -2*clm*64+4*clmblkdst; /* Next dest pointer */ \
rdi += -2*clm*XMM_SCD8;	/* Restore premult/sine/cosine pointer */ \
\
loops_reset; \
for(unsigned int loopA = 14*2*2-1; loopA; ) { \
	do{ \
		do{ \
			for(unsigned int loopD = 2*clm; loopD; loopD--) { \
				r8_sg4cl_eight_complex_unfft8(rsi, 64, blkdst, 2*blkdst, rcx, 64, clmblkdst, 2*clmblkdst, rdi-128); \
				rdi += XMM_SCD8;		/* Next premult/sin/cos data */ \
			} \
			rsi += -2*clm*64+4*blkdst;	/* Next source pointer */ \
			rcx += -2*clm*64+4*clmblkdst; /* Next dest pointer */ \
			rdi += -2*clm*XMM_SCD8;	/* Restore premult/sin/cos data pointer */ \
		}while(--loopA & 1); \
		rcx += -2*4*clmblkdst+clmblkdst8; \
	}while(loopA & (2*2-1)); \
	rdi += 2*clm*XMM_SCD8;	/* Next premult/sin/cos pointer */ \
} \
rsi = (uintptr_t)g->scratch_area;	/* Restore source pointer */ \
rdx = (rdi-128);		/* Next real data sin/cos pointer */ \
end_timer(16); \
\
/* Do inverse FFT levels 6,7 */ \
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
/* Do 104 four_complex_unfft operations */ \
/*	distance between fft data elements is 8K */ \
/*	do 52 macros each processing 32 data values */ \
\
start_timer(17); \
rdi = (rdx+8*2*clm*XMM_SCND2R); /* Complex data sin/cos ptr */ \
loops_init(56*2*clm);		/* 56 iterations of 2*clm */ \
\
r4_x4cl_wpn_eight_reals_unfft_preload; \
for(unsigned int loopA = 4; loopA; loopA--) { \
	for(unsigned int loopB = 2*clm; loopB; loopB--) { \
		r4_x4cl_wpn_eight_reals_unfft(rsi, 64, 4*clmblkdst, clmblkdst8, rdi, 4*2*clm*XMM_SCND2, rdx, 4*2*clm*XMM_SCND2R); \
		rdi += XMM_SCND2;		/* Next sine/cosine pointer */ \
		rdx += XMM_SCND2R;		/* Next sine/cosine pointer */ \
	} \
	rsi += -2*clm*64+clmblkdst; /* Next source pointer */ \
} \
rsi += -4*clmblkdst+2*clmblkdst8; /* Next source pointer */ \
rdi += -4*2*clm*XMM_SCND2;	/* Restore sine/cosine pointer */ \
\
loops_reset; \
r4_x4cl_wpn_four_complex_djbunfft_preload; \
for(unsigned int loopA = 13; loopA; loopA--) { \
	for(unsigned int loopB = 4; loopB; loopB--) { \
		for(unsigned int loopC = 2*clm; loopC; loopC--) { \
			r4_x4cl_wpn_four_complex_djbunfft(rsi, 64, 4*clmblkdst, clmblkdst8, rdi, 4*2*clm*XMM_SCND2); \
			rdi += XMM_SCND2;		/* Next clm's sin/cos data */ \
		} \
		rsi += -2*clm*64+clmblkdst; /* Next source pointer */ \
	} \
	rsi += -4*clmblkdst+2*clmblkdst8; /* Next source pointer */ \
	rdi += -4*2*clm*XMM_SCND2;	/* Restore sine/cosine pointer */ \
} \
rsi += -14*2*clmblkdst8;	/* Restore source pointer */ \
end_timer(17); \
\
/* Do inverse FFT levels 1,2,3,4,5 (actually 4.807 levels) */ \
/* On input the 64-byte cache lines hold these data values: */ \
/*	0K	16K	32K	48K */ \
/*	1K	... */ \
/*	... */ \
/*	15K */ \
/*	64K	... */ \
/*	... */ \
/* On output the 64-byte cache lines hold these data values: */ \
/*	0K	224K	448K	672K */ \
/*	1K	... */ \
/*	... */ \
/*	223K */ \
\
/* Do 128 28_reals_last_unfft operations */ \
/*	distance between fft data elements is 32K */ \
/*	do 16 macros each processing 112 data values */ \
\
start_timer(20); \
rdi = (uintptr_t)g->sincos2;/* Load fixed sin/cos postmultipliers */ \
set_data_prefetch_ptrs;		/* Init prefetching pointers */ \
loops_init(16*2*clm);		/* 16 iterations of 2*clm */ \
r7_x14cl_28_reals_last_unfft_preload; \
for(unsigned int loopA = 2; loopA; loopA--) { \
	for(unsigned int loopB = 8*2*clm; loopB; ) { \
		xtouch(xptr(rcx+rbp+(clm-1)*128));	/* Preload TLB */ \
		do{ \
			r7_x14cl_28_reals_last_unfft(rsi, 64, 2*clmblkdst8, rdi, 16*XMM_SCD13); \
		}while(--loopB & (2*clm-1)); \
		newxloop1quadpre(clm, c1b, c1a, clmblkdst, rdi, XMM_SCD13); /* 2*clm iterations */ \
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
copy_scratch_data_to_fft(224, clm); \
goto	endlp;			/* Skip FFT code if flag not set */ \
\
/* Fire up auxiliary threads to do pass 1 forward FFT */ \
\
/* pass1: */pass1_forward_fft_setup(b0b); \
\
/* Copy data from the FFT data area back to the scratch area */ \
\
b0b:	copy_fft_data_to_scratch(224, clm); \
\
/* Do FFT levels 1,2,3,4,5 (actually 4.807 levels) */ \
/* */ \
/* On input the 64-byte cache lines hold these data values: */ \
/*	0K	224K	448K	672K */ \
/*	1K	... */ \
/*	... */ \
/*	223K */ \
/* On output the 64-byte cache lines hold these data values: */ \
/*	0K	32K	64K	96K */ \
/*	1K	... */ \
/*	... */ \
/*	31K	... */ \
/*	128K	... */ \
/*	... */ \
\
/* Do 128 28_reals_first_fft operations */ \
/*	distance between fft data elements is 32K */ \
/*	do 32 iterations each processing 56 data values */ \
\
yes_fft:start_timer(22); \
rsi = (uintptr_t)g->scratch_area;	/* Get address of scratch area */ \
rdi = (uintptr_t)g->sincos2;/* Load fixed sin/cos postmultipliers */ \
set_data_prefetch_ptrs;		/* Init prefetching pointers */ \
pfing(rcx += (64-8)*blkdst);	/* We've already prefetched 64 blocks, but */ \
;					/* code below will fetch 8 blocks too many */ \
loops_init(32*2*clm);		/* 32 iters of 2*clm */ \
r7_x7cl_28_reals_first_fft_scratch_preload; \
for(unsigned int loopA = 4; loopA; loopA--) { \
	for(unsigned int loopB = 8*2*clm; loopB; ) { \
		xtouch(xptr(rcx+rbp+(clm-1)*128));	/* Preload TLB */ \
		do{ \
			r7_x7cl_28_reals_first_fft_scratch(rsi, 64, 4*clmblkdst8, rdi); \
		}while(--loopB & (2*clm-1)); \
		newxloop1doublepre(clm, b1b, b1a, clmblkdst, rdi, XMM_SCD13); /* 2*clm iters prefetching 4*clm cache lines */ \
	} \
	rsi += -8*clmblkdst+clmblkdst8; \
} \
rsi += -4*clmblkdst8;	/* Restore source pointer */ \
end_timer(22); \
\
/* Do FFT levels 6,7 */ \
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
/* Do 8 eight_reals_fft operations */ \
/* Do 8 four_complex_fft operations */ \
/*	distance between fft data elements is 8K */ \
/*	do 8 macros each processing 32 data values */ \
/* Do 96 four_complex_fft operations */ \
/*	distance between fft data elements is 8K */ \
/*	do 48 macros each processing 32 data values */ \
\
start_timer(25); \
rdx = (uintptr_t)g->premult_addr;/* Address of premultiplier/sin/cos table */ \
rdx += 2*clm*XMM_SCD4 + 14*2*clm*XMM_SCD8; /* Real data sin/cos ptr */ \
rdi = (rdx+8*2*clm*XMM_SCND2R); /* Complex data sin/cos ptr */ \
loops_init(56*2*clm);		/* 7 iterations of 8 iters of 2*clm */ \
r4_x4cl_wpn_eight_reals_four_complex_djbfft_preload; \
for(unsigned int loopA = 8*2*clm; loopA; ) { \
	xtouch(xptr(rcx+rbp+(clm-1)*128));	/* Preload TLB */ \
	do{ \
		r4_x4cl_wpn_eight_reals_four_complex_djbfft(rsi, 64, clmblkdst8, 2*clmblkdst8, rdi, rdx); \
		rdi += XMM_SCND2;		/* Next sine/cosine pointer */ \
		rdx += XMM_SCND2R;		/* Next sine/cosine pointer */ \
	}while(--loopA & (2*clm-1)); \
	newxloop1(clm, b6b, clmblkdst);	/* 2*clm iterations */ \
} \
rsi += -8*clmblkdst+4*clmblkdst8; \
rdi += -8*2*clm*XMM_SCND2;	/* Restore sin/cos ptr */ \
\
loops_reset; \
r4_x4cl_wpn_four_complex_djbfft_preload; \
for(unsigned int loopA = 6; loopA; loopA--) { \
	for(unsigned int loopB = 8*2*clm; loopB; ) { \
		xtouch(xptr(rcx+rbp+(clm-1)*128));	/* Preload TLB */ \
		do{ \
			r4_x4cl_wpn_four_complex_djbfft(rsi, 64, clmblkdst8, 2*clmblkdst8, rdi); \
			rdi += XMM_SCND2;		/* Next sine/cosine pointer */ \
		}while(--loopB & (2*clm-1)); \
		newxloop1doublepre(clm, b7b, b7a, clmblkdst); /* 2*clm iterations */ \
	} \
	rsi += -8*clmblkdst+4*clmblkdst8; /* Next source pointer */ \
	rdi += -8*2*clm*XMM_SCND2;	/* Restore sin/cos ptr */ \
} \
rsi += -7*4*clmblkdst8;	/* Restore source pointer */ \
end_timer(25); \
\
/* Do FFT levels 8,9,10 */ \
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
/* Do 54 eight_complex_fft8 operations */ \
/*	distance between fft data elements is 1K */ \
/*	do 27 macros each processing 64 data values */ \
\
start_timer(26); \
rdi += -8*2*clm*XMM_SCND2R-14*2*clm*XMM_SCD8+128; /* Complex data sin/cos ptr and use -128 */ \
;					/* addressing to generate 1 byte offsets in radix-8 macro */ \
rdx = (rdi-128-2*clm*XMM_SCD4); /* Real data sin/cos ptr */ \
rcx = (uintptr_t)g->data_addr;/* Destination pointer */ \
loops_init_prefetch(28*2*clm, 128, 1, rbp); /* 28 iters of 2*clm */ \
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
for(unsigned int loopA = 14*2-1; loopA; ) { \
	do{ \
		for(unsigned int loopC = 2*clm; loopC; loopC--) { \
			r8_sg8cl_eight_complex_fft8(rsi, 64, clmblkdst, 2*clmblkdst, 4*clmblkdst, rcx, 64, blkdst, 2*blkdst, 4*blkdst, rdi-128); \
			rdi += XMM_SCD8;		/* Next premult/sin/cos ptr */ \
		} \
		rsi += -2*clm*64+clmblkdst8; /* Next source pointer */ \
		rcx += -2*clm*64+8*blkdst;	/* Next dest pointer */ \
		rdi += -2*clm*XMM_SCD8;	/* Restore premult/sine/cosine pointer */ \
	}while(--loopA & 1); \
	rdi += 2*clm*XMM_SCD8;	/* Next premult/sine/cosine pointer */ \
} \
end_timer(26); \
\
/* Work on next 4*clm sets of 896 values */ \
\
endlp:	xpass1_get_next_block(pass2, c0b, b0b, _wpn); /* Get next block and jump */ \
/* done: */ xfft_footer_rn; \
\
\
\
/* Perform a pass 1 FFT of length 1024 doubles. */ \
/* Memory examples are for a 1M FFT with 10 levels done in pass 2. */ \
/* Sin/cos data needed: 8*2*clm*XMM_SCND2 + 8*2*clm*XMM_SCND2R + 2*clm*XMM_SCD4 + 16*2*clm*XMM_SCD8 = 94*2*clm 64-byte cache lines */ \
\

#define r4dwpn_pass1sc1024(procname, xpass2name, pass2_levels, clm) /* (fftname, pass2_levels, clmarg) */ \
\
/* set_FFT_constants(fftname, 0, 0, pass2_levels, clmarg); */ \
xfft_header_rn(procname, xpass2name, pass2_levels, clm); \
\
/* Do pass 2 (the last 10 FFT levels if this is a 1M FFT) */ \
\
/* pass2: */start_timer(0); \
xpass2name(g); \
end_timer(0); \
if(g->ffttype == 1)		/* We're done if FFTing only */ \
	return; /* goto done; */ \
\
/* Do the inverse FFT */ \
\
pass1_inverse_fft_setup(c0b, procname); \
\
/* Do inverse FFT levels 8,9,10 */ \
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
/* Do 63 eight_complex_unfft8 operations */ \
/*	distance between fft data elements is 1K */ \
/*	do 63 macros each processing 32 data values */ \
\
c0b:	start_timer(16); \
rsi = (uintptr_t)g->data_addr;/* Get source address */ \
rcx = (uintptr_t)g->scratch_area;	/* Get address of scratch area */ \
rdx = (uintptr_t)g->premult_addr;/* Real data sin/cos pointer */ \
rdi = (rdx+2*clm*XMM_SCD4+128); /* Complex data sin/cos ptr and use -128 */ \
;					/* addressing to generate 1 byte offsets in radix-8 macro */ \
pfing(rbp = (rdx+(94-32)*2*clm*64)); /* Finish prefetching the sin/cos data */ \
loops_init_prefetch(64*2*clm, 64, 2, rbp); /* 64*2*clm iterations */ \
/*bug	xtouchpages rbp, 64*2*clm * 64 / 2 ;; Prefetch TLBs */ \
for(unsigned int loopA = 2*clm; loopA; loopA--) { \
	r8_sg4cl_sixteen_reals_unfft8(rsi, 64, blkdst, 2*blkdst, rcx, 64, clmblkdst, 2*clmblkdst, rdi-128, rdx); \
	rdi += XMM_SCD8;		/* Next premult/sin/cos data */ \
	rdx += XMM_SCD4;		/* Next sine/cosine pointer */ \
} \
rsi += -2*clm*64+4*blkdst;	/* Next source pointer */ \
rcx += -2*clm*64+4*clmblkdst; /* Next dest pointer */ \
rdi += -2*clm*XMM_SCD8;	/* Restore premult/sine/cosine pointer */ \
\
loops_reset; \
/* c9a: */ xtouch(xptr(rbp+4096));		/* Prefetch TLB */ \
for(unsigned int loopA = 16*2*2-1; loopA; ) { \
	do{ \
		do{ \
			for(unsigned int loopD = 2*clm; loopD; loopD--) { \
				r8_sg4cl_eight_complex_unfft8(rsi, 64, blkdst, 2*blkdst, rcx, 64, clmblkdst, 2*clmblkdst, rdi-128); \
				rdi += XMM_SCD8;		/* Next premult/sin/cos data */ \
			} \
			rsi += -2*clm*64+4*blkdst;	/* Next source pointer */ \
			rcx += -2*clm*64+4*clmblkdst; /* Next dest pointer */ \
			rdi += -2*clm*XMM_SCD8;	/* Restore premult/sin/cos data pointer */ \
		}while(--loopA & 1); \
		rcx += -2*4*clmblkdst+clmblkdst8; \
	}while(loopA & (2*2-1)); \
	rdi += 2*clm*XMM_SCD8;	/* Next premult/sin/cos pointer */ \
} \
rsi = (uintptr_t)g->scratch_area;	/* Restore source pointer */ \
rdx = (rdi-128);		/* Next real sin/cos data pointer */ \
end_timer(16); \
\
/* Do inverse FFT levels 6,7 */ \
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
/* Do 120 four_complex_unfft operations */ \
/*	distance between fft data elements is 8K */ \
/*	do 60 macros each processing 32 data values */ \
\
start_timer(17); \
rdi = (rdx+8*2*clm*XMM_SCND2R); /* Complex data sin/cos ptr */ \
loops_init(64*2*clm);		/* 64 iterations of 2*clm */ \
\
r4_x4cl_wpn_eight_reals_unfft_preload; \
for(unsigned int loopA = 4; loopA; loopA--) { \
	for(unsigned int loopB = 2*clm; loopB; loopB--) { \
		r4_x4cl_wpn_eight_reals_unfft(rsi, 64, 4*clmblkdst, clmblkdst8, rdi, 4*2*clm*XMM_SCND2, rdx, 4*2*clm*XMM_SCND2R); \
		rdi += XMM_SCND2;		/* Next sine/cosine pointer */ \
		rdx += XMM_SCND2R;		/* Next sine/cosine pointer */ \
	} \
	rsi += -2*clm*64+clmblkdst; /* Next source pointer */ \
} \
rsi += -4*clmblkdst+2*clmblkdst8; /* Next source pointer */ \
rdi += -4*2*clm*XMM_SCND2;	/* Restore sine/cosine pointer */ \
\
loops_reset; \
r4_x4cl_wpn_four_complex_djbunfft_preload; \
for(unsigned int loopA = 15; loopA; loopA--) { \
	for(unsigned int loopB = 4; loopB; loopB--) { \
		for(unsigned int loopC = 2*clm; loopC; loopC--) { \
			r4_x4cl_wpn_four_complex_djbunfft(rsi, 64, 4*clmblkdst, clmblkdst8, rdi, 4*2*clm*XMM_SCND2); \
			rdi += XMM_SCND2;		/* Next clm's sin/cos data */ \
		} \
		rsi += -2*clm*64+clmblkdst; /* Next source pointer */ \
	} \
	rsi += -4*clmblkdst+2*clmblkdst8; /* Next source pointer */ \
	rdi += -4*2*clm*XMM_SCND2;	/* Restore sine/cosine pointer */ \
} \
rsi += -16*2*clmblkdst8;	/* Restore source pointer */ \
end_timer(17); \
\
/* Do inverse FFT levels 4,5 */ \
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
/*	256K	... */ \
/*	... */ \
\
/* Do 32 eight_reals_unfft operations */ \
/*	distance between fft data elements is 32K */ \
/*	do 16 macros each processing 32 data values */ \
/* Do 96 four_complex_unfft operations */ \
/*	distance between fft data elements is 32K */ \
/*	do 48 macros each processing 32 data values */ \
\
start_timer(18); \
rdx = (uintptr_t)g->sincos2;/* Load fixed sin/cos data ptr */ \
rdx += 128*XMM_SCD3;	/* Real data sin/cos ptr */ \
rdi = (rdx+32*XMM_SCD2);	/* Complex data sin/cos ptr */ \
loops_init(64*2*clm);		/* 64 iterations of 2*clm */ \
\
r4_x4cl_eight_reals_unfft_preload; \
for(unsigned int loopA = 2; loopA; loopA--) { \
	for(unsigned int loopB = 8; loopB; loopB--) { \
		for(unsigned int loopC = 2*clm; loopC; loopC--) { \
			r4_x4cl_eight_reals_unfft(rsi, 64, 2*clmblkdst8, 4*clmblkdst8, rdi, 16*XMM_SCD2, rdx, 16*XMM_SCD2); \
		} \
		rsi += -2*clm*64+clmblkdst; /* Next source pointer */ \
		rdi += XMM_SCD2;		/* Next sine/cosine pointer */ \
		rdx += XMM_SCD2;		/* Next sine/cosine pointer */ \
	} \
	rsi += -8*clmblkdst+clmblkdst8; /* Next source pointer */ \
} \
rsi += -2*clmblkdst8+8*clmblkdst8; /* Next source pointer */ \
rdi += -2*8*XMM_SCD2;	/* Restore sine/cosine pointer */ \
\
loops_reset; \
r4_x4cl_four_complex_djbunfft_preload; \
for(unsigned int loopA = 3; loopA; loopA--) { \
	for(unsigned int loopB = 2; loopB; loopB--) { \
		for(unsigned int loopC = 8; loopC; loopC--) { \
			for(unsigned int loopD = 2*clm; loopD; loopD--) { \
				r4_x4cl_four_complex_djbunfft(rsi, 64, 2*clmblkdst8, 4*clmblkdst8, rdi, 16*XMM_SCD2); \
			} \
			rsi += -2*clm*64+clmblkdst; /* Next source pointer */ \
			rdi += XMM_SCD2;		/* Next clm's sin/cos data */ \
		} \
		rsi += -8*clmblkdst+clmblkdst8; /* Next source pointer */ \
	} \
	rsi += -2*clmblkdst8+8*clmblkdst8; /* Next source pointer */ \
	rdi += -2*8*XMM_SCD2;	/* Restore sine/cosine pointer */ \
} \
rsi += -4*8*clmblkdst8;	/* Restore source pointer */ \
end_timer(18); \
\
/* Do inverse FFT levels 1,2,3 */ \
/* On input the 64-byte cache lines hold these data values: */ \
/*	0K	64K	128K	192K */ \
/*	1K	... */ \
/*	... */ \
/*	63K */ \
/*	256K	... */ \
/*	... */ \
/* On output the 64-byte cache lines hold these data values: */ \
/*	0K	256K	512K	768K */ \
/*	1K	... */ \
/*	... */ \
/*	255K */ \
\
/* Do 128 eight_reals_last_unfft operations */ \
/*	distance between fft data elements is 128K */ \
/*	do 64 macros each processing 32 data values */ \
\
start_timer(20); \
rdi = (uintptr_t)g->sincos2;/* Load fixed sin/cos postmultipliers */ \
set_data_prefetch_ptrs;		/* Init prefetching pointers */ \
loops_init(64*2*clm);		/* 64 iterations of 2*clm */ \
r4_x4cl_eight_reals_last_unfft_preload; \
for(unsigned int loopA = 8; loopA; loopA--) { \
	for(unsigned int loopB = 8*2*clm; loopB; ) { \
		xtouch(xptr(rcx+rbp+(clm-1)*128));	/* Preload TLB */ \
		do{ \
			r4_x4cl_eight_reals_last_unfft(rsi, 64, 8*clmblkdst8, 16*clmblkdst8, rdi, 64*XMM_SCD3); \
		}while(--loopB & (2*clm-1)); \
		newxloop1(clm, c1b, clmblkdst, rdi, XMM_SCD3); /* 2*clm iterations */ \
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
copy_scratch_data_to_fft(256, clm); \
goto	endlp;			/* Skip FFT code if flag not set */ \
\
/* Fire up auxiliary threads to do pass 1 forward FFT */ \
\
/* pass1: */pass1_forward_fft_setup(b0b); \
\
/* Copy data from the FFT data area back to the scratch area */ \
\
b0b:	copy_fft_data_to_scratch(256, clm); \
\
/* Do FFT levels 1,2,3 */ \
/* */ \
/* On input the 64-byte cache lines hold these data values: */ \
/*	0K	256K	512K	768K */ \
/*	1K	... */ \
/*	... */ \
/*	255K */ \
/* On output the 64-byte cache lines hold these data values: */ \
/*	0K	128K	256K	384K */ \
/*	1K	... */ \
/*	... */ \
/*	127K	... */ \
/*	512K	... */ \
/*	... */ \
\
/* Do 128 eight_reals_first_fft operations */ \
/*	distance between fft data elements is 128K */ \
/*	do 128 iterations each processing 16 data values */ \
\
yes_fft:start_timer(22); \
rsi = (uintptr_t)g->scratch_area;	/* Get address of scratch area */ \
rdi = (uintptr_t)g->sincos2;/* Load fixed sin/cos postmultipliers */ \
set_data_prefetch_ptrs;		/* Init prefetching pointers */ \
pfing(rcx += 64*blkdst);		/* We've already prefetched 64 blocks */ \
loops_init(128*2*clm);		/* 16 iters of 4 iters of 2*2*clm */ \
r4_x2cl_eight_reals_first_fft_scratch_preload; \
for(unsigned int loopA = 16; loopA; loopA--) { \
	for(unsigned int loopB = 4*2*2*clm; loopB; ) { \
		xtouch(xptr(rcx+rbp+(clm-1)*128));	/* Preload TLB */ \
		do{ \
			do{ \
				r4_x2cl_eight_reals_first_fft_scratch(rsi, 64, 16*clmblkdst8, rdi); \
			}while(--loopB & (2*clm-1)); \
			newxloop2(clm, b1b, clmblkdst, rdi, XMM_SCD3); /* 2*2*clm iters prefetching 2*clm cache lines */ \
		}while(loopB & (2*2*clm-1)); \
	} \
	rsi += -4*2*clmblkdst+clmblkdst8; \
} \
rsi += -16*clmblkdst8;	/* Restore source pointer */ \
rdx = rdi;		/* Next real data sin/cos ptr */ \
end_timer(22); \
\
/* Do FFT levels 4,5 */ \
/* */ \
/* On input the 64-byte cache lines hold these data values: */ \
/*	0K	128K	256K	384K */ \
/*	1K	... */ \
/*	... */ \
/*	127K	... */ \
/*	512K	... */ \
/*	... */ \
/* On output the 64-byte cache lines hold these data values: */ \
/*	0K	32K	64K	96K */ \
/*	1K	... */ \
/*	... */ \
/*	31K	... */ \
/*	128K	... */ \
/*	... */ \
\
/* Do 32 eight_reals_fft operations */ \
/* Do 32 four_complex_fft operations */ \
/*	distance between fft data elements is 32K */ \
/*	do 32 macros each processing 32 data values */ \
/* Do 64 four_complex_fft operations */ \
/*	distance between fft data elements is 32K */ \
/*	do 32 macros each processing 32 data values */ \
\
start_timer(24); \
rdi = (rdx+32*XMM_SCD2);	/* Complex data sin/cos ptr */ \
loops_init(64*2*clm);		/* 8 iters of 8 of 2*clm */ \
r4_x4cl_eight_reals_four_complex_djbfft_preload; \
for(unsigned int loopA = 4; loopA; loopA--) { \
	for(unsigned int loopB = 8*2*clm; loopB; ) { \
		xtouch(xptr(rcx+rbp+(clm-1)*128));	/* Preload TLB */ \
		do{ \
			r4_x4cl_eight_reals_four_complex_djbfft(rsi, 64, 4*clmblkdst8, 8*clmblkdst8, rdi, rdx); \
		}while(--loopB & (2*clm-1)); \
		newxloop1(clm, b4b, clmblkdst, rdi, XMM_SCD2, rdx, XMM_SCD2); /* 2*clm iterations */ \
	} \
	rsi += -8*clmblkdst+clmblkdst8; \
} \
rsi += -4*clmblkdst8+16*clmblkdst8; /* Next source pointer */ \
rdi += -32*XMM_SCD2;	/* Restore sin/cos ptr */ \
\
loops_reset; \
r4_x4cl_four_complex_djbfft_preload; \
for(unsigned int loopA = 4; loopA; loopA--) { \
	for(unsigned int loopB = 8*2*clm; loopB; ) { \
		xtouch(xptr(rcx+rbp+(clm-1)*128));	/* Preload TLB */ \
		do{ \
			r4_x4cl_four_complex_djbfft(rsi, 64, 4*clmblkdst8, 8*clmblkdst8, rdi); \
		}while(--loopB & (2*clm-1)); \
		newxloop1(clm, b5b, clmblkdst, rdi, XMM_SCD2); /* 2*clm iterations */ \
	} \
	rsi += -8*clmblkdst+clmblkdst8; \
} \
rsi += -4*clmblkdst8-16*clmblkdst8; /* Restore source pointer */ \
end_timer(24); \
\
/* Do FFT levels 6,7 */ \
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
/* Do 8 eight_reals_fft operations */ \
/* Do 8 four_complex_fft operations */ \
/*	distance between fft data elements is 8K */ \
/*	do 8 macros each processing 32 data values */ \
/* Do 112 four_complex_fft operations */ \
/*	distance between fft data elements is 8K */ \
/*	do 56 macros each processing 32 data values */ \
\
start_timer(25); \
rdx = (uintptr_t)g->premult_addr;/* Address of premultiplier/sin/cos table */ \
rdx += 2*clm*XMM_SCD4 + 16*2*clm*XMM_SCD8; /* Normalization ptr */ \
rdi = (rdx+8*2*clm*XMM_SCND2R); /* Complex data sin/cos ptr */ \
loops_init(64*2*clm);		/* 8 iterations of 8 iters of 2*clm */ \
r4_x4cl_wpn_eight_reals_four_complex_djbfft_preload; \
for(unsigned int loopA = 8*2*clm; loopA; ) { \
	xtouch(xptr(rcx+rbp+(clm-1)*128));	/* Preload TLB */ \
	do{ \
		r4_x4cl_wpn_eight_reals_four_complex_djbfft(rsi, 64, clmblkdst8, 2*clmblkdst8, rdi, rdx); \
		rdi += XMM_SCND2;		/* Next sine/cosine pointer */ \
		rdx += XMM_SCND2R;		/* Next sine/cosine pointer */ \
	}while(--loopA & (2*clm-1)); \
	newxloop1(clm, b6b, clmblkdst);	/* 2*clm iterations */ \
} \
rsi += -8*clmblkdst+4*clmblkdst8; \
rdi += -8*2*clm*XMM_SCND2;	/* Restore sin/cos ptr */ \
\
loops_reset; \
r4_x4cl_wpn_four_complex_djbfft_preload; \
for(unsigned int loopA = 7; loopA; loopA--) { \
	for(unsigned int loopB = 8*2*clm; loopB; ) { \
		xtouch(xptr(rcx+rbp+(clm-1)*128));	/* Preload TLB */ \
		do{ \
			r4_x4cl_wpn_four_complex_djbfft(rsi, 64, clmblkdst8, 2*clmblkdst8, rdi); \
			rdi += XMM_SCND2;		/* Next sine/cosine pointer */ \
		}while(--loopB & (2*clm-1)); \
		newxloop1(clm, b7b, clmblkdst);	/* 2*clm iterations */ \
	} \
	rsi += -8*clmblkdst+4*clmblkdst8; /* Next source pointer */ \
	rdi += -8*2*clm*XMM_SCND2;	/* Restore sin/cos ptr */ \
} \
rsi += -8*4*clmblkdst8;	/* Restore source pointer */ \
end_timer(25); \
\
/* Do FFT levels 8,9,10 */ \
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
/* Do 62 eight_complex_fft8 operations */ \
/*	distance between fft data elements is 1K */ \
/*	do 31 macros each processing 64 data values */ \
\
start_timer(26); \
rdi += -8*2*clm*XMM_SCND2R-16*2*clm*XMM_SCD8+128; /* Complex data sin/cos ptr and use -128 */ \
;					/* addressing to generate 1 byte offsets in radix-8 macro */ \
rdx = (rdi-128-2*clm*XMM_SCD4); /* Real data sin/cos ptr */ \
rcx = (uintptr_t)g->data_addr;/* Destination pointer */ \
loops_init_prefetch(32*2*clm, 128, 1, rbp); /* 32 iters of 2*clm */ \
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
for(unsigned int loopA = 16*2-1; loopA; ) { \
	do{ \
		for(unsigned int loopC = 2*clm; loopC; loopC--) { \
			r8_sg8cl_eight_complex_fft8(rsi, 64, clmblkdst, 2*clmblkdst, 4*clmblkdst, rcx, 64, blkdst, 2*blkdst, 4*blkdst, rdi-128); \
			rdi += XMM_SCD8;		/* Next premult/sin/cos ptr */ \
		} \
		rsi += -2*clm*64+clmblkdst8; /* Next source pointer */ \
		rcx += -2*clm*64+8*blkdst;	/* Next dest pointer */ \
		rdi += -2*clm*XMM_SCD8;	/* Restore premult/sine/cosine pointer */ \
	}while(--loopA & 1); \
	rdi += 2*clm*XMM_SCD8;	/* Next premult/sine/cosine pointer */ \
} \
end_timer(26); \
\
/* Work on next 4*clm sets of 1024 values */ \
\
endlp:	xpass1_get_next_block(pass2, c0b, b0b, _wpn); /* Get next block and jump */ \
/* done: */ xfft_footer_rn; \
\
\
\
/* Perform a pass 1 FFT of length 1280 doubles. */ \
/* Memory examples are for a 1280K FFT with 10 levels done in pass 2. */ \
/* Sin/cos data needed: 16*2*clm*XMM_SCD2 + 16*2*clm*XMM_SCD2 + 4*2*clm*XMM_SCND2 + 4*2*clm*XMM_SCND2R + */ \
/*			2*clm*XMM_SCD2 + 10*2*clm*XMM_SCD4 = 67*clm 64-byte cache lines */ \
\

#define r4dwpn_pass1sc1280(procname, xpass2name, pass2_levels, clm) /* (fftname, pass2_levels, clmarg) */ \
\
/* set_FFT_constants(fftname, 0, 0, pass2_levels, clmarg); */ \
xfft_header_rn(procname, xpass2name, pass2_levels, clm); \
\
/* Do pass 2 (the last 10 FFT levels if this is a 1280K FFT) */ \
\
/* pass2: */start_timer(0); \
xpass2name(g); \
end_timer(0); \
if(g->ffttype == 1)		/* We're done if FFTing only */ \
	return; /* goto done; */ \
\
/* Do the inverse FFT */ \
\
pass1_inverse_fft_setup(c0b, procname); \
\
/* Do inverse FFT levels 10,11 */ \
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
/* Do 159 four_complex_unfft4 operations */ \
/*	distance between fft data elements is 1K */ \
/*	do 159 macros each processing 16 data values */ \
\
c0b:	start_timer(16); \
rsi = (uintptr_t)g->data_addr;/* Get source address */ \
rcx = (uintptr_t)g->scratch_area;	/* Get address of scratch area */ \
rdx = (uintptr_t)g->premult_addr;/* Real data sin/cos pointer */ \
rdi = (rdx+2*clm*XMM_SCD2); /* Complex data sin/cos ptr */ \
pfing(rbp = (rdx+(67-40)*2*clm*64)); /* Finish prefetching the sin/cos data */ \
loops_init_prefetch(160*2*clm, 64, 4, rbp); /* 160*2*clm iterations */ \
/*bug	xtouchpages rbp, 160*2*clm * 64 / 4 ;; Prefetch TLBs */ \
for(unsigned int loopA = 2*clm; loopA; loopA--) { \
	r4_sg2cl_eight_reals_unfft4(rsi, 64, blkdst, rcx, 64, clmblkdst, rdi, rdx); \
	rdi += XMM_SCD4;		/* Next premult/sin/cos data */ \
	rdx += XMM_SCD2;		/* Next sine/cosine pointer */ \
} \
rsi += -2*clm*64+2*blkdst;	/* Next source pointer */ \
rcx += -2*clm*64+2*clmblkdst; /* Next dest pointer */ \
rdi += -2*clm*XMM_SCD4;	/* Restore premult/sine/cosine pointer */ \
\
loops_reset; \
for(unsigned int loopA = 10*4*4-1; loopA; ) { \
	do{ \
		do{ \
			for(unsigned int loopD = 2*clm; loopD; loopD--) { \
				r4_sg2cl_four_complex_unfft4(rsi, 64, blkdst, rcx, 64, clmblkdst, rdi); \
				rdi += XMM_SCD4;		/* Next premult/sin/cos data */ \
			} \
			rsi += -2*clm*64+2*blkdst;	/* Next source pointer */ \
			rcx += -2*clm*64+2*clmblkdst; /* Next dest pointer */ \
			rdi += -2*clm*XMM_SCD4;	/* Restore premult/sin/cos data pointer */ \
		}while(--loopA & 3); \
		rcx += -4*2*clmblkdst+clmblkdst8; \
	}while(loopA & (4*4-1)); \
	rdi += 2*clm*XMM_SCD4;	/* Next premult/sin/cos pointer */ \
} \
rsi = (uintptr_t)g->scratch_area;	/* Restore source pointer */ \
rdx = rdi;		/* Next real data sin/cos pointer */ \
end_timer(16); \
\
/* Do inverse FFT levels 8,9 */ \
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
/* Do 156 four_complex_unfft operations */ \
/*	distance between fft data elements is 4K */ \
/*	do 78 macros each processing 32 data values */ \
\
start_timer(17); \
rdi = (rdx+4*2*clm*XMM_SCND2R); /* Complex data sin/cos ptr */ \
loops_init(80*2*clm);		/* 80*2*clm iterations */ \
\
r4_x4cl_wpn_eight_reals_unfft_preload; \
for(unsigned int loopA = 2; loopA; loopA--) { \
	for(unsigned int loopB = 2*clm; loopB; loopB--) { \
		r4_x4cl_wpn_eight_reals_unfft(rsi, 64, 2*clmblkdst, 4*clmblkdst, rdi, 2*2*clm*XMM_SCND2, rdx, 2*2*clm*XMM_SCND2R); \
		rdi += XMM_SCND2;		/* Next sine/cosine pointer */ \
		rdx += XMM_SCND2R;		/* Next sine/cosine pointer */ \
	} \
	rsi += -2*clm*64+clmblkdst; /* Next source pointer */ \
} \
rsi += -2*clmblkdst+clmblkdst8; /* Next source pointer */ \
rdi += -2*2*clm*XMM_SCND2;	/* Restore sine/cosine pointer */ \
\
loops_reset; \
r4_x4cl_wpn_four_complex_djbunfft_preload; \
for(unsigned int loopA = 39; loopA; loopA--) { \
	for(unsigned int loopB = 2; loopB; loopB--) { \
		for(unsigned int loopC = 2*clm; loopC; loopC--) { \
			r4_x4cl_wpn_four_complex_djbunfft(rsi, 64, 2*clmblkdst, 4*clmblkdst, rdi, 2*2*clm*XMM_SCND2); \
			rdi += XMM_SCND2;		/* Next clm's sin/cos data */ \
		} \
		rsi += -2*clm*64+clmblkdst; /* Next source pointer */ \
	} \
	rsi += -2*clmblkdst+clmblkdst8; /* Next source pointer */ \
	rdi += -2*2*clm*XMM_SCND2;	/* Restore sine/cosine pointer */ \
} \
rsi += -40*clmblkdst8;	/* Restore source pointer */ \
rdx = (rdi+4*2*clm*XMM_SCND2); /* Next real data sin/cos ptr */ \
end_timer(17); \
\
/* Do inverse FFT levels 6,7 */ \
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
/*	128K	... */ \
/*	... */ \
\
/* Do 16 eight_reals_unfft operations */ \
/*	distance between fft data elements is 16K */ \
/*	do 8 macros each processing 32 data values */ \
/* Do 144 four_complex_unfft operations */ \
/*	distance between fft data elements is 16K */ \
/*	do 72 macros each processing 32 data values */ \
\
start_timer(18); \
rdi = (rdx+16*2*clm*XMM_SCD2); /* Complex data sin/cos ptr */ \
loops_init(80*2*clm);		/* 80*2*clm iterations */ \
\
r4_x4cl_eight_reals_unfft_preload; \
for(unsigned int loopA = 8; loopA; loopA--) { \
	for(unsigned int loopB = 2*clm; loopB; loopB--) { \
		r4_x4cl_eight_reals_unfft(rsi, 64, clmblkdst8, 2*clmblkdst8, rdi, 8*2*clm*XMM_SCD2, rdx, 8*2*clm*XMM_SCD2); \
		rdi += XMM_SCD2;		/* Next sine/cosine pointer */ \
		rdx += XMM_SCD2;		/* Next sine/cosine pointer */ \
	} \
	rsi += -2*clm*64+clmblkdst; /* Next source pointer */ \
} \
rsi += -8*clmblkdst+4*clmblkdst8; /* Next source pointer */ \
rdi += -8*2*clm*XMM_SCD2;	/* Restore sine/cosine pointer */ \
\
loops_reset; \
r4_x4cl_four_complex_djbunfft_preload; \
for(unsigned int loopA = 9; loopA; loopA--) { \
	for(unsigned int loopB = 8; loopB; loopB--) { \
		for(unsigned int loopC = 2*clm; loopC; loopC--) { \
			r4_x4cl_four_complex_djbunfft(rsi, 64, clmblkdst8, 2*clmblkdst8, rdi, 8*2*clm*XMM_SCD2); \
			rdi += XMM_SCD2;		/* Next clm's sin/cos data */ \
		} \
		rsi += -2*clm*64+clmblkdst; /* Next source pointer */ \
	} \
	rsi += -8*clmblkdst+4*clmblkdst8; /* Next source pointer */ \
	rdi += -8*2*clm*XMM_SCD2;	/* Restore sine/cosine pointer */ \
} \
rsi += -10*4*clmblkdst8;	/* Restore source pointer */ \
end_timer(18); \
\
/* Do inverse FFT levels 1,2,3,4,5 (actually 4.32 levels) */ \
/* On input the 64-byte cache lines hold these data values: */ \
/*	0K	32K	64K	96K */ \
/*	1K	... */ \
/*	... */ \
/*	31K */ \
/*	128K	... */ \
/*	... */ \
/* On output the 64-byte cache lines hold these data values: */ \
/*	0K	320K	640K	960K */ \
/*	1K	... */ \
/*	... */ \
/*	319K */ \
\
/* Do 256 20_reals_last_unfft operations */ \
/*	distance between fft data elements is 64K */ \
/*	do 32 macros each processing 80 data values */ \
\
start_timer(20); \
rdi = (uintptr_t)g->sincos2;/* Load fixed sin/cos postmultipliers */ \
set_data_prefetch_ptrs;		/* Init prefetching pointers */ \
loops_init(32*2*clm);		/* 32 iterations of 2*clm */ \
r5_x10cl_20_reals_last_unfft_preload; \
for(unsigned int loopA = 4; loopA; loopA--) { \
	for(unsigned int loopB = 8*2*clm; loopB; ) { \
		xtouch(xptr(rcx+rbp+(clm-1)*128));	/* Preload TLB */ \
		do{ \
			r5_x10cl_20_reals_last_unfft(rsi, 64, 4*clmblkdst8, rdi, 32*XMM_SCD9); \
		}while(--loopB & (2*clm-1)); \
		newxloop1(clm, c1b, clmblkdst, rdi, XMM_SCD9); /* 2*clm iterations */ \
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
copy_scratch_data_to_fft(320, clm); \
goto	endlp;			/* Skip FFT code if flag not set */ \
\
/* Fire up auxiliary threads to do pass 1 forward FFT */ \
\
/* pass1: */pass1_forward_fft_setup(b0b); \
\
/* Copy data from the FFT data area back to the scratch area */ \
\
b0b:	copy_fft_data_to_scratch(320, clm); \
\
/* Do FFT levels 1,2,3,4,5 (actually 4.32 levels) */ \
/* */ \
/* On input the 64-byte cache lines hold these data values: */ \
/*	0K	320K	640K	960K */ \
/*	1K	... */ \
/*	... */ \
/*	319K */ \
/* On output the 64-byte cache lines hold these data values: */ \
/*	0K	64K	128K	196K */ \
/*	1K	... */ \
/*	... */ \
/*	63K	... */ \
/*	256K	... */ \
/*	... */ \
\
/* Do 256 20_reals_first_fft operations */ \
/*	distance between fft data elements is 64K */ \
/*	do 64 iterations each processing 40 data values */ \
\
yes_fft:start_timer(22); \
rsi = (uintptr_t)g->scratch_area;	/* Get address of scratch area */ \
rdi = (uintptr_t)g->sincos2;/* Load fixed sin/cos postmultipliers */ \
set_data_prefetch_ptrs;		/* Init prefetching pointers */ \
pfing(rcx += (32-0)*blkdst);	/* We've already prefetched 32 blocks, perfect! */ \
loops_init(64*2*clm);		/* 64 iters of 2*clm */ \
r5_x5cl_20_reals_first_fft_scratch_preload; \
for(unsigned int loopA = 8; loopA; loopA--) { \
	for(unsigned int loopB = 8*2*clm; loopB; ) { \
		xtouch(xptr(rcx+rbp+(clm-1)*128));	/* Preload TLB */ \
		do{ \
			r5_x5cl_20_reals_first_fft_scratch(rsi, 64, 8*clmblkdst8, rdi); \
		}while(--loopB & (2*clm-1)); \
		newxloop1doublepre(clm, b1b, b1a, clmblkdst, rdi, XMM_SCD9); /* 2*clm iters prefetching 4*clm cache lines */ \
	} \
	rsi += -8*clmblkdst+clmblkdst8; \
} \
rsi += -8*clmblkdst8;	/* Restore source pointer */ \
end_timer(22); \
\
/* Do FFT levels 6,7 */ \
/* */ \
/* On input the 64-byte cache lines hold these data values: */ \
/*	0K	64K	128K	196K */ \
/*	1K	... */ \
/*	... */ \
/*	63K	... */ \
/*	256K	... */ \
/*	... */ \
/* On output the 64-byte cache lines hold these data values: */ \
/*	0K	16K	32K	48K */ \
/*	1K	... */ \
/*	... */ \
/*	15K	... */ \
/*	64K	... */ \
/*	... */ \
\
/* Do 16 eight_reals_fft operations */ \
/* Do 16 four_complex_fft operations */ \
/*	distance between fft data elements is 16K */ \
/*	do 16 macros each processing 32 data values */ \
/* Do 128 four_complex_fft operations */ \
/*	distance between fft data elements is 16K */ \
/*	do 64 macros each processing 32 data values */ \
\
start_timer(24); \
rdx = (uintptr_t)g->premult_addr;/* Address of premultiplier/sin/cos table */ \
rdx += 4*2*clm*XMM_SCND2 + 4*2*clm*XMM_SCND2R + 2*clm*XMM_SCD2 + 10*2*clm*XMM_SCD4; /* Real data sin/cos ptr */ \
rdi = (rdx+16*2*clm*XMM_SCD2); /* Next sin/cos ptr */ \
loops_init(80*2*clm);		/* 5 iterations of 16 iters of 2*clm */ \
r4_x4cl_eight_reals_four_complex_djbfft_preload; \
for(unsigned int loopA = 2; loopA; loopA--) { \
	for(unsigned int loopB = 8*2*clm; loopB; ) { \
		xtouch(xptr(rcx+rbp+(clm-1)*128));	/* Preload TLB */ \
		do{ \
			r4_x4cl_eight_reals_four_complex_djbfft(rsi, 64, 2*clmblkdst8, 4*clmblkdst8, rdi, rdx); \
			rdi += XMM_SCD2;		/* Next sine/cosine pointer */ \
			rdx += XMM_SCD2;		/* Next sine/cosine pointer */ \
		}while(--loopB & (2*clm-1)); \
		newxloop1(clm, b4b, clmblkdst);	/* 2*clm iterations */ \
	} \
	rsi += -8*clmblkdst+clmblkdst8; \
} \
rsi += -2*clmblkdst8+8*clmblkdst8; \
rdi += -16*2*clm*XMM_SCD2;	/* Restore sin/cos ptr */ \
\
loops_reset; \
r4_x4cl_four_complex_djbfft_preload; \
for(unsigned int loopA = 4; loopA; loopA--) { \
	for(unsigned int loopB = 2; loopB; loopB--) { \
		for(unsigned int loopC = 8*2*clm; loopC; ) { \
			xtouch(xptr(rcx+rbp+(clm-1)*128));	/* Preload TLB */ \
			do{ \
				r4_x4cl_four_complex_djbfft(rsi, 64, 2*clmblkdst8, 4*clmblkdst8, rdi); \
				rdi += XMM_SCD2;		/* Next sine/cosine pointer */ \
			}while(--loopC & (2*clm-1)); \
			newxloop1(clm, b5b, clmblkdst);	/* 2*clm iterations */ \
		} \
		rsi += -8*clmblkdst+clmblkdst8; /* Next source pointer */ \
	} \
	rsi += -2*clmblkdst8+8*clmblkdst8; /* Next source pointer */ \
	rdi += -16*2*clm*XMM_SCD2;	/* Restore sin/cos ptr */ \
} \
rsi += -5*8*clmblkdst8;	/* Restore source pointer */ \
end_timer(24); \
\
/* Do FFT levels 8,9 */ \
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
/* Do 152 four_complex_fft operations */ \
/*	distance between fft data elements is 4K */ \
/*	do 76 macros each processing 32 data values */ \
\
start_timer(25); \
rdi += -16*2*clm*XMM_SCD2-4*2*clm*XMM_SCND2; /* Complex data sin/cos ptr */ \
rdx = (rdi-4*2*clm*XMM_SCND2R); /* Real data sin/cos ptr */ \
loops_init(80*2*clm);		/* 20 iterations of 4 iters of 2*clm */ \
r4_x4cl_wpn_eight_reals_four_complex_djbfft_preload; \
for(unsigned int loopA = 4*2*clm; loopA; ) { \
	xtouch(xptr(rcx+rbp+(clm-1)*128));	/* Preload TLB */ \
	do{ \
		r4_x4cl_wpn_eight_reals_four_complex_djbfft(rsi, 64, 4*clmblkdst, clmblkdst8, rdi, rdx); \
		rdi += XMM_SCND2;		/* Next sine/cosine pointer */ \
		rdx += XMM_SCND2R;		/* Next sine/cosine pointer */ \
	}while(--loopA & (2*clm-1)); \
	newxloop1(clm, b6b, clmblkdst);	/* 2*clm iterations */ \
} \
rsi += -4*clmblkdst+2*clmblkdst8; \
rdi += -4*2*clm*XMM_SCND2;	/* Restore sin/cos ptr */ \
\
loops_reset; \
r4_x4cl_wpn_four_complex_djbfft_preload; \
for(unsigned int loopA = 19; loopA; loopA--) { \
	for(unsigned int loopB = 4*2*clm; loopB; ) { \
		xtouch(xptr(rcx+rbp+(clm-1)*128));	/* Preload TLB */ \
		do{ \
			r4_x4cl_wpn_four_complex_djbfft(rsi, 64, 4*clmblkdst, clmblkdst8, rdi); \
			rdi += XMM_SCND2;		/* Next sine/cosine pointer */ \
		}while(--loopB & (2*clm-1)); \
		newxloop1(clm, b7b, clmblkdst);	/* 2*clm iterations */ \
	} \
	rsi += -4*clmblkdst+2*clmblkdst8; /* Next source pointer */ \
	rdi += -4*2*clm*XMM_SCND2;	/* Restore sin/cos ptr */ \
} \
rsi += -20*2*clmblkdst8;	/* Restore source pointer */ \
end_timer(25); \
\
/* Do FFT levels 10,11 */ \
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
/* Do 1 four_complex_fft4 operations */ \
/*	distance between fft data elements is 1K */ \
/*	do 1 macros each processing 32 data values */ \
/* Do 158 four_complex_fft4 operations */ \
/*	distance between fft data elements is 1K */ \
/*	do 79 macros each processing 32 data values */ \
\
start_timer(26); \
rdi += -4*2*clm*XMM_SCND2R-10*2*clm*XMM_SCD4; /* Complex data sin/cos ptr */ \
rdx = (rdi-2*clm*XMM_SCD2); /* Real data sin/cos ptr */ \
rcx = (uintptr_t)g->data_addr;/* Destination pointer */ \
loops_init_prefetch(80*2*clm, 64, 2, rbp); /* 80 iters of 2*clm */ \
pfing(rbp = (uintptr_t)g->premult_prefetch);/* Load premult/sin/cos prefetch pointer */ \
r4_sg4cl_eight_reals_four_complex_fft4_preload; \
for(unsigned int loopA = 2*clm; loopA; loopA--) { \
	r4_sg4cl_eight_reals_four_complex_fft4(rsi, 64, clmblkdst, 2*clmblkdst, rcx, 64, blkdst, 2*blkdst, rdi, rdx); \
	rdi += XMM_SCD4;		/* Next sine/cosine pointer */ \
	rdx += XMM_SCD2;		/* Next sine/cosine pointer */ \
} \
rsi += -2*clm*64+4*clmblkdst; /* Next source pointer */ \
rcx += -2*clm*64+4*blkdst;	/* Next dest pointer */ \
rdi += -2*clm*XMM_SCD4;	/* Restore premult/sine/cosine pointer */ \
\
loops_reset; \
for(unsigned int loopA = 10*4*2-1; loopA; ) { \
	do{ \
		do{ \
			for(unsigned int loopD = 2*clm; loopD; loopD--) { \
				r4_sg4cl_four_complex_fft4(rsi, 64, clmblkdst, 2*clmblkdst, rcx, 64, blkdst, 2*blkdst, rdi); \
				rdi += XMM_SCD4;		/* Next premult/sin/cos ptr */ \
			} \
			rsi += -2*clm*64+4*clmblkdst; /* Next source pointer */ \
			rcx += -2*clm*64+4*blkdst;	/* Next dest pointer */ \
			rdi += -2*clm*XMM_SCD4;	/* Restore premult/sine/cosine pointer */ \
		}while(--loopA & 1); \
		rsi += -2*4*clmblkdst+clmblkdst8; /* Next source pointer */ \
	}while(loopA & (4*2-1)); \
	rdi += 2*clm*XMM_SCD4;	/* Next premult/sine/cosine pointer */ \
} \
end_timer(26); \
\
/* Work on next 4*clm sets of 1280 values */ \
\
endlp:	xpass1_get_next_block(pass2, c0b, b0b, _wpn); /* Get next block and jump */ \
/* done: */ xfft_footer_rn; \
\
\
\
/* Perform a pass 1 FFT of length 1536 doubles. */ \
/* Memory examples are for a 1536K FFT with 10 levels done in pass 2. */ \
/* Sin/cos data needed: 16*2*clm*XMM_SCD1 + 16*2*clm*XMM_SCD1 + 4*2*clm*XMM_SCND2 + 4*2*clm*XMM_SCD2R + */ \
/*			2*clm*XMM_SCD2 + 16*2*clm*XMM_SCD4 = 63*2*clm 64-byte cache lines */ \
\

#define r4dwpn_pass1sc1536(procname, xpass2name, pass2_levels, clm) /* (fftname, pass2_levels, clmarg) */ \
\
/* set_FFT_constants(fftname, 0, 0, pass2_levels, clmarg); */ \
xfft_header_rn(procname, xpass2name, pass2_levels, clm); \
\
/* Do pass 2 (the last 10 FFT levels if this is a 1536K FFT) */ \
\
/* pass2: */start_timer(0); \
xpass2name(g); \
end_timer(0); \
if(g->ffttype == 1)		/* We're done if FFTing only */ \
	return; /* goto done; */ \
\
/* Do the inverse FFT */ \
\
pass1_inverse_fft_setup(c0b, procname); \
\
/* Do inverse FFT levels 10,11 */ \
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
/* Do 191 four_complex_unfft4 operations */ \
/*	distance between fft data elements is 1K */ \
/*	do 191 macros each processing 16 data values */ \
\
c0b:	start_timer(16); \
rsi = (uintptr_t)g->data_addr;/* Get source address */ \
rcx = (uintptr_t)g->scratch_area;	/* Get address of scratch area */ \
rdx = (uintptr_t)g->premult_addr;/* Real data sin/cos pointer */ \
rdi = (rdx+2*clm*XMM_SCD2); /* Complex data sin/cos ptr */ \
pfing(rbp = (rdx+(63-24)*2*clm*64)); /* Finish prefetching the sin/cos data */ \
loops_init_prefetch(192*2*clm, 64, 8, rbp); /* 192*2*clm iterations */ \
/*bug	xtouchpages rbp, 192*2*clm * 64 / 8 ;; Prefetch TLBs */ \
for(unsigned int loopA = 2*clm; loopA; loopA--) { \
	r4_sg2cl_eight_reals_unfft4(rsi, 64, blkdst, rcx, 64, clmblkdst, rdi, rdx); \
	rdi += XMM_SCD4;		/* Next premult/sin/cos data */ \
	rdx += XMM_SCD2;		/* Next sine/cosine pointer */ \
} \
rsi += -2*clm*64+2*blkdst;	/* Next source pointer */ \
rcx += -2*clm*64+2*clmblkdst; /* Next dest pointer */ \
rdi += -2*clm*XMM_SCD4;	/* Restore premult/sine/cosine pointer */ \
\
loops_reset; \
for(unsigned int loopA = 16*3*4-1; loopA; ) { \
	for(unsigned int loopB = 3; loopB; loopB--) { \
		do{ \
			for(unsigned int loopD = 2*clm; loopD; loopD--) { \
				r4_sg2cl_four_complex_unfft4(rsi, 64, blkdst, rcx, 64, clmblkdst, rdi); \
				rdi += XMM_SCD4;		/* Next premult/sin/cos data */ \
			} \
			rsi += -2*clm*64+2*blkdst;	/* Next source pointer */ \
			rcx += -2*clm*64+2*clmblkdst; /* Next dest pointer */ \
			rdi += -2*clm*XMM_SCD4;	/* Restore premult/sin/cos data pointer */ \
		}while(--loopA & 3); \
		rcx += -4*2*clmblkdst+clmblkdst8; \
	} \
	rdi += 2*clm*XMM_SCD4;	/* Next premult/sin/cos pointer */ \
} \
rsi = (uintptr_t)g->scratch_area;	/* Restore source pointer */ \
rdx = rdi;		/* Next real data sin/cos pointer */ \
end_timer(16); \
\
/* Do inverse FFT levels 8,9 */ \
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
/* Do 188 four_complex_unfft operations */ \
/*	distance between fft data elements is 4K */ \
/*	do 94 macros each processing 32 data values */ \
\
start_timer(17); \
rdi = (rdx+4*2*clm*XMM_SCND2R); /* Real data sin/cos ptr */ \
loops_init(96*2*clm);		/* 96*2*clm iterations */ \
\
r4_x4cl_wpn_eight_reals_unfft_preload; \
for(unsigned int loopA = 2; loopA; loopA--) { \
	for(unsigned int loopB = 2*clm; loopB; loopB--) { \
		r4_x4cl_wpn_eight_reals_unfft(rsi, 64, 2*clmblkdst, 4*clmblkdst, rdi, 2*2*clm*XMM_SCND2, rdx, 2*2*clm*XMM_SCND2R); \
		rdi += XMM_SCND2;		/* Next sine/cosine pointer */ \
		rdx += XMM_SCND2R;		/* Next sine/cosine pointer */ \
	} \
	rsi += -2*clm*64+clmblkdst; /* Next source pointer */ \
} \
rsi += -2*clmblkdst+clmblkdst8; /* Next source pointer */ \
rdi += -2*2*clm*XMM_SCND2;	/* Restore sine/cosine pointer */ \
\
loops_reset; \
r4_x4cl_wpn_four_complex_djbunfft_preload; \
for(unsigned int loopA = 47; loopA; loopA--) { \
	for(unsigned int loopB = 2; loopB; loopB--) { \
		for(unsigned int loopC = 2*clm; loopC; loopC--) { \
			r4_x4cl_wpn_four_complex_djbunfft(rsi, 64, 2*clmblkdst, 4*clmblkdst, rdi, 2*2*clm*XMM_SCND2); \
			rdi += XMM_SCND2;		/* Next clm's sin/cos data */ \
		} \
		rsi += -2*clm*64+clmblkdst; /* Next source pointer */ \
	} \
	rsi += -2*clmblkdst+clmblkdst8; /* Next source pointer */ \
	rdi += -2*2*clm*XMM_SCND2;	/* Restore sine/cosine pointer */ \
} \
rsi += -48*clmblkdst8;	/* Restore source pointer */ \
rdx = (rdi+4*2*clm*XMM_SCND2); /* Next real data sin/cos ptr */ \
end_timer(17); \
\
/* Do inverse FFT levels 6,7 (actually 1.585 levels) */ \
/* On input the 64-byte cache lines hold these data values: */ \
/*	0K	8K	16K	24K */ \
/*	1K	... */ \
/*	... */ \
/*	7K */ \
/*	32K	... */ \
/*	... */ \
/* On output the 64-byte cache lines hold these data values: */ \
/*	0K	24K	48K	72K */ \
/*	1K	... */ \
/*	... */ \
/*	23K */ \
/*	96K	... */ \
/*	... */ \
\
/* Do 16 six_reals_unfft operations */ \
/*	distance between fft data elements is 16K */ \
/*	do 8 macros each processing 24 data values */ \
/* Do 240 three_complex_unfft operations */ \
/*	distance between fft data elements is 16K */ \
/*	do 120 macros each processing 24 data values */ \
\
start_timer(18); \
rdi = (rdx+16*2*clm*XMM_SCD1); /* Real data sin/cos ptr */ \
loops_init(128*2*clm);		/* 128 iterations of 2*clm */ \
for(unsigned int loopA = 8; loopA; loopA--) { \
	for(unsigned int loopB = 2*clm; loopB; loopB--) { \
		r3_x3cl_six_reals_unfft(rsi, 64, clmblkdst8, rdi, 8*2*clm*XMM_SCD1, rdx, 8*2*clm*XMM_SCD1); \
		rdi += XMM_SCD1;		/* Next sine/cosine pointer */ \
		rdx += XMM_SCD1;		/* Next sine/cosine pointer */ \
	} \
	rsi += -2*clm*64+clmblkdst; /* Next source pointer */ \
} \
rsi += -8*clmblkdst+3*clmblkdst8; /* Next source pointer */ \
rdi += -8*2*clm*XMM_SCD1;	/* Restore sine/cosine pointer */ \
\
loops_reset; \
r3_x3cl_three_complex_djbunfft_preload; \
for(unsigned int loopA = 15; loopA; loopA--) { \
	for(unsigned int loopB = 8; loopB; loopB--) { \
		for(unsigned int loopC = 2*clm; loopC; loopC--) { \
			r3_x3cl_three_complex_djbunfft(rsi, 64, clmblkdst8, rdi, 8*2*clm*XMM_SCD1); \
			rdi += XMM_SCD1;		/* Next clm's sin/cos data */ \
		} \
		rsi += -2*clm*64+clmblkdst; /* Next source pointer */ \
	} \
	rsi += -8*clmblkdst+3*clmblkdst8; /* Next source pointer */ \
	rdi += -8*2*clm*XMM_SCD1;	/* Restore sine/cosine pointer */ \
} \
rsi += -16*3*clmblkdst8;	/* Restore source pointer */ \
end_timer(18); \
\
/* Do inverse FFT levels 4,5 */ \
/* On input the 64-byte cache lines hold these data values: */ \
/*	0K	24K	48K	72K */ \
/*	1K	... */ \
/*	... */ \
/*	23K */ \
/*	96K	... */ \
/*	... */ \
/* On output the 64-byte cache lines hold these data values: */ \
/*	0K	96K	192K	288K */ \
/*	1K	... */ \
/*	... */ \
/*	95K */ \
/*	384K	... */ \
/*	... */ \
\
/* Do 48 eight_reals_unfft operations */ \
/*	distance between fft data elements is 48K */ \
/*	do 24 macros each processing 32 data values */ \
/* Do 144 four_complex_unfft operations */ \
/*	distance between fft data elements is 48K */ \
/*	do 72 macros each processing 32 data values */ \
\
start_timer(19); \
rdx = (uintptr_t)g->sincos2;/* Load fixed sin/cos data ptr */ \
rdx += 192*XMM_SCD3;	/* Real data sin/cos ptr */ \
rdi = (rdx+48*XMM_SCD2);	/* Complex data sin/cos ptr */ \
\
loops_init(96*2*clm);		/* 24*2*clm iterations */ \
r4_x4cl_eight_reals_unfft_preload; \
for(unsigned int loopA = 3; loopA; loopA--) { \
	for(unsigned int loopB = 8; loopB; loopB--) { \
		for(unsigned int loopC = 2*clm; loopC; loopC--) { \
			r4_x4cl_eight_reals_unfft(rsi, 64, 3*clmblkdst8, 6*clmblkdst8, rdi, 24*XMM_SCD2, rdx, 24*XMM_SCD2); \
		} \
		rsi += -2*clm*64+clmblkdst; /* Next source pointer */ \
		rdi += XMM_SCD2;		/* Next sine/cosine pointer */ \
		rdx += XMM_SCD2;		/* Next sine/cosine pointer */ \
	} \
	rsi += -8*clmblkdst+clmblkdst8; /* Next source pointer */ \
} \
rsi += -3*clmblkdst8+12*clmblkdst8; /* Next source pointer */ \
rdi += -24*XMM_SCD2;	/* Restore sine/cosine pointer */ \
\
loops_reset; \
r4_x4cl_four_complex_djbunfft_preload; \
for(unsigned int loopA = 3; loopA; loopA--) { \
	for(unsigned int loopB = 3; loopB; loopB--) { \
		for(unsigned int loopC = 8; loopC; loopC--) { \
			for(unsigned int loopD = 2*clm; loopD; loopD--) { \
				r4_x4cl_four_complex_djbunfft(rsi, 64, 3*clmblkdst8, 6*clmblkdst8, rdi, 24*XMM_SCD2); \
			} \
			rsi += -2*clm*64+clmblkdst; /* Next source pointer */ \
			rdi += XMM_SCD2;		/* Next clm's sin/cos data */ \
		} \
		rsi += -8*clmblkdst+clmblkdst8; /* Next source pointer */ \
	} \
	rsi += -3*clmblkdst8+12*clmblkdst8; /* Next source pointer */ \
	rdi += -24*XMM_SCD2;	/* Restore sine/cosine pointer */ \
} \
rsi += -4*12*clmblkdst8;	/* Restore source pointer */ \
end_timer(19); \
\
/* Do inverse FFT levels 1,2,3 */ \
/* On input the 64-byte cache lines hold these data values: */ \
/*	0K	96K	192K	288K */ \
/*	1K	... */ \
/*	... */ \
/*	95K */ \
/*	384K	... */ \
/*	... */ \
/* On output the 64-byte cache lines hold these data values: */ \
/*	0K	384K	768K	1152K */ \
/*	1K	... */ \
/*	... */ \
/*	383K */ \
\
/* Do 192 eight_reals_last_unfft operations */ \
/*	distance between fft data elements is 192K */ \
/*	do 96 macros each processing 32 data values */ \
\
start_timer(20); \
rdi = (uintptr_t)g->sincos2;/* Load fixed sin/cos postmultipliers */ \
loops_init(96*2*clm);		/* 96 iterations of 2*clm */ \
r4_x4cl_eight_reals_last_unfft_preload; \
for(unsigned int loopA = 12; loopA; loopA--) { \
	for(unsigned int loopB = 8; loopB; loopB--) { \
		for(unsigned int loopC = 2*clm; loopC; loopC--) { \
			r4_x4cl_eight_reals_last_unfft(rsi, 64, 12*clmblkdst8, 24*clmblkdst8, rdi, 96*XMM_SCD3); \
		} \
		newxloop1nopre(clm, c1b, clmblkdst, rdi, XMM_SCD3); /* 2*clm iterations */ \
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
copy_scratch_data_to_fft(384, clm); \
goto	endlp;			/* Skip FFT code if flag not set */ \
\
/* Fire up auxiliary threads to do pass 1 forward FFT */ \
\
/* pass1: */pass1_forward_fft_setup(b0b); \
\
/* Copy data from the FFT data area back to the scratch area */ \
\
b0b:	copy_fft_data_to_scratch(384, clm); \
\
/* Do FFT levels 1,2,3 */ \
/* */ \
/* On input the 64-byte cache lines hold these data values: */ \
/*	0K	384K	768K	1152K */ \
/*	1K	... */ \
/*	... */ \
/*	383K */ \
/* On output the 64-byte cache lines hold these data values: */ \
/*	0K	192K	384K	576K */ \
/*	1K	... */ \
/*	... */ \
/*	191K	... */ \
/*	768K	... */ \
/*	... */ \
\
/* Do 192 eight_reals_first_fft operations */ \
/*	distance between fft data elements is 192K */ \
/*	do 192 iterations each processing 16 data values */ \
\
yes_fft:start_timer(22); \
rsi = (uintptr_t)g->scratch_area;	/* Get address of scratch area */ \
rdi = (uintptr_t)g->sincos2;/* Load fixed sin/cos postmultipliers */ \
set_data_prefetch_ptrs;		/* Init prefetching pointers */ \
loops_init(192*2*clm);		/* 24 iters of 4 iters of 2*2*clm */ \
r4_x2cl_eight_reals_first_fft_scratch_preload; \
for(unsigned int loopA = 24; loopA; loopA--) { \
	for(unsigned int loopB = 4*2*2*clm; loopB; ) { \
		xtouch(xptr(rcx+rbp+(clm-1)*128));	/* Preload TLB */ \
		do{ \
			do{ \
				r4_x2cl_eight_reals_first_fft_scratch(rsi, 64, 24*clmblkdst8, rdi); \
			}while(--loopB & (2*clm-1)); \
			newxloop2(clm, b1b, clmblkdst, rdi, XMM_SCD3); /* 2*2*clm iters prefetching 2*clm cache lines */ \
		}while(loopB & (2*2*clm-1)); \
	} \
	rsi += -4*2*clmblkdst+clmblkdst8; \
} \
rsi += -24*clmblkdst8;	/* Restore source pointer */ \
rdx = rdi;		/* Next real data sin/cos ptr */ \
pfing(rcx += -4*blkdst);		/* We've just prefetched 96 blocks, correct for */ \
;					/* code below that prefetches 4 blocks too many */ \
end_timer(22); \
\
/* Do FFT levels 4,5 */ \
/* */ \
/* On input the 64-byte cache lines hold these data values: */ \
/*	0K	192K	384K	576K */ \
/*	1K	... */ \
/*	... */ \
/*	191K	... */ \
/*	768K	... */ \
/*	... */ \
/* On output the 64-byte cache lines hold these data values: */ \
/*	0K	48K	96K	144K */ \
/*	1K	... */ \
/*	... */ \
/*	47K	... */ \
/*	192K	... */ \
/*	... */ \
\
/* Do 48 eight_reals_fft operations */ \
/* Do 48 four_complex_fft operations */ \
/*	distance between fft data elements is 48K */ \
/*	do 48 macros each processing 32 data values */ \
/* Do 96 four_complex_fft operations */ \
/*	distance between fft data elements is 48K */ \
/*	do 48 macros each processing 32 data values */ \
\
start_timer(23); \
rdi = (rdx+48*XMM_SCD2);	/* Complex data sin/cos ptr */ \
loops_init(48*2*clm);		/* 6 iters of 8 of 2*clm */ \
r4_x4cl_eight_reals_four_complex_djbfft_preload; \
for(unsigned int loopA = 6; loopA; loopA--) { \
	for(unsigned int loopB = 4*2*2*clm; loopB; ) { \
		xtouch(xptr(rcx+rbp+(clm-1)*128));	/* Preload TLB */ \
		do{ \
			do{ \
				r4_x4cl_eight_reals_four_complex_djbfft(rsi, 64, 6*clmblkdst8, 12*clmblkdst8, rdi, rdx); \
			}while(--loopB & (2*clm-1)); \
			newxloop2(clm, b2b, clmblkdst, rdi, XMM_SCD2, rdx, XMM_SCD2); /* 2*2*clm iterations */ \
		}while(loopB & (2*2*clm-1)); \
	} \
	rsi += -8*clmblkdst+clmblkdst8; \
} \
rsi += -6*clmblkdst8+24*clmblkdst8; \
rdi += -48*XMM_SCD2;	/* Restore sin/cos ptr */ \
\
loops_init(48*2*clm);		/* 6 iters of 8 of 2*clm */ \
r4_x4cl_four_complex_djbfft_preload; \
for(unsigned int loopA = 6; loopA; loopA--) { \
	for(unsigned int loopB = 8*2*clm; loopB; ) { \
		xtouch(xptr(rcx+rbp+(clm-1)*128));	/* Preload TLB */ \
		do{ \
			r4_x4cl_four_complex_djbfft(rsi, 64, 6*clmblkdst8, 12*clmblkdst8, rdi); \
		}while(--loopB & (2*clm-1)); \
		newxloop1(clm, b3b, clmblkdst, rdi, XMM_SCD2); /* 2*clm iterations */ \
	} \
	rsi += -8*clmblkdst+clmblkdst8; /* Next source pointer */ \
} \
rsi += -6*clmblkdst8-24*clmblkdst8; /* Next source pointer */ \
end_timer(23); \
\
/* Do FFT levels 6,7 (actually 1.585 levels) */ \
/* */ \
/* On input the 64-byte cache lines hold these data values: */ \
/*	0K	48K	96K	144K */ \
/*	1K	... */ \
/*	... */ \
/*	47K	... */ \
/*	192K	... */ \
/*	... */ \
/* On output the 64-byte cache lines hold these data values: */ \
/*	0K	16K	32K	48K */ \
/*	1K	... */ \
/*	... */ \
/*	15K	... */ \
/*	64K	... */ \
/*	... */ \
\
/* Do 16 six_reals_fft operations */ \
/* Do 16 three_complex_fft operations */ \
/*	distance between fft data elements is 16K */ \
/*	do 16 macros each processing 24 data values */ \
/* Do 224 three_complex_fft operations */ \
/*	distance between fft data elements is 16K */ \
/*	do 112 macros each processing 24 data values */ \
\
start_timer(24); \
rdx = (uintptr_t)g->premult_addr;/* Address of premultiplier/sin/cos table */ \
rdx += 4*2*clm*XMM_SCND2 + 4*2*clm*XMM_SCND2R + 2*clm*XMM_SCD2 + 16*2*clm*XMM_SCD4; /* Real data sin/cos ptr */ \
rdi = (rdx+16*2*clm*XMM_SCD1); /* Complex data sin/cos ptr */ \
loops_init(128*2*clm);		/* 8 iters of 2 iters of 8 of 2*clm */ \
r3_x3cl_six_reals_three_complex_djbfft_preload; \
for(unsigned int loopA = 2; loopA; loopA--) { \
	for(unsigned int loopB = 8*2*clm; loopB; ) { \
		xtouch(xptr(rcx+rbp+(clm-1)*128));	/* Preload TLB */ \
		do{ \
			r3_x3cl_six_reals_three_complex_djbfft(rsi, 64, 2*clmblkdst8, rdi, rdx); \
			rdi += XMM_SCD1;		/* Next sine/cosine pointer */ \
			rdx += XMM_SCD1;		/* Next sine/cosine pointer */ \
		}while(--loopB & (2*clm-1)); \
		newxloop1(clm, b4b, clmblkdst);	/* 2*clm iterations */ \
	} \
	rsi += -8*clmblkdst+clmblkdst8; \
} \
rsi += -2*clmblkdst8+6*clmblkdst8; /* Next source pointer */ \
rdi += -16*2*clm*XMM_SCD1;	/* Restore sin/cos ptr */ \
\
loops_reset; \
r3_x3cl_three_complex_djbfft_preload; \
for(unsigned int loopA = 7; loopA; loopA--) { \
	for(unsigned int loopB = 2; loopB; loopB--) { \
		for(unsigned int loopC = 8*2*clm; loopC; ) { \
			xtouch(xptr(rcx+rbp+(clm-1)*128));	/* Preload TLB */ \
			do{ \
				r3_x3cl_three_complex_djbfft(rsi, 64, 2*clmblkdst8, rdi); \
				rdi += XMM_SCD1;		/* Next sine/cosine pointer */ \
			}while(--loopC & (2*clm-1)); \
			newxloop1(clm, b5b, clmblkdst);	/* 2*clm iterations */ \
		} \
		rsi += -8*clmblkdst+clmblkdst8; \
	} \
	rsi += -2*clmblkdst8+6*clmblkdst8; /* Restore source pointer */ \
	rdi += -16*2*clm*XMM_SCD1;	/* Restore sin/cos ptr */ \
} \
rsi += -8*6*clmblkdst8;	/* Restore source pointer */ \
end_timer(24); \
\
/* Do FFT levels 8,9 */ \
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
/* Do 184 four_complex_fft operations */ \
/*	distance between fft data elements is 4K */ \
/*	do 92 macros each processing 32 data values */ \
\
start_timer(25); \
rdi += -16*2*clm*XMM_SCD1-4*2*clm*XMM_SCND2; /* Complex data sin/cos ptr */ \
rdx = (rdi-4*2*clm*XMM_SCND2R); /* Real data sin/cos ptr */ \
loops_init(96*2*clm);		/* 24 iterations of 4 iters of 2*clm */ \
r4_x4cl_wpn_eight_reals_four_complex_djbfft_preload; \
for(unsigned int loopA = 4; loopA; loopA--) { \
	for(unsigned int loopB = 2*clm; loopB; loopB--) { \
		r4_x4cl_wpn_eight_reals_four_complex_djbfft(rsi, 64, 4*clmblkdst, clmblkdst8, rdi, rdx); \
		rdi += XMM_SCND2;		/* Next sine/cosine pointer */ \
		rdx += XMM_SCND2R;		/* Next sine/cosine pointer */ \
	} \
	newxloop1nopre(clm, b6b, clmblkdst);	/* 2*clm iterations */ \
} \
rsi += -4*clmblkdst+2*clmblkdst8; \
rdi += -4*2*clm*XMM_SCND2;	/* Restore sin/cos ptr */ \
\
loops_reset; \
r4_x4cl_wpn_four_complex_djbfft_preload; \
for(unsigned int loopA = 23; loopA; loopA--) { \
	for(unsigned int loopB = 4*2*clm; loopB; ) { \
		xtouch(xptr(rcx+rbp+(clm-1)*128));	/* Preload TLB */ \
		do{ \
			r4_x4cl_wpn_four_complex_djbfft(rsi, 64, 4*clmblkdst, clmblkdst8, rdi); \
			rdi += XMM_SCND2;		/* Next sine/cosine pointer */ \
		}while(--loopB & (2*clm-1)); \
		newxloop1(clm, b7b, clmblkdst);	/* 2*clm iterations */ \
	} \
	rsi += -4*clmblkdst+2*clmblkdst8; /* Next source pointer */ \
	rdi += -4*2*clm*XMM_SCND2;	/* Restore sin/cos ptr */ \
} \
rsi += -24*2*clmblkdst8;	/* Restore source pointer */ \
end_timer(25); \
\
/* Do FFT levels 10,11 */ \
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
/* Do 1 four_complex_fft4 operations */ \
/*	distance between fft data elements is 1K */ \
/*	do 1 macros each processing 32 data values */ \
/* Do 190 four_complex_fft4 operations */ \
/*	distance between fft data elements is 1K */ \
/*	do 95 macros each processing 32 data values */ \
\
start_timer(26); \
rdi += -4*2*clm*XMM_SCND2R-16*2*clm*XMM_SCD4; /* Complex data sin/cos ptr */ \
rdx = (rdi-2*clm*XMM_SCD2); /* Real data sin/cos ptr */ \
rcx = (uintptr_t)g->data_addr;/* Destination pointer */ \
loops_init_prefetch(96*2*clm, 64, 2, rbp); /* 96 iters of 2*clm */ \
pfing(rbp = (uintptr_t)g->premult_prefetch);/* Load premult/sin/cos prefetch pointer */ \
r4_sg4cl_eight_reals_four_complex_fft4_preload; \
for(unsigned int loopA = 2*clm; loopA; loopA--) { \
	r4_sg4cl_eight_reals_four_complex_fft4(rsi, 64, clmblkdst, 2*clmblkdst, rcx, 64, blkdst, 2*blkdst, rdi, rdx); \
	rdi += XMM_SCD4;		/* Next sine/cosine pointer */ \
	rdx += XMM_SCD2;		/* Next sine/cosine pointer */ \
} \
rsi += -2*clm*64+4*clmblkdst; /* Next source pointer */ \
rcx += -2*clm*64+4*blkdst;	/* Next dest pointer */ \
rdi += -2*clm*XMM_SCD4;	/* Restore premult/sine/cosine pointer */ \
\
loops_reset; \
for(unsigned int loopA = 16*3*2-1; loopA; ) { \
	for(unsigned int loopB = 3; loopB; loopB--) { \
		do{ \
			for(unsigned int loopD = 2*clm; loopD; loopD--) { \
				r4_sg4cl_four_complex_fft4(rsi, 64, clmblkdst, 2*clmblkdst, rcx, 64, blkdst, 2*blkdst, rdi); \
				rdi += XMM_SCD4;		/* Next premult/sin/cos ptr */ \
			} \
			rsi += -2*clm*64+4*clmblkdst; /* Next source pointer */ \
			rcx += -2*clm*64+4*blkdst;	/* Next dest pointer */ \
			rdi += -2*clm*XMM_SCD4;	/* Restore premult/sine/cosine pointer */ \
		}while(--loopA & 1); \
		rsi += -2*4*clmblkdst+clmblkdst8; /* Next source pointer */ \
	} \
	rdi += 2*clm*XMM_SCD4;	/* Next premult/sine/cosine pointer */ \
} \
end_timer(26); \
\
/* Work on next 4*clm sets of 1536 values */ \
\
endlp:	xpass1_get_next_block(pass2, c0b, b0b, _wpn); /* Get next block and jump */ \
/* done: */ xfft_footer_rn; \
\
\
\
/* Perform a pass 1 FFT of length 1792 doubles. */ \
/* Memory examples are for a 1792K FFT with 10 levels done in pass 2. */ \
/* Sin/cos data needed: 16*2*clm*XMM_SCND2 + 16*2*clm*XMM_SCND2R + 4*2*clm*XMM_SCD2 + 4*2*clm*XMM_SCD2 + */ \
/*			2*clm*XMM_SCD2 + 14*2*clm*XMM_SCD4 = 93*2*clm 64-byte cache lines */ \
\

#define r4dwpn_pass1sc1792(procname, xpass2name, pass2_levels, clm) /* (fftname, pass2_levels, clmarg) */ \
\
/* set_FFT_constants(fftname, 0, 0, pass2_levels, clmarg); */ \
xfft_header_rn(procname, xpass2name, pass2_levels, clm); \
\
/* Do pass 2 (the last 10 FFT levels if this is a 1792K FFT) */ \
\
/* pass2: */start_timer(0); \
xpass2name(g); \
end_timer(0); \
if(g->ffttype == 1)		/* We're done if FFTing only */ \
	return; /* goto done; */ \
\
/* Do the inverse FFT */ \
\
pass1_inverse_fft_setup(c0b, procname); \
\
/* Do inverse FFT levels 10,11 */ \
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
/* Do 223 four_complex_unfft4 operations */ \
/*	distance between fft data elements is 1K */ \
/*	do 223 macros each processing 16 data values */ \
\
c0b:	start_timer(16); \
rsi = (uintptr_t)g->data_addr;/* Get source address */ \
rcx = (uintptr_t)g->scratch_area;	/* Get address of scratch area */ \
rdx = (uintptr_t)g->premult_addr;/* Real data sin/cos pointer */ \
rdi = (rdx+2*clm*XMM_SCD2); /* Complex data sin/cos ptr */ \
pfing(rbp = (rdx+(93-56)*2*clm*64)); /* Finish prefetching the sin/cos data */ \
loops_init_prefetch(224*2*clm, 64, 4, rbp); /* 224*2*clm iterations */ \
/*bug	xtouchpages rbp, 224*2*clm * 64 / 4 ;; Prefetch TLBs */ \
for(unsigned int loopA = 2*clm; loopA; loopA--) { \
	r4_sg2cl_eight_reals_unfft4(rsi, 64, blkdst, rcx, 64, clmblkdst, rdi, rdx); \
	rdi += XMM_SCD4;		/* Next premult/sin/cos data */ \
	rdx += XMM_SCD2;		/* Next sine/cosine pointer */ \
} \
rsi += -2*clm*64+2*blkdst;	/* Next source pointer */ \
rcx += -2*clm*64+2*clmblkdst; /* Next dest pointer */ \
rdi += -2*clm*XMM_SCD4;	/* Restore premult/sine/cosine pointer */ \
\
loops_reset; \
for(unsigned int loopA = 14*4*4-1; loopA; ) { \
	do{ \
		do{ \
			for(unsigned int loopD = 2*clm; loopD; loopD--) { \
				r4_sg2cl_four_complex_unfft4(rsi, 64, blkdst, rcx, 64, clmblkdst, rdi); \
				rdi += XMM_SCD4;		/* Next premult/sin/cos data */ \
			} \
			rsi += -2*clm*64+2*blkdst;	/* Next source pointer */ \
			rcx += -2*clm*64+2*clmblkdst; /* Next dest pointer */ \
			rdi += -2*clm*XMM_SCD4;	/* Restore premult/sin/cos data pointer */ \
		}while(--loopA & 3); \
		rcx += -4*2*clmblkdst+clmblkdst8; \
	}while(loopA & (4*4-1)); \
	rdi += 2*clm*XMM_SCD4;	/* Next premult/sin/cos pointer */ \
} \
rsi = (uintptr_t)g->scratch_area;	/* Restore source pointer */ \
rdx = rdi;		/* Next real data sin/cos pointer */ \
end_timer(16); \
\
/* Do inverse FFT levels 8,9 */ \
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
/* Do 220 four_complex_unfft operations */ \
/*	distance between fft data elements is 4K */ \
/*	do 110 macros each processing 32 data values */ \
\
start_timer(17); \
rdi = (rdx+4*2*clm*XMM_SCD2); /* Complex data sin/cos ptr */ \
loops_init(112*2*clm);		/* 112*2*clm iterations */ \
\
r4_x4cl_eight_reals_unfft_preload; \
for(unsigned int loopA = 2; loopA; loopA--) { \
	for(unsigned int loopB = 2*clm; loopB; loopB--) { \
		r4_x4cl_eight_reals_unfft(rsi, 64, 2*clmblkdst, 4*clmblkdst, rdi, 2*2*clm*XMM_SCD2, rdx, 2*2*clm*XMM_SCD2); \
		rdi += XMM_SCD2;		/* Next sine/cosine pointer */ \
		rdx += XMM_SCD2;		/* Next sine/cosine pointer */ \
	} \
	rsi += -2*clm*64+clmblkdst; /* Next source pointer */ \
} \
rsi += -2*clmblkdst+clmblkdst8; /* Next source pointer */ \
rdi += -2*2*clm*XMM_SCD2;	/* Restore sine/cosine pointer */ \
\
loops_reset; \
r4_x4cl_four_complex_djbunfft_preload; \
for(unsigned int loopA = 55; loopA; loopA--) { \
	for(unsigned int loopB = 2; loopB; loopB--) { \
		for(unsigned int loopC = 2*clm; loopC; loopC--) { \
			r4_x4cl_four_complex_djbunfft(rsi, 64, 2*clmblkdst, 4*clmblkdst, rdi, 2*2*clm*XMM_SCD2); \
			rdi += XMM_SCD2;		/* Next clm's sin/cos data */ \
		} \
		rsi += -2*clm*64+clmblkdst; /* Next source pointer */ \
	} \
	rsi += -2*clmblkdst+clmblkdst8; /* Next source pointer */ \
	rdi += -2*2*clm*XMM_SCD2;	/* Restore sine/cosine pointer */ \
} \
rsi += -56*clmblkdst8;	/* Restore source pointer */ \
rdx = (rdi+4*2*clm*XMM_SCD2); /* Next real data sin/cos ptr */ \
end_timer(17); \
\
/* Do inverse FFT levels 6,7 */ \
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
/*	128K	... */ \
/*	... */ \
\
/* Do 16 eight_reals_unfft operations */ \
/*	distance between fft data elements is 16K */ \
/*	do 8 macros each processing 32 data values */ \
/* Do 208 four_complex_unfft operations */ \
/*	distance between fft data elements is 16K */ \
/*	do 104 macros each processing 32 data values */ \
\
start_timer(18); \
rdi = (rdx+16*2*clm*XMM_SCND2R); /* Complex data sin/cos ptr */ \
loops_init(112*2*clm);		/* 112*2*clm iterations */ \
\
r4_x4cl_wpn_eight_reals_unfft_preload; \
for(unsigned int loopA = 8; loopA; loopA--) { \
	for(unsigned int loopB = 2*clm; loopB; loopB--) { \
		r4_x4cl_wpn_eight_reals_unfft(rsi, 64, clmblkdst8, 2*clmblkdst8, rdi, 8*2*clm*XMM_SCND2, rdx, 8*2*clm*XMM_SCND2R); \
		rdi += XMM_SCND2;		/* Next sine/cosine pointer */ \
		rdx += XMM_SCND2R;		/* Next sine/cosine pointer */ \
	} \
	rsi += -2*clm*64+clmblkdst; /* Next source pointer */ \
} \
rsi += -8*clmblkdst+4*clmblkdst8; /* Next source pointer */ \
rdi += -8*2*clm*XMM_SCND2;	/* Restore sine/cosine pointer */ \
\
loops_reset; \
r4_x4cl_wpn_four_complex_djbunfft_preload; \
for(unsigned int loopA = 13; loopA; loopA--) { \
	for(unsigned int loopB = 8; loopB; loopB--) { \
		for(unsigned int loopC = 2*clm; loopC; loopC--) { \
			r4_x4cl_wpn_four_complex_djbunfft(rsi, 64, clmblkdst8, 2*clmblkdst8, rdi, 8*2*clm*XMM_SCND2); \
			rdi += XMM_SCND2;		/* Next clm's sin/cos data */ \
		} \
		rsi += -2*clm*64+clmblkdst; /* Next source pointer */ \
	} \
	rsi += -8*clmblkdst+4*clmblkdst8; /* Next source pointer */ \
	rdi += -8*2*clm*XMM_SCND2;	/* Restore sine/cosine pointer */ \
} \
rsi += -14*4*clmblkdst8;	/* Restore source pointer */ \
end_timer(18); \
\
/* Do inverse FFT levels 1,2,3,4,5 (actually 4.807 levels) */ \
/* On input the 64-byte cache lines hold these data values: */ \
/*	0K	32K	64K	96K */ \
/*	1K	... */ \
/*	... */ \
/*	31K */ \
/*	128K	... */ \
/*	... */ \
/* On output the 64-byte cache lines hold these data values: */ \
/*	0K	448K	896K	1344K */ \
/*	1K	... */ \
/*	... */ \
/*	447K */ \
\
/* Do 256 28_reals_last_unfft operations */ \
/*	distance between fft data elements is 64K */ \
/*	do 32 macros each processing 112 data values */ \
\
start_timer(20); \
rdi = (uintptr_t)g->sincos2;/* Load fixed sin/cos postmultipliers */ \
set_data_prefetch_ptrs;		/* Init prefetching pointers */ \
loops_init(32*2*clm);		/* 32 iterations of 2*clm */ \
r7_x14cl_28_reals_last_unfft_preload; \
for(unsigned int loopA = 4; loopA; loopA--) { \
	for(unsigned int loopB = 8*2*clm; loopB; ) { \
		xtouch(xptr(rcx+rbp+(clm-1)*128));	/* Preload TLB */ \
		do{ \
			r7_x14cl_28_reals_last_unfft(rsi, 64, 4*clmblkdst8, rdi, 32*XMM_SCD13); \
		}while(--loopB & (2*clm-1)); \
		newxloop1quadpre(clm, c1b, c1a, clmblkdst, rdi, XMM_SCD13); /* 2*clm iterations */ \
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
copy_scratch_data_to_fft(448, clm); \
goto	endlp;			/* Skip FFT code if flag not set */ \
\
/* Fire up auxiliary threads to do pass 1 forward FFT */ \
\
/* pass1: */pass1_forward_fft_setup(b0b); \
\
/* Copy data from the FFT data area back to the scratch area */ \
\
b0b:	copy_fft_data_to_scratch(448, clm); \
\
/* Do FFT levels 1,2,3,4,5 (actually 4.807 levels) */ \
/* */ \
/* On input the 64-byte cache lines hold these data values: */ \
/*	0K	448K	896K	1344K */ \
/*	1K	... */ \
/*	... */ \
/*	447K */ \
/* On output the 64-byte cache lines hold these data values: */ \
/*	0K	64K	128K	192K */ \
/*	1K	... */ \
/*	... */ \
/*	63K	... */ \
/*	256K	... */ \
/*	... */ \
\
/* Do 256 28_reals_first_fft operations */ \
/*	distance between fft data elements is 64K */ \
/*	do 64 iterations each processing 56 data values */ \
\
yes_fft:start_timer(22); \
rsi = (uintptr_t)g->scratch_area;	/* Get address of scratch area */ \
rdi = (uintptr_t)g->sincos2;/* Load fixed sin/cos postmultipliers */ \
set_data_prefetch_ptrs;		/* Init prefetching pointers */ \
pfing(rcx += (128-12)*blkdst);	/* We've already prefetched 128 blocks, but */ \
;					/* code below will fetch 12 blocks too many */ \
loops_init(64*2*clm);		/* 64 iters of 2*clm */ \
r7_x7cl_28_reals_first_fft_scratch_preload; \
for(unsigned int loopA = 8; loopA; loopA--) { \
	for(unsigned int loopB = 8*2*clm; loopB; ) { \
		xtouch(xptr(rcx+rbp+(clm-1)*128));	/* Preload TLB */ \
		do{ \
			r7_x7cl_28_reals_first_fft_scratch(rsi, 64, 8*clmblkdst8, rdi); \
		}while(--loopB & (2*clm-1)); \
		newxloop1doublepre(clm, b1b, b1a, clmblkdst, rdi, XMM_SCD13); /* 2*clm iters prefetching 4*clm cache lines */ \
	} \
	rsi += -8*clmblkdst+clmblkdst8; \
} \
rsi += -8*clmblkdst8;	/* Restore source pointer */ \
end_timer(22); \
\
/* Do FFT levels 6,7 */ \
/* */ \
/* On input the 64-byte cache lines hold these data values: */ \
/*	0K	64K	128K	192K */ \
/*	1K	... */ \
/*	... */ \
/*	63K	... */ \
/*	256K	... */ \
/*	... */ \
/* On output the 64-byte cache lines hold these data values: */ \
/*	0K	16K	32K	48K */ \
/*	1K	... */ \
/*	... */ \
/*	15K	... */ \
/*	64K	... */ \
/*	... */ \
\
/* Do 16 eight_reals_fft operations */ \
/* Do 16 four_complex_fft operations */ \
/*	distance between fft data elements is 16K */ \
/*	do 16 macros each processing 32 data values */ \
/* Do 192 four_complex_fft operations */ \
/*	distance between fft data elements is 16K */ \
/*	do 96 macros each processing 32 data values */ \
\
start_timer(24); \
rdx = (uintptr_t)g->premult_addr;/* Address of premultiplier/sin/cos table */ \
rdx += 4*2*clm*XMM_SCD2 + 4*2*clm*XMM_SCD2 + 2*clm*XMM_SCD2 + 14*2*clm*XMM_SCD4; /* Real data sin/cos ptr */ \
rdi = (rdx+16*2*clm*XMM_SCND2R); /* Complex data sin/cos ptr */ \
loops_init(112*2*clm);		/* 7 iterations of 16 iters of 2*clm */ \
r4_x4cl_wpn_eight_reals_four_complex_djbfft_preload; \
for(unsigned int loopA = 2; loopA; loopA--) { \
	for(unsigned int loopB = 8; loopB; loopB--) { \
		for(unsigned int loopC = 2*clm; loopC; loopC--) { \
			r4_x4cl_wpn_eight_reals_four_complex_djbfft(rsi, 64, 2*clmblkdst8, 4*clmblkdst8, rdi, rdx); \
			rdi += XMM_SCND2;		/* Next sine/cosine pointer */ \
			rdx += XMM_SCND2R;		/* Next sine/cosine pointer */ \
		} \
		newxloop1nopre(clm, b4b, clmblkdst);	 /* 2*clm iterations */ \
	} \
	rsi += -8*clmblkdst+clmblkdst8; \
} \
rsi += -2*clmblkdst8+8*clmblkdst8; \
rdi += -16*2*clm*XMM_SCND2; /* Restore sin/cos ptr */ \
\
loops_reset; \
r4_x4cl_wpn_four_complex_djbfft_preload; \
for(unsigned int loopA = 6; loopA; loopA--) { \
	for(unsigned int loopB = 2; loopB; loopB--) { \
		for(unsigned int loopC = 8*2*clm; loopC; ) { \
			xtouch(xptr(rcx+rbp+(clm-1)*128));	/* Preload TLB */ \
			do{ \
				r4_x4cl_wpn_four_complex_djbfft(rsi, 64, 2*clmblkdst8, 4*clmblkdst8, rdi); \
				rdi += XMM_SCND2;		/* Next sine/cosine pointer */ \
			}while(--loopC & (2*clm-1)); \
			newxloop1(clm, b5b, clmblkdst);	/* 2*clm iterations */ \
		} \
		rsi += -8*clmblkdst+clmblkdst8; /* Next source pointer */ \
	} \
	rsi += -2*clmblkdst8+8*clmblkdst8; /* Next source pointer */ \
	rdi += -16*2*clm*XMM_SCND2; /* Restore sin/cos ptr */ \
} \
rsi += -7*8*clmblkdst8;	/* Restore source pointer */ \
end_timer(24); \
\
/* Do FFT levels 8,9 */ \
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
/* Do 216 four_complex_fft operations */ \
/*	distance between fft data elements is 4K */ \
/*	do 108 macros each processing 32 data values */ \
\
start_timer(25); \
rdi += -16*2*clm*XMM_SCND2R-4*2*clm*XMM_SCD2; /* Complex data sin/cos ptr */ \
rdx = (rdi-4*2*clm*XMM_SCD2); /* Real data sin/cos ptr */ \
loops_init(112*2*clm);		/* 28 iterations of 4 iters of 2*clm */ \
r4_x4cl_eight_reals_four_complex_djbfft_preload; \
for(unsigned int loopA = 4; loopA; loopA--) { \
	for(unsigned int loopB = 2*clm; loopB; loopB--) { \
		r4_x4cl_eight_reals_four_complex_djbfft(rsi, 64, 4*clmblkdst, clmblkdst8, rdi, rdx); \
		rdi += XMM_SCD2;		/* Next sine/cosine pointer */ \
		rdx += XMM_SCD2;		/* Next sine/cosine pointer */ \
	} \
	newxloop1nopre(clm, b6b, clmblkdst);	/* 2*clm iterations */ \
} \
rsi += -4*clmblkdst+2*clmblkdst8; \
rdi += -4*2*clm*XMM_SCD2;	/* Restore sin/cos ptr */ \
\
loops_reset; \
r4_x4cl_four_complex_djbfft_preload; \
for(unsigned int loopA = 27; loopA; loopA--) { \
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
rsi += -28*2*clmblkdst8;	/* Restore source pointer */ \
end_timer(25); \
\
/* Do FFT levels 10,11 */ \
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
/* Do 1 four_complex_fft4 operations */ \
/*	distance between fft data elements is 1K */ \
/*	do 1 macros each processing 32 data values */ \
/* Do 222 four_complex_fft4 operations */ \
/*	distance between fft data elements is 1K */ \
/*	do 111 macros each processing 32 data values */ \
\
start_timer(26); \
rdi += -4*2*clm*XMM_SCD2-14*2*clm*XMM_SCD4; /* Complex data sin/cos ptr */ \
rdx = (rdi-2*clm*XMM_SCD2); /* Real data sin/cos ptr */ \
rcx = (uintptr_t)g->data_addr;/* Destination pointer */ \
loops_init_prefetch(112*2*clm, 64, 2, rbp); /* 112 iters of 2*clm */ \
pfing(rbp = (uintptr_t)g->premult_prefetch);/* Load premult/sin/cos prefetch pointer */ \
r4_sg4cl_eight_reals_four_complex_fft4_preload; \
for(unsigned int loopA = 2*clm; loopA; loopA--) { \
	r4_sg4cl_eight_reals_four_complex_fft4(rsi, 64, clmblkdst, 2*clmblkdst, rcx, 64, blkdst, 2*blkdst, rdi, rdx); \
	rdi += XMM_SCD4;		/* Next sine/cosine pointer */ \
	rdx += XMM_SCD2;		/* Next sine/cosine pointer */ \
} \
rsi += -2*clm*64+4*clmblkdst; /* Next source pointer */ \
rcx += -2*clm*64+4*blkdst;	/* Next dest pointer */ \
rdi += -2*clm*XMM_SCD4;	/* Restore premult/sine/cosine pointer */ \
\
loops_reset; \
for(unsigned int loopA = 14*4*2-1; loopA; ) { \
	do{ \
		do{ \
			for(unsigned int loopD = 2*clm; loopD; loopD--) { \
				r4_sg4cl_four_complex_fft4(rsi, 64, clmblkdst, 2*clmblkdst, rcx, 64, blkdst, 2*blkdst, rdi); \
				rdi += XMM_SCD4;		/* Next premult/sin/cos ptr */ \
			} \
			rsi += -2*clm*64+4*clmblkdst; /* Next source pointer */ \
			rcx += -2*clm*64+4*blkdst;	/* Next dest pointer */ \
			rdi += -2*clm*XMM_SCD4;	/* Restore premult/sine/cosine pointer */ \
		}while(--loopA & 1); \
		rsi += -2*4*clmblkdst+clmblkdst8; /* Next source pointer */ \
	}while(loopA & (4*2-1)); \
	rdi += 2*clm*XMM_SCD4;	/* Next premult/sine/cosine pointer */ \
} \
end_timer(26); \
\
/* Work on next 4*clm sets of 1792 values */ \
\
endlp:	xpass1_get_next_block(pass2, c0b, b0b, _wpn); /* Get next block and jump */ \
/* done: */ xfft_footer_rn; \
\
\
\
/* Perform a pass 1 FFT of length 2048 doubles. */ \
/* Memory examples are for a 2048K FFT with 10 levels done in pass 2. */ \
/* Sin/cos data needed: 16*2*clm*XMM_SCND2 + 16*2*clm*XMM_SCND2R + 4*2*clm*XMM_SCD2 + 4*2*clm*XMM_SCD2 + */ \
/*			2*clm*XMM_SCD2 + 16*2*clm*XMM_SCD4 = 97*2*clm 64-byte cache lines */ \
\

#define r4dwpn_pass1sc2048(procname, xpass2name, pass2_levels, clm) /* (fftname, pass2_levels, clmarg) */ \
\
/* set_FFT_constants(fftname, 0, 0, pass2_levels, clmarg); */ \
xfft_header_rn(procname, xpass2name, pass2_levels, clm); \
\
/* Do pass 2 (the last 10 FFT levels if this is a 2048K FFT) */ \
\
/* pass2: */start_timer(0); \
xpass2name(g); \
end_timer(0); \
if(g->ffttype == 1)		/* We're done if FFTing only */ \
	return; /* goto done; */ \
\
/* Do the inverse FFT */ \
\
pass1_inverse_fft_setup(c0b, procname); \
\
/* Do inverse FFT levels 10.11 */ \
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
/* Do 255 four_complex_unfft4 operations */ \
/*	distance between fft data elements is 1K */ \
/*	do 255 macros each processing 16 data values */ \
\
c0b:	start_timer(16); \
rsi = (uintptr_t)g->data_addr;/* Get source address */ \
rcx = (uintptr_t)g->scratch_area;	/* Get address of scratch area */ \
rdx = (uintptr_t)g->premult_addr;/* Real data sin/cos pointer */ \
rdi = (rdx+2*clm*XMM_SCD2); /* Complex data sin/cos ptr */ \
pfing(rbp = (rdx+(97-64)*2*clm*64)); /* Finish prefetching the sin/cos data */ \
loops_init_prefetch_touch(256*2*clm, 64, 4, rbp, 4096, c9a); /* 256*2*clm iterations */ \
xtouch(xptr(rbp));			/* Prefetch TLB */ \
for(unsigned int loopA = 2*clm; loopA; loopA--) { \
	r4_sg2cl_eight_reals_unfft4(rsi, 64, blkdst, rcx, 64, clmblkdst, rdi, rdx); \
	rdi += XMM_SCD4;		/* Next premult/sin/cos data */ \
	rdx += XMM_SCD2;		/* Next sine/cosine pointer */ \
} \
rsi += -2*clm*64+2*blkdst;	/* Next source pointer */ \
rcx += -2*clm*64+2*clmblkdst; /* Next dest pointer */ \
rdi += -2*clm*XMM_SCD4;	/* Restore premult/sine/cosine pointer */ \
\
loops_reset; \
/* c9a: */ xtouch(xptr(rbp+4096));		/* Prefetch TLB */ \
for(unsigned int loopA = 16*4*4-1; loopA; ) { \
	do{ \
		do { \
			for(unsigned int loopD = 2*clm; loopD; loopD--) { \
				r4_sg2cl_four_complex_unfft4(rsi, 64, blkdst, rcx, 64, clmblkdst, rdi); \
				rdi += XMM_SCD4;		/* Next premult/sin/cos data */ \
			} \
			rsi += -2*clm*64+2*blkdst;	/* Next source pointer */ \
			rcx += -2*clm*64+2*clmblkdst; /* Next dest pointer */ \
			rdi += -2*clm*XMM_SCD4;	/* Restore premult/sin/cos data pointer */ \
		}while(--loopA & 3); \
		rcx += -4*2*clmblkdst+clmblkdst8; \
	}while(loopA & (4*4-1)); \
	rdi += 2*clm*XMM_SCD4;	/* Next premult/sin/cos pointer */ \
} \
rsi = (uintptr_t)g->scratch_area;	/* Restore source pointer */ \
rdx = rdi;		/* Next real data sin/cos pointer */ \
end_timer(16); \
\
/* Do inverse FFT levels 8,9 */ \
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
/* Do 252 four_complex_unfft operations */ \
/*	distance between fft data elements is 4K */ \
/*	do 126 macros each processing 32 data values */ \
\
start_timer(17); \
rdi = (rdx+4*2*clm*XMM_SCD2); /* Complex data sin/cos ptr */ \
loops_init(128*2*clm);		/* 128*2*clm iterations */ \
\
r4_x4cl_eight_reals_unfft_preload; \
for(unsigned int loopA = 2; loopA; loopA--) { \
	for(unsigned int loopB = 2*clm; loopB; loopB--) { \
		r4_x4cl_eight_reals_unfft(rsi, 64, 2*clmblkdst, 4*clmblkdst, rdi, 2*2*clm*XMM_SCD2, rdx, 2*2*clm*XMM_SCD2); \
		rdi += XMM_SCD2;		/* Next sine/cosine pointer */ \
		rdx += XMM_SCD2;		/* Next sine/cosine pointer */ \
	} \
	rsi += -2*clm*64+clmblkdst; /* Next source pointer */ \
} \
rsi += -2*clmblkdst+clmblkdst8; /* Next source pointer */ \
rdi += -2*2*clm*XMM_SCD2;	/* Restore sine/cosine pointer */ \
\
loops_reset; \
r4_x4cl_four_complex_djbunfft_preload; \
for(unsigned int loopA = 63; loopA; loopA--) { \
	for(unsigned int loopB = 2; loopB; loopB--) { \
		for(unsigned int loopC = 2*clm; loopC; loopC--) { \
			r4_x4cl_four_complex_djbunfft(rsi, 64, 2*clmblkdst, 4*clmblkdst, rdi, 2*2*clm*XMM_SCD2); \
			rdi += XMM_SCD2;		/* Next clm's sin/cos data */ \
		} \
		rsi += -2*clm*64+clmblkdst; /* Next source pointer */ \
	} \
	rsi += -2*clmblkdst+clmblkdst8; /* Next source pointer */ \
	rdi += -2*2*clm*XMM_SCD2;	/* Restore sine/cosine pointer */ \
} \
rsi += -64*clmblkdst8;	/* Restore source pointer */ \
rdx = (rdi+4*2*clm*XMM_SCD2); /* Next real data sin/cos ptr */ \
end_timer(17); \
\
/* Do inverse FFT levels 6,7 */ \
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
/*	128K	... */ \
/*	... */ \
\
/* Do 16 eight_reals_unfft operations */ \
/*	distance between fft data elements is 16K */ \
/*	do 8 macros each processing 32 data values */ \
/* Do 240 four_complex_unfft operations */ \
/*	distance between fft data elements is 16K */ \
/*	do 120 macros each processing 32 data values */ \
\
start_timer(18); \
rdi = (rdx+16*2*clm*XMM_SCND2R); /* Complex data sin/cos ptr */ \
loops_init(128*2*clm);		/* 128 iterations of 2*clm */ \
r4_x4cl_wpn_eight_reals_unfft_preload; \
for(unsigned int loopA = 8; loopA; loopA--) { \
	for(unsigned int loopB = 2*clm; loopB; loopB--) { \
		r4_x4cl_wpn_eight_reals_unfft(rsi, 64, clmblkdst8, 2*clmblkdst8, rdi, 8*2*clm*XMM_SCND2, rdx, 8*2*clm*XMM_SCND2R); \
		rdi += XMM_SCND2;		/* Next sine/cosine pointer */ \
		rdx += XMM_SCND2R;		/* Next sine/cosine pointer */ \
	} \
	rsi += -2*clm*64+clmblkdst; /* Next source pointer */ \
} \
rsi += -8*clmblkdst+4*clmblkdst8; /* Next source pointer */ \
rdi += -8*2*clm*XMM_SCND2;	/* Restore sine/cosine pointer */ \
\
loops_reset; \
r4_x4cl_wpn_four_complex_djbunfft_preload; \
for(unsigned int loopA = 15; loopA; loopA--) { \
	for(unsigned int loopB = 8; loopB; loopB--) { \
		for(unsigned int loopC = 2*clm; loopC; loopC--) { \
			r4_x4cl_wpn_four_complex_djbunfft(rsi, 64, clmblkdst8, 2*clmblkdst8, rdi, 8*2*clm*XMM_SCND2); \
			rdi += XMM_SCND2;		/* Next clm's sin/cos data */ \
		} \
		rsi += -2*clm*64+clmblkdst; /* Next source pointer */ \
	} \
	rsi += -8*clmblkdst+4*clmblkdst8; /* Next source pointer */ \
	rdi += -8*2*clm*XMM_SCND2;	/* Restore sine/cosine pointer */ \
} \
rsi += -16*4*clmblkdst8;	/* Restore source pointer */ \
end_timer(18); \
\
/* Do inverse FFT levels 4,5 */ \
/* On input the 64-byte cache lines hold these data values: */ \
/*	0K	32K	64K	96K */ \
/*	1K	... */ \
/*	... */ \
/*	31K */ \
/*	128K	... */ \
/*	... */ \
/* On output the 64-byte cache lines hold these data values: */ \
/*	0K	128K	256K	384K */ \
/*	1K	... */ \
/*	... */ \
/*	127K */ \
/*	512K	... */ \
/*	... */ \
\
/* Do 64 eight_reals_unfft operations */ \
/*	distance between fft data elements is 64K */ \
/*	do 32 macros each processing 32 data values */ \
/* Do 192 four_complex_unfft operations */ \
/*	distance between fft data elements is 64K */ \
/*	do 96 macros each processing 32 data values */ \
\
start_timer(19); \
rdx = (uintptr_t)g->sincos2;/* Load fixed sin/cos data ptr */ \
rdx += 256*XMM_SCD3;	/* Real data sin/cos ptr */ \
rdi = (rdx+64*XMM_SCD2);	/* Complex data sin/cos ptr */ \
loops_init(128*2*clm);		/* 128 iterations of 2*clm */ \
\
r4_x4cl_eight_reals_unfft_preload; \
for(unsigned int loopA = 4; loopA; loopA--) { \
	for(unsigned int loopB = 8; loopB; loopB--) { \
		for(unsigned int loopC = 2*clm; loopC; loopC--) { \
			r4_x4cl_eight_reals_unfft(rsi, 64, 4*clmblkdst8, 8*clmblkdst8, rdi, 32*XMM_SCD2, rdx, 32*XMM_SCD2); \
		} \
		rsi += -2*clm*64+clmblkdst; /* Next source pointer */ \
		rdi += XMM_SCD2;		/* Next sine/cosine pointer */ \
		rdx += XMM_SCD2;		/* Next sine/cosine pointer */ \
	} \
	rsi += -8*clmblkdst+clmblkdst8; /* Next source pointer */ \
} \
rsi += -4*clmblkdst8+16*clmblkdst8; /* Next source pointer */ \
rdi += -32*XMM_SCD2;	/* Restore sine/cosine pointer */ \
\
loops_reset; \
r4_x4cl_four_complex_djbunfft_preload; \
for(unsigned int loopA = 3; loopA; loopA--) { \
	for(unsigned int loopB = 4; loopB; loopB--) { \
		for(unsigned int loopC = 8; loopC; loopC--) { \
			for(unsigned int loopD = 2*clm; loopD; loopD--) { \
				r4_x4cl_four_complex_djbunfft(rsi, 64, 4*clmblkdst8, 8*clmblkdst8, rdi, 32*XMM_SCD2); \
			} \
			rsi += -2*clm*64+clmblkdst; /* Next source pointer */ \
			rdi += XMM_SCD2;		/* Next clm's sin/cos data */ \
		} \
		rsi += -8*clmblkdst+clmblkdst8; /* Next source pointer */ \
	} \
	rsi += -4*clmblkdst8+16*clmblkdst8; /* Next source pointer */ \
	rdi += -32*XMM_SCD2;	/* Restore sine/cosine pointer */ \
} \
rsi += -4*16*clmblkdst8;	/* Restore source pointer */ \
end_timer(19); \
\
/* Do inverse FFT levels 1,2,3 */ \
/* On input the 64-byte cache lines hold these data values: */ \
/*	0K	128K	256K	384K */ \
/*	1K	... */ \
/*	... */ \
/*	127K */ \
/*	512K	... */ \
/*	... */ \
/* On output the 64-byte cache lines hold these data values: */ \
/*	0K	512K	1024K	1536K */ \
/*	1K	... */ \
/*	... */ \
/*	511K */ \
\
/* Do 256 eight_reals_last_unfft operations */ \
/*	distance between fft data elements is 256K */ \
/*	do 128 macros each processing 32 data values */ \
\
start_timer(20); \
rdi = (uintptr_t)g->sincos2;/* Load fixed sin/cos postmultipliers */ \
loops_init(128*2*clm);		/* 128 iterations of 2*clm */ \
r4_x4cl_eight_reals_last_unfft_preload; \
for(unsigned int loopA = 16; loopA; loopA--) { \
	for(unsigned int loopB = 8; loopB; loopB--) { \
		for(unsigned int loopC = 2*clm; loopC; loopC--) { \
			r4_x4cl_eight_reals_last_unfft(rsi, 64, 16*clmblkdst8, 32*clmblkdst8, rdi, 128*XMM_SCD3); \
		} \
		newxloop1nopre(clm, c1b, clmblkdst, rdi, XMM_SCD3); /* 2*clm iterations */ \
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
copy_scratch_data_to_fft(512, clm); \
goto	endlp;			/* Skip FFT code if flag not set */ \
\
/* Fire up auxiliary threads to do pass 1 forward FFT */ \
\
/* pass1: */pass1_forward_fft_setup(b0b); \
\
/* Copy data from the FFT data area back to the scratch area */ \
\
b0b:	copy_fft_data_to_scratch(512, clm); \
\
/* Do FFT levels 1,2,3 */ \
/* */ \
/* On input the 64-byte cache lines hold these data values: */ \
/*	0K	512K	1024K	1536K */ \
/*	1K	... */ \
/*	... */ \
/*	511K */ \
/* On output the 64-byte cache lines hold these data values: */ \
/*	0K	128K	256K	384K */ \
/*	1K	... */ \
/*	... */ \
/*	127K */ \
\
/* Do 256 eight_reals_first_fft operations */ \
/*	distance between fft data elements is 256K */ \
/*	do 256 iterations each processing 16 data values */ \
\
yes_fft:start_timer(22); \
rsi = (uintptr_t)g->scratch_area;	/* Get address of scratch area */ \
rdi = (uintptr_t)g->sincos2;/* Load fixed sin/cos postmultipliers */ \
set_data_prefetch_ptrs;		/* Init prefetching pointers */ \
loops_init(256*2*clm);		/* 32 iters of 4 iters of 2*2*clm */ \
r4_x2cl_eight_reals_first_fft_scratch_preload; \
for(unsigned int loopA = 32; loopA; loopA--) { \
	for(unsigned int loopB = 4*2*2*clm; loopB; ) { \
		xtouch(xptr(rcx+rbp+(clm-1)*128));	/* Preload TLB */ \
		do{ \
			do{ \
				r4_x2cl_eight_reals_first_fft_scratch(rsi, 64, 32*clmblkdst8, rdi); \
			}while(--loopB & (2*clm-1)); \
			newxloop2(clm, b1b, clmblkdst, rdi, XMM_SCD3); /* 2*2*clm iters prefetching 2*clm cache lines */ \
		}while(loopB & (2*2*clm-1)); \
	} \
	rsi += -4*2*clmblkdst+clmblkdst8; \
} \
rsi += -32*clmblkdst8;	/* Restore source pointer */ \
rdx = rdi;		/* Next real data sin/cos ptr */ \
end_timer(22); \
\
/* Do FFT levels 4,5 */ \
/* */ \
/* On input the 64-byte cache lines hold these data values: */ \
/*	0K	128K	256K	384K */ \
/*	1K	... */ \
/*	... */ \
/*	127K */ \
/* On output the 64-byte cache lines hold these data values: */ \
/*	0K	64K	128K	192K */ \
/*	1K	... */ \
/*	... */ \
/*	63K	... */ \
/*	256K	... */ \
/*	... */ \
\
/* Do 64 eight_reals_fft operations */ \
/* Do 64 four_complex_fft operations */ \
/*	distance between fft data elements is 64K */ \
/*	do 64 macros each processing 32 data values */ \
/* Do 128 four_complex_fft operations */ \
/*	distance between fft data elements is 16K */ \
/*	do 64 macros each processing 32 data values */ \
\
start_timer(23); \
rdi = (rdx+64*XMM_SCD2);	/* Complex data sin/cos ptr */ \
loops_init(128*2*clm);		/* 2 sets of 8 iters of 8 of 2*clm */ \
r4_x4cl_eight_reals_four_complex_djbfft_preload; \
for(unsigned int loopA = 8; loopA; loopA--) { \
	for(unsigned int loopB = 8*2*clm; loopB; ) { \
		xtouch(xptr(rcx+rbp+(clm-1)*128));	/* Preload TLB */ \
		do{ \
			r4_x4cl_eight_reals_four_complex_djbfft(rsi, 64, 8*clmblkdst8, 16*clmblkdst8, rdi, rdx); \
		}while(--loopB & (2*clm-1)); \
		newxloop1(clm, b2b, clmblkdst, rdi, XMM_SCD2, rdx, XMM_SCD2); /* 2*clm iterations */ \
	} \
	rsi += -8*clmblkdst+clmblkdst8; \
} \
rsi += -8*clmblkdst8+32*clmblkdst8; /* Next source pointer */ \
rdi += -64*XMM_SCD2;	/* Restore sin/cos ptr */ \
\
loops_reset; \
r4_x4cl_four_complex_djbfft_preload; \
for(unsigned int loopA = 8; loopA; loopA--) { \
	for(unsigned int loopB = 8*2*clm; loopB; ) { \
		xtouch(xptr(rcx+rbp+(clm-1)*128));	/* Preload TLB */ \
		do{ \
			r4_x4cl_four_complex_djbfft(rsi, 64, 8*clmblkdst8, 16*clmblkdst8, rdi); \
		}while(--loopB & (2*clm-1)); \
		newxloop1(clm, b3b, clmblkdst, rdi, XMM_SCD2);	/* 2*clm iterations */ \
	} \
	rsi += -8*clmblkdst+clmblkdst8; \
} \
rsi += -8*clmblkdst8-32*clmblkdst8; /* Restore source pointer */ \
end_timer(23); \
\
/* Do FFT levels 6,7 */ \
/* */ \
/* On input the 64-byte cache lines hold these data values: */ \
/*	0K	64K	128K	192K */ \
/*	1K	... */ \
/*	... */ \
/*	63K	... */ \
/*	256K	... */ \
/*	... */ \
/* On output the 64-byte cache lines hold these data values: */ \
/*	0K	16K	32K	48K */ \
/*	1K	... */ \
/*	... */ \
/*	15K	... */ \
/*	64K	... */ \
/*	... */ \
\
/* Do 16 eight_reals_fft operations */ \
/* Do 16 four_complex_fft operations */ \
/*	distance between fft data elements is 16K */ \
/*	do 16 macros each processing 32 data values */ \
/* Do 224 four_complex_fft operations */ \
/*	distance between fft data elements is 16K */ \
/*	do 112 macros each processing 32 data values */ \
\
start_timer(24); \
rdx = (uintptr_t)g->premult_addr;/* Address of premultiplier/sin/cos table */ \
rdx += 4*2*clm*XMM_SCD2 + 4*2*clm*XMM_SCD2 + 2*clm*XMM_SCD2 + 16*2*clm*XMM_SCD4; /* Real data sin/cos ptr */ \
rdi = (rdx+16*2*clm*XMM_SCND2R); /* Complex data sin/cos ptr */ \
loops_init(128*2*clm);		/* 8 iters of 16 of 2*clm */ \
r4_x4cl_wpn_eight_reals_four_complex_djbfft_preload; \
for(unsigned int loopA = 2; loopA; loopA--) { \
	for(unsigned int loopB = 8*2*clm; loopB; ) { \
		xtouch(xptr(rcx+rbp+(clm-1)*128));	/* Preload TLB */ \
		do{ \
			r4_x4cl_wpn_eight_reals_four_complex_djbfft(rsi, 64, 2*clmblkdst8, 4*clmblkdst8, rdi, rdx); \
			rdi += XMM_SCND2;		/* Next sine/cosine pointer */ \
			rdx += XMM_SCND2R;		/* Next sine/cosine pointer */ \
		}while(--loopB & (2*clm-1)); \
		newxloop1(clm, b4b, clmblkdst);	/* 2*clm iterations */ \
	} \
	rsi += -4*2*clmblkdst+clmblkdst8; \
} \
rsi += -2*clmblkdst8+8*clmblkdst8; /* Next source pointer */ \
rdi += -16*2*clm*XMM_SCND2;	/* Restore sin/cos ptr */ \
\
loops_reset; \
r4_x4cl_wpn_four_complex_djbfft_preload; \
for(unsigned int loopA = 7; loopA; loopA--) { \
	for(unsigned int loopB = 2; loopB; loopB--) { \
		for(unsigned int loopC = 8*2*clm; loopC; ) { \
			xtouch(xptr(rcx+rbp+(clm-1)*128));	/* Preload TLB */ \
			do{ \
				r4_x4cl_wpn_four_complex_djbfft(rsi, 64, 2*clmblkdst8, 4*clmblkdst8, rdi); \
				rdi += XMM_SCND2;		/* Next sine/cosine pointer */ \
			}while(--loopC & (2*clm-1)); \
			newxloop1(clm, b5b, clmblkdst);	/* 2*clm iterations */ \
		} \
		rsi += -4*2*clmblkdst+clmblkdst8; \
	} \
	rsi += -2*clmblkdst8+8*clmblkdst8; /* Next source pointer */ \
	rdi += -16*2*clm*XMM_SCND2; /* Restore sin/cos ptr */ \
} \
rsi += -8*8*clmblkdst8;	/* Restore source pointer */ \
end_timer(24); \
\
/* Do FFT levels 8,9 */ \
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
/* Do 248 four_complex_fft operations */ \
/*	distance between fft data elements is 4K */ \
/*	do 124 macros each processing 32 data values */ \
\
start_timer(25); \
rdi += -16*2*clm*XMM_SCND2R-4*2*clm*XMM_SCD2; /* Complex data sin/cos ptr */ \
rdx = (rdi-4*2*clm*XMM_SCD2); /* Real data sin/cos ptr */ \
loops_init(128*2*clm);		/* 32 iterations of 4 iters of 2*clm */ \
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
for(unsigned int loopA = 31; loopA; loopA--) { \
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
rsi += -32*2*clmblkdst8;	/* Restore source pointer */ \
end_timer(25); \
\
/* Do FFT levels 10,11 */ \
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
/* Do 1 four_complex_fft4 operations */ \
/*	distance between fft data elements is 1K */ \
/*	do 1 macros each processing 32 data values */ \
/* Do 254 four_complex_fft4 operations */ \
/*	distance between fft data elements is 1K */ \
/*	do 127 macros each processing 32 data values */ \
\
start_timer(26); \
rdi += -4*2*clm*XMM_SCD2-16*2*clm*XMM_SCD4; /* Complex data sin/cos ptr */ \
rdx = (rdi-2*clm*XMM_SCD2); /* Real data sin/cos ptr */ \
rcx = (uintptr_t)g->data_addr;/* Destination pointer */ \
loops_init_prefetch(128*2*clm, 64, 2, rbp); /* 128 iters of 2*clm */ \
pfing(rbp = (uintptr_t)g->premult_prefetch);/* Load premult/sin/cos prefetch pointer */ \
r4_sg4cl_eight_reals_four_complex_fft4_preload; \
for(unsigned int loopA = 2*clm; loopA; loopA--) { \
	r4_sg4cl_eight_reals_four_complex_fft4(rsi, 64, clmblkdst, 2*clmblkdst, rcx, 64, blkdst, 2*blkdst, rdi, rdx); \
	rdi += XMM_SCD4;		/* Next sine/cosine pointer */ \
	rdx += XMM_SCD2;		/* Next sine/cosine pointer */ \
} \
rsi += -2*clm*64+4*clmblkdst; /* Next source pointer */ \
rcx += -2*clm*64+4*blkdst;	/* Next dest pointer */ \
rdi += -2*clm*XMM_SCD4;	/* Restore premult/sine/cosine pointer */ \
\
loops_reset; \
for(unsigned int loopA = 16*4*2-1; loopA; ) { \
	do{ \
		do{ \
			for(unsigned int loopD = 2*clm; loopD; loopD--) { \
				r4_sg4cl_four_complex_fft4(rsi, 64, clmblkdst, 2*clmblkdst, rcx, 64, blkdst, 2*blkdst, rdi); \
				rdi += XMM_SCD4;		/* Next premult/sin/cos ptr */ \
			} \
			rsi += -2*clm*64+4*clmblkdst; /* Next source pointer */ \
			rcx += -2*clm*64+4*blkdst;	/* Next dest pointer */ \
			rdi += -2*clm*XMM_SCD4;	/* Restore premult/sine/cosine pointer */ \
		}while(--loopA & 1); \
		rsi += -2*4*clmblkdst+clmblkdst8; /* Next source pointer */ \
	}while(loopA & (4*2-1)); \
	rdi += 2*clm*XMM_SCD4;	/* Next premult/sine/cosine pointer */ \
} \
end_timer(26); \
\
/* Work on next 4*clm sets of 2048 values */ \
\
endlp:	xpass1_get_next_block(pass2, c0b, b0b, _wpn); /* Get next block and jump */ \
/* done: */ xfft_footer_rn; \
\
\
\
/* Perform a pass 1 FFT of length 2560 doubles. */ \
/* Memory examples are for a 2560K FFT with 10 levels done in pass 2. */ \
/* Sin/cos data needed: 32*2*clm*XMM_SCD2 + 32*2*clm*XMM_SCD2 + 8*2*clm*XMM_SCND2 + 8*2*clm*XMM_SCND2R + */ \
/*			2*clm*XMM_SCD4 + 10*2*clm*XMM_SCD8 = 134*2*clm 64-byte cache lines */ \
\

#define r4dwpn_pass1sc2560(procname, xpass2name, pass2_levels, clm) /* (fftname, pass2_levels, clmarg) */ \
\
/* set_FFT_constants(fftname, 0, 0, pass2_levels, clmarg); */ \
xfft_header_rn(procname, xpass2name, pass2_levels, clm); \
\
/* Do pass 2 (the last 10 FFT levels if this is a 2560K FFT) */ \
\
/* pass2: */start_timer(0); \
xpass2name(g); \
end_timer(0); \
if(g->ffttype == 1)		/* We're done if FFTing only */ \
	return; /* goto done; */ \
\
/* Do the inverse FFT */ \
\
pass1_inverse_fft_setup(c0b, procname); \
\
/* Do inverse FFT levels 10,11,12 */ \
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
/* Do 159 eight_complex_unfft8 operations */ \
/*	distance between fft data elements is 1K */ \
/*	do 159 macros each processing 32 data values */ \
\
c0b:	start_timer(16); \
rsi = (uintptr_t)g->data_addr;/* Get source address */ \
rcx = (uintptr_t)g->scratch_area;	/* Get address of scratch area */ \
rdx = (uintptr_t)g->premult_addr;/* Real data sin/cos pointer */ \
rdi = (rdx+2*clm*XMM_SCD4+128); /* Complex data sin/cos ptr and use -128 */ \
;					/* addressing to generate 1 byte offsets in radix-8 macro */ \
pfing(rbp = (rdx+(134-40-20)*2*clm*64)); /* Finish prefetching the sin/cos data */ \
loops_init_prefetch(160*2*clm, 64, 4, rbp); /* 160*2*clm iterations */ \
/*bug	xtouchpages rbp, 160*2*clm * 64 / 4 ;; Prefetch TLBs */ \
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
for(unsigned int loopA = 10*8*2-1; loopA; ) { \
	do{ \
		do{ \
			for(unsigned int loopD = 2*clm; loopD; loopD--) { \
				r8_sg4cl_eight_complex_unfft8(rsi, 64, blkdst, 2*blkdst, rcx, 64, clmblkdst, 2*clmblkdst, rdi-128); \
				rdi += XMM_SCD8;		/* Next premult/sin/cos data */ \
			} \
			rsi += -2*clm*64+4*blkdst;	/* Next source pointer */ \
			rcx += -2*clm*64+4*clmblkdst; /* Next dest pointer */ \
			rdi += -2*clm*XMM_SCD8;	/* Restore premult/sin/cos data pointer */ \
		}while(--loopA & 1); \
		rcx += -2*4*clmblkdst+clmblkdst8; \
	}while(loopA & (8*2-1)); \
	rdi += 2*clm*XMM_SCD8;	/* Next premult/sin/cos pointer */ \
} \
rsi = (uintptr_t)g->scratch_area;	/* Restore source pointer */ \
rdx = (rdi-128);		/* Next real data sin/cos pointer */ \
end_timer(16); \
\
/* Do inverse FFT levels 8,9 */ \
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
/* Do 312 four_complex_unfft operations */ \
/*	distance between fft data elements is 8K */ \
/*	do 156 macros each processing 32 data values */ \
\
start_timer(17); \
rdi = (rdx+8*2*clm*XMM_SCND2R); /* Complex data sin/cos ptr */ \
loops_init_prefetch(160*2*clm, 64, 8, rbp); /* 160 iterations of 2*clm */ \
/*	xtouchpages rbp, 160*2*clm * 64 / 8 ;; Prefetch TLBs */ \
r4_x4cl_wpn_eight_reals_unfft_preload; \
for(unsigned int loopA = 4; loopA; loopA--) { \
	for(unsigned int loopB = 2*clm; loopB; loopB--) { \
		r4_x4cl_wpn_eight_reals_unfft(rsi, 64, 4*clmblkdst, clmblkdst8, rdi, 4*2*clm*XMM_SCND2, rdx, 4*2*clm*XMM_SCND2R); \
		rdi += XMM_SCND2;		/* Next sine/cosine pointer */ \
		rdx += XMM_SCND2R;		/* Next sine/cosine pointer */ \
	} \
	rsi += -2*clm*64+clmblkdst; /* Next source pointer */ \
} \
rsi += -4*clmblkdst+2*clmblkdst8; /* Next source pointer */ \
rdi += -4*2*clm*XMM_SCND2;	/* Restore sine/cosine pointer */ \
\
loops_reset; \
r4_x4cl_wpn_four_complex_djbunfft_preload; \
for(unsigned int loopA = 39; loopA; loopA--) { \
	for(unsigned int loopB = 4; loopB; loopB--) { \
		for(unsigned int loopC = 2*clm; loopC; loopC--) { \
			r4_x4cl_wpn_four_complex_djbunfft(rsi, 64, 4*clmblkdst, clmblkdst8, rdi, 4*2*clm*XMM_SCND2); \
			rdi += XMM_SCND2;		/* Next clm's sin/cos data */ \
		} \
		rsi += -2*clm*64+clmblkdst; /* Next source pointer */ \
	} \
	rsi += -4*clmblkdst+2*clmblkdst8; /* Next source pointer */ \
	rdi += -4*2*clm*XMM_SCND2;	/* Restore sine/cosine pointer */ \
} \
rsi += -40*2*clmblkdst8;	/* Restore source pointer */ \
rdx = (rdi+8*2*clm*XMM_SCND2); /* Next real data sin/cos ptr */ \
end_timer(17); \
\
/* Do inverse FFT levels 6,7 */ \
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
/*	256K	... */ \
/*	... */ \
\
/* Do 32 eight_reals_unfft operations */ \
/*	distance between fft data elements is 32K */ \
/*	do 16 macros each processing 32 data values */ \
/* Do 288 four_complex_unfft operations */ \
/*	distance between fft data elements is 32K */ \
/*	do 144 macros each processing 32 data values */ \
\
start_timer(18); \
rdi = (rdx+32*2*clm*XMM_SCD2); /* Complex data sin/cos ptr */ \
loops_init(160*2*clm);		/* 160 iterations of 2*clm */ \
r4_x4cl_eight_reals_unfft_preload; \
for(unsigned int loopA = 2; loopA; loopA--) { \
	for(unsigned int loopB = 8; loopB; loopB--) { \
		for(unsigned int loopC = 2*clm; loopC; loopC--) { \
			r4_x4cl_eight_reals_unfft(rsi, 64, 2*clmblkdst8, 4*clmblkdst8, rdi, 16*2*clm*XMM_SCD2, rdx, 16*2*clm*XMM_SCD2); \
			rdi += XMM_SCD2;		/* Next sine/cosine pointer */ \
			rdx += XMM_SCD2;		/* Next sine/cosine pointer */ \
		} \
		rsi += -2*clm*64+clmblkdst; /* Next source pointer */ \
	} \
	rsi += -8*clmblkdst+clmblkdst8; /* Next source pointer */ \
} \
rsi += -2*clmblkdst8+8*clmblkdst8; /* Next source pointer */ \
rdi += -16*2*clm*XMM_SCD2;	/* Restore sine/cosine pointer */ \
\
loops_reset; \
r4_x4cl_four_complex_djbunfft_preload; \
for(unsigned int loopA = 9; loopA; loopA--) { \
	for(unsigned int loopB = 2; loopB; loopB--) { \
		for(unsigned int loopC = 8; loopC; loopC--) { \
			for(unsigned int loopD = 2*clm; loopD; loopD--) { \
				r4_x4cl_four_complex_djbunfft(rsi, 64, 2*clmblkdst8, 4*clmblkdst8, rdi, 16*2*clm*XMM_SCD2); \
				rdi += XMM_SCD2;		/* Next clm's sin/cos data */ \
			} \
			rsi += -2*clm*64+clmblkdst; /* Next source pointer */ \
		} \
		rsi += -8*clmblkdst+clmblkdst8; /* Next source pointer */ \
	} \
	rsi += -2*clmblkdst8+8*clmblkdst8; /* Next source pointer */ \
	rdi += -16*2*clm*XMM_SCD2;	/* Restore sine/cosine pointer */ \
} \
rsi += -10*8*clmblkdst8;	/* Restore source pointer */ \
end_timer(18); \
\
/* Do inverse FFT levels 1,2,3,4,5 (actually 4.32 levels) */ \
/* On input the 64-byte cache lines hold these data values: */ \
/*	0K	64K	128K	192K */ \
/*	1K	... */ \
/*	... */ \
/*	63K */ \
/*	256K	... */ \
/*	... */ \
/* On output the 64-byte cache lines hold these data values: */ \
/*	0K	640K	1280K	1920K */ \
/*	1K	... */ \
/*	... */ \
/*	639K */ \
\
/* Do 512 20_reals_last_unfft operations */ \
/*	distance between fft data elements is 128K */ \
/*	do 64 macros each processing 80 data values */ \
\
start_timer(20); \
rdi = (uintptr_t)g->sincos2;/* Load fixed sin/cos postmultipliers */ \
set_data_prefetch_ptrs;		/* Init prefetching pointers */ \
loops_init(64*2*clm);		/* 64 iterations of 2*clm */ \
r5_x10cl_20_reals_last_unfft_preload; \
for(unsigned int loopA = 8; loopA; loopA--) { \
	for(unsigned int loopB = 8*2*clm; loopB; ) { \
		xtouch(xptr(rcx+rbp+(clm-1)*128));	/* Preload TLB */ \
		do{ \
			r5_x10cl_20_reals_last_unfft(rsi, 64, 8*clmblkdst8, rdi, 64*XMM_SCD9); \
		}while(--loopB & (2*clm-1)); \
		newxloop1(clm, c1b, clmblkdst, rdi, XMM_SCD9); /* 2*clm iterations */ \
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
copy_scratch_data_to_fft(640, clm); \
goto	endlp;			/* Skip FFT code if flag not set */ \
\
/* Fire up auxiliary threads to do pass 1 forward FFT */ \
\
/* pass1: */pass1_forward_fft_setup(b0b); \
\
/* Copy data from the FFT data area back to the scratch area */ \
\
b0b:	copy_fft_data_to_scratch(640, clm); \
\
/* Do FFT levels 1,2,3,4,5 (actually 4.32 levels) */ \
/* */ \
/* On input the 64-byte cache lines hold these data values: */ \
/*	0K	640K	1280K	1920K */ \
/*	1K	... */ \
/*	... */ \
/*	639K */ \
/* On output the 64-byte cache lines hold these data values: */ \
/*	0K	128K	256K	384K */ \
/*	1K	... */ \
/*	... */ \
/*	127K	... */ \
/*	512K	... */ \
/*	... */ \
\
/* Do 512 20_reals_first_fft operations */ \
/*	distance between fft data elements is 128K */ \
/*	do 128 iterations each processing 40 data values */ \
\
yes_fft:start_timer(22); \
rsi = (uintptr_t)g->scratch_area;	/* Get address of scratch area */ \
rdi = (uintptr_t)g->sincos2;/* Load fixed sin/cos postmultipliers */ \
set_data_prefetch_ptrs;		/* Init prefetching pointers */ \
pfing(rcx += (64-0)*blkdst);	/* We've already prefetched 64 blocks, perfect! */ \
loops_init(128*2*clm);		/* 128 iters of 2*clm */ \
r5_x5cl_20_reals_first_fft_scratch_preload; \
for(unsigned int loopA = 16; loopA; loopA--) { \
	for(unsigned int loopB = 8*2*clm; loopB; ) { \
		xtouch(xptr(rcx+rbp+(clm-1)*128));	/* Preload TLB */ \
		do{ \
			r5_x5cl_20_reals_first_fft_scratch(rsi, 64, 16*clmblkdst8, rdi); \
		}while(--loopB & (2*clm-1)); \
		newxloop1doublepre(clm, b1b, b1a, clmblkdst, rdi, XMM_SCD9); /* 2*clm iters prefetching 4*clm cache lines */ \
	} \
	rsi += -8*clmblkdst+clmblkdst8; \
} \
rsi += -16*clmblkdst8;	/* Restore source pointer */ \
end_timer(22); \
\
/* Do FFT levels 6,7 */ \
/* */ \
/* On input the 64-byte cache lines hold these data values: */ \
/*	0K	128K	256K	384K */ \
/*	1K	... */ \
/*	... */ \
/*	127K	... */ \
/*	512K	... */ \
/*	... */ \
/* On output the 64-byte cache lines hold these data values: */ \
/*	0K	32K	64K	96K */ \
/*	1K	... */ \
/*	... */ \
/*	31K	... */ \
/*	128K	... */ \
/*	... */ \
\
/* Do 32 eight_reals_fft operations */ \
/* Do 32 four_complex_fft operations */ \
/*	distance between fft data elements is 32K */ \
/*	do 32 macros each processing 32 data values */ \
/* Do 256 four_complex_fft operations */ \
/*	distance between fft data elements is 32K */ \
/*	do 128 macros each processing 32 data values */ \
\
start_timer(24); \
rdx = (uintptr_t)g->premult_addr;/* Address of premultiplier/sin/cos table */ \
rdx += 8*2*clm*XMM_SCND2 + 8*2*clm*XMM_SCND2R + 2*clm*XMM_SCD4 + 10*2*clm*XMM_SCD8; /* Real data sin/cos ptr */ \
rdi = (rdx+32*2*clm*XMM_SCD2); /* Complex data sin/cos ptr */ \
loops_init(160*2*clm);		/* 5 iterations of 32 iters of 2*clm */ \
r4_x4cl_eight_reals_four_complex_djbfft_preload; \
for(unsigned int loopA = 4; loopA; loopA--) { \
	for(unsigned int loopB = 8*2*clm; loopB; ) { \
		xtouch(xptr(rcx+rbp+(clm-1)*128));	/* Preload TLB */ \
		do{ \
			r4_x4cl_eight_reals_four_complex_djbfft(rsi, 64, 4*clmblkdst8, 8*clmblkdst8, rdi, rdx); \
			rdi += XMM_SCD2;		/* Next sine/cosine pointer */ \
			rdx += XMM_SCD2;		/* Next sine/cosine pointer */ \
		}while(--loopB & (2*clm-1)); \
		newxloop1(clm, b4b, clmblkdst);	/* 2*clm iterations */ \
	} \
	rsi += -8*clmblkdst+clmblkdst8; \
} \
rsi += -4*clmblkdst8+16*clmblkdst8; \
rdi += -32*2*clm*XMM_SCD2;	/* Restore sin/cos ptr */ \
\
loops_reset; \
r4_x4cl_four_complex_djbfft_preload; \
for(unsigned int loopA = 4; loopA; loopA--) { \
	for(unsigned int loopB = 4; loopB; loopB--) { \
		for(unsigned int loopC = 8*2*clm; loopC; ) { \
			xtouch(xptr(rcx+rbp+(clm-1)*128));	/* Preload TLB */ \
			do{ \
				r4_x4cl_four_complex_djbfft(rsi, 64, 4*clmblkdst8, 8*clmblkdst8, rdi); \
				rdi += XMM_SCD2;		/* Next sine/cosine pointer */ \
			}while(--loopC & (2*clm-1)); \
			newxloop1(clm, b5b, clmblkdst);	/* 2*clm iterations */ \
		} \
		rsi += -8*clmblkdst+clmblkdst8; /* Next source pointer */ \
	} \
	rsi += -4*clmblkdst8+16*clmblkdst8; /* Next source pointer */ \
	rdi += -32*2*clm*XMM_SCD2;	/* Restore sin/cos ptr */ \
} \
rsi += -5*16*clmblkdst8;	/* Restore source pointer */ \
end_timer(24); \
\
/* Do FFT levels 8,9 */ \
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
/* Do 8 eight_reals_fft operations */ \
/* Do 8 four_complex_fft operations */ \
/*	distance between fft data elements is 8K */ \
/*	do 8 macros each processing 32 data values */ \
/* Do 308 four_complex_fft operations */ \
/*	distance between fft data elements is 8K */ \
/*	do 152 macros each processing 32 data values */ \
\
start_timer(25); \
rdi += -32*2*clm*XMM_SCD2-8*2*clm*XMM_SCND2; /* Complex data sin/cos ptr */ \
rdx = (rdi-8*2*clm*XMM_SCND2R); /* Real data sin/cos ptr */ \
loops_init(160*2*clm);		/* 20 iterations of 8 iters of 2*clm */ \
r4_x4cl_wpn_eight_reals_four_complex_djbfft_preload; \
for(unsigned int loopA = 8*2*clm; loopA; ) { \
	xtouch(xptr(rcx+rbp+(clm-1)*128));	/* Preload TLB */ \
	do{ \
		r4_x4cl_wpn_eight_reals_four_complex_djbfft(rsi, 64, clmblkdst8, 2*clmblkdst8, rdi, rdx); \
		rdi += XMM_SCND2;		/* Next sine/cosine pointer */ \
		rdx += XMM_SCND2R;		/* Next sine/cosine pointer */ \
	}while(--loopA & (2*clm-1)); \
	newxloop1(clm, b6b, clmblkdst);	/* 2*clm iterations */ \
} \
rsi += -8*clmblkdst+4*clmblkdst8; \
rdi += -8*2*clm*XMM_SCND2;	/* Restore sin/cos ptr */ \
\
loops_reset; \
r4_x4cl_wpn_four_complex_djbfft_preload; \
for(unsigned int loopA = 19; loopA; loopA--) { \
	for(unsigned int loopB = 8*2*clm; loopB; ) { \
		xtouch(xptr(rcx+rbp+(clm-1)*128));	/* Preload TLB */ \
		do{ \
			r4_x4cl_wpn_four_complex_djbfft(rsi, 64, clmblkdst8, 2*clmblkdst8, rdi); \
			rdi += XMM_SCND2;		/* Next sine/cosine pointer */ \
		}while(--loopB & (2*clm-1)); \
		newxloop1(clm, b7b, clmblkdst);	/* 2*clm iterations */ \
	} \
	rsi += -8*clmblkdst+4*clmblkdst8; /* Next source pointer */ \
	rdi += -8*2*clm*XMM_SCND2;	/* Restore sin/cos ptr */ \
} \
rsi += -20*4*clmblkdst8;	/* Restore source pointer */ \
end_timer(25); \
\
/* Do FFT levels 10,11,12 */ \
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
/* Do 158 eight_complex_fft8 operations */ \
/*	distance between fft data elements is 1K */ \
/*	do 79 macros each processing 64 data values */ \
\
start_timer(26); \
rdi += -8*2*clm*XMM_SCND2R-10*2*clm*XMM_SCD8+128; /* Complex data premult/sin/cos ptr and use -128 */ \
;					/* addressing to generate 1 byte offsets in radix-8 macro */ \
rdx = (rdi-128-2*clm*XMM_SCD4); /* Real data sin/cos ptr */ \
rcx = (uintptr_t)g->data_addr;/* Destination pointer */ \
loops_init_prefetch(80*2*clm, 64, 1, rbp); /* 80 iters of 2*clm */ \
pfing(rbp = (uintptr_t)g->premult_prefetch);/* Load premult/sin/cos prefetch pointer */ \
r8_sg8cl_sixteen_reals_eight_complex_fft8_preload; \
for(unsigned int loopA = 2*clm; loopA; loopA--) { \
	r8_sg8cl_sixteen_reals_eight_complex_fft8(rsi, 64, clmblkdst, 2*clmblkdst, 4*clmblkdst, rcx, 64, blkdst, 2*blkdst, 4*blkdst, rdi-128, rdx); \
	rdi += XMM_SCD8;		/* Next sine/cosine pointer */ \
	rdx += XMM_SCD4;		/* Next sine/cosine pointer */ \
} \
rsi += -2*clm*64+clmblkdst8; /* Next source pointer */ \
rcx += -2*clm*64+8*blkdst;	/* Next dest pointer */ \
rdi += -2*clm*XMM_SCD8;	/* Restore sin/cos pointer */ \
\
loops_reset; \
for(unsigned int loopA = 10*8-1; loopA; ) { \
	do{ \
		for(unsigned int loopC = 2*clm; loopC; loopC--) { \
			r8_sg8cl_eight_complex_fft8(rsi, 64, clmblkdst, 2*clmblkdst, 4*clmblkdst, rcx, 64, blkdst, 2*blkdst, 4*blkdst, rdi-128); \
			rdi += XMM_SCD8;		/* Next premult/sin/cos ptr */ \
		} \
		rsi += -2*clm*64+clmblkdst8; /* Next source pointer */ \
		rcx += -2*clm*64+8*blkdst;	/* Next dest pointer */ \
		rdi += -2*clm*XMM_SCD8;	/* Restore premult/sine/cosine pointer */ \
	}while(--loopA & 7); \
	rdi += 2*clm*XMM_SCD8;	/* Next premult/sine/cosine pointer */ \
} \
end_timer(26); \
\
/* Work on next 4*clm sets of 2560 values */ \
\
endlp:	xpass1_get_next_block(pass2, c0b, b0b, _wpn); /* Get next block and jump */ \
/* done: */ xfft_footer_rn; \
\
\
\
/* Perform a pass 1 FFT of length 3072 doubles. */ \
/* Memory examples are for a 3072K FFT with 10 levels done in pass 2. */ \
/* Sin/cos data needed: 32*2*clm*XMM_SCD1 + 32*2*clm*XMM_SCD1 + 8*2*clm*XMM_SCND2 + 8*2*clm*XMM_SCND2R + */ \
/*			2*clm*XMM_SCD4 + 16*2*clm*XMM_SCD8 = 126*2*clm 64-byte cache lines */ \
\

#define r4dwpn_pass1sc3072(procname, xpass2name, pass2_levels, clm) /* (fftname, pass2_levels, clmarg) */ \
\
/* set_FFT_constants(fftname, 0, 0, pass2_levels, clmarg); */ \
xfft_header_rn(procname, xpass2name, pass2_levels, clm); \
\
/* Do pass 2 (the last 10 FFT levels if this is a 3072K FFT) */ \
\
/* pass2: */start_timer(0); \
xpass2name(g); \
end_timer(0); \
if(g->ffttype == 1)		/* We're done if FFTing only */ \
	return; /* goto done; */ \
\
/* Do the inverse FFT */ \
\
pass1_inverse_fft_setup(c0b, procname); \
\
/* Do inverse FFT levels 10,11,12 */ \
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
/* Do 191 eight_complex_unfft8 operations */ \
/*	distance between fft data elements is 1K */ \
/*	do 191 macros each processing 32 data values */ \
\
c0b:	start_timer(16); \
rsi = (uintptr_t)g->data_addr;/* Get source address */ \
rcx = (uintptr_t)g->scratch_area;	/* Get address of scratch area */ \
rdx = (uintptr_t)g->premult_addr;/* Real data sin/cos pointer */ \
rdi = (rdx+2*clm*XMM_SCD4+128); /* Complex data sin/cos ptr and use -128 */ \
;					/* addressing to generate 1 byte offsets in radix-8 macro */ \
pfing(rbp = (rdx+(126-48)*2*clm*64)); /* Finish prefetching the sin/cos data */ \
loops_init_prefetch(192*2*clm, 64, 4, rbp); /* 192*2*clm iterations */ \
/*bug	xtouchpages rbp, 192*2*clm * 64 / 4 ;; Prefetch TLBs */ \
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
for(unsigned int loopA = 16*6*2-1; loopA; ) { \
	for(unsigned int loopB = 6; loopB; loopB--) { \
		do{ \
			for(unsigned int loopD = 2*clm; loopD; loopD--) { \
				r8_sg4cl_eight_complex_unfft8(rsi, 64, blkdst, 2*blkdst, rcx, 64, clmblkdst, 2*clmblkdst, rdi-128); \
				rdi += XMM_SCD8;		/* Next premult/sin/cos data */ \
			} \
			rsi += -2*clm*64+4*blkdst;	/* Next source pointer */ \
			rcx += -2*clm*64+4*clmblkdst; /* Next dest pointer */ \
			rdi += -2*clm*XMM_SCD8;	/* Restore premult/sin/cos data pointer */ \
		}while(--loopA & 1); \
		rcx += -2*4*clmblkdst+clmblkdst8; \
	} \
	rdi += 2*clm*XMM_SCD8;	/* Next premult/sin/cos pointer */ \
} \
rsi = (uintptr_t)g->scratch_area;	/* Restore source pointer */ \
rdx = (rdi-128);		/* Next real data sin/cos pointer */ \
end_timer(16); \
\
/* Do inverse FFT levels 8,9 */ \
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
/* Do 376 four_complex_unfft operations */ \
/*	distance between fft data elements is 8K */ \
/*	do 188 macros each processing 32 data values */ \
\
start_timer(17); \
rdi = (rdx+8*2*clm*XMM_SCND2R); /* Complex data sin/cos ptr */ \
loops_init(192*2*clm);		/* 192 iterations of 2*clm */ \
\
r4_x4cl_wpn_eight_reals_unfft_preload; \
for(unsigned int loopA = 4; loopA; loopA--) { \
	for(unsigned int loopB = 2*clm; loopB; loopB--) { \
		r4_x4cl_wpn_eight_reals_unfft(rsi, 64, 4*clmblkdst, clmblkdst8, rdi, 4*2*clm*XMM_SCND2, rdx, 4*2*clm*XMM_SCND2R); \
		rdi += XMM_SCND2;		/* Next sine/cosine pointer */ \
		rdx += XMM_SCND2R;		/* Next sine/cosine pointer */ \
	} \
	rsi += -2*clm*64+clmblkdst; /* Next source pointer */ \
} \
rsi += -4*clmblkdst+2*clmblkdst8; /* Next source pointer */ \
rdi += -4*2*clm*XMM_SCND2;	/* Restore sine/cosine pointer */ \
\
loops_reset; \
r4_x4cl_wpn_four_complex_djbunfft_preload; \
for(unsigned int loopA = 47; loopA; loopA--) { \
	for(unsigned int loopB = 4; loopB; loopB--) { \
		for(unsigned int loopC = 2*clm; loopC; loopC--) { \
			r4_x4cl_wpn_four_complex_djbunfft(rsi, 64, 4*clmblkdst, clmblkdst8, rdi, 4*2*clm*XMM_SCND2); \
			rdi += XMM_SCND2;		/* Next clm's sin/cos data */ \
		} \
		rsi += -2*clm*64+clmblkdst; /* Next source pointer */ \
	} \
	rsi += -4*clmblkdst+2*clmblkdst8; /* Next source pointer */ \
	rdi += -4*2*clm*XMM_SCND2;	/* Restore sine/cosine pointer */ \
} \
rsi += -48*2*clmblkdst8;	/* Restore source pointer */ \
rdx = (rdi+8*2*clm*XMM_SCND2); /* Next real data sin/cos ptr */ \
end_timer(17); \
\
/* Do inverse FFT levels 6,7 (actually 1.585 levels) */ \
/* On input the 64-byte cache lines hold these data values: */ \
/*	0K	16K	32K	48K */ \
/*	1K	... */ \
/*	... */ \
/*	15K */ \
/*	64K	... */ \
/*	... */ \
/* On output the 64-byte cache lines hold these data values: */ \
/*	0K	48K	96K	144K */ \
/*	1K	... */ \
/*	... */ \
/*	47K */ \
/*	192K	... */ \
/*	... */ \
\
/* Do 32 six_reals_unfft operations */ \
/*	distance between fft data elements is 32K */ \
/*	do 16 macros each processing 32 data values */ \
/* Do 480 three_complex_unfft operations */ \
/*	distance between fft data elements is 32K */ \
/*	do 240 macros each processing 32 data values */ \
\
start_timer(18); \
rdi = (rdx+32*2*clm*XMM_SCD1); /* Complex data sin/cos ptr */ \
loops_init(256*2*clm);		/* 256 iterations of 2*clm */ \
\
for(unsigned int loopA = 2; loopA; loopA--) { \
	for(unsigned int loopB = 8; loopB; loopB--) { \
		for(unsigned int loopC = 2*clm; loopC; loopC--) { \
			r3_x3cl_six_reals_unfft(rsi, 64, 2*clmblkdst8, rdi, 16*2*clm*XMM_SCD1, rdx, 16*2*clm*XMM_SCD1); \
			rdi += XMM_SCD1;		/* Next sine/cosine pointer */ \
			rdx += XMM_SCD1;		/* Next sine/cosine pointer */ \
		} \
		rsi += -2*clm*64+clmblkdst; /* Next source pointer */ \
	} \
	rsi += -8*clmblkdst+clmblkdst8; /* Next source pointer */ \
} \
rsi += -2*clmblkdst8+6*clmblkdst8; /* Next source pointer */ \
rdi += -16*2*clm*XMM_SCD1; /* Restore sine/cosine pointer */ \
\
loops_reset; \
r3_x3cl_three_complex_djbunfft_preload; \
for(unsigned int loopA = 15; loopA; loopA--) { \
	for(unsigned int loopB = 2; loopB; loopB--) { \
		for(unsigned int loopC = 8; loopC; loopC--) { \
			for(unsigned int loopD = 2*clm; loopD; loopD--) { \
				r3_x3cl_three_complex_djbunfft(rsi, 64, 2*clmblkdst8, rdi, 16*2*clm*XMM_SCD1); \
				rdi += XMM_SCD1;		/* Next clm's sin/cos data */ \
			} \
			rsi += -2*clm*64+clmblkdst; /* Next source pointer */ \
		} \
		rsi += -8*clmblkdst+clmblkdst8; /* Next source pointer */ \
	} \
	rsi += -2*clmblkdst8+6*clmblkdst8; /* Next source pointer */ \
	rdi += -16*2*clm*XMM_SCD1; /* Restore sine/cosine pointer */ \
} \
rsi += -16*6*clmblkdst8;	/* Restore source pointer */ \
end_timer(18); \
\
/* Do inverse FFT levels 4,5 */ \
/* On input the 64-byte cache lines hold these data values: */ \
/*	0K	48K	96K	144K */ \
/*	1K	... */ \
/*	... */ \
/*	47K */ \
/*	192K	... */ \
/*	... */ \
/* On output the 64-byte cache lines hold these data values: */ \
/*	0K	192K	384K	576K */ \
/*	1K	... */ \
/*	... */ \
/*	191K */ \
/*	768K	... */ \
/*	... */ \
\
/* Do 96 eight_reals_unfft operations */ \
/*	distance between fft data elements is 96K */ \
/*	do 48 macros each processing 32 data values */ \
/* Do 288 four_complex_unfft operations */ \
/*	distance between fft data elements is 96K */ \
/*	do 144 macros each processing 32 data values */ \
\
start_timer(19); \
rdx = (uintptr_t)g->sincos2;/* Load fixed sin/cos data ptr */ \
rdx += 384*XMM_SCD3;	/* Real data sin/cos ptr */ \
rdi = (rdx+96*XMM_SCD2);	/* Complex data sin/cos ptr */ \
loops_init(192*2*clm);		/* 192 iterations of 2*clm */ \
\
r4_x4cl_eight_reals_unfft_preload; \
for(unsigned int loopA = 6; loopA; loopA--) { \
	for(unsigned int loopB = 8; loopB; loopB--) { \
		for(unsigned int loopC = 2*clm; loopC; loopC--) { \
			r4_x4cl_eight_reals_unfft(rsi, 64, 6*clmblkdst8, 12*clmblkdst8, rdi, 48*XMM_SCD2, rdx, 48*XMM_SCD2); \
		} \
		rsi += -2*clm*64+clmblkdst; /* Next source pointer */ \
		rdi += XMM_SCD2;		/* Next sine/cosine pointer */ \
		rdx += XMM_SCD2;		/* Next sine/cosine pointer */ \
	} \
	rsi += -8*clmblkdst+clmblkdst8; /* Next source pointer */ \
} \
rsi += -6*clmblkdst8+24*clmblkdst8; /* Next source pointer */ \
rdi += -48*XMM_SCD2;	/* Restore sine/cosine pointer */ \
\
loops_reset; \
r4_x4cl_four_complex_djbunfft_preload; \
for(unsigned int loopA = 3; loopA; loopA--) { \
	for(unsigned int loopB = 6; loopB; loopB--) { \
		for(unsigned int loopC = 8; loopC; loopC--) { \
			for(unsigned int loopD = 2*clm; loopD; loopD--) { \
				r4_x4cl_four_complex_djbunfft(rsi, 64, 6*clmblkdst8, 12*clmblkdst8, rdi, 48*XMM_SCD2); \
			} \
			rsi += -2*clm*64+clmblkdst; /* Next source pointer */ \
			rdi += XMM_SCD2;		/* Next clm's sin/cos data */ \
		} \
		rsi += -8*clmblkdst+clmblkdst8; /* Next source pointer */ \
	} \
	rsi += -6*clmblkdst8+24*clmblkdst8; /* Next source pointer */ \
	rdi += -48*XMM_SCD2;	/* Restore sine/cosine pointer */ \
} \
rsi += -4*24*clmblkdst8;	/* Restore source pointer */ \
end_timer(19); \
\
/* Do inverse FFT levels 1,2,3 */ \
/* On input the 64-byte cache lines hold these data values: */ \
/*	0K	192K	384K	576K */ \
/*	1K	... */ \
/*	... */ \
/*	191K */ \
/*	768K	... */ \
/*	... */ \
/* On output the 64-byte cache lines hold these data values: */ \
/*	0K	768K	1536K	2304K */ \
/*	1K	... */ \
/*	... */ \
/*	767K */ \
\
/* Do 384 eight_reals_last_unfft operations */ \
/*	distance between fft data elements is 384K */ \
/*	do 192 macros each processing 32 data values */ \
\
start_timer(20); \
rdi = (uintptr_t)g->sincos2;/* Load fixed sin/cos postmultipliers */ \
loops_init(192*2*clm);		/* 192 iterations of 2*clm */ \
r4_x4cl_eight_reals_last_unfft_preload; \
for(unsigned int loopA = 24; loopA; loopA--) { \
	for(unsigned int loopB = 8; loopB; loopB--) { \
		for(unsigned int loopC = 2*clm; loopC; loopC--) { \
			r4_x4cl_eight_reals_last_unfft(rsi, 64, 24*clmblkdst8, 48*clmblkdst8, rdi, 192*XMM_SCD3); \
		} \
		newxloop1nopre(clm, c1b, clmblkdst, rdi, XMM_SCD3); /* 2*clm iterations */ \
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
copy_scratch_data_to_fft(768, clm); \
goto	endlp;			/* Skip FFT code if flag not set */ \
\
/* Fire up auxiliary threads to do pass 1 forward FFT */ \
\
/* pass1: */pass1_forward_fft_setup(b0b); \
\
/* Copy data from the FFT data area back to the scratch area */ \
\
b0b:	copy_fft_data_to_scratch(768, clm); \
\
/* Do FFT levels 1,2,3 */ \
/* */ \
/* On input the 64-byte cache lines hold these data values: */ \
/*	0K	768K	1536K	2304K */ \
/*	1K	... */ \
/*	... */ \
/*	767K */ \
/* On output the 64-byte cache lines hold these data values: */ \
/*	0K	384K	768K	1152K */ \
/*	1K	... */ \
/*	... */ \
/*	383K	... */ \
/*	1536K	... */ \
/*	... */ \
\
/* Do 384 eight_reals_first_fft operations */ \
/*	distance between fft data elements is 384K */ \
/*	do 384 iterations each processing 16 data values */ \
\
yes_fft:start_timer(22); \
rsi = (uintptr_t)g->scratch_area;	/* Get address of scratch area */ \
rdi = (uintptr_t)g->sincos2;/* Load fixed sin/cos postmultipliers */ \
set_data_prefetch_ptrs;		/* Init prefetching pointers */ \
loops_init(384*2*clm);		/* 48 iters of 4 iters of 2*2*clm */ \
r4_x2cl_eight_reals_first_fft_scratch_preload; \
for(unsigned int loopA = 48; loopA; loopA--) { \
	for(unsigned int loopB = 4*2*2*clm; loopB; ) { \
		xtouch(xptr(rcx+rbp+(clm-1)*128));	/* Preload TLB */ \
		do{ \
			do{ \
				r4_x2cl_eight_reals_first_fft_scratch(rsi, 64, 48*clmblkdst8, rdi); \
			}while(--loopB & (2*clm-1)); \
			newxloop2(clm, b1b, clmblkdst, rdi, XMM_SCD3); /* 2*2*clm iters prefetching 2*clm cache lines */ \
		}while(loopB & (2*2*clm-1)); \
	} \
	rsi += -4*2*clmblkdst+clmblkdst8; \
} \
rsi += -48*clmblkdst8;	/* Restore source pointer */ \
rdx = rdi;		/* Next real data sin/cos ptr */ \
pfing(rcx += -8*blkdst);		/* We just prefetched 192 blocks, correct for */ \
;					/* code below that will prefetch 8 too many */ \
end_timer(22); \
\
/* Do FFT levels 4,5 */ \
/* */ \
/* On input the 64-byte cache lines hold these data values: */ \
/*	0K	384K	768K	1152K */ \
/*	1K	... */ \
/*	... */ \
/*	383K	... */ \
/*	1536K	... */ \
/*	... */ \
/* On output the 64-byte cache lines hold these data values: */ \
/*	0K	96K	192K	288K */ \
/*	1K	... */ \
/*	... */ \
/*	95K	... */ \
/*	384K	... */ \
/*	... */ \
\
/* Do 96 eight_reals_fft operations */ \
/* Do 96 four_complex_fft operations */ \
/*	distance between fft data elements is 96K */ \
/*	do 96 macros each processing 32 data values */ \
/* Do 192 four_complex_fft operations */ \
/*	distance between fft data elements is 96K */ \
/*	do 96 macros each processing 32 data values */ \
\
start_timer(24); \
rdi = (rdx+96*XMM_SCD2);	/* Complex data sin/cos ptr */ \
loops_init(96*2*clm);		/* 12 iterations of 8 iters of 2*clm */ \
r4_x4cl_eight_reals_four_complex_djbfft_preload; \
for(unsigned int loopA = 12; loopA; loopA--) { \
	for(unsigned int loopB = 4*2*2*clm; loopB; ) { \
		xtouch(xptr(rcx+rbp+(clm-1)*128));	/* Preload TLB */ \
		do{ \
			do{ \
				r4_x4cl_eight_reals_four_complex_djbfft(rsi, 64, 12*clmblkdst8, 24*clmblkdst8, rdi, rdx); \
			}while(--loopB & (2*clm-1)); \
			newxloop2(clm, b2b, clmblkdst, rdi, XMM_SCD2, rdx, XMM_SCD2); /* 2*2*clm iterations */ \
		}while(loopB & (2*2*clm-1)); \
	} \
	rsi += -8*clmblkdst+clmblkdst8; \
} \
rsi += -12*clmblkdst8+48*clmblkdst8; \
rdi += -96*XMM_SCD2;	/* Restore sin/cos ptr */ \
\
loops_init(96*2*clm);		/* 12 iterations of 8 iters of 2*clm */ \
r4_x4cl_four_complex_djbfft_preload; \
for(unsigned int loopA = 12; loopA; loopA--) { \
	for(unsigned int loopB = 8*2*clm; loopB; ) { \
		xtouch(xptr(rcx+rbp+(clm-1)*128));	/* Preload TLB */ \
		do{ \
			r4_x4cl_four_complex_djbfft(rsi, 64, 12*clmblkdst8, 24*clmblkdst8, rdi); \
		}while(--loopB & (2*clm-1)); \
		newxloop1(clm, b3b, clmblkdst, rdi, XMM_SCD2);	/* 2*clm iterations */ \
	} \
	rsi += -8*clmblkdst+clmblkdst8; /* Next source pointer */ \
} \
rsi += -12*clmblkdst8-48*clmblkdst8; /* Restore source pointer */ \
end_timer(24); \
\
/* Do FFT levels 6,7 (actually 1.585 levels) */ \
/* */ \
/* On input the 64-byte cache lines hold these data values: */ \
/*	0K	96K	192K	288K */ \
/*	1K	... */ \
/*	... */ \
/*	95K	... */ \
/*	384K	... */ \
/*	... */ \
/* On output the 64-byte cache lines hold these data values: */ \
/*	0K	32K	64K	96K */ \
/*	1K	... */ \
/*	... */ \
/*	31K	... */ \
/*	128K	... */ \
/*	... */ \
\
/* Do 32 six_reals_fft operations */ \
/* Do 32 three_complex_fft operations */ \
/*	distance between fft data elements is 32K */ \
/*	do 32 macros each processing 32 data values */ \
/* Do 448 three_complex_fft operations */ \
/*	distance between fft data elements is 32K */ \
/*	do 224 macros each processing 32 data values */ \
\
start_timer(23); \
rdx = (uintptr_t)g->premult_addr;/* Address of premultiplier/sin/cos table */ \
rdx += 8*2*clm*XMM_SCND2 + 8*2*clm*XMM_SCND2R + 2*clm*XMM_SCD4 + 16*2*clm*XMM_SCD8; /* Real data sin/cos ptr */ \
rdi = (rdx+32*2*clm*XMM_SCD1); /* Real data sin/cos ptr */ \
loops_init(256*2*clm);		/* 2 sets of 16 iters of 8 of 2*clm */ \
r3_x3cl_six_reals_three_complex_djbfft_preload; \
for(unsigned int loopA = 4; loopA; loopA--) { \
	for(unsigned int loopB = 8*2*clm; loopB; ) { \
		xtouch(xptr(rcx+rbp+(clm-1)*128));	/* Preload TLB */ \
		do{ \
			r3_x3cl_six_reals_three_complex_djbfft(rsi, 64, 4*clmblkdst8, rdi, rdx); \
			rdi += XMM_SCD1;		/* Next sine/cosine pointer */ \
			rdx += XMM_SCD1;		/* Next sine/cosine pointer */ \
		}while(--loopB & (2*clm-1)); \
		newxloop1(clm, b4b, clmblkdst);	/* 2*clm iterations */ \
	} \
	rsi += -8*clmblkdst+clmblkdst8; \
} \
rsi += -4*clmblkdst8+12*clmblkdst8; /* Next source pointer */ \
rdi += -32*2*clm*XMM_SCD1; /* Restore sin/cos ptr */ \
\
loops_reset; \
r3_x3cl_three_complex_djbfft_preload; \
for(unsigned int loopA = 7; loopA; loopA--) { \
	for(unsigned int loopB = 4; loopB; loopB--) { \
		for(unsigned int loopC = 8*2*clm; loopC; ) { \
			xtouch(xptr(rcx+rbp+(clm-1)*128));	/* Preload TLB */ \
			do{ \
				r3_x3cl_three_complex_djbfft(rsi, 64, 4*clmblkdst8, rdi); \
				rdi += XMM_SCD1;		/* Next sine/cosine pointer */ \
			}while(--loopC & (2*clm-1)); \
			newxloop1(clm, b5b, clmblkdst);	/* 2*clm iterations */ \
		} \
		rsi += -8*clmblkdst+clmblkdst8; \
	} \
	rsi += -4*clmblkdst8+12*clmblkdst8; /* Next source pointer */ \
	rdi += -32*2*clm*XMM_SCD1; /* Restore sin/cos ptr */ \
} \
rsi += -8*12*clmblkdst8;	/* Restore source pointer */ \
end_timer(23); \
\
/* Do FFT levels 8,9 */ \
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
/* Do 8 eight_reals_fft operations */ \
/* Do 8 four_complex_fft operations */ \
/*	distance between fft data elements is 8K */ \
/*	do 8 macros each processing 32 data values */ \
/* Do 368 four_complex_fft operations */ \
/*	distance between fft data elements is 8K */ \
/*	do 184 macros each processing 32 data values */ \
\
start_timer(25); \
rdi += -32*2*clm*XMM_SCD1-8*2*clm*XMM_SCND2; /* Complex data sin/cos ptr */ \
rdx = (rdi-8*2*clm*XMM_SCND2R); /* Real data sin/cos ptr */ \
loops_init(192*2*clm);		/* 24 iterations of 8 iters of 2*clm */ \
r4_x4cl_wpn_eight_reals_four_complex_djbfft_preload; \
for(unsigned int loopA = 8*2*clm; loopA; ) { \
	xtouch(xptr(rcx+rbp+(clm-1)*128));	/* Preload TLB */ \
	do{ \
		r4_x4cl_wpn_eight_reals_four_complex_djbfft(rsi, 64, clmblkdst8, 2*clmblkdst8, rdi, rdx); \
		rdi += XMM_SCND2;		/* Next sine/cosine pointer */ \
		rdx += XMM_SCND2R;		/* Next sine/cosine pointer */ \
	}while(--loopA & (2*clm-1)); \
	newxloop1nopre(clm, b6b, clmblkdst);	/* 2*clm iterations */ \
} \
rsi += -8*clmblkdst+4*clmblkdst8; \
rdi += -8*2*clm*XMM_SCND2;	/* Restore sin/cos ptr */ \
\
loops_reset; \
r4_x4cl_wpn_four_complex_djbfft_preload; \
for(unsigned int loopA = 23; loopA; loopA--) { \
	for(unsigned int loopB = 8*2*clm; loopB; ) { \
		xtouch(xptr(rcx+rbp+(clm-1)*128));	/* Preload TLB */ \
		do{ \
			r4_x4cl_wpn_four_complex_djbfft(rsi, 64, clmblkdst8, 2*clmblkdst8, rdi); \
			rdi += XMM_SCND2;		/* Next sine/cosine pointer */ \
		}while(--loopB & (2*clm-1)); \
		newxloop1(clm, b7b, clmblkdst);	/* 2*clm iterations */ \
	} \
	rsi += -8*clmblkdst+4*clmblkdst8; /* Next source pointer */ \
	rdi += -8*2*clm*XMM_SCND2;	/* Restore sin/cos ptr */ \
} \
rsi += -24*4*clmblkdst8;	/* Restore source pointer */ \
end_timer(25); \
\
/* Do FFT levels 10,11,12 */ \
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
/* Do 190 eight_complex_fft8 operations */ \
/*	distance between fft data elements is 1K */ \
/*	do 95 macros each processing 64 data values */ \
\
start_timer(26); \
rdi += -8*2*clm*XMM_SCND2R-16*2*clm*XMM_SCD8+128; /* Complex data sin/cos ptr and use -128 */ \
;					/* addressing to generate 1 byte offsets in radix-8 macro */ \
rdx = (rdi-128-2*clm*XMM_SCD4); /* Real data sin/cos ptr */ \
rcx = (uintptr_t)g->data_addr;/* Destination pointer */ \
loops_init_prefetch(96*2*clm, 64, 1, rbp); /* 96 iters of 2*clm */ \
pfing(rbp = (uintptr_t)g->premult_prefetch);/* Load premult/sin/cos prefetch pointer */ \
r8_sg8cl_sixteen_reals_eight_complex_fft8_preload; \
for(unsigned int loopA = 2*clm; loopA; loopA--) { \
	r8_sg8cl_sixteen_reals_eight_complex_fft8(rsi, 64, clmblkdst, 2*clmblkdst, 4*clmblkdst, rcx, 64, blkdst, 2*blkdst, 4*blkdst, rdi-128, rdx); \
	rdi += XMM_SCD8;		/* Next sine/cosine pointer */ \
	rdx += XMM_SCD4;		/* Next sine/cosine pointer */ \
} \
rsi += -2*clm*64+clmblkdst8; /* Next source pointer */ \
rcx += -2*clm*64+8*blkdst;	/* Next dest pointer */ \
rdi += -2*clm*XMM_SCD8;	/* Restore sin/cos pointer */ \
\
loops_reset; \
for(unsigned int loopA = 16*(6+2)-1; loopA; ) { \
	loopA -= 2; \
	do{ /* loop 5 or 6 times */ \
		for(unsigned int loopC = 2*clm; loopC; loopC--) { \
			r8_sg8cl_eight_complex_fft8(rsi, 64, clmblkdst, 2*clmblkdst, 4*clmblkdst, rcx, 64, blkdst, 2*blkdst, 4*blkdst, rdi-128); \
			rdi += XMM_SCD8;		/* Next premult/sin/cos ptr */ \
		} \
		rsi += -2*clm*64+clmblkdst8; /* Next source pointer */ \
		rcx += -2*clm*64+8*blkdst;	/* Next dest pointer */ \
		rdi += -2*clm*XMM_SCD8;	/* Restore premult/sine/cosine pointer */ \
	}while(--loopA & 7); \
	rdi += 2*clm*XMM_SCD8;	/* Next premult/sine/cosine pointer */ \
} \
end_timer(26); \
\
/* Work on next 4*clm sets of 3072 values */ \
\
endlp:	xpass1_get_next_block(pass2, c0b, b0b, _wpn); /* Get next block and jump */ \
/* done: */ xfft_footer_rn; \
\
\
\
/* Perform a pass 1 FFT of length 3584 doubles. */ \
/* Memory examples are for a 3584K FFT with 10 levels done in pass 2. */ \
/* Sin/cos data needed: 32*2*clm*XMM_SCD2 + 32*2*clm*XMM_SCD2 + 8*2*clm*XMM_SCND2 + 8*2*clm*XMM_SCND2R + */ \
/*			2*clm*XMM_SCD4 + 14*2*clm*XMM_SCD8 = 150*2*clm 64-byte cache lines */ \
\

#define r4dwpn_pass1sc3584(procname, xpass2name, pass2_levels, clm) /* (fftname, pass2_levels, clmarg) */ \
\
/* set_FFT_constants(fftname, 0, 0, pass2_levels, clmarg); */ \
xfft_header_rn(procname, xpass2name, pass2_levels, clm); \
\
/* Do pass 2 (the last 10 FFT levels if this is a 3584K FFT) */ \
\
/* pass2: */start_timer(0); \
xpass2name(g); \
end_timer(0); \
if(g->ffttype == 1)		/* We're done if FFTing only */ \
	return; /* goto done; */ \
\
/* Do the inverse FFT */ \
\
pass1_inverse_fft_setup(c0b, procname); \
\
/* Do inverse FFT levels 10,11,12 */ \
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
/* Do 223 eight_complex_unfft8 operations */ \
/*	distance between fft data elements is 1K */ \
/*	do 223 macros each processing 32 data values */ \
\
c0b:	start_timer(16); \
rsi = (uintptr_t)g->data_addr;/* Get source address */ \
rcx = (uintptr_t)g->scratch_area;	/* Get address of scratch area */ \
rdx = (uintptr_t)g->premult_addr;/* Real data sin/cos pointer */ \
rdi = (rdx+2*clm*XMM_SCD4+128); /* Complex data sin/cos ptr and use -128 */ \
;					/* addressing to generate 1 byte offsets in radix-8 macro */ \
pfing(rbp = (rdx+(150-56)*2*clm*64)); /* Finish prefetching the sin/cos data */ \
loops_init_prefetch(224*2*clm, 64, 4, rbp); /* 224*2*clm iterations */ \
/*bug	xtouchpages rbp, 224*2*clm * 64 / 4 ;; Prefetch TLBs */ \
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
for(unsigned int loopA = 14*8*2-1; loopA; ) { \
	do{ \
		do{ \
			for(unsigned int loopD = 2*clm; loopD; loopD--) { \
				r8_sg4cl_eight_complex_unfft8(rsi, 64, blkdst, 2*blkdst, rcx, 64, clmblkdst, 2*clmblkdst, rdi-128); \
				rdi += XMM_SCD8;		/* Next premult/sin/cos data */ \
			} \
			rsi += -2*clm*64+4*blkdst;	/* Next source pointer */ \
			rcx += -2*clm*64+4*clmblkdst; /* Next dest pointer */ \
			rdi += -2*clm*XMM_SCD8;	/* Restore premult/sin/cos data pointer */ \
		}while(--loopA & 1); \
		rcx += -2*4*clmblkdst+clmblkdst8; \
	}while(loopA & (8*2-1)); \
	rdi += 2*clm*XMM_SCD8;	/* Next premult/sin/cos pointer */ \
} \
rsi = (uintptr_t)g->scratch_area;	/* Restore source pointer */ \
rdx = (rdi-128);		/* Next real data sin/cos pointer */ \
end_timer(16); \
\
/* Do inverse FFT levels 8,9 */ \
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
/* Do 440 four_complex_unfft operations */ \
/*	distance between fft data elements is 8K */ \
/*	do 220 macros each processing 32 data values */ \
\
start_timer(17); \
rdi = (rdx+8*2*clm*XMM_SCND2R); /* Complex data sin/cos ptr */ \
loops_init(224*2*clm);		/* 224 iterations of 2*clm */ \
r4_x4cl_wpn_eight_reals_unfft_preload; \
for(unsigned int loopA = 4; loopA; loopA--) { \
	for(unsigned int loopB = 2*clm; loopB; loopB--) { \
		r4_x4cl_wpn_eight_reals_unfft(rsi, 64, 4*clmblkdst, clmblkdst8, rdi, 4*2*clm*XMM_SCND2, rdx, 4*2*clm*XMM_SCND2R); \
		rdi += XMM_SCND2;		/* Next sine/cosine pointer */ \
		rdx += XMM_SCND2R;		/* Next sine/cosine pointer */ \
	} \
	rsi += -2*clm*64+clmblkdst; /* Next source pointer */ \
} \
rsi += -4*clmblkdst+2*clmblkdst8; /* Next source pointer */ \
rdi += -4*2*clm*XMM_SCND2;	/* Restore sine/cosine pointer */ \
\
loops_reset; \
r4_x4cl_wpn_four_complex_djbunfft_preload; \
for(unsigned int loopA = 55; loopA; loopA--) { \
	for(unsigned int loopB = 4; loopB; loopB--) { \
		for(unsigned int loopC = 2*clm; loopC; loopC--) { \
			r4_x4cl_wpn_four_complex_djbunfft(rsi, 64, 4*clmblkdst, clmblkdst8, rdi, 4*2*clm*XMM_SCND2); \
			rdi += XMM_SCND2;		/* Next clm's sin/cos data */ \
		} \
		rsi += -2*clm*64+clmblkdst; /* Next source pointer */ \
	} \
	rsi += -4*clmblkdst+2*clmblkdst8; /* Next source pointer */ \
	rdi += -4*2*clm*XMM_SCND2;	/* Restore sine/cosine pointer */ \
} \
rsi += -56*2*clmblkdst8;	/* Restore source pointer */ \
rdx = (rdi+8*2*clm*XMM_SCND2); /* Next real data sin/cos ptr */ \
end_timer(17); \
\
/* Do inverse FFT levels 6,7 */ \
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
/*	256K	... */ \
/*	... */ \
\
/* Do 32 eight_reals_unfft operations */ \
/*	distance between fft data elements is 32K */ \
/*	do 16 macros each processing 32 data values */ \
/* Do 416 four_complex_unfft operations */ \
/*	distance between fft data elements is 32K */ \
/*	do 208 macros each processing 32 data values */ \
\
start_timer(18); \
rdi = (rdx+32*2*clm*XMM_SCD2); /* Complex data sin/cos ptr */ \
loops_init(224*2*clm);		/* 224 iterations of 2*clm */ \
r4_x4cl_eight_reals_unfft_preload; \
for(unsigned int loopA = 2; loopA; loopA--) { \
	for(unsigned int loopB = 8; loopB; loopB--) { \
		for(unsigned int loopC = 2*clm; loopC; loopC--) { \
			r4_x4cl_eight_reals_unfft(rsi, 64, 2*clmblkdst8, 4*clmblkdst8, rdi, 16*2*clm*XMM_SCD2, rdx, 16*2*clm*XMM_SCD2); \
			rdi += XMM_SCD2;		/* Next sine/cosine pointer */ \
			rdx += XMM_SCD2;		/* Next sine/cosine pointer */ \
		} \
		rsi += -2*clm*64+clmblkdst; /* Next source pointer */ \
	} \
	rsi += -8*clmblkdst+clmblkdst8; /* Next source pointer */ \
} \
rsi += -2*clmblkdst8+8*clmblkdst8; /* Next source pointer */ \
rdi += -16*2*clm*XMM_SCD2;	/* Restore sine/cosine pointer */ \
\
loops_reset; \
r4_x4cl_four_complex_djbunfft_preload; \
for(unsigned int loopA = 13; loopA; loopA--) { \
	for(unsigned int loopB = 2; loopB; loopB--) { \
		for(unsigned int loopC = 8; loopC; loopC--) { \
			for(unsigned int loopD = 2*clm; loopD; loopD--) { \
				r4_x4cl_four_complex_djbunfft(rsi, 64, 2*clmblkdst8, 4*clmblkdst8, rdi, 16*2*clm*XMM_SCD2); \
				rdi += XMM_SCD2;		/* Next clm's sin/cos data */ \
			} \
			rsi += -2*clm*64+clmblkdst; /* Next source pointer */ \
		} \
		rsi += -8*clmblkdst+clmblkdst8; /* Next source pointer */ \
	} \
	rsi += -2*clmblkdst8+8*clmblkdst8; /* Next source pointer */ \
	rdi += -16*2*clm*XMM_SCD2;	/* Restore sine/cosine pointer */ \
} \
rsi += -14*8*clmblkdst8;	/* Restore source pointer */ \
end_timer(18); \
\
/* Do inverse FFT levels 1,2,3,4,5 (actually 4.807 levels) */ \
/* On input the 64-byte cache lines hold these data values: */ \
/*	0K	64K	128K	192K */ \
/*	1K	... */ \
/*	... */ \
/*	63K */ \
/*	256K	... */ \
/*	... */ \
/* On output the 64-byte cache lines hold these data values: */ \
/*	0K	896K	1792K	2688K */ \
/*	1K	... */ \
/*	... */ \
/*	895K */ \
\
/* Do 512 28_reals_last_unfft operations */ \
/*	distance between fft data elements is 128K */ \
/*	do 64 macros each processing 112 data values */ \
\
start_timer(20); \
rdi = (uintptr_t)g->sincos2;/* Load fixed sin/cos postmultipliers */ \
set_data_prefetch_ptrs;		/* Init prefetching pointers */ \
loops_init(64*2*clm);		/* 64iterations of 2*clm */ \
r7_x14cl_28_reals_last_unfft_preload; \
for(unsigned int loopA = 8; loopA; loopA--) { \
	for(unsigned int loopB = 8*2*clm; loopB; ) { \
		xtouch(xptr(rcx+rbp+(clm-1)*128));	/* Preload TLB */ \
		do{ \
			r7_x14cl_28_reals_last_unfft(rsi, 64, 8*clmblkdst8, rdi, 64*XMM_SCD13); \
		}while(--loopB & (2*clm-1)); \
		newxloop1quadpre(clm, c1b, c1a, clmblkdst, rdi, XMM_SCD13); /* 2*clm iterations */ \
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
copy_scratch_data_to_fft(896, clm); \
goto	endlp;			/* Skip FFT code if flag not set */ \
\
/* Fire up auxiliary threads to do pass 1 forward FFT */ \
\
/* pass1: */pass1_forward_fft_setup(b0b); \
\
/* Copy data from the FFT data area back to the scratch area */ \
\
b0b:	copy_fft_data_to_scratch(896, clm); \
\
/* Do FFT levels 1,2,3,4,5 (actually 4.807 levels) */ \
/* */ \
/* On input the 64-byte cache lines hold these data values: */ \
/*	0K	896K	1792K	2688K */ \
/*	1K	... */ \
/*	... */ \
/*	895K */ \
/* On output the 64-byte cache lines hold these data values: */ \
/*	0K	128K	256K	384K */ \
/*	1K	... */ \
/*	... */ \
/*	127K	... */ \
/*	512K	... */ \
/*	... */ \
\
/* Do 512 28_reals_first_fft operations */ \
/*	distance between fft data elements is 128K */ \
/*	do 128 iterations each processing 56 data values */ \
\
yes_fft:start_timer(22); \
rsi = (uintptr_t)g->scratch_area;	/* Get address of scratch area */ \
rdi = (uintptr_t)g->sincos2;/* Load fixed sin/cos postmultipliers */ \
set_data_prefetch_ptrs;		/* Init prefetching pointers */ \
pfing(rcx += (256-24)*blkdst);	/* We've already prefetched 256 blocks, but */ \
;					/* code below will fetch 24 blocks too many */ \
loops_init(128*2*clm);		/* 128 iters of 2*clm */ \
r7_x7cl_28_reals_first_fft_scratch_preload; \
for(unsigned int loopA = 16; loopA; loopA--) { \
	for(unsigned int loopB = 8*2*clm; loopB; ) { \
		xtouch(xptr(rcx+rbp+(clm-1)*128));	/* Preload TLB */ \
		do{ \
			r7_x7cl_28_reals_first_fft_scratch(rsi, 64, 16*clmblkdst8, rdi); \
		}while(--loopB & (2*clm-1)); \
		newxloop1doublepre(clm, b1b, b1a, clmblkdst, rdi, XMM_SCD13); /* 2*clm iters prefetching 4*clm cache lines */ \
	} \
	rsi += -8*clmblkdst+clmblkdst8; \
} \
rsi += -16*clmblkdst8;	/* Restore source pointer */ \
end_timer(22); \
\
/* Do FFT levels 6,7 */ \
/* */ \
/* On input the 64-byte cache lines hold these data values: */ \
/*	0K	128K	256K	384K */ \
/*	1K	... */ \
/*	... */ \
/*	127K	... */ \
/*	512K	... */ \
/*	... */ \
/* On output the 64-byte cache lines hold these data values: */ \
/*	0K	32K	64K	96K */ \
/*	1K	... */ \
/*	... */ \
/*	31K	... */ \
/*	128K	... */ \
/*	... */ \
\
/* Do 32 eight_reals_fft operations */ \
/* Do 32 four_complex_fft operations */ \
/*	distance between fft data elements is 32K */ \
/*	do 32 macros each processing 32 data values */ \
/* Do 384 four_complex_fft operations */ \
/*	distance between fft data elements is 32K */ \
/*	do 192 macros each processing 32 data values */ \
\
start_timer(24); \
rdx = (uintptr_t)g->premult_addr;/* Address of premultiplier/sin/cos table */ \
rdx += 8*2*clm*XMM_SCND2 + 8*2*clm*XMM_SCND2R + 2*clm*XMM_SCD4 + 14*2*clm*XMM_SCD8; /* Real data sin/cos ptr */ \
rdi = (rdx+32*2*clm*XMM_SCD2); /* Complex data sin/cos ptr */ \
loops_init(224*2*clm);		/* 7 iterations of 32 iters of 2*clm */ \
r4_x4cl_eight_reals_four_complex_djbfft_preload; \
for(unsigned int loopA = 4; loopA; loopA--) { \
	for(unsigned int loopB = 8*2*clm; loopB; ) { \
		xtouch(xptr(rcx+rbp+(clm-1)*128));	/* Preload TLB */ \
		do{ \
			r4_x4cl_eight_reals_four_complex_djbfft(rsi, 64, 4*clmblkdst8, 8*clmblkdst8, rdi, rdx); \
			rdi += XMM_SCD2;		/* Next sine/cosine pointer */ \
			rdx += XMM_SCD2;		/* Next sine/cosine pointer */ \
		}while(--loopB & (2*clm-1)); \
		newxloop1nopre(clm, b4b, clmblkdst);	/* 2*clm iterations */ \
	} \
	rsi += -8*clmblkdst+clmblkdst8; \
} \
rsi += -4*clmblkdst8+16*clmblkdst8; \
rdi += -32*2*clm*XMM_SCD2;	/* Restore sin/cos ptr */ \
\
loops_reset; \
r4_x4cl_four_complex_djbfft_preload; \
for(unsigned int loopA = 6; loopA; loopA--) { \
	for(unsigned int loopB = 4; loopB; loopB--) { \
		for(unsigned int loopC = 8*2*clm; loopC; ) { \
			xtouch(xptr(rcx+rbp+(clm-1)*128));	/* Preload TLB */ \
			do{ \
				r4_x4cl_four_complex_djbfft(rsi, 64, 4*clmblkdst8, 8*clmblkdst8, rdi); \
				rdi += XMM_SCD2;		/* Next sine/cosine pointer */ \
			}while(--loopC & (2*clm-1)); \
			newxloop1(clm, b5b, clmblkdst);	/* 2*clm iterations */ \
		} \
		rsi += -8*clmblkdst+clmblkdst8; /* Next source pointer */ \
	} \
	rsi += -4*clmblkdst8+16*clmblkdst8; /* Next source pointer */ \
	rdi += -32*2*clm*XMM_SCD2;	/* Restore sin/cos ptr */ \
} \
rsi += -7*16*clmblkdst8;	/* Restore source pointer */ \
end_timer(24); \
\
/* Do FFT levels 8,9 */ \
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
/* Do 8 eight_reals_fft operations */ \
/* Do 8 four_complex_fft operations */ \
/*	distance between fft data elements is 8K */ \
/*	do 8 macros each processing 32 data values */ \
/* Do 432 four_complex_fft operations */ \
/*	distance between fft data elements is 8K */ \
/*	do 216 macros each processing 32 data values */ \
\
start_timer(25); \
rdi += -32*2*clm*XMM_SCD2-8*2*clm*XMM_SCND2; /* Complex data sin/cos ptr */ \
rdx = (rdi-8*2*clm*XMM_SCND2R); /* Real data sin/cos ptr */ \
loops_init(224*2*clm);		/* 28 iterations of 8 iters of 2*clm */ \
r4_x4cl_wpn_eight_reals_four_complex_djbfft_preload; \
for(unsigned int loopA = 8*2*clm; loopA; ) { \
	xtouch(xptr(rcx+rbp+(clm-1)*128));	/* Preload TLB */ \
	do{ \
		r4_x4cl_wpn_eight_reals_four_complex_djbfft(rsi, 64, clmblkdst8, 2*clmblkdst8, rdi, rdx); \
		rdi += XMM_SCND2;		/* Next sine/cosine pointer */ \
		rdx += XMM_SCND2R;		/* Next sine/cosine pointer */ \
	}while(--loopA & (2*clm-1)); \
	newxloop1nopre(clm, b6b, clmblkdst);	/* 2*clm iterations */ \
} \
rsi += -8*clmblkdst+4*clmblkdst8; \
rdi += -8*2*clm*XMM_SCND2;	/* Restore sin/cos ptr */ \
\
loops_reset; \
r4_x4cl_wpn_four_complex_djbfft_preload; \
for(unsigned int loopA = 27; loopA; loopA--) { \
	for(unsigned int loopB = 8*2*clm; loopB; ) { \
		xtouch(xptr(rcx+rbp+(clm-1)*128));	/* Preload TLB */ \
		do{ \
			r4_x4cl_wpn_four_complex_djbfft(rsi, 64, clmblkdst8, 2*clmblkdst8, rdi); \
			rdi += XMM_SCND2;		/* Next sine/cosine pointer */ \
		}while(--loopB & (2*clm-1)); \
		newxloop1(clm, b7b, clmblkdst);	/* 2*clm iterations */ \
	} \
	rsi += -8*clmblkdst+4*clmblkdst8; /* Next source pointer */ \
	rdi += -8*2*clm*XMM_SCND2;	/* Restore sin/cos ptr */ \
} \
rsi += -28*4*clmblkdst8;	/* Restore source pointer */ \
end_timer(25); \
\
/* Do FFT levels 10,11,12 */ \
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
/* Do 222 eight_complex_fft8 operations */ \
/*	distance between fft data elements is 1K */ \
/*	do 111 macros each processing 64 data values */ \
\
start_timer(26); \
rdi += -8*2*clm*XMM_SCND2R-14*2*clm*XMM_SCD8+128; /* Complex data premult/sin/cos ptr and use -128 */ \
;					/* addressing to generate 1 byte offsets in radix-8 macro */ \
rdx = (rdi-128-2*clm*XMM_SCD4); /* Real data sin/cos ptr */ \
rcx = (uintptr_t)g->data_addr;/* Destination pointer */ \
loops_init_prefetch(112*2*clm, 64, 1, rbp); /* 112 iters of 2*clm */ \
pfing(rbp = (uintptr_t)g->premult_prefetch);/* Load premult/sin/cos prefetch pointer */ \
r8_sg8cl_sixteen_reals_eight_complex_fft8_preload; \
for(unsigned int loopA = 2*clm; loopA; loopA--) { \
	r8_sg8cl_sixteen_reals_eight_complex_fft8(rsi, 64, clmblkdst, 2*clmblkdst, 4*clmblkdst, rcx, 64, blkdst, 2*blkdst, 4*blkdst, rdi-128, rdx); \
	rdi += XMM_SCD8;		/* Next sine/cosine pointer */ \
	rdx += XMM_SCD4;		/* Next sine/cosine pointer */ \
} \
rsi += -2*clm*64+clmblkdst8; /* Next source pointer */ \
rcx += -2*clm*64+8*blkdst;	/* Next dest pointer */ \
rdi += -2*clm*XMM_SCD8;	/* Restore sin/cos pointer */ \
\
loops_reset; \
for(unsigned int loopA = 14*8-1; loopA; ) { \
	do{ \
		for(unsigned int loopC = 2*clm; loopC; loopC--) { \
			r8_sg8cl_eight_complex_fft8(rsi, 64, clmblkdst, 2*clmblkdst, 4*clmblkdst, rcx, 64, blkdst, 2*blkdst, 4*blkdst, rdi-128); \
			rdi += XMM_SCD8;		/* Next premult/sin/cos ptr */ \
		} \
		rsi += -2*clm*64+clmblkdst8; /* Next source pointer */ \
		rcx += -2*clm*64+8*blkdst;	/* Next dest pointer */ \
		rdi += -2*clm*XMM_SCD8;	/* Restore premult/sine/cosine pointer */ \
	}while(--loopA & 7); \
	rdi += 2*clm*XMM_SCD8;	/* Next premult/sine/cosine pointer */ \
} \
end_timer(26); \
\
/* Work on next 4*clm sets of 3584 values */ \
\
endlp:	xpass1_get_next_block(pass2, c0b, b0b, _wpn); /* Get next block and jump */ \
/* done: */ xfft_footer_rn; \
\
\
\
/* Perform a pass 1 FFT of length 4096 doubles. */ \
/* Memory examples are for a 4M FFT with 10 levels done in pass 2. */ \
/* Sin/cos data needed: 32*2*clm*XMM_SCD2 + 32*2*clm*XMM_SCD2 + 8*2*clm*XMM_SCND2 + 8*2*clm*XMM_SCND2R + */ \
/*			2*clm*XMM_SCD4 + 16*2*clm*XMM_SCD8 = 158*2*clm 64-byte cache lines */ \
\

#define r4dwpn_pass1sc4096(procname, xpass2name, pass2_levels, clm) /* (fftname, pass2_levels, clmarg) */ \
\
/* set_FFT_constants(fftname, 0, 0, pass2_levels, clmarg); */ \
xfft_header_rn(procname, xpass2name, pass2_levels, clm); \
\
/* Do pass 2 (the last 10 FFT levels if this is a 4M FFT) */ \
\
/* pass2: */start_timer(0); \
xpass2name(g); \
end_timer(0); \
if(g->ffttype == 1)		/* We're done if FFTing only */ \
	return; /* goto done; */ \
\
/* Do the inverse FFT */ \
\
pass1_inverse_fft_setup(c0b, procname); \
\
/* Do inverse FFT levels 10,11,12 */ \
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
/* Do 255 eight_complex_unfft8 operations */ \
/*	distance between fft data elements is 1K */ \
/*	do 255 macros each processing 32 data values */ \
\
c0b:	start_timer(16); \
rsi = (uintptr_t)g->data_addr;/* Get source address */ \
rcx = (uintptr_t)g->scratch_area;	/* Get address of scratch area */ \
rdx = (uintptr_t)g->premult_addr;/* Real data sin/cos pointer */ \
rdi = (rdx+2*clm*XMM_SCD4+128); /* Complex data sin/cos ptr and use -128 */ \
;					/* addressing to generate 1 byte offsets in radix-8 macro */ \
pfing(rbp = (rdx+(158-32)*2*clm*64)); /* Finish prefetching the sin/cos data */ \
loops_init_prefetch(256*2*clm, 64, 8, rbp); /* 256*2*clm iterations */ \
/*bug	xtouchpages rbp, 256*2*clm * 64 / 8 ;; Prefetch TLBs */ \
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
/* c9a: */ xtouch(xptr(rbp+4096));		/* Prefetch TLB */ \
for(unsigned int loopA = 16*8*2-1; loopA; ) { \
	do{ \
		do{ \
			for(unsigned int loopD = 2*clm; loopD; loopD--) { \
				r8_sg4cl_eight_complex_unfft8(rsi, 64, blkdst, 2*blkdst, rcx, 64, clmblkdst, 2*clmblkdst, rdi-128); \
				rdi += XMM_SCD8;		/* Next premult/sin/cos data */ \
			} \
			rsi += -2*clm*64+4*blkdst;	/* Next source pointer */ \
			rcx += -2*clm*64+4*clmblkdst; /* Next dest pointer */ \
			rdi += -2*clm*XMM_SCD8;	/* Restore premult/sin/cos data pointer */ \
		}while(--loopA & 1); \
		rcx += -2*4*clmblkdst+clmblkdst8; \
	}while(loopA & (8*2-1)); \
	rdi += 2*clm*XMM_SCD8;	/* Next premult/sin/cos pointer */ \
} \
rsi = (uintptr_t)g->scratch_area;	/* Restore source pointer */ \
rdx = (rdi-128);		/* Next real data sin/cos pointer */ \
end_timer(16); \
\
/* Do inverse FFT levels 8,9 */ \
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
/* Do 504 four_complex_unfft operations */ \
/*	distance between fft data elements is 8K */ \
/*	do 252 macros each processing 32 data values */ \
\
start_timer(17); \
rdi = (rdx+8*2*clm*XMM_SCND2R); /* Complex data sin/cos ptr */ \
loops_init(256*2*clm);		/* 256 iterations of 2*clm */ \
r4_x4cl_wpn_eight_reals_unfft_preload; \
for(unsigned int loopA = 4; loopA; loopA--) { \
	for(unsigned int loopB = 2*clm; loopB; loopB--) { \
		r4_x4cl_wpn_eight_reals_unfft(rsi, 64, 4*clmblkdst, clmblkdst8, rdi, 4*2*clm*XMM_SCND2, rdx, 4*2*clm*XMM_SCND2R); \
		rdi += XMM_SCND2;		/* Next sine/cosine pointer */ \
		rdx += XMM_SCND2R;		/* Next sine/cosine pointer */ \
	} \
	rsi += -2*clm*64+clmblkdst; /* Next source pointer */ \
} \
rsi += -4*clmblkdst+2*clmblkdst8; /* Next source pointer */ \
rdi += -4*2*clm*XMM_SCND2;	/* Restore sine/cosine pointer */ \
\
loops_reset; \
r4_x4cl_wpn_four_complex_djbunfft_preload; \
for(unsigned int loopA = 63; loopA; loopA--) { \
	for(unsigned int loopB = 4; loopB; loopB--) { \
		for(unsigned int loopC = 2*clm; loopC; loopC--) { \
			r4_x4cl_wpn_four_complex_djbunfft(rsi, 64, 4*clmblkdst, clmblkdst8, rdi, 4*2*clm*XMM_SCND2); \
			rdi += XMM_SCND2;		/* Next clm's sin/cos data */ \
		} \
		rsi += -2*clm*64+clmblkdst; /* Next source pointer */ \
	} \
	rsi += -4*clmblkdst+2*clmblkdst8; /* Next source pointer */ \
	rdi += -4*2*clm*XMM_SCND2;	/* Restore sine/cosine pointer */ \
} \
rsi += -64*2*clmblkdst8;	/* Restore source pointer */ \
rdx = (rdi+8*2*clm*XMM_SCND2); /* Next real data sin/cos ptr */ \
end_timer(17); \
\
/* Do inverse FFT levels 6,7 */ \
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
/*	256K	... */ \
/*	... */ \
\
/* Do 32 eight_reals_unfft operations */ \
/*	distance between fft data elements is 32K */ \
/*	do 16 macros each processing 32 data values */ \
/* Do 480 four_complex_unfft operations */ \
/*	distance between fft data elements is 32K */ \
/*	do 240 macros each processing 32 data values */ \
\
start_timer(18); \
rdi = (rdx+32*2*clm*XMM_SCD2); /* Complex data sin/cos ptr */ \
loops_init(256*2*clm);		/* 256 iterations of 2*clm */ \
\
r4_x4cl_eight_reals_unfft_preload; \
for(unsigned int loopA = 2; loopA; loopA--) { \
	for(unsigned int loopB = 8; loopB; loopB--) { \
		for(unsigned int loopC = 2*clm; loopC; loopC--) { \
			r4_x4cl_eight_reals_unfft(rsi, 64, 2*clmblkdst8, 4*clmblkdst8, rdi, 16*2*clm*XMM_SCD2, rdx, 16*2*clm*XMM_SCD2); \
			rdi += XMM_SCD2;		/* Next sine/cosine pointer */ \
			rdx += XMM_SCD2;		/* Next sine/cosine pointer */ \
		} \
		rsi += -2*clm*64+clmblkdst; /* Next source pointer */ \
	} \
	rsi += -8*clmblkdst+clmblkdst8; /* Next source pointer */ \
} \
rsi += -2*clmblkdst8+8*clmblkdst8; /* Next source pointer */ \
rdi += -16*2*clm*XMM_SCD2; /* Restore sine/cosine pointer */ \
\
loops_reset; \
r4_x4cl_four_complex_djbunfft_preload; \
for(unsigned int loopA = 15; loopA; loopA--) { \
	for(unsigned int loopB = 2; loopB; loopB--) { \
		for(unsigned int loopC = 8; loopC; loopC--) { \
			for(unsigned int loopD = 2*clm; loopD; loopD--) { \
				r4_x4cl_four_complex_djbunfft(rsi, 64, 2*clmblkdst8, 4*clmblkdst8, rdi, 16*2*clm*XMM_SCD2); \
				rdi += XMM_SCD2;		/* Next clm's sin/cos data */ \
			} \
			rsi += -2*clm*64+clmblkdst; /* Next source pointer */ \
		} \
		rsi += -8*clmblkdst+clmblkdst8; /* Next source pointer */ \
	} \
	rsi += -2*clmblkdst8+8*clmblkdst8; /* Next source pointer */ \
	rdi += -16*2*clm*XMM_SCD2; /* Restore sine/cosine pointer */ \
} \
rsi += -16*8*clmblkdst8;	/* Restore source pointer */ \
end_timer(18); \
\
/* Do inverse FFT levels 4,5 */ \
/* On input the 64-byte cache lines hold these data values: */ \
/*	0K	64K	128K	192K */ \
/*	1K	... */ \
/*	... */ \
/*	63K */ \
/*	256K	... */ \
/*	... */ \
/* On output the 64-byte cache lines hold these data values: */ \
/*	0K	256K	512K	768K */ \
/*	1K	... */ \
/*	... */ \
/*	255K */ \
/*	1024K	... */ \
/*	... */ \
\
/* Do 128 eight_reals_unfft operations */ \
/*	distance between fft data elements is 128K */ \
/*	do 64 macros each processing 32 data values */ \
/* Do 384 four_complex_unfft operations */ \
/*	distance between fft data elements is 128K */ \
/*	do 192 macros each processing 32 data values */ \
\
start_timer(19); \
rdx = (uintptr_t)g->sincos2;/* Load fixed sin/cos data ptr */ \
rdx += 512*XMM_SCD3;	/* Real data sin/cos ptr */ \
rdi = (rdx+128*XMM_SCD2);	/* Complex data sin/cos ptr */ \
loops_init(256*2*clm);		/* 256 iterations of 2*clm */ \
\
r4_x4cl_eight_reals_unfft_preload; \
for(unsigned int loopA = 8; loopA; loopA--) { \
	for(unsigned int loopB = 8; loopB; loopB--) { \
		for(unsigned int loopC = 2*clm; loopC; loopC--) { \
			r4_x4cl_eight_reals_unfft(rsi, 64, 8*clmblkdst8, 16*clmblkdst8, rdi, 64*XMM_SCD2, rdx, 64*XMM_SCD2); \
		} \
		rsi += -2*clm*64+clmblkdst; /* Next source pointer */ \
		rdi += XMM_SCD2;		/* Next sine/cosine pointer */ \
		rdx += XMM_SCD2;		/* Next sine/cosine pointer */ \
	} \
	rsi += -8*clmblkdst+clmblkdst8; /* Next source pointer */ \
} \
rsi += -8*clmblkdst8+32*clmblkdst8; /* Next source pointer */ \
rdi += -64*XMM_SCD2;	/* Restore sine/cosine pointer */ \
\
loops_reset; \
r4_x4cl_four_complex_djbunfft_preload; \
for(unsigned int loopA = 3; loopA; loopA--) { \
	for(unsigned int loopB = 8; loopB; loopB--) { \
		for(unsigned int loopC = 8; loopC; loopC--) { \
			for(unsigned int loopD = 2*clm; loopD; loopD--) { \
				r4_x4cl_four_complex_djbunfft(rsi, 64, 8*clmblkdst8, 16*clmblkdst8, rdi, 64*XMM_SCD2); \
			} \
			rsi += -2*clm*64+clmblkdst; /* Next source pointer */ \
			rdi += XMM_SCD2;		/* Next sin/cos ptr */ \
		} \
		rsi += -8*clmblkdst+clmblkdst8; /* Next source pointer */ \
	} \
	rsi += -8*clmblkdst8+32*clmblkdst8; /* Next source pointer */ \
	rdi += -64*XMM_SCD2;	/* Restore sine/cosine pointer */ \
} \
rsi += -4*32*clmblkdst8;	/* Restore source pointer */ \
end_timer(19); \
\
/* Do inverse FFT levels 1,2,3 */ \
/* On input the 64-byte cache lines hold these data values: */ \
/*	0K	256K	512K	768K */ \
/*	1K	... */ \
/*	... */ \
/*	255K */ \
/*	1024K	... */ \
/*	... */ \
/* On output the 64-byte cache lines hold these data values: */ \
/*	0K	1024K	2048K	3072K */ \
/*	1K	... */ \
/*	... */ \
/*	1023K */ \
\
/* Do 512 eight_reals_last_unfft operations */ \
/*	distance between fft data elements is 512K */ \
/*	do 256 macros each processing 32 data values */ \
\
start_timer(20); \
rdi = (uintptr_t)g->sincos2;/* Load fixed sin/cos postmultipliers */ \
loops_init(256*2*clm);		/* 256 iterations of 2*clm */ \
r4_x4cl_eight_reals_last_unfft_preload; \
for(unsigned int loopA = 32; loopA; loopA--) { \
	for(unsigned int loopB = 8; loopB; loopB--) { \
		for(unsigned int loopC = 2*clm; loopC; loopC--) { \
			r4_x4cl_eight_reals_last_unfft(rsi, 64, 32*clmblkdst8, 64*clmblkdst8, rdi, 256*XMM_SCD3); \
		} \
		newxloop1nopre(clm, c1b, clmblkdst, rdi, XMM_SCD3); /* 2*clm iterations */ \
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
copy_scratch_data_to_fft(1024, clm); \
goto	endlp;			/* Skip FFT code if flag not set */ \
\
/* Fire up auxiliary threads to do pass 1 forward FFT */ \
\
/* pass1: */pass1_forward_fft_setup(b0b); \
\
/* Copy data from the FFT data area back to the scratch area */ \
\
b0b:	copy_fft_data_to_scratch(1024, clm); \
\
/* Do FFT levels 1,2,3 */ \
/* */ \
/* On input the 64-byte cache lines hold these data values: */ \
/*	0K	1024K	2048K	3072K */ \
/*	1K	... */ \
/*	... */ \
/*	1023K */ \
/* On output the 64-byte cache lines hold these data values: */ \
/*	0K	512K	1024K	1536K */ \
/*	1K	... */ \
/*	... */ \
/*	511K	... */ \
/*	2048K	... */ \
/*	... */ \
\
/* Do 512 eight_reals_first_fft operations */ \
/*	distance between fft data elements is 512K */ \
/*	do 512 iterations each processing 16 data values */ \
\
yes_fft:start_timer(22); \
rsi = (uintptr_t)g->scratch_area;	/* Get address of scratch area */ \
rdi = (uintptr_t)g->sincos2;/* Load fixed sin/cos postmultipliers */ \
set_data_prefetch_ptrs;		/* Init prefetching pointers */ \
loops_init(512*2*clm);		/* 64 iters of 4 iters of 2*2*clm */ \
r4_x2cl_eight_reals_first_fft_scratch_preload; \
for(unsigned int loopA = 64; loopA; loopA--) { \
	for(unsigned int loopB = 4*2*2*clm; loopB; ) { \
		xtouch(xptr(rcx+rbp+(clm-1)*128));	/* Preload TLB */ \
		do{ \
			do{ \
				r4_x2cl_eight_reals_first_fft_scratch(rsi, 64, 64*clmblkdst8, rdi); \
			}while(--loopB & (2*clm-1)); \
			newxloop2(clm, b1b, clmblkdst, rdi, XMM_SCD3); /* 2*2*clm iters prefetching 2*clm cache lines */ \
		}while(loopB & (2*2*clm-1)); \
	} \
	rsi += -4*2*clmblkdst+clmblkdst8; \
} \
rsi += -64*clmblkdst8;	/* Restore source pointer */ \
rdx = rdi;		/* Next real data sin/cos ptr */ \
end_timer(22); \
\
/* Do FFT levels 4,5 */ \
/* */ \
/* On input the 64-byte cache lines hold these data values: */ \
/*	0K	512K	1024K	1536K */ \
/*	1K	... */ \
/*	... */ \
/*	511K	... */ \
/*	2048K	... */ \
/*	... */ \
/* On output the 64-byte cache lines hold these data values: */ \
/*	0K	128K	256K	384K */ \
/*	1K	... */ \
/*	... */ \
/*	127K	... */ \
/*	512K	... */ \
/*	... */ \
\
/* Do 128 eight_reals_fft operations */ \
/* Do 128 four_complex_fft operations */ \
/*	distance between fft data elements is 128K */ \
/*	do 128 macros each processing 32 data values */ \
/* Do 256 four_complex_fft operations */ \
/*	distance between fft data elements is 128K */ \
/*	do 128 macros each processing 32 data values */ \
\
start_timer(23); \
rdi = (rdx+128*XMM_SCD2);	/* Complex data sin/cos ptr */ \
loops_init(128*2*clm);		/* 16 iters of 8 of 2*clm */ \
r4_x4cl_eight_reals_four_complex_djbfft_preload; \
for(unsigned int loopA = 16; loopA; loopA--) { \
	for(unsigned int loopB = 8*2*clm; loopB; ) { \
		xtouch(xptr(rcx+rbp+(clm-1)*128));	/* Preload TLB */ \
		do{ \
			r4_x4cl_eight_reals_four_complex_djbfft(rsi, 64, 16*clmblkdst8, 32*clmblkdst8, rdi, rdx); \
		}while(--loopB & (2*clm-1)); \
		newxloop1(clm, b2b, clmblkdst, rdi, XMM_SCD2, rdx, XMM_SCD2); /* 2*clm iterations */ \
	} \
	rsi += -8*clmblkdst+clmblkdst8; \
} \
rsi += -16*clmblkdst8+64*clmblkdst8; /* Next source pointer */ \
rdi += -128*XMM_SCD2;	/* Restore sin/cos ptr */ \
\
loops_init(128*2*clm);		/* 16 iters of 8 of 2*clm */ \
r4_x4cl_four_complex_djbfft_preload; \
for(unsigned int loopA = 16; loopA; loopA--) { \
	for(unsigned int loopB = 8*2*clm; loopB; ) { \
		xtouch(xptr(rcx+rbp+(clm-1)*128));	/* Preload TLB */ \
		do{ \
			r4_x4cl_four_complex_djbfft(rsi, 64, 16*clmblkdst8, 32*clmblkdst8, rdi); \
		}while(--loopB & (2*clm-1)); \
		newxloop1(clm, b3b, clmblkdst, rdi, XMM_SCD2);	/* 2*clm iterations */ \
	} \
	rsi += -8*clmblkdst+clmblkdst8; \
} \
rsi += -16*clmblkdst8-64*clmblkdst8; /* Restore source pointer */ \
end_timer(23); \
\
/* Do FFT levels 6,7 */ \
/* */ \
/* On input the 64-byte cache lines hold these data values: */ \
/*	0K	128K	256K	384K */ \
/*	1K	... */ \
/*	... */ \
/*	127K	... */ \
/*	512K	... */ \
/*	... */ \
/* On output the 64-byte cache lines hold these data values: */ \
/*	0K	32K	64K	96K */ \
/*	1K	... */ \
/*	... */ \
/*	31K	... */ \
/*	128K	... */ \
/*	... */ \
\
/* Do 32 eight_reals_fft operations */ \
/* Do 32 four_complex_fft operations */ \
/*	distance between fft data elements is 32K */ \
/*	do 32 macros each processing 32 data values */ \
/* Do 448 four_complex_fft operations */ \
/*	distance between fft data elements is 32K */ \
/*	do 224 macros each processing 32 data values */ \
\
start_timer(24); \
rdx = (uintptr_t)g->premult_addr;/* Address of premultiplier/sin/cos table */ \
rdx += 8*2*clm*XMM_SCND2 + 8*2*clm*XMM_SCND2R + 2*clm*XMM_SCD4 + 16*2*clm*XMM_SCD8; /* Real data sin/cos ptr */ \
rdi = (rdx+32*2*clm*XMM_SCD2); /* Real data sin/cos ptr */ \
loops_init(256*2*clm);		/* 8 iters of 32 of 2*clm */ \
r4_x4cl_eight_reals_four_complex_djbfft_preload; \
for(unsigned int loopA = 4; loopA; loopA--) { \
	for(unsigned int loopB = 8*2*clm; loopB; ) { \
		xtouch(xptr(rcx+rbp+(clm-1)*128));	/* Preload TLB */ \
		do{ \
			r4_x4cl_eight_reals_four_complex_djbfft(rsi, 64, 4*clmblkdst8, 8*clmblkdst8, rdi, rdx); \
			rdi += XMM_SCD2;		/* Next sine/cosine pointer */ \
			rdx += XMM_SCD2;		/* Next sine/cosine pointer */ \
		}while(--loopB & (2*clm-1)); \
		newxloop1(clm, b4b, clmblkdst);	/* 2*clm iterations */ \
	} \
	rsi += -8*clmblkdst+clmblkdst8; \
} \
rsi += -4*clmblkdst8+16*clmblkdst8; /* Next source pointer */ \
rdi += -32*2*clm*XMM_SCD2;	/* Restore sin/cos ptr */ \
\
loops_reset; \
r4_x4cl_four_complex_djbfft_preload; \
for(unsigned int loopA = 7; loopA; loopA--) { \
	for(unsigned int loopB = 4; loopB; loopB--) { \
		for(unsigned int loopC = 8*2*clm; loopC; ) { \
			xtouch(xptr(rcx+rbp+(clm-1)*128));	/* Preload TLB */ \
			do{ \
				r4_x4cl_four_complex_djbfft(rsi, 64, 4*clmblkdst8, 8*clmblkdst8, rdi); \
				rdi += XMM_SCD2;		/* Next sine/cosine pointer */ \
			}while(--loopC & (2*clm-1)); \
			newxloop1(clm, b5b, clmblkdst);	/* 2*clm iterations */ \
		} \
		rsi += -8*clmblkdst+clmblkdst8; \
	} \
	rsi += -4*clmblkdst8+16*clmblkdst8; /* Next source pointer */ \
	rdi += -32*2*clm*XMM_SCD2;	/* Restore sin/cos ptr */ \
} \
rsi += -8*16*clmblkdst8;	/* Restore source pointer */ \
end_timer(24); \
\
/* Do FFT levels 8,9 */ \
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
/* Do 8 eight_reals_fft operations */ \
/* Do 8 four_complex_fft operations */ \
/*	distance between fft data elements is 8K */ \
/*	do 8 macros each processing 32 data values */ \
/* Do 496 four_complex_fft operations */ \
/*	distance between fft data elements is 8K */ \
/*	do 248 macros each processing 32 data values */ \
\
start_timer(25); \
rdi += -32*2*clm*XMM_SCD2-8*2*clm*XMM_SCND2; /* Complex data sin/cos ptr */ \
rdx = (rdi-8*2*clm*XMM_SCND2R); /* Real data sin/cos ptr */ \
loops_init(256*2*clm);		/* 32 iterations of 8 iters of 2*clm */ \
r4_x4cl_wpn_eight_reals_four_complex_djbfft_preload; \
for(unsigned int loopA = 8*2*clm; loopA; ) { \
	xtouch(xptr(rcx+rbp+(clm-1)*128));	/* Preload TLB */ \
	do{ \
		r4_x4cl_wpn_eight_reals_four_complex_djbfft(rsi, 64, clmblkdst8, 2*clmblkdst8, rdi, rdx); \
		rdi += XMM_SCND2;		/* Next sine/cosine pointer */ \
		rdx += XMM_SCND2R;		/* Next sine/cosine pointer */ \
	}while(--loopA & (2*clm-1)); \
	newxloop1(clm, b6b, clmblkdst);	/* 2*clm iterations */ \
} \
rsi += -8*clmblkdst+4*clmblkdst8; \
rdi += -8*2*clm*XMM_SCND2;	/* Restore sin/cos ptr */ \
\
loops_reset; \
r4_x4cl_wpn_four_complex_djbfft_preload; \
for(unsigned int loopA = 31; loopA; loopA--) { \
	for(unsigned int loopB = 8*2*clm; loopB; ) { \
		xtouch(xptr(rcx+rbp+(clm-1)*128));	/* Preload TLB */ \
		do{ \
			r4_x4cl_wpn_four_complex_djbfft(rsi, 64, clmblkdst8, 2*clmblkdst8, rdi); \
			rdi += XMM_SCND2;		/* Next sine/cosine pointer */ \
		}while(--loopB & (2*clm-1)); \
		newxloop1(clm, b7b, clmblkdst);	/* 2*clm iterations */ \
	} \
	rsi += -8*clmblkdst+4*clmblkdst8; /* Next source pointer */ \
	rdi += -8*2*clm*XMM_SCND2;	/* Restore sin/cos ptr */ \
} \
rsi += -32*4*clmblkdst8;	/* Restore source pointer */ \
end_timer(25); \
\
/* Do FFT levels 10,11,12 */ \
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
/* Do 254 eight_complex_fft8 operations */ \
/*	distance between fft data elements is 1K */ \
/*	do 127 macros each processing 64 data values */ \
\
start_timer(26); \
rdi += -8*2*clm*XMM_SCND2R-16*2*clm*XMM_SCD8+128; /* Complex data sin/cos ptr and use -128 */ \
;					/* addressing to generate 1 byte offsets in radix-8 macro */ \
rdx = (rdi-128-2*clm*XMM_SCD4); /* Real data sin/cos ptr */ \
rcx = (uintptr_t)g->data_addr;/* Destination pointer */ \
loops_init_prefetch(128*2*clm, 64, 1, rbp); /* 128 iters of 2*clm */ \
pfing(rbp = (uintptr_t)g->premult_prefetch);/* Load premult/sin/cos prefetch pointer */ \
r8_sg8cl_sixteen_reals_eight_complex_fft8_preload; \
for(unsigned int loopA = 2*clm; loopA; loopA--) { \
	r8_sg8cl_sixteen_reals_eight_complex_fft8(rsi, 64, clmblkdst, 2*clmblkdst, 4*clmblkdst, rcx, 64, blkdst, 2*blkdst, 4*blkdst, rdi-128, rdx); \
	rdi += XMM_SCD8;		/* Next sine/cosine pointer */ \
	rdx += XMM_SCD4;		/* Next sine/cosine pointer */ \
} \
rsi += -2*clm*64+clmblkdst8; /* Next source pointer */ \
rcx += -2*clm*64+8*blkdst;	/* Next dest pointer */ \
rdi += -2*clm*XMM_SCD8;	/* Restore sin/cos pointer */ \
\
loops_reset; \
for(unsigned int loopA = 16*8-1; loopA; ) { \
	do{ \
		for(unsigned int loopC = 2*clm; loopC; loopC--) { \
			r8_sg8cl_eight_complex_fft8(rsi, 64, clmblkdst, 2*clmblkdst, 4*clmblkdst, rcx, 64, blkdst, 2*blkdst, 4*blkdst, rdi-128); \
			rdi += XMM_SCD8;		/* Next premult/sin/cos ptr */ \
		} \
		rsi += -2*clm*64+clmblkdst8; /* Next source pointer */ \
		rcx += -2*clm*64+8*blkdst;	/* Next dest pointer */ \
		rdi += -2*clm*XMM_SCD8;	/* Restore premult/sine/cosine pointer */ \
	}while(--loopA & 7); \
	rdi += 2*clm*XMM_SCD8;	/* Next premult/sine/cosine pointer */ \
} \
end_timer(26); \
\
/* Work on next 4*clm sets of 4096 values */ \
\
endlp:	xpass1_get_next_block(pass2, c0b, b0b, _wpn); /* Get next block and jump */ \
/* done: */ xfft_footer_rn; \
\
\
\
/* Perform a pass 1 negacyclic FFT of length 128 doubles. */ \
/* Memory examples are for a 128K FFT with 10 levels done in pass 2. */ \
/* Sin/cos data needed: 4*2*clm*XMM_SCND2 + 4*2*clm*XMM_SCD4 = 16*2*clm 64-byte cache lines */ \
\

#define r4dwpn_pass1sc128ac(procname, xpass2name, pass2_levels, clm) /* (fftname, pass2_levels, clmarg) */ \
\
/* set_FFT_constants(fftname, 0, 1, pass2_levels, clmarg); */ \
xfft_header_rn(procname, xpass2name, pass2_levels, clm); \
\
/* Do pass 2 (the last 10 FFT levels if this is a 128K FFT) */ \
\
/* pass2: */start_timer(0); \
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
/* Do 16 four_complex_unfft4 operations */ \
/*	distance between fft data elements is 1K */ \
/*	do 16 macros each processing 16 data values */ \
\
c0b:	start_timer(16); \
rsi = (uintptr_t)g->data_addr;/* Get source address */ \
rcx = (uintptr_t)g->scratch_area;	/* Get address of scratch area */ \
rdi = (uintptr_t)g->premult_addr;/* Load premult/sin/cos pointer */ \
pfing(rbp = (rdi+(16-8)*2*clm*64)); /* Finish prefetching the sin/cos/premult data */ \
loops_init_prefetch(16*2*clm, 64, 2, rbp); /* 16*2*clm iterations */ \
/*bug	xtouchpages rbp, 16*2*clm * 64 / 2 ;; Prefetch TLBs */ \
for(unsigned int loopA = 4; loopA; loopA--) { \
	for(unsigned int loopB = 4; loopB; loopB--) { \
		for(unsigned int loopC = 2*clm; loopC; loopC--) { \
			r4_sg2cl_four_complex_unfft4(rsi, 64, blkdst, rcx, 64, clmblkdst, rdi); \
			rdi += XMM_SCD4;		/* Next premult/sin/cos data */ \
		} \
		rsi += -2*clm*64+2*blkdst;	/* Next source pointer */ \
		rcx += -2*clm*64+2*clmblkdst; /* Next dest pointer */ \
		rdi += -2*clm*XMM_SCD4;	/* Restore premult/sin/cos data pointer */ \
	} \
	rcx += -4*2*clmblkdst+clmblkdst8; \
	rdi += 2*clm*XMM_SCD4;	/* Next premult/sin/cos pointer */ \
} \
rsi = (uintptr_t)g->scratch_area;	/* Restore source pointer */ \
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
r4_x4cl_wpn_four_complex_djbunfft_preload; \
for(unsigned int loopA = 4; loopA; loopA--) { \
	for(unsigned int loopB = 2*2*clm; loopB; ) { \
		xtouch(xptr(rcx+rbp+(clm-1)*128));	/* Preload TLB */ \
		do{ \
			r4_x4cl_wpn_four_complex_djbunfft(rsi, 64, 2*clmblkdst, 4*clmblkdst, rdi, 2*2*clm*XMM_SCND2); \
			rdi += XMM_SCND2;		/* Next clm's sin/cos data */ \
		}while(--loopB & (2*clm-1)); \
		newxloop1(clm, c7b, clmblkdst);	/* 2*clm iterations */ \
	} \
	rsi += -2*clmblkdst+clmblkdst8; /* Next source pointer */ \
	rdi += -2*2*clm*XMM_SCND2;	/* Restore sine/cosine pointer */ \
} \
rsi += -4*clmblkdst8;	/* Restore source pointer */ \
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
rdx = (uintptr_t)g->sincos1;/* Load fixed negacyclic premultipliers */ \
rdi = (uintptr_t)g->sincos2;/* Load fixed sin/cos postmultipliers */ \
loops_init(8*2*clm);		/* 8 iterations of 2*clm */ \
for(unsigned int loopA = 8*2*clm; loopA; ) { \
	xtouch(xptr(rcx+rbp+(clm-1)*128));	/* Preload TLB */ \
	do{ \
		r4_x4cl_four_complex_last_djbunfft(rsi, 64, clmblkdst8, 2*clmblkdst8, rdi, 8*XMM_SCD2, rdx, 8*XMM_SCD4); \
	}while(--loopA & (2*clm-1)); \
	newxloop1(clm, c1b, clmblkdst, rdi, XMM_SCD2, rdx, XMM_SCD4); /* 2*clm iterations */ \
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
/* pass1: */pass1_forward_fft_setup(b0b); \
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
rdx = (uintptr_t)g->sincos1;/* Address of fixed negacyclic premultipliers */ \
rdi = (uintptr_t)g->sincos2;/* Load fixed sin/cos postmultipliers */ \
set_data_prefetch_ptrs;		/* Init prefetching pointers */ \
pfing(rcx += 16*blkdst);		/* We've already prefetched 16 blocks */ \
loops_init(16*2*clm);		/* 2 iters of 4 iters of 2*2*clm */ \
for(unsigned int loopA = 2; loopA; loopA--) { \
	for(unsigned int loopB = 4*2*2*clm; loopB; ) { \
		xtouch(xptr(rcx+rbp+(clm-1)*128));	/* Preload TLB */ \
		do{ \
			do{ \
				r4_x2cl_four_complex_first_djbfft_scratch(rsi, 64, 2*clmblkdst8, rdi, rdx); \
			}while(--loopB & (2*clm-1)); \
			newxloop2(clm, b1b, clmblkdst, rdi, XMM_SCD2, rdx, XMM_SCD4); /* 2*2*clm iters prefetching 2*clm cache lines */ \
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
/*	do 8 macros each processing 8 data values */ \
\
start_timer(25); \
rdi = (uintptr_t)g->premult_addr;/* Address of premultiplier/sin/cos table */ \
rdi += 4*2*clm*XMM_SCD4; \
loops_init(8*2*clm);		/* 2 iterations of 4 iters of 2*clm */ \
r4_x4cl_wpn_four_complex_djbfft_preload; \
for(unsigned int loopA = 2; loopA; loopA--) { \
	for(unsigned int loopB = 4*2*clm; loopB; ) { \
		xtouch(xptr(rcx+rbp+(clm-1)*128));	/* Preload TLB */ \
		do{ \
			r4_x4cl_wpn_four_complex_djbfft(rsi, 64, 4*clmblkdst, clmblkdst8, rdi); \
			rdi += XMM_SCND2;		/* Next sine/cosine pointer */ \
		}while(--loopB & (2*clm-1)); \
		newxloop1(clm, b7b, clmblkdst);	/* 2*clm iterations */ \
	} \
	rsi += -4*clmblkdst+2*clmblkdst8; /* Next source pointer */ \
	rdi += -4*2*clm*XMM_SCND2;	/* Restore sin/cos ptr */ \
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
rdi += -4*2*clm*XMM_SCD4;	/* Next premult/sin/cos ptr */ \
rcx = (uintptr_t)g->data_addr;/* Destination pointer */ \
loops_init_prefetch(8*2*clm, 64, 1, rbp); /* 8 iters of 2*clm */ \
pfing(rbp = (uintptr_t)g->premult_prefetch);/* Load premult/sin/cos prefetch pointer */ \
for(unsigned int loopA = 4; loopA; loopA--) { \
	for(unsigned int loopB = 2; loopB; loopB--) { \
		for(unsigned int loopC = 2*clm; loopC; loopC--) { \
			r4_sg4cl_four_complex_fft4(rsi, 64, clmblkdst, 2*clmblkdst, rcx, 64, blkdst, 2*blkdst, rdi); \
			rdi += XMM_SCD4;		/* Next premult/sin/cos ptr */ \
		} \
		rsi += -2*clm*64+4*clmblkdst; /* Next source pointer */ \
		rcx += -2*clm*64+4*blkdst;	/* Next dest pointer */ \
		rdi += -2*clm*XMM_SCD4;	/* Restore premult/sine/cosine pointer */ \
	} \
	rsi += -2*4*clmblkdst+clmblkdst8; /* Next source pointer */ \
	rdi += 2*clm*XMM_SCD4;	/* Next premult/sine/cosine pointer */ \
} \
end_timer(26); \
\
/* Work on next 4*clm sets of 128 values */ \
\
endlp:	xpass1_get_next_block(pass2, c0b, b0b, _wpn); /* Get next block and jump */ \
/* done: */ xfft_footer_rn; \
\
\
\
/* Perform a pass 1 negacyclic FFT of length 256 doubles. */ \
/* Memory examples are for a 256K FFT with 10 levels done in pass 2. */ \
/* Sin/cos data needed: 8*2*clm*XMM_SCND2 + 4*2*clm*XMM_SCD8 = 32*2*clm 64-byte cache lines */ \
\

#define r4dwpn_pass1sc256ac(procname, xpass2name, pass2_levels, clm) /* (fftname, pass2_levels, clmarg) */ \
\
/* set_FFT_constants(fftname, 0, 1, pass2_levels, clmarg); */ \
xfft_header_rn(procname, xpass2name, pass2_levels, clm); \
\
/* Do pass 2 (the last 10 FFT levels if this is a 256K FFT) */ \
\
/* pass2: */start_timer(0); \
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
pfing(rbp = (rdi+(32-16)*2*clm*64)); /* Finish prefetching the sin/cos/premult data */ \
loops_init_prefetch(16*2*clm, 64, 1, rbp); /* 16*2*clm iterations */ \
/*bug	xtouchpages rbp, 16*2*clm * 128 / 1 ;; Prefetch TLBs */ \
rdi += 128;		/* Generates 1 byte offsets in radix-8 macro */ \
for(unsigned int loopA = 4; loopA; loopA--) { \
	for(unsigned int loopB = 2; loopB; loopB--) { \
		for(unsigned int loopC = 2; loopC; loopC--) { \
			for(unsigned int loopD = 2*clm; loopD; loopD--) { \
				r8_sg4cl_eight_complex_unfft8(rsi, 64, blkdst, 2*blkdst, rcx, 64, clmblkdst, 2*clmblkdst, rdi-128); \
				rdi += XMM_SCD8;		/* Next premult/sin/cos data */ \
			} \
			rsi += -2*clm*64+4*blkdst;	/* Next source pointer */ \
			rcx += -2*clm*64+4*clmblkdst; /* Next dest pointer */ \
			rdi += -2*clm*XMM_SCD8;	/* Restore premult/sin/cos data pointer */ \
		} \
		rcx += -2*4*clmblkdst+clmblkdst8; \
	} \
	rdi += 2*clm*XMM_SCD8;	/* Next premult/sin/cos pointer */ \
} \
rsi = (uintptr_t)g->scratch_area;	/* Restore source pointer */ \
rdi += -128;		/* Next sin/cos data pointer */ \
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
r4_x4cl_wpn_four_complex_djbunfft_preload; \
for(unsigned int loopA = 4; loopA; loopA--) { \
	for(unsigned int loopB = 4*2*clm; loopB; ) { \
		xtouch(xptr(rcx+rbp+(clm-1)*128));	/* Preload TLB */ \
		do{ \
			r4_x4cl_wpn_four_complex_djbunfft(rsi, 64, 4*clmblkdst, clmblkdst8, rdi, 4*2*clm*XMM_SCND2); \
			rdi += XMM_SCND2;		/* Next clm's sin/cos data */ \
		}while(--loopB & (2*clm-1)); \
		newxloop1(clm, c7b, clmblkdst);	/* Loop 2*clm times */ \
	} \
	rsi += -4*clmblkdst+2*clmblkdst8; /* Next source pointer */ \
	rdi += -4*2*clm*XMM_SCND2;	/* Restore sine/cosine pointer */ \
} \
rsi += -4*2*clmblkdst8;	/* Restore source pointer */ \
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
rdx = (uintptr_t)g->sincos1;/* Load fixed negacyclic premultipliers */ \
rdi = (uintptr_t)g->sincos2;/* Load fixed sin/cos postmultipliers */ \
loops_init(16*2*clm);		/* 16 iterations of 2*clm */ \
for(unsigned int loopA = 2; loopA; loopA--) { \
	for(unsigned int loopB = 8*2*clm; loopB; ) { \
		xtouch(xptr(rcx+rbp+(clm-1)*128));	/* Preload TLB */ \
		do{ \
			r4_x4cl_four_complex_last_djbunfft(rsi, 64, 2*clmblkdst8, 4*clmblkdst8, rdi, 16*XMM_SCD2, rdx, 16*XMM_SCD4); \
		}while(--loopB & (2*clm-1)); \
		newxloop1(clm, c1b, clmblkdst, rdi, XMM_SCD2, rdx, XMM_SCD4); /* 2*clm iterations */ \
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
/* pass1: */pass1_forward_fft_setup(b0b); \
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
rdx = (uintptr_t)g->sincos1;/* Address of fixed negacyclic premultipliers */ \
rdi = (uintptr_t)g->sincos2;/* Load fixed sin/cos postmultipliers */ \
set_data_prefetch_ptrs;		/* Init prefetching pointers */ \
pfing(rcx += 32*blkdst);		/* We've already prefetched 32 blocks */ \
loops_init(32*2*clm);		/* 4 iters of 4 iters of 2*2*clm */ \
for(unsigned int loopA = 4; loopA; loopA--) { \
	for(unsigned int loopB = 4*2*2*clm; loopB; ) { \
		xtouch(xptr(rcx+rbp+(clm-1)*128));	/* Preload TLB */ \
		do{ \
			do{ \
				r4_x2cl_four_complex_first_djbfft_scratch(rsi, 64, 4*clmblkdst8, rdi, rdx); \
			}while(--loopB & (2*clm-1)); \
			newxloop2(clm, b1b, clmblkdst, rdi, XMM_SCD2, rdx, XMM_SCD4); /* 2*2*clm iters prefetching 2*clm cache lines */ \
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
/* Do 32 four_complex_fft operations */ \
/*	distance between fft data elements is 8K */ \
/*	do 16 macros each processing 32 data values */ \
\
start_timer(25); \
rdi = (uintptr_t)g->premult_addr;/* Address of premultiplier/sin/cos table */ \
rdi += 4*2*clm*XMM_SCD8; \
loops_init(16*2*clm);		/* 2 iterations of 8 iters of 2*clm */ \
r4_x4cl_wpn_four_complex_djbfft_preload; \
for(unsigned int loopA = 2; loopA; loopA--) { \
	for(unsigned int loopB = 8*2*clm; loopB; ) { \
		xtouch(xptr(rcx+rbp+(clm-1)*128));	/* Preload TLB */ \
		do{ \
			r4_x4cl_wpn_four_complex_djbfft(rsi, 64, clmblkdst8, 2*clmblkdst8, rdi); \
			rdi += XMM_SCND2;		/* Next sine/cosine pointer */ \
		}while(--loopB & (2*clm-1)); \
		newxloop1(clm, b7b, clmblkdst);	/* 2*clm iterations */ \
	} \
	rsi += -8*clmblkdst+4*clmblkdst8; /* Next source pointer */ \
	rdi += -8*2*clm*XMM_SCND2;	/* Restore sin/cos ptr */ \
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
rdi += -4*2*clm*XMM_SCD8+128; /* Next premult/sin/cos ptr (gen 1 byte offsets in radix-8 macro) */ \
rcx = (uintptr_t)g->data_addr;/* Destination pointer */ \
loops_init_prefetch(8*2*clm, 128, 1, rbp); /* 8 iters of 2*clm */ \
pfing(rbp = (uintptr_t)g->premult_prefetch);/* Load premult/sin/cos prefetch pointer */ \
for(unsigned int loopA = 4; loopA; loopA--) { \
	for(unsigned int loopB = 2; loopB; loopB--) { \
		for(unsigned int loopC = 2*clm; loopC; loopC--) { \
			r8_sg8cl_eight_complex_fft8(rsi, 64, clmblkdst, 2*clmblkdst, 4*clmblkdst, rcx, 64, blkdst, 2*blkdst, 4*blkdst, rdi-128); \
			rdi += XMM_SCD8;		/* Next premult/sin/cos ptr */ \
		} \
		rsi += -2*clm*64+clmblkdst8; /* Next source pointer */ \
		rcx += -2*clm*64+8*blkdst;	/* Next dest pointer */ \
		rdi += -2*clm*XMM_SCD8;	/* Restore premult/sine/cosine pointer */ \
	} \
	rdi += 2*clm*XMM_SCD8;	/* Next premult/sine/cosine pointer */ \
} \
end_timer(26); \
\
/* Work on next 4*clm sets of 256 values */ \
\
endlp:	xpass1_get_next_block(pass2, c0b, b0b, _wpn); /* Get next block and jump */ \
/* done: */ xfft_footer_rn; \
\
\
\
/* Perform a pass 1 negacyclic FFT of length 384 doubles. */ \
/* Memory examples are for a 384K FFT with 10 levels done in pass 2. */ \
/* Sin/cos data needed: 16*2*clm*XMM_SCD1 + 4*2*clm*XMM_SCND2 + 4*2*clm*XMM_SCD4 = 24*2*clm 64-byte cache lines */ \
\

#define r4dwpn_pass1sc384ac(procname, xpass2name, pass2_levels, clm) /* (fftname, pass2_levels, clmarg) */ \
\
/* set_FFT_constants(fftname, 0, 1, pass2_levels, clmarg); */ \
xfft_header_rn(procname, xpass2name, pass2_levels, clm); \
\
/* Do pass 2 (the last 10 FFT levels if this is a 384K FFT) */ \
\
/* pass2: */start_timer(0); \
xpass2name(g); \
end_timer(0); \
if(g->ffttype == 1)		/* We're done if FFTing only */ \
	return; /* goto done; */ \
\
/* Do the inverse FFT */ \
\
pass1_inverse_fft_setup(c0b, procname); \
\
/* Do inverse FFT levels 8,9 */ \
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
/* Do 48 four_complex_unfft4 operations */ \
/*	distance between fft data elements is 1K */ \
/*	do 48 macros each processing 16 data values */ \
\
c0b:	start_timer(16); \
rsi = (uintptr_t)g->data_addr;/* Get source address */ \
rcx = (uintptr_t)g->scratch_area;	/* Get address of scratch area */ \
rdi = (uintptr_t)g->premult_addr;/* Load premult/sin/cos pointer */ \
pfing(rbp = (rdi+(24-12)*2*clm*64)); /* Finish prefetching the sin/cos/premult data */ \
loops_init_prefetch(48*2*clm, 64, 4, rbp); /* 48*2*clm iterations */ \
/*bug	xtouchpages rbp, 48*2*clm * 64 / 4 ;; Prefetch TLBs */ \
for(unsigned int loopA = 4; loopA; loopA--) { \
	for(unsigned int loopB = 3; loopB; loopB--) { \
		for(unsigned int loopC = 4; loopC; loopC--) { \
			for(unsigned int loopD = 2*clm; loopD; loopD--) { \
				r4_sg2cl_four_complex_unfft4(rsi, 64, blkdst, rcx, 64, clmblkdst, rdi); \
				rdi += XMM_SCD4;		/* Next premult/sin/cos data */ \
			} \
			rsi += -2*clm*64+2*blkdst;	/* Next source pointer */ \
			rcx += -2*clm*64+2*clmblkdst; /* Next dest pointer */ \
			rdi += -2*clm*XMM_SCD4;	/* Restore premult/sin/cos data pointer */ \
		} \
		rcx += -4*2*clmblkdst+clmblkdst8; \
	} \
	rdi += 2*clm*XMM_SCD4;	/* Next premult/sin/cos pointer */ \
} \
rsi = (uintptr_t)g->scratch_area;	/* Restore source pointer */ \
end_timer(16); \
\
/* Do inverse FFT levels 6,7 */ \
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
/* Do 48 four_complex_unfft operations */ \
/*	distance between fft data elements is 4K */ \
/*	do 24 macros each processing 32 data values */ \
\
start_timer(17); \
loops_init(24*2*clm);		/* 24*2*clm iterations */ \
r4_x4cl_wpn_four_complex_djbunfft_preload; \
for(unsigned int loopA = 12; loopA; loopA--) { \
	for(unsigned int loopB = 2; loopB; loopB--) { \
		for(unsigned int loopC = 2*clm; loopC; loopC--) { \
			r4_x4cl_wpn_four_complex_djbunfft(rsi, 64, 2*clmblkdst, 4*clmblkdst, rdi, 2*2*clm*XMM_SCND2); \
			rdi += XMM_SCND2;		/* Next clm's sin/cos data */ \
		} \
		rsi += -2*clm*64+clmblkdst; /* Next source pointer */ \
	} \
	rsi += -2*clmblkdst+clmblkdst8; /* Next source pointer */ \
	rdi += -2*2*clm*XMM_SCND2;	/* Restore sine/cosine pointer */ \
} \
rsi += -12*clmblkdst8;	/* Restore source pointer */ \
rdi += 4*2*clm*XMM_SCND2;	/* Next sin/cos data pointer */ \
end_timer(17); \
\
/* Do inverse FFT levels 4,5 (actually 1.585 levels) */ \
/* On input the 64-byte cache lines hold these data values: */ \
/*	0K	8K	16K	24K */ \
/*	1K	... */ \
/*	... */ \
/*	7K */ \
/*	32K	... */ \
/*	... */ \
/* On output the 64-byte cache lines hold these data values: */ \
/*	0K	24K	48K	72K */ \
/*	1K	... */ \
/*	... */ \
/*	23K */ \
/*	96K	... */ \
/*	... */ \
\
/* Do 64 three_complex_unfft operations */ \
/*	distance between fft data elements is 16K */ \
/*	do 32 macros each processing 24 data values */ \
\
start_timer(18); \
loops_init(32*2*clm);		/* 32 iterations of 2*clm */ \
r3_x3cl_three_complex_djbunfft_preload; \
for(unsigned int loopA = 4; loopA; loopA--) { \
	for(unsigned int loopB = 8; loopB; loopB--) { \
		for(unsigned int loopC = 2*clm; loopC; loopC--) { \
			r3_x3cl_three_complex_djbunfft(rsi, 64, clmblkdst8, rdi, 8*2*clm*XMM_SCD1); \
			rdi += XMM_SCD1;		/* Next clm's sin/cos data */ \
		} \
		rsi += -2*clm*64+clmblkdst; /* Next source pointer */ \
	} \
	rsi += -8*clmblkdst+3*clmblkdst8; /* Next source pointer */ \
	rdi += -8*2*clm*XMM_SCD1;	/* Restore sine/cosine pointer */ \
} \
rsi += -4*3*clmblkdst8;	/* Restore source pointer */ \
end_timer(18); \
\
/* Do inverse FFT levels 1,2,3 */ \
/* On input the 64-byte cache lines hold these data values: */ \
/*	0K	24K	48K	72K */ \
/*	1K	... */ \
/*	... */ \
/*	23K */ \
/*	96K	... */ \
/*	... */ \
/* On output the 64-byte cache lines hold these data values: */ \
/*	0K	96K	192K	288K */ \
/*	1K	... */ \
/*	... */ \
/*	95K */ \
\
/* Do 48 four_complex_unfft operations */ \
/*	distance between fft data elements is 48K */ \
/*	do 24 macros each processing 32 data values */ \
\
start_timer(20); \
rdx = (uintptr_t)g->sincos1;/* Load fixed negacyclic premultipliers */ \
rdi = (uintptr_t)g->sincos2;/* Load fixed sin/cos postmultipliers */ \
set_data_prefetch_ptrs;		/* Init prefetching pointers */ \
loops_init(24*2*clm);		/* 24 iterations of 2*clm */ \
for(unsigned int loopA = 3; loopA; loopA--) { \
	for(unsigned int loopB = 8*2*clm; loopB; ) { \
		xtouch(xptr(rcx+rbp+(clm-1)*128));	/* Preload TLB */ \
		do{ \
			r4_x4cl_four_complex_last_djbunfft(rsi, 64, 3*clmblkdst8, 6*clmblkdst8, rdi, 24*XMM_SCD2, rdx, 24*XMM_SCD4); \
		}while(--loopB & (2*clm-1)); \
		newxloop1(clm, c1b, clmblkdst, rdi, XMM_SCD2, rdx, XMM_SCD4); /* 2*clm iterations */ \
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
copy_scratch_data_to_fft(96, clm); \
goto	endlp;			/* Skip FFT code if flag not set */ \
\
/* Fire up auxiliary threads to do pass 1 forward FFT */ \
\
/* pass1: */pass1_forward_fft_setup(b0b); \
\
/* Copy data from the FFT data area back to the scratch area */ \
\
b0b:	copy_fft_data_to_scratch(96, clm); \
\
/* Do FFT levels 1,2,3 */ \
/* */ \
/* On input the 64-byte cache lines hold these data values: */ \
/*	0K	96K	192K	288K */ \
/*	1K	... */ \
/*	... */ \
/*	95K */ \
/* On output the 64-byte cache lines hold these data values: */ \
/*	0K	48K	96K	144K */ \
/*	1K	... */ \
/*	... */ \
/*	47K	... */ \
/*	192K	... */ \
/*	... */ \
\
/* Do 48 four_complex_first_fft operations */ \
/*	distance between fft data elements is 48K */ \
/*	do 48 macros each processing 16 data values */ \
\
yes_fft:start_timer(22); \
rsi = (uintptr_t)g->scratch_area;	/* Get address of scratch area */ \
rdx = (uintptr_t)g->sincos1;/* Address of fixed negacyclic premultipliers */ \
rdi = (uintptr_t)g->sincos2;/* Load fixed sin/cos postmultipliers */ \
set_data_prefetch_ptrs;		/* Init prefetching pointers */ \
pfing(rcx += (24-8)*blkdst);	/* We've already prefetched 24 blocks, but */ \
;					/* code below prefetches 8 blocks too many */ \
loops_init(48*2*clm);		/* 6 iters of 4 iters of 2*2*clm */ \
for(unsigned int loopA = 6; loopA; loopA--) { \
	for(unsigned int loopB = 4*2*2*clm; loopB; ) { \
		xtouch(xptr(rcx+rbp+(clm-1)*128));	/* Preload TLB */ \
		do{ \
			do{ \
				r4_x2cl_four_complex_first_djbfft_scratch(rsi, 64, 6*clmblkdst8, rdi, rdx); \
			}while(--loopB & (2*clm-1)); \
			newxloop2(clm, b1b, clmblkdst, rdi, XMM_SCD2, rdx, XMM_SCD4); /* 2*2*clm iters prefetching 2*clm cache lines */ \
		}while(loopB & (2*2*clm-1)); \
	} \
	rsi += -4*2*clmblkdst+clmblkdst8; \
} \
rsi += -6*clmblkdst8;	/* Restore source pointer */ \
end_timer(22); \
\
/* Do FFT levels 4,5 (actually 1.585 levels) */ \
/* */ \
/* On input the 64-byte cache lines hold these data values: */ \
/*	0K	48K	96K	144K */ \
/*	1K	... */ \
/*	... */ \
/*	47K	... */ \
/*	192K	... */ \
/*	... */ \
/* On output the 64-byte cache lines hold these data values: */ \
/*	0K	16K	32K	48K */ \
/*	1K	... */ \
/*	... */ \
/*	15K	... */ \
/*	64K	... */ \
/*	... */ \
\
/* Do 64 three_complex_fft operations */ \
/*	distance between fft data elements is 16K */ \
/*	do 32 macros each processing 24 data values */ \
\
start_timer(24); \
rdi = (uintptr_t)g->premult_addr;/* Address of premultiplier/sin/cos table */ \
rdi += 4*2*clm*XMM_SCND2 + 4*2*clm*XMM_SCD4; \
loops_init(32*2*clm);		/* 2 iters of 2 iters of 8 of 2*clm */ \
r3_x3cl_three_complex_djbfft_preload; \
for(unsigned int loopA = 2; loopA; loopA--) { \
	for(unsigned int loopB = 2; loopB; loopB--) { \
		for(unsigned int loopC = 8*2*clm; loopC; ) { \
			xtouch(xptr(rcx+rbp+(clm-1)*128));	/* Preload TLB */ \
			do{ \
				r3_x3cl_three_complex_djbfft(rsi, 64, 2*clmblkdst8, rdi); \
				rdi += XMM_SCD1;		/* Next sine/cosine pointer */ \
			}while(--loopC & (2*clm-1)); \
			newxloop1(clm, b5b, clmblkdst);	/* 2*clm iterations */ \
		} \
		rsi += -8*clmblkdst+clmblkdst8; \
	} \
	rsi += -2*clmblkdst8+6*clmblkdst8; /* Next source pointer */ \
	rdi += -16*2*clm*XMM_SCD1;	/* Restore sin/cos ptr */ \
} \
rsi += -2*6*clmblkdst8;	/* Restore source pointer */ \
end_timer(24); \
\
/* Do FFT levels 6,7 */ \
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
/* Do 48 four_complex_fft operations */ \
/*	distance between fft data elements is 4K */ \
/*	do 24 macros each processing 32 data values */ \
\
start_timer(25); \
rdi += -4*2*clm*XMM_SCND2;	/* Next sin/cos ptr */ \
loops_init(24*2*clm);		/* 6 iterations of 4 iters of 2*clm */ \
r4_x4cl_wpn_four_complex_djbfft_preload; \
for(unsigned int loopA = 6; loopA; loopA--) { \
	for(unsigned int loopB = 4*2*clm; loopB; ) { \
		xtouch(xptr(rcx+rbp+(clm-1)*128));	/* Preload TLB */ \
		do{ \
			r4_x4cl_wpn_four_complex_djbfft(rsi, 64, 4*clmblkdst, clmblkdst8, rdi); \
			rdi += XMM_SCND2;		/* Next sine/cosine pointer */ \
		}while(--loopB & (2*clm-1)); \
		newxloop1(clm, b7b, clmblkdst);	/* 2*clm iterations */ \
	} \
	rsi += -4*clmblkdst+2*clmblkdst8; /* Next source pointer */ \
	rdi += -4*2*clm*XMM_SCND2;	/* Restore sin/cos ptr */ \
} \
rsi += -6*2*clmblkdst8;	/* Restore source pointer */ \
end_timer(25); \
\
/* Do FFT levels 8,9 */ \
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
/* Do 48 four_complex_fft operations */ \
/*	distance between fft data elements is 1K */ \
/*	do 24 macros each processing 32 data values */ \
\
start_timer(26); \
rdi += -4*2*clm*XMM_SCD4;	/* Next premult/sin/cos ptr */ \
rcx = (uintptr_t)g->data_addr;/* Destination pointer */ \
loops_init_prefetch(24*2*clm, 64, 2, rbp); /* 24 iters of 2*clm */ \
pfing(rbp = (uintptr_t)g->premult_prefetch);/* Load premult/sin/cos prefetch pointer */ \
for(unsigned int loopA = 4; loopA; loopA--) { \
	for(unsigned int loopB = 3; loopB; loopB--) { \
		for(unsigned int loopC = 2; loopC; loopC--) { \
			for(unsigned int loopD = 2*clm; loopD; loopD--) { \
				r4_sg4cl_four_complex_fft4(rsi, 64, clmblkdst, 2*clmblkdst, rcx, 64, blkdst, 2*blkdst, rdi); \
				rdi += XMM_SCD4;		/* Next premult/sin/cos ptr */ \
			} \
			rsi += -2*clm*64+4*clmblkdst; /* Next source pointer */ \
			rcx += -2*clm*64+4*blkdst;	/* Next dest pointer */ \
			rdi += -2*clm*XMM_SCD4;	/* Restore premult/sine/cosine pointer */ \
		} \
		rsi += -2*4*clmblkdst+clmblkdst8; /* Next source pointer */ \
	} \
	rdi += 2*clm*XMM_SCD4;	/* Next premult/sine/cosine pointer */ \
} \
end_timer(26); \
\
/* Work on next 4*clm sets of 384 values */ \
\
endlp:	xpass1_get_next_block(pass2, c0b, b0b, _wpn); /* Get next block and jump */ \
/* done: */ xfft_footer_rn; \
\
\
\
/* Perform a pass 1 negacyclic FFT of length 512 doubles. */ \
/* Memory examples are for a 512K FFT with 10 levels done in pass 2. */ \
/* Sin/cos data needed: 16*2*clm*XMM_SCD2 + 4*2*clm*XMM_SCND2 + 4*2*clm*XMM_SCD4 = 32*2*clm 64-byte cache lines */ \
\

#define r4dwpn_pass1sc512ac(procname, xpass2name, pass2_levels, clm) /* (fftname, pass2_levels, clmarg) */ \
\
/* set_FFT_constants(fftname, 0, 1, pass2_levels, clmarg); */ \
xfft_header_rn(procname, xpass2name, pass2_levels, clm); \
\
/* Do pass 2 (the last 10 FFT levels if this is a 512K FFT) */ \
\
/* pass2: */start_timer(0); \
xpass2name(g); \
end_timer(0); \
if(g->ffttype == 1)		/* We're done if FFTing only */ \
	return; /* goto done; */ \
\
/* Do the inverse FFT */ \
\
pass1_inverse_fft_setup(c0b, procname); \
\
/* Do inverse FFT levels 8,9 */ \
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
/* Do 64 four_complex_unfft4 operations */ \
/*	distance between fft data elements is 1K */ \
/*	do 64 macros each processing 16 data values */ \
\
c0b:	start_timer(16); \
rsi = (uintptr_t)g->data_addr;/* Get source address */ \
rcx = (uintptr_t)g->scratch_area;	/* Get address of scratch area */ \
rdi = (uintptr_t)g->premult_addr;/* Load premult/sin/cos pointer */ \
pfing(rbp = (rdi+(32-16)*2*clm*64)); /* Finish prefetching the sin/cos/premult data */ \
loops_init_prefetch(64*2*clm, 64, 4, rbp); /* 64*2*clm iterations */ \
/*bug	xtouchpages rbp, 64*2*clm * 64 / 4 ;; Prefetch TLBs */ \
for(unsigned int loopA = 4; loopA; loopA--) { \
	for(unsigned int loopB = 4; loopB; loopB--) { \
		for(unsigned int loopC = 4; loopC; loopC--) { \
			for(unsigned int loopD = 2*clm; loopD; loopD--) { \
				r4_sg2cl_four_complex_unfft4(rsi, 64, blkdst, rcx, 64, clmblkdst, rdi); \
				rdi += XMM_SCD4;		/* Next premult/sin/cos data */ \
			} \
			rsi += -2*clm*64+2*blkdst;	/* Next source pointer */ \
			rcx += -2*clm*64+2*clmblkdst; /* Next dest pointer */ \
			rdi += -2*clm*XMM_SCD4;	/* Restore premult/sin/cos data pointer */ \
		} \
		rcx += -4*2*clmblkdst+clmblkdst8; \
	} \
	rdi += 2*clm*XMM_SCD4;	/* Next premult/sin/cos pointer */ \
} \
rsi = (uintptr_t)g->scratch_area;	/* Restore source pointer */ \
end_timer(16); \
\
/* Do inverse FFT levels 6,7 */ \
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
/* Do 64 four_complex_unfft operations */ \
/*	distance between fft data elements is 4K */ \
/*	do 32 macros each processing 32 data values */ \
\
start_timer(17); \
loops_init(32*2*clm);		/* 32*2*clm iterations */ \
r4_x4cl_wpn_four_complex_djbunfft_preload; \
for(unsigned int loopA = 16; loopA; loopA--) { \
	for(unsigned int loopB = 2; loopB; loopB--) { \
		for(unsigned int loopC = 2*clm; loopC; loopC--) { \
			r4_x4cl_wpn_four_complex_djbunfft(rsi, 64, 2*clmblkdst, 4*clmblkdst, rdi, 2*2*clm*XMM_SCND2); \
			rdi += XMM_SCND2;		/* Next clm's sin/cos data */ \
		} \
		rsi += -2*clm*64+clmblkdst; /* Next source pointer */ \
	} \
	rsi += -2*clmblkdst+clmblkdst8; /* Next source pointer */ \
	rdi += -2*2*clm*XMM_SCND2;	/* Restore sine/cosine pointer */ \
} \
rsi += -16*clmblkdst8;	/* Restore source pointer */ \
rdi += 4*2*clm*XMM_SCND2;	/* Next sin/cos data pointer */ \
end_timer(17); \
\
/* Do inverse FFT levels 4,5 */ \
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
/*	128K	... */ \
/*	... */ \
\
/* Do 64 four_complex_unfft operations */ \
/*	distance between fft data elements is 16K */ \
/*	do 32 macros each processing 32 data values */ \
\
start_timer(18); \
loops_init(32*2*clm);		/* 32 iterations of 2*clm */ \
r4_x4cl_four_complex_djbunfft_preload; \
for(unsigned int loopA = 4; loopA; loopA--) { \
	for(unsigned int loopB = 8; loopB; loopB--) { \
		for(unsigned int loopC = 2*clm; loopC; loopC--) { \
			r4_x4cl_four_complex_djbunfft(rsi, 64, clmblkdst8, 2*clmblkdst8, rdi, 8*2*clm*XMM_SCD2); \
			rdi += XMM_SCD2;		/* Next clm's sin/cos data */ \
		} \
		rsi += -2*clm*64+clmblkdst; /* Next source pointer */ \
	} \
	rsi += -8*clmblkdst+4*clmblkdst8; /* Next source pointer */ \
	rdi += -8*2*clm*XMM_SCD2;	/* Restore sine/cosine pointer */ \
} \
rsi += -4*4*clmblkdst8;	/* Restore source pointer */ \
end_timer(18); \
\
/* Do inverse FFT levels 1,2,3 */ \
/* On input the 64-byte cache lines hold these data values: */ \
/*	0K	32K	64K	96K */ \
/*	1K	... */ \
/*	... */ \
/*	31K */ \
/*	128K	... */ \
/*	... */ \
/* On output the 64-byte cache lines hold these data values: */ \
/*	0K	128K	256K	384K */ \
/*	1K	... */ \
/*	... */ \
/*	127K */ \
\
/* Do 64 four_complex_unfft operations */ \
/*	distance between fft data elements is 64K */ \
/*	do 32 macros each processing 32 data values */ \
\
start_timer(20); \
rdx = (uintptr_t)g->sincos1;/* Load fixed negacyclic premultipliers */ \
rdi = (uintptr_t)g->sincos2;/* Load fixed sin/cos postmultipliers */ \
set_data_prefetch_ptrs;		/* Init prefetching pointers */ \
loops_init(32*2*clm);		/* 32 iterations of 2*clm */ \
for(unsigned int loopA = 4; loopA; loopA--) { \
	for(unsigned int loopB = 8*2*clm; loopB; ) { \
		xtouch(xptr(rcx+rbp+(clm-1)*128));	/* Preload TLB */ \
		do{ \
			r4_x4cl_four_complex_last_djbunfft(rsi, 64, 4*clmblkdst8, 8*clmblkdst8, rdi, 32*XMM_SCD2, rdx, 32*XMM_SCD4); \
		}while(--loopB & (2*clm-1)); \
		newxloop1(clm, c1b, clmblkdst, rdi, XMM_SCD2, rdx, XMM_SCD4); /* 2*clm iterations */ \
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
copy_scratch_data_to_fft(128, clm); \
goto	endlp;			/* Skip FFT code if flag not set */ \
\
/* Fire up auxiliary threads to do pass 1 forward FFT */ \
\
/* pass1: */pass1_forward_fft_setup(b0b); \
\
/* Copy data from the FFT data area back to the scratch area */ \
\
b0b:	copy_fft_data_to_scratch(128, clm); \
\
/* Do FFT levels 1,2,3 */ \
/* */ \
/* On input the 64-byte cache lines hold these data values: */ \
/*	0K	128K	256K	384K */ \
/*	1K	... */ \
/*	... */ \
/*	127K */ \
/* On output the 64-byte cache lines hold these data values: */ \
/*	0K	64K	128K	192K */ \
/*	1K	... */ \
/*	... */ \
/*	63K	... */ \
/*	256K	... */ \
/*	... */ \
\
/* Do 64 four_complex_first_fft operations */ \
/*	distance between fft data elements is 64K */ \
/*	do 64 iterations each processing 16 data values */ \
\
yes_fft:start_timer(22); \
rsi = (uintptr_t)g->scratch_area;	/* Get address of scratch area */ \
rdx = (uintptr_t)g->sincos1;/* Address of fixed negacyclic premultipliers */ \
rdi = (uintptr_t)g->sincos2;/* Load fixed sin/cos postmultipliers */ \
set_data_prefetch_ptrs;		/* Init prefetching pointers */ \
pfing(rcx += 32*blkdst);		/* We've already prefetched 32 blocks */ \
loops_init(64*2*clm);		/* 8 iters of 4 iters of 2*2*clm */ \
for(unsigned int loopA = 8; loopA; loopA--) { \
	for(unsigned int loopB = 4*2*2*clm; loopB; ) { \
		xtouch(xptr(rcx+rbp+(clm-1)*128));	/* Preload TLB */ \
		do{ \
			do{ \
				r4_x2cl_four_complex_first_djbfft_scratch(rsi, 64, 8*clmblkdst8, rdi, rdx); \
			}while(--loopB & (2*clm-1)); \
			newxloop2(clm, b1b, clmblkdst, rdi, XMM_SCD2, rdx, XMM_SCD4); /* 2*2*clm iters prefetching 2*clm cache lines */ \
		}while(loopB & (2*2*clm-1)); \
	} \
	rsi += -4*2*clmblkdst+clmblkdst8; \
} \
rsi += -8*clmblkdst8;	/* Restore source pointer */ \
end_timer(22); \
\
/* Do FFT levels 4,5 */ \
/* */ \
/* On input the 64-byte cache lines hold these data values: */ \
/*	0K	64K	128K	192K */ \
/*	1K	... */ \
/*	... */ \
/*	63K	... */ \
/*	256K	... */ \
/*	... */ \
/* On output the 64-byte cache lines hold these data values: */ \
/*	0K	16K	32K	48K */ \
/*	1K	... */ \
/*	... */ \
/*	15K	... */ \
/*	64K	... */ \
/*	... */ \
\
/* Do 64 four_complex_fft operations */ \
/*	distance between fft data elements is 16K */ \
/*	do 32 macros each processing 32 data values */ \
\
start_timer(24); \
rdi = (uintptr_t)g->premult_addr;/* Address of premultiplier/sin/cos table */ \
rdi += 4*2*clm*XMM_SCND2 + 4*2*clm*XMM_SCD4; \
loops_init(32*2*clm);		/* 2 iters of 2 iters of 8 of 2*clm */ \
r4_x4cl_four_complex_djbfft_preload; \
for(unsigned int loopA = 2; loopA; loopA--) { \
	for(unsigned int loopB = 2; loopB; loopB--) { \
		for(unsigned int loopC = 8*2*clm; loopC; ) { \
			xtouch(xptr(rcx+rbp+(clm-1)*128));	/* Preload TLB */ \
			do{ \
				r4_x4cl_four_complex_djbfft(rsi, 64, 2*clmblkdst8, 4*clmblkdst8, rdi); \
				rdi += XMM_SCD2;		/* Next sine/cosine pointer */ \
			}while(--loopC & (2*clm-1)); \
			newxloop1(clm, b5b, clmblkdst);	/* 2*clm iterations */ \
		} \
		rsi += -8*clmblkdst+clmblkdst8; \
	} \
	rsi += -2*clmblkdst8+8*clmblkdst8; /* Next source pointer */ \
	rdi += -16*2*clm*XMM_SCD2;	/* Restore sin/cos ptr */ \
} \
rsi += -2*8*clmblkdst8;	/* Restore source pointer */ \
end_timer(24); \
\
/* Do FFT levels 6,7 */ \
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
/* Do 64 four_complex_fft operations */ \
/*	distance between fft data elements is 4K */ \
/*	do 32 macros each processing 32 data values */ \
\
start_timer(25); \
rdi += -4*2*clm*XMM_SCND2;	/* Next sin/cos ptr */ \
loops_init(32*2*clm);		/* 8 iterations of 4 iters of 2*clm */ \
r4_x4cl_wpn_four_complex_djbfft_preload; \
for(unsigned int loopA = 8; loopA; loopA--) { \
	for(unsigned int loopB = 4*2*clm; loopB; ) { \
		xtouch(xptr(rcx+rbp+(clm-1)*128));	/* Preload TLB */ \
		do{ \
			r4_x4cl_wpn_four_complex_djbfft(rsi, 64, 4*clmblkdst, clmblkdst8, rdi); \
			rdi += XMM_SCND2;		/* Next sine/cosine pointer */ \
		}while(--loopB & (2*clm-1)); \
		newxloop1(clm, b7b, clmblkdst);	/* 2*clm iterations */ \
	} \
	rsi += -4*clmblkdst+2*clmblkdst8; /* Next source pointer */ \
	rdi += -4*2*clm*XMM_SCND2;	/* Restore sin/cos ptr */ \
} \
rsi += -8*2*clmblkdst8;	/* Restore source pointer */ \
end_timer(25); \
\
/* Do FFT levels 8,9 */ \
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
/* Do 64 four_complex_fft operations */ \
/*	distance between fft data elements is 1K */ \
/*	do 32 macros each processing 32 data values */ \
\
start_timer(26); \
rdi += -4*2*clm*XMM_SCD4;	/* Next premult/sin/cos ptr */ \
rcx = (uintptr_t)g->data_addr;/* Destination pointer */ \
loops_init_prefetch(32*2*clm, 64, 2, rbp); /* 32 iters of 2*clm */ \
pfing(rbp = (uintptr_t)g->premult_prefetch);/* Load premult/sin/cos prefetch pointer */ \
for(unsigned int loopA = 4; loopA; loopA--) { \
	for(unsigned int loopB = 4; loopB; loopB--) { \
		for(unsigned int loopC = 2; loopC; loopC--) { \
			for(unsigned int loopD = 2*clm; loopD; loopD--) { \
				r4_sg4cl_four_complex_fft4(rsi, 64, clmblkdst, 2*clmblkdst, rcx, 64, blkdst, 2*blkdst, rdi); \
				rdi += XMM_SCD4;		/* Next premult/sin/cos ptr */ \
			} \
			rsi += -2*clm*64+4*clmblkdst; /* Next source pointer */ \
			rcx += -2*clm*64+4*blkdst;	/* Next dest pointer */ \
			rdi += -2*clm*XMM_SCD4;	/* Restore premult/sine/cosine pointer */ \
		} \
		rsi += -2*4*clmblkdst+clmblkdst8; /* Next source pointer */ \
	} \
	rdi += 2*clm*XMM_SCD4;	/* Next premult/sine/cosine pointer */ \
} \
end_timer(26); \
\
/* Work on next 4*clm sets of 512 values */ \
\
endlp:	xpass1_get_next_block(pass2, c0b, b0b, _wpn); /* Get next block and jump */ \
/* done: */ xfft_footer_rn; \
\
\
\
/* Perform a pass 1 negacyclic FFT of length 640 doubles. */ \
/* Memory examples are for a 640K FFT with 10 levels done in pass 2. */ \
/* Sin/cos data needed: 16*2*clm*XMM_SCD2 + 4*2*clm*XMM_SCND2 + 4*2*clm*XMM_SCD4 = 32*2*clm 64-byte cache lines */ \
\

#define r4dwpn_pass1sc640ac(procname, xpass2name, pass2_levels, clm) /* (fftname, pass2_levels, clmarg) */ \
\
/* set_FFT_constants(fftname, 0, 1, pass2_levels, clmarg); */ \
xfft_header_rn(procname, xpass2name, pass2_levels, clm); \
\
/* Do pass 2 (the last 10 FFT levels if this is a 640K FFT) */ \
\
/* pass2: */start_timer(0); \
xpass2name(g); \
end_timer(0); \
if(g->ffttype == 1)		/* We're done if FFTing only */ \
	return; /* goto done; */ \
\
/* Do the inverse FFT */ \
\
pass1_inverse_fft_setup(c0b, procname); \
\
/* Do inverse FFT levels 8,9 */ \
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
/* Do 80 four_complex_unfft4 operations */ \
/*	distance between fft data elements is 1K */ \
/*	do 80 macros each processing 16 data values */ \
\
c0b:	start_timer(16); \
rsi = (uintptr_t)g->data_addr;/* Get source address */ \
rcx = (uintptr_t)g->scratch_area;	/* Get address of scratch area */ \
rdi = (uintptr_t)g->premult_addr;/* Load premult/sin/cos pointer */ \
pfing(rbp = (rdi+(32-20)*2*clm*64)); /* Finish prefetching the sin/cos/premult data */ \
loops_init_prefetch(80*2*clm, 64, 4, rbp); /* 80*2*clm iterations */ \
/*bug	xtouchpages rbp, 80*2*clm * 64 / 4 ;; Prefetch TLBs */ \
for(unsigned int loopA = 4; loopA; loopA--) { \
	for(unsigned int loopB = 5; loopB; loopB--) { \
		for(unsigned int loopC = 4; loopC; loopC--) { \
			for(unsigned int loopD = 2*clm; loopD; loopD--) { \
				r4_sg2cl_four_complex_unfft4(rsi, 64, blkdst, rcx, 64, clmblkdst, rdi); \
				rdi += XMM_SCD4;		/* Next premult/sin/cos data */ \
			} \
			rsi += -2*clm*64+2*blkdst;	/* Next source pointer */ \
			rcx += -2*clm*64+2*clmblkdst; /* Next dest pointer */ \
			rdi += -2*clm*XMM_SCD4;	/* Restore premult/sin/cos data pointer */ \
		} \
		rcx += -4*2*clmblkdst+clmblkdst8; \
	} \
	rdi += 2*clm*XMM_SCD4;	/* Next premult/sin/cos pointer */ \
} \
rsi = (uintptr_t)g->scratch_area;	/* Restore source pointer */ \
end_timer(16); \
\
/* Do inverse FFT levels 6,7 */ \
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
/* Do 80 four_complex_unfft operations */ \
/*	distance between fft data elements is 4K */ \
/*	do 40 macros each processing 32 data values */ \
\
start_timer(17); \
loops_init(40*2*clm);		/* 40*2*clm iterations */ \
r4_x4cl_wpn_four_complex_djbunfft_preload; \
for(unsigned int loopA = 20; loopA; loopA--) { \
	for(unsigned int loopB = 2; loopB; loopB--) { \
		for(unsigned int loopC = 2*clm; loopC; loopC--) { \
			r4_x4cl_wpn_four_complex_djbunfft(rsi, 64, 2*clmblkdst, 4*clmblkdst, rdi, 2*2*clm*XMM_SCND2); \
			rdi += XMM_SCND2;		/* Next clm's sin/cos data */ \
		} \
		rsi += -2*clm*64+clmblkdst; /* Next source pointer */ \
	} \
	rsi += -2*clmblkdst+clmblkdst8; /* Next source pointer */ \
	rdi += -2*2*clm*XMM_SCND2;	/* Restore sine/cosine pointer */ \
} \
rsi += -20*clmblkdst8;	/* Restore source pointer */ \
rdi += 4*2*clm*XMM_SCND2;	/* Next sin/cos data pointer */ \
end_timer(17); \
\
/* Do inverse FFT levels 4,5 (actually 2.322 levels) */ \
/* On input the 64-byte cache lines hold these data values: */ \
/*	0K	8K	16K	24K */ \
/*	1K	... */ \
/*	... */ \
/*	7K */ \
/*	32K	... */ \
/*	... */ \
/* On output the 64-byte cache lines hold these data values: */ \
/*	0K	40K	80K	120K */ \
/*	1K	... */ \
/*	... */ \
/*	39K */ \
/*	160K	... */ \
/*	... */ \
\
/* Do 64 five_complex_unfft operations */ \
/*	distance between fft data elements is 16K */ \
/*	do 32 macros each processing 40 data values */ \
\
start_timer(18); \
loops_init(32*2*clm);		/* 32 iterations of 2*clm */ \
r5_x5cl_five_complex_djbunfft_preload; \
for(unsigned int loopA = 4; loopA; loopA--) { \
	for(unsigned int loopB = 8; loopB; loopB--) { \
		for(unsigned int loopC = 2*clm; loopC; loopC--) { \
			r5_x5cl_five_complex_djbunfft(rsi, 64, clmblkdst8, rdi, 8*2*clm*XMM_SCD2); \
			rdi += XMM_SCD2;		/* Next clm's sin/cos data */ \
		} \
		rsi += -2*clm*64+clmblkdst; /* Next source pointer */ \
	} \
	rsi += -8*clmblkdst+5*clmblkdst8; /* Next source pointer */ \
	rdi += -8*2*clm*XMM_SCD2;	/* Restore sine/cosine pointer */ \
} \
rsi += -4*5*clmblkdst8;	/* Restore source pointer */ \
end_timer(18); \
\
/* Do inverse FFT levels 1,2,3 */ \
/* On input the 64-byte cache lines hold these data values: */ \
/*	0K	40K	80K	120K */ \
/*	1K	... */ \
/*	... */ \
/*	39K */ \
/*	160K	... */ \
/*	... */ \
/* On output the 64-byte cache lines hold these data values: */ \
/*	0K	160K	320K	480K */ \
/*	1K	... */ \
/*	... */ \
/*	159K */ \
\
/* Do 80 four_complex_unfft operations */ \
/*	distance between fft data elements is 80K */ \
/*	do 40 macros each processing 32 data values */ \
\
start_timer(20); \
rdx = (uintptr_t)g->sincos1;/* Load fixed negacyclic premultipliers */ \
rdi = (uintptr_t)g->sincos2;/* Load fixed sin/cos postmultipliers */ \
set_data_prefetch_ptrs;		/* Init prefetching pointers */ \
loops_init(40*2*clm);		/* 40 iterations of 2*clm */ \
for(unsigned int loopA = 5; loopA; loopA--) { \
	for(unsigned int loopB = 4*2*2*clm; loopB; ) { \
		xtouch(xptr(rcx+rbp+(clm-1)*128));	/* Preload TLB */ \
		do{ \
			do{ \
				r4_x4cl_four_complex_last_djbunfft(rsi, 64, 5*clmblkdst8, 10*clmblkdst8, rdi, 40*XMM_SCD2, rdx, 40*XMM_SCD4); \
			}while(--loopB & (2*clm-1)); \
			newxloop2(clm, c1b, clmblkdst, rdi, XMM_SCD2, rdx, XMM_SCD4); /* 2*2*clm iterations */ \
		}while(loopB & (2*2*clm-1)); \
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
copy_scratch_data_to_fft(160, clm); \
goto	endlp;			/* Skip FFT code if flag not set */ \
\
/* Fire up auxiliary threads to do pass 1 forward FFT */ \
\
/* pass1: */pass1_forward_fft_setup(b0b); \
\
/* Copy data from the FFT data area back to the scratch area */ \
\
b0b:	copy_fft_data_to_scratch(160, clm); \
\
/* Do FFT levels 1,2,3 */ \
/* */ \
/* On input the 64-byte cache lines hold these data values: */ \
/*	0K	160K	320K	480K */ \
/*	1K	... */ \
/*	... */ \
/*	159K */ \
/* On output the 64-byte cache lines hold these data values: */ \
/*	0K	80K	160K	240K */ \
/*	1K	... */ \
/*	... */ \
/*	79K	... */ \
/*	320K	... */ \
/*	... */ \
\
/* Do 80 four_complex_first_fft operations */ \
/*	distance between fft data elements is 80K */ \
/*	do 80 iterations each processing 16 data values */ \
\
yes_fft:start_timer(22); \
rsi = (uintptr_t)g->scratch_area;	/* Get address of scratch area */ \
rdx = (uintptr_t)g->sincos1;/* Address of fixed negacyclic premultipliers */ \
rdi = (uintptr_t)g->sincos2;/* Load fixed sin/cos postmultipliers */ \
set_data_prefetch_ptrs;		/* Init prefetching pointers */ \
pfing(rcx += (20-4)*blkdst);	/* We've already prefetched 20 blocks, but */ \
;					/* code below prefetches 4 blocks too many */ \
loops_init(80*2*clm);		/* 10 iters of 4 iters of 2*2*clm */ \
for(unsigned int loopA = 10; loopA; loopA--) { \
	for(unsigned int loopB = 4*2*2*clm; loopB; ) { \
		xtouch(xptr(rcx+rbp+(clm-1)*128));	/* Preload TLB */ \
		do{ \
			do{ \
				r4_x2cl_four_complex_first_djbfft_scratch(rsi, 64, 10*clmblkdst8, rdi, rdx); \
			}while(--loopB & (2*clm-1)); \
			newxloop2(clm, b1b, clmblkdst, rdi, XMM_SCD2, rdx, XMM_SCD4); /* 2*2*clm iters prefetching 2*clm cache lines */ \
		}while(loopB & (2*2*clm-1)); \
	} \
	rsi += -4*2*clmblkdst+clmblkdst8; \
} \
rsi += -10*clmblkdst8;	/* Restore source pointer */ \
end_timer(22); \
\
/* Do FFT levels 4,5 (actually 2.322 levels) */ \
/* */ \
/* On input the 64-byte cache lines hold these data values: */ \
/*	0K	80K	160K	240K */ \
/*	1K	... */ \
/*	... */ \
/*	79K	... */ \
/*	320K	... */ \
/*	... */ \
/* On output the 64-byte cache lines hold these data values: */ \
/*	0K	16K	32K	48K */ \
/*	1K	... */ \
/*	... */ \
/*	15K	... */ \
/*	64K	... */ \
/*	... */ \
\
/* Do 64 five_complex_fft operations */ \
/*	distance between fft data elements is 16K */ \
/*	do 32 macros each processing 40 data values */ \
\
start_timer(24); \
rdi = (uintptr_t)g->premult_addr;/* Address of premultiplier/sin/cos table */ \
rdi += 4*2*clm*XMM_SCND2 + 4*2*clm*XMM_SCD4; \
loops_init(32*2*clm);		/* 2 iters of 2 iters of 8 of 2*clm */ \
r5_x5cl_five_complex_djbfft_preload; \
for(unsigned int loopA = 2; loopA; loopA--) { \
	for(unsigned int loopB = 2; loopB; loopB--) { \
		for(unsigned int loopC = 8*2*clm; loopC; ) { \
			xtouch(xptr(rcx+rbp+(clm-1)*128));	/* Preload TLB */ \
			do{ \
				r5_x5cl_five_complex_djbfft(rsi, 64, 2*clmblkdst8, rdi); \
				rdi += XMM_SCD2;		/* Next sine/cosine pointer */ \
			}while(--loopC & (2*clm-1)); \
			newxloop1doublepre(clm, b5b, b5a, clmblkdst); /* 2*clm iterations */ \
		} \
		rsi += -8*clmblkdst+clmblkdst8; \
	} \
	rsi += -2*clmblkdst8+10*clmblkdst8; /* Next source pointer */ \
	rdi += -16*2*clm*XMM_SCD2;	/* Restore sin/cos ptr */ \
} \
rsi += -2*10*clmblkdst8;	/* Restore source pointer */ \
end_timer(24); \
\
/* Do FFT levels 6,7 */ \
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
/* Do 80 four_complex_fft operations */ \
/*	distance between fft data elements is 4K */ \
/*	do 40 macros each processing 32 data values */ \
\
start_timer(25); \
rdi += -4*2*clm*XMM_SCND2;	/* Next sin/cos ptr */ \
loops_init(40*2*clm);		/* 10 iterations of 4 iters of 2*clm */ \
r4_x4cl_wpn_four_complex_djbfft_preload; \
for(unsigned int loopA = 10; loopA; loopA--) { \
	for(unsigned int loopB = 4*2*clm; loopB; ) { \
		xtouch(xptr(rcx+rbp+(clm-1)*128));	/* Preload TLB */ \
		do{ \
			r4_x4cl_wpn_four_complex_djbfft(rsi, 64, 4*clmblkdst, clmblkdst8, rdi); \
			rdi += XMM_SCND2;		/* Next sine/cosine pointer */ \
		}while(--loopB & (2*clm-1)); \
		newxloop1(clm, b7b, clmblkdst);	/* 2*clm iterations */ \
	} \
	rsi += -4*clmblkdst+2*clmblkdst8; /* Next source pointer */ \
	rdi += -4*2*clm*XMM_SCND2;	/* Restore sin/cos ptr */ \
} \
rsi += -10*2*clmblkdst8;	/* Restore source pointer */ \
end_timer(25); \
\
/* Do FFT levels 8,9 */ \
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
/* Do 80 four_complex_fft operations */ \
/*	distance between fft data elements is 1K */ \
/*	do 40 macros each processing 32 data values */ \
\
start_timer(26); \
rdi += -4*2*clm*XMM_SCD4;	/* Next premult/sin/cos ptr */ \
rcx = (uintptr_t)g->data_addr;/* Destination pointer */ \
loops_init_prefetch(40*2*clm, 64, 2, rbp); /* 40 iters of 2*clm */ \
pfing(rbp = (uintptr_t)g->premult_prefetch);/* Load premult/sin/cos prefetch pointer */ \
for(unsigned int loopA = 4; loopA; loopA--) { \
	for(unsigned int loopB = 5; loopB; loopB--) { \
		for(unsigned int loopC = 2; loopC; loopC--) { \
			for(unsigned int loopD = 2*clm; loopD; loopD--) { \
				r4_sg4cl_four_complex_fft4(rsi, 64, clmblkdst, 2*clmblkdst, rcx, 64, blkdst, 2*blkdst, rdi); \
				rdi += XMM_SCD4;		/* Next premult/sin/cos ptr */ \
			} \
			rsi += -2*clm*64+4*clmblkdst; /* Next source pointer */ \
			rcx += -2*clm*64+4*blkdst;	/* Next dest pointer */ \
			rdi += -2*clm*XMM_SCD4;	/* Restore premult/sine/cosine pointer */ \
		} \
		rsi += -2*4*clmblkdst+clmblkdst8; /* Next source pointer */ \
	} \
	rdi += 2*clm*XMM_SCD4;	/* Next premult/sine/cosine pointer */ \
} \
end_timer(26); \
\
/* Work on next 4*clm sets of 640 values */ \
\
endlp:	xpass1_get_next_block(pass2, c0b, b0b, _wpn); /* Get next block and jump */ \
/* done: */ xfft_footer_rn; \
\
\
\
/* Perform a pass 1 negacyclic FFT of length 768 doubles. */ \
/* Memory examples are for a 768K FFT with 10 levels done in pass 2. */ \
/* Sin/cos data needed: 32*2*clm*XMM_SCD1 + 8*2*clm*XMM_SCND2 + 4*2*clm*XMM_SCD8 = 48*2*clm 64-byte cache lines */ \
\

#define r4dwpn_pass1sc768ac(procname, xpass2name, pass2_levels, clm) /* (fftname, pass2_levels, clmarg) */ \
\
/* set_FFT_constants(fftname, 0, 1, pass2_levels, clmarg); */ \
xfft_header_rn(procname, xpass2name, pass2_levels, clm); \
\
/* Do pass 2 (the last 10 FFT levels if this is a 768K FFT) */ \
\
/* pass2: */start_timer(0); \
xpass2name(g); \
end_timer(0); \
if(g->ffttype == 1)		/* We're done if FFTing only */ \
	return; /* goto done; */ \
\
/* Do the inverse FFT */ \
\
pass1_inverse_fft_setup(c0b, procname); \
\
/* Do inverse FFT levels 8,9,10 */ \
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
/* Do 48 eight_complex_unfft8 operations */ \
/*	distance between fft data elements is 1K */ \
/*	do 48 macros each processing 32 data values */ \
\
c0b:	start_timer(16); \
rsi = (uintptr_t)g->data_addr;/* Get source address */ \
rcx = (uintptr_t)g->scratch_area;	/* Get address of scratch area */ \
rdi = (uintptr_t)g->premult_addr;/* Load premult/sin/cos pointer */ \
pfing(rbp = (rdi+(48-24)*2*clm*64)); /* Finish prefetching the sin/cos/premult data */ \
loops_init_prefetch(48*2*clm, 64, 2, rbp); /* 48*2*clm iterations */ \
/*bug	xtouchpages rbp, 48*2*clm * 64 / 2 ;; Prefetch TLBs */ \
rdi += 128;		/* Generates 1 byte offsets in radix-8 macro */ \
for(unsigned int loopA = 4; loopA; loopA--) { \
	for(unsigned int loopB = 6; loopB; loopB--) { \
		for(unsigned int loopC = 2; loopC; loopC--) { \
			for(unsigned int loopD = 2*clm; loopD; loopD--) { \
				r8_sg4cl_eight_complex_unfft8(rsi, 64, blkdst, 2*blkdst, rcx, 64, clmblkdst, 2*clmblkdst, rdi-128); \
				rdi += XMM_SCD8;		/* Next premult/sin/cos data */ \
			} \
			rsi += -2*clm*64+4*blkdst;	/* Next source pointer */ \
			rcx += -2*clm*64+4*clmblkdst; /* Next dest pointer */ \
			rdi += -2*clm*XMM_SCD8;	/* Restore premult/sin/cos data pointer */ \
		} \
		rcx += -2*4*clmblkdst+clmblkdst8; \
	} \
	rdi += 2*clm*XMM_SCD8;	/* Next premult/sin/cos pointer */ \
} \
rsi = (uintptr_t)g->scratch_area;	/* Restore source pointer */ \
rdi += -128;		/* Next sin/cos data pointer */ \
end_timer(16); \
\
/* Do inverse FFT levels 6,7 */ \
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
/* Do 96 four_complex_unfft operations */ \
/*	distance between fft data elements is 8K */ \
/*	do 48 macros each processing 32 data values */ \
\
start_timer(17); \
loops_init(48*2*clm);		/* 48 iterations of 2*clm */ \
r4_x4cl_wpn_four_complex_djbunfft_preload; \
for(unsigned int loopA = 12; loopA; loopA--) { \
	for(unsigned int loopB = 4; loopB; loopB--) { \
		for(unsigned int loopC = 2*clm; loopC; loopC--) { \
			r4_x4cl_wpn_four_complex_djbunfft(rsi, 64, 4*clmblkdst, clmblkdst8, rdi, 4*2*clm*XMM_SCND2); \
			rdi += XMM_SCND2;		/* Next clm's sin/cos data */ \
		} \
		rsi += -2*clm*64+clmblkdst; /* Next source pointer */ \
	} \
	rsi += -4*clmblkdst+2*clmblkdst8; /* Next source pointer */ \
	rdi += -4*2*clm*XMM_SCND2;	/* Restore sine/cosine pointer */ \
} \
rsi += -12*2*clmblkdst8;	/* Restore source pointer */ \
rdi += 8*2*clm*XMM_SCND2;	/* Next sin/cos data pointer */ \
end_timer(17); \
\
/* Do inverse FFT levels 4,5 (actually 1.585 levels) */ \
/* On input the 64-byte cache lines hold these data values: */ \
/*	0K	16K	32K	48K */ \
/*	1K	... */ \
/*	... */ \
/*	15K */ \
/*	64K	... */ \
/*	... */ \
/* On output the 64-byte cache lines hold these data values: */ \
/*	0K	48K	96K	144K */ \
/*	1K	... */ \
/*	... */ \
/*	47K */ \
/*	192K	... */ \
/*	... */ \
\
/* Do 128 three_complex_unfft operations */ \
/*	distance between fft data elements is 32K */ \
/*	do 64 macros each processing 24 data values */ \
\
start_timer(18); \
loops_init(64*2*clm);		/* 64 iterations of 2*clm */ \
r3_x3cl_three_complex_djbunfft_preload; \
for(unsigned int loopA = 4; loopA; loopA--) { \
	for(unsigned int loopB = 2; loopB; loopB--) { \
		for(unsigned int loopC = 8; loopC; loopC--) { \
			for(unsigned int loopD = 2*clm; loopD; loopD--) { \
				r3_x3cl_three_complex_djbunfft(rsi, 64, 2*clmblkdst8, rdi, 16*2*clm*XMM_SCD1); \
				rdi += XMM_SCD1;		/* Next clm's sin/cos data */ \
			} \
			rsi += -2*clm*64+clmblkdst; /* Next source pointer */ \
		} \
		rsi += -8*clmblkdst+clmblkdst8; /* Next source pointer */ \
	} \
	rsi += -2*clmblkdst8+6*clmblkdst8; /* Next source pointer */ \
	rdi += -2*8*2*clm*XMM_SCD1; /* Restore sine/cosine pointer */ \
} \
rsi += -4*6*clmblkdst8;	/* Restore source pointer */ \
end_timer(18); \
\
/* Do inverse FFT levels 1,2,3 */ \
/* On input the 64-byte cache lines hold these data values: */ \
/*	0K	48K	96K	144K */ \
/*	1K	... */ \
/*	... */ \
/*	47K */ \
/*	192K	... */ \
/*	... */ \
/* On output the 64-byte cache lines hold these data values: */ \
/*	0K	192K	384K	576K */ \
/*	1K	... */ \
/*	... */ \
/*	191K */ \
\
/* Do 96 four_complex_unfft operations */ \
/*	distance between fft data elements is 96K */ \
/*	do 48 macros each processing 32 data values */ \
\
start_timer(20); \
rdx = (uintptr_t)g->sincos1;/* Load fixed negacyclic premultipliers */ \
rdi = (uintptr_t)g->sincos2;/* Load fixed sin/cos postmultipliers */ \
set_data_prefetch_ptrs;		/* Init prefetching pointers */ \
loops_init(48*2*clm);		/* 48 iterations of 2*clm */ \
for(unsigned int loopA = 6; loopA; loopA--) { \
	for(unsigned int loopB = 8*2*clm; loopB; ) { \
		xtouch(xptr(rcx+rbp+(clm-1)*128));	/* Preload TLB */ \
		do{ \
			r4_x4cl_four_complex_last_djbunfft(rsi, 64, 6*clmblkdst8, 12*clmblkdst8, rdi, 48*XMM_SCD2, rdx, 48*XMM_SCD4); \
		}while(--loopB & (2*clm-1)); \
		newxloop1(clm, c1b, clmblkdst, rdi, XMM_SCD2, rdx, XMM_SCD4); /* 2*clm iterations */ \
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
copy_scratch_data_to_fft(192, clm); \
goto	endlp;			/* Skip FFT code if flag not set */ \
\
/* Fire up auxiliary threads to do pass 1 forward FFT */ \
\
/* pass1: */pass1_forward_fft_setup(b0b); \
\
/* Copy data from the FFT data area back to the scratch area */ \
\
b0b:	copy_fft_data_to_scratch(192, clm); \
\
/* Do FFT levels 1,2,3 */ \
/* */ \
/* On input the 64-byte cache lines hold these data values: */ \
/*	0K	192K	384K	576K */ \
/*	1K	... */ \
/*	... */ \
/*	191K */ \
/* On output the 64-byte cache lines hold these data values: */ \
/*	0K	96K	192K	288K */ \
/*	1K	... */ \
/*	... */ \
/*	95K	... */ \
/*	384K	... */ \
/*	... */ \
\
/* Do 96 four_complex_first_fft operations */ \
/*	distance between fft data elements is 96K */ \
/*	do 96 iterations each processing 16 data values */ \
\
yes_fft:start_timer(22); \
rsi = (uintptr_t)g->scratch_area;	/* Get address of scratch area */ \
rdx = (uintptr_t)g->sincos1;/* Address of fixed negacyclic premultipliers */ \
rdi = (uintptr_t)g->sincos2;/* Load fixed sin/cos postmultipliers */ \
set_data_prefetch_ptrs;		/* Init prefetching pointers */ \
pfing(rcx += (48-16)*blkdst);	/* We've already prefetched 48 blocks, but */ \
;					/* code below prefetches 16 blocks too many */ \
loops_init(96*2*clm);		/* 12 iters of 4 iters of 2*2*clm */ \
for(unsigned int loopA = 12; loopA; loopA--) { \
	for(unsigned int loopB = 4*2*2*clm; loopB; ) { \
		xtouch(xptr(rcx+rbp+(clm-1)*128));	/* Preload TLB */ \
		do{ \
			do{ \
				r4_x2cl_four_complex_first_djbfft_scratch(rsi, 64, 12*clmblkdst8, rdi, rdx); \
			}while(--loopB & (2*clm-1)); \
			newxloop2(clm, b1b, clmblkdst, rdi, XMM_SCD2, rdx, XMM_SCD4); /* 2*2*clm iters prefetching 2*clm cache lines */ \
		}while(loopB & (2*2*clm-1)); \
	} \
	rsi += -4*2*clmblkdst+clmblkdst8; \
} \
rsi += -12*clmblkdst8;	/* Restore source pointer */ \
end_timer(22); \
\
/* Do FFT levels 4,5 (actually 1.585 levels) */ \
/* */ \
/* On input the 64-byte cache lines hold these data values: */ \
/*	0K	96K	192K	288K */ \
/*	1K	... */ \
/*	... */ \
/*	95K	... */ \
/*	384K	... */ \
/*	... */ \
/* On output the 64-byte cache lines hold these data values: */ \
/*	0K	32K	64K	96K */ \
/*	1K	... */ \
/*	... */ \
/*	31K	... */ \
/*	128K	... */ \
/*	... */ \
\
/* Do 128 three_complex_fft operations */ \
/*	distance between fft data elements is 32K */ \
/*	do 64 macros each processing 24 data values */ \
\
start_timer(24); \
rdi = (uintptr_t)g->premult_addr;/* Address of premultiplier/sin/cos table */ \
rdi += 8*2*clm*XMM_SCND2 + 4*2*clm*XMM_SCD8; \
loops_init(64*2*clm);		/* 2 iters of 4 iters of 8 of 2*clm */ \
r3_x3cl_three_complex_djbfft_preload; \
for(unsigned int loopA = 2; loopA; loopA--) { \
	for(unsigned int loopB = 4; loopB; loopB--) { \
		for(unsigned int loopC = 8*2*clm; loopC; ) { \
			xtouch(xptr(rcx+rbp+(clm-1)*128));	/* Preload TLB */ \
			do{ \
				r3_x3cl_three_complex_djbfft(rsi, 64, 4*clmblkdst8, rdi); \
				rdi += XMM_SCD1;		/* Next sine/cosine pointer */ \
			}while(--loopC & (2*clm-1)); \
			newxloop1(clm, b5b, clmblkdst);	/* 2*clm iterations */ \
		} \
		rsi += -8*clmblkdst+clmblkdst8; \
	} \
	rsi += -4*clmblkdst8+12*clmblkdst8; /* Next source pointer */ \
	rdi += -32*2*clm*XMM_SCD1;	/* Restore sin/cos ptr */ \
} \
rsi += -2*12*clmblkdst8;	/* Restore source pointer */ \
end_timer(24); \
\
/* Do FFT levels 6,7 */ \
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
/* Do 96 four_complex_fft operations */ \
/*	distance between fft data elements is 8K */ \
/*	do 48 macros each processing 32 data values */ \
\
start_timer(25); \
rdi += -8*2*clm*XMM_SCND2;	/* Next sin/cos ptr */ \
loops_init(48*2*clm);		/* 6 iterations of 8 iters of 2*clm */ \
r4_x4cl_wpn_four_complex_djbfft_preload; \
for(unsigned int loopA = 6; loopA; loopA--) { \
	for(unsigned int loopB = 8*2*clm; loopB; ) { \
		xtouch(xptr(rcx+rbp+(clm-1)*128));	/* Preload TLB */ \
		do{ \
			r4_x4cl_wpn_four_complex_djbfft(rsi, 64, clmblkdst8, 2*clmblkdst8, rdi); \
			rdi += XMM_SCND2;		/* Next sine/cosine pointer */ \
		}while(--loopB & (2*clm-1)); \
		newxloop1(clm, b7b, clmblkdst);	/* 2*clm iterations */ \
	} \
	rsi += -8*clmblkdst+4*clmblkdst8; /* Next source pointer */ \
	rdi += -8*2*clm*XMM_SCND2;	/* Restore sin/cos ptr */ \
} \
rsi += -6*4*clmblkdst8;	/* Restore source pointer */ \
end_timer(25); \
\
/* Do FFT levels 8,9,10 */ \
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
/* Do 48 eight_complex_fft operations */ \
/*	distance between fft data elements is 1K */ \
/*	do 24 macros each processing 64 data values */ \
\
start_timer(26); \
rdi += -4*2*clm*XMM_SCD8+128; /* Next premult/sin/cos ptr (gen 1 byte offsets in radix-8 macro) */ \
rcx = (uintptr_t)g->data_addr;/* Destination pointer */ \
loops_init_prefetch(24*2*clm, 64, 1, rbp); /* 24 iters of 2*clm */ \
pfing(rbp = (uintptr_t)g->premult_prefetch);/* Load premult/sin/cos prefetch pointer */ \
for(unsigned int loopA = 4; loopA; loopA--) { \
	for(unsigned int loopB = 6; loopB; loopB--) { \
		for(unsigned int loopC = 2*clm; loopC; loopC--) { \
			r8_sg8cl_eight_complex_fft8(rsi, 64, clmblkdst, 2*clmblkdst, 4*clmblkdst, rcx, 64, blkdst, 2*blkdst, 4*blkdst, rdi-128); \
			rdi += XMM_SCD8;		/* Next premult/sin/cos ptr */ \
		} \
		rsi += -2*clm*64+clmblkdst8; /* Next source pointer */ \
		rcx += -2*clm*64+8*blkdst;	/* Next dest pointer */ \
		rdi += -2*clm*XMM_SCD8;	/* Restore premult/sine/cosine pointer */ \
	} \
	rdi += 2*clm*XMM_SCD8;	/* Next premult/sine/cosine pointer */ \
} \
end_timer(26); \
\
/* Work on next 4*clm sets of 768 values */ \
\
endlp:	xpass1_get_next_block(pass2, c0b, b0b, _wpn); /* Get next block and jump */ \
/* done: */ xfft_footer_rn; \
\
\
\
/* Perform a pass 1 negacyclic FFT of length 1024 doubles. */ \
/* Memory examples are for a 1M FFT with 10 levels done in pass 2. */ \
/* Prefetch 32*2*clm*XMM_SCD2 + 8*2*clm*XMM_SCND2 + 4*2*clm*XMM_SCD8 = 64*2*clm 64-byte cache lines */ \
\

#define r4dwpn_pass1sc1024ac(procname, xpass2name, pass2_levels, clm) /* (fftname, pass2_levels, clmarg) */ \
\
/* set_FFT_constants(fftname, 0, 1, pass2_levels, clmarg); */ \
xfft_header_rn(procname, xpass2name, pass2_levels, clm); \
\
/* Do pass 2 (the last 10 FFT levels if this is a 1M FFT) */ \
\
/* pass2: */start_timer(0); \
xpass2name(g); \
end_timer(0); \
if(g->ffttype == 1)		/* We're done if FFTing only */ \
	return; /* goto done; */ \
\
/* Do the inverse FFT */ \
\
pass1_inverse_fft_setup(c0b, procname); \
\
/* Do inverse FFT levels 8,9,10 */ \
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
/* Do 64 eight_complex_unfft8 operations */ \
/*	distance between fft data elements is 1K */ \
/*	do 64 macros each processing 32 data values */ \
\
c0b:	start_timer(16); \
rsi = (uintptr_t)g->data_addr;/* Get source address */ \
rcx = (uintptr_t)g->scratch_area;	/* Get address of scratch area */ \
rdi = (uintptr_t)g->premult_addr;/* Load premult/sin/cos pointer */ \
pfing(rbp = (rdi+(64-32)*2*clm*64)); /* Finish prefetching the sin/cos/premult data */ \
loops_init_prefetch(64*2*clm, 64, 2, rbp); /* 64*2*clm iterations */ \
/*bug	xtouchpages rbp, 64*2*clm * 64 / 2 ;; Prefetch TLBs */ \
rdi += 128;		/* Generates 1 byte offsets in radix-8 macro */ \
for(unsigned int loopA = 4; loopA; loopA--) { \
	for(unsigned int loopB = 8; loopB; loopB--) { \
		for(unsigned int loopC = 2; loopC; loopC--) { \
			for(unsigned int loopD = 2*clm; loopD; loopD--) { \
				r8_sg4cl_eight_complex_unfft8(rsi, 64, blkdst, 2*blkdst, rcx, 64, clmblkdst, 2*clmblkdst, rdi-128); \
				rdi += XMM_SCD8;		/* Next premult/sin/cos data */ \
			} \
			rsi += -2*clm*64+4*blkdst;	/* Next source pointer */ \
			rcx += -2*clm*64+4*clmblkdst; /* Next dest pointer */ \
			rdi += -2*clm*XMM_SCD8;	/* Restore premult/sin/cos data pointer */ \
		} \
		rcx += -2*4*clmblkdst+clmblkdst8; \
	} \
	rdi += 2*clm*XMM_SCD8;	/* Next premult/sin/cos pointer */ \
} \
rsi = (uintptr_t)g->scratch_area;	/* Restore source pointer */ \
rdi += -128;		/* Next sin/cos data pointer */ \
end_timer(16); \
\
/* Do inverse FFT levels 6,7 */ \
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
/* Do 128 four_complex_unfft operations */ \
/*	distance between fft data elements is 8K */ \
/*	do 64 macros each processing 32 data values */ \
\
start_timer(17); \
loops_init(64*2*clm);		/* 64 iterations of 2*clm */ \
r4_x4cl_wpn_four_complex_djbunfft_preload; \
for(unsigned int loopA = 16; loopA; loopA--) { \
	for(unsigned int loopB = 4; loopB; loopB--) { \
		for(unsigned int loopC = 2*clm; loopC; loopC--) { \
			r4_x4cl_wpn_four_complex_djbunfft(rsi, 64, 4*clmblkdst, clmblkdst8, rdi, 4*2*clm*XMM_SCND2); \
			rdi += XMM_SCND2;		/* Next clm's sin/cos data */ \
		} \
		rsi += -2*clm*64+clmblkdst; /* Next source pointer */ \
	} \
	rsi += -4*clmblkdst+2*clmblkdst8; /* Next source pointer */ \
	rdi += -4*2*clm*XMM_SCND2;	/* Restore sine/cosine pointer */ \
} \
rsi += -16*2*clmblkdst8;	/* Restore source pointer */ \
rdi += 8*2*clm*XMM_SCND2;	/* Next sin/cos data pointer */ \
end_timer(17); \
\
/* Do inverse FFT levels 4,5 */ \
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
/*	256K	... */ \
/*	... */ \
\
/* Do 128 four_complex_unfft operations */ \
/*	distance between fft data elements is 32K */ \
/*	do 64 macros each processing 32 data values */ \
\
start_timer(18); \
loops_init(64*2*clm);		/* 64 iterations of 2*clm */ \
r4_x4cl_four_complex_djbunfft_preload; \
for(unsigned int loopA = 4; loopA; loopA--) { \
	for(unsigned int loopB = 2; loopB; loopB--) { \
		for(unsigned int loopC = 8; loopC; loopC--) { \
			for(unsigned int loopD = 2*clm; loopD; loopD--) { \
				r4_x4cl_four_complex_djbunfft(rsi, 64, 2*clmblkdst8, 4*clmblkdst8, rdi, 16*2*clm*XMM_SCD2); \
				rdi += XMM_SCD2;		/* Next clm's sin/cos data */ \
			} \
			rsi += -2*clm*64+clmblkdst; /* Next source pointer */ \
		} \
		rsi += -8*clmblkdst+clmblkdst8; /* Next source pointer */ \
	} \
	rsi += -2*clmblkdst8+8*clmblkdst8; /* Next source pointer */ \
	rdi += -2*8*2*clm*XMM_SCD2; /* Restore sine/cosine pointer */ \
} \
rsi += -4*8*clmblkdst8;	/* Restore source pointer */ \
end_timer(18); \
\
/* Do inverse FFT levels 1,2,3 */ \
/* On input the 64-byte cache lines hold these data values: */ \
/*	0K	64K	128K	192K */ \
/*	1K	... */ \
/*	... */ \
/*	63K */ \
/*	256K	... */ \
/*	... */ \
/* On output the 64-byte cache lines hold these data values: */ \
/*	0K	256K	512K	768K */ \
/*	1K	... */ \
/*	... */ \
/*	255K */ \
\
/* Do 128 four_complex_unfft operations */ \
/*	distance between fft data elements is 128K */ \
/*	do 64 macros each processing 32 data values */ \
\
start_timer(20); \
rdx = (uintptr_t)g->sincos1;/* Load fixed negacyclic premultipliers */ \
rdi = (uintptr_t)g->sincos2;/* Load fixed sin/cos postmultipliers */ \
set_data_prefetch_ptrs;		/* Init prefetching pointers */ \
loops_init(64*2*clm);		/* 64 iterations of 2*clm */ \
for(unsigned int loopA = 8; loopA; loopA--) { \
	for(unsigned int loopB = 8*2*clm; loopB; ) { \
		xtouch(xptr(rcx+rbp+(clm-1)*128));	/* Preload TLB */ \
		do{ \
			r4_x4cl_four_complex_last_djbunfft(rsi, 64, 8*clmblkdst8, 16*clmblkdst8, rdi, 64*XMM_SCD2, rdx, 64*XMM_SCD4); \
		}while(--loopB & (2*clm-1)); \
		newxloop1(clm, c1b, clmblkdst, rdi, XMM_SCD2, rdx, XMM_SCD4); /* 2*clm iterations */ \
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
copy_scratch_data_to_fft(256, clm); \
goto	endlp;			/* Skip FFT code if flag not set */ \
\
/* Fire up auxiliary threads to do pass 1 forward FFT */ \
\
/* pass1: */pass1_forward_fft_setup(b0b); \
\
/* Copy data from the FFT data area back to the scratch area */ \
\
b0b:	copy_fft_data_to_scratch(256, clm); \
\
/* Do FFT levels 1,2,3 */ \
/* */ \
/* On input the 64-byte cache lines hold these data values: */ \
/*	0K	256K	512K	768K */ \
/*	1K	... */ \
/*	... */ \
/*	255K */ \
/* On output the 64-byte cache lines hold these data values: */ \
/*	0K	128K	256K	384K */ \
/*	1K	... */ \
/*	... */ \
/*	127K	... */ \
/*	512K	... */ \
/*	... */ \
\
/* Do 128 four_complex_first_fft operations */ \
/*	distance between fft data elements is 128K */ \
/*	do 128 iterations each processing 16 data values */ \
\
yes_fft:start_timer(22); \
rsi = (uintptr_t)g->scratch_area;	/* Get address of scratch area */ \
rdx = (uintptr_t)g->sincos1;/* Address of fixed negacyclic premultipliers */ \
rdi = (uintptr_t)g->sincos2;/* Load fixed sin/cos postmultipliers */ \
set_data_prefetch_ptrs;		/* Init prefetching pointers */ \
pfing(rcx += 64*blkdst);		/* We've already prefetched 64 blocks */ \
loops_init(128*2*clm);		/* 16 iters of 4 iters of 2*2*clm */ \
for(unsigned int loopA = 16; loopA; loopA--) { \
	for(unsigned int loopB = 4*2*2*clm; loopB; ) { \
		xtouch(xptr(rcx+rbp+(clm-1)*128));	/* Preload TLB */ \
		do{ \
			do{ \
				r4_x2cl_four_complex_first_djbfft_scratch(rsi, 64, 16*clmblkdst8, rdi, rdx); \
			}while(--loopB & (2*clm-1)); \
			newxloop2(clm, b1b, clmblkdst, rdi, XMM_SCD2, rdx, XMM_SCD4); /* 2*2*clm iters prefetching 2*clm cache lines */ \
		}while(loopB & (2*2*clm-1)); \
	} \
	rsi += -4*2*clmblkdst+clmblkdst8; \
} \
rsi += -16*clmblkdst8;	/* Restore source pointer */ \
end_timer(22); \
\
/* Do FFT levels 4,5 */ \
/* */ \
/* On input the 64-byte cache lines hold these data values: */ \
/*	0K	128K	256K	384K */ \
/*	1K	... */ \
/*	... */ \
/*	127K	... */ \
/*	512K	... */ \
/*	... */ \
/* On output the 64-byte cache lines hold these data values: */ \
/*	0K	32K	64K	96K */ \
/*	1K	... */ \
/*	... */ \
/*	31K	... */ \
/*	128K	... */ \
/*	... */ \
\
/* Do 128 four_complex_fft operations */ \
/*	distance between fft data elements is 32K */ \
/*	do 64 macros each processing 32 data values */ \
\
start_timer(24); \
rdi = (uintptr_t)g->premult_addr;/* Address of premultiplier/sin/cos table */ \
rdi += 8*2*clm*XMM_SCND2 + 4*2*clm*XMM_SCD8; \
loops_init(64*2*clm);		/* 2 iters of 4 iters of 8 of 2*clm */ \
r4_x4cl_four_complex_djbfft_preload; \
for(unsigned int loopA = 2; loopA; loopA--) { \
	for(unsigned int loopB = 4; loopB; loopB--) { \
		for(unsigned int loopC = 8*2*clm; loopC; ) { \
			xtouch(xptr(rcx+rbp+(clm-1)*128));	/* Preload TLB */ \
			do{ \
				r4_x4cl_four_complex_djbfft(rsi, 64, 4*clmblkdst8, 8*clmblkdst8, rdi); \
				rdi += XMM_SCD2;		/* Next sine/cosine pointer */ \
			}while(--loopC & (2*clm-1)); \
			newxloop1(clm, b5b, clmblkdst);	/* 2*clm iterations */ \
		} \
		rsi += -8*clmblkdst+clmblkdst8; \
	} \
	rsi += -4*clmblkdst8+16*clmblkdst8; /* Next source pointer */ \
	rdi += -32*2*clm*XMM_SCD2;	/* Restore sin/cos ptr */ \
} \
rsi += -2*16*clmblkdst8;	/* Restore source pointer */ \
end_timer(24); \
\
/* Do FFT levels 6,7 */ \
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
/* Do 128 four_complex_fft operations */ \
/*	distance between fft data elements is 8K */ \
/*	do 64 macros each processing 32 data values */ \
\
start_timer(25); \
rdi += -8*2*clm*XMM_SCND2;	/* Next sin/cos ptr */ \
loops_init(64*2*clm);		/* 8 iterations of 8 iters of 2*clm */ \
r4_x4cl_wpn_four_complex_djbfft_preload; \
for(unsigned int loopA = 8; loopA; loopA--) { \
	for(unsigned int loopB = 8*2*clm; loopB; ) { \
		xtouch(xptr(rcx+rbp+(clm-1)*128));	/* Preload TLB */ \
		do{ \
			r4_x4cl_wpn_four_complex_djbfft(rsi, 64, clmblkdst8, 2*clmblkdst8, rdi); \
			rdi += XMM_SCND2;		/* Next sine/cosine pointer */ \
		}while(--loopB & (2*clm-1)); \
		newxloop1(clm, b7b, clmblkdst);	/* 2*clm iterations */ \
	} \
	rsi += -8*clmblkdst+4*clmblkdst8; /* Next source pointer */ \
	rdi += -8*2*clm*XMM_SCND2;	/* Restore sin/cos ptr */ \
} \
rsi += -8*4*clmblkdst8;	/* Restore source pointer */ \
end_timer(25); \
\
/* Do FFT levels 8,9,10 */ \
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
rdi += -4*2*clm*XMM_SCD8+128; /* Next premult/sin/cos ptr (gen 1 byte offsets in radix-8 macro) */ \
rcx = (uintptr_t)g->data_addr;/* Destination pointer */ \
loops_init_prefetch(32*2*clm, 64, 1, rbp); /* 32 iters of 2*clm */ \
pfing(rbp = (uintptr_t)g->premult_prefetch);/* Load premult/sin/cos prefetch pointer */ \
for(unsigned int loopA = 4; loopA; loopA--) { \
	for(unsigned int loopB = 8; loopB; loopB--) { \
		for(unsigned int loopC = 2*clm; loopC; loopC--) { \
			r8_sg8cl_eight_complex_fft8(rsi, 64, clmblkdst, 2*clmblkdst, 4*clmblkdst, rcx, 64, blkdst, 2*blkdst, 4*blkdst, rdi-128); \
			rdi += XMM_SCD8;		/* Next premult/sin/cos ptr */ \
		} \
		rsi += -2*clm*64+clmblkdst8; /* Next source pointer */ \
		rcx += -2*clm*64+8*blkdst;	/* Next dest pointer */ \
		rdi += -2*clm*XMM_SCD8;	/* Restore premult/sine/cosine pointer */ \
	} \
	rdi += 2*clm*XMM_SCD8;	/* Next premult/sine/cosine pointer */ \
} \
end_timer(26); \
\
/* Work on next 4*clm sets of 1024 values */ \
\
endlp:	xpass1_get_next_block(pass2, c0b, b0b, _wpn); /* Get next block and jump */ \
/* done: */ xfft_footer_rn; \
\
\
\
/* Perform a pass 1 negacyclic FFT of length 1280 doubles. */ \
/* Memory examples are for a 1280K FFT with 10 levels done in pass 2. */ \
/* Prefetch 32*2*clm*XMM_SCD2 + 8*2*clm*XMM_SCND2 + 4*2*clm*XMM_SCD8 = 64*2*clm 64-byte cache lines */ \
\

#define r4dwpn_pass1sc1280ac(procname, xpass2name, pass2_levels, clm) /* (fftname, pass2_levels, clmarg) */ \
\
/* set_FFT_constants(fftname, 0, 1, pass2_levels, clmarg); */ \
xfft_header_rn(procname, xpass2name, pass2_levels, clm); \
\
/* Do pass 2 (the last 10 FFT levels if this is a 1280K FFT) */ \
\
/* pass2: */start_timer(0); \
xpass2name(g); \
end_timer(0); \
if(g->ffttype == 1)		/* We're done if FFTing only */ \
	return; /* goto done; */ \
\
/* Do the inverse FFT */ \
\
pass1_inverse_fft_setup(c0b, procname); \
\
/* Do inverse FFT levels 8,9,10 */ \
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
/* Do 80 eight_complex_unfft8 operations */ \
/*	distance between fft data elements is 1K */ \
/*	do 80 macros each processing 32 data values */ \
\
c0b:	start_timer(16); \
rsi = (uintptr_t)g->data_addr;/* Get source address */ \
rcx = (uintptr_t)g->scratch_area;	/* Get address of scratch area */ \
rdi = (uintptr_t)g->premult_addr;/* Load premult/sin/cos pointer */ \
pfing(rbp = (rdi+(64-40)*2*clm*64)); /* Finish prefetching the sin/cos/premult data */ \
loops_init_prefetch(80*2*clm, 64, 2, rbp); /* 80*2*clm iterations */ \
/*bug	xtouchpages rbp, 80*2*clm * 64 / 2 ;; Prefetch TLBs */ \
rdi += 128;		/* Generates 1 byte offsets in radix-8 macro */ \
for(unsigned int loopA = 4; loopA; loopA--) { \
	for(unsigned int loopB = 10; loopB; loopB--) { \
		for(unsigned int loopC = 2; loopC; loopC--) { \
			for(unsigned int loopD = 2*clm; loopD; loopD--) { \
				r8_sg4cl_eight_complex_unfft8(rsi, 64, blkdst, 2*blkdst, rcx, 64, clmblkdst, 2*clmblkdst, rdi-128); \
				rdi += XMM_SCD8;		/* Next premult/sin/cos data */ \
			} \
			rsi += -2*clm*64+4*blkdst;	/* Next source pointer */ \
			rcx += -2*clm*64+4*clmblkdst; /* Next dest pointer */ \
			rdi += -2*clm*XMM_SCD8;	/* Restore premult/sin/cos data pointer */ \
		} \
		rcx += -2*4*clmblkdst+clmblkdst8; \
	} \
	rdi += 2*clm*XMM_SCD8;	/* Next premult/sin/cos pointer */ \
} \
rsi = (uintptr_t)g->scratch_area;	/* Restore source pointer */ \
rdi += -128;		/* Next sin/cos data pointer */ \
end_timer(16); \
\
/* Do inverse FFT levels 6,7 */ \
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
/* Do 160 four_complex_unfft operations */ \
/*	distance between fft data elements is 8K */ \
/*	do 80 macros each processing 32 data values */ \
\
start_timer(17); \
loops_init(80*2*clm);		/* 80 iterations of 2*clm */ \
r4_x4cl_wpn_four_complex_djbunfft_preload; \
for(unsigned int loopA = 20; loopA; loopA--) { \
	for(unsigned int loopB = 4; loopB; loopB--) { \
		for(unsigned int loopC = 2*clm; loopC; loopC--) { \
			r4_x4cl_wpn_four_complex_djbunfft(rsi, 64, 4*clmblkdst, clmblkdst8, rdi, 4*2*clm*XMM_SCND2); \
			rdi += XMM_SCND2;		/* Next clm's sin/cos data */ \
		} \
		rsi += -2*clm*64+clmblkdst; /* Next source pointer */ \
	} \
	rsi += -4*clmblkdst+2*clmblkdst8; /* Next source pointer */ \
	rdi += -4*2*clm*XMM_SCND2;	/* Restore sine/cosine pointer */ \
} \
rsi += -20*2*clmblkdst8;	/* Restore source pointer */ \
rdi += 8*2*clm*XMM_SCND2;	/* Next sin/cos data pointer */ \
end_timer(17); \
\
/* Do inverse FFT levels 4,5 (actually 2.32 levels) */ \
/* On input the 64-byte cache lines hold these data values: */ \
/*	0K	16K	32K	48K */ \
/*	1K	... */ \
/*	... */ \
/*	15K */ \
/*	64K	... */ \
/*	... */ \
/* On output the 64-byte cache lines hold these data values: */ \
/*	0K	80K	160K	240K */ \
/*	1K	... */ \
/*	... */ \
/*	79K */ \
/*	320K	... */ \
/*	... */ \
\
/* Do 128 five_complex_unfft operations */ \
/*	distance between fft data elements is 32K */ \
/*	do 64 macros each processing 40 data values */ \
\
start_timer(18); \
loops_init(64*2*clm);		/* 64 iterations of 2*clm */ \
r5_x5cl_five_complex_djbunfft_preload; \
for(unsigned int loopA = 4; loopA; loopA--) { \
	for(unsigned int loopB = 2; loopB; loopB--) { \
		for(unsigned int loopC = 8; loopC; loopC--) { \
			for(unsigned int loopD = 2*clm; loopD; loopD--) { \
				r5_x5cl_five_complex_djbunfft(rsi, 64, 2*clmblkdst8, rdi, 16*2*clm*XMM_SCD2); \
				rdi += XMM_SCD2;		/* Next clm's sin/cos data */ \
			} \
			rsi += -2*clm*64+clmblkdst; /* Next source pointer */ \
		} \
		rsi += -8*clmblkdst+clmblkdst8; /* Next source pointer */ \
	} \
	rsi += -2*clmblkdst8+10*clmblkdst8; /* Next source pointer */ \
	rdi += -2*8*2*clm*XMM_SCD2; /* Restore sine/cosine pointer */ \
} \
rsi += -4*10*clmblkdst8;	/* Restore source pointer */ \
end_timer(18); \
\
/* Do inverse FFT levels 1,2,3 */ \
/* On input the 64-byte cache lines hold these data values: */ \
/*	0K	80K	160K	240K */ \
/*	1K	... */ \
/*	... */ \
/*	79K */ \
/*	320K	... */ \
/*	... */ \
/* On output the 64-byte cache lines hold these data values: */ \
/*	0K	320K	640K	960K */ \
/*	1K	... */ \
/*	... */ \
/*	319K */ \
\
/* Do 160 four_complex_unfft operations */ \
/*	distance between fft data elements is 160K */ \
/*	do 80 macros each processing 32 data values */ \
\
start_timer(20); \
rdx = (uintptr_t)g->sincos1;/* Load fixed negacyclic premultipliers */ \
rdi = (uintptr_t)g->sincos2;/* Load fixed sin/cos postmultipliers */ \
set_data_prefetch_ptrs;		/* Init prefetching pointers */ \
loops_init(80*2*clm);		/* 80 iterations of 2*clm */ \
for(unsigned int loopA = 10; loopA; loopA--) { \
	for(unsigned int loopB = 4*2*2*clm; loopB; ) { \
		xtouch(xptr(rcx+rbp+(clm-1)*128));	/* Preload TLB */ \
		do{ \
			do{ \
				r4_x4cl_four_complex_last_djbunfft(rsi, 64, 10*clmblkdst8, 20*clmblkdst8, rdi, 80*XMM_SCD2, rdx, 80*XMM_SCD4); \
			}while(--loopB & (2*clm-1)); \
			newxloop2(clm, c1b, clmblkdst, rdi, XMM_SCD2, rdx, XMM_SCD4); /* 2*2*clm iterations */ \
		}while(loopB & (2*2*clm-1)); \
	} \
	rsi += -4*2*clmblkdst+clmblkdst8; /* Next source pointer */ \
} \
end_timer(20); \
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
/* pass1: */pass1_forward_fft_setup(b0b); \
\
/* Copy data from the FFT data area back to the scratch area */ \
\
b0b:	copy_fft_data_to_scratch(320, clm); \
\
/* Do FFT levels 1,2,3 */ \
/* */ \
/* On input the 64-byte cache lines hold these data values: */ \
/*	0K	320K	640K	960K */ \
/*	1K	... */ \
/*	... */ \
/*	319K */ \
/* On output the 64-byte cache lines hold these data values: */ \
/*	0K	160K	320K	480K */ \
/*	1K	... */ \
/*	... */ \
/*	159K	... */ \
/*	640K	... */ \
/*	... */ \
\
/* Do 160 four_complex_first_fft operations */ \
/*	distance between fft data elements is 160K */ \
/*	do 160 iterations each processing 16 data values */ \
\
yes_fft:start_timer(22); \
rsi = (uintptr_t)g->scratch_area;	/* Get address of scratch area */ \
rdx = (uintptr_t)g->sincos1;/* Address of fixed negacyclic premultipliers */ \
rdi = (uintptr_t)g->sincos2;/* Load fixed sin/cos postmultipliers */ \
set_data_prefetch_ptrs;		/* Init prefetching pointers */ \
pfing(rcx += (40-8)*blkdst);	/* We've already prefetched 40 blocks, but */ \
;					/* code below prefetches 8 blocks too many */ \
loops_init(160*2*clm);		/* 20 iters of 4 iters of 2*2*clm */ \
for(unsigned int loopA = 20; loopA; loopA--) { \
	for(unsigned int loopB = 4*2*2*clm; loopB; ) { \
		xtouch(xptr(rcx+rbp+(clm-1)*128));	/* Preload TLB */ \
		do{ \
			do{ \
				r4_x2cl_four_complex_first_djbfft_scratch(rsi, 64, 20*clmblkdst8, rdi, rdx); \
			}while(--loopB & (2*clm-1)); \
			newxloop2(clm, b1b, clmblkdst, rdi, XMM_SCD2, rdx, XMM_SCD4); /* 2*2*clm iters prefetching 2*clm cache lines */ \
		}while(loopB & (2*2*clm-1)); \
	} \
	rsi += -4*2*clmblkdst+clmblkdst8; \
} \
rsi += -20*clmblkdst8;	/* Restore source pointer */ \
end_timer(22); \
\
/* Do FFT levels 4,5 (actually 2.32 levels) */ \
/* */ \
/* On input the 64-byte cache lines hold these data values: */ \
/*	0K	160K	320K	480K */ \
/*	1K	... */ \
/*	... */ \
/*	159K	... */ \
/*	640K	... */ \
/*	... */ \
/* On output the 64-byte cache lines hold these data values: */ \
/*	0K	32K	64K	96K */ \
/*	1K	... */ \
/*	... */ \
/*	31K	... */ \
/*	128K	... */ \
/*	... */ \
\
/* Do 128 five_complex_fft operations */ \
/*	distance between fft data elements is 40K */ \
/*	do 64 macros each processing 40 data values */ \
\
start_timer(24); \
rdi = (uintptr_t)g->premult_addr;/* Address of premultiplier/sin/cos table */ \
rdi += 8*2*clm*XMM_SCND2 + 4*2*clm*XMM_SCD8; \
loops_init(64*2*clm);		/* 2 iters of 4 iters of 8 of 2*clm */ \
r5_x5cl_five_complex_djbfft_preload; \
for(unsigned int loopA = 2; loopA; loopA--) { \
	for(unsigned int loopB = 4; loopB; loopB--) { \
		for(unsigned int loopC = 8*2*clm; loopC; ) { \
			xtouch(xptr(rcx+rbp+(clm-1)*128));	/* Preload TLB */ \
			do{ \
				r5_x5cl_five_complex_djbfft(rsi, 64, 4*clmblkdst8, rdi); \
				rdi += XMM_SCD2;		/* Next sine/cosine pointer */ \
			}while(--loopC & (2*clm-1)); \
			newxloop1doublepre(clm, b5b, b5a, clmblkdst); /* 2*clm iterations */ \
		} \
		rsi += -8*clmblkdst+clmblkdst8; \
	} \
	rsi += -4*clmblkdst8+20*clmblkdst8; /* Next source pointer */ \
	rdi += -32*2*clm*XMM_SCD2;	/* Restore sin/cos ptr */ \
} \
rsi += -2*20*clmblkdst8;	/* Restore source pointer */ \
end_timer(24); \
\
/* Do FFT levels 6,7 */ \
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
/* Do 160 four_complex_fft operations */ \
/*	distance between fft data elements is 8K */ \
/*	do 80 macros each processing 32 data values */ \
\
start_timer(25); \
rdi += -8*2*clm*XMM_SCND2;	/* Next sin/cos ptr */ \
loops_init(80*2*clm);		/* 20 iterations of 2 iters of 2*2*clm */ \
r4_x4cl_wpn_four_complex_djbfft_preload; \
for(unsigned int loopA = 10; loopA; loopA--) { \
	for(unsigned int loopB = 8*2*clm; loopB; ) { \
		xtouch(xptr(rcx+rbp+(clm-1)*128));	/* Preload TLB */ \
		do{ \
			r4_x4cl_wpn_four_complex_djbfft(rsi, 64, clmblkdst8, 2*clmblkdst8, rdi); \
			rdi += XMM_SCND2;		/* Next sine/cosine pointer */ \
		}while(--loopB & (2*clm-1)); \
		newxloop1(clm, b7b, clmblkdst);	/* 2*clm iterations */ \
	} \
	rsi += -8*clmblkdst+4*clmblkdst8; /* Next source pointer */ \
	rdi += -8*2*clm*XMM_SCND2;	/* Restore sin/cos ptr */ \
} \
rsi += -10*4*clmblkdst8;	/* Restore source pointer */ \
end_timer(25); \
\
/* Do FFT levels 8,9,10 */ \
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
/* Do 80 eight_complex_fft operations */ \
/*	distance between fft data elements is 1K */ \
/*	do 40 macros each processing 64 data values */ \
\
start_timer(26); \
rdi += -4*2*clm*XMM_SCD8+128; /* Next premult/sin/cos ptr (gen 1 byte offsets in radix-8 macro) */ \
rcx = (uintptr_t)g->data_addr;/* Destination pointer */ \
loops_init_prefetch(40*2*clm, 64, 1, rbp); /* 40 iters of 2*clm */ \
pfing(rbp = (uintptr_t)g->premult_prefetch);/* Load premult/sin/cos prefetch pointer */ \
for(unsigned int loopA = 4; loopA; loopA--) { \
	for(unsigned int loopB = 10; loopB; loopB--) { \
		for(unsigned int loopC = 2*clm; loopC; loopC--) { \
			r8_sg8cl_eight_complex_fft8(rsi, 64, clmblkdst, 2*clmblkdst, 4*clmblkdst, rcx, 64, blkdst, 2*blkdst, 4*blkdst, rdi-128); \
			rdi += XMM_SCD8;		/* Next premult/sin/cos ptr */ \
		} \
		rsi += -2*clm*64+clmblkdst8; /* Next source pointer */ \
		rcx += -2*clm*64+8*blkdst;	/* Next dest pointer */ \
		rdi += -2*clm*XMM_SCD8;	/* Restore premult/sine/cosine pointer */ \
	} \
	rdi += 2*clm*XMM_SCD8;	/* Next premult/sine/cosine pointer */ \
} \
end_timer(26); \
\
/* Work on next 4*clm sets of 1280 values */ \
\
endlp:	xpass1_get_next_block(pass2, c0b, b0b, _wpn); /* Get next block and jump */ \
/* done: */ xfft_footer_rn; \
\
\
\
/* Perform a pass 1 negacyclic FFT of length 1536 doubles. */ \
/* Memory examples are for a 1536K FFT with 10 levels done in pass 2. */ \
/* Sin/cos data needed: 16*2*clm*XMM_SCD1 + 4*2*clm*XMM_SCND2 + 16*2*clm*XMM_SCD4 = 48*2*clm 64-byte cache lines */ \
\

#define r4dwpn_pass1sc1536ac(procname, xpass2name, pass2_levels, clm) /* (fftname, pass2_levels, clmarg) */ \
\
/* set_FFT_constants(fftname, 0, 1, pass2_levels, clmarg); */ \
xfft_header_rn(procname, xpass2name, pass2_levels, clm); \
\
/* Do pass 2 (the last 10 FFT levels if this is a 1536K FFT) */ \
\
/* pass2: */start_timer(0); \
xpass2name(g); \
end_timer(0); \
if(g->ffttype == 1)		/* We're done if FFTing only */ \
	return; /* goto done; */ \
\
/* Do the inverse FFT */ \
\
pass1_inverse_fft_setup(c0b, procname); \
\
/* Do inverse FFT levels 10,11 */ \
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
/* Do 192 four_complex_unfft4 operations */ \
/*	distance between fft data elements is 1K */ \
/*	do 192 macros each processing 16 data values */ \
\
c0b:	start_timer(16); \
rsi = (uintptr_t)g->data_addr;/* Get source address */ \
rcx = (uintptr_t)g->scratch_area;	/* Get address of scratch area */ \
rdi = (uintptr_t)g->premult_addr;/* Load premult/sin/cos pointer */ \
pfing(rbp = (rdi+(48-24)*2*clm*64)); /* Finish prefetching the sin/cos/premult data */ \
loops_init_prefetch(192*2*clm, 64, 8, rbp); /* 192*2*clm iterations */ \
/*bug	xtouchpages rbp, 192*2*clm * 64 / 8 ;; Prefetch TLBs */ \
for(unsigned int loopA = 16; loopA; loopA--) { \
	for(unsigned int loopB = 3; loopB; loopB--) { \
		for(unsigned int loopC = 4; loopC; loopC--) { \
			for(unsigned int loopD = 2*clm; loopD; loopD--) { \
				r4_sg2cl_four_complex_unfft4(rsi, 64, blkdst, rcx, 64, clmblkdst, rdi); \
				rdi += XMM_SCD4;		/* Next premult/sin/cos data */ \
			} \
			rsi += -2*clm*64+2*blkdst;	/* Next source pointer */ \
			rcx += -2*clm*64+2*clmblkdst; /* Next dest pointer */ \
			rdi += -2*clm*XMM_SCD4;	/* Restore premult/sin/cos data pointer */ \
		} \
		rcx += -4*2*clmblkdst+clmblkdst8; \
	} \
	rdi += 2*clm*XMM_SCD4;	/* Next premult/sin/cos pointer */ \
} \
rsi = (uintptr_t)g->scratch_area;	/* Restore source pointer */ \
end_timer(16); \
\
/* Do inverse FFT levels 8,9 */ \
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
/* Do 192 four_complex_unfft operations */ \
/*	distance between fft data elements is 4K */ \
/*	do 96 macros each processing 32 data values */ \
\
start_timer(17); \
loops_init(96*2*clm);		/* 96*2*clm iterations */ \
r4_x4cl_wpn_four_complex_djbunfft_preload; \
for(unsigned int loopA = 48; loopA; loopA--) { \
	for(unsigned int loopB = 2; loopB; loopB--) { \
		for(unsigned int loopC = 2*clm; loopC; loopC--) { \
			r4_x4cl_wpn_four_complex_djbunfft(rsi, 64, 2*clmblkdst, 4*clmblkdst, rdi, 2*2*clm*XMM_SCND2); \
			rdi += XMM_SCND2;		/* Next clm's sin/cos data */ \
		} \
		rsi += -2*clm*64+clmblkdst; /* Next source pointer */ \
	} \
	rsi += -2*clmblkdst+clmblkdst8; /* Next source pointer */ \
	rdi += -2*2*clm*XMM_SCND2;	/* Restore sine/cosine pointer */ \
} \
rsi += -48*clmblkdst8;	/* Restore source pointer */ \
rdi += 4*2*clm*XMM_SCND2;	/* Next sin/cos data pointer */ \
end_timer(17); \
\
/* Do inverse FFT levels 6,7 (actually 1.585 levels) */ \
/* On input the 64-byte cache lines hold these data values: */ \
/*	0K	8K	16K	24K */ \
/*	1K	... */ \
/*	... */ \
/*	7K */ \
/*	32K	... */ \
/*	... */ \
/* On output the 64-byte cache lines hold these data values: */ \
/*	0K	24K	48K	72K */ \
/*	1K	... */ \
/*	... */ \
/*	23K */ \
/*	96K	... */ \
/*	... */ \
\
/* Do 256 three_complex_unfft operations */ \
/*	distance between fft data elements is 16K */ \
/*	do 128 macros each processing 24 data values */ \
\
start_timer(18); \
loops_init(128*2*clm);		/* 128 iterations of 2*clm */ \
r3_x3cl_three_complex_djbunfft_preload; \
for(unsigned int loopA = 16; loopA; loopA--) { \
	for(unsigned int loopB = 8; loopB; loopB--) { \
		for(unsigned int loopC = 2*clm; loopC; loopC--) { \
			r3_x3cl_three_complex_djbunfft(rsi, 64, clmblkdst8, rdi, 8*2*clm*XMM_SCD1); \
			rdi += XMM_SCD1;		/* Next clm's sin/cos data */ \
		} \
		rsi += -2*clm*64+clmblkdst; /* Next source pointer */ \
	} \
	rsi += -8*clmblkdst+3*clmblkdst8; /* Next source pointer */ \
	rdi += -8*2*clm*XMM_SCD1;	/* Restore sine/cosine pointer */ \
} \
rsi += -16*3*clmblkdst8;	/* Restore source pointer */ \
end_timer(18); \
\
/* Do inverse FFT levels 4,5 */ \
/* On input the 64-byte cache lines hold these data values: */ \
/*	0K	24K	48K	72K */ \
/*	1K	... */ \
/*	... */ \
/*	23K */ \
/*	96K	... */ \
/*	... */ \
/* On output the 64-byte cache lines hold these data values: */ \
/*	0K	96K	192K	288K */ \
/*	1K	... */ \
/*	... */ \
/*	23K */ \
/*	96K	... */ \
/*	... */ \
\
/* Do 192 four_complex_unfft operations */ \
/*	distance between fft data elements is 48K */ \
/*	do 96 macros each processing 32 data values */ \
\
start_timer(19); \
rdi = (uintptr_t)g->sincos2;/* Load fixed sin/cos data ptr */ \
rdi += 192*XMM_SCD2;	/* Sin/cos ptr */ \
loops_init(96*2*clm);		/* 96*2*clm iterations */ \
r4_x4cl_four_complex_djbunfft_preload; \
for(unsigned int loopA = 4; loopA; loopA--) { \
	for(unsigned int loopB = 3; loopB; loopB--) { \
		for(unsigned int loopC = 8; loopC; loopC--) { \
			for(unsigned int loopD = 2*clm; loopD; loopD--) { \
				r4_x4cl_four_complex_djbunfft(rsi, 64, 3*clmblkdst8, 6*clmblkdst8, rdi, 24*XMM_SCD2); \
			} \
			rsi += -2*clm*64+clmblkdst; /* Next source pointer */ \
			rdi += XMM_SCD2;		/* Next clm's sin/cos data */ \
		} \
		rsi += -8*clmblkdst+clmblkdst8; /* Next source pointer */ \
	} \
	rsi += -3*clmblkdst8+12*clmblkdst8; /* Next source pointer */ \
	rdi += -24*XMM_SCD2;	/* Restore sine/cosine pointer */ \
} \
rsi += -4*12*clmblkdst8;	/* Restore source pointer */ \
end_timer(19); \
\
/* Do inverse FFT levels 1,2,3 */ \
/* On input the 64-byte cache lines hold these data values: */ \
/*	0K	96K	192K	288K */ \
/*	1K	... */ \
/*	... */ \
/*	23K */ \
/*	96K	... */ \
/*	... */ \
/* On output the 64-byte cache lines hold these data values: */ \
/*	0K	384K	768K	1152K */ \
/*	1K	... */ \
/*	... */ \
/*	383K */ \
\
/* Do 192 four_complex_unfft operations */ \
/*	distance between fft data elements is 192K */ \
/*	do 96 macros each processing 32 data values */ \
\
start_timer(20); \
rdx = (uintptr_t)g->sincos1;/* Load fixed negacyclic premultipliers */ \
rdi = (uintptr_t)g->sincos2;/* Load fixed sin/cos postmultipliers */ \
loops_init(96*2*clm);		/* 96 iterations of 2*clm */ \
for(unsigned int loopA = 12; loopA; loopA--) { \
	for(unsigned int loopB = 8; loopB; loopB--) { \
		for(unsigned int loopC = 2*clm; loopC; loopC--) { \
			r4_x4cl_four_complex_last_djbunfft(rsi, 64, 12*clmblkdst8, 24*clmblkdst8, rdi, 96*XMM_SCD2, rdx, 96*XMM_SCD4); \
		} \
		newxloop1nopre(clm, c1b, clmblkdst, rdi, XMM_SCD2, rdx, XMM_SCD4); /* 2*clm iterations */ \
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
copy_scratch_data_to_fft(384, clm); \
goto	endlp;			/* Skip FFT code if flag not set */ \
\
/* Fire up auxiliary threads to do pass 1 forward FFT */ \
\
/* pass1: */pass1_forward_fft_setup(b0b); \
\
/* Copy data from the FFT data area back to the scratch area */ \
\
b0b:	copy_fft_data_to_scratch(384, clm); \
\
/* Do FFT levels 1,2,3 */ \
/* */ \
/* On input the 64-byte cache lines hold these data values: */ \
/*	0K	384K	768K	1152K */ \
/*	1K	... */ \
/*	... */ \
/*	383K */ \
/* On output the 64-byte cache lines hold these data values: */ \
/*	0K	192K	384K	576K */ \
/*	1K	... */ \
/*	... */ \
/*	191K	... */ \
/*	768K	... */ \
/*	... */ \
\
/* Do 192 four_complex_first_fft operations */ \
/*	distance between fft data elements is 192K */ \
/*	do 192 macros each processing 16 data values */ \
\
yes_fft:start_timer(22); \
rsi = (uintptr_t)g->scratch_area;	/* Get address of scratch area */ \
rdx = (uintptr_t)g->sincos1;/* Address of fixed negacyclic premultipliers */ \
rdi = (uintptr_t)g->sincos2;/* Load fixed sin/cos postmultipliers */ \
set_data_prefetch_ptrs;		/* Init prefetching pointers */ \
loops_init(192*2*clm);		/* 24 iters of 4 iters of 2*2*clm */ \
for(unsigned int loopA = 24; loopA; loopA--) { \
	for(unsigned int loopB = 4*2*2*clm; loopB; ) { \
		xtouch(xptr(rcx+rbp+(clm-1)*128));	/* Preload TLB */ \
		do{ \
			do{ \
				r4_x2cl_four_complex_first_djbfft_scratch(rsi, 64, 24*clmblkdst8, rdi, rdx); \
			}while(--loopB & (2*clm-1)); \
			newxloop2(clm, b1b, clmblkdst, rdi, XMM_SCD2, rdx, XMM_SCD4); /* 2*2*clm iters prefetching 2*clm cache lines */ \
		}while(loopB & (2*2*clm-1)); \
	} \
	rsi += -4*2*clmblkdst+clmblkdst8; \
} \
rsi += -24*clmblkdst8;	/* Restore source pointer */ \
pfing(rcx += -32*blkdst);		/* We just prefetched 96 blocks, correct for */ \
;					/* code below which prefetches 32 blocks too many */ \
end_timer(22); \
\
/* Do FFT levels 4,5 */ \
/* */ \
/* On input the 64-byte cache lines hold these data values: */ \
/*	0K	192K	384K	576K */ \
/*	1K	... */ \
/*	... */ \
/*	191K	... */ \
/*	768K	... */ \
/*	... */ \
/* On output the 64-byte cache lines hold these data values: */ \
/*	0K	48K	96K	144K */ \
/*	1K	... */ \
/*	... */ \
/*	47K	... */ \
/*	192K	... */ \
/*	... */ \
\
/* Do 192 four_complex_fft operations */ \
/*	distance between fft data elements is 48K */ \
/*	do 96 macros each processing 32 data values */ \
\
start_timer(23); \
loops_init(96*2*clm);		/* 6 iterations of 16 iters of 2*clm */ \
r4_x4cl_four_complex_djbfft_preload; \
for(unsigned int loopA = 2; loopA; loopA--) { \
	for(unsigned int loopB = 6; loopB; loopB--) { \
		for(unsigned int loopC = 8*2*clm; loopC; ) { \
			xtouch(xptr(rcx+rbp+(clm-1)*128));	/* Preload TLB */ \
			do{ \
				r4_x4cl_four_complex_djbfft(rsi, 64, 6*clmblkdst8, 12*clmblkdst8, rdi); \
			}while(--loopC & (2*clm-1)); \
			newxloop1(clm, b3b, clmblkdst, rdi, XMM_SCD2);	/* 2*clm iterations */ \
		} \
		rsi += -8*clmblkdst+clmblkdst8; /* Next source pointer */ \
	} \
	rsi += -6*clmblkdst8+24*clmblkdst8; /* Next source pointer */ \
	rdi += -48*XMM_SCD2;	/* Restore sin/cos ptr */ \
} \
rsi += -2*24*clmblkdst8;	/* Restore source pointer */ \
end_timer(23); \
\
/* Do FFT levels 6,7 (actually 1.585 levels) */ \
/* */ \
/* On input the 64-byte cache lines hold these data values: */ \
/*	0K	48K	96K	144K */ \
/*	1K	... */ \
/*	... */ \
/*	47K	... */ \
/*	192K	... */ \
/*	... */ \
/* On output the 64-byte cache lines hold these data values: */ \
/*	0K	16K	32K	48K */ \
/*	1K	... */ \
/*	... */ \
/*	15K	... */ \
/*	64K	... */ \
/*	... */ \
\
/* Do 256 three_complex_fft operations */ \
/*	distance between fft data elements is 16K */ \
/*	do 128 macros each processing 24 data values */ \
\
start_timer(24); \
rdi = (uintptr_t)g->premult_addr;/* Address of premultiplier/sin/cos table */ \
rdi += 4*2*clm*XMM_SCND2 + 16*2*clm*XMM_SCD4; \
loops_init(128*2*clm);		/* 2 iters of 8 iters of 8 of 2*clm */ \
r3_x3cl_three_complex_djbfft_preload; \
for(unsigned int loopA = 8; loopA; loopA--) { \
	for(unsigned int loopB = 2; loopB; loopB--) { \
		for(unsigned int loopC = 8*2*clm; loopC; ) { \
			xtouch(xptr(rcx+rbp+(clm-1)*128));	/* Preload TLB */ \
			do{ \
				r3_x3cl_three_complex_djbfft(rsi, 64, 2*clmblkdst8, rdi); \
				rdi += XMM_SCD1;		/* Next sine/cosine pointer */ \
			}while(--loopC & (2*clm-1)); \
			newxloop1(clm, b5b, clmblkdst);	/* 2*clm iterations */ \
		} \
		rsi += -8*clmblkdst+clmblkdst8; \
	} \
	rsi += -2*clmblkdst8+6*clmblkdst8; /* Next source pointer */ \
	rdi += -16*2*clm*XMM_SCD1;	/* Restore sin/cos ptr */ \
} \
rsi += -8*6*clmblkdst8;	/* Restore source pointer */ \
end_timer(24); \
\
/* Do FFT levels 8,9 */ \
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
/* Do 192 four_complex_fft operations */ \
/*	distance between fft data elements is 4K */ \
/*	do 96 macros each processing 32 data values */ \
\
start_timer(25); \
rdi += -4*2*clm*XMM_SCND2;	/* Next sin/cos ptr */ \
loops_init(96*2*clm);		/* 24 iterations of 4 iters of 2*clm */ \
r4_x4cl_wpn_four_complex_djbfft_preload; \
for(unsigned int loopA = 24; loopA; loopA--) { \
	for(unsigned int loopB = 4*2*clm; loopB; ) { \
		xtouch(xptr(rcx+rbp+(clm-1)*128));	/* Preload TLB */ \
		do{ \
			r4_x4cl_wpn_four_complex_djbfft(rsi, 64, 4*clmblkdst, clmblkdst8, rdi); \
			rdi += XMM_SCND2;		/* Next sine/cosine pointer */ \
		}while(--loopB & (2*clm-1)); \
		newxloop1(clm, b7b, clmblkdst);	/* 2*clm iterations */ \
	} \
	rsi += -4*clmblkdst+2*clmblkdst8; /* Next source pointer */ \
	rdi += -4*2*clm*XMM_SCND2;	/* Restore sin/cos ptr */ \
} \
rsi += -24*2*clmblkdst8;	/* Restore source pointer */ \
end_timer(25); \
\
/* Do FFT levels 10,11 */ \
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
/* Do 192 four_complex_fft operations */ \
/*	distance between fft data elements is 1K */ \
/*	do 96 macros each processing 32 data values */ \
\
start_timer(26); \
rdi += -16*2*clm*XMM_SCD4;	/* Next premult/sin/cos ptr */ \
rcx = (uintptr_t)g->data_addr;/* Destination pointer */ \
loops_init_prefetch(96*2*clm, 64, 4, rbp); /* 96 iters of 2*clm */ \
pfing(rbp = (uintptr_t)g->premult_prefetch);/* Load premult/sin/cos prefetch pointer */ \
for(unsigned int loopA = 16; loopA; loopA--) { \
	for(unsigned int loopB = 3; loopB; loopB--) { \
		for(unsigned int loopC = 2; loopC; loopC--) { \
			for(unsigned int loopD = 2*clm; loopD; loopD--) { \
				r4_sg4cl_four_complex_fft4(rsi, 64, clmblkdst, 2*clmblkdst, rcx, 64, blkdst, 2*blkdst, rdi); \
				rdi += XMM_SCD4;		/* Next premult/sin/cos ptr */ \
			} \
			rsi += -2*clm*64+4*clmblkdst; /* Next source pointer */ \
			rcx += -2*clm*64+4*blkdst;	/* Next dest pointer */ \
			rdi += -2*clm*XMM_SCD4;	/* Restore premult/sine/cosine pointer */ \
		} \
		rsi += -2*4*clmblkdst+clmblkdst8; /* Next source pointer */ \
	} \
	rdi += 2*clm*XMM_SCD4;	/* Next premult/sine/cosine pointer */ \
} \
end_timer(26); \
\
/* Work on next 4*clm sets of 1536 values */ \
\
endlp:	xpass1_get_next_block(pass2, c0b, b0b, _wpn); /* Get next block and jump */ \
/* done: */ xfft_footer_rn; \
\
\
\
/* Perform a pass 1 negacyclic FFT of length 2048 doubles. */ \
/* Memory examples are for a 2048K FFT with 10 levels done in pass 2. */ \
/* Sin/cos data needed: 16*2*clm*XMM_SCND2 + 4*2*clm*XMM_SCD2 + 16*2*clm*XMM_SCD4 = 68*2*clm 64-byte cache lines */ \
\

#define r4dwpn_pass1sc2048ac(procname, xpass2name, pass2_levels, clm) /* (fftname, pass2_levels, clmarg) */ \
\
/* set_FFT_constants(fftname, 0, 1, pass2_levels, clmarg); */ \
xfft_header_rn(procname, xpass2name, pass2_levels, clm); \
\
/* Do pass 2 (the last 10 FFT levels if this is a 2M FFT) */ \
\
/* pass2: */start_timer(0); \
xpass2name(g); \
end_timer(0); \
if(g->ffttype == 1)		/* We're done if FFTing only */ \
	return; /* goto done; */ \
\
/* Do the inverse FFT */ \
\
pass1_inverse_fft_setup(c0b, procname); \
\
/* Do inverse FFT levels 10,11 */ \
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
/* Do 256 four_complex_unfft4 operations */ \
/*	distance between fft data elements is 1K */ \
/*	do 256 macros each processing 16 data values */ \
\
c0b:	start_timer(16); \
rsi = (uintptr_t)g->data_addr;/* Get source address */ \
rcx = (uintptr_t)g->scratch_area;	/* Get address of scratch area */ \
rdi = (uintptr_t)g->premult_addr;/* Load premult/sin/cos pointer */ \
pfing(rbp = (rdi+(68-32)*2*clm*64)); /* Finish prefetching the sin/cos/premult data */ \
loops_init_prefetch(256*2*clm, 64, 8, rbp); /* 256*2*clm iterations */ \
/*bug	xtouchpages rbp, 256*2*clm * 64 / 8 ;; Prefetch TLBs */ \
for(unsigned int loopA = 16; loopA; loopA--) { \
	for(unsigned int loopB = 4; loopB; loopB--) { \
		for(unsigned int loopC = 4; loopC; loopC--) { \
			for(unsigned int loopD = 2*clm; loopD; loopD--) { \
				r4_sg2cl_four_complex_unfft4(rsi, 64, blkdst, rcx, 64, clmblkdst, rdi); \
				rdi += XMM_SCD4;		/* Next premult/sin/cos data */ \
			} \
			rsi += -2*clm*64+2*blkdst;	/* Next source pointer */ \
			rcx += -2*clm*64+2*clmblkdst; /* Next dest pointer */ \
			rdi += -2*clm*XMM_SCD4;	/* Restore premult/sin/cos data pointer */ \
		} \
		rcx += -4*2*clmblkdst+clmblkdst8; \
	} \
	rdi += 2*clm*XMM_SCD4;	/* Next premult/sin/cos pointer */ \
} \
rsi = (uintptr_t)g->scratch_area;	/* Restore source pointer */ \
end_timer(16); \
\
/* Do inverse FFT levels 8,9 */ \
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
/* Do 256 four_complex_unfft operations */ \
/*	distance between fft data elements is 4K */ \
/*	do 128 macros each processing 32 data values */ \
\
start_timer(17); \
loops_init(128*2*clm);		/* 128*2*clm iterations */ \
r4_x4cl_four_complex_djbunfft_preload; \
for(unsigned int loopA = 64; loopA; loopA--) { \
	for(unsigned int loopB = 2; loopB; loopB--) { \
		for(unsigned int loopC = 2*clm; loopC; loopC--) { \
			r4_x4cl_four_complex_djbunfft(rsi, 64, 2*clmblkdst, 4*clmblkdst, rdi, 2*2*clm*XMM_SCD2); \
			rdi += XMM_SCD2;		/* Next clm's sin/cos data */ \
		} \
		rsi += -2*clm*64+clmblkdst; /* Next source pointer */ \
	} \
	rsi += -2*clmblkdst+clmblkdst8; /* Next source pointer */ \
	rdi += -2*2*clm*XMM_SCD2;	/* Restore sine/cosine pointer */ \
} \
rsi += -64*clmblkdst8;	/* Restore source pointer */ \
rdi += 4*2*clm*XMM_SCD2;	/* Next sin/cos data pointer */ \
end_timer(17); \
\
/* Do inverse FFT levels 6,7 */ \
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
/*	128K	... */ \
/*	... */ \
\
/* Do 256 four_complex_unfft operations */ \
/*	distance between fft data elements is 16K */ \
/*	do 128 macros each processing 32 data values */ \
\
start_timer(18); \
loops_init(128*2*clm);		/* 128 iterations of 2*clm */ \
r4_x4cl_wpn_four_complex_djbunfft_preload; \
for(unsigned int loopA = 16; loopA; loopA--) { \
	for(unsigned int loopB = 8; loopB; loopB--) { \
		for(unsigned int loopC = 2*clm; loopC; loopC--) { \
			r4_x4cl_wpn_four_complex_djbunfft(rsi, 64, clmblkdst8, 2*clmblkdst8, rdi, 8*2*clm*XMM_SCND2); \
			rdi += XMM_SCND2;		/* Next clm's sin/cos data */ \
		} \
		rsi += -2*clm*64+clmblkdst; /* Next source pointer */ \
	} \
	rsi += -8*clmblkdst+4*clmblkdst8; /* Next source pointer */ \
	rdi += -8*2*clm*XMM_SCND2;	/* Restore sine/cosine pointer */ \
} \
rsi += -16*4*clmblkdst8;	/* Restore source pointer */ \
end_timer(18); \
\
/* Do inverse FFT levels 4,5 */ \
/* On input the 64-byte cache lines hold these data values: */ \
/*	0K	32K	64K	96K */ \
/*	1K	... */ \
/*	... */ \
/*	31K */ \
/*	128K	... */ \
/*	... */ \
/* On output the 64-byte cache lines hold these data values: */ \
/*	0K	128K	256K	384K */ \
/*	1K	... */ \
/*	... */ \
/*	127K */ \
/*	512K	... */ \
/*	... */ \
\
/* Do 256 four_complex_unfft operations */ \
/*	distance between fft data elements is 64K */ \
/*	do 128 macros each processing 32 data values */ \
\
start_timer(19); \
rdi = (uintptr_t)g->sincos2;/* Load fixed sin/cos data ptr */ \
rdi += 256*XMM_SCD2;	/* Sin/cos ptr */ \
loops_init(128*2*clm);		/* 128 iterations of 2*clm */ \
r4_x4cl_four_complex_djbunfft_preload; \
for(unsigned int loopA = 4; loopA; loopA--) { \
	for(unsigned int loopB = 4; loopB; loopB--) { \
		for(unsigned int loopC = 8; loopC; loopC--) { \
			for(unsigned int loopD = 2*clm; loopD; loopD--) { \
				r4_x4cl_four_complex_djbunfft(rsi, 64, 4*clmblkdst8, 8*clmblkdst8, rdi, 32*XMM_SCD2); \
			} \
			rsi += -2*clm*64+clmblkdst; /* Next source pointer */ \
			rdi += XMM_SCD2;		/* Next sin/cos data */ \
		} \
		rsi += -8*clmblkdst+clmblkdst8; /* Next source pointer */ \
	} \
	rsi += -4*clmblkdst8+16*clmblkdst8; /* Next source pointer */ \
	rdi += -32*XMM_SCD2;	/* Restore sine/cosine pointer */ \
} \
rsi += -4*16*clmblkdst8;	/* Restore source pointer */ \
end_timer(19); \
\
/* Do inverse FFT levels 1,2,3 */ \
/* On input the 64-byte cache lines hold these data values: */ \
/*	0K	128K	256K	384K */ \
/*	1K	... */ \
/*	... */ \
/*	127K */ \
/*	512K	... */ \
/*	... */ \
/* On output the 64-byte cache lines hold these data values: */ \
/*	0K	512K	1024K	1536K */ \
/*	1K	... */ \
/*	... */ \
/*	511K */ \
\
/* Do 256 four_complex_unfft operations */ \
/*	distance between fft data elements is 256K */ \
/*	do 128 macros each processing 32 data values */ \
\
start_timer(20); \
rdx = (uintptr_t)g->sincos1;/* Load fixed negacyclic premultipliers */ \
rdi = (uintptr_t)g->sincos2;/* Load fixed sin/cos postmultipliers */ \
loops_init(128*2*clm);		/* 128 iterations of 2*clm */ \
for(unsigned int loopA = 16; loopA; loopA--) { \
	for(unsigned int loopB = 8; loopB; loopB--) { \
		for(unsigned int loopC = 2*clm; loopC; loopC--) { \
			r4_x4cl_four_complex_last_djbunfft(rsi, 64, 16*clmblkdst8, 32*clmblkdst8, rdi, 128*XMM_SCD2, rdx, 128*XMM_SCD4); \
		} \
		newxloop1nopre(clm, c1b, clmblkdst, rdi, XMM_SCD2, rdx, XMM_SCD4); /* 2*clm iterations */ \
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
copy_scratch_data_to_fft(512, clm); \
goto	endlp;			/* Skip FFT code if flag not set */ \
\
/* Fire up auxiliary threads to do pass 1 forward FFT */ \
\
/* pass1: */pass1_forward_fft_setup(b0b); \
\
/* Copy data from the FFT data area back to the scratch area */ \
\
b0b:	copy_fft_data_to_scratch(512, clm); \
\
/* Do FFT levels 1,2,3 */ \
/* */ \
/* On input the 64-byte cache lines hold these data values: */ \
/*	0K	512K	1024K	1536K */ \
/*	1K	... */ \
/*	... */ \
/*	511K */ \
/* On output the 64-byte cache lines hold these data values: */ \
/*	0K	256K	512K	768K */ \
/*	1K	... */ \
/*	... */ \
/*	255K	... */ \
/*	1024K	... */ \
/*	... */ \
\
/* Do 256 four_complex_first_fft operations */ \
/*	distance between fft data elements is 256K */ \
/*	do 256 iterations each processing 16 data values */ \
\
yes_fft:start_timer(22); \
rsi = (uintptr_t)g->scratch_area;	/* Get address of scratch area */ \
rdx = (uintptr_t)g->sincos1;/* Address of fixed negacyclic premultipliers */ \
rdi = (uintptr_t)g->sincos2;/* Load fixed sin/cos postmultipliers */ \
set_data_prefetch_ptrs;		/* Init prefetching pointers */ \
loops_init(256*2*clm);		/* 32 iters of 4 iters of 2*2*clm */ \
for(unsigned int loopA = 32; loopA; loopA--) { \
	for(unsigned int loopB = 4*2*2*clm; loopB; ) { \
		xtouch(xptr(rcx+rbp+(clm-1)*128));	/* Preload TLB */ \
		do{ \
			do{ \
				r4_x2cl_four_complex_first_djbfft_scratch(rsi, 64, 32*clmblkdst8, rdi, rdx); \
			}while(--loopB & (2*clm-1)); \
			newxloop2(clm, b1b, clmblkdst, rdi, XMM_SCD2, rdx, XMM_SCD4); /* 2*2*clm iters prefetching 2*clm cache lines */ \
		}while(loopB & (2*2*clm-1)); \
	} \
	rsi += -4*2*clmblkdst+clmblkdst8; \
} \
rsi += -32*clmblkdst8;	/* Restore source pointer */ \
end_timer(22); \
\
/* Do FFT levels 4,5 */ \
/* */ \
/* On input the 64-byte cache lines hold these data values: */ \
/*	0K	256K	512K	768K */ \
/*	1K	... */ \
/*	... */ \
/*	255K	... */ \
/*	1024K	... */ \
/*	... */ \
/* On output the 64-byte cache lines hold these data values: */ \
/*	0K	64K	128K	192K */ \
/*	1K	... */ \
/*	... */ \
/*	63K	... */ \
/*	256K	... */ \
/*	... */ \
\
/* Do 256 four_complex_fft operations */ \
/*	distance between fft data elements is 64K */ \
/*	do 128 macros each processing 32 data values */ \
\
start_timer(23); \
loops_init(128*2*clm);		/* 2 iters of 8 iters of 8 of 2*clm */ \
r4_x4cl_four_complex_djbfft_preload; \
for(unsigned int loopA = 2; loopA; loopA--) { \
	for(unsigned int loopB = 8; loopB; loopB--) { \
		for(unsigned int loopC = 8*2*clm; loopC; ) { \
			xtouch(xptr(rcx+rbp+(clm-1)*128));	/* Preload TLB */ \
			do{ \
				r4_x4cl_four_complex_djbfft(rsi, 64, 8*clmblkdst8, 16*clmblkdst8, rdi); \
			}while(--loopC & (2*clm-1)); \
			newxloop1(clm, b3b, clmblkdst, rdi, XMM_SCD2);	/* 2*clm iterations */ \
		} \
		rsi += -8*clmblkdst+clmblkdst8; \
	} \
	rsi += -8*clmblkdst8+32*clmblkdst8; /* Next source pointer */ \
	rdi += -64*XMM_SCD2;	/* Restore sin/cos ptr */ \
} \
rsi += -2*32*clmblkdst8;	/* Restore source pointer */ \
end_timer(23); \
\
/* Do FFT levels 6,7 */ \
/* */ \
/* On input the 64-byte cache lines hold these data values: */ \
/*	0K	64K	128K	192K */ \
/*	1K	... */ \
/*	... */ \
/*	63K	... */ \
/*	256K	... */ \
/*	... */ \
/* On output the 64-byte cache lines hold these data values: */ \
/*	0K	16K	32K	48K */ \
/*	1K	... */ \
/*	... */ \
/*	15K	... */ \
/*	64K	... */ \
/*	... */ \
\
/* Do 256 four_complex_fft operations */ \
/*	distance between fft data elements is 16K */ \
/*	do 128 macros each processing 32 data values */ \
\
start_timer(24); \
rdi = (uintptr_t)g->premult_addr;/* Address of premultiplier/sin/cos table */ \
rdi += 4*2*clm*XMM_SCD2 + 16*2*clm*XMM_SCD4; \
loops_init(128*2*clm);		/* 8 iters of 2 iters of 8 of 2*clm */ \
r4_x4cl_wpn_four_complex_djbfft_preload; \
for(unsigned int loopA = 8; loopA; loopA--) { \
	for(unsigned int loopB = 2; loopB; loopB--) { \
		for(unsigned int loopC = 8*2*clm; loopC; ) { \
			xtouch(xptr(rcx+rbp+(clm-1)*128));	/* Preload TLB */ \
			do{ \
				r4_x4cl_wpn_four_complex_djbfft(rsi, 64, 2*clmblkdst8, 4*clmblkdst8, rdi); \
				rdi += XMM_SCND2;		/* Next sine/cosine pointer */ \
			}while(--loopC & (2*clm-1)); \
			newxloop1(clm, b5b, clmblkdst);	/* 2*clm iterations */ \
		} \
		rsi += -8*clmblkdst+clmblkdst8; \
	} \
	rsi += -2*clmblkdst8+8*clmblkdst8; /* Next source pointer */ \
	rdi += -16*2*clm*XMM_SCND2; /* Restore sin/cos ptr */ \
} \
rsi += -8*8*clmblkdst8;	/* Restore source pointer */ \
end_timer(24); \
\
/* Do FFT levels 8,9 */ \
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
/* Do 256 four_complex_fft operations */ \
/*	distance between fft data elements is 4K */ \
/*	do 128 macros each processing 32 data values */ \
\
start_timer(25); \
rdi += -4*2*clm*XMM_SCD2;	/* Next sin/cos ptr */ \
loops_init(128*2*clm);		/* 32 iterations of 4 iters of 2*clm */ \
r4_x4cl_four_complex_djbfft_preload; \
for(unsigned int loopA = 32; loopA; loopA--) { \
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
rsi += -32*2*clmblkdst8;	/* Restore source pointer */ \
end_timer(25); \
\
/* Do FFT levels 10,11 */ \
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
/* Do 64 four_complex_fft operations */ \
/*	distance between fft data elements is 1K */ \
/*	do 32 macros each processing 32 data values */ \
\
start_timer(26); \
rdi += -16*2*clm*XMM_SCD4;	/* Next premult/sin/cos ptr */ \
rcx = (uintptr_t)g->data_addr;/* Destination pointer */ \
loops_init_prefetch(128*2*clm, 64, 2, rbp); /* 128 iters of 2*clm */ \
pfing(rbp = (uintptr_t)g->premult_prefetch);/* Load premult/sin/cos prefetch pointer */ \
for(unsigned int loopA = 16; loopA; loopA--) { \
	for(unsigned int loopB = 4; loopB; loopB--) { \
		for(unsigned int loopC = 2; loopC; loopC--) { \
			for(unsigned int loopD = 2*clm; loopD; loopD--) { \
				r4_sg4cl_four_complex_fft4(rsi, 64, clmblkdst, 2*clmblkdst, rcx, 64, blkdst, 2*blkdst, rdi); \
				rdi += XMM_SCD4;		/* Next premult/sin/cos ptr */ \
			} \
			rsi += -2*clm*64+4*clmblkdst; /* Next source pointer */ \
			rcx += -2*clm*64+4*blkdst;	/* Next dest pointer */ \
			rdi += -2*clm*XMM_SCD4;	/* Restore premult/sine/cosine pointer */ \
		} \
		rsi += -2*4*clmblkdst+clmblkdst8; /* Next source pointer */ \
	} \
	rdi += 2*clm*XMM_SCD4;	/* Next premult/sine/cosine pointer */ \
} \
end_timer(26); \
\
/* Work on next 4*clm sets of 2048 values */ \
\
endlp:	xpass1_get_next_block(pass2, c0b, b0b, _wpn); /* Get next block and jump */ \
/* done: */ xfft_footer_rn; \
\
\
\
/* Perform a pass 1 negacyclic FFT of length 2560 doubles. */ \
/* Memory examples are for a 2560K FFT with 10 levels done in pass 2. */ \
/* Sin/cos data needed: 16*2*clm*XMM_SCD2 + 4*2*clm*XMM_SCND2 + 16*2*clm*XMM_SCD4 = 56*2*clm 64-byte cache lines */ \
\

#define r4dwpn_pass1sc2560ac(procname, xpass2name, pass2_levels, clm) /* (fftname, pass2_levels, clmarg) */ \
\
/* set_FFT_constants(fftname, 0, 1, pass2_levels, clmarg); */ \
xfft_header_rn(procname, xpass2name, pass2_levels, clm); \
\
/* Do pass 2 (the last 10 FFT levels if this is a 2560K FFT) */ \
\
/* pass2: */start_timer(0); \
xpass2name(g); \
end_timer(0); \
if(g->ffttype == 1)		/* We're done if FFTing only */ \
	return; /* goto done; */ \
\
/* Do the inverse FFT */ \
\
pass1_inverse_fft_setup(c0b, procname); \
\
/* Do inverse FFT levels 10,11 */ \
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
/* Do 320 four_complex_unfft4 operations */ \
/*	distance between fft data elements is 1K */ \
/*	do 320 macros each processing 16 data values */ \
\
c0b:	start_timer(16); \
rsi = (uintptr_t)g->data_addr;/* Get source address */ \
rcx = (uintptr_t)g->scratch_area;	/* Get address of scratch area */ \
rdi = (uintptr_t)g->premult_addr;/* Load premult/sin/cos pointer */ \
pfing(rbp = (rdi+(56-20)*2*clm*64)); /* Finish prefetching the sin/cos/premult data */ \
loops_init_prefetch(320*2*clm, 64, 16, rbp); /* 320*2*clm iterations */ \
/*bug	xtouchpages rbp, 320*2*clm * 64 / 16 ;; Prefetch TLBs */ \
for(unsigned int loopA = 16; loopA; loopA--) { \
	for(unsigned int loopB = 5; loopB; loopB--) { \
		for(unsigned int loopC = 4; loopC; loopC--) { \
			for(unsigned int loopD = 2*clm; loopD; loopD--) { \
				r4_sg2cl_four_complex_unfft4(rsi, 64, blkdst, rcx, 64, clmblkdst, rdi); \
				rdi += XMM_SCD4;		/* Next premult/sin/cos data */ \
			} \
			rsi += -2*clm*64+2*blkdst;	/* Next source pointer */ \
			rcx += -2*clm*64+2*clmblkdst; /* Next dest pointer */ \
			rdi += -2*clm*XMM_SCD4;	/* Restore premult/sin/cos data pointer */ \
		} \
		rcx += -4*2*clmblkdst+clmblkdst8; \
	} \
	rdi += 2*clm*XMM_SCD4;	/* Next premult/sin/cos pointer */ \
} \
rsi = (uintptr_t)g->scratch_area;	/* Restore source pointer */ \
end_timer(16); \
\
/* Do inverse FFT levels 8,9 */ \
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
/* Do 320 four_complex_unfft operations */ \
/*	distance between fft data elements is 4K */ \
/*	do 160 macros each processing 32 data values */ \
\
start_timer(17); \
loops_init(160*2*clm);		/* 160*2*clm iterations */ \
r4_x4cl_wpn_four_complex_djbunfft_preload; \
for(unsigned int loopA = 80; loopA; loopA--) { \
	for(unsigned int loopB = 2; loopB; loopB--) { \
		for(unsigned int loopC = 2*clm; loopC; loopC--) { \
			r4_x4cl_wpn_four_complex_djbunfft(rsi, 64, 2*clmblkdst, 4*clmblkdst, rdi, 2*2*clm*XMM_SCND2); \
			rdi += XMM_SCND2;		/* Next clm's sin/cos data */ \
		} \
		rsi += -2*clm*64+clmblkdst; /* Next source pointer */ \
	} \
	rsi += -2*clmblkdst+clmblkdst8; /* Next source pointer */ \
	rdi += -2*2*clm*XMM_SCND2;	/* Restore sine/cosine pointer */ \
} \
rsi += -80*clmblkdst8;	/* Restore source pointer */ \
rdi += 4*2*clm*XMM_SCND2;	/* Next sin/cos data pointer */ \
end_timer(17); \
\
/* Do inverse FFT levels 6,7 (actually 2.322 levels) */ \
/* On input the 64-byte cache lines hold these data values: */ \
/*	0K	8K	16K	24K */ \
/*	1K	... */ \
/*	... */ \
/*	7K */ \
/*	32K	... */ \
/*	... */ \
/* On output the 64-byte cache lines hold these data values: */ \
/*	0K	40K	80K	120K */ \
/*	1K	... */ \
/*	... */ \
/*	39K */ \
/*	160K	... */ \
/*	... */ \
\
/* Do 256 five_complex_unfft operations */ \
/*	distance between fft data elements is 16K */ \
/*	do 128 macros each processing 40 data values */ \
\
start_timer(18); \
loops_init(128*2*clm);		/* 128 iterations of 2*clm */ \
r5_x5cl_five_complex_djbunfft_preload; \
for(unsigned int loopA = 16; loopA; loopA--) { \
	for(unsigned int loopB = 8; loopB; loopB--) { \
		for(unsigned int loopC = 2*clm; loopC; loopC--) { \
			r5_x5cl_five_complex_djbunfft(rsi, 64, clmblkdst8, rdi, 8*2*clm*XMM_SCD2); \
			rdi += XMM_SCD2;		/* Next sin/cos ptr */ \
		} \
		rsi += -2*clm*64+clmblkdst; /* Next source pointer */ \
	} \
	rsi += -8*clmblkdst+5*clmblkdst8; /* Next source pointer */ \
	rdi += -8*2*clm*XMM_SCD2;	/* Restore sine/cosine pointer */ \
} \
rsi += -16*5*clmblkdst8;	/* Restore source pointer */ \
end_timer(18); \
\
/* Do inverse FFT levels 4,5 */ \
/* On input the 64-byte cache lines hold these data values: */ \
/*	0K	40K	80K	120K */ \
/*	1K	... */ \
/*	... */ \
/*	39K */ \
/*	160K	... */ \
/*	... */ \
/* On output the 64-byte cache lines hold these data values: */ \
/*	0K	160K	320K	480K */ \
/*	1K	... */ \
/*	... */ \
/*	159K */ \
/*	640K	... */ \
/*	... */ \
\
/* Do 320 four_complex_unfft operations */ \
/*	distance between fft data elements is 80K */ \
/*	do 160 macros each processing 32 data values */ \
\
start_timer(19); \
rdi = (uintptr_t)g->sincos2;/* Load fixed sin/cos data ptr */ \
rdi += 320*XMM_SCD2;	/* Sin/cos ptr */ \
loops_init(160*2*clm);		/* 160*2*clm iterations */ \
r4_x4cl_four_complex_djbunfft_preload; \
for(unsigned int loopA = 4; loopA; loopA--) { \
	for(unsigned int loopB = 5; loopB; loopB--) { \
		for(unsigned int loopC = 8; loopC; loopC--) { \
			for(unsigned int loopD = 2*clm; loopD; loopD--) { \
				r4_x4cl_four_complex_djbunfft(rsi, 64, 5*clmblkdst8, 10*clmblkdst8, rdi, 40*XMM_SCD2); \
			} \
			rsi += -2*clm*64+clmblkdst; /* Next source pointer */ \
			rdi += XMM_SCD2;		/* Next sin/cos ptr */ \
		} \
		rsi += -8*clmblkdst+clmblkdst8; /* Next source pointer */ \
	} \
	rsi += -5*clmblkdst8+20*clmblkdst8; /* Next source pointer */ \
	rdi += -40*XMM_SCD2;	/* Restore sine/cosine pointer */ \
} \
rsi += -4*20*clmblkdst8;	/* Restore source pointer */ \
end_timer(19); \
\
/* Do inverse FFT levels 1,2,3 */ \
/* On input the 64-byte cache lines hold these data values: */ \
/*	0K	160K	320K	480K */ \
/*	1K	... */ \
/*	... */ \
/*	159K */ \
/*	640K	... */ \
/*	... */ \
/* On output the 64-byte cache lines hold these data values: */ \
/*	0K	640K	1280K	1920K */ \
/*	1K	... */ \
/*	... */ \
/*	639K */ \
\
/* Do 320 four_complex_unfft operations */ \
/*	distance between fft data elements is 320K */ \
/*	do 160 macros each processing 32 data values */ \
\
start_timer(20); \
rdx = (uintptr_t)g->sincos1;/* Load fixed negacyclic premultipliers */ \
rdi = (uintptr_t)g->sincos2;/* Load fixed sin/cos postmultipliers */ \
loops_init(160*2*clm);		/* 160 iterations of 2*clm */ \
for(unsigned int loopA = 20; loopA; loopA--) { \
	for(unsigned int loopB = 8; loopB; loopB--) { \
		for(unsigned int loopC = 2*clm; loopC; loopC--) { \
			r4_x4cl_four_complex_last_djbunfft(rsi, 64, 20*clmblkdst8, 40*clmblkdst8, rdi, 160*XMM_SCD2, rdx, 160*XMM_SCD4); \
		} \
		newxloop1(clm, c1b, clmblkdst, rdi, XMM_SCD2, rdx, XMM_SCD4); /* 2*clm iterations */ \
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
copy_scratch_data_to_fft(640, clm); \
goto	endlp;			/* Skip FFT code if flag not set */ \
\
/* Fire up auxiliary threads to do pass 1 forward FFT */ \
\
/* pass1: */pass1_forward_fft_setup(b0b); \
\
/* Copy data from the FFT data area back to the scratch area */ \
\
b0b:	copy_fft_data_to_scratch(640, clm); \
\
/* Do FFT levels 1,2,3 */ \
/* */ \
/* On input the 64-byte cache lines hold these data values: */ \
/*	0K	640K	1280K	1920K */ \
/*	1K	... */ \
/*	... */ \
/*	639K */ \
/* On output the 64-byte cache lines hold these data values: */ \
/*	0K	320K	640K	960K */ \
/*	1K	... */ \
/*	... */ \
/*	319K	... */ \
/*	1280K	... */ \
/*	... */ \
\
/* Do 320 four_complex_first_fft operations */ \
/*	distance between fft data elements is 320K */ \
/*	do 320 iterations each processing 16 data values */ \
\
yes_fft:start_timer(22); \
rsi = (uintptr_t)g->scratch_area;	/* Get address of scratch area */ \
rdx = (uintptr_t)g->sincos1;/* Address of fixed negacyclic premultipliers */ \
rdi = (uintptr_t)g->sincos2;/* Load fixed sin/cos postmultipliers */ \
set_data_prefetch_ptrs;		/* Init prefetching pointers */ \
loops_init(320*2*clm);		/* 40 iters of 2 iters of 4*2*clm */ \
for(unsigned int loopA = 40; loopA; loopA--) { \
	for(unsigned int loopB = 2*4*2*clm; loopB; ) { \
		xtouch(xptr(rcx+rbp+(clm-1)*128));	/* Preload TLB */ \
		do{ \
			do{ \
				r4_x2cl_four_complex_first_djbfft_scratch(rsi, 64, 40*clmblkdst8, rdi, rdx); \
			}while(--loopB & (2*clm-1)); \
			newxloop4(clm, b1b, clmblkdst, rdi, XMM_SCD2, rdx, XMM_SCD4); /* 4*2*clm iters prefetching 2*clm cache lines */ \
		}while(loopB & (4*2*clm-1)); \
	} \
	rsi += -2*4*clmblkdst+clmblkdst8; \
} \
rsi += -40*clmblkdst8;	/* Restore source pointer */ \
pfing(rcx += -16*blkdst);		/* We just prefetched 80 blocks, correct for */ \
;					/* code below which will prefetch 16 blocks too many */ \
end_timer(22); \
\
/* Do FFT levels 4,5 */ \
/* */ \
/* On input the 64-byte cache lines hold these data values: */ \
/*	0K	320K	640K	960K */ \
/*	1K	... */ \
/*	... */ \
/*	319K	... */ \
/*	1280K	... */ \
/*	... */ \
/* On output the 64-byte cache lines hold these data values: */ \
/*	0K	80K	160K	240K */ \
/*	1K	... */ \
/*	... */ \
/*	79K	... */ \
/*	320K	... */ \
/*	... */ \
\
/* Do 320 four_complex_fft operations */ \
/*	distance between fft data elements is 80K */ \
/*	do 160 macros each processing 32 data values */ \
\
start_timer(23); \
loops_init(160*2*clm);		/* 10 iterations of 16 iters of 2*clm */ \
r4_x4cl_four_complex_djbfft_preload; \
for(unsigned int loopA = 2; loopA; loopA--) { \
	for(unsigned int loopB = 10; loopB; loopB--) { \
		for(unsigned int loopC = 8*2*clm; loopC; ) { \
			xtouch(xptr(rcx+rbp+(clm-1)*128));	/* Preload TLB */ \
			do{ \
				r4_x4cl_four_complex_djbfft(rsi, 64, 10*clmblkdst8, 20*clmblkdst8, rdi); \
			}while(--loopC & (2*clm-1)); \
			newxloop1(clm, b3b, clmblkdst, rdi, XMM_SCD2);	/* 2*clm iterations */ \
		} \
		rsi += -8*clmblkdst+clmblkdst8; /* Next source pointer */ \
	} \
	rsi += -10*clmblkdst8+40*clmblkdst8; /* Next source pointer */ \
	rdi += -80*XMM_SCD2;	/* Restore sin/cos ptr */ \
} \
rsi += -2*40*clmblkdst8;	/* Restore source pointer */ \
end_timer(23); \
\
/* Do FFT levels 6,7 (actually 2.322 levels) */ \
/* */ \
/* On input the 64-byte cache lines hold these data values: */ \
/*	0K	80K	160K	240K */ \
/*	1K	... */ \
/*	... */ \
/*	79K	... */ \
/*	320K	... */ \
/*	... */ \
/* On output the 64-byte cache lines hold these data values: */ \
/*	0K	16K	32K	48K */ \
/*	1K	... */ \
/*	... */ \
/*	15K	... */ \
/*	64K	... */ \
/*	... */ \
\
/* Do 256 five_complex_fft operations */ \
/*	distance between fft data elements is 16K */ \
/*	do 128 macros each processing 40 data values */ \
\
start_timer(24); \
rdi = (uintptr_t)g->premult_addr;/* Address of premultiplier/sin/cos table */ \
rdi += 4*2*clm*XMM_SCND2 + 16*2*clm*XMM_SCD4; \
loops_init(128*2*clm);		/* 2 iters of 8 iters of 8 of 2*clm */ \
r5_x5cl_five_complex_djbfft_preload; \
for(unsigned int loopA = 8; loopA; loopA--) { \
	for(unsigned int loopB = 2; loopB; loopB--) { \
		for(unsigned int loopC = 8*2*clm; loopC; ) { \
			xtouch(xptr(rcx+rbp+(clm-1)*128));	/* Preload TLB */ \
			do{ \
				r5_x5cl_five_complex_djbfft(rsi, 64, 2*clmblkdst8, rdi); \
				rdi += XMM_SCD2;		/* Next sine/cosine pointer */ \
			}while(--loopC & (2*clm-1)); \
			newxloop1doublepre(clm, b5b, b5a, clmblkdst); /* 2*clm iterations */ \
		} \
		rsi += -8*clmblkdst+clmblkdst8; \
	} \
	rsi += -2*clmblkdst8+10*clmblkdst8; /* Next source pointer */ \
	rdi += -16*2*clm*XMM_SCD2;	/* Restore sin/cos ptr */ \
} \
rsi += -8*10*clmblkdst8;	/* Restore source pointer */ \
end_timer(24); \
\
/* Do FFT levels 8,9 */ \
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
/* Do 320 four_complex_fft operations */ \
/*	distance between fft data elements is 4K */ \
/*	do 160 macros each processing 32 data values */ \
\
start_timer(25); \
rdi += -4*2*clm*XMM_SCND2;	/* Next sin/cos ptr */ \
loops_init(160*2*clm);		/* 40 iterations of 4 iters of 2*clm */ \
r4_x4cl_wpn_four_complex_djbfft_preload; \
for(unsigned int loopA = 40; loopA; loopA--) { \
	for(unsigned int loopB = 4*2*clm; loopB; ) { \
		xtouch(xptr(rcx+rbp+(clm-1)*128));	/* Preload TLB */ \
		do{ \
			r4_x4cl_wpn_four_complex_djbfft(rsi, 64, 4*clmblkdst, clmblkdst8, rdi); \
			rdi += XMM_SCND2;		/* Next sine/cosine pointer */ \
		}while(--loopB & (2*clm-1)); \
		newxloop1(clm, b7b, clmblkdst);	/* 2*clm iterations */ \
	} \
	rsi += -4*clmblkdst+2*clmblkdst8; /* Next source pointer */ \
	rdi += -4*2*clm*XMM_SCND2;	/* Restore sin/cos ptr */ \
} \
rsi += -40*2*clmblkdst8;	/* Restore source pointer */ \
end_timer(25); \
\
/* Do FFT levels 10,11 */ \
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
/* Do 320 four_complex_fft operations */ \
/*	distance between fft data elements is 1K */ \
/*	do 160 macros each processing 32 data values */ \
\
start_timer(26); \
rdi += -16*2*clm*XMM_SCD4;	/* Next premult/sin/cos ptr */ \
rcx = (uintptr_t)g->data_addr;/* Destination pointer */ \
loops_init_prefetch(160*2*clm, 64, 4, rbp); /* 160 iters of 2*clm */ \
pfing(rbp = (uintptr_t)g->premult_prefetch);/* Load premult/sin/cos prefetch pointer */ \
for(unsigned int loopA = 16; loopA; loopA--) { \
	for(unsigned int loopB = 5; loopB; loopB--) { \
		for(unsigned int loopC = 2; loopC; loopC--) { \
			for(unsigned int loopD = 2*clm; loopD; loopD--) { \
				r4_sg4cl_four_complex_fft4(rsi, 64, clmblkdst, 2*clmblkdst, rcx, 64, blkdst, 2*blkdst, rdi); \
				rdi += XMM_SCD4;		/* Next premult/sin/cos ptr */ \
			} \
			rsi += -2*clm*64+4*clmblkdst; /* Next source pointer */ \
			rcx += -2*clm*64+4*blkdst;	/* Next dest pointer */ \
			rdi += -2*clm*XMM_SCD4;	/* Restore premult/sine/cosine pointer */ \
		} \
		rsi += -2*4*clmblkdst+clmblkdst8; /* Next source pointer */ \
	} \
	rdi += 2*clm*XMM_SCD4;	/* Next premult/sine/cosine pointer */ \
} \
end_timer(26); \
\
/* Work on next 4*clm sets of 2560 values */ \
\
endlp:	xpass1_get_next_block(pass2, c0b, b0b, _wpn); /* Get next block and jump */ \
/* done: */ xfft_footer_rn; \
\
\
\
/* Perform a pass 1 negacyclic FFT of length 3072 doubles. */ \
/* Memory examples are for a 3M FFT with 10 levels done in pass 2. */ \
/* Sin/cos data needed: 32*2*clm*XMM_SCD1 + 8*2*clm*XMM_SCND2 + 16*2*clm*XMM_SCD8 = 96*2*clm 64-byte cache lines */ \
\

#define r4dwpn_pass1sc3072ac(procname, xpass2name, pass2_levels, clm) /* (fftname, pass2_levels, clmarg) */ \
\
/* set_FFT_constants(fftname, 0, 1, pass2_levels, clmarg); */ \
xfft_header_rn(procname, xpass2name, pass2_levels, clm); \
\
/* Do pass 2 (the last 10 FFT levels if this is a 3M FFT) */ \
\
/* pass2: */start_timer(0); \
xpass2name(g); \
end_timer(0); \
if(g->ffttype == 1)		/* We're done if FFTing only */ \
	return; /* goto done; */ \
\
/* Do the inverse FFT */ \
\
pass1_inverse_fft_setup(c0b, procname); \
\
/* Do inverse FFT levels 10,11,12 */ \
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
/* Do 192 eight_complex_unfft8 operations */ \
/*	distance between fft data elements is 1K */ \
/*	do 192 macros each processing 32 data values */ \
\
c0b:	start_timer(16); \
rsi = (uintptr_t)g->data_addr;/* Get source address */ \
rcx = (uintptr_t)g->scratch_area;	/* Get address of scratch area */ \
rdi = (uintptr_t)g->premult_addr;/* Load premult/sin/cos pointer */ \
pfing(rbp = (rdi+(96-48)*2*clm*64)); /* Finish prefetching the sin/cos/premult data */ \
loops_init_prefetch(192*2*clm, 64, 4, rbp); /* 192*2*clm iterations */ \
/*bug	xtouchpages rbp, 192*2*clm * 64 / 4 ;; Prefetch TLBs */ \
rdi += 128;		/* Generates 1 byte offsets in radix-8 macro */ \
for(unsigned int loopA = 16; loopA; loopA--) { \
	for(unsigned int loopB = 6; loopB; loopB--) { \
		for(unsigned int loopC = 2; loopC; loopC--) { \
			for(unsigned int loopD = 2*clm; loopD; loopD--) { \
				r8_sg4cl_eight_complex_unfft8(rsi, 64, blkdst, 2*blkdst, rcx, 64, clmblkdst, 2*clmblkdst, rdi-128); \
				rdi += XMM_SCD8;		/* Next premult/sin/cos data */ \
			} \
			rsi += -2*clm*64+4*blkdst;	/* Next source pointer */ \
			rcx += -2*clm*64+4*clmblkdst; /* Next dest pointer */ \
			rdi += -2*clm*XMM_SCD8;	/* Restore premult/sin/cos data pointer */ \
		} \
		rcx += -2*4*clmblkdst+clmblkdst8; \
	} \
	rdi += 2*clm*XMM_SCD8;	/* Next premult/sin/cos pointer */ \
} \
rsi = (uintptr_t)g->scratch_area;	/* Restore source pointer */ \
rdi += -128;		/* Next sin/cos data pointer */ \
end_timer(16); \
\
/* Do inverse FFT levels 8,9 */ \
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
/* Do 384 four_complex_unfft operations */ \
/*	distance between fft data elements is 8K */ \
/*	do 192 macros each processing 32 data values */ \
\
start_timer(17); \
loops_init(192*2*clm);		/* 192*2*clm iterations */ \
r4_x4cl_wpn_four_complex_djbunfft_preload; \
for(unsigned int loopA = 48; loopA; loopA--) { \
	for(unsigned int loopB = 4; loopB; loopB--) { \
		for(unsigned int loopC = 2*clm; loopC; loopC--) { \
			r4_x4cl_wpn_four_complex_djbunfft(rsi, 64, 4*clmblkdst, clmblkdst8, rdi, 4*2*clm*XMM_SCND2); \
			rdi += XMM_SCND2;		/* Next clm's sin/cos data */ \
		} \
		rsi += -2*clm*64+clmblkdst; /* Next source pointer */ \
	} \
	rsi += -4*clmblkdst+2*clmblkdst8; /* Next source pointer */ \
	rdi += -4*2*clm*XMM_SCND2;	/* Restore sine/cosine pointer */ \
} \
rsi += -48*2*clmblkdst8;	/* Restore source pointer */ \
rdi += 8*2*clm*XMM_SCND2;	/* Next sin/cos data pointer */ \
end_timer(17); \
\
/* Do inverse FFT levels 6,7 (actually 1.585 levels) */ \
/* On input the 64-byte cache lines hold these data values: */ \
/*	0K	16K	32K	48K */ \
/*	1K	... */ \
/*	... */ \
/*	15K */ \
/*	64K	... */ \
/*	... */ \
/* On output the 64-byte cache lines hold these data values: */ \
/*	0K	48K	96K	144K */ \
/*	1K	... */ \
/*	... */ \
/*	47K */ \
/*	192K	... */ \
/*	... */ \
\
/* Do 512 three_complex_unfft operations */ \
/*	distance between fft data elements is 32K */ \
/*	do 256 macros each processing 24 data values */ \
\
start_timer(18); \
loops_init(256*2*clm);		/* 256 iterations of 2*clm */ \
r3_x3cl_three_complex_djbunfft_preload; \
for(unsigned int loopA = 16; loopA; loopA--) { \
	for(unsigned int loopB = 2; loopB; loopB--) { \
		for(unsigned int loopC = 8; loopC; loopC--) { \
			for(unsigned int loopD = 2*clm; loopD; loopD--) { \
				r3_x3cl_three_complex_djbunfft(rsi, 64, 2*clmblkdst8, rdi, 16*2*clm*XMM_SCD1); \
				rdi += XMM_SCD1;		/* Next sin/cos ptr */ \
			} \
			rsi += -2*clm*64+clmblkdst; /* Next source pointer */ \
		} \
		rsi += -8*clmblkdst+clmblkdst8; /* Next source pointer */ \
	} \
	rsi += -2*clmblkdst8+6*clmblkdst8; /* Next source pointer */ \
	rdi += -16*2*clm*XMM_SCD1; /* Restore sine/cosine pointer */ \
} \
rsi += -16*6*clmblkdst8;	/* Restore source pointer */ \
end_timer(18); \
\
/* Do inverse FFT levels 4,5 */ \
/* On input the 64-byte cache lines hold these data values: */ \
/*	0K	48K	96K	144K */ \
/*	1K	... */ \
/*	... */ \
/*	47K */ \
/*	192K	... */ \
/*	... */ \
/* On output the 64-byte cache lines hold these data values: */ \
/*	0K	192K	384K	576K */ \
/*	1K	... */ \
/*	... */ \
/*	191K */ \
/*	768K	... */ \
/*	... */ \
\
/* Do 384 four_complex_unfft operations */ \
/*	distance between fft data elements is 96K */ \
/*	do 192 macros each processing 32 data values */ \
\
start_timer(19); \
rdi = (uintptr_t)g->sincos2;/* Load fixed sin/cos data ptr */ \
rdi += 384*XMM_SCD2;	/* Sin/cos ptr */ \
loops_init(192*2*clm);		/* 192 iterations of 2*clm */ \
r4_x4cl_four_complex_djbunfft_preload; \
for(unsigned int loopA = 4; loopA; loopA--) { \
	for(unsigned int loopB = 6; loopB; loopB--) { \
		for(unsigned int loopC = 8; loopC; loopC--) { \
			for(unsigned int loopD = 2*clm; loopD; loopD--) { \
				r4_x4cl_four_complex_djbunfft(rsi, 64, 6*clmblkdst8, 12*clmblkdst8, rdi, 48*XMM_SCD2); \
			} \
			rsi += -2*clm*64+clmblkdst; /* Next source pointer */ \
			rdi += XMM_SCD2;		/* Next sin/cos ptr */ \
		} \
		rsi += -8*clmblkdst+clmblkdst8; /* Next source pointer */ \
	} \
	rsi += -6*clmblkdst8+24*clmblkdst8; /* Next source pointer */ \
	rdi += -48*XMM_SCD2;	/* Restore sine/cosine pointer */ \
} \
rsi += -4*24*clmblkdst8;	/* Restore source pointer */ \
end_timer(19); \
\
/* Do inverse FFT levels 1,2,3 */ \
/* On input the 64-byte cache lines hold these data values: */ \
/*	0K	192K	384K	576K */ \
/*	1K	... */ \
/*	... */ \
/*	191K */ \
/*	768K	... */ \
/*	... */ \
/* On output the 64-byte cache lines hold these data values: */ \
/*	0K	768K	1536K	2304K */ \
/*	1K	... */ \
/*	... */ \
/*	767K */ \
\
/* Do 384 four_complex_unfft operations */ \
/*	distance between fft data elements is 384K */ \
/*	do 192 macros each processing 32 data values */ \
\
start_timer(20); \
rdx = (uintptr_t)g->sincos1;/* Load fixed negacyclic premultipliers */ \
rdi = (uintptr_t)g->sincos2;/* Load fixed sin/cos postmultipliers */ \
loops_init(192*2*clm);		/* 192 iterations of 2*clm */ \
for(unsigned int loopA = 24; loopA; loopA--) { \
	for(unsigned int loopB = 8; loopB; loopB--) { \
		for(unsigned int loopC = 2*clm; loopC; loopC--) { \
			r4_x4cl_four_complex_last_djbunfft(rsi, 64, 24*clmblkdst8, 48*clmblkdst8, rdi, 192*XMM_SCD2, rdx, 192*XMM_SCD4); \
		} \
		newxloop1(clm, c1b, clmblkdst, rdi, XMM_SCD2, rdx, XMM_SCD4); /* 2*clm iterations */ \
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
copy_scratch_data_to_fft(768, clm); \
goto	endlp;			/* Skip FFT code if flag not set */ \
\
/* Fire up auxiliary threads to do pass 1 forward FFT */ \
\
/* pass1: */pass1_forward_fft_setup(b0b); \
\
/* Copy data from the FFT data area back to the scratch area */ \
\
b0b:	copy_fft_data_to_scratch(768, clm); \
\
/* Do FFT levels 1,2,3 */ \
/* */ \
/* On input the 64-byte cache lines hold these data values: */ \
/*	0K	768K	1536K	2304K */ \
/*	1K	... */ \
/*	... */ \
/*	767K */ \
/* On output the 64-byte cache lines hold these data values: */ \
/*	0K	384K	768K	1152K */ \
/*	1K	... */ \
/*	... */ \
/*	383K	... */ \
/*	1536K	... */ \
/*	... */ \
\
/* Do 384 four_complex_first_fft operations */ \
/*	distance between fft data elements is 384K */ \
/*	do 384 iterations each processing 16 data values */ \
\
yes_fft:start_timer(22); \
rsi = (uintptr_t)g->scratch_area;	/* Get address of scratch area */ \
rdx = (uintptr_t)g->sincos1;/* Address of fixed negacyclic premultipliers */ \
rdi = (uintptr_t)g->sincos2;/* Load fixed sin/cos postmultipliers */ \
set_data_prefetch_ptrs;		/* Init prefetching pointers */ \
loops_init(384*2*clm);		/* 48 iters of 4 iters of 2*2*clm */ \
for(unsigned int loopA = 48; loopA; loopA--) { \
	for(unsigned int loopB = 4*2*2*clm; loopB; ) { \
		xtouch(xptr(rcx+rbp+(clm-1)*128));	/* Preload TLB */ \
		do{ \
			do{ \
				r4_x2cl_four_complex_first_djbfft_scratch(rsi, 64, 48*clmblkdst8, rdi, rdx); \
			}while(--loopB & (2*clm-1)); \
			newxloop2(clm, b1b, clmblkdst, rdi, XMM_SCD2, rdx, XMM_SCD4); /* 2*2*clm iters prefetching 2*clm cache lines */ \
		}while(loopB & (2*2*clm-1)); \
	} \
	rsi += -4*2*clmblkdst+clmblkdst8; \
} \
rsi += -48*clmblkdst8;	/* Restore source pointer */ \
pfing(rcx += -64*blkdst);		/* We just prefetched 192 blocks, correct for code */ \
;					/* below which will prefetch 64 blocks too many */ \
end_timer(22); \
\
/* Do FFT levels 4,5 */ \
/* */ \
/* On input the 64-byte cache lines hold these data values: */ \
/*	0K	384K	768K	1152K */ \
/*	1K	... */ \
/*	... */ \
/*	383K	... */ \
/*	1536K	... */ \
/*	... */ \
/* On output the 64-byte cache lines hold these data values: */ \
/*	0K	96K	192K	288K */ \
/*	1K	... */ \
/*	... */ \
/*	95K	... */ \
/*	384K	... */ \
/*	... */ \
\
/* Do 384 four_complex_fft operations */ \
/*	distance between fft data elements is 96K */ \
/*	do 192 macros each processing 32 data values */ \
\
start_timer(23); \
loops_init(192*2*clm);		/* 6 iterations of 32 iters of 2*clm */ \
r4_x4cl_four_complex_djbfft_preload; \
for(unsigned int loopA = 2; loopA; loopA--) { \
	for(unsigned int loopB = 12; loopB; loopB--) { \
		for(unsigned int loopC = 8*2*clm; loopC; ) { \
			xtouch(xptr(rcx+rbp+(clm-1)*128));	/* Preload TLB */ \
			do{ \
				r4_x4cl_four_complex_djbfft(rsi, 64, 12*clmblkdst8, 24*clmblkdst8, rdi); \
			}while(--loopC & (2*clm-1)); \
			newxloop1(clm, b3b, clmblkdst, rdi, XMM_SCD2);	/* 2*clm iterations */ \
		} \
		rsi += -8*clmblkdst+clmblkdst8; /* Next source pointer */ \
	} \
	rsi += -12*clmblkdst8+48*clmblkdst8; /* Next source pointer */ \
	rdi += -96*XMM_SCD2;	/* Restore sin/cos ptr */ \
} \
rsi += -2*48*clmblkdst8;	/* Restore source pointer */ \
end_timer(23); \
\
/* Do FFT levels 6,7 (actually 1.585 levels) */ \
/* */ \
/* On input the 64-byte cache lines hold these data values: */ \
/*	0K	96K	192K	288K */ \
/*	1K	... */ \
/*	... */ \
/*	95K	... */ \
/*	384K	... */ \
/*	... */ \
/* On output the 64-byte cache lines hold these data values: */ \
/*	0K	32K	64K	96K */ \
/*	1K	... */ \
/*	... */ \
/*	31K	... */ \
/*	128K	... */ \
/*	... */ \
\
/* Do 512 three_complex_fft operations */ \
/*	distance between fft data elements is 32K */ \
/*	do 256 macros each processing 24 data values */ \
\
start_timer(24); \
rdi = (uintptr_t)g->premult_addr;/* Address of premultiplier/sin/cos table */ \
rdi += 8*2*clm*XMM_SCND2 + 16*2*clm*XMM_SCD8; \
loops_init(256*2*clm);		/* 2 iters of 16 iters of 8 of 2*clm */ \
r3_x3cl_three_complex_djbfft_preload; \
for(unsigned int loopA = 8; loopA; loopA--) { \
	for(unsigned int loopB = 4; loopB; loopB--) { \
		for(unsigned int loopC = 8*2*clm; loopC; ) { \
			xtouch(xptr(rcx+rbp+(clm-1)*128));	/* Preload TLB */ \
			do{ \
				r3_x3cl_three_complex_djbfft(rsi, 64, 4*clmblkdst8, rdi); \
				rdi += XMM_SCD1;		/* Next sine/cosine pointer */ \
			}while(--loopC & (2*clm-1)); \
			newxloop1(clm, b5b, clmblkdst);	/* 2*clm iterations */ \
		} \
		rsi += -8*clmblkdst+clmblkdst8; \
	} \
	rsi += -4*clmblkdst8+12*clmblkdst8; /* Next source pointer */ \
	rdi += -32*2*clm*XMM_SCD1; /* Restore sin/cos ptr */ \
} \
rsi += -8*12*clmblkdst8;	/* Restore source pointer */ \
end_timer(24); \
\
/* Do FFT levels 8,9 */ \
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
/* Do 384 four_complex_fft operations */ \
/*	distance between fft data elements is 8K */ \
/*	do 192 macros each processing 32 data values */ \
\
start_timer(25); \
rdi += -8*2*clm*XMM_SCND2;	/* Next sin/cos ptr */ \
loops_init(192*2*clm);		/* 24 iterations of 8 iters of 2*clm */ \
r4_x4cl_wpn_four_complex_djbfft_preload; \
for(unsigned int loopA = 24; loopA; loopA--) { \
	for(unsigned int loopB = 8*2*clm; loopB; ) { \
		xtouch(xptr(rcx+rbp+(clm-1)*128));	/* Preload TLB */ \
		do{ \
			r4_x4cl_wpn_four_complex_djbfft(rsi, 64, clmblkdst8, 2*clmblkdst8, rdi); \
			rdi += XMM_SCND2;		/* Next sine/cosine pointer */ \
		}while(--loopB & (2*clm-1)); \
		newxloop1(clm, b7b, clmblkdst);	/* 2*clm iterations */ \
	} \
	rsi += -8*clmblkdst+4*clmblkdst8; /* Next source pointer */ \
	rdi += -8*2*clm*XMM_SCND2;	/* Restore sin/cos ptr */ \
} \
rsi += -24*4*clmblkdst8;	/* Restore source pointer */ \
end_timer(25); \
\
/* Do FFT levels 10,11,12 */ \
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
/* Do 192 eight_complex_fft operations */ \
/*	distance between fft data elements is 1K */ \
/*	do 96 macros each processing 64 data values */ \
\
start_timer(26); \
rdi += -16*2*clm*XMM_SCD8+128; /* Next premult/sin/cos ptr (gen 1 byte offsets in radix-8 macro) */ \
rcx = (uintptr_t)g->data_addr;/* Destination pointer */ \
loops_init_prefetch(96*2*clm, 64, 2, rbp); /* 96 iters of 2*clm */ \
pfing(rbp = (uintptr_t)g->premult_prefetch);/* Load premult/sin/cos prefetch pointer */ \
for(unsigned int loopA = 16; loopA; loopA--) { \
	for(unsigned int loopB = 6; loopB; loopB--) { \
		for(unsigned int loopC = 2*clm; loopC; loopC--) { \
			r8_sg8cl_eight_complex_fft8(rsi, 64, clmblkdst, 2*clmblkdst, 4*clmblkdst, rcx, 64, blkdst, 2*blkdst, 4*blkdst, rdi-128); \
			rdi += XMM_SCD8;		/* Next premult/sin/cos ptr */ \
		} \
		rsi += -2*clm*64+clmblkdst8; /* Next source pointer */ \
		rcx += -2*clm*64+8*blkdst;	/* Next dest pointer */ \
		rdi += -2*clm*XMM_SCD8;	/* Restore premult/sine/cosine pointer */ \
	} \
	rdi += 2*clm*XMM_SCD8;	/* Next premult/sine/cosine pointer */ \
} \
end_timer(26); \
\
/* Work on next 4*clm sets of 3072 values */ \
\
endlp:	xpass1_get_next_block(pass2, c0b, b0b, _wpn); /* Get next block and jump */ \
/* done: */ xfft_footer_rn; \
\
\
\
/* Perform a pass 1 negacyclic FFT of length 4096 doubles. */ \
/* Memory examples are for a 4M FFT with 10 levels done in pass 2. */ \
/* Sin/cos data needed: 32*2*clm*XMM_SCD2 + 8*2*clm*XMM_SCND2 + 16*2*clm*XMM_SCD8 = 112*2*clm 64-byte cache lines */ \
\

#define r4dwpn_pass1sc4096ac(procname, xpass2name, pass2_levels, clm) /* (fftname, pass2_levels, clmarg) */ \
\
/* set_FFT_constants(fftname, 0, 1, pass2_levels, clmarg); */ \
xfft_header_rn(procname, xpass2name, pass2_levels, clm); \
\
/* Do pass 2 (the last 10 FFT levels if this is a 4M FFT) */ \
\
/* pass2: */start_timer(0); \
xpass2name(g); \
end_timer(0); \
if(g->ffttype == 1)		/* We're done if FFTing only */ \
	return; /* goto done; */ \
\
/* Do the inverse FFT */ \
\
pass1_inverse_fft_setup(c0b, procname); \
\
/* Do inverse FFT levels 10,11,12 */ \
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
/* Do 256 eight_complex_unfft8 operations */ \
/*	distance between fft data elements is 1K */ \
/*	do 256 macros each processing 32 data values */ \
\
c0b:	start_timer(16); \
rsi = (uintptr_t)g->data_addr;/* Get source address */ \
rcx = (uintptr_t)g->scratch_area;	/* Get address of scratch area */ \
rdi = (uintptr_t)g->premult_addr;/* Load premult/sin/cos pointer */ \
pfing(rbp = (rdi+(112-64)*2*clm*64)); /* Finish prefetching the sin/cos/premult data */ \
loops_init_prefetch(256*2*clm, 64, 4, rbp); /* 256*2*clm iterations */ \
/*bug	xtouchpages rbp, 256*2*clm * 64 / 4 ;; Prefetch TLBs */ \
rdi += 128;		/* Generates 1 byte offsets in radix-8 macro */ \
for(unsigned int loopA = 16; loopA; loopA--) { \
	for(unsigned int loopB = 8; loopB; loopB--) { \
		for(unsigned int loopC = 2; loopC; loopC--) { \
			for(unsigned int loopD = 2*clm; loopD; loopD--) { \
				r8_sg4cl_eight_complex_unfft8(rsi, 64, blkdst, 2*blkdst, rcx, 64, clmblkdst, 2*clmblkdst, rdi-128); \
				rdi += XMM_SCD8;		/* Next premult/sin/cos data */ \
			} \
			rsi += -2*clm*64+4*blkdst;	/* Next source pointer */ \
			rcx += -2*clm*64+4*clmblkdst; /* Next dest pointer */ \
			rdi += -2*clm*XMM_SCD8;	/* Restore premult/sin/cos data pointer */ \
		} \
		rcx += -2*4*clmblkdst+clmblkdst8; \
	} \
	rdi += 2*clm*XMM_SCD8;	/* Next premult/sin/cos pointer */ \
} \
rsi = (uintptr_t)g->scratch_area;	/* Restore source pointer */ \
rdi += -128;		/* Next sin/cos data pointer */ \
end_timer(16); \
\
/* Do inverse FFT levels 8,9 */ \
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
/* Do 512 four_complex_unfft operations */ \
/*	distance between fft data elements is 8K */ \
/*	do 256 macros each processing 32 data values */ \
\
start_timer(17); \
loops_init(256*2*clm);		/* 256*2*clm iterations */ \
r4_x4cl_wpn_four_complex_djbunfft_preload; \
for(unsigned int loopA = 64; loopA; loopA--) { \
	for(unsigned int loopB = 4; loopB; loopB--) { \
		for(unsigned int loopC = 2*clm; loopC; loopC--) { \
			r4_x4cl_wpn_four_complex_djbunfft(rsi, 64, 4*clmblkdst, clmblkdst8, rdi, 4*2*clm*XMM_SCND2); \
			rdi += XMM_SCND2;		/* Next clm's sin/cos data */ \
		} \
		rsi += -2*clm*64+clmblkdst; /* Next source pointer */ \
	} \
	rsi += -4*clmblkdst+2*clmblkdst8; /* Next source pointer */ \
	rdi += -4*2*clm*XMM_SCND2;	/* Restore sine/cosine pointer */ \
} \
rsi += -64*2*clmblkdst8;	/* Restore source pointer */ \
rdi += 8*2*clm*XMM_SCND2;	/* Next sin/cos data pointer */ \
end_timer(17); \
\
/* Do inverse FFT levels 6,7 */ \
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
/*	256K	... */ \
/*	... */ \
\
/* Do 512 four_complex_unfft operations */ \
/*	distance between fft data elements is 32K */ \
/*	do 256 macros each processing 32 data values */ \
\
start_timer(18); \
loops_init(256*2*clm);		/* 256 iterations of 2*clm */ \
r4_x4cl_four_complex_djbunfft_preload; \
for(unsigned int loopA = 16; loopA; loopA--) { \
	for(unsigned int loopB = 2; loopB; loopB--) { \
		for(unsigned int loopC = 8; loopC; loopC--) { \
			for(unsigned int loopD = 2*clm; loopD; loopD--) { \
				r4_x4cl_four_complex_djbunfft(rsi, 64, 2*clmblkdst8, 4*clmblkdst8, rdi, 16*2*clm*XMM_SCD2); \
				rdi += XMM_SCD2;		/* Next clm's sin/cos data */ \
			} \
			rsi += -2*clm*64+clmblkdst; /* Next source pointer */ \
		} \
		rsi += -8*clmblkdst+clmblkdst8; /* Next source pointer */ \
	} \
	rsi += -2*clmblkdst8+8*clmblkdst8; /* Next source pointer */ \
	rdi += -2*8*2*clm*XMM_SCD2; /* Restore sine/cosine pointer */ \
} \
rsi += -16*8*clmblkdst8;	/* Restore source pointer */ \
end_timer(18); \
\
/* Do inverse FFT levels 4,5 */ \
/* On input the 64-byte cache lines hold these data values: */ \
/*	0K	64K	128K	192K */ \
/*	1K	... */ \
/*	... */ \
/*	63K */ \
/*	256K	... */ \
/*	... */ \
/* On output the 64-byte cache lines hold these data values: */ \
/*	0K	256K	512K	768K */ \
/*	1K	... */ \
/*	... */ \
/*	255K */ \
/*	1024K	... */ \
/*	... */ \
\
/* Do 512 four_complex_unfft operations */ \
/*	distance between fft data elements is 128K */ \
/*	do 256 macros each processing 32 data values */ \
\
start_timer(19); \
rdi = (uintptr_t)g->sincos2;/* Load fixed sin/cos data ptr */ \
rdi += 512*XMM_SCD2;	/* Sin/cos ptr */ \
loops_init(256*2*clm);		/* 256 iterations of 2*clm */ \
r4_x4cl_four_complex_djbunfft_preload; \
for(unsigned int loopA = 4; loopA; loopA--) { \
	for(unsigned int loopB = 8; loopB; loopB--) { \
		for(unsigned int loopC = 8; loopC; loopC--) { \
			for(unsigned int loopD = 2*clm; loopD; loopD--) { \
				r4_x4cl_four_complex_djbunfft(rsi, 64, 8*clmblkdst8, 16*clmblkdst8, rdi, 64*XMM_SCD2); \
			} \
			rsi += -2*clm*64+clmblkdst; /* Next source pointer */ \
			rdi += XMM_SCD2;		/* Next clm's sin/cos data */ \
		} \
		rsi += -8*clmblkdst+clmblkdst8; /* Next source pointer */ \
	} \
	rsi += -8*clmblkdst8+32*clmblkdst8; /* Next source pointer */ \
	rdi += -64*XMM_SCD2;	/* Restore sine/cosine pointer */ \
} \
rsi += -4*32*clmblkdst8;	/* Restore source pointer */ \
end_timer(19); \
\
/* Do inverse FFT levels 1,2,3 */ \
/* On input the 64-byte cache lines hold these data values: */ \
/*	0K	256K	512K	768K */ \
/*	1K	... */ \
/*	... */ \
/*	255K */ \
/*	1024K	... */ \
/*	... */ \
/* On output the 64-byte cache lines hold these data values: */ \
/*	0K	1024K	2048K	3072K */ \
/*	1K	... */ \
/*	... */ \
/*	1023K */ \
\
/* Do 512 four_complex_unfft operations */ \
/*	distance between fft data elements is 512K */ \
/*	do 256 macros each processing 32 data values */ \
\
start_timer(20); \
rdx = (uintptr_t)g->sincos1;/* Load fixed negacyclic premultipliers */ \
rdi = (uintptr_t)g->sincos2;/* Load fixed sin/cos postmultipliers */ \
loops_init(256*2*clm);		/* 256 iterations of 2*clm */ \
for(unsigned int loopA = 32; loopA; loopA--) { \
	for(unsigned int loopB = 8; loopB; loopB--) { \
		for(unsigned int loopC = 2*clm; loopC; loopC--) { \
			r4_x4cl_four_complex_last_djbunfft(rsi, 64, 32*clmblkdst8, 64*clmblkdst8, rdi, 256*XMM_SCD2, rdx, 256*XMM_SCD4); \
		} \
		newxloop1(clm, c1b, clmblkdst, rdi, XMM_SCD2, rdx, XMM_SCD4); /* 2*clm iterations */ \
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
copy_scratch_data_to_fft(1024, clm); \
goto	endlp;			/* Skip FFT code if flag not set */ \
\
/* Fire up auxiliary threads to do pass 1 forward FFT */ \
\
/* pass1: */pass1_forward_fft_setup(b0b); \
\
/* Copy data from the FFT data area back to the scratch area */ \
\
b0b:	copy_fft_data_to_scratch(1024, clm); \
\
/* Do FFT levels 1,2,3 */ \
/* */ \
/* On input the 64-byte cache lines hold these data values: */ \
/*	0K	1024K	2048K	3072K */ \
/*	1K	... */ \
/*	... */ \
/*	1023K */ \
/* On output the 64-byte cache lines hold these data values: */ \
/*	0K	512K	1024K	1536K */ \
/*	1K	... */ \
/*	... */ \
/*	511K	... */ \
/*	2048K	... */ \
/*	... */ \
\
/* Do 512 four_complex_first_fft operations */ \
/*	distance between fft data elements is 512K */ \
/*	do 512 iterations each processing 16 data values */ \
\
yes_fft:start_timer(22); \
rsi = (uintptr_t)g->scratch_area;	/* Get address of scratch area */ \
rdx = (uintptr_t)g->sincos1;/* Address of fixed negacyclic premultipliers */ \
rdi = (uintptr_t)g->sincos2;/* Load fixed sin/cos postmultipliers */ \
set_data_prefetch_ptrs;		/* Init prefetching pointers */ \
loops_init(512*2*clm);		/* 64 iters of 4 iters of 2*2*clm */ \
for(unsigned int loopA = 64; loopA; loopA--) { \
	for(unsigned int loopB = 4*2*2*clm; loopB; ) { \
		xtouch(xptr(rcx+rbp+(clm-1)*128));	/* Preload TLB */ \
		do{ \
			do{ \
				r4_x2cl_four_complex_first_djbfft_scratch(rsi, 64, 64*clmblkdst8, rdi, rdx); \
			}while(--loopB & (2*clm-1)); \
			newxloop2(clm, b1b, clmblkdst, rdi, XMM_SCD2, rdx, XMM_SCD4); /* 2*2*clm iters prefetching 2*clm cache lines */ \
		}while(loopB & (2*2*clm-1)); \
	} \
	rsi += -4*2*clmblkdst+clmblkdst8; \
} \
rsi += -64*clmblkdst8;	/* Restore source pointer */ \
end_timer(22); \
\
/* Do FFT levels 4,5 */ \
/* */ \
/* On input the 64-byte cache lines hold these data values: */ \
/*	0K	512K	1024K	1536K */ \
/*	1K	... */ \
/*	... */ \
/*	511K	... */ \
/*	2048K	... */ \
/*	... */ \
/* On output the 64-byte cache lines hold these data values: */ \
/*	0K	128K	256K	384K */ \
/*	1K	... */ \
/*	... */ \
/*	127K	... */ \
/*	512K	... */ \
/*	... */ \
\
/* Do 512 four_complex_fft operations */ \
/*	distance between fft data elements is 128K */ \
/*	do 256 macros each processing 32 data values */ \
\
start_timer(23); \
loops_init(256*2*clm);		/* 2 iters of 16 iters of 8 of 2*clm */ \
r4_x4cl_four_complex_djbfft_preload; \
for(unsigned int loopA = 2; loopA; loopA--) { \
	for(unsigned int loopB = 16; loopB; loopB--) { \
		for(unsigned int loopC = 8*2*clm; loopC; ) { \
			xtouch(xptr(rcx+rbp+(clm-1)*128));	/* Preload TLB */ \
			do{ \
				r4_x4cl_four_complex_djbfft(rsi, 64, 16*clmblkdst8, 32*clmblkdst8, rdi); \
			}while(--loopC & (2*clm-1)); \
			newxloop1(clm, b3b, clmblkdst, rdi, XMM_SCD2);	/* 2*clm iterations */ \
		} \
		rsi += -8*clmblkdst+clmblkdst8; \
	} \
	rsi += -16*clmblkdst8+64*clmblkdst8; /* Next source pointer */ \
	rdi += -128*XMM_SCD2;	/* Restore sin/cos ptr */ \
} \
rsi += -2*64*clmblkdst8;	/* Restore source pointer */ \
end_timer(23); \
\
/* Do FFT levels 6,7 */ \
/* */ \
/* On input the 64-byte cache lines hold these data values: */ \
/*	0K	128K	256K	384K */ \
/*	1K	... */ \
/*	... */ \
/*	127K	... */ \
/*	512K	... */ \
/*	... */ \
/* On output the 64-byte cache lines hold these data values: */ \
/*	0K	32K	64K	96K */ \
/*	1K	... */ \
/*	... */ \
/*	31K	... */ \
/*	128K	... */ \
/*	... */ \
\
/* Do 512 four_complex_fft operations */ \
/*	distance between fft data elements is 32K */ \
/*	do 256 macros each processing 32 data values */ \
\
start_timer(24); \
rdi = (uintptr_t)g->premult_addr;/* Address of premultiplier/sin/cos table */ \
rdi += 8*2*clm*XMM_SCND2 + 16*2*clm*XMM_SCD8; \
loops_init(256*2*clm);		/* 8 iters of 4 iters of 8 of 2*clm */ \
r4_x4cl_four_complex_djbfft_preload; \
for(unsigned int loopA = 8; loopA; loopA--) { \
	for(unsigned int loopB = 4; loopB; loopB--) { \
		for(unsigned int loopC = 8*2*clm; loopC; ) { \
			xtouch(xptr(rcx+rbp+(clm-1)*128));	/* Preload TLB */ \
			do{ \
				r4_x4cl_four_complex_djbfft(rsi, 64, 4*clmblkdst8, 8*clmblkdst8, rdi); \
				rdi += XMM_SCD2;		/* Next sine/cosine pointer */ \
			}while(--loopC & (2*clm-1)); \
			newxloop1(clm, b5b, clmblkdst);	/* 2*clm iterations */ \
		} \
		rsi += -8*clmblkdst+clmblkdst8; \
	} \
	rsi += -4*clmblkdst8+16*clmblkdst8; /* Next source pointer */ \
	rdi += -32*2*clm*XMM_SCD2;	/* Restore sin/cos ptr */ \
} \
rsi += -8*16*clmblkdst8;	/* Restore source pointer */ \
end_timer(24); \
\
/* Do FFT levels 8,9 */ \
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
/* Do 512 four_complex_fft operations */ \
/*	distance between fft data elements is 8K */ \
/*	do 256 macros each processing 32 data values */ \
\
start_timer(25); \
rdi += -8*2*clm*XMM_SCND2;	/* Next sin/cos ptr */ \
loops_init(256*2*clm);		/* 32 iterations of 8 iters of 2*clm */ \
r4_x4cl_wpn_four_complex_djbfft_preload; \
for(unsigned int loopA = 32; loopA; loopA--) { \
	for(unsigned int loopB = 8*2*clm; loopB; ) { \
		xtouch(xptr(rcx+rbp+(clm-1)*128));	/* Preload TLB */ \
		do{ \
			r4_x4cl_wpn_four_complex_djbfft(rsi, 64, clmblkdst8, 2*clmblkdst8, rdi); \
			rdi += XMM_SCND2;		/* Next sine/cosine pointer */ \
		}while(--loopB & (2*clm-1)); \
		newxloop1(clm, b7b, clmblkdst);	/* 2*clm iterations */ \
	} \
	rsi += -8*clmblkdst+4*clmblkdst8; /* Next source pointer */ \
	rdi += -8*2*clm*XMM_SCND2;	/* Restore sin/cos ptr */ \
} \
rsi += -32*4*clmblkdst8;	/* Restore source pointer */ \
end_timer(25); \
\
/* Do FFT levels 10,11,12 */ \
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
/* Do 256 eight_complex_fft operations */ \
/*	distance between fft data elements is 1K */ \
/*	do 128 macros each processing 64 data values */ \
\
start_timer(26); \
rdi += -16*2*clm*XMM_SCD8+128; /* Next premult/sin/cos ptr (gen 1 byte offsets in radix-8 macro) */ \
rcx = (uintptr_t)g->data_addr;/* Destination pointer */ \
loops_init_prefetch(128*2*clm, 64, 2, rbp); /* 128 iters of 2*clm */ \
pfing(rbp = (uintptr_t)g->premult_prefetch);/* Load premult/sin/cos prefetch pointer */ \
for(unsigned int loopA = 16; loopA; loopA--) { \
	for(unsigned int loopB = 8; loopB; loopB--) { \
		for(unsigned int loopC = 2*clm; loopC; loopC--) { \
			r8_sg8cl_eight_complex_fft8(rsi, 64, clmblkdst, 2*clmblkdst, 4*clmblkdst, rcx, 64, blkdst, 2*blkdst, 4*blkdst, rdi-128); \
			rdi += XMM_SCD8;		/* Next premult/sin/cos ptr */ \
		} \
		rsi += -2*clm*64+clmblkdst8; /* Next source pointer */ \
		rcx += -2*clm*64+8*blkdst;	/* Next dest pointer */ \
		rdi += -2*clm*XMM_SCD8;	/* Restore premult/sine/cosine pointer */ \
	} \
	rdi += 2*clm*XMM_SCD8;	/* Next premult/sine/cosine pointer */ \
} \
end_timer(26); \
\
/* Work on next 4*clm sets of 4096 values */ \
\
endlp:	xpass1_get_next_block(pass2, c0b, b0b, _wpn); /* Get next block and jump */ \
/* done: */ xfft_footer_rn; \
\
\
\
/* Perform a pass 1 negacyclic FFT of length 5120 doubles. */ \
/* Memory examples are for a 5M FFT with 10 levels done in pass 2. */ \
/* Sin/cos data needed: 32*2*clm*XMM_SCD2 + 8*2*clm*XMM_SCND2 + 16*2*clm*XMM_SCD8 = 112*2*clm 64-byte cache lines */ \
\

#define r4dwpn_pass1sc5120ac(procname, xpass2name, pass2_levels, clm) /* (fftname, pass2_levels, clmarg) */ \
\
/* set_FFT_constants(fftname, 0, 1, pass2_levels, clmarg); */ \
xfft_header_rn(procname, xpass2name, pass2_levels, clm); \
\
/* Do pass 2 (the last 10 FFT levels if this is a 5M FFT) */ \
\
/* pass2: */start_timer(0); \
xpass2name(g); \
end_timer(0); \
if(g->ffttype == 1)		/* We're done if FFTing only */ \
	return; /* goto done; */ \
\
/* Do the inverse FFT */ \
\
pass1_inverse_fft_setup(c0b, procname); \
\
/* Do inverse FFT levels 10,11,12 */ \
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
/* Do 320 eight_complex_unfft8 operations */ \
/*	distance between fft data elements is 1K */ \
/*	do 320 macros each processing 32 data values */ \
\
c0b:	start_timer(16); \
rsi = (uintptr_t)g->data_addr;/* Get source address */ \
rcx = (uintptr_t)g->scratch_area;	/* Get address of scratch area */ \
rdi = (uintptr_t)g->premult_addr;/* Load premult/sin/cos pointer */ \
pfing(rbp = (rdi+(112-40)*2*clm*64)); /* Finish prefetching the sin/cos/premult data */ \
loops_init_prefetch(320*2*clm, 64, 8, rbp); /* 320*2*clm iterations */ \
/*bug	xtouchpages rbp, 320*2*clm * 64 / 8 ;; Prefetch TLBs */ \
rdi += 128;		/* Generates 1 byte offsets in radix-8 macro */ \
for(unsigned int loopA = 16; loopA; loopA--) { \
	for(unsigned int loopB = 10; loopB; loopB--) { \
		for(unsigned int loopC = 2; loopC; loopC--) { \
			for(unsigned int loopD = 2*clm; loopD; loopD--) { \
				r8_sg4cl_eight_complex_unfft8(rsi, 64, blkdst, 2*blkdst, rcx, 64, clmblkdst, 2*clmblkdst, rdi-128); \
				rdi += XMM_SCD8;		/* Next premult/sin/cos data */ \
			} \
			rsi += -2*clm*64+4*blkdst;	/* Next source pointer */ \
			rcx += -2*clm*64+4*clmblkdst; /* Next dest pointer */ \
			rdi += -2*clm*XMM_SCD8;	/* Restore premult/sin/cos data pointer */ \
		} \
		rcx += -2*4*clmblkdst+clmblkdst8; \
	} \
	rdi += 2*clm*XMM_SCD8;	/* Next premult/sin/cos pointer */ \
} \
rsi = (uintptr_t)g->scratch_area;	/* Restore source pointer */ \
rdi += -128;		/* Next sin/cos data pointer */ \
end_timer(16); \
\
/* Do inverse FFT levels 8,9 */ \
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
/* Do 640 four_complex_unfft operations */ \
/*	distance between fft data elements is 8K */ \
/*	do 320 macros each processing 32 data values */ \
\
start_timer(17); \
loops_init(320*2*clm);		/* 320 iterations of 2*clm */ \
r4_x4cl_wpn_four_complex_djbunfft_preload; \
for(unsigned int loopA = 80; loopA; loopA--) { \
	for(unsigned int loopB = 4; loopB; loopB--) { \
		for(unsigned int loopC = 2*clm; loopC; loopC--) { \
			r4_x4cl_wpn_four_complex_djbunfft(rsi, 64, 4*clmblkdst, clmblkdst8, rdi, 4*2*clm*XMM_SCND2); \
			rdi += XMM_SCND2;		/* Next clm's sin/cos data */ \
		} \
		rsi += -2*clm*64+clmblkdst; /* Next source pointer */ \
	} \
	rsi += -4*clmblkdst+2*clmblkdst8; /* Next source pointer */ \
	rdi += -4*2*clm*XMM_SCND2;	/* Restore sine/cosine pointer */ \
} \
rsi += -80*2*clmblkdst8;	/* Restore source pointer */ \
rdi += 8*2*clm*XMM_SCND2;	/* Next sin/cos data pointer */ \
end_timer(17); \
\
/* Do inverse FFT levels 6,7 (actually 2.32 levels) */ \
/* On input the 64-byte cache lines hold these data values: */ \
/*	0K	16K	32K	48K */ \
/*	1K	... */ \
/*	... */ \
/*	15K */ \
/*	64K	... */ \
/*	... */ \
/* On output the 64-byte cache lines hold these data values: */ \
/*	0K	80K	160K	240K */ \
/*	1K	... */ \
/*	... */ \
/*	79K */ \
/*	320K	... */ \
/*	... */ \
\
/* Do 512 five_complex_unfft operations */ \
/*	distance between fft data elements is 32K */ \
/*	do 256 macros each processing 40 data values */ \
\
start_timer(18); \
loops_init(256*2*clm);		/* 256 iterations of 2*clm */ \
r5_x5cl_five_complex_djbunfft_preload; \
for(unsigned int loopA = 16; loopA; loopA--) { \
	for(unsigned int loopB = 2; loopB; loopB--) { \
		for(unsigned int loopC = 8; loopC; loopC--) { \
			for(unsigned int loopD = 2*clm; loopD; loopD--) { \
				r5_x5cl_five_complex_djbunfft(rsi, 64, 2*clmblkdst8, rdi, 16*2*clm*XMM_SCD2); \
				rdi += XMM_SCD2;		/* Next clm's sin/cos data */ \
			} \
			rsi += -2*clm*64+clmblkdst; /* Next source pointer */ \
		} \
		rsi += -8*clmblkdst+clmblkdst8; /* Next source pointer */ \
	} \
	rsi += -2*clmblkdst8+10*clmblkdst8; /* Next source pointer */ \
	rdi += -16*2*clm*XMM_SCD2; /* Restore sine/cosine pointer */ \
} \
rsi += -16*10*clmblkdst8;	/* Restore source pointer */ \
end_timer(18); \
\
/* Do inverse FFT levels 4,5 */ \
/* On input the 64-byte cache lines hold these data values: */ \
/*	0K	80K	160K	240K */ \
/*	1K	... */ \
/*	... */ \
/*	79K */ \
/*	320K	... */ \
/*	... */ \
/* On output the 64-byte cache lines hold these data values: */ \
/*	0K	320K	640K	960K */ \
/*	1K	... */ \
/*	... */ \
/*	319K */ \
/*	1280K	... */ \
/*	... */ \
\
/* Do 640 four_complex_unfft operations */ \
/*	distance between fft data elements is 160K */ \
/*	do 320 macros each processing 32 data values */ \
\
start_timer(19); \
rdi = (uintptr_t)g->sincos2;/* Load fixed sin/cos data ptr */ \
rdi += 640*XMM_SCD2;	/* Sin/cos ptr */ \
loops_init(320*2*clm);		/* 320 iterations of 2*clm */ \
r4_x4cl_four_complex_djbunfft_preload; \
for(unsigned int loopA = 4; loopA; loopA--) { \
	for(unsigned int loopB = 10; loopB; loopB--) { \
		for(unsigned int loopC = 8; loopC; loopC--) { \
			for(unsigned int loopD = 2*clm; loopD; loopD--) { \
				r4_x4cl_four_complex_djbunfft(rsi, 64, 10*clmblkdst8, 20*clmblkdst8, rdi, 80*XMM_SCD2); \
			} \
			rsi += -2*clm*64+clmblkdst; /* Next source pointer */ \
			rdi += XMM_SCD2;		/* Next sin/cos ptr */ \
		} \
		rsi += -8*clmblkdst+clmblkdst8; /* Next source pointer */ \
	} \
	rsi += -10*clmblkdst8+40*clmblkdst8; /* Next source pointer */ \
	rdi += -80*XMM_SCD2;	/* Restore sine/cosine pointer */ \
} \
rsi += -4*40*clmblkdst8;	/* Restore source pointer */ \
end_timer(19); \
\
/* Do inverse FFT levels 1,2,3 */ \
/* On input the 64-byte cache lines hold these data values: */ \
/*	0K	320K	640K	960K */ \
/*	1K	... */ \
/*	... */ \
/*	319K */ \
/*	1280K	... */ \
/*	... */ \
/* On output the 64-byte cache lines hold these data values: */ \
/*	0K	1280K	2560K	3840K */ \
/*	1K	... */ \
/*	... */ \
/*	1279K */ \
\
/* Do 640 four_complex_unfft operations */ \
/*	distance between fft data elements is 640K */ \
/*	do 320 macros each processing 32 data values */ \
\
start_timer(20); \
rdx = (uintptr_t)g->sincos1;/* Load fixed negacyclic premultipliers */ \
rdi = (uintptr_t)g->sincos2;/* Load fixed sin/cos postmultipliers */ \
loops_init(320*2*clm);		/* 320 iterations of 2*clm */ \
for(unsigned int loopA = 40; loopA; loopA--) { \
	for(unsigned int loopB = 8; loopB; loopB--) { \
		for(unsigned int loopC = 2*clm; loopC; loopC--) { \
			r4_x4cl_four_complex_last_djbunfft(rsi, 64, 40*clmblkdst8, 80*clmblkdst8, rdi, 320*XMM_SCD2, rdx, 320*XMM_SCD4); \
		} \
		newxloop1(clm, c1b, clmblkdst, rdi, XMM_SCD2, rdx, XMM_SCD4); /* 2*clm iterations */ \
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
copy_scratch_data_to_fft(1280, clm); \
goto	endlp;			/* Skip FFT code if flag not set */ \
\
/* Fire up auxiliary threads to do pass 1 forward FFT */ \
\
/* pass1: */pass1_forward_fft_setup(b0b); \
\
/* Copy data from the FFT data area back to the scratch area */ \
\
b0b:	copy_fft_data_to_scratch(1280, clm); \
\
/* Do FFT levels 1,2,3 */ \
/* */ \
/* On input the 64-byte cache lines hold these data values: */ \
/*	0K	1280K	2560K	3840K */ \
/*	1K	... */ \
/*	... */ \
/*	1279K */ \
/* On output the 64-byte cache lines hold these data values: */ \
/*	0K	640K	1280K	1920K */ \
/*	1K	... */ \
/*	... */ \
/*	639K	... */ \
/*	2560K	... */ \
/*	... */ \
\
/* Do 640 four_complex_first_fft operations */ \
/*	distance between fft data elements is 640K */ \
/*	do 640 iterations each processing 16 data values */ \
\
yes_fft:start_timer(22); \
rsi = (uintptr_t)g->scratch_area;	/* Get address of scratch area */ \
rdx = (uintptr_t)g->sincos1;/* Address of fixed negacyclic premultipliers */ \
rdi = (uintptr_t)g->sincos2;/* Load fixed sin/cos postmultipliers */ \
set_data_prefetch_ptrs;		/* Init prefetching pointers */ \
loops_init(640*2*clm);		/* 80 iters of 2 iters of 4*2*clm */ \
for(unsigned int loopA = 80; loopA; loopA--) { \
	for(unsigned int loopB = 2*4*2*clm; loopB; ) { \
		xtouch(xptr(rcx+rbp+(clm-1)*128));	/* Preload TLB */ \
		do{ \
			do{ \
				r4_x2cl_four_complex_first_djbfft_scratch(rsi, 64, 80*clmblkdst8, rdi, rdx); \
			}while(--loopB & (2*clm-1)); \
			newxloop4(clm, b1b, clmblkdst, rdi, XMM_SCD2, rdx, XMM_SCD4); /* 4*2*clm iters prefetching 2*clm cache lines */ \
		}while(loopB & (4*2*clm-1)); \
	} \
	rsi += -2*4*clmblkdst+clmblkdst8; \
} \
rsi += -80*clmblkdst8;	/* Restore source pointer */ \
pfing(rcx += -32*blkdst);		/* We just prefetched 160 blocks, corect for code */ \
;					/* below which will prefetch 32 blocks too many */ \
end_timer(22); \
\
/* Do FFT levels 4,5 */ \
/* */ \
/* On input the 64-byte cache lines hold these data values: */ \
/*	0K	640K	1280K	1920K */ \
/*	1K	... */ \
/*	... */ \
/*	639K	... */ \
/*	2560K	... */ \
/*	... */ \
/* On output the 64-byte cache lines hold these data values: */ \
/*	0K	160K	320K	480K */ \
/*	1K	... */ \
/*	... */ \
/*	159K	... */ \
/*	640K	... */ \
/*	... */ \
\
/* Do 640 four_complex_fft operations */ \
/*	distance between fft data elements is 160K */ \
/*	do 320 macros each processing 32 data values */ \
\
start_timer(23); \
loops_init(320*2*clm);		/* 10 iterations of 32 iters of 2*clm */ \
r4_x4cl_four_complex_djbfft_preload; \
for(unsigned int loopA = 2; loopA; loopA--) { \
	for(unsigned int loopB = 20; loopB; loopB--) { \
		for(unsigned int loopC = 8*2*clm; loopC; ) { \
			xtouch(xptr(rcx+rbp+(clm-1)*128));	/* Preload TLB */ \
			do{ \
				r4_x4cl_four_complex_djbfft(rsi, 64, 20*clmblkdst8, 40*clmblkdst8, rdi); \
			}while(--loopC & (2*clm-1)); \
			newxloop1(clm, b3b, clmblkdst, rdi, XMM_SCD2);	/* 2*clm iterations */ \
		} \
		rsi += -8*clmblkdst+clmblkdst8; /* Next source pointer */ \
	} \
	rsi += -20*clmblkdst8+80*clmblkdst8; /* Next source pointer */ \
	rdi += -160*XMM_SCD2;	/* Restore sin/cos ptr */ \
} \
rsi += -2*80*clmblkdst8;	/* Restore source pointer */ \
end_timer(23); \
\
/* Do FFT levels 6,7 (actually 2.32 levels) */ \
/* */ \
/* On input the 64-byte cache lines hold these data values: */ \
/*	0K	160K	320K	480K */ \
/*	1K	... */ \
/*	... */ \
/*	159K	... */ \
/*	640K	... */ \
/*	... */ \
/* On output the 64-byte cache lines hold these data values: */ \
/*	0K	32K	64K	96K */ \
/*	1K	... */ \
/*	... */ \
/*	31K	... */ \
/*	128K	... */ \
/*	... */ \
\
/* Do 512 five_complex_fft operations */ \
/*	distance between fft data elements is 32K */ \
/*	do 256 macros each processing 40 data values */ \
\
start_timer(24); \
rdi = (uintptr_t)g->premult_addr;/* Address of premultiplier/sin/cos table */ \
rdi += 8*2*clm*XMM_SCND2 + 16*2*clm*XMM_SCD8; \
loops_init(256*2*clm);		/* 2 iters of 16 iters of 8 of 2*clm */ \
r5_x5cl_five_complex_djbfft_preload; \
for(unsigned int loopA = 8; loopA; loopA--) { \
	for(unsigned int loopB = 4; loopB; loopB--) { \
		for(unsigned int loopC = 8*2*clm; loopC; ) { \
			xtouch(xptr(rcx+rbp+(clm-1)*128));	/* Preload TLB */ \
			do{ \
				r5_x5cl_five_complex_djbfft(rsi, 64, 4*clmblkdst8, rdi); \
				rdi += XMM_SCD2;		/* Next sine/cosine pointer */ \
			}while(--loopC & (2*clm-1)); \
			newxloop1doublepre(clm, b5b, b5a, clmblkdst); /* 2*clm iterations */ \
		} \
		rsi += -8*clmblkdst+clmblkdst8; \
	} \
	rsi += -4*clmblkdst8+20*clmblkdst8; /* Next source pointer */ \
	rdi += -32*2*clm*XMM_SCD2; /* Restore sin/cos ptr */ \
} \
rsi += -8*20*clmblkdst8;	/* Restore source pointer */ \
end_timer(24); \
\
/* Do FFT levels 8,9 */ \
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
/* Do 640 four_complex_fft operations */ \
/*	distance between fft data elements is 8K */ \
/*	do 320 macros each processing 32 data values */ \
\
start_timer(25); \
rdi += -8*2*clm*XMM_SCND2;	/* Next sin/cos ptr */ \
loops_init(320*2*clm);		/* 40 iterations of 8 iters of 2*clm */ \
r4_x4cl_wpn_four_complex_djbfft_preload; \
for(unsigned int loopA = 40; loopA; loopA--) { \
	for(unsigned int loopB = 8*2*clm; loopB; ) { \
		xtouch(xptr(rcx+rbp+(clm-1)*128));	/* Preload TLB */ \
		do{ \
			r4_x4cl_wpn_four_complex_djbfft(rsi, 64, clmblkdst8, 2*clmblkdst8, rdi); \
			rdi += XMM_SCND2;		/* Next sine/cosine pointer */ \
		}while(--loopB & (2*clm-1)); \
		newxloop1(clm, b7b, clmblkdst);	/* 2*clm iterations */ \
	} \
	rsi += -8*clmblkdst+4*clmblkdst8; /* Next source pointer */ \
	rdi += -8*2*clm*XMM_SCND2;	/* Restore sin/cos ptr */ \
} \
rsi += -40*4*clmblkdst8;	/* Restore source pointer */ \
end_timer(25); \
\
/* Do FFT levels 10,11,12 */ \
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
/* Do 320 eight_complex_fft operations */ \
/*	distance between fft data elements is 1K */ \
/*	do 160 macros each processing 64 data values */ \
\
start_timer(26); \
rdi += -16*2*clm*XMM_SCD8+128; /* Next premult/sin/cos ptr (gen 1 byte offsets in radix-8 macro) */ \
rcx = (uintptr_t)g->data_addr;/* Destination pointer */ \
loops_init_prefetch(160*2*clm, 64, 2, rbp); /* 160 iters of 2*clm */ \
pfing(rbp = (uintptr_t)g->premult_prefetch);/* Load premult/sin/cos prefetch pointer */ \
for(unsigned int loopA = 16; loopA; loopA--) { \
	for(unsigned int loopB = 10; loopB; loopB--) { \
		for(unsigned int loopC = 2*clm; loopC; loopC--) { \
			r8_sg8cl_eight_complex_fft8(rsi, 64, clmblkdst, 2*clmblkdst, 4*clmblkdst, rcx, 64, blkdst, 2*blkdst, 4*blkdst, rdi-128); \
			rdi += XMM_SCD8;		/* Next premult/sin/cos ptr */ \
		} \
		rsi += -2*clm*64+clmblkdst8; /* Next source pointer */ \
		rcx += -2*clm*64+8*blkdst;	/* Next dest pointer */ \
		rdi += -2*clm*XMM_SCD8;	/* Restore premult/sine/cosine pointer */ \
	} \
	rdi += 2*clm*XMM_SCD8;	/* Next premult/sine/cosine pointer */ \
} \
end_timer(26); \
\
/* Work on next 4*clm sets of 5120 values */ \
\
endlp:	xpass1_get_next_block(pass2, c0b, b0b, _wpn); /* Get next block and jump */ \
/* done: */ xfft_footer_rn; \
\
