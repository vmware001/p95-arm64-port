#pragma once

typedef void (*NORMRTN_pass1)(struct gwasm_data *__restrict g);

#ifndef FROM_ASM
void xgw_carries(struct gwasm_data *__restrict g);
#endif

// xfft_1_ret and xfft_3_ret: used in one pass

#define xfft_1_ret

#if 0 // To asm code.
#define xfft_3_ret { \
register struct gwasm_data *r11_arg asm("r11") = g; \
register vec2f64 xmm7_arg asm("xmm7") = xmm7; \
__asm__ volatile ( \
	".intel_syntax noprefix\n" \
	"call qword ptr [r11 + 0x18]\n" \
	".att_syntax\n" \
	: : /*"r" (g->NORMRTN),*/ "r" (r11_arg), "x" (xmm7_arg) \
	: "rax", "rcx", "rdx", "rbx", "rsi", "rdi", "rbp", \
	"r8", "r9", "r10", "memory",  \
	"xmm0", "xmm1", "xmm2", "xmm3", "xmm4", "xmm5", "xmm6", "xmm7", \
	"xmm8", "xmm9", "xmm10", "xmm11", "xmm12", "xmm13", "xmm14", "xmm15" ); \
}
#else  // To converted C code.
#define xfft_3_ret (*(NORMRTN_pass1)g->NORMRTN)(g);
#endif

#define set_data_prefetch_ptrs

// Used in "hg" (home-grown) FFTs
#define XMM_SCD	96			// Sizeof an SSE2 sin/cos table entry
#define XMM_PMD 128			// Sizeof an entry in premultiplier table

// Used in r4 (traditional radix-4) FFTs
#define XMM_SCD1  32	// Sizeof an SSE2 small sin/cos table entry
#define XMM_SCD2  64	// Sizeof an SSE2 DJB sin/cos table entry
#define XMM_SCD3  96	// Sizeof an eight_reals sin/cos table entry
#define XMM_SCD4 128	// Sizeof an entry in premultiplier table
#define XMM_SCD8 256	// Sizeof an radix-8 with premultipliers table entry
#define XMM_SCD9 288	// Sizeof a real radix-20 with premultipliers table entry
#define XMM_SCD13 416	// Sizeof a real radix-28 with premultipliers table entry
#define XMM_HSCD1 16	// Sizeof a half sin/cos table entry used in real FFTs in pass 2
#define XMM_HSCD2 32	// Sizeof a half sin/cos table entry used in real FFTs in pass 2

// Used in r4dwpn (radix-4 delay with partial normalization) FFTs
#define XMM_SCND2 128	// Sizeof a complex radix-4 or radix-5 sin/cos/normalization table entry
#define XMM_SCND2R 96	// Sizeof a real radix-4 or radix-5 sin/cos/normalization table entry
#define XMM_GMD 160	// Sizeof a grp two-to-phi/two-to-minus-phi multiplier block in a r4dwpn FFT

#define dist128 (128*64+128)			/* Distance between 128 cache lines */



#if 0
#define xpass1_get_next_block(pass2, c0b, b0b, suffix) \
do{ \
	int a = g->pass1_get_next_block(g); \
	\
	if(a < 1) goto c0b; /* Do another inverse FFT/norm/FFT */ \
	if(a == 1) goto b0b; /* Do another forward FFT */ \
	if(a < 4) break; /* Pass 1 main thread done or pass 1 auxiliary thread done */ \
	if(a == 4) goto pass2; /* Start pass 2 */ \
	\
	/* Split the carries into high and low carries (so they do not exceed the */ \
	/* desired number of bits per FFT word.  Then add the carries back to the */ \
	/* FFT data. */ \
	\
	start_timer(29); \
	/* call */ xgw_carries##suffix(g); \
	end_timer(29); \
}while(1); \
\
/* All done, end timer and fall through to exit code */ \
\
end_timer(1)
#endif


#define TEST2EXEC0(A,B) A
#define TEST2EXEC1(A,B) B

#define pass1_normalize(scratch, do_forward_fft) { \
	pass1_pre_carries(g); \
	start_timer(28); \
	rsi = (uintptr_t) TEST2EXEC##scratch(g->data_addr, g->scratch_area); \
	((c_func1_uptr)g->NORMRTN)(g, rsi); \
	end_timer(28); \
	uintptr_t a = pass1_post_carries(g); \
	if(a & 0xFF) \
		goto do_forward_fft; \
}





/* */ \
/* Macros to copy data to and from scratch area */ \
/* */ \
\

#define copy_scratch_data_to_fft(cachelines, clm) \
start_timer(30); \
rsi = (uintptr_t)g->scratch_area;	/* Get address of scratch area */ \
rcx = (uintptr_t)g->data_addr;/* Load source address */ \
for(unsigned int loopZ = cachelines/8; loopZ; loopZ--) { \
	for(unsigned int loopY = 8; loopY; loopY--) { \
		for(unsigned int loopX = 2*clm; loopX; loopX--) { \
			xmm0 = xptr(rsi);		/* Copy a cache line */ \
			xmm1 = xptr(rsi+16); \
			xmm2 = xptr(rsi+32); \
			xmm3 = xptr(rsi+48); \
			xptr(rcx) = xmm0; \
			xptr(rcx+16) = xmm1; \
			xptr(rcx+32) = xmm2; \
			xptr(rcx+48) = xmm3; \
			rsi += 64;			/* Next source pointer */ \
			rcx += 64;			/* Next dest pointer */ \
		} \
		rsi += -2*clm*64+clmblkdst;/* Next source pointer */ \
		rcx += -2*clm*64+blkdst;	/* Next dest pointer */ \
	} \
	rsi += -8*clmblkdst+clmblkdst8; \
} \
end_timer(30); \
\
\

#define copy_scratch_no_pad_data_to_fft(cachelines, clm) \
start_timer(30); \
rsi = (uintptr_t)g->scratch_area;	/* Get address of scratch area */ \
rcx = (uintptr_t)g->data_addr;/* Load source address */ \
for(unsigned int loopZ = cachelines; loopZ; loopZ--) { \
	for(unsigned int loopY = 2*clm; loopY; loopY--) { \
		xmm0 = xptr(rsi);		/* Copy a cache line */ \
		xmm1 = xptr(rsi+16); \
		xmm2 = xptr(rsi+32); \
		xmm3 = xptr(rsi+48); \
		xptr(rcx) = xmm0; \
		xptr(rcx+16) = xmm1; \
		xptr(rcx+32) = xmm2; \
		xptr(rcx+48) = xmm3; \
		rsi += 64;			/* Next source pointer */ \
		rcx += 64;			/* Next dest pointer */ \
	} \
	rsi += -2*clm*64+clmblkdst;/* Next source pointer */ \
	rcx += -2*clm*64+blkdst;	/* Next dest pointer */ \
} \
end_timer(30); \
\
\

#define copy_fft_data_to_scratch(cachelines, clm) \
start_timer(31); \
rsi = (uintptr_t)g->data_addr;/* Load source address */ \
rsi += g->DIST_TO_FFTSRCARG; \
rcx = (uintptr_t)g->scratch_area;	/* Get address of scratch area */ \
for(unsigned int loopZ = cachelines/8; loopZ; loopZ--) { \
	for(unsigned int loopY = 8; loopY; loopY--) { \
		for(unsigned int loopX = 2*clm; loopX; loopX--) { \
			xmm0 = xptr(rsi);		/* Copy a cache line */ \
			xmm1 = xptr(rsi+16); \
			xmm2 = xptr(rsi+32); \
			xmm3 = xptr(rsi+48); \
			xptr(rcx) = xmm0; \
			xptr(rcx+16) = xmm1; \
			xptr(rcx+32) = xmm2; \
			xptr(rcx+48) = xmm3; \
			rsi += 64;			/* Next source pointer */ \
			rcx += 64;			/* Next dest pointer */ \
		} \
		rsi += -2*clm*64+blkdst;	/* Next source pointer */ \
		rcx += -2*clm*64+clmblkdst;/* Next dest pointer */ \
	} \
	rcx += -8*clmblkdst+clmblkdst8; \
} \
end_timer(31); \
\
\

#define copy_fft_data_to_scratch_no_pad(cachelines, clm) \
start_timer(31); \
rsi = (uintptr_t)g->data_addr;/* Load source address */ \
rsi += g->DIST_TO_FFTSRCARG; \
rcx = (uintptr_t)g->scratch_area;	/* Get address of scratch area */ \
for(unsigned int loopZ = cachelines; loopZ; loopZ--) { \
	for(unsigned int loopY = 2*clm; loopY; loopY--) { \
		xmm0 = xptr(rsi);		/* Copy a cache line */ \
		xmm1 = xptr(rsi+16); \
		xmm2 = xptr(rsi+32); \
		xmm3 = xptr(rsi+48); \
		xptr(rcx) = xmm0; \
		xptr(rcx+16) = xmm1; \
		xptr(rcx+32) = xmm2; \
		xptr(rcx+48) = xmm3; \
		rsi += 64;			/* Next source pointer */ \
		rcx += 64;			/* Next dest pointer */ \
	} \
	rsi += -2*clm*64+blkdst;	/* Next source pointer */ \
	rcx += -2*clm*64+clmblkdst;/* Next dest pointer */ \
} \
end_timer(31); \
\

