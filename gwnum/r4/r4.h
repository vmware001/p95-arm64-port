#pragma once
/* Copyright 2009-2024 - Mersenne Research, Inc.  All rights reserved */ \
/* Author:  George Woltman */ \
/* Email: woltman@alum.mit.edu */ \
/* */ \
\

/* */ \
/* */ \
/* All new macros for version 26 of gwnum where we do a very traditional, primarily */ \
/* radix-4, FFT.  The forward FFT macros multiply by the sin/cos values at the end */ \
/* of the macro and the inverse FFTs multiply by the sin/cos values at the start of */ \
/* the macro.  We also implement the Daniel J. Bernstein (DJB) "exponent-1" idea to */ \
/* save sin/cos memory. */ \
/* */ \
/* */ \
\
/* In forward FFTs, */ \
/*    input values are R1+R5i, R2+R6i, R3+R7i, R4+R8i */ \
/*    output values are R1+R2i, R3+R4i, R5+R6i, R7+R8i */ \
\
/* In inverse FFTs, */ \
/*    input values are R1+R2i, R3+R4i, R5+R6i, R7+R8i */ \
/*    output values are R1+R5i, R2+R6i, R3+R7i, R4+R8i */ \
\
/* */ \
/* ************************************* four-complex-djbfft variants ****************************************** */ \
/* */ \
/* Macros to do Daniel J. Bernstein's exponent-1 butterflies.  The */ \
/* difference with a standard four-complex-fft is in the postmultiply step. */ \
/* Instead of multiplying by w^2x, w^x, w^3x, we multiply by w^2x, w^x, w^-x. */ \
/* */ \
\

#define r4_x4cl_four_complex_djbfft_preload \
	r4_x4cl_2sc_four_complex_djbfft_preload;

/* Special version used in pass 2 when first levels are radix-3. */

#define r4_x4cl_2sc_four_complex_djbfft_preload \
	r4_x4c_2sc_djbfft_mem_preload; \
;	/*r4_x4c_2sc_djbfft_partial_mem_preload -- assumed to be identical to r4_x4c_2sc_djbfft_mem_preload */ \


#ifdef LOW_CPU_REGISTER_COUNT
#define r4_x4cl_four_complex_djbfft(srcreg,srcinc,d1,d2,screg) { \
	vec2f64 *sc1a = (vec2f64*)(screg); \
	r4_x4cl_r4sc_four_complex_djbfft(srcreg,srcinc,d1,d2,sc0a[0],sc0a[1],sc0a[2],sc0a[3]); \
}
#define r4_x4cl_2sc_four_complex_djbfft(srcreg,srcinc,d1,d2,screg1,screg2) { \
	vec2f64 *sc1a = (vec2f64*)(screg1); \
	vec2f64 *sc2a = (vec2f64*)(screg2); \
	r4_x4cl_r4sc_four_complex_djbfft(srcreg,srcinc,d1,d2,sc1a[0],sc1a[1],sc2a[0],sc2a[1]); \
}
#else
#define r4_x4cl_four_complex_djbfft(srcreg,srcinc,d1,d2,screg) { \
	vec2f64 *sc1a = (vec2f64*)(screg); \
	vec2f64 sc0r = sc1a[0], sc1r = sc1a[1]; \
	vec2f64 sc2r = sc1a[2], sc3r = sc1a[3]; \
	r4_x4cl_r4sc_four_complex_djbfft(srcreg,srcinc,d1,d2,sc0r,sc1r,sc2r,sc3r); \
}
#define r4_x4cl_2sc_four_complex_djbfft(srcreg,srcinc,d1,d2,screg1,screg2) { \
	vec2f64 *sc1a = (vec2f64*)(screg1); \
	vec2f64 *sc2a = (vec2f64*)(screg2); \
	vec2f64 sc0r = sc1a[0], sc1r = sc1a[1]; \
	vec2f64 sc2r = sc2a[0], sc3r = sc2a[1]; \
	r4_x4cl_r4sc_four_complex_djbfft(srcreg,srcinc,d1,d2,sc0r,sc1r,sc2r,sc3r); \
}
#endif

#define r4_x4cl_r4sc_four_complex_djbfft(srcreg,srcinc,d1,d2,screg1,screg2,screg3,screg4) { \
	vec2f64 *sd0 = (vec2f64*)(srcreg); \
	vec2f64 *sd1 = (vec2f64*)(srcreg+d1); \
	vec2f64 *sd2 = (vec2f64*)(srcreg+d2); \
	vec2f64 *sd3 = (vec2f64*)(srcreg+d2+d1); \
	r4_x4c_r4sc_djbfft_mem(sd0[0],sd1[0],sd2[0],sd3[0], sd0[1],sd1[1],sd2[1],sd3[1], screg1,screg2,screg3,screg4,srcreg+srcinc,d1,sd0[0],sd0[1]); \
	/* Saved R1 */ \
	/* Saved I1 */ \
	xmm0 = sd0[2];	/* R1 */ \
	xmm5 = sd0[3];	/* R5 */ \
	sd0[2] = xmm7;	/* Save R2 */ \
	sd0[3] = xmm6;	/* Save I2 */ \
	xmm7 = sd1[2];	/* R2 */ \
	xmm6 = sd1[3];	/* R6 */ \
	sd1[0] = xmm3;	/* Save R3 */ \
	sd1[1] = xmm1;	/* Save I3 */ \
	sd1[2] = xmm2;	/* Save R4 */ \
	sd1[3] = xmm4;	/* Save I4 */ \
	r4_x4c_r4sc_djbfft_partial_mem(xmm0,xmm7,xmm3,xmm2, xmm5,xmm6,xmm1,xmm4, sd2[2],sd3[2],sd2[3],sd3[3], screg1,screg2,screg3,screg4,srcreg+srcinc+d2,d1, sd2[0],sd2[1]); \
	/* Saved R1 */ \
	/* Saved I1 */ \
	sd2[2] = xmm4;	/* Save R2 */ \
	sd2[3] = xmm1;	/* Save I2 */ \
	sd3[0] = xmm2;	/* Save R3 */ \
	sd3[1] = xmm7;	/* Save I3 */ \
	sd3[2] = xmm0;	/* Save R4 */ \
	sd3[3] = xmm5;	/* Save I4 */ \
	srcreg += srcinc; \
}

/* Used in last levels of an r4 FFT pass 1 (not in r4delay FFT).  No swizzling. */
#if 0 //IFDEF UNUSED

#define r4_g4cl_four_complex_djbfft_preload \
	r4_x4c_djbfft_mem_preload; \
;	/*r4_x4c_djbfft_partial_mem_preload -- assumed same as r4_x4c_djbfft_mem_preload */ \
\

#define r4_g4cl_four_complex_djbfft(srcreg,srcinc,d1,d2,dstreg,dstinc,e1,e2,screg) \
	xprefetch(u8ptr(srcreg+srcinc)); \
	r4_x4c_djbfft_mem(xptr(srcreg),xptr(srcreg+d1),xptr(srcreg+d2),xptr(srcreg+d2+d1),xptr(srcreg+16),xptr(srcreg+d1+16),xptr(srcreg+d2+16),xptr(srcreg+d2+d1+16),screg,0,dstreg+dstinc,e1,xptr(dstreg),xptr(dstreg+16)); \
	xprefetch(u8ptr(srcreg+srcinc+d1)); \
/*	xstore	[dstreg], xmm0		;; Save R1 */ \
/*	xstore	[dstreg+16], xmm0	;; Save I1 */ \
	xmm0 = xptr(srcreg+32);	/* R1 */ \
	xmm5 = xptr(srcreg+48);	/* R5 */ \
	xptr(dstreg+32) = xmm7;	/* Save R2 */ \
	xptr(dstreg+48) = xmm6;	/* Save I2 */ \
	xmm7 = xptr(srcreg+d1+32);	/* R2 */ \
	xmm6 = xptr(srcreg+d1+48);	/* R6 */ \
	xptr(dstreg+e1) = xmm3;	/* Save R3 */ \
	xptr(dstreg+e1+16) = xmm1;	/* Save I3 */ \
	xptr(dstreg+e1+32) = xmm2;	/* Save R4 */ \
	xptr(dstreg+e1+48) = xmm4;	/* Save I4 */ \
	xprefetch(u8ptr(srcreg+srcinc+d2)); \
	r4_x4c_djbfft_partial_mem(xmm0,xmm7,xmm3,xmm2,xmm5,xmm6,xmm1,xmm4,xptr(srcreg+d2+32),xptr(srcreg+d2+d1+32),xptr(srcreg+d2+48),xptr(srcreg+d2+d1+48),screg,0,dstreg+dstinc+e2,e1,xptr(dstreg+e2),xptr(dstreg+e2+16); \
	xprefetch(u8ptr(srcreg+srcinc+d2+d1))); \
	srcreg += srcinc; \
/*	xstore	[dstreg+e2], xmm0	;; Save R1 */ \
/*	xstore	[dstreg+e2+16], xmm0	;; Save I1 */ \
	xptr(dstreg+e2+32) = xmm4;	/* Save R2 */ \
	xptr(dstreg+e2+48) = xmm1;	/* Save I2 */ \
	xptr(dstreg+e2+e1) = xmm2;	/* Save R3 */ \
	xptr(dstreg+e2+e1+16) = xmm7;	/* Save I3 */ \
	xptr(dstreg+e2+e1+32) = xmm0;	/* Save R4 */ \
	xptr(dstreg+e2+e1+48) = xmm5;	/* Save I4 */ \
	dstreg += dstinc; \

#endif
\
/* Used in last levels of an r4 FFT pass 1 (not in r4delay FFT).  Swizzling. */ \
\

#define r4_sg4cl_four_complex_djbfft_preload \
	r4_x4c_djbfft_mem_preload; \
;	/*r4_x4c_djbfft_partial_mem_preload -- assumed to be identical to r4_x4c_djbfft_mem_preload */ \
\
\

#define r4_sg4cl_four_complex_djbfft(srcreg,srcinc,d1,d2,dstreg,dstinc,e1,e2,screg) { \
	vec2f64 *sd0 = (vec2f64*)(srcreg); \
	vec2f64 *sd1 = (vec2f64*)(srcreg+d1); \
	vec2f64 *sd2 = (vec2f64*)(srcreg+d2); \
	vec2f64 *sd3 = (vec2f64*)(srcreg+d2+d1); \
	vec2f64 *de0 = (vec2f64*)(dstreg); \
	vec2f64 *de1 = (vec2f64*)(dstreg+e1); \
	vec2f64 *de2 = (vec2f64*)(dstreg+e2); \
	vec2f64 *de3 = (vec2f64*)(dstreg+e2+e1); \
	xprefetch(u8ptr(srcreg+srcinc)); \
	r4_x4c_djbfft_mem(sd0[0],sd1[0],sd2[0],sd3[0],sd0[1],sd1[1],sd2[1],sd3[1], screg,0,dstreg+dstinc,e1, de0[0],de0[2]); \
	xprefetch(u8ptr(srcreg+srcinc+d1)); \
	shuffle_store_partial(de0[0], de0[1], xmm0, xmm7);				/* Save R1,R2 */ \
	shuffle_store_partial(de0[2], de0[3], xmm5, xmm6);			/* Save I1,I2 */ \
	xmm0 = sd0[2];	/* R1 */ \
	xmm7 = sd1[2];	/* R2 */ \
	xmm5 = sd0[3];	/* R5 */ \
	shuffle_store_with_temp(de1[0], de1[1], xmm3, xmm2, xmm6);		/* Save R3,R4 */ \
	shuffle_store_with_temp(de1[2], de1[3], xmm1, xmm4, xmm6);	/* Save I3,I4 */ \
	xmm6 = sd1[3];	/* R6 */ \
	xprefetch(u8ptr(srcreg+srcinc+d2)); \
	r4_x4c_djbfft_partial_mem(xmm0,xmm7,xmm3,xmm2,xmm5,xmm6,xmm1,xmm4, sd2[2],sd3[2],sd2[3],sd3[3], screg,0,dstreg+dstinc+e2,e1, de2[0],de2[2]); \
	xprefetch(u8ptr(srcreg+srcinc+d2+d1)); \
	srcreg += srcinc; \
	shuffle_store_partial(de2[0], de2[1], xmm3, xmm4);			/* Save R1,R2 */ \
	shuffle_store_partial(de2[2], de2[3], xmm6, xmm1);		/* Save I1,I2 */ \
	shuffle_store_with_temp(de3[0], de3[1], xmm2, xmm0, xmm3);	/* Save R3,R4 */ \
	shuffle_store_with_temp(de3[2], de3[3], xmm7, xmm5, xmm3);	/* Save I3,I4 */ \
	dstreg += dstinc; \
}

/* Used in first levels of pass 2 when last levels of pass 1 did not swizzle */
#if 0 //IFDEF UNUSED

#define r4_s2cl_four_complex_djbfft(srcreg,srcinc,d1,screg) \
	shuffle_load(xmm0, xmm2, xptr(srcreg+rbx), xptr(srcreg+32+rbx)); /* R1,R3 */ \
	shuffle_load(xmm1, xmm3, xptr(srcreg+d1+rbx), xptr(srcreg+d1+32+rbx)); /* R2,R4 */ \
	xprefetch(u8ptr(srcreg+srcinc+rbx)); \
	shuffle_load(xmm4, xmm6, xptr(srcreg+16+rbx), xptr(srcreg+48+rbx)); /* R5,R7 */ \
	shuffle_load(xmm5, xmm7, xptr(srcreg+d1+16+rbx), xptr(srcreg+d1+48+rbx)); /* R6,R8 */ \
	xprefetch(u8ptr(srcreg+srcinc+d1+rbx)); \
	r4_x4c_djbfft(xmm0,xmm1,xmm2,xmm3,xmm4,xmm5,xmm6,xmm7,screg,0,srcreg+srcinc,d1,xptr(srcreg),xptr(srcreg+16)); \
/*	xstore	[srcreg], xmm0		;; Save R1 */ \
/*	xstore	[srcreg+16], xmm0	;; Save I1 */ \
	xptr(srcreg+32) = xmm7;	/* Save R2 */ \
	xptr(srcreg+48) = xmm6;	/* Save I2 */ \
	xptr(srcreg+d1) = xmm3;	/* Save R3 */ \
	xptr(srcreg+d1+16) = xmm1;	/* Save I3 */ \
	xptr(srcreg+d1+32) = xmm0;	/* Save R4 */ \
	xptr(srcreg+d1+48) = xmm4;	/* Save I4 */ \
	srcreg += srcinc; \

#endif
\
/* Used in first levels of pass 2 when last levels of pass 1 swizzled. */ \
\

#define r4_f2cl_four_complex_djbfft_preload \
	r4_x4c_djbfft_mem_preload; \
\

#define r4_f2cl_four_complex_djbfft(srcreg,srcinc,d1,screg) { \
	uintptr_t src_rbx = srcreg+rbx;\
	vec2f64 *sd0 = (vec2f64*)(srcreg); \
	vec2f64 *sd1 = (vec2f64*)(srcreg+d1); \
	vec2f64 *bd0 = (vec2f64*)(src_rbx); \
	vec2f64 *bd1 = (vec2f64*)(src_rbx+d1); \
	r4_x4c_djbfft_mem(bd0[0],bd1[0],bd0[1],bd1[1],bd0[2],bd1[2],bd0[3],bd1[3], screg,0,src_rbx+srcinc,d1, sd0[0],sd0[1]); \
/*	xstore	[srcreg], xmm0		;; Save R1 */ \
/*	xstore	[srcreg+16], xmm0	;; Save I1 */ \
	sd0[2] = xmm7;	/* Save R2 */ \
	sd0[3] = xmm6;	/* Save I2 */ \
	sd1[0] = xmm3;	/* Save R3 */ \
	sd1[1] = xmm1;	/* Save I3 */ \
	sd1[2] = xmm2;	/* Save R4 */ \
	sd1[3] = xmm4;	/* Save I4 */ \
	srcreg += srcinc; \
}

/* Used in later levels of pass 2. */ \
\

#define r4_nf2cl_four_complex_djbfft_preload \
	r4_x4c_djbfft_mem_preload; \

#define r4_nf2cl_four_complex_djbfft(srcreg,srcinc,d1,screg) { \
	vec2f64 *sd0 = (vec2f64*)(srcreg); \
	vec2f64 *sd1 = (vec2f64*)(srcreg+d1); \
	r4_x4c_djbfft_mem(sd0[0],sd1[0],sd0[1],sd1[1],sd0[2],sd1[2],sd0[3],sd1[3], screg,0,srcreg+srcinc,d1, sd0[0],sd0[1]); \
/*	xstore	[srcreg], xmm0		;; Save R1 */ \
/*	xstore	[srcreg+16], xmm0	;; Save I1 */ \
	sd0[2] = xmm7;	/* Save R2 */ \
	sd0[3] = xmm6;	/* Save I2 */ \
	sd1[0] = xmm3;	/* Save R3 */ \
	sd1[1] = xmm1;	/* Save I3 */ \
	sd1[2] = xmm2;	/* Save R4 */ \
	sd1[3] = xmm4;	/* Save I4 */ \
	srcreg += srcinc; \
}

/* Used in pass 2 after radix-3 levels.  Uses two sin/cos ptrs instead of one. */ \


#define r4_x2cl_2sc_four_complex_djbfft_preload \
	r4_x4c_2sc_djbfft_mem_preload; \

#define r4_x2cl_2sc_four_complex_djbfft(srcreg,srcinc,d1,screg1,screg2) { \
	vec2f64 *sd0 = (vec2f64*)(srcreg); \
	vec2f64 *sd1 = (vec2f64*)(srcreg+d1); \
	r4_x4c_2sc_djbfft_mem(sd0[0],sd1[0],sd0[1],sd1[1],sd0[2],sd1[2],sd0[3],sd1[3], screg1,screg2,srcreg+srcinc,d1, sd0[0],sd0[1]); \
/*	xstore	[srcreg], xmm0		;; Save R1 */ \
/*	xstore	[srcreg+16], xmm0	;; Save I1 */ \
	sd0[2] = xmm7;	/* Save R2 */ \
	sd0[3] = xmm6;	/* Save I2 */ \
	sd1[0] = xmm3;	/* Save R3 */ \
	sd1[1] = xmm1;	/* Save I3 */ \
	sd1[2] = xmm2;	/* Save R4 */ \
	sd1[3] = xmm4;	/* Save I4 */ \
	srcreg += srcinc; \
}

/* Used in pass 2 after radix-5 levels. */ \
\

#define r4_x2cl_four_complex_djbfft_preload \
	r4_x4c_2sc_djbfft_mem_preload; \


#define r4_x2cl_four_complex_djbfft(srcreg,srcinc,d1,screg1) { \
	vec2f64 *sd0 = (vec2f64*)(srcreg); \
	vec2f64 *sd1 = (vec2f64*)(srcreg+d1); \
	r4_x4c_2sc_djbfft_mem(sd0[0],sd1[0],sd0[1],sd1[1],sd0[2],sd1[2],sd0[3],sd1[3], screg1,screg1+32,srcreg+srcinc,d1, sd0[0],sd0[1]); \
/*	xstore	[srcreg], xmm0		;; Save R1 */ \
/*	xstore	[srcreg+16], xmm0	;; Save I1 */ \
	sd0[2] = xmm7;	/* Save R2 */ \
	sd0[3] = xmm6;	/* Save I2 */ \
	sd1[0] = xmm3;	/* Save R3 */ \
	sd1[1] = xmm1;	/* Save I3 */ \
	sd1[2] = xmm2;	/* Save R4 */ \
	sd1[3] = xmm4;	/* Save I4 */ \
	srcreg += srcinc; \
}

/* The common macros to get the four complex djbfft job done */ \
\

#define r4_x4c_djbfft_mem_preload \
	r4_x4c_2sc_djbfft_mem_preload; \

#define r4_x4c_djbfft_mem(R1,R2,R3,R4,R5,R6,R7,R8,screg,off,pre1,pre2,dst1,dst2) { \
	vec2f64 *scv = (vec2f64*)(screg+off); \
	r4_x4c_r4sc_djbfft_mem(R1,R2,R3,R4,R5,R6,R7,R8,scv[0],scv[1],scv[2],scv[3],pre1,pre2,dst1,dst2); \
}


#define r4_x4c_2sc_djbfft_mem_preload \

#define r4_x4c_2sc_djbfft_mem(R1,R2,R3,R4,R5,R6,R7,R8, screg1,screg2, pre1,pre2,dst1,dst2) \
	r4_x4c_r4sc_djbfft_mem(R1,R2,R3,R4,R5,R6,R7,R8, xptr(screg1),xptr(screg1+16),xptr(screg2),xptr(screg2+16), pre1,pre2,dst1,dst2)

#define r4_x4c_r4sc_djbfft_mem(R1,R2,R3,R4,R5,R6,R7,R8,SC0,SC1,SC2,SC3,pre1,pre2,dst1,dst2) \
	xmm0 = R1;			/* R1 */ \
	xmm2 = R3;			/* R3 */ \
	xmm2 += xmm0;		/* R3 = R1 + R3 (new R1) */ \
\
	xmm1 = R2;			/* R2 */ \
	xmm3 = R4;			/* R4 */ \
	xmm3 += xmm1;		/* R4 = R2 + R4 (new R2) */ \
\
	xmm4 = R5;			/* I1 */ \
	xmm6 = R7;			/* I3 */ \
	xmm6 += xmm4;		/* I3 = I1 + I3 (new I1) */ \
\
	xmm0 -= R3;			/* R1 = R1 - R3 (new R3) */ \
	xmm1 -= R4;			/* R2 = R2 - R4 (new R4) */ \
\
	xprefetchw(u8ptr(pre1)); \
\
	xmm5 = xmm3;		/* Copy R2 */ \
	xmm3 += xmm2;		/* R2 = R1 + R2 (final R1) */ \
	xmm2 -= xmm5;		/* R1 = R1 - R2 (final R2) */ \
\
	xmm5 = R6;			/* I2 */ \
	xmm7 = R8;			/* I4 */ \
	xmm7 += xmm5;		/* I4 = I2 + I4 (new I2) */ \
\
	xmm4 -= R7;			/* I1 = I1 - I3 (new I3) */ \
	xmm5 -= R8;			/* I2 = I2 - I4 (new I4) */ \
\
	dst1 = xmm3;		/* Save R1 */ \
\
	xmm3 = xmm7;		/* Copy I2 */ \
	xmm7 += xmm6;		/* I2 = I1 + I2 (final I1) */ \
\
	xmm6 -= xmm3;		/* I1 = I1 - I2 (final I2) */ \
\
	dst2 = xmm7;		/* Save I1 */ \
\
	xmm3 = xmm1; \
	xmm1 += xmm4;		/* R4 = I3 + R4 (final I3) */ \
\
	xmm7 = xmm0;		/* Copy R3 */ \
	xmm0 -= xmm5;		/* R3 = R3 - I4 (final R3) */ \
\
	xmm4 -= xmm3;		/* I3 = I3 - R4 (final I4) */ \
	xmm5 += xmm7;		/* I4 = R3 + I4 (final R4) */ \
\
	xprefetchw(u8ptr(pre1+pre2)); \
\
	xmm7 = xmm2 * SC3 - xmm6;	/* A2 = R2 * cosine/sine - I2 */ \
	xmm3 = xmm0 * SC1 - xmm1;	/* A3 = R3 * cosine/sine - I3 */ \
	xmm6 = xmm6 * SC3 + xmm2;	/* B2 = I2 * cosine/sine + R2 */ \
	xmm2 = xmm5 * SC1 + xmm4;	/* A4 = R4 * cosine/sine + I4 */ \
	xmm1 = xmm1 * SC1 + xmm0;	/* B3 = I3 * cosine/sine + R3 */ \
	xmm4 = xmm4 * SC1 - xmm5;	/* B4 = I4 * cosine/sine - R4 */ \
\
	xmm7 *= SC2;		/* A2 = A2 * sine (final R2) */ \
	xmm6 *= SC2;		/* B2 = B2 * sine (final I2) */ \
\
	xmm3 *= SC0;		/* A3 = A3 * sine (final R3) */ \
	xmm1 *= SC0;		/* B3 = B3 * sine (final I3) */ \
\
	xmm2 *= SC0;		/* A4 = A4 * sine (final R4) */ \
	xmm4 *= SC0;		/* B4 = B4 * sine (final I4) */ \
\
\

#define r4_x4c_djbfft_partial_mem_preload \
	r4_x4c_2sc_djbfft_partial_mem_preload; \
\
\

#define r4_x4c_djbfft_partial_mem(r1,r2,r3,r4,r5,r6,r7,r8,mem3,mem4,mem7,mem8,screg,off,pre1,pre2,dst1,dst2) { \
	vec2f64 *scv = (vec2f64*)(screg+off); \
	r4_x4c_r4sc_djbfft_partial_mem(r1,r2,r3,r4,r5,r6,r7,r8,mem3,mem4,mem7,mem8,scv[0],scv[1],scv[2],scv[3],pre1,pre2,dst1,dst2); \
}


#define r4_x4c_2sc_djbfft_partial_mem_preload \
\
\

#define r4_x4c_2sc_djbfft_partial_mem(r1,r2,r3,r4,r5,r6,r7,r8,mem3,mem4,mem7,mem8,screg1,screg2,pre1,pre2,dst1,dst2) \
	r4_x4c_r4sc_djbfft_partial_mem(r1,r2,r3,r4,r5,r6,r7,r8,mem3,mem4,mem7,mem8,xptr(screg1),xptr(screg1+16),xptr(screg2),xptr(screg2+16),pre1,pre2,dst1,dst2)

#define r4_x4c_r4sc_djbfft_partial_mem(r1,r2,r3,r4,r5,r6,r7,r8,mem3,mem4,mem7,mem8,SC0,SC1,SC2,SC3,pre1,pre2,dst1,dst2) \
	r7 = r3 = mem3; \
	r3 += r1;			/* R3 = R1 + R3 (new R1) */ \
	r1 -= r7;			/* R1 = R1 - R3 (new R3) */ \
\
	r7 = r4 = mem4; \
	r4 += r2;			/* R4 = R2 + R4 (new R2) */ \
	r2 -= r7;			/* R2 = R2 - R4 (new R4) */ \
\
	r7 = r3; \
	r3 -= r4;			/* R1 = R1 - R2 (final R2) */ \
	r4 += r7;			/* R2 = R1 + R2 (final R1) */ \
\
	dst1 = r4;			/* Save R1 */ \
\
	r7 = mem7; \
	r8 = mem8; \
\
	r4 = r6; \
	r6 -= r8;			/* I2 = I2 - I4 (new I4) */ \
	r8 += r4;			/* I4 = I2 + I4 (new I2) */ \
\
	r4 = r5; \
	r5 -= r7;			/* I1 = I1 - I3 (new I3) */ \
	r7 += r4;			/* I3 = I1 + I3 (new I1) */ \
\
	r4 = r1; \
	r1 -= r6;			/* R3 = R3 - I4 (final R3) */ \
	r6 += r4;			/* I4 = R3 + I4 (final R4) */ \
\
	r4 = r2; \
	r2 += r5;			/* R4 = I3 + R4 (final I3) */ \
	r5 -= r4;			/* I3 = I3 - R4 (final I4) */ \
\
	xprefetchw(u8ptr(pre1)); \
\
	r4 = r7; \
	r7 -= r8;			/* I1 = I1 - I2 (final I2) */ \
	r8 += r4;			/* I2 = I1 + I2 (final I1) */ \
\
	dst2 = r8;			/* Save I1 */ \
\
	r4 = r1 * SC1 - r2;	/* A3 = R3 * cosine/sine - I3 */ \
	r2 = r2 * SC1 + r1;	/* B3 = I3 * cosine/sine + R3 */ \
	r8 = r3 * SC3 - r7;	/* A2 = R2 * cosine/sine - I2 */ \
	r7 = r7 * SC3 + r3;	/* B2 = I2 * cosine/sine + R2 */ \
	r1 = r6 * SC1 + r5;	/* A4 = R4 * cosine/sine + I4 */ \
	r5 = r5 * SC1 - r6;	/* B4 = I4 * cosine/sine - R4 */ \
\
	xprefetchw(u8ptr(pre1+pre2)); \
	r4 *= SC0;			/* A3 = A3 * sine (final R3) */ \
	r2 *= SC0;			/* B3 = B3 * sine (final I3) */ \
	r8 *= SC2;			/* A2 = A2 * sine (final R2) */ \
	r7 *= SC2;			/* B2 = B2 * sine (final I2) */ \
	r1 *= SC0;			/* A4 = A4 * sine (final R4) */ \
	r5 *= SC0;			/* B4 = B4 * sine (final I4) */ \
\

#if 0 //IFDEF UNUSED

#define r4_x4c_djbfft(r1,r2,r3,r4,r5,r6,r7,r8,screg,off,pre1,pre2,dst1,dst2) \
	r1 -= r3;			/* R1 = R1 - R3 (new R3) */ \
	r3 += r3; /* multwo */ \
	r3 += r1;			/* R3 = R1 + R3 (new R1) */ \
\
	r2 -= r4;			/* R2 = R2 - R4 (new R4) */ \
	r4 += r4; /* multwo */			/* R4 = R4 * 2 */ \
	r4 += r2;			/* R4 = R2 + R4 (new R2) */ \
\
	r3 -= r4;			/* R1 = R1 - R2 (final R2) */ \
	r4 += r4; /* multwo */			/* R2 = R2 * 2 */ \
	r4 += r3;			/* R2 = R1 + R2 (final R1) */ \
\
	dst1 = r4;		/* Save R1 */ \
\
	r4 = r6; \
	r6 -= r8;			/* I2 = I2 - I4 (new I4) */ \
	r8 += r4;			/* I4 = I2 + I4 (new I2) */ \
\
	r4 = r5; \
	r5 -= r7;			/* I1 = I1 - I3 (new I3) */ \
	r7 += r4;			/* I3 = I1 + I3 (new I1) */ \
\
	r4 = r1; \
	r1 -= r6;			/* R3 = R3 - I4 (final R3) */ \
	r6 += r4;			/* I4 = R3 + I4 (final R4) */ \
\
	r4 = r2; \
	r2 += r5;			/* R4 = I3 + R4 (final I3) */ \
	r5 -= r4;			/* I3 = I3 - R4 (final I4) */ \
\
	xprefetchw(u8ptr(pre1)); \
\
	r4 = r7; \
	r7 -= r8;			/* I1 = I1 - I2 (final I2) */ \
	r8 += r4;			/* I2 = I1 + I2 (final I1) */ \
\
	r4 = xptr(screg+off+0+16);	/* cosine/sine */ \
	r4 *= r1;			/* A3 = R3 * cosine/sine */ \
\
	dst2 = r8;		/* Save I1 */ \
\
	r8 = xptr(screg+off+32+16);	/* cosine/sine */ \
	r8 *= r3;			/* A2 = R2 * cosine/sine */ \
\
	r4 -= r2;			/* A3 = A3 - I3 */ \
	r2 *= xptr(screg+off+0+16);	/* B3 = I3 * cosine/sine */ \
	r2 += r1;			/* B3 = B3 + R3 */ \
\
	r1 = xptr(screg+off+0+16);	/* cosine/sine */ \
	r1 *= r6;			/* A4 = R4 * cosine/sine */ \
\
	r8 -= r7;			/* A2 = A2 - I2 */ \
	r7 *= xptr(screg+off+32+16);	/* B2 = I2 * cosine/sine */ \
\
	xprefetchw(u8ptr(pre1+pre2)); \
\
	r1 += r5;			/* A4 = A4 + I4 */ \
	r5 *= xptr(screg+off+0+16);	/* B4 = I4 * cosine/sine */ \
\
	r7 += r3;			/* B2 = B2 + R2 */ \
	r5 -= r6;			/* B4 = B4 - R4 */ \
\
	r3 = xptr(screg+off+0);	/* sine */ \
	r4 *= r3;			/* A3 = A3 * sine (final R3) */ \
	r2 *= r3;			/* B3 = B3 * sine (final I3) */ \
	r8 *= xptr(screg+off+32);	/* A2 = A2 * sine (final R2) */ \
	r1 *= r3;			/* A4 = A4 * sine (final R4) */ \
	r7 *= xptr(screg+off+32);	/* B2 = B2 * sine (final I2) */ \
	r5 *= r3;			/* B4 = B4 * sine (final I4) */ \

#endif
\
/* 32-bit AMD K8 optimized versions of the above macros */ \
\
/* These macros could be improved further if we didn't keep output register compatibility */ \
/* with the Core2/P4 macros.  That is, we'd have to rewrite all the callers of these */ \
/* macros like we do for the Core2/P4 64-bit macros.  One optimization revolves around */ \
/* noting that subpd/multwo/addpd can be changed to addpd/multwo/subpd giving us greater */ \
/* flexibility in ordering the computation of intermediate results. */ \

#if 0 // (@INSTR(,%xarch,<K8>) NE 0);
\
/* K8 Cheat sheet for scheduling dependency chains */ \
/*	      12345678901234567890123456789012345678901234567890123456789012345678901234567890 */ \
/*r13(r3)     AAAA */ \
/*i24(i4)       AAAA */ \
/*i13(i3)         AAAA */ \
/*r24(r4)           AAAA */ \
/*r13(r1)             AAAA */ \
/*r24(r2)               AAAA */ \
/*i13(i1)                 AAAA */ \
/*i24(i2)                   AAAA */ \
/*mR3(depR3-I4)               AAAA */ \
/*mI4(depI3-R4)	                AAAA			 */ \
/*mR2(depR1-R2)	                  AAAA */ \
/*mI2(depI1-I2) 	            AAAA */ \
/*mR4(depR3+I4)               MMMM    AAAA */ \
/*mI3(depI3+R4)	                MMMM    AAAA				 */ \
/*mR1(depR1+R2)	                  MMMM    AAAA */ \
/*mI1(depI1+I2)	                    MMMM    AAAA */ \
/*A3		                      MMMM      AAAAMMMM */ \
/*A2                                    MMMM  AAAAMMMM */ \
/*B2		                          MMMM    AAAAMMMM */ \
/*B3		                            MMMM      AAAAMMMM */ \
/*A4                                          MMMM  AAAAMMMM */ \
/*B4		                                MMMM    AAAAMMMM */ \
\

#define r4_x4c_2sc_djbfft_mem(R1,R2,R3,R4,R5,R6,R7,R8,screg1,screg2,pre1,pre2,dst1,dst2) \
	xmm0 = R1;		/* R1 */ \
	xmm2 = R3;		/* R3 */ \
	xmm0 -= xmm2;		/* R1 = R1 - R3 (new R3)		;; 1-4 */ \
\
	xmm5 = R6;		/* I2 */ \
	xmm7 = R8;		/* I4 */ \
	xmm5 -= xmm7;		/* I2 = I2 - I4 (new I4)		;; 3-6 */ \
\
	xmm4 = R5;		/* I1 */ \
	xmm6 = R7;		/* I3 */ \
	xmm4 -= xmm6;		/* I1 = I1 - I3 (new I3)		;; 5-8 */ \
\
	xmm1 = R2;		/* R2 */ \
	xmm3 = R4;		/* R4 */ \
	xmm1 -= xmm3;		/* R2 = R2 - R4 (new R4)		;; 7-10 */ \
\
	xmm2 += R1;		/* R3 = R1 + R3 (new R1)		;; 9-12 */ \
\
	xmm3 += R2;		/* R4 = R2 + R4 (new R2)		;; 11-14 */ \
\
	xmm6 += R5;		/* I3 = I1 + I3 (new I1)		;; 13-16 */ \
\
	xmm7 += R6;		/* I4 = I2 + I4 (new I2)		;; 15-18 */ \
\
	xprefetchw(u8ptr(pre1)); \
\
	xmm0 -= xmm5;		/* R3 = R3 - I4 (final R3)		;; 17-20 */ \
	xmm5 += xmm5; /* multwo */;			/* I4 = I4 * 2				;; 17-20 */ \
\
	xmm4 -= xmm1;		/* I3 = I3 - R4 (final I4)		;; 19-22 */ \
	xmm1 += xmm1; /* multwo */;			/* R4 = R4 * 2				;; 19-22 */ \
\
	xmm2 -= xmm3;		/* R1 = R1 - R2 (final R2)		;; 21-24 */ \
	xmm3 += xmm3; /* multwo */;			/* R2 = R2 * 2				;; 21-24 */ \
\
	xmm6 -= xmm7;		/* I1 = I1 - I2 (final I2)		;; 23-26 */ \
	xmm7 += xmm7; /* multwo */;			/* I2 = I2 * 2				;; 23-26 */ \
\
	xmm5 += xmm0;		/* I4 = R3 + I4 (final R4)		;; 25-28 */ \
	xmm1 += xmm4;		/* R4 = I3 + R4 (final I3)		;; 27-30 */ \
	xprefetchw(u8ptr(pre1+pre2)); \
	xmm3 += xmm2;		/* R2 = R1 + R2 (final R1)		;; 29-32 */ \
	xmm7 += xmm6;		/* I2 = I1 + I2 (final I1)		;; 31-34 */ \
\
	dst1 = xmm3;		/* Save R1				;; 33 */ \
\
	xmm3 = xptr(screg1+16);	/* cosine/sine */ \
	xmm3 *= xmm0;		/* A3 = R3 * cosine/sine		;; 25-28 */ \
\
	dst2 = xmm7;		/* Save I1				;; 35 */ \
\
	xmm7 = xptr(screg2+16);	/* cosine/sine */ \
	xmm7 *= xmm2;		/* A2 = R2 * cosine/sine		;; 27-30 */ \
\
	xmm7 -= xmm6;		/* A2 = A2 - I2				;; 33-36 */ \
	xmm6 *= xptr(screg2+16);	/* B2 = I2 * cosine/sine		;; 29-32 */ \
\
	xmm3 -= xmm1;		/* A3 = A3 - I3				;; 35-38 */ \
	xmm1 *= xptr(screg1+16);	/* B3 = I3 * cosine/sine		;; 31-34 */ \
\
	xmm6 += xmm2;		/* B2 = B2 + R2				;; 37-40 */ \
\
	xmm2 = xptr(screg1+16);	/* cosine/sine */ \
	xmm2 *= xmm5;		/* A4 = R4 * cosine/sine		;; 33-36 */ \
\
	xmm2 += xmm4;		/* A4 = A4 + I4				;; 39-42 */ \
	xmm4 *= xptr(screg1+16);	/* B4 = I4 * cosine/sine		;; 35-38 */ \
\
	xmm7 *= xptr(screg2);		/* A2 = A2 * sine (final R2)		;; 37-40 */ \
\
	xmm6 *= xptr(screg2);		/* B2 = B2 * sine (final I2)		;; 41-44 */ \
\
	xmm1 += xmm0;		/* B3 = B3 + R3				;; 41-44 */ \
	xmm0 = xptr(screg1);		/* Sine */ \
	xmm3 *= xmm0;		/* A3 = A3 * sine (final R3)		;; 39-42 */ \
\
	xmm4 -= xmm5;		/* B4 = B4 - R4				;; 43-46 */ \
	xmm2 *= xmm0;		/* A4 = A4 * sine (final R4)		;; 43-46 */ \
\
	xmm1 *= xmm0;		/* B3 = B3 * sine (final I3)		;; 45-48 */ \
	xmm4 *= xmm0;		/* B4 = B4 * sine (final I4)		;; 47-50 */ \
\
\
/* K8 Cheat sheet for scheduling dependency chains */ \
/*	      12345678901234567890123456789012345678901234567890123456789012345678901234567890 */ \
/*i13(i1)     AAAA */ \
/*i13(i3)       AAAA */ \
/*r24(r2)         AAAA */ \
/*r24(r4)           AAAA */ \
/*r13(r1)             AAAA */ \
/*r13(r3)               AAAA */ \
/*i24(i2)                 AAAA */ \
/*i24(i4)                   AAAA */ \
/*mI4(depI3-R4)	              AAAA			 */ \
/*mR3(depR3-I4)                 AAAA */ \
/*mR2(depR1-R2)	                  AAAA */ \
/*mI2(depI1-I2) 	            AAAA */ \
/*mI3(depI3+R4)	              MMMM    AAAA				 */ \
/*mR4(depR3+I4)                 MMMM    AAAA */ \
/*mR1(depR1+R2)	                  MMMM    AAAA */ \
/*mI1(depI1+I2)	                    MMMM    AAAA */ \
/*A3		                      MMMM    AAAAMMMM */ \
/*A2                                    MMMM      AAAAMMMM */ \
/*B3		                          MMMM  AAAAMMMM */ \
/*A4		                            MMMM    AAAAMMMM */ \
/*B2                                          MMMM    AAAAMMMM */ \
/*B4		                                MMMM    AAAAMMMM */ \
\

#define r4_x4c_2sc_djbfft_partial_mem(r1,r2,r3,r4,r5,r6,r7,r8,mem3,mem4,mem7,mem8,screg1,screg2,pre1,pre2,dst1,dst2) \
	r7 = mem7;						/* K8 */ \
	r7 += r5;			/* I3 = I1 + I3 (new I1)	;; 1-4 */ \
	r5 -= mem7;		/* I1 = I1 - I3 (new I3)	;; 3-6 */ \
\
	r4 = mem4; \
	r4 += r2;			/* R4 = R2 + R4 (new R2)	;; 5-8 */ \
	r2 -= mem4;		/* R2 = R2 - R4 (new R4)	;; 7-10 */ \
\
	r3 = mem3; \
	r3 += r1;			/* R3 = R1 + R3 (new R1)	;; 9-12 */ \
	r1 -= mem3;		/* R1 = R1 - R3 (new R3)	;; 11-14 */ \
\
	r8 = mem8; \
	r8 += r6;			/* I4 = I2 + I4 (new I2)	;; 13-16 */ \
	r6 -= mem8;		/* I2 = I2 - I4 (new I4)	;; 15-18 */ \
\
	xprefetchw(u8ptr(pre1)); \
\
	r5 -= r2;			/* I3 = I3 - R4 (final I4)	;; 17-20 */ \
	r2 += r2; /* multwo */			/* R4 = R4 * 2			;; 17-20 */ \
\
	r1 -= r6;			/* R3 = R3 - I4 (final R3)	;; 19-22 */ \
	r6 += r6; /* multwo */			/* I4 = I4 * 2			;; 19-22 */ \
\
	r3 -= r4;			/* R1 = R1 - R2 (final R2)	;; 21-24 */ \
	r4 += r4; /* multwo */			/* R2 = R2 * 2			;; 21-24 */ \
\
	r7 -= r8;			/* I1 = I1 - I2 (final I2)	;; 23-26 */ \
	r8 += r8; /* multwo */			/* I2 = I2 * 2			;; 23-26 */ \
\
	r6 += r1;			/* I4 = R3 + I4 (final R4)	;; 25-28 */ \
	r2 += r5;			/* R4 = I3 + R4 (final I3)	;; 27-30 */ \
	xprefetchw(u8ptr(pre1+pre2)); \
	r4 += r3;			/* R2 = R1 + R2 (final R1)	;; 29-32 */ \
	r8 += r7;			/* I2 = I1 + I2 (final I1)	;; 31-34 */ \
\
	dst1 = r4;		/* Save R1 */ \
\
	r4 = xptr(screg1+16);		/* cosine/sine */ \
	r4 *= r1;			/* A3 = R3 * cosine/sine	;; 25-28 */ \
\
	dst2 = r8;		/* Save I1 */ \
\
	r8 = xptr(screg2+16);		/* cosine/sine */ \
	r8 *= r3;			/* A2 = R2 * cosine/sine	;; 27-30 */ \
\
	r4 -= r2;			/* A3 = A3 - I3			;; 33-36 */ \
	r2 *= xptr(screg1+16);		/* B3 = I3 * cosine/sine	;; 29-32 */ \
\
	r2 += r1;			/* B3 = B3 + R3			;; 35-38 */ \
\
	r1 = xptr(screg1+16);		/* cosine/sine */ \
	r1 *= r6;			/* A4 = R4 * cosine/sine	;; 31-34 */ \
\
	r8 -= r7;			/* A2 = A2 - I2			;; 37-40 */ \
	r7 *= xptr(screg2+16);		/* B2 = I2 * cosine/sine	;; 33-36 */ \
\
	r1 += r5;			/* A4 = A4 + I4			;; 39-42 */ \
	r5 *= xptr(screg1+16);		/* B4 = I4 * cosine/sine	;; 35-38 */ \
\
	r7 += r3;			/* B2 = B2 + R2			;; 41-44 */ \
	r3 = xptr(screg1);		/* sine */ \
	r4 *= r3;			/* A3 = A3 * sine (final R3)	;; 37-40 */ \
\
	r2 *= r3;			/* B3 = B3 * sine (final I3)	;; 39-42 */ \
\
	r8 *= xptr(screg2);		/* A2 = A2 * sine (final R2)	;; 41-44 */ \
\
	r5 -= r6;			/* B4 = B4 - R4			;; 43-46 */ \
	r1 *= r3;			/* A4 = A4 * sine (final R4)	;; 43-46 */ \
\
	r7 *= xptr(screg2);		/* B2 = B2 * sine (final I2)	;; 45-48 */ \
\
	r5 *= r3;			/* B4 = B4 * sine (final I4)	;; 47-50 */ \
\
\

#endif
\
/* 64-bit Intel and AMD K10 implementations of the above - use the extra XMM registers */ \

#if 0 // (@INSTR(,%xarch,<CORE>) NE 0) OR (@INSTR(,%xarch,<P4>) NE 0) OR (@INSTR(,%xarch,<BLEND>) NE 0) OR (@INSTR(,%xarch,<K10>) NE 0);
#if 0 //DEF X86_64;
\
/* Theoretical best case is 44 clocks on a Core 2.  Now at 47.2 clocks. */ \
\

#define r4_x4cl_2sc_four_complex_djbfft_preload \
\
\

#define r4_x4cl_2sc_four_complex_djbfft(srcreg,srcinc,d1,d2,screg1,screg2) \
	xmm0 = xptr(srcreg);		/* R1 */ \
	xmm2 = xptr(srcreg+d2);	/* R3 */ \
	xmm5 = xptr(srcreg+d1+16);	/* I2 */ \
	xmm7 = xptr(srcreg+d2+d1+16);	/* I4 */ \
	xmm1 = xptr(srcreg+d1);	/* R2 */ \
	xmm3 = xptr(srcreg+d2+d1);	/* R4 */ \
	xmm4 = xptr(srcreg+16);	/* I1 */ \
	xmm6 = xptr(srcreg+d2+16);	/* I3 */ \
\
	xmm8 = xmm0;		/* Copy R1				; 1-3 */ \
	xmm0 -= xmm2;		/* R1 = R1 - R3 (new R3)		; 1-3 */ \
\
	xmm9 = xmm5;		/* Copy I2				; 2-4 */ \
	xmm5 -= xmm7;		/* I2 = I2 - I4 (new I4)		; 2-4 */ \
\
	xmm10 = xmm1;		/* Copy R2				; 3-5 */ \
	xmm1 -= xmm3;		/* R2 = R2 - R4 (new R4)		; 3-5 */ \
\
	xmm11 = xmm4;		/* Copy I1				; 4-6 */ \
	xmm4 -= xmm6;		/* I1 = I1 - I3 (new I3)		; 4-6	avail 12+ */ \
	xmm12 = xptr(screg1+16);	/* cosine/sine				; 4 */ \
\
	xmm2 += xmm8;		/* R3 = R1 + R3 (new R1)		; 5-7	avail 8,13+ */ \
	xmm13 = xmm12;							/* 5-7 */ \
\
	xmm3 += xmm10;		/* R4 = R2 + R4 (new R2)		; 6-8	avail 8,10,14+ */ \
	xmm14 = xmm12;							/* 6-8 */ \
\
	xmm7 += xmm9;		/* I4 = I2 + I4 (new I2)		; 7-9	avail 8.10.9,15 */ \
	xmm8 = xmm5;		/* Copy I4				; 7-9 */ \
\
	xmm6 += xmm11;		/* I3 = I1 + I3 (new I1)		; 8-10	avail 10.9,15,11 */ \
	xmm15 = xmm12;							/* 8-10 */ \
\
	xmm5 += xmm0;		/* I4 = R3 + I4 (final R4)		; 9-11	avail 10.9,11 */ \
	xmm9 = xmm4;		/* Copy I3				; 9-11 */ \
	xmm10 = xptr(screg2+16);	/* cosine/sine				; 9 */ \
\
	xmm0 -= xmm8;		/* R3 = R3 - I4 (final R3)		; 10-12	avail 11,8 */ \
\
	xprefetchw(u8ptr(srcreg+srcinc)); \
\
	xmm4 -= xmm1;		/* I3 = I3 - R4 (final I4)		; 11-13	avail 11,8 */ \
	xmm11 = xmm10;							/* 11-13 */ \
\
	xmm1 += xmm9;		/* R4 = I3 + R4 (final I3)		; 12-14	avail 8,9 */ \
	xmm12 *= xmm5;		/* A4 = R4 * cosine/sine		; 12-16 */ \
	xmm8 = xmm2;		/* Copy R1				; 12-14 */ \
\
	xmm2 -= xmm3;		/* R1 = R1 - R2 (final R2)		; 13-15	avail 9 */ \
	xmm13 *= xmm0;		/* A3 = R3 * cosine/sine		; 13-17 */ \
	xmm9 = xmm6;		/* Copy I1				; 13-15 */ \
\
	xmm6 -= xmm7;		/* I1 = I1 - I2 (final I2)		; 14-16	avail none */ \
	xmm14 *= xmm4;		/* B4 = I4 * cosine/sine		; 14-18 */ \
\
	xmm3 += xmm8;		/* R2 = R1 + R2 (final R1)		; 15-17	avail 8 */ \
	xmm15 *= xmm1;		/* B3 = I3 * cosine/sine		; 15-19 */ \
	xmm8 = xptr(screg1);		/* Sine					; 15 */ \
\
	xmm7 += xmm9;		/* I2 = I1 + I2 (final I1)		; 16-18	avail 9 */ \
	xmm10 *= xmm2;		/* A2 = R2 * cosine/sine		; 16-20 */ \
	xmm9 = xptr(screg2);							/* 16 */ \
\
	xmm12 += xmm4;		/* A4 = A4 + I4				; 17-19	avail 4 */ \
	xmm11 *= xmm6;		/* B2 = I2 * cosine/sine		; 17-21 */ \
	xmm4 = xptr(srcreg+32);	/*#2 R1					; 17 */ \
\
	xmm13 -= xmm1;		/* A3 = A3 - I3				; 18-20	avail 1 */ \
	xmm1 = xptr(srcreg+d2+32);	/* next R3 */ \
	xptr(srcreg) = xmm3;		/* Save R1				; 18	avail 3 */ \
\
	xprefetchw(u8ptr(srcreg+srcinc+d1)); \
\
	xmm14 -= xmm5;		/* B4 = B4 - R4				; 19-21	avail 3,5 */ \
	xmm3 = xptr(srcreg+d1+48);	/*#2 I2					; 19 */ \
	xptr(srcreg+16) = xmm7;	/* Save I1				; 19	avail 5,7 */ \
\
	xmm15 += xmm0;		/* B3 = B3 + R3				; 20-22	avail 5,7,0 */ \
	xmm12 *= xmm8;		/* A4 = A4 * sine (final R4)		; 20-24 */ \
	xmm5 = xptr(srcreg+d2+d1+48);	/*#2 I4 */ \
\
	xmm10 -= xmm6;		/* A2 = A2 - I2				; 21-23	avail 7,0,6 */ \
	xmm13 *= xmm8;		/* A3 = A3 * sine (final R3)		; 21-25 */ \
	xmm7 = xptr(srcreg+d1+32);	/*#2 R2 */ \
\
	xmm11 += xmm2;		/* B2 = B2 + R2				; 22-24	avail 0,6,2 */ \
	xmm14 *= xmm8;		/* B4 = B4 * sine (final I4)		; 22-26 */ \
	xmm0 = xptr(srcreg+d2+d1+32);	/*#2 R4 */ \
\
	xmm6 = xmm4;		/*#2 Copy R1				; 23-25	avail 2 */ \
	xmm4 -= xmm1;		/*#2 R1 = R1 - R3 (new R3)		; 23-25 */ \
	xmm15 *= xmm8;		/* B3 = B3 * sine (final I3)		; 23-27	avail 2,8 */ \
\
	xmm2 = xmm3;		/*#2 Copy I2				; 24-26	avail 8 */ \
	xmm3 -= xmm5;		/*#2 I2 = I2 - I4 (new I4)		; 24-26 */ \
	xmm10 *= xmm9;		/* A2 = A2 * sine (final R2)		; 24-28 */ \
\
	xmm8 = xmm7;		/*#2 Copy R2				; 25-27	avail none */ \
	xmm7 -= xmm0;		/*#2 R2 = R2 - R4 (new R4)		; 25-27 */ \
	xmm11 *= xmm9;		/* B2 = B2 * sine (final I2)		; 25-29	avail 9 */ \
	xmm9 = xptr(srcreg+48);	/*#2 I1					; 25 */ \
	xptr(srcreg+d1+32) = xmm12;	/* Save R4				; 25	avail 12 */ \
\
	xmm12 = xptr(srcreg+d2+48);	/*#2 I3					; 26	avail none */ \
	xptr(srcreg+d1) = xmm13;	/* Save R3				; 26	avail 13 */ \
	xmm13 = xmm9;		/*#2 Copy I1				; 26-28	avail none */ \
	xmm9 -= xmm12;		/*#2 I1 = I1 - I3 (new I3)		; 26-28 */ \
\
	xmm1 += xmm6;		/*#2 R3 = R1 + R3 (new R1)		; 27-29	avail 6 */ \
	xptr(srcreg+d1+48) = xmm14;	/* Save I4				; 27	avail 6,14 */ \
	xmm14 = xptr(screg1+16);	/*#2 cosine/sine			; 27	avail 6 */ \
\
	xmm0 += xmm8;		/*#2 R4 = R2 + R4 (new R2)		; 28-30	avail 6,8 */ \
	xptr(srcreg+d1+16) = xmm15;	/* Save I3				; 28	avail 6,8,15 */ \
	xmm15 = xmm14;							/* 28-30	avail 6,8 */ \
\
	xmm5 += xmm2;		/*#2 I4 = I2 + I4 (new I2)		; 29-31	avail 6,8,2 */ \
	xptr(srcreg+32) = xmm10;	/* Save R2				; 29	avail 6,8,2,10 */ \
	xmm6 = xmm3;		/*#2 Copy I4				; 29-31	avail 8,2,10 */ \
\
	xmm12 += xmm13;		/*#2 I3 = I1 + I3 (new I1)		; 30-32	avail 8,2,10,13 */ \
	xptr(srcreg+48) = xmm11;	/* Save I2				; 30	avail 8,2,10,13,11 */ \
	xmm10 = xmm14;							/* 30-32	avail 8,2,13,11 */ \
\
	xprefetchw(u8ptr(srcreg+srcinc+d2)); \
\
	xmm3 += xmm4;		/*#2 I4 = R3 + I4 (final R4)		; 31-33 */ \
	xmm11 = xmm14;							/* 31-33	avail 8,2,13 */ \
	xmm2 = xmm9;		/*#2 Copy I3				; 31-33	avail 8,13 */ \
\
	xmm4 -= xmm6;		/*#2 R3 = R3 - I4 (final R3)		; 32-34	avail 8,13,6 */ \
	xmm8 = xptr(screg2+16);	/*#2 cosine/sine			; 32	avail 13,6 */ \
\
	xmm9 -= xmm7;		/*#2 I3 = I3 - R4 (final I4)		; 33-35 */ \
	xmm13 = xmm8;							/* 33-35	avail 6 */ \
\
	xmm7 += xmm2;		/*#2 R4 = I3 + R4 (final I3)		; 34-36	avail 6,2 */ \
	xmm14 *= xmm3;		/*#2 A4 = R4 * cosine/sine		; 34-38 */ \
	xmm6 = xmm1;		/*#2 Copy R1				; 34-36	avail 2 */ \
\
	xmm1 -= xmm0;		/*#2 R1 = R1 - R2 (final R2)		; 35-37 */ \
	xmm15 *= xmm4;		/*#2 A3 = R3 * cosine/sine		; 35-39 */ \
	xmm2 = xmm12;		/*#2 Copy I1				; 35-37	avail none */ \
\
	xmm12 -= xmm5;		/*#2 I1 = I1 - I2 (final I2)		; 36-38 */ \
	xmm10 *= xmm9;		/*#2 B4 = I4 * cosine/sine		; 36-40 */ \
\
	xmm0 += xmm6;		/*#2 R2 = R1 + R2 (final R1)		; 37-39	avail 6 */ \
	xmm11 *= xmm7;		/*#2 B3 = I3 * cosine/sine		; 37-41 */ \
	xmm6 = xptr(screg1);		/*#2 Sine				; 37	avail none */ \
\
	xmm5 += xmm2;		/*#2 I2 = I1 + I2 (final I1)		; 38-40	avail 2 */ \
	xmm8 *= xmm1;		/*#2 A2 = R2 * cosine/sine		; 38-42 */ \
	xmm2 = xptr(screg2);		/*#2 Sine				; 38	avail none */ \
\
	xmm14 += xmm9;		/*#2 A4 = A4 + I4			; 39-41	avail 9 */ \
	xmm13 *= xmm12;		/*#2 B2 = I2 * cosine/sine		; 39-43 */ \
\
	xmm15 -= xmm7;		/*#2 A3 = A3 - I3			; 40-42	avail 9,7 */ \
	xptr(srcreg+d2) = xmm0;	/*#2 Save R1				; 40	avail 9,7,0 */ \
\
	xprefetchw(u8ptr(srcreg+srcinc+d2+d1)); \
\
	xmm10 -= xmm3;		/*#2 B4 = B4 - R4			; 41-43	avail 9,7,0,3 */ \
	xptr(srcreg+d2+16) = xmm5;	/*#2 Save I1				; 41	avail 9,7,0,3,5 */ \
\
	xmm11 += xmm4;		/*#2 B3 = B3 + R3			; 42-44	avail 9,7,0,3,5,4 */ \
	xmm14 *= xmm6;		/*#2 A4 = A4 * sine (final R4)		; 42-46 */ \
\
	xmm8 -= xmm12;		/*#2 A2 = A2 - I2			; 43-45	avail 9,7,0,3,5,4,12 */ \
	xmm15 *= xmm6;		/*#2 A3 = A3 * sine (final R3)		; 43-47 */ \
\
	xmm13 += xmm1;		/*#2 B2 = B2 + R2			; 44-46	avail 9,7,0,3,5,4,12,1 */ \
	xmm10 *= xmm6;		/*#2 B4 = B4 * sine (final I4)		; 44-48 */ \
\
	xmm11 *= xmm6;		/*#2 B3 = B3 * sine (final I3)		; 45-49 */ \
	xmm8 *= xmm2;		/*#2 A2 = A2 * sine (final R2)		; 46-50 */ \
	xmm13 *= xmm2;		/*#2 B2 = B2 * sine (final I2)		; 47-51 */ \
\
	xptr(srcreg+d2+d1+32) = xmm14; /*#2 Save R4				; 47 */ \
	xptr(srcreg+d2+d1) = xmm15;	/*#2 Save R3				; 48 */ \
	xptr(srcreg+d2+d1+48) = xmm10; /*#2 Save I4				; 49 */ \
	xptr(srcreg+d2+d1+16) = xmm11; /*#2 Save I3				; 50 */ \
	xptr(srcreg+d2+32) = xmm8;	/*#2 Save R2				; 51 */ \
	xptr(srcreg+d2+48) = xmm13;	/*#2 Save I2				; 52 */ \
\
	srcreg += srcinc; \
\
\

#define r4_f2cl_four_complex_djbfft_preload \
\
\

#define r4_f2cl_four_complex_djbfft(srcreg,srcinc,d1,screg) \
	r4_x4c_djbfft_mem64(xptr(srcreg+rbx),xptr(srcreg+d1+rbx),xptr(srcreg+16+rbx),xptr(srcreg+d1+16+rbx),xptr(srcreg+32+rbx),xptr(srcreg+d1+32+rbx),xptr(srcreg+48+rbx),xptr(srcreg+d1+48+rbx),screg,0,srcreg+srcinc+rbx,d1,xptr(srcreg),xptr(srcreg+16)); \
	xptr(srcreg+d1+32) = xmm12;	/* Save R4 */ \
	xptr(srcreg+d1) = xmm13;	/* Save R3 */ \
	xptr(srcreg+d1+48) = xmm14;	/* Save I4 */ \
	xptr(srcreg+d1+16) = xmm15;	/* Save I3 */ \
	xptr(srcreg+32) = xmm10;	/* Save R2 */ \
	xptr(srcreg+48) = xmm11;	/* Save I2 */ \
	srcreg += srcinc; \
\
\

#define r4_nf2cl_four_complex_djbfft_preload \
\
\

#define r4_nf2cl_four_complex_djbfft(srcreg,srcinc,d1,screg) \
	r4_x4c_djbfft_mem64(xptr(srcreg),xptr(srcreg+d1),xptr(srcreg+16),xptr(srcreg+d1+16),xptr(srcreg+32),xptr(srcreg+d1+32),xptr(srcreg+48),xptr(srcreg+d1+48),screg,0,srcreg+srcinc,d1,xptr(srcreg),xptr(srcreg+16)); \
	xptr(srcreg+d1+32) = xmm12;	/* Save R4 */ \
	xptr(srcreg+d1) = xmm13;	/* Save R3 */ \
	xptr(srcreg+d1+48) = xmm14;	/* Save I4 */ \
	xptr(srcreg+d1+16) = xmm15;	/* Save I3 */ \
	xptr(srcreg+32) = xmm10;	/* Save R2 */ \
	xptr(srcreg+48) = xmm11;	/* Save I2 */ \
	srcreg += srcinc; \
\
\

#define r4_x2cl_2sc_four_complex_djbfft_preload \
\
\

#define r4_x2cl_2sc_four_complex_djbfft(srcreg,srcinc,d1,screg1,screg2) \
	r4_x4c_2sc_djbfft_mem64(xptr(srcreg),xptr(srcreg+d1),xptr(srcreg+16),xptr(srcreg+d1+16),xptr(srcreg+32),xptr(srcreg+d1+32),xptr(srcreg+48),xptr(srcreg+d1+48),screg1,screg2,srcreg+srcinc,d1,xptr(srcreg),xptr(srcreg+16)); \
	xptr(srcreg+d1+32) = xmm12;	/* Save R4 */ \
	xptr(srcreg+d1) = xmm13;	/* Save R3 */ \
	xptr(srcreg+d1+48) = xmm14;	/* Save I4 */ \
	xptr(srcreg+d1+16) = xmm15;	/* Save I3 */ \
	xptr(srcreg+32) = xmm10;	/* Save R2 */ \
	xptr(srcreg+48) = xmm11;	/* Save I2 */ \
	srcreg += srcinc; \
\
\

#define r4_x2cl_four_complex_djbfft_preload \
\
\

#define r4_x2cl_four_complex_djbfft(srcreg,srcinc,d1,screg1) \
	r4_x4c_2sc_djbfft_mem64(xptr(srcreg),xptr(srcreg+d1),xptr(srcreg+16),xptr(srcreg+d1+16),xptr(srcreg+32),xptr(srcreg+d1+32),xptr(srcreg+48),xptr(srcreg+d1+48),screg1,screg1+32,srcreg+srcinc,d1,xptr(srcreg),xptr(srcreg+16)); \
	xptr(srcreg+d1+32) = xmm12;	/* Save R4 */ \
	xptr(srcreg+d1) = xmm13;	/* Save R3 */ \
	xptr(srcreg+d1+48) = xmm14;	/* Save I4 */ \
	xptr(srcreg+d1+16) = xmm15;	/* Save I3 */ \
	xptr(srcreg+32) = xmm10;	/* Save R2 */ \
	xptr(srcreg+48) = xmm11;	/* Save I2 */ \
	srcreg += srcinc; \
\
\
/* Core 2 optimal is 22 clocks, currently at 24.4 clocks. */ \
\

#define r4_x4c_djbfft_mem64(R1,R2,R3,R4,R5,R6,R7,R8,screg,off,pre1,pre2,dst1,dst2) \
	r4_x4c_2sc_djbfft_mem64(R1,R2,R3,R4,R5,R6,R7,R8,screg+off,screg+off+32,pre1,pre2,dst1,dst2); \
\

#define r4_x4c_2sc_djbfft_mem64(R1,R2,R3,R4,R5,R6,R7,R8,screg1,screg2,pre1,pre2,dst1,dst2) \
	xmm0 = R1;		/* R1 */ \
	xmm2 = R3;		/* R3 */ \
	xmm8 = xmm0;		/* Copy R1				; 1-3 */ \
	xmm0 -= xmm2;		/* R1 = R1 - R3 (new R3)		; 1-3 */ \
	xmm5 = R6;		/* I2 */ \
	xmm7 = R8;		/* I4 */ \
	xmm9 = xmm5;		/* Copy I2				; 2-4 */ \
	xmm5 -= xmm7;		/* I2 = I2 - I4 (new I4)		; 2-4 */ \
	xmm1 = R2;		/* R2 */ \
	xmm3 = R4;		/* R4 */ \
	xmm10 = xmm1;		/* Copy R2				; 3-5 */ \
	xmm1 -= xmm3;		/* R2 = R2 - R4 (new R4)		; 3-5 */ \
	xmm4 = R5;		/* I1 */ \
	xmm6 = R7;		/* I3 */ \
	xmm11 = xmm4;		/* Copy I1				; 4-6 */ \
	xmm4 -= xmm6;		/* I1 = I1 - I3 (new I3)		; 4-6 */ \
	xmm12 = xptr(screg1+16);	/* cosine/sine				; 4 */ \
\
	xmm2 += xmm8;		/* R3 = R1 + R3 (new R1)		; 5-7 */ \
	xmm13 = xmm12;							/* 5-7 */ \
	xmm3 += xmm10;		/* R4 = R2 + R4 (new R2)		; 6-8 */ \
	xmm14 = xmm12;							/* 6-8 */ \
	xmm7 += xmm9;		/* I4 = I2 + I4 (new I2)		; 7-9 */ \
	xmm8 = xmm5;		/* Copy I4				; 7-9 */ \
	xmm6 += xmm11;		/* I3 = I1 + I3 (new I1)		; 8-10 */ \
	xmm15 = xmm12;							/* 8-10 */ \
\
	xprefetchw(u8ptr(pre1)); \
\
	xmm5 += xmm0;		/* I4 = R3 + I4 (final R4)		; 9-11 */ \
	xmm9 = xmm4;		/* Copy I3				; 9-11 */ \
	xmm10 = xptr(screg2+16);	/* cosine/sine				; 9 */ \
\
	xmm0 -= xmm8;		/* R3 = R3 - I4 (final R3)		; 10-12 */ \
\
	xmm4 -= xmm1;		/* I3 = I3 - R4 (final I4)		; 11-13 */ \
	xmm11 = xmm10;							/* 11-13 */ \
\
	xmm1 += xmm9;		/* R4 = I3 + R4 (final I3)		; 12-14 */ \
	xmm12 *= xmm5;		/* A4 = R4 * cosine/sine		; 12-16 */ \
	xmm8 = xmm2;		/* Copy R1				; 12-14 */ \
\
	xmm2 -= xmm3;		/* R1 = R1 - R2 (final R2)		; 13-15 */ \
	xmm13 *= xmm0;		/* A3 = R3 * cosine/sine		; 13-17 */ \
	xmm9 = xmm6;		/* Copy I1				; 13-15 */ \
\
	xmm6 -= xmm7;		/* I1 = I1 - I2 (final I2)		; 14-16 */ \
	xmm14 *= xmm4;		/* B4 = I4 * cosine/sine		; 14-18 */ \
\
	xmm3 += xmm8;		/* R2 = R1 + R2 (final R1)		; 15-17 */ \
	xmm15 *= xmm1;		/* B3 = I3 * cosine/sine		; 15-19 */ \
\
	xmm7 += xmm9;		/* I2 = I1 + I2 (final I1)		; 16-18 */ \
	xmm10 *= xmm2;		/* A2 = R2 * cosine/sine		; 16-20 */ \
	xmm8 = xptr(screg1);		/* Sine					; 16 */ \
\
	xmm12 += xmm4;		/* A4 = A4 + I4				; 17-19 */ \
	xmm11 *= xmm6;		/* B2 = I2 * cosine/sine		; 17-21 */ \
	xmm9 = xptr(screg2);							/* 17 */ \
\
	xmm13 -= xmm1;		/* A3 = A3 - I3				; 18-20 */ \
	dst1 = xmm3;		/* Save R1				; 18 */ \
\
	xprefetchw(u8ptr(pre1+pre2)); \
\
	xmm14 -= xmm5;		/* B4 = B4 - R4				; 19-21 */ \
	dst2 = xmm7;		/* Save I1				; 19 */ \
\
	xmm15 += xmm0;		/* B3 = B3 + R3				; 20-22 */ \
	xmm12 *= xmm8;		/* A4 = A4 * sine (final R4)		; 20-24 */ \
\
	xmm10 -= xmm6;		/* A2 = A2 - I2				; 21-23 */ \
	xmm13 *= xmm8;		/* A3 = A3 * sine (final R3)		; 21-25 */ \
\
	xmm11 += xmm2;		/* B2 = B2 + R2				; 22-24 */ \
	xmm14 *= xmm8;		/* B4 = B4 * sine (final I4)		; 22-26 */ \
\
	xmm15 *= xmm8;		/* B3 = B3 * sine (final I3)		; 23-27 */ \
	xmm10 *= xmm9;		/* A2 = A2 * sine (final R2)		; 24-28 */ \
	xmm11 *= xmm9;		/* B2 = B2 * sine (final I2)		; 25-29 */ \
\

#endif
#endif
\
/* 64-bit AMD K8 optimized versions of the above macros.  Same as 32-bit macros with some */ \
/* constants preloaded.  We can probably do better than that. */ \

#if 0 // (@INSTR(,%xarch,<K8>) NE 0);
#if 0 //DEF X86_64;
\

#define r4_x4c_2sc_djbfft_mem_preload \
	xmm15 = XMM_TWO; \
\
\

#define r4_x4c_2sc_djbfft_mem(R1,R2,R3,R4,R5,R6,R7,R8,screg1,screg2,pre1,pre2,dst1,dst2) \
	xmm0 = R1;		/* R1 */ \
	xmm2 = R3;		/* R3 */ \
	xmm0 -= xmm2;		/* R1 = R1 - R3 (new R3)		;; 1-4 */ \
\
	xmm5 = R6;		/* I2 */ \
	xmm7 = R8;		/* I4 */ \
	xmm5 -= xmm7;		/* I2 = I2 - I4 (new I4)		;; 3-6 */ \
\
	xmm4 = R5;		/* I1 */ \
	xmm6 = R7;		/* I3 */ \
	xmm4 -= xmm6;		/* I1 = I1 - I3 (new I3)		;; 5-8 */ \
\
	xmm1 = R2;		/* R2 */ \
	xmm3 = R4;		/* R4 */ \
	xmm1 -= xmm3;		/* R2 = R2 - R4 (new R4)		;; 7-10 */ \
\
	xmm2 += R1;		/* R3 = R1 + R3 (new R1)		;; 9-12 */ \
\
	xmm3 += R2;		/* R4 = R2 + R4 (new R2)		;; 11-14 */ \
\
	xmm6 += R5;		/* I3 = I1 + I3 (new I1)		;; 13-16 */ \
\
	xmm7 += R6;		/* I4 = I2 + I4 (new I2)		;; 15-18 */ \
\
	xprefetchw(u8ptr(pre1)); \
\
	xmm0 -= xmm5;		/* R3 = R3 - I4 (final R3)		;; 17-20 */ \
	xmm5 *= xmm15;		/* I4 = I4 * 2				;; 17-20 */ \
\
	xmm4 -= xmm1;		/* I3 = I3 - R4 (final I4)		;; 19-22 */ \
	xmm1 *= xmm15;		/* R4 = R4 * 2				;; 19-22 */ \
\
	xmm2 -= xmm3;		/* R1 = R1 - R2 (final R2)		;; 21-24 */ \
	xmm3 *= xmm15;		/* R2 = R2 * 2				;; 21-24 */ \
\
	xmm6 -= xmm7;		/* I1 = I1 - I2 (final I2)		;; 23-26 */ \
	xmm7 *= xmm15;		/* I2 = I2 * 2				;; 23-26 */ \
\
	xmm5 += xmm0;		/* I4 = R3 + I4 (final R4)		;; 25-28 */ \
	xmm1 += xmm4;		/* R4 = I3 + R4 (final I3)		;; 27-30 */ \
	xprefetchw(u8ptr(pre1+pre2)); \
	xmm3 += xmm2;		/* R2 = R1 + R2 (final R1)		;; 29-32 */ \
	xmm7 += xmm6;		/* I2 = I1 + I2 (final I1)		;; 31-34 */ \
\
	dst1 = xmm3;		/* Save R1				;; 33 */ \
\
	xmm3 = xptr(screg1+16);	/* cosine/sine */ \
	xmm3 *= xmm0;		/* A3 = R3 * cosine/sine		;; 25-28 */ \
\
	dst2 = xmm7;		/* Save I1				;; 35 */ \
\
	xmm7 = xptr(screg2+16);	/* cosine/sine */ \
	xmm7 *= xmm2;		/* A2 = R2 * cosine/sine		;; 27-30 */ \
\
	xmm7 -= xmm6;		/* A2 = A2 - I2				;; 33-36 */ \
	xmm6 *= xptr(screg2+16);	/* B2 = I2 * cosine/sine		;; 29-32 */ \
\
	xmm3 -= xmm1;		/* A3 = A3 - I3				;; 35-38 */ \
	xmm1 *= xptr(screg1+16);	/* B3 = I3 * cosine/sine		;; 31-34 */ \
\
	xmm6 += xmm2;		/* B2 = B2 + R2				;; 37-40 */ \
\
	xmm2 = xptr(screg1+16);	/* cosine/sine */ \
	xmm2 *= xmm5;		/* A4 = R4 * cosine/sine		;; 33-36 */ \
\
	xmm2 += xmm4;		/* A4 = A4 + I4				;; 39-42 */ \
	xmm4 *= xptr(screg1+16);	/* B4 = I4 * cosine/sine		;; 35-38 */ \
\
	xmm7 *= xptr(screg2);		/* A2 = A2 * sine (final R2)		;; 37-40 */ \
\
	xmm6 *= xptr(screg2);		/* B2 = B2 * sine (final I2)		;; 41-44 */ \
\
	xmm1 += xmm0;		/* B3 = B3 + R3				;; 41-44 */ \
	xmm0 = xptr(screg1);		/* Sine */ \
	xmm3 *= xmm0;		/* A3 = A3 * sine (final R3)		;; 39-42 */ \
\
	xmm4 -= xmm5;		/* B4 = B4 - R4				;; 43-46 */ \
	xmm2 *= xmm0;		/* A4 = A4 * sine (final R4)		;; 43-46 */ \
\
	xmm1 *= xmm0;		/* B3 = B3 * sine (final I3)		;; 45-48 */ \
	xmm4 *= xmm0;		/* B4 = B4 * sine (final I4)		;; 47-50 */ \
\
\

#define r4_x4c_2sc_djbfft_partial_mem_preload \
	xmm15 = XMM_TWO; \
\
\

#define r4_x4c_2sc_djbfft_partial_mem(r1,r2,r3,r4,r5,r6,r7,r8,mem3,mem4,mem7,mem8,screg1,screg2,pre1,pre2,dst1,dst2) \
	r7 = mem7;						/* K8 */ \
	r7 += r5;			/* I3 = I1 + I3 (new I1)	;; 1-4 */ \
	r5 -= mem7;		/* I1 = I1 - I3 (new I3)	;; 3-6 */ \
\
	r4 = mem4; \
	r4 += r2;			/* R4 = R2 + R4 (new R2)	;; 5-8 */ \
	r2 -= mem4;		/* R2 = R2 - R4 (new R4)	;; 7-10 */ \
\
	r3 = mem3; \
	r3 += r1;			/* R3 = R1 + R3 (new R1)	;; 9-12 */ \
	r1 -= mem3;		/* R1 = R1 - R3 (new R3)	;; 11-14 */ \
\
	r8 = mem8; \
	r8 += r6;			/* I4 = I2 + I4 (new I2)	;; 13-16 */ \
	r6 -= mem8;		/* I2 = I2 - I4 (new I4)	;; 15-18 */ \
\
	xprefetchw(u8ptr(pre1)); \
\
	r5 -= r2;			/* I3 = I3 - R4 (final I4)	;; 17-20 */ \
	r2 *= xmm15;		/* R4 = R4 * 2			;; 17-20 */ \
\
	r1 -= r6;			/* R3 = R3 - I4 (final R3)	;; 19-22 */ \
	r6 *= xmm15;		/* I4 = I4 * 2			;; 19-22 */ \
\
	r3 -= r4;			/* R1 = R1 - R2 (final R2)	;; 21-24 */ \
	r4 *= xmm15;		/* R2 = R2 * 2			;; 21-24 */ \
\
	r7 -= r8;			/* I1 = I1 - I2 (final I2)	;; 23-26 */ \
	r8 *= xmm15;		/* I2 = I2 * 2			;; 23-26 */ \
\
	r6 += r1;			/* I4 = R3 + I4 (final R4)	;; 25-28 */ \
	r2 += r5;			/* R4 = I3 + R4 (final I3)	;; 27-30 */ \
	xprefetchw(u8ptr(pre1+pre2)); \
	r4 += r3;			/* R2 = R1 + R2 (final R1)	;; 29-32 */ \
	r8 += r7;			/* I2 = I1 + I2 (final I1)	;; 31-34 */ \
\
	dst1 = r4;		/* Save R1 */ \
\
	r4 = xptr(screg1+16);		/* cosine/sine */ \
	r4 *= r1;			/* A3 = R3 * cosine/sine	;; 25-28 */ \
\
	dst2 = r8;		/* Save I1 */ \
\
	r8 = xptr(screg2+16);		/* cosine/sine */ \
	r8 *= r3;			/* A2 = R2 * cosine/sine	;; 27-30 */ \
\
	r4 -= r2;			/* A3 = A3 - I3			;; 33-36 */ \
	r2 *= xptr(screg1+16);		/* B3 = I3 * cosine/sine	;; 29-32 */ \
\
	r2 += r1;			/* B3 = B3 + R3			;; 35-38 */ \
\
	r1 = xptr(screg1+16);		/* cosine/sine */ \
	r1 *= r6;			/* A4 = R4 * cosine/sine	;; 31-34 */ \
\
	r8 -= r7;			/* A2 = A2 - I2			;; 37-40 */ \
	r7 *= xptr(screg2+16);		/* B2 = I2 * cosine/sine	;; 33-36 */ \
\
	r1 += r5;			/* A4 = A4 + I4			;; 39-42 */ \
	r5 *= xptr(screg1+16);		/* B4 = I4 * cosine/sine	;; 35-38 */ \
\
	r7 += r3;			/* B2 = B2 + R2			;; 41-44 */ \
	r3 = xptr(screg1);		/* sine */ \
	r4 *= r3;			/* A3 = A3 * sine (final R3)	;; 37-40 */ \
\
	r2 *= r3;			/* B3 = B3 * sine (final I3)	;; 39-42 */ \
\
	r8 *= xptr(screg2);		/* A2 = A2 * sine (final R2)	;; 41-44 */ \
\
	r5 -= r6;			/* B4 = B4 - R4			;; 43-46 */ \
	r1 *= r3;			/* A4 = A4 * sine (final R4)	;; 43-46 */ \
\
	r7 *= xptr(screg2);		/* B2 = B2 * sine (final I2)	;; 45-48 */ \
\
	r5 *= r3;			/* B4 = B4 * sine (final I4)	;; 47-50 */ \
\
\

#endif
#endif
\
/* */ \
/* ************************************* four-complex-djbunfft variants ****************************************** */ \
/* */ \
\

#define r4_x4cl_four_complex_djbunfft_preload \
	r4_x4c_djbunfft_mem_preload; \
;	/*r4_x4c_djbunfft_partial_mem -- assume same as r4_x4c_djbunfft_mem_preload */ \
\
\

#define r4_x4cl_four_complex_djbunfft(srcreg,srcinc,d1,d2,screg_a,scoff) { \
	vec2f64 *sd0 = (vec2f64*)(srcreg); \
	vec2f64 *sd1 = (vec2f64*)(srcreg+d1); \
	vec2f64 *sd2 = (vec2f64*)(srcreg+d2); \
	vec2f64 *sd3 = (vec2f64*)(srcreg+d2+d1); \
	char *screg = (char*)(screg_a); /* improves clang compiler */ \
\
	r4_x4c_djbunfft_mem(sd0[1],sd0[3],sd1[1],sd1[3],sd2[1],sd2[3],sd3[1],sd3[3], sd1[1],screg,scoff,srcreg+srcinc+d2,d1); \
/*	xstore	[srcreg+d1+16], xmm6	;; Save R3 */ \
	vec2f64 xma0,xma1,xma2,xma3,xma4,xma5,xma6,xma7; \
	xma2 = sd1[0];	/* R3 */ \
	sd1[0] = xmm7;	/* Save R1 */ \
	xma3 = sd1[2];	/* R4 */ \
	sd1[2] = xmm4;	/* Save I1 */ \
	sd1[3] = xmm3;	/* Save I3 */ \
	xma6 = sd3[0];	/* R7 */ \
	xma7 = sd3[2];	/* R8 */ \
	sd3[0] = xmm5;	/* Save R2 */ \
	sd3[1] = xmm0;	/* Save R4 */ \
	sd3[2] = xmm2;	/* Save I2 */ \
	sd3[3] = xmm1;	/* Save I4 */ \
	r4_x4c_djbunfft_partial_mem(xma0,xma1,xma2,xma3,xma4,xma5,xma6,xma7, sd0[0],sd0[2],sd2[0],sd2[2],sd0[0], screg,0,srcreg+srcinc,d1); \
/*	xstore	[srcreg], xma1		;; Save R1 */ \
	sd0[1] = xma2;	/* Save R3 */ \
	sd0[2] = xma4;	/* Save I1 */ \
	sd0[3] = xma0;	/* Save I3 */ \
	sd2[0] = xma3;	/* Save R2 */ \
	sd2[1] = xma5;	/* Save R4 */ \
	sd2[2] = xma6;	/* Save I2 */ \
	sd2[3] = xma7;	/* Save I4 */ \
	srcreg += srcinc; \
}

/* Used in pass 2 when there are radix-3 levels.  Uses 2 sin/cos ptrs */ \


#define r4_x4cl_2sc_four_complex_djbunfft_preload \
	r4_x4c_2sc_djbunfft_mem_preload; \
;	/*r4_x4c_2sc_djbunfft_partial_mem_preload -- assume same as r4_x4c_2sc_djbunfft_mem_preload */ \
\
\

#define r4_x4cl_2sc_four_complex_djbunfft(srcreg,srcinc,d1,d2,screg1,scoff1,screg2,scoff2) { \
	vec2f64 *sd0 = (vec2f64*)(srcreg); \
	vec2f64 *sd1 = (vec2f64*)(srcreg+d1); \
	vec2f64 *sd2 = (vec2f64*)(srcreg+d2); \
	vec2f64 *sd3 = (vec2f64*)(srcreg+d2+d1); \
	r4_x4c_2sc_djbunfft_mem(sd0[1],sd0[3],sd1[1],sd1[3],sd2[1],sd2[3],sd3[1],sd3[3], sd1[1],screg1+scoff1,screg2+scoff2,srcreg+srcinc+d2,d1); \
/*	xstore	[srcreg+d1+16], xmm6	;; Save R3 */ \
	vec2f64 xma0,xma1,xma2,xma3,xma4,xma5,xma6,xma7; \
	xma2 = sd1[0];	/* R3 */ \
	sd1[0] = xmm7;	/* Save R1 */ \
	xma3 = sd1[2];	/* R4 */ \
	sd1[2] = xmm4;	/* Save I1 */ \
	sd1[3] = xmm3;	/* Save I3 */ \
	xma6 = sd3[0];	/* R7 */ \
	xma7 = sd3[2];	/* R8 */ \
	sd3[0] = xmm5;	/* Save R2 */ \
	sd3[1] = xmm0;	/* Save R4 */ \
	sd3[2] = xmm2;	/* Save I2 */ \
	sd3[3] = xmm1;	/* Save I4 */ \
	r4_x4c_2sc_djbunfft_partial_mem(xma0,xma1,xma2,xma3,xma4,xma5,xma6,xma7, sd0[0],sd0[2],sd2[0],sd2[2],sd0[0], screg1,screg2,srcreg+srcinc,d1); \
/*	xstore	[srcreg], xma1		;; Save R1 */ \
	sd0[1] = xma2;	/* Save R3 */ \
	sd0[2] = xma4;	/* Save I1 */ \
	sd0[3] = xma0;	/* Save I3 */ \
	sd2[0] = xma3;	/* Save R2 */ \
	sd2[1] = xma5;	/* Save R4 */ \
	sd2[2] = xma6;	/* Save I2 */ \
	sd2[3] = xma7;	/* Save I4 */ \
	srcreg += srcinc; \
}

/* Used in last levels of an r4 FFT pass 1 (not in r4delay FFT).  No swizzling. */
#if 0 //IFDEF UNUSED

#define r4_g2cl_four_complex_djbunfft_preload \
	r4_x4c_djbunfft_mem_preload; \
\

#define r4_g2cl_four_complex_djbunfft(srcreg,srcinc,d1,dstreg,dstinc,e1,screg) \
	r4_x4c_djbunfft_mem(xptr(srcreg),xptr(srcreg+16),xptr(srcreg+32),xptr(srcreg+48),xptr(srcreg+d1),xptr(srcreg+d1+16),xptr(srcreg+d1+32),xptr(srcreg+d1+48),xptr(dstreg+16),screg,0,dstreg+dstinc,e1); \
	srcreg += srcinc; \
	xptr(dstreg) = xmm7;		/* Save R1 */ \
/*	xstore	[dstreg+16], xmm6	;; Save R3 */ \
	xptr(dstreg+32) = xmm4;	/* Save I1 */ \
	xptr(dstreg+48) = xmm3;	/* Save I3 */ \
	xptr(dstreg+e1) = xmm5;	/* Save R2 */ \
	xptr(dstreg+e1+16) = xmm0;	/* Save R4 */ \
	xptr(dstreg+e1+32) = xmm2;	/* Save I2 */ \
	xptr(dstreg+e1+48) = xmm1;	/* Save I4 */ \
	dstreg += dstinc; \
\
\
/* Used in first levels of pass 2.  Swizzling. */ \

#define r4_s4cl_four_complex_djbunfft(srcreg,srcinc,d1,d2,screg,scdist) \
	r4_x4c_djbunfft_mem_shuffle(xptr(srcreg+16),xptr(srcreg+48),xptr(srcreg+d1+16),xptr(srcreg+d1+48),xptr(srcreg+d2+16),xptr(srcreg+d2+48),xptr(srcreg+d2+d1+16),xptr(srcreg+d2+d1+48),xptr(srcreg+d1+16),xptr(srcreg+d1+48),screg,scdist,srcreg+srcinc+d2,d1); \
/*	shuffle_store [srcreg+d1+16], [srcreg+d1+48], xmm2, xmm6 ;; Save I1, I3 */ \
	xmm2 = xptr(srcreg+d1);	/* R3 */ \
	xmm6 = xptr(srcreg+d1+32);	/* R4 */ \
	shuffle_store(xptr(srcreg+d1), xptr(srcreg+d1+32), xmm5, xmm3); /* Save R1, R3 */ \
	xmm5 = xptr(srcreg+d2+d1);	/* R7 */ \
	xmm3 = xptr(srcreg+d2+d1+32);	/* R8 */ \
	shuffle_store(xptr(srcreg+d2+d1), xptr(srcreg+d2+d1+32), xmm7, xmm0); /* Save R2, R4 */ \
	shuffle_store(xptr(srcreg+d2+d1+16), xptr(srcreg+d2+d1+48), xmm4, xmm1); /* Save I2, I4 */ \
	r4_x4c_djbunfft_partial_mem_shuffle(xmm7,xmm0,xmm2,xmm6,xmm4,xmm1,xmm5,xmm3,xptr(srcreg),xptr(srcreg+32),xptr(srcreg+d2),xptr(srcreg+d2+32),xptr(srcreg),xptr(srcreg+32),screg,0,srcreg+srcinc,d1); \
/*	shuffle_store [srcreg], [srcreg+32], xmm0, xmm2 ;; Save R1, R3 */ \
	shuffle_store(xptr(srcreg+16), xptr(srcreg+48), xmm4, xmm7); /* Save I1, I3 */ \
	shuffle_store(xptr(srcreg+d2), xptr(srcreg+d2+32), xmm6, xmm1); /* Save R2, R4 */ \
	shuffle_store(xptr(srcreg+d2+16), xptr(srcreg+d2+48), xmm5, xmm3); /* Save I2, I4 */ \
	srcreg += srcinc; \

#endif
\
/* Used in last levels of an r4 FFT pass 1 (not in r4delay FFT).  Swizzling. */ \

#define r4_sg2cl_four_complex_djbunfft(srcreg,srcinc,d1,dstreg,dstinc,e1,screg) { \
	vec2f64 *sd0 = (vec2f64*)(srcreg); \
	vec2f64 *sd1 = (vec2f64*)(srcreg+d1); \
	vec2f64 *de0 = (vec2f64*)(dstreg); \
	vec2f64 *de1 = (vec2f64*)(dstreg+e1); \
	vec2f64 *sc = (vec2f64*)(screg); \
\
	vec2f64 temp1; \
	shuffle_load_with_temp(xmm7, xmm2, sd1[0], sd1[1], xmm1); /* R3,R4 */ \
	xmm1 = sc[1];		/* cosine/sine */ \
	xmm3 = xmm7;		/* Copy R3 */ \
	xmm0 = xmm2;		/* Copy R4 */ \
\
	shuffle_load_with_temp(xmm5, xmm4, sd1[2], sd1[3], xmm6); /* I3,I4 */ \
	xmm7 = xmm7 * xmm1 + xmm5;	/* A3 = R3 * cosine/sine + I3 */ \
	xmm2 = xmm2 * xmm1 - xmm4;	/* A4 = R4 * cosine/sine - I4 */ \
	xmm5 = xmm5 * xmm1 - xmm3;	/* B3 = I3 * cosine/sine - R3 */ \
	xmm4 = xmm4 * xmm1 + xmm0;	/* B4 = I4 * cosine/sine + R4 */ \
\
	shuffle_load_with_temp(xmm0, xmm1, sd0[0], sd0[1], xmm6);	/* R1,R2 */ \
\
	xmm3 = sc[0];		/* Sine */ \
	xmm7 *= xmm3;		/* A3 = A3 * sine (new R3) */ \
	xmm5 *= xmm3;		/* B3 = B3 * sine (new I3) */ \
	xmm2 *= xmm3;		/* A4 = A4 * sine (new R4) */ \
	xmm4 *= xmm3;		/* B4 = B4 * sine (new I4) */ \
\
	temp1 = xmm5;		/* Temporarily save I3 */ \
\
	shuffle_load(xmm5, xmm3, sd0[2], sd0[3]); /* I1,I2 */ \
	xmm6 = xmm1 * sc[3] + xmm3;	/* A2 = R2 * cosine/sine + I2 */ \
	xmm3 = xmm3 * sc[3] - xmm1;	/* B2 = I2 * cosine/sine - R2 */ \
	xmm1 = sc[2]; \
	xmm6 *= xmm1;		/* A2 = A2 * sine (new R2) */ \
	xmm3 *= xmm1;		/* B2 = B2 * sine (new I2) */ \
\
	xprefetchw(u8ptr(dstreg+dstinc)); \
	srcreg += srcinc; \
\
	xmm1 = xmm0;		/* Copy R1 */ \
	xmm0 -= xmm6;		/* R1 = R1 - R2 (new R2) */ \
	xmm6 += xmm1;		/* R2 = R1 + R2 (new R1) */ \
\
	xmm1 = xmm2;		/* Copy R4 */ \
	xmm2 -= xmm7;		/* R4 = R4 - R3 (new I4) */ \
	xmm7 += xmm1;		/* R3 = R4 + R3 (new R3) */ \
\
	xmm1 = xmm5;		/* Copy I1 */ \
	xmm5 -= xmm3;		/* I1 = I1 - I2 (new I2) */ \
	xmm3 += xmm1;		/* I2 = I1 + I2 (new I1) */ \
\
	xmm1 = xmm7;		/* Copy R3 */ \
	xmm7 += xmm6;		/* R3 = R1 + R3 (final R1) */ \
	xmm6 -= xmm1;		/* R1 = R1 - R3 (final R3) */ \
\
	xmm1 = temp1;		/* Reload I3 */ \
	de0[0] = xmm7;		/* Save R1 */ \
\
	xmm7 = xmm1;		/* Copy I3 */ \
	xmm1 -= xmm4;		/* I3 = I3 - I4 (new R4) */ \
	xmm4 += xmm7;		/* I4 = I3 + I4 (new I3) */ \
\
	xprefetchw(u8ptr(dstreg+dstinc+e1)); \
\
	xmm7 = xmm0;		/* Copy R2 */ \
	xmm0 -= xmm1;		/* R2 = R2 - R4 (final R4) */ \
	xmm1 += xmm7;		/* R4 = R2 + R4 (final R2) */ \
\
	xmm7 = xmm5;		/* Copy I2 */ \
	xmm5 -= xmm2;		/* I2 = I2 - I4 (final I4) */ \
	xmm2 += xmm7;		/* I4 = I2 + I4 (final I2) */ \
\
	xmm7 = xmm3;		/* Copy I1 */ \
	xmm3 -= xmm4;		/* I1 = I1 - I3 (final I3) */ \
	xmm4 += xmm7;		/* I3 = I1 + I3 (final I1) */ \
\
	de0[1] = xmm6;	/* Save R3 */ \
	de0[2] = xmm4;	/* Save I1 */ \
	de0[3] = xmm3;	/* Save I3 */ \
	de1[0] = xmm1;	/* Save R2 */ \
	de1[1] = xmm0;	/* Save R4 */ \
	de1[2] = xmm2;	/* Save I2 */ \
	de1[3] = xmm5;	/* Save I4 */ \
	dstreg += dstinc; \
}


#define r4_x4c_djbunfft_mem_preload \
	r4_x4c_2sc_djbunfft_mem_preload; \
\
\

#define r4_x4c_djbunfft_mem(mem1,mem2,mem3,mem4,mem5,mem6,mem7,mem8,dst3,screg,off,pre1,pre2) \
	r4_x4c_2sc_djbunfft_mem(mem1,mem2,mem3,mem4,mem5,mem6,mem7,mem8,dst3,screg+off,screg+off+32,pre1,pre2); \
\
\

#define r4_x4c_2sc_djbunfft_mem_preload \
\
\

#define r4_x4c_2sc_djbunfft_mem_PART1(mem1,mem2,mem3,mem4,mem5,mem6,mem7,mem8,screg1,screg2,pre1,pre2) { \
	vec2f64 *sc1 = (vec2f64*)(screg1); \
	vec2f64 *sc2 = (vec2f64*)(screg2); \
	xmm3 = sc2[1];		/* cosine/sine */ \
	xmm0 = mem4;		/* I2 */ \
	xmm1 = mem3;		/* I2 */ \
	xmm6 = xmm1 * xmm3 + xmm0;	/* A2 = R2 * cosine/sine + I2 */ \
	xmm3 = xmm3 * xmm0 - xmm1;	/* B2 = I2 * cosine/sine - R2 */ \
\
	xmm4 = sc1[1];		/* cosine/sine */ \
	xmm2 = mem5;		/* R3 */ \
	xmm5 = mem6;		/* I3 */ \
	xmm7 = xmm2 * xmm4 + xmm5;		/* A3 = R3 * cosine/sine + I3 */ \
	xmm5 = xmm5 * xmm4 - xmm2;		/* B3 = I3 * cosine/sine - R3 */ \
\
	xmm1 = mem7;		/* R4 */ \
	xmm0 = mem8;		/* I4 */ \
	xmm2 = xmm1 * xmm4 - xmm0;		/* A4 = R4 * cosine/sine - I4 */ \
	xmm4 = xmm4 * xmm0 + xmm1;		/* B4 = I4 * cosine/sine + R4 */ \
\
	xmm1 = sc2[0];		/* Sine */ \
	xmm6 *= xmm1;		/* A2 = A2 * sine (new R2) */ \
	xmm3 *= xmm1;		/* B2 = B2 * sine (new I2) */ \
	xmm1 = sc1[0];		/* Sine */ \
	xmm7 *= xmm1;		/* A3 = A3 * sine (new R3) */ \
	xmm5 *= xmm1;		/* B3 = B3 * sine (new I3) */ \
	xmm2 *= xmm1;		/* A4 = A4 * sine (new R4) */ \
	xmm4 *= xmm1;		/* B4 = B4 * sine (new I4) */ \
\
	xprefetchw(u8ptr(pre1)); \
\
	xmm0 = mem1;		/* R1 */ \
	xmm0 -= xmm6;		/* R1 = R1 - R2 (new R2) */ \
	xmm6 += mem1;		/* R2 = R1 + R2 (new R1) */ \
\
	xmm1 = xmm2;		/* Copy R4 */ \
	xmm2 -= xmm7;		/* R4 = R4 - R3 (new I4) */ \
	xmm7 += xmm1;		/* R3 = R4 + R3 (new R3) */ \
\
	xmm1 = xmm5;		/* Copy I3 */ \
	xmm5 -= xmm4;		/* I3 = I3 - I4 (new R4) */ \
	xmm4 += xmm1;		/* I4 = I3 + I4 (new I3) */ \
\
	xmm1 = xmm6;		/* Copy R1 */ \
	xmm6 -= xmm7;		/* R1 = R1 - R3 (final R3) */ \
	xmm7 += xmm1;		/* R3 = R1 + R3 (final R1) */ \
\
	xmm1 = mem2;		/* I1 */ \
	xmm1 -= xmm3;		/* I1 = I1 - I2 (new I2) */ \
	xmm3 += mem2;		/* I2 = I1 + I2 (new I1) */ \
	xprefetchw(u8ptr(pre1+pre2)); \
}

#if 0 // unused
#define r4_x4c_2sc_djbunfft_mem_nodst3(mem1,mem2,mem3,mem4,mem5,mem6,mem7,mem8,screg1,screg2,pre1,pre2) \
	r4_x4c_2sc_djbunfft_mem_PART1(mem1,mem2,mem3,mem4,mem5,mem6,mem7,mem8,screg1,screg2,pre1,pre2); \
	xmm0 -= xmm5;		/* R2 = R2 - R4 (final R4) */ \
	xmm5 += xmm5; /* multwo */;			/* R4 = R4 * 2 */ \
	xmm5 += xmm0;		/* R4 = R2 + R4 (final R2) */ \
\
	xmm1 -= xmm2;		/* I2 = I2 - I4 (final I4) */ \
	xmm2 += xmm2; /* multwo */;			/* I4 = I4 * 2 */ \
	xmm2 += xmm1;		/* I4 = I2 + I4 (final I2) */ \
\
	xmm3 -= xmm4;		/* I1 = I1 - I3 (final I3) */ \
	xmm4 += xmm4; /* multwo */;			/* I3 = I3 * 2 */ \
	xmm4 += xmm3;		/* I3 = I1 + I3 (final I1) */
#endif

#define r4_x4c_2sc_djbunfft_mem(mem1,mem2,mem3,mem4,mem5,mem6,mem7,mem8,dst3,screg1,screg2,pre1,pre2) \
	r4_x4c_2sc_djbunfft_mem_PART1(mem1,mem2,mem3,mem4,mem5,mem6,mem7,mem8,screg1,screg2,pre1,pre2); \
	dst3 = xmm6; \
\
	xmm6 = xmm0;		/* Copy R2 */ \
	xmm0 -= xmm5;		/* R2 = R2 - R4 (final R4) */ \
	xmm5 += xmm6;		/* R4 = R2 + R4 (final R2) */ \
\
	xmm6 = xmm1;		/* Copy I2 */ \
	xmm1 -= xmm2;		/* I2 = I2 - I4 (final I4) */ \
	xmm2 += xmm6;		/* I4 = I2 + I4 (final I2) */ \
\
	xmm6 = xmm3;		/* Copy I1 */ \
	xmm3 -= xmm4;		/* I1 = I1 - I3 (final I3) */ \
	xmm4 += xmm6;		/* I3 = I1 + I3 (final I1) */ \


#define r4_x4c_djbunfft_partial_mem_preload \
	r4_x4c_2sc_djbunfft_partial_mem_preload; \
\
\

#define r4_x4c_djbunfft_partial_mem(r1,r2,r3,r4,r5,r6,r7,r8,mem1,mem2,mem5,mem6,dst1,screg,off,pre1,pre2) \
	r4_x4c_2sc_djbunfft_partial_mem(r1,r2,r3,r4,r5,r6,r7,r8,mem1,mem2,mem5,mem6,dst1,screg+off,screg+off+32,pre1,pre2); \
\
\

#define r4_x4c_2sc_djbunfft_partial_mem_preload \
\
\

#define r4_x4c_2sc_djbunfft_partial_mem_PART1(r1,r2,r3,r4,r5,r6,r7,r8,mem1,mem2,mem5,mem6,screg1,screg2,pre1,pre2) { \
	vec2f64 *sc1 = (vec2f64*)(screg1); \
	vec2f64 *sc2 = (vec2f64*)(screg2); \
\
	r1 = sc2[1];		/* cosine/sine */ \
	r2 = r3;			/* Copy R2 */ \
	r3 = r3 * r1 + r4;	/* A2 = R2 * cosine/sine + I2 */ \
	r1 = r1 * r4 - r2;	/* B2 = I2 * cosine/sine - R2 */ \
\
	r4 = sc1[1];		/* cosine/sine */ \
	r6 = r7;			/* Copy R4 */ \
	r7 = r7 * r4 - r8;	/* A4 = R4 * cosine/sine - I4 */ \
	r5 = r4 * r8 + r6;	/* B4 = I4 * cosine/sine + R4 */ \
\
	r6 = mem5;			/* R3 */ \
	r8 = mem6;			/* I3 */ \
	r2 = r6 * r4 + r8;	/* A3 = R3 * cosine/sine + I3 */ \
	r4 = r4 * r8 - r6;	/* B3 = I3 * cosine/sine - R3 */ \
\
	r8 = sc2[0];		/* Sine */ \
	r3 *= r8;			/* A2 = A2 * sine (new R2) */ \
	r1 *= r8;			/* B2 = B2 * sine (new I2) */ \
	r8 = sc1[0];		/* Sine */ \
	r7 *= r8;			/* A4 = A4 * sine (new R4) */ \
	r5 *= r8;			/* B4 = B4 * sine (new I4) */ \
	r2 *= r8;			/* A3 = A3 * sine (new R3) */ \
	r4 *= r8;			/* B3 = B3 * sine (new I3) */ \
\
	xprefetchw(u8ptr(pre1)); \
\
	r6 = mem1;			/* R1 */ \
	r6 -= r3;			/* R1 = R1 - R2 (new R2) */ \
	r3 += mem1;			/* R2 = R1 + R2 (new R1) */ \
\
	r8 = r7;			/* Copy R4 */ \
	r7 -= r2;			/* R4 = R4 - R3 (new I4) */ \
	r2 += r8;			/* R3 = R4 + R3 (new R3) */ \
\
	r8 = r4;			/* Copy I3 */ \
	r4 -= r5;			/* I3 = I3 - I4 (new R4) */ \
	r5 += r8;			/* I4 = I3 + I4 (new I3) */ \
\
	r8 = r3;			/* Copy R1 */ \
	r3 -= r2;			/* R1 = R1 - R3 (final R3) */ \
	r2 += r8;			/* R3 = R1 + R3 (final R1) */ \
\
	r8 = mem2;			/* I1 */ \
	r8 -= r1;			/* I1 = I1 - I2 (new I2) */ \
	r1 += mem2;			/* I2 = I1 + I2 (new I1) */ \
\
	xprefetchw(u8ptr(pre1+pre2)); \
}
#if 0 // unused
#define r4_x4c_2sc_djbunfft_partial_mem_nodst1(r1,r2,r3,r4,r5,r6,r7,r8,mem1,mem2,mem5,mem6,screg1,screg2,pre1,pre2) \
	r4_x4c_2sc_djbunfft_partial_mem_PART1(r1,r2,r3,r4,r5,r6,r7,r8,mem1,mem2,mem5,mem6,screg1,screg2,pre1,pre2) \
	r6 -= r4;			/* R2 = R2 - R4 (final R4) */ \
	r4 += r4; /* multwo */			/* R4 = R4 * 2 */ \
	r4 += r6;			/* R4 = R2 + R4 (final R2) */ \
\
	r8 -= r7;			/* I2 = I2 - I4 (final I4) */ \
	r7 += r7; /* multwo */			/* I4 = I4 * 2 */ \
	r7 += r8;			/* I4 = I2 + I4 (final I2) */ \
\
	r1 -= r5;			/* I1 = I1 - I3 (final I3) */ \
	r5 += r5; /* multwo */			/* I3 = I3 * 2 */ \
	r5 += r1;			/* I3 = I1 + I3 (final I1) */
#endif

#define r4_x4c_2sc_djbunfft_partial_mem(r1,r2,r3,r4,r5,r6,r7,r8, mem1,mem2,mem5,mem6, dst1,screg1,screg2,pre1,pre2) \
	r4_x4c_2sc_djbunfft_partial_mem_PART1( r1,r2,r3,r4,r5,r6,r7,r8, mem1,mem2,mem5,mem6, screg1,screg2,pre1,pre2) \
	dst1 = r2; \
\
	r2 = r6;			/* Copy R2 */ \
	r6 -= r4;			/* R2 = R2 - R4 (final R4) */ \
	r4 += r2;			/* R4 = R2 + R4 (final R2) */ \
\
	r2 = r8;			/* Copy I2 */ \
	r8 -= r7;			/* I2 = I2 - I4 (final I4) */ \
	r7 += r2;			/* I4 = I2 + I4 (final I2) */ \
\
	r2 = r1;			/* Copy I1 */ \
	r1 -= r5;			/* I1 = I1 - I3 (final I3) */ \
	r5 += r2;			/* I3 = I1 + I3 (final I1) */ \


#if 0 //IFDEF UNUSED

#define r4_x4c_djbunfft_mem_shuffle(mem1,mem2,mem3,mem4,mem5,mem6,mem7,mem8,dst5,dst7,screg,off,pre1,pre2) \
	xmm3 = xptr(screg+off+32+16);	/* cosine/sine */ \
	xmm6 = mem4;		/* I2 */ \
	xmm6 *= xmm3;		/* B2 = I2 * cosine/sine */ \
	xmm0 = mem3;		/* R2 */ \
	xmm3 *= xmm0;		/* A2 = R2 * cosine/sine */ \
\
	xmm5 = xptr(screg+off+0+16);	/* cosine/sine */ \
	xmm7 = mem6;		/* I3 */ \
	xmm7 *= xmm5;		/* B3 = I3 * cosine/sine */ \
	xmm1 = mem5;		/* R3 */ \
	xmm4 = xmm5;		/* cosine/sine */ \
	xmm5 *= xmm1;		/* A3 = R3 * cosine/sine */ \
\
	xmm2 = mem8;		/* I4 */ \
	xmm2 *= xmm4;		/* B4 = I4 * cosine/sine */ \
	xmm6 -= xmm0;		/* B2 = B2 - R2 */ \
	xmm0 = mem7;		/* R4 */ \
	xmm4 *= xmm0;		/* A4 = R4 * cosine/sine */ \
\
	xmm7 -= xmm1;		/* B3 = B3 - R3 */ \
	xmm2 += xmm0;		/* B4 = B4 + R4 */ \
	xmm3 += mem4;		/* A2 = A2 + I2 */ \
	xmm5 += mem6;		/* A3 = A3 + I3 */ \
	xmm4 -= mem8;		/* A4 = A4 - I4 */ \
\
	xmm6 *= xptr(screg+off+32);	/* B2 = B2 * sine (new I2) */ \
	xmm1 = xptr(screg+off+0);	/* Sine */ \
	xmm7 *= xmm1;		/* B3 = B3 * sine (new I3) */ \
	xmm2 *= xmm1;		/* B4 = B4 * sine (new I4) */ \
	xmm3 *= xptr(screg+off+32);	/* A2 = A2 * sine (new R2) */ \
	xmm5 *= xmm1;		/* A3 = A3 * sine (new R3) */ \
	xmm4 *= xmm1;		/* A4 = A4 * sine (new R4) */ \
\
	xprefetchw(u8ptr(pre1)); \
\
	xmm1 = mem2;		/* I1 */ \
	xmm1 -= xmm6;		/* I1 = I1 - I2 (new I2) */ \
	xmm6 += mem2;		/* I2 = I1 + I2 (new I1) */ \
\
	xmm0 = xmm7;		/* Copy I3 */ \
	xmm7 -= xmm2;		/* I3 = I3 - I4 (new R4) */ \
	xmm2 += xmm0;		/* I4 = I3 + I4 (new I3) */ \
\
	xmm0 = xmm4;		/* Copy R4 */ \
	xmm4 -= xmm5;		/* R4 = R4 - R3 (new I4) */ \
	xmm5 += xmm0;		/* R3 = R4 + R3 (new R3) */ \
\
	xmm0 = xmm6;		/* Copy I1 */ \
	xmm6 -= xmm2;		/* I1 = I1 - I3 (final I3) */ \
	xmm2 += xmm0;		/* I3 = I1 + I3 (final I1) */ \
\
	xmm0 = mem1;		/* R1 */ \
	xmm0 -= xmm3;		/* R1 = R1 - R2 (new R2) */ \
	xmm3 += mem1;		/* R2 = R1 + R2 (new R1) */ \
\
	xprefetchw(u8ptr(pre1+pre2)); \
\
	shuffle_store(dst5, dst7, xmm2, xmm6); /* Save I1, I3 */ \
\
	xmm6 = xmm1;		/* Copy I2 */ \
	xmm1 -= xmm4;		/* I2 = I2 - I4 (final I4) */ \
	xmm4 += xmm6;		/* I4 = I2 + I4 (final I2) */ \
\
	xmm2 = xmm3;		/* Copy R1 */ \
	xmm3 -= xmm5;		/* R1 = R1 - R3 (final R3) */ \
	xmm5 += xmm2;		/* R3 = R1 + R3 (final R1) */ \
\
	xmm6 = xmm0;		/* Copy R2 */ \
	xmm0 -= xmm7;		/* R2 = R2 - R4 (final R4) */ \
	xmm7 += xmm6;		/* R4 = R2 + R4 (final R2) */ \
\
\

#define r4_x4c_djbunfft_partial_mem_shuffle(r1,r2,r3,r4,r5,r6,r7,r8,mem1,mem2,mem5,mem6,dst1,dst3,screg,off,pre1,pre2) \
	r1 = xptr(screg+off+32+16);	/* cosine/sine */ \
	r2 = r3;			/* Copy R2 */ \
	r3 *= r1;			/* A2 = R2 * cosine/sine */ \
	r1 *= r4;			/* B2 = I2 * cosine/sine */ \
\
	r5 = xptr(screg+off+0+16);	/* cosine/sine */ \
	r6 = r7;			/* Copy R4 */ \
	r7 *= r5;			/* A4 = R4 * cosine/sine */ \
	r3 += r4;			/* A2 = A2 + I2 */ \
	r4 = r5;			/* cosine/sine */ \
	r5 *= r8;			/* B4 = I4 * cosine/sine */ \
	r1 -= r2;			/* B2 = B2 - R2 */ \
\
	r2 = mem5;		/* R3 */ \
	r2 *= r4;			/* A3 = R3 * cosine/sine */ \
	r4 *= mem6;		/* B3 = I3 * cosine/sine */ \
\
	r7 -= r8;			/* A4 = A4 - I4 */ \
	r5 += r6;			/* B4 = B4 + R4 */ \
	r2 += mem6;		/* A3 = A3 + I3 */ \
	r4 -= mem5;		/* B3 = B3 - R3 */ \
\
	r3 *= xptr(screg+off+32);	/* A2 = A2 * sine (new R2) */ \
	r1 *= xptr(screg+off+32);	/* B2 = B2 * sine (new I2) */ \
	r8 = xptr(screg+off+0);	/* Sine */ \
	r7 *= r8;			/* A4 = A4 * sine (new R4) */ \
	r5 *= r8;			/* B4 = B4 * sine (new I4) */ \
	r2 *= r8;			/* A3 = A3 * sine (new R3) */ \
	r4 *= r8;			/* B3 = B3 * sine (new I3) */ \
\
	xprefetchw(u8ptr(pre1)); \
\
	r6 = mem1;		/* R1 */ \
	r6 -= r3;			/* R1 = R1 - R2 (new R2) */ \
	r3 += mem1;		/* R2 = R1 + R2 (new R1) */ \
\
	r8 = r7;			/* Copy R4 */ \
	r7 -= r2;			/* R4 = R4 - R3 (new I4) */ \
	r2 += r8;			/* R3 = R4 + R3 (new R3) */ \
\
	r8 = r4;			/* Copy I3 */ \
	r4 -= r5;			/* I3 = I3 - I4 (new R4) */ \
	r5 += r8;			/* I4 = I3 + I4 (new I3) */ \
\
	r8 = r3;			/* Copy R1 */ \
	r3 -= r2;			/* R1 = R1 - R3 (final R3) */ \
	r2 += r8;			/* R3 = R1 + R3 (final R1) */ \
\
	r8 = mem2;		/* I1 */ \
	r8 -= r1;			/* I1 = I1 - I2 (new I2) */ \
	r1 += mem2;		/* I2 = I1 + I2 (new I1) */ \
\
	xprefetchw(u8ptr(pre1+pre2)); \
\
	shuffle_store(dst1, dst3, r2, r3);	 /* Save R1 and R3 */ \
\
	r2 = r6;			/* Copy R2 */ \
	r6 -= r4;			/* R2 = R2 - R4 (final R4) */ \
	r4 += r2;			/* R4 = R2 + R4 (final R2) */ \
\
	r3 = r8;			/* Copy I2 */ \
	r8 -= r7;			/* I2 = I2 - I4 (final I4) */ \
	r7 += r3;			/* I4 = I2 + I4 (final I2) */ \
\
	r2 = r1;			/* Copy I1 */ \
	r1 -= r5;			/* I1 = I1 - I3 (final I3) */ \
	r5 += r2;			/* I3 = I1 + I3 (final I1) */ \

#endif

#if 0 //IFDEF UNUSED

#define r4_x4c_djbunfft(r1,r2,r3,r4,r5,r6,r7,r8,mem7,mem8,dst1,screg,off,pre1,pre2) \
	r7 = xptr(screg+off+32+16);	/* cosine/sine */ \
	r7 *= r3;			/* A2 = R2 * cosine/sine */ \
	r7 += r4;			/* A2 = A2 + I2 */ \
\
	r8 = xptr(screg+off+0+16);	/* cosine/sine */ \
	r8 *= r5;			/* A3 = R3 * cosine/sine */ \
	r8 += r6;			/* A3 = A3 + I3 */ \
\
	r4 *= xptr(screg+off+32+16);	/* B2 = I2 * cosine/sine */ \
	r4 -= r3;			/* B2 = B2 - R2 */ \
\
	r3 = mem7;		/* R4 */ \
	r3 *= xptr(screg+off+0+16);	/* A4 = R4 * cosine/sine */ \
	r3 -= mem8;		/* A4 = A4 - I4 */ \
\
	r7 *= xptr(screg+off+32);	/* A2 = A2 * sine (new R2) */ \
	r8 *= xptr(screg+off+0);	/* A3 = A3 * sine (new R3) */ \
	r3 *= xptr(screg+off+0);	/* A4 = A4 * sine (new R4) */ \
\
	xprefetchw(u8ptr(pre1)); \
\
	r6 *= xptr(screg+off+0+16);	/* B3 = I3 * cosine/sine */ \
	r6 -= r5;			/* B3 = B3 - R3 */ \
\
	r5 = r7;			/* R2 = R2 * 2 */ \
	r7 += r1;			/* R2 = R1 + R2 (new R1) */ \
	r1 -= r5;			/* R1 = R1 - R2 (new R2) */ \
\
	r5 = r8; 			/* R3 = R3 * 2 */ \
	r8 += r3;			/* R3 = R4 + R3 (new R3) */ \
	r3 -= r5;			/* R4 = R4 - R3 (new I4) */ \
\
	r5 = r8; \
	r8 += r7;			/* R3 = R1 + R3 (final R1) */ \
	r7 -= r5;			/* R1 = R1 - R3 (final R3) */ \
\
	r5 = mem8;		/* I4 */ \
	r5 *= xptr(screg+off+0+16);	/* B4 = I4 * cosine/sine */ \
	r5 += mem7;		/* B4 = B4 + R4 */ \
\
	r4 *= xptr(screg+off+32);	/* B2 = B2 * sine (new I2) */ \
	r6 *= xptr(screg+off+0);	/* B3 = B3 * sine (new I3) */ \
	r5 *= xptr(screg+off+0);	/* B4 = B4 * sine (new I4) */ \
\
	xprefetchw(u8ptr(pre1+pre2)); \
\
	IFB <dst1>; \
\
	r2 -= r4;			/* I1 = I1 - I2 (new I2) */ \
	r4 += r4; /* multwo */			/* I2 = I2 * 2 */ \
	r4 += r2;			/* I2 = I1 + I2 (new I1) */ \
\
	r6 -= r5;			/* I3 = I3 - I4 (new R4) */ \
	r5 += r5; /* multwo */			/* I4 = I4 * 2 */ \
	r5 += r6;			/* I4 = I3 + I4 (new I3) */ \
\
	r1 -= r6;			/* R2 = R2 - R4 (final R4) */ \
	r6 += r6; /* multwo */			/* R4 = R4 * 2 */ \
	r6 += r1;			/* R4 = R2 + R4 (final R2) */ \
\
	r2 -= r3;			/* I2 = I2 - I4 (final I4) */ \
	r3 += r3; /* multwo */			/* I4 = I4 * 2 */ \
	r3 += r2;			/* I4 = I2 + I4 (final I2) */ \
\
	r4 -= r5;			/* I1 = I1 - I3 (final I3) */ \
	r5 += r5; /* multwo */			/* I3 = I3 * 2 */ \
	r5 += r4;			/* I3 = I1 + I3 (final I1) */ \
\
	ELSE; \
\
	dst1 = r8; \
\
	r8 = r2; \
	r2 -= r4;			/* I1 = I1 - I2 (new I2) */ \
	r4 += r8;			/* I2 = I1 + I2 (new I1) */ \
\
	r8 = r6; \
	r6 -= r5;			/* I3 = I3 - I4 (new R4) */ \
	r5 += r8;			/* I4 = I3 + I4 (new I3) */ \
\
	r8 = r1; \
	r1 -= r6;			/* R2 = R2 - R4 (final R4) */ \
	r6 += r8;			/* R4 = R2 + R4 (final R2) */ \
\
	r8 = r2; \
	r2 -= r3;			/* I2 = I2 - I4 (final I4) */ \
	r3 += r8;			/* I4 = I2 + I4 (final I2) */ \
\
	r8 = r4; \
	r4 -= r5;			/* I1 = I1 - I3 (final I3) */ \
	r5 += r8;			/* I3 = I1 + I3 (final I1) */ \
\
	ENDIF; \

#endif
\
\
/* 32-bit AMD optimized versions of the above macros */ \

#if 0 // (@INSTR(,%xarch,<K8>) NE 0);
\
/* K8 Cheat sheet for scheduling dependency chains */ \
/*	      12345678901234567890123456789012345678901234567890123456789012345678901234567890 */ \
/*A2	      MMMMAAAA    MMMM */ \
/*B2	        MMMMAAAA    MMMM */ \
/*A4              MMMMAAAA */ \
/*A3                MMMMAAAA */ \
/*B3		      MMMMAAAA */ \
/*B4		        MMMMAAAA */ \
\
/*i4(r4-r3)                   AAAAMMMM */ \
/*r3(r4+r3)                   MMMMAAAAMMMM */ \
/*r4(i3-i4)                     AAAAMMMM */ \
/*i3(i3+i4)                     MMMMAAAAMMMM */ \
/*i2(i1-i2)                           AAAA */ \
/*r2(r1-r2)                             AAAA */ \
/*r1(r1+r2)                               AAAA */ \
/*i1(i1+i2)                                 AAAA */ \
\
/*I4(I2-I4)	                              AAAA */ \
/*I2(I2+I4)			              MMMM    AAAA */ \
/*R4(R2-R4)	                                AAAA			 */ \
/*R2(R2+R4)			                MMMM    AAAA */ \
/*R3(R1-R3)	                                  AAAA */ \
/*R1(R1+R3)	                                  MMMM    AAAA */ \
/*I3(I1-I3)	                                    AAAA				 */ \
/*I1(I1+I3)                                         MMMM    AAAA */ \
\

#define r4_x4c_2sc_djbunfft_mem(mem1,mem2,mem3,mem4,mem5,mem6,mem7,mem8,dst3,screg1,screg2,pre1,pre2) \
	xmm3 = xptr(screg2+16);	/* cosine/sine			;; K8 */ \
	xmm6 = mem3;		/* R2 */ \
	xmm6 *= xmm3;		/* A2 = R2 * cosine/sine	;; 1-4 */ \
	xmm0 = mem4;		/* I2 */ \
	xmm3 *= xmm0;		/* B2 = I2 * cosine/sine	;; 3-6 */ \
\
	xmm5 = xptr(screg1+16);	/* cosine/sine */ \
	xmm2 = mem7;		/* R4 */ \
	xmm2 *= xmm5;		/* A4 = R4 * cosine/sine	;; 5-8 */ \
	xmm6 += xmm0;		/* A2 = A2 + I2			;; 5-8 */ \
\
	xmm7 = mem5;		/* R3 */ \
	xmm7 *= xmm5;		/* A3 = R3 * cosine/sine	;; 7-10 */ \
	xmm3 -= mem3;		/* B2 = B2 - R2			;; 7-10 */ \
\
	xmm1 = mem6;		/* I3 */ \
	xmm5 *= xmm1;		/* B3 = I3 * cosine/sine	;; 9-12 */ \
	xmm4 = mem8;		/* I4 */ \
	xmm2 -= xmm4;		/* A4 = A4 - I4			;; 9-12 */ \
\
	xmm4 *= xptr(screg1+16);	/* B4 = I4 * cosine/sine	;; 11-14 */ \
	xmm7 += xmm1;		/* A3 = A3 + I3			;; 11-14 */ \
\
	xmm6 *= xptr(screg2);		/* A2 = A2 * sine (new R2)	;; 13-16 */ \
	xmm5 -= mem5;		/* B3 = B3 - R3			;; 13-16 */ \
\
	xmm3 *= xptr(screg2);		/* B2 = B2 * sine (new I2)	;; 15-18 */ \
	xmm4 += mem7;		/* B4 = B4 + R4			;; 15-18 */ \
\
	xprefetchw(u8ptr(pre1)); \
\
	xmm2 -= xmm7;		/* R4 = R4 - R3 (new I4)	;; 17-20 */ \
	xmm7 += xmm7; /* multwo */;			/* R3 = R3 * 2			;; 17-20 */ \
\
	xmm5 -= xmm4;		/* I3 = I3 - I4 (new R4)	;; 19-22 */ \
	xmm4 += xmm4; /* multwo */;			/* I4 = I4 * 2			;; 19-22 */ \
\
	xmm7 += xmm2;		/* R3 = R4 + R3 (new R3)	;; 21-24 */ \
	xmm0 = xptr(screg1);		/* Sine */ \
	xmm2 *= xmm0;		/* new I4 = new I4 * sine	;; 21-24 */ \
\
	xmm4 += xmm5;		/* I4 = I3 + I4 (new I3)	;; 23-26 */ \
	xmm5 *= xmm0;		/* new R4 = new R4 * sine	;; 23-26 */ \
\
	xmm1 = mem2;		/* I1 */ \
	xmm1 -= xmm3;		/* I1 = I1 - I2 (new I2)	;; 25-28 */ \
	xmm7 *= xmm0;		/* new R3 = new R3 * sine	;; 25-28 */ \
\
	xmm4 *= xmm0;		/* new I3 = new I3 * sine	;; 27-30 */ \
	xmm0 = mem1;		/* R1 */ \
	xmm0 -= xmm6;		/* R1 = R1 - R2 (new R2)	;; 27-30 */ \
\
	xmm6 += mem1;		/* R2 = R1 + R2 (new R1)	;; 29-32 */ \
\
	xmm3 += mem2;		/* I2 = I1 + I2 (new I1)	;; 31-34 */ \
\
	xprefetchw(u8ptr(pre1+pre2)); \
\
	xmm1 -= xmm2;		/* I2 = I2 - I4 (final I4)	;; 33-36 */ \
	xmm2 += xmm2; /* multwo */;			/* I4 = I4 * 2			;; 33-36 */ \
\
	xmm0 -= xmm5;		/* R2 = R2 - R4 (final R4)	;; 35-38 */ \
	xmm5 += xmm5; /* multwo */;			/* R4 = R4 * 2			;; 35-38 */ \
\
	xmm6 -= xmm7;		/* R1 = R1 - R3 (final R3)	;; 37-40 */ \
	xmm7 += xmm7; /* multwo */;			/* R3 = R3 * 2			;; 37-40 */ \
\
	xmm3 -= xmm4;		/* I1 = I1 - I3 (final I3)	;; 39-42 */ \
	xmm4 += xmm4; /* multwo */;			/* I3 = I3 * 2			;; 39-42 */ \
\
	xmm2 += xmm1;		/* I4 = I2 + I4 (final I2)	;; 41-44 */ \
	IFNB <dst3>; \
	dst3 = xmm6;						/* 41 */ \
	ENDIF; \
\
	xmm5 += xmm0;		/* R4 = R2 + R4 (final R2)	;; 43-46 */ \
	xmm7 += xmm6;		/* R3 = R1 + R3 (final R1)	;; 45-48 */ \
	xmm4 += xmm3;		/* I3 = I1 + I3 (final I1)	;; 47-50 */ \
\
\
\
/* K8 Cheat sheet for scheduling dependency chains */ \
/*	      12345678901234567890123456789012345678901234567890123456789012345678901234567890 */ \
/*B2	      MMMMAAAA    MMMM */ \
/*A2	        MMMMAAAA    MMMM */ \
/*B4		  MMMMAAAA */ \
/*A4                MMMMAAAA */ \
/*A3                  MMMMAAAA */ \
/*B3		        MMMMAAAA */ \
\
/*i4(r4-r3)                   AAAAMMMM */ \
/*r3(r4+r3)                   MMMMAAAAMMMM */ \
/*r4(i3-i4)                     AAAAMMMM */ \
/*i3(i3+i4)                     MMMMAAAAMMMM */ \
/*i2(i1-i2)                           AAAA */ \
/*r2(r1-r2)                             AAAA */ \
/*r1(r1+r2)                               AAAA */ \
/*i1(i1+i2)                                 AAAA */ \
\
/*I4(I2-I4)	                              AAAA */ \
/*I2(I2+I4)			              MMMM    AAAA */ \
/*R4(R2-R4)	                                AAAA			 */ \
/*R2(R2+R4)			                MMMM    AAAA */ \
/*R3(R1-R3)	                                  AAAA */ \
/*R1(R1+R3)	                                  MMMM    AAAA */ \
/*I3(I1-I3)	                                    AAAA				 */ \
/*I1(I1+I3)                                         MMMM    AAAA */ \
\

#define r4_x4c_2sc_djbunfft_partial_mem(r1,r2,r3,r4,r5,r6,r7,r8,mem1,mem2,mem5,mem6,dst1,screg1,screg2,pre1,pre2) \
	r1 = xptr(screg2+16);		/* cosine/sine */ \
	r1 *= r4;			/* B2 = I2 * cosine/sine	;; 1-4 */ \
\
	r1 -= r3;			/* B2 = B2 - R2			;; 5-8 */ \
\
	r3 *= xptr(screg2+16);		/* A2 = R2 * cosine/sine	;; 3-6 */ \
\
	r5 = xptr(screg1+16);		/* cosine/sine */ \
	r5 *= r8;			/* B4 = I4 * cosine/sine	;; 5-8 */ \
\
	r3 += r4;			/* A2 = A2 + I2			;; 7-10 */ \
\
	r5 += r7;			/* B4 = B4 + R4			;; 9-12 */ \
\
	r4 = xptr(screg1+16);		/* cosine/sine */ \
	r7 *= r4;			/* A4 = R4 * cosine/sine	;; 7-10 */ \
\
	r2 = mem5;		/* R3 */ \
	r2 *= r4;			/* A3 = R3 * cosine/sine	;; 9-12 */ \
\
	r7 -= r8;			/* A4 = A4 - I4			;; 11-14 */ \
	r4 *= mem6;		/* B3 = I3 * cosine/sine	;; 11-14 */ \
\
	r2 += mem6;		/* A3 = A3 + I3			;; 13-16 */ \
	r1 *= xptr(screg2);		/* B2 = B2 * sine (new I2)	;; 13-16 */ \
\
	r4 -= mem5;		/* B3 = B3 - R3			;; 15-18 */ \
	r3 *= xptr(screg2);		/* A2 = A2 * sine (new R2)	;; 15-18 */ \
\
	xprefetchw(u8ptr(pre1)); \
\
	r7 -= r2;			/* R4 = R4 - R3 (new I4)	;; 17-20 */ \
	r2 += r2; /* multwo */			/* R3 = R3 * 2			;; 17-20 */ \
\
	r4 -= r5;			/* I3 = I3 - I4 (new R4)	;; 19-22 */ \
	r5 += r5; /* multwo */			/* I4 = I4 * 2			;; 19-22 */ \
\
	r2 += r7;			/* R3 = R4 + R3 (new R3)	;; 21-24 */ \
	r6 = xptr(screg1);		/* Sine */ \
	r7 *= r6;			/* new I4 = new I4 * sine	;; 21-24 */ \
\
	r5 += r4;			/* I4 = I3 + I4 (new I3)	;; 23-26 */ \
	r4 *= r6;			/* new R4 = new R4 * sine	;; 23-26 */ \
\
	r8 = mem2;		/* I1 */ \
	r8 -= r1;			/* I1 = I1 - I2 (new I2)	;; 25-28 */ \
	r2 *= r6;			/* new R3 = new R3 * sine	;; 25-28 */ \
\
	r5 *= r6;			/* new I3 = new I3 * sine	;; 27-30 */ \
	r6 = mem1;		/* R1 */ \
	r6 -= r3;			/* R1 = R1 - R2 (new R2)	;; 27-30 */ \
\
	r3 += mem1;		/* R2 = R1 + R2 (new R1)	;; 29-32 */ \
\
	xprefetchw(u8ptr(pre1+pre2)); \
\
	r1 += mem2;		/* I2 = I1 + I2 (new I1)	;; 31-34 */ \
\
	r8 -= r7;			/* I2 = I2 - I4 (final I4)	;; 33-36 */ \
	r7 += r7; /* multwo */			/* I4 = I4 * 2			;; 33-36 */ \
\
	r6 -= r4;			/* R2 = R2 - R4 (final R4)	;; 35-38 */ \
	r4 += r4; /* multwo */			/* R4 = R4 * 2			;; 35-38 */ \
\
	r3 -= r2;			/* R1 = R1 - R3 (final R3)	;; 37-40 */ \
	r2 += r2; /* multwo */			/* R3 = R3 * 2			;; 37-40 */ \
\
	r1 -= r5;			/* I1 = I1 - I3 (final I3)	;; 39-42 */ \
	r5 += r5; /* multwo */			/* I3 = I3 * 2			;; 39-42 */ \
\
	r7 += r8;			/* I4 = I2 + I4 (final I2)	;; 41-44 */ \
	r4 += r6;			/* R4 = R2 + R4 (final R2)	;; 43-46 */ \
	r2 += r3;			/* R3 = R1 + R3 (final R1)	;; 45-48 */ \
	r5 += r1;			/* I3 = I1 + I3 (final I1)	;; 47-50 */ \
\
	IFNB <dst1>; \
	dst1 = r2; \
	ENDIF; \
\
\

#endif
\
/* 64-bit Intel and K10 implementations of the above - use the extra XMM registers */ \
/* For some reason this implementation is worse on a Pentium 4. */ \

#if 0 // X86_64 && (@INSTR(,%xarch,<CORE>) NE 0) OR (@INSTR(,%xarch,<BLEND>) NE 0) OR (@INSTR(,%xarch,<K10>) NE 0);
\
/* Theoretical best case is 44 clocks on a Core 2.  Now at 51.4 clocks. */ \
#undef r4_x4cl_four_complex_djbunfft
#define r4_x4cl_four_complex_djbunfft(srcreg,srcinc,d1,d2,screg,scoff) { \
	vec2f64 *sd0 = (vec2f64*)(srcreg); \
	vec2f64 *sd1 = (vec2f64*)(srcreg+d1); \
	vec2f64 *sd2 = (vec2f64*)(srcreg+d2); \
	vec2f64 *sd3 = (vec2f64*)(srcreg+d2+d1); \
	vec2f64 xmm8,xmm9,xmm10,xmm11,xmm12,xmm13,xmm14,xmm15; \
	xmm0 = sd3[1];	/* R4 */ \
	xmm1 = xptr(screg+scoff+16); /* cosine/sine */ \
	xmm2 = xmm0;		/* Copy R4				; 1-3 */ \
	xmm0 *= xmm1;		/* A4 = R4 * cosine/sine		; 1-5 */ \
\
	xmm5 = sd1[3];	/* I2 */ \
	xmm6 = xptr(screg+scoff+32+16); /* cosine/sine */ \
	xmm7 = xmm5;		/* Copy I2				; 2-4 */ \
	xmm5 *= xmm6;		/* B2 = I2 * cosine/sine		; 2-6 */ \
\
	xmm3 = sd2[1];	/* R3 */ \
	xmm4 = xmm3;		/* Copy R3				; 3-5 */ \
	xmm3 *= xmm1;		/* A3 = R3 * cosine/sine		; 3-7 */ \
\
	xmm8 = sd1[1];	/* R2 */ \
	xmm9 = xmm8;		/* Copy R2				; 4-6 */ \
	xmm8 *= xmm6;		/* A2 = R2 * cosine/sine		; 4-8	avail 6,10+ */ \
\
	xmm6 = sd2[3];	/* I3 */ \
	xmm10 = xmm6;		/* Copy I3				; 5-7	avail 11+ */ \
	xmm6 *= xmm1;		/* B3 = I3 * cosine/sine		; 5-9 */ \
\
	xmm11 = sd3[3];/* I4 */ \
	xmm0 -= xmm11;		/* A4 = A4 - I4 (new R4 / sine)		; 6-8 */ \
	xmm11 *= xmm1;		/* B4 = I4 * cosine/sine		; 6-10	avail 1,12+ */ \
\
	xmm5 -= xmm9;		/* B2 = B2 - R2				; 7-9	avail 1,9,12+ */ \
	xmm1 = xptr(screg+scoff+32);	/* sine					; 7  */ \
\
	xmm3 += xmm10;		/* A3 = A3 + I3 (new R3 / sine)		; 8-10	avail 9,10,12+ */ \
\
	xmm8 += xmm7;		/* A2 = A2 + I2				; 9-11	avail 9,10,7,12+ */ \
	xmm7 = sd0[3];	/* I1					; 9 */ \
\
	xmm6 -= xmm4;		/* B3 = B3 - R3 (new I3 / sine)		; 10-12	avail 9,10,4,12+ */ \
	xmm5 *= xmm1;		/* B2 = B2 * sine (new I2)		; 10-14 */ \
	xmm4 = xmm0;		/* Copy R4 / sine			; 10-12	avail 9,10,12+ */ \
\
	xmm11 += xmm2;		/* B4 = B4 + R4 (new I4 / sine)		; 11-13	avail 9,10,2,12+ */ \
\
	xmm0 -= xmm3;		/* R4 = R4 - R3 (newer I4 / sine)	; 12-14 */ \
	xmm8 *= xmm1;		/* A2 = A2 * sine (new R2)		; 12-16	avail 9,10,2,1,12+ */ \
	xmm1 = xptr(screg+scoff);	/* sine					; 12	avail 9,10,2,12+ */ \
\
	xmm3 += xmm4;		/* R3 = R4 + R3 (newer R3 / sine)	; 13-15	avail 9,10,2,4,12+ */ \
	xmm4 = xmm6;		/* Copy I3 / sine			; 13-15 */ \
	xmm2 = sd0[1];	/* R1					; 13	avail 9,10,12+ */ \
\
	xprefetchw(u8ptr(srcreg+srcinc+d2)); \
\
	xmm6 -= xmm11;		/* I3 = I3 - I4 (newer R4 / sine)	; 14-16 */ \
	xmm10 = xmm7;		/* Copy I1				; 14-16 avail 9,12+ */ \
	xmm9 = sd3[0];	/*#2 R4					; 14	avail 12+ */ \
\
	xmm7 -= xmm5;		/* I1 = I1 - I2 (newer I2)		; 15-17 */ \
	xmm0 *= xmm1;		/* newer I4 * sine			; 15-19 */ \
	xmm12 = xptr(screg+16);	/*#2 cosine/sine			; 15	avail 13+ */ \
\
	xmm11 += xmm4;		/* I4 = I3 + I4 (newer I3 / sine)	; 16-18 avail 4,13+ */ \
	xmm3 *= xmm1;		/* newer R3 * sine			; 16-20 */ \
	xmm4 = xmm2;		/* Copy R1				; 16-18	avail 13+ */ \
\
	xmm5 += xmm10;		/* I2 = I1 + I2 (newer I1)		; 17-19	avail 10,13+ */ \
	xmm6 *= xmm1;		/* newer R4 * sine			; 17-21 */ \
	xmm13 = sd1[2];	/*#2 I2					; 17	avail 10,14+ */ \
\
	xmm2 += xmm8;		/* R2 = R1 + R2 (newer R1)		; 18-20 */ \
	xmm14 = xptr(screg+32+16);	/*#2 cosine/sine			; 18	avail 10,15 */ \
\
	xmm4 -= xmm8;		/* R1 = R1 - R2 (newer R2)		; 19-21	avail 10,8,15 */ \
	xmm11 *= xmm1;		/* newer I3 * sine			; 19-23	avail 10,8,1,15 */ \
	xmm1 = xmm7;		/* Copy I2				; 19-21	avail 10,8,15 */ \
\
	xprefetchw(u8ptr(srcreg+srcinc+d2+d1)); \
\
	xmm7 -= xmm0;		/* I2 = I2 - I4 (final I4)		; 20-22	avail 10,8,15 storable 7 */ \
	xmm15 = xmm9;		/*#2 Copy R4				; 20-22	avail 10,8 storable 7 */ \
\
	xmm8 = xmm2;		/* Copy R1				; 21-23	avail 10 storable 7 */ \
	xmm2 -= xmm3;		/* R1 = R1 - R3 (final R3)		; 21-23	avail 10 storable 7,2 */ \
\
	xmm0 += xmm1;		/* I4 = I2 + I4 (final I2)		; 22-24	avail 10,1 storable 7,2,0 */ \
	xmm1 = xmm4;		/* Copy R2				; 22-24	avail 10 storable 7,2,0 */ \
\
	xmm4 -= xmm6;		/* R2 = R2 - R4 (final R4)		; 23-25	avail 10 storable 7,2,0,4 */ \
	xmm9 *= xmm12;		/*#2 A4 = R4 * cosine/sine		; 23-27 */ \
	xmm10 = xmm13;		/*#2 Copy I2				; 23-25	avail none storable 7,2,0,4 */ \
	sd3[3] = xmm7;	/* Save I4				; 23	avail 7 storable 2,0,4 */ \
\
	xmm3 += xmm8;		/* R3 = R1 + R3 (final R1)		; 24-26	avail 7,8 storable 2,0,4,3 */ \
	xmm13 *= xmm14;		/*#2 B2 = I2 * cosine/sine		; 24-28 */ \
	xmm7 = xmm5;		/* Copy I1				; 24-26	avail 8 storable 2,0,4,3 */ \
	sd1[1] = xmm2;	/* Save R3				; 24	avail 8,2 storable 0,4,3 */ \
\
	xmm6 += xmm1;		/* R4 = R2 + R4 (final R2)		; 25-27	avail 8,2,1 storable 0,4,3,6 */ \
	xmm2 = sd2[0];	/*#2 R3					; 25	avail 8,1 storable 0,4,3,6 */ \
	xmm8 = xmm2;		/*#2 Copy R3				; 25-27	avail 1 storable 0,4,3,6 */ \
	xmm2 *= xmm12;		/*#2 A3 = R3 * cosine/sine		; 25-29 */ \
\
	xmm5 -= xmm11;		/* I1 = I1 - I3 (final I3)		; 26-28	avail 1 storable 0,4,3,6,5 */ \
	xmm1 = sd1[0];	/*#2 R2					; 26	avail none storable 0,4,3,6,5 */ \
	sd3[1] = xmm4;	/* Save R4				; 26	avail 4 storable 0,3,6,5 */ \
	xmm4 = xmm1;		/*#2 Copy R2				; 26-28	avail none storable 0,3,6,5 */ \
	xmm1 *= xmm14;		/*#2 A2 = R2 * cosine/sine		; 26-30	avail 14 storable 0,3,6,5 */ \
\
	xmm11 += xmm7;		/* I3 = I1 + I3 (final I1)		; 27-29	avail 14,7 storable 0,3,6,5,11 */ \
	xmm14 = sd2[2];	/*#2 I3					; 27	avail 7 storable 0,3,6,5,11 */ \
	xmm7 = xmm14;		/*#2 Copy I3				; 27-29	avail none storable 0,3,6,5,11 */ \
	xmm14 *= xmm12;		/*#2 B3 = I3 * cosine/sine		; 27-31 */ \
	sd1[0] = xmm3;	/* Save R1				; 27	avail 3 storable 0,6,5,11 */ \
\
	xmm3 = sd3[2]; /*#2 I4					; 28	avail none storable 0,6,5,11 */ \
	xmm9 -= xmm3;		/*#2 A4 = A4 - I4 (new R4 / sine)	; 28-30 */ \
	xmm3 *= xmm12;		/*#2 B4 = I4 * cosine/sine		; 28-32	avail 12 storable 0,6,5,11 */ \
	sd3[0] = xmm6;	/* Save R2				; 28	avail 12,6 storable 0,5,11 */ \
\
	xmm13 -= xmm4;		/*#2 B2 = B2 - R2			; 29-31	avail 12,6,4 storable 0,5,11 */ \
	xmm6 = xptr(screg+32);	/*#2 sine				; 29	avail 12,4 storable 0,5,11 */ \
	sd1[3] = xmm5;	/* Save I3				; 29	avail 12,4,5 storable 0,11 */ \
\
	xmm2 += xmm7;		/*#2 A3 = A3 + I3 (new R3 / sine)	; 30-32	avail 12,4,5,7 storable 0,11 */ \
	xmm5 = sd0[2];	/*#2 I1					; 30	avail 12,4,7 storable 0,11 */ \
	sd1[2] = xmm11;	/* Save I1				; 30	avail 12,4,7,11 storable 0 */ \
\
	xmm1 += xmm10;		/*#2 A2 = A2 + I2			; 31-33	avail 12,4,7,11,10 storable 0 */ \
	sd3[2] = xmm0;	/* Save I2				; 25	avail 12,4,7,11,10,0 */ \
\
	xmm14 -= xmm8;		/*#2 B3 = B3 - R3 (new I3 / sine)	; 32-34	avail 12,4,7,11,10,0,8 */ \
	xmm13 *= xmm6;		/*#2 B2 = B2 * sine (new I2)		; 32-36 */ \
	xmm4 = xmm9;		/*#2 Copy R4 / sine			; 32-34	avail 12,7,11,10,0,8 */ \
\
	xmm3 += xmm15;		/*#2 B4 = B4 + R4 (new I4 / sine)	; 33-35	avail 12,7,11,10,0,8,15 */ \
\
	xmm9 -= xmm2;		/*#2 R4 = R4 - R3 (newer I4 / sine)	; 34-36 */ \
	xmm1 *= xmm6;		/*#2 A2 = A2 * sine (new R2)		; 34-38	avail 12,7,11,10,0,8,15,6 */ \
	xmm6 = xptr(screg);		/*#2 sine				; 34	avail 12,7,11,10,0,8,15 */ \
\
	xmm2 += xmm4;		/*#2 R3 = R4 + R3 (newer R3 / sine)	; 35-37	avail 12,7,11,10,0,8,15,4 */ \
	xmm4 = xmm14;		/*#2 Copy I3 / sine			; 35-37	avail 12,7,11,10,0,8,15 */ \
	xmm7 = sd0[0];		/*#2 R1					; 35	avail 12,11,10,0,8,15 */ \
\
	xprefetchw(u8ptr(srcreg+srcinc)); \
\
	xmm14 -= xmm3;		/*#2 I3 = I3 - I4 (newer R4 / sine)	; 36-38 */ \
	xmm0 = xmm5;		/*#2 Copy I1				; 36-38	avail 12,11,10,8,15 */ \
\
	xmm5 -= xmm13;		/*#2 I1 = I1 - I2 (newer I2)		; 37-39 */ \
	xmm9 *= xmm6;		/*#2 newer I4 * sine			; 37-41 */ \
\
	xmm3 += xmm4;		/*#2 I4 = I3 + I4 (newer I3 / sine)	; 38-40	avail 12,11,10,8,15,0 */ \
	xmm2 *= xmm6;		/*#2 newer R3 * sine			; 38-42 */ \
	xmm4 = xmm7;		/*#2 Copy R1				; 38-40	avail 12,11,10,8,15 */ \
\
	xmm13 += xmm0;		/*#2 I2 = I1 + I2 (newer I1)		; 39-41	avail 12,11,10,8,15,0 */ \
	xmm14 *= xmm6;		/*#2 newer R4 * sine			; 39-43 */ \
\
	xmm7 += xmm1;		/*#2 R2 = R1 + R2 (newer R1)		; 40-42 */ \
	xmm0 = xmm5;		/*#2 Copy I2				; 40-42	avail 12,11,10,8,15 */ \
\
	xmm4 -= xmm1;		/*#2 R1 = R1 - R2 (newer R2)		; 41-43	avail 12,11,10,8,15,1 */ \
	xmm3 *= xmm6;		/*#2 newer I3 * sine			; 41-45	avail 12,11,10,8,15,1,6 */ \
\
	xprefetchw(u8ptr(srcreg+srcinc+d1)); \
\
	xmm5 -= xmm9;		/*#2 I2 = I2 - I4 (final I4)		; 42-44	avail 12,11,10,8,15,6 storable 5 */ \
	xmm1 = xmm13;		/*#2 Copy I1				; 42-44	avail 12,11,10,8,15,6 storable 5 */ \
\
	xmm6 = xmm7;		/*#2 Copy R1				; 43-45	avail 12,11,10,8,15 storable 5 */ \
	xmm7 -= xmm2;		/*#2 R1 = R1 - R3 (final R3)		; 43-45	avail 12,11,10,8,15 storable 5,7 */ \
\
	xmm9 += xmm0;		/*#2 I4 = I2 + I4 (final I2)		; 44-46	avail 12,11,10,8,15,0 storable 5,7,9 */ \
	xmm0 = xmm4;		/*#2 Copy R2				; 44-46	avail 12,11,10,8,15 storable 5,7,9 */ \
\
	xmm4 -= xmm14;		/*#2 R2 = R2 - R4 (final R4)		; 45-47	avail 12,11,10,8,15 storable 5,7,9,4 */ \
	sd2[3] = xmm5;	/*#2 Save I4				; 45	avail 12,11,10,8,15,5 storable 7,9,4 */ \
\
	xmm2 += xmm6;		/*#2 R3 = R1 + R3 (final R1)		; 46-48	avail 12,11,10,8,15,5,6 storable 7,9,4,2 */ \
	sd0[1] = xmm7;	/*#2 Save R3				; 46	avail 12,11,10,8,15,5,6,7 storable 9,4,2 */ \
\
	xmm14 += xmm0;		/*#2 R4 = R2 + R4 (final R2)		; 47-49	avail 12,11,10,8,15,5,6,7,0 storable 9,4,2,14 */ \
	sd2[2] = xmm9;	/*#2 Save I2				; 47	avail 12,11,10,8,15,5,6,7,0,9 storable 4,2,14 */ \
\
	xmm13 -= xmm3;		/*#2 I1 = I1 - I3 (final I3)		; 48-50	avail 12,11,10,8,15,5,6,7,0,9 storable 4,2,14,13 */ \
	sd2[1] = xmm4;	/*#2 Save R4				; 48	avail 12,11,10,8,15,5,6,7,0,9,4 storable 2,14,13 */ \
\
	xmm3 += xmm1;		/*#2 I3 = I1 + I3 (final I1)		; 49-51	avail 12,11,10,8,15,5,6,7,0,9,4,1 storable 2,14,13,3 */ \
	sd0[0] = xmm2;		/*#2 Save R1				; 49	avail 12,11,10,8,15,5,6,7,0,9,4,1,2 storable 14,13,3 */ \
\
	sd2[0] = xmm14;	/*#2 Save R2				; 50 */ \
	sd0[3] = xmm13;	/*#2 Save I3				; 51 */ \
	sd0[2] = xmm3;	/*#2 Save I1				; 52 */ \
\
	srcreg += srcinc; \
}
#endif

/* 64-bit AMD K8 optimized versions of the above macros.  The same as the 32-bit */ \
/* version with some constants preloaded.  We could probably do better. */ \

#if 0 // (@INSTR(,%xarch,<K8>) NE 0);
#if 0 //DEF X86_64;
\

#define r4_x4c_2sc_djbunfft_mem_preload \
	xmm15 = XMM_TWO; \
\
\

#define r4_x4c_2sc_djbunfft_mem(mem1,mem2,mem3,mem4,mem5,mem6,mem7,mem8,dst3,screg1,screg2,pre1,pre2) \
	xmm3 = xptr(screg2+16);	/* cosine/sine			;; K8 */ \
	xmm6 = mem3;		/* R2 */ \
	xmm6 *= xmm3;		/* A2 = R2 * cosine/sine	;; 1-4 */ \
	xmm0 = mem4;		/* I2 */ \
	xmm3 *= xmm0;		/* B2 = I2 * cosine/sine	;; 3-6 */ \
\
	xmm5 = xptr(screg1+16);	/* cosine/sine */ \
	xmm2 = mem7;		/* R4 */ \
	xmm2 *= xmm5;		/* A4 = R4 * cosine/sine	;; 5-8 */ \
	xmm6 += xmm0;		/* A2 = A2 + I2			;; 5-8 */ \
\
	xmm7 = mem5;		/* R3 */ \
	xmm7 *= xmm5;		/* A3 = R3 * cosine/sine	;; 7-10 */ \
	xmm3 -= mem3;		/* B2 = B2 - R2			;; 7-10 */ \
\
	xmm1 = mem6;		/* I3 */ \
	xmm5 *= xmm1;		/* B3 = I3 * cosine/sine	;; 9-12 */ \
	xmm4 = mem8;		/* I4 */ \
	xmm2 -= xmm4;		/* A4 = A4 - I4			;; 9-12 */ \
\
	xmm4 *= xptr(screg1+16);	/* B4 = I4 * cosine/sine	;; 11-14 */ \
	xmm7 += xmm1;		/* A3 = A3 + I3			;; 11-14 */ \
\
	xmm6 *= xptr(screg2);		/* A2 = A2 * sine (new R2)	;; 13-16 */ \
	xmm5 -= mem5;		/* B3 = B3 - R3			;; 13-16 */ \
\
	xmm3 *= xptr(screg2);		/* B2 = B2 * sine (new I2)	;; 15-18 */ \
	xmm4 += mem7;		/* B4 = B4 + R4			;; 15-18 */ \
\
	xprefetchw(u8ptr(pre1)); \
\
	xmm2 -= xmm7;		/* R4 = R4 - R3 (new I4)	;; 17-20 */ \
	xmm7 *= xmm15;		/* R3 = R3 * 2			;; 17-20 */ \
\
	xmm5 -= xmm4;		/* I3 = I3 - I4 (new R4)	;; 19-22 */ \
	xmm4 *= xmm15;		/* I4 = I4 * 2			;; 19-22 */ \
\
	xmm7 += xmm2;		/* R3 = R4 + R3 (new R3)	;; 21-24 */ \
	xmm0 = xptr(screg1);		/* Sine */ \
	xmm2 *= xmm0;		/* new I4 = new I4 * sine	;; 21-24 */ \
\
	xmm4 += xmm5;		/* I4 = I3 + I4 (new I3)	;; 23-26 */ \
	xmm5 *= xmm0;		/* new R4 = new R4 * sine	;; 23-26 */ \
\
	xmm1 = mem2;		/* I1 */ \
	xmm1 -= xmm3;		/* I1 = I1 - I2 (new I2)	;; 25-28 */ \
	xmm7 *= xmm0;		/* new R3 = new R3 * sine	;; 25-28 */ \
\
	xmm4 *= xmm0;		/* new I3 = new I3 * sine	;; 27-30 */ \
	xmm0 = mem1;		/* R1 */ \
	xmm0 -= xmm6;		/* R1 = R1 - R2 (new R2)	;; 27-30 */ \
\
	xmm6 += mem1;		/* R2 = R1 + R2 (new R1)	;; 29-32 */ \
\
	xmm3 += mem2;		/* I2 = I1 + I2 (new I1)	;; 31-34 */ \
\
	xprefetchw(u8ptr(pre1+pre2)); \
\
	xmm1 -= xmm2;		/* I2 = I2 - I4 (final I4)	;; 33-36 */ \
	xmm2 *= xmm15;		/* I4 = I4 * 2			;; 33-36 */ \
\
	xmm0 -= xmm5;		/* R2 = R2 - R4 (final R4)	;; 35-38 */ \
	xmm5 *= xmm15;		/* R4 = R4 * 2			;; 35-38 */ \
\
	xmm6 -= xmm7;		/* R1 = R1 - R3 (final R3)	;; 37-40 */ \
	xmm7 *= xmm15;		/* R3 = R3 * 2			;; 37-40 */ \
\
	xmm3 -= xmm4;		/* I1 = I1 - I3 (final I3)	;; 39-42 */ \
	xmm4 *= xmm15;		/* I3 = I3 * 2			;; 39-42 */ \
\
	xmm2 += xmm1;		/* I4 = I2 + I4 (final I2)	;; 41-44 */ \
	IFNB <dst3>; \
	dst3 = xmm6;						/* 41 */ \
	ENDIF; \
\
	xmm5 += xmm0;		/* R4 = R2 + R4 (final R2)	;; 43-46 */ \
	xmm7 += xmm6;		/* R3 = R1 + R3 (final R1)	;; 45-48 */ \
	xmm4 += xmm3;		/* I3 = I1 + I3 (final I1)	;; 47-50 */ \
\
\
\

#define r4_x4c_2sc_djbunfft_partial_mem_preload \
	xmm15 = XMM_TWO; \
\
\

#define r4_x4c_2sc_djbunfft_partial_mem(r1,r2,r3,r4,r5,r6,r7,r8,mem1,mem2,mem5,mem6,dst1,screg1,screg2,pre1,pre2) \
	r1 = xptr(screg2+16);		/* cosine/sine */ \
	r1 *= r4;			/* B2 = I2 * cosine/sine	;; 1-4 */ \
\
	r1 -= r3;			/* B2 = B2 - R2			;; 5-8 */ \
\
	r3 *= xptr(screg2+16);		/* A2 = R2 * cosine/sine	;; 3-6 */ \
\
	r5 = xptr(screg1+16);		/* cosine/sine */ \
	r5 *= r8;			/* B4 = I4 * cosine/sine	;; 5-8 */ \
\
	r3 += r4;			/* A2 = A2 + I2			;; 7-10 */ \
\
	r5 += r7;			/* B4 = B4 + R4			;; 9-12 */ \
\
	r4 = xptr(screg1+16);		/* cosine/sine */ \
	r7 *= r4;			/* A4 = R4 * cosine/sine	;; 7-10 */ \
\
	r2 = mem5;		/* R3 */ \
	r2 *= r4;			/* A3 = R3 * cosine/sine	;; 9-12 */ \
\
	r7 -= r8;			/* A4 = A4 - I4			;; 11-14 */ \
	r4 *= mem6;		/* B3 = I3 * cosine/sine	;; 11-14 */ \
\
	r2 += mem6;		/* A3 = A3 + I3			;; 13-16 */ \
	r1 *= xptr(screg2);		/* B2 = B2 * sine (new I2)	;; 13-16 */ \
\
	r4 -= mem5;		/* B3 = B3 - R3			;; 15-18 */ \
	r3 *= xptr(screg2);		/* A2 = A2 * sine (new R2)	;; 15-18 */ \
\
	xprefetchw(u8ptr(pre1)); \
\
	r7 -= r2;			/* R4 = R4 - R3 (new I4)	;; 17-20 */ \
	r2 *= xmm15;		/* R3 = R3 * 2			;; 17-20 */ \
\
	r4 -= r5;			/* I3 = I3 - I4 (new R4)	;; 19-22 */ \
	r5 *= xmm15;		/* I4 = I4 * 2			;; 19-22 */ \
\
	r2 += r7;			/* R3 = R4 + R3 (new R3)	;; 21-24 */ \
	r6 = xptr(screg1);		/* Sine */ \
	r7 *= r6;			/* new I4 = new I4 * sine	;; 21-24 */ \
\
	r5 += r4;			/* I4 = I3 + I4 (new I3)	;; 23-26 */ \
	r4 *= r6;			/* new R4 = new R4 * sine	;; 23-26 */ \
\
	r8 = mem2;		/* I1 */ \
	r8 -= r1;			/* I1 = I1 - I2 (new I2)	;; 25-28 */ \
	r2 *= r6;			/* new R3 = new R3 * sine	;; 25-28 */ \
\
	r5 *= r6;			/* new I3 = new I3 * sine	;; 27-30 */ \
	r6 = mem1;		/* R1 */ \
	r6 -= r3;			/* R1 = R1 - R2 (new R2)	;; 27-30 */ \
\
	r3 += mem1;		/* R2 = R1 + R2 (new R1)	;; 29-32 */ \
\
	xprefetchw(u8ptr(pre1+pre2)); \
\
	r1 += mem2;		/* I2 = I1 + I2 (new I1)	;; 31-34 */ \
\
	r8 -= r7;			/* I2 = I2 - I4 (final I4)	;; 33-36 */ \
	r7 *= xmm15;		/* I4 = I4 * 2			;; 33-36 */ \
\
	r6 -= r4;			/* R2 = R2 - R4 (final R4)	;; 35-38 */ \
	r4 *= xmm15;		/* R4 = R4 * 2			;; 35-38 */ \
\
	r3 -= r2;			/* R1 = R1 - R3 (final R3)	;; 37-40 */ \
	r2 *= xmm15;		/* R3 = R3 * 2			;; 37-40 */ \
\
	r1 -= r5;			/* I1 = I1 - I3 (final I3)	;; 39-42 */ \
	r5 *= xmm15;		/* I3 = I3 * 2			;; 39-42 */ \
\
	r7 += r8;			/* I4 = I2 + I4 (final I2)	;; 41-44 */ \
	r4 += r6;			/* R4 = R2 + R4 (final R2)	;; 43-46 */ \
	r2 += r3;			/* R3 = R1 + R3 (final R1)	;; 45-48 */ \
	r5 += r1;			/* I3 = I1 + I3 (final I1)	;; 47-50 */ \
\
	IFNB <dst1>; \
	dst1 = r2; \
	ENDIF; \
\
\

#endif
#endif
\
/* */ \
/* ************************************* four-complex-first-fft variants ****************************************** */ \
/* */ \
\
/* This code applies the roots-of-minus-1 premultipliers in a negacyclic FFT. */ \
/* It also applies the three sin/cos multipliers after a radix-4 butterfly.  We save memory */ \
/* by splitting the roots-of-minus-1 premultipliers such that every macro uses the  */ \
/* same premultiplier data and we have 4 sin/cos postmultipliers.  Every sin/cos postmultiplier, */ \
/* a very big table anyway, is multiplied by the other part of the split roots-of-minus-1. */ \
/* The common premuliplier data is 1, .924+.383i, SQRTHALF+SQRTHALFi, .383+.924i. */ \
/* This scheme is used by our simple radix-4 DJB "r4" FFTs. */ \
\
/*r4_x2cl_four_complex_first_fft4 MACRO srcreg,srcinc,d1 */ \
/*	r4_x2cl_four_complex_first_fft4_cmn srcreg,rbx,srcinc,d1,rdi */ \
/*	ENDM */ \

/* r4_x2cl_four_complex_first_fft4_cmn unused */

#define r4_x2cl_four_complex_first_fft4_scratch(srcreg,srcinc,d1,screg) { \
	vec2f64 *sd0 = (vec2f64*)(srcreg); \
	vec2f64 *sd1 = (vec2f64*)(srcreg+d1); \
	vec2f64 *sc = (vec2f64*)(screg); \
	vec2f64 temp2; \
	xmm7 = CONST2_P924; \
	xmm1 = CONST2_P383; \
	vec2f64 s0 = sd1[0], s1 = sd1[1], s2 = sd1[2], s3 = sd1[3]; \
\
	xmm4 = s0;				/* R2 */ \
	xmm0 = s2;				/* R6 (I2) */ \
	xmm4 *= xmm7;			/* A2 = R2 * .924 */ \
	xmm4 -= (xmm0 * xmm1);	/* A2 = A2 - (I2 * .383) (new R2) */ \
\
	xmm6 = s0;				/* R2 */ \
	xmm2 = s2;				/* R6 (I2) */ \
	xmm2 *= xmm7;			/* D2 = I2 * .924 */ \
	xmm6 = xmm6 * xmm1 + xmm2;	/* B2 = R2 * .383 + D2 (new I2) */ \
\
	xmm5 = s1;				/* R4 */ \
	xmm0 = s3;				/* R8 (I4) */ \
	xmm5 *= xmm1;			/* A4 = R4 * .383 */ \
	xmm5 -= (xmm0 * xmm7);	/* A4 = A4 - (I4 * .924) (new R4) */ \
\
	xmm1 *= s3;				/* D4 = I4 * .383 */ \
	xmm7 = xmm7 * s1 + xmm1;	/* B4 = R4 * .924 + D4 (new I4) */ \
\
	xprefetchw(u8ptr(srcreg+srcinc)); \
\
	s0 = sd0[0]; s1 = sd0[1]; s2 = sd0[2]; s3 = sd0[3]; \
\
	xmm1 = s1;				/* R3 */ \
	xmm3 = s3;				/* R7 (I3) */ \
	xmm1 -= xmm3;			/* A3 = R3 - I3 */ \
	xmm3 += s1;				/* B3 = R3 + I3 */ \
	xmm1 *= XMM_SQRTHALF2;		/* A3 = A3 * SQRTHALF (new R3) */ \
	xmm3 *= XMM_SQRTHALF2;		/* B3 = B3 * SQRTHALF (new I3) */ \
\
	xprefetchw(u8ptr(srcreg+srcinc+d1)); \
\
	xmm0 = xmm4;			/* Copy R2 */ \
	xmm4 -= xmm5;			/* R2 = R2 - R4 (new R4) */ \
	xmm5 += xmm0;			/* R4 = R2 + R4 (new R2) */ \
\
	xmm2 = xmm6;			/* Copy I2 */ \
	xmm6 -= xmm7;			/* I2 = I2 - I4 (new I4) */ \
	xmm7 += xmm2;			/* I4 = I2 + I4 (new I2) */ \
\
	xmm0 = s0;				/* R1 */ \
	xmm0 -= xmm1;			/* R1 = R1 - R3 (new R3) */ \
	xmm1 += s0;				/* R3 = R1 + R3 (new R1) */ \
\
	xmm2 = s2;				/* R5 (I1) */ \
	xmm2 -= xmm3;			/* I1 = I1 - I3 (new I3) */ \
	xmm3 += s2;				/* I3 = I1 + I3 (new I1) */ \
\
	temp2 = xmm0; \
	xmm0 -= xmm6;			/* R3 = R3 - I4 (new R3) */ \
	xmm6 += temp2;			/* I4 = R3 + I4 (new R4) */ \
	temp2 = xmm2; \
	xmm2 -= xmm4;			/* I3 = I3 - R4 (new I4) */ \
	xmm4 += temp2;			/* R4 = I3 + R4 (new I3) */ \
	temp2 = xmm1; \
	xmm1 -= xmm5;			/* R1 = R1 - R2 (new R2) */ \
	xmm5 += temp2;			/* R2 = R1 + R2 (new R1) */ \
	temp2 = xmm3; \
	xmm3 -= xmm7;			/* I1 = I1 - I2 (new I2) */ \
	xmm7 += temp2;			/* I2 = I1 + I2 (new I1) */ \
\
	vec2f64 saved_R1 = xmm5;	/* Save R1 */ \
	vec2f64 saved_I1 = xmm7;	/* Save I1 */ \
\
	s1 = sc[3];					/* cosine/sine */ \
	xmm5 = s1 * xmm0 - xmm4;	/* A3 = R3 * cosine/sine - I3 */ \
	xmm4 = xmm4 * s1 + xmm0;	/* B3 = I3 * cosine/sine + R3 */ \
\
	s2 = sc[5];					/* cosine/sine */ \
	xmm7 = s2 * xmm1 - xmm3;	/* A2 = R2 * cosine/sine - I2 */ \
	xmm3 = xmm3 * s2 + xmm1;	/* B2 = I2 * cosine/sine + R2 */ \
\
	s3 = sc[7];					/* cosine/sine */ \
	xmm0 = s3 * xmm6 - xmm2;	/* A4 = R4 * cosine/sine - I4 */ \
	xmm2 = xmm2 * s3 + xmm6;	/* B4 = I4 * cosine/sine + R4 */ \
\
	s0 = sc[1];					/* cosine/sine */ \
	xmm6 = saved_I1;			/* Restore I1 */ \
	xmm1 = s0 * saved_R1 - xmm6;	/* A1 = R1 * cosine/sine - I1 */ \
	xmm6 = xmm6 * s0 + saved_R1;	/* B1 = I1 * cosine/sine + R1 */ \
\
	s1 = sc[2]; \
	xmm5 *= s1;				/* A3 = A3 * sine (final R3) */ \
	xmm4 *= s1;				/* B3 = B3 * sine (final I3) */ \
	s2 = sc[4]; \
	xmm7 *= s2;				/* A2 = A2 * sine (final R2) */ \
	xmm3 *= s2;				/* B2 = B2 * sine (final I2) */ \
	s3 = sc[6]; \
	xmm0 *= s3;				/* A4 = A4 * sine (final R4) */ \
	xmm2 *= s3;				/* B4 = B4 * sine (final I4) */ \
	s0 = sc[0]; \
	xmm1 *= s0;				/* A1 = A1 * sine (final R1) */ \
	xmm6 *= s0;				/* B1 = B1 * sine (final I1) */ \
\
	sd1[0] = xmm5;			/* Save R3 */ \
	sd1[1] = xmm4;			/* Save I3 */ \
	sd1[2] = xmm0;			/* Save R4 */ \
	sd1[3] = xmm2;			/* Save I4 */ \
	sd0[0] = xmm1;			/* Save R1 */ \
	sd0[1] = xmm6;			/* Save I1 */ \
	sd0[2] = xmm7;			/* Save R2 */ \
	sd0[3] = xmm3;			/* Save I2 */ \
	srcreg += srcinc; \
}


/* This code applies the roots-of-minus-1 premultipliers in a negacyclic FFT. */ \
/* It also applies the three sin/cos multipliers after the first radix-4 butterfly. */ \
/* We save memory by splitting the roots-of-minus-1 premultipliers such that every */ \
/* macro uses the same premultiplier data.  The rest of the premultiplier is applied */ \
/* at the end of pass 1.  This scheme is used by the r4delay FFTs. */ \
\
/*r4_x2cl_four_complex_first_djbfft MACRO srcreg,srcinc,d1,screg,pmreg */ \
/*	r4_x2cl_four_complex_first_djbfft_cmn srcreg,rbx,srcinc,d1,screg,pmreg */ \
/*	ENDM */ \

/* r4_x2cl_four_complex_first_djbfft_cmn1 unused */

#define r4_x2cl_four_complex_first_djbfft_scratch(srcreg,srcinc,d1,screg,pmreg) { \
	vec2f64 *sd0 = (vec2f64*)(srcreg); \
	vec2f64 *sd1 = (vec2f64*)(srcreg+d1); \
	vec2f64 *sc = (vec2f64*)(screg); \
	vec2f64 *pm = (vec2f64*)(pmreg); \
	vec2f64 s1, s2, temp2; \
	s1 = pm[1]; \
	xmm1 = sd0[0];				/* R1 */ \
	xmm2 = sd0[2];				/* R5 (I1) */ \
	xmm7 = s1 * xmm1 - xmm2;	/* A1 = R1 * cosine/sine - I1 */ \
	xmm2 = xmm2 * s1 + xmm1;	/* B1 = I1 * cosine/sine + R1 */ \
\
	s1 = pm[3]; \
	xmm3 = sd1[0];				/* R2 */ \
	xmm4 = sd1[2];				/* R6 (I2) */ \
	xmm6 = s1 * xmm3 - xmm4;	/* A2 = R2 * cosine/sine - I2 */ \
	xmm4 = xmm4 * s1 + xmm3;	/* B2 = I2 * cosine/sine + R2 */ \
\
	s1 = pm[5]; \
	xmm0 = sd0[1];				/* R3 */ \
	xmm1 = sd0[3];				/* R7 (I3) */ \
	xmm5 = s1 * xmm0 - xmm1;	/* A3 = R3 * cosine/sine - I3 */ \
	xmm1 = xmm1 * s1 + xmm0;	/* B3 = I3 * cosine/sine + R3 */ \
\
	s1 = pm[7]; \
	xmm0 = sd1[3];				/* R8 (I4) */ \
	xmm3 = s1 * sd1[1] - xmm0;	/* A4 = R4 * cosine/sine - I4 */ \
	xmm0 = xmm0 * s1 + sd1[1];	/* B4 = I4 * cosine/sine + R4 */ \
\
	s1 = pm[0]; \
	xmm7 *= s1;				/* A1 = A1 * sine */ \
	xmm2 *= s1;				/* B1 = B1 * sine */ \
	s1 = pm[2]; \
	xmm6 *= s1;				/* A2 = A2 * sine */ \
	xmm4 *= s1;				/* B2 = B2 * sine */ \
	s1 = pm[4]; \
	xmm5 *= s1;				/* A3 = A3 * sine */ \
	xmm1 *= s1;				/* B3 = B3 * sine */ \
	s1 = pm[6]; \
	xmm3 *= s1;				/* A4 = A4 * sine */ \
	xmm0 *= s1;				/* B4 = B4 * sine */ \
\
	xprefetchw(u8ptr(srcreg+srcinc)); \
	xprefetchw(u8ptr(srcreg+srcinc+d1)); \
\
	temp2 = xmm6; \
	xmm6 -= xmm3;			/* R2 = R2 - R4 (new R4) */ \
	xmm3 += temp2;			/* R4 = R2 + R4 (new R2) */ \
\
	temp2 = xmm4; \
	xmm4 -= xmm0;			/* I2 = I2 - I4 (new I4) */ \
	xmm0 += temp2;			/* I4 = I2 + I4 (new I2) */ \
\
	temp2 = xmm7; \
	xmm7 -= xmm5;			/* R1 = R1 - R3 (new R3) */ \
	xmm5 += temp2;			/* R3 = R1 + R3 (new R1) */ \
\
	temp2 = xmm2; \
	xmm2 -= xmm1;			/* I1 = I1 - I3 (new I3) */ \
	xmm1 += temp2;			/* I3 = I1 + I3 (new I1) */ \
\
	temp2 = xmm7; \
	xmm7 -= xmm4;			/* R3 = R3 - I4 (new R3) */ \
	xmm4 += temp2;			/* I4 = R3 + I4 (new R4) */ \
	temp2 = xmm2; \
	xmm2 -= xmm6;			/* I3 = I3 - R4 (new I4) */ \
	xmm6 += temp2;			/* R4 = I3 + R4 (new I3) */ \
	temp2 = xmm5; \
	xmm5 -= xmm3;			/* R1 = R1 - R2 (new R2) */ \
	xmm3 += temp2;			/* R2 = R1 + R2 (new R1) */ \
	temp2 = xmm1; \
	xmm1 -= xmm0;			/* I1 = I1 - I2 (new I2) */ \
	xmm0 += temp2;			/* I2 = I1 + I2 (new I1) */ \
\
	sd0[0] = xmm3;			/* Save R1 */ \
	sd0[1] = xmm0;			/* Save I1 */ \
\
	s1 = sc[1]; \
	s2 = sc[3]; \
	xmm3 = s1 * xmm7 - xmm6;	/* A3 = R3 * cosine/sine - I3 */ \
	xmm6 = xmm6 * s1 + xmm7;	/* B3 = I3 * cosine/sine + R3 */ \
\
	xmm0 = s2 * xmm5 - xmm1;	/* A2 = R2 * cosine/sine - I2 */ \
	xmm1 = xmm1 * s2 + xmm5;	/* B2 = I2 * cosine/sine + R2 */ \
\
	xmm7 = s1 * xmm4 + xmm2;	/* A4 = R4 * cosine/sine + I4 */ \
	xmm2 = xmm2 * s1 - xmm4;	/* B4 = I4 * cosine/sine - R4 */ \
\
	xmm5 = sc[0]; \
	xmm3 *= xmm5;			/* A3 = A3 * sine (final R3) */ \
	xmm6 *= xmm5;			/* B3 = B3 * sine (final I3) */ \
	s1 = sc[2]; \
	xmm0 *= s1;				/* A2 = A2 * sine (final R2) */ \
	xmm1 *= s1;				/* B2 = B2 * sine (final I2) */ \
	xmm7 *= xmm5;			/* A4 = A4 * sine (final R4) */ \
	xmm2 *= xmm5;			/* B4 = B4 * sine (final I4) */ \
\
	sd1[0] = xmm3;			/* Save R3 */ \
	sd1[1] = xmm6;			/* Save I3 */ \
	sd1[2] = xmm7;			/* Save R4 */ \
	sd1[3] = xmm2;			/* Save I4 */ \
	sd0[2] = xmm0;			/* Save R2 */ \
	sd0[3] = xmm1;			/* Save I2 */ \
	srcreg += srcinc; \
}

#if 0 //DEF X86_64;
\
/*r4_x2cl_four_complex_first_djbfft MACRO srcreg,srcinc,d1,screg,pmreg */ \
/*	r4_x2cl_four_complex_first_djbfft_cmn srcreg,rbx,srcinc,d1,screg,pmreg */ \
/*	ENDM */ \

#define r4_x2cl_four_complex_first_djbfft_scratch(srcreg,srcinc,d1,screg,pmreg) \
	r4_x2cl_four_complex_first_djbfft_cmn(srcreg,0,srcinc,d1,screg,pmreg); \
\

#define r4_x2cl_four_complex_first_djbfft_cmn(srcreg,off,srcinc,d1,screg,pmreg) \
\
	xmm0 = xptr(srcreg+off+d1);		/* R2 */ \
	xmm1 = xptr(pmreg+32+16);		/* cosine/sine 2 */ \
	xmm2 = xmm0;			/* Copy R2			; 1-3 */ \
	xmm0 *= xmm1;			/* A2 = R2 * cosine/sine	; 1-5 */ \
\
	xmm3 = xptr(srcreg+off+d1+16);	/* R4 */ \
	xmm4 = xptr(pmreg+96+16);		/* cosine/sine 4 */ \
	xmm5 = xmm3;			/* Copy R4			; 2-4 */ \
	xmm3 *= xmm4;			/* A4 = R4 * cosine/sine	; 2-6 */ \
\
	xmm6 = xptr(srcreg+off+d1+32);	/* I2 */ \
	xmm1 *= xmm6;			/* B2 = I2 * cosine/sine	; 3-7 */ \
\
	xmm7 = xptr(srcreg+off+d1+48);	/* I4 */ \
	xmm4 *= xmm7;			/* B4 = I4 * cosine/sine	; 4-8 */ \
\
	xmm8 = xptr(srcreg+off);		/* R1 */ \
	xmm9 = xptr(pmreg+0+16);		/* cosine/sine 1 */ \
	xmm10 = xmm8;			/* Copy R1			; 5-7 */ \
	xmm8 *= xmm9;			/* A1 = R1 * cosine/sine	; 5-9 */ \
\
	xmm0 -= xmm6;			/* A2 = A2 - I2			; 6-8 */ \
	xmm11 = xptr(srcreg+off+16);		/* R3 */ \
	xmm12 = xptr(pmreg+64+16);		/* cosine/sine 3 */ \
	xmm13 = xmm11;			/* Copy R3			; 6-8 */ \
	xmm11 *= xmm12;			/* A3 = R3 * cosine/sine	; 6-10 */ \
\
	xmm3 -= xmm7;			/* A4 = A4 - I4			; 7-9 */ \
	xmm14 = xptr(srcreg+off+32);		/* I1 */ \
	xmm9 *= xmm14;			/* B1 = I1 * cosine/sine	; 7-11 */ \
\
	xmm1 += xmm2;			/* B2 = B2 + R2			; 8-10 */ \
	xmm15 = xptr(srcreg+off+48);		/* I3 */ \
	xmm12 *= xmm15;			/* B3 = I3 * cosine/sine	; 8-12	avail 6,7,2 */ \
\
	xmm4 += xmm5;			/* B4 = B4 + R4			; 9-11	avail 6,7,2,5 */ \
	xmm6 = xptr(pmreg+32);		/* sine 2			; 9	avail 7,2,5 */ \
	xmm0 *= xmm6;			/* A2 = A2 * sine		; 9-13 */ \
\
	xmm8 -= xmm14;			/* A1 = A1 - I1			; 10-12	avail 7,2,5,14 */ \
	xmm7 = xptr(pmreg+96);		/* sine 4			; 10	avail 2,5,14 */ \
	xmm3 *= xmm7;			/* A4 = A4 * sine		; 10-14 */ \
\
	xmm11 -= xmm15;			/* A3 = A3 - I3			; 11-13	avail 2,5,14,15 */ \
	xmm1 *= xmm6;			/* B2 = B2 * sine		; 11-15	avail 2,5,14,15,6 */ \
	xmm2 = xptr(pmreg+0);			/* sine 1			; 11	avail 5,14,15,6 */ \
\
	xmm9 += xmm10;			/* B1 = B1 + R1			; 12-14	avail 5,14,15,6,10 */ \
	xmm4 *= xmm7;			/* B4 = B4 * sine		; 12-16	avail 5,14,15,6,10,7 */ \
	xmm5 = xptr(pmreg+64);		/* sine 3			; 12	avail 14,15,6,10,7 */ \
\
	xmm12 += xmm13;			/* B3 = B3 + R3			; 13-15	avail 14,15,6,10,7,13 */ \
	xmm8 *= xmm2;			/* A1 = A1 * sine		; 13-17 */ \
\
	xmm11 *= xmm5;			/* A3 = A3 * sine		; 14-18 */ \
	xmm6 = xmm0;			/* Copy R2			; 14 on Core i7, 14-16 on Core 2 avail 14,15,10,7,13 */ \
\
	xprefetchw(u8ptr(srcreg+srcinc)); \
\
	xmm0 -= xmm3;			/* R2 = R2 - R4 (new R4)	; 15-17 */ \
	xmm9 *= xmm2;			/* B1 = B1 * sine		; 15-19	avail 14,15,10,7,13,2 */ \
\
	xmm3 += xmm6;			/* R4 = R2 + R4 (new R2)	; 16-18	avail 14,15,10,7,13,2,6 */ \
	xmm12 *= xmm5;			/* B3 = B3 * sine		; 16-20	avail 14,15,10,7,13,2,6,5 */ \
\
	xmm7 = xmm1;			/* Copy I2			; 16 on Core i7, 16-18 on Core 2 */ \
	xmm1 -= xmm4;			/* I2 = I2 - I4 (new I4)	; 17-19 */ \
	xmm4 += xmm7;			/* I4 = I2 + I4 (new I2)	; 18-20 */ \
\
	xmm7 = xmm8;			/* Copy R1			; 18 on Core i7, 16-18 on Core 2 */ \
	xmm8 -= xmm11;			/* R1 = R1 - R3 (new R3)	; 19-21 */ \
	xmm11 += xmm7;			/* R3 = R1 + R3 (new R1)	; 20-22 */ \
\
	xmm7 = xmm9;			/* Copy I1			; 20 on Core i7, 16-18 on Core 2 */ \
\
	xmm9 -= xmm12;			/* I1 = I1 - I3 (new I3)	; 21-23 */ \
	xmm2 = xmm1;			/* Copy I4			; 21-23 (20) avail 14,15,10,13,6,5 */ \
\
	xprefetchw(u8ptr(srcreg+srcinc+d1)); \
\
	xmm12 += xmm7;			/* I3 = I1 + I3 (new I1)	; 22-24 avail 14,15,10,13,6,5,7 */ \
	xmm6 = xptr(screg+0+16);		/* cosine/sine 3&4		; 22	avail 14,15,10,13,5,7 */ \
\
	xmm1 += xmm8;			/* I4 = R3 + I4 (newer R4)	; 23-25 */ \
	xmm7 = xmm0;			/* Copy R4			; 23-25 (18) avail 14,15,10,13,5 */ \
\
	xmm8 -= xmm2;			/* R3 = R3 - I4 (newer R3)	; 24-26 avail 14,15,10,13,5,2 */ \
	xmm2 = xmm11;			/* Copy R1			; 24-26 (23) avail 14,15,10,13,5 */ \
\
	xmm0 += xmm9;			/* R4 = I3 + R4 (newer I3)	; 25-27 */ \
	xmm10 = xmm12;			/* Copy I1			; 27-29 (25) avail 15,13,5 */ \
\
	xmm9 -= xmm7;			/* I3 = I3 - R4 (newer I4)	; 26-28 avail 14,15,13,5,7 */ \
	xmm14 = xmm1;			/* Copy R4			; 26-28 avail 15,13,5,7 */ \
	xmm1 *= xmm6;			/* A4 = R4 * cosine/sine	; 26-30 */ \
\
	xmm11 -= xmm3;			/* R1 = R1 - R2 (newer R2)	; 27-29 */ \
	xmm15 = xmm8;			/* Copy R3			; 27-29 avail 13,5,7 */ \
	xmm8 *= xmm6;			/* A3 = R3 * cosine/sine	; 27-31 */ \
\
	xmm12 -= xmm4;			/* I1 = I1 - I2 (newer I2)	; 28-30 */ \
	xmm13 = xmm0;			/* Copy I3			; 28-30 avail 5,7 */ \
	xmm0 *= xmm6;			/* B3 = I3 * cosine/sine	; 28-32 */ \
\
	xmm3 += xmm2;			/* R2 = R1 + R2 (final R1)	; 29-31	avail 5,7,2 */ \
	xmm6 *= xmm9;			/* B4 = I4 * cosine/sine	; 29-33 */ \
	xmm5 = xptr(screg+32+16);		/* cosine/sine 2		; 29	avail 7,2 */ \
\
	xmm4 += xmm10;			/* I2 = I1 + I2 (final I1)	; 30-32	avail 7,2,10 */ \
	xmm7 = xmm11;			/* Copy R2			; 30-32	avail 2,10 */ \
	xmm11 *= xmm5;			/* A2 = R2 * cosine/sine	; 30-34 */ \
\
	xmm1 += xmm9;			/* A4 = A4 + I4			; 31-33	avail 2,10,9 */ \
	xmm5 *= xmm12;			/* B2 = I2 * cosine/sine	; 31-35 */ \
	xmm2 = xptr(screg+0);			/* sine 3&4			; 31	avail 10,9 */ \
\
	xmm8 -= xmm13;			/* A3 = A3 - I3			; 32-34	avail 10,9,13 */ \
	xptr(srcreg) = xmm3;			/* Save R1			; 32	avail 10,9,13,3 */ \
\
	xmm0 += xmm15;			/* B3 = B3 + R3			; 33-35	avail 10,9,13,3,15 */ \
	xptr(srcreg+16) = xmm4;		/* Save I1			; 33	avail 10,9,13,3,15,4 */ \
\
	xmm6 -= xmm14;			/* B4 = B4 - R4			; 34-36	avail 10,9,13,3,15,4,14 */ \
	xmm1 *= xmm2;			/* A4 = A4 * sine (final R4)	; 34-38 */ \
	xmm3 = xptr(screg+32);		/* sine 2			; 34	avail 10,9,13,15,4,14 */ \
\
	xmm11 -= xmm12;			/* A2 = A2 - I2			; 35-37	avail 10,9,13,15,4,14,12 */ \
	xmm8 *= xmm2;			/* A3 = A3 * sine (final R3)	; 35-39 */ \
\
	xmm5 += xmm7;			/* B2 = B2 + R2			; 36-38	avail 10,9,13,15,4,14,12,7 */ \
	xmm0 *= xmm2;			/* B3 = B3 * sine (final I3)	; 36-40 */ \
\
	xmm6 *= xmm2;			/* B4 = B4 * sine (final I4)	; 37-41	avail 10,9,13,15,4,14,12,7,2 */ \
\
	xmm11 *= xmm3;			/* A2 = A2 * sine (final R2)	; 38-42 */ \
\
	xmm5 *= xmm3;			/* B2 = B2 * sine (final I2)	; 39-43	avail 10,9,13,15,4,14,12,7,2,3 */ \
\
	xptr(srcreg+d1+32) = xmm1;		/* Save R4			; 39 */ \
	xptr(srcreg+d1) = xmm8;		/* Save R3			; 40 */ \
	xptr(srcreg+d1+16) = xmm0;		/* Save I3			; 41 */ \
	xptr(srcreg+d1+48) = xmm6;		/* Save I4			; 42 */ \
	xptr(srcreg+32) = xmm11;		/* Save R2			; 43 */ \
	xptr(srcreg+48) = xmm5;		/* Save I2			; 44 */ \
	srcreg += srcinc; \
\

#endif
\
/* */ \
/* ************************************* four-complex-last-unfft variants ****************************************** */ \
/* */ \
\
/* This code applies the sin/cos multipliers before a radix-4 butterfly. */ \
/* Then it applies the premultipliers since the negacyclic inverse FFT is complete. */ \
/* We save memory by splitting the roots-of-minus-1 premultipliers such that every */ \
/* macro uses the same premultiplier data.  Every sin/cos postmultiplier, a very big */ \
/* table anyway, are then multiplied by the other part of the split roots-of-minus-1. */ \
/* The common postmuliplier data is 1, .924-.383i, SQRTHALF-SQRTHALFi, .383-.924i. */ \
\

#define r4_x4cl_four_complex_last_unfft4(srcreg,srcinc,d1,d2,screg,off) { \
	vec2f64 *sd0 = (vec2f64*)(srcreg); \
	vec2f64 *sd1 = (vec2f64*)(srcreg+d1); \
	vec2f64 *sd2 = (vec2f64*)(srcreg+d2); \
	vec2f64 *sd3 = (vec2f64*)(srcreg+d2+d1); \
	vec2f64 xma0,xma1,xma2,xma3,xma4,xma5,xma6,xma7; \
	xmm0 = sd0[0];	/* R1 */ \
	xmm1 = sd0[2];	/* R2 */ \
	xmm4 = sd2[0];	/* R5 */ \
	xmm5 = sd2[2];	/* R6 */ \
	r4_x4c_unfft_postmult(xmm0,xmm1,xmm2,xmm3,xmm4,xmm5,xmm6,xmm7,sd1[0],sd1[2],sd3[0],sd3[2],sd0[0],sd0[2],screg,0,srcreg+srcinc,d1); \
/*	xstore	[srcreg], xmm		;; Save R1 */ \
/*	xstore	[srcreg+32], xmm	;; Save R5 */ \
	xma0 = sd0[1];	/* R1 */ \
	xma1 = sd0[3];	/* R2 */ \
	sd0[1] = xmm2;	/* Save R3 */ \
	sd0[3] = xmm3;	/* Save R7 */ \
	xma4 = sd2[1];	/* R5 */ \
	xma5 = sd2[3];	/* R6 */ \
	sd2[0] = xmm7;	/* Save R2 */ \
	sd2[1] = xmm4;	/* Save R4 */ \
	sd2[2] = xmm0;	/* Save R6 */ \
	sd2[3] = xmm1;	/* Save R8 */ \
	r4_x4c_unfft_postmult(xma0,xma1,xma2,xma3,xma4,xma5,xma6,xma7,sd1[1],sd1[3],sd3[1],sd3[3],sd1[0],sd1[2],screg,off,srcreg+srcinc+d2,d1); \
/*	xstore	[srcreg+d1], xmb	;; Save R1 */ \
/*	xstore	[srcreg+d1+32], xmb	;; Save R5 */ \
	sd1[1] = xma2;	/* Save R3 */ \
	sd1[3] = xma3;	/* Save R7 */ \
	sd3[0] = xma7;	/* Save R2 */ \
	sd3[1] = xma4;	/* Save R4 */ \
	sd3[2] = xma0;	/* Save R6 */ \
	sd3[3] = xma1;	/* Save R8 */ \
	srcreg += srcinc; \
}


#define r4_x4c_unfft_postmult(r1, r2, r3, r4, r5, r6, r7, r8, mem3, mem4, mem7, mem8, dest1, dest2, screg, off, pre1, pre2) { \
	vec2f64 s1, m1, temp; \
	vec2f64 *sc = (vec2f64*)(screg+off); \
\
	s1 = sc[1];			/* cosine/sine */ \
	r7 = s1 * r1 + r2;	/* A1 = R1 * cosine/sine + I1 */ \
	r2 = r2 * s1 - r1;	/* B1 = I1 * cosine/sine - R1 */ \
\
	s1 = r8 = sc[3];	/* cosine/sine */ \
	r8 = s1 * r5 + r6;	/* A3 = R3 * cosine/sine + I3 */ \
	r6 = r6 * s1 - r5;	/* B3 = I3 * cosine/sine - R3 */ \
\
	m1 = mem3; \
	r4 = mem4;			/* R4 (I2) */ \
	s1 = sc[5];			/* cosine/sine */ \
	r3 = s1 * m1 + r4;	/* A2 = R2 * cosine/sine + I2 */ \
	r4 = r4 * s1 - m1;	/* B2 = I2 * cosine/sine - R2 */ \
\
	m1 = mem7; \
	r5 = mem8;			/* I4 */ \
	s1 = sc[7];			/* cosine/sine */ \
	r1 = s1 * m1 + r5;	/* A4 = R4 * cosine/sine + I4 */ \
	r5 = r5 * s1 - m1;	/* B4 = I4 * cosine/sine - R4 */ \
\
	s1 = sc[0]; \
	r7 *= s1;			/* A1 = A1 * sine (new R1) */ \
	r2 *= s1;			/* B1 = B1 * sine (new I1) */ \
\
	s1 = sc[2]; \
	r8 *= s1;			/* A3 = A3 * sine (new R3) */ \
	r6 *= s1;			/* B3 = B3 * sine (new I3) */ \
\
	s1 = sc[4]; \
	r3 *= s1;			/* A2 = A2 * sine (new R2) */ \
	r4 *= s1;			/* B2 = B2 * sine (new I2) */ \
\
	s1 = sc[6]; \
	r1 *= s1;			/* A4 = A4 * sine (new R4) */ \
	r5 *= s1;			/* B4 = B4 * sine (new I4) */ \
\
	xprefetchw(u8ptr(pre1)); \
\
	temp = r7; \
	r7 -= r3;			/* R1 = R1 - R2 (new R2) */ \
	r3 += temp;			/* R2 = R1 + R2 (new R1) */ \
\
	temp = r2; \
	r2 -= r4;			/* I1 = I1 - I2 (new I2) */ \
	r4 += temp;			/* I2 = I1 + I2 (new I1) */ \
\
	temp = r1; \
	r1 -= r8;			/* R4 = R4 - R3 (new I4) */ \
	r8 += temp;			/* R3 = R4 + R3 (new R3) */ \
\
	temp = r6; \
	r6 -= r5;			/* I3 = I3 - I4 (new R4) */ \
	r5 += temp;			/* I4 = I3 + I4 (new I3) */ \
\
	xprefetchw(u8ptr(pre1+pre2)); \
\
	temp = r7; \
	r7 -= r6;			/* R2 = R2 - R4 (new R4) */ \
	r6 += temp;			/* R4 = R2 + R4 (new R2) */ \
\
	temp = r2; \
	r2 -= r1;			/* I2 = I2 - I4 (new I4) */ \
	r1 += temp;			/* I4 = I2 + I4 (new I2) */ \
\
	temp = r3; \
	r3 -= r8;			/* R1 = R1 - R3 (new R3) */ \
	r8 += temp;			/* R3 = R1 + R3 (final R1) */ \
\
	temp = r4; \
	r4 -= r5;			/* I1 = I1 - I3 (new I3) */ \
	r5 += temp;			/* I3 = I1 + I3 (final I1) */ \
\
	dest1 = r8;			/* Save R1 */ \
	dest2 = r5;			/* Save I1 */ \
\
	r3 *= XMM_SQRTHALF2;	/* A3 = R3 * SQRTHALF */ \
	r4 *= XMM_SQRTHALF2;	/* B3 = I3 * SQRTHALF */ \
	r8 = r4;				/* Copy B3 */ \
	r4 -= r3;				/* B3 = B3 - A3 (final I3) */ \
	r3 += r8;				/* A3 = A3 + B3	(final R3) */ \
\
	r8 = r6 * CONST2_P924;	/* A2 = R2 * .924 */ \
	r5 = r7 * CONST2_P383;	/* A4 = R4 * .383 */ \
	r8 += r1 * CONST2_P383;	/* A2 = A2 + I2 * .383 (final R2) */ \
	r5 += r2 * CONST2_P924;	/* A4 = A4 + I4 * .924	(final R4) */ \
	r1 *= CONST2_P924;		/* D2 = I2 * .924 */ \
	r2 *= CONST2_P383;		/* D4 = I4 * .383 */ \
	r1 -= r6 * CONST2_P383;	/* D2 = D2 - R2 * .383 (final I2) */ \
	r2 -= r7 * CONST2_P924;	/* D4 = D4 - R4 * .924	(final I4) */ \
}

/* This code applies the sin/cos multipliers before a radix-4 butterfly. */ \
/* After the butterfly, it applies the negacyclic premultipliers since the inverse FFT is complete. */ \
\

#define r4_x4cl_four_complex_last_djbunfft(srcreg,srcinc,d1,d2,screg,off,pmreg,pmoff) { \
	vec2f64 *sd0 = (vec2f64*)(srcreg); \
	vec2f64 *sd1 = (vec2f64*)(srcreg+d1); \
	vec2f64 *sd2 = (vec2f64*)(srcreg+d2); \
	vec2f64 *sd3 = (vec2f64*)(srcreg+d2+d1); \
	vec2f64 xma0,xma1,xma2,xma3,xma4,xma5,xma6,xma7; \
	xmm2 = sd1[1];	/* R2 */ \
	xmm6 = sd3[1];	/* R4 */ \
	xmm3 = sd1[3];	/* I2 */ \
	xmm7 = sd3[3];	/* I4 */ \
	r4_x4c_djbunfft_postmult(xmm0,xmm1,xmm2,xmm3,xmm4,xmm5,xmm6,xmm7, sd0[1],sd0[3],sd2[1],sd2[3],sd1[1],sd1[3], screg,off,pmreg,pmoff,srcreg+srcinc+d2,d1); \
	xma2 = sd1[0];	/* R2 */ \
	sd1[0] = xmm5;	/* Save R1 */ \
	sd1[1] = xmm2;	/* Save R3 */ \
	xma3 = sd1[2];	/* I2 */ \
	sd1[2] = xmm7;	/* Save I1 */ \
	sd1[3] = xmm3;	/* Save I3 */ \
	xma6 = sd3[0];	/* R4 */ \
	xma7 = sd3[2];	/* I4 */ \
	sd3[0] = xmm1;	/* Save R2 */ \
	sd3[1] = xmm4;	/* Save R4 */ \
	sd3[2] = xmm6;	/* Save I2 */ \
	sd3[3] = xmm0;	/* Save I4 */ \
	r4_x4c_djbunfft_postmult(xma0,xma1,xma2,xma3,xma4,xma5,xma6,xma7, sd0[0],sd0[2],sd2[0],sd2[2],sd0[0],sd0[2], screg,0,pmreg,0,srcreg+srcinc,d1); \
	sd0[0] = xma5;	/* Save R1 */ \
	sd0[1] = xma2;	/* Save R3 */ \
	sd0[2] = xma7;	/* Save I1 */ \
	sd0[3] = xma3;	/* Save I3 */ \
	sd2[0] = xma1;	/* Save R2 */ \
	sd2[1] = xma4;	/* Save R4 */ \
	sd2[2] = xma6;	/* Save I2 */ \
	sd2[3] = xma0;	/* Save I4 */ \
	srcreg += srcinc; \
}


#define r4_x4c_djbunfft_postmult(r1,r2,r3,r4,r5,r6,r7,r8,mem1,mem2,mem5,mem6,dest1,dest2,screg,off,pmreg,pmoff,pre1,pre2) { \
	vec2f64 s0,s1,s2,s3; \
	vec2f64 *sc = (vec2f64*)(screg+off); \
	vec2f64 *pm = (vec2f64*)(pmreg+pmoff); \
	\
	s2 = r5 = sc[3];		/* cosine/sine */ \
	r5 = r3 * s2 + r4;		/* A2 = R2 * cosine/sine + I2 */ \
	r4 = r4 * s2 - r3;		/* B2 = I2 * cosine/sine - R2 */ \
\
	r6 = r7;				/* Copy R4 */ \
	r1 = sc[1];				/* cosine/sine */ \
	r7 = r7 * r1 - r8;		/* A4 = R4 * cosine/sine - I4 */ \
	r8 = r8 * r1 + r6;		/* B4 = I4 * cosine/sine + R4 */ \
\
	s0 = mem5; \
	s1 = mem6; \
	r2 = s0 * r1 + s1;		/* A3 = R3 * cosine/sine + I3 */ \
	r1 = r1 * s1 - s0;		/* B3 = I3 * cosine/sine - R3 */ \
\
	s2 = sc[2]; \
	r5 *= s2;				/* A2 = A2 * sine (new R2) */ \
	r4 *= s2;				/* B2 = B2 * sine (new I2) */ \
\
	r3 = sc[0]; \
	r7 *= r3;				/* A4 = A4 * sine (new R4) */ \
	r8 *= r3;				/* B4 = B4 * sine (new I4) */ \
\
	r2 *= r3;				/* A3 = A3 * sine (new R3) */ \
	r1 *= r3;				/* B3 = B3 * sine (new I3) */ \
\
	xprefetchw(u8ptr(pre1)); \
\
	r6 = mem1;				/* R1 */ \
	r6 -= r5;				/* R1 = R1 - R2 (new R2) */ \
	r5 += mem1;				/* R2 = R1 + R2 (new R1) */ \
	r3 = r7;				/* Copy R4 */ \
	r7 -= r2;				/* R4 = R4 - R3 (new I4) */ \
	r2 += r3;				/* R3 = R4 + R3 (new R3) */ \
	r3 = r1;				/* Copy I3 */ \
	r1 -= r8;				/* I3 = I3 - I4 (new R4) */ \
	r8 += r3;				/* I4 = I3 + I4 (new I3) */ \
\
	xprefetchw(u8ptr(pre1+pre2)); \
\
	r3 = r6;				/* Copy R2 */ \
	r6 -= r1;				/* R2 = R2 - R4 (new R4) */ \
	r1 += r3;				/* R4 = R2 + R4 (new R2) */ \
\
	r3 = mem2;				/* I1 */ \
	r3 -= r4;				/* I1 = I1 - I2 (new I2) */ \
	r4 += mem2;				/* I2 = I1 + I2 (new I1) */ \
\
	dest1 = r6;				/* Save R4 */ \
\
	r6 = r5;				/* Copy R1 */ \
	r5 -= r2;				/* R1 = R1 - R3 (new R3) */ \
	r2 += r6;				/* R3 = R1 + R3 (new R1) */ \
\
	r6 = r7;				/* Copy I4 */ \
	r7 += r3;				/* I4 = I2 + I4 (new I2) */ \
	r3 -= r6;				/* I2 = I2 - I4 (new I4) */ \
\
	r6 = r4;				/* Copy I1 */ \
	r4 -= r8;				/* I1 = I1 - I3 (new I3) */ \
	r8 += r6;				/* I3 = I1 + I3 (new I1) */ \
\
	dest2 = r3;				/* Save I4 */ \
\
	s0 = pm[1]; \
	s2 = pm[5]; \
	s1 = pm[3]; \
	s3 = pm[7]; \
\
	r6 = s0 * r2 + r8;		/* A1 = R1 * cosine/sine + I1 */ \
	r8 = r8 * s0 - r2;		/* B1 = I1 * cosine/sine - R1 */ \
	r3 = s2 * r5 + r4;		/* A3 = R3 * cosine/sine + I3 */ \
	r4 = r4 * s2 - r5;		/* B3 = I3 * cosine/sine - R3 */ \
	r2 = s1 * r1 + r7;		/* A2 = R2 * cosine/sine + I2 */ \
	r7 = r7 * s1 - r1;		/* B2 = I2 * cosine/sine - R2 */ \
	r5 = dest1 * s3 + dest2;	/* A4 = R4 * cosine/sine + I4 */ \
	r1 = dest2 * s3 - dest1;	/* B4 = I4 * cosine/sine - R4 */ \
\
	s0 = pm[0]; \
	s2 = pm[4]; \
	s1 = pm[2]; \
	s3 = pm[6]; \
	r6 *= s0;				/* A1 = A1 * sine */ \
	r8 *= s0;				/* B1 = B1 * sine */ \
	r3 *= s2;				/* A3 = A3 * sine */ \
	r4 *= s2;				/* B3 = B3 * sine */ \
	r2 *= s1;				/* A2 = A2 * sine */ \
	r7 *= s1;				/* B2 = B2 * sine */ \
	r5 *= s3;				/* A4 = A4 * sine */ \
	r1 *= s3;				/* B4 = B4 * sine */ \
}

#if 0 //DEF X86_64;
\

#define r4_x4cl_four_complex_last_djbunfft(srcreg,srcinc,d1,d2,screg,scoff,pmreg,pmoff) \
\
	xmm0 = xptr(srcreg+d1+16);		/* R2 */ \
	xmm1 = xptr(screg+scoff+32+16);	/* cosine/sine 2 */ \
	xmm2 = xmm0;			/* Copy R2 */ \
	xmm0 *= xmm1;			/* A2 = R2 * cosine/sine	; 1-5 */ \
	xmm3 = xptr(srcreg+d1+48);		/* I2 */ \
	xmm1 *= xmm3;			/* B2 = I2 * cosine/sine	; 2-6 */ \
\
	xmm4 = xptr(srcreg+d2+d1+16);		/* R4 */ \
	xmm5 = xptr(screg+scoff+0+16);	/* cosine/sine 3&4 */ \
	xmm6 = xmm4;			/* Copy R4 */ \
	xmm4 *= xmm5;			/* A4 = R4 * cosine/sine	; 3-7 */ \
	xmm7 = xptr(srcreg+d2+d1+48);		/* I4 */ \
	xmm8 = xmm7;			/* Copy I4 */ \
	xmm7 *= xmm5;			/* B4 = I4 * cosine/sine	; 4-8 */ \
\
	xmm9 = xptr(srcreg+d2+16);		/* R3 */ \
	xmm10 = xmm9;			/* Copy R3 */ \
	xmm9 *= xmm5;			/* A3 = R3 * cosine/sine	; 5-7 */ \
\
	xmm0 += xmm3;			/* A2 = A2 + I2			; 6-8	avail 3,11+ */ \
	xmm11 = xptr(srcreg+d2+48);		/* I3				; 6	avail 3,12+ */ \
	xmm5 *= xmm11;			/* B3 = I3 * cosine/sine	; 6-10 */ \
\
	xmm1 -= xmm2;			/* B2 = B2 - R2			; 7-9	avail 2,3,12+ */ \
	xmm12 = xptr(srcreg+d1);		/*#2 R2 */ \
	xmm13 = xptr(screg+32+16);		/*#2 cosine/sine 2 */ \
	xmm14 = xmm12;			/*#2 Copy R2 */ \
	xmm12 *= xmm13;			/*#2 A2 = R2 * cosine/sine	; 7-11	avail 2,3,15 */ \
\
	xmm4 -= xmm8;			/* A4 = A4 - I4			; 8-10	avail 2,3,8,15 */ \
	xmm15 = xptr(srcreg+d1+32);		/*#2 I2				; 8	avail 2,3,8 */ \
	xmm13 *= xmm15;			/*#2 B2 = I2 * cosine/sine	; 8-12 */ \
\
	xmm7 += xmm6;			/* B4 = B4 + R4			; 9-11	avail 2,3,8,6 */ \
	xmm2 = xptr(screg+scoff+32);		/* sine 2			; 9	avail 3,8,6 */ \
	xmm0 *= xmm2;			/* A2 = A2 * sine (R2)		; 9-13 */ \
\
	xmm9 += xmm11;			/* A3 = A3 + I3			; 10-12	avail 3,8,6,11 */ \
	xmm1 *= xmm2;			/* B2 = B2 * sine (I2)		; 10-14	avail 3,8,6,11,2 */ \
	xmm3 = xptr(screg+scoff+0);		/* sine 3&4			; 10	avail 8,6,11,2 */ \
\
	xmm5 -= xmm10;			/* B3 = B3 - R3			; 11-13	avail 8,6,11,2,10 */ \
	xmm6 = xptr(srcreg+16);		/* R1				; 11	avail 8,11,2,10 */ \
	xmm4 *= xmm3;			/* A4 = A4 * sine (R4)		; 11-15 */ \
\
	xmm12 += xmm15;			/*#2 A2 = A2 + I2		; 12-14	avail 8,11,2,10,15 */ \
	xmm7 *= xmm3;			/* B4 = B4 * sine (I4)		; 12-16 */ \
	xmm2 = xmm6;			/* Copy R1			; 12-14	avail 8,11,10,15 */ \
	xmm8 = xptr(screg+32);		/*#2 sine 2			; 12	avail 11,10,15 */ \
\
	xmm13 -= xmm14;			/*#2 B2 = B2 - R2		; 13-15	avail 11,10,15,14 */ \
	xmm9 *= xmm3;			/* A3 = A3 * sine (R3)		; 13-17 */ \
	xmm11 = xptr(srcreg+48);		/* I1				; 13	avail 10,15,14 */ \
\
	xmm6 -= xmm0;			/* R1 = R1 - R2 (new R2)	; 14-16 */ \
	xmm5 *= xmm3;			/* B3 = B3 * sine (I3)		; 14-18	avail 10,15,14,3 */ \
	xmm3 = xmm11;			/* Copy I1			; 14-16	avail 10,15,14 */ \
\
	xmm0 += xmm2;			/* R2 = R1 + R2 (new R1)	; 15-17	avail 10,15,14,2 */ \
	xmm12 *= xmm8;			/*#2 A2 = A2 * sine (R2)	; 15-19 (can safely delay) */ \
\
	xmm11 -= xmm1;			/* I1 = I1 - I2 (new I2)	; 16-18 */ \
	xmm13 *= xmm8;			/*#2 B2 = B2 * sine (I2)	; 16-20	(can safely delay) avail 10,15,14,2,8 */ \
	xmm2 = xmm4;			/* Copy R4			; 16-18	avail 10,15,14,8 */ \
\
	xprefetchw(u8ptr(srcreg+srcinc+d2)); \
\
	xmm1 += xmm3;			/* I2 = I1 + I2 (new I1)	; 17-19	avail 10,15,14,8,3 */ \
\
	xmm4 -= xmm9;			/* R4 = R4 - R3 (new I4)	; 18-20 */ \
	xmm3 = xmm7;			/* Copy I4			; 18-20 (17) avail 10,15,14,8 */ \
\
	xmm9 += xmm2;			/* R3 = R4 + R3 (new R3)	; 19-21	avail 10,15,14,8,2 */ \
\
	xmm7 += xmm5;			/* I4 = I3 + I4 (new I3)	; 20-22 */ \
	xmm2 = xmm11;			/* Copy I2			; 20-22 (19) avail 10,15,14,8 */ \
\
	xmm5 -= xmm3;			/* I3 = I3 - I4 (new R4)	; 21-23	avail 10,15,14,8,3 */ \
	xmm10 = xptr(pmreg+pmoff+96+16);	/* cosine/sine 4		; 21	avail 15,14,8,3 */ \
\
	xprefetchw(u8ptr(srcreg+srcinc+d2+d1)); \
\
	xmm11 -= xmm4;			/* I2 = I2 - I4 (newer I4)	; 22-24 */ \
	xmm3 = xmm6;			/* Copy R2			; 22-24 (17) avail 15,14,8 */ \
	xmm15 = xptr(pmreg+pmoff+32+16);	/* cosine/sine 2		; 22	avail 14,8 */ \
\
	xmm4 += xmm2;			/* I4 = I2 + I4 (newer I2)	; 23-25	avail 14,8,2 */ \
	xmm2 = xmm0;			/* Copy R1			; 23-25 (18) avail 14,8 */ \
\
	xmm6 -= xmm5;			/* R2 = R2 - R4 (newer R4)	; 24-26 */ \
	xmm8 = xmm1;			/* Copy I1			; 24-26 (20) avail 14 */ \
\
	xmm5 += xmm3;			/* R4 = R2 + R4 (newer R2)	; 25-27	avail 14,3 */ \
	xmm14 = xmm11;			/* Copy I4			; 25-27	avail 3 */ \
	xmm11 *= xmm10;			/* B4 = I4 * cosine/sine	; 25-29 */ \
\
	xmm0 -= xmm9;			/* R1 = R1 - R3 (newer R3)	; 26-28 */ \
	xmm3 = xmm4;			/* Copy I2			; 26-28	avail none */ \
	xmm4 *= xmm15;			/* B2 = I2 * cosine/sine	; 26-30 */ \
\
	xmm1 -= xmm7;			/* I1 = I1 - I3 (newer I3)	; 27-29 */ \
	xmm10 *= xmm6;			/* A4 = R4 * cosine/sine	; 27-31 */ \
\
	xmm9 += xmm2;			/* R3 = R1 + R3 (newer R1)	; 28-30	avail 2 */ \
	xmm15 *= xmm5;			/* A2 = R2 * cosine/sine	; 28-32 */ \
	xmm2 = xptr(pmreg+pmoff+64+16);	/* cosine/sine 3		; 28	avail none */ \
\
	xmm7 += xmm8;			/* I3 = I1 + I3 (newer I1)	; 29-31	avail 8 */ \
	xmm8 = xmm0;			/* Copy R3			; 29-31	avail none */ \
	xmm0 *= xmm2;			/* A3 = R3 * cosine/sine	; 29-33 */ \
\
	xmm11 -= xmm6;			/* B4 = B4 - R4			; 30-32	avail 6 */ \
	xmm2 *= xmm1;			/* B3 = I3 * cosine/sine	; 30-34 */ \
	xmm6 = xptr(pmreg+pmoff+0+16);	/* cosine/sine 1		; 30	avail none */ \
\
	xmm4 -= xmm5;			/* B2 = B2 - R2			; 31-33	avail 5 */ \
	xmm5 = xmm9;			/* Copy R1			; 31-33	avail none */ \
	xmm9 *= xmm6;			/* A1 = R1 * cosine/sine	; 31-35 */ \
\
	xmm10 += xmm14;			/* A4 = A4 + I4			; 32-34	avail 14 */ \
	xmm6 *= xmm7;			/* B1 = I1 * cosine/sine	; 32-36 */ \
	xmm14 = xptr(pmreg+pmoff+96);		/* sine 4			; 32	avail none */ \
\
	xmm15 += xmm3;			/* A2 = A2 + I2			; 33-35	avail 3 */ \
	xmm11 *= xmm14;			/* B4 = B4 * sine (final I4)	; 33-37	avail 3 storable 11 */ \
	xmm3 = xptr(pmreg+pmoff+32);		/* sine 2			; 33	avail none storable 11 */ \
\
	xmm0 += xmm1;			/* A3 = A3 + I3			; 34-36	avail 1 storable 11 */ \
	xmm4 *= xmm3;			/* B2 = B2 * sine (final I2)	; 34-38	avail 1 storable 11,4 */ \
	xmm1 = xptr(srcreg+d2+d1);		/*#2 R4				; 34	avail none storable 11,4 */ \
\
	xmm2 -= xmm8;			/* B3 = B3 - R3			; 35-37	avail 8 storable 11,4 */ \
	xmm10 *= xmm14;			/* A4 = A4 * sine (final R4)	; 35-39	avail 8,14 storable 11,4,10 */ \
	xmm8 = xptr(screg+0+16);		/*#2 cosine/sine 3&4		; 35	avail 14 storable 11,4,10 */ \
\
	xmm9 += xmm7;			/* A1 = A1 + I1			; 36-38	avail 14,7 storable 11,4,10 */ \
	xmm15 *= xmm3;			/* A2 = A2 * sine (final R2)	; 36-40	avail 14,7,3 storable 11,4,10,15 */ \
	xmm7 = xptr(srcreg+d2);		/*#2 R3				; 36	avail 14,3 storable 11,4,10,15 */ \
\
	xmm6 -= xmm5;			/* B1 = B1 - R1			; 37-39	avail 14,3,5 storable 11,4,10,15 */ \
	xmm3 = xmm1;			/*#2 Copy R4			; 37-39	avail 14,5 storable 11,4,10,15 */ \
	xmm1 *= xmm8;			/*#2 A4 = R4 * cosine/sine	; 37-41 */ \
	xmm14 = xptr(srcreg+d2+d1+32);	/*#2 I4				; 37	avail 5 storable 11,4,10,15 */ \
\
	xmm5 = xmm7;			/*#2 Copy R3			; 38-40	avail none storable 11,4,10,15 */ \
	xmm7 *= xmm8;			/*#2 A3 = R3 * cosine/sine	; 38-42 */ \
	xptr(srcreg+d2+d1+48) = xmm11;	/* Save I4			; 38	avail 11 storable 4,10,15 */ \
\
	xmm11 = xmm14;			/*#2 Copy I4			; 39-41	avail none storable 4,10,15 */ \
	xmm14 *= xmm8;			/*#2 B4 = I4 * cosine/sine	; 39-43 */ \
	xptr(srcreg+d2+d1+32) = xmm4;		/* Save I2			; 39 (read in at clock 36) avail 4 storable 10,15 */ \
\
	xmm4 = xptr(srcreg+d2+32);		/*#2 I3				; 40	avail none storable 10,15 */ \
	xmm8 *= xmm4;			/*#2 B3 = I3 * cosine/sine	; 40-44 */ \
	xptr(srcreg+d2+d1+16) = xmm10;	/* Save R4			; 40	avail 10 storable 15 */ \
\
	xmm10 = xptr(pmreg+pmoff+64);		/* sine 3			; 41	avail none storable 15 */ \
	xmm0 *= xmm10;			/* A3 = A3 * sine (final R3)	; 41-45	avail none storable 15,0 */ \
	xptr(srcreg+d2+d1) = xmm15;		/* Save R2			; 41 (read in at clock 34) avail 15 storable 0 */ \
\
	xmm1 -= xmm11;			/*#2 A4 = A4 - I4		; 42-44	avail 15,11 storable 0 */ \
	xmm2 *= xmm10;			/* B3 = B3 * sine (final I3)	; 42-46	avail 15,11,10 storable 0.2 */ \
	xmm11 = xptr(pmreg+pmoff+0);		/* sine 1			; 42	avail 15,10 storable 0,2 */ \
\
	xmm7 += xmm4;			/*#2 A3 = A3 + I3		; 43-45	avail 15,10,4 storable 0,2 */ \
	xmm9 *= xmm11;			/* A1 = A1 * sine (final R1)	; 43-47	avail 15,10,4 storable 0,2,9 */ \
	xmm4 = xptr(screg+0);			/*#2 sine 3&4			; 43	avail 15,10 storable 0,2,9 */ \
\
	xmm14 += xmm3;			/*#2 B4 = B4 + R4		; 44-46	avail 15,10,3 storable 0,2,9 */ \
	xmm6 *= xmm11;			/* B1 = B1 * sine (final I1)	; 44-48	avail 15,10,3,11 storable 0,2,9,6 */ \
	xmm11 = xmm12;			/*#2 Copy R2			; 44-46 (20) avail 15,10,3 storable 0,2,9,6 */ \
	xmm3 = xptr(srcreg);			/*#2 R1				; 44	avail 15,10 storable 0,2,9,6 */ \
\
	xmm8 -= xmm5;			/*#2 B3 = B3 - R3		; 45-47	avail 15,10,5 storable 0,2,9,6 */ \
	xmm1 *= xmm4;			/*#2 A4 = A4 * sine (R4)	; 45-49 */ \
	xmm10 = xptr(srcreg+32);		/*#2 I1				; 45	avail 15,5 storable 0,2,9,6 */ \
\
	xmm12 += xmm3;			/*#2 R2 = R1 + R2 (new R1)	; 46-48 */ \
	xmm7 *= xmm4;			/*#2 A3 = A3 * sine (R3)	; 46-50 */ \
	xmm5 = xmm13;			/*#2 Copy I2			; 46-48 (21) avail 15 storable 0,2,9,6 */ \
	xptr(srcreg+d1+16) = xmm0;		/* Save R3			; 46	avail 15,0 storable 2,9,6 */ \
\
	xmm3 -= xmm11;			/*#2 R1 = R1 - R2 (new R2)	; 47-49	avail 15,0,11 storable 2,9,6 */ \
	xmm14 *= xmm4;			/*#2 B4 = B4 * sine (I4)	; 47-51 */ \
	xptr(srcreg+d1+48) = xmm2;		/* Save I3			; 47	avail 15,0,11,2 storable 9,6 */ \
\
	xmm13 += xmm10;			/*#2 I2 = I1 + I2 (new I1)	; 48-50 */ \
	xmm8 *= xmm4;			/*#2 B3 = B3 * sine (I3)	; 48-52	avail 15,0,11,2,4 storable 9,6 */ \
	xptr(srcreg+d1) = xmm9;		/* Save R1			; 48 (read in at clock 7) avail 15,0,11,2,4,9 storable 6 */ \
\
	xmm10 -= xmm5;			/*#2 I1 = I1 - I2 (new I2)	; 49-51	avail 15,0,11,2,4,9,5 storable 6 */ \
	xptr(srcreg+d1+32) = xmm6;		/* Save I1			; 49 (read at clock 8) avail 15,0,11,2,4,9,5,6 */ \
\
	xprefetchw(u8ptr(srcreg+srcinc)); \
\
	xmm0 = xmm1;			/*#2 Copy R4			; 50-52 (50 on Core i7) avail 15,11,2,4,9,5,6 */ \
;										/* 50 STALL */ \
\
	xmm1 -= xmm7;			/*#2 R4 = R4 - R3 (new I4)	; 51-52 */ \
\
	xmm7 += xmm0;			/*#2 R3 = R4 + R3 (new R3)	; 52-53	avail 15,11,2,4,9,5,6,0 */ \
	xmm2 = xmm14;			/*#2 Copy I4			; 52-54 (52 on Core i7) avail 15,11,4,9,5,6,0 */ \
\
	xmm14 += xmm8;			/*#2 I4 = I3 + I4 (new I3)	; 53-55 */ \
	xmm4 = xmm10;			/*#2 Copy I2			; 53-55 (52) avail 15,11,9,5,6,0 */ \
\
	xmm8 -= xmm2;			/*#2 I3 = I3 - I4 (new R4)	; 54-56	avail 15,11,9,5,6,0,2 */ \
	xmm6 = xptr(pmreg+96+16);		/*#2 cosine/sine 4		; 54	avail 15,11,9,5,0,2 */ \
\
	xprefetchw(u8ptr(srcreg+srcinc+d1)); \
\
	xmm10 -= xmm1;			/*#2 I2 = I2 - I4 (newer I4)	; 55-57 */ \
	xmm5 = xmm3;			/*#2 Copy R2			; 55-57 (50) avail 15,11,9,0,2 */ \
	xmm0 = xptr(pmreg+32+16);		/*#2 cosine/sine 2		; 55	avail 15,11,9,2 */ \
\
	xmm1 += xmm4;			/*#2 I4 = I2 + I4 (newer I2)	; 56-58	avail 15,11,9,2,4 */ \
	xmm9 = xmm12;			/*#2 Copy R1			; 56-58 (49) avail 15,11,2,4 */ \
\
	xmm3 -= xmm8;			/*#2 R2 = R2 - R4 (newer R4)	; 57-59 */ \
	xmm11 = xmm13;			/*#2 Copy I1			; 59-61 (52) avail 15,2,4 */ \
\
	xmm8 += xmm5;			/*#2 R4 = R2 + R4 (newer R2)	; 58-60	avail 15,2,4,5 */ \
	xmm2 = xmm10;			/*#2 Copy I4			; 58-60	avail 15,4,5 */ \
	xmm10 *= xmm6;			/*#2 B4 = I4 * cosine/sine	; 58-62 */ \
\
	xmm12 -= xmm7;			/*#2 R1 = R1 - R3 (newer R3)	; 59-61 */ \
	xmm4 = xmm1;			/*#2 Copy I2			; 59-61	avail 15,5 */ \
	xmm1 *= xmm0;			/*#2 B2 = I2 * cosine/sine	; 59-63 */ \
\
	xmm13 -= xmm14;			/*#2 I1 = I1 - I3 (newer I3)	; 60-62 */ \
	xmm6 *= xmm3;			/*#2 A4 = R4 * cosine/sine	; 60-64 */ \
	xmm5 = xptr(pmreg+64+16);		/*#2 cosine/sine 3		; 60	avail 15 */ \
\
	xmm7 += xmm9;			/*#2 R3 = R1 + R3 (newer R1)	; 61-63	avail 15,9 */ \
	xmm0 *= xmm8;			/*#2 A2 = R2 * cosine/sine	; 61-65 */ \
	xmm15 = xptr(pmreg+0+16);		/*#2 cosine/sine 1		; 61	avail 9 */ \
\
	xmm14 += xmm11;			/*#2 I3 = I1 + I3 (newer I1)	; 62-64	avail 9,11 */ \
	xmm9 = xmm12;			/*#2 Copy R3			; 62-64	avail 11 */ \
	xmm12 *= xmm5;			/*#2 A3 = R3 * cosine/sine	; 62-66 */ \
\
	xmm10 -= xmm3;			/*#2 B4 = B4 - R4		; 63-65	avail 11,3 */ \
	xmm5 *= xmm13;			/*#2 B3 = I3 * cosine/sine	; 63-67 */ \
	xmm11 = xptr(pmreg+96);		/*#2 sine 4			; 63	avail 3 */ \
\
	xmm1 -= xmm8;			/*#2 B2 = B2 - R2		; 64-66	avail 3,8 */ \
	xmm3 = xmm7;			/*#2 Copy R1			; 64-66	avail 8 */ \
	xmm7 *= xmm15;			/*#2 A1 = R1 * cosine/sine	; 64-68 */ \
\
	xmm6 += xmm2;			/*#2 A4 = A4 + I4		; 65-67	avail 8,2 */ \
	xmm15 *= xmm14;			/*#2 B1 = I1 * cosine/sine	; 65-69 */ \
	xmm2 = xptr(pmreg+32);		/*#2 sine 2			; 65	avail 8 */ \
\
	xmm0 += xmm4;			/*#2 A2 = A2 + I2		; 66-68	avail 8,4 */ \
	xmm10 *= xmm11;			/*#2 B4 = B4 * sine (final I4)	; 66-70	avail 8,4 storable 10 */ \
\
	xmm12 += xmm13;			/*#2 A3 = A3 + I3		; 67-69	avail 8,4,13 storable 10 */ \
	xmm1 *= xmm2;			/*#2 B2 = B2 * sine (final I2)	; 67-71	avail 8,4,13 storable 10,1 */ \
	xmm4 = xptr(pmreg+64);		/*#2 sine 3			; 67	avail 8,13 storable 10,1 */ \
\
	xmm5 -= xmm9;			/*#2 B3 = B3 - R3		; 68-70	avail 8,13,9 storable 10,1 */ \
	xmm6 *= xmm11;			/*#2 A4 = A4 * sine (final R4)	; 68-72	avail 8,13,9,11 storable 10,1,6 */ \
	xmm8 = xptr(pmreg+0);			/*#2 sine 1			; 68 */ \
\
	xmm7 += xmm14;			/*#2 A1 = A1 + I1		; 69-71	avail 13,9,11,14 storable 10,1,6 */ \
	xmm0 *= xmm2;			/*#2 A2 = A2 * sine (final R2)	; 69-73	avail 13,9,11,14,2 storable 10,1,6,0 */ \
\
	xmm15 -= xmm3;			/*#2 B1 = B1 - R1		; 70-72	avail 13,9,11,14,2,3 storable 10,1,6,0 */ \
	xmm12 *= xmm4;			/*#2 A3 = A3 * sine (final R3)	; 70-74	avail 13,9,11,14,2,3 storable 10,1,6,0,12 */ \
\
	xmm5 *= xmm4;			/*#2 B3 = B3 * sine (final I3)	; 71-75	avail 13,9,11,14,2,3,4 storable 10,1,6,0,12,5 */ \
	xptr(srcreg+d2+48) = xmm10;		/*#2 Save I4			; 71	avail 13,9,11,14,2,3,4,10 storable 1,6,0,12,5 */ \
\
	xmm7 *= xmm8;			/*#2 A1 = A1 * sine (final R1)	; 72-76	avail 13,9,11,14,2,3,4,10 storable 1,6,0,12,5,7 */ \
	xptr(srcreg+d2+32) = xmm1;		/*#2 Save I2			; 72	avail 13,9,11,14,2,3,4,10,1 storable 6,0,12,5,7 */ \
\
	xmm15 *= xmm8;			/*#2 B1 = B1 * sine (final I1)	; 73-77	avail 13,9,11,14,2,3,4,10,1,8 storable 6,0,12,5,7,15 */ \
	xptr(srcreg+d2+16) = xmm6;		/*#2 Save R4			; 73 */ \
\
	xptr(srcreg+d2) = xmm0;		/*#2 Save R2			; 74 */ \
	xptr(srcreg+16) = xmm12;		/*#2 Save R3			; 75 */ \
	xptr(srcreg+48) = xmm5;		/*#2 Save I3			; 76 */ \
	xptr(srcreg) = xmm7;			/*#2 Save R1			; 77 */ \
	xptr(srcreg+32) = xmm15;		/*#2 Save I1			; 78 */ \
\
	srcreg += srcinc; \
\

#endif
\
/* */ \
/* ******************************* four-complex-with-partial-normalization variants ************************************* */ \
/* */ \
/* These macros are used in pass 1 of r4dwpn two pass FFTs.  They are like the standard four-complex */ \
/* DJBFFT macros except that a normalization multiplier has been pre-applied to the sine multiplier. */ \
/* Consequently, the forward FFT and inverse FFT use different sine multipliers. */ \
/* Also, a normalization multiplier must be applied to the final R1/I1 value. */ \
/* */ \
\

#define r4_x4cl_wpn_four_complex_djbfft_preload \
	r4_x4c_wpn_djbfft_mem_preload; \
;	/* r4_x4c_wpn_djbfft_partial_mem_preload -- assume same as r4_x4c_wpn_djbfft_mem_preload */ \

#ifdef LOW_CPU_REGISTER_COUNT
#define r4_x4cl_wpn_four_complex_djbfft(srcreg,srcinc,d1,d2,screg) { \
	vec2f64 *sd0 = (vec2f64*)(srcreg); \
	vec2f64 *sd1 = (vec2f64*)(srcreg+d1); \
	vec2f64 *sd2 = (vec2f64*)(srcreg+d2); \
	vec2f64 *sd3 = (vec2f64*)(srcreg+d2+d1); \
	vec2f64 *sc = (vec2f64*)(screg); \
	r4_x4c_wpn_djbfft_mem4(sd0[0],sd1[0],sd2[0],sd3[0],sd0[1],sd1[1],sd2[1],sd3[1], sc[2],sc[3],sc[5],sc[6], sc[0], srcreg+srcinc,d1,sd0[0],sd0[1]); \
/*	xstore	[srcreg], xmm0		;; Save R1 */ \
/*	xstore	[srcreg+16], xmm0	;; Save I1 */ \
	xmm0 = sd0[2];	/* R1 */ \
	xmm5 = sd0[3];	/* R5 */ \
	sd0[2] = xmm7;	/* Save R2 */ \
	sd0[3] = xmm6;	/* Save I2 */ \
	xmm7 = sd1[2];	/* R2 */ \
	xmm6 = sd1[3];	/* R6 */ \
	sd1[0] = xmm3;	/* Save R3 */ \
	sd1[1] = xmm1;	/* Save I3 */ \
	sd1[2] = xmm2;	/* Save R4 */ \
	sd1[3] = xmm4;	/* Save I4 */ \
	r4_x4c_wpn_djbfft_partial_mem4(xmm0,xmm7,xmm3,xmm2,xmm5,xmm6,xmm1,xmm4, sd2[2],sd3[2],sd2[3],sd3[3], sc[2],sc[3],sc[5],sc[6], sc[0], srcreg+srcinc+d2,d1,sd2[0],sd2[1]); \
/*	xstore	[srcreg+d2], xmm0	;; Save R1 */ \
/*	xstore	[srcreg+d2+16], xmm0	;; Save I1 */ \
	sd2[2] = xmm4;	/* Save R2 */ \
	sd2[3] = xmm1;	/* Save I2 */ \
	sd3[0] = xmm2;	/* Save R3 */ \
	sd3[1] = xmm7;	/* Save I3 */ \
	sd3[2] = xmm0;	/* Save R4 */ \
	sd3[3] = xmm5;	/* Save I4 */ \
	srcreg += srcinc; \
}
#else
#define r4_x4cl_wpn_four_complex_djbfft(srcreg,srcinc,d1,d2,screg) { \
	vec2f64 *sd0 = (vec2f64*)(srcreg); \
	vec2f64 *sd1 = (vec2f64*)(srcreg+d1); \
	vec2f64 *sd2 = (vec2f64*)(srcreg+d2); \
	vec2f64 *sd3 = (vec2f64*)(srcreg+d2+d1); \
	vec2f64 *sc = (vec2f64*)(screg); \
	vec2f64 normreg = sc[0]; \
	vec2f64 sc0 = sc[2]; \
	vec2f64 sc1 = sc[3]; \
	vec2f64 sc2 = sc[5]; \
	vec2f64 sc3 = sc[6]; \
	r4_x4c_wpn_djbfft_mem4(sd0[0],sd1[0],sd2[0],sd3[0],sd0[1],sd1[1],sd2[1],sd3[1], sc0,sc1,sc2,sc3, normreg,srcreg+srcinc,d1,sd0[0],sd0[1]); \
/*	xstore	[srcreg], xmm0		;; Save R1 */ \
/*	xstore	[srcreg+16], xmm0	;; Save I1 */ \
	xmm0 = sd0[2];	/* R1 */ \
	xmm5 = sd0[3];	/* R5 */ \
	sd0[2] = xmm7;	/* Save R2 */ \
	sd0[3] = xmm6;	/* Save I2 */ \
	xmm7 = sd1[2];	/* R2 */ \
	xmm6 = sd1[3];	/* R6 */ \
	sd1[0] = xmm3;	/* Save R3 */ \
	sd1[1] = xmm1;	/* Save I3 */ \
	sd1[2] = xmm2;	/* Save R4 */ \
	sd1[3] = xmm4;	/* Save I4 */ \
	r4_x4c_wpn_djbfft_partial_mem4(xmm0,xmm7,xmm3,xmm2,xmm5,xmm6,xmm1,xmm4, sd2[2],sd3[2],sd2[3],sd3[3], sc0,sc1,sc2,sc3, normreg,srcreg+srcinc+d2,d1,sd2[0],sd2[1]); \
/*	xstore	[srcreg+d2], xmm0	;; Save R1 */ \
/*	xstore	[srcreg+d2+16], xmm0	;; Save I1 */ \
	sd3[0] = xmm2;	/* Save R3 */ \
	sd3[1] = xmm7;	/* Save I3 */ \
	sd3[2] = xmm0;	/* Save R4 */ \
	sd3[3] = xmm5;	/* Save I4 */ \
	sd2[2] = xmm4;	/* Save R2 */ \
	sd2[3] = xmm1;	/* Save I2 */ \
	srcreg += srcinc; \
}
#endif

#define r4_x4c_wpn_djbfft_mem_preload

//#define r4_x4c_wpn_djbfft_mem(R1,R2,R3,R4,R5,R6,R7,R8,screg1,screg2,normreg,pre1,pre2,dst1,dst2)

#define r4_x4c_wpn_djbfft_mem4(R1,R2,R3,R4,R5,R6,R7,R8,SC0,SC1,SC2,SC3,normreg,pre1,pre2,dst1,dst2) { \
	vec2f64 s0,s1; \
	xmm0 = R1;			/* R1 */ \
	xmm2 = R3;			/* R3 */ \
	xmm2 += xmm0;		/* R3 = R1 + R3 (new R1) */ \
\
	xmm1 = R2;			/* R2 */ \
	xmm3 = R4;			/* R4 */ \
	xmm3 += xmm1;		/* R4 = R2 + R4 (new R2) */ \
\
	xmm4 = R5;			/* I1 */ \
	xmm6 = R7;			/* I3 */ \
	xmm6 += xmm4;		/* I3 = I1 + I3 (new I1) */ \
\
	xmm0 -= R3;			/* R1 = R1 - R3 (new R3) */ \
	xmm1 -= R4;			/* R2 = R2 - R4 (new R4) */ \
\
	xprefetchw(u8ptr(pre1)); \
\
	xmm5 = xmm3;		/* Copy R2 */ \
	xmm3 += xmm2;		/* R2 = R1 + R2 (final R1) */ \
	xmm2 -= xmm5;		/* R1 = R1 - R2 (final R2) */ \
\
	xmm5 = R6;			/* I2 */ \
	xmm7 = R8;			/* I4 */ \
	xmm7 += xmm5;		/* I4 = I2 + I4 (new I2) */ \
\
	xmm4 -= R7;			/* I1 = I1 - I3 (new I3) */ \
	xmm5 -= R8;			/* I2 = I2 - I4 (new I4) */ \
\
	s0 = normreg; \
	xmm3 *= s0;			/* Apply normalization multiplier to R1 */ \
	dst1 = xmm3;		/* Save R1 */ \
\
	xmm3 = xmm7;		/* Copy I2 */ \
	xmm7 += xmm6;		/* I2 = I1 + I2 (final I1) */ \
\
	xmm6 -= xmm3;		/* I1 = I1 - I2 (final I2) */ \
\
	xmm7 *= s0;			/* Apply normalization multiplier to I1 */ \
	dst2 = xmm7;		/* Save I1 */ \
\
	xmm3 = xmm1; \
	xmm1 += xmm4;		/* R4 = I3 + R4 (final I3) */ \
\
	xmm7 = xmm0;		/* Copy R3 */ \
	xmm0 -= xmm5;		/* R3 = R3 - I4 (final R3) */ \
\
	xmm4 -= xmm3;		/* I3 = I3 - R4 (final I4) */ \
	xmm5 += xmm7;		/* I4 = R3 + I4 (final R4) */ \
\
	xprefetchw(u8ptr(pre1+pre2)); \
\
	s1 = SC3; \
	xmm7 = s1 * xmm2 - xmm6;	/* A2 = R2 * cosine/sine - I2 */ \
	xmm6 = xmm6 * s1 + xmm2;	/* B2 = I2 * cosine/sine + R2 */ \
	s0 = SC1; \
	xmm3 = s0 * xmm0 - xmm1;	/* A3 = R3 * cosine/sine - I3 */ \
	xmm1 = xmm1 * s0 + xmm0;	/* B3 = I3 * cosine/sine + R3 */ \
	xmm2 = s0 * xmm5 + xmm4;	/* A4 = R4 * cosine/sine + I4 */ \
	xmm4 = xmm4 * s0 - xmm5;	/* B4 = I4 * cosine/sine - R4 */ \
\
	s1 = SC2; \
	xmm5 = SC0;			/* Sine */ \
	xmm7 *= s1;			/* A2 = A2 * sine (final R2) */ \
	xmm6 *= s1;			/* B2 = B2 * sine (final I2) */ \
	xmm3 *= xmm5;		/* A3 = A3 * sine (final R3) */ \
	xmm1 *= xmm5;		/* B3 = B3 * sine (final I3) */ \
	xmm2 *= xmm5;		/* A4 = A4 * sine (final R4) */ \
	xmm4 *= xmm5;		/* B4 = B4 * sine (final I4) */ \
}


#define r4_x4c_wpn_djbfft_partial_mem_preload

/*#define r4_x4c_wpn_djbfft_partial_mem(r1,r2,r3,r4,r5,r6,r7,r8, mem3,mem4,mem7,mem8, screg1,screg2,normreg, pre1,pre2,dst1,dst2) \
	r4_x4c_wpn_djbfft_partial_mem4(r1,r2,r3,r4,r5,r6,r7,r8, mem3,mem4,mem7,mem8, xptr(screg1),xptr(screg1+16),xptr(screg2),xptr(screg2+16), xptr(normreg), pre1,pre2,dst1,dst2)
*/
#define r4_x4c_wpn_djbfft_partial_mem4(r1,r2,r3,r4,r5,r6,r7,r8, mem3,mem4,mem7,mem8, SC0,SC1,SC2,SC3, normreg, pre1,pre2,dst1,dst2) { \
	vec2f64 s1,s2; \
	s1 = r3 = mem3; \
	r3 += r1;			/* R3 = R1 + R3 (new R1) */ \
	r1 -= s1;			/* R1 = R1 - R3 (new R3) */ \
\
	s1 = r4 = mem4; \
	r4 += r2;			/* R4 = R2 + R4 (new R2) */ \
	r2 -= s1;			/* R2 = R2 - R4 (new R4) */ \
\
	r7 = r3; \
	r3 -= r4;			/* R1 = R1 - R2 (final R2) */ \
	r4 += r7;			/* R2 = R1 + R2 (final R1) */ \
\
	r4 *= normreg;		/* Apply normalization multiplier to R1 */ \
	dst1 = r4;			/* Save R1 */ \
\
	r7 = mem7; \
	r8 = mem8; \
\
	r4 = r6; \
	r6 -= r8;			/* I2 = I2 - I4 (new I4) */ \
	r8 += r4;			/* I4 = I2 + I4 (new I2) */ \
\
	r4 = r5; \
	r5 -= r7;			/* I1 = I1 - I3 (new I3) */ \
	r7 += r4;			/* I3 = I1 + I3 (new I1) */ \
\
	r4 = r1; \
	r1 -= r6;			/* R3 = R3 - I4 (final R3) */ \
	r6 += r4;			/* I4 = R3 + I4 (final R4) */ \
\
	r4 = r2; \
	r2 += r5;			/* R4 = I3 + R4 (final I3) */ \
	r5 -= r4;			/* I3 = I3 - R4 (final I4) */ \
\
	r4 = r7; \
	r7 -= r8;			/* I1 = I1 - I2 (final I2) */ \
	r8 += r4;			/* I2 = I1 + I2 (final I1) */ \
\
	xprefetchw(u8ptr(pre1)); \
\
\
	r8 *= normreg;		/* Apply normalization multiplier to I1 */ \
	dst2 = r8;			/* Save I1 */ \
\
	s1 = SC1; \
	r4 = s1 * r1 - r2;	/* A3 = R3 * cosine/sine - I3 */ \
	r2 = r2 * s1 + r1;	/* B3 = I3 * cosine/sine + R3 */ \
	r1 = s1 * r6 + r5;	/* A4 = R4 * cosine/sine + I4 */ \
	r5 = r5 * s1 - r6;	/* B4 = I4 * cosine/sine - R4 */ \
	s2 = SC3; \
	r8 = s2 * r3 - r7;	/* A2 = R2 * cosine/sine - I2 */ \
	r7 = r7 * s2 + r3;	/* B2 = I2 * cosine/sine + R2 */ \
\
	xprefetchw(u8ptr(pre1+pre2)); \
\
	r3 = SC0;			/* sine */ \
	s2 = SC2; \
	r4 *= r3;			/* A3 = A3 * sine (final R3) */ \
	r2 *= r3;			/* B3 = B3 * sine (final I3) */ \
	r1 *= r3;			/* A4 = A4 * sine (final R4) */ \
	r5 *= r3;			/* B4 = B4 * sine (final I4) */ \
	r8 *= s2;			/* A2 = A2 * sine (final R2) */ \
	r7 *= s2;			/* B2 = B2 * sine (final I2) */ \
}

/* 32-bit AMD K8 optimized versions of the above macros */ \

#if 0 // (@INSTR(,%xarch,<K8>) NE 0);
\

#define r4_x4c_wpn_djbfft_mem(R1,R2,R3,R4,R5,R6,R7,R8,screg1,screg2,normreg,pre1,pre2,dst1,dst2) \
	xmm0 = R1;		/* R1 */ \
	xmm2 = R3;		/* R3 */ \
	xmm2 += xmm0;		/* R3 = R1 + R3 (new R1) */ \
\
	xmm1 = R2;		/* R2 */ \
	xmm3 = R4;		/* R4 */ \
	xmm3 += xmm1;		/* R4 = R2 + R4 (new R2) */ \
\
	xmm4 = R5;		/* I1 */ \
	xmm6 = R7;		/* I3 */ \
	xmm6 += xmm4;		/* I3 = I1 + I3 (new I1) */ \
\
	xmm0 -= R3;		/* R1 = R1 - R3 (new R3) */ \
	xmm1 -= R4;		/* R2 = R2 - R4 (new R4) */ \
\
	xprefetchw(u8ptr(pre1)); \
\
	xmm3 += xmm2;		/* R2 = R1 + R2 (final R1) */ \
	xmm2 += xmm2; /* multwo */; \
	xmm2 -= xmm3;		/* R1 = R1 - R2 (final R2) */ \
\
	xmm5 = R6;		/* I2 */ \
	xmm7 = R8;		/* I4 */ \
	xmm7 += xmm5;		/* I4 = I2 + I4 (new I2) */ \
\
	xmm4 -= R7;		/* I1 = I1 - I3 (new I3) */ \
	xmm5 -= R8;		/* I2 = I2 - I4 (new I4) */ \
\
	xmm3 *= xptr(normreg);		/* Apply normalization multiplier to R1 */ \
	dst1 = xmm3;		/* Save R1 */ \
\
	xmm7 += xmm6;		/* I2 = I1 + I2 (final I1) */ \
	xmm6 += xmm6; /* multwo */; \
	xmm6 -= xmm7;		/* I1 = I1 - I2 (final I2) */ \
\
	xmm7 *= xptr(normreg);		/* Apply normalization multiplier to I1 */ \
	dst2 = xmm7;		/* Save I1 */ \
\
	xmm1 += xmm4;		/* R4 = I3 + R4 (final I3) */ \
	xmm4 += xmm4; /* multwo */; \
\
	xmm0 -= xmm5;		/* R3 = R3 - I4 (final R3) */ \
	xmm5 += xmm5; /* multwo */; \
	xmm4 -= xmm1;		/* I3 = I3 - R4 (final I4) */ \
	xmm5 += xmm0;		/* I4 = R3 + I4 (final R4) */ \
\
	xprefetchw(u8ptr(pre1+pre2)); \
\
	xmm7 = xptr(screg2+16);	/* cosine/sine */ \
	xmm7 *= xmm2;		/* A2 = R2 * cosine/sine */ \
	xmm7 -= xmm6;		/* A2 = A2 - I2 */ \
\
	xmm3 = xptr(screg1+16);	/* cosine/sine */ \
	xmm3 *= xmm0;		/* A3 = R3 * cosine/sine */ \
	xmm3 -= xmm1;		/* A3 = A3 - I3 */ \
\
	xmm6 *= xptr(screg2+16);	/* B2 = I2 * cosine/sine */ \
	xmm6 += xmm2;		/* B2 = B2 + R2 */ \
\
	xmm2 = xptr(screg1+16);	/* cosine/sine */ \
	xmm2 *= xmm5;		/* A4 = R4 * cosine/sine */ \
	xmm2 += xmm4;		/* A4 = A4 + I4 */ \
\
	xmm1 *= xptr(screg1+16);	/* B3 = I3 * cosine/sine */ \
	xmm1 += xmm0;		/* B3 = B3 + R3 */ \
\
	xmm4 *= xptr(screg1+16);	/* B4 = I4 * cosine/sine */ \
	xmm4 -= xmm5;		/* B4 = B4 - R4 */ \
\
	xmm7 *= xptr(screg2);		/* A2 = A2 * sine (final R2) */ \
	xmm6 *= xptr(screg2);		/* B2 = B2 * sine (final I2) */ \
\
	xmm5 = xptr(screg1);		/* Sine */ \
	xmm3 *= xmm5;		/* A3 = A3 * sine (final R3) */ \
	xmm1 *= xmm5;		/* B3 = B3 * sine (final I3) */ \
\
	xmm2 *= xmm5;		/* A4 = A4 * sine (final R4) */ \
	xmm4 *= xmm5;		/* B4 = B4 * sine (final I4) */ \
\
\

#define r4_x4c_wpn_djbfft_partial_mem(r1,r2,r3,r4,r5,r6,r7,r8,mem3,mem4,mem7,mem8,screg1,screg2,normreg,pre1,pre2,dst1,dst2) \
	r3 = mem3; \
	r3 += r1;			/* R3 = R1 + R3 (new R1) */ \
	r1 -= mem3;		/* R1 = R1 - R3 (new R3) */ \
\
	r4 = mem4; \
	r4 += r2;			/* R4 = R2 + R4 (new R2) */ \
	r2 -= mem4;		/* R2 = R2 - R4 (new R4) */ \
\
	r3 -= r4;			/* R1 = R1 - R2 (final R2) */ \
	r4 += r4; /* multwo */ \
	r4 += r3;			/* R2 = R1 + R2 (final R1) */ \
\
	r4 *= xptr(normreg);		/* Apply normalization multiplier to R1 */ \
	dst1 = r4;		/* Save R1 */ \
\
	r7 = mem7; \
	r8 = mem8; \
\
	r6 -= r8;			/* I2 = I2 - I4 (new I4) */ \
	r8 += r8; /* multwo */ \
	r8 += r6;			/* I4 = I2 + I4 (new I2) */ \
\
	r5 -= r7;			/* I1 = I1 - I3 (new I3) */ \
	r7 += r7; /* multwo */ \
	r7 += r5;			/* I3 = I1 + I3 (new I1) */ \
\
	r1 -= r6;			/* R3 = R3 - I4 (final R3) */ \
	r6 += r6; /* multwo */ \
	r6 += r1;			/* I4 = R3 + I4 (final R4) */ \
\
	r2 += r5;			/* R4 = I3 + R4 (final I3) */ \
	r5 += r5; /* multwo */ \
	r5 -= r2;			/* I3 = I3 - R4 (final I4) */ \
\
	xprefetchw(u8ptr(pre1)); \
\
	r7 -= r8;			/* I1 = I1 - I2 (final I2) */ \
	r8 += r8; /* multwo */ \
	r8 += r7;			/* I2 = I1 + I2 (final I1) */ \
\
	r4 = xptr(screg1+16);		/* cosine/sine */ \
	r4 *= r1;			/* A3 = R3 * cosine/sine */ \
\
	r8 *= xptr(normreg);		/* Apply normalization multiplier to I1 */ \
	dst2 = r8;		/* Save I1 */ \
\
	r8 = xptr(screg2+16);		/* cosine/sine */ \
	r8 *= r3;			/* A2 = R2 * cosine/sine */ \
\
	r4 -= r2;			/* A3 = A3 - I3 */ \
	r2 *= xptr(screg1+16);		/* B3 = I3 * cosine/sine */ \
	r2 += r1;			/* B3 = B3 + R3 */ \
\
	r1 = xptr(screg1+16);		/* cosine/sine */ \
	r1 *= r6;			/* A4 = R4 * cosine/sine */ \
\
	r8 -= r7;			/* A2 = A2 - I2 */ \
	r7 *= xptr(screg2+16);		/* B2 = I2 * cosine/sine */ \
\
	xprefetchw(u8ptr(pre1+pre2)); \
\
	r1 += r5;			/* A4 = A4 + I4 */ \
	r5 *= xptr(screg1+16);		/* B4 = I4 * cosine/sine */ \
\
	r7 += r3;			/* B2 = B2 + R2 */ \
	r5 -= r6;			/* B4 = B4 - R4 */ \
\
	r3 = xptr(screg1);		/* sine */ \
	r4 *= r3;			/* A3 = A3 * sine (final R3) */ \
	r2 *= r3;			/* B3 = B3 * sine (final I3) */ \
	r8 *= xptr(screg2);		/* A2 = A2 * sine (final R2) */ \
	r1 *= r3;			/* A4 = A4 * sine (final R4) */ \
	r7 *= xptr(screg2);		/* B2 = B2 * sine (final I2) */ \
	r5 *= r3;			/* B4 = B4 * sine (final I4) */ \
\

#endif
\
/* 64-bit Intel implementations of the above - use the extra XMM registers */ \

#if 0 // (@INSTR(,%xarch,<CORE>) NE 0) OR (@INSTR(,%xarch,<P4>) NE 0) OR (@INSTR(,%xarch,<BLEND>) NE 0) OR (@INSTR(,%xarch,<K10>) NE 0);
#if 0 //DEF X86_64;
\
/* Theoretical best case is 44 clocks on a Core 2.  Now at ?? clocks. */ \
\

#define r4_x4cl_wpn_four_complex_djbfft(srcreg,srcinc,d1,d2,screg) \
	xmm0 = xptr(srcreg);		/* R1 */ \
	xmm2 = xptr(srcreg+d2);	/* R3 */ \
	xmm5 = xptr(srcreg+d1+16);	/* I2 */ \
	xmm7 = xptr(srcreg+d2+d1+16);	/* I4 */ \
	xmm1 = xptr(srcreg+d1);	/* R2 */ \
	xmm3 = xptr(srcreg+d2+d1);	/* R4 */ \
	xmm4 = xptr(srcreg+16);	/* I1 */ \
	xmm6 = xptr(srcreg+d2+16);	/* I3 */ \
\
	xmm8 = xmm0;		/* Copy R1				; 1-3 */ \
	xmm0 -= xmm2;		/* R1 = R1 - R3 (new R3)		; 1-3 */ \
\
	xmm9 = xmm5;		/* Copy I2				; 2-4 */ \
	xmm5 -= xmm7;		/* I2 = I2 - I4 (new I4)		; 2-4 */ \
\
	xmm10 = xmm1;		/* Copy R2				; 3-5 */ \
	xmm1 -= xmm3;		/* R2 = R2 - R4 (new R4)		; 3-5 */ \
\
	xmm11 = xmm4;		/* Copy I1				; 4-6 */ \
	xmm4 -= xmm6;		/* I1 = I1 - I3 (new I3)		; 4-6	avail 12+ */ \
	xmm12 = xptr(screg+32+16);	/* cosine/sine				; 4 */ \
\
	xmm2 += xmm8;		/* R3 = R1 + R3 (new R1)		; 5-7	avail 8,13+ */ \
	xmm13 = xmm12;							/* 5-7 */ \
\
	xmm3 += xmm10;		/* R4 = R2 + R4 (new R2)		; 6-8	avail 8,10,14+ */ \
	xmm14 = xmm12;							/* 6-8 */ \
\
	xmm7 += xmm9;		/* I4 = I2 + I4 (new I2)		; 7-9	avail 8.10.9,15 */ \
	xmm8 = xmm5;		/* Copy I4				; 7-9 */ \
\
	xmm6 += xmm11;		/* I3 = I1 + I3 (new I1)		; 8-10	avail 10.9,15,11 */ \
	xmm15 = xmm12;							/* 8-10 */ \
\
	xmm5 += xmm0;		/* I4 = R3 + I4 (final R4)		; 9-11	avail 10.9,11 */ \
	xmm9 = xmm4;		/* Copy I3				; 9-11 */ \
	xmm10 = xptr(screg+80+16);	/* cosine/sine				; 9 */ \
\
	xmm0 -= xmm8;		/* R3 = R3 - I4 (final R3)		; 10-12	avail 11,8 */ \
\
	xprefetchw(u8ptr(srcreg+srcinc)); \
\
	xmm4 -= xmm1;		/* I3 = I3 - R4 (final I4)		; 11-13	avail 11,8 */ \
	xmm11 = xmm10;							/* 11-13 */ \
\
	xmm1 += xmm9;		/* R4 = I3 + R4 (final I3)		; 12-14	avail 8,9 */ \
	xmm12 *= xmm5;		/* A4 = R4 * cosine/sine		; 12-16 */ \
	xmm8 = xmm2;		/* Copy R1				; 12-14 */ \
\
	xmm2 -= xmm3;		/* R1 = R1 - R2 (final R2)		; 13-15	avail 9 */ \
	xmm13 *= xmm0;		/* A3 = R3 * cosine/sine		; 13-17 */ \
	xmm9 = xmm6;		/* Copy I1				; 13-15 */ \
\
	xmm6 -= xmm7;		/* I1 = I1 - I2 (final I2)		; 14-16	avail none */ \
	xmm14 *= xmm4;		/* B4 = I4 * cosine/sine		; 14-18 */ \
\
	xmm3 += xmm8;		/* R2 = R1 + R2 (final R1)		; 15-17	avail 8 */ \
	xmm15 *= xmm1;		/* B3 = I3 * cosine/sine		; 15-19 */ \
	xmm8 = xptr(screg+32);	/* Sine * normalization value		; 15 */ \
\
	xmm7 += xmm9;		/* I2 = I1 + I2 (final I1)		; 16-18	avail 9 */ \
	xmm10 *= xmm2;		/* A2 = R2 * cosine/sine		; 16-20 */ \
	xmm9 = xptr(screg+80);	/* Sine * normalization value		; 16 */ \
\
	xmm12 += xmm4;		/* A4 = A4 + I4				; 17-19	avail 4 */ \
	xmm11 *= xmm6;		/* B2 = I2 * cosine/sine		; 17-21 */ \
	xmm4 = xptr(srcreg+32);	/*#2 R1					; 17 */ \
\
	xmm13 -= xmm1;		/* A3 = A3 - I3				; 18-20	avail 1 */ \
	xmm3 *= xptr(screg);		/* Apply normalization to final R1	; 18-22 */ \
	xmm1 = xptr(srcreg+d2+32);	/* next R3 */ \
\
	xprefetchw(u8ptr(srcreg+srcinc+d1)); \
\
	xmm14 -= xmm5;		/* B4 = B4 - R4				; 19-21	avail 5 */ \
	xmm7 *= xptr(screg);		/* Apply normalization to final I1	; 19-23 */ \
	xptr(srcreg) = xmm3;		/* Save R1				; 23	avail 3,5 */ \
	xmm3 = xptr(srcreg+d1+48);	/*#2 I2					; 19 */ \
\
	xmm15 += xmm0;		/* B3 = B3 + R3				; 20-22	avail 5,0 */ \
	xmm12 *= xmm8;		/* A4 = A4 * sine (final R4)		; 20-24 */ \
	xmm5 = xptr(srcreg+d2+d1+48);	/*#2 I4 */ \
\
	xmm10 -= xmm6;		/* A2 = A2 - I2				; 21-23	avail 0,6 */ \
	xmm13 *= xmm8;		/* A3 = A3 * sine (final R3)		; 21-25 */ \
	xptr(srcreg+16) = xmm7;	/* Save I1				; 24	avail 0,6,7 */ \
	xmm7 = xptr(srcreg+d1+32);	/*#2 R2 */ \
\
	xmm11 += xmm2;		/* B2 = B2 + R2				; 22-24	avail 0,6,2 */ \
	xmm14 *= xmm8;		/* B4 = B4 * sine (final I4)		; 22-26 */ \
	xmm0 = xptr(srcreg+d2+d1+32);	/*#2 R4 */ \
\
	xmm6 = xmm4;		/*#2 Copy R1				; 23-25	avail 2 */ \
	xmm4 -= xmm1;		/*#2 R1 = R1 - R3 (new R3)		; 23-25 */ \
	xmm15 *= xmm8;		/* B3 = B3 * sine (final I3)		; 23-27	avail 2,8 */ \
\
	xmm2 = xmm3;		/*#2 Copy I2				; 24-26	avail 8 */ \
	xmm3 -= xmm5;		/*#2 I2 = I2 - I4 (new I4)		; 24-26 */ \
	xmm10 *= xmm9;		/* A2 = A2 * sine (final R2)		; 24-28 */ \
\
	xmm8 = xmm7;		/*#2 Copy R2				; 25-27	avail none */ \
	xmm7 -= xmm0;		/*#2 R2 = R2 - R4 (new R4)		; 25-27 */ \
	xmm11 *= xmm9;		/* B2 = B2 * sine (final I2)		; 25-29	avail 9 */ \
	xmm9 = xptr(srcreg+48);	/*#2 I1					; 25 */ \
	xptr(srcreg+d1+32) = xmm12;	/* Save R4				; 25	avail 12 */ \
\
	xmm12 = xptr(srcreg+d2+48);	/*#2 I3					; 26	avail none */ \
	xptr(srcreg+d1) = xmm13;	/* Save R3				; 26	avail 13 */ \
	xmm13 = xmm9;		/*#2 Copy I1				; 26-28	avail none */ \
	xmm9 -= xmm12;		/*#2 I1 = I1 - I3 (new I3)		; 26-28 */ \
\
	xmm1 += xmm6;		/*#2 R3 = R1 + R3 (new R1)		; 27-29	avail 6 */ \
	xptr(srcreg+d1+48) = xmm14;	/* Save I4				; 27	avail 6,14 */ \
	xmm14 = xptr(screg+32+16);	/*#2 cosine/sine			; 27	avail 6 */ \
\
	xmm0 += xmm8;		/*#2 R4 = R2 + R4 (new R2)		; 28-30	avail 6,8 */ \
	xptr(srcreg+d1+16) = xmm15;	/* Save I3				; 28	avail 6,8,15 */ \
	xmm15 = xmm14;							/* 28-30	avail 6,8 */ \
\
	xmm5 += xmm2;		/*#2 I4 = I2 + I4 (new I2)		; 29-31	avail 6,8,2 */ \
	xptr(srcreg+32) = xmm10;	/* Save R2				; 29	avail 6,8,2,10 */ \
	xmm6 = xmm3;		/*#2 Copy I4				; 29-31	avail 8,2,10 */ \
\
	xmm12 += xmm13;		/*#2 I3 = I1 + I3 (new I1)		; 30-32	avail 8,2,10,13 */ \
	xptr(srcreg+48) = xmm11;	/* Save I2				; 30	avail 8,2,10,13,11 */ \
	xmm10 = xmm14;							/* 30-32	avail 8,2,13,11 */ \
\
	xprefetchw(u8ptr(srcreg+srcinc+d2)); \
\
	xmm3 += xmm4;		/*#2 I4 = R3 + I4 (final R4)		; 31-33 */ \
	xmm11 = xmm14;							/* 31-33	avail 8,2,13 */ \
	xmm2 = xmm9;		/*#2 Copy I3				; 31-33	avail 8,13 */ \
\
	xmm4 -= xmm6;		/*#2 R3 = R3 - I4 (final R3)		; 32-34	avail 8,13,6 */ \
	xmm8 = xptr(screg+80+16);	/*#2 cosine/sine			; 32	avail 13,6 */ \
\
	xmm9 -= xmm7;		/*#2 I3 = I3 - R4 (final I4)		; 33-35 */ \
	xmm13 = xmm8;							/* 33-35	avail 6 */ \
\
	xmm7 += xmm2;		/*#2 R4 = I3 + R4 (final I3)		; 34-36	avail 6,2 */ \
	xmm14 *= xmm3;		/*#2 A4 = R4 * cosine/sine		; 34-38 */ \
	xmm6 = xmm1;		/*#2 Copy R1				; 34-36	avail 2 */ \
\
	xmm1 -= xmm0;		/*#2 R1 = R1 - R2 (final R2)		; 35-37 */ \
	xmm15 *= xmm4;		/*#2 A3 = R3 * cosine/sine		; 35-39 */ \
	xmm2 = xmm12;		/*#2 Copy I1				; 35-37	avail none */ \
\
	xmm12 -= xmm5;		/*#2 I1 = I1 - I2 (final I2)		; 36-38 */ \
	xmm10 *= xmm9;		/*#2 B4 = I4 * cosine/sine		; 36-40 */ \
\
	xmm0 += xmm6;		/*#2 R2 = R1 + R2 (final R1)		; 37-39	avail 6 */ \
	xmm11 *= xmm7;		/*#2 B3 = I3 * cosine/sine		; 37-41 */ \
	xmm6 = xptr(screg+32);	/*#2 Sine * normalization value		; 37	avail none */ \
\
	xmm5 += xmm2;		/*#2 I2 = I1 + I2 (final I1)		; 38-40	avail 2 */ \
	xmm8 *= xmm1;		/*#2 A2 = R2 * cosine/sine		; 38-42 */ \
	xmm2 = xptr(screg+80);	/*#2 Sine * normalization value		; 38	avail none */ \
\
	xmm14 += xmm9;		/*#2 A4 = A4 + I4			; 39-41	avail 9 */ \
	xmm13 *= xmm12;		/*#2 B2 = I2 * cosine/sine		; 39-43 */ \
\
	xmm15 -= xmm7;		/*#2 A3 = A3 - I3			; 40-42	avail 9,7 */ \
	xmm0 *= xptr(screg);		/*#2 Apply normalization to final R1	; 40-44 */ \
\
	xprefetchw(u8ptr(srcreg+srcinc+d2+d1)); \
\
	xmm10 -= xmm3;		/*#2 B4 = B4 - R4			; 41-43	avail 9,7,3 */ \
	xmm5 *= xptr(screg);		/*#2 Apply normalization to final I1	; 41-45 */ \
\
	xmm11 += xmm4;		/*#2 B3 = B3 + R3			; 42-44	avail 9,7,3,4 */ \
	xmm14 *= xmm6;		/*#2 A4 = A4 * sine (final R4)		; 42-46 */ \
\
	xmm8 -= xmm12;		/*#2 A2 = A2 - I2			; 43-45	avail 9,7,3,4,12 */ \
	xmm15 *= xmm6;		/*#2 A3 = A3 * sine (final R3)		; 43-47 */ \
\
	xmm13 += xmm1;		/*#2 B2 = B2 + R2			; 44-46	avail 9,7,3,4,12,1 */ \
	xmm10 *= xmm6;		/*#2 B4 = B4 * sine (final I4)		; 44-48 */ \
\
	xmm11 *= xmm6;		/*#2 B3 = B3 * sine (final I3)		; 45-49 */ \
	xptr(srcreg+d2) = xmm0;	/*#2 Save R1				; 45 */ \
\
	xmm8 *= xmm2;		/*#2 A2 = A2 * sine (final R2)		; 46-50 */ \
	xptr(srcreg+d2+16) = xmm5;	/*#2 Save I1				; 46 */ \
\
	xmm13 *= xmm2;		/*#2 B2 = B2 * sine (final I2)		; 47-51 */ \
	xptr(srcreg+d2+d1+32) = xmm14; /*#2 Save R4				; 47 */ \
\
	xptr(srcreg+d2+d1) = xmm15;	/*#2 Save R3				; 48 */ \
	xptr(srcreg+d2+d1+48) = xmm10; /*#2 Save I4				; 49 */ \
	xptr(srcreg+d2+d1+16) = xmm11; /*#2 Save I3				; 50 */ \
	xptr(srcreg+d2+32) = xmm8;	/*#2 Save R2				; 51 */ \
	xptr(srcreg+d2+48) = xmm13;	/*#2 Save I2				; 52 */ \
\
	srcreg += srcinc; \

#endif
#endif
\
/* 64-bit AMD K8 optimized versions of the above macros */ \

#if 0 // (@INSTR(,%xarch,<K8>) NE 0);
#if 0 //DEF X86_64;
\

#define r4_x4c_wpn_djbfft_mem_preload \
	xmm15 = XMM_TWO; \
\
\

#define r4_x4c_wpn_djbfft_mem(R1,R2,R3,R4,R5,R6,R7,R8,screg1,screg2,normreg,pre1,pre2,dst1,dst2) \
	xmm0 = R1;		/* R1 */ \
	xmm2 = R3;		/* R3 */ \
	xmm2 += xmm0;		/* R3 = R1 + R3 (new R1) */ \
\
	xmm1 = R2;		/* R2 */ \
	xmm3 = R4;		/* R4 */ \
	xmm3 += xmm1;		/* R4 = R2 + R4 (new R2) */ \
\
	xmm4 = R5;		/* I1 */ \
	xmm6 = R7;		/* I3 */ \
	xmm6 += xmm4;		/* I3 = I1 + I3 (new I1) */ \
\
	xmm0 -= R3;		/* R1 = R1 - R3 (new R3) */ \
	xmm1 -= R4;		/* R2 = R2 - R4 (new R4) */ \
\
	xprefetchw(u8ptr(pre1)); \
\
	xmm3 += xmm2;		/* R2 = R1 + R2 (final R1) */ \
	xmm2 *= xmm15; \
	xmm2 -= xmm3;		/* R1 = R1 - R2 (final R2) */ \
\
	xmm5 = R6;		/* I2 */ \
	xmm7 = R8;		/* I4 */ \
	xmm7 += xmm5;		/* I4 = I2 + I4 (new I2) */ \
\
	xmm4 -= R7;		/* I1 = I1 - I3 (new I3) */ \
	xmm5 -= R8;		/* I2 = I2 - I4 (new I4) */ \
\
	xmm3 *= xptr(normreg);		/* Apply normalization multiplier to R1 */ \
	dst1 = xmm3;		/* Save R1 */ \
\
	xmm7 += xmm6;		/* I2 = I1 + I2 (final I1) */ \
	xmm6 *= xmm15; \
	xmm6 -= xmm7;		/* I1 = I1 - I2 (final I2) */ \
\
	xmm7 *= xptr(normreg);		/* Apply normalization multiplier to I1 */ \
	dst2 = xmm7;		/* Save I1 */ \
\
	xmm1 += xmm4;		/* R4 = I3 + R4 (final I3) */ \
	xmm4 *= xmm15; \
\
	xmm0 -= xmm5;		/* R3 = R3 - I4 (final R3) */ \
	xmm5 *= xmm15; \
	xmm4 -= xmm1;		/* I3 = I3 - R4 (final I4) */ \
	xmm5 += xmm0;		/* I4 = R3 + I4 (final R4) */ \
\
	xprefetchw(u8ptr(pre1+pre2)); \
\
	xmm7 = xptr(screg2+16);	/* cosine/sine */ \
	xmm7 *= xmm2;		/* A2 = R2 * cosine/sine */ \
	xmm7 -= xmm6;		/* A2 = A2 - I2 */ \
\
	xmm3 = xptr(screg1+16);	/* cosine/sine */ \
	xmm3 *= xmm0;		/* A3 = R3 * cosine/sine */ \
	xmm3 -= xmm1;		/* A3 = A3 - I3 */ \
\
	xmm6 *= xptr(screg2+16);	/* B2 = I2 * cosine/sine */ \
	xmm6 += xmm2;		/* B2 = B2 + R2 */ \
\
	xmm2 = xptr(screg1+16);	/* cosine/sine */ \
	xmm2 *= xmm5;		/* A4 = R4 * cosine/sine */ \
	xmm2 += xmm4;		/* A4 = A4 + I4 */ \
\
	xmm1 *= xptr(screg1+16);	/* B3 = I3 * cosine/sine */ \
	xmm1 += xmm0;		/* B3 = B3 + R3 */ \
\
	xmm4 *= xptr(screg1+16);	/* B4 = I4 * cosine/sine */ \
	xmm4 -= xmm5;		/* B4 = B4 - R4 */ \
\
	xmm7 *= xptr(screg2);		/* A2 = A2 * sine (final R2) */ \
	xmm6 *= xptr(screg2);		/* B2 = B2 * sine (final I2) */ \
\
	xmm5 = xptr(screg1);		/* Sine */ \
	xmm3 *= xmm5;		/* A3 = A3 * sine (final R3) */ \
	xmm1 *= xmm5;		/* B3 = B3 * sine (final I3) */ \
\
	xmm2 *= xmm5;		/* A4 = A4 * sine (final R4) */ \
	xmm4 *= xmm5;		/* B4 = B4 * sine (final I4) */ \
\
\

#define r4_x4c_wpn_djbfft_partial_mem_preload \
	xmm15 = XMM_TWO; \
\
\

#define r4_x4c_wpn_djbfft_partial_mem(r1,r2,r3,r4,r5,r6,r7,r8,mem3,mem4,mem7,mem8,screg1,screg2,normreg,pre1,pre2,dst1,dst2) \
	r3 = mem3; \
	r3 += r1;			/* R3 = R1 + R3 (new R1) */ \
	r1 -= mem3;		/* R1 = R1 - R3 (new R3) */ \
\
	r4 = mem4; \
	r4 += r2;			/* R4 = R2 + R4 (new R2) */ \
	r2 -= mem4;		/* R2 = R2 - R4 (new R4) */ \
\
	r3 -= r4;			/* R1 = R1 - R2 (final R2) */ \
	r4 *= xmm15; \
	r4 += r3;			/* R2 = R1 + R2 (final R1) */ \
\
	r4 *= xptr(normreg);		/* Apply normalization multiplier to R1 */ \
	dst1 = r4;		/* Save R1 */ \
\
	r7 = mem7; \
	r8 = mem8; \
\
	r6 -= r8;			/* I2 = I2 - I4 (new I4) */ \
	r8 *= xmm15; \
	r8 += r6;			/* I4 = I2 + I4 (new I2) */ \
\
	r5 -= r7;			/* I1 = I1 - I3 (new I3) */ \
	r7 *= xmm15; \
	r7 += r5;			/* I3 = I1 + I3 (new I1) */ \
\
	r1 -= r6;			/* R3 = R3 - I4 (final R3) */ \
	r6 *= xmm15; \
	r6 += r1;			/* I4 = R3 + I4 (final R4) */ \
\
	r2 += r5;			/* R4 = I3 + R4 (final I3) */ \
	r5 *= xmm15; \
	r5 -= r2;			/* I3 = I3 - R4 (final I4) */ \
\
	xprefetchw(u8ptr(pre1)); \
\
	r7 -= r8;			/* I1 = I1 - I2 (final I2) */ \
	r8 *= xmm15; \
	r8 += r7;			/* I2 = I1 + I2 (final I1) */ \
\
	r4 = xptr(screg1+16);		/* cosine/sine */ \
	r4 *= r1;			/* A3 = R3 * cosine/sine */ \
\
	r8 *= xptr(normreg);		/* Apply normalization multiplier to I1 */ \
	dst2 = r8;		/* Save I1 */ \
\
	r8 = xptr(screg2+16);		/* cosine/sine */ \
	r8 *= r3;			/* A2 = R2 * cosine/sine */ \
\
	r4 -= r2;			/* A3 = A3 - I3 */ \
	r2 *= xptr(screg1+16);		/* B3 = I3 * cosine/sine */ \
	r2 += r1;			/* B3 = B3 + R3 */ \
\
	r1 = xptr(screg1+16);		/* cosine/sine */ \
	r1 *= r6;			/* A4 = R4 * cosine/sine */ \
\
	r8 -= r7;			/* A2 = A2 - I2 */ \
	r7 *= xptr(screg2+16);		/* B2 = I2 * cosine/sine */ \
\
	xprefetchw(u8ptr(pre1+pre2)); \
\
	r1 += r5;			/* A4 = A4 + I4 */ \
	r5 *= xptr(screg1+16);		/* B4 = I4 * cosine/sine */ \
\
	r7 += r3;			/* B2 = B2 + R2 */ \
	r5 -= r6;			/* B4 = B4 - R4 */ \
\
	r3 = xptr(screg1);		/* sine */ \
	r4 *= r3;			/* A3 = A3 * sine (final R3) */ \
	r2 *= r3;			/* B3 = B3 * sine (final I3) */ \
	r8 *= xptr(screg2);		/* A2 = A2 * sine (final R2) */ \
	r1 *= r3;			/* A4 = A4 * sine (final R4) */ \
	r7 *= xptr(screg2);		/* B2 = B2 * sine (final I2) */ \
	r5 *= r3;			/* B4 = B4 * sine (final I4) */ \
\

#endif
#endif
\

#define r4_x4cl_wpn_four_complex_djbunfft_preload \
	r4_x4c_wpn_djbunfft_mem_preload; \
;	/*r4_x4c_wpn_djbunfft_partial_mem_preload -- assume same as r4_x4c_wpn_djbunfft_mem_preload */ \
\
\

#define r4_x4cl_wpn_four_complex_djbunfft(srcreg,srcinc,d1,d2,screg_a,scoff) { \
	vec2f64 *sd0 = (vec2f64*)(srcreg); \
	vec2f64 *sd1 = (vec2f64*)(srcreg+d1); \
	vec2f64 *sd2 = (vec2f64*)(srcreg+d2); \
	vec2f64 *sd3 = (vec2f64*)(srcreg+d2+d1); \
	char *screg = (char*)(screg_a); /* improves clang compiler */ \
	vec2f64 xma0,xma1,xma2,xma3,xma4,xma5,xma6,xma7; \
\
	r4_x4c_wpn_djbunfft_mem(sd0[1],sd0[3],sd1[1],sd1[3],sd2[1],sd2[3],sd3[1],sd3[3], sd1[1], screg+32+scoff,screg+scoff,srcreg+srcinc+d2,d1); \
/*	xstore	[srcreg+d1+16], xmm6	;; Save R3 */ \
	xma2 = sd1[0];	/* R3 */ \
	sd1[0] = xmm7;	/* Save R1 */ \
	xma3 = sd1[2];	/* R4 */ \
	sd1[2] = xmm4;	/* Save I1 */ \
	sd1[3] = xmm3;	/* Save I3 */ \
	xma6 = sd3[0];	/* R7 */ \
	xma7 = sd3[2];	/* R8 */ \
	sd3[0] = xmm5;	/* Save R2 */ \
	sd3[1] = xmm0;	/* Save R4 */ \
	sd3[2] = xmm2;	/* Save I2 */ \
	sd3[3] = xmm1;	/* Save I4 */ \
	r4_x4c_wpn_djbunfft_partial_mem(xma0,xma1,xma2,xma3,xma4,xma5,xma6,xma7, sd0[0],sd0[2],sd2[0],sd2[2],sd0[0], screg+32,screg,srcreg+srcinc,d1); \
/*	xstore	[srcreg], xma1		;; Save R1 */ \
	sd0[1] = xma2;	/* Save R3 */ \
	sd0[2] = xma4;	/* Save I1 */ \
	sd0[3] = xma0;	/* Save I3 */ \
	sd2[0] = xma3;	/* Save R2 */ \
	sd2[1] = xma5;	/* Save R4 */ \
	sd2[2] = xma6;	/* Save I2 */ \
	sd2[3] = xma7;	/* Save I4 */ \
	srcreg += srcinc; \
}


#define r4_x4c_wpn_djbunfft_mem_preload \
\
\

#define r4_x4c_wpn_djbunfft_mem(mem1,mem2,mem3m,mem4,mem5m,mem6,mem7m,mem8,dst3,screg,normreg,pre1,pre2) { \
	vec2f64 mem3 = mem3m, mem5 = mem5m, mem7 = mem7m; \
	xmm3 = xptr(screg+48+16);	/* cosine/sine */ \
	xmm7 = mem3;		/* R2 */ \
	xmm0 = mem4;		/* I2 */ \
	xmm6 = xmm7 * xmm3 + xmm0;	/* A2 = R2 * cosine/sine + I2 */ \
	xmm3 = xmm3 * xmm0 - xmm7;	/* B2 = I2 * cosine/sine - R2 */ \
\
	xmm4 = xptr(screg+16);	/* cosine/sine */ \
	xmm2 = mem5;		/* R3 */ \
	xmm1 = mem6;		/* I3 */ \
	xmm7 = xmm2 * xmm4 + xmm1;	/* A3 = R3 * cosine/sine + I3 */ \
	xmm5 = xmm4 * xmm1 - xmm2;		/* B3 = I3 * cosine/sine - R3 */ \
\
	xmm1 = mem7;		/* R4 */ \
	xmm0 = mem8;		/* I4 */ \
	xmm2 = xmm1 * xmm4 - xmm0;	/* A4 = R4 * cosine/sine - I4 */ \
	xmm4 = xmm4 * xmm0 + xmm1;	/* B4 = I4 * cosine/sine + R4 */ \
\
	xmm0 = mem1;		/* R1 */ \
	xmm0 *= xptr(normreg+16);	/* R1 * normalization_inverse */ \
\
	xmm1 = xptr(screg+48+32); \
	xmm6 *= xmm1;		/* A2 = A2 * sine*normalization_inverse (new R2) */ \
	xmm3 *= xmm1;		/* B2 = B2 * sine*normalization_inverse (new I2) */ \
	xmm1 = xptr(screg+32);	/* Sine * normalization_inverse */ \
	xmm7 *= xmm1;		/* A3 = A3 * sine*normalization_inverse (new R3) */ \
	xmm5 *= xmm1;		/* B3 = B3 * sine*normalization_inverse (new I3) */ \
	xmm2 *= xmm1;		/* A4 = A4 * sine*normalization_inverse (new R4) */ \
	xmm4 *= xmm1;		/* B4 = B4 * sine*normalization_inverse (new I4) */ \
\
	xprefetchw(u8ptr(pre1)); \
\
	xmm1 = xmm0;		/* Copy R1 */ \
	xmm0 -= xmm6;		/* R1 = R1 - R2 (new R2) */ \
	xmm6 += xmm1;		/* R2 = R1 + R2 (new R1) */ \
\
	xmm1 = xmm2;		/* Copy R4 */ \
	xmm2 -= xmm7;		/* R4 = R4 - R3 (new I4) */ \
	xmm7 += xmm1;		/* R3 = R4 + R3 (new R3) */ \
\
	xmm1 = xmm5;		/* Copy I3 */ \
	xmm5 -= xmm4;		/* I3 = I3 - I4 (new R4) */ \
	xmm4 += xmm1;		/* I4 = I3 + I4 (new I3) */ \
\
	xmm1 = xmm6;		/* Copy R1 */ \
	xmm6 -= xmm7;		/* R1 = R1 - R3 (final R3) */ \
	xmm7 += xmm1;		/* R3 = R1 + R3 (final R1) */ \
\
	xmm1 = mem2;		/* I1 */ \
	xmm1 *= xptr(normreg+16);	/* I1 * normalization_inverse */ \
\
	dst3 = xmm6; \
\
	xmm6 = xmm1;		/* Copy I1 */ \
	xmm1 -= xmm3;		/* I1 = I1 - I2 (new I2) */ \
	xmm3 += xmm6;		/* I2 = I1 + I2 (new I1) */ \
\
	xprefetchw(u8ptr(pre1+pre2)); \
\
	xmm6 = xmm0;		/* Copy R2 */ \
	xmm0 -= xmm5;		/* R2 = R2 - R4 (final R4) */ \
	xmm5 += xmm6;		/* R4 = R2 + R4 (final R2) */ \
\
	xmm6 = xmm1;		/* Copy I2 */ \
	xmm1 -= xmm2;		/* I2 = I2 - I4 (final I4) */ \
	xmm2 += xmm6;		/* I4 = I2 + I4 (final I2) */ \
\
	xmm6 = xmm3;		/* Copy I1 */ \
	xmm3 -= xmm4;		/* I1 = I1 - I3 (final I3) */ \
	xmm4 += xmm6;		/* I3 = I1 + I3 (final I1) */ \
}
\

#define r4_x4c_wpn_djbunfft_partial_mem_preload \
\
\

#define r4_x4c_wpn_djbunfft_partial_mem(r1,r2,r3,r4,r5,r6,r7,r8,mem1,mem2,mem5m,mem6m,dst1,screg,normreg,pre1,pre2) { \
	r1 = xptr(screg+48+16);	/* cosine/sine */ \
	r2 = r3;			/* Copy R2 */ \
	r3 = r3 * r1 + r4;	/* A2 = R2 * cosine/sine + I2 */ \
	r1 = r1 * r4 - r2;	/* B2 = I2 * cosine/sine - R2 */ \
\
	r4 = xptr(screg+16);	/* cosine/sine */ \
	r6 = r7;			/* Copy R4 */ \
	r7 = r7 * r4 - r8;	/* A4 = R4 * cosine/sine - I4 */ \
	r5 = r4 * r8 + r6;	/* B4 = I4 * cosine/sine + R4 */ \
\
	vec2f64 mem5 = mem5m, mem6 = mem6m; \
	r2 = r4 * mem5 + mem6;	/* A3 = R3 * cosine/sine + I3 */ \
	r4 = r4 * mem6 - mem5;	/* B3 = I3 * cosine/sine - R3 */ \
\
	r6 = mem1;			/* R1 */ \
	r6 *= xptr(normreg+16);	/* R1 * normalization_inverse */ \
\
	r8 = xptr(screg+48+32); \
	r3 *= r8;			/* A2 = A2 * sine*normalization_inverse (new R2) */ \
	r1 *= r8;			/* B2 = B2 * sine*normalization_inverse (new I2) */ \
	r8 = xptr(screg+32);		/* Sine * normalization_inverse */ \
	r7 *= r8;			/* A4 = A4 * sine*normalization_inverse (new R4) */ \
	r5 *= r8;			/* B4 = B4 * sine*normalization_inverse (new I4) */ \
	r2 *= r8;			/* A3 = A3 * sine*normalization_inverse (new R3) */ \
	r4 *= r8;			/* B3 = B3 * sine*normalization_inverse (new I3) */ \
\
	xprefetchw(u8ptr(pre1)); \
\
	r8 = r6;			/* Copy R1 */ \
	r6 -= r3;			/* R1 = R1 - R2 (new R2) */ \
	r3 += r8;			/* R2 = R1 + R2 (new R1) */ \
\
	r8 = r7;			/* Copy R4 */ \
	r7 -= r2;			/* R4 = R4 - R3 (new I4) */ \
	r2 += r8;			/* R3 = R4 + R3 (new R3) */ \
\
	r8 = r4;			/* Copy I3 */ \
	r4 -= r5;			/* I3 = I3 - I4 (new R4) */ \
	r5 += r8;			/* I4 = I3 + I4 (new I3) */ \
\
	r8 = r3;			/* Copy R1 */ \
	r3 -= r2;			/* R1 = R1 - R3 (final R3) */ \
	r2 += r8;			/* R3 = R1 + R3 (final R1) */ \
\
	r8 = mem2;			/* I1 */ \
	r8 *= xptr(normreg+16);	/* I1 * normalization_inverse */ \
\
	dst1 = r2; \
\
	r2 = r8;			/* Copy I1 */ \
	r8 -= r1;			/* I1 = I1 - I2 (new I2) */ \
	r1 += r2;			/* I2 = I1 + I2 (new I1) */ \
\
	xprefetchw(u8ptr(pre1+pre2)); \
\
	r2 = r6;			/* Copy R2 */ \
	r6 -= r4;			/* R2 = R2 - R4 (final R4) */ \
	r4 += r2;			/* R4 = R2 + R4 (final R2) */ \
\
	r2 = r8;			/* Copy I2 */ \
	r8 -= r7;			/* I2 = I2 - I4 (final I4) */ \
	r7 += r2;			/* I4 = I2 + I4 (final I2) */ \
\
	r2 = r1;			/* Copy I1 */ \
	r1 -= r5;			/* I1 = I1 - I3 (final I3) */ \
	r5 += r2;			/* I3 = I1 + I3 (final I1) */ \
}

/* 32-bit AMD K8 optimized versions of the above macros */ \

#if 0 // (@INSTR(,%xarch,<K8>) NE 0);
\
/* K8 Cheat sheet for scheduling dependency chains */ \
/*	      12345678901234567890123456789012345678901234567890123456789012345678901234567890 */ \
/*A2	      MMMMAAAA    MMMM */ \
/*B2	        MMMMAAAA    MMMM */ \
/*A4              MMMMAAAA */ \
/*A3                MMMMAAAA */ \
/*B3		      MMMMAAAA */ \
/*B4		        MMMMAAAA */ \
\
/*i1			          MMMM */ \
/*r1		                    MMMM */ \
\
/*i4(r4-r3)                   AAAA        MMMM */ \
/*r3(r4+r3)                   MMMMAAAA        MMMM */ \
/*r4(i3-i4)                     AAAA        MMMM */ \
/*i3(i3+i4)                     MMMMAAAA        MMMM */ \
/*i2(i1-i2)                           AAAA */ \
/*i1(i1+i2)                           MMMMAAAA */ \
/*r2(r1-r2)                             AAAA */ \
/*r1(r1+r2)                             MMMMAAAA */ \
\
/*I4(I2-I4)	                              AAAA */ \
/*I2(I2+I4)			                  MMMMAAAA */ \
/*R4(R2-R4)	                                AAAA			 */ \
/*R2(R2+R4)			                    MMMMAAAA */ \
/*R3(R1-R3)	                                  AAAA */ \
/*R1(R1+R3)	                                      MMMMAAAA */ \
/*I3(I1-I3)	                                    AAAA				 */ \
/*I1(I1+I3)                                             MMMMAAAA */ \
\

#define r4_x4c_wpn_djbunfft_mem(mem1,mem2,mem3,mem4,mem5,mem6,mem7,mem8,dst3,screg,normreg,pre1,pre2) \
	xmm3 = xptr(screg+48+16);	/* cosine/sine */ \
	xmm6 = mem3;		/* R2 */ \
	xmm6 *= xmm3;		/* A2 = R2 * cosine/sine	; 1-4 */ \
\
	xmm0 = mem4;		/* I2 */ \
	xmm3 *= xmm0;		/* B2 = I2 * cosine/sine	; 3-6 */ \
\
	xmm4 = xptr(screg+16);	/* cosine/sine */ \
	xmm2 = mem7;		/* R4 */ \
	xmm2 *= xmm4;		/* A4 = R4 * cosine/sine	; 5-8 */ \
	xmm6 += xmm0;		/* A2 = A2 + I2			; 5-8 */ \
\
	xmm5 = xptr(screg+16);	/* cosine/sine */ \
	xmm7 = mem5;		/* R3 */ \
	xmm7 *= xmm5;		/* A3 = R3 * cosine/sine	; 7-10 */ \
	xmm3 -= mem3;		/* B2 = B2 - R2			; 7-10 */ \
\
	xmm1 = mem6;		/* I3 */ \
	xmm5 *= xmm1;		/* B3 = I3 * cosine/sine	; 9-12 */ \
	xmm0 = mem8;		/* I4 */ \
	xmm2 -= xmm0;		/* A4 = A4 - I4			; 9-12 */ \
\
	xmm4 *= xmm0;		/* B4 = I4 * cosine/sine	; 11-14 */ \
	xmm7 += xmm1;		/* A3 = A3 + I3			; 11-14 */ \
\
	xmm5 -= mem5;		/* B3 = B3 - R3			; 13-16 */ \
	xmm6 *= xptr(screg+48+32);	/* A2 = A2 * sine*normalization_inverse (new R2) ; 13-16 */ \
\
	xmm4 += mem7;		/* B4 = B4 + R4			; 15-18 */ \
	xmm3 *= xptr(screg+48+32);	/* B2 = B2 * sine*normalization_inverse (new I2) ; 15-18 */ \
\
	xmm2 -= xmm7;		/* R4 = R4 - R3 (new I4)	; 17-20 */ \
	xmm7 += xmm7; /* multwo */;							/* 17-20 */ \
\
	xprefetchw(u8ptr(pre1)); \
\
	xmm5 -= xmm4;		/* I3 = I3 - I4 (new R4)	; 19-22 */ \
	xmm4 += xmm4; /* multwo */;							/* 19-22 */ \
\
	xmm1 = mem2;		/* I1 */ \
	xmm7 += xmm2;		/* R3 = R4 + R3 (new R3)	; 21-24 */ \
	xmm1 *= xptr(normreg+16);	/* I1 * normalization_inverse	; 21-24 */ \
\
	xmm0 = mem1;		/* R1 */ \
	xmm4 += xmm5;		/* I4 = I3 + I4 (new I3)	; 23-26 */ \
	xmm0 *= xptr(normreg+16);	/* R1 * normalization_inverse	; 23-26 */ \
\
	xmm1 -= xmm3;		/* I1 = I1 - I2 (new I2)	; 25-28 */ \
	xmm3 += xmm3; /* multwo */;							/* 25-28 */ \
\
	xmm0 -= xmm6;		/* R1 = R1 - R2 (new R2)	; 27-30 */ \
	xmm6 += xmm6; /* multwo */;							/* 27-30 */ \
\
	xprefetchw(u8ptr(pre1+pre2)); \
\
	xmm3 += xmm1;		/* I2 = I1 + I2 (new I1)	; 29-32 */ \
	xmm2 *= xptr(screg+32);	/* I4 = I4 * sine*normalization_inverse */ \
\
	xmm6 += xmm0;		/* R2 = R1 + R2 (new R1)	; 31-34 */ \
	xmm5 *= xptr(screg+32);	/* R4 = R4 * sine*normalization_inverse */ \
\
	xmm1 -= xmm2;		/* I2 = I2 - I4 (final I4)	; 33-36 */ \
	xmm7 *= xptr(screg+32);	/* R3 = R3 * sine*normalization_inverse */ \
\
	xmm0 -= xmm5;		/* R2 = R2 - R4 (final R4)	; 35-38 */ \
	xmm4 *= xptr(screg+32);	/* I3 = I3 * sine*normalization_inverse */ \
\
	xmm6 -= xmm7;		/* R1 = R1 - R3 (final R3)	; 37-40 */ \
	xmm2 += xmm2; /* multwo */; \
\
	xmm3 -= xmm4;		/* I1 = I1 - I3 (final I3)	; 39-42 */ \
	xmm5 += xmm5; /* multwo */; \
\
	xmm2 += xmm1;		/* I4 = I2 + I4 (final I2)	; 41-44 */ \
	xmm7 += xmm7; /* multwo */; \
\
	xmm5 += xmm0;		/* R4 = R2 + R4 (final R2)	; 43-46 */ \
	xmm4 += xmm4; /* multwo */; \
\
	xmm7 += xmm6;		/* R3 = R1 + R3 (final R1)	; 45-48 */ \
	dst3 = xmm6;						/* 41 */ \
\
	xmm4 += xmm3;		/* I3 = I1 + I3 (final I1)	; 47-50 */ \
\
\

#define r4_x4c_wpn_djbunfft_partial_mem(r1,r2,r3,r4,r5,r6,r7,r8,mem1,mem2,mem5,mem6,dst1,screg,normreg,pre1,pre2) \
	r1 = xptr(screg+48+16);	/* cosine/sine */ \
	r1 *= r4;			/* B2 = I2 * cosine/sine	; 1-4 */ \
\
	r1 -= r3;			/* B2 = B2 - R2			; 5-8 */ \
	r3 *= xptr(screg+48+16);	/* A2 = R2 * cosine/sine	; 3-6 */ \
	r5 = xptr(screg+16);		/* cosine/sine */ \
\
	r2 = mem5;		/* R3 */ \
	r5 *= r8;			/* B4 = I4 * cosine/sine	; 5-8 */ \
\
	r3 += r4;			/* A2 = A2 + I2			; 7-10 */ \
	r4 = xptr(screg+16);		/* cosine/sine */ \
	r2 *= r4;			/* A3 = R3 * cosine/sine	; 7-10 */ \
\
	r5 += r7;			/* B4 = B4 + R4			; 9-12 */ \
	r7 *= xptr(screg+16);		/* A4 = R4 * cosine/sine	; 9-12 */ \
\
	r2 += mem6;		/* A3 = A3 + I3			; 11-14 */ \
	r4 *= mem6;		/* B3 = I3 * cosine/sine	; 11-14 */ \
\
	r7 -= r8;			/* A4 = A4 - I4			; 13-16 */ \
	r1 *= xptr(screg+48+32);	/* B2 = B2 * sine*normalization_inverse (new I2) */ \
\
	r4 -= mem5;		/* B3 = B3 - R3			; 15-18 */ \
	r3 *= xptr(screg+48+32);	/* A2 = A2 * sine*normalization_inverse (new R2) */ \
\
	r7 -= r2;			/* R4 = R4 - R3 (new I4)	; 17-20 */ \
	r2 += r2; /* multwo */							/* 17-20 */ \
\
	xprefetchw(u8ptr(pre1)); \
\
	r4 -= r5;			/* I3 = I3 - I4 (new R4)	; 19-22 */ \
	r5 += r5; /* multwo */							/* 19-22 */ \
\
	r8 = mem2;		/* I1 */ \
	r2 += r7;			/* R3 = R4 + R3 (new R3)	; 21-24 */ \
	r8 *= xptr(normreg+16);	/* I1 * normalization_inverse	; 21-24 */ \
\
	r6 = mem1;		/* R1 */ \
	r5 += r4;			/* I4 = I3 + I4 (new I3)	; 23-26 */ \
	r6 *= xptr(normreg+16);	/* R1 * normalization_inverse	; 23-26 */ \
\
	r8 -= r1;			/* I1 = I1 - I2 (new I2)	; 25-28 */ \
	r1 += r1; /* multwo */							/* 25-28 */ \
\
	xprefetchw(u8ptr(pre1+pre2)); \
\
	r6 -= r3;			/* R1 = R1 - R2 (new R2)	; 27-30 */ \
	r3 += r3; /* multwo */							/* 27-30 */ \
\
	r1 += r8;			/* I2 = I1 + I2 (new I1)	; 29-32 */ \
	r7 *= xptr(screg+32);		/* I4 = I4 * sine*normalization_inverse */ \
\
	r3 += r6;			/* R2 = R1 + R2 (new R1)	; 31-34 */ \
	r4 *= xptr(screg+32);		/* R4 = R4 * sine*normalization_inverse */ \
\
	r8 -= r7;			/* I2 = I2 - I4 (final I4)	; 33-36 */ \
	r2 *= xptr(screg+32);		/* R3 = R3 * sine*normalization_inverse */ \
\
	r6 -= r4;			/* R2 = R2 - R4 (final R4)	; 35-38 */ \
	r5 *= xptr(screg+32);		/* I3 = I3 * sine*normalization_inverse */ \
\
	r3 -= r2;			/* R1 = R1 - R3 (final R3)	; 37-40 */ \
	r7 += r7; /* multwo */ \
\
	r1 -= r5;			/* I1 = I1 - I3 (final I3)	; 39-42 */ \
	r4 += r4; /* multwo */ \
\
	r7 += r8;			/* I4 = I2 + I4 (final I2)	; 41-44 */ \
	r2 += r2; /* multwo */ \
\
	r4 += r6;			/* R4 = R2 + R4 (final R2)	; 43-46 */ \
	r5 += r5; /* multwo */ \
\
	r2 += r3;			/* R3 = R1 + R3 (final R1)	; 45-48 */ \
\
	r5 += r1;			/* I3 = I1 + I3 (final I1)	; 47-50 */ \
\
	dst1 = r2; \
\

#endif
\
/* 64-bit Intel and AMD K10 implementations of the above - use the extra XMM registers */ \

#if 0 // (@INSTR(,%xarch,<CORE>) NE 0) OR (@INSTR(,%xarch,<P4>) NE 0) OR (@INSTR(,%xarch,<BLEND>) NE 0) OR (@INSTR(,%xarch,<K10>) NE 0);
#if 0 //DEF X86_64;
\
/* Theoretical best case is 44 clocks on a Core 2.  Now at ?? clocks. */ \

#define r4_x4cl_wpn_four_complex_djbunfft(srcreg,srcinc,d1,d2,screg,scoff) \
	xmm0 = xptr(srcreg+d2+d1+16);	/* R4 */ \
	xmm1 = xptr(screg+scoff+32+16); /* cosine/sine */ \
	xmm2 = xmm0;		/* Copy R4				; 1-3 */ \
	xmm0 *= xmm1;		/* A4 = R4 * cosine/sine		; 1-5 */ \
\
	xmm5 = xptr(srcreg+d1+48);	/* I2 */ \
	xmm6 = xptr(screg+scoff+80+16); /* cosine/sine */ \
	xmm7 = xmm5;		/* Copy I2				; 2-4 */ \
	xmm5 *= xmm6;		/* B2 = I2 * cosine/sine		; 2-6 */ \
\
	xmm3 = xptr(srcreg+d2+16);	/* R3 */ \
	xmm4 = xmm3;		/* Copy R3				; 3-5 */ \
	xmm3 *= xmm1;		/* A3 = R3 * cosine/sine		; 3-7 */ \
\
	xmm8 = xptr(srcreg+d1+16);	/* R2 */ \
	xmm9 = xmm8;		/* Copy R2				; 4-6 */ \
	xmm8 *= xmm6;		/* A2 = R2 * cosine/sine		; 4-8	avail 6,10+ */ \
\
	xmm6 = xptr(srcreg+d2+48);	/* I3 */ \
	xmm10 = xmm6;		/* Copy I3				; 5-7	avail 11+ */ \
	xmm6 *= xmm1;		/* B3 = I3 * cosine/sine		; 5-9 */ \
\
	xmm11 = xptr(srcreg+d2+d1+48);/* I4 */ \
	xmm0 -= xmm11;		/* A4 = A4 - I4 (new R4 / sine)		; 6-8 */ \
	xmm11 *= xmm1;		/* B4 = I4 * cosine/sine		; 6-10	avail 1,12+ */ \
\
	xmm5 -= xmm9;		/* B2 = B2 - R2				; 7-9	avail 1,9,12+ */ \
	xmm1 = xptr(screg+scoff+80+32); /* sine * normalization_inverse	; 7 */ \
\
	xmm3 += xmm10;		/* A3 = A3 + I3 (new R3 / sine)		; 8-10	avail 9,10,12+ */ \
\
	xmm8 += xmm7;		/* A2 = A2 + I2				; 9-11	avail 9,10,7,12+ */ \
	xmm7 = xptr(srcreg+48);	/* I1					; 9 */ \
	xmm7 *= xptr(screg+scoff+16);	/* I1 * normalization_inverse		; 9-13 */ \
\
	xmm6 -= xmm4;		/* B3 = B3 - R3 (new I3 / sine)		; 10-12	avail 9,10,4,12+ */ \
	xmm5 *= xmm1;		/* B2 = B2 * sine (new I2)		; 10-14 */ \
	xmm4 = xmm0;		/* Copy R4 / sine			; 10-12	avail 9,10,12+ */ \
\
	xmm11 += xmm2;		/* B4 = B4 + R4 (new I4 / sine)		; 11-13	avail 9,10,2,12+ */ \
	xmm2 = xptr(srcreg+16);	/* R1					; 11	avail 9,10,12+ */ \
	xmm2 *= xptr(screg+scoff+16);	/* R1 * normalization_inverse		; 11-15 */ \
\
	xmm0 -= xmm3;		/* R4 = R4 - R3 (newer I4 / sine)	; 12-14 */ \
	xmm8 *= xmm1;		/* A2 = A2 * sine (new R2)		; 12-16	avail 9,10,1,12+ */ \
	xmm1 = xptr(screg+scoff+32+32); /* sine * normalization_inverse	; 12	avail 9,10,12+ */ \
\
	xmm3 += xmm4;		/* R3 = R4 + R3 (newer R3 / sine)	; 13-15	avail 9,10,4,12+ */ \
	xmm4 = xmm6;		/* Copy I3 / sine			; 13-15 */ \
\
	xprefetchw(u8ptr(srcreg+srcinc+d2)); \
\
	xmm6 -= xmm11;		/* I3 = I3 - I4 (newer R4 / sine)	; 14-16 */ \
	xmm10 = xmm7;		/* Copy I1				; 14-16 avail 9,12+ */ \
	xmm9 = xptr(srcreg+d2+d1);	/*#2 R4					; 14	avail 12+ */ \
\
	xmm7 -= xmm5;		/* I1 = I1 - I2 (newer I2)		; 15-17 */ \
	xmm0 *= xmm1;		/* newer I4 * sine			; 15-19 */ \
	xmm12 = xptr(screg+32+16);	/*#2 cosine/sine			; 15	avail 13+ */ \
\
	xmm11 += xmm4;		/* I4 = I3 + I4 (newer I3 / sine)	; 16-18 avail 4,13+ */ \
	xmm3 *= xmm1;		/* newer R3 * sine			; 16-20 */ \
	xmm4 = xmm2;		/* Copy R1				; 16-18	avail 13+ */ \
\
	xmm5 += xmm10;		/* I2 = I1 + I2 (newer I1)		; 17-19	avail 10,13+ */ \
	xmm6 *= xmm1;		/* newer R4 * sine			; 17-21 */ \
	xmm13 = xptr(srcreg+d1+32);	/*#2 I2					; 17	avail 10,14+ */ \
\
	xmm2 += xmm8;		/* R2 = R1 + R2 (newer R1)		; 18-20 */ \
	xmm14 = xptr(screg+80+16);	/*#2 cosine/sine			; 18	avail 10,15 */ \
\
	xmm4 -= xmm8;		/* R1 = R1 - R2 (newer R2)		; 19-21	avail 10,8,15 */ \
	xmm11 *= xmm1;		/* newer I3 * sine			; 19-23	avail 10,8,1,15 */ \
	xmm1 = xmm7;		/* Copy I2				; 19-21	avail 10,8,15 */ \
\
	xprefetchw(u8ptr(srcreg+srcinc+d2+d1)); \
\
	xmm7 -= xmm0;		/* I2 = I2 - I4 (final I4)		; 20-22	avail 10,8,15 storable 7 */ \
	xmm15 = xmm9;		/*#2 Copy R4				; 20-22	avail 10,8 storable 7 */ \
\
	xmm8 = xmm2;		/* Copy R1				; 21-23	avail 10 storable 7 */ \
	xmm2 -= xmm3;		/* R1 = R1 - R3 (final R3)		; 21-23	avail 10 storable 7,2 */ \
\
	xmm0 += xmm1;		/* I4 = I2 + I4 (final I2)		; 22-24	avail 10,1 storable 7,2,0 */ \
	xmm1 = xmm4;		/* Copy R2				; 22-24	avail 10 storable 7,2,0 */ \
\
	xmm4 -= xmm6;		/* R2 = R2 - R4 (final R4)		; 23-25	avail 10 storable 7,2,0,4 */ \
	xmm9 *= xmm12;		/*#2 A4 = R4 * cosine/sine		; 23-27 */ \
	xmm10 = xmm13;		/*#2 Copy I2				; 23-25	avail none storable 7,2,0,4 */ \
	xptr(srcreg+d2+d1+48) = xmm7;	/* Save I4				; 23	avail 7 storable 2,0,4 */ \
\
	xmm3 += xmm8;		/* R3 = R1 + R3 (final R1)		; 24-26	avail 7,8 storable 2,0,4,3 */ \
	xmm13 *= xmm14;		/*#2 B2 = I2 * cosine/sine		; 24-28 */ \
	xmm7 = xmm5;		/* Copy I1				; 24-26	avail 8 storable 2,0,4,3 */ \
	xptr(srcreg+d1+16) = xmm2;	/* Save R3				; 24	avail 8,2 storable 0,4,3 */ \
\
	xmm6 += xmm1;		/* R4 = R2 + R4 (final R2)		; 25-27	avail 8,2,1 storable 0,4,3,6 */ \
	xmm2 = xptr(srcreg+d2);	/*#2 R3					; 25	avail 8,1 storable 0,4,3,6 */ \
	xmm8 = xmm2;		/*#2 Copy R3				; 25-27	avail 1 storable 0,4,3,6 */ \
	xmm2 *= xmm12;		/*#2 A3 = R3 * cosine/sine		; 25-29 */ \
\
	xmm5 -= xmm11;		/* I1 = I1 - I3 (final I3)		; 26-28	avail 1 storable 0,4,3,6,5 */ \
	xmm1 = xptr(srcreg+d1);	/*#2 R2					; 26	avail none storable 0,4,3,6,5 */ \
	xptr(srcreg+d2+d1+16) = xmm4;	/* Save R4				; 26	avail 4 storable 0,3,6,5 */ \
	xmm4 = xmm1;		/*#2 Copy R2				; 26-28	avail none storable 0,3,6,5 */ \
	xmm1 *= xmm14;		/*#2 A2 = R2 * cosine/sine		; 26-30	avail 14 storable 0,3,6,5 */ \
\
	xmm11 += xmm7;		/* I3 = I1 + I3 (final I1)		; 27-29	avail 14,7 storable 0,3,6,5,11 */ \
	xmm14 = xptr(srcreg+d2+32);	/*#2 I3					; 27	avail 7 storable 0,3,6,5,11 */ \
	xmm7 = xmm14;		/*#2 Copy I3				; 27-29	avail none storable 0,3,6,5,11 */ \
	xmm14 *= xmm12;		/*#2 B3 = I3 * cosine/sine		; 27-31 */ \
	xptr(srcreg+d1) = xmm3;	/* Save R1				; 27	avail 3 storable 0,6,5,11 */ \
\
	xmm3 = xptr(srcreg+d2+d1+32); /*#2 I4					; 28	avail none storable 0,6,5,11 */ \
	xmm9 -= xmm3;		/*#2 A4 = A4 - I4 (new R4 / sine)	; 28-30 */ \
	xmm3 *= xmm12;		/*#2 B4 = I4 * cosine/sine		; 28-32	avail 12 storable 0,6,5,11 */ \
	xptr(srcreg+d2+d1) = xmm6;	/* Save R2				; 28	avail 12,6 storable 0,5,11 */ \
\
	xmm13 -= xmm4;		/*#2 B2 = B2 - R2			; 29-31	avail 12,6,4 storable 0,5,11 */ \
	xmm6 = xptr(screg+80+32);	/*#2 sine * normalization_inverse	; 29	avail 12,4 storable 0,5,11 */ \
	xptr(srcreg+d1+48) = xmm5;	/* Save I3				; 29	avail 12,4,5 storable 0,11 */ \
\
	xmm2 += xmm7;		/*#2 A3 = A3 + I3 (new R3 / sine)	; 30-32	avail 12,4,5,7 storable 0,11 */ \
	xmm5 = xptr(srcreg+32);	/*#2 I1					; 30	avail 12,4,7 storable 0,11 */ \
	xmm5 *= xptr(screg+16);	/*#2 I1 * normalization_inverse		; 30-34 */ \
	xptr(srcreg+d1+32) = xmm11;	/* Save I1				; 30	avail 12,4,7,11 storable 0 */ \
\
	xmm1 += xmm10;		/*#2 A2 = A2 + I2			; 31-33	avail 12,4,7,11,10 storable 0 */ \
	xptr(srcreg+d2+d1+32) = xmm0;	/* Save I2				; 25	avail 12,4,7,11,10,0 */ \
	xmm7 = xptr(srcreg);		/*#2 R1					; 31	avail 12,4,11,10,0 */ \
\
	xmm14 -= xmm8;		/*#2 B3 = B3 - R3 (new I3 / sine)	; 32-34	avail 12,4,11,10,0,8 */ \
	xmm13 *= xmm6;		/*#2 B2 = B2 * sine (new I2)		; 32-36 */ \
	xmm4 = xmm9;		/*#2 Copy R4 / sine			; 32-34	avail 12,11,10,0,8 */ \
\
	xmm3 += xmm15;		/*#2 B4 = B4 + R4 (new I4 / sine)	; 33-35	avail 12,11,10,0,8,15 */ \
	xmm7 *= xptr(screg+16);	/*#2 R1 * normalization_inverse		; 33-37 */ \
\
	xmm9 -= xmm2;		/*#2 R4 = R4 - R3 (newer I4 / sine)	; 34-36 */ \
	xmm1 *= xmm6;		/*#2 A2 = A2 * sine (new R2)		; 34-38	avail 12,11,10,0,8,15,6 */ \
	xmm6 = xptr(screg+32++32);	/*#2 sine * normalization_inverse	; 34	avail 12,11,10,0,8,15 */ \
\
	xmm2 += xmm4;		/*#2 R3 = R4 + R3 (newer R3 / sine)	; 35-37	avail 12,11,10,0,8,15,4 */ \
	xmm4 = xmm14;		/*#2 Copy I3 / sine			; 35-37	avail 12,11,10,0,8,15 */ \
\
	xprefetchw(u8ptr(srcreg+srcinc)); \
\
	xmm14 -= xmm3;		/*#2 I3 = I3 - I4 (newer R4 / sine)	; 36-38 */ \
	xmm0 = xmm5;		/*#2 Copy I1				; 36-38	avail 12,11,10,8,15 */ \
\
	xmm5 -= xmm13;		/*#2 I1 = I1 - I2 (newer I2)		; 37-39 */ \
	xmm9 *= xmm6;		/*#2 newer I4 * sine			; 37-41 */ \
\
	xmm3 += xmm4;		/*#2 I4 = I3 + I4 (newer I3 / sine)	; 38-40	avail 12,11,10,8,15,0 */ \
	xmm2 *= xmm6;		/*#2 newer R3 * sine			; 38-42 */ \
	xmm4 = xmm7;		/*#2 Copy R1				; 38-40	avail 12,11,10,8,15 */ \
\
	xmm13 += xmm0;		/*#2 I2 = I1 + I2 (newer I1)		; 39-41	avail 12,11,10,8,15,0 */ \
	xmm14 *= xmm6;		/*#2 newer R4 * sine			; 39-43 */ \
\
	xmm7 += xmm1;		/*#2 R2 = R1 + R2 (newer R1)		; 40-42 */ \
	xmm0 = xmm5;		/*#2 Copy I2				; 40-42	avail 12,11,10,8,15 */ \
\
	xmm4 -= xmm1;		/*#2 R1 = R1 - R2 (newer R2)		; 41-43	avail 12,11,10,8,15,1 */ \
	xmm3 *= xmm6;		/*#2 newer I3 * sine			; 41-45	avail 12,11,10,8,15,1,6 */ \
\
	xprefetchw(u8ptr(srcreg+srcinc+d1)); \
\
	xmm5 -= xmm9;		/*#2 I2 = I2 - I4 (final I4)		; 42-44	avail 12,11,10,8,15,6 storable 5 */ \
	xmm1 = xmm13;		/*#2 Copy I1				; 42-44	avail 12,11,10,8,15,6 storable 5 */ \
\
	xmm6 = xmm7;		/*#2 Copy R1				; 43-45	avail 12,11,10,8,15 storable 5 */ \
	xmm7 -= xmm2;		/*#2 R1 = R1 - R3 (final R3)		; 43-45	avail 12,11,10,8,15 storable 5,7 */ \
\
	xmm9 += xmm0;		/*#2 I4 = I2 + I4 (final I2)		; 44-46	avail 12,11,10,8,15,0 storable 5,7,9 */ \
	xmm0 = xmm4;		/*#2 Copy R2				; 44-46	avail 12,11,10,8,15 storable 5,7,9 */ \
\
	xmm4 -= xmm14;		/*#2 R2 = R2 - R4 (final R4)		; 45-47	avail 12,11,10,8,15 storable 5,7,9,4 */ \
	xptr(srcreg+d2+48) = xmm5;	/*#2 Save I4				; 45	avail 12,11,10,8,15,5 storable 7,9,4 */ \
\
	xmm2 += xmm6;		/*#2 R3 = R1 + R3 (final R1)		; 46-48	avail 12,11,10,8,15,5,6 storable 7,9,4,2 */ \
	xptr(srcreg+16) = xmm7;	/*#2 Save R3				; 46	avail 12,11,10,8,15,5,6,7 storable 9,4,2 */ \
\
	xmm14 += xmm0;		/*#2 R4 = R2 + R4 (final R2)		; 47-49	avail 12,11,10,8,15,5,6,7,0 storable 9,4,2,14 */ \
	xptr(srcreg+d2+32) = xmm9;	/*#2 Save I2				; 47	avail 12,11,10,8,15,5,6,7,0,9 storable 4,2,14 */ \
\
	xmm13 -= xmm3;		/*#2 I1 = I1 - I3 (final I3)		; 48-50	avail 12,11,10,8,15,5,6,7,0,9 storable 4,2,14,13 */ \
	xptr(srcreg+d2+16) = xmm4;	/*#2 Save R4				; 48	avail 12,11,10,8,15,5,6,7,0,9,4 storable 2,14,13 */ \
\
	xmm3 += xmm1;		/*#2 I3 = I1 + I3 (final I1)		; 49-51	avail 12,11,10,8,15,5,6,7,0,9,4,1 storable 2,14,13,3 */ \
	xptr(srcreg) = xmm2;		/*#2 Save R1				; 49	avail 12,11,10,8,15,5,6,7,0,9,4,1,2 storable 14,13,3 */ \
\
	xptr(srcreg+d2) = xmm14;	/*#2 Save R2				; 50 */ \
	xptr(srcreg+48) = xmm13;	/*#2 Save I3				; 51 */ \
	xptr(srcreg+32) = xmm3;	/*#2 Save I1				; 52 */ \
\
	srcreg += srcinc; \

#endif
#endif
\
/* 64-bit AMD K8 optimized versions of the above macros */ \

#if 0 // (@INSTR(,%xarch,<K8>) NE 0);
#if 0 //DEF X86_64;
\

#define r4_x4c_wpn_djbunfft_mem_preload \
	xmm15 = XMM_TWO; \
\
\

#define r4_x4c_wpn_djbunfft_mem(mem1,mem2,mem3,mem4,mem5,mem6,mem7,mem8,dst3,screg,normreg,pre1,pre2) \
	xmm3 = xptr(screg+48+16);	/* cosine/sine */ \
	xmm6 = mem3;		/* R2 */ \
	xmm6 *= xmm3;		/* A2 = R2 * cosine/sine	; 1-4 */ \
\
	xmm0 = mem4;		/* I2 */ \
	xmm3 *= xmm0;		/* B2 = I2 * cosine/sine	; 3-6 */ \
\
	xmm4 = xptr(screg+16);	/* cosine/sine */ \
	xmm2 = mem7;		/* R4 */ \
	xmm2 *= xmm4;		/* A4 = R4 * cosine/sine	; 5-8 */ \
	xmm6 += xmm0;		/* A2 = A2 + I2			; 5-8 */ \
\
	xmm5 = xptr(screg+16);	/* cosine/sine */ \
	xmm7 = mem5;		/* R3 */ \
	xmm7 *= xmm5;		/* A3 = R3 * cosine/sine	; 7-10 */ \
	xmm3 -= mem3;		/* B2 = B2 - R2			; 7-10 */ \
\
	xmm1 = mem6;		/* I3 */ \
	xmm5 *= xmm1;		/* B3 = I3 * cosine/sine	; 9-12 */ \
	xmm0 = mem8;		/* I4 */ \
	xmm2 -= xmm0;		/* A4 = A4 - I4			; 9-12 */ \
\
	xmm4 *= xmm0;		/* B4 = I4 * cosine/sine	; 11-14 */ \
	xmm7 += xmm1;		/* A3 = A3 + I3			; 11-14 */ \
\
	xmm5 -= mem5;		/* B3 = B3 - R3			; 13-16 */ \
	xmm6 *= xptr(screg+48+32);	/* A2 = A2 * sine*normalization_inverse (new R2) ; 13-16 */ \
\
	xmm4 += mem7;		/* B4 = B4 + R4			; 15-18 */ \
	xmm3 *= xptr(screg+48+32);	/* B2 = B2 * sine*normalization_inverse (new I2) ; 15-18 */ \
\
	xmm2 -= xmm7;		/* R4 = R4 - R3 (new I4)	; 17-20 */ \
	xmm7 *= xmm15;						/* 17-20 */ \
\
	xprefetchw(u8ptr(pre1)); \
\
	xmm5 -= xmm4;		/* I3 = I3 - I4 (new R4)	; 19-22 */ \
	xmm4 *= xmm15;						/* 19-22 */ \
\
	xmm1 = mem2;		/* I1 */ \
	xmm7 += xmm2;		/* R3 = R4 + R3 (new R3)	; 21-24 */ \
	xmm1 *= xptr(normreg+16);	/* I1 * normalization_inverse	; 21-24 */ \
\
	xmm0 = mem1;		/* R1 */ \
	xmm4 += xmm5;		/* I4 = I3 + I4 (new I3)	; 23-26 */ \
	xmm0 *= xptr(normreg+16);	/* R1 * normalization_inverse	; 23-26 */ \
\
	xmm1 -= xmm3;		/* I1 = I1 - I2 (new I2)	; 25-28 */ \
	xmm3 *= xmm15;						/* 25-28 */ \
\
	xmm0 -= xmm6;		/* R1 = R1 - R2 (new R2)	; 27-30 */ \
	xmm6 *= xmm15;						/* 27-30 */ \
\
	xprefetchw(u8ptr(pre1+pre2)); \
\
	xmm3 += xmm1;		/* I2 = I1 + I2 (new I1)	; 29-32 */ \
	xmm2 *= xptr(screg+32);	/* I4 = I4 * sine*normalization_inverse */ \
\
	xmm6 += xmm0;		/* R2 = R1 + R2 (new R1)	; 31-34 */ \
	xmm5 *= xptr(screg+32);	/* R4 = R4 * sine*normalization_inverse */ \
\
	xmm1 -= xmm2;		/* I2 = I2 - I4 (final I4)	; 33-36 */ \
	xmm7 *= xptr(screg+32);	/* R3 = R3 * sine*normalization_inverse */ \
\
	xmm0 -= xmm5;		/* R2 = R2 - R4 (final R4)	; 35-38 */ \
	xmm4 *= xptr(screg+32);	/* I3 = I3 * sine*normalization_inverse */ \
\
	xmm6 -= xmm7;		/* R1 = R1 - R3 (final R3)	; 37-40 */ \
	xmm2 *= xmm15; \
\
	xmm3 -= xmm4;		/* I1 = I1 - I3 (final I3)	; 39-42 */ \
	xmm5 *= xmm15; \
\
	xmm2 += xmm1;		/* I4 = I2 + I4 (final I2)	; 41-44 */ \
	xmm7 *= xmm15; \
\
	xmm5 += xmm0;		/* R4 = R2 + R4 (final R2)	; 43-46 */ \
	xmm4 *= xmm15; \
\
	xmm7 += xmm6;		/* R3 = R1 + R3 (final R1)	; 45-48 */ \
	dst3 = xmm6;						/* 41 */ \
\
	xmm4 += xmm3;		/* I3 = I1 + I3 (final I1)	; 47-50 */ \
\
\

#define r4_x4c_wpn_djbunfft_partial_mem_preload \
	xmm15 = XMM_TWO; \
\
\

#define r4_x4c_wpn_djbunfft_partial_mem(r1,r2,r3,r4,r5,r6,r7,r8,mem1,mem2,mem5,mem6,dst1,screg,normreg,pre1,pre2) \
	r1 = xptr(screg+48+16);	/* cosine/sine */ \
	r1 *= r4;			/* B2 = I2 * cosine/sine	; 1-4 */ \
\
	r1 -= r3;			/* B2 = B2 - R2			; 5-8 */ \
	r3 *= xptr(screg+48+16);	/* A2 = R2 * cosine/sine	; 3-6 */ \
	r5 = xptr(screg+16);		/* cosine/sine */ \
\
	r2 = mem5;		/* R3 */ \
	r5 *= r8;			/* B4 = I4 * cosine/sine	; 5-8 */ \
\
	r3 += r4;			/* A2 = A2 + I2			; 7-10 */ \
	r4 = xptr(screg+16);		/* cosine/sine */ \
	r2 *= r4;			/* A3 = R3 * cosine/sine	; 7-10 */ \
\
	r5 += r7;			/* B4 = B4 + R4			; 9-12 */ \
	r7 *= xptr(screg+16);		/* A4 = R4 * cosine/sine	; 9-12 */ \
\
	r2 += mem6;		/* A3 = A3 + I3			; 11-14 */ \
	r4 *= mem6;		/* B3 = I3 * cosine/sine	; 11-14 */ \
\
	r7 -= r8;			/* A4 = A4 - I4			; 13-16 */ \
	r1 *= xptr(screg+48+32);	/* B2 = B2 * sine*normalization_inverse (new I2) */ \
\
	r4 -= mem5;		/* B3 = B3 - R3			; 15-18 */ \
	r3 *= xptr(screg+48+32);	/* A2 = A2 * sine*normalization_inverse (new R2) */ \
\
	r7 -= r2;			/* R4 = R4 - R3 (new I4)	; 17-20 */ \
	r2 *= xmm15;						/* 17-20 */ \
\
	xprefetchw(u8ptr(pre1)); \
\
	r4 -= r5;			/* I3 = I3 - I4 (new R4)	; 19-22 */ \
	r5 *= xmm15;						/* 19-22 */ \
\
	r8 = mem2;		/* I1 */ \
	r2 += r7;			/* R3 = R4 + R3 (new R3)	; 21-24 */ \
	r8 *= xptr(normreg+16);	/* I1 * normalization_inverse	; 21-24 */ \
\
	r6 = mem1;		/* R1 */ \
	r5 += r4;			/* I4 = I3 + I4 (new I3)	; 23-26 */ \
	r6 *= xptr(normreg+16);	/* R1 * normalization_inverse	; 23-26 */ \
\
	r8 -= r1;			/* I1 = I1 - I2 (new I2)	; 25-28 */ \
	r1 *= xmm15;						/* 25-28 */ \
\
	xprefetchw(u8ptr(pre1+pre2)); \
\
	r6 -= r3;			/* R1 = R1 - R2 (new R2)	; 27-30 */ \
	r3 *= xmm15;						/* 27-30 */ \
\
	r1 += r8;			/* I2 = I1 + I2 (new I1)	; 29-32 */ \
	r7 *= xptr(screg+32);		/* I4 = I4 * sine*normalization_inverse */ \
\
	r3 += r6;			/* R2 = R1 + R2 (new R1)	; 31-34 */ \
	r4 *= xptr(screg+32);		/* R4 = R4 * sine*normalization_inverse */ \
\
	r8 -= r7;			/* I2 = I2 - I4 (final I4)	; 33-36 */ \
	r2 *= xptr(screg+32);		/* R3 = R3 * sine*normalization_inverse */ \
\
	r6 -= r4;			/* R2 = R2 - R4 (final R4)	; 35-38 */ \
	r5 *= xptr(screg+32);		/* I3 = I3 * sine*normalization_inverse */ \
\
	r3 -= r2;			/* R1 = R1 - R3 (final R3)	; 37-40 */ \
	r7 *= xmm15; \
\
	r1 -= r5;			/* I1 = I1 - I3 (final I3)	; 39-42 */ \
	r4 *= xmm15; \
\
	r7 += r8;			/* I4 = I2 + I4 (final I2)	; 41-44 */ \
	r2 *= xmm15; \
\
	r4 += r6;			/* R4 = R2 + R4 (final R2)	; 43-46 */ \
	r5 *= xmm15; \
\
	r2 += r3;			/* R3 = R1 + R3 (final R1)	; 45-48 */ \
\
	r5 += r1;			/* I3 = I1 + I3 (final I1)	; 47-50 */ \
\
	dst1 = r2; \
\

#endif
#endif
\
\
/* */ \
/* ************************************* four-complex-with_square and variants ****************************************** */ \
/* */ \
/* These macros are used in the last levels of pass 2 in two pass FFTs. */ \
/* */ \
\
/* */ \
/* The last two levels of the forward FFT are performed. */ \
/* No sin/cos multipliers are needed. */ \
/* */ \
\

#define r4_x4cl_four_complex_fft_final(srcreg,srcinc,d1,d2) { \
	vec2f64 *sd0 = (vec2f64*)(srcreg); \
	vec2f64 *sd1 = (vec2f64*)(srcreg+d1); \
	vec2f64 *sd2 = (vec2f64*)(srcreg+d2); \
	vec2f64 *sd3 = (vec2f64*)(srcreg+d2+d1); \
	r4_x4c_simple_fft_mem(sd0[0],sd1[0],sd2[0],sd3[0],sd0[1],sd1[1],sd2[1],sd3[1], srcreg+srcinc,d1); \
	sd0[0] = xmm3;	/* Save R1 */ \
	sd0[1] = xmm7;	/* Save R2 */ \
	xmm3 = sd0[2];	/* R1 */ \
	xmm7 = sd0[3];	/* R5 */ \
	sd0[2] = xmm2;	/* Save R3 */ \
	sd0[3] = xmm6;	/* Save R4 */ \
	xmm2 = sd1[2];	/* R2 */ \
	xmm6 = sd1[3];	/* R6 */ \
	sd1[0] = xmm0;	/* Save R5 */ \
	sd1[1] = xmm1;	/* Save R6 */ \
	sd1[2] = xmm5;	/* Save R7 */ \
	sd1[3] = xmm4;	/* Save R8 */ \
	r4_x4c_simple_fft_partial_mem(xmm3,xmm2,xmm0,xmm5,xmm7,xmm6,xmm1,xmm4, sd2[2],sd3[2],sd2[3],sd3[3],srcreg+srcinc+d2,d1); \
	sd2[0] = xmm5;	/* Save R1 */ \
	sd2[1] = xmm4;	/* Save R2 */ \
	sd2[2] = xmm0;	/* Save R3 */ \
	sd2[3] = xmm1;	/* Save R4 */ \
	sd3[0] = xmm3;	/* Save R5 */ \
	sd3[1] = xmm2;	/* Save R6 */ \
	sd3[2] = xmm6;	/* Save R7 */ \
	sd3[3] = xmm7;	/* Save R8 */ \
	srcreg += srcinc; \
}

/* */ \
/* The last two levels of the forward FFT are performed, point-wise */ \
/* squaring, and first two levels of the inverse FFT are performed. */ \
/* No sin/cos multipliers are needed. */ \
/* */ \
\

#define r4_x4cl_four_complex_with_square_preload \
\
\

#define r4_x4cl_four_complex_with_square(srcreg,srcinc,d1,d2) { \
	vec2f64 *sd0 = (vec2f64*)(srcreg); \
	vec2f64 *sd1 = (vec2f64*)(srcreg+d1); \
	vec2f64 *sd2 = (vec2f64*)(srcreg+d2); \
	vec2f64 *sd3 = (vec2f64*)(srcreg+d2+d1); \
	vec2f64 temp; \
	r4_x4c_simple_fft_mem(sd0[0],sd1[0],sd2[0],sd3[0],sd0[1],sd1[1],sd2[1],sd3[1], srcreg+srcinc,d1); \
	temp = xmm0;			/* R5 */ \
	xp_complex_square(xmm3, xmm7, xmm0);	/* Square R1, R2 */ \
	xp_complex_square(xmm2, xmm6, xmm0);	/* Square R3, R4 */ \
	xp_complex_square(xmm5, xmm4, xmm0);	/* Square R7, R8 */ \
	xmm0 = temp; \
	temp = xmm4; \
	xp_complex_square(xmm0, xmm1, xmm4);	/* Square R5, R6 */ \
	r4_x4c_simple_unfft(xmm3,xmm7,xmm2,xmm6,xmm0,xmm1,xmm5,xmm4, temp,,,sd0[0]); \
/*	xstore	[srcreg], xmm0			;; Save R1 */ \
	sd0[1] = xmm2;		/* Save R3 */ \
	xmm0 = sd0[2];		/* R1 */ \
	xmm2 = sd0[3];		/* R5 */ \
	sd0[2] = xmm4;		/* Save R5 */ \
	sd0[3] = xmm6;		/* Save R7 */ \
	xmm4 = sd1[2];		/* R2 */ \
	xmm6 = sd1[3];		/* R6 */ \
	sd1[0] = xmm1;		/* Save R2 */ \
	sd1[1] = xmm3;		/* Save R4 */ \
	sd1[2] = xmm5;		/* Save R6 */ \
	sd1[3] = xmm7;		/* Save R8 */ \
	r4_x4c_simple_fft_partial_mem(xmm0,xmm4,xmm1,xmm3,xmm2,xmm6,xmm5,xmm7, sd2[2],sd3[2],sd2[3],sd3[3],srcreg+srcinc+d2,d1); \
	temp = xmm0;		/* R5 */ \
	xp_complex_square(xmm3, xmm7, xmm0);	/* Square R1, R2 */ \
	xp_complex_square(xmm1, xmm5, xmm0);	/* Square R3, R4 */ \
	xp_complex_square(xmm6, xmm2, xmm0);	/* Square R7, R8 */ \
	xmm0 = temp;		/* R5 */ \
	temp = xmm2;		/* R8 */ \
	xp_complex_square(xmm0, xmm4, xmm2);	/* Square R5, R6 */ \
	r4_x4c_simple_unfft(xmm3,xmm7,xmm1,xmm5,xmm0,xmm4,xmm6,xmm2, temp,,,sd2[0]); \
/*	xstore	[srcreg+d2], xmm0		;; Save R1 */ \
	sd2[1] = xmm1;		/* Save R3 */ \
	sd2[2] = xmm2;		/* Save R5 */ \
	sd2[3] = xmm5;		/* Save R7 */ \
	sd3[0] = xmm4;		/* Save R2 */ \
	sd3[1] = xmm3;		/* Save R4 */ \
	sd3[2] = xmm6;		/* Save R6 */ \
	sd3[3] = xmm7;		/* Save R8 */ \
	srcreg += srcinc; \
}

/* */ \
/* The last two levels of the forward FFT are performed, point-wise */ \
/* multiplication, and first two levels of the inverse FFT are performed. */ \
/* No sin/cos multipliers are needed. */ \
/* */ \
\

#define r4_x4cl_four_complex_with_mult(srcreg,srcinc,d1,d2) { \
	uintptr_t src_rbp = srcreg + rbp; \
	vec2f64 *sd0 = (vec2f64*)(srcreg); \
	vec2f64 *sd1 = (vec2f64*)(srcreg+d1); \
	vec2f64 *sd2 = (vec2f64*)(srcreg+d2); \
	vec2f64 *sd3 = (vec2f64*)(srcreg+d2+d1); \
	vec2f64 *pd0 = (vec2f64*)(src_rbp); \
	vec2f64 *pd1 = (vec2f64*)(src_rbp+d1); \
	vec2f64 *pd2 = (vec2f64*)(src_rbp+d2); \
	vec2f64 *pd3 = (vec2f64*)(src_rbp+d2+d1); \
	vec2f64 temp; \
	r4_x4c_simple_fft_mem(sd0[0],sd1[0],sd2[0],sd3[0],sd0[1],sd1[1],sd2[1],sd3[1], srcreg+srcinc,d1); \
	xp4c_mulf_r(xmm3,xmm7,xmm2,xmm6,xmm0,xmm1,xmm5,xmm4, pd0[0],pd0[1],pd0[2],pd0[3],pd1[0],pd1[1],pd1[2],pd1[3]); \
	temp = xmm4; \
	r4_x4c_simple_unfft(xmm3,xmm7,xmm2,xmm6,xmm0,xmm1,xmm5,xmm4,temp,,,sd0[0]); \
/*	xstore	[srcreg], xmm0			;; Save R1 */ \
	sd0[1] = xmm2;		/* Save R3 */ \
	xmm0 = sd0[2];		/* R1 */ \
	xmm2 = sd0[3];		/* R5 */ \
	sd0[2] = xmm4;		/* Save R5 */ \
	sd0[3] = xmm6;		/* Save R7 */ \
	xmm4 = sd1[2];		/* R2 */ \
	xmm6 = sd1[3];		/* R6 */ \
	sd1[0] = xmm1;		/* Save R2 */ \
	sd1[1] = xmm3;		/* Save R4 */ \
	sd1[2] = xmm5;		/* Save R6 */ \
	sd1[3] = xmm7;		/* Save R8 */ \
	r4_x4c_simple_fft_partial_mem(xmm0,xmm4,xmm1,xmm3,xmm2,xmm6,xmm5,xmm7, sd2[2],sd3[2],sd2[3],sd3[3],srcreg+srcinc+d2,d1); \
	xp4c_mulf_r(xmm3,xmm7,xmm1,xmm5,xmm0,xmm4,xmm6,xmm2, pd2[0],pd2[1],pd2[2],pd2[3],pd3[0],pd3[1],pd3[2],pd3[3]); \
	temp = xmm2; \
	r4_x4c_simple_unfft(xmm3,xmm7,xmm1,xmm5,xmm0,xmm4,xmm6,xmm2,temp,,,sd2[0]); \
/*	xstore	[srcreg+d2], xmm0		;; Save R1 */ \
	sd2[1] = xmm1;		/* Save R3 */ \
	sd2[2] = xmm2;		/* Save R5 */ \
	sd2[3] = xmm5;		/* Save R7 */ \
	sd3[0] = xmm4;		/* Save R2 */ \
	sd3[1] = xmm3;		/* Save R4 */ \
	sd3[2] = xmm6;		/* Save R6 */ \
	sd3[3] = xmm7;		/* Save R8 */ \
	srcreg += srcinc; \
}

/* */ \
/* Point-wise multiplication and first two levels of the inverse FFT are performed. */ \
/* No sin/cos multipliers are needed. */ \
/* */ \
\

#define r4_x4cl_four_complex_with_mulf(srcreg,srcinc,d1,d2) { \
	uintptr_t src_rbx = srcreg + rbx; \
	uintptr_t src_rbp = srcreg + rbp; \
	vec2f64 *sd0 = (vec2f64*)(srcreg); \
	vec2f64 *sd1 = (vec2f64*)(srcreg+d1); \
	vec2f64 *sd2 = (vec2f64*)(srcreg+d2); \
	vec2f64 *sd3 = (vec2f64*)(srcreg+d2+d1); \
	vec2f64 *bd0 = (vec2f64*)(src_rbx); \
	vec2f64 *bd1 = (vec2f64*)(src_rbx+d1); \
	vec2f64 *bd2 = (vec2f64*)(src_rbx+d2); \
	vec2f64 *bd3 = (vec2f64*)(src_rbx+d2+d1); \
	vec2f64 *pd0 = (vec2f64*)(src_rbp); \
	vec2f64 *pd1 = (vec2f64*)(src_rbp+d1); \
	vec2f64 *pd2 = (vec2f64*)(src_rbp+d2); \
	vec2f64 *pd3 = (vec2f64*)(src_rbp+d2+d1); \
	vec2f64 temp; \
	xmm3 = bd0[0];		/* R1 */ \
	xmm7 = bd0[1];		/* R2 */ \
	xmm2 = bd0[2];		/* R3 */ \
	xmm6 = bd0[3];		/* R4 */ \
	xmm0 = bd1[0];		/* R5 */ \
	xmm1 = bd1[1];		/* R6 */ \
	xmm5 = bd1[2];		/* R7 */ \
	xmm4 = bd1[3];		/* R8 */ \
	xp4c_mulf_r(xmm3,xmm7,xmm2,xmm6,xmm0,xmm1,xmm5,xmm4, pd0[0],pd0[1],pd0[2],pd0[3],pd1[0],pd1[1],pd1[2],pd1[3]); \
	temp = xmm4; \
	r4_x4c_simple_unfft(xmm3,xmm7,xmm2,xmm6,xmm0,xmm1,xmm5,xmm4, temp,srcreg+srcinc,d1,sd0[0]); \
/*	xstore	[srcreg], xmm0			;; Save R1 */ \
	sd0[1] = xmm2;		/* Save R3 */ \
	sd0[2] = xmm4;		/* Save R5 */ \
	sd0[3] = xmm6;		/* Save R7 */ \
	sd1[0] = xmm1;		/* Save R2 */ \
	sd1[1] = xmm3;		/* Save R4 */ \
	sd1[2] = xmm5;		/* Save R6 */ \
	sd1[3] = xmm7;		/* Save R8 */ \
	xmm3 = bd2[0];		/* R1 */ \
	xmm7 = bd2[1];		/* R2 */ \
	xmm1 = bd2[2];		/* R3 */ \
	xmm5 = bd2[3];		/* R4 */ \
	xmm0 = bd3[0];		/* R5 */ \
	xmm4 = bd3[1];		/* R6 */ \
	xmm6 = bd3[2];		/* R7 */ \
	xmm2 = bd3[3];		/* R8 */ \
	xp4c_mulf_r(xmm3,xmm7,xmm1,xmm5,xmm0,xmm4,xmm6,xmm2, pd2[0],pd2[1],pd2[2],pd2[3],pd3[0],pd3[1],pd3[2],pd3[3]); \
	temp = xmm2; \
	r4_x4c_simple_unfft(xmm3,xmm7,xmm1,xmm5,xmm0,xmm4,xmm6,xmm2, temp,srcreg+srcinc+d2,d1,sd2[0]); \
/*	xstore	[srcreg+d2], xmm0		;; Save R1 */ \
	sd2[1] = xmm1;		/* Save R3 */ \
	sd2[2] = xmm2;		/* Save R5 */ \
	sd2[3] = xmm5;		/* Save R7 */ \
	sd3[0] = xmm4;		/* Save R2 */ \
	sd3[1] = xmm3;		/* Save R4 */ \
	sd3[2] = xmm6;		/* Save R6 */ \
	sd3[3] = xmm7;		/* Save R8 */ \
	srcreg += srcinc; \
}



#define r4_x4c_simple_fft_mem(R1,R2,R3,R4,R5,R6,R7,R8,pre1,pre2) { \
	vec2f64 temp2; \
	xmm0 = R1;			/* R1 */ \
	xmm2 = R3;			/* R3 */ \
	xmm5 = xmm0;		/* Copy R1 */ \
	xmm0 -= xmm2;		/* R1 = R1 - R3 (new R3) */ \
	xmm2 += xmm5;		/* R3 = R1 + R3 (new R1) */ \
\
	xmm4 = R5;			/* I1 */ \
	xmm6 = R7;			/* I3 */ \
	xmm5 = xmm4;		/* Copy I1 */ \
	xmm4 -= xmm6;		/* I1 = I1 - I3 (new I3) */ \
	xmm6 += xmm5;		/* I3 = I1 + I3 (new I1) */ \
\
	xprefetchw(u8ptr(pre1)); \
\
	xmm1 = R2;			/* R2 */ \
	xmm3 = R4;			/* R4 */ \
	xmm5 = xmm1;		/* Copy R2 */ \
	xmm1 -= xmm3;		/* R2 = R2 - R4 (new R4) */ \
	xmm3 += xmm5;		/* R4 = R2 + R4 (new R2) */ \
\
	xmm5 = R6;			/* I2 */ \
	xmm7 = R8;			/* I4 */ \
	xmm5 -= xmm7;		/* I2 = I2 - I4 (new I4) */ \
	xmm7 += R6;			/* I4 = I2 + I4 (new I2) */ \
\
	xprefetchw(u8ptr(pre1+pre2)); \
\
	temp2 = xmm4; \
	xmm4 -= xmm1;		/* I3 = I3 - R4 (final I4) */ \
	xmm1 += temp2;		/* R4 = I3 + R4 (final I3) */ \
\
	temp2 = xmm0; \
	xmm0 -= xmm5;		/* R3 = R3 - I4 (final R3) */ \
	xmm5 += temp2;		/* I4 = R3 + I4 (final R4) */ \
\
	temp2 = xmm2; \
	xmm2 -= xmm3;		/* R1 = R1 - R2 (final R2) */ \
	xmm3 += temp2;		/* R2 = R1 + R2 (final R1) */ \
\
	temp2 = xmm6; \
	xmm6 -= xmm7;		/* I1 = I1 - I2 (final I2) */ \
	xmm7 += temp2;		/* I2 = I1 + I2 (final I1) */ \
}


#define r4_x4c_simple_fft_partial_mem(r1, r2, r3, r4, r5, r6, r7, r8, mem3, mem4, mem7, mem8, pre1, pre2) { \
	vec2f64 temp; \
	r3 = mem3;			/* R3 */ \
	r8 = r3;			/* Copy R3 */ \
	r3 += r1;			/* R3 = R1 + R3 (new R1) */ \
	r1 -= r8;			/* R1 = R1 - R3 (new R3) */ \
\
	r7 = mem7;			/* I3 */ \
	r8 = r7;			/* Copy I3 */ \
	r7 += r5;			/* I3 = I1 + I3 (new I1) */ \
	r5 -= r8;			/* I1 = I1 - I3 (new I3) */ \
\
	xprefetchw(u8ptr(pre1)); \
\
	r4 = mem4;			/* R4 */ \
	r8 = r4;			/* Copy R4 */ \
	r4 += r2;			/* R4 = R2 + R4 (new R2) */ \
	r2 -= r8;			/* R2 = R2 - R4 (new R4) */ \
\
	vec2f64 m8 = mem8; \
	r8 = m8;			/* I4 */ \
	r8 += r6;			/* I4 = I2 + I4 (new I2) */ \
	r6 -= m8;			/* I2 = I2 - I4 (new I4) */ \
\
	temp = r5; \
	r5 -= r2;			/* I3 = I3 - R4 (final I4) */ \
	r2 += temp;			/* R4 = I3 + R4 (final I3) */ \
\
	xprefetchw(u8ptr(pre1+pre2)); \
\
	temp = r1; \
	r1 -= r6;			/* R3 = R3 - I4 (final R3) */ \
	r6 += temp;			/* I4 = R3 + I4 (final R4) */ \
\
	temp = r3; \
	r3 -= r4;			/* R1 = R1 - R2 (final R2) */ \
	r4 += temp;			/* R2 = R1 + R2 (final R1) */ \
\
	temp = r7; \
	r7 -= r8;			/* I1 = I1 - I2 (final I2) */ \
	r8 += temp;			/* I2 = I1 + I2 (final I1) */ \
}


#define r4_x4c_simple_unfft(r1, r2, r3, r4, r5, r6, r7, r8, mem8, pre1, pre2, dest1) { \
	r8 = r1;			/* R1 */ \
	r1 -= r3;			/* R1 = R1 - R2 (new R2) */ \
	r3 += r8;			/* R2 = R1 + R2 (new R1) */ \
\
	r8 = r7;			/* R4 */ \
	r7 -= r5;			/* R4 = R4 - R3 (new I4) */ \
	r5 += r8;			/* R3 = R4 + R3 (new R3) */ \
\
	r8 = r2;			/* I2 */ \
	r2 -= r4;			/* I1 = I1 - I2 (new I2) */ \
	r4 += r8;			/* I2 = I1 + I2 (new I1) */ \
\
	r8 = r3; \
	r3 -= r5;			/* R1 = R1 - R3 (final R3) */ \
	r5 += r8;			/* R3 = R1 + R3 (final R1) */ \
\
	/*IFNB <pre1>;*/ \
	xprefetchw(u8ptr(pre1)); \
	/*ENDIF;*/ \
\
	vec2f64 m8 = mem8; \
	r8 = m8;			/* I4 */ \
	r8 += r6;			/* I4 = I3 + I4 (new I3) */ \
	r6 -= m8;			/* I3 = I3 - I4 (new R4) */ \
\
	dest1 = r5; \
\
	r5 = r1;			/* Copy R2 */ \
	r1 -= r6;			/* R2 = R2 - R4 (final R4) */ \
	r6 += r5;			/* R4 = R2 + R4 (final R2) */ \
\
	r5 = r2;			/* Copy I2 */ \
	r2 -= r7;			/* I2 = I2 - I4 (final I4) */ \
	r7 += r5;			/* I4 = I2 + I4 (final I2) */ \
\
	/*IFNB <pre1>;*/ \
	xprefetchw(u8ptr(pre1+pre2)); \
	/*ENDIF;*/ \
\
	r5 = r4;			/* Copy I1 */ \
	r4 -= r8;			/* I1 = I1 - I3 (final I3) */ \
	r8 += r5;			/* I3 = I1 + I3 (final I1) */ \
}

/* 32-bit AMD K8 optimized versions of the above macros */ \

#if 0 // (@INSTR(,%xarch,<K8>) NE 0);
\

#define r4_x4c_simple_fft_mem(R1,R2,R3,R4,R5,R6,R7,R8,pre1,pre2) \
	xmm0 = R1;		/* R1 */ \
	xmm2 = R3;		/* R3 */ \
	xmm0 -= xmm2;		/* R1 = R1 - R3 (new R3) */ \
	xmm2 += R1;		/* R3 = R1 + R3 (new R1) */ \
\
	xmm4 = R5;		/* I1 */ \
	xmm6 = R7;		/* I3 */ \
	xmm4 -= xmm6;		/* I1 = I1 - I3 (new I3) */ \
	xmm6 += R5;		/* I3 = I1 + I3 (new I1) */ \
\
	xprefetchw(u8ptr(pre1)); \
\
	xmm1 = R2;		/* R2 */ \
	xmm3 = R4;		/* R4 */ \
	xmm1 -= xmm3;		/* R2 = R2 - R4 (new R4) */ \
	xmm3 += R2;		/* R4 = R2 + R4 (new R2) */ \
\
	xmm5 = R6;		/* I2 */ \
	xmm7 = R8;		/* I4 */ \
	xmm5 -= xmm7;		/* I2 = I2 - I4 (new I4) */ \
	xmm7 += R6;		/* I4 = I2 + I4 (new I2) */ \
\
	xprefetchw(u8ptr(pre1+pre2)); \
\
	xmm4 -= xmm1;		/* I3 = I3 - R4 (final I4) */ \
	xmm1 += xmm1; /* multwo */;			/* R4 = R4 * 2 */ \
\
	xmm0 -= xmm5;		/* R3 = R3 - I4 (final R3) */ \
	xmm5 += xmm5; /* multwo */;			/* I4 = I4 * 2 */ \
\
	xmm2 -= xmm3;		/* R1 = R1 - R2 (final R2) */ \
	xmm3 += xmm3; /* multwo */;			/* R2 = R2 * 2 */ \
\
	xmm6 -= xmm7;		/* I1 = I1 - I2 (final I2) */ \
	xmm7 += xmm7; /* multwo */;			/* I2 = I2 * 2 */ \
\
	xmm1 += xmm4;		/* R4 = I3 + R4 (final I3) */ \
	xmm5 += xmm0;		/* I4 = R3 + I4 (final R4) */ \
	xmm3 += xmm2;		/* R2 = R1 + R2 (final R1) */ \
	xmm7 += xmm6;		/* I2 = I1 + I2 (final I1) */ \
\
\

#define r4_x4c_simple_fft_partial_mem(r1, r2, r3, r4, r5, r6, r7, r8, mem3, mem4, mem7, mem8, pre1, pre2) \
	r3 = mem3;		/* R3 */ \
	r3 += r1;			/* R3 = R1 + R3 (new R1) */ \
	r1 -= mem3;		/* R1 = R1 - R3 (new R3) */ \
\
	r7 = mem7;		/* I3 */ \
	r7 += r5;			/* I3 = I1 + I3 (new I1) */ \
	r5 -= mem7;		/* I1 = I1 - I3 (new I3) */ \
\
	xprefetchw(u8ptr(pre1)); \
\
	r4 = mem4;		/* R4 */ \
	r4 += r2;			/* R4 = R2 + R4 (new R2) */ \
	r2 -= mem4;		/* R2 = R2 - R4 (new R4) */ \
\
	r8 = mem8;		/* I4 */ \
	r8 += r6;			/* I4 = I2 + I4 (new I2) */ \
	r6 -= mem8;		/* I2 = I2 - I4 (new I4) */ \
\
	xprefetchw(u8ptr(pre1+pre2)); \
\
	r5 -= r2;			/* I3 = I3 - R4 (final I4) */ \
	r2 += r2; /* multwo */			/* R4 = R4 * 2 */ \
\
	r1 -= r6;			/* R3 = R3 - I4 (final R3) */ \
	r6 += r6; /* multwo */			/* I4 = I4 * 2 */ \
\
	r3 -= r4;			/* R1 = R1 - R2 (final R2) */ \
	r4 += r4; /* multwo */			/* R2 = R2 * 2 */ \
\
	r7 -= r8;			/* I1 = I1 - I2 (final I2) */ \
	r8 += r8; /* multwo */			/* I2 = I2 * 2 */ \
\
	r2 += r5;			/* R4 = I3 + R4 (final I3) */ \
	r6 += r1;			/* I4 = R3 + I4 (final R4) */ \
	r4 += r3;			/* R2 = R1 + R2 (final R1) */ \
	r8 += r7;			/* I2 = I1 + I2 (final I1) */ \
\
\

#define r4_x4c_simple_unfft(r1, r2, r3, r4, r5, r6, r7, r8, mem8, pre1, pre2, dest1) \
	r1 -= r3;			/* R1 = R1 - R2 (new R2) */ \
	r3 += r3; /* multwo */			/* R2 = R2 * 2 */ \
	r3 += r1;			/* R2 = R1 + R2 (new R1) */ \
\
	r7 -= r5;			/* R4 = R4 - R3 (new I4) */ \
	r5 += r5; /* multwo */			/* R3 = R3 * 2 */ \
	r5 += r7;			/* R3 = R4 + R3 (new R3) */ \
\
	r2 -= r4;			/* I1 = I1 - I2 (new I2) */ \
	r4 += r4; /* multwo */			/* I2 = I2 * 2 */ \
	r4 += r2;			/* I2 = I1 + I2 (new I1) */ \
\
	r8 = mem8;		/* I4 */ \
	r8 += r6;			/* I4 = I3 + I4 (new I3) */ \
	r6 -= mem8;		/* I3 = I3 - I4 (new R4) */ \
\
	IFNB <pre1>; \
	xprefetchw(u8ptr(pre1)); \
	ENDIF; \
\
	r3 -= r5;			/* R1 = R1 - R3 (final R3) */ \
	r5 += r5; /* multwo */			/* R3 = R3 * 2 */ \
\
	r1 -= r6;			/* R2 = R2 - R4 (final R4) */ \
	r6 += r6; /* multwo */			/* R4 = R4 * 2 */ \
\
	r2 -= r7;			/* I2 = I2 - I4 (final I4) */ \
	r7 += r7; /* multwo */			/* R4 = R4 * 2 */ \
\
	r4 -= r8;			/* I1 = I1 - I3 (final I3) */ \
	r8 += r8; /* multwo */			/* I3 = I3 * 2 */ \
\
	IFNB <pre1>; \
	xprefetchw(u8ptr(pre1+pre2)); \
	ENDIF; \
\
	r5 += r3;			/* R3 = R1 + R3 (final R1) */ \
	r6 += r1;			/* R4 = R2 + R4 (final R2) */ \
	dest1 = r5; \
	r7 += r2;			/* I4 = I2 + I4 (final I2) */ \
	r8 += r4;			/* I3 = I1 + I3 (final I1) */ \
\

#endif
\
/* 64-bit Intel and AMD K10 optimized versions of the above macros */ \

#if 0 // (@INSTR(,%xarch,<CORE>) NE 0) OR (@INSTR(,%xarch,<P4>) NE 0) OR (@INSTR(,%xarch,<K10>) NE 0) OR (@INSTR(,%xarch,<BLEND>) NE 0);
#if 0 //DEF X86_64;
\

#define r4_x4c_simple_fft_mem(R1,R2,R3,R4,R5,R6,R7,R8,pre1,pre2) \
	xmm0 = R1;		/* R1 */ \
	xmm2 = R3;		/* R3 */ \
	xmm5 = xmm0;		/* Copy R1 */ \
	xmm0 -= xmm2;		/* R1 = R1 - R3 (new R3) */ \
	xmm2 += xmm5;		/* R3 = R1 + R3 (new R1) */ \
\
	xmm4 = R5;		/* I1 */ \
	xmm6 = R7;		/* I3 */ \
	xmm5 = xmm4;		/* Copy I1 */ \
	xmm4 -= xmm6;		/* I1 = I1 - I3 (new I3) */ \
	xmm6 += xmm5;		/* I3 = I1 + I3 (new I1) */ \
\
	xprefetchw(u8ptr(pre1)); \
\
	xmm1 = R2;		/* R2 */ \
	xmm3 = R4;		/* R4 */ \
	xmm5 = xmm1;		/* Copy R2 */ \
	xmm1 -= xmm3;		/* R2 = R2 - R4 (new R4) */ \
	xmm3 += xmm5;		/* R4 = R2 + R4 (new R2) */ \
\
	xmm5 = R6;		/* I2 */ \
	xmm7 = R8;		/* I4 */ \
	xmm8 = xmm5;		/* Copy I2 */ \
	xmm5 -= xmm7;		/* I2 = I2 - I4 (new I4) */ \
	xmm7 += xmm8;		/* I4 = I2 + I4 (new I2) */ \
\
	xprefetchw(u8ptr(pre1+pre2)); \
\
	xmm8 = xmm4;		/* Copy I3 */ \
	xmm4 -= xmm1;		/* I3 = I3 - R4 (final I4) */ \
	xmm1 += xmm8;		/* R4 = I3 + R4 (final I3) */ \
\
	xmm8 = xmm0;		/* Copy R3 */ \
	xmm0 -= xmm5;		/* R3 = R3 - I4 (final R3) */ \
	xmm5 += xmm8;		/* I4 = R3 + I4 (final R4) */ \
\
	xmm8 = xmm2;		/* Copy R1 */ \
	xmm2 -= xmm3;		/* R1 = R1 - R2 (final R2) */ \
	xmm3 += xmm8;		/* R2 = R1 + R2 (final R1) */ \
\
	xmm8 = xmm6;		/* Copy I1 */ \
	xmm6 -= xmm7;		/* I1 = I1 - I2 (final I2) */ \
	xmm7 += xmm8;		/* I2 = I1 + I2 (final I1) */ \
\
\

#define r4_x4c_simple_fft_partial_mem(r1, r2, r3, r4, r5, r6, r7, r8, mem3, mem4, mem7, mem8, pre1, pre2) \
	r3 = mem3;		/* R3 */ \
	r8 = r3;			/* Copy R3 */ \
	r3 += r1;			/* R3 = R1 + R3 (new R1) */ \
	r1 -= r8;			/* R1 = R1 - R3 (new R3) */ \
\
	r7 = mem7;		/* I3 */ \
	r8 = r7;			/* Copy I3 */ \
	r7 += r5;			/* I3 = I1 + I3 (new I1) */ \
	r5 -= r8;			/* I1 = I1 - I3 (new I3) */ \
\
	xprefetchw(u8ptr(pre1)); \
\
	r4 = mem4;		/* R4 */ \
	r8 = r4;			/* Copy R4 */ \
	r4 += r2;			/* R4 = R2 + R4 (new R2) */ \
	r2 -= r8;			/* R2 = R2 - R4 (new R4) */ \
\
	r8 = mem8;		/* I4 */ \
	xmm8 = r8;		/* Copy I4 */ \
	r8 += r6;			/* I4 = I2 + I4 (new I2) */ \
	r6 -= xmm8;		/* I2 = I2 - I4 (new I4) */ \
\
	xmm8 = r5;		/* Copy I3 */ \
	r5 -= r2;			/* I3 = I3 - R4 (final I4) */ \
	r2 += xmm8;		/* R4 = I3 + R4 (final I3) */ \
\
	xprefetchw(u8ptr(pre1+pre2)); \
\
	xmm8 = r1;		/* Copy R3 */ \
	r1 -= r6;			/* R3 = R3 - I4 (final R3) */ \
	r6 += xmm8;		/* I4 = R3 + I4 (final R4) */ \
\
	xmm8 = r3;		/* Copy R1 */ \
	r3 -= r4;			/* R1 = R1 - R2 (final R2) */ \
	r4 += xmm8;		/* R2 = R1 + R2 (final R1) */ \
\
	xmm8 = r7;		/* Copy I1 */ \
	r7 -= r8;			/* I1 = I1 - I2 (final I2) */ \
	r8 += xmm8;		/* I2 = I1 + I2 (final I1) */ \
\
\

#define r4_x4c_simple_unfft(r1, r2, r3, r4, r5, r6, r7, r8, mem8, pre1, pre2, dest1) \
	r8 = r1;			/* R1 */ \
	r1 -= r3;			/* R1 = R1 - R2 (new R2) */ \
	r3 += r8;			/* R2 = R1 + R2 (new R1) */ \
\
	r8 = r7;			/* R4 */ \
	r7 -= r5;			/* R4 = R4 - R3 (new I4) */ \
	r5 += r8;			/* R3 = R4 + R3 (new R3) */ \
\
	r8 = r2;			/* I2 */ \
	r2 -= r4;			/* I1 = I1 - I2 (new I2) */ \
	r4 += r8;			/* I2 = I1 + I2 (new I1) */ \
\
	r8 = r3; \
	r3 -= r5;			/* R1 = R1 - R3 (final R3) */ \
	r5 += r8;			/* R3 = R1 + R3 (final R1) */ \
\
	IFNB <pre1>; \
	xprefetchw(u8ptr(pre1)); \
	ENDIF; \
\
	r8 = mem8;		/* I4 */ \
	xmm8 = r8;		/* Copy I4 */ \
	r8 += r6;			/* I4 = I3 + I4 (new I3) */ \
	r6 -= xmm8;		/* I3 = I3 - I4 (new R4) */ \
\
	dest1 = r5; \
\
	r5 = r1;			/* Copy R2 */ \
	r1 -= r6;			/* R2 = R2 - R4 (final R4) */ \
	r6 += r5;			/* R4 = R2 + R4 (final R2) */ \
\
	r5 = r2;			/* Copy I2 */ \
	r2 -= r7;			/* I2 = I2 - I4 (final I4) */ \
	r7 += r5;			/* I4 = I2 + I4 (final I2) */ \
\
	IFNB <pre1>; \
	xprefetchw(u8ptr(pre1+pre2)); \
	ENDIF; \
\
	r5 = r4;			/* Copy I1 */ \
	r4 -= r8;			/* I1 = I1 - I3 (final I3) */ \
	r8 += r5;			/* I3 = I1 + I3 (final I1) */ \
\

#endif
#endif
\
/* 64-bit AMD K8 optimized versions of the above macros.  Derived from Intel 64-bit */ \
/* version -- could probably be optimized further. */ \

#if 0 // (@INSTR(,%xarch,<K8>) NE 0);
#if 0 //DEF X86_64;
\

#define r4_x4cl_four_complex_with_square_preload \
	xmm15 = XMM_TWO; \
\
\

#define r4_x4cl_four_complex_with_square(srcreg,srcinc,d1,d2) \
	xmm0 = xptr(srcreg+16);	/* I1 */ \
	xmm2 = xptr(srcreg+d2+16);	/* I3 */ \
	xmm0 -= xmm2;		/* I1 = I1 - I3 (new I3)		; 1-3 */ \
\
	xmm3 = xptr(srcreg+d1);	/* R2 */ \
	xmm5 = xptr(srcreg+d2+d1);	/* R4 */ \
	xmm3 -= xmm5;		/* R2 = R2 - R4 (new R4)		; 2-4 */ \
\
	xmm6 = xptr(srcreg);		/* R1 */ \
	xmm8 = xptr(srcreg+d2);	/* R3 */ \
	xmm6 -= xmm8;		/* R1 = R1 - R3 (new R3)		; 3-5 */ \
\
	xmm9 = xptr(srcreg+d1+16);	/* I2 */ \
	xmm11 = xptr(srcreg+d2+d1+16);/* I4 */ \
	xmm9 -= xmm11;		/* I2 = I2 - I4 (new I4)		; 4-6 */ \
\
	xmm0 -= xmm3;		/* I3 = I3 - R4 (final I4)		; 5-7 */ \
	xmm3 *= xmm15;		/* R4 = R4 * 2 */ \
	xmm3 += xmm0;		/* R4 = I3 + R4 (final I3)		; 6-8 */ \
\
	xmm6 -= xmm9;		/* R3 = R3 - I4 (final R3)		; 7-9 */ \
	xmm9 *= xmm15;		/* I4 = I4 * 2 */ \
	xmm9 += xmm6;		/* I4 = R3 + I4 (final R4)		; 8-10 */ \
\
	xmm13 = xmm0;		/* Copy I4 */ \
	xmm0 *= xmm0;		/* I4 * I4				; 8-12 */ \
\
	xmm2 += xptr(srcreg+16);	/* I3 = I1 + I3 (new I1)		; 9-11	avail 1,14,12 */ \
	xmm1 = xmm3;		/* Copy I3 */ \
	xmm3 *= xmm3;		/* I3 * I3				; 9-13	avail 14,12 */ \
\
	xmm11 += xptr(srcreg+d1+16);	/* I4 = I2 + I4 (new I2)		; 10-12	avail 14,12,10 */ \
	xmm1 *= xmm6;		/* I3 * R3				; 10-14 */ \
\
	xmm8 += xptr(srcreg);		/* R3 = R1 + R3 (new R1)		; 11-13	avail 14,12,10,7 */ \
	xmm13 *= xmm9;		/* I4 * R4				; 11-15 */ \
\
	xmm5 += xptr(srcreg+d1);	/* R4 = R2 + R4 (new R2)		; 12-14	avail 14,12,10,7,4 */ \
	xmm6 *= xmm6;		/* R3 * R3				; 12-16 */ \
\
	xprefetchw(u8ptr(srcreg+srcinc)); \
\
	xmm2 -= xmm11;		/* I1 = I1 - I2 (final I2)		; 13-15	avail 14,12,10,7 */ \
	xmm11 *= xmm15;		/* I2 = I2 * 2 */ \
	xmm9 *= xmm9;		/* R4 * R4				; 13-17 */ \
\
	xmm11 += xmm2;		/* I2 = I1 + I2 (final I1)		; 14-16	avail 14,12,10,7,4 */ \
	xmm10 = xptr(srcreg+48);	/*#2 I1 */ \
\
	xmm8 -= xmm5;		/* R1 = R1 - R2 (final R2)		; 15-17	avail 14,12,7 */ \
	xmm5 *= xmm15;		/* R2 = R2 * 2 */ \
	xmm1 *= xmm15;		/* I3 * R3 * 2 (new I3)			; 15-19 */ \
\
	xmm5 += xmm8;		/* R2 = R1 + R2 (final R1)		; 16-18	avail 14,12,7,4 */ \
	xmm4 = xmm2;		/* Copy I2 */ \
	xmm2 *= xmm2;		/* I2 * I2				; 16-20	avail 14,12,7 */ \
\
	xmm6 -= xmm3;		/* R3^2 - I3^2 (new R3)			; 17-19	avail 14,12,7,3 */ \
	xmm3 = xmm11;		/* Copy I1 */ \
	xmm11 *= xmm11;		/* I1 * I1				; 17-21	avail 14,12,7 */ \
\
	xmm9 -= xmm0;		/* R4^2 - I4^2 (new R4)			; 18-20	avail 14,12,7,0 */ \
	xmm4 *= xmm8;		/* I2 * R2				; 18-22 */ \
	xmm14 = xptr(srcreg+d2+48);	/*#2 I3 */ \
\
	xprefetchw(u8ptr(srcreg+srcinc+d1)); \
\
	xmm0 = xmm10;		/*#2 Copy I1 */ \
	xmm10 -= xmm14;		/*#2 I1 = I1 - I3 (new I3)		; 19-21	avail 12,7 */ \
	xmm3 *= xmm5;		/* I1 * R1				; 19-23 */ \
\
	xmm14 += xmm0;		/*#2 I3 = I1 + I3 (new I1)		; 20-22	avail 12,7,0 */ \
	xmm8 *= xmm8;		/* R2 * R2				; 20-24 */ \
	xmm12 = xptr(srcreg+d1+32);	/*#2 R2 */ \
\
	xmm9 -= xmm6;		/* R4 = R4 - R3 (new I4)		; 21-23	avail 7 */ \
	xmm6 *= xmm15;		/* R3 = R3 * 2 */ \
	xmm5 *= xmm5;		/* R1 * R1				; 21-25 */ \
\
	xmm6 += xmm9;		/* R3 = R4 + R3 (new R3)		; 22-24	avail 7,0 */ \
	xmm13 *= xmm15;		/* I4 * R4 * 2 (new I4)			; 22-26 (16) */ \
	xmm7 = xptr(srcreg+32);	/*#2 R1 */ \
\
	xmm12 -= xptr(srcreg+d2+d1+32);/*#2 R2 = R2 - R4 (new R4)		; 23-25 */ \
	xmm4 *= xmm15;		/* I2 * R2 * 2 (new I2)			; 23-27 */ \
	xmm0 = xptr(srcreg+d2+32);	/*#2 R3 */ \
\
	xmm7 -= xmm0;		/*#2 R1 = R1 - R3 (new R3)		; 24-26	avail none */ \
	xmm3 *= xmm15;		/* I1 * R1 * 2 (new I1)			; 24-28 */ \
\
	xmm8 -= xmm2;		/* R2^2 - I2^2 (new R2)			; 25-27	avail 2 */ \
\
	xprefetchw(u8ptr(srcreg+srcinc+d2)); \
\
	xmm5 -= xmm11;		/* R1^2 - I1^2 (new R1)			; 26-28	avail 2,11 */ \
\
	xmm13 += xmm1;		/* I4 = I3 + I4 (new I3)		; 27-29 */ \
	xmm1 *= xmm15;		/* I4 = I4 * 2 */ \
	xmm1 -= xmm13;		/* I3 = I3 - I4 (new R4)		; 28-30 */ \
	xmm11 = xptr(srcreg+d1+48);	/*#2 I2						avail 2 */ \
\
	xmm8 += xmm5;		/* R2 = R1 + R2 (new R1)		; 29-31 */ \
	xmm5 *= xmm15; \
	xmm5 -= xmm8;		/* R1 = R1 - R2 (new R2)		; 30-32 */ \
\
	xmm3 -= xmm4;		/* I1 = I1 - I2 (new I2)		; 31-33 */ \
	xmm4 *= xmm15;		/* I2 = I2 * 2 */ \
	xmm4 += xmm3;		/* I2 = I1 + I2 (new I1)		; 32-34 */ \
\
	xmm8 -= xmm6;		/* R1 = R1 - R3 (final R3)		; 33-35 */ \
	xmm6 *= xmm15;		/* R3 = R3 * 2 */ \
	xmm6 += xmm8;		/* R3 = R1 + R3 (final R1)		; 34-36	avail 2 storable 8,6 */ \
\
	xmm2 = xmm5;		/* Copy R2 */ \
	xmm5 -= xmm1;		/* R2 = R2 - R4 (final R4)		; 35-37	avail none storable 8,6,5 */ \
\
	xptr(srcreg+16) = xmm8;	/* Save R3				; 36 */ \
	xmm8 = xptr(srcreg+d2+d1+48);	/*#2 I4 */ \
	xptr(srcreg) = xmm6;		/* Save R1				; 37 */ \
	xmm6 = xmm11;		/*#2 Copy I2 */ \
	xmm11 -= xmm8;		/*#2 I2 = I2 - I4 (new I4)		; 36-38	avail none storable 5 */ \
\
	xptr(srcreg+d1+16) = xmm5;	/* Save R4				; 38 */ \
	xmm10 -= xmm12;		/*#2 I3 = I3 - R4 (final I4)		; 37-39	avail none */ \
	xmm12 *= xmm15;		/*#2 R4 = R4 * 2 */ \
\
	xmm12 += xmm10;		/*#2 R4 = I3 + R4 (final I3)		; 38-40	avail 5 */ \
\
	xprefetchw(u8ptr(srcreg+srcinc+d2+d1)); \
\
	xmm7 -= xmm11;		/*#2 R3 = R3 - I4 (final R3)		; 39-41	avail none */ \
	xmm11 *= xmm15;		/*#2 I4 = I4 * 2 */ \
\
	xmm11 += xmm7;		/*#2 I4 = R3 + I4 (final R4)		; 40-42 avail 5 */ \
	xmm5 = xmm10;		/*#2 Copy I4 */ \
	xmm10 *= xmm10;		/*#2 I4 * I4				; 40-44 avail none */ \
\
	xmm1 += xmm2;		/* R4 = R2 + R4 (final R2)		; 41-43	avail 2 storable 1 */ \
	xmm2 = xmm12;		/*#2 Copy I3 */ \
	xmm12 *= xmm12;		/*#2 I3 * I3				; 41-45	storable 1 */ \
\
	xmm6 += xmm8;		/*#2 I4 = I2 + I4 (new I2)		; 42-44	avail 8 storable 1 */ \
	xmm2 *= xmm7;		/*#2 I3 * R3				; 42-46 */ \
\
	xmm0 += xptr(srcreg+32);	/*#2 R3 = R1 + R3 (new R1)		; 43-45 */ \
	xmm5 *= xmm11;		/*#2 I4 * R4				; 43-47 */ \
\
	xmm8 = xptr(srcreg+d1+32);	/*#2 R2 */ \
	xmm8 += xptr(srcreg+d2+d1+32);	/*#2 R4 = R2 + R4 (new R2)		; 44-46 */ \
	xmm7 *= xmm7;		/*#2 R3 * R3				; 44-48 */ \
	xptr(srcreg+d1) = xmm1;	/* Save R2				; 44	avail 1 */ \
\
	xmm14 -= xmm6;		/*#2 I1 = I1 - I2 (final I2)		; 45-47	avail none */ \
	xmm6 *= xmm15;		/*#2 I2 = I2 * 2 */ \
	xmm11 *= xmm11;		/*#2 R4 * R4				; 45-49 */ \
\
	xmm6 += xmm14;		/*#2 I2 = I1 + I2 (final I1)		; 46-48	avail 1 */ \
\
	xmm0 -= xmm8;		/*#2 R1 = R1 - R2 (final R2)		; 47-49	avail none */ \
	xmm8 *= xmm15;		/*#2 R2 = R2 * 2 */ \
	xmm2 *= xmm15;		/*#2 I3 * R3 * 2 (new I3)		; 47-51 */ \
\
	xmm8 += xmm0;		/*#2 R2 = R1 + R2 (final R1)		; 48-50	avail 1 */ \
	xmm1 = xmm14;		/*#2 Copy I2 */ \
	xmm14 *= xmm14;		/*#2 I2 * I2				; 48-52	avail none */ \
\
	xmm7 -= xmm12;		/*#2 R3^2 - I3^2 (new R3)		; 49-51	avail 12 */ \
	xmm12 = xmm6;		/*#2 Copy I1 */ \
	xmm6 *= xmm6;		/*#2 I1 * I1				; 49-53	avail none */ \
\
	xmm11 -= xmm10;		/*#2 R4^2 - I4^2 (new R4)		; 50-52	avail 10 */ \
	xmm1 *= xmm0;		/*#2 I2 * R2				; 50-54 */ \
\
	xmm3 -= xmm9;		/* I2 = I2 - I4 (final I4)		; 51-53	avail none storable 3 */ \
	xmm9 *= xmm15;		/* I4 = I4 * 2 */ \
	xmm12 *= xmm8;		/*#2 I1 * R1				; 51-55 */ \
\
	xmm9 += xmm3;		/* I4 = I2 + I4 (final I2)		; 52-54	avail 10 storable 3,9 */ \
	xmm0 *= xmm0;		/*#2 R2 * R2				; 52-56 */ \
\
	xmm11 -= xmm7;		/*#2 R4 = R4 - R3 (new I4)		; 53-55	avail none storable 3,9 */ \
	xmm7 *= xmm15;		/*#2 R3 = R3 * 2 */ \
	xmm8 *= xmm8;		/*#2 R1 * R1				; 53-57 */ \
\
	xmm7 += xmm11;		/*#2 R3 = R4 + R3 (new R3)		; 54-56	avail 10 storable 3,9 */ \
	xmm5 *= xmm15;		/*#2 I4 * R4 * 2 (new I4)		; 54-58 */ \
	xptr(srcreg+d1+48) = xmm3;	/* Save I4				; 54	avail 10,3 storable 9 */ \
\
	xmm4 -= xmm13;		/* I1 = I1 - I3 (final I3)		; 55-57	avail 10 storable 9,4 */ \
	xmm13 *= xmm15;		/* I3 = I3 * 2 */ \
	xmm1 *= xmm15;		/*#2 I2 * R2 * 2 (new I2)		; 55-59 */ \
	xptr(srcreg+d1+32) = xmm9;	/* Save I2				; 55	avail 10,9 storable 4 */ \
\
	xmm13 += xmm4;		/* I3 = I1 + I3 (final I1)		; 56-58	avail 10,9,3 storable 4,13 */ \
	xmm12 *= xmm15;		/*#2 I1 * R1 * 2 (new I1)		; 56-60 */ \
\
	xmm0 -= xmm14;		/*#2 R2^2 - I2^2 (new R2)		; 57-59	avail 10,9,3,14 storable 4,13 */ \
\
	xmm8 -= xmm6;		/*#2 R1^2 - I1^2 (new R1)		; 58-60	avail 10,9,3,14,6 storable 4,13 */ \
	xptr(srcreg+48) = xmm4;	/* Save I3				; 58	avail 10,9,3,14,6,4 storable 13 */ \
\
	xmm5 += xmm2;		/*#2 I4 = I3 + I4 (new I3)		; 59-61 */ \
	xmm2 *= xmm15;		/*#2 I4 = I4 * 2 */ \
	xmm2 -= xmm5;		/*#2 I3 = I3 - I4 (new R4)		; 60-62 */ \
	xptr(srcreg+32) = xmm13;	/* Save I1				; 59	avail 10,9,3,14,6,4,13 */ \
\
	xmm0 += xmm8;		/*#2 R2 = R1 + R2 (new R1)		; 61-63 */ \
	xmm8 *= xmm15;		/*#2 R2 = R2 * 2 */ \
	xmm8 -= xmm0;		/*#2 R1 = R1 - R2 (new R2)		; 62-64 */ \
\
	xmm12 -= xmm1;		/*#2 I1 = I1 - I2 (new I2)		; 63-65 */ \
	xmm1 *= xmm15;		/*#2 I2 = I2 * 2 */ \
	xmm1 += xmm12;		/*#2 I2 = I1 + I2 (new I1)		; 64-66 */ \
\
	xmm0 -= xmm7;		/*#2 R1 = R1 - R3 (final R3)		; 65-67 */ \
	xmm7 *= xmm15;		/*#2 R3 = R3 * 2 */ \
	xmm7 += xmm0;		/*#2 R3 = R1 + R3 (final R1)		; 66-68	avail 10,9,3,14,6,4,13 storable 0,7 */ \
\
	xmm8 -= xmm2;		/*#2 R2 = R2 - R4 (final R4)		; 67-69 */ \
	xmm2 *= xmm15;		/*#2 R4 = R4 * 2 */ \
	xmm2 += xmm8;		/*#2 R4 = R2 + R4 (final R2)		; 68-70	avail 10,9,3,14,6,4,13 storable 0,7,8,2 */ \
	xptr(srcreg+d2+16) = xmm0;	/*#2 Save R3				; 68 */ \
\
	xmm12 -= xmm11;		/*#2 I2 = I2 - I4 (final I4)		; 69-71	avail 10,9,14,6,4,13,0 storable 7,8,2,12 */ \
	xmm11 *= xmm15;		/*#2 I4 = I4 * 2 */ \
	xptr(srcreg+d2) = xmm7;	/*#2 Save R1				; 69 */ \
	xmm11 += xmm12;		/*#2 I4 = I2 + I4 (final I2)		; 70-72	avail 10,9,3,14,6,4,13,0,7 storable 8,2,12,11 */ \
	xptr(srcreg+d2+d1+16) = xmm8;	/*#2 Save R4				; 70 */ \
\
	xmm1 -= xmm5;		/*#2 I1 = I1 - I3 (final I3)		; 71-73 */ \
	xmm5 *= xmm15;		/*#2 I3 = I3 * 2 */ \
	xptr(srcreg+d2+d1) = xmm2;	/*#2 Save R2				; 71 */ \
	xmm5 += xmm1;		/*#2 I3 = I1 + I3 (final I1)		; 72-74	avail 10,9,3,14,6,4,13,0,7,8,2 storable,12,11,1,5 */ \
	xptr(srcreg+d2+d1+48) = xmm12; /*#2 Save I4				; 72 */ \
\
	xptr(srcreg+d2+d1+32) = xmm11; /*#2 Save I2				; 73 */ \
	xptr(srcreg+d2+48) = xmm1;	/*#2 Save I3				; 74 */ \
	xptr(srcreg+d2+32) = xmm5;	/*#2 Save I1				; 75 */ \
\
	srcreg += srcinc; \
\

#endif
#endif
\
/* 64-bit Intel and AMD K10 optimized versions of the above macros. */ \

#if 0 // (@INSTR(,%xarch,<CORE>) NE 0) OR (@INSTR(,%xarch,<P4>) NE 0) OR (@INSTR(,%xarch,<K10>) NE 0) OR (@INSTR(,%xarch,<BLEND>) NE 0);
#if 0 //DEF X86_64;
\

#define r4_x4cl_four_complex_fft_final(srcreg,srcinc,d1,d2) \
	xmm15 = xptr(srcreg);		/* R1 */ \
	xmm2 = xptr(srcreg+d2);	/* R3 */ \
	xmm0 = xmm15;		/* Copy R1 */ \
	xmm15 -= xmm2;		/* R1 = R1 - R3 (new R3)		; 1-3 */ \
	xmm2 += xmm0;		/* R3 = R1 + R3 (new R1)		; 2-4 */ \
\
	xmm4 = xptr(srcreg+16);	/* I1 */ \
	xmm6 = xptr(srcreg+d2+16);	/* I3 */ \
	xmm0 = xmm4;		/* Copy I1 */ \
	xmm4 -= xmm6;		/* I1 = I1 - I3 (new I3)		; 3-5 */ \
	xmm6 += xmm0;		/* I3 = I1 + I3 (new I1)		; 4-6 */ \
\
	xprefetchw(u8ptr(srcreg+srcinc)); \
\
	xmm1 = xptr(srcreg+d1);	/* R2 */ \
	xmm3 = xptr(srcreg+d2+d1);	/* R4 */ \
	xmm0 = xmm1;		/* Copy R2 */ \
	xmm1 -= xmm3;		/* R2 = R2 - R4 (new R4)		; 5-7 */ \
	xmm3 += xmm0;		/* R4 = R2 + R4 (new R2)		; 6-8 */ \
\
	xmm5 = xptr(srcreg+d1+16);	/* I2 */ \
	xmm7 = xptr(srcreg+d2+d1+16);	/* I4 */ \
	xmm0 = xmm5;		/* Copy I2 */ \
	xmm5 -= xmm7;		/* I2 = I2 - I4 (new I4)		; 7-9 */ \
	xmm7 += xmm0;		/* I4 = I2 + I4 (new I2)		; 8-10 */ \
\
	xprefetchw(u8ptr(srcreg+srcinc+d1)); \
\
	xmm0 = xmm4;		/* Copy I3 */ \
	xmm4 -= xmm1;		/* I3 = I3 - R4 (final I4)		; 9-11 */ \
	xmm1 += xmm0;		/* R4 = I3 + R4 (final I3)		; 10-12 */ \
\
	xmm12 = xptr(srcreg+d2+d1+48);/* I4 */ \
	xmm13 = xptr(srcreg+d1+48);	/* R6 */ \
	xmm0 = xmm12;		/* Copy I4 */ \
	xmm12 += xmm13;		/* I4 = I2 + I4 (new I2)		; 11-13 */ \
	xmm13 -= xmm0;		/* I2 = I2 - I4 (new I4)		; 12-14 */ \
	xptr(srcreg+d1+48) = xmm4;	/* Save I4				; 12 */ \
\
	xmm0 = xmm2;		/* Copy R1 */ \
	xmm2 -= xmm3;		/* R1 = R1 - R2 (final R2)		; 13-15 */ \
	xptr(srcreg+d1+16) = xmm1;	/* Save I3				; 13 */ \
	xmm3 += xmm0;		/* R2 = R1 + R2 (final R1)		; 14-16 */ \
\
	xmm4 = xptr(srcreg+d2+32);	/* R3 */ \
	xmm1 = xptr(srcreg+32);	/* R1 */ \
	xmm0 = xmm4;		/* Copy R3 */ \
	xmm4 += xmm1;		/* R3 = R1 + R3 (new R1)		; 15-17 */ \
	xmm1 -= xmm0;		/* R1 = R1 - R3 (new R3)		; 16-18 */ \
	xptr(srcreg+32) = xmm2;	/* Save R2				; 16 */ \
\
	xprefetchw(u8ptr(srcreg+srcinc+d2)); \
\
	xmm0 = xmm6;		/* Copy I1 */ \
	xmm6 -= xmm7;		/* I1 = I1 - I2 (final I2)		; 17-19 */ \
	xptr(srcreg) = xmm3;		/* Save R1				; 17 */ \
	xmm7 += xmm0;		/* I2 = I1 + I2 (final I1)		; 18-20 */ \
\
	xmm2 = xptr(srcreg+d2+48);	/* I3 */ \
	xmm3 = xptr(srcreg+48);	/* R5 */ \
	xmm0 = xmm2;		/* Copy I3 */ \
	xmm2 += xmm3;		/* I3 = I1 + I3 (new I1)		; 19-21 */ \
	xmm3 -= xmm0;		/* I1 = I1 - I3 (new I3)		; 20-22 */ \
	xptr(srcreg+48) = xmm6;	/* Save I2				; 20 */ \
\
	xmm0 = xmm15;		/* Copy R3 */ \
	xmm15 -= xmm5;		/* R3 = R3 - I4 (final R3)		; 21-23 */ \
	xptr(srcreg+16) = xmm7;	/* Save I1				; 21 */ \
	xmm5 += xmm0;		/* I4 = R3 + I4 (final R4)		; 22-24 */ \
\
	xmm6 = xptr(srcreg+d2+d1+32);	/* R4 */ \
	xmm7 = xptr(srcreg+d1+32);	/* R2 */ \
	xmm0 = xmm6;		/* Copy R4 */ \
	xmm6 += xmm7;		/* R4 = R2 + R4 (new R2)		; 23-25 */ \
	xmm7 -= xmm0;		/* R2 = R2 - R4 (new R4)		; 24-26 */ \
	xptr(srcreg+d1) = xmm15;	/* Save R3				; 24 */ \
\
	xmm0 = xmm3;		/* Copy I3 */ \
	xmm3 -= xmm7;		/* I3 = I3 - R4 (final I4)		; 25-27 */ \
	xmm7 += xmm0;		/* R4 = I3 + R4 (final I3)		; 26-28 */ \
	xptr(srcreg+d1+32) = xmm5;	/* Save R4				; 25 */ \
\
	xprefetchw(u8ptr(srcreg+srcinc+d2+d1)); \
\
	xmm0 = xmm1;		/* Copy R3 */ \
	xmm1 -= xmm13;		/* R3 = R3 - I4 (final R3)		; 27-29 */ \
	xmm13 += xmm0;		/* I4 = R3 + I4 (final R4)		; 28-30 */ \
	xptr(srcreg+d2+d1+48) = xmm3;	/* Save I4				; 28 */ \
\
	xmm0 = xmm4;		/* Copy R1 */ \
	xmm4 -= xmm6;		/* R1 = R1 - R2 (final R2)		; 29-31 */ \
	xptr(srcreg+d2+d1+16) = xmm7;	/* Save I3				; 29 */ \
	xmm6 += xmm0;		/* R2 = R1 + R2 (final R1)		; 30-32 */ \
	xptr(srcreg+d2+d1) = xmm1;	/* Save R3				; 30 */ \
\
	xmm0 = xmm2;		/* Copy I1 */ \
	xmm2 -= xmm12;		/* I1 = I1 - I2 (final I2)		; 31-33 */ \
	xptr(srcreg+d2+d1+32) = xmm13; /* Save R4				; 31 */ \
	xmm12 += xmm0;		/* I2 = I1 + I2 (final I1)		; 32-34 */ \
\
	xptr(srcreg+d2+32) = xmm4;	/* Save R2				; 32 */ \
	xptr(srcreg+d2) = xmm6;	/* Save R1				; 33 */ \
	xptr(srcreg+d2+48) = xmm2;	/* Save I2				; 34 */ \
	xptr(srcreg+d2+16) = xmm12;	/* Save I1				; 35 */ \
	srcreg += srcinc; \
\
\

#define r4_x4cl_four_complex_with_square_preload \
	xmm15 = XMM_TWO; \
\
\

#define r4_x4cl_four_complex_with_square(srcreg,srcinc,d1,d2) \
	xmm0 = xptr(srcreg+16);	/* I1 */ \
	xmm1 = xptr(srcreg+d2+16);	/* I3 */ \
	xmm2 = xmm0;		/* Copy I1 */ \
	xmm0 -= xmm1;		/* I1 = I1 - I3 (new I3)		; 1-3 */ \
\
	xmm3 = xptr(srcreg+d1);	/* R2 */ \
	xmm4 = xptr(srcreg+d2+d1);	/* R4 */ \
	xmm5 = xmm3;		/* Copy R2 */ \
	xmm3 -= xmm4;		/* R2 = R2 - R4 (new R4)		; 2-4 */ \
\
	xmm6 = xptr(srcreg);		/* R1 */ \
	xmm7 = xptr(srcreg+d2);	/* R3 */ \
	xmm8 = xmm6;		/* Copy R1 */ \
	xmm6 -= xmm7;		/* R1 = R1 - R3 (new R3)		; 3-5 */ \
\
	xmm9 = xptr(srcreg+d1+16);	/* I2 */ \
	xmm10 = xptr(srcreg+d2+d1+16);/* I4 */ \
	xmm11 = xmm9;		/* Copy I2 */ \
	xmm9 -= xmm10;		/* I2 = I2 - I4 (new I4)		; 4-6 */ \
\
	xmm12 = xmm0;		/* Copy I3 */ \
	xmm0 -= xmm3;		/* I3 = I3 - R4 (final I4)		; 5-7 */ \
\
	xmm3 += xmm12;		/* R4 = I3 + R4 (final I3)		; 6-8 */ \
\
	xmm12 = xmm6;		/* Copy R3 */ \
	xmm6 -= xmm9;		/* R3 = R3 - I4 (final R3)		; 7-9 */ \
\
	xmm9 += xmm12;		/* I4 = R3 + I4 (final R4)		; 8-10 */ \
	xmm13 = xmm0;		/* Copy I4 */ \
	xmm0 *= xmm0;		/* I4 * I4				; 8-12 */ \
\
	xmm2 += xmm1;		/* I3 = I1 + I3 (new I1)		; 9-11	avail 1,14,12 */ \
	xmm1 = xmm3;		/* Copy I3 */ \
	xmm3 *= xmm3;		/* I3 * I3				; 9-13	avail 14,12 */ \
\
	xmm11 += xmm10;		/* I4 = I2 + I4 (new I2)		; 10-12	avail 14,12,10 */ \
	xmm1 *= xmm6;		/* I3 * R3				; 10-14 */ \
\
	xmm8 += xmm7;		/* R3 = R1 + R3 (new R1)		; 11-13	avail 14,12,10,7 */ \
	xmm13 *= xmm9;		/* I4 * R4				; 11-15 */ \
\
	xmm5 += xmm4;		/* R4 = R2 + R4 (new R2)		; 12-14	avail 14,12,10,7,4 */ \
	xmm6 *= xmm6;		/* R3 * R3				; 12-16 */ \
\
	xprefetchw(u8ptr(srcreg+srcinc)); \
\
	xmm4 = xmm2;		/* Copy I1 */ \
	xmm2 -= xmm11;		/* I1 = I1 - I2 (final I2)		; 13-15	avail 14,12,10,7 */ \
	xmm9 *= xmm9;		/* R4 * R4				; 13-17 */ \
\
	xmm11 += xmm4;		/* I2 = I1 + I2 (final I1)		; 14-16	avail 14,12,10,7,4 */ \
	xmm10 = xptr(srcreg+48);	/*#2 I1 */ \
\
	xmm4 = xmm8;		/* Copy R1 */ \
	xmm8 -= xmm5;		/* R1 = R1 - R2 (final R2)		; 15-17	avail 14,12,7 */ \
	xmm1 *= xmm15;		/* I3 * R3 * 2 (new I3)			; 15-19 */ \
\
	xmm5 += xmm4;		/* R2 = R1 + R2 (final R1)		; 16-18	avail 14,12,7,4 */ \
	xmm4 = xmm2;		/* Copy I2 */ \
	xmm2 *= xmm2;		/* I2 * I2				; 16-20	avail 14,12,7 */ \
\
	xmm6 -= xmm3;		/* R3^2 - I3^2 (new R3)			; 17-19	avail 14,12,7,3 */ \
	xmm3 = xmm11;		/* Copy I1 */ \
	xmm11 *= xmm11;		/* I1 * I1				; 17-21	avail 14,12,7 */ \
\
	xmm9 -= xmm0;		/* R4^2 - I4^2 (new R4)			; 18-20	avail 14,12,7,0 */ \
	xmm4 *= xmm8;		/* I2 * R2				; 18-22 */ \
	xmm14 = xptr(srcreg+d2+48);	/*#2 I3 */ \
\
	xprefetchw(u8ptr(srcreg+srcinc+d1)); \
\
	xmm0 = xmm10;		/*#2 Copy I1 */ \
	xmm10 -= xmm14;		/*#2 I1 = I1 - I3 (new I3)		; 19-21	avail 12,7 */ \
	xmm3 *= xmm5;		/* I1 * R1				; 19-23 */ \
\
	xmm14 += xmm0;		/*#2 I3 = I1 + I3 (new I1)		; 20-22	avail 12,7,0 */ \
	xmm8 *= xmm8;		/* R2 * R2				; 20-24 */ \
	xmm12 = xptr(srcreg+d1+32);	/*#2 R2 */ \
\
	xmm0 = xmm9;		/* Copy R4 */ \
	xmm9 -= xmm6;		/* R4 = R4 - R3 (new I4)		; 21-23	avail 7 */ \
	xmm5 *= xmm5;		/* R1 * R1				; 21-25 */ \
\
	xmm6 += xmm0;		/* R3 = R4 + R3 (new R3)		; 22-24	avail 7,0 */ \
	xmm13 *= xmm15;		/* I4 * R4 * 2 (new I4)			; 22-26 (16) */ \
	xmm7 = xptr(srcreg+32);	/*#2 R1 */ \
\
	xmm12 -= xptr(srcreg+d2+d1+32);/*#2 R2 = R2 - R4 (new R4)		; 23-25 */ \
	xmm4 *= xmm15;		/* I2 * R2 * 2 (new I2)			; 23-27 */ \
	xmm0 = xptr(srcreg+d2+32);	/*#2 R3 */ \
\
	xmm7 -= xmm0;		/*#2 R1 = R1 - R3 (new R3)		; 24-26	avail none */ \
	xmm3 *= xmm15;		/* I1 * R1 * 2 (new I1)			; 24-28 */ \
\
	xmm8 -= xmm2;		/* R2^2 - I2^2 (new R2)			; 25-27	avail 2 */ \
\
	xprefetchw(u8ptr(srcreg+srcinc+d2)); \
\
	xmm5 -= xmm11;		/* R1^2 - I1^2 (new R1)			; 26-28	avail 2,11 */ \
\
	xmm2 = xmm13;		/* Copy I4 */ \
	xmm13 += xmm1;		/* I4 = I3 + I4 (new I3)		; 27-29 */ \
	xmm1 -= xmm2;		/* I3 = I3 - I4 (new R4)		; 28-30 */ \
	xmm11 = xptr(srcreg+d1+48);	/*#2 I2						avail 2 */ \
\
	xmm2 = xmm8;		/* Copy R2 */ \
	xmm8 += xmm5;		/* R2 = R1 + R2 (new R1)		; 29-31 */ \
	xmm5 -= xmm2;		/* R1 = R1 - R2 (new R2)		; 30-32 */ \
\
	xmm2 = xmm3;		/* Copy I1 */ \
	xmm3 -= xmm4;		/* I1 = I1 - I2 (new I2)		; 31-33 */ \
	xmm4 += xmm2;		/* I2 = I1 + I2 (new I1)		; 32-34 */ \
\
	xmm2 = xmm8;		/* Copy R1 */ \
	xmm8 -= xmm6;		/* R1 = R1 - R3 (final R3)		; 33-35 */ \
	xmm6 += xmm2;		/* R3 = R1 + R3 (final R1)		; 34-36	avail 2 storable 8,6 */ \
\
	xmm2 = xmm5;		/* Copy R2 */ \
	xmm5 -= xmm1;		/* R2 = R2 - R4 (final R4)		; 35-37	avail none storable 8,6,5 */ \
\
	xptr(srcreg+16) = xmm8;	/* Save R3				; 36 */ \
	xmm8 = xptr(srcreg+d2+d1+48);	/*#2 I4 */ \
	xptr(srcreg) = xmm6;		/* Save R1				; 37 */ \
	xmm6 = xmm11;		/*#2 Copy I2 */ \
	xmm11 -= xmm8;		/*#2 I2 = I2 - I4 (new I4)		; 36-38	avail none storable 5 */ \
\
	xptr(srcreg+d1+16) = xmm5;	/* Save R4				; 38 */ \
	xmm5 = xmm10;		/*#2 Copy I3 */ \
	xmm10 -= xmm12;		/*#2 I3 = I3 - R4 (final I4)		; 37-39	avail none */ \
\
	xmm12 += xmm5;		/*#2 R4 = I3 + R4 (final I3)		; 38-40	avail 5 */ \
\
	xprefetchw(u8ptr(srcreg+srcinc+d2+d1)); \
\
	xmm5 = xmm7;		/*#2 Copy R3 */ \
	xmm7 -= xmm11;		/*#2 R3 = R3 - I4 (final R3)		; 39-41	avail none */ \
\
	xmm11 += xmm5;		/*#2 I4 = R3 + I4 (final R4)		; 40-42 avail 5 */ \
	xmm5 = xmm10;		/*#2 Copy I4 */ \
	xmm10 *= xmm10;		/*#2 I4 * I4				; 40-44 avail none */ \
\
	xmm1 += xmm2;		/* R4 = R2 + R4 (final R2)		; 41-43	avail 2 storable 1 */ \
	xmm2 = xmm12;		/*#2 Copy I3 */ \
	xmm12 *= xmm12;		/*#2 I3 * I3				; 41-45	storable 1 */ \
\
	xmm6 += xmm8;		/*#2 I4 = I2 + I4 (new I2)		; 42-44	avail 8 storable 1 */ \
	xmm2 *= xmm7;		/*#2 I3 * R3				; 42-46 */ \
\
	xmm0 += xptr(srcreg+32);	/*#2 R3 = R1 + R3 (new R1)		; 43-45 */ \
	xmm5 *= xmm11;		/*#2 I4 * R4				; 43-47 */ \
\
	xmm8 = xptr(srcreg+d1+32);	/*#2 R2 */ \
	xmm8 += xptr(srcreg+d2+d1+32);	/*#2 R4 = R2 + R4 (new R2)		; 44-46 */ \
	xmm7 *= xmm7;		/*#2 R3 * R3				; 44-48 */ \
	xptr(srcreg+d1) = xmm1;	/* Save R2				; 44	avail 1 */ \
\
	xmm1 = xmm14;		/*#2 Copy I1 */ \
	xmm14 -= xmm6;		/*#2 I1 = I1 - I2 (final I2)		; 45-47	avail none */ \
	xmm11 *= xmm11;		/*#2 R4 * R4				; 45-49 */ \
\
	xmm6 += xmm1;		/*#2 I2 = I1 + I2 (final I1)		; 46-48	avail 1 */ \
\
	xmm1 = xmm0;		/*#2 Copy R1 */ \
	xmm0 -= xmm8;		/*#2 R1 = R1 - R2 (final R2)		; 47-49	avail none */ \
	xmm2 *= xmm15;		/*#2 I3 * R3 * 2 (new I3)		; 47-51 */ \
\
	xmm8 += xmm1;		/*#2 R2 = R1 + R2 (final R1)		; 48-50	avail 1 */ \
	xmm1 = xmm14;		/*#2 Copy I2 */ \
	xmm14 *= xmm14;		/*#2 I2 * I2				; 48-52	avail none */ \
\
	xmm7 -= xmm12;		/*#2 R3^2 - I3^2 (new R3)		; 49-51	avail 12 */ \
	xmm12 = xmm6;		/*#2 Copy I1 */ \
	xmm6 *= xmm6;		/*#2 I1 * I1				; 49-53	avail none */ \
\
	xmm11 -= xmm10;		/*#2 R4^2 - I4^2 (new R4)		; 50-52	avail 10 */ \
	xmm1 *= xmm0;		/*#2 I2 * R2				; 50-54 */ \
\
	xmm10 = xmm3;		/* Copy I2 */ \
	xmm3 -= xmm9;		/* I2 = I2 - I4 (final I4)		; 51-53	avail none storable 3 */ \
	xmm12 *= xmm8;		/*#2 I1 * R1				; 51-55 */ \
\
	xmm9 += xmm10;		/* I4 = I2 + I4 (final I2)		; 52-54	avail 10 storable 3,9 */ \
	xmm0 *= xmm0;		/*#2 R2 * R2				; 52-56 */ \
\
	xmm10 = xmm11;		/*#2 Copy R4 */ \
	xmm11 -= xmm7;		/*#2 R4 = R4 - R3 (new I4)		; 53-55	avail none storable 3,9 */ \
	xmm8 *= xmm8;		/*#2 R1 * R1				; 53-57 */ \
\
	xmm7 += xmm10;		/*#2 R3 = R4 + R3 (new R3)		; 54-56	avail 10 storable 3,9 */ \
	xmm5 *= xmm15;		/*#2 I4 * R4 * 2 (new I4)		; 54-58 */ \
	xptr(srcreg+d1+48) = xmm3;	/* Save I4				; 54	avail 10,3 storable 9 */ \
\
	xmm3 = xmm4;		/* Copy I1 */ \
	xmm4 -= xmm13;		/* I1 = I1 - I3 (final I3)		; 55-57	avail 10 storable 9,4 */ \
	xmm1 *= xmm15;		/*#2 I2 * R2 * 2 (new I2)		; 55-59 */ \
	xptr(srcreg+d1+32) = xmm9;	/* Save I2				; 55	avail 10,9 storable 4 */ \
\
	xmm13 += xmm3;		/* I3 = I1 + I3 (final I1)		; 56-58	avail 10,9,3 storable 4,13 */ \
	xmm12 *= xmm15;		/*#2 I1 * R1 * 2 (new I1)		; 56-60 */ \
\
	xmm0 -= xmm14;		/*#2 R2^2 - I2^2 (new R2)		; 57-59	avail 10,9,3,14 storable 4,13 */ \
\
	xmm8 -= xmm6;		/*#2 R1^2 - I1^2 (new R1)		; 58-60	avail 10,9,3,14,6 storable 4,13 */ \
	xptr(srcreg+48) = xmm4;	/* Save I3				; 58	avail 10,9,3,14,6,4 storable 13 */ \
\
	xmm3 = xmm5;		/*#2 Copy I4 */ \
	xmm5 += xmm2;		/*#2 I4 = I3 + I4 (new I3)		; 59-61 */ \
	xmm2 -= xmm3;		/*#2 I3 = I3 - I4 (new R4)		; 60-62 */ \
	xptr(srcreg+32) = xmm13;	/* Save I1				; 59	avail 10,9,3,14,6,4,13 */ \
\
	xmm3 = xmm0;		/*#2 Copy R2 */ \
	xmm0 += xmm8;		/*#2 R2 = R1 + R2 (new R1)		; 61-63 */ \
	xmm8 -= xmm3;		/*#2 R1 = R1 - R2 (new R2)		; 62-64 */ \
\
	xmm3 = xmm12;		/*#2 Copy I1 */ \
	xmm12 -= xmm1;		/*#2 I1 = I1 - I2 (new I2)		; 63-65 */ \
	xmm1 += xmm3;		/*#2 I2 = I1 + I2 (new I1)		; 64-66 */ \
\
	xmm3 = xmm0;		/*#2 Copy R1 */ \
	xmm0 -= xmm7;		/*#2 R1 = R1 - R3 (final R3)		; 65-67 */ \
	xmm7 += xmm3;		/*#2 R3 = R1 + R3 (final R1)		; 66-68	avail 10,9,3,14,6,4,13 storable 0,7 */ \
\
	xmm3 = xmm8;		/*#2 Copy R2 */ \
	xmm8 -= xmm2;		/*#2 R2 = R2 - R4 (final R4)		; 67-69 */ \
	xmm2 += xmm3;		/*#2 R4 = R2 + R4 (final R2)		; 68-70	avail 10,9,3,14,6,4,13 storable 0,7,8,2 */ \
	xptr(srcreg+d2+16) = xmm0;	/*#2 Save R3				; 68 */ \
\
	xmm3 = xmm12;		/*#2 Copy I2 */ \
	xmm12 -= xmm11;		/*#2 I2 = I2 - I4 (final I4)		; 69-71	avail 10,9,14,6,4,13,0 storable 7,8,2,12 */ \
	xptr(srcreg+d2) = xmm7;	/*#2 Save R1				; 69 */ \
	xmm11 += xmm3;		/*#2 I4 = I2 + I4 (final I2)		; 70-72	avail 10,9,3,14,6,4,13,0,7 storable 8,2,12,11 */ \
	xptr(srcreg+d2+d1+16) = xmm8;	/*#2 Save R4				; 70 */ \
\
	xmm3 = xmm1;		/*#2 Copy I1 */ \
	xmm1 -= xmm5;		/*#2 I1 = I1 - I3 (final I3)		; 71-73 */ \
	xptr(srcreg+d2+d1) = xmm2;	/*#2 Save R2				; 71 */ \
	xmm5 += xmm3;		/*#2 I3 = I1 + I3 (final I1)		; 72-74	avail 10,9,3,14,6,4,13,0,7,8,2 storable,12,11,1,5 */ \
	xptr(srcreg+d2+d1+48) = xmm12; /*#2 Save I4				; 72 */ \
\
	xptr(srcreg+d2+d1+32) = xmm11; /*#2 Save I2				; 73 */ \
	xptr(srcreg+d2+48) = xmm1;	/*#2 Save I3				; 74 */ \
	xptr(srcreg+d2+32) = xmm5;	/*#2 Save I1				; 75 */ \
\
	srcreg += srcinc; \
\
\

#define r4_x4cl_four_complex_with_mult(srcreg,srcinc,d1,d2) \
	xmm1 = xptr(srcreg+16);		/* I1 */ \
	xmm2 = xptr(srcreg+d2+16);		/* I3 */ \
	xmm3 = xmm1;			/* Copy I1 */ \
	xmm1 -= xmm2;			/* I1 = I1 - I3 (new I3)	; 1-3 */ \
\
	xmm4 = xptr(srcreg+d1);		/* R2 */ \
	xmm5 = xptr(srcreg+d2+d1);		/* R4 */ \
	xmm6 = xmm4;			/* Copy R2 */ \
	xmm4 -= xmm5;			/* R2 = R2 - R4 (new R4)	; 2-4 */ \
\
	xmm7 = xptr(srcreg);			/* R1 */ \
	xmm8 = xptr(srcreg+d2);		/* R3 */ \
	xmm9 = xmm7;			/* Copy R1 */ \
	xmm7 -= xmm8;			/* R1 = R1 - R3 (new R3)	; 3-5 */ \
\
	xmm10 = xptr(srcreg+d1+16);		/* I2 */ \
	xmm11 = xptr(srcreg+d2+d1+16);	/* I4 */ \
	xmm12 = xmm10;			/* Copy I2 */ \
	xmm10 -= xmm11;			/* I2 = I2 - I4 (new I4)	; 4-6	avail 0,13+ */ \
\
	xmm0 = xmm1;			/* Copy I3 */ \
	xmm1 += xmm4;			/* R4 = I3 + R4 (final I3)	; 5-7 */ \
\
	xmm0 -= xmm4;			/* I3 = I3 - R4 (final I4)	; 6-8	avail 4,13+ */ \
	xmm13 = xptr(srcreg+d1+16+rbp);	/* Load I3-from-mem */ \
\
	xmm4 = xmm7;			/* Copy R3 */ \
	xmm7 -= xmm10;			/* R3 = R3 - I4 (final R3)	; 7-9 */ \
	xmm14 = xptr(srcreg+d1+rbp);		/* Load R3-from-mem */ \
\
	xmm10 += xmm4;			/* I4 = R3 + I4 (final R4)	; 8-10 */ \
	xmm15 = xmm1;			/* Copy I3 */ \
	xmm1 *= xmm13;			/* I3 * I3-from-mem		; 8-12	avail 4 */ \
\
	xmm6 += xmm5;			/* R4 = R2 + R4 (new R2)	; 9-11	avail 4,5 */ \
	xmm15 *= xmm14;			/* I3 * R3-from-mem		; 9-13 */ \
\
	xmm9 += xmm8;			/* R3 = R1 + R3 (new R1)	; 10-12	avail 4,5,8 */ \
	xmm14 *= xmm7;			/* R3 * R3-from-mem		; 10-14 */ \
	xmm5 = xptr(srcreg+d1+48+rbp);	/* Load I4-from-mem */ \
\
	xmm12 += xmm11;			/* I4 = I2 + I4 (new I2)	; 11-13	avail 4,8,11 */ \
	xmm7 *= xmm13;			/* R3 * I3-from-mem		; 11-15	avail 4,8,11,13 */ \
	xmm11 = xptr(srcreg+d1+32+rbp);	/* Load R4-from-mem */ \
\
	xmm3 += xmm2;			/* I3 = I1 + I3 (new I1)	; 12-14	avail 4,8,13,2 */ \
	xmm8 = xmm0;			/* Copy I4 */ \
	xmm0 *= xmm5;			/* I4 * I4-from-mem		; 12-16	avail 4,13,2 */ \
\
	xmm4 = xmm9;			/* Copy R1 */ \
	xmm9 -= xmm6;			/* R1 = R1 - R2 (final R2)	; 13-15	avail 13,2 */ \
	xmm8 *= xmm11;			/* I4 * R4-from-mem		; 13-17 */ \
\
	xmm6 += xmm4;			/* R2 = R1 + R2 (final R1)	; 14-16	avail 4,13,2 */ \
	xmm11 *= xmm10;			/* R4 * R4-from-mem		; 14-18 */ \
	xmm13 = xptr(srcreg+32+rbp);		/* Load R2-from-mem */ \
\
	xmm4 = xmm3;			/* Copy I1 */ \
	xmm3 -= xmm12;			/* I1 = I1 - I2 (final I2)	; 15-17	avail 2 */ \
	xmm10 *= xmm5;			/* R4 * I4-from-mem		; 15-19	avail 2,5 */ \
\
	xprefetchw(u8ptr(srcreg+srcinc)); \
\
	xmm12 += xmm4;			/* I2 = I1 + I2 (final I1)	; 16-18	avail 4,2,5 */ \
	xmm2 = xmm9;			/* Copy R2 */ \
	xmm9 *= xmm13;			/* R2 * R2-from-mem		; 16-20	avail 4,5 */ \
\
	xmm14 -= xmm1;			/* R3*R3-I3*I3 (new R3)		; 17-19	avail 4,5,1 */ \
	xmm5 = xptr(srcreg+48+rbp);		/* Load I2-from-mem */ \
	xmm2 *= xmm5;			/* R2 * I2-from-mem		; 17-21	avail 4,1 */ \
\
	xmm7 += xmm15;			/* R3*I3+I3*R3 (new I3)		; 18-20	avail 4,1,15 */ \
	xmm5 *= xmm3;			/* I2 * I2-from-mem		; 18-22 */ \
	xmm1 = xptr(srcreg+rbp);		/* Load R1-from-mem */ \
\
	xmm11 -= xmm0;			/* R4*R4-I4*I4 (new R4)		; 19-21	avail 4,15,0 */ \
	xmm3 *= xmm13;			/* I2 * R2-from-mem		; 19-23	avail 4,15,0,13 */ \
	xmm15 = xptr(srcreg+48);		/*#2 I1 */ \
\
	xmm10 += xmm8;			/* R4*I4+I4*R4 (new I4)		; 20-22	avail 4,0,13,8 */ \
	xmm0 = xmm6;			/* Copy R1 */ \
	xmm6 *= xmm1;			/* R1 * R1-from-mem		; 20-24	avail 4,13,8 */ \
\
	xprefetchw(u8ptr(srcreg+srcinc+d1)); \
\
	xmm13 = xptr(srcreg+d2+48);		/*#2 I3 */ \
	xmm8 = xmm15;			/*#2 Copy I1 */ \
	xmm15 -= xmm13;			/*#2 I1 = I1 - I3 (new I3)	; 21-23	avail 4 */ \
	xmm4 = xptr(srcreg+16+rbp);		/* Load I1-from-mem */ \
	xmm0 *= xmm4;			/* R1 * I1-from-mem		; 21-25	avail none */ \
\
	xmm8 += xmm13;			/*#2 I3 = I1 + I3 (new I1)	; 22-24 avail 13 */ \
	xmm4 *= xmm12;			/* I1 * I1-from-mem		; 22-26 */ \
\
	xmm9 -= xmm5;			/* R2*R2-I2*I2 (new R2)		; 23-25 avail 13,5 */ \
	xmm1 *= xmm12;			/* I1 * R1-from-mem		; 23-27 avail 13,5,12 */ \
\
	xmm2 += xmm3;			/* R2*I2+I2*R2 (new I2)		; 24-26 avail 13,5,12,3 */ \
\
	xmm5 = xmm11;			/* Copy R4 */ \
	xmm11 -= xmm14;			/* R4 = R4 - R3 (new I4)	; 25-27 avail 13,12,3 */ \
\
	xmm14 += xmm5;			/* R3 = R4 + R3 (new R3)	; 26-28 avail 5,12,3 */ \
	xmm13 = xptr(srcreg+d1+32);		/*#2 R2 */ \
\
	xmm6 -= xmm4;			/* R1*R1-I1*I1 (new R1)		; 27-29 avail 5,12,3,4 */ \
	xmm12 = xptr(srcreg+d2+d1+32);	/*#2 R4 */ \
\
	xmm0 += xmm1;			/* R1*I1+I1*R1 (new I1)		; 28-30 avail 5,3,4,1 */ \
	xmm3 = xptr(srcreg+32);		/*#2 R1 */ \
\
	xmm5 = xmm10;			/* Copy I4 */ \
	xmm10 += xmm7;			/* I4 = I3 + I4 (new I3)	; 29-31 avail 4,1 */ \
\
	xmm7 -= xmm5;			/* I3 = I3 - I4 (new R4)	; 30-32 avail 5,4,1 */ \
	xmm4 = xptr(srcreg+d2+32);		/*#2 R3 */ \
\
	xmm5 = xmm6;			/* Copy R1 */ \
	xmm6 -= xmm9;			/* R1 = R1 - R2 (new R2)	; 31-33	avail 1 */ \
\
	xmm9 += xmm5;			/* R2 = R1 + R2 (new R1)	; 32-34	avail 5,1 */ \
\
	xmm5 = xmm0;			/* Copy I1 */ \
	xmm0 -= xmm2;			/* I1 = I1 - I2 (new I2)	; 33-35	avail 1 */ \
\
	xmm2 += xmm5;			/* I2 = I1 + I2 (new I1)	; 34-36	avail 5,1 */ \
\
	xmm5 = xmm9;			/* Copy R1 */ \
	xmm9 -= xmm14;			/* R1 = R1 - R3 (final R3)	; 35-37	avail 1 storable 9 */ \
\
	xmm14 += xmm5;			/* R3 = R1 + R3 (final R1)	; 36-38	avail 5,1 storable 9,14 */ \
\
	xmm5 = xmm6;			/* Copy R2 */ \
	xmm6 -= xmm7;			/* R2 = R2 - R4 (final R4)	; 37-39	avail 1 storable 9,14,6 */ \
\
	xmm7 += xmm5;			/* R4 = R2 + R4 (final R2)	; 38-40	avail 5,1 storable 9,14,6,7 */ \
	xptr(srcreg+16) = xmm9;		/* Save R3			; 38 */ \
\
	xmm1 = xmm13;			/*#2 Copy R2 */ \
	xmm13 -= xmm12;			/*#2 R2 = R2 - R4 (new R4)	; 39-41	avail 5,9 storable 14,6,7 */ \
	xptr(srcreg) = xmm14;			/* Save R1			; 39 */ \
	xmm14 = xptr(srcreg+d1+48);		/*#2 I2 */ \
\
	xmm9 = xmm3;			/*#2 Copy R1 */ \
	xmm3 -= xmm4;			/*#2 R1 = R1 - R3 (new R3)	; 40-42	avail 5 storable 6,7 */ \
	xptr(srcreg+d1+16) = xmm6;		/* Save R4			; 40 */ \
\
	xmm5 = xptr(srcreg+d2+d1+48);		/*#2 I4 */ \
	xmm6 = xmm14;			/*#2 Copy I2 */ \
	xmm14 -= xmm5;			/*#2 I2 = I2 - I4 (new I4)	; 41-43	avail none storable 7 */ \
	xptr(srcreg+d1) = xmm7;		/* Save R2			; 41	avail 7 */ \
\
	xmm7 = xmm15;			/*#2 Copy I3 */ \
	xmm15 += xmm13;			/*#2 R4 = I3 + R4 (final I3)	; 42-44	avail none */ \
\
	xmm7 -= xmm13;			/*#2 I3 = I3 - R4 (final I4)	; 43-45	avail 13 */ \
\
	xmm13 = xmm3;			/*#2 Copy R3 */ \
	xmm3 -= xmm14;			/*#2 R3 = R3 - I4 (final R3)	; 44-46 */ \
\
	xmm14 += xmm13;			/*#2 I4 = R3 + I4 (final R4)	; 45-47 */ \
	xmm13 = xptr(srcreg+d2+d1+16+rbp);	/*#2 Load I3-from-mem */ \
\
	xmm1 += xmm12;			/*#2 R4 = R2 + R4 (new R2)	; 46-48	avail 12 */ \
	xmm12 = xmm15;			/*#2 Copy I3 */ \
	xmm15 *= xmm13;			/*#2 I3 * I3-from-mem		; 46-50	avail none */ \
\
	xmm4 += xmm9;			/*#2 R3 = R1 + R3 (new R1)	; 47-49	avail 9 */ \
	xmm9 = xptr(srcreg+d2+d1+rbp);	/*#2 Load R3-from-mem */ \
	xmm12 *= xmm9;			/*#2 I3 * R3-from-mem		; 47-51	avail none */ \
\
	xmm6 += xmm5;			/*#2 I4 = I2 + I4 (new I2)	; 48-50	avail 5 */ \
	xmm9 *= xmm3;			/*#2 R3 * R3-from-mem		; 48-52 */ \
\
	xmm5 = xmm0;			/* Copy I2 */ \
	xmm0 -= xmm11;			/* I2 = I2 - I4 (final I4)	; 49-51	avail none storable 0 */ \
	xmm3 *= xmm13;			/*#2 R3 * I3-from-mem		; 49-53	avail 13 storable 0 */ \
	xmm13 = xptr(srcreg+d2+d1+48+rbp);	/*#2 Load I4-from-mem */ \
\
	xmm11 += xmm5;			/* I4 = I2 + I4 (final I2)	; 50-52	avail 5 storable 0,11 */ \
	xmm5 = xmm7;			/*#2 Copy I4 */ \
	xmm7 *= xmm13;			/*#2 I4 * I4-from-mem		; 50-54	avail none storable 0,11 */ \
\
	xptr(srcreg+d1+48) = xmm0;		/* Save I4			; 52 */ \
	xmm0 = xmm4;			/*#2 Copy R1 */ \
	xmm4 -= xmm1;			/*#2 R1 = R1 - R2 (final R2)	; 51-53	avail none storable 11 */ \
	xptr(srcreg+d1+32) = xmm11;		/* Save I2			; 53 */ \
	xmm11 = xptr(srcreg+d2+d1+32+rbp);	/*#2 Load R4-from-mem */ \
	xmm5 *= xmm11;			/*#2 I4 * R4-from-mem		; 51-55	avail none */ \
\
	xmm1 += xmm0;			/*#2 R2 = R1 + R2 (final R1)	; 52-54	avail 0 */ \
	xmm11 *= xmm14;			/*#2 R4 * R4-from-mem		; 52-56 */ \
\
	xprefetchw(u8ptr(srcreg+srcinc+d2)); \
\
	xmm0 = xmm8;			/*#2 Copy I1 */ \
	xmm8 -= xmm6;			/*#2 I1 = I1 - I2 (final I2)	; 53-55	avail none */ \
	xmm14 *= xmm13;			/*#2 R4 * I4-from-mem		; 53-57	avail 13 */ \
\
	xmm6 += xmm0;			/*#2 I2 = I1 + I2 (final I1)	; 54-56	avail 13,0 */ \
	xmm13 = xptr(srcreg+d2+32+rbp);	/*#2 Load R2-from-mem */ \
	xmm0 = xmm4;			/*#2 Copy R2 */ \
	xmm4 *= xmm13;			/*#2 R2 * R2-from-mem		; 54-58	avail none */ \
\
	xmm9 -= xmm15;			/*#2 R3*R3-I3*I3 (new R3)	; 55-57	avail 15 */ \
	xmm15 = xptr(srcreg+d2+48+rbp);	/*#2 Load I2-from-mem */ \
	xmm0 *= xmm15;			/*#2 R2 * I2-from-mem		; 55-59	avail none */ \
\
	xmm3 += xmm12;			/*#2 R3*I3+I3*R3 (new I3)	; 56-58	avail 12 */ \
	xmm15 *= xmm8;			/*#2 I2 * I2-from-mem		; 56-60 */ \
\
	xmm11 -= xmm7;			/*#2 R4*R4-I4*I4 (new R4)	; 57-59	avail 12,7 */ \
	xmm8 *= xmm13;			/*#2 I2 * R2-from-mem		; 57-61	avail 12,7,13 */ \
	xmm7 = xptr(srcreg+d2+rbp);		/*#2 Load R1-from-mem */ \
\
	xprefetchw(u8ptr(srcreg+srcinc+d2+d1)); \
\
	xmm14 += xmm5;			/*#2 R4*I4+I4*R4 (new I4)	; 58-60	avail 12,13,5 */ \
	xmm13 = xmm1;			/*#2 Copy R1 */ \
	xmm1 *= xmm7;			/*#2 R1 * R1-from-mem		; 58-62	avail 12,5 */ \
\
	xmm5 = xmm2;			/* Copy I1 */ \
	xmm2 -= xmm10;			/* I1 = I1 - I3 (final I3)	; 59-61	avail 12 storable 2 */ \
	xmm12 = xptr(srcreg+d2+16+rbp);	/*#2 Load I1-from-mem	 */ \
	xmm13 *= xmm12;			/*#2 R1 * I1-from-mem		; 59-63	avail none */ \
\
	xmm10 += xmm5;			/* I3 = I1 + I3 (final I1)	; 60-62	avail 5 storable 2,10 */ \
	xmm12 *= xmm6;			/*#2 I1 * I1-from-mem		; 60-64 */ \
\
	xmm4 -= xmm15;			/*#2 R2*R2-I2*I2 (new R2)	; 61-63	avail 5,15 storable 2,10 */ \
	xmm6 *= xmm7;			/*#2 I1 * R1-from-mem		; 61-65	avail 5,15,7 storable 2,10 */ \
\
	xmm0 += xmm8;			/*#2 R2*I2+I2*R2 (new I2)	; 62-64	avail 5,15,7,8 storable 2,10 */ \
	xptr(srcreg+48) = xmm2;		/* Save I3			; 62 */ \
\
	xmm5 = xmm11;			/*#2 Copy R4 */ \
	xmm11 -= xmm9;			/*#2 R4 = R4 - R3 (new I4)	; 63-65	avail 15,7,8,2 storable 10 */ \
	xmm9 += xmm5;			/*#2 R3 = R4 + R3 (new R3)	; 64-66	avail 5,15,7,8,2 storable 10 */ \
	xptr(srcreg+32) = xmm10;		/* Save I1			; 63 */ \
\
	xmm1 -= xmm12;			/*#2 R1*R1-I1*I1 (new R1)	; 65-67	avail 5,15,7,8,2,10,12 */ \
\
	xmm13 += xmm6;			/*#2 R1*I1+I1*R1 (new I1)	; 66-68	avail 5,15,7,8,2,10,12,6 */ \
\
	xmm5 = xmm14;			/*#2 Copy I4 */ \
	xmm14 += xmm3;			/*#2 I4 = I3 + I4 (new I3)	; 67-69 */ \
	xmm3 -= xmm5;			/*#2 I3 = I3 - I4 (new R4)	; 68-70 */ \
\
	xmm5 = xmm1;			/*#2 Copy R1 */ \
	xmm1 -= xmm4;			/*#2 R1 = R1 - R2 (new R2)	; 69-71 */ \
	xmm4 += xmm5;			/*#2 R2 = R1 + R2 (new R1)	; 70-72 */ \
\
	xmm5 = xmm13;			/*#2 Copy I1 */ \
	xmm13 -= xmm0;			/*#2 I1 = I1 - I2 (new I2)	; 71-73 */ \
	xmm0 += xmm5;			/*#2 I2 = I1 + I2 (new I1)	; 72-74 */ \
\
	xmm5 = xmm1;			/*#2 Copy R2 */ \
	xmm1 -= xmm3;			/*#2 R2 = R2 - R4 (final R4)	; 73-75 */ \
	xmm3 += xmm5;			/*#2 R4 = R2 + R4 (final R2)	; 74-76 */ \
\
	xmm5 = xmm4;			/*#2 Copy R1 */ \
	xmm4 -= xmm9;			/*#2 R1 = R1 - R3 (final R3)	; 75-77 */ \
	xmm9 += xmm5;			/*#2 R3 = R1 + R3 (final R1)	; 76-78 */ \
	xptr(srcreg+d2+d1+16) = xmm1;		/*#2 Save R4			; 76 */ \
\
	xmm5 = xmm13;			/*#2 Copy I2 */ \
	xmm13 -= xmm11;			/*#2 I2 = I2 - I4 (final I4)	; 77-79 */ \
	xptr(srcreg+d2+d1) = xmm3;		/*#2 Save R2			; 77 */ \
	xmm11 += xmm5;			/*#2 I4 = I2 + I4 (final I2)	; 78-80 */ \
	xptr(srcreg+d2+16) = xmm4;		/*#2 Save R3			; 78 */ \
\
	xmm5 = xmm0;			/*#2 Copy I1 */ \
	xmm0 -= xmm14;			/*#2 I1 = I1 - I3 (final I3)	; 79-81 */ \
	xptr(srcreg+d2) = xmm9;		/*#2 Save R1			; 79 */ \
	xmm14 += xmm5;			/*#2 I3 = I1 + I3 (final I1)	; 80-82 */ \
	xptr(srcreg+d2+d1+48) = xmm13;	/*#2 Save I4			; 80 */ \
\
	xptr(srcreg+d2+d1+32) = xmm11;	/*#2 Save I2			; 81 */ \
	xptr(srcreg+d2+48) = xmm0;		/*#2 Save I3			; 82 */ \
	xptr(srcreg+d2+32) = xmm14;		/*#2 Save I1			; 83 */ \
	srcreg += srcinc; \
\
\

#define r4_x4cl_four_complex_with_mulf(srcreg,srcinc,d1,d2) { \
	vec2f64 xmm8,xmm9,xmm10,xmm11,xmm12,xmm13,xmm14,xmm15; \
	xmm1 = xptr(srcreg+d1+16+rbx);	/* I3 */ \
	xmm2 = xptr(srcreg+d1+16+rbp);	/* Load I3-other */ \
	xmm3 = xmm1;			/* Copy I3 */ \
	xmm1 *= xmm2;			/* I3 * I3-other		; 1-5 */ \
\
	xmm4 = xptr(srcreg+d1+rbp);		/* Load R3-other */ \
	xmm3 *= xmm4;			/* I3 * R3-other		; 2-6 */ \
\
	xmm5 = xptr(srcreg+d1+rbx);		/* R3 */ \
	xmm4 *= xmm5;			/* R3 * R3-other		; 3-7 */ \
\
	xmm5 *= xmm2;			/* R3 * I3-other		; 4-8 */ \
	xmm6 = xptr(srcreg+d1+48+rbx);	/* I4 */ \
\
	xmm7 = xptr(srcreg+d1+48+rbp);	/* Load I4-other */ \
	xmm2 = xmm6;			/* Copy I4 */ \
	xmm6 *= xmm7;			/* I4 * I4-other		; 5-9 */ \
\
	xmm0 = xptr(srcreg+d1+32+rbp);	/* Load R4-other */ \
	xmm2 *= xmm0;			/* I4 * R4-other		; 6-10 */ \
\
	xmm8 = xptr(srcreg+d1+32+rbx);	/* R4 */ \
	xmm0 *= xmm8;			/* R4 * R4-other		; 7-11 */ \
\
	xmm4 -= xmm1;			/* R3*R3-I3*I3 (new R3)		; 8-10	avail 1,9+ */ \
	xmm8 *= xmm7;			/* R4 * I4-other		; 8-12	avail 1,7,9+ */ \
	xmm9 = xptr(srcreg+32+rbx);		/* R2 */ \
\
	xmm5 += xmm3;			/* R3*I3+I3*R3 (new I3)		; 9-11	avail 1,7,3,9+ */ \
	xmm10 = xptr(srcreg+32+rbp);		/* Load R2-other */ \
	xmm11 = xmm9;			/* Copy R2 */ \
	xmm9 *= xmm10;			/* R2 * R2-other		; 9-13	avail 1,7,3,12+ */ \
\
	xprefetchw(u8ptr(srcreg+srcinc)); \
\
	xmm3 = xptr(srcreg+48+rbp);		/* Load I2-other */ \
	xmm11 *= xmm3;			/* R2 * I2-other		; 10-14	avail 1,7,12+ */ \
\
	xmm7 = xptr(srcreg+48+rbx);		/* I2 */ \
	xmm3 *= xmm7;			/* I2 * I2-other		; 11-15	avail 1,12+ */ \
\
	xmm0 -= xmm6;			/* R4*R4-I4*I4 (new R4)		; 12-14	avail 1,6,12+ */ \
	xmm7 *= xmm10;			/* I2 * R2-other		; 12-16	avail 1,6,10,12+ */ \
	xmm6 = xptr(srcreg+rbx);		/* R1 */ \
\
	xmm8 += xmm2;			/* R4*I4+I4*R4 (new I4)		; 13-15	avail 1,2,10,12+ */ \
	xmm1 = xptr(srcreg+rbp);		/* Load R1-other */ \
	xmm2 = xmm6;			/* Copy R1 */ \
	xmm6 *= xmm1;			/* R1 * R1-other		; 13-17	avail 10,12+ */ \
\
	xprefetchw(u8ptr(srcreg+srcinc+d1)); \
\
	xmm10 = xptr(srcreg+16+rbp);		/* Load I1-other */ \
	xmm2 *= xmm10;			/* R1 * I1-other		; 14-18	avail 12+ */ \
\
	xmm12 = xptr(srcreg+16+rbx);		/* I1 */ \
	xmm10 *= xmm12;			/* I1 * I1-other		; 15-19	avail 13+ */ \
\
	xmm9 -= xmm3;			/* R2*R2-I2*I2 (new R2)		; 16-18	avail 3,13+ */ \
	xmm12 *= xmm1;			/* I1 * R1-other		; 16-20	avail 3,1,13+ */ \
\
	xmm11 += xmm7;			/* R2*I2+I2*R2 (new I2)		; 17-19 avail 3,1,7,13+ */ \
\
	xmm7 = xmm0;			/* Copy R4 */ \
	xmm0 -= xmm4;			/* R4 = R4 - R3 (new I4)	; 18-20 */ \
	xmm4 += xmm7;			/* R3 = R4 + R3 (new R3)	; 19-21 */ \
\
	xmm6 -= xmm10;			/* R1*R1-I1*I1 (new R1)		; 20-22 avail 3,1,7,10,13+ */ \
	xmm13 = xptr(srcreg+d2+d1+16+rbx);	/*#2 I3 */ \
\
	xmm2 += xmm12;			/* R1*I1+I1*R1 (new I1)		; 21-23 avail 3,1,7,10,12,14+ */ \
	xmm14 = xptr(srcreg+d2+d1+16+rbp);	/*#2 Load I3-other */ \
\
	xmm7 = xmm8;			/* Copy I4 */ \
	xmm8 += xmm5;			/* I4 = I3 + I4 (new I3)	; 22-24 avail 3,1,7,10,12,15 */ \
\
	xmm5 -= xmm7;			/* I3 = I3 - I4 (new R4)	; 23-25 */ \
	xmm10 = xptr(srcreg+d2+d1+rbp);	/*#2 Load R3-other */ \
\
	xprefetchw(u8ptr(srcreg+srcinc+d2)); \
\
	xmm7 = xmm6;			/* Copy R1 */ \
	xmm6 -= xmm9;			/* R1 = R1 - R2 (new R2)	; 24-26 */ \
	xmm1 = xptr(srcreg+d2+d1+rbx);	/*#2 R3 */ \
\
	xmm9 += xmm7;			/* R2 = R1 + R2 (new R1)	; 25-27 */ \
	xmm12 = xmm13;			/*#2 Copy I3 */ \
	xmm13 *= xmm14;			/*#2 I3 * I3-other		; 25-29 avail 3,7,15 */ \
\
	xmm7 = xmm2;			/* Copy I1 */ \
	xmm2 -= xmm11;			/* I1 = I1 - I2 (new I2)	; 26-28 avail 3,15 */ \
	xmm12 *= xmm10;			/*#2 I3 * R3-other		; 26-30 avail 3,15 */ \
\
	xprefetchw(u8ptr(srcreg+srcinc+d2+d1)); \
\
	xmm11 += xmm7;			/* I2 = I1 + I2 (new I1)	; 27-29 avail 3,7,15 */ \
	xmm10 *= xmm1;			/*#2 R3 * R3-other		; 27-31 */ \
	xmm15 = xptr(srcreg+d2+d1+48+rbx);	/*#2 I4 */ \
\
	xmm7 = xmm9;			/* Copy R1 */ \
	xmm9 -= xmm4;			/* R1 = R1 - R3 (final R3)	; 28-30 avail 3 storable 9 */ \
	xmm1 *= xmm14;			/*#2 R3 * I3-other		; 28-32 avail 3,14 storable 9 */ \
\
	xmm4 += xmm7;			/* R3 = R1 + R3 (final R1)	; 29-31 avail 3,7,14 storable 9,4 */ \
	xmm3 = xptr(srcreg+d2+d1+48+rbp);	/*#2 Load I4-other */ \
	xmm14 = xmm15;			/*#2 Copy I4 */ \
	xmm15 *= xmm3;			/*#2 I4 * I4-other		; 29-33 avail 7 storable 9,4 */ \
\
	xmm7 = xmm6;			/* Copy R2 */ \
	xmm6 -= xmm5;			/* R2 = R2 - R4 (final R4)	; 30-32 avail none storable 9,4,6 */ \
	xptr(srcreg+16) = xmm9;		/* Save R3			; 31 */ \
	xmm9 = xptr(srcreg+d2+d1+32+rbp);	/*#2 Load R4-other */ \
	xmm14 *= xmm9;			/*#2 I4 * R4-other		; 30-34	avail none storable 4,6 */ \
\
	xmm5 += xmm7;			/* R4 = R2 + R4 (final R2)	; 31-33	avail 7 storable 4,6,5 */ \
	xmm7 = xptr(srcreg+d2+d1+32+rbx);	/*#2 R4 */ \
	xmm9 *= xmm7;			/*#2 R4 * R4-other		; 31-35	avail none storable 4,6,5 */ \
\
	xmm10 -= xmm13;			/*#2 R3*R3-I3*I3 (new R3)	; 32-34	avail 13 storable 4,6,5 */ \
	xmm7 *= xmm3;			/*#2 R4 * I4-other		; 32-36	avail 13,3 storable 4,6,5 */ \
	xmm13 = xptr(srcreg+d2+32+rbx);	/*#2 R2 */ \
	xptr(srcreg) = xmm4;			/* Save R1			; 32	avail 3,4 storable 6,5 */ \
\
	xmm1 += xmm12;			/*#2 R3*I3+I3*R3 (new I3)	; 33-35	avail 3,4,12 storable 6,5 */ \
	xmm3 = xptr(srcreg+d2+32+rbp);	/*#2 Load R2-other */ \
	xmm4 = xmm13;			/*#2 Copy R2 */ \
	xmm13 *= xmm3;			/*#2 R2 * R2-other		; 33-37	avail 12 storable 6,5 */ \
	xptr(srcreg+d1+16) = xmm6;		/* Save R4			; 33	avail 6,12 storable 5 */ \
\
	xmm6 = xmm2;			/* Copy I2 */ \
	xmm2 -= xmm0;			/* I2 = I2 - I4 (final I4)	; 34-36	avail 12 storable 5,2 */ \
	xmm12 = xptr(srcreg+d2+48+rbp);	/*#2 Load I2-other */ \
	xmm4 *= xmm12;			/*#2 R2 * I2-other		; 34-38	avail none storable 5,2 */ \
	xptr(srcreg+d1) = xmm5;		/* Save R2			; 34	avail 5 storable 2 */ \
\
	xmm0 += xmm6;			/* I4 = I2 + I4 (final I2)	; 35-37	avail 5,6 storable 2,0 */ \
	xmm5 = xptr(srcreg+d2+48+rbx);	/*#2 I2 */ \
	xmm12 *= xmm5;			/*#2 I2 * I2-other		; 35-39	avail 6 storable 2,0 */ \
\
	xmm9 -= xmm15;			/*#2 R4*R4-I4*I4 (new R4)	; 36-38	avail 6,15 storable 2,0 */ \
	xmm5 *= xmm3;			/*#2 I2 * R2-other		; 36-40	avail 3,6,15 storable 2,0 */ \
	xmm6 = xptr(srcreg+d2+rbx);		/*#2 R1 */ \
\
	xmm7 += xmm14;			/*#2 R4*I4+I4*R4 (new I4)	; 37-39	avail 3,15,14 storable 2,0 */ \
	xmm15 = xptr(srcreg+d2+rbp);		/*#2 Load R1-other */ \
	xmm3 = xmm6;			/*#2 Copy R1 */ \
	xmm6 *= xmm15;			/*#2 R1 * R1-other		; 37-41	avail 14 storable 2,0 */ \
	xptr(srcreg+d1+48) = xmm2;		/* Save I4			; 37	avail 2,14 storable 0 */ \
\
	xmm2 = xmm11;			/* Copy I1 */ \
	xmm11 -= xmm8;			/* I1 = I1 - I3 (final I3)	; 38-40	avail 14 storable 0,11 */ \
	xmm14 = xptr(srcreg+d2+16+rbp);	/*#2 Load I1-other	 */ \
	xmm3 *= xmm14;			/*#2 R1 * I1-other		; 38-42	avail none storable 0,11 */ \
	xptr(srcreg+d1+32) = xmm0;		/* Save I2			; 38	avail 0 storable 11 */ \
\
	xmm8 += xmm2;			/* I3 = I1 + I3 (final I1)	; 39-41	avail 0,2 storable 11,8 */ \
	xmm0 = xptr(srcreg+d2+16+rbx);	/*#2 I1 */ \
	xmm14 *= xmm0;			/*#2 I1 * I1-other		; 39-43	avail 2 storable 11,8 */ \
\
	xmm13 -= xmm12;			/*#2 R2*R2-I2*I2 (new R2)	; 40-42	avail 2,12 storable 11,8 */ \
	xmm0 *= xmm15;			/*#2 I1 * R1-other		; 40-44	avail 2,12,15 storable 11,8 */ \
\
	xmm4 += xmm5;			/*#2 R2*I2+I2*R2 (new I2)	; 41-43	avail 2,12,15,5 storable 11,8 */ \
	xptr(srcreg+48) = xmm11;		/* Save I3			; 41	avail 2,12,15,5,11 storable 8 */ \
\
	xmm5 = xmm9;			/*#2 Copy R4 */ \
	xmm9 -= xmm10;			/*#2 R4 = R4 - R3 (new I4)	; 42-44 */ \
	xmm10 += xmm5;			/*#2 R3 = R4 + R3 (new R3)	; 43-45 */ \
	xptr(srcreg+32) = xmm8;		/* Save I1			; 42	avail 2,12,15,5,11,8 */ \
\
	xmm6 -= xmm14;			/*#2 R1*R1-I1*I1 (new R1)	; 44-46	avail 2,12,15,5,11,8,14 */ \
\
	xmm3 += xmm0;			/*#2 R1*I1+I1*R1 (new I1)	; 45-47	avail 2,12,15,5,11,8,14,0 */ \
\
	xmm0 = xmm7;			/*#2 Copy I4 */ \
	xmm7 += xmm1;			/*#2 I4 = I3 + I4 (new I3)	; 46-48 */ \
	xmm1 -= xmm0;			/*#2 I3 = I3 - I4 (new R4)	; 47-49 */ \
\
	xmm0 = xmm6;			/*#2 Copy R1 */ \
	xmm6 -= xmm13;			/*#2 R1 = R1 - R2 (new R2)	; 48-50 */ \
	xmm13 += xmm0;			/*#2 R2 = R1 + R2 (new R1)	; 49-51 */ \
\
	xmm0 = xmm3;			/*#2 Copy I1 */ \
	xmm3 -= xmm4;			/*#2 I1 = I1 - I2 (new I2)	; 50-52 */ \
	xmm4 += xmm0;			/*#2 I2 = I1 + I2 (new I1)	; 51-53 */ \
\
	xmm0 = xmm6;			/*#2 Copy R2 */ \
	xmm6 -= xmm1;			/*#2 R2 = R2 - R4 (final R4)	; 52-54 */ \
	xmm1 += xmm0;			/*#2 R4 = R2 + R4 (final R2)	; 53-55 */ \
\
	xmm0 = xmm13;			/*#2 Copy R1 */ \
	xmm13 -= xmm10;			/*#2 R1 = R1 - R3 (final R3)	; 54-56 */ \
	xmm10 += xmm0;			/*#2 R3 = R1 + R3 (final R1)	; 55-57 */ \
	xptr(srcreg+d2+d1+16) = xmm6;		/*#2 Save R4			; 55 */ \
\
	xmm0 = xmm3;			/*#2 Copy I2 */ \
	xmm3 -= xmm9;			/*#2 I2 = I2 - I4 (final I4)	; 56-58 */ \
	xptr(srcreg+d2+d1) = xmm1;		/*#2 Save R2			; 56 */ \
	xmm9 += xmm0;			/*#2 I4 = I2 + I4 (final I2)	; 57-59 */ \
	xptr(srcreg+d2+16) = xmm13;		/*#2 Save R3			; 57 */ \
\
	xmm0 = xmm4;			/*#2 Copy I1 */ \
	xmm4 -= xmm7;			/*#2 I1 = I1 - I3 (final I3)	; 58-60 */ \
	xptr(srcreg+d2) = xmm10;		/*#2 Save R1			; 58 */ \
	xmm7 += xmm0;			/*#2 I3 = I1 + I3 (final I1)	; 59-61 */ \
	xptr(srcreg+d2+d1+48) = xmm3;		/*#2 Save I4			; 59 */ \
\
	xptr(srcreg+d2+d1+32) = xmm9;		/*#2 Save I2			; 60 */ \
	xptr(srcreg+d2+48) = xmm4;		/*#2 Save I3			; 61 */ \
	xptr(srcreg+d2+32) = xmm7;		/*#2 Save I1			; 62 */ \
	srcreg += srcinc; \
}
#endif
#endif

/* */ \
/* ************************************* four-complex-fft4 variants ****************************************** */ \
/* */ \
/* These macros are used in the last levels of pass 1.  Four sin/cos multipliers are needed to */ \
/* finish off the partial sin/cos multiplies that were done in the first levels of pass 1. */ \
/* FFTs of type r4delay do this to reduce memory usage at the cost of some extra */ \
/* complex multiplies. */ \
\
\
/* */ \
/* In the split premultiplier case, we apply part of the roots of -1 at the */ \
/* end of the first pass.  Thus we have 4 sin/cos/premultipliers instead */ \
/* of the usual 3. */ \
/* */ \
\
/* Used in last levels of pass 1 (split premultiplier and delay cases).  No swizzling. */
#if 0 //IFDEF UNUSED

#define r4_g4cl_four_complex_fft4(srcreg,srcinc,d1,d2,dstreg,dstinc,e1,e2,screg) \
	xprefetch(u8ptr(srcreg+srcinc)); \
	r4_x4c_fft4_mem(xptr(srcreg),xptr(srcreg+d1),xptr(srcreg+d2),xptr(srcreg+d2+d1),xptr(srcreg+16),xptr(srcreg+d1+16),xptr(srcreg+d2+16),xptr(srcreg+d2+d1+16),xptr(dstreg),xptr(dstreg+16),screg,0,dstreg+dstinc,e1); \
	xprefetch(u8ptr(srcreg+srcinc+d1)); \
	xptr(dstreg) = xmm2;		/* Save R1 */ \
	xptr(dstreg+16) = xmm5;	/* Save I1 */ \
	xmm2 = xptr(srcreg+32);	/* R1 */ \
	xmm5 = xptr(srcreg+48);	/* R5 */ \
	xptr(dstreg+32) = xmm7;	/* Save R2 */ \
	xptr(dstreg+48) = xmm6;	/* Save I2 */ \
	xmm7 = xptr(srcreg+d1+32);	/* R2 */ \
	xmm6 = xptr(srcreg+d1+48);	/* R6 */ \
	xptr(dstreg+e1) = xmm3;	/* Save R3 */ \
	xptr(dstreg+e1+16) = xmm1;	/* Save I3 */ \
	xptr(dstreg+e1+32) = xmm0;	/* Save R4 */ \
	xptr(dstreg+e1+48) = xmm4;	/* Save I4 */ \
	xprefetch(u8ptr(srcreg+srcinc+d2)); \
	r4_x4c_fft4_partial_mem(xmm2,xmm7,xmm3,xmm0,xmm5,xmm6,xmm1,xmm4,xptr(srcreg+d2+32),xptr(srcreg+d2+d1+32),xptr(srcreg+d2+48),xptr(srcreg+d2+d1+48),xptr(dstreg+e2),xptr(dstreg+e2+16),screg,0,dstreg+dstinc+e2,e1); \
	xprefetch(u8ptr(srcreg+srcinc+d2+d1)); \
	srcreg += srcinc; \
	xptr(dstreg+e2) = xmm3;	/* Save R1 */ \
	xptr(dstreg+e2+16) = xmm6;	/* Save I1 */ \
	xptr(dstreg+e2+32) = xmm4;	/* Save R2 */ \
	xptr(dstreg+e2+48) = xmm1;	/* Save I2 */ \
	xptr(dstreg+e2+e1) = xmm0;	/* Save R3 */ \
	xptr(dstreg+e2+e1+16) = xmm7;	/* Save I3 */ \
	xptr(dstreg+e2+e1+32) = xmm2;	/* Save R4 */ \
	xptr(dstreg+e2+e1+48) = xmm5;	/* Save I4 */ \
	dstreg += dstinc; \
\
\

#define r4_x4c_fft4_partial_mem(r1,r2,r3,r4,r5,r6,r7,r8,mem3,mem4,mem7,mem8,dest1,dest2,screg,off,pre1,pre2) \
	r3 = mem3;			/* R3 */ \
	r3 += r1;				/* R3 = R1 + R3 (new R1) */ \
	r1 -= mem3;			/* R1 = R1 - R3 (new R3) */ \
\
	r4 = mem4;			/* R4 */ \
	r4 += r2;				/* R4 = R2 + R4 (new R2) */ \
	r2 -= mem4;			/* R2 = R2 - R4 (new R4) */ \
\
	r7 = r3; \
	r3 -= r4;				/* R1 = R1 - R2 (final R2) */ \
	r4 += r7;				/* R2 = R1 + R2 (final R1) */ \
\
	r7 = mem7;			/* I3 */ \
	r7 += r5;				/* I3 = I1 + I3 (new I1) */ \
	r5 -= mem7;			/* I1 = I1 - I3 (new I3) */ \
\
	r8 = mem8;			/* I4 */ \
	r8 += r6;				/* I4 = I2 + I4 (new I2) */ \
	r6 -= mem8;			/* I2 = I2 - I4 (new I4) */ \
\
	dest1 = r4;			/* Save R1 */ \
\
	r4 = r5; \
	r5 -= r2;				/* I3 = I3 - R4 (final I4) */ \
	r2 += r4;				/* R4 = I3 + R4 (final I3) */ \
\
	r4 = r1; \
	r1 -= r6;				/* R3 = R3 - I4 (final R3) */ \
	r6 += r4;				/* I4 = R3 + I4 (final R4) */ \
\
	xprefetchw(u8ptr(pre1)); \
\
	r4 = r7; \
	r7 -= r8;				/* I1 = I1 - I2 (final I2) */ \
	r8 += r4;				/* I2 = I1 + I2 (final I1) */ \
\
	r4 = xptr(screg+off+32+16);		/* cosine/sine */ \
	r4 *= r1;				/* A3 = R3 * cosine/sine */ \
	r4 -= r2;				/* A3 = A3 - I3 */ \
	r2 *= xptr(screg+off+32+16);		/* B3 = I3 * cosine/sine */ \
	r2 += r1;				/* B3 = B3 + R3 */ \
\
	dest2 = r8;			/* Save I1 */ \
\
	r8 = xptr(screg+off+64+16);		/* cosine/sine */ \
	r8 *= r3;				/* A2 = R2 * cosine/sine */ \
	r8 -= r7;				/* A2 = A2 - I2 */ \
	r7 *= xptr(screg+off+64+16);		/* B2 = I2 * cosine/sine */ \
	r7 += r3;				/* B2 = B2 + R2 */ \
\
	xprefetchw(u8ptr(pre1+pre2)); \
\
	r1 = xptr(screg+off+96+16);		/* cosine/sine */ \
	r1 *= r6;				/* A4 = R4 * cosine/sine */ \
	r1 -= r5;				/* A4 = A4 - I4 */ \
	r5 *= xptr(screg+off+96+16);		/* B4 = I4 * cosine/sine */ \
	r5 += r6;				/* B4 = B4 + R4 */ \
\
	r3 = xptr(screg+off+32); \
	r4 *= r3;				/* A3 = A3 * sine (final R3) */ \
	r2 *= r3;				/* B3 = B3 * sine (final I3) */ \
	r6 = xptr(screg+off+64); \
	r8 *= r6;				/* A2 = A2 * sine (final R2) */ \
	r7 *= r6;				/* B2 = B2 * sine (final I2) */ \
\
	r3 = xptr(screg+off+0+16);		/* cosine/sine */ \
	r3 *= dest1;			/* A1 = R1 * cosine/sine */ \
	r6 = dest2;			/* Restore I1 */ \
	r3 -= r6;				/* A1 = A1 - I1 */ \
	r6 *= xptr(screg+off+0+16);		/* B1 = I1 * cosine/sine */ \
	r6 += dest1;			/* B1 = B1 + R1 */ \
\
	r1 *= xptr(screg+off+96);		/* A4 = A4 * sine (final R4) */ \
	r5 *= xptr(screg+off+96);		/* B4 = B4 * sine (final I4) */ \
	r3 *= xptr(screg+off+0);		/* A1 = A1 * sine (final R1) */ \
	r6 *= xptr(screg+off+0);		/* B1 = B1 * sine (final I1) */ \

#endif
\
/* Used in last levels of pass 1 (split premultiplier and delay cases).  Swizzling. */ \

#ifdef LOW_CPU_REGISTER_COUNT
#define r4_sg4cl_four_complex_fft4(srcreg,srcinc,d1,d2,dstreg,dstinc,e1,e2,screg) { \
	vec2f64 *sd0 = (vec2f64*)(srcreg); \
	vec2f64 *sd1 = (vec2f64*)(srcreg+d1); \
	vec2f64 *sd2 = (vec2f64*)(srcreg+d2); \
	vec2f64 *sd3 = (vec2f64*)(srcreg+d2+d1); \
	vec2f64 *de0 = (vec2f64*)(dstreg); \
	vec2f64 *de1 = (vec2f64*)(dstreg+e1); \
	vec2f64 *de2 = (vec2f64*)(dstreg+e2); \
	vec2f64 *de3 = (vec2f64*)(dstreg+e2+e1); \
	vec2f64 *sc = (vec2f64*)(screg); \
	xprefetch(u8ptr(srcreg+srcinc)); \
	r4_x4c_fft4_mem8(sd0[0],sd1[0],sd2[0],sd3[0],sd0[1],sd1[1],sd2[1],sd3[1], de0[0],de0[1], sc[0],sc[1],sc[2],sc[3],sc[4],sc[5],sc[6],sc[7], dstreg+dstinc,e1); \
	xprefetch(u8ptr(srcreg+srcinc+d1)); \
	shuffle_store(de0[0], de0[1], xmm2, xmm7);					/* Save R1,R2 */ \
	shuffle_store_with_temp(de0[2], de0[3], xmm5, xmm6, xmm7);		/* Save I1,I2 */ \
	shuffle_store_with_temp(de1[0], de1[1], xmm3, xmm0, xmm7);		/* Save R3,R4 */ \
	shuffle_store_with_temp(de1[2], de1[3], xmm1, xmm4, xmm7);	/* Save I3,I4 */ \
	xprefetch(u8ptr(srcreg+srcinc+d2)); \
	r4_x4c_fft4_mem8(sd0[2],sd1[2],sd2[2],sd3[2],sd0[3],sd1[3],sd2[3],sd3[3], de2[0],de2[1], sc[0],sc[1],sc[2],sc[3],sc[4],sc[5],sc[6],sc[7], dstreg+dstinc+e2,e1); \
	xprefetch(u8ptr(srcreg+srcinc+d2+d1)); \
	srcreg += srcinc; \
	shuffle_store(de2[0], de2[1], xmm2, xmm7);				/* Save R1,R2 */ \
	shuffle_store_with_temp(de2[2], de2[3], xmm5, xmm6, xmm7);	/* Save I1,I2 */ \
	shuffle_store_with_temp(de3[0], de3[1], xmm3, xmm0, xmm7);	/* Save R3,R4 */ \
	shuffle_store_with_temp(de3[2], de3[3], xmm1, xmm4, xmm7);	/* Save I3,I4 */ \
	dstreg += dstinc; \
}
#else
#define r4_sg4cl_four_complex_fft4(srcreg,srcinc,d1,d2,dstreg,dstinc,e1,e2,screg) { \
	vec2f64 *sd0 = (vec2f64*)(srcreg); \
	vec2f64 *sd1 = (vec2f64*)(srcreg+d1); \
	vec2f64 *sd2 = (vec2f64*)(srcreg+d2); \
	vec2f64 *sd3 = (vec2f64*)(srcreg+d2+d1); \
	vec2f64 *de0 = (vec2f64*)(dstreg); \
	vec2f64 *de1 = (vec2f64*)(dstreg+e1); \
	vec2f64 *de2 = (vec2f64*)(dstreg+e2); \
	vec2f64 *de3 = (vec2f64*)(dstreg+e2+e1); \
	vec2f64 *sc = (vec2f64*)(screg); \
	vec2f64 sc0 = sc[0], sc1 = sc[1], sc2 = sc[2], sc3 = sc[3]; \
	vec2f64 sc4 = sc[4], sc5 = sc[5], sc6 = sc[6], sc7 = sc[7]; \
	xprefetch(u8ptr(srcreg+srcinc)); \
	r4_x4c_fft4_mem8(sd0[0],sd1[0],sd2[0],sd3[0],sd0[1],sd1[1],sd2[1],sd3[1], de0[0],de0[1], sc0,sc1,sc2,sc3,sc4,sc5,sc6,sc7, dstreg+dstinc,e1); \
	xprefetch(u8ptr(srcreg+srcinc+d1)); \
	shuffle_store(de0[0], de0[1], xmm2, xmm7);					/* Save R1,R2 */ \
	shuffle_store_with_temp(de0[2], de0[3], xmm5, xmm6, xmm7);		/* Save I1,I2 */ \
	shuffle_store_with_temp(de1[0], de1[1], xmm3, xmm0, xmm7);		/* Save R3,R4 */ \
	shuffle_store_with_temp(de1[2], de1[3], xmm1, xmm4, xmm7);	/* Save I3,I4 */ \
	xprefetch(u8ptr(srcreg+srcinc+d2)); \
	r4_x4c_fft4_mem8(sd0[2],sd1[2],sd2[2],sd3[2],sd0[3],sd1[3],sd2[3],sd3[3], de2[0],de2[1], sc0,sc1,sc2,sc3,sc4,sc5,sc6,sc7, dstreg+dstinc+e2,e1); \
	xprefetch(u8ptr(srcreg+srcinc+d2+d1)); \
	srcreg += srcinc; \
	shuffle_store(de2[0], de2[1], xmm2, xmm7);				/* Save R1,R2 */ \
	shuffle_store_with_temp(de2[2], de2[3], xmm5, xmm6, xmm7);	/* Save I1,I2 */ \
	shuffle_store_with_temp(de3[0], de3[1], xmm3, xmm0, xmm7);	/* Save R3,R4 */ \
	shuffle_store_with_temp(de3[2], de3[3], xmm1, xmm4, xmm7);	/* Save I3,I4 */ \
	dstreg += dstinc; \
}
#endif

#define r4_x4c_fft4_mem8(mem1,mem2,mem3,mem4,mem5,mem6,mem7,mem8,unused1,unused2,SC0,SC1,SC2,SC3,SC4,SC5,SC6,SC7,pre1,pre2) { \
	vec2f64 m, r1, i1; \
	xmm2 = mem3;			/* R3 */ \
	xmm0 = mem1;			/* R1 */ \
	xmm5 = xmm2;			/* Copy R3 */ \
	xmm2 += xmm0;			/* R3 = R1 + R3 (new R1) */ \
	xmm0 -= xmm5;			/* R1 = R1 - R3 (new R3) */ \
\
	xmm3 = mem4;			/* R4 */ \
	xmm1 = mem2;			/* R2 */ \
	xmm5 = xmm3;			/* Copy R4 */ \
	xmm3 += xmm1;			/* R4 = R2 + R4 (new R2) */ \
	xmm1 -= xmm5;			/* R2 = R2 - R4 (new R4) */ \
\
	xmm6 = xmm2; \
	xmm2 -= xmm3;			/* R1 = R1 - R2 (final R2) */ \
	xmm3 += xmm6;			/* R2 = R1 + R2 (final R1) */ \
\
	xmm6 = mem7;			/* I3 */ \
	xmm4 = mem5;			/* I1 */ \
	xmm5 = xmm6;			/* Copy I3 */ \
	xmm6 += xmm4;			/* I3 = I1 + I3 (new I1) */ \
	xmm4 -= xmm5;			/* I1 = I1 - I3 (new I3) */ \
\
	xmm5 = mem6;			/* I2 */ \
	xmm7 = mem8;			/* I4 */ \
	xmm7 += xmm5;			/* I4 = I2 + I4 (new I2) */ \
	xmm5 -= mem8;			/* I2 = I2 - I4 (new I4) */ \
\
	r1 = xmm3;				/* Save R1 */ \
\
	xmm3 = xmm4; \
	xmm4 -= xmm1;			/* I3 = I3 - R4 (final I4) */ \
	xmm1 += xmm3;			/* R4 = I3 + R4 (final I3) */ \
\
	xmm3 = xmm0; \
	xmm0 -= xmm5;			/* R3 = R3 - I4 (final R3) */ \
	xmm5 += xmm3;			/* I4 = R3 + I4 (final R4) */ \
\
	xprefetchw(u8ptr(pre1)); \
\
	xmm3 = xmm6; \
	xmm6 -= xmm7;			/* I1 = I1 - I2 (final I2) */ \
	xmm7 += xmm3;			/* I2 = I1 + I2 (final I1) */ \
\
	m = SC3;				/* cosine/sine */ \
	xmm3 = m * xmm0 - xmm1;	/* A3 = R3 * cosine/sine - I3 */ \
	xmm1 = xmm1 * m + xmm0;	/* B3 = I3 * cosine/sine + R3 */ \
\
	i1 = xmm7;				/* Save I1 */ \
\
	m = SC5;				/* cosine/sine */ \
	xmm7 = m * xmm2 - xmm6;	/* A2 = R2 * cosine/sine - I2 */ \
	xmm6 = xmm6 * m + xmm2;	/* B2 = I2 * cosine/sine + R2 */ \
\
	xprefetchw(u8ptr(pre1+pre2)); \
\
	m = SC7;				/* cosine/sine */ \
	xmm0 = m * xmm5 - xmm4;	/* A4 = R4 * cosine/sine - I4 */ \
	xmm4 = xmm4 * m + xmm5;	/* B4 = I4 * cosine/sine + R4 */ \
\
	xmm2 = SC2; \
	xmm3 *= xmm2;			/* A3 = A3 * sine (final R3) */ \
	xmm1 *= xmm2;			/* B3 = B3 * sine (final I3) */ \
	xmm5 = SC4; \
	xmm7 *= xmm5;			/* A2 = A2 * sine (final R2) */ \
	xmm6 *= xmm5;			/* B2 = B2 * sine (final I2) */ \
\
	m = SC1;				/* cosine/sine */ \
	xmm2 = m * r1 - i1;		/* A1 = R1 * cosine/sine - I1 */ \
	xmm5 = i1 * m + r1;		/* B1 = I1 * cosine/sine + R1 */ \
\
	m = SC6; \
	xmm0 *= m;				/* A4 = A4 * sine (final R4) */ \
	xmm4 *= m;				/* B4 = B4 * sine (final I4) */ \
	m = SC0; \
	xmm2 *= m;				/* A1 = A1 * sine (final R1) */ \
	xmm5 *= m;				/* B1 = B1 * sine (final I1) */ \
}

/* 32-bit AMD K8 optimized versions of the above macros.  Derived from Intel 32-bit */ \
/* version -- could probably be optimized further. */ \

#if 0 // (@INSTR(,%xarch,<K8>) NE 0);
\

#define r4_x4c_fft4_mem(mem1,mem2,mem3,mem4,mem5,mem6,mem7,mem8,dest1,dest2,screg,off,pre1,pre2) \
	xmm2 = mem3;			/* R3 */ \
	xmm0 = mem1;			/* R1 */ \
	xmm2 += xmm0;			/* R3 = R1 + R3 (new R1) */ \
	xmm0 -= mem3;			/* R1 = R1 - R3 (new R3) */ \
\
	xmm3 = mem4;			/* R4 */ \
	xmm1 = mem2;			/* R2 */ \
	xmm3 += xmm1;			/* R4 = R2 + R4 (new R2) */ \
	xmm1 -= mem4;			/* R2 = R2 - R4 (new R4) */ \
\
	xmm6 = mem7;			/* I3 */ \
	xmm4 = mem5;			/* I1 */ \
	xmm6 += xmm4;			/* I3 = I1 + I3 (new I1) */ \
	xmm4 -= mem7;			/* I1 = I1 - I3 (new I3) */ \
\
	xmm5 = mem6;			/* I2 */ \
	xmm7 = mem8;			/* I4 */ \
	xmm7 += xmm5;			/* I4 = I2 + I4 (new I2) */ \
	xmm5 -= mem8;			/* I2 = I2 - I4 (new I4) */ \
\
	xprefetchw(u8ptr(pre1)); \
\
	xmm2 -= xmm3;			/* R1 = R1 - R2 (final R2) */ \
	xmm3 += xmm3; /* multwo */;				/* R2 = R2 * 2 */ \
\
	xmm4 -= xmm1;			/* I3 = I3 - R4 (final I4) */ \
	xmm1 += xmm1; /* multwo */;				/* R4 = R4 * 2 */ \
\
	xmm6 -= xmm7;			/* I1 = I1 - I2 (final I2) */ \
	xmm7 += xmm7; /* multwo */;				/* I2 = I2 * 2 */ \
\
	xmm0 -= xmm5;			/* R3 = R3 - I4 (final R3) */ \
	xmm5 += xmm5; /* multwo */;				/* I4 = I4 * 2 */ \
\
	xprefetchw(u8ptr(pre1+pre2)); \
\
	xmm3 += xmm2;			/* R2 = R1 + R2 (final R1) */ \
	xmm1 += xmm4;			/* R4 = I3 + R4 (final I3) */ \
	xmm7 += xmm6;			/* I2 = I1 + I2 (final I1) */ \
	xmm5 += xmm0;			/* I4 = R3 + I4 (final R4) */ \
\
	dest1 = xmm3;			/* Save R1 */ \
\
	xmm3 = xptr(screg+off+32+16);		/* cosine/sine */ \
	xmm3 *= xmm0;			/* A3 = R3 * cosine/sine */ \
\
	dest2 = xmm7;			/* Save I1 */ \
\
	xmm7 = xptr(screg+off+64+16);		/* cosine/sine */ \
	xmm7 *= xmm2;			/* A2 = R2 * cosine/sine */ \
\
	xmm3 -= xmm1;			/* A3 = A3 - I3 */ \
	xmm1 *= xptr(screg+off+32+16);		/* B3 = I3 * cosine/sine */ \
\
	xmm7 -= xmm6;			/* A2 = A2 - I2 */ \
	xmm6 *= xptr(screg+off+64+16);		/* B2 = I2 * cosine/sine */ \
\
	xmm1 += xmm0;			/* B3 = B3 + R3 */ \
\
	xmm0 = xptr(screg+off+96+16);		/* cosine/sine */ \
	xmm0 *= xmm5;			/* A4 = R4 * cosine/sine */ \
\
	xmm6 += xmm2;			/* B2 = B2 + R2 */ \
\
	xmm2 = xptr(screg+off+0+16);		/* cosine/sine */ \
	xmm2 *= dest1;			/* A1 = R1 * cosine/sine */ \
\
	xmm0 -= xmm4;			/* A4 = A4 - I4 */ \
	xmm4 *= xptr(screg+off+96+16);		/* B4 = I4 * cosine/sine */ \
\
	xmm4 += xmm5;			/* B4 = B4 + R4 */ \
\
	xmm5 = dest2;			/* Restore I1 */ \
	xmm2 -= xmm5;			/* A1 = A1 - I1 */ \
	xmm5 *= xptr(screg+off+0+16);		/* B1 = I1 * cosine/sine */ \
\
	xmm5 += dest1;			/* B1 = B1 + R1 */ \
\
	xmm3 *= xptr(screg+off+32);		/* A3 = A3 * sine (final R3) */ \
	xmm1 *= xptr(screg+off+32);		/* B3 = B3 * sine (final I3) */ \
	xmm7 *= xptr(screg+off+64);		/* A2 = A2 * sine (final R2) */ \
	xmm6 *= xptr(screg+off+64);		/* B2 = B2 * sine (final I2) */ \
	xmm0 *= xptr(screg+off+96);		/* A4 = A4 * sine (final R4) */ \
	xmm4 *= xptr(screg+off+96);		/* B4 = B4 * sine (final I4) */ \
	xmm2 *= xptr(screg+off+0);		/* A1 = A1 * sine (final R1) */ \
	xmm5 *= xptr(screg+off+0);		/* B1 = B1 * sine (final I1) */ \
\

#endif
\
/* 64-bit versions of the above macros. */ \

#if 0 //DEF X86_64;
\

#define r4_sg4cl_four_complex_fft4(srcreg,srcinc,d1,d2,dstreg,dstinc,e1,e2,screg) \
	xmm0 = xptr(srcreg);			/* R1 */ \
	xmm1 = xptr(srcreg+d2);		/* R3 */ \
	xmm2 = xmm0;			/* Copy R1			; 1-3 */ \
	xmm0 += xmm1;			/* R3 = R1 + R3 (new R1)	; 1-3 */ \
	xmm3 = xptr(srcreg+d1);		/* R2 */ \
	xmm4 = xptr(srcreg+d2+d1);		/* R4 */ \
	xmm5 = xmm3;			/* Copy R2			; 2-4 */ \
	xmm3 += xmm4;			/* R4 = R2 + R4 (new R2)	; 2-4 */ \
	xmm6 = xptr(srcreg+16);		/* I1 */ \
	xmm7 = xptr(srcreg+d2+16);		/* I3 */ \
	xmm8 = xmm6;			/* Copy I1			; 3-5 */ \
	xmm6 += xmm7;			/* I3 = I1 + I3 (new I1)	; 3-5 */ \
	xmm9 = xptr(srcreg+d1+16);		/* I2 */ \
	xmm10 = xptr(srcreg+d2+d1+16);	/* I4 */ \
	xmm11 = xmm9;			/* Copy I2			; 4-6 */ \
	xmm9 += xmm10;			/* I4 = I2 + I4 (new I2)	; 4-6	avail 12+ */ \
\
	xprefetch(u8ptr(srcreg+srcinc)); \
\
	xmm2 -= xmm1;			/* R1 = R1 - R3 (new R3)	; 5-7	avail 1,12+ */ \
\
	xmm5 -= xmm4;			/* R2 = R2 - R4 (new R4)	; 6-8	avail 1,4,12+ */ \
\
	xprefetchw(u8ptr(dstreg+dstinc)); \
\
	xmm8 -= xmm7;			/* I1 = I1 - I3 (new I3)	; 7-9	avail 1,4,7,12+ */ \
	xmm1 = xmm0;			/* Copy R1			; 7-9 (4) avail 4,7,12+ */ \
\
	xmm11 -= xmm10;			/* I2 = I2 - I4 (new I4)	; 8-10	avail 4,7,10,12+ */ \
	xmm12 = xptr(screg+64+16);		/* cosine/sine 2		; 8	avail 4,7,10,13+ */ \
\
	xmm0 -= xmm3;			/* R1 = R1 - R2 (final R2)	; 9-11 */ \
	xmm4 = xmm6;			/* Copy I1			; 9-11 (6) avail 7,10,13+ */ \
	xmm13 = xptr(screg+0+16);		/* cosine/sine 1		; 9 */ \
\
	xmm3 += xmm1;			/* R2 = R1 + R2 (final R1)	; 10-12 avail 7,10,1,14+ */ \
	xmm7 = xmm8;			/* Copy I3			; 10-12 avail 10,1,14+ */ \
\
	xprefetch(u8ptr(srcreg+srcinc+d1)); \
\
	xmm6 -= xmm9;			/* I1 = I1 - I2 (final I2)	; 11-13 */ \
	xmm1 = xmm2;			/* Copy R3			; 11-13 (8) avail 10,14+ */ \
\
	xmm9 += xmm4;			/* I2 = I1 + I2 (final I1)	; 12-14 avail 10,4,14+ */ \
	xmm10 = xmm0;			/* Copy R2			; 12-14 avail 4,14+ */ \
	xmm0 *= xmm12;			/* A2 = R2 * cosine/sine	; 12-16 */ \
\
	xprefetchw(u8ptr(dstreg+dstinc+e1)); \
\
	xmm8 -= xmm5;			/* I3 = I3 - R4 (final I4)	; 13-15 */ \
	xmm4 = xmm3;			/* Copy R1			; 13-15 avail 14+ */ \
	xmm3 *= xmm13;			/* A1 = R1 * cosine/sine	; 13-17 */ \
\
	xmm5 += xmm7;			/* R4 = I3 + R4 (final I3)	; 14-16 avail 7,14+ */ \
	xmm12 *= xmm6;			/* B2 = I2 * cosine/sine	; 14-18 */ \
	xmm14 = xptr(screg+96+16);		/* cosine/sine 4		; 14	avail 7,15 */ \
\
	xmm2 -= xmm11;			/* R3 = R3 - I4 (final R3)	; 15-17 */ \
	xmm13 *= xmm9;			/* B1 = I1 * cosine/sine	; 15-19 */ \
	xmm15 = xptr(screg+32+16);		/* cosine/sine 3		; 15	avail 7 */ \
\
	xmm11 += xmm1;			/* I4 = R3 + I4 (final R4)	; 16-18 avail 7,1 */ \
	xmm7 = xmm8;			/* Copy I4			; 16-18 avail 1 */ \
	xmm8 *= xmm14;			/* B4 = I4 * cosine/sine	; 16-20 */ \
\
	xprefetch(u8ptr(srcreg+srcinc+d2)); \
\
	xmm0 -= xmm6;			/* A2 = A2 - I2			; 17-19	avail 1,6 */ \
	xmm1 = xmm5;			/* Copy I3			; 17-19	avail 6 */ \
	xmm5 *= xmm15;			/* B3 = I3 * cosine/sine	; 17-21 */ \
\
	xmm3 -= xmm9;			/* A1 = A1 - I1			; 18-20	avail 6,9 */ \
	xmm15 *= xmm2;			/* A3 = R3 * cosine/sine	; 18-22 */ \
	xmm6 = xptr(screg+64);		/* sine 2			; 18	avail 9 */ \
\
	xmm12 += xmm10;			/* B2 = B2 + R2			; 19-21	avail 9,10 */ \
	xmm14 *= xmm11;			/* A4 = R4 * cosine/sine	; 19-23 */ \
	xmm9 = xptr(screg);			/* sine 1			; 19	avail 10 */ \
\
	xmm13 += xmm4;			/* B1 = B1 + R1			; 20-22	avail 10,4 */ \
	xmm0 *= xmm6;			/* A2 = A2 * sine (final R2)	; 20-24 */ \
	xmm4 = xptr(screg+96);		/* sine 4			; 20	avail 10,11 */ \
\
	xmm8 += xmm11;			/* B4 = B4 + R4			; 21-23	avail 10,11 */ \
	xmm3 *= xmm9;			/* A1 = A1 * sine (final R1)	; 21-25 */ \
	xmm11 = xptr(screg+32);		/* sine 3			; 21	avail 10 */ \
\
	xmm5 += xmm2;			/* B3 = B3 + R3			; 22-24	avail 10,2 */ \
	xmm12 *= xmm6;			/* B2 = B2 * sine (final I2)	; 22-26	avail 10,2,6 */ \
	xmm10 = xptr(srcreg+32);		/*#2 R1				; 22	avail 2,6 */ \
\
	xmm15 -= xmm1;			/* A3 = A3 - I3			; 23-25	avail 2,6,1 */ \
	xmm13 *= xmm9;			/* B1 = B1 * sine (final I1)	; 23-27	avail 2,6,1,9 */ \
	xmm6 = xmm10;			/*#2 Copy R1			; 1-3	avail 2,1,9 */ \
	xmm2 = xptr(srcreg+d2+32);		/*#2 R3				; 23	avail 1,9 */ \
\
	xmm14 -= xmm7;			/* A4 = A4 - I4			; 24-26	avail 1,9,7 */ \
	xmm8 *= xmm4;			/* B4 = B4 * sine (final I4)	; 24-28 */ \
	xmm1 = xptr(srcreg+d1+32);		/*#2 R2				; 24	avail 9,7 */ \
	xmm9 = xmm1;			/*#2 Copy R2			; 2-4	avail 7 */ \
\
	xmm10 += xmm2;			/*#2 R3 = R1 + R3 (new R1)	; 1-3 */ \
	xmm5 *= xmm11;			/* B3 = B3 * sine (final I3)	; 25-29 */ \
	xmm7 = xptr(srcreg+d2+d1+32);		/*#2 R4				; 25	avail none */ \
\
	xmm1 += xmm7;			/*#2 R4 = R2 + R4 (new R2)	; 2-4 */ \
	xmm15 *= xmm11;			/* A3 = A3 * sine (final R3)	; 26-30	avail 11 */ \
\
	shuffle_store_with_temp(xptr(dstreg), xptr(dstreg+16), xmm3, xmm0, xmm11);	/* Save R1,R2 in clocks 26-30  avail 3,0,11 */ \
\
	xmm0 = xptr(srcreg+48);		/*#2 I1 */ \
	xmm3 = xptr(srcreg+d2+48);		/*#2 I3 */ \
	xmm11 = xmm0;			/*#2 Copy I1			; 3-5 */ \
	xmm0 += xmm3;			/*#2 I3 = I1 + I3 (new I1)	; 3-5 */ \
	xmm14 *= xmm4;			/* A4 = A4 * sine (final R4)	; 27-31	avail 4 */ \
\
	shuffle_store_with_temp(xptr(dstreg+32), xptr(dstreg+48), xmm13, xmm12, xmm4);	/* Save I1,I2 in clocks 28-32  avail 13,12,4 */ \
\
	xmm4 = xptr(srcreg+d1+48);		/*#2 I2 */ \
	xmm12 = xptr(srcreg+d2+d1+48);	/*#2 I4 */ \
	xmm13 = xmm4;			/*#2 Copy I2			; 4-6 */ \
	xmm4 += xmm12;			/*#2 I4 = I2 + I4 (new I2)	; 4-6	avail none */ \
\
	xprefetchw(u8ptr(dstreg+dstinc+e2)); \
\
	xmm6 -= xmm2;			/*#2 R1 = R1 - R3 (new R3)	; 5-7	avail 2 */ \
\
	xmm9 -= xmm7;			/*#2 R2 = R2 - R4 (new R4)	; 6-8	avail 2,7 */ \
\
	shuffle_store_with_temp(xptr(dstreg+e1+32), xptr(dstreg+e1+48), xmm5, xmm8, xmm2); /* Save I3,I4 in clocks 30-34 */ \
\
	xmm11 -= xmm3;			/*#2 I1 = I1 - I3 (new I3)	; 7-9 */ \
	xmm2 = xmm10;			/*#2 Copy R1			; 7-9 (4) */ \
\
	xprefetch(u8ptr(srcreg+srcinc+d2+d1)); \
\
	xmm13 -= xmm12;			/*#2 I2 = I2 - I4 (new I4)	; 8-10 */ \
	xmm5 = xptr(screg+64+16);		/*#2 cosine/sine 2		; 8 */ \
\
	shuffle_store_with_temp(xptr(dstreg+e1), xptr(dstreg+e1+16), xmm15, xmm14, xmm3);	/* Save R3,R4 in clocks 32-36 */ \
\
	xmm10 -= xmm1;			/*#2 R1 = R1 - R2 (final R2)	; 9-11 */ \
	xmm7 = xmm0;			/*#2 Copy I1			; 9-11 (6) */ \
	xmm8 = xptr(screg+0+16);		/*#2 cosine/sine 1		; 9 */ \
\
	xmm1 += xmm2;			/*#2 R2 = R1 + R2 (final R1)	; 10-12 */ \
	xmm3 = xmm11;			/*#2 Copy I3			; 10-12 */ \
\
	xprefetchw(u8ptr(dstreg+dstinc+e2+e1)); \
\
	xmm0 -= xmm4;			/*#2 I1 = I1 - I2 (final I2)	; 11-13 */ \
	xmm2 = xmm6;			/*#2 Copy R3			; 11-13 (8) */ \
\
	xmm4 += xmm7;			/*#2 I2 = I1 + I2 (final I1)	; 12-14 */ \
	xmm12 = xmm10;			/*#2 Copy R2			; 12-14 */ \
	xmm10 *= xmm5;			/*#2 A2 = R2 * cosine/sine	; 12-16 */ \
\
	xmm11 -= xmm9;			/*#2 I3 = I3 - R4 (final I4)	; 13-15 */ \
	xmm7 = xmm1;			/*#2 Copy R1			; 13-15 */ \
	xmm1 *= xmm8;			/*#2 A1 = R1 * cosine/sine	; 13-17 */ \
\
	xmm9 += xmm3;			/*#2 R4 = I3 + R4 (final I3)	; 14-16 */ \
	xmm5 *= xmm0;			/*#2 B2 = I2 * cosine/sine	; 14-18 */ \
	xmm14 = xptr(screg+96+16);		/*#2 cosine/sine 4		; 14 */ \
\
	xmm6 -= xmm13;			/*#2 R3 = R3 - I4 (final R3)	; 15-17 */ \
	xmm8 *= xmm4;			/*#2 B1 = I1 * cosine/sine	; 15-19 */ \
	xmm15 = xptr(screg+32+16);		/*#2 cosine/sine 3		; 15 */ \
\
	xmm13 += xmm2;			/*#2 I4 = R3 + I4 (final R4)	; 16-18 */ \
	xmm3 = xmm11;			/*#2 Copy I4			; 16-18 */ \
	xmm11 *= xmm14;			/*#2 B4 = I4 * cosine/sine	; 16-20 */ \
\
	xmm10 -= xmm0;			/*#2 A2 = A2 - I2		; 17-19 */ \
	xmm2 = xmm9;			/*#2 Copy I3			; 17-19 */ \
	xmm9 *= xmm15;			/*#2 B3 = I3 * cosine/sine	; 17-21 */ \
\
	xmm1 -= xmm4;			/*#2 A1 = A1 - I1		; 18-20 */ \
	xmm15 *= xmm6;			/*#2 A3 = R3 * cosine/sine	; 18-22 */ \
	xmm0 = xptr(screg+64);		/*#2 sine 2			; 18 */ \
\
	xmm5 += xmm12;			/*#2 B2 = B2 + R2		; 19-21 */ \
	xmm14 *= xmm13;			/*#2 A4 = R4 * cosine/sine	; 19-23 */ \
	xmm4 = xptr(screg);			/*#2 sine 1			; 19 */ \
\
	xmm8 += xmm7;			/*#2 B1 = B1 + R1		; 20-22 */ \
	xmm10 *= xmm0;			/*#2 A2 = A2 * sine (final R2)	; 20-24 */ \
\
	srcreg += srcinc; \
\
	xmm11 += xmm13;			/*#2 B4 = B4 + R4		; 21-23 */ \
	xmm1 *= xmm4;			/*#2 A1 = A1 * sine (final R1)	; 21-25 */ \
	xmm7 = xptr(screg+96);		/*#2 sine 4			; 21 */ \
\
	xmm9 += xmm6;			/*#2 B3 = B3 + R3		; 22-24 */ \
	xmm5 *= xmm0;			/*#2 B2 = B2 * sine (final I2)	; 22-26 */ \
	xmm6 = xptr(screg+32);		/*#2 sine 3			; 22 */ \
\
	xmm15 -= xmm2;			/*#2 A3 = A3 - I3		; 23-25 */ \
	xmm8 *= xmm4;			/*#2 B1 = B1 * sine (final I1)	; 23-27 */ \
\
	xmm14 -= xmm3;			/*#2 A4 = A4 - I4		; 24-26 */ \
	xmm11 *= xmm7;			/*#2 B4 = B4 * sine (final I4)	; 24-28 */ \
\
	xmm9 *= xmm6;			/*#2 B3 = B3 * sine (final I3)	; 25-29 */ \
\
	xmm15 *= xmm6;			/*#2 A3 = A3 * sine (final R3)	; 26-30 */ \
\
	xmm14 *= xmm7;			/*#2 A4 = A4 * sine (final R4)	; 27-31 */ \
\
	shuffle_store_with_temp(xptr(dstreg+e2), xptr(dstreg+e2+16), xmm1, xmm10, xmm2);		/*#2 Save R1,R2 */ \
	shuffle_store_with_temp(xptr(dstreg+e2+32), xptr(dstreg+e2+48), xmm8, xmm5, xmm2);	/*#2 Save I1,I2 */ \
	shuffle_store_with_temp(xptr(dstreg+e2+e1), xptr(dstreg+e2+e1+16), xmm15, xmm14, xmm2);	/*#2 Save R3,R4 */ \
	shuffle_store_with_temp(xptr(dstreg+e2+e1+32), xptr(dstreg+e2+e1+48), xmm9, xmm11, xmm2);	/*#2 Save I3,I4 */ \
	dstreg += dstinc; \
\

#endif
\
/* */ \
/* ************************************* four-complex-unfft4 variants ****************************************** */ \
/* */ \
\
/* Used in last levels of pass 1 (split premultiplier and delay cases).  No swizzling. */
#if 0 //IFDEF UNUSED

#define r4_g2cl_four_complex_unfft4(srcreg,srcinc,d1,dstreg,dstinc,e1,screg) \
	xmm0 = xptr(srcreg);		/* R1 */ \
	xmm1 = xptr(srcreg+16);	/* I1 */ \
	xmm4 = xptr(srcreg+d1);	/* R3 */ \
	xmm5 = xptr(srcreg+d1+16);	/* I3 */ \
	r4_x4c_unfft4_partial_mem(xmm0,xmm1,xmm2,xmm3,xmm4,xmm5,xmm6,xmm7,xptr(srcreg+32),xptr(srcreg+48),xptr(srcreg+d1+32),xptr(srcreg+d1+48),xptr(dstreg),xptr(dstreg+32),screg,0,dstreg+dstinc,e1); \
	srcreg += srcinc; \
	xptr(dstreg) = xmm7;		/* Save R1 */ \
	xptr(dstreg+16) = xmm2;	/* Save R3 */ \
	xptr(dstreg+32) = xmm4;	/* Save I1 */ \
	xptr(dstreg+48) = xmm3;	/* Save I3 */ \
	xptr(dstreg+e1) = xmm5;	/* Save R2 */ \
	xptr(dstreg+e1+16) = xmm6;	/* Save R4 */ \
	xptr(dstreg+e1+32) = xmm0;	/* Save I2 */ \
	xptr(dstreg+e1+48) = xmm1;	/* Save I4 */ \
	dstreg += dstinc; \
\
\

#define r4_x4c_unfft4_partial_mem(r1,r2,r3,r4,r5,r6,r7,r8,mem3,mem4,mem7,mem8,dest1,dest2,screg,off,pre1,pre2) \
	r7 = xptr(screg+off+0+16);		/* cosine/sine */ \
	r7 *= r1;				/* A1 = R1 * cosine/sine */ \
	r7 += r2;				/* A1 = A1 + I1 */ \
	r2 *= xptr(screg+off+0+16);		/* B1 = I1 * cosine/sine */ \
	r2 -= r1;				/* B1 = B1 - R1 */ \
\
	r8 = xptr(screg+off+32+16);		/* cosine/sine */ \
	r8 *= r5;				/* A3 = R3 * cosine/sine */ \
	r8 += r6;				/* A3 = A3 + I3 */ \
	r6 *= xptr(screg+off+32+16);		/* B3 = I3 * cosine/sine */ \
	r6 -= r5;				/* B3 = B3 - R3 */ \
\
	r3 = xptr(screg+off+64+16);		/* cosine/sine */ \
	r3 *= mem3;			/* A2 = R2 * cosine/sine */ \
	r4 = mem4;			/* R4 (I2) */ \
	r3 += r4;				/* A2 = A2 + I2 */ \
	r4 *= xptr(screg+off+64+16);		/* B2 = I2 * cosine/sine */ \
	r4 -= mem3;			/* B2 = B2 - R2 */ \
\
	r1 = xptr(screg+off+96+16);		/* cosine/sine */ \
	r1 *= mem7;			/* A4 = R4 * cosine/sine */ \
	r5 = mem8;			/* I4 */ \
	r1 += r5;				/* A4 = A4 + I4 */ \
	r5 *= xptr(screg+off+96+16);		/* B4 = I4 * cosine/sine */ \
	r5 -= mem7;			/* B4 = B4 - R4 */ \
\
	r7 *= xptr(screg+off+0);		/* A1 = A1 * sine (new R1) */ \
	r2 *= xptr(screg+off+0);		/* B1 = B1 * sine (new I1) */ \
\
	r8 *= xptr(screg+off+32);		/* A3 = A3 * sine (new R3) */ \
	r6 *= xptr(screg+off+32);		/* B3 = B3 * sine (new I3) */ \
\
	r3 *= xptr(screg+off+64);		/* A2 = A2 * sine (new R2) */ \
	r4 *= xptr(screg+off+64);		/* B2 = B2 * sine (new I2) */ \
\
	r1 *= xptr(screg+off+96);		/* A4 = A4 * sine (new R4) */ \
	r5 *= xptr(screg+off+96);		/* B4 = B4 * sine (new I4) */ \
\
	xprefetchw(u8ptr(pre1)); \
\
	r7 -= r3;				/* R1 = R1 - R2 (new R2) */ \
	r3 += r3; /* multwo */				/* R2 = R2 * 2 */ \
	r3 += r7;				/* R2 = R1 + R2 (new R1) */ \
\
	r2 -= r4;				/* I1 = I1 - I2 (new I2) */ \
	r4 += r4; /* multwo */				/* I2 = I2 * 2 */ \
	r4 += r2;				/* I2 = I1 + I2 (new I1) */ \
\
	r1 -= r8;				/* R4 = R4 - R3 (new I4) */ \
	r8 += r8; /* multwo */				/* R3 = R3 * 2 */ \
	r8 += r1;				/* R3 = R4 + R3 (new R3) */ \
\
	r6 -= r5;				/* I3 = I3 - I4 (new R4) */ \
	r5 += r5; /* multwo */				/* I4 = I4 * 2 */ \
	r5 += r6;				/* I4 = I3 + I4 (new I3) */ \
\
	xprefetchw(u8ptr(pre1+pre2)); \
\
	r7 -= r6;				/* R2 = R2 - R4 (final R4) */ \
	r6 += r6; /* multwo */				/* R4 = R4 * 2 */ \
	r6 += r7;				/* R4 = R2 + R4 (final R2) */ \
\
	r2 -= r1;				/* I2 = I2 - I4 (final I4) */ \
	r1 += r1; /* multwo */				/* I4 = I4 * 2 */ \
	r1 += r2;				/* I4 = I2 + I4 (final I2) */ \
\
	r3 -= r8;				/* R1 = R1 - R3 (final R3) */ \
	r8 += r8; /* multwo */				/* R3 = R3 * 2 */ \
	r8 += r3;				/* R3 = R1 + R3 (final R1) */ \
\
	r4 -= r5;				/* I1 = I1 - I3 (final I3) */ \
	r5 += r5; /* multwo */				/* I3 = I3 * 2 */ \
	r5 += r4;				/* I3 = I1 + I3 (final I1) */ \

#endif
\
/* Used in last levels of pass 1 (r4delay case).  Swizzling. */ \

#define r4_sg2cl_four_complex_unfft4(srcreg,srcinc,d1,dstreg,dstinc,e1,screg) { \
	vec2f64 *sd0 = (vec2f64*)(srcreg); \
	vec2f64 *sd1 = (vec2f64*)(srcreg+d1); \
	vec2f64 *de0 = (vec2f64*)(dstreg); \
	vec2f64 *de1 = (vec2f64*)(dstreg+e1); \
	vec2f64 *sc = (vec2f64*)(screg); \
\
	shuffle_load_with_temp(xmm0, xmm4, sd0[0], sd0[1], xmm7);		/* R1,R2 */ \
	shuffle_load_with_temp(xmm1, xmm3, sd0[2], sd0[3], xmm7);	/* I1,I2 */ \
\
	xmm6 = xmm0 * sc[1] + xmm1;	/* A1 = R1 * cosine/sine + I1 */ \
	xmm2 = xmm4 * sc[5] + xmm3;	/* A2 = R2 * cosine/sine + I2 */ \
	xmm1 = xmm1 * sc[1] - xmm0;	/* B1 = I1 * cosine/sine - R1 */ \
	xmm3 = xmm3 * sc[5] - xmm4;	/* B2 = I2 * cosine/sine - R2 */ \
\
	xmm6 *= sc[0];			/* A1 = A1 * sine (new R1) */ \
	xmm2 *= sc[4];			/* A2 = A2 * sine (new R2) */ \
\
	vec2f64 saved_R1 = xmm6;	/* Save new R1 temporarily */ \
	vec2f64 saved_R2 = xmm2;	/* Save new R2 temporarily */ \
\
	shuffle_load_with_temp(xmm4, xmm5, sd1[0], sd1[1], xmm7);	/* R3,R4 */ \
	shuffle_load_with_temp(xmm6, xmm2, sd1[2], sd1[3], xmm7);	/* I3,I4 */ \
\
	xmm7 = xmm4 * sc[3] + xmm6;	/* A3 = R3 * cosine/sine + I3 */ \
	xmm0 = xmm5 * sc[7] + xmm2;	/* A4 = R4 * cosine/sine + I4 */ \
	xmm6 = xmm6 * sc[3] - xmm4;	/* B3 = I3 * cosine/sine - R3 */ \
	xmm2 = xmm2 * sc[7] - xmm5;	/* B4 = I4 * cosine/sine - R4 */ \
\
	xmm1 *= sc[0];			/* B1 = B1 * sine (new I1) */ \
	xmm3 *= sc[4];			/* B2 = B2 * sine (new I2) */ \
	xmm7 *= sc[2];			/* A3 = A3 * sine (new R3) */ \
	xmm0 *= sc[6];			/* A4 = A4 * sine (new R4) */ \
	xmm6 *= sc[2];			/* B3 = B3 * sine (new I3) */ \
	xmm2 *= sc[6];			/* B4 = B4 * sine (new I4) */ \
\
	srcreg += srcinc; \
\
	xmm4 = xmm1;			/* Copy I1 */ \
	xmm1 -= xmm3;			/* I1 = I1 - I2 (new I2) */ \
	xmm3 += xmm4;			/* I2 = I1 + I2 (new I1) */ \
\
	xprefetchw(u8ptr(dstreg+dstinc)); \
\
	xmm4 = xmm0;			/* Copy R4 */ \
	xmm0 -= xmm7;			/* R4 = R4 - R3 (new I4) */ \
	xmm7 += xmm4;			/* R3 = R4 + R3 (new R3) */ \
\
	xprefetchw(u8ptr(dstreg+dstinc+e1)); \
\
	xmm4 = xmm6;			/* Copy I3 */ \
	xmm6 -= xmm2;			/* I3 = I3 - I4 (new R4) */ \
	xmm2 += xmm4;			/* I4 = I3 + I4 (new I3) */ \
\
	xprefetch(u8ptr(srcreg)); \
\
	xmm4 = xmm1;			/* Copy I2 */ \
	xmm1 -= xmm0;			/* I2 = I2 - I4 (final I4) */ \
	xmm0 += xmm4;			/* I4 = I2 + I4 (final I2) */ \
\
	xprefetch(u8ptr(srcreg+d1)); \
\
	xmm4 = saved_R1;		/* Reload new R1 */ \
	xmm5 = saved_R2;		/* Reload new R2 */ \
	xmm4 -= xmm5;			/* R1 = R1 - R2 (new R2) */ \
	xmm5 += saved_R1;		/* R2 = R1 + R2 (new R1) */ \
\
	de1[3] = xmm1;			/* Save I4 */ \
	de1[2] = xmm0;			/* Save I2 */ \
\
	xmm1 = xmm3;			/* Copy I1 */ \
	xmm3 -= xmm2;			/* I1 = I1 - I3 (final I3) */ \
	xmm2 += xmm1;			/* I3 = I1 + I3 (final I1) */ \
\
	xmm1 = xmm4;			/* Copy R2 */ \
	xmm4 -= xmm6;			/* R2 = R2 - R4 (final R4) */ \
	xmm6 += xmm1;			/* R4 = R2 + R4 (final R2) */ \
\
	xmm1 = xmm5;			/* Copy R1 */ \
	xmm5 -= xmm7;			/* R1 = R1 - R3 (final R3) */ \
	xmm7 += xmm1;			/* R3 = R1 + R3 (final R1) */ \
\
	de0[2] = xmm2;			/* Save I1 */ \
	de0[3] = xmm3;			/* Save I3 */ \
	de1[0] = xmm6;			/* Save R2 */ \
	de1[1] = xmm4;			/* Save R4 */ \
	de0[0] = xmm7;			/* Save R1 */ \
	de0[1] = xmm5;			/* Save R3 */ \
	dstreg += dstinc; \
}

#if 0 //DEF X86_64;
\

#define r4_sg2cl_four_complex_unfft4(srcreg,srcinc,d1,dstreg,dstinc,e1,screg) \
	shuffle_load_with_temp(xmm1, xmm2, xptr(srcreg), xptr(srcreg+16), xmm0);		/* R1,R2 */ \
\
	xmm3 = xptr(screg+0+16);		/* cosine/sine 1		; 1 */ \
	xmm15 = xmm1;			/* Copy R1			; 1-3 */ \
	xmm1 *= xmm3;			/* A1 = R1 * cosine/sine	; 1-5 */ \
\
	xmm4 = xptr(screg+64+16);		/* cosine/sine 2		; 2 */ \
	xmm5 = xmm2;			/* Copy R2			; 2-4 */ \
	xmm2 *= xmm4;			/* A2 = R2 * cosine/sine	; 2-6 */ \
\
	shuffle_load_with_temp(xmm6, xmm7, xptr(srcreg+32), xptr(srcreg+48), xmm0);	/* I1,I2 */ \
\
	xmm3 *= xmm6;			/* B1 = I1 * cosine/sine	; 3-7 */ \
	xmm4 *= xmm7;			/* B2 = I2 * cosine/sine	; 4-8 */ \
\
	shuffle_load_with_temp(xmm8, xmm9, xptr(srcreg+d1), xptr(srcreg+d1+16), xmm0);	/* R3,R4 */ \
\
	xmm10 = xptr(screg+32+16);		/* cosine/sine 3		; 5 */ \
	xmm11 = xmm8;			/* Copy R3			; 5-7 */ \
	xmm8 *= xmm10;			/* A3 = R3 * cosine/sine	; 5-9 */ \
\
	xmm1 += xmm6;			/* A1 = A1 + I1			; 6-8	avail 0,6,12-14 */ \
	xmm12 = xptr(screg+96+16);		/* cosine/sine 4		; 6	avail 0,6,13,14 */ \
	xmm13 = xmm9;			/* Copy R4			; 6-8	avail 0,6,14 */ \
	xmm9 *= xmm12;			/* A4 = R4 * cosine/sine	; 6-10 */ \
\
	shuffle_load_with_temp(xmm6, xmm14, xptr(srcreg+d1+32), xptr(srcreg+d1+48), xmm0); /* I3,I4 */ \
\
	xmm2 += xmm7;			/* A2 = A2 + I2			; 7-9	avail 0,7 */ \
	xmm10 *= xmm6;			/* B3 = I3 * cosine/sine	; 7-11 */ \
	xmm0 = xptr(screg+0);			/* sine 1			; 7	avail 7 */ \
\
	xmm3 -= xmm15;			/* B1 = B1 - R1			; 8-10	avail 7,15 */ \
	xmm12 *= xmm14;			/* B4 = I4 * cosine/sine	; 8-12 */ \
	xmm7 = xptr(screg+64);		/* sine 2			; 8	avail 15 */ \
\
	xmm4 -= xmm5;			/* B2 = B2 - R2			; 9-11	avail 15,5 */ \
	xmm1 *= xmm0;			/* A1 = A1 * sine (new R1)	; 9-13 */ \
\
	xmm8 += xmm6;			/* A3 = A3 + I3			; 10-12	avail 15,5,6 */ \
	xmm2 *= xmm7;			/* A2 = A2 * sine (new R2)	; 10-14 */ \
	xmm5 = xptr(screg+32);		/* sine 3			; 10	avail 15,6 */ \
\
	xmm9 += xmm14;			/* A4 = A4 + I4			; 11-13	avail 15,6,14 */ \
	xmm3 *= xmm0;			/* B1 = B1 * sine (new I1)	; 11-15	avail 15,6,14,0 */ \
	xmm6 = xptr(screg+96);		/* sine 4			; 11	avail 15,14,0 */ \
\
	xmm10 -= xmm11;			/* B3 = B3 - R3			; 12-14	avail 15,14,0,11 */ \
	xmm4 *= xmm7;			/* B2 = B2 * sine (new I2)	; 12-16	avail 15,14,0,11,7 */ \
\
	xmm12 -= xmm13;			/* B4 = B4 - R4			; 13-15	avail 15,14,0,11,7,13 */ \
	xmm8 *= xmm5;			/* A3 = A3 * sine (new R3)	; 13-17 */ \
\
	xmm0 = xmm1;			/* Copy R1			; 14-16 (14 on a Core i7) avail 15,14,11,7,13 */ \
	xmm9 *= xmm6;			/* A4 = A4 * sine (new R4)	; 14-18 */ \
	srcreg += srcinc; \
\
	xmm1 -= xmm2;			/* R1 = R1 - R2 (newer R2)	; 15-17 */ \
	xmm10 *= xmm5;			/* B3 = B3 * sine (new I3)	; 15-19	avail 15,14,11,7,13,5 */ \
\
	xprefetchw(u8ptr(dstreg+dstinc)); \
\
	xmm2 += xmm0;			/* R2 = R1 + R2 (newer R1)	; 17-19 (16-18 on a Core i7) avail 15,14,11,7,13,5,0 */ \
	xmm12 *= xmm6;			/* B4 = B4 * sine (new I4)	; 16-20	avail 15,14,11,7,13,5,0,12 */ \
	xmm0 = xmm3;			/* Copy I1			; 16-18	avail 15,14,11,7,13,5,12 */ \
\
	xmm3 -= xmm4;			/* I1 = I1 - I2 (newer I2)	; 18-20 */ \
	xmm5 = xmm8;			/* Copy R3			; 18-20	avail 15,14,11,7,13,12 */ \
\
	xmm4 += xmm0;			/* I2 = I1 + I2 (newer I1)	; 19-21	avail 15,14,11,7,13,12,0 */ \
\
	xprefetchw(u8ptr(dstreg+dstinc+e1)); \
\
	xmm8 += xmm9;			/* R3 = R4 + R3 (newer R3)	; 20-22 */ \
	xmm0 = xmm10;			/* Copy I3			; 20-22	avail 15,14,11,7,13,12 */ \
\
	xmm9 -= xmm5;			/* R4 = R4 - R3 (newer I4)	; 21-23	avail 15,14,11,7,13,12,5 */ \
	xmm7 = xmm2;			/* Copy R1			; 21-23 (20) avail 15,14,11,13,12,5 */ \
\
	xmm10 -= xmm12;			/* I3 = I3 - I4 (newer R4)	; 22-24 */ \
	xmm5 = xmm3;			/* Copy I2			; 22-24 (21) avail 15,14,11,13,12 */ \
\
	xprefetch(u8ptr(srcreg)); \
\
	xmm12 += xmm0;			/* I4 = I3 + I4 (newer I3)	; 23-25	avail 15,14,11,13,12,0 */ \
	xmm11 = xmm1;			/* Copy R2			; 23-25 (18) avail 15,14,13,12,0 */ \
\
	xmm2 -= xmm8;			/* R1 = R1 - R3 (final R3)	; 24-26 */ \
	xmm0 = xmm4;			/* Copy I1			; 24-26 (22) avail 15,14,13,12 */ \
\
	xmm8 += xmm7;			/* R3 = R1 + R3 (final R1)	; 25-27	avail 15,14,13,12,7 */ \
\
	xprefetch(u8ptr(srcreg+d1)); \
\
	xmm3 -= xmm9;			/* I2 = I2 - I4 (final I4)	; 26-28 */ \
\
	xmm9 += xmm5;			/* I4 = I2 + I4 (final I2)	; 27-29	avail 15,14,13,12,7,5 */ \
	xptr(dstreg+16) = xmm2;		/* Save R3			; 27 */ \
\
	xmm1 -= xmm10;			/* R2 = R2 - R4 (final R4)	; 28-30 */ \
	xptr(dstreg) = xmm8;			/* Save R1			; 28 */ \
\
	xmm10 += xmm11;			/* R4 = R2 + R4 (final R2)	; 29-31	avail 15,14,13,12,7,5,11 */ \
	xptr(dstreg+e1+48) = xmm3;		/* Save I4			; 29 */ \
\
	xmm4 -= xmm12;			/* I1 = I1 - I3 (final I3)	; 30-32 */ \
	xptr(dstreg+e1+32) = xmm9;		/* Save I2			; 30 */ \
\
	xmm12 += xmm0;			/* I3 = I1 + I3 (final I1)	; 31-33	avail 15,14,13,12,7,5,11,0 */ \
	xptr(dstreg+e1+16) = xmm1;		/* Save R4			; 31 */ \
\
	xptr(dstreg+e1) = xmm10;		/* Save R2			; 32 */ \
	xptr(dstreg+48) = xmm4;		/* Save I3			; 33 */ \
	xptr(dstreg+32) = xmm12;		/* Save I1			; 34 */ \
	dstreg += dstinc; \

#endif
\
/* */ \
/* ************************************* eight-reals-fft variants ****************************************** */ \
/* */ \
\
/* These macros operate on eight reals doing 2 and 3/4 levels of the FFT and applying */ \
/* the sin/cos multipliers afterwards.  The output is 2 reals (only 2 levels of FFT done) */ \
/* and 3 complex numbers (3 levels of FFT performed).  These macros take a screg */ \
/* that points to twiddles w^n, w^2n, and w^5n. */ \
\
/*r4_x2cl_eight_reals_first_fft_preload MACRO */ \
/*	r4_x2cl_eight_reals_first_fft_cmn_preload */ \
/*	ENDM */ \
/*r4_x2cl_eight_reals_first_fft MACRO srcreg,srcinc,d1,screg */ \
/*	r4_x2cl_eight_reals_first_fft_cmn srcreg,rbx,srcinc,d1,screg */ \
/*	ENDM */ \
\
/* r4_x2cl_eight_reals_first_fft_cmn unused */

#define r4_x2cl_eight_reals_first_fft_scratch_preload \
	r4_x8r_fft_mem_preload; \

#define r4_x2cl_eight_reals_first_fft_scratch(srcreg,srcinc,d1,screg) { \
	vec2f64 *sd0 = (vec2f64*)(srcreg); \
	vec2f64 *sd1 = (vec2f64*)(srcreg+d1); \
	r4_x8r_fft_mem(sd0[0],sd1[0],sd0[1],sd1[1],sd0[2],sd1[2],sd0[3],sd1[3], screg,screg+32,screg+64,srcreg+srcinc,d1, sd0[0],sd0[1]); \
/*	xstore	[srcreg], xmm0		;; Save real value #1 */ \
/*	xstore	[srcreg+16], xmm0	;; Save real value #2 */ \
	sd0[2] = xmm7;	/* Save R2 */ \
	sd0[3] = xmm5;	/* Save I2 */ \
	sd1[0] = xmm4;	/* Save R3 */ \
	sd1[1] = xmm3;	/* Save I3 */ \
	sd1[2] = xmm6;	/* Save R4 */ \
	sd1[3] = xmm2;	/* Save I4 */ \
	srcreg += srcinc; \
}

/* Common macro for eight-reals-FFT doing 2 and 3/4 levels. */ \
\

#define r4_x8r_fft_mem_preload \
\

#define r4_x8r_fft_mem(mem1,mem2,mem3,mem4,mem5,mem6,mem7,mem8,screg1,screg2,screg5,pre1,pre2,dst1,dst2) { \
	vec2f64 t; \
	xmm3 = mem4;		/* R4 */ \
	xmm7 = mem8;		/* R8 */ \
	xmm0 = xmm3; \
	xmm3 -= xmm7;		/* new R8 = R4 - R8 */ \
	xmm7 += xmm0;		/* new R4 = R4 + R8 */ \
\
	xmm1 = mem2;		/* R2 */ \
	xmm5 = mem6;		/* R6 */ \
	xmm4 = xmm1; \
	xmm1 -= xmm5;		/* new R6 = R2 - R6 */ \
	xmm5 += xmm4;		/* new R2 = R2 + R6 */ \
\
	xmm3 *= XMM_SQRTHALF2;	/* R8 = R8 * square root */ \
	xmm1 *= XMM_SQRTHALF2;	/* R6 = R6 * square root */ \
\
	xmm0 = mem1;		/* R1 */ \
	xmm4 = mem5;		/* R5 */ \
	xmm2 = xmm0; \
	xmm0 -= xmm4;		/* new R5 = R1 - R5 */ \
	xmm4 += xmm2;		/* new R1 = R1 + R5 */ \
\
	xmm2 = xmm5;		/* Copy R2 */ \
	xmm5 -= xmm7;		/* R2 = R2 - R4 (final I2) */ \
	xmm7 += xmm2;		/* R4 = R2 + R4 (final I1, a.k.a 2nd real result) */ \
\
	xmm6 = xmm1;		/* Copy R6 */ \
	xmm1 -= xmm3;		/* R6 = R6 - R8 (Real part) */ \
	xmm3 += xmm6;		/* R8 = R6 + R8 (Imaginary part) */ \
\
	xmm6 = mem7;		/* R7 */ \
	xmm6 += mem3;		/* new R3 = R3 + R7 */ \
\
	xmm2 = xmm4;		/* Copy R1 */ \
	xmm4 -= xmm6;		/* R1 = R1 - R3 (final R2) */ \
	xmm6 += xmm2;		/* R3 = R1 + R3 (final R1) */ \
\
	xmm2 = mem3;		/* R3 */ \
	xmm2 -= mem7;		/* new R7 = R3 - R7 */ \
\
	xprefetchw(u8ptr(pre1)); \
\
	dst2 = xmm7;		/* Save I1 */ \
	dst1 = xmm6;		/* Save R1 */ \
\
	xmm7 = xmm0;		/* Copy R5 */ \
	xmm0 -= xmm1;		/* R5 = R5 - R6 (final R4) */ \
	xmm1 += xmm7;		/* R6 = R5 + R6 (final R3) */ \
\
	xmm6 = xmm2;		/* Copy R7 */ \
	xmm2 -= xmm3;		/* R7 = R7 - R8 (final I4) */ \
	xmm3 += xmm6;		/* R8 = R7 + R8 (final I3) */ \
\
	xprefetchw(u8ptr(pre1+pre2)); \
\
	t = xptr(screg2+16);	/* cosine/sine for w^2n */ \
	xmm7 = t * xmm4 - xmm5;	/* A2 = R2 * cosine/sine - I2 */ \
	xmm5 = t * xmm5 + xmm4;	/* B2 = I2 * cosine/sine + R2 */ \
\
	t = xptr(screg5+16);	/* cosine/sine for w^5n */ \
	xmm6 = t * xmm0 - xmm2;	/* A4 = R4 * cosine/sine - I4 */ \
	xmm2 = t * xmm2 + xmm0;	/* B4 = I4 * cosine/sine + R4 */ \
\
	t = xptr(screg1+16);	/* cosine/sine for w^n */ \
	xmm4 = t * xmm1 - xmm3;	/* A3 = R3 * cosine/sine - I3 */ \
	xmm3 = t * xmm3 + xmm1;	/* B3 = I3 * cosine/sine + R3 */ \
\
	xmm0 = xptr(screg2); \
	xmm7 *= xmm0;		/* A2 = A2 * sine (final R2) */ \
	xmm5 *= xmm0;		/* B2 = B2 * sine (final I2) */ \
	xmm1 = xptr(screg5); \
	xmm6 *= xmm1;		/* A4 = A4 * sine (final R4) */ \
	xmm2 *= xmm1;		/* B4 = B4 * sine (final I4) */ \
	xmm0 = xptr(screg1); \
	xmm4 *= xmm0;		/* A3 = A3 * sine (final R3) */ \
	xmm3 *= xmm0;		/* B3 = B3 * sine (final I3) */ \
}

/* 32-bit AMD K8 optimized versions of the above macros */

#if 0 // (@INSTR(,%xarch,<K8>) NE 0);
\
/*R8	AAAA */ \
/*R4	  AAAA */ \
/*R6	    AAAA */ \
/*R2	      AAAA */ \
/*R8*sqrt    MMMM */ \
/*R6*sqrt        MMMM */ \
/*R1	        AAAA */ \
/*R5	          AAAA */ \
/*FI2(R2-R4)	    AAAA */ \
/*FI1(R2+R4)	    MMMMAAAA */ \
/*R6(R6-R8)	      AAAA */ \
/*R8(R6+R8)	      MMMMAAAA */ \
/*R3		            AAAA */ \
/*R7			      AAAA */ \
/*FR2(R1-R3)		        AAAA */ \
/*FR1(R1+R3)		        MMMMAAAA */ \
/*FR4(R5-R6)			  AAAA */ \
/*FR3(R5+R6)			  MMMM  AAAA */ \
/*FI4(R7-R8)			      AAAA */ \
/*FI3(R7+R8)			      MMMMAAAA */ \
/*A2	                            MMMMMMMMAAAA   */ \
/*A4				          MMMMAAAA  MMMM */ \
/*B2	                    MMMM            MMMMAAAA     */ \
/*B4					      MMMMAAAAMMMM */ \
/*A3					        MMMMAAAAMMMM */ \
/*B3					          MMMMAAAAMMMM */ \
\

#define r4_x8r_fft_mem(mem1,mem2,mem3,mem4,mem5,mem6,mem7,mem8,screg1,screg2,screg5,pre1,pre2,dst1,dst2) \
	xmm3 = mem4;		/* R4				; K8 */ \
	xmm7 = mem8;		/* R8 */ \
	xmm3 -= xmm7;		/* new R8 = R4 - R8		; 1-4 */ \
	xmm7 += mem4;		/* new R4 = R4 + R8		; 3-6 */ \
\
	xmm1 = mem2;		/* R2 */ \
	xmm5 = mem6;		/* R6 */ \
	xmm1 -= xmm5;		/* new R6 = R2 - R6		; 5-8 */ \
	xmm5 += mem2;		/* new R2 = R2 + R6		; 7-10 */ \
\
	xmm3 *= XMM_SQRTHALF2;	/* R8 = R8 * square root	; 5-8 */ \
\
	xmm0 = mem1;		/* R1 */ \
	xmm4 = mem5;		/* R5 */ \
	xmm0 -= xmm4;		/* new R5 = R1 - R5		; 9-12 */ \
	xmm4 += mem1;		/* new R1 = R1 + R5		; 11-14 */ \
\
	xmm1 *= XMM_SQRTHALF2;	/* R6 = R6 * square root	; 9-12 */ \
\
	xmm5 -= xmm7;		/* R2 = R2 - R4 (final I2)	; 13-16 */ \
	xmm7 += xmm7; /* multwo */;			/* R4 = R4 * 2			; 13-16 */ \
	xmm1 -= xmm3;		/* R6 = R6 - R8 (Real part)	; 15-18 */ \
	xmm3 += xmm3; /* multwo */;			/* R8 = R8 * 2			; 15-18 */ \
	xmm7 += xmm5;		/* R4 = R2 + R4 (final I1)	; 17-20 */ \
	xmm3 += xmm1;		/* R8 = R6 + R8 (Imaginary part); 19-22 */ \
\
	xmm6 = mem7;		/* R7 */ \
	xmm6 += mem3;		/* new R3 = R3 + R7		; 21-24 */ \
\
	xmm2 = mem3;		/* R3 */ \
	xmm2 -= mem7;		/* new R7 = R3 - R7		; 23-26 */ \
	xmm5 *= xptr(screg2);		/* I2 = I2 * sine		; 23-26 */ \
\
	dst2 = xmm7;		/* Save I1			; 21 */ \
\
	xmm4 -= xmm6;		/* R1 = R1 - R3 (final R2)	; 25-28 */ \
	xmm6 += xmm6; /* multwo */;			/* R3 = R3 * 2			; 25-28 */ \
\
	xprefetchw(u8ptr(pre1)); \
\
	xmm0 -= xmm1;		/* R5 = R5 - R6 (final R4)	; 27-30 */ \
	xmm1 += xmm1; /* multwo */;			/* R6 = R6 * 2			; 27-30 */ \
\
	xmm6 += xmm4;		/* R3 = R1 + R3 (final R1)	; 29-32 */ \
	xmm4 *= xptr(screg2);		/* R2 = R2 * sine		; 29-32 */ \
\
	xprefetchw(u8ptr(pre1+pre2)); \
\
	xmm2 -= xmm3;		/* R7 = R7 - R8 (final I4)	; 31-34 */ \
	xmm3 += xmm3; /* multwo */;			/* R8 = R8 * 2			; 31-34 */ \
	xmm7 = xptr(screg2+16);	/* cosine/sine for w^2n */ \
\
	xmm1 += xmm0;		/* R6 = R5 + R6 (final R3)	; 33-36 */ \
	xmm7 *= xmm4;		/* A2 = R2 * cosine/sine	; 33-36 */ \
	dst1 = xmm6;		/* Save R1			; 33 */ \
\
	xmm3 += xmm2;		/* R8 = R7 + R8 (final I3)	; 35-38 */ \
	xmm6 = xptr(screg5+16);	/* cosine/sine for w^5n */ \
	xmm6 *= xmm0;		/* A4 = R4 * cosine/sine	; 35-38 */ \
\
	xmm7 -= xmm5;		/* A2 = A2 - I2			; 37-40 */ \
	xmm5 *= xptr(screg2+16);	/* B2 = I2 * cosine/sine	; 37-40 */ \
\
	xmm6 -= xmm2;		/* A4 = A4 - I4			; 39-42 */ \
	xmm2 *= xptr(screg5+16);	/* B4 = I4 * cosine/sine	; 39-42 */ \
\
	xmm5 += xmm4;		/* B2 = B2 + R2			; 41-44 */ \
	xmm4 = xptr(screg1+16);	/* cosine/sine for w^n */ \
	xmm4 *= xmm1;		/* A3 = R3 * cosine/sine	; 41-44 */ \
\
	xmm2 += xmm0;		/* B4 = B4 + R4			; 43-46 */ \
\
	xmm4 -= xmm3;		/* A3 = A3 - I3			; 45-48 */ \
	xmm3 *= xptr(screg1+16);	/* B3 = I3 * cosine/sine	; 43-46 */ \
	xmm6 *= xptr(screg5);		/* A4 = A4 * sine (final R4)	; 45-48 */ \
\
	xmm3 += xmm1;		/* B3 = B3 + R3			; 47-50 */ \
	xmm2 *= xptr(screg5);		/* B4 = B4 * sine (final I4)	; 47-50 */ \
	xmm4 *= xptr(screg1);		/* A3 = A3 * sine (final R3)	; 49-52 */ \
	xmm3 *= xptr(screg1);		/* B3 = B3 * sine (final I3)	; 51-54 */ \
\

#endif
\
/* 64-bit Intel and AMD K10 optimized versions of the above macros */ \

#if 0 // (@INSTR(,%xarch,<CORE>) NE 0) OR (@INSTR(,%xarch,<P4>) NE 0) OR (@INSTR(,%xarch,<K10>) NE 0) OR (@INSTR(,%xarch,<BLEND>) NE 0);
#if 0 //DEF X86_64;
\
/* Optimal Core 2 timing is 24 clocks, currently at 27 clocks */ \
\

#define r4_x8r_fft_mem_preload \
	xmm15 = XMM_SQRTHALF2; \
\
\

#define r4_x8r_fft_mem(mem1,mem2,mem3,mem4,mem5,mem6,mem7,mem8,screg1,screg2,screg5,pre1,pre2,dst1,dst2) \
	xmm14 = mem4;		/* R4 */ \
	xmm11 = mem8;		/* R8 */ \
	xmm0 = xmm14;		/* Copy R4 */ \
	xmm14 -= xmm11;		/* new R8 = R4 - R8			; 1-3 */ \
\
	xmm4 = mem2;		/* R2 */ \
	xmm5 = mem6;		/* R6 */ \
	xmm1 = xmm4;		/* Copy R2 */ \
	xmm4 -= xmm5;		/* new R6 = R2 - R6			; 2-4 */ \
\
	xmm11 += xmm0;		/* new R4 = R4 + R8			; 3-5 */ \
\
	xmm5 += xmm1;		/* new R2 = R2 + R6			; 4-6 */ \
	xmm14 *= xmm15;		/* R8 = R8 * square root		; 4-8 */ \
\
	xmm0 = mem1;		/* R1 */ \
	xmm1 = mem5;		/* R5 */ \
	xmm6 = xmm0;		/* Copy R1 */ \
	xmm0 += xmm1;		/* new R1 = R1 + R5			; 5-7 */ \
	xmm4 *= xmm15;		/* R6 = R6 * square root		; 5-9 */ \
\
	xmm13 = mem3;		/* R3 */ \
	xmm8 = mem7;		/* R7 */ \
	xmm9 = xmm13;		/* Copy R3 */ \
	xmm13 += xmm8; 		/* new R3 = R3 + R7			; 6-8 */ \
\
	xmm6 -= xmm1;		/* new R5 = R1 - R5			; 7-9 */ \
\
	xmm1 = xmm5;		/* Copy R2 */ \
	xmm5 -= xmm11;		/* R2 = R2 - R4 (final I2)		; 8-10 */ \
\
	xmm10 = xmm0;		/* Copy R1 */ \
	xmm0 -= xmm13;		/* R1 = R1 - R3 (final R2)		; 9-11 */ \
\
	xmm7 = xmm4;		/* Copy R6 */ \
	xmm4 -= xmm14;		/* R6 = R6 - R8 (Real part)		; 10-12 */ \
\
	xmm14 += xmm7;		/* R8 = R6 + R8 (Imaginary part)	; 11-13 */ \
	xmm7 = xptr(screg2+16);	/* cosine/sine for w^2n */ \
\
	xmm9 -= xmm8;		/* new R7 = R3 - R7			; 12-14 */ \
	xmm2 = xptr(screg5+16);	/* cosine/sine for w^5n */ \
\
	xmm8 = xmm6;		/* Copy R5 */ \
	xmm6 -= xmm4;		/* R5 = R5 - R6 (final R4)		; 13-15 */ \
	xmm3 = xptr(screg1+16);	/* cosine/sine for w^n */ \
\
	xprefetchw(u8ptr(pre1)); \
\
	xmm4 += xmm8;		/* R6 = R5 + R6 (final R3)		; 14-16 */ \
	xmm12 = xmm5;		/* Copy I2 */ \
	xmm5 *= xmm7;		/* B2 = I2 * cosine/sine		; 14-18 */ \
\
	xmm8 = xmm9;		/* Copy R7 */ \
	xmm9 -= xmm14;		/* R7 = R7 - R8 (final I4)		; 15-17 */ \
	xmm7 *= xmm0;		/* A2 = R2 * cosine/sine		; 15-19 */ \
\
	xmm14 += xmm8;		/* R8 = R7 + R8 (final I3)		; 16-18 */ \
	xmm8 = xmm6;		/* Copy R4 */ \
	xmm6 *= xmm2;		/* A4 = R4 * cosine/sine		; 16-20 */ \
\
	xmm11 += xmm1;		/* R4 = R2 + R4 (final I1, a.k.a 2nd real result) ; 17-19 */ \
	xmm1 = xmm4;		/* Copy R3 */ \
	xmm4 *= xmm3;		/* A3 = R3 * cosine/sine		; 17-21 */ \
\
	xmm13 += xmm10;		/* R3 = R1 + R3 (final R1)		; 18-20 */ \
	xmm2 *= xmm9;		/* B4 = I4 * cosine/sine		; 18-22 */ \
	xmm10 = xptr(screg2);		/* sine */ \
\
	xmm5 += xmm0;		/* B2 = B2 + R2				; 19-21 */ \
	xmm3 *= xmm14;		/* B3 = I3 * cosine/sine		; 19-23 */ \
	xmm0 = xptr(screg5); \
\
	xprefetchw(u8ptr(pre1+pre2)); \
\
	xmm7 -= xmm12;		/* A2 = A2 - I2				; 20-22 */ \
	xmm12 = xptr(screg1); \
	dst2 = xmm11;		/* Save I1				; 20 */ \
\
	xmm6 -= xmm9;		/* A4 = A4 - I4				; 21-23 */ \
	dst1 = xmm13;		/* Save R1				; 21 */ \
\
	xmm4 -= xmm14;		/* A3 = A3 - I3				; 22-24 */ \
	xmm5 *= xmm10;		/* B2 = B2 * sine (final I2)		; 22-26 */ \
\
	xmm2 += xmm8;		/* B4 = B4 + R4				; 23-25 */ \
	xmm7 *= xmm10;		/* A2 = A2 * sine (final R2)		; 23-27 */ \
\
	xmm3 += xmm1;		/* B3 = B3 + R3				; 24-26 */ \
	xmm6 *= xmm0;		/* A4 = A4 * sine (final R4)		; 24-28 */ \
\
	xmm4 *= xmm12;		/* A3 = A3 * sine (final R3)		; 25-29 */ \
	xmm2 *= xmm0;		/* B4 = B4 * sine (final I4)		; 26-30 */ \
	xmm3 *= xmm12;		/* B3 = B3 * sine (final I3)		; 27-31 */ \
\

#endif
#endif
\
/* 64-bit AMD K8 optimized versions of the above macros.  Same as 32-bit version with */ \
/* some constants preloaded. */ \

#if 0 // (@INSTR(,%xarch,<K8>) NE 0);
#if 0 //DEF X86_64;
\

#define r4_x8r_fft_mem_preload \
	xmm14 = XMM_SQRTHALF2; \
	xmm15 = XMM_TWO; \
\
\

#define r4_x8r_fft_mem(mem1,mem2,mem3,mem4,mem5,mem6,mem7,mem8,screg1,screg2,screg5,pre1,pre2,dst1,dst2) \
	xmm3 = mem4;		/* R4				; K8 */ \
	xmm7 = mem8;		/* R8 */ \
	xmm3 -= xmm7;		/* new R8 = R4 - R8		; 1-4 */ \
	xmm7 += mem4;		/* new R4 = R4 + R8		; 3-6 */ \
\
	xmm1 = mem2;		/* R2 */ \
	xmm5 = mem6;		/* R6 */ \
	xmm1 -= xmm5;		/* new R6 = R2 - R6		; 5-8 */ \
	xmm5 += mem2;		/* new R2 = R2 + R6		; 7-10 */ \
\
	xmm3 *= xmm14;		/* R8 = R8 * square root	; 5-8 */ \
\
	xmm0 = mem1;		/* R1 */ \
	xmm4 = mem5;		/* R5 */ \
	xmm0 -= xmm4;		/* new R5 = R1 - R5		; 9-12 */ \
	xmm4 += mem1;		/* new R1 = R1 + R5		; 11-14 */ \
\
	xmm1 *= xmm14;		/* R6 = R6 * square root	; 9-12 */ \
\
	xmm5 -= xmm7;		/* R2 = R2 - R4 (final I2)	; 13-16 */ \
	xmm7 *= xmm15;		/* R4 = R4 * 2			; 13-16 */ \
	xmm1 -= xmm3;		/* R6 = R6 - R8 (Real part)	; 15-18 */ \
	xmm3 *= xmm15;		/* R8 = R8 * 2			; 15-18 */ \
	xmm7 += xmm5;		/* R4 = R2 + R4 (final I1)	; 17-20 */ \
	xmm3 += xmm1;		/* R8 = R6 + R8 (Imaginary part); 19-22 */ \
\
	xmm6 = mem7;		/* R7 */ \
	xmm6 += mem3;		/* new R3 = R3 + R7		; 21-24 */ \
\
	xmm2 = mem3;		/* R3 */ \
	xmm2 -= mem7;		/* new R7 = R3 - R7		; 23-26 */ \
	xmm5 *= xptr(screg2);		/* I2 = I2 * sine		; 23-26 */ \
\
	dst2 = xmm7;		/* Save I1			; 21 */ \
\
	xmm4 -= xmm6;		/* R1 = R1 - R3 (final R2)	; 25-28 */ \
	xmm6 *= xmm15;		/* R3 = R3 * 2			; 25-28 */ \
\
	xprefetchw(u8ptr(pre1)); \
\
	xmm0 -= xmm1;		/* R5 = R5 - R6 (final R4)	; 27-30 */ \
	xmm1 *= xmm15;		/* R6 = R6 * 2			; 27-30 */ \
\
	xmm6 += xmm4;		/* R3 = R1 + R3 (final R1)	; 29-32 */ \
	xmm4 *= xptr(screg2);		/* R2 = R2 * sine		; 29-32 */ \
\
	xprefetchw(u8ptr(pre1+pre2)); \
\
	xmm2 -= xmm3;		/* R7 = R7 - R8 (final I4)	; 31-34 */ \
	xmm3 *= xmm15;		/* R8 = R8 * 2			; 31-34 */ \
	xmm7 = xptr(screg2+16);	/* cosine/sine for w^2n */ \
\
	xmm1 += xmm0;		/* R6 = R5 + R6 (final R3)	; 33-36 */ \
	xmm7 *= xmm4;		/* A2 = R2 * cosine/sine	; 33-36 */ \
	dst1 = xmm6;		/* Save R1			; 33 */ \
\
	xmm3 += xmm2;		/* R8 = R7 + R8 (final I3)	; 35-38 */ \
	xmm6 = xptr(screg5+16);	/* cosine/sine for w^5n */ \
	xmm6 *= xmm0;		/* A4 = R4 * cosine/sine	; 35-38 */ \
\
	xmm7 -= xmm5;		/* A2 = A2 - I2			; 37-40 */ \
	xmm5 *= xptr(screg2+16);	/* B2 = I2 * cosine/sine	; 37-40 */ \
\
	xmm6 -= xmm2;		/* A4 = A4 - I4			; 39-42 */ \
	xmm2 *= xptr(screg5+16);	/* B4 = I4 * cosine/sine	; 39-42 */ \
\
	xmm5 += xmm4;		/* B2 = B2 + R2			; 41-44 */ \
	xmm4 = xptr(screg1+16);	/* cosine/sine for w^n */ \
	xmm4 *= xmm1;		/* A3 = R3 * cosine/sine	; 41-44 */ \
\
	xmm2 += xmm0;		/* B4 = B4 + R4			; 43-46 */ \
\
	xmm4 -= xmm3;		/* A3 = A3 - I3			; 45-48 */ \
	xmm3 *= xptr(screg1+16);	/* B3 = I3 * cosine/sine	; 43-46 */ \
	xmm6 *= xptr(screg5);		/* A4 = A4 * sine (final R4)	; 45-48 */ \
\
	xmm3 += xmm1;		/* B3 = B3 + R3			; 47-50 */ \
	xmm2 *= xptr(screg5);		/* B4 = B4 * sine (final I4)	; 47-50 */ \
	xmm4 *= xptr(screg1);		/* A3 = A3 * sine (final R3)	; 49-52 */ \
	xmm3 *= xptr(screg1);		/* B3 = B3 * sine (final I3)	; 51-54 */ \
\

#endif
#endif
\
/* */ \
/* ************************************* eight-reals-unfft variants ****************************************** */ \
/* */ \
\
/* These macros produce eight reals after doing 2 and 3/4 levels of the inverse FFT applying */ \
/* the sin/cos multipliers beforehand.  The input is 2 reals (only 2 levels of inverse FFT done) */ \
/* and 3 complex numbers (3 levels of inverse FFT performed).  These macros take a screg */ \
/* that points to twiddles w^n, w^2n, and w^5n. */ \
\

#define r4_x4cl_eight_reals_last_unfft_preload \
	r4_x4cl_eight_reals_unfft_cmn_preload; \
\
\

#define r4_x4cl_eight_reals_last_unfft(srcreg,srcinc,d1,d2,screg,scoff) \
	r4_x4cl_eight_reals_unfft_cmn(srcreg,srcinc,d1,d2,screg,scoff,screg+32,scoff,screg+64,scoff); \
\
\
/* This macro is similar to r4_x4cl_eight_reals_last_unfft. */ \
/* The difference is this macro takes two screg pointers. */ \
/* One points to w^n and w^5n, the other points to w^2n which is */ \
/* also the screg used by the four-complex-FFT. */ \
\

#define r4_x4cl_eight_reals_unfft_preload \
	r4_x4cl_eight_reals_unfft_cmn_preload; \
\
\

#define r4_x4cl_eight_reals_unfft(srcreg,srcinc,d1,d2,screg1,scoff1,screg2,scoff2) \
	r4_x4cl_eight_reals_unfft_cmn(srcreg,srcinc,d1,d2,screg2,scoff2,screg1,scoff1,screg2+32,scoff2); \
\
\
/* Common eight reals unfft macro */ \
\

#define r4_x4cl_eight_reals_unfft_cmn_preload \
\
\

#define r4_x4cl_eight_reals_unfft_cmn(srcreg,srcinc,d1,d2,screg1,scoff1,screg2,scoff2,screg5,scoff5) { \
	vec2f64 *sd0 = (vec2f64*)(srcreg); \
	vec2f64 *sd1 = (vec2f64*)(srcreg+d1); \
	vec2f64 *sd2 = (vec2f64*)(srcreg+d2); \
	vec2f64 *sd3 = (vec2f64*)(srcreg+d2+d1); \
	vec2f64 xma0,xma1,xma2,xma3,xma4,xma5,xma6,xma7; \
	r4_x8r_unfft_mem(sd0[1],sd0[3],sd1[1],sd1[3],sd2[1],sd2[3],sd3[1],sd3[3], screg1+scoff1,screg2+scoff2, screg5+scoff5,srcreg+srcinc+d2,d1, sd3[1],sd3[3]); \
/*	xstore	[srcreg+d3+16], xmm2	;; Save R4 */ \
/*	xstore	[srcreg+d3+48], xmm5	;; Save R8 */ \
	xma6 = sd3[0];	/* Load R4 */ \
	xma7 = sd3[2];	/* Load I4 */ \
	sd3[0] = xmm4;	/* Save R2 */ \
	sd3[2] = xmm3;	/* Save R6 */ \
	xma2 = sd1[0];	/* Load R2 */ \
	xma3 = sd1[2];	/* Load I2 */ \
	sd1[0] = xmm1;	/* Save R1 */ \
	sd1[1] = xmm7;	/* Save R3 */ \
	sd1[2] = xmm0;	/* Save R5 */ \
	sd1[3] = xmm6;	/* Save R7 */ \
	r4_x8r_unfft_partial_mem(xma0,xma1,xma2,xma3,xma4,xma5,xma6,xma7, sd0[0],sd0[2],sd2[0],sd2[2], screg1,screg2,screg5,srcreg+srcinc,d1, sd2[1],sd2[3]); \
	sd2[0] = xma4;	/* Save R2 */ \
	sd2[2] = xma3;	/* Save R6 */ \
	sd0[0] = xma1;		/* Save R1 */ \
	sd0[1] = xma7;	/* Save R3 */ \
	sd0[2] = xma0;	/* Save R5 */ \
	sd0[3] = xma6;	/* Save R7 */ \
/*	xstore	[srcreg+d2+16], xma2	;; Save R4 */ \
/*	xstore	[srcreg+d2+48], xma5	;; Save R8 */ \
	srcreg += srcinc; \
}

/* Common macro to do the 2 and 3/4 inverse FFT producing eight reals. */ \
\

#define r4_x8r_unfft_mem(mem1,mem2,mem3,mem4,mem5,mem6,mem7,mem8,screg1,screg2,screg5,pre1,pre2,dst4,dst8) \
	xmm6 = xptr(screg5+16);		/* cosine/sine for w^5n */ \
	xmm0 = mem7;				/* R4 */ \
	xmm7 = mem8;				/* I4 */ \
	xmm1 = xmm6 * xmm0 + xmm7;	/* A4 = R4 * cosine/sine + I4 */ \
	xmm7 = xmm6 * xmm7 - xmm0;	/* B4 = I4 * cosine/sine - R4 */ \
\
	xmm2 = xptr(screg1+16);		/* cosine/sine for w^n */ \
	xmm5 = mem5;				/* R3 */ \
	xmm6 = mem6;				/* I3 */ \
	xmm4 = xmm2 * xmm5 + xmm6;	/* A3 = R3 * cosine/sine + I3 */ \
	xmm5 = xmm2 * xmm6 - xmm5;	/* B3 = I3 * cosine/sine - R3 */ \
\
	xmm6 = xptr(screg2+16);		/* cosine/sine for w^2n */ \
	xmm2 = mem3;				/* R2 */ \
	xmm3 = mem4;				/* I2 */ \
	xmm0 = xmm6 * xmm2 + xmm3;	/* A2 = R2 * cosine/sine + I2 */ \
	xmm3 = xmm6 * xmm3 - xmm2;	/* B2 = I2 * cosine/sine - R2 */ \
\
	xprefetchw(u8ptr(pre1)); \
\
	xmm6 = xptr(screg5); \
	xmm1 *= xmm6;		/* new R7 = A4 * sine (final R4) */ \
	xmm7 *= xmm6;		/* new R8 = B4 * sine (final I4) */ \
	xmm2 = xptr(screg1); \
	xmm4 *= xmm2;		/* new R5 = A3 * sine (final R3) */ \
	xmm5 *= xmm2;		/* new R6 = B3 * sine (final I3) */ \
	xmm6 = xptr(screg2); \
	xmm0 *= xmm6;		/* new R3 = A2 * sine (final R2) */ \
	xmm3 *= xmm6;		/* new R4 = B2 * sine (final I2) */ \
\
	xprefetchw(u8ptr(pre1+pre2)); \
\
	xmm6 = xmm4;		/* Copy R5 */ \
	xmm4 -= xmm1;		/* new R6 = R5 - R7 */ \
	xmm1 += xmm6;		/* new R5 = R5 + R7 */ \
\
	xmm2 = xmm5;		/* Copy R6 */ \
	xmm5 -= xmm7;		/* new R8 = R6 - R8 */ \
	xmm7 += xmm2;		/* new R7 = R6 + R8 */ \
\
	xmm6 = xmm4;		/* Copy R6 */ \
	xmm4 += xmm5;		/* R6 = R6 + R8 */ \
	xmm5 -= xmm6;		/* R8 = R8 - R6 */ \
\
	xmm4 *= XMM_SQRTHALF2;	/* R6 = R6 * square root of 1/2 */ \
	xmm5 *= XMM_SQRTHALF2;	/* R8 = R8 * square root of 1/2 */ \
\
	xmm2 = mem2;		/* R2 */ \
	xmm2 -= xmm3;		/* R2 = R2 - R4 (new R4) */ \
	xmm3 += mem2;		/* R4 = R2 + R4 (new R2) */ \
\
	xmm6 = xmm2;		/* Copy R4 */ \
	xmm2 -= xmm5;		/* R4 = R4 - R8 (final R8) */ \
	xmm5 += xmm6;		/* R8 = R4 + R8 (final R4) */ \
\
	xmm6 = xmm3;		/* Copy R2 */ \
	xmm3 -= xmm4;		/* R2 = R2 - R6 (final R6) */ \
	xmm4 += xmm6;		/* R6 = R2 + R6 (final R2) */ \
\
	xmm6 = mem1;		/* R1 */ \
	xmm6 -= xmm0;		/* R1 = R1 - R3 (new R3) */ \
	xmm0 += mem1;		/* R3 = R1 + R3 (new R1) */ \
\
	dst8 = xmm2;		/* Save R8 */ \
	dst4 = xmm5;		/* Save R4 */ \
\
	xmm2 = xmm6;		/* Copy R3 */ \
	xmm6 -= xmm7;		/* R3 = R3 - R7 (final R7) */ \
	xmm7 += xmm2;		/* R7 = R3 + R7 (final R3) */ \
\
	xmm5 = xmm0;		/* Copy R1 */ \
	xmm0 -= xmm1;		/* R1 = R1 - R5 (final R5) */ \
	xmm1 += xmm5;		/* R5 = R1 + R5 (final R1) */ \
\
\

#define r4_x8r_unfft_partial_mem(r1,r2,r3,r4,r5,r6,r7,r8,mem1,mem2,mem5,mem6,screg1,screg2,screg5,pre1,pre2,dst4,dst8) \
	r1 = xptr(screg5+16);	/* cosine/sine for w^5n */ \
	r2 = r1 * r7 + r8;	/* A4 = R4 * cosine/sine + I4 */ \
	r8 = r1 * r8 - r7;	/* B4 = I4 * cosine/sine - R4 */ \
\
	r1 = xptr(screg1+16);		/* cosine/sine for w^n */ \
	r7 = mem5;			/* R3 */ \
	r6 = mem6;			/* I3 */ \
	r5 = r1 * r7 + r6;	/* A3 = R3 * cosine/sine + I3 */ \
	r6 = r1 * r6 - r7;	/* B3 = I3 * cosine/sine - R3 */ \
\
	r7 = xptr(screg2+16);	/* cosine/sine for w^2n */ \
	r1 = r7 * r3 + r4;	/* A2 = R2 * cosine/sine + I2 */ \
	r4 = r7 * r4 - r3;	/* B2 = I2 * cosine/sine - R2 */ \
\
	xprefetchw(u8ptr(pre1)); \
\
	r7 = xptr(screg5); \
	r2 *= r7;			/* new R7 = A4 * sine (final R4) */ \
	r8 *= r7;			/* new R8 = B4 * sine (final I4) */ \
	r3 = xptr(screg1); \
	r5 *= r3;			/* new R5 = A3 * sine (final R3) */ \
	r6 *= r3;			/* new R6 = B3 * sine (final I3) */ \
	r7 = xptr(screg2); \
	r1 *= r7;			/* new R3 = A2 * sine (final R2) */ \
	r4 *= r7;			/* new R4 = B2 * sine (final I2) */ \
\
	xprefetchw(u8ptr(pre1+pre2)); \
\
	r7 = r5;			/* Copy R5 */ \
	r5 -= r2;			/* new R6 = R5 - R7 */ \
	r2 += r7;			/* new R5 = R5 + R7 */ \
\
	r3 = r6;			/* Copy R6 */ \
	r6 -= r8;			/* new R8 = R6 - R8 */ \
	r8 += r3;			/* new R7 = R6 + R8 */ \
\
	r7 = r5;			/* Copy R6 */ \
	r5 += r6;			/* R6 = R6 + R8 */ \
	r6 -= r7;			/* R8 = R8 - R6 */ \
\
	r5 *= XMM_SQRTHALF2;	/* R6 = R6 * square root of 1/2 */ \
	r6 *= XMM_SQRTHALF2;	/* R8 = R8 * square root of 1/2 */ \
\
	r3 = mem2;		/* R2 */ \
	r3 -= r4;			/* R2 = R2 - R4 (new R4) */ \
	r4 += mem2;		/* R4 = R2 + R4 (new R2) */ \
\
	r7 = r3;			/* Copy R4 */ \
	r3 -= r6;			/* R4 = R4 - R8 (final R8) */ \
	r6 += r7;			/* R8 = R4 + R8 (final R4) */ \
\
	r7 = r4;			/* Copy R2 */ \
	r4 -= r5;			/* R2 = R2 - R6 (final R6) */ \
	r5 += r7;			/* R6 = R2 + R6 (final R2) */ \
\
	r7 = mem1;		/* R1 */ \
	r7 -= r1;			/* R1 = R1 - R3 (new R3) */ \
	r1 += mem1;		/* R3 = R1 + R3 (new R1) */ \
\
	dst8 = r3;		/* Save R8 */ \
	dst4 = r6;		/* Save R4 */ \
\
	r3 = r7;			/* Copy R3 */ \
	r7 -= r8;			/* R3 = R3 - R7 (final R7) */ \
	r8 += r3;			/* R7 = R3 + R7 (final R3) */ \
\
	r6 = r1;			/* Copy R1 */ \
	r1 -= r2;			/* R1 = R1 - R5 (final R5) */ \
	r2 += r6;			/* R5 = R1 + R5 (final R1) */ \
\
\
/* This macro is similar to r4_x4cl_eight_reals_unfft. */ \
/* The difference is this macro has different source and destination. */ \
\
/* Used in the last levels of an r4 FFT pass 1 (swizzling): */ \

#define r4_sg2cl_eight_reals_djbunfft(srcreg,srcinc,d1,dstreg,dstinc,e1,screg1,screg2) { \
	vec2f64 *sd0 = (vec2f64*)(srcreg); \
	vec2f64 *sd1 = (vec2f64*)(srcreg+d1); \
	vec2f64 *de0 = (vec2f64*)(dstreg); \
	vec2f64 *de1 = (vec2f64*)(dstreg+e1); \
	r4_s8r_unfft4(sd0[0],sd0[2],sd0[1],sd0[3],sd1[0],sd1[2],sd1[1],sd1[3],screg2,screg1,screg2+32,dstreg+dstinc,e1,de1[1],de1[3]); \
	srcreg += srcinc; \
/*	xstore	[dstreg+e1+48], xmm2	;; Save R8 */ \
/*	xstore	[dstreg+e1+16], xmm5	;; Save R4 */ \
	de1[2] = xmm3;	/* Save R6 */ \
	de1[0] = xmm4;	/* Save R2 */ \
	de0[3] = xmm6;	/* Save R7 */ \
	de0[2] = xmm0;	/* Save R5 */ \
	de0[1] = xmm7;	/* Save R3 */ \
	de0[0] = xmm1;	/* Save R1 */ \
	dstreg += dstinc; \
}

/* The macros below take an XMM_SCD4 pointer.  We know that the first */ \
/* XMM_SCD4 entry is 1, thus we skip over it to get to the w^1 entry we desire. */ \
\
/* This macro is used in the last levels of an r4delay FFT pass 1 (no swizzling): */
#if 0 //IFDEF UNUSED

#define r4_g2cl_eight_reals_unfft4(srcreg,srcinc,d1,dstreg,dstinc,e1,screg1,screg2) \
	r4_x8r_unfft_mem(xptr(srcreg),xptr(srcreg+16),xptr(srcreg+32),xptr(srcreg+48),xptr(srcreg+d1),xptr(srcreg+d1+16),xptr(srcreg+d1+32),xptr(srcreg+d1+48),screg2,screg1+32,screg2+32,dstreg+dstinc,e1,xptr(dstreg+e1+16),xptr(dstreg+e1+48)); \
	srcreg += srcinc; \
	xptr(dstreg+16) = xmm7;	/* Save R3 */ \
/*	xstore	[dstreg+e1+16], xmm2	;; Save R4 */ \
	xptr(dstreg+48) = xmm6;	/* Save R7 */ \
/*	xstore	[dstreg+e1+48], xmm5	;; Save R8 */ \
	xptr(dstreg) = xmm1;		/* Save R1 */ \
	xptr(dstreg+e1) = xmm4;	/* Save R2 */ \
	xptr(dstreg+32) = xmm0;	/* Save R5 */ \
	xptr(dstreg+e1+32) = xmm3;	/* Save R6 */ \
	dstreg += dstinc; \

#endif
\
/* Used in the last levels of an r4delay FFT pass 1 (swizzling): */ \

#define r4_sg2cl_eight_reals_unfft4(srcreg,srcinc,d1,dstreg,dstinc,e1,screg1,screg2) { \
	vec2f64 *sd0 = (vec2f64*)(srcreg); \
	vec2f64 *sd1 = (vec2f64*)(srcreg+d1); \
	vec2f64 *de0 = (vec2f64*)(dstreg); \
	vec2f64 *de1 = (vec2f64*)(dstreg+e1); \
	r4_s8r_unfft4(sd0[0],sd0[2],sd0[1],sd0[3],sd1[0],sd1[2],sd1[1],sd1[3], screg2,screg1+32,screg2+32,dstreg+dstinc,e1,de1[1],de1[3]); \
	srcreg += srcinc; \
/*	xstore	[dstreg+e1+48], xmm2	;; Save R8 */ \
/*	xstore	[dstreg+e1+16], xmm5	;; Save R4 */ \
	de1[2] = xmm3;	/* Save R6 */ \
	de1[0] = xmm4;	/* Save R2 */ \
	de0[3] = xmm6;	/* Save R7 */ \
	de0[2] = xmm0;	/* Save R5 */ \
	de0[1] = xmm7;	/* Save R3 */ \
	de0[0] = xmm1;	/* Save R1 */ \
	dstreg += dstinc; \
}


#define r4_s8r_unfft4(memr1_1,memr1_2,memr2,memi2,memr3,memi3,memr4,memi4,screg1,screg2,screg5,pre1,pre2,dst4,dst8) { \
	shuffle_load_with_temp(xmm0, xmm3, memr3, memr4, xmm6); /* R3,R4 */ \
	shuffle_load_with_temp(xmm5, xmm7, memi3, memi4, xmm6); /* I3,I4 */ \
	xmm6 = xptr(screg5+16);	/* cosine/sine for w^5n */ \
	xmm2 = xptr(screg1+16);	/* cosine/sine for w^n */ \
	xmm1 = xmm6 * xmm3 + xmm7;	/* A4 = R4 * cosine/sine + I4 */ \
	xmm4 = xmm2 * xmm0 + xmm5;	/* A3 = R3 * cosine/sine + I3 */ \
	xmm7 = xmm6 * xmm7 - xmm3;	/* B4 = I4 * cosine/sine - R4 */ \
	xmm5 = xmm2 * xmm5 - xmm0;	/* B3 = I3 * cosine/sine - R3 */ \
\
	shuffle_load_with_temp(xmm2, xmm6, memr1_1, memr2, xmm0); /* Real value #1,R2 */ \
	vec2f64 temp1 = xmm2;		/* Temporarily save real value #1 */ \
	shuffle_load_with_temp(xmm2, xmm3, memr1_2, memi2, xmm0); /* Real value #2,I2 */ \
	vec2f64 s1 = xptr(screg2+16); \
	xmm0 = s1 * xmm6 + xmm3;	/* A2 = R2 * cosine/sine for w^2n + I2 */ \
	xmm3 = s1 * xmm3 - xmm6;	/* B2 = I2 * cosine/sine - R2 */ \
\
	xprefetchw(u8ptr(pre1)); \
\
	xmm6 = xptr(screg5); \
	xmm1 *= xmm6;		/* new R7 = A4 * sine (final R4) */ \
	xmm7 *= xmm6;		/* new R8 = B4 * sine (final I4) */ \
	xmm6 = xptr(screg1); \
	xmm4 *= xmm6;		/* new R5 = A3 * sine (final R3) */ \
	xmm5 *= xmm6;		/* new R6 = B3 * sine (final I3) */ \
	xmm6 = xptr(screg2); \
	xmm0 *= xmm6;		/* new R3 = A2 * sine (final R2) */ \
	xmm3 *= xmm6;		/* new R4 = B2 * sine (final I2) */ \
\
	xprefetchw(u8ptr(pre1+pre2)); \
\
	xmm6 = xmm4;		/* Copy R5 */ \
	xmm4 -= xmm1;		/* new R6 = R5 - R7 */ \
	xmm1 += xmm6;		/* new R5 = R5 + R7 */ \
\
	xmm6 = xmm5;		/* Copy R6 */ \
	xmm5 -= xmm7;		/* new R8 = R6 - R8 */ \
	xmm7 += xmm6;		/* new R7 = R6 + R8 */ \
\
	xmm6 = xmm4;		/* Copy R6 */ \
	xmm4 += xmm5;		/* R6 = R6 + R8 */ \
	xmm5 -= xmm6;		/* R8 = R8 - R6 */ \
\
	xmm4 *= XMM_SQRTHALF2;	/* R6 = R6 * square root of 1/2 */ \
	xmm5 *= XMM_SQRTHALF2;	/* R8 = R8 * square root of 1/2 */ \
\
	xmm6 = xmm2;		/* Copy R2 */ \
	xmm2 -= xmm3;		/* R2 = R2 - R4 (new R4) */ \
	xmm3 += xmm6;		/* R4 = R2 + R4 (new R2) */ \
\
	xmm6 = xmm2;		/* Copy R4 */ \
	xmm2 -= xmm5;		/* R4 = R4 - R8 (final R8) */ \
	xmm5 += xmm6;		/* R8 = R4 + R8 (final R4) */ \
\
	xmm6 = xmm3;		/* Copy R2 */ \
	xmm3 -= xmm4;		/* R2 = R2 - R6 (final R6) */ \
	xmm4 += xmm6;		/* R6 = R2 + R6 (final R2) */ \
\
	xmm6 = temp1;		/* R1 */ \
	xmm6 -= xmm0;		/* R1 = R1 - R3 (new R3) */ \
	xmm0 += temp1;		/* R3 = R1 + R3 (new R1) */ \
\
	dst8 = xmm2;		/* Save R8 */ \
	dst4 = xmm5;		/* Save R4 */ \
\
	xmm2 = xmm6;		/* Copy R3 */ \
	xmm6 -= xmm7;		/* R3 = R3 - R7 (final R7) */ \
	xmm7 += xmm2;		/* R7 = R3 + R7 (final R3) */ \
\
	xmm5 = xmm0;		/* Copy R1 */ \
	xmm0 -= xmm1;		/* R1 = R1 - R5 (final R5) */ \
	xmm1 += xmm5;		/* R5 = R1 + R5 (final R1) */ \
}

/* 32-bit AMD K8 optimized versions of the above macros */ \

#if 0 // (@INSTR(,%xarch,<K8>) NE 0);
\

#define r4_x8r_unfft_mem(mem1,mem2,mem3,mem4,mem5,mem6,mem7,mem8,screg1,screg2,screg5,pre1,pre2,dst4,dst8) \
	xmm6 = xptr(screg5+16);	/* cosine/sine for w^5n */ \
	xmm1 = mem7;		/* R4 */ \
	xmm1 *= xmm6;		/* A4 = R4 * cosine/sine */ \
	xmm7 = mem8;		/* I4 */ \
	xmm1 += xmm7;		/* A4 = A4 + I4 */ \
	xmm7 *= xmm6;		/* B4 = I4 * cosine/sine */ \
\
	xmm2 = xptr(screg1+16);	/* cosine/sine for w^n */ \
	xmm4 = mem5;		/* R3 */ \
	xmm4 *= xmm2;		/* A3 = R3 * cosine/sine for w^n */ \
	xmm5 = mem6;		/* I3 */ \
	xmm4 += xmm5;		/* A3 = A3 + I3 */ \
	xmm5 *= xmm2;		/* B3 = I3 * cosine/sine */ \
\
	xmm6 = xptr(screg2+16);	/* cosine/sine for w^2n */ \
	xmm0 = mem3;		/* R2 */ \
	xmm0 *= xmm6;		/* A2 = R2 * cosine/sine */ \
	xmm3 = mem4;		/* I2 */ \
	xmm0 += xmm3;		/* A2 = A2 + I2 */ \
	xmm3 *= xmm6;		/* B2 = I2 * cosine/sine */ \
\
	xmm7 -= mem7;		/* B4 = B4 - R4 */ \
	xmm5 -= mem5;		/* B3 = B3 - R3 */ \
	xmm3 -= mem3;		/* B2 = B2 - R2 */ \
\
	xprefetchw(u8ptr(pre1)); \
\
	xmm6 = xptr(screg5); \
	xmm1 *= xmm6;		/* new R7 = A4 * sine (final R4) */ \
	xmm7 *= xmm6;		/* new R8 = B4 * sine (final I4) */ \
	xmm2 = xptr(screg1); \
	xmm4 *= xmm2;		/* new R5 = A3 * sine (final R3) */ \
	xmm5 *= xmm2;		/* new R6 = B3 * sine (final I3) */ \
	xmm6 = xptr(screg2); \
	xmm0 *= xmm6;		/* new R3 = A2 * sine (final R2) */ \
	xmm3 *= xmm6;		/* new R4 = B2 * sine (final I2) */ \
\
	xprefetchw(u8ptr(pre1+pre2)); \
\
	xmm4 -= xmm1;		/* new R6 = R5 - R7 */ \
	xmm1 += xmm1; /* multwo */; \
	xmm5 -= xmm7;		/* new R8 = R6 - R8 */ \
	xmm7 += xmm7; /* multwo */; \
	xmm1 += xmm4;		/* new R5 = R5 + R7 */ \
	xmm7 += xmm5;		/* new R7 = R6 + R8 */ \
\
	xmm4 += xmm5;		/* R6 = R6 + R8 */ \
	xmm5 += xmm5; /* multwo */; \
	xmm5 -= xmm4;		/* R8 = R8 - R6 */ \
\
	xmm4 *= XMM_SQRTHALF2;	/* R6 = R6 * square root of 1/2 */ \
	xmm5 *= XMM_SQRTHALF2;	/* R8 = R8 * square root of 1/2 */ \
\
	xmm2 = mem2;		/* R2 */ \
	xmm2 -= xmm3;		/* R2 = R2 - R4 (new R4) */ \
	xmm3 += mem2;		/* R4 = R2 + R4 (new R2) */ \
\
	xmm2 -= xmm5;		/* R4 = R4 - R8 (final R8) */ \
	xmm5 += xmm5; /* multwo */; \
	xmm3 -= xmm4;		/* R2 = R2 - R6 (final R6) */ \
	xmm4 += xmm4; /* multwo */; \
	xmm5 += xmm2;		/* R8 = R4 + R8 (final R4) */ \
	xmm4 += xmm3;		/* R6 = R2 + R6 (final R2) */ \
\
	xmm6 = mem1;		/* R1 */ \
	xmm6 -= xmm0;		/* R1 = R1 - R3 (new R3) */ \
	xmm0 += mem1;		/* R3 = R1 + R3 (new R1) */ \
\
	dst8 = xmm2;		/* Save R8 */ \
	dst4 = xmm5;		/* Save R4 */ \
\
	xmm6 -= xmm7;		/* R3 = R3 - R7 (final R7) */ \
	xmm7 += xmm7; /* multwo */; \
	xmm0 -= xmm1;		/* R1 = R1 - R5 (final R5) */ \
	xmm1 += xmm1; /* multwo */; \
	xmm7 += xmm6;		/* R7 = R3 + R7 (final R3) */ \
	xmm1 += xmm0;		/* R5 = R1 + R5 (final R1) */ \
\
\

#define r4_x8r_unfft_partial_mem(r1,r2,r3,r4,r5,r6,r7,r8,mem1,mem2,mem5,mem6,screg1,screg2,screg5,pre1,pre2,dst4,dst8) \
	r2 = xptr(screg5+16);		/* cosine/sine for w^5n */ \
	r2 *= r7;			/* A4 = R4 * cosine/sine */ \
	r2 += r8;			/* A4 = A4 + I4 */ \
	r8 *= xptr(screg5+16);		/* B4 = I4 * cosine/sine */ \
	r8 -= r7;			/* B4 = B4 - R4 */ \
\
	r1 = xptr(screg1+16);		/* cosine/sine for w^n */ \
	r5 = mem5;		/* R3 */ \
	r5 *= r1;			/* A3 = R3 * cosine/sine for w^n */ \
	r6 = mem6;		/* I3 */ \
	r5 += r6;			/* A3 = A3 + I3 */ \
	r6 *= r1;			/* B3 = I3 * cosine/sine */ \
	r6 -= mem5;		/* B3 = B3 - R3 */ \
\
	r1 = xptr(screg2+16);		/* cosine/sine for w^2n */ \
	r1 *= r3;			/* A2 = R2 * cosine/sine */ \
	r1 += r4;			/* A2 = A2 + I2 */ \
	r4 *= xptr(screg2+16);		/* B2 = I2 * cosine/sine */ \
	r4 -= r3;			/* B2 = B2 - R2 */ \
\
	xprefetchw(u8ptr(pre1)); \
\
	r7 = xptr(screg5); \
	r2 *= r7;			/* new R7 = A4 * sine (final R4) */ \
	r8 *= r7;			/* new R8 = B4 * sine (final I4) */ \
	r3 = xptr(screg1); \
	r5 *= r3;			/* new R5 = A3 * sine (final R3) */ \
	r6 *= r3;			/* new R6 = B3 * sine (final I3) */ \
	r7 = xptr(screg2); \
	r1 *= r7;			/* new R3 = A2 * sine (final R2) */ \
	r4 *= r7;			/* new R4 = B2 * sine (final I2) */ \
\
	xprefetchw(u8ptr(pre1+pre2)); \
\
	r5 -= r2;			/* new R6 = R5 - R7 */ \
	r2 += r2; /* multwo */ \
	r6 -= r8;			/* new R8 = R6 - R8 */ \
	r8 += r8; /* multwo */ \
	r2 += r5;			/* new R5 = R5 + R7 */ \
	r8 += r6;			/* new R7 = R6 + R8 */ \
\
	r5 += r6;			/* R6 = R6 + R8 */ \
	r6 += r6; /* multwo */ \
	r6 -= r5;			/* R8 = R8 - R6 */ \
\
	r5 *= XMM_SQRTHALF2;	/* R6 = R6 * square root of 1/2 */ \
	r6 *= XMM_SQRTHALF2;	/* R8 = R8 * square root of 1/2 */ \
\
	r3 = mem2;		/* R2 */ \
	r3 -= r4;			/* R2 = R2 - R4 (new R4) */ \
	r4 += mem2;		/* R4 = R2 + R4 (new R2) */ \
\
	r3 -= r6;			/* R4 = R4 - R8 (final R8) */ \
	r6 += r6; /* multwo */ \
	r4 -= r5;			/* R2 = R2 - R6 (final R6) */ \
	r5 += r5; /* multwo */ \
	r6 += r3;			/* R8 = R4 + R8 (final R4) */ \
	r5 += r4;			/* R6 = R2 + R6 (final R2) */ \
\
	r7 = mem1;		/* R1 */ \
	r7 -= r1;			/* R1 = R1 - R3 (new R3) */ \
	r1 += mem1;		/* R3 = R1 + R3 (new R1) */ \
\
	dst8 = r3;		/* Save R8 */ \
	dst4 = r6;		/* Save R4 */ \
\
	r7 -= r8;			/* R3 = R3 - R7 (final R7) */ \
	r8 += r8; /* multwo */ \
	r1 -= r2;			/* R1 = R1 - R5 (final R5) */ \
	r2 += r2; /* multwo */ \
	r8 += r7;			/* R7 = R3 + R7 (final R3) */ \
	r2 += r1;			/* R5 = R1 + R5 (final R1) */ \
\

#endif
\
/* 64-bit Intel and AMD K10 optimized versions of the above macros */ \

#if 0 // (@INSTR(,%xarch,<CORE>) NE 0) OR (@INSTR(,%xarch,<P4>) NE 0) OR (@INSTR(,%xarch,<K10>) NE 0) OR (@INSTR(,%xarch,<BLEND>) NE 0);
#if 0 //DEF X86_64;
\

#define r4_x4cl_eight_reals_unfft_cmn_preload \
	xmm15 = XMM_SQRTHALF2; \
\
\
/* Theoretical best Core 2 timing is 48 clocks.  The empty add slot at clock 14 means */ \
/* this macro should run in 49 clocks.  Actual timing is 58 clocks. */ \
\

#define r4_x4cl_eight_reals_unfft_cmn(srcreg,srcinc,d1,d2,screg1,scoff1,screg2,scoff2,screg5,scoff5) \
	xmm0 = xptr(srcreg+d2+d1+16);	/* R4 */ \
	xmm1 = xptr(screg5+scoff5+16);/* cosine/sine for w^5n */ \
	xmm2 = xmm0;		/* Copy R4 */ \
	xmm0 *= xmm1;		/* A4 = R4 * cosine/sine		; 1-5 */ \
\
	xmm3 = xptr(srcreg+d2+16);	/* R3 */ \
	xmm4 = xptr(screg1+scoff1+16);/* cosine/sine for w^n */ \
	xmm5 = xmm3;		/* Copy R3 */ \
	xmm3 *= xmm4;		/* A3 = R3 * cosine/sine for w^n	; 2-6 */ \
\
	xmm6 = xptr(srcreg+d2+d1+48);	/* I4 */ \
	xmm1 *= xmm6;		/* B4 = I4 * cosine/sine		; 3-7 */ \
\
	xmm7 = xptr(srcreg+d2+48);	/* I3 */ \
	xmm4 *= xmm7;		/* B3 = I3 * cosine/sine		; 4-8 */ \
\
	xmm8 = xptr(srcreg+d2+d1);	/*#2 R4 */ \
	xmm9 = xptr(screg5+16);	/*#2 cosine/sine for w^5n */ \
	xmm10 = xmm8;		/*#2 Copy R4 */ \
	xmm8 *= xmm9;		/*#2 A4 = R4 * cosine/sine		; 5-9 */ \
\
	xmm0 += xmm6;		/* A4 = A4 + I4				; 6-8	avail 6,11+ */ \
	xmm11 = xptr(srcreg+d2);	/*#2 R3 */ \
	xmm12 = xptr(screg1+16);	/*#2 cosine/sine for w^n */ \
	xmm13 = xmm11;		/*#2 Copy R3 */ \
	xmm11 *= xmm12;		/*#2 A3 = R3 * cosine/sine for w^n	; 6-10	avail 6,14 */ \
\
	xmm3 += xmm7;		/* A3 = A3 + I3				; 7-9	avail 6,7,14 */ \
	xmm14 = xptr(srcreg+d2+d1+32);/*#2 I4 */ \
	xmm9 *= xmm14;		/*#2 B4 = I4 * cosine/sine		; 7-11	avail 6,7 */ \
\
	xmm1 -= xmm2;		/* B4 = B4 - R4				; 8-10	avail 6,7,2 */ \
	xmm6 = xptr(srcreg+d2+32);	/*#2 I3 */ \
	xmm12 *= xmm6;		/*#2 B3 = I3 * cosine/sine		; 8-12	avail 7,2 */ \
\
	xmm4 -= xmm5;		/* B3 = B3 - R3				; 9-11	avail 7,2,5 */ \
	xmm7 = xptr(screg5+scoff5);	/* sine */ \
	xmm0 *= xmm7;		/* new R7 = A4 * sine (final R4)	; 9-13	avail 2,5 */ \
\
	xmm8 += xmm14;		/*#2 A4 = A4 + I4			; 10-12	avail 2,5,14 */ \
	xmm2 = xptr(screg1+scoff1);	/* sine */ \
	xmm3 *= xmm2;		/* new R5 = A3 * sine (final R3)	; 10-14	avail 5,14 */ \
\
	xmm11 += xmm6;		/*#2 A3 = A3 + I3			; 11-13	avail 5,14,6 */ \
	xmm1 *= xmm7;		/* new R8 = B4 * sine (final I4)	; 11-15	avail 5,14,6,7 */ \
\
	xmm9 -= xmm10;		/*#2 B4 = B4 - R4			; 12-14	avail 5,14,6,7,10 */ \
	xmm4 *= xmm2;		/* new R6 = B3 * sine (final I3)	; 12-16	avail 5,14,6,7,10,2 */ \
\
	xmm12 -= xmm13;		/*#2 B3 = B3 - R3			; 13-15	avail 5,14,6,7,10,2,13 */ \
	xmm5 = xptr(srcreg+d1+16);	/* R2 */ \
	xmm14 = xptr(screg2+scoff2+16);/* cosine/sine for w^2n */ \
	xmm6 = xmm5;		/* Copy R2 */ \
	xmm5 *= xmm14;		/* A2 = R2 * cosine/sine		; 13-17	avail 7,10,2,13 */ \
\
	xmm7 = xptr(srcreg+d1+48);	/* I2 */ \
	xmm14 *= xmm7;		/* B2 = I2 * cosine/sine		; 14-18	avail 10,2,13 */ \
\
	xprefetchw(u8ptr(srcreg+srcinc)); \
\
	xmm10 = xmm3;		/* Copy R5 */ \
	xmm3 -= xmm0;		/* new R6 = R5 - R7			; 15-17	avail 2,13 */ \
	xmm2 = xptr(screg5);		/*#2 sine */ \
	xmm8 *= xmm2;		/*#2 new R7 = A4 * sine (final R4)	; 15-19	avail 13 */ \
\
	xmm10 += xmm0;		/* new R5 = R5 + R7			; 16-18	avail 13,0 */ \
	xmm13 = xptr(screg1);		/*#2 sine */ \
	xmm11 *= xmm13;		/*#2 new R5 = A3 * sine (final R3)	; 16-20	avail 0 */ \
\
	xmm0 = xmm4;		/* Copy R6 */ \
	xmm4 -= xmm1;		/* new R8 = R6 - R8			; 17-19 */ \
	xmm9 *= xmm2;		/*#2 new R8 = B4 * sine (final I4)	; 17-21 avail 2 */ \
\
	xmm5 += xmm7;		/* A2 = A2 + I2				; 18-20	avail 2,7 */ \
	xmm12 *= xmm13;		/*#2 new R6 = B3 * sine (final I3)	; 18-22	avail 2,7,13 */ \
\
	xmm0 += xmm1;		/* new R7 = R6 + R8			; 19-21	avail 2,7,13,1 */ \
	xmm3 *= xmm15;		/* R6 = R6 * square root of 1/2		; 19-23 */ \
\
	xmm14 -= xmm6;		/* B2 = B2 - R2				; 20-22	avail 2,7,13,1,6 */ \
	xmm4 *= xmm15;		/* R8 = R8 * square root of 1/2		; 20-24 */ \
\
	xprefetchw(u8ptr(srcreg+srcinc+d1)); \
\
	xmm2 = xmm11;		/*#2 Copy R5 */ \
	xmm11 -= xmm8;		/*#2 new R6 = R5 - R7			; 21-23	avail 7,13,1,6 */ \
\
	xmm2 += xmm8;		/*#2 new R5 = R5 + R7			; 22-24	avail 7,13,1,6,8 */ \
	xmm7 = xptr(screg2+scoff2);	/* sine */ \
	xmm5 *= xmm7;		/* new R3 = A2 * sine (final R2)	; 22-26	avail 13,1,6,8 */ \
\
	xmm13 = xmm12;		/*#2 Copy R6 */ \
	xmm12 -= xmm9;		/*#2 new R8 = R6 - R8			; 23-25	avail 1,6,8 */ \
	xmm14 *= xmm7;		/* new R4 = B2 * sine (final I2)	; 23-27	avail 1,6,8,7 */ \
\
	xmm13 += xmm9;		/*#2 new R7 = R6 + R8			; 24-26	avail 1,6,8,7,9 */ \
	xmm6 = xptr(srcreg+16);	/* R1 */ \
\
	xmm1 = xmm3;		/* Copy R6 */ \
	xmm3 += xmm4;		/* R6 = R6 + R8				; 25-27	avail 8,7,9 */ \
	xmm7 = xptr(srcreg+48);	/* R2 */ \
\
	xmm4 -= xmm1;		/* R8 = R8 - R6				; 26-28	avail 8,9,1 */ \
\
	xprefetchw(u8ptr(srcreg+srcinc+d2));	/*#2 */ \
\
	xmm8 = xmm6;		/* Copy R1 */ \
	xmm6 -= xmm5;		/* R1 = R1 - R3 (new R3)		; 27-29	avail 9,1 */ \
	xmm8 += xmm5;		/* R3 = R1 + R3 (new R1)		; 28-30	avail 9,1,5 */ \
\
	xmm9 = xmm7;		/* Copy R2 */ \
	xmm7 -= xmm14;		/* R2 = R2 - R4 (new R4)		; 29-31	avail 1,5 */ \
\
	xmm14 += xmm9;		/* R4 = R2 + R4 (new R2)		; 30-32	avail 1,5,9 */ \
	xmm1 = xptr(srcreg+d1);	/*#2 R2 */ \
	xmm5 = xptr(screg2+16);	/*#2 cosine/sine for w^2n */ \
	xmm1 *= xmm5;		/*#2 A2 = R2 * cosine/sine		; 30-34	avail 9 */ \
\
	xmm9 = xmm6;		/* Copy R3 */ \
	xmm6 -= xmm0;		/* R3 = R3 - R7 (final R7)		; 31-33	avail none storable 6 */ \
	xptr(srcreg+d1+48) = xmm6;	/* Save R7				; 34	avail 6 */ \
	xmm6 = xptr(srcreg+d1+32);	/*#2 I2 */ \
	xmm5 *= xmm6;		/*#2 B2 = I2 * cosine/sine		; 31-35	avail none */ \
\
	xmm0 += xmm9;		/* R7 = R3 + R7 (final R3)		; 32-34	avail 9 storable 0 */ \
\
	xmm9 = xmm8;		/* Copy R1 */ \
	xmm8 -= xmm10;		/* R1 = R1 - R5 (final R5)		; 33-35	avail none storable 0,8 */ \
	xmm10 += xmm9;		/* R5 = R1 + R5 (final R1)		; 34-36	avail 9 storable 0,8,10 */ \
\
	xmm1 += xmm6;		/*#2 A2 = A2 + I2			; 35-37	avail 9,6 storable 0,8,10 */ \
	xmm11 *= xmm15;		/*#2 R6 = R6 * square root of 1/2	; 35-39 */ \
	xptr(srcreg+d1+16) = xmm0;	/* Save R3				; 35	avail 9,6,0 storable 8,10 */ \
\
	xmm5 -= xptr(srcreg+d1);	/*#2 B2 = B2 - R2			; 36-38 */ \
	xmm12 *= xmm15;		/*#2 R8 = R8 * square root of 1/2	; 36-40 */ \
	xptr(srcreg+d1+32) = xmm8;	/* Save R5				; 36	avail 9,6,0,8 storable 10 */ \
\
	xmm9 = xmm7;		/* Copy R4 */ \
	xmm7 -= xmm4;		/* R4 = R4 - R8 (final R8)		; 37-39	avail 6,0,8 storable 10,7 */ \
	xptr(srcreg+d1) = xmm10;	/* Save R1				; 37	avail 6,0,8,10 storable 7 */ \
\
	xprefetchw(u8ptr(srcreg+srcinc+d2+d1));/*#2 */ \
\
	xmm4 += xmm9;		/* R8 = R4 + R8 (final R4)		; 38-40	avail 6,0,8,10,9 storable 7,4 */ \
	xmm6 = xptr(screg2);		/*#2 sine */ \
	xmm1 *= xmm6;		/*#2 new R3 = A2 * sine (final R2)	; 38-42	avail 0,8,10,9 storable 7,4 */ \
\
	xmm0 = xmm14;		/* Copy R2 */ \
	xmm14 -= xmm3;		/* R2 = R2 - R6 (final R6)		; 39-41	avail 8,10,9 storable 7,4,14 */ \
	xmm5 *= xmm6;		/*#2 new R4 = B2 * sine (final I2)	; 39-43	avail 8,10,9,6 storable 7,4,14 */ \
\
	xmm3 += xmm0;		/* R6 = R2 + R6 (final R2)		; 40-42	avail 8,10,9,6,0 storable 7,4,14,3 */ \
	xmm10 = xptr(srcreg);		/*#2 R1 */ \
	xptr(srcreg+d2+d1+48) = xmm7;	/* Save R8				; 40	avail 8,9,6,0,7 storable 4,14,3 */ \
\
	xmm8 = xmm11;		/*#2 Copy R6 */ \
	xmm11 += xmm12;		/*#2 R6 = R6 + R8			; 41-43	avail 9,6,0,7 storable 4,14,3 */ \
	xptr(srcreg+d2+d1+16) = xmm4;	/* Save R4				; 41	avail 9,6,0,7,4 storable 14,3 */ \
\
	xmm12 -= xmm8;		/*#2 R8 = R8 - R6			; 42-44	avail 9,6,0,7,4,8 storable 14,3 */ \
	xmm6 = xptr(srcreg+32);	/*#2 R2 */ \
	xptr(srcreg+d2+d1+32) = xmm14;/* Save R6				; 42	avail 9,0,7,4,8,14 storable 3 */ \
\
	xmm9 = xmm10;		/*#2 Copy R1 */ \
	xmm10 -= xmm1;		/*#2 R1 = R1 - R3 (new R3)		; 43-45	avail 0,7,4,8,14 storable 3 */ \
	xptr(srcreg+d2+d1) = xmm3;	/* Save R2				; 43	avail 0,7,4,8,14,3 */ \
\
	xmm1 += xmm9;		/*#2 R3 = R1 + R3 (new R1)		; 44-46	avail 0,7,4,8,14,3,9 */ \
\
	xmm0 = xmm6;		/*#2 Copy R2 */ \
	xmm6 -= xmm5;		/*#2 R2 = R2 - R4 (new R4)		; 45-47	avail 7,4,8,14,3,9 */ \
\
	xmm5 += xmm0;		/*#2 R4 = R2 + R4 (new R2)		; 46-48	avail 7,4,8,14,3,9,0 */ \
\
	xmm7 = xmm10;		/*#2 Copy R3 */ \
	xmm10 -= xmm13;		/*#2 R3 = R3 - R7 (final R7)		; 47-49	avail 4,8,14,3,9,0 storable 10 */ \
	xmm13 += xmm7;		/*#2 R7 = R3 + R7 (final R3)		; 48-50	avail 4,8,14,3,9,0,7 storable 10,13 */ \
\
	xmm4 = xmm1;		/*#2 Copy R1 */ \
	xmm1 -= xmm2;		/*#2 R1 = R1 - R5 (final R5)		; 49-51	avail 8,14,3,9,0,7 storable 10,13,1 */ \
\
	xmm2 += xmm4;		/*#2 R5 = R1 + R5 (final R1)		; 50-52	avail 8,14,3,9,0,7,4 storable 10,13,1,2 */ \
	xptr(srcreg+48) = xmm10;	/*#2 Save R7				; 50	avail 8,14,3,9,0,7,4,10 storable 13,1,2 */ \
\
	xmm8 = xmm6;		/*#2 Copy R4 */ \
	xmm6 -= xmm12;		/*#2 R4 = R4 - R8 (final R8)		; 51-53	avail 14,3,9,0,7,4,10 storable 13,1,2,6 */ \
	xptr(srcreg+16) = xmm13;	/*#2 Save R3				; 51	avail 14,3,9,0,7,4,10,13 storable 1,2,6 */ \
\
	xmm12 += xmm8;		/*#2 R8 = R4 + R8 (final R4)		; 52-54	avail 14,3,9,0,7,4,10,13,8 storable 1,2,6,12 */ \
	xptr(srcreg+32) = xmm1;	/*#2 Save R5				; 52	avail 14,3,9,0,7,4,10,13,8,1 storable 2,6,12 */ \
\
	xmm8 = xmm5;		/*#2 Copy R2 */ \
	xmm5 -= xmm11;		/*#2 R2 = R2 - R6 (final R6)		; 53-55	avail 14,3,9,0,7,4,10,13,1 storable 2,6,12,5 */ \
	xptr(srcreg) = xmm2;		/*#2 Save R1				; 53	avail 14,3,9,0,7,4,10,13,1,2 storable 6,12,5 */ \
\
	xmm11 += xmm8;		/*#2 R6 = R2 + R6 (final R2)		; 54-56	avail 14,3,9,0,7,4,10,13,1,2,8 storable 6,12,5,11 */ \
	xptr(srcreg+d2+48) = xmm6;	/*#2 Save R8				; 54 */ \
	xptr(srcreg+d2+16) = xmm12;	/*#2 Save R4				; 55 */ \
	xptr(srcreg+d2+32) = xmm5;	/*#2 Save R6				; 56 */ \
	xptr(srcreg+d2) = xmm11;	/*#2 Save R2				; 57 */ \
\
	srcreg += srcinc; \
\

#endif
#endif
\
/* 64-bit AMD K8 optimized versions of the above macros.  Same as 32-bit version with */ \
/* some constants preloaded. NOT REALLY ANY BETTER!!  We should try rewritting with */ \
/* combining r4_x8r_unfft_mem and r4_x8r_unfft_partial_mem into one macro.  Should */ \
/* yield significant savings. */ \

#if 0 // (@INSTR(,%xarch,<K8>) NE 0);
#if 0 //DEF X86_64;
\

#define r4_x4cl_eight_reals_unfft_cmn_preload \
	xmm14 = XMM_SQRTHALF2; \
	xmm15 = XMM_TWO; \
\
\

#define r4_x8r_unfft_mem(mem1,mem2,mem3,mem4,mem5,mem6,mem7,mem8,screg1,screg2,screg5,pre1,pre2,dst4,dst8) \
	xmm6 = xptr(screg5+16);	/* cosine/sine for w^5n */ \
	xmm1 = mem7;		/* R4 */ \
	xmm1 *= xmm6;		/* A4 = R4 * cosine/sine */ \
	xmm7 = mem8;		/* I4 */ \
	xmm1 += xmm7;		/* A4 = A4 + I4 */ \
	xmm7 *= xmm6;		/* B4 = I4 * cosine/sine */ \
\
	xmm2 = xptr(screg1+16);	/* cosine/sine for w^n */ \
	xmm4 = mem5;		/* R3 */ \
	xmm4 *= xmm2;		/* A3 = R3 * cosine/sine for w^n */ \
	xmm5 = mem6;		/* I3 */ \
	xmm4 += xmm5;		/* A3 = A3 + I3 */ \
	xmm5 *= xmm2;		/* B3 = I3 * cosine/sine */ \
\
	xmm6 = xptr(screg2+16);	/* cosine/sine for w^2n */ \
	xmm0 = mem3;		/* R2 */ \
	xmm0 *= xmm6;		/* A2 = R2 * cosine/sine */ \
	xmm3 = mem4;		/* I2 */ \
	xmm0 += xmm3;		/* A2 = A2 + I2 */ \
	xmm3 *= xmm6;		/* B2 = I2 * cosine/sine */ \
\
	xmm7 -= mem7;		/* B4 = B4 - R4 */ \
	xmm5 -= mem5;		/* B3 = B3 - R3 */ \
	xmm3 -= mem3;		/* B2 = B2 - R2 */ \
\
	xprefetchw(u8ptr(pre1)); \
\
	xmm6 = xptr(screg5); \
	xmm1 *= xmm6;		/* new R7 = A4 * sine (final R4) */ \
	xmm7 *= xmm6;		/* new R8 = B4 * sine (final I4) */ \
	xmm2 = xptr(screg1); \
	xmm4 *= xmm2;		/* new R5 = A3 * sine (final R3) */ \
	xmm5 *= xmm2;		/* new R6 = B3 * sine (final I3) */ \
	xmm6 = xptr(screg2); \
	xmm0 *= xmm6;		/* new R3 = A2 * sine (final R2) */ \
	xmm3 *= xmm6;		/* new R4 = B2 * sine (final I2) */ \
\
	xprefetchw(u8ptr(pre1+pre2)); \
\
	xmm4 -= xmm1;		/* new R6 = R5 - R7 */ \
	xmm1 *= xmm15;		/* R7 = R7 * 2 */ \
	xmm5 -= xmm7;		/* new R8 = R6 - R8 */ \
	xmm7 *= xmm15;		/* R8 = R8 * 2 */ \
	xmm1 += xmm4;		/* new R5 = R5 + R7 */ \
	xmm7 += xmm5;		/* new R7 = R6 + R8 */ \
\
	xmm4 += xmm5;		/* R6 = R6 + R8 */ \
	xmm5 *= xmm15;		/* R8 = R8 * 2 */ \
	xmm5 -= xmm4;		/* R8 = R8 - R6 */ \
\
	xmm4 *= xmm14;		/* R6 = R6 * square root of 1/2 */ \
	xmm5 *= xmm14;		/* R8 = R8 * square root of 1/2 */ \
\
	xmm2 = mem2;		/* R2 */ \
	xmm2 -= xmm3;		/* R2 = R2 - R4 (new R4) */ \
	xmm3 += mem2;		/* R4 = R2 + R4 (new R2) */ \
\
	xmm2 -= xmm5;		/* R4 = R4 - R8 (final R8) */ \
	xmm5 *= xmm15;		/* R8 = R8 * 2 */ \
	xmm3 -= xmm4;		/* R2 = R2 - R6 (final R6) */ \
	xmm4 *= xmm15;		/* R6 = R6 * 2 */ \
	xmm5 += xmm2;		/* R8 = R4 + R8 (final R4) */ \
	xmm4 += xmm3;		/* R6 = R2 + R6 (final R2) */ \
\
	xmm6 = mem1;		/* R1 */ \
	xmm6 -= xmm0;		/* R1 = R1 - R3 (new R3) */ \
	xmm0 += mem1;		/* R3 = R1 + R3 (new R1) */ \
\
	dst8 = xmm2;		/* Save R8 */ \
	dst4 = xmm5;		/* Save R4 */ \
\
	xmm6 -= xmm7;		/* R3 = R3 - R7 (final R7) */ \
	xmm7 *= xmm15;		/* R7 = R7 * 2 */ \
	xmm0 -= xmm1;		/* R1 = R1 - R5 (final R5) */ \
	xmm1 *= xmm15;		/* R5 = R5 * 2 */ \
	xmm7 += xmm6;		/* R7 = R3 + R7 (final R3) */ \
	xmm1 += xmm0;		/* R5 = R1 + R5 (final R1) */ \
\
\

#define r4_x8r_unfft_partial_mem(r1,r2,r3,r4,r5,r6,r7,r8,mem1,mem2,mem5,mem6,screg1,screg2,screg5,pre1,pre2,dst4,dst8) \
	r2 = xptr(screg5+16);		/* cosine/sine for w^5n */ \
	r2 *= r7;			/* A4 = R4 * cosine/sine */ \
	r2 += r8;			/* A4 = A4 + I4 */ \
	r8 *= xptr(screg5+16);		/* B4 = I4 * cosine/sine */ \
	r8 -= r7;			/* B4 = B4 - R4 */ \
\
	r1 = xptr(screg1+16);		/* cosine/sine for w^n */ \
	r5 = mem5;		/* R3 */ \
	r5 *= r1;			/* A3 = R3 * cosine/sine for w^n */ \
	r6 = mem6;		/* I3 */ \
	r5 += r6;			/* A3 = A3 + I3 */ \
	r6 *= r1;			/* B3 = I3 * cosine/sine */ \
	r6 -= mem5;		/* B3 = B3 - R3 */ \
\
	r1 = xptr(screg2+16);		/* cosine/sine for w^2n */ \
	r1 *= r3;			/* A2 = R2 * cosine/sine */ \
	r1 += r4;			/* A2 = A2 + I2 */ \
	r4 *= xptr(screg2+16);		/* B2 = I2 * cosine/sine */ \
	r4 -= r3;			/* B2 = B2 - R2 */ \
\
	xprefetchw(u8ptr(pre1)); \
\
	r7 = xptr(screg5); \
	r2 *= r7;			/* new R7 = A4 * sine (final R4) */ \
	r8 *= r7;			/* new R8 = B4 * sine (final I4) */ \
	r3 = xptr(screg1); \
	r5 *= r3;			/* new R5 = A3 * sine (final R3) */ \
	r6 *= r3;			/* new R6 = B3 * sine (final I3) */ \
	r7 = xptr(screg2); \
	r1 *= r7;			/* new R3 = A2 * sine (final R2) */ \
	r4 *= r7;			/* new R4 = B2 * sine (final I2) */ \
\
	xprefetchw(u8ptr(pre1+pre2)); \
\
	r5 -= r2;			/* new R6 = R5 - R7 */ \
	r2 *= xmm15;		/* R7 = R7 * 2 */ \
	r6 -= r8;			/* new R8 = R6 - R8 */ \
	r8 *= xmm15;		/* R8 = R8 * 2 */ \
	r2 += r5;			/* new R5 = R5 + R7 */ \
	r8 += r6;			/* new R7 = R6 + R8 */ \
\
	r5 += r6;			/* R6 = R6 + R8 */ \
	r6 *= xmm15;		/* R8 = R8 * 2 */ \
	r6 -= r5;			/* R8 = R8 - R6 */ \
\
	r5 *= xmm14;		/* R6 = R6 * square root of 1/2 */ \
	r6 *= xmm14;		/* R8 = R8 * square root of 1/2 */ \
\
	r3 = mem2;		/* R2 */ \
	r3 -= r4;			/* R2 = R2 - R4 (new R4) */ \
	r4 += mem2;		/* R4 = R2 + R4 (new R2) */ \
\
	r3 -= r6;			/* R4 = R4 - R8 (final R8) */ \
	r6 *= xmm15;		/* R8 = R8 * 2 */ \
	r4 -= r5;			/* R2 = R2 - R6 (final R6) */ \
	r5 *= xmm15;		/* R6 = R6 * 2 */ \
	r6 += r3;			/* R8 = R4 + R8 (final R4) */ \
	r5 += r4;			/* R6 = R2 + R6 (final R2) */ \
\
	r7 = mem1;		/* R1 */ \
	r7 -= r1;			/* R1 = R1 - R3 (new R3) */ \
	r1 += mem1;		/* R3 = R1 + R3 (new R1) */ \
\
	dst8 = r3;		/* Save R8 */ \
	dst4 = r6;		/* Save R4 */ \
\
	r7 -= r8;			/* R3 = R3 - R7 (final R7) */ \
	r8 *= xmm15;		/* R7 = R7 * 2 */ \
	r1 -= r2;			/* R1 = R1 - R5 (final R5) */ \
	r2 *= xmm15;		/* R5 = R5 * 2 */ \
	r8 += r7;			/* R7 = R3 + R7 (final R3) */ \
	r2 += r1;			/* R5 = R1 + R5 (final R1) */ \
\

#endif
#endif
\
/* The with-partial-normalization variant */ \
\

#define r4_x4cl_wpn_eight_reals_unfft_preload \
\
\

#define r4_x4cl_wpn_eight_reals_unfft(srcreg,srcinc,d1,d2,screg1,scoff1,screg2,scoff2) { \
	vec2f64 *sd0 = (vec2f64*)(srcreg); \
	vec2f64 *sd1 = (vec2f64*)(srcreg+d1); \
	vec2f64 *sd2 = (vec2f64*)(srcreg+d2); \
	vec2f64 *sd3 = (vec2f64*)(srcreg+d2+d1); \
	vec2f64 xma0,xma1,xma2,xma3,xma4,xma5,xma6,xma7; \
	/* tmp = srcinc+d2; */ \
	r4_x8r_wpn_unfft_mem(sd0[1],sd0[3],sd1[1],sd1[3],sd2[1],sd2[3],sd3[1],sd3[3],screg2+scoff2,screg1+scoff1+32,screg2+scoff2+48,screg1+scoff1,srcreg+srcinc+d2,d1,sd3[1],sd3[3]); \
/*	xstore	[srcreg+d3+16], xmm2	;; Save R4 */ \
/*	xstore	[srcreg+d3+48], xmm5	;; Save R8 */ \
	xma6 = sd3[0];	/* Load R4 */ \
	xma7 = sd3[2];	/* Load I4 */ \
	sd3[0] = xmm4;	/* Save R2 */ \
	sd3[2] = xmm3;	/* Save R6 */ \
	xma2 = sd1[0];	/* Load R2 */ \
	xma3 = sd1[2];	/* Load I2 */ \
	sd1[0] = xmm1;	/* Save R1 */ \
	sd1[1] = xmm7;	/* Save R3 */ \
	sd1[2] = xmm0;	/* Save R5 */ \
	sd1[3] = xmm6;	/* Save R7 */ \
	r4_x8r_wpn_unfft_partial_mem(xma0,xma1,xma2,xma3,xma4,xma5,xma6,xma7,sd0[0],sd0[2],sd2[0],sd2[2],screg2,screg1+32,screg2+48,screg1,srcreg+srcinc,d1,sd2[1],sd2[3]); \
/*	xstore	[srcreg+d2+16], xma2	;; Save R4 */ \
/*	xstore	[srcreg+d2+48], xma5	;; Save R8 */ \
	sd2[0] = xma4;	/* Save R2 */ \
	sd2[2] = xma3;	/* Save R6 */ \
	sd0[0] = xma1;	/* Save R1 */ \
	sd0[1] = xma7;	/* Save R3 */ \
	sd0[2] = xma0;	/* Save R5 */ \
	sd0[3] = xma6;	/* Save R7 */ \
	srcreg += srcinc; \
}


#define r4_x8r_wpn_unfft_mem(mem1,mem2,mem3,mem4,mem5,mem6,mem7,mem8,screg1,screg2,screg5,normreg,pre1,pre2,dst4,dst8) \
	xmm6 = xptr(screg5+16);	/* cosine/sine for w^5n */ \
	xmm0 = mem7;		/* R4 */ \
	xmm7 = mem8;		/* I4 */ \
	xmm1 = xmm6 * xmm0 + xmm7;	/* A4 = R4 * cosine/sine + I4 */ \
	xmm7 = xmm6 * xmm7 - xmm0;	/* B4 = I4 * cosine/sine - R4 */ \
\
	xmm2 = xptr(screg1+16);	/* cosine/sine for w^n */ \
	xmm0 = mem5;		/* R3 */ \
	xmm5 = mem6;		/* I3 */ \
	xmm4 = xmm2 * xmm0 + xmm5;	/* A3 = R3 * cosine/sine for w^n + I3 */ \
	xmm5 = xmm2 * xmm5 - xmm0;	/* B3 = I3 * cosine/sine - R3 */ \
\
	xmm6 = xptr(screg2+16);	/* cosine/sine for w^2n */ \
	xmm2 = mem3;		/* R2 */ \
	xmm3 = mem4;		/* I2 */ \
	xmm0 = xmm6 * xmm2 + xmm3;	/* A2 = R2 * cosine/sine + I2 */ \
	xmm3 = xmm6 * xmm3 - xmm2;	/* B2 = I2 * cosine/sine - R2 */ \
\
	xprefetchw(u8ptr(pre1)); \
\
	xmm6 = xptr(screg5+32);	/* Sine * normalization_inverse */ \
	xmm1 *= xmm6;		/* new R7 = A4 * sine (final R4) */ \
	xmm7 *= xmm6;		/* new R8 = B4 * sine (final I4) */ \
	xmm2 = xptr(screg1+32);	/* Sine * normalization_inverse */ \
	xmm4 *= xmm2;		/* new R5 = A3 * sine (final R3) */ \
	xmm5 *= xmm2;		/* new R6 = B3 * sine (final I3) */ \
	xmm6 = xptr(screg2+32);	/* Sine * normalization_inverse */ \
	xmm0 *= xmm6;		/* new R3 = A2 * sine (final R2) */ \
	xmm3 *= xmm6;		/* new R4 = B2 * sine (final I2) */ \
\
	xmm2 = mem2;		/* R2 */ \
	xmm2 *= xptr(normreg+16);	/* R2 * normalization_inverse */ \
\
	xprefetchw(u8ptr(pre1+pre2)); \
\
	xmm6 = xmm4;		/* Copy R5 */ \
	xmm4 -= xmm1;		/* new R6 = R5 - R7 */ \
	xmm1 += xmm6;		/* new R5 = R5 + R7 */ \
\
	xmm6 = xmm5;		/* Copy R6 */ \
	xmm5 -= xmm7;		/* new R8 = R6 - R8 */ \
	xmm7 += xmm6;		/* new R7 = R6 + R8 */ \
\
	xmm6 = xmm4;		/* Copy R6 */ \
	xmm4 += xmm5;		/* R6 = R6 + R8 */ \
	xmm5 -= xmm6;		/* R8 = R8 - R6 */ \
\
	xmm4 *= XMM_SQRTHALF2;	/* R6 = R6 * square root of 1/2 */ \
	xmm5 *= XMM_SQRTHALF2;	/* R8 = R8 * square root of 1/2 */ \
\
	xmm6 = xmm2;		/* Copy R2 */ \
	xmm2 -= xmm3;		/* R2 = R2 - R4 (new R4) */ \
	xmm3 += xmm6;		/* R4 = R2 + R4 (new R2) */ \
\
	xmm6 = xmm2;		/* Copy R4 */ \
	xmm2 -= xmm5;		/* R4 = R4 - R8 (final R8) */ \
	xmm5 += xmm6;		/* R8 = R4 + R8 (final R4) */ \
\
	xmm6 = xmm3;		/* Copy R2 */ \
	xmm3 -= xmm4;		/* R2 = R2 - R6 (final R6) */ \
	xmm4 += xmm6;		/* R6 = R2 + R6 (final R2) */ \
\
	xmm6 = mem1;		/* R1 */ \
	xmm6 *= xptr(normreg+16);	/* R1 * normalization_inverse */ \
\
	dst8 = xmm2;		/* Save R8 */ \
\
	xmm2 = xmm6;		/* Copy R1 */ \
	xmm6 -= xmm0;		/* R1 = R1 - R3 (new R3) */ \
	xmm0 += xmm2;		/* R3 = R1 + R3 (new R1) */ \
\
	dst4 = xmm5;		/* Save R4 */ \
\
	xmm2 = xmm6;		/* Copy R3 */ \
	xmm6 -= xmm7;		/* R3 = R3 - R7 (final R7) */ \
	xmm7 += xmm2;		/* R7 = R3 + R7 (final R3) */ \
\
	xmm5 = xmm0;		/* Copy R1 */ \
	xmm0 -= xmm1;		/* R1 = R1 - R5 (final R5) */ \
	xmm1 += xmm5;		/* R5 = R1 + R5 (final R1) */ \
\
\

#define r4_x8r_wpn_unfft_partial_mem(r1,r2,r3,r4,r5,r6,r7,r8,mem1,mem2,mem5,mem6,screg1,screg2,screg5,normreg,pre1,pre2,dst4,dst8) \
	r1 = xptr(screg5+16);	/* cosine/sine for w^5n */ \
	r2 = r1 * r7 + r8;	/* A4 = R4 * cosine/sine + I4 */ \
	r8 = r1 * r8 - r7;	/* B4 = I4 * cosine/sine - R4 */ \
\
	r1 = xptr(screg1+16);	/* cosine/sine for w^n */ \
	r7 = mem5;			/* R3 */ \
	r6 = mem6;			/* I3 */ \
	r5 = r1 * r7 + r6;	/* A3 = R3 * cosine/sine + I3 */ \
	r6 = r1 * r6 - r7;	/* B3 = I3 * cosine/sine - R3 */ \
\
	r7 = xptr(screg2+16);	/* cosine/sine for w^2n */ \
	r1 = r7 * r3 + r4;	/* A2 = R2 * cosine/sine + I2 */ \
	r4 = r7 * r4 - r3;	/* B2 = I2 * cosine/sine - R2 */ \
\
	xprefetchw(u8ptr(pre1)); \
\
	r7 = xptr(screg5+32);		/* Sine * normalization_inverse */ \
	r2 *= r7;			/* new R7 = A4 * sine (final R4) */ \
	r8 *= r7;			/* new R8 = B4 * sine (final I4) */ \
	r3 = xptr(screg1+32);		/* Sine * normalization_inverse */ \
	r5 *= r3;			/* new R5 = A3 * sine (final R3) */ \
	r6 *= r3;			/* new R6 = B3 * sine (final I3) */ \
	r7 = xptr(screg2+32);		/* Sine * normalization_inverse */ \
	r1 *= r7;			/* new R3 = A2 * sine (final R2) */ \
	r4 *= r7;			/* new R4 = B2 * sine (final I2) */ \
\
	r3 = mem2;		/* R2 */ \
	r3 *= xptr(normreg+16);	/* R2 * normalization_inverse */ \
\
	xprefetchw(u8ptr(pre1+pre2)); \
\
	r7 = r5;			/* Copy R5 */ \
	r5 -= r2;			/* new R6 = R5 - R7 */ \
	r2 += r7;			/* new R5 = R5 + R7 */ \
\
	r7 = r6;			/* Copy R6 */ \
	r6 -= r8;			/* new R8 = R6 - R8 */ \
	r8 += r7;			/* new R7 = R6 + R8 */ \
\
	r7 = r5;			/* Copy R6 */ \
	r5 += r6;			/* R6 = R6 + R8 */ \
	r6 -= r7;			/* R8 = R8 - R6 */ \
\
	r5 *= XMM_SQRTHALF2;	/* R6 = R6 * square root of 1/2 */ \
	r6 *= XMM_SQRTHALF2;	/* R8 = R8 * square root of 1/2 */ \
\
	r7 = r3;			/* Copy R2 */ \
	r3 -= r4;			/* R2 = R2 - R4 (new R4) */ \
	r4 += r7;			/* R4 = R2 + R4 (new R2) */ \
\
	r7 = r3;			/* Copy R4 */ \
	r3 -= r6;			/* R4 = R4 - R8 (final R8) */ \
	r6 += r7;			/* R8 = R4 + R8 (final R4) */ \
\
	r7 = mem1;		/* R1 */ \
	r7 *= xptr(normreg+16);	/* R1 * normalization_inverse */ \
\
	dst8 = r3;		/* Save R8 */ \
\
	r3 = r4;			/* Copy R2 */ \
	r4 -= r5;			/* R2 = R2 - R6 (final R6) */ \
	r5 += r3;			/* R6 = R2 + R6 (final R2) */ \
\
	r3 = r7;			/* Copy R1 */ \
	r7 -= r1;			/* R1 = R1 - R3 (new R3) */ \
	r1 += r3;			/* R3 = R1 + R3 (new R1) */ \
\
	dst4 = r6;		/* Save R4 */ \
\
	r3 = r7;			/* Copy R3 */ \
	r7 -= r8;			/* R3 = R3 - R7 (final R7) */ \
	r8 += r3;			/* R7 = R3 + R7 (final R3) */ \
\
	r6 = r1;			/* Copy R1 */ \
	r1 -= r2;			/* R1 = R1 - R5 (final R5) */ \
	r2 += r6;			/* R5 = R1 + R5 (final R1) */ \
\

#if 0 //DEF X86_64;
\
/* Theoretical best Core 2 timing is 48 clocks.  The empty add slot at clock 14 means */ \
/* this macro should run in 49 clocks.  Actual timing is ?? clocks. */ \
\

#define r4_x4cl_wpn_eight_reals_unfft_preload \
	xmm15 = XMM_SQRTHALF2; \
\
\

#define r4_x4cl_wpn_eight_reals_unfft(srcreg,srcinc,d1,d2,screg1,scoff1,screg2,scoff2) \
	xmm0 = xptr(srcreg+d2+d1+16);	/* R4 */ \
	xmm1 = xptr(screg2+scoff2+48+16);/* cosine/sine for w^5n */ \
	xmm2 = xmm0;		/* Copy R4 */ \
	xmm0 *= xmm1;		/* A4 = R4 * cosine/sine		; 1-5 */ \
\
	xmm3 = xptr(srcreg+d2+16);	/* R3 */ \
	xmm4 = xptr(screg2+scoff2+16);/* cosine/sine for w^n */ \
	xmm5 = xmm3;		/* Copy R3 */ \
	xmm3 *= xmm4;		/* A3 = R3 * cosine/sine for w^n	; 2-6 */ \
\
	xmm6 = xptr(srcreg+d2+d1+48);	/* I4 */ \
	xmm1 *= xmm6;		/* B4 = I4 * cosine/sine		; 3-7 */ \
\
	xmm7 = xptr(srcreg+d2+48);	/* I3 */ \
	xmm4 *= xmm7;		/* B3 = I3 * cosine/sine		; 4-8 */ \
\
	xmm8 = xptr(srcreg+d2+d1);	/*#2 R4 */ \
	xmm9 = xptr(screg2+48+16);	/*#2 cosine/sine for w^5n */ \
	xmm10 = xmm8;		/*#2 Copy R4 */ \
	xmm8 *= xmm9;		/*#2 A4 = R4 * cosine/sine		; 5-9 */ \
\
	xmm0 += xmm6;		/* A4 = A4 + I4				; 6-8	avail 6,11+ */ \
	xmm11 = xptr(srcreg+d2);	/*#2 R3 */ \
	xmm12 = xptr(screg2+16);	/*#2 cosine/sine for w^n */ \
	xmm13 = xmm11;		/*#2 Copy R3 */ \
	xmm11 *= xmm12;		/*#2 A3 = R3 * cosine/sine for w^n	; 6-10	avail 6,14 */ \
\
	xmm3 += xmm7;		/* A3 = A3 + I3				; 7-9	avail 6,7,14 */ \
	xmm14 = xptr(srcreg+d2+d1+32);/*#2 I4 */ \
	xmm9 *= xmm14;		/*#2 B4 = I4 * cosine/sine		; 7-11	avail 6,7 */ \
\
	xmm1 -= xmm2;		/* B4 = B4 - R4				; 8-10	avail 6,7,2 */ \
	xmm6 = xptr(srcreg+d2+32);	/*#2 I3 */ \
	xmm12 *= xmm6;		/*#2 B3 = I3 * cosine/sine		; 8-12	avail 7,2 */ \
\
	xmm4 -= xmm5;		/* B3 = B3 - R3				; 9-11	avail 7,2,5 */ \
	xmm7 = xptr(screg2+scoff2+48+32); /* sine * normalization_inverse */ \
	xmm0 *= xmm7;		/* new R7 = A4 * sine (final R4)	; 9-13	avail 2,5 */ \
\
	xmm8 += xmm14;		/*#2 A4 = A4 + I4			; 10-12	avail 2,5,14 */ \
	xmm2 = xptr(screg2+scoff2+32); /* sine * normalization_inverse */ \
	xmm3 *= xmm2;		/* new R5 = A3 * sine (final R3)	; 10-14	avail 5,14 */ \
\
	xmm11 += xmm6;		/*#2 A3 = A3 + I3			; 11-13	avail 5,14,6 */ \
	xmm1 *= xmm7;		/* new R8 = B4 * sine (final I4)	; 11-15	avail 5,14,6,7 */ \
\
	xmm9 -= xmm10;		/*#2 B4 = B4 - R4			; 12-14	avail 5,14,6,7,10 */ \
	xmm4 *= xmm2;		/* new R6 = B3 * sine (final I3)	; 12-16	avail 5,14,6,7,10,2 */ \
\
	xmm12 -= xmm13;		/*#2 B3 = B3 - R3			; 13-15	avail 5,14,6,7,10,2,13 */ \
	xmm5 = xptr(srcreg+d1+16);	/* R2 */ \
	xmm14 = xptr(screg1+scoff1+32+16);/* cosine/sine for w^2n */ \
	xmm6 = xmm5;		/* Copy R2 */ \
	xmm5 *= xmm14;		/* A2 = R2 * cosine/sine		; 13-17	avail 7,10,2,13 */ \
\
	xmm7 = xptr(srcreg+d1+48);	/* I2 */ \
	xmm14 *= xmm7;		/* B2 = I2 * cosine/sine		; 14-18	avail 10,2,13 */ \
\
	xprefetchw(u8ptr(srcreg+srcinc)); \
\
	xmm10 = xmm3;		/* Copy R5 */ \
	xmm3 -= xmm0;		/* new R6 = R5 - R7			; 15-17	avail 2,13 */ \
	xmm2 = xptr(screg2+48+32);	/*#2 sine * normalization_inverse */ \
	xmm8 *= xmm2;		/*#2 new R7 = A4 * sine (final R4)	; 15-19	avail 13 */ \
\
	xmm10 += xmm0;		/* new R5 = R5 + R7			; 16-18	avail 13,0 */ \
	xmm13 = xptr(screg2+32);	/*#2 sine * normalization_inverse */ \
	xmm11 *= xmm13;		/*#2 new R5 = A3 * sine (final R3)	; 16-20	avail 0 */ \
\
	xmm0 = xmm4;		/* Copy R6 */ \
	xmm4 -= xmm1;		/* new R8 = R6 - R8			; 17-19 */ \
	xmm9 *= xmm2;		/*#2 new R8 = B4 * sine (final I4)	; 17-21 avail 2 */ \
\
	xmm5 += xmm7;		/* A2 = A2 + I2				; 18-20	avail 2,7 */ \
	xmm12 *= xmm13;		/*#2 new R6 = B3 * sine (final I3)	; 18-22	avail 2,7,13 */ \
\
	xmm0 += xmm1;		/* new R7 = R6 + R8			; 19-21	avail 2,7,13,1 */ \
	xmm3 *= xmm15;		/* R6 = R6 * square root of 1/2		; 19-23 */ \
\
	xmm14 -= xmm6;		/* B2 = B2 - R2				; 20-22	avail 2,7,13,1,6 */ \
	xmm4 *= xmm15;		/* R8 = R8 * square root of 1/2		; 20-24 */ \
	xmm6 = xptr(srcreg+16);	/* R1					; 20 */ \
\
	xprefetchw(u8ptr(srcreg+srcinc+d1)); \
\
	xmm2 = xmm11;		/*#2 Copy R5 */ \
	xmm11 -= xmm8;		/*#2 new R6 = R5 - R7			; 21-23	avail 7,13,1 */ \
	xmm6 *= xptr(screg1+scoff1+16);/* R1 * normalization_inverse		; 21-25 */ \
\
	xmm2 += xmm8;		/*#2 new R5 = R5 + R7			; 22-24	avail 7,13,1,8 */ \
	xmm7 = xptr(screg1+scoff1+32+32); /* sine * normalization_inverse */ \
	xmm5 *= xmm7;		/* new R3 = A2 * sine (final R2)	; 22-26	avail 13,1,8 */ \
\
	xmm13 = xmm12;		/*#2 Copy R6 */ \
	xmm12 -= xmm9;		/*#2 new R8 = R6 - R8			; 23-25	avail 1,8 */ \
	xmm14 *= xmm7;		/* new R4 = B2 * sine (final I2)	; 23-27	avail 1,8,7 */ \
	xmm7 = xptr(srcreg+48);	/* R2 */ \
\
	xmm13 += xmm9;		/*#2 new R7 = R6 + R8			; 24-26	avail 1,8,9 */ \
	xmm7 *= xptr(screg1+scoff1+16);/* R2 * normalization_inverse		; 24-28 */ \
\
	xmm1 = xmm3;		/* Copy R6 */ \
	xmm3 += xmm4;		/* R6 = R6 + R8				; 25-27	avail 8,7,9 */ \
\
	xmm4 -= xmm1;		/* R8 = R8 - R6				; 26-28	avail 8,9,1 */ \
\
	xprefetchw(u8ptr(srcreg+srcinc+d2));	/*#2 */ \
\
	xmm8 = xmm6;		/* Copy R1 */ \
	xmm6 -= xmm5;		/* R1 = R1 - R3 (new R3)		; 27-29	avail 9,1 */ \
	xmm8 += xmm5;		/* R3 = R1 + R3 (new R1)		; 28-30	avail 9,1,5 */ \
\
	xmm9 = xmm7;		/* Copy R2 */ \
	xmm7 -= xmm14;		/* R2 = R2 - R4 (new R4)		; 29-31	avail 1,5 */ \
\
	xmm14 += xmm9;		/* R4 = R2 + R4 (new R2)		; 30-32	avail 1,5,9 */ \
	xmm1 = xptr(srcreg+d1);	/*#2 R2 */ \
	xmm5 = xptr(screg1+32+16);	/*#2 cosine/sine for w^2n */ \
	xmm1 *= xmm5;		/*#2 A2 = R2 * cosine/sine		; 30-34	avail 9 */ \
\
	xmm9 = xmm6;		/* Copy R3 */ \
	xmm6 -= xmm0;		/* R3 = R3 - R7 (final R7)		; 31-33	avail none storable 6 */ \
	xptr(srcreg+d1+48) = xmm6;	/* Save R7				; 34	avail 6 */ \
	xmm6 = xptr(srcreg+d1+32);	/*#2 I2 */ \
	xmm5 *= xmm6;		/*#2 B2 = I2 * cosine/sine		; 31-35	avail none */ \
\
	xmm0 += xmm9;		/* R7 = R3 + R7 (final R3)		; 32-34	avail 9 storable 0 */ \
\
	xmm9 = xmm8;		/* Copy R1 */ \
	xmm8 -= xmm10;		/* R1 = R1 - R5 (final R5)		; 33-35	avail none storable 0,8 */ \
	xmm10 += xmm9;		/* R5 = R1 + R5 (final R1)		; 34-36	avail 9 storable 0,8,10 */ \
\
	xmm1 += xmm6;		/*#2 A2 = A2 + I2			; 35-37	avail 9,6 storable 0,8,10 */ \
	xmm11 *= xmm15;		/*#2 R6 = R6 * square root of 1/2	; 35-39 */ \
	xptr(srcreg+d1+16) = xmm0;	/* Save R3				; 35	avail 9,6,0 storable 8,10 */ \
\
	xmm5 -= xptr(srcreg+d1);	/*#2 B2 = B2 - R2			; 36-38 */ \
	xmm12 *= xmm15;		/*#2 R8 = R8 * square root of 1/2	; 36-40 */ \
	xptr(srcreg+d1+32) = xmm8;	/* Save R5				; 36	avail 9,6,0,8 storable 10 */ \
\
	xmm9 = xmm7;		/* Copy R4 */ \
	xmm7 -= xmm4;		/* R4 = R4 - R8 (final R8)		; 37-39	avail 6,0,8 storable 10,7 */ \
	xptr(srcreg+d1) = xmm10;	/* Save R1				; 37	avail 6,0,8,10 storable 7 */ \
	xmm10 = xptr(srcreg);		/*#2 R1					; 37 */ \
	xmm10 *= xptr(screg1+16);	/*#2 R1 * normalization_inverse		; 37-41 */ \
\
	xprefetchw(u8ptr(srcreg+srcinc+d2+d1));/*#2 */ \
\
	xmm4 += xmm9;		/* R8 = R4 + R8 (final R4)		; 38-40	avail 6,0,8,9 storable 7,4 */ \
	xmm6 = xptr(screg1+32+32);	/*#2 sine * normalization_inverse */ \
	xmm1 *= xmm6;		/*#2 new R3 = A2 * sine (final R2)	; 38-42	avail 0,8,9 storable 7,4 */ \
\
	xmm0 = xmm14;		/* Copy R2 */ \
	xmm14 -= xmm3;		/* R2 = R2 - R6 (final R6)		; 39-41	avail 8,9 storable 7,4,14 */ \
	xmm5 *= xmm6;		/*#2 new R4 = B2 * sine (final I2)	; 39-43	avail 8,9,6 storable 7,4,14 */ \
	xmm6 = xptr(srcreg+32);	/*#2 R2					; 39 */ \
\
	xmm3 += xmm0;		/* R6 = R2 + R6 (final R2)		; 40-42	avail 8,9,0 storable 7,4,14,3 */ \
	xmm6 *= xptr(screg1+16);	/*#2 R2 * normalization_inverse		; 40-44 */ \
	xptr(srcreg+d2+d1+48) = xmm7;	/* Save R8				; 40	avail 8,9,0,7 storable 4,14,3 */ \
\
	xmm8 = xmm11;		/*#2 Copy R6 */ \
	xmm11 += xmm12;		/*#2 R6 = R6 + R8			; 41-43	avail 9,0,7 storable 4,14,3 */ \
	xptr(srcreg+d2+d1+16) = xmm4;	/* Save R4				; 41	avail 9,0,7,4 storable 14,3 */ \
\
	xmm12 -= xmm8;		/*#2 R8 = R8 - R6			; 42-44	avail 9,0,7,4,8 storable 14,3 */ \
	xptr(srcreg+d2+d1+32) = xmm14;/* Save R6				; 42	avail 9,0,7,4,8,14 storable 3 */ \
\
	xmm9 = xmm10;		/*#2 Copy R1 */ \
	xmm10 -= xmm1;		/*#2 R1 = R1 - R3 (new R3)		; 43-45	avail 0,7,4,8,14 storable 3 */ \
	xptr(srcreg+d2+d1) = xmm3;	/* Save R2				; 43	avail 0,7,4,8,14,3 */ \
\
	xmm1 += xmm9;		/*#2 R3 = R1 + R3 (new R1)		; 44-46	avail 0,7,4,8,14,3,9 */ \
\
	xmm0 = xmm6;		/*#2 Copy R2 */ \
	xmm6 -= xmm5;		/*#2 R2 = R2 - R4 (new R4)		; 45-47	avail 7,4,8,14,3,9 */ \
\
	xmm5 += xmm0;		/*#2 R4 = R2 + R4 (new R2)		; 46-48	avail 7,4,8,14,3,9,0 */ \
\
	xmm7 = xmm10;		/*#2 Copy R3 */ \
	xmm10 -= xmm13;		/*#2 R3 = R3 - R7 (final R7)		; 47-49	avail 4,8,14,3,9,0 storable 10 */ \
	xmm13 += xmm7;		/*#2 R7 = R3 + R7 (final R3)		; 48-50	avail 4,8,14,3,9,0,7 storable 10,13 */ \
\
	xmm4 = xmm1;		/*#2 Copy R1 */ \
	xmm1 -= xmm2;		/*#2 R1 = R1 - R5 (final R5)		; 49-51	avail 8,14,3,9,0,7 storable 10,13,1 */ \
\
	xmm2 += xmm4;		/*#2 R5 = R1 + R5 (final R1)		; 50-52	avail 8,14,3,9,0,7,4 storable 10,13,1,2 */ \
	xptr(srcreg+48) = xmm10;	/*#2 Save R7				; 50	avail 8,14,3,9,0,7,4,10 storable 13,1,2 */ \
\
	xmm8 = xmm6;		/*#2 Copy R4 */ \
	xmm6 -= xmm12;		/*#2 R4 = R4 - R8 (final R8)		; 51-53	avail 14,3,9,0,7,4,10 storable 13,1,2,6 */ \
	xptr(srcreg+16) = xmm13;	/*#2 Save R3				; 51	avail 14,3,9,0,7,4,10,13 storable 1,2,6 */ \
\
	xmm12 += xmm8;		/*#2 R8 = R4 + R8 (final R4)		; 52-54	avail 14,3,9,0,7,4,10,13,8 storable 1,2,6,12 */ \
	xptr(srcreg+32) = xmm1;	/*#2 Save R5				; 52	avail 14,3,9,0,7,4,10,13,8,1 storable 2,6,12 */ \
\
	xmm8 = xmm5;		/*#2 Copy R2 */ \
	xmm5 -= xmm11;		/*#2 R2 = R2 - R6 (final R6)		; 53-55	avail 14,3,9,0,7,4,10,13,1 storable 2,6,12,5 */ \
	xptr(srcreg) = xmm2;		/*#2 Save R1				; 53	avail 14,3,9,0,7,4,10,13,1,2 storable 6,12,5 */ \
\
	xmm11 += xmm8;		/*#2 R6 = R2 + R6 (final R2)		; 54-56	avail 14,3,9,0,7,4,10,13,1,2,8 storable 6,12,5,11 */ \
	xptr(srcreg+d2+48) = xmm6;	/*#2 Save R8				; 54 */ \
	xptr(srcreg+d2+16) = xmm12;	/*#2 Save R4				; 55 */ \
	xptr(srcreg+d2+32) = xmm5;	/*#2 Save R6				; 56 */ \
	xptr(srcreg+d2) = xmm11;	/*#2 Save R2				; 57 */ \
\
	srcreg += srcinc; \
\

#endif
\
/* */ \
/* ************************************* eight-reals-four-complex-djbfft variants ****************************************** */ \
/* */ \
/* Because of the way we store FFT data elements in cache lines during the pass 1 of the forward FFT, a cache line */ \
/* will contain data needing an eight-reals-fft as well as data needing a four-complex-fft. */ \
/* */ \
\
/* Macro to do an eight_reals_fft and a four_complex_djbfft. */ \
/* The eight-reals macro is very similar to the eight-reals-first-FFT. */ \
/* Both do 2 and 3/4 levels of the FFT.  One difference this */ \
/* macro takes two screg pointers.  One points to w^n and w^5n, the other */ \
/* points to w^2n which is also the screg used by the four-complex-FFT. */ \
\

#define r4_x4cl_eight_reals_four_complex_djbfft_preload \
	r4_x8r_fft_mem_preload; \
;	/*r4_x4c_djbfft_partial_mem_preload -- must be compatible with r4_x8r_fft_mem_preload!! */ \
\
\

#define r4_x4cl_eight_reals_four_complex_djbfft(srcreg,srcinc,d1,d2,screg1,screg2) { \
	vec2f64 *sd0 = (vec2f64*)(srcreg); \
	vec2f64 *sd1 = (vec2f64*)(srcreg+d1); \
	vec2f64 *sd2 = (vec2f64*)(srcreg+d2); \
	vec2f64 *sd3 = (vec2f64*)(srcreg+d2+d1); \
	r4_x8r_fft_mem(sd0[0],sd1[0],sd2[0],sd3[0],sd0[1],sd1[1],sd2[1],sd3[1],screg2,screg1,screg2+32,srcreg+srcinc,d1,sd0[0],sd0[1]); \
/*	xstore	[srcreg], xmm0		;; Save real value #1 */ \
/*	xstore	[srcreg+16], xmm0	;; Save real value #2 */ \
	xmm0 = sd0[2];	/* R1 */ \
	xmm1 = sd0[3];	/* R5 */ \
	sd0[2] = xmm7;	/* Save R2 */ \
	sd0[3] = xmm5;	/* Save I2 */ \
	xmm7 = sd1[2];	/* R2 */ \
	xmm5 = sd1[3];	/* R6 */ \
	sd1[0] = xmm4;	/* Save R3 */ \
	sd1[1] = xmm3;	/* Save I3 */ \
	sd1[2] = xmm6;	/* Save R4 */ \
	sd1[3] = xmm2;	/* Save I4 */ \
	r4_x4c_djbfft_partial_mem(xmm0,xmm7,xmm3,xmm2,xmm1,xmm5,xmm4,xmm6,sd2[2],sd3[2],sd2[3],sd3[3],screg1,0,srcreg+srcinc+d2,d1,sd2[0],sd2[1]); \
/*	xstore	[srcreg+d2], xmm0	;; Save real value #1 */ \
/*	xstore	[srcreg+d2+16], xmm0	;; Save real value #2 */ \
	sd2[2] = xmm6;	/* Save R2 */ \
	sd2[3] = xmm4;	/* Save I2 */ \
	sd3[0] = xmm2;	/* Save R3 */ \
	sd3[1] = xmm7;	/* Save I3 */ \
	sd3[2] = xmm0;	/* Save R4 */ \
	sd3[3] = xmm1;	/* Save I4 */ \
	srcreg += srcinc; \
}

/* Macro to do an eight_reals_fft and a four_complex_djbfft. */ \
/* This macro is used in the last levels of an r4 FFT pass 1 (no swizzling) */ \
\

#define r4_sg4cl_eight_reals_four_complex_djbfft_preload \
	r4_x8r_fft_mem_preload; \
;	/*r4_x4c_djbfft_partial_mem_preload -- must be compatible with r4_x8r_fft_mem_preload!! */ \
\
\

#define r4_sg4cl_eight_reals_four_complex_djbfft(srcreg,srcinc,d1,d2,dstreg,dstinc,e1,e2,screg1,screg2) { \
	vec2f64 *sd0 = (vec2f64*)(srcreg); \
	vec2f64 *sd1 = (vec2f64*)(srcreg+d1); \
	vec2f64 *sd2 = (vec2f64*)(srcreg+d2); \
	vec2f64 *sd3 = (vec2f64*)(srcreg+d2+d1); \
	vec2f64 *de0 = (vec2f64*)(dstreg); \
	vec2f64 *de1 = (vec2f64*)(dstreg+e1); \
	vec2f64 *de2 = (vec2f64*)(dstreg+e2); \
	vec2f64 *de3 = (vec2f64*)(dstreg+e2+e1); \
	r4_x8r_fft_mem(sd0[0],sd1[0],sd2[0],sd3[0],sd0[1],sd1[1],sd2[1],sd3[1], screg2,screg1,screg2+32,dstreg+dstinc,e1, de0[0],de0[2]); \
	shuffle_store_partial(de0[0], de0[1], xmm0, xmm7);				/* Save real value #1,R2 */ \
	shuffle_store_partial(de0[2], de0[3], xmm1, xmm5);			/* Save real value #2,I2 */ \
	xmm0 = sd0[2];	/* R1 */ \
	xmm5 = sd0[3];	/* R5 */ \
	shuffle_store_with_temp(de1[0], de1[1], xmm4, xmm6, xmm1);		/* Save R3,R4 */ \
	shuffle_store_with_temp(de1[2], de1[3], xmm3, xmm2, xmm1);	/* Save I3,I4 */ \
	xmm7 = sd1[2];	/* R2 */ \
	xmm1 = sd1[3];	/* R6 */ \
	xprefetch(u8ptr(srcreg+srcinc+d2)); \
	r4_x4c_djbfft_partial_mem(xmm0,xmm7,xmm3,xmm2,xmm5,xmm1,xmm6,xmm4, sd2[2],sd3[2],sd2[3],sd3[3],screg1,0,dstreg+dstinc+e2,e1, de2[0],de2[2]); \
	xprefetch(u8ptr(srcreg+srcinc+d2+d1)); \
	srcreg += srcinc; \
	shuffle_store_partial(de2[0], de2[1], xmm3, xmm4);			/* Save R1,R2 */ \
	shuffle_store_partial(de2[2], de2[3], xmm1, xmm6);		/* Save I1,I2 */ \
	shuffle_store_with_temp(de3[0], de3[1], xmm2, xmm0, xmm3);	/* Save R3,R4 */ \
	shuffle_store_with_temp(de3[2], de3[3], xmm7, xmm5, xmm3);	/* Save I3,I4 */ \
	dstreg += dstinc; \
}


/* Macro to do an eight_reals_fft and a four_complex_fft4 at the end of pass 1. */ \
/* We could do a four_complex_djbfft since the 4 twiddle factors belong to the */ \
/* first grouping.  However, this will only work if the corresponding unfft does */ \
/* a djbunfft.  This macro is very similar to r4_x4cl_eight_reals_four_complex_djbfft. */ \

/* This macro is used in the last levels of an r4delay pass 1 (no swizzling) */

#if 0 //IFDEF UNUSED
#define r4_g4cl_eight_reals_four_complex_fft4_preload \
	r4_x8r_fft_mem_preload; \
;	/*r4_x4c_djbfft_partial_mem_preload -- must be compatible with r4_x8r_fft_mem_preload!! */ \
\

#define r4_g4cl_eight_reals_four_complex_fft4(srcreg,srcinc,d1,d2,dstreg,dstinc,e1,e2,screg1,screg2) \
	r4_x8r_fft_mem(xptr(srcreg),xptr(srcreg+d1),xptr(srcreg+d2),xptr(srcreg+d2+d1),xptr(srcreg+16),xptr(srcreg+d1+16),xptr(srcreg+d2+16),xptr(srcreg+d2+d1+16),screg2,screg1+32,screg2+32,dstreg+dstinc,e1,xptr(dstreg),xptr(dstreg+16)); \
/*	xstore	[dstreg], xmm0		;; Save real value #1 */ \
/*	xstore	[dstreg+16], xmm0	;; Save real value #2 */ \
	xmm0 = xptr(srcreg+32);	/* R1 */ \
	xmm1 = xptr(srcreg+48);	/* R5 */ \
	xptr(dstreg+32) = xmm7;	/* Save R2 */ \
	xptr(dstreg+48) = xmm5;	/* Save I2 */ \
	xmm7 = xptr(srcreg+d1+32);	/* R2 */ \
	xmm5 = xptr(srcreg+d1+48);	/* R6 */ \
	xptr(dstreg+e1) = xmm4;	/* Save R3 */ \
	xptr(dstreg+e1+16) = xmm3;	/* Save I3 */ \
	xptr(dstreg+e1+32) = xmm6;	/* Save R4 */ \
	xptr(dstreg+e1+48) = xmm2;	/* Save I4 */ \
	r4_x4c_fft_partial_mem(xmm0,xmm7,xmm3,xmm2,xmm1,xmm5,xmm4,xmm6,xptr(srcreg+d2+32),xptr(srcreg+d2+d1+32),xptr(srcreg+d2+48),xptr(srcreg+d2+d1+48),screg1,32,dstreg+dstinc+e2,e1,xptr(dstreg+e2),xptr(dstreg+e2+16)); \
	srcreg += srcinc; \
/*	xstore	[dstreg+e2], xmm0	;; Save R1 */ \
/*	xstore	[dstreg+e2+16], xmm0	;; Save I1 */ \
	xptr(dstreg+e2+32) = xmm6;	/* Save R2 */ \
	xptr(dstreg+e2+48) = xmm4;	/* Save I2 */ \
	xptr(dstreg+e2+e1) = xmm2;	/* Save R3 */ \
	xptr(dstreg+e2+e1+16) = xmm7;	/* Save I3 */ \
	xptr(dstreg+e2+e1+32) = xmm0;	/* Save R4 */ \
	xptr(dstreg+e2+e1+48) = xmm1;	/* Save I4 */ \
	dstreg += dstinc; \

#endif
\
/* Used in the last levels of an r4delay pass 1 (swizzling) */ \

#define r4_sg4cl_eight_reals_four_complex_fft4_preload \
	r4_x8r_fft_mem_preload; \
;	/*r4_x4c_djbfft_mem_preload -- must be compatible with r4_x8r_fft_mem_preload!! */ \
\

#define r4_sg4cl_eight_reals_four_complex_fft4(srcreg,srcinc,d1,d2,dstreg,dstinc,e1,e2,screg1,screg2) { \
	vec2f64 *sd0 = (vec2f64*)(srcreg); \
	vec2f64 *sd1 = (vec2f64*)(srcreg+d1); \
	vec2f64 *sd2 = (vec2f64*)(srcreg+d2); \
	vec2f64 *sd3 = (vec2f64*)(srcreg+d2+d1); \
	vec2f64 *de0 = (vec2f64*)(dstreg); \
	vec2f64 *de1 = (vec2f64*)(dstreg+e1); \
	vec2f64 *de2 = (vec2f64*)(dstreg+e2); \
	vec2f64 *de3 = (vec2f64*)(dstreg+e2+e1); \
	r4_x8r_fft_mem(sd0[0],sd1[0],sd2[0],sd3[0],sd0[1],sd1[1],sd2[1],sd3[1], screg2,screg1+32,screg2+32,dstreg+dstinc,e1, de0[0],de0[2]); \
	shuffle_store_partial(de0[0], de0[1], xmm0, xmm7);				/* Save real value #1,R2 */ \
	shuffle_store_partial(de0[2], de0[3], xmm1, xmm5);			/* Save real value #2,I2 */ \
	shuffle_store_with_temp(de1[0], de1[1], xmm4, xmm6, xmm7);		/* Save R3,R4 */ \
	shuffle_store_with_temp(de1[2], de1[3], xmm3, xmm2, xmm7);	/* Save I3,I4 */ \
	xprefetch(u8ptr(srcreg+srcinc+d2)); \
	r4_x4c_fft_mem(sd0[2],sd1[2],sd2[2],sd3[2],sd0[3],sd1[3],sd2[3],sd3[3], screg1+32,screg1+64,screg1+96,dstreg+dstinc+e2,e1, de2[0],de2[2]); \
	xprefetch(u8ptr(srcreg+srcinc+d2+d1)); \
	srcreg += srcinc; \
	shuffle_store_partial(de2[0], de2[1], xmm2, xmm7);			/* Save R1,R2 */ \
	shuffle_store_partial(de2[2], de2[3], xmm5, xmm6);		/* Save I1,I2 */ \
	shuffle_store_with_temp(de3[0], de3[1], xmm3, xmm0, xmm7);	/* Save R3,R4 */ \
	shuffle_store_with_temp(de3[2], de3[3], xmm1, xmm4, xmm7);	/* Save I3,I4 */ \
	dstreg += dstinc; \
}

/* The with-partial-normalization variant */ \
\

#define r4_x4cl_wpn_eight_reals_four_complex_djbfft_preload \
	r4_x8r_wpn_fft_mem_preload; \
;	/*r4_x4c_wpn_djbfft_partial_mem -- assumed to be a subset of r4_x8r_wpn_fft_mem_preload */ \
\
\

#define r4_x4cl_wpn_eight_reals_four_complex_djbfft(srcreg,srcinc,d1,d2,screg1,screg2) { \
	vec2f64 *sd0 = (vec2f64*)(srcreg); \
	vec2f64 *sd1 = (vec2f64*)(srcreg+d1); \
	vec2f64 *sd2 = (vec2f64*)(srcreg+d2); \
	vec2f64 *sd3 = (vec2f64*)(srcreg+d2+d1); \
	vec2f64 *sc1 = (vec2f64*)(screg1); \
	vec2f64 *sc2 = (vec2f64*)(screg2); \
	vec2f64 sc1_0 = sc1[0]; \
	vec2f64 sc1_2 = sc1[2]; \
	vec2f64 sc1_3 = sc1[3]; \
	r4_x8r_wpn_fft_mem4(sd0[0],sd1[0],sd2[0],sd3[0],sd0[1],sd1[1],sd2[1],sd3[1], sc2[0],sc2[1],sc1_2,sc1_3,sc2[3],sc2[4],sc1_0, srcreg+srcinc,d1, sd0[0],sd0[1]); \
/*	xstore	[srcreg], xmm0		;; Save real value #1 */ \
/*	xstore	[srcreg+16], xmm0	;; Save real value #2 */ \
	xmm0 = sd0[2];	/* R1 */ \
	xmm1 = sd0[3];	/* R5 */ \
	sd0[2] = xmm7;	/* Save R2 */ \
	sd0[3] = xmm5;	/* Save I2 */ \
	xmm7 = sd1[2];	/* R2 */ \
	xmm5 = sd1[3];	/* R6 */ \
	sd1[0] = xmm4;	/* Save R3 */ \
	sd1[1] = xmm3;	/* Save I3 */ \
	sd1[2] = xmm6;	/* Save R4 */ \
	sd1[3] = xmm2;	/* Save I4 */ \
	r4_x4c_wpn_djbfft_partial_mem4(xmm0,xmm7,xmm3,xmm2,xmm1,xmm5,xmm4,xmm6, sd2[2],sd3[2],sd2[3],sd3[3], sc1_2,sc1_3,sc1[5],sc1[6],sc1_0, srcreg+srcinc+d2,d1, sd2[0],sd2[1]); \
/*	xstore	[srcreg+d2], xmm0	;; Save real value #1 */ \
/*	xstore	[srcreg+d2+16], xmm0	;; Save real value #2 */ \
	sd2[2] = xmm6;	/* Save R2 */ \
	sd2[3] = xmm4;	/* Save I2 */ \
	sd3[0] = xmm2;	/* Save R3 */ \
	sd3[1] = xmm7;	/* Save I3 */ \
	sd3[2] = xmm0;	/* Save R4 */ \
	sd3[3] = xmm1;	/* Save I4 */ \
	srcreg += srcinc; \
}


#define r4_x8r_wpn_fft_mem_preload

#define r4_x8r_wpn_fft_mem4(mem1,mem2,mem3,mem4,mem5,mem6,mem7,mem8,SC0,SC1,SC2,SC3,SC4,SC5,normreg,pre1,pre2,dst1,dst2) \
	xmm3 = mem4;		/* R4 */ \
	xmm7 = mem8;		/* R8 */ \
	xmm0 = xmm3; \
	xmm3 -= xmm7;		/* new R8 = R4 - R8 */ \
	xmm7 += xmm0;		/* new R4 = R4 + R8 */ \
\
	xmm1 = mem2;		/* R2 */ \
	xmm5 = mem6;		/* R6 */ \
	xmm4 = xmm1; \
	xmm1 -= xmm5;		/* new R6 = R2 - R6 */ \
	xmm5 += xmm4;		/* new R2 = R2 + R6 */ \
\
	xmm3 *= XMM_SQRTHALF2;	/* R8 = R8 * square root */ \
	xmm1 *= XMM_SQRTHALF2;	/* R6 = R6 * square root */ \
\
	xmm0 = mem1;		/* R1 */ \
	xmm4 = mem5;		/* R5 */ \
	xmm2 = xmm0; \
	xmm0 -= xmm4;		/* new R5 = R1 - R5 */ \
	xmm4 += xmm2;		/* new R1 = R1 + R5 */ \
\
	xmm2 = xmm5;		/* Copy R2 */ \
	xmm5 -= xmm7;		/* R2 = R2 - R4 (final I2) */ \
	xmm7 += xmm2;		/* R4 = R2 + R4 (final I1, a.k.a 2nd real result) */ \
\
	xmm6 = xmm1;		/* Copy R6 */ \
	xmm1 -= xmm3;		/* R6 = R6 - R8 (Real part) */ \
	xmm3 += xmm6;		/* R8 = R6 + R8 (Imaginary part) */ \
\
	xmm6 = mem7;		/* R7 */ \
	xmm6 += mem3;		/* new R3 = R3 + R7 */ \
\
	xmm2 = xmm4;		/* Copy R1 */ \
	xmm4 -= xmm6;		/* R1 = R1 - R3 (final R2) */ \
	xmm6 += xmm2;		/* R3 = R1 + R3 (final R1) */ \
\
	xmm2 = mem3;		/* R3 */ \
	xmm2 -= mem7;		/* new R7 = R3 - R7 */ \
\
	xprefetchw(u8ptr(pre1)); \
\
	xmm6 *= normreg;	/* R1 * normalization value */ \
	dst1 = xmm6;		/* Save R1 */ \
	xmm7 *= normreg;	/* I1 * normalization value */ \
	dst2 = xmm7;		/* Save I1 */ \
\
	xmm7 = xmm0;		/* Copy R5 */ \
	xmm0 -= xmm1;		/* R5 = R5 - R6 (final R4) */ \
	xmm1 += xmm7;		/* R6 = R5 + R6 (final R3) */ \
\
	xmm6 = xmm2;		/* Copy R7 */ \
	xmm2 -= xmm3;		/* R7 = R7 - R8 (final I4) */ \
	xmm3 += xmm6;		/* R8 = R7 + R8 (final I3) */ \
\
	xprefetchw(u8ptr(pre1+pre2)); \
\
	xmm6 = SC3;			/* cosine/sine for w^2n */ \
	xmm7 = xmm4 * xmm6 - xmm5;	/* A2 = R2 * cosine/sine - I2 */ \
	xmm5 = xmm5 * xmm6 + xmm4;	/* B2 = I2 * cosine/sine + R2 */ \
\
	xmm4 = SC5;			/* cosine/sine for w^5n */ \
	xmm6 = xmm0 * xmm4 - xmm2;	/* A4 = R4 * cosine/sine - I4 */ \
	xmm2 = xmm2 * xmm4 + xmm0;	/* B4 = I4 * cosine/sine + R4 */ \
\
	xmm0 = SC1;			/* cosine/sine for w^n */ \
	xmm4 = xmm1 * xmm0 - xmm3;	/* A3 = R3 * cosine/sine - I3 */ \
	xmm3 = xmm3 * xmm0 + xmm1;	/* B3 = I3 * cosine/sine + R3 */ \
\
	xmm0 = SC2; \
	xmm7 *= xmm0;		/* A2 = A2 * sine (final R2) */ \
	xmm5 *= xmm0;		/* B2 = B2 * sine (final I2) */ \
	xmm1 = SC4; \
	xmm6 *= xmm1;		/* A4 = A4 * sine (final R4) */ \
	xmm2 *= xmm1;		/* B4 = B4 * sine (final I4) */ \
	xmm0 = SC0; \
	xmm4 *= xmm0;		/* A3 = A3 * sine (final R3) */ \
	xmm3 *= xmm0;		/* B3 = B3 * sine (final I3) */ \
\
\
/* 32-bit AMD K8 version */ \

#if 0 // (@INSTR(,%xarch,<K8>) NE 0);
\

#define r4_x8r_wpn_fft_mem(mem1,mem2,mem3,mem4,mem5,mem6,mem7,mem8,screg1,screg2,screg5,normreg,pre1,pre2,dst1,dst2) \
	xmm3 = mem4;		/* R4 */ \
	xmm7 = mem8;		/* R8 */ \
	xmm3 -= xmm7;		/* new R8 = R4 - R8 */ \
	xmm7 += mem4;		/* new R4 = R4 + R8 */ \
\
	xmm1 = mem2;		/* R2 */ \
	xmm5 = mem6;		/* R6 */ \
	xmm1 -= xmm5;		/* new R6 = R2 - R6 */ \
	xmm5 += mem2;		/* new R2 = R2 + R6 */ \
\
	xmm3 *= XMM_SQRTHALF2;	/* R8 = R8 * square root */ \
	xmm1 *= XMM_SQRTHALF2;	/* R6 = R6 * square root */ \
\
	xmm0 = mem1;		/* R1 */ \
	xmm4 = mem5;		/* R5 */ \
	xmm0 -= xmm4;		/* new R5 = R1 - R5 */ \
	xmm4 += mem1;		/* new R1 = R1 + R5 */ \
\
	xmm5 -= xmm7;		/* R2 = R2 - R4 (final I2) */ \
	xmm7 += xmm7; /* multwo */; \
	xmm7 += xmm5;		/* R4 = R2 + R4 (final I1, a.k.a 2nd real result) */ \
\
	xmm1 -= xmm3;		/* R6 = R6 - R8 (Real part) */ \
	xmm3 += xmm3; /* multwo */; \
	xmm3 += xmm1;		/* R8 = R6 + R8 (Imaginary part) */ \
\
	xmm6 = mem7;		/* R7 */ \
	xmm6 += mem3;		/* new R3 = R3 + R7 */ \
\
	xmm4 -= xmm6;		/* R1 = R1 - R3 (final R2) */ \
	xmm6 += xmm6; /* multwo */; \
	xmm6 += xmm4;		/* R3 = R1 + R3 (final R1) */ \
\
	xmm2 = mem3;		/* R3 */ \
	xmm2 -= mem7;		/* new R7 = R3 - R7 */ \
\
	xprefetchw(u8ptr(pre1)); \
\
	xmm6 *= xptr(normreg);		/* R1 * normalization value */ \
	dst1 = xmm6;		/* Save R1 */ \
	xmm7 *= xptr(normreg);		/* I1 * normalization value */ \
	dst2 = xmm7;		/* Save I1 */ \
\
	xmm0 -= xmm1;		/* R5 = R5 - R6 (final R4) */ \
	xmm1 += xmm1; /* multwo */; \
	xmm1 += xmm0;		/* R6 = R5 + R6 (final R3) */ \
\
	xmm2 -= xmm3;		/* R7 = R7 - R8 (final I4) */ \
	xmm3 += xmm3; /* multwo */; \
	xmm3 += xmm2;		/* R8 = R7 + R8 (final I3) */ \
\
	xprefetchw(u8ptr(pre1+pre2)); \
\
	xmm7 = xptr(screg2+16);	/* cosine/sine for w^2n */ \
	xmm7 *= xmm4;		/* A2 = R2 * cosine/sine */ \
	xmm7 -= xmm5;		/* A2 = A2 - I2 */ \
	xmm5 *= xptr(screg2+16);	/* B2 = I2 * cosine/sine */ \
	xmm5 += xmm4;		/* B2 = B2 + R2 */ \
\
	xmm6 = xptr(screg5+16);	/* cosine/sine for w^5n */ \
	xmm6 *= xmm0;		/* A4 = R4 * cosine/sine */ \
	xmm6 -= xmm2;		/* A4 = A4 - I4 */ \
	xmm2 *= xptr(screg5+16);	/* B4 = I4 * cosine/sine */ \
	xmm2 += xmm0;		/* B4 = B4 + R4 */ \
\
	xmm4 = xptr(screg1+16);	/* cosine/sine for w^n */ \
	xmm4 *= xmm1;		/* A3 = R3 * cosine/sine */ \
	xmm4 -= xmm3;		/* A3 = A3 - I3 */ \
	xmm3 *= xptr(screg1+16);	/* B3 = I3 * cosine/sine */ \
	xmm3 += xmm1;		/* B3 = B3 + R3 */ \
\
	xmm0 = xptr(screg2); \
	xmm7 *= xmm0;		/* A2 = A2 * sine (final R2) */ \
	xmm5 *= xmm0;		/* B2 = B2 * sine (final I2) */ \
	xmm1 = xptr(screg5); \
	xmm6 *= xmm1;		/* A4 = A4 * sine (final R4) */ \
	xmm2 *= xmm1;		/* B4 = B4 * sine (final I4) */ \
	xmm0 = xptr(screg1); \
	xmm4 *= xmm0;		/* A3 = A3 * sine (final R3) */ \
	xmm3 *= xmm0;		/* B3 = B3 * sine (final I3) */ \
\

#endif
\
/* 64-bit Intel version. */ \

#if 0 // (@INSTR(,%xarch,<CORE>) NE 0) OR (@INSTR(,%xarch,<P4>) NE 0) OR (@INSTR(,%xarch,<K10>) NE 0) OR (@INSTR(,%xarch,<BLEND>) NE 0);
#if 0 //DEF X86_64;
\
/* Optimal Core 2 timing is 24 clocks, currently at ?? clocks */ \
\

#define r4_x8r_wpn_fft_mem_preload \
	xmm15 = XMM_SQRTHALF2; \
\
\

#define r4_x8r_wpn_fft_mem(mem1,mem2,mem3,mem4,mem5,mem6,mem7,mem8,screg1,screg2,screg5,normreg,pre1,pre2,dst1,dst2) \
	xmm14 = mem4;		/* R4 */ \
	xmm11 = mem8;		/* R8 */ \
	xmm0 = xmm14;		/* Copy R4 */ \
	xmm14 -= xmm11;		/* new R8 = R4 - R8			; 1-3 */ \
\
	xmm4 = mem2;		/* R2 */ \
	xmm5 = mem6;		/* R6 */ \
	xmm1 = xmm4;		/* Copy R2 */ \
	xmm4 -= xmm5;		/* new R6 = R2 - R6			; 2-4 */ \
\
	xmm11 += xmm0;		/* new R4 = R4 + R8			; 3-5 */ \
\
	xmm5 += xmm1;		/* new R2 = R2 + R6			; 4-6 */ \
	xmm14 *= xmm15;		/* R8 = R8 * square root 1/2		; 4-8 */ \
\
	xmm0 = mem1;		/* R1 */ \
	xmm1 = mem5;		/* R5 */ \
	xmm6 = xmm0;		/* Copy R1 */ \
	xmm0 += xmm1;		/* new R1 = R1 + R5			; 5-7 */ \
	xmm4 *= xmm15;		/* R6 = R6 * square root 1/2		; 5-9 */ \
\
	xmm13 = mem3;		/* R3 */ \
	xmm8 = mem7;		/* R7 */ \
	xmm9 = xmm13;		/* Copy R3 */ \
	xmm13 += xmm8; 		/* new R3 = R3 + R7			; 6-8 */ \
\
	xmm6 -= xmm1;		/* new R5 = R1 - R5			; 7-9 */ \
\
	xmm1 = xmm5;		/* Copy R2 */ \
	xmm5 -= xmm11;		/* R2 = R2 - R4 (final I2)		; 8-10 */ \
\
	xmm10 = xmm0;		/* Copy R1 */ \
	xmm0 -= xmm13;		/* R1 = R1 - R3 (final R2)		; 9-11 */ \
\
	xmm7 = xmm4;		/* Copy R6 */ \
	xmm4 -= xmm14;		/* R6 = R6 - R8 (Real part)		; 10-12 */ \
\
	xmm14 += xmm7;		/* R8 = R6 + R8 (Imaginary part)	; 11-13 */ \
	xmm7 = xptr(screg2+16);	/* cosine/sine for w^2n */ \
\
	xmm9 -= xmm8;		/* new R7 = R3 - R7			; 12-14 */ \
	xmm2 = xptr(screg5+16);	/* cosine/sine for w^5n */ \
\
	xmm8 = xmm6;		/* Copy R5 */ \
	xmm6 -= xmm4;		/* R5 = R5 - R6 (final R4)		; 13-15 */ \
	xmm3 = xptr(screg1+16);	/* cosine/sine for w^n */ \
\
	xprefetchw(u8ptr(pre1)); \
\
	xmm4 += xmm8;		/* R6 = R5 + R6 (final R3)		; 14-16 */ \
	xmm12 = xmm5;		/* Copy I2 */ \
	xmm5 *= xmm7;		/* B2 = I2 * cosine/sine		; 14-18 */ \
\
	xmm8 = xmm9;		/* Copy R7 */ \
	xmm9 -= xmm14;		/* R7 = R7 - R8 (final I4)		; 15-17 */ \
	xmm7 *= xmm0;		/* A2 = R2 * cosine/sine		; 15-19 */ \
\
	xmm14 += xmm8;		/* R8 = R7 + R8 (final I3)		; 16-18 */ \
	xmm8 = xmm6;		/* Copy R4 */ \
	xmm6 *= xmm2;		/* A4 = R4 * cosine/sine		; 16-20 */ \
\
	xmm11 += xmm1;		/* R4 = R2 + R4 (final I1, a.k.a 2nd real result) ; 17-19 */ \
	xmm1 = xmm4;		/* Copy R3 */ \
	xmm4 *= xmm3;		/* A3 = R3 * cosine/sine		; 17-21 */ \
\
	xmm13 += xmm10;		/* R3 = R1 + R3 (final R1)		; 18-20 */ \
	xmm2 *= xmm9;		/* B4 = I4 * cosine/sine		; 18-22 */ \
	xmm10 = xptr(screg2);		/* Sine * normalization value */ \
\
	xmm5 += xmm0;		/* B2 = B2 + R2				; 19-21 */ \
	xmm3 *= xmm14;		/* B3 = I3 * cosine/sine		; 19-23 */ \
	xmm0 = xptr(screg5);		/* Sine * normalization value */ \
\
	xprefetchw(u8ptr(pre1+pre2)); \
\
	xmm7 -= xmm12;		/* A2 = A2 - I2				; 20-22 */ \
	xmm11 *= xptr(normreg);	/* I1 * normalization value		; 20-24 */ \
	xmm12 = xptr(screg1);		/* Sine * normalization value */ \
\
	xmm6 -= xmm9;		/* A4 = A4 - I4				; 21-23 */ \
	xmm13 *= xptr(normreg);	/* R1 * normalization value		; 21-25 */ \
\
	xmm4 -= xmm14;		/* A3 = A3 - I3				; 22-24 */ \
	xmm5 *= xmm10;		/* B2 = B2 * sine (final I2)		; 22-26 */ \
\
	xmm2 += xmm8;		/* B4 = B4 + R4				; 23-25 */ \
	xmm7 *= xmm10;		/* A2 = A2 * sine (final R2)		; 23-27 */ \
\
	xmm3 += xmm1;		/* B3 = B3 + R3				; 24-26 */ \
	xmm6 *= xmm0;		/* A4 = A4 * sine (final R4)		; 24-28 */ \
\
	xmm4 *= xmm12;		/* A3 = A3 * sine (final R3)		; 25-29 */ \
	dst2 = xmm11;		/* Save I1				; 25 */ \
\
	xmm2 *= xmm0;		/* B4 = B4 * sine (final I4)		; 26-30 */ \
	dst1 = xmm13;		/* Save R1				; 26 */ \
\
	xmm3 *= xmm12;		/* B3 = B3 * sine (final I3)		; 27-31 */ \
\

#endif
#endif
\
/* 64-bit AMD K8 optimized version.  This is the 32-bit K8 version, but with XMM_TWO and */ \
/* XMM_SQRTHALF preloaded. */ \

#if 0 // (@INSTR(,%xarch,<K8>) NE 0);
#if 0 //DEF X86_64;

/* WARNING:  this preload must be compatible with r4_x4c_wpn_djbfft_partial_mem */ \

#define r4_x8r_wpn_fft_mem_preload \
	xmm14 = XMM_SQRTHALF2; \
	xmm15 = XMM_TWO; \
\
\

#define r4_x8r_wpn_fft_mem(mem1,mem2,mem3,mem4,mem5,mem6,mem7,mem8,screg1,screg2,screg5,normreg,pre1,pre2,dst1,dst2) \
	xmm3 = mem4;		/* R4 */ \
	xmm7 = mem8;		/* R8 */ \
	xmm3 -= xmm7;		/* new R8 = R4 - R8 */ \
	xmm7 += mem4;		/* new R4 = R4 + R8 */ \
\
	xmm1 = mem2;		/* R2 */ \
	xmm5 = mem6;		/* R6 */ \
	xmm1 -= xmm5;		/* new R6 = R2 - R6 */ \
	xmm5 += mem2;		/* new R2 = R2 + R6 */ \
\
	xmm3 *= xmm14;		/* R8 = R8 * square root */ \
	xmm1 *= xmm14;		/* R6 = R6 * square root */ \
\
	xmm0 = mem1;		/* R1 */ \
	xmm4 = mem5;		/* R5 */ \
	xmm0 -= xmm4;		/* new R5 = R1 - R5 */ \
	xmm4 += mem1;		/* new R1 = R1 + R5 */ \
\
	xmm5 -= xmm7;		/* R2 = R2 - R4 (final I2) */ \
	xmm7 *= xmm15; \
	xmm7 += xmm5;		/* R4 = R2 + R4 (final I1, a.k.a 2nd real result) */ \
\
	xmm1 -= xmm3;		/* R6 = R6 - R8 (Real part) */ \
	xmm3 *= xmm15; \
	xmm3 += xmm1;		/* R8 = R6 + R8 (Imaginary part) */ \
\
	xmm6 = mem7;		/* R7 */ \
	xmm6 += mem3;		/* new R3 = R3 + R7 */ \
\
	xmm4 -= xmm6;		/* R1 = R1 - R3 (final R2) */ \
	xmm6 *= xmm15; \
	xmm6 += xmm4;		/* R3 = R1 + R3 (final R1) */ \
\
	xmm2 = mem3;		/* R3 */ \
	xmm2 -= mem7;		/* new R7 = R3 - R7 */ \
\
	xprefetchw(u8ptr(pre1)); \
\
	xmm6 *= xptr(normreg);		/* R1 * normalization value */ \
	dst1 = xmm6;		/* Save R1 */ \
	xmm7 *= xptr(normreg);		/* I1 * normalization value */ \
	dst2 = xmm7;		/* Save I1 */ \
\
	xmm0 -= xmm1;		/* R5 = R5 - R6 (final R4) */ \
	xmm1 *= xmm15; \
	xmm1 += xmm0;		/* R6 = R5 + R6 (final R3) */ \
\
	xmm2 -= xmm3;		/* R7 = R7 - R8 (final I4) */ \
	xmm3 *= xmm15; \
	xmm3 += xmm2;		/* R8 = R7 + R8 (final I3) */ \
\
	xprefetchw(u8ptr(pre1+pre2)); \
\
	xmm7 = xptr(screg2+16);	/* cosine/sine for w^2n */ \
	xmm7 *= xmm4;		/* A2 = R2 * cosine/sine */ \
	xmm7 -= xmm5;		/* A2 = A2 - I2 */ \
	xmm5 *= xptr(screg2+16);	/* B2 = I2 * cosine/sine */ \
	xmm5 += xmm4;		/* B2 = B2 + R2 */ \
\
	xmm6 = xptr(screg5+16);	/* cosine/sine for w^5n */ \
	xmm6 *= xmm0;		/* A4 = R4 * cosine/sine */ \
	xmm6 -= xmm2;		/* A4 = A4 - I4 */ \
	xmm2 *= xptr(screg5+16);	/* B4 = I4 * cosine/sine */ \
	xmm2 += xmm0;		/* B4 = B4 + R4 */ \
\
	xmm4 = xptr(screg1+16);	/* cosine/sine for w^n */ \
	xmm4 *= xmm1;		/* A3 = R3 * cosine/sine */ \
	xmm4 -= xmm3;		/* A3 = A3 - I3 */ \
	xmm3 *= xptr(screg1+16);	/* B3 = I3 * cosine/sine */ \
	xmm3 += xmm1;		/* B3 = B3 + R3 */ \
\
	xmm0 = xptr(screg2); \
	xmm7 *= xmm0;		/* A2 = A2 * sine (final R2) */ \
	xmm5 *= xmm0;		/* B2 = B2 * sine (final I2) */ \
	xmm1 = xptr(screg5); \
	xmm6 *= xmm1;		/* A4 = A4 * sine (final R4) */ \
	xmm2 *= xmm1;		/* B4 = B4 * sine (final I4) */ \
	xmm0 = xptr(screg1); \
	xmm4 *= xmm0;		/* A3 = A3 * sine (final R3) */ \
	xmm3 *= xmm0;		/* B3 = B3 * sine (final I3) */ \
\

#endif
#endif

/* */ \
/* ************************************* half-eight-reals-four-complex-djbfft variants ****************************************** */ \
/* */ \
/* Because of the way we store FFT data elements in cache lines during the pass 2 of the forward FFT, a cache line */ \
/* will contain data needing an eight-reals-fft in the low word only as well as data needing a four-complex-ffts. */ \
/* */ \
\
/* Macro to swizzle and do an eight_reals_fft and a four_complex_djbfft in first levels of pass 2. */ \
/* The eight-reals macro and one of the four-complex only use half the XMM  */ \
/* register.  This isn't very efficient, but this macro isn't called a whole lot. */ \
\
/* This is used in the first levels of pass 2 if pass 2 does the swizzling */
#if 0 //IFDEF UNUSED

#define r4_s2cl_eight_reals_four_complex_djbfft(srcreg,srcinc,d1,screg1,screg2) \
	shuffle_load(xmm0, xmm2, xptr(srcreg+rbx), xptr(srcreg+32+rbx)); /* R1,R3 */ \
	shuffle_load(xmm1, xmm3, xptr(srcreg+d1+rbx), xptr(srcreg+d1+32+rbx)); /* R2,R4 */ \
	shuffle_load(xmm4, xmm6, xptr(srcreg+16+rbx), xptr(srcreg+48+rbx)); /* R5,R7 */ \
	shuffle_load(xmm5, xmm7, xptr(srcreg+d1+16+rbx), xptr(srcreg+d1+48+rbx)); /* R6,R8 */ \
	xptr(srcreg) = xmm0; \
	xptr(srcreg+16) = xmm1; \
	xptr(srcreg+32) = xmm2; \
	xptr(srcreg+48) = xmm3; \
	xptr(srcreg+d1) = xmm4; \
	xptr(srcreg+d1+16) = xmm5; \
	xptr(srcreg+d1+32) = xmm6; \
	xptr(srcreg+d1+48) = xmm7; \
	r4_h8r_h4c_djbfft_mem(xptr(srcreg),xptr(srcreg+16),xptr(srcreg+32),xptr(srcreg+48),xptr(srcreg+d1),xptr(srcreg+d1+16),xptr(srcreg+d1+32),xptr(srcreg+d1+48),screg1,screg2,srcreg+srcinc,d1,xptr(srcreg),xptr(srcreg+16)); \
/*	xstore	[srcreg], xmm0		;; Save real value #1 */ \
/*	xstore	[srcreg+16], xmm0	;; Save real value #2 */ \
	xptr(srcreg+32) = xmm7;	/* Save R2 */ \
	xptr(srcreg+48) = xmm6;	/* Save I2 */ \
	xptr(srcreg+d1) = xmm3;	/* Save R3 */ \
	xptr(srcreg+d1+16) = xmm1;	/* Save I3 */ \
	xptr(srcreg+d1+32) = xmm0;	/* Save R4 */ \
	xptr(srcreg+d1+48) = xmm4;	/* Save I4 */ \
	srcreg += srcinc; \

#endif
\
/* This is used in the first levels of pass 2 if pass 1 does the swizzling */ \

#define r4_fh2cl_eight_reals_four_complex_djbfft(srcreg,srcinc,d1,screg1,screg2) { \
	uintptr_t src_rbx = srcreg+rbx; \
	vec2f64 *sd0 = (vec2f64*)(srcreg); \
	vec2f64 *sd1 = (vec2f64*)(srcreg+d1); \
	vec2f64 *bd0 = (vec2f64*)(src_rbx); \
	vec2f64 *bd1 = (vec2f64*)(src_rbx+d1); \
	r4_h8r_h4c_djbfft_mem(bd0[0],bd1[0],bd0[1],bd1[1],bd0[2],bd1[2],bd0[3],bd1[3], screg1,screg2,srcreg+srcinc,d1, sd0[0],sd0[1]); \
/*	xstore	[srcreg], xmm0		;; Save real value #1 */ \
/*	xstore	[srcreg+16], xmm0	;; Save real value #2 */ \
	sd0[2] = xmm7;	/* Save R2 */ \
	sd0[3] = xmm6;	/* Save I2 */ \
	sd1[0] = xmm3;	/* Save R3 */ \
	sd1[1] = xmm1;	/* Save I3 */ \
	sd1[2] = xmm0;	/* Save R4 */ \
	sd1[3] = xmm4;	/* Save I4 */ \
	srcreg += srcinc; \
}

/* This is used in the later levels of pass 2 if first levels used a radix-3 or radix-5. */ \

#define r4_h2cl_eight_reals_four_complex_djbfft(srcreg,srcinc,d1,screg1,screg2) { \
	vec2f64 *sd0 = (vec2f64*)(srcreg); \
	vec2f64 *sd1 = (vec2f64*)(srcreg+d1); \
	r4_h8r_h4c_djbfft_mem(sd0[0],sd1[0],sd0[1],sd1[1],sd0[2],sd1[2],sd0[3],sd1[3], screg1,screg2,srcreg+srcinc,d1, sd0[0],sd0[1]); \
/*	xstore	[srcreg], xmm0		;; Save real value #1 */ \
/*	xstore	[srcreg+16], xmm0	;; Save real value #2 */ \
	sd0[2] = xmm7;	/* Save R2 */ \
	sd0[3] = xmm6;	/* Save I2 */ \
	sd1[0] = xmm3;	/* Save R3 */ \
	sd1[1] = xmm1;	/* Save I3 */ \
	sd1[2] = xmm0;	/* Save R4 */ \
	sd1[3] = xmm4;	/* Save I4 */ \
	srcreg += srcinc; \
}

/* Like r4_h2cl_eight_reals_four_complex_djbfft except it uses 2 sin/cos ptrs for the four-complex djbfft. */ \
/* This is used in the later levels of pass 2 if first levels used a radix-3. */ \

#define r4_h2cl_2sc_eight_reals_four_complex_djbfft(srcreg,srcinc,d1,screg1,screg2,screg3) { \
\
	vec2f64 *sd0 = (vec2f64*)(srcreg); \
	vec2f64 *sd1 = (vec2f64*)(srcreg+d1); \
	r4_h8r_h4c_2sc_djbfft_mem(sd0[0],sd1[0],sd0[1],sd1[1],sd0[2],sd1[2],sd0[3],sd1[3], screg1,screg2,screg3,srcreg+srcinc,d1, sd0[0],sd0[1]); \
/*	xstore	[srcreg], xmm0		;; Save real value #1 */ \
/*	xstore	[srcreg+16], xmm0	;; Save real value #2 */ \
	sd0[2] = xmm7;	/* Save R2 */ \
	sd0[3] = xmm6;	/* Save I2 */ \
	sd1[0] = xmm3;	/* Save R3 */ \
	sd1[1] = xmm1;	/* Save I3 */ \
	sd1[2] = xmm0;	/* Save R4 */ \
	sd1[3] = xmm4;	/* Save I4 */ \
	srcreg += srcinc; \
}

/* Macro to do an eight_reals_fft and three four_complex_djbfft in pass 2. */ \
/* The eight-reals macro and one of the four-complex only use half the XMM  */ \
/* register.  This isn't very efficient, but this macro isn't called a whole lot. */ \
\

#define r4_h4cl_eight_reals_four_complex_djbfft_preload \
	r4_x4c_djbfft_partial_mem_preload; \
\
\

#define r4_h4cl_eight_reals_four_complex_djbfft(srcreg,srcinc,d1,d2,screg1,screg2) { \
	vec2f64 *sd0 = (vec2f64*)(srcreg); \
	vec2f64 *sd1 = (vec2f64*)(srcreg+d1); \
	vec2f64 *sd2 = (vec2f64*)(srcreg+d2); \
	vec2f64 *sd3 = (vec2f64*)(srcreg+d2+d1); \
	r4_h8r_h4c_djbfft_mem(sd0[0],sd1[0],sd2[0],sd3[0],sd0[1],sd1[1],sd2[1],sd3[1], screg1,screg2,srcreg+srcinc,d1,sd0[0],sd0[1]); \
/*	xstore	[srcreg], xmm0		;; Save real value #1 */ \
/*	xstore	[srcreg+16], xmm0	;; Save real value #2 */ \
	sd1[0] = xmm3;	/* Save R3 */ \
	sd1[1] = xmm1;	/* Save I3 */ \
	xmm1 = sd0[2];	/* R1 */ \
	xmm3 = sd0[3];	/* R5 */ \
	sd0[2] = xmm7;	/* Save R2 */ \
	sd0[3] = xmm6;	/* Save I2 */ \
	xmm7 = sd1[2];	/* R2 */ \
	xmm6 = sd1[3];	/* R6 */ \
	sd1[2] = xmm0;	/* Save R4 */ \
	sd1[3] = xmm4;	/* Save I4 */ \
	r4_x4c_djbfft_partial_mem(xmm1,xmm7,xmm0,xmm2,xmm3,xmm6,xmm4,xmm5, sd2[2],sd3[2],sd2[3],sd3[3], screg1,0,srcreg+srcinc+d2,d1, sd2[0],sd2[1]); \
/*	xstore	[srcreg+d2], xmm0	;; Save R1 */ \
/*	xstore	[srcreg+d2+16], xmm0	;; Save I1 */ \
	sd2[2] = xmm5;	/* Save R2 */ \
	sd2[3] = xmm4;	/* Save I2 */ \
	sd3[0] = xmm2;	/* Save R3 */ \
	sd3[1] = xmm7;	/* Save I3 */ \
	sd3[2] = xmm1;	/* Save R4 */ \
	sd3[3] = xmm3;	/* Save I4 */ \
	srcreg += srcinc; \
}

/* Like r4_h4cl_eight_reals_four_complex_djbfft but takes 2 sin/cos ptrs for the 4 complex djbfft. */ \
/* Used in pass 2 when the first levels are a three-complex FFT. */ \
\

#define r4_h4cl_2sc_eight_reals_four_complex_djbfft_preload \
	r4_x4c_2sc_djbfft_partial_mem_preload; \
\
\

#define r4_h4cl_2sc_eight_reals_four_complex_djbfft(srcreg,srcinc,d1,d2,screg1,screg2,screg3) { \
	vec2f64 *sd0 = (vec2f64*)(srcreg); \
	vec2f64 *sd1 = (vec2f64*)(srcreg+d1); \
	vec2f64 *sd2 = (vec2f64*)(srcreg+d2); \
	vec2f64 *sd3 = (vec2f64*)(srcreg+d2+d1); \
	r4_h8r_h4c_2sc_djbfft_mem(sd0[0],sd1[0],sd2[0],sd3[0],sd0[1],sd1[1],sd2[1],sd3[1], screg1,screg2,screg3,srcreg+srcinc,d1, sd0[0],sd0[1]); \
/*	xstore	[srcreg], xmm0		;; Save real value #1 */ \
/*	xstore	[srcreg+16], xmm0	;; Save real value #2 */ \
	sd1[0] = xmm3;	/* Save R3 */ \
	sd1[1] = xmm1;	/* Save I3 */ \
	xmm1 = sd0[2];	/* R1 */ \
	xmm3 = sd0[3];	/* R5 */ \
	sd0[2] = xmm7;	/* Save R2 */ \
	sd0[3] = xmm6;	/* Save I2 */ \
	xmm7 = sd1[2];	/* R2 */ \
	xmm6 = sd1[3];	/* R6 */ \
	sd1[2] = xmm0;	/* Save R4 */ \
	sd1[3] = xmm4;	/* Save I4 */ \
	r4_x4c_2sc_djbfft_partial_mem(xmm1,xmm7,xmm0,xmm2,xmm3,xmm6,xmm4,xmm5, sd2[2],sd3[2],sd2[3],sd3[3], screg1,screg2,srcreg+srcinc+d2,d1, sd2[0],sd2[1]); \
/*	xstore	[srcreg+d2], xmm0	;; Save R1 */ \
/*	xstore	[srcreg+d2+16], xmm0	;; Save I1 */ \
	sd2[2] = xmm5;	/* Save R2 */ \
	sd2[3] = xmm4;	/* Save I2 */ \
	sd3[0] = xmm2;	/* Save R3 */ \
	sd3[1] = xmm7;	/* Save I3 */ \
	sd3[2] = xmm1;	/* Save R4 */ \
	sd3[3] = xmm3;	/* Save I4 */ \
	srcreg += srcinc; \
}


/* Does an r4_x8r_fft_mem on the low word of the xmm register */ \
/* Does an r4_x4c_djbfft_mem on the high word of the xmm register */ \

#define r4_h8r_h4c_djbfft_mem(mem1,mem2,mem3,mem4,mem5,mem6,mem7,mem8,screg1,screg2,pre1,pre2,dst1,dst2) \
	r4_h8r_h4c_2sc_djbfft_mem(mem1,mem2,mem3,mem4,mem5,mem6,mem7,mem8,screg1,screg1+32,screg2,pre1,pre2,dst1,dst2); \
\

#if 1
#define r4_h8r_h4c_2sc_djbfft_mem(mem1,mem2,mem3,mem4,mem5,mem6,mem7,mem8,screg1,screg2,screg3,pre1,pre2,dst1,dst2) { \
	vec2f64 M1=mem1, M2=mem2, M3=mem3, M4=mem4, M5=mem5, M6=mem6, M7=mem7, M8=mem8; \
	xprefetchw(u8ptr(pre1)); \
	xprefetchw(u8ptr(pre1+pre2)); \
\
	/* [0]: Do the eight reals part */ \
	/* [1]: Do the four-complex part */ \
\
	xmm1[0] = M4[0];		/* R4 */ \
	xmm7[0] = M8[0];		/* R8 */ \
	xmm1[0] -= xmm7[0];		/* new R8 = R4 - R8 */ \
	xmm7[0] += M4[0];		/* new R4 = R4 + R8 */ \
	xmm5[0] = M2[0];		/* R2 */ \
	xmm6[0] = M6[0];		/* R6 */ \
	xmm5[0] -= xmm6[0];		/* new R6 = R2 - R6 */ \
	xmm6[0] += M2[0];		/* new R2 = R2 + R6 */ \
\
	xmm1[0] *= XMM_SQRTHALF1;	/* R8 = R8 * square root */ \
	xmm5[0] *= XMM_SQRTHALF1;	/* R6 = R6 * square root */ \
\
	xmm2[0] = M1[0];		/* R1 */ \
	xmm3[0] = M5[0];		/* R5 */ \
	xmm2[0] -= xmm3[0];		/* new R5 = R1 - R5 */ \
	xmm3[0] += M1[0];		/* new R1 = R1 + R5 */ \
\
	xmm0[0] = xmm5[0];		/* Copy R6 */ \
	xmm5[0] -= xmm1[0];		/* R6 = R6 - R8 (Real part) */ \
	xmm1[0] += xmm0[0];		/* R8 = R6 + R8 (Imaginary part) */ \
\
	xmm0[0] = M7[0];		/* R7 */ \
	xmm0[0] += M3[0];		/* new R3 = R3 + R7 */ \
\
	xmm5[1] = M1[1];		/* R1 */ \
	xmm3[1] = M3[1];		/* R3 */ \
	xmm3[1] += xmm5[1];		/* R3 = R1 + R3 (new R1) */ \
	xmm5[1] -= M3[1];		/* R1 = R1 - R3 (new R3) */ \
\
	xmm1[1] = M2[1];		/* R2 */ \
	xmm0[1] = M4[1];		/* R4 */ \
	xmm0[1] += xmm1[1];		/* R4 = R2 + R4 (new R2) */ \
	xmm1[1] -= M4[1];		/* R2 = R2 - R4 (new R4) */ \
\
	xmm4 = xmm3;			/* Copy R1 */ \
	xmm3 -= xmm0;			/* R1R1 = R1R1 - R3R2 (final R2R2) */ \
	xmm0 += xmm4;			/* R3R2 = R1R1 + R3R2 (final R1R1) */ \
\
	dst1 = xmm0;			/* Save R1 R1 */ \
\
	xmm2[1] = M6[1];		/* I2 */ \
	xmm7[1] = M8[1];		/* I4 */ \
	xmm2[1] -= xmm7[1];		/* I2 = I2 - I4 (new I4) */ \
	xmm7[1] += M6[1];		/* I4 = I2 + I4 (new I2) */ \
\
	xmm4[0] = M3[0];		/* R3 */ \
	xmm4[0] -= M7[0];		/* new R7 = R3 - R7 */ \
\
	xmm4[1] = M5[1];		/* I1 */ \
	xmm6[1] = M7[1];		/* I3 */ \
	xmm4[1] -= xmm6[1];		/* I1 = I1 - I3 (new I3) */ \
	xmm6[1] += M5[1];		/* I3 = I1 + I3 (new I1) */ \
\
	xmm0[0] = xmm5[0]; \
	xmm5[0] += xmm2[0];		/* R6 = R5 + R6 (final R3) */ \
	xmm2[0] -= xmm0[0];		/* R5 = R5 - R6 (final R4) */ \
\
	xmm0[1] = xmm5[1]; \
	xmm5[1] -= xmm2[1];		/* R3 = R3 - I4 (final R3) */ \
	xmm2[1] += xmm0[1];		/* I4 = R3 + I4 (final R4) */ \
\
	xmm0 = xmm6; \
	xmm6 -= xmm7;			/* R2,I1 = R2,I1 - R4,I2 (final I2,I2) */ \
	xmm7 += xmm0;			/* R4,I2 = R2,I1 + R4,I2 (final I1,I1, a.k.a 2nd real result) */ \
\
	dst2 = xmm7;			/* Save I1 I1 */ \
\
\
	xmm0 = xmm4; \
	xmm4 -= xmm1;			/* R7,I3 = R7,I3 - R8,R4 (final I4,I4) */ \
	xmm1 += xmm0;			/* R8,R4 = R7,I3 + R8,R4 (final I3,I3) */ \
\
	double sc1_16 = f64ptr(screg1+16); /* cosine/sine for w^2n */ \
	double sc2_16 = f64ptr(screg2+16); \
	vec2f64 sc; \
\
	sc[0] = sc1_16; sc[1] = sc2_16; \
	xmm7 = xmm3 * sc - xmm6;	/* A2 = R2 * cosine/sine - I2 */ \
	xmm6 = xmm6 * sc + xmm3;	/* B2 = I2 * cosine/sine + R2 */ \
\
	vec2f64 sc3_0 = xptr(screg3); /* cosine/sine for w^n */ \
\
	sc[0] = sc3_0[1]; sc[1] = sc1_16; \
	xmm3 = xmm5 * sc - xmm1;	/* A3 = R3 * cosine/sine - I3 */ \
	xmm1 = xmm1 * sc + xmm5;	/* B3 = I3 * cosine/sine + R3 */ \
\
	vec2f64 sc3_16 = -xptr(screg3+16); /* cosine/sine for w^5n (negative) */ \
\
	sc[0] = sc3_16[1]; sc[1] = sc1_16; \
	xmm0 = xmm2 * sc + xmm4;	/* A4 = R4 * cosine/sine (-,+) I4 */ \
	xmm4 = xmm4 * sc - xmm2;	/* B4 = I4 * cosine/sine (+,-) R4 */ \
\
	xmm5[0] = sc3_0[0]; xmm5[1] = f64ptr(screg1); \
	xmm3 *= xmm5;			/* A3 = A3 * sine (final R3) */ \
	xmm1 *= xmm5;			/* B3 = B3 * sine (final I3) */ \
\
	xmm5[0] = sc3_16[0]; \
	xmm0 *= xmm5;			/* A4 = A4 * sine (final R4) */ \
	xmm4 *= xmm5;			/* B4 = B4 * sine (final I4) */ \
\
	xmm5[0] = xmm5[1]; xmm5[1] = f64ptr(screg2); \
	xmm7 *= xmm5;			/* A2 = A2 * sine (final R2) */ \
	xmm6 *= xmm5;			/* B2 = B2 * sine (final I2) */ \
}
#else
#define r4_h8r_h4c_2sc_djbfft_mem(mem1,mem2,mem3,mem4,mem5,mem6,mem7,mem8,screg1,screg2,screg3,pre1,pre2,dst1,dst2) { \
	double f0,f1,f2,f3,f4,f5,f6,f7; \
;	/* Do the four-complex part */ \
\
	f0 = mem1[1];		/* R1 */ \
	f2 = mem3[1];		/* R3 */ \
	f2 += f0;		/* R3 = R1 + R3 (new R1) */ \
	f0 -= mem3[1];		/* R1 = R1 - R3 (new R3) */ \
\
	f1 = mem2[1];		/* R2 */ \
	f3 = mem4[1];		/* R4 */ \
	f3 += f1;		/* R4 = R2 + R4 (new R2) */ \
	f1 -= mem4[1];		/* R2 = R2 - R4 (new R4) */ \
\
	f6 = f2; \
	f2 -= f3;		/* R1 = R1 - R2 (final R2) */ \
	f3 += f6;		/* R2 = R1 + R2 (final R1) */ \
\
	dst1[1] = f3;		/* Save R1 */ \
\
	f5 = mem6[1];		/* I2 */ \
	f7 = mem8[1];		/* I4 */ \
	f5 -= f7;		/* I2 = I2 - I4 (new I4) */ \
	f7 += mem6[1];		/* I4 = I2 + I4 (new I2) */ \
\
	f4 = mem5[1];		/* I1 */ \
	f6 = mem7[1];		/* I3 */ \
	f4 -= f6;		/* I1 = I1 - I3 (new I3) */ \
	f6 += mem5[1];		/* I3 = I1 + I3 (new I1) */ \
\
	f3 = f0; \
	f0 -= f5;		/* R3 = R3 - I4 (final R3) */ \
	f5 += f3;		/* I4 = R3 + I4 (final R4) */ \
\
	f3 = f1; \
	f1 += f4;		/* R4 = I3 + R4 (final I3) */ \
	f4 -= f3;		/* I3 = I3 - R4 (final I4) */ \
\
	f3 = f6; \
	f6 -= f7;		/* I1 = I1 - I2 (final I2) */ \
	f7 += f3;		/* I2 = I1 + I2 (final I1) */ \
\
	f3 = f64ptr(screg1+16);	/* cosine/sine */ \
	f3 *= f0;		/* A3 = R3 * cosine/sine */ \
\
	dst2[1] = f7;		/* Save I1 */ \
\
	f7 = f64ptr(screg2+16);	/* cosine/sine */ \
	f7 *= f2;		/* A2 = R2 * cosine/sine */ \
\
	f3 -= f1;		/* A3 = A3 - I3 */ \
	f1 *= f64ptr(screg1+16);	/* B3 = I3 * cosine/sine */ \
	f1 += f0;		/* B3 = B3 + R3 */ \
\
	f0 = f64ptr(screg1+16);	/* cosine/sine */ \
	f0 *= f5;		/* A4 = R4 * cosine/sine */ \
\
	f7 -= f6;		/* A2 = A2 - I2 */ \
	f6 *= f64ptr(screg2+16);	/* B2 = I2 * cosine/sine */ \
\
	f0 += f4;		/* A4 = A4 + I4 */ \
	f4 *= f64ptr(screg1+16);	/* B4 = I4 * cosine/sine */ \
\
	f6 += f2;		/* B2 = B2 + R2 */ \
	f4 -= f5;		/* B4 = B4 - R4 */ \
\
	f2 = f64ptr(screg1);	/* sine */ \
	f3 *= f2;		/* A3 = A3 * sine (final R3) */ \
	f1 *= f2;		/* B3 = B3 * sine (final I3) */ \
	f7 *= f64ptr(screg2);	/* A2 = A2 * sine (final R2) */ \
	f0 *= f2;		/* A4 = A4 * sine (final R4) */ \
	f6 *= f64ptr(screg2);	/* B2 = B2 * sine (final I2) */ \
	f4 *= f2;		/* B4 = B4 * sine (final I4) */ \
\
;	/* Do the eight reals part */ \
\
	xmm3[1] = f3;		/* Copy to high part of XMM register */ \
	xmm1[1] = f1; \
	xmm7[1] = f7; \
	xmm0[1] = f0; \
	xmm6[1] = f6; \
	xmm4[1] = f4; \
\
	f1 = mem4[0];		/* R4 */ \
	f7 = mem8[0];		/* R8 */ \
	f1 -= f7;		/* new R8 = R4 - R8 */ \
	f7 += mem4[0];		/* new R4 = R4 + R8 */ \
	f5 = mem2[0];		/* R2 */ \
	f6 = mem6[0];		/* R6 */ \
	f5 -= f6;		/* new R6 = R2 - R6 */ \
	f6 += mem2[0];		/* new R2 = R2 + R6 */ \
\
	/*vec2f64 t51; \
	t51[0] = f1; t51[1] = f5; \
	t51 *= XMM_SQRTHALF2; \
	f1 = t51[0]; f5 = t51[1];*/ \
	f1 *= XMM_SQRTHALF1;	/* R8 = R8 * square root */ \
	f5 *= XMM_SQRTHALF1;	/* R6 = R6 * square root */ \
\
	f2 = mem1[0];		/* R1 */ \
	f3 = mem5[0];		/* R5 */ \
	f2 -= f3;		/* new R5 = R1 - R5 */ \
	f3 += mem1[0];		/* new R1 = R1 + R5 */ \
\
	f4 = f6;		/* Copy R2 */ \
	f6 -= f7;		/* R2 = R2 - R4 (final I2) */ \
	f7 += f4;		/* R4 = R2 + R4 (final I1, a.k.a 2nd real result) */ \
\
	f0 = f5;		/* Copy R6 */ \
	f5 -= f1;		/* R6 = R6 - R8 (Real part) */ \
	f1 += f0;		/* R8 = R6 + R8 (Imaginary part) */ \
\
	f0 = mem7[0];		/* R7 */ \
	f0 += mem3[0];		/* new R3 = R3 + R7 */ \
\
	f4 = f3;		/* Copy R1 */ \
	f3 -= f0;		/* R1 = R1 - R3 (final R2) */ \
	f0 += f4;		/* R3 = R1 + R3 (final R1) */ \
\
	f4 = mem3[0];		/* R3 */ \
	f4 -= mem7[0];		/* new R7 = R3 - R7 */ \
\
	xprefetchw(u8ptr(pre1)); \
\
	dst1[0] = f0;		/* Save R1 */ \
	dst2[0] = f7;		/* Save I1 */ \
\
	f7 = f2;		/* Copy R5 */ \
	f2 -= f5;		/* R5 = R5 - R6 (final R4) */ \
	f5 += f7;		/* R6 = R5 + R6 (final R3) */ \
\
	f0 = f4;		/* Copy R7 */ \
	f4 -= f1;		/* R7 = R7 - R8 (final I4) */ \
	f1 += f0;		/* R8 = R7 + R8 (final I3) */ \
\
	xprefetchw(u8ptr(pre1+pre2)); \
\
	f7 = f64ptr(screg1+16);	/* cosine/sine for w^2n */ \
	f7 *= f3;		/* A2 = R2 * cosine/sine */ \
	f7 -= f6;		/* A2 = A2 - I2 */ \
	f6 *= f64ptr(screg1+16);	/* B2 = I2 * cosine/sine */ \
	f6 += f3;		/* B2 = B2 + R2 */ \
\
	f0 = f64ptr(screg3+24);	/* cosine/sine for w^5n */ \
	f0 *= f2;		/* A4 = R4 * cosine/sine */ \
	f0 -= f4;		/* A4 = A4 - I4 */ \
	f4 *= f64ptr(screg3+24);	/* B4 = I4 * cosine/sine */ \
	f4 += f2;		/* B4 = B4 + R4 */ \
\
	f3 = f64ptr(screg3+8);	/* cosine/sine for w^n */ \
	f3 *= f5;		/* A3 = R3 * cosine/sine */ \
	f3 -= f1;		/* A3 = A3 - I3 */ \
	f1 *= f64ptr(screg3+8);	/* B3 = I3 * cosine/sine */ \
	f1 += f5;		/* B3 = B3 + R3 */ \
\
	f7 *= f64ptr(screg1);	/* A2 = A2 * sine (final R2) */ \
	f6 *= f64ptr(screg1);	/* B2 = B2 * sine (final I2) */ \
	f0 *= f64ptr(screg3+16);	/* A4 = A4 * sine (final R4) */ \
	f4 *= f64ptr(screg3+16);	/* B4 = B4 * sine (final I4) */ \
	f3 *= f64ptr(screg3);	/* A3 = A3 * sine (final R3) */ \
	f1 *= f64ptr(screg3);	/* B3 = B3 * sine (final I3) */ \
\
	xmm0[0] = f0; \
	xmm1[0] = f1; \
	xmm3[0] = f3; \
	xmm4[0] = f4; \
	xmm6[0] = f6; \
	xmm7[0] = f7; \
}
#endif

/* */ \
/* ************************************* half-eight-reals-four-complex-djbunfft variants ****************************************** */ \
/* */ \
\
/* Macro to do an eight_reals_unfft and a four_complex_djbunfft  */ \
/* and swizzle the results in the last levels of pass 2. */ \
/* The eight-reals operation is done in the lower half of the XMM */ \
/* register.  This isn't very efficient, but this macro isn't called a whole lot. */ \

#if 0 //IFDEF UNUSED

#define r4_s4cl_eight_reals_four_complex_djbunfft(srcreg,srcinc,d1,d2,screg1,scoff1,screg2,scoff2) { \
	const uintptr_t d3 = d2 + d1; \
	r4_h8r_h4c_djbunfft_mem(xptr(srcreg+16),xptr(srcreg+48),xptr(srcreg+d1+16),xptr(srcreg+d1+48),xptr(srcreg+d2+16),xptr(srcreg+d2+48),xptr(srcreg+d3+16),xptr(srcreg+d3+48),screg1,scoff1,screg2,scoff2,srcreg+srcinc+d2,d1,xptr(g->u.xmm.XMM_TMP1),xptr(g->u.xmm.XMM_TMP2)); \
	xmm2 = xptr(g->u.xmm.XMM_TMP1);		/* Load R4 */ \
	xmm5 = xptr(g->u.xmm.XMM_TMP2);		/* Load R8 */ \
	shuffle_store(xptr(srcreg+d1+16), xptr(srcreg+d1+48), xmm0, xmm6); /* Save R5, R7 */ \
	shuffle_store(xptr(srcreg+d3+16), xptr(srcreg+d3+48), xmm3, xmm5); /* Save R6, R8 */ \
	xmm0 = xptr(srcreg+d1);	/* Load R2 */ \
	xmm6 = xptr(srcreg+d1+32);	/* Load I2 */ \
	xmm3 = xptr(srcreg+d3);	/* Load R4 */ \
	xmm5 = xptr(srcreg+d3+32);	/* Load I4 */ \
	shuffle_store(xptr(srcreg+d1), xptr(srcreg+d1+32), xmm1, xmm7); /* Save R1, R3 */ \
	shuffle_store(xptr(srcreg+d3), xptr(srcreg+d3+32), xmm4, xmm2); /* Save R2, R4 */ \
	xptr(g->u.xmm.XMM_TMP1) = xmm0; \
	xptr(g->u.xmm.XMM_TMP2) = xmm6; \
	xptr(g->u.xmm.XMM_TMP3) = xmm3; \
	xptr(g->u.xmm.XMM_TMP4) = xmm5; \
	r4_h8r_h4c_djbunfft_mem(xptr(srcreg),xptr(srcreg+32),xptr(g->u.xmm.XMM_TMP1),xptr(g->u.xmm.XMM_TMP2),xptr(srcreg+d2),xptr(srcreg+d2+32),xptr(g->u.xmm.XMM_TMP3),xptr(g->u.xmm.XMM_TMP4),screg1,0,screg2,0,srcreg+srcinc,d1,xptr(g->u.xmm.XMM_TMP1),xptr(g->u.xmm.XMM_TMP2)); \
	xmm2 = xptr(g->u.xmm.XMM_TMP1);		/* Load R4 */ \
	xmm5 = xptr(g->u.xmm.XMM_TMP2);		/* Load R8 */ \
	shuffle_store(xptr(srcreg), xptr(srcreg+32), xmm1, xmm7); /* Save R1, R3 */ \
	shuffle_store(xptr(srcreg+16), xptr(srcreg+48), xmm0, xmm6); /* Save R5, R7 */ \
	shuffle_store(xptr(srcreg+d2), xptr(srcreg+d2+32), xmm4, xmm2); /* Save R2, R4 */ \
	shuffle_store(xptr(srcreg+d2+16), xptr(srcreg+d2+48), xmm3, xmm5); /* Save R6, R8 */ \
	srcreg += srcinc;
}
#endif
\
/* Macro to do an eight_reals_unfft and a four_complex_djbunfft in pass 2. */ \
/* The eight-reals operation is done in the lower half of the XMM */ \
/* register.  This isn't very efficient, but this macro isn't called a whole lot. */ \
\

#define r4_h4cl_eight_reals_four_complex_djbunfft(srcreg,srcinc,d1,d2,screg1,scoff1,screg2,scoff2) { \
	vec2f64 *sd0 = (vec2f64*)(srcreg); \
	vec2f64 *sd1 = (vec2f64*)(srcreg+d1); \
	vec2f64 *sd2 = (vec2f64*)(srcreg+d2); \
	vec2f64 *sd3 = (vec2f64*)(srcreg+d2+d1); \
	vec2f64 TMP1, TMP2, TMP3, TMP4; \
	r4_h8r_h4c_djbunfft_mem(sd0[1],sd0[3],sd1[1],sd1[3], sd2[1],sd2[3],sd3[1],sd3[3], screg1,scoff1,screg2,scoff2, srcreg+srcinc+d2,d1,sd3[1],sd3[3]); \
/*	xstore	[srcreg+d3+16], xmm5	;; Save R4 */ \
/*	xstore	[srcreg+d3+48], xmm2	;; Save R8 */ \
	TMP3 = sd3[0];	/* Load R4 */ \
	TMP4 = sd3[2];	/* Load I4 */ \
	sd3[0] = xmm4;	/* Save R2 */ \
	sd3[2] = xmm3;	/* Save R6 */ \
	TMP1 = sd1[0];	/* Load R2 */ \
	TMP2 = sd1[2];	/* Load I2 */ \
	sd1[0] = xmm1;	/* Save R1 */ \
	sd1[1] = xmm7;	/* Save R3 */ \
	sd1[2] = xmm0;	/* Save R5 */ \
	sd1[3] = xmm6;	/* Save R7 */ \
	r4_h8r_h4c_djbunfft_mem(sd0[0],sd0[2],TMP1,TMP2, sd2[0],sd2[2],TMP3,TMP4, screg1,0,screg2,0, srcreg+srcinc,d1,sd2[1],sd2[3]); \
	sd0[0] = xmm1;		/* Save R1 */ \
	sd2[0] = xmm4;	/* Save R2 */ \
	sd0[1] = xmm7;	/* Save R3 */ \
/*	xstore	[srcreg+d2+16], xmm5	;; Save R4 */ \
	sd0[2] = xmm0;	/* Save R5 */ \
	sd2[2] = xmm3;	/* Save R6 */ \
	sd0[3] = xmm6;	/* Save R7 */ \
/*	xstore	[srcreg+d2+48], xmm2	;; Save R8 */ \
	srcreg += srcinc; \
}

/* Like r4_h4cl_eight_reals_four_complex_djbunfft but uses 2 sin/cos ptrs. */ \
/* Used in pass 2 when the first levels are radix-3. */ \

#define r4_h4cl_2sc_eight_reals_four_complex_djbunfft(srcreg,srcinc,d1,d2,screg1,scoff1,screg2,scoff2,screg3,scoff3) { \
	vec2f64 *sd0 = (vec2f64*)(srcreg); \
	vec2f64 *sd1 = (vec2f64*)(srcreg+d1); \
	vec2f64 *sd2 = (vec2f64*)(srcreg+d2); \
	vec2f64 *sd3 = (vec2f64*)(srcreg+d2+d1); \
	vec2f64 TMP1, TMP2, TMP3, TMP4; \
	r4_h8r_h4c_2sc_djbunfft_mem(sd0[1],sd0[3],sd1[1],sd1[3], sd2[1],sd2[3],sd3[1],sd3[3], screg1+scoff1,screg2+scoff2,screg3+scoff3, srcreg+srcinc+d2,d1,sd3[1],sd3[3]); \
/*	xstore	[srcreg+d3+16], xmm5	;; Save R4 */ \
/*	xstore	[srcreg+d3+48], xmm2	;; Save R8 */ \
	TMP3 = sd3[0];	/* Load R4 */ \
	TMP4 = sd3[2];	/* Load I4 */ \
	sd3[0] = xmm4;	/* Save R2 */ \
	sd3[2] = xmm3;	/* Save R6 */ \
	TMP1 = sd1[0];	/* Load R2 */ \
	TMP2 = sd1[2];	/* Load I2 */ \
	sd1[0] = xmm1;	/* Save R1 */ \
	sd1[1] = xmm7;	/* Save R3 */ \
	sd1[2] = xmm0;	/* Save R5 */ \
	sd1[3] = xmm6;	/* Save R7 */ \
	r4_h8r_h4c_2sc_djbunfft_mem(sd0[0],sd0[2],TMP1,TMP2, sd2[0],sd2[2],TMP3,TMP4, screg1,screg2,screg3,srcreg+srcinc,d1, sd2[1],sd2[3]); \
	sd0[0] = xmm1;		/* Save R1 */ \
	sd2[0] = xmm4;	/* Save R2 */ \
	sd0[1] = xmm7;	/* Save R3 */ \
/*	xstore	[srcreg+d2+16], xmm5	;; Save R4 */ \
	sd0[2] = xmm0;	/* Save R5 */ \
	sd2[2] = xmm3;	/* Save R6 */ \
	sd0[3] = xmm6;	/* Save R7 */ \
/*	xstore	[srcreg+d2+48], xmm2	;; Save R8 */ \
	srcreg += srcinc; \
}


#define r4_h8r_h4c_djbunfft_mem(mem1,mem2,mem3,mem4,mem5,mem6,mem7,mem8, screg1,scoff1,screg2,scoff2, pre1,pre2,dst4,dst8) \
	r4_h8r_h4c_2sc_djbunfft_mem(mem1,mem2,mem3,mem4,mem5,mem6,mem7,mem8, screg1+scoff1,screg1+scoff1+32,screg2+scoff2, pre1,pre2,dst4,dst8); \
\

#if 1
#define r4_h8r_h4c_2sc_djbunfft_mem(mem1,mem2,mem3,mem4,mem5,mem6,mem7,mem8, screg1,screg2,screg3, pre1,pre2,dst4,dst8) { \
	vec2f64 M1=mem1, M2=mem2, M3=mem3, M4=mem4, M5=mem5, M6=mem6, M7=mem7, M8=mem8; \
	xprefetchw(u8ptr(pre1)); \
	xprefetchw(u8ptr(pre1+pre2)); \
\
	double sc10 = f64ptr(screg1), sc12 = f64ptr(screg1+16); /* cosine/sine for w^2n */ \
	double sc20 = f64ptr(screg2), sc22 = f64ptr(screg2+16); \
	vec2f64 sc30 = xptr(screg3); /* cosine/sine for w^n */ \
	vec2f64 sc32 = -xptr(screg3+16); /* cosine/sine for w^5n (negative) */ \
	vec2f64 tmp; \
\
	/* [0]: Do the eight-reals part */ \
	/* [1]: Do the four complex part */ \
\
	vec2f64 sc; \
	sc[0] = sc12; sc[1] = sc22; \
	xmm0 = M3 * sc + M4;	/* A2 = R2 * cosine/sine + I2 */ \
	xmm3 = M4 * sc - M3;	/* B2 = I2 * cosine/sine - R2 */ \
\
	sc[0] = sc32[1]; sc[1] = sc12; \
	xmm1 = M7 * sc - M8;	/* A4 = R4 * cosine/sine (+,-) I4 */ \
	xmm7 = M8 * sc + M7;	/* B4 = I4 * cosine/sine (-,+) R4 */ \
\
	sc[0] = sc30[1]; sc[1] = sc12; \
	xmm4 = M5 * sc + M6;	/* A3 = R3 * cosine/sine + I3 */ \
	xmm5 = M6 * sc - M5;	/* B3 = I3 * cosine/sine - R3 */ \
\
	sc[0] = sc32[0]; sc[1] = sc10; \
	xmm1 *= sc;				/* new R7,A4 = A4 * sine */ \
	xmm7 *= sc;				/* new R8,B4 = B4 * sine */ \
	sc[0] = sc30[0]; sc[1] = sc10; \
	xmm4 *= sc;				/* new R5,A3 = A3 * sine */ \
	xmm5 *= sc;				/* new R6,B3 = B3 * sine */ \
	sc[0] = sc10; sc[1] = sc20; \
	xmm0 *= sc;				/* new R3,A2 = A2 * sine */ \
	xmm3 *= sc;				/* new R4,B2 = B2 * sine */ \
\
	tmp = xmm5;				/* Copy R6,I3 */ \
	xmm5 -= xmm7;			/* new R8,R4 = R6,I3 - R8,I4 */ \
	xmm7 += tmp;			/* new R7,I3 = R6,I3 + R8,I4 */ \
\
	xmm6 = M1;				/* R1 */ \
	xmm6 -= xmm0;			/* R1,R1 = R1,R2 - R3,R2 (new R3,R2) */ \
	xmm0 += M1;				/* R3,R2 = R1,R1 + R3,R2 (new R1,R1) */ \
\
	xmm2 = M2;				/* I1,R1 */ \
	xmm2 -= xmm3;			/* I1,R2 = I1,R2 - I2,R4 (new I2,R4) */ \
	xmm3 += M2;				/* I2,R4 = I1,R2 + I2,R4 (new I1,R2) */ \
\
	tmp[0] = xmm4[0];		/* Copy R5 */ \
	xmm4[0] -= xmm1[0];		/* new R6 = R5 - R7 */ \
	xmm1[0] += tmp[0];		/* new R5 = R5 + R7 */ \
\
	tmp[0] = xmm4[0];		/* Copy R6 */ \
	xmm4[0] += xmm5[0];		/* R6 = R6 + R8 */ \
	xmm5[0] -= tmp[0];		/* R8 = R8 - R6 */ \
\
	xmm4[0] *= XMM_SQRTHALF1;	/* R6 = R6 * square root of 1/2 */ \
	xmm5[0] *= XMM_SQRTHALF1;	/* R8 = R8 * square root of 1/2 */ \
\
	tmp[1] = xmm1[1];		/* Copy R4 */ \
	xmm1[1] -= xmm4[1];		/* R4 = R4 - R3 (new I4) */ \
	xmm4[1] += tmp[1];		/* R3 = R4 + R3 (new R3) */ \
\
	tmp[1] = xmm0[1];		/* Copy R1 */ \
	xmm0[1] -= xmm4[1];		/* R1 = R1 - R3 (final R3) */ \
	xmm4[1] += tmp[1];		/* R3 = R1 + R3 (final R1) */ \
\
	tmp[1] = xmm6[1];		/* Copy R2 */ \
	xmm6[1] -= xmm5[1];		/* R2 = R2 - R4 (final R4) */ \
	xmm5[1] += tmp[1];		/* R4 = R2 + R4 (final R2) */ \
\
	tmp[0] = xmm2[0];		/* Copy R4 */ \
	xmm2[0] -= xmm5[0];		/* R4 = R4 - R8 (final R8) */ \
	xmm5[0] += tmp[0];		/* R8 = R4 + R8 (final R4) */ \
\
	tmp[0] = xmm3[0];		/* Copy R2 */ \
	xmm3[0] -= xmm4[0];		/* R2 = R2 - R6 (final R6) */ \
	xmm4[0] += tmp[0];		/* R6 = R2 + R6 (final R2) */ \
\
	vec2f64 tmp1 = {xmm5[0], xmm6[1]}; \
	dst4 = tmp1;			/* Save R4 */ \
\
	tmp[0] = xmm6[0];		/* Copy R3 */ \
	xmm6[0] -= xmm7[0];		/* R3 = R3 - R7 (final R7) */ \
	xmm7[0] += tmp[0];		/* R7 = R3 + R7 (final R3) */ \
\
	tmp[0] = xmm0[0];		/* Copy R1 */ \
	xmm0[0] -= xmm1[0];		/* R1 = R1 - R5 (final R5) */ \
	xmm1[0] += tmp[0];		/* R5 = R1 + R5 (final R1) */ \
\
	tmp[1] = xmm2[1];		/* Copy I2 */ \
	xmm2[1] -= xmm1[1];		/* I2 = I2 - I4 (final I4) */ \
	xmm1[1] += tmp[1];		/* I4 = I2 + I4 (final I2) */ \
\
	tmp[1] = xmm3[1];		/* Copy I1 */ \
	xmm3[1] -= xmm7[1];		/* I1 = I1 - I3 (final I3) */ \
	xmm7[1] += tmp[1];		/* I3 = I1 + I3 (final I1) */ \
\
	xmm6[1] = xmm3[1]; \
	xmm3[1] = xmm1[1]; \
	xmm1[1] = xmm4[1]; \
	xmm4[1] = xmm5[1]; \
	xmm5[1] = xmm0[1]; xmm0[1] = xmm7[1]; xmm7[1] = xmm5[1]; \
\
	dst8 = xmm2;		/* Save R8 */ \
}
#else
#define r4_h8r_h4c_2sc_djbunfft_mem(mem1,mem2,mem3,mem4,mem5,mem6,mem7,mem8, screg1,screg2,screg3, pre1,pre2,dst4,dst8) { \
	double f0,f1,f2,f3,f4,f5,f6,f7; \
\
	/* Do the four complex part */ \
\
	f7 = f64ptr2(screg2,16);	/* cosine/sine */ \
	f7 *= mem3[1];		/* A2 = R2 * cosine/sine */ \
	f6 = mem4[1];		/* I2 */ \
	f7 += f6;		/* A2 = A2 + I2 */ \
	f6 *= f64ptr2(screg2,16);	/* B2 = I2 * cosine/sine */ \
\
	f3 = f64ptr2(screg1,16);	/* cosine/sine */ \
	f4 = f3;		/* cosine/sine */ \
	f3 *= mem7[1];		/* A4 = R4 * cosine/sine */ \
	f0 = mem8[1];		/* I4 */ \
	f3 -= f0;		/* A4 = A4 - I4 */ \
	f0 *= f4;		/* B4 = I4 * cosine/sine */ \
	f6 -= mem3[1];		/* B2 = B2 - R2 */ \
\
	f1 = mem5[1];		/* R3 */ \
	f1 *= f4;		/* A3 = R3 * cosine/sine */ \
	f4 *= mem6[1];		/* B3 = I3 * cosine/sine */ \
\
	f0 += mem7[1];		/* B4 = B4 + R4 */ \
	f1 += mem6[1];		/* A3 = A3 + I3 */ \
	f4 -= mem5[1];		/* B3 = B3 - R3 */ \
\
	f7 *= f64ptr(screg2);	/* A2 = A2 * sine (new R2) */ \
	f6 *= f64ptr(screg2);	/* B2 = B2 * sine (new I2) */ \
	f2 = f64ptr(screg1);	/* Sine */ \
	f3 *= f2;		/* A4 = A4 * sine (new R4) */ \
	f0 *= f2;		/* B4 = B4 * sine (new I4) */ \
	f1 *= f2;		/* A3 = A3 * sine (new R3) */ \
	f4 *= f2;		/* B3 = B3 * sine (new I3) */ \
\
	f5 = mem1[1];		/* R1 */ \
	f5 -= f7;		/* R1 = R1 - R2 (new R2) */ \
	f7 += mem1[1];		/* R2 = R1 + R2 (new R1) */ \
\
	f2 = f3;		/* Copy R4 */ \
	f3 -= f1;		/* R4 = R4 - R3 (new I4) */ \
	f1 += f2;		/* R3 = R4 + R3 (new R3) */ \
\
	f2 = f4;		/* Copy I3 */ \
	f4 -= f0;		/* I3 = I3 - I4 (new R4) */ \
	f0 += f2;		/* I4 = I3 + I4 (new I3) */ \
\
	f2 = f7;		/* Copy R1 */ \
	f7 -= f1;		/* R1 = R1 - R3 (final R3) */ \
	f1 += f2;		/* R3 = R1 + R3 (final R1) */ \
\
	f2 = f5;		/* Copy R2 */ \
	f5 -= f4;		/* R2 = R2 - R4 (final R4) */ \
	f4 += f2;		/* R4 = R2 + R4 (final R2) */ \
\
	f2 = mem2[1];		/* I1 */ \
	f2 -= f6;		/* I1 = I1 - I2 (new I2) */ \
	f6 += mem2[1];		/* I2 = I1 + I2 (new I1) */ \
\
	dst4[1] = f5; \
\
	f5 = f2;		/* Copy I2 */ \
	f2 -= f3;		/* I2 = I2 - I4 (final I4) */ \
	f3 += f5;		/* I4 = I2 + I4 (final I2) */ \
\
	f5 = f6;		/* Copy I1 */ \
	f6 -= f0;		/* I1 = I1 - I3 (final I3) */ \
	f0 += f5;		/* I3 = I1 + I3 (final I1) */ \
\
	dst8[1] = f2; \
\
;	/* Do the eight-reals part */ \
\
	xmm0[1] = f0; \
	xmm1[1] = f1; \
	xmm3[1] = f3; \
	xmm4[1] = f4; \
	xmm6[1] = f6; \
	xmm7[1] = f7; \
\
	f1 = f64ptr2(screg3,24);	/* cosine/sine for w^5n */ \
	f1 *= mem7[0];		/* A4 = R4 * cosine/sine */ \
	f7 = mem8[0];		/* I4 */ \
	f1 += f7;		/* A4 = A4 + I4 */ \
	f7 *= f64ptr2(screg3,24);	/* B4 = I4 * cosine/sine */ \
	f7 -= mem7[0];		/* B4 = B4 - R4 */ \
\
	f4 = f64ptr2(screg3,8);	/* cosine/sine for w^n */ \
	f4 *= mem5[0];		/* A3 = R3 * cosine/sine for w^n */ \
	f5 = mem6[0];		/* I3 */ \
	f4 += f5;		/* A3 = A3 + I3 */ \
	f5 *= f64ptr2(screg3,8);	/* B3 = I3 * cosine/sine */ \
	f5 -= mem5[0];		/* B3 = B3 - R3 */ \
\
	f0 = f64ptr2(screg1,16);	/* cosine/sine for w^2n */ \
	f0 *= mem3[0];		/* A2 = R2 * cosine/sine */ \
	f3 = mem4[0];		/* I2 */ \
	f0 += f3;		/* A2 = A2 + I2 */ \
	f3 *= f64ptr2(screg1,16);	/* B2 = I2 * cosine/sine */ \
	f3 -= mem3[0];		/* B2 = B2 - R2 */ \
\
	xprefetchw(u8ptr(pre1)); \
\
	f1 *= f64ptr2(screg3,16);	/* new R7 = A4 * sine */ \
	f7 *= f64ptr2(screg3,16);	/* new R8 = B4 * sine */ \
	f4 *= f64ptr(screg3);	/* new R5 = A3 * sine */ \
	f5 *= f64ptr(screg3);	/* new R6 = B3 * sine */ \
	f0 *= f64ptr(screg1);	/* new R3 = A2 * sine */ \
	f3 *= f64ptr(screg1);	/* new R4 = B2 * sine */ \
\
	xprefetchw(u8ptr(pre1+pre2)); \
\
	f6 = f4;		/* Copy R5 */ \
	f4 -= f1;		/* new R6 = R5 - R7 */ \
	f1 += f6;		/* new R5 = R5 + R7 */ \
\
	f2 = f5;		/* Copy R6 */ \
	f5 -= f7;		/* new R8 = R6 - R8 */ \
	f7 += f2;		/* new R7 = R6 + R8 */ \
\
	f6 = f4;		/* Copy R6 */ \
	f4 += f5;		/* R6 = R6 + R8 */ \
	f5 -= f6;		/* R8 = R8 - R6 */ \
\
	/
	/*vec2f64 t51; \
	t51[0] = f4; t51[1] = f5; \
	t51 *= XMM_SQRTHALF2; \
	f4 = t51[0]; f5 = t51[1];*/ \
	f4 *= XMM_SQRTHALF1;	/* R6 = R6 * square root of 1/2 */ \
	f5 *= XMM_SQRTHALF1;	/* R8 = R8 * square root of 1/2 */ \
\
	f2 = mem2[0];		/* R2 */ \
	f2 -= f3;		/* R2 = R2 - R4 (new R4) */ \
	f3 += mem2[0];		/* R4 = R2 + R4 (new R2) */ \
\
	f6 = f2;		/* Copy R4 */ \
	f2 -= f5;		/* R4 = R4 - R8 (final R8) */ \
	f5 += f6;		/* R8 = R4 + R8 (final R4) */ \
\
	f6 = f3;		/* Copy R2 */ \
	f3 -= f4;		/* R2 = R2 - R6 (final R6) */ \
	f4 += f6;		/* R6 = R2 + R6 (final R2) */ \
\
	f6 = mem1[0];		/* R1 */ \
	f6 -= f0;		/* R1 = R1 - R3 (new R3) */ \
	f0 += mem1[0];		/* R3 = R1 + R3 (new R1) */ \
\
	dst8[0] = f2;		/* Save R8 */ \
	dst4[0] = f5;		/* Save R4 */ \
\
	f2 = f6;		/* Copy R3 */ \
	f6 -= f7;		/* R3 = R3 - R7 (final R7) */ \
	f7 += f2;		/* R7 = R3 + R7 (final R3) */ \
\
	f5 = f0;		/* Copy R1 */ \
	f0 -= f1;		/* R1 = R1 - R5 (final R5) */ \
	f1 += f5;		/* R5 = R1 + R5 (final R1) */ \
	xmm0[0] = f0; \
	xmm1[0] = f1; \
	xmm3[0] = f3; \
	xmm4[0] = f4; \
	xmm6[0] = f6; \
	xmm7[0] = f7; \
}
#endif

/* */ \
/* ********************************* half-eight-reals-four-complex-fft-with-square variants *************************************** */ \
/* */ \
\
/* Macro to do an eight_reals_fft and three four_complex_djbfft in the final levels of pass 2. */ \
/* The eight-reals macro and one of the four-complex only use half the XMM  */ \
/* register.  This isn't very efficient, but this macro is called only once. */ \
\

#define r4_h4cl_eight_reals_four_complex_fft_final(srcreg,srcinc,d1,d2) { \
	vec2f64 *sd0 = (vec2f64*)(srcreg); \
	vec2f64 *sd1 = (vec2f64*)(srcreg+d1); \
	vec2f64 *sd2 = (vec2f64*)(srcreg+d2); \
	vec2f64 *sd3 = (vec2f64*)(srcreg+d2+d1); \
	r4_h8r_h4c_simple_fft_mem(sd0[0],sd1[0],sd2[0],sd3[0],sd0[1],sd1[1],sd2[1],sd3[1], sd0[0],sd0[1]); \
/*	xstore	[srcreg], xmm0		;; Save R1 */ \
/*	xstore	[srcreg+16], xmm0	;; Save I1 */ \
	xmm3 = sd0[2];	/* R1 */ \
	xmm7 = sd0[3];	/* R5 */ \
	sd0[2] = xmm2;	/* Save R2 */ \
	sd0[3] = xmm6;	/* Save I2 */ \
	xmm2 = sd1[2];	/* R2 */ \
	xmm6 = sd1[3];	/* R6 */ \
	sd1[0] = xmm0;	/* Save R3 */ \
	sd1[1] = xmm1;	/* Save I3 */ \
	sd1[2] = xmm5;	/* Save R4 */ \
	sd1[3] = xmm4;	/* Save I4 */ \
	r4_x4c_simple_fft_partial_mem(xmm3,xmm2,xmm0,xmm5,xmm7,xmm6,xmm1,xmm4, sd2[2],sd3[2],sd2[3],sd3[3],srcreg+srcinc+d2,d1); \
	sd2[0] = xmm5;	/* Save R1 */ \
	sd2[1] = xmm4;	/* Save R2 */ \
	sd2[2] = xmm0;	/* Save R3 */ \
	sd2[3] = xmm1;	/* Save R4 */ \
	sd3[0] = xmm3;	/* Save R5 */ \
	sd3[1] = xmm2;	/* Save R6 */ \
	sd3[2] = xmm6;	/* Save R7 */ \
	sd3[3] = xmm7;	/* Save R8 */ \
	srcreg += srcinc; \
}


#define r4_h4cl_eight_reals_four_complex_with_square(srcreg,srcinc,d1,d2) { \
	vec2f64 *sd0 = (vec2f64*)(srcreg); \
	vec2f64 *sd1 = (vec2f64*)(srcreg+d1); \
	vec2f64 *sd2 = (vec2f64*)(srcreg+d2); \
	vec2f64 *sd3 = (vec2f64*)(srcreg+d2+d1); \
	r4_h8r_h4c_simple_fft_mem(sd0[0],sd1[0],sd2[0],sd3[0],sd0[1],sd1[1],sd2[1],sd3[1], sd0[0],sd0[1]); \
	xp_complex_square(xmm2, xmm6, xmm7);	/* Square R2, I2 */ \
	sd1[0] = xmm2; \
	sd1[1] = xmm6; \
	xp_complex_square(xmm0, xmm1, xmm7);	/* Square R3, I3 */ \
	sd2[0] = xmm0; \
	sd2[1] = xmm1; \
	xp_complex_square(xmm5, xmm4, xmm7);	/* Square R4, I4 */ \
	sd3[0] = xmm5; \
	sd3[1] = xmm4; \
	xmm0[0] = f64ptr(srcreg);		/* R1 */ \
	xmm1[0] = f64ptr(srcreg+16);		/* I1 */ \
	xmm0[0] *= xmm0[0];			/* Square R1, I1 low (actually they are both real) */ \
	xmm1[0] *= xmm1[0]; \
	f64ptr(srcreg-16) = xmm0[0];		/* Save product of sum of FFT values */ \
	f64ptr(srcreg) = xmm0[0]; \
	f64ptr(srcreg+16) = xmm1[0]; \
	xmm0[0] = f64ptr(srcreg+8);		/* R1 */ \
	xmm1[0] = f64ptr(srcreg+24);		/* I1 */ \
	xs_complex_square(xmm0[0], xmm1[0], xmm2[0]);	/* Square R1, I1 high */ \
	f64ptr(srcreg+8) = xmm0[0]; \
	f64ptr(srcreg+24) = xmm1[0]; \
	r4_h8r_h4c_simple_unfft_mem(sd0[0],sd0[1],sd1[0],sd1[1],sd2[0],sd2[1],sd3[0],sd3[1], srcreg+srcinc,d1,sd1[0]); \
	sd0[0] = xmm1;		/* Save R1 */ \
	sd0[1] = xmm7;		/* Save R3 */ \
	xmm1 = sd0[2];		/* R1 */ \
	xmm7 = sd0[3];		/* R5 */ \
	sd0[2] = xmm0;		/* Save R5 */ \
	sd0[3] = xmm6;		/* Save R7 */ \
	xmm0 = sd1[2];		/* R2 */ \
	xmm6 = sd1[3];		/* R6 */ \
/*	xstore	[srcreg+d1], xmm4		;; Save R2 */ \
	sd1[1] = xmm5;		/* Save R4 */ \
	sd1[2] = xmm3;		/* Save R6 */ \
	sd1[3] = xmm2;		/* Save R8 */ \
	r4_x4c_simple_fft_partial_mem(xmm1,xmm0,xmm2,xmm3,xmm7,xmm6,xmm5,xmm4,sd2[2],sd3[2],sd2[3],sd3[3],srcreg+srcinc+d2,d1); \
	sd2[0] = xmm1;		/* R5 */ \
	xp_complex_square(xmm3, xmm4, xmm1);	/* Square R1, R2 */ \
	xp_complex_square(xmm2, xmm5, xmm1);	/* Square R3, R4 */ \
	xp_complex_square(xmm6, xmm7, xmm1);	/* Square R7, R8 */ \
	xmm1 = sd2[0];		/* R5 */ \
	sd2[0] = xmm7;		/* R8 */ \
	xp_complex_square(xmm1, xmm0, xmm7);	/* Square R5, R6 */ \
	r4_x4c_simple_unfft(xmm3,xmm4,xmm2,xmm5,xmm1,xmm0,xmm6,xmm7,sd2[0],srcreg+srcinc+d2,d1,sd2[0]); \
/*	xstore	[srcreg+d2], xmm1		;; Save R1 */ \
	sd2[1] = xmm2;		/* Save R3 */ \
	sd2[2] = xmm7;		/* Save R5 */ \
	sd2[3] = xmm5;		/* Save R7 */ \
	sd3[0] = xmm0;		/* Save R2 */ \
	sd3[1] = xmm3;		/* Save R4 */ \
	sd3[2] = xmm6;		/* Save R6 */ \
	sd3[3] = xmm4;		/* Save R8 */ \
	srcreg += srcinc; \
}


#define r4_h4cl_eight_reals_four_complex_with_mult(srcreg,srcinc,d1,d2) { \
	uintptr_t src_rbp = srcreg+rbp; \
	vec2f64 *sd0 = (vec2f64*)(srcreg); \
	vec2f64 *sd1 = (vec2f64*)(srcreg+d1); \
	vec2f64 *sd2 = (vec2f64*)(srcreg+d2); \
	vec2f64 *sd3 = (vec2f64*)(srcreg+d2+d1); \
	vec2f64 *pd0 = (vec2f64*)(src_rbp); \
	vec2f64 *pd1 = (vec2f64*)(src_rbp+d1); \
	vec2f64 *pd2 = (vec2f64*)(src_rbp+d2); \
	vec2f64 *pd3 = (vec2f64*)(src_rbp+d2+d1); \
	r4_h8r_h4c_simple_fft_mem(sd0[0],sd1[0],sd2[0],sd3[0],sd0[1],sd1[1],sd2[1],sd3[1], sd0[0],sd0[1]); \
	xmm3 = sd0[0]; \
	xmm7 = sd0[1]; \
	xp4c_mulf_r(xmm3,xmm7,xmm2,xmm6,xmm0,xmm1,xmm5,xmm4, pd0[0],pd0[1],pd0[2],pd0[3],pd1[0],pd1[1],pd1[2],pd1[3]); \
	sd1[0] = xmm2; \
	sd1[1] = xmm6; \
	sd2[0] = xmm0; \
	sd2[1] = xmm1; \
	sd3[0] = xmm5; \
	sd3[1] = xmm4; \
	xmm3[0] = f64ptr(srcreg); \
	xmm3[0] *= f64ptr(src_rbp);		/* Multiply R1 low */ \
	f64ptr(srcreg-16) = xmm3[0];		/* Save product of sum of FFT values */ \
	xmm7[0] = f64ptr(srcreg+16); \
	xmm7[0] *= f64ptr(src_rbp+16);	/* Multiply I1 low (actually a real value) */ \
	sd0[0] = xmm3; \
	sd0[1] = xmm7; \
	r4_h8r_h4c_simple_unfft_mem(sd0[0],sd0[1],sd1[0],sd1[1],sd2[0],sd2[1],sd3[0],sd3[1], srcreg+srcinc,d1,sd1[0]); \
	sd0[0] = xmm1;			/* Save R1 */ \
	sd0[1] = xmm7;		/* Save R3 */ \
	xmm1 = sd0[2];		/* R1 */ \
	xmm7 = sd0[3];		/* R5 */ \
	sd0[2] = xmm0;		/* Save R5 */ \
	sd0[3] = xmm6;		/* Save R7 */ \
	xmm0 = sd1[2];		/* R2 */ \
	xmm6 = sd1[3];		/* R6 */ \
/*	xstore	[srcreg+d1], xmm4		;; Save R2 */ \
	sd1[1] = xmm5;		/* Save R4 */ \
	sd1[2] = xmm3;		/* Save R6 */ \
	sd1[3] = xmm2;		/* Save R8 */ \
	r4_x4c_simple_fft_partial_mem(xmm1,xmm0,xmm2,xmm3,xmm7,xmm6,xmm5,xmm4,sd2[2],sd3[2],sd2[3],sd3[3],srcreg+srcinc+d2,d1); \
	xp4c_mulf_r(xmm3,xmm4,xmm2,xmm5,xmm1,xmm0,xmm6,xmm7, pd2[0],pd2[1],pd2[2],pd2[3],pd3[0],pd3[1],pd3[2],pd3[3]); \
	sd2[0] = xmm7; \
	r4_x4c_simple_unfft(xmm3,xmm4,xmm2,xmm5,xmm1,xmm0,xmm6,xmm7,sd2[0],srcreg+srcinc+d2,d1,sd2[0]); \
/*	xstore	[srcreg+d2], xmm1		;; Save R1 */ \
	sd2[1] = xmm2;		/* Save R3 */ \
	sd2[2] = xmm7;		/* Save R5 */ \
	sd2[3] = xmm5;		/* Save R7 */ \
	sd3[0] = xmm0;		/* Save R2 */ \
	sd3[1] = xmm3;		/* Save R4 */ \
	sd3[2] = xmm6;		/* Save R6 */ \
	sd3[3] = xmm4;		/* Save R8 */ \
	srcreg += srcinc; \
}


#define r4_h4cl_eight_reals_four_complex_with_mulf(srcreg,srcinc,d1,d2) { \
	uintptr_t src_rbx = srcreg+rbx; \
	uintptr_t src_rbp = srcreg+rbp; \
	vec2f64 *sd0 = (vec2f64*)(srcreg); \
	vec2f64 *sd1 = (vec2f64*)(srcreg+d1); \
	vec2f64 *sd2 = (vec2f64*)(srcreg+d2); \
	vec2f64 *sd3 = (vec2f64*)(srcreg+d2+d1); \
	vec2f64 *bd0 = (vec2f64*)(src_rbx); \
	vec2f64 *bd1 = (vec2f64*)(src_rbx+d1); \
	vec2f64 *bd2 = (vec2f64*)(src_rbx+d2); \
	vec2f64 *bd3 = (vec2f64*)(src_rbx+d2+d1); \
	vec2f64 *pd0 = (vec2f64*)(src_rbp); \
	vec2f64 *pd1 = (vec2f64*)(src_rbp+d1); \
	vec2f64 *pd2 = (vec2f64*)(src_rbp+d2); \
	vec2f64 *pd3 = (vec2f64*)(src_rbp+d2+d1); \
	xmm3 = bd0[0];		/* R1 */ \
	xmm7 = bd0[1];		/* R2 */ \
	xmm2 = bd0[2];		/* R3 */ \
	xmm6 = bd0[3];		/* R4 */ \
	xmm0 = bd1[0];		/* R5 */ \
	xmm1 = bd1[1];		/* R6 */ \
	xmm5 = bd1[2];		/* R7 */ \
	xmm4 = bd1[3];		/* R8 */ \
	xp4c_mulf_r(xmm3,xmm7,xmm2,xmm6,xmm0,xmm1,xmm5,xmm4, pd0[0],pd0[1],pd0[2],pd0[3],pd1[0],pd1[1],pd1[2],pd1[3]); \
	sd0[2] = xmm2; \
	sd0[3] = xmm6; \
	sd1[0] = xmm0; \
	sd1[1] = xmm1; \
	sd1[2] = xmm5; \
	sd1[3] = xmm4; \
	xmm3[0] = f64ptr(src_rbx); \
	xmm3[0] *= f64ptr(src_rbp);		/* Multiply R1 low */ \
	f64ptr(srcreg-16) = xmm3[0];		/* Save product of sum of FFT values */ \
	xmm7[0] = f64ptr(src_rbx+16); \
	xmm7[0] *= f64ptr(src_rbp+16);	/* Multiply I1 low (actually a real value) */ \
	sd0[0] = xmm3; \
	sd0[1] = xmm7; \
	r4_h8r_h4c_simple_unfft_mem(sd0[0],sd0[1],sd0[2],sd0[3],sd1[0],sd1[1],sd1[2],sd1[3],srcreg+srcinc,d1,sd1[0]); \
	sd0[0] = xmm1;			/* Save R1 */ \
	sd0[1] = xmm7;		/* Save R3 */ \
	sd0[2] = xmm0;		/* Save R5 */ \
	sd0[3] = xmm6;		/* Save R7 */ \
/*	xstore	[srcreg+d1], xmm4		;; Save R2 */ \
	sd1[1] = xmm5;		/* Save R4 */ \
	sd1[2] = xmm3;		/* Save R6 */ \
	sd1[3] = xmm2;		/* Save R8 */ \
	xmm3 = bd2[0];		/* R1 */ \
	xmm7 = bd2[1];		/* R2 */ \
	xmm1 = bd2[2];		/* R3 */ \
	xmm5 = bd2[3];		/* R4 */ \
	xmm0 = bd3[0];		/* R5 */ \
	xmm4 = bd3[1];		/* R6 */ \
	xmm6 = bd3[2];		/* R7 */ \
	xmm2 = bd3[3];		/* R8 */ \
	xp4c_mulf_r(xmm3,xmm7,xmm1,xmm5,xmm0,xmm4,xmm6,xmm2, pd2[0],pd2[1],pd2[2],pd2[3],pd3[0],pd3[1],pd3[2],pd3[3]); \
	sd2[0] = xmm2; \
	r4_x4c_simple_unfft(xmm3,xmm7,xmm1,xmm5,xmm0,xmm4,xmm6,xmm2,sd2[0],srcreg+srcinc+d2,d1,sd2[0]); \
/*	xstore	[srcreg+d2], xmm0		;; Save R1 */ \
	sd2[1] = xmm1;		/* Save R3 */ \
	sd2[2] = xmm2;		/* Save R5 */ \
	sd2[3] = xmm5;		/* Save R7 */ \
	sd3[0] = xmm4;		/* Save R2 */ \
	sd3[1] = xmm3;		/* Save R4 */ \
	sd3[2] = xmm6;		/* Save R6 */ \
	sd3[3] = xmm7;		/* Save R8 */ \
	srcreg += srcinc; \
}

/* Does an r4_x8r_fft_mem on the low word of the xmm register */ \
/* Does an r4_x4c_djbfft_mem on the high word of the xmm register */ \

#if 1
#define r4_h8r_h4c_simple_fft_mem(mem1,mem2,mem3,mem4,mem5,mem6,mem7,mem8,dst1,dst2) { \
	vec2f64 M1=mem1, M2=mem2, M3=mem3, M4=mem4, M5=mem5, M6=mem6, M7=mem7, M8=mem8; \
\
	/*[0]: Do the eight reals part */ \
	/*[1]: Do the four-complex part */ \
\
	xmm1[0] = M4[0];		/* R4 */ \
	xmm3[0] = M8[0];		/* R8 */ \
	xmm1[0] -= xmm3[0];		/* new R8 = R4 - R8 */ \
	xmm3[0] += M4[0];		/* new R4 = R4 + R8 */ \
	xmm0[0] = M2[0];		/* R2 */ \
	xmm6[0] = M6[0];		/* R6 */ \
	xmm0[0] -= xmm6[0];		/* new R6 = R2 - R6 */ \
	xmm6[0] += M2[0];		/* new R2 = R2 + R6 */ \
\
	/*vec2xmm6[0]4 t51; \
	t51[0] = xmm1[0]; t51[1] = xmm0[0]; \
	t51 *= XMM_SQRTHALF2; \
	xmm1[0] = t51[0]; xmm0[0] = t51[1];*/ \
	xmm1[0] *= XMM_SQRTHALF1;	/* R8 = R8 * square root */ \
	xmm0[0] *= XMM_SQRTHALF1;	/* R6 = R6 * square root */ \
\
	xmm5[0] = M1[0];		/* R1 */ \
	xmm2[0] = M5[0];		/* R5 */ \
	xmm5[0] -= xmm2[0];		/* new R5 = R1 - R5 */ \
	xmm2[0] += M1[0];		/* new R1 = R1 + R5 */ \
\
	xmm4[0] = xmm6[0];		/* Copy R2 */ \
	xmm6[0] -= xmm3[0];		/* R2 = R2 - R4 (final I2) */ \
	xmm3[0] += xmm4[0];		/* R4 = R2 + R4 (new I1) */ \
\
	xmm7[0] = xmm0[0];		/* Copy R6 */ \
	xmm0[0] -= xmm1[0];		/* R6 = R6 - R8 (Real part) */ \
	xmm1[0] += xmm7[0];		/* R8 = R6 + R8 (Imaginary part) */ \
\
	xmm7[0] = M7[0];		/* R7 */ \
	xmm7[0] += M3[0];		/* new R3 = R3 + R7 */ \
\
	xmm4[0] = xmm2[0];		/* Copy R1 */ \
	xmm2[0] -= xmm7[0];		/* R1 = R1 - R3 (final R2) */ \
	xmm7[0] += xmm4[0];		/* R3 = R1 + R3 (new R1) */ \
\
	xmm4[0] = xmm7[0];		/* Copy R1 */ \
	xmm7[0] -= xmm3[0];		/* R1 = R1 - R2 (final I1) */ \
	xmm3[0] += xmm4[0];		/* R2 = R1 + R2 (final R1) */ \
\
	xmm4[0] = M3[0];		/* R3 */ \
	xmm4[0] -= M7[0];		/* new R7 = R3 - R7 */ \
\
\
	xmm5[1] = M1[1];		/* R1 */ \
	xmm2[1] = M3[1];		/* R3 */ \
	xmm2[1] += xmm5[1];		/* R3 = R1 + R3 (new R1) */ \
	xmm5[1] -= M3[1];		/* R1 = R1 - R3 (new R3) */ \
\
	xmm1[1] = M2[1];		/* R2 */ \
	xmm3[1] = M4[1];		/* R4 */ \
	xmm3[1] += xmm1[1];		/* R4 = R2 + R4 (new R2) */ \
	xmm1[1] -= M4[1];		/* R2 = R2 - R4 (new R4) */ \
\
	xmm6[1] = xmm2[1]; \
	xmm2[1] -= xmm3[1];		/* R1 = R1 - R2 (final R2) */ \
	xmm3[1] += xmm6[1];		/* R2 = R1 + R2 (final R1) */ \
\
	xmm0[1] = M6[1];		/* I2 */ \
	xmm7[1] = M8[1];		/* I4 */ \
	xmm0[1] = xmm7[1] - xmm0[1]; /* I2 = I2 - I4 (new negative I4) */ \
	xmm7[1] += M6[1];		/* I4 = I2 + I4 (new I2) */ \
\
	xmm4[1] = M5[1];		/* I1 */ \
	xmm6[1] = M7[1];		/* I3 */ \
	xmm4[1] -= xmm6[1];		/* I1 = I1 - I3 (new I3) */ \
	xmm6[1] += M5[1];		/* I3 = I1 + I3 (new I1) */ \
\
	dst1 = xmm3;			/* Save R1 R1 */ \
\
	xmm3[1] = xmm6[1]; \
	xmm6[1] -= xmm7[1];		/* I1 = I1 - I2 (final I2) */ \
	xmm7[1] += xmm3[1];		/* I2 = I1 + I2 (final I1) */ \
\
	xmm3 = xmm5;			/* Copy R5 R3 */ \
	/* xmm5[1]: R3 */  \
	/* xmm0[1]: negative I4 */  \
	xmm5 = xmm3 - xmm0;		/* R5,I4 = R5,R3 - R6,-I4 (final R4,R4) */ \
	xmm0 = xmm0 + xmm3;		/* R6,R3 = R5,R3 + R6,-I4 (final R3,R3) */ \
\
	xmm3 = xmm4;			/* Copy R8I4 */ \
	xmm4 -= xmm1;			/* R7I3 = R7I3 - R8R4 (final I4I4) */ \
	xmm1 += xmm3;			/* R8R4 = R7I3 + R8R4 (final I3I3) */ \
\
	dst2 = xmm7;			/* Save I1 I1 */ \
}
#else
#define r4_h8r_h4c_simple_fft_mem(mem1,mem2,mem3,mem4,mem5,mem6,mem7,mem8,dst1,dst2) { \
	double f0,f1,f2,f3,f4,f5,f6,f7; \
\
	/* Do the four-complex part */ \
\
	f0 = mem1[1];		/* R1 */ \
	f2 = mem3[1];		/* R3 */ \
	f2 += f0;		/* R3 = R1 + R3 (new R1) */ \
	f0 -= mem3[1];		/* R1 = R1 - R3 (new R3) */ \
\
	f1 = mem2[1];		/* R2 */ \
	f3 = mem4[1];		/* R4 */ \
	f3 += f1;		/* R4 = R2 + R4 (new R2) */ \
	f1 -= mem4[1];		/* R2 = R2 - R4 (new R4) */ \
\
	f6 = f2; \
	f2 -= f3;		/* R1 = R1 - R2 (final R2) */ \
	f3 += f6;		/* R2 = R1 + R2 (final R1) */ \
\
	dst1[1] = f3;		/* Save R1 */ \
\
	f5 = mem6[1];		/* I2 */ \
	f7 = mem8[1];		/* I4 */ \
	f5 -= f7;		/* I2 = I2 - I4 (new I4) */ \
	f7 += mem6[1];		/* I4 = I2 + I4 (new I2) */ \
\
	f4 = mem5[1];		/* I1 */ \
	f6 = mem7[1];		/* I3 */ \
	f4 -= f6;		/* I1 = I1 - I3 (new I3) */ \
	f6 += mem5[1];		/* I3 = I1 + I3 (new I1) */ \
\
	f3 = f0; \
	f0 -= f5;		/* R3 = R3 - I4 (final R3) */ \
	f5 += f3;		/* I4 = R3 + I4 (final R4) */ \
\
	f3 = f1; \
	f1 += f4;		/* R4 = I3 + R4 (final I3) */ \
	f4 -= f3;		/* I3 = I3 - R4 (final I4) */ \
\
	f3 = f6; \
	f6 -= f7;		/* I1 = I1 - I2 (final I2) */ \
	f7 += f3;		/* I2 = I1 + I2 (final I1) */ \
\
	dst2[1] = f7;		/* Save I1 */ \
\
;	/* Do the eight reals part */ \
\
	xmm2[1] = f2;		/* Copy to high part of XMM register */ \
	xmm6[1] = f6; \
	xmm0[1] = f0; \
	xmm1[1] = f1; \
	xmm5[1] = f5; \
	xmm4[1] = f4; \
\
	f1 = mem4[0];		/* R4 */ \
	f7 = mem8[0];		/* R8 */ \
	f1 -= f7;		/* new R8 = R4 - R8 */ \
	f7 += mem4[0];		/* new R4 = R4 + R8 */ \
	f0 = mem2[0];		/* R2 */ \
	f6 = mem6[0];		/* R6 */ \
	f0 -= f6;		/* new R6 = R2 - R6 */ \
	f6 += mem2[0];		/* new R2 = R2 + R6 */ \
\
	/*vec2f64 t51; \
	t51[0] = f1; t51[1] = f0; \
	t51 *= XMM_SQRTHALF2; \
	f1 = t51[0]; f0 = t51[1];*/ \
	f1 *= XMM_SQRTHALF1;	/* R8 = R8 * square root */ \
	f0 *= XMM_SQRTHALF1;	/* R6 = R6 * square root */ \
\
	f5 = mem1[0];		/* R1 */ \
	f2 = mem5[0];		/* R5 */ \
	f5 -= f2;		/* new R5 = R1 - R5 */ \
	f2 += mem1[0];		/* new R1 = R1 + R5 */ \
\
	f4 = f6;		/* Copy R2 */ \
	f6 -= f7;		/* R2 = R2 - R4 (final I2) */ \
	f7 += f4;		/* R4 = R2 + R4 (new I1) */ \
\
	f3 = f0;		/* Copy R6 */ \
	f0 -= f1;		/* R6 = R6 - R8 (Real part) */ \
	f1 += f3;		/* R8 = R6 + R8 (Imaginary part) */ \
\
	f3 = mem7[0];		/* R7 */ \
	f3 += mem3[0];		/* new R3 = R3 + R7 */ \
\
	f4 = f2;		/* Copy R1 */ \
	f2 -= f3;		/* R1 = R1 - R3 (final R2) */ \
	f3 += f4;		/* R3 = R1 + R3 (new R1) */ \
\
	f4 = f3;		/* Copy R1 */ \
	f3 -= f7;		/* R1 = R1 - R2 (final I1) */ \
	f7 += f4;		/* R2 = R1 + R2 (final R1) */ \
\
	f4 = mem3[0];		/* R3 */ \
	f4 -= mem7[0];		/* new R7 = R3 - R7 */ \
\
	dst1[0] = f7;		/* Save R1 */ \
	dst2[0] = f3;		/* Save I1 */ \
\
	f7 = f5;		/* Copy R5 */ \
	f5 -= f0;		/* R5 = R5 - R6 (final R4) */ \
	f0 += f7;		/* R6 = R5 + R6 (final R3) */ \
\
	f3 = f4;		/* Copy R7 */ \
	f4 -= f1;		/* R7 = R7 - R8 (final I4) */ \
	f1 += f3;		/* R8 = R7 + R8 (final I3) */ \
	xmm0[0] = f0; \
	xmm1[0] = f1; \
	xmm2[0] = f2; \
	xmm4[0] = f4; \
	xmm5[0] = f5; \
	xmm6[0] = f6; \

#endif

#if 1
#define r4_h8r_h4c_simple_unfft_mem(mem1,mem2,mem3,mem4,mem5,mem6,mem7,mem8,pre1,pre2,dst2) { \
	vec2f64 M1=mem1, M2=mem2, M3=mem3, M4=mem4, M5=mem5, M6=mem6, M7=mem7, M8=mem8; \
	xprefetchw(u8ptr(pre1)); \
	xprefetchw(u8ptr(pre1+pre2)); \
\
;	/* Do the four complex part */ \
\
	xmm4 = M5;		/* R5 R3 */ \
	xmm1 = M7;		/* R7 R4 */ \
	xmm4 -= xmm1;		/* R5R3 = R5R3 - R7R4 (new R6 negI4) */ \
	xmm1 += M5;		/* R7R4 = R5R3 + R7R4 (new R3 R3) */ \
	xmm3[1] = xmm4[1]; \
\
	xmm5 = M6;		/* R6 I3 */ \
	xmm7 = M8;		/* R8 I4 */ \
	xmm5 -= xmm7;		/* R6I3 = R6I3 - R8I4 (new R8R4) */ \
	xmm7 += M6;		/* R8I4 = R6I3 + R8I4 (new R7I3) */ \
	xmm4[1] = xmm5[1]; \
	xmm0[1] = xmm7[1]; \
\
	xmm5[1] = M1[1];		/* R1 */ \
	xmm7[1] = M3[1];		/* R2 */ \
	xmm5[1] -= xmm7[1];		/* R1 = R1 - R2 (new R2) */ \
	xmm7[1] += M1[1];		/* R2 = R1 + R2 (new R1) */ \
\
	xmm2[1] = xmm7[1];		/* Copy R1 */ \
	xmm7[1] -= xmm1[1];		/* R1 = R1 - R3 (final R3) */ \
	xmm1[1] += xmm2[1];		/* R3 = R1 + R3 (final R1) */ \
\
\
	xmm2[1] = xmm5[1];		/* Copy R2 */ \
	xmm5[1] -= xmm4[1];		/* R2 = R2 - R4 (final R4) */ \
	xmm4[1] += xmm2[1];		/* R4 = R2 + R4 (final R2) */ \
\
	xmm2[1] = M2[1];		/* I1 */ \
	xmm6[1] = M4[1];		/* I2 */ \
	xmm2[1] -= xmm6[1];		/* I1 = I1 - I2 (new I2) */ \
	xmm6[1] += M2[1];		/* I2 = I1 + I2 (new I1) */ \
\
;	/* Do the eight-reals part */ \
\
	xmm6[0] = xmm4[0];		/* Copy R6 */ \
	xmm4[0] += xmm5[0];		/* R6 = R6 + R8 */ \
	xmm5[0] -= xmm6[0];		/* R8 = R8 - R6 */ \
\
	/*vec2xmm6[0]4 t51; \
	t51[0] = xmm4[0]; t51[1] = xmm5[0]; \
	t51 *= XMM_SQRTHALF2; \
	xmm4[0] = t51[0]; xmm5[0] = t51[1];*/ \
	xmm4[0] *= XMM_SQRTHALF1;	/* R6 = R6 * square root of 1/2 */ \
	xmm5[0] *= XMM_SQRTHALF1;	/* R8 = R8 * square root of 1/2 */ \
\
	xmm2[0] = M1[0];		/* R1 */ \
	xmm6[0] = M2[0];		/* R2 */ \
	xmm2[0] -= xmm6[0];		/* new R2 = R1 - R2 */ \
	xmm2[0] *= 0.5;			/* Mul R2 by HALF */ \
	xmm6[0] += xmm2[0];		/* new R1 = R1 + R2 */ \
\
	xmm3[0] = M4[0];		/* R4 */ \
	xmm3[0] += xmm2[0];		/* R4 = R2 + R4 (new R2) */ \
	xmm2[0] -= M4[0];		/* R2 = R2 - R4 (new R4) */ \
\
	xmm0[0] = xmm2[0];		/* Copy R4 */ \
	xmm2[0] -= xmm5[0];		/* R4 = R4 - R8 (final R8) */ \
	xmm5[0] += xmm0[0];		/* R8 = R4 + R8 (final R4) */ \
\
	xmm0[0] = xmm3[0];		/* Copy R2 */ \
	xmm3[0] -= xmm4[0];		/* R2 = R2 - R6 (final R6) */ \
	xmm4[0] += xmm0[0];		/* R6 = R2 + R6 (final R2) */ \
\
	xmm0[0] = M3[0];		/* R3 */ \
	xmm0[0] += xmm6[0];		/* R3 = R1 + R3 (new R1) */ \
	xmm6[0] -= M3[0];		/* R1 = R1 - R3 (new R3) */ \
\
	dst2 = xmm4;		/* Save R4 */ \
\
	xmm4[0] = xmm6[0];		/* Copy R3 */ \
	xmm6[0] -= xmm7[0];		/* R3 = R3 - R7 (final R7) */ \
	xmm7[0] += xmm4[0];		/* R7 = R3 + R7 (final R3) */ \
\
	xmm4[0] = xmm0[0];		/* Copy R1 */ \
	xmm0[0] -= xmm1[0];		/* R1 = R1 - R5 (final R5) */ \
	xmm1[0] += xmm4[0];		/* R5 = R1 + R5 (final R1) */ \
\
	xmm4[1] = xmm2[1];		/* Copy I2 */ \
	xmm2[1] += xmm3[1];		/* I2 = I2 - negI4 (final I4) */ \
	xmm3[1] = xmm4[1] - xmm3[1]; /* I4 = I2 + negI4 (final I2) */ \
\
	xmm4[1] = xmm6[1];		/* Copy I1 */ \
	xmm6[1] -= xmm0[1];		/* I1 = I1 - I3 (final I3) */ \
	xmm0[1] += xmm4[1];		/* I3 = I1 + I3 (final I1) */ \
}
#else
#define r4_h8r_h4c_simple_unfft_mem(mem1,mem2,mem3,mem4,mem5,mem6,mem7,mem8,pre1,pre2,dst2) { \
	double f0,f1,f2,f3,f4,f5,f6,f7; \
\
;	/* Do the four complex part */ \
\
	f5 = mem1[1];		/* R1 */ \
	f7 = mem3[1];		/* R2 */ \
	f5 -= f7;		/* R1 = R1 - R2 (new R2) */ \
	f7 += mem1[1];		/* R2 = R1 + R2 (new R1) */ \
\
	f3 = mem7[1];		/* R4 */ \
	f1 = mem5[1];		/* R3 */ \
	f3 -= f1;		/* R4 = R4 - R3 (new I4) */ \
	f1 += mem7[1];		/* R3 = R4 + R3 (new R3) */ \
\
	xprefetchw(u8ptr(pre1)); \
\
	f4 = mem6[1];		/* I3 */ \
	f0 = mem8[1];		/* I4 */ \
	f4 -= f0;		/* I3 = I3 - I4 (new R4) */ \
	f0 += mem6[1];		/* I4 = I3 + I4 (new I3) */ \
\
	f2 = f7;		/* Copy R1 */ \
	f7 -= f1;		/* R1 = R1 - R3 (final R3) */ \
	f1 += f2;		/* R3 = R1 + R3 (final R1) */ \
\
	xprefetchw(u8ptr(pre1+pre2)); \
\
	f2 = f5;		/* Copy R2 */ \
	f5 -= f4;		/* R2 = R2 - R4 (final R4) */ \
	f4 += f2;		/* R4 = R2 + R4 (final R2) */ \
\
	f2 = mem2[1];		/* I1 */ \
	f6 = mem4[1];		/* I2 */ \
	f2 -= f6;		/* I1 = I1 - I2 (new I2) */ \
	f6 += mem2[1];		/* I2 = I1 + I2 (new I1) */ \
\
	dst2[1] = f4; \
\
	f4 = f2;		/* Copy I2 */ \
	f2 -= f3;		/* I2 = I2 - I4 (final I4) */ \
	f3 += f4;		/* I4 = I2 + I4 (final I2) */ \
\
	f4 = f6;		/* Copy I1 */ \
	f6 -= f0;		/* I1 = I1 - I3 (final I3) */ \
	f0 += f4;		/* I3 = I1 + I3 (final I1) */ \
\
;	/* Do the eight-reals part */ \
\
	xmm0[1] = f0; \
	xmm1[1] = f1; \
	xmm2[1] = f2; \
	xmm3[1] = f3; \
	xmm5[1] = f5; \
	xmm6[1] = f6; \
	xmm7[1] = f7; \
\
	f4 = mem5[0];		/* R5 */ \
	f1 = mem7[0];		/* R7 */ \
	f4 -= f1;		/* new R6 = R5 - R7 */ \
	f1 += mem5[0];		/* new R5 = R5 + R7 */ \
\
	f5 = mem6[0];		/* R6 */ \
	f7 = mem8[0];		/* R8 */ \
	f5 -= f7;		/* new R8 = R6 - R8 */ \
	f7 += mem6[0];		/* new R7 = R6 + R8 */ \
\
	f6 = f4;		/* Copy R6 */ \
	f4 += f5;		/* R6 = R6 + R8 */ \
	f5 -= f6;		/* R8 = R8 - R6 */ \
\
	/*vec2f64 t51; \
	t51[0] = f4; t51[1] = f5; \
	t51 *= XMM_SQRTHALF2; \
	f4 = t51[0]; f5 = t51[1];*/ \
	f4 *= XMM_SQRTHALF1;	/* R6 = R6 * square root of 1/2 */ \
	f5 *= XMM_SQRTHALF1;	/* R8 = R8 * square root of 1/2 */ \
\
	f2 = mem1[0];		/* R1 */ \
	f6 = mem2[0];		/* R2 */ \
	f2 -= f6;		/* new R2 = R1 - R2 */ \
	f2 *= 0.5;			/* Mul R2 by HALF */ \
	f6 += f2;		/* new R1 = R1 + R2 */ \
\
	f3 = mem4[0];		/* R4 */ \
	f3 += f2;		/* R4 = R2 + R4 (new R2) */ \
	f2 -= mem4[0];		/* R2 = R2 - R4 (new R4) */ \
\
	f0 = f2;		/* Copy R4 */ \
	f2 -= f5;		/* R4 = R4 - R8 (final R8) */ \
	f5 += f0;		/* R8 = R4 + R8 (final R4) */ \
\
	f0 = f3;		/* Copy R2 */ \
	f3 -= f4;		/* R2 = R2 - R6 (final R6) */ \
	f4 += f0;		/* R6 = R2 + R6 (final R2) */ \
\
	f0 = mem3[0];		/* R3 */ \
	f0 += f6;		/* R3 = R1 + R3 (new R1) */ \
	f6 -= mem3[0];		/* R1 = R1 - R3 (new R3) */ \
\
	dst2[0] = f4;		/* Save R4 */ \
\
	f4 = f6;		/* Copy R3 */ \
	f6 -= f7;		/* R3 = R3 - R7 (final R7) */ \
	f7 += f4;		/* R7 = R3 + R7 (final R3) */ \
\
	f4 = f0;		/* Copy R1 */ \
	f0 -= f1;		/* R1 = R1 - R5 (final R5) */ \
	f1 += f4;		/* R5 = R1 + R5 (final R1) */ \
	xmm0[0] = f0; \
	xmm1[0] = f1; \
	xmm2[0] = f2; \
	xmm3[0] = f3; \
	xmm5[0] = f5; \
	xmm6[0] = f6; \
	xmm7[0] = f7; \
}
#endif

/* */ \
/* ************************************* four-complex-fft variants ****************************************** */ \
/* */ \

#if 0 //IFDEF UNUSED
\
/* Macros to operate on 4 64-byte cache lines doing 4 four_complex_ffts */ \
/* with sin/cos complex multiplies after the butterfly. */ \
\

#define r4_x4cl_four_complex_fft(srcreg,srcinc,d1,d2,screg) \
	r4_x4c_fft_mem(xptr(srcreg),xptr(srcreg+d1),xptr(srcreg+d2),xptr(srcreg+d2+d1),xptr(srcreg+16),xptr(srcreg+d1+16),xptr(srcreg+d2+16),xptr(srcreg+d2+d1+16),screg+0,screg+32,screg+64,srcreg+srcinc,d1,xptr(srcreg),xptr(srcreg+16)); \
/*	xstore	[srcreg], xmm2		;; Save R1 */ \
/*	xstore	[srcreg+16], xmm5	;; Save R2 */ \
	xmm2 = xptr(srcreg+32);	/* R1 */ \
	xmm5 = xptr(srcreg+48);	/* R5 */ \
	xptr(srcreg+32) = xmm7;	/* Save R3 */ \
	xptr(srcreg+48) = xmm6;	/* Save R4 */ \
	xptr(srcreg+d1) = xmm3;	/* Save R5 */ \
	xptr(srcreg+d1+16) = xmm1;	/* Save R6 */ \
	xmm7 = xptr(srcreg+d1+32);	/* R2 */ \
	xmm6 = xptr(srcreg+d1+48);	/* R6 */ \
	xptr(srcreg+d1+32) = xmm0;	/* Save R7 */ \
	xptr(srcreg+d1+48) = xmm4;	/* Save R8 */ \
	xmm3 = xptr(srcreg+d2+32);	/* R3 */ \
	xmm1 = xptr(srcreg+d2+48);	/* R7 */ \
	xmm0 = xptr(srcreg+d2+d1+32);	/* R4 */ \
	xmm4 = xptr(srcreg+d2+d1+48);	/* R8 */ \
	r4_x4c_fft(xmm2, xmm7, xmm3, xmm0, xmm5, xmm6, xmm1, xmm4, screg, 0, srcreg+srcinc+d2, d1, xptr(srcreg+d2), xptr(srcreg+d2+16)); \
/*	xstore	[srcreg+d2], xmm3	;; Save R1 */ \
/*	xstore	[srcreg+d2+16], xmm6	;; Save R2 */ \
	xptr(srcreg+d2+32) = xmm4;	/* Save R3 */ \
	xptr(srcreg+d2+48) = xmm1;	/* Save R4 */ \
	xptr(srcreg+d2+d1) = xmm0;	/* Save R5 */ \
	xptr(srcreg+d2+d1+16) = xmm7;	/* Save R6 */ \
	xptr(srcreg+d2+d1+32) = xmm2;	/* Save R7 */ \
	xptr(srcreg+d2+d1+48) = xmm5;	/* Save R8 */ \
	srcreg += srcinc; \
\
\

#define r4_g4cl_four_complex_fft(srcreg,srcinc,d1,d2,dstreg,dstinc,e1,e2) \
	xprefetch(u8ptr(srcreg+srcinc)); \
	r4_x4c_fft_mem(xptr(srcreg),xptr(srcreg+d1),xptr(srcreg+d2),xptr(srcreg+d2+d1),xptr(srcreg+16),xptr(srcreg+d1+16),xptr(srcreg+d2+16),xptr(srcreg+d2+d1+16),rdi+0,rdi+32,rdi+64,dstreg+dstinc,e1,xptr(dstreg),xptr(dstreg+16)); \
	xprefetch(u8ptr(srcreg+srcinc+d1)); \
/*	xstore	[dstreg], xmm2		;; Save R1 */ \
/*	xstore	[dstreg+16], xmm5	;; Save R2 */ \
	xptr(dstreg+32) = xmm7;	/* Save R3 */ \
	xptr(dstreg+48) = xmm6;	/* Save R4 */ \
	xptr(dstreg+e1) = xmm3;	/* Save R5 */ \
	xptr(dstreg+e1+16) = xmm1;	/* Save R6 */ \
	xptr(dstreg+e1+32) = xmm0;	/* Save R7 */ \
	xptr(dstreg+e1+48) = xmm4;	/* Save R8 */ \
	xprefetchw(u8ptr(srcreg+srcinc+d2)); \
	r4_x4c_fft_mem(xptr(srcreg+32),xptr(srcreg+d1+32),xptr(srcreg+d2+32),xptr(srcreg+d2+d1+32),xptr(srcreg+48),xptr(srcreg+d1+48),xptr(srcreg+d2+48),xptr(srcreg+d2+d1+48),rdi+0,rdi+32,rdi+64,dstreg+dstinc+e2,e1,xptr(dstreg+e2),xptr(dstreg+e2+16)); \
	xprefetchw(u8ptr(srcreg+srcinc+d2+d1)); \
	srcreg += srcinc; \
/*	xstore	[dstreg+e2], xmm2	;; Save R1 */ \
/*	xstore	[dstreg+e2+16], xmm5	;; Save R2 */ \
	xptr(dstreg+e2+32) = xmm7;	/* Save R3 */ \
	xptr(dstreg+e2+48) = xmm6;	/* Save R4 */ \
	xptr(dstreg+e2+e1) = xmm3;	/* Save R5 */ \
	xptr(dstreg+e2+e1+16) = xmm1;	/* Save R6 */ \
	xptr(dstreg+e2+e1+32) = xmm0;	/* Save R7 */ \
	xptr(dstreg+e2+e1+48) = xmm4;	/* Save R8 */ \
	dstreg += dstinc; \

#endif
\

#define r4_x4c_fft_mem(R1,R2,R3,R4,R5,R6,R7,R8,screg1,screg2,screg3,pre1,pre2,dst1,dst2) \
	vec2f64 temp2; \
	xmm0 = R1;			/* R1 */ \
	xmm2 = R3;			/* R3 */ \
	xmm0 -= xmm2;		/* R1 = R1 - R3 (new R3) */ \
	xmm2 += R1;			/* R3 = R1 + R3 (new R1) */ \
\
	xmm4 = R5;			/* I1 */ \
	xmm6 = R7;			/* I3 */ \
	xmm4 -= xmm6;		/* I1 = I1 - I3 (new I3) */ \
	xmm6 += R5;			/* I3 = I1 + I3 (new I1) */ \
\
	xmm1 = R2;			/* R2 */ \
	xmm3 = R4;			/* R4 */ \
	xmm1 -= xmm3;		/* R2 = R2 - R4 (new R4) */ \
	xmm3 += R2;			/* R4 = R2 + R4 (new R2) */ \
\
	xmm5 = R6;			/* I2 */ \
	xmm7 = R8;			/* I4 */ \
	xmm5 -= xmm7;		/* I2 = I2 - I4 (new I4) */ \
	xmm7 += R6;			/* I4 = I2 + I4 (new I2) */ \
\
	temp2 = xmm4; \
	xmm4 -= xmm1;		/* I3 = I3 - R4 (final I4) */ \
	xmm1 += temp2;		/* R4 = I3 + R4 (final I3) */ \
\
	temp2 = xmm0; \
	xmm0 -= xmm5;		/* R3 = R3 - I4 (final R3) */ \
	xmm5 += temp2;		/* I4 = R3 + I4 (final R4) */ \
\
	xprefetchw(u8ptr(pre1)); \
\
	temp2 = xmm2; \
	xmm2 -= xmm3;		/* R1 = R1 - R2 (final R2) */ \
	xmm3 += temp2;		/* R2 = R1 + R2 (final R1) */ \
\
	temp2 = xmm6; \
	xmm6 -= xmm7;		/* I1 = I1 - I2 (final I2) */ \
	xmm7 += temp2;		/* I2 = I1 + I2 (final I1) */ \
\
	dst1 = xmm3;		/* Save R1 */ \
	dst2 = xmm7;		/* Save I1 */ \
\
	xmm7 = xptr(screg1+16);	/* cosine/sine */ \
	xmm3 = xmm0 * xmm7 - xmm1;	/* A3 = R3 * cosine/sine - I3 */ \
	xmm1 = xmm1 * xmm7 + xmm0;	/* B3 = I3 * cosine/sine + R3 */ \
\
	xmm0 = xptr(screg2+16);	/* cosine/sine */ \
	xmm7 = xmm2 * xmm0 - xmm6;	/* A2 = R2 * cosine/sine - I2 */ \
	xmm6 = xmm6 * xmm0 + xmm2;	/* B2 = I2 * cosine/sine + R2 */ \
\
	xprefetchw(u8ptr(pre1+pre2)); \
\
	xmm2 = xptr(screg3+16);	/* cosine/sine */ \
	xmm0 = xmm5 * xmm2 - xmm4;	/* A4 = R4 * cosine/sine - I4 */ \
	xmm4 = xmm4 * xmm2 + xmm5;	/* B4 = I4 * cosine/sine + R4 */ \
\
	xmm2 = xptr(screg1); \
	xmm3 *= xmm2;			/* A3 = A3 * sine (final R3) */ \
	xmm1 *= xmm2;			/* B3 = B3 * sine (final I3) */ \
	xmm2 = xptr(screg2); \
	xmm7 *= xmm2;			/* A2 = A2 * sine (final R2) */ \
	xmm6 *= xmm2;			/* B2 = B2 * sine (final I2) */ \
	xmm2 = xptr(screg3); \
	xmm0 *= xmm2;			/* A4 = A4 * sine (final R4) */ \
	xmm4 *= xmm2;			/* B4 = B4 * sine (final I4) */ \
\
\

#if 0 // unused
#define r4_x4c_fft_partial_mem(r1, r2, r3, r4, r5, r6, r7, r8, mem3, mem4, mem7, mem8, screg, off, pre1, pre2, dst1, dst2) \
	r3 = mem3;		/* R3 */ \
	r3 += r1;			/* R3 = R1 + R3 (new R1) */ \
	r1 -= mem3;		/* R1 = R1 - R3 (new R3) */ \
\
	r7 = mem7;		/* I3 */ \
	r7 += r5;			/* I3 = I1 + I3 (new I1) */ \
	r5 -= mem7;		/* I1 = I1 - I3 (new I3) */ \
\
	r4 = mem4;		/* R4 */ \
	r4 += r2;			/* R4 = R2 + R4 (new R2) */ \
	r2 -= mem4;		/* R2 = R2 - R4 (new R4) */ \
\
	r8 = mem8;		/* I4 */ \
	r8 += r6;			/* I4 = I2 + I4 (new I2) */ \
	r6 -= mem8;		/* I2 = I2 - I4 (new I4) */ \
\
	r5 -= r2;			/* I3 = I3 - R4 (final I4) */ \
	r2 += r2; /* multwo */			/* R4 = R4 * 2 */ \
	r2 += r5;			/* R4 = I3 + R4 (final I3) */ \
\
	r1 -= r6;			/* R3 = R3 - I4 (final R3) */ \
	r6 += r6; /* multwo */			/* I4 = I4 * 2 */ \
	r6 += r1;			/* I4 = R3 + I4 (final R4) */ \
\
	xprefetchw(u8ptr(pre1)); \
\
	r3 -= r4;			/* R1 = R1 - R2 (final R2) */ \
	r4 += r4; /* multwo */			/* R2 = R2 * 2 */ \
	r4 += r3;			/* R2 = R1 + R2 (final R1) */ \
\
	r7 -= r8;			/* I1 = I1 - I2 (final I2) */ \
	r8 += r8; /* multwo */			/* I2 = I2 * 2 */ \
	r8 += r7;			/* I2 = I1 + I2 (final I1) */ \
\
	dst1 = r4;		/* Save R1 */ \
	dst2 = r8;		/* Save I1 */ \
\
	r4 = xptr(screg+off+0+16);	/* cosine/sine */ \
	r4 *= r1;			/* A3 = R3 * cosine/sine */ \
	r4 -= r2;			/* A3 = A3 - I3 */ \
	r2 *= xptr(screg+off+0+16);	/* B3 = I3 * cosine/sine */ \
	r2 += r1;			/* B3 = B3 + R3 */ \
\
	r8 = xptr(screg+off+32+16);	/* cosine/sine */ \
	r8 *= r3;			/* A2 = R2 * cosine/sine */ \
	r8 -= r7;			/* A2 = A2 - I2 */ \
	r7 *= xptr(screg+off+32+16);	/* B2 = I2 * cosine/sine */ \
	r7 += r3;			/* B2 = B2 + R2 */ \
\
	xprefetchw(u8ptr(pre1+pre2)); \
\
	r1 = xptr(screg+off+64+16);	/* cosine/sine */ \
	r1 *= r6;			/* A4 = R4 * cosine/sine */ \
	r1 -= r5;			/* A4 = A4 - I4 */ \
	r5 *= xptr(screg+off+64+16);	/* B4 = I4 * cosine/sine */ \
	r5 += r6;			/* B4 = B4 + R4 */ \
\
	r4 *= xptr(screg+off+0);	/* A3 = A3 * sine (final R3) */ \
	r2 *= xptr(screg+off+0);	/* B3 = B3 * sine (final I3) */ \
	r8 *= xptr(screg+off+32);	/* A2 = A2 * sine (final R2) */ \
	r7 *= xptr(screg+off+32);	/* B2 = B2 * sine (final I2) */ \
	r1 *= xptr(screg+off+64);	/* A4 = A4 * sine (final R4) */ \
	r5 *= xptr(screg+off+64);	/* B4 = B4 * sine (final I4) */ \

#endif
