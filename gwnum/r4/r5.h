#pragma once

/* Copyright 2009-2016 - Mersenne Research, Inc.  All rights reserved */ \
/* Author:  George Woltman */ \
/* Email: woltman@alum.mit.edu */ \
/* */ \
\
/* */ \
/* */ \
/* All new macros for version 26 of gwnum.  Do a radix-5 step in an FFT. */ \
/* The forward FFT macros multiply by the sin/cos values at the end of the macro */ \
/* and the inverse FFTs multiply by the sin/cos values at the start of the macro. */ \
/* We also implement the Daniel J. Bernstein (DJB) "exponent-1" idea to save sin/cos */ \
/* memory. */ \
/* */ \
/* */ \
\
/* */ \
/* ************************************* five-complex-djbfft variants ****************************************** */ \
/* */ \
\

#define r5_x5cl_five_complex_djbfft_preload \
\
\

#define r5_x5cl_five_complex_djbfft(srcreg,srcinc,d1,screg) { \
	vec2f64 *sd0 = (vec2f64*)(srcreg); \
	vec2f64 *sd1 = (vec2f64*)(srcreg+(d1)); \
	vec2f64 *sd2 = (vec2f64*)(srcreg+2*(d1)); \
	vec2f64 *sd3 = (vec2f64*)(srcreg+3*(d1)); \
	vec2f64 *sd4 = (vec2f64*)(srcreg+4*(d1)); \
	vec2f64 TMP5, TMP6, TMP7, TMP8; \
	r5_x5c_djbfft_mem(sd0[0],sd0[1],sd1[0],sd1[1],sd2[0],sd2[1],sd3[0],sd3[1],sd4[0],sd4[1],screg,screg+32,sd0[0],sd0[1],sd2[0]); \
	sd1[0] = xmm3;		/* Save R3 */ \
	sd1[1] = xmm4;		/* Save I3 */ \
/*	xstore	[srcreg+d2], xmm7		;; Save R5 */ \
	sd2[1] = xmm6;		/* Save I5 */ \
	xmm3 = sd0[2];		/* R1 */ \
	xmm4 = sd0[3];		/* I1 */ \
	xmm7 = sd1[2];		/* R2 */ \
	xmm6 = sd1[3];		/* I2 */ \
	sd0[2] = xmm1;		/* Save R2 */ \
	sd0[3] = xmm0;		/* Save I2 */ \
	sd1[2] = xmm5;		/* Save R4 */ \
	sd1[3] = xmm2;		/* Save I4 */ \
	TMP5 = xmm3; \
	TMP6 = xmm4; \
	TMP7 = xmm7; \
	TMP8 = xmm6; \
	r5_x5c_djbfft_mem(TMP5,TMP6,TMP7,TMP8,sd2[2],sd2[3],sd3[2],sd3[3],sd4[2],sd4[3],screg,screg+32,sd2[2],sd2[3],sd4[2]); \
	sd3[0] = xmm1;		/* Save R2 */ \
	sd3[1] = xmm0;		/* Save I2 */ \
	sd3[2] = xmm3;		/* Save R3 */ \
	sd3[3] = xmm4;		/* Save I3 */ \
	sd4[0] = xmm5;		/* Save R4 */ \
	sd4[1] = xmm2;		/* Save I4 */ \
/*	xstore	[srcreg+d4+32], xmm7		;; Save R5 */ \
	sd4[3] = xmm6;		/* Save I5 */ \
	srcreg += srcinc; \
}
\
/* Used in first levels of pass 2.  No swizzling. */ \
\

#define r5_f5cl_five_complex_djbfft_preload \
	r5_f5cl_five_complex_djbfft_common_preload; \
\
\

#define r5_f5cl_five_complex_djbfft(srcreg,srcinc,d1,screg,scoff) \
	r5_f5cl_five_complex_djbfft_common(srcreg,rbx,srcinc,d1,screg,scoff,screg+32,scoff); \
\
\
/* Used in pass 2 where the memory layout is different and the two */ \
/* five-complex-djbffts use different sin/cos data. */ \
\

#define r5_nf5cl_five_complex_djbfft_preload \
	r5_f5cl_five_complex_djbfft_common_preload; \
\
\

#define r5_nf5cl_five_complex_djbfft(srcreg,srcinc,d1,screg,scoff) \
	r5_f5cl_five_complex_djbfft_common(srcreg,0,srcinc,d1,screg,scoff,screg+32,scoff); \
\
\
/* Used in pass 2 where the memory layout is different and the two */ \
/* five-complex-djbffts use different sin/cos data.  Uses two sin/cos pointers */ \
/* because the first levels of pass 2 was radix-3 which only has one sin/cos value. */ \
\

#define r5_nf5cl_2sc_five_complex_djbfft_preload \
	r5_f5cl_five_complex_djbfft_common_preload; \
\
\

#define r5_nf5cl_2sc_five_complex_djbfft(srcreg,srcinc,d1,screg1,scoff1,screg2,scoff2) \
	r5_f5cl_five_complex_djbfft_common(srcreg,0,srcinc,d1,screg1,scoff1,screg2,scoff2); \
\
\
/* The common implementation of five_complex_djbfft */ \
\

#define r5_f5cl_five_complex_djbfft_common_preload \
\
\

#define r5_f5cl_five_complex_djbfft_common(srcreg,srcoff,srcinc,d1,screg1,scoff1,screg2,scoff2) { \
	uintptr_t src_off = (uintptr_t)(srcreg + srcoff); \
	vec2f64 *sd0 = (vec2f64*)(srcreg); \
	vec2f64 *sd1 = (vec2f64*)(srcreg+(d1)); \
	vec2f64 *sd2 = (vec2f64*)(srcreg+2*(d1)); \
	vec2f64 *sd3 = (vec2f64*)(srcreg+3*(d1)); \
	vec2f64 *sd4 = (vec2f64*)(srcreg+4*(d1)); \
	vec2f64 *bd0 = (vec2f64*)(src_off); \
	vec2f64 *bd1 = (vec2f64*)(src_off+(d1)); \
	vec2f64 *bd2 = (vec2f64*)(src_off+2*(d1)); \
	vec2f64 *bd3 = (vec2f64*)(src_off+3*(d1)); \
	vec2f64 *bd4 = (vec2f64*)(src_off+4*(d1)); \
	/*uintptr_t d2 = 2*d1; \
	uintptr_t d3 = 3*d1; \
	uintptr_t d4 = 4*d1; \
	uintptr_t n1 = d2+32; \
	uintptr_t n2 = d4+32; \
	uintptr_t n3 = d1+16; \
	uintptr_t n4 = d1+48; \
	uintptr_t n5 = d3+16; \
	uintptr_t n6 = d3+48;*/ \
	vec2f64 TMP5,TMP6,TMP7,TMP8; \
	r5_x5c_djbfft_mem(bd0[0],bd0[2],bd2[0],bd2[2],bd4[0], bd4[2],bd1[1],bd1[3],bd3[1],bd3[3], screg1,screg2, sd0[0],sd0[2],sd4[0]); \
	sd2[0] = xmm3;		/* Save R3 */ \
	sd2[2] = xmm4;		/* Save I3 */ \
/*	xstore	[srcreg+d4], xmm7		;; Save R5 */ \
	sd4[2] = xmm6;		/* Save I5 */ \
	TMP5 = bd1[0];	/* R1 */ \
	TMP6 = bd1[2];	/* I1 */ \
	TMP7 = bd3[0];	/* R2 */ \
	TMP8 = bd3[2];	/* I2 */ \
	sd1[0] = xmm1;		/* Save R2 */ \
	sd1[2] = xmm0;		/* Save I2 */ \
	sd3[0] = xmm5;		/* Save R4 */ \
	sd3[2] = xmm2;		/* Save I4 */ \
	r5_x5c_djbfft_mem(TMP5,TMP6,TMP7,TMP8,bd0[1], bd0[3],bd2[1],bd2[3],bd4[1],bd4[3], screg1+scoff1,screg2+scoff2, sd0[1],sd0[3],sd4[1]); \
	sd1[1] = xmm1;		/* Save R2 */ \
	sd1[3] = xmm0;		/* Save I2 */ \
	sd2[1] = xmm3;		/* Save R3 */ \
	sd2[3] = xmm4;		/* Save I3 */ \
	sd3[1] = xmm5;		/* Save R4 */ \
	sd3[3] = xmm2;		/* Save I4 */ \
/*	xstore	[srcreg+d4+16], xmm7		;; Save R5 */ \
	sd4[3] = xmm6;		/* Save I5 */ \
	srcreg += srcinc; \
}

/* Do a 5-complex FFT.  A 5-complex FFT is: */ \
/* r25a=r2+r5 */ \
/* r34a=r3+r4 */ \
/* i25s=i2-i5 */ \
/* i34s=i3-i4 */ \
/* outr(0) = r1 + r25a + r34a */ \
/* t1=cos2*r25a + cos4*r34a + r1 */ \
/* t2=sin2*i25s + sin4*i34s */ \
/* outr(1)=t1-t2 */ \
/* outr(4)=t1+t2 */ \
/* t3=cos4*r25a + cos2*r34a + r1 */ \
/* t4=sin4*i25s - sin2*i34s */ \
/* outr(2)=t3-t4 */ \
/* outr(3)=t3+t4 */ \
/* r25s=r2-r5 */ \
/* r34s=r3-r4 */ \
/* i25a=i2+i5 */ \
/* i34a=i3+i4 */ \
/* outi(0)=i1+i25a+i34a */ \
/* t5=cos2*i25a + cos4*i34a + i1 */ \
/* t6=sin2*r25s + sin4*r34s */ \
/* outi(1)=t5+t6 */ \
/* outi(4)=t5-t6 */ \
/* t7=cos4*i25a + cos2*i34a + i1 */ \
/* t8=sin4*r25s - sin2*r34s */ \
/* outi(2)=t7+t8 */ \
/* outi(3)=t7-t8 */ \
/* Where cos2 = cos 2*pi/5 = 0.309, sin2 = 0.951, */ \
/* cos4 =-0.809, sin4 = 0.588 */ \
/* Finally, multiply 4 of the 5 results by twiddle factors. */ \
\

#define r5_x5c_djbfft_mem(r1,i1,r2,i2,r3,i3,r4,i4,r5,i5,screg1,screg2,dstr1,dsti1,dstr5) { \
	vec2f64 TMP1, TMP2, TMP3, TMP4; \
	vec2f64 r25s, r34s, i25a, i34a; \
	xmm4 = r2; \
	xmm5 = r5; \
	xmm0 = xmm4 + xmm5;		/* r25a=r2+r5 */ \
	r25s = xmm4 - xmm5;		/* i25s=r2-r5 */ \
	xmm4 = r3; \
	xmm5 = r4; \
	xmm1 = xmm4 + xmm5;		/* r34a=r3+r4 */ \
	r34s = xmm4 - xmm5;		/* i25s=r3-r4 */ \
	xmm4 = i2; \
	xmm5 = i5; \
	xmm2 = xmm4 - xmm5;		/* i25s=i2-i5 */ \
	i25a = xmm4 + xmm5;		/* i25a=i2+i5 */ \
	xmm4 = i3; \
	xmm5 = i4; \
	xmm3 = xmm4 - xmm5;		/* i34s=i3-i4 */ \
	i34a = xmm4 + xmm5;		/* i34a=i3+i4 */ \
\
	xmm4 = xmm0 * CONST2_P309;	/* cos2*r25a */ \
	xmm4 += xmm1 * CONST2_M809;	/* cos2*r25a + cos4*r34a */ \
	xmm5 = xmm0 * CONST2_M809;	/* cos4*r25a */ \
	xmm5 += xmm1 * CONST2_P309;	/* cos4*r25a+cos2*r34a */ \
	xmm6 = xmm2 * CONST2_P951;	/* sin2*i25s */ \
	xmm6 += xmm3 * CONST2_P588;	/* t2=sin2*i25s + sin4*i34s */ \
	xmm2 = xmm2 * CONST2_P588;	/* sin4*i25s */ \
	xmm2 -= xmm3 * CONST2_P951;	/* t4=sin4*i25s-sin2*i34s */ \
	xmm0 += xmm1;			/* r25a + r34a */ \
	xmm7 = r1; \
	xmm4 += xmm7;			/* t1=cos2*r25a + cos4*r34a + r1 */ \
	xmm0 += xmm7;			/* outr(0) = r1 + r25a + r34a */ \
	xmm5 += xmm7;			/* t3=cos4*r25a+cos2*r34a+r1 */ \
	xmm7 = xmm4; \
	xmm4 -= xmm6;			/* outr(1)=t1-t2 */ \
	xmm6 += xmm7;			/* outr(4)=t1+t2 */ \
	xmm7 = xmm5; \
	xmm5 -= xmm2;			/* outr(2)=t3-t4 */ \
	xmm2 += xmm7;			/* outr(3)=t3+t4 */ \
\
	TMP1 = xmm4;			/* Save new r2 */ \
	TMP2 = xmm5;			/* Save new r3 */ \
	TMP3 = xmm2;			/* Save new r4 */ \
	TMP4 = xmm6;			/* Save new r5 */ \
\
	dstr1 = xmm0;			/* WARNING: dstr1 may be alias to r3 */ \
\
	xmm2 = i25a;			/* i25a=i2+i5 */ \
	xmm3 = i34a;			/* i34a=i3+i4 */ \
	xmm4 = r25s;			/* r25s=r2-r5 */ \
	xmm1 = r34s;			/* r34s=r3-r4 */ \
\
	xmm5 = xmm2; \
	xmm5 += xmm3;			/* i25a+i34a */ \
	xmm6 = xmm2 * CONST2_P309;	/* cos2*i25a */ \
	xmm6 += xmm3 * CONST2_M809;	/* cos2*i25a + cos4*i34a */ \
	xmm0 = xmm4 * CONST2_P951;	/* sin2*r25s */ \
	xmm0 += xmm1 * CONST2_P588;	/* t6=sin2*r25s + sin4*r34s */ \
	xmm2 = xmm2 * CONST2_M809;	/* cos4*i25a */ \
	xmm2 += xmm3 * CONST2_P309;	/* cos4*i25a + cos2*i34a */ \
	xmm4 = xmm4 * CONST2_P588;	/* sin4*r25s */ \
	xmm4 -= xmm1 * CONST2_P951;	/* t8=sin4*r25s - sin2*r34s */ \
	xmm3 = i1; \
	xmm6 += xmm3;			/* t5=cos2*i25a + cos4*i34a + i1 */ \
	xmm2 += xmm3;			/* t7=cos4*i25a + cos2*i34a + i1 */ \
\
	xmm7 = xmm6; \
	xmm6 -= xmm0;			/* outi(4)=t5-t6 */ \
	xmm0 += xmm7;			/* outi(1)=t5+t6 */ \
	xmm1 = xmm2; \
	xmm2 -= xmm4;			/* outi(3)=t7-t8 */ \
	xmm4 += xmm1;			/* outi(2)=t7+t8 */ \
	xmm5 += xmm3;			/* outi(0)=i1+i25a+i34a */ \
\
	xmm3 = xptr(screg1+16); \
	xmm7 = TMP4 * xmm3 + xmm6;	/* A5 = R5 * cosine/sine + I5 */ \
	xmm1 = TMP1 * xmm3 - xmm0;	/* A2 = R2 * cosine/sine - I2 */ \
	xmm6 = xmm6 * xmm3 - TMP4;	/* B5 = I5 * cosine/sine - R5 */ \
	xmm0 = xmm0 * xmm3 + TMP1;	/* B2 = I2 * cosine/sine + R2 */ \
	TMP1 = xptr(screg1); \
	xmm7 *= TMP1;			/* A5 = A5 * sine (new R5) */ \
	xmm1 *= TMP1;			/* A2 = A2 * sine (new R2) */ \
	xmm6 *= TMP1;			/* B5 = B5 * sine (new I5) */ \
	xmm0 *= TMP1;			/* B2 = B2 * sine (new I2) */ \
\
	TMP1 = xptr(screg2+16); \
	TMP4 = xptr(screg2); \
	dsti1 = xmm5; \
	dstr5 = xmm7; \
	xmm3 = TMP2 * TMP1 - xmm4;	/* A3 = R3 * cosine/sine - I3 */ \
	xmm5 = TMP3 * TMP1 + xmm2;	/* A4 = R4 * cosine/sine + I4 */ \
	xmm4 = xmm4 * TMP1 + TMP2;	/* B3 = I3 * cosine/sine + R3 */ \
	xmm2 = xmm2 * TMP1 - TMP3;	/* B4 = I4 * cosine/sine - R4 */ \
	xmm3 *= TMP4;			/* A3 = A3 * sine (new R3) */ \
	xmm5 *= TMP4;			/* A4 = A4 * sine (new R4) */ \
	xmm4 *= TMP4;			/* B3 = B3 * sine (new I3) */ \
	xmm2 *= TMP4;			/* B4 = B4 * sine (new I4) */ \
}

#if 0 //DEF X86_64;
\

#define r5_x5cl_five_complex_djbfft_preload \
	xmm13 = CONST2_P309; \
	xmm14 = CONST2_M809; \
	xmm15 = CONST2_P951; \
\
\

#define r5_x5cl_five_complex_djbfft(srcreg,srcinc,d1,screg) \
	xmm0 = xptr(srcreg+d1+32);		/* Load R2 */ \
	xmm1 = xptr(srcreg+4*d1+32);		/* Load R5 */ \
	xmm2 = xmm0;			/* Copy R2 */ \
	xmm0 += xmm1;			/* r25a=r2+r5				; 1-3 */ \
\
	xmm3 = xptr(srcreg+2*d1+32);		/* Load R3 */ \
	xmm4 = xptr(srcreg+3*d1+32);		/* Load R4 */ \
	xmm5 = xmm3;			/* Copy R3 */ \
	xmm3 += xmm4;			/* r34a=r3+r4				; 2-4 */ \
\
	xmm2 -= xmm1;			/* r25s=r2-r5				; 3-5		avail 1,6-12 */ \
\
	xmm5 -= xmm4;			/* r34s=r3-r4				; 4-6		avail 1,4,6-12 */ \
	xmm6 = xmm13; \
	xmm6 *= xmm0;			/* cos2*r25a				; 4-8		avail 1,4,7-12 */ \
\
	xmm7 = xmm0;			/* Copy r25a */ \
	xmm0 += xmm3;			/* r25a + r34a				; 5-7		avail 1,4,8-12 */ \
	xmm8 = xmm14; \
	xmm8 *= xmm3;			/* cos4*r34a				; 5-9		avail 1,4,9-12 */ \
\
	xmm9 = xptr(srcreg+d1+48);		/* Load I2 */ \
	xmm10 = xptr(srcreg+4*d1+48);		/* Load I5 */ \
	xmm11 = xmm9;			/* Copy I2 */ \
	xmm9 += xmm10;			/* i25a=i2+i5				; 6-8		avail 1,4,12 */ \
	xmm7 *= xmm14;			/* cos4*r25a				; 6-10 */ \
\
	xmm11 -= xmm10;			/* i25s=i2-i5				; 7-9		avail 1,4,10,12 */ \
	xmm3 *= xmm13;			/* cos2*r34a				; 7-11 */ \
\
	xmm12 = xptr(srcreg+2*d1+48);		/* Load I3 */ \
	xmm1 = xptr(srcreg+3*d1+48);		/* Load I4 */ \
	xmm4 = xmm12;			/* Copy I3 */ \
	xmm12 += xmm1;			/* i34a=i3+i4				; 8-10		avail 10 */ \
	xmm10 = xmm15; \
	xmm10 *= xmm2;			/* sin2*r25s				; 8-12		avail none */ \
\
	xmm4 -= xmm1;			/* i34s=i3-i4				; 9-11		avail 1 */ \
	xmm1 = CONST2_P588; \
	xmm1 *= xmm5;			/* sin4*r34s				; 9-13		avail none */ \
\
	xmm6 += xmm8;			/* cos2*r25a + cos4*r34a		; 10-12		avail 8 */ \
	xmm2 *= CONST2_P588;			/* sin4*r25s				; 10-14 */ \
\
	xmm8 = xptr(srcreg+32);		/* Load R1 */ \
	xmm0 += xmm8;			/* outr(0) = r1+r25a+r34a (Final R1)	; 11-13		avail none storable 0 */ \
	xmm5 *= xmm15;			/* sin2*r34s				; 11-15 */ \
\
	xmm7 += xmm3;			/* cos4*r25a + cos2*r34a		; 12-14		avail 3 storable 0 */ \
	xmm3 = xmm13; \
	xmm3 *= xmm9;			/* cos2*i25a				; 12-16		avail none storable 0 */ \
\
	xmm6 += xmm8;			/* t1=cos2*r25a + cos4*r34a + r1	; 13-15 */ \
	xptr(srcreg+2*d1+32) = xmm0;		/* Save R1 */ \
	xmm0 = xmm14; \
	xmm0 *= xmm12;			/* cos4*i34a				; 13-17		avail none */ \
\
	xmm10 += xmm1;			/* t6=sin2*r25s + sin4*r34s		; 14-16		avail 1 */ \
	xmm1 = xmm9;			/* Copy i25a */ \
	xmm9 *= xmm14;			/* cos4*i25a				; 14-18		avail none */ \
\
	xmm1 += xmm12;			/* i25a+i34a				; 15-17(11) */ \
	xmm12 *= xmm13;			/* cos2*i34a				; 15-19 */ \
\
	xmm7 += xmm8;			/* t3=cos4*r25a+cos2*r34a+r1		; 16-18		avail 8 */ \
	xmm8 = xmm15; \
	xmm8 *= xmm11;			/* sin2*i25s				; 16-20		avail none */ \
\
	xmm2 -= xmm5;			/* t8=sin4*r25s - sin2*r34s		; 17-19(16)	avail 5 */ \
	xmm5 = CONST2_P588; \
	xmm5 *= xmm4;			/* sin4*i34s				; 17-21		avail none */ \
\
	xmm3 += xmm0;			/* cos2*i25a + cos4*i34a		; 18-20		avail 0 */ \
	xmm11 *= CONST2_P588;			/* sin4*i25s				; 18-22 */ \
\
	xmm0 = xptr(srcreg+48);		/* Load I1 */ \
	xmm1 += xmm0;			/* outi(0)=i1+i25a+i34a (Final I1)	; 19-21(18)	avail none storable 1 */ \
	xmm4 *= xmm15;			/* sin2*i34s				; 19-23 */ \
\
	xmm9 += xmm12;			/* cos4*i25a + cos2*i34a		; 20-22		avail 12 storable 1 */ \
\
	xmm3 += xmm0;			/* t5=cos2*i25a + cos4*i34a + i1	; 21-23 */ \
\
	xmm8 += xmm5;			/* t2=sin2*i25s + sin4*i34s		; 22-24		avail 12,5 storable 1 */ \
	xptr(srcreg+2*d1+48) = xmm1;		/* Save I1 */ \
\
	xmm9 += xmm0;			/* t7=cos4*i25a + cos2*i34a + i1	; 23-25		avail 12,5,1,0 */ \
\
	xmm11 -= xmm4;			/* t4=sin4*i25s-sin2*i34s		; 24-26		avail 12,5,1,0,4 */ \
\
	xmm12 = xmm10;			/* Copy t6 */ \
	xmm10 += xmm3;			/* outi(1)=t5+t6 (new I2)		; 25-27 */ \
\
	xmm3 -= xmm12;			/* outi(4)=t5-t6 (new I5)		; 26-28 */ \
\
	xmm5 = xmm6;			/* Copy t1 */ \
	xmm6 -= xmm8;			/* outr(1)=t1-t2 (new R2)		; 27-29 */ \
\
	xmm8 += xmm5;			/* outr(4)=t1+t2 (new R5)		; 28-30 */ \
	xmm1 = xptr(screg+16);		/* cosine/sine */ \
	xmm0 = xmm10;			/* Copy I2 */ \
	xmm10 *= xmm1;			/* B2 = I2 * cosine/sine		; 28-32		avail 4,12,5 */ \
\
	xmm4 = xmm2;			/* Copy t8 */ \
	xmm2 += xmm9;			/* outi(2)=t7+t8 (new I3)		; 29-31 */ \
	xmm12 = xmm3;			/* Copy I5 */ \
	xmm3 *= xmm1;			/* B5 = I5 * cosine/sine		; 29-33		avail 5 */ \
	xmm9 -= xmm4;			/* outi(3)=t7-t8 (new I4)		; 30-32		avail 5,4 */ \
	xmm5 = xmm6;			/* Copy R2 */ \
	xmm6 *= xmm1;			/* A2 = R2 * cosine/sine		; 30-34		avail 4 */ \
	xmm4 = xmm7;			/* Copy t3 */ \
	xmm7 -= xmm11;			/* outr(2)=t3-t4 (new R3)		; 31-33 */ \
	xmm1 *= xmm8;			/* A5 = R5 * cosine/sine		; 31-35 */ \
	xmm11 += xmm4;			/* outr(3)=t3+t4 (new R4)		; 32-34 */ \
\
	xmm4 = xmm2;			/* Copy I3 */ \
	xmm10 += xmm5;			/* B2 = B2 + R2				; 33-35 */ \
	xmm5 = xptr(screg+32+16);		/* cosine/sine */ \
	xmm2 *= xmm5;			/* B3 = I3 * cosine/sine		; 32-36		avail none */ \
\
	xmm3 -= xmm8;			/* B5 = B5 - R5				; 34-36 */ \
	xmm8 = xmm9;			/* Copy I4 */ \
	xmm9 *= xmm5;			/* B4 = I4 * cosine/sine		; 33-37 */ \
\
	xmm6 -= xmm0;			/* A2 = A2 - I2				; 35-37 */ \
	xmm0 = xmm7;			/* Copy R3 */ \
	xmm7 *= xmm5;	 		/* A3 = R3 * cosine/sine		; 34-38 */ \
\
	xmm5 *= xmm11;			/* A4 = R4 * cosine/sine		; 35-39 */ \
\
	xmm1 += xmm12;			/* A5 = A5 + I5				; 36-38		avail 12 */ \
	xmm12 = xptr(screg);			/* sine */ \
	xmm10 *= xmm12;			/* B2 = B2 * sine (final I2)		; 36-40		avail none */ \
\
	xmm2 += xmm0;			/* B3 = B3 + R3				; 37-39 */ \
	xmm3 *= xmm12;			/* B5 = B5 * sine (final I5)		; 37-41 */ \
	xmm0 = xptr(screg+32); \
\
	xmm9 -= xmm11;			/* B4 = B4 - R4				; 38-40 */ \
	xmm6 *= xmm12;			/* A2 = A2 * sine (final R2)		; 38-42 */ \
	xmm11 = xptr(srcreg+d1);		/* Load R2 */ \
\
	xmm7 -= xmm4;			/* A3 = A3 - I3				; 39-41 */ \
	xmm1 *= xmm12;			/* A5 = A5 * sine (final R5)		; 39-43		avail 12 */ \
	xmm4 = xptr(srcreg+4*d1);		/* Load R5 */ \
\
	xmm5 += xmm8;			/* A4 = A4 + I4				; 40-42		avail 5,12 */ \
	xmm2 *= xmm0;			/* B3 = B3 * sine (final I3)		; 40-44 */ \
	xmm12 = xmm11;			/* Copy R2						avail 5 */ \
\
	xmm11 += xmm4;			/* r25a=r2+r5				; 1-3 */ \
	xmm9 *= xmm0;			/* B4 = B4 * sine (final I4)		; 41-45 */ \
	xmm8 = xptr(srcreg+2*d1);		/* Load R3						avail none */ \
\
	xptr(srcreg+4*d1+48) = xmm3;		/* Save I5				; 42 */ \
	xmm3 = xmm8;			/* Copy R3 */ \
	xptr(srcreg+4*d1+32) = xmm1;		/* Save R5				; 44 */ \
	xmm1 = xptr(srcreg+3*d1);		/* Load R4 */ \
	xmm8 += xmm1;			/* r34a=r3+r4				; 2-4 */ \
	xmm7 *= xmm0;			/* A3 = A3 * sine (final R3)		; 42-46 */ \
\
	xmm12 -= xmm4;			/* r25s=r2-r5				; 3-5 */ \
	xmm5 *= xmm0;			/* A4 = A4 * sine (final R4)		; 43-47		avail 8 */ \
\
	xmm3 -= xmm1;			/* r34s=r3-r4				; 4-6 */ \
	xmm0 = xmm13; \
	xmm0 *= xmm11;			/* cos2*r25a				; 4-8 */ \
\
	xptr(srcreg+3*d1) = xmm6;		/* Save R2				; 43 */ \
	xmm6 = xmm11;			/* Copy r25a */ \
	xmm11 += xmm8;			/* r25a + r34a				; 5-7 */ \
	xptr(srcreg+3*d1+48) = xmm2;		/* Save I3				; 45 */ \
	xmm2 = xmm14; \
	xmm2 *= xmm8;			/* cos4*r34a				; 5-9 */ \
\
	xptr(srcreg+3*d1+32) = xmm7;		/* Save R3				; 47 */ \
	xmm7 = xptr(srcreg+d1+16);		/* Load I2 */ \
	xptr(srcreg+4*d1) = xmm5;		/* Save R4				; 48 */ \
	xmm5 = xptr(srcreg+4*d1+16);		/* Load I5 */ \
	xptr(srcreg+4*d1+16) = xmm9;		/* Save I4				; 46 */ \
	xmm9 = xmm7;			/* Copy I2 */ \
	xmm7 += xmm5;			/* i25a=i2+i5				; 6-8 */ \
	xmm6 *= xmm14;			/* cos4*r25a				; 6-10 */ \
\
	xmm9 -= xmm5;			/* i25s=i2-i5				; 7-9 */ \
	xmm8 *= xmm13;			/* cos2*r34a				; 7-11 */ \
\
	xmm4 = xptr(srcreg+3*d1+16);		/* Load I4 */ \
	xptr(srcreg+3*d1+16) = xmm10;		/* Save I2				; 41 */ \
	xmm10 = xptr(srcreg+2*d1+16);		/* Load I3 */ \
	xmm1 = xmm10;			/* Copy I3 */ \
	xmm10 += xmm4;			/* i34a=i3+i4				; 8-10 */ \
	xmm5 = xmm15; \
	xmm5 *= xmm12;			/* sin2*r25s				; 8-12 */ \
\
	xmm1 -= xmm4;			/* i34s=i3-i4				; 9-11 */ \
	xmm4 = CONST2_P588; \
	xmm4 *= xmm3;			/* sin4*r34s				; 9-13 */ \
\
	xmm0 += xmm2;			/* cos2*r25a + cos4*r34a		; 10-12 */ \
	xmm12 *= CONST2_P588;			/* sin4*r25s				; 10-14 */ \
\
	xmm2 = xptr(srcreg);			/* Load R1 */ \
	xmm11 += xmm2;			/* outr(0) = r1+r25a+r34a (Final R1)	; 11-13 */ \
	xmm3 *= xmm15;			/* sin2*r34s				; 11-15 */ \
\
	xmm6 += xmm8;			/* cos4*r25a + cos2*r34a		; 12-14 */ \
	xmm8 = xmm13; \
	xmm8 *= xmm7;			/* cos2*i25a				; 12-16 */ \
\
	xmm0 += xmm2;			/* t1=cos2*r25a + cos4*r34a + r1	; 13-15 */ \
	xptr(srcreg) = xmm11;			/* Save R1 */ \
	xmm11 = xmm14; \
	xmm11 *= xmm10;			/* cos4*i34a				; 13-17 */ \
\
	xmm5 += xmm4;			/* t6=sin2*r25s + sin4*r34s		; 14-16 */ \
	xmm4 = xmm7;			/* Copy i25a */ \
	xmm7 *= xmm14;			/* cos4*i25a				; 14-18 */ \
\
	xmm4 += xmm10;			/* i25a+i34a				; 15-17(11) */ \
	xmm10 *= xmm13;			/* cos2*i34a				; 15-19 */ \
\
	xmm6 += xmm2;			/* t3=cos4*r25a+cos2*r34a+r1		; 16-18 */ \
	xmm2 = xmm15; \
	xmm2 *= xmm9;			/* sin2*i25s				; 16-20 */ \
\
	xmm12 -= xmm3;			/* t8=sin4*r25s - sin2*r34s		; 17-19(16) */ \
	xmm3 = CONST2_P588; \
	xmm3 *= xmm1;			/* sin4*i34s				; 17-21 */ \
\
	xmm8 += xmm11;			/* cos2*i25a + cos4*i34a		; 18-20 */ \
	xmm9 *= CONST2_P588;			/* sin4*i25s				; 18-22 */ \
\
	xmm11 = xptr(srcreg+16);		/* Load I1 */ \
	xmm4 += xmm11;			/* outi(0)=i1+i25a+i34a (Final I1)	; 19-21(18) */ \
	xmm1 *= xmm15;			/* sin2*i34s				; 19-23 */ \
\
	xmm7 += xmm10;			/* cos4*i25a + cos2*i34a		; 20-22 */ \
\
	xmm8 += xmm11;			/* t5=cos2*i25a + cos4*i34a + i1	; 21-23 */ \
\
	xmm2 += xmm3;			/* t2=sin2*i25s + sin4*i34s		; 22-24 */ \
	xptr(srcreg+16) = xmm4;		/* Save I1 */ \
\
	xmm7 += xmm11;			/* t7=cos4*i25a + cos2*i34a + i1	; 23-25 */ \
\
	xmm9 -= xmm1;			/* t4=sin4*i25s-sin2*i34s		; 24-26 */ \
\
	xmm10 = xmm5;			/* Copy t6 */ \
	xmm5 += xmm8;			/* outi(1)=t5+t6 (new I2)		; 25-27 */ \
\
	xmm8 -= xmm10;			/* outi(4)=t5-t6 (new I5)		; 26-28 */ \
\
	xmm3 = xmm0;			/* Copy t1 */ \
	xmm0 -= xmm2;			/* outr(1)=t1-t2 (new R2)		; 27-29 */ \
\
	xmm2 += xmm3;			/* outr(4)=t1+t2 (new R5)		; 28-30 */ \
	xmm4 = xptr(screg+16);		/* cosine/sine */ \
	xmm11 = xmm5;			/* Copy I2 */ \
	xmm5 *= xmm4;			/* B2 = I2 * cosine/sine		; 28-32 */ \
\
	xmm1 = xmm12;			/* Copy t8 */ \
	xmm12 += xmm7;			/* outi(2)=t7+t8 (new I3)		; 29-31 */ \
	xmm10 = xmm8;			/* Copy I5 */ \
	xmm8 *= xmm4;			/* B5 = I5 * cosine/sine		; 29-33 */ \
	xmm7 -= xmm1;			/* outi(3)=t7-t8 (new I4)		; 30-32 */ \
	xmm3 = xmm0;			/* Copy R2 */ \
	xmm0 *= xmm4;			/* A2 = R2 * cosine/sine		; 30-34 */ \
	xmm1 = xmm6;			/* Copy t3 */ \
	xmm6 -= xmm9;			/* outr(2)=t3-t4 (new R3)		; 31-33 */ \
	xmm4 *= xmm2;			/* A5 = R5 * cosine/sine		; 31-35 */ \
	xmm9 += xmm1;			/* outr(3)=t3+t4 (new R4)		; 32-34 */ \
\
	xmm1 = xmm12;			/* Copy I3 */ \
	xmm5 += xmm3;			/* B2 = B2 + R2				; 33-35 */ \
	xmm3 = xptr(screg+32+16);		/* cosine/sine */ \
	xmm12 *= xmm3;			/* B3 = I3 * cosine/sine		; 32-36 */ \
\
	xmm8 -= xmm2;			/* B5 = B5 - R5				; 34-36 */ \
	xmm2 = xmm7;			/* Copy I4 */ \
	xmm7 *= xmm3;			/* B4 = I4 * cosine/sine		; 33-37 */ \
\
	xmm0 -= xmm11;			/* A2 = A2 - I2				; 35-37 */ \
	xmm11 = xmm6;			/* Copy R3 */ \
	xmm6 *= xmm3;	 		/* A3 = R3 * cosine/sine		; 34-38 */ \
\
	xmm3 *= xmm9;			/* A4 = R4 * cosine/sine		; 35-39 */ \
\
	xmm4 += xmm10;			/* A5 = A5 + I5				; 36-38 */ \
	xmm10 = xptr(screg);			/* sine */ \
	xmm5 *= xmm10;			/* B2 = B2 * sine (final I2)		; 36-40 */ \
\
	xmm12 += xmm11;			/* B3 = B3 + R3				; 37-39 */ \
	xmm8 *= xmm10;			/* B5 = B5 * sine (final I5)		; 37-41 */ \
	xmm7 -= xmm9;			/* B4 = B4 - R4				; 38-40 */ \
	xmm0 *= xmm10;			/* A2 = A2 * sine (final R2)		; 38-42 */ \
	xmm6 -= xmm1;			/* A3 = A3 - I3				; 39-41 */ \
	xmm4 *= xmm10;			/* A5 = A5 * sine (final R5)		; 39-43 */ \
	xmm3 += xmm2;			/* A4 = A4 + I4				; 40-42 */ \
\
	xmm11 = xptr(screg+32); \
	xmm12 *= xmm11;			/* B3 = B3 * sine (final I3)		; 40-44 */ \
	xmm7 *= xmm11;			/* B4 = B4 * sine (final I4)		; 41-45 */ \
	xptr(srcreg+48) = xmm5;		/* Save I2				; 41 */ \
	xmm6 *= xmm11;			/* A3 = A3 * sine (final R3)		; 42-46 */ \
	xptr(srcreg+2*d1+16) = xmm8;		/* Save I5				; 42 */ \
	xmm3 *= xmm11;			/* A4 = A4 * sine (final R4)		; 43-47 */ \
	xptr(srcreg+32) = xmm0;		/* Save R2				; 43 */ \
\
	xptr(srcreg+2*d1) = xmm4;		/* Save R5 */ \
	xptr(srcreg+d1+16) = xmm12;		/* Save I3 */ \
	xptr(srcreg+d1+48) = xmm7;		/* Save I4 */ \
	xptr(srcreg+d1) = xmm6;		/* Save R3 */ \
	xptr(srcreg+d1+32) = xmm3;		/* Save R4 */ \
\
	srcreg += srcinc; \
\
\

#define r5_f5cl_five_complex_djbfft_common_preload \
	xmm13 = CONST2_P309; \
	xmm14 = CONST2_M809; \
	xmm15 = CONST2_P951; \
\
\

#define r5_f5cl_five_complex_djbfft_common(srcreg,srcoff,srcinc,d1,screg1,scoff1,screg2,scoff2) \
	xmm0 = xptr(srcreg+srcoff+3*d1);	/* Load R2 */ \
	xmm1 = xptr(srcreg+srcoff+4*d1+16);	/* Load R5 */ \
	xmm2 = xmm0;			/* Copy R2 */ \
	xmm0 += xmm1;			/* r25a=r2+r5				; 1-3 */ \
\
	xmm3 = xptr(srcreg+srcoff+16);	/* Load R3 */ \
	xmm4 = xptr(srcreg+srcoff+2*d1+16);	/* Load R4 */ \
	xmm5 = xmm3;			/* Copy R3 */ \
	xmm3 += xmm4;			/* r34a=r3+r4				; 2-4 */ \
\
	xmm2 -= xmm1;			/* r25s=r2-r5				; 3-5		avail 1,6-12 */ \
\
	xmm5 -= xmm4;			/* r34s=r3-r4				; 4-6		avail 1,4,6-12 */ \
	xmm6 = xmm13; \
	xmm6 *= xmm0;			/* cos2*r25a				; 4-8		avail 1,4,7-12 */ \
\
	xmm7 = xmm0;			/* Copy r25a */ \
	xmm0 += xmm3;			/* r25a + r34a				; 5-7		avail 1,4,8-12 */ \
	xmm8 = xmm14; \
	xmm8 *= xmm3;			/* cos4*r34a				; 5-9		avail 1,4,9-12 */ \
\
	xmm9 = xptr(srcreg+srcoff+3*d1+32);	/* Load I2 */ \
	xmm10 = xptr(srcreg+srcoff+4*d1+48);	/* Load I5 */ \
	xmm11 = xmm9;			/* Copy I2 */ \
	xmm9 += xmm10;			/* i25a=i2+i5				; 6-8		avail 1,4,12 */ \
	xmm7 *= xmm14;			/* cos4*r25a				; 6-10 */ \
\
	xmm11 -= xmm10;			/* i25s=i2-i5				; 7-9		avail 1,4,10,12 */ \
	xmm3 *= xmm13;			/* cos2*r34a				; 7-11 */ \
\
	xmm12 = xptr(srcreg+srcoff+48);	/* Load I3 */ \
	xmm1 = xptr(srcreg+srcoff+2*d1+48);	/* Load I4 */ \
	xmm4 = xmm12;			/* Copy I3 */ \
	xmm12 += xmm1;			/* i34a=i3+i4				; 8-10		avail 10 */ \
	xmm10 = xmm15; \
	xmm10 *= xmm2;			/* sin2*r25s				; 8-12		avail none */ \
\
	xmm4 -= xmm1;			/* i34s=i3-i4				; 9-11		avail 1 */ \
	xmm1 = CONST2_P588; \
	xmm1 *= xmm5;			/* sin4*r34s				; 9-13		avail none */ \
\
	xmm6 += xmm8;			/* cos2*r25a + cos4*r34a		; 10-12		avail 8 */ \
	xmm2 *= CONST2_P588;			/* sin4*r25s				; 10-14 */ \
\
	xmm8 = xptr(srcreg+srcoff+d1);	/* Load R1 */ \
	xmm0 += xmm8;			/* outr(0) = r1+r25a+r34a (Final R1)	; 11-13		avail none storable 0 */ \
	xmm5 *= xmm15;			/* sin2*r34s				; 11-15 */ \
\
	xmm7 += xmm3;			/* cos4*r25a + cos2*r34a		; 12-14		avail 3 storable 0 */ \
	xmm3 = xmm13; \
	xmm3 *= xmm9;			/* cos2*i25a				; 12-16		avail none storable 0 */ \
\
	xmm6 += xmm8;			/* t1=cos2*r25a + cos4*r34a + r1	; 13-15 */ \
	xptr(srcreg+16) = xmm0;		/* Save R1 */ \
	xmm0 = xmm14; \
	xmm0 *= xmm12;			/* cos4*i34a				; 13-17		avail none */ \
\
	xmm10 += xmm1;			/* t6=sin2*r25s + sin4*r34s		; 14-16		avail 1 */ \
	xmm1 = xmm9;			/* Copy i25a */ \
	xmm9 *= xmm14;			/* cos4*i25a				; 14-18		avail none */ \
\
	xmm1 += xmm12;			/* i25a+i34a				; 15-17(11) */ \
	xmm12 *= xmm13;			/* cos2*i34a				; 15-19 */ \
\
	xmm7 += xmm8;			/* t3=cos4*r25a+cos2*r34a+r1		; 16-18		avail 8 */ \
	xmm8 = xmm15; \
	xmm8 *= xmm11;			/* sin2*i25s				; 16-20		avail none */ \
\
	xmm2 -= xmm5;			/* t8=sin4*r25s - sin2*r34s		; 17-19(16)	avail 5 */ \
	xmm5 = CONST2_P588; \
	xmm5 *= xmm4;			/* sin4*i34s				; 17-21		avail none */ \
\
	xmm3 += xmm0;			/* cos2*i25a + cos4*i34a		; 18-20		avail 0 */ \
	xmm11 *= CONST2_P588;			/* sin4*i25s				; 18-22 */ \
\
	xmm0 = xptr(srcreg+srcoff+d1+32);	/* Load I1 */ \
	xmm1 += xmm0;			/* outi(0)=i1+i25a+i34a (Final I1)	; 19-21(18)	avail none storable 1 */ \
	xmm4 *= xmm15;			/* sin2*i34s				; 19-23 */ \
\
	xmm9 += xmm12;			/* cos4*i25a + cos2*i34a		; 20-22		avail 12 storable 1 */ \
\
	xmm3 += xmm0;			/* t5=cos2*i25a + cos4*i34a + i1	; 21-23 */ \
\
	xmm8 += xmm5;			/* t2=sin2*i25s + sin4*i34s		; 22-24		avail 12,5 storable 1 */ \
	xptr(srcreg+48) = xmm1;		/* Save I1 */ \
\
	xmm9 += xmm0;			/* t7=cos4*i25a + cos2*i34a + i1	; 23-25		avail 12,5,1,0 */ \
\
	xmm11 -= xmm4;			/* t4=sin4*i25s-sin2*i34s		; 24-26		avail 12,5,1,0,4 */ \
\
	xmm12 = xmm10;			/* Copy t6 */ \
	xmm10 += xmm3;			/* outi(1)=t5+t6 (new I2)		; 25-27 */ \
\
	xmm3 -= xmm12;			/* outi(4)=t5-t6 (new I5)		; 26-28 */ \
\
	xmm5 = xmm6;			/* Copy t1 */ \
	xmm6 -= xmm8;			/* outr(1)=t1-t2 (new R2)		; 27-29 */ \
\
	xmm8 += xmm5;			/* outr(4)=t1+t2 (new R5)		; 28-30 */ \
	xmm1 = xptr(screg1+scoff1+16);	/* cosine/sine */ \
	xmm0 = xmm10;			/* Copy I2 */ \
	xmm10 *= xmm1;			/* B2 = I2 * cosine/sine		; 28-32		avail 4,12,5 */ \
\
	xmm4 = xmm2;			/* Copy t8 */ \
	xmm2 += xmm9;			/* outi(2)=t7+t8 (new I3)		; 29-31 */ \
	xmm12 = xmm3;			/* Copy I5 */ \
	xmm3 *= xmm1;			/* B5 = I5 * cosine/sine		; 29-33		avail 5 */ \
	xmm9 -= xmm4;			/* outi(3)=t7-t8 (new I4)		; 30-32		avail 5,4 */ \
	xmm5 = xmm6;			/* Copy R2 */ \
	xmm6 *= xmm1;			/* A2 = R2 * cosine/sine		; 30-34		avail 4 */ \
	xmm4 = xmm7;			/* Copy t3 */ \
	xmm7 -= xmm11;			/* outr(2)=t3-t4 (new R3)		; 31-33 */ \
	xmm1 *= xmm8;			/* A5 = R5 * cosine/sine		; 31-35 */ \
	xmm11 += xmm4;			/* outr(3)=t3+t4 (new R4)		; 32-34 */ \
\
	xmm4 = xmm2;			/* Copy I3 */ \
	xmm10 += xmm5;			/* B2 = B2 + R2				; 33-35 */ \
	xmm5 = xptr(screg2+scoff2+16);	/* cosine/sine */ \
	xmm2 *= xmm5;			/* B3 = I3 * cosine/sine		; 32-36		avail none */ \
\
	xmm3 -= xmm8;			/* B5 = B5 - R5				; 34-36 */ \
	xmm8 = xmm9;			/* Copy I4 */ \
	xmm9 *= xmm5;			/* B4 = I4 * cosine/sine		; 33-37 */ \
\
	xmm6 -= xmm0;			/* A2 = A2 - I2				; 35-37 */ \
	xmm0 = xmm7;			/* Copy R3 */ \
	xmm7 *= xmm5;	 		/* A3 = R3 * cosine/sine		; 34-38 */ \
\
	xmm5 *= xmm11;			/* A4 = R4 * cosine/sine		; 35-39 */ \
\
	xmm1 += xmm12;			/* A5 = A5 + I5				; 36-38		avail 12 */ \
	xmm12 = xptr(screg1+scoff1);		/* sine */ \
	xmm10 *= xmm12;			/* B2 = B2 * sine (final I2)		; 36-40		avail none */ \
\
	xmm2 += xmm0;			/* B3 = B3 + R3				; 37-39 */ \
	xmm3 *= xmm12;			/* B5 = B5 * sine (final I5)		; 37-41 */ \
	xmm0 = xptr(screg2+scoff2); \
\
	xmm9 -= xmm11;			/* B4 = B4 - R4				; 38-40 */ \
	xmm6 *= xmm12;			/* A2 = A2 * sine (final R2)		; 38-42 */ \
	xmm11 = xptr(srcreg+srcoff+2*d1);	/* Load R2 */ \
\
	xmm7 -= xmm4;			/* A3 = A3 - I3				; 39-41 */ \
	xmm1 *= xmm12;			/* A5 = A5 * sine (final R5)		; 39-43		avail 12 */ \
	xmm4 = xptr(srcreg+srcoff+3*d1+16);	/* Load R5 */ \
\
	xmm5 += xmm8;			/* A4 = A4 + I4				; 40-42		avail 5,12 */ \
	xmm2 *= xmm0;			/* B3 = B3 * sine (final I3)		; 40-44 */ \
	xmm12 = xmm11;			/* Copy R2						avail 5 */ \
\
	xmm11 += xmm4;			/* r25a=r2+r5				; 1-3 */ \
	xmm9 *= xmm0;			/* B4 = B4 * sine (final I4)		; 41-45 */ \
	xmm8 = xptr(srcreg+srcoff+4*d1);	/* Load R3						avail none */ \
\
	xptr(srcreg+4*d1+48) = xmm3;		/* Save I5				; 42 */ \
	xmm3 = xmm8;			/* Copy R3 */ \
	xptr(srcreg+4*d1+16) = xmm1;		/* Save R5				; 44 */ \
	xmm1 = xptr(srcreg+srcoff+d1+16);	/* Load R4 */ \
	xmm8 += xmm1;			/* r34a=r3+r4				; 2-4 */ \
	xmm7 *= xmm0;			/* A3 = A3 * sine (final R3)		; 42-46 */ \
\
	xmm12 -= xmm4;			/* r25s=r2-r5				; 3-5 */ \
	xmm5 *= xmm0;			/* A4 = A4 * sine (final R4)		; 43-47		avail 8 */ \
\
	xmm3 -= xmm1;			/* r34s=r3-r4				; 4-6 */ \
	xmm0 = xmm13; \
	xmm0 *= xmm11;			/* cos2*r25a				; 4-8 */ \
\
	xptr(srcreg+d1+16) = xmm6;		/* Save R2				; 43 */ \
	xmm6 = xmm11;			/* Copy r25a */ \
	xmm11 += xmm8;			/* r25a + r34a				; 5-7 */ \
	xptr(srcreg+2*d1+48) = xmm2;		/* Save I3				; 45 */ \
	xmm2 = xmm14; \
	xmm2 *= xmm8;			/* cos4*r34a				; 5-9 */ \
\
	xptr(srcreg+2*d1+16) = xmm7;		/* Save R3				; 47 */ \
	xmm7 = xptr(srcreg+srcoff+2*d1+32);	/* Load I2 */ \
	xptr(srcreg+3*d1+16) = xmm5;		/* Save R4				; 48 */ \
	xmm5 = xptr(srcreg+srcoff+3*d1+48);	/* Load I5 */ \
	xptr(srcreg+3*d1+48) = xmm9;		/* Save I4				; 46 */ \
	xmm9 = xmm7;			/* Copy I2 */ \
	xmm7 += xmm5;			/* i25a=i2+i5				; 6-8 */ \
	xmm6 *= xmm14;			/* cos4*r25a				; 6-10 */ \
\
	xmm9 -= xmm5;			/* i25s=i2-i5				; 7-9 */ \
	xmm8 *= xmm13;			/* cos2*r34a				; 7-11 */ \
\
	xmm4 = xptr(srcreg+srcoff+d1+48);	/* Load I4 */ \
	xptr(srcreg+d1+48) = xmm10;		/* Save I2				; 41 */ \
	xmm10 = xptr(srcreg+srcoff+4*d1+32);	/* Load I3 */ \
	xmm1 = xmm10;			/* Copy I3 */ \
	xmm10 += xmm4;			/* i34a=i3+i4				; 8-10 */ \
	xmm5 = xmm15; \
	xmm5 *= xmm12;			/* sin2*r25s				; 8-12 */ \
\
	xmm1 -= xmm4;			/* i34s=i3-i4				; 9-11 */ \
	xmm4 = CONST2_P588; \
	xmm4 *= xmm3;			/* sin4*r34s				; 9-13 */ \
\
	xmm0 += xmm2;			/* cos2*r25a + cos4*r34a		; 10-12 */ \
	xmm12 *= CONST2_P588;			/* sin4*r25s				; 10-14 */ \
\
	xmm2 = xptr(srcreg+srcoff);		/* Load R1 */ \
	xmm11 += xmm2;			/* outr(0) = r1+r25a+r34a (Final R1)	; 11-13 */ \
	xmm3 *= xmm15;			/* sin2*r34s				; 11-15 */ \
\
	xmm6 += xmm8;			/* cos4*r25a + cos2*r34a		; 12-14 */ \
	xmm8 = xmm13; \
	xmm8 *= xmm7;			/* cos2*i25a				; 12-16 */ \
\
	xmm0 += xmm2;			/* t1=cos2*r25a + cos4*r34a + r1	; 13-15 */ \
	xptr(srcreg) = xmm11;			/* Save R1 */ \
	xmm11 = xmm14; \
	xmm11 *= xmm10;			/* cos4*i34a				; 13-17 */ \
\
	xmm5 += xmm4;			/* t6=sin2*r25s + sin4*r34s		; 14-16 */ \
	xmm4 = xmm7;			/* Copy i25a */ \
	xmm7 *= xmm14;			/* cos4*i25a				; 14-18 */ \
\
	xmm4 += xmm10;			/* i25a+i34a				; 15-17(11) */ \
	xmm10 *= xmm13;			/* cos2*i34a				; 15-19 */ \
\
	xmm6 += xmm2;			/* t3=cos4*r25a+cos2*r34a+r1		; 16-18 */ \
	xmm2 = xmm15; \
	xmm2 *= xmm9;			/* sin2*i25s				; 16-20 */ \
\
	xmm12 -= xmm3;			/* t8=sin4*r25s - sin2*r34s		; 17-19(16) */ \
	xmm3 = CONST2_P588; \
	xmm3 *= xmm1;			/* sin4*i34s				; 17-21 */ \
\
	xmm8 += xmm11;			/* cos2*i25a + cos4*i34a		; 18-20 */ \
	xmm9 *= CONST2_P588;			/* sin4*i25s				; 18-22 */ \
\
	xmm11 = xptr(srcreg+srcoff+32);	/* Load I1 */ \
	xmm4 += xmm11;			/* outi(0)=i1+i25a+i34a (Final I1)	; 19-21(18) */ \
	xmm1 *= xmm15;			/* sin2*i34s				; 19-23 */ \
\
	xmm7 += xmm10;			/* cos4*i25a + cos2*i34a		; 20-22 */ \
\
	xmm8 += xmm11;			/* t5=cos2*i25a + cos4*i34a + i1	; 21-23 */ \
\
	xmm2 += xmm3;			/* t2=sin2*i25s + sin4*i34s		; 22-24 */ \
	xptr(srcreg+32) = xmm4;		/* Save I1 */ \
\
	xmm7 += xmm11;			/* t7=cos4*i25a + cos2*i34a + i1	; 23-25 */ \
\
	xmm9 -= xmm1;			/* t4=sin4*i25s-sin2*i34s		; 24-26 */ \
\
	xmm10 = xmm5;			/* Copy t6 */ \
	xmm5 += xmm8;			/* outi(1)=t5+t6 (new I2)		; 25-27 */ \
\
	xmm8 -= xmm10;			/* outi(4)=t5-t6 (new I5)		; 26-28 */ \
\
	xmm3 = xmm0;			/* Copy t1 */ \
	xmm0 -= xmm2;			/* outr(1)=t1-t2 (new R2)		; 27-29 */ \
\
	xmm2 += xmm3;			/* outr(4)=t1+t2 (new R5)		; 28-30 */ \
	xmm4 = xptr(screg1+16);		/* cosine/sine */ \
	xmm11 = xmm5;			/* Copy I2 */ \
	xmm5 *= xmm4;			/* B2 = I2 * cosine/sine		; 28-32 */ \
\
	xmm1 = xmm12;			/* Copy t8 */ \
	xmm12 += xmm7;			/* outi(2)=t7+t8 (new I3)		; 29-31 */ \
	xmm10 = xmm8;			/* Copy I5 */ \
	xmm8 *= xmm4;			/* B5 = I5 * cosine/sine		; 29-33 */ \
	xmm7 -= xmm1;			/* outi(3)=t7-t8 (new I4)		; 30-32 */ \
	xmm3 = xmm0;			/* Copy R2 */ \
	xmm0 *= xmm4;			/* A2 = R2 * cosine/sine		; 30-34 */ \
	xmm1 = xmm6;			/* Copy t3 */ \
	xmm6 -= xmm9;			/* outr(2)=t3-t4 (new R3)		; 31-33 */ \
	xmm4 *= xmm2;			/* A5 = R5 * cosine/sine		; 31-35 */ \
	xmm9 += xmm1;			/* outr(3)=t3+t4 (new R4)		; 32-34 */ \
\
	xmm1 = xmm12;			/* Copy I3 */ \
	xmm5 += xmm3;			/* B2 = B2 + R2				; 33-35 */ \
	xmm3 = xptr(screg2+16);		/* cosine/sine */ \
	xmm12 *= xmm3;			/* B3 = I3 * cosine/sine		; 32-36 */ \
\
	xmm8 -= xmm2;			/* B5 = B5 - R5				; 34-36 */ \
	xmm2 = xmm7;			/* Copy I4 */ \
	xmm7 *= xmm3;			/* B4 = I4 * cosine/sine		; 33-37 */ \
\
	xmm0 -= xmm11;			/* A2 = A2 - I2				; 35-37 */ \
	xmm11 = xmm6;			/* Copy R3 */ \
	xmm6 *= xmm3;	 		/* A3 = R3 * cosine/sine		; 34-38 */ \
\
	xmm3 *= xmm9;			/* A4 = R4 * cosine/sine		; 35-39 */ \
\
	xmm4 += xmm10;			/* A5 = A5 + I5				; 36-38 */ \
	xmm10 = xptr(screg1);			/* sine */ \
	xmm5 *= xmm10;			/* B2 = B2 * sine (final I2)		; 36-40 */ \
\
	xmm12 += xmm11;			/* B3 = B3 + R3				; 37-39 */ \
	xmm8 *= xmm10;			/* B5 = B5 * sine (final I5)		; 37-41 */ \
	xmm7 -= xmm9;			/* B4 = B4 - R4				; 38-40 */ \
	xmm0 *= xmm10;			/* A2 = A2 * sine (final R2)		; 38-42 */ \
	xmm6 -= xmm1;			/* A3 = A3 - I3				; 39-41 */ \
	xmm4 *= xmm10;			/* A5 = A5 * sine (final R5)		; 39-43 */ \
	xmm3 += xmm2;			/* A4 = A4 + I4				; 40-42 */ \
\
	xmm11 = xptr(screg2); \
	xmm12 *= xmm11;			/* B3 = B3 * sine (final I3)		; 40-44 */ \
	xmm7 *= xmm11;			/* B4 = B4 * sine (final I4)		; 41-45 */ \
	xptr(srcreg+d1+32) = xmm5;		/* Save I2				; 41 */ \
	xmm6 *= xmm11;			/* A3 = A3 * sine (final R3)		; 42-46 */ \
	xptr(srcreg+4*d1+32) = xmm8;		/* Save I5				; 42 */ \
	xmm3 *= xmm11;			/* A4 = A4 * sine (final R4)		; 43-47 */ \
	xptr(srcreg+d1) = xmm0;		/* Save R2				; 43 */ \
\
	xptr(srcreg+4*d1) = xmm4;		/* Save R5 */ \
	xptr(srcreg+2*d1+32) = xmm12;		/* Save I3 */ \
	xptr(srcreg+3*d1+32) = xmm7;		/* Save I4 */ \
	xptr(srcreg+2*d1) = xmm6;		/* Save R3 */ \
	xptr(srcreg+3*d1) = xmm3;		/* Save R4 */ \
\
	srcreg += srcinc; \

#endif
\
/* */ \
/* ************************************* five-complex-djbunfft variants ****************************************** */ \
/* */ \
\

#define r5_x5cl_five_complex_djbunfft_preload \
	r5_x5cl_2sc_five_complex_djbunfft_preload; \
\
\

#define r5_x5cl_five_complex_djbunfft(srcreg,srcinc,d1,screg,scoff) \
	r5_x5cl_2sc_five_complex_djbunfft(srcreg,srcinc,d1,screg,scoff,screg+32,scoff); \
\
\
/* Same as above except it uses two sin/cos pointers because the first */ \
/* levels of pass 2 was radix-3 which only has one sin/cos value. */ \
\

#define r5_x5cl_2sc_five_complex_djbunfft_preload \
\
\

#define r5_x5cl_2sc_five_complex_djbunfft(srcreg,srcinc,d1,screg1,scoff1,screg2,scoff2) { \
	vec2f64 *sd0 = (vec2f64*)(srcreg); \
	vec2f64 *sd1 = (vec2f64*)(srcreg+(d1)); \
	vec2f64 *sd2 = (vec2f64*)(srcreg+2*(d1)); \
	vec2f64 *sd3 = (vec2f64*)(srcreg+3*(d1)); \
	vec2f64 *sd4 = (vec2f64*)(srcreg+4*(d1)); \
	vec2f64 TEMP5, TEMP6, TEMP7, TEMP8;\
	TEMP5 = sd1[1]; \
	TEMP6 = sd3[1]; \
	r5_x5c_djbunfft_mem(sd0[0],sd0[2],sd1[0],sd1[2],sd2[0],sd2[2],sd3[0],sd3[2],sd4[0],sd4[2], screg1,screg2, sd0[0],sd2[0],sd4[0],sd1[1],sd3[1]); \
	sd0[2] = xmm5;	/* Save I1 */ \
	sd2[2] = xmm6;	/* Save I2 */ \
	sd4[2] = xmm2;	/* Save I3 */ \
	TEMP7 = sd1[3];	/* I2 */ \
	sd1[3] = xmm0;	/* Save I4 */ \
	TEMP8 = sd3[3];	/* I4 */ \
	sd3[3] = xmm4;	/* Save I5 */ \
	r5_x5c_djbunfft_mem(sd0[1],sd0[3],TEMP5,TEMP7,sd2[1],sd2[3],TEMP6,TEMP8,sd4[1],sd4[3], screg1+scoff1,screg2+scoff2,sd1[0],sd3[0],sd0[1],sd2[1],sd4[1]); \
	sd1[2] = xmm5;	/* Save I1 */ \
	sd3[2] = xmm6;	/* Save I2 */ \
	sd0[3] = xmm2;	/* Save I3 */ \
	sd2[3] = xmm0;	/* Save I4 */ \
	sd4[3] = xmm4;	/* Save I5 */ \
	srcreg += srcinc; \
}

/* Do a 5-complex inverse FFT. */ \
/* First we apply twiddle factors to 4 of the 5 input numbers. */ \
/* A 5-complex inverse FFT is like the forward FFT except all the  */ \
/* sin values are negated. */ \
/* Assumes imaginary part of source values is 32 bytes after real part (to workaround MASM line-length limitations) */ \
/* (was: r5_x5c_djbunfft_mem MACRO r1,i1,r2,i2,r3,i3,r4,i4,r5,i5,screg1,screg2,dstr1,dstr2,dstr3,dstr4,dstr5) */ \

#define r5_x5c_djbunfft_mem(R1,I1,R2,I2,R3,I3,R4,I4,R5,I5,screg1,screg2,dstr1,dstr2,dstr3,dstr4,dstr5) { \
	vec2f64 TMP1, TMP2, TMP3, TMP4; \
	xmm1 = R2; \
	xmm0 = xptr(screg1+16); \
	xmm7 = R5; \
	xmm4 = I2; \
	xmm6 = I5; \
	xmm1 = xmm1 * xmm0 + xmm4;	/* A2 = R2 * cosine/sine + I2 */ \
	xmm7 = xmm7 * xmm0 - xmm6;	/* A5 = R5 * cosine/sine - I5 */ \
	xmm4 = xmm4 * xmm0 - R2;	/* B2 = I2 * cosine/sine - R2 */ \
	xmm6 = xmm6 * xmm0 + R5;	/* B5 = I5 * cosine/sine + R5 */ \
	xmm3 = xptr(screg1); \
	xmm1 *= xmm3;			/* A2 = A2 * sine (new R2) */ \
	xmm4 *= xmm3;			/* B2 = B2 * sine (new I2) */ \
	xmm7 *= xmm3;			/* A5 = A5 * sine (new R5) */ \
	xmm6 *= xmm3;			/* B5 = B5 * sine (new I5) */ \
\
	xmm2 = xmm1;			/* Copy R2 */ \
	xmm1 -= xmm7;			/* r25s=r2-r5 */ \
	xmm7 += xmm2;			/* r25a=r2+r5 */ \
	xmm2 = xmm4;			/* Copy I2 */ \
	xmm4 -= xmm6;			/* i25s=i2-i5 */ \
	xmm6 += xmm2;			/* i25a=i2+i5 */ \
\
	TMP1 = xmm1;			/* Save r2-r5 */ \
	TMP2 = xmm6;			/* Save i2+i5 */ \
\
	vec2f64 r4r = R4; \
	vec2f64 r3r = R3; \
	xmm1 = xptr(screg2+16); \
	xmm2 = I4; \
	xmm0 = I3; \
	xmm5 = r4r * xmm1 - xmm2;	/* A4 = R4 * cosine/sine - I4 */ \
	xmm3 = r3r * xmm1 + xmm0;	/* A3 = R3 * cosine/sine + I3 */ \
	xmm2 = xmm2 * xmm1 + r4r;	/* B4 = I4 * cosine/sine + R4 */ \
	xmm0 = xmm0 * xmm1 - r3r;	/* B3 = I3 * cosine/sine - R3 */ \
	xmm6 = xptr(screg2); \
	xmm5 *= xmm6;			/* A4 = A4 * sine (new R4) */ \
	xmm2 *= xmm6;			/* B4 = B4 * sine (new I4) */ \
	xmm3 *= xmm6;			/* A3 = A3 * sine (new R3) */ \
	xmm0 *= xmm6;			/* B3 = B3 * sine (new I3) */ \
\
	xmm1 = xmm3;			/* Copy R3 */ \
	xmm3 -= xmm5;			/* r34s=r3-r4 */ \
	xmm5 += xmm1;			/* r34a=r3+r4 */ \
	TMP3 = xmm3; \
	xmm6 = xmm0;			/* Copy I3 */ \
	xmm0 -= xmm2;			/* i34s=i3-i4 */ \
	xmm2 += xmm6;			/* i34a=i3+i4 */ \
	TMP4 = xmm2; \
\
	xmm1 = xmm7 * CONST2_P309;	/* cos2*r25a */ \
	xmm1 += xmm5 * CONST2_M809;	/* cos2*r25a + cos4*r34a */ \
	xmm3 = xmm4 * CONST2_P951;	/* sin2*i25s */ \
	xmm3 += xmm0 * CONST2_P588;	/* t2=sin2*i25s + sin4*i34s */ \
	xmm2 = xmm7 * CONST2_M809;	/* cos4*r25a */ \
	xmm2 += xmm5 * CONST2_P309;			/* cos4*r25a+cos2*r34a */ \
	xmm4 = xmm4 * CONST2_P588;	/* sin4*i25s */ \
	xmm4 -= xmm0 * CONST2_P951;			/* t4=sin4*i25s-sin2*i34s */ \
	xmm7 += xmm5;			/* r25a + r34a */ \
	xmm6 = R1; \
	xmm1 += xmm6;			/* t1=cos2*r25a + cos4*r34a + r1 */ \
	xmm7 += xmm6;			/* outr(0) = r1 + r25a + r34a */ \
	xmm2 += xmm6;			/* t3=cos4*r25a+cos2*r34a+r1 */ \
	xmm6 = xmm1; \
	xmm1 -= xmm3;			/* outr(4)=t1-t2 */ \
	xmm3 += xmm6;			/* outr(1)=t1+t2 */ \
	xmm5 = xmm2; \
	xmm2 -= xmm4;			/* outr(3)=t3-t4 */ \
	xmm4 += xmm5;			/* outr(2)=t3+t4 */ \
\
	dstr1 = xmm7; \
	dstr2 = xmm3;			/* Save new r2 */ \
	dstr3 = xmm4;			/* Save new r3 */ \
	dstr4 = xmm2;			/* Save new r4 */ \
	dstr5 = xmm1;			/* Save new r5 */ \
\
	xmm0 = TMP1;			/* r25s=r2-r5 */ \
	xmm2 = TMP2;			/* i25a=i2+i5 */ \
	xmm1 = TMP3;			/* r34s=r3-r4 */ \
	xmm3 = TMP4;			/* i34a=i3+i4 */ \
\
	xmm5 = xmm2; \
	xmm5 += xmm3;			/* i25a+i34a */ \
	xmm6 = xmm2 * CONST2_P309;	/* cos2*i25a */ \
	xmm6 += xmm3 * CONST2_M809;	/* cos2*i25a + cos4*i34a */ \
	xmm2 = xmm2 * CONST2_M809;	/* cos4*i25a */ \
	xmm2 += xmm3 * CONST2_P309;	/* cos4*i25a + cos2*i34a */ \
	xmm4 = xmm0 * CONST2_P951;	/* sin2*r25s */ \
	xmm4 += xmm1 * CONST2_P588;	/* t6=sin2*r25s + sin4*r34s */ \
	xmm0 = xmm0 * CONST2_P588;	/* sin4*r25s */ \
	xmm0 -= xmm1 * CONST2_P951;	/* t8=sin4*r25s - sin2*r34s */ \
	xmm3 = I1; /*i1 */ \
	xmm6 += xmm3;			/* t5=cos2*i25a + cos4*i34a + i1 */ \
	xmm2 += xmm3;			/* t7=cos4*i25a + cos2*i34a + i1 */ \
\
	xmm7 = xmm6; \
	xmm6 -= xmm4;			/* outi(1)=t5-t6 */ \
	xmm4 += xmm7;			/* outi(4)=t5+t6 */ \
	xmm1 = xmm2; \
	xmm2 -= xmm0;			/* outi(2)=t7-t8 */ \
	xmm0 += xmm1;			/* outi(3)=t7+t8 */ \
	xmm5 += xmm3;			/* outi(0)=i1+i25a+i34a */ \
}

#if 0 //DEF X86_64;


#define r5_x5cl_2sc_five_complex_djbunfft_preload \
	xmm14 = CONST2_P309; \
	xmm15 = CONST2_M809; \
\
\

#define r5_x5cl_2sc_five_complex_djbunfft(srcreg,srcinc,d1,screg1,scoff1,screg2,scoff2) \
	xmm0 = xptr(srcreg+d1);		/* Load R2 */ \
	xmm1 = xptr(screg1+16);		/* cosine/sine */ \
	xmm2 = xmm0;			/* Copy R2 */ \
	xmm0 *= xmm1;			/* A2 = R2 * cosine/sine		; 1-5 */ \
	xmm3 = xptr(srcreg+4*d1);		/* Load R5 */ \
	xmm4 = xmm1;			/* Copy cosine/sine */ \
	xmm1 *= xmm3;			/* A5 = R5 * cosine/sine		; 2-6 */ \
	xmm5 = xptr(srcreg+3*d1);		/* Load R4 */ \
	xmm6 = xptr(screg2+16);		/* cosine/sine */ \
	xmm7 = xmm5;			/* Copy R4 */ \
	xmm5 *= xmm6;			/* A4 = R4 * cosine/sine		; 3-7 */ \
	xmm8 = xptr(srcreg+2*d1);		/* Load R3 */ \
	xmm9 = xmm6;			/* Copy cosine/sine */ \
	xmm6 *= xmm8;			/* A3 = R3 * cosine/sine		; 4-8 */ \
	xmm10 = xptr(srcreg+d1+32);		/* Load I2 */ \
	xmm11 = xmm4;			/* Copy cosine/sine */ \
	xmm4 *= xmm10;			/* B2 = I2 * cosine/sine		; 5-9 */ \
	xmm0 += xmm10;			/* A2 = A2 + I2				; 6-8			avail 12,10 */ \
	xmm12 = xptr(srcreg+4*d1+32);		/* Load I5 */ \
	xmm11 *= xmm12;			/* B5 = I5 * cosine/sine		; 6-10			avail 10 */ \
	xmm1 -= xmm12;			/* A5 = A5 - I5				; 7-9			avail 10,12 */ \
	xmm10 = xptr(srcreg+3*d1+32);		/* Load I4 */ \
	xmm12 = xmm9;			/* Copy cosine/sine */ \
	xmm9 *= xmm10;			/* B4 = I4 * cosine/sine		; 7-11			avail none */ \
	xmm5 -= xmm10;			/* A4 = A4 - I4				; 8-10 */ \
	xmm10 = xptr(srcreg+2*d1+32);		/* Load I3 */ \
	xmm12 *= xmm10;			/* B3 = I3 * cosine/sine		; 8-12 */ \
	xmm6 += xmm10;			/* A3 = A3 + I3				; 9-11 */ \
	xmm10 = xptr(screg1);			/* sine */ \
	xmm0 *= xmm10;			/* A2 = A2 * sine (new R2)		; 9-13 */ \
	xmm4 -= xmm2;			/* B2 = B2 - R2				; 10-12			avail 2 */ \
	xmm1 *= xmm10;			/* A5 = A5 * sine (new R5)		; 10-14 */ \
	xmm11 += xmm3;			/* B5 = B5 + R5				; 11-13			avail 2,3 */ \
	xmm2 = xptr(screg2);			/* sine */ \
	xmm5 *= xmm2;			/* A4 = A4 * sine (new R4)		; 11-15			avail 3 */ \
	xmm9 += xmm7;			/* B4 = B4 + R4				; 12-14			avail 3,7 */ \
	xmm6 *= xmm2;			/* A3 = A3 * sine (new R3)		; 12-16 */ \
	xmm12 -= xmm8;			/* B3 = B3 - R3				; 13-15			avail 3,7,8 */ \
	xmm4 *= xmm10;			/* B2 = B2 * sine (new I2)		; 13-17 */ \
	xmm11 *= xmm10;			/* B5 = B5 * sine (new I5)		; 14-18			avail 3,7,8,10 */ \
	xmm3 = xmm0;			/* Copy R2 */ \
	xmm0 += xmm1;			/* r25a=r2+r5				; 15-17			avail 7,8,10 */ \
	xmm9 *= xmm2;			/* B4 = B4 * sine (new I4)		; 15-19 */ \
	xmm3 -= xmm1;			/* r25s=r2-r5				; 16-18			avail 7,8,10,1 */ \
	xmm12 *= xmm2;			/* B3 = B3 * sine (new I3)		; 16-20			avail 7,8,10,1,2 */ \
	xmm7 = xmm5;			/* Copy R4 */ \
	xmm5 += xmm6;			/* r34a=r3+r4				; 17-19 */ \
	xmm6 -= xmm7;			/* r34s=r3-r4				; 18-20			avail 8,10,1,2,7 */ \
	xmm8 = xmm14; \
	xmm8 *= xmm0;			/* cos2*r25a				; 18-22			avail 10,1,2,7 */ \
	xmm10 = xmm4;			/* Copy I2 */ \
	xmm4 += xmm11;			/* i25a=i2+i5				; 19-21			avail 1,2,7 */ \
	xmm1 = xmm15; \
	xmm1 *= xmm0;			/* cos4*r25a				; 19-23			avail 2,7 */ \
	xmm10 -= xmm11;			/* i25s=i2-i5				; 20-22			avail 2,7,11 */ \
	xmm2 = xmm15; \
	xmm2 *= xmm5;			/* cos4*r34a				; 20-24			avail 7,11 */ \
	xmm0 += xmm5;			/* r25a + r34a				; 21-23 */ \
	xmm5 *= xmm14;			/* cos2*r34a				; 21-25 */ \
	xmm7 = xmm9;			/* Copy I4 */ \
	xmm9 += xmm12;			/* i34a=i3+i4				; 22-24			avail 11 */ \
	xmm11 = CONST2_P951; \
	xmm11 *= xmm3;			/* sin2*r25s				; 22-26			avail none */ \
	xmm12 -= xmm7;			/* i34s=i3-i4				; 23-25			avail 7 */ \
	xmm3 *= CONST2_P588;			/* sin4*r25s				; 23-27 */ \
\
	xmm13 = xptr(srcreg);			/* Load R1  */ \
	xmm0 += xmm13;			/* outr(0) = r1 + r25a + r34a		; 24-26 */ \
	xmm7 = CONST2_P588; \
	xmm7 *= xmm6;			/* sin4*r34s				; 24-28			avail none storable 0 */ \
	xmm8 += xmm2;			/* cos2*r25a + cos4*r34a		; 25-27			avail 2 storable 0 */ \
	xmm6 *= CONST2_P951;			/* sin2*r34s				; 25-29 */ \
	xmm1 += xmm5;			/* cos4*r25a+cos2*r34a			; 26-28			avail 2,5 storable 0 */ \
	xmm2 = xmm14; \
	xmm2 *= xmm4;			/* cos2*i25a				; 26-30			avail 5 storable 0 */ \
	xmm5 = xmm4;			/* Copy i25a */ \
	xmm4 += xmm9;			/* i25a+i34a				; 27-29 */ \
	xptr(srcreg) = xmm0;			/* Save R1				; 27			avail 0 */ \
	xmm0 = xmm15; \
	xmm0 *= xmm9;			/* cos4*i34a				; 27-31			avail none */ \
	xmm8 += xmm13;			/* t1=cos2*r25a + cos4*r34a + r1	; 28-30 */ \
	xmm5 *= xmm15;			/* cos4*i25a				; 28-32 */ \
	xmm1 += xmm13;			/* t3=cos4*r25a+cos2*r34a+r1		; 29-31			avail 13 */ \
	xmm9 *= xmm14;			/* cos2*i34a				; 29-33 */ \
	xmm11 += xmm7;			/* t6=sin2*r25s + sin4*r34s		; 30-32 */ \
	xmm7 = CONST2_P951; \
	xmm7 *= xmm10;			/* sin2*i25s				; 30-34 */ \
	xmm3 -= xmm6;			/* t8=sin4*r25s - sin2*r34s		; 31-33 */ \
	xmm6 = CONST2_P588; \
	xmm6 *= xmm12;			/* sin4*i34s				; 31-35 */ \
	xmm2 += xmm0;			/* cos2*i25a + cos4*i34a		; 32-34			avail 0 */ \
	xmm10 *= CONST2_P588;			/* sin4*i25s				; 32-36 */ \
\
	xmm0 = xptr(srcreg+32);		/* Load I1 */ \
	xmm4 += xmm0;			/* outi(0)=i1+i25a+i34a			; 33-35 (30) */ \
	xmm12 *= CONST2_P951;			/* sin2*i34s				; 33-37 */ \
	xmm5 += xmm9;			/* cos4*i25a + cos2*i34a		; 34-36 */ \
	xmm9 = xptr(srcreg+d1+16);		/* Load next R2 */ \
	xmm2 += xmm0;			/* t5=cos2*i25a + cos4*i34a + i1	; 35-37 */ \
	xptr(srcreg+32) = xmm4;		/* Save I1				; 36 */ \
	xmm7 += xmm6;			/* t2=sin2*i25s + sin4*i34s		; 36-38 */ \
	xmm6 = xptr(srcreg+3*d1+16);		/* Load next R4 */ \
	xmm5 += xmm0;			/* t7=cos4*i25a + cos2*i34a + i1	; 37-39			avail 4,0 */ \
	xmm13 = xptr(srcreg+d1+48);		/* Load next I2 */ \
	xmm10 -= xmm12;			/* t4=sin4*i25s-sin2*i34s		; 38-40			avail 4,0,12 */ \
	xmm4 = xptr(srcreg+3*d1+48);		/* Load next I4 */ \
\
	xmm12 = xmm11;			/* Copy t6 */ \
	xmm11 += xmm2;			/* outi(4)=t5+t6			; 39-41 */ \
	xmm0 = xptr(screg1+scoff1+16);	/* next cosine/sine */ \
	xptr(srcreg+3*d1+48) = xmm11;		/* Save I5				; 42 */ \
	xmm2 -= xmm12;			/* outi(1)=t5-t6			; 40-42 */ \
	xmm12 = xmm8;			/* Copy t1 */ \
	xptr(srcreg+2*d1+32) = xmm2;		/* Save I2				; 43 */ \
	xmm8 -= xmm7;			/* outr(4)=t1-t2			; 41-43 */ \
	xmm2 = xptr(srcreg+4*d1+16);		/* Load next R5 */ \
	xptr(srcreg+3*d1+16) = xmm8;		/* Save R5				; 44 */ \
\
	xmm7 += xmm12;			/* outr(1)=t1+t2			; 42-44 */ \
	xmm11 = xmm9;			/* Copy R2 */ \
	xmm9 *= xmm0;			/* A2 = R2 * cosine/sine		; 1-5 (42-46) */ \
	xptr(srcreg+2*d1) = xmm7;		/* Save R2				; 45 */ \
\
	xmm12 = xmm3;			/* Copy t8 */ \
	xmm3 += xmm5;			/* outi(3)=t7+t8			; 43-45 */ \
	xmm8 = xmm0;			/* Copy cosine/sine */ \
	xmm0 *= xmm2;			/* A5 = R5 * cosine/sine		; 2-6 (43-47) */ \
	xptr(srcreg+d1+48) = xmm3;		/* Save I4				; 46 */ \
\
	xmm5 -= xmm12;			/* outi(2)=t7-t8			; 44-46 */ \
	xmm7 = xptr(screg2+scoff2+16);	/* next cosine/sine */ \
	xmm3 = xmm6;			/* Copy R4 */ \
	xmm6 *= xmm7;			/* A4 = R4 * cosine/sine		; 3-7 (44-48) */ \
	xptr(srcreg+4*d1+32) = xmm5;		/* Save	I3				; 47 */ \
\
	xmm12 = xmm1;			/* Copy t3 */ \
	xmm1 -= xmm10;			/* outr(3)=t3-t4			; 45-47 */ \
	xmm5 = xptr(srcreg+2*d1+16);		/* Load R3 */ \
	xptr(srcreg+d1+16) = xmm1;		/* Save R4				; 48 */ \
	xmm1 = xmm7;			/* Copy cosine/sine */ \
	xmm7 *= xmm5;			/* A3 = R3 * cosine/sine		; 4-8 (45-49) */ \
\
	xmm10 += xmm12;			/* outr(2)=t3+t4			; 46-48 */ \
	xmm12 = xmm8;			/* Copy cosine/sine */ \
	xmm8 *= xmm13;			/* B2 = I2 * cosine/sine		; 5-9 (46-50) */ \
	xptr(srcreg+4*d1) = xmm10;		/* Save R3				; 49 */ \
\
	xmm9 += xmm13;			/* A2 = A2 + I2				; 6-8 */ \
	xmm10 = xptr(srcreg+4*d1+48);		/* Load I5 */ \
	xmm12 *= xmm10;			/* B5 = I5 * cosine/sine		; 6-10 */ \
	xmm0 -= xmm10;			/* A5 = A5 - I5				; 7-9 */ \
	xmm10 = xmm1;			/* Copy cosine/sine */ \
	xmm1 *= xmm4;			/* B4 = I4 * cosine/sine		; 7-11 */ \
	xmm6 -= xmm4;			/* A4 = A4 - I4				; 8-10 */ \
	xmm4 = xptr(srcreg+2*d1+48);		/* Load I3 */ \
	xmm10 *= xmm4;			/* B3 = I3 * cosine/sine		; 8-12 */ \
	xmm7 += xmm4;			/* A3 = A3 + I3				; 9-11 */ \
	xmm4 = xptr(screg1+scoff1);		/* sine */ \
	xmm9 *= xmm4;			/* A2 = A2 * sine (new R2)		; 9-13 */ \
	xmm8 -= xmm11;			/* B2 = B2 - R2				; 10-12 */ \
	xmm0 *= xmm4;			/* A5 = A5 * sine (new R5)		; 10-14 */ \
	xmm12 += xmm2;			/* B5 = B5 + R5				; 11-13 */ \
	xmm11 = xptr(screg2+scoff2);		/* sine */ \
	xmm6 *= xmm11;			/* A4 = A4 * sine (new R4)		; 11-15 */ \
	xmm1 += xmm3;			/* B4 = B4 + R4				; 12-14 */ \
	xmm7 *= xmm11;			/* A3 = A3 * sine (new R3)		; 12-16 */ \
	xmm10 -= xmm5;			/* B3 = B3 - R3				; 13-15 */ \
	xmm8 *= xmm4;			/* B2 = B2 * sine (new I2)		; 13-17 */ \
	xmm12 *= xmm4;			/* B5 = B5 * sine (new I5)		; 14-18 */ \
	xmm2 = xmm9;			/* Copy R2 */ \
	xmm9 += xmm0;			/* r25a=r2+r5				; 15-17 */ \
	xmm1 *= xmm11;			/* B4 = B4 * sine (new I4)		; 15-19 */ \
	xmm2 -= xmm0;			/* r25s=r2-r5				; 16-18 */ \
	xmm10 *= xmm11;			/* B3 = B3 * sine (new I3)		; 16-20 */ \
	xmm3 = xmm6;			/* Copy R4 */ \
	xmm6 += xmm7;			/* r34a=r3+r4				; 17-19 */ \
	xmm7 -= xmm3;			/* r34s=r3-r4				; 18-20 */ \
	xmm5 = xmm14; \
	xmm5 *= xmm9;			/* cos2*r25a				; 18-22 */ \
	xmm4 = xmm8;			/* Copy I2 */ \
	xmm8 += xmm12;			/* i25a=i2+i5				; 19-21 */ \
	xmm0 = xmm15; \
	xmm0 *= xmm9;			/* cos4*r25a				; 19-23 */ \
	xmm4 -= xmm12;			/* i25s=i2-i5				; 20-22 */ \
	xmm11 = xmm15; \
	xmm11 *= xmm6;			/* cos4*r34a				; 20-24 */ \
	xmm9 += xmm6;			/* r25a + r34a				; 21-23 */ \
	xmm6 *= xmm14;			/* cos2*r34a				; 21-25 */ \
	xmm3 = xmm1;			/* Copy I4 */ \
	xmm1 += xmm10;			/* i34a=i3+i4				; 22-24 */ \
	xmm12 = CONST2_P951; \
	xmm12 *= xmm2;			/* sin2*r25s				; 22-26 */ \
	xmm10 -= xmm3;			/* i34s=i3-i4				; 23-25 */ \
	xmm2 *= CONST2_P588;			/* sin4*r25s				; 23-27 */ \
\
	xmm13 = xptr(srcreg+16);		/* Load R1  */ \
	xmm9 += xmm13;			/* outr(0) = r1 + r25a + r34a		; 24-26 */ \
	xmm3 = CONST2_P588; \
	xmm3 *= xmm7;			/* sin4*r34s				; 24-28 */ \
	xmm5 += xmm11;			/* cos2*r25a + cos4*r34a		; 25-27 */ \
	xmm7 *= CONST2_P951;			/* sin2*r34s				; 25-29 */ \
	xmm0 += xmm6;			/* cos4*r25a+cos2*r34a			; 26-28 */ \
	xmm11 = xmm14; \
	xmm11 *= xmm8;			/* cos2*i25a				; 26-30 */ \
	xmm6 = xmm8;			/* Copy i25a */ \
	xmm8 += xmm1;			/* i25a+i34a				; 27-29 */ \
	xptr(srcreg+d1) = xmm9;		/* Save R1				; 27 */ \
	xmm9 = xmm15; \
	xmm9 *= xmm1;			/* cos4*i34a				; 27-31 */ \
	xmm5 += xmm13;			/* t1=cos2*r25a + cos4*r34a + r1	; 28-30 */ \
	xmm6 *= xmm15;			/* cos4*i25a				; 28-32 */ \
	xmm0 += xmm13;			/* t3=cos4*r25a+cos2*r34a+r1		; 29-31 */ \
	xmm1 *= xmm14;			/* cos2*i34a				; 29-33 */ \
	xmm12 += xmm3;			/* t6=sin2*r25s + sin4*r34s		; 30-32 */ \
	xmm3 = CONST2_P951; \
	xmm3 *= xmm4;			/* sin2*i25s				; 30-34 */ \
	xmm2 -= xmm7;			/* t8=sin4*r25s - sin2*r34s		; 31-33 */ \
	xmm7 = CONST2_P588; \
	xmm7 *= xmm10;			/* sin4*i34s				; 31-35 */ \
	xmm11 += xmm9;			/* cos2*i25a + cos4*i34a		; 32-34 */ \
	xmm4 *= CONST2_P588;			/* sin4*i25s				; 32-36 */ \
\
	xmm9 = xptr(srcreg+48);		/* Load I1 */ \
	xmm8 += xmm9;			/* outi(0)=i1+i25a+i34a			; 33-35 (30) */ \
	xmm10 *= CONST2_P951;			/* sin2*i34s				; 33-37 */ \
	xmm6 += xmm1;			/* cos4*i25a + cos2*i34a		; 34-36 */ \
	xmm11 += xmm9;			/* t5=cos2*i25a + cos4*i34a + i1	; 35-37 */ \
	xmm3 += xmm7;			/* t2=sin2*i25s + sin4*i34s		; 36-38 */ \
	xptr(srcreg+d1+32) = xmm8;		/* Save I1				; 36 */ \
	xmm6 += xmm9;			/* t7=cos4*i25a + cos2*i34a + i1	; 37-39 */ \
	xmm4 -= xmm10;			/* t4=sin4*i25s-sin2*i34s		; 38-40 */ \
\
	xmm1 = xmm12;			/* Copy t6 */ \
	xmm12 += xmm11;			/* outi(4)=t5+t6 */ \
	xmm11 -= xmm1;			/* outi(1)=t5-t6 */ \
	xmm7 = xmm5;			/* Copy t1 */ \
	xmm5 -= xmm3;			/* outr(4)=t1-t2 */ \
	xmm3 += xmm7;			/* outr(1)=t1+t2 */ \
	xmm8 = xmm2;			/* Copy t8 */ \
	xmm2 += xmm6;			/* outi(3)=t7+t8 */ \
	xmm6 -= xmm8;			/* outi(2)=t7-t8 */ \
	xmm9 = xmm0;			/* Copy t3 */ \
	xmm0 -= xmm4;			/* outr(3)=t3-t4 */ \
	xmm4 += xmm9;			/* outr(2)=t3+t4 */ \
\
	xptr(srcreg+4*d1+48) = xmm12;		/* Save I5 */ \
	xptr(srcreg+3*d1+32) = xmm11;		/* Save I2 */ \
	xptr(srcreg+4*d1+16) = xmm5;		/* Save R5 */ \
	xptr(srcreg+3*d1) = xmm3;		/* Save R2 */ \
	xptr(srcreg+2*d1+48) = xmm2;		/* Save I4 */ \
	xptr(srcreg+48) = xmm6;		/* Save	I3 */ \
	xptr(srcreg+2*d1+16) = xmm0;		/* Save R4 */ \
	xptr(srcreg+16) = xmm4;		/* Save R3 */ \
\
	srcreg += srcinc; \
\

#endif
\
/* */ \
/* ************************************* ten-reals-fft variants ****************************************** */ \
/* */ \
\

#if 0 // unused
#define r5_x5cl_ten_reals_five_complex_djbfft(srcreg,srcinc,d1,screg1,screg2) \
	d2=2*d1; \
	d3=3*d1; \
	d4=4*d1; \
\
	r5_x5c_djbfft_mem(xptr(srcreg+32),xptr(srcreg+48),xptr(srcreg+d1+32),xptr(srcreg+d1+48),xptr(srcreg+d2+32),xptr(srcreg+d2+48),xptr(srcreg+d3+32),xptr(srcreg+d3+48),xptr(srcreg+d4+32),xptr(srcreg+d4+48),screg1,screg1+32,xptr(srcreg+d2+32),xptr(srcreg+d2+48),xptr(srcreg+d4+32)); \
\
	xmm7 = xptr(srcreg+d3);		/* Load R4 */ \
	TMP4 = xmm7; \
	xmm7 = xptr(srcreg+d3+16);		/* Load R9 */ \
	TMP5 = xmm7; \
	xmm7 = xptr(srcreg+d4);		/* Load R5 */ \
	TMP6 = xmm7; \
	xmm7 = xptr(srcreg+d4+16);		/* Load R10 */ \
	TMP7 = xmm7; \
\
	xptr(srcreg+d3) = xmm1;		/* Save R2 */ \
	xptr(srcreg+d3+16) = xmm0;		/* Save I2 */ \
	xptr(srcreg+d3+32) = xmm3;		/* Save R3 */ \
	xptr(srcreg+d3+48) = xmm4;		/* Save I3 */ \
	xptr(srcreg+d4) = xmm5;		/* Save R4 */ \
	xptr(srcreg+d4+16) = xmm2;		/* Save I4 */ \
/*	xstore	[srcreg+d4+32], xmm7		;; Save R5 */ \
	xptr(srcreg+d4+48) = xmm6;		/* Save I5 */ \
\
	r5_x10r_fft_mem(xptr(srcreg),xptr(srcreg+d1),xptr(srcreg+d2),TMP4,TMP6,xptr(srcreg+16),xptr(srcreg+d1+16),xptr(srcreg+d2+16),TMP5,TMP7,screg1,screg2,xptr(srcreg),xptr(srcreg+16)); \
	xptr(srcreg+32) = xmm2;		/* Save R2 */ \
	xptr(srcreg+48) = xmm4;		/* Save I2 */ \
	xptr(srcreg+d1) = xmm1;		/* Save R3 */ \
	xptr(srcreg+d1+16) = xmm0;		/* Save I3 */ \
	xptr(srcreg+d1+32) = xmm5;		/* Save R4 */ \
	xptr(srcreg+d1+48) = xmm6;		/* Save I4 */ \
	xptr(srcreg+d2) = xmm3;		/* Save R5 */ \
	xptr(srcreg+d2+16) = xmm7;		/* Save I5 */ \
\
	srcreg += srcinc; \
\
\
/* To calculate a 10-reals FFT (in a shorthand notation): */ \
/* r1 + r2 + r3 + r4 + r5 + r6 + r7 + r8 + r9 + r10	*  w^0000000000 */ \
/* r1 + r2 + r3 + r4 + r5 + r6 + r7 + r8 + r9 + r10	*  w^0123456789 */ \
/* r1 + r2 + r3 + r4 + r5 + r6 + r7 + r8 + r9 + r10	*  w^0246802468 */ \
/* r1 + r2 + r3 + r4 + r5 + r6 + r7 + r8 + r9 + r10	*  w^0369258147 */ \
/* r1 + r2 + r3 + r4 + r5 + r6 + r7 + r8 + r9 + r10	*  w^0482604826 */ \
/* r1 + r2 + r3 + r4 + r5 + r6 + r7 + r8 + r9 + r10	*  w^0505050505 */ \
/* r1 + r2 + r3 + r4 + r5 + r6 + r7 + r8 + r9 + r10	*  w^0628406284 */ \
/* r1 + r2 + r3 + r4 + r5 + r6 + r7 + r8 + r9 + r10	*  w^0741852963 */ \
/* r1 + r2 + r3 + r4 + r5 + r6 + r7 + r8 + r9 + r10	*  w^0864208642 */ \
/* r1 + r2 + r3 + r4 + r5 + r6 + r7 + r8 + r9 + r10	*  w^0987654321 */ \
/* Noting that w^5 = -1 and that Hermetian symmetry means we won't need */ \
/* to calculate the last 5 rows: */ \
/* r1 + r2 + r3 + r4 + r5 + r6 + r7 + r8 + r9 + r10	*  w^0000000000 */ \
/* r1 + r2 + r3 + r4 + r5 - r6 - r7 - r8 - r9 - r10	*  w^0123401234 */ \
/* r1 + r2 + r3 - r4 - r5 + r6 + r7 + r8 - r9 - r10	*  w^0241302413 */ \
/* r1 + r2 - r3 - r4 + r5 - r6 - r7 + r8 + r9 - r10	*  w^0314203142 */ \
/* r1 + r2 - r3 + r4 - r5 + r6 + r7 - r8 + r9 - r10	*  w^0432104321 */ \
/* Reorganize into odds and evens, that is */ \
/* r1 + r3 + r5 + r7 + r9  * powers + (r2 + r4 + r6 + r8 + r10) * powers */ \
/* Thus: */ \
/* r1 + r3 + r5 + r7 + r9	*  w^00000	+ r2 + r4 + r6 + r8 + r10	*  w^00000 */ \
/* r1 + r3 + r5 - r7 - r9	*  w^02413	+ r2 + r4 - r6 - r8 - r10	*  w^13024 */ \
/* r1 + r3 - r5 + r7 - r9	*  w^04321	+ r2 - r4 + r6 + r8 - r10	*  w^21043 */ \
/* r1 - r3 + r5 - r7 + r9	*  w^01234	+ r2 - r4 - r6 + r8 - r10	*  w^34012 */ \
/* r1 - r3 - r5 + r7 + r9	*  w^03142	+ r2 + r4 + r6 - r8 - r10	*  w^42031 */ \
/* Apply the sin/cos values: */ \
/* w^1/10 = .809 + .588i */ \
/* w^2/10 = .309 + .951i */ \
/* w^3/10 = -.309 + .951i */ \
/* w^4/10 = -.809 + .588i */ \
/* reals: */ \
/* r1 + r3 + r5 + r7 + r9			+ r2 + r4 + r6 + r8 + r10 */ \
/* r1 + .309r3 - .809r5 - .809r7 + .309r9	+ .809r2 - .309r4 - r6 - .309r8 + .809r10 */ \
/* r1 - .809r3 + .309r5 + .309r7 - .809r9	+ .309r2 - .809r4 + r6 - .809r8 + .309r10 */ \
/* r1 - .809r3 + .309r5 + .309r7 - .809r9	+ -.309r2 + .809r4 - r6 + .809r8 - .309r10 */ \
/* r1 + .309r3 - .809r5 - .809r7 + .309r9	+ -.809r2 + .309r4 + r6 + .309r8 - .809r10 */ \
/* imaginarys: */ \
/* 0						+ 0 */ \
/*  + .951r3 + .588r5 - .588r7 - .951r9		+ .588r2 + .951r4 - .951r8 - .588r10 */ \
/*  + .588r3 - .951r5 + .951r7 - .588r9		+ .951r2 - .588r4 + .588r8 - .951r10 */ \
/*  - .588r3 + .951r5 - .951r7 + .588r9		+ .951r2 - .588r4 + .588r8 - .951r10 */ \
/*  - .951r3 - .588r5 + .588r7 + .951r9		+ .588r2 + .951r4 - .951r8 - .588r10 */ \
/* Further simplifying: */ \
/* reals: */ \
/* r1 + r3 + r5 + r7 + r9		+ r2 + r4 + r6 + r8 + r10 */ \
/* r1 + .309(r3+r9) - .809(r5+r7)	+ .809(r2+r10) - .309(r4+r8) - r6 */ \
/* r1 - .809(r3+r9) + .309(r5+r7)	+ .309(r2+r10) - .809(r4+r8) + r6 */ \
/* r1 - .809(r3+r9) + .309(r5+r7)	+ -.309(r2+r10) + .809(r4+r8) - r6 */ \
/* r1 + .309(r3+r9) - .809(r5+r7)	+ -.809(r2+r10) + .309(r4+r8) + r6 */ \
/* imaginarys: */ \
/* 0					+ 0 */ \
/*  + .951(r3-r9) + .588(r5-r7)		+ .588(r2-r10) + .951(r4-r8) */ \
/*  + .588(r3-r9) - .951(r5-r7)		+ .951(r2-r10) - .588(r4-r8) */ \
/*  - .588(r3-r9) + .951(r5-r7)		+ .951(r2-r10) - .588(r4-r8) */ \
/*  - .951(r3-r9) - .588(r5-r7)		+ .588(r2-r10) + .951(r4-r8) */ \
\

#define r5_x10r_fft_mem(mem1,mem2,mem3,mem4,mem5,mem6,mem7,mem8,mem9,mem10,screg1,screg2,dst1,dst2) { \
	vec2f64 TMP1, TMP2; \
	xmm0 = mem3;		/* R3 */ \
	xmm0 += mem9;		/* R3+R9 */ \
	xmm1 = mem5;		/* R5 */ \
	xmm1 += mem7;		/* R5+R7 */ \
	xmm4 = CONST2_P309; \
	xmm4 *= xmm0;		/* new oddR2 = .309*(R3+R9) */ \
	xmm7 = mem1;		/* R1 */ \
	xmm4 += xmm7;		/* new oddR2 += R1 */ \
	xmm5 = CONST2_M809; \
	xmm5 *= xmm0;		/* new oddR3 = -.809*(R3+R9) */ \
	xmm0 += xmm7;		/* new oddR1 = R1+R3+R9 */ \
	xmm5 += xmm7;		/* new oddR3 += R1 */ \
	xmm7 = CONST2_M809; \
	xmm7 *= xmm1;		/* -.809*(R5+R7) */ \
	xmm0 += xmm1;		/* new oddR1 = R1+R3+R5+R7+R9 */ \
	xmm1 *= CONST2_P309;		/* .309*(R5+R7) */ \
	xmm4 += xmm7;		/* new oddR2 += -.809*(R5+R7) */ \
	xmm5 += xmm1;		/* new oddR3 += .309*(R5+R7) */ \
	dst1 = xmm0;		/* Save oddR1 */ \
\
	xmm0 = mem2;		/* R2 */ \
	xmm0 += mem10;		/* R2+R10 */ \
	xmm1 = mem4;		/* R4 */ \
	xmm1 += mem8;		/* R4+R8 */ \
	xmm2 = CONST2_P309; \
	xmm2 *= xmm0;		/* new evenR3 = .309*(R2+R10) */ \
	xmm7 = mem6;		/* R6 */ \
	xmm2 += xmm7;		/* new evenR3 += R6 */ \
	xmm3 = CONST2_M809; \
	xmm3 *= xmm0;		/* new evenR5 = -.809*(R2+R10) */ \
	xmm0 += xmm7;		/* new evenR1 = R6+R2+R10 */ \
	xmm3 += xmm7;		/* new evenR5 += R6 */ \
	xmm7 = CONST2_M809; \
	xmm7 *= xmm1;		/* -.809*(R4+R8) */ \
	xmm0 += xmm1;		/* new evenR1 = R6+R2+R4+R8+R10 */ \
	xmm1 *= CONST2_P309;		/* .309*(R4+R8) */ \
	xmm2 += xmm7;		/* new evenR3 += -.809*(R4+R8) */ \
	xmm3 += xmm1;		/* new evenR5 += .309*(R4+R8) */ \
	dst2 = xmm0;		/* Save evenR1 */ \
\
	xmm6 = xmm4;		/* Copy oddR2 */ \
	xmm4 -= xmm3;		/* New R2 = oddR2 - evenR5 */ \
	xmm3 += xmm6;		/* New R5 = oddR2 + evenR5 */ \
\
	xmm6 = xmm5;		/* Copy oddR3 */ \
	xmm5 -= xmm2;		/* New R4 = oddR3 - evenR3 */ \
	xmm2 += xmm6;		/* New R3 = oddR3 + evenR3 */ \
\
	TMP1 = xmm4;		/* Temporarily save new R2 */ \
	TMP2 = xmm2;		/* Temporarily save new R3 */ \
\
	xmm0 = mem3;		/* R3 */ \
	xmm0 -= mem9;		/* R3-R9 */ \
	xmm1 = mem5;		/* R5 */ \
	xmm1 -= mem7;		/* R5-R7 */ \
	xmm4 = CONST2_P951; \
	xmm4 *= xmm0;		/* new oddI2 = .951*(R3-R9) */ \
	xmm6 = CONST2_P588; \
	xmm0 *= xmm6;		/* new oddI3 = .588*(R3-R9) */ \
	xmm6 *= xmm1;		/* .588*(R5-R7) */ \
	xmm1 *= CONST2_P951;		/* .951*(R5-R7) */ \
	xmm4 += xmm6;		/* new oddI2 += .588*(R5-R7) */ \
	xmm0 -= xmm1;		/* new oddI3 -= .951*(R5-R7) */ \
\
	xmm7 = mem2;		/* R2 */ \
	xmm7 -= mem10;		/* R2-R10 */ \
	xmm1 = mem4;		/* R4 */ \
	xmm1 -= mem8;		/* R4-R8 */ \
	xmm6 = CONST2_P951; \
	xmm6 *= xmm7;		/* new evenI3 = .951*(R2-R10) */ \
	xmm2 = CONST2_P588; \
	xmm7 *= xmm2;		/* new evenI2 = .588*(R2-R10) */ \
	xmm2 *= xmm1;		/* .588*(R4-R8) */ \
	xmm1 *= CONST2_P951;		/* .951*(R4-R8) */ \
	xmm6 -= xmm2;		/* new evenI3 -= .588*(R4-R8) */ \
	xmm7 += xmm1;		/* new evenI2 += .951*(R4-R8) */ \
\
	xmm1 = xmm7;		/* Copy evenI2 */ \
	xmm7 -= xmm4;		/* New I5 = evenI2 - oddI2 */ \
	xmm4 += xmm1;		/* New I2 = evenI2 + oddI2 */ \
\
	xmm1 = xmm6;		/* Copy evenI3 */ \
	xmm6 -= xmm0;		/* New I4 = evenI3 - oddI3 */ \
	xmm0 += xmm1;		/* New I3 = evenI3 + oddI3 */ \
\
;	/* apply twiddles */ \
\
	xmm1 = xptr(screg1+48);	/* cosine/sine for w^4 */ \
	xmm2 = xmm3;		/* Copy R5 */ \
	xmm3 *= xmm1;		/* A5 = R5 * cosine/sine */ \
	xmm3 -= xmm7;		/* A5 = A5 - I5 */ \
	xmm7 *= xmm1;		/* B5 = I5 * cosine/sine */ \
	xmm7 += xmm2;		/* B5 = B5 + R5 */ \
\
	xmm1 = xptr(screg2+48);	/* cosine/sine for w^3 */ \
	xmm2 = xmm5;		/* Copy R4 */ \
	xmm5 *= xmm1;		/* A4 = R4 * cosine/sine */ \
	xmm5 -= xmm6;		/* A4 = A4 - I4 */ \
	xmm6 *= xmm1;		/* B4 = I4 * cosine/sine */ \
	xmm6 += xmm2;		/* B4 = B4 + R4 */ \
\
	xmm2 = xptr(screg1+32); \
	xmm3 *= xmm2;		/* A5 = A5 * sine (new R5) */ \
	xmm7 *= xmm2;		/* B5 = B5 * sine (new I5) */ \
	xmm2 = xptr(screg2+32); \
	xmm5 *= xmm2;		/* A4 = A4 * sine (new R4) */ \
	xmm6 *= xmm2;		/* B4 = B4 * sine (new I4) */ \
\
	xmm1 = xptr(screg1+16);	/* cosine/sine for w^2 */ \
	xmm1 *= TMP2;		/* A3 = R3 * cosine/sine */ \
	xmm1 -= xmm0;		/* A3 = A3 - I3 */ \
	xmm0 *= xptr(screg1+16);	/* B3 = I3 * cosine/sine */ \
	xmm0 += TMP2;		/* B3 = B3 + R3 */ \
\
	xmm2 = xptr(screg2+16);	/* cosine/sine for w^1 */ \
	xmm2 *= TMP1;		/* A2 = R2 * cosine/sine */ \
	xmm2 -= xmm4;		/* A2 = A2 - I2 */ \
	xmm4 *= xptr(screg2+16);	/* B2 = I2 * cosine/sine */ \
	xmm4 += TMP1;		/* B2 = B2 + R2 */ \
\
	xmm1 *= xptr(screg1);		/* A3 = A3 * sine (new R3) */ \
	xmm0 *= xptr(screg1);		/* B3 = B3 * sine (new I3) */ \
	xmm2 *= xptr(screg2);		/* A2 = A2 * sine (new R2) */ \
	xmm4 *= xptr(screg2);		/* B2 = B2 * sine (new I2) */ \

#endif

/* This is used in the first levels of pass 2 if pass 1 does the swizzling */ \
/* The ten-reals macro and the five-complex share an XMM register. */ \
/* This isn't very efficient, but this macro isn't called a whole lot. */ \

#define r5_fh5cl_ten_reals_five_complex_djbfft(srcreg,srcinc,d1,screg1,scoff1,screg2,scoff2) { \
	uintptr_t src_rbx = srcreg+rbx; \
	vec2f64 *sd0 = (vec2f64*)(srcreg); \
	vec2f64 *sd1 = (vec2f64*)(srcreg+(d1)); \
	vec2f64 *sd2 = (vec2f64*)(srcreg+2*(d1)); \
	vec2f64 *sd3 = (vec2f64*)(srcreg+3*(d1)); \
	vec2f64 *sd4 = (vec2f64*)(srcreg+4*(d1)); \
	vec2f64 *bd0 = (vec2f64*)(src_rbx); \
	vec2f64 *bd1 = (vec2f64*)(src_rbx+(d1)); \
	vec2f64 *bd2 = (vec2f64*)(src_rbx+2*(d1)); \
	vec2f64 *bd3 = (vec2f64*)(src_rbx+3*(d1)); \
	vec2f64 *bd4 = (vec2f64*)(src_rbx+4*(d1)); \
	/*const uintptr_t n0 = d1; \
	const uintptr_t n1 = d2+16; \
	const uintptr_t n2 = d4+16; \
	const uintptr_t n3 = d1+32; \
	const uintptr_t n4 = d3+32; \
	const uintptr_t n5 = d2+48; \
	const uintptr_t n6 = d4+48;*/ \
	const uintptr_t s1 = scoff1; \
	const uintptr_t s2 = scoff2; \
	r5_h10r_h5c_djbfft_mem(bd1[0],bd3[0],bd0[1],bd2[1],bd4[1], bd1[2],bd3[2],bd0[3],bd2[3],bd4[3], screg1+s1,screg1+s1+32,screg2+s2, sd0[1],sd0[3],sd4[1]); \
/*	xstore	[srcreg+16], xmm2		;; Save R1 #1/R1 */ \
/*	xstore	[srcreg+48], xmm5		;; Save R1 #2/I1 */ \
	sd2[1] = xmm3;		/* Save R3 */ \
	sd2[3] = xmm4;		/* Save I3 */ \
	vec2f64 TMP5 = bd1[1];	/* R4 */ \
	sd1[1] = xmm1;		/* Save R2 */ \
	vec2f64 TMP6 = bd1[3];	/* R9 */ \
	sd1[3] = xmm0;		/* Save I2 */ \
	vec2f64 TMP7 = bd3[1];	/* R5 */ \
	sd3[1] = xmm5;		/* Save R4 */ \
	vec2f64 TMP8 = bd3[3];	/* R10 */ \
	sd3[3] = xmm2;		/* Save I4 */ \
/*	xstore	[srcreg+4*d1+16], xmm7		;; Save R5 */ \
	sd4[3] = xmm6;		/* Save I5 */ \
	r5_h10r_h5c_djbfft_mem(bd0[0],bd2[0],bd4[0],TMP5,TMP7, bd0[2],sd2[2],sd4[2],TMP6,TMP8, screg1,screg1+32,screg2,sd0[0],sd0[2],sd4[0]); \
/*	xstore	[srcreg], xmm2			;; Save R1 #1/R1 */ \
/*	xstore	[srcreg+32], xmm5		;; Save R1 #2/I1 */ \
	sd1[0] = xmm1;		/* Save R2 */ \
	sd1[2] = xmm0;		/* Save I2 */ \
	sd2[0] = xmm3;		/* Save R3 */ \
	sd2[2] = xmm4;		/* Save I3 */ \
	sd3[0] = xmm5;		/* Save R4 */ \
	sd3[2] = xmm2;		/* Save I4 */ \
/*	xstore	[srcreg+4*d1], xmm7		;; Save R5 */ \
	sd4[2] = xmm6;		/* Save I5 */ \
	srcreg += srcinc; \
}

/* This is used in the later radix-5 levels of pass 2. */ \
/* The ten-reals macro and the five-complex share an XMM register. */ \
/* This isn't very efficient, but this macro isn't called a whole lot. */ \

#define r5_h5cl_ten_reals_five_complex_djbfft(srcreg,srcinc,d1,screg1,scoff1,screg2,scoff2) { \
	vec2f64 *sd0 = (vec2f64*)(srcreg); \
	vec2f64 *sd1 = (vec2f64*)(srcreg+(d1)); \
	vec2f64 *sd2 = (vec2f64*)(srcreg+2*(d1)); \
	vec2f64 *sd3 = (vec2f64*)(srcreg+3*(d1)); \
	vec2f64 *sd4 = (vec2f64*)(srcreg+4*(d1)); \
	vec2f64 TMP5, TMP6, TMP7, TMP8; \
	/*const uintptr_t n0 = d1; \
	const uintptr_t n1 = d2+16; \
	const uintptr_t n2 = d4+16; \
	const uintptr_t n3 = d1+32; \
	const uintptr_t n4 = d3+32; \
	const uintptr_t n5 = d2+48; \
	const uintptr_t n6 = d4+48;*/ \
	const uintptr_t s1 = scoff1; \
	const uintptr_t s2 = scoff2; \
	r5_h10r_h5c_djbfft_mem(sd1[0],sd3[0],sd0[1],sd2[1],sd4[1], sd1[2],sd3[2],sd0[3],sd2[3],sd4[3], screg1+s1,screg1+s1+32,screg2+s2, sd0[1],sd0[3],sd4[1]); \
/*	xstore	[srcreg+16], xmm2		;; Save R1 #1/R1 */ \
/*	xstore	[srcreg+48], xmm5		;; Save R1 #2/I1 */ \
/*	xstore	[srcreg+4*d1+16], xmm7		;; Save R5 */ \
	sd4[3] = xmm6;		/* Save I5 */ \
	sd2[1] = xmm3;		/* Save R3 */ \
	sd2[3] = xmm4;		/* Save I3 */ \
	TMP5 = sd1[1];		/* R4 */ \
	sd1[1] = xmm1;		/* Save R2 */ \
	TMP6 = sd1[3];		/* R9 */ \
	sd1[3] = xmm0;		/* Save I2 */ \
	TMP7 = sd3[1];		/* R5 */ \
	sd3[1] = xmm5;		/* Save R4 */ \
	TMP8 = sd3[3];		/* R10 */ \
	sd3[3] = xmm2;		/* Save I4 */ \
	r5_h10r_h5c_djbfft_mem(sd0[0],sd2[0],sd4[0],TMP5,TMP7, sd0[2],sd2[2],sd4[2],TMP6,TMP8, screg1,screg1+32,screg2, sd0[0],sd0[2],sd4[0]); \
/*	xstore	[srcreg], xmm2			;; Save R1 #1/R1 */ \
/*	xstore	[srcreg+32], xmm5		;; Save R1 #2/I1 */ \
	sd1[0] = xmm1;		/* Save R2 */ \
	sd1[2] = xmm0;		/* Save I2 */ \
	sd2[0] = xmm3;		/* Save R3 */ \
	sd2[2] = xmm4;		/* Save I3 */ \
	sd3[0] = xmm5;		/* Save R4 */ \
	sd3[2] = xmm2;		/* Save I4 */ \
/*	xstore	[srcreg+4*d1], xmm7		;; Save R5 */ \
	sd4[2] = xmm6;		/* Save I5 */ \
	srcreg += srcinc; \
}

/* This is used in the later radix-5 levels of pass 2.  Uses two sin/cos */ \
/* ptrs for case where first levels in pass 2 were radix-3. */ \
/* The ten-reals macro and the five-complex share an XMM register. */ \
/* This isn't very efficient, but this macro isn't called a whole lot. */ \

#define r5_h5cl_2sc_ten_reals_five_complex_djbfft(srcreg,srcinc,d1,screg1,s1,screg2,s2,screg3,s3) { \
	vec2f64 *sd0 = (vec2f64*)(srcreg); \
	vec2f64 *sd1 = (vec2f64*)(srcreg+(d1)); \
	vec2f64 *sd2 = (vec2f64*)(srcreg+2*(d1)); \
	vec2f64 *sd3 = (vec2f64*)(srcreg+3*(d1)); \
	vec2f64 *sd4 = (vec2f64*)(srcreg+4*(d1)); \
	vec2f64 TMP5, TMP6, TMP7, TMP8; \
	/*s1 = scoff1;*/ \
	/*s2 = scoff2;*/ \
	/*s3 = scoff3;*/ \
	r5_h10r_h5c_djbfft_mem(sd1[0],sd3[0],sd0[1],sd2[1],sd4[1], sd1[2],sd3[2],sd0[3],sd2[3],sd4[3], screg1+s1,screg2+s2,screg3+s3, sd0[1],sd0[3],sd4[1]); \
/*	xstore	[srcreg+16], xmm2		;; Save R1 #1/R1 */ \
/*	xstore	[srcreg+48], xmm5		;; Save R1 #2/I1 */ \
	sd2[1] = xmm3;		/* Save R3 */ \
	sd2[3] = xmm4;		/* Save I3 */ \
	TMP5 = sd1[1];		/* R4 */ \
	sd1[1] = xmm1;		/* Save R2 */ \
	TMP6 = sd1[3];		/* R9 */ \
	sd1[3] = xmm0;		/* Save I2 */ \
	TMP7 = sd3[1];		/* R5 */ \
	sd3[1] = xmm5;		/* Save R4 */ \
	TMP8 = sd3[3];		/* R10 */ \
	sd3[3] = xmm2;		/* Save I4 */ \
/*	xstore	[srcreg+4*d1+16], xmm7		;; Save R5 */ \
	sd4[3] = xmm6;		/* Save I5 */ \
	r5_h10r_h5c_djbfft_mem(sd0[0],sd2[0],sd4[0],TMP5,TMP7, sd0[2],sd2[2],sd4[2],TMP6,TMP8, screg1,screg2,screg3, sd0[0],sd0[2],sd4[0]); \
/*	xstore	[srcreg], xmm2			;; Save R1 #1/R1 */ \
/*	xstore	[srcreg+32], xmm5		;; Save R1 #2/I1 */ \
	sd1[0] = xmm1;		/* Save R2 */ \
	sd1[2] = xmm0;		/* Save I2 */ \
	sd2[0] = xmm3;		/* Save R3 */ \
	sd2[2] = xmm4;		/* Save I3 */ \
	sd3[0] = xmm5;		/* Save R4 */ \
	sd3[2] = xmm2;		/* Save I4 */ \
/*	xstore	[srcreg+4*d1], xmm7		;; Save R5 */ \
	sd4[2] = xmm6;		/* Save I5 */ \
	srcreg += srcinc; \
}


#define r5_h10r_h5c_djbfft_mem(mem1p,mem2p,mem3p,mem4p,mem5p,mem6p,mem7p,mem8p,mem9p,mem10p, screg1,screg2,screg3,dst1,dst2,dstr5) { \
	double TMP1, TMP2, TMP3, TMP4;\
	double f0,f1,f2,f3,f4,f5,f6,f7; \
	vec2f64 mem2 = mem2p; \
	vec2f64 mem5 = mem5p; \
	vec2f64 mem3 = mem3p; \
	vec2f64 mem4 = mem4p; \
	vec2f64 mem7 = mem7p; \
	vec2f64 mem10 = mem10p; \
	vec2f64 mem8 = mem8p; \
	vec2f64 mem9 = mem9p; \
	vec2f64 mem1 = mem1p; \
	vec2f64 mem6 = mem6p; \
	/* Do the five complex part */ \
	f0 = mem2[1]; \
	f0 += mem5[1];			/* r25a=r2+r5 */ \
	f1 = mem3[1]; \
	f1 += mem4[1];			/* r34a=r3+r4 */ \
	f2 = mem7[1]; \
	f2 -= mem10[1];			/* i25s=i2-i5 */ \
	f3 = mem8[1]; \
	f3 -= mem9[1];			/* i34s=i3-i4 */ \
	f4 = f0 * CONST_P309;	/* cos2*r25a */ \
	f4 += f1 * CONST_M809;	/* cos2*r25a + cos4*r34a */ \
	f6 = f2 * CONST_P951;	/* sin2*i25s */ \
	f6 += f3 * CONST_P588;	/* t2=sin2*i25s + sin4*i34s */ \
	f5 = f0 * CONST_M809;	/* cos4*r25a */ \
	f5 += f1 * CONST_P309;	/* cos4*r25a+cos2*r34a */ \
	f2 = f2 * CONST_P588;	/* sin4*i25s */ \
	f2 -= f3 * CONST_P951;	/* t4=sin4*i25s-sin2*i34s */ \
	f0 = f0 + f1;			/* r25a + r34a */ \
	f7 = mem1[1]; \
	f4 += f7;				/* t1=cos2*r25a + cos4*r34a + r1 */ \
	f0 += f7;				/* outr(0) = r1 + r25a + r34a */ \
	f5 += f7;				/* t3=cos4*r25a+cos2*r34a+r1 */ \
	f7 = f4; \
	f4 -= f6;				/* outr(1)=t1-t2 */ \
	f6 += f7;				/* outr(4)=t1+t2 */ \
	f7 = f5; \
	f5 -= f2;				/* outr(2)=t3-t4 */ \
	f2 += f7;				/* outr(3)=t3+t4 */ \
\
	TMP1 = f4;				/* Save new r2 */ \
	TMP2 = f5;				/* Save new r3 */ \
	TMP3 = f2;				/* Save new r4 */ \
	TMP4 = f6;				/* Save new r5 */ \
\
	f2 = mem7[1]; \
	f2 += mem10[1];			/* i25a=i2+i5 */ \
	f3 = mem8[1]; \
	f3 += mem9[1];			/* i34a=i3+i4 */ \
	f4 = mem2[1]; \
	f4 -= mem5[1];			/* r25s=r2-r5 */ \
	f1 = mem3[1]; \
	f1 -= mem4[1];			/* r34s=r3-r4 */ \
\
	dst1[1] = f0;			/* WARNING: dst1 may be alias to mem3 */ \
\
	f5 = f2; \
	f5 += f3;				/* i25a+i34a */ \
	f6 = f2 * CONST_P309;	/* cos2*i25a */ \
	f6 += f3 * CONST_M809;	/* cos2*i25a + cos4*i34a */ \
	f2 = f2 * CONST_M809;	/* cos4*i25a */ \
	f2 += f3 * CONST_P309;	/* cos4*i25a + cos2*i34a */ \
	f0 = f4 * CONST_P951;	/* sin2*r25s */ \
	f0 += f1 * CONST_P588;	/* t6=sin2*r25s + sin4*r34s */ \
	f4 = f4 * CONST_P588;	/* sin4*r25s */ \
	f4 -= f1 * CONST_P951;	/* t8=sin4*r25s - sin2*r34s */ \
	f3 = mem6[1]; \
	f6 += f3;				/* t5=cos2*i25a + cos4*i34a + i1 */ \
	f2 += f3;				/* t7=cos4*i25a + cos2*i34a + i1 */ \
\
	f7 = f6; \
	f6 -= f0;				/* outi(4)=t5-t6 */ \
	f0 += f7;				/* outi(1)=t5+t6 */ \
	f1 = f2; \
	f2 -= f4;				/* outi(3)=t7-t8 */ \
	f4 += f1;				/* outi(2)=t7+t8 */ \
	f5 += f3;				/* outi(0)=i1+i25a+i34a */ \
\
	f3 = f64ptr(screg1+16);	/* Cosine/sine for w^n */ \
	f7 = TMP4 * f3 + f6;	/* A5 = R5 * cosine/sine + I5 */ \
	f1 = TMP1 * f3 - f0;	/* A2 = R2 * cosine/sine - I2 */ \
	f6 = f6 * f3 - TMP4;	/* B5 = I5 * cosine/sine - R5 */ \
	f0 = f0 * f3 + TMP1;	/* B2 = I2 * cosine/sine + R2 */ \
	dst2[1] = f5; \
	f5 = f64ptr(screg1); \
	f7 *= f5;				/* A5 = A5 * sine (new R5) */ \
	dstr5[1] = f7; \
	f1 *= f5;				/* A2 = A2 * sine (new R2) */ \
	f6 *= f5;				/* B5 = B5 * sine (new I5) */ \
	f0 *= f5;				/* B2 = B2 * sine (new I2) */ \
\
	f7 = f64ptr(screg2+16);		/* Cosine/sine for w^2n */ \
	f3 = TMP2 * f7 - f4;	/* A3 = R3 * cosine/sine - I3 */ \
	f5 = TMP3 * f7 + f2;	/* A4 = R4 * cosine/sine + I4 */ \
	f4 = f4 * f7 + TMP2;	/* B3 = I3 * cosine/sine + R3 */ \
	f2 = f2 * f7 - TMP3;	/* B4 = I4 * cosine/sine - R4 */ \
	f7 = f64ptr(screg2); \
	f3 *= f7;				/* A3 = A3 * sine (new R3) */ \
	f5 *= f7;				/* A4 = A4 * sine (new R4) */ \
	f4 *= f7;				/* B3 = B3 * sine (new I3) */ \
	f2 *= f7;				/* B4 = B4 * sine (new I4) */ \
\
;	/* Do the ten-reals part */ \
	xmm1[1] = f1;			/* Copy to high part of XMM register */ \
	xmm6[1] = f6; \
	xmm0[1] = f0; \
	xmm3[1] = f3; \
	xmm5[1] = f5; \
	xmm4[1] = f4; \
	xmm2[1] = f2; \
\
	f3 = mem3[0];			/* R3 */ \
	f3 += mem9[0];			/* R3+R9 */ \
	f2 = mem5[0];			/* R5 */ \
	f2 += mem7[0];			/* R5+R7 */ \
	f4 = mem1[0];			/* R1 */ \
	f0 = f3 * CONST_P309 + f4;	/* new oddR2 = .309*(R3+R9) + R1 */ \
	f5 = f3 * CONST_M809 + f4;	/* new oddR3 = -.809*(R3+R9) + R1 */ \
	f3 += f4;				/* new oddR1 = R1+R3+R9 */ \
	f0 += f2 * CONST_M809;	/* new oddR2 += -.809*(R5+R7) */ \
	f3 += f2;				/* new oddR1 = R1+R3+R5+R7+R9 */ \
	f5 += f2 * CONST_P309;	/* new oddR3 += .309*(R5+R7) */ \
\
	f6 = mem2[0];			/* R2 */ \
	f6 += mem10[0];			/* R2+R10 */ \
	f2 = mem4[0];			/* R4 */ \
	f2 += mem8[0];			/* R4+R8 */ \
	f4 = mem6[0];			/* R6 */ \
	f1 = f6 * CONST_P309 + f4;	/* new evenR3 = .309*(R2+R10) + R6 */ \
	f7 = f6 * CONST_M809 + f4;	/* new evenR5 = -.809*(R2+R10) + R6 */ \
	f6 += f4;				/* new evenR1 = R6+R2+R10 */ \
	f1 += f2 * CONST_M809;	/* new evenR3 += -.809*(R4+R8) */ \
	f6 += f2;				/* new evenR1 = R6+R2+R4+R8+R10 */ \
	f7 += f2 *= CONST_P309;	/* new evenR5 += .309*(R4+R8) */ \
\
	f2 = f0;				/* Copy oddR2 */ \
	f0 -= f7;				/* New R2 = oddR2 - evenR5 */ \
	f7 += f2;				/* New R5 = oddR2 + evenR5 */ \
\
	f2 = f5;				/* Copy oddR3 */ \
	f5 -= f1;				/* New R4 = oddR3 - evenR3 */ \
	f1 += f2;				/* New R3 = oddR3 + evenR3 */ \
\
	TMP1 = f0;				/* Temporarily save new R2 */ \
	TMP2 = f1;				/* Temporarily save new R3 */ \
\
	f4 = mem3[0];			/* R3 */ \
	f4 -= mem9[0];			/* R3-R9 */ \
	dst1[0] = f3;			/* Save oddR1, WARNING: dst1 may be alias to mem1 or mem3 */ \
	f3 = mem5[0];			/* R5 */ \
	f3 -= mem7[0];			/* R5-R7 */ \
	f0 = f4 * CONST_P951;	/* new oddI2 = .951*(R3-R9) */ \
	f4 *= CONST_P588;		/* new oddI3 = .588*(R3-R9) */ \
	f0 += f3 * CONST_P588;	/* new oddI2 += .588*(R5-R7) */ \
	f4 -= f3 * CONST_P951;	/* new oddI3 -= .951*(R5-R7) */ \
\
	f3 = mem4[0];			/* R4 */ \
	f3 -= mem8[0];			/* R4-R8 */ \
	dst2[0] = f6;			/* Save evenR1, WARNING: dst2 may be alias to mem6 or mem8 */ \
	f6 = mem2[0];			/* R2 */ \
	f6 -= mem10[0];			/* R2-R10 */ \
	f2 = f6 * CONST_P951;	/* new evenI3 = .951*(R2-R10) */ \
	f6 *= CONST_P588;		/* new evenI2 = .588*(R2-R10) */ \
	f2 -= f3 * CONST_P588;	/* new evenI3 -= .588*(R4-R8) */ \
	f6 += f3 * CONST_P951;	/* new evenI2 += .951*(R4-R8) */ \
\
	f3 = f6;				/* Copy evenI2 */ \
	f6 -= f0;				/* New I5 = evenI2 - oddI2 */ \
	f0 += f3;				/* New I2 = evenI2 + oddI2 */ \
\
	f3 = f2;				/* Copy evenI3 */ \
	f2 -= f4;				/* New I4 = evenI3 - oddI3 */ \
	f4 += f3;				/* New I3 = evenI3 + oddI3 */ \
\
;	/* apply twiddles */ \
\
	f3 = f64ptr(screg2+16);	/* cosine/sine for w^4 */ \
	f1 = f7;				/* Copy R5 */ \
	f7 = f7 * f3 - f6;		/* A5 = R5 * cosine/sine - I5 */ \
	f6 = f6 * f3 + f1;		/* B5 = I5 * cosine/sine + R5 */ \
\
	f3 = f64ptr(screg3+24);	/* cosine/sine for w^3 */ \
	f1 = f5;				/* Copy R4 */ \
	f5 = f5 * f3 - f2;		/* A4 = R4 * cosine/sine - I4 */ \
	f2 = f2 * f3 + f1;		/* B4 = I4 * cosine/sine + R4 */ \
\
	f1 = f64ptr(screg2); \
	f7 *= f1;				/* A5 = A5 * sine (new R5) */ \
	f6 *= f1;				/* B5 = B5 * sine (new I5) */ \
	f1 = f64ptr(screg3+16); \
	f5 *= f1;				/* A4 = A4 * sine (new R4) */ \
	f2 *= f1;				/* B4 = B4 * sine (new I4) */ \
\
	dstr5[0] = f7;			/* Save R5 */ \
\
	f1 = f64ptr(screg1+16);	/* cosine/sine for w^2 */ \
	f3 = TMP2 * f1 - f4;	/* A3 = R3 * cosine/sine - I3 */ \
	f4 = f4 * f1 + TMP2;	/* B3 = I3 * cosine/sine + R3 */ \
\
	f7 = f64ptr(screg3+8);	/* cosine/sine for w^1 */ \
	f1 = f7 * TMP1 - f0;	/* A2 = R2 * cosine/sine - I2 */ \
	f0 = f0 * f7 + TMP1;	/* B2 = I2 * cosine/sine + R2 */ \
\
	f3 *= f64ptr(screg1);	/* A3 = A3 * sine (new R3) */ \
	f4 *= f64ptr(screg1);	/* B3 = B3 * sine (new I3) */ \
	f1 *= f64ptr(screg3);	/* A2 = A2 * sine (new R2) */ \
	f0 *= f64ptr(screg3);	/* B2 = B2 * sine (new I2) */ \
	xmm1[0] = f1; \
	xmm6[0] = f6; \
	xmm0[0] = f0; \
	xmm3[0] = f3; \
	xmm5[0] = f5; \
	xmm4[0] = f4; \
	xmm2[0] = f2; \
}


/* */ \
/* */ \
/* ************************************* ten-reals-unfft variants ****************************************** */ \
/* */ \


#if 0 // unused
#define r5_x5cl_ten_reals_unfft(srcreg,srcinc,d1,screg1,scoff1,screg2,scoff2) { \
	vec2f64 TMP1, TMP2, TMP3, TMP4; \
	d2 = 2*d1; \
	d3 = 3*d1; \
	d4 = 4*d1; \
	s1 = scoff1; \
	s2 = scoff2; \
	n1 = d1; \
	n2 = d1+32; \
	n4 = d2+32; \
	n6 = d3+32; \
	n8 = d4+32; \
	x10u(xptr(srcreg),xptr(srcreg+32),xptr(srcreg+n1),xptr(srcreg+n2),xptr(srcreg+d2),xptr(srcreg+n4),xptr(srcreg+d3),xptr(srcreg+n6),xptr(srcreg+d4),xptr(srcreg+n8),screg1,screg2,xptr(srcreg),xptr(srcreg+d4),xptr(srcreg+32),xptr(srcreg+n4)); \
	xmm1 = xptr(srcreg+d1+16); \
	TMP1 = xmm1; \
	xmm1 = xptr(srcreg+3*d1+16); \
	TMP2 = xmm1; \
	xmm1 = xptr(srcreg+d1+48); \
	TMP3 = xmm1; \
	xmm1 = xptr(srcreg+3*d1+48); \
	TMP4 = xmm1; \
/*	xstore	[srcreg], xmm4		;; Save R1 */ \
	xptr(srcreg+2*d1) = xmm0;	/* Save R2 */ \
/*	xstore	[srcreg+4*d1], xmm5	;; Save R3 */ \
	xptr(srcreg+d1+16) = xmm6;	/* Save R4 */ \
	xptr(srcreg+3*d1+16) = xmm7;	/* Save R5 */ \
/*	xstore	[srcreg+32], xmm1	;; Save R6 */ \
/*	xstore	[srcreg+2*d1+32], xmm0	;; Save R7 */ \
	xptr(srcreg+4*d1+32) = xmm3;	/* Save R8 */ \
	xptr(srcreg+d1+48) = xmm4;	/* Save R9 */ \
	xptr(srcreg+3*d1+48) = xmm2;	/* Save R10 */ \
	n1 = d1; \
	n2 = d2+16; \
	n3 = d2+48; \
	n4 = d4+16; \
	n5 = d4+48; \
	n6 = d1+32; \
	n7 = d3+32; \
	x10u(xptr(srcreg+16),xptr(srcreg+48),TMP1,TMP3,xptr(srcreg+n2),xptr(srcreg+n3),TMP2,TMP4,xptr(srcreg+n4),xptr(srcreg+n5),screg1+s1,screg2+s2,xptr(srcreg+n1),xptr(srcreg+16),xptr(srcreg+n6),xptr(srcreg+n7)); \
/*	xstore	[srcreg+d1], xmm4	;; Save R1 */ \
	xptr(srcreg+3*d1) = xmm0;	/* Save R2 */ \
/*	xstore	[srcreg+16], xmm5	;; Save R3 */ \
	xptr(srcreg+2*d1+16) = xmm6;	/* Save R4 */ \
	xptr(srcreg+4*d1+16) = xmm7;	/* Save R5 */ \
/*	xstore	[srcreg+d1+32], xmm3	;; Save R6 */ \
/*	xstore	[srcreg+3*d1+32], xmm3	;; Save R7 */ \
	xptr(srcreg+48) = xmm3;	/* Save R8 */ \
	xptr(srcreg+2*d1+48) = xmm4;	/* Save R9 */ \
	xptr(srcreg+4*d1+48) = xmm2;	/* Save R10 */ \
	srcreg += srcinc; \
}


/* To calculate a 10-reals unFFT (in a shorthand notation): */ \
/* c1 + c2 + c3 + c4 + c5 + c6 + c7 + c8 + c9 + c10	*  w^0000000000 */ \
/* c1 + c2 + c3 + c4 + c5 + c6 + c7 + c8 + c9 + c10	*  w^0123456789 */ \
/* c1 + c2 + c3 + c4 + c5 + c6 + c7 + c8 + c9 + c10	*  w^0246802468 */ \
/* c1 + c2 + c3 + c4 + c5 + c6 + c7 + c8 + c9 + c10	*  w^0369258147 */ \
/* c1 + c2 + c3 + c4 + c5 + c6 + c7 + c8 + c9 + c10	*  w^0482604826 */ \
/* c1 + c2 + c3 + c4 + c5 + c6 + c7 + c8 + c9 + c10	*  w^0505050505 */ \
/* c1 + c2 + c3 + c4 + c5 + c6 + c7 + c8 + c9 + c10	*  w^0628406284 */ \
/* c1 + c2 + c3 + c4 + c5 + c6 + c7 + c8 + c9 + c10	*  w^0741852963 */ \
/* c1 + c2 + c3 + c4 + c5 + c6 + c7 + c8 + c9 + c10	*  w^0864208642 */ \
/* c1 + c2 + c3 + c4 + c5 + c6 + c7 + c8 + c9 + c10	*  w^0987654321 */ \
/* Noting that w^5 = -1 */ \
/* c1 + c2 + c3 + c4 + c5 + c6 + c7 + c8 + c9 + c10	*  w^0000000000 */ \
/* c1 + c2 + c3 + c4 + c5 - c6 - c7 - c8 - c9 - c10	*  w^0123401234 */ \
/* c1 + c2 + c3 - c4 - c5 + c6 + c7 + c8 - c9 - c10	*  w^0241302413 */ \
/* c1 + c2 - c3 - c4 + c5 - c6 - c7 + c8 + c9 - c10	*  w^0314203142 */ \
/* c1 + c2 - c3 + c4 - c5 + c6 + c7 - c8 + c9 - c10	*  w^0432104321 */ \
/* c1 - c2 + c3 - c4 + c5 - c6 + c7 - c8 + c9 - c10	*  w^0000000000 */ \
/* c1 - c2 + c3 - c4 + c5 + c6 - c7 + c8 - c9 + c10	*  w^0123401234 */ \
/* c1 - c2 + c3 + c4 - c5 - c6 + c7 - c8 - c9 + c10	*  w^0241302413 */ \
/* c1 - c2 - c3 + c4 + c5 + c6 - c7 - c8 + c9 + c10	*  w^0314203142 */ \
/* c1 - c2 - c3 - c4 - c5 - c6 + c7 + c8 + c9 + c10	*  w^0432104321 */ \
/* incoming is:	r1_1 = c1r + c6r */ \
/*		c2 = c2r + c2i */ \
/*		c3 = c3r + c3i */ \
/*		c4 = c4r + c4i */ \
/*		c5 = c5r + c5i */ \
/*		r1_2 = c1r - c6r */ \
/*		c7 = c5r - c5i	(implied) */ \
/*		c8 = c4r - c4i	(implied) */ \
/*		c9 = c3r - c3i	(implied) */ \
/*		c10 = c2r - c2i	(implied) */ \
/* And noticing the signs of the real and imaginary parts of the sin/cos values: */ \
/* w^1/10 = .809 - .588i */ \
/* w^2/10 = .309 - .951i */ \
/* w^3/10 = -.309 - .951i */ \
/* w^4/10 = -.809 - .588i */ \
/* We get reals: */ \
/* r1_1 + 2c2 + 2c3 + 2c4 + 2c5	*  w^00000 */ \
/* r1_2 + 2c2 + 2c3 + 2c4 + 2c5	*  w^01234 */ \
/* r1_1 + 2c2 + 2c3 - 2c4 - 2c5	*  w^02413 */ \
/* r1_2 + 2c2 - 2c3 - 2c4 + 2c5	*  w^03142 */ \
/* r1_1 + 2c2 - 2c3 + 2c4 - 2c5	*  w^04321 */ \
/* r1_2 - 2c2 + 2c3 - 2c4 + 2c5	*  w^00000 */ \
/* r1_1 - 2c2 + 2c3 - 2c4 + 2c5	*  w^01234 */ \
/* r1_2 - 2c2 + 2c3 + 2c4 - 2c5	*  w^02413 */ \
/* r1_1 - 2c2 - 2c3 + 2c4 + 2c5	*  w^03142 */ \
/* r1_2 - 2c2 - 2c3 - 2c4 - 2c5	*  w^04321 */ \
/* Now drop the multiplication by 2 (the actual r1_1 and r1_2 inputs are already doubled) */ \
/* and expand the sin/cos multipliers: */ \
/* r1_1 + c2r + c3r + c4r + c5r */ \
/* r1_2 + .809c2r + .588c2i + .309c3r + .951c3i - .309c4r + .951c4i - .809c5r + .588c5i */ \
/* r1_1 + .309c2r + .951c2i - .809c3r + .588c3i - .809c4r - .588c4i + .309c5r - .951c5i */ \
/* r1_2 - .309c2r + .951c2i - .809c3r - .588c3i + .809c4r - .588c4i + .309c5r + .951c5i */ \
/* r1_1 - .809c2r + .588c2i + .309c3r - .951c3i + .309c4r + .951c4i - .809c5r - .588c5i */ \
/* r1_2 - c2r + c3r - c4r + c5r */ \
/* r1_1 - .809c2r - .588c2i + .309c3r + .951c3i + .309c4r - .951c4i - .809c5r + .588c5i */ \
/* r1_2 - .309c2r - .951c2i - .809c3r + .588c3i + .809c4r + .588c4i + .309c5r - .951c5i */ \
/* r1_1 + .309c2r - .951c2i - .809c3r - .588c3i - .809c4r + .588c4i + .309c5r + .951c5i */ \
/* r1_2 + .809c2r - .588c2i + .309c3r - .951c3i - .309c4r - .951c4i - .809c5r - .588c5i */ \
/* Simplify: */ \
/* r1_1 + c2r + c3r + c4r + c5r */ \
/* r1_2 + .809(c2r-c5r) + .588(c2i+c5i) + .309(c3r-c4r) + .951(c3i+c4i) */ \
/* r1_1 + .309(c2r+c5r) + .951(c2i-c5i) - .809(c3r+c4r) + .588(c3i-c4i) */ \
/* r1_2 - .309(c2r-c5r) + .951(c2i+c5i) - .809(c3r-c4r) - .588(c3i+c4i) */ \
/* r1_1 - .809(c2r+c5r) + .588(c2i-c5i) + .309(c3r+c4r) - .951(c3i-c4i) */ \
/* r1_2 - c2r + c3r - c4r + c5r */ \
/* r1_1 - .809(c2r+c5r) - .588(c2i-c5i) + .309(c3r+c4r) + .951(c3i-c4i) */ \
/* r1_2 - .309(c2r-c5r) - .951(c2i+c5i) - .809(c3r-c4r) + .588(c3i+c4i) */ \
/* r1_1 + .309(c2r+c5r) - .951(c2i-c5i) - .809(c3r+c4r) - .588(c3i-c4i) */ \
/* r1_2 + .809(c2r-c5r) - .588(c2i+c5i) + .309(c3r-c4r) - .951(c3i+c4i) */ \
\
/*r5_x10r_unfft_mem MACRO memr1_1,memr1_2,memr2,memi2,memr3,memi3,memr4,memi4,memr5,memi5,screg1,screg2,dst1,dst3,dst6,dst7 */ \

#define x10u(memr1_1,memr1_2,memr2,memi2,memr3,memi3,memr4,memi4,memr5,memi5,screg1,screg2,dst1,dst3,dst6,dst7) \
	xmm2 = memr2;		/* R2 */ \
	xmm2 *= xptr(screg2+16);	/* A2 = R2 * cosine/sine */ \
	xmm4 = memr3;		/* R3 */ \
	xmm4 *= xptr(screg1+16);	/* A3 = R3 * cosine/sine */ \
	xmm3 = memi2;		/* I2 */ \
	xmm2 += xmm3;		/* A2 = A2 + I2 */ \
	xmm5 = memi3;		/* I3 */ \
	xmm4 += xmm5;		/* A3 = A3 + I3 */ \
	xmm3 *= xptr(screg2+16);	/* B2 = I2 * cosine/sine */ \
	xmm5 *= xptr(screg1+16);	/* B3 = I3 * cosine/sine */ \
	xmm3 -= memr2;		/* B2 = B2 - R2 */ \
	xmm5 -= memr3;		/* B3 = B3 - R3 */ \
	xmm2 *= xptr(screg2);		/* A2 = A2 * sine (new R2) */ \
	xmm4 *= xptr(screg1);		/* A3 = A3 * sine (new R3) */ \
	xmm3 *= xptr(screg2);		/* B2 = B2 * sine (new I2) */ \
	xmm5 *= xptr(screg1);		/* B3 = B3 * sine (new I3) */ \
	memr2 = xmm2; \
	memi2 = xmm3; \
	memr3 = xmm4; \
	memi3 = xmm5; \
\
	xmm2 = memr4;		/* R4 */ \
	xmm2 *= xptr(screg2+48);	/* A4 = R4 * cosine/sine */ \
	xmm4 = memr5;		/* R5 */ \
	xmm4 *= xptr(screg1+48);	/* A5 = R5 * cosine/sine */ \
	xmm3 = memi4;		/* I4 */ \
	xmm2 += xmm3;		/* A4 = A4 + I4 */ \
	xmm5 = memi5;		/* I5 */ \
	xmm4 += xmm5;		/* A5 = A5 + I5 */ \
	xmm3 *= xptr(screg2+48);	/* B4 = I4 * cosine/sine */ \
	xmm5 *= xptr(screg1+48);	/* B5 = I5 * cosine/sine */ \
	xmm3 -= memr4;		/* B4 = B4 - R4 */ \
	xmm5 -= memr5;		/* B5 = B5 - R5 */ \
	xmm2 *= xptr(screg2+32);	/* A4 = A4 * sine (new R4) */ \
	xmm4 *= xptr(screg1+32);	/* A5 = A5 * sine (new R5) */ \
	xmm3 *= xptr(screg2+32);	/* B4 = B4 * sine (new I4) */ \
	xmm5 *= xptr(screg1+32);	/* B5 = B5 * sine (new I5) */ \
	memr4 = xmm2; \
	memi4 = xmm3; \
	memr5 = xmm4; \
	memi5 = xmm5; \
\
	xmm0 = memr2;		/* R2 */ \
	xmm0 += memr5;		/* R2+R5 */ \
	xmm1 = memr3;		/* R3 */ \
	xmm1 += memr4;		/* R3+R4 */ \
	xmm4 = CONST2_P309; \
	xmm4 *= xmm0;		/* new R3 = .309*(R2+R5) */ \
	xmm7 = memr1_1;		/* R1_1 */ \
	xmm4 += xmm7;		/* new R3 += R1_1 */ \
	xmm5 = CONST2_M809; \
	xmm5 *= xmm0;		/* new R5 = -.809*(R2+R5) */ \
	xmm0 += xmm7;		/* new R1 = R1_1+R2+R5 */ \
	xmm5 += xmm7;		/* new R5 += R1_1 */ \
	xmm7 = CONST2_M809; \
	xmm7 *= xmm1;		/* -.809*(R3+R4) */ \
	xmm0 += xmm1;		/* final R1 = R1_1+R2+R5+R3+R4 */ \
	xmm1 *= CONST2_P309;		/* .309*(R3+R4) */ \
	xmm4 += xmm7;		/* new R3 += -.809*(R3+R4) */ \
	xmm5 += xmm1;		/* new R5 += .309*(R3+R4) */ \
	dst1 = xmm0;		/* Save final R1 */ \
\
	xmm7 = memi2;		/* I2 */ \
	xmm7 -= memi5;		/* I2-I5 */ \
	xmm1 = memi3;		/* I3 */ \
	xmm1 -= memi4;		/* I3-I4 */ \
	xmm6 = CONST2_P951; \
	xmm6 *= xmm7;		/* tmp1 = .951*(I2-I5) */ \
	xmm2 = CONST2_P588; \
	xmm7 *= xmm2;		/* tmp2 = .588*(I2-I5) */ \
	xmm2 *= xmm1;		/* .588*(I3-I4) */ \
	xmm1 *= CONST2_P951;		/* .951*(I3-I4) */ \
	xmm6 += xmm2;		/* tmp1 += .588*(I3-I4) */ \
	xmm7 -= xmm1;		/* tmp2 -= .951*(I3-I4) */ \
\
	xmm1 = xmm4;		/* Copy new R3 */ \
	xmm4 -= xmm6;		/* final R9 = new R3 - tmp1 */ \
	xmm6 += xmm1;		/* final R3 = new R3 + tmp1 */ \
\
	xmm1 = xmm5;		/* Copy new R5 */ \
	xmm5 -= xmm7;		/* final R7 = new R5 - tmp2 */ \
	xmm7 += xmm1;		/* final R5 = new R5 + tmp2 */ \
\
	xmm0 = memr3;		/* R3 */ \
	xmm0 -= memr4;		/* R3-R4 */ \
	xmm1 = memr2;		/* R2 */ \
	xmm1 -= memr5;		/* R2-R5 */ \
	xmm2 = CONST2_P309; \
	xmm2 *= xmm0;		/* new R2 = .309*(R3-R4) */ \
	dst3 = xmm6;		/* Save final R3 */ \
	xmm6 = memr1_2;		/* R1_2 */ \
	xmm2 += xmm6;		/* new R2 += R1_2 */ \
	xmm3 = CONST2_M809; \
	xmm3 *= xmm0;		/* new R4 = -.809*(R3-R4) */ \
	xmm0 += xmm6;		/* new R6 = R1_2+(R3-R4) */ \
	xmm3 += xmm6;		/* new R4 += R1_2 */ \
	xmm6 = CONST2_M809; \
	xmm6 *= xmm1;		/* -.809*(R2-R5) */ \
	xmm0 -= xmm1;		/* new R6 = R6+(R3-R4)-(R2-R5) */ \
	xmm1 *= CONST2_P309;		/* .309*(R2-R5) */ \
	xmm2 -= xmm6;		/* new R2 -= -.809*(R2-R5) */ \
	xmm3 -= xmm1;		/* new R4 -= .309*(R2-R5) */ \
	dst6 = xmm0;		/* Save final R6 */ \
\
	xmm0 = memi2;		/* I2 */ \
	xmm0 += memi5;		/* I2+I5 */ \
	xmm1 = memi3;		/* I3 */ \
	xmm1 += memi4;		/* I3+I4 */ \
	xmm6 = CONST2_P951; \
	xmm6 *= xmm0;		/* tmp2 = .951*(I2+I5) */ \
	dst7 = xmm5;		/* Save final R7 */ \
	xmm5 = CONST2_P588; \
	xmm0 *= xmm5;		/* tmp1 = .588*(I2+I5) */ \
	xmm5 *= xmm1;		/* .588*(I3+I4) */ \
	xmm1 *= CONST2_P951;		/* .951*(I3+I4) */ \
	xmm6 -= xmm5;		/* tmp2 -= .588*(I3+I4) */ \
	xmm0 += xmm1;		/* tmp1 += .951*(I3+I4) */ \
\
	xmm1 = xmm2;		/* Copy new R2 */ \
	xmm2 -= xmm0;		/* Final R10 = new R2 - tmp1 */ \
	xmm0 += xmm1;		/* Final R2 = new R2 + tmp1 */ \
\
	xmm1 = xmm3;		/* Copy new R4 */ \
	xmm3 -= xmm6;		/* Final R8 = new R4 - tmp2 */ \
	xmm6 += xmm1;		/* Final R4 = new R4 + tmp2 */ \

#endif
/* Macro to do an ten_reals_unfft and a five_complex_djbunfft in pass 2. */ \
/* The ten-reals operation is done in the lower half of the XMM */ \
/* register.  This isn't very efficient, but this macro isn't called a whole lot. */ \
\

#define r5_h5cl_ten_reals_five_complex_djbunfft(srcreg,srcinc,d1,screg1,scoff1,screg2,scoff2) { \
	vec2f64 *sd0 = (vec2f64*)(srcreg); \
	vec2f64 *sd1 = (vec2f64*)(srcreg+(d1)); \
	vec2f64 *sd2 = (vec2f64*)(srcreg+2*(d1)); \
	vec2f64 *sd3 = (vec2f64*)(srcreg+3*(d1)); \
	vec2f64 *sd4 = (vec2f64*)(srcreg+4*(d1)); \
	xmm6 = sd1[0];	/* Load R2 */ \
	xmm7 = sd1[2];	/* Load I2 */ \
	vec2f64 TMP5 = xmm6; \
	vec2f64 TMP6 = xmm7; \
	xmm6 = sd3[0];	/* Load R4 */ \
	xmm7 = sd3[2];	/* Load I4 */ \
	vec2f64 TMP7 = xmm6; \
	vec2f64 TMP8 = xmm7; \
	/*const uintptr_t n1 = d1+16; \
	const uintptr_t n2 = d1+48; \
	const uintptr_t n3 = d2+16; \
	const uintptr_t n4 = d2+48; \
	const uintptr_t n5 = d3+16; \
	const uintptr_t n6 = d3+48; \
	const uintptr_t n7 = d4+16; \
	const uintptr_t n8 = d4+48; \
	const uintptr_t n0 = d1; \
	const uintptr_t n9 = d3; \
	const uintptr_t n10 = d1+32;*/ \
	r5_h10r_h5c_djbunfft_mem(sd0[1],sd0[3],sd1[1],sd1[3],sd2[1], sd2[3],sd3[1],sd3[3],sd4[1],sd4[3], screg1+scoff1,screg1+scoff1+32,screg2+scoff2, sd1[0],sd3[0],sd0[1],sd2[1],sd4[1],sd1[2]); \
/*	xstore	[srcreg+d1], xmm4	;; Save R1 */ \
/*	xstore	[srcreg+3*d1], xmm5	;; Save R2 */ \
/*	xstore	[srcreg+16], xmm3	;; Save R3 */ \
/*	xstore	[srcreg+2*d1+16], xmm6	;; Save R4 */ \
/*	xstore	[srcreg+4*d1+16], xmm0	;; Save R5 */ \
/*	xstore	[srcreg+d1+32], xmm2	;; Save R6 */ \
	sd3[2] = xmm6;	/* Save R7 */ \
	sd0[3] = xmm2;	/* Save R8 */ \
	sd2[3] = xmm0;	/* Save R9 */ \
	sd4[3] = xmm4;	/* Save R10 */ \
	r5_h10r_h5c_djbunfft_mem(sd0[0],sd0[2],TMP5,TMP6,sd2[0], sd2[2],TMP7,TMP8,sd4[0],sd4[2], screg1,screg1+32,screg2, sd0[0],sd2[0],sd4[0],sd1[1],sd3[1],sd0[2]); \
/*	xstore	[srcreg], xmm4		;; Save R1 */ \
/*	xstore	[srcreg+2*d1], xmm5	;; Save R2 */ \
/*	xstore	[srcreg+4*d1], xmm3	;; Save R3 */ \
/*	xstore	[srcreg+d1+16], xmm6	;; Save R4 */ \
/*	xstore	[srcreg+3*d1+16], xmm0	;; Save R5 */ \
/*	xstore	[srcreg+32], xmm2	;; Save R6 */ \
	sd2[2] = xmm6;	/* Save R7 */ \
	sd4[2] = xmm2;	/* Save R8 */ \
	sd1[3] = xmm0;	/* Save R9 */ \
	sd3[3] = xmm4;	/* Save R10 */ \
	srcreg += srcinc; \
}

/* Like r5_h5cl_ten_reals_five_complex_djbunfft but uses 2 sin/cos ptrs */ \
/* in case where first levels of pass 2 were radix-3. */ \
\

#define r5_h5cl_2sc_ten_reals_five_complex_djbunfft(srcreg,srcinc,d1,screg1,scoff1,screg2,scoff2,screg3,scoff3) { \
	vec2f64 *sd0 = (vec2f64*)(srcreg); \
	vec2f64 *sd1 = (vec2f64*)(srcreg+(d1)); \
	vec2f64 *sd2 = (vec2f64*)(srcreg+2*(d1)); \
	vec2f64 *sd3 = (vec2f64*)(srcreg+3*(d1)); \
	vec2f64 *sd4 = (vec2f64*)(srcreg+4*(d1)); \
	vec2f64 TMP5, TMP6, TMP7, TMP8; \
	xmm6 = sd1[0];	/* Load R2 */ \
	xmm7 = sd1[2];	/* Load I2 */ \
	TMP5 = xmm6; \
	TMP6 = xmm7; \
	xmm6 = sd3[0];	/* Load R4 */ \
	xmm7 = sd3[2];	/* Load I4 */ \
	TMP7 = xmm6; \
	TMP8 = xmm7; \
	r5_h10r_h5c_djbunfft_mem(sd0[1],sd0[3],sd1[1],sd1[3],sd2[1], sd2[3],sd3[1],sd3[3],sd4[1],sd4[3], screg1+scoff1,screg2+scoff2,screg3+scoff3, sd1[0],sd3[0],sd0[1],sd2[1],sd4[1],sd1[2]); \
/*	xstore	[srcreg+d1], xmm4	;; Save R1 */ \
/*	xstore	[srcreg+3*d1], xmm5	;; Save R2 */ \
/*	xstore	[srcreg+16], xmm3	;; Save R3 */ \
/*	xstore	[srcreg+2*d1+16], xmm6	;; Save R4 */ \
/*	xstore	[srcreg+4*d1+16], xmm0	;; Save R5 */ \
/*	xstore	[srcreg+d1+32], xmm2	;; Save R6 */ \
	sd3[2] = xmm6;	/* Save R7 */ \
	sd0[3] = xmm2;	/* Save R8 */ \
	sd2[3] = xmm0;	/* Save R9 */ \
	sd4[3] = xmm4;	/* Save R10 */ \
	r5_h10r_h5c_djbunfft_mem(sd0[0],sd0[2],TMP5,TMP6,sd2[0], sd2[2],TMP7,TMP8,sd4[0],sd4[2], screg1,screg2,screg3, sd0[0],sd2[0],sd4[0],sd1[1],sd3[1],sd0[2]); \
/*	xstore	[srcreg], xmm4		;; Save R1 */ \
/*	xstore	[srcreg+2*d1], xmm5	;; Save R2 */ \
/*	xstore	[srcreg+4*d1], xmm3	;; Save R3 */ \
/*	xstore	[srcreg+d1+16], xmm6	;; Save R4 */ \
/*	xstore	[srcreg+3*d1+16], xmm0	;; Save R5 */ \
/*	xstore	[srcreg+32], xmm2	;; Save R6 */ \
	sd2[2] = xmm6;	/* Save R7 */ \
	sd4[2] = xmm2;	/* Save R8 */ \
	sd1[3] = xmm0;	/* Save R9 */ \
	sd3[3] = xmm4;	/* Save R10 */ \
	srcreg += srcinc; \
}


#define r5_h10r_h5c_djbunfft_mem(memr1p,memi1p,memr2,memi2,memr3p,memi3p,memr4,memi4,memr5p,memi5p, screg1,screg2,screg3, dstr1,dstr2,dstr3,dstr4,dstr5,dstr6) { \
	double f0,f1,f2,f3,f4,f5,f6,f7; \
	double TMP1, TMP2, TMP3, TMP4; \
	double TMP2R, TMP2I, TMP3R, TMP3I;\
	double TMP4R, TMP4I, TMP5R, TMP5I;\
	vec2f64 memr3 = memr3p; \
	vec2f64 memi3 = memi3p; \
	vec2f64 memr5 = memr5p; \
	vec2f64 memi5 = memi5p; \
	/* Do the five complex part */ \
\
	f0 = f64ptr(screg1+16);	/* Load cosine/sine for w^n */ \
	f2 = memr2[1];			/* Load R2 */ \
	f4 = memi2[1];			/* I2 */ \
	f5 = memr5[1];			/* Load R5 */ \
	f6 = memi5[1];			/* I5 */ \
	f1 = f2 * f0 + f4;			/* A2 = R2 * cosine/sine+ I2 */ \
	f7 = f5 * f0 - f6;		/* A5 = R5 * cosine/sine - I5 */ \
	f4 = f4 * f0 - f2;		/* B2 = I2 * cosine/sine - R2 */ \
	f6 = f6 * f0 + f5;		/* B5 = I5 * cosine/sine + R5 */ \
	f3 = f64ptr(screg1); \
	f1 *= f3;				/* A2 = A2 * sine (new R2) */ \
	f4 *= f3;				/* B2 = B2 * sine (new I2) */ \
	f7 *= f3;				/* A5 = A5 * sine (new R5) */ \
	f6 *= f3;				/* B5 = B5 * sine (new I5) */ \
\
	f2 = f1;				/* Copy R2 */ \
	f1 -= f7;				/* r25s=r2-r5 */ \
	f7 += f2;				/* r25a=r2+r5 */ \
	f2 = f4;				/* Copy I2 */ \
	f4 -= f6;				/* i25s=i2-i5 */ \
	f6 += f2;				/* i25a=i2+i5 */ \
\
	TMP1 = f1;				/* Save r2-r5 */ \
	TMP2 = f6;				/* Save i2+i5 */ \
\
	f5 = memr4[1];			/* Load R4 */ \
	f3 = memr3[1];			/* Load R3 */ \
	f2 = memi4[1];			/* I4 */ \
	f0 = memi3[1];			/* I3 */ \
	f1 = f64ptr(screg2+16);	/* Load cosine/sine for w^2n */ \
	f5 = f5 * f1 - f2;		/* A4 = R4 * cosine/sine - I4 */ \
	f3 = f3 * f1 + f0;		/* A3 = R3 * cosine/sine + I3 */ \
	f2 = f2 * f1 + memr4[1];	/* B4 = I4 * cosine/sine */ \
	f0 = f0 * f1 - memr3[1];	/* B3 = I3 * cosine/sine */ \
	f6 = f64ptr(screg2); \
	f5 *= f6;				/* A4 = A4 * sine (new R4) */ \
	f2 *= f6;				/* B4 = B4 * sine (new I4) */ \
	f3 *= f6;				/* A3 = A3 * sine (new R3) */ \
	f0 *= f6;				/* B3 = B3 * sine (new I3) */ \
\
	f1 = f3;				/* Copy R3 */ \
	f3 -= f5;				/* r34s=r3-r4 */ \
	f5 += f1;				/* r34a=r3+r4 */ \
	TMP3 = f3; \
	f6 = f0;				/* Copy I3 */ \
	f0 -= f2;				/* i34s=i3-i4 */ \
	f2 += f6;				/* i34a=i3+i4 */ \
	TMP4 = f2; \
\
	f1 = f7 * CONST_P309;	/* cos2*r25a */ \
	f1 += f5 * CONST_M809;	/* cos2*r25a + cos4*r34a */ \
	f3 = f4 * CONST_P951;	/* sin2*i25s */ \
	f3 += f0 * CONST_P588;	/* t2=sin2*i25s + sin4*i34s */ \
	f2 = f7 * CONST_M809;	/* cos4*r25a */ \
	f2 += f5 * CONST_P309;	/* cos4*r25a+cos2*r34a */ \
	f4 = f4 * CONST_P588;		/* sin4*i25s */ \
	f4 -= f0 * CONST_P951;	/* t4=sin4*i25s-sin2*i34s */ \
	f7 += f5;				/* r25a + r34a */ \
	vec2f64 memr1 = memr1p; \
	f6 = memr1[1]; \
	f1 += f6;				/* t1=cos2*r25a + cos4*r34a + r1 */ \
	f7 += f6;				/* outr(0) = r1 + r25a + r34a */ \
	f2 += f6;				/* t3=cos4*r25a+cos2*r34a+r1 */ \
	f6 = f1; \
	f1 -= f3;				/* outr(4)=t1-t2 */ \
	f3 += f6;				/* outr(1)=t1+t2 */ \
	f5 = f2; \
	f2 -= f4;				/* outr(3)=t3-t4 */ \
	f4 += f5;				/* outr(2)=t3+t4 */ \
\
	dstr1[1] = f7; \
	dstr2[1] = f3;			/* Save new r2 */ \
	dstr3[1] = f4;			/* Save new r3 */ \
	dstr4[1] = f2;			/* Save new r4 */ \
	dstr5[1] = f1;			/* Save new r5 */ \
\
	f0 = TMP1;				/* r25s=r2-r5 */ \
	f2 = TMP2;				/* i25a=i2+i5 */ \
	f1 = TMP3;				/* r34s=r3-r4 */ \
	f3 = TMP4;				/* i34a=i3+i4 */ \
\
	f5 = f2; \
	f5 += f3;				/* i25a+i34a */ \
	f6 = f2 * CONST_P309;	/* cos2*i25a */ \
	f6 += f3 * CONST_M809;	/* cos2*i25a + cos4*i34a */ \
	f2 = f2 * CONST_M809;	/* cos4*i25a */ \
	f2 += f3 * CONST_P309;	/* cos4*i25a + cos2*i34a */ \
	f4 = f0 * CONST_P951;	/* sin2*r25s */ \
	f4 += f1 * CONST_P588;	/* t6=sin2*r25s + sin4*r34s */ \
	f0 = f0 * CONST_P588;	/* sin4*r25s */ \
	f0 -= f1 * CONST_P951;	/* t8=sin4*r25s - sin2*r34s */ \
	vec2f64 memi1 = memi1p; \
	f3 = memi1[1];			/* I1 */ \
	f6 += f3;				/* t5=cos2*i25a + cos4*i34a + i1 */ \
	f2 += f3;				/* t7=cos4*i25a + cos2*i34a + i1 */ \
\
	f7 = f6; \
	f6 -= f4;				/* outi(1)=t5-t6 */ \
	f4 += f7;				/* outi(4)=t5+t6 */ \
	f1 = f2; \
	f2 -= f0;				/* outi(2)=t7-t8 */ \
	f0 += f1;				/* outi(3)=t7+t8 */ \
	f5 += f3;				/* outi(0)=i1+i25a+i34a */ \
\
;	/* Do the ten reals part */ \
\
	xmm6[1] = f6; \
	xmm4[1] = f4; \
	xmm2[1] = f2; \
	xmm0[1] = f0; \
	xmm5[1] = f5; \
\
	vec2f64 sc30 = xptr(screg3); \
	f5 = f64ptr(screg1+16); \
	f1 = memr2[0];				/* R2 */ \
	f0 = memr3[0];				/* R3 */ \
	f3 = memi2[0];				/* I2 */ \
	f7 = memi3[0];				/* I3 */ \
	f1 = f1 * sc30[1] + f3;		/* A2 = R2 * cosine/sine + I2 */ \
	f0 = memr3[0] * f5 + f7;	/* A3 = R3 * cosine/sine + I3 */ \
	f3 = f3 * sc30[1] - memr2[0];	/* B2 = I2 * cosine/sine - R2 */ \
	f7 = f7 * f5 - memr3[0];	/* B3 = I3 * cosine/sine - R3 */ \
	\
	f5 = f64ptr(screg1); \
	f1 *= sc30[0];				/* A2 = A2 * sine (new R2) */ \
	f0 *= f5;					/* A3 = A3 * sine (new R3) */ \
	f3 *= sc30[0];				/* B2 = B2 * sine (new I2) */ \
	f7 *= f5;					/* B3 = B3 * sine (new I3) */ \
	TMP2R = f1; \
	TMP2I = f3; \
	TMP3R = f0; \
	TMP3I = f7; \
\
	vec2f64 sc31 = xptr(screg3+16); \
	f5 = f64ptr(screg2+16); \
	f1 = memr4[0];				/* R4 */ \
	f0 = memr5[0];				/* R5 */ \
	f3 = memi4[0];				/* I4 */ \
	f1 = f1 * sc31[1] + f3;		/* A4 = R4 * cosine/sine + I4 */ \
	f7 = memi5[0];				/* I5 */ \
	f0 = f0 * f5 + f7;			/* A5 = R5 * cosine/sine + I5 */ \
	f3 = f3 * sc31[1] - memr4[0];	/* B4 = I4 * cosine/sine - R4 */ \
	f7 = f7 * f5 - memr5[0];	/* B5 = I5 * cosine/sine - R5 */ \
	f5 = f64ptr(screg2); \
	f1 *= sc31[0];				/* A4 = A4 * sine (new R4) */ \
	f0 *= f5;					/* A5 = A5 * sine (new R5) */ \
	f3 *= sc31[0];				/* B4 = B4 * sine (new I4) */ \
	f7 *= f5;					/* B5 = B5 * sine (new I5) */ \
	TMP4R = f1; \
	TMP4I = f3; \
	TMP5R = f0; \
	TMP5I = f7; \
\
	f7 =  TMP2R;			/* R2 */ \
	f7 += TMP5R;			/* R2+R5 */ \
	f1 =  TMP3R;			/* R3 */ \
	f1 += TMP4R;			/* R3+R4 */ \
	f3 = memr1[0];			/* R1_1 */ \
	f0 = f7 * CONST_P309 + f3;	/* new R3 = .309*(R2+R5) + R1_1 */ \
	f6 = f7 * CONST_M809 + f3;	/* new R5 = -.809*(R2+R5) + R1_1 */ \
	f7 += f3;				/* new R1 = R1_1+R2+R5 */ \
	f0 += f1 * CONST_M809;	/* new R3 += -.809*(R3+R4) */ \
	f6 += f1 * CONST_P309;	/* new R5 += .309*(R3+R4) */ \
	f7 += f1;				/* final R1 = R1_1+R2+R5+R3+R4 */ \
	dstr1[0] = f7;			/* Save final R1 */ \
\
	f3  = TMP2I;			/* I2 */ \
	f3 -= TMP5I;			/* I2-I5 */ \
	f1  = TMP3I;			/* I3 */ \
	f1 -= TMP4I;			/* I3-I4 */ \
	f7 = f3 * CONST_P951;	/* tmp1 = .951*(I2-I5) */ \
	f7 += f1 * CONST_P588;	/* tmp1 += .588*(I3-I4) */ \
	f3 = f3 * CONST_P588;	/* tmp2 = .588*(I2-I5) */ \
	f3 -= f1 * CONST_P951;	/* tmp2 -= .951*(I3-I4) */ \
\
	f1 = f0;				/* Copy new R3 */ \
	f0 -= f7;				/* final R9 = new R3 - tmp1 */ \
	f7 += f1;				/* final R3 = new R3 + tmp1 */ \
\
	f1 = f6;				/* Copy new R5 */ \
	f6 -= f3;				/* final R7 = new R5 - tmp2 */ \
	f3 += f1;				/* final R5 = new R5 + tmp2 */ \
\
	f5 =  TMP3R;			/* R3 */ \
	f5 -= TMP4R;			/* R3-R4 */ \
	f1 =  TMP2R;			/* R2 */ \
	f1 -= TMP5R;			/* R2-R5 */ \
	dstr3[0] = f7;		/* Save final R3 */ \
	f7 = memi1[0];	/* R1_2 */ \
	f4 = f5 * CONST_P309 + f7;	/* new R2 = .309*(R3-R4) + R1_2 */ \
	f2 = f5 * CONST_M809 + f7;	/* new R4 = -.809*(R3-R4) + R1_2 */ \
	f5 += f7;				/* new R6 = R1_2+(R3-R4) */ \
	f4 -= f1 * CONST_M809;	/* new R2 -= -.809*(R2-R5) */ \
	f2 -= f1 * CONST_P309;	/* new R4 -= .309*(R2-R5) */ \
	f5 -= f1;				/* final R6 = R6+(R3-R4)-(R2-R5) */ \
\
	dstr5[0] = f3;			/* Save final R5 */ \
	xmm5[0] = f5; \
	dstr6 = xmm5;			/* Save final R6 */ \
\
	f3  = TMP2I;	/* I2 */ \
	f3 += TMP5I;	/* I2+I5 */ \
	f1  = TMP3I;	/* I3 */ \
	f1 += TMP4I;	/* I3+I4 */ \
	f7 = f3 * CONST_P951;	/* tmp2 = .951*(I2+I5) */ \
	f7 -= f1 * CONST_P588;	/* tmp2 -= .588*(I3+I4) */ \
	f3 = f3 * CONST_P588;	/* tmp1 = .588*(I2+I5) */ \
	f3 += f1 * CONST_P951;	/* tmp1 += .951*(I3+I4) */ \
\
	f1 = f4;				/* Copy new R2 */ \
	f4 -= f3;				/* Final R10 = new R2 - tmp1 */ \
	f3 += f1;				/* Final R2 = new R2 + tmp1 */ \
\
	f1 = f2;				/* Copy new R4 */ \
	f2 -= f7;				/* Final R8 = new R4 - tmp2 */ \
	f7 += f1;				/* Final R4 = new R4 + tmp2 */ \
\
	dstr2[0] = f3;			/* Save final R2 */ \
	dstr4[0] = f7;			/* Save final R4 */ \
\
	xmm6[0] = f6; \
	xmm4[0] = f4; \
	xmm2[0] = f2; \
	xmm0[0] = f0; \
	/* xmm5[0] = f5; unused */ \
}


/* */ \
/* ************************************* 20-reals-first-fft variants ****************************************** */ \
/* */ \
\
/* This should in theory be faster than an 8-real step 1 followed by a 5-complex (or 10-real) step 2. */ \
/* To see this, count the adds and muls to process 80 reals as either */ \
/*    1)  10 * eight-reals, 2 ten-reals and 6 five-complex, or */ \
/*    2)  4 * twenty-reals, 1/2 * (one eight-reals and 9 four-complex) */ \
\
/* These macros operate on twenty reals doing 4.32 levels of the FFT.  The output is */ \
/* 2 reals and 9 complex numbers. */ \
\
/*r5_x5cl_20_reals_first_fft_preload MACRO */ \
/*	r5_x5cl_20_reals_first_fft_cmn_preload */ \
/*	ENDM */ \
/*r5_x5cl_20_reals_first_fft MACRO srcreg,srcinc,d1,screg */ \
/*	r5_x5cl_20_reals_first_fft_cmn srcreg,rbx,srcinc,d1,screg */ \
/*	ENDM */ \

#define r5_x5cl_20_reals_first_fft_scratch_preload \
	r5_x5cl_20_reals_first_fft_cmn_preload; \
\

#define r5_x5cl_20_reals_first_fft_scratch(srcreg,srcinc,d1,screg) \
	r5_x5cl_20_reals_first_fft_cmn(srcreg,srcinc,d1,screg); \
\
\
/* To calculate a 20-reals FFT, we calculate 20 complex values in a brute force way (using a shorthand notation): */ \
/* r1 + r2 + ... + r20	*  w^0000000000... */ \
/* r1 + r2 + ... + r20	*  w^0123456789A... */ \
/* r1 + r2 + ... + r20	*  w^02468ACE.... */ \
/*    ... */ \
/* r1 + r2 + ... + r20	*  w^...A987654321 */ \
/* Note that Hermetian symmetry means we won't need to calculate the last 10 complex values. */ \
/* */ \
/* The sin/cos values (w = 20th root of unity) are: */ \
/* w^1 = .951 + .309i */ \
/* w^2 = .809 + .588i */ \
/* w^3 = .588 + .809i */ \
/* w^4 = .309 + .951i */ \
/* w^5 = 0 + 1i */ \
/* w^6 = -.309 + .951i */ \
/* w^7 = -.588 + .809i */ \
/* w^8 = -.809 + .588i */ \
/* w^9 = -.951 + .309i */ \
/* w^10 = -1 */ \
/* */ \
/* Applying the sin/cos values above (and noting that combining r2 and r20, r3 and r19, etc. will simplify calculations): */ \
/* reals: */ \
/* r1     +(r2+r20)     +(r3+r19)     +(r4+r18)     +(r5+r17) + (r6+r16)     +(r7+r15)     +(r8+r14)     +(r9+r13)     +(r10+r12) + r11 */ \
/* r1 +.951(r2+r20) +.809(r3+r19) +.588(r4+r18) +.309(r5+r17)            -.309(r7+r15) -.588(r8+r14) -.809(r9+r13) -.951(r10+r12) - r11 */ \
/* r1 +.809(r2+r20) +.309(r3+r19) -.309(r4+r18) -.809(r5+r17) - (r6+r16) -.809(r7+r15) -.309(r8+r14) +.309(r9+r13) +.809(r10+r12) + r11 */ \
/* r1 +.588(r2+r20) -.309(r3+r19) -.951(r4+r18) -.809(r5+r17)            +.809(r7+r15) +.951(r8+r14) +.309(r9+r13) -.588(r10+r12) - r11 */ \
/* r1 +.309(r2+r20) -.809(r3+r19) -.809(r4+r18) +.309(r5+r17) + (r6+r16) +.309(r7+r15) -.809(r8+r14) -.809(r9+r13) +.309(r10+r12) + r11 */ \
/* r1                   -(r3+r19)                   +(r5+r17)                -(r7+r15)                   +(r9+r13)                - r11 */ \
/* r1 -.309(r2+r20) -.809(r3+r19) +.809(r4+r18) +.309(r5+r17) - (r6+r16) +.309(r7+r15) +.809(r8+r14) -.809(r9+r13) -.309(r10+r12) + r11 */ \
/* r1 -.588(r2+r20) -.309(r3+r19) +.951(r4+r18) -.809(r5+r17)            +.809(r7+r15) -.951(r8+r14) +.309(r9+r13) +.588(r10+r12) - r11 */ \
/* r1 -.809(r2+r20) +.309(r3+r19) +.309(r4+r18) -.809(r5+r17) + (r6+r16) -.809(r7+r15) +.309(r8+r14) +.309(r9+r13) -.809(r10+r12) + r11 */ \
/* r1 -.951(r2+r20) +.809(r3+r19) -.588(r4+r18) +.309(r5+r17)            -.309(r7+r15) +.588(r8+r14) -.809(r9+r13) +.951(r10+r12) - r11 */ \
/* r1     -(r2+r20)     +(r3+r19)     -(r4+r18)     +(r5+r17) - (r6+r16)     +(r7+r15)     -(r8+r14)     +(r9+r13)     -(r10+r12) + r11 */ \
/* */ \
/* imaginarys: */ \
/* 0 */ \
/* +.309(r2-r20) +.588(r3-r19) +.809(r4-r18) +.951(r5-r17) + (r6-r16) +.951(r7-r15) +.809(r8-r14) +.588(r9-r13) +.309(r10-r12) */ \
/* +.588(r2-r20) +.951(r3-r19) +.951(r4-r18) +.588(r5-r17)            -.588(r7-r15) -.951(r8-r14) -.951(r9-r13) -.588(r10-r12) */ \
/* +.809(r2-r20) +.951(r3-r19) +.309(r4-r18) -.588(r5-r17) - (r6-r16) -.588(r7-r15) +.309(r8-r14) +.951(r9-r13) +.809(r10-r12) */ \
/* +.951(r2-r20) +.588(r3-r19) -.588(r4-r18) -.951(r5-r17)            +.951(r7-r15) +.588(r8-r14) -.588(r9-r13) -.951(r10-r12) */ \
/*      (r2-r20)                   -(r4-r18)               + (r6-r16)                   -(r8-r14)                   +(r10-r12) */ \
/* +.951(r2-r20) -.588(r3-r19) -.588(r4-r18) +.975(r5-r17)            -.951(r7-r15) +.588(r8-r14) +.588(r9-r13) -.951(r10-r12) */ \
/* +.809(r2-r20) -.951(r3-r19) +.309(r4-r18) +.588(r5-r17) - (r6-r16) +.588(r7-r15) +.309(r8-r14) -.951(r9-r13) +.809(r10-r12) */ \
/* +.588(r2-r20) -.951(r3-r19) +.951(r4-r18) -.588(r5-r17)            +.588(r7-r15) -.951(r8-r14) +.951(r9-r13) -.588(r10-r12) */ \
/* +.309(r2-r20) -.588(r3-r19) +.809(r4-r18) -.951(r5-r17) + (r6-r16) -.951(r7-r15) +.809(r8-r14) -.588(r9-r13) +.309(r10-r12) */ \
/* 0 */ \
/* */ \
/* There are many more symmetries we can take advantage of.   For example, the (r2+/-r20) column */ \
/* always has the same multiplier as the (r10+/-r12) column.  This is true for all the "even" columns. */ \
/* Also the computations for the 2nd row are very similar to the computations for the 10th row, */ \
/* the 3rd row are similar to the 9th, etc.  Finally, note that for the odd columns, there are */ \
/* only two multipliers to apply and can be combined with every fourth column. */ \
/* */ \
/* Lastly, output would normally be 9 complex and 2 reals but the users of this routine */ \
/* expect us to "back up" the 2 reals by one level.  That is: */ \
/*	real #1A:  r1 + r3+r19 + r5+r17 + ... */ \
/*	real #1B:  r2+r20 + r4+r18 + ... */ \
\
/* Store intermediate results in g->u.xmm.XMM_COL_MULTS (an 8KB buffer used in normalization) */ \
\

#if 1
#define r5_x5cl_20_reals_first_fft_cmn_preload \

#define r5_x5cl_20_reals_first_fft_cmn(srcreg,srcinc,d1,screg) { \
	vec2f64 *sd0 = (vec2f64*)(srcreg); \
	vec2f64 *sd1 = (vec2f64*)(srcreg+d1*1); \
	vec2f64 *sd2 = (vec2f64*)(srcreg+d1*2); \
	vec2f64 *sd3 = (vec2f64*)(srcreg+d1*3); \
	vec2f64 *sd4 = (vec2f64*)(srcreg+d1*4); \
	vec2f64 *sc = (vec2f64*)(screg); \
	vec2f64 t0,t1,t2,t3,t4,t5,t6,t7,t8,t9; \
	vec2f64 t10,t11,t12,t13,t14,t15,t16,t17; \
\
;	/* Do the odd columns for the real results */ \
\
	xmm7 = sd4[0];			/* r5 */ \
	xmm0 = sd1[3];			/* r17 */ \
	vec2f64 r5s17 = xmm7 - xmm0; \
	xmm7 += xmm0;			/* r5+r17 */ \
	xmm2 = sd0[0];			/* r1 */ \
	xmm0 = xmm7 + xmm2;		/* r1+(r5+r17) */ \
	xmm1 = xmm2 + xmm7 * CONST2_P309;	/* r1+.309(r5+r17) */ \
	xmm2 = xmm2 - xmm7 * CONST2_P809;	/* r1-.809(r5+r17) */ \
\
	xmm4 = sd3[1];			/* r9 */ \
	xmm5 = sd2[2];			/* r13 */ \
	vec2f64 r9s13 = xmm4 - xmm5; \
	xmm4 += xmm5;			/* r9+r13 */ \
	xmm0 += xmm4;			/* r1+(r5+r17)+(r9+r13) */ \
	xmm1 -= xmm4 * CONST2_P809;	/* r1+.309(r5+r17)-.809(r9+r13) */ \
	xmm2 += xmm4 * CONST2_P309;	/* r1-.809(r5+r17)+.309(r9+r13) */ \
\
	xmm7 = sd2[0];			/* r3 */ \
	xmm4 = sd3[3];			/* r19 */ \
	vec2f64 r3s19 = xmm7 - xmm4; \
	xmm7 += xmm4;			/* r3+r19 */ \
	xmm4 = sd0[2];			/* r11 */ \
	xmm5 = xmm7 + xmm4;		/* (r3+r19)+r11 */ \
	xmm6 = xmm7 * CONST2_P809 - xmm4;	/* .809(r3+r19)-r11 */ \
	xmm7 = xmm7 * CONST2_P309 + xmm4;	/* .309(r3+r19)+r11 */ \
\
	xmm4 = sd1[1];			/* r7 */ \
	xmm3 = sd4[2];			/* r15 */ \
	vec2f64 r7s15 = xmm4 - xmm3; \
	xmm4 += xmm3;			/* r7+r15 */ \
	xmm5 += xmm4;			/* (r3+r19)+(r7+r15)+r11 */ \
	xmm6 -= xmm4 * CONST2_P309;	/* .809(r3+r19)-.309(r7+r15)-r11 */ \
	xmm7 -= xmm4 * CONST2_P809;	/* .309(r3+r19)-.809(r7+r15)+r11 */ \
\
	xmm4 = xmm0;			/* Copy 1-mod-4-cols row #1 */ \
	xmm0 -= xmm5;			/* Real odd-cols row #6 (final real #6) */ \
	xmm5 += xmm4;			/* Real odd-cols row #1 (final real #1A) */ \
\
	xmm3 = xmm1;			/* Copy 1-mod-4-cols row #2 */ \
	xmm1 -= xmm6;			/* Real odd-cols row #5 */ \
	xmm6 += xmm3;			/* Real odd-cols row #2 */ \
\
	xmm4 = xmm2;			/* Copy 1-mod-4-cols row #3 */ \
	xmm2 -= xmm7;			/* Real odd-cols row #4 */ \
	xmm7 += xmm4;			/* Real odd-cols row #3 */ \
\
	t4 =xmm0;				/* Real #6 */ \
	sd0[0] = xmm5;			/* Final real #1A */ \
	t3 = xmm1;				/* Real odd-cols row #5 */ \
	t0 = xmm6;				/* Real odd-cols row #2 */ \
	t2 = xmm2;				/* Real odd-cols row #4 */ \
	t1 = xmm7;				/* Real odd-cols row #3 */ \
\
;	/* Do the even columns for the real results */ \
\
	xmm7 = sd1[0];			/* r2 */ \
	xmm2 = sd4[3];			/* r20 */ \
	vec2f64 r2s20 = xmm7 - xmm2; \
	xmm7 += xmm2;			/* r2+r20 */ \
	xmm0 = sd4[1];			/* r10 */ \
	xmm2 = sd1[2];			/* r12 */ \
	vec2f64 r10s12 = xmm0 - xmm2; \
	xmm0 += xmm2;			/* r10+r12 */ \
	xmm2 = xmm7; \
	xmm7 -= xmm0;			/* (r2+r20)-(r10+r12) */ \
	xmm0 += xmm2;			/* (r2+r20)+(r10+r12) */ \
\
\
	xmm4 = sd3[0];			/* r4 */ \
	xmm2 = sd2[3];			/* r18 */ \
	vec2f64 r4s18 = xmm4 - xmm2; \
	xmm4 += xmm2;			/* r4+r18 */ \
	xmm1 = sd2[1];			/* r8 */ \
	xmm2 = sd3[2];			/* r14 */ \
	vec2f64 r8s14 = xmm1 - xmm2; \
	xmm1 += xmm2;			/* r8+r14 */ \
	xmm2 = xmm4; \
	xmm4 -= xmm1;			/* (r4+r18)-(r8+r14) */ \
	xmm1 += xmm2;			/* (r4+r18)+(r8+r14) */ \
\
	xmm6 = xmm7 * CONST2_P951;	/* .951((r2+r20)-(r10+r12)) */ \
	xmm6 += xmm4 * CONST2_P588;	/* .951((r2+r20)-(r10+r12))+.588((r4+r18)-(r8+r14)) */ \
	xmm7 = xmm7 * CONST2_P588;	/* .588((r2+r20)-(r10+r12)) */ \
	xmm7 -= xmm4 * CONST2_P951;	/* .588((r2+r20)-(r10+r12))-.951((r4+r18)-(r8+r14)) */ \
\
	t5 = xmm6;				/* Save real even-cols row #2 */ \
	t6 = xmm7;				/* Save real even-cols row #4 */ \
\
	xmm4 = xmm0 * CONST2_P809;	/* .809((r2+r20)+(r10+r12)) */ \
	xmm5 = xmm0 * CONST2_P309;	/* .309((r2+r20)+(r10+r12)) */ \
\
	xmm0 += xmm1;				/* ((r2+r20)+(r10+r12))+((r4+r18)+(r8+r14)) */ \
	xmm3 *= xmm1 * CONST2_P309;	/* .309((r4+r18)+(r8+r14)) */ \
	xmm4 -= xmm1 * CONST2_P309;	/* .809((r2+r20)+(r10+r12))-.309((r4+r18)+(r8+r14)) */ \
	xmm5 -= xmm1 * CONST2_P809;	/* .309((r2+r20)+(r10+r12))-.809((r4+r18)+(r8+r14)) */ \
\
	xmm2 = sd0[1];			/* r6 */ \
	xmm3 = sd0[3];			/* r16 */ \
	xmm2 -= xmm3;			/* r6-r16 */ \
	xmm3 += sd0[1];			/* r6+r16 */ \
	xmm0 += xmm3;			/* ((r2+r20)+(r10+r12))+((r4+r18)+(r8+r14))+(r6+r16) */ \
	xmm4 -= xmm3;			/* .809((r2+r20)+(r10+r12))-.309((r4+r18)+(r8+r14))-(r6+r16) */ \
	xmm5 += xmm3;			/* .309((r2+r20)+(r10+r12))-.809((r4+r18)+(r8+r14))+(r6+r16) */ \
\
	sd0[1] = xmm0;			/* Save final real #1B (real even-cols row #1) */ \
	t7 = xmm4;				/* Save real even-cols row #3 */ \
	t8 = xmm5;				/* Save real even-cols row #5 */ \
\
;	/* Do the even columns for the imaginary results */ \
\
	xmm0 = r2s20;			/* r2-r20 */ \
	xmm7 = r10s12;			/* r10-r12 */ \
	xmm3 = xmm7; \
	xmm7 += xmm0;			/* (r2-r20)+(r10-r12) */ \
	xmm0 -= xmm3;			/* (r2-r20)-(r10-r12) */ \
\
	xmm4 = xmm7 + xmm2;			/* ((r2-r20)+(r10-r12))+(r6-r16) */ \
	xmm6 = xmm7 * CONST2_P309 + xmm2;	/* .309((r2-r20)+(r10-r12))+(r6-r16) */ \
	xmm7 = xmm7 * CONST2_P809 - xmm2;	/* .809((r2-r20)+(r10-r12))-(r6-r16) */ \
\
	xmm1 = r4s18;			/* r4-r18 */ \
	xmm3 = r8s14;			/* r8-r14 */ \
	xmm2 = xmm3; \
	xmm3 += xmm1;			/* (r4-r18)+(r8-r14) */ \
	xmm1 -= xmm2;			/* (r4-r18)-(r8-r14) */ \
\
	xmm4 -= xmm3;			/* ((r2-r20)+(r10-r12))-((r4-r18)+(r8-r14))+(r6-r16) */ \
	xmm6 += xmm3 * CONST2_P809;	/* .309((r2-r20)+(r10-r12))+.809((r4-r18)+(r8-r14))+(r6-r16) */ \
	xmm7 += xmm3 * CONST2_P309;	/* .809((r2-r20)+(r10-r12))+.309((r4-r18)+(r8-r14))-(r6-r16) */ \
\
	t9 = xmm4;	/* Save imag row #6 */ \
	t10 = xmm6;	/* Save imag even-cols row #2 */ \
	t11 = xmm7;	/* Save imag even-cols row #4 */ \
\
	xmm5 = xmm1 * CONST2_P951;	/* .951((r4-r18)-(r8-r14)) */ \
	xmm5 += xmm0 * CONST2_P588;	/* .588((r2-r20)-(r10-r12))+.951((r4-r18)-(r8-r14)) */ \
	xmm0 = xmm0 * CONST2_P951;	/* .951((r2-r20)-(r10-r12)) */ \
	xmm0 -= xmm1 * CONST2_P588;	/* .951((r2-r20)-(r10-r12))-.588((r4-r18)-(r8-r14)) */ \
\
	t12 = xmm5;	/* Save imag even-cols row #3 */ \
	t13 = xmm0;	/* Save imag even-cols row #5 */ \
\
;	/* Do the odd columns for the imag results */ \
\
	xmm1 = r5s17;			/* r5-r17 */ \
	xmm7 = r9s13;			/* r9-r13 */ \
	xmm0 = xmm1 * CONST2_P951;	/* .951(r5-r17) */ \
	xmm0 += xmm7 * CONST2_P588;			/* .951(r5-r17)+.588(r9-r13) */ \
	xmm1 = xmm1 * CONST2_P588;	/* .588(r5-r17) */ \
	xmm1 -= xmm7 * CONST2_P951;			/* .588(r5-r17)-.951(r9-r13) */ \
\
	xmm7 = r3s19;			/* r3-r19 */ \
	xmm6 = xmm7 * CONST2_P588;	/* .588(r3-r19) */ \
	xmm7 = xmm7 * CONST2_P951;	/* .951(r3-r19) */ \
\
	xmm2 = r7s15;			/* r7-r15 */ \
	xmm6 += xmm2 * CONST2_P951;	/* .588(r3-r19)+.951(r7-r15) */ \
	xmm7 -= xmm2 * CONST2_P588;	/* .951(r3-r19)-.588(r7-r15) */ \
\
	xmm4 = xmm0;			/* Copy 1-mod-4 imag odd-cols row #2 */ \
	xmm0 += xmm6;			/* Imag odd-cols row #2 */ \
	xmm6 -= xmm4;			/* Imag odd-cols row #5 */ \
\
	xmm5 = xmm1;			/* Copy 1-mod-4 imag odd-cols row #4 */ \
	xmm1 += xmm7;			/* Imag odd-cols row #3 */ \
	xmm7 -= xmm5;			/* Imag odd-cols row #4 */ \
\
	t14 = xmm0;	/* Imag odd-cols row #2 */ \
	t17 = xmm6;	/* Imag odd-cols row #5 */ \
	t15 = xmm1;	/* Imag odd-cols row #3 */ \
	t16 = xmm7;	/* Imag odd-cols row #4 */ \
\
;	/* Now combine the even and odd columns then do the post-multiply by twiddle factors. */ \
\
	xmm0 = t0;				/* Real odd-cols row #2 */ \
	xmm1 = t5;				/* Real even-cols row #2 */ \
	xmm0 -= xmm1;			/* Real #10 */ \
	xmm1 += t0;				/* Real #2 */ \
	xmm2 = t10;				/* Imag even-cols row #2 */ \
	xmm3 = t14;				/* Imag odd-cols row #2 */ \
	xmm2 -= xmm3;			/* Imag #10 */ \
	xmm3 += t10;			/* Imag #2 */ \
\
	xmm5 = sc[8*2+1]; \
	xmm7 = xmm0;			/* Copy R10 */ \
	xmm0 = xmm7 * xmm5 - xmm2;	/* A10 = R10 * cosine/sine - I10 */ \
	xmm2 = xmm2 * xmm5 + xmm7;	/* B10 = I10 * cosine/sine + R10 */ \
	xmm0 *= sc[8*2];		/* A10 = A10 * sine (final R10) */ \
	xmm2 *= sc[8*2];		/* B10 = B10 * sine (final I10) */ \
	sd4[2] = xmm0;			/* Save final R10 */ \
	sd4[3] = xmm2;			/* Save final I10 */ \
\
	xmm4 = sc[1]; \
	xmm7 = sc[0]; \
	xmm6 = xmm1;			/* Copy R2 */ \
	xmm1 = xmm6 * xmm4 - xmm3;			/* A2 = R2 * cosine/sine - I2 */ \
	xmm3 = xmm3 * xmm4 + xmm6;			/* B2 = I2 * cosine/sine + R2 */ \
	xmm1 *= sc[0];			/* A2 = A2 * sine (final R2) */ \
	xmm3 *= sc[0];			/* B2 = B2 * sine (final I2) */ \
\
	sd0[2] = xmm1;			/* Save final R2 */ \
	sd0[3] = xmm3;			/* Save final I2 */ \
\
	xmm0 = t1;				/* Real odd-cols row #3 */ \
	xmm1 = t7;				/* Real even-cols row #3 */ \
	xmm0 -= xmm1;			/* Real #9 */ \
	xmm1 += t1;				/* Real #3 */ \
	xmm2 = t12;				/* Imag even-cols row #3 */ \
	xmm3 = t15;				/* Imag odd-cols row #3 */ \
	xmm2 -= xmm3;			/* Imag #9 */ \
	xmm3 += t12;			/* Imag #3 */ \
\
	xmm5 = sc[7*2+1];		/* cosine/sine */ \
	xmm7 = xmm0;			/* Copy R9 */ \
	xmm0 = xmm7 * xmm5 - xmm2;	/* A9 = R9 * cosine/sine - I9 */ \
	xmm2 = xmm2 * xmm5 + xmm7;	/* B9 = I9 * cosine/sine + R9 */ \
	xmm0 *= sc[7*2];		/* A9 = A9 * sine (final R9) */ \
	xmm2 *= sc[7*2];		/* B9 = B9 * sine (final I9) */ \
\
	xmm4 = sc[2+1];			/* cosine/sine */ \
	xmm6 = xmm1;			/* Copy R3 */ \
	xmm1 = xmm6 * xmm4 - xmm3;	/* A3 = R3 * cosine/sine - I3 */ \
	xmm3 = xmm3 * xmm4 + xmm6;	/* B3 = I3 * cosine/sine + R3 */ \
	xmm1 *= sc[2];			/* A3 = A3 * sine (final R3) */ \
	xmm3 *= sc[2];			/* B3 = B3 * sine (final I3) */ \
\
	sd4[0] = xmm0;			/* Save final R9 */ \
	sd4[1] = xmm2;			/* Save final I9 */ \
	sd1[0] = xmm1;			/* Save final R3 */ \
	sd1[1] = xmm3;			/* Save final I3 */ \
\
	xmm0 = t2;				/* Real odd-cols row #4 */ \
	xmm1 = t6;				/* Real even-cols row #4 */ \
	xmm0 -= xmm1;			/* Real #8 */ \
	xmm1 += t2;				/* Real #4 */ \
	xmm2 = t11;				/* Imag even-cols row #4 */ \
	xmm3 = t16;				/* Imag odd-cols row #4 */ \
	xmm2 -= xmm3;			/* Imag #8 */ \
	xmm3 += t11;			/* Imag #4 */ \
\
	xmm5 = sc[6*2+1];		/* cosine/sine */ \
	xmm7 = xmm0;			/* Copy R8 */ \
	xmm0 = xmm7 * xmm5 - xmm2;	/* A8 = R8 * cosine/sine - I8 */ \
	xmm2 = xmm2 * xmm5 + xmm7;	/* B8 = I8 * cosine/sine + R8 */ \
	xmm0 *= sc[6*2];		/* A8 = A8 * sine (final R8) */ \
	xmm2 *= sc[6*2];		/* B8 = B8 * sine (final I8) */ \
\
	xmm4 = sc[2*2+1];		/* cosine/sine */ \
	xmm6 = xmm1;			/* Copy R4 */ \
	xmm1 = xmm6 * xmm4 - xmm3;	/* A4 = R4 * cosine/sine - I4 */ \
	xmm3 = xmm3 * xmm4 + xmm6;	/* B4 = I4 * cosine/sine + R4 */ \
	xmm1 *= sc[2*2];		/* A4 = A4 * sine (final R4) */ \
	xmm3 *= sc[2*2];		/* B4 = B4 * sine (final I4) */ \
\
	sd3[2] = xmm0;			/* Save final R8 */ \
	sd3[3] = xmm2;			/* Save final I8 */ \
	sd1[2] = xmm1;			/* Save final R4 */ \
	sd1[3] = xmm3;			/* Save final I4 */ \
\
	xmm0 = t3;				/* Real odd-cols row #5 */ \
	xmm1 = t8;				/* Real even-cols row #5 */ \
	xmm0 -= xmm1;			/* Real #7 */ \
	xmm1 += t3;				/* Real #5 */ \
	xmm2 = t13;				/* Imag even-cols row #5 */ \
	xmm3 = t17;				/* Imag odd-cols row #5 */ \
	xmm2 -= xmm3;			/* Imag #7 */ \
	xmm3 += t13;			/* Imag #5 */ \
\
	xmm4 = sc[5*2+1];		/* cosine/sine */ \
	xmm6 = xmm0;			/* Copy R7 */ \
	xmm0 = xmm6 * xmm4 - xmm2;	/* A7 = R7 * cosine/sine - I7 */ \
	xmm2 = xmm2 * xmm4 + xmm6;	/* B7 = I7 * cosine/sine + R7 */ \
	xmm0 *= sc[5*2];		/* A7 = A7 * sine (final R7) */ \
	xmm2 *= sc[5*2];		/* B7 = B7 * sine (final I7) */ \
\
	xmm5 = sc[3*2+1];		/* cosine/sine */ \
	xmm7 = xmm1;			/* Copy R5 */ \
	xmm1 = xmm7 * xmm5 - xmm3;	/* A5 = R5 * cosine/sine - I5 */ \
	xmm3 = xmm3 * xmm5 + xmm7;	/* B5 = I5 * cosine/sine + R5 */ \
	xmm1 *= sc[3*2];		/* A5 = A5 * sine (final R5) */ \
	xmm3 *= sc[3*2];		/* B5 = B5 * sine (final I5) */ \
\
	sd3[0] = xmm0;			/* Save final R7 */ \
	sd3[1] = xmm2;			/* Save final I7 */ \
\
	xmm6 = t4;				/* Real #6 */ \
	xmm2 = t9;				/* Imag #6 */ \
	xmm4 = sc[4*2+1];		/* cosine/sine */ \
	xmm0 = xmm6 * xmm4 - xmm2;	/* A6 = R6 * cosine/sine - I6 */ \
	xmm2 = xmm2 * xmm4 + xmm6;	/* B6 = I6 * cosine/sine + R6 */ \
	xmm0 *= sc[4*2];		/* A6 = A6 * sine (final R6) */ \
	xmm2 *= sc[4*2];		/* B6 = B6 * sine (final I6) */ \
\
	sd2[0] = xmm1;			/* Save final R5 */ \
	sd2[1] = xmm3;			/* Save final I5 */ \
	sd2[2] = xmm0;			/* Save final R6 */ \
	sd2[3] = xmm2;			/* Save final I6 */ \
\
	srcreg += srcinc; \
}

/* 64-bit version of the above using more registers for better instruction scheduling */
#else //DEF X86_64;

#define r5_x5cl_20_reals_first_fft_cmn_preload \
/*	xmm15 = CONST2_P309;*/ \
/*	xmm14 = CONST2_P809;*/ \

#define r5_x5cl_20_reals_first_fft_cmn(srcreg,srcinc,d1,screg) { \
	vec2f64 *sd0 = (vec2f64*)(srcreg); \
	vec2f64 *sd1 = (vec2f64*)(srcreg+d1*1); \
	vec2f64 *sd2 = (vec2f64*)(srcreg+d1*2); \
	vec2f64 *sd3 = (vec2f64*)(srcreg+d1*3); \
	vec2f64 *sd4 = (vec2f64*)(srcreg+d1*4); \
	vec2f64 *sc = (vec2f64*)(screg); \
	vec2f64 xmm8,xmm9,xmm10,xmm11,xmm12,xmm13,xmm14; \
	vec2f64 t0,t1,t2,t3,t4,t5,t6,t7,t8,t9; \
	vec2f64 t10,t11,t12,t13,t14,t15; \
\
;	/* Do the odd columns */ \
\
	xmm0 = sd4[0];			/* r5 */ \
	xmm1 = sd1[3];			/* r17 */ \
	xmm2 = xmm0;			/* Copy r5 */ \
	xmm0 += xmm1;			/* r5+r17				; 1-3 */ \
\
	xmm3 = sd3[1];			/* r9 */ \
	xmm4 = sd2[2];			/* r13 */ \
	xmm5 = xmm3;			/* Copy r9 */ \
	xmm3 += xmm4;			/* r9+r13				; 2-4 */ \
\
	xmm2 -= xmm1;			/* r5-r17				; 3-5			avail 1 */ \
\
	xmm5 -= xmm4;			/* r9-r13				; 4-6			avail 1,4 */ \
	xmm6 = CONST2_P309; \
	xmm6 *= xmm0;			/* .309(r5+r17)				; 4-8 */ \
\
	xmm7 = sd2[0];			/* r3 */ \
	xmm8 = sd3[3];			/* r19 */ \
	xmm9 = xmm7;			/* Copy r3 */ \
	xmm7 += xmm8;			/* r3+r19				; 5-7 */ \
	xmm10 = CONST2_P809; \
	xmm10 *= xmm3;			/* .809(r9+r13)				; 5-9 */ \
\
	xmm11 = sd1[1];			/* r7 */ \
	xmm12 = sd4[2];			/* r15 */ \
	xmm1 = xmm11;			/* Copy r7 */ \
	xmm11 += xmm12;			/* r7+r15				; 6-8			avail 4 */ \
	xmm4 = CONST2_P809; \
	xmm4 *= xmm0;			/* .809(r5+r17)				; 6-10			avail none */ \
\
	xmm9 -= xmm8;			/* r3-r19				; 7-9 */ \
	xmm8 = CONST2_P309; \
	xmm8 *= xmm3;			/* .309(r9+r13)				; 7-11 */ \
\
	xmm1 -= xmm12;			/* r7-r15				; 8-10 */ \
	xmm12 = CONST2_P809; \
	xmm12 *= xmm7;			/* .809(r3+r19)				; 8-12 */ \
\
	xmm0 += xmm3;			/* (r5+r17)+(r9+r13)			; 9-11 */ \
	xmm3 = CONST2_P309; \
	xmm3 *= xmm11;			/* .309(r7+r15)				; 9-13 */ \
\
	xmm6 -= xmm10;			/* .309(r5+r17)-.809(r9+r13)		; 10-12 */ \
	xmm10 = CONST2_P309; \
	xmm10 *= xmm7;			/* .309(r3+r19)				; 10-14 */ \
\
	xmm7 += xmm11;			/* (r3+r19)+(r7+r15)			; 11-13 */ \
	xmm11 *= CONST2_P809;			/* .809(r7+r15)				; 11-15 */ \
\
	xmm8 -= xmm4;			/* -.809(r5+r17)+.309(r9+r13)		; 12-14			avail 4 */ \
\
	xmm4 = sd0[0];			/* r1 */ \
	xmm0 += xmm4;			/* r1+(r5+r17)+(r9+r13)			; 13-15			avail none */ \
\
	xmm12 -= xmm3;			/* .809(r3+r19)-.309(r7+r15)		; 14-16			avail 3 */ \
\
	xmm3 = sd0[2];			/* r11 */ \
	xmm7 += xmm3;			/* (r3+r19)+(r7+r15)+r11		; 15-17			avail none */ \
\
	xmm10 -= xmm11;			/* .309(r3+r19)-.809(r7+r15)		; 16-18 */ \
	xmm11 = CONST2_P951; \
	xmm11 *= xmm2;			/* .951(r5-r17)				; 16-20 (12) */ \
\
	xmm6 += xmm4;			/* r1+.309(r5+r17)-.809(r9+r13)		; 17-19 */ \
	xmm2 *= CONST2_P588;			/* .588(r5-r17)				; 17-21 (13) */ \
\
	xmm8 += xmm4;			/* r1-.809(r5+r17)+.309(r9+r13)		; 18-20 */ \
	xmm4 = CONST2_P588; \
	xmm4 *= xmm5;			/* .588(r9-r13)				; 18-22 (14) */ \
\
	xmm12 -= xmm3;			/* .809(r3+r19)-.309(r7+r15)-r11	; 19-21 */ \
	xmm5 *= CONST2_P951;			/* .951(r9-r13)				; 19-23 (15) */ \
\
	xmm10 += xmm3;			/* .309(r3+r19)-.809(r7+r15)+r11	; 20-22			avail 3 */ \
\
	xmm3 = xmm0;			/* Copy 1-mod-4-cols row #1 */ \
	xmm0 -= xmm7;			/* Real odd-cols row #6 (final real #6)	; 21-23	(18)		avail none storable 0 */ \
	xmm7 += xmm3;			/* Real odd-cols row #1 (final real #1A); 22-24	(19)		avail 3 storable 0,7 */ \
	xmm3 = CONST2_P588; \
	xmm3 *= xmm9;			/* .588(r3-r19)				; 22-26 (16)		avail none storable 0,7 */ \
\
	xmm11 += xmm4;			/* .951(r5-r17)+.588(r9-r13)		; 23-25 */ \
	xmm4 = CONST2_P951; \
	xmm4 *= xmm1;			/* .951(r7-r15)				; 23-27 (18) */ \
\
	xmm2 -= xmm5;			/* .588(r5-r17)-.951(r9-r13)		; 24-26			avail 5 storable 0,7 */ \
	xmm9 *= CONST2_P951;	/* .951(r3-r19)				; 24-28 (17) */ \
	t4 =xmm0;				/* Real #6				; 24			avail 5,0 storable 7 */ \
\
	xmm5 = xmm6;			/* Copy 1-mod-4-cols row #2 */ \
	xmm6 -= xmm12;			/* Real odd-cols row #5			; 25-27	(22)		avail 0 storable 7,6 */ \
	xmm1 *= CONST2_P588;	/* .588(r7-r15)				; 25-29 (19) */ \
	sd0[0] = xmm7;			/* Final real #1A			; 25			avail 0,7 storable 6 */ \
\
	xmm12 += xmm5;			/* Real odd-cols row #2			; 26-28			avail 5,0,7 storable 6,12 */ \
	xmm0 = sd1[0];			/* r2 */ \
\
	xmm5 = xmm8;			/* Copy 1-mod-4-cols row #3 */ \
	xmm8 -= xmm10;			/* Real odd-cols row #4			; 27-29			avail 7 storable 6,12,8 */ \
	xmm7 = sd4[3];			/* r20 */ \
\
	xmm10 += xmm5;			/* Real odd-cols row #3			; 28-30			avail 5 storable 6,12,8,10 */ \
	t3 = xmm6;				/* Real odd-cols row #5			; 28			avail 5,6 storable 12,8,10 */ \
	xmm6 = sd4[1];			/* r10 */ \
\
	xmm3 += xmm4;			/* .588(r3-r19)+.951(r7-r15)		; 29-31			avail 5,4 storable 12,8,10 */ \
	t0 = xmm12;				/* Real odd-cols row #2			; 29			avail 5,4,12 storable 8,10 */ \
	xmm4 = sd1[2];			/* r12 */ \
\
	xmm9 -= xmm1;			/* .951(r3-r19)-.588(r7-r15)		; 30-32			avail 5,1,12 storable 8,10 */ \
	t2 = xmm8;				/* Real odd-cols row #4			; 30			avail 5,1,12,8 storable 10 */ \
	xmm12 = sd3[0];			/* r4 */ \
\
;	/* Start the even columns (the first part will overlap with the last few odd column calculations) */ \
\
	xmm5 = xmm0;			/* Copy r2 */ \
	xmm0 += xmm7;			/* r2+r20							; 1-3	avail 1,8 storable 10 */ \
	t1 = xmm10;				/* Real odd-cols row #3			; 31			 */ \
\
	xmm1 = xmm6;			/* Copy r10 */ \
	xmm6 += xmm4;			/* r10+r12							; 2-4	avail 8,10 */ \
\
	xmm10 = xmm11;			/* Copy 1-mod-4 imag odd-cols row #2 */ \
	xmm11 += xmm3;			/* Imag odd-cols row #2			; 33-35 (32) */ \
	xmm8 = sd2[3];			/* r18 */ \
\
	xmm3 -= xmm10;			/* Imag odd-cols row #5			; 34-36 */ \
\
	xmm10 = xmm2;			/* Copy 1-mod-4 imag odd-cols row #4 */ \
	xmm2 += xmm9;			/* Imag odd-cols row #3			; 35-37 (33) */ \
\
	xmm9 -= xmm10;			/* Imag odd-cols row #4			; 36-38 */ \
	t14 = xmm11;			/* Imag odd-cols row #2			; 36				avail 10,11 */ \
	xmm11 = sd2[1];			/* r8 */ \
\
	xmm10 = xmm12;			/* Copy r4 */ \
	xmm12 += xmm8;			/* r4+r18							; 3-5 */ \
	t12 = xmm3;				/* Imag odd-cols row #5			; 37				avail 3 */ \
	xmm3 = sd3[2];			/* r14 */ \
\
	t15 = xmm2;				/* Imag odd-cols row #3			; 38 */ \
	xmm2 = xmm11;			/* Copy r8 */ \
	xmm11 += xmm3;			/* r8+r14							; 4-6	avail none */ \
\
	t10 = xmm9;				/* Imag odd-cols row #4			; 39 */ \
	xmm9 = xmm0;			/* Copy r2+r20 */ \
	xmm0 -= xmm6;			/* (r2+r20)-(r10+r12)						; 5-7 */ \
\
	xmm6 += xmm9;			/* (r2+r20)+(r10+r12)						; 6-8 */ \
\
	xmm9 = xmm12;			/* Copy r4+r18 */ \
	xmm12 -= xmm11;			/* (r4+r18)-(r8+r14)						; 7-9 */ \
\
	xmm11 += xmm9;			/* (r4+r18)+(r8+r14)						; 8-10 */ \
	xmm9 = CONST2_P951; \
	xmm9 *= xmm0;			/* .951((r2+r20)-(r10+r12))					; 8-12 */ \
\
	xmm5 -= xmm7;			/* r2-r20							; 9-11	avail 7 */ \
	xmm0 *= CONST2_P588;			/* .588((r2+r20)-(r10+r12))					; 9-13 */ \
\
	xmm1 -= xmm4;			/* r10-r12							; 10-12 */ \
	xmm4 = CONST2_P588; \
	xmm4 *= xmm12;			/* .588((r4+r18)-(r8+r14))					; 10-14 */ \
\
	xmm10 -= xmm8;			/* r4-r18							; 11-13	avail 7,8 */ \
	xmm12 *= CONST2_P951;			/* .951((r4+r18)-(r8+r14))					; 11-15 */ \
\
	xmm2 -= xmm3;			/* r8-r14							; 12-14 */ \
	xmm7 = CONST2_P809; \
	xmm7 *= xmm6;			/* .809((r2+r20)+(r10+r12))					; 12-16	avail 8,3 */ \
\
	xmm8 = xmm5;			/* Copy r2-r20 */ \
	xmm5 += xmm1;			/* (r2-r20)+(r10-r12)						; 13-15 */ \
	xmm3 = CONST2_P309; \
	xmm3 *= xmm11;			/* .309((r4+r18)+(r8+r14))					; 13-17	avail none */ \
\
	xmm8 -= xmm1;			/* (r2-r20)-(r10-r12)						; 14-16 */ \
	xmm1 = CONST2_P309; \
	xmm1 *= xmm6;			/* .309((r2+r20)+(r10+r12))					; 14-18 */ \
\
	xmm6 += xmm11;			/* ((r2+r20)+(r10+r12))+((r4+r18)+(r8+r14))			; 15-17 */ \
	xmm11 *= CONST2_P809;			/* .809((r4+r18)+(r8+r14))					; 15-19 */ \
\
	xmm9 += xmm4;			/* .951((r2+r20)-(r10+r12))+.588((r4+r18)-(r8+r14))		; 16-18 */ \
	xmm4 = CONST2_P309; \
	xmm4 *= xmm5;			/* .309((r2-r20)+(r10-r12))					; 16-20	avail none storable 9 */ \
\
	xmm0 -= xmm12;			/* .588((r2+r20)-(r10+r12))-.951((r4+r18)-(r8+r14))		; 17-19 */ \
	xmm12 = CONST2_P809; \
	xmm12 *= xmm5;			/* .809((r2-r20)+(r10-r12))					; 17-21	avail none storable 9,0 */ \
\
	t5 = xmm9;				/* Save real even-cols row #2					; 19 */ \
	xmm9 = xmm10;			/* Copy r4-r18 */ \
	xmm10 += xmm2;			/* (r4-r18)+(r8-r14)						; 18-20	avail none storable 0 */ \
\
	xmm9 -= xmm2;			/* (r4-r18)-(r8-r14)						; 19-21	avail 2 storable 0 */ \
	xmm2 = sd0[1];			/* r6 */ \
\
	xmm7 -= xmm3;			/* .809((r2+r20)+(r10+r12))-.309((r4+r18)+(r8+r14))		; 20-22	avail 3 storable 0 */ \
	xmm3 = sd0[3];			/* r16 */ \
	t6 = xmm0;				/* Save real even-cols row #4					; 20	avail 0 */ \
\
	xmm0 = xmm2;			/* Copy r6 */ \
	xmm2 += xmm3;			/* r6+r16							; 21-23 avail none */ \
\
	xmm0 -= xmm3;			/* r6-r16							; 22-24	avail 3 */ \
	xmm3 = CONST2_P809; \
	xmm3 *= xmm10;			/* .809((r4-r18)+(r8-r14))					; 22-26 (21) avail none */ \
\
	xmm1 -= xmm11;			/* .309((r2+r20)+(r10+r12))-.809((r4+r18)+(r8+r14))		; 23-25 */ \
	xmm11 = CONST2_P309; \
	xmm11 *= xmm10;			/* .309((r4-r18)+(r8-r14))					; 23-27 */ \
\
	xmm6 += xmm2;			/* ((r2+r20)+(r10+r12))+((r4+r18)+(r8+r14))+(r6+r16)		; 24-26	storable 6 */ \
	xmm7 -= xmm2;			/* .809((r2+r20)+(r10+r12))-.309((r4+r18)+(r8+r14))-(r6+r16)	; 25-27	storable 6,7 */ \
	xmm1 += xmm2;			/* .309((r2+r20)+(r10+r12))-.809((r4+r18)+(r8+r14))+(r6+r16)	; 26-28	avail 2 storable 6,7,1 */ \
\
	xmm4 += xmm3;			/* .309((r2-r20)+(r10-r12))+.809((r4-r18)+(r8-r14))		; 27-29 */ \
	xmm3 = CONST2_P588; \
	xmm3 *= xmm8;			/* .588((r2-r20)-(r10-r12))					; 27-31 (18) */ \
	sd0[1] = xmm6;			/* Save final real #1B (real even-cols row #1)			; 27	avail 2,6 storable 7,1 */ \
	xmm6 = t0;				/* Real odd-cols row #2 */ \
\
	xmm12 += xmm11;			/* .809((r2-r20)+(r10-r12))+.309((r4-r18)+(r8-r14))		; 28-30	avail 2,11 storable 7,1 */ \
	xmm8 *= CONST2_P951;	/* .951((r2-r20)-(r10-r12))					; 28-32	(19) */ \
	xmm11 = t5;				/* Real even-cols row #2 */ \
	t7 = xmm7;				/* Save real even-cols row #3					; 28	avail 2,7 storable 1 */ \
\
	xmm5 -= xmm10;			/* ((r2-r20)+(r10-r12))-((r4-r18)+(r8-r14))			; 29-31	avail 2,7,10 storable 1 */ \
	xmm2 = CONST2_P951; \
	xmm2 *= xmm9;			/* .951((r4-r18)-(r8-r14))					; 29-33 (22) */ \
	xmm7 = t14;				/* Imag odd-cols row #2 */ \
	t8 = xmm1;				/* Save real even-cols row #5					; 29	avail 10,1 */ \
\
	xmm4 += xmm0;			/* .309((r2-r20)+(r10-r12))+.809((r4-r18)+(r8-r14))+(r6-r16)	; 30-32	avail 10,1 storable 4 */ \
	xmm9 *= CONST2_P588;	/* .588((r4-r18)-(r8-r14))					; 30-34 */ \
\
;	/* Now combine the even and odd columns then do the post-multiply by twiddle factors (start as we finish up even cols) */ \
\
	xmm10 = xmm6;			/* Copy real odd-cols row #2 */ \
	xmm6 -= xmm11;			/* Real #10				; 1-3 */ \
	xmm11 += xmm10;			/* Real #2				; 2-4 */ \
	xmm1 = sc[8*2+1];		/* cosine/sine */ \
	xmm10 = xmm4;			/* Copy imag even-cols row #2 */ \
	xmm4 -= xmm7;			/* Imag #10				; 3-5 */ \
	xmm7 += xmm10;			/* Imag #2				; 4-6	avail 10 */ \
	xmm10 = t1;				/* Real odd-cols row #3 */ \
\
	xmm3 += xmm2;			/* .588((r2-r20)-(r10-r12))+.951((r4-r18)-(r8-r14))		; 35-37 avail 2 storable 3 */ \
	xmm2 = t7;				/* Real even-cols row #3 */ \
	xmm8 -= xmm9;			/* .951((r2-r20)-(r10-r12))-.588((r4-r18)-(r8-r14))		; 36-38 avail 9 storable 3,8 */ \
	xmm9 = t15;				/* Imag odd-cols row #3 */ \
	xmm12 -= xmm0;			/* .809((r2-r20)+(r10-r12))+.309((r4-r18)+(r8-r14))-(r6-r16)	; 37-39 avail none storable 3,8,12 */ \
	xmm5 += xmm0;			/* ((r2-r20)+(r10-r12))-((r4-r18)+(r8-r14))+(r6-r16)		; 38-40 avail 0 storable 3,8,12,5 */ \
\
	xmm0 = xmm6;			/* Copy R10 */ \
	xmm6 *= xmm1;			/* A10 = R10 * cosine/sine		; 8-12 (4) avail none storable 3,8,12,5 */ \
\
	t13 = xmm8;				/* Save imag even-cols row #5					; 39 */ \
	xmm8 = xmm10;			/* Copy real odd-cols row #3 */ \
	xmm10 -= xmm2;			/* Real #9				; 9-11 */ \
	xmm1 *= xmm4;			/* B10 = I10 * cosine/sine		; 9-13 (6) avail none storable 3,12,5 */ \
\
	xmm2 += xmm8;			/* Real #3				; 10-12 */ \
	xmm8 = sc[1];			/* cosine/sine */ \
	t11 = xmm12;			/* Save imag even-cols row #4					; 40 */ \
	xmm12 = xmm11;			/* Copy R2 */ \
	xmm11 *= xmm8;			/* A2 = R2 * cosine/sine		; 10-14 (5) avail none storable 3,5 */ \
\
	t9 = xmm5;				/* Save imag row #6						; 41 */ \
	xmm5 = xmm3;		 	/* Copy imag even-cols row #3 */ \
	xmm3 -= xmm9;			/* Imag #9				; 11-13 */ \
	xmm8 *= xmm7;			/* B2 = I2 * cosine/sine		; 11-15	avail none */ \
\
	xmm9 += xmm5;			/* Imag #3				; 12-14	avail 5 */ \
	xmm5 = sc[7*2+1];		/* cosine/sine */ \
	xmm13 = xmm10;			/* Copy R9 */ \
	xmm10 *= xmm5;			/* A9 = R9 * cosine/sine		; 12-16	avail none */ \
\
	xmm6 -= xmm4;			/* A10 = A10 - I10			; 13-15	avail 4 */ \
	xmm4 = sc[2+1];			/* cosine/sine */ \
	xmm14 = xmm2;			/* Copy R3 */ \
	xmm2 *= xmm4;			/* A3 = R3 * cosine/sine		; 13-17 avail none */ \
\
	xmm1 += xmm0;			/* B10 = B10 + R10			; 14-16	avail 0 */ \
	xmm5 *= xmm3;			/* B9 = I9 * cosine/sine		; 14-18 */ \
	xmm0 = sc[8*2];			/* sine */ \
\
	xmm11 -= xmm7;			/* A2 = A2 - I2				; 15-17	avail 7 */ \
	xmm4 *= xmm9;			/* B3 = I3 * cosine/sine		; 15-19 */ \
	xmm7 = sc[0];			/* sine */ \
\
	xmm8 += xmm12;			/* B2 = B2 + R2				; 16-18 */ \
	xmm12 = t2;				/* Real odd-cols row #4 */ \
	xmm6 *= xmm0;			/* A10 = A10 * sine (final R10)		; 16-20	avail none storable 6 */ \
\
	xmm10 -= xmm3;			/* A9 = A9 - I9				; 17-19	avail 3 storable 6 */ \
	xmm1 *= xmm0;			/* B10 = B10 * sine (final I10)		; 17-21	avail 3,0 storable 6,1 */ \
	xmm3 = t6;				/* Real even-cols row #4 */ \
\
	xmm2 -= xmm9;			/* A3 = A3 - I3				; 18-20	avail 0,9 storable 6,1 */ \
	xmm11 *= xmm7;			/* A2 = A2 * sine (final R2)		; 18-22	avail 0,9 storable 6,1,11 */ \
	xmm0 = sc[7*2];			/* sine */ \
\
	xmm5 += xmm13;			/* B9 = B9 + R9				; 19-21	avail 9,13 storable 6,1,11 */ \
	xmm8 *= xmm7;			/* B2 = B2 * sine (final I2)		; 19-23	avail 9,13,7 storable 6,1,11,8 */ \
	xmm7 = sc[2];			/* sine */ \
\
	xmm4 += xmm14;			/* B3 = B3 + R3				; 20-22	avail 9,13,14 storable 6,1,11,8 */ \
	xmm10 *= xmm0;			/* A9 = A9 * sine (final R9)		; 20-24	avail 9,13,14 storable 6,1,11,8,10 */ \
	xmm9 = t11;				/* Imag even-cols row #4 */ \
\
	xmm13 = xmm12;			/* Copy real odd-cols row #4 */ \
	xmm12 -= xmm3;			/* Real #8				; 21-23	avail 14 storable 6,1,11,8,10 */ \
	xmm2 *= xmm7;			/* A3 = A3 * sine (final R3)		; 21-25 */ \
	xmm14 = t10;			/* Imag odd-cols row #4 */ \
	sd4[2] = xmm6;			/* Save final R10			; 21	avail 6 storable 1,11,8,10,2 */ \
\
	xmm6 = xmm9;			/* Copy imag even-cols row #4 */ \
	xmm9 -= xmm14;			/* Imag #8				; 22-24 */ \
	xmm5 *= xmm0;			/* B9 = B9 * sine (final I9)		; 22-26	avail 0 storable 1,11,8,10,2,5 */ \
	xmm0 = sc[6*2+1];		/* cosine/sine */ \
	sd4[3] = xmm1;			/* Save final I10			; 22	avail 1 storable 11,8,10,2,5 */ \
\
	xmm3 += xmm13;			/* Real #4				; 23-25	avail 1,13 storable 11,8,10,2,5 */ \
	xmm4 *= xmm7;			/* B3 = B3 * sine (final I3)		; 23-27	avail 1,13,7 storable 11,8,10,2,5,4 */ \
	xmm1 = t3;				/* Real odd-cols row #5 */ \
	sd0[2] = xmm11;			/* Save final R2			; 23 */ \
\
	xmm14 += xmm6;			/* Imag #4				; 24-26	avail 13,7,11,6 storable 8,10,2,5,4 */ \
	xmm13 = xmm12;			/* Copy R8 */ \
	xmm12 *= xmm0;			/* A8 = R8 * cosine/sine		; 24-28	avail 7,11,6 storable 8,10,2,5,4 */ \
	xmm7 = t8;				/* Real even-cols row #5 */ \
	sd0[3] = xmm8;			/* Save final I2			; 24	avail 11,6,8 storable 10,2,5,4 */ \
\
	xmm8 = xmm1;			/* Copy real odd-cols row #5 */ \
	xmm1 -= xmm7;			/* Real #7				; 25-27	avail 11,6 storable 10,2,5,4 */ \
	xmm0 *= xmm9;			/* B8 = I8 * cosine/sine		; 25-29 */ \
	sd4[0] = xmm10;			/* Save final R9			; 25	avail 11,6,10 storable 2,5,4 */ \
\
	xmm7 += xmm8;			/* Real #5				; 26-28	avail 11,6,10,8 storable 2,5,4 */ \
	xmm8 = sc[2*2+1];		/* cosine/sine */ \
	xmm10 = xmm3;			/* Copy R4 */ \
	xmm3 *= xmm8;			/* A4 = R4 * cosine/sine		; 26-30	avail 11,6 storable 2,5,4 */ \
	sd1[0] = xmm2;			/* Save final R3			; 26	avail 11,6,2 storable 5,4 */ \
\
	xmm11 = t13;			/* Imag even-cols row #5 */ \
	xmm6 = t12;				/* Imag odd-cols row #5 */ \
	xmm2 = xmm11;			/* Copy imag even-cols row #5 */ \
	xmm11 -= xmm6;			/* Imag #7				; 27-29	avail none storable 5,4 */ \
	xmm8 *= xmm14;			/* B4 = I4 * cosine/sine		; 27-31 */ \
	sd4[1] = xmm5;			/* Save final I9			; 27	avail 5 storable 4 */ \
\
	xmm6 += xmm2;			/* Imag #5				; 28-30	avail 5,2 storable 4 */ \
	xmm5 = t4;				/* Real #6 */ \
	xmm5 *= sc[4*2+1];		/* A6 = R6 * cosine/sine		; 28-32 */ \
	sd1[1] = xmm4;			/* Save final I3			; 28	avail 2,4 */ \
\
	xmm12 -= xmm9;			/* A8 = A8 - I8				; 29-31	avail 2,4,9 */ \
	xmm2 = sc[5*2+1];		/* cosine/sine */ \
	xmm4 = xmm1;			/* Copy R7 */ \
	xmm1 *= xmm2;			/* A7 = R7 * cosine/sine		; 29-33	avail 9 */ \
\
	xmm0 += xmm13;			/* B8 = B8 + R8				; 30-32	avail 9,13 */ \
	xmm2 *= xmm11;			/* B7 = I7 * cosine/sine		; 30-34 */ \
\
	xmm3 -= xmm14;			/* A4 = A4 - I4				; 31-33	avail 9,13,14 */ \
	xmm9 = sc[3*2+1];		/* cosine/sine */ \
	xmm13 = xmm7;			/* Copy R5 */ \
	xmm7 *= xmm9;			/* A5 = R5 * cosine/sine		; 31-35	avail 14 */ \
\
	xmm8 += xmm10;			/* B4 = B4 + R4				; 32-34	avail 14,10 */ \
	xmm9 *= xmm6;			/* B5 = I5 * cosine/sine		; 32-35 */ \
	xmm14 = t9;				/* Imag #6 */ \
\
	xmm5 -= xmm14;			/* A6 = A6 - I6				; 33-35	avail 10 */ \
	xmm14 *= sc[4*2+1];		/* B6 = I6 * cosine/sine		; 33-36 */ \
	xmm10 = sc[6*2];		/* sine */ \
\
	xmm1 -= xmm11;			/* A7 = A7 - I7				; 34-36	avail 10,11 */ \
	xmm12 *= xmm10;			/* A8 = A8 * sine (final R8)		; 34-38	avail 11 storable 12 */ \
	xmm11 = sc[2*2];		/* sine */ \
\
	xmm2 += xmm4;			/* B7 = B7 + R7				; 35-37	avail 4 storable 12 */ \
	xmm0 *= xmm10;			/* B8 = B8 * sine (final I8)		; 35-39	avail 4,10 storable 12,0 */ \
	xmm4 = sc[5*2];			/* sine */ \
\
	xmm7 -= xmm6;			/* A5 = A5 - I5				; 36-38	avail 10,6 storable 12,0 */ \
	xmm3 *= xmm11;			/* A4 = A4 * sine (final R4)		; 36-40	avail 10,6 storable 12,0,3 */ \
	xmm6 = sc[4*2];			/* sine */ \
\
	xmm9 += xmm13;			/* B5 = B5 + R5				; 37-39	avail 10,13 storable 12,0,3 */ \
	xmm8 *= xmm11;			/* B4 = B4 * sine (final I4)		; 37-41	avail 10,13,11 storable 12,0,3,8 */ \
\
	xmm14 += t4;			/* B6 = B6 + R6				; 38-40 */ \
	xmm5 *= xmm6;			/* A6 = A6 * sine (final R6)		; 38-42	avail 10,13,11 storable 12,0,3,8,5 */ \
\
	xmm1 *= xmm4;			/* A7 = A7 * sine (final R7)		; 39-43	avail 10,13,11 storable 12,0,3,8,5,1 */ \
	sd3[2] = xmm12;			/* Save final R8			; 39	avail 10,13,11,12 storable 0,3,8,5,1 */ \
\
	xmm2 *= xmm4;			/* B7 = B7 * sine (final I7)		; 40-44	avail 10,13,11,12,4 storable 0,3,8,5,1,2 */ \
	sd3[3] = xmm0;			/* Save final I8			; 40	avail 10,13,11,12,4,0 storable 3,8,5,1,2 */ \
\
	xmm10 = sc[3*2];		/* sine  */ \
	xmm7 *= xmm10;			/* A5 = A5 * sine (final R5)		; 41-45	avail 13,11,12,4,0 storable 3,8,5,1,2,7 */ \
	sd1[2] = xmm3;			/* Save final R4			; 41	avail 13,11,12,4,0,3 storable 8,5,1,2,7 */ \
\
	xmm9 *= xmm10;			/* B5 = B5 * sine (final I5)		; 42-46	avail 13,11,12,4,0,3,10 storable 8,5,1,2,7,9 */ \
	sd1[3] = xmm8;			/* Save final I4			; 42	avail 13,11,12,4,0,3,10,8 storable 5,1,2,7,9 */ \
\
	xmm14 *= xmm6;			/* B6 = B6 * sine (final I6)		; 43-47	avail 13,11,12,4,0,3,10,8,6 storable 5,1,2,7,9,14 */ \
	sd2[2] = xmm5;			/* Save final R6			; 43	avail 13,11,12,4,0,3,10,8,6,5 storable 1.2,7,9,14 */ \
\
	sd3[0] = xmm1;			/* Save final R7			; 44 */ \
	sd3[1] = xmm2;			/* Save final I7			; 45 */ \
	sd2[0] = xmm7;			/* Save final R5			; 46 */ \
	sd2[1] = xmm9;			/* Save final I5			; 47 */ \
	sd2[3] = xmm14;			/* Save final I6			; 48 */ \
	/*xmm14 = CONST2_P809;*/ \
\
	srcreg += srcinc; \
}

#endif

/* */ \
/* ************************************* 20-reals-last-unfft variants ****************************************** */ \
/* */ \
\
/* These macros produce 20 reals after doing 4.32 levels of the inverse FFT applying */ \
/* the sin/cos multipliers beforehand.  The input is 2 real and 9 complex numbers. */ \
\
/* To calculate a 20-reals inverse FFT, we calculate 20 real values from 20 complex inputs in a brute force way. */ \
/* First we note that the 20 complex values are computed from the 9 complex and 2 real inputs using Hermetian symmetry, thusly: */ \
/* c1 = r1A + 0*i */ \
/* c2 = r2 + i2*i */ \
/* ... */ \
/* c10 = r10 + i10*i */ \
/* c11 = r1B + 0*i */ \
/* c12 = r10 - i10*i */ \
/* ... */ \
/* c20 = r2 - i2*i  */ \
/* */ \
/* The brute force calculations are: */ \
/* */ \
/* c1 + c2 + ... + c20	*  w^-0000000000... */ \
/* c1 + c2 + ... + c20	*  w^-0123456789A... */ \
/* c1 + c2 + ... + c20	*  w^-02468ACE.... */ \
/*    ... */ \
/* c1 + c2 + ... + c20	*  w^-...A987654321 */ \
/* */ \
/* The sin/cos values (w = 20th root of unity) are: */ \
/* w^-1 = .951 - .309i */ \
/* w^-2 = .809 - .588i */ \
/* w^-3 = .588 - .809i */ \
/* w^-4 = .309 - .951i */ \
/* w^-5 = 0 - 1i */ \
/* w^-6 = -.309 - .951i */ \
/* w^-7 = -.588 - .809i */ \
/* w^-8 = -.809 - .588i */ \
/* w^-9 = -.951 - .309i */ \
/* w^-10 = -1 */ \
/* */ \
/* Applying the sin/cos values above, taking advantage of symmetry, and ignoring a lot of multiplies by 2: */ \
/* r1     +(r2+r10)     +(r3+r9)     +(r4+r8)     +(r5+r7) + r6 + r11 */ \
/* r1 +.951(r2-r10) +.809(r3-r9) +.588(r4-r8) +.309(r5-r7)      - r11 +.309(i2+i10) +.588*(i3+i9) +.809(i4+i8) +.951(i5+i7) + i6 */ \
/* r1 +.809(r2+r10) +.309(r3+r9) -.309(r4+r8) -.809(r5+r7) - r6 + r11 +.588(i2-i10) +.951*(i3-i9) +.951(i4-i8) +.588(i5-i7) */ \
/* r1 +.588(r2-r10) -.309(r3-r9) -.951(r4-r8) -.809(r5-r7)      - r11 +.809(i2+i10) +.951*(i3+i9) +.309(i4+i8) -.588(i5+i7) - i6 */ \
/* r1 +.309(r2+r10) -.809(r3+r9) -.809(r4+r8) +.309(r5+r7) + r6 + r11 +.951(i2-i10) +.588*(i3-i9) -.588(i4-i8) -.951(i5-i7) */ \
/* r1                   -(r3-r9)                  +(r5-r7)      - r11     +(i2+i10)                   -(i4+i8)              + i6 */ \
/* r1 -.309(r2+r10) -.809(r3+r9) +.809(r4+r8) +.309(r5+r7) - r6 + r11 +.951(i2-i10) -.588*(i3-i9) -.588(i4-i8) +.951(i5-i7) */ \
/* r1 -.588(r2-r10) -.309(r3-r9) +.951(r4-r8) -.809(r5-r7)      - r11 +.809(i2+i10) -.951*(i3+i9) +.309(i4+i8) +.588(i5+i7) - i6 */ \
/* r1 -.809(r2+r10) +.309(r3+r9) +.309(r4+r8) -.809(r5+r7) + r6 + r11 +.588(i2-i10) -.951*(i3-i9) +.951(i4-i8) -.588(i5-i7) */ \
/* r1 -.951(r2-r10) +.809(r3-r9) -.588(r4-r8) +.309(r5-r7)      - r11 +.309(i2+i10) -.588*(i3+i9) +.809(i4+i8) -.951(i5+i7) + i6 */ \
/* r1     -(r2+r10)     +(r3+r9)     -(r4+r8)     +(r5+r7) - r6 + r11 */ \
/* r1 -.951(r2-r10) +.809(r3-r9) -.588(r4-r8) +.309(r5-r7)      - r11 -.309(i2+i10) +.588*(i3+i9) -.809(i4+i8) +.951(i5+i7) */ \
/* ... r13 thru r20 are the same as r8 through r1 but with the sign of the imaginary component changed. */ \
/* */ \
/* Also remember that due to the funny way we do things reals input r1A = r1+r11 and r1B = r1-11 */ \
\
/* Store intermediate results in g->u.xmm.XMM_COL_MULTS (an 8KB buffer used in normalization) */ \
\

#define r5_x10cl_20_reals_last_unfft_preload \
\
\

#define r5_x10cl_20_reals_last_unfft(srcreg,srcinc,d1,screg,scoff) \
	vec2f64 temp1[19]; \
	vec2f64 temp2[19]; \
\
	/* Using pointers improves clang ARM64 compiler. */ \
	vec2f64 *sd0 = (vec2f64*)(srcreg); \
	vec2f64 *sd1 = (vec2f64*)(srcreg+d1*1); \
	vec2f64 *sd2 = (vec2f64*)(srcreg+d1*2); \
	vec2f64 *sd3 = (vec2f64*)(srcreg+d1*3); \
	vec2f64 *sd4 = (vec2f64*)(srcreg+d1*4); \
	vec2f64 *sd5 = (vec2f64*)(srcreg+d1*5); \
	vec2f64 *sd6 = (vec2f64*)(srcreg+d1*6); \
	vec2f64 *sd7 = (vec2f64*)(srcreg+d1*7); \
	vec2f64 *sd8 = (vec2f64*)(srcreg+d1*8); \
	vec2f64 *sd9 = (vec2f64*)(srcreg+d1*9); \
\
	vec2f64 *sc1 = (vec2f64*)(screg); \
	vec2f64 *sc2 = (vec2f64*)(screg+scoff); \
\
;	/* Apply the 9 twiddle factors */ \
\
\
	r5_x9c_twiddle_PART8(sc1,0,temp1) \
	r5_x9c_twiddle_PART(sc1,0,temp1,sd1,sd9,0,16) \
	r5_x9c_twiddle_PART(sc1,0,temp1,sd3,sd7,4,12) \
	r5_x9c_twiddle_PART(sc2,1,temp2,sd1,sd9,0,16) \
	r5_x9c_twiddle_PART(sc2,1,temp2,sd3,sd7,4,12) \
	r5_x9c_twiddle_PART(sc2,1,temp2,sd4,sd6,6,10) \
	r5_x9c_twiddle_PART(sc2,1,temp2,sd2,sd8,2,14) \
	r5_x9c_twiddle_PART8(sc2,1,temp2) \
\
;	/* Do the 20 reals inverse FFT */ \
\
	r5_x20r_unfft(sd1,sd3,sd5,sd7,sd9,sd0[1],sd0[3],temp2); \
	r5_x9c_twiddle_PART(sc1,0,temp1,sd4,sd6,6,10) \
	r5_x9c_twiddle_PART(sc1,0,temp1,sd2,sd8,2,14) \
	r5_x20r_unfft(sd0,sd2,sd4,sd6,sd8,sd0[0],sd0[2],temp1); \
\
	srcreg += srcinc; \
\
\
/* Apply the 9 twiddle factors.  Also do the first +/- on the reals so that */ \
/* we have a perfect balance of adds and subtracts. */ \

#define r5_x9c_twiddle_PART(SC,SN,tmpreg,SD1,SD2,SCN1,SCN2) \
	xmm0 = SC[SCN1+1];		/* cosine/sine */ \
	xmm4 = SC[SCN2+1];		/* cosine/sine */ \
	xmm1 = SD1[0+SN];		/* R5 */ \
	xmm5 = SD2[0+SN];		/* R7 */ \
	xmm3 = SD1[2+SN];		/* I5 */ \
	xmm7 = SD2[2+SN];		/* I7 */ \
	xmm2 = xmm1 * xmm0 + xmm3;	/* A5 = R5 * cosine/sine + I5 */ \
	xmm6 = xmm5 * xmm4 + xmm7;	/* A7 = R7 * cosine/sine + I7 */ \
	xmm3 = xmm3 * xmm0 - xmm1;	/* B5 = I5 * cosine/sine - R5 */ \
	xmm7 = xmm7 * xmm4 - xmm5;	/* B7 = I7 * cosine/sine - R7 */ \
	xmm1 = SC[SCN1];		/* sine */ \
	xmm5 = SC[SCN2];		/* sine */ \
	xmm2 *= xmm1;			/* R5 = A5 * sine */ \
	xmm6 *= xmm5;			/* R7 = A7 * sine */ \
	xmm3 *= xmm1;			/* I5 = B5 * sine */ \
	xmm7 *= xmm5;			/* I7 = B7 * sine */ \
	xmm1 = xmm2;			/* Copy R5 */ \
	xmm2 += xmm6;			/* R5+R7 */ \
	xmm1 -= xmm6;			/* R5-R7 */ \
	xmm5 = xmm3;			/* Copy I5 */ \
	xmm3 += xmm7;			/* I5+I7 */ \
	xmm5 -= xmm7;			/* I5-I7 */ \
	tmpreg[SCN2] = xmm2;	/* Save R5+R7 */ \
	tmpreg[SCN1] = xmm1;	/* Save R5-R7 */ \
	tmpreg[SCN2+1] = xmm3;	/* Save I5+I7 */ \
	tmpreg[SCN1+1] = xmm5;	/* Save I5-I7 */ \

#define r5_x9c_twiddle_PART8(SC,SN,tmpreg) \
	xmm4 = sd5[0+SN];		/* R6 */ \
	xmm5 = sd5[2+SN];		/* I6 */ \
	xmm0 = SC[4*2+1];		/* cosine/sine */ \
	xmm1 = SC[4*2];			/* sine */ \
	xmm2 = xmm4 * xmm0 + xmm5;	/* A6 = R6 * cosine/sine + I6 */ \
	xmm3 = xmm5 * xmm0 - xmm4;	/* B6 = I6 * cosine/sine - R6 */ \
	xmm2 *= xmm1;			/* R6 = A6 * sine */ \
	xmm3 *= xmm1;			/* I6 = B6 * sine */ \
	tmpreg[8] = xmm2;		/* Save R6 */ \
	tmpreg[9] = xmm3;		/* Save I6 */ \



#define r5_x20r_unfft(sd0,sd2,sd4,sd6,sd8,r1A,r1B,tmpreg) \
\
;	/* Calculate odd columns derived from real inputs (even rows) */ \
\
	xmm5 = tmpreg[2];			/* r3-r9 */ \
	xmm4 = r1B;					/* r1-r11 */ \
	xmm6 = xmm5 * CONST2_P809 + xmm4;	/* r1+.809(r3-r9)-r11 */ \
	xmm7 = xmm4 - xmm5 * CONST2_P309;	/* r1-.309(r3-r9)-r11 */ \
	xmm0 = xmm4 - xmm5;			/* r1-(r3-r9)-r11 */ \
\
	xmm1 = tmpreg[6];			/* r5-r7 */ \
	xmm6 += xmm1 * CONST2_P309;	/* r1+.809(r3-r9)+.309(r5-r7)-r11 */ \
	xmm7 -= xmm1 * CONST2_P809;	/* r1-.309(r3-r9)-.809(r5-r7)-r11 */ \
	xmm1 += xmm0;				/* r1-(r3-r9)+(r5-r7)-r11 */ \
\
	tmpreg[18] = xmm1;			/* Save odd-real-cols row #6 (also is real-cols row #6) */ \
\
;	/* Calculate even columns derived from real inputs (even rows) */ \
;	/* From above, odd-real-col rols rows #2,4 are in xmm6, xmm7 */ \
\
	xmm2 = tmpreg[0];			/* r2-r10 */ \
	xmm0 = tmpreg[4];			/* r4-r8 */ \
	xmm1 = xmm2 * CONST2_P951;	/* .951(r2-r10) */ \
	xmm2 = xmm2 * CONST2_P588;	/* .588(r2-r10) */ \
	xmm1 += xmm0 * CONST2_P588;	/* .951(r2-r10)+.588(r4-r8) */ \
	xmm2 -= xmm0 * CONST2_P951;	/* .588(r2-r10)-.951(r4-r8) */ \
\
;	/* Combine even and odd columns (even rows) */ \
\
	xmm0 = xmm6;			/* Copy odd-real-cols row #2 */ \
	xmm6 -= xmm1;			/* real-cols row #10 (odd#2 - even#2) */ \
	xmm1 += xmm0;			/* real-cols row #2 (odd#2 + even#2) */ \
\
	xmm4 = xmm7;			/* Copy odd-real-cols row #4 */ \
	xmm7 -= xmm2;			/* real-cols row #8 (odd#4 - even#4) */ \
	xmm2 += xmm4;			/* real-cols row #4 (odd#4 + even#4) */ \
\
	tmpreg[2] = xmm6;		/* Save real-cols row #10 */ \
	tmpreg[0] = xmm1;			/* Save real-cols row #2 */ \
	tmpreg[6] = xmm7;		/* Save real-cols row #8 */ \
	tmpreg[4] = xmm2;		/* Save real-cols row #4 */ \
\
;	/* Calculate even columns derived from real inputs (odd rows) */ \
\
	xmm4 = tmpreg[8];			/* r6 */ \
	xmm5 = tmpreg[16];			/* r2+r10 */ \
\
	xmm0 = xmm5 + xmm4;			/* (r2+10)+r6 */ \
	xmm6 = xmm5 * CONST2_P809 - xmm4;	/* .809(r2+r10)-r6 */ \
	xmm7 = xmm5 * CONST2_P309 + xmm4;	/* .309(r2+r10)+r6 */ \
\
	xmm1 = tmpreg[12];			/* r4+r8 */ \
	xmm0 += xmm1;				/* (r2+r10)+(r4+r8)+r6 */ \
	xmm6 -= xmm1 * CONST2_P309;	/* .809(r2+r10)-.309(r4+r8)-r6 */ \
	xmm7 -= xmm1 * CONST2_P809;	/* .309(r2+r10)-.809(r4+r8)+r6 */ \
\
;	/* Calculate odd columns derived from real inputs (odd rows) */ \
;	/* From above, even-real-cols row #1,3,5 are in xmm0,xmm6,xmm7 */ \
\
	xmm5 = tmpreg[14];			/* r3+r9 */ \
	xmm3 = r1A;					/* r1+r11 */ \
	xmm1 = xmm5 + xmm3;			/* r1+(r3+r9)+r11 */ \
	xmm2 = xmm5 * CONST2_P309 + xmm3;	/* r1+.309(r3+r9)+r11 */ \
	xmm3 -= xmm5 * CONST2_P809;	/* r1-.809(r3+r9)+r11 */ \
\
	xmm4 = tmpreg[10];			/* r5+r7 */ \
	xmm1 += xmm4;				/* r1+(r3+r9)+(r5+r7)+r11 */ \
	xmm2 -= xmm4 * CONST2_P809;	/* r1+.309(r3+r9)-.809(r5+r7)+r11 */ \
	xmm3 += xmm4 * CONST2_P309;	/* r1-.809(r3+r9)+.309(r5+r7)+r11 */ \
\
;	/* Combine even and odd columns (odd rows) */ \
\
	xmm4 = xmm0;			/* Copy even-real-cols row #1 */ \
	xmm0 += xmm1;			/* real-cols row #1 (and final R1) */ \
	xmm1 -= xmm4;			/* real-cols row #11 (and final R11) */ \
	sd0[0] = xmm0;			/* Save final R1 */ \
	sd0[2] = xmm1;			/* Save final R11 */ \
\
	xmm5 = xmm6;			/* Copy even-real-cols row #3 */ \
	xmm6 += xmm2;			/* real-cols row #3 */ \
	xmm2 -= xmm5;			/* real-cols row #9 */ \
	tmpreg[8] = xmm6;		/* Save real-cols row #3 */ \
	tmpreg[10] = xmm2;		/* Save real-cols row #9 */ \
\
	xmm4 = xmm7;			/* Copy even-real-cols row #5 */ \
	xmm7 += xmm3;			/* real-cols row #5 */ \
	xmm3 -= xmm4;	 		/* real-cols row #7 */ \
	tmpreg[12] = xmm7;		/* Save real-cols row #5 */ \
	tmpreg[14] = xmm3;		/* Save real-cols row #7 */ \
\
;	/* Calculate even columns derived from imaginary inputs (even rows) */ \
\
	xmm3 = tmpreg[17];			/* i2+i10 */ \
	xmm4 = tmpreg[9];			/* i6 */ \
	xmm1 = tmpreg[13];			/* i4+i8 */ \
	xmm0 = xmm3 + xmm4;			/* (i2+i10)+i6 */ \
	xmm6 = xmm3 * CONST2_P309 + xmm4;	/* .309(i2+i10)+i6 */ \
	xmm7 = xmm3 * CONST2_P809 - xmm4;	/* .809(i2+i10)-i6 */ \
	xmm0 -= xmm1;				/* (i2+i10)-(i4+i8)+i6 */ \
	xmm6 += xmm1 * CONST2_P809;	/* .309(i2+i10)+.809(i4+i8)+i6 */ \
	xmm7 += xmm1 * CONST2_P309;	/* .809(i2+i10)+.309(i4+i8)-i6 */ \
\
;	/* Combine real and imaginary data for row #6 */ \
\
	xmm4 = tmpreg[18];		/* Load real-cols row #6 */ \
	xmm4 -= xmm0;			/* final R16 */ \
	xmm0 += tmpreg[18];		/* final R6 */ \
	sd0[3] = xmm4;		/* Save R16 */ \
	sd0[1] = xmm0;		/* Save R6 */ \
\
;	/* Calculate odd columns derived from imaginary inputs (even rows) */ \
;	/* From above, even-imag-cols row #2,4 are in xmm6, xmm7 */ \
\
	xmm4 = tmpreg[15];			/* i3+i9 */ \
	xmm0 = tmpreg[11];			/* i5+i7 */ \
	xmm3 = xmm4 * CONST2_P588;	/* .588(i3+i9) */ \
	xmm4 = xmm4 * CONST2_P951;	/* .951(i3+i9) */ \
	xmm3 += xmm0 * CONST2_P951;	/* .588(i3+i9)+.951(i5+i7) */ \
	xmm4 -= xmm0 * CONST2_P588;	/* .951(i3+i9)-.588(i5+i7) */ \
\
;	/* Combine even and odd columns, then real and imag data (even rows) */ \
\
	xmm0 = xmm6;			/* Copy even-imag-cols row #2 */ \
	xmm6 -= xmm3;			/* imag-cols row #10 (even#2 - odd#2) */ \
	xmm3 += xmm0;			/* imag-cols row #2 (even#2 + odd#2) */ \
	xmm1 = xmm7;			/* Copy even-imag-cols row #4 */ \
	xmm7 -= xmm4;			/* imag-cols row #8 (even#4 - odd#4) */ \
	xmm4 += xmm1;			/* imag-cols row #4 (even#4 + odd#4) */ \
\
	xmm0 = tmpreg[2];		/* Load real-cols row #10 */ \
	xmm0 -= xmm6;			/* final R12 */ \
	xmm6 += tmpreg[2];		/* final R10 */ \
	xmm1 = tmpreg[0];		/* Load real-cols row #2 */ \
	xmm1 -= xmm3;			/* final R20 */ \
	xmm3 += tmpreg[0];		/* final R2 */ \
	sd2[2] = xmm0;			/* Save R12 */ \
	sd8[1] = xmm6;			/* Save R10 */ \
	sd8[3] = xmm1;			/* Save R20 */ \
	sd2[0] = xmm3;			/* Save R2 */ \
\
	xmm0 = tmpreg[6];		/* Load real-cols row #8 */ \
	xmm0 -= xmm7;			/* final R14 */ \
	xmm7 += tmpreg[6];		/* final R8 */ \
	xmm1 = tmpreg[4];		/* Load real-cols row #4 */ \
	xmm1 -= xmm4;			/* final R18 */ \
	xmm4 += tmpreg[4];		/* final R4 */ \
	sd6[2] = xmm0;			/* Save R14 */ \
	sd4[1] = xmm7;			/* Save R8 */ \
	sd4[3] = xmm1;			/* Save R18 */ \
	sd6[0] = xmm4;			/* Save R4 */ \
\
;	/* Calculate even columns derived from imaginary inputs (odd rows) */ \
\
	xmm7 = tmpreg[1];			/* i2-i10 */ \
	xmm4 = tmpreg[5];			/* i4-i8 */ \
	xmm6 = xmm7 * CONST2_P588;	/* .588(i2-i10) */ \
	xmm7 = xmm7 * CONST2_P951;	/* .951(i2-i10) */ \
	xmm6 += xmm4 * CONST2_P951;	/* .588(i2-i10)+.951(i4-i8) */ \
	xmm7 -= xmm4 * CONST2_P588;	/* .951(i2-i10)-.588(i4-i8) */ \
\
;	/* Calculate odd columns derived from imaginary inputs (odd rows) */ \
;	/* From above, even-imag-cols row #3,5 are in xmm6,xmm7 */ \
\
	xmm4 = tmpreg[3];			/* i3-i9 */ \
	xmm0 = tmpreg[7];			/* i5-i7 */ \
	xmm3 = xmm4 * CONST2_P951;	/* .951(i3-i9) */ \
	xmm4 = xmm4 * CONST2_P588;	/* .588(i3-i9) */ \
	xmm3 += xmm0 * CONST2_P588;	/* .951(i3-i9)+.588(i5-i7) */ \
	xmm4 -= xmm0 * CONST2_P951;	/* .588(i3-i9)-.951(i5-i7) */ \
\
;	/* Combine even and odd columns, then real and imag data (odd rows) */ \
\
	xmm0 = xmm6;			/* Copy even-imag-cols row #3 */ \
	xmm6 -= xmm3;			/* imag-cols row #9 (even#3 - odd#3) */ \
	xmm3 += xmm0;			/* imag-cols row #3 (even#3 + odd#3) */ \
	xmm1 = xmm7;			/* Copy even-imag-cols row #5 */ \
	xmm7 -= xmm4;			/* imag-cols row #7 (even#5 - odd#5) */ \
	xmm4 += xmm1;			/* imag-cols row #5 (even#5 + odd#5) */ \
\
	xmm0 = tmpreg[10];		/* Load real-cols row #9 */ \
	xmm0 -= xmm6;			/* final R13 */ \
	xmm6 += tmpreg[10];		/* final R9 */ \
	xmm1 = tmpreg[8];		/* Load real-cols row #3 */ \
	xmm1 -= xmm3;			/* final R19 */ \
	xmm3 += tmpreg[8];		/* final R3 */ \
	sd4[2] = xmm0;			/* Save R13 */ \
	sd6[1] = xmm6;			/* Save R9 */ \
	sd6[3] = xmm1;			/* Save R19 */ \
	sd4[0] = xmm3;			/* Save R3 */ \
\
	xmm0 = tmpreg[14];		/* Load real-cols row #7 */ \
	xmm0 -= xmm7;			/* final R15 */ \
	xmm7 += tmpreg[14];		/* final R7 */ \
	xmm1 = tmpreg[12];		/* Load real-cols row #5 */ \
	xmm1 -= xmm4;			/* final R17 */ \
	xmm4 += tmpreg[12];		/* final R5 */ \
	sd8[2] = xmm0;			/* Save R15 */ \
	sd2[1] = xmm7;			/* Save R7 */ \
	sd2[3] = xmm1;			/* Save R17 */ \
	sd8[0] = xmm4;			/* Save R5 */ \

#if 0 //DEF X86_64;
\

#define r5_x10cl_20_reals_last_unfft_preload \
	xmm15 = CONST2_P309; \
	xmm14 = CONST2_P809; \
\
\
/* Core 2 optimal is 204 clocks, currently at 239 clocks. */ \
\

#define r5_x10cl_20_reals_last_unfft(srcreg,srcinc,d1,screg,scoff) \
\
;	/* Apply the 9 twiddle factors to both sets of data */ \
\
	xmm0 = xptr(srcreg+d1);		/* R2 */ \
	xmm1 = xptr(screg+16);		/* cosine/sine */ \
	xmm2 = xmm0;			/* Copy R2 */ \
	xmm0 *= xmm1;			/* A2 = R2 * cosine/sine		; 1-5 */ \
	xmm3 = xptr(srcreg+9*d1);		/* R10 */ \
	xmm4 = xptr(screg+8*32+16);		/* cosine/sine */ \
	xmm5 = xmm3;			/* Copy R10 */ \
	xmm3 *= xmm4;			/* A10 = R10 * cosine/sine		; 2-6 */ \
	xmm6 = xptr(srcreg+d1+32);		/* I2 */ \
	xmm1 *= xmm6;			/* B2 = I2 * cosine/sine		; 3-7 */ \
	xmm7 = xptr(srcreg+9*d1+32);		/* I10 */ \
	xmm4 *= xmm7;			/* B10 = I10 * cosine/sine		; 4-8 */ \
	xmm8 = xptr(srcreg+5*d1);		/* R6 */ \
	xmm8 *= xptr(screg+4*32+16);		/* A6 = R6 * cosine/sine		; 5-9 */ \
\
	xmm0 += xmm6;			/* A2 = A2 + I2				; 6-8	avail 6 */ \
	xmm9 = xptr(srcreg+2*d1);		/* R3 */ \
	xmm10 = xptr(screg+32+16);		/* cosine/sine */ \
	xmm11 = xmm9;			/* Copy R3 */ \
	xmm9 *= xmm10;			/* A3 = R3 * cosine/sine		; 6-10 */ \
\
	xmm3 += xmm7;			/* A10 = A10 + I10			; 7-9	avail 6,7 */ \
	xmm12 = xptr(srcreg+8*d1);		/* R9 */ \
	xmm13 = xptr(screg+7*32+16);		/* cosine/sine */ \
	xmm6 = xmm12;			/* Copy R9 */ \
	xmm12 *= xmm13;			/* A9 = R9 * cosine/sine		; 7-11	avail 7 */ \
\
	xmm1 -= xmm2;			/* B2 = B2 - R2				; 8-10	avail 7,2 */ \
	xmm7 = xptr(srcreg+2*d1+32);		/* I3 */ \
	xmm10 *= xmm7;			/* B3 = I3 * cosine/sine		; 8-12	avail 2 */ \
\
	xmm4 -= xmm5;			/* B10 = B10 - R10			; 9-11	avail 2,5 */ \
	xmm2 = xptr(srcreg+8*d1+32);		/* I9 */ \
	xmm13 *= xmm2;			/* B9 = I9 * cosine/sine		; 9-13	avail 5 */ \
\
	xmm8 += xptr(srcreg+5*d1+32);		/* A6 = A6 + I6				; 10-12 */ \
	xmm5 = xptr(screg);			/* sine */ \
	xmm0 *= xmm5;			/* R2 = A2 * sine			; 10-14	avail none */ \
\
	xmm9 += xmm7;			/* A3 = A3 + I3				; 11-13 */ \
	xmm7 = xptr(screg+8*32);		/* sine */ \
	xmm3 *= xmm7;			/* R10 = A10 * sine			; 11-15 */ \
\
	xmm12 += xmm2;			/* A9 = A9 + I9				; 12-14	avail 2 */ \
	xmm1 *= xmm5;			/* I2 = B2 * sine			; 12-16	avail 2,5 */ \
	xmm2 = xptr(screg+32);		/* sine */ \
\
	xmm10 -= xmm11;			/* B3 = B3 - R3				; 13-15	avail 5,11 */ \
	xmm4 *= xmm7;			/* I10 = B10 * sine			; 13-17	avail 5,11,7 */ \
	xmm5 = xptr(screg+7*32);		/* sine */ \
\
	xmm13 -= xmm6;			/* B9 = B9 - R9				; 14-16	avail 11,7,6 */ \
	xmm9 *= xmm2;			/* R3 = A3 * sine			; 14-18	avail 11,7,6 */ \
	xmm11 = xptr(srcreg+3*d1);		/* R4 */ \
\
	xmm12 *= xmm5;			/* R9 = A9 * sine			; 15-19 */ \
	xmm7 = xptr(screg+2*32+16);		/* cosine/sine */ \
\
	xmm6 = xmm0;			/* Copy R2 */ \
	xmm0 += xmm3;			/* R2+R10				; 16-18	avail none storable 0 */ \
	xmm10 *= xmm2;			/* I3 = B3 * sine			; 16-20	avail 2 storable 0 */ \
	xmm2 = xptr(srcreg+7*d1);		/* R8 */ \
\
	xmm6 -= xmm3;			/* R2-R10				; 17-19	avail 3 storable 0,6 */ \
	xmm13 *= xmm5;			/* I9 = B9 * sine			; 17-21	avail 3,5 storable 0,6 */ \
\
	xmm8 *= xptr(screg+4*32);		/* R6 = A6 * sine			; 18-22	avail 3,5 storable 0,6,8 */ \
	xmm3 = xmm1;			/* Copy I2 */ \
	xmm1 += xmm4;			/* I2+I10				; 18-20	avail 5 storable 0,6,8,1 */ \
\
	xmm3 -= xmm4;			/* I2-I10				; 19-21	avail 5,4 storable 0,6,8,1,3 */ \
	xmm5 = xmm11;			/* Copy R4 */ \
	xmm11 *= xmm7;			/* A4 = R4 * cosine/sine		; 19-23	avail 4 storable 0,6,8,1,3 */ \
	xptr2(g->u.xmm.XMM_COL_MULTS,8*32) = xmm0;	/* Save R2+R10				; 19	avail 4,0 storable 6,8,1,3 */ \
\
	xmm0 = xmm9;			/* Copy R3 */ \
	xmm9 += xmm12;			/* R3+R9				; 20-22	avail 4 storable 6,8,1,3,9 */ \
	xmm4 = xptr(screg+6*32+16);		/* cosine/sine */ \
	t0 = xmm6;		/* Save R2-R10				; 20	avail 6 storable 8,1,3,9 */ \
	xmm6 = xmm2;			/* Copy R8 */ \
	xmm2 *= xmm4;			/* A8 = R8 * cosine/sine		; 20-24	avail none storable 8,1,3,9 */ \
\
	xmm0 -= xmm12;			/* R3-R9				; 21-23	avail 12 storable 8,1,3,9,0 */ \
	xmm12 = xptr(srcreg+3*d1+32);		/* I4 */ \
	xmm7 *= xmm12;			/* B4 = I4 * cosine/sine		; 21-25 */ \
	xptr2(g->u.xmm.XMM_COL_MULTS,8*32+16) = xmm1;	/* Save I2+I10				; 21	avail 1 storable 8,3,9,0 */ \
\
	xmm1 = xmm10;			/* Copy I3 */ \
	xmm10 += xmm13;			/* I3+I9				; 22-24	avail none storable 8,3,9,0,10 */ \
	t1 = xmm3;		/* Save I2-I10				; 22	avail 3 storable 8,9,0,10 */ \
	xmm3 = xptr(srcreg+7*d1+32);		/* I8 */ \
	xmm4 *= xmm3;			/* B8 = I8 * cosine/sine		; 22-26	avail none storable 8,9,0,10 */ \
\
	xmm1 -= xmm13;			/* I3-I9				; 23-25	avail 13 storable 8,9,0,10,1 */ \
	xmm13 = xptr(srcreg+5*d1+32);		/* I6 */ \
	xmm13 *= xptr(screg+4*32+16);		/* B6 = I6 * cosine/sine		; 23-27	avail none storable 8,9,0,10,1 */ \
	xptr2(g->u.xmm.XMM_COL_MULTS,4*32) = xmm8;	/* Save R6				; 23 */ \
	xptr2(g->u.xmm.XMM_COL_MULTS,7*32) = xmm9;	/* Save R3+R9				; 23	avail 8,9 storable 0,10,1 */ \
\
	xmm11 += xmm12;			/* A4 = A4 + I4				; 24-26	avail 8,9,12 storable 0,10,1 */ \
	xmm8 = xptr(srcreg+4*d1);		/* R5 */ \
	xmm9 = xptr(screg+3*32+16);		/* cosine/sine */ \
	xmm12 = xmm8;			/* Copy R5 */ \
	xmm8 *= xmm9;			/* A5 = R5 * cosine/sine		; 24-28	avail none storable 0,10,1 */ \
	t2 = xmm0;		/* Save R3-R9				; 24	avail 0 storable 10,1 */ \
\
	xmm2 += xmm3;			/* A8 = A8 + I8				; 25-27	avail 0,3 storable 10,1 */ \
	xmm0 = xptr(srcreg+6*d1);		/* R7 */ \
	xmm3 = xptr(screg+5*32+16);		/* cosine/sine */ \
	xptr2(g->u.xmm.XMM_COL_MULTS,7*32+16) = xmm10;	/* Save I3+I9				; 25 */ \
	xmm10 = xmm0;			/* Copy R7 */ \
	xmm0 *= xmm3;			/* A7 = R7 * cosine/sine		; 25-29	avail none storable 1 */ \
\
	xmm7 -= xmm5;			/* B4 = B4 - R4				; 26-28	avail 5 storable 1 */ \
	xmm5 = xptr(srcreg+4*d1+32);		/* I5 */ \
	xmm9 *= xmm5;			/* B5 = I5 * cosine/sine		; 26-30	avail none storable 1 */ \
	xptr2(g->u.xmm.XMM_COL_MULTS,32+16) = xmm1;	/* Save I3-I9				; 26	avail 1 */ \
\
	xmm4 -= xmm6;			/* B8 = B8 - R8				; 27-29	avail 1,6 */ \
	xmm1 = xptr(srcreg+6*d1+32);		/* I7 */ \
	xmm3 *= xmm1;			/* B7 = I7 * cosine/sine		; 27-31	avail 6 */ \
\
	xmm13 -= xptr(srcreg+5*d1);		/* B6 = B6 - R6				; 28-30 */ \
	xmm6 = xptr(screg+2*32);		/* sine */ \
	xmm11 *= xmm6;			/* R4 = A4 * sine			; 28-32 avail none */ \
\
	xmm8 += xmm5;			/* A5 = A5 + I5				; 29-31 */ \
	xmm5 = xptr(screg+6*32);		/* sine */ \
	xmm2 *= xmm5;			/* R8 = A8 * sine			; 29-33 */ \
\
	xmm0 += xmm1;			/* A7 = A7 + I7				; 30-32	avail 1 */ \
	xmm7 *= xmm6;			/* I4 = B4 * sine			; 30-34	avail 1,6 */ \
	xmm1 = xptr(screg+3*32);		/* sine */ \
\
	xmm9 -= xmm12;			/* B5 = B5 - R5				; 31-33	avail 6,12 */ \
	xmm4 *= xmm5;			/* I8 = B8 * sine			; 31-35	avail 6,12,5 */ \
	xmm6 = xptr(screg+5*32);		/* sine */ \
\
	xmm3 -= xmm10;			/* B7 = B7 - R7				; 32-34	avail 6,12,5,10 */ \
	xmm8 *= xmm1;			/* R5 = A5 * sine			; 32-36 */ \
	xmm5 = xptr(srcreg+16+d1);		/* R2 */ \
\
	xmm0 *= xmm6;			/* R7 = A7 * sine			; 33-37	avail 12,10 */ \
	xmm10 = xptr(screg+scoff+16);		/* cosine/sine */ \
\
	xmm12 = xmm11;			/* Copy R4 */ \
	xmm11 += xmm2;			/* R4+R8				; 34-36	avail none storable 11 */ \
	xmm9 *= xmm1;			/* I5 = B5 * sine			; 34-38	avail 1 storable 11 */ \
	xmm1 = xptr(srcreg+16+9*d1);		/* R10 */ \
\
	xmm12 -= xmm2;			/* R4-R8				; 35-37	avail 2 storable 11,12 */ \
	xmm3 *= xmm6;			/* I7 = B7 * sine			; 35-39	avail 2,6 storable 11,12 */ \
	xmm2 = xptr(screg+scoff+8*32+16);	/* cosine/sine */ \
\
	xmm6 = xmm7;			/* Copy I4 */ \
	xmm7 += xmm4;			/* I4+I8				; 36-38	avail none storable 11,12,7 */ \
	xmm13 *= xptr(screg+4*32);		/* I6 = B6 * sine			; 36-40	avail none storable 11,12,7,13 */ \
\
	xmm6 -= xmm4;			/* I4-I8				; 37-39	avail 4 storable 11,12,7,13,6 */ \
	xmm4 = xmm5;			/* Copy R2 */ \
	xmm5 *= xmm10;			/* A2 = R2 * cosine/sine		; 1-5 */ \
	xptr2(g->u.xmm.XMM_COL_MULTS,6*32) = xmm11;	/* Save R4+R8				; 37	avail 11 storable 12,7,13,6 */ \
\
	xmm11 = xmm8;			/* Copy R5 */ \
	xmm8 += xmm0;			/* R5+R7				; 38-40	avail none storable 12,7,13,6,8 */ \
	xptr2(g->u.xmm.XMM_COL_MULTS,2*32) = xmm12;	/* Save R4-R8				; 38	avail 12 storable 7,13,6,8 */ \
	xmm12 = xmm1;			/* Copy R10 */ \
	xmm1 *= xmm2;			/* A10 = R10 * cosine/sine		; 2-6 */ \
\
	xmm11 -= xmm0;			/* R5-R7				; 39-41	avail 0 storable 7,13,6,8,11 */ \
	xmm0 = xptr(srcreg+16+d1+32);		/* I2 */ \
	xmm10 *= xmm0;			/* B2 = I2 * cosine/sine		; 3-7 */ \
	xptr2(g->u.xmm.XMM_COL_MULTS,6*32+16) = xmm7;	/* Save I4+I8				; 39	avail 7 storable 13,6,8,11 */ \
\
	xmm7 = xmm9;			/* Copy I5 */ \
	xmm9 += xmm3;			/* I5+I7				; 40-42	avail none storable 13,6,8,11,9 */ \
	xptr2(g->u.xmm.XMM_COL_MULTS,2*32+16) = xmm6;	/* Save I4-I8				; 40	avail 6 storable 13,8,11,9 */ \
	xmm6 = xptr(srcreg+16+9*d1+32);	/* I10 */ \
	xmm2 *= xmm6;			/* B10 = I10 * cosine/sine		; 4-8 */ \
\
	xmm7 -= xmm3;			/* I5-I7				; 41-43	avail 3 storable 13,8,11,9,7 */ \
	xmm3 = xptr(srcreg+16+5*d1);		/* R6 */ \
	xmm3 *= xptr(screg+scoff+4*32+16);	/* A6 = R6 * cosine/sine		; 5-9 */ \
	xptr2(g->u.xmm.XMM_COL_MULTS,4*32+16) = xmm13;	/* Save I6				; 41 */ \
	xptr2(g->u.xmm.XMM_COL_MULTS,5*32) = xmm8;	/* Save R5+R7				; 41 */ \
\
	xmm5 += xmm0;			/* A2 = A2 + I2				; 6-8 */ \
	xmm13 = xptr(srcreg+16+2*d1);		/* R3 */ \
	xmm8 = xptr(screg+scoff+32+16);	/* cosine/sine */ \
	xptr2(g->u.xmm.XMM_COL_MULTS,3*32) = xmm11;	/* Save R5-R7				; 42 */ \
	xmm11 = xmm13;			/* Copy R3 */ \
	xmm13 *= xmm8;			/* A3 = R3 * cosine/sine		; 6-10 */ \
\
	xmm1 += xmm6;			/* A10 = A10 + I10			; 7-9 */ \
	xptr2(g->u.xmm.XMM_COL_MULTS,5*32+16) = xmm9;	/* Save I5+I7				; 43 */ \
	xmm9 = xptr(srcreg+16+8*d1);		/* R9 */ \
	xptr2(g->u.xmm.XMM_COL_MULTS,3*32+16) = xmm7;	/* Save I5-I7				; 44 */ \
	xmm7 = xptr(screg+scoff+7*32+16);	/* cosine/sine */ \
	xmm0 = xmm9;			/* Copy R9 */ \
	xmm9 *= xmm7;			/* A9 = R9 * cosine/sine		; 7-11 */ \
\
	xmm10 -= xmm4;			/* B2 = B2 - R2				; 8-10 */ \
	xmm6 = xptr(srcreg+16+2*d1+32);	/* I3 */ \
	xmm8 *= xmm6;			/* B3 = I3 * cosine/sine		; 8-12 */ \
\
	xmm2 -= xmm12;			/* B10 = B10 - R10			; 9-11 */ \
	xmm4 = xptr(srcreg+16+8*d1+32);	/* I9 */ \
	xmm7 *= xmm4;			/* B9 = I9 * cosine/sine		; 9-13 */ \
\
	xmm3 += xptr(srcreg+16+5*d1+32);	/* A6 = A6 + I6				; 10-12 */ \
	xmm12 = xptr(screg+scoff);		/* sine */ \
	xmm5 *= xmm12;			/* R2 = A2 * sine			; 10-14 */ \
\
	xmm13 += xmm6;			/* A3 = A3 + I3				; 11-13 */ \
	xmm6 = xptr(screg+scoff+8*32);	/* sine */ \
	xmm1 *= xmm6;			/* R10 = A10 * sine			; 11-15 */ \
\
	xmm9 += xmm4;			/* A9 = A9 + I9				; 12-14 */ \
	xmm10 *= xmm12;			/* I2 = B2 * sine			; 12-16 */ \
	xmm4 = xptr(screg+scoff+32);		/* sine */ \
\
	xmm8 -= xmm11;			/* B3 = B3 - R3				; 13-15 */ \
	xmm2 *= xmm6;			/* I10 = B10 * sine			; 13-17 */ \
	xmm12 = xptr(screg+scoff+7*32);	/* sine */ \
\
	xmm7 -= xmm0;			/* B9 = B9 - R9				; 14-16 */ \
	xmm13 *= xmm4;			/* R3 = A3 * sine			; 14-18 */ \
	xmm11 = xptr(srcreg+16+3*d1);		/* R4 */ \
\
	xmm9 *= xmm12;			/* R9 = A9 * sine			; 15-19 */ \
	xmm6 = xptr(screg+scoff+2*32+16);	/* cosine/sine */ \
\
	xmm0 = xmm5;			/* Copy R2 */ \
	xmm5 += xmm1;			/* R2+R10				; 16-18 */ \
	xmm8 *= xmm4;			/* I3 = B3 * sine			; 16-20 */ \
	xmm4 = xptr(srcreg+16+7*d1);		/* R8 */ \
\
	xmm0 -= xmm1;			/* R2-R10				; 17-19 */ \
	xmm7 *= xmm12;			/* I9 = B9 * sine			; 17-21 */ \
\
	xmm3 *= xptr(screg+scoff+4*32);	/* R6 = A6 * sine			; 18-22 */ \
	xmm1 = xmm10;			/* Copy I2 */ \
	xmm10 += xmm2;			/* I2+I10				; 18-20 */ \
\
	xmm1 -= xmm2;			/* I2-I10				; 19-21 */ \
	xmm12 = xmm11;			/* Copy R4 */ \
	xmm11 *= xmm6;			/* A4 = R4 * cosine/sine		; 19-23 */ \
	xptr2(g->u.xmm.XMM_COL_MULTS,10*32+8*32) = xmm5; /* Save R2+R10				; 19 */ \
\
	xmm5 = xmm13;			/* Copy R3 */ \
	xmm13 += xmm9;			/* R3+R9				; 20-22 */ \
	xmm2 = xptr(screg+scoff+6*32+16);	/* cosine/sine */ \
	xptr2(g->u.xmm.XMM_COL_MULTS,10*32+0) = xmm0;	/* Save R2-R10				; 20 */ \
	xmm0 = xmm4;			/* Copy R8 */ \
	xmm4 *= xmm2;			/* A8 = R8 * cosine/sine		; 20-24 */ \
\
	xmm5 -= xmm9;			/* R3-R9				; 21-23 */ \
	xmm9 = xptr(srcreg+16+3*d1+32);	/* I4 */ \
	xmm6 *= xmm9;			/* B4 = I4 * cosine/sine		; 21-25 */ \
	xptr2(g->u.xmm.XMM_COL_MULTS,10*32+8*32+16) = xmm10; /* Save I2+I10			; 21 */ \
\
	xmm10 = xmm8;			/* Copy I3 */ \
	xmm8 += xmm7;			/* I3+I9				; 22-24 */ \
	xptr2(g->u.xmm.XMM_COL_MULTS,10*32+16) = xmm1;	/* Save I2-I10				; 22 */ \
	xmm1 = xptr(srcreg+16+7*d1+32);	/* I8 */ \
	xmm2 *= xmm1;			/* B8 = I8 * cosine/sine		; 22-26 */ \
\
	xmm10 -= xmm7;			/* I3-I9				; 23-25 */ \
	xmm7 = xptr(srcreg+16+5*d1+32);	/* I6 */ \
	xmm7 *= xptr(screg+scoff+4*32+16);	/* B6 = I6 * cosine/sine		; 23-27 */ \
	xptr2(g->u.xmm.XMM_COL_MULTS,10*32+4*32) = xmm3; /* Save R6				; 23 */ \
	xptr2(g->u.xmm.XMM_COL_MULTS,10*32+7*32) = xmm13; /* Save R3+R9				; 23 */ \
\
	xmm11 += xmm9;			/* A4 = A4 + I4				; 24-26 */ \
	xmm3 = xptr(srcreg+16+4*d1);		/* R5 */ \
	xmm13 = xptr(screg+scoff+3*32+16);	/* cosine/sine */ \
	xmm9 = xmm3;			/* Copy R5 */ \
	xmm3 *= xmm13;			/* A5 = R5 * cosine/sine		; 24-28 */ \
	xptr2(g->u.xmm.XMM_COL_MULTS,10*32+32) = xmm5;	/* Save R3-R9				; 24 */ \
\
	xmm4 += xmm1;			/* A8 = A8 + I8				; 25-27 */ \
	xmm5 = xptr(srcreg+16+6*d1);		/* R7 */ \
	xmm1 = xptr(screg+scoff+5*32+16);	/* cosine/sine */ \
	xptr2(g->u.xmm.XMM_COL_MULTS,10*32+7*32+16) = xmm8; /* Save I3+I9			; 25 */ \
	xmm8 = xmm5;			/* Copy R7 */ \
	xmm5 *= xmm1;			/* A7 = R7 * cosine/sine		; 25-29 */ \
\
	xmm6 -= xmm12;			/* B4 = B4 - R4				; 26-28 */ \
	xmm12 = xptr(srcreg+16+4*d1+32);	/* I5 */ \
	xmm13 *= xmm12;			/* B5 = I5 * cosine/sine		; 26-30 */ \
	xptr2(g->u.xmm.XMM_COL_MULTS,10*32+32+16) = xmm10; /* Save I3-I9			; 26 */ \
\
	xmm2 -= xmm0;			/* B8 = B8 - R8				; 27-29 */ \
	xmm10 = xptr(srcreg+16+6*d1+32);	/* I7 */ \
	xmm1 *= xmm10;			/* B7 = I7 * cosine/sine		; 27-31 */ \
\
	xmm7 -= xptr(srcreg+16+5*d1);		/* B6 = B6 - R6				; 28-30 */ \
	xmm0 = xptr(screg+scoff+2*32);	/* sine */ \
	xmm11 *= xmm0;			/* R4 = A4 * sine			; 28-32 */ \
\
	xmm3 += xmm12;			/* A5 = A5 + I5				; 29-31 */ \
	xmm12 = xptr(screg+scoff+6*32);	/* sine */ \
	xmm4 *= xmm12;			/* R8 = A8 * sine			; 29-33 */ \
\
	xmm5 += xmm10;			/* A7 = A7 + I7				; 30-32 */ \
	xmm6 *= xmm0;			/* I4 = B4 * sine			; 30-34 */ \
	xmm10 = xptr(screg+scoff+3*32);	/* sine */ \
\
	xmm13 -= xmm9;			/* B5 = B5 - R5				; 31-33 */ \
	xmm2 *= xmm12;			/* I8 = B8 * sine			; 31-35 */ \
	xmm0 = xptr(screg+scoff+5*32);	/* sine */ \
\
	xmm1 -= xmm8;			/* B7 = B7 - R7				; 32-34 */ \
	xmm3 *= xmm10;			/* R5 = A5 * sine			; 32-36 */ \
\
	xmm5 *= xmm0;			/* R7 = A7 * sine			; 33-37 */ \
\
	xmm9 = xmm11;			/* Copy R4 */ \
	xmm11 += xmm4;			/* R4+R8				; 34-36 */ \
	xmm13 *= xmm10;			/* I5 = B5 * sine			; 34-38 */ \
\
	xmm9 -= xmm4;			/* R4-R8				; 35-37 */ \
	xmm1 *= xmm0;			/* I7 = B7 * sine			; 35-39 */ \
\
	xmm0 = xmm6;			/* Copy I4 */ \
	xmm6 += xmm2;			/* I4+I8				; 36-38 */ \
	xmm7 *= xptr(screg+scoff+4*32);	/* I6 = B6 * sine			; 36-40 */ \
\
	xmm0 -= xmm2;			/* I4-I8				; 37-39 */ \
	xptr2(g->u.xmm.XMM_COL_MULTS,10*32+6*32) = xmm11; /* Save R4+R8				; 37 */ \
\
	xmm11 = xmm3;			/* Copy R5 */ \
	xmm3 += xmm5;			/* R5+R7				; 38-40 */ \
	xptr2(g->u.xmm.XMM_COL_MULTS,10*32+2*32) = xmm9; /* Save R4-R8				; 38 */ \
\
	xmm11 -= xmm5;			/* R5-R7				; 39-41 */ \
	xptr2(g->u.xmm.XMM_COL_MULTS,10*32+6*32+16) = xmm6; /* Save I4+I8			; 39 */ \
\
	xmm6 = xmm13;			/* Copy I5 */ \
	xmm13 += xmm1;			/* I5+I7				; 40-42 */ \
	xptr2(g->u.xmm.XMM_COL_MULTS,10*32+2*32+16) = xmm0; /* Save I4-I8			; 40 */ \
\
	xmm6 -= xmm1;			/* I5-I7				; 41-43 */ \
	xptr2(g->u.xmm.XMM_COL_MULTS,10*32+4*32+16) = xmm7; /* Save I6				; 41 */ \
	xptr2(g->u.xmm.XMM_COL_MULTS,10*32+5*32) = xmm3; /* Save R5+R7				; 41 */ \
/*	xstore	g->u.xmm.XMM_COL_MULTS[10*32][3*32], xmm11 ;; Save R5-R7				; 42 */ \
	xptr2(g->u.xmm.XMM_COL_MULTS,10*32+5*32+16) = xmm13; /* Save I5+I7			; 43 */ \
	xptr2(g->u.xmm.XMM_COL_MULTS,10*32+3*32+16) = xmm6; /* Save I5-I7			; 44 */ \
\
;	/* Do the 20 reals inverse FFT */ \
\
	r5_x20r_unfft(srcreg+d1,d1,xptr(srcreg+16),xptr(srcreg+48),xptr2(g->u.xmm.XMM_COL_MULTS,10*32),1); \
	r5_x20r_unfft(srcreg,d1,xptr(srcreg),xptr(srcreg+32),g->u.xmm.XMM_COL_MULTS,0); \
\
	srcreg += srcinc; \
\
\
/* NOTE: On first call, g->u.xmm.XMM_COL_MULTS[10*32][3*32] is known to be in xmm11. */ \
/* On second call, xmm13 holds XMM_P951. */ \
\

#define r5_x20r_unfft(srcreg,d1,r1A,r1B,tmpreg,firstcall) \
\
;	/* Calculate odd columns derived from real inputs (even rows) */ \
\
	xmm0 = tmpreg[2];		/* r3-r9 */ \
	xmm1 = xmm14; \
	xmm1 *= xmm0;			/* .809(r3-r9)				; 1-5 */ \
	xmm2 = xmm15; \
	xmm2 *= xmm0;			/* .309(r3-r9)				; 2-6 */ \
\
	IF firstcall EQ 0; \
	xmm11 = tmpreg[6];		/* r5-r7 */ \
	ENDIF; \
	xmm4 = xmm15; \
	xmm4 *= xmm11;			/* .309(r5-r7)				; 3-7 */ \
	xmm5 = xmm14; \
	xmm5 *= xmm11;			/* .809(r5-r7)				; 4-8 */ \
\
	xmm6 = r1B;			/* r1-r11 */ \
	xmm0 -= xmm6;			/* -(r1-(r3-r9)-r11)			; 5-7 (1) */ \
	xmm7 = tmpreg[0];			/* r2-r10 */ \
	IF firstcall EQ 1; \
	xmm13 = CONST2_P951; \
	ENDIF; \
	xmm8 = xmm13; \
	xmm8 *= xmm7;			/* .951(r2-r10)				; 5-9 */ \
\
	xmm1 += xmm6;			/* r1+.809(r3-r9)-r11			; 6-8 */ \
	xmm7 *= CONST2_P588;			/* .588(r2-r10)				; 6-10 */ \
\
	xmm6 -= xmm2;			/* r1-.309(r3-r9)-r11			; 7-9	avail 2 */ \
	xmm9 = tmpreg[4];		/* r4-r8 */ \
	xmm10 = CONST2_P588; \
	xmm10 *= xmm9;			/* .588(r4-r8)				; 7-11 */ \
\
	xmm11 -= xmm0;			/* r1-(r3-r9)+(r5-r7)-r11		; 8-10	avail 2,0 */ \
	xmm9 *= xmm13;			/* .951(r4-r8)				; 8-12 */ \
\
	xmm1 += xmm4;			/* r1+.809(r3-r9)+.309(r5-r7)-r11	; 9-11 (8) avail 2,0,4 */ \
	xmm3 = tmpreg[17];		/* i2+i10 */ \
	xmm12 = xmm15; \
	xmm12 *= xmm3;			/* .309(i2+i10)				; 9-13 */ \
\
	xmm6 -= xmm5;			/* r1-.309(r3-r9)-.809(r5-r7)-r11	; 10-12 (9) avail 2,0,4,5 */ \
	xmm2 = xmm14; \
	xmm2 *= xmm3;			/* .809(i2+i10)				; 10-14	avail 0,4,5 */ \
\
	xmm0 = tmpreg[9];		/* i6 */ \
	xmm3 += xmm0;			/* (i2+i10)+i6				; 11-13	avail 4,5 */ \
	xmm4 = tmpreg[13];		/* i4+i8 */ \
	xmm5 = xmm14; \
	xmm5 *= xmm4;			/* .809(i4+i8)				; 11-15 avail none */ \
\
	xmm8 += xmm10;			/* .951(r2-r10)+.588(r4-r8)		; 12-14	avail 10 */ \
	xmm10 = xmm15; \
	xmm10 *= xmm4;			/* .309(i4+i8)				; 12-16	avail none */ \
\
	xmm7 -= xmm9;			/* .588(r2-r10)-.951(r4-r8)		; 13-15 avail 9 */ \
\
	xmm12 += xmm0;			/* .309(i2+i10)+i6			; 14-16 */ \
\
	xmm2 -= xmm0;			/* .809(i2+i10)-i6			; 15-17	avail 9,0 */ \
	xmm9 = tmpreg[15];		/* i3+i9 */ \
	xmm0 = CONST2_P588; \
	xmm0 *= xmm9;			/* .588(i3+i9)				; 15-19 (13) avail none */ \
\
	xmm3 -= xmm4;			/* (i2+i10)-(i4+i8)+i6			; 16-18 avail 4 */ \
	xmm9 *= xmm13;			/* .951(i3+i9)				; 16-20 */ \
\
	xmm12 += xmm5;			/* .309(i2+i10)+.809(i4+i8)+i6		; 17-19	avail 4,5 */ \
	xmm4 = tmpreg[11];		/* i5+i7 */ \
	xmm5 = xmm13; \
	xmm5 *= xmm4;			/* .951(i5+i7)				; 17-21	avail none */ \
\
	xmm2 += xmm10;			/* .809(i2+i10)+.309(i4+i8)-i6		; 18-20	avail 10 */ \
	xmm4 *= CONST2_P588;			/* .588(i5+i7)				; 18-22 */ \
\
	xmm10 = xmm1;			/* Copy odd-real-cols row #2 */ \
	xmm1 -= xmm8;			/* real-cols row #10 (odd#2 - even#2)	; 19-21 (15) */ \
	xmm8 += xmm10;			/* real-cols row #2 (odd#2 + even#2)	; 20-22 (16) */ \
\
	xmm10 = xmm6;			/* Copy odd-real-cols row #4 */ \
	xmm6 -= xmm7;			/* real-cols row #8 (odd#4 - even#4)	; 21-23 (16) */ \
	xmm7 += xmm10;			/* real-cols row #4 (odd#4 + even#4)	; 22-24 (17) */ \
\
	xmm0 += xmm5;			/* .588(i3+i9)+.951(i5+i7)		; 23-25 (22) avail 10,5 */ \
	xmm10 = tmpreg[16];		/* r2+r10 */ \
\
	xmm9 -= xmm4;			/* .951(i3+i9)-.588(i5+i7)		; 24-26 (23) avail 5,4 */ \
\
	xmm4 = xmm11;			/* Copy  real-cols row #6 */ \
	xmm11 -= xmm3;			/* final R16				; 25-27 (19) */ \
\
	xmm3 += xmm4;			/* final R6				; 26-28 */ \
	xmm5 = xmm14; \
\
	xmm4 = xmm12;			/* Copy even-imag-cols row #2 */ \
	xmm12 -= xmm0;			/* imag-cols row #10 (even#2 - odd#2)	; 27-29 */ \
\
	xmm0 += xmm4;			/* imag-cols row #2 (even#2 + odd#2)	; 28-30 */ \
	xptr(srcreg+48) = xmm11;		/* Save R16				; 28	avail 3,4 */ \
	xmm11 = xmm15; \
\
	xmm4 = xmm2;			/* Copy even-imag-cols row #4 */ \
	xmm2 -= xmm9;			/* imag-cols row #8 (even#4 - odd#4)	; 29-31 */ \
	xptr(srcreg+16) = xmm3;		/* Save R6				; 29 */ \
\
	xmm9 += xmm4;			/* imag-cols row #4 (even#4 + odd#4)	; 30-32	avail 11,4 */ \
	xmm3 = tmpreg[12];		/* r4+r8 */ \
\
	xmm4 = xmm1;			/* Copy real-cols row #10 */ \
	xmm1 -= xmm12;			/* final R12				; 31-33 */ \
\
	xmm12 += xmm4;			/* final R10				; 32-34 */ \
\
	xmm4 = xmm8;			/* Copy real-cols row #2 */ \
	xmm8 -= xmm0;			/* final R20				; 33-35 */ \
\
	xmm0 += xmm4;			/* final R2				; 34-36 */ \
	xptr(srcreg+2*d1+32) = xmm1;		/* Save R12				; 34	avail 1,4 */ \
	xmm1 = xmm15; \
\
	xmm4 = xmm6;			/* Copy real-cols row #8 */ \
	xmm6 -= xmm2;			/* final R14				; 35-37 */ \
	xmm5 *= xmm10;			/* .809(r2+r10)				; 1-5 */ \
	xptr(srcreg+8*d1+16) = xmm12;		/* Save R10				; 35 */ \
\
	xmm2 += xmm4;			/* final R8				; 36-38 */ \
	xmm11 *= xmm10;			/* .309(r2+r10)				; 2-6 */ \
	xptr(srcreg+8*d1+48) = xmm8;		/* Save R20				; 36	avail 12,8,4 */ \
\
	xmm4 = xmm7;			/* Load real-cols row #4 */ \
	xmm7 -= xmm9;			/* final R18				; 37-39 */ \
	xmm1 *= xmm3;			/* .309(r4+r8)				; 3-7 */ \
	xptr(srcreg+2*d1) = xmm0;		/* Save R2				; 37 */ \
\
	xmm9 += xmm4;			/* final R4				; 38-40 */ \
	xmm12 = xmm14; \
	xmm12 *= xmm3;			/* .809(r4+r8)				; 4-8 */ \
	xptr(srcreg+6*d1+32) = xmm6;		/* Save R14				; 38	avail 8,0,6,4 */ \
\
	xmm8 = tmpreg[8];		/* r6 */ \
	xmm10 += xmm8;			/* (r2+r10)+r6				; 5-7 */ \
	xmm0 = tmpreg[14];		/* r3+r9 */ \
	xmm6 = xmm15; \
	xmm6 *= xmm0;			/* .309(r3+r9)				; 5-9 */ \
	xptr(srcreg+4*d1+16) = xmm2;		/* Save R8				; 39 */ \
\
	xmm5 -= xmm8;			/* .809(r2+r10)-r6			; 6-8 */ \
	xmm4 = xmm14; \
	xmm4 *= xmm0;			/* .809(r3+r9)				; 6-10 */ \
	xptr(srcreg+4*d1+48) = xmm7;		/* Save R18				; 40 */ \
\
	xmm11 += xmm8;			/* .309(r2+r10)+r6			; 7-9	avail 8 */ \
	xmm2 = tmpreg[10];		/* r5+r7 */ \
	xmm7 = xmm14; \
	xmm7 *= xmm2;			/* .809(r5+r7)				; 7-11 */ \
	xptr(srcreg+6*d1) = xmm9;		/* Save R4				; 41 */ \
\
	xmm10 += xmm3;			/* (r2+r10)+(r4+r8)+r6			; 8-10	avail 8,11 */ \
	xmm9 = xmm15; \
	xmm9 *= xmm2;			/* .309(r5+r7)				; 8-12 */ \
\
	xmm5 -= xmm1;			/* .809(r2+r10)-.309(r4+r8)-r6		; 9-11	avail 8,11,1 */ \
	xmm8 = tmpreg[1];		/* i2-i10 */ \
	xmm3 = CONST2_P588; \
	xmm3 *= xmm8;			/* .588(i2-i10)				; 9-13	avail 1 */ \
\
	xmm11 -= xmm12;			/* .309(r2+r10)-.809(r4+r8)+r6		; 10-12	avail 1,12 */ \
	xmm8 *= xmm13;			/* .951(i2-i10)				; 10-14 */ \
\
	xmm0 += xmm2;			/* (r3+r9)+(r5+r7)			; 11-13	avail 1,12,2 */ \
	xmm1 = tmpreg[5];		/* i4-i8 */ \
	xmm12 = xmm13; \
	xmm12 *= xmm1;			/* .951(i4-i8)				; 11-15	avail 2 */ \
\
	xmm6 -= xmm7;			/* .309(r3+r9)-.809(r5+r7)		; 12-14	avail 2,7 */ \
	xmm1 *= CONST2_P588;			/* .588(i4-i8)				; 12-16 */ \
\
	xmm9 -= xmm4;			/* -.809(r3+r9)+.309(r5+r7)		; 13-15	avail 2,7,4 */ \
	xmm2 = tmpreg[3];		/* i3-i9 */ \
	xmm7 = xmm13; \
	xmm7 *= xmm2;			/* .951(i3-i9)				; 13-17	avail 4 */ \
\
	xmm4 = r1A;			/* r1+r11 */ \
	xmm0 += xmm4;			/* r1+(r3+r9)+(r5+r7)+r11		; 14-16	avail none */ \
	xmm2 *= CONST2_P588;			/* .588(i3-i9)				; 14-18 */ \
\
	xmm6 += xmm4;			/* r1+.309(r3+r9)-.809(r5+r7)+r11	; 15-17 */ \
	xmm9 += xmm4;			/* r1-.809(r3+r9)+.309(r5+r7)+r11	; 16-18 avail 4 */ \
\
	xmm3 += xmm12;			/* .588(i2-i10)+.951(i4-i8)		; 17-19	avail 4,12 */ \
	xmm4 = tmpreg[7];		/* i5-i7 */ \
	xmm12 = CONST2_P588; \
	xmm12 *= xmm4;			/* .588(i5-i7)				; 17-21 (15) avail none */ \
\
	xmm8 -= xmm1;			/* .951(i2-i10)-.588(i4-i8)		; 18-20	avail 1 */ \
	xmm4 *= xmm13;			/* .951(i5-i7)				; 18-22 (16) */ \
\
	xmm1 = xmm10;			/* Copy even-real-cols row #1 */ \
	xmm10 += xmm0;			/* real-cols row #1 (and final R1)	; 19-21 */ \
	xmm0 -= xmm1;			/* real-cols row #11 (and final R11)	; 20-22 */ \
\
	xmm1 = xmm5;			/* Copy even-real-cols row #3 */ \
	xmm5 += xmm6;			/* real-cols row #3			; 21-23 */ \
	xmm6 -= xmm1;			/* real-cols row #9			; 22-24 */ \
	xptr(srcreg) = xmm10;			/* Save final R1			; 22	avail 10,1 */ \
\
	xmm7 += xmm12;			/* .951(i3-i9)+.588(i5-i7)		; 23-25	avail 10,12,1 */ \
	xptr(srcreg+32) = xmm0;		/* Save final R11			; 23	avail 10,12,0,1 */ \
	xmm2 -= xmm4;			/* .588(i3-i9)-.951(i5-i7)		; 24-26	avail 10,12,0,4,1 */ \
\
	xmm1 = xmm11;			/* Copy even-real-cols row #5 */ \
	xmm11 += xmm9;			/* real-cols row #5			; 25-27 */ \
	xmm9 -= xmm1;	 		/* real-cols row #7			; 26-28 */ \
\
	xmm1 = xmm3;			/* Copy even-imag-cols row #3 */ \
	xmm3 -= xmm7;			/* imag-cols row #9 (even#3 - odd#3)	; 27-29 */ \
	xmm7 += xmm1;			/* imag-cols row #3 (even#3 + odd#3)	; 28-30 */ \
	xmm1 = xmm8;			/* Copy even-imag-cols row #5 */ \
	xmm8 -= xmm2;			/* imag-cols row #7 (even#5 - odd#5)	; 29-31 */ \
	xmm2 += xmm1;			/* imag-cols row #5 (even#5 + odd#5)	; 30-32 */ \
\
	xmm1 = xmm6;			/* Copy real-cols row #9 */ \
	xmm6 -= xmm3;			/* final R13				; 31-33 */ \
	xmm3 += xmm1;			/* final R9				; 32-34 */ \
	xmm1 = xmm5;			/* Copy real-cols row #3 */ \
	xmm5 -= xmm7;			/* final R19				; 33-35 */ \
	xmm7 += xmm1;			/* final R3				; 34-36 */ \
	xptr(srcreg+4*d1+32) = xmm6;		/* Save R13				; 34 */ \
\
	xmm1 = xmm9;			/* Copy real-cols row #7 */ \
	xmm9 -= xmm8;			/* final R15				; 35-37 */ \
	xptr(srcreg+6*d1+16) = xmm3;		/* Save R9				; 35 */ \
	xmm8 += xmm1;			/* final R7				; 36-38 */ \
	xptr(srcreg+6*d1+48) = xmm5;		/* Save R19				; 36 */ \
\
	xmm1 = xmm11;			/* Copy real-cols row #5 */ \
	xmm11 -= xmm2;			/* final R17				; 37-39 */ \
	xptr(srcreg+4*d1) = xmm7;		/* Save R3				; 37 */ \
	xmm2 += xmm1;			/* final R5				; 38-40 */ \
	xptr(srcreg+8*d1+32) = xmm9;		/* Save R15				; 38 */ \
	xptr(srcreg+2*d1+16) = xmm8;		/* Save R7				; 39 */ \
	xptr(srcreg+2*d1+48) = xmm11;		/* Save R17				; 40 */ \
	xptr(srcreg+8*d1) = xmm2;		/* Save R5				; 41 */ \
\

#endif
