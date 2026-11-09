#pragma once

/* Copyright 2009-2016 - Mersenne Research, Inc.  All rights reserved */ \
/* Author:  George Woltman */ \
/* Email: woltman@alum.mit.edu */ \
/* */ \
\
/* */ \
/* */ \
/* All new macros for version 26 of gwnum.  Do a radix-3 step in an FFT. */ \
/* The forward FFT macros multiply by the sin/cos values at the end of the macro */ \
/* and the inverse FFTs multiply by the sin/cos values at the start of the macro. */ \
/* We also implement the Daniel J. Bernstein (DJB) "exponent-1" idea to save sin/cos */ \
/* memory. */ \
/* */ \
/* */ \
\
/* */ \
/* ************************************* three-complex-djbfft variants ****************************************** */ \
/* */ \
\

#define r3_x3cl_three_complex_djbfft_preload \
	r3_x3c_djbfft_partial_mem_preload; \
\
\

#define r3_x3cl_three_complex_djbfft(srcreg,srcinc,d1,screg) { \
	vec2f64 *sd0 = (vec2f64*)(srcreg); \
	vec2f64 *sd1 = (vec2f64*)(srcreg+(d1)); \
	vec2f64 *sd2 = (vec2f64*)(srcreg+2*(d1)); \
	xmm2 = sd2[2];		/* R3 */ \
	xmm5 = sd2[3];		/* I3 */ \
	r3_x3c_djbfft_partial_mem(xmm0,xmm1,xmm2,xmm3,xmm4,xmm5,xmm6,xmm7, sd0[2],sd1[2],sd0[3],sd1[3], screg,srcreg+srcinc,d1); \
	xmm6 = sd2[0];		/* R3 */ \
	xmm7 = sd2[1];		/* I3 */ \
	sd1[2] = xmm2;		/* Save R1 */ \
	sd1[3] = xmm5;		/* Save I1 */ \
	sd2[0] = xmm0;		/* Save R2 */ \
	sd2[1] = xmm1;		/* Save I2 */ \
	sd2[2] = xmm4;		/* Save R3 */ \
	sd2[3] = xmm3;		/* Save I3 */ \
	r3_x3c_djbfft_partial_mem(xmm0,xmm1,xmm6,xmm3,xmm4,xmm7,xmm2,xmm5, sd0[0],sd1[0],sd0[1],sd1[1], screg,srcreg+srcinc+2*d1,0); \
	sd0[0] = xmm6;			/* Save R1 */ \
	sd0[1] = xmm7;		/* Save I1 */ \
	sd0[2] = xmm0;		/* Save R2 */ \
	sd0[3] = xmm1;		/* Save I2 */ \
	sd1[0] = xmm4;		/* Save R3 */ \
	sd1[1] = xmm3;		/* Save I3 */ \
	srcreg += srcinc; \
}

/* Used in first levels of pass 2.  No swizzling. */ \
\

#define r3_f3cl_three_complex_djbfft_preload \
	r3_f3cl_three_complex_djbfft_common_preload; \
\
\

#define r3_f3cl_three_complex_djbfft(srcreg,srcinc,d1,screg,scoff) \
	r3_f3cl_three_complex_djbfft_common(srcreg,rbx,srcinc,d1,screg,scoff); \
\
\
/* Used in pass 2 where the memory layout is different and the two */ \
/* three-complex-djbffts use different sin/cos data. */ \
\

#define r3_nf3cl_three_complex_djbfft_preload \
	r3_f3cl_three_complex_djbfft_common_preload; \
\
\

#define r3_nf3cl_three_complex_djbfft(srcreg,srcinc,d1,screg,scoff) \
	r3_f3cl_three_complex_djbfft_common(srcreg,0,srcinc,d1,screg,scoff); \
\
\
/* The common 3-complex djbfft code. */ \
\

#define r3_f3cl_three_complex_djbfft_common_preload \
	r3_x3c_djbfft_partial_mem_preload; \
\
\

#define r3_f3cl_three_complex_djbfft_common(srcreg,srcoff,srcinc,d1,screg_a,scoff) { \
	uintptr_t src_off = (uintptr_t)(srcreg + srcoff); \
	vec2f64 *sd0 = (vec2f64*)(srcreg); \
	vec2f64 *sd1 = (vec2f64*)(srcreg+(d1)); \
	vec2f64 *sd2 = (vec2f64*)(srcreg+2*(d1)); \
	vec2f64 *bd0 = (vec2f64*)(src_off); \
	vec2f64 *bd1 = (vec2f64*)(src_off+(d1)); \
	vec2f64 *bd2 = (vec2f64*)(src_off+2*(d1)); \
	char *screg = (char*)(screg_a); /* improves clang compiler */ \
	xmm2 = bd2[1];	/* R3 */ \
	xmm5 = bd2[3];	/* I3 */ \
	r3_x3c_djbfft_partial_mem(xmm0,xmm1,xmm2,xmm3,xmm4,xmm5,xmm6,xmm7, bd1[0],bd0[1],bd1[2],bd0[3], screg+scoff,srcreg+srcinc,d1); \
	xmm6 = bd1[1];	/* R3 */ \
	xmm7 = bd1[3];	/* I3 */ \
	sd0[1] = xmm2;		/* Save R1 */ \
	sd0[3] = xmm5;		/* Save I1 */ \
	sd1[1] = xmm0;		/* Save R2 */ \
	sd1[3] = xmm1;		/* Save I2 */ \
	sd2[1] = xmm4;		/* Save R3 */ \
	sd2[3] = xmm3;		/* Save I3 */ \
	r3_x3c_djbfft_partial_mem(xmm0,xmm1,xmm6,xmm3,xmm4,xmm7,xmm2,xmm5, bd0[0],bd2[0],bd0[2],bd2[2], screg,srcreg+srcinc+2*d1,0); \
	sd0[0] = xmm6;			/* Save R1 */ \
	sd0[2] = xmm7;		/* Save I1 */ \
	sd1[0] = xmm0;		/* Save R2 */ \
	sd1[2] = xmm1;		/* Save I2 */ \
	sd2[0] = xmm4;		/* Save R3 */ \
	sd2[2] = xmm3;		/* Save I3 */ \
	srcreg += srcinc; \
}

/* Do a 3-complex FFT.  The input values are R1+R4i, R2+R5i, R3+R6i */ \
/* A 3-complex FFT is: */ \
/* Res1:  (R1+R2+R3) + (I1+I2+I3)i */ \
/* Res2:  (R1-.5R2-.866I2-.5R3+.866I3) + (I1-.5I2+.866R2-.5I3-.866R3)i */ \
/* Res3:  (R1-.5R2+.866I2-.5R3-.866I3) + (I1-.5I2-.866R2-.5I3+.866R3)i */ \
/* Then multiply 2 of the 3 results by twiddle factors. */ \
\

#define r3_x3c_djbfft_partial_mem_preload \
\
\

#define r3_x3c_djbfft_partial_mem(r1,r2,r3,r4,r5,r6,t1,t2,memr1,memr2,memi1,memi2,screg,pre1,pre2) \
	r2 = memr2; \
	t1 = r3; \
	r3 += r2;			/* R2 + R3 */ \
	r5 = memi2; \
	t2 = r6; \
	r6 += r5;			/* I2 + I3 */ \
	r2 -= t1;			/* R2 - R3 */ \
	t1 = r3;			/* 0.5 * (R2 + R3) */ \
	r5 -= t2;			/* I2 - I3 */ \
	t2 = r6;			/* 0.5 * (I2 + I3) */ \
	r1 = memr1; \
	r3 += r1;			/* R1 + R2 + R3 (final R1) */ \
	r2 *= CONST2_P866;	/* 0.866 * (R2 - R3) */ \
	r1 -= t1 * 0.5;		/* (R1-.5R2-.5R3) */ \
	r5 *= CONST2_P866;	/* 0.866 * (I2 - I3) */ \
	r4 = memi1; \
	r6 += r4;			/* I1 + I2 + I3 (final I1) */ \
	r4 -= t2 * 0.5;		/* (I1-.5I2-.5I3) */ \
\
	xprefetchw(u8ptr(pre1)); \
\
	t1 = r1; \
	r1 -= r5;			/* Final R2 */ \
	t2 = r4; \
	r4 -= r2;			/* Final I3 */ \
	r5 += t1;			/* Final R3 */ \
	r2 += t2;			/* Final I2 */ \
\
	/*  IF 0 pre2 NE 0; \
		xprefetchw(u8ptr(pre1+pre2)); \
		ENDIF */ \
\
	t1 = r1;			/* Copy R2 */ \
	t2 = xptr(screg+16); \
	r1 = t1 * t2 - r2;	/* A2 = R2 * cosine/sine - I2 */ \
	r2 = r2 * t2 + t1;	/* B2 = I2 * cosine/sine + R2 */ \
\
	t1 = r5;			/* Copy R3 */ \
	r5 = r5 * t2 + r4;	/* A3 = R3 * cosine/sine  + I3 */ \
	r4 = r4 * t2 - t1;	/* B3 = I3 * cosine/sine - R3 */ \
\
	t2 = xptr(screg); \
	r1 *= t2;			/* A2 = A2 * sine (new R2) */ \
	r2 *= t2;			/* B2 = B2 * sine (new I2) */ \
	r5 *= t2;			/* A3 = A3 * sine (new R3) */ \
	r4 *= t2;			/* B3 = B3 * sine (new I3) */ \
\
\
/* AMD K8 version of some of the above macros */ \

#if 0 // (@INSTR(,%xarch,<K8>) NE 0);


#define r3_x3c_djbfft_partial_mem(r1,r2,r3,r4,r5,r6,t1,t2,memr1,memr2,memi1,memi2,screg,pre1,pre2) \
	r2 = memr2; \
	r2 -= r3;			/* R2 - R3 */ \
	r3 += memr2;		/* R2 + R3 */ \
	r5 = memi2; \
	r5 -= r6;			/* I2 - I3 */ \
	r6 += memi2;		/* I2 + I3 */ \
\
	r1 = 0.5; \
	r1 *= r3;			/* 0.5 * (R2 + R3) */ \
	r4 = 0.5; \
	r4 *= r6;			/* 0.5 * (I2 + I3) */ \
\
	t1 = memr1; \
	r3 += t1;			/* R1 + R2 + R3 (final R1) */ \
\
	r2 *= CONST2_P866;		/* 0.866 * (R2 - R3) */ \
	t1 -= r1;			/* (R1-.5R2-.5R3) */ \
	r5 *= CONST2_P866;		/* 0.866 * (I2 - I3) */ \
	t2 = memi1; \
	r6 += t2;			/* I1 + I2 + I3 (final I1) */ \
	t2 -= r4;			/* (I1-.5I2-.5I3) */ \
\
	xprefetchw(u8ptr(pre1)); \
\
	t1 -= r5;			/* Final R2 */ \
	r5 += r5; /* multwo */ \
	t2 -= r2;			/* Final I3 */ \
	r2 += r2; /* multwo */ \
	r5 += t1;			/* Final R3 */ \
	r2 += t2;			/* Final I2 */ \
\
	IF 0 pre2 NE 0; \
	xprefetchw(u8ptr(pre1+pre2)); \
	ENDIF \
\
	r1 = xptr(screg+16); \
	r1 *= t1;			/* A2 = R2 * cosine/sine */ \
\
	r4 = xptr(screg+16); \
	r4 *= t2;			/* B3 = I3 * cosine/sine */ \
\
	r1 -= r2;			/* A2 = A2 - I2 */ \
	r2 *= xptr(screg+16);		/* B2 = I2 * cosine/sine */ \
\
	r4 -= r5;			/* B3 = B3 - R3 */ \
	r5 *= xptr(screg+16);		/* A3 = R3 * cosine/sine */ \
\
	r2 += t1;			/* B2 = B2 + R2 */ \
	r5 += t2;			/* A3 = A3 + I3 */ \
\
	t2 = xptr(screg); \
	r1 *= t2;			/* A2 = A2 * sine (new R2) */ \
	r4 *= t2;			/* B3 = B3 * sine (new I3) */ \
	r2 *= t2;			/* B2 = B2 * sine (new I2) */ \
	r5 *= t2;			/* A3 = A3 * sine (new R3) */ \
\

#endif

/* 64-bit version of some of the above macros - use the extra registers */ \

#if 0 // (@INSTR(,%xarch,<CORE>) NE 0) OR (@INSTR(,%xarch,<BLEND>) NE 0);
#if 0  // IFDEF X86_64;

/* Core 2 optimal is 32 clocks.  Currently at 35.5 clocks. */ \
\

#define r3_x3cl_three_complex_djbfft_preload \
	xmm14 = 0.5; \
	xmm15 = CONST2_P866; \
\
\

#define r3_x3cl_three_complex_djbfft(srcreg,srcinc,d1,screg) \
	xmm0 = xptr(srcreg+d1+32);		/* R2 */ \
	xmm1 = xptr(srcreg+2*d1+32);		/* R3 */ \
	xmm2 = xmm0;			/* Copy R2 */ \
	xmm0 += xmm1;			/* R2 + R3			; 1-3 */ \
\
	xmm3 = xptr(srcreg+d1+48);		/* I2 */ \
	xmm4 = xptr(srcreg+2*d1+48);		/* I3 */ \
	xmm5 = xmm3;			/* Copy I2 */ \
	xmm3 += xmm4;			/* I2 + I3			; 2-4 */ \
\
	xmm2 -= xmm1;			/* R2 - R3			; 3-5 */ \
	xmm1 = xmm14;			/* 0.5 */ \
\
	xmm5 -= xmm4;			/* I2 - I3			; 4-6			avail 1,4,7-13 */ \
	xmm1 *= xmm0;			/* 0.5 * (R2 + R3)		; 4-8			avail 4,7-13 */ \
\
	xmm4 = xptr(srcreg+d1);		/*#2 R2 */ \
	xmm6 = xptr(srcreg+2*d1);		/*#2 R3 */ \
	xmm7 = xmm4;			/*#2 Copy R2						avail 8-13 */ \
	xmm4 += xmm6;			/*#2 R2 + R3			; 5-7 */ \
	xmm8 = xmm14;			/* 0.5 */ \
	xmm8 *= xmm3;			/* 0.5 * (I2 + I3)		; 5-9			avail 9-13 */ \
\
	xmm9 = xptr(srcreg+d1+16);		/*#2 I2 */ \
	xmm10 = xptr(srcreg+2*d1+16);		/*#2 I3 */ \
	xmm11 = xmm9;			/*#2 Copy I2 */ \
	xmm9 += xmm10;			/*#2 I2 + I3			; 6-8			avail 12-13 */ \
	xmm2 *= xmm15;			/* 0.866 * (R2 - R3)		; 6-10 */ \
\
	xmm7 -= xmm6;			/*#2 R2 - R3			; 7-9			avail 6,12-13 */ \
	xmm5 *= xmm15;			/* 0.866 * (I2 - I3)		; 7-11 */ \
	xmm6 = xmm14;			/* 0.5 */ \
\
	xmm11 -= xmm10;			/*#2 I2 - I3			; 8-10			avail 6,10,12-13 */ \
	xmm6 *= xmm4;			/*#2 0.5 * (R2 + R3)		; 8-12			avail 10,12-13 */ \
	xmm12 = xmm14;			/* 0.5 */ \
\
	xmm10 = xptr(srcreg+32);		/* R1 */ \
	xmm0 += xmm10;			/* R1 + R2 + R3 (final R1)	; 9-11			avail 12-13 storable 0 */ \
	xmm12 *= xmm9;			/*#2 0.5 * (I2 + I3)		; 9-13			avail 13 storable 0 */ \
\
	xmm10 -= xmm1;			/* (R1-.5R2-.5R3)		; 10-12			avail 1,13 storable 0 */ \
	xmm7 *= xmm15;			/*#2 0.866 * (R2 - R3)		; 10-14 */ \
\
	xmm1 = xptr(srcreg+48);		/* I1 */ \
	xmm3 += xmm1;			/* I1 + I2 + I3 (final I1)	; 11-13			avail 13 storable 0,3 */ \
	xmm11 *= xmm15;			/*#2 0.866 * (I2 - I3)		; 11-15 */ \
\
	xmm1 -= xmm8;			/* (I1-.5I2-.5I3)		; 12-14			avail 8,13 storable 0,3 */ \
	xmm8 = xptr(srcreg);			/*#2 R1 */ \
	xptr(srcreg+d1+32) = xmm0;		/* Save R1			; 12			avail 0,13 storable 3 */ \
\
	xmm4 += xmm8;			/*#2 R1 + R2 + R3 (final R1)	; 13-15			avail 0,13 storable 3,4 */ \
	xmm13 = xptr(screg+16);		/* cosine/sine						avail 0 storable 3,4 */ \
\
	xmm8 -= xmm6;			/*#2 (R1-.5R2-.5R3)		; 14-16			avail 0,6 storable 3,4 */ \
	xmm0 = xmm13;			/* Copy cosine/sine					avail 6 storable 3,4 */ \
	xmm6 = xptr(srcreg+16);		/*#2 I1 */ \
	xptr(srcreg+d1+48) = xmm3;		/* Save I1			; 14			avail 3 storable 4 */ \
\
	xmm9 += xmm6;			/*#2 I1 + I2 + I3 (final I1)	; 15-17			avail 3 storable 4,9 */ \
	xmm3 = xmm13;			/* Copy cosine/sine					avail none storable 4,9 */ \
\
	xmm6 -= xmm12;			/*#2 (I1-.5I2-.5I3)		; 16-18			avail 12 storable 4,9 */ \
	xptr(srcreg) = xmm4;			/*#2 Save R1			; 16			avail 4,12 storable 9 */ \
\
	xmm12 = xmm10;			/* Copy R1-.5R2-.5R3 */ \
	xmm10 -= xmm5;			/* Final R2 (R-I)		; 17-19			avail 4 storable 9 */ \
\
	xmm5 += xmm12;			/* Final R3 (R+I)		; 18-20 */ \
	xptr(srcreg+16) = xmm9;		/*#2 Save I1			; 18			avail 4,9,12 */ \
\
	xmm12 = xmm1;			/* Copy I1-.5I2-.5I3					avail 4,9 */ \
	xmm1 += xmm2;			/* Final I2 (I+R)		; 19-21 */ \
\
	xmm12 -= xmm2;			/* Final I3 (I-R)		; 20-22			avail 2,4,9 */ \
	xmm2 = xmm13;			/* Copy cosine/sine */ \
	xmm0 *= xmm10;			/* A2 = R2 * cosine/sine	; 20-24  		avail 4,9 */ \
\
	xmm4 = xmm8;			/*#2 Copy R1-.5R2-.5R3					avail 9 */ \
	xmm8 -= xmm11;			/*#2 Final R2 (R-I)		; 21-23 */ \
	xmm3 *= xmm5;			/* A3 = R3 * cosine/sine	; 21-25 */ \
\
	xmm11 += xmm4;			/*#2 Final R3 (R+I)		; 22-24  		avail 4,9 */ \
	xmm9 = xmm13;			/* Copy cosine/sine					avail 4 */ \
	xmm2 *= xmm1;			/* B2 = I2 * cosine/sine	; 22-26 */ \
\
	xmm4 = xmm6;			/*#2 Copy I1-.5I2-.5I3					avail none */ \
	xmm6 += xmm7;			/*#2 Final I2 (I+R)		; 23-25 */ \
	xmm9 *= xmm12;			/* B3 = I3 * cosine/sine	; 23-27 */ \
\
	xmm4 -= xmm7;			/*#2 Final I3 (I-R)		; 24-26			avail 7 */ \
	xmm7 = xmm13;			/* Copy cosine/sine					avail none */ \
	xmm13 *= xmm8;			/*#2 A2 = R2 * cosine/sine	; 24-28 */ \
\
	xmm0 -= xmm1;			/* A2 = A2 - I2			; 25-27			avail 1 */ \
	xmm1 = xmm7;			/* Copy cosine/sine					avail none */ \
	xmm7 *= xmm11;			/*#2 A3 = R3 * cosine/sine	; 25-29 */ \
\
	xmm3 += xmm12;			/* A3 = A3 + I3			; 26-28			avail 12 */ \
	xmm12 = xmm1;			/* Copy cosine/sine					avail none */ \
	xmm1 *= xmm6;			/*#2 B2 = I2 * cosine/sine	; 26-30 */ \
\
	xmm2 += xmm10;			/* B2 = B2 + R2			; 27-29			avail 10 */ \
	xmm12 *= xmm4;			/*#2 B3 = I3 * cosine/sine	; 27-31 */ \
\
	xmm9 -= xmm5;			/* B3 = B3 - R3			; 28-30			avail 5,10 */ \
	xmm10 = xptr(screg); \
	xmm0 *= xmm10;			/* A2 = A2 * sine (new R2)	; 28-32 */ \
\
	xmm13 -= xmm6;			/*#2 A2 = A2 - I2		; 29-31 */ \
	xmm3 *= xmm10;			/* A3 = A3 * sine (new R3)	; 29-33 */ \
	xmm7 += xmm4;			/*#2 A3 = A3 + I3		; 30-32 */ \
	xmm2 *= xmm10;			/* B2 = B2 * sine (new I2)	; 30-34 */ \
	xmm1 += xmm8;			/*#2 B2 = B2 + R2		; 31-33 */ \
	xmm9 *= xmm10;			/* B3 = B3 * sine (new I3)	; 31-35 */ \
	xmm12 -= xmm11;			/*#2 B3 = B3 - R3		; 32-34 */ \
	xmm13 *= xmm10;			/*#2 A2 = A2 * sine (new R2)	; 32-36 */ \
	xmm7 *= xmm10;			/*#2 A3 = A3 * sine (new R3)	; 33-37 */ \
	xptr(srcreg+2*d1) = xmm0;		/* Save R2			; 33 */ \
	xmm1 *= xmm10;			/*#2 B2 = B2 * sine (new I2)	; 34-38 */ \
	xptr(srcreg+2*d1+32) = xmm3;		/* Save R3			; 34 */ \
	xmm12 *= xmm10;			/*#2 B3 = B3 * sine (new I3)	; 35-39 */ \
\
	xptr(srcreg+2*d1+16) = xmm2;		/* Save I2 */ \
	xptr(srcreg+2*d1+48) = xmm9;		/* Save I3 */ \
	xptr(srcreg+32) = xmm13;		/* Save R2#2 */ \
	xptr(srcreg+d1) = xmm7;		/* Save R3#2 */ \
	xptr(srcreg+48) = xmm1;		/* Save I2#2 */ \
	xptr(srcreg+d1+16) = xmm12;		/* Save I3#2 */ \
	srcreg += srcinc; \
\
\

#define r3_f3cl_three_complex_djbfft_common_preload \
	xmm14 = 0.5; \
	xmm15 = CONST2_P866; \
\
\

#define r3_f3cl_three_complex_djbfft_common(srcreg,srcoff,srcinc,d1,screg,scoff) \
	xmm0 = xptr(srcreg+srcoff+2*d1);	/* R2 */ \
	xmm1 = xptr(srcreg+srcoff+d1+16);	/* R3 */ \
	xmm2 = xmm0;			/* Copy R2 */ \
	xmm0 += xmm1;			/* R2 + R3			; 1-3 */ \
\
	xmm3 = xptr(srcreg+srcoff+2*d1+32);	/* I2 */ \
	xmm4 = xptr(srcreg+srcoff+d1+48);	/* I3 */ \
	xmm5 = xmm3;			/* Copy I2 */ \
	xmm3 += xmm4;			/* I2 + I3			; 2-4 */ \
\
	xmm2 -= xmm1;			/* R2 - R3			; 3-5 */ \
	xmm1 = xmm14;			/* 0.5 */ \
\
	xmm5 -= xmm4;			/* I2 - I3			; 4-6			avail 1,4,7-13 */ \
	xmm1 *= xmm0;			/* 0.5 * (R2 + R3)		; 4-8			avail 4,7-13 */ \
\
	xmm4 = xptr(srcreg+srcoff+16);	/*#2 R2 */ \
	xmm6 = xptr(srcreg+srcoff+2*d1+16);	/*#2 R3 */ \
	xmm7 = xmm4;			/*#2 Copy R2						avail 8-13 */ \
	xmm4 += xmm6;			/*#2 R2 + R3			; 5-7 */ \
	xmm8 = xmm14;			/* 0.5 */ \
	xmm8 *= xmm3;			/* 0.5 * (I2 + I3)		; 5-9			avail 9-13 */ \
\
	xmm9 = xptr(srcreg+srcoff+48);	/*#2 I2 */ \
	xmm10 = xptr(srcreg+srcoff+2*d1+48);	/*#2 I3 */ \
	xmm11 = xmm9;			/*#2 Copy I2 */ \
	xmm9 += xmm10;			/*#2 I2 + I3			; 6-8			avail 12-13 */ \
	xmm2 *= xmm15;			/* 0.866 * (R2 - R3)		; 6-10 */ \
\
	xmm7 -= xmm6;			/*#2 R2 - R3			; 7-9			avail 6,12-13 */ \
	xmm5 *= xmm15;			/* 0.866 * (I2 - I3)		; 7-11 */ \
	xmm6 = xmm14;			/* 0.5 */ \
\
	xmm11 -= xmm10;			/*#2 I2 - I3			; 8-10			avail 6,10,12-13 */ \
	xmm6 *= xmm4;			/*#2 0.5 * (R2 + R3)		; 8-12			avail 10,12-13 */ \
	xmm12 = xmm14;			/* 0.5 */ \
\
	xmm10 = xptr(srcreg+srcoff);		/* R1 */ \
	xmm0 += xmm10;			/* R1 + R2 + R3 (final R1)	; 9-11			avail 12-13 storable 0 */ \
	xmm12 *= xmm9;			/*#2 0.5 * (I2 + I3)		; 9-13			avail 13 storable 0 */ \
\
	xmm10 -= xmm1;			/* (R1-.5R2-.5R3)		; 10-12			avail 1,13 storable 0 */ \
	xmm7 *= xmm15;			/*#2 0.866 * (R2 - R3)		; 10-14 */ \
\
	xmm1 = xptr(srcreg+srcoff+32);	/* I1 */ \
	xmm3 += xmm1;			/* I1 + I2 + I3 (final I1)	; 11-13			avail 13 storable 0,3 */ \
	xmm11 *= xmm15;			/*#2 0.866 * (I2 - I3)		; 11-15 */ \
\
	xmm1 -= xmm8;			/* (I1-.5I2-.5I3)		; 12-14			avail 8,13 storable 0,3 */ \
	xmm8 = xptr(srcreg+srcoff+d1);	/*#2 R1 */ \
	xptr(srcreg) = xmm0;			/* Save R1			; 12			avail 0,13 storable 3 */ \
\
	xmm4 += xmm8;			/*#2 R1 + R2 + R3 (final R1)	; 13-15			avail 0,13 storable 3,4 */ \
	xmm13 = xptr(screg+16);		/* cosine/sine						avail 0 storable 3,4 */ \
\
	xmm8 -= xmm6;			/*#2 (R1-.5R2-.5R3)		; 14-16			avail 0,6 storable 3,4 */ \
	xmm0 = xmm13;			/* Copy cosine/sine					avail 6 storable 3,4 */ \
	xmm6 = xptr(srcreg+srcoff+d1+32);	/*#2 I1 */ \
	xptr(srcreg+32) = xmm3;		/* Save I1			; 14			avail 3 storable 4 */ \
\
	xmm9 += xmm6;			/*#2 I1 + I2 + I3 (final I1)	; 15-17			avail 3 storable 4,9 */ \
	xmm3 = xmm13;			/* Copy cosine/sine					avail none storable 4,9 */ \
\
	xmm6 -= xmm12;			/*#2 (I1-.5I2-.5I3)		; 16-18			avail 12 storable 4,9 */ \
	xptr(srcreg+16) = xmm4;		/*#2 Save R1			; 16			avail 4,12 storable 9 */ \
\
	xmm12 = xmm10;			/* Copy R1-.5R2-.5R3 */ \
	xmm10 -= xmm5;			/* Final R2 (R-I)		; 17-19			avail 4 storable 9 */ \
\
	xmm5 += xmm12;			/* Final R3 (R+I)		; 18-20 */ \
	xptr(srcreg+48) = xmm9;		/*#2 Save I1			; 18			avail 4,9,12 */ \
\
	xmm12 = xmm1;			/* Copy I1-.5I2-.5I3					avail 4,9 */ \
	xmm1 += xmm2;			/* Final I2 (I+R)		; 19-21 */ \
\
	xmm12 -= xmm2;			/* Final I3 (I-R)		; 20-22			avail 2,4,9 */ \
	xmm2 = xmm13;			/* Copy cosine/sine */ \
	xmm0 *= xmm10;			/* A2 = R2 * cosine/sine	; 20-24  		avail 4,9 */ \
\
	xmm4 = xmm8;			/*#2 Copy R1-.5R2-.5R3					avail 9 */ \
	xmm8 -= xmm11;			/*#2 Final R2 (R-I)		; 21-23 */ \
	xmm3 *= xmm5;			/* A3 = R3 * cosine/sine	; 21-25 */ \
\
	xmm11 += xmm4;			/*#2 Final R3 (R+I)		; 22-24  		avail 4,9 */ \
	xmm9 = xptr(screg+scoff+16);		/* Copy cosine/sine					avail 4 */ \
	xmm2 *= xmm1;			/* B2 = I2 * cosine/sine	; 22-26 */ \
\
	xmm4 = xmm6;			/*#2 Copy I1-.5I2-.5I3					avail none */ \
	xmm6 += xmm7;			/*#2 Final I2 (I+R)		; 23-25 */ \
	xmm13 *= xmm12;			/* B3 = I3 * cosine/sine	; 23-27 */ \
\
	xmm4 -= xmm7;			/*#2 Final I3 (I-R)		; 24-26			avail 7 */ \
	xmm7 = xmm9;			/* Copy cosine/sine					avail none */ \
	xmm9 *= xmm8;			/*#2 A2 = R2 * cosine/sine	; 24-28 */ \
\
	xmm0 -= xmm1;			/* A2 = A2 - I2			; 25-27			avail 1 */ \
	xmm1 = xmm7;			/* Copy cosine/sine					avail none */ \
	xmm7 *= xmm11;			/*#2 A3 = R3 * cosine/sine	; 25-29 */ \
\
	xmm3 += xmm12;			/* A3 = A3 + I3			; 26-28			avail 12 */ \
	xmm12 = xmm1;			/* Copy cosine/sine					avail none */ \
	xmm1 *= xmm6;			/*#2 B2 = I2 * cosine/sine	; 26-30 */ \
\
	xmm2 += xmm10;			/* B2 = B2 + R2			; 27-29			avail 10 */ \
	xmm12 *= xmm4;			/*#2 B3 = I3 * cosine/sine	; 27-31 */ \
\
	xmm13 -= xmm5;			/* B3 = B3 - R3			; 28-30			avail 5,10 */ \
	xmm10 = xptr(screg); \
	xmm0 *= xmm10;			/* A2 = A2 * sine (new R2)	; 28-32 */ \
\
	xmm9 -= xmm6;			/*#2 A2 = A2 - I2		; 29-31 */ \
	xmm3 *= xmm10;			/* A3 = A3 * sine (new R3)	; 29-33 */ \
	xmm7 += xmm4;			/*#2 A3 = A3 + I3		; 30-32 */ \
	xmm2 *= xmm10;			/* B2 = B2 * sine (new I2)	; 30-34 */ \
	xmm1 += xmm8;			/*#2 B2 = B2 + R2		; 31-33 */ \
	xmm13 *= xmm10;			/* B3 = B3 * sine (new I3)	; 31-35 */ \
	xmm10 = xptr(screg+scoff); \
	xmm12 -= xmm11;			/*#2 B3 = B3 - R3		; 32-34 */ \
	xmm9 *= xmm10;			/*#2 A2 = A2 * sine (new R2)	; 32-36 */ \
	xmm7 *= xmm10;			/*#2 A3 = A3 * sine (new R3)	; 33-37 */ \
	xptr(srcreg+d1) = xmm0;		/* Save R2			; 33 */ \
	xmm1 *= xmm10;			/*#2 B2 = B2 * sine (new I2)	; 34-38 */ \
	xptr(srcreg+2*d1) = xmm3;		/* Save R3			; 34 */ \
	xmm12 *= xmm10;			/*#2 B3 = B3 * sine (new I3)	; 35-39 */ \
\
	xptr(srcreg+d1+32) = xmm2;		/* Save I2 */ \
	xptr(srcreg+2*d1+32) = xmm13;		/* Save I3 */ \
	xptr(srcreg+d1+16) = xmm9;		/* Save R2#2 */ \
	xptr(srcreg+2*d1+16) = xmm7;		/* Save R3#2 */ \
	xptr(srcreg+d1+48) = xmm1;		/* Save I2#2 */ \
	xptr(srcreg+2*d1+48) = xmm12;		/* Save I3#2 */ \
	srcreg += srcinc; \

#endif
#endif

/* The Pentium-4 and AMD K10 64-bit version.  This is just the 32-bit version */ \
/* with constants preloaded. */ \

#if 0 // (@INSTR(,%xarch,<P4>) NE 0) OR (@INSTR(,%xarch,<K10>) NE 0);
#if 0 // IFDEF X86_64;

/* Identical to the 32-bit version but uses registers for HALF and P866.  Very likely */ \
/* this could be optimized further. */ \
\

#define r3_x3c_djbfft_partial_mem_preload \
	xmm14 = 0.5; \
	xmm15 = CONST2_P866; \
\
\

#define r3_x3c_djbfft_partial_mem(r1,r2,r3,r4,r5,r6,t1,t2,memr1,memr2,memi1,memi2,screg,pre1,pre2) \
	r2 = memr2; \
	t1 = r3; \
	r3 += r2;			/* R2 + R3 */ \
	r5 = memi2; \
	t2 = r6; \
	r6 += r5;			/* I2 + I3 */ \
	r2 -= t1;			/* R2 - R3 */ \
	t1 = xmm14; \
	t1 *= r3;			/* 0.5 * (R2 + R3) */ \
	r5 -= t2;			/* I2 - I3 */ \
	t2 = xmm14; \
	t2 *= r6;			/* 0.5 * (I2 + I3) */ \
	r1 = memr1; \
	r3 += r1;			/* R1 + R2 + R3 (final R1) */ \
	r2 *= xmm15;		/* 0.866 * (R2 - R3) */ \
	r1 -= t1;			/* (R1-.5R2-.5R3) */ \
	r5 *= xmm15;		/* 0.866 * (I2 - I3) */ \
	r4 = memi1; \
	r6 += r4;			/* I1 + I2 + I3 (final I1) */ \
	r4 -= t2;			/* (I1-.5I2-.5I3) */ \
\
	xprefetchw(u8ptr(pre1)); \
\
	t1 = r1; \
	r1 -= r5;			/* Final R2 */ \
	t2 = r4; \
	r4 -= r2;			/* Final I3 */ \
	r5 += t1;			/* Final R3 */ \
	r2 += t2;			/* Final I2 */ \
\
	if 0 // pre2 NE 0; \
	xprefetchw(u8ptr(pre1+pre2)); \
	endif \
\
	t1 = r1;			/* Copy R2 */ \
	t2 = xptr(screg+16); \
	r1 *= t2;			/* A2 = R2 * cosine/sine */ \
	r1 -= r2;			/* A2 = A2 - I2 */ \
	r2 *= t2;			/* B2 = I2 * cosine/sine */ \
	r2 += t1;			/* B2 = B2 + R2 */ \
\
	t1 = r5;			/* Copy R3 */ \
	r5 *= t2;			/* A3 = R3 * cosine/sine */ \
	r5 += r4;			/* A3 = A3 + I3 */ \
	r4 *= t2;			/* B3 = I3 * cosine/sine */ \
	r4 -= t1;			/* B3 = B3 - R3 */ \
\
	t2 = xptr(screg); \
	r1 *= t2;			/* A2 = A2 * sine (new R2) */ \
	r2 *= t2;			/* B2 = B2 * sine (new I2) */ \
	r5 *= t2;			/* A3 = A3 * sine (new R3) */ \
	r4 *= t2;			/* B3 = B3 * sine (new I3) */ \

#endif
#endif

/* 64-bit AMD K8 version.  Same as the 32-bit version -- can be improved */

#if 0 // (@INSTR(,%xarch,<K8>) NE 0);
#if 0 // IFDEF X86_64;

#define r3_x3c_djbfft_partial_mem_preload \
	xmm14 = CONST2_P866; \
	xmm15 = XMM_TWO; \
\
\

#define r3_x3c_djbfft_partial_mem(r1,r2,r3,r4,r5,r6,t1,t2,memr1,memr2,memi1,memi2,screg,pre1,pre2) \
	r2 = memr2; \
	r2 -= r3;			/* R2 - R3 */ \
	r3 += memr2;		/* R2 + R3 */ \
	r5 = memi2; \
	r5 -= r6;			/* I2 - I3 */ \
	r6 += memi2;		/* I2 + I3 */ \
\
	r1 = 0.5; \
	r1 *= r3;			/* 0.5 * (R2 + R3) */ \
	r4 = 0.5; \
	r4 *= r6;			/* 0.5 * (I2 + I3) */ \
\
	t1 = memr1; \
	r3 += t1;			/* R1 + R2 + R3 (final R1) */ \
\
	r2 *= xmm14;		/* 0.866 * (R2 - R3) */ \
	t1 -= r1;			/* (R1-.5R2-.5R3) */ \
	r5 *= xmm14;		/* 0.866 * (I2 - I3) */ \
	t2 = memi1; \
	r6 += t2;			/* I1 + I2 + I3 (final I1) */ \
	t2 -= r4;			/* (I1-.5I2-.5I3) */ \
\
	xprefetchw(u8ptr(pre1)); \
\
	t1 -= r5;			/* Final R2 */ \
	r5 *= xmm15; \
	t2 -= r2;			/* Final I3 */ \
	r2 *= xmm15; \
	r5 += t1;			/* Final R3 */ \
	r2 += t2;			/* Final I2 */ \
\
	if 0 // pre2 NE 0; \
	xprefetchw(u8ptr(pre1+pre2)); \
	endif \
\
	r1 = xptr(screg+16); \
	r1 *= t1;			/* A2 = R2 * cosine/sine */ \
\
	r4 = xptr(screg+16); \
	r4 *= t2;			/* B3 = I3 * cosine/sine */ \
\
	r1 -= r2;			/* A2 = A2 - I2 */ \
	r2 *= xptr(screg+16);		/* B2 = I2 * cosine/sine */ \
\
	r4 -= r5;			/* B3 = B3 - R3 */ \
	r5 *= xptr(screg+16);		/* A3 = R3 * cosine/sine */ \
\
	r2 += t1;			/* B2 = B2 + R2 */ \
	r5 += t2;			/* A3 = A3 + I3 */ \
\
	t2 = xptr(screg); \
	r1 *= t2;			/* A2 = A2 * sine (new R2) */ \
	r4 *= t2;			/* B3 = B3 * sine (new I3) */ \
	r2 *= t2;			/* B2 = B2 * sine (new I2) */ \
	r5 *= t2;			/* A3 = A3 * sine (new R3) */ \


#endif
#endif

/* */ \
/* ************************************* three-complex-djbunfft variants ****************************************** */ \
/* */ \
\

#define r3_x3cl_three_complex_djbunfft_preload \
	r3_x3c_djbunfft_partial_mem_preload; \
\
\

#define r3_x3cl_three_complex_djbunfft(srcreg,srcinc,d1,screg_a,scoff) { \
	vec2f64 *sd0 = (vec2f64*)(srcreg); \
	vec2f64 *sd1 = (vec2f64*)(srcreg+(d1)); \
	vec2f64 *sd2 = (vec2f64*)(srcreg+2*(d1)); \
	char *screg = (char*)(screg_a); /* improves clang compiler */ \
	xmm2 = sd1[0];	/* R2 */ \
	xmm3 = sd1[2];	/* I2 */ \
	r3_x3c_djbunfft_partial_mem(xmm0,xmm1,xmm2,xmm3,xmm4,xmm5,xmm6,xmm7, sd0[0],sd0[2],sd2[0],sd2[2], screg,0,srcreg+srcinc,d1); \
	xmm6 = sd1[1];	/* R2 */ \
	xmm7 = sd1[3];	/* I2 */ \
	sd0[0] = xmm4;		/* Save R1 */ \
	sd0[2] = xmm5;	/* Save I1 */ \
	sd2[0] = xmm3;	/* Save R2 */ \
	sd2[2] = xmm1;	/* Save I2 */ \
	sd1[1] = xmm0;	/* Save R3 */ \
	sd1[3] = xmm2;	/* Save I3 */ \
	r3_x3c_djbunfft_partial_mem(xmm0,xmm1,xmm6,xmm7,xmm4,xmm5,xmm2,xmm3, sd0[1],sd0[3],sd2[1],sd2[3], screg,scoff,srcreg+srcinc+2*d1,0); \
	sd1[0] = xmm4;	/* Save R1 */ \
	sd1[2] = xmm5;	/* Save I1 */ \
	sd0[1] = xmm7;	/* Save R2 */ \
	sd0[3] = xmm1;	/* Save I2 */ \
	sd2[1] = xmm0;	/* Save R3 */ \
	sd2[3] = xmm6;	/* Save I3 */ \
	srcreg += srcinc; \
}

/* Do a 3-complex inverse FFT.  The input values are R1+R2i, R3+R4i, R5+R6i */ \
/* First we apply twiddle factors to 2 of the 3 input numbers. */ \
/* A 3-complex inverse FFT is: */ \
/* Res1:  (R1+R2+R3) + (I1+I2+I3)i */ \
/* Res2:  (R1-.5R2+.866I2-.5R3-.866I3) + (I1-.5I2-.866R2-.5I3+.866R3)i */ \
/* Res3:  (R1-.5R2-.866I2-.5R3+.866I3) + (I1-.5I2+.866R2-.5I3-.866R3)i */ \
\

#define r3_x3c_djbunfft_partial_mem_preload \
\
\

#define r3_x3c_djbunfft_partial_mem(r1,r2,r3,r4,r5,r6,t1,t2,memr1,memi1,memr3,memi3,screg,scoff,pre1,pre2) \
	t1 = r3;			/* Copy R2 */ \
	t2 = xptr(screg+scoff+16); \
	r5 = memr3; \
	r3 = r3 * t2 + r4;	/* A2 = R2 * cosine/sine + I2 */ \
	r6 = memi3; \
	r5 = r5 * t2 - r6;	/* A3 = R3 * cosine/sine - I3 */ \
	r4 = r4 * t2 - t1;	/* B2 = I2 * cosine/sine - R2 */ \
	r6 = r6 * t2 + memr3;	/* B3 = I3 * cosine/sine + R3 */ \
	t2 = xptr(screg+scoff); \
	r3 *= t2;			/* A2 = A2 * sine (final R2) */ \
	r5 *= t2;			/* A3 = A3 * sine (final R3) */ \
	r4 *= t2;			/* B2 = B2 * sine (final I2) */ \
	r6 *= t2;			/* B3 = B3 * sine (final I3) */ \
\
	xprefetchw(u8ptr(pre1)); \
\
	t1 = r3; \
	r3 -= r5;			/* R2 - R3 */ \
	r5 += t1;			/* R2 + R3 */ \
	r3 *= CONST2_P866;	/* 0.866 * (R2 - R3) */ \
	t2 = r4; \
	r4 -= r6;			/* I2 - I3 */ \
\
	/*if 0 // pre2 NE 0; \
	xprefetchw(u8ptr(pre1+pre2)); \
	endif */ \
\
	r6 += t2;			/* I2 + I3 */ \
	r4 *= CONST2_P866;	/* 0.866 * (I2 - I3) */ \
	t1 = r5;			/* 0.5 * (R2 + R3) */ \
	r1 = memr1; \
	r5 += r1;			/* R1 + R2 + R3 (final R1) */ \
	r1 -= t1 * 0.5;		/* (R1-.5R2-.5R3) */ \
	t1 = r6;			/* 0.5 * (I2 + I3) */ \
	r2 = memi1; \
	r6 += r2;			/* I1 + I2 + I3 (final I1) */ \
	r2 -= t1 * 0.5;		/* (I1-.5I2-.5I3) */ \
	t1 = r1; \
	r1 -= r4;			/* Final R3 */ \
	r4 += t1;			/* Final R2 */ \
	t2 = r2; \
	r2 -= r3;			/* Final I2 */ \
	r3 += t2;			/* Final I3 */ \
\
\
/* AMD K8 version of some of the above macros.  Could be improved slightly by eliminating the  */ \
/* first xcopy.  Half the time this is called the preloaded r3 value is still in memory. */ \
/* Or we eliminate the initial xcopy (do a cosine/sine load instead) but it will */ \
/* change the output registers. */ \

#if 0 // (@INSTR(,%xarch,<K8>) NE 0);

#define r3_x3c_djbunfft_partial_mem(r1,r2,r3,r4,r5,r6,t1,t2,memr1,memi1,memr3,memi3,screg,scoff,pre1,pre2) \
	t1 = r3;			/* Copy R2 */ \
	t2 = xptr(screg+scoff+16); \
	r3 *= t2;			/* A2 = R2 * cosine/sine */ \
	r5 = memr3; \
	r5 *= t2;			/* A3 = R3 * cosine/sine */ \
	r3 += r4;			/* A2 = A2 + I2 */ \
	r6 = memi3; \
	r5 -= r6;			/* A3 = A3 - I3 */ \
	r4 *= t2;			/* B2 = I2 * cosine/sine */ \
	r6 *= t2;			/* B3 = I3 * cosine/sine */ \
	r4 -= t1;			/* B2 = B2 - R2 */ \
	r6 += memr3;		/* B3 = B3 + R3 */ \
	t2 = xptr(screg+scoff); \
	r3 *= t2;			/* A2 = A2 * sine (final R2) */ \
	r5 *= t2;			/* A3 = A3 * sine (final R3) */ \
	r4 *= t2;			/* B2 = B2 * sine (final I2) */ \
	r6 *= t2;			/* B3 = B3 * sine (final I3) */ \
\
	xprefetchw(u8ptr(pre1)); \
\
	r3 -= r5;			/* R2 - R3 */ \
	r5 += r5; /* multwo */; \
	r4 -= r6;			/* I2 - I3 */ \
	r6 += r6; /* multwo */; \
	r5 += r3;			/* R2 + R3 */ \
	r3 *= CONST2_P866;		/* 0.866 * (R2 - R3) */ \
	r6 += r4;			/* I2 + I3 */ \
\
	if 0 // pre2 NE 0; \
	xprefetchw(u8ptr(pre1+pre2)); \
	endif \
\
	r4 *= CONST2_P866;		/* 0.866 * (I2 - I3) */ \
	t1 = 0.5; \
	t1 *= r5;			/* 0.5 * (R2 + R3) */ \
	r1 = memr1; \
	r5 += r1;			/* R1 + R2 + R3 (final R1) */ \
	r1 -= t1;			/* (R1-.5R2-.5R3) */ \
	t1 = 0.5; \
	t1 *= r6;			/* 0.5 * (I2 + I3) */ \
	r2 = memi1; \
	r6 += r2;			/* I1 + I2 + I3 (final I1) */ \
	r2 -= t1;			/* (I1-.5I2-.5I3) */ \
\
	r1 -= r4;			/* Final R3 */ \
	r4 += f4; /* multwo */; \
	r2 -= r3;			/* Final I2 */ \
	r3 += r3; /* multwo */; \
	r4 += r1;			/* Final R2 */ \
	r3 += r2;			/* Final I3 */ \


#endif

/* 64-bit version of some of the above macros - use the extra registers */ \

#if 0 // (@INSTR(,%xarch,<CORE>) NE 0) OR (@INSTR(,%xarch,<K10>) NE 0) OR (@INSTR(,%xarch,<BLEND>) NE 0);
#if 0 // IFDEF X86_64;

/* Core 2 optimal is 32 clocks.  Currently at 39 clocks. */ \

#define r3_x3cl_three_complex_djbunfft_preload \
	xmm14 = 0.5; \
	xmm15 = CONST2_P866; \
\
\

#define r3_x3cl_three_complex_djbunfft(srcreg,srcinc,d1,screg,scoff) \
	xmm0 = xptr(srcreg+d1);	/* R2 */ \
	xmm1 = xmm0;		/* Copy R2 */ \
	xmm13 = xptr(screg+16);	/* cosine/sine */ \
	xmm0 *= xmm13;		/* A2 = R2 * cosine/sine		; 1-5 */ \
\
	xmm2 = xptr(srcreg+2*d1);	/* R3 */ \
	xmm3 = xmm2;		/* Copy R3 */ \
	xmm2 *= xmm13;		/* A3 = R3 * cosine/sine		; 2-6 */ \
\
	xmm4 = xptr(srcreg+d1+32);	/* I2 */ \
	xmm5 = xmm4;		/* Copy I2 */ \
	xmm4 *= xmm13;		/* B2 = I2 * cosine/sine		; 3-7 */ \
\
	xmm6 = xptr(srcreg+2*d1+32);	/* I3 */ \
	xmm13 *= xmm6;		/* B3 = I3 * cosine/sine		; 4-8 */ \
\
	xmm7 = xptr(srcreg+d1+16);	/*#2 R2 */ \
	xmm8 = xmm7;		/*#2 Copy R2 */ \
	xmm12 = xptr(screg+scoff+16);	/*#2 cosine/sine */ \
	xmm7 *= xmm12;		/*#2 A2 = R2 * cosine/sine		; 5-9 */ \
\
	xmm0 += xmm5;		/* A2 = A2 + I2				; 6-8			avail 5,9-11 */ \
	xmm9 = xptr(srcreg+2*d1+16);	/*#2 R3 */ \
	xmm10 = xmm9;		/*#2 Copy R3							avail 5,11 */ \
	xmm9 *= xmm12;		/*#2 A3 = R3 * cosine/sine		; 6-10 */ \
\
	xmm2 -= xmm6;		/* A3 = A3 - I3				; 7-9			avail 5,6,11 */ \
	xmm5 = xptr(srcreg+d1+48);	/*#2 I2 */ \
	xmm6 = xmm5;		/*#2 Copy I2							avail 11 */ \
	xmm5 *= xmm12;		/*#2 B2 = I2 * cosine/sine		; 7-11 */ \
\
	xmm4 -= xmm1;		/* B2 = B2 - R2				; 8-10			avail 1,11 */ \
	xmm11 = xptr(srcreg+2*d1+48);	/*#2 I3								avail 1 */ \
	xmm12 *= xmm11;		/*#2 B3 = I3 * cosine/sine		; 8-12 */ \
\
	xmm13 += xmm3;		/* B3 = B3 + R3				; 9-11			avail 1,3 */ \
	xmm1 = xptr(screg);		/* sine								avail 3 */ \
\
	xmm0 *= xmm1;		/* A2 = A2 * sine (new R2)		; 9-13 */ \
	xmm7 += xmm6;		/*#2 A2 = A2 + I2 (new R2 / sine)	; 10-12			avail 3,6 */ \
	xmm2 *= xmm1;		/* A3 = A3 * sine (new R3)		; 10-14 */ \
\
	xprefetchw(u8ptr(srcreg+srcinc)); \
\
	xmm9 -= xmm11;		/*#2 A3 = A3 - I3 (new R3 / sine)	; 11-13			avail 3,6,11 */ \
	xmm4 *= xmm1;		/* B2 = B2 * sine (new I2)		; 11-15 */ \
	xmm5 -= xmm8;		/*#2 B2 = B2 - R2 (new I2 / sine)	; 12-14			avail 3,6,8,11 */ \
	xmm13 *= xmm1;		/* B3 = B3 * sine (new I3)		; 12-16			avail 1,3,6,8,11 */ \
	xmm12 += xmm10;		/*#2 B3 = B3 + R3 (new I3 / sine)	; 13-15			avail 1,3,6,8,10,11 */ \
\
	xprefetchw(u8ptr(srcreg+srcinc+d1)); \
\
	xmm3 = xmm7;		/*#2 Copy R2 / sine */ \
	xmm7 += xmm9;		/*#2 (R2 + R3) / sine			; 14-16			avail 1,6,8,10,11 */ \
	xmm8 = xmm0;		/* Copy R2 */ \
	xmm0 += xmm2;		/* R2 + R3				; 15-17			avail 1,6,10,11 */ \
	xmm6 = xmm5;		/*#2 Copy I2 / sine */ \
\
	xprefetchw(u8ptr(srcreg+srcinc+2*d1)); \
\
	xmm5 += xmm12;		/*#2 (I2 + I3) / sine			; 16-18			avail 1,10,11 */ \
	xmm10 = xmm4;		/* Copy I2 */ \
\
	xmm4 += xmm13;		/* I2 + I3				; 17-19			avail 1,11 */ \
	xmm11 = xptr(screg+scoff);	/*#2 sine							avail 1 */ \
	xmm7 *= xmm11;		/*#2 R2 + R3				; 17-21 */ \
\
	xmm3 -= xmm9;		/*#2 (R2 - R3) /sine			; 18-20			avail 1,9 */ \
	xmm1 = xmm14;		/* 0.5 */ \
	xmm1 *= xmm0;		/* 0.5 * (R2 + R3)			; 18-22			avail 9 */ \
\
	xmm6 -= xmm12;		/*#2 (I2 - I3) / sine			; 19-21			avail 9,12 */ \
	xmm5 *= xmm11;		/*#2 I2 + I3				; 19-23 */ \
\
	xmm8 -= xmm2;		/* R2 - R3				; 20-22			avail 2,9,12 */ \
	xmm12 = xmm14;		/* 0.5 */ \
	xmm12 *= xmm4;		/* 0.5 * (I2 + I3)			; 20-24			avail 2,9 */ \
\
	xmm10 -= xmm13;		/* I2 - I3				; 21-23			avail 2,9,13 */ \
	xmm3 *= xmm11;		/*#2 R2 - R3				; 21-25 */ \
\
	xmm9 = xptr(srcreg);		/* R1 */ \
	xmm0 += xmm9;		/* R1 + R2 + R3 (final R1)		; 22-24			avail 2,13 storable 0 */ \
	xmm6 *= xmm11;		/*#2 I2 - I3				; 22-26			avail 2,11,13 storable 0 */ \
\
	xmm9 -= xmm1;		/* (R1-.5R2-.5R3)			; 23-25			avail 1,2,11,13 storable 0 */ \
	xmm13 = xmm14;		/*#2 0.5 */ \
	xmm13 *= xmm7;		/*#2 0.5 * (R2 + R3)			; 23-27			avail 1,2,11 storable 0 */ \
\
	xmm2 = xptr(srcreg+32);	/* I1 */ \
	xmm4 += xmm2;		/* I1 + I2 + I3 (final I1)		; 24-26			avail 1,11 storable 0,4 */ \
	xmm11 = xmm14;		/*#2 0.5 */ \
	xmm11 *= xmm5;		/*#2 0.5 * (I2 + I3)			; 24-28			avail 1 storable 0,4 */ \
\
	xmm2 -= xmm12;		/* (I1-.5I2-.5I3)			; 25-27			avail 1,12 storable 0,4 */ \
	xmm10 *= xmm15;		/* 0.866 * (I2 - I3)			; 25-29 */ \
\
	xmm1 = xptr(srcreg+16);	/*#2 R1 */ \
	xmm7 += xmm1;		/*#2 R1 + R2 + R3 (final R1)		; 26-28			avail 12 storable 0,4,7 */ \
	xmm8 *= xmm15;		/* 0.866 * (R2 - R3)			; 26-30 */ \
\
	xmm12 = xptr(srcreg+48);	/*#2 I1 */ \
	xmm5 += xmm12;		/*#2 I1 + I2 + I3 (final I1)		; 27-29			avail none storable 0.4.5.7 */ \
	xmm6 *= xmm15;		/*#2 0.866 * (I2 - I3)			; 27-31 */ \
\
	xmm1 -= xmm13;		/*#2 (R1-.5R2-.5R3)			; 28-30			avail 13 storable 0.4.5.7 */ \
	xmm3 *= xmm15;		/*#2 0.866 * (R2 - R3)			; 28-32 */ \
	xptr(srcreg) = xmm0;		/* Save R1				; 25 */ \
\
	xmm12 -= xmm11;		/*#2 (I1-.5I2-.5I3)			; 29-31			avail 11,13 */ \
	xptr(srcreg+32) = xmm4;	/* Save I1				; 27 */ \
\
	xmm13 = xmm9;		/* Copy R1-.5R2-.5R3 */ \
	xmm9 -= xmm10;		/* Final R3 (R-I)			; 30-32 */ \
	xptr(srcreg+d1) = xmm7;	/*#2 Save R1				; 29 */ \
	xmm10 += xmm13;		/* Final R2 (R+I)			; 31-33 */ \
	xptr(srcreg+d1+32) = xmm5;	/*#2 Save I1				; 30 */ \
\
	xmm11 = xmm2;		/* Copy I1-.5I2-.5I3 */ \
	xmm2 -= xmm8;		/* Final I2 (I-R)			; 32-34 */ \
	xmm8 += xmm11;		/* Final I3 (I+R)			; 33-35 */ \
	xptr(srcreg+d1+16) = xmm9;	/* Save R3				; 33 */ \
\
	xmm0 = xmm1;		/*#2 Copy R1-.5R2-.5R3 */ \
	xmm1 -= xmm6;		/*#2 Final R3 (R-I)			; 34-36 */ \
	xptr(srcreg+2*d1) = xmm10;	/* Save R2				; 34 */ \
	xmm6 += xmm0;		/*#2 Final R2 (R+I)			; 35-37 */ \
	xptr(srcreg+2*d1+32) = xmm2;	/* Save I2				; 35 */ \
\
	xmm4 = xmm12;		/*#2 Copy I1-.5I2-.5I3 */ \
	xmm12 -= xmm3;		/*#2 Final I2 (I-R)			; 36-38 */ \
	xptr(srcreg+d1+48) = xmm8;	/* Save I3				; 36 */ \
	xmm3 += xmm4;		/*#2 Final I3 (I+R)			; 37-39 */ \
\
	xptr(srcreg+2*d1+16) = xmm1;	/*#2 Save R3				; 37 */ \
	xptr(srcreg+16) = xmm6;	/*#2 Save R2 */ \
	xptr(srcreg+48) = xmm12;	/*#2 Save I2 */ \
	xptr(srcreg+2*d1+48) = xmm3;	/*#2 Save I3 */ \
	srcreg += srcinc; \

#endif
#endif

/* The Pentium-4 64-bit version */ \

#if 0 // @INSTR(,%xarch,<P4>) EQ 1;
#if 0 // IFDEF X86_64;

/* Same as 32-bit version but HALF and P866 in registers.  Very likely this */ \
/* could be optimized further. */ \
\

#define r3_x3c_djbunfft_partial_mem_preload \
	xmm14 = 0.5; \
	xmm15 = CONST2_P866; \
\
\

#define r3_x3c_djbunfft_partial_mem(r1,r2,r3,r4,r5,r6,t1,t2,memr1,memi1,memr3,memi3,screg,scoff,pre1,pre2) \
	t1 = r3;			/* Copy R2 */ \
	t2 = xptr(screg+scoff+16); \
	r3 *= t2;			/* A2 = R2 * cosine/sine */ \
	r5 = memr3; \
	r5 *= t2;			/* A3 = R3 * cosine/sine */ \
	r3 += r4;			/* A2 = A2 + I2 */ \
	r6 = memi3; \
	r5 -= r6;			/* A3 = A3 - I3 */ \
	r4 *= t2;			/* B2 = I2 * cosine/sine */ \
	r6 *= t2;			/* B3 = I3 * cosine/sine */ \
	r4 -= t1;			/* B2 = B2 - R2 */ \
	r6 += memr3;		/* B3 = B3 + R3 */ \
	t2 = xptr(screg+scoff); \
	r3 *= t2;			/* A2 = A2 * sine (final R2) */ \
	r5 *= t2;			/* A3 = A3 * sine (final R3) */ \
	r4 *= t2;			/* B2 = B2 * sine (final I2) */ \
	r6 *= t2;			/* B3 = B3 * sine (final I3) */ \
\
	t1 = r3; \
	r3 -= r5;			/* R2 - R3 */ \
	r5 += t1;			/* R2 + R3 */ \
	r3 *= xmm15;		/* 0.866 * (R2 - R3) */ \
	t2 = r4; \
	r4 -= r6;			/* I2 - I3 */ \
	r6 += t2;			/* I2 + I3 */ \
	r4 *= xmm15;		/* 0.866 * (I2 - I3) */ \
	t1 = xmm14; \
	t1 *= r5;			/* 0.5 * (R2 + R3) */ \
	r1 = memr1; \
	r5 += r1;			/* R1 + R2 + R3 (final R1) */ \
	r1 -= t1;			/* (R1-.5R2-.5R3) */ \
	t1 = xmm14; \
	t1 *= r6;			/* 0.5 * (I2 + I3) */ \
	r2 = memi1; \
	r6 += r2;			/* I1 + I2 + I3 (final I1) */ \
	r2 -= t1;			/* (I1-.5I2-.5I3) */ \
	t1 = r1; \
	r1 -= r4;			/* Final R3 */ \
	r4 += t1;			/* Final R2 */ \
	t2 = r2; \
	r2 -= r3;			/* Final I2 */ \
	r3 += t2;			/* Final I3 */ \

#endif
#endif

/* 64-bit AMD K8 version.  same as the 32-bit version but with some preloaded constants. */ \
/* Probably could be improved. */ \

#if 0 // (@INSTR(,%xarch,<K8>) NE 0);
#if 0 // X86_64;


#define r3_x3c_djbunfft_partial_mem_preload \
	xmm14 = CONST2_P866; \
	xmm15 = XMM_TWO; \
\
\

#define r3_x3c_djbunfft_partial_mem(r1,r2,r3,r4,r5,r6,t1,t2,memr1,memi1,memr3,memi3,screg,scoff,pre1,pre2) \
	t1 = r3;			/* Copy R2 */ \
	t2 = xptr(screg+scoff+16); \
	r3 *= t2;			/* A2 = R2 * cosine/sine */ \
	r5 = memr3; \
	r5 *= t2;			/* A3 = R3 * cosine/sine */ \
	r3 += r4;			/* A2 = A2 + I2 */ \
	r6 = memi3; \
	r5 -= r6;			/* A3 = A3 - I3 */ \
	r4 *= t2;			/* B2 = I2 * cosine/sine */ \
	r6 *= t2;			/* B3 = I3 * cosine/sine */ \
	r4 -= t1;			/* B2 = B2 - R2 */ \
	r6 += memr3;		/* B3 = B3 + R3 */ \
	t2 = xptr(screg+scoff); \
	r3 *= t2;			/* A2 = A2 * sine (final R2) */ \
	r5 *= t2;			/* A3 = A3 * sine (final R3) */ \
	r4 *= t2;			/* B2 = B2 * sine (final I2) */ \
	r6 *= t2;			/* B3 = B3 * sine (final I3) */ \
\
	xprefetchw(u8ptr(pre1)); \
\
	r3 -= r5;			/* R2 - R3 */ \
	r5 *= xmm15; \
	r4 -= r6;			/* I2 - I3 */ \
	r6 *= xmm15; \
	r5 += r3;			/* R2 + R3 */ \
	r3 *= xmm14;		/* 0.866 * (R2 - R3) */ \
	r6 += r4;			/* I2 + I3 */ \
\
	if 0 // pre2 NE 0; \
	xprefetchw(u8ptr(pre1+pre2)); \
	endif \
\
	r4 *= xmm14;		/* 0.866 * (I2 - I3) */ \
	t1 = 0.5; \
	t1 *= r5;			/* 0.5 * (R2 + R3) */ \
	r1 = memr1; \
	r5 += r1;			/* R1 + R2 + R3 (final R1) */ \
	r1 -= t1;			/* (R1-.5R2-.5R3) */ \
	t1 = 0.5; \
	t1 *= r6;			/* 0.5 * (I2 + I3) */ \
	r2 = memi1; \
	r6 += r2;			/* I1 + I2 + I3 (final I1) */ \
	r2 -= t1;			/* (I1-.5I2-.5I3) */ \
\
	r1 -= r4;			/* Final R3 */ \
	r4 *= xmm15; \
	r2 -= r3;			/* Final I2 */ \
	r3 *= xmm15; \
	r4 += r1;			/* Final R2 */ \
	r3 += r2;			/* Final I3 */ \


#endif
#endif

/* */ \
/* ************************************* six-reals-fft variants ****************************************** */ \
/* */ \
\

#define r3_x3cl_six_reals_three_complex_djbfft_preload \
	r3_x3c_djbfft_partial_mem_preload; \
\
\

#define r3_x3cl_six_reals_three_complex_djbfft(srcreg,srcinc,d1,screg1,screg2) { \
	vec2f64 *sd0 = (vec2f64*)(srcreg); \
	vec2f64 *sd1 = (vec2f64*)(srcreg+(d1)); \
	vec2f64 *sd2 = (vec2f64*)(srcreg+2*(d1)); \
	xmm2 = sd2[2];		/* R3 */ \
	xmm5 = sd2[3];		/* I3 */ \
	r3_x3c_djbfft_partial_mem(xmm0,xmm1,xmm2,xmm3,xmm4,xmm5,xmm6,xmm7,sd0[2],sd1[2],sd0[3],sd1[3],screg1,srcreg+srcinc,d1); \
	xmm6 = sd2[0];		/* R3 */ \
	xmm7 = sd2[1];		/* R6 */ \
	sd1[2] = xmm2;		/* Save R1 */ \
	sd1[3] = xmm5;		/* Save I1 */ \
	sd2[0] = xmm0;		/* Save R2 */ \
	sd2[1] = xmm1;		/* Save I2 */ \
	sd2[2] = xmm4;		/* Save R3 */ \
	sd2[3] = xmm3;		/* Save I3 */ \
	r3_x6r_fft_partial_mem(xmm0,xmm1,xmm6,xmm3,xmm4,xmm7,xmm2,xmm5,sd0[0],sd1[0],sd0[1],sd1[1],screg1,screg2,srcreg+srcinc+2*d1,0); \
	sd0[0] = xmm2;			/* Save R1 #1 */ \
	sd0[1] = xmm7;		/* Save R1 #2 */ \
	sd0[2] = xmm1;		/* Save R2 */ \
	sd0[3] = xmm6;		/* Save I2 */ \
	sd1[0] = xmm0;		/* Save R3 */ \
	sd1[1] = xmm5;		/* Save I3 */ \
	srcreg += srcinc; \
}

/* R1 #1 = R1 + R3 + R5 */ \
/* R1 #2 = R2 + R4 + R6 */ \
/* R2 = R1 - R4 + 0.5 * (R2 - R3 - R5 + R6) */ \
/* I2 = 0.866 * (R2 + R3 - R5 - R6) */ \
/* R3 = R1 + R4 - 0.5 * (R2 + R3 + R5 + R6) */ \
/* I3 = 0.866 * (R2 - R3 + R5 - R6) */

#if 0 // UNUSED
#define r3_x6r_fft_mem(mem1,mem2,mem3,mem4,mem5,mem6,screg1,screg2) \
	untested_macro; \
\
	xmm2 = mem5;		/* R5 */ \
	xmm6 = mem3;		/* R3 */ \
	xmm2 += xmm6;		/* R3 + R5 */ \
	xmm1 = mem6;		/* R6 */ \
	xmm7 = mem2;		/* R2 */ \
	xmm1 += xmm7;		/* R2 + R6 */ \
	xmm6 -= mem5;		/* R3 - R5 */ \
	xmm7 -= mem6;		/* R2 - R6 */ \
\
	xmm4 = 0.5; \
	xmm4 *= xmm2;		/* 0.5 * (R3 + R5) */ \
	xmm0 = mem1;		/* R1 */ \
	xmm2 += xmm0;		/* final R1 #1 = R1 + R3 + R5 */ \
	xmm5 = 0.5; \
	xmm5 *= xmm1;		/* 0.5 * (R2 + R6) */ \
	xmm3 = mem4;		/* R4 */ \
	xmm1 += xmm3;		/* final R1 #2 = R2 + R4 + R6 */ \
	xmm6 *= CONST2_P866;		/* new I3 = 0.866 * (R3 - R5) */ \
	xmm0 -= xmm4;		/* new R2 = R1 - 0.5 * (R3 + R5) */ \
	xmm7 *= CONST2_P866;		/* new I2 = 0.866 * (R2 - R6) */ \
	xmm5 -= xmm3;		/* new R3 = 0.5 * (R2 + R6) - R4 */ \
\
	xmm4 = xmm0;		/* Copy R2 */ \
	xmm0 -= xmm5;		/* R2 = R2 - R3 (final R3) */ \
	xmm5 += xmm4;		/* R3 = R2 + R3 (final R2) */ \
	xmm3 = xmm7;		/* Copy I2 */ \
	xmm7 -= xmm6;		/* I2 = I2 - I3 (final I3) */ \
	xmm6 += xmm3;		/* I3 = I2 + I3 (final I2) */ \
\
	xmm3 = xptr(screg1+16);	/* cosine/sine for w^2 */ \
	xmm4 = xmm0;		/* Copy R3 */ \
	xmm0 *= xmm3;		/* A3 = R3 * cosine/sine */ \
	xmm0 -= xmm7;		/* A3 = A3 - I3 */ \
	xmm7 *= xmm3;		/* B3 = I3 * cosine/sine */ \
	xmm7 += xmm4;		/* B3 = B3 + R3 */ \
\
	xmm3 = xptr(screg2+16);	/* cosine/sine for w^1 */ \
	xmm4 = xmm5;		/* Copy R2 */ \
	xmm5 *= xmm3;		/* A2 = R2 * cosine/sine */ \
	xmm5 -= xmm6;		/* A2 = A2 - I2 */ \
	xmm6 *= xmm3;		/* B2 = I2 * cosine/sine */ \
	xmm6 += xmm4;		/* B2 = B2 + R2 */ \
\
	xmm3 = xptr(screg1); \
	xmm0 *= xmm3;		/* A3 = A3 * sine (new R3) */ \
	xmm7 *= xmm3;		/* B3 = B3 * sine (new I3) */ \
	xmm4 = xptr(screg2); \
	xmm5 *= xmm4;		/* A2 = A2 * sine (new R2) */ \
	xmm6 *= xmm4;		/* B2 = B2 * sine (new I2) */ \

#endif

/* This is used in the first levels of pass 2 if pass 1 does the swizzling */ \
/* The six-reals macro and the three-complex share an XMM register. */ \
/* This isn't very efficient, but this macro isn't called a whole lot. */ \

#define r3_fh3cl_six_reals_three_complex_djbfft(srcreg,srcinc,d1,screg1,scoff1,screg2,scoff2) { \
	uintptr_t src_rbx = srcreg+rbx; \
	vec2f64 *sd0 = (vec2f64*)(srcreg); \
	vec2f64 *sd1 = (vec2f64*)(srcreg+(d1)); \
	vec2f64 *sd2 = (vec2f64*)(srcreg+2*(d1)); \
	vec2f64 *bd0 = (vec2f64*)(src_rbx); \
	vec2f64 *bd1 = (vec2f64*)(src_rbx+(d1)); \
	vec2f64 *bd2 = (vec2f64*)(src_rbx+2*(d1)); \
	vec2f64 TMP1, TMP2; \
	r3_h6r_h3c_djbfft_mem(bd1[0],bd0[1],bd2[1],bd1[2],bd0[3],bd2[3], screg1+scoff1,screg2+scoff2,srcreg+srcinc,d1); \
	TMP1 = bd1[1];		/* R3 */ \
	TMP2 = bd1[3];		/* R6 */ \
	sd0[1] = xmm2;		/* Save R1 #1/R1 */ \
	sd0[3] = xmm5;		/* Save R1 #2/I1 */ \
	sd1[1] = xmm0;		/* Save R2 */ \
	sd1[3] = xmm1;		/* Save I2 */ \
	sd2[1] = xmm4;		/* Save R3 */ \
	sd2[3] = xmm3;		/* Save I3 */ \
	r3_h6r_h3c_djbfft_mem(bd0[0],bd2[0],TMP1,bd0[2],bd2[2],TMP2, screg1,screg2,srcreg+srcinc+2*d1,0); \
	sd0[0] = xmm2;			/* Save R1 #1/R1 */ \
	sd0[2] = xmm5;		/* Save R1 #2/I1 */ \
	sd1[0] = xmm0;		/* Save R2 */ \
	sd1[2] = xmm1;		/* Save I2 */ \
	sd2[0] = xmm4;		/* Save R3 */ \
	sd2[2] = xmm3;		/* Save I3 */ \
	srcreg += srcinc; \
}

/* This is used in the later radix-3 levels of pass 2 */ \
/* The six-reals macro and the three-complex share an XMM register. */ \
/* This isn't very efficient, but this macro isn't called a whole lot. */ \

#define r3_h3cl_six_reals_three_complex_djbfft(srcreg,srcinc,d1,screg1,scoff1,screg2,scoff2) { \
	vec2f64 *sd0 = (vec2f64*)(srcreg); \
	vec2f64 *sd1 = (vec2f64*)(srcreg+(d1)); \
	vec2f64 *sd2 = (vec2f64*)(srcreg+2*(d1)); \
	vec2f64 TMP1, TMP2; \
	r3_h6r_h3c_djbfft_mem(sd1[0],sd0[1],sd2[1],sd1[2],sd0[3],sd2[3], screg1+scoff1,screg2+scoff2,srcreg+srcinc,d1); \
	TMP1 = sd1[1];		/* R3 */ \
	TMP2 = sd1[3];		/* R6 */ \
	sd0[1] = xmm2;		/* Save R1 #1/R1 */ \
	sd0[3] = xmm5;		/* Save R1 #2/I1 */ \
	sd1[1] = xmm0;		/* Save R2 */ \
	sd1[3] = xmm1;		/* Save I2 */ \
	sd2[1] = xmm4;		/* Save R3 */ \
	sd2[3] = xmm3;		/* Save I3 */ \
	r3_h6r_h3c_djbfft_mem(sd0[0],sd2[0],TMP1,sd0[2],sd2[2],TMP2, screg1,screg2,srcreg+srcinc+2*d1,0); \
	sd0[0] = xmm2;			/* Save R1 #1/R1 */ \
	sd0[2] = xmm5;		/* Save R1 #2/I1 */ \
	sd1[0] = xmm0;		/* Save R2 */ \
	sd1[2] = xmm1;		/* Save I2 */ \
	sd2[0] = xmm4;		/* Save R3 */ \
	sd2[2] = xmm3;		/* Save I3 */ \
	srcreg += srcinc; \
}


#define r3_x6r_fft_partial_mem(r1,r2,r3,r4,r5,r6,t1,t2,mem1,mem2,mem4,mem5,screg1,screg2,pre1,pre2) \
	t1 = mem5;			/* R5 */ \
	t1 += r3;			/* R3 + R5 */ \
	r3 -= mem5;			/* R3 - R5 */ \
	t2 = mem2;			/* R2 */ \
	t2 -= r6;			/* R2 - R6 */ \
	r6 += mem2;			/* R2 + R6 */ \
\
	r5 = t1;			/* 0.5 * (R3 + R5) */ \
	r1 = mem1;			/* R1 */ \
	t1 += r1;			/* final R1 #1 = R1 + R3 + R5 */ \
	r2 = r6;			/* 0.5 * (R2 + R6) */ \
	r4 = mem4;			/* R4 */ \
	r6 += r4;			/* final R1 #2 = R2 + R4 + R6 */ \
	r3 *= CONST2_P866;	/* new I3 = 0.866 * (R3 - R5) */ \
	r1 -= r5 * 0.5;		/* new R2 = R1 - 0.5 * (R3 + R5) */ \
	t2 *= CONST2_P866;	/* new I2 = 0.866 * (R2 - R6) */ \
	r2 = r2 * 0.5 - r4;	/* new R3 = 0.5 * (R2 + R6) - R4 */ \
\
	xprefetchw(u8ptr(pre1)); \
\
	r5 = r1;			/* Copy R2 */ \
	r1 -= r2;			/* R2 = R2 - R3 (final R3) */ \
	r2 += r5;			/* R3 = R2 + R3 (final R2) */ \
	r4 = t2;			/* Copy I2 */ \
	t2 -= r3;			/* I2 = I2 - I3 (final I3) */ \
	r3 += r4;			/* I3 = I2 + I3 (final I2) */ \
\
	/*if 0 // pre2 NE 0; \
	xprefetchw(u8ptr(pre1+pre2)); \
	endif */ \
\
	r4 = xptr(screg1+16);		/* cosine/sine for w^2 */ \
	r5 = r1;			/* Copy R3 */ \
	r1 = r1 * r4 - t2;	/* A3 = R3 * cosine/sine - I3 */ \
	t2 = t2 * r4 + r5;	/* B3 = I3 * cosine/sine + R3 */ \
\
	r4 = xptr(screg2+16);		/* cosine/sine for w^1 */ \
	r5 = r2;			/* Copy R2 */ \
	r2 = r2 * r4 - r3;	/* A2 = R2 * cosine/sine - I2 */ \
	r3 = r3 * r4 + r5;	/* B2 = I2 * cosine/sine + R2 */ \
\
	r4 = xptr(screg1); \
	r1 *= r4;			/* A3 = A3 * sine (new R3) */ \
	t2 *= r4;			/* B3 = B3 * sine (new I3) */ \
	r5 = xptr(screg2); \
	r2 *= r5;			/* A2 = A2 * sine (new R2) */ \
	r3 *= r5;			/* B2 = B2 * sine (new I2) */ \
\
\

#define r3_h6r_h3c_djbfft_mem(mem1,mem2,mem3,mem4,mem5,mem6,screg1,screg2,pre1,pre2) { \
	double f0,f1,f2,f3,f4,f5,f6,f7; \
	/* Do the three complex part */ \
	f1 = mem2[1];		/* R2 */ \
	f2 = mem3[1];		/* R3 */ \
	f1 -= f2;			/* R2 - R3 */ \
	f4 = mem5[1];		/* I2 */ \
	f5 = mem6[1];		/* I3 */ \
	f4 -= f5;			/* I2 - I3 */ \
	f2 += mem2[1];		/* R2 + R3 */ \
	f1 *= CONST_P866;	/* 0.866 * (R2 - R3) */ \
	f5 += mem5[1];		/* I2 + I3 */ \
	f4 *= CONST_P866;	/* 0.866 * (I2 - I3) */ \
	f6 = f2;			/* 0.5 * (R2 + R3) */ \
	f0 = mem1[1]; \
	f2 += f0;			/* R1 + R2 + R3 (final R1) */ \
	f0 -= f6 * 0.5;		/* (R1-.5R2-.5R3) */ \
	f6 = f5;			/* 0.5 * (I2 + I3) */ \
	f3 = mem4[1];		/* I1 */ \
	f5 += f3;			/* I1 + I2 + I3 (final I1) */ \
	f3 -= f6 * 0.5;		/* (I1-.5I2-.5I3) */ \
	f6 = f0; \
	f0 -= f4;			/* Final R2 */ \
	f7 = f3; \
	f3 -= f1;			/* Final I3 */ \
	f4 += f6;			/* Final R3 */ \
	f1 += f7;			/* Final I2 */ \
\
	xprefetchw(u8ptr(pre1)); \
\
	f6 = f0;			/* Copy R2 */ \
	f7 = f64ptr(screg1+16); \
	f0 = f0 * f7 - f1;	/* A2 = R2 * cosine/sine - I2 */ \
	f1 = f1 * f7 + f6;	/* B2 = I2 * cosine/sine+ R2 */ \
\
	f6 = f4;			/* Copy R3 */ \
	f4 = f4 * f7 + f3;			/* A3 = R3 * cosine/sine + I3 */ \
	f3 = f3 * f7 - f6;			/* B3 = I3 * cosine/sine - R3 */ \
\
	/*if 0 // pre2 NE 0; \
	xprefetchw(u8ptr(pre1+pre2)); \
	endif*/ \
\
	f7 = f64ptr(screg1); \
	f0 *= f7;			/* A2 = A2 * sine (new R2) */ \
	f1 *= f7;			/* B2 = B2 * sine (new I2) */ \
	f4 *= f7;			/* A3 = A3 * sine (new R3) */ \
	f3 *= f7;			/* B3 = B3 * sine (new I3) */ \
\
;	/* Do the six-reals part */ \
	xmm2[1] = f2;		/* Copy to high part of XMM register */ \
	xmm5[1] = f5; \
	xmm0[1] = f0; \
	xmm1[1] = f1; \
	xmm4[1] = f4; \
	xmm3[1] = f3; \
\
	f2 = mem5[0];		/* R5 */ \
	f1 = mem3[0];		/* R3 */ \
	f2 += f1;			/* R3 + R5 */ \
	f1 -= mem5[0];		/* R3 - R5 */ \
	f3 = mem2[0];		/* R2 */ \
	f5 = mem6[0];		/* R6 */ \
	f3 -= f5;			/* R2 - R6 */ \
	f5 += mem2[0];		/* R2 + R6 */ \
\
	f7 = f2;			/* 0.5 * (R3 + R5) */ \
	f4 = mem1[0];		/* R1 */ \
	f2 += f4;			/* final R1 #1 = R1 + R3 + R5 */ \
	f0 = f5;			/* 0.5 * (R2 + R6) */ \
	f6 = mem4[0];		/* R4 */ \
	f5 += f6;			/* final R1 #2 = R2 + R4 + R6 */ \
	f1 *= CONST_P866;	/* new I3 = 0.866 * (R3 - R5) */ \
	f4 -= f7 * 0.5;		/* new R2 = R1 - 0.5 * (R3 + R5) */ \
	f3 *= CONST_P866;	/* new I2 = 0.866 * (R2 - R6) */ \
	f0 = f0 * 0.5 - f6;	/* new R3 = 0.5 * (R2 + R6) - R4 */ \
\
	f7 = f4;			/* Copy R2 */ \
	f4 -= f0;			/* R2 = R2 - R3 (final R3) */ \
	f0 += f7;			/* R3 = R2 + R3 (final R2) */ \
	f6 = f3;			/* Copy I2 */ \
	f3 -= f1;			/* I2 = I2 - I3 (final I3) */ \
	f1 += f6;			/* I3 = I2 + I3 (final I2) */ \
\
	f6 = f64ptr(screg1+16);	/* cosine/sine for w^2 */ \
	f7 = f4;			/* Copy R3 */ \
	f4 = f4 * f6 - f3;			/* A3 = R3 * cosine/sine - I3 */ \
	f3 = f3 * f6 + f7;			/* B3 = I3 * cosine/sine + R3 */ \
\
	f6 = f64ptr(screg2+8);	/* cosine/sine for w^1 */ \
	f7 = f0;			/* Copy R2 */ \
	f0 = f0 * f6 - f1;			/* A2 = R2 * cosine/sine - I2 */ \
	f1 = f1 * f6 + f7;			/* B2 = I2 * cosine/sine + R2 */ \
\
	f6 = f64ptr(screg1); \
	f4 *= f6;			/* A3 = A3 * sine (new R3) */ \
	f3 *= f6;			/* B3 = B3 * sine (new I3) */ \
	f7 = f64ptr(screg2); \
	f0 *= f7;			/* A2 = A2 * sine (new R2) */ \
	f1 *= f7;			/* B2 = B2 * sine (new I2) */ \
	xmm0[0] = f0; \
	xmm1[0] = f1; \
	xmm2[0] = f2; \
	xmm3[0] = f3; \
	xmm4[0] = f4; \
	xmm5[0] = f5; \
}

/* */ \
/* ************************************* six-reals-unfft variants ****************************************** */ \
/* */ \
\

#define r3_x3cl_six_reals_unfft(srcreg,srcinc,d1,screg1,scoff1,screg2,scoff2) { \
	vec2f64 *sd0 = (vec2f64*)(srcreg); \
	vec2f64 *sd1 = (vec2f64*)(srcreg+(d1)); \
	vec2f64 *sd2 = (vec2f64*)(srcreg+2*(d1)); \
	xmm2 = sd1[0];	/* R2 */ \
	xmm3 = sd1[2];	/* I2 */ \
	r3_x6r_unfft_partial_mem(xmm0,xmm1,xmm2,xmm3,xmm4,xmm5,xmm6,xmm7, sd0[0],sd0[2],sd2[0],sd2[2], screg1,0,screg2,0,srcreg+srcinc,d1); \
	sd0[0] = xmm4;		/* Save R1 */ \
	sd0[2] = xmm1;	/* Save R4 */ \
	sd2[0] = xmm5;	/* Save R2 */ \
	sd2[2] = xmm0;	/* Save R5 */ \
	xmm0 = sd1[1];	/* R2 */ \
	xmm1 = sd1[3];	/* I2 */ \
	sd1[1] = xmm3;	/* Save R3 */ \
	sd1[3] = xmm7;	/* Save R6 */ \
	r3_x6r_unfft_partial_mem(xmm2,xmm3,xmm0,xmm1,xmm4,xmm5,xmm6,xmm7, sd0[1],sd0[3],sd2[1],sd2[3], screg1,scoff1,screg2,scoff2,srcreg+srcinc+2*d1,0); \
	sd0[1] = xmm5;	/* Save R2 */ \
	sd0[3] = xmm2;	/* Save R5 */ \
	sd1[0] = xmm4;	/* Save R1 */ \
	sd1[2] = xmm3;	/* Save R4 */ \
	sd2[1] = xmm1;	/* Save R3 */ \
	sd2[3] = xmm7;	/* Save R6 */ \
	srcreg += srcinc; \
}

/* R1 = R1#1 + (R2 + R3) */ \
/* R2 = R1#2 + 0.5 * (R2 - R3) + 0.866 * (I2 + I3) */ \
/* R3 = R1#1 - 0.5 * (R2 + R3) + 0.866 * (I2 - I3) */ \
/* R4 = R1#2 - (R2 - R3) */ \
/* R5 = R1#1 - 0.5 * (R2 + R3) - 0.866 * (I2 - I3) */ \
/* R6 = R1#2 + 0.5 * (R2 - R3) - 0.866 * (I2 + I3) */ \
\

#define r3_x6r_unfft_partial_mem(r1,r2,r3,r4,r5,r6,t1,t2,memr1_1,memr1_2,memr3,memi3,screg1,scoff1,screg2,scoff2,pre1,pre2) { \
	vec2f64 s1, s2; \
	s2 = xptr(screg2+scoff2+16); \
	s1 = xptr(screg1+scoff1+16); \
	t1 = r3;			/* Copy R2 */ \
	r5 = memr3;			/* R3 */ \
	r3 = r3 * s2 + r4;	/* A2 = R2 * cosine/sine + I2 */ \
	r6 = memi3;			/* I3 */ \
	r5 = r5 * s1 + r6;	/* A3 = R3 * cosine/sine + I3 */ \
	r4 = r4 * s2 - t1;	/* B2 = I2 * cosine/sine - R2 */ \
	r6 = r6 * s1 - memr3;	/* B3 = I3 * cosine/sine - R3 */ \
	s2 = xptr(screg2+scoff2); \
	s1 = xptr(screg1+scoff1); \
	r3 *= s2;	/* A2 = A2 * sine (new R2) */ \
	r5 *= s1;	/* A3 = A3 * sine (new R3) */ \
	r4 *= s2;	/* B2 = B2 * sine (new I2) */ \
	r6 *= s1;	/* B3 = B3 * sine (new I3) */ \
\
	xprefetchw(u8ptr(pre1)); \
\
	t1 = r3;			/* Copy R2 */ \
	r3 -= r5;			/* R2 - R3 */ \
	r5 += t1;			/* R2 + R3 */ \
	t2 = r4;			/* Copy I2 */ \
	r4 -= r6;			/* I2 - I3 */ \
	r6 += t2;			/* I2 + I3 */ \
\
	/*if 0 // pre2 NE 0; \
	xprefetchw(u8ptr(pre1+pre2)); \
	#endif*/ \
\
	t2 = r3;			/* 0.5 * (R2 - R3) */ \
	t1 = r5;			/* 0.5 * (R2 + R3) */ \
	r4 *= CONST2_P866;		/* 0.866 * (I2 - I3) */ \
	r6 *= CONST2_P866;		/* 0.866 * (I2 + I3) */ \
	r2 = memr1_2;		/* R1#2 */ \
	t2 = t2 * 0.5 + r2;		/* R1#2 + 0.5 * (R2 - R3) */ \
	r1 = memr1_1;		/* R1#1 */ \
	r5 += r1;			/* final R1 = R1#1 + (R2 + R3) */ \
	r1 -= t1 * 0.5;			/* R1#1 - 0.5 * (R2 + R3) */ \
	r2 -= r3;			/* final R4 = R1#2 - (R2 - R3) */ \
	t1 = r6;			/* Copy 0.866 * (I2 + I3) */ \
	r6 += t2;			/* final R2 = R1#2 + 0.5 * (R2 - R3) + 0.866 * (I2 + I3) */ \
	t2 -= t1;			/* final R6 = R1#2 + 0.5 * (R2 - R3) - 0.866 * (I2 + I3) */ \
	r3 = r4;			/* Copy 0.866 * (I2 - I3) */ \
	r4 += r1;			/* final R3 = R1#1 - 0.5 * (R2 + R3) + 0.866 * (I2 - I3) */ \
	r1 -= r3;			/* final R5 = R1#1 - 0.5 * (R2 + R3) - 0.866 * (I2 - I3) */ \
}

/* Macro to do an six_reals_unfft and a three_complex_djbunfft in pass 2. */ \
/* The six-reals operation is done in the lower half of the XMM */ \
/* register.  This isn't very efficient, but this macro isn't called a whole lot. */ \
\

#define r3_h3cl_six_reals_three_complex_djbunfft(srcreg,srcinc,d1,screg1,scoff1,screg2,scoff2) { \
	vec2f64 *sd0 = (vec2f64*)(srcreg); \
	vec2f64 *sd1 = (vec2f64*)(srcreg+(d1)); \
	vec2f64 *sd2 = (vec2f64*)(srcreg+2*(d1)); \
	vec2f64 TMP1, TMP2; \
	r3_h6r_h3c_djbunfft_mem(sd0[1],sd0[3],sd1[1],sd1[3],sd2[1],sd2[3], screg1+scoff1,screg2+scoff2,srcreg+srcinc,d1); \
	xmm6 = sd1[0];	/* Load R2 */ \
	xmm7 = sd1[2];	/* Load I2 */ \
	TMP1 = xmm6; \
	TMP2 = xmm7; \
	sd1[0] = xmm4;	/* Save R1 */ \
	sd0[1] = xmm3;	/* Save R2 */ \
	sd2[1] = xmm0;	/* Save R3 */ \
	sd1[2] = xmm5;	/* Save R4 */ \
	sd0[3] = xmm1;	/* Save R5 */ \
	sd2[3] = xmm2;	/* Save R6 */ \
	r3_h6r_h3c_djbunfft_mem(sd0[0],sd0[2],TMP1,TMP2,sd2[0],sd2[2], screg1,screg2,srcreg+srcinc+2*d1,0); \
	sd0[0] = xmm4;		/* Save R1 */ \
	sd2[0] = xmm3;	/* Save R2 */ \
	sd1[1] = xmm0;	/* Save R3 */ \
	sd0[2] = xmm5;	/* Save R4 */ \
	sd2[2] = xmm1;	/* Save R5 */ \
	sd1[3] = xmm2;	/* Save R6 */ \
	srcreg += srcinc; \
}



#define r3_h6r_h3c_djbunfft_mem(memr1_1,memr1_2,memr2,memi2,memr3,memi3,screg1,screg2,pre1,pre2) { \
	double f0,f1,f2,f3,f4,f5,f6,f7; \
	double s10 = f64ptr(screg1); \
	double s11 = f64ptr(screg1+16); \
\
	/* Do the three complex part */ \
	f2 = memr2[1];		/* R2 */ \
	f6 = f2;			/* Copy R2 */ \
	f7 = s11;	/* cosine/sine */ \
	f4 = memr3[1];		/* R3 */ \
	f3 = memi2[1];		/* I2 */ \
	f5 = memi3[1];		/* I3 */ \
	f2 = f2 * f7 + f3;	/* A2 = R2 * cosine/sine + I2 */ \
	f4 = f4 * f7 - f5;	/* A3 = R3 * cosine/sine - I3 */ \
	f3 = f3 * f7 - f6;	/* B2 = I2 * cosine/sine - R2 */ \
	f5 = f5 * f7 + memr3[1];	/* B3 = I3 * cosine/sine + R3 */ \
	f7 = s10;	/* sine */ \
	f2 *= f7;			/* A2 = A2 * sine (final R2) */ \
	f4 *= f7;			/* A3 = A3 * sine (final R3) */ \
	f3 *= f7;			/* B2 = B2 * sine (final I2) */ \
	f5 *= f7;			/* B3 = B3 * sine (final I3) */ \
\
	xprefetchw(u8ptr(pre1)); \
\
	f6 = f2;			/* Copy R2 */ \
	f2 -= f4;			/* R2 - R3 */ \
	f4 += f6;			/* R2 + R3 */ \
	f2 *= CONST_P866;	/* 0.866 * (R2 - R3) */ \
	f7 = f3;			/* Copy I2 */ \
	f3 -= f5;			/* I2 - I3 */ \
	f5 += f7;			/* I2 + I3 */ \
	f3 *= CONST_P866;	/* 0.866 * (I2 - I3) */ \
	f6 = f4;			/* 0.5 * (R2 + R3) */ \
	f0 = memr1_1[1];	/* R1 */ \
	f4 += f0;			/* R1 + R2 + R3 (final R1) */ \
	f0 -= f6 * 0.5;		/* (R1-.5R2-.5R3) */ \
	f6 = f5;			/* 0.5 * (I2 + I3) */ \
	f1 = memr1_2[1];	/* I1 */ \
	f5 += f1;			/* I1 + I2 + I3 (final I1) */ \
	f1 -= f6 * 0.5;		/* (I1-.5I2-.5I3) */ \
	f6 = f0; \
	f0 -= f3;			/* Final R3 */ \
	f3 += f6;			/* Final R2 */ \
	f7 = f1; \
	f1 -= f2;			/* Final I2 */ \
	f2 += f7;			/* Final I3 */ \
\
	/*if 0 // pre2 NE 0; \
	xprefetchw(u8ptr(pre1+pre2)); \
	#endif*/ \
\
	/* Do the six reals part */ \
\
	xmm4[1] = f4; \
	xmm3[1] = f3; \
	xmm0[1] = f0; \
	xmm5[1] = f5; \
	xmm1[1] = f1; \
	xmm2[1] = f2; \
\
	vec2f64 s2 = xptr(screg2); \
	f6 = memr2[0];		/* R2 */ \
	f7 = f6;			/* Copy R2 */ \
	f4 = memr3[0];		/* R3 */ \
	f0 = memi2[0];		/* R2 */ \
	f3 = memi3[0];		/* I3 */ \
	f6 = f6 * s2[1] + f0;	/* A2 = R2 * cosine/sine + I2 */ \
	f4 = f4 * s11 + f3;		/* A3 = R3 * cosine/sine + I3 */ \
	f0 = f0 * s2[1] - f7;	/* B2 = I2 * cosine/sine - R2 */ \
	f3 = f3 * s11 - memr3[0];	/* B3 = I3 * cosine/sine - R3 */ \
	f6 *= s2[0];		/* A2 = A2 * sine (new R2) */ \
	f4 *= s10;			/* A3 = A3 * sine (new R3) */ \
	f0 *= s2[0];		/* B2 = B2 * sine (new I2) */ \
	f3 *= s10;			/* B3 = B3 * sine (new I3) */ \
\
	f7 = f6;			/* Copy R2 */ \
	f6 -= f4;			/* R2 - R3 */ \
	f4 += f7;			/* R2 + R3 */ \
	f2 = f0;			/* Copy I2 */ \
	f0 -= f3;			/* I2 - I3 */ \
	f3 += f2;			/* I2 + I3 */ \
	f2 = f6;			/* 0.5 * (R2 - R3) */ \
	f7 = f4;			/* 0.5 * (R2 + R3) */ \
	f0 *= CONST_P866;	/* 0.866 * (I2 - I3) */ \
	f3 *= CONST_P866;	/* 0.866 * (I2 + I3) */ \
	f5 = memr1_2[0];	/* R1#2 */ \
	f2 = f2 * 0.5 + f5;	/* R1#2 + 0.5 * (R2 - R3) */ \
	f1 = memr1_1[0];	/* R1#1 */ \
	f4 += f1;			/* final R1 = R1#1 + (R2 + R3) */ \
	f1 -= f7 * 0.5;		/* R1#1 - 0.5 * (R2 + R3) */ \
	f5 -= f6;			/* final R4 = R1#2 - (R2 - R3) */ \
	f7 = f3;			/* Copy 0.866 * (I2 + I3) */ \
	f3 += f2;			/* final R2 = R1#2 + 0.5 * (R2 - R3) + 0.866 * (I2 + I3) */ \
	f2 -= f7;			/* final R6 = R1#2 + 0.5 * (R2 - R3) - 0.866 * (I2 + I3) */ \
	f6 = f0;			/* Copy 0.866 * (I2 - I3) */ \
	f0 += f1;			/* final R3 = R1#1 - 0.5 * (R2 + R3) + 0.866 * (I2 - I3) */ \
	f1 -= f6;			/* final R5 = R1#1 - 0.5 * (R2 + R3) - 0.866 * (I2 - I3) */ \
	xmm0[0] = f0; \
	xmm1[0] = f1; \
	xmm2[0] = f2; \
	xmm3[0] = f3; \
	xmm4[0] = f4; \
	xmm5[0] = f5; \
}
