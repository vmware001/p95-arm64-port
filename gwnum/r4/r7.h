#pragma once

/* Copyright 2009-2010 - Mersenne Research, Inc.  All rights reserved */ \
/* Author:  George Woltman */ \
/* Email: woltman@alum.mit.edu */ \
/* */ \
\
/* */ \
/* Macros for a radix-7 first step in a real FFT. */ \
/* */ \
\
/* */ \
/* ************************************* 28-reals-first-fft variants ****************************************** */ \
/* */ \
\
/* These macros operate on 28 reals doing 4.807 levels of the FFT, applying */ \
/* the sin/cos multipliers afterwards.  The output is 2 reals and 13 complex numbers. */ \
\
/*r7_x7cl_28_reals_first_fft_preload MACRO */ \
/*	r7_x7cl_28_reals_first_fft_cmn_preload */ \
/*	ENDM */ \
/* */ \
/*r7_x7cl_28_reals_first_fft MACRO srcreg,srcinc,d1,screg */ \
/*	r7_x7cl_28_reals_first_fft_cmn srcreg,rbx,srcinc,d1,screg */ \
/*	ENDM */ \
\

#define r7_x7cl_28_reals_first_fft_scratch_preload \
	r7_x7cl_28_reals_first_fft_cmn_preload; \
\
\

#define r7_x7cl_28_reals_first_fft_scratch(srcreg,srcinc,d1,screg) \
	r7_x7cl_28_reals_first_fft_cmn(srcreg,srcinc,d1,screg); \
\
\
/* To calculate a 28-reals FFT, we calculate 28 complex values in a brute force way (using a shorthand notation): */ \
/* r1 + r2 + ... + r28	*  w^0000000000... */ \
/* r1 + r2 + ... + r28	*  w^0123456789A... */ \
/* r1 + r2 + ... + r28	*  w^02468ACE.... */ \
/*    ... */ \
/* r1 + r2 + ... + r28	*  w^...A987654321 */ \
/* Note that Hermetian symmetry means we won't need to calculate the last 14 complex values. */ \
/* */ \
/* The sin/cos values (w = 28th root of unity) are: */ \
/* w^1 = .975 + .223i */ \
/* w^2 = .901 + .434i */ \
/* w^3 = .782 + .623i */ \
/* w^4 = .623 + .782i */ \
/* w^5 = .434 + .901i */ \
/* w^6 = .223 + .975i */ \
/* w^7 = 0 + 1i */ \
/* w^8 = -.223 + .975i */ \
/* w^9 = -.434 + .901i */ \
/* w^10 = -.623 + .782i */ \
/* w^11 = -.782 + .623i */ \
/* w^12 = -.901 + .434i */ \
/* w^13 = -.975 + .223i */ \
/* w^14 = -1 */ \
/* */ \
/* Applying the sin/cos values above (and noting that combining r2 and r28, r3 and r27, etc. will simplify calculations): */ \
/* reals: */ \
/* r1     +(r2+r28)     +(r3+r27)     +(r4+r26)     +(r5+r25)     +(r6+r24)     +(r7+r23) + (r8+r22)     +(r9+r21)     +(r10+r20)     +(r11+r19)     +(r12+r18)     +(r13+r17)     +(r14+r16) + r15 */ \
/* r1 +.975(r2+r28) +.901(r3+r27) +.782(r4+r26) +.623(r5+r25) +.434(r6+r24) +.223(r7+r23)            -.223(r9+r21) -.434(r10+r20) -.623(r11+r19) -.782(r12+r18) -.901(r13+r17) -.975(r14+r16) - r15 */ \
/* r1 +.901(r2+r28) +.623(r3+r27) +.223(r4+r26) -.223(r5+r25) -.623(r6+r24) -.901(r7+r23) - (r8+r22) -.901(r9+r21) -.623(r10+r20) -.223(r11+r19) +.223(r12+r18) +.623(r13+r17) +.901(r14+r16) + r15 */ \
/* r1 +.782(r2+r28) +.223(r3+r27) -.434(r4+r26) -.901(r5+r25) -.975(r6+r24) -.623(r7+r23)            +.623(r9+r21) +.975(r10+r20) +.901(r11+r19) +.434(r12+r18) -.223(r13+r17) -.782(r14+r16) - r15 */ \
/* r1 +.623(r2+r28) -.223(r3+r27) -.901(r4+r26) -.901(r5+r25) -.223(r6+r24) +.623(r7+r23) + (r8+r22) +.623(r9+r21) -.223(r10+r20) -.901(r11+r19) -.901(r12+r18) -.223(r13+r17) +.623(r14+r16) + r15 */ \
/* r1 +.434(r2+r28) -.623(r3+r27) -.975(r4+r26) -.223(r5+r25) +.782(r6+r24) +.901(r7+r23)            -.901(r9+r21) -.782(r10+r20) +.223(r11+r19) +.975(r12+r18) +.623(r13+r17) -.434(r14+r16) - r15 */ \
/* r1 +.223(r2+r28) -.901(r3+r27) -.623(r4+r26) +.623(r5+r25) +.901(r6+r24) -.223(r7+r23) - (r8+r22) -.223(r9+r21) +.901(r10+r20) +.623(r11+r19) -.623(r12+r18) -.901(r13+r17) +.223(r14+r16) + r15 */ \
/* r1                   -(r3+r27)                   +(r5+r25)                   -(r7+r23)                +(r9+r21)                    -(r11+r19)                    +(r13+r17)                - r15 */ \
/* r1 -.223(r2+r28) -.901(r3+r27) +.623(r4+r26) +.623(r5+r25) -.901(r6+r24) -.223(r7+r23) + (r8+r22) -.223(r9+r21) -.901(r10+r20) +.623(r11+r19) +.623(r12+r18) -.901(r13+r17) -.223(r14+r16) + r15 */ \
/* r1 -.434(r2+r28) -.623(r3+r27) +.975(r4+r26) -.223(r5+r25) -.782(r6+r24) +.901(r7+r23)            -.901(r9+r21) +.782(r10+r20) +.223(r11+r19) -.975(r12+r18) +.623(r13+r17) +.434(r14+r16) - r15 */ \
/* r1 -.623(r2+r28) -.223(r3+r27) +.901(r4+r26) -.901(r5+r25) +.223(r6+r24) +.623(r7+r23) - (r8+r22) +.623(r9+r21) +.223(r10+r20) -.901(r11+r19) +.901(r12+r18) -.223(r13+r17) -.623(r14+r16) + r15 */ \
/* r1 -.782(r2+r28) +.223(r3+r27) +.434(r4+r26) -.901(r5+r25) +.975(r6+r24) -.623(r7+r23)            +.623(r9+r21) -.975(r10+r20) +.901(r11+r19) -.434(r12+r18) -.223(r13+r17) +.782(r14+r16) - r15 */ \
/* r1 -.901(r2+r28) +.623(r3+r27) -.223(r4+r26) -.223(r5+r25) +.623(r6+r24) -.901(r7+r23) + (r8+r22) -.901(r9+r21) +.623(r10+r20) -.223(r11+r19) -.223(r12+r18) +.623(r13+r17) -.901(r14+r16) + r15 */ \
/* r1 -.975(r2+r28) +.901(r3+r27) -.782(r4+r26) +.623(r5+r25) -.434(r6+r24) +.223(r7+r23)            -.223(r9+r21) +.434(r10+r20) -.623(r11+r19) +.782(r12+r18) -.901(r13+r17) +.975(r14+r16) - r15 */ \
/* r1     -(r2+r28)     +(r3+r27)     -(r4+r26)     +(r5+r25)     -(r6+r24)     +(r7+r23) - (r8+r22)     +(r9+r21)     -(r10+r20)     +(r11+r19)     -(r12+r18)     +(r13+r17)     -(r14+r16) + r15 */ \
/* */ \
/* imaginarys: */ \
/* 0 */ \
/* +.223(r2-r28) +.434(r3-r27) +.623(r4-r26) +.782(r5-r25) +.901(r6-r24) +.975(r7-r23) + (r8-r22) +.975(r9-r21) +.901(r10-r20) +.782(r11-r19) +.623(r12-r18) +.434(r13-r17) +.223(r14-r16) */ \
/* +.434(r2-r28) +.782(r3-r27) +.975(r4-r26) +.975(r5-r25) +.782(r6-r24) +.434(r7-r23)            -.434(r9-r21) -.782(r10-r20) -.975(r11-r19) -.975(r12-r18) -.782(r13-r17) -.434(r14-r16) */ \
/* +.623(r2-r28) +.975(r3-r27) +.901(r4-r26) +.434(r5-r25) -.223(r6-r24) -.782(r7-r23) - (r8-r22) -.782(r9-r21) -.223(r10-r20) +.434(r11-r19) +.901(r12-r18) +.975(r13-r17) +.623(r14-r16) */ \
/* +.782(r2-r28) +.975(r3-r27) +.434(r4-r26) -.434(r5-r25) -.975(r6-r24) -.782(r7-r23)            +.782(r9-r21) +.975(r10-r20) +.434(r11-r19) -.434(r12-r18) -.975(r13-r17) -.782(r14-r16) */ \
/* +.901(r2-r28) +.782(r3-r27) -.223(r4-r26) -.975(r5-r25) -.623(r6-r24) +.434(r7-r23) + (r8-r22) +.434(r9-r21) -.623(r10-r20) -.975(r11-r19) -.223(r12-r18) +.782(r13-r17) +.901(r14-r16) */ \
/* +.975(r2-r28) +.434(r3-r27) -.782(r4-r26) -.782(r5-r25) +.434(r6-r24) +.975(r7-r23)            -.975(r9-r21) -.434(r10-r20) +.782(r11-r19) +.782(r12-r18) -.434(r13-r17) -.975(r14-r16) */ \
/*      (r2-r28)                   -(r4-r26)                   +(r6-r24)               - (r8-r22)                   +(r10-r20)                    -(r12-r18)                    +(r14-r16) */ \
/* +.975(r2-r28) -.434(r3-r27) -.782(r4-r26) +.782(r5-r25) +.434(r6-r24) -.975(r7-r23)            +.975(r9-r21) -.434(r10-r20) -.782(r11-r19) +.782(r12-r18) +.434(r13-r17) -.975(r14-r16) */ \
/* +.901(r2-r28) -.782(r3-r27) -.223(r4-r26) +.975(r5-r25) -.623(r6-r24) -.434(r7-r23) + (r8-r22) -.434(r9-r21) -.623(r10-r20) +.975(r11-r19) -.223(r12-r18) -.782(r13-r17) +.901(r14-r16) */ \
/* +.782(r2-r28) -.975(r3-r27) +.434(r4-r26) +.434(r5-r25) -.975(r6-r24) +.782(r7-r23)            -.782(r9-r21) +.975(r10-r20) -.434(r11-r19) -.434(r12-r18) +.975(r13-r17) -.782(r14-r16) */ \
/* +.623(r2-r28) -.975(r3-r27) +.901(r4-r26) -.434(r5-r25) -.223(r6-r24) +.782(r7-r23) - (r8-r22) +.782(r9-r21) -.223(r10-r20) -.434(r11-r19) +.901(r12-r18) -.975(r13-r17) +.623(r14-r16) */ \
/* +.434(r2-r28) -.782(r3-r27) +.975(r4-r26) -.975(r5-r25) +.782(r6-r24) -.434(r7-r23)            +.434(r9-r21) -.782(r10-r20) +.975(r11-r19) -.975(r12-r18) +.782(r13-r17) -.434(r14-r16) */ \
/* +.223(r2-r28) -.434(r3-r27) +.623(r4-r26) -.782(r5-r25) +.901(r6-r24) -.975(r7-r23) + (r8-r22) -.975(r9-r21) +.901(r10-r20) -.782(r11-r19) +.623(r12-r18) -.434(r13-r17) +.223(r14-r16) */ \
/* 0 */ \
/* */ \
/* There are many more symmetries we can take advantage of.   For example, the (r2+/-r28) column */ \
/* always has the same multiplier as the (r14+/-r16) column.  This is true for all the "even" columns. */ \
/* Also the computations for the 2nd row are very similar to the computations for the 14th row, */ \
/* the 3rd row are similar to the 13th, etc.  Finally, note that for the odd columns, there are */ \
/* only three multipliers to apply and can be combined with every fourth column. */ \
/* */ \
/* Lastly, output would normally be 13 complex and 2 reals.  but the users of this routine */ \
/* expect us to "back up" the 2 reals by one level.  That is: */ \
/*	real #1A:  r1 + r3+r27 + r5+r25 + ... */ \
/*	real #1B:  r2+r28 + r4+r26 + ... */ \
\
/* Store intermediate results in g->u.xmm.XMM_COL_MULTS (an 8KB buffer used in normalization) */ \
\

#define r7_x7cl_28_reals_first_fft_cmn_preload


#define r7_x7cl_28_reals_first_fft_cmn(srcreg,srcinc,d1,screg) { \
	vec2f64 t0,t1,t2,t3,t4,t5,t6,t7,t8,t9; \
	vec2f64 t10,t11,t12,t13,t14,t15,t16,t17,t18,t19; \
	vec2f64 t20,t21,t22,t23,t24,t25; \
\
	vec2f64 *sd0 = (vec2f64*)(srcreg); \
	vec2f64 *sd1 = (vec2f64*)(srcreg+d1*1); \
	vec2f64 *sd2 = (vec2f64*)(srcreg+d1*2); \
	vec2f64 *sd3 = (vec2f64*)(srcreg+d1*3); \
	vec2f64 *sd4 = (vec2f64*)(srcreg+d1*4); \
	vec2f64 *sd5 = (vec2f64*)(srcreg+d1*5); \
	vec2f64 *sd6 = (vec2f64*)(srcreg+d1*6); \
	vec2f64 *sc = (vec2f64*)(screg); \
	\
	/* Do the odd columns for the real results */ \
\
	xmm4 = sd4[0];			/* r5 */ \
	xmm4 += sd3[3];			/* r5+r25 */ \
	xmm3 = sd0[0];			/* r1 */ \
	xmm5 = sd1[1];			/* r9 */ \
	xmm5 += sd6[2];			/* r9+r21 */ \
	xmm6 = sd5[1];			/* r13 */ \
	xmm6 += sd2[2];			/* r13+r17 */ \
	xmm0 = xmm3 + xmm4;					/* r1+(r5+r25) */ \
	xmm1 = xmm3 + xmm4 * CONST2_P623;	/* r1+.623(r5+r25) */ \
	xmm2 = xmm3 - xmm4 * CONST2_P223;	/* r1-.223(r5+r25) */ \
	xmm3 = xmm3 - xmm4 * CONST2_P901;	/* r1-.901(r5+r25) */ \
	xmm0 += xmm5;				/* r1+(r5+r25)+(r9+r21) */ \
	xmm1 -= xmm5 * CONST2_P223;	/* r1+.623(r5+r25)-.223(r9+r21) */ \
	xmm2 -= xmm5 * CONST2_P901;	/* r1-.223(r5+r25)-.901(r9+r21) */ \
	xmm3 += xmm5 * CONST2_P623;	/* r1-.901(r5+r25)+.623(r9+r21) */ \
	xmm0 += xmm6;				/* r1+(r5+r25)+(r9+r21)+(r13+r17) */ \
	xmm1 -= xmm6 * CONST2_P901;	/* r1+.623(r5+r25)-.223(r9+r21)-.901(r13+r17) */ \
	xmm2 += xmm6 * CONST2_P623;	/* r1-.223(r5+r25)-.901(r9+r21)+.623(r13+r17) */ \
	xmm3 -= xmm6 * CONST2_P223;	/* r1-.901(r5+r25)+.623(r9+r21)-.223(r13+r17) */ \
\
	t0 = xmm0; \
	t1 = xmm1; \
	t2 = xmm2; \
	t3 = xmm3; \
\
	xmm4 = sd2[0];			/* r3 */ \
	xmm4 += sd5[3];			/* r3+r27 */ \
	xmm3 = sd0[2];			/* r15 */ \
	xmm5 = sd6[0];			/* r7 */ \
	xmm5 += sd1[3];			/* r7+r23 */ \
	xmm6 = sd3[1];			/* r11 */ \
	xmm6 += sd4[2];			/* r11+r19 */ \
	xmm0 = xmm4 + xmm3;					/* (r3+r27)+r15 */ \
	xmm1 = xmm4 * CONST2_P901 - xmm3;	/* .901(r3+r27)-r15 */ \
	xmm2 = xmm4 * CONST2_P623 + xmm3;	/* .623(r3+r27)+r15 */ \
	xmm3 = xmm4 * CONST2_P223 - xmm3;	/* .223(r3+r27)-r15 */ \
	xmm0 += xmm5;				/* (r3+r27)+(r7+r23)+r15 */ \
	xmm1 += xmm5 * CONST2_P223;	/* .901(r3+r27)+.223(r7+r23)-r15 */ \
	xmm2 -= xmm5 * CONST2_P901;	/* .623(r3+r27)-.901(r7+r23)+r15 */ \
	xmm3 -= xmm5 * CONST2_P623;	/* .223(r3+r27)-.623(r7+r23)-r15 */ \
	xmm0 += xmm6;				/* (r3+r27)+(r7+r23)+(r11+r19)+r15 */ \
	xmm1 -= xmm6 * CONST2_P623;	/* .901(r3+r27)+.223(r7+r23)-.623(r11+r19)-r15 */ \
	xmm2 -= xmm6 * CONST2_P223;	/* .623(r3+r27)-.901(r7+r23)-.223(r11+r19)+r15 */ \
	xmm3 += xmm6 * CONST2_P901;	/* .223(r3+r27)-.623(r7+r23)+.901(r11+r19)-r15 */ \
\
	xmm4 = t0; \
	xmm4 -= xmm0;			/* Real odd-cols row #8 (final real #8) */ \
	xmm0 += t0;				/* Real odd-cols row #1 (final real #1A) */ \
\
	xmm5 = t1; \
	xmm5 -= xmm1;			/* Real odd-cols row #7 */ \
	xmm1 += t1;				/* Real odd-cols row #2 */ \
\
	xmm6 = t2; \
	xmm6 -= xmm2;			/* Real odd-cols row #6 */ \
	xmm2 += t2;				/* Real odd-cols row #3 */ \
\
	xmm7 = t3; \
	xmm7 -= xmm3;			/* Real odd-cols row #5 */ \
	xmm3 += t3;				/* Real odd-cols row #4 */ \
\
	t6 = xmm4;				/* Real #8 */ \
	t5 = xmm5;				/* Real odd-cols row #7 */ \
	t0 = xmm1;				/* Real odd-cols row #2 */ \
	t4 = xmm6;				/* Real odd-cols row #6 */ \
	t1 = xmm2;				/* Real odd-cols row #3 */ \
	t3 = xmm7;				/* Real odd-cols row #5 */ \
	t2 = xmm3;				/* Real odd-cols row #4 */ \
\
;	/* Do the even columns for the real results */ \
\
	xmm5 = sd1[0];			/* r2 */ \
	xmm5 += sd6[3];			/* r2+r28 */ \
	xmm6 = sd6[1];			/* r14 */ \
	xmm6 += sd1[2];			/* r14+r16 */ \
	sd0[0] = xmm0;			/* Final real #1A */ \
	xmm7 = xmm5 - xmm6;		/* (r2+r28)-(r14+r16) */ \
	xmm0 = xmm5 + xmm6;		/* (r2+r28)+(r14+r16) */ \
\
	xmm4 = sd3[0];			/* r4 */ \
	xmm4 += sd4[3];			/* r4+r26 */ \
	xmm1 = sd4[1];			/* r12 */ \
	xmm1 += sd3[2];			/* r12+r18 */ \
	xmm3 = xmm4; \
	xmm4 -= xmm1;			/* (r4+r26)-(r12+r18) */ \
	xmm1 += xmm3;			/* (r4+r26)+(r12+r18) */ \
\
	xmm5 = xmm7 * CONST2_P975;	/* .975((r2+r28)-(r14+r16)) */ \
	xmm6 = xmm7 * CONST2_P782;	/* .782((r2+r28)-(r14+r16)) */ \
	xmm7 = xmm7 * CONST2_P434;	/* .434((r2+r28)-(r14+r16)) */ \
	xmm5 += xmm4 * CONST2_P782;	/* .975((r2+r28)-(r14+r16))+.782((r4+r26)-(r12+r18)) */ \
	xmm6 -= xmm4 * CONST2_P434;	/* .782((r2+r28)-(r14+r16))-.434((r4+r26)-(r12+r18)) */ \
	xmm7 -= xmm4 * CONST2_P975;	/* .434((r2+r28)-(r14+r16))-.975((r4+r26)-(r12+r18)) */ \
\
	xmm4 = sd5[0];			/* r6 */ \
	xmm4 += sd2[3];			/* r6+r24 */ \
	xmm2 = sd2[1];			/* r10 */ \
	xmm2 += sd5[2];			/* r10+r20 */ \
	xmm3 = xmm4; \
	xmm4 -= xmm2;			/* (r6+r24)-(r10+r20) */ \
	xmm2 += xmm3;			/* (r6+r24)+(r10+r20) */ \
\
	xmm5 += xmm4 * CONST2_P434;	/* .975((r2+r28)-(r14+r16))+.782((r4+r26)-(r12+r18))+.434((r6+r24)-(r10+r20)) */ \
	xmm6 -= xmm4 * CONST2_P975;	/* .782((r2+r28)-(r14+r16))-.434((r4+r26)-(r12+r18))-.975((r6+r24)-(r10+r20)) */ \
	xmm7 += xmm4 * CONST2_P782;	/* .434((r2+r28)-(r14+r16))-.975((r4+r26)-(r12+r18))+.782((r6+r24)-(r10+r20)) */ \
\
	t7 = xmm5;	/* Save real even-cols row #2 */ \
	t8 = xmm6;	/* Save real even-cols row #4 */ \
	t9 = xmm7;	/* Save real even-cols row #6 */ \
\
	xmm3 = sd0[1];			/* r8 */ \
	xmm7 = sd0[3];			/* r22 */ \
	xmm3 -= xmm7;			/* r8-r22 */ \
	xmm7 += sd0[1];			/* r8+r22 */ \
\
	xmm4 = xmm0 + xmm7;					/*     ((r2+r28)+(r14+r16))+(r8+r22) */ \
	xmm5 = xmm0 * CONST2_P901 - xmm7;	/* .901((r2+r28)+(r14+r16))-(r8+r22) */ \
	xmm6 = xmm0 * CONST2_P623 + xmm7;	/* .623((r2+r28)+(r14+r16))+(r8+r22) */ \
	xmm7 = xmm0 * CONST2_P223 - xmm7;	/* .223((r2+r28)+(r14+r16))-(r8+r22) */ \
	xmm4 += xmm1;				/*     ((r2+r28)+(r14+r16))+    ((r4+r26)+(r12+r18))+(r8+r22) */ \
	xmm5 += xmm1 * CONST2_P223;	/* .901((r2+r28)+(r14+r16))+.223((r4+r26)+(r12+r18))-(r8+r22) */ \
	xmm6 -= xmm1 * CONST2_P901;	/* .623((r2+r28)+(r14+r16))-.901((r4+r26)+(r12+r18))+(r8+r22) */ \
	xmm7 -= xmm1 * CONST2_P623;	/* .223((r2+r28)+(r14+r16))-.623((r4+r26)+(r12+r18))-(r8+r22) */ \
	xmm4 += xmm2;				/*     ((r2+r28)+(r14+r16))+    ((r4+r26)+(r12+r18))+    ((r6+r24)+(r10+r20)+(r8+r22)) */ \
	xmm5 -= xmm2 * CONST2_P623;	/* .901((r2+r28)+(r14+r16))+.223((r4+r26)+(r12+r18))-.623((r6+r24)+(r10+r20))-(r8+r22) */ \
	xmm6 -= xmm2 * CONST2_P223;	/* .623((r2+r28)+(r14+r16))-.901((r4+r26)+(r12+r18))-.223((r6+r24)+(r10+r20))+(r8+r22) */ \
	xmm7 += xmm2 * CONST2_P901;	/* .223((r2+r28)+(r14+r16))-.623((r4+r26)+(r12+r18))+.901((r6+r24)+(r10+r20))-(r8+r22) */ \
\
\
	t10 = xmm5;				/* Save real even-cols row #3 */ \
	t11 = xmm6;				/* Save real even-cols row #5 */ \
	t12 = xmm7;				/* Save real even-cols row #7 */ \
\
;	/* Do the even columns for the imaginary results */ \
\
	xmm0 = sd1[0];			/* r2 */ \
	xmm0 -= sd6[3];			/* r2-r28 */ \
	xmm7 = sd6[1];			/* r14 */ \
	xmm7 -= sd1[2];			/* r14-r16 */ \
	xmm2 = xmm7; \
	xmm7 += xmm0;			/* (r2-r28)+(r14-r16) */ \
	xmm0 -= xmm2;			/* (r2-r28)-(r14-r16) */ \
	sd0[1] = xmm4;			/* Save final real #1B (real even-cols row #1) */ \
\
	xmm4 = xmm7 - xmm3;					/* ((r2-r28)+(r14-r16))-(r8-r22) */ \
	xmm5 = xmm7 * CONST2_P223 + xmm3;	/* .223((r2-r28)+(r14-r16))+(r8-r22) */ \
	xmm6 = xmm7 * CONST2_P623 - xmm3;	/* .623((r2-r28)+(r14-r16))-(r8-r22) */ \
	xmm7 = xmm7 * CONST2_P901 + xmm3;	/* .901((r2-r28)+(r14-r16))+(r8-r22) */ \
\
	xmm1 = sd3[0];			/* r4 */ \
	xmm1 -= sd4[3];			/* r4-r26 */ \
	xmm3 = sd4[1];			/* r12 */ \
	xmm3 -= sd3[2];			/* r12-r18 */ \
	xmm2 = xmm3; \
	xmm3 += xmm1;			/* (r4-r26)+(r12-r18) */ \
	xmm1 -= xmm2;			/* (r4-r26)-(r12-r18) */ \
\
	xmm4 -= xmm3;				/* ((r2-r28)+(r14-r16))-((r4-r26)+(r12-r18))-(r8-r22) */ \
	xmm5 += xmm3 * CONST2_P623;	/* .223((r2-r28)+(r14-r16))+.623((r4-r26)+(r12-r18))+(r8-r22) */ \
	xmm6 += xmm3 * CONST2_P901;	/* .623((r2-r28)+(r14-r16))+.901((r4-r26)+(r12-r18))-(r8-r22) */ \
	xmm7 -= xmm3 * CONST2_P223;	/* .901((r2-r28)+(r14-r16))-.223((r4-r26)+(r12-r18))+(r8-r22) */ \
\
	xmm2 = sd5[0];			/* r6 */ \
	xmm2 -= sd2[3];			/* r6-r24 */ \
	xmm3 = sd2[1];			/* r10 */ \
	xmm3 -= sd5[2];			/* r10-r20 */ \
	xmm2 -= xmm3;			/* (r6-r24)-(r10-r20) */ \
	xmm3 += xmm3; \
	xmm3 += xmm2;			/* (r6-r24)+(r10-r20) */ \
\
	xmm4 += xmm3;			/* ((r2-r28)+(r14-r16))-((r4-r26)+(r12-r18))+((r6-r24)+(r10-r20))-(r8-r22) */ \
	xmm5 += xmm3 * CONST2_P901;	/* .223((r2-r28)+(r14-r16))+.623((r4-r26)+(r12-r18))+.901((r6-r24)+(r10-r20))+(r8-r22) */ \
	xmm6 -= xmm3 * CONST2_P223;	/* .623((r2-r28)+(r14-r16))+.901((r4-r26)+(r12-r18))-.223((r6-r24)+(r10-r20))-(r8-r22) */ \
	xmm7 -= xmm3 * CONST2_P623;	/* .901((r2-r28)+(r14-r16))-.223((r4-r26)+(r12-r18))-.623((r6-r24)+(r10-r20))+(r8-r22) */ \
\
	t13 = xmm4;				/* Save imag row #8 */ \
	t14 = xmm5;				/* Save imag even-cols row #2 */ \
	t15 = xmm6;				/* Save imag even-cols row #4 */ \
	t16 = xmm7;				/* Save imag even-cols row #6 */ \
\
	xmm5 = xmm0 * CONST2_P434;	/* .434((r2-r28)-(r14-r16)) */ \
	xmm6 = xmm0 * CONST2_P782;	/* .782((r2-r28)-(r14-r16)) */ \
	xmm0 = xmm0 * CONST2_P975;	/* .975((r2-r28)-(r14-r16)) */ \
	xmm5 += xmm1 * CONST2_P975;	/* .434((r2-r28)-(r14-r16))+.975((r4-r26)-(r12-r18)) */ \
	xmm6 += xmm1 * CONST2_P434;	/* .782((r2-r28)-(r14-r16))+.434((r4-r26)-(r12-r18)) */ \
	xmm0 -= xmm1 * CONST2_P782;	/* .975((r2-r28)-(r14-r16))-.782((r4-r26)-(r12-r18)) */ \
	xmm5 += xmm2 * CONST2_P782;	/* .434((r2-r28)-(r14-r16))+.975((r4-r26)-(r12-r18))+.782((r6-r24)-(r10-r20)) */ \
	xmm6 -= xmm2 * CONST2_P975;	/* .782((r2-r28)-(r14-r16))+.434((r4-r26)-(r12-r18))-.975((r6-r24)-(r10-r20)) */ \
	xmm0 += xmm2 * CONST2_P434;	/* .975((r2-r28)-(r14-r16))-.782((r4-r26)-(r12-r18))+.434((r6-r24)-(r10-r20)) */ \
\
	t17 = xmm5;				/* Save imag even-cols row #3 */ \
	t18 = xmm6;				/* Save imag even-cols row #5 */ \
	t19 = xmm0;				/* Save imag even-cols row #7 */ \
\
	xmm2 = sd4[0];			/* r5 */ \
	xmm2 -= sd3[3];			/* r5-r25 */ \
	xmm6 = sd1[1];			/* r9 */ \
	xmm6 -= sd6[2];			/* r9-r21 */ \
	xmm7 = sd5[1];			/* r13 */ \
	xmm7 -= sd2[2];			/* r13-r17 */ \
\
	xmm0 = t6;				/* Real #8 */ \
	xmm1 = t13;				/* Imag #8 */ \
	xmm5 = sc[6*2+1];		/* cosine/sine */ \
	xmm0 = t6 * xmm5 - xmm1;	/* A8 = R8 * cosine/sine - I8 */ \
	xmm1 = xmm1 * xmm5 + t6;	/* B8 = I8 * cosine/sine + R8 */ \
	xmm0 *= sc[6*2];		/* A8 = A8 * sine (final R8) */ \
	xmm1 *= sc[6*2];		/* B8 = B8 * sine (final I8) */ \
	sd3[2] = xmm0;			/* Save final R8 */ \
	sd3[3] = xmm1;			/* Save final I8 */ \
\
;	/* Do the odd columns for the imag results */ \
\
	xmm0 = xmm2 * CONST2_P782;	/* .782(r5-r25) */ \
	xmm1 = xmm2 * CONST2_P975;	/* .975(r5-r25) */ \
	xmm2 = xmm2 * CONST2_P434;	/* .434(r5-r25) */ \
	xmm0 += xmm6 * CONST2_P975;	/* .782(r5-r25)+.975(r9-r21) */ \
	xmm1 -= xmm6 * CONST2_P434;	/* .975(r5-r25)-.434(r9-r21) */ \
	xmm2 -= xmm6 * CONST2_P782;	/* .434(r5-r25)-.782(r9-r21) */ \
	xmm0 += xmm7 * CONST2_P434;	/* .782(r5-r25)+.975(r9-r21)+.434(r13-r17) */ \
	xmm1 -= xmm7 * CONST2_P782;	/* .975(r5-r25)-.434(r9-r21)-.782(r13-r17) */ \
	xmm2 += xmm7 * CONST2_P975;	/* .434(r5-r25)-.782(r9-r21)+.975(r13-r17) */ \
\
	t20 = xmm0; \
	t21 = xmm1; \
	t22 = xmm2; \
\
	xmm2 = sd2[0];			/* r3 */ \
	xmm2 -= sd5[3];			/* r3-r27 */ \
	xmm6 = sd6[0];			/* r7 */ \
	xmm6 -= sd1[3];			/* r7-r23 */ \
	xmm7 = sd3[1];			/* r11 */ \
	xmm7 -= sd4[2];			/* r11-r19 */ \
	xmm0 = xmm2 * CONST2_P434;	/* .434(r3-r27) */ \
	xmm1 = xmm2 * CONST2_P782;	/* .782(r3-r27) */ \
	xmm2 = xmm2 * CONST2_P975;	/* .975(r3-r27) */ \
	xmm0 += xmm6 * CONST2_P975;	/* .434(r3-r27)+.975(r7-r23) */ \
	xmm1 += xmm6 * CONST2_P434;	/* .782(r3-r27)+.434(r7-r23) */ \
	xmm2 -= xmm6 * CONST2_P782;	/* .975(r3-r27)-.782(r7-r23) */ \
	xmm0 += xmm7 * CONST2_P782;	/* .434(r3-r27)+.975(r7-r23)+.782(r11-r19) */ \
	xmm1 -= xmm7 * CONST2_P975;	/* .782(r3-r27)+.434(r7-r23)-.975(r11-r19) */ \
	xmm2 += xmm7 * CONST2_P434;	/* .975(r3-r27)-.782(r7-r23)+.434(r11-r19) */ \
\
	xmm5 = t20; \
	xmm5 += xmm0;			/* Imag odd-cols row #2 */ \
	xmm0 -= t20;			/* Imag odd-cols row #7 */ \
\
	xmm6 = t21; \
	xmm6 += xmm1;			/* Imag odd-cols row #3 */ \
	xmm1 -= t21;			/* Imag odd-cols row #6 */ \
\
	xmm7 = t22; \
	xmm7 += xmm2;			/* Imag odd-cols row #4 */ \
	xmm2 -= t22;			/* Imag odd-cols row #5 */ \
\
	t20 = xmm5;				/* Imag odd-cols row #2 */ \
	t25 = xmm0;				/* Imag odd-cols row #7 */ \
	t21 = xmm6;				/* Imag odd-cols row #3 */ \
	t24 = xmm1;				/* Imag odd-cols row #6 */ \
	t22 = xmm7;				/* Imag odd-cols row #4 */ \
	t23 = xmm2;				/* Imag odd-cols row #5 */ \
\
;	/* Now combine the even and odd columns then do the post-multiply by twiddle factors. */ \
\
	xmm7 = t0;				/* Real odd-cols row #2 */ \
	xmm6 = t7;				/* Real even-cols row #2 */ \
	xmm7 -= xmm6;			/* Real #14 */ \
	xmm6 += t0;				/* Real #2 */ \
	xmm2 = t14;				/* Imag even-cols row #2 */ \
	xmm5 = t20;				/* Imag odd-cols row #2 */ \
	xmm2 -= xmm5;			/* Imag #14 */ \
	xmm5 += t14;			/* Imag #2 */ \
\
	xmm3 = sc[12*2+1];		/* cosine/sine */ \
	xmm4 = sc[1];			/* cosine/sine */ \
	xmm0 = xmm7 * xmm3 - xmm2;	/* A14 = R14 * cosine/sine - I14 */ \
	xmm1 = xmm6 * xmm4 - xmm5;	/* A2 = R2 * cosine/sine - I2 */ \
	xmm2 = xmm2 * xmm3 + xmm7;	/* B14 = I14 * cosine/sine + R14 */ \
	xmm5 = xmm5 * xmm4 + xmm6;	/* B2 = I2 * cosine/sine + R2 */ \
	xmm0 *= sc[12*2];		/* A14 = A14 * sine (final R14) */ \
	xmm2 *= sc[12*2];		/* B14 = B14 * sine (final I14) */ \
	xmm1 *= sc[0];			/* A2 = A2 * sine (final R2) */ \
	xmm5 *= sc[0];			/* B2 = B2 * sine (final I2) */ \
	sd6[2] = xmm0;			/* Save final R14 */ \
	sd6[3] = xmm2;			/* Save final I14 */ \
	sd0[2] = xmm1;			/* Save final R2 */ \
	sd0[3] = xmm5;			/* Save final I2 */ \
\
	xmm7 = t1;				/* Real odd-cols row #3 */ \
	xmm6 = t10;				/* Real even-cols row #3 */ \
	xmm7 -= xmm6;			/* Real #13 */ \
	xmm6 += t1;				/* Real #3 */ \
	xmm2 = t17;				/* Imag even-cols row #3 */ \
	xmm3 = t21;				/* Imag odd-cols row #3 */ \
	xmm2 -= xmm3;			/* Imag #13 */ \
	xmm3 += t17;			/* Imag #3 */ \
\
	xmm5 = sc[11*2+1];		/* cosine/sine */ \
	xmm4 = sc[2+1];			/* cosine/sine */ \
	xmm0 = xmm7 * xmm5 - xmm2;	/* A13 = R13 * cosine/sine - I13 */ \
	xmm1 = xmm6 * xmm4 - xmm3;	/* A3 = R3 * cosine/sine - I3 */ \
	xmm2 = xmm2 * xmm5 + xmm7;	/* B13 = I13 * cosine/sine + R13 */ \
	xmm3 = xmm3 * xmm4 + xmm6;	/* B3 = I3 * cosine/sine + R3 */ \
	xmm0 *= sc[11*2];		/* A13 = A13 * sine (final R13) */ \
	xmm2 *= sc[11*2];		/* B13 = B13 * sine (final I13) */ \
	xmm1 *= sc[2];			/* A3 = A3 * sine (final R3) */ \
	xmm3 *= sc[2];			/* B3 = B3 * sine (final I3) */ \
	sd6[0] = xmm0;			/* Save final R13 */ \
	sd6[1] = xmm2;			/* Save final I13 */ \
	sd1[0] = xmm1;			/* Save final R3 */ \
	sd1[1] = xmm3;			/* Save final I3 */ \
\
	xmm7 = t2;				/* Real odd-cols row #4 */ \
	xmm6 = t8;				/* Real even-cols row #4 */ \
	xmm7 -= xmm6;			/* Real #12 */ \
	xmm6 += t2;				/* Real #4 */ \
	xmm2 = t15;				/* Imag even-cols row #4 */ \
	xmm3 = t22;				/* Imag odd-cols row #4 */ \
	xmm2 -= xmm3;			/* Imag #12 */ \
	xmm3 += t15;			/* Imag #4 */ \
\
	xmm5 = sc[10*2+1];		/* cosine/sine */ \
	xmm4 = sc[2*2+1];		/* cosine/sine */ \
	xmm0 = xmm7 * xmm5 - xmm2;	/* A12 = R12 * cosine/sine - I12 */ \
	xmm1 = xmm6 * xmm4 - xmm3;	/* A4 = R4 * cosine/sine - I4 */ \
	xmm2 = xmm2 * xmm5 + xmm7;	/* B12 = I12 * cosine/sine + R12 */ \
	xmm3 = xmm3 * xmm4 + xmm6;	/* B4 = I4 * cosine/sine + R4 */ \
	xmm0 *= sc[10*2];		/* A12 = A12 * sine (final R12) */ \
	xmm2 *= sc[10*2];		/* B12 = B12 * sine (final I12) */ \
	xmm1 *= sc[2*2];		/* A4 = A4 * sine (final R4) */ \
	xmm3 *= sc[2*2];		/* B4 = B4 * sine (final I4) */ \
	sd5[2] = xmm0;			/* Save final R12 */ \
	sd5[3] = xmm2;			/* Save final I12 */ \
	sd1[2] = xmm1;			/* Save final R4 */ \
	sd1[3] = xmm3;			/* Save final I4 */ \
\
	xmm7 = t3;				/* Real odd-cols row #5 */ \
	xmm6 = t11;				/* Real even-cols row #5 */ \
	xmm7 -= xmm6;			/* Real #11 */ \
	xmm6 += t3;				/* Real #5 */ \
	xmm2 = t18;				/* Imag even-cols row #5 */ \
	xmm3 = t23;				/* Imag odd-cols row #5 */ \
	xmm2 -= xmm3;			/* Imag #11 */ \
	xmm3 += t18;			/* Imag #5 */ \
\
	xmm5 = sc[9*2+1];		/* cosine/sine */ \
	xmm4 = sc[3*2+1];		/* cosine/sine */ \
	xmm0 = xmm7 * xmm5 - xmm2;	/* A11 = R11 * cosine/sine - I11 */ \
	xmm1 = xmm6 * xmm4 - xmm3;	/* A5 = R5 * cosine/sine - I5 */ \
	xmm2 = xmm2 * xmm5 + xmm7;	/* B11 = I11 * cosine/sine + R11 */ \
	xmm3 = xmm3 * xmm4 + xmm6;	/* B5 = I5 * cosine/sine + R5 */ \
	xmm0 *= sc[9*2];		/* A11 = A11 * sine (final R11) */ \
	xmm2 *= sc[9*2];		/* B11 = B11 * sine (final I11) */ \
	xmm1 *= sc[3*2];		/* A5 = A5 * sine (final R5) */ \
	xmm3 *= sc[3*2];		/* B5 = B5 * sine (final I5) */ \
	sd5[0] = xmm0;			/* Save final R11 */ \
	sd5[1] = xmm2;			/* Save final I11 */ \
	sd2[0] = xmm1;			/* Save final R5 */ \
	sd2[1] = xmm3;			/* Save final I5 */ \
\
	xmm7 = t4;				/* Real odd-cols row #6 */ \
	xmm6 = t9;				/* Real even-cols row #6 */ \
	xmm7 -= xmm6;			/* Real #10 */ \
	xmm6 += t4;				/* Real #6 */ \
	xmm2 = t16;				/* Imag even-cols row #6 */ \
	xmm3 = t24;				/* Imag odd-cols row #6 */ \
	xmm2 -= xmm3;			/* Imag #10 */ \
	xmm3 += t16;			/* Imag #6 */ \
\
	xmm5 = sc[8*2+1];		/* cosine/sine */ \
	xmm4 = sc[4*2+1];		/* cosine/sine */ \
	xmm0 = xmm7 * xmm5 - xmm2;	/* A10 = R10 * cosine/sine - I10 */ \
	xmm1 = xmm6 * xmm4 - xmm3;	/* A6 = R6 * cosine/sine - I6 */ \
	xmm2 = xmm2 * xmm5 + xmm7;	/* B10 = I10 * cosine/sine + R10 */ \
	xmm3 = xmm3 * xmm4 + xmm6;	/* B6 = I6 * cosine/sine + R6 */ \
	xmm0 *= sc[8*2];		/* A10 = A10 * sine (final R10) */ \
	xmm2 *= sc[8*2];		/* B10 = B10 * sine (final I10) */ \
	xmm1 *= sc[4*2];		/* A6 = A6 * sine (final R6) */ \
	xmm3 *= sc[4*2];		/* B6 = B6 * sine (final I6) */ \
	sd4[2] = xmm0;			/* Save final R10 */ \
	sd4[3] = xmm2;			/* Save final I10 */ \
	sd2[2] = xmm1;			/* Save final R6 */ \
	sd2[3] = xmm3;			/* Save final I6 */ \
\
	xmm7 = t5;				/* Real odd-cols row #7 */ \
	xmm6 = t12;				/* Real even-cols row #7 */ \
	xmm7 -= xmm6;			/* Real #9 */ \
	xmm6 += t5;				/* Real #7 */ \
	xmm2 = t19;				/* Imag even-cols row #7 */ \
	xmm3 = t25;				/* Imag odd-cols row #7 */ \
	xmm2 -= xmm3;			/* Imag #9 */ \
	xmm3 += t19;			/* Imag #7 */ \
\
	xmm5 = sc[7*2+1];		/* cosine/sine */ \
	xmm4 = sc[5*2+1];		/* cosine/sine */ \
	xmm0 = xmm7 * xmm5 - xmm2;	/* A9 = R9 * cosine/sine - I9 */ \
	xmm1 = xmm6 * xmm4 - xmm3;	/* A7 = R7 * cosine/sine - I7 */ \
	xmm2 = xmm2 * xmm5 + xmm7;	/* B9 = I9 * cosine/sine + R9 */ \
	xmm3 = xmm3 * xmm4 + xmm6;	/* B7 = I7 * cosine/sine + R7 */ \
	xmm0 *= sc[7*2];		/* A9 = A9 * sine (final R9) */ \
	xmm2 *= sc[7*2];		/* B9 = B9 * sine (final I9) */ \
	xmm1 *= sc[5*2];		/* A7 = A7 * sine (final R7) */ \
	xmm3 *= sc[5*2];		/* B7 = B7 * sine (final I7) */ \
	sd4[0] = xmm0;			/* Save final R9 */ \
	sd4[1] = xmm2;			/* Save final I9 */ \
	sd3[0] = xmm1;			/* Save final R7 */ \
	sd3[1] = xmm3;			/* Save final I7 */ \
\
	srcreg += srcinc; \
}


/* */ \
/* ************************************* 28-reals-last-unfft variants ****************************************** */ \
/* */ \
\
/* These macros produce 28 reals after doing 4.807 levels of the inverse FFT applying */ \
/* the sin/cos multipliers beforehand.  The input is 2 real and 13 complex numbers. */ \
\
/* To calculate a 28-reals inverse FFT, we calculate 28 real values from 28 complex inputs in a brute force way. */ \
/* First we note that the 28 complex values are computed from the 13 complex and 2 real inputs using Hermetian symmetry, thusly: */ \
/* c1 = r1A + 0*i */ \
/* c2 = r2 + i2*i */ \
/* ... */ \
/* c14 = r14 + i14*i */ \
/* c15 = r1B + 0*i */ \
/* c16 = r14 - i14*i */ \
/* ... */ \
/* c28 = r2 - i2*i  */ \
/* */ \
/* The brute force calculations are: */ \
/* */ \
/* c1 + c2 + ... + c28	*  w^-0000000000... */ \
/* c1 + c2 + ... + c28	*  w^-0123456789A... */ \
/* c1 + c2 + ... + c28	*  w^-02468ACE.... */ \
/*    ... */ \
/* c1 + c2 + ... + c28	*  w^-...A987654321 */ \
/* */ \
/* The sin/cos values (w = 28th root of unity) are: */ \
/* w^-1 = .975 - .223i */ \
/* w^-2 = .901 - .434i */ \
/* w^-3 = .782 - .623i */ \
/* w^-4 = .623 - .782i */ \
/* w^-5 = .434 - .901i */ \
/* w^-6 = .223 - .975i */ \
/* w^-7 = 0 - 1i */ \
/* w^-8 = -.223 - .975i */ \
/* w^-9 = -.434 - .901i */ \
/* w^-10 = -.623 - .782i */ \
/* w^-11 = -.782 - .623i */ \
/* w^-12 = -.901 - .434i */ \
/* w^-13 = -.975 - .223i */ \
/* w^-14 = -1 */ \
/* */ \
/* Applying the sin/cos values above, taking advantage of symmetry, and ignoring a lot of multiplies by 2: */ \
/* r1     +(r2+r14)     +(r3+r13)     +(r4+r12)     +(r5+r11)     +(r6+r10)     +(r7+r9) + r8 + r15 */ \
/* r1 +.975(r2-r14) +.901(r3-r13) +.782(r4-r12) +.623(r5-r11) +.434(r6-r10) +.223(r7-r9)      - r15 +.223(i2+i14) +.434*(i3+i13) +.623(i4+i12) +.782(i5+i11) +.901(i6+i10) +.975(i7+i9) + i8 */ \
/* r1 +.901(r2+r14) +.623(r3+r13) +.223(r4+r12) -.223(r5+r11) -.623(r6+r10) -.901(r7+r9) - r8 + r15 +.434(i2-i14) +.782*(i3-i13) +.975(i4-i12) +.975(i5-i11) +.782(i6-i10) +.434(i7-i9) */ \
/* r1 +.782(r2-r14) +.223(r3-r13) -.434(r4-r12) -.901(r5-r11) -.975(r6-r10) -.623(r7-r9)      - r15 +.623(i2+i14) +.975*(i3+i13) +.901(i4+i12) +.434(i5+i11) -.223(i6+i10) -.782(i7+i9) - i8 */ \
/* r1 +.623(r2+r14) -.223(r3+r13) -.901(r4+r12) -.901(r5+r11) -.223(r6+r10) +.623(r7+r9) + r8 + r15 +.782(i2-i14) +.975*(i3-i13) +.434(i4-i12) -.434(i5-i11) -.975(i6-i10) -.782(i7-i9) */ \
/* r1 +.434(r2-r14) -.623(r3-r13) -.975(r4-r12) -.223(r5-r11) +.782(r6-r10) +.901(r7-r9)      - r15 +.901(i2+i14) +.782*(i3+i13) -.223(i4+i12) -.975(i5+i11) -.623(i6+i10) +.434(i7+i9) + i8 */ \
/* r1 +.223(r2+r14) -.901(r3+r13) -.623(r4+r12) +.623(r5+r11) +.901(r6+r10) -.223(r7+r9) - r8 + r15 +.975(i2-i14) +.434*(i3-i13) -.782(i4-i12) -.782(i5-i11) +.434(i6-i10) +.975(i7-i9) */ \
/* r1                   -(r3-r13)                   +(r5-r11)                   -(r7-r9)      - r15     +(i2+i14)                    -(i4+i12)                   +(i6+i10)              - i8 */ \
/* r1 -.223(r2+r14) -.901(r3+r13) +.623(r4+r12) +.623(r5+r11) -.901(r6+r10) -.223(r7+r9) + r8 + r15 +.975(i2-i14) -.434*(i3-i13) -.782(i4-i12) +.782(i5-i11) +.434(i6-i10) -.975(i7-i9) */ \
/* r1 -.434(r2-r14) -.623(r3-r13) +.975(r4-r12) -.223(r5-r11) -.782(r6-r10) +.901(r7-r9)      - r15 +.901(i2+i14) -.782*(i3+i13) -.223(i4+i12) +.975(i5+i11) -.623(i6+i10) -.434(i7+i9) + i8 */ \
/* r1 -.623(r2+r14) -.223(r3+r13) +.901(r4+r12) -.901(r5+r11) +.223(r6+r10) +.623(r7+r9) - r8 + r15 +.782(i2-i14) -.975*(i3-i13) +.434(i4-i12) +.434(i5-i11) -.975(i6-i10) +.782(i7-i9) */ \
/* r1 -.782(r2-r14) +.223(r3-r13) +.434(r4-r12) -.901(r5-r11) +.975(r6-r10) -.623(r7-r9)      - r15 +.623(i2+i14) -.975*(i3+i13) +.901(i4+i12) -.434(i5+i11) -.223(i6+i10) +.782(i7+i9) - i8 */ \
/* r1 -.901(r2+r14) +.623(r3+r13) -.223(r4+r12) -.223(r5+r11) +.623(r6+r10) -.901(r7+r9) + r8 + r15 +.434(i2-i14) -.782*(i3-i13) +.975(i4-i12) -.975(i5-i11) +.782(i6-i10) -.434(i7-i9) */ \
/* r1 -.975(r2-r14) +.901(r3-r13) -.782(r4-r12) +.623(r5-r11) -.434(r6-r10) +.223(r7-r9)      - r15 +.223(i2+i14) -.434*(i3+i13) +.623(i4+i12) -.782(i5+i11) +.901(i6+i10) -.975(i7+i9) + i8 */ \
/* r1     -(r2+r14)     +(r3+r13)     -(r4+r12)     +(r5+r11)     -(r6+r10)     +(r7+r9) - r8 + r15 */ \
/* r1 -.975(r2-r14) +.901(r3-r13) -.782(r4-r12) +.623(r5-r11) -.434(r6-r10) +.223(r7-r9)      - r15 -.223(i2-i14) +.434*(i3-i13) -.623(i4-i12) +.782(i5-i11) -.901(i6-i10) +.975(i7-i9) - i8 */ \
/* ... r17 thru r28 are the same as r12 through r1 but with the sign of the imaginary component changed. */ \
/* */ \
/* Also remember that due to the funny way we do things reals input r1A = r1+r15 and r1B = r1-15 */ \
\
/* Store intermediate results in g->u.xmm.XMM_COL_MULTS (an 8KB buffer used in normalization) */ \
\

#define r7_x14cl_28_reals_last_unfft_preload

#define r7_x14cl_28_reals_last_unfft(srcreg,srcinc,d1,screg,scoff) { \
	vec2f64 temp1[27]; \
	vec2f64 temp2[27]; \
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
	vec2f64 *sd10 = (vec2f64*)(srcreg+d1*10); \
	vec2f64 *sd11 = (vec2f64*)(srcreg+d1*11); \
	vec2f64 *sd12 = (vec2f64*)(srcreg+d1*12); \
	vec2f64 *sd13 = (vec2f64*)(srcreg+d1*13); \
\
	/* Apply the 13 twiddle factors */ \
\
	vec2f64 *sc1 = (vec2f64*)(screg); \
	vec2f64 *sc2 = (vec2f64*)(screg+scoff); \
\
	r7_x13c_twiddle_PART7(sc1,0,temp1) \
	r7_x13c_twiddle_PART(sc1,0,temp1,sd1,sd13,0,24) \
	r7_x13c_twiddle_PART(sc1,0,temp1,sd3,sd11,4,20) \
	r7_x13c_twiddle_PART(sc1,0,temp1,sd5,sd9, 8,16) \
	r7_x13c_twiddle_PART(sc2,1,temp2,sd2,sd12,2,22) \
	r7_x13c_twiddle_PART(sc2,1,temp2,sd4,sd10,6,18) \
	r7_x13c_twiddle_PART(sc2,1,temp2,sd6,sd8,10,14) \
	r7_x13c_twiddle_PART(sc2,1,temp2,sd1,sd13,0,24) \
	r7_x13c_twiddle_PART(sc2,1,temp2,sd3,sd11,4,20) \
	r7_x13c_twiddle_PART(sc2,1,temp2,sd5,sd9, 8,16) \
	r7_x13c_twiddle_PART7(sc2,1,temp2) \
\
	/* Do the 28 reals inverse FFT */ \
\
	r7_x28r_unfft(sd1,sd3,sd5,sd7,sd9,sd11,sd13,sd0[1],sd0[3],temp2); \
	r7_x13c_twiddle_PART(sc1,0,temp1,sd2,sd12,2,22) \
	r7_x13c_twiddle_PART(sc1,0,temp1,sd4,sd10,6,18) \
	r7_x13c_twiddle_PART(sc1,0,temp1,sd6,sd8,10,14) \
	r7_x28r_unfft(sd0,sd2,sd4,sd6,sd8,sd10,sd12,sd0[0],sd0[2],temp1); \
\
	srcreg += srcinc; \
}

/* Apply the 13 twiddle factors.  Also do the first +/- on the reals so that */ \
/* we have a perfect balance of adds and subtracts. */ \

#define r7_x13c_twiddle_PART(SC,SN,tmpreg,SD1,SD2,SCN1,SCN2) \
	xmm1 = SD1[0+SN];		/* R7 */ \
	xmm3 = SD1[2+SN];		/* I7 */ \
	xmm0 = SC[SCN1+1];		/* cosine/sine */ \
	xmm5 = SD2[0+SN];		/* R9 */ \
	xmm7 = SD2[2+SN];		/* I9 */ \
	xmm4 = SC[SCN2+1];		/* cosine/sine */ \
	xmm2 = xmm1 * xmm0 + xmm3;	/* A7 = R7 * cosine/sine + I7 */ \
	xmm6 = xmm5 * xmm4 + xmm7;	/* A9 = R9 * cosine/sine + I9 */ \
	xmm3 = xmm3 * xmm0 - xmm1;	/* B7 = I7 * cosine/sine - R7 */ \
	xmm7 = xmm7 * xmm4 - xmm5;	/* B9 = I9 * cosine/sine - R9 */ \
	xmm1 = SC[SCN1];		/* sine */ \
	xmm5 = SC[SCN2];		/* sine */ \
	xmm2 *= xmm1;			/* R7 = A7 * sine */ \
	xmm6 *= xmm5;			/* R9 = A9 * sine */ \
	xmm3 *= xmm1;			/* I7 = B7 * sine */ \
	xmm7 *= xmm5;			/* I9 = B9 * sine */ \
	xmm1 = xmm2;			/* Copy R7 */ \
	xmm2 += xmm6;			/* R7+R9 */ \
	xmm1 -= xmm6;			/* R7-R9 */ \
	xmm5 = xmm3;			/* Copy I7 */ \
	xmm3 += xmm7;			/* I7+I9 */ \
	xmm5 -= xmm7;			/* I7-I9 */ \
	tmpreg[SCN2] = xmm2;	/* Save R7+R9 */ \
	tmpreg[SCN1] = xmm1;	/* Save R7-R9 */ \
	tmpreg[SCN2+1] = xmm3;	/* Save I7+I9 */ \
	tmpreg[SCN1+1] = xmm5;	/* Save I7-I9 */ \

#define r7_x13c_twiddle_PART7(SC,SN,tmpreg) \
	xmm4 = sd7[0+SN];		/* R8 */ \
	xmm3 = sd7[2+SN];		/* I8 */ \
	xmm0 = SC[13];			/* cosine/sine */ \
	xmm1 = SC[12];			/* sine */ \
	xmm2 = xmm4 * xmm0 + xmm3;	/* A8 = R8 * cosine/sine  + I8 */ \
	xmm3 = xmm3 * xmm0 - xmm4;	/* B8 = I8 * cosine/sine - R8 */ \
	tmpreg[12] = xmm2 * xmm1;	/* R8 = A8 * sine */ \
	tmpreg[13] = xmm3 * xmm1;	/* I8 = B8 * sine */ \



#define r7_x28r_unfft(sd0,sd2,sd4,sd6,sd8,sd10,sd12,r1A,r1B,tmpreg) \
\
;	/* Calculate odd columns derived from real inputs (even rows) */ \
\
	xmm1 = tmpreg[2];			/* r3-r13 */ \
	xmm6 = r1B;					/* r1-r15 */ \
	xmm7 = xmm6 - xmm1;					/* r1-(r3-r13)-r15 */ \
	xmm4 = xmm6 + xmm1 * CONST2_P901;	/* r1+.901(r3-r13)-r15 */ \
	xmm5 = xmm6 + xmm1 * CONST2_P223;	/* r1+.223(r3-r13)-r15 */ \
	xmm0 = xmm6 - xmm1 * CONST2_P623;	/* r1-.623(r3-r13)-r15 */ \
\
	xmm1 = tmpreg[6];			/* r5-r11 */ \
	xmm7 += xmm1;				/* r1-(r3-r13)+(r5-r11)-r15 */ \
	xmm4 += xmm1 * CONST2_P623;	/* r1+.901(r3-r13)+.623(r5-r11)-r15 */ \
	xmm5 -= xmm1 * CONST2_P901;	/* r1+.223(r3-r13)-.901(r5-r11)-r15 */ \
	xmm0 -= xmm1 * CONST2_P223;	/* r1-.623(r3-r13)-.223(r5-r11)-r15 */ \
\
	xmm1 = tmpreg[10];			/* r7-r9 */ \
	xmm7 -= xmm1;				/* r1-(r3-r13)+(r5-r11)-(r7-r9)-r15 */ \
	xmm4 += xmm1 * CONST2_P223;	/* r1+.901(r3-r13)+.623(r5-r11)+.223(r7-r9)-r15 */ \
	xmm5 -= xmm1 * CONST2_P623;	/* r1+.223(r3-r13)-.901(r5-r11)-.623(r7-r9)-r15 */ \
	xmm0 += xmm1 * CONST2_P901;	/* r1-.623(r3-r13)-.223(r5-r11)+.901(r7-r9)-r15 */ \
\
	tmpreg[26] = xmm7;		/* Save odd-real-cols row #8 (also is real-cols row #8) */ \
\
;	/* Calculate even columns derived from real inputs (even rows) */ \
;	/* From above, odd-real-col rols rows #2,4,6 are in xmm4, xmm5, xmm0 */ \
\
	xmm3 = tmpreg[0];			/* r2-r14 */ \
	xmm1 = xmm3 * CONST2_P975;	/* .975(r2-r14) */ \
	xmm2 = xmm3 * CONST2_P782;	/* .782(r2-r14) */ \
	xmm3 = xmm3 * CONST2_P434;	/* .434(r2-r14) */ \
\
	xmm6 = tmpreg[4];			/* r4-r12 */ \
	xmm1 += xmm6 * CONST2_P782;	/* .975(r2-r14)+.782(r4-r12) */ \
	xmm2 -= xmm6 * CONST2_P434;	/* .782(r2-r14)-.434(r4-r12) */ \
	xmm3 -= xmm6 * CONST2_P975;	/* .434(r2-r14)-.975(r4-r12) */ \
\
	xmm6 = tmpreg[8];			/* r6-r10 */ \
	xmm1 += xmm6 * CONST2_P434;	/* .975(r2-r14)+.782(r4-r12)+.434(r6-r10) */ \
	xmm2 -= xmm6 * CONST2_P975;	/* .782(r2-r14)-.434(r4-r12)-.975(r6-r10) */ \
	xmm3 += xmm6 * CONST2_P782;	/* .434(r2-r14)-.975(r4-r12)+.782(r6-r10) */ \
\
;	/* Combine even and odd columns (even rows) */ \
\
	xmm7 = xmm4;			/* Copy odd-real-cols row #2 */ \
	xmm4 -= xmm1;			/* real-cols row #14 (odd#2 - even#2) */ \
	xmm1 += xmm7;			/* real-cols row #2 (odd#2 + even#2) */ \
\
	xmm6 = xmm5;			/* Copy odd-real-cols row #4 */ \
	xmm5 -= xmm2;			/* real-cols row #12 (odd#4 - even#4) */ \
	xmm2 += xmm6;			/* real-cols row #4 (odd#4 + even#4) */ \
\
	xmm7 = xmm0;			/* Copy odd-real-cols row #6 */ \
	xmm0 -= xmm3;			/* real-cols row #10 (odd#6 - even#6) */ \
	xmm3 += xmm7;			/* real-cols row #6 (odd#6 + even#6) */ \
\
	tmpreg[2] = xmm4;		/* Save real-cols row #14 */ \
	tmpreg[0] = xmm1;			/* Save real-cols row #2 */ \
	tmpreg[6] = xmm5;		/* Save real-cols row #12 */ \
	tmpreg[4] = xmm2;		/* Save real-cols row #4 */ \
	tmpreg[10] = xmm0;		/* Save real-cols row #10 */ \
	tmpreg[8] = xmm3;		/* Save real-cols row #6 */ \
\
;	/* Calculate even columns derived from real inputs (odd rows) */ \
\
	xmm2 = tmpreg[24];			/* r2+r14 */ \
	xmm4 = tmpreg[12];			/* r8 */ \
	xmm0 = xmm2 + xmm4;					/*     (r2+r14)+r8 */ \
	xmm5 = xmm2 * CONST2_P901 - xmm4;	/* .901(r2+r14)-r8 */ \
	xmm6 = xmm2 * CONST2_P623 + xmm4;	/* .623(r2+r14)+r8 */ \
	xmm7 = xmm2 * CONST2_P223 - xmm4;	/* .223(r2+r14)-r8 */ \
\
	xmm1 = tmpreg[20];		/* r4+r12 */ \
	xmm0 += xmm1;				/*     (r2+r14)+    (r4+r12)+r8 */ \
	xmm5 += xmm1 * CONST2_P223;	/* .901(r2+r14)+.223(r4+r12)-r8 */ \
	xmm6 -= xmm1 * CONST2_P901;	/* .623(r2+r14)-.901(r4+r12)+r8 */ \
	xmm7 -= xmm1 * CONST2_P623;	/* .223(r2+r14)-.623(r4+r12)-r8 */ \
\
	xmm2 = tmpreg[16];			/* r6+r10 */ \
	xmm0 += xmm2;				/*     (r2+r14)+    (r4+r12)+    (r6+r10)+r8 */ \
	xmm5 -= xmm2 * CONST2_P623;	/* .901(r2+r14)+.223(r4+r12)-.623(r6+r10)-r8 */ \
	xmm6 -= xmm2 * CONST2_P223;	/* .623(r2+r14)-.901(r4+r12)-.223(r6+r10)+r8 */ \
	xmm7 += xmm2 * CONST2_P901;	/* .223(r2+r14)-.623(r4+r12)+.901(r6+r10)-r8 */ \
\
\
	tmpreg[12] = xmm0;		/* Save even-real-cols row #1	;; We could save a few loads and stores */ \
	tmpreg[16] = xmm5;		/* Save even-real-cols row #3	;; if two of these registers were left */ \
	tmpreg[20] = xmm6;		/* Save even-real-cols row #5	;; unchanged through the next section */ \
	tmpreg[24] = xmm7;		/* Save even-real-cols row #7 */ \
\
;	/* Calculate odd columns derived from real inputs (odd rows) */ \
\
	xmm3 = tmpreg[22];			/* r3+r13 */ \
	xmm7 = r1A;					/* r1+r15 */ \
	xmm0 = xmm7 + xmm3;					/* r1+    (r3+r13)+r15 */ \
	xmm5 = xmm7 + xmm3 * CONST2_P623;	/* r1+.623(r3+r13)+r15 */ \
	xmm4 = xmm7 - xmm3 * CONST2_P223;	/* r1-.223(r3+r13)+r15 */ \
	xmm7 = xmm7 - xmm3 * CONST2_P901;	/* r1-.901(r3+r13)+r15 */ \
	\
	xmm1 = tmpreg[18];			/* r5+r11 */ \
	xmm0 += xmm1;				/* r1+(r3+r13)+(r5+r11)+r15 */ \
	xmm5 -= xmm1 * CONST2_P223;	/* r1+.623(r3+r13)-.223(r5+r11)+r15 */ \
	xmm4 -= xmm1 * CONST2_P901;	/* r1-.223(r3+r13)-.901(r5+r11)+r15 */ \
	xmm7 += xmm1 * CONST2_P623;	/* r1-.901(r3+r13)+.623(r5+r11)+r15 */ \
\
	xmm2 = tmpreg[14];			/* r7+r9 */ \
	xmm0 += xmm2;				/* r1+(r3+r13)+(r5+r11)+(r7+r9)+r15 */ \
	xmm5 -= xmm2 * CONST2_P901;	/* r1+.623(r3+r13)-.223(r5+r11)-.901(r7+r9)+r15 */ \
	xmm4 += xmm2 * CONST2_P623;	/* r1-.223(r3+r13)-.901(r5+r11)+.623(r7+r9)+r15 */ \
	xmm7 -= xmm2 * CONST2_P223;	/* r1-.901(r3+r13)+.623(r5+r11)-.223(r7+r9)+r15 */ \
\
;	/* Combine even and odd columns (odd rows) */ \
\
	xmm1 = tmpreg[12];		/* even-real-cols row #1 */ \
	xmm1 += xmm0;			/* real-cols row #1 (and final R1) */ \
	xmm0 -= tmpreg[12];		/* real-cols row #15 (and final R15) */ \
	sd0[0] = xmm1;			/* Save final R1 */ \
	sd0[2] = xmm0;			/* Save final R15 */ \
\
	xmm2 = tmpreg[16];		/* even-real-cols row #3 */ \
	xmm2 += xmm5;			/* real-cols row #3 */ \
	xmm5 -= tmpreg[16];		/* real-cols row #13 */ \
	tmpreg[12] = xmm2;		/* Save real-cols row #3 */ \
	tmpreg[14] = xmm5;		/* Save real-cols row #13 */ \
\
	xmm3 = tmpreg[20];		/* even-real-cols row #5 */ \
	xmm3 += xmm4;			/* real-cols row #5 */ \
	xmm4 -= tmpreg[20]; 		/* real-cols row #11 */ \
	tmpreg[16] = xmm3;		/* Save real-cols row #5 */ \
	tmpreg[18] = xmm4;		/* Save real-cols row #11 */ \
\
	xmm6 = tmpreg[24];		/* even-real-cols row #7 */ \
	xmm6 += xmm7;			/* real-cols row #7 */ \
	xmm7 -= tmpreg[24];		/* real-cols row #9 */ \
	tmpreg[20] = xmm6;		/* Save real-cols row #7 */ \
	tmpreg[24] = xmm7;		/* Save real-cols row #9 */ \
\
;	/* Calculate even columns derived from imaginary inputs (even rows) */ \
\
	xmm3 = tmpreg[25];			/* i2+i14 */ \
	xmm4 = tmpreg[13];			/* i8 */ \
	xmm0 = xmm3 - xmm4;			/* (i2+i14)-i8 */ \
	xmm5 = xmm3 * CONST2_P223 + xmm4;	/* .223(i2+i14)+i8 */ \
	xmm6 = xmm3 * CONST2_P623 - xmm4;	/* .623(i2+i14)-i8 */ \
	xmm7 = xmm3 * CONST2_P901 + xmm4;	/* .901(i2+i14)+i8 */ \
\
	xmm1 = tmpreg[21];			/* i4+i12 */ \
	xmm0 -= xmm1;				/* (i2+i14)-(i4+i12)-i8 */ \
	xmm5 += xmm1 * CONST2_P623;	/* .223(i2+i14)+.623(i4+i12)+i8 */ \
	xmm6 += xmm1 * CONST2_P901;	/* .623(i2+i14)+.901(i4+i12)-i8 */ \
	xmm7 -= xmm1 * CONST2_P223;	/* .901(i2+i14)-.223(i4+i12)+i8 */ \
\
	xmm2 = tmpreg[17];			/* i6+i10 */ \
	xmm0 += xmm2;				/* (i2+i14)-(i4+i12)+(i6+i10)-i8 */ \
	xmm5 += xmm2 * CONST2_P901;	/* .223(i2+i14)+.623(i4+i12)+.901(i6+i10)+i8 */ \
	xmm6 -= xmm2 * CONST2_P223;	/* .623(i2+i14)+.901(i4+i12)-.223(i6+i10)-i8 */ \
	xmm7 -= xmm2 * CONST2_P623;	/* .901(i2+i14)-.223(i4+i12)-.623(i6+i10)+i8 */ \
\
;	/* Combine real and imaginary data for row #8 */ \
\
	xmm4 = tmpreg[26];		/* Load real-cols row #8 */ \
	xmm4 -= xmm0;			/* final R22 */ \
	xmm0 += tmpreg[26];		/* final R8 */ \
	sd0[3] = xmm4;			/* Save R22 */ \
	sd0[1] = xmm0;			/* Save R8 */ \
\
;	/* Calculate odd columns derived from imaginary inputs (even rows) */ \
;	/* From above, even-imag-cols row #2,4,6 are in xmm5, xmm6, xmm7 */ \
\
	xmm4 = tmpreg[23];			/* i3+i13 */ \
	xmm2 = xmm4 * CONST2_P434;	/* .434(i3+i13) */ \
	xmm3 = xmm4 * CONST2_P975;	/* .975(i3+i13) */ \
	xmm4 = xmm4 * CONST2_P782;	/* .782(i3+i13) */ \
\
	xmm0 = tmpreg[19];			/* i5+i11 */ \
	xmm2 += xmm0 * CONST2_P782;	/* .434(i3+i13)+.782(i5+i11) */ \
	xmm3 += xmm0 * CONST2_P434;	/* .975(i3+i13)+.434(i5+i11) */ \
	xmm4 -= xmm0 * CONST2_P975;	/* .782(i3+i13)-.975(i5+i11) */ \
\
	xmm0 = tmpreg[15];		/* i7+i9 */ \
	xmm2 += xmm0 * CONST2_P975;	/* .434(i3+i13)+.782(i5+i11)+.975(i7+i9) */ \
	xmm3 -= xmm0 * CONST2_P782;	/* .975(i3+i13)+.434(i5+i11)-.782(i7+i9) */ \
	xmm4 += xmm0 * CONST2_P434;	/* .782(i3+i13)-.975(i5+i11)+.434(i7+i9) */ \
\
;	/* Combine even and odd columns, then real and imag data (even rows) */ \
\
	xmm0 = xmm5;			/* Copy even-imag-cols row #2 */ \
	xmm5 -= xmm2;			/* imag-cols row #14 (even#2 - odd#2) */ \
	xmm2 += xmm0;			/* imag-cols row #2 (even#2 + odd#2) */ \
	xmm1 = xmm6;			/* Copy even-imag-cols row #4 */ \
	xmm6 -= xmm3;			/* imag-cols row #12 (even#4 - odd#4) */ \
	xmm3 += xmm1;			/* imag-cols row #4 (even#4 + odd#4) */ \
	xmm0 = xmm7;			/* Copy even-imag-cols row #6 */ \
	xmm7 -= xmm4;			/* imag-cols row #10 (even#6 - odd#6) */ \
	xmm4 += xmm0;			/* imag-cols row #6 (even#6 + odd#6) */ \
\
	xmm0 = tmpreg[2];		/* Load real-cols row #14 */ \
	xmm0 -= xmm5;			/* final R16 */ \
	xmm5 += tmpreg[2];		/* final R14 */ \
	xmm1 = tmpreg[0];		/* Load real-cols row #2 */ \
	xmm1 -= xmm2;			/* final R28 */ \
	xmm2 += tmpreg[0];		/* final R2 */ \
	sd2[2] = xmm0;			/* Save R16 */ \
	sd12[1] = xmm5;			/* Save R14 */ \
	sd12[3] = xmm1;			/* Save R28 */ \
	sd2[0] = xmm2;			/* Save R2 */ \
\
	xmm0 = tmpreg[6];		/* Load real-cols row #12 */ \
	xmm0 -= xmm6;			/* final R18 */ \
	xmm6 += tmpreg[6];		/* final R12 */ \
	xmm1 = tmpreg[4];		/* Load real-cols row #4 */ \
	xmm1 -= xmm3;			/* final R26 */ \
	xmm3 += tmpreg[4];		/* final R4 */ \
	sd6[2] = xmm0;			/* Save R18 */ \
	sd8[1] = xmm6;			/* Save R12 */ \
	sd8[3] = xmm1;			/* Save R26 */ \
	sd6[0] = xmm3;			/* Save R4 */ \
\
	xmm0 = tmpreg[10];		/* Load real-cols row #10 */ \
	xmm0 -= xmm7;			/* final R20 */ \
	xmm7 += tmpreg[10];		/* final R10 */ \
	xmm1 = tmpreg[8];		/* Load real-cols row #6 */ \
	xmm1 -= xmm4;			/* final R24 */ \
	xmm4 += tmpreg[8];		/* final R6 */ \
	sd10[2] = xmm0;			/* Save R20 */ \
	sd4[1] = xmm7;			/* Save R10 */ \
	sd4[3] = xmm1;			/* Save R24 */ \
	sd10[0] = xmm4;			/* Save R6 */ \
\
;	/* Calculate even columns derived from imaginary inputs (odd rows) */ \
\
	xmm7 = tmpreg[1];			/* i2-i14 */ \
	xmm5 = xmm7 * CONST2_P434;	/* .434(i2-i14) */ \
	xmm6 = xmm7 * CONST2_P782;	/* .782(i2-i14) */ \
	xmm7 = xmm7 * CONST2_P975;	/* .975(i2-i14) */ \
\
	xmm4 = tmpreg[5];			/* i4-i12 */ \
	xmm5 += xmm4 * CONST2_P975;	/* .434(i2-i14)+.975(i4-i12) */ \
	xmm6 += xmm4 * CONST2_P434;	/* .782(i2-i14)+.434(i4-i12) */ \
	xmm7 -= xmm4 * CONST2_P782;	/* .975(i2-i14)-.782(i4-i12) */ \
\
	xmm4 = tmpreg[9];			/* i6-i10 */ \
	xmm5 += xmm4 * CONST2_P782;	/* .434(i2-i14)+.975(i4-i12)+.782(i6-i10) */ \
	xmm6 -= xmm4 * CONST2_P975;	/* .782(i2-i14)+.434(i4-i12)-.975(i6-i10) */ \
	xmm7 += xmm4 * CONST2_P434;	/* .975(i2-i14)-.782(i4-i12)+.434(i6-i10) */ \
\
;	/* Calculate odd columns derived from imaginary inputs (odd rows) */ \
;	/* From above, even-imag-cols row #3,5,7 are in xmm5,xmm6,xmm7 */ \
\
	xmm4 = tmpreg[3];			/* i3-i13 */ \
	xmm2 = xmm4 * CONST2_P782;	/* .782(i3-i13) */ \
	xmm3 = xmm4 * CONST2_P975;	/* .975(i3-i13) */ \
	xmm4 = xmm4 * CONST2_P434;	/* .434(i3-i13) */ \
\
	xmm0 = tmpreg[7];			/* i5-i11 */ \
	xmm2 += xmm0 * CONST2_P975;	/* .782(i3-i13)+.975(i5-i11) */ \
	xmm3 -= xmm0 * CONST2_P434;	/* .975(i3-i13)-.434(i5-i11) */ \
	xmm4 -= xmm0 * CONST2_P782;	/* .434(i3-i13)-.782(i5-i11) */ \
\
	xmm0 = tmpreg[11];			/* i7-i9 */ \
	xmm2 += xmm0 * CONST2_P434;	/* .782(i3-i13)+.975(i5-i11)+.434(i7-i9) */ \
	xmm3 -= xmm0 * CONST2_P782;	/* .975(i3-i13)-.434(i5-i11)-.782(i7-i9) */ \
	xmm4 += xmm0 * CONST2_P975;	/* .434(i3-i13)-.782(i5-i11)+.975(i7-i9) */ \
\
;	/* Combine even and odd columns, then real and imag data (odd rows) */ \
\
	xmm0 = xmm5;			/* Copy even-imag-cols row #3 */ \
	xmm5 -= xmm2;			/* imag-cols row #13 (even#3 - odd#3) */ \
	xmm2 += xmm0;			/* imag-cols row #3 (even#3 + odd#3) */ \
	xmm1 = xmm6;			/* Copy even-imag-cols row #5 */ \
	xmm6 -= xmm3;			/* imag-cols row #11 (even#5 - odd#5) */ \
	xmm3 += xmm1;			/* imag-cols row #5 (even#5 + odd#5) */ \
	xmm0 = xmm7;			/* Copy even-imag-cols row #7 */ \
	xmm7 -= xmm4;			/* imag-cols row #9 (even#7 - odd#7) */ \
	xmm4 += xmm0;			/* imag-cols row #7 (even#7 + odd#7) */ \
\
	xmm0 = tmpreg[14];		/* Load real-cols row #13 */ \
	xmm0 -= xmm5;			/* final R17 */ \
	xmm5 += tmpreg[14];		/* final R13 */ \
	xmm1 = tmpreg[12];		/* Load real-cols row #3 */ \
	xmm1 -= xmm2;			/* final R27 */ \
	xmm2 += tmpreg[12];		/* final R3 */ \
	sd4[2] = xmm0;			/* Save R17 */ \
	sd10[1] = xmm5;			/* Save R13 */ \
	sd10[3] = xmm1;			/* Save R27 */ \
	sd4[0] = xmm2;			/* Save R3 */ \
\
	xmm0 = tmpreg[18];		/* Load real-cols row #11 */ \
	xmm0 -= xmm6;			/* final R19 */ \
	xmm6 += tmpreg[18];		/* final R11 */ \
	xmm1 = tmpreg[16];		/* Load real-cols row #5 */ \
	xmm1 -= xmm3;			/* final R25 */ \
	xmm3 += tmpreg[16];		/* final R5 */ \
	sd8[2] = xmm0;			/* Save R19 */ \
	sd6[1] = xmm6;			/* Save R11 */ \
	sd6[3] = xmm1;			/* Save R25 */ \
	sd8[0] = xmm3;			/* Save R5 */ \
\
	xmm0 = tmpreg[24];		/* Load real-cols row #9 */ \
	xmm0 -= xmm7;			/* final R21 */ \
	xmm7 += tmpreg[24];		/* final R9 */ \
	xmm1 = tmpreg[20];		/* Load real-cols row #7 */ \
	xmm1 -= xmm4;			/* final R23 */ \
	xmm4 += tmpreg[20];		/* final R7 */ \
	sd12[2] = xmm0;			/* Save R21 */ \
	sd2[1] = xmm7;			/* Save R9 */ \
	sd2[3] = xmm1;			/* Save R23 */ \
	sd12[0] = xmm4;			/* Save R7 */
