#pragma once
	//
#include "hg64.h"


#define s2cl_eight_reals_first_fft(srcreg,srcinc,d1); \
	shuffle_load(xmm0,xmm2,xptr(srcreg+rbx),xptr(srcreg+16+rbx)); /* R1,R3 */ \
	shuffle_load(xmm1,xmm3,xptr(srcreg+d1+rbx),xptr(srcreg+d1+16+rbx)); /* R2,R4 */ \
	shuffle_load(xmm4,xmm6,xptr(srcreg+32+rbx),xptr(srcreg+48+rbx)); /* R5,R7 */ \
	shuffle_load(xmm5,xmm7,xptr(srcreg+d1+32+rbx),xptr(srcreg+d1+48+rbx)); /* R6,R8 */ \
	x8r_fft; \
	xptr(srcreg) = xmm7; \
	xptr(srcreg+16) = xmm6; \
	xptr(srcreg+32) = xmm4; \
	xptr(srcreg+48) = xmm5; \
	xptr(srcreg+d1) = xmm1; \
	xptr(srcreg+d1+16) = xmm3; \
	xptr(srcreg+d1+32) = xmm0; \
	xptr(srcreg+d1+48) = xmm2; \
	srcreg += srcinc; \

#define x2cl_eight_reals_first_fft(srcreg,srcinc,d1); \
	xmm0 = xptr(srcreg+rbx); \
	xmm1 = xptr(srcreg+d1+rbx); \
	xmm2 = xptr(srcreg+16+rbx); \
	xmm3 = xptr(srcreg+d1+16+rbx); \
	xmm4 = xptr(srcreg+32+rbx); \
	xmm5 = xptr(srcreg+d1+32+rbx); \
	xmm6 = xptr(srcreg+48+rbx); \
	xmm7 = xptr(srcreg+d1+48+rbx); \
	x8r_fft; \
	xptr(srcreg) = xmm7; \
	xptr(srcreg+16) = xmm6; \
	xptr(srcreg+32) = xmm4; \
	xptr(srcreg+48) = xmm5; \
	xptr(srcreg+d1) = xmm1; \
	xptr(srcreg+d1+16) = xmm3; \
	xptr(srcreg+d1+32) = xmm0; \
	xptr(srcreg+d1+48) = xmm2; \
	srcreg += srcinc; \

#define x2cl_eight_reals_fft(srcreg,srcinc,d1); \
	xmm0 = xptr(srcreg); \
	xmm1 = xptr(srcreg+d1); \
	xmm2 = xptr(srcreg+16); \
	xmm3 = xptr(srcreg+d1+16); \
	xmm4 = xptr(srcreg+32); \
	xmm5 = xptr(srcreg+d1+32); \
	xmm6 = xptr(srcreg+48); \
	xmm7 = xptr(srcreg+d1+48); \
	x8r_fft; \
	xptr(srcreg) = xmm7; \
	xptr(srcreg+16) = xmm6; \
	xptr(srcreg+32) = xmm4; \
	xptr(srcreg+48) = xmm5; \
	xptr(srcreg+d1) = xmm1; \
	xptr(srcreg+d1+16) = xmm3; \
	xptr(srcreg+d1+32) = xmm0; \
	xptr(srcreg+d1+48) = xmm2; \
	srcreg += srcinc; \

#define  g2cl_eight_reals_first_fft(srcreg,srcinc,d1,dstreg,dstinc,e1); \
	xmm0 = xptr(srcreg+rbx); \
	xmm1 = xptr(srcreg+d1+rbx); \
	xmm2 = xptr(srcreg+16+rbx); \
	xmm3 = xptr(srcreg+d1+16+rbx); \
	xmm4 = xptr(srcreg+32+rbx); \
	xmm5 = xptr(srcreg+d1+32+rbx); \
	xmm6 = xptr(srcreg+48+rbx); \
	xmm7 = xptr(srcreg+d1+48+rbx); \
	srcreg += srcinc; \
	x8r_fft; \
	xptr(dstreg) = xmm7; \
	xptr(dstreg+16) = xmm6; \
	xptr(dstreg+32) = xmm4; \
	xptr(dstreg+48) = xmm5; \
	xptr(dstreg+e1) = xmm1; \
	xptr(dstreg+e1+16) = xmm3; \
	xptr(dstreg+e1+32) = xmm0; \
	xptr(dstreg+e1+48) = xmm2; \
	dstreg += dstinc; \



#define x4c_unfft(r1, r2, r3, r4, r5, r6, r7, r8) { \
	vec2f64 tmp1,tmp2,tmp3,tmp4; \
	r1 -= r3;			/* new R2 = R1 - R2 */ \
	r3 *= 2; \
	r3 += r1;			/* new R1 = R1 + R2 */ \
	r6 -= r8;			/* new R4 = I3 - I4 */ \
	r8 *= 2; \
	r8 += r6;			/* new I3 = I3 + I4 */ \
	r2 -= r4;			/* new I2 = I1 - I2 */ \
	r4 *= 2; \
	r4 += r2;			/* new I1 = I1 + I2 */ \
	r7 -= r5;			/* new I4 = R4 - R3 */ \
	r5 *= 2; \
	r5 += r7;			/* new R3 = R3 + R4 */ \
	r1 -= r6;			/* R2 = R2 - R4 (new R4) */ \
	r6 *= 2;			/* R4 = R4 * 2 */ \
	r6 += r1;			/* R4 = R2 + R4 (new R2) */ \
	r2 -= r7;			/* I2 = I2 - I4 (new I4) */ \
	r7 *= 2;			/* I4 = I4 * 2 */ \
	r7 += r2;			/* I4 = I2 + I4 (new I2) */ \
	r3 -= r5;			/* R1 = R1 - R3 (new R3) */ \
	tmp1 = r6;		/* Save new R2 */ \
	r2 *= xptr(rdi+64);		/* B4 = new I4 * sine */ \
	tmp2 = r3;		/* Save new R3 */ \
	r1 *= xptr(rdi+64);		/* A4 = new R4 * sine */ \
	r6 *= xptr(rdi+0+16);		/* A2 = new R2 * cosine/sine */ \
	tmp3 = r2;		/* Save B4 */ \
	r3 *= xptr(rdi+32+16);		/* A3 = new R3 * cosine/sine */ \
	r4 -= r8;			/* I1 = I1 - I3 (new I3) */ \
	r2 *= xptr(rdi+64+16);		/* C4 = B4 * cosine/sine */ \
	r8 *= 2;			/* I3 = I3 * 2 */ \
	tmp4 = r4;		/* Save I1 */ \
	r6 += r7;			/* A2 = A2 + new I2 */ \
	r7 *= xptr(rdi+0+16);		/* B2 = new I2 * cosine/sine */ \
	r3 += r4;			/* A3 = A3 + new I3 */ \
	r4 *= xptr(rdi+32+16);		/* B3 = new I3 * cosine/sine */ \
	r5 *= 2;			/* R3 = R3 * 2 */ \
	r2 -= r1;			/* C4 = C4 - A4 (final I4) */ \
	r1 *= xptr(rdi+64+16);		/* A4 = A4 * cosine/sine */ \
	r7 -= tmp1;		/* B2 = B2 - new R2 */ \
	r6 *= xptr(rdi+0);		/* A2 = A2 * sine (final R2) */ \
	r4 -= tmp2;		/* B3 = B3 - new R3 */ \
	r3 *= xptr(rdi+32);		/* A3 = A3 * sine (final R3) */ \
	r8 += tmp4;		/* I3 = I1 + I3 (new & final I1) */ \
	r7 *= xptr(rdi+0);		/* B2 = B2 * sine (final I2) */ \
	r1 += tmp3;		/* A4 = B4 + A4 (final R4) */ \
	r4 *= xptr(rdi+32);		/* B3 = B3 * sine (final I3) */ \
	r5 += tmp2;		/* R3 = R1 + R3 (new & final R1) */ \
}


#define s2cl_four_complex_with_mulf(srcreg,srcinc,d1) \
	xmm3 = xptr(srcreg+rbx);	/* R1 */ \
	xmm7 = xptr(srcreg+16+rbx);	/* R2 */ \
	xmm2 = xptr(srcreg+32+rbx);	/* R3 */ \
	xmm6 = xptr(srcreg+48+rbx);	/* R4 */ \
	xmm0 = xptr(srcreg+d1+rbx);	/* R5 */ \
	xmm1 = xptr(srcreg+d1+16+rbx);/* R6 */ \
	xmm5 = xptr(srcreg+d1+32+rbx);/* R7 */ \
	xmm4 = xptr(srcreg+d1+48+rbx);/* R8 */ \
	\
	xp4c_mulf(xmm3, xmm7, xmm2, xmm6, xmm0, xmm1, xmm5, xmm4, srcreg, srcreg+16, srcreg+32, srcreg+48, srcreg+d1, srcreg+d1+16, srcreg+d1+32, srcreg+d1+48); \
	\
	x4c_unfft(xmm3, xmm7, xmm2, xmm6, xmm0, xmm1, xmm5, xmm4); \
	\
	shuffle_store(xptr(srcreg), xptr(srcreg+d1), xmm0, xmm1); /* Save R1 and R2 */ \
	shuffle_store(xptr(srcreg+16), xptr(srcreg+d1+16), xmm2, xmm3); /* Save R3 and R4 */ \
	shuffle_store(xptr(srcreg+32), xptr(srcreg+d1+32), xmm4, xmm5); /* Save R5 and R6 */ \
	shuffle_store(xptr(srcreg+48), xptr(srcreg+d1+48), xmm6, xmm7); /* Save R7 and R8 */ \
	\
	srcreg += srcinc; \


#define s2cl_eight_reals_fft_2_final(srcreg,srcinc,d1); \
	xmm0[0] = f64ptr(srcreg);	/* R1 */ \
	xmm1[0] = f64ptr(srcreg+8);	/* R2 */ \
	xmm2[0] = f64ptr(srcreg+d1);	/* R3 */ \
	xmm3[0] = f64ptr(srcreg+d1+8);	/* R4 */ \
	xmm4[0] = f64ptr(srcreg+16);	/* R5 */ \
	xmm5[0] = f64ptr(srcreg+24);	/* R6 */ \
	xmm6[0] = f64ptr(srcreg+d1+16);	/* R7 */ \
	xmm7[0] = f64ptr(srcreg+d1+24);	/* R8 */ \
	xs8r_fft(xmm0[0], xmm1[0], xmm2[0], xmm3[0], xmm4[0], xmm5[0], xmm6[0], xmm7[0]); \
	f64ptr(srcreg) = xmm3[0];	/* Save R1 */ \
	f64ptr(srcreg+8) = xmm2[0];	/* Save R2 */ \
	f64ptr(srcreg+16) = xmm0[0];	/* Save R3 */ \
	f64ptr(srcreg+24) = xmm1[0];	/* Save R4 */ \
	xmm0[0] = f64ptr(srcreg+32);	/* R1 */ \
	xmm1[0] = f64ptr(srcreg+40);	/* R2 */ \
	xmm2[0] = f64ptr(srcreg+48);	/* R5 */ \
	xmm3[0] = f64ptr(srcreg+56);	/* R6 */ \
	f64ptr(srcreg+32) = xmm5[0];	/* Save R5 */ \
	f64ptr(srcreg+40) = xmm7[0];	/* Save R6 */ \
	f64ptr(srcreg+48) = xmm4[0];	/* Save R7 */ \
	f64ptr(srcreg+56) = xmm6[0];	/* Save R8 */ \
	xmm4[0] = f64ptr(srcreg+d1+32);	/* R3 */ \
	xmm5[0] = f64ptr(srcreg+d1+40);	/* R4 */ \
	xmm6[0] = f64ptr(srcreg+d1+48);	/* R7 */ \
	xs4c_fft(xmm0[0], xmm1[0], xmm4[0], xmm5[0], xmm2[0], xmm3[0], xmm6[0], xmm7[0], srcreg+d1+56); \
	f64ptr(srcreg+d1) = xmm5[0];	/* Save R1 */ \
	f64ptr(srcreg+d1+8) = xmm1[0];	/* Save R2 */ \
	f64ptr(srcreg+d1+16) = xmm7[0];	/* Save R3 */ \
	f64ptr(srcreg+d1+24) = xmm6[0];	/* Save R4 */ \
	f64ptr(srcreg+d1+32) = xmm0[0];	/* Save R5 */ \
	f64ptr(srcreg+d1+40) = xmm4[0];	/* Save R6 */ \
	f64ptr(srcreg+d1+48) = xmm3[0];	/* Save R7 */ \
	f64ptr(srcreg+d1+56) = xmm2[0];	/* Save R8 */ \
	srcreg += srcinc; \


/* Do an eight_reals_with_square_2 on 8 doubles and */
/* a four_complex_with_square_2 on 8 doubles */
#define s2cl_eight_reals_with_square_2(srcreg,srcinc,d1); \
	xmult7	srcreg, srcreg; \
\
	xmm0[0] = f64ptr(srcreg);	/* R1 */ \
	xmm1[0] = f64ptr(srcreg+8);	/* R2 */ \
	xmm2[0] = f64ptr(srcreg+d1);	/* R3 */ \
	xmm3[0] = f64ptr(srcreg+d1+8);	/* R4 */ \
	xmm4[0] = f64ptr(srcreg+16);	/* R5 */ \
	xmm5[0] = f64ptr(srcreg+24);	/* R6 */ \
	xmm6[0] = f64ptr(srcreg+d1+16);	/* R7 */ \
	xmm7[0] = f64ptr(srcreg+d1+24);	/* R8 */ \
	xs8r_fft xmm0, xmm1, xmm2, xmm3, xmm4, xmm5, xmm6, xmm7; \
\
	xmm3[0] *= xmm3[0];		/* R1 = R1 * R1 */ \
	xmm2[0] *= xmm2[0];		/* R2 = R2 * R2 */ \
	f64ptr(rsi-16) = xmm3[0];	/* Save product of sum of FFT values */ \
	xs_complex_square(xmm0[0], xmm1[0], xmm3[0]);	/* Square R3, R4 */ \
	xs_complex_square(xmm5[0], xmm7[0], xmm3[0]);	/* Square R5, R6 */ \
	xs_complex_square(xmm4[0], xmm6[0], xmm3[0]);	/* Square R7, R8 */ \
	xmm3[0] = f64ptr(rsi-16);	/* Restore xmm3 */ \
\
	xs8r_unfft(xmm3[0], xmm2[0], xmm0[0], xmm1[0], xmm5[0], xmm7[0], xmm4[0], xmm6[0]); \
\
	f64ptr(srcreg) = xmm0[0];	/* Save R1 */ \
	f64ptr(srcreg+8) = xmm1[0];	/* Save R2 */ \
	f64ptr(srcreg+16) = xmm2[0];	/* Save R3 */ \
	f64ptr(srcreg+24) = xmm3[0];	/* Save R4 */ \
	xmm0[0] = f64ptr(srcreg+32);	/* R1 */ \
	xmm1[0] = f64ptr(srcreg+40);	/* R2 */ \
	xmm2[0] = f64ptr(srcreg+48);	/* R5 */ \
	xmm3[0] = f64ptr(srcreg+56);	/* R6 */ \
	f64ptr(srcreg+32) = xmm4[0];	/* Save R5 */ \
	f64ptr(srcreg+40) = xmm5[0];	/* Save R6 */ \
	f64ptr(srcreg+48) = xmm6[0];	/* Save R7 */ \
	f64ptr(srcreg+56) = xmm7[0];	/* Save R8 */ \
	xmm4[0] = f64ptr(srcreg+d1+32);	/* R3 */ \
	xmm5[0] = f64ptr(srcreg+d1+40);	/* R4 */ \
	xmm6[0] = f64ptr(srcreg+d1+48);	/* R7 */ \
\
	xs4c_fft(xmm0[0], xmm1[0], xmm4[0], xmm5[0], xmm2[0], xmm3[0], xmm6[0], xmm7[0], srcreg+d1+56); \
\
	g->u.xmm.XMM_TMP1 = xmm0[0]; \
	xs_complex_square xmm5, xmm1, xmm0;	/* Square R1, R2 */ \
	xs_complex_square xmm7, xmm6, xmm0;	/* Square R3, R4 */ \
	xmm0[0] = g->u.xmm.XMM_TMP1; \
	g->u.xmm.XMM_TMP1 = xmm5[0]; \
	xs_complex_square xmm0, xmm4, xmm5;	/* Square R5, R6 */ \
	xs_complex_square xmm3, xmm2, xmm5;	/* Square R7, R8 */ \
	xmm5[0] = g->u.xmm.XMM_TMP1; \
\
	xs4c_unfft(xmm5[0], xmm1[0], xmm7[0], xmm6[0], xmm0[0], xmm4[0], xmm3[0], xmm2[0]); \
\
	f64ptr(srcreg+d1) = xmm0[0];	/* Save R1 */ \
	f64ptr(srcreg+d1+8) = xmm4[0];	/* Save R2 */ \
	f64ptr(srcreg+d1+16) = xmm7[0];	/* Save R3 */ \
	f64ptr(srcreg+d1+24) = xmm5[0];	/* Save R4 */ \
	f64ptr(srcreg+d1+32) = xmm2[0];	/* Save R5 */ \
	f64ptr(srcreg+d1+40) = xmm3[0];	/* Save R6 */ \
	f64ptr(srcreg+d1+48) = xmm6[0];	/* Save R7 */ \
	f64ptr(srcreg+d1+56) = xmm1[0];	/* Save R8 */ \
	srcreg += srcinc; \

\
/* Do an eight_reals_with_mult_2 on 8 doubles and */ \
/* a four_complex_with_mult on 8 doubles */ \
#define s2cl_eight_reals_with_mult_2(srcreg,srcinc,d1); \
	xmult7	srcreg, srcreg+rbp; \
\
	xmm0[0] = f64ptr(srcreg);	/* R1 */ \
	xmm1[0] = f64ptr(srcreg+8);	/* R2 */ \
	xmm2[0] = f64ptr(srcreg+d1);	/* R3 */ \
	xmm3[0] = f64ptr(srcreg+d1+8);	/* R4 */ \
	xmm4[0] = f64ptr(srcreg+16);	/* R5 */ \
	xmm5[0] = f64ptr(srcreg+24);	/* R6 */ \
	xmm6[0] = f64ptr(srcreg+d1+16);	/* R7 */ \
	xmm7[0] = f64ptr(srcreg+d1+24);	/* R8 */ \
	xs8r_fft(xmm0[0], xmm1[0], xmm2[0], xmm3[0], xmm4[0], xmm5[0], xmm6[0], xmm7[0]); \
\
	xs8r_mulf(xmm3[0], xmm2[0], xmm0[0], xmm1[0], xmm5[0], xmm7[0], xmm4[0], xmm6[0], srcreg, srcreg+8, srcreg+16, srcreg+24, srcreg+32, srcreg+40, srcreg+48, srcreg+56); \
\
	xs8r_unfft(xmm3[0], xmm2[0], xmm0[0], xmm1[0], xmm5[0], xmm7[0], xmm4[0], xmm6[0]); \
\
	f64ptr(srcreg) = xmm0[0];	/* Save R1 */ \
	f64ptr(srcreg+8) = xmm1[0];	/* Save R2 */ \
	f64ptr(srcreg+16) = xmm2[0];	/* Save R3 */ \
	f64ptr(srcreg+24) = xmm3[0];	/* Save R4 */ \
	xmm0[0] = f64ptr(srcreg+32);	/* R1 */ \
	xmm1[0] = f64ptr(srcreg+40);	/* R2 */ \
	xmm2[0] = f64ptr(srcreg+48);	/* R5 */ \
	xmm3[0] = f64ptr(srcreg+56);	/* R6 */ \
	f64ptr(srcreg+32) = xmm4[0];	/* Save R5 */ \
	f64ptr(srcreg+40) = xmm5[0];	/* Save R6 */ \
	f64ptr(srcreg+48) = xmm6[0];	/* Save R7 */ \
	f64ptr(srcreg+56) = xmm7[0];	/* Save R8 */ \
	xmm4[0] = f64ptr(srcreg+d1+32);	/* R3 */ \
	xmm5[0] = f64ptr(srcreg+d1+40);	/* R4 */ \
	xmm6[0] = f64ptr(srcreg+d1+48);	/* R7 */ \
\
	xs4c_fft(xmm0[0], xmm1[0], xmm4[0], xmm5[0], xmm2[0], xmm3[0], xmm6[0], xmm7[0], srcreg+d1+56); \
\
	xs4c_mulf(xmm5[0], xmm1[0], xmm7[0], xmm6[0], xmm0[0], xmm4[0], xmm3[0], xmm2[0], srcreg+d1, srcreg+d1+8, srcreg+d1+16, srcreg+d1+24, srcreg+d1+32, srcreg+d1+40, srcreg+d1+48, srcreg+d1+56); \
\
	xs4c_unfft(xmm5[0], xmm1[0], xmm7[0], xmm6[0], xmm0[0], xmm4[0], xmm3[0], xmm2[0]); \
\
	f64ptr(srcreg+d1) = xmm0[0];	/* Save R1 */ \
	f64ptr(srcreg+d1+8) = xmm4[0];	/* Save R2 */ \
	f64ptr(srcreg+d1+16) = xmm7[0];	/* Save R3 */ \
	f64ptr(srcreg+d1+24) = xmm5[0];	/* Save R4 */ \
	f64ptr(srcreg+d1+32) = xmm2[0];	/* Save R5 */ \
	f64ptr(srcreg+d1+40) = xmm3[0];	/* Save R6 */ \
	f64ptr(srcreg+d1+48) = xmm6[0];	/* Save R7 */ \
	f64ptr(srcreg+d1+56) = xmm1[0];	/* Save R8 */ \
	srcreg += srcinc; \


/* Do an eight_reals_with_mulf_2 on 8 doubles and */ \
/* a four_complex_with_mulf_2 on 8 doubles */
#define s2cl_eight_reals_with_mulf_2(srcreg,srcinc,d1); \
	xmult7	srcreg+rbx, srcreg+rbp; \
\
	xmm3[0] = f64ptr(srcreg+rbx);	/* R1 */ \
	xmm2[0] = f64ptr(srcreg+8+rbx);	/* R2 */ \
	xmm0[0] = f64ptr(srcreg+16+rbx);/* R3 */ \
	xmm1[0] = f64ptr(srcreg+24+rbx);/* R4 */ \
	xmm5[0] = f64ptr(srcreg+32+rbx);/* R5 */ \
	xmm7[0] = f64ptr(srcreg+40+rbx);/* R6 */ \
	xmm4[0] = f64ptr(srcreg+48+rbx);/* R7 */ \
	xmm6[0] = f64ptr(srcreg+56+rbx);/* R8 */ \
\
	xs8r_mulf(xmm3[0], xmm2[0], xmm0[0], xmm1[0], xmm5[0], xmm7[0], xmm4[0], xmm6[0], srcreg, srcreg+8, srcreg+16, srcreg+24, srcreg+32, srcreg+40, srcreg+48, srcreg+56); \
\
	xs8r_unfft(xmm3[0], xmm2[0], xmm0[0], xmm1[0], xmm5[0], xmm7[0], xmm4[0], xmm6[0]); \
\
	f64ptr(srcreg) = xmm0[0];	/* Save R1 */ \
	f64ptr(srcreg+8) = xmm1[0];	/* Save R2 */ \
	f64ptr(srcreg+16) = xmm2[0];	/* Save R3 */ \
	f64ptr(srcreg+24) = xmm3[0];	/* Save R4 */ \
	f64ptr(srcreg+32) = xmm4[0];	/* Save R5 */ \
	f64ptr(srcreg+40) = xmm5[0];	/* Save R6 */ \
	f64ptr(srcreg+48) = xmm6[0];	/* Save R7 */ \
	f64ptr(srcreg+56) = xmm7[0];	/* Save R8 */ \
\
	xmm5[0] = f64ptr(srcreg+d1+rbx);/* R1 */ \
	xmm7[0] = f64ptr(srcreg+d1+8+rbx);/* R2 */ \
	xmm4[0] = f64ptr(srcreg+d1+16+rbx);/* R5 */ \
	xmm6[0] = f64ptr(srcreg+d1+24+rbx);/* R6 */ \
	xmm0[0] = f64ptr(srcreg+d1+32+rbx);/* R3 */ \
	xmm1[0] = f64ptr(srcreg+d1+40+rbx);/* R4 */ \
	xmm3[0] = f64ptr(srcreg+d1+48+rbx);/* R7 */ \
	xmm2[0] = f64ptr(srcreg+d1+56+rbx);/* R8 */ \
\
	xs4c_mulf(xmm5[0], xmm7[0], xmm4[0], xmm6[0], xmm0[0], xmm1[0], xmm3[0], xmm2[0], srcreg+d1, srcreg+d1+8, srcreg+d1+16, srcreg+d1+24, srcreg+d1+32, srcreg+d1+40, srcreg+d1+48, srcreg+d1+56); \
\
	xs4c_unfft(xmm5[0], xmm7[0], xmm4[0], xmm6[0], xmm0[0], xmm1[0], xmm3[0], xmm2[0]); \
\
	f64ptr(srcreg+d1) = xmm0[0];	/* Save R1 */ \
	f64ptr(srcreg+d1+8) = xmm1[0];	/* Save R2 */ \
	f64ptr(srcreg+d1+16) = xmm4[0];	/* Save R3 */ \
	f64ptr(srcreg+d1+24) = xmm5[0];	/* Save R4 */ \
	f64ptr(srcreg+d1+32) = xmm2[0];	/* Save R5 */ \
	f64ptr(srcreg+d1+40) = xmm3[0];	/* Save R6 */ \
	f64ptr(srcreg+d1+48) = xmm6[0];	/* Save R7 */ \
	f64ptr(srcreg+d1+56) = xmm7[0];	/* Save R8 */ \
	srcreg += srcinc; \
\


#define xs4c_fft(r1, r2, r3, r4, r5, r6, r7, r8, mem8); \
	r8 = f64ptr(rdi+32+24);	/* cosine/sine */ \
	r8 *= r3;			/* A3 = R3 * cosine/sine	;1-6 */ \
	r8 -= r7;			/* A3 = A3 - I3			;8-11 */ \
	r7 *= f64ptr(rdi+32+24);	/* B3 = I3 * cosine/sine	;3-8 */ \
	r7 += r3;			/* B3 = B3 + R3			;10-13 */ \
	r3 = f64ptr(rdi+0+24);	/* cosine/sine */ \
	r3 *= r2;			/* A2 = R2 * cosine/sine	;5-10 */ \
	r3 -= r6;			/* A2 = A2 - I2			;12-15 */ \
	r6 *= f64ptr(rdi+0+24);	/* B2 = I2 * cosine/sine	;9-14 */ \
	r6 += r2;			/* B2 = B2 + R2			;16-19 */ \
	r2 = f64ptr(rdi+64+24);	/* cosine/sine */ \
	r2 *= mem8;		/* B4 = I4 * cosine/sine	;11-16 */ \
	r2 += r4;			/* B4 = B4 + R4			;18-21 */ \
	r4 *= f64ptr(rdi+64+24);	/* A4 = R4 * cosine/sine	;7-12 */ \
	r4 -= mem8;		/* A4 = A4 - I4			;14-17 */ \
	r8 *= f64ptr(rdi+32+8);	/* A3 = A3 * sine (new R3)	;13-18 */ \
	r7 *= f64ptr(rdi+32+8);	/* B3 = B3 * sine (new I3)	;15-20 */ \
	r3 *= f64ptr(rdi+0+8);		/* A2 = A2 * sine (new R2)	;17-22 */ \
	r4 *= f64ptr(rdi+64+8);	/* A4 = A4 * sine (new R4)	;19-24 */ \
	 r1 -= r8;			/* R1 = R1 - R3 (new R3)	;20-23 */ \
	 r8 *= 2; \
	r6 *= f64ptr(rdi+0+8);		/* B2 = B2 * sine (new I2)	;21-26 */ \
	 r5 -= r7;			/* I1 = I1 - I3 (new I3)	;22-25 */ \
	 r7 *= 2; \
	r2 *= f64ptr(rdi+64+8);	/* B4 = B4 * sine (new I4)	;23-28 */ \
	 r8 += r1;			/* R3 = R1 + R3 (new R1)	;24-27 */ \
	 r3 -= r4;			/* R2 = R2 - R4 (new R4)	;26-29 */ \
	 r4 *= 2;			/* R4 = R4 * 2			;27-32 */ \
	 r7 += r5;			/* I3 = I1 + I3 (new I1)	;28-31 */ \
	 r6 -= r2;			/* I2 = I2 - I4 (new I4)	;30-33 */ \
	 r2 *= 2;			/* I4 = I4 * 2			;31-36 */ \
	r5 -= r3;			/* I3 = I3 - R4 (final I4)	;32-35 */ \
	 r4 += r3;			/* R4 = R2 + R4 (new R2)	;34-37 */ \
	r3 *= 2;			/* R4 = R4 * 2			;35-40 */ \
	 r2 += r6;			/* I4 = I2 + I4 (new I2)	;36-39 */ \
	r1 -= r6;			/* R3 = R3 - I4 (final R3)	;38-41 */ \
	r6 *= 2;			/* I4 = I4 * 2			;39-44 */ \
	r8 -= r4;			/* R1 = R1 - R2 (final R2)	;40-43 */ \
	r4 *= 2;			/* R2 = R2 * 2			;41-46 */ \
	r7 -= r2;			/* I1 = I1 - I2 (final I2)	;42-45 */ \
	r2 *= 2;			/* I2 = I2 * 2			;43-48 */ \
	r3 += r5;			/* R4 = I3 + R4 (final I3)	;44-47 */ \
	r6 += r1;			/* I4 = R3 + I4 (final R4)	;46-49 */ \
	r4 += r8;			/* R2 = R1 + R2 (final R1)	;48-51 */ \
	r2 += r7;			/* I2 = I1 + I2 (final I1)	;50-53 */ \

#define s2cl_four_complex_fft_final(srcreg,srcinc,d1) \
	shuffle_load(xmm0,xmm1,u64ptr(srcreg),u64ptr(srcreg+32)); /* R1,R2 */ \
	shuffle_load(xmm2,xmm3,u64ptr(srcreg+d1),u64ptr(srcreg+d1+32)); /* R3,R4 */ \
	shuffle_load(xmm4,xmm5,u64ptr(srcreg+16),u64ptr(srcreg+48)); /* R5,R6 */ \
	shuffle_load(xmm6,xmm7,u64ptr(srcreg+d1+16),u64ptr(srcreg+d1+48)); /* R7,R8 */ \
	xptr(g->u.xmm.XMM_TMP1) = xmm7; \
	x4c_fft(xmm0, xmm1, xmm2, xmm3, xmm4, xmm5, xmm6, xmm7, xptr(g->u.xmm.XMM_TMP1), rdi, 0, srcreg+srcinc, d1, srcreg+d1+48, srcreg+d1); \
	xptr(srcreg) = xmm3;		/* Save R1 */ \
	xptr(srcreg+16) = xmm1;	/* Save R2 */ \
	xptr(srcreg+32) = xmm7;	/* Save R3 */ \
	xptr(srcreg+48) = xmm6;	/* Save R4 */ \
/*	xstore	[srcreg+d1], xmm0	;; Save R5 */ \
	xptr(srcreg+d1+16) = xmm2;	/* Save R6 */ \
	xptr(srcreg+d1+32) = xmm5;	/* Save R7 */ \
/*	xstore	[srcreg+d1+48], xmm4	;; Save R8 */ \
	srcreg += srcinc; \
