#pragma once

/* Copyright 2009-2010 - Mersenne Research, Inc.  All rights reserved */ \
/* Author:  George Woltman */ \
/* Email: woltman@alum.mit.edu */ \
/* */ \
\
/* */ \
/* */ \
/* All new macros for version 26 of gwnum where we do a very traditional, primarily */ \
/* radix-4, FFT.  These macros do a radix-2 step for FFTs with an odd number of levels. */ \
/* The forward FFT macros multiply by the sin/cos values at the end */ \
/* of the macro and the inverse FFTs multiply by the sin/cos values at the start of */ \
/* the macro.  We also implement the Daniel J. Bernstein (DJB) "exponent-1" idea to */ \
/* save sin/cos memory. */ \
/* */ \
/* */ \
\
/* Macro to operate on 2 64-byte cache lines.  It does 2 two_complex_ffts. */ \
/* That is, two radix-2 butterflies followed by a complex multiply using sin/cos data. */ \
\

#define r2_x2cl_two_complex_fft(srcreg,srcinc,d1,screg) \
	xmm2 = xptr(srcreg);		/* R1 */ \
	xmm0 = xptr(srcreg+d1);	/* R2 */ \
	xmm1 = xmm2;		/* Copy R1 */ \
	xmm2 -= xmm0;		/* R1 = R1 - R2 (new R2) */ \
	xmm0 += xmm1;		/* R2 = R1 + R2 (new and final R1) */ \
\
	xmm3 = xptr(srcreg+16);	/* R3 (I1) */ \
	xmm1 = xptr(srcreg+d1+16);	/* R4 (I2) */ \
	xmm4 = xmm3;		/* Copy I1 */ \
	xmm3 -= xmm1;		/* I1 = I1 - I2 (new I2) */ \
	xmm1 += xmm4;		/* I2 = I1 + I2 (new and final I1) */ \
\
	xmm6 = xptr(srcreg+32);	/* R1_2 */ \
	xmm4 = xptr(srcreg+d1+32);	/* R2_2 */ \
	xmm7 = xmm6;		/* Copy R1_2 */ \
	xmm6 -= xmm4;		/* R1_2 = R1_2-R2_2 (new R2_2) */ \
	xmm4 += xmm7;		/* R2_2 = R1_2+R2_2 (new and final R1_2) */ \
\
	xptr(srcreg) = xmm0;		/* Save R1 */ \
	xptr(srcreg+16) = xmm1;	/* Save R2 */ \
\
	xprefetchw(u8ptr(srcreg+srcinc)); \
\
	xmm7 = xptr(srcreg+48);	/* R3_2 */ \
	xmm5 = xptr(srcreg+d1+48);	/* R4_2 */ \
	xmm0 = xmm7;		/* Copy I1_2 */ \
	xmm7 -= xmm5;		/* I1_2 = I1_2-I2_2 (new I2_2) */ \
	xmm5 += xmm0;		/* I2_2 = I1_2+I2_2 (new and final I1_2) */ \
\
	xmm0 = xptr(screg+16);	/* cosine/sine */ \
	xmm1 = xmm0;		/* cosine/sine */ \
	xmm0 *= xmm2;		/* A2 = R2 * cosine/sine */ \
	xmm1 *= xmm3;		/* B2 = I2 * cosine/sine */ \
\
	xptr(srcreg+d1) = xmm4;	/* Save R1_2 */ \
	xptr(srcreg+d1+16) = xmm5;	/* Save R2_2 */ \
\
	xmm4 = xptr(screg+16);	/* cosine/sine */ \
	xmm5 = xmm4;		/* cosine/sine */ \
	xmm4 *= xmm6;		/* A2_2 = R2_2 * cosine/sine */ \
	xmm5 *= xmm7;		/* B2_2 = I2_2 * cosine/sine */ \
\
	xprefetchw(u8ptr(srcreg+srcinc+d1)); \
\
	xmm0 -= xmm3;		/* A2 = A2 - I2 */ \
	xmm1 += xmm2;		/* B2 = B2 + R2 */ \
\
	xmm4 -= xmm7;		/* A2_2 = A2_2 - I2_2 */ \
	xmm5 += xmm6;		/* B2_2 = B2_2 + R2_2 */ \
\
	xmm3 = xptr(screg);		/* Load sine */ \
	xmm0 *= xmm3;		/* A2 = A2 * sine (final R2) */ \
	xmm1 *= xmm3;		/* B2 = B2 * sine (final I2) */ \
\
	xmm4 *= xmm3;		/* A2_2 = A2_2 * sine (final R2_2) */ \
	xmm5 *= xmm3;		/* B2_2 = B2_2 * sine (final I2_2) */ \
\
	xptr(srcreg+32) = xmm0;	/* Save R3 */ \
	xptr(srcreg+48) = xmm1;	/* Save R4 */ \
	xptr(srcreg+d1+32) = xmm4;	/* Save R3_2 */ \
	xptr(srcreg+d1+48) = xmm5;	/* Save R4_2 */ \
\
	srcreg += srcinc; \
\
\
\
/* Macro to operate on 2 64-byte cache lines.  It does 2 two_complex_unffts */ \
/* That is, two complex multiplies using sin/cos data followed by radix-2 butterflies. */ \
\

#define r2_x2cl_two_complex_unfft(srcreg,srcinc,d1,screg,scdist) \
	xmm2 = xptr(srcreg+d1);		/* R3 (R2) */ \
	xmm3 = xptr(srcreg+d1+32);		/* R4 (I2) */ \
\
	xmm4 = xptr(screg+16);		/* cosine/sine */ \
	xmm4 *= xmm2;			/* A2 = R2 * cosine/sine */ \
	xmm4 += xmm3;			/* A2 = A2 + I2 */ \
	xmm5 = xptr(screg+16);		/* cosine/sine */ \
	xmm5 *= xmm3;			/* B2 = I2 * cosine/sine */ \
	xmm5 -= xmm2;			/* B2 = B2 - R2 */ \
	xmm3 = xptr(screg);			/* sine */ \
	xmm4 *= xmm3;			/* A2 = A2 * sine (new R2) */ \
	xmm5 *= xmm3;			/* B2 = B2 * sine (new I2) */ \
\
	xprefetchw(u8ptr(srcreg+srcinc)); \
\
	xmm6 = xptr(srcreg+d1+16);		/* R3 (R2_2) */ \
	xmm7 = xptr(srcreg+d1+48);		/* R4 (I2_2) */ \
	xmm2 = xptr(screg+scdist+16);		/* cosine/sine */ \
	xmm2 *= xmm6;			/* A2_2 = R2_2 * cosine/sine */ \
	xmm2 += xmm7;			/* A2_2 = A2_2 + I2_2 */ \
	xmm3 = xptr(screg+scdist+16);		/* cosine/sine */ \
	xmm3 *= xmm7;			/* B2_2 = I2_2 * cosine/sine */ \
	xmm3 -= xmm6;			/* B2_2 = B2_2 - R2 */ \
	xmm7 = xptr(screg+scdist);		/* sine */ \
	xmm2 *= xmm7;			/* A2_2 = A2_2 * sine (new R2_2) */ \
	xmm3 *= xmm7;			/* B2_2 = B2_2 * sine (new I2_2) */ \
\
	xprefetchw(u8ptr(srcreg+srcinc+d1)); \
\
	xmm0 = xptr(srcreg);			/* R1 (R1) */ \
	xmm0 -= xmm4;			/* final R2 = R1 - R2 */ \
	xmm4 += xptr(srcreg);			/* final R1 = R1 + R2 */ \
\
	xmm1 = xptr(srcreg+32);		/* R2 (I1) */ \
	xmm1 -= xmm5;			/* final I2 = I1 - I2 */ \
	xmm5 += xptr(srcreg+32);		/* final I1 = I1 + I2 */ \
\
	xmm6 = xptr(srcreg+16);		/* R1 (R1_2) */ \
	xmm6 -= xmm2;			/* final R2_2 = R1_2 - R2_2 */ \
	xmm2 += xptr(srcreg+16);		/* final R1_2 = R1_2 + R2_2 */ \
\
	xmm7 = xptr(srcreg+48);		/* R2 (I1_2) */ \
	xmm7 -= xmm3;			/* final I2_2 = I1_2 - I2_2 */ \
	xmm3 += xptr(srcreg+48);		/* final I1_2 = I1_2 + I2_2 */ \
\
	xptr(srcreg) = xmm4;			/* Save R1 */ \
	xptr(srcreg+16) = xmm0;		/* Save R2 */ \
	xptr(srcreg+32) = xmm5;		/* Save R3 */ \
	xptr(srcreg+48) = xmm1;		/* Save R4 */ \
\
	xptr(srcreg+d1) = xmm2;		/* Save R1_2 */ \
	xptr(srcreg+d1+16) = xmm6;		/* Save R2_2 */ \
	xptr(srcreg+d1+32) = xmm3;		/* Save R3_2 */ \
	xptr(srcreg+d1+48) = xmm7;		/* Save R4_2 */ \
	srcreg += srcinc; \
\
