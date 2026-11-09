#pragma once
	//

#define s2cl_eight_reals_first_fft(srcreg,srcinc,d1) { \
	uintptr_t src_rbx = srcreg+rbx; \
	shuffle_load(xmm0,xmm2,xptr(src_rbx),xptr(src_rbx+16)); /* R1,R3 */ \
	shuffle_load(xmm1,xmm3,xptr(src_rbx+d1),xptr(src_rbx+d1+16)); /* R2,R4 */ \
	shuffle_load(xmm4,xmm6,xptr(src_rbx+32),xptr(src_rbx+48)); /* R5,R7 */ \
	shuffle_load(xmm5,xmm7,xptr(src_rbx+d1+32),xptr(src_rbx+d1+48)); /* R6,R8 */ \
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
}

#define x2cl_eight_reals_first_fft(srcreg,srcinc,d1) { \
	uintptr_t src_rbx = srcreg+rbx; \
	xmm0 = xptr(src_rbx); \
	xmm1 = xptr(src_rbx+d1); \
	xmm2 = xptr(src_rbx+16); \
	xmm3 = xptr(src_rbx+d1+16); \
	xmm4 = xptr(src_rbx+32); \
	xmm5 = xptr(src_rbx+d1+32); \
	xmm6 = xptr(src_rbx+48); \
	xmm7 = xptr(src_rbx+d1+48); \
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
}

#define x2cl_eight_reals_fft(srcreg,srcinc,d1) \
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
\

#define g2cl_eight_reals_first_fft(srcreg,srcinc,d1,dstreg,dstinc,e1) { \
	uintptr_t src_rbx = srcreg+rbx; \
	xmm0 = xptr(src_rbx); \
	xmm1 = xptr(src_rbx+d1); \
	xmm2 = xptr(src_rbx+16); \
	xmm3 = xptr(src_rbx+d1+16); \
	xmm4 = xptr(src_rbx+32); \
	xmm5 = xptr(src_rbx+d1+32); \
	xmm6 = xptr(src_rbx+48); \
	xmm7 = xptr(src_rbx+d1+48); \
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
}

#if 1
#define x8r_fft { \
	vec2f64 temp2; \
\
	temp2 = xmm1; \
	xmm1 -= xmm5;		/* new R6 = R2 - R6 */ \
	xmm5 += temp2;		/* new R2 = R2 + R6 */ \
\
	temp2 = xmm3; \
	xmm3 -= xmm7;		/* new R8 = R4 - R8 */ \
	xmm7 += temp2;		/* new R4 = R4 + R8 */ \
\
	temp2 = xmm0; \
	xmm0 -= xmm4;		/* new R5 = R1 - R5 */ \
	xmm4 += temp2;		/* new R1 = R1 + R5 */ \
\
	temp2 = xmm2; \
	xmm2 -= xmm6;		/* new R7 = R3 - R7 */ \
	xmm6 += temp2;		/* new R3 = R3 + R7 */ \
\
	xmm1 *= XMM_SQRTHALF2;	/* R6 = R6 * square root */ \
	xmm3 *= XMM_SQRTHALF2;	/* R8 = R8 * square root */ \
\
	temp2 = xmm5; \
	xmm5 -= xmm7;		/* R2 = R2 - R4 (new & final R4) */ \
	xmm7 += temp2;		/* R4 = R2 + R4 (new R2) */ \
\
	temp2 = xmm4; \
	xmm4 -= xmm6;		/* R1 = R1 - R3 (new & final R3) */ \
	xmm6 += temp2;		/* R3 = R1 + R3 (new R1) */ \
\
	temp2 = xmm1; \
	xmm1 -= xmm3;		/* R6 = R6 - R8 (Real part) */ \
	xmm3 += temp2;		/* R8 = R6 + R8 (Imaginary part) */ \
\
	temp2 = xmm0; \
	xmm0 -= xmm1;		/* R5 = R5 - R6 (final R7) */ \
	xmm1 += temp2;		/* R6 = R5 + R6 (final R5) */ \
\
	temp2 = xmm2; \
	xmm2 -= xmm3;		/* R7 = R7 - R8 (final R8) */ \
	xmm3 += temp2;		/* R8 = R7 + R8 (final R6) */ \
\
	temp2 = xmm6; \
	xmm6 -= xmm7;		/* R1 = R1 - R2 (final R2) */ \
	xmm7 += temp2;		/* R2 = R1 + R2 (final R1) */ \
}
#else
#define x8r_fft \
	xmm3 -= xmm7;		/* new R8 = R4 - R8 */ \
	xmm7 += xmm7; /* multwo */ \
	xmm7 += xmm3;		/* new R4 = R4 + R8 */ \
	xmm1 -= xmm5;		/* new R6 = R2 - R6 */ \
	xmm5 += xmm5; /* multwo */ \
	xmm5 += xmm1;		/* new R2 = R2 + R6 */ \
	 xmm3 *= XMM_SQRTHALF2;	/* R8 = R8 * square root */ \
	 xmm1 *= XMM_SQRTHALF2;	/* R6 = R6 * square root */ \
	xmm0 -= xmm4;		/* new R5 = R1 - R5 */ \
	xmm4 += xmm4; /* multwo */ \
	xmm4 += xmm0;		/* new R1 = R1 + R5 */ \
	 xmm5 -= xmm7;		/* R2 = R2 - R4 (new & final R4) */ \
	 xmm7 += xmm7; /* multwo */			/* R4 = R4 * 2 */ \
	xmm2 -= xmm6;		/* new R7 = R3 - R7 */ \
	xmm6 += xmm6; /* multwo */ \
	xmm6 += xmm2;		/* new R3 = R3 + R7 */ \
	 xmm1 -= xmm3;		/* R6 = R6 - R8 (Real part) */ \
	 xmm3 += xmm3; /* multwo */			/* R8 = R8 * 2 */ \
	 xmm4 -= xmm6;		/* R1 = R1 - R3 (new & final R3) */ \
	 xmm6 += xmm6; /* multwo */			/* R3 = R3 * 2 */ \
	 xmm7 += xmm5;		/* R4 = R2 + R4 (new R2) */ \
	xmm3 += xmm1;		/* R8 = R6 + R8 (Imaginary part) */ \
	xmm0 -= xmm1;		/* R5 = R5 - R6 (final R7) */ \
	xmm1 += xmm1; /* multwo */			/* R6 = R6 * 2 */ \
	 xmm6 += xmm4;		/* R3 = R1 + R3 (new R1) */ \
	xmm2 -= xmm3;		/* R7 = R7 - R8 (final R8) */ \
	xmm3 += xmm3; /* multwo */			/* R8 = R8 * 2 */ \
	xmm6 -= xmm7;		/* R1 = R1 - R2 (final R2) */ \
	xmm7 += xmm7; /* multwo */			/* R2 = R2 * 2 */ \
	xmm1 += xmm0;		/* R6 = R5 + R6 (final R5) */ \
	xmm3 += xmm2;		/* R8 = R7 + R8 (final R6) */ \
	xmm7 += xmm6;		/* R2 = R1 + R2 (final R1) */
#endif

/* Macro to operate on 4 64-byte cache lines.  It does the last */ \
/* three inverse FFT levels of a one pass FFT. */ \

#define x4cl_eight_reals_last_unfft(srcreg,srcinc,d1,d2) \
	xmm0 = xptr(srcreg);		/* R1 */ \
	xmm1 = xptr(srcreg+32);	/* R2 */ \
	xmm2 = xptr(srcreg+d1);	/* R3 */ \
	xmm3 = xptr(srcreg+d1+32);	/* R4 */ \
	xmm4 = xptr(srcreg+d2);	/* R5 */ \
	xmm5 = xptr(srcreg+d2+32);	/* R6 */ \
	xmm6 = xptr(srcreg+d2+d1);	/* R7 */ \
	xmm7 = xptr(srcreg+d2+d1+32);	/* R8 */ \
	x8r_unfft(xmm0, xmm1, xmm2, xmm3, xmm4, xmm5, xmm6, xmm7); \
	xptr(srcreg) = xmm6;		/* Save R1 */ \
	xptr(srcreg+d2) = xmm4;	/* Save R2 */ \
	xptr(srcreg+32) = xmm2;	/* Save R5 */ \
	xptr(srcreg+d2+32) = xmm3;	/* Save R6 */ \
	xmm6 = xptr(srcreg+16);	/* R1 */ \
	xmm4 = xptr(srcreg+48);	/* R2 */ \
	xmm2 = xptr(srcreg+d2+16);	/* R5 */ \
	xmm3 = xptr(srcreg+d2+48);	/* R6 */ \
	xptr(srcreg+16) = xmm7;	/* Save R3 */ \
	xptr(srcreg+d2+16) = xmm5;	/* Save R4 */ \
	xptr(srcreg+48) = xmm1;	/* Save R7 */ \
	xptr(srcreg+d2+48) = xmm0;	/* Save R8 */ \
	xmm7 = xptr(srcreg+d1+16);	/* R3 */ \
	xmm5 = xptr(srcreg+d1+48);	/* R4 */ \
	xmm1 = xptr(srcreg+d2+d1+16);	/* R7 */ \
	xmm0 = xptr(srcreg+d2+d1+48);	/* R8 */ \
	x8r_unfft(xmm6, xmm4, xmm7, xmm5, xmm2, xmm3, xmm1, xmm0); \
	xptr(srcreg+d1) = xmm1;	/* Save R1 */ \
	xptr(srcreg+d2+d1) = xmm2;	/* Save R2 */ \
	xptr(srcreg+d1+16) = xmm0;	/* Save R3 */ \
	xptr(srcreg+d2+d1+16) = xmm3;	/* Save R4 */ \
	xptr(srcreg+d1+32) = xmm7;	/* Save R5 */ \
	xptr(srcreg+d2+d1+32) = xmm5;	/* Save R6 */ \
	xptr(srcreg+d1+48) = xmm4;	/* Save R7 */ \
	xptr(srcreg+d2+d1+48) = xmm6;	/* Save R8 */ \
	srcreg += srcinc; \
\

#if 1
#define x8r_unfft(r1, r2, r3, r4, r5, r6, r7, r8) { \
	vec2f64 temp; \
\
	temp = r5; \
	r5 -= r7;			/* new R6 = R5 - R7		;5-8 */ \
	r7 += temp;			/* new R5 = R5 + R7		;7-10 */ \
\
	temp = r6; \
	r6 -= r8;			/* new R8 = R6 - R8		;1-4 */ \
	r8 += temp;			/* new R7 = R6 + R8		;3-6 */ \
\
	temp = r1; \
	r1 -= r2;			/* new R2 = R1 - R2		;9-12 */ \
	r2 += temp;			/* new R1 = R1 + R2		;11-14 */ \
\
	temp = r6; \
	r6 -= r5;			/* R8 = R8 - R6			;13-16 */ \
	r5 += temp;			/* R6 = R6 + R8			;15-18 */ \
\
	temp = r1; \
	r1 -= r4;			/* R2 = R2 - R4 (new R4)	;17-20 */ \
	r4 += temp;			/* R4 = R2 + R4 (new R2)	;27-30 */ \
\
	temp = r2; \
	r2 -= r3;			/* R1 = R1 - R3 (new R3)	;19-22 */ \
	r3 += temp;			/* R3 = R1 + R3 (new R1)	; */ \
\
	r5 *= XMM_SQRTHALF2;	/* R6 = R6 * square root of 1/2	;22-27 */ \
	r6 *= XMM_SQRTHALF2;	/* R8 = R8 * square root of 1/2	;18-23 */ \
\
	temp = r2; \
	r2 -= r8;			/* R3 = R3 - R7 (final R7)	; */ \
	r8 += temp;			/* R7 = R3 + R7 (final R3)	; */ \
\
	temp = r4; \
	r4 -= r5;			/* R2 = R2 - R6 (final R6)	; */ \
	r5 += temp;			/* R6 = R2 + R6 (final R2)	; */ \
\
	temp = r1; \
	r1 -= r6;			/* newR4 = newR4-newR8(final R8); */ \
	r6 += temp;			/* R8 = R4 + R8 (final R4)	; */ \
\
	temp = r3; \
	r3 -= r7;			/* R1 = R1 - R5 (final R5)	; */ \
	r7 += temp;			/* R5 = R1 + R5 (final R1)	; */ \
}
#else
#define x8r_unfft(r1, r2, r3, r4, r5, r6, r7, r8) \
	r6 -= r8;			/* new R8 = R6 - R8		;1-4 */ \
	r8 += r8; /* multwo */ \
	r8 += r6;			/* new R7 = R6 + R8		;3-6 */ \
	r5 -= r7;			/* new R6 = R5 - R7		;5-8 */ \
	r7 += r7; /* multwo */ \
	r7 += r5;			/* new R5 = R5 + R7		;7-10 */ \
	r1 -= r2;			/* new R2 = R1 - R2		;9-12 */ \
	r2 += r2; /* multwo */ \
	r2 += r1;			/* new R1 = R1 + R2		;11-14 */ \
	r6 -= r5;			/* R8 = R8 - R6			;13-16 */ \
	r5 += r5; /* multwo */ \
	r5 += r6;			/* R6 = R6 + R8			;15-18 */ \
	r1 -= r4;			/* R2 = R2 - R4 (new R4)	;17-20 */ \
	r6 *= XMM_SQRTHALF2;	/* R8 = R8 * square root of 1/2	;18-23 */ \
	r4 += r4; /* multwo */			/* R4 = R4 * 2			;20-25 */ \
	r5 *= XMM_SQRTHALF2;	/* R6 = R6 * square root of 1/2	;22-27 */ \
	r2 -= r3;			/* R1 = R1 - R3 (new R3)	;19-22 */ \
	r3 += r3; /* multwo */			/* R3 = R3 * 2			;24-29 */ \
	r4 += r1;			/* R4 = R2 + R4 (new R2)	;27-30 */ \
	r1 -= r6;			/* newR4 = newR4-newR8(final R8); */ \
	r6 += r6; /* multwo */			/* R8 = R8 * 2			; */ \
	r3 += r2;			/* R3 = R1 + R3 (new R1)	; */ \
	r2 -= r8;			/* R3 = R3 - R7 (final R7)	; */ \
	r8 += r8; /* multwo */			/* R7 = R7 * 2			; */ \
	r4 -= r5;			/* R2 = R2 - R6 (final R6)	; */ \
	r5 += r5; /* multwo */			/* R6 = R6 * 2			; */ \
	r3 -= r7;			/* R1 = R1 - R5 (final R5)	; */ \
	r7 += r7; /* multwo */			/* R5 = R5 * 2			; */ \
	r6 += r1;			/* R8 = R4 + R8 (final R4)	; */ \
	r8 += r2;			/* R7 = R3 + R7 (final R3)	; */ \
	r5 += r4;			/* R6 = R2 + R6 (final R2)	; */ \
	r7 += r3;			/* R5 = R1 + R5 (final R1)	; */
#endif


#define xfive_reals_fft_preload \
\

#define s5cl_five_reals_first_fft(srcreg,srcinc,d1) { \
	uintptr_t src_rbx = srcreg+rbx; \
	shuffle_load(xmm0,xmm7,xptr(src_rbx),xptr(src_rbx+16)); /* R1,R1 */ \
	xptr(srcreg+16) =xmm7;		/* Save it */ \
	shuffle_load(xmm1,xmm2,xptr(src_rbx+3*d1),xptr(src_rbx+3*d1+16)); /*R2,R3 */ \
	shuffle_load(xmm3,xmm4,xptr(src_rbx+d1+32),xptr(src_rbx+d1+48)); /* R4,R5 */ \
	x5r_fft(xmm0, xmm1, xmm2, xmm3, xmm4, xmm5, xmm6, xmm7); \
	xptr(srcreg) = xmm0;			/* Save R1 */ \
	shuffle_load(xmm0,xmm2,xptr(src_rbx+d1),xptr(src_rbx+d1+16)); /* R1,R2 */ \
	xptr(srcreg+d1) = xmm7;		/* Save R2 */ \
	xptr(srcreg+d1+16) = xmm1;		/* Save R3 */ \
	xptr(srcreg+d1+32) = xmm5;		/* Save R4 */ \
	xptr(srcreg+d1+48) = xmm6;		/* Save R5 */ \
	shuffle_load(xmm1,xmm7,xptr(src_rbx+32),xptr(src_rbx+48)); /* R3,R3 */ \
	xptr(srcreg+48) =xmm7;		/* Save it */ \
	shuffle_load(xmm3,xmm4,xptr(src_rbx+3*d1+32),xptr(src_rbx+3*d1+48));/*R4,R5 */ \
	x5r_fft(xmm0, xmm2, xmm1, xmm3, xmm4, xmm5, xmm6, xmm7); \
	xptr(srcreg+32) = xmm0;		/* Save R1 */ \
	xptr(srcreg+3*d1) = xmm7;		/* Save R2 */ \
	xptr(srcreg+3*d1+16) = xmm2;		/* Save R3 */ \
	xptr(srcreg+3*d1+32) = xmm5;		/* Save R4 */ \
	xptr(srcreg+3*d1+48) = xmm6;		/* Save R5 */ \
\
	xmm0 = xptr(srcreg+16);		/* R1 */ \
	shuffle_load(xmm1,xmm2,xptr(src_rbx+4*d1),xptr(src_rbx+4*d1+16));/*R2,R3 */ \
	shuffle_load(xmm3,xmm4,xptr(src_rbx+2*d1+32),xptr(src_rbx+2*d1+48));/*R4,R5 */ \
	x5r_fft(xmm0, xmm1, xmm2, xmm3, xmm4, xmm5, xmm6, xmm7); \
	xptr(srcreg+16) = xmm0;		/* Save R1 */ \
	shuffle_load(xmm0,xmm2,xptr(src_rbx+2*d1),xptr(src_rbx+2*d1+16));/*R1,R2 */ \
	xptr(srcreg+2*d1) = xmm7;		/* Save R2 */ \
	xptr(srcreg+2*d1+16) = xmm1;		/* Save R3 */ \
	xptr(srcreg+2*d1+32) = xmm5;		/* Save R4 */ \
	xptr(srcreg+2*d1+48) = xmm6;		/* Save R5 */ \
	xmm1 = xptr(srcreg+48);		/* R3 */ \
	shuffle_load(xmm3,xmm4,xptr(src_rbx+4*d1+32),xptr(src_rbx+4*d1+48));/*R4,R5 */ \
	x5r_fft(xmm0, xmm2, xmm1, xmm3, xmm4, xmm5, xmm6, xmm7); \
	xptr(srcreg+48) = xmm0;		/* Save R1 */ \
	xptr(srcreg+4*d1) = xmm7;		/* Save R2 */ \
	xptr(srcreg+4*d1+16) = xmm2;		/* Save R3 */ \
	xptr(srcreg+4*d1+32) = xmm5;		/* Save R4 */ \
	xptr(srcreg+4*d1+48) = xmm6;		/* Save R5 */ \
	srcreg += srcinc; \
}

#define x5cl_five_reals_first_fft(srcreg,srcinc,d1) { \
	uintptr_t src_rbx = srcreg+rbx; \
	xmm0 = xptr(src_rbx); \
	xmm1 = xptr(src_rbx+3*d1); \
	xmm2 = xptr(src_rbx+3*d1+16); \
	xmm3 = xptr(src_rbx+d1+32); \
	xmm4 = xptr(src_rbx+d1+48); \
	x5r_fft(xmm0, xmm1, xmm2, xmm3, xmm4, xmm5, xmm6, xmm7); \
	xmm2 = xptr(src_rbx+d1);	/* Load R1 */ \
	xmm3 = xptr(src_rbx+d1+16); /* Load R2 */ \
	xptr(srcreg) = xmm0;		/* Save R1 */ \
	xptr(srcreg+d1) = xmm7;	/* Save R2 */ \
	xptr(srcreg+d1+16) = xmm1;	/* Save R3 */ \
	xptr(srcreg+d1+32) = xmm5;	/* Save R4 */ \
	xptr(srcreg+d1+48) = xmm6;	/* Save R5 */ \
	xmm0 = xptr(src_rbx+32);	/* Load R3 */ \
	xmm1 = xptr(src_rbx+3*d1+32); /* Load R4 */ \
	xmm4 = xptr(src_rbx+3*d1+48); /* Load R5 */ \
	x5r_fft(xmm2, xmm3, xmm0, xmm1, xmm4, xmm5, xmm6, xmm7); \
	xptr(srcreg+32) = xmm2;	/* Save R1 */ \
	xptr(srcreg+3*d1) = xmm7;	/* Save R2 */ \
	xptr(srcreg+3*d1+16) = xmm3;	/* Save R3 */ \
	xptr(srcreg+3*d1+32) = xmm5;	/* Save R4 */ \
	xptr(srcreg+3*d1+48) = xmm6;	/* Save R5 */ \
\
	xmm0 = xptr(src_rbx+16); \
	xmm1 = xptr(src_rbx+4*d1); \
	xmm2 = xptr(src_rbx+4*d1+16); \
	xmm3 = xptr(src_rbx+2*d1+32); \
	xmm4 = xptr(src_rbx+2*d1+48); \
	x5r_fft(xmm0, xmm1, xmm2, xmm3, xmm4, xmm5, xmm6, xmm7); \
	xmm2 = xptr(src_rbx+2*d1); /* Load R1 */ \
	xmm3 = xptr(src_rbx+2*d1+16); /* Load R2 */ \
	xptr(srcreg+16) = xmm0;	/* Save R1 */ \
	xptr(srcreg+2*d1) = xmm7;	/* Save R2 */ \
	xptr(srcreg+2*d1+16) = xmm1;	/* Save R3 */ \
	xptr(srcreg+2*d1+32) = xmm5;	/* Save R4 */ \
	xptr(srcreg+2*d1+48) = xmm6;	/* Save R5 */ \
	xmm0 = xptr(src_rbx+48);	/* Load R3 */ \
	xmm1 = xptr(src_rbx+4*d1+32); /* Load R4 */ \
	xmm4 = xptr(src_rbx+4*d1+48); /* Load R5 */ \
	x5r_fft(xmm2, xmm3, xmm0, xmm1, xmm4, xmm5, xmm6, xmm7); \
	xptr(srcreg+48) = xmm2;	/* Save R1 */ \
	xptr(srcreg+4*d1) = xmm7;	/* Save R2 */ \
	xptr(srcreg+4*d1+16) = xmm3;	/* Save R3 */ \
	xptr(srcreg+4*d1+32) = xmm5;	/* Save R4 */ \
	xptr(srcreg+4*d1+48) = xmm6;	/* Save R5 */ \
	srcreg += srcinc; \
}

#define x5cl_five_reals_fft(srcreg,srcinc,d1) \
	xmm0 = xptr(srcreg); \
	xmm1 = xptr(srcreg+3*d1); \
	xmm2 = xptr(srcreg+3*d1+16); \
	xmm3 = xptr(srcreg+d1+32); \
	xmm4 = xptr(srcreg+d1+48); \
	x5r_fft(xmm0, xmm1, xmm2, xmm3, xmm4, xmm5, xmm6, xmm7); \
	xmm2 = xptr(srcreg+d1);	/* Load R1 */ \
	xmm3 = xptr(srcreg+d1+16);	/* Load R2 */ \
	xptr(srcreg) = xmm0;		/* Save R1 */ \
	xptr(srcreg+d1) = xmm7;	/* Save R2 */ \
	xptr(srcreg+d1+16) = xmm1;	/* Save R3 */ \
	xptr(srcreg+d1+32) = xmm5;	/* Save R4 */ \
	xptr(srcreg+d1+48) = xmm6;	/* Save R5 */ \
	xmm0 = xptr(srcreg+32);	/* Load R3 */ \
	xmm1 = xptr(srcreg+3*d1+32);	/* Load R4 */ \
	xmm4 = xptr(srcreg+3*d1+48);	/* Load R5 */ \
	x5r_fft(xmm2, xmm3, xmm0, xmm1, xmm4, xmm5, xmm6, xmm7); \
	xptr(srcreg+32) = xmm2;	/* Save R1 */ \
	xptr(srcreg+3*d1) = xmm7;	/* Save R2 */ \
	xptr(srcreg+3*d1+16) = xmm3;	/* Save R3 */ \
	xptr(srcreg+3*d1+32) = xmm5;	/* Save R4 */ \
	xptr(srcreg+3*d1+48) = xmm6;	/* Save R5 */ \
\
	xmm0 = xptr(srcreg+16); \
	xmm1 = xptr(srcreg+4*d1); \
	xmm2 = xptr(srcreg+4*d1+16); \
	xmm3 = xptr(srcreg+2*d1+32); \
	xmm4 = xptr(srcreg+2*d1+48); \
	x5r_fft(xmm0, xmm1, xmm2, xmm3, xmm4, xmm5, xmm6, xmm7); \
	xmm2 = xptr(srcreg+2*d1);	/* Load R1 */ \
	xmm3 = xptr(srcreg+2*d1+16);	/* Load R2 */ \
	xptr(srcreg+16) = xmm0;	/* Save R1 */ \
	xptr(srcreg+2*d1) = xmm7;	/* Save R2 */ \
	xptr(srcreg+2*d1+16) = xmm1;	/* Save R3 */ \
	xptr(srcreg+2*d1+32) = xmm5;	/* Save R4 */ \
	xptr(srcreg+2*d1+48) = xmm6;	/* Save R5 */ \
	xmm0 = xptr(srcreg+48);	/* Load R3 */ \
	xmm1 = xptr(srcreg+4*d1+32);	/* Load R4 */ \
	xmm4 = xptr(srcreg+4*d1+48);	/* Load R5 */ \
	x5r_fft(xmm2, xmm3, xmm0, xmm1, xmm4, xmm5, xmm6, xmm7); \
	xptr(srcreg+48) = xmm2;	/* Save R1 */ \
	xptr(srcreg+4*d1) = xmm7;	/* Save R2 */ \
	xptr(srcreg+4*d1+16) = xmm3;	/* Save R3 */ \
	xptr(srcreg+4*d1+32) = xmm5;	/* Save R4 */ \
	xptr(srcreg+4*d1+48) = xmm6;	/* Save R5 */ \
	srcreg += srcinc; \
\

#define g5cl_five_reals_first_fft(srcreg,srcinc,d1,dstreg,dstinc,e1) { \
	uintptr_t src_rbx = srcreg+rbx; \
	xmm0 = xptr(src_rbx); \
	xmm1 = xptr(src_rbx+3*d1); \
	xmm2 = xptr(src_rbx+3*d1+16); \
	xmm3 = xptr(src_rbx+d1+32); \
	xmm4 = xptr(src_rbx+d1+48); \
	x5r_fft(xmm0, xmm1, xmm2, xmm3, xmm4, xmm5, xmm6, xmm7); \
	xmm2 = xptr(src_rbx+d1);	/* Load R1 */ \
	xmm3 = xptr(src_rbx+d1+16); /* Load R2 */ \
	xptr(dstreg) = xmm0;		/* Save R1 */ \
	xptr(dstreg+e1) = xmm7;	/* Save R2 */ \
	xptr(dstreg+e1+16) = xmm1;	/* Save R3 */ \
	xptr(dstreg+e1+32) = xmm5;	/* Save R4 */ \
	xptr(dstreg+e1+48) = xmm6;	/* Save R5 */ \
	xmm0 = xptr(src_rbx+32);	/* Load R3 */ \
	xmm1 = xptr(src_rbx+3*d1+32); /* Load R4 */ \
	xmm4 = xptr(src_rbx+3*d1+48); /* Load R5 */ \
	x5r_fft(xmm2, xmm3, xmm0, xmm1, xmm4, xmm5, xmm6, xmm7); \
	xptr(dstreg+32) = xmm2;	/* Save R1 */ \
	xptr(dstreg+3*e1) = xmm7;	/* Save R2 */ \
	xptr(dstreg+3*e1+16) = xmm3;	/* Save R3 */ \
	xptr(dstreg+3*e1+32) = xmm5;	/* Save R4 */ \
	xptr(dstreg+3*e1+48) = xmm6;	/* Save R5 */ \
\
	xmm0 = xptr(src_rbx+16); \
	xmm1 = xptr(src_rbx+4*d1); \
	xmm2 = xptr(src_rbx+4*d1+16); \
	xmm3 = xptr(src_rbx+2*d1+32); \
	xmm4 = xptr(src_rbx+2*d1+48); \
	x5r_fft(xmm0, xmm1, xmm2, xmm3, xmm4, xmm5, xmm6, xmm7); \
	xmm2 = xptr(src_rbx+2*d1); /* Load R1 */ \
	xmm3 = xptr(src_rbx+2*d1+16); /* Load R2 */ \
	xptr(dstreg+16) = xmm0;	/* Save R1 */ \
	xptr(dstreg+2*e1) = xmm7;	/* Save R2 */ \
	xptr(dstreg+2*e1+16) = xmm1;	/* Save R3 */ \
	xptr(dstreg+2*e1+32) = xmm5;	/* Save R4 */ \
	xptr(dstreg+2*e1+48) = xmm6;	/* Save R5 */ \
	xmm0 = xptr(src_rbx+48);	/* Load R3 */ \
	xmm1 = xptr(src_rbx+4*d1+32); /* Load R4 */ \
	xmm4 = xptr(src_rbx+4*d1+48); /* Load R5 */ \
	srcreg += srcinc; \
	x5r_fft(xmm2, xmm3, xmm0, xmm1, xmm4, xmm5, xmm6, xmm7); \
	xptr(dstreg+48) = xmm2;	/* Save R1 */ \
	xptr(dstreg+4*e1) = xmm7;	/* Save R2 */ \
	xptr(dstreg+4*e1+16) = xmm3;	/* Save R3 */ \
	xptr(dstreg+4*e1+32) = xmm5;	/* Save R4 */ \
	xptr(dstreg+4*e1+48) = xmm6;	/* Save R5 */ \
	dstreg += dstinc; \
}

#define x5r_fft(r1, r2, r3, r4, r5, t1, t2, t3) \
	t1 = r5;			/* 0-5 Copy R5 */ \
	r5 += r2;			/* 1-4 T1 = R2 + R5 */ \
	t2 = r4;			/* 2-7 Copy R4 */ \
	r4 += r3;			/* 3-5 T2 = R3 + R4 */ \
	t3 = r1;			/* 4-9 newR2 = R1 */ \
	r2 -= t1;			/* 6-9 T3 = R2 - R5 */ \
	t1 = r1;			/* 7-12 newR3 = R1 */ \
	r3 -= t2;			/* 8-11 T4 = R3 - R4 */ \
	t2 = CONST2_P618;		/* 9-14 const1 (.588/.951) */ \
	r1 += r5;			/* 10-13 newR1 = R1 + T1 */ \
	r5 *= CONST2_P309;		/* 11-16 T1 = T1 * .309 */ \
	r2 *= CONST2_P951;		/* 13-18 T3 = T3 * .951 (new I2) */ \
	r1 += r4;			/* 14-17 newR1 = newR1 + T2 */ \
	r3 *= CONST2_P588;		/* 15-20 T4 = T4 * .588 */ \
	t3 += r5;			/* 17-20 newR2 = newR2 + T1 */ \
	r4 *= CONST2_M809;		/* 18-23 T2 = T2 * -.809 */ \
	r5 *= CONST2_M262;		/* 20-25 T1 = T1 * (-.809/.309) */ \
	t2 *= r2;			/* 22-27 T3 = T3 * (.588/.951) */ \
	r2 += r3;			/* 23-26 newI2 = newI2 + T4 */ \
	r3 *= CONST2_M162;		/* 24-29 T4 = T4 * (-.951/.588) */ \
	t3 += r4;			/* 25-28 newR2 = newR2 + T2 */ \
	r4 *= CONST2_M382;		/* 26-31 T2 = T2 * (.309/-.809) */ \
	t1 += r5;			/* 27-30 newR3 = newR3 + T1 */ \
	t2 += r3;			/* 30-33 T3 = T3 + T4 (final I3) */ \
	t1 += r4;			/* 32-35 newR3 = newR3 + T2 */ \
\
\

#define xfive_reals_unfft_preload \
\

#define x5cl_five_reals_last_unfft(srcreg,srcinc,d1) \
	xmm0 = xptr(srcreg);		/* R1 */ \
	xmm1 = xptr(srcreg+d1);	/* R2 */ \
	xmm2 = xptr(srcreg+d1+32);	/* R3 */ \
	xmm3 = xptr(srcreg+3*d1);	/* R4 */ \
	xmm4 = xptr(srcreg+3*d1+32);	/* R5 */ \
	x5r_unfft(xmm0, xmm1, xmm2, xmm3, xmm4, xmm5, xmm6, xmm7, xptr(srcreg)); \
	xmm1 = xptr(srcreg+d1+16);	/* R2 */ \
	xmm2 = xptr(srcreg+d1+48);	/* R3 */ \
	xptr(srcreg+d1+32) = xmm3;	/* Save R4 */ \
	xptr(srcreg+d1+48) = xmm7;	/* Save R5 */ \
	xmm3 = xptr(srcreg+3*d1+16);	/* R4 */ \
	xmm4 = xptr(srcreg+3*d1+48);	/* R5 */ \
	xptr(srcreg+3*d1) = xmm5;	/* Save R2 */ \
	xptr(srcreg+3*d1+16) = xmm6;	/* Save R3 */ \
	xmm0 = xptr(srcreg+32);	/* R1 */ \
	x5r_unfft(xmm0, xmm1, xmm2, xmm3, xmm4, xmm5, xmm6, xmm7, xptr(srcreg+d1)); \
	xptr(srcreg+d1+16) = xmm5;	/* Save R2 */ \
	xptr(srcreg+32) = xmm6;	/* Save R3 */ \
	xptr(srcreg+3*d1+32) = xmm3;	/* Save R4 */ \
	xptr(srcreg+3*d1+48) = xmm7;	/* Save R5 */ \
\
	xmm0 = xptr(srcreg+16);	/* R1 */ \
	xmm1 = xptr(srcreg+2*d1);	/* R2 */ \
	xmm2 = xptr(srcreg+2*d1+32);	/* R3 */ \
	xmm3 = xptr(srcreg+4*d1);	/* R4 */ \
	xmm4 = xptr(srcreg+4*d1+32);	/* R5 */ \
	x5r_unfft(xmm0, xmm1, xmm2, xmm3, xmm4, xmm5, xmm6, xmm7, xptr(srcreg+16)); \
	xmm1 = xptr(srcreg+2*d1+16);	/* R2 */ \
	xmm2 = xptr(srcreg+2*d1+48);	/* R3 */ \
	xptr(srcreg+2*d1+32) = xmm3;	/* Save R4 */ \
	xptr(srcreg+2*d1+48) = xmm7;	/* Save R5 */ \
	xmm3 = xptr(srcreg+4*d1+16);	/* R4 */ \
	xmm4 = xptr(srcreg+4*d1+48);	/* R5 */ \
	xptr(srcreg+4*d1) = xmm5;	/* Save R2 */ \
	xptr(srcreg+4*d1+16) = xmm6;	/* Save R3 */ \
	xmm0 = xptr(srcreg+48);	/* R1 */ \
	x5r_unfft(xmm0, xmm1, xmm2, xmm3, xmm4, xmm5, xmm6, xmm7, xptr(srcreg+2*d1)); \
	xptr(srcreg+2*d1+16) = xmm5;	/* Save R2 */ \
	xptr(srcreg+48) = xmm6;	/* Save R3 */ \
	xptr(srcreg+4*d1+32) = xmm3;	/* Save R4 */ \
	xptr(srcreg+4*d1+48) = xmm7;	/* Save R5 */ \
	srcreg += srcinc; \
\

#define x5r_unfft(r1, r2, r3, r4, r5, t1, t2, t3, mem1) \
	t1 = CONST2_P309;		/* Load .309 */ \
	t1 *= r2;			/* 1-6 R2*.309 */ \
	t2 = CONST2_M809;		/* Load -.809 */ \
	t2 *= r2;			/* 3-8 R2*-.809 */ \
	r2 += r4;			/* 4-7 R2+R3 */ \
	t3 = CONST2_M809;		/* Load .309  BUG in comment */ \
	t3 *= r4;			/* 5-10 R3*-.809 */ \
	r2 += r1;			/* 6-9 R1+R2+R3 (final R1) */ \
	r4 *= CONST2_P309;		/* 7-12 R3*.309 */ \
	t1 += r1;			/* 8-11 R1 + R2*.309 */ \
	mem1 = r2;		/* Save final R1 */ \
	r2 = CONST2_P951;		/* Load 0.951 */ \
	r2 *= r3;	 		/* 9-14 I2*.951 */ \
	t2 += r1;			/* 10-13 R1 + R2*-.809 */ \
	r1 = CONST2_P588;		/* Load 0.588 */ \
	r1 *= r5;			/* 11-16 I3*.588 */ \
	t1 += t3;			/* 12-15 R1 + R2*.309 - R3*.809 */ \
	r3 *= CONST2_P588;		/* 13-18 I2*.588 */ \
	t2 += r4;			/* 14-17 R1 - R2*.809 + R3*.309 */ \
	r5 *= CONST2_P951;		/* 15-20 I3*-.951 */ \
	t3 = t1;			/* 16-21 R1 + R2*.309 - R3*.809 */ \
	r2 += r1;			/* 17-20 I2*.951 + I3*.588 */ \
	r4 = t2;			/* 18-23 R1 - R2*.809 + R3*.309 */ \
	r3 -= r5;			/* 21-24 I2*.588 - I3*.951 */ \
	t1 += r2;			/* 23-26 final R2 */ \
	t3 -= r2;			/* 25-28 final R5 */ \
	t2 += r3;			/* 27-30 final R3 */ \
	r4 -= r3;			/* 29-31 final R4 */ \
\
\

#define xsix_reals_fft_preload \
\

#define s3cl_six_reals_first_fft(srcreg,srcinc,d1) { \
	uintptr_t src_rbx = srcreg+rbx; \
	low_load(xmm0, xptr(src_rbx), xptr(src_rbx+16)); /* R1 */ \
	low_load(xmm3, xptr(src_rbx+32), xptr(src_rbx+48)); /* R4 */ \
	high_load(xmm2, xptr(src_rbx+d1), xptr(src_rbx+d1+16)); /* R3 */ \
	high_load(xmm5, xptr(src_rbx+d1+32), xptr(src_rbx+d1+48)); /* R6 */ \
	low_load(xmm1, xptr(src_rbx+2*d1), xptr(src_rbx+2*d1+16)); /* R2 */ \
	low_load(xmm4, xptr(src_rbx+2*d1+32), xptr(src_rbx+2*d1+48)); /* R5 */ \
	x6r_fft(xmm0, xmm1, xmm2, xmm3, xmm4, xmm5, xmm6, xmm7); \
	high_load(xmm3, xptr(src_rbx), xptr(src_rbx+16)); /* R2 */ \
	high_load(xmm4, xptr(src_rbx+32), xptr(src_rbx+48)); /* R5 */ \
	xptr(srcreg) = xmm1;			/* Save R1 */ \
	xptr(srcreg+32) = xmm2;		/* Save R2 */ \
	low_load(xmm1, xptr(src_rbx+d1), xptr(src_rbx+d1+16)); /* R1 */ \
	low_load(xmm2, xptr(src_rbx+d1+32), xptr(src_rbx+d1+48));	/* R4 */ \
	xptr(srcreg+d1) = xmm5;		/* Save R3 */ \
	xptr(srcreg+d1+16) = xmm7;		/* Save R4 */ \
	xptr(srcreg+d1+32) = xmm0;		/* Save R5 */ \
	xptr(srcreg+d1+48) = xmm6;		/* Save R6 */ \
	high_load(xmm5, xptr(src_rbx+2*d1), xptr(src_rbx+2*d1+16)); /* R3 */ \
	high_load(xmm7, xptr(src_rbx+2*d1+32), xptr(src_rbx+2*d1+48)); /* R6 */ \
	x6r_fft(xmm1, xmm3, xmm5, xmm2, xmm4, xmm7, xmm0, xmm6); \
	xptr(srcreg+16) = xmm3;		/* Save R1 */ \
	xptr(srcreg+48) = xmm5;		/* Save R2 */ \
	xptr(srcreg+2*d1) = xmm7;		/* Save R3 */ \
	xptr(srcreg+2*d1+16) = xmm6;		/* Save R4 */ \
	xptr(srcreg+2*d1+32) = xmm1;		/* Save R5 */ \
	xptr(srcreg+2*d1+48) = xmm0;		/* Save R6 */ \
	srcreg += srcinc; \
}

#define x3cl_six_reals_first_fft(srcreg,srcinc,d1) { \
	uintptr_t src_rbx = srcreg+rbx; \
	xmm0 = xptr(src_rbx); \
	xmm3 = xptr(src_rbx+32); \
	xmm2 = xptr(src_rbx+d1+16); \
	xmm5 = xptr(src_rbx+d1+48); \
	xmm1 = xptr(src_rbx+2*d1); \
	xmm4 = xptr(src_rbx+2*d1+32); \
	x6r_fft(xmm0, xmm1, xmm2, xmm3, xmm4, xmm5, xmm6, xmm7); \
	xptr(srcreg) = xmm1; \
	xptr(srcreg+32) = xmm2; \
	xptr(srcreg+d1+16) = xmm7; \
	xptr(srcreg+d1+48) = xmm6; \
	xmm1 = xptr(src_rbx+d1);		/* R1 */ \
	xmm2 = xptr(src_rbx+d1+32);	/* R4 */ \
	xptr(srcreg+d1) = xmm5; \
	xptr(srcreg+d1+32) = xmm0; \
	xmm3 = xptr(src_rbx+16);		/* R2 */ \
	xmm4 = xptr(src_rbx+48);		/* R5 */ \
	xmm5 = xptr(src_rbx+2*d1+16);	/* R3 */ \
	xmm6 = xptr(src_rbx+2*d1+48);	/* R6 */ \
	x6r_fft(xmm1, xmm3, xmm5, xmm2, xmm4, xmm6, xmm0, xmm7); \
	xptr(srcreg+16) = xmm3; \
	xptr(srcreg+48) = xmm5; \
	xptr(srcreg+2*d1) = xmm6; \
	xptr(srcreg+2*d1+16) = xmm7; \
	xptr(srcreg+2*d1+32) = xmm1; \
	xptr(srcreg+2*d1+48) = xmm0; \
	srcreg += srcinc; \
}

#define x3cl_six_reals_fft(srcreg,srcinc,d1) \
	xmm0 = xptr(srcreg); \
	xmm3 = xptr(srcreg+32); \
	xmm2 = xptr(srcreg+d1+16); \
	xmm5 = xptr(srcreg+d1+48); \
	xmm1 = xptr(srcreg+2*d1); \
	xmm4 = xptr(srcreg+2*d1+32); \
	x6r_fft(xmm0, xmm1, xmm2, xmm3, xmm4, xmm5, xmm6, xmm7); \
	xptr(srcreg) = xmm1; \
	xptr(srcreg+32) = xmm2; \
	xptr(srcreg+d1+16) = xmm7; \
	xptr(srcreg+d1+48) = xmm6; \
	xmm1 = xptr(srcreg+d1);		/* R1 */ \
	xmm2 = xptr(srcreg+d1+32);		/* R4 */ \
	xptr(srcreg+d1) = xmm5; \
	xptr(srcreg+d1+32) = xmm0; \
	xmm3 = xptr(srcreg+16);		/* R2 */ \
	xmm4 = xptr(srcreg+48);		/* R5 */ \
	xmm5 = xptr(srcreg+2*d1+16);		/* R3 */ \
	xmm6 = xptr(srcreg+2*d1+48);		/* R6 */ \
	x6r_fft(xmm1, xmm3, xmm5, xmm2, xmm4, xmm6, xmm0, xmm7); \
	xptr(srcreg+16) = xmm3; \
	xptr(srcreg+48) = xmm5; \
	xptr(srcreg+2*d1) = xmm6; \
	xptr(srcreg+2*d1+16) = xmm7; \
	xptr(srcreg+2*d1+32) = xmm1; \
	xptr(srcreg+2*d1+48) = xmm0; \
	srcreg += srcinc; \
\

#define g3cl_six_reals_first_fft(srcreg,srcinc,d1,dstreg,dstinc,e1) { \
	uintptr_t src_rbx = srcreg+rbx; \
	xmm0 = xptr(src_rbx); \
	xmm3 = xptr(src_rbx+32); \
	xmm2 = xptr(src_rbx+d1+16); \
	xmm5 = xptr(src_rbx+d1+48); \
	xmm1 = xptr(src_rbx+2*d1); \
	xmm4 = xptr(src_rbx+2*d1+32); \
	x6r_fft(xmm0, xmm1, xmm2, xmm3, xmm4, xmm5, xmm6, xmm7); \
	xptr(dstreg) = xmm1; \
	xptr(dstreg+32) = xmm2; \
	xptr(dstreg+e1+16) = xmm7; \
	xptr(dstreg+e1+48) = xmm6; \
	xmm1 = xptr(src_rbx+d1);		/* R1 */ \
	xmm2 = xptr(src_rbx+d1+32);	/* R4 */ \
	xptr(dstreg+e1) = xmm5; \
	xptr(dstreg+e1+32) = xmm0; \
	xmm3 = xptr(src_rbx+16);		/* R2 */ \
	xmm4 = xptr(src_rbx+48);		/* R5 */ \
	xmm5 = xptr(src_rbx+2*d1+16);	/* R3 */ \
	xmm6 = xptr(src_rbx+2*d1+48);	/* R6 */ \
	srcreg += srcinc; \
	x6r_fft(xmm1, xmm3, xmm5, xmm2, xmm4, xmm6, xmm0, xmm7); \
	xptr(dstreg+16) = xmm3; \
	xptr(dstreg+48) = xmm5; \
	xptr(dstreg+2*e1) = xmm6; \
	xptr(dstreg+2*e1+16) = xmm7; \
	xptr(dstreg+2*e1+32) = xmm1; \
	xptr(dstreg+2*e1+48) = xmm0; \
	dstreg += dstinc; \
}
/* Simplifying the pseudo code from pfa.mac yields: */ \
/* new R1 = R1 + R3 + R5 */ \
/* new R2 = R1 - 0.5 * (R3 + R5) */ \
/* new R3 = R2 + R4 + R6 */ \
/* new I2 = 0.866 * (R3 - R5) */ \
/* new R4 = 0.5 * (R2 + R6) - R4 */ \
/* new I4 = 0.866 * (R2 - R6) */ \
/* R1 + R3 (final R1) */ \
/* R1 - R3 (final R2) */ \
/* R2 + R4 (final R3) */ \
/* R2 - R4 (final R5) */ \
/* I2 + I4 (final R4) */ \
/* I2 - I4 (final R6) */ \

#define x6r_fft(r1, r2, r3, r4, r5, r6, t1, t2) \
	t1 = r3; \
	r3 += r5;			/* T4 = R3 + R5 */ \
	t2 = r2; \
	r2 += r6;			/* T2 = R2 + R6 */ \
	t1 -= r5;			/* T3 = R3 - R5 */ \
	t2 -= r6;			/* T1 = R2 - R6 */ \
	r5 = r3 * 0.5;			/* 0.5 * (R3 + R5) */ \
	r3 += r1;			/* new R1 = R1 + R3 + R5 */ \
	r6 = r2 * 0.5;			/* 0.5 * (R2 + R6) */ \
	r2 += r4;			/* new R3 = R2 + R4 + R6 */ \
	t1 *= CONST2_P866;		/* new I2 = 0.866 * (R3 - R5) */ \
	r1 -= r5;			/* new R2 = R1 - 0.5 * (R3 + R5) */ \
	t2 *= CONST2_P866;		/* new I4 = 0.866 * (R2 - R6) */ \
	r6 -= r4;			/* new R4 = 0.5 * (R2 + R6) - R4 */ \
\
	r5 = r3; \
	r3 -= r2;			/* R1 = R1 - R3 (final R2) */ \
	r2 += r5;			/* R3 = R1 + R3 (final R1) */ \
	r4 = t1; \
	t1 -= t2;			/* I2 = I2 - I4 (final R6) */ \
	t2 += r4;			/* I4 = I2 + I4 (final R4) */ \
	r5 = r1; \
	r1 -= r6;			/* R2 = R2 - R4 (final R5) */ \
	r6 += r5;			/* R4 = R2 + R4 (final R3) */ \
\
\

#define xsix_reals_unfft_preload \
\

#define x3cl_six_reals_last_unfft(srcreg,srcinc,d1) \
	xmm0 = xptr(srcreg);		/* R1 */ \
	xmm1 = xptr(srcreg+32);	/* R2 */ \
	xmm2 = xptr(srcreg+d1);	/* R3 */ \
	xmm3 = xptr(srcreg+d1+32);	/* R4 */ \
	xmm4 = xptr(srcreg+2*d1);	/* R5 */ \
	xmm5 = xptr(srcreg+2*d1+32);	/* R6 */ \
	x6r_unfft(xmm0, xmm1, xmm2, xmm3, xmm4, xmm5, xmm6, xmm7); \
	xptr(srcreg) = xmm4;		/* Save R1 */ \
	xptr(srcreg+32) = xmm0;	/* Save R4 */ \
	xptr(srcreg+2*d1) = xmm3;	/* Save R2 */ \
	xptr(srcreg+2*d1+32) = xmm1;	/* Save R5 */ \
	xmm2 = xptr(srcreg+d1+16);	/* R3 */ \
	xmm3 = xptr(srcreg+d1+48);	/* R4 */ \
	xptr(srcreg+d1+16) = xmm5;	/* Save R3 */ \
	xptr(srcreg+d1+48) = xmm7;	/* Save R6 */ \
	xmm0 = xptr(srcreg+16);	/* R1 */ \
	xmm1 = xptr(srcreg+48);	/* R2 */ \
	xmm4 = xptr(srcreg+2*d1+16);	/* R5 */ \
	xmm5 = xptr(srcreg+2*d1+48);	/* R6 */ \
	x6r_unfft(xmm0, xmm1, xmm2, xmm3, xmm4, xmm5, xmm6, xmm7); \
	xptr(srcreg+16) = xmm3;	/* Save R2 */ \
	xptr(srcreg+48) = xmm1;	/* Save R5 */ \
	xptr(srcreg+d1) = xmm4;	/* Save R1 */ \
	xptr(srcreg+d1+32) = xmm0;	/* Save R4 */ \
	xptr(srcreg+2*d1+16) = xmm5;	/* Save R3 */ \
	xptr(srcreg+2*d1+48) = xmm7;	/* Save R6 */ \
	srcreg += srcinc; \
\
/* Simplifying the pseudo code in pfa.mac we get: */ \
/* R1 + R2 (new R1) */ \
/* R1 - R2 (new R3) */ \
/* R3 + R5 (new R2) */ \
/* R3 - R5 (new R4) */ \
/* R4 + R6 (new I2) */ \
/* R4 - R6 (new I4) */ \
/* final R1 = R1 + R2 */ \
/* final R3 = R1 - 0.5 * R2 + 0.866 * I2 */ \
/* final R5 = R1 - 0.5 * R2 - 0.866 * I2 */ \
/* final R2 = R3 + 0.5 * R4 + 0.866 * I4 */ \
/* final R4 = R3 - R4 */ \
/* final R6 = R3 + 0.5 * R4 - 0.866 * I4 */ \

#define x6r_unfft(r1, r2, r3, r4, r5, r6, t1, t2) \
	t1 = r3; \
	r3 -= r5;			/* R3 - R5 (new R4) */ \
	r5 += t1;			/* R3 + R5 (new R2) */ \
	t2 = r4; \
	r4 -= r6;			/* R4 - R6 (new I4) */ \
	r6 += t2;			/* R4 + R6 (new I2) */ \
	t2 = r3 * 0.5;			/* 0.5 * R4 */ \
	t1 = r1; \
	r1 -= r2;			/* R1 - R2 (new R3) */ \
	r2 += t1;			/* R1 + R2 (new R1) */ \
	t1 = r5 * 0.5;			/* 0.5 * R2 */ \
	r4 *= CONST2_P866;		/* 0.866 * I4 */ \
	r6 *= CONST2_P866;		/* 0.866 * I2 */ \
	t2 += r1;			/* R3 + 0.5 * R4 */ \
	r5 += r2;			/* final R1 = R1 + R2 */ \
	r2 -= t1;			/* R1 - 0.5 * R2 */ \
	r1 -= r3;			/* final R4 = R3 - R4 */ \
	t1 = r4;			/* Copy 0.866 * I4 */ \
	r4 += t2;			/* final R2 = R3 + 0.5 * R4 + 0.866 * I4 */ \
	t2 -= t1;			/* final R6 = R3 + 0.5 * R4 - 0.866 * I4 */ \
	r3 = r6;			/* Copy 0.866 * I2 */ \
	r6 += r2;			/* final R3 = R1 - 0.5 * R2 + 0.866 * I2 */ \
	r2 -= r3;			/* final R5 = R1 - 0.5 * R2 - 0.866 * I2 */ \
\
\

#define xseven_reals_fft_preload \
\

#define s7cl_seven_reals_first_fft(srcreg,srcinc,d1) { \
	uintptr_t src_rbx = srcreg+rbx; \
	shuffle_load(xmm0,xmm7,xptr(src_rbx),xptr(src_rbx+16)); /* R1,R1 */ \
	xptr(srcreg+16) =xmm7;		/* Save it */ \
	xmm1 = xptr(src_rbx+d1+16); \
	/*movlpd*/ xmm1[0] = f64ptr(src_rbx+d1+8);	/* R2 */ \
	shuffle_load(xmm2,xmm3,xptr(src_rbx+5*d1),xptr(src_rbx+5*d1+16)); /*R3,R4 */ \
	low_load(xmm4, xptr(src_rbx+d1+32), xptr(src_rbx+d1+48));	/* R5 */ \
	shuffle_load(xmm5,xmm6,xptr(src_rbx+3*d1+32),xptr(src_rbx+3*d1+48)); /*R6,R7 */ \
	x7r_fft(xmm0, xmm1, xmm2, xmm3, xmm4, xmm5, xmm6, xmm7, xptr(srcreg)); \
	low_load(xmm2, xptr(src_rbx+d1), xptr(src_rbx+d1+16)); /* R1 */ \
	high_load(xmm6, xptr(src_rbx+d1+32), xptr(src_rbx+d1+48)); /* R5 */ \
	xptr(srcreg+d1) = xmm0;		/* Save R2 */ \
	xptr(srcreg+d1+32) = xmm1;		/* Save R3 */ \
	shuffle_load(xmm1,xmm0,xptr(src_rbx+3*d1),xptr(src_rbx+3*d1+16)); /*R2,R3 */ \
	xptr(srcreg+3*d1) = xmm4;		/* Save R4 */ \
	xptr(srcreg+3*d1+16) = xmm5;		/* Save R5 */ \
	xptr(srcreg+3*d1+32) = xmm7;		/* Save R6 */ \
	xptr(srcreg+3*d1+48) = xmm3;		/* Save R7 */ \
	shuffle_load(xmm3,xmm7,xptr(src_rbx+32),xptr(src_rbx+48)); /* R4,R4 */ \
	xptr(srcreg+48) =xmm7;		/* Save it */ \
	shuffle_load(xmm5,xmm4,xptr(src_rbx+5*d1+32),xptr(src_rbx+5*d1+48));/*R6,R7 */ \
	x7r_fft(xmm2, xmm1, xmm0, xmm3, xmm6, xmm5, xmm4, xmm7, xptr(srcreg+32)); \
	xptr(srcreg+d1+16) = xmm2;		/* Save R2 */ \
	xptr(srcreg+d1+48) = xmm1;		/* Save R3 */ \
	xptr(srcreg+5*d1) = xmm6;		/* Save R4 */ \
	xptr(srcreg+5*d1+16) = xmm5;		/* Save R5 */ \
	xptr(srcreg+5*d1+32) = xmm7;		/* Save R6 */ \
	xptr(srcreg+5*d1+48) = xmm3;		/* Save R7 */ \
\
	xmm0 = xptr(srcreg+16);		/* R1 */ \
	high_load(xmm1, xptr(src_rbx+2*d1), xptr(src_rbx+2*d1+16)); /* R2 */ \
	shuffle_load(xmm2,xmm3,xptr(src_rbx+6*d1),xptr(src_rbx+6*d1+16));/*R3,R4 */ \
	low_load(xmm4, xptr(src_rbx+2*d1+32), xptr(src_rbx+2*d1+48)); /* R5 */ \
	shuffle_load(xmm5,xmm6,xptr(src_rbx+4*d1+32),xptr(src_rbx+4*d1+48));/*R6,R7 */ \
	x7r_fft(xmm0, xmm1, xmm2, xmm3, xmm4, xmm5, xmm6, xmm7, xptr(srcreg+16)); \
	low_load(xmm2, xptr(src_rbx+2*d1), xptr(src_rbx+2*d1+16)); /* R1 */ \
	high_load(xmm6, xptr(src_rbx+2*d1+32), xptr(src_rbx+2*d1+48)); \
	xptr(srcreg+2*d1) = xmm0;		/* Save R2 */ \
	xptr(srcreg+2*d1+32) = xmm1;		/* Save R3 */ \
	shuffle_load(xmm1,xmm0,xptr(src_rbx+4*d1),xptr(src_rbx+4*d1+16));/*R2,R3 */ \
	xptr(srcreg+4*d1) = xmm4;		/* Save R4 */ \
	xptr(srcreg+4*d1+16) = xmm5;		/* Save R5 */ \
	xptr(srcreg+4*d1+32) = xmm7;		/* Save R6 */ \
	xptr(srcreg+4*d1+48) = xmm3;		/* Save R7 */ \
	xmm3 = xptr(srcreg+48);		/* R4 */ \
	shuffle_load(xmm5,xmm4,xptr(src_rbx+6*d1+32),xptr(src_rbx+6*d1+48));/*R6,R7 */ \
	x7r_fft(xmm2, xmm1, xmm0, xmm3, xmm6, xmm5, xmm4, xmm7, xptr(srcreg+48)); \
	xptr(srcreg+2*d1+16) = xmm2;		/* Save R2 */ \
	xptr(srcreg+2*d1+48) = xmm1;		/* Save R3 */ \
	xptr(srcreg+6*d1) = xmm6;		/* Save R4 */ \
	xptr(srcreg+6*d1+16) = xmm5;		/* Save R5 */ \
	xptr(srcreg+6*d1+32) = xmm7;		/* Save R6 */ \
	xptr(srcreg+6*d1+48) = xmm3;		/* Save R7 */ \
	srcreg += srcinc; \
}
/* 215.55 clocks */ \

#define x7cl_seven_reals_first_fft(srcreg,srcinc,d1) { \
	uintptr_t src_rbx = srcreg+rbx; \
	xmm0 = xptr(src_rbx);	/* Load R1 */ \
	xmm1 = xptr(src_rbx+d1+16);/* Load R2 */ \
	xmm2 = xptr(src_rbx+5*d1);/* Load R3 */ \
	x7r_fft_mem(xptr(src_rbx+5*d1+16), xptr(src_rbx+d1+32), xptr(src_rbx+3*d1+32), xptr(src_rbx+3*d1+48), xptr(srcreg), xptr(srcreg+3*d1+32), 0); \
	xptr(srcreg+d1+32) = xmm2;	/* Save I2 */ \
	xmm2 = xptr(src_rbx+3*d1+16);/* Load R3 */ \
	xptr(srcreg+3*d1+16) = xmm1;	/* Save I3 */ \
	xmm1 = xptr(src_rbx+3*d1);/* Load R2 */ \
	xptr(srcreg+3*d1) = xmm0;	/* Save R3 */ \
	xmm0 = xptr(src_rbx+d1);	/* Load R1 */ \
	xptr(srcreg+d1) = xmm4;	/* Save R2 */ \
	xptr(srcreg+3*d1+48) = xmm3;	/* Save I4 */ \
	x7r_fft_mem(xptr(src_rbx+32), xptr(src_rbx+d1+48), xptr(src_rbx+5*d1+32), xptr(src_rbx+5*d1+48), xptr(srcreg+32), xptr(srcreg+5*d1+32), 1); \
	xptr(srcreg+d1+16) = xmm4;	/* Save R2 */ \
	xptr(srcreg+d1+48) = xmm2;	/* Save I2 */ \
	xptr(srcreg+5*d1) = xmm0;	/* Save R3 */ \
	xptr(srcreg+5*d1+16) = xmm1;	/* Save I3 */ \
	xptr(srcreg+5*d1+48) = xmm3;	/* Save I4 */ \
\
	xmm0 = xptr(src_rbx+16);	/* Load R1 */ \
	xmm1 = xptr(src_rbx+2*d1+16);/* Load R2 */ \
	xmm2 = xptr(src_rbx+6*d1);/* Load R3 */ \
	x7r_fft_mem(xptr(src_rbx+6*d1+16), xptr(src_rbx+2*d1+32), xptr(src_rbx+4*d1+32), xptr(src_rbx+4*d1+48), xptr(srcreg+16), xptr(srcreg+4*d1+32), 0); \
	xptr(srcreg+2*d1+32) = xmm2;	/* Save I2 */ \
	xmm2 = xptr(src_rbx+4*d1+16);/* Load R3 */ \
	xptr(srcreg+4*d1+16) = xmm1;	/* Save I3 */ \
	xmm1 = xptr(src_rbx+4*d1);/* Load R2 */ \
	xptr(srcreg+4*d1) = xmm0;	/* Save R3 */ \
	xmm0 = xptr(src_rbx+2*d1);/* Load R1 */ \
	xptr(srcreg+2*d1) = xmm4;	/* Save R2 */ \
	xptr(srcreg+4*d1+48) = xmm3;	/* Save I4 */ \
	x7r_fft_mem(xptr(src_rbx+48), xptr(src_rbx+2*d1+48), xptr(src_rbx+6*d1+32), xptr(src_rbx+6*d1+48), xptr(srcreg+48), xptr(srcreg+6*d1+32), 1); \
	xptr(srcreg+2*d1+16) = xmm4;	/* Save R2 */ \
	xptr(srcreg+2*d1+48) = xmm2;	/* Save I2 */ \
	xptr(srcreg+6*d1) = xmm0;	/* Save R3 */ \
	xptr(srcreg+6*d1+16) = xmm1;	/* Save I3 */ \
	xptr(srcreg+6*d1+48) = xmm3;	/* Save I4 */ \
	srcreg += srcinc; \
}

#define x7cl_seven_reals_fft(srcreg,srcinc,d1) \
	xmm0 = xptr(srcreg);		/* Load R1 */ \
	xmm1 = xptr(srcreg+d1+16);	/* Load R2 */ \
	xmm2 = xptr(srcreg+5*d1);	/* Load R3 */ \
	x7r_fft_mem(xptr(srcreg+5*d1+16), xptr(srcreg+d1+32), xptr(srcreg+3*d1+32), xptr(srcreg+3*d1+48), xptr(srcreg), xptr(srcreg+3*d1+32), 0); \
	xptr(srcreg+d1+32) = xmm2;	/* Save I2 */ \
	xmm2 = xptr(srcreg+3*d1+16);	/* Load R3 */ \
	xptr(srcreg+3*d1+16) = xmm1;	/* Save I3 */ \
	xmm1 = xptr(srcreg+3*d1);	/* Load R2 */ \
	xptr(srcreg+3*d1) = xmm0;	/* Save R3 */ \
	xmm0 = xptr(srcreg+d1);	/* Load R1 */ \
	xptr(srcreg+d1) = xmm4;	/* Save R2 */ \
	xptr(srcreg+3*d1+48) = xmm3;	/* Save I4 */ \
	x7r_fft_mem(xptr(srcreg+32), xptr(srcreg+d1+48), xptr(srcreg+5*d1+32), xptr(srcreg+5*d1+48), xptr(srcreg+32), xptr(srcreg+5*d1+32), 1); \
	xptr(srcreg+d1+16) = xmm4;	/* Save R2 */ \
	xptr(srcreg+d1+48) = xmm2;	/* Save I2 */ \
	xptr(srcreg+5*d1) = xmm0;	/* Save R3 */ \
	xptr(srcreg+5*d1+16) = xmm1;	/* Save I3 */ \
	xptr(srcreg+5*d1+48) = xmm3;	/* Save I4 */ \
\
	xmm0 = xptr(srcreg+16);	/* Load R1 */ \
	xmm1 = xptr(srcreg+2*d1+16);	/* Load R2 */ \
	xmm2 = xptr(srcreg+6*d1);	/* Load R3 */ \
	x7r_fft_mem(xptr(srcreg+6*d1+16), xptr(srcreg+2*d1+32), xptr(srcreg+4*d1+32), xptr(srcreg+4*d1+48), xptr(srcreg+16), xptr(srcreg+4*d1+32), 0); \
	xptr(srcreg+2*d1+32) = xmm2;	/* Save I2 */ \
	xmm2 = xptr(srcreg+4*d1+16);	/* Load R3 */ \
	xptr(srcreg+4*d1+16) = xmm1;	/* Save I3 */ \
	xmm1 = xptr(srcreg+4*d1);	/* Load R2 */ \
	xptr(srcreg+4*d1) = xmm0;	/* Save R3 */ \
	xmm0 = xptr(srcreg+2*d1);	/* Load R1 */ \
	xptr(srcreg+2*d1) = xmm4;	/* Save R2 */ \
	xptr(srcreg+4*d1+48) = xmm3;	/* Save I4 */ \
	x7r_fft_mem(xptr(srcreg+48), xptr(srcreg+2*d1+48), xptr(srcreg+6*d1+32), xptr(srcreg+6*d1+48), xptr(srcreg+48), xptr(srcreg+6*d1+32), 1); \
	xptr(srcreg+2*d1+16) = xmm4;	/* Save R2 */ \
	xptr(srcreg+2*d1+48) = xmm2;	/* Save I2 */ \
	xptr(srcreg+6*d1) = xmm0;	/* Save R3 */ \
	xptr(srcreg+6*d1+16) = xmm1;	/* Save I3 */ \
	xptr(srcreg+6*d1+48) = xmm3;	/* Save I4 */ \
	srcreg += srcinc; \
\

#define g7cl_seven_reals_first_fft(srcreg,srcinc,d1,dstreg,dstinc,e1) { \
	uintptr_t src_rbx = srcreg+rbx; \
	xmm0 = xptr(src_rbx);	/* Load R1 */ \
	xmm1 = xptr(src_rbx+d1+16);/* Load R2 */ \
	xmm2 = xptr(src_rbx+5*d1);/* Load R3 */ \
	x7r_fft_mem(xptr(src_rbx+5*d1+16), xptr(src_rbx+d1+32), xptr(src_rbx+3*d1+32), xptr(src_rbx+3*d1+48), xptr(dstreg), xptr(dstreg+3*e1+32), 0); \
	xptr(dstreg+e1+32) = xmm2;	/* Save I2 */ \
	xmm2 = xptr(src_rbx+3*d1+16);/* Load R3 */ \
	xptr(dstreg+3*e1+16) = xmm1;	/* Save I3 */ \
	xmm1 = xptr(src_rbx+3*d1);/* Load R2 */ \
	xptr(dstreg+3*e1) = xmm0;	/* Save R3 */ \
	xmm0 = xptr(src_rbx+d1);	/* Load R1 */ \
	xptr(dstreg+e1) = xmm4;	/* Save R2 */ \
	xptr(dstreg+3*e1+48) = xmm3;	/* Save I4 */ \
	x7r_fft_mem(xptr(src_rbx+32), xptr(src_rbx+d1+48), xptr(src_rbx+5*d1+32), xptr(src_rbx+5*d1+48), xptr(dstreg+32), xptr(dstreg+5*e1+32), 1); \
	xptr(dstreg+e1+16) = xmm4;	/* Save R2 */ \
	xptr(dstreg+e1+48) = xmm2;	/* Save I2 */ \
	xptr(dstreg+5*e1) = xmm0;	/* Save R3 */ \
	xptr(dstreg+5*e1+16) = xmm1;	/* Save I3 */ \
	xptr(dstreg+5*e1+48) = xmm3;	/* Save I4 */ \
\
	xmm0 = xptr(src_rbx+16);	/* Load R1 */ \
	xmm1 = xptr(src_rbx+2*d1+16);/* Load R2 */ \
	xmm2 = xptr(src_rbx+6*d1);/* Load R3 */ \
	x7r_fft_mem(xptr(src_rbx+6*d1+16), xptr(src_rbx+2*d1+32), xptr(src_rbx+4*d1+32), xptr(src_rbx+4*d1+48), xptr(dstreg+16), xptr(dstreg+4*d1+32), 0); \
	xptr(dstreg+2*e1+32) = xmm2;	/* Save I2 */ \
	xmm2 = xptr(src_rbx+4*d1+16);/* Load R3 */ \
	xptr(dstreg+4*e1+16) = xmm1;	/* Save I3 */ \
	xmm1 = xptr(src_rbx+4*d1);/* Load R2 */ \
	xptr(dstreg+4*e1) = xmm0;	/* Save R3 */ \
	xmm0 = xptr(src_rbx+2*d1);/* Load R1 */ \
	xptr(dstreg+2*e1) = xmm4;	/* Save R2 */ \
	xptr(dstreg+4*e1+48) = xmm3;	/* Save I4 */ \
	x7r_fft_mem(xptr(src_rbx+48), xptr(src_rbx+2*d1+48), xptr(src_rbx+6*d1+32), xptr(src_rbx+6*d1+48), xptr(dstreg+48), xptr(dstreg+6*d1+32), 1); \
	srcreg += srcinc; \
	xptr(dstreg+2*e1+16) = xmm4;	/* Save R2 */ \
	xptr(dstreg+2*e1+48) = xmm2;	/* Save I2 */ \
	xptr(dstreg+6*e1) = xmm0;	/* Save R3 */ \
	xptr(dstreg+6*e1+16) = xmm1;	/* Save I3 */ \
	xptr(dstreg+6*e1+48) = xmm3;	/* Save I4 */ \
	dstreg += dstinc; \
}
/* xmm0, xmm1, xmm2 are preloaded with m1, m2, m3 */ \
/* destr1 may be same address as m4 (in that case set m4_conflict to one) */ \
/* destr4 may be same address as m6 */ \

#define x7r_fft_mem(m4, m5, m6, m7, destr1, destr4, m4_conflict) { \
	xmm3 = m7;		/*	T1 = R7 */ \
	xmm3 += xmm1;		/*1-4	T1 = R2+R7 */ \
	xmm4 = m6;		/*	T2 = R6 */ \
	xmm4 += xmm2;		/*3-6	T2 = R3+R6 */ \
	xmm5 = m4;		/*	T3 = R4 */ \
	xmm5 += m5;		/*5-8	T3 = R4+R5 */ \
	xmm6 = CONST2_P623; \
	xmm6 *= xmm3;		/*6-11	T1 = T1 * .623 */ \
	xmm3 += xmm0;		/*7-10  R1+T1 */ \
	xmm7 = CONST2_P623; \
	xmm7 *= xmm4;		/*8-13	T2 = T2 * .623 */ \
	xmm4 += xmm5;		/*9-12	T2+T3 */ \
	xmm5 *= CONST2_P623;		/*10-15	T3 = T3 * .623 */ \
	xmm3 += xmm4;		/*14-17	R1+T1+T2+T3 (final R1) */ \
\
	xmm4 = xmm0;		/*	newR2 = R1 */ \
	xmm4 += xmm6;		/*12-15	newR2 = R1 + T1 */ \
	xmm6 *= CONST2_M358;		/*13-18	T1 = T1 * (-.223/.623) */ \
	\
	/* IF m4_conflict EQ 1; */ \
	vec2f64 tmp = xmm3;		/*	Save R1 */ \
	/* ELSE; */ \
	/* destr1 = xmm3 */;		/*	Save R1 */ \
	/* ENDIF; */ \
	xmm3 = xmm0;		/*	newR4 = R1 */ \
	xmm0 += xmm5;		/*16-19	newR3 = R1 + T3 */ \
	xmm5 *= CONST2_M358;		/*17-22	T3 = T3 * (-.223/.623) */ \
	xmm3 += xmm7;		/*14-17	newR4 = R1 + T2 */ \
	xmm7 *= CONST2_M358;		/*15-20	T2 = T2 * (-.223/.623) */ \
	xmm0 += xmm6;		/*20-23	newR3 = newR3 + T1 */ \
	xmm6 *= CONST2_P404;		/*21-26	T1 = T1 * (-.901/-.223) */ \
	xmm3 += xmm5;		/*24-27	newR4 = newR4 + T3 */ \
	xmm5 *= CONST2_P404;		/*25-30	T3 = T3 * (-.901/-.223) */ \
	xmm4 += xmm7;		/*22-25	newR2 = newR2 + T2 */ \
	xmm7 *= CONST2_P404;		/*23-28	T2 = T2 * (-.901/-.223) */ \
	xmm3 += xmm6;		/*30-33	newR4 = newR4 + T1 (final R4) */ \
	xmm0 += xmm7;		/*34-37	newR3 = newR3 + T2 (final R3) */ \
	xmm4 += xmm5;		/*36-39	newR2 = newR2 + T3 (final R2) */ \
\
	xmm1 -= m7;		/*28-31	S1 = R2-R7 */ \
	xmm1 *= CONST2_P975;		/*33-38	S1 = S1 * .975 */ \
	xmm2 -= m6;		/*26-29	S2 = R3-R6 */ \
	xmm2 *= CONST2_P975;		/*31-36	S2 = S2 * .975 */ \
	destr4 = xmm3;		/*	Save R4 */ \
	xmm3 = m4;		/*	S3 = R4 */ \
	xmm3 -= m5;		/*32-35	S3 = R4-R5 */ \
	xmm3 *= CONST2_P975;		/*37-42	S3 = S3 * .975 */ \
	xmm5 = CONST2_P445;		/*	(.434/.975) */ \
	xmm5 *= xmm1;		/*41-46	S1 = S1 * (.434/.975), newI3=S1 */ \
	xmm6 = CONST2_P445;		/*	(.434/.975) */ \
	xmm6 *= xmm2;		/*39-44	S2 = S2 * (.434/.975), newI2=S2 */ \
	xmm7 = CONST2_P445;		/*	(.434/.975) */ \
	xmm7 *= xmm3;		/*43-48	S3 = S3 * (.434/.975), newI4=S3 */ \
	xmm1 -= xmm6;		/*45-48	newI3 = newI3 - S2 */ \
	xmm6 *= CONST2_P180;		/*46-51	S2 = S2 * (.782/.434) */ \
	xmm3 += xmm5;		/*47-50	newI4 = newI4 + S1 */ \
	xmm5 *= CONST2_P180;		/*48-53	S1 = S1 * (.782/.434) */ \
	xmm2 += xmm7;		/*49-52	newI2 = newI2 + S3 */ \
	xmm7 *= CONST2_P180;		/*50-55	S3 = S3 * (.782/.434) */ \
	xmm3 -= xmm6;		/*52-55	newI4 = newI4 - S2 (final I4) */ \
	xmm2 += xmm5;		/*54-57	newI2 = newI2 + S1 (final I2) */ \
	xmm1 -= xmm7;		/*56-59	newI3 = newI3 - S3 (final I3) */ \
	/* IF m4_conflict EQ 1; */ \
	xmm6 = tmp;		/*	Reload final R1 */ \
	destr1 = xmm6;		/*	Save final R1 */ \
	/* ENDIF; */ \
}

#define x7r_fft(r1, r2, r3, r4, r5, r6, r7, t1, memr1) { \
	vec2f64 tmp; \
	t1 = r2; \
	r2 -= r7;			/*	R2-R7 */ \
	r7 += t1;			/* T1 = R2+R7 */ \
	t1 = r3; \
	r3 -= r6;			/*	R3-R6 */ \
	r6 += t1;			/* T2 = R3+R6 */ \
	t1 = r4; \
	r4 -= r5;			/*	R4-R5 */ \
	r5 += t1;			/* T3 = R4+R5 */ \
	t1 = r1;			/* R1 */ \
	t1 += r7;			/* R1+T1 */ \
	t1 += r6;			/* R1+T1+T2 */ \
	t1 += r5;			/* R1+T1+T2+T3 (final R1) */ \
	memr1 = t1; \
	r7 *= CONST2_P623;		/* T1 = T1 * .623 */ \
	r6 *= CONST2_P623;		/* T2 = T2 * .623 */ \
	r5 *= CONST2_P623;		/* T3 = T3 * .623 */ \
	tmp = r2; \
	r2 = r1; \
	t1 = r1; \
	r1 += r7;			/* newR2 = R1 + T1 */ \
	r2 += r5;			/* newR3 = R1 + T3 */ \
	t1 += r6;			/* newR4 = R1 + T2 */ \
	r7 *= CONST2_M358;		/* T1 = T1 * (-.223/.623) */ \
	r6 *= CONST2_M358;		/* T2 = T2 * (-.223/.623) */ \
	r5 *= CONST2_M358;		/* T3 = T3 * (-.223/.623) */ \
	r1 += r6;			/* newR2 = newR2 + T2 */ \
	r2 += r7;			/* newR3 = newR3 + T1 */ \
	t1 += r5;			/* newR4 = newR4 + T3 */ \
	r7 *= CONST2_P404;		/* T1 = T1 * (-.901/-.223) */ \
	r6 *= CONST2_P404;		/* T2 = T2 * (-.901/-.223) */ \
	r5 *= CONST2_P404;		/* T3 = T3 * (-.901/-.223) */ \
	r1 += r5;			/* newR2 = newR2 + T3 (final R2) */ \
	r2 += r6;			/* newR3 = newR3 + T2 (final R3) */ \
	t1 += r7;			/* newR4 = newR4 + T1 (final R4) */ \
	r7 = tmp;		/* T1 = R2-R7 */ \
	r7 *= CONST2_P975;		/* T1 = T1 * .975 */ \
	r3 *= CONST2_P975;		/* T2 = T2 * .975 */ \
	r4 *= CONST2_P975;		/* T3 = T3 * .975 */ \
	tmp = r2;		/* final R3 */ \
	r2 = r3;			/* newI2 = T2 */ \
	r6 = r7;			/* newI3 = T1 */ \
	r5 = r4;			/* newI4 = T3 */ \
	r7 *= CONST2_P445;		/* T1 = T1 * (.434/.975) */ \
	r3 *= CONST2_P445;		/* T2 = T2 * (.434/.975) */ \
	r5 *= CONST2_P445;		/* T3 = T3 * (.434/.975) */ \
	r2 += r5;			/* newI2 = newI2 + T3 */ \
	r6 -= r3;			/* newI3 = newI3 - T2 */ \
	r4 += r7;			/* newI4 = newI4 + T1 */ \
	r7 *= CONST2_P180;		/* T1 = T1 * (.782/.434) */ \
	r3 *= CONST2_P180;		/* T2 = T2 * (.782/.434) */ \
	r5 *= CONST2_P180;		/* T3 = T3 * (.782/.434) */ \
	r2 += r7;			/* newI2 = newI2 + T1 (final I2) */ \
	r6 -= r5;			/* newI3 = newI3 - T3 (final I3) */ \
	r4 -= r3;			/* newI4 = newI4 - T2 (final I4) */ \
	r5 = tmp; 		/* final R3 */ \
}

#define xseven_reals_unfft_preload \
\
/* 198.15 clocks */ \

#define x7cl_seven_reals_last_unfft(srcreg,srcinc,d1) \
	xmm0 = xptr(srcreg);		/* R1 */ \
	xmm1 = xptr(srcreg+d1);	/* R2 */ \
	xmm3 = xptr(srcreg+3*d1);	/* R4 */ \
	xmm4 = xptr(srcreg+3*d1+32);	/* R5 */ \
	xmm5 = xptr(srcreg+5*d1);	/* R6 */ \
	x7r_unfft_mem(xmm0, xmm1, xmm2, xmm3, xmm4, xmm5, xmm6, xmm7, xptr(srcreg+d1+32), xptr(srcreg+5*d1+32), xptr(srcreg)); \
	xptr(srcreg+5*d1) = xmm1;	/* Save R3 */ \
	xmm1 = xptr(srcreg+d1+16);	/* R2 */ \
	xptr(srcreg+d1+16) = xmm5;	/* Save R2 */ \
	xptr(srcreg+d1+32) = xmm0;	/* Save R5 */ \
	xmm5 = xptr(srcreg+5*d1+16);	/* R6 */ \
	xptr(srcreg+5*d1+16) = xmm3;	/* Save R4 */ \
	xptr(srcreg+3*d1+32) = xmm7;	/* Save R6 */ \
	xmm4 = xptr(srcreg+3*d1+48);	/* R5 */ \
	xptr(srcreg+3*d1+48) = xmm2;	/* Save R7 */ \
	xmm0 = xptr(srcreg+32);	/* R1 */ \
	xmm3 = xptr(srcreg+3*d1+16);	/* R4 */ \
	x7r_unfft_mem(xmm0, xmm1, xmm2, xmm3, xmm4, xmm5, xmm6, xmm7, xptr(srcreg+d1+48), xptr(srcreg+5*d1+48), xptr(srcreg+d1)); \
	xptr(srcreg+3*d1) = xmm5;	/* Save R2 */ \
	xptr(srcreg+3*d1+16) = xmm1;	/* Save R3 */ \
	xptr(srcreg+32) = xmm3;	/* Save R4 */ \
	xptr(srcreg+d1+48) = xmm0;	/* Save R5 */ \
	xptr(srcreg+5*d1+32) = xmm7;	/* Save R6 */ \
	xptr(srcreg+5*d1+48) = xmm2;	/* Save R7 */ \
\
	xmm0 = xptr(srcreg+16);	/* R1 */ \
	xmm1 = xptr(srcreg+2*d1);	/* R2 */ \
	xmm3 = xptr(srcreg+4*d1);	/* R4 */ \
	xmm4 = xptr(srcreg+4*d1+32);	/* R5 */ \
	xmm5 = xptr(srcreg+6*d1);	/* R6 */ \
	x7r_unfft_mem(xmm0, xmm1, xmm2, xmm3, xmm4, xmm5, xmm6, xmm7, xptr(srcreg+2*d1+32), xptr(srcreg+6*d1+32), xptr(srcreg+16)); \
	xptr(srcreg+6*d1) = xmm1;	/* Save R3 */ \
	xmm1 = xptr(srcreg+2*d1+16);	/* R2 */ \
	xptr(srcreg+2*d1+16) = xmm5;	/* Save R2 */ \
	xptr(srcreg+2*d1+32) = xmm0;	/* Save R5 */ \
	xmm5 = xptr(srcreg+6*d1+16);	/* R6 */ \
	xptr(srcreg+6*d1+16) = xmm3;	/* Save R4 */ \
	xptr(srcreg+4*d1+32) = xmm7;	/* Save R6 */ \
	xmm4 = xptr(srcreg+4*d1+48);	/* R5 */ \
	xptr(srcreg+4*d1+48) = xmm2;	/* Save R7 */ \
	xmm0 = xptr(srcreg+48);	/* R1 */ \
	xmm3 = xptr(srcreg+4*d1+16);	/* R4 */ \
	x7r_unfft_mem(xmm0, xmm1, xmm2, xmm3, xmm4, xmm5, xmm6, xmm7, xptr(srcreg+2*d1+48), xptr(srcreg+6*d1+48), xptr(srcreg+2*d1)); \
	xptr(srcreg+4*d1) = xmm5;	/* Save R2 */ \
	xptr(srcreg+4*d1+16) = xmm1;	/* Save R3 */ \
	xptr(srcreg+48) = xmm3;	/* Save R4 */ \
	xptr(srcreg+2*d1+48) = xmm0;	/* Save R5 */ \
	xptr(srcreg+6*d1+32) = xmm7;	/* Save R6 */ \
	xptr(srcreg+6*d1+48) = xmm2;	/* Save R7 */ \
	srcreg += srcinc; \
\
/* All but r3 and r7 must be pre-loaded */ \

#define x7r_unfft_mem(r1, r2, r3, r4, r5, r6, r7, t1, memr3, memr7, outmemr1) { \
	t1 = r1;			/* R1 */ \
	t1 += r2;			/* R1 + R2 */ \
	t1 += r4;			/* R1 + R2 + R3 */ \
	t1 += r6;			/* R1 + R2 + R3 + R4 (final R1) */ \
	outmemr1 = t1;		/* Save final R1 */ \
\
	r7 = r1;			/* A2 = R1 */ \
	t1 = r1;			/* A3 = R1 */ \
	r2 *= CONST2_P623;		/* S2 = R2 * .623 */ \
	r4 *= CONST2_P623;		/* S3 = R3 * .623 */ \
	r6 *= CONST2_P623;		/* S4 = R4 * .623 */ \
	r7 += r2;			/* A2 = A2 + S2 */ \
	t1 += r6;			/* A3 = A3 + S4 */ \
	r1 += r4;			/* A4 = A4 + S3 */ \
	r2 *= CONST2_M358;		/* S2 = S2 * (-.223/.623) */ \
	r4 *= CONST2_M358;		/* S3 = S3 * (-.223/.623) */ \
	r6 *= CONST2_M358;		/* S4 = S4 * (-.223/.623) */ \
	r7 += r4;			/* A2 = A2 + S3 */ \
	t1 += r2;			/* A3 = A3 + S2 */ \
	r1 += r6;			/* A4 = A4 + S4 */ \
	r2 *= CONST2_P404;		/* S2 = S2 * (-.901/-.223) */ \
	r4 *= CONST2_P404;		/* S3 = S3 * (-.901/-.223) */ \
	r6 *= CONST2_P404;		/* S4 = S4 * (-.901/-.223) */ \
	r7 += r6;			/* A2 = A2 + S4 */ \
	t1 += r4;			/* A3 = A3 + S3 */ \
	r1 += r2;			/* A4 = A4 + S2 */ \
\
	vec2f64 tmp2 = r7;		/* Save A2 */ \
	r3 = memr3;		/* Load I2 */ \
	r7 = memr7;		/* Load I3 */ \
	r3 *= CONST2_P975;		/* T2 = I2*.975 */ \
	r5 *= CONST2_P975;		/* T3 = I3*.975 */ \
	r7 *= CONST2_P975;		/* T4 = I4*.975 */ \
	r6 = r5;			/* B2 = T3 */ \
	r2 = r3;			/* B3 = T2 */ \
	r4 = r7;			/* B4 = T4 */ \
	r3 *= CONST2_P445;		/* T2 = T2 * (.434/.975) */ \
	r5 *= CONST2_P445;		/* T3 = T3 * (.434/.975) */ \
	r7 *= CONST2_P445;		/* T4 = T4 * (.434/.975) */ \
	r6 += r7;			/* B2 = B2 + T4 */ \
	r2 -= r5;			/* B3 = B3 - T3 */ \
	r4 += r3;			/* B4 = B4 + T2 */ \
	r3 *= CONST2_P180;		/* T2 = T2 * (.782/.434) */ \
	r5 *= CONST2_P180;		/* T3 = T3 * (.782/.434) */ \
	r7 *= CONST2_P180;		/* T4 = T4 * (.782/.434) */ \
	r6 += r3;			/* B2 = B2 + T2 */ \
	r3 = tmp2;		/* Reload A2 */ \
	r2 -= r7;			/* B3 = B3 - T4 */ \
	r4 -= r5;			/* B4 = B4 - T3 */ \
\
	r3 -= r6;			/* A2 = A2 - B2 (final R7) */ \
	r6 += tmp2;		/* B2 = A2 + B2 (final R2) */ \
	r5 = t1; \
	t1 -= r2;			/* A3 = A3 - B3 (final R6) */ \
	r2 += r5;			/* B3 = A3 + B3 (final R3) */ \
	r5 = r1; \
	r1 -= r4;			/* A4 = A4 - B4 (final R5) */ \
	r4 += r5;			/* B4 = A4 + B4 (final R4) */ \
}

#define x7r_unfft(r1, r2, r3, r4, r5, r6, r7, t1, memr1) { \
	t1 = r1;			/* R1 */ \
	t1 += r2;			/* R1 + R2 */ \
	t1 += r4;			/* R1 + R2 + R3 */ \
	t1 += r6;			/* R1 + R2 + R3 + R4 (final R1) */ \
	r3 *= CONST2_P975;		/* T2 = I2*.975 */ \
	r5 *= CONST2_P975;		/* T3 = I3*.975 */ \
	r7 *= CONST2_P975;		/* T4 = I4*.975 */ \
	memr1 = t1;		/* Save final R1 */ \
	vec2f64 tmp2 = r2;		/* Save R2 */ \
	vec2f64 tmp3 = r4;		/* Save R3 */ \
	t1 = r5;			/* B2 = T3 */ \
	r2 = r3;			/* B3 = T2 */ \
	r4 = r7;			/* B4 = T4 */ \
	r3 *= CONST2_P445;		/* T2 = T2 * (.434/.975) */ \
	r5 *= CONST2_P445;		/* T3 = T3 * (.434/.975) */ \
	r7 *= CONST2_P445;		/* T4 = T4 * (.434/.975) */ \
	t1 += r7;			/* B2 = B2 + T4 */ \
	r2 -= r5;			/* B3 = B3 - T3 */ \
	r4 += r3;			/* B4 = B4 + T2 */ \
	r3 *= CONST2_P180;		/* T2 = T2 * (.782/.434) */ \
	r5 *= CONST2_P180;		/* T3 = T3 * (.782/.434) */ \
	r7 *= CONST2_P180;		/* T4 = T4 * (.782/.434) */ \
	t1 += r3;			/* B2 = B2 + T2 */ \
	r2 -= r7;			/* B3 = B3 - T4 */ \
	r4 -= r5;			/* B4 = B4 - T3 */ \
	r3 = tmp2;		/* Reload R2 */ \
	r5 = tmp3;		/* Reload R3 */ \
	tmp2 = t1;		/* Save B2 */ \
	r7 = r1;			/* A2 = R1 */ \
	t1 = r1;			/* A3 = R1 */ \
	r3 *= CONST2_P623;		/* S2 = R2 * .623 */ \
	r5 *= CONST2_P623;		/* S3 = R3 * .623 */ \
	r6 *= CONST2_P623;		/* S4 = R4 * .623 */ \
	r7 += r3;			/* A2 = A2 + S2 */ \
	t1 += r6;			/* A3 = A3 + S4 */ \
	r1 += r5;			/* A4 = A4 + S3 */ \
	r3 *= CONST2_M358;		/* S2 = S2 * (-.223/.623) */ \
	r5 *= CONST2_M358;		/* S3 = S3 * (-.223/.623) */ \
	r6 *= CONST2_M358;		/* S4 = S4 * (-.223/.623) */ \
	r7 += r5;			/* A2 = A2 + S3 */ \
	t1 += r3;			/* A3 = A3 + S2 */ \
	r1 += r6;			/* A4 = A4 + S4 */ \
	r3 *= CONST2_P404;		/* S2 = S2 * (-.901/-.223) */ \
	r5 *= CONST2_P404;		/* S3 = S3 * (-.901/-.223) */ \
	r6 *= CONST2_P404;		/* S4 = S4 * (-.901/-.223) */ \
	r7 += r6;			/* A2 = A2 + S4 */ \
	t1 += r5;			/* A3 = A3 + S3 */ \
	r1 += r3;			/* A4 = A4 + S2 */ \
	r3 = tmp2;		/* Reload B2 */ \
	r5 = r7; \
	r7 -= r3;			/* A2 = A2 - B2 (final R7) */ \
	r3 += r5;			/* B2 = A2 + B2 (final R2) */ \
	r5 = t1; \
	t1 -= r2;			/* A3 = A3 - B3 (final R6) */ \
	r2 += r5;			/* B3 = A3 + B3 (final R3) */ \
	r5 = r1; \
	r1 -= r4;			/* A4 = A4 - B4 (final R5) */ \
	r4 += r5;			/* B4 = A4 + B4 (final R4) */ \
}
\
\
\
\
/************************************************************************ */ \
/* These macros process the interior FFT levels. */ \
/************************************************************************ */ \
\
\
/* Macro to operate on 4 64-byte cache lines.  It does 4 two_four_reals_fft, */ \
/* 4 two_two_complex_fft_2, 4 four_complex_fft in a one pass FFT.  The */ \
/* x2cl versions is used in the PFA-6 case of the two pass FFT. */ \

#define x4cl_eight_reals_fft_2(srcreg,srcinc,d1,d2,screg) \
	xmm0 = xptr(srcreg);		/* R1 */ \
	xmm4 = xptr(srcreg+16);	/* R5 */ \
	xmm1 = xptr(srcreg+d1);	/* R2 */ \
	xmm2 = xptr(srcreg+d2);	/* R3 */ \
	xmm6 = xptr(srcreg+d2+16);	/* R7 */ \
	xmm3 = xptr(srcreg+d2+d1);	/* R4 */ \
	x8r2_fft(xmm0, xmm1, xmm2, xmm3, xmm4, xmm5, xmm6, xmm7, xptr(srcreg+d1+16), xptr(srcreg+d2+d1+16), xptr(srcreg)); \
	xmm3 = xptr(srcreg+32);	/* R1 */ \
	xptr(srcreg+16) = xmm2;	/* Save R2 */ \
	xmm2 = xptr(srcreg+48);	/* R5 */ \
	xptr(srcreg+32) = xmm0;	/* Save R3 */ \
	xmm0 = xptr(srcreg+d1+32);	/* R2 */ \
	xptr(srcreg+48) = xmm1;	/* Save R4 */ \
	xmm1 = xptr(srcreg+d1+48);	/* R6 */ \
	xptr(srcreg+d1) = xmm5;	/* Save R5 */ \
	xmm5 = xptr(srcreg+d2+32);	/* R3 */ \
	xptr(srcreg+d1+16) = xmm7;	/* Save R6 */ \
	xmm7 = xptr(srcreg+d2+48);	/* R7 */ \
	xptr(srcreg+d1+32) = xmm4;	/* Save R7 */ \
	xmm4 = xptr(srcreg+d2+d1+32);	/* R4 */ \
	xptr(srcreg+d1+48) = xmm6;	/* Save R8 */ \
	x4c_fft_dst12(xmm3, xmm0, xmm5, xmm4, xmm2, xmm1, xmm7, xmm6, xptr(srcreg+d2+d1+48), screg, 0, srcreg+srcinc+d2, d1, xptr(srcreg+d2+d1+48), xptr(srcreg+d2+d1)); \
	xptr(srcreg+d2) = xmm4;	/* Save R1 */ \
	xptr(srcreg+d2+16) = xmm0;	/* Save R2 */ \
	xptr(srcreg+d2+32) = xmm6;	/* Save R3 */ \
	xptr(srcreg+d2+48) = xmm7;	/* Save R4 */ \
/*	xstore	[srcreg+d2+d1], xmm3	;; Save R5 */ \
	xptr(srcreg+d2+d1+16) = xmm5;	/* Save R6 */ \
	xptr(srcreg+d2+d1+32) = xmm1;	/* Save R7 */ \
/*	xstore	[srcreg+d2+d1+48], xmm2	;; Save R8 */ \
	srcreg += srcinc; \
\

#define g4cl_eight_reals_fft_2(srcreg,srcinc,d1,d2,dstreg,dstinc,e1,e2) \
	xmm0 = xptr(srcreg);		/* R1 */ \
	xmm4 = xptr(srcreg+16);	/* R5 */ \
	xmm1 = xptr(srcreg+d1);	/* R2 */ \
	xmm2 = xptr(srcreg+d2);	/* R3 */ \
	xmm6 = xptr(srcreg+d2+16);	/* R7 */ \
	xmm3 = xptr(srcreg+d2+d1);	/* R4 */ \
	x8r2_fft(xmm0, xmm1, xmm2, xmm3, xmm4, xmm5, xmm6, xmm7, xptr(srcreg+d1+16), xptr(srcreg+d2+d1+16), xptr(dstreg)); \
	xmm3 = xptr(srcreg+32);	/* R1 */ \
	xptr(dstreg+16) = xmm2;	/* Save R2 */ \
	xmm2 = xptr(srcreg+48);	/* R5 */ \
	xptr(dstreg+32) = xmm0;	/* Save R3 */ \
	xmm0 = xptr(srcreg+d1+32);	/* R2 */ \
	xptr(dstreg+48) = xmm1;	/* Save R4 */ \
	xmm1 = xptr(srcreg+d1+48);	/* R6 */ \
	xptr(dstreg+e1) = xmm5;	/* Save R5 */ \
	xmm5 = xptr(srcreg+d2+32);	/* R3 */ \
	xptr(dstreg+e1+16) = xmm7;	/* Save R6 */ \
	xmm7 = xptr(srcreg+d2+48);	/* R7 */ \
	xptr(dstreg+e1+32) = xmm4;	/* Save R7 */ \
	xmm4 = xptr(srcreg+d2+d1+32);	/* R4 */ \
	xptr(dstreg+e1+48) = xmm6;	/* Save R8 */ \
	x4c_fft_dst12(xmm3, xmm0, xmm5, xmm4, xmm2, xmm1, xmm7, xmm6, xptr(srcreg+d2+d1+48), rdi, 0, dstreg+dstinc+e2, e1, xptr(dstreg+e2+e1+48), xptr(dstreg+e2+e1)); \
	srcreg += srcinc; \
	xptr(dstreg+e2) = xmm4;	/* Save R1 */ \
	xptr(dstreg+e2+16) = xmm0;	/* Save R2 */ \
	xptr(dstreg+e2+32) = xmm6;	/* Save R3 */ \
	xptr(dstreg+e2+48) = xmm7;	/* Save R4 */ \
/*	xstore	[dstreg+e2+e1], xmm3	;; Save R5 */ \
	xptr(dstreg+e2+e1+16) = xmm5;	/* Save R6 */ \
	xptr(dstreg+e2+e1+32) = xmm1;	/* Save R7 */ \
/*	xstore	[dstreg+e2+e1+48], xmm2	;; Save R8 */ \
	dstreg += dstinc; \
\

#define x2cl_eight_reals_fft_2(srcreg,srcinc,d1) \
	xmm0 = xptr(srcreg);		/* R1 */ \
	xmm1 = xptr(srcreg+d1);	/* R2 */ \
	xmm2 = xptr(srcreg+16);	/* R3 */ \
	xmm3 = xptr(srcreg+d1+16);	/* R4 */ \
	xmm4 = xptr(srcreg+32);	/* R5 */ \
	xmm6 = xptr(srcreg+48);	/* R7 */ \
	x8r2_fft(xmm0, xmm1, xmm2, xmm3, xmm4, xmm5, xmm6, xmm7, xptr(srcreg+d1+32), xptr(srcreg+d1+48), xptr(srcreg)); \
	xptr(srcreg+16) = xmm2;	/* Save R2 */ \
	xptr(srcreg+32) = xmm0;	/* Save R3 */ \
	xptr(srcreg+48) = xmm1;	/* Save R4 */ \
	xptr(srcreg+d1) = xmm5;	/* Save R5 */ \
	xptr(srcreg+d1+16) = xmm7;	/* Save R6 */ \
	xptr(srcreg+d1+32) = xmm4;	/* Save R7 */ \
	xptr(srcreg+d1+48) = xmm6;	/* Save R8 */ \
	srcreg += srcinc; \
\

#define x8r2_fft(r1, r2, r3, r4, r5, r6, r7, r8, mem6, mem8, dest1) \
	r6 = r1; \
	r1 -= r3;			/* new R3 = R1 - R3 (final R3) */ \
	r3 += r6;			/* new R1 = R1 + R3 */ \
	r8 = r2; \
	r2 -= r4;			/* new R4 = R2 - R4 (final R4) */ \
	r4 += r8;			/* new R2 = R2 + R4 */ \
	r6 = r3; \
	r3 -= r4;			/* R1 = R1 - R2 (final R2) */ \
	r4 += r6;			/* R2 = R1 + R2 (final R1) */ \
	r8 = mem8; \
	r6 = mem6; \
	r6 -= r8;			/* R2 - I2 */ \
	r8 += mem6;		/* R2 + I2 */ \
	r6 *= XMM_SQRTHALF2;	/* newR2 */ \
	dest1 = r4; \
	r8 *= XMM_SQRTHALF2;	/* newI2 */ \
	r4 = r7; \
	r7 -= r8;			/* I1 = I1 - I2 (new I2) */ \
	r8 += r4;			/* I2 = I1 + I2 (new I1) */ \
	r4 = r5; \
	r5 -= r6;			/* R1 = R1 - R2 (new R2) */ \
	r6 += r4;			/* R2 = R1 + R2 (new R1) */ \
\
\
\
/* Macro to operate on 4 64-byte cache lines.  It does 4 two_four_reals_unfft, */ \
/* 4 two_two_complex_unfft_2, 4 four_complex_unfft in a one pass FFT. */ \

#define x4cl_eight_reals_unfft_2(srcreg,srcinc,d1,d2) \
	xmm0 = xptr(srcreg);		/* R1 */ \
	xmm1 = xptr(srcreg+32);	/* R2 */ \
	xmm2 = xptr(srcreg+d1);	/* R3 */ \
	xmm3 = xptr(srcreg+d1+32);	/* R4 */ \
	xmm4 = xptr(srcreg+d2);	/* R5 */ \
	xmm5 = xptr(srcreg+d2+32);	/* R6 */ \
	xmm6 = xptr(srcreg+d2+d1);	/* R7 */ \
	xprefetchw(u8ptr(srcreg+srcinc)); \
	xprefetchw(u8ptr(srcreg+srcinc+d1)); \
	x8r2_unfft(xmm0, xmm1, xmm2, xmm3, xmm4, xmm5, xmm6, xmm7, xptr(srcreg+d2+d1+32), xptr(srcreg)); \
	xmm2 = xptr(srcreg+16);	/* R1 */ \
	xptr(srcreg+d2) = xmm3;	/* Save R2 */ \
	xmm3 = xptr(srcreg+d2+16);	/* R5 */ \
	xptr(srcreg+16) = xmm1;	/* Save R3 */ \
	xmm1 = xptr(srcreg+48);	/* R2 */ \
	xptr(srcreg+d2+16) = xmm0;	/* Save R4 */ \
	xmm0 = xptr(srcreg+d2+48);	/* R6 */ \
	xptr(srcreg+32) = xmm6;	/* Save R5 */ \
	xmm6 = xptr(srcreg+d1+16);	/* R3 */ \
	xptr(srcreg+d2+32) = xmm4;	/* Save R6 */ \
	xmm4 = xptr(srcreg+d1+48);	/* R4 */ \
	xptr(srcreg+48) = xmm7;	/* Save R7 */ \
	xmm7 = xptr(srcreg+d2+d1+16);	/* R7 */ \
	xptr(srcreg+d2+48) = xmm5;	/* Save R8 */ \
	xprefetchw(u8ptr(srcreg+srcinc+d2)); \
	xprefetchw(u8ptr(srcreg+srcinc+d2+d1)); \
	x8r2_unfft(xmm2, xmm1, xmm6, xmm4, xmm3, xmm0, xmm7, xmm5, xptr(srcreg+d2+d1+48), xptr(srcreg+d1)); \
	xptr(srcreg+d2+d1) = xmm4;	/* Save R2 */ \
	xptr(srcreg+d1+16) = xmm1;	/* Save R3 */ \
	xptr(srcreg+d2+d1+16) = xmm2;	/* Save R4 */ \
	xptr(srcreg+d1+32) = xmm7;	/* Save R5 */ \
	xptr(srcreg+d2+d1+32) = xmm3;	/* Save R6 */ \
	xptr(srcreg+d1+48) = xmm5;	/* Save R7 */ \
	xptr(srcreg+d2+d1+48) = xmm0;	/* Save R8 */ \
	srcreg += srcinc; \
\

#define x8r2_unfft(r1, r2, r3, r4, r5, r6, r7, r8, mem8, dest1) \
	r8 = r1; \
	r1 -= r2;			/* new R2 = R1 - R2 */ \
	r2 += r8;			/* new R1 = R1 + R2 */ \
	r8 = r1; \
	r1 -= r4;			/* R2 = R2 - R4 (final R4) */ \
	r4 += r8;			/* R4 = R2 + R4 (final R2) */ \
	r8 = r2; \
	r2 -= r3;			/* R1 = R1 - R3 (final R3) */ \
	r3 += r8;			/* R3 = R1 + R3 (final R1) */ \
	r8 = mem8; \
	r8 += r6;			/* new I1 = I1 + I2 */ \
	dest1 = r3; \
	r6 -= mem8;		/* new I2 = I1 - I2 */ \
	r3 = r5; \
	r5 -= r7;			/* new R2 = R1 - R2 */ \
	r7 += r3;			/* new R1 = R1 + R2 */ \
	r3 = XMM_SQRTHALF2; \
	r6 *= r3;			/* B2 = I2 * sine */ \
	r5 *= r3;			/* A2 = R2 * sine */ \
	r3 = r6;			/* Save B2 (C2 = B2) */ \
	r6 -= r5;			/* C2 = C2 - A2 (new I2) */ \
	r5 += r3;			/* A2 = A2 + B2 (new R2) */ \
\
\
\

#define x1cl_half_eight_reals_fft_2(srcreg,srcinc) \
	xmm0 = xptr(srcreg);		/* R1 */ \
	xmm1 = xptr(srcreg+16);	/* R2 */ \
	xmm2 = xptr(srcreg+32);	/* R3 */ \
	xmm3 = xptr(srcreg+48);	/* R4 */ \
	x8r2_half_fft(xmm0, xmm1, xmm2, xmm3, xmm4, xmm5, xmm6); \
	xptr(srcreg) = xmm3;		/* Save R1 */ \
	xptr(srcreg+16) = xmm2;	/* Save R2 */ \
	xptr(srcreg+32) = xmm0;	/* Save R3 */ \
	xptr(srcreg+48) = xmm1;	/* Save R4 */ \
	srcreg += srcinc; \
\

#define x8r2_half_fft(r1, r2, r3, r4, t1, t2, t3) \
	t1 = r1; \
	r1 -= r3;			/* new R3 = R1 - R3 (final R3) */ \
	r3 += t1;			/* new R1 = R1 + R3 */ \
	t2 = r2; \
	r2 -= r4;			/* new R4 = R2 - R4 (final R4) */ \
	r4 += t2;			/* new R2 = R2 + R4 */ \
	t3 = r3; \
	r3 -= r4;			/* R1 = R1 - R2 (final R2) */ \
	r4 += t3;			/* R2 = R1 + R2 (final R1) */ \
\
\

#define x2cl_half_eight_reals_unfft_2(srcreg,srcinc,d1) \
	xmm0 = xptr(srcreg);		/* R1 */ \
	xmm1 = xptr(srcreg+32);	/* R2 */ \
	xmm2 = xptr(srcreg+d1);	/* R3 */ \
	xmm3 = xptr(srcreg+d1+32);	/* R4 */ \
	x8r2_half_unfft(xmm0, xmm1, xmm2, xmm3, xmm4, xmm5, xmm6); \
	xmm4 = xptr(srcreg+16);	/* R1 */ \
	xmm5 = xptr(srcreg+48);	/* R2 */ \
	xmm6 = xptr(srcreg+d1+16);	/* R3 */ \
	xmm7 = xptr(srcreg+d1+48);	/* R4 */ \
	xptr(srcreg) = xmm2;		/* Save R1 */ \
	xptr(srcreg+16) = xmm3;	/* Save R2 */ \
	xptr(srcreg+32) = xmm1;	/* Save R3 */ \
	xptr(srcreg+48) = xmm0;	/* Save R4 */ \
	x8r2_half_unfft(xmm4, xmm5, xmm6, xmm7, xmm2, xmm3, xmm1); \
	xptr(srcreg+d1) = xmm6;	/* Save R1 */ \
	xptr(srcreg+d1+16) = xmm7;	/* Save R2 */ \
	xptr(srcreg+d1+32) = xmm5;	/* Save R3 */ \
	xptr(srcreg+d1+48) = xmm4;	/* Save R4 */ \
	srcreg += srcinc; \
\

#define x1cl_half_eight_reals_unfft_2(srcreg,srcinc) \
	xmm0 = xptr(srcreg);		/* R1 */ \
	xmm1 = xptr(srcreg+16);	/* R2 */ \
	xmm2 = xptr(srcreg+32);	/* R3 */ \
	xmm3 = xptr(srcreg+48);	/* R4 */ \
	x8r2_half_unfft(xmm0, xmm1, xmm2, xmm3, xmm4, xmm5, xmm6); \
	xptr(srcreg) = xmm2;		/* Save R1 */ \
	xptr(srcreg+16) = xmm3;	/* Save R2 */ \
	xptr(srcreg+32) = xmm1;	/* Save R3 */ \
	xptr(srcreg+48) = xmm0;	/* Save R4 */ \
	srcreg += srcinc; \
\

#define x8r2_half_unfft(r1, r2, r3, r4, t1, t2, t3) \
	t1 = r1; \
	r1 -= r2;			/* new R2 = R1 - R2 */ \
	r2 += t1;			/* new R1 = R1 + R2 */ \
	t2 = r1; \
	r1 -= r4;			/* R2 = R2 - R4 (final R4) */ \
	r4 += t2;			/* R4 = R2 + R4 (final R2) */ \
	t3 = r2; \
	r2 -= r3;			/* R1 = R1 - R3 (final R3) */ \
	r3 += t3;			/* R3 = R1 + R3 (final R1) */ \
\
\
/* Macro to operate on 2 64-byte cache lines.  It does 2 two_four_reals_unfft, */ \
/* 2 two_two_complex_unfft_2, 2 four_complex_unfft in a one pass FFT. */ \

#define x2cl_eight_reals_unfft_2(srcreg,srcinc,d1) \
	xmm0 = xptr(srcreg);		/* R1 */ \
	xmm1 = xptr(srcreg+16);	/* R2 */ \
	xmm2 = xptr(srcreg+32);	/* R3 */ \
	xmm3 = xptr(srcreg+48);	/* R4 */ \
	xmm4 = xptr(srcreg+d1);	/* R5 */ \
	xmm5 = xptr(srcreg+d1+16);	/* R6 */ \
	xmm6 = xptr(srcreg+d1+32);	/* R7 */ \
	x8r2_unfft(xmm0, xmm1, xmm2, xmm3, xmm4, xmm5, xmm6, xmm7, xptr(srcreg+d1+48), xptr(srcreg)); \
	xptr(srcreg+d1) = xmm3;	/* Save R2 */ \
	xptr(srcreg+16) = xmm1;	/* Save R3 */ \
	xptr(srcreg+d1+16) = xmm0;	/* Save R4 */ \
	xptr(srcreg+32) = xmm6;	/* Save R5 */ \
	xptr(srcreg+d1+32) = xmm4;	/* Save R6 */ \
	xptr(srcreg+48) = xmm7;	/* Save R7 */ \
	xptr(srcreg+d1+48) = xmm5;	/* Save R8 */ \
	srcreg += srcinc; \
\

#define g2cl_eight_reals_unfft_2(srcreg,srcinc,d1,dstreg,dstinc,e1) \
	xmm0 = xptr(srcreg);		/* R1 */ \
	xmm1 = xptr(srcreg+16);	/* R2 */ \
	xmm2 = xptr(srcreg+32);	/* R3 */ \
	xmm3 = xptr(srcreg+48);	/* R4 */ \
	xmm4 = xptr(srcreg+d1);	/* R5 */ \
	xmm5 = xptr(srcreg+d1+16);	/* R6 */ \
	xmm6 = xptr(srcreg+d1+32);	/* R7 */ \
	x8r2_unfft(xmm0, xmm1, xmm2, xmm3, xmm4, xmm5, xmm6, xmm7, xptr(srcreg+d1+48), xptr(dstreg)); \
	srcreg += srcinc; \
	xptr(dstreg+e1) = xmm3;	/* Save R2 */ \
	xptr(dstreg+16) = xmm1;	/* Save R3 */ \
	xptr(dstreg+e1+16) = xmm0;	/* Save R4 */ \
	xptr(dstreg+32) = xmm6;	/* Save R5 */ \
	xptr(dstreg+e1+32) = xmm4;	/* Save R6 */ \
	xptr(dstreg+48) = xmm7;	/* Save R7 */ \
	xptr(dstreg+e1+48) = xmm5;	/* Save R8 */ \
	dstreg += dstinc; \
\
\
\
/* Macro to operate on 4 64-byte cache lines.  It does 4 four_complex_ffts */ \
/* in a one pass FFT. */ \
/* 98.35 clocks */ \

#define x4cl_four_complex_fft(srcreg,srcinc,d1,d2,screg) \
	x4c_fft_mem_dst2(xptr(srcreg),xptr(srcreg+d1),xptr(srcreg+d2),xptr(srcreg+d2+d1),xptr(srcreg+16),xptr(srcreg+d1+16),xptr(srcreg+d2+16),xptr(srcreg+d2+d1+16),screg,0,srcreg+srcinc,d1,xptr(srcreg+d1)); \
	xptr(srcreg) = xmm5;		/* Save R1 */ \
	xmm5 = xptr(srcreg+32);	/* R1 */ \
	xptr(srcreg+16) = xmm7;	/* Save R2 */ \
	xmm7 = xptr(srcreg+48);	/* R5 */ \
	xptr(srcreg+32) = xmm1;	/* Save R3 */ \
	xmm1 = xptr(srcreg+d1+32);	/* R2 */ \
	xptr(srcreg+48) = xmm3;	/* Save R4 */ \
	xmm3 = xptr(srcreg+d1+48);	/* R6 */ \
/*	xstore	[srcreg+d1], xmm2	;; Save R5 */ \
	xmm2 = xptr(srcreg+d2+32);	/* R3 */ \
	xptr(srcreg+d1+16) = xmm4;	/* Save R6 */ \
	xmm4 = xptr(srcreg+d2+48);	/* R7 */ \
	xptr(srcreg+d1+32) = xmm6;	/* Save R7 */ \
	xmm6 = xptr(srcreg+d2+d1+32);	/* R4 */ \
	xptr(srcreg+d1+48) = xmm0;	/* Save R8 */ \
	x4c_fft_dst12(xmm5, xmm1, xmm2, xmm6, xmm7, xmm3, xmm4, xmm0, xptr(srcreg+d2+d1+48), screg, XMM_SCD, srcreg+srcinc+d2, d1, xptr(srcreg+d2+d1+48), xptr(srcreg+d2+d1)); \
	xptr(srcreg+d2) = xmm6;	/* Save R1 */ \
	xptr(srcreg+d2+16) = xmm1;	/* Save R2 */ \
	xptr(srcreg+d2+32) = xmm0;	/* Save R3 */ \
	xptr(srcreg+d2+48) = xmm4;	/* Save R4 */ \
/*	xstore	[srcreg+d2+d1], xmm5	;; Save R5 */ \
	xptr(srcreg+d2+d1+16) = xmm2;	/* Save R6 */ \
	xptr(srcreg+d2+d1+32) = xmm3;	/* Save R7 */ \
/*	xstore	[srcreg+d2+d1+48], xmm7	;; Save R8 */ \
	srcreg += srcinc; \
\

#define g4cl_four_complex_fft(srcreg,srcinc,d1,d2,dstreg,dstinc,e1,e2) \
	xprefetchw(u8ptr(srcreg+srcinc)); \
	x4c_fft_mem_dst12(xptr(srcreg),xptr(srcreg+d1),xptr(srcreg+d2),xptr(srcreg+d2+d1),xptr(srcreg+16),xptr(srcreg+d1+16),xptr(srcreg+d2+16),xptr(srcreg+d2+d1+16),rdi,0,dstreg+dstinc,e1,xptr(dstreg+e1+48),xptr(dstreg+e1)); \
	xprefetchw(u8ptr(srcreg+srcinc+d1)); \
	xptr(dstreg) = xmm5;		/* Save R1 */ \
	xptr(dstreg+16) = xmm7;	/* Save R2 */ \
	xptr(dstreg+32) = xmm1;	/* Save R3 */ \
	xptr(dstreg+48) = xmm3;	/* Save R4 */ \
/*	xstore	[dstreg+e1], xmm2	;; Save R5 */ \
	xptr(dstreg+e1+16) = xmm4;	/* Save R6 */ \
	xptr(dstreg+e1+32) = xmm6;	/* Save R7 */ \
/*	xstore	[dstreg+e1+48], xmm0	;; Save R8 */ \
	xprefetchw(u8ptr(srcreg+srcinc+d2)); \
	x4c_fft_mem_dst12(xptr(srcreg+32),xptr(srcreg+d1+32),xptr(srcreg+d2+32),xptr(srcreg+d2+d1+32),xptr(srcreg+48),xptr(srcreg+d1+48),xptr(srcreg+d2+48),xptr(srcreg+d2+d1+48),rdi,XMM_SCD,dstreg+dstinc+e2,e1,xptr(dstreg+e2+e1+48),xptr(dstreg+e2+e1)); \
	xprefetchw(u8ptr(srcreg+srcinc+d2+d1)); \
	srcreg += srcinc; \
	xptr(dstreg+e2) = xmm5;	/* Save R1 */ \
	xptr(dstreg+e2+16) = xmm7;	/* Save R2 */ \
	xptr(dstreg+e2+32) = xmm1;	/* Save R3 */ \
	xptr(dstreg+e2+48) = xmm3;	/* Save R4 */ \
/*	xstore	[dstreg+e2+e1], xmm2	;; Save R5 */ \
	xptr(dstreg+e2+e1+16) = xmm4;	/* Save R6 */ \
	xptr(dstreg+e2+e1+32) = xmm6;	/* Save R7 */ \
/*	xstore	[dstreg+e2+e1+48], xmm0	;; Save R8 */ \
	dstreg += dstinc; \
\
/* 50.35 clocks */ \

#define x2cl_four_complex_fft(srcreg,srcinc,d1) \
	x4c_fft_mem_dst12(xptr(srcreg), xptr(srcreg+d1), xptr(srcreg+16), xptr(srcreg+d1+16), xptr(srcreg+32), xptr(srcreg+d1+32), xptr(srcreg+48), xptr(srcreg+d1+48), rdi,0,srcreg+srcinc,d1,xptr(srcreg+d1+48),xptr(srcreg+d1)); \
	xptr(srcreg) = xmm5;		/* Save R1 */ \
	xptr(srcreg+16) = xmm7;	/* Save R2 */ \
	xptr(srcreg+32) = xmm1;	/* Save R3 */ \
	xptr(srcreg+48) = xmm3;	/* Save R4 */ \
/*	xstore	[srcreg+d1], xmm2	;; Save R5 */ \
	xptr(srcreg+d1+16) = xmm4;	/* Save R6 */ \
	xptr(srcreg+d1+32) = xmm6;	/* Save R7 */ \
/*	xstore	[srcreg+d1+48], xmm0	;; Save R8 */ \
	srcreg += srcinc; \
\

#if 1
#define x4c_fft_dst12(r1, r2, r3, r4, r5, r6, r7, r8, mem8, screg, off, pre1, pre2, dst1, dst2) \
	x4c_fft_PART1(r1, r2, r3, r4, r5, r6, r7, r8, mem8, screg, off, pre1, pre2); \
	dst1 = r5; \
	dst2 = r1;
#define x4c_fft_dst1(r1, r2, r3, r4, r5, r6, r7, r8, mem8, screg, off, pre1, pre2, dst1) \
	x4c_fft_PART1(r1, r2, r3, r4, r5, r6, r7, r8, mem8, screg, off, pre1, pre2); \
	dst1 = r5;
#define x4c_fft_dst2(r1, r2, r3, r4, r5, r6, r7, r8, mem8, screg, off, pre1, pre2, dst2) \
	x4c_fft_PART1(r1, r2, r3, r4, r5, r6, r7, r8, mem8, screg, off, pre1, pre2); \
	dst2 = r1;
#define x4c_fft_nodst(r1, r2, r3, r4, r5, r6, r7, r8, mem8, screg, off, pre1, pre2) \
	x4c_fft_PART1(r1, r2, r3, r4, r5, r6, r7, r8, mem8, screg, off, pre1, pre2); \

#define x4c_fft_PART1(r1, r2, r3, r4, r5, r6, r7, r8, mem8, screg, off, pre1, pre2) { \
	vec2f64 temp; \
	temp = xptr(screg+off+32+16);	/* cosine/sine */ \
	r8 = temp * r3;		/* A3 = R3 * cosine/sine	;1-6 */ \
	r8 -= r7;			/* A3 = A3 - I3			;8-11 */ \
	r7 *= temp;			/* B3 = I3 * cosine/sine	;3-8 */ \
	r7 += r3;			/* B3 = B3 + R3			;10-13 */ \
	temp = xptr(screg+off+0+16);	/* cosine/sine */ \
	r3 = temp * r2;			/* A2 = R2 * cosine/sine	;5-10 */ \
	r3 -= r6;			/* A2 = A2 - I2			;12-15 */ \
	r6 *= temp;	/* B2 = I2 * cosine/sine	;9-14 */ \
	r6 += r2;			/* B2 = B2 + R2			;16-19 */ \
	temp = xptr(screg+off+64+16);	/* cosine/sine */ \
	r2 = temp * mem8;		/* B4 = I4 * cosine/sine	;11-16 */ \
	r2 += r4;			/* B4 = B4 + R4			;18-21 */ \
	r4 *= temp;			/* A4 = R4 * cosine/sine	;7-12 */ \
	r4 -= mem8;			/* A4 = A4 - I4			;14-17 */ \
	temp = xptr(screg+off+32); \
	r8 *= temp;					/* A3 = A3 * sine (new R3)	;13-18 */ \
	r7 *= temp;					/* B3 = B3 * sine (new I3)	;15-20 */ \
	r4 *= xptr(screg+off+64);	/* A4 = A4 * sine (new R4)	;19-24 */ \
	temp = xptr(screg+off); \
	r3 *= temp;					/* A2 = A2 * sine (new R2)	;17-22 */ \
	r6 *= temp;					/* B2 = B2 * sine (new I2)	;21-26 */ \
	xprefetchw(u8ptr(pre1)); \
\
	temp = r1; \
	r1 -= r8;			/* R1 = R1 - R3 (new R3)	;20-23 */ \
	r8 += temp;			/* R3 = R1 + R3 (new R1)	;24-27 */ \
	temp = r5; \
	r5 -= r7;			/* I1 = I1 - I3 (new I3)	;22-25 */ \
	r7 += temp;			/* I3 = I1 + I3 (new I1)	;28-31 */ \
\
	r2 *= xptr(screg+off+64);	/* B4 = B4 * sine (new I4)	;23-28 */ \
	xprefetchw(u8ptr(pre1+pre2)); \
	temp = r3; \
	r3 -= r4;			/* R2 = R2 - R4 (new R4)	;26-29 */ \
	r4 += temp;			/* R4 = R2 + R4 (new R2)	;34-37 */ \
\
	temp = r6; \
	r6 -= r2;			/* I2 = I2 - I4 (new I4)	;30-33 */ \
	r2 += temp;			/* I4 = I2 + I4 (new I2)	;36-39 */ \
\
	temp = r5; \
	r5 -= r3;			/* I3 = I3 - R4 (final I4)	;32-35 */ \
	r3 += temp;			/* R4 = I3 + R4 (final I3)	;44-47 */ \
	temp = r1; \
	r1 -= r6;			/* R3 = R3 - I4 (final R3)	;38-41 */ \
	r6 += temp;			/* I4 = R3 + I4 (final R4)	;46-49 */ \
	temp = r8; \
	r8 -= r4;			/* R1 = R1 - R2 (final R2)	;40-43 */ \
	r4 += temp;			/* R2 = R1 + R2 (final R1)	;48-51 */ \
	temp = r7; \
	r7 -= r2;			/* I1 = I1 - I2 (final I2)	;42-45 */ \
	r2 += temp;			/* I2 = I1 + I2 (final I1)	;50-53 */ \
}

#else
#define x4c_fft_dst12(r1, r2, r3, r4, r5, r6, r7, r8, mem8, screg, off, pre1, pre2, dst1, dst2) \
	x4c_fft__PART1(r1, r2, r3, r4, r5, r6, r7, r8, mem8, screg, off, pre1, pre2) \
	x4c_fft_PART2_dst1(r1, r2, r3, r4, r5, r6, r7, r8, mem8, screg, off, dst1, exec, dst2)
#define x4c_fft_dst1(r1, r2, r3, r4, r5, r6, r7, r8, mem8, screg, off, pre1, pre2, dst1) \
	x4c_fft__PART1(r1, r2, r3, r4, r5, r6, r7, r8, mem8, screg, off, pre1, pre2) \
	x4c_fft_PART2_dst1(r1, r2, r3, r4, r5, r6, r7, r8, mem8, screg, off, dst1, noexec, 0)
#define x4c_fft_dst2(r1, r2, r3, r4, r5, r6, r7, r8, mem8, screg, off, pre1, pre2, dst2) \
	x4c_fft__PART1(r1, r2, r3, r4, r5, r6, r7, r8, mem8, screg, off, pre1, pre2) \
	x4c_fft_PART2_nodst1(r1, r2, r3, r4, r5, r6, r7, r8, mem8, screg, off, exec, dst2)
#define x4c_fft_nodst(r1, r2, r3, r4, r5, r6, r7, r8, mem8, screg, off, pre1, pre2) \
	x4c_fft__PART1(r1, r2, r3, r4, r5, r6, r7, r8, mem8, screg, off, pre1, pre2) \
	x4c_fft_PART2_nodst1(r1, r2, r3, r4, r5, r6, r7, r8, mem8, screg, off, noexec, 0)

#define x4c_fft__PART1(r1, r2, r3, r4, r5, r6, r7, r8, mem8, screg, off, pre1, pre2) \
	r8 = xptr(screg+off+32+16);	/* cosine/sine */ \
	r8 *= r3;			/* A3 = R3 * cosine/sine	;1-6 */ \
	r8 -= r7;			/* A3 = A3 - I3			;8-11 */ \
	r7 *= xptr(screg+off+32+16);	/* B3 = I3 * cosine/sine	;3-8 */ \
	r7 += r3;			/* B3 = B3 + R3			;10-13 */ \
	r3 = xptr(screg+off+0+16);	/* cosine/sine */ \
	r3 *= r2;			/* A2 = R2 * cosine/sine	;5-10 */ \
	r3 -= r6;			/* A2 = A2 - I2			;12-15 */ \
	r6 *= xptr(screg+off+0+16);	/* B2 = I2 * cosine/sine	;9-14 */ \
	r6 += r2;			/* B2 = B2 + R2			;16-19 */ \
	r2 = xptr(screg+off+64+16);	/* cosine/sine */ \
	r2 *= mem8;		/* B4 = I4 * cosine/sine	;11-16 */ \
	r2 += r4;			/* B4 = B4 + R4			;18-21 */ \
	r4 *= xptr(screg+off+64+16);	/* A4 = R4 * cosine/sine	;7-12 */ \
	r4 -= mem8;		/* A4 = A4 - I4			;14-17 */ \
	r8 *= xptr(screg+off+32);	/* A3 = A3 * sine (new R3)	;13-18 */ \
	r7 *= xptr(screg+off+32);	/* B3 = B3 * sine (new I3)	;15-20 */ \
	r3 *= xptr(screg+off+0);	/* A2 = A2 * sine (new R2)	;17-22 */ \
	r4 *= xptr(screg+off+64);	/* A4 = A4 * sine (new R4)	;19-24 */ \
	xprefetchw(u8ptr(pre1)); \
	 r1 -= r8;			/* R1 = R1 - R3 (new R3)	;20-23 */ \
	 r8 += r8; /* multwo */ \
	r6 *= xptr(screg+off+0);	/* B2 = B2 * sine (new I2)	;21-26 */ \
	 r5 -= r7;			/* I1 = I1 - I3 (new I3)	;22-25 */ \
	 r7 += r7; /* multwo */ \
	r2 *= xptr(screg+off+64);	/* B4 = B4 * sine (new I4)	;23-28 */ \
	xprefetchw(u8ptr(pre1+pre2)); \
	 r8 += r1;			/* R3 = R1 + R3 (new R1)	;24-27 */ \
	 r3 -= r4;			/* R2 = R2 - R4 (new R4)	;26-29 */ \
	 r4 += r4; /* multwo */			/* R4 = R4 * 2			;27-32 */ \
	 r7 += r5;			/* I3 = I1 + I3 (new I1)	;28-31 */ \
	 r6 -= r2;			/* I2 = I2 - I4 (new I4)	;30-33 */ \
	 r2 += r2; /* multwo */			/* I4 = I4 * 2			;31-36 */ \
	r5 -= r3;			/* I3 = I3 - R4 (final I4)	;32-35 */
#define x4c_fft_PART2_dst1(r1, r2, r3, r4, r5, r6, r7, r8, mem8, screg, off, dst1, dst2w, dst2) \
	dst1 = r5; \
	 r4 += r3;			/* R4 = R2 + R4 (new R2)	;34-37 */ \
	r3 += r3; /* multwo */			/* R4 = R4 * 2			;35-40 */ \
	 r2 += r6;			/* I4 = I2 + I4 (new I2)	;36-39 */ \
	r3 += r5;			/* R4 = I3 + R4 (final I3)	;44-47 */ \
	r5 = r1; \
	r1 -= r6;			/* R3 = R3 - I4 (final R3)	;38-41 */ \
dst2w(dst2 = r1); \
	r6 += r5;			/* I4 = R3 + I4 (final R4)	;46-49 */ \
	r5 = r8; \
	r8 -= r4;			/* R1 = R1 - R2 (final R2)	;40-43 */ \
	r4 += r5;			/* R2 = R1 + R2 (final R1)	;48-51 */ \
	r5 = r7; \
	r7 -= r2;			/* I1 = I1 - I2 (final I2)	;42-45 */ \
	r2 += r5;			/* I2 = I1 + I2 (final I1)	;50-53 */
#define x4c_fft_PART2_nodst1(r1, r2, r3, r4, r5, r6, r7, r8, mem8, screg, off, dst2w, dst2) \
	 r4 += r3;			/* R4 = R2 + R4 (new R2)	;34-37 */ \
	r3 += r3; /* multwo */			/* R4 = R4 * 2			;35-40 */ \
	 r2 += r6;			/* I4 = I2 + I4 (new I2)	;36-39 */ \
	r1 -= r6;			/* R3 = R3 - I4 (final R3)	;38-41 */ \
dst2w(dst2 = r1); \
	r6 += r6; /* multwo */			/* I4 = I4 * 2			;39-44 */ \
	r8 -= r4;			/* R1 = R1 - R2 (final R2)	;40-43 */ \
	r4 += r4; /* multwo */			/* R2 = R2 * 2			;41-46 */ \
	r7 -= r2;			/* I1 = I1 - I2 (final I2)	;42-45 */ \
	r2 += r2; /* multwo */			/* I2 = I2 * 2			;43-48 */ \
	r3 += r5;			/* R4 = I3 + R4 (final I3)	;44-47 */ \
	r6 += r1;			/* I4 = R3 + I4 (final R4)	;46-49 */ \
	r4 += r8;			/* R2 = R1 + R2 (final R1)	;48-51 */ \
	r2 += r7;			/* I2 = I1 + I2 (final I1)	;50-53 */
#endif
#if 1
#define x4c_fft_mem_dst12(R1,R2,R3,R4,R5,R6,R7,R8,screg,off,pre1,pre2,dst1,dst2) \
	x4c_fft_mem_PART1(R1,R2,R3,R4,R5,R6,R7,R8,screg,off,pre1,pre2) \
	dst1 = xmm0; \
	dst2 = xmm2; \

#define x4c_fft_mem_dst1(R1,R2,R3,R4,R5,R6,R7,R8,screg,off,pre1,pre2,dst1) \
	x4c_fft_mem_PART1(R1,R2,R3,R4,R5,R6,R7,R8,screg,off,pre1,pre2) \
	dst1 = xmm0; \

#define x4c_fft_mem_dst2(R1,R2,R3,R4,R5,R6,R7,R8,screg,off,pre1,pre2,dst2) \
	x4c_fft_mem_PART1(R1,R2,R3,R4,R5,R6,R7,R8,screg,off,pre1,pre2) \
	dst2 = xmm2; \

#define x4c_fft_mem_nodst(R1,R2,R3,R4,R5,R6,R7,R8,screg,off,pre1,pre2) \
	x4c_fft_mem_PART1(R1,R2,R3,R4,R5,R6,R7,R8,screg,off,pre1,pre2) \

#define x4c_fft_mem_PART1(R1,R2,R3,R4,R5,R6,R7,R8,screg,off,pre1,pre2) { \
	vec2f64 temp2; \
	xmm0 = R3;		/* R3 */ \
	xmm1 = xptr(screg+off+32+16);	/* cosine/sine */ \
	xmm3 = xmm1; \
	xmm1 *= xmm0;		/* A3 = R3 * cosine/sine	;1-6 */ \
	xmm2 = R7;		/* I3 */ \
	xmm3 *= xmm2;		/* B3 = I3 * cosine/sine	;3-8 */ \
	xmm4 = R2;		/* R2 */ \
	xmm6 = xptr(screg+off+0+16);	/* cosine/sine */ \
	xmm4 *= xmm6;		/* A2 = R2 * cosine/sine	;5-10 */ \
	xmm5 = R4;		/* R4 */ \
	xmm7 = xptr(screg+off+64+16);	/* cosine/sine */ \
	xmm5 *= xmm7;		/* A4 = R4 * cosine/sine	;7-12 */ \
	xmm1 -= xmm2;		/* A3 = A3 - I3			;8-11 */ \
	xmm2 = R6;		/* I2 */ \
	xmm6 *= xmm2;		/* B2 = I2 * cosine/sine	;9-14 */ \
	xmm3 += xmm0;		/* B3 = B3 + R3			;10-13 */ \
	xmm0 = R8;		/* I4 */ \
	xmm7 *= xmm0;		/* B4 = I4 * cosine/sine	;11-16 */ \
	xmm4 -= xmm2;		/* A2 = A2 - I2			;12-15 */ \
	xmm2 = xptr(screg+off+32);	/* sine */ \
	xmm1 *= xmm2;		/* A3 = A3 * sine (new R3)	;13-18 */ \
	xmm5 -= xmm0;		/* A4 = A4 - I4			;14-17 */ \
	xmm3 *= xmm2;		/* B3 = B3 * sine (new I3)	;15-20 */ \
	xmm6 += R2;		/* B2 = B2 + R2			;16-19 */ \
	xmm0 = xptr(screg+off+0);	/* sine */ \
	xmm4 *= xmm0;		/* A2 = A2 * sine (new R2)	;17-22 */ \
	xprefetchw(u8ptr(pre1)); \
	xmm7 += R4;		/* B4 = B4 + R4			;18-21 */ \
	xmm5 *= xptr(screg+off+64);	/* A4 = A4 * sine (new R4)	;19-24 */ \
	 xmm2 = R1;		/* R1 */ \
	 xmm2 -= xmm1;		/* R1 = R1 - R3 (new R3)	;20-23 */ \
	xmm6 *= xmm0;		/* B2 = B2 * sine (new I2)	;21-26 */ \
	 xmm0 = R5;		/* I1 */ \
	 xmm0 -= xmm3;		/* I1 = I1 - I3 (new I3)	;22-25 */ \
	xmm7 *= xptr(screg+off+64);	/* B4 = B4 * sine (new I4)	;23-28 */ \
	 xmm1 += R1;		/* R3 = R1 + R3 (new R1)	;24-27 */ \
	xprefetchw(u8ptr(pre1+pre2)); \
	xmm3 += R5;		/* I3 = I1 + I3 (new I1)	;28-31 */ \
\
	temp2 = xmm4; \
	xmm4 -= xmm5;		/* R2 = R2 - R4 (new R4)	;26-29 */ \
	xmm5 += temp2;		/* R4 = R2 + R4 (new R2)	;34-37 */ \
	temp2 = xmm6; \
	xmm6 -= xmm7;		/* I2 = I2 - I4 (new I4)	;30-33 */ \
	xmm7 += temp2;		/* I4 = I2 + I4 (new I2)	;36-39 */ \
\
	temp2 = xmm0; \
	xmm0 -= xmm4;		/* I3 = I3 - R4 (final I4)	;32-35 */ \
	xmm4 += temp2;		/* R4 = I3 + R4 (final I3)	;44-47 */ \
	temp2 = xmm2; \
	xmm2 -= xmm6;		/* R3 = R3 - I4 (final R3)	;38-41 */ \
	xmm6 += temp2;		/* I4 = R3 + I4 (final R4)	;46-49 */ \
	temp2 = xmm1; \
	xmm1 -= xmm5;		/* R1 = R1 - R2 (final R2)	;40-43 */ \
	xmm5 += temp2;		/* R2 = R1 + R2 (final R1)	;48-51 */ \
	temp2 = xmm3; \
	xmm3 -= xmm7;		/* I1 = I1 - I2 (final I2)	;42-45 */ \
	xmm7 += temp2;		/* I2 = I1 + I2 (final I1)	;50-53 */ \
}

#else
#define x4c_fft_mem_dst12(R1,R2,R3,R4,R5,R6,R7,R8,screg,off,pre1,pre2,dst1,dst2) \
	x4c_fft_mem_PART1(R1,R2,R3,R4,R5,R6,R7,R8,screg,off,pre1,pre2) \
	x4c_fft_mem_PART2_dst1(dst1, exec, dst2)
#define x4c_fft_mem_dst1(R1,R2,R3,R4,R5,R6,R7,R8,screg,off,pre1,pre2,dst1) \
	x4c_fft_mem_PART1(R1,R2,R3,R4,R5,R6,R7,R8,screg,off,pre1,pre2) \
	x4c_fft_mem_PART2_dst1(dst1, exec, dst2)
#define x4c_fft_mem_dst2(R1,R2,R3,R4,R5,R6,R7,R8,screg,off,pre1,pre2,dst2) \
	x4c_fft_mem_PART1(R1,R2,R3,R4,R5,R6,R7,R8,screg,off,pre1,pre2) \
	x4c_fft_mem_PART2_nodst1 \
	x4c_fft_mem_PART3_nodst1_dst2(dst2)
#define x4c_fft_mem_nodst(R1,R2,R3,R4,R5,R6,R7,R8,screg,off,pre1,pre2) \
	x4c_fft_mem_PART1(R1,R2,R3,R4,R5,R6,R7,R8,screg,off,pre1,pre2) \
	x4c_fft_mem_PART2_nodst1 \
	x4c_fft_mem_PART3_nodst1_nodst2

#define x4c_fft_mem_PART1(R1,R2,R3,R4,R5,R6,R7,R8,screg,off,pre1,pre2) \
	xmm0 = R3;		/* R3 */ \
	xmm1 = xptr(screg+off+32+16);	/* cosine/sine */ \
	xmm1 *= xmm0;		/* A3 = R3 * cosine/sine	;1-6 */ \
	xmm2 = R7;		/* I3 */ \
	xmm3 = xptr(screg+off+32+16);	/* cosine/sine */ \
	xmm3 *= xmm2;		/* B3 = I3 * cosine/sine	;3-8 */ \
	xmm4 = R2;		/* R2 */ \
	xmm6 = xptr(screg+off+0+16);	/* cosine/sine */ \
	xmm4 *= xmm6;		/* A2 = R2 * cosine/sine	;5-10 */ \
	xmm5 = R4;		/* R4 */ \
	xmm7 = xptr(screg+off+64+16);	/* cosine/sine */ \
	xmm5 *= xmm7;		/* A4 = R4 * cosine/sine	;7-12 */ \
	xmm1 -= xmm2;		/* A3 = A3 - I3			;8-11 */ \
	xmm2 = R6;		/* I2 */ \
	xmm6 *= xmm2;		/* B2 = I2 * cosine/sine	;9-14 */ \
	xmm3 += xmm0;		/* B3 = B3 + R3			;10-13 */ \
	xmm0 = R8;		/* I4 */ \
	xmm7 *= xmm0;		/* B4 = I4 * cosine/sine	;11-16 */ \
	xmm4 -= xmm2;		/* A2 = A2 - I2			;12-15 */ \
	xmm2 = xptr(screg+off+32);	/* sine */ \
	xmm1 *= xmm2;		/* A3 = A3 * sine (new R3)	;13-18 */ \
	xmm5 -= xmm0;		/* A4 = A4 - I4			;14-17 */ \
	xmm3 *= xmm2;		/* B3 = B3 * sine (new I3)	;15-20 */ \
	xmm6 += R2;		/* B2 = B2 + R2			;16-19 */ \
	xmm0 = xptr(screg+off+0);	/* sine */ \
	xmm4 *= xmm0;		/* A2 = A2 * sine (new R2)	;17-22 */ \
	xprefetchw(u8ptr(pre1)); \
	xmm7 += R4;		/* B4 = B4 + R4			;18-21 */ \
	xmm5 *= xptr(screg+off+64);	/* A4 = A4 * sine (new R4)	;19-24 */ \
	 xmm2 = R1;		/* R1 */ \
	 xmm2 -= xmm1;		/* R1 = R1 - R3 (new R3)	;20-23 */ \
	xmm6 *= xmm0;		/* B2 = B2 * sine (new I2)	;21-26 */ \
	 xmm0 = R5;		/* I1 */ \
	 xmm0 -= xmm3;		/* I1 = I1 - I3 (new I3)	;22-25 */ \
	xmm7 *= xptr(screg+off+64);	/* B4 = B4 * sine (new I4)	;23-28 */ \
	 xmm1 += R1;		/* R3 = R1 + R3 (new R1)	;24-27 */ \
	xprefetchw(u8ptr(pre1+pre2)); \
	 xmm4 -= xmm5;		/* R2 = R2 - R4 (new R4)	;26-29 */ \
	 xmm5 += xmm5; /* multwo */			/* R4 = R4 * 2			;27-32 */ \
	 xmm3 += R5;		/* I3 = I1 + I3 (new I1)	;28-31 */ \
	 xmm6 -= xmm7;		/* I2 = I2 - I4 (new I4)	;30-33 */ \
	 xmm7 += xmm7; /* multwo */			/* I4 = I4 * 2			;31-36 */ \
	xmm0 -= xmm4;		/* I3 = I3 - R4 (final I4)	;32-35 */

#define x4c_fft_mem_PART2_dst1(dst1, dst2w, dst2) \
	dst1 = xmm0; \
	 xmm5 += xmm4;		/* R4 = R2 + R4 (new R2)	;34-37 */ \
	xmm4 += xmm4; /* multwo */			/* R4 = R4 * 2			;35-40 */ \
	 xmm7 += xmm6;		/* I4 = I2 + I4 (new I2)	;36-39 */ \
	xmm4 += xmm0;		/* R4 = I3 + R4 (final I3)	;44-47 */ \
	xmm0 = xmm2; \
	xmm2 -= xmm6;		/* R3 = R3 - I4 (final R3)	;38-41 */ \
dst2w(dst2 = xmm2); \
	xmm6 += xmm0;		/* I4 = R3 + I4 (final R4)	;46-49 */ \
	xmm0 = xmm1; \
	xmm1 -= xmm5;		/* R1 = R1 - R2 (final R2)	;40-43 */ \
	xmm5 += xmm0;		/* R2 = R1 + R2 (final R1)	;48-51 */ \
	xmm0 = xmm3; \
	xmm3 -= xmm7;		/* I1 = I1 - I2 (final I2)	;42-45 */ \
	xmm7 += xmm0;		/* I2 = I1 + I2 (final I1)	;50-53 */

#define x4c_fft_mem_PART2_nodst1 \
	 xmm5 += xmm4;		/* R4 = R2 + R4 (new R2)	;34-37 */ \
	xmm4 += xmm4; /* multwo */			/* R4 = R4 * 2			;35-40 */ \
	 xmm7 += xmm6;		/* I4 = I2 + I4 (new I2)	;36-39 */ \
	xmm2 -= xmm6;		/* R3 = R3 - I4 (final R3)	;38-41 */
#define x4c_fft_mem_PART3_nodst1_dst2(dst2) \
	dst2 = xmm2; \
	xmm6 += xmm6; /* multwo */			/* I4 = I4 * 2			;39-44 */ \
	xmm6 += xmm2;		/* I4 = R3 + I4 (final R4)	;46-49 */ \
	xmm2 = xmm1; \
	xmm1 -= xmm5;		/* R1 = R1 - R2 (final R2)	;40-43 */ \
	xmm5 += xmm2;		/* R2 = R1 + R2 (final R1)	;48-51 */ \
	xmm2 = xmm3; \
	xmm3 -= xmm7;		/* I1 = I1 - I2 (final I2)	;42-45 */ \
	xmm4 += xmm0;		/* R4 = I3 + R4 (final I3)	;44-47 */ \
	xmm7 += xmm2;		/* I2 = I1 + I2 (final I1)	;50-53 */
#define x4c_fft_mem_PART3_nodst1_nodst2 \
	xmm6 += xmm6; /* multwo */			/* I4 = I4 * 2			;39-44 */ \
	xmm1 -= xmm5;		/* R1 = R1 - R2 (final R2)	;40-43 */ \
	xmm5 += xmm5; /* multwo */			/* R2 = R2 * 2			;41-46 */ \
	xmm3 -= xmm7;		/* I1 = I1 - I2 (final I2)	;42-45 */ \
	xmm7 += xmm7; /* multwo */			/* I2 = I2 * 2			;43-48 */ \
	xmm4 += xmm0;		/* R4 = I3 + R4 (final I3)	;44-47 */ \
	xmm6 += xmm2;		/* I4 = R3 + I4 (final R4)	;46-49 */ \
	xmm5 += xmm1;		/* R2 = R1 + R2 (final R1)	;48-51 */ \
	xmm7 += xmm3;		/* I2 = I1 + I2 (final I1)	;50-53 */
#endif


/* Macro to operate on 4 64-byte cache lines.  It does 4 four_complex_unffts */ \
/* in a one pass FFT. */ \
/* 98 clocks */ \

#define x4cl_four_complex_unfft(srcreg,srcinc,d1,d2,screg) \
	xmm0 = xptr(srcreg);		/* R1 */ \
	xmm1 = xptr(srcreg+32);	/* R2 */ \
	xmm2 = xptr(srcreg+d1);	/* R3 */ \
	xmm3 = xptr(srcreg+d1+32);	/* R4 */ \
	xmm4 = xptr(srcreg+d2);	/* R5 */ \
	xmm5 = xptr(srcreg+d2+32);	/* R6 */ \
	best_x4c_unfft(xmm0, xmm1, xmm2, xmm3, xmm4, xmm5, xmm6, xmm7, xptr(srcreg+d2+d1), xptr(srcreg+d2+d1+32), xptr(srcreg), xptr(srcreg+32), screg, 0, srcreg+srcinc, d1); \
	xmm2 = xptr(srcreg+16);	/* R1 */ \
	xmm3 = xptr(srcreg+d2+16);	/* R5 */ \
	xptr(srcreg+d2) = xmm1;	/* Save R2 */ \
	xptr(srcreg+d2+16) = xmm4;	/* Save R4 */ \
	xmm1 = xptr(srcreg+48);	/* R2 */ \
	xmm4 = xptr(srcreg+d2+48);	/* R6 */ \
	xptr(srcreg+d2+32) = xmm0;	/* Save R6 */ \
	xptr(srcreg+d2+48) = xmm7;	/* Save R8 */ \
	xptr(srcreg+16) = xmm5;	/* Save R3 */ \
	xptr(srcreg+48) = xmm6;	/* Save R7 */ \
	xmm0 = xptr(srcreg+d1+16);	/* R3 */ \
	xmm5 = xptr(srcreg+d1+48);	/* R4 */ \
	best_x4c_unfft(xmm2, xmm1, xmm0, xmm5, xmm3, xmm4, xmm6, xmm7, xptr(srcreg+d2+d1+16), xptr(srcreg+d2+d1+48), xptr(srcreg+d1), xptr(srcreg+d1+32), screg, 0, srcreg+srcinc+d2, d1); \
	xptr(srcreg+d2+d1) = xmm1;	/* Save R2 */ \
	xptr(srcreg+d2+d1+16) = xmm3;	/* Save R4 */ \
	xptr(srcreg+d2+d1+32) = xmm2;	/* Save R6 */ \
	xptr(srcreg+d2+d1+48) = xmm7;	/* Save R8 */ \
	xptr(srcreg+d1+16) = xmm4;	/* Save R3 */ \
	xptr(srcreg+d1+48) = xmm6;	/* Save R7 */ \
	srcreg += srcinc; \
\

#define best_x4c_unfft(r1, r2, r3, r4, r5, r6, r7, r8, mem7, mem8, dest1, dest2, screg, off, ... /*, pre1, pre2*/ ) { \
	vec2f64 temp; \
	temp = r1; \
	r1 -= r3;			/* new R2 = R1 - R2 */ \
	r3 += temp;			/* new R1 = R1 + R2 */ \
	temp = mem8; \
	r8 = temp + r6;		/* new I3 = I3 + I4 */ \
	r6 -= temp;		/* new R4 = I3 - I4 */ \
	temp = r2; \
	r2 -= r4;			/* new I2 = I1 - I2 */ \
	r4 += temp;			/* new I1 = I1 + I2 */ \
	temp = mem7; \
	r7 = temp - r5;		/* new I4 = R4 - R3 */ \
	r5 += temp;			/* new R3 = R3 + R4 */ \
	/*IFNB <pre1> \
	xprefetchw(u8ptr(pre1)); \
	ENDIF*/ \
	temp = r1; \
	r1 -= r6;			/* R2 = R2 - R4 (new R4) */ \
	r6 += temp;			/* R4 = R2 + R4 (new R2) */ \
	temp = r2; \
	r2 -= r7;			/* I2 = I2 - I4 (new I4) */ \
	r7 += temp;			/* I4 = I2 + I4 (new I2) */ \
	temp = r3; \
	r3 -= r5;			/* R1 = R1 - R3 (new R3) */ \
	r5 += temp;			/* R3 = R1 + R3 (new & final R1) */ \
	/*IFNB <pre1> \
	xprefetchw(u8ptr(pre1+pre2)); \
	ENDIF*/ \
	dest1 = r5;		/* Save final R1 */ \
	r5 = r4; \
	r4 -= r8;			/* I1 = I1 - I3 (new I3) */ \
	r8 += r5;			/* I3 = I1 + I3 (new & final I1) */ \
	dest2 = r8;		/* Save final I1 */ \
	r5 = xptr(screg+off+64+16);	/* cosine/sine */ \
	r5 *= r1;			/* A4 = new R4 * cosine/sine */ \
	r8 = xptr(screg+off+64+16);	/* cosine/sine */ \
	r8 *= r2;			/* B4 = new I4 * cosine/sine */ \
	r5 += r2;			/* A4 = A4 + new I4 */ \
	r8 -= r1;			/* B4 = B4 - new R4 */ \
	r5 *= xptr(screg+off+64);	/* A4 = A4 * sine (final R4) */ \
	r8 *= xptr(screg+off+64);	/* B4 = B4 * sine (final I4) */ \
	r2 = xptr(screg+off+0+16);	/* cosine/sine */ \
	r2 *= r6;			/* A2 = new R2 * cosine/sine */ \
	r1 = xptr(screg+off+0+16);	/* cosine/sine */ \
	r1 *= r7;			/* B2 = new I2 * cosine/sine */ \
	r2 += r7;			/* A2 = A2 + new I2 */ \
	r1 -= r6;			/* B2 = B2 - new R2 */ \
	r6 = xptr(screg+off+32+16);	/* cosine/sine */ \
	r6 *= r3;			/* A3 = new R3 * cosine/sine */ \
	r7 = xptr(screg+off+32+16);	/* cosine/sine */ \
	r7 *= r4;			/* B3 = new I3 * cosine/sine */ \
	r6 += r4;			/* A3 = A3 + new I3 */ \
	r7 -= r3;			/* B3 = B3 - new R3 */ \
	r2 *= xptr(screg+off+0);	/* A2 = A2 * sine (final R2) */ \
	r1 *= xptr(screg+off+0);	/* B2 = B2 * sine (final I2) */ \
	r6 *= xptr(screg+off+32);	/* A3 = A3 * sine (final R3) */ \
	r7 *= xptr(screg+off+32);	/* B3 = B3 * sine (final I3) */ \
}

#define new_x4c_unfft(r1, r2, r3, r4, r5, r6, r7, r8, dest1, off) { \
	vec2f64 temp; \
	temp = r1; \
	r1 -= r3;			/* new R2 = R1 - R2 */ \
	r3 += temp;			/* new R1 = R1 + R2 */ \
	temp = r6; \
	r6 -= r8;			/* new R4 = I3 - I4 */ \
	r8 += temp;			/* new I3 = I3 + I4 */ \
	temp = r2; \
	r2 -= r4;			/* new I2 = I1 - I2 */ \
	r4 += temp;			/* new I1 = I1 + I2 */ \
	temp = r7; \
	r7 -= r5;			/* new I4 = R4 - R3 */ \
	r5 += temp;			/* new R3 = R3 + R4 */ \
	temp = r1; \
	r1 -= r6;			/* R2 = R2 - R4 (new R4) */ \
	r6 += temp;			/* R4 = R2 + R4 (new R2) */ \
	temp = r2; \
	r2 -= r7;			/* I2 = I2 - I4 (new I4) */ \
	r7 += temp;			/* I4 = I2 + I4 (new I2) */ \
	temp = r3; \
	r3 -= r5;			/* R1 = R1 - R3 (new R3) */ \
	r5 += temp;			/* R3 = R1 + R3 (new & final R1) */ \
	dest1 = r5;		/* Save final R1 */ \
	r2 *= xptr(rdi+off+64);	/* B4 = new I4 * sine */ \
	r1 *= xptr(rdi+off+64);	/* A4 = new R4 * sine */ \
	r5 = xptr(rdi+off+64+16);	/* cosine/sine */ \
	r5 *= r2;			/* C4 = B4 * cosine/sine */ \
	r5 -= r1;			/* C4 = C4 - A4 (final I4) */ \
	r1 *= xptr(rdi+off+64+16);	/* A4 = A4 * cosine/sine */ \
	r1 += r2;			/* A4 = B4 + A4 (final R4) */ \
	temp = r4; \
	r4 -= r8;			/* I1 = I1 - I3 (new I3) */ \
	r8 += temp;			/* I3 = I1 + I3 (new & final I1) */ \
	r2 = xptr(rdi+off+0+16);	/* cosine/sine */ \
	r2 *= r6;			/* A2 = new R2 * cosine/sine */ \
	r2 += r7;			/* A2 = A2 + new I2 */ \
	r7 *= xptr(rdi+off+0+16);	/* B2 = new I2 * cosine/sine */ \
	r7 -= r6;			/* B2 = B2 - new R2 */ \
	r6 = xptr(rdi+off+32+16);	/* cosine/sine */ \
	r6 *= r3;			/* A3 = new R3 * cosine/sine */ \
	r6 += r4;			/* A3 = A3 + new I3 */ \
	r4 *= xptr(rdi+off+32+16);	/* B3 = new I3 * cosine/sine */ \
	r4 -= r3;			/* B3 = B3 - new R3 */ \
	r2 *= xptr(rdi+off+0);		/* A2 = A2 * sine (final R2) */ \
	r7 *= xptr(rdi+off+0);		/* B2 = B2 * sine (final I2) */ \
	r6 *= xptr(rdi+off+32);	/* A3 = A3 * sine (final R3) */ \
	r4 *= xptr(rdi+off+32);	/* B3 = B3 * sine (final I3) */ \
}

#define x4c_unfft(r1, r2, r3, r4, r5, r6, r7, r8) { \
	vec2f64 tmp1, tmp2, tmp3, tmp4; \
	tmp4 = r1; \
	r1 -= r3;			/* new R2 = R1 - R2 */ \
	r3 += tmp4;			/* new R1 = R1 + R2 */ \
	tmp4 = r6; \
	r6 -= r8;			/* new R4 = I3 - I4 */ \
	r8 += tmp4;			/* new I3 = I3 + I4 */ \
	tmp4 = r2; \
	r2 -= r4;			/* new I2 = I1 - I2 */ \
	r4 += tmp4;			/* new I1 = I1 + I2 */ \
	tmp4 = r7; \
	r7 -= r5;			/* new I4 = R4 - R3 */ \
	r5 += tmp4;			/* new R3 = R3 + R4 */ \
	tmp4 = r1; \
	r1 -= r6;			/* R2 = R2 - R4 (new R4) */ \
	r6 += tmp4;			/* R4 = R2 + R4 (new R2) */ \
	tmp1 = r6;		/* Save new R2 */ \
	tmp4 = r2; \
	r2 -= r7;			/* I2 = I2 - I4 (new I4) */ \
	r7 += tmp4;			/* I4 = I2 + I4 (new I2) */ \
	tmp4 = r3; \
	r3 -= r5;			/* R1 = R1 - R3 (new R3) */ \
	r5 += tmp4;		/* R3 = R1 + R3 (new & final R1) */ \
	tmp2 = r3;		/* Save new R3 */ \
	r2 *= xptr(rdi+64);		/* B4 = new I4 * sine */ \
	r1 *= xptr(rdi+64);		/* A4 = new R4 * sine */ \
	r6 *= xptr(rdi+0+16);		/* A2 = new R2 * cosine/sine */ \
	tmp3 = r2;		/* Save B4 */ \
	r3 *= xptr(rdi+32+16);		/* A3 = new R3 * cosine/sine */ \
	r2 *= xptr(rdi+64+16);		/* C4 = B4 * cosine/sine */ \
	tmp4 = r4; \
	r4 -= r8;			/* I1 = I1 - I3 (new I3) */ \
	r8 += tmp4;		/* I3 = I1 + I3 (new & final I1) */ \
	r6 += r7;			/* A2 = A2 + new I2 */ \
	r7 *= xptr(rdi+0+16);		/* B2 = new I2 * cosine/sine */ \
	r3 += r4;			/* A3 = A3 + new I3 */ \
	r4 *= xptr(rdi+32+16);		/* B3 = new I3 * cosine/sine */ \
	r2 -= r1;			/* C4 = C4 - A4 (final I4) */ \
	r1 *= xptr(rdi+64+16);		/* A4 = A4 * cosine/sine */ \
	r7 -= tmp1;		/* B2 = B2 - new R2 */ \
	r6 *= xptr(rdi+0);		/* A2 = A2 * sine (final R2) */ \
	r4 -= tmp2;		/* B3 = B3 - new R3 */ \
	r3 *= xptr(rdi+32);		/* A3 = A3 * sine (final R3) */ \
	r7 *= xptr(rdi+0);		/* B2 = B2 * sine (final I2) */ \
	r1 += tmp3;		/* A4 = B4 + A4 (final R4) */ \
	r4 *= xptr(rdi+32);		/* B3 = B3 * sine (final I3) */ \
}

/* Macro to operate on 2 64-byte cache lines.  It does 2 four_complex_unffts */
/* in a two pass FFT. */

#define x2cl_four_complex_unfft(srcreg,srcinc,d1) \
	xmm0 = xptr(srcreg);		/* R1 */ \
	xmm1 = xptr(srcreg+16);	/* R2 */ \
	xmm2 = xptr(srcreg+32);	/* R3 */ \
	xmm3 = xptr(srcreg+48);	/* R4 */ \
	xmm4 = xptr(srcreg+d1);	/* R5 */ \
	xmm5 = xptr(srcreg+d1+16);	/* R6 */ \
	best_x4c_unfft(xmm0, xmm1, xmm2, xmm3, xmm4, xmm5, xmm6, xmm7, xptr(srcreg+d1+32), xptr(srcreg+d1+48), xptr(srcreg), xptr(srcreg+32), rdi, 0, srcreg+srcinc, d1); \
	xptr(srcreg+d1) = xmm1;	/* Save R2 */ \
	xptr(srcreg+d1+16) = xmm4;	/* Save R4 */ \
	xptr(srcreg+d1+32) = xmm0;	/* Save R6 */ \
	xptr(srcreg+d1+48) = xmm7;	/* Save R8 */ \
	xptr(srcreg+16) = xmm5;	/* Save R3 */ \
	xptr(srcreg+48) = xmm6;	/* Save R7 */ \
	srcreg += srcinc; \
\

#define g2cl_four_complex_unfft(srcreg,srcinc,d1,dstreg,dstinc,e1) \
	xmm0 = xptr(srcreg);		/* R1 */ \
	xmm1 = xptr(srcreg+16);	/* R2 */ \
	xmm2 = xptr(srcreg+32);	/* R3 */ \
	xmm3 = xptr(srcreg+48);	/* R4 */ \
	xmm4 = xptr(srcreg+d1);	/* R5 */ \
	xmm5 = xptr(srcreg+d1+16);	/* R6 */ \
	best_x4c_unfft(xmm0, xmm1, xmm2, xmm3, xmm4, xmm5, xmm6, xmm7, xptr(srcreg+d1+32), xptr(srcreg+d1+48), xptr(dstreg), xptr(dstreg+32), rdi, 0, dstreg+dstinc, e1); \
	srcreg += srcinc; \
	xptr(dstreg+e1) = xmm1;	/* Save R2 */ \
	xptr(dstreg+e1+16) = xmm4;	/* Save R4 */ \
	xptr(dstreg+e1+32) = xmm0;	/* Save R6 */ \
	xptr(dstreg+e1+48) = xmm7;	/* Save R8 */ \
	xptr(dstreg+16) = xmm5;	/* Save R3 */ \
	xptr(dstreg+48) = xmm6;	/* Save R7 */ \
	dstreg += dstinc; \
\
\
\
/* Macro to operate on 2 64-byte cache line.  It does does one FFT level */ \
/* given 8 reals. */ \

#define x2cl_four_reals_fft_1(srcreg,srcinc,d1) \
	xmm0 = xptr(srcreg);		/* R1 */ \
	xmm1 = xptr(srcreg+16);	/* R2 */ \
	xmm2 = xptr(srcreg+32);	/* R3 */ \
	xmm3 = xptr(srcreg+48);	/* R4 */ \
	half_x4r_fft(xmm0, xmm1, xmm2, xmm3, xmm4); \
	xptr(srcreg) = xmm1;		/* Save R1 */ \
	xptr(srcreg+16) = xmm0;	/* Save R2 */ \
	xptr(srcreg+32) = xmm2;	/* Save R3 */ \
	xptr(srcreg+48) = xmm3;	/* Save R4 */ \
	xmm4 = xptr(srcreg+d1);	/* R1 */ \
	xmm5 = xptr(srcreg+d1+16);	/* R2 */ \
	xmm6 = xptr(srcreg+d1+32);	/* R3 */ \
	xmm7 = xptr(srcreg+d1+48);	/* R4 */ \
	half_x4r_fft(xmm4, xmm5, xmm6, xmm7, xmm0); \
	xptr(srcreg+d1) = xmm5;	/* Save R1 */ \
	xptr(srcreg+d1+16) = xmm4;	/* Save R2 */ \
	xptr(srcreg+d1+32) = xmm6;	/* Save R3 */ \
	xptr(srcreg+d1+48) = xmm7;	/* Save R4 */ \
	srcreg += srcinc; \
\

#define half_x4r_fft(r1, r2, r3, r4, t1) \
	t1 = r1; \
	r1 -= r2;			/* R1 - R2 (final R2) */ \
	r2 += t1;			/* R1 + R2 (final R1) */ \
;					/* Nop R3 */ \
;					/* Nop R4 */ \
\
\
/* Macro to operate on 2 64-byte cache lines.  It does does one FFT level */ \
/* given 4 reals and 2 complex values. */ \

#define x2cl_eight_reals_fft_1(srcreg,srcinc,d1) \
	xmm0 = xptr(srcreg);		/* R1 */ \
	xmm1 = xptr(srcreg+d1);	/* R2 */ \
	xmm4 = xptr(srcreg+32);	/* R5 */ \
	xmm5 = xptr(srcreg+d1+32);	/* R6 */ \
	xmm6 = xptr(srcreg+48);	/* R7 */ \
	xmm7 = xptr(srcreg+d1+48);	/* R8 */ \
	x4r_fft(xmm0, xmm1, xmm2, xmm3, xmm4, xmm5, xmm6, xmm7); \
	xmm2 = xptr(srcreg+16);	/* R3 */ \
	xmm3 = xptr(srcreg+d1+16);	/* R4 */ \
	xptr(srcreg) = xmm1;		/* Save R1 */ \
	xptr(srcreg+16) = xmm0;	/* Save R2 */ \
	xptr(srcreg+32) = xmm2;	/* Save R3 */ \
	xptr(srcreg+48) = xmm3;	/* Save R4 */ \
	xptr(srcreg+d1) = xmm5;	/* Save R5 */ \
	xptr(srcreg+d1+16) = xmm7;	/* Save R6 */ \
	xptr(srcreg+d1+32) = xmm4;	/* Save R7 */ \
	xptr(srcreg+d1+48) = xmm6;	/* Save R8 */ \
	srcreg += srcinc; \
\

#define s2cl_eight_reals_fft_1(srcreg,srcinc,d1) { \
	uintptr_t src_rbx = srcreg+rbx; \
	shuffle_load(xmm0,xmm1,xptr(src_rbx),xptr(src_rbx+d1)); /* R1,R2 */ \
	shuffle_load(xmm4,xmm5,xptr(src_rbx+32),xptr(src_rbx+d1+32)); /* R5,R6 */ \
	shuffle_load(xmm6,xmm7,xptr(src_rbx+48),xptr(src_rbx+d1+48)); /* R7,R8 */ \
	x4r_fft(xmm0, xmm1, xmm2, xmm3, xmm4, xmm5, xmm6, xmm7); \
	shuffle_load(xmm2,xmm3,xptr(src_rbx+16),xptr(src_rbx+d1+16)); /* R3,R4 */ \
	xptr(srcreg) = xmm1;		/* Save R1 */ \
	xptr(srcreg+16) = xmm0;	/* Save R2 */ \
	xptr(srcreg+32) = xmm2;	/* Save R3 */ \
	xptr(srcreg+48) = xmm3;	/* Save R4 */ \
	xptr(srcreg+d1) = xmm5;	/* Save R5 */ \
	xptr(srcreg+d1+16) = xmm7;	/* Save R6 */ \
	xptr(srcreg+d1+32) = xmm4;	/* Save R7 */ \
	xptr(srcreg+d1+48) = xmm6;	/* Save R8 */ \
	srcreg += srcinc; \
}

#define x4r_fft(r1, r2, t3, t4, r5, r6, r7, r8) \
	t3 = r1; \
	r1 -= r2;			/* R1 - R2 (final R2) */ \
	r2 += t3;			/* R1 + R2 (final R1) */ \
;					/* Nop R3 */ \
;					/* Nop R4 */ \
	t4 = r6; \
	r6 -= r8;			/* R2 - I2 */ \
	r8 += t4;			/* R2 + I2 */ \
	t3 = XMM_SQRTHALF2; \
	r6 *= t3;			/* newR2 */ \
	r8 *= t3;			/* newI2 */ \
	t4 = r7; \
	r7 -= r8;			/* I1 = I1 - I2 (new I2) */ \
	r8 += t4;			/* I2 = I1 + I2 (new I1) */ \
	t3 = r5; \
	r5 -= r6;			/* R1 = R1 - R2 (new R2) */ \
	r6 += t3;			/* R2 = R1 + R2 (new R1) */ \
\
\
\
/* Macro to operate on 2 64-byte cache lines.  It does one level */ \
/* of inverse FFT. */ \

#define x2cl_eight_reals_unfft_1(srcreg,srcinc,d1) \
	xmm0 = xptr(srcreg);		/* R1 */ \
	xmm1 = xptr(srcreg+32);	/* R2 */ \
	xmm4 = xptr(srcreg+16);	/* R1 */ \
	xmm5 = xptr(srcreg+48);	/* R2 */ \
	xmm2 = xptr(srcreg+d1);	/* R3 */ \
	xmm3 = xptr(srcreg+d1+32);	/* R4 */ \
	xmm6 = xptr(srcreg+d1+16);	/* R3 */ \
/*	xload	xmm7, [srcreg+d1+48]	;; R4 */ \
	xmm7 = xmm0; \
	xmm0 -= xmm1;		/* R1 - R2 (final R2) */ \
	xmm1 += xmm7;		/* R1 + R2 (final R1) */ \
	xmm7 = xmm4; \
	xmm4 -= xmm5;		/* R1 - R2 (final R2) */ \
	xmm5 += xmm7;		/* R1 + R2 (final R1) */ \
	xptr(srcreg) = xmm1;		/* R1 */ \
	xptr(srcreg+16) = xmm0;	/* R2 */ \
	xptr(srcreg+32) = xmm2;	/* R3 */ \
	xptr(srcreg+48) = xmm3;	/* R4 */ \
	xptr(srcreg+d1) = xmm5;	/* R5 */ \
	xptr(srcreg+d1+16) = xmm4;	/* R6 */ \
	xptr(srcreg+d1+32) = xmm6;	/* R7 */ \
/*	xstore	[srcreg+d1+48], xmm7	;; R8 */ \
	srcreg += srcinc; \
\

#define x1cl_eight_reals_unfft_1(srcreg,srcinc) \
	xmm0 = xptr(srcreg);		/* R1 */ \
	xmm1 = xptr(srcreg+16);	/* R2 */ \
	xmm2 = xmm0; \
	xmm0 -= xmm1;		/* R1 - R2 (final R2) */ \
	xmm1 += xmm2;		/* R1 + R2 (final R1) */ \
	xptr(srcreg) = xmm1;		/* R1 */ \
	xptr(srcreg+16) = xmm0;	/* R2 */ \
	srcreg += srcinc; \
\
\
/* Macro to operate on 4 64-byte cache lines.  It does one level */ \
/* of inverse FFT. */ \

#define s4cl_eight_reals_unfft_1(srcreg,srcinc,d1,d2) \
	xmm0 = xptr(srcreg);		/* R1 */ \
	xmm1 = xptr(srcreg+32);	/* R2 */ \
	xmm4 = xptr(srcreg+d2);	/* R5 */ \
	xmm5 = xptr(srcreg+d2+32);	/* R6 */ \
	xmm6 = xptr(srcreg+d2+d1);	/* R7 */ \
	xmm7 = xptr(srcreg+d2+d1+32);	/* R8 */ \
	x8r1_unfft(xmm0, xmm1, xmm2, xmm3, xmm4, xmm5, xmm6, xmm7, xptr(srcreg+d1), xptr(srcreg+d1+32)); \
	shuffle_store(xptr(srcreg), xptr(srcreg+d2), xmm1, xmm0);	/* Save R1 and R2 */ \
	xmm0 = xptr(srcreg+16);	/* R1 */ \
	xmm1 = xptr(srcreg+d2+16);	/* R5 */ \
	shuffle_store(xptr(srcreg+16), xptr(srcreg+d2+16), xmm2, xmm3); /* Save R3 and R4 */ \
	xmm2 = xptr(srcreg+48);	/* R2 */ \
	xmm3 = xptr(srcreg+d2+48);	/* R6 */ \
	shuffle_store(xptr(srcreg+32), xptr(srcreg+d2+32), xmm6, xmm4); /* Save R5 and R6 */ \
	shuffle_store(xptr(srcreg+48), xptr(srcreg+d2+48), xmm7, xmm5); /* Save R7 and R8 */ \
	xmm6 = xptr(srcreg+d2+d1+16);	/* R7 */ \
	xmm7 = xptr(srcreg+d2+d1+48);	/* R8 */ \
	x8r1_unfft(xmm0, xmm2, xmm4, xmm5, xmm1, xmm3, xmm6, xmm7, xptr(srcreg+d1+16), xptr(srcreg+d1+48)); \
	shuffle_store(xptr(srcreg+d1), xptr(srcreg+d2+d1), xmm2, xmm0); /* Save R1 and R2 */ \
	shuffle_store(xptr(srcreg+d1+16), xptr(srcreg+d2+d1+16), xmm4, xmm5); /* R3, R4 */ \
	shuffle_store(xptr(srcreg+d1+32), xptr(srcreg+d2+d1+32), xmm6, xmm1); /* R5, R6 */ \
	shuffle_store(xptr(srcreg+d1+48), xptr(srcreg+d2+d1+48), xmm7, xmm3); /* R7, R8 */ \
	srcreg += srcinc; \
\

#define x8r1_unfft(r1, r2, r3, r4, r5, r6, r7, r8, mem3, mem4) \
	r3 = r1; \
	r1 -= r2;			/* R1 - R2 (final R2) */ \
	r2 += r3;			/* R1 + R2 (final R1) */ \
	r4 = r6; \
	r6 -= r8;			/* new I2 = I1 - I2 */ \
	r8 += r4;			/* new I1 = I1 + I2 */ \
	r3 = r5; \
	r5 -= r7;			/* new R2 = R1 - R2 */ \
	r7 += r3;			/* new R1 = R1 + R2 */ \
	r4 = XMM_SQRTHALF2; \
	r6 *= r4;			/* B2 = I2 * sine */ \
	r5 *= r4;			/* A2 = R2 * sine */ \
	r3 = r6;			/* Save B2 (C2 = B2) */ \
	r6 -= r5;			/* C2 = C2 - A2 (new I2) */ \
	r5 += r3;			/* A2 = A2 + B2 (new R2) */ \
	r3 = mem3; \
	r4 = mem4; \
\
\
/* Macro to operate on 2 64-byte cache lines.  It does 2 two_complex_ffts */ \
/* in a one pass FFT. */ \
/* 30.45 clocks */ \

#define x2cl_two_complex_fft(srcreg,srcinc,d1,screg) \
	xmm2 = xptr(srcreg+d1);	/* R2 */ \
	xmm0 = xptr(screg+48);		/* cosine/sine */ \
	xmm0 *= xmm2;		/* A2 = R2 * cosine/sine	;1-6 */ \
	xmm3 = xptr(srcreg+d1+16);	/* R4 */ \
	xmm1 = xptr(screg+48);	/* cosine/sine */ \
	xmm1 *= xmm3;		/* B2 = I2 * cosine/sine	;3-8 */ \
	xmm6 = xptr(srcreg+d1+32);	/* R2_2 */ \
	xmm4 = xptr(screg+XMM_SCD+48);/* cosine/sine */ \
	xmm4 *= xmm6;		/* A2_2 = R2_2 * cosine/sine	;5-10 */ \
	xmm7 = xptr(srcreg+d1+48);	/* R4_2 */ \
	xmm5 = xptr(screg+XMM_SCD+48);/* cosine/sine */ \
	xmm5 *= xmm7;		/* B2_2 = I2_2 * cosine/sine	;7-12 */ \
	xmm0 -= xmm3;		/* A2 = A2 - I2			;8-11 */ \
	xmm1 += xmm2;		/* B2 = B2 + R2			;10-13 */ \
	xprefetchw(u8ptr(srcreg+srcinc)); \
	xmm4 -= xmm7;		/* A2_2 = A2_2 - I2_2		;12-15 */ \
	xmm0 *= xptr(screg+32);	/* A2 = A2 * sine (new R2)	;13-18 */ \
	xmm5 += xmm6;		/* B2_2 = B2_2 + R2_2		;14-17 */ \
	xmm1 *= xptr(screg+32);	/* B2 = B2 * sine (new I2)	;15-20 */ \
	xmm4 *= xptr(screg+XMM_SCD+32);/* A2_2 = A2_2 * sine (new R2_2);17-22 */ \
	xmm2 = xptr(srcreg);		/* R1 */ \
	xprefetchw(u8ptr(srcreg+srcinc+d1)); \
	xmm2 -= xmm0;		/* R1 = R1 - R2 (final R2)	;18-21 */ \
	xmm5 *= xptr(screg+XMM_SCD+32);/* B2_2 = B2_2 * sine (new I2_2);19-24 */ \
	xmm0 += xptr(srcreg);		/* R2 = R1 + R2 (final R1)	;20-23 */ \
	xmm3 = xptr(srcreg+16);	/* R3 */ \
	xmm3 -= xmm1;		/* I1 = I1 - I2 (final I2)	;22-25 */ \
	xmm1 += xptr(srcreg+16);	/* I2 = I1 + I2 (final I1)	;24-27 */ \
	xmm6 = xptr(srcreg+32);	/* R1_2 */ \
	xmm6 -= xmm4;		/* R1_2 = R1_2-R2_2 (final R2_2);26-29 */ \
	xmm4 += xptr(srcreg+32);	/* R2_2 = R1_2+R2_2 (final R1_2);28-31 */ \
	xmm7 = xptr(srcreg+48);	/* R3_2 */ \
	xmm7 -= xmm5;		/* I1_2 = I1_2-I2_2 (final I2_2);30-33 */ \
	xmm5 += xptr(srcreg+48);	/* I2_2 = I1_2+I2_2 (final I1_2);32-35 */ \
	xptr(srcreg) = xmm0;		/* Save R1 */ \
	xptr(srcreg+16) = xmm1;	/* Save R2 */ \
	xptr(srcreg+32) = xmm2;	/* Save R3 */ \
	xptr(srcreg+48) = xmm3;	/* Save R4 */ \
	xptr(srcreg+d1) = xmm4;	/* Save R1_2 */ \
	xptr(srcreg+d1+16) = xmm5;	/* Save R2_2 */ \
	xptr(srcreg+d1+32) = xmm6;	/* Save R3_2 */ \
	xptr(srcreg+d1+48) = xmm7;	/* Save R4_2 */ \
	srcreg += srcinc; \
\
\
/* Like the above except that it uses a different memory layout and uses */ \
/* only one sin/cos set.  See xfft48p for an example use. */ \

#define x2cl_two_complex_fft_in_place(srcreg,srcinc,d1,screg) \
	xmm2 = xptr(srcreg+16);	/* R2 */ \
	xmm0 = xptr(screg+48);		/* cosine/sine */ \
	xmm0 *= xmm2;		/* A2 = R2 * cosine/sine	;1-6 */ \
	xmm3 = xptr(srcreg+48);	/* R4 */ \
	xmm1 = xptr(screg+48);	/* cosine/sine */ \
	xmm1 *= xmm3;		/* B2 = I2 * cosine/sine	;3-8 */ \
	xmm6 = xptr(srcreg+d1+16);	/* R2_2 */ \
	xmm4 = xptr(screg+48);	/* cosine/sine */ \
	xmm4 *= xmm6;		/* A2_2 = R2_2 * cosine/sine	;5-10 */ \
	xmm7 = xptr(srcreg+d1+48);	/* R4_2 */ \
	xmm5 = xptr(screg+48);	/* cosine/sine */ \
	xmm5 *= xmm7;		/* B2_2 = I2_2 * cosine/sine	;7-12 */ \
	xmm0 -= xmm3;		/* A2 = A2 - I2			;8-11 */ \
	xmm1 += xmm2;		/* B2 = B2 + R2			;10-13 */ \
	xmm4 -= xmm7;		/* A2_2 = A2_2 - I2_2		;12-15 */ \
	xmm0 *= xptr(screg+32);	/* A2 = A2 * sine (new R2)	;13-18 */ \
	xmm5 += xmm6;		/* B2_2 = B2_2 + R2_2		;14-17 */ \
	xmm1 *= xptr(screg+32);	/* B2 = B2 * sine (new I2)	;15-20 */ \
	xmm4 *= xptr(screg+32);	/* A2_2 = A2_2 * sine (new R2_2);17-22 */ \
	xmm2 = xptr(srcreg);		/* R1 */ \
	xmm2 -= xmm0;		/* R1 = R1 - R2 (final R2)	;18-21 */ \
	xmm5 *= xptr(screg+32);	/* B2_2 = B2_2 * sine (new I2_2);19-24 */ \
	xmm0 += xptr(srcreg);		/* R2 = R1 + R2 (final R1)	;20-23 */ \
	xmm3 = xptr(srcreg+32);	/* R3 */ \
	xmm3 -= xmm1;		/* I1 = I1 - I2 (final I2)	;22-25 */ \
	xmm1 += xptr(srcreg+32);	/* I2 = I1 + I2 (final I1)	;24-27 */ \
	xmm6 = xptr(srcreg+d1);	/* R1_2 */ \
	xmm6 -= xmm4;		/* R1_2 = R1_2-R2_2 (final R2_2);26-29 */ \
	xmm4 += xptr(srcreg+d1);	/* R2_2 = R1_2+R2_2 (final R1_2);28-31 */ \
	xmm7 = xptr(srcreg+d1+32);	/* R3_2 */ \
	xmm7 -= xmm5;		/* I1_2 = I1_2-I2_2 (final I2_2);30-33 */ \
	xmm5 += xptr(srcreg+d1+32);	/* I2_2 = I1_2+I2_2 (final I1_2);32-35 */ \
	xptr(srcreg) = xmm0;		/* Save R1 */ \
	xptr(srcreg+16) = xmm1;	/* Save R2 */ \
	xptr(srcreg+32) = xmm2;	/* Save R3 */ \
	xptr(srcreg+48) = xmm3;	/* Save R4 */ \
	xptr(srcreg+d1) = xmm4;	/* Save R1_2 */ \
	xptr(srcreg+d1+16) = xmm5;	/* Save R2_2 */ \
	xptr(srcreg+d1+32) = xmm6;	/* Save R3_2 */ \
	xptr(srcreg+d1+48) = xmm7;	/* Save R4_2 */ \
	srcreg += srcinc; \
\
\
/* Macro to operate on 2 64-byte cache lines.  It does 2 two_complex_unffts */ \
/* in a one pass FFT. */ \
/* 31.45 clocks */ \

#define x2cl_two_complex_unfft(srcreg,srcinc,d1) \
	xmm0 = xptr(srcreg);		/* R1 */ \
	xmm1 = xptr(srcreg+32);	/* R2 */ \
	xmm2 = xptr(srcreg+d1);	/* R3 */ \
	xmm3 = xptr(srcreg+d1+32);	/* R4 */ \
	x2c_unfft(xmm0, xmm1, xmm2, xmm3, xmm4, xmm5); \
	xptr(srcreg) = xmm2;		/* Save R1 */ \
	xmm2 = xptr(srcreg+16);	/* R1 */ \
	xptr(srcreg+16) = xmm4;	/* Save R2 */ \
	xmm4 = xptr(srcreg+48);	/* R2 */ \
	xptr(srcreg+32) = xmm3;	/* Save R3 */ \
	xmm3 = xptr(srcreg+d1+16);	/* R3 */ \
	xptr(srcreg+48) = xmm5;	/* Save R4 */ \
	xmm5 = xptr(srcreg+d1+48);	/* R4 */ \
	x2c_unfft(xmm2, xmm4, xmm3, xmm5, xmm0, xmm1); \
	xptr(srcreg+d1) = xmm3;	/* Save R1 */ \
	xptr(srcreg+d1+16) = xmm0;	/* Save R2 */ \
	xptr(srcreg+d1+32) = xmm5;	/* Save R3 */ \
	xptr(srcreg+d1+48) = xmm1;	/* Save R4 */ \
	srcreg += srcinc; \

// BUG: asm have unused d1: x1cl_two_complex_unfft(srcreg,srcinc,d1)
#define x1cl_two_complex_unfft(srcreg,srcinc) \
	xmm0 = xptr(srcreg);		/* R1 */ \
	xmm1 = xptr(srcreg+16);	/* R2 */ \
	xmm2 = xptr(srcreg+32);	/* R3 */ \
	xmm3 = xptr(srcreg+48);	/* R4 */ \
	x2c_unfft(xmm0, xmm1, xmm2, xmm3, xmm4, xmm5); \
	xptr(srcreg) = xmm2;		/* Save R1 */ \
	xptr(srcreg+16) = xmm4;	/* Save R2 */ \
	xptr(srcreg+32) = xmm3;	/* Save R3 */ \
	xptr(srcreg+48) = xmm5;	/* Save R4 */ \
	srcreg += srcinc; \
\

#define x2c_unfft(r1, r2, r3, r4, t1, t2) \
	t1 = r1; \
	r1 -= r3;			/* new R2 = R1 - R2 */ \
	r3 += t1;			/* new R1 = R1 + R2 */ \
	t2 = r2; \
	r2 -= r4;			/* new I2 = I1 - I2 */ \
	r4 += t2;			/* new I1 = I1 + I2 */ \
	t1 = xptr(rdi+32+16);		/* cosine/sine */ \
	t1 *= r1;			/* A2 = new R2 * cosine/sine */ \
	t1 += r2;			/* A2 = A2 + new I2 */ \
	t2 = xptr(rdi+32+16);		/* cosine/sine */ \
	t2 *= r2;			/* B2 = new I2 * cosine/sine */ \
	t2 -= r1;			/* B2 = B2 - new R2 */ \
	r1 = xptr(rdi+32);		/* sine */ \
	t1 *= r1;			/* A2 = A2 * sine (final R2) */ \
	t2 *= r1;			/* B2 = B2 * sine (final I2) */ \

#ifdef SHOULD_BE_FASTER_BUT_ISNT
/* 28.35 clocks, but is slower in practice */
#define x2cl_two_complex_unfft(srcreg,srcinc,d1) \
	xmm4 = xptr(srcreg);		/* R1 */ \
	xmm2 = xptr(srcreg+d1);	/* R2 */ \
	xmm4 -= xmm2;		/* new R2 = R1 - R2		;1-4 */ \
	xmm5 = xptr(srcreg+32);	/* I1 */ \
	xmm5 -= xptr(srcreg+d1+32);	/* new I2 = I1 - I2		;3-6 */ \
	xmm6 = xptr(srcreg+16);	/* R1_2 */ \
	xmm6 -= xptr(srcreg+d1+16);	/* new R2_2 = R1_2 - R2_2	;5-8 */ \
	xmm7 = xptr(rdi+48);		/* cosine/sine */ \
	xmm0 = xmm4;		/* save new R2 */ \
	xmm4 *= xmm7;		/* A2 = new R2 * cosine/sine	;6-11 */ \
	xmm3 = xptr(srcreg+48);	/* I1_2 */ \
	xmm3 -= xptr(srcreg+d1+48);	/* new I2_2 = I1_2 - I2_2	;7-10 */ \
	xmm1 = xmm5;		/* save new I2 */ \
	xmm5 *= xmm7;		/* B2 = new I2 * cosine/sine	;8-13 */ \
	xmm4 += xmm1;		/* A2 = A2 + new I2		;13-16 */ \
	xmm1 = xptr(srcreg);		/* R1 */ \
	xmm1 += xmm2;		/* new R1 = R1 + R2		;9-12 */ \
	xmm2 = xmm6;		/* save new R2_2 */ \
	xmm6 *= xmm7;		/* A2_2 = new R2_2 * cosine/sine;10-15 */ \
	xmm5 -= xmm0;		/* B2 = B2 - new R2		;15-18 */ \
	xmm0 = xptr(srcreg+32);	/* I1 */ \
	xmm0 += xptr(srcreg+d1+32);	/* new I1 = I1 + I2		;11-14 */ \
	xmm7 *= xmm3;		/* B2_2 = new I2_2 * cosine/sine;12-17 */ \
	xmm6 += xmm3;		/* A2_2 = A2_2 + new I2_2	;17-20 */ \
	xmm7 -= xmm2;		/* B2_2 = B2_2 - new R2_2	;19-22 */ \
	xmm3 = xptr(srcreg+16);	/* R1_2 */ \
	xmm3 += xptr(srcreg+d1+16);	/* new R1_2 = R1_2 + R2_2	;21-24 */ \
	xptr(srcreg) = xmm1;		/* Save R1 */ \
	xmm1 = xptr(rdi+32);		/* sine */ \
	xmm4 *= xmm1;		/* A2 = A2 * sine (final R2)	;22-27 */ \
	xmm2 = xptr(srcreg+48);	/* I1_2 */ \
	xmm2 += xptr(srcreg+d1+48);	/* new I1_2 = I1_2 + I2_2	;23-26 */ \
	xmm5 *= xmm1;		/* B2 = B2 * sine (final I2)	;24-29 */ \
	xmm6 *= xmm1;		/* A2_2 = A2_2*sine (final R2_2);26-31 */ \
	xmm7 *= xmm1;		/* B2_2 = B2_2*sine (final I2_2);28-33 */ \
	xptr(srcreg+16) = xmm4;	/* Save R2 */ \
	xptr(srcreg+32) = xmm0;	/* Save R3 */ \
	xptr(srcreg+48) = xmm5;	/* Save R4 */ \
	xptr(srcreg+d1) = xmm3;	/* Save R1_2 */ \
	xptr(srcreg+d1+16) = xmm6;	/* Save R2_2 */ \
	xptr(srcreg+d1+32) = xmm2;	/* Save R3_2 */ \
	xptr(srcreg+d1+48) = xmm7;	/* Save R4_2 */ \
	srcreg += srcinc;
#endif
\
\
/************************************************************************ */ \
/* */ \
/* These macros do the last two FFT levels. */ \
/* */ \
/************************************************************************ */ \
\
/* Macros for one pass FFTs and the real data in a two pass FFT. */ \
/* These macros process two cache lines containing the first 16 FFT values. */ \
/* The first 8 values are real data, then next 8 values represent 4 complex */ \
/* numbers.  Note that the high word of each XMM register contains the */ \
/* complex data.  Also note that these macros are called only once and */ \
/* we can be a bit sloppy with the optimization. */ \
\
/* Do an eight_reals_fft_2 on 8 doubles and a four_complex_fft_2 on 8 doubles */ \

#define s2cl_eight_reals_fft_2_final(srcreg,srcinc,d1) \
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
	xs4c_fft(xmm0[0], xmm1[0], xmm4[0], xmm5[0], xmm2[0], xmm3[0], xmm6[0], xmm7[0], f64ptr(srcreg+d1+56)); \
	f64ptr(srcreg+d1) = xmm5[0];	/* Save R1 */ \
	f64ptr(srcreg+d1+8) = xmm1[0];	/* Save R2 */ \
	f64ptr(srcreg+d1+16) = xmm7[0];	/* Save R3 */ \
	f64ptr(srcreg+d1+24) = xmm6[0];	/* Save R4 */ \
	f64ptr(srcreg+d1+32) = xmm0[0];	/* Save R5 */ \
	f64ptr(srcreg+d1+40) = xmm4[0];	/* Save R6 */ \
	f64ptr(srcreg+d1+48) = xmm3[0];	/* Save R7 */ \
	f64ptr(srcreg+d1+56) = xmm2[0];	/* Save R8 */ \
	srcreg += srcinc; \
\
\
/* Do an eight_reals_with_square_2 on 8 doubles and */ \
/* a four_complex_with_square_2 on 8 doubles */ \

#define s2cl_eight_reals_with_square_2(srcreg,srcinc,d1) { \
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
	xs4c_fft(xmm0[0], xmm1[0], xmm4[0], xmm5[0], xmm2[0], xmm3[0], xmm6[0], xmm7[0], f64ptr(srcreg+d1+56)); \
\
	double tmp = xmm0[0]; \
	xs_complex_square(xmm5[0], xmm1[0], xmm0[0]);	/* Square R1, R2 */ \
	xs_complex_square(xmm7[0], xmm6[0], xmm0[0]);	/* Square R3, R4 */ \
	xmm0[0] = tmp; \
	tmp = xmm5[0]; \
	xs_complex_square(xmm0[0], xmm4[0], xmm5[0]);	/* Square R5, R6 */ \
	xs_complex_square(xmm3[0], xmm2[0], xmm5[0]);	/* Square R7, R8 */ \
	xmm5[0] = tmp; \
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
}

/* Do an eight_reals_with_mult_2 on 8 doubles and */
/* a four_complex_with_mult on 8 doubles */

#define s2cl_eight_reals_with_mult_2(srcreg,srcinc,d1) \
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
	xs8r_mulf(xmm3[0], xmm2[0], xmm0[0], xmm1[0], xmm5[0], xmm7[0], xmm4[0], xmm6[0], (srcreg), (srcreg+8), (srcreg+16), (srcreg+24), (srcreg+32), (srcreg+40), (srcreg+48), (srcreg+56)); \
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
	xs4c_fft(xmm0[0], xmm1[0], xmm4[0], xmm5[0], xmm2[0], xmm3[0], xmm6[0], xmm7[0], f64ptr(srcreg+d1+56)); \
\
	xs4c_mulf(xmm5[0], xmm1[0], xmm7[0], xmm6[0], xmm0[0], xmm4[0], xmm3[0], xmm2[0], (srcreg+d1), (srcreg+d1+8), (srcreg+d1+16), (srcreg+d1+24), (srcreg+d1+32), (srcreg+d1+40), (srcreg+d1+48), (srcreg+d1+56)); \
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
\
/* Do an eight_reals_with_mulf_2 on 8 doubles and */ \
/* a four_complex_with_mulf_2 on 8 doubles */ \

#define s2cl_eight_reals_with_mulf_2(srcreg,srcinc,d1) \
	uintptr_t src_rbx = srcreg+rbx; \
\
	xmm3[0] = f64ptr(src_rbx);	/* R1 */ \
	xmm2[0] = f64ptr(src_rbx+8);	/* R2 */ \
	xmm0[0] = f64ptr(src_rbx+16);/* R3 */ \
	xmm1[0] = f64ptr(src_rbx+24);/* R4 */ \
	xmm5[0] = f64ptr(src_rbx+32);/* R5 */ \
	xmm7[0] = f64ptr(src_rbx+40);/* R6 */ \
	xmm4[0] = f64ptr(src_rbx+48);/* R7 */ \
	xmm6[0] = f64ptr(src_rbx+56);/* R8 */ \
\
	xs8r_mulf(xmm3[0], xmm2[0], xmm0[0], xmm1[0], xmm5[0], xmm7[0], xmm4[0], xmm6[0], (srcreg), (srcreg+8), (srcreg+16), (srcreg+24), (srcreg+32), (srcreg+40), (srcreg+48), (srcreg+56)); \
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
	xmm5[0] = f64ptr(src_rbx+d1);/* R1 */ \
	xmm7[0] = f64ptr(src_rbx+d1+8);/* R2 */ \
	xmm4[0] = f64ptr(src_rbx+d1+16);/* R5 */ \
	xmm6[0] = f64ptr(src_rbx+d1+24);/* R6 */ \
	xmm0[0] = f64ptr(src_rbx+d1+32);/* R3 */ \
	xmm1[0] = f64ptr(src_rbx+d1+40);/* R4 */ \
	xmm3[0] = f64ptr(src_rbx+d1+48);/* R7 */ \
	xmm2[0] = f64ptr(src_rbx+d1+56);/* R8 */ \
\
	xs4c_mulf(xmm5[0], xmm7[0], xmm4[0], xmm6[0], xmm0[0], xmm1[0], xmm3[0], xmm2[0], (srcreg+d1), (srcreg+d1+8), (srcreg+d1+16), (srcreg+d1+24), (srcreg+d1+32), (srcreg+d1+40), (srcreg+d1+48), (srcreg+d1+56)); \
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
\
\
/* Help macros for s2cl_eight_reals_fft_2_final, s2cl_eight_reals_with_square_2, */ \
/* s2cl_eight_reals_with_mult_2, s2cl_eight_reals_with_mulf_2 */ \
\
/* Perform last two fft levels, results returned in */ \
/* R1=xmm3,R2=xmm2,R3=xmm0,R4=xmm1,R5=xmm5,R6=xmm7,R7=xmm4,R8=xmm6 */ \

#define xs8r_fft(r0, r1, r2, r3, r4, r5, r6, r7) { \
	double temp; \
	r5 *= XMM_SQRTHALF1;	/* R6 = R6 * square root of 1/2 */ \
	r7 *= XMM_SQRTHALF1;	/* R8 = R8 * square root of 1/2 */ \
	temp = r0; \
	r0 -= r2;			/* new R3 = R1 - R3 (final R3) */ \
	r2 += temp;			/* new R1 = R1 + R3 */ \
	temp = r1; \
	r1 -= r3;			/* new R4 = R2 - R4 (final R4) */ \
	r3 += temp;			/* new R2 = R2 + R4 */ \
	temp = r5; \
	r5 -= r7;			/* R6 = R6 - R8 (Real part) */ \
	r7 += temp;			/* R8 = R6 + R8 (Imaginary part) */ \
	temp = r2; \
	r2 -= r3;			/* R1 = R1 - R2 (final R2) */ \
	r3 += temp;			/* R2 = R1 + R2 (final R1) */ \
	temp = r4; \
	r4 -= r5;			/* R5 = R5 - R6 (final R7) */ \
	r5 += temp;			/* R6 = R5 + R6 (final R5) */ \
	temp = r6; \
	r6 -= r7;			/* R7 = R7 - R8 (final R8) */ \
	r7 += temp;			/* R8 = R7 + R8 (final R6) */ \
}

/* Perform the multiply step.  Multiply registers with values from memory. */ \
/* Input and output registers are: */ \
/* R1=xmm3,R2=xmm2,R3=xmm4,R4=xmm5,R5=xmm0,R6=xmm1,R7=xmm6,R8=xmm7 */ \

#define xs8r_mulf(r1, r2, r3, r4, r5, r6, r7, r8, m1, m2, m3, m4, m5, m6, m7, m8) { \
	r1 *= f64ptr2(m1,rbp);		/* R11 */ \
	r2 *= f64ptr2(m2,rbp);		/* R22 */ \
	f64ptr(rsi-16) = r1;		/* Save xmm3, Save product of sum of FFT values */ \
	double r1b, r2b; /* XMM_TMP1 = r2 */ /* Save xmm2 */ \
	xs_complex_mult(r3, r4, f64ptr2(m3,rbp), f64ptr2(m4,rbp), r2b, r1b); \
	xs_complex_mult(r5, r6, f64ptr2(m5,rbp), f64ptr2(m6,rbp), r2b, r1b); \
	xs_complex_mult(r7, r8, f64ptr2(m7,rbp), f64ptr2(m8,rbp), r2b, r1b); \
	/* r2 = XMM_TMP1; */		/* Restore xmm2 */ \
	/* r1 = f64ptr(rsi-16); */	/* Restore xmm3 */ \
}
\
/* Perform 2 levels of inverse FFT.  Input registers are: */ \
/* R1=xmm3,R2=xmm2,R3=xmm4,R4=xmm5,R5=xmm0,R6=xmm1,R7=xmm6,R8=xmm7 */ \
/* Output registers are: */ \
/* R1=xmm4,R2=xmm5,R3=xmm2,R4=xmm3,R5=xmm6,R6=xmm0,R7=xmm7,R8=xmm1 */ \

#define xs8r_unfft(r1, r2, r3, r4, r5, r6, r7, r8) \
	double temp; \
	r1 -= r2;			/* R1 = R1 - R2 (new R2) */ \
	r1 *= 0.5;			/* Mul R1 by HALF */ \
	r2 += r1;			/* R2 = R1 + R2 (new R1) */ \
\
	temp = r5; \
	r5 -= r7;			/* R5 = R5 - R7 (new R6) */ \
	r7 += temp;			/* R7 = R5 + R7 (new R5) */ \
\
	temp = r6; \
	r6 -= r8;			/* R6 = R6 - R8 (new R8) */ \
	r8 += temp;			/* R8 = R6 + R8 (new R7) */ \
\
	temp = r6; \
	r6 -= r5;			/* R8 = R8 - R6 */ \
	r5 += temp;			/* R6 = R6 + R8 */ \
	r5 *= XMM_SQRTHALF1;	/* R6 = R6 * square root of 1/2 */ \
	r6 *= XMM_SQRTHALF1;	/* R8 = R8 * square root of 1/2 */ \
\
	temp = r2; \
	r2 -= r3;			/* R1 = R1 - R3 (new R3) */ \
	r3 += temp;			/* R3 = R1 + R3 (new R1) */ \
\
	temp = r1; \
	r1 -= r4;			/* R2 = R2 - R4 (new R4) */ \
	r4 += temp;			/* R4 = R2 + R4 (new R2) */ \
\
\
/* Perform last 2 levels of FFT. */ \
/* ***Optimization - these are sin/cos values special??? */ \

#define xs4c_fft(r1, r2, r3, r4, r5, r6, r7, r8, mem8) { \
	double temp, sine; \
	sine = f64ptr(rdi+32+24);	/* cosine/sine */ \
	r8 = sine * r3;			/* A3 = R3 * cosine/sine	;1-6 */ \
	r8 -= r7;			/* A3 = A3 - I3			;8-11 */ \
	r7 *= sine;	/* B3 = I3 * cosine/sine	;3-8 */ \
	r7 += r3;			/* B3 = B3 + R3			;10-13 */ \
	sine = f64ptr(rdi+0+24);	/* cosine/sine */ \
	r3 = sine * r2;			/* A2 = R2 * cosine/sine	;5-10 */ \
	r3 -= r6;			/* A2 = A2 - I2			;12-15 */ \
	r6 *= sine;	/* B2 = I2 * cosine/sine	;9-14 */ \
	r6 += r2;			/* B2 = B2 + R2			;16-19 */ \
	sine = f64ptr(rdi+64+24);	/* cosine/sine */ \
	temp = mem8; \
	r2 = sine * temp;		/* B4 = I4 * cosine/sine	;11-16 */ \
	r2 += r4;			/* B4 = B4 + R4			;18-21 */ \
	r4 *= sine;	/* A4 = R4 * cosine/sine	;7-12 */ \
	r4 -= temp;		/* A4 = A4 - I4			;14-17 */ \
	sine = f64ptr(rdi+32+8); \
	r8 *= sine;	/* A3 = A3 * sine (new R3)	;13-18 */ \
	r7 *= sine;	/* B3 = B3 * sine (new I3)	;15-20 */ \
	sine = f64ptr(rdi+0+8); \
	r3 *= sine;		/* A2 = A2 * sine (new R2)	;17-22 */ \
	r6 *= sine;		/* B2 = B2 * sine (new I2)	;21-26 */ \
	sine = f64ptr(rdi+64+8); \
	r4 *= sine;	/* A4 = A4 * sine (new R4)	;19-24 */ \
	r2 *= sine;	/* B4 = B4 * sine (new I4)	;23-28 */ \
	temp = r1; \
	r1 -= r8;			/* R1 = R1 - R3 (new R3)	;20-23 */ \
	r8 += temp;			/* R3 = R1 + R3 (new R1)	;24-27 */ \
	temp = r5; \
	r5 -= r7;			/* I1 = I1 - I3 (new I3)	;22-25 */ \
	r7 += temp;			/* I3 = I1 + I3 (new I1)	;28-31 */ \
	temp = r3; \
	r3 -= r4;			/* R2 = R2 - R4 (new R4)	;26-29 */ \
	r4 += temp;			/* R4 = R2 + R4 (new R2)	;34-37 */ \
	temp = r6; \
	r6 -= r2;			/* I2 = I2 - I4 (new I4)	;30-33 */ \
	r2 += temp;			/* I4 = I2 + I4 (new I2)	;36-39 */ \
\
	temp = r5; \
	r5 -= r3;			/* I3 = I3 - R4 (final I4)	;32-35 */ \
	r3 += temp;			/* R4 = I3 + R4 (final I3)	;44-47 */ \
	temp = r1; \
	r1 -= r6;			/* R3 = R3 - I4 (final R3)	;38-41 */ \
	r6 += temp;			/* I4 = R3 + I4 (final R4)	;46-49 */ \
	temp = r8; \
	r8 -= r4;			/* R1 = R1 - R2 (final R2)	;40-43 */ \
	r4 += temp;			/* R2 = R1 + R2 (final R1)	;48-51 */ \
	temp = r7; \
	r7 -= r2;			/* I1 = I1 - I2 (final I2)	;42-45 */ \
	r2 += temp;			/* I2 = I1 + I2 (final I1)	;50-53 */ \
}

/* Perform the first two inverse FFT levels.  Output registers are: */ \
/* R1=xmm2,I1=xmm0,R2=xmm4,I2=xmm6,R3=xmm1,I3=xmm3,R4=xmm5,I4=xmm7 */ \

#define xs4c_unfft(r1, r2, r3, r4, r5, r6, r7, r8) { \
	double tmp1,tmp2,tmp3,tmp4; \
	tmp4 = r1; \
	r1 -= r3;			/* new R2 = R1 - R2 */ \
	r3 += tmp4;			/* new R1 = R1 + R2 */ \
	tmp4 = r6; \
	r6 -= r8;			/* new R4 = I3 - I4 */ \
	r8 += tmp4;			/* new I3 = I3 + I4 */ \
	tmp4 = r2; \
	r2 -= r4;			/* new I2 = I1 - I2 */ \
	r4 += tmp4;			/* new I1 = I1 + I2 */ \
	tmp4 = r7; \
	r7 -= r5;			/* new I4 = R4 - R3 */ \
	r5 += tmp4;			/* new R3 = R3 + R4 */ \
\
	tmp4 = r1; \
	r1 -= r6;			/* R2 = R2 - R4 (new R4) */ \
	r6 += tmp4;			/* R4 = R2 + R4 (new R2) */ \
	tmp1 = r6;			/* Save new R2 */ \
	tmp4 = r2; \
	r2 -= r7;			/* I2 = I2 - I4 (new I4) */ \
	r7 += tmp4;			/* I4 = I2 + I4 (new I2) */ \
	tmp4 = r3; \
	r3 -= r5;			/* R1 = R1 - R3 (new R3) */ \
	r5 += tmp4;			/* R3 = R1 + R3 (new & final R1) */ \
	tmp4 = r4; \
	r4 -= r8;			/* I1 = I1 - I3 (new I3) */ \
	r8 += tmp4;			/* I3 = I1 + I3 (new & final I1) */ \
\
	double *s = (double*)(rdi); \
	tmp2 = r3;		/* Save new R3 */ \
	r2 *= s[9];		/* B4 = new I4 * sine */ \
	tmp3 = r2;		/* Save B4 */ \
	r1 *= s[9];		/* A4 = new R4 * sine */ \
	r6 *= s[3];		/* A2 = new R2 * cosine/sine */ \
	r3 *= s[7];		/* A3 = new R3 * cosine/sine */ \
	r2 *= s[11];	/* C4 = B4 * cosine/sine */ \
	r6 += r7;		/* A2 = A2 + new I2 */ \
	r7 *= s[3];		/* B2 = new I2 * cosine/sine */ \
	r3 += r4;		/* A3 = A3 + new I3 */ \
	r4 *= s[7];		/* B3 = new I3 * cosine/sine */ \
	r2 -= r1;		/* C4 = C4 - A4 (final I4) */ \
	r1 *= s[11];	/* A4 = A4 * cosine/sine */ \
	r7 -= tmp1;		/* B2 = B2 - new R2 */ \
	r6 *= s[1];		/* A2 = A2 * sine (final R2) */ \
	r4 -= tmp2;		/* B3 = B3 - new R3 */ \
	r3 *= s[5];		/* A3 = A3 * sine (final R3) */ \
	r7 *= s[1];		/* B2 = B2 * sine (final I2) */ \
	r1 += tmp3;		/* A4 = B4 + A4 (final R4) */ \
	r4 *= s[5];		/* B3 = B3 * sine (final I3) */ \
}
\
\
/* These macros process two cache lines containing 16 FFT values. */ \
/* This data represents 8 complex numbers. */ \
/* These macros are called frequently and should be optimized. */ \
\

#define s2cl_four_complex_fft_final(srcreg,srcinc,d1) \
	shuffle_load(xmm0,xmm1,xptr(srcreg),xptr(srcreg+32)); /* R1,R2 */ \
	shuffle_load(xmm2,xmm3,xptr(srcreg+d1),xptr(srcreg+d1+32)); /* R3,R4 */ \
	shuffle_load(xmm4,xmm5,xptr(srcreg+16),xptr(srcreg+48)); /* R5,R6 */ \
	shuffle_load(xmm6,xmm7,xptr(srcreg+d1+16),xptr(srcreg+d1+48)); /* R7,R8 */ \
	vec2f64 tmp7 = xmm7; \
	x4c_fft_dst12(xmm0, xmm1, xmm2, xmm3, xmm4, xmm5, xmm6, xmm7, tmp7, rdi, 0, srcreg+srcinc, d1, xptr(srcreg+d1+48), xptr(srcreg+d1)); \
	xptr(srcreg) = xmm3;		/* Save R1 */ \
	xptr(srcreg+16) = xmm1;	/* Save R2 */ \
	xptr(srcreg+32) = xmm7;	/* Save R3 */ \
	xptr(srcreg+48) = xmm6;	/* Save R4 */ \
/*	xstore	[srcreg+d1], xmm0	;; Save R5 */ \
	xptr(srcreg+d1+16) = xmm2;	/* Save R6 */ \
	xptr(srcreg+d1+32) = xmm5;	/* Save R7 */ \
/*	xstore	[srcreg+d1+48], xmm4	;; Save R8 */ \
	srcreg += srcinc; \
\
\

#define s2cl_four_complex_with_square(srcreg,srcinc,d1) { \
	vec2f64 tmp1; \
	shuffle_load(xmm0,xmm1,xptr(srcreg),xptr(srcreg+32)); /* R1,R2 */ \
	shuffle_load(xmm2,xmm3,xptr(srcreg+d1),xptr(srcreg+d1+32)); /* R3,R4 */ \
	shuffle_load(xmm4,xmm5,xptr(srcreg+16),xptr(srcreg+48)); /* R5,R6 */ \
	shuffle_load(xmm6,xmm7,xptr(srcreg+d1+16),xptr(srcreg+d1+48)); /* R7,R8 */ \
	tmp1 = xmm7; \
	x4c_fft_dst2(xmm0, xmm1, xmm2, xmm3, xmm4, xmm5, xmm6, xmm7, tmp1, rdi, 0, srcreg+srcinc, d1, tmp1); \
\
/*	xstore	XMM_TMP1, xmm0 */ \
	xp_complex_square(xmm3, xmm1, xmm0);	/* Square R1, R2 */ \
	xp_complex_square(xmm7, xmm6, xmm0);	/* Square R3, R4 */ \
	xmm0 = tmp1; \
	tmp1 = xmm3; \
	xp_complex_square(xmm0, xmm2, xmm3);	/* Square R5, R6 */ \
	xp_complex_square(xmm5, xmm4, xmm3);	/* Square R7, R8 */ \
	xmm3 = tmp1; \
\
	x4c_unfft(xmm3, xmm1, xmm7, xmm6, xmm0, xmm2, xmm5, xmm4); \
\
	shuffle_store(xptr(srcreg), xptr(srcreg+d1), xmm0, xmm2);	/* Save R1, R2 */ \
	shuffle_store(xptr(srcreg+16), xptr(srcreg+d1+16), xmm7, xmm3); /* Save R3, R4 */ \
	shuffle_store(xptr(srcreg+32), xptr(srcreg+d1+32), xmm4, xmm5); /* Save R5, R6 */ \
	shuffle_store(xptr(srcreg+48), xptr(srcreg+d1+48), xmm6, xmm1); /* Save R7, R8 */ \
\
	srcreg += srcinc; \
}


#define s2cl_four_complex_with_mult(srcreg,srcinc,d1) { \
	shuffle_load(xmm0,xmm1,xptr(srcreg),xptr(srcreg+32)); /* R1,R2 */ \
	shuffle_load(xmm2,xmm3,xptr(srcreg+d1),xptr(srcreg+d1+32)); /* R3,R4 */ \
	shuffle_load(xmm4,xmm5,xptr(srcreg+16),xptr(srcreg+48)); /* R5,R6 */ \
	shuffle_load(xmm6,xmm7,xptr(srcreg+d1+16),xptr(srcreg+d1+48)); /* R7,R8 */ \
	vec2f64 tmp1 = xmm7; \
	x4c_fft_nodst(xmm0, xmm1, xmm2, xmm3, xmm4, xmm5, xmm6, xmm7, tmp1, rdi, 0, srcreg+srcinc, d1); \
\
	xp4c_mulf(xmm3, xmm1, xmm7, xmm6, xmm0, xmm2, xmm5, xmm4, (srcreg), (srcreg+16), (srcreg+32), (srcreg+48), (srcreg+d1), (srcreg+d1+16), (srcreg+d1+32), (srcreg+d1+48)); \
\
	x4c_unfft(xmm3, xmm1, xmm7, xmm6, xmm0, xmm2, xmm5, xmm4); \
\
	shuffle_store(xptr(srcreg), xptr(srcreg+d1), xmm0, xmm2); /* Save R1 and R2 */ \
	shuffle_store(xptr(srcreg+16), xptr(srcreg+d1+16), xmm7, xmm3); /* Save R3 and R4 */ \
	shuffle_store(xptr(srcreg+32), xptr(srcreg+d1+32), xmm4, xmm5); /* Save R5 and R6 */ \
	shuffle_store(xptr(srcreg+48), xptr(srcreg+d1+48), xmm6, xmm1); /* Save R7 and R8 */ \
\
	srcreg += srcinc; \
}


#define s2cl_four_complex_with_mulf(srcreg,srcinc,d1) { \
	uintptr_t src_rbx = srcreg+rbx; \
	xmm3 = xptr(src_rbx);	/* R1 */ \
	xmm7 = xptr(src_rbx+16);	/* R2 */ \
	xmm2 = xptr(src_rbx+32);	/* R3 */ \
	xmm6 = xptr(src_rbx+48);	/* R4 */ \
	xmm0 = xptr(src_rbx+d1);	/* R5 */ \
	xmm1 = xptr(src_rbx+d1+16);/* R6 */ \
	xmm5 = xptr(src_rbx+d1+32);/* R7 */ \
	xmm4 = xptr(src_rbx+d1+48);/* R8 */ \
\
	xp4c_mulf(xmm3, xmm7, xmm2, xmm6, xmm0, xmm1, xmm5, xmm4, (srcreg), (srcreg+16), (srcreg+32), (srcreg+48), (srcreg+d1), (srcreg+d1+16), (srcreg+d1+32), (srcreg+d1+48)); \
\
	x4c_unfft(xmm3, xmm7, xmm2, xmm6, xmm0, xmm1, xmm5, xmm4); \
\
	shuffle_store(xptr(srcreg), xptr(srcreg+d1), xmm0, xmm1); /* Save R1 and R2 */ \
	shuffle_store(xptr(srcreg+16), xptr(srcreg+d1+16), xmm2, xmm3); /* Save R3 and R4 */ \
	shuffle_store(xptr(srcreg+32), xptr(srcreg+d1+32), xmm4, xmm5); /* Save R5 and R6 */ \
	shuffle_store(xptr(srcreg+48), xptr(srcreg+d1+48), xmm6, xmm7); /* Save R7 and R8 */ \
\
	srcreg += srcinc; \
}

/* These macros are used in the complex sections of two pass FFTs. */ \
/* These macros process four cache lines containing 32 FFT values. */ \
/* This data represents 16 complex numbers. */ \
/* These macros are called frequently and should be optimized. */ \
\

#define x4cl_four_complex_fft_final(srcreg,srcinc,d1,d2) \
	x4cl_four_complex_fft(srcreg,srcinc,d1,d2,rdi); \
\
\
/* 234.25 */ \

#define x4cl_four_complex_with_square(srcreg,srcinc,d1,d2) { \
	vec2f64 tmp1, tmp2, tmp3, tmp4, tmp5, tmp6, tmp7, tmp8; \
	x4c_fft_mem_dst2(xptr(srcreg), xptr(srcreg+d1), xptr(srcreg+d2), xptr(srcreg+d2+d1), xptr(srcreg+16), xptr(srcreg+d1+16), xptr(srcreg+d2+16), xptr(srcreg+d2+d1+16), rdi, 0, srcreg+srcinc, d1, tmp5); \
/*	tmp5 = xmm2;*/			/* R5 */ \
	tmp6 = xmm4;			/* R6 */ \
	xp_complex_square(xmm6, xmm0, xmm2);	/* Square R7, R8 */ \
	xp_complex_square(xmm1, xmm3, xmm4);	/* Square R3, R4 */ \
	xp_complex_square(xmm5, xmm7, xmm2);	/* Square R1, R2 */ \
	tmp7 = xmm6; \
	tmp8 = xmm0; \
	tmp3 = xmm1; \
	tmp4 = xmm3; \
	tmp1 = xmm5; \
	tmp2 = xmm7; \
\
	x4c_fft_mem_dst2(xptr(srcreg+32), xptr(srcreg+d1+32), xptr(srcreg+d2+32), xptr(srcreg+d2+d1+32), xptr(srcreg+48), xptr(srcreg+d1+48), xptr(srcreg+d2+48), xptr(srcreg+d2+d1+48), rdi,XMM_SCD,srcreg+srcinc+d2,d1, xptr(srcreg+d2+d1)); \
/*	xptr(srcreg+d2+d1) = xmm2;*/		/* R5 */ \
	xptr(srcreg+d2+d1+16) = xmm4;		/* R6 */ \
	xp_complex_square(xmm6, xmm0, xmm2);	/* Square R7, R8 */ \
	xp_complex_square(xmm1, xmm3, xmm4);	/* Square R3, R4 */ \
	xp_complex_square(xmm5, xmm7, xmm2);	/* Square R1, R2 */ \
	xptr(srcreg+d2+d1+32) = xmm6;		/* R7 */ \
	xptr(srcreg+d2+d1+48) = xmm0;		/* R8 */ \
	xptr(srcreg+d2+32) = xmm1;		/* R3 */ \
	xptr(srcreg+d2+48) = xmm3;		/* R4 */ \
	xptr(srcreg+d2) = xmm5;		/* R1 */ \
	xptr(srcreg+d2+16) = xmm7;		/* R2 */ \
\
	xmm4 = tmp5;		/* R5 */ \
	xmm5 = tmp6;		/* R6 */ \
	xp_complex_square(xmm4, xmm5, xmm6); /* Square R5, R6 */ \
	xmm0 = tmp1;		/* R1 */ \
	xmm2 = tmp3;		/* R3 */ \
	xmm1 = tmp2;		/* R2 */ \
	xmm3 = tmp4;		/* R4 */ \
	best_x4c_unfft(xmm0, xmm1, xmm2, xmm3, xmm4, xmm5, xmm6, xmm7, tmp7, tmp8, xptr(srcreg), xptr(srcreg+32), rdi, 0); \
	xptr(srcreg+16) = xmm5;	/* Save R3 */ \
	xptr(srcreg+48) = xmm6;	/* Save R7 */ \
	xptr(srcreg+d1) = xmm1;	/* Save R2 */ \
	xptr(srcreg+d1+16) = xmm4;	/* Save R4 */ \
	xptr(srcreg+d1+32) = xmm0;	/* Save R6 */ \
	xptr(srcreg+d1+48) = xmm7;	/* Save R8 */ \
\
	xmm4 = xptr(srcreg+d2+d1);	/* R5 */ \
	xmm5 = xptr(srcreg+d2+d1+16);	/* R6 */ \
	xp_complex_square(xmm4, xmm5, xmm6); /* Square R5, R6 */ \
	xmm0 = xptr(srcreg+d2);	/* R1 */ \
	xmm2 = xptr(srcreg+d2+32);	/* R3 */ \
	xmm1 = xptr(srcreg+d2+16);	/* R2 */ \
	xmm3 = xptr(srcreg+d2+48);	/* R4 */ \
	best_x4c_unfft(xmm0, xmm1, xmm2, xmm3, xmm4, xmm5, xmm6, xmm7, xptr(srcreg+d2+d1+32), xptr(srcreg+d2+d1+48), xptr(srcreg+d2), xptr(srcreg+d2+32), rdi, XMM_SCD); \
	xptr(srcreg+d2+16) = xmm5;	/* Save R3 */ \
	xptr(srcreg+d2+48) = xmm6;	/* Save R7 */ \
	xptr(srcreg+d2+d1) = xmm1;	/* Save R2 */ \
	xptr(srcreg+d2+d1+16) = xmm4;	/* Save R4 */ \
	xptr(srcreg+d2+d1+32) = xmm0;	/* Save R6 */ \
	xptr(srcreg+d2+d1+48) = xmm7;	/* Save R8 */ \
\
	srcreg += srcinc; \
}


#define x4cl_four_complex_with_mult(srcreg,srcinc,d1,d2) \
	x4c_fft_mem_nodst(xptr(srcreg), xptr(srcreg+d1), xptr(srcreg+d2), xptr(srcreg+d2+d1), xptr(srcreg+16), xptr(srcreg+d1+16), xptr(srcreg+d2+16), xptr(srcreg+d2+d1+16), rdi,0,srcreg+srcinc,d1); \
\
	xp4c_mulf(xmm5, xmm7, xmm1, xmm3, xmm2, xmm4, xmm6, xmm0, (srcreg), (srcreg+16), (srcreg+32), (srcreg+48), (srcreg+d1), (srcreg+d1+16), (srcreg+d1+32), (srcreg+d1+48)); \
\
	new_x4c_unfft(xmm5, xmm7, xmm1, xmm3, xmm2, xmm4, xmm6, xmm0, xptr(srcreg), 0); \
	xmm1 = xptr(srcreg+32);	/* R1 */ \
	xptr(srcreg+d1) = xmm7;	/* Save R2 */ \
	xmm7 = xptr(srcreg+48);	/* R5 */ \
	xptr(srcreg+16) = xmm4;	/* Save R3 */ \
	xmm4 = xptr(srcreg+d1+32);	/* R2 */ \
	xptr(srcreg+d1+16) = xmm5;	/* Save R4 */ \
	xmm5 = xptr(srcreg+d1+48);	/* R6 */ \
	xptr(srcreg+32) = xmm0;	/* Save R5 */ \
	xmm0 = xptr(srcreg+d2+32);	/* R3 */ \
	xptr(srcreg+d1+32) = xmm6;	/* Save R6 */ \
	xmm6 = xptr(srcreg+d2+48);	/* R7 */ \
	xptr(srcreg+48) = xmm3;	/* Save R7 */ \
	xmm3 = xptr(srcreg+d2+d1+32);	/* R4 */ \
	xptr(srcreg+d1+48) = xmm2;	/* Save R8 */ \
\
	x4c_fft_nodst(xmm1, xmm4, xmm0, xmm3, xmm7, xmm5, xmm6, xmm2, xptr(srcreg+d2+d1+48), rdi, XMM_SCD, srcreg+srcinc+d2, d1); \
\
	xp4c_mulf(xmm3, xmm4, xmm2, xmm6, xmm1, xmm0, xmm5, xmm7, (srcreg+d2), (srcreg+d2+16), (srcreg+d2+32), (srcreg+d2+48), (srcreg+d2+d1), (srcreg+d2+d1+16), (srcreg+d2+d1+32), (srcreg+d2+d1+48)); \
\
	new_x4c_unfft(xmm3, xmm4, xmm2, xmm6, xmm1, xmm0, xmm5, xmm7, xptr(srcreg+d2), XMM_SCD); \
\
	xptr(srcreg+d2+d1) = xmm4;	/* Save R2 */ \
	xptr(srcreg+d2+16) = xmm0;	/* Save R3 */ \
	xptr(srcreg+d2+d1+16) = xmm3;	/* Save R4 */ \
	xptr(srcreg+d2+32) = xmm7;	/* Save R5 */ \
	xptr(srcreg+d2+d1+32) = xmm5;	/* Save R6 */ \
	xptr(srcreg+d2+48) = xmm6;	/* Save R7 */ \
	xptr(srcreg+d2+d1+48) = xmm1;	/* Save R8 */ \
	srcreg += srcinc; \
\
\

#define x4cl_four_complex_with_mulf(srcreg,srcinc,d1,d2) { \
	uintptr_t src_rbx = srcreg+rbx; \
	xmm5 = xptr(src_rbx);	/* R1 */ \
	xmm7 = xptr(src_rbx+16);	/* R2 */ \
	xmm1 = xptr(src_rbx+32);	/* R3 */ \
	xmm3 = xptr(src_rbx+48);	/* R4 */ \
	xmm2 = xptr(src_rbx+d1);	/* R5 */ \
	xmm4 = xptr(src_rbx+d1+16);/* R6 */ \
	xmm6 = xptr(src_rbx+d1+32);/* R7 */ \
	xmm0 = xptr(src_rbx+d1+48);/* R8 */ \
\
	xp4c_mulf(xmm5, xmm7, xmm1, xmm3, xmm2, xmm4, xmm6, xmm0, (srcreg), (srcreg+16), (srcreg+32), (srcreg+48), (srcreg+d1), (srcreg+d1+16), (srcreg+d1+32), (srcreg+d1+48)); \
\
	new_x4c_unfft(xmm5, xmm7, xmm1, xmm3, xmm2, xmm4, xmm6, xmm0, xptr(srcreg), 0); \
	xptr(srcreg+d1) = xmm7;	/* Save R2 */ \
	xptr(srcreg+16) = xmm4;	/* Save R3 */ \
	xptr(srcreg+d1+16) = xmm5;	/* Save R4 */ \
	xptr(srcreg+32) = xmm0;	/* Save R5 */ \
	xptr(srcreg+d1+32) = xmm6;	/* Save R6 */ \
	xptr(srcreg+48) = xmm3;	/* Save R7 */ \
	xptr(srcreg+d1+48) = xmm2;	/* Save R8 */ \
\
	xmm3 = xptr(src_rbx+d2);	/* R1 */ \
	xmm4 = xptr(src_rbx+d2+16);/* R2 */ \
	xmm2 = xptr(src_rbx+d2+32);/* R3 */ \
	xmm6 = xptr(src_rbx+d2+48);/* R4 */ \
	xmm1 = xptr(src_rbx+d2+d1);/* R5 */ \
	xmm0 = xptr(src_rbx+d2+d1+16);/* R6 */ \
	xmm5 = xptr(src_rbx+d2+d1+32);/* R7 */ \
	xmm7 = xptr(src_rbx+d2+d1+48);/* R8 */ \
\
	xp4c_mulf(xmm3, xmm4, xmm2, xmm6, xmm1, xmm0, xmm5, xmm7, (srcreg+d2), (srcreg+d2+16), (srcreg+d2+32), (srcreg+d2+48), (srcreg+d2+d1), (srcreg+d2+d1+16), (srcreg+d2+d1+32), (srcreg+d2+d1+48)); \
\
	new_x4c_unfft(xmm3, xmm4, xmm2, xmm6, xmm1, xmm0, xmm5, xmm7, xptr(srcreg+d2), XMM_SCD); \
\
	xptr(srcreg+d2+d1) = xmm4;	/* Save R2 */ \
	xptr(srcreg+d2+16) = xmm0;	/* Save R3 */ \
	xptr(srcreg+d2+d1+16) = xmm3;	/* Save R4 */ \
	xptr(srcreg+d2+32) = xmm7;	/* Save R5 */ \
	xptr(srcreg+d2+d1+32) = xmm5;	/* Save R6 */ \
	xptr(srcreg+d2+48) = xmm6;	/* Save R7 */ \
	xptr(srcreg+d2+d1+48) = xmm1;	/* Save R8 */ \
	srcreg += srcinc; \
}


/******************************************************************** */ \
/* These macros are used in the premultiplier step of two pass FFTs */ \
/******************************************************************** */ \
\
/* 84 clocks */ \

#define s2cl_four_complex_gpm_fft(srcreg,srcinc,d1) { \
	uintptr_t src_rbx = srcreg+rbx; \
	shuffle_load(xmm0, xmm1, xptr(src_rbx), xptr(src_rbx+32)); /* R1,R3 */ \
	xmm6 = xmm0;		/* Save R1 */ \
	xmm0 *= xptr(rdi+16);		/* A1 = R1 * premul_real/premul_imag */ \
	xmm7 = xmm1;		/* Save R3 */ \
	xmm1 *= xptr(rdi+80);		/* A3 = R3 * premul_real/premul_imag */ \
\
	shuffle_load(xmm2, xmm3, xptr(src_rbx+16), xptr(src_rbx+48)); /* R5,R7 */ \
	xprefetch(u8ptr(src_rbx+srcinc)); \
	xmm0 -= xmm2;		/* A1 = A1 - I1 */ \
	xmm2 *= xptr(rdi+16);		/* B1 = I1 * premul_real/premul_imag */ \
	xmm1 -= xmm3;		/* A3 = A3 - I3 */ \
	xmm3 *= xptr(rdi+80);		/* B3 = I3 * premul_real/premul_imag */ \
	xmm2 += xmm6;		/* B1 = B1 + R1 */ \
	xmm0 *= xptr(rdi);		/* A1 = A1 * premul_imag (new R1) */ \
	xmm3 += xmm7;		/* B3 = B3 + R3 */ \
	xmm2 *= xptr(rdi);		/* B1 = B1 * premul_imag (new I1) */ \
\
	shuffle_load(xmm4, xmm5, xptr(src_rbx+d1), xptr(src_rbx+d1+32)); /* R2,R4 */ \
	xprefetch(u8ptr(src_rbx+srcinc+d1)); \
	xptr(srcreg) = xmm0;		/* Save new R1 */ \
	xmm0 = xmm4;		/* Save R2 */ \
	xmm4 *= xptr(rdi+48);		/* A2 = R2 * premul_real/premul_imag */ \
	xmm1 *= xptr(rdi+64);		/* A3 = A3 * premul_imag (new R3) */ \
	xmm3 *= xptr(rdi+64);		/* B3 = B3 * premul_imag (new I3) */ \
\
	shuffle_load(xmm6,xmm7,xptr(src_rbx+d1+16),xptr(src_rbx+d1+48)); /* R6,R8 */ \
	xprefetchw(u8ptr(srcreg+srcinc)); \
	xmm4 -= xmm6;		/* A2 = A2 - I2 */ \
	xmm6 *= xptr(rdi+48);		/* B2 = I2 * premul_real/premul_imag */ \
	xmm6 += xmm0;		/* B2 = B2 + R2 */ \
	xmm0 = xmm5;		/* Save R4 */ \
	xmm5 *= xptr(rdi+112);		/* A4 = R4 * premul_real/premul_imag */ \
	xmm5 -= xmm7;		/* A4 = A4 - I4 */ \
	xmm7 *= xptr(rdi+112);		/* B4 = I4 * premul_real/premul_imag */ \
	xmm4 *= xptr(rdi+32);		/* A2 = A2 * premul_imag (new R2) */ \
	xmm7 += xmm0;		/* B4 = B4 + R4 */ \
	xmm6 *= xptr(rdi+32);		/* B2 = B2 * premul_imag (new I2) */ \
	xmm5 *= xptr(rdi+96);		/* A4 = A4 * premul_imag (new R4) */ \
	xmm7 *= xptr(rdi+96);		/* B4 = B4 * premul_imag (new I4) */ \
\
	xprefetchw(u8ptr(srcreg+srcinc+d1)); \
	 xmm0 = xmm2; \
	 xmm2 -= xmm3;		/* I1 = I1 - I3 (new I3) */ \
	 xmm3 += xmm0;		/* I3 = I1 + I3 (new I1) */ \
	 xmm0 = xmm4; \
	 xmm4 -= xmm5;		/* R2 = R2 - R4 (new R4) */ \
	 xmm5 += xmm0;		/* R4 = R2 + R4 (new R2) */ \
	 xmm0 = xmm6; \
	 xmm6 -= xmm7;		/* I2 = I2 - I4 (new I4) */ \
	 xmm7 += xmm0;		/* I4 = I2 + I4 (new I2) */ \
	xmm0 = xmm2; \
	xmm2 -= xmm4;		/* I3 = I3 - R4 (final I4) */ \
	xptr(srcreg+d1+48) = xmm2; \
	xmm4 += xmm0;		/* R4 = I3 + R4 (final I3) */ \
	 xmm0 = xptr(srcreg);		/* Reload new R1 */ \
	 xmm0 -= xmm1;		/* R1 = R1 - R3 (new R3) */ \
	 xmm1 += xptr(srcreg);		/* R3 = R1 + R3 (new R1) */ \
	xmm2 = xmm3; \
	xmm3 -= xmm7;		/* I1 = I1 - I2 (final I2) */ \
	xptr(srcreg+48) = xmm3; \
	xmm7 += xmm2;		/* I2 = I1 + I2 (final I1) */ \
	xmm2 = xmm0; \
	xmm0 -= xmm6;		/* R3 = R3 - I4 (final R3) */ \
	xmm3 = xmm1; \
	xmm1 -= xmm5;		/* R1 = R1 - R2 (final R2) */ \
	xmm6 += xmm2;		/* I4 = R3 + I4 (final R4) */ \
	xmm5 += xmm3;		/* R2 = R1 + R2 (final R1) */ \
	xptr(srcreg+d1) = xmm0; \
	xptr(srcreg+d1+16) = xmm4; \
	xptr(srcreg+d1+32) = xmm6; \
	xptr(srcreg) = xmm5; \
	xptr(srcreg+16) = xmm7; \
	xptr(srcreg+32) = xmm1; \
	srcreg += srcinc; \
}

/* This code is identical to the above except it is used in the single pass negacyclic FFTs and has a somewhat different memory layout. */ \

#define s2cl_four_complex_first_fft(srcreg,srcinc,d1) { \
	uintptr_t src_rbx = srcreg+rbx; \
	vec2f64 temp1, s1; \
	shuffle_load(xmm0,xmm1,xptr(src_rbx),xptr(src_rbx+16)); /* R1,R3 */ \
	shuffle_load(xmm2,xmm3,xptr(src_rbx+32),xptr(src_rbx+48)); /* R5,R7 */ \
	shuffle_load(xmm4,xmm5,xptr(src_rbx+d1),xptr(src_rbx+d1+16)); /* R2,R4 */ \
	shuffle_load(xmm6,xmm7,xptr(src_rbx+d1+32),xptr(src_rbx+d1+48)); /* R6,R8 */ \
\
	temp1 = xmm5;		/* Save R4 */ \
\
	xmm5 = xmm0;		/* Save R1 */ \
	s1 = xptr(rdi+16); \
	xmm0 *= s1;		/* A1 = R1 * premul_real/premul_imag */ \
	xmm0 -= xmm2;		/* A1 = A1 - I1 */ \
	xmm2 *= s1;		/* B1 = I1 * premul_real/premul_imag */ \
	xmm2 += xmm5;		/* B1 = B1 + R1 */ \
\
	xmm5 = xmm1;		/* Save R3 */ \
	s1 = xptr(rdi+80); \
	xmm1 *= s1;		/* A3 = R3 * premul_real/premul_imag */ \
	xmm1 -= xmm3;		/* A3 = A3 - I3 */ \
	xmm3 *= s1;		/* B3 = I3 * premul_real/premul_imag */ \
	xmm3 += xmm5;		/* B3 = B3 + R3 */ \
	s1 = xptr(rdi); \
	xmm0 *= s1;		/* A1 = A1 * premul_imag (new R1) */ \
	xmm2 *= s1;		/* B1 = B1 * premul_imag (new I1) */ \
\
	xmm5 = xmm4;		/* Save R2 */ \
	s1 = xptr(rdi+48); \
	xmm4 *= s1;		/* A2 = R2 * premul_real/premul_imag */ \
	xmm4 -= xmm6;		/* A2 = A2 - I2 */ \
	xmm6 *= s1;		/* B2 = I2 * premul_real/premul_imag */ \
	xmm6 += xmm5;		/* B2 = B2 + R2 */ \
	s1 = xptr(rdi+64); \
	xmm1 *= s1;		/* A3 = A3 * premul_imag (new R3) */ \
	xmm3 *= s1;		/* B3 = B3 * premul_imag (new I3) */ \
\
	xmm5 = temp1;		/* Reload R4 */ \
	s1 = xptr(rdi+112); \
	xmm5 *= s1;		/* A4 = R4 * premul_real/premul_imag */ \
	xmm5 -= xmm7;		/* A4 = A4 - I4 */ \
	xmm7 *= s1;		/* B4 = I4 * premul_real/premul_imag */ \
	xmm7 += temp1;		/* B4 = B4 + R4 */ \
	s1 = xptr(rdi+32); \
	xmm4 *= s1;		/* A2 = A2 * premul_imag (new R2) */ \
	xmm6 *= s1;		/* B2 = B2 * premul_imag (new I2) */ \
	s1 = xptr(rdi+96); \
	xmm5 *= s1;		/* A4 = A4 * premul_imag (new R4) */ \
	xmm7 *= s1;		/* B4 = B4 * premul_imag (new I4) */ \
\
	temp1 = xmm0; \
	xmm0 -= xmm1;		/* R1 = R1 - R3 (new R3) */ \
	xmm1 += temp1;		/* R3 = R1 + R3 (new R1) */ \
	temp1 = xmm2; \
	xmm2 -= xmm3;		/* I1 = I1 - I3 (new I3) */ \
	xmm3 += temp1;		/* I3 = I1 + I3 (new I1) */ \
	temp1 = xmm4; \
	xmm4 -= xmm5;		/* R2 = R2 - R4 (new R4) */ \
	xmm5 += temp1;		/* R4 = R2 + R4 (new R2) */ \
	temp1 = xmm6; \
	xmm6 -= xmm7;		/* I2 = I2 - I4 (new I4) */ \
	xmm7 += temp1;		/* I4 = I2 + I4 (new I2) */ \
\
	temp1 = xmm0; \
	xmm0 -= xmm6;		/* R3 = R3 - I4 (new R3) */ \
	xmm6 += temp1;		/* I4 = R3 + I4 (new R4) */ \
	temp1 = xmm2; \
	xmm2 -= xmm4;		/* I3 = I3 - R4 (new I4) */ \
	xmm4 += temp1;		/* R4 = I3 + R4 (new I3) */ \
	temp1 = xmm1; \
	xmm1 -= xmm5;		/* R1 = R1 - R2 (new R2) */ \
	xmm5 += temp1;		/* R2 = R1 + R2 (new R1) */ \
	temp1 = xmm3; \
	xmm3 -= xmm7;		/* I1 = I1 - I2 (new I2) */ \
	xmm7 += temp1;		/* I2 = I1 + I2 (new I1) */ \
\
	xptr(srcreg+d1) = xmm0; \
	xptr(srcreg+d1+16) = xmm4; \
	xptr(srcreg+d1+32) = xmm6; \
	xptr(srcreg+d1+48) = xmm2; \
	xptr(srcreg) = xmm5; \
	xptr(srcreg+16) = xmm7; \
	xptr(srcreg+32) = xmm1; \
	xptr(srcreg+48) = xmm3; \
	srcreg += srcinc; \
}
/* This code is identical to the above except it is used in the two pass negacyclic FFTs and has a somewhat different memory layout. */ \

#define x2cl_four_complex_first_fft(srcreg,srcinc,d1) \
	x2cl_four_complex_first_fft_cmn(srcreg,srcinc,d1,rbx); \
\

#define x2cl_four_complex_first_fft_scratch(srcreg,srcinc,d1) \
	x2cl_four_complex_first_fft_cmn(srcreg,srcinc,d1,0); \
\

#define x2cl_four_complex_first_fft_cmn(srcreg,srcinc,d1,off) { \
	vec2f64 temp1, s1; \
	xmm0 = xptr(srcreg+off);		/* R1 */ \
	xmm1 = xptr(srcreg+off+16);		/* R3 */ \
	xmm2 = xptr(srcreg+off+32);		/* R5 */ \
	xmm3 = xptr(srcreg+off+48);		/* R7 */ \
	xmm4 = xptr(srcreg+off+d1);		/* R2 */ \
	xmm5 = xptr(srcreg+off+d1+16);	/* R4 */ \
	xmm6 = xptr(srcreg+off+d1+32);	/* R6 */ \
	xmm7 = xptr(srcreg+off+d1+48);	/* R8 */ \
\
	temp1 = xmm5;		/* Save R4 */ \
\
	xmm5 = xmm0;		/* Save R1 */ \
	s1 = xptr(rdi+16); \
	xmm0 *= s1;		/* A1 = R1 * premul_real/premul_imag */ \
	xmm0 -= xmm2;		/* A1 = A1 - I1 */ \
	xmm2 *= s1;		/* B1 = I1 * premul_real/premul_imag */ \
	xmm2 += xmm5;		/* B1 = B1 + R1 */ \
\
	xmm5 = xmm1;		/* Save R3 */ \
	s1 = xptr(rdi+80); \
	xmm1 *= s1;		/* A3 = R3 * premul_real/premul_imag */ \
	xmm1 -= xmm3;		/* A3 = A3 - I3 */ \
	xmm3 *= s1;		/* B3 = I3 * premul_real/premul_imag */ \
	xmm3 += xmm5;		/* B3 = B3 + R3 */ \
	s1 = xptr(rdi); \
	xmm0 *= s1;		/* A1 = A1 * premul_imag (new R1) */ \
	xmm2 *= s1;		/* B1 = B1 * premul_imag (new I1) */ \
\
	xmm5 = xmm4;		/* Save R2 */ \
	s1 = xptr(rdi+48); \
	xmm4 *= s1;		/* A2 = R2 * premul_real/premul_imag */ \
	xmm4 -= xmm6;		/* A2 = A2 - I2 */ \
	xmm6 *= s1;		/* B2 = I2 * premul_real/premul_imag */ \
	xmm6 += xmm5;		/* B2 = B2 + R2 */ \
	s1 = xptr(rdi+64); \
	xmm1 *= s1;		/* A3 = A3 * premul_imag (new R3) */ \
	xmm3 *= s1;		/* B3 = B3 * premul_imag (new I3) */ \
\
	xmm5 = temp1;		/* Reload R4 */ \
	s1 = xptr(rdi+112); \
	xmm5 *= s1;		/* A4 = R4 * premul_real/premul_imag */ \
	xmm5 -= xmm7;		/* A4 = A4 - I4 */ \
	xmm7 *= s1;		/* B4 = I4 * premul_real/premul_imag */ \
	xmm7 += temp1;		/* B4 = B4 + R4 */ \
	s1 = xptr(rdi+32); \
	xmm4 *= s1;		/* A2 = A2 * premul_imag (new R2) */ \
	xmm6 *= s1;		/* B2 = B2 * premul_imag (new I2) */ \
	s1 = xptr(rdi+96); \
	xmm5 *= s1;		/* A4 = A4 * premul_imag (new R4) */ \
	xmm7 *= s1;		/* B4 = B4 * premul_imag (new I4) */ \
\
	temp1 = xmm0; \
	xmm0 -= xmm1;		/* R1 = R1 - R3 (new R3) */ \
	xmm1 += temp1;		/* R3 = R1 + R3 (new R1) */ \
	temp1 = xmm2; \
	xmm2 -= xmm3;		/* I1 = I1 - I3 (new I3) */ \
	xmm3 += temp1;		/* I3 = I1 + I3 (new I1) */ \
	temp1 = xmm4; \
	xmm4 -= xmm5;		/* R2 = R2 - R4 (new R4) */ \
	xmm5 += temp1;		/* R4 = R2 + R4 (new R2) */ \
	temp1 = xmm6; \
	xmm6 -= xmm7;		/* I2 = I2 - I4 (new I4) */ \
	xmm7 += temp1;		/* I4 = I2 + I4 (new I2) */ \
\
	temp1 = xmm0; \
	xmm0 -= xmm6;		/* R3 = R3 - I4 (new R3) */ \
	xmm6 += temp1;		/* I4 = R3 + I4 (new R4) */ \
	temp1 = xmm2; \
	xmm2 -= xmm4;		/* I3 = I3 - R4 (new I4) */ \
	xmm4 += temp1;		/* R4 = I3 + R4 (new I3) */ \
	temp1 = xmm1; \
	xmm1 -= xmm5;		/* R1 = R1 - R2 (new R2) */ \
	xmm5 += temp1;		/* R2 = R1 + R2 (new R1) */ \
	temp1 = xmm3; \
	xmm3 -= xmm7;		/* I1 = I1 - I2 (new I2) */ \
	xmm7 += temp1;		/* I2 = I1 + I2 (new I1) */ \
\
	xptr(srcreg+d1) = xmm0; \
	xptr(srcreg+d1+16) = xmm4; \
	xptr(srcreg+d1+32) = xmm6; \
	xptr(srcreg+d1+48) = xmm2; \
	xptr(srcreg) = xmm5; \
	xptr(srcreg+16) = xmm7; \
	xptr(srcreg+32) = xmm1; \
	xptr(srcreg+48) = xmm3; \
	srcreg += srcinc; \
}


/* 126.35 clocks */

#define s4cl_four_complex_gpm_unfft(srcreg,srcinc,d1,d2,off) \
	xmm0 = xptr(srcreg);		/* R1 */ \
	xmm1 = xptr(srcreg+32);	/* R2 */ \
	xmm4 = xptr(srcreg+d2);	/* R5 */ \
	xmm5 = xptr(srcreg+d2+32);	/* R6 */ \
	x4gpm_unfft(xmm0,xmm1,xmm2,xmm3,xmm4,xmm5,xmm6,xmm7,xptr(srcreg+d1),xptr(srcreg+d1+32),xptr(srcreg+d2+d1),xptr(srcreg+d2+d1+32),xptr(srcreg),xptr(srcreg+32),0,srcreg+srcinc,d1); \
	shuffle_store(xptr(srcreg), xptr(srcreg+32), xmm5, xmm4); /* Save R1 and R3 */ \
	xmm5 = xptr(srcreg+16);	/* R1 */ \
	xmm4 = xptr(srcreg+48);	/* R2 */ \
	shuffle_store(xptr(srcreg+16), xptr(srcreg+48), xmm6, xmm7); /* Save R5 and R7 */ \
	shuffle_store(xptr(srcreg+d2), xptr(srcreg+d2+32), xmm0, xmm2); /* Save R2 and R4 */ \
	xmm0 = xptr(srcreg+d2+16);	/* R5 */ \
	xmm2 = xptr(srcreg+d2+48);	/* R6 */ \
	shuffle_store(xptr(srcreg+d2+16), xptr(srcreg+d2+48), xmm1, xmm3); /* Save R6, R8 */ \
	x4gpm_unfft(xmm5,xmm4,xmm6,xmm7,xmm0,xmm2,xmm1,xmm3,xptr(srcreg+d1+16),xptr(srcreg+d1+48),xptr(srcreg+d2+d1+16),xptr(srcreg+d2+d1+48),xptr(srcreg+d1),xptr(srcreg+d1+32),off,srcreg+srcinc+d2,d1); \
	shuffle_store(xptr(srcreg+d1), xptr(srcreg+d1+32), xmm2, xmm0); /* Save R1 and R3 */ \
	shuffle_store(xptr(srcreg+d1+16), xptr(srcreg+d1+48), xmm1, xmm3); /* Save R5, R7 */ \
	shuffle_store(xptr(srcreg+d2+d1), xptr(srcreg+d2+d1+32), xmm5, xmm6); /* R2, R4 */ \
	shuffle_store(xptr(srcreg+d2+d1+16), xptr(srcreg+d2+d1+48), xmm4, xmm7); /* R6,R8 */ \
\
	srcreg += srcinc; \
\
/* Similar to above but used in negacyclic FFT case. */ \

#define x4cl_four_complex_last_unfft(srcreg,srcinc,d1,d2,off) \
	xmm0 = xptr(srcreg);		/* R1 */ \
	xmm1 = xptr(srcreg+32);	/* R2 */ \
	xmm4 = xptr(srcreg+d2);	/* R5 */ \
	xmm5 = xptr(srcreg+d2+32);	/* R6 */ \
	x4gpm_unfft(xmm0,xmm1,xmm2,xmm3,xmm4,xmm5,xmm6,xmm7,xptr(srcreg+d1),xptr(srcreg+d1+32),xptr(srcreg+d1+d2),xptr(srcreg+d1+d2+32),xptr(srcreg),xptr(srcreg+32),0,srcreg+srcinc,d1); \
	xptr(srcreg) = xmm5;		/* Save R1 */ \
	xptr(srcreg+32) = xmm6;	/* Save R5 */ \
	xmm5 = xptr(srcreg+16);	/* R1 */ \
	xmm6 = xptr(srcreg+48);	/* R2 */ \
	xptr(srcreg+16) = xmm4;	/* Save R3 */ \
	xptr(srcreg+48) = xmm7;	/* Save R7 */ \
	xptr(srcreg+d2) = xmm0;	/* Save R2 */ \
	xptr(srcreg+d2+32) = xmm1;	/* Save R6 */ \
	xmm4 = xptr(srcreg+d2+16);	/* R5 */ \
	xmm0 = xptr(srcreg+d2+48);	/* R6 */ \
	xptr(srcreg+d2+16) = xmm2;	/* Save R4 */ \
	xptr(srcreg+d2+48) = xmm3;	/* Save R8 */ \
	x4gpm_unfft(xmm5,xmm6,xmm2,xmm3,xmm4,xmm0,xmm1,xmm7,xptr(srcreg+d1+16),xptr(srcreg+d1+48),xptr(srcreg+d1+d2+16),xptr(srcreg+d1+d2+48),xptr(srcreg+d1),xptr(srcreg+d1+32),off,srcreg+srcinc+d2,d1); \
	xptr(srcreg+d1) = xmm0;	/* Save R1 */ \
	xptr(srcreg+d1+16) = xmm4;	/* Save R3 */ \
	xptr(srcreg+d1+32) = xmm1;	/* Save R5 */ \
	xptr(srcreg+d1+48) = xmm7;	/* Save R7 */ \
	xptr(srcreg+d2+d1) = xmm5;	/* Save R2 */ \
	xptr(srcreg+d2+d1+16) = xmm2;	/* Save R4 */ \
	xptr(srcreg+d2+d1+32) = xmm6;	/* Save R6 */ \
	xptr(srcreg+d2+d1+48) = xmm3;	/* Save R8 */ \
	srcreg += srcinc; \
\

#define x4gpm_unfft(r1,r2,r3,r4,r5,r6,r7,r8,mem3,mem4,mem7,mem8,dest1,dest2,off,pre1,pre2) \
	x4c_unfft4_cmn(r1,r2,r3,r4,r5,r6,r7,r8,mem3,mem4,mem7,mem8,dest1,dest2,off,off+32,off+64,off+96,pre1,pre2); \
\
\
/* New cpm code that uses more premultipliers but saves lots of xsincos_complex */ \
/* data. */ \

#define x4cl_four_complex_cpm_fft(srcreg,srcinc,d1,d2,off) \
	xmm0 = xptr(srcreg);		/* R1 */ \
	xmm4 = xptr(srcreg+16);	/* R5 */ \
	xmm1 = xptr(srcreg+d1);	/* R2 */ \
	xmm5 = xptr(srcreg+d1+16);	/* R6 */ \
	xmm2 = xptr(srcreg+d2);	/* R3 */ \
	xmm6 = xptr(srcreg+d2+16);	/* R7 */ \
\
	x4c_fft4_cmn(xmm0,xmm1,xmm2,xmm3,xmm4,xmm5,xmm6,xmm7,xptr(srcreg+d2+d1),xptr(srcreg+d2+d1+16),xptr(srcreg+d1),0,32,64,96,srcreg+srcinc,d1); \
\
	xptr(srcreg) = xmm1;		/* Save R1 */ \
	xptr(srcreg+16) = xmm5;	/* Save I1 */ \
	xmm1 = xptr(srcreg+32);	/* R1 */ \
	xmm5 = xptr(srcreg+48);	/* R5 */ \
	xptr(srcreg+32) = xmm0;	/* Save R2 */ \
	xptr(srcreg+48) = xmm4;	/* Save I2 */ \
	xmm0 = xptr(srcreg+d1+32);	/* R2 */ \
	xmm4 = xptr(srcreg+d1+48);	/* R6 */ \
	xptr(srcreg+d1+16) = xmm2;	/* Save I3 */ \
	xptr(srcreg+d1+32) = xmm6;	/* Save R4 */ \
	xptr(srcreg+d1+48) = xmm7;	/* Save I4 */ \
	xmm2 = xptr(srcreg+d2+32);	/* R3 */ \
	xmm3 = xptr(srcreg+d2+48);	/* R7 */ \
\
	x4c_fft4_cmn(xmm1,xmm0,xmm2,xmm6,xmm5,xmm4,xmm3,xmm7,xptr(srcreg+d2+d1+32),xptr(srcreg+d2+d1+48),xptr(srcreg+d2+d1),off,off+32,off+64,off+96,srcreg+srcinc+d2,d1); \
\
	xptr(srcreg+d2) = xmm0;	/* Save R1 */ \
	xptr(srcreg+d2+16) = xmm4;	/* Save I1 */ \
	xptr(srcreg+d2+32) = xmm1;	/* Save R2 */ \
	xptr(srcreg+d2+48) = xmm5;	/* Save I2 */ \
	xptr(srcreg+d2+d1+16) = xmm2;	/* Save I3 */ \
	xptr(srcreg+d2+d1+32) = xmm3;	/* Save R4 */ \
	xptr(srcreg+d2+d1+48) = xmm7;	/* Save I4 */ \
	srcreg += srcinc; \
\
\
/* 114.05 clocks */ \

#define x4c_fft4_cmn(r1,r2,r3,r4,r5,r6,r7,r8,mem4,mem8,dest3,off1,off2,off3,off4,pre1,pre2) { \
	vec2f64 temp; \
	temp = xptr(rdi+off1+16); \
	r4 = temp;	/* premul_real/premul_imag */ \
	r4 *= r1;			/* A1 = R1 * premul_real/premul_imag */ \
	r8 = temp;	/* premul_real/premul_imag */ \
	r8 *= r5;			/* B1 = I1 * premul_real/premul_imag */ \
	r4 -= r5;			/* A1 = A1 - I1 */ \
	r8 += r1;			/* B1 = B1 + R1 */ \
	xprefetchw(u8ptr(pre1)); \
\
	temp = xptr(rdi+off3+16); \
	r1 = temp;	/* premul_real/premul_imag */ \
	r1 *= r3;			/* A3 = R3 * premul_real/premul_imag */ \
	r5 = temp;	/* premul_real/premul_imag */ \
	r5 *= r7;			/* B3 = I3 * premul_real/premul_imag */ \
	r1 -= r7;			/* A3 = A3 - I3 */ \
	r5 += r3;			/* B3 = B3 + R3 */ \
	temp = xptr(rdi+off1); \
	r4 *= temp;		/* A1 = A1 * premul_imag (new R1) */ \
	r8 *= temp;		/* B1 = B1 * premul_imag (new I1) */ \
\
	temp = xptr(rdi+off2+16); \
	r3 = temp;	/* premul_real/premul_imag */ \
	r3 *= r2;		 	/* A2 = R2 * premul_real/premul_imag */ \
	r7 = temp;	/* premul_real/premul_imag */ \
	r7 *= r6;			/* B2 = I2 * premul_real/premul_imag */ \
	r3 -= r6;			/* A2 = A2 - I2 */ \
	r7 += r2;			/* B2 = B2 + R2 */ \
	temp = xptr(rdi+off3); \
	r1 *= temp;		/* A3 = A3 * premul_imag (new R3) */ \
	r5 *= temp;		/* B3 = B3 * premul_imag (new I3) */ \
\
	temp = xptr(rdi+off4+16); \
	r2 = temp;	/* premul_real/premul_imag */ \
	r2 *= mem4;	 	/* A4 = R4 * premul_real/premul_imag */ \
	r6 = temp;	/* premul_real/premul_imag */ \
	r6 *= mem8;		/* B4 = I4 * premul_real/premul_imag */ \
	r2 -= mem8;		/* A4 = A4 - I4 */ \
	r6 += mem4;		/* B4 = B4 + R4 */ \
	temp = xptr(rdi+off2); \
	r3 *= temp;		/* A2 = A2 * premul_imag (new R2) */ \
	r7 *= temp;		/* B2 = B2 * premul_imag (new I2) */ \
	temp = xptr(rdi+off4); \
	r2 *= temp;		/* A4 = A4 * premul_imag (new R4) */ \
	r6 *= temp;		/* B4 = B4 * premul_imag (new I4) */ \
	xprefetchw(u8ptr(pre1+pre2)); \
\
	temp = r4; \
	r4 -= r1;			/* R1 = R1 - R3 (new R3) */ \
	r1 += temp;			/* R3 = R1 + R3 (new R1) */ \
	temp = r8; \
	r8 -= r5;			/* I1 = I1 - I3 (new I3) */ \
	r5 += temp;			/* I3 = I1 + I3 (new I1) */ \
	temp = r3; \
	r3 -= r2;			/* R2 = R2 - R4 (new R4) */ \
	r2 += temp;			/* R4 = R2 + R4 (new R2) */ \
	temp = r7; \
	r7 -= r6;			/* I2 = I2 - I4 (new I4) */ \
	r6 += temp;			/* I4 = I2 + I4 (new I2) */ \
\
	temp = r4; \
	r4 -= r7;			/* R3 = R3 - I4 (final R3) */ \
	r7 += temp;			/* I4 = R3 + I4 (final R4) */ \
	dest3 = r4; \
	r4 = r8; \
	r8 -= r3;			/* I3 = I3 - R4 (final I4) */ \
	r3 += r4;			/* R4 = I3 + R4 (final I3) */ \
	r4 = r1; \
	r1 -= r2;			/* R1 = R1 - R2 (final R2) */ \
	r2 += r4;			/* R2 = R1 + R2 (final R1) */ \
	r4 = r5; \
	r5 -= r6;			/* I1 = I1 - I2 (final I2) */ \
	r6 += r4;			/* I2 = I1 + I2 (final I1) */ \
}


#define x4cl_four_complex_cpm_unfft(srcreg,srcinc,d1,d2) \
	xmm0 = xptr(srcreg);		/* R1 */ \
	xmm1 = xptr(srcreg+32);	/* R2 */ \
	xmm4 = xptr(srcreg+d2);	/* R5 */ \
	xmm5 = xptr(srcreg+d2+32);	/* R6 */ \
	x4c_cpm_unfft(xmm0,xmm1,xmm2,xmm3,xmm4,xmm5,xmm6,xmm7,xptr(srcreg+d1),xptr(srcreg+d1+32),xptr(srcreg+d1+d2),xptr(srcreg+d1+d2+32),xptr(srcreg),xptr(srcreg+32),xptr(srcreg+d2),xptr(srcreg+d2+32),srcreg+srcinc,d1); \
	xmm5 = xptr(srcreg+16);	/* R1 */ \
	xmm6 = xptr(srcreg+48);	/* R2 */ \
	xmm0 = xptr(srcreg+d2+16);	/* R5 */ \
	xmm1 = xptr(srcreg+d2+48);	/* R6 */ \
	xptr(srcreg+16) = xmm4;	/* Save R3 */ \
	xptr(srcreg+48) = xmm7;	/* Save I3 */ \
	xptr(srcreg+d2+16) = xmm2;	/* Save R4 */ \
	xptr(srcreg+d2+48) = xmm3;	/* Save I4 */ \
	x4c_cpm_unfft(xmm5,xmm6,xmm4,xmm7,xmm0,xmm1,xmm2,xmm3,xptr(srcreg+d1+16),xptr(srcreg+d1+48),xptr(srcreg+d1+d2+16),xptr(srcreg+d1+d2+48),xptr(srcreg+d1),xptr(srcreg+d1+32),xptr(srcreg+d1+d2),xptr(srcreg+d1+d2+32),srcreg+srcinc+d2,d1); \
	xptr(srcreg+d1+16) = xmm0;	/* Save R3 */ \
	xptr(srcreg+d1+48) = xmm3;	/* Save I3 */ \
	xptr(srcreg+d2+d1+16) = xmm4;	/* Save R4 */ \
	xptr(srcreg+d2+d1+48) = xmm7;	/* Save I4 */ \
	srcreg += srcinc; \
\
\

#define x4c_cpm_unfft(r1,r2,r3,r4,r5,r6,r7,r8,mem3,mem4,mem7,mem8,dest1,dest2,dest3,dest4,pre1,pre2) \
	x4c_unfft4_cmn(r1,r2,r3,r4,r5,r6,r7,r8,mem3,mem4,mem7,mem8,dest1,dest2,0,32,64,96,pre1,pre2); \
	dest1 = r6;		/* Save R1 */ \
	dest2 = r7;		/* Save I1 */ \
	dest3 = r1;		/* Save R2 */ \
	dest4 = r2;		/* Save I2 */ \
\
\

#define x4c_unfft4_cmn(r1,r2,r3,r4,r5,r6,r7,r8,mem3,mem4,mem7,mem8,dest1,dest2,off1,off2,off3,off4,pre1,pre2) { \
	vec2f64 temp; \
	r7 = r2; \
	r3 = mem3; \
	r4 = mem4; \
	r2 -= r4;			/* new I2 = I1 - I2 */ \
	r8 = mem8; \
	r8 += r6;			/* new I3 = I3 + I4 */ \
	r4 += r7;			/* new I1 = I1 + I2 */ \
	r7 = r1; \
	r1 -= r3;			/* new R2 = R1 - R2 */ \
	r3 += r7;			/* new R1 = R1 + R2 */ \
	r7 = mem7; \
	r7 -= r5;			/* new I4 = R4 - R3 */ \
	r5 += mem7;		/* new R3 = R4 + R3 */ \
	r6 -= mem8;		/* new R4 = I3 - I4 */ \
\
	xprefetchw(u8ptr(pre1)); \
\
	temp = r4; \
	r4 -= r8;			/* I1 = I1 - I3 (new I3) */ \
	r8 += temp;			/* I3 = I1 + I3 (new I1) */ \
	temp = r3; \
	r3 -= r5;			/* R1 = R1 - R3 (new R3) */ \
	r5 += temp;			/* R3 = R1 + R3 (new R1) */ \
	temp = r2; \
	r2 -= r7;			/* I2 = I2 - I4 (new I4) */ \
	r7 += temp;			/* I4 = I2 + I4 (new I2) */ \
	temp = r1; \
	r1 -= r6;			/* R2 = R2 - R4 (new R4) */ \
	r6 += temp;			/* R4 = R2 + R4 (new R2) */ \
\
	xprefetchw(u8ptr(pre1+pre2)); \
	dest1 = r5;		/* save intermediate R1 */ \
	dest2 = r8;		/* save intermediate I1 */ \
	temp = xptr(rdi+off3+16); \
	r5 = temp;	/* pre_real/pre_imag */ \
	r5 *= r3;			/* A3 = R3 * pre_real/pre_imag */ \
	r8 = temp;	/* pre_real/pre_imag */ \
	r8 *= r4;			/* B3 = I3 * pre_real/pre_imag */ \
	r5 += r4;			/* A3 = A3 + I3 */ \
	r8 -= r3;			/* B3 = B3 - R3 */ \
	temp = xptr(rdi+off4+16); \
	r3 = temp;	/* pre_real/pre_imag */ \
	r3 *= r1;			/* A4 = R4 * pre_real/pre_imag */ \
	r4 = temp;	/* pre_real/pre_imag */ \
	r4 *= r2;			/* B4 = I4 * pre_real/pre_imag */ \
	r3 += r2;			/* A4 = A4 + I4 */ \
	r4 -= r1;			/* B4 = B4 - R4 */ \
	temp = xptr(rdi+off2+16); \
	r1 = temp;	/* pre_real/pre_imag */ \
	r1 *= r6;			/* A2 = R2 * pre_real/pre_imag */ \
	r2 = temp;	/* pre_real/pre_imag */ \
	r2 *= r7;			/* B2 = I2 * pre_real/pre_imag */ \
	temp = xptr(rdi+off4); \
	r3 *= temp;		/* A4 = A4 * pre_imag (final R4) */ \
	r4 *= temp;		/* B4 = B4 * pre_imag (final I4) */ \
	r1 += r7;			/* A2 = A2 + I2 */ \
	r2 -= r6;			/* B2 = B2 - R2 */ \
	temp = xptr(rdi+off3); \
	r5 *= temp;		/* A3 = A3 * pre_imag (final R3) */ \
	r8 *= temp;		/* B3 = B3 * pre_imag (final I3) */ \
	temp = xptr(rdi+off1+16); \
	r6 = temp;	/* pre_real/pre_imag */ \
	r6 *= dest1;		/* A1 = R1 * pre_real/pre_imag */ \
	r7 = temp;	/* pre_real/pre_imag */ \
	r7 *= dest2;		/* B1 = I1 * pre_real/pre_imag */ \
	r6 += dest2;		/* A1 = A1 + I1 */ \
	r7 -= dest1;		/* B1 = B1 - R1 */ \
	temp = xptr(rdi+off2); \
	r1 *= temp;		/* A2 = A2 * pre_imag (final R2) */ \
	r2 *= temp;		/* B2 = B2 * pre_imag (final I2) */ \
	temp = xptr(rdi+off1); \
	r6 *= temp;		/* A1 = A1 * pre_imag (final R1) */ \
	r7 *= temp;		/* B1 = B1 * pre_imag (final I1) */ \
}


#define x3cl_three_complex_first_fft(srcreg,srcinc,d1) { \
	uintptr_t src_rbx = srcreg+rbx; \
	xmm0 = xptr(src_rbx);		/* R1 */ \
	xmm3 = xptr(src_rbx+32);		/* R4 */ \
	xmm2 = xptr(src_rbx+d1+16);	/* R3 */ \
	xmm5 = xptr(src_rbx+d1+48);	/* R6 */ \
	xmm1 = xptr(src_rbx+2*d1);	/* R2 */ \
	xmm4 = xptr(src_rbx+2*d1+32);	/* R5 */ \
	x3c_premult(xmm0, xmm1, xmm2, xmm3, xmm4, xmm5, xmm6, xmm7, 0); \
	x3c_fft(xmm0, xmm1, xmm2, xmm3, xmm4, xmm5, xmm6, xmm7); \
	xmm6 = xptr(src_rbx+d1);		/* R1 */ \
	xmm7 = xptr(src_rbx+d1+32);	/* R4 */ \
	xptr(srcreg) = xmm2;			/* Save R1 */ \
	xptr(srcreg+32) = xmm5;		/* Save R2 */ \
	xptr(srcreg+d1) = xmm0;		/* Save R3 */ \
	xptr(srcreg+d1+16) = xmm1;		/* Save R4 */ \
	xptr(srcreg+d1+32) = xmm4;		/* Save R5 */ \
	xptr(srcreg+d1+48) = xmm3;		/* Save R6 */ \
	xmm1 = xptr(src_rbx+16);		/* R2 */ \
	xmm4 = xptr(src_rbx+48);		/* R5 */ \
	xmm2 = xptr(src_rbx+2*d1+16);	/* R3 */ \
	xmm5 = xptr(src_rbx+2*d1+48);	/* R6 */ \
	x3c_premult(xmm6, xmm1, xmm2, xmm7, xmm4, xmm5, xmm0, xmm3, 96); \
	x3c_fft(xmm6, xmm1, xmm2, xmm7, xmm4, xmm5, xmm0, xmm3); \
	xptr(srcreg+16) = xmm2;		/* Save R1 */ \
	xptr(srcreg+48) = xmm5;		/* Save R2 */ \
	xptr(srcreg+2*d1) = xmm6;		/* Save R3 */ \
	xptr(srcreg+2*d1+16) = xmm1;		/* Save R4 */ \
	xptr(srcreg+2*d1+32) = xmm4;		/* Save R5 */ \
	xptr(srcreg+2*d1+48) = xmm7;		/* Save R6 */ \
	srcreg += srcinc; \
}

/* Same as above except we don't use the rbx offset as the FFT data has */ \
/* already been copied to the scratch area using the rbx offset. */ \

#define x3cl_three_complex_fft(srcreg,srcinc,d1) \
	xmm0 = xptr(srcreg);			/* R1 */ \
	xmm3 = xptr(srcreg+32);		/* R4 */ \
	xmm2 = xptr(srcreg+d1+16);		/* R3 */ \
	xmm5 = xptr(srcreg+d1+48);		/* R6 */ \
	xmm1 = xptr(srcreg+2*d1);		/* R2 */ \
	xmm4 = xptr(srcreg+2*d1+32);		/* R5 */ \
	x3c_premult(xmm0, xmm1, xmm2, xmm3, xmm4, xmm5, xmm6, xmm7, 0); \
	x3c_fft(xmm0, xmm1, xmm2, xmm3, xmm4, xmm5, xmm6, xmm7); \
	xmm6 = xptr(srcreg+d1);		/* R1 */ \
	xmm7 = xptr(srcreg+d1+32);		/* R4 */ \
	xptr(srcreg) = xmm2;			/* Save R1 */ \
	xptr(srcreg+32) = xmm5;		/* Save R2 */ \
	xptr(srcreg+d1) = xmm0;		/* Save R3 */ \
	xptr(srcreg+d1+16) = xmm1;		/* Save R4 */ \
	xptr(srcreg+d1+32) = xmm4;		/* Save R5 */ \
	xptr(srcreg+d1+48) = xmm3;		/* Save R6 */ \
	xmm1 = xptr(srcreg+16);		/* R2 */ \
	xmm4 = xptr(srcreg+48);		/* R5 */ \
	xmm2 = xptr(srcreg+2*d1+16);		/* R3 */ \
	xmm5 = xptr(srcreg+2*d1+48);		/* R6 */ \
	x3c_premult(xmm6, xmm1, xmm2, xmm7, xmm4, xmm5, xmm0, xmm3, 96); \
	x3c_fft(xmm6, xmm1, xmm2, xmm7, xmm4, xmm5, xmm0, xmm3); \
	xptr(srcreg+16) = xmm2;		/* Save R1 */ \
	xptr(srcreg+48) = xmm5;		/* Save R2 */ \
	xptr(srcreg+2*d1) = xmm6;		/* Save R3 */ \
	xptr(srcreg+2*d1+16) = xmm1;		/* Save R4 */ \
	xptr(srcreg+2*d1+32) = xmm4;		/* Save R5 */ \
	xptr(srcreg+2*d1+48) = xmm7;		/* Save R6 */ \
	srcreg += srcinc; \
\
\
/* This code is identical to the above except it is used in the single pass negacyclic FFTs and swizzles the inputs. */ \

#define s3cl_three_complex_first_fft(srcreg,srcinc,d1) { \
	uintptr_t src_rbx = srcreg+rbx; \
	low_load(xmm0, xptr(src_rbx), xptr(src_rbx+16)); /* R1 */ \
	low_load(xmm3, xptr(src_rbx+32), xptr(src_rbx+48)); /* R4 */ \
	high_load(xmm2, xptr(src_rbx+d1), xptr(src_rbx+d1+16)); /* R3 */ \
	high_load(xmm5, xptr(src_rbx+d1+32), xptr(src_rbx+d1+48)); /* R6 */ \
	low_load(xmm1, xptr(src_rbx+2*d1), xptr(src_rbx+2*d1+16)); /* R2 */ \
	low_load(xmm4, xptr(src_rbx+2*d1+32), xptr(src_rbx+2*d1+48)); /* R5 */ \
	x3c_premult(xmm0, xmm1, xmm2, xmm3, xmm4, xmm5, xmm6, xmm7, 0); \
	x3c_fft(xmm0, xmm1, xmm2, xmm3, xmm4, xmm5, xmm6, xmm7); \
	high_load(xmm6, xptr(src_rbx), xptr(src_rbx+16)); /* R2 */ \
	high_load(xmm7, xptr(src_rbx+32), xptr(src_rbx+48)); /* R5 */ \
	xptr(srcreg) = xmm2;			/* Save R1 */ \
	xptr(srcreg+32) = xmm5;		/* Save R2 */ \
	low_load(xmm2, xptr(src_rbx+d1), xptr(src_rbx+d1+16)); /* R1 */ \
	low_load(xmm5, xptr(src_rbx+d1+32), xptr(src_rbx+d1+48));	/* R4 */ \
	xptr(srcreg+d1) = xmm0;		/* Save R3 */ \
	xptr(srcreg+d1+16) = xmm1;		/* Save R4 */ \
	xptr(srcreg+d1+32) = xmm4;		/* Save R5 */ \
	xptr(srcreg+d1+48) = xmm3;		/* Save R6 */ \
	high_load(xmm0, xptr(src_rbx+2*d1), xptr(src_rbx+2*d1+16)); /* R3 */ \
	high_load(xmm1, xptr(src_rbx+2*d1+32), xptr(src_rbx+2*d1+48)); /* R6 */ \
	x3c_premult(xmm2, xmm6, xmm0, xmm5, xmm7, xmm1, xmm3, xmm4, 96); \
	x3c_fft(xmm2, xmm6, xmm0, xmm5, xmm7, xmm1, xmm3, xmm4); \
	xptr(srcreg+16) = xmm0;		/* Save R1 */ \
	xptr(srcreg+48) = xmm1;		/* Save R2 */ \
	xptr(srcreg+2*d1) = xmm2;		/* Save R3 */ \
	xptr(srcreg+2*d1+16) = xmm6;		/* Save R4 */ \
	xptr(srcreg+2*d1+32) = xmm7;		/* Save R5 */ \
	xptr(srcreg+2*d1+48) = xmm5;		/* Save R6 */ \
	srcreg += srcinc; \
}


#define x3c_premult(r1, r2, r3, r4, r5, r6, t1, t2, off) \
	t1 = r3;			/* Copy R3 */ \
	r3 *= xptr(rdi+off+80);	/* A3 = R3 * premul_real/premul_imag */ \
	r3 -= r6;			/* A3 = A3 - I3 */ \
	r6 *= xptr(rdi+off+80);	/* B3 = I3 * premul_real/premul_imag */ \
	r6 += t1;			/* B3 = B3 + R3 */ \
\
	t1 = r2;			/* Copy R2 */ \
	r2 *= xptr(rdi+off+48);	/* A2 = R2 * premul_real/premul_imag */ \
	r2 -= r5;			/* A2 = A2 - I2 */ \
	r5 *= xptr(rdi+off+48);	/* B2 = I2 * premul_real/premul_imag */ \
	r5 += t1;			/* B2 = B2 + R2 */ \
\
	r3 *= xptr(rdi+off+64);	/* A3 = A3 * premul_imag (new R3) */ \
	r6 *= xptr(rdi+off+64);	/* B3 = B3 * premul_imag (new I3) */ \
\
	t1 = r1;			/* Copy R1 */ \
	r1 *= xptr(rdi+off+16);	/* A1 = R1 * premul_real/premul_imag */ \
	r1 -= r4;			/* A1 = A1 - I1 */ \
	r4 *= xptr(rdi+off+16);	/* B1 = I1 * premul_real/premul_imag */ \
	r4 += t1;			/* B1 = B1 + R1 */ \
\
	r2 *= xptr(rdi+off+32);	/* A2 = A2 * premul_imag (new R2) */ \
	r5 *= xptr(rdi+off+32);	/* B2 = B2 * premul_imag (new I2) */ \
\
	r1 *= xptr(rdi+off);		/* A1 = A1 * premul_imag (new R1) */ \
	r4 *= xptr(rdi+off);		/* B1 = B1 * premul_imag (new I1) */ \
\
\
/* Do a 3-complex FFT.  The input values are R1+R4i, R2+R5i, R3+R6i */ \
/* A 3-complex FFT is: */ \
/* Res1:  (R1+R2+R3) + (I1+I2+I3)i */ \
/* Res2:  (R1-.5R2-.866I2-.5R3+.866I3) + (I1-.5I2+.866R2-.5I3-.866R3)i */ \
/* Res3:  (R1-.5R2+.866I2-.5R3-.866I3) + (I1-.5I2-.866R2-.5I3+.866R3)i */ \

#define x3c_fft(r1, r2, r3, r4, r5, r6, t1, t2) \
	t1 = r2; \
	r2 -= r3;			/* R2 - R3 */ \
	t2 = r5; \
	r5 -= r6;			/* I2 - I3 */ \
	r3 += t1;			/* R2 + R3 */ \
	r2 *= CONST2_P866;		/* 0.866 * (R2 - R3) */ \
	r6 += t2;			/* I2 + I3 */ \
	r5 *= CONST2_P866;		/* 0.866 * (I2 - I3) */ \
	t1 = r3 * 0.5;			/* 0.5 * (R2 + R3) */ \
	r3 += r1;			/* R1 + R2 + R3 (final R1) */ \
	r1 -= t1;			/* (R1-.5R2-.5R3) */ \
	t1 = r6 * 0.5;			/* 0.5 * (I2 + I3) */ \
	r6 += r4;			/* I1 + I2 + I3 (final I1) */ \
	r4 -= t1;			/* (I1-.5I2-.5I3) */ \
	t1 = r1; \
	r1 -= r5;			/* Final R2 */ \
	t2 = r4; \
	r4 -= r2;			/* Final I3 */ \
	r5 += t1;			/* Final R3 */ \
	r2 += t2;			/* Final I2 */ \
\
\

#define x3cl_three_complex_last_unfft(srcreg,srcinc,d1) \
	xmm0 = xptr(srcreg);		/* R1 */ \
	xmm1 = xptr(srcreg+32);	/* R2 */ \
	xmm2 = xptr(srcreg+d1);	/* R3 */ \
	xmm3 = xptr(srcreg+d1+32);	/* R4 */ \
	xmm4 = xptr(srcreg+2*d1);	/* R5 */ \
	xmm5 = xptr(srcreg+2*d1+32);	/* R6 */ \
	x3c_unfft(xmm0, xmm1, xmm2, xmm3, xmm4, xmm5, xmm6, xmm7); \
	x3c_postmult(xmm4, xmm5, xmm3, xmm1, xmm0, xmm2, xmm6, xmm7, 0); \
	xptr(srcreg) = xmm3;		/* Save R1 */ \
	xptr(srcreg+32) = xmm1;	/* Save R4 */ \
	xptr(srcreg+2*d1) = xmm0;	/* Save R2 */ \
	xptr(srcreg+2*d1+32) = xmm2;	/* Save R5 */ \
	xmm2 = xptr(srcreg+d1+16);	/* R3 */ \
	xmm3 = xptr(srcreg+d1+48);	/* R4 */ \
	xptr(srcreg+d1+16) = xmm6;	/* Save R3 */ \
	xptr(srcreg+d1+48) = xmm7;	/* Save R6 */ \
	xmm0 = xptr(srcreg+16);	/* R1 */ \
	xmm1 = xptr(srcreg+48);	/* R2 */ \
	xmm4 = xptr(srcreg+2*d1+16);	/* R5 */ \
	xmm5 = xptr(srcreg+2*d1+48);	/* R6 */ \
	x3c_unfft(xmm0, xmm1, xmm2, xmm3, xmm4, xmm5, xmm6, xmm7); \
	x3c_postmult(xmm4, xmm5, xmm3, xmm1, xmm0, xmm2, xmm6, xmm7, 96); \
	xptr(srcreg+16) = xmm0;	/* Save R2 */ \
	xptr(srcreg+48) = xmm2;	/* Save R5 */ \
	xptr(srcreg+d1) = xmm3;	/* Save R1 */ \
	xptr(srcreg+d1+32) = xmm1;	/* Save R4 */ \
	xptr(srcreg+2*d1+16) = xmm6;	/* Save R3 */ \
	xptr(srcreg+2*d1+48) = xmm7;	/* Save R6 */ \
	srcreg += srcinc; \
\
\
/* Do a 3-complex inverse FFT.  The input values are R1+R2i, R3+R4i, R5+R6i */ \
/* A 3-complex inverse FFT is: */ \
/* Res1:  (R1+R2+R3) + (I1+I2+I3)i */ \
/* Res2:  (R1-.5R2+.866I2-.5R3-.866I3) + (I1-.5I2-.866R2-.5I3+.866R3)i */ \
/* Res3:  (R1-.5R2-.866I2-.5R3+.866I3) + (I1-.5I2+.866R2-.5I3-.866R3)i */ \

#define x3c_unfft(r1, r2, r3, r4, r5, r6, t1, t2) \
	t1 = r3; \
	r3 -= r5;			/* R2 - R3 */ \
	r5 += t1;			/* R2 + R3 */ \
	r3 *= CONST2_P866;		/* 0.866 * (R2 - R3) */ \
	t2 = r4; \
	r4 -= r6;			/* I2 - I3 */ \
	r6 += t2;			/* I2 + I3 */ \
	r4 *= CONST2_P866;		/* 0.866 * (I2 - I3) */ \
	t1 = r5 * 0.5;			/* 0.5 * (R2 + R3) */ \
	r5 += r1;			/* R1 + R2 + R3 (final R1) */ \
	r1 -= t1;			/* (R1-.5R2-.5R3) */ \
	t1 = r6 * 0.5;			/* 0.5 * (I2 + I3) */ \
	r6 += r2;			/* I1 + I2 + I3 (final I1) */ \
	r2 -= t1;			/* (I1-.5I2-.5I3) */ \
	t1 = r1; \
	r1 -= r4;			/* Final R3 */ \
	r4 += t1;			/* Final R2 */ \
	t2 = r2; \
	r2 -= r3;			/* Final I2 */ \
	r3 += t2;			/* Final I3 */ \
\
\

#define x3c_postmult(r1, r2, r3, r4, r5, r6, t1, t2, off) \
	t1 = xptr(rdi+off+80);	/* pre_real/pre_imag */ \
	t1 *= r5;			/* A3 = R3 * pre_real/pre_imag */ \
	t2 = xptr(rdi+off+80);	/* pre_real/pre_imag */ \
	t2 *= r6;			/* B3 = I3 * pre_real/pre_imag */ \
	t1 += r6;			/* A3 = A3 + I3 */ \
	t2 -= r5;			/* B3 = B3 - R3 */ \
\
	r5 = xptr(rdi+off+48);	/* pre_real/pre_imag */ \
	r5 *= r3;			/* A2 = R2 * pre_real/pre_imag */ \
	r6 = xptr(rdi+off+48);	/* pre_real/pre_imag */ \
	r6 *= r4;			/* B2 = I2 * pre_real/pre_imag */ \
	r5 += r4;			/* A2 = A2 + I2 */ \
	r6 -= r3;			/* B2 = B2 - R2 */ \
	t1 *= xptr(rdi+off+64);	/* A3 = A3 * pre_imag (final R3) */ \
	t2 *= xptr(rdi+off+64);	/* B3 = B3 * pre_imag (final I3) */ \
\
	r3 = xptr(rdi+off+16);	/* pre_real/pre_imag */ \
	r3 *= r1;			/* A1 = R1 * pre_real/pre_imag */ \
	r4 = xptr(rdi+off+16);	/* pre_real/pre_imag */ \
	r4 *= r2;			/* B1 = I1 * pre_real/pre_imag */ \
	r3 += r2;			/* A1 = A1 + I1 */ \
	r4 -= r1;			/* B1 = B1 - R1 */ \
	r5 *= xptr(rdi+off+32);	/* A2 = A2 * pre_imag (final R2) */ \
	r6 *= xptr(rdi+off+32);	/* B2 = B2 * pre_imag (final I2) */ \
	r3 *= xptr(rdi+off);		/* A1 = A1 * pre_imag (final R1) */ \
	r4 *= xptr(rdi+off);		/* B1 = B1 * pre_imag (final I1) */ \
\

