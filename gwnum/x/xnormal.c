#include <stdlib.h>
#include <ctype.h>

#include "../gwnum.h"
#include "../gwcommon.h"
#include "../unravel.h"
#include "../gwtables.h"
#include "../xarch.h"
#include "../xbasics.h"
#include "xnormal.h"

// xfunc572 used in 6 macros:
// xnorm012_1d_zpad
// xnorm012_2d_zpad_part1_cmn
// xnorm012_wpn_zpad_part1_cmn
// xnorm_smallmul_1d_zpad_cleanup
// xnorm_smallmul_2d_fft_zpad
// xnorm_smallmul_wpn_fft_zpad

#if 0
void xfunc572(struct gwasm_data *__restrict g) {
	double f0, f1, f2, f4, f5; \
;	/* Divide the zpad data by k.  Store the integer part in XMM_TMP */ \
;	/* and the remainder in ZPAD0.  Later we will wrap the integer part */ \
;	/* down to the bottom of the FFT data area (and multiply by -c). */ \
;	/* And we will store the remainder in the upper half of the FFT */ \
;	/* data area. */ \
\
;	/* Note there are three cases to handle.  K is smaller than a big word. */ \
;	/* K is between one and 2 big words in size.  And K is more than */ \
;	/* 2 big words in size. */ \
\
	if(g->ZPAD_TYPE >= 2) {		/* Are we dealing with case 1,2,or 3 */ \
	if(g->ZPAD_TYPE > 2) {
		/* jl smallk; */			/* One word case */ \
		/* je mediumk; */			/* Two word case */ \
	\
		/* This case does the divide by k where k is three words */ \
		f0 = g->ZPAD0_6[6];/* Load zpad word (high bits) */ \
		f1 = g->ZPAD0_6[5];/* Load zpad word (middle bits) */ \
		f2 = g->ZPAD0_6[4];/* Load zpad word (low bits) */ \
		f4 = g->ZPAD_INVERSE_K6;	/* Load shifted 1/k */ \
		f4 *= f0;		/* Mul ZPAD by shifted 1/k */ \
		f4 = roundf64(f4);	/* Round to integer	 */ \
		f5 = g->ZPAD_K6_HI;	/* Load high bits of k */ \
		f5 *= f4; \
		f0 -= f5;		/* Calculate high bits of remainder */ \
		f5 = g->ZPAD_K6_MID;	/* Load middle bits of k */ \
		f5 *= f4; \
		f1 -= f5;		/* Calculate middle bits of remainder */ \
		f5 = g->ZPAD_K6_LO;	/* Load low bits of k */ \
		f5 *= f4; \
		f2 -= f5;		/* Calculate low bits of remainder */ \
		g->u.xmm.XMM_TMP5[0] = f4;	/* Save word of zpad / k */ \
	\
		f0 *= g->ZPAD_SHIFT6;	/* Shift previous zpad word */ \
		f1 += f0;		/* Add to create new high zpad bits */ \
		f0 = g->ZPAD0_6[3];/* Load zpad word (new low bits) */ \
		f5 = g->ZPAD_SHIFT5;	/* Combine high and medium bits */ \
		f5 *= f1; \
		f5 += f2; \
		f4 = g->ZPAD_INVERSE_K5;	/* Load shifted 1/k */ \
		f4 *= f5;		/* Mul ZPAD by shifted 1/k */ \
		f4 = roundf64(f4);/* Round to integer */ \
		f5 = g->ZPAD_K5_HI;	/* Load high bits of k */ \
		f5 *= f4; \
		f1 -= f5;		/* Calculate high bits of remainder */ \
		f5 = g->ZPAD_K5_MID;	/* Load middle bits of k */ \
		f5 *= f4; \
		f2 -= f5;		/* Calculate middle bits of remainder */ \
		f5 = g->ZPAD_K5_LO;	/* Load low bits of k */ \
		f5 *= f4; \
		f0 -= f5;		/* Calculate low bits of remainder */ \
		g->u.xmm.XMM_TMP4[0] = f4;	/* Save word of zpad / k */ \
	\
		f1 *= g->ZPAD_SHIFT5;	/* Shift previous zpad word */ \
		f2 += f1;		/* Add to create new high zpad bits */ \
		f1 = g->ZPAD0_6[2];/* Load zpad word (new low bits) */ \
		f5 = g->ZPAD_SHIFT4;	/* Combine high and medium bits */ \
		f5 *= f2; \
		f5 += f0; \
		f4 = g->ZPAD_INVERSE_K4;	/* Load shifted 1/k */ \
		f4 *= f5;		/* Mul ZPAD by shifted 1/k */ \
		f4 = roundf64(f4);/* Round to integer */ \
		f5 = g->ZPAD_K4_HI;	/* Load high bits of k */ \
		f5 *= f4; \
		f2 -= f5;		/* Calculate high bits of remainder */ \
		f5 = g->ZPAD_K4_MID;	/* Load middle bits of k */ \
		f5 *= f4; \
		f0 -= f5;		/* Calculate middle bits of remainder */ \
		f5 = g->ZPAD_K4_LO;	/* Load low bits of k */ \
		f5 *= f4; \
		f1 -= f5;		/* Calculate low bits of remainder */ \
		g->u.xmm.XMM_TMP3[0] = f4;	/* Save word of zpad / k */ \
	\
		f2 *= g->ZPAD_SHIFT4;	/* Shift previous zpad word */ \
		f0 += f2;		/* Add to create new high zpad bits */ \
		f2 = g->ZPAD0_6[1];/* Load zpad word (new low bits) */ \
		f5 = g->ZPAD_SHIFT3;	/* Combine high and medium bits */ \
		f5 *= f0; \
		f5 += f1; \
		f4 = g->ZPAD_INVERSE_K3;	/* Load shifted 1/k */ \
		f4 *= f5;		/* Mul ZPAD by shifted 1/k */ \
		f4 = roundf64(f4);/* Round to integer */ \
		f5 = g->ZPAD_K3_HI;	/* Load high bits of k */ \
		f5 *= f4; \
		f0 -= f5;		/* Calculate high bits of remainder */ \
		f5 = g->ZPAD_K3_MID;	/* Load middle bits of k */ \
		f5 *= f4; \
		f1 -= f5;		/* Calculate middle bits of remainder */ \
		f5 = g->ZPAD_K3_LO;	/* Load low bits of k */ \
		f5 *= f4; \
		f2 -= f5;		/* Calculate low bits of remainder */ \
		g->u.xmm.XMM_TMP2[0] = f4;	/* Save word of zpad / k */ \
	\
		f0 *= g->ZPAD_SHIFT3;	/* Shift previous zpad word */ \
		f1 += f0;		/* Add to create new high zpad bits */ \
		f0 = g->ZPAD0_6[0];/* Load zpad word (new low bits) */ \
		f5 = g->ZPAD_SHIFT2;	/* Combine high and medium bits */ \
		f5 *= f1; \
		f5 += f2; \
		f4 = g->ZPAD_INVERSE_K2;	/* Load shifted 1/k */ \
		f4 *= f5;		/* Mul ZPAD by shifted 1/k */ \
		f4 = roundf64(f4);/* Round to integer */ \
		f5 = g->ZPAD_K2_HI;	/* Load high bits of k */ \
		f5 *= f4; \
		f1 -= f5;		/* Calculate high bits of remainder */ \
		f5 = g->ZPAD_K2_MID;	/* Load middle bits of k */ \
		f5 *= f4; \
		f2 -= f5;		/* Calculate middle bits of remainder */ \
		f5 = g->ZPAD_K2_LO;	/* Load low bits of k */ \
		f5 *= f4; \
		f0 -= f5;		/* Calculate low bits of remainder */ \
		g->u.xmm.XMM_TMP1[0] = f4;	/* Save word of zpad / k */ \
	\
		f1 *= g->ZPAD_SHIFT2;	/* Shift previous zpad word */ \
		f2 += f1;		/* Add to create new high zpad bits */ \
		f2 *= g->ZPAD_SHIFT1;	/* Shift previous zpad word */ \
		f0 += f2;		/* Add to create new high zpad bits */ \
		g->ZPAD0_6[0] = f0;		/* Save remainder of zpad / k */ \
	\
		f1 = 0;		/* Zero words that other cases set */ \
		g->u.xmm.XMM_TMP6[0] = f1; \
	}else{ \
		/* This case does the divide by k where k is two words */ \
		f0 = g->ZPAD0_6[6];/* Load zpad word (high bits) */ \
		f1 = g->ZPAD0_6[5];/* Load zpad word (low bits) */ \
		f4 = g->ZPAD_INVERSE_K6;	/* Load shifted 1/k */ \
		f4 *= f0;		/* Mul ZPAD by shifted 1/k */ \
		f4 = roundf64(f4);	/* Round to integer	 */ \
		f5 = g->ZPAD_K6_HI;	/* Load high bits of k */ \
		f5 *= f4; \
		f0 -= f5;		/* Calculate high bits of remainder */ \
		f5 = g->ZPAD_K6_LO;	/* Load low bits of k */ \
		f5 *= f4; \
		f1 -= f5;		/* Calculate low bits of remainder */ \
		g->u.xmm.XMM_TMP6[0] = f4;	/* Save word of zpad / k */ \
	\
		f0 *= g->ZPAD_SHIFT6;	/* Shift previous zpad word */ \
		f0 += f1;		/* Add to create new high zpad bits */ \
		f1 = g->ZPAD0_6[4];/* Load zpad word (new low bits) */ \
		f4 = g->ZPAD_INVERSE_K5;	/* Load shifted 1/k */ \
		f4 *= f0;		/* Mul ZPAD by shifted 1/k */ \
		f4 = roundf64(f4);/* Round to integer */ \
		f5 = g->ZPAD_K5_HI;	/* Load high bits of k */ \
		f5 *= f4; \
		f0 -= f5;		/* Calculate high bits of remainder */ \
		f5 = g->ZPAD_K5_LO;	/* Load low bits of k */ \
		f5 *= f4; \
		f1 -= f5;		/* Calculate low bits of remainder */ \
		g->u.xmm.XMM_TMP5[0] = f4;	/* Save word of zpad / k */ \
	\
		f0 *= g->ZPAD_SHIFT5;	/* Shift previous zpad word */ \
		f0 += f1;		/* Add to create new high zpad bits */ \
		f1 = g->ZPAD0_6[3];/* Load zpad word (new low bits) */ \
		f4 = g->ZPAD_INVERSE_K4;	/* Load shifted 1/k */ \
		f4 *= f0;		/* Mul ZPAD by shifted 1/k */ \
		f4 = roundf64(f4);/* Round to integer */ \
		f5 = g->ZPAD_K4_HI;	/* Load high bits of k */ \
		f5 *= f4; \
		f0 -= f5;		/* Calculate high bits of remainder */ \
		f5 = g->ZPAD_K4_LO;	/* Load low bits of k */ \
		f5 *= f4; \
		f1 -= f5;		/* Calculate low bits of remainder */ \
		g->u.xmm.XMM_TMP4[0] = f4;	/* Save word of zpad / k */ \
	\
		f0 *= g->ZPAD_SHIFT4;	/* Shift previous zpad word */ \
		f0 += f1;		/* Add to create new high zpad bits */ \
		f1 = g->ZPAD0_6[2];/* Load zpad word (new low bits) */ \
		f4 = g->ZPAD_INVERSE_K3;	/* Load shifted 1/k */ \
		f4 *= f0;		/* Mul ZPAD by shifted 1/k */ \
		f4 = roundf64(f4);/* Round to integer */ \
		f5 = g->ZPAD_K3_HI;	/* Load high bits of k */ \
		f5 *= f4; \
		f0 -= f5;		/* Calculate high bits of remainder */ \
		f5 = g->ZPAD_K3_LO;	/* Load low bits of k */ \
		f5 *= f4; \
		f1 -= f5;		/* Calculate low bits of remainder */ \
		g->u.xmm.XMM_TMP3[0] = f4;	/* Save word of zpad / k */ \
	\
		f0 *= g->ZPAD_SHIFT3;	/* Shift previous zpad word */ \
		f0 += f1;		/* Add to create new high zpad bits */ \
		f1 = g->ZPAD0_6[1];/* Load zpad word (new low bits) */ \
		f4 = g->ZPAD_INVERSE_K2;	/* Load shifted 1/k */ \
		f4 *= f0;		/* Mul ZPAD by shifted 1/k */ \
		f4 = roundf64(f4);/* Round to integer */ \
		f5 = g->ZPAD_K2_HI;	/* Load high bits of k */ \
		f5 *= f4; \
		f0 -= f5;		/* Calculate high bits of remainder */ \
		f5 = g->ZPAD_K2_LO;	/* Load low bits of k */ \
		f5 *= f4; \
		f1 -= f5;		/* Calculate low bits of remainder */ \
		g->u.xmm.XMM_TMP2[0] = f4;	/* Save word of zpad / k */ \
	\
		f0 *= g->ZPAD_SHIFT2;	/* Shift previous zpad word */ \
		f0 += f1;		/* Add to create new high zpad bits */ \
		f1 = g->ZPAD0_6[0];/* Load zpad word (new low bits) */ \
		f4 = g->ZPAD_INVERSE_K1;	/* Load shifted 1/k */ \
		f4 *= f0;		/* Mul ZPAD by shifted 1/k */ \
		f4 = roundf64(f4);/* Round to integer */ \
		f5 = g->ZPAD_K1_HI;	/* Load high bits of k */ \
		f5 *= f4; \
		f0 -= f5;		/* Calculate high bits of remainder */ \
		f5 = g->ZPAD_K1_LO;	/* Load low bits of k */ \
		f5 *= f4; \
		f1 -= f5;		/* Calculate low bits of remainder */ \
		g->u.xmm.XMM_TMP1[0] = f4;	/* Save word of zpad / k */ \
	\
		f0 *= g->ZPAD_SHIFT1;	/* Shift previous zpad word */ \
		f0 += f1;		/* Add to create new high zpad bits */ \
		g->ZPAD0_6[0] = f0;		/* Save remainder of zpad / k */ \
	}
	}else{ \
		/* This case does the divide by k where k is one word */ \
		/* Assume ZPAD5 and ZPAD6 are zero. */ \
		f0 = g->ZPAD0_6[4];/* Load zpad data */ \
		f4 = g->ZPAD_INVERSE_K1;	/* Load by 1/k */ \
		f4 *= f0;		/* Mul ZPAD data by 1/k */ \
		f4 = roundf64(f4);	/* Round to integer	 */ \
		g->u.xmm.XMM_TMP5[0] = f4;	/* Save integer part */ \
		f4 *= g->ZPAD_K1_LO;	/* Compute remainder */ \
		f0 -= f4; \
	\
		f0 *= g->ZPAD_SHIFT4;	/* Shift previous zpad word */ \
		f0 += g->ZPAD0_6[3];/* Add in zpad data */ \
		f4 = g->ZPAD_INVERSE_K1;	/* Load by 1/k */ \
		f4 *= f0;		/* Mul ZPAD data by 1/k */ \
		f4 = roundf64(f4);	/* Round to integer	 */ \
		g->u.xmm.XMM_TMP4[0] = f4;	/* Save integer part */ \
		f4 *= g->ZPAD_K1_LO;	/* Compute remainder */ \
		f0 -= f4; \
	\
		f0 *= g->ZPAD_SHIFT3;	/* Shift previous zpad word */ \
		f0 += g->ZPAD0_6[2];/* Add in zpad data */ \
		f4 = g->ZPAD_INVERSE_K1;	/* Load by 1/k */ \
		f4 *= f0;		/* Mul ZPAD data by 1/k */ \
		f4 = roundf64(f4);	/* Round to integer	 */ \
		g->u.xmm.XMM_TMP3[0] = f4;	/* Save integer part */ \
		f4 *= g->ZPAD_K1_LO;	/* Compute remainder */ \
		f0 -= f4; \
	\
		f0 *= g->ZPAD_SHIFT2;	/* Shift previous zpad word */ \
		f0 += g->ZPAD0_6[1];/* Add in zpad data */ \
		f4 = g->ZPAD_INVERSE_K1;	/* Load by 1/k */ \
		f4 *= f0;		/* Mul ZPAD data by 1/k */ \
		f4 = roundf64(f4);	/* Round to integer	 */ \
		g->u.xmm.XMM_TMP2[0] = f4;	/* Save integer part */ \
		f4 *= g->ZPAD_K1_LO;	/* Compute remainder */ \
		f0 -= f4; \
	\
		f0 *= g->ZPAD_SHIFT1;	/* Shift previous zpad word */ \
		f0 += g->ZPAD0_6[0];/* Add in zpad data */ \
		f4 = g->ZPAD_INVERSE_K1;	/* Load by 1/k */ \
		f4 *= f0;		/* Mul ZPAD data by 1/k */ \
		f4 = roundf64(f4);	/* Round to integer	 */ \

		g->u.xmm.XMM_TMP1[0] = f4;	/* Save integer part */ \
		f4 *= g->ZPAD_K1_LO;	/* Compute remainder */ \
		f0 -= f4; \
		g->ZPAD0_6[0] = f0;		/* Save remainder */ \
	\
		f1 = 0;		/* Zero words that other cases set */ \
		g->u.xmm.XMM_TMP6[0] = f1; \
	};
}

#else
void xfunc572(struct gwasm_data *__restrict g) {
	double f0,f1,f2,f4,f5;
	vec2f64 TMP12, TMP34, TMP56;
	char ZPAD_TYPE = g->ZPAD_TYPE;
	if(g->ZPAD_TYPE >= 2) {
	if(g->ZPAD_TYPE > 2) {
		f0 = g->ZPAD0_6[6];
		f1 = g->ZPAD0_6[5];
		f2 = g->ZPAD0_6[4];

		f4 = f0 * g->ZPAD_INVERSE_K6;
		f4 = roundf64(f4);
		f5 = f4 * g->ZPAD_K6_HI; f0 -= f5;
		f5 = f4 * g->ZPAD_K6_MID; f1 -= f5;
		f5 = f4 * g->ZPAD_K6_LO; f2 -= f5;
		TMP56[0] = f4;
		TMP56[1] = 0.0;
		g->u.xmm.XMM_TMP3 = TMP56;

		f0 *= g->ZPAD_SHIFT6;
		f1 += f0;
		f0 = g->ZPAD0_6[3];
		f5 = f1 * g->ZPAD_SHIFT5;
		f5 += f2;
		f4 = f5 * g->ZPAD_INVERSE_K5;
		f4 = roundf64(f4);
		f5 = f4 * g->ZPAD_K5_HI; f1 -= f5;
		f5 = f4 * g->ZPAD_K5_MID; f2 -= f5;
		f5 = f4 * g->ZPAD_K5_LO; f0 -= f5;
		TMP34[1] = f4; /* XMM_TMP4 */

		f1 *= g->ZPAD_SHIFT5;
		f2 += f1;
		f1 = g->ZPAD0_6[2];
		f5 = f2 * g->ZPAD_SHIFT4;
		f5 += f0;
		f4 = f5 * g->ZPAD_INVERSE_K4;
		f4 = roundf64(f4);
		f5 = f4 * g->ZPAD_K4_HI; f2 -= f5;
		f5 = f4 * g->ZPAD_K4_MID; f0 -= f5;
		f5 = f4 * g->ZPAD_K4_LO; f1 -= f5;
		TMP34[0] = f4;
		g->u.xmm.XMM_TMP2 = TMP34;
		;
		f2 *= g->ZPAD_SHIFT4;
		f0 += f2;
		f2 = g->ZPAD0_6[1];
		f5 = f0 * g->ZPAD_SHIFT3;
		f5 += f1;
		f4 = f5 * g->ZPAD_INVERSE_K3;
		f4 = roundf64(f4);
		f5 = f4 * g->ZPAD_K3_HI; f0 -= f5;
		f5 = f4 * g->ZPAD_K3_MID; f1 -= f5;
		f5 = f4 * g->ZPAD_K3_LO; f2 -= f5;
		TMP12[1] = f4;

		f0 *= g->ZPAD_SHIFT3;
		f1 += f0;
		f0 = g->ZPAD0_6[0];
		f5 = f1 * g->ZPAD_SHIFT2;
		f5 += f2;
		f4 = f5 * g->ZPAD_INVERSE_K2;
		f4 = roundf64(f4);
		f5 = f4 * g->ZPAD_K2_HI; f1 -= f5;
		f5 = f4 * g->ZPAD_K2_MID; f2 -= f5;
		f5 = f4 * g->ZPAD_K2_LO; f0 -= f5;
		TMP12[0] = f4;
		g->u.xmm.XMM_TMP1 = TMP12;

		f1 *= g->ZPAD_SHIFT2;
		f2 += f1;
		f2 *= g->ZPAD_SHIFT1;
		f0 += f2;
		g->ZPAD0_6[0] = f0;
		;
	}else{
		f0 = g->ZPAD0_6[6];
		f1 = g->ZPAD0_6[5];
		f4 = f0 * g->ZPAD_INVERSE_K6;
		f4 = roundf64(f4);
		f5 = f4 * g->ZPAD_K6_HI; f0 -= f5;
		f5 = f4 * g->ZPAD_K6_LO; f1 -= f5;
		TMP56[1] = f4;

		f0 *= g->ZPAD_SHIFT6;
		f0 += f1;
		f1 = g->ZPAD0_6[4];
		f4 = f0 * g->ZPAD_INVERSE_K5;
		f4 = roundf64(f4);
		f5 = f4 * g->ZPAD_K5_HI; f0 -= f5;
		f5 = f4 * g->ZPAD_K5_LO; f1 -= f5;
		TMP56[0] = f4;
		g->u.xmm.XMM_TMP3 = TMP56;

		f0 *= g->ZPAD_SHIFT5;
		f0 += f1;
		f1 = g->ZPAD0_6[3];
		f4 = f0 * g->ZPAD_INVERSE_K4;
		f4 = roundf64(f4);
		f5 = f4 * g->ZPAD_K4_HI; f0 -= f5;
		f5 = f4 * g->ZPAD_K4_LO; f1 -= f5;
		TMP34[1] = f4;

		f0 *= g->ZPAD_SHIFT4;
		f0 += f1;
		f1 = g->ZPAD0_6[2];
		f4 = f0 * g->ZPAD_INVERSE_K3;
		f4 = roundf64(f4);
		f5 = f4 * g->ZPAD_K3_HI; f0 -= f5;
		f5 = f4 * g->ZPAD_K3_LO; f1 -= f5;
		TMP34[0] = f4;
		g->u.xmm.XMM_TMP2 = TMP34;

		f0 *= g->ZPAD_SHIFT3;
		f0 += f1;
		f1 = g->ZPAD0_6[1];
		f4 = f0 * g->ZPAD_INVERSE_K2;
		f4 = roundf64(f4);
		f5 = f4 * g->ZPAD_K2_HI; f0 -= f5;
		f5 = f4 * g->ZPAD_K2_LO; f1 -= f5;
		TMP12[1] = f4;

		f0 *= g->ZPAD_SHIFT2;
		f0 += f1;
		f1 = g->ZPAD0_6[0];
		f4 = f0 * g->ZPAD_INVERSE_K1;
		f4 = roundf64(f4);
		f5 = f4 * g->ZPAD_K1_HI; f0 -= f5;
		f5 = f4 * g->ZPAD_K1_LO; f1 -= f5;
		TMP12[0] = f4;
		g->u.xmm.XMM_TMP1 = TMP12;

		f0 *= g->ZPAD_SHIFT1;
		f0 += f1;
		g->ZPAD0_6[0] = f0;
	}
	}else{
		double ZPAD_INVERSE_K1 = g->ZPAD_INVERSE_K1;
		double ZPAD_K1_LO = g->ZPAD_K1_LO;
		f0 = g->ZPAD0_6[4];
		f4 = f0 * ZPAD_INVERSE_K1;
		f4 = roundf64(f4);
		TMP56[0] = f4;
		TMP56[1] = 0.0;
		f4 *= ZPAD_K1_LO;
		f0 -= f4;
		g->u.xmm.XMM_TMP3 = TMP56;

		f0 *= g->ZPAD_SHIFT4;
		f0 += g->ZPAD0_6[3];
		f4 = f0 * ZPAD_INVERSE_K1;
		f4 = roundf64(f4);
		TMP34[1] = f4;
		f4 *= ZPAD_K1_LO;
		f0 -= f4;

		f0 *= g->ZPAD_SHIFT3;
		f0 += g->ZPAD0_6[2];
		f4 = f0 * ZPAD_INVERSE_K1;
		f4 = roundf64(f4);
		TMP34[0] = f4;
		f4 *= ZPAD_K1_LO;
		f0 -= f4;
		g->u.xmm.XMM_TMP2 = TMP34;

		f0 *= g->ZPAD_SHIFT2;
		f0 += g->ZPAD0_6[1];
		f4 = f0 * ZPAD_INVERSE_K1;
		f4 = roundf64(f4);
		TMP12[1] = f4;
		f4 *= ZPAD_K1_LO;
		f0 -= f4;

		f0 *= g->ZPAD_SHIFT1;
		f0 += g->ZPAD0_6[0];
		f4 = f0 * ZPAD_INVERSE_K1;
		f4 = roundf64(f4);
		TMP12[0] = f4;
		f4 *= ZPAD_K1_LO;
		f0 -= f4;
		g->u.xmm.XMM_TMP1 = TMP12;
		g->ZPAD0_6[0] = f0;		/* Save remainder */
	};
}
#endif
