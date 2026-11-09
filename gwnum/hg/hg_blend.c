#include <stdlib.h>
#include <ctype.h>

#include "../gwnum.h"
#include "../gwcommon.h"
#include "../unravel.h"
#include "../gwtables.h"
#include "../xarch.h"
#include "../xbasics.h"
#include "../xmult.h"
#include "hg.h"
#include "hgonepass.h"
#include <stdio.h>

//             func. macro.
// negacyclic:   _ac    p

#define INITVAR1 \
	uintptr_t rbx,rsi,rdi; \
	vec2f64 xmm0,xmm1,xmm2,xmm3,xmm4,xmm5,xmm6,xmm7; \
	rsi = (uintptr_t)g->DESTARG; \
	rbx = (uintptr_t)g->DIST_TO_FFTSRCARG;

#define N452(A) #A
#define N_TO_STRING(A) N452(A)
#define P1 static char Z1=0; if(Z1 == 0) \
	{Z1=1; puts("hg_blend.c " N_TO_STRING(__LINE__));}

#if 1
void xfft_hg_32_op_BLEND(struct gwasm_data *__restrict g) {INITVAR1; xfft32(0);}
void xfft_hg_48_op_BLEND(struct gwasm_data *__restrict g) {INITVAR1; xfft48(0);}
void xfft_hg_64_op_BLEND(struct gwasm_data *__restrict g) {INITVAR1; xfft64(0);}
void xfft_hg_80_op_BLEND(struct gwasm_data *__restrict g) {INITVAR1; xfft80(0);}
void xfft_hg_96_op_BLEND(struct gwasm_data *__restrict g) {INITVAR1; xfft96(0);}
void xfft_hg_112_op_BLEND(struct gwasm_data *__restrict g) {INITVAR1; xfft112(0);}
void xfft_hg_128_op_BLEND(struct gwasm_data *__restrict g) {INITVAR1; xfft128(0);}
void xfft_hg_160_op_BLEND(struct gwasm_data *__restrict g) {INITVAR1; xfft160(0);}
void xfft_hg_192_op_BLEND(struct gwasm_data *__restrict g) {INITVAR1; xfft192(0);}
void xfft_hg_224_op_BLEND(struct gwasm_data *__restrict g) {INITVAR1; xfft224(0);}
void xfft_hg_256_op_BLEND(struct gwasm_data *__restrict g) {INITVAR1; xfft256(0);}
void xfft_hg_320_op_BLEND(struct gwasm_data *__restrict g) {INITVAR1; xfft320(0);}
void xfft_hg_384_op_BLEND(struct gwasm_data *__restrict g) {INITVAR1; xfft384(0);}
void xfft_hg_448_op_BLEND(struct gwasm_data *__restrict g) {INITVAR1; xfft448(0);}
void xfft_hg_512_op_BLEND(struct gwasm_data *__restrict g) {INITVAR1; xfft512(0);}
void xfft_hg_640_op_BLEND(struct gwasm_data *__restrict g) {INITVAR1; xfft640(0);}
void xfft_hg_768_op_BLEND(struct gwasm_data *__restrict g) {INITVAR1; xfft768(0);}
void xfft_hg_896_op_BLEND(struct gwasm_data *__restrict g) {INITVAR1; xfft896(0);}
void xfft_hg_1024_op_BLEND(struct gwasm_data *__restrict g) {INITVAR1; xfft1024(0);}
void xfft_hg_1280_op_BLEND(struct gwasm_data *__restrict g) {INITVAR1; xfft1280(0);}
void xfft_hg_1536_op_BLEND(struct gwasm_data *__restrict g) {INITVAR1; xfft1536(0);}
void xfft_hg_1792_op_BLEND(struct gwasm_data *__restrict g) {INITVAR1; xfft1792(0);}
void xfft_hg_2048_op_BLEND(struct gwasm_data *__restrict g) {INITVAR1; xfft2048(0);}
void xfft_hg_2560_op_BLEND(struct gwasm_data *__restrict g) {INITVAR1; xfft2560(0);}
void xfft_hg_3072_op_BLEND(struct gwasm_data *__restrict g) {INITVAR1; xfft3072(0);}
void xfft_hg_3584_op_BLEND(struct gwasm_data *__restrict g) {INITVAR1; xfft3584(0);}
void xfft_hg_4096_op_BLEND(struct gwasm_data *__restrict g) {INITVAR1; xfft4096(0);}
void xfft_hg_5120_op_BLEND(struct gwasm_data *__restrict g) {INITVAR1; xfft5120(0);}
void xfft_hg_6144_op_BLEND(struct gwasm_data *__restrict g) {INITVAR1; xfft6144(0);}
//void xfft_hg_7168_op_BLEND(struct gwasm_data *__restrict g) {INITVAR1; xfft7168(0);}
void xfft_hg_8192_op_BLEND(struct gwasm_data *__restrict g) {INITVAR1; xfft8192(0);}
#endif
#if 1
void xfft_hg_32_op_ac_BLEND(struct gwasm_data *__restrict g) {INITVAR1; xfft32p(0);}
void xfft_hg_48_op_ac_BLEND(struct gwasm_data *__restrict g) {INITVAR1; xfft48p(0);}
void xfft_hg_64_op_ac_BLEND(struct gwasm_data *__restrict g) {INITVAR1; xfft64p(0);}
void xfft_hg_96_op_ac_BLEND(struct gwasm_data *__restrict g) {INITVAR1; xfft96p(0);}
void xfft_hg_128_op_ac_BLEND(struct gwasm_data *__restrict g) {INITVAR1; xfft128p(0);}
void xfft_hg_192_op_ac_BLEND(struct gwasm_data *__restrict g) {INITVAR1; xfft192p(0);}
void xfft_hg_256_op_ac_BLEND(struct gwasm_data *__restrict g) {INITVAR1; xfft256p(0);}
void xfft_hg_384_op_ac_BLEND(struct gwasm_data *__restrict g) {INITVAR1; xfft384p(0);}
void xfft_hg_512_op_ac_BLEND(struct gwasm_data *__restrict g) {INITVAR1; xfft512p(0);}
void xfft_hg_768_op_ac_BLEND(struct gwasm_data *__restrict g) {INITVAR1; xfft768p(0);}
void xfft_hg_1024_op_ac_BLEND(struct gwasm_data *__restrict g) {INITVAR1; xfft1024p(0);}
void xfft_hg_1536_op_ac_BLEND(struct gwasm_data *__restrict g) {INITVAR1; xfft1536p(0);}
void xfft_hg_2048_op_ac_BLEND(struct gwasm_data *__restrict g) {INITVAR1; xfft2048p(0);}
void xfft_hg_3072_op_ac_BLEND(struct gwasm_data *__restrict g) {INITVAR1; xfft3072p(0);}
void xfft_hg_4096_op_ac_BLEND(struct gwasm_data *__restrict g) {INITVAR1; xfft4096p(0);}
void xfft_hg_6144_op_ac_BLEND(struct gwasm_data *__restrict g) {INITVAR1; xfft6144p(0);}
//void xfft_hg_8192_op_ac_BLEND(struct gwasm_data *__restrict g) {INITVAR1; xfft8192p(0);}
#endif




// xfft_hg_32_BLEND

//
	//
//
