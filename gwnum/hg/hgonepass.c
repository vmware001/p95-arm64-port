#include <stdlib.h>
#include <ctype.h>
#include <stdio.h>

#include "../gwnum.h"

#include "../gwcommon.h"
#include "../unravel.h"
#include "../gwtables.h"
#include "../xarch.h"
#include "../xbasics.h"
#include "../xmult.h"
#include "hg.h"

void hgonepass_xmiddle_1(struct gwasm_data *__restrict g) {
	//puts("xmiddle_1");
	uintptr_t rsi, rdi;
	vec2f64 xmm0,xmm1,xmm2,xmm3,xmm4,xmm5,xmm6,xmm7;
	rsi = (uintptr_t)g->DESTARG; /* Restore source pointer */
	rdi = (uintptr_t)g->sincos1; /* Load sin/cos pointer */
	s2cl_eight_reals_fft_2_final(rsi, 2*64, 64);
	rdi += XMM_SCD;		/* Next sine/cosine pointer */
	unsigned int loop = g->count1;
	do{
		s2cl_four_complex_fft_final(rsi, 2*64, 64);
		rdi += XMM_SCD;		/* Next sine/cosine pointer */
	}while(--loop);
}
void hgonepass_xmiddle_1p(struct gwasm_data *__restrict g) {
	//puts("xmiddle_1p");
	uintptr_t rsi, rdi;
	vec2f64 xmm0,xmm1,xmm2,xmm3,xmm4,xmm5,xmm6,xmm7;
	rsi = (uintptr_t)g->DESTARG; /* Restore source pointer */
	rdi = (uintptr_t)g->sincos1; /* Load sin/cos pointer */
	unsigned int loop = g->count1;
	do{
		s2cl_four_complex_fft_final(rsi, 2*64, 64);
		rdi += XMM_SCD;		/* Next sine/cosine pointer */
	}while(--loop);
}
void hgonepass_xmiddle_2(struct gwasm_data *__restrict g) {
	//puts("xmiddle_2");
	uintptr_t rsi, rdi;
	vec2f64 xmm0,xmm1,xmm2,xmm3,xmm4,xmm5,xmm6,xmm7;
	rsi = (uintptr_t)g->DESTARG; /* Restore source pointer */
	rdi = (uintptr_t)g->sincos1; /* Load sin/cos pointer */
	s2cl_eight_reals_with_square_2(rsi, 2*64, 64);
	/* 319/383/447/511 iterations */
	rdi += XMM_SCD;		/* Next sine/cosine pointer */
	unsigned int loop = g->count1;
	do{
		/*printf("loop %u", loop);*/
		s2cl_four_complex_with_square(rsi, 2*64, 64);
		rdi += XMM_SCD;		/* Next sine/cosine pointer */
	}while(--loop);
}
void hgonepass_xmiddle_2p(struct gwasm_data *__restrict g) {
	//puts("xmiddle_2p");
	uintptr_t rsi, rdi;
	vec2f64 xmm0,xmm1,xmm2,xmm3,xmm4,xmm5,xmm6,xmm7;
	rsi = (uintptr_t)g->DESTARG; /* Restore source pointer */
	rdi = (uintptr_t)g->sincos1; /* Load sin/cos pointer */
	unsigned int loop = g->count1;
	do{
		s2cl_four_complex_with_square(rsi, 2*64, 64);
		rdi += XMM_SCD;		/* Next sine/cosine pointer */
	}while(--loop);
}
void hgonepass_xmiddle_3(struct gwasm_data *__restrict g) {
	//puts("xmiddle_3");
	uintptr_t rsi, rdi, rbp;
	vec2f64 xmm0,xmm1,xmm2,xmm3,xmm4,xmm5,xmm6,xmm7;
	rsi = (uintptr_t)g->DESTARG; /* Restore source pointer */
	rdi = (uintptr_t)g->sincos1; /* Load sin/cos pointer */
	rbp = (uintptr_t)g->DIST_TO_MULSRCARG;
	s2cl_eight_reals_with_mult_2(rsi, 2*64, 64);
	/* 319/383/447/511 iterations */
	rdi += XMM_SCD;		/* Next sine/cosine pointer */
	unsigned int loop = g->count1;
	do{
		s2cl_four_complex_with_mult(rsi, 2*64, 64);
		rdi += XMM_SCD;		/* Next sine/cosine pointer */
	}while(--loop);
}
void hgonepass_xmiddle_3p(struct gwasm_data *__restrict g) {
	//puts("xmiddle_3p");
	uintptr_t rsi, rdi, rbp;
	vec2f64 xmm0,xmm1,xmm2,xmm3,xmm4,xmm5,xmm6,xmm7;
	rsi = (uintptr_t)g->DESTARG; /* Restore source pointer */
	rdi = (uintptr_t)g->sincos1; /* Load sin/cos pointer */
	rbp = (uintptr_t)g->DIST_TO_MULSRCARG;
	unsigned int loop = g->count1;
	do{
		s2cl_four_complex_with_mult(rsi, 2*64, 64);
		rdi += XMM_SCD;		/* Next sine/cosine pointer */
	}while(--loop);
}
void hgonepass_xmiddle_4(struct gwasm_data *__restrict g) {
	//\puts("xmiddle_4");
	uintptr_t rsi, rdi, rbp, rbx;
	vec2f64 xmm0,xmm1,xmm2,xmm3,xmm4,xmm5,xmm6,xmm7;
	rsi = (uintptr_t)g->DESTARG; /* Restore source pointer */
	rdi = (uintptr_t)g->sincos1; /* Load sin/cos pointer */
	rbp = (uintptr_t)g->DIST_TO_MULSRCARG;
	rbx = (uintptr_t)g->DIST_TO_FFTSRCARG;
	s2cl_eight_reals_with_mulf_2(rsi, 2*64, 64);
	rdi += XMM_SCD;		/* Next sine/cosine pointer */
	unsigned int loop = g->count1;
	do{
		s2cl_four_complex_with_mulf(rsi, 2*64, 64);
		rdi += XMM_SCD;		/* Next sine/cosine pointer */
	}while(--loop);
}
void hgonepass_xmiddle_4p(struct gwasm_data *__restrict g) {
	//puts("xmiddle_4p");
	uintptr_t rsi, rdi, rbp, rbx;
	vec2f64 xmm0,xmm1,xmm2,xmm3,xmm4,xmm5,xmm6,xmm7;
	rsi = (uintptr_t)g->DESTARG; /* Restore source pointer */
	rdi = (uintptr_t)g->sincos1; /* Load sin/cos pointer */
	rbp = (uintptr_t)g->DIST_TO_MULSRCARG;
	rbx = (uintptr_t)g->DIST_TO_FFTSRCARG;
	unsigned int loop = g->count1;
	do{
		s2cl_four_complex_with_mulf(rsi, 2*64, 64);
		rdi += XMM_SCD;		/* Next sine/cosine pointer */
	}while(--loop);
}
