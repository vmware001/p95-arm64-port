#pragma once
#include <stdio.h>

// double define to allow it's inputs to expand correctly.
#define CONCAT2(A, B) A##B
#define CONCAT(A, B) CONCAT2(A, B)

#ifndef ARCH1
#define ARCH1 CORE
#endif


#define A_OR_Bnonoexecnoexec(C) C
#define A_OR_Bnonoexecexec(C) C
#define A_OR_Bnoexecnoexec(C)
#define A_OR_Bnoexecexec(C) C
#define A_OR_Bexecnoexec(C) C
#define A_OR_Bexecexec(C) C
#define A_OR_B(A, B, C) A_OR_B##A##B(C)

#define inorm(lab, ttp, echk, const1, base2, sse4) \
void CONCAT(lab, ARCH1)(struct gwasm_data *__restrict g, uintptr_t rsi) { \
	uintptr_t rax, rbx, rdi, rbp; \
A_OR_B##base2##ttp(uintptr_t rcx); \
A_OR_B##base2##ttp(uintptr_t r8); \
A_OR_B##base2##ttp(uintptr_t r9); \
ttp(uintptr_t rdx); \
uintptr_t loopcount1, loopcount2; \
	unsigned int loopcount3; \
	vec2f64 xmm0, xmm1, xmm2, xmm3, xmm4, xmm5, xmm7; \
	vec2f64 xmm8, xmm9, xmm10, xmm11, xmm12, xmm13; \
ttp(vec2f64 xmm14); \
ttp(vec2f64 xmm15); \
	uintptr_t saved_rsi = rsi; \
echk(vec2f64 xmm6 = g->u.xmm.XMM_MAXERR); \
	rbx = (uintptr_t) g->norm_ptr2;  \
ttp(rax = (uint32_t)g->cache_line_multiplier); \
	rdi = ((uintptr_t)g->u.xmm.XMM_COL_MULTS) + 128; /* LEA */ \
ttp(setlp:) \
	xnorm_2d_setup(ttp); \
ttp(rdi += 512); \
ttp(rbx += 32); \
ttp(rax -= 1); \
ttp(if(rax) goto setlp); \
	\
ttp(rdx = (uintptr_t)g->norm_grp_mults); \
	rbp = (uintptr_t)g->carries; \
	rdi = (uintptr_t)g->norm_ptr1; \
	rax = (uint32_t)g->addcount1; \
	loopcount2 = rax; \
	loopcount3 = 0; \
ttp(rax = u8ptr(rdi+0)); \
ttp(rcx = u8ptr(rdi+1)); \
ttp(r8 = u8ptr(rdi+2)); /*64b*/ \
ttp(r9 = u8ptr(rdi+3)); /*64b*/ \
no##ttp(base2(rax = 0)); \
no##ttp(base2(rcx = 0)); \
no##ttp(base2(r8 = 0)); /*64b*/ \
no##ttp(base2(r9 = 0)); /*64b*/ \
	do{ \
		rbx = (uint32_t)g->cache_line_multiplier; \
		loopcount1 = rbx; \
		rbx = ((uintptr_t)g->u.xmm.XMM_COL_MULTS); /* LEA */ \
		/* L2prefetch128 [rdx+128] */ \
		vec2f64 carry0 = xptr(rbp+0*16); \
		vec2f64 carry1 = xptr(rbp+1*16); \
		vec2f64 carry2 = xptr(rbp+2*16); \
		vec2f64 carry3 = xptr(rbp+3*16); \
		do{ /*ilp1:*/ \
			/* xprefetchw [rsi+64] */ \
			xnorm_2d(ttp, echk, const1, base2, sse4, carry0, carry1, carry2, carry3); \
			rsi += 64; \
		ttp(rbx += 512); \
		ttp(rdi += 4); \
			loopcount1--; \
		}while(loopcount1); \
		xptr(rbp+0*16) = carry0; \
		xptr(rbp+1*16) = carry1; \
		xptr(rbp+2*16) = carry2; \
		xptr(rbp+3*16) = carry3; \
		loopcount2--; \
		if(loopcount2 == 0) \
			break; \
		rsi += g->normblkdst; \
		rbp += 64; \
	ttp(rdx += 128); \
		loopcount3 += 0x80000000u / 4u; \
		if(loopcount3 == 0) \
			rsi += g->normblkdst8; /* Add 128 every 8 clmblkdsts */ \
	}while(1); \
	echk(g->u.xmm.XMM_MAXERR = xmm6); \
	\
	if(g->this_block == g->last_pass1_block) { \
		rsi = saved_rsi; \
		xnorm_top_carry; \
	} \
}

#define zpnorm(lab, ttp, echk, const1, base2, sse4, khi, c1, cm1) \
void CONCAT(lab, ARCH1)(struct gwasm_data *__restrict g, uintptr_t rsi) { \
	uintptr_t rax, rcx, rbx, rdi, rbp; \
ttp(uintptr_t rdx); \
	uintptr_t loopcount1z, loopcount2z; \
	unsigned int loopcount3z; \
	vec2f64 xmm0, xmm1, xmm2, xmm3, xmm4, xmm5, \
		xmm8, xmm9, xmm10, xmm11, xmm12, xmm13; \
\
	/* Handled in C code by pass1_pre_carries */ \
	/*zpad_sub7(g);*/		/* Subtract 7 ZPAD words from lowest FFT words */ \
echk(vec2f64 xmm6 = g->u.xmm.XMM_MAXERR); \
	rbx = (uintptr_t) g->norm_ptr2; \
ttp(rax = (uint32_t)g->cache_line_multiplier); \
	rdi = ((uintptr_t)g->u.xmm.XMM_COL_MULTS) + 128; /* LEA */ \
ttp(setlp:) \
	xnorm_2d_setup(ttp); \
ttp(rdi += 512); \
ttp(rbx += 32); \
ttp(rax -= 1); \
ttp(if(rax) goto setlp); \
	\
ttp(rdx = (uintptr_t)g->norm_grp_mults); \
	rbp = (uintptr_t)g->carries; \
	rdi = (uintptr_t)g->norm_ptr1; \
	rax = (uint32_t)g->addcount1; \
	loopcount2z = rax; \
	loopcount3z = 0; \
ttp(rax = u8ptr(rdi+0)); \
no##ttp(rax = 0); \
no##ttp(rcx = 0); \
	do{ /*ilp0:*/ \
		rbx = (uint32_t)g->cache_line_multiplier; \
		loopcount1z = rbx; \
		rbx = ((uintptr_t)g->u.xmm.XMM_COL_MULTS); /* LEA */ \
/* L2prefetch128 [rdx+128] */ \
		vec2f64 carry0 = xptr(rbp+0*16); \
		vec2f64 carry1 = xptr(rbp+1*16); \
		vec2f64 carry2 = xptr(rbp+2*16); \
		vec2f64 carry3 = xptr(rbp+3*16); \
		do{ /*ilp1:*/ \
			/* xprefetchw [rsi+64] */ \
			xnorm_2d_zpad(ttp, echk, const1, base2, sse4, khi, c1, cm1, carry0, carry1, carry2, carry3); \
			rsi += 64; \
			ttp(rbx += 512); \
			ttp(rdi += 4); \
			loopcount1z--; \
		}while(loopcount1z); \
		xptr(rbp+0*16) = carry0; \
		xptr(rbp+1*16) = carry1; \
		xptr(rbp+2*16) = carry2; \
		xptr(rbp+3*16) = carry3; \
		loopcount2z--; \
		if(loopcount2z == 0) \
			break; \
		rsi += g->normblkdst; \
		rbp += 64; \
	ttp(rdx += 128); \
		loopcount3z += 0x80000000u / 4u; \
		if(loopcount3z == 0) \
			rsi += g->normblkdst8; /* Add 128 every 8 clmblkdsts */ \
	}while(1); \
echk(g->u.xmm.XMM_MAXERR = xmm6); \
}
