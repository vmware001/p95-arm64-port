#pragma once

#ifndef ARCH1
#define ARCH1 CORE
#endif




// From: pass1_normalize, in xmult.h; which is from hgpass1.h

#define inorm(lab, ttp, echk, const1, base2, sse4) \
void CONCAT(lab, ARCH1)(struct gwasm_data *__restrict g, uintptr_t rsi) { \
	uintptr_t rdi, rbp; \
	ttp(uintptr_t rdx); \
	ttp(uintptr_t rbx); \
	unsigned int blk8_counter; \
	vec2f64 xmm0, xmm1, xmm2, xmm3, xmm4, xmm5, xmm7; \
	vec2f64 xmm8, xmm9, xmm10, xmm11, xmm12, xmm13; \
	base2(const1(no##echk(vec2f64 xmm6, xmm15))); \
	base2(vec2f64 xmm14); \
	uintptr_t saved_rsi = rsi; \
	xnorm_wpn_preload(ttp, echk, const1, base2, sse4); \
echk(vec2f64 xmm6 = g->u.xmm.XMM_MAXERR); \
ttp(rdx = (uintptr_t)g->norm_grp_mults); \
	rbp = (uintptr_t)g->carries; \
	rdi = (uintptr_t)g->norm_ptr1; \
	blk8_counter = 0; \
	uint32_t loopcount3 = g->count3; \
	no##ttp(loopcount3 *= g->count2;) \
ttp(rbx = u16ptr(rdi)); \
	do{ /*ilp0:*/ \
		ttp(uint32_t loopcount2 = g->count2); \
		ttp(do){ /*ilp1:*/ \
			uint32_t loopcount1 = g->cache_line_multiplier; \
			vec2f64 carry0 = xptr(rbp+0*16); \
			vec2f64 carry1 = xptr(rbp+1*16); \
			vec2f64 carry2 = xptr(rbp+2*16); \
			vec2f64 carry3 = xptr(rbp+3*16); \
			do{ /*ilp2:*/ \
				/*xprefetchw [rsi+64]*/; \
				xnorm_wpn(ttp, echk, const1, base2, sse4, carry0, carry1, carry2, carry3); \
				rsi += 64; \
				rdi += 2; \
				loopcount1--; \
			}while(loopcount1); \
			xptr(rbp+0*16) = carry0; \
			xptr(rbp+1*16) = carry1; \
			xptr(rbp+2*16) = carry2; \
			xptr(rbp+3*16) = carry3; \
			rsi += g->normblkdst; /* Skip gap in blkdst or clmblkdst */ \
			rbp += 64; /* Next set of carries */ \
			blk8_counter += 0x20000000u; /* 2^32 / 8 */ \
			if(blk8_counter == 0) { \
				rsi += g->normblkdst8; \
			} \
		}ttp(while(--loopcount2)); \
		ttp(rdx += XMM_GMD*4); \
		loopcount3--; \
	}while(loopcount3); \
echk(g->u.xmm.XMM_MAXERR = xmm6); \
	\
	if(g->this_block == g->last_pass1_block) { \
		rsi = saved_rsi; \
		xnorm_top_carry; \
	} \
}

#define zpnorm(lab, ttp, echk, const1, base2, sse4, khi, c1, cm1) \
void CONCAT(lab, ARCH1)(struct gwasm_data *__restrict g, uintptr_t rsi) { \
	uintptr_t rbp; \
	ttp(uintptr_t rdx, rbx, rdi); \
	unsigned int blk8_counterz; \
	vec2f64 xmm0, xmm1, xmm2, xmm3, xmm4, xmm5, xmm7; \
	vec2f64 xmm8, xmm9, xmm10, xmm11, xmm12, xmm15; \
\
	/* Handled in C code by pass1_pre_carries */ \
	/*zpad_sub7(g);*/		/* Subtract 7 ZPAD words from lowest FFT words */ \
	no##const1(no##c1(no##cm1(vec2f64 xmm14))); \
	no##const1(    c1(vec2f64 xmm14)); \
	const1(vec2f64 xmm14); \
	base2(no##echk(sse4(khi(vec2f64 xmm6)))); \
	no##base2(no##echk(khi(vec2f64 xmm6))); \
\
	xnorm_wpn_zpad_preload(ttp, echk, const1, base2, sse4, khi, c1, cm1); \
echk(vec2f64 xmm6 = g->u.xmm.XMM_MAXERR); \
ttp(rdx = (uintptr_t)g->norm_grp_mults); \
	rbp = (uintptr_t)g->carries; \
ttp(rdi = (uintptr_t)g->norm_ptr1); \
	blk8_counterz = 0; \
	uint32_t loopcount3z = g->count3; \
	no##ttp(loopcount3z *= g->count2;) \
ttp(rbx = u16ptr(rdi)); \
	do{ /*ilp0:*/ \
		ttp(uint32_t loopcount2z = g->count2); \
		ttp(do){ /*ilp1:*/ \
			uint32_t loopcount1z = g->cache_line_multiplier; \
			vec2f64 carry0 = xptr(rbp+0*16); \
			vec2f64 carry1 = xptr(rbp+1*16); \
			vec2f64 carry2 = xptr(rbp+2*16); \
			vec2f64 carry3 = xptr(rbp+3*16); \
			do{ /*ilp2:*/ \
				/*xprefetchw [rsi+64]*/; \
				xnorm_wpn_zpad(ttp, echk, const1, base2, sse4, khi, c1, cm1, carry0, carry1, carry2, carry3); \
				rsi += 64; \
				ttp(rdi += 2); \
				loopcount1z--; \
			}while(loopcount1z); \
			xptr(rbp+0*16) = carry0; \
			xptr(rbp+1*16) = carry1; \
			xptr(rbp+2*16) = carry2; \
			xptr(rbp+3*16) = carry3; \
			rsi += g->normblkdst; \
			rbp += 64; \
			blk8_counterz += 0x20000000u; /* 2^32 / 8 */ \
			if(blk8_counterz == 0) { \
				rsi += g->normblkdst8; \
			} \
		}ttp(while(--loopcount2z)); \
		ttp(rdx += XMM_GMD*4); \
		loopcount3z--; \
	}while(loopcount3z); \
echk(g->u.xmm.XMM_MAXERR = xmm6); \
}
