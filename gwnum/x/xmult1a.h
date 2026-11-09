#pragma once
//#include <stdio.h>

#define inorm(lab, ttp, echk, const1, base2, sse4) \
void lab##BLEND(struct gwasm_data *__restrict g) { \
	uintptr_t rdx, rsi, rbp, tmp1; \
	vec2f64 xmm0, xmm1, xmm2, xmm3, xmm4, xmm5; \
	rsi = (uintptr_t)g->DESTARG; \
	/*static int TT=0; if(!TT) {TT=1; puts("inorm1 " #lab","#ttp","#echk","#const1","#base2","#sse4);}*/ \
	f64ptr(rsi+g->ADDIN_OFFSET) += g->ADDIN_VALUE; \
	xmm2 = XMM_BIGVAL2; \
	xmm3 = xmm2; \
echk(vec2f64 xmm6 = {0.0, 0.0}); \
echk(xmm6[0] = g->MAXERR); \
	rbp = (uintptr_t)g->norm_col_mults; \
ttp(uintptr_t rdi = (uintptr_t)g->norm_biglit_array); \
	/*rax = 0;*/ \
	rdx = (uint32_t)g->normcount1; /* 32 bit unsigned int */ \
	uintptr_t saved_reg3 = rsi; \
ttp(uintptr_t saved_reg2 = rdi); \
ttp(uintptr_t saved_reg1 = rbp); \
	do{ \
		unsigned int ebx = rdx & 0x7FF; \
		do{ \
			xnorm_1d(ttp, echk, const1, base2, sse4); \
			rsi += 64; \
			ttp(rbp += 128); \
			ttp(rdi += 4); \
			ebx -= 1; \
		}while(ebx != 0); \
		tmp1 = rsi; rsi = saved_reg3; saved_reg3 = tmp1; \
		ttp(tmp1 = rdi); ttp(rdi = saved_reg2); ttp(saved_reg2 = tmp1); \
		ttp(tmp1 = rbp); ttp(rbp = saved_reg1); ttp(saved_reg1 = tmp1); \
		xnorm012_1d_mid(ttp, base2); \
		rsi = saved_reg3; \
		ttp(rdi = saved_reg2); \
		ttp(rbp = saved_reg1); \
		rdx >>= 11; \
	}while(rdx != 0); \
	\
/*echk(printf("----- e2 %f %f\n", xmm6[0], xmm6[1]));*/ \
/* printf("err %f\n", xmm6[0], xmm6[1]); */ \
echk(g->MAXERR = vec2reducemax(xmm6)); \
	rsi = (uintptr_t)g->DESTARG; \
	uintptr_t rdi2 = (uintptr_t)g->ADDIN_OFFSET; \
	/* Unlike AVX, FMA, and AVX512 FFTs, one pass SSE2 FFT's ADDIN_OFFSET can be different than POST_ADDIN_OFFSET */ \
	/* Look at if this is the 1st, 2nd, 3rd, or 4th double in a group of 4 doubles */ \
	rdx = rdi2 & 0b11000; \
	if(rdx != 0 && rdx != 0b11000) \
		rdi2 ^= 0b11000; \
	f64ptr(rsi+rdi2) += g->POSTADDIN_VALUE; \
	inorm_end_##base2(g, xmm2, xmm3); \
}
#define inorm_end_noexec non2dn
#define inorm_end_exec idn

#define zpnorm(lab, ttp, echk, const1, base2, sse4, khi, c1, cm1) \
void lab##BLEND(struct gwasm_data *__restrict g) { \
	/*static int TT=0; if(!TT) {TT=1; puts("zpnorm1 " #lab","#ttp","#echk","#const1"," #base2","#sse4","#khi","#c1","#cm1);}*/ \
	uintptr_t rdx, rbx, rsi, rdi, rbp, tmp1; \
	uintptr_t saved_reg3, saved_reg2, saved_reg1; \
	vec2f64 xmm0, xmm1, xmm2, xmm3, xmm4, xmm5; \
	zpad_sub7(g);		/* Subtract 7 ZPAD words from lowest FFT words */ \
	rsi = (uintptr_t)g->DESTARG;		/* Addr of multiplied number */ \
	xmm2 = XMM_BIGVAL2; \
	xmm3[1] = xmm3[0] = 0.0; \
	echk(vec2f64 xmm6 = {0.0, 0.0}); \
	echk(xmm6[0] = g->MAXERR); \
	rbp = (uintptr_t)g->norm_col_mults;	/* Addr of the multipliers */ \
	rdi = (uintptr_t)g->norm_biglit_array;	/* Addr of the big/little flags array */ \
	/*rax = 0;*/		/* Clear big/lit flag */ \
	rbx = (uint32_t)g->normcount1;		/* Load loop counter */ \
	do{ \
		saved_reg1 = rbx;		/* Save loop counter */ \
		rbx &= 0x7FF;		/* Grab 11 bits of the counter */ \
		saved_reg2 = rdi;		/* remember edi for xnorm012_1d_mid */ \
		saved_reg3 = rsi;		/* remember esi for xnorm012_1d_mid */ \
		rdx = rbp;		/* remember ebp for xnorm012_1d_mid */ \
		do{ \
			xnorm_1d_zpad(ttp, echk, const1, base2, sse4, khi, c1, cm1); /* Normalize 8 values */ \
			rsi += 64;			/* Next cache line */ \
			ttp(rbp += 128);		/* Next set of 8 multipliers */ \
			ttp(rdi += 4);			/* Next big/little flags */ \
			rbx -= 1;			/* Test loop counter */ \
		}while(rbx != 0);			/* Loop til done */ \
		rbx = saved_reg3;		/* Restore FFT data addr */ \
		tmp1 = rdi; rdi = saved_reg2; saved_reg2 = tmp1; \
		xnorm012_1d_mid_zpad(const1, base2); /* Rotate carries and add in carries */ \
		rdi = saved_reg2;		/* Restore big/lit pointer */ \
		rbx = saved_reg1;		/* Restore loop counter */ \
		rbx >>= 11;			/* Get next loop amount */ \
	}while(rbx != 0); \
	\
/*echk(printf("--- err %f %f\n", xmm6[0], xmm6[1]));*/ \
echk(g->MAXERR = vec2reducemax(xmm6)); \
	zpnorm_end_##base2##const1(g, xmm2, xmm3); \
}

#define zpnorm_end_noexecexec non2zpcdn
#define zpnorm_end_noexecnoexec non2zpdn
#define zpnorm_end_execexec zpcdn
#define zpnorm_end_execnoexec zpdn




void non2zpdn(struct gwasm_data *__restrict g, vec2f64 xmm2, vec2f64 xmm3);
void non2zpcdn(struct gwasm_data *__restrict g, vec2f64 xmm2, vec2f64 xmm3);
void zpcdn(struct gwasm_data *__restrict g, vec2f64 xmm2, vec2f64 xmm3);
void zpdn(struct gwasm_data *__restrict g, vec2f64 xmm2, vec2f64 xmm3);
void non2dn(struct gwasm_data *__restrict g, vec2f64 xmm2, vec2f64 xmm3);
void zdn(struct gwasm_data *__restrict g, vec2f64 xmm2, vec2f64 xmm3);
void idn(struct gwasm_data *__restrict g, vec2f64 xmm2, vec2f64 xmm3);
