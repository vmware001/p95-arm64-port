#include <stdlib.h>
#include <ctype.h>
#include <stdio.h>

#include "../gwnum.h"
#include "../gwcommon.h"
#include "../unravel.h"
#include "../gwtables.h"
#include "../xarch.h"
#include "../xbasics.h"
#include "xnormal.h"



typedef double __attribute__ ((vector_size (64), aligned(16))) vec8f64b;
typedef double __attribute__ ((vector_size (32), aligned(16))) vec4f64b;


void gwxaddq2(struct gwasm_data *__restrict g)
{
	char *rcx = g->SRCARG;
	char *rdx = g->SRC2ARG;
	char *rsi = g->DESTARG;
	uintptr_t rbx = g->addcount1;
	do{
		unsigned int a = g->normval4;
		do{
			unsigned int carry;
			do{
				*(vec8f64b*)rsi = *(vec8f64b*)rdx + *(vec8f64b*)rcx;
				rcx += 64;
				rdx += 64;
				rsi += 64;
				a = __builtin_addc(a, 0x80000000u / 64u, 0, &carry);
			}while(!carry);
			rcx += 128;
			rdx += 128;
			rsi += 128;
			a--;
		}while(a);
		rcx += g->pass2gapsize;
		rdx += g->pass2gapsize;
		rsi += g->pass2gapsize;
		rbx--;
	}while(rbx);
}
//DEST = SRC2 - SRC
void gwxsubq2(struct gwasm_data *__restrict g)
{
	char *rcx = g->SRCARG;
	char *rdx = g->SRC2ARG;
	char *rsi = g->DESTARG;
	uintptr_t rbx = g->addcount1;
	do{
		unsigned int a = g->normval4;
		do{
			unsigned int carry;
			do{
				*(vec8f64b*)rsi = *(vec8f64b*)rdx - *(vec8f64b*)rcx;
				rcx += 64;
				rdx += 64;
				rsi += 64;
				a = __builtin_addc(a, 0x80000000u / 64u, 0, &carry);
			}while(!carry);
			rcx += 128;
			rdx += 128;
			rsi += 128;
			a--;
		}while(a);
		rcx += g->pass2gapsize;
		rdx += g->pass2gapsize;
		rsi += g->pass2gapsize;
		rbx--;
	}while(rbx);
}

// DEST = SRC + SRC2, DEST2 = SRC - SRC2
void gwxaddsubq2(struct gwasm_data *__restrict g)
{
	char *rcx = g->SRCARG;
	char *rdx = g->SRC2ARG;
	char *rsi = g->DESTARG;
	char *rbp = g->DEST2ARG;
	uintptr_t rbx = g->addcount1;
	do{
		unsigned int a = g->normval4;
		do{
			unsigned int carry;
			do{
				vec4f64b x = *(vec4f64b*)rcx;
				vec4f64b y = *(vec4f64b*)rdx;
				*(vec4f64b*)rsi = x + y;
				*(vec4f64b*)rbp = x - y;
				vec4f64b x2 = ((vec4f64b*)rcx)[1];
				vec4f64b y2 = ((vec4f64b*)rdx)[1];
				((vec4f64b*)rsi)[1] = x2 + y2;
				((vec4f64b*)rbp)[1] = x2 - y2;
				rcx += 64;
				rdx += 64;
				rsi += 64;
				a = __builtin_addc(a, 0x80000000u / 64u, 0, &carry);
			}while(!carry);
			rcx += 128;
			rdx += 128;
			rsi += 128;
			a--;
		}while(a);
		rcx += g->pass2gapsize;
		rdx += g->pass2gapsize;
		rsi += g->pass2gapsize;
		rbx--;
	}while(rbx);
}



static void final_carries_2(struct gwasm_data *__restrict g, uintptr_t rsi, vec2f64 xmm6, vec2f64 xmm7) {
	uintptr_t rax,rbx,rcx,rdi,rbp;
	rbp = (uintptr_t)g->norm_grp_mults;	// Addr of the group multipliers
	rbx = (uintptr_t)g->norm_col_mults;	// Addr of the column multipliers
	rdi = (uintptr_t)g->norm_biglit_array;	// Addr of the big/little flags array

	if(g->ZERO_PADDED_FFT == 0)
	{
		xnorm_top_carry_cmn(rsi, xmm7[0], 2);
		if(g->B_IS_2 == 0) {
			xnorm_smallmul_2d_fft(noexec);
		}else{
			xnorm_smallmul_2d_fft(exec);
		};
	}else{
		if(g->B_IS_2 == 0) {
			xnorm_smallmul_2d_fft_zpad(noexec);
		}else{
			xnorm_smallmul_2d_fft_zpad(exec);
		};

	}
}

void gwxadd2(struct gwasm_data *__restrict g) {
	vec2f64 xmm0,xmm1,xmm2,xmm3,xmm4,xmm5,xmm6,xmm7;
	vec2f64 t1,t2,t3,t4;
	uintptr_t rax, rcx, rdx, rbx, rbp, rdi, rsi;
	unsigned int loopcount1, loopcount2, loopcount3, loopcount4, loopcount5;
	uintptr_t saved_reg;
	rcx = (uintptr_t)g->SRCARG;
	rdx = (uintptr_t)g->SRC2ARG;
	rsi = (uintptr_t)g->DESTARG;
	rbp = (uintptr_t)g->norm_grp_mults;
	rdi = (uintptr_t)g->norm_biglit_array;
	xmm7 = XMM_BIGVAL2;
	t1 = xmm7;
	t2 = xmm7;
	t3 = xmm7;
	t4 = xmm7;
	loopcount1 = g->count3;

	do{ // asec:
		loopcount2 = loopcount1 & 0x7FF;
		g->norm_ptr1 = (void*)rsi;
		g->norm_ptr2 = (void*)rbp;
		do{

			loopcount3 = g->normval4;
			rbx = (uintptr_t)g->norm_col_mults;
			do{
				if(g->B_IS_2 == 0)
				{
					if(g->RATIONAL_FFT != 0){
						loopcount4 = 128;
						do{
							xnorm_op_2d(addpd, noexec, noexec, saved_reg, t1,t2,t3,t4);
						}while(--loopcount4);
					}else{
						loopcount4 = g->normval1;
						do{
							loopcount5 = g->cache_line_multiplier;
							do{
								xnorm_op_2d(addpd, exec, noexec, saved_reg, t1,t2,t3,t4);
							}while(--loopcount5);
							rdi += g->normval2;
						}while(--loopcount4);
					}
				}else{
					if(g->RATIONAL_FFT != 0){
						loopcount4 = 128;
						do{
							xnorm_op_2d(addpd, noexec, exec, saved_reg, t1,t2,t3,t4);
						}while(--loopcount4);
					}else{
						loopcount4 = g->normval1;
						do{
							loopcount5 = g->cache_line_multiplier;
							do{
								xnorm_op_2d(addpd, exec, exec, saved_reg, t1,t2,t3,t4);
							}while(--loopcount5);
							rdi += g->normval2;
						}while(--loopcount4);
					}
				} // chunk done
				rcx += 128;
				rdx += 128;
				rsi += 128;
				loopcount3--;
			}while(loopcount3);
			rcx += g->pass2gapsize;
			rdx += g->pass2gapsize;
			rsi += g->pass2gapsize;
			// Block done
			rsi -= g->pass1blkdst;
			xnorm_op_2d_blk(rsi, rbp, t1,t2,t3,t4);
			rsi += g->pass1blkdst;
			if(g->RATIONAL_FFT == 0) {
				rbp += 128;
				rdi += g->normval3;
			}
			loopcount2--;
		}while(loopcount2);
		rax = (uintptr_t)g->norm_ptr1;
		rbx = (uintptr_t)g->norm_ptr2;
		xnorm_op_2d_sec(t1,t2,t3,t4);
		loopcount1 >>= 11;
	}while(loopcount1 != 0);
	rsi = (uintptr_t)g->DESTARG;
	xmm6 = t1;
	xmm7 = t3;
	final_carries_2(g, rsi, xmm6, xmm7);
}

void gwxsub2(struct gwasm_data *__restrict g) {
	vec2f64 xmm0,xmm1,xmm2,xmm3,xmm4,xmm5,xmm6,xmm7;
	vec2f64 t1,t2,t3,t4;
	uintptr_t rax, rcx, rdx, rbx, rbp, rdi, rsi;
	unsigned int loopcount1, loopcount2, loopcount3, loopcount4, loopcount5;
	uintptr_t saved_reg;
	rcx = (uintptr_t)g->SRCARG;
	rdx = (uintptr_t)g->SRC2ARG;
	rsi = (uintptr_t)g->DESTARG;
	rbp = (uintptr_t)g->norm_grp_mults;
	rdi = (uintptr_t)g->norm_biglit_array;
	xmm7 = XMM_BIGVAL2;
	t1 = xmm7;
	t2 = xmm7;
	t3 = xmm7;
	t4 = xmm7;
	loopcount1 = g->count3;

	do{ // asec:
		loopcount2 = loopcount1 & 0x7FF;
		g->norm_ptr1 = (void*)rsi;
		g->norm_ptr2 = (void*)rbp;
		do{

			loopcount3 = g->normval4;
			rbx = (uintptr_t)g->norm_col_mults;
			do{
				if(g->B_IS_2 == 0)
				{
					if(g->RATIONAL_FFT != 0){
						loopcount4 = 128;
						do{
							xnorm_op_2d(subpd, noexec, noexec, saved_reg, t1,t2,t3,t4);
						}while(--loopcount4);
					}else{
						loopcount4 = g->normval1;
						do{
							loopcount5 = g->cache_line_multiplier;
							do{
								xnorm_op_2d(subpd, exec, noexec, saved_reg, t1,t2,t3,t4);
							}while(--loopcount5);
							rdi += g->normval2;
						}while(--loopcount4);
					}
				}else{
					if(g->RATIONAL_FFT != 0){
						loopcount4 = 128;
						do{
							xnorm_op_2d(subpd, noexec, exec, saved_reg, t1,t2,t3,t4);
						}while(--loopcount4);
					}else{
						loopcount4 = g->normval1;
						do{
							loopcount5 = g->cache_line_multiplier;
							do{
								xnorm_op_2d(subpd, exec, exec, saved_reg, t1,t2,t3,t4);
							}while(--loopcount5);
							rdi += g->normval2;
						}while(--loopcount4);
					}
				} // chunk done
				rcx += 128;
				rdx += 128;
				rsi += 128;
				loopcount3--;
			}while(loopcount3);
			rcx += g->pass2gapsize;
			rdx += g->pass2gapsize;
			rsi += g->pass2gapsize;
			// Block done
			rsi -= g->pass1blkdst;
			xnorm_op_2d_blk(rsi, rbp, t1,t2,t3,t4);
			rsi += g->pass1blkdst;
			if(g->RATIONAL_FFT == 0) {
				rbp += 128;
				rdi += g->normval3;
			}
			loopcount2--;
		}while(loopcount2);
		rax = (uintptr_t)g->norm_ptr1;
		rbx = (uintptr_t)g->norm_ptr2;
		xnorm_op_2d_sec(t1,t2,t3,t4);
		loopcount1 >>= 11;
	}while(loopcount1 != 0);
	rsi = (uintptr_t)g->DESTARG;
	xmm6 = t1;
	xmm7 = t3;
	final_carries_2(g, rsi, xmm6, xmm7);
}


#if 0 // Alternative method, slower
void gwxaddsub2(struct gwasm_data *__restrict g) {
	void *dest1arg = g->DESTARG;
	gwxadd2(g);
	g->DESTARG = g->DEST2ARG;
	gwxsub2(g);
	g->DESTARG = dest1arg;
	return;
}
#else // Primary method
void gwxaddsub2(struct gwasm_data *__restrict g) {
	vec2f64 xmm0,xmm1,xmm2,xmm3,xmm4,xmm5,xmm6,xmm7;
	vec2f64 t1,t2,t3,t4,t5,t6,t7,t8;
	uintptr_t rax, rcx, rdx, rbx, rbp, rdi, rsi;
	unsigned int loopcount1, loopcount2, loopcount3, loopcount4, loopcount5;
	uintptr_t saved_grp_ptr, saved_dest1_ptr, saved_dest2_ptr, saved_reg;
	rcx = (uintptr_t)g->SRCARG;
	rdx = (uintptr_t)g->SRC2ARG;
	rsi = (uintptr_t)g->DESTARG;
	rbp = (uintptr_t)g->DEST2ARG;
	rbx = (uintptr_t)g->norm_grp_mults;
	rdi = (uintptr_t)g->norm_biglit_array;
	xmm7 = XMM_BIGVAL2;
	t1 = xmm7;
	t2 = xmm7;
	t3 = xmm7;
	t4 = xmm7;
	t5 = xmm7;
	t6 = xmm7;
	t7 = xmm7;
	t8 = xmm7;
	loopcount1 = g->count3;

	do{ // asec:
		loopcount2 = loopcount1 & 0x7FF;
		saved_dest1_ptr = rsi;
		saved_dest2_ptr = rbp;
		saved_grp_ptr = rbx;
		do{
			loopcount3 = g->normval4;
			rax = (uintptr_t)g->norm_col_mults;
			do{
				if(g->B_IS_2 == 0)
				{
					if(g->RATIONAL_FFT != 0){
						loopcount4 = 128;
						do{
							xnorm_addsub_2d(noexec, noexec, saved_reg, t1,t2,t3,t4,t5,t6,t7,t8);
						}while(--loopcount4);
					}else{
						loopcount4 = g->normval1;
						do{
							loopcount5 = g->cache_line_multiplier;
							do{
								xnorm_addsub_2d(exec, noexec, saved_reg, t1,t2,t3,t4,t5,t6,t7,t8);
							}while(--loopcount5);
							rdi += g->normval2;
						}while(--loopcount4);
					}
				}else{
					if(g->RATIONAL_FFT != 0){
						loopcount4 = 128;
						do{
							xnorm_addsub_2d(noexec, exec, saved_reg, t1,t2,t3,t4,t5,t6,t7,t8);
						}while(--loopcount4);
					}else{
						loopcount4 = g->normval1;
						do{
							loopcount5 = g->cache_line_multiplier;
							do{
								xnorm_addsub_2d(exec, exec, saved_reg, t1,t2,t3,t4,t5,t6,t7,t8);
							}while(--loopcount5);
							rdi += g->normval2;
						}while(--loopcount4);
					}
				} // chunk done
				rcx += 128;
				rdx += 128;
				rsi += 128;
				rbp += 128;
				loopcount3--;
			}while(loopcount3);
			rcx += g->pass2gapsize;
			rdx += g->pass2gapsize;
			rsi += g->pass2gapsize;
			rbp += g->pass2gapsize;

			// Block done
			rsi -= g->pass1blkdst;
			xnorm_op_2d_blk(rsi, rbx, t1,t2,t3,t4);
			rbp -= g->pass1blkdst;
			xnorm_op_2d_blk(rbp, rbx, t5,t6,t7,t8);
			rsi += g->pass1blkdst;
			rbp += g->pass1blkdst;
			if(g->RATIONAL_FFT == 0) {
				rbx += 128;
				rdi += g->normval3;
			}
			loopcount2--;
		}while(loopcount2);
		saved_reg = rbx;
		rax = saved_dest1_ptr;
		rbx = saved_grp_ptr;
		xnorm_op_2d_sec(t1,t2,t3,t4);
		rax = saved_dest2_ptr;
		xnorm_op_2d_sec(t5,t6,t7,t8);
		rbx = saved_reg;
		loopcount1 >>= 11;
	}while(loopcount1);
	// xptr(g->u.xmm.XMM_TMP8) = xptr(g->u.xmm.XMM_TMP5);
	rsi = (uintptr_t)g->DESTARG;
	xmm6 = t1;
	xmm7 = t3;
	final_carries_2(g, rsi, xmm6, xmm7);
	rsi = (uintptr_t)g->DEST2ARG;
	xmm6 = t5;
	xmm7 = t7;
	final_carries_2(g, rsi, xmm6, xmm7);
}
#endif




void gwxmuls2(struct gwasm_data *__restrict g) {
	vec2f64 xmm0,xmm1,xmm2,xmm3,xmm4,xmm5,xmm6,xmm7;
	vec2f64 t1,t2,t3,t4,t5;
	uintptr_t rax, rcx, rbx, rbp, rdi, rsi, tmp1;
	unsigned int loopcount1, loopcount2, loopcount3, loopcount4, loopcount5;
	uintptr_t saved_sec_biglit, saved_blk_start, saved_blk_biglit;
	rsi = (uintptr_t)g->DESTARG;
	rbp = (uintptr_t)g->norm_grp_mults;
	rdi = (uintptr_t)g->norm_biglit_array;
	xmm0[1] = xmm0[0] = g->DBLARG;
	if(g->RATIONAL_FFT == 0)
	{
		xmm0 *= g->u.xmm.XMM_NORM012_FF;
	}
	t5 = xmm0;


	xmm7 = XMM_BIGVAL2;
	t1 = xmm7;
	t2 = xmm7;
	t3 = xmm7;
	t4 = xmm7;
	loopcount1 = g->count3;

	do{ // asec:
		loopcount2 = loopcount1 & 0x7FF;
		g->norm_ptr1 = (void*)rsi;
		g->norm_ptr2 = (void*)rbp;
		saved_sec_biglit = rdi;
		do{
			saved_blk_start = rsi;
			saved_blk_biglit = rdi;
			loopcount3 = g->normval4;
			rbx = (uintptr_t)g->norm_col_mults;
			do{
				if(g->B_IS_2 == 0)
				{
					if(g->RATIONAL_FFT != 0){
						loopcount4 = 128;
						do{
							xnorm_smallmul_2d(noexec, noexec, t1,t2,t3,t4,t5);
						}while(--loopcount4);
					}else{
						loopcount4 = g->normval1;
						do{
							loopcount5 = g->cache_line_multiplier;
							do{
								xnorm_smallmul_2d(exec, noexec, t1,t2,t3,t4,t5);
							}while(--loopcount5);
							rdi += g->normval2;
						}while(--loopcount4);
					}
				}else{
					if(g->RATIONAL_FFT != 0){
						loopcount4 = 128;
						do{
							xnorm_smallmul_2d(noexec, exec, t1,t2,t3,t4,t5);
						}while(--loopcount4);
					}else{
						loopcount4 = g->normval1;
						do{
							loopcount5 = g->cache_line_multiplier;
							do{
								xnorm_smallmul_2d(exec, exec, t1,t2,t3,t4,t5);
							}while(--loopcount5);
							rdi += g->normval2;
						}while(--loopcount4);
					}
				} // Chunk done
				rsi += 128;
				loopcount3--;
			}while(loopcount3);
			// Block done
			tmp1 = rsi; rsi = saved_blk_start; saved_blk_start = tmp1;
			tmp1 = rdi; rdi = saved_blk_biglit; saved_blk_biglit = tmp1;
			rbx = (uintptr_t)g->norm_col_mults;
			if(g->B_IS_2 == 0){
				xnorm_smallmul_2d_blk(noexec, t1,t2,t3,t4);
			}else{
				xnorm_smallmul_2d_blk(exec, t1,t2,t3,t4);
			}
			rsi = saved_blk_start;
			rdi = saved_blk_biglit;
			rsi += g->pass2gapsize;
			if(g->RATIONAL_FFT == 0) {
				rdi += g->normval3;
				rbp += 128;
			}
			loopcount2--;
		}while(loopcount2);
		// Section done
		tmp1 = rsi; rsi = (intptr_t)g->norm_ptr1; g->norm_ptr1 = (void*)tmp1;
		tmp1 = rbp; rbp = (intptr_t)g->norm_ptr2; g->norm_ptr2 = (void*)tmp1;
		tmp1 = rdi; rdi = saved_sec_biglit; saved_sec_biglit = tmp1;
		rbx = (uintptr_t)g->norm_col_mults;
		if(g->B_IS_2 == 0){
			xnorm_smallmul_2d_sec(noexec, t1,t2,t3,t4);
		}else{
			xnorm_smallmul_2d_sec(exec, t1,t2,t3,t4);
		}
		rsi = (intptr_t)g->norm_ptr1;
		rbp = (intptr_t)g->norm_ptr2;
		rdi = saved_sec_biglit;
		loopcount1 >>= 11;
	}while(loopcount1 != 0);
	rsi = (uintptr_t)g->DESTARG;
	xmm6 = t1;
	xmm7 = t3;
	final_carries_2(g, rsi, xmm6, xmm7);
}

