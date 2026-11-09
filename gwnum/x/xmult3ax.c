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
#include "xnormal.h"

static void final_carries_3(struct gwasm_data *__restrict g, uintptr_t rsi, vec2f64 xmm6, vec2f64 xmm7) {
	uintptr_t rax,rcx,rdi,rbp;
	rbp = (uintptr_t)g->norm_grp_mults;	// Addr of the group multipliers
	//rbx = (uintptr_t)g->norm_col_mults;	// Addr of the column multipliers
	rdi = (uintptr_t)g->norm_biglit_array;	// Addr of the big/little flags array

	if(g->ZERO_PADDED_FFT == 0)
	{
		xnorm_top_carry_cmn(rsi, xmm7[0], 2);
		if(g->B_IS_2 == 0) {
			xnorm_smallmul_wpn_fft(noexec);
		}else{
			xnorm_smallmul_wpn_fft(exec);
		}
	}else{
		if(g->B_IS_2 == 0) {
			xnorm_smallmul_wpn_fft_zpad(noexec);
		}else{
			xnorm_smallmul_wpn_fft_zpad(exec);
		}
	}
}

void gwxadd3(struct gwasm_data *__restrict g) {
	vec2f64 xmm0,xmm1,xmm2,xmm3,xmm4,xmm5,xmm6,xmm7;
	vec2f64 t1,t2,t3,t4;
	uintptr_t rax, rcx, rdx, rbx, rbp, rdi, rsi, tmp1;
	unsigned int loopcount1, loopcount2, loopcount3, loopcount4, loopcount5;
	uintptr_t saved_blk_start, saved_blk_biglit;
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
	do{
		loopcount2 = g->count2;
		do{
			saved_blk_start = rsi;
			saved_blk_biglit = rdi;
			loopcount3 = g->normval4;
			do{
				if(g->B_IS_2 == 0)
				{
					if(g->RATIONAL_FFT != 0){
						loopcount4 = 128;
						do{
							xnorm_op_wpn(addpd, noexec, noexec, t1,t2,t3,t4);
						}while(--loopcount4);
					}else{
						loopcount4 = g->normval1;
						do{
							loopcount5 = g->cache_line_multiplier;
							do{
								xnorm_op_wpn(addpd, exec, noexec, t1,t2,t3,t4);
							}while(--loopcount5);
							rdi += g->normval2;
						}while(--loopcount4);
					}
				}else{
					if(g->RATIONAL_FFT != 0){
						loopcount4 = 128;
						do{
							xnorm_op_wpn(addpd, noexec, exec, t1,t2,t3,t4);
						}while(--loopcount4);
					}else{
						loopcount4 = g->normval1;
						do{
							loopcount5 = g->cache_line_multiplier;
							do{
								xnorm_op_wpn(addpd, exec, exec, t1,t2,t3,t4);
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
			tmp1 = rsi; rsi = saved_blk_start; saved_blk_start = tmp1;
			tmp1 = rdi; rdi = saved_blk_biglit; saved_blk_biglit = tmp1;
			xnorm_op_wpn_blk(rsi, rbp, t1,t2,t3,t4);
			rsi = saved_blk_start;
			rdi = saved_blk_biglit;
			if(g->RATIONAL_FFT == 0) {
				rdi += g->normval3;
			}
			loopcount2--;
		}while(loopcount2);
		if(g->RATIONAL_FFT == 0) {
			rbp += 4*XMM_GMD;
		}
		loopcount1--;
	}while(loopcount1);
	rax = (uintptr_t)g->DESTARG;
	rbx = (uintptr_t)g->norm_grp_mults;
	rdi = (uintptr_t)g->norm_biglit_array;
	xnorm_op_wpn_sec(t1,t2,t3,t4);
	rsi = (uintptr_t)g->DESTARG;
	xmm6 = t1;
	xmm7 = t3;
	final_carries_3(g, rsi, xmm6, xmm7);
}

void gwxsub3(struct gwasm_data *__restrict g) {
	vec2f64 xmm0,xmm1,xmm2,xmm3,xmm4,xmm5,xmm6,xmm7;
	vec2f64 t1,t2,t3,t4;
	uintptr_t rax, rcx, rdx, rbx, rbp, rdi, rsi, tmp1;
	unsigned int loopcount1, loopcount2, loopcount3, loopcount4, loopcount5;
	uintptr_t saved_blk_start, saved_blk_biglit;
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
	do{
		loopcount2 = g->count2;
		do{
			saved_blk_start = rsi;
			saved_blk_biglit = rdi;
			loopcount3 = g->normval4;
			do{
				if(g->B_IS_2 == 0)
				{
					if(g->RATIONAL_FFT != 0){
						loopcount4 = 128;
						do{
							xnorm_op_wpn(subpd, noexec, noexec, t1,t2,t3,t4);
						}while(--loopcount4);
					}else{
						loopcount4 = g->normval1;
						do{
							loopcount5 = g->cache_line_multiplier;
							do{
								xnorm_op_wpn(subpd, exec, noexec, t1,t2,t3,t4);
							}while(--loopcount5);
							rdi += g->normval2;
						}while(--loopcount4);
					}
				}else{
					if(g->RATIONAL_FFT != 0){
						loopcount4 = 128;
						do{
							xnorm_op_wpn(subpd, noexec, exec, t1,t2,t3,t4);
						}while(--loopcount4);
					}else{
						loopcount4 = g->normval1;
						do{
							loopcount5 = g->cache_line_multiplier;
							do{
								xnorm_op_wpn(subpd, exec, exec, t1,t2,t3,t4);
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
			tmp1 = rsi; rsi = saved_blk_start; saved_blk_start = tmp1;
			tmp1 = rdi; rdi = saved_blk_biglit; saved_blk_biglit = tmp1;
			xnorm_op_wpn_blk(rsi, rbp, t1,t2,t3,t4);
			rsi = saved_blk_start;
			rdi = saved_blk_biglit;
			if(g->RATIONAL_FFT == 0) {
				rdi += g->normval3;
			}
			loopcount2--;
		}while(loopcount2);
		if(g->RATIONAL_FFT == 0) {
			rbp += 4*XMM_GMD;
		}
		loopcount1--;
	}while(loopcount1);
	rax = (uintptr_t)g->DESTARG;
	rbx = (uintptr_t)g->norm_grp_mults;
	rdi = (uintptr_t)g->norm_biglit_array;
	xnorm_op_wpn_sec(t1,t2,t3,t4);
	rsi = (uintptr_t)g->DESTARG;
	xmm6 = t1;
	xmm7 = t3;
	final_carries_3(g, rsi, xmm6, xmm7);
}

void gwxaddsub3(struct gwasm_data *__restrict g) {
	vec2f64 xmm0,xmm1,xmm2,xmm3,xmm4,xmm5,xmm6,xmm7;
	vec2f64 t1,t2,t3,t4,t5,t6,t7,t8;
	uintptr_t rax, rcx, rdx, rbx, rbp, rdi, rsi, tmp1;
	unsigned int loopcount1, loopcount2, loopcount3, loopcount4, loopcount5;
	uintptr_t saved_blk_start, saved_blk_start2, saved_blk_biglit;
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
	do{
		loopcount2 = g->count2;
		do{
			saved_blk_start = rsi;
			saved_blk_start2 = rbp;
			saved_blk_biglit = rdi;
			loopcount3 = g->normval4;
			do{
				if(g->B_IS_2 == 0)
				{
					if(g->RATIONAL_FFT != 0){
						loopcount4 = 128;
						do{
							xnorm_addsub_wpn(noexec, noexec, t1,t2,t3,t4,t5,t6,t7,t8);
						}while(--loopcount4);
					}else{
						loopcount4 = g->normval1;
						do{
							loopcount5 = g->cache_line_multiplier;
							do{
								xnorm_addsub_wpn(exec, noexec, t1,t2,t3,t4,t5,t6,t7,t8);
							}while(--loopcount5);
							rdi += g->normval2;
						}while(--loopcount4);
					}
				}else{
					if(g->RATIONAL_FFT != 0){
						loopcount4 = 128;
						do{
							xnorm_addsub_wpn(noexec, exec, t1,t2,t3,t4,t5,t6,t7,t8);
						}while(--loopcount4);
					}else{
						loopcount4 = g->normval1;
						do{
							loopcount5 = g->cache_line_multiplier;
							do{
								xnorm_addsub_wpn(exec, exec, t1,t2,t3,t4,t5,t6,t7,t8);
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
			tmp1 = rsi; rsi = saved_blk_start; saved_blk_start = tmp1;
			tmp1 = rbp; rbp = saved_blk_start2; saved_blk_start2 = tmp1;
			tmp1 = rdi; rdi = saved_blk_biglit; saved_blk_biglit = tmp1;
			xnorm_op_wpn_blk(rsi, rbx, t1,t2,t3,t4);
			xnorm_op_wpn_blk(rbp, rbx, t5,t6,t7,t8);
			rsi = saved_blk_start;
			rbp = saved_blk_start2;
			rdi = saved_blk_biglit;
			if(g->RATIONAL_FFT == 0) {
				rdi += g->normval3;
			}
			loopcount2--;
		}while(loopcount2);
		if(g->RATIONAL_FFT == 0) {
			rbx += 4*XMM_GMD;
		}
		loopcount1--;
	}while(loopcount1);
	rax = (uintptr_t)g->DESTARG;
	rbx = (uintptr_t)g->norm_grp_mults;
	rdi = (uintptr_t)g->norm_biglit_array;
	xnorm_op_wpn_sec(t1,t2,t3,t4);
	rax = (uintptr_t)g->DEST2ARG;
	xnorm_op_wpn_sec(t5,t6,t7,t8);
	rsi = (uintptr_t)g->DESTARG;
	//xptr(g->u.xmm.XMM_TMP8) = xptr(g->u.xmm.XMM_TMP5);
	xmm6 = t1;
	xmm7 = t3;
	final_carries_3(g, rsi, xmm6, xmm7);
	rsi = (uintptr_t)g->DEST2ARG;
	xmm6 = t5;
	xmm7 = t7;
	final_carries_3(g, rsi, xmm6, xmm7);
}


void gwxmuls3(struct gwasm_data *__restrict g) {
	vec2f64 xmm0,xmm1,xmm2,xmm3,xmm4,xmm5,xmm6,xmm7;
	vec2f64 t1,t2,t3,t4,t5;
	uintptr_t rax, rcx, rbp, rdi, rsi, tmp1;
	unsigned int loopcount1, loopcount2, loopcount3, loopcount4, loopcount5;
	uintptr_t saved_blk_start, saved_blk_biglit;
	rsi = (uintptr_t)g->DESTARG;
	rbp = (uintptr_t)g->norm_grp_mults;
	rdi = (uintptr_t)g->norm_biglit_array;
	xmm0[1] = xmm0[0] = g->DBLARG;
	t5 = xmm0;
	xmm7 = XMM_BIGVAL2;
	t1 = xmm7;
	t2 = xmm7;
	t3 = xmm7;
	t4 = xmm7;
	loopcount1 = g->count3;
	do{
		loopcount2 = g->count2;
		do{
			saved_blk_start = rsi;
			saved_blk_biglit = rdi;
			loopcount3 = g->normval4;
			do{
				if(g->B_IS_2 == 0)
				{
					if(g->RATIONAL_FFT != 0){
						loopcount4 = 128;
						do{
							xnorm_smallmul_wpn(noexec, noexec, t1,t2,t3,t4,t5);
						}while(--loopcount4);
					}else{
						loopcount4 = g->normval1;
						do{
							loopcount5 = g->cache_line_multiplier;
							do{
								xnorm_smallmul_wpn(exec, noexec, t1,t2,t3,t4,t5);
							}while(--loopcount5);
							rdi += g->normval2;
						}while(--loopcount4);
					}
				}else{
					if(g->RATIONAL_FFT != 0){
						loopcount4 = 128;
						do{
							xnorm_smallmul_wpn(noexec, exec, t1,t2,t3,t4,t5);
						}while(--loopcount4);
					}else{
						loopcount4 = g->normval1;
						do{
							loopcount5 = g->cache_line_multiplier;
							do{
								xnorm_smallmul_wpn(exec, exec, t1,t2,t3,t4,t5);
							}while(--loopcount5);
							rdi += g->normval2;
						}while(--loopcount4);
					}
				} // chunk done
				rsi += 128;
				loopcount3--;
			}while(loopcount3);
			rsi += g->pass2gapsize;
			// Block done
			tmp1 = rsi; rsi = saved_blk_start; saved_blk_start = tmp1;
			tmp1 = rdi; rdi = saved_blk_biglit; saved_blk_biglit = tmp1;
			if(g->B_IS_2 == 0) {
				xnorm_smallmul_wpn_blk(noexec, t1,t2,t3,t4);
			}else{
				xnorm_smallmul_wpn_blk(exec, t1,t2,t3,t4);
			}
			rsi = saved_blk_start;
			rdi = saved_blk_biglit;
			if(g->RATIONAL_FFT == 0) {
				rdi += g->normval3;
			}
			loopcount2--;
		}while(loopcount2);
		if(g->RATIONAL_FFT == 0) {
			rbp += 4*XMM_GMD;
		}
		loopcount1--;
	}while(loopcount1);
	rsi = (uintptr_t)g->DESTARG;
	rbp = (uintptr_t)g->norm_grp_mults;
	rdi = (uintptr_t)g->norm_biglit_array;
	if(g->B_IS_2 == 0) {
		xnorm_smallmul_wpn_sec(noexec, t1,t2,t3,t4);
	}else{
		xnorm_smallmul_wpn_sec(exec, t1,t2,t3,t4);
	}
	rsi = (uintptr_t)g->DESTARG;
	xmm6 = t1;
	xmm7 = t3;
	final_carries_3(g, rsi, xmm6, xmm7);
}
