#include <stdlib.h>
#include <ctype.h>

#include "../gwnum.h"
#include "../gwcommon.h"
#include "../unravel.h"
#include "../gwtables.h"
#include "../xarch.h"
#include "../xbasics.h"
#include "xnormal.h"
#include <stdio.h>

typedef double __attribute__ ((vector_size (64), aligned(16))) vec8f64b;
typedef double __attribute__ ((vector_size (32), aligned(16))) vec4f64b;


void gwxaddq1(struct gwasm_data *__restrict g)
{
	vec8f64b *c = g->SRCARG;
	vec8f64b *d = g->SRC2ARG;
	vec8f64b *si = g->DESTARG;
	int a = g->addcount1;
	do{
		si[0] = d[0] + c[0];
		c++; d++; si++;
		a--;
	}while(a != 0);
}
//DEST = SRC2 - SRC
void gwxsubq1(struct gwasm_data *__restrict g)
{
	vec8f64b *c = g->SRCARG;
	vec8f64b *d = g->SRC2ARG;
	vec8f64b *si = g->DESTARG;
	int a = g->addcount1;
	do{
		si[0] = d[0] - c[0];
		c++; d++; si++;
		a--;
	}while(a != 0);
}

// DEST = SRC + SRC2, DEST2 = SRC - SRC2
void gwxaddsubq1(struct gwasm_data *__restrict g)
{
	vec4f64b *c = g->SRCARG;
	vec4f64b *d = g->SRC2ARG;
	vec4f64b *si = g->DESTARG;
	vec4f64b *bp = g->DEST2ARG;
	int a = g->addcount1;
	do{
		vec4f64b x = c[0];
		vec4f64b y = d[0];
		si[0] = x + y;
		bp[0] = x - y;
		vec4f64b x2 = c[1];
		vec4f64b y2 = d[1];
		si[1] = x2 + y2;
		bp[1] = x2 - y2;
		c+=2; d+=2; si+=2; bp+=2;
		a--;
	}while(a != 0);
}

typedef unsigned long long __attribute__ ((vector_size (64), aligned(16))) vec8u64b;
void gwxcopy4kb(struct gwasm_data *__restrict g)
{
	vec8u64b *s1 = g->SRCARG;
	vec8u64b *s2 = g->SRC2ARG;
	vec8u64b *di = g->DESTARG;
	int a = 4096/64; /* Count of 64 byte chunks in 4KB */
	do{
		di[0] = s1[0] & s2[0];
		s1++; s2++; di++;
		a--;
	}while(a != 0);
}





#if 0
#define final_carries_1 { \
	vec2f64 *bp = g->norm_col_mults; \
	vec2f64 *di = g->norm_biglit_array; \
	if(g->B_IS_2 == 0) { \
		xnorm_smallmul_1d_cleanup(noexec); \
	}else{ \
		xnorm_smallmul_1d_cleanup(exec); \
	} \
}
#else
static void final_carries_1(struct gwasm_data *__restrict g, uintptr_t rsi, vec2f64 xmm2, vec2f64 xmm3) {
	uintptr_t rdi,rbp;
	rbp = (uintptr_t)g->norm_col_mults;
	rdi = (uintptr_t)g->norm_biglit_array;	// Addr of the big/little flags array
#if 0
	if(g->ZERO_PADDED_FFT == 0) {	/* Zero-padded FFT? */
		xnorm_top_carry_1d;	/* No, do a very standard carry */
		if(g->B_IS_2 == 0) {
			xnorm012_1d(noexec, noexec);
		}else{
			xnorm012_1d(noexec, exec);
		}
	}else{			/* Yes, do special zpad carry */
		if(g->B_IS_2 == 0) {
			xnorm_smallmul_1d_zpad_cleanup(noexec);
		}else{
			xnorm_smallmul_1d_zpad_cleanup(exec);
		}
	}
#else
	if(g->B_IS_2 == 0) {
		xnorm_smallmul_1d_cleanup(noexec);
	}else{
		xnorm_smallmul_1d_cleanup(exec);
	};
#endif
}
#endif

#if 0
void gwxadd1(struct gwasm_data *__restrict g)
{
	vec2f64 *c = g->SRCARG;
	vec2f64 *d = g->SRC2ARG;
	vec2f64 *si = g->DESTARG;
	vec2f64 xmm0,xmm1,xmm2,xmm3,xmm4,xmm5;
	xmm2 = *(vec2f64*)&g->u.xmm.XMM_BIGVAL;
	xmm3 = xmm2;
	int loopcount1 = g->normcount1;
	vec2f64 *bp = g->norm_col_mults;
	unsigned char *di_u8 = g->norm_biglit_array;
	vec2f64 *saved_dest_ptr = si;
	vec2f64 *saved_col_ptr = bp;
//nadd0:
	do{
		int loopcount2 = loopcount1 & 0x7FF;
		if(g->B_IS_2 == 0)
		{
			if(g->RATIONAL_FFT != 0){
				do{
					xnorm_op_1d(addpd, noexec, noexec);
				}while(--loopcount2);
			}else{
				do{
					xnorm_op_1d(addpd, exec, noexec);
				}while(--loopcount2);
			}
		}else{
			if(g->RATIONAL_FFT != 0){
				do{
					xnorm_op_1d(addpd, noexec, exec);
				}while(--loopcount2);
			}else{
				do{
					xnorm_op_1d(addpd, exec, exec);
				}while(--loopcount2);
			}
		}
		vec2f64 *rbx = saved_col_ptr;
		vec2f64 *rax = saved_dest_ptr;
		xnorm_op_1d_mid_cleanup;
		loopcount1 >>= 11;
	}while(loopcount1 != 0);
	vec2f64 *rbx = saved_col_ptr;
	vec2f64 *rax = saved_dest_ptr;
	vec2f64 *rsi = g->DESTARG;
	final_carries_1;
}
#endif


void gwxadd1(struct gwasm_data *__restrict g) {
	vec2f64 xmm0,xmm1,xmm2,xmm3,xmm4,xmm5;
	uintptr_t rax, rcx, rdx, rbx, rbp, rdi, rsi;
	uintptr_t saved_dest_ptr, saved_col_ptr;
	unsigned int loopcount1, loopcount2;

	rcx = (uintptr_t)g->SRCARG;
	rdx = (uintptr_t)g->SRC2ARG;
	rsi = (uintptr_t)g->DESTARG;
	xmm2 = XMM_BIGVAL2;	// Start process with no carry
	xmm3 = xmm2;
	loopcount1 = g->normcount1;
	rbp = (uintptr_t)g->norm_col_mults;	// Address of the multipliers
	rdi = (uintptr_t)g->norm_biglit_array;	// Addr of the big/little flags array
	do{
		loopcount2 = loopcount1 & 0x7FF;		// Grab 11 bits of the counter
		saved_dest_ptr = rsi;	// remember rsi for xnorm_op_1d_mid
		saved_col_ptr = rbp;	// remember rbp for xnorm_op_1d_mid
		rax = 0;		// Clear big/lit flags
		rbx = 0;
		if(g->B_IS_2 == 0)
		{
			if(g->RATIONAL_FFT != 0){
				do{
					xnorm_op_1d(addpd, noexec, noexec);
				}while(--loopcount2);
			}else{
				do{
					xnorm_op_1d(addpd, exec, noexec);
				}while(--loopcount2);
			}
		}else{
			if(g->RATIONAL_FFT != 0){
				do{
					xnorm_op_1d(addpd, noexec, exec);
				}while(--loopcount2);
			}else{
				do{
					xnorm_op_1d(addpd, exec, exec);
				}while(--loopcount2);
			}
		}
		rbx = saved_col_ptr;	// Restore multipliers pointer
		rax = saved_dest_ptr;	// Restore dest pointer
		xnorm_op_1d_mid_cleanup;		// Rotate carries and add in carries
		loopcount1 = loopcount1 >> 11;			// Get next loop amount
	}while(loopcount1 != 0);

	rsi = (uintptr_t)(g->DESTARG);
	final_carries_1(g, rsi, xmm2, xmm3);
}

// gwxsub1 is copied from gwxadd1, then change addpd to subpd
void gwxsub1(struct gwasm_data *__restrict g) {
	vec2f64 xmm0,xmm1,xmm2,xmm3,xmm4,xmm5;
	uintptr_t rax, rcx, rdx, rbx, rbp, rdi, rsi;
	uintptr_t saved_dest_ptr, saved_col_ptr;
	unsigned int loopcount1, loopcount2;

	rcx = (uintptr_t)g->SRCARG;
	rdx = (uintptr_t)g->SRC2ARG;
	rsi = (uintptr_t)g->DESTARG;
	xmm2 = XMM_BIGVAL2;	// Start process with no carry
	xmm3 = xmm2;
	loopcount1 = g->normcount1;
	rbp = (uintptr_t)g->norm_col_mults;	// Address of the multipliers
	rdi = (uintptr_t)g->norm_biglit_array;	// Addr of the big/little flags array
	do{
		loopcount2 = loopcount1 & 0x7FF;		// Grab 11 bits of the counter
		saved_dest_ptr = rsi;	// remember rsi for xnorm_op_1d_mid
		saved_col_ptr = rbp;	// remember rbp for xnorm_op_1d_mid
		rax = 0;		// Clear big/lit flags
		rbx = 0;
		if(g->B_IS_2 == 0)
		{
			if(g->RATIONAL_FFT != 0){
				do{
					xnorm_op_1d(subpd, noexec, noexec);
				}while(--loopcount2);
			}else{
				do{
					xnorm_op_1d(subpd, exec, noexec);
				}while(--loopcount2);
			}
		}else{
			if(g->RATIONAL_FFT != 0){
				do{
					xnorm_op_1d(subpd, noexec, exec);
				}while(--loopcount2);
			}else{
				do{
					xnorm_op_1d(subpd, exec, exec);
				}while(--loopcount2);
			}
		}
		rbx = saved_col_ptr;	// Restore multipliers pointer
		rax = saved_dest_ptr;	// Restore dest pointer
		xnorm_op_1d_mid_cleanup;		// Rotate carries and add in carries
		loopcount1 = loopcount1 >> 11;			// Get next loop amount
	}while(loopcount1 != 0);

	rsi = (uintptr_t)(g->DESTARG);
	final_carries_1(g, rsi, xmm2, xmm3);
}

void gwxaddsub1(struct gwasm_data *__restrict g) {
	vec2f64 xmm0,xmm1,xmm2,xmm3,xmm4,xmm5,xmm6,xmm7;
	uintptr_t rax, rcx, rdx, rbx, rbp, rdi, rsi;
	uintptr_t saved_dest1_ptr, saved_dest2_ptr, saved_col_ptr;
	unsigned int loopcount1, loopcount2;

	rcx = (uintptr_t)g->SRCARG;
	rdx = (uintptr_t)g->SRC2ARG;
	rsi = (uintptr_t)g->DESTARG;
	rbp = (uintptr_t)g->DEST2ARG;
	xmm2 = XMM_BIGVAL2;	// Start process with no carry
	xmm3 = xmm2;
	xmm6 = xmm2;
	xmm7 = xmm2;
	loopcount1 = g->normcount1;		// Load loop counter
	rbx = (uintptr_t)g->norm_col_mults;	// Address of the multipliers
	rdi = (uintptr_t)g->norm_biglit_array;	// Addr of the big/little flags array
	do{
		loopcount2 = loopcount1 & 0x7FF;		// Grab 11 bits of the counter
		mov(	saved_dest1_ptr, rsi);
		mov(	saved_dest2_ptr, rbp);
		mov(	saved_col_ptr, rbx);
		rax = 0;		// Clear big/lit flags
		if(g->B_IS_2 == 0)
		{
			if(g->RATIONAL_FFT != 0){
				do{
					xnorm_addsub_1d(noexec, noexec);
				}while(--loopcount2);
			}else{
				do{
					xnorm_addsub_1d(exec, noexec);
				}while(--loopcount2);
			}
		}else{
			if(g->RATIONAL_FFT != 0){
				do{
					xnorm_addsub_1d(noexec, exec);
				}while(--loopcount2);
			}else{
				do{
					xnorm_addsub_1d(exec, exec);
				}while(--loopcount2);
			}
		}
		{uintptr_t t1 = rbx; rbx = saved_col_ptr; saved_col_ptr = t1;}; // xchg
		xnorm_addsub_1d_mid_cleanup(saved_dest1_ptr, saved_dest2_ptr);
		//rax = saved_dest1_ptr;
		//xnorm_op_1d_mid_cleanup(xmm2, xmm3);
		//rax = saved_dest2_ptr;
		//xnorm_op_1d_mid_cleanup(xmm6, xmm7);
		rbx = saved_col_ptr;
		loopcount1 = loopcount1 >> 11;			// Get next loop amount
	}while(loopcount1 != 0);

	rsi = (uintptr_t)(g->DESTARG);
	final_carries_1(g, rsi, xmm2, xmm3);
	rsi = (uintptr_t)(g->DEST2ARG);
	final_carries_1(g, rsi, xmm6, xmm7);
}



// gwxsub1 is copied from gwxadd1, then change addpd to subpd
void gwxmuls1(struct gwasm_data *__restrict g) {
	vec2f64 xmm0,xmm1,xmm2,xmm3	,xmm4,xmm5,xmm7;
	uintptr_t rbp, rdi, rsi, tmp1;
	uintptr_t saved_dest_ptr, saved_col_ptr, saved_biglit_ptr;
	unsigned int loopcount1, loopcount2;

	rsi = (uintptr_t)g->DESTARG;
	xmm7[1] = xmm7[0] = g->DBLARG;
	if(g->RATIONAL_FFT == 0)
		xmm7 *= g->u.xmm.XMM_NORM012_FF;
	xmm2 = XMM_BIGVAL2;    // Start process with no carry
	xmm3 = xmm2;

	loopcount1 = g->normcount1;
	rbp = (uintptr_t)g->norm_col_mults;	// Address of the multipliers
	rdi = (uintptr_t)g->norm_biglit_array;	// Addr of the big/little flags array
	saved_dest_ptr = rsi; // remember for xnorm_smallmul_1d_mid
	saved_col_ptr = rbp;
	saved_biglit_ptr = rdi;
	if(g->B_IS_2 == 0){
		do{
			loopcount2 = loopcount1 & 0x7FF;		// Grab 11 bits of the counter
			//saved_dest_ptr = rsi;	// remember rsi for xnorm_op_1d_mid
			//saved_col_ptr = rbp;	// remember rbp for xnorm_op_1d_mid
			if(g->RATIONAL_FFT != 0){
				do{
					xnorm_smallmul_1d(noexec, noexec);
				}while(--loopcount2);
			}else{
				do{
					xnorm_smallmul_1d(exec, noexec);
				}while(--loopcount2);
			}
			tmp1 = rsi; rsi = saved_dest_ptr; saved_dest_ptr = tmp1;
			tmp1 = rdi; rdi = saved_biglit_ptr; saved_biglit_ptr = tmp1;
			tmp1 = rbp; rbp = saved_col_ptr;; saved_col_ptr = tmp1;
			xnorm_smallmul_1d_mid_cleanup(noexec);		// Rotate carries and add in carries
			rsi = saved_dest_ptr;
			rdi = saved_biglit_ptr;
			rbp = saved_col_ptr;
			loopcount1 = loopcount1 >> 11;			// Get next loop amount
		}while(loopcount1 != 0);
	}else{
		do{
			loopcount2 = loopcount1 & 0x7FF;		// Grab 11 bits of the counter
			//saved_dest_ptr = rsi;	// remember rsi for xnorm_op_1d_mid
			//saved_col_ptr = rbp;	// remember rbp for xnorm_op_1d_mid
			if(g->RATIONAL_FFT != 0){
				do{
					xnorm_smallmul_1d(noexec, exec);
				}while(--loopcount2);
			}else{
				do{
					xnorm_smallmul_1d(exec, exec);
				}while(--loopcount2);
			}
			tmp1 = rsi; rsi = saved_dest_ptr; saved_dest_ptr = tmp1;
			tmp1 = rdi; rdi = saved_biglit_ptr; saved_biglit_ptr = tmp1;
			tmp1 = rbp; rbp = saved_col_ptr;; saved_col_ptr = tmp1;
			xnorm_smallmul_1d_mid_cleanup(exec);		// Rotate carries and add in carries
			rsi = saved_dest_ptr;
			rdi = saved_biglit_ptr;
			rbp = saved_col_ptr;
			loopcount1 = loopcount1 >> 11;			// Get next loop amount
		}while(loopcount1 != 0);
	}

	rsi = (uintptr_t)(g->DESTARG);
	final_carries_1(g, rsi, xmm2, xmm3);
}
