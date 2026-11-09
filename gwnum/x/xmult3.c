/* Copyright 2001-2012 Mersenne Research, Inc.  All rights reserved */
/* Author:  George Woltman */
/* Email: woltman@alum.mit.edu */
/* */
/* These routine implement some common cleanup code for r4dwpn (r4delay with partial normalization) FFTs */
/* */

#include <stdlib.h>
#include <ctype.h>
#include <stdio.h>

#include "../gwnum.h"
#include "../gwcommon.h"
#include "../unravel.h"
#include "../gwtables.h"
//#define FROM_ASM
#include "../xarch.h"
#include "../xmult.h"
#include "../xbasics.h"
#include "xnormal.h"


/****************************************** */
/* Routine for finishing off a r4dwpn FFT */
/****************************************** */

/* Split the accumulated carries into two carries - a high carry and a */
/* low carry.  Handle both the with and without two-to-phi array cases. */
/* Add these carries back into the FFT data. */

//biglit_incr	EQU	PPTR u64ptr(rsp+first_local);
//grp_incr	EQU	PPTR u64ptr(rsp+first_local+SZPTR);
//loopcount1	EQU	DPTR u64ptr(rsp+first_local+2*SZPTR);
//loopcount2	EQU	DPTR u64ptr(rsp+first_local+2*SZPTR+4);

#ifdef FROM_ASM
void xgw_carries_wpn()
#else
void xgw_carries_wpn(struct gwasm_data *__restrict g)
#endif
{
#ifdef FROM_ASM
	register struct gwasm_data *__restrict g __asm("r11");
	__asm volatile ("" : "=r"(g));  // Let compiler know r11 is set
	//puts("xgw_carries_wpn");
	//__asm volatile ("" : "+r"(g));  // force back into register r11
#endif
	uintptr_t rax, rcx, rdx, rbx, rbp, rsi, rdi;
	vec2f64 xmm0,xmm1,xmm2,xmm3,xmm4,xmm5,xmm6;
	uintptr_t biglit_incr;
	uintptr_t grp_incr, loopcount1, loopcount2;

	rax = g->cache_line_multiplier;/* Cache lines in each pass1 loop */
	rax <<= 1;			/* Compute biglit increment */
	rdx = 4*XMM_GMD;		/* Compute grp increment */
	if(g->RATIONAL_FFT != 0)	/* Don't bump these two pointers */
	{							/* for rational FFTs */
		rax = 0;		/* Zero biglit_incr */
		rdx = 0;		/* Zero grp_incr */
	}
	biglit_incr = rax;	/* Save computed increments */
	grp_incr = rdx;
	rsi = (uintptr_t)g->carries;		/* Addr of the carries */
	if(g->ZERO_PADDED_FFT == 0) {	/* Special case the zero padded FFT case */
		xnorm012_wpn_part1;
		rbp = (uintptr_t)g->data_addr;		/* Addr of the FFT data */
		rdi = (uintptr_t)g->norm_ptr1;		/* Addr of the big/little flags array */
		rdx = (uintptr_t)g->norm_grp_mults;	/* Addr of the group multipliers */
		loopcount2 = g->count3;		/* Load count of grp multipliers */
		do {
			loopcount1 = g->count2;	/* Load wpn count */
			do{
				if(g->B_IS_2 == 0) { /* Is b = 2? */
					xnorm012_wpn(noexec);		/* Split carries for one cache line */
				}else{ /* Yes, do simpler roundings */
					xnorm012_wpn(exec);		/* Split carries for one cache line */
				}
				rsi += 64;			/* Next carries pointer */
				rbp += g->pass1blkdst;	/* Next FFT data pointer */
				rdi += biglit_incr;	/* Next big/little flags pointer */
			}while(--loopcount1);
			rdx += grp_incr;		/* Next group multiplier */
		}while(--loopcount2);
	}else{
		xnorm012_wpn_zpad_part1;
		rbp = (uintptr_t)g->data_addr;		/* Addr of the FFT data */
		rdi = (uintptr_t)g->norm_ptr1;		/* Addr of the big/little flags array */
		rdx = (uintptr_t)g->norm_grp_mults;	/* Addr of the group multipliers */
		loopcount2 = g->count3; /* Load count of grp multipliers */
		do{
			loopcount1 = g->count2;	/* Load wpn count */
			do{
				xnorm012_wpn_zpad;		/* Split carries for one cache line */
				rsi += 64;			/* Next carries pointer */
				rbp += g->pass1blkdst;	/* Next FFT data pointer */
				rdi += biglit_incr;	/* Next big/little flags pointer */
			}while(--loopcount1);
			rdx += grp_incr;		/* Next group multiplier */
		}while(--loopcount2);
	}
}

