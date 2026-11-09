#include <stdlib.h>
#include <ctype.h>
#include <stdio.h>

#include "../gwnum.h"
#include "../gwtables.h"
#include "../gwcommon.h"
#include "../unravel.h"
//#define FROM_ASM
#include "../xarch.h"
#include "../xbasics.h"
#include "../xmult.h"
#include "xnormal.h"

#if 0
#include <xmmintrin.h>
#endif

typedef void (*func_gwasm_data)(struct gwasm_data *__restrict g);
//void xpass2_r4_8_levels_t(struct gwasm_data *__restrict g);
//void xpass2_hg_8_levels_t_P4(struct gwasm_data *__restrict g);
//void xpass2_r4_8_levels(struct gwasm_data *__restrict g);

#if 0
void pass1_aux_entry_point(struct gwasm_data *__restrict g) {
	//__sync_synchronize(); // ARM64 may need this to avoid multithread errors.
#if 0
	unsigned int mxcsr = _mm_getcsr();
	mxcsr &= ~(_MM_MASK_DIV_ZERO | _MM_MASK_OVERFLOW | _MM_MASK_INVALID);
	_mm_setcsr(mxcsr);
#endif

	//static void *p; if(p != g->thread_work_routine) {
	//	p = g->thread_work_routine;
	//	printf("pass1_aux_entry_point %p\n", p);
	//}
	((func_gwasm_data)g->thread_work_routine)(g);
}
#endif
#if 0
void pass2_aux_entry_point(struct gwasm_data *__restrict g) {
	//__sync_synchronize(); // ARM64 may need this to avoid multithread errors.
#if 0
	unsigned int mxcsr = _mm_getcsr();
	mxcsr &= ~(_MM_MASK_DIV_ZERO | _MM_MASK_OVERFLOW | _MM_MASK_INVALID);
	_mm_setcsr(mxcsr);
#endif

	//static void *p; if(p != g->thread_work_routine) {
	//	p = g->thread_work_routine;
	//	printf("pass2_aux_entry_point %p\n", p);
	//}

	//xpass2_r4_8_levels_t(g);
	//xpass2_hg_8_levels_t_P4(g);
	((func_gwasm_data)g->thread_work_routine)(g);
}
#endif

#if 1
// used in xmult.mac line 320, then in r4delay which is unused
// used in hgpass1
#ifdef FROM_ASM
__attribute__((preserve_all)) void xgw_carries()
#else
void xgw_carries(struct gwasm_data *__restrict g)
#endif
{
#ifdef FROM_ASM
	register struct gwasm_data *__restrict g __asm("r11");
	__asm volatile ("" : "=r"(g)); // Let compiler know r11 is set.
	//puts("xgw_carries");
#endif

	uintptr_t rax,rcx,rdx,rbx,rsi,rdi,rbp;
	vec2f64 xmm0,xmm1,xmm2,xmm3,xmm4,xmm5,xmm6,xmm7;
	uintptr_t loopcount1;
	rsi = (uintptr_t)g->carries;		/* Addr of the carries */
	if(g->ZERO_PADDED_FFT == 0) { /* Special case the zero padded FFT case */
		xnorm012_2d_part1;
		rbp = (uintptr_t)g->data_addr;		/* Addr of the FFT data */
		rdi = (uintptr_t)g->norm_ptr1;		/* Addr of the big/little flags array */
		rdx = (uintptr_t)g->norm_grp_mults;	/* Addr of the group multipliers */
		rbx = (uintptr_t)g->norm_ptr2;		/* Addr of the column multipliers */
		rax = g->addcount1;		/* Load count of cache lines in carries array */
		loopcount1 = rax;		/* Save for later */
		if(g->B_IS_2 == 0)		/* Is b = 2? */
		{					/* Yes, do simpler roundings */
			do{
				xnorm012_2d(noexec);		/* Split carries for one cache line */
				rax = g->cache_line_multiplier;/* Cache lines in each pass1 loop */
				rsi += 64;			/* Next carries pointer */
				rbp += g->pass1blkdst;	/* Next FFT data pointer */
				if(g->RATIONAL_FFT == 0) {
					rdx += 128;		/* Next group multiplier */
					rdi = rdi+rax*4;	/* Next big/little flags pointer */
				}
				loopcount1--;		/* Test loop counter */
			}while(loopcount1 != 0);			/* Next carry row */			/* Done */
		}else{
			do{
				xnorm012_2d(exec);		/* Split carries for one cache line */
				rax = g->cache_line_multiplier;/* Cache lines in each pass1 loop */
				rsi += 64;			/* Next carries pointer */
				rbp += g->pass1blkdst;	/* Next FFT data pointer */
				if(g->RATIONAL_FFT == 0) {
					rdx += 128;		/* Next group multiplier */
					rdi = rdi+rax*4;	/* Next big/little flags pointer */
				}
				loopcount1--;		/* Test loop counter */
			}while(loopcount1 != 0);			/* Next carry row */
		}
	}else{			/* Done */
		xnorm012_2d_zpad_part1;
		rbp = (uintptr_t)g->data_addr;		/* Addr of the FFT data */
		rdi = (uintptr_t)g->norm_ptr1;		/* Addr of the big/little flags array */
		rdx = (uintptr_t)g->norm_grp_mults;	/* Addr of the group multipliers */
		rbx = (uintptr_t)g->norm_ptr2;		/* Addr of the column multipliers */
		rax = g->addcount1;		/* Load count of cache lines in carries array */
		loopcount1 = rax;		/* Save for later */
		do{
			xnorm012_2d_zpad;		/* Split carries for one cache line */
			rax = g->cache_line_multiplier;/* Cache lines in each pass1 loop */
			rsi += 64;			/* Next carries pointer */
			rbp += g->pass1blkdst;	/* Next FFT data pointer */
			if(g->RATIONAL_FFT == 0) {
				rdx += 128;		/* Next group multiplier */
				rdi = rdi+rax*4;	/* Next big/little flags pointer */
			}
			loopcount1--;
		}while(loopcount1 != 0);			/* Next carry row in section */
	}
#ifdef FROM_ASM
	__asm volatile ("" : "+r"(g));  // force back into register r11
#endif
}
#endif
