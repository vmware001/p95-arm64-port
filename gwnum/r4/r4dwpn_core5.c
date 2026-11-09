#include <stdlib.h>
#include <ctype.h>

#include "../gwnum.h"
#include "../gwcommon.h"
#include "../unravel.h"
#include "../gwtables.h"
#include "../xarch.h"
#include "../xbasics.h"
#include "../xmult.h"
#include "../gwnum_func.h"

#include "r3.h"
#include "r4.h"
#include "r5.h"
#include "r7.h"
#include "r8.h"
#include "r4pass2.h"
#include "r4dwpnpass1sc.h"

#if 0  // test with pass2 assembly
void xpass2_r4_14_levels_asm(struct gwasm_data *__restrict g) {
	if(g->ffttype <= 9)
	{
		printf("ffttype %i, NEGACYCLIC_FFT = %i\n", g->ffttype, g->NEGACYCLIC_FFT);
		xpass2_r4_14_levels(g);
		return;
	}

	register struct gwasm_data *r11_arg __asm("r11") = g;
	__asm volatile (
		/* ".intel_syntax noprefix\n" */
		"call xpass2_r4dwpn_14_levels_CORE\n"
		/* ".att_syntax\n" */
		: : /*"r" (g->NORMRTN),*/ "r" (r11_arg)
		: "rax", "rcx", "rdx", "rbx", "rsi", "rdi", "rbp",
		"r8", "r9", "r10", "memory",
		"xmm0", "xmm1", "xmm2", "xmm3", "xmm4", "xmm5", "xmm6", "xmm7",
		"xmm8", "xmm9", "xmm10", "xmm11", "xmm12", "xmm13", "xmm14", "xmm15" );
}
#endif

r4dwpn_pass1sc896(xfft_r4dwpn_8960K_10240_4_CORE, xpass2_r4_10240_levels, 10240, 4);
r4dwpn_pass1sc1024(xfft_r4dwpn_9M_9216_4_CORE, xpass2_r4_9216_levels, 9216, 4);
r4dwpn_pass1sc768(xfft_r4dwpn_9600K_12800_4_CORE, xpass2_r4_12800_levels, 12800, 4);
r4dwpn_pass1sc1024(xfft_r4dwpn_10M_10240_4_CORE, xpass2_r4_10240_levels, 10240, 4);
r4dwpn_pass1sc3584(xfft_r4dwpn_10752K_3072_4_CORE, xpass2_r4_3072_levels, 3072, 4);
r4dwpn_pass1sc896(xfft_r4dwpn_11200K_12800_4_CORE, xpass2_r4_12800_levels, 12800, 4);
r4dwpn_pass1sc2560(xfft_r4dwpn_11520K_4608_4_CORE, xpass2_r4_4608_levels, 4608, 4);
r4dwpn_pass1sc768(xfft_r4dwpn_12M_14_4_CORE, xpass2_r4_14_levels, 14, 4);
r4dwpn_pass1sc1024(xfft_r4dwpn_12800K_12800_4_CORE, xpass2_r4_12800_levels, 12800, 4);
r4dwpn_pass1sc3584(xfft_r4dwpn_13440K_3840_4_CORE, xpass2_r4_3840_levels, 3840, 4);
r4dwpn_pass1sc3072(xfft_r4dwpn_13824K_4608_4_CORE, xpass2_r4_4608_levels, 4608, 4);
r4dwpn_pass1sc3584(xfft_r4dwpn_14M_12_4_CORE, xpass2_r4_12_levels, 12, 4);
r4dwpn_pass1sc1024(xfft_r4dwpn_15M_15360_4_CORE, xpass2_r4_15360_levels, 15360, 4);
r4dwpn_pass1sc640(xfft_r4dwpn_16000K_25600_4_CORE, xpass2_r4_25600_levels, 25600, 4);
r4dwpn_pass1sc4096(xfft_r4dwpn_16M_12_4_CORE, xpass2_r4_12_levels, 12, 4);
r4dwpn_pass1sc3584(xfft_r4dwpn_17920K_5120_4_CORE, xpass2_r4_5120_levels, 5120, 4);
r4dwpn_pass1sc2048(xfft_r4dwpn_18M_9216_4_CORE, xpass2_r4_9216_levels, 9216, 4);
r4dwpn_pass1sc4096(xfft_r4dwpn_18M_4608_4_CORE, xpass2_r4_4608_levels, 4608, 4);
r4dwpn_pass1sc768(xfft_r4dwpn_19200K_25600_4_CORE, xpass2_r4_25600_levels, 25600, 4);
r4dwpn_pass1sc3072(xfft_r4dwpn_19200K_6400_4_CORE, xpass2_r4_6400_levels, 6400, 4);
r4dwpn_pass1sc4096(xfft_r4dwpn_20M_5120_4_CORE, xpass2_r4_5120_levels, 5120, 4);
r4dwpn_pass1sc3584(xfft_r4dwpn_21M_6144_4_CORE, xpass2_r4_6144_levels, 6144, 4);
r4dwpn_pass1sc896(xfft_r4dwpn_22400K_25600_4_CORE, xpass2_r4_25600_levels, 25600, 4);
r4dwpn_pass1sc1536(xfft_r4dwpn_23040K_15360_4_CORE, xpass2_r4_15360_levels, 15360, 4);
r4dwpn_pass1sc4096(xfft_r4dwpn_24M_6144_4_CORE, xpass2_r4_6144_levels, 6144, 4);
r4dwpn_pass1sc1024(xfft_r4dwpn_25M_25600_4_CORE, xpass2_r4_25600_levels, 25600, 4);
r4dwpn_pass1sc3584(xfft_r4dwpn_26880K_7680_4_CORE, xpass2_r4_7680_levels, 7680, 4);
r4dwpn_pass1sc3072(xfft_r4dwpn_27M_9216_4_CORE, xpass2_r4_9216_levels, 9216, 4);
r4dwpn_pass1sc3584(xfft_r4dwpn_28M_13_4_CORE, xpass2_r4_13_levels, 13, 4);
r4dwpn_pass1sc4096(xfft_r4dwpn_30M_7680_4_CORE, xpass2_r4_7680_levels, 7680, 4);
r4dwpn_pass1sc1280(xfft_r4dwpn_32000K_25600_4_CORE, xpass2_r4_25600_levels, 25600, 4);
r4dwpn_pass1sc4096(xfft_r4dwpn_32M_13_4_CORE, xpass2_r4_13_levels, 13, 4);
r4dwpn_pass1sc3072(xfft_r4dwpn_75M_25600_4_CORE, xpass2_r4_25600_levels, 25600, 4);
