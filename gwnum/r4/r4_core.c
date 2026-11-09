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
#include "r4pass1sc.h"

#if 0  // test with pass2 assembly
void xpass2_r4_8_levels_asm(struct gwasm_data *__restrict g) {
	//if(g->ffttype == 99)
	{
		//printf("g->NEGACYCLIC_FFT = %i\n", g->NEGACYCLIC_FFT);
	//	xpass2_r4_8_levels(g);
	//	return;
	}

	register struct gwasm_data *r11_arg __asm("r11") = g;
	__asm volatile (
		/* ".intel_syntax noprefix\n" */
		"call xpass2_r4_8_levels_CORE\n"
		/* ".att_syntax\n" */
		: : /*"r" (g->NORMRTN),*/ "r" (r11_arg)
		: "rax", "rcx", "rdx", "rbx", "rsi", "rdi", "rbp",
		"r8", "r9", "r10", "memory",
		"xmm0", "xmm1", "xmm2", "xmm3", "xmm4", "xmm5", "xmm6", "xmm7",
		"xmm8", "xmm9", "xmm10", "xmm11", "xmm12", "xmm13", "xmm14", "xmm15" );
}
#endif


r4_pass1sc32(xfft_r4_8K_np_8_4_CORE, xpass2_r4_8_levels, 8, 4);
r4_pass1sc64(xfft_r4_16K_np_8_4_CORE, xpass2_r4_8_levels, 8, 4);
r4_pass1sc80(xfft_r4_20K_8_4_CORE, xpass2_r4_8_levels, 8, 4);
r4_pass1sc80(xfft_r4_20K_np_8_4_CORE, xpass2_r4_8_levels, 8, 4);
r4_pass1sc96(xfft_r4_24K_8_4_CORE, xpass2_r4_8_levels, 8, 4);
r4_pass1sc96(xfft_r4_24K_np_8_4_CORE, xpass2_r4_8_levels, 8, 4);
r4_pass1sc32(xfft_r4_24K_np_768_4_CORE, xpass2_r4_768_levels, 768, 4);
r4_pass1sc32(xfft_r4_40K_1280_4_CORE, xpass2_r4_1280_levels, 1280, 4);
r4_pass1sc32(xfft_r4_40K_np_1280_4_CORE, xpass2_r4_1280_levels, 1280, 4);
r4_pass1sc64(xfft_r4_48K_np_768_4_CORE, xpass2_r4_768_levels, 768, 4);
r4_pass1sc224(xfft_r4_56K_8_4_CORE, xpass2_r4_8_levels, 8, 4);
r4_pass1sc96(xfft_r4_72K_np_768_4_CORE, xpass2_r4_768_levels, 768, 4);
r4_pass1sc96(xfft_r4_72K_768_2_CORE, xpass2_r4_768_levels, 768, 2);
r4_pass1sc64(xfft_r4_80K_np_1280_4_CORE, xpass2_r4_1280_levels, 1280, 4);
r4_pass1sc112(xfft_r4_84K_768_4_CORE, xpass2_r4_768_levels, 768, 4);
r4_pass1sc112(xfft_r4_84K_np_768_4_CORE, xpass2_r4_768_levels, 768, 4);
r4_pass1sc96(xfft_r4_144K_1536_4_CORE, xpass2_r4_1536_levels, 1536, 4);

r4_pass1sc32ac(xfft_r4_8K_ac_np_8_4_CORE, xpass2_r4_8_levels, 8, 4);
r4_pass1sc64ac(xfft_r4_16K_ac_np_8_4_CORE, xpass2_r4_8_levels, 8, 4);
r4_pass1sc96ac(xfft_r4_24K_ac_np_8_4_CORE, xpass2_r4_8_levels, 8, 4);
r4_pass1sc32ac(xfft_r4_24K_ac_np_768_4_CORE, xpass2_r4_768_levels, 768, 4);
r4_pass1sc32ac(xfft_r4_40K_ac_1280_4_CORE, xpass2_r4_1280_levels, 1280, 4);
r4_pass1sc32ac(xfft_r4_40K_ac_np_1280_4_CORE, xpass2_r4_1280_levels, 1280, 4);
r4_pass1sc64ac(xfft_r4_48K_ac_np_768_4_CORE, xpass2_r4_768_levels, 768, 4);
r4_pass1sc64ac(xfft_r4_48K_ac_768_4_CORE, xpass2_r4_768_levels, 768, 4);
r4_pass1sc96ac(xfft_r4_72K_ac_768_4_CORE, xpass2_r4_768_levels, 768, 4);
r4_pass1sc96ac(xfft_r4_72K_ac_np_768_4_CORE, xpass2_r4_768_levels, 768, 4);
r4_pass1sc96ac(xfft_r4_72K_ac_768_2_CORE, xpass2_r4_768_levels, 768, 2);
r4_pass1sc64ac(xfft_r4_80K_ac_np_1280_4_CORE, xpass2_r4_1280_levels, 1280, 4);
r4_pass1sc64ac(xfft_r4_80K_ac_1280_2_CORE, xpass2_r4_1280_levels, 1280, 2);
r4_pass1sc96ac(xfft_r4_144K_ac_1536_4_CORE, xpass2_r4_1536_levels, 1536, 4);

#if 0
r4_pass1sc32(xfft_r4_8K_8_4_CORE, xpass2_r4_8_levels, 8, 4);
r4_pass1sc64(xfft_r4_16K_8_4_CORE, xpass2_r4_8_levels, 8, 4);
r4_pass1sc80(xfft_r4_20K_8_4_CORE, xpass2_r4_8_levels, 8, 4);
r4_pass1sc96(xfft_r4_24K_8_4_CORE, xpass2_r4_8_levels, 8, 4);
r4_pass1sc112(xfft_r4_28K_8_4_CORE, xpass2_r4_8_levels, 8, 4);
r4_pass1sc128(xfft_r4_32K_8_4_CORE, xpass2_r4_8_levels, 8, 4);
r4_pass1sc224(xfft_r4_56K_8_4_CORE, xpass2_r4_8_levels, 8, 4);
r4_pass1sc256(xfft_r4_64K_8_4_CORE, xpass2_r4_8_levels, 8, 4);
r4_pass1sc32(xfft_r4_24K_768_4_CORE, xpass2_r4_768_levels, 768, 4);
r4_pass1sc64(xfft_r4_48K_768_4_CORE, xpass2_r4_768_levels, 768, 4);
r4_pass1sc80(xfft_r4_60K_768_4_CORE, xpass2_r4_768_levels, 768, 4);
r4_pass1sc96(xfft_r4_72K_768_4_CORE, xpass2_r4_768_levels, 768, 4);
r4_pass1sc112(xfft_r4_84K_768_4_CORE, xpass2_r4_768_levels, 768, 4);
r4_pass1sc128(xfft_r4_96K_768_4_CORE, xpass2_r4_768_levels, 768, 4);
r4_pass1sc224(xfft_r4_168K_768_4_CORE, xpass2_r4_768_levels, 768, 4);
r4_pass1sc256(xfft_r4_192K_768_4_CORE, xpass2_r4_768_levels, 768, 4);
r4_pass1sc32(xfft_r4_32K_10_4_CORE, xpass2_r4_10_levels, 10, 4);
r4_pass1sc64(xfft_r4_64K_10_4_CORE, xpass2_r4_10_levels, 10, 4);
r4_pass1sc80(xfft_r4_80K_10_4_CORE, xpass2_r4_10_levels, 10, 4);
r4_pass1sc96(xfft_r4_96K_10_4_CORE, xpass2_r4_10_levels, 10, 4);
r4_pass1sc112(xfft_r4_112K_10_4_CORE, xpass2_r4_10_levels, 10, 4);
r4_pass1sc128(xfft_r4_128K_10_4_CORE, xpass2_r4_10_levels, 10, 4);
r4_pass1sc224(xfft_r4_224K_10_4_CORE, xpass2_r4_10_levels, 10, 4);
r4_pass1sc256(xfft_r4_256K_10_4_CORE, xpass2_r4_10_levels, 10, 4);
r4_pass1sc32(xfft_r4_40K_1280_4_CORE, xpass2_r4_1280_levels, 1280, 4);
r4_pass1sc64(xfft_r4_80K_1280_4_CORE, xpass2_r4_1280_levels, 1280, 4);
r4_pass1sc80(xfft_r4_100K_1280_4_CORE, xpass2_r4_1280_levels, 1280, 4);
r4_pass1sc96(xfft_r4_120K_1280_4_CORE, xpass2_r4_1280_levels, 1280, 4);
r4_pass1sc112(xfft_r4_140K_1280_4_CORE, xpass2_r4_1280_levels, 1280, 4);
r4_pass1sc128(xfft_r4_160K_1280_4_CORE, xpass2_r4_1280_levels, 1280, 4);
r4_pass1sc224(xfft_r4_280K_1280_4_CORE, xpass2_r4_1280_levels, 1280, 4);
r4_pass1sc256(xfft_r4_320K_1280_4_CORE, xpass2_r4_1280_levels, 1280, 4);
r4_pass1sc32(xfft_r4_48K_1536_4_CORE, xpass2_r4_1536_levels, 1536, 4);
r4_pass1sc64(xfft_r4_96K_1536_4_CORE, xpass2_r4_1536_levels, 1536, 4);
r4_pass1sc80(xfft_r4_120K_1536_4_CORE, xpass2_r4_1536_levels, 1536, 4);
r4_pass1sc96(xfft_r4_144K_1536_4_CORE, xpass2_r4_1536_levels, 1536, 4);
r4_pass1sc112(xfft_r4_168K_1536_4_CORE, xpass2_r4_1536_levels, 1536, 4);
r4_pass1sc128(xfft_r4_192K_1536_4_CORE, xpass2_r4_1536_levels, 1536, 4);
r4_pass1sc224(xfft_r4_336K_1536_4_CORE, xpass2_r4_1536_levels, 1536, 4);
r4_pass1sc256(xfft_r4_384K_1536_4_CORE, xpass2_r4_1536_levels, 1536, 4);

r4_pass1sc32ac(xfft_r4_8K_ac_8_4_CORE, xpass2_r4_8_levels, 8, 4);
r4_pass1sc64ac(xfft_r4_16K_ac_8_4_CORE, xpass2_r4_8_levels, 8, 4);
r4_pass1sc96ac(xfft_r4_24K_ac_8_4_CORE, xpass2_r4_8_levels, 8, 4);
r4_pass1sc128ac(xfft_r4_32K_ac_8_4_CORE, xpass2_r4_8_levels, 8, 4);
r4_pass1sc256ac(xfft_r4_64K_ac_8_4_CORE, xpass2_r4_8_levels, 8, 4);
r4_pass1sc32ac(xfft_r4_24K_ac_768_4_CORE, xpass2_r4_768_levels, 768, 4);
r4_pass1sc64ac(xfft_r4_48K_ac_768_4_CORE, xpass2_r4_768_levels, 768, 4);
r4_pass1sc96ac(xfft_r4_72K_ac_768_4_CORE, xpass2_r4_768_levels, 768, 4);
r4_pass1sc128ac(xfft_r4_96K_ac_768_4_CORE, xpass2_r4_768_levels, 768, 4);
r4_pass1sc256ac(xfft_r4_192K_ac_768_4_CORE, xpass2_r4_768_levels, 768, 4);
r4_pass1sc32ac(xfft_r4_32K_ac_10_4_CORE, xpass2_r4_10_levels, 10, 4);
r4_pass1sc64ac(xfft_r4_64K_ac_10_4_CORE, xpass2_r4_10_levels, 10, 4);
r4_pass1sc96ac(xfft_r4_96K_ac_10_4_CORE, xpass2_r4_10_levels, 10, 4);
r4_pass1sc128ac(xfft_r4_128K_ac_10_4_CORE, xpass2_r4_10_levels, 10, 4);
r4_pass1sc256ac(xfft_r4_256K_ac_10_4_CORE, xpass2_r4_10_levels, 10, 4);
r4_pass1sc32ac(xfft_r4_40K_ac_1280_4_CORE, xpass2_r4_1280_levels, 1280, 4);
r4_pass1sc64ac(xfft_r4_80K_ac_1280_4_CORE, xpass2_r4_1280_levels, 1280, 4);
r4_pass1sc96ac(xfft_r4_120K_ac_1280_4_CORE, xpass2_r4_1280_levels, 1280, 4);
r4_pass1sc128ac(xfft_r4_160K_ac_1280_4_CORE, xpass2_r4_1280_levels, 1280, 4);
r4_pass1sc256ac(xfft_r4_320K_ac_1280_4_CORE, xpass2_r4_1280_levels, 1280, 4);
r4_pass1sc32ac(xfft_r4_48K_ac_1536_4_CORE, xpass2_r4_1536_levels, 1536, 4);
r4_pass1sc64ac(xfft_r4_96K_ac_1536_4_CORE, xpass2_r4_1536_levels, 1536, 4);
r4_pass1sc96ac(xfft_r4_144K_ac_1536_4_CORE, xpass2_r4_1536_levels, 1536, 4);
r4_pass1sc128ac(xfft_r4_192K_ac_1536_4_CORE, xpass2_r4_1536_levels, 1536, 4);
r4_pass1sc256ac(xfft_r4_384K_ac_1536_4_CORE, xpass2_r4_1536_levels, 1536, 4);

#endif
