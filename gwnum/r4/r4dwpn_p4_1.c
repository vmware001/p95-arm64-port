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
void xpass2_r4_11_levels_asm(struct gwasm_data *__restrict g) {
    if(g->ffttype == 4)
    {
        //printf("g->NEGACYCLIC_FFT = %i\n", g->NEGACYCLIC_FFT);
        xpass2_r4_11_levels(g);
        return;
    }

    register struct gwasm_data *r11_arg __asm("r11") = g;
    __asm volatile (
        /* ".intel_syntax noprefix\n" */
        "call xpass2_r4dwpn_11_levels_CORE\n"
        /* ".att_syntax\n" */
        : : /*"r" (g->NORMRTN),*/ "r" (r11_arg)
        : "rax", "rcx", "rdx", "rbx", "rsi", "rdi", "rbp",
        "r8", "r9", "r10", "memory",
        "xmm0", "xmm1", "xmm2", "xmm3", "xmm4", "xmm5", "xmm6", "xmm7",
        "xmm8", "xmm9", "xmm10", "xmm11", "xmm12", "xmm13", "xmm14", "xmm15" );
}
#endif


#if 1
r4dwpn_pass1sc128(xfft_r4dwpn_32K_8_4_P4, xpass2_r4_8_levels, 8, 4);
r4dwpn_pass1sc256(xfft_r4dwpn_64K_8_4_P4, xpass2_r4_8_levels, 8, 4);
r4dwpn_pass1sc320(xfft_r4dwpn_80K_8_4_P4, xpass2_r4_8_levels, 8, 4);
r4dwpn_pass1sc448(xfft_r4dwpn_112K_8_4_P4, xpass2_r4_8_levels, 8, 4);
r4dwpn_pass1sc640(xfft_r4dwpn_160K_8_4_P4, xpass2_r4_8_levels, 8, 4);
r4dwpn_pass1sc896(xfft_r4dwpn_224K_8_4_P4, xpass2_r4_8_levels, 8, 4);
r4dwpn_pass1sc320(xfft_r4dwpn_240K_768_4_P4, xpass2_r4_768_levels, 768, 4);
r4dwpn_pass1sc384(xfft_r4dwpn_288K_768_4_P4, xpass2_r4_768_levels, 768, 4);
r4dwpn_pass1sc256(xfft_r4dwpn_320K_1280_4_P4, xpass2_r4_1280_levels, 1280, 4);
r4dwpn_pass1sc320(xfft_r4dwpn_320K_10_4_P4, xpass2_r4_10_levels, 10, 4);
r4dwpn_pass1sc448(xfft_r4dwpn_336K_768_4_P4, xpass2_r4_768_levels, 768, 4);
r4dwpn_pass1sc256(xfft_r4dwpn_384K_1536_4_P4, xpass2_r4_1536_levels, 1536, 4);
r4dwpn_pass1sc320(xfft_r4dwpn_400K_1280_4_P4, xpass2_r4_1280_levels, 1280, 4);
r4dwpn_pass1sc448(xfft_r4dwpn_448K_10_4_P4, xpass2_r4_10_levels, 10, 4);
r4dwpn_pass1sc640(xfft_r4dwpn_480K_768_4_P4, xpass2_r4_768_levels, 768, 4);
r4dwpn_pass1sc256(xfft_r4dwpn_512K_11_4_P4, xpass2_r4_11_levels, 11, 4);
r4dwpn_pass1sc448(xfft_r4dwpn_560K_1280_4_P4, xpass2_r4_1280_levels, 1280, 4);
r4dwpn_pass1sc384(xfft_r4dwpn_576K_1536_4_P4, xpass2_r4_1536_levels, 1536, 4);
r4dwpn_pass1sc640(xfft_r4dwpn_640K_10_4_P4, xpass2_r4_10_levels, 10, 4);
r4dwpn_pass1sc896(xfft_r4dwpn_672K_768_4_P4, xpass2_r4_768_levels, 768, 4);
r4dwpn_pass1sc320(xfft_r4dwpn_720K_2304_4_P4, xpass2_r4_2304_levels, 2304, 4);
r4dwpn_pass1sc640(xfft_r4dwpn_800K_1280_4_P4, xpass2_r4_1280_levels, 1280, 4);
r4dwpn_pass1sc896(xfft_r4dwpn_896K_10_4_P4, xpass2_r4_10_levels, 10, 4);
r4dwpn_pass1sc640(xfft_r4dwpn_960K_1536_4_P4, xpass2_r4_1536_levels, 1536, 4);
r4dwpn_pass1sc256(xfft_r4dwpn_1M_12_4_P4, xpass2_r4_12_levels, 12, 4);
#endif

