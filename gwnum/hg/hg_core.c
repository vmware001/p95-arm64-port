#include <stdlib.h>
#include <ctype.h>

#include "../gwnum.h"
#include "../gwcommon.h"
#include "../unravel.h"
#include "../gwtables.h"
#include "../xarch.h"
#include "../xbasics.h"
#include "../xmult.h"
#include "hg.h"
#include "hgonepass.h"
#include "hgpass1.h"
#include "hgpass1sc.h"
#include "hgpass2.h"

#include "../gwnum_func.h"
#include <stdio.h>


// function name, symbol name
// complex / negacyclic: _ac
// pass1 in-place: _ip
// pass1sc: no in-place, no _ip
// _np: no prefetch


#if 0  // test with pass2 assembly
void xpass2_hg_8_levels_asm(struct gwasm_data *__restrict g) {
    if(g->ffttype == 1)
    {
        printf("g->NEGACYCLIC_FFT = %i\n", g->NEGACYCLIC_FFT);
        xpass2_hg_8_levels(g);
        return;
    }

    register struct gwasm_data *r11_arg __asm("r11") = g;
    __asm volatile (
    /* ".intel_syntax noprefix\n" */
    "call xpass2_hg_8_levels_P4\n"
    /* ".att_syntax\n" */
    : : /*"r" (g->NORMRTN),*/ "r" (r11_arg)
    : "rax", "rcx", "rdx", "rbx", "rsi", "rdi", "rbp",
    "r8", "r9", "r10", "memory",
    "xmm0", "xmm1", "xmm2", "xmm3", "xmm4", "xmm5", "xmm6", "xmm7",
    "xmm8", "xmm9", "xmm10", "xmm11", "xmm12", "xmm13", "xmm14", "xmm15" );
}
#endif

hg_pass1levels6pfa5(xfft_hg_10K_ip_8_4_CORE, xpass2_hg_8_levels, 8, 4);
hg_pass1levels6pfa6(xfft_hg_12K_ip_8_4_CORE, xpass2_hg_8_levels, 8, 4);
hg_pass1levels6pfa7(xfft_hg_14K_ip_8_4_CORE, xpass2_hg_8_levels, 8, 4);
//hg_pass1levels6pfa8(xfft_hg_16K_ip_8_4_CORE, xpass2_hg_8_levels, 8, 4);
//hg_pass1levels7pfa5(xfft_hg_20K_ip_8_4_CORE, xpass2_hg_8_levels, 8, 4);
//hg_pass1levels7pfa6(xfft_hg_24K_ip_8_4_CORE, xpass2_hg_8_levels, 8, 4);
hg_pass1levels7pfa7(xfft_hg_28K_ip_8_4_CORE, xpass2_hg_8_levels, 8, 4);
hg_pass1levels7pfa8(xfft_hg_32K_ip_8_4_CORE, xpass2_hg_8_levels, 8, 4);
hg_pass1levels5pfa6(xfft_hg_48K_ip_11_4_CORE, xpass2_hg_11_levels, 11, 4);

//hg_pass1levels5pfa6(xfft_hg_48K_11_4_CORE, xpass2_hg_11_levels, 11, 4);
hg_pass1sclevels9pfa6(xfft_hg_1536K_12_4_CORE, xpass2_hg_12_levels, 12, 4);

hg_pass1levels6complex3(xfft_hg_12K_ac_ip_8_4_CORE, xpass2_hg_8_levels, 8, 4);
//hg_pass1levels6complex4(xfft_hg_16K_ac_ip_8_4_CORE, xpass2_hg_8_levels, 8, 4);
//hg_pass1levels7complex3(xfft_hg_24K_ac_ip_8_4_CORE, xpass2_hg_8_levels, 8, 4);
//hg_pass1levels7complex4(xfft_hg_32K_ac_ip_8_4_CORE, xpass2_hg_8_levels, 8, 4);


// xfft_hg_768K_12_2_CORE
