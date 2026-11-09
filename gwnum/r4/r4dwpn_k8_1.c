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

r4dwpn_pass1sc128(xfft_r4dwpn_32K_8_2_K8, xpass2_r4_8_levels, 8, 2);
r4dwpn_pass1sc256(xfft_r4dwpn_64K_8_1_K8, xpass2_r4_8_levels, 8, 1);
r4dwpn_pass1sc320(xfft_r4dwpn_80K_8_1_K8, xpass2_r4_8_levels, 8, 1);
r4dwpn_pass1sc128(xfft_r4dwpn_96K_768_2_K8, xpass2_r4_768_levels, 768, 2);
r4dwpn_pass1sc448(xfft_r4dwpn_112K_8_2_K8, xpass2_r4_8_levels, 8, 2);
r4dwpn_pass1sc128(xfft_r4dwpn_128K_10_2_K8, xpass2_r4_10_levels, 10, 2);
r4dwpn_pass1sc128(xfft_r4dwpn_160K_1280_2_K8, xpass2_r4_1280_levels, 1280, 2);
r4dwpn_pass1sc256(xfft_r4dwpn_192K_768_2_K8, xpass2_r4_768_levels, 768, 2);
r4dwpn_pass1sc896(xfft_r4dwpn_224K_8_1_K8, xpass2_r4_8_levels, 8, 1);
r4dwpn_pass1sc320(xfft_r4dwpn_240K_768_2_K8, xpass2_r4_768_levels, 768, 2);
r4dwpn_pass1sc256(xfft_r4dwpn_256K_10_2_K8, xpass2_r4_10_levels, 10, 2);
r4dwpn_pass1sc384(xfft_r4dwpn_288K_768_2_K8, xpass2_r4_768_levels, 768, 2);
r4dwpn_pass1sc256(xfft_r4dwpn_320K_1280_1_K8, xpass2_r4_1280_levels, 1280, 1);
r4dwpn_pass1sc448(xfft_r4dwpn_336K_768_2_K8, xpass2_r4_768_levels, 768, 2);
r4dwpn_pass1sc384(xfft_r4dwpn_384K_10_2_K8, xpass2_r4_10_levels, 10, 2);
r4dwpn_pass1sc320(xfft_r4dwpn_400K_1280_1_K8, xpass2_r4_1280_levels, 1280, 1);
r4dwpn_pass1sc448(xfft_r4dwpn_448K_10_2_K8, xpass2_r4_10_levels, 10, 2);
r4dwpn_pass1sc384(xfft_r4dwpn_480K_1280_2_K8, xpass2_r4_1280_levels, 1280, 2);
r4dwpn_pass1sc512(xfft_r4dwpn_512K_10_2_K8, xpass2_r4_10_levels, 10, 2);
r4dwpn_pass1sc256(xfft_r4dwpn_576K_2304_4_K8, xpass2_r4_2304_levels, 2304, 4);
r4dwpn_pass1sc640(xfft_r4dwpn_640K_10_1_K8, xpass2_r4_10_levels, 10, 1);
r4dwpn_pass1sc896(xfft_r4dwpn_672K_768_2_K8, xpass2_r4_768_levels, 768, 2);
r4dwpn_pass1sc256(xfft_r4dwpn_768K_3072_1_K8, xpass2_r4_3072_levels, 3072, 1);
r4dwpn_pass1sc640(xfft_r4dwpn_800K_1280_1_K8, xpass2_r4_1280_levels, 1280, 1);
r4dwpn_pass1sc384(xfft_r4dwpn_864K_2304_2_K8, xpass2_r4_2304_levels, 2304, 2);
r4dwpn_pass1sc896(xfft_r4dwpn_896K_10_1_K8, xpass2_r4_10_levels, 10, 1);
r4dwpn_pass1sc256(xfft_r4dwpn_960K_3840_2_K8, xpass2_r4_3840_levels, 3840, 2);
r4dwpn_pass1sc256(xfft_r4dwpn_1M_12_2_K8, xpass2_r4_12_levels, 12, 2);
