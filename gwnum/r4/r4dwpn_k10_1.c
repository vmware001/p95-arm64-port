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

r4dwpn_pass1sc128(xfft_r4dwpn_32K_np_8_4_K10, xpass2_r4_8_levels, 8, 4);
r4dwpn_pass1sc256(xfft_r4dwpn_64K_8_2_K10, xpass2_r4_8_levels, 8, 2);
r4dwpn_pass1sc320(xfft_r4dwpn_80K_8_2_K10, xpass2_r4_8_levels, 8, 2);
r4dwpn_pass1sc128(xfft_r4dwpn_96K_768_4_K10, xpass2_r4_768_levels, 768, 4);
r4dwpn_pass1sc448(xfft_r4dwpn_112K_8_4_K10, xpass2_r4_8_levels, 8, 4);
r4dwpn_pass1sc128(xfft_r4dwpn_128K_10_4_K10, xpass2_r4_10_levels, 10, 4);
r4dwpn_pass1sc128(xfft_r4dwpn_160K_1280_2_K10, xpass2_r4_1280_levels, 1280, 2);
r4dwpn_pass1sc256(xfft_r4dwpn_192K_768_4_K10, xpass2_r4_768_levels, 768, 4);
r4dwpn_pass1sc896(xfft_r4dwpn_224K_8_4_K10, xpass2_r4_8_levels, 8, 4);
r4dwpn_pass1sc320(xfft_r4dwpn_240K_768_4_K10, xpass2_r4_768_levels, 768, 4);
r4dwpn_pass1sc256(xfft_r4dwpn_256K_10_4_K10, xpass2_r4_10_levels, 10, 4);
r4dwpn_pass1sc384(xfft_r4dwpn_288K_768_4_K10, xpass2_r4_768_levels, 768, 4);
r4dwpn_pass1sc320(xfft_r4dwpn_320K_10_1_K10, xpass2_r4_10_levels, 10, 1);
r4dwpn_pass1sc448(xfft_r4dwpn_336K_768_4_K10, xpass2_r4_768_levels, 768, 4);
r4dwpn_pass1sc384(xfft_r4dwpn_384K_10_4_K10, xpass2_r4_10_levels, 10, 4);
r4dwpn_pass1sc320(xfft_r4dwpn_400K_1280_1_K10, xpass2_r4_1280_levels, 1280, 1);
r4dwpn_pass1sc448(xfft_r4dwpn_448K_10_4_K10, xpass2_r4_10_levels, 10, 4);
r4dwpn_pass1sc384(xfft_r4dwpn_480K_1280_4_K10, xpass2_r4_1280_levels, 1280, 4);
r4dwpn_pass1sc512(xfft_r4dwpn_512K_10_4_K10, xpass2_r4_10_levels, 10, 4);
r4dwpn_pass1sc448(xfft_r4dwpn_560K_1280_4_K10, xpass2_r4_1280_levels, 1280, 4);
r4dwpn_pass1sc384(xfft_r4dwpn_576K_1536_4_K10, xpass2_r4_1536_levels, 1536, 4);
r4dwpn_pass1sc640(xfft_r4dwpn_640K_10_4_K10, xpass2_r4_10_levels, 10, 4);
r4dwpn_pass1sc448(xfft_r4dwpn_672K_1536_4_K10, xpass2_r4_1536_levels, 1536, 4);
r4dwpn_pass1sc320(xfft_r4dwpn_720K_2304_4_K10, xpass2_r4_2304_levels, 2304, 4);
r4dwpn_pass1sc768(xfft_r4dwpn_768K_10_4_K10, xpass2_r4_10_levels, 10, 4);
r4dwpn_pass1sc640(xfft_r4dwpn_800K_1280_4_K10, xpass2_r4_1280_levels, 1280, 4);
r4dwpn_pass1sc384(xfft_r4dwpn_864K_2304_4_K10, xpass2_r4_2304_levels, 2304, 4);
r4dwpn_pass1sc896(xfft_r4dwpn_896K_10_4_K10, xpass2_r4_10_levels, 10, 4);
r4dwpn_pass1sc640(xfft_r4dwpn_960K_1536_4_K10, xpass2_r4_1536_levels, 1536, 4);
r4dwpn_pass1sc256(xfft_r4dwpn_1M_12_4_K10, xpass2_r4_12_levels, 12, 4);
