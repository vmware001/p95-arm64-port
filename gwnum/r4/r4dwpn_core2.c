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

r4dwpn_pass1sc448(xfft_r4dwpn_560K_1280_4_CORE, xpass2_r4_1280_levels, 1280, 4);
r4dwpn_pass1sc256(xfft_r4dwpn_576K_2304_4_CORE, xpass2_r4_2304_levels, 2304, 4);
r4dwpn_pass1sc384(xfft_r4dwpn_576K_1536_2_CORE, xpass2_r4_1536_levels, 1536, 2);
r4dwpn_pass1sc768(xfft_r4dwpn_576K_768_4_CORE, xpass2_r4_768_levels, 768, 4);
r4dwpn_pass1sc256(xfft_r4dwpn_640K_2560_4_CORE, xpass2_r4_2560_levels, 2560, 4);
r4dwpn_pass1sc640(xfft_r4dwpn_640K_10_4_CORE, xpass2_r4_10_levels, 10, 4);
r4dwpn_pass1sc640(xfft_r4dwpn_640K_10_1_CORE, xpass2_r4_10_levels, 10, 1);
r4dwpn_pass1sc448(xfft_r4dwpn_672K_1536_4_CORE, xpass2_r4_1536_levels, 1536, 4);
r4dwpn_pass1sc448(xfft_r4dwpn_672K_1536_1_CORE, xpass2_r4_1536_levels, 1536, 1);
r4dwpn_pass1sc896(xfft_r4dwpn_672K_768_4_CORE, xpass2_r4_768_levels, 768, 4);
r4dwpn_pass1sc320(xfft_r4dwpn_720K_2304_4_CORE, xpass2_r4_2304_levels, 2304, 4);
r4dwpn_pass1sc256(xfft_r4dwpn_768K_3072_4_CORE, xpass2_r4_3072_levels, 3072, 4);
r4dwpn_pass1sc256(xfft_r4dwpn_768K_3072_2_CORE, xpass2_r4_3072_levels, 3072, 2);
r4dwpn_pass1sc1024(xfft_r4dwpn_768K_768_4_CORE, xpass2_r4_768_levels, 768, 4);
r4dwpn_pass1sc320(xfft_r4dwpn_800K_2560_4_CORE, xpass2_r4_2560_levels, 2560, 4);
r4dwpn_pass1sc640(xfft_r4dwpn_800K_1280_4_CORE, xpass2_r4_1280_levels, 1280, 4);
r4dwpn_pass1sc640(xfft_r4dwpn_800K_1280_1_CORE, xpass2_r4_1280_levels, 1280, 1);
r4dwpn_pass1sc384(xfft_r4dwpn_864K_2304_4_CORE, xpass2_r4_2304_levels, 2304, 4);
r4dwpn_pass1sc448(xfft_r4dwpn_896K_11_2_CORE, xpass2_r4_11_levels, 11, 2);
r4dwpn_pass1sc896(xfft_r4dwpn_896K_10_4_CORE, xpass2_r4_10_levels, 10, 4);
r4dwpn_pass1sc256(xfft_r4dwpn_960K_3840_4_CORE, xpass2_r4_3840_levels, 3840, 4);
r4dwpn_pass1sc256(xfft_r4dwpn_960K_3840_2_CORE, xpass2_r4_3840_levels, 3840, 2);
r4dwpn_pass1sc1280(xfft_r4dwpn_960K_768_4_CORE, xpass2_r4_768_levels, 768, 4);
r4dwpn_pass1sc256(xfft_r4dwpn_1M_12_4_CORE, xpass2_r4_12_levels, 12, 4);
r4dwpn_pass1sc256(xfft_r4dwpn_1M_12_2_CORE, xpass2_r4_12_levels, 12, 2);
r4dwpn_pass1sc1024(xfft_r4dwpn_1M_10_4_CORE, xpass2_r4_10_levels, 10, 4);
r4dwpn_pass1sc448(xfft_r4dwpn_1120K_2560_4_CORE, xpass2_r4_2560_levels, 2560, 4);
r4dwpn_pass1sc256(xfft_r4dwpn_1152K_4608_4_CORE, xpass2_r4_4608_levels, 4608, 4);
r4dwpn_pass1sc384(xfft_r4dwpn_1152K_3072_4_CORE, xpass2_r4_3072_levels, 3072, 4);
r4dwpn_pass1sc384(xfft_r4dwpn_1152K_3072_2_CORE, xpass2_r4_3072_levels, 3072, 2);
r4dwpn_pass1sc320(xfft_r4dwpn_1200K_3840_4_CORE, xpass2_r4_3840_levels, 3840, 4);
r4dwpn_pass1sc320(xfft_r4dwpn_1200K_3840_2_CORE, xpass2_r4_3840_levels, 3840, 2);
r4dwpn_pass1sc256(xfft_r4dwpn_1280K_5120_4_CORE, xpass2_r4_5120_levels, 5120, 4);
r4dwpn_pass1sc320(xfft_r4dwpn_1280K_12_4_CORE, xpass2_r4_12_levels, 12, 4);
r4dwpn_pass1sc320(xfft_r4dwpn_1280K_12_2_CORE, xpass2_r4_12_levels, 12, 2);
r4dwpn_pass1sc640(xfft_r4dwpn_1280K_11_2_CORE, xpass2_r4_11_levels, 11, 2);
