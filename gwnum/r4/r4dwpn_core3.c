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

r4dwpn_pass1sc448(xfft_r4dwpn_1344K_3072_4_CORE, xpass2_r4_3072_levels, 3072, 4);
r4dwpn_pass1sc448(xfft_r4dwpn_1344K_3072_1_CORE, xpass2_r4_3072_levels, 3072, 1);
r4dwpn_pass1sc896(xfft_r4dwpn_1344K_1536_4_CORE, xpass2_r4_1536_levels, 1536, 4);
r4dwpn_pass1sc320(xfft_r4dwpn_1440K_4608_4_CORE, xpass2_r4_4608_levels, 4608, 4);
r4dwpn_pass1sc384(xfft_r4dwpn_1440K_3840_2_CORE, xpass2_r4_3840_levels, 3840, 2);
r4dwpn_pass1sc640(xfft_r4dwpn_1440K_2304_4_CORE, xpass2_r4_2304_levels, 2304, 4);
r4dwpn_pass1sc384(xfft_r4dwpn_1536K_12_2_CORE, xpass2_r4_12_levels, 12, 2);
r4dwpn_pass1sc512(xfft_r4dwpn_1536K_3072_4_CORE, xpass2_r4_3072_levels, 3072, 4);
r4dwpn_pass1sc768(xfft_r4dwpn_1536K_11_4_CORE, xpass2_r4_11_levels, 11, 4);
r4dwpn_pass1sc256(xfft_r4dwpn_1600K_6400_4_CORE, xpass2_r4_6400_levels, 6400, 4);
r4dwpn_pass1sc320(xfft_r4dwpn_1600K_5120_4_CORE, xpass2_r4_5120_levels, 5120, 4);
r4dwpn_pass1sc320(xfft_r4dwpn_1600K_5120_2_CORE, xpass2_r4_5120_levels, 5120, 2);
r4dwpn_pass1sc448(xfft_r4dwpn_1680K_3840_4_CORE, xpass2_r4_3840_levels, 3840, 4);
r4dwpn_pass1sc384(xfft_r4dwpn_1728K_4608_4_CORE, xpass2_r4_4608_levels, 4608, 4);
r4dwpn_pass1sc768(xfft_r4dwpn_1728K_2304_4_CORE, xpass2_r4_2304_levels, 2304, 4);
r4dwpn_pass1sc448(xfft_r4dwpn_1792K_12_1_CORE, xpass2_r4_12_levels, 12, 1);
r4dwpn_pass1sc896(xfft_r4dwpn_1792K_11_4_CORE, xpass2_r4_11_levels, 11, 4);
r4dwpn_pass1sc384(xfft_r4dwpn_1920K_5120_1_CORE, xpass2_r4_5120_levels, 5120, 1);
r4dwpn_pass1sc640(xfft_r4dwpn_1920K_3072_4_CORE, xpass2_r4_3072_levels, 3072, 4);
r4dwpn_pass1sc512(xfft_r4dwpn_2M_12_4_CORE, xpass2_r4_12_levels, 12, 4);
r4dwpn_pass1sc512(xfft_r4dwpn_2M_12_2_CORE, xpass2_r4_12_levels, 12, 2);
r4dwpn_pass1sc1024(xfft_r4dwpn_2M_11_4_CORE, xpass2_r4_11_levels, 11, 4);
r4dwpn_pass1sc448(xfft_r4dwpn_2240K_5120_4_CORE, xpass2_r4_5120_levels, 5120, 4);
r4dwpn_pass1sc896(xfft_r4dwpn_2240K_2560_4_CORE, xpass2_r4_2560_levels, 2560, 4);
r4dwpn_pass1sc768(xfft_r4dwpn_2304K_3072_4_CORE, xpass2_r4_3072_levels, 3072, 4);
r4dwpn_pass1sc768(xfft_r4dwpn_2304K_3072_1_CORE, xpass2_r4_3072_levels, 3072, 1);
r4dwpn_pass1sc384(xfft_r4dwpn_2400K_6400_4_CORE, xpass2_r4_6400_levels, 6400, 4);
r4dwpn_pass1sc640(xfft_r4dwpn_2400K_3840_4_CORE, xpass2_r4_3840_levels, 3840, 4);
r4dwpn_pass1sc512(xfft_r4dwpn_2560K_5120_4_CORE, xpass2_r4_5120_levels, 5120, 4);
r4dwpn_pass1sc640(xfft_r4dwpn_2560K_12_4_CORE, xpass2_r4_12_levels, 12, 4);
r4dwpn_pass1sc640(xfft_r4dwpn_2560K_12_1_CORE, xpass2_r4_12_levels, 12, 1);
r4dwpn_pass1sc1280(xfft_r4dwpn_2560K_11_4_CORE, xpass2_r4_11_levels, 11, 4);
r4dwpn_pass1sc896(xfft_r4dwpn_2688K_3072_4_CORE, xpass2_r4_3072_levels, 3072, 4);
r4dwpn_pass1sc896(xfft_r4dwpn_2688K_3072_1_CORE, xpass2_r4_3072_levels, 3072, 1);
r4dwpn_pass1sc448(xfft_r4dwpn_2800K_6400_4_CORE, xpass2_r4_6400_levels, 6400, 4);
r4dwpn_pass1sc448(xfft_r4dwpn_2800K_6400_2_CORE, xpass2_r4_6400_levels, 6400, 2);
r4dwpn_pass1sc320(xfft_r4dwpn_2880K_9216_4_CORE, xpass2_r4_9216_levels, 9216, 4);
r4dwpn_pass1sc640(xfft_r4dwpn_2880K_4608_2_CORE, xpass2_r4_4608_levels, 4608, 2);
r4dwpn_pass1sc768(xfft_r4dwpn_2880K_3840_4_CORE, xpass2_r4_3840_levels, 3840, 4);
r4dwpn_pass1sc768(xfft_r4dwpn_3M_12_1_CORE, xpass2_r4_12_levels, 12, 1);
r4dwpn_pass1sc1024(xfft_r4dwpn_3M_3072_4_CORE, xpass2_r4_3072_levels, 3072, 4);
