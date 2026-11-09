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

r4dwpn_pass1sc384(xfft_r4dwpn_1152K_3072_1_K8, xpass2_r4_3072_levels, 3072, 1);
r4dwpn_pass1sc320(xfft_r4dwpn_1200K_3840_2_K8, xpass2_r4_3840_levels, 3840, 2);
r4dwpn_pass1sc256(xfft_r4dwpn_1280K_5120_2_K8, xpass2_r4_5120_levels, 5120, 2);
r4dwpn_pass1sc448(xfft_r4dwpn_1344K_3072_1_K8, xpass2_r4_3072_levels, 3072, 1);
r4dwpn_pass1sc384(xfft_r4dwpn_1440K_3840_2_K8, xpass2_r4_3840_levels, 3840, 2);
r4dwpn_pass1sc384(xfft_r4dwpn_1536K_12_2_K8, xpass2_r4_12_levels, 12, 2);
r4dwpn_pass1sc320(xfft_r4dwpn_1600K_5120_2_K8, xpass2_r4_5120_levels, 5120, 2);
r4dwpn_pass1sc768(xfft_r4dwpn_1728K_2304_2_K8, xpass2_r4_2304_levels, 2304, 2);
r4dwpn_pass1sc448(xfft_r4dwpn_1792K_12_2_K8, xpass2_r4_12_levels, 12, 2);
r4dwpn_pass1sc384(xfft_r4dwpn_1920K_5120_2_K8, xpass2_r4_5120_levels, 5120, 2);
r4dwpn_pass1sc512(xfft_r4dwpn_2M_12_4_K8, xpass2_r4_12_levels, 12, 4);
r4dwpn_pass1sc1024(xfft_r4dwpn_2304K_2304_4_K8, xpass2_r4_2304_levels, 2304, 4);
r4dwpn_pass1sc384(xfft_r4dwpn_2400K_6400_1_K8, xpass2_r4_6400_levels, 6400, 1);
r4dwpn_pass1sc640(xfft_r4dwpn_2560K_12_2_K8, xpass2_r4_12_levels, 12, 2);
r4dwpn_pass1sc896(xfft_r4dwpn_2688K_3072_1_K8, xpass2_r4_3072_levels, 3072, 1);
r4dwpn_pass1sc768(xfft_r4dwpn_2880K_3840_2_K8, xpass2_r4_3840_levels, 3840, 2);
r4dwpn_pass1sc768(xfft_r4dwpn_3M_12_1_K8, xpass2_r4_12_levels, 12, 1);
r4dwpn_pass1sc640(xfft_r4dwpn_3200K_5120_1_K8, xpass2_r4_5120_levels, 5120, 1);
r4dwpn_pass1sc1536(xfft_r4dwpn_3456K_2304_4_K8, xpass2_r4_2304_levels, 2304, 4);
r4dwpn_pass1sc896(xfft_r4dwpn_3584K_12_1_K8, xpass2_r4_12_levels, 12, 1);
r4dwpn_pass1sc768(xfft_r4dwpn_3840K_5120_1_K8, xpass2_r4_5120_levels, 5120, 1);
r4dwpn_pass1sc1024(xfft_r4dwpn_4M_12_2_K8, xpass2_r4_12_levels, 12, 2);
r4dwpn_pass1sc512(xfft_r4dwpn_4608K_9216_4_K8, xpass2_r4_9216_levels, 9216, 4);
r4dwpn_pass1sc768(xfft_r4dwpn_4800K_6400_4_K8, xpass2_r4_6400_levels, 6400, 4);
r4dwpn_pass1sc1024(xfft_r4dwpn_5M_5120_4_K8, xpass2_r4_5120_levels, 5120, 4);
r4dwpn_pass1sc448(xfft_r4dwpn_5376K_12288_4_K8, xpass2_r4_12288_levels, 12288, 4);
r4dwpn_pass1sc1536(xfft_r4dwpn_5760K_3840_4_K8, xpass2_r4_3840_levels, 3840, 4);
r4dwpn_pass1sc1536(xfft_r4dwpn_6M_12_4_K8, xpass2_r4_12_levels, 12, 4);
