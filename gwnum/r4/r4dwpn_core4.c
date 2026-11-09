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

r4dwpn_pass1sc640(xfft_r4dwpn_3200K_5120_4_CORE, xpass2_r4_5120_levels, 5120, 4);
r4dwpn_pass1sc640(xfft_r4dwpn_3200K_5120_2_CORE, xpass2_r4_5120_levels, 5120, 2);
r4dwpn_pass1sc640(xfft_r4dwpn_3200K_5120_1_CORE, xpass2_r4_5120_levels, 5120, 1);
r4dwpn_pass1sc896(xfft_r4dwpn_3360K_3840_4_CORE, xpass2_r4_3840_levels, 3840, 4);
r4dwpn_pass1sc384(xfft_r4dwpn_3456K_9216_4_CORE, xpass2_r4_9216_levels, 9216, 4);
r4dwpn_pass1sc768(xfft_r4dwpn_3456K_4608_4_CORE, xpass2_r4_4608_levels, 4608, 4);
r4dwpn_pass1sc896(xfft_r4dwpn_3584K_12_4_CORE, xpass2_r4_12_levels, 12, 4);
r4dwpn_pass1sc896(xfft_r4dwpn_3584K_12_2_CORE, xpass2_r4_12_levels, 12, 2);
r4dwpn_pass1sc3584(xfft_r4dwpn_3584K_10_4_CORE, xpass2_r4_10_levels, 10, 4);
r4dwpn_pass1sc768(xfft_r4dwpn_3840K_5120_4_CORE, xpass2_r4_5120_levels, 5120, 4);
r4dwpn_pass1sc768(xfft_r4dwpn_3840K_5120_1_CORE, xpass2_r4_5120_levels, 5120, 1);
r4dwpn_pass1sc1024(xfft_r4dwpn_3840K_3840_4_CORE, xpass2_r4_3840_levels, 3840, 4);
r4dwpn_pass1sc1024(xfft_r4dwpn_4M_12_4_CORE, xpass2_r4_12_levels, 12, 4);
r4dwpn_pass1sc1024(xfft_r4dwpn_4M_12_1_CORE, xpass2_r4_12_levels, 12, 1);
r4dwpn_pass1sc896(xfft_r4dwpn_4480K_5120_4_CORE, xpass2_r4_5120_levels, 5120, 4);
r4dwpn_pass1sc512(xfft_r4dwpn_4608K_9216_4_CORE, xpass2_r4_9216_levels, 9216, 4);
r4dwpn_pass1sc1024(xfft_r4dwpn_4608K_4608_4_CORE, xpass2_r4_4608_levels, 4608, 4);
r4dwpn_pass1sc768(xfft_r4dwpn_4800K_6400_4_CORE, xpass2_r4_6400_levels, 6400, 4);
r4dwpn_pass1sc640(xfft_r4dwpn_5M_13_4_CORE, xpass2_r4_13_levels, 13, 4);
r4dwpn_pass1sc1024(xfft_r4dwpn_5M_5120_4_CORE, xpass2_r4_5120_levels, 5120, 4);
r4dwpn_pass1sc896(xfft_r4dwpn_5376K_6144_4_CORE, xpass2_r4_6144_levels, 6144, 4);
r4dwpn_pass1sc1792(xfft_r4dwpn_5376K_3072_4_CORE, xpass2_r4_3072_levels, 3072, 4);
r4dwpn_pass1sc896(xfft_r4dwpn_5600K_6400_4_CORE, xpass2_r4_6400_levels, 6400, 4);
r4dwpn_pass1sc640(xfft_r4dwpn_5760K_9216_4_CORE, xpass2_r4_9216_levels, 9216, 4);
r4dwpn_pass1sc1280(xfft_r4dwpn_5760K_4608_4_CORE, xpass2_r4_4608_levels, 4608, 4);
r4dwpn_pass1sc1024(xfft_r4dwpn_6M_6144_4_CORE, xpass2_r4_6144_levels, 6144, 4);
r4dwpn_pass1sc1024(xfft_r4dwpn_6400K_6400_4_CORE, xpass2_r4_6400_levels, 6400, 4);
r4dwpn_pass1sc896(xfft_r4dwpn_6720K_7680_4_CORE, xpass2_r4_7680_levels, 7680, 4);
r4dwpn_pass1sc1792(xfft_r4dwpn_6720K_3840_4_CORE, xpass2_r4_3840_levels, 3840, 4);
r4dwpn_pass1sc768(xfft_r4dwpn_6912K_9216_4_CORE, xpass2_r4_9216_levels, 9216, 4);
r4dwpn_pass1sc1792(xfft_r4dwpn_7M_12_4_CORE, xpass2_r4_12_levels, 12, 4);
r4dwpn_pass1sc3584(xfft_r4dwpn_7M_11_4_CORE, xpass2_r4_11_levels, 11, 4);
r4dwpn_pass1sc1024(xfft_r4dwpn_7680K_7680_4_CORE, xpass2_r4_7680_levels, 7680, 4);
r4dwpn_pass1sc2560(xfft_r4dwpn_7680K_3072_4_CORE, xpass2_r4_3072_levels, 3072, 4);
r4dwpn_pass1sc640(xfft_r4dwpn_8000K_12800_4_CORE, xpass2_r4_12800_levels, 12800, 4);
r4dwpn_pass1sc1024(xfft_r4dwpn_8M_13_4_CORE, xpass2_r4_13_levels, 13, 4);
