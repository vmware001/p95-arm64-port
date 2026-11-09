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

r4dwpn_pass1sc448(xfft_r4dwpn_1120K_2560_4_P4, xpass2_r4_2560_levels, 2560, 4);
r4dwpn_pass1sc256(xfft_r4dwpn_1152K_4608_4_P4, xpass2_r4_4608_levels, 4608, 4);
r4dwpn_pass1sc320(xfft_r4dwpn_1200K_3840_4_P4, xpass2_r4_3840_levels, 3840, 4);
r4dwpn_pass1sc640(xfft_r4dwpn_1280K_11_4_P4, xpass2_r4_11_levels, 11, 4);
r4dwpn_pass1sc896(xfft_r4dwpn_1344K_1536_4_P4, xpass2_r4_1536_levels, 1536, 4);
r4dwpn_pass1sc320(xfft_r4dwpn_1440K_4608_4_P4, xpass2_r4_4608_levels, 4608, 4);
r4dwpn_pass1sc640(xfft_r4dwpn_1600K_2560_4_P4, xpass2_r4_2560_levels, 2560, 4);
r4dwpn_pass1sc448(xfft_r4dwpn_1680K_3840_4_P4, xpass2_r4_3840_levels, 3840, 4);
r4dwpn_pass1sc640(xfft_r4dwpn_1920K_3072_4_P4, xpass2_r4_3072_levels, 3072, 4);
r4dwpn_pass1sc320(xfft_r4dwpn_2000K_6400_4_P4, xpass2_r4_6400_levels, 6400, 4);
r4dwpn_pass1sc256(xfft_r4dwpn_2M_13_4_P4, xpass2_r4_13_levels, 13, 4);
r4dwpn_pass1sc256(xfft_r4dwpn_2304K_9216_4_P4, xpass2_r4_9216_levels, 9216, 4);
r4dwpn_pass1sc640(xfft_r4dwpn_2400K_3840_4_P4, xpass2_r4_3840_levels, 3840, 4);
r4dwpn_pass1sc640(xfft_r4dwpn_2560K_12_4_P4, xpass2_r4_12_levels, 12, 4);
r4dwpn_pass1sc896(xfft_r4dwpn_2688K_3072_4_P4, xpass2_r4_3072_levels, 3072, 4);
r4dwpn_pass1sc640(xfft_r4dwpn_2880K_4608_4_P4, xpass2_r4_4608_levels, 4608, 4);
r4dwpn_pass1sc768(xfft_r4dwpn_3M_12_4_P4, xpass2_r4_12_levels, 12, 4);
