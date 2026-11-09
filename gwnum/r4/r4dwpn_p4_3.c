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

#include "r4.h"
#include "r5.h"
#include "r7.h"
#include "r8.h"
#include "r4pass2.h"
#include "r4dwpnpass1sc.h"

r4dwpn_pass1sc640(xfft_r4dwpn_3200K_5120_4_P4, xpass2_r4_5120_levels, 5120, 4);
r4dwpn_pass1sc896(xfft_r4dwpn_3360K_3840_4_P4, xpass2_r4_3840_levels, 3840, 4);
r4dwpn_pass1sc896(xfft_r4dwpn_3584K_12_4_P4, xpass2_r4_12_levels, 12, 4);
r4dwpn_pass1sc640(xfft_r4dwpn_3840K_6144_4_P4, xpass2_r4_6144_levels, 6144, 4);
r4dwpn_pass1sc640(xfft_r4dwpn_4000K_6400_4_P4, xpass2_r4_6400_levels, 6400, 4);
r4dwpn_pass1sc1024(xfft_r4dwpn_4M_12_4_P4, xpass2_r4_12_levels, 12, 4);
r4dwpn_pass1sc896(xfft_r4dwpn_4480K_5120_4_P4, xpass2_r4_5120_levels, 5120, 4);
r4dwpn_pass1sc1024(xfft_r4dwpn_4608K_4608_4_P4, xpass2_r4_4608_levels, 4608, 4);
r4dwpn_pass1sc640(xfft_r4dwpn_4800K_7680_4_P4, xpass2_r4_7680_levels, 7680, 4);
r4dwpn_pass1sc640(xfft_r4dwpn_5M_13_4_P4, xpass2_r4_13_levels, 13, 4);
r4dwpn_pass1sc896(xfft_r4dwpn_5376K_6144_4_P4, xpass2_r4_6144_levels, 6144, 4);
r4dwpn_pass1sc896(xfft_r4dwpn_5600K_6400_4_P4, xpass2_r4_6400_levels, 6400, 4);
r4dwpn_pass1sc640(xfft_r4dwpn_5760K_9216_4_P4, xpass2_r4_9216_levels, 9216, 4);
r4dwpn_pass1sc1024(xfft_r4dwpn_6M_6144_4_P4, xpass2_r4_6144_levels, 6144, 4);
