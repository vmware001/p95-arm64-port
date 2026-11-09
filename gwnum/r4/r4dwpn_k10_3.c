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

r4dwpn_pass1sc640(xfft_r4dwpn_3200K_5120_4_K10, xpass2_r4_5120_levels, 5120, 4);
r4dwpn_pass1sc768(xfft_r4dwpn_3456K_4608_4_K10, xpass2_r4_4608_levels, 4608, 4);
r4dwpn_pass1sc896(xfft_r4dwpn_3584K_12_4_K10, xpass2_r4_12_levels, 12, 4);
r4dwpn_pass1sc768(xfft_r4dwpn_3840K_5120_4_K10, xpass2_r4_5120_levels, 5120, 4);
r4dwpn_pass1sc1024(xfft_r4dwpn_4M_12_4_K10, xpass2_r4_12_levels, 12, 4);
r4dwpn_pass1sc896(xfft_r4dwpn_4480K_5120_4_K10, xpass2_r4_5120_levels, 5120, 4);
r4dwpn_pass1sc768(xfft_r4dwpn_4608K_6144_4_K10, xpass2_r4_6144_levels, 6144, 4);
r4dwpn_pass1sc768(xfft_r4dwpn_4800K_6400_4_K10, xpass2_r4_6400_levels, 6400, 4);
r4dwpn_pass1sc1280(xfft_r4dwpn_5M_12_4_K10, xpass2_r4_12_levels, 12, 4);
r4dwpn_pass1sc896(xfft_r4dwpn_5376K_6144_4_K10, xpass2_r4_6144_levels, 6144, 4);
r4dwpn_pass1sc896(xfft_r4dwpn_5600K_6400_4_K10, xpass2_r4_6400_levels, 6400, 4);
r4dwpn_pass1sc1280(xfft_r4dwpn_5760K_4608_4_K10, xpass2_r4_4608_levels, 4608, 4);
r4dwpn_pass1sc1536(xfft_r4dwpn_6M_12_4_K10, xpass2_r4_12_levels, 12, 4);
r4dwpn_pass1sc1280(xfft_r4dwpn_6400K_5120_4_K10, xpass2_r4_5120_levels, 5120, 4);
r4dwpn_pass1sc896(xfft_r4dwpn_6720K_7680_4_K10, xpass2_r4_7680_levels, 7680, 4);
r4dwpn_pass1sc768(xfft_r4dwpn_6912K_9216_4_K10, xpass2_r4_9216_levels, 9216, 4);
r4dwpn_pass1sc1792(xfft_r4dwpn_7M_12_4_K10, xpass2_r4_12_levels, 12, 4);
r4dwpn_pass1sc1280(xfft_r4dwpn_7680K_6144_4_K10, xpass2_r4_6144_levels, 6144, 4);
r4dwpn_pass1sc1536(xfft_r4dwpn_18M_12288_4_K10, xpass2_r4_12288_levels, 12288, 4);
r4dwpn_pass1sc1280(xfft_r4dwpn_19200K_15360_4_K10, xpass2_r4_15360_levels, 15360, 4);
