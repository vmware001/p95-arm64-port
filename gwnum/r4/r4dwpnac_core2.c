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
#include "r8.h"
#include "r4pass2.h"
#include "r4dwpnpass1sc.h"

r4dwpn_pass1sc768ac(xfft_r4dwpn_2304K_ac_3072_4_CORE, xpass2_r4_3072_levels, 3072, 4);
r4dwpn_pass1sc768ac(xfft_r4dwpn_2304K_ac_3072_1_CORE, xpass2_r4_3072_levels, 3072, 1);
r4dwpn_pass1sc384ac(xfft_r4dwpn_2400K_ac_6400_4_CORE, xpass2_r4_6400_levels, 6400, 4);
r4dwpn_pass1sc384ac(xfft_r4dwpn_2400K_ac_6400_2_CORE, xpass2_r4_6400_levels, 6400, 2);
r4dwpn_pass1sc512ac(xfft_r4dwpn_2560K_ac_5120_4_CORE, xpass2_r4_5120_levels, 5120, 4);
r4dwpn_pass1sc512ac(xfft_r4dwpn_2560K_ac_5120_2_CORE, xpass2_r4_5120_levels, 5120, 2);
r4dwpn_pass1sc1280ac(xfft_r4dwpn_2560K_ac_11_4_CORE, xpass2_r4_11_levels, 11, 4);
r4dwpn_pass1sc768ac(xfft_r4dwpn_2880K_ac_3840_4_CORE, xpass2_r4_3840_levels, 3840, 4);
r4dwpn_pass1sc768ac(xfft_r4dwpn_2880K_ac_3840_1_CORE, xpass2_r4_3840_levels, 3840, 1);
r4dwpn_pass1sc768ac(xfft_r4dwpn_3M_ac_12_4_CORE, xpass2_r4_12_levels, 12, 4);
r4dwpn_pass1sc768ac(xfft_r4dwpn_3M_ac_12_1_CORE, xpass2_r4_12_levels, 12, 1);
r4dwpn_pass1sc512ac(xfft_r4dwpn_3200K_ac_6400_4_CORE, xpass2_r4_6400_levels, 6400, 4);
r4dwpn_pass1sc640ac(xfft_r4dwpn_3200K_ac_5120_4_CORE, xpass2_r4_5120_levels, 5120, 4);
r4dwpn_pass1sc384ac(xfft_r4dwpn_3456K_ac_9216_2_CORE, xpass2_r4_9216_levels, 9216, 2);
r4dwpn_pass1sc768ac(xfft_r4dwpn_3456K_ac_4608_4_CORE, xpass2_r4_4608_levels, 4608, 4);
r4dwpn_pass1sc768ac(xfft_r4dwpn_3840K_ac_5120_4_CORE, xpass2_r4_5120_levels, 5120, 4);
r4dwpn_pass1sc768ac(xfft_r4dwpn_3840K_ac_5120_1_CORE, xpass2_r4_5120_levels, 5120, 1);
r4dwpn_pass1sc256ac(xfft_r4dwpn_4M_ac_14_4_CORE, xpass2_r4_14_levels, 14, 4);
r4dwpn_pass1sc1024ac(xfft_r4dwpn_4M_ac_12_4_CORE, xpass2_r4_12_levels, 12, 4);
r4dwpn_pass1sc1024ac(xfft_r4dwpn_4M_ac_12_2_CORE, xpass2_r4_12_levels, 12, 2);
r4dwpn_pass1sc512ac(xfft_r4dwpn_4608K_ac_9216_4_CORE, xpass2_r4_9216_levels, 9216, 4);
r4dwpn_pass1sc1024ac(xfft_r4dwpn_4608K_ac_4608_4_CORE, xpass2_r4_4608_levels, 4608, 4);
r4dwpn_pass1sc768ac(xfft_r4dwpn_4800K_ac_6400_4_CORE, xpass2_r4_6400_levels, 6400, 4);
r4dwpn_pass1sc1024ac(xfft_r4dwpn_5M_ac_5120_4_CORE, xpass2_r4_5120_levels, 5120, 4);
r4dwpn_pass1sc768ac(xfft_r4dwpn_5760K_ac_7680_4_CORE, xpass2_r4_7680_levels, 7680, 4);
r4dwpn_pass1sc1280ac(xfft_r4dwpn_5760K_ac_4608_4_CORE, xpass2_r4_4608_levels, 4608, 4);
r4dwpn_pass1sc384ac(xfft_r4dwpn_6M_ac_14_4_CORE, xpass2_r4_14_levels, 14, 4);
r4dwpn_pass1sc1024ac(xfft_r4dwpn_6M_ac_6144_4_CORE, xpass2_r4_6144_levels, 6144, 4);
r4dwpn_pass1sc1024ac(xfft_r4dwpn_6400K_ac_6400_4_CORE, xpass2_r4_6400_levels, 6400, 4);
r4dwpn_pass1sc768ac(xfft_r4dwpn_6912K_ac_9216_4_CORE, xpass2_r4_9216_levels, 9216, 4);
r4dwpn_pass1sc768ac(xfft_r4dwpn_7680K_ac_10240_4_CORE, xpass2_r4_10240_levels, 10240, 4);
r4dwpn_pass1sc1024ac(xfft_r4dwpn_7680K_ac_7680_4_CORE, xpass2_r4_7680_levels, 7680, 4);
r4dwpn_pass1sc640ac(xfft_r4dwpn_8000K_ac_12800_4_CORE, xpass2_r4_12800_levels, 12800, 4);
r4dwpn_pass1sc1024ac(xfft_r4dwpn_8M_ac_13_4_CORE, xpass2_r4_13_levels, 13, 4);
r4dwpn_pass1sc1024ac(xfft_r4dwpn_9M_ac_9216_4_CORE, xpass2_r4_9216_levels, 9216, 4);
r4dwpn_pass1sc768ac(xfft_r4dwpn_9600K_ac_12800_4_CORE, xpass2_r4_12800_levels, 12800, 4);
r4dwpn_pass1sc1024ac(xfft_r4dwpn_10M_ac_10240_4_CORE, xpass2_r4_10240_levels, 10240, 4);
r4dwpn_pass1sc2560ac(xfft_r4dwpn_11520K_ac_4608_4_CORE, xpass2_r4_4608_levels, 4608, 4);
r4dwpn_pass1sc768ac(xfft_r4dwpn_12M_ac_14_4_CORE, xpass2_r4_14_levels, 14, 4);
r4dwpn_pass1sc1024ac(xfft_r4dwpn_12800K_ac_12800_4_CORE, xpass2_r4_12800_levels, 12800, 4);
r4dwpn_pass1sc3072ac(xfft_r4dwpn_13824K_ac_4608_4_CORE, xpass2_r4_4608_levels, 4608, 4);
r4dwpn_pass1sc1024ac(xfft_r4dwpn_15M_ac_15360_4_CORE, xpass2_r4_15360_levels, 15360, 4);
r4dwpn_pass1sc640ac(xfft_r4dwpn_16000K_ac_25600_4_CORE, xpass2_r4_25600_levels, 25600, 4);
r4dwpn_pass1sc4096ac(xfft_r4dwpn_16M_ac_12_4_CORE, xpass2_r4_12_levels, 12, 4);
r4dwpn_pass1sc4096ac(xfft_r4dwpn_18M_ac_4608_4_CORE, xpass2_r4_4608_levels, 4608, 4);
r4dwpn_pass1sc768ac(xfft_r4dwpn_19200K_ac_25600_4_CORE, xpass2_r4_25600_levels, 25600, 4);
r4dwpn_pass1sc4096ac(xfft_r4dwpn_20M_ac_5120_4_CORE, xpass2_r4_5120_levels, 5120, 4);
r4dwpn_pass1sc1536ac(xfft_r4dwpn_23040K_ac_15360_4_CORE, xpass2_r4_15360_levels, 15360, 4);
r4dwpn_pass1sc4096ac(xfft_r4dwpn_24M_ac_6144_4_CORE, xpass2_r4_6144_levels, 6144, 4);
r4dwpn_pass1sc1024ac(xfft_r4dwpn_25M_ac_25600_4_CORE, xpass2_r4_25600_levels, 25600, 4);
r4dwpn_pass1sc3072ac(xfft_r4dwpn_27M_ac_9216_4_CORE, xpass2_r4_9216_levels, 9216, 4);
r4dwpn_pass1sc4096ac(xfft_r4dwpn_30M_ac_7680_4_CORE, xpass2_r4_7680_levels, 7680, 4);
r4dwpn_pass1sc1280ac(xfft_r4dwpn_32000K_ac_25600_4_CORE, xpass2_r4_25600_levels, 25600, 4);
r4dwpn_pass1sc4096ac(xfft_r4dwpn_32M_ac_13_4_CORE, xpass2_r4_13_levels, 13, 4);
