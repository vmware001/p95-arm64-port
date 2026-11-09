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

r4dwpn_pass1sc128ac(xfft_r4dwpn_32K_ac_np_8_4_CORE, xpass2_r4_8_levels, 8, 4);
r4dwpn_pass1sc128ac(xfft_r4dwpn_32K_ac_8_2_CORE, xpass2_r4_8_levels, 8, 2);
r4dwpn_pass1sc256ac(xfft_r4dwpn_64K_ac_8_4_CORE, xpass2_r4_8_levels, 8, 4);
r4dwpn_pass1sc256ac(xfft_r4dwpn_64K_ac_np_8_4_CORE, xpass2_r4_8_levels, 8, 4);
r4dwpn_pass1sc256ac(xfft_r4dwpn_64K_ac_8_1_CORE, xpass2_r4_8_levels, 8, 1);
r4dwpn_pass1sc128ac(xfft_r4dwpn_96K_ac_768_4_CORE, xpass2_r4_768_levels, 768, 4);
r4dwpn_pass1sc128ac(xfft_r4dwpn_96K_ac_np_768_4_CORE, xpass2_r4_768_levels, 768, 4);
r4dwpn_pass1sc128ac(xfft_r4dwpn_96K_ac_768_2_CORE, xpass2_r4_768_levels, 768, 2);
r4dwpn_pass1sc384ac(xfft_r4dwpn_96K_ac_8_4_CORE, xpass2_r4_8_levels, 8, 4);
r4dwpn_pass1sc128ac(xfft_r4dwpn_128K_ac_10_4_CORE, xpass2_r4_10_levels, 10, 4);
r4dwpn_pass1sc128ac(xfft_r4dwpn_128K_ac_np_10_4_CORE, xpass2_r4_10_levels, 10, 4);
r4dwpn_pass1sc128ac(xfft_r4dwpn_128K_ac_10_2_CORE, xpass2_r4_10_levels, 10, 2);
r4dwpn_pass1sc128ac(xfft_r4dwpn_160K_ac_1280_4_CORE, xpass2_r4_1280_levels, 1280, 4);
r4dwpn_pass1sc128ac(xfft_r4dwpn_160K_ac_1280_2_CORE, xpass2_r4_1280_levels, 1280, 2);
r4dwpn_pass1sc256ac(xfft_r4dwpn_192K_ac_768_4_CORE, xpass2_r4_768_levels, 768, 4);
r4dwpn_pass1sc256ac(xfft_r4dwpn_192K_ac_768_1_CORE, xpass2_r4_768_levels, 768, 1);
r4dwpn_pass1sc768ac(xfft_r4dwpn_192K_ac_8_4_CORE, xpass2_r4_8_levels, 8, 4);
r4dwpn_pass1sc256ac(xfft_r4dwpn_256K_ac_10_4_CORE, xpass2_r4_10_levels, 10, 4);
r4dwpn_pass1sc256ac(xfft_r4dwpn_256K_ac_10_2_CORE, xpass2_r4_10_levels, 10, 2);
r4dwpn_pass1sc1024ac(xfft_r4dwpn_256K_ac_8_4_CORE, xpass2_r4_8_levels, 8, 4);
r4dwpn_pass1sc384ac(xfft_r4dwpn_288K_ac_768_4_CORE, xpass2_r4_768_levels, 768, 4);
r4dwpn_pass1sc384ac(xfft_r4dwpn_288K_ac_768_2_CORE, xpass2_r4_768_levels, 768, 2);
r4dwpn_pass1sc256ac(xfft_r4dwpn_320K_ac_1280_4_CORE, xpass2_r4_1280_levels, 1280, 4);
r4dwpn_pass1sc256ac(xfft_r4dwpn_320K_ac_1280_2_CORE, xpass2_r4_1280_levels, 1280, 2);
r4dwpn_pass1sc256ac(xfft_r4dwpn_384K_ac_1536_4_CORE, xpass2_r4_1536_levels, 1536, 4);
r4dwpn_pass1sc384ac(xfft_r4dwpn_384K_ac_10_4_CORE, xpass2_r4_10_levels, 10, 4);
r4dwpn_pass1sc384ac(xfft_r4dwpn_384K_ac_10_2_CORE, xpass2_r4_10_levels, 10, 2);
r4dwpn_pass1sc384ac(xfft_r4dwpn_480K_ac_1280_4_CORE, xpass2_r4_1280_levels, 1280, 4);
r4dwpn_pass1sc384ac(xfft_r4dwpn_480K_ac_1280_1_CORE, xpass2_r4_1280_levels, 1280, 1);
r4dwpn_pass1sc256ac(xfft_r4dwpn_512K_ac_11_4_CORE, xpass2_r4_11_levels, 11, 4);
r4dwpn_pass1sc512ac(xfft_r4dwpn_512K_ac_10_4_CORE, xpass2_r4_10_levels, 10, 4);
r4dwpn_pass1sc256ac(xfft_r4dwpn_576K_ac_2304_4_CORE, xpass2_r4_2304_levels, 2304, 4);
r4dwpn_pass1sc384ac(xfft_r4dwpn_576K_ac_1536_4_CORE, xpass2_r4_1536_levels, 1536, 4);
r4dwpn_pass1sc768ac(xfft_r4dwpn_576K_ac_768_4_CORE, xpass2_r4_768_levels, 768, 4);
r4dwpn_pass1sc256ac(xfft_r4dwpn_640K_ac_2560_4_CORE, xpass2_r4_2560_levels, 2560, 4);
r4dwpn_pass1sc512ac(xfft_r4dwpn_640K_ac_1280_4_CORE, xpass2_r4_1280_levels, 1280, 4);
r4dwpn_pass1sc256ac(xfft_r4dwpn_768K_ac_3072_4_CORE, xpass2_r4_3072_levels, 3072, 4);
r4dwpn_pass1sc384ac(xfft_r4dwpn_768K_ac_11_2_CORE, xpass2_r4_11_levels, 11, 2);
r4dwpn_pass1sc768ac(xfft_r4dwpn_768K_ac_10_4_CORE, xpass2_r4_10_levels, 10, 4);
r4dwpn_pass1sc640ac(xfft_r4dwpn_800K_ac_1280_4_CORE, xpass2_r4_1280_levels, 1280, 4);
r4dwpn_pass1sc640ac(xfft_r4dwpn_800K_ac_1280_1_CORE, xpass2_r4_1280_levels, 1280, 1);
r4dwpn_pass1sc384ac(xfft_r4dwpn_864K_ac_2304_4_CORE, xpass2_r4_2304_levels, 2304, 4);
r4dwpn_pass1sc256ac(xfft_r4dwpn_960K_ac_3840_4_CORE, xpass2_r4_3840_levels, 3840, 4);
r4dwpn_pass1sc256ac(xfft_r4dwpn_960K_ac_3840_2_CORE, xpass2_r4_3840_levels, 3840, 2);
r4dwpn_pass1sc768ac(xfft_r4dwpn_960K_ac_1280_4_CORE, xpass2_r4_1280_levels, 1280, 4);
r4dwpn_pass1sc256ac(xfft_r4dwpn_1M_ac_12_4_CORE, xpass2_r4_12_levels, 12, 4);
r4dwpn_pass1sc256ac(xfft_r4dwpn_1M_ac_12_2_CORE, xpass2_r4_12_levels, 12, 2);
r4dwpn_pass1sc256ac(xfft_r4dwpn_1152K_ac_4608_4_CORE, xpass2_r4_4608_levels, 4608, 4);
r4dwpn_pass1sc384ac(xfft_r4dwpn_1152K_ac_3072_4_CORE, xpass2_r4_3072_levels, 3072, 4);
r4dwpn_pass1sc384ac(xfft_r4dwpn_1152K_ac_3072_2_CORE, xpass2_r4_3072_levels, 3072, 2);
r4dwpn_pass1sc256ac(xfft_r4dwpn_1280K_ac_5120_4_CORE, xpass2_r4_5120_levels, 5120, 4);
r4dwpn_pass1sc256ac(xfft_r4dwpn_1280K_ac_5120_2_CORE, xpass2_r4_5120_levels, 5120, 2);
r4dwpn_pass1sc384ac(xfft_r4dwpn_1440K_ac_3840_4_CORE, xpass2_r4_3840_levels, 3840, 4);
r4dwpn_pass1sc384ac(xfft_r4dwpn_1440K_ac_3840_2_CORE, xpass2_r4_3840_levels, 3840, 2);
r4dwpn_pass1sc384ac(xfft_r4dwpn_1536K_ac_12_1_CORE, xpass2_r4_12_levels, 12, 1);
r4dwpn_pass1sc512ac(xfft_r4dwpn_1536K_ac_3072_4_CORE, xpass2_r4_3072_levels, 3072, 4);
r4dwpn_pass1sc768ac(xfft_r4dwpn_1536K_ac_11_4_CORE, xpass2_r4_11_levels, 11, 4);
r4dwpn_pass1sc256ac(xfft_r4dwpn_1600K_ac_6400_4_CORE, xpass2_r4_6400_levels, 6400, 4);
r4dwpn_pass1sc256ac(xfft_r4dwpn_1600K_ac_6400_2_CORE, xpass2_r4_6400_levels, 6400, 2);
r4dwpn_pass1sc384ac(xfft_r4dwpn_1728K_ac_4608_4_CORE, xpass2_r4_4608_levels, 4608, 4);
r4dwpn_pass1sc768ac(xfft_r4dwpn_1728K_ac_2304_4_CORE, xpass2_r4_2304_levels, 2304, 4);
r4dwpn_pass1sc384ac(xfft_r4dwpn_1920K_ac_5120_4_CORE, xpass2_r4_5120_levels, 5120, 4);
r4dwpn_pass1sc384ac(xfft_r4dwpn_1920K_ac_5120_1_CORE, xpass2_r4_5120_levels, 5120, 1);
r4dwpn_pass1sc512ac(xfft_r4dwpn_2M_ac_12_4_CORE, xpass2_r4_12_levels, 12, 4);
r4dwpn_pass1sc512ac(xfft_r4dwpn_2M_ac_12_2_CORE, xpass2_r4_12_levels, 12, 2);
r4dwpn_pass1sc1024ac(xfft_r4dwpn_2M_ac_11_4_CORE, xpass2_r4_11_levels, 11, 4);
