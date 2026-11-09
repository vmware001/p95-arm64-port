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
#include "r4pass1sc.h"

//void xpass2_r4_1280_levels_asm(struct gwasm_data *__restrict g);

r4_pass1sc32(xfft_r4_8K_np_8_4_P4, xpass2_r4_8_levels, 8, 4);
r4_pass1sc64(xfft_r4_16K_8_4_P4, xpass2_r4_8_levels, 8, 4);
r4_pass1sc64(xfft_r4_16K_np_8_4_P4, xpass2_r4_8_levels, 8, 4);
r4_pass1sc80(xfft_r4_20K_np_8_4_P4, xpass2_r4_8_levels, 8, 4);
r4_pass1sc112(xfft_r4_28K_8_4_P4, xpass2_r4_8_levels, 8, 4);
r4_pass1sc32(xfft_r4_40K_np_1280_4_P4, xpass2_r4_1280_levels, 1280, 4);
r4_pass1sc64(xfft_r4_48K_768_4_P4, xpass2_r4_768_levels, 768, 4);
r4_pass1sc224(xfft_r4_56K_8_4_P4, xpass2_r4_8_levels, 8, 4);
r4_pass1sc112(xfft_r4_84K_768_4_P4, xpass2_r4_768_levels, 768, 4);

r4_pass1sc32ac(xfft_r4_8K_ac_np_8_4_P4, xpass2_r4_8_levels, 8, 4);
r4_pass1sc32ac(xfft_r4_40K_ac_np_1280_4_P4, xpass2_r4_1280_levels, 1280, 4);
r4_pass1sc64ac(xfft_r4_48K_ac_np_768_4_P4, xpass2_r4_768_levels, 768, 4);
r4_pass1sc64ac(xfft_r4_80K_ac_1280_4_P4, xpass2_r4_1280_levels, 1280, 4);


#if 0
r4_pass1sc32(xfft_r4_8K_8_4_P4, xpass2_r4_8_levels, 8, 4);
r4_pass1sc64(xfft_r4_16K_8_4_P4, xpass2_r4_8_levels, 8, 4);
r4_pass1sc80(xfft_r4_20K_8_4_P4, xpass2_r4_8_levels, 8, 4);
r4_pass1sc96(xfft_r4_24K_8_4_P4, xpass2_r4_8_levels, 8, 4);
r4_pass1sc112(xfft_r4_28K_8_4_P4, xpass2_r4_8_levels, 8, 4);
r4_pass1sc128(xfft_r4_32K_8_4_P4, xpass2_r4_8_levels, 8, 4);
r4_pass1sc224(xfft_r4_56K_8_4_P4, xpass2_r4_8_levels, 8, 4);
r4_pass1sc256(xfft_r4_64K_8_4_P4, xpass2_r4_8_levels, 8, 4);
r4_pass1sc32(xfft_r4_24K_768_4_P4, xpass2_r4_768_levels, 768, 4);
r4_pass1sc64(xfft_r4_48K_768_4_P4, xpass2_r4_768_levels, 768, 4);
r4_pass1sc80(xfft_r4_60K_768_4_P4, xpass2_r4_768_levels, 768, 4);
r4_pass1sc96(xfft_r4_72K_768_4_P4, xpass2_r4_768_levels, 768, 4);
r4_pass1sc112(xfft_r4_84K_768_4_P4, xpass2_r4_768_levels, 768, 4);
r4_pass1sc128(xfft_r4_96K_768_4_P4, xpass2_r4_768_levels, 768, 4);
r4_pass1sc224(xfft_r4_168K_768_4_P4, xpass2_r4_768_levels, 768, 4);
r4_pass1sc256(xfft_r4_192K_768_4_P4, xpass2_r4_768_levels, 768, 4);
r4_pass1sc32(xfft_r4_32K_10_4_P4, xpass2_r4_10_levels, 10, 4);
r4_pass1sc64(xfft_r4_64K_10_4_P4, xpass2_r4_10_levels, 10, 4);
r4_pass1sc80(xfft_r4_80K_10_4_P4, xpass2_r4_10_levels, 10, 4);
r4_pass1sc96(xfft_r4_96K_10_4_P4, xpass2_r4_10_levels, 10, 4);
r4_pass1sc112(xfft_r4_112K_10_4_P4, xpass2_r4_10_levels, 10, 4);
r4_pass1sc128(xfft_r4_128K_10_4_P4, xpass2_r4_10_levels, 10, 4);
r4_pass1sc224(xfft_r4_224K_10_4_P4, xpass2_r4_10_levels, 10, 4);
r4_pass1sc256(xfft_r4_256K_10_4_P4, xpass2_r4_10_levels, 10, 4);
r4_pass1sc32(xfft_r4_40K_1280_4_P4, xpass2_r4_1280_levels, 1280, 4);
r4_pass1sc64(xfft_r4_80K_1280_4_P4, xpass2_r4_1280_levels, 1280, 4);
r4_pass1sc80(xfft_r4_100K_1280_4_P4, xpass2_r4_1280_levels, 1280, 4);
r4_pass1sc96(xfft_r4_120K_1280_4_P4, xpass2_r4_1280_levels, 1280, 4);
r4_pass1sc112(xfft_r4_140K_1280_4_P4, xpass2_r4_1280_levels, 1280, 4);
r4_pass1sc128(xfft_r4_160K_1280_4_P4, xpass2_r4_1280_levels, 1280, 4);
r4_pass1sc224(xfft_r4_280K_1280_4_P4, xpass2_r4_1280_levels, 1280, 4);
r4_pass1sc256(xfft_r4_320K_1280_4_P4, xpass2_r4_1280_levels, 1280, 4);
r4_pass1sc32(xfft_r4_48K_1536_4_P4, xpass2_r4_1536_levels, 1536, 4);
r4_pass1sc64(xfft_r4_96K_1536_4_P4, xpass2_r4_1536_levels, 1536, 4);
r4_pass1sc80(xfft_r4_120K_1536_4_P4, xpass2_r4_1536_levels, 1536, 4);
r4_pass1sc96(xfft_r4_144K_1536_4_P4, xpass2_r4_1536_levels, 1536, 4);
r4_pass1sc112(xfft_r4_168K_1536_4_P4, xpass2_r4_1536_levels, 1536, 4);
r4_pass1sc128(xfft_r4_192K_1536_4_P4, xpass2_r4_1536_levels, 1536, 4);
r4_pass1sc224(xfft_r4_336K_1536_4_P4, xpass2_r4_1536_levels, 1536, 4);
r4_pass1sc256(xfft_r4_384K_1536_4_P4, xpass2_r4_1536_levels, 1536, 4);

r4_pass1sc32ac(xfft_r4_8K_ac_8_4_P4, xpass2_r4_8_levels, 8, 4);
r4_pass1sc64ac(xfft_r4_16K_ac_8_4_P4, xpass2_r4_8_levels, 8, 4);
r4_pass1sc96ac(xfft_r4_24K_ac_8_4_P4, xpass2_r4_8_levels, 8, 4);
r4_pass1sc128ac(xfft_r4_32K_ac_8_4_P4, xpass2_r4_8_levels, 8, 4);
r4_pass1sc256ac(xfft_r4_64K_ac_8_4_P4, xpass2_r4_8_levels, 8, 4);
r4_pass1sc32ac(xfft_r4_24K_ac_768_4_P4, xpass2_r4_768_levels, 768, 4);
r4_pass1sc64ac(xfft_r4_48K_ac_768_4_P4, xpass2_r4_768_levels, 768, 4);
r4_pass1sc96ac(xfft_r4_72K_ac_768_4_P4, xpass2_r4_768_levels, 768, 4);
r4_pass1sc128ac(xfft_r4_96K_ac_768_4_P4, xpass2_r4_768_levels, 768, 4);
r4_pass1sc256ac(xfft_r4_192K_ac_768_4_P4, xpass2_r4_768_levels, 768, 4);
r4_pass1sc32ac(xfft_r4_32K_ac_10_4_P4, xpass2_r4_10_levels, 10, 4);
r4_pass1sc64ac(xfft_r4_64K_ac_10_4_P4, xpass2_r4_10_levels, 10, 4);
r4_pass1sc96ac(xfft_r4_96K_ac_10_4_P4, xpass2_r4_10_levels, 10, 4);
r4_pass1sc128ac(xfft_r4_128K_ac_10_4_P4, xpass2_r4_10_levels, 10, 4);
r4_pass1sc256ac(xfft_r4_256K_ac_10_4_P4, xpass2_r4_10_levels, 10, 4);
r4_pass1sc32ac(xfft_r4_40K_ac_1280_4_P4, xpass2_r4_1280_levels, 1280, 4);
r4_pass1sc64ac(xfft_r4_80K_ac_1280_4_P4, xpass2_r4_1280_levels, 1280, 4);
r4_pass1sc96ac(xfft_r4_120K_ac_1280_4_P4, xpass2_r4_1280_levels, 1280, 4);
r4_pass1sc128ac(xfft_r4_160K_ac_1280_4_P4, xpass2_r4_1280_levels, 1280, 4);
r4_pass1sc256ac(xfft_r4_320K_ac_1280_4_P4, xpass2_r4_1280_levels, 1280, 4);
r4_pass1sc32ac(xfft_r4_48K_ac_1536_4_P4, xpass2_r4_1536_levels, 1536, 4);
r4_pass1sc64ac(xfft_r4_96K_ac_1536_4_P4, xpass2_r4_1536_levels, 1536, 4);
r4_pass1sc96ac(xfft_r4_144K_ac_1536_4_P4, xpass2_r4_1536_levels, 1536, 4);
r4_pass1sc128ac(xfft_r4_192K_ac_1536_4_P4, xpass2_r4_1536_levels, 1536, 4);
r4_pass1sc256ac(xfft_r4_384K_ac_1536_4_P4, xpass2_r4_1536_levels, 1536, 4);

#endif
