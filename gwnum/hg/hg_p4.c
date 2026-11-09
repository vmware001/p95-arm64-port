#include <stdlib.h>
#include <ctype.h>

#include "../gwnum.h"
#include "../gwcommon.h"
#include "../unravel.h"
#include "../gwtables.h"
#include "../xarch.h"
#include "../xbasics.h"
#include "../xmult.h"
#include "hg.h"
#include "hgonepass.h"
#include "hgpass1.h"
#include "hgpass1sc.h"
#include "hgpass2.h"

#include "../gwnum_func.h"



// function name, symbol name
// complex / negacyclic: _ac
// pass1 in-place: _ip
// pass1sc: no in-place, no _ip
// _np: no prefetch
#if 0
hg_pass1levels6pfa5(xfft_hg_10K_ip_8_4_P4, xpass2_hg_8_levels, 8, 4);
hg_pass1levels6pfa6(xfft_hg_12K_ip_8_4_P4, xpass2_hg_8_levels, 8, 4);
hg_pass1levels6pfa7(xfft_hg_14K_ip_8_4_P4, xpass2_hg_8_levels, 8, 4);
hg_pass1levels6pfa8(xfft_hg_16K_ip_8_4_P4, xpass2_hg_8_levels, 8, 4);
hg_pass1levels7pfa5(xfft_hg_20K_ip_8_4_P4, xpass2_hg_8_levels, 8, 4);
hg_pass1levels7pfa6(xfft_hg_24K_ip_8_4_P4, xpass2_hg_8_levels, 8, 4);
hg_pass1levels7pfa7(xfft_hg_28K_ip_8_4_P4, xpass2_hg_8_levels, 8, 4);
hg_pass1levels7pfa8(xfft_hg_32K_ip_8_4_P4, xpass2_hg_8_levels, 8, 4);

hg_pass1levels6complex3(xfft_hg_12K_ac_ip_8_4_P4, xpass2_hg_8_levels, 8, 4);
hg_pass1levels6complex4(xfft_hg_16K_ac_ip_8_4_P4, xpass2_hg_8_levels, 8, 4);
hg_pass1levels7complex3(xfft_hg_24K_ac_ip_8_4_P4, xpass2_hg_8_levels, 8, 4);
hg_pass1levels7complex4(xfft_hg_32K_ac_ip_8_4_P4, xpass2_hg_8_levels, 8, 4);


hg_pass1levels5pfa5(xfft_hg_40K_ip_11_4_P4, xpass2_hg_11_levels, 11, 4);
hg_pass1levels5pfa6(xfft_hg_48K_ip_11_4_P4, xpass2_hg_11_levels, 11, 4);
hg_pass1levels5pfa7(xfft_hg_56K_ip_11_4_P4, xpass2_hg_11_levels, 11, 4);
hg_pass1levels5pfa8(xfft_hg_64K_ip_11_4_P4, xpass2_hg_11_levels, 11, 4);
hg_pass1levels6pfa5(xfft_hg_80K_ip_11_4_P4, xpass2_hg_11_levels, 11, 4);
hg_pass1levels6pfa6(xfft_hg_96K_ip_11_4_P4, xpass2_hg_11_levels, 11, 4);
hg_pass1levels6pfa7(xfft_hg_112K_ip_11_4_P4, xpass2_hg_11_levels, 11, 4);
hg_pass1levels6pfa8(xfft_hg_128K_ip_11_4_P4, xpass2_hg_11_levels, 11, 4);


hg_pass1levels5complex3(xfft_hg_48K_ac_ip_11_4_P4, xpass2_hg_11_levels, 11, 4);
hg_pass1levels5complex4(xfft_hg_64K_ac_ip_11_4_P4, xpass2_hg_11_levels, 11, 4);
hg_pass1levels6complex3(xfft_hg_96K_ac_ip_11_4_P4, xpass2_hg_11_levels, 11, 4);
hg_pass1levels6complex4(xfft_hg_128K_ac_ip_11_4_P4, xpass2_hg_11_levels, 11, 4);

#endif



#if 0
hg_pass1sclevels8pfa5(xfft_hg_160K_10_4_P4, xpass2_hg_10_levels, 10, 4);
hg_pass1sclevels8pfa5(xfft_hg_160K_10_2_P4, xpass2_hg_10_levels, 10, 2);
hg_pass1sclevels8pfa6(xfft_hg_192K_10_4_P4, xpass2_hg_10_levels, 10, 4);
hg_pass1sclevels8pfa6(xfft_hg_192K_10_2_P4, xpass2_hg_10_levels, 10, 2);
hg_pass1sclevels8pfa6(xfft_hg_192K_10_1_P4, xpass2_hg_10_levels, 10, 1);
hg_pass1sclevels8pfa7(xfft_hg_224K_10_4_P4, xpass2_hg_10_levels, 10, 4);
hg_pass1sclevels8pfa7(xfft_hg_224K_10_2_P4, xpass2_hg_10_levels, 10, 2);
hg_pass1sclevels8pfa8(xfft_hg_256K_10_4_P4, xpass2_hg_10_levels, 10, 4);

hg_pass1sclevels9pfa5(xfft_hg_320K_10_4_P4, xpass2_hg_10_levels, 10, 4);
hg_pass1sclevels9pfa6(xfft_hg_384K_10_4_P4, xpass2_hg_10_levels, 10, 4);
hg_pass1sclevels9pfa7(xfft_hg_448K_10_4_P4, xpass2_hg_10_levels, 10, 4);
hg_pass1sclevels9pfa8(xfft_hg_512K_10_4_P4, xpass2_hg_10_levels, 10, 4);


hg_pass1sclevels10pfa5(xfft_hg_640K_10_4_P4, xpass2_hg_10_levels, 10, 4);
hg_pass1sclevels10pfa6(xfft_hg_768K_10_4_P4, xpass2_hg_10_levels, 10, 4);
hg_pass1sclevels10pfa7(xfft_hg_896K_10_4_P4, xpass2_hg_10_levels, 10, 4);
hg_pass1sclevels10pfa8(xfft_hg_1024K_10_4_P4, xpass2_hg_10_levels, 10, 4);
hg_pass1sclevels11pfa5(xfft_hg_640K_11_4_P4, xpass2_hg_11_levels, 11, 4);
hg_pass1sclevels11pfa6(xfft_hg_768K_11_4_P4, xpass2_hg_11_levels, 11, 4);
hg_pass1sclevels11pfa7(xfft_hg_896K_11_4_P4, xpass2_hg_11_levels, 11, 4);
hg_pass1sclevels11pfa8(xfft_hg_1024K_11_4_P4, xpass2_hg_11_levels, 11, 4);
hg_pass1sclevels12pfa5(xfft_hg_640K_12_4_P4, xpass2_hg_12_levels, 12, 4);
hg_pass1sclevels12pfa6(xfft_hg_768K_12_4_P4, xpass2_hg_12_levels, 12, 4);
hg_pass1sclevels12pfa7(xfft_hg_896K_12_4_P4, xpass2_hg_12_levels, 12, 4);
hg_pass1sclevels12pfa8(xfft_hg_1024K_12_4_P4, xpass2_hg_12_levels, 12, 4);

hg_pass1sclevels8complex3(testfunc01, xpass2_hg_8_levels, 8, 4);
hg_pass1sclevels9complex3(testfunc02, xpass2_hg_8_levels, 9, 4);
hg_pass1sclevels10complex3(testfunc03, xpass2_hg_8_levels, 10, 4);
hg_pass1sclevels11complex3(testfunc04, xpass2_hg_8_levels, 11, 4);
hg_pass1sclevels12complex3(testfunc05, xpass2_hg_8_levels, 12, 4);
hg_pass1sclevels8complex4(testfunc41, xpass2_hg_8_levels, 8, 4);
hg_pass1sclevels9complex4(testfunc42, xpass2_hg_8_levels, 9, 4);
hg_pass1sclevels10complex4(testfunc43, xpass2_hg_8_levels, 10, 4);
hg_pass1sclevels11complex4(testfunc44, xpass2_hg_8_levels, 11, 4);
hg_pass1sclevels12complex4(testfunc45, xpass2_hg_8_levels, 12, 4);

#endif


#if 0
hg_pass1levels6pfa5(xfft_hg_10K_ip_8_4_P4, xpass2_hg_8_levels, 8, 4);
hg_pass1levels6pfa6(xfft_hg_12K_ip_8_4_P4, xpass2_hg_8_levels, 8, 4);
hg_pass1levels6pfa7(xfft_hg_14K_ip_8_4_P4, xpass2_hg_8_levels, 8, 4);
hg_pass1levels6pfa8(xfft_hg_16K_ip_8_4_P4, xpass2_hg_8_levels, 8, 4);
hg_pass1levels7pfa5(xfft_hg_20K_ip_8_4_P4, xpass2_hg_8_levels, 8, 4);
hg_pass1levels7pfa6(xfft_hg_24K_ip_8_4_P4, xpass2_hg_8_levels, 8, 4);
hg_pass1levels7pfa7(xfft_hg_28K_ip_8_4_P4, xpass2_hg_8_levels, 8, 4);
hg_pass1levels7pfa8(xfft_hg_32K_ip_8_4_P4, xpass2_hg_8_levels, 8, 4);

hg_pass1levels5pfa5(xfft_hg_40K_ip_8_4_P4, xpass2_hg_11_levels, 11, 4);
hg_pass1levels5pfa6(xfft_hg_48K_ip_8_4_P4, xpass2_hg_11_levels, 11, 4);
hg_pass1levels5pfa7(xfft_hg_56K_ip_8_4_P4, xpass2_hg_11_levels, 11, 4);
hg_pass1levels5pfa8(xfft_hg_64K_ip_8_4_P4, xpass2_hg_11_levels, 11, 4);
hg_pass1levels6pfa5(xfft_hg_80K_ip_8_4_P4, xpass2_hg_11_levels, 11, 4);
hg_pass1levels6pfa6(xfft_hg_96K_ip_8_4_P4, xpass2_hg_11_levels, 11, 4);
hg_pass1levels6pfa7(xfft_hg_112K_ip_8_4_P4, xpass2_hg_11_levels, 11, 4);
hg_pass1levels6pfa8(xfft_hg_128K_ip_8_4_P4, xpass2_hg_11_levels, 11, 4);

hg_pass1sclevels8pfa5(xfft_hg_160K_10_4_P4, xpass2_hg_10_levels, 10, 4);
hg_pass1sclevels8pfa6(xfft_hg_192K_10_4_P4, xpass2_hg_10_levels, 10, 4);
hg_pass1sclevels8pfa7(xfft_hg_224K_10_4_P4, xpass2_hg_10_levels, 10, 4);
hg_pass1sclevels8pfa8(xfft_hg_256K_10_4_P4, xpass2_hg_10_levels, 10, 4);
hg_pass1sclevels9pfa5(xfft_hg_320K_10_4_P4, xpass2_hg_10_levels, 10, 4);
hg_pass1sclevels9pfa6(xfft_hg_384K_10_4_P4, xpass2_hg_10_levels, 10, 4);
hg_pass1sclevels9pfa7(xfft_hg_448K_10_4_P4, xpass2_hg_10_levels, 10, 4);
hg_pass1sclevels9pfa8(xfft_hg_512K_10_4_P4, xpass2_hg_10_levels, 10, 4);

hg_pass1sclevels9pfa5(xfft_hg_640K_11_4_P4, xpass2_hg_11_levels, 11, 4);
hg_pass1sclevels9pfa6(xfft_hg_768K_11_4_P4, xpass2_hg_11_levels, 11, 4);
hg_pass1sclevels9pfa7(xfft_hg_896K_11_4_P4, xpass2_hg_11_levels, 11, 4);
hg_pass1sclevels9pfa8(xfft_hg_1024K_11_4_P4, xpass2_hg_11_levels, 11, 4);
hg_pass1sclevels10pfa5(xfft_hg_1280K_11_2_P4, xpass2_hg_11_levels, 11, 2);
hg_pass1sclevels10pfa6(xfft_hg_1546K_11_1_P4, xpass2_hg_11_levels, 11, 1);
hg_pass1sclevels10pfa7(xfft_hg_1792K_11_1_P4, xpass2_hg_11_levels, 11, 1);

hg_pass1sclevels8pfa5(xfft_hg_640K_12_4_P4, xpass2_hg_12_levels, 12, 4);
hg_pass1sclevels8pfa6(xfft_hg_768K_12_4_P4, xpass2_hg_12_levels, 12, 4);
hg_pass1sclevels8pfa7(xfft_hg_896K_12_4_P4, xpass2_hg_12_levels, 12, 4);
hg_pass1sclevels8pfa8(xfft_hg_1024K_12_4_P4, xpass2_hg_12_levels, 12, 4);


#endif
//void xpass2_hg_8_levels_asm(struct gwasm_data *__restrict g);

hg_pass1levels6pfa5(xfft_hg_10K_ip_8_4_P4, xpass2_hg_8_levels, 8, 4);
hg_pass1levels6pfa6(xfft_hg_12K_ip_8_4_P4, xpass2_hg_8_levels, 8, 4);
hg_pass1levels6pfa7(xfft_hg_14K_ip_8_4_P4, xpass2_hg_8_levels, 8, 4);
//hg_pass1levels6pfa8(xfft_hg_16K_ip_8_4_P4, xpass2_hg_8_levels, 8, 4);
//hg_pass1levels7pfa5(xfft_hg_20K_ip_8_4_P4, xpass2_hg_8_levels, 8, 4);
hg_pass1levels7pfa6(xfft_hg_24K_ip_8_4_P4, xpass2_hg_8_levels, 8, 4);
//hg_pass1levels7pfa7(xfft_hg_28K_ip_8_4_P4, xpass2_hg_8_levels, 8, 4);
hg_pass1levels7pfa8(xfft_hg_32K_ip_8_4_P4, xpass2_hg_8_levels, 8, 4);

hg_pass1sclevels8pfa6(xfft_hg_768K_12_2_P4, xpass2_hg_12_levels, 12, 2);

hg_pass1levels6complex3(xfft_hg_12K_ac_ip_8_4_P4, xpass2_hg_8_levels, 8, 4);
//hg_pass1levels6complex4(xfft_hg_16K_ac_ip_8_4_P4, xpass2_hg_8_levels, 8, 4);
hg_pass1levels7complex3(xfft_hg_24K_ac_ip_8_4_P4, xpass2_hg_8_levels, 8, 4);
hg_pass1levels7complex4(xfft_hg_32K_ac_ip_8_4_P4, xpass2_hg_8_levels, 8, 4);


// xfft_hg_768K_12_2_P4
