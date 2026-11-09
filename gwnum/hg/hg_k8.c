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


hg_pass1levels6pfa5(xfft_hg_10K_ip_8_4_K8, xpass2_hg_8_levels, 8, 4);
hg_pass1levels6pfa6(xfft_hg_12K_ip_8_4_K8, xpass2_hg_8_levels, 8, 4);
//hg_pass1levels6pfa7(xfft_hg_14K_ip_8_4_K8, xpass2_hg_8_levels, 8, 4);
//hg_pass1levels6pfa8(xfft_hg_16K_ip_8_4_K8, xpass2_hg_8_levels, 8, 4);
//hg_pass1levels7pfa5(xfft_hg_20K_ip_8_4_K8, xpass2_hg_8_levels, 8, 4);
//hg_pass1levels7pfa6(xfft_hg_24K_ip_8_4_K8, xpass2_hg_8_levels, 8, 4);
//hg_pass1levels7pfa7(xfft_hg_28K_ip_8_4_K8, xpass2_hg_8_levels, 8, 4);
//hg_pass1levels7pfa8(xfft_hg_32K_ip_8_4_K8, xpass2_hg_8_levels, 8, 4);

hg_pass1levels6complex3(xfft_hg_12K_ac_ip_8_4_K8, xpass2_hg_8_levels, 8, 4);
//hg_pass1levels6complex4(xfft_hg_16K_ac_ip_8_4_K8, xpass2_hg_8_levels, 8, 4);
//hg_pass1levels7complex3(xfft_hg_24K_ac_ip_8_4_K8, xpass2_hg_8_levels, 8, 4);
//hg_pass1levels7complex4(xfft_hg_32K_ac_ip_8_4_K8, xpass2_hg_8_levels, 8, 4);


// xfft_hg_768K_12_2_K8
