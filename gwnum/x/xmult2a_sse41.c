#include <stdlib.h>
#include <ctype.h>

#include "../gwnum.h"
#include "../gwtables.h"
#include "../gwnum_func.h"
#include "../gwcommon.h"
#include "../unravel.h"
#include "../xarch.h"
#include "../xbasics.h"
#include "xnormal.h"

#include "xmult2a.h"

inorm(xr2s4, noexec, noexec, noexec, exec, exec);
inorm(xr2es4, noexec, exec, noexec, exec, exec);
inorm(xr2cs4, noexec, noexec, exec, exec, exec);
inorm(xr2ecs4, noexec, exec, exec, exec, exec);
inorm(xi2s4, exec, noexec, noexec, exec, exec);
inorm(xi2es4, exec, exec, noexec, exec, exec);
inorm(xi2cs4, exec, noexec, exec, exec, exec);
inorm(xi2ecs4, exec, exec, exec, exec, exec);
zpnorm(xr2zps4, noexec, noexec, noexec, exec, exec, exec, noexec, noexec);
zpnorm(xr2zps4c1, noexec, noexec, noexec, exec, exec, exec, exec, noexec);
zpnorm(xr2zps4cm1, noexec, noexec, noexec, exec, exec, exec, noexec, exec);
zpnorm(xr2zpes4, noexec, exec, noexec, exec, exec, exec, noexec, noexec);
zpnorm(xr2zpes4c1, noexec, exec, noexec, exec, exec, exec, exec, noexec);
zpnorm(xr2zpes4cm1, noexec, exec, noexec, exec, exec, exec, noexec, exec);
zpnorm(xr2zpcs4, noexec, noexec, exec, exec, exec, exec, noexec, noexec);
zpnorm(xr2zpecs4, noexec, exec, exec, exec, exec, exec, noexec, noexec);
zpnorm(xi2zps4, exec, noexec, noexec, exec, exec, exec, noexec, noexec);
zpnorm(xi2zps4c1, exec, noexec, noexec, exec, exec, exec, exec, noexec);
zpnorm(xi2zps4cm1, exec, noexec, noexec, exec, exec, exec, noexec, exec);
zpnorm(xi2zpes4, exec, exec, noexec, exec, exec, exec, noexec, noexec);
zpnorm(xi2zpes4c1, exec, exec, noexec, exec, exec, exec, exec, noexec);
zpnorm(xi2zpes4cm1, exec, exec, noexec, exec, exec, exec, noexec, exec);
zpnorm(xi2zpcs4, exec, noexec, exec, exec, exec, exec, noexec, noexec);
zpnorm(xi2zpecs4, exec, exec, exec, exec, exec, exec, noexec, noexec);
zpnorm(xr2zps4k, noexec, noexec, noexec, exec, exec, noexec, noexec, noexec);
zpnorm(xr2zps4kc1, noexec, noexec, noexec, exec, exec, noexec, exec, noexec);
zpnorm(xr2zps4kcm1, noexec, noexec, noexec, exec, exec, noexec, noexec, exec);
zpnorm(xr2zpes4k, noexec, exec, noexec, exec, exec, noexec, noexec, noexec);
zpnorm(xr2zpes4kc1, noexec, exec, noexec, exec, exec, noexec, exec, noexec);
zpnorm(xr2zpes4kcm1, noexec, exec, noexec, exec, exec, noexec, noexec, exec);
zpnorm(xr2zpcs4k, noexec, noexec, exec, exec, exec, noexec, noexec, noexec);
zpnorm(xr2zpecs4k, noexec, exec, exec, exec, exec, noexec, noexec, noexec);
zpnorm(xi2zps4k, exec, noexec, noexec, exec, exec, noexec, noexec, noexec);
zpnorm(xi2zps4kc1, exec, noexec, noexec, exec, exec, noexec, exec, noexec);
zpnorm(xi2zps4kcm1, exec, noexec, noexec, exec, exec, noexec, noexec, exec);
zpnorm(xi2zpes4k, exec, exec, noexec, exec, exec, noexec, noexec, noexec);
zpnorm(xi2zpes4kc1, exec, exec, noexec, exec, exec, noexec, exec, noexec);
zpnorm(xi2zpes4kcm1, exec, exec, noexec, exec, exec, noexec, noexec, exec);
zpnorm(xi2zpcs4k, exec, noexec, exec, exec, exec, noexec, noexec, noexec);
zpnorm(xi2zpecs4k, exec, exec, exec, exec, exec, noexec, noexec, noexec);

inorm(xr2bs4, noexec, noexec, noexec, noexec, exec);
inorm(xr2ebs4, noexec, exec, noexec, noexec, exec);
inorm(xr2cbs4, noexec, noexec, exec, noexec, exec);
inorm(xr2ecbs4, noexec, exec, exec, noexec, exec);
inorm(xi2bs4, exec, noexec, noexec, noexec, exec);
inorm(xi2ebs4, exec, exec, noexec, noexec, exec);
inorm(xi2cbs4, exec, noexec, exec, noexec, exec);
inorm(xi2ecbs4, exec, exec, exec, noexec, exec);
zpnorm(xr2zpbs4, noexec, noexec, noexec, noexec, exec, exec, noexec, noexec);
zpnorm(xr2zpbs4c1, noexec, noexec, noexec, noexec, exec, exec, exec, noexec);
zpnorm(xr2zpbs4cm1, noexec, noexec, noexec, noexec, exec, exec, noexec, exec);
zpnorm(xr2zpebs4, noexec, exec, noexec, noexec, exec, exec, noexec, noexec);
zpnorm(xr2zpebs4c1, noexec, exec, noexec, noexec, exec, exec, exec, noexec);
zpnorm(xr2zpebs4cm1, noexec, exec, noexec, noexec, exec, exec, noexec, exec);
zpnorm(xr2zpcbs4, noexec, noexec, exec, noexec, exec, exec, noexec, noexec);
zpnorm(xr2zpecbs4, noexec, exec, exec, noexec, exec, exec, noexec, noexec);
zpnorm(xi2zpbs4, exec, noexec, noexec, noexec, exec, exec, noexec, noexec);
zpnorm(xi2zpbs4c1, exec, noexec, noexec, noexec, exec, exec, exec, noexec);
zpnorm(xi2zpbs4cm1, exec, noexec, noexec, noexec, exec, exec, noexec, exec);
zpnorm(xi2zpebs4, exec, exec, noexec, noexec, exec, exec, noexec, noexec);
zpnorm(xi2zpebs4c1, exec, exec, noexec, noexec, exec, exec, exec, noexec);
zpnorm(xi2zpebs4cm1, exec, exec, noexec, noexec, exec, exec, noexec, exec);
zpnorm(xi2zpcbs4, exec, noexec, exec, noexec, exec, exec, noexec, noexec);
zpnorm(xi2zpecbs4, exec, exec, exec, noexec, exec, exec, noexec, noexec);
zpnorm(xr2zpbs4k, noexec, noexec, noexec, noexec, exec, noexec, noexec, noexec);
zpnorm(xr2zpbs4kc1, noexec, noexec, noexec, noexec, exec, noexec, exec, noexec);
zpnorm(xr2zpbs4kcm1, noexec, noexec, noexec, noexec, exec, noexec, noexec, exec);
zpnorm(xr2zpebs4k, noexec, exec, noexec, noexec, exec, noexec, noexec, noexec);
zpnorm(xr2zpebs4kc1, noexec, exec, noexec, noexec, exec, noexec, exec, noexec);
zpnorm(xr2zpebs4kcm1, noexec, exec, noexec, noexec, exec, noexec, noexec, exec);
zpnorm(xr2zpcbs4k, noexec, noexec, exec, noexec, exec, noexec, noexec, noexec);
zpnorm(xr2zpecbs4k, noexec, exec, exec, noexec, exec, noexec, noexec, noexec);
zpnorm(xi2zpbs4k, exec, noexec, noexec, noexec, exec, noexec, noexec, noexec);
zpnorm(xi2zpbs4kc1, exec, noexec, noexec, noexec, exec, noexec, exec, noexec);
zpnorm(xi2zpbs4kcm1, exec, noexec, noexec, noexec, exec, noexec, noexec, exec);
zpnorm(xi2zpebs4k, exec, exec, noexec, noexec, exec, noexec, noexec, noexec);
zpnorm(xi2zpebs4kc1, exec, exec, noexec, noexec, exec, noexec, exec, noexec);
zpnorm(xi2zpebs4kcm1, exec, exec, noexec, noexec, exec, noexec, noexec, exec);
zpnorm(xi2zpcbs4k, exec, noexec, exec, noexec, exec, noexec, noexec, noexec);
zpnorm(xi2zpecbs4k, exec, exec, exec, noexec, exec, noexec, noexec, noexec);
