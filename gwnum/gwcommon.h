/*----------------------------------------------------------------------
| gwcommon.h
|
| This file contains handy #defines that I use in all my projects
| 
|  Copyright 2005-2021 Mersenne Research, Inc.  All rights reserved.
+---------------------------------------------------------------------*/

#ifndef _GWCOMMON_H
#define _GWCOMMON_H

#ifndef TRUE
#define TRUE	1
#endif
#ifndef FALSE
#define FALSE	0
#endif

/* PORTING (ARM64 Linux): gwnum_c1 only listed __arm64__ (Apple clang).  Linux */
/* aarch64 compilers (clang/gcc target aarch64-linux-*) define __aarch64__, so  */
/* without this the build would select the old 32-bit/x87 code paths (NEW_STACK_SIZE, */
/* FFT(1) emulation, and other layout/behaviour).  We only need the "64-bit mode"  */
/* semantics of this macro, not actual x86 code (all SSE2/SSE4.1 kernels are C).   */
#if defined(__x86_64__) || defined(_M_X64) || defined(__arm64__) || defined(__aarch64__) || defined(_M_ARM64)
#define X86_64
#endif

/* In many cases where the C code is interfacing with the assembly code */
/* we must declare variables that are exactly 32 or 64 bits wide.  This is the */
/* portable way to do this.  Finally supported in Visual Studio 2013. */

#include "inttypes.h"
#ifndef __cplusplus
#include "stdbool.h"
#endif

/* Define the ASSERT macro I use while debugging */

#include <assert.h>
#ifdef GDEBUG
#define ASSERTG		assert
#else
#define ASSERTG(a)
#endif
#define GWASSERT  	assert

#endif
