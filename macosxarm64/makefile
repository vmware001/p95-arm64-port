# ===========================================================================
# mprime for macOS on Apple Silicon (arm64)
#
# Native build (no Rosetta).  The FFT/multiply/normalization kernels are the
# gwnum 31.7 "assembly converted to C" ARM64 port (128-bit vectors emulated
# with arm_neon.h) - the same code used by the Linux arm64 port in this tree.
# -DX86_64 is kept from the upstream macOS makefile: in the shared sources it
# only means "64-bit build" (64-bit factoring path, 64-bit benchmark strings),
# it is not x86-specific.  -DPORT=10 is the upstream MacOSX 64-bit port id.
#
# Requirements (installed on the Mac, e.g. with Homebrew):
#   brew install gmp hwloc boost curl
# Xcode command line tools provide clang, libc++, IOKit and CoreFoundation.
#
# Build:
#   make -C macosxarm64            # builds gwnum.a/polymult.a first, then mprime
#   ./macosxarm64/build-macos-arm64.sh   # same, with checks and a smoke test
#
# Notes:
#   * GCC->clang: the gwnum C port uses GCC vector extensions, which Apple
#     clang supports; arm_neon.h is provided by clang.
#   * APPLE_CPU tunes the scheduler/ISA level (default apple-m1; try apple-m2,
#     apple-m3, apple-m4, or APPLE_CPU= for a generic armv8 build).
#   * factor64 is the portable C replacement for the x86 factor64.asm
#     (factor64_c.c in this directory).
# ===========================================================================

ENVP = MACOSX_DEPLOYMENT_TARGET=11.0

APPLE_CPU ?= apple-m1

BREW   ?= brew
BREW_PREFIX := $(shell $(BREW) --prefix 2>/dev/null)
GMP_PREFIX   ?= $(shell $(BREW) --prefix gmp 2>/dev/null)
HWLOC_PREFIX ?= $(shell $(BREW) --prefix hwloc 2>/dev/null)
CURL_PREFIX  ?= $(shell $(BREW) --prefix curl 2>/dev/null)
BOOST_PREFIX ?= $(shell $(BREW) --prefix boost 2>/dev/null)

ARCHFLAGS = -arch arm64 -mmacosx-version-min=11.0
ifneq ($(APPLE_CPU),)
ARCHFLAGS += -mcpu=$(APPLE_CPU)
endif

INCS = -I.. -I../gwnum -I../linux \
       $(if $(GMP_PREFIX),-I$(GMP_PREFIX)/include) \
       $(if $(HWLOC_PREFIX),-I$(HWLOC_PREFIX)/include) \
       $(if $(CURL_PREFIX),-I$(CURL_PREFIX)/include) \
       $(if $(BOOST_PREFIX),-I$(BOOST_PREFIX)/include)

CFLAGS   = $(ARCHFLAGS) $(INCS) -std=gnu99 -O3 -ffp-contract=on -MMD -MP -DX86_64 -DCOMMAND_LINE_MPRIME -DPORT=10 -Wno-unused-result
CPPFLAGS = $(ARCHFLAGS) $(INCS) -std=c++11 -O3 -ffp-contract=on -MMD -MP -DX86_64 -DCOMMAND_LINE_MPRIME -DPORT=10 -Wno-unused-result

LFLAGS = $(ARCHFLAGS) -Wl,-no_pie

GMP_LIB   = $(if $(GMP_PREFIX),-L$(GMP_PREFIX)/lib) -lgmp
HWLOC_LIB = $(if $(HWLOC_PREFIX),-L$(HWLOC_PREFIX)/lib) -lhwloc
CURL_LIB  = $(if $(CURL_PREFIX),-L$(CURL_PREFIX)/lib) -lcurl

GWNUM_LIBDIR = $(CURDIR)
LIBS = $(GWNUM_LIBDIR)/gwnum.a $(GWNUM_LIBDIR)/polymult.a \
       $(HWLOC_LIB) $(GMP_LIB) $(CURL_LIB) \
       -lm -lpthread -framework IOKit -framework CoreFoundation -lc++

CC  = clang
CPP = clang++

FACTOR_SRC = $(if $(wildcard factor64_c.c),factor64_c.c,factor64_stub.c)
FACTOROBJ  = $(FACTOR_SRC:.c=.o)

OBJS = prime.o menu.o cJSON.o ecm.o exponentiate.o pair.o pm1prob.o
EXE  = mprime

##############################################################################

all: $(EXE)

# The two gwnum archives are built from the shared ARM64 C port.  "force"
# makes the sub-make run on every build (it is cheap when nothing changed);
# without it the archives would never be rebuilt once they exist.
$(GWNUM_LIBDIR)/gwnum.a $(GWNUM_LIBDIR)/polymult.a: force
	$(MAKE) -C ../gwnum -f Makefile.arm64 ARCH=macos APPLE_CPU="$(APPLE_CPU)" \
		OBJDIR="$(CURDIR)/obj-gwnum" LIBDIR="$(CURDIR)"

force:

$(EXE): security $(OBJS) $(FACTOROBJ) $(GWNUM_LIBDIR)/gwnum.a $(GWNUM_LIBDIR)/polymult.a
	$(ENVP) $(CPP) $(LFLAGS) -o $(EXE) $(OBJS) $(FACTOROBJ) $(LIBS)

# Upstream ships empty "trial version" stubs for the closed-source security
# code; the shared makefiles create them with touch.  Do the same, portably.
security:
	-touch ../security.h ../security.c ../secure5.c

cJSON.o: ../cJSON.c
	$(ENVP) $(CC) $(CFLAGS) -c ../cJSON.c -o cJSON.o

ecm.o: ../ecm.cpp
	$(ENVP) $(CPP) $(CPPFLAGS) -c ../ecm.cpp -o ecm.o

exponentiate.o: ../exponentiate.c
	$(ENVP) $(CC) $(CFLAGS) -c ../exponentiate.c -o exponentiate.o

pair.o: ../pair.cpp
	$(ENVP) $(CPP) $(CPPFLAGS) -c ../pair.cpp -o pair.o

pm1prob.o: ../pm1prob.c
	$(ENVP) $(CC) $(CFLAGS) -c ../pm1prob.c -o pm1prob.o

%.o: %.c
	$(ENVP) $(CC) $(CFLAGS) -c $< -o $@

clean:
	-rm -f $(EXE) $(OBJS) $(FACTOROBJ) $(OBJS:.o=.d) $(FACTOROBJ:.o=.d)
	-rm -rf obj-gwnum gwnum.a polymult.a

distclean: clean
	-rm -rf obj-gwnum

# Automatic header/source dependencies (-MMD -MP).  prime.c #includes
# commona.c/commonb.c/commonc.c and friends, so tracking them is essential:
# without it editing ../commonc.c only relinks while prime.o stays stale.
-include $(OBJS:.o=.d) $(FACTOROBJ:.o=.d)

.PHONY: all clean distclean security force
