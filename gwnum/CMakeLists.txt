cmake_minimum_required(VERSION 3.10)

project(gwnum_c1_lib LANGUAGES C CXX)
#add_compile_options("-ffp-model=strict" "-")


add_library(gwnum_c1_lib
cpuid.c
radix.c
gwbench.c
gwdbldbl.cpp
gwini.c
gwthread.cpp
gwtables.c
gwutil.c
gwnum.c
giants.c
#main.c
#main1.s

multtable2.c
mult_gwinfo.c
gianthlp.c
xarch.c


r4/r4_p4.c
r4/r4_core.c
r4/r4_k8.c
r4/r4_k10.c
r4/r4dwpn_core1.c
r4/r4dwpn_core2.c
r4/r4dwpn_core3.c
r4/r4dwpn_core4.c
r4/r4dwpn_core5.c
r4/r4dwpnac_core1.c
r4/r4dwpnac_core2.c
r4/r4dwpn_k10_1.c
r4/r4dwpn_k10_2.c
r4/r4dwpn_k10_3.c
r4/r4dwpnac_k10_1.c
r4/r4dwpn_k8_1.c
r4/r4dwpn_k8_2.c
r4/r4dwpnac_k8_1.c
r4/r4dwpn_p4_1.c
r4/r4dwpn_p4_2.c
r4/r4dwpn_p4_3.c
r4/r4dwpnac_p4_1.c
r4/r4pass2_8.c
r4/r4pass2_768.c
r4/r4pass2_10.c
r4/r4pass2_1280.c
r4/r4pass2_1536.c
r4/r4pass2_11.c
r4/r4pass2_2304.c
r4/r4pass2_2560.c
r4/r4pass2_3072.c
r4/r4pass2_3840.c
r4/r4pass2_12.c
r4/r4pass2_4608.c
r4/r4pass2_5120.c
r4/r4pass2_6144.c
r4/r4pass2_6400.c
r4/r4pass2_7680.c
r4/r4pass2_13.c
r4/r4pass2_9216.c
r4/r4pass2_10240.c
r4/r4pass2_12288.c
r4/r4pass2_12800.c
r4/r4pass2_15360.c
r4/r4pass2_14.c
#r4/r4pass2_20480.c
r4/r4pass2_25600.c


hg/hg_blend.c
hg/hg_core.c
hg/hg_p4.c
hg/hg_k8.c
hg/hg_k10.c
hg/hgonepass.c
hg/hgpass2_8.c
#hg/hgpass2_10.c
hg/hgpass2_11.c
hg/hgpass2_12.c
#hg/hgpass2_13.c

x/xmult1a.c
x/xmult2.c
x/xmult2a.c
x/xmult3.c
x/xmult3a.c

x/xmult1ax.c
x/xmult2ax.c
x/xmult3ax.c
x/xnormal.c
)


if(CMAKE_SYSTEM_PROCESSOR MATCHES "x86_64|amd64|i386|i686|x86")
target_sources(gwnum_c1_lib PRIVATE
x/xmult2a_k8.c
x/xmult1a_sse41.c
x/xmult2a_sse41.c
x/xmult2a_k8_sse41.c
x/xmult3a_k8.c
x/xmult3a_sse41.c
x/xmult3a_k8_sse41.c
cpuid_xsave.c
)

set_source_files_properties(cpuid_xsave.c PROPERTIES COMPILE_FLAGS "-mxsave")
set_source_files_properties(xmult1ax.c PROPERTIES COMPILE_FLAGS "-march=x86-64-v2 -fmacro-backtrace-limit=0")
set_source_files_properties(xmult1a_sse41.c PROPERTIES COMPILE_FLAGS "-msse4.1")
set_source_files_properties(xmult2a_sse41.c PROPERTIES COMPILE_FLAGS "-msse4.1")
set_source_files_properties(xmult3a_sse41.c PROPERTIES COMPILE_FLAGS "-msse4.1")
set_source_files_properties(xmult2a_k8_sse41.c PROPERTIES COMPILE_FLAGS "-msse4.1")
set_source_files_properties(xmult3a_k8_sse41.c PROPERTIES COMPILE_FLAGS "-msse4.1")
#set_source_files_properties(main.c PROPERTIES COMPILE_FLAGS "-msse4.1")
#set_source_files_properties(gianthlp.c PROPERTIES COMPILE_FLAGS "-Os")
endif()

# -fno-fast-math -frounding-math -fsignaling-nans -ffp-contract=off
# -mcmodel=small or large
# -fdirect-access-external-data

# -ffp-contract=fast have a chance of round problem due to extra precision.

#target_compile_options(gwnum_c1_lib PRIVATE
#"-ffp-contract=on"
#"-ftrapping-math"
#"-frounding-math"
#"-fno-fast-math"
# "-fno-pie"
#)

if (CMAKE_CXX_COMPILER_ID MATCHES "Clang")
	#message(STATUS "Using Clang compiler")
	if(CMAKE_SYSTEM_PROCESSOR MATCHES "x86_64|amd64|i386|i686|x86|arm64|aarch64")
		#target_compile_options(gwnum_c1_lib PRIVATE "-ffp-exception-behavior=maytrap")
	endif()
else()
	# Avoid gcc default fp-contract=fast with FMA (ARM64), possible round errors.
	include(CheckCXXCompilerFlag)
	check_cxx_compiler_flag("-ffp-contract=on" HAVE_FP_CONTRACT)
	if(HAVE_FP_CONTRACT)
		target_compile_options(gwnum_c1_lib PRIVATE -ffp-contract=on)
	endif()
endif()

#target_link_options(gwnum_c1_lib PRIVATE "-no-pie")
#target_link_libraries(gwnum_c1_lib PRIVATE /mnt/2/1/code/gwnum/gwnum/linux64/gwnum.a)


install(TARGETS gwnum_c1_lib RUNTIME DESTINATION bin)
target_include_directories(gwnum_c1_lib PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/qd)
