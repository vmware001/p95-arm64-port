#pragma once
//#define __restrict

// double define to allow it's inputs to expand correctly.
#define CONCAT2(A, B) A##B
#define CONCAT(A, B) CONCAT2(A, B)

// vec2f64 in gwtables.h
//typedef double __attribute__ ((vector_size (16))) vec2f64;

#if 0 // const in C
__attribute__((visibility("hidden"))) extern const vec2f64 XMM_BIGVAL2;
__attribute__((visibility("hidden"))) extern const vec2f64 XMM_SQRTHALF2;
#else
static const vec2f64 XMM_BIGVAL2 = { 0x1.8p+52, 0x1.8p+52 };
static const vec2f64 XMM_SQRTHALF2 = {
0.7071067811865475244,
0.7071067811865475244 };
#endif
#define XMM_BIGVAL1 (XMM_BIGVAL2[0])
#define XMM_SQRTHALF1 (XMM_SQRTHALF2[0])
// g->u.xmm.XMM_BIGVAL

#define addpd(A,B) A=A+B
#define subpd(A,B) A=A-B
#define mulpd(A,B) A=A*B
#define divpd(A,B) A=A/B

#define addsd(A,B) A=A+B
#define subsd(A,B) A=A-B
#define mulsd(A,B) A=A*B
#define divsd(A,B) A=A/B

#define mov(A,B) A = B
#define movsd(A,B) A = B
#define xload(A,B) A = B

// need tmp, A and B can be the same.
//#define shufpd(A,B,I) {vec2f64 t = B; A[0] = A[I & 1]; A[1] = t[(I >> 1) & 1]}
#define shufpd(A,B,I) A = __builtin_shufflevector(A, B, I&1, (I>>1)&1|2);

#define unpcklpd(A,B) A[1] = (B)[0];
#define unpckhpd(A,B) A[0] = A[1]; A[1] = (B)[1];
#define unpcklo unpcklpd
#define unpckhi unpckhpd

#define low_load(A,C,D) { A[0] = C[0]; A[1] = D[0]; }
#define high_load(A,C,D) { A[0] = C[1]; A[1] = D[1]; }

#define shuffle_store(M1,M2,A,B) { \
	M1 = __builtin_shufflevector(A,B,0,2); \
	M2 = __builtin_shufflevector(A,B,1,3); \
}
#define shuffle_store_with_temp(M1,M2,A,B,C) shuffle_store(M1,M2,A,B)

#define shuffle_store_partial(M1,M2,A,B) \
	A = M1; \
	shuffle_store(M1,M2,A,B);

#define shuffle_load(A,B,M1,M2) { \
	A = __builtin_shufflevector(M1,M2,0,2); \
	B = __builtin_shufflevector(M1,M2,1,3); \
}

#define shuffle_load_with_temp(A,B,C,D,E) shuffle_load(A,B,C,D)

// x86  maxsd A,B is (A>B?A:B) ; when any input is NAN, then return B.
// arm64 fmax A,B is fmax(A,B) ; when any one input are NAN, then return other number.
// __builtin_elementwise_max is fast on ARM64, but slow multi-instruction on x86.

#if defined(__x86_64) || defined(__i386__)
	#define vec2max __builtin_ia32_maxpd
	#define vec2min __builtin_ia32_minpd
	static inline double vec2reducemax(vec2f64 a) {return a[0]>a[1]?a[0]:a[1];}
	static inline double vec2reducemin(vec2f64 a) {return a[0]<a[1]?a[0]:a[1];}
#elif defined(__has_builtin) && __has_builtin(__builtin_elementwise_max) \
		&& (!defined(__arm__) || (defined(__ARM_ARCH) && __ARM_ARCH >= 8))
	#define vec2reducemax __builtin_reduce_max
	#define vec2reducemin __builtin_reduce_min
	#define vec2max __builtin_elementwise_max
	#define vec2min __builtin_elementwise_min
#elif defined(__arm64__)
	#define vec2reducemax vmaxvq_f64
	#define vec2reducemin vminvq_f64
	#define vec2max vmaxq_f64
	#define vec2min vminq_f64
#else
	static inline vec2f64 vec2max(vec2f64 a, vec2f64 b)
		{vec2f64 r = {a[0]>b[0]?a[0]:b[0], a[1]>b[1]?a[1]:b[1]}; return r;}
	static inline vec2f64 vec2min(vec2f64 a, vec2f64 b)
		{vec2f64 r = {a[0]<b[0]?a[0]:b[0], a[1]<b[1]?a[1]:b[1]}; return r;}
	static inline double vec2reducemax(vec2f64 a) {return a[0]>a[1]?a[0]:a[1];}
	static inline double vec2reducemin(vec2f64 a) {return a[0]<a[1]?a[0]:a[1];}
#endif

#define maxpd(A,B) A = vec2max(A,B)
#define minpd(A,B) A = vec2min(A,B)

#define u8ptr(A) *(unsigned char*)(A)
#define u16ptr(A) *(unsigned short*)(A)
//#define u8ptr2(A,B) *((unsigned char*)(A) + (unsigned char*)(B))
#define xptr(A) (*(vec2f64*)(A))
#define f64ptr(A) (*(double*)(A))
#define loadu8(A,B) A = (*(unsigned char*)(B))
#define mov64(A,B) A = (uintptr_t)(B)

#define xptr2(A,B) (*(vec2f64*)((char*)(A) + (B)))
#define f64ptr2(A,B) (*(double*)((char*)(A) + (B)))

// Some CPU don't have fast roundeven.

#if defined(__arm__) && defined(__ARM_ARCH) && __ARM_ARCH >= 8 \
		&& defined(__clang_major__) && __clang_major__ <= 21
	// 32 bit, -march=armv8-a or -mcpu=cortex-a53
	// Old clang have roundeven compile error backend bug.
	#if 1
		#define roundf64 __builtin_elementwise_nearbyint
		#define roundx __builtin_elementwise_nearbyint
	#else
	static inline double roundf64(double a) {
		__asm__("vrintn.f64 %[out], %[in]" : [out]"=w"(a) : [in]"w"(a));
		return a;
	}
	static inline vec2f64 roundx(vec2f64 v) {
		v[0] = roundf64(v[0]);
		v[1] = roundf64(v[1]);
		return v;
	}
	#endif
#elif defined(__has_builtin) && __has_builtin(__builtin_elementwise_roundeven) \
		&& (!defined(__SSE2__) || defined(__SSE4_1__)) \
		&& (!defined(__arm__) || (defined(__ARM_ARCH) && __ARM_ARCH >= 8))
	#define roundf64 __builtin_elementwise_roundeven
	#define roundx __builtin_elementwise_roundeven
#elif defined(__SSE4_1__)
	/* PORTING NOTE: gwnum_c1 had a typo here ("roundpd(X,8)"); it only shows */
	/* up on compilers that lack __builtin_elementwise_roundeven (e.g. GCC),  */
	/* which is why it went unnoticed.  VROUNDPD imm8=8 => round-to-nearest-  */
	/* even using the immediate, suppress precision exception.                 */
	#define roundf64(A) __builtin_roundeven(A)
	#define roundx(A) __builtin_ia32_roundpd(A,8)
#elif defined(__arm64__)
	#define roundf64(A) __builtin_roundeven(A)
	#define roundx(A) vrndnq_f64(A)
#else
	#define roundf64(A) (((A) + XMM_BIGVAL1) - XMM_BIGVAL1)
	#define roundx(A) (((A) + XMM_BIGVAL2) - XMM_BIGVAL2)
#endif

#if defined(__has_builtin) && __has_builtin(__builtin_elementwise_abs)
	#define absf64 __builtin_elementwise_abs
	#define absx __builtin_elementwise_abs
#elif defined(__arm64__)
	#include <arm_neon.h>
	#define absf64 __builtin_fabs
	#define absx vabsq_f64
#else
	#include <math.h>
	#define absf64 fabs
	static inline vec2f64 absx(vec2f64 a) {a[0] = fabs(a[0]); a[1] = fabs(a[1]); return a;}
#endif


// Const numbers from gwdbldblc.cpp, gwasm_constants

static const double CONST_P951 = 0x1.e6f0e134454ffp-1;   // 0.95105651629515353
static const double CONST_P618 = 0x1.3c6ef372fe95p-1;    // 0.61803398874989490
static const double CONST_P309 = 0x1.3c6ef372fe95p-2;    // 0.30901699437494745
static const double CONST_M262 = -0x1.4f1bbcdcbfa53p+1;  // -2.61803398874989446
static const double CONST_P588 = 0x1.2cf2304755a5ep-1;   // 0.58778525229247314
static const double CONST_M162 = -0x1.9e3779b97f4a8p+0;  // -1.61803398874989490
static const double CONST_M809 = -0x1.9e3779b97f4a8p-1;  // -0.80901699437494745
static const double CONST_M382 = -0x1.8722191a02d61p-2;  // -0.38196601125010515
static const double CONST_P866 = 0x1.bb67ae8584caap-1;   // 0.86602540378443860
static const double CONST_P433 = 0x1.bb67ae8584caap-2;   // 0.43301270189221930
static const double CONST_P577 = 0x1.279a74590331cp-1;   // 0.57735026918962573
static const double CONST_P975 = 0x1.f329c0558e969p-1;   // 0.97492791218182362
static const double CONST_P445 = 0x1.c7b90e3024582p-2;   // 0.44504186791262879
//#define CONST_P434_P975 CONST_P445
static const double CONST_P180 = 0x1.cd4bca9cb5c71p+0;   // 1.80193773580483829
static const double CONST_P623 = 0x1.3f3a0e28bedd1p-1;   // 0.62348980185873348
static const double CONST_M358 = -0x1.6d761c42b2c41p-2;  // -0.35689586789220945
static const double CONST_P404 = 0x1.032176315d291p+2;   // 4.04891733952230570
static const double CONST_M223 = -0x1.c7b90e3024582p-3;  // -0.22252093395631439
static const double CONST_M901 = -0x1.cd4bca9cb5c71p-1;  // -0.90096886790241915
static const double CONST_M691 = -0x1.6250a3698fe64p-1;  // -0.69202147163009586
static const double CONST_P924 = 0x1.d906bcf328d46p-1;   // 0.92387953251128674
static const double CONST_P383 = 0x1.87de2a6aea963p-2;   // 0.38268343236508978
static const double CONST_P782 = 0x1.904c37505de4bp-1;   // 0.78183148246802980
static const double CONST_P434 = 0x1.bc4c04d71abc1p-2;   // 0.43388373911755812
static const double CONST_P975_P434 = 0x1.1f9d07145f6e9p+1;   // 2.24697960371746719
static const double CONST_P782_P434 = 0x1.cd4bca9cb5c71p+0;   // 1.80193773580483829
static const double CONST_P259 = 0x1.0907dc193069p-2;    // 0.25881904510252074
static const double CONST_P966 = 0x1.ee8dd4748bf15p-1;   // 0.96592582628906831
static const double CONST_P259_P707 = 0x1.76cf5d0b09954p-2;   // 0.36602540378443860
static const double CONST_P966_P707 = 0x1.5db3d742c2655p+0;   // 1.36602540378443860
static const double CONST_P924_P383 = 0x1.3504f333f9de6p+1;   // 2.41421356237309492
static const double CONST_P981_P195 = 0x1.41bfee2424771p+2;   // 5.02733949212584807
static const double CONST_P981 = 0x1.f6297cff75cbp-1;    // 0.98078528040323043
static const double CONST_P195 = 0x1.8f8b83c69a60bp-3;   // 0.19509032201612828
static const double CONST_P831_P556 = 0x1.7f218e25a7462p+0;   // 1.49660576266548917
static const double CONST_P831 = 0x1.a9b66290ea1a3p-1;   // 0.83146961230254524
static const double CONST_P556 = 0x1.1c73b39ae68c8p-1;   // 0.55557023301960218
static const double CONST_P556_P195 = 0x1.6c835e79946a3p+1;   // 2.84775906502257348
static const double CONST_P901_P975 = 0x1.d928bdedb0faap-1;   // 0.92413896109109328
static const double CONST_P623_P975 = 0x1.476fb0b30acecp-1;   // 0.63952400384496633
static const double CONST_P223_P975 = 0x1.d37150903fb55p-3;   // 0.22824347439014994
static const double CONST_P1_P975 = 0x1.0695615eb6194p+0;   // 1.02571686327255396
static const double CONST_P782_P975 = 0x1.9a9795396b8e2p-1;   // 0.80193773580483829
static const double CONST_P809 = -CONST_M809;
static const double CONST_P901 = -CONST_M901;
static const double CONST_P223 = -CONST_M223;

#define CONST2_P951 __extension__(vec2f64){CONST_P951, CONST_P951}
#define CONST2_P618 __extension__(vec2f64){CONST_P618, CONST_P618}
#define CONST2_P309 __extension__(vec2f64){CONST_P309, CONST_P309}
#define CONST2_M262 __extension__(vec2f64){CONST_M262, CONST_M262}
#define CONST2_P588 __extension__(vec2f64){CONST_P588, CONST_P588}
#define CONST2_M162 __extension__(vec2f64){CONST_M162, CONST_M162}
#define CONST2_M809 __extension__(vec2f64){CONST_M809, CONST_M809}
#define CONST2_M382 __extension__(vec2f64){CONST_M382, CONST_M382}
#define CONST2_P866 __extension__(vec2f64){CONST_P866, CONST_P866}
#define CONST2_P433 __extension__(vec2f64){CONST_P433, CONST_P433}
#define CONST2_P577 __extension__(vec2f64){CONST_P577, CONST_P577}
#define CONST2_P975 __extension__(vec2f64){CONST_P975, CONST_P975}
//#define CONST2_P434_P975 __extension__(vec2f64){CONST_P434_P975, CONST_P434_P975}
#define CONST2_P445 __extension__(vec2f64){CONST_P445, CONST_P445}
#define CONST2_P180 __extension__(vec2f64){CONST_P180, CONST_P180}
#define CONST2_P623 __extension__(vec2f64){CONST_P623, CONST_P623}
#define CONST2_P358 __extension__(vec2f64){CONST_P358, CONST_P358}
#define CONST2_M358 __extension__(vec2f64){CONST_M358, CONST_M358}
#define CONST2_P404 __extension__(vec2f64){CONST_P404, CONST_P404}
#define CONST2_M223 __extension__(vec2f64){CONST_M223, CONST_M223}
#define CONST2_M901 __extension__(vec2f64){CONST_M901, CONST_M901}
#define CONST2_M691 __extension__(vec2f64){CONST_M691, CONST_M691}
#define CONST2_P924 __extension__(vec2f64){CONST_P924, CONST_P924}
#define CONST2_P383 __extension__(vec2f64){CONST_P383, CONST_P383}
#define CONST2_P782 __extension__(vec2f64){CONST_P782, CONST_P782}
#define CONST2_P434 __extension__(vec2f64){CONST_P434, CONST_P434}
#define CONST2_P975_P434 __extension__(vec2f64){CONST_P975_P434, CONST_P975_P434}
#define CONST2_P782_P434 __extension__(vec2f64){CONST_P782_P434, CONST_P782_P434}
#define CONST2_P259 __extension__(vec2f64){CONST_P259, CONST_P259}
#define CONST2_P966 __extension__(vec2f64){CONST_P966, CONST_P966}
#define CONST2_P259_P707 __extension__(vec2f64){CONST_P259_P707, CONST_P259_P707}
#define CONST2_P966_P707 __extension__(vec2f64){CONST_P966_P707, CONST_P966_P707}
#define CONST2_P924_P383 __extension__(vec2f64){CONST_P924_P383, CONST_P924_P383}
#define CONST2_P981_P195 __extension__(vec2f64){CONST_P981_P195, CONST_P981_P195}
#define CONST2_P981 __extension__(vec2f64){CONST_P981, CONST_P981}
#define CONST2_P195 __extension__(vec2f64){CONST_P195, CONST_P195}
#define CONST2_P831_P556 __extension__(vec2f64){CONST_P831_P556, CONST_P831_P556}
#define CONST2_P831 __extension__(vec2f64){CONST_P831, CONST_P831}
#define CONST2_P556 __extension__(vec2f64){CONST_P556, CONST_P556}
#define CONST2_P556_P195 __extension__(vec2f64){CONST_P556_P195, CONST_P556_P195}
#define CONST2_P901_P975 __extension__(vec2f64){CONST_P901_P975, CONST_P901_P975}
#define CONST2_P623_P975 __extension__(vec2f64){CONST_P623_P975, CONST_P623_P975}
#define CONST2_P223_P975 __extension__(vec2f64){CONST_P223_P975, CONST_P223_P975}
#define CONST2_P1_P975 __extension__(vec2f64){CONST_P1_P975, CONST_P1_P975}
#define CONST2_P782_P975 __extension__(vec2f64){CONST_P782_P975, CONST_P782_P975}

#define CONST2_P809 __extension__(vec2f64){CONST_P809, CONST_P809}
#define CONST2_P901 __extension__(vec2f64){CONST_P901, CONST_P901}
#define CONST2_P223 __extension__(vec2f64){CONST_P223, CONST_P223}

