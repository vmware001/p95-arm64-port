// edtest.cpp - independent check of gwnum Edwards point arithmetic.
//
// Copies of production ed_dbl/ed_add/ed_extend/ed_check from ecm.cpp run
// against random inputs at a given modulus.  All inputs and outputs are
// printed as hex; verify.py recomputes everything with plain big-int math.
// If the C output disagrees with Python, the gwnum/ECM implementation is
// wrong at this FFT size - no Prime95 logic involved.
//
// Build (aarch64, links only gwnum + gmp):
//   zig c++ -target aarch64-linux-musl -I<gwnum> -O2 -std=c++11 -c edtest.cpp
//   zig c++ -target aarch64-linux-musl -static -s -o edtest edtest.o gwnum.a libgmp.a -lm -lpthread

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include <signal.h>
#include <ucontext.h>
#include <unistd.h>

static void segv_handler (int sig, siginfo_t *si, void *uc0) {
	ucontext_t *uc = (ucontext_t *) uc0;
	uint64_t *pc = (uint64_t *) uc->uc_mcontext.pc;
	uint64_t *fp = (uint64_t *) uc->uc_mcontext.regs[29];
	uint64_t *lr = (uint64_t *) uc->uc_mcontext.regs[30];
	char buf[512];
	int n = snprintf (buf, sizeof buf, "SEGV addr=%p pc=%p lr=%p sp=%p fp=%p\n",
			  si->si_addr, pc, lr, uc->uc_mcontext.sp, fp);
	write (2, buf, n);
	for (int i = 0; i < 16 && fp; i++) {
		uint64_t *ra = (uint64_t *) fp[1];
		uint64_t *nfp = (uint64_t *) fp[0];
		if (!ra) break;
		n = snprintf (buf, sizeof buf, "  frame %d: fp=%p ret=%p\n", i, fp, ra);
		write (2, buf, n);
		if (nfp <= fp) break;
		fp = nfp;
	}
	_exit (139);
}

static void install_segv (void) {
	struct sigaction sa;
	memset (&sa, 0, sizeof sa);
	sa.sa_sigaction = segv_handler;
	sa.sa_flags = SA_SIGINFO;
	sigaction (SIGSEGV, &sa, NULL);
}

#include "gwnum.h"
#include "gwcommon.h"

/* ---- stubs for prime95 services the copied code touches ---- */
#define INI_FILE "prime.txt"
extern "C" int IniGetInt (const char *file, const char *key, int dflt);

/* ---- options for ed_add/ed_dbl (from ecm.cpp) ---- */
#define	ED_XYZ			0x1
#define ED_RESULT_FOR_ADD	0x2
#define ED_RESULT_FOR_DBL	0x4
#define ED_STARTNEXTFFT		0x8
#define ED_FFT_S1		0x10
#define ED_FFT_S1Z		0x20
#define ED_FFT_S1T		0x40
#define ED_SUBTRACT		0x80
#define ED_FORCE_EXTEND		0x100
#define ED_NO_GWALLOC		0x200

/* ---- minimal ecmhandle (copied functions only touch these fields) ---- */
struct ed {
	gwnum	x, y, z, t;
	bool	alternate_format;
};
struct ecmhandle {
	gwhandle gwdata;
	gwnum	ed_a;		/* NULL = untwisted */
	gwnum	ed_half;
	gwnum	ed_d;		/* curve parameter (only ed_check uses it) */
	giant	N;
};

/* ================= copied from ecm.cpp (production code, unmodified math) ================= */

/* Create the value FFT(1/2) */
static void create_ed_half (ecmhandle *ecmdata) {
	ecmdata->ed_half = gwalloc (&ecmdata->gwdata);
	gshiftright (1, ecmdata->N);
	gianttogw (&ecmdata->gwdata, ecmdata->N, ecmdata->ed_half);
	gwsmalladd (&ecmdata->gwdata, 1, ecmdata->ed_half);
	gwfft (&ecmdata->gwdata, ecmdata->ed_half, ecmdata->ed_half);
	gshiftleft (1, ecmdata->N);
	iaddg (1, ecmdata->N);
}

static void ed_extend (ecmhandle *ecmdata, struct ed *e, int options) {
	int mul_options = (options & ED_STARTNEXTFFT) ? GWMUL_STARTNEXTFFT : 0;
	ASSERTG (!e->alternate_format);
	if (!(options & ED_FORCE_EXTEND) && e->t != NULL) return;
	if (!(options & ED_FORCE_EXTEND)) e->t = gwalloc (&ecmdata->gwdata);
	gwmul3 (&ecmdata->gwdata, e->x, e->y, e->t, GWMUL_FFT_S12 | mul_options);	// t = xy
	if (e->z != NULL) {
		gwmul3 (&ecmdata->gwdata, e->x, e->z, e->x, GWMUL_FFT_S2 | mul_options);	// x = xz
		gwmul3 (&ecmdata->gwdata, e->y, e->z, e->y, GWMUL_FFT_S2 | mul_options);	// y = yz
		gwsquare2 (&ecmdata->gwdata, e->z, e->z, mul_options);			// z = z^2
	}
}

static void ed_add (ecmhandle *ecmdata, struct ed *in1, struct ed *in2, struct ed *out, int options) {
	gwnum	tmp1, tmp2;
	bool	safe11 = mul_safe (&ecmdata->gwdata, 1, 1);
	int	subtract = options & ED_SUBTRACT;
	int	extended = options & ED_RESULT_FOR_ADD;
	int	mul_S1Z_options = (options & ED_FFT_S1Z) ? GWMUL_FFT_S1 : 0;
	int	mul_S1T_options = (options & ED_FFT_S1T) ? GWMUL_FFT_S1 : 0;
	int	mul_options = (options & ED_STARTNEXTFFT) ? GWMUL_STARTNEXTFFT : 0;

	ASSERTG (!in1->alternate_format && !in2->alternate_format);
	ASSERTG (ecmdata->ed_a == NULL || FFT_state (ecmdata->ed_a) == FULLY_FFTed);
	ASSERTG (!(options & ED_NO_GWALLOC) || (in1 != out && in1->t != NULL && out->t != NULL));

	if (in1->t == NULL) ed_extend (ecmdata, in1, ED_STARTNEXTFFT);
	if (in2->t == NULL) ed_extend (ecmdata, in2, ED_STARTNEXTFFT);

	gwnum in1z = in1->z;
	gwnum in2z = in2->z;
	if (out->z == NULL) out->z = gwalloc (&ecmdata->gwdata);
	if (out->t != NULL) tmp1 = out->t;
	else tmp1 = gwalloc (&ecmdata->gwdata);
	if (options & ED_NO_GWALLOC) tmp2 = in1->t;
	else tmp2 = gwalloc (&ecmdata->gwdata);

	gwnum D;
	if (in2z == NULL) D = in1->t;
	else gwmul3 (&ecmdata->gwdata, in1->t, in2z, D = tmp2, mul_S1T_options | GWMUL_FFT_S2 | GWMUL_STARTNEXTFFT_IF(safe11 || !extended)); // D = T1 * Z2
	gwnum C;
	if (in1z == NULL) C = in2->t;
	else gwmul3 (&ecmdata->gwdata, in1z, in2->t, C = out->z, mul_S1Z_options | GWMUL_FFT_S2 | GWMUL_STARTNEXTFFT_IF(safe11 || !extended)); // C = Z1 * T2
	gwaddsub4o (&ecmdata->gwdata, D, C, tmp2, out->z, GWADD_FORCE_NORMALIZE);	/* ARM64: delayed-norm H breaks the 16K multiply - force normalize */
        if (subtract) {
		gwswap (tmp2, out->z);
		gwmulmuladd5 (&ecmdata->gwdata, in1->x, in2->y, in1->y, in2->x, tmp1, GWMUL_FFT_S1234 | GWMUL_STARTNEXTFFT);	// F
		if (ecmdata->ed_a == NULL)
			gwmulmulsub5 (&ecmdata->gwdata, in1->y, in2->y, in1->x, in2->x, out->y, GWMUL_STARTNEXTFFT);	// G
		else {
			gwmul3 (&ecmdata->gwdata, in1->x, in2->x, out->x, GWMUL_STARTNEXTFFT);
			gwmulmulsub5 (&ecmdata->gwdata, in1->y, in2->y, ecmdata->ed_a, out->x, out->y, GWMUL_STARTNEXTFFT);
		}
        } else {
		gwmulmulsub5 (&ecmdata->gwdata, in1->x, in2->y, in1->y, in2->x, tmp1, GWMUL_FFT_S1234 | GWMUL_STARTNEXTFFT);	// F = X1*Y2 - Y1*X2
		if (ecmdata->ed_a == NULL)
			gwmulmuladd5 (&ecmdata->gwdata, in1->y, in2->y, in1->x, in2->x, out->y, GWMUL_STARTNEXTFFT);	// G
		else {
			gwmul3 (&ecmdata->gwdata, in1->x, in2->x, out->x, GWMUL_FFT_S1 | GWMUL_STARTNEXTFFT);
			gwmulmuladd5 (&ecmdata->gwdata, in1->y, in2->y, ecmdata->ed_a, out->x, out->y, GWMUL_STARTNEXTFFT);
		}
	}
	if (!extended) {
		gwmul3 (&ecmdata->gwdata, tmp2, tmp1, tmp2, GWMUL_FFT_S2 | mul_options);			// X3 = E * F
		gwmul3 (&ecmdata->gwdata, tmp1, out->y, tmp1, GWMUL_FFT_S2 | mul_options);			// Z3 = F * G
		gwmul3 (&ecmdata->gwdata, out->y, out->z, out->z, mul_options);					// Y3 = G * H
		gwswap (out->z, out->y);
		gwswap (tmp2, out->x);
		gwswap (tmp1, out->z);
		out->t = NULL;
	} else {
		gwmul3 (&ecmdata->gwdata, tmp2, tmp1, out->x, GWMUL_FFT_S12 | mul_options);			// X3 = E * F
		gwmul3 (&ecmdata->gwdata, tmp1, out->y, tmp1, GWMUL_FFT_S2 | mul_options);			// Z3 = F * G
		gwmul3 (&ecmdata->gwdata, out->y, out->z, out->y, GWMUL_FFT_S2 | mul_options);			// Y3 = G * H
		gwmul3 (&ecmdata->gwdata, tmp2, out->z, tmp2, mul_options);					// T3 = E * H
		gwswap (tmp1, out->z);
		out->t = tmp2, tmp2 = NULL;
	}

	if (options & ED_NO_GWALLOC) {
		if (!extended) out->t = tmp2;
		in1->t = tmp1;
	} else {
		gwfree (&ecmdata->gwdata, tmp1);
		gwfree (&ecmdata->gwdata, tmp2);
	}
}

static void ed_dbl (ecmhandle *ecmdata, struct ed *in, struct ed *out, int options) {
	gwnum	tmp1, tmp2;
	int	extended = options & ED_RESULT_FOR_ADD;
	int	mul_S1_options = (options & ED_FFT_S1) ? GWMUL_FFT_S1 : 0;
	int	mul_options = (options & ED_STARTNEXTFFT) ? GWMUL_STARTNEXTFFT : 0;
	int	can_alternate_format = !(options & ED_XYZ) && !extended;

	ASSERTG (ecmdata->ed_a == NULL || FFT_state (ecmdata->ed_a) == FULLY_FFTed);
	ASSERTG (!(options & ED_NO_GWALLOC) || (in != out && in->t != NULL && out->t != NULL && extended));

	static	int algo5 = 2;
	static	int algo5a = 2;
	static	int algo7 = 2;
	if (algo5 == 2) {
		algo5 = IniGetInt (INI_FILE, "EdDblAlgo5", 1);
		algo5a = IniGetInt (INI_FILE, "EdDblAlgo5a", 1);
		algo7 = IniGetInt (INI_FILE, "EdDblAlgo7", 1);
	}

	gwnum inz = in->z;
	if (out->z == NULL) out->z = gwalloc (&ecmdata->gwdata);
	if (out->t != NULL) tmp1 = out->t, out->t = NULL;
	else tmp1 = gwalloc (&ecmdata->gwdata);
	if (options & ED_NO_GWALLOC) tmp2 = in->t;
	else tmp2 = gwalloc (&ecmdata->gwdata);
	gwsetmulbyconst (&ecmdata->gwdata, 2);

	if (algo5 && squaresquareadd_safe (&ecmdata->gwdata, 0, 0, 0, 0)) {
	    if (inz == NULL) dbltogw (&ecmdata->gwdata, 2.0, out->z);
	    else gwsquare2 (&ecmdata->gwdata, inz, out->z, mul_S1_options | GWMUL_MULBYCONST | GWMUL_STARTNEXTFFT);	// C = 2Z^2
	    if (ecmdata->ed_a == NULL) {
		gwsquare2 (&ecmdata->gwdata, in->x, tmp2, GWMUL_FFT_S1 | GWMUL_MULBYCONST | GWMUL_STARTNEXTFFT);	// G-H = 2X^2
		gwmul3 (&ecmdata->gwdata, in->x, in->y, tmp1, GWMUL_FFT_S2 | GWMUL_MULBYCONST | GWMUL_STARTNEXTFFT);	// E = 2XY
		gwmulmuladd5 (&ecmdata->gwdata, in->y, in->y, in->x, in->x, out->y, GWMUL_STARTNEXTFFT);		// G = Y^2 + X^2
	    } else {
		gwsquare2 (&ecmdata->gwdata, in->x, tmp2, GWMUL_FFT_S1 | GWMUL_STARTNEXTFFT);
		gwmul3 (&ecmdata->gwdata, in->x, in->y, out->x, GWMUL_FFT_S2 | GWMUL_MULBYCONST | GWMUL_STARTNEXTFFT);
		gwmulmuladd5 (&ecmdata->gwdata, ecmdata->ed_a, tmp2, in->y, in->y, out->y, GWMUL_FFT_S2 | GWMUL_STARTNEXTFFT);
		gwmul3 (&ecmdata->gwdata, ecmdata->ed_a, tmp2, tmp2, GWMUL_MULBYCONST | GWMUL_STARTNEXTFFT);
		gwswap (tmp1, out->x);
	    }
	    if (!extended) {
		gwsubmul4 (&ecmdata->gwdata, out->y, tmp2, out->y, tmp2, GWMUL_FFT_S13 | mul_options);		// Y3 = H * G
		gwsubmul4 (&ecmdata->gwdata, out->z, out->y, out->y, out->x, GWMUL_FFT_S1 | mul_options);		// Z3 = F * G
		gwsubmul4 (&ecmdata->gwdata, out->z, out->y, tmp1, tmp1, mul_options);				// X3 = F * E
		gwswap (tmp2, out->y);
		gwswap (out->x, out->z);
		gwswap (tmp1, out->x);
	    } else {
		gwsubmul4 (&ecmdata->gwdata, out->y, tmp2, out->y, out->x, GWMUL_FFT_S123 | mul_options);	// Y3 = H * G
		gwsubmul4 (&ecmdata->gwdata, out->y, tmp2, tmp1, tmp2, GWMUL_FFT_S3 | mul_options);		// T3 = H * E
		gwsubmul4 (&ecmdata->gwdata, out->z, out->y, tmp1, tmp1, GWMUL_FFT_S1 | mul_options);		// X3 = F * E
		gwsubmul4 (&ecmdata->gwdata, out->z, out->y, out->y, out->z, mul_options);			// Z3 = F * G
		gwswap (out->x, out->y);
		gwswap (tmp2, out->t);
		gwswap (tmp1, out->x);
	    }
	}

	else if (algo5a && squaremuladd_safe (&ecmdata->gwdata, 0, 0, 0, 0)) {
	    if (ecmdata->ed_half == NULL) create_ed_half (ecmdata);
	    if (inz == NULL) dbltogw (&ecmdata->gwdata, 2.0, out->z);
	    else gwsquare2 (&ecmdata->gwdata, inz, out->z, mul_S1_options | GWMUL_MULBYCONST | GWMUL_STARTNEXTFFT);	// C = 2Z^2
	    if (ecmdata->ed_a == NULL) {
		gwsquare2 (&ecmdata->gwdata, in->x, tmp2, GWMUL_FFT_S1 | GWMUL_MULBYCONST | GWMUL_STARTNEXTFFT);	// G-H = 2X^2
		gwmul3 (&ecmdata->gwdata, in->x, in->y, tmp1, GWMUL_FFT_S2 | GWMUL_MULBYCONST | GWMUL_STARTNEXTFFT);	// E = 2XY
		gwmulmuladd5 (&ecmdata->gwdata, ecmdata->ed_half, tmp2, in->y, in->y, out->y, GWMUL_STARTNEXTFFT);	// G = Y^2 + 1/2 * 2X^2
	    } else {
		gwsquare2 (&ecmdata->gwdata, in->x, tmp2, GWMUL_FFT_S1 | GWMUL_STARTNEXTFFT);
		gwmul3 (&ecmdata->gwdata, in->x, in->y, out->x, GWMUL_FFT_S2 | GWMUL_MULBYCONST | GWMUL_STARTNEXTFFT);
		gwmul3 (&ecmdata->gwdata, ecmdata->ed_a, tmp2, tmp2, GWMUL_MULBYCONST | GWMUL_STARTNEXTFFT);
		gwmulmuladd5 (&ecmdata->gwdata, ecmdata->ed_half, tmp2, in->y, in->y, out->y, GWMUL_FFT_S2 | GWMUL_STARTNEXTFFT);
		gwswap (tmp1, out->x);
	    }
	    if (!extended) {
		gwsubmul4 (&ecmdata->gwdata, out->y, tmp2, out->y, tmp2, GWMUL_FFT_S13 | mul_options);
		gwsubmul4 (&ecmdata->gwdata, out->z, out->y, out->y, out->x, GWMUL_FFT_S1 | mul_options);
		gwsubmul4 (&ecmdata->gwdata, out->z, out->y, tmp1, tmp1, mul_options);
		gwswap (tmp2, out->y);
		gwswap (out->x, out->z);
		gwswap (tmp1, out->x);
	    } else {
		gwsubmul4 (&ecmdata->gwdata, out->y, tmp2, out->y, out->x, GWMUL_FFT_S123 | mul_options);
		gwsubmul4 (&ecmdata->gwdata, out->y, tmp2, tmp1, tmp2, GWMUL_FFT_S3 | mul_options);
		gwsubmul4 (&ecmdata->gwdata, out->z, out->y, tmp1, tmp1, GWMUL_FFT_S1 | mul_options);
		gwsubmul4 (&ecmdata->gwdata, out->z, out->y, out->y, out->z, mul_options);
		gwswap (out->x, out->y);
		gwswap (tmp2, out->t);
		gwswap (tmp1, out->x);
	    }
	}

	else if (algo7 && squareadd_safe (&ecmdata->gwdata, 0, 0, 0) && ecmdata->ed_a == NULL && inz != NULL && (in->alternate_format || can_alternate_format)) {
	    if (ecmdata->ed_half == NULL) create_ed_half (ecmdata);

	    if (!in->alternate_format) {
		gwsquare2 (&ecmdata->gwdata, inz, out->z, mul_S1_options | GWMUL_MULBYCONST | GWMUL_STARTNEXTFFT);	// C = 2Z^2
		gwsquare2 (&ecmdata->gwdata, in->x, tmp2, GWMUL_FFT_S1);						// A = X^2
		gwmul3 (&ecmdata->gwdata, in->x, in->y, out->x, GWMUL_FFT_S2 | GWMUL_MULBYCONST | GWMUL_STARTNEXTFFT);	// E = 2XY
		gwsquare2 (&ecmdata->gwdata, in->y, out->y, 0);							// B = Y^2
		gwaddsub4o (&ecmdata->gwdata, out->y, tmp2, tmp2, out->y, GWADD_FORCE_NORMALIZE);			// G, H
		gwsubmul4 (&ecmdata->gwdata, out->z, tmp2, tmp2, tmp1, GWMUL_FFT_S12 | GWMUL_MULBYCONST | GWMUL_STARTNEXTFFT);	// 2Z3 = 2 * F * G
		gwsubmul4 (&ecmdata->gwdata, out->z, tmp2, out->x, out->x, GWMUL_MULBYCONST | GWMUL_STARTNEXTFFT);	// 2X3 = 2 * F * E
		gwmulmuladd5 (&ecmdata->gwdata, tmp2, out->y, ecmdata->ed_half, out->x, out->y, GWMUL_FFT_S4 | GWMUL_STARTNEXTFFT); // Y3+X3
		gwswap (tmp1, out->z);
		out->alternate_format = TRUE;
	    }

	    else {
		gwsquare2 (&ecmdata->gwdata, in->y, tmp2, GWMUL_FFT_S1 | GWMUL_MULBYCONST | GWMUL_STARTNEXTFFT);	// 2M = 2(Y+X)^2
		gwsubmul4 (&ecmdata->gwdata, in->y, in->x, in->y, out->y, GWMUL_FFT_S1 | GWMUL_MULBYCONST | GWMUL_STARTNEXTFFT); // 2H
		gwmuladd4 (&ecmdata->gwdata, in->x, in->x, out->y, out->x, GWMUL_FFT_S3 | GWMUL_STARTNEXTFFT);		// 2G
		if (can_alternate_format) {
		    gwmulsub4 (&ecmdata->gwdata, inz, inz, out->x, out->z, GWMUL_FFT_S3 | GWMUL_MULBYCONST | GWMUL_STARTNEXTFFT); // 4F
		    gwmul3 (&ecmdata->gwdata, out->z, out->x, tmp1, GWMUL_FFT_S1 | GWMUL_STARTNEXTFFT);			// 8Z3
		    gwsubmul4 (&ecmdata->gwdata, tmp2, out->x, out->z, tmp2, GWMUL_STARTNEXTFFT);				// 8X3
		    gwmulmuladd5 (&ecmdata->gwdata, out->y, out->x, ecmdata->ed_half, tmp2, out->y, GWMUL_FFT_S4 | GWMUL_STARTNEXTFFT); // 4(Y3+X3)
		    gwswap (tmp1, out->z);
		    gwswap (tmp2, out->x);
		} else {
		    gwmulsub4 (&ecmdata->gwdata, inz, inz, out->x, out->z, GWMUL_FFT_S3 | GWMUL_STARTNEXTFFT);		// 2F
		    gwmul3 (&ecmdata->gwdata, out->z, out->x, tmp1, GWMUL_FFT_S1 | mul_options);				// 4Z3
		    gwsubmul4 (&ecmdata->gwdata, tmp2, out->x, out->z, out->z, GWMUL_FFT_S1 | mul_options);			// 4X3
		    if (extended) gwsubmul4 (&ecmdata->gwdata, tmp2, out->x, out->y, tmp2, GWMUL_FFT_S3 | mul_options);	// 4T3
		    gwmul3 (&ecmdata->gwdata, out->y, out->x, out->y, mul_options);						// 4Y3
		    gwswap (out->z, out->x);
		    gwswap (tmp1, out->z);
		    if (extended) gwswap (tmp2, out->t);
		    out->alternate_format = FALSE;
		}
	    }
	}

	else {													// Algorithm 4
	    bool safe21 = mul_safe (&ecmdata->gwdata, 2, 1);
	    if (inz == NULL) dbltogw (&ecmdata->gwdata, 2.0, out->z);
	    else gwsquare2 (&ecmdata->gwdata, inz, out->z, mul_S1_options | GWMUL_MULBYCONST | GWMUL_STARTNEXTFFT);	// C = 2Z^2
	    if (ecmdata->ed_a == NULL) {
		gwsquare2 (&ecmdata->gwdata, in->x, tmp2, GWMUL_FFT_S1 | GWMUL_STARTNEXTFFT_IF(safe21));		// A = X^2
	    } else {
		gwsquare2 (&ecmdata->gwdata, in->x, tmp2, GWMUL_FFT_S1 | GWMUL_STARTNEXTFFT);
		gwmul3 (&ecmdata->gwdata, ecmdata->ed_a, tmp2, tmp2, GWMUL_STARTNEXTFFT_IF(safe21));
	    }
	    gwmul3 (&ecmdata->gwdata, in->x, in->y, out->x, GWMUL_FFT_S2 | GWMUL_MULBYCONST | GWMUL_STARTNEXTFFT);	// E = 2XY
	    gwsquare2 (&ecmdata->gwdata, in->y, out->y, GWMUL_STARTNEXTFFT_IF(safe21));					// B = Y^2
	    gwaddsub4o (&ecmdata->gwdata, out->y, tmp2, tmp2, out->y, GWADD_DELAYNORM_IF(safe21));			// G, H
	    if (!extended) {
		gwmul3 (&ecmdata->gwdata, tmp2, out->y, out->y, GWMUL_FFT_S1 | mul_options);				// Y3 = G * H
		gwsubmul4 (&ecmdata->gwdata, out->z, tmp2, tmp2, tmp1, GWMUL_FFT_S1 | mul_options);			// Z3 = F * G
		gwsubmul4 (&ecmdata->gwdata, out->z, tmp2, out->x, out->x, mul_options);				// X3 = F * E
		gwswap (tmp1, out->z);
	    } else {
		gwmul3 (&ecmdata->gwdata, tmp2, out->y, tmp1, GWMUL_FFT_S12 | mul_options);				// Y3 = G * H
		gwmul3 (&ecmdata->gwdata, out->x, out->y, out->y, GWMUL_FFT_S1 | mul_options);				// T3 = E * H
		gwsubmul4 (&ecmdata->gwdata, out->z, tmp2, out->x, out->x, GWMUL_FFT_S1 | mul_options);			// X3 = F * E
		gwsubmul4 (&ecmdata->gwdata, out->z, tmp2, tmp2, out->z, mul_options);					// Z3 = F * G
		gwswap (out->y, out->t);
		gwswap (tmp1, out->y);
	    }
	}

	if (options & ED_NO_GWALLOC)
		in->t = (tmp1 != NULL ? tmp1 : tmp2);
	else {
		gwfree (&ecmdata->gwdata, tmp1);
		gwfree (&ecmdata->gwdata, tmp2);
	}
}

static bool ed_check (ecmhandle *ecmdata, struct ed *in) {
	gwnum	t0, t1, t2, t3;
	bool	result;

	ASSERTG (!in->alternate_format);

	t0 = gwalloc (&ecmdata->gwdata);
	t1 = gwalloc (&ecmdata->gwdata);
	t2 = gwalloc (&ecmdata->gwdata);
	t3 = gwalloc (&ecmdata->gwdata);

	gwsquare2 (&ecmdata->gwdata, in->x, t0, GWMUL_PRESERVE_S1);
	gwsquare2 (&ecmdata->gwdata, in->y, t1, GWMUL_PRESERVE_S1);
	if (in->z != NULL) gwsquare2 (&ecmdata->gwdata, in->z, t2, GWMUL_PRESERVE_S1);
	else dbltogw (&ecmdata->gwdata, 1.0, t2);
	gwmul3 (&ecmdata->gwdata, t0, t1, t3, 0);
	gwmul3 (&ecmdata->gwdata, t3, ecmdata->ed_d, t3, 0);
	if (ecmdata->ed_a != NULL) gwmul3 (&ecmdata->gwdata, t0, ecmdata->ed_a, t0, 0);
	gwadd3 (&ecmdata->gwdata, t0, t1, t0);
	gwmul3 (&ecmdata->gwdata, t0, t2, t0, 0);
	gwsquare2 (&ecmdata->gwdata, t2, t2, 0);
	gwadd3 (&ecmdata->gwdata, t2, t3, t2);
	gwsub3 (&ecmdata->gwdata, t0, t2, t0);
	giant tmp = popg (&ecmdata->gwdata.gdata, ((unsigned long) ecmdata->gwdata.bit_length >> 5) + 5);
	gwtogiant (&ecmdata->gwdata, t0, tmp);
	modg (ecmdata->N, tmp);
	result = (tmp->sign == 0);
	pushg (&ecmdata->gwdata.gdata, 1);

	gwfree (&ecmdata->gwdata, t0);
	gwfree (&ecmdata->gwdata, t1);
	gwfree (&ecmdata->gwdata, t2);
	gwfree (&ecmdata->gwdata, t3);

	return (result);
}

/* ================= test harness ================= */

static uint64_t lcg_state = 0x123456789abcdef0ULL;
static uint32_t lcg (void) { lcg_state = lcg_state * 6364136223846793005ULL + 1442695040888963407ULL; return (uint32_t)(lcg_state >> 33); }

static void print_giant (const char *tag, giant g) {
	int words = (g->sign < 0) ? -g->sign : g->sign;
	printf ("%s ", tag);
	for (int i = words - 1; i >= 0; i--) printf ("%08x", g->n[i]);
	printf ("\n");
}

static gwhandle *gwp;	/* for helpers */
static void print_gwnum (const char *tag, gwnum a) {
	giant g = popg (&gwp->gdata, ((int)gwp->bit_length >> 5) + 5);
	gwtogiant (gwp, a, g);
	print_giant (tag, g);
	pushg (&gwp->gdata, 1);
}

/* deterministic random residue < N (words limbs, top masked) */
static void rand_giant (giant g, int words, uint32_t topmask) {
	g->sign = words;
	for (int i = 0; i < words; i++) g->n[i] = lcg ();
	g->n[words-1] &= topmask;
	if (g->n[words-1] == 0) g->n[words-1] = 1;
}

static gwnum rand_gwnum (int words, uint32_t topmask) {
	giant g = popg (&gwp->gdata, words + 2);
	rand_giant (g, words, topmask);
	gwnum a = gwalloc (gwp);
	gianttogw (gwp, g, a);
	pushg (&gwp->gdata, 1);
	return a;
}

int main (int argc, char **argv) {
	install_segv ();
	unsigned long ns[2];
	int nn = 0;
	for (int i = 1; i < argc; i++) ns[nn++] = atol (argv[i]);
	if (nn == 0) { ns[0] = 264169; ns[1] = 1277; nn = 2; }

	for (int t = 0; t < nn; t++) {
		unsigned long n = ns[t];
		int words = (int)(n >> 5) + 1;
		uint32_t topmask = (n & 31) ? ((1u << (n & 31)) - 1) : 0xffffffffu;

		struct ecmhandle ed;
		memset (&ed, 0, sizeof (ed));
		gwp = &ed.gwdata;
		gwinit (&ed.gwdata);
		gwset_maxmulbyconst (&ed.gwdata, 2);	/* ed_dbl multiplies by 2 */
		gwset_bench_cores (&ed.gwdata, 1);	/* as in production ecm setup */
		gwset_bench_workers (&ed.gwdata, 1);
		gwset_num_threads (&ed.gwdata, 1);
		gwset_safety_margin (&ed.gwdata, 0.0);
		gwset_larger_fftlen_count (&ed.gwdata, 0);
		gwset_minimum_fftlen (&ed.gwdata, 0);
		if (gwsetup (&ed.gwdata, 1.0, 2.0, (double) n, -1.0)) { printf ("gwsetup failed\n"); return 1; }

		/* N = 2^n - 1 */
		giant Ng = popg (&ed.gwdata.gdata, words + 1);
		Ng->sign = (int)((n + 31) / 32);
		for (int i = 0; i < Ng->sign; i++) Ng->n[i] = 0xffffffff;
		Ng->n[Ng->sign-1] = topmask;
		ed.N = Ng;
		/* curve parameter d (random nonzero, used only by ed_check) */
		ed.ed_a = NULL;
		ed.ed_half = NULL;
		ed.ed_d = rand_gwnum (words, topmask);

		printf ("=== n=%lu fft=%d ===\n", n, (int) gwfftlen (&ed.gwdata));

		/* basic op sanity: c = a*b mod N via plain gwmul3 */
		gwnum a = rand_gwnum (words, topmask);
		gwnum b = rand_gwnum (words, topmask);
		gwnum c = gwalloc (&ed.gwdata);
		gwmul3 (&ed.gwdata, a, b, c, 0);
		print_gwnum ("BASIC_a", a); print_gwnum ("BASIC_b", b); print_gwnum ("BASIC_ab", c);

		/* isolated primitive probes at this FFT size */
		{
			gwnum s1 = rand_gwnum (words, topmask), s2 = rand_gwnum (words, topmask);
			gwnum s3 = rand_gwnum (words, topmask), s4 = rand_gwnum (words, topmask);
			gwnum d1 = gwalloc (&ed.gwdata), d2 = gwalloc (&ed.gwdata);
			gwmulmuladd5 (&ed.gwdata, s1, s2, s3, s4, d1, 0);
			print_gwnum ("MMA5_a", s1); print_gwnum ("MMA5_b", s2); print_gwnum ("MMA5_c", s3);
			print_gwnum ("MMA5_d", s4); print_gwnum ("MMA5_abcd", d1);
			gwmulmuladd5 (&ed.gwdata, s1, s2, s3, s4, d2, GWMUL_STARTNEXTFFT);
			print_gwnum ("MMA5S_abcd", d2);
			gwnum e1 = gwalloc (&ed.gwdata);
			gwsubmul4 (&ed.gwdata, s1, s2, s3, e1, 0);
			print_gwnum ("SM4_amb", e1);	/* (s1-s2)*s3 */
			gwnum f1 = gwalloc (&ed.gwdata), f2 = gwalloc (&ed.gwdata);
			gwaddsub4o (&ed.gwdata, s1, s2, f1, f2, GWADD_DELAYNORM_IF (0));
			print_gwnum ("AS4_add", f1); print_gwnum ("AS4_sub", f2);
		}

		/* ---- ed_dbl, projective output with XYZ (save/stage-end path) ---- */
		{
			struct ed in, out; memset (&in, 0, sizeof in); memset (&out, 0, sizeof out);
			in.x = rand_gwnum (words, topmask); in.y = rand_gwnum (words, topmask); in.z = rand_gwnum (words, topmask);
			out.x = gwalloc (&ed.gwdata); out.y = gwalloc (&ed.gwdata);	/* out coords must pre-exist like the accumulator */
			print_gwnum ("DBL_XYZ_x", in.x); print_gwnum ("DBL_XYZ_y", in.y); print_gwnum ("DBL_XYZ_z", in.z);
			ed_dbl (&ed, &in, &out, ED_RESULT_FOR_DBL | ED_XYZ);
			print_gwnum ("DBL_XYZ_X3", out.x); print_gwnum ("DBL_XYZ_Y3", out.y); print_gwnum ("DBL_XYZ_Z3", out.z);
		}

		/* ---- ed_dbl, extended output + STARTNEXTFFT (mid-loop prep path) ---- */
		{
			struct ed in, out; memset (&in, 0, sizeof in); memset (&out, 0, sizeof out);
			in.x = rand_gwnum (words, topmask); in.y = rand_gwnum (words, topmask); in.z = rand_gwnum (words, topmask);
			out.x = gwalloc (&ed.gwdata); out.y = gwalloc (&ed.gwdata); out.t = gwalloc (&ed.gwdata);
			print_gwnum ("DBLF_x", in.x); print_gwnum ("DBLF_y", in.y); print_gwnum ("DBLF_z", in.z);
			ed_dbl (&ed, &in, &out, ED_RESULT_FOR_ADD | ED_STARTNEXTFFT);
			print_gwnum ("DBLF_X3", out.x); print_gwnum ("DBLF_Y3", out.y); print_gwnum ("DBLF_Z3", out.z); print_gwnum ("DBLF_T3", out.t);
		}

		/* ---- ed_add on two extended points, XYZ output (save-path add) ---- */
		{
			struct ed p1, p2, out;
			memset (&p1, 0, sizeof p1); memset (&p2, 0, sizeof p2); memset (&out, 0, sizeof out);
			p1.x = rand_gwnum (words, topmask); p1.y = rand_gwnum (words, topmask); p1.z = rand_gwnum (words, topmask);
			p2.x = rand_gwnum (words, topmask); p2.y = rand_gwnum (words, topmask); p2.z = rand_gwnum (words, topmask);
			ed_extend (&ed, &p1, 0); ed_extend (&ed, &p2, 0);
			out.x = gwalloc (&ed.gwdata); out.y = gwalloc (&ed.gwdata);
			print_gwnum ("ADDX_x1", p1.x); print_gwnum ("ADDX_y1", p1.y); print_gwnum ("ADDX_t1", p1.t); print_gwnum ("ADDX_z1", p1.z);
			print_gwnum ("ADDX_x2", p2.x); print_gwnum ("ADDX_y2", p2.y); print_gwnum ("ADDX_t2", p2.t); print_gwnum ("ADDX_z2", p2.z);
			ed_add (&ed, &p1, &p2, &out, ED_RESULT_FOR_DBL | ED_XYZ);
			print_gwnum ("ADDX_X3", out.x); print_gwnum ("ADDX_Y3", out.y); print_gwnum ("ADDX_Z3", out.z);
		}

		/* ---- same add with STARTNEXTFFT (mid-loop add) ---- */
		{
			struct ed p1, p2, out;
			memset (&p1, 0, sizeof p1); memset (&p2, 0, sizeof p2); memset (&out, 0, sizeof out);
			p1.x = rand_gwnum (words, topmask); p1.y = rand_gwnum (words, topmask); p1.z = rand_gwnum (words, topmask);
			p2.x = rand_gwnum (words, topmask); p2.y = rand_gwnum (words, topmask); p2.z = rand_gwnum (words, topmask);
			ed_extend (&ed, &p1, 0); ed_extend (&ed, &p2, 0);
			out.x = gwalloc (&ed.gwdata); out.y = gwalloc (&ed.gwdata); out.t = gwalloc (&ed.gwdata);
			print_gwnum ("ADDS_x1", p1.x); print_gwnum ("ADDS_y1", p1.y); print_gwnum ("ADDS_t1", p1.t); print_gwnum ("ADDS_z1", p1.z);
			print_gwnum ("ADDS_x2", p2.x); print_gwnum ("ADDS_y2", p2.y); print_gwnum ("ADDS_t2", p2.t); print_gwnum ("ADDS_z2", p2.z);
			ed_add (&ed, &p1, &p2, &out, ED_RESULT_FOR_DBL | ED_STARTNEXTFFT);
			print_gwnum ("ADDS_X3", out.x); print_gwnum ("ADDS_Y3", out.y); print_gwnum ("ADDS_Z3", out.z);
		}

		gwdone (&ed.gwdata);
		fflush (stdout);
	}
	return 0;
}
