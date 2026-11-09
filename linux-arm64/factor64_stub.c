/* factor64_stub.c -- LINK-TEST ONLY weak stubs for the 64-bit trial-factoring
 * assembly entry points that upstream linux64 gets from factor64.asm
 * (prebuilt as linux64/factor64.o, x86-64 only).
 *
 * The real ARM64 implementation lives in factor64_c.c (portable C port of
 * factor64.asm, owned by the factor64-c work item).  This stub exists only so
 * that the mprime aarch64 binary can be linked and started before that file
 * lands.  The Makefile.arm64 picks factor64_c.c automatically once it exists
 * and never links this file then.
 *
 * The stubs are intentionally weak: if factor64_c.c is linked in as well, its
 * strong definitions win.
 *
 * WARNING: with this stub linked, trial factoring silently reports "no factor
 * found" for every range (factor64_sieve does nothing, factor64_tf finds
 * nothing).  Every LL/PRP/P-1 result is unaffected (they do not use this
 * code), but do NOT trust Factor=... worktodo results from a binary built
 * with this stub.  A loud warning is printed so such a binary is easy to spot.
 */

#include <stdio.h>

/* Keep the struct definitions in sync with commonb.c (struct facasm_data,
 * 64-bit version).  We only need the opaque pointer here; using the real
 * struct is not required for the stubs, so declare the parameter as void*
 * compatible with the prototype in commonb.c via the common header typedefs.
 * commonb.c declares:
 *	EXTERNC int factor64_pass_setup (struct facasm_data *);
 *	EXTERNC int factor64_small (struct facasm_data *);
 *	EXTERNC int factor64_sieve (struct facasm_data *);
 *	EXTERNC int factor64_tf (struct facasm_data *);
 * We use a forward declaration so no layout is duplicated here.
 */

struct facasm_data;

static void factor64_stub_warn (void)
{
	static int warned = 0;
	if (!warned) {
		warned = 1;
		fprintf (stderr, "WARNING: mprime built with factor64_stub.c "
				 "(no trial-factoring implementation linked); "
				 "Factor= results are NOT valid!\n");
	}
}

__attribute__((weak))
int factor64_pass_setup (struct facasm_data *asm_data)
{
	(void) asm_data;
	factor64_stub_warn ();
	return (0);
}

__attribute__((weak))
int factor64_small (struct facasm_data *asm_data)
{
	(void) asm_data;
	factor64_stub_warn ();
	return (0);		/* 0 = no more small factors found */
}

__attribute__((weak))
int factor64_sieve (struct facasm_data *asm_data)
{
	(void) asm_data;
	factor64_stub_warn ();
	/* The asm sieves a 12KB bit array; a zeroed "all candidate" sieve means
	   every chunk is TFed as a no-op.  commonb.c tracks progress via
	   facdist12K, not via the sieve contents, so this terminates. */
	return (0);
}

__attribute__((weak))
int factor64_tf (struct facasm_data *asm_data)
{
	(void) asm_data;
	factor64_stub_warn ();
	return (0);		/* 0 = no factor found in this sieve area */
}
