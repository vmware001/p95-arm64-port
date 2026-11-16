# contrib/validation — proving an optimised build still computes correctly

Any change to the gwnum FFT kernels (compiler flags, hand-written assembled
kernels, data-layout changes, ...) must pass **all three layers** below before
it is allowed anywhere near the farm. Speed numbers alone mean nothing: a
broken FFT kernel still produces plausible-looking ECM "no factor" results
while silently wasting weeks of computation.

The current production build (`-O3 -ffp-contract=on -mcpu=cortex_a53`, the
"mc4" binary) passes all three layers.

## Layer 1 — bit-exact arithmetic vs GMP (`gwnum_exact_strict.c`)

gwnum's FFT multiply is **exact by design** (the FFT length is chosen so the
roundoff is far below the carry threshold). Therefore for any inputs a,b the
recovered integer must equal the true value *exactly* - a single differing bit
means the optimisation is broken.

The test performs, for every iteration, four independent checks and compares
each result limb-by-limb against GMP:

| check | operation | expectation |
|---|---|---|
| `mul`    | `gwmul3(a,b)`            | `a*b mod 2^n-1` |
| `add`    | `gwadd3o(a,b)`           | `a+b mod 2^n-1` |
| `sqr`    | `gwmul3(a,a)`            | `a*a mod 2^n-1` |
| `fftmul` | both inputs pre-`gwfft`ed | `a*b mod 2^n-1` |

Inputs cycle through random values, deliberately correlated worst cases
(`a = b = 2^n-2`, `a = b = 1`) and values near the modulus, because that is
where FFT roundoff is largest. The PRNG is seeded deterministically so runs
are reproducible.

Build (host):

```
zig cc -target aarch64-linux-musl -O2 -I<gwnum dir> -I<_deps>/arm64/include \
        -c gwnum_exact_strict.c -o exact_strict.o
zig cc -target aarch64-linux-musl -static -no-pie -s -o exact_strict \
        exact_strict.o libgwnum.a libpolymult.a libgmp.a -lm -lpthread
```

Run (on the box):

```
./exact_strict 1277 400      # small: 64-double FFT
./exact_strict 11213 100     # 512
./exact_strict 264169 15     # production size: 16K FFT (the farm's hot path)
```

Expected: every line ends in `ALL_EXACT`, exit code 0.

## Layer 2 — end-to-end PRP against known historical answers (`mprime_known_results.sh`)

Runs the *complete* client (FFT + carries + Gerbicz check + result reporting)
on small exponents whose answers have been public for decades:

| exponent | expected answer | note |
|---|---|---|
| M1277  | composite `"status":"C"` | smallest Mersenne not fully factored |
| M11213 | prime `"status":"P"` | Mersenne prime, 1963 |
| M19937 | prime `"status":"P"` | Mersenne prime, 1971 |

It creates an isolated scratch dir (`/data/local/tmp/val`), copies the binary
under test, disables PrimeNet, runs the three tests and prints the status
counts. Expected: exactly `1 x C` and `2 x P`, no mprime error/roundoff lines.
The log also shows mprime's own Gerbicz error checks passing.

```
sh mprime_known_results.sh [/path/to/mprime]
```

## Layer 3 — reproduce the existing production results

Two complementary checks:

1. **Digest equality between the old and the new build** on the real workload
   sizes (16K FFT for M264169, 400K for M7831391): run the same number of
   full-size multiplies with a deterministic input and compare a hash over the
   result gwnum's raw FFT doubles / limbs. Any difference means the new build
   is not doing the same arithmetic. (`_work/asm/fftbench_asm.c` from the
   assembly experiment implements this; a digest must be stable across repeated
   runs of the same binary before it can be trusted.)
2. **Re-run a completed ECM curve** from `results.json.txt` with the new build
   and confirm the same outcome is produced (stage 1+2 complete, no factor, no
   FFT-error retries and no "roundoff > 0.40" warnings in run.log).

## Practical notes

- Correctness tests are cheap: all of layer 1 takes < 3 s on the box, layer 2
  a few seconds, so run them after *every* rebuild, not just before deployment.
- `giantstruct` is `{int maxsize; int sign; uint32_t *n;}` — the limbs live at
  `g->n`, **not** at an offset of 8 bytes into the struct. Writing to
  `(char*)g + 8` corrupts the `n` pointer and segfaults (this was a real bug in
  the first version of the harness).
- Timing measurements must be taken with the farm paused
  (`kill -STOP $(cat mprime.pid)`, then `kill -CONT` **in all cases**), and
  `taskset` on toybox Android takes a bitmask (`taskset 4` = CPU2), not `-c`.
