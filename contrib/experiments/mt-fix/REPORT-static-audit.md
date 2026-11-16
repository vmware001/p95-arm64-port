# MT corruption in gwnum 31.7 ARM64 (asm→C port) — static audit

Task: `task-5`. Scope: READ-ONLY audit of `p95v3107b02.arm64.source/gwnum/` (git root
`D:\data\code\p95\p95v3107b02.arm64.source`) to explain: 8-thread-per-worker ECM stage 1 on
M264169/M291143 (16K FFT, Pass1=64, Pass2=256, clm=4) intermittently fails the end-of-stage-1
`ed_check` ("Point no longer on the Edwards curve"), while 1-thread (exact), 4-thread, and
8-thread PRP-squares runs are clean.

Method: (a) exhaustive scan of the converted kernel trees for shared/static state and for reads of
MT-scoped `asm_data` fields; (b) field-by-field cross-check of every `g->FIELD` read in compiled
converted sources against `auxiliary_thread`'s copy list + `internal_gwsetup` static init; (c)
line-by-line comparison of the converted normalisation/carry macros against their `.mac` originals;
(d) architecture review of the gwnum MT state machine (`gwnum.c`) for x86-TSO assumptions.

## Verdict (short)

No defect was found in the *converted x86→C kernels* that is thread-unsafe or MT-only wrong; the
converted code is a faithful translation and touches no shared state. The strongest, evidence-backed
candidate is in the **upstream MT framework, not the conversion**: the pass-1 carry-section state
machine (`gwnum/gwnum.c`) synchronises cross-thread fields with plain loads/stores and relies on
**x86-TSO store order** for correctness. That ordering is not guaranteed on ARM64 and the port only
added acquire/release to `gwthread.cpp` atomics, not to these plain fields. The failure mode matches
the observation exactly: intermittent, MT-only, increases with thread count, and produces a
*numerically clean* (roundoff-invisible) wrong product.

## Ranked suspects

### S1 (HIGH) — ARM64 memory-ordering hole in the pass-1 carry-section handshake
`gwnum/gwnum.c:6313-6319` (writer), `gwnum/gwnum.c:6366-6377` (reader), `gwnum/gwnum.c:5684-5696`
(section split writer), `gwnum/gwnum.h:923-931` (plain `int` fields).

Sequence per state-1 (inverse-FFT+normalise) pass of every multiply:
1. Thread *i* normalises the first `num_postfft_blocks` blocks of its section (writing FFT data),
   then does a **plain store** `sections[i].carry_in_blocks_finished = TRUE` and signals an event:
   `gwnum.c:6316-6319`.
2. Thread *j* (owner of the previous section) must not add its carries into those blocks until they
   are normalised. Its fast-path test is two **plain loads**:
   `while (last_block == sections[carry_out].start_block && !sections[carry_out].carry_in_blocks_finished)`
   — `gwnum.c:6368-6369`.
3. Section splitting re-uses a thread's section slot and relies on store order: it writes
   `start_block/next_block/last_block/carry_out_section`, then `carry_in_blocks_finished = FALSE`,
   with the explicit comment `// Must be done AFTER changing start_block!` (`gwnum.c:5685-5689`) and
   the invariant comment at `gwnum.c:6344-6345` ("SECTION SPLITTING CODE MUST CHANGE THESE VALUES IN
   A PARTICULAR ORDER SUCH THAT THIS CHECK ALWAYS WORKS").

Why it is safe on x86 and broken here: on x86-64 stores from one CPU become visible to other CPUs in
program order, so (a) if *j* sees the flag TRUE it also sees the preceding FFT-data writes, and (b) if
*j* sees the new `start_block` it also sees `carry_in_blocks_finished == FALSE`. On ARM64 neither
store→store nor store→load order is guaranteed for plain accesses. Therefore *j* can
- see `carry_in_blocks_finished == TRUE` from thread *i*'s **previous, already-finished** section
  (slot reuse at 5685-5689) while `start_block` is already the new split section's start; or
- see the current section's flag TRUE before *i*'s normalisation writes are visible;
and then **skip the wait and run `xgw_carries` over blocks that still hold pre-normalisation
inverse-FFT values** (`g->data_addr` blocks, `x/gw_carries` = `x/xmult2.c:xgw_carries`). It writes
back "stale data + carry", and the very next thing that happens is the state-3 *post-FFT forward FFT*
on exactly those blocks (`gwnum.c:6415-6430`) — no normalisation/roundoff check runs on them again in
this operation. Result: the product is wrong by sub-ULP amounts, i.e. the classic "wrong but
numerically clean" result that only trips `ed_check` much later.
- Note `gwevent_signal`/`gwevent_wait` do provide release/acquire (pthread mutex), but only on the
  *slow* path. The fast-path flag read takes no lock, so it has no acquire.

Why the observations match:
- MT-only (plain field accesses are per-operation, not per-thread state).
- Intermittent: window is the store-visibility gap; also a logic race on slot reuse.
- Worse with 8 threads than 4: `min_split_size = num_postfft_blocks + 2*clm = 8+8+8 = 24`
  (`gwnum.c:5597-5606`) while `num_pass1_blocks = PASS2_SIZE/2 = 128` (`gwnum.c:6595-6597`), so at
  most ~5 sections exist; with 8 threads some helpers finish early and *steal/split* other live
  sections mid-operation, i.e. exactly the slot-reuse/reorder path above. With 4 threads (0,24,56,88)
  there is little or no stealing.
- ECM vs PRP: ECM stage 1 performs ≈10 `gwmul3`-class FFT passes per ladder bit (`ecm.cpp:2006-2096`,
  `ed_add`/`ed_dup`), PRP does 1 `gwsquare2` per bit → ~10× more handshakes per unit time, so the
  same low per-FFT probability surfaces first in ECM. (PRP with `GWMUL_STARTNEXTFFT` also runs the
  handshake, so it is not immune — only less exposed.)
- "Clean" corruption: the reader replaces normalised values by pre-normalised ones (difference <1 per
  word) and adds a carry; forward-FFT then spreads a ~1-ULP difference over the spectrum, which is
  the same magnitude as ordinary roundoff, so gwnum's maxerr never reaches 0.5.

Confidence: MEDIUM-HIGH for the mechanism being a real ARM64 port bug; MEDIUM that it is *the*
observed corruption (needs the experiment below).

Cheapest experiment to confirm/refute:
1. **TSan build on x86 native** (fastest, no ARM needed): build gwnum with
   `CC=g++ CFLAGS/CXXFLAGS="-fsanitize=thread -O1 -g"` using `gwnum/Makefile.arm64 ARCH=x86`,
   link `_work/mt/mtstress.c` (differential ST-vs-8T harness already present) with the same flags,
   run `mtstress --threads 8`. Expect TSan reports on `gwnum.c:6317` vs `:6368-6369`
   (`carry_in_blocks_finished`) and on `gwnum.c:5685` vs `:6368` (`start_block`). (TSan on x86
   flags unsynchronised plain accesses even though x86 hardware usually hides them.)
2. **Barrier patch on ARM64** (2-line experiment, do not ship): in `gwnum.c` make
   the flag accesses atomic: writer lines 5689 and 6317 →
   `__atomic_store_n(&sections[..].carry_in_blocks_finished, v, __ATOMIC_RELEASE);`; reader line
   6369 → `__atomic_load_n(..., __ATOMIC_ACQUIRE)` (and read the flag *before* `start_block`).
   Rerun the 8-thread ECM reproducer; if corruption disappears, S1 is confirmed.
3. **Probability amplifier** (no patch): make the window deterministic with an instrumented debug
   build — insert `__sync_synchronize()`-free `nanosleep`/`sched_yield` right after line 5689 and
   right after 6317 when an env var is set, and/or add a delay after the first-block normalisation.
   On ARM64 the unmodified code should then corrupt within seconds.

### S2 (MED-HIGH) — same class, other plain fields of the section state machine
`gwnum/gwnum.c:5615-5627` (reads other sections' `section_state`, `last_block`, `next_block` while
choosing a section to split), `5702-5704` (writes own `next_block`, `section_state` **after** the
lock at 5696 was released), `6327-6331` (`next_block`/`last_block` read before taking the section
lock), `6396-6435` (`section_state` transitions). All are plain `int`s (`gwnum.h:923-931`); the only
locking is `change_in_progress` (`gwnum.c:5567-5578`) and most readers never take it. On ARM64 a
splitter can act on an inconsistent snapshot (e.g. `next_block` already bumped but `section_state`
not yet visible) and split a section at the wrong point, or a thread can see `section_state==1`
while the owner has already consumed the last block. Fix is the same discipline
(release on the state writes, acquire after the lock/at the checks), so the S1 experiment covers it.

### S3 (MED) — `MAXERR`/`XMM_MAXERR` handshake only copied in pass-1 state 1
`gwnum/gwnum.c:5947-5960` copies `XMM_MAXERR` to helpers only when `pass1_state==1`; helpers never
re-seed it for pass-2 or `PASS1_STATE_MULTITHREAD_OP`, and the main thread's `XMM_MAXERR` is never
reset between operations (only merged in `pass1_get_next_block`, `6253-6262`). If a helper ever reads
`XMM_MAXERR` in another state (the pass-1/2 normalisers do read it: `x/xmult2a.h:35`, `x/xnormal.h`,
`x/xnormal.c`) it inherits a stale/garbage maximum; this cannot change the arithmetic, but it can mask
or fabricate roundoff reporting. Relevant because `f195b22` "neutralised" a spurious `MAXERR=0.5`
readout on ARM64 instead of root-causing it — worth re-checking with the MT reproducer whether the
real maxerr is being hidden. Not believed to be the ed_check cause.

### S4 (LOW) — `asm_data` snapshot taken once per helper
`gwnum/gwnum.c:5869-5888`: a helper `memcpy`s the whole main-thread `asm_data` on its first
activation and only refreshes a whitelist per call (`5900`, `5920-5961`). I cross-checked the
converted kernels' complete read set against this whitelist (see "checked and ruled out") and found
no per-operation field that is missed *for the SSE2 paths*. The residual risk is ordering: the
`memcpy` reads fields the main thread writes just before `signal_auxiliary_threads`; the event
mutex does give acquire/release, so this is believed sound, but it is a place to look if S1 is
refuted.

### S5 (LOW) — ARM64 CPU flags / X86_64 macro interactions
The port defines `X86_64` on aarch64 (`gwnum/gwcommon.h:19-24`). I audited every `#ifdef X86_64`
in `gwnum.c` (lines 30, 286, 335, 348, 5044, 5253, 10322, 10462, 10728, 10760, 10879, 11129, 11398,
11666): all "64-bit asm feature" paths are additionally gated on `CPU_AVX*|CPU_FMA3`, which cpuid.c
forces off on ARM64, so the AVX-only shortcuts (two-destination FFT, fast `gwunfft`, direct
`gwmulmuladd5`, FFT(1)-is-one) are inert; `NEW_STACK_SIZE=0` matches the allocation code. No
MT-specific consequence found. (Also verified `sse2_explode9K8` maps K8→CORE on ARM64,
`gwnum.c:348-352`, consistent with cpuid.c leaving CPU_3DNOW_PREFETCH/CPU_SSE41 clear.)

## Checked and ruled out (evidence)

1. **No new shared/static state in converted kernels.** Full scan of `r4/ hg/ x/` (`*.c` + `*.h`):
   no function-local `static` (only commented-out debug lines `x/xmult1a.h:9,59`; the `P1()` debug
   macro `r4/r4pass1sc.h:48-49` is compiled out because `DEBUG524` is undefined); no writable
   file-scope objects other than `static const` trig tables (`x/xarch.h:166-262`) and
   `x/xarch.h:15-18`; only `static` helper functions (`x/xmult1ax.c:94`, `x/xmult2ax.c:117`,
   `x/xmult3ax.c:14`, `x/xmult1a.h` macro-generated).
2. **Converted kernels never touch `gwdata`, `parent_gwdata` or `asm_data->gwdata`.**
   `grep -n 'gwdata|parent_gwdata' gwnum/{r4,hg,x}/` → 0 hits. All per-thread state is reached only
   through the per-thread `struct gwasm_data *g` parameter. (Task focus area 2 is thereby ruled out
   for the kernels; the `asm_data->gwdata` use is confined to the C callbacks in `gwnum.c`.)
3. **`thread_num` is invisible to the kernels.** It appears only in `gwnum.c`
   (`auxiliary_thread`, `pass1_state1_*`, `pass1_get_next_block_mt`, `pass1_post_carries`); neither
   the converted C nor the `.mac` originals reference it (grep `thread_num` in `*.mac` → 0). The
   main/helper distinction in the asm path is entirely via the `PASS1_EXIT_THREAD=3` return code of
   `pass1_get_next_block_mt`; the C macro reproduces that: `r4/r4pass1sc.h:56-66`,
   `r4/r4dwpnpass1sc.h:56-66`, `hg/hgpass1.h:57-67` vs `xmult.mac:303-327`
   (`jl c0b / je b0b / jl done / je pass2 / else xgw_carries → again`). The C form
   (`if (branch1&2) return; if (branch1!=5) break;`) is equivalent for codes 0..5; codes 0/4 fall
   through the macro body to the outer `do{...}while(1)` dispatch (`r4/r4pass1sc.h:36-38`) exactly
   like the asm `goto c0b`/fall-through-to-pass2.
4. **MT section indexing is not done in the kernels.** `start_block`/`carry_out_section`/
   `pass1_carry_sections` never appear in `r4/ hg/ x/`; `norm_ptr1`/`norm_ptr2`/`carries`/
   `data_addr`/`premult_addr` are absolute pointers computed in `gwnum.c`
   (`pass1_state1_norm_addrs` 5478-5494, `pass1_state1_carry_addrs` 5496-5514,
   `pass1_state1_assign_*` 5582-5751) and consumed verbatim (`x/xmult2a.h:36-49`,
   `x/xmult2.c:78-84`). `this_block` is used by the kernels only in the two legal asm places:
   `this_block==0` carry wrap/sign (`x/xnormal.h:1860`, `xnormal.mac:2531`) and
   `this_block==last_pass1_block` top-carry (`x/xmult2a.h:93-96`, `xmult2a.asm:100-104`). So "offset
   relative to section instead of absolute" (task focus area 3) is ruled out for the converted code.
5. **`xgw_carries` faithful to asm.** `x/xmult2.c:66-136` vs `xmult2.asm:110-166`: carries pointer
   (`g->carries`), `rp = data_addr`, per-row `rsi+=64 / rbp+=pass1blkdst / rdi+=clm*4 / rdx+=128`,
   `RATIONAL_FFT`/`B_IS_2` branches, zpad variant and `loopcount1=addcount1` all match.
6. **`xnorm012_2d` + `_part1` faithful.** `x/xnormal.h:1858-2308` vs `xnormal.mac:2528-2983`:
   part-1 `count3` 3-section shuffle loop, negate-last-carry (`-=BIGVAL; *=-C; +=BIGVAL`),
   wrap into carries[0]/[2], the `SPREAD_CARRY_OVER_EXTRA_WORDS` 2-word vs 6-word paths,
   `BIGLIT_INCR2/4` offsets, and the final 4-carry clear all match instruction-for-instruction
   (including rounding via `roundf64` ≡ the add/sub-BIGVAL trick).
7. **`inorm`/`xnorm_2d` (the pass-1 normaliser) faithful.** `x/xmult2a.h:21-97` vs
   `xmult2a.asm:45-108`: register roles are the asm's (`rsi`=FFT data, `rbp`=carries,
   `rdi`=norm_ptr1/biglit, `rbx`=norm_ptr2/col-mults copy); `setlp` precompute loop, per-row
   `rsi+=normblkdst`, `rbp+=64`, `rdx+=128`, 8-iteration `normblkdst8` cadence, `XMM_MAXERR` save,
   and the top-carry call gated on `this_block==last_pass1_block` all match.
8. **SSE2 add/sub/addsub/smallmul are single-threaded on ARM64** — important scoping result:
   `gwadd3o/gwsub3o/gwaddsub4o` only call `multithread_op` when `cpu_flags & (CPU_AVX512F|CPU_AVX)`
   (`gwnum.c:12736-12793`, `12883-12940`, `13033-13104`), `gwsmallmul` likewise
   (`gwnum.c:13316-13330`), and `cpuid.c` forces ARM64 to SSE..SSE42. So ECM's point operations are
   single-threaded; the only MT-exposed code for this workload is `gwfft/gwmul3/gwsquare2`
   (pass1+pass2) and `gwcopy`. This eliminates a large candidate space (multithread-op helper setup,
   `do_multithread_op_work`, 4KB/AVX-512 block partitioning, `gwy3_apply_carries` path).
9. **Helper copy list completeness.** I extracted *all* `g->FIELD`/`g->u.xmm.FIELD` references from
   every compiled converted source (142 + 31 distinct names) and compared them with
   `auxiliary_thread`'s copy list (`gwnum.c:5900`, `5920-5961`) and with `internal_gwsetup`'s static
   initialisation (`gwnum.c:3340-5240`). Every kernel-read field is either (a) static after
   `internal_gwsetup` (and therefore covered by the one-time `memcpy` at 5870: `NEGACYCLIC_FFT`,
   `B_IS_2`, `RATIONAL_FFT`, `ZERO_PADDED_FFT`, `SPREAD_CARRY_OVER_EXTRA_WORDS`, `TOP_CARRY_NEEDS_
   ADJUSTING`, `ADDIN_*` constants, `XMM_LIMIT_*`, `XMM_K_*`, `XMM_MINUS_C`, `XMM_NORM012_FF`,
   `XMM_T*_FUDGE`, `sincos*`, `xsincos_complex`, `plus1_premults`, `norm_grp_mults`,
   `norm_col_mults`, `norm_biglit_array`, `count*`, `normval*`, `pass1blkdst`, `pass2blkdst`,
   `pass2gapsize`, `normblkdst*`, `cache_line_multiplier`, `last_pass1_block`, callbacks), or
   (b) written per call for that thread (`DESTARG` 5900; `DIST_TO_*`, `SRC2ARG`, `SRC3ARG`,
   `DEST2ARG`, `ffttype`, `mul4_opcode`, `thread_work_routine` 5920-5927; state-1 set
   `NORMRTN`, `const_fft`, `add_sub_smallmul_op`, `ADDIN_*`, `XMM_MULCONST*`, `XMM_MAXERR`
   5928-5960; `DBLARG` for multithread-op at 5893), or (c) per-block pointers/counters set by the
   block-assignment routines on the thread's own `asm_data` (`this_block`, `next_block`,
   `data_addr`, `data_prefetch`, `premult_addr`, `premult_prefetch`, `norm_ptr1`, `norm_ptr2`,
   `carries`, `scratch_area`, `thread_num`). No gap for the SSE2 pass1/pass2 paths, so "helpers use
   the main thread's `g`/stale per-thread field" (focus area 4) is ruled out.
10. **Delayed-norm/`addcount` 16K area does not apply here.** For M264169 (k=1, b=2, c=-1)
    `TOP_CARRY_NEEDS_ADJUSTING = (k>1.0 && !ZPAD) = FALSE` (`gwnum.c:4931`), so
    `xnorm_top_carry*` is a no-op for this workload (`x/xnormal.h:1733-1734`); and
    `SPREAD_CARRY_OVER_EXTRA_WORDS` only picks the same 8-block post-FFT window for clm=8.
    The delayed-normalisation bug fixed at caller level in `72bb02a` (ECM `ed_add`: force-normalise
    E/H because an unnormalised multiply input is mishandled at 16K) is a *single-thread* defect that
    is already worked around in `ecm.cpp`; nothing in this audit explains an MT-only variant of it.
11. **`gwthread.cpp`/`gwthread.h` atomics are sound for ARM64** (`std::atomic` with
    acquire/release, `__yield` spin, pthread-cond manual-reset event emulation that mirrors the
    Windows semantics, including `signalled_count` to avoid lost/spurious wakeups).

## Recommended next actions (cheapest first)

1. Let the existing `_work/mt/mtstress.exe` (differential ST vs 8T, x86 build of the same converted
   C) finish `_work/mt/soak8-n264169.log` / `soak8-n291143.log`. Two outcomes:
   - mismatch on x86 → the bug is in shared C, and a local gdb/TSan session is possible (fast path);
   - clean on x86 → strongly supports an ARM64 memory-model-only cause (S1/S2) and/or a
     compiler-on-aarch64 effect.
2. Run the TSan x86 build as in S1 experiment 1 and record the reported pairs of lines.
3. On the ARM64 box, apply the release/acquire patch of S1 (two writers + one reader; note the
   required order: store `start_block/next_block/last_block` first, then RELEASE-store the flag;
   reader: ACQUIRE-load the flag, then read `start_block`) and rerun the 8-thread ECM.
4. If S1/S2 are refuted, instrument `xnorm012_2d`/`xnorm_2d` to hash every block's data right after
   the pass-1 normalise and right after the state-3 post-FFT, plus a per-block "carry added" log, to
   localise the first bad block to (normalise | carry-add | post-FFT | pass2).

## Live evidence at time of writing

`_work/mt/mtstress.exe` (x86-64 native build of the same converted C; ARMs flags faked to
`0x3F2F`) is already running the decisive differential test, started 01:30:
- `_work/mt/soak8-n264169.log`: phase-1 single-thread reference of 1,000,000 iterations of
  ops `0..5` (`gwmul3`/`gwmulmuladd5`/`gwsubmul4`/chained-FFT variants) finished in 557 s at
  1795 it/s, "type-1 FFT length 16K, Pass1=64, Pass2=256, clm=4"; phase-2 8-thread run started.
- `_work/mt/soak8-n291143.log`: same, reference 484 s, phase 2 started.
If either reports a digest mismatch, the bug is in the shared C (converted kernels *or*
gwnum.c) and is debuggable locally; a clean x86 soak supports an ARM64-specific memory-model
cause (S1/S2). Note the x86 host is not a memory-model proxy: x86 TSO hides S1/S2 by
construction, so a clean x86 soak is *expected* under the S1/S2 hypothesis.

## Caveats

- This is a static audit; nothing was executed and no gwnum source was modified.
- The converted-kernel conclusions are based on spot line-by-line comparison of the macros on the
  SSE2-r4-two-pass legal path (r4 `xfft_r4_16K_*_CORE` → `r4/r4pass1sc.h`, `xr2*`/`gw_carries` from
  `x/xmult2a.c`/`x/xmult2.c`, `xnormal.h`); the hg/ and dwpn/ variants were only grep-audited (they
  are not used by this 16K workload).
- S1/S2 are *upstream gwnum.c* issues, not asm→C conversion issues, but they are the only MT-only
  correctness hazards found that match the "intermittent + clean wrong result + worse with more
  threads" signature.
