# Package simplification — design

Date: 2026-10-06. Branch: `package-simplification` (from main `e59ebd6`). One branch, one PR.

Companion file: `2026-10-06-package-simplification-items.md`, the approved item catalog. It
holds each item's full change text, evidence and test impact, keyed by the IDs this spec uses
(`types/T5`, `enum/S1`, `time/TT02`, `loc/TY1`, …). This spec sets the goals, the order, the
gates and the decisions; the catalog supplies the detail.

## 1. Goal

Make EnzymeRates smaller, simpler and faster to develop on, without losing documentation,
coverage or correctness.

Success means all of the following:

1. `src/` shrinks by the approved items: about 4,070 lines net from 10,857. Comments and
   docstrings of surviving code stay.
2. The test suite takes at most ~550 s of test time (baseline 1,104–1,152 s), measured with the
   timed runner of §9.
3. No coverage is lost. Every deleted test was signed off by Denis or tested only code this
   branch deletes. Output stays pristine.
4. Four invariants hold unchanged: the `rate_equation` performance contract (0 allocations,
   < 120 ns), Canonical Step Form, the `name(p, m)` chokepoint, and Pass-2 Wegscheider absorption.
5. Every behaviour-preserving commit reproduces the equivalence harness (§4) byte for byte.

## 2. Scope

In scope: the approved source items, the performance fix (§5), the approved test items, the bug
fixes Denis approved, and the doc edits these force.

Out of scope, each a follow-up PR:

- Making per-regulator multiplicities work (`A(1, 2)` in `@enzyme_reaction`). Today the syntax
  parses and is stored in `RegulatorMults.allowed_multiplicities` and the Sig, but no move or
  derivation reads it: regulatory sites always take the catalytic multiplicity. This branch
  leaves the syntax and the field untouched.
- Moving dead-end decorations out of `init_mechanisms`.
- The K-type `n = 1` L phantom.

## 3. Baseline measurements

Full suite on main: 19m15s wall, 1,152 s of test time, 64,694 tests, 3.9 GB peak memory.

An instrumented run split each top-level testset into compile, GC and runtime: 1,104 s total,
514 s compile (47 %), 205 s GC (19 %). The largest testsets:

| Testset | Time | Compile | GC |
|---|---|---|---|
| `_hyperbolic_catalysis` matches the derived denominator | 258 s | 0 s | 115 s |
| Enzyme Derivation Tests | 81 s | 62 s | 4 s |
| init division-freeness (bi_bi_pp) | 73 s | 73 s | 2 s |
| compile-budget (subprocesses) | 72 s | — | — |
| catalytic moves on bi-bi to depth 2 | 70 s | 0 s | 39 s |
| rate-eq dedup-key partition stability | 56 s | 52 s | 2 s |
| `init_mechanisms` on ter-ter within 150 s | 45 s | 0 s | 11 s |
| Fitting | 44 s | 6 s | 3 s |

A profile of the hyperbolic sweep found the dominant cost. About 65 % of enumeration samples and
60 % of derivation samples sit in `_free_enz_set`. `name(p, m)` calls `_rep_step`, which rebuilds
the mechanism's free-enzyme set (three passes over all steps, rendering species names to Strings)
on every call. `_assemble_constraints` reaches it twice per lookup through its `step_name`
closure. A throwaway spike that memoized the set per mechanism cut the ping-pong depth-2
enumeration from 86 s to 30 s (122 GB → 45 GB allocated) with identical output.

## 4. Equivalence harness

Before any source change, record a golden snapshot; every behaviour-preserving commit must
reproduce it.

- **Reactions:** uni-uni, uni-bi, ordered bi-bi, ping-pong-capable bi-bi, one ter-ter seed set,
  and one reaction with an allosteric regulator and a competitive inhibitor.
- **Recorded per reaction:** `init_mechanisms`, `expand_mechanisms` to depth 2 (depth 1 for
  ter-ter), and `seed_mechanisms` where it applies. Record each mechanism's
  `string(compile_mechanism(m))` in enumeration order, its `rate_equation_string`, and its
  `fitted_params`.
- **Storage:** a scratch script and its output files, outside the repo. The harness is a refactor
  tool, not a test; the suite keeps its own pins.
- **Approved behaviour changes:** each lands in its own commit. That commit regenerates the
  snapshot, and its message states what changed and how many mechanisms or strings moved.
- **Enumeration order:** behaviour-preserving commits keep it. A commit that changes order only,
  with the same set, must say so and be approved.

## 5. Performance fix

### 5a. Lazy per-mechanism cache (approved option)

- `Mechanism` and `AllostericMechanism` gain one non-semantic field: a lazily filled cache
  (a mutable container) holding the free-enzyme set and each kinetic group's naming
  representative. The field stays out of `==`, `hash`, the Sig encoder and `show`. Canonical Step
  Form, dedup and compiled types are therefore unchanged.
- `Species` stores its rendered name at construction; `name(::Species)` returns it. The stored
  name stays out of `==` and `hash`.
- `_rep_step` reads the cached representative. `name(p, m)` keeps its signature, so the
  chokepoint and its AST-walker test stay as they are.
- `_assemble_constraints` renders each step's names once (thermo/T7).
- Gates: the harness, the perf test, and the full suite. Re-time the suite right after this
  commit and report per-testset numbers before starting the test commits.

The generated `==`/`hash` loop (types/T5) must skip the cache fields. Land T5 after this commit,
or write the skip into T5.

### 5b. Per-type JIT

Helpers reached from `rate_equation_string`, `parameters`, `fitted_params` and `_kcat_forward`
re-specialize on every `EnzymeMechanism{Sig}` type. That is why "dedup-key partition" spends 52
of its 56 s compiling. Lift to the concrete mechanism at the entry point and mark the singleton
argument `@nospecialize` where the body does not need the type (time/TT35). `rate_equation`'s
generated body still compiles per type, as its performance contract requires. Gate: measured
compile time before and after on 50 fresh bi-bi mechanisms, plus the harness and the perf test.

## 6. Order of work

One or more commits per phase. Run the full suite before every commit; iterate with focused
runs.

1. **Harness** (§4), recorded on the untouched branch.
2. **Performance fix** (§5a), then re-time. Then §5b.
3. **Test-speed items that need no source change** (§8), so later phases run on a faster suite.
4. **Pure deletions** of dead and test-only source code, with the tests whose only subject is
   that code.
5. **Consolidations**, in the order the cross-file critic set (catalog, "Cross-file critic"):
   types foundations (T5 → T13 → T14 → T6 → T12 → T17 → T4 → T3 → T7 → T8 → T9 → T10 → T11),
   then thermo and derivation, then enumeration (including S1, S4, F), then DSL (including D7),
   then identify and fitting. Items that share call sites land together, as the critic's
   resolutions say: types/T17 with enum/R1, R2 and deriv/V6–V8; deriv/V11 with enum/RC;
   thermo/T7 before deriv/V14–V16; thermo/T8 with deriv/V14 and V18; deriv/V18 owns the body
   builder (thermo/T9 shrinks accordingly).
6. **Approved changes that alter recorded outputs** (enumeration sets or order, rate strings,
   parameter names, `show` text, CSV or log text), one commit each, each regenerating the harness
   and re-pinning its tests (§7). Approved removals of internal or unused methods that leave the
   recorded outputs unchanged (types/T7, T10, deriv/V27, identify/V2) may land in phase 5.
7. **Finish:** Aqua and JET, the docs build and doctests, docs edits for removed or renamed
   internals (`developer.md`, `enumeration_engine.md`, the fitting tutorial), version 0.9.0
   (breaking: `parameters(::AllostericEnzymeMechanism, Full)` is removed), and a final timed run
   for the PR description.

Use group-prefixed IDs in commits and the plan (`types/T1`, `thermo/T1`), because bare IDs
collide across groups.

## 7. Decisions

From Denis's sign-off page (2026-10-06).

**Rewrites approved** (each gated on the harness): enum/S1 `_catalytic_topologies` as one
walker; dsl/D7 DSL step parser emits Exprs directly; enum/S4 dead-end builder as one
connectivity fill; thermo/T2 dependent-parameter solve as column-order RREF; enum/F delete
`_seed_candidate_screen`, only if timed ter-ter and bi-bi ping-pong `init_mechanisms` runs are at
most ~5 % slower without it. S4 also asserts that no two distinct species in a mechanism render
the same name (the critic's `:ENADH` case).

**Bug fixes approved:** types/T20, dsl/D3, D9, D12, D14, D15, D18, deriv/V23, V24, V25,
enum/S5, D2, identify/V3, B8, C2. deriv/V24 follows TDD: write the failing numerical test first,
and fix only if it fails.

**Output and API changes approved:** types/T7, T10, T21 with enum/A9, T22 (three Parameter
types; update the chokepoint test's type regex), T23 (drop the chain display), thermo/T11 (one
name sort, which replaces deriv/V26), deriv/V27, identify/B6, V2, dsl/D8, D10, D19, enum/A5.

**Declined.** These items stay as they are:

- thermo/T10: pivot priority keeps today's rule; correct the `_step_priority` docstring to match
  the code.
- identify/FP1: keep `FittingProblem(::Mechanism, …)`.
- identify/FR2: keep `lb`/`ub` (they set the optimizer box).
- identify/FR3: keep `abstol`/`reltol`/`callback`.
- identify/B7: keep the `fit_inherited` column.
- enum/A10: keep `_requires_hyperbolic_catalysis` as the extension hook.
- types/T19 with dsl/D4: the follow-up PR of §2.
- identify/B9: moot.

## 8. Test plan

Denis approved three levers: cheaper with the same coverage, deletion of subsumed tests, and
lower statistical thresholds. No opt-in slow tier.

**Approved time items** (catalog prefix `time/`; a merged `loc/` item lands with its partner):
TT01 with loc/T19, TT02, TT03 with loc/CB2, TT05 with loc/FI1, TT06, TT10, TT11, TT12, TT13 with
loc/ID1, TT14, TT17 with loc/ID3 and loc/ID4, TT18, TT20 with loc/ID11, TT22, TT23, TT25, TT27
with loc/ID2, TT29 with loc/T22, TT31 with loc/T16, TT32 with loc/T11, TT33 with loc/T20, TT34,
TT35, TT36. Items that need no sign-off because coverage is unchanged are also in: TT04, TT08,
TT09, TT15, TT16, TT21, TT24, TT26, TT28, TT30. Declined: TT07 with loc/DR1 (the QSSA and
stiff-ODE oracles and the OrdinaryDiffEqFIRK test dependency stay) and TT19. loc/DR16 depends on
loc/DR1, so it drops out.

**Approved code-size items:** every `loc/` item in the catalog, plus dsl/D21, identify/T2, T3,
CV2, CV3 and deriv/V7.

**Projection** (seconds of test time, from measurements and spikes):

| Step | Saving | Suite |
|---|---|---|
| Baseline | | 1,104 |
| §5a cache | ~250 | ~850 |
| §5b per-type JIT | ~40 | ~810 |
| Division-freeness on polynomials (TT02) | ~60 | |
| Compile-budget without the duplicate ter-ter subprocess (TT03) | ~47 | |
| Rank oracle without derivation or compilation (TT06, TT09) | ~37 | |
| Hyperbolic sweep: drop bi-bi, memoize per state (TT01, TT04, TT14) | ~35 | |
| BBO fits at maxtime 0.5 s (TT05) | ~25 | |
| Remaining items | ~60 | ~550 |

**Rules:**

- An item marked "measure first" in the catalog lands only if its measurement holds. Example:
  TT12 lands only if the allosteric identify run still recovers the mechanism at
  `max_param_count` 6.
- A test rewritten from compiling `rate_equation` to checking symbolic polynomials must assert
  the same property.
- A test whose only subject is deleted source code goes with that code.
- Enumeration-engine tests keep the CLAUDE.md rules: inline fixtures, exact child sets and
  counts, chemistry written the enumerator's way.

**Checkpoint.** After §5a, re-time the suite. If the projection falls short, the next levers are
the in-process ter-ter `init_mechanisms` pin (45 s) and the ter-ter seed and split tests
(~27 s). Any test change beyond the approved list goes to Denis first.

## 9. Verification

- **Each commit:** the full suite, green and pristine, plus the harness for behaviour-preserving
  commits.
- **Performance:** `test_rate_equation_performance` stays green in every commit, and §5 is
  re-timed as stated there.
- **Timing:** a scratch runner prints every nested testset's time (`@testset verbose=true`
  around each file) and the per-file compile and GC split
  (`Base.cumulative_compile_time_ns`). The PR description carries the before and after tables.
- **Lines:** `wc -l src/*.jl test/*.jl` before and after, in the PR description.
- **Docs:** the docs build and doctests at the end; every renamed or deleted internal that a
  docs page names is updated in the commit that changes it.

## 10. Risks

- **Canonical Step Form** (types/T13, T14, T16, T21): each lands alone, gated on the Canonical
  Step Form tests and the harness.
- **Pivot solve** (thermo/T2): the critic re-derived the equivalence argument. The harness
  checks it on every recorded mechanism.
- **Pass-2 absorption** (deriv/V9, V14, V16, thermo/T7, T8): gated on the exact strings of the
  non-competitive-inhibitor and non-essential-activator families, which use the absorption.
- **Cache field** (§5a): a stale or semantic cache would break dedup. The cache derives only from
  `steps(m)` at construction and stays out of `==`, `hash` and the Sig. The harness and the
  dedup tests check it.
- **Large rewrites** (S1, S4, D7): each lands in its own commit with the harness, so a regression
  bisects to one rewrite.
