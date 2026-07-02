# LDH HPC Run — Bug Fixes and Search-Engine Improvements

**Date:** 2026-07-01
**Branch:** several-bug-fixes (design); implementation branching decided in the plan.
**Source:** `docs/ldh_hpc_results/2026_06_27_results` (a 2.4 GB LDH `identify_rate_equation`
run, `EnzymeRates` at `main`, 2026-06-27) and `docs/ldh_hpc_results/identify_ldh.jl`.

## Purpose

Use the LDH run to fix real bugs and improve `identify_rate_equation`. Denis flagged six
areas; analysis (100,787 fits distilled, plus reproduction on current code) resolved them
into **five code fixes** (§5a, §5b, §1, §2, and the trivial §3), **one report** (§6), and
**§4** — dedup is not broken, so the canonical-key upgrade is deferred and LOOCV's existing
`eq_hash` dedup is kept (a confirming test only). This spec designs all of them for a single
implementation effort.

## Findings (evidence)

The run selected a sensible 8-parameter `AllostericEnzymeMechanism` (train loss 0.0119),
so the pipeline works — but it was grossly inefficient and shed most of its candidates to
a codegen crash.

| Signal | Value |
|---|---|
| Total fits | 100,787 |
| Distinct rate equations (`eq_hash`) | 39,970 (**2.5× fit redundancy**) |
| Equations re-fit across ≥2 iteration files | 16,630 |
| Global best train loss | 0.00932, first reached at **iteration 8** |
| Iterations run | **22** (labelled "target n_params" up to 26) despite `max_param_count=13` |
| Errored fits | 54,702 — **UndefVarError 54,674**, `_kcat_forward` "no components" **28** |
| Non-Success retcodes | 310 (0.31%), **100% `MaxTime`** |

All error and inefficiency claims were reproduced on current code (branch
`several-bug-fixes`), not just read from the run.

## Design

### §5a — Inactive-state parameter elision (`UndefVarError`) — correctness bug

**Root cause.** The allosteric derivation *defines* inactive-state (I) symbols and
*references* them at separate sites that choose the I-symbol set by different criteria, so
`DEFINED(I) ⊉ REFERENCED(I)` and a referenced symbol dangles. `_build_dep_assignments`
(`src/rate_eq_derivation.jl:1470`) predicts I-symbols from group structure
(`i_names_set`, `:1501-1518`, with an `i_dead && tag !== :NonequalAI` gate and an
`if !i_dead` synthesized-dep block); `_allosteric_num_den_exprs` (`:1581`) references
whatever survives into the retained `Q_I` (and `N_I` when the I-cycle is live). PR #54
retained `Q_I` as enzyme mass but reasoned only about deps touching an `:OnlyA` symbol;
it missed the cases below. Three trigger paths, one defect:

1. **`i_dead` Case-B synthesized deps** (majority): a dep whose Haldane RHS references a
   `:NonequalAI` symbol gets an I-name only inside `if !i_dead`, so when `i_dead` the
   I-name never enters `i_names_set`, is never assigned, yet survives into `Q_I`
   (e.g. `k_I_ELactateNAD_to_ENADHPyruvate`, `K_I_Lactate_ENAD`).
2. **`i_dead` phantom**: an emitted `:NonequalAI` dep assignment whose RHS references an
   independent I-counterpart that the phantom-survival filter correctly dropped from the
   destructure — the assignment is dead but references an undefined symbol.
3. **non-`i_dead` group-rep/polynomial mismatch** (~6.5%, e.g. `kon_I_NAD_EPyruvate`):
   `i_names_set` collects only the representative step's I-name per group, but `Q_I`
   references per-enzyme-form I-binding names that don't fold to the rep.

**Fix.** Derive the I-symbol set from the *actual* reference polynomial, not from group
structure. Add one shared helper `_i_state_referenced_syms(am) → S_I = { I-state symbols
in den_i_poly } ∪ ({ in num_i_poly } when the I-cycle is live)`, factoring `den_i_poly` /
`num_i_poly` out of `_allosteric_num_den_exprs` so definition and reference read the same
polynomials. Then:

- **`_build_dep_assignments`**: replace the group-structure `i_names_set` construction
  with `i_names_set = S_I`; keep the existing catalytic-dep loop and its `a_only → 0`
  zeroing. This assigns surviving Case-B deps (path 1) and excludes non-surviving phantom
  names (path 2).
- **`_dependent_param_exprs`**: gate `indep_I_list` and `dep_I` on `S_I` **plus a
  single-pass RHS closure** — include any independent I-name appearing in the substituted
  RHS of an emitted dep. (Closure terminates in one pass because dep RHSes reference only
  independent params — the Gaussian-elimination invariant; a surviving dep's Haldane RHS
  can reference an independent I-name not itself in `Q_I`, so closure is required —
  verified 10/40 mechanisms still dangle without it.)
- **Single source:** `S_I` and the closure live in one helper consumed by both sites so the
  assignment set and the destructure set cannot desync (the top implementation-regression
  risk).
- **Do not** fold regulator (`Kreg`) mirrors into `S_I`; `S_I` is catalytic-only. Keep
  reg-site definitions driven by reg-site structure exactly as today.
- **Companion:** reconcile `_synthesized_dep_i_names` (`:1263`, which early-returns
  `Symbol[]` when any `:OnlyA` group is present) with `S_I` so `parameters(Full)`
  enumerates the newly-emitted synth-dep I-names and the `name(p, m)` chokepoint AST-walker
  test (`test/test_types.jl:1577-1644`) stays consistent.

**Correctness / perf.** RHS expressions are untouched, so Haldane/Wegscheider algebra and
Canonical Step Form are preserved. The change is a **no-op for currently-working
mechanisms** (verified 0/120 indep-set mismatch), so `MECHANISM_TEST_SPECS` and the
0-alloc / sub-100 ns contract are unaffected; broken mechanisms gain flat scalar
assignments only. Prototype verified: **0/40 previously-failing mechanisms dangle** with
the fix.

**Tests.** **Add the actual failing LDH mechanisms as permanent `MECHANISM_TEST_SPECS`
entries** (`test/mechanism_definitions_for_test_enzyme_derivation.jl`), transcribed via
`@allosteric_mechanism_src` — **one per trigger path**: Case-B reverse `k_I_*_to_*` (e.g.
the mechanism dangling on `k_I_ELactateNAD_to_ENADHPyruvate`), Case-B binding `K_I_*` (e.g.
`K_I_Lactate_ENAD`), `i_dead` phantom `kon_I_*`, and non-`i_dead` multi-substrate
`kon_I_NAD_EPyruvate`. Their topologies are already captured from the reproduction, so
transcription is mechanical. As `MECHANISM_TEST_SPECS` entries they inherit the full
derivation battery (Full/Reduced parameters, analytical rate, 0-alloc/perf), so a future
regression fails the build. Plus: an invariant test (`DEFINED(I) ⊇ REFERENCED(I)`: every
free parameter symbol in the generated body is destructured or assigned) over that battery;
each fixture asserts a finite numeric `rate_equation` equal to the analytical MWC value; a
runnable-transcript test (`rate_equation_string` LHS-covers every v-line symbol); and a
0-alloc / sub-100 ns assertion on one of the new LDH `i_dead` fixtures (the existing perf
test uses PFK-1/HK, which do not exercise this path).

### §5b — `_kcat_forward` "no components" and `NaN` at products=0 — derivation bug

**Root cause (investigated directly, overturning the initial "zero-kcat" reading).** The 28
mechanisms are **not** zero-kcat. Their forward rate at saturating substrates as
products→0⁺ converges to a finite nonzero limit (13.78, 4.98, 4.96 for three samples);
only at *exactly* products=0 does it evaluate to `NaN`. The derived numerator and
denominator share an **uncancelled common factor involving products** — the denominator
has no constant (free-enzyme `1`) term, so its reference form is product-bound. Two
symptoms, one cause: (a) `_kcat_forward`'s substrate-only-monomial matcher finds nothing
and hits `isempty(a_keys)` (`:882`); (b) `rate_equation` returns `0/0 = NaN` at products=0,
which can corrupt fits on any zero-product data point. This is a division-free-derivation /
free-enzyme-reference artifact and overlaps the existing division-free work.

Returning `kcat = 0.0` (an earlier proposal) is **wrong**: it would stamp kcat=0 on an
enzyme whose kcat is ~14 and corrupt `scale_k_to_kcat` rescaling.

**Fix.** Root-cause the common factor in the numerator/denominator derivation and cancel
it (restore the free-enzyme reference / extend the conc-GCD reduction so the denominator
regains its constant term). This simultaneously restores the substrate-only numerator
monomial (so `_kcat_forward` works unchanged) and removes the products=0 `NaN`. The exact
cancellation site is not yet pinned; implementation begins by locating where the
free-enzyme-per-segment reference leaves a shared product factor for these topologies. The
change is compile-time (derivation) only, so the runtime perf contract is not at risk, but
it is in load-bearing code and requires full-suite TDD.

**Acceptance criteria.** For the 28 mechanisms: `_kcat_forward` returns a finite value
equal to the numerical grid-peak forward turnover (matching the products→0⁺ limit), and
`rate_equation` returns a finite value at products=0. No change to any currently-working
mechanism's rendered equation or kcat.

**Tests.** **Add ≥1 of the 28 as a `MECHANISM_TEST_SPECS` entry** (via
`@allosteric_mechanism_src`) with an `analytical_kcat_fn`, so kcat correctness is a permanent
guard: assert `_kcat_forward` equals the numerical grid-peak forward turnover and
`rate_equation` is finite at products=0. Plus a guard that no valid mechanism's denominator
lacks a constant term, and `rescale_parameter_values` round-trips
(`_kcat_forward(result) ≈ scale_k_to_kcat`).

### §1 — Parsimony filter references only `c-1`

**Root cause.** `src/identify_rate_equation.jl:579-581` sets `parsimony_cutoff =
loss_parsimony_threshold * best_loss_by_count[c-1]`, referencing only count `c-1`. Two
failure modes: (a) **count gaps** — an expansion move can add >1 parameter, so a tier can be
populated at `c` while `c-1` was never fit; then `haskey(...,c-1)` is false and the
parsimony filter is *silently disabled* at that count; (b) **non-monotone** `best_loss_by_count`
— if `best(c-1) > best(c-3)`, the `c-1` reference is looser than the true best simpler model.

**Fix.** Extract a pure helper:

```julia
function _parsimony_cutoff(best_loss_by_count::Dict{Int,Float64}, c::Int,
                           loss_parsimony_threshold::Float64)
    prev = [best_loss_by_count[k] for k in keys(best_loss_by_count) if k < c]
    isempty(prev) && return nothing
    loss_parsimony_threshold * minimum(prev)
end
```

Replace the inline expression with a call to it. `min` over all counts `< c` is ≤
`best(c-1)`, so the cutoff is strictly tighter-or-equal (more pruning; also dampens §2). Do
**not** add `loss_abs_threshold` to the parsimony term (the additive term guards the
*relative* cutoff against `best(c)→0`; the parsimony term is deliberately strict, and
`min_beam_width` remains the anti-collapse floor). Behavior-neutral on this monotone,
gapless LDH run; corrects gapped/non-monotone runs.

**Tests.** `_parsimony_cutoff`: `nothing` at the base tier; min over all `<c` not just
`c-1`; count-gap case (`c-1` absent) returns a cutoff rather than `nothing`; non-monotone
case picks the true minimum; wired into `_select_beam`, the new cutoff admits ≤ the old.

### §2 — Search runs to full structural exhaustion (fit-dedup only)

**Root cause.** Termination is only `isempty(frontier)` (`:566,612`); `_process_batch` drops
`n > max_param_count` before fitting, so the frontier only ever holds `n ≤ 13` entries, and
the loop runs until the entire reachable ≤13 structural space is enumerated and expanded.
`target` (`:613`) is a bare iteration counter (it reached 26 as a label artifact; no fit
exceeded 13). The frontier keeps all structurally-distinct mechanisms with **no `eq_hash`
dedup** (`:488`), so the 2.5× equation-identical mechanisms all re-fit; 16,630 re-fit across
iterations because `unique!` is per-batch.

*Confirmed against the data: all 100,787 rows are independent fits, none copied.* The only
existing dedup is `unique!` on mechanism structs, within a batch — never by `eq_hash`, never
across iterations. Repeated `eq_hash`es show per-row loss variation at the CMA-ES noise floor
(e.g. `bfe5523cdb919538` appears **128 times with 127 distinct loss values**, spread ~2e-13),
which is only possible if each was fit independently; a copy would be bit-identical. So that
one equation was fit 128 separate times — exactly the waste fit-dedup removes.

**Fix — fit-dedup by `eq_hash` (only).**  `eq_hash` is already computed pre-fit
(`:456-458`). Split `_process_batch` into two passes with a cross-iteration
`memo::Dict{UInt64, fitresult}` threaded from `_beam_search`:

- **PASS-1** (`pmap`): compile + cap-check + render each mechanism → `(mech, n, em_type,
  eq_text, eq_hash)` or `FitFailure`. No fit.
- **PASS-2** (master): fit **one representative per unseen `eq_hash`** at full `n_restarts`;
  store the **raw pre-rescale** fit (log-space optimum + loss + retcode) in `memo`.

  *How a duplicate's parameters are inferred (no re-fit).* Same `eq_hash` ⟹ identical reduced
  rate equation ⟹ identical `fitted_params` **names and order** (the destructuring line
  `(; …) = params` is part of the dedup key, and `fitted_params` *is* that independent set).
  The optimizer's objective — `loss!` over the shared rate function and the same data — is
  therefore byte-identical for a duplicate, so an independent fit would target the same
  optimum. We reuse the representative's raw optimum directly: the duplicate inherits the
  same fitted-parameter values by name, no new optimization.

  *Why rescale is re-run per mechanism.* `fit_rate_equation` returns params **rescaled** so
  that `_kcat_forward(m, result) ≈ scale_k_to_kcat` — it multiplies the SS rate constants by
  `scale_k_to_kcat / _kcat_forward(m, raw)`. Both `_kcat_forward` and `_ss_rate_constant_names`
  are computed from the mechanism **structure** (`steps(am)`, `_all_i_state_parameters`,
  `_raw_symbolic_rate_polys_allosteric`), i.e. per-structure `@generated`/enumeration code —
  *not* from the reduced equation. Mathematically kcat is a property of the shared reduced
  function, so the two structures *should* rescale identically; but we don't want the reported
  values to depend on that assumption (§5b shows this exact code can be structure-sensitive).
  So we apply **each mechanism's own** `rescale_parameter_values` to the shared raw fit. This
  guarantees every row's reported params equal what a standalone fit-then-rescale of that
  mechanism would produce, and it costs one `_kcat_forward` evaluation per row (no
  optimization). `loss` is computed **pre-rescale** (`fitting.jl:256-261`), so it is
  `eq_hash`-invariant and **selection is unaffected** regardless — the rescale only touches
  reported parameter *values*.

  *Observability.* Add a boolean column **`fit_inherited`** to each result row (and the CSVs):
  `false` for the representative that was actually fit, `true` for a row whose fit was reused
  from an earlier same-`eq_hash` fit. This lets a run be audited for how much fitting the
  dedup saved and which rows share an optimum.

Savings: 100,787 → 39,970 fits (~60% fewer), fitting the dominant cost (`n_restarts=20` ×
CMA-ES at `maxtime=60 s`).

**Scope of the guarantee.** Fit-dedup is **structural-coverage-neutral** — the same *set* of
distinct equations is fit, each once at full budget — but **not selection-invariant**: it
replaces "min over K duplicate fits" with a single full-budget fit, removing a multiplicity
bias that favored equations with more structural realizations. This can shift selection
among near-ties, within the search's existing unseeded-`randn` nondeterminism
(`fitting.jl:253`). The durable invariant is *"same set of distinct equations, each fit once
at full `n_restarts` budget."*

**Rejected / deferred.**
- **No-improvement termination:** dropped (Denis's call, confirmed by verification). On LDH
  it saves only ~7% (`max_stale=3`) / 13% (`max_stale=2`) because the expensive plateau
  iterations must be *observed* before staleness is detectable, and it overlaps fit-dedup.
  Its value is a robustness bound for pathological runs, not LDH throughput. Record this
  rationale; revisit only if a runaway non-LDH run motivates it.
- **Expansion-dedup by `eq_hash`: rejected.** `expand_mechanisms` is a pure function of
  *structure*, so two same-`eq_hash` mechanisms yield *different* children (e.g. a split
  kinetic group vs a Wegscheider-tied folded group), and deduping there drops reachable
  equations.
- **`min_beam_width=50`: keep.** It amplifies the plateau, not the peak; MWC recovery needs
  ≥10; lowering it does nothing at the peak.

**Risks.** PASS-1 still compiles every structurally-distinct mechanism (redundant *compiles*
remain; compile ≪ fit — a structural compile-memo keyed by `hash(mech)` is a possible
follow-up, YAGNI for now). PASS-2 needs the compiled singleton to fit the representative:
recompile the representative in PASS-2 (simplest; compile ≪ fit) rather than serialize the
large `Sig` type across workers. Post-§5, a representative whose fit *throws* marks all its
`eq_hash` duplicates failed (all-or-nothing per equation) — acceptable, within the
nondeterminism envelope; add a regression test.

**Tests.** Fit-dedup within a batch (two distinct structures, same `eq_hash` → fit invoked
once via a counting stub optimizer; both `BatchEntry`s carry identical loss/retcode; each
carries its *own* rescaled params); cross-batch memo hit (no new fit); rescaling-reuse
soundness (two same-`eq_hash` structures with different full enumerations get correctly
different rescaled params from the same raw fit); an assertion that `best_loss_by_count` and
`cv_pool` equal the non-dedup path for the set of distinct equations. (Note: same `eq_hash`
⟹ identical `fitted_params` names *and order*, since the destructuring line is retained in
the dedup key — the earlier "different tuple order" test scenario is unreachable; test the
real invariant instead.)

### §3 — Misleading progress message

**Root cause.** `:602-605` prints `target n_params=$target`, a counter decoupled from fitted
counts (reached 26 vs a real ceiling of 13).

**Fix.** Report the observed child `n_params` range from `child_entries`:

```julia
child_np = [e.n_params for e in child_entries]
np_label = isempty(child_np) ? "n/a" :
    minimum(child_np) == maximum(child_np) ? string(minimum(child_np)) :
    "$(minimum(child_np))-$(maximum(child_np))"
# "Iteration $iteration (child n_params $np_label): $(length(parents)) parents → …"
```

`child_entries` may be empty when a batch is all-failures (the line is guarded by
`!isempty(child_entries) || !isempty(child_failures)`), so handle empty → `"n/a"`. `target`
stays as the loop-control variable.

**Tests.** Range-label helper (`empty → "n/a"`, single → `"8"`, spread → `"7-13"`); the
emitted line no longer contains `target n_params=`; a beam-search assertion that no progress
line reports `n_params` exceeding `max_param_count`.

### §4 — Equation dedup: keep textual `eq_hash`; canonical key deferred

**Finding.** `_rate_eq_dedup_key` (`:305`) is textual — it strips `#` and
`ANNOTATION_SUBSTITUTED` lines and hashes the rest. It is **correct** for exact re-renders
(verified: stripped key reproduces the stored `eq_hash`); it is **not broken**. Its only
limitation is that it is blind to reparametrized-equivalent equations: different enumerated
mechanisms can derive the same reduced rate *function* yet render with different
independent-K symbols / Wegscheider elimination choices (Canonical Step Form is
per-mechanism, not cross-mechanism), so functional twins get distinct hashes. This has a
real cost — the run's n=8 CV bucket scored ~5 reparametrizations of one function, so some
genuinely-different n=8 models were never cross-validated (loss-proxy twin ratio ~1.8–3.7×
depending on tolerance).

**Decision (Denis, 2026-07-01): defer the canonical key.** A canonical rational-function key
would collapse the twins, but its compute payoff is **unmeasurable in advance** (a loss-based
estimate cannot separate true algebraic twins from data-indistinguishable different
functions — the count runs smoothly from ~1.3× at tol 1e-12 to ~8.8× at 1e-4 with no
plateau), and keying fit-dedup on it would risk a false-merge corrupting selection. It is
also substantial algebraic work intersecting the Canonical Step Form / Haldane guard.
Not worth bundling into this effort. Documented as a known limitation; revisit as its own PR
if selection breadth becomes a priority.

**LOOCV dedup — already present, keep it.** Candidate selection for LOOCV already dedups by
(textual) `eq_hash` in two places: `_offer_cv!` (`:510-524`) keeps the top `n_cv_candidates`
**distinct-`eq_hash`** entries per param count, and the `_cv_model_selection` picker
(`:871-888`) skips repeat hashes when choosing candidates. `§2`'s fit-dedup preserves this
(duplicates still produce `BatchEntry` rows that feed `_offer_cv!`). So no change is needed
for LOOCV dedup; confirm it with a test and add a one-line comment at `_rate_eq_dedup_key`
noting the exact-render-only scope (functional twins may still occupy distinct candidate
slots — the deferred limitation above).

**Tests.** Assert `_offer_cv!` and the `_cv_model_selection` picker keep at most one row per
`eq_hash` per param count (guards the existing LOOCV dedup against regression).

### §6 — Non-Success retcodes (report, not a bug)

310 of 100,787 fits (0.31%) returned non-Success, and every one was `MaxTime` — the 60 s
per-fit budget, concentrated in the high-parameter allosteric iterations where the landscape
is hardest. The fits produced usable parameters; these are not errors. Callers who want more
time can raise the `maxtime` keyword (default 60.0 s, forwarded to `Optimization.solve`).
Add a doctest/kwarg-forwarding test guarding that `maxtime` reaches `Optimization.solve`.

## Testing strategy

TDD per §. Prefer deterministic unit tests over the pure helpers (`_parsimony_cutoff`,
range-label, fit-dedup memo with a counting stub optimizer, LOOCV `eq_hash`-uniqueness) over
full search runs. Add a `DEFINED ⊇ REFERENCED` invariant test for the allosteric body, and a
0-alloc / sub-100 ns perf assertion on a fixed LDH `i_dead` mechanism (the existing perf
test uses PFK-1/HK, which do not exercise the fixed path). Every existing test must stay
green; tests that pin specific LDH winners will change once §5 admits the shed candidates —
update them to reflect corrected behavior, never to paper over it.

## Risks and non-goals

- **§5a/§5b touch the load-bearing derivation.** The runtime perf contract holds because
  both are compile-time; enforce with a perf assertion on the fixed path.
- **§5 admits ~54.7k previously-shed candidates**, so a re-run's selected equation will
  differ from 2026-06-27 (intended). §5 fixing errors *increases* the raw candidate count,
  which is why §2 fit-dedup is included.
- **64-bit `eq_hash` collision** (birthday bound ~4e-11 over ~40k equations) is treated as
  negligible in §2; a collision would silently merge two equations' fits. Not guarded.
- **§4 deferral is a known limitation, not a fix:** functional twins still occupy distinct
  LOOCV candidate slots, so the CV pool can under-sample genuinely-different models at a param
  count (the n=8 "~5 copies" effect persists). Accepted for now.
- **Non-goals:** the canonical rate-function dedup key (§4, deferred), no-improvement
  termination (deferred), expansion-dedup (rejected), a structural compile-memo (YAGNI),
  lowering `min_beam_width`, and re-litigating PR #54's decision to retain `Q_I` for dead
  I-states.

## Implementation phasing

1. **§5a** (UndefVar S_I fix) — highest-value correctness fix; unblocks the shed candidates.
2. **§5b** (kcat/`NaN` common-factor) — derivation root-cause; independent of §5a.
3. **§1** (parsimony helper) — small, isolated.
4. **§3** (progress message) — trivial.
5. **§2** (fit-dedup + `fit_inherited` column) — after §5, since §5 changes the candidate
   volume it operates on; includes the LOOCV `eq_hash`-uniqueness confirming test (§4).
6. **§6** — report + `maxtime` kwarg-forwarding test.

## Open questions

- **§5b:** exact cancellation site in the free-enzyme-per-segment / conc-GCD reduction —
  resolved during implementation; is there overlap to reuse from the
  `rate-equation-no-concentration-division` branch?
- **§5a:** confirm the `_synthesized_dep_i_names` / `parameters(Full)` reconciliation passes
  the `name(p, m)` chokepoint AST-walker test.
- **docs/ldh_hpc_results (2.4 GB) is untracked** — it should not be committed to git as-is;
  confirm whether to keep it local, drop it, or track a small subset.
