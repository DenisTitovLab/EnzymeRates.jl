# Allosteric multi-substrate fixes: implementation spec

Status: design approved; reproductions confirmed; kcat semantics verified; ready
for implementation (TDD).

This spec turns the two investigation notes into an implementation plan and adds a
third fix that surfaced during reproduction. Read first:

- [`allosteric_dead_istate_undefined_param_bug.md`](../../design_notes/allosteric_dead_istate_undefined_param_bug.md) — Family A root cause.
- [`allosteric_multisubstrate_errors_and_fix_plan.md`](../../design_notes/allosteric_multisubstrate_errors_and_fix_plan.md) — Families A and B at LDH scale.

## Scope: three fixes, one PR

The LDH run (`NADH + Pyruvate ⇌ Lactate + NAD`, `oligomeric_state = 4`,
`scale_k_to_kcat = 1.0`) sheds ~140k allosteric candidates to two error families.
Both are real bugs in valid MWC mechanisms. A third fix makes such failures
debuggable.

- **Fix A** — dead-inactive-state rate equations reference undefined I-state
  parameters (correctness crash).
- **Fix B** — `_kcat_forward` rejects mechanisms with more than one
  saturating-substrate pattern (feature gap).
- **Fix C** — CSV failure rows discard the failing mechanism's structure
  (debuggability).

A and B are distinct bugs that share a trigger class: multi-substrate / oligomeric
allosteric mechanisms, which the uni-uni test suite never exercises. C is
independent; it would have removed the friction we hit while reproducing A.

## Reproductions (confirmed)

Both families reproduce from the run's own `docs/Cluster_results/2026_06_23_results`
CSVs, using the run's loss filter (keep mechanisms within `1.3 ×` the minimum loss
at a tier), one round of `expand_mechanisms`, and the allosteric children. Each
mechanism round-trips from its `mechanism_type` Sig string via
`Mechanism(em)` / `AllostericMechanism(em)`, and the bug fires from structure
alone — no fit, no real parameter values.

- **Family B**: filter `equation_search_iteration_1.csv` (22 parents), expand
  (254 children, 133 allosteric), call `_kcat_forward`. 54 children trigger the
  assert, with `a_keys` counts of 2 and 3. Smallest: the LDH homotetramer,
  7 fitted parameters.
- **Family A**: filter `equation_search_iteration_2.csv`, expand, keep allosteric
  dead-inactive-state children, build `rate_equation`. A child with 8 fitted
  parameters leaves `k_I_ELactateNAD_to_ENADHPyruvate` undefined — the exact
  symbol the LDH run reported. Its catalytic step carries `:OnlyA`, so the
  inactive state is dead.

Family A first appears at search depth 3 and Family B at depth 2, which is why a
single expansion of the iteration-1 (B) and iteration-2 (A) survivors reaches
each. Error rows in the CSVs are *not* round-trippable: `_failure_row` records the
bare `"EnzymeRates.AllostericMechanism"` and an empty `rate_equation`. Fix C
removes that limitation.

## Fix A — keep inactive-state assignments

`_allosteric_num_den_exprs` zeroes the inactive-state numerator when the inactive
state is dead but keeps the inactive-state denominator `Q_I` (the `L * den_I`
term, `rate_eq_derivation.jl:1679–1683`). `Q_I` references I-state dependent
parameters. Three callers then drop every I-state assignment:

- `_build_allosteric_rate_body` — line 1699.
- `rate_equation_string` — line 1758.
- `_kcat_forward` — line 933 (combined with `i_pattern_dead`).

**Change:** stop eliding. Replace each `_i_state_dead(...) ? Expr[] : i_assignments_`
with `i_assignments_`.

**Why it is safe:** `_build_dep_assignments` already assigns `0` to any I-state
dependent whose right-hand side references an `:OnlyA` symbol
(`rate_eq_derivation.jl:1596–1597`). Keeping the assignments therefore defines
every symbol `Q_I` uses, with the same values the live-inactive-state path already
emits. The reverse-catalytic mirror (`k_I_EP_to_ES`, `k_I_ELactateNAD_to_ENADHPyruvate`)
appears in `Q_I` only as an added term, never as a divisor, so its `0` value
introduces no `Inf`.

The justifying comments (lines 1696–1698 and 931–932) become false once the
elision is gone — they claim the assignments are dead. Remove them.

The `_kcat_forward` site (line 933) keeps its separate `i_pattern_dead` handling,
which Fix B restructures (below).

## Fix B — max kcat over all saturating patterns

`_kcat_forward` for `AllostericEnzymeMechanism` picks `met_key = a_keys[1]` and
asserts a single saturating pattern (`rate_eq_derivation.jl:883–887`):

```julia
length(a_keys) == 1 || error("…multiple saturating-substrate kcat components…")
```

Multi-substrate oligomeric MWC mechanisms routinely produce several. The
non-allosteric `_kcat_forward` already handles this (lines 788–818): it builds one
candidate per pattern and returns `max(...)`.

**Change:** loop the per-pattern construction (the `met_key`-dependent block,
lines 887–1016) over every `met_key ∈ a_keys`. Collect each pattern's
regulator-corner expressions, then return `max` across patterns and corners.
Keep the empty-`a_keys` error (lines 880–882); remove the `length(a_keys) == 1`
assert. `i_pattern_dead` stays per-pattern inside the loop, and Fix A keeps the
assignments each pattern needs.

**`max` is the correct operation, not a heuristic.** `_kcat_forward` returns the
*peak achievable forward turnover* — the supremum of `v/E_total` over all substrate
levels, substrate ratios, and effector states, at products = 0. Each matched
monomial `n_m/d_m` is the turnover in one saturation regime; `max` selects the
best. This was verified end-to-end on a random-order steady-state bi-bi (which has
`A·B, A²·B, A·B²` in the numerator): `_kcat_forward` equalled the numerical grid
peak of the forward rate to 10 digits. The `max` also self-corrects two ways:
substrate-inhibition monomials are denominator-only and never match; effector
inhibition is handled by the regulator corner-`max` (it picks the inhibitor-off
corner — confirmed to equal the `K_reg→∞` neutralized value on a V-type PFK case).
**Product-containing matched monomials must be excluded** from the candidate set:
`_kcat_forward` is evaluated at products = 0, so a met_key carrying a product is
outside the evaluation domain. For non-allosteric `CatN = 1` they happen never to
win, but for allosteric oligomers (`CatN > 1`) the per-pattern closed form can make
a product-containing candidate exceed the true peak (observed 1.30× overestimate),
so the allosteric `_kcat_forward` filters `a_keys` to substrate-only patterns. The
non-allosteric `_kcat_forward` already does the `max` (lines 788–818); Fix B brings
the allosteric path to parity and adds the products-filter.

**Docstring.** State this contract on the allosteric `_kcat_forward`: kcat is the
peak productive forward turnover (max over saturating patterns and regulator
corners); denominator-only substrate-inhibition regimes and product-present
regimes do not inflate it.

## Fix C — record the failing mechanism

`_failure_row` (`identify_rate_equation.jl:380`) records
`mechanism_type = string(typeof(f.mech))`, which collapses a concrete
`AllostericMechanism` to its bare type name and loses all structure.

**Change:** record the round-trippable parametric Sig string, degrading
gracefully:

```julia
mechanism_type = try
    string(typeof(compile_mechanism(f.mech)))
catch
    string(typeof(f.mech))   # compile itself failed — keep today's behavior
end,
```

Leave `rate_equation = missing`: `rate_equation_string` can itself raise, and the
Sig string already makes the mechanism reproducible. No `FitFailure` change —
`f.mech` is already the concrete mechanism.

## Considered and rejected: role-aware regulator corners

We examined whether the corner-`max` mishandles a species that is both a substrate
and an allosteric inhibitor (e.g. PFK / ATP). For a V-type inhibitor the strict
`[S]→∞` limit is the self-inhibited rate, which the corner-`max` does *not* return
— it returns the inhibitor-neutralized value (verified: `max` = `K_reg→∞` to 4
digits). That is **correct** for kcat-as-normalization-anchor: kcat is the peak
productive turnover, not the self-inhibited tail. So no "pin substrate-ligand
corners on / product-ligand corners off" rework is needed. The only untreated edge
is a species that is both a **product** and an allosteric **activator** (`max`
would switch it on though products = 0); it is rare and left as a documented
limitation, not fixed here.

Anchoring kcat to a real *measured* value (rather than the analytic peak) is a
separate enhancement — see
[`kcat_rescale_at_measured_concentrations.md`](../../design_notes/kcat_rescale_at_measured_concentrations.md).

## Tests (TDD: write failing tests first)

Add to `test/test_rate_eq_derivation.jl` (A, B) and a fitting/identify test for C.

- **A, uni-uni (minimal):** `{E+P→EP (SS), E+S→ES (RE), ES→EP (SS,:OnlyA)}` for
  `S → P` with competitive inhibitor `R`. Assert `rate_equation` evaluates and
  every symbol in `rate_equation_string` is defined (destructure or constraint
  left-hand side).
- **A, bi-substrate:** a hand-written `@allosteric_mechanism` over a two-substrate,
  two-product reaction with an `:OnlyA` ternary→ternary catalytic step
  (`oligomeric_state = 2` suffices), mirroring the reproduced LDH structure. Same
  assertion; confirms the multi-symbol case.
- **B:** a bi-substrate `oligomeric_state ≥ 2` MWC mechanism whose `_kcat_forward`
  yields `length(a_keys) > 1`. Assert `_kcat_forward` returns a finite value **and
  equals the numerical peak of `rate_equation` over a substrate grid at products =
  0** (the peak-productive-turnover contract), not just that it runs.
- **B, kcat semantics (non-allosteric guard):** a random-order steady-state bi-bi
  (numerator `A·B, A²·B, A·B²`). Assert `_kcat_forward` equals the numerical
  grid-peak forward rate. Guards the `max` contract that Fix B mirrors.
- **C:** force a `FitFailure` and assert its `_failure_row.mechanism_type`
  round-trips via `Core.eval` + `Mechanism`/`AllostericMechanism` to a mechanism
  equal to the original.

The reproduced fixtures' round-trippable Sig strings are saved under the session
scratchpad for reference when authoring the hand-written test mechanisms.

## Guardrails

- **Performance contract.** `rate_equation` must stay 0-allocation and sub-100ns
  for every mechanism in `MECHANISM_TEST_SPECS`, including the MWC dimer
  (spec 24B). Fix A changes only dead-inactive-state bodies. Confirm spec 24B is
  not dead-inactive-state (classical MWC keeps both states catalytic, so Fix A
  leaves it untouched), and run `test_performance` regardless. Fixes B and C never
  touch `rate_equation`.
- **Validation.** Re-run a bounded LDH-like search (bi-bi, `oligomeric_state = 2`
  or `4`, low `max_param_count`) and confirm both error families vanish from the
  result CSVs. Spot-check that previously erroring mechanisms now produce finite
  losses and kcat values.
- **Full suite.** Run `Pkg.test()` before committing.
