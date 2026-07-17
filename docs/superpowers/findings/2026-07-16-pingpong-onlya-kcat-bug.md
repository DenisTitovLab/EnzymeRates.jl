# Ping-pong `:OnlyA` `_kcat_forward` bug — analysis

## Symptom

A PFK-P HPC run errored on 112 of 691 initial mechanisms, every one with the same
message:

```
_kcat_forward: AllostericEnzymeMechanism produced no kcat components —
saturating-substrate pattern not found in numerator
```

Iteration 1 (5936 children) errored on **zero**. All 112 failures are at the base
tier. This is *after* the `:OnlyA` thermodynamic guard (PR #70) landed — that guard
is working (the original run errored 118 of 359; the expansion tier is now clean),
but a distinct class survives it.

## What all 112 have in common

- **All are ping-pong** — a covalent enzyme intermediate (`E(; residual = …)`), so
  the mechanism graph has **two empty-bound forms**: free `E` and the covalent
  intermediate.
- **All are accepted by the guard** (`_onlya_haldane_violation` returns `nothing`),
  and correctly so — see "The mechanisms are valid" below.

By `:OnlyA` placement:

| class | count | `:OnlyA` on |
|---|---|---|
| V-system | 14 | one chemical (half-reaction) step, no binding |
| balanced K-system | 14 | one substrate binding **and** one product binding |
| one-sided + chem | 84 | one chemical step + one one-sided binding (42 substrate, 42 product) |

## The mechanisms are valid, not traps

The two one-sided cases were checked against an independent n=1 two-conformation
mass-action ground truth (formulation 1: only free `E` flips between the A and I
conformations; the covalent intermediate does not flip). Both give `v = 0` at the
equilibrium metabolite ratio (`1.7e-16`, `7.6e-17`) and `v ≠ 0` off-equilibrium —
they are thermodynamically consistent and genuinely catalyze.

The earlier "irreversible-feed trap" hypothesis was **wrong**. The 4th field of a
`Step` is `is_equilibrium` (rapid-equilibrium vs steady-state), **not**
reversibility — every step is microscopically reversible, so the reverse of a
"feed" step simply holds equilibrium mass with zero net flux. There is no trap.
The V-system (14) and balanced K-system (14) classes are valid by inspection; the
guard is right to accept all of them.

## Root cause: ping-pong breaks the free-enzyme `d_free` normalization

The allosteric combine normalizes each conformation to its free-enzyme
spanning-tree weight `d_free`. That machinery assumes **one** free-enzyme form. A
ping-pong has **two** empty-bound forms, and `_reachable_from_free` seeds both as
"free," so the single-free-enzyme normalization is ill-defined. It fails two ways:

- **`d_free_I = 0`** — the inactive graph fragments into an isolated, zero-mass
  covalent island (e.g. `{E_cov_I, E(F16BP)_cov_I}`). Dividing by it gives `NaN`;
  `rate_equation` returns `NaN`. (56 of the 84 one-sided cases return `NaN`.)
- **`d_free_I` carries a product** — the cross-weight then multiplies the A-state
  numerator by that product concentration, so every saturating-substrate group key
  acquires a product factor. `_kcat_forward` evaluates kcat at products = 0 and
  filters out any product-bearing key, so `a_keys` empties and it throws
  "no kcat components." (The remaining cases, and the kcat crash generally.) This
  is how the throw arises mechanically, but it does not mean `_kcat_forward` is
  wrong — see "Options" below.

These are not the same defect. `d_free_I = 0` (err1) is a genuine normalization
bug: dividing by a zero-mass reference form that should never have been treated as
free. A product-bearing `d_free_I` (err2) is not a bug — it correctly reflects a
real absorbing trap in the mechanism, detailed under "Options" below. Neither is a
thermodynamic-validity problem or a guard gap; the guard correctly accepts all of
these valid mechanisms.

## Relationship to the guard and the redesign

`_onlya_haldane_violation` is a Haldane-*satisfiability* check, and it is **complete
for that contract** here — not leaking. In both mechanisms the sole catalytic cycle
carries the F6P `:OnlyA` binding *and* the `:OnlyA` chemical step, so the cycle is
genuinely satisfiable (the inactive chemical step's rate-constant ratio
`k_I_f/k_I_r` vanishes from the rate law and is free to absorb F6P's affinity
divergence). Dropping the `:OnlyA` chemical edge leaves a spanning tree — zero
constraint rows — so the guard's sign test has nothing to inspect and returns
`nothing`, which is the *correct* verdict. These are valid mechanisms; tightening
the guard cannot and should not catch them.

This diagnosis originally motivated proposing the solve-then-limit redesign as
the fix — see
`docs/superpowers/specs/2026-07-16-mwc-solve-then-limit-derivation-design.md` /
`docs/superpowers/plans/2026-07-16-mwc-solve-then-limit-derivation.md`. That
redesign was **declined**; both failure modes were instead resolved (or found
not to need resolving) locally. See "Options" below for the measured outcome
and `docs/superpowers/specs/2026-07-16-mwc-derivation-targeted-fixes-design.md`
for the full measurements behind declining the redesign.

err1's `d_free_I = 0` annihilates the A-numerator (`×0^n`) → `rate_equation` is
`NaN` structurally (all draws, not just at equilibrium); err2's product-bearing
`d_free_I` keeps `rate_equation` finite (the injected factor cancels between
numerator and denominator) but still empties the kcat `a_keys`. Two different
causes, two different symptoms — see "Options" below.

## Options — measured outcome

The three options above (band-aid sentinel / derivation-time degeneracy check /
solve-then-limit redesign) were superseded once err1 and err2 were actually
measured on `mwc-targeted-fixes`. Neither needed a new degeneracy check or a
sentinel; one was a real bug and got fixed, the other was not a bug.

**Option 3 (redesign) — declined.** The solve-then-limit rewrite proposed in
`docs/superpowers/specs/2026-07-16-mwc-solve-then-limit-derivation-design.md`
was independently reproduced claim-by-claim and declined; the full measurement
lives in
`docs/superpowers/specs/2026-07-16-mwc-derivation-targeted-fixes-design.md`.
The `:OnlyA` limit equals graph deletion on every constructable mechanism, so
the rewrite could not have changed either failure mode here — both were
resolved (or found not to need resolving) with local fixes instead.

**err1 (`d_free_I = 0` / NaN) — fixed.** `_reachable_from_free`'s seed treated
every empty-`bound` form as free, so a ping-pong's covalent intermediate (empty
`bound`, non-empty `residual`) seeded itself as a second free root. Under
formulation 1 only free enzyme flips, so a component free `E` cannot reach
holds no inactive mass and must be stranded — the seed must also require an
empty `residual`. Tightening it (commit `29eee7e`) fixes it: `d_free_I` 0→1,
`rate_equation` NaN→0.0966, `_kcat_forward` crash→ok. No regression: 1803/1803
derivation tests pass, golden reference byte-identical.

**err2 (`_kcat_forward` crash) — not a bug; this doc's root cause for it was
wrong.** In the inactive conformation the F6P binding and the
`E(F6P)→E(ADP)` iso step are `:OnlyA` and therefore deleted, so at `F16BP = 0`
the form `E_res_I` has no route out: it cannot flip (formulation 1 flips only
free `E`; `E_res` carries a residual) and it cannot react. It is an
**absorbing trap** — all enzyme drains into it and `v = 0` genuinely, not a
normalization artifact. `_kcat_forward` evaluates forward turnover at
products = 0, and the true forward turnover there really is zero, so
"no kcat components" is a **true report**. There is nothing to fix in
`_kcat_forward`.

Supporting measurement: `v → 0` as `F16BP → 0` for any fixed `L > 0` (2.33e-7
at `F16BP = 1e-6`; 2.33e-10 at `F16BP = 1e-9`), while `v = 0.331` at `L = 0`.
`v` is smooth in `L` at fixed `F16BP` — this is *not* a discontinuity in `L`;
the `F16BP → 0` and `L → 0` limits do not commute.

Option 2's obvious local fix — making `_kcat_forward` group on un-normalized
polynomials — was tried and **breaks an existing gate**:
`test/test_rate_eq_derivation.jl:2209` ("Fix B", live I-state) goes to 7.7%
kcat error (0.2457 vs. grid-peak 0.2661). The free-enzyme normalization is a
common factor *within* a conformation but not across the `L`-weighted A/I
combine, so stripping it from the group key changes which polynomial terms
match, and the fix is not clean.

**Open question for the repo owner** (not resolved here): whether err2-class
mechanisms should be (a) accepted with a non-fatal kcat sentinel, (b) rejected
at construction as degenerate, or (c) left as-is, since the behavior is
correct and the HPC crash is enumeration meeting a genuinely degenerate
mechanism.

## Evidence

Scripts under `docs/superpowers/prototypes/2026-07-16-mwc-solve-then-limit/`:
`RC_gt.jl` (the ground-truth validity check), `RC_probe.jl` (the two-free-form
graph fragmentation), `RC_dfree.jl` (the `d_free_I` = 0 / product-bearing values).
The 112-mechanism categorization and per-class `rate_equation` finiteness were
computed directly from the run's `initial_mechanisms.csv`.
