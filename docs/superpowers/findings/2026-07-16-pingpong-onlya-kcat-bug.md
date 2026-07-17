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
  "no kcat components." (The remaining cases, and the kcat crash generally.)

Both are the same defect — the free-enzyme normalization has no well-defined
single reference form in a ping-pong. It is *not* a thermodynamic-validity problem
and *not* a guard gap; the guard correctly accepts these valid mechanisms, and the
derivation then mishandles them.

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

The failure is entirely downstream, in the free-enzyme normalization — which is
exactly what the solve-then-limit redesign rebuilds (per-state normalization done
symbolically, ping-pong bi-bi already named as the plan's normalization test case).
Spec/plan: `docs/superpowers/specs/2026-07-16-mwc-solve-then-limit-derivation-design.md`,
`docs/superpowers/plans/2026-07-16-mwc-solve-then-limit-derivation.md`.

Note err1 and err2 fail *differently* from the one root cause: err1's `d_free_I = 0`
annihilates the A-numerator (`×0^n`) → `rate_equation` is `NaN` structurally (all
draws, not just at equilibrium); err2's product-bearing `d_free_I` keeps
`rate_equation` finite (the injected factor cancels between numerator and
denominator) but still empties the kcat `a_keys`. Same defect, two symptoms.

## Options

1. **Immediate HPC unblock (band-aid):** make `_kcat_forward` return a sentinel
   ("no forward kcat") instead of throwing. The 112 stop crashing — valid
   mechanisms rescale approximately (only the reported parameter scale is off, not
   the fitted loss or model selection); the `NaN`-`rate_equation` cases fail to fit
   gracefully and lose in the beam. The run completes. Does not fix correctness.
2. **Derivation-time degeneracy check (provably-complete crash-stopper):** the
   exact defect condition is `d_free_I` not being a metabolite-free monomial (`0`
   or metabolite-bearing). Detecting that where `d_free_I` is computed and
   rejecting/skipping the mechanism is provably complete — unlike an
   enumeration-time structural pre-filter keyed on `:OnlyA` *bindings*, which is
   strictly narrower (a mechanism with an `:OnlyA` chemical step closing the sole
   cycle and all-`:EqualAI` bindings hits the same degenerate `d_free_I` with no
   `:OnlyA` binding to key on). Still *discards* valid mechanisms rather than
   deriving them, but it stops the crash without mislabeling the cause.
3. **Redesign (the real fix):** the solve-then-limit derivation removes the
   fragmenting single-free-form normalization; ping-pong normalizes correctly by
   construction, so these mechanisms *derive* instead of being skipped.

Not recommended: a guard extension or an enumeration-time `:OnlyA`-binding
pre-filter — both are narrower than the actual `d_free_I` degeneracy condition and
would reject valid mechanisms.

## Evidence

Scripts under `docs/superpowers/prototypes/2026-07-16-mwc-solve-then-limit/`:
`RC_gt.jl` (the ground-truth validity check), `RC_probe.jl` (the two-free-form
graph fragmentation), `RC_dfree.jl` (the `d_free_I` = 0 / product-bearing values).
The 112-mechanism categorization and per-class `rate_equation` finiteness were
computed directly from the run's `initial_mechanisms.csv`.
