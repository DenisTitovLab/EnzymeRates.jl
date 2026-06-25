# Bug: allosteric dead-inactive-state rate equation references an undefined I-state parameter

Status: confirmed and root-caused; not yet fixed.

## Symptom

For certain allosteric (MWC) mechanisms, evaluating the rate equation raises

```
UndefVarError: k_I_EP_to_ES not defined in EnzymeRates
```

and `rate_equation_string` for the same mechanism prints a `v = …` line that uses
`k_I_EP_to_ES` without ever defining it (it is neither in the `params`
destructure nor in any Wegscheider/Haldane constraint). It surfaced while fitting
allosteric candidates during `identify_rate_equation`.

`k_I_EP_to_ES` is the **inactive (I) state** mirror of the reverse catalytic rate
constant `k_EP_to_ES` (the `EP → ES` step). It is a *dependent* parameter (the
I-state mirror of an active-state constant that is itself Haldane-dependent), and
it is supposed to be emitted as a constraint assignment in the rate-equation body
— but in this case it is not.

## Reproduction

The reaction `S → P` with a competitive inhibitor `R`. Of 85 allosteric
mechanisms enumerated (init + 3 expansion rounds), **6 trigger the bug**. They all
share a V-type-like shape — the inactive state binds substrate but cannot turn it
over (the catalytic step is tagged `:OnlyA`):

```
SS  E + P → EP
RE  E + S → ES
SS  ES   → EP      ← catalytic step, inactive-state-dead (:OnlyA)
```

`rate_equation_string` for one of them:

```
(; K_A_S_E, k_ES_to_EP, koff_A_P_E, kon_A_P_E, K_I_S_E, L, Keq, E_total) = params
(; S, P) = concs
# Haldane constraints:
k_EP_to_ES = (1 / Keq) * (1 / K_A_S_E) * k_ES_to_EP * koff_A_P_E * (1 / kon_A_P_E)
v = E_total * (k_ES_to_EP * koff_A_P_E * S / K_A_S_E - k_EP_to_ES * kon_A_P_E * P)
  / ( k_EP_to_ES + koff_A_P_E + … + kon_A_P_E * P
      + L * (k_ES_to_EP * S / K_I_S_E + S * k_I_EP_to_ES / K_I_S_E + k_I_EP_to_ES) )
                                              ^^^^^^^^^^^^^                ^^^^^^^^^^^^^
                                              undefined: in v, but no constraint defines it
```

The active-state reverse constant `k_EP_to_ES` *is* defined (the Haldane line).
Its I-state mirror `k_I_EP_to_ES`, which appears only in the `L * (…)` term, is
not.

## Root cause

The `L * (…)` term is the **inactive-state denominator** `Q_I` — the MWC
"enzyme mass" contribution of the inactive conformation. The bug is a wrong
assumption in the dead-inactive-state optimization in
`src/rate_eq_derivation.jl`.

When the inactive state cannot catalyze (`_i_state_dead`, caused by an `:OnlyA`
catalytic group), `_allosteric_num_den_exprs` (≈ lines 1612–1688):

- zeroes the inactive-state **numerator**: `N_I = 0` (line 1649), so the
  `L * num_I` term is dropped;
- but **keeps** the inactive-state **denominator** `Q_I` (line 1651) and returns
  `den_A + L * den_I` (line 1683), because the inactive conformation still holds
  enzyme even though it cannot turn over. `Q_I` is built by renaming the
  active-state denominator polynomial with `rename_I` (line 1646), which includes
  the synthesized I-state name `k_I_EP_to_ES` (added by `_add_case_b_renames!`,
  line 1631). So `k_I_EP_to_ES` genuinely appears in the kept denominator.

Then both body builders drop the inactive-state assignments entirely when the
inactive state is dead:

- `_build_allosteric_rate_body` (line 1699):
  `i_assignments = _i_state_dead(M()) ? Expr[] : i_assignments_`
- `rate_equation_string` (line 1758): the identical line.

with the justification (comment at ≈ 1696–1698): *"When the I-state cycle is
broken, i_assignments … become dead code — they're only referenced from the
L*num_I branch, which is now elided."*

**That justification is false.** The inactive-state assignments are *also*
referenced by `Q_I` in the kept `L * den_I` denominator term. `_build_dep_assignments`
does correctly build the assignment `k_I_EP_to_ES = (1/Keq) * … * k_I_ES_to_EP`
(lines 1592–1602) — it is then thrown away by the caller. With the assignment
gone and `k_I_EP_to_ES` being a *dependent* (so it is in `merged_dep`, never in
`merged_indep`/the `params` destructure — see `_dependent_param_exprs` line 1390),
the symbol is left undefined in the generated body.

In short: **the dead-inactive-state path zeroes the I-state numerator and elides
all I-state parameter assignments, but keeps the I-state denominator `Q_I`, which
still depends on those assignments.**

### Why it is narrow

All three conditions must hold:

1. The mechanism is allosteric with a **catalytically-dead inactive state**
   (`:OnlyA` catalytic group — the V-type-like case where the inactive
   conformation binds but cannot turn over). This is what elides
   `i_assignments`.
2. The active-state **reverse catalytic** constant is Haldane-dependent (always
   true for a reversible catalytic step constrained by `Keq`), so its I-state
   mirror is a *dependent* `k_I_*` rather than a destructured independent.
3. That I-state dependent **survives into `Q_I`** (the inactive-state
   denominator), so it is referenced despite the numerator being zeroed.

Independent I-state constants (`:NonequalAI`, in the `params` list) and `:EqualAI`
mirrors (which share the active-state symbol and so are destructured) are fine —
only the Haldane-dependent reverse-catalytic I-mirror is left dangling. Hence the
symptom is always this one symbol, `k_I_EP_to_ES`.

## What it is *not*

It is not a missing entry in the parameter list: `k_I_EP_to_ES` is correctly a
dependent parameter, not an independent one, so its absence from the `params`
destructure is correct. Adding it to the independent set would create a phantom
fittable parameter. The defect is purely that its *defining constraint* is elided
while its *use* (in `Q_I`) is retained.

## Fix direction

Do not blanket-elide `i_assignments` when the inactive state is dead. The correct
behavior is to keep exactly the inactive-state assignments that the retained
`Q_I` (denominator) still references, and drop only those used solely by the
removed `L * num_I` numerator. Two reasonable options:

- Stop eliding `i_assignments` in the `i_dead` branch (the assignments are cheap;
  the optimization was eliminating now-genuinely-dead numerator-only bindings, but
  the denominator ones are not dead). Verify no truly-dead binding then references
  a zeroed symbol.
- Or, when `i_dead`, compute the set of symbols actually referenced by `Q_I`
  (`den_I`) and keep only those `i_assignments`.

The same elision appears in two places that must stay in sync —
`_build_allosteric_rate_body` (line 1699) and `rate_equation_string` (line 1758) —
so the fix must touch both (or factor the decision into one shared helper).

A regression test should assert that `rate_equation` evaluates (and
`rate_equation_string` defines every symbol it uses) for a V-type / dead-inactive-state
allosteric mechanism whose reverse catalytic constant is Haldane-dependent — e.g.
the `S → P` (+R) mechanism `{E+P→EP (SS), E+S→ES (RE), ES→EP (SS, :OnlyA)}` above.
