# Allosteric multi-substrate derivation errors (LDH run) and a combined fix plan

Status: investigated and root-caused; fix proposed, not implemented.

## Context

A real `identify_rate_equation` run on LDH (`NADH + Pyruvate ⇌ Lactate + NAD`,
`oligomeric_state = 4`, `scale_k_to_kcat = 1.0`) produces a large allosteric
(MWC) search space and a large number of **error rows** in the result CSVs.
Counting the error messages across `equation_search_iteration_*.csv`, two
distinct families dominate, each affecting tens of thousands of candidate fits
(the search reaches >1M rows by iteration 6, of which ~140k carry an error):

- **Family A — `UndefVarError: k_I_… / K_I_… / koff_I_…` not defined in EnzymeRates**
  (e.g. `k_I_ELactateNAD_to_ENADHPyruvate`, `K_I_Pyruvate_E`, `koff_I_Pyruvate_E`).
- **Family B — `_kcat_forward: AllostericEnzymeMechanism with multiple
  saturating-substrate kcat components (N found) is unsupported`** (observed
  `N = 2,3,5,6,7,8,12,16`).

Both fire only for **allosteric, multi-substrate / oligomeric** mechanisms — the
class that a bi-bi homotetramer generates in bulk but that the simpler reactions
in the test suite (mostly uni-uni MWC) never exercise. A shallow local
enumeration (`bi-bi`, `oligomeric_state = 2`, first 150 allosteric mechanisms)
triggers neither, consistent with the errors appearing only at search depth 3+ in
the LDH run.

## Family A: undefined I-state parameters (the `k_I_EP_to_ES` bug, at scale)

This is the **same root cause** documented in
[`allosteric_dead_istate_undefined_param_bug.md`](allosteric_dead_istate_undefined_param_bug.md).
Briefly: when the inactive state cannot catalyze (`_i_state_dead`, caused by an
`:OnlyA` catalytic group — the V-type-like case), the rate-equation builders zero
the inactive-state **numerator** but keep the inactive-state **denominator**
`Q_I` ("enzyme mass"), then **elide all inactive-state dependent-parameter
assignments** on the false assumption that they are only used by the now-removed
numerator. `Q_I` still references those I-state dependents, so the body uses them
undefined.

On a single-substrate reaction this surfaced as exactly one symbol
(`k_I_EP_to_ES`, the I-mirror of the Haldane-dependent reverse catalytic
constant). On LDH's **bi-substrate** reaction the inactive-state denominator has
many more I-state dependents — binding constants (`K_I_Pyruvate_E`), off-rates
(`koff_I_Pyruvate_E`), and the catalytic reverse (`k_I_ELactateNAD_to_ENADHPyruvate`)
— so the same bug manifests as many undefined symbols.

The faulty `i_dead ? Expr[] : i_assignments_` elision appears in **three**
functions in `src/rate_eq_derivation.jl`:

- `_build_allosteric_rate_body` — line ~1699 (the runtime `rate_equation` crash);
- `rate_equation_string` — line ~1758 (the string shows the undefined symbol);
- `_kcat_forward` (allosteric) — line ~933 (paired with zeroing `A_I`/`B_I`; needs
  review — when `i_state_dead` but the saturating pattern is *present* in the
  inactive state, `B_I = den_k_I^CatN` is non-zero and can reference an elided
  I-dependent, the same defect).

## Family B: `_kcat_forward` multiple saturating-substrate components

`_kcat_forward` computes the forward turnover number (used only when
`scale_k_to_kcat` is set — the LDH run sets it to `1.0`) by grouping the rate
equation's numerator/denominator by metabolite "saturation pattern" and reading
off `num_k / den_k` at the saturating pattern.

The allosteric implementation (`src/rate_eq_derivation.jl`, `_kcat_forward` for
`AllostericEnzymeMechanism`, lines ~833–930) **assumes exactly one** saturating
pattern:

```julia
a_keys = sort!([k for k in keys(num_A_groups) if haskey(den_A_groups, k)])
…
length(a_keys) == 1 ||
    error("_kcat_forward: AllostericEnzymeMechanism with multiple " *
          "saturating-substrate kcat components ($(length(a_keys)) found) " *
          "is unsupported")
```

with the comment "single component for mechanisms exercised here; assert keeps
that constraint visible." Multi-substrate oligomeric MWC mechanisms routinely
produce more than one saturating-substrate pattern (alternative catalytic routes,
and the way the numerator factorizes with `catalytic_multiplicity > 1`), so the
assert fires. This is a **deliberately-unsupported case**, not a latent crash —
but it kills the candidate's fit just the same.

Note the contrast with the **non-allosteric** `_kcat_forward` (same file, lines
~788–818), which already handles multiple candidates correctly:

```julia
candidates = [:($nk / $dk) for (nk, dk) in components]
result = length(candidates) == 1 ? candidates[1] : Expr(:call, :max, candidates...)
```

i.e. kcat = max forward turnover over the alternative saturating paths — exactly
the generalization the allosteric path is missing.

## Are A and B related?

**They are two distinct bugs, not one.** Different functions, different root
causes: A is a wrong dead-code elimination in the rate-equation body (a
correctness bug that crashes valid mechanisms); B is an unimplemented case in the
kcat extraction (a feature gap that refuses to fit otherwise-valid mechanisms).

They share a **theme and a trigger class**: the MWC allosteric derivation was
built and tested on uni-uni mechanisms (single substrate, single saturating
component, simple inactive state), and both gaps surface specifically on the
richer **multi-substrate / oligomeric** mechanisms that the LDH bi-bi homotetramer
generates. They also share an incidental code smell: the same `i_dead` elision
pattern (Family A) is copy-pasted across three functions, one of which is the very
`_kcat_forward` that carries Family B. Fixing them together is natural because
they live in the same file and the same allosteric code, and because both must be
fixed before the LDH search can complete without shedding most of its allosteric
candidates.

## Combined fix plan

Reproduce first (TDD), then fix both in one PR in `src/rate_eq_derivation.jl`,
each with a regression test built from a minimal multi-substrate / dead-I-state
allosteric mechanism.

### Step 0 — reproductions / failing tests

- **A:** the uni-uni mechanism `{E+P→EP (SS), E+S→ES (RE), ES→EP (SS, :OnlyA)}`
  (S→P, +R) already reproduces the undefined I-state dependent. Add a bi-substrate
  dead-I-state mechanism too. Assertion: `rate_equation` evaluates and
  `rate_equation_string` defines every symbol it uses.
- **B:** find a bi-bi `oligomeric_state ≥ 2` mechanism whose `_kcat_forward`
  produces `length(a_keys) > 1` (search the deeper enumeration; the LDH run shows
  they are common at depth 3+). Assertion: `_kcat_forward` returns a finite value.

### Step 1 — Fix A (undefined I-state dependents)

Replace the blanket `i_dead ? Expr[] : i_assignments_` elision with "keep exactly
the inactive-state dependent assignments that the retained `Q_I` (denominator)
still references." Concretely:

1. Factor the decision into one shared helper (the elision is currently duplicated
   in `_build_allosteric_rate_body`, `rate_equation_string`, and `_kcat_forward` —
   they must not drift).
2. When `i_dead`, compute the set of symbols referenced by the kept inactive-state
   denominator expression and keep only those `i_assignments` (drop the genuinely
   numerator-only ones). Simplest correct version: stop eliding entirely and rely
   on the assignments being harmless when unused — but verify no assignment then
   references a zeroed/`:OnlyA` symbol (those are already set to `0` at
   lines ~1596–1597, so this is safe).
3. Re-check `_kcat_forward`'s `B_I` branch: when `i_state_dead` and the saturating
   pattern survives in the inactive state, `B_I` references I-state deps and needs
   the same kept assignments.

### Step 2 — Fix B (multiple saturating-substrate kcat components)

Generalize the allosteric `_kcat_forward` to mirror the non-allosteric one: build
one kcat candidate per saturating pattern `met_key ∈ a_keys` (each using the
existing per-pattern `A/B` construction at lines ~887–930), and return the `max`
over patterns (and, as today, over regulator corners). Semantics: kcat = maximum
forward turnover at saturation across alternative catalytic routes, consistent
with the non-allosteric docstring ("returns the max over all paths"). Remove the
`length(a_keys) == 1 || error(...)` assertion (keep the genuinely-empty-`a_keys`
error at lines ~880–882).

### Risks / verification

- **Performance contract.** `rate_equation` must stay 0-alloc / sub-100ns for the
  non-allosteric `MECHANISM_TEST_SPECS`. Fix A only adds assignments on
  *allosteric, dead-I-state* mechanisms; confirm allosteric mechanisms are not part
  of that perf contract (they are larger and not in `MECHANISM_TEST_SPECS`), and
  that the non-allosteric body is untouched.
- **Validation.** After both fixes, re-run a bounded LDH-like search (e.g. bi-bi
  `oligomeric_state = 2` or `4`, low `max_param_count`) and confirm the two error
  families no longer appear in the result CSVs. Spot-check that the previously
  erroring mechanisms now produce finite losses and sensible kcat values.
- **Scope check.** These mechanisms *should* be fittable (they are valid MWC
  mechanisms), so fixing the derivation is correct. An alternative — not generating
  some of them — is a separate enumeration-design question and would also interact
  with the redundancy/identifiability discussion in the loss-dominance note.
