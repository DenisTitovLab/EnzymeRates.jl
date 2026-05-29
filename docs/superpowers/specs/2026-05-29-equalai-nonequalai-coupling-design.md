# EqualAI × NonequalAI Coupling — Contained Fix

**Date:** 2026-05-29 (revised after deeper analysis)
**Branch:** `refactor-to-concrete-types-instead-of-symbols`
**Status:** Approved design, pending implementation
**Companion:** `2026-05-29-direction-symmetry-constraint-resolution.md` (the
larger, post-refactor follow-up this fix consciously defers to — a general,
all-mechanism rewrite of dependent-parameter removal, not just allosteric).

> **Revision note.** An earlier version of this doc proposed rejecting mixed
> EqualAI × NonequalAI configurations via a rank/nullspace validator and
> rewriting the enumeration. That design was based on a **false premise** —
> reading the actual hand-verified mechanisms showed the configurations it
> would reject are physically valid. The corrected analysis and the much
> smaller real fix are below. The rank/symmetry analysis that came out of
> the investigation is preserved in the companion doc for a future PR.

## 1. The real bug

On the structural-parameter-names branch, four allosteric tests fail:
- **PK Constraints** — `n_haldane_constraints = 1` (should be 2),
  `n_mirror_constraints = 0` (should be 4).
- **PK Haldane Equilibrium**.
- **m_mixed p_eq** — rate ≈ 2.615 at chemical equilibrium (should be 0).
- **Allosteric Analytical Rate**.

All four share one root cause in the synth-dep machinery
(`src/rate_eq_derivation.jl`). When an `:EqualAI` catalytic group's
*dependent* parameter (e.g. the catalytic reverse `k2r`/`k5r`, eliminated by
the Haldane) has a dep-expression RHS that references a `:NonequalAI`
symbol, the inactive-state assignment needs a **distinct** synthesized
inactive name (`k2r_T`). But:

```
rename_T[k] = name(_flip_to_inactive(_param_for_symbol(am, k)), am)
```

`_flip_to_inactive` correctly returns an `:EqualAI` parameter **unchanged**
(per spec). So `rename_T[k] == k` — a self-map. Then `dep_T[k]` *overwrites*
the correct `dep_R[k]` in `merge(dep_R, dep_T)`, and the rate-equation body
emits two conflicting assignments for the same symbol → wrong rate.

The old positional naming hid this: it suffixed `_T` to make a new symbol
(`:k2r_T`), so the overwrite never happened. The bug is **downstream code
assuming `_flip_to_inactive`'s result always differs from its input** — not
`_flip_to_inactive` itself, which is correct.

## 2. The corrected model (why these configs are VALID, not rejectable)

The A and I states share the same enzyme-form graph; thermodynamic cycles
pin **ratio-type** quantities per state. Each cycle eliminates one
*dependent* parameter, computed per state. A dependent parameter is **free
to differ between states** — it is derived, not a shared user value.

So an `:EqualAI` SS group whose **forward** rate is shared can still have a
**reverse** rate that legitimately differs between states (the reverse is
the dependent param the Haldane derives). That is exactly PK
(`PEP :NonequalAI`, catalysis `:EqualAI`: `k5r` and `k5r_T` differ, derived
from `K1` and `K1_T`) and `m_mixed` (`:NonequalAI` S-binding + `:EqualAI`
catalysis, expecting rate 0 at equilibrium via a synthesized `k2r_T`). Both
are correct mechanisms the suite *expects to work*.

Therefore: **do not reject these configurations, and do not change the
enumeration.** They are valid. The single defect is that the new
structural-name code computes them wrong. (Full derivation —
dependent-param absorption, the rank/degeneracy test, and the deeper
direction-symmetry principle — is in the companion doc.)

## 3. The fix

In the synth-dep / dep-assignment machinery, when building the
inactive-state rename map: if a dependent parameter `k` is `:EqualAI` (so
`_flip_to_inactive` returns it unchanged) **but its dep-expression RHS
references at least one `:NonequalAI` symbol**, synthesize a **distinct
inactive-state name** for `k` instead of self-mapping. This restores the
old positional code's behavior (a separate `k_..._T` symbol) without
touching `_flip_to_inactive`'s semantics.

- Apply uniformly across the synth-dep call sites that build `rename_T`
  (`_dependent_param_exprs(::AllostericEnzymeMechanism)`,
  `_build_dep_assignments`, `_synthesized_dep_t_names`, and the
  `rate_equation_string` path).
- **Verify Haldane *and* Wegscheider coverage.** The same overwrite can
  occur for a Wegscheider-derived dependent binding `K` (random-order
  loops), not only the Haldane-derived reverse rate. Confirm whether the
  current `rename_T` loop already walks Wegscheider deps; if it only handles
  Haldane, extend the fix to Wegscheider deps too. (Establish by inspecting
  `_dependent_param_exprs_allosteric` and the dep-set the `rename_T` loop
  iterates.)

## 4. Tests (true TDD — pre-derive expected values independently)

- **PK Constraints / Haldane Equilibrium:** with the fix, PK must show
  `n_haldane = 2`, `n_mirror = 4`, and zero net rate at equilibrium.
  Re-derive these counts from the dep machinery on the corrected PK output,
  matching implementation to truth — no blind golden edits.
- **m_mixed:** the existing `@test isapprox(rate_eq, 0.0; atol=1e-10)` must
  pass. **Keep m_mixed as-is** (it is a valid single-NonequalAI mechanism);
  the earlier plan to convert it to a `@test_throws` or all-NonequalAI was
  wrong.
- **Allosteric Analytical Rate:** passes once the synth-dep emits the
  correct distinct inactive assignment.
- **Wegscheider regression:** add a random-order mechanism where a
  Wegscheider-derived dependent `K` is `:EqualAI` with a `:NonequalAI`
  partner in the loop; assert correct (zero-at-equilibrium / matching
  analytical) behavior. Pre-derive the expected dependent-`K` per-state
  values with independent linear algebra, not the code under test.
- **No test deletion**, no `@testset` / `MECHANISM_TEST_SPECS` removal.
- Leave the deferred **`canonical-hash partition stability`** failure
  (21 vs 23) alone — parent session's territory.

## 5. Explicitly NOT in this fix

- **No rank/nullspace validator** and **no constructor rejection** — the
  configs are valid.
- **No enumeration change** (`_expand_change_allo_state` etc.) — single
  NonequalAI flips are valid (absorbed by the dependent reverse).
- **No representation change**, **no `√`**, **no symmetric resolution** —
  all deferred to the companion-doc follow-up PR.
- Do **not** touch the structural-naming chokepoint (`_state_tag`,
  `_render_binding`, `_render_iso`) or `_flip_to_inactive`'s semantics.
- The genuinely-degenerate **lone-NonequalAI pure-RE Wegscheider** case
  (an over-parametrization, not a crash) is left for the follow-up — see
  `2026-05-29-nonequalai-rank-validity.md`; it does not block suite-green.
- **`parameters(m, Full)` has a pre-existing duplicate** for Case-B shapes: the
  base enumeration (`_all_i_state_parameters`) already emits the I-state reverse
  (e.g. `:k_I_EP_to_ES`) as a `Krev(rep, :I)` struct, and `_synthesized_dep_t_names`
  emits the same name — so `parameters(Full)` is not `allunique`. This is not
  introduced by this fix (the self-map merely changed which symbol duplicated),
  only affects `parameters(Full)` whose consumer (`identify_rate_equation`
  canonicalizer) is still pending, and is entangled with the parent refactor's
  parameters-API design (whether `_all_i_state_parameters` should emit an
  I-mirror for an `:EqualAI` group whose reverse is a dependent Haldane
  constant). Deferred to the parent refactor / a follow-up; not chased here.

## 6. Done when

- The 4 failing allosteric tests pass; full suite green except the deferred
  hash-partition failure.
- No new failures in unrelated testsets.
- CLAUDE.md "Allosteric state taxonomy" gains a short note: an `:EqualAI`
  group's *dependent* parameter (e.g. SS reverse rate) may legitimately
  differ between A/I states when the Haldane/Wegscheider RHS references a
  `:NonequalAI` symbol; the synth-dep synthesizes a distinct inactive name
  for it. Points to the companion doc for the planned symmetric resolution.
