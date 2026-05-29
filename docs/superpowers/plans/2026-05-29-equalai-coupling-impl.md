# EqualAI × NonequalAI Coupling — Contained Fix Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Repair the synth-dep overwrite so that an `:EqualAI` *dependent*
parameter whose Haldane/Wegscheider RHS references a `:NonequalAI` symbol
receives a **distinct** synthesized inactive-state name instead of self-mapping
and overwriting its active-state assignment — making PK and `m_mixed` compute
correctly.

**Architecture:** Four call sites in `src/rate_eq_derivation.jl` build a
`rename_T` map with an identical "Pass 2" loop:
`rename_T[k] = name(_flip_to_inactive(_param_for_symbol(am, k)), am)`. For an
`:EqualAI` dependent param, `_flip_to_inactive` correctly returns the input
unchanged, so `rename_T[k] == k` (self-map). Downstream, `dep_T[k]` then
overwrites `dep_R[k]` in `merge(dep_R, dep_T)` (and the polynomial T-rename is a
no-op), producing the wrong rate. Fix: a shared helper that, when the flip is a
no-op, constructs the explicit `:I`-state variant so the inactive name is
distinct. Apply the **same** helper at all four sites so the synthesized inactive
names are *identical* across the dep-expr, dep-assignment, and polynomial paths
(consistency is the correctness property).

**Tech Stack:** Julia, `@generated` rate-equation derivation, the
`name(p, m)` parameter-naming chokepoint, `Pkg.test()`.

**Scope guardrails (do NOT cross — see the spec):**
- Do **not** touch the parameter-naming chokepoint internals (`_state_tag`,
  `_render_binding`, `_render_iso`) or `_flip_to_inactive`'s semantics.
- Do **not** change the fifth flip site `src/mechanism_enumeration.jl:2000`
  (`_build_name_map`) — it feeds the **deferred** canonical-hash partition test
  (21 vs 23), the parent session's territory.
- Do **not** add a rank validator, change the enumeration, introduce `√`, or any
  representation change — those are the follow-up PRs (see
  `docs/superpowers/specs/2026-05-29-direction-symmetry-constraint-resolution.md`
  and `2026-05-29-nonequalai-rank-validity.md`).
- No test deletion; no `@testset` / `MECHANISM_TEST_SPECS` removal.

**Spec:** `docs/superpowers/specs/2026-05-29-equalai-nonequalai-coupling-design.md`

---

## File Structure

- `src/types.jl` — add `_force_inactive(p)` next to `_flip_to_inactive`
  (~line 1455): the inactive-state variant of a parameter *regardless* of its
  allosteric tag.
- `src/rate_eq_derivation.jl` — add `_dep_inactive_name(am, k)` and
  `_add_case_b_renames!(rename_T, deps, am)`; replace the four identical Pass-2
  loops (lines ~766, ~1277, ~1424, ~1543) with calls to the helper.
- `test/test_rate_eq_derivation.jl` — keep `m_mixed` and the PK goldens
  (they become the regression tests); add one Wegscheider regression test
  (random-order mechanism, rate ≈ 0 at chemical equilibrium).
- `.claude/CLAUDE.md` — one note under "Allosteric state taxonomy".

---

## Task 1: Establish the failing baseline

**Files:**
- Test: `test/test_rate_eq_derivation.jl` (existing)

- [ ] **Step 1: Run the failing allosteric tests and capture output**

Run:
```bash
julia --project -e 'using EnzymeRates, Test; include("test/mechanism_definitions_for_test_enzyme_derivation.jl"); include("test/test_rate_eq_derivation.jl")'
```
Expected: the `Allosteric edge cases` testset shows `m_mixed` rate ≈ 2.615 at
equilibrium (should be 0); the `PK` testset shows `n_haldane_constraints = 1`
(expected 2) and `n_mirror_constraints = 0` (expected 4); "Allosteric Analytical
Rate" / "PK Haldane Equilibrium" fail. Record the exact numbers — these are the
TDD targets.

- [ ] **Step 2: No commit** (baseline only).

---

## Task 2: Add `_force_inactive` (tag-agnostic inactive variant)

**Files:**
- Modify: `src/types.jl` (next to `_flip_to_inactive`, ~line 1455)
- Test: `test/test_types.jl`

- [ ] **Step 1: Write the failing test**

In `test/test_types.jl`, add (use a real `Step`/`Krev` per the patterns already
in that file; the assertion is the point):
```julia
@testset "_force_inactive forces :I regardless of tag" begin
    # An :EqualAI parameter has no :I variant under _flip_to_inactive
    # (returns itself); _force_inactive must return the explicit :I variant.
    s = EnzymeRates.Step(/* an SS step from an existing test fixture */)
    p_eq = EnzymeRates.Krev(s, :EqualAI)
    @test EnzymeRates._flip_to_inactive(p_eq) === p_eq          # unchanged
    @test EnzymeRates._force_inactive(p_eq) == EnzymeRates.Krev(s, :I)
    p_a = EnzymeRates.Krev(s, :A)
    @test EnzymeRates._force_inactive(p_a) == EnzymeRates.Krev(s, :I)
end
```
(Reuse a `Step` already constructed elsewhere in `test/test_types.jl`; do not
invent fields.)

- [ ] **Step 2: Run test to verify it fails**

Run: `julia --project -e 'using Pkg; Pkg.test()'` (or include just
`test/test_types.jl`).
Expected: FAIL — `_force_inactive` not defined.

- [ ] **Step 3: Implement `_force_inactive`**

In `src/types.jl`, immediately after the `_flip_to_inactive` definitions:
```julia
# Inactive-state variant of a parameter REGARDLESS of its allosteric tag.
# Unlike `_flip_to_inactive` (which returns an `:EqualAI`/`:None` param
# unchanged), this forces the `:I` state, used to give a *dependent* `:EqualAI`
# parameter a distinct inactive name when the Haldane/Wegscheider relation
# makes it differ between states.
_force_inactive(p::P) where {P <: Union{Kd, Kiso, Kon, Koff, Kfor, Krev}} =
    P(p.step, :I)
_force_inactive(p::Kreg) = Kreg(p.site, p.ligand, :I)
```

- [ ] **Step 4: Run test to verify it passes**

Run: `julia --project -e 'using Pkg; Pkg.test()'` (or just `test/test_types.jl`).
Expected: PASS.

- [ ] **Step 5: Commit**

```bash
git add src/types.jl test/test_types.jl
git commit -m "feat: _force_inactive yields the :I variant regardless of allo tag

Co-Authored-By: Claude Opus 4.8 (1M context) <noreply@anthropic.com>"
```

---

## Task 3: Add the shared rename helper and verify a distinct inactive name

**Files:**
- Modify: `src/rate_eq_derivation.jl` (add helpers near the synth-dep code, e.g.
  just above `_dependent_param_exprs(::Type{AllostericEnzymeMechanism…})`,
  ~line 1235)

- [ ] **Step 1: Write the failing test**

In `test/test_rate_eq_derivation.jl`, add a focused unit test using the existing
`m_mixed`'s catalytic mechanism (the `:EqualAI` catalysis reverse `k2r` is the
dep whose RHS references the `:NonequalAI` S-binding K):
```julia
@testset "_dep_inactive_name distinct for promoted EqualAI dep" begin
    am = EnzymeRates.AllostericMechanism(m_mixed)   # m_mixed defined above in this file
    dep_R, _ = EnzymeRates._dependent_param_exprs_allosteric(am)
    nonequalai = Set(EnzymeRates.name(p_R, am)
                     for (p_R, _) in EnzymeRates._I_rename_parameters(am))
    # find an EqualAI dep whose RHS references a NonequalAI symbol
    k = first(k for (k, v) in dep_R
              if EnzymeRates._expr_references_any(v, nonequalai)
                 && !(k in nonequalai))
    @test EnzymeRates._dep_inactive_name(am, k) != k
end
```

- [ ] **Step 2: Run test to verify it fails**

Run: `julia --project -e 'using Pkg; Pkg.test()'`
Expected: FAIL — `_dep_inactive_name` not defined.

- [ ] **Step 3: Implement the helpers**

In `src/rate_eq_derivation.jl`:
```julia
# Distinct inactive-state name for a *dependent* parameter being promoted to
# per-state (Case B: an `:EqualAI` dep whose Haldane/Wegscheider RHS references
# a `:NonequalAI` symbol). For a `:NonequalAI`/`:A` dep, `_flip_to_inactive`
# already yields a distinct `:I` name; for an `:EqualAI` dep it is a no-op, so
# fall back to the forced `:I` variant to avoid a self-map.
function _dep_inactive_name(am, k::Symbol)
    p = _param_for_symbol(am, k)
    nm = name(_flip_to_inactive(p), am)
    nm == k ? name(_force_inactive(p), am) : nm
end

# Pass 2 of the I-rename construction, shared across the four synth-dep sites so
# the synthesized inactive names are identical everywhere. Returns the Pass-1
# key set (callers that need to distinguish synthesized entries use it).
function _add_case_b_renames!(rename_T::Dict{Symbol, Symbol}, deps, am)
    renamed_set = Set{Symbol}(keys(rename_T))
    for (k, v) in deps
        haskey(rename_T, k) && continue
        _expr_references_any(v, renamed_set) || continue
        rename_T[k] = _dep_inactive_name(am, k)
    end
    return renamed_set
end
```

- [ ] **Step 4: Run test to verify it passes**

Run: `julia --project -e 'using Pkg; Pkg.test()'`
Expected: PASS for the new unit test. (Other allosteric tests may still fail —
the sites are not wired up yet; that is Task 4.)

- [ ] **Step 5: Verify the chokepoint renders the `:I` variant distinctly**

Sanity-check interactively that `name(_force_inactive(p), am)` produces an
inactive-tokened symbol distinct from `name(p, am)` for an `:EqualAI` `p`
(the structural chokepoint renders `:I` with a state token, `:EqualAI` without).
If it does NOT (e.g. the chokepoint consults the group's `cat_allo_state` and
ignores `p.state`), STOP and report — the helper needs a different inactive-name
source and the chokepoint is the parent's territory.

- [ ] **Step 6: Commit**

```bash
git add src/rate_eq_derivation.jl test/test_rate_eq_derivation.jl
git commit -m "feat: shared synth-dep I-rename helper with distinct EqualAI promotion name

Co-Authored-By: Claude Opus 4.8 (1M context) <noreply@anthropic.com>"
```

---

## Task 4: Wire all four synth-dep sites to the helper

**Files:**
- Modify: `src/rate_eq_derivation.jl` lines ~766, ~1277, ~1424, ~1543

Each site currently has this block (variable name is `dep_R_all` at 766/1277/1543
and `dep_R` at 1424):
```julia
    renamed_set = Set{Symbol}(keys(rename_T))
    for (k, v) in dep_R_all
        haskey(rename_T, k) && continue
        _expr_references_any(v, renamed_set) || continue
        rename_T[k] = name(_flip_to_inactive(_param_for_symbol(am, k)), am)
    end
```

- [ ] **Step 1: Replace the block at each of the four sites**

Replace each occurrence with (matching the local deps variable name):
```julia
    renamed_set = _add_case_b_renames!(rename_T, dep_R_all, am)
```
(Use `dep_R` at the `_build_dep_assignments` site, line ~1424, where the local is
named `dep_R`.) Confirm `renamed_set` is still in scope for the later uses in
`_build_dep_assignments` (lines ~1463) and `_dependent_param_exprs` — the helper
returns it, so the assignment preserves it.

- [ ] **Step 2: Run the targeted allosteric tests**

Run:
```bash
julia --project -e 'using EnzymeRates, Test; include("test/mechanism_definitions_for_test_enzyme_derivation.jl"); include("test/test_rate_eq_derivation.jl")'
```
Expected: `m_mixed` rate ≈ 0 at equilibrium; PK `n_haldane_constraints = 2`,
`n_mirror_constraints = 4`; "PK Haldane Equilibrium" and "Allosteric Analytical
Rate" PASS.

- [ ] **Step 3: If any still fail, STOP and diagnose** (do not pile on fixes).
Check that all four sites produce identical synthesized names (a mismatch
between the polynomial rename at 766/1543 and the dep map at 1277/1424 is the
likely culprit) — they should, since all four now call the same helper.

- [ ] **Step 4: Commit**

```bash
git add src/rate_eq_derivation.jl
git commit -m "fix: synth-dep emits distinct inactive name for EqualAI deps coupled to NonequalAI

Co-Authored-By: Claude Opus 4.8 (1M context) <noreply@anthropic.com>"
```

---

## Task 5: Wegscheider regression test

**Files:**
- Test: `test/test_rate_eq_derivation.jl` (in the `Allosteric edge cases`
  testset, after `m_mixed`)

The four targeted tests are Haldane-side. The **"Random-order Bi-Bi"** spec
mechanism (`test/mechanism_definitions_for_test_enzyme_derivation.jl` ~line 504,
`n_wegscheider=1`, `n_mirror=0`) has a genuine independent Wegscheider cycle (the
ungrouped binding square `E→EA→EAB→EB→E`). Made allosteric with one binding
group `:NonequalAI`, its Wegscheider-dependent `:EqualAI` binding K hits the same
Case-B path on a **Wegscheider** cycle (not a Haldane cycle), so the contained
fix must reach it. The contained fix *computes* this by promoting that dependent
binding K (it does NOT reject — rejection is the rank-algorithm follow-up; see
`docs/superpowers/specs/2026-05-29-nonequalai-rank-validity.md`). Assert the
mechanism-agnostic invariant — **zero net rate at chemical equilibrium** — which
needs no full analytical pre-derivation.

- [ ] **Step 1: Write the test**

```julia
# Random-order Bi-Bi: the ungrouped binding square is a genuine Wegscheider
# cycle. One binding is :NonequalAI, so the Wegscheider-dependent :EqualAI
# binding K is promoted (Case B on a WEGSCHEIDER cycle). Any thermodynamically
# consistent mechanism MUST give zero net rate at chemical equilibrium.
# (This config is over-parametrized — the rank-validity follow-up will REJECT
#  it; here we only verify the contained fix computes it consistently.)
cm_ro = @enzyme_mechanism begin
    substrates: A, B
    products:   P, Q
    steps: begin
        E + A <--> E(A)
        E + B <--> E(B)
        E(A) + B <--> E(A, B)
        E(B) + A <--> E(A, B)
        E(A, B) <--> E(P, Q)
        E(P, Q) <--> E(Q) + P
        E(Q) <--> E + Q
    end
end
m_ro = EnzymeRates.AllostericEnzymeMechanism(
    cm_ro,
    (2, (:NonequalAI, :EqualAI, :EqualAI, :EqualAI, :EqualAI, :EqualAI, :EqualAI)),
    (((:I,), 2, (:NonequalAI,)),),
)
# Overall A + B ⇌ P + Q, so Keq = P·Q/(A·B).
Keq_ro = 4.0
A_eq, B_eq = 1.5, 2.0
P_eq, Q_eq = 3.0, (Keq_ro * A_eq * B_eq / 3.0)
p_ro = (/* fill independent params from rate_equation_string(m_ro); set Keq=Keq_ro, E_total=1.0, L=2.0, reg K's positive */)
@test isapprox(
    rate_equation(m_ro, (A=A_eq, B=B_eq, P=P_eq, Q=Q_eq, I=0.5), p_ro),
    0.0; atol=1e-9)
```
Before writing the literal `p_ro`, obtain the exact independent-parameter names
by printing `EnzymeRates.rate_equation_string(m_ro)` and `fitted_params(m_ro)`;
set each to a positive value and `Keq` to `Keq_ro`. The expected value (0) is
pre-derived from thermodynamics, independent of the code under test.

- [ ] **Step 2: Confirm it exercises a Wegscheider Case-B promotion**

Verify the random-order loop actually yields a Wegscheider-derived `:EqualAI`
dep referencing the `:NonequalAI` binding: check that
`_dependent_param_exprs_allosteric(AllostericMechanism(m_ro))` has a dep whose
RHS references the NonequalAI symbol and whose own group is `:EqualAI`. If the
elimination happens to pivot differently and no such dep exists, adjust which
binding is `:NonequalAI` (or the step ordering) until it does, then re-derive
`p_ro`. Note this in the test comment.

- [ ] **Step 3: Run the test**

Run the include line from Task 1.
Expected: PASS (rate ≈ 0). (Optionally confirm it would have FAILED pre-fix by
stashing Task 4 — informative but not required.)

- [ ] **Step 4: Commit**

```bash
git add test/test_rate_eq_derivation.jl
git commit -m "test: Wegscheider-cycle EqualAI×NonequalAI gives zero rate at equilibrium

Co-Authored-By: Claude Opus 4.8 (1M context) <noreply@anthropic.com>"
```

---

## Task 6: Full suite green (except the deferred hash partition)

**Files:** none (verification)

- [ ] **Step 1: Run the full suite**

Run: `julia --project -e 'using Pkg; Pkg.test()'`
Expected: all pass **except** the `canonical-hash partition stability` test
(deferred). Note: fixing `_dependent_param_exprs` changes the dep-expr set the
canonical hash consumes, so that deferred test's failure numbers may *shift*
(e.g. no longer exactly 21 vs 23). That is expected — **do not chase it**.

- [ ] **Step 2: Confirm no NON-deferred test newly fails**

Diff the failing-test list against the known-deferred set
(`canonical-hash partition stability` only). If any other test fails, STOP and
diagnose the root cause before proceeding (do not adjust unrelated goldens).

- [ ] **Step 3: No commit** (verification only).

---

## Task 7: Document the constraint in CLAUDE.md

**Files:**
- Modify: `.claude/CLAUDE.md` — "Allosteric state taxonomy" section.

- [ ] **Step 1: Add the note**

Append to that section:
```markdown
- An `:EqualAI` group's *dependent* parameter (e.g. an SS reverse rate, or a
  Wegscheider-dependent binding K) may legitimately differ between A and I
  states when its Haldane/Wegscheider RHS references a `:NonequalAI` symbol.
  The synth-dep machinery synthesizes a **distinct** inactive-state name for
  such a dependent param (shared `_dep_inactive_name` / `_add_case_b_renames!`
  helper) — it is NOT a tag violation, since the dependent value is derived,
  not user-shared. The planned direction-symmetric resolution of this
  (speeds shared, ratios derived) and the NonequalAI degeneracy rejection are
  follow-up PRs (see docs/superpowers/specs/2026-05-29-direction-symmetry-
  constraint-resolution.md and 2026-05-29-nonequalai-rank-validity.md).
```

- [ ] **Step 2: Commit**

```bash
git add .claude/CLAUDE.md
git commit -m "docs: CLAUDE.md note on EqualAI dependent params differing between A/I

Co-Authored-By: Claude Opus 4.8 (1M context) <noreply@anthropic.com>"
```

---

## Self-Review checklist (done while writing)

- **Spec coverage:** synth-dep overwrite repair (Tasks 2–4); Haldane + Wegscheider
  (one helper covers both — all four loops iterate full dep set; Task 5 proves
  Wegscheider); m_mixed kept and asserted 0 (Task 4/5); PK goldens 2/4 (Task 4);
  no validator/enumeration/√ (scope guardrails); CLAUDE.md note (Task 7). ✓
- **Placeholders:** the only intentionally-deferred literals are the `Step`
  fixture in Task 2 and `p_ro` in Task 5, both with explicit instructions to
  source real values from existing fixtures / `rate_equation_string` — not
  free-floating TODOs. ✓
- **Name consistency:** `_force_inactive` (types.jl), `_dep_inactive_name` /
  `_add_case_b_renames!` (rate_eq_derivation.jl) used identically across Tasks
  3–4. ✓
