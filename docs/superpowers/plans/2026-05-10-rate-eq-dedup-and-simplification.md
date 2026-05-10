# Rate-equation dedup + derivation simplification: implementation plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Collapse three sources of rate-equation `eq_hash` duplication (factoring variants, Ka↔Kd inversion, split-with-tie ≡ pre-merged) and simplify the derivation pipeline by ~640 net body lines, while preserving all rate-equation correctness.

**Architecture:** Four sequential commits in one PR. Commit 0 removes the parameter-identifiability code (independent prep). Commit 1 makes polynomial construction use Kd directly (Source B). Commit 2 always emits expanded polynomials and deletes the entire algebraic-factoring family + types + tests (Source A). Commit 3 absorbs single-symbol Wegscheider RE ties into the kinetic-group rename map and adds section-labeled output + canonicalizer normalization (Source C).

**Tech Stack:** Julia 1.9+, EnzymeRates.jl package at `/home/denis.linux/.julia/dev/EnzymeRates`. Tests use Julia's `Test`, `CSV`, `DataFrames`. Project has `Aqua` and `JET` quality gates.

**Reference spec:** `docs/superpowers/specs/2026-05-09-rate-eq-dedup-and-simplification-design.md`

**Reference investigation:** `dedup_investigation/investigation_eq_hash_duplication.md`. The 4-mechanism LDH n=7 cluster used by the dedup tests is at rows 22 (eq_hash `831e36af`), 27 (`9c7141ac`), 31 (`89f33d51`), 36 (`b362dd75`) of `dedup_investigation/cv_results.csv`.

---

## Pre-flight

### TDD development loop

Cold `Pkg.test()` runs are slow because of precompilation. Use a long-lived Julia REPL with Revise.jl for fast iteration:

```bash
cd /home/denis.linux/.julia/dev/EnzymeRates
julia --project
```

Then in the REPL:

```julia
using Revise
using EnzymeRates
using Test
include("test/mechanism_definitions_for_test_enzyme_derivation.jl")
include("test/test_eq_hash_dedup.jl")  # after creating in Task 0.8
```

After source edits, just re-`include` the affected test file — Revise picks up source changes without REPL restart.

For final validation use `Pkg.test()`:

```bash
julia --project -e 'using Pkg; Pkg.test()'
```

### Branch / worktree

The current branch is `improve-deduplication-of-ident-equations`. All commits go on this branch.

### Test data dependency

The dedup tests in `test/test_eq_hash_dedup.jl` and the CSV-replay test depend on:
- `dedup_investigation/cv_results.csv`
- `dedup_investigation/params_estimate_{5,6,7,8}.csv`

These files exist in the repo (untracked under `dedup_investigation/`). Do NOT add them to git. Do NOT delete them. Tests guard for missing files with skip-if-not-present.

---

## Commit 0 — Remove parameter-identifiability code + write dedup test scaffold

**Goal of commit:** Pure deletion of `structural_identifiability_deficit` and supporting machinery + write the eq-hash dedup test file (initially red across all sources). Identifiability removal is independent of the dedup work but lands first because it unblocks Commit 2's helper deletions.

**Files:**
- Modify: `src/EnzymeRates.jl`
- Modify: `src/rate_eq_derivation.jl`
- Modify: `src/identify_rate_equation.jl`
- Modify: `test/test_rate_eq_derivation.jl`
- Modify: `test/test_identify_rate_equation.jl`
- Modify: `test/mechanism_definitions_for_test_enzyme_derivation.jl`
- Modify: `.claude/CLAUDE.md`
- Create: `test/test_eq_hash_dedup.jl`
- Modify: `test/runtests.jl`

### Task 0.1: Drop the public export of `structural_identifiability_deficit`

**Files:**
- Modify: `src/EnzymeRates.jl:16`

- [ ] **Step 1: Open `src/EnzymeRates.jl` and remove the export line for `structural_identifiability_deficit`.**

Find the line:

```julia
export structural_identifiability_deficit
```

Delete it. Save.

- [ ] **Step 2: Reload in REPL, confirm package still loads.**

```julia
include("src/EnzymeRates.jl")  # or just rely on Revise
```

Expected: no errors. The function still exists internally; just not exported.

### Task 0.2: Delete `structural_identifiability_deficit` for `EnzymeMechanism`

**Files:**
- Modify: `src/rate_eq_derivation.jl:780–800`

- [ ] **Step 1: Locate and delete the function and its docstring.**

Around lines 780–800 of `src/rate_eq_derivation.jl`:

```julia
"""
    structural_identifiability_deficit(m::EnzymeMechanism) → Int

Number of excess parameters beyond what is structurally identifiable from kinetic data.
"""
@generated function structural_identifiability_deficit(::M) where {M <: EnzymeMechanism}
    _, indep = _dependent_param_exprs(M)
    n_k = length(indep)
    num_fs, denom_terms = _raw_symbolic_rate_polys(M)
    num = _expand_factored_sigma(num_fs)
    den = _expand_to_poly(denom_terms)
    mm(mono) = sort([
        s => e for (s, e) in mono
        if !is_k_parameter(s) && s != :E_total && s != :Keq
    ])
    n_num = length(unique(mm(k) for (k, _) in num))
    n_denom = length(unique(mm(k) for (k, _) in den))
    n_k - (n_num - 1) - (n_denom - 1)
end
```

Delete the docstring + function. Also delete the `# ─── Structural Identifiability ───` section header comment immediately above if it bounds only this function.

- [ ] **Step 2: Reload, confirm no parse errors.**

In REPL: trigger Revise reload (any edit to a tracked file triggers it). Confirm `EnzymeRates` namespace is still callable.

### Task 0.3: Delete `structural_identifiability_deficit` for `AllostericEnzymeMechanism` and `_count_allosteric_rate_monomials`

**Files:**
- Modify: `src/rate_eq_derivation.jl:1672–1779` (file end)

- [ ] **Step 1: Locate and delete both functions.**

The block runs from line 1672 to the end of file (1779). The `structural_identifiability_deficit` overload starts at 1674 and ends at 1682. The `_count_allosteric_rate_monomials` function starts at 1690 and runs to 1778. Together with the section header comment at 1672, this is ~110 lines.

Delete from the section comment header at line 1672 through the closing `end` at line 1778. If the file then ends with stray blank lines, leave those — `git diff` will show only meaningful deletions.

- [ ] **Step 2: Reload and verify.**

### Task 0.4: Remove the identifiability test helper and its caller

**Files:**
- Modify: `test/test_rate_eq_derivation.jl:475–484, 751`

- [ ] **Step 1: Delete the `test_identifiability` function.**

Around line 475:

```julia
function test_identifiability(spec::MechanismTestSpec)
    m = spec.mechanism
    @testset "Identifiability" begin
        # over-counts identifiable degrees of freedom for factored
        # forms; we only check the boolean direction with a relaxed
        # is_identifiable check — the magnitude is not biophysically
        # meaningful at its current shape.
        @test (structural_identifiability_deficit(m) <= 0) ==
              spec.expected_is_identifiable
    end
end
```

Delete the function.

- [ ] **Step 2: Delete the call site at line 751.**

Around line 751 — the line `test_identifiability(spec)` — delete it.

- [ ] **Step 3: Reword file-header comment at line 2.**

Around line 2:

```julia
# Validates structure, constraints, identifiability,
```

Drop "identifiability," from the list.

### Task 0.5: Remove the `expected_is_identifiable` field and all 37 fixture occurrences

**Files:**
- Modify: `test/mechanism_definitions_for_test_enzyme_derivation.jl`

- [ ] **Step 1: Delete the field declaration in the `MechanismTestSpec` struct.**

Around line 32:

```julia
    expected_is_identifiable::Bool
```

Delete the line.

- [ ] **Step 2: Delete all 37 `expected_is_identifiable=...` fixture lines.**

Use a search-and-delete pass. Count occurrences first:

```bash
grep -c "expected_is_identifiable=" test/mechanism_definitions_for_test_enzyme_derivation.jl
```

Expected: 38 (37 fixture-call-site lines + 1 field declaration deleted in Step 1; if Step 1 already ran, expect 37).

Then delete each occurrence. They're typically on their own line within a `MechanismTestSpec(...)` constructor call:

```julia
            expected_is_identifiable=true,
```

Delete each entire line. After deletion, re-run the count — should be 0.

```bash
grep -c "expected_is_identifiable=" test/mechanism_definitions_for_test_enzyme_derivation.jl
```

Expected: 0.

- [ ] **Step 3: Verify the file still parses.**

```bash
julia --project -e 'include("test/mechanism_definitions_for_test_enzyme_derivation.jl")'
```

Expected: no parse errors. (This won't run any tests; just loads the fixture definitions.)

### Task 0.6: Reword stale identifiability and Phase G comments

**Files:**
- Modify: `src/identify_rate_equation.jl:147, 246`
- Modify: `src/rate_eq_derivation.jl:22`
- Modify: `test/test_identify_rate_equation.jl:54`

- [ ] **Step 1: Reword `identify_rate_equation.jl:246` defensive-lookup comment.**

Around line 245–251 — the existing comment reads:

```julia
    # Defensive lookup: a fitted key may not appear in the body
    # (e.g., a structurally-unidentifiable ghost param on a
    # zeroed `:NonequalRT` path), in which case `spec_name_map`
    # has no entry. Fall back to the spec key itself in cached_params
    # if both maps lack the canonical token; if even that misses,
    # use NaN as a sentinel that downstream loss/CV will surface.
```

Replace with:

```julia
    # Defensive lookup: a fitted key may not appear in the body
    # (e.g., a parameter on a zeroed `:NonequalRT` path), in which
    # case `spec_name_map` has no entry. Fall back to the spec key
    # itself in cached_params if both maps lack the canonical token;
    # if even that misses, use NaN as a sentinel that downstream
    # loss/CV will surface.
```

- [ ] **Step 2: Reword `identify_rate_equation.jl:147` Phase G reference.**

Around line 146–149 — find:

```julia
`parameters(m, Full)` is defined for both `EnzymeMechanism` and
`AllostericEnzymeMechanism` (Phase G.0). Allosteric coverage
includes T-state names, regulator-site names, and the allosteric
coupling `L` automatically.
```

Replace with:

```julia
`parameters(m, Full)` is defined for both `EnzymeMechanism` and
`AllostericEnzymeMechanism`. Allosteric coverage includes T-state
names, regulator-site names, and the allosteric coupling `L`
automatically.
```

- [ ] **Step 3: Reword `rate_eq_derivation.jl:22` Phase G reference.**

Around lines 20–25 — find:

```julia
  allosteric Full mode is used as a name source by Phase G's
  rate-equation canonicalizer; no `rate_equation` method is
  defined for `(::AllostericEnzymeMechanism, ::FullMode)`, so
  this mode is for canonicalization, not runtime evaluation.
```

Replace with:

```julia
  allosteric Full mode is used as a name source by the
  rate-equation canonicalizer in `identify_rate_equation`; no
  `rate_equation` method is defined for
  `(::AllostericEnzymeMechanism, ::FullMode)`, so this mode is for
  canonicalization, not runtime evaluation.
```

- [ ] **Step 4: Reword `test/test_identify_rate_equation.jl:54` ghost-param comment.**

Around line 54 — find a comment reading approximately:

```julia
# 5 identifiable params + 1 ghost (k3f_T)
```

Replace with:

```julia
# 5 fitted params + 1 zeroed-path param (k3f_T)
```

### Task 0.7: Update CLAUDE.md to drop identifiability references

**Files:**
- Modify: `.claude/CLAUDE.md:291`

- [ ] **Step 1: Find the rate_eq_derivation.jl bullet.**

Around line 291:

```markdown
- `src/rate_eq_derivation.jl` — King-Altman/Cha rate equation derivation via `@generated` functions; parameters API; identifiability checks; kcat computation (`_is_ss_rate_constant`, `_kcat_components`, `_kcat_forward`); `rescale_parameter_values`; AllostericEnzymeMechanism MWC rate equation assembly (`_build_allosteric_rate_body`, `_count_allosteric_rate_monomials`, `rate_equation_string`, `structural_identifiability_deficit`); helpers `_onlyT_syms`, `_onlyR_syms`, `_T_rename`, `_build_kinetic_rename_map`, `_build_dep_assignments`.
```

Replace with:

```markdown
- `src/rate_eq_derivation.jl` — King-Altman/Cha rate equation derivation via `@generated` functions; parameters API; kcat computation (`_is_ss_rate_constant`, `_kcat_components`, `_kcat_forward`); `rescale_parameter_values`; AllostericEnzymeMechanism MWC rate equation assembly (`_build_allosteric_rate_body`, `rate_equation_string`); helpers `_onlyR_syms`, `_T_rename`, `_build_kinetic_rename_map`, `_build_dep_assignments`.
```

(Drops "identifiability checks", `_count_allosteric_rate_monomials`, `structural_identifiability_deficit`, and the non-existent `_onlyT_syms`.)

### Task 0.8: Create `test/test_eq_hash_dedup.jl` with TDD scaffold

**Files:**
- Create: `test/test_eq_hash_dedup.jl`
- Modify: `test/runtests.jl`

- [ ] **Step 1: Create the test file with all four testsets (initially RED for sources not yet fixed).**

Create `test/test_eq_hash_dedup.jl`:

```julia
# ABOUTME: Tests that algebraically-equivalent rate equations from different
# ABOUTME: mechanism specs collapse to one eq_hash (Sources A, B, C dedup).

using CSV, DataFrames, Test
using EnzymeRates
using EnzymeRates: _canonical_rate_eq_hash

const _CV_PATH = joinpath(@__DIR__, "..", "dedup_investigation", "cv_results.csv")

if !isfile(_CV_PATH)
    @warn "Skipping eq_hash dedup tests — $(_CV_PATH) not present."
else
    const _CV = CSV.read(_CV_PATH, DataFrame)

    # Reconstruct an EnzymeMechanism instance from a CSV row's mechanism_type column.
    _mech(row_idx) = eval(Meta.parse(_CV[row_idx, :mechanism_type]))()

    @testset "eq_hash dedup" begin
        # Sanity: confirm the row indices match the investigation's eq_hashes.
        @test _CV[22, :eq_hash] == "831e36af"
        @test _CV[27, :eq_hash] == "9c7141ac"
        @test _CV[31, :eq_hash] == "89f33d51"
        @test _CV[36, :eq_hash] == "b362dd75"

        @testset "Source A: factoring variants of same polynomial" begin
            # Mech 2 vs Mech 1 (rows 27 vs 22): denominator factored differently.
            @test _canonical_rate_eq_hash(_mech(27)) ==
                  _canonical_rate_eq_hash(_mech(22))
        end

        @testset "Source B: Ka↔Kd inversion artifact" begin
            # Mech 1 vs Mech 3 (rows 22 vs 31): both define K8, only Mech 1 uses it.
            @test _canonical_rate_eq_hash(_mech(22)) ==
                  _canonical_rate_eq_hash(_mech(31))
            # Mech 2 vs Mech 4 (rows 27 vs 36): Mech 4 has dead K12 = 1/(1/K2).
            @test _canonical_rate_eq_hash(_mech(27)) ==
                  _canonical_rate_eq_hash(_mech(36))
        end

        @testset "Source C: split-with-tie ≡ pre-merged" begin
            # All four cluster into one canonical equation.
            for j in (27, 31, 36)
                @test _canonical_rate_eq_hash(_mech(22)) ==
                      _canonical_rate_eq_hash(_mech(j))
            end
        end

        @testset "Section labels render correctly" begin
            α = _mech(22)
            s = rate_equation_string(α)
            @test occursin("# User defined constraints:", s)
            @test occursin("(substituted into v)", s)

            β = _mech(31)
            s = rate_equation_string(β)
            @test occursin("# Wegscheider constraints:", s)
            @test occursin("(substituted into v)", s)
        end
    end
end
```

- [ ] **Step 2: Add include line to `test/runtests.jl`.**

Insert immediately after the existing `include("test_identify_rate_equation.jl")` line, so dedup tests run after the regular rate-eq derivation tests:

```julia
    include("test_eq_hash_dedup.jl")
```

The full block should now read:

```julia
@testset "EnzymeRates.jl" begin
    include("test_accessors.jl")
    include("test_types.jl")
    include("test_dsl.jl")
    include("test_sym_poly.jl")
    include("test_rate_eq_derivation.jl")
    include("test_fitting.jl")
    include("test_mechanism_enumeration.jl")
    include("test_identify_rate_equation.jl")
    include("test_eq_hash_dedup.jl")
    include("test_readme_runs.jl")
    include("test_aqua_jet.jl")
end
```

- [ ] **Step 3: Run the new test file in isolation; confirm Sources A, B, C are red and the sanity tests pass.**

```julia
include("test/test_eq_hash_dedup.jl")
```

Expected:
- Sanity tests (4 row-index assertions) pass.
- Source A test: FAIL (4 hashes are distinct under current code).
- Source B tests: FAIL.
- Source C tests: FAIL.
- Section-label tests: FAIL (no section headers in current output).

This is correct — these tests will turn green incrementally over Commits 1–3.

### Task 0.9: Run the full test suite; confirm all non-dedup tests still pass

- [ ] **Step 1: Run full test suite.**

```bash
julia --project -e 'using Pkg; Pkg.test()'
```

Expected: all existing tests pass. The only new failures are in `test_eq_hash_dedup.jl` (Sources A, B, C — expected red until later commits; sanity tests should be green).

If any *other* test fails, investigate before continuing — Commit 0 is supposed to be a pure-deletion commit and should not break anything outside the deleted helpers.

### Task 0.10: Commit 0

- [ ] **Step 1: Stage and commit.**

```bash
git add src/EnzymeRates.jl src/rate_eq_derivation.jl src/identify_rate_equation.jl \
        test/test_rate_eq_derivation.jl test/test_identify_rate_equation.jl \
        test/mechanism_definitions_for_test_enzyme_derivation.jl \
        test/test_eq_hash_dedup.jl test/runtests.jl \
        .claude/CLAUDE.md
git commit -m "$(cat <<'EOF'
Remove parameter-identifiability code; add dedup test scaffold

structural_identifiability_deficit's current implementation is not
correct for the use cases that matter and will be redesigned later.
Removing it now unblocks the polynomial-helper deletions in the
following commits and avoids constraining the redesign.

Adds test/test_eq_hash_dedup.jl with TDD tests for Sources A, B, C of
eq_hash duplication. The Source-A/B/C testsets are intentionally red at
this commit and turn green over the following three commits.

Co-Authored-By: Claude Opus 4.7 (1M context) <noreply@anthropic.com>
EOF
)"
```

---

## Commit 1 — Kd-by-construction (Source B fix)

**Goal of commit:** Polynomials store binding K's directly with the dissociation-constant convention. Eliminate the Ka↔Kd inversion layer and all `1/(1/X)` artifacts. Source-B dedup tests turn green.

**Files:**
- Modify: `src/rate_eq_derivation.jl` (`_compute_alpha`, `_kcat_forward` for both type families, `_raw_rate_expr_and_symbols`, `_rate_v_line`)
- Modify: `src/sym_poly_for_rate_eq_derivation.jl` (`_poly_to_expr`, `_factored_sigma_to_expr`, `_factored_poly_to_expr`, `_denom_terms_to_expr`, `to_rate_expr` — drop `inverted_params` parameter)
- Modify: `src/thermodynamic_constr_for_rate_eq_derivation.jl` (`_dependent_param_exprs`, drop `_apply_kd_inversion`, `_constraint_expr_strings`, `_build_rate_body`)
- Modify: existing rate-equation-string fixtures that show `1/(1/X)` artifacts

### Task 1.1: Audit fixture files for `1/(1/X)` patterns and predict the post-Step-1 expected strings

**Files:**
- Read: `test/mechanism_definitions_for_test_enzyme_derivation.jl`
- Read: `test/test_rate_eq_derivation.jl`

- [ ] **Step 1: Scan fixtures for the inversion artifact.**

```bash
grep -n "1 / (1 /" test/mechanism_definitions_for_test_enzyme_derivation.jl test/test_rate_eq_derivation.jl
```

For each match, the post-Step-1 expected text drops the `1 / (1 /` and trailing `)`. Example transformation:

- Before: `K8 = 1 / (1 / K4)`
- After: `K8 = K4`

- Before: `(1 / Keq) * (1 / (1 / K1)) * (1 / K2) * k3f`
- After: `(1 / Keq) * K1 * (1 / K2) * k3f`

- Before: `(1 + S / K1)` where the polynomial monomials currently print as `S * (1/K1)` due to inversion
- After: same `(1 + S / K1)` — the printed form was already Kd-style; only Haldane RHS expressions and the factored bases change

- [ ] **Step 2: Make a notes file or list of the exact fixtures and substitutions to perform in Task 1.2.**

This is reading-only; no edits yet. The substitutions are mechanical text-edits applied in the next task.

### Task 1.2: Update fixtures with post-Step-1 expected strings (TDD: failing tests come first)

**Files:**
- Modify: `test/mechanism_definitions_for_test_enzyme_derivation.jl`
- Modify: `test/test_rate_eq_derivation.jl:1142–1166` (byte-identical allosteric fixture)

- [ ] **Step 1: For each fixture identified in Task 1.1, rewrite the expected string.**

Apply the mechanical substitution:
- Replace `1 / (1 / K_i)` → `K_i` everywhere, for each `K_i` symbol.
- Replace `(1 / (1 / K_i))` → `K_i` (with surrounding parens dropped if reasonable).

For the byte-identical allosteric fixture at `test/test_rate_eq_derivation.jl:1160–1164`:

Before:

```julia
    expected = raw"""(; K1, K2, k3f, K1_T, K2_T, k3f_T, K_R_reg1, K_R_T_reg1, L, Keq, E_total) = params
(; S, P, R) = concs
k3r = (1 / Keq) * (1 / (1 / K1)) * (1 / K2) * k3f
k3r_T = (1 / Keq) * (1 / (1 / K1_T)) * (1 / K2_T) * k3f_T
v = E_total * (...) / (...)"""
```

After:

```julia
    expected = raw"""(; K1, K2, k3f, K1_T, K2_T, k3f_T, K_R_reg1, K_R_T_reg1, L, Keq, E_total) = params
(; S, P, R) = concs
k3r = (1 / Keq) * K1 * (1 / K2) * k3f
k3r_T = (1 / Keq) * K1_T * (1 / K2_T) * k3f_T
v = E_total * (...) / (...)"""
```

(Keep the v-line as-is — Step 2 will rewrite that part.)

- [ ] **Step 2: Run the affected fixture tests to confirm they FAIL under current code.**

In REPL:

```julia
include("test/test_rate_eq_derivation.jl")
```

Expected: the byte-identical allosteric fixture and any other fixture you updated should now FAIL (current code produces the old `1/(1/K)` form). The Source-B dedup tests in `test_eq_hash_dedup.jl` are also still failing.

If a fixture updated in Step 1 passes immediately, you wrote the wrong post-Step-1 expectation; revisit.

### Task 1.3: Flip `_compute_alpha` to put binding K's in `alpha_den`

**Files:**
- Modify: `src/rate_eq_derivation.jl:149–179`

- [ ] **Step 1: Locate `_compute_alpha` and find the binding-step accumulator.**

Around lines 168–177:

```julia
                K = poly_sym(Symbol("K$idx"))
                if i_f == cur && j_f ∉ visited
                    alpha_num[j_f] = poly_mul(poly_mul(alpha_num[cur], K), mp(m_l))
                    alpha_den[j_f] = poly_mul(alpha_den[cur], mp(m_r))
                    push!(visited, j_f); push!(queue, j_f)
                elseif j_f == cur && i_f ∉ visited
                    alpha_num[i_f] = poly_mul(alpha_num[cur], mp(m_r))
                    alpha_den[i_f] = poly_mul(poly_mul(alpha_den[cur], K), mp(m_l))
                    push!(visited, i_f); push!(queue, i_f)
                end
```

- [ ] **Step 2: Move K from `alpha_num` to `alpha_den` for binding-direction steps; keep iso steps unchanged.**

Reasoning: alpha is the relative concentration `[E_form]/[E_ref]`. For a binding step `E + S ⇌ ES`, Kd convention gives `[ES]/[E] = [S]/K_d`. Forward traversal (cur=E, target=ES) puts `S` in `alpha_num` and `K_d` in `alpha_den`; backward traversal (cur=ES, target=E) puts `K_d` in `alpha_num` and `S` in `alpha_den`. For iso steps (no metabolite), K stays in `alpha_num` for the forward-traversed form (Ka convention unchanged from original code).

Replace the binding/iso code block at lines ~168–177 of `_compute_alpha` with the version below:

```julia
                K = poly_sym(Symbol("K$idx"))
                # iso = both sides have no metabolite (pure isomerization)
                is_iso = isempty(m_l) && isempty(m_r)
                if i_f == cur && j_f ∉ visited
                    if is_iso
                        # iso: K (Ka) in alpha_num for the more-isomerized form
                        alpha_num[j_f] = poly_mul(alpha_num[cur], K)
                        alpha_den[j_f] = alpha_den[cur]
                    else
                        # binding step: Kd in alpha_den (Kd convention)
                        alpha_num[j_f] = poly_mul(alpha_num[cur], mp(m_l))
                        alpha_den[j_f] = poly_mul(poly_mul(alpha_den[cur], K), mp(m_r))
                    end
                    push!(visited, j_f); push!(queue, j_f)
                elseif j_f == cur && i_f ∉ visited
                    if is_iso
                        alpha_num[i_f] = alpha_num[cur]
                        alpha_den[i_f] = poly_mul(alpha_den[cur], K)
                    else
                        # backward binding: Kd into alpha_num, m_l (the bound met) into alpha_den
                        alpha_num[i_f] = poly_mul(alpha_num[cur], K)
                        alpha_den[i_f] = poly_mul(alpha_den[cur], mp(m_l))
                    end
                    push!(visited, i_f); push!(queue, i_f)
                end
```

- [ ] **Step 3: Reload, run a single-mechanism numerical-equivalence test as a smoke check.**

In REPL — pick a simple mechanism from the fixtures (e.g., the first uni-uni RE one) and check that `rate_equation` still gives the right number for randomly-sampled params:

```julia
include("test/mechanism_definitions_for_test_enzyme_derivation.jl")
m = first_simple_re_spec.mechanism  # use whichever spec name exists for a uni-uni RE mech
# Sample params, compute rate_equation, compare to ode_steady_state_flux:
# (full smoke check is run by Pkg.test() at end of commit)
```

If `_compute_alpha` is the only change so far and Haldane elimination still does Ka, the rate equation will be NUMERICALLY WRONG by a factor of K^n until Task 1.4 is complete. That's expected — defer numerical validation to after Task 1.4.

### Task 1.4: Sign-flip A-matrix in `_dependent_param_exprs` for binding K columns

**Files:**
- Modify: `src/thermodynamic_constr_for_rate_eq_derivation.jl:198–212`

- [ ] **Step 1: Locate the A-matrix construction loop.**

Around lines 198–212:

```julia
    # Translate cycle-incidence columns into the merged-parameter A matrix.
    # Non-representative steps' columns are folded into their representative
    # via the rename map; this is mathematically equivalent to a kinetic-group
    # equality constraint (K_idx = K_rep, k_idx_f = k_rep_f, ...).
    A = zeros(Rational{BigInt}, nc, n_vars)
    rhs = Rational{BigInt}.(rhs_coeffs)
    for i in 1:nc, j in 1:nsteps
        C[i, j] == 0 && continue
        if eq_steps[j]
            sym = Symbol("K$j")
            sym = get(rename, sym, sym)
            A[i, sym_col[sym]] += C[i, j]
        else
            kf = Symbol("k$(j)f"); kr = Symbol("k$(j)r")
            kf = get(rename, kf, kf); kr = get(rename, kr, kr)
            A[i, sym_col[kf]] += C[i, j]
            A[i, sym_col[kr]] -= C[i, j]
        end
    end
```

- [ ] **Step 2: Add a binding-K set and apply the sign-flip.**

Replace the loop with:

```julia
    # Binding K's are Kd in the polynomial; cycle products use 1/Kd, so
    # binding-K column entries get a sign flip on top of the cycle incidence.
    enz_set_for_binding = Set(enz_names)
    binding_K_set = Set{Symbol}()
    for (j, (lhs, _, _, _)) in enumerate(rxns)
        eq_steps[j] || continue
        any(s ∉ enz_set_for_binding for s in lhs) || continue
        sym = Symbol("K$j")
        push!(binding_K_set, get(rename, sym, sym))
    end

    A = zeros(Rational{BigInt}, nc, n_vars)
    rhs = Rational{BigInt}.(rhs_coeffs)
    for i in 1:nc, j in 1:nsteps
        C[i, j] == 0 && continue
        if eq_steps[j]
            sym = Symbol("K$j")
            sym = get(rename, sym, sym)
            sign_factor = sym in binding_K_set ? -1 : 1
            A[i, sym_col[sym]] += sign_factor * C[i, j]
        else
            kf = Symbol("k$(j)f"); kr = Symbol("k$(j)r")
            kf = get(rename, kf, kf); kr = get(rename, kr, kr)
            A[i, sym_col[kf]] += C[i, j]
            A[i, sym_col[kr]] -= C[i, j]
        end
    end
```

### Task 1.5: Drop `_apply_kd_inversion` and its three callers

**Files:**
- Modify: `src/thermodynamic_constr_for_rate_eq_derivation.jl:292–318, 348`
- Modify: `src/rate_eq_derivation.jl:895, 1467`

- [ ] **Step 1: Delete `_apply_kd_inversion`.**

In `src/thermodynamic_constr_for_rate_eq_derivation.jl`, around lines 292–306, delete the function and its docstring:

```julia
"""Apply K→1/K inversion to Haldane dep_exprs.
When a dependent param is itself a binding K, its RHS is wrapped in `inv_fn`
to compensate for the implicit LHS inversion (Ka→Kd)."""
function _apply_kd_inversion(dep_exprs, M::Type{<:EnzymeMechanism}, inv_fn)
    binding_Ks = Set(_binding_K_symbols(M))
    isempty(binding_Ks) && return dep_exprs
    inv_subs = Dict(K => inv_fn(K) for K in binding_Ks)
    Dict(
        k => begin
            rhs = substitute_params_expr(v, inv_subs)
            k in binding_Ks ? inv_fn(rhs) : rhs
        end
        for (k, v) in dep_exprs
    )
end
```

- [ ] **Step 2: Drop the call in `_constraint_expr_strings`.**

Around line 312:

```julia
function _constraint_expr_strings(M::Type{<:EnzymeMechanism})
    lines = String[]
    dep_exprs, _ = _dependent_param_exprs(M)
    if !isempty(dep_exprs)
        dep_exprs = _apply_kd_inversion(dep_exprs, M, K -> :(1 / $K))
        for (sym, expr) in sort(collect(dep_exprs); by=p -> string(p[1]))
            push!(lines, "$sym = $(string(expr))")
        end
    end
    lines
end
```

Replace with:

```julia
function _constraint_expr_strings(M::Type{<:EnzymeMechanism})
    lines = String[]
    dep_exprs, _ = _dependent_param_exprs(M)
    for (sym, expr) in sort(collect(dep_exprs); by=p -> string(p[1]))
        push!(lines, "$sym = $(string(expr))")
    end
    lines
end
```

- [ ] **Step 3: Drop the call in `_build_rate_body(M, ::Type{ReducedMode})`.**

Around line 348:

```julia
function _build_rate_body(M, ::Type{ReducedMode})
    expr, _, conc_syms = _raw_rate_expr_and_symbols(M)
    dep_exprs, indep = _dependent_param_exprs(M)
    dep_exprs = _apply_kd_inversion(dep_exprs, M, K -> :(inv($K)))
    hw_params = (indep..., :Keq, :E_total)
    ...
```

Remove the `dep_exprs = _apply_kd_inversion(...)` line. The remaining body unchanged.

- [ ] **Step 4: Drop the call in EnzymeMechanism `_kcat_forward`.**

In `src/rate_eq_derivation.jl` around line 895:

```julia
@generated function _kcat_forward(
    ::M, params::NamedTuple,
) where {M <: EnzymeMechanism}
    components = _kcat_components(M)
    dep_exprs, indep = _dependent_param_exprs(M)
    dep_exprs = _apply_kd_inversion(dep_exprs, M, K -> :(inv($K)))
    hw_params = (indep..., :Keq)
    ...
```

Remove the `_apply_kd_inversion` line.

Also drop the binding-K substitution block immediately below (lines ~899–907):

```julia
    # Apply Kd inversion to component expressions: raw polys use Ka,
    # but params store Kd for binding K's
    binding_Ks = Set(_binding_K_symbols(M))
    kd_subs = Dict(K => :(inv($K)) for K in binding_Ks)
    candidates = [
        :($(substitute_params_expr(nk, kd_subs)) /
          $(substitute_params_expr(dk, kd_subs)))
        for (nk, dk) in components
    ]
```

Replace with:

```julia
    candidates = [:($nk / $dk) for (nk, dk) in components]
```

- [ ] **Step 5: Drop the call in `_build_dep_assignments`.**

Around line 1467:

```julia
function _build_dep_assignments(
    M_type::Type{<:AllostericEnzymeMechanism}, inv_fn,
)
    m = M_type()
    CM = typeof(catalytic_mechanism(m))

    dep_R, indep_R = _dependent_param_exprs(CM)
    dep_R_kd = _apply_kd_inversion(dep_R, CM, inv_fn)
    sorted_deps = sort(collect(dep_R_kd); by=first)
    ...
```

Replace with:

```julia
function _build_dep_assignments(
    M_type::Type{<:AllostericEnzymeMechanism},
)
    m = M_type()
    CM = typeof(catalytic_mechanism(m))

    dep_R, indep_R = _dependent_param_exprs(CM)
    sorted_deps = sort(collect(dep_R); by=first)
    ...
```

(Drop the `inv_fn` parameter and the `_apply_kd_inversion` call. `dep_R_kd` rename → `dep_R` everywhere it's used downstream.)

- [ ] **Step 6: Update the two callers of `_build_dep_assignments` to drop the `inv_fn` argument.**

```bash
grep -n "_build_dep_assignments" /home/denis.linux/.julia/dev/EnzymeRates/src/rate_eq_derivation.jl
```

Two call sites: in `_build_allosteric_rate_body` (around line 1618) and in allosteric `rate_equation_string` (around line 1655). Each currently looks like:

```julia
    r_assignments, t_assignments_ = _build_dep_assignments(M_type, K -> :(inv($K)))
```

or

```julia
    r_assignments, t_assignments_ = _build_dep_assignments(M, K -> :(1 / $K))
```

Drop the second argument:

```julia
    r_assignments, t_assignments_ = _build_dep_assignments(M_type)
```

### Task 1.6: Drop `inverted_params` parameter from `_poly_to_expr` and the factored-expr helpers

**Files:**
- Modify: `src/sym_poly_for_rate_eq_derivation.jl:150–191, 375–443`

- [ ] **Step 1: Update `_poly_to_expr` signature and drop the `inverted_params` branch.**

Around line 150:

```julia
function _poly_to_expr(p::POLY, param_syms::Set{Symbol}, conc_syms::Set{Symbol},
                       inverted_params::Set{Symbol}=Set{Symbol}())
    isempty(p) && return 0
    pos, neg = Any[], Any[]
    sorted = sort(...)
    for (mono, coeff) in sorted
        ...
        for (s, e) in sorted_mono
            if s in inverted_params
                tgt, ex = e > 0 ? (df, e) : (nf, -e)
            else
                tgt, ex = e > 0 ? (nf, e) : (df, -e)
            end
            ex != 0 && push!(tgt, ex == 1 ? s : :($s ^ $ex))
        end
        ...
```

Replace with the simpler form (drop the `inverted_params` parameter and always use the second branch):

```julia
function _poly_to_expr(p::POLY, param_syms::Set{Symbol}, conc_syms::Set{Symbol})
    isempty(p) && return 0
    pos, neg = Any[], Any[]
    sorted = sort(
        collect(p);
        by=x -> (
            sum(e for (s,e) in x[1] if s ∉ param_syms; init=0),
            x[2] < 0,
            Tuple(string(s) for (s, _) in x[1]),
        ),
    )
    for (mono, coeff) in sorted
        nf, df = Any[], Any[]
        abs_c = abs(coeff)
        cn = Int(numerator(abs_c))
        cd = Int(denominator(abs_c))
        cn != 1 && push!(nf, cn)
        cd != 1 && push!(df, cd)
        sorted_mono = sort(
            mono;
            by=sp -> (sp.first in param_syms ? 0 : 1, string(sp.first)),
        )
        for (s, e) in sorted_mono
            tgt, ex = e > 0 ? (nf, e) : (df, -e)
            ex != 0 && push!(tgt, ex == 1 ? s : :($s ^ $ex))
        end
        num_part = isempty(nf) ? 1 : _nest_binary(:*, nf)
        term = isempty(df) ? num_part : :($num_part / $(_nest_binary(:*, df)))
        coeff > 0 ? push!(pos, term) : push!(neg, term)
    end
    pe = isempty(pos) ? nothing : _nest_binary(:+, pos)
    ne = isempty(neg) ? nothing : _nest_binary(:+, neg)
    pe !== nothing && ne !== nothing && return :($pe - $ne)
    pe !== nothing && return pe
    ne !== nothing && return :(- $ne)
    return 0
end
```

- [ ] **Step 2: Drop the `inverted_params` parameter from the four factored-expr helpers.**

In `src/sym_poly_for_rate_eq_derivation.jl` around lines 375, 400, 424, 451 — for each of `_factored_poly_to_expr`, `_factored_sigma_to_expr`, `_denom_terms_to_expr`, `to_rate_expr`: drop the `inverted_params::Set{Symbol}` parameter and remove `, inverted_params` from internal `_poly_to_expr` calls. (These functions are entirely deleted in Commit 2; this Step is a temporary tidying so the package compiles in the interim. If you'd rather skip the interim cleanup and accept a compile error until Commit 2 finishes, you can leave the parameters and just stop passing them — Julia will dispatch on default args.)

Easier alternative: change the default argument value from `Set{Symbol}()` to nothing-meaningful and stop passing from callers. They go away in Commit 2 anyway. Pick whichever keeps the package compiling.

- [ ] **Step 3: Drop `_poly_to_expr` overload at `rate_eq_derivation.jl:1516`.**

Around line 1516 there's a small overload:

```julia
function _poly_to_expr(p::POLY, param_syms, conc_syms, inv_set)
    fs = FactoredSigma([poly_one()], [FactoredPoly([p], [1])])
    _factored_sigma_to_expr(fs, param_syms, conc_syms, inv_set)
end
```

This is the four-arg version that wraps a POLY as FactoredSigma. After Step 2 it's gone anyway, but for now: remove the `inv_set` parameter or delete the overload entirely (since the three-arg `_poly_to_expr` in `sym_poly_for_rate_eq_derivation.jl` covers POLY directly).

### Task 1.7: Drop `_binding_K_symbols` non-allosteric callers and `inv_set` plumbing

**Files:**
- Modify: `src/rate_eq_derivation.jl:702–711, 749–763`

- [ ] **Step 1: Update `_raw_rate_expr_and_symbols` to drop `inv_set`.**

Around lines 702–711:

```julia
function _raw_rate_expr_and_symbols(M::Type{<:EnzymeMechanism})
    num, denom_terms = _raw_symbolic_rate_polys(M)
    m = M()
    param_syms = Set{Symbol}(_raw_param_symbols(m))
    conc_syms = Set{Symbol}(metabolites(m))
    inv_set = Set(_binding_K_symbols(M))
    expr = to_rate_expr(num, denom_terms, param_syms, conc_syms, inv_set)
    all_params = _sorted_raw_param_symbols(M)
    return expr, all_params, metabolites(m)
end
```

Replace with:

```julia
function _raw_rate_expr_and_symbols(M::Type{<:EnzymeMechanism})
    num, denom_terms = _raw_symbolic_rate_polys(M)
    m = M()
    param_syms = Set{Symbol}(_raw_param_symbols(m))
    conc_syms = Set{Symbol}(metabolites(m))
    expr = to_rate_expr(num, denom_terms, param_syms, conc_syms)
    all_params = _sorted_raw_param_symbols(M)
    return expr, all_params, metabolites(m)
end
```

- [ ] **Step 2: Update `_rate_v_line` similarly.**

Around lines 749–763:

```julia
function _rate_v_line(M::Type{<:EnzymeMechanism})
    num_fs, denom_terms = _raw_symbolic_rate_polys(M)
    m = M()
    ps = Set{Symbol}(_raw_param_symbols(m))
    cs = Set{Symbol}(metabolites(m))
    inv = Set(_binding_K_symbols(M))
    num_str = _expr_to_string(
        _factored_sigma_to_expr(num_fs, ps, cs, inv),
    )
    den_str = _expr_to_string(
        _denom_terms_to_expr(denom_terms, ps, cs, inv),
    )
    "v = E_total * ($num_str) / ($den_str)"
end
```

Drop the `inv = Set(_binding_K_symbols(M))` line and the trailing `, inv` arguments to the factored-expr converters. Keep the rest of the function intact for now — Commit 2 simplifies it further. Resulting form:

```julia
function _rate_v_line(M::Type{<:EnzymeMechanism})
    num_fs, denom_terms = _raw_symbolic_rate_polys(M)
    m = M()
    ps = Set{Symbol}(_raw_param_symbols(m))
    cs = Set{Symbol}(metabolites(m))
    num_str = _expr_to_string(_factored_sigma_to_expr(num_fs, ps, cs))
    den_str = _expr_to_string(_denom_terms_to_expr(denom_terms, ps, cs))
    "v = E_total * ($num_str) / ($den_str)"
end
```

- [ ] **Step 3: Audit remaining callers of `_binding_K_symbols`.**

```bash
grep -n "_binding_K_symbols" src/rate_eq_derivation.jl
```

After Tasks 1.5–1.7, callers should be only the allosteric branch (`_allosteric_num_den_exprs` lines around 1535, 1541; `_kcat_forward` allosteric lines around 973–974). Keep those — Commit 2 reworks the allosteric branch.

Also in factor-poly code around line 491 — that's deleted in Commit 2.

If any non-allosteric, non-factor-poly caller of `_binding_K_symbols` remains, you missed a substitution above; fix before continuing.

### Task 1.8: Drop binding-K `inv($K)` substitutions in allosteric `_kcat_forward`

**Files:**
- Modify: `src/rate_eq_derivation.jl:973–991`

- [ ] **Step 1: Delete the substitute_params_expr block.**

Around line 973:

```julia
    # Apply Kd inversion: raw polys use Ka convention, params use Kd
    binding_Ks_R = Set(_binding_K_symbols(CM))
    binding_Ks_T = Set(get(rename_T, K, K) for K in binding_Ks_R)
    num_k_R_expr = substitute_params_expr(
        raw_num_k_R, Dict(K => :(inv($K)) for K in binding_Ks_R))
    den_k_R_expr = substitute_params_expr(
        raw_den_k_R, Dict(K => :(inv($K)) for K in binding_Ks_R))
```

(and the corresponding T-state block immediately below it)

Replace with direct use:

```julia
    num_k_R_expr = raw_num_k_R
    den_k_R_expr = raw_den_k_R
    # T-state expressions also use raw forms now
```

(and similarly drop the T-state substitution wrapping). The polynomial is already in Kd form after Task 1.3 + 1.4.

### Task 1.9: Run dedup tests; confirm Source-B turns green

- [ ] **Step 1: In REPL, re-include the dedup test file.**

```julia
include("test/test_eq_hash_dedup.jl")
```

Expected:
- Sanity tests: green.
- **Source A: still RED** (factoring variants haven't been collapsed yet).
- **Source B: GREEN** ← this commit's target.
- **Source C: still RED** (split-with-tie not absorbed yet).
- Section labels: still RED.

If Source B is still red, the polynomial-construction or Haldane sign-flip is wrong. Diff a Source-B mechanism's `rate_equation_string` output against what you predicted; investigate.

### Task 1.10: Run full test suite

- [ ] **Step 1: Full Pkg.test() run.**

```bash
julia --project -e 'using Pkg; Pkg.test()'
```

Expected: all tests pass *except* Source A, Source C, and section-label tests in `test_eq_hash_dedup.jl` (still red as designed). All hand-updated fixtures from Task 1.2 should now match.

If a non-dedup-related test fails, investigate: most likely the polynomial-flip changes the printed form of some fixture you missed in Task 1.2.

### Task 1.11: Commit 1

- [ ] **Step 1: Stage and commit.**

```bash
git add src/rate_eq_derivation.jl src/sym_poly_for_rate_eq_derivation.jl \
        src/thermodynamic_constr_for_rate_eq_derivation.jl \
        test/mechanism_definitions_for_test_enzyme_derivation.jl \
        test/test_rate_eq_derivation.jl
git commit -m "$(cat <<'EOF'
Kd-by-construction: polynomials store binding K as Kd directly

Source B fix: eliminate the Ka↔Kd inversion layer that produced
1/(1/X) artifacts in Haldane closures. Polynomial construction in
_compute_alpha now puts binding K's in alpha_den (Kd convention);
_dependent_param_exprs sign-flips A-matrix entries for binding K
columns so dep_exprs come out directly in Kd form. _apply_kd_inversion
and the inverted_params/inv_set plumbing are removed.

Source-B dedup tests turn green.

Co-Authored-By: Claude Opus 4.7 (1M context) <noreply@anthropic.com>
EOF
)"
```

---

## Commit 2 — Always-expanded emission + helper deletion (Source A fix)

**Goal of commit:** Drop algebraic factoring of printed numerator/denominator. Expand fully via `_poly_to_expr`. Allosteric MWC outer factoring stays; inner conformation polynomials expand. Delete the `FactoredPoly`/`FactoredSigma`/`DenomTerm` type family and unreachable expansion helpers. Source-A dedup tests turn green.

**Files:**
- Modify: `src/rate_eq_derivation.jl` (drop `_factor_poly`, `_try_algebraic_factor_sigma`, `_try_poly_power`, `_haldane_equality_substitutions`; rework `_rate_v_line`, `_raw_symbolic_rate_polys`, `_kcat_components`, allosteric `_kcat_forward`, `_allosteric_num_den_exprs`)
- Modify: `src/sym_poly_for_rate_eq_derivation.jl` (drop factored types, `_expand_*`, `_poly_power`, `_try_poly_exact_div`, `_estimate_expanded_term_count`, factored-expr helpers, `unfactored_denom_term`)
- Delete: `test/test_sym_poly.jl`
- Modify: `test/runtests.jl` (drop the include)
- Modify: `test/mechanism_definitions_for_test_enzyme_derivation.jl` (update 20 factored-form fixtures)
- Modify: `test/test_rate_eq_derivation.jl` (regenerate byte-identical allosteric fixture)

### Task 2.1: Hand-expand the 20 `expected_factored_*` fixtures (TDD: failing tests come first)

**Files:**
- Modify: `test/mechanism_definitions_for_test_enzyme_derivation.jl`

- [ ] **Step 1: List all fixtures with factored expectations.**

```bash
grep -n "expected_factored_num=\|expected_factored_denom=" test/mechanism_definitions_for_test_enzyme_derivation.jl
```

Expected: 40 lines (20 num + 20 denom). Note line numbers.

- [ ] **Step 2: For each fixture, expand the factored form to a flat sum.**

The factored form looks like:

```julia
expected_factored_num=
    "k_cat * (S - K_eq * P) / K_S",
expected_factored_denom=
    "(1 + S / K_S) * (1 + P / K_P)",
```

Expand the denominator: `(1 + S / K_S) * (1 + P / K_P) = 1 + S / K_S + P / K_P + S * P / (K_S * K_P)`.

Format the expanded string in the same printing convention as `_poly_to_expr` produces:
- Terms separated by ` + ` for positive coefficients, ` - ` for negative.
- Sort terms by (metabolite-degree, lex order of symbol names) — the existing `_poly_to_expr` sort key.
- Coefficient on its own if rational; e.g., `2 * X` not `X * 2`.
- Use `met / K` for binding-K monomials; `K` symbol stays in the denominator of the term.

For the example above, the expanded denominator becomes:

```julia
expected_factored_denom=
    "1 + P / K_P + S / K_S + P * S / (K_P * K_S)",
```

Sort lex: P before S; "P / K_P" then "S / K_S" then "P * S / (K_P * K_S)" (degree-2 last).

Numerator: `k_cat * (S - K_eq * P) / K_S = k_cat * S / K_S - k_cat * K_eq * P / K_S`. Expanded form follows the same sort.

For each of the 20 fixtures, do this transformation. This is mechanical but takes time; allocate ~30 min.

If a fixture's expanded form would exceed visual readability (e.g., 50+ terms), accept that — the test compares strings byte-for-byte regardless of length.

- [ ] **Step 3: Run the affected fixture tests; confirm they FAIL under current code.**

```julia
include("test/test_rate_eq_derivation.jl")
```

Expected: the 20 `expected_factored_*` tests now FAIL because the current code still emits the factored form. This is the TDD discipline — Source A tests in `test_eq_hash_dedup.jl` are also still red.

### Task 2.2: Drop `_factor_poly`, `_try_algebraic_factor_sigma`, `_try_poly_power`

**Files:**
- Modify: `src/rate_eq_derivation.jl:200–427`

- [ ] **Step 1: Delete the three functions.**

In `src/rate_eq_derivation.jl`, around lines 200–427, delete:
- `_try_algebraic_factor_sigma` (~120 lines, line ~210)
- `_try_poly_power` (~20 lines, line ~343)
- `_factor_poly` (~22 lines, line ~405)

Also delete the `# ─── Algebraic Regulator Factoring ─────────` and `# ─── Raw Rate Equation Derivation ───` section comments if they bound only the deleted functions.

### Task 2.2.5: Identify pre-Haldane-merged fixtures so Pkg.test() doesn't surface them mid-commit

**Files:**
- Read: `test/mechanism_definitions_for_test_enzyme_derivation.jl`
- Read: `test/test_rate_eq_derivation.jl`

The `_haldane_equality_substitutions` helper (deleted in Task 2.3 below) currently merges pairs like `k8r → k7r` when both resolve to the same Haldane expression. Removing it before Commit 3's Pass 2 absorption lands means the "merged" fixture form will fail until Commit 3.

- [ ] **Step 1: Grep for fixtures that show this merge pattern.**

```bash
grep -nE "k[0-9]+r" test/mechanism_definitions_for_test_enzyme_derivation.jl test/test_rate_eq_derivation.jl | head -50
```

For each fixture whose Haldane line collapses two `k_r` symbols to one, identify which `k_r` got eliminated. The pattern in current output is: only one of `{k_ar, k_br}` appears in the rate equation v-line; the constraint line writes the surviving symbol.

- [ ] **Step 2: Update affected fixtures to the un-merged form for Commit 2.**

For each, expand the constraint line so both `k_ar` and `k_br` are individually defined (each with its full Haldane RHS — they'll be numerically equal but textually distinct after `_haldane_equality_substitutions` is removed).

This is mechanical but takes judgment: if you can't tell from the fixture which pair was being merged, run the current code on the mechanism in REPL and observe the output first, then update.

- [ ] **Step 3: Verify Pkg.test() at the end of Commit 2 has no remaining surprises from un-merged Haldanes.**

This step is exercised in Task 2.15 below. Pre-handling here prevents that step from cascading into many fixture updates.

- [ ] **Step 4: Note that Commit 3's Pass 2 absorption will re-merge these and the fixtures will need a SECOND update in Task 3.7.**

This is the cost of the interim state. Expected.

### Task 2.3: Drop `_haldane_equality_substitutions`

**Files:**
- Modify: `src/rate_eq_derivation.jl:380–395, 543–550`

- [ ] **Step 1: Delete `_haldane_equality_substitutions`.**

Around line 380:

```julia
"""
Build substitution pairs merging Haldane-derived parameters that have
identical expressions after user constraints. ...
"""
function _haldane_equality_substitutions(dep_exprs)
    ...
end
```

Delete the function.

- [ ] **Step 2: Drop the call site in `_raw_symbolic_rate_polys`.**

Around lines 541–549:

```julia
    # Merge Haldane-derived equal parameters (e.g., k8r→k7r when
    # both resolve to the same thermodynamic expression)
    haldane_subs = _haldane_equality_substitutions(dep_exprs)
    if !isempty(haldane_subs)
        hsub_map = Dict{Symbol, Symbol}(haldane_subs)
        num = _rename_symbols(num, hsub_map)
        denom_terms = [_rename_symbols(dt, hsub_map)
                       for dt in denom_terms]
    end
```

Delete this entire block. Source C in Commit 3 handles this case uniformly via Pass 2 of `_build_kinetic_rename_map`.

### Task 2.4: Change `_raw_symbolic_rate_polys` to return flat POLYs

**Files:**
- Modify: `src/rate_eq_derivation.jl:432–571`

- [ ] **Step 1: Locate the function.**

The function builds `num` (POLY) and `denom_terms::Vector{DenomTerm}`. It needs to return `(num::POLY, den::POLY)` instead.

The denominator is built as `denom_terms[g] = DenomTerm(_factor_poly(csigma, ...), D[g])`. After dropping `_factor_poly`, replace with the flat polynomial: each group's contribution to the denominator is `csigma * D[g]`. Sum across groups for the flat denominator.

- [ ] **Step 2: Rewrite the denominator-construction block.**

Around lines 491–517:

```julia
    # Factor sigma for each RE group
    binding_Ks = Set{Symbol}(
        Symbol("K$i")
        for (i, (lhs, rhs, _, _)) in enumerate(rxns)
        if eq_steps[i] && any(s ∉ enz_set for s in lhs) &&
           all(s ∈ enz_set for s in rhs)
    )
    normalize = G == 1 && sigma_den[1] != poly_one()
    denom_terms = DenomTerm[]
    for g in 1:G
        raw_sigma = if normalize
            reduce(
                poly_add,
                (_poly_div_mono(alpha_num[i], alpha_den[i])
                 for i in groups[g]),
            )
        else
            sigma_num[g]
        end
        csigma = _rename_symbols(raw_sigma, rename_map)
        push!(denom_terms, DenomTerm(
            _factor_poly(
                csigma, rxns, eq_steps, enz_set, rename_map;
                binding_Ks,
            ),
            D[g],
        ))
    end
```

Replace with:

```julia
    normalize = G == 1 && sigma_den[1] != poly_one()
    den = poly_zero()
    for g in 1:G
        raw_sigma = if normalize
            reduce(
                poly_add,
                (_poly_div_mono(alpha_num[i], alpha_den[i])
                 for i in groups[g]),
            )
        else
            sigma_num[g]
        end
        csigma = _rename_symbols(raw_sigma, rename_map)
        den = poly_add(den, poly_mul(csigma, D[g]))
    end
```

- [ ] **Step 3: Rewrite the post-derivation cleanup block.**

Around lines 525–571:

```julia
    abs_nu = abs(nu_ref)
    if abs_nu != 1
        for (i, dt) in enumerate(denom_terms)
            denom_terms[i] = DenomTerm(
                dt.sigma,
                poly_mul(dt.cofactor, poly_const(abs_nu)),
            )
        end
    end

    # Apply kinetic-group renaming (K2 → K1 etc.) to numerator
    num = _rename_symbols(num, rename_map)
    denom_terms = [_rename_symbols(dt, rename_map) for dt in denom_terms]

    [Haldane equality block already deleted in Task 2.3]

    n_terms = (length(num) +
               _estimate_expanded_term_count(denom_terms))
    if n_terms > MAX_RATE_EQUATION_TERMS
        error(...)
    end

    # Factor numerator (only if it reduces display terms)
    num_fs = _factor_poly(
        num, rxns, eq_steps, enz_set, rename_map;
        binding_Ks, check_benefit=true,
    )

    num_fs, denom_terms
end
```

Replace with:

```julia
    abs_nu = abs(nu_ref)
    abs_nu != 1 && (den = poly_mul(den, poly_const(abs_nu)))

    # Apply kinetic-group renaming (K2 → K1 etc.) to numerator and denominator
    num = _rename_symbols(num, rename_map)
    den = _rename_symbols(den, rename_map)

    n_terms = length(num) + length(den)
    if n_terms > MAX_RATE_EQUATION_TERMS
        error(
            "Rate equation for this mechanism has $n_terms polynomial " *
            "terms (limit: $MAX_RATE_EQUATION_TERMS). Equations this " *
            "large take a very long time to compile and are unlikely " *
            "to be practically useful for parameter fitting.",
        )
    end

    num, den
end
```

- [ ] **Step 4: Update the call site type annotations.**

Look at the function's documentation comments and any return-type annotations earlier in the file. Update to reflect `(num::POLY, den::POLY)` rather than `(num_fs::FactoredSigma, denom_terms::Vector{DenomTerm})`.

### Task 2.5: Simplify `_rate_v_line` to flat-POLY emission

**Files:**
- Modify: `src/rate_eq_derivation.jl:749–763`

- [ ] **Step 1: Replace with the four-line version.**

Around lines 749–763:

```julia
function _rate_v_line(M::Type{<:EnzymeMechanism})
    num, den = _raw_symbolic_rate_polys(M)
    m = M()
    ps = Set{Symbol}(_raw_param_symbols(m))
    cs = Set{Symbol}(metabolites(m))
    "v = E_total * ($(_expr_to_string(_poly_to_expr(num, ps, cs)))) / " *
        "($(_expr_to_string(_poly_to_expr(den, ps, cs))))"
end
```

### Task 2.6: Simplify `_raw_rate_expr_and_symbols` and `to_rate_expr`

**Files:**
- Modify: `src/rate_eq_derivation.jl:702–711`
- Modify: `src/sym_poly_for_rate_eq_derivation.jl:447–462`

- [ ] **Step 1: Update `_raw_rate_expr_and_symbols`.**

Around lines 702–711:

```julia
function _raw_rate_expr_and_symbols(M::Type{<:EnzymeMechanism})
    num, denom_terms = _raw_symbolic_rate_polys(M)
    m = M()
    param_syms = Set{Symbol}(_raw_param_symbols(m))
    conc_syms = Set{Symbol}(metabolites(m))
    expr = to_rate_expr(num, denom_terms, param_syms, conc_syms)
    all_params = _sorted_raw_param_symbols(M)
    return expr, all_params, metabolites(m)
end
```

The variable name `denom_terms` is now misleading (it's a flat POLY). Rename to `den`:

```julia
function _raw_rate_expr_and_symbols(M::Type{<:EnzymeMechanism})
    num, den = _raw_symbolic_rate_polys(M)
    m = M()
    param_syms = Set{Symbol}(_raw_param_symbols(m))
    conc_syms = Set{Symbol}(metabolites(m))
    expr = to_rate_expr(num, den, param_syms, conc_syms)
    all_params = _sorted_raw_param_symbols(M)
    return expr, all_params, metabolites(m)
end
```

- [ ] **Step 2: Replace `to_rate_expr` with the flat-POLY version.**

In `src/sym_poly_for_rate_eq_derivation.jl` around lines 447–462:

```julia
function to_rate_expr(
    num::Union{POLY, FactoredSigma},
    denom_terms::Vector{DenomTerm},
    param_syms::Set{Symbol}, conc_syms::Set{Symbol},
    inverted_params::Set{Symbol}=Set{Symbol}(),
)
    num_expr = num isa POLY ?
        _poly_to_expr(num, param_syms, conc_syms, inverted_params) :
        _factored_sigma_to_expr(
            num, param_syms, conc_syms, inverted_params,
        )
    den_expr = _denom_terms_to_expr(
        denom_terms, param_syms, conc_syms, inverted_params,
    )
    :(E_total * ($num_expr) / ($den_expr))
end
```

Replace with:

```julia
function to_rate_expr(
    num::POLY, den::POLY,
    param_syms::Set{Symbol}, conc_syms::Set{Symbol},
)
    num_expr = _poly_to_expr(num, param_syms, conc_syms)
    den_expr = _poly_to_expr(den, param_syms, conc_syms)
    :(E_total * ($num_expr) / ($den_expr))
end
```

### Task 2.7: Update `_kcat_components` to consume flat POLYs

**Files:**
- Modify: `src/rate_eq_derivation.jl:858–875`

- [ ] **Step 1: Drop the `_expand_*` calls.**

Around line 858:

```julia
function _kcat_components(M::Type{<:EnzymeMechanism})
    num_fs, denom_terms = _raw_symbolic_rate_polys(M)
    num = _expand_factored_sigma(num_fs)
    den = _expand_to_poly(denom_terms)
    num_groups, den_groups = _kcat_groups_from_polys(num, den)
    ...
```

Replace with:

```julia
function _kcat_components(M::Type{<:EnzymeMechanism})
    num, den = _raw_symbolic_rate_polys(M)
    num_groups, den_groups = _kcat_groups_from_polys(num, den)
    ...
```

### Task 2.8: Update allosteric `_kcat_forward` to consume flat POLYs

**Files:**
- Modify: `src/rate_eq_derivation.jl:929–991`

- [ ] **Step 1: Replace the `_expand_*` calls with direct flat-POLY use.**

Around line 943:

```julia
    num_fs, denom_terms = _raw_symbolic_rate_polys(CM)
    r_only_syms = _onlyR_syms(m)
    rename_T = _T_rename(m)
    num_R_poly = _expand_factored_sigma(num_fs)
    den_R_poly = _expand_to_poly(denom_terms)
    num_T_poly = _rename_symbols(
        _zero_symbols_in_poly(_expand_factored_sigma(num_fs), r_only_syms),
        rename_T)
    den_T_poly = _rename_symbols(
        _zero_symbols_in_poly(_expand_to_poly(denom_terms), r_only_syms),
        rename_T)
```

Replace with:

```julia
    num_R_poly, den_R_poly = _raw_symbolic_rate_polys(CM)
    r_only_syms = _onlyR_syms(m)
    rename_T = _T_rename(m)
    num_T_poly = _rename_symbols(
        _zero_symbols_in_poly(num_R_poly, r_only_syms), rename_T)
    den_T_poly = _rename_symbols(
        _zero_symbols_in_poly(den_R_poly, r_only_syms), rename_T)
```

### Task 2.9: Update `_allosteric_num_den_exprs` to consume flat POLYs

**Files:**
- Modify: `src/rate_eq_derivation.jl:1526–1611`

- [ ] **Step 1: Replace factored-expr emission with `_poly_to_expr`.**

Around lines 1526–1572:

```julia
function _allosteric_num_den_exprs(M_type::Type{<:AllostericEnzymeMechanism})
    m = M_type()
    CM = typeof(catalytic_mechanism(m))
    CatN = catalytic_multiplicity(m)
    RS = regulatory_sites(m)

    num_fs, denom_terms = _raw_symbolic_rate_polys(CM)
    cat_params = Set{Symbol}(_raw_param_symbols(CM()))
    cat_mets = Set{Symbol}(metabolites(CM()))
    binding_Ks_r = Set(_binding_K_symbols(CM))

    r_only_syms = _onlyR_syms(m)
    rename_T = _T_rename(m)
    binding_Ks_t = Set(get(rename_T, K, K) for K in binding_Ks_r)

    N_R = _factored_sigma_to_expr(num_fs, cat_params, cat_mets, binding_Ks_r)
    Q_R = _denom_terms_to_expr(denom_terms, cat_params, cat_mets, binding_Ks_r)

    if isempty(r_only_syms)
        N_T = substitute_params_expr(
            _factored_sigma_to_expr(num_fs, cat_params, cat_mets, binding_Ks_r),
            rename_T)
        Q_T = substitute_params_expr(
            _denom_terms_to_expr(denom_terms, cat_params, cat_mets, binding_Ks_r),
            rename_T)
    else
        num_t_poly = _rename_symbols(
            _zero_symbols_in_poly(_expand_factored_sigma(num_fs), r_only_syms),
            rename_T)
        den_t_poly = _rename_symbols(
            _zero_symbols_in_poly(_expand_to_poly(denom_terms), r_only_syms),
            rename_T)
        N_T = _t_state_dead(m) ? 0 :
              _poly_to_expr(num_t_poly, cat_params, cat_mets, binding_Ks_t)
        Q_T = _poly_to_expr(den_t_poly, cat_params, cat_mets, binding_Ks_t)
    end
    ...
```

Replace with:

```julia
function _allosteric_num_den_exprs(M_type::Type{<:AllostericEnzymeMechanism})
    m = M_type()
    CM = typeof(catalytic_mechanism(m))
    CatN = catalytic_multiplicity(m)
    RS = regulatory_sites(m)

    num_R_poly, den_R_poly = _raw_symbolic_rate_polys(CM)
    cat_params = Set{Symbol}(_raw_param_symbols(CM()))
    cat_mets = Set{Symbol}(metabolites(CM()))

    r_only_syms = _onlyR_syms(m)
    rename_T = _T_rename(m)

    N_R = _poly_to_expr(num_R_poly, cat_params, cat_mets)
    Q_R = _poly_to_expr(den_R_poly, cat_params, cat_mets)

    if isempty(r_only_syms)
        N_T = substitute_params_expr(N_R, rename_T)
        Q_T = substitute_params_expr(Q_R, rename_T)
    else
        num_t_poly = _rename_symbols(
            _zero_symbols_in_poly(num_R_poly, r_only_syms), rename_T)
        den_t_poly = _rename_symbols(
            _zero_symbols_in_poly(den_R_poly, r_only_syms), rename_T)
        N_T = _t_state_dead(m) ? 0 :
              _poly_to_expr(num_t_poly, cat_params, cat_mets)
        Q_T = _poly_to_expr(den_t_poly, cat_params, cat_mets)
    end
    ...
```

(Keep the rest of the function — `make_num_term`, `make_den_term`, `reg_Q_R`, etc. — unchanged. Drop only the `binding_Ks_r`/`binding_Ks_t` declarations.)

### Task 2.10: Drop `FactoredPoly`, `FactoredSigma`, `DenomTerm` types and their helpers

**Files:**
- Modify: `src/sym_poly_for_rate_eq_derivation.jl:344–525`

- [ ] **Step 1: Delete the type definitions and `unfactored_denom_term`.**

Around lines 344–370. Delete:

```julia
# ─── Factored denominator types ──────────────────────────────

"""Product of POLY factors with integer exponents: ..."""
struct FactoredPoly
    factors::Vector{POLY}
    exponents::Vector{Int}
end

struct FactoredSigma
    coefficients::Vector{POLY}
    products::Vector{FactoredPoly}
end

struct DenomTerm
    sigma::FactoredSigma
    cofactor::POLY
end

function unfactored_denom_term(sigma_num::POLY, cofactor::POLY)
    DenomTerm(...)
end
```

- [ ] **Step 2: Delete the factored-expr helpers.**

Around lines 372–443. Delete:
- `_factored_poly_to_expr`
- `_factored_sigma_to_expr`
- `_denom_terms_to_expr`

- [ ] **Step 3: Delete the factored `_rename_symbols` overloads.**

Around lines 466–487. Delete the three overloads for `FactoredPoly`, `FactoredSigma`, `DenomTerm`.

- [ ] **Step 4: Delete the expansion helpers.**

Around lines 489–544. Delete:
- `_expand_factored_poly`
- `_expand_factored_sigma`
- `_expand_to_poly`
- `_estimate_expanded_term_count`

### Task 2.11: Drop `_poly_power` and `_try_poly_exact_div`

**Files:**
- Modify: `src/sym_poly_for_rate_eq_derivation.jl:48–119`

- [ ] **Step 1: Delete `_poly_power`.**

Around lines 47–55. Delete:

```julia
"""Raise a POLY to a non-negative integer power via repeated multiplication."""
function _poly_power(p::POLY, n::Int)
    n == 0 && return poly_one()
    result = poly_one()
    for _ in 1:n
        result = poly_mul(result, p)
    end
    result
end
```

- [ ] **Step 2: Delete `_try_poly_exact_div`.**

Around lines 64–116. Delete the function (it's the largest single function in the file at ~50 lines).

### Task 2.12: Delete `test/test_sym_poly.jl` and remove its include

**Files:**
- Delete: `test/test_sym_poly.jl`
- Modify: `test/runtests.jl`

- [ ] **Step 1: Delete the test file.**

```bash
git rm test/test_sym_poly.jl
```

- [ ] **Step 2: Remove the include from `test/runtests.jl`.**

In `test/runtests.jl`, delete the line:

```julia
    include("test_sym_poly.jl")
```

### Task 2.13: Re-derive byte-identical allosteric fixture

**Files:**
- Modify: `test/test_rate_eq_derivation.jl:1142–1166`

- [ ] **Step 1: Run the test to capture the current output.**

```julia
using EnzymeRates: rate_equation_string
include("test/test_rate_eq_derivation.jl")
```

The byte-identical allosteric test will fail. Before fixing, capture what the new code produces:

```julia
rxn_allo = @enzyme_reaction begin
    substrates: S[C]
    products:   P[C]
    allosteric_regulators: R
    oligomeric_state: 2
end
init = EnzymeRates.init_mechanisms(rxn_allo)
base = first(init)
used_groups = sort!(collect(Set(s.kinetic_group for s in base.steps)))
spec = EnzymeRates.AllostericMechanismSpec(
    base, 2, [[:R]], [2],
    Dict(g => :NonequalRT for g in used_groups),
    Dict(:R => :NonequalRT),
    base.n_fit_params_estimate + 5)
m_allo = EnzymeRates.AllostericEnzymeMechanism(spec)
println(rate_equation_string(m_allo))
```

- [ ] **Step 2: Verify numerical equivalence with the existing `test_rate_equation_string` infrastructure.**

The fixture's `expected` is byte-exact, but the same testset has a numerical-equivalence assertion via `_eval_rate_string`. Before committing the regenerated string, manually verify the output is mathematically equivalent — same v/(num/den), Kd convention, fully expanded inner polynomials.

- [ ] **Step 3: Replace the `expected` raw string with the captured output.**

Update lines 1160–1164 of `test/test_rate_eq_derivation.jl` to embed the captured string verbatim.

### Task 2.14: Run dedup tests; confirm Source A turns green

- [ ] **Step 1: Run the dedup test file.**

```julia
include("test/test_eq_hash_dedup.jl")
```

Expected:
- Sanity: green.
- Source A: **GREEN** ← this commit's target.
- Source B: green (from Commit 1).
- Source C: still RED.
- Section labels: still RED.

### Task 2.15: Run full test suite

- [ ] **Step 1: Full Pkg.test() run.**

```bash
julia --project -e 'using Pkg; Pkg.test()'
```

Expected: all tests pass except Source C and section-label tests in `test_eq_hash_dedup.jl`. The 20 hand-expanded fixtures should match. The byte-identical allosteric fixture should match (after Task 2.13).

If a non-dedup test fails: most likely a fixture you missed, or the polynomial-renaming order differs in `_rename_symbols`. Investigate; do not bypass.

### Task 2.15.5: Update CLAUDE.md line 290 (sym_poly bullet)

**Files:**
- Modify: `.claude/CLAUDE.md:290`

- [ ] **Step 1: Find the current bullet.**

Around line 290:

```markdown
- `src/sym_poly_for_rate_eq_derivation.jl` — Symbolic polynomial algebra (`POLY`/`FactoredSigma`/`FactoredPoly`/`DenomTerm`); `_rename_symbols`, `_zero_symbols_in_poly` for MWC allosteric-state-driven substitution.
```

- [ ] **Step 2: Replace with a description that no longer mentions the deleted types.**

```markdown
- `src/sym_poly_for_rate_eq_derivation.jl` — Symbolic polynomial algebra (`POLY` = `Dict{MONO, Rational{Int}}`); arithmetic ops (`poly_add`, `poly_mul`, `sym_det`, etc.); `_rename_symbols` and `_zero_symbols_in_poly` for MWC allosteric-state-driven substitution; `_poly_to_expr` and `_expr_to_string` for emission.
```

### Task 2.16: Commit 2

- [ ] **Step 1: Stage and commit.**

```bash
git add src/rate_eq_derivation.jl src/sym_poly_for_rate_eq_derivation.jl \
        test/test_rate_eq_derivation.jl \
        test/mechanism_definitions_for_test_enzyme_derivation.jl \
        test/runtests.jl
git rm test/test_sym_poly.jl
git commit -m "$(cat <<'EOF'
Always-expanded rate equation emission; drop factored-poly type family

Source A fix: drop algebraic factoring of printed numerator/denominator
for non-allosteric and allosteric inner conformation polynomials.
_raw_symbolic_rate_polys now returns (num::POLY, den::POLY) directly.

Removes _factor_poly, _try_algebraic_factor_sigma, _try_poly_power,
_haldane_equality_substitutions, FactoredPoly/FactoredSigma/DenomTerm
types and their _rename_symbols overloads, _factored_*_to_expr,
_expand_factored_*, _expand_to_poly, _estimate_expanded_term_count,
_poly_power, _try_poly_exact_div. Deletes the entire test/test_sym_poly.jl
(was 299 lines testing only the now-removed types).

The MWC outer factoring (E_total * (N_R + L*N_T) / (Q_R^N + L*Q_T^N))
is structurally required and stays; only the inner conformation
polynomials expand.

Source-A dedup tests turn green.

Co-Authored-By: Claude Opus 4.7 (1M context) <noreply@anthropic.com>
EOF
)"
```

---

## Commit 3 — Section labels + Wegscheider absorption + canonicalizer normalization (Source C fix)

**Goal of commit:** Single-symbol Wegscheider RE ties get absorbed into the kinetic-group rename map at polynomial level; `rate_equation_string` emits constraints under three section headers; eq_hash canonicalizer normalizes single-symbol equality lines so Source C duplicates collapse to one hash. Source-C and section-label dedup tests turn green; CSV-replay test passes.

**Files:**
- Modify: `src/thermodynamic_constr_for_rate_eq_derivation.jl` (extract `_dependent_param_exprs_kernel`)
- Modify: `src/rate_eq_derivation.jl` (Pass 2 of `_build_kinetic_rename_map`, sectioned `rate_equation_string`, `ANNOTATION_SUBSTITUTED`)
- Modify: `src/identify_rate_equation.jl` (canonicalizer block normalization)
- Create: `test/test_dedup_csv_replay.jl`
- Modify: `test/runtests.jl` (add include for CSV-replay test)
- Modify: existing fixtures that show constraint lines (add `(substituted into v)` annotation + section headers)

### Task 3.1: Extract `_dependent_param_exprs_kernel`

**Files:**
- Modify: `src/thermodynamic_constr_for_rate_eq_derivation.jl:161–290`

- [ ] **Step 1: Refactor `_dependent_param_exprs` into a kernel + wrapper pair.**

The current function takes `M::Type{<:EnzymeMechanism}` and uses `_build_kinetic_rename_map(M())` internally (line 192). Extract a kernel that accepts the rename map as a parameter, so we can call it twice with different renames.

Replace `_dependent_param_exprs(M)` with:

```julia
"""
Internal kernel for Gaussian-elimination-based Wegscheider/Haldane resolution.
Takes the rename map as a parameter so callers can supply either the
user-defined kinetic-group rename (Pass 1 only) or the full rename map
including absorbed single-symbol Wegscheider ties (Pass 1 + Pass 2).
"""
function _dependent_param_exprs_kernel(
    M::Type{<:EnzymeMechanism},
    rename::AbstractDict{Symbol, Symbol},
)
    m = M()
    rxns = reactions(m)
    eq_steps = equilibrium_steps(m)
    enz_names = enzyme_forms(m)
    enz_set = Set(enz_names)

    # ... existing body of _dependent_param_exprs starting after line 192,
    # but using the `rename` parameter instead of calling _build_kinetic_rename_map ...
end

"""
    _dependent_param_exprs(M::Type{<:EnzymeMechanism}) → (dep_exprs, indep)

Wrapper that uses the full rename map (user-defined kinetic groups +
absorbed single-symbol Wegscheider RE ties).
"""
function _dependent_param_exprs(M::Type{<:EnzymeMechanism})
    rename = _build_kinetic_rename_map(M)
    _dependent_param_exprs_kernel(M, rename)
end
```

The body of the kernel is identical to the current `_dependent_param_exprs` body except the line `rename = _build_kinetic_rename_map(m)` is replaced by reading the `rename` parameter.

- [ ] **Step 2: Verify nothing else changed.**

```bash
julia --project -e 'using EnzymeRates; m = EnzymeRates.Tests; nothing'
# Or just rely on Revise reload + run an existing test
```

The wrapper preserves the existing public interface; nothing else should break.

### Task 3.1.5: Inline the 2-arg overloads of `_enzyme_incidence_matrix` and `_thermodynamic_constraints`

**Files:**
- Modify: `src/thermodynamic_constr_for_rate_eq_derivation.jl:35–122`

The current chain has 4 functions across 2 indirection layers. The 1-arg `_enzyme_incidence_matrix(M::Type)` overload (lines 47–53) has zero callers as currently written. The multi-arg `_enzyme_incidence_matrix` and `_thermodynamic_constraints` versions can be inlined into a single `_thermodynamic_constraints(M::Type)` function. Saves ~30 lines.

- [ ] **Step 1: Delete the 1-arg `_enzyme_incidence_matrix(M::Type)` overload.**

Around lines 47–53:

```julia
function _enzyme_incidence_matrix(M::Type{<:EnzymeMechanism})
    m = M()
    enz_names = enzyme_forms(m)
    _enzyme_incidence_matrix(
        enz_names, Set(enz_names), reactions(m),
    )
end
```

Delete the whole block. Confirm zero callers via `grep -n "_enzyme_incidence_matrix(M" src/`.

- [ ] **Step 2: Inline the multi-arg `_enzyme_incidence_matrix` body into `_thermodynamic_constraints`.**

Around lines 35–45 (multi-arg `_enzyme_incidence_matrix`) and 91–112 (multi-arg `_thermodynamic_constraints`), merge by inlining the matrix construction into the multi-arg `_thermodynamic_constraints` body. Keep the helper as a closure or local block — either is fine.

- [ ] **Step 3: Inline the multi-arg `_thermodynamic_constraints` body into the 1-arg version (which is the only caller).**

Around lines 114–122. The 1-arg version becomes the single function; the multi-arg version goes away.

After this Step the file has:
- `_integer_nullspace` (one function, unchanged)
- `_thermodynamic_constraints(M::Type{<:EnzymeMechanism})` (one function, ~50 lines, including inlined matrix construction)
- `_classify_cycle` (one function, unchanged)

- [ ] **Step 4: Reload, run a small mechanism test as a smoke check.**

```julia
include("test/mechanism_definitions_for_test_enzyme_derivation.jl")
m = first(specs).mechanism  # any spec
EnzymeRates._dependent_param_exprs(typeof(m))
```

Expected: same return values as before. The inlining is purely structural; no behavior changes.

### Task 3.2: Add Pass 2 to `_build_kinetic_rename_map`

**Files:**
- Modify: `src/rate_eq_derivation.jl:90–108`

- [ ] **Step 1: Add the second pass that absorbs single-symbol Wegscheider RE ties.**

Replace the existing function:

```julia
function _build_kinetic_rename_map(m::EnzymeMechanism)
    rename = Dict{Symbol, Symbol}()
    eq = equilibrium_steps(m)
    for g in kinetic_groups(m)
        idxs = steps_in_group(m, g)
        length(idxs) == 1 && continue
        rep = first(idxs)
        for idx in idxs
            idx == rep && continue
            if eq[idx]
                rename[Symbol("K$idx")] = Symbol("K$rep")
            else
                rename[Symbol("k$(idx)f")] = Symbol("k$(rep)f")
                rename[Symbol("k$(idx)r")] = Symbol("k$(rep)r")
            end
        end
    end
    rename
end
```

With:

```julia
"""Predicate: is this Symbol an RE binding K (i.e., `K{digits}` form)?"""
function _is_re_K(sym::Symbol)
    s = string(sym)
    length(s) >= 2 && s[1] == 'K' && all(isdigit, s[2:end])
end

function _build_kinetic_rename_map(M::Type{<:EnzymeMechanism})
    m = M()
    rename = Dict{Symbol, Symbol}()
    eq = equilibrium_steps(m)

    # Pass 1: user-defined kinetic-group merges
    for g in kinetic_groups(m)
        idxs = steps_in_group(m, g)
        length(idxs) == 1 && continue
        rep = first(idxs)
        for idx in idxs
            idx == rep && continue
            if eq[idx]
                rename[Symbol("K$idx")] = Symbol("K$rep")
            else
                rename[Symbol("k$(idx)f")] = Symbol("k$(rep)f")
                rename[Symbol("k$(idx)r")] = Symbol("k$(rep)r")
            end
        end
    end

    # Pass 2: single-symbol Wegscheider RE ties from raw dep_exprs (Pass 1 only)
    dep_raw, _ = _dependent_param_exprs_kernel(M, rename)
    for (lhs, rhs) in dep_raw
        rhs isa Symbol || continue
        _is_re_K(lhs) && _is_re_K(rhs) || continue
        rename[lhs] = get(rename, rhs, rhs)
    end

    rename
end

# Backward-compatible 1-arg overload taking a mechanism instance
_build_kinetic_rename_map(m::EnzymeMechanism) = _build_kinetic_rename_map(typeof(m))
```

- [ ] **Step 2: Audit callers of `_build_kinetic_rename_map`.**

```bash
grep -n "_build_kinetic_rename_map" src/*.jl test/*.jl
```

Update any caller that passes an instance to pass `typeof(m)` instead — or rely on the backward-compatible overload above. Either is fine; preference is for the new type-based call signature.

### Task 3.3: Add `ANNOTATION_SUBSTITUTED` constant and restructure `rate_equation_string` for sections

**Files:**
- Modify: `src/rate_eq_derivation.jl` (around `rate_equation_string` for `EnzymeMechanism` and `AllostericEnzymeMechanism`)

- [ ] **Step 1: Add the constant near the top of `rate_eq_derivation.jl`.**

After the `_AnyMechanism` alias (line 3):

```julia
"""Suffix appended to single-symbol equality lines whose LHS got folded
into the kinetic-group rename map. Both display sites (User defined
kinetic-group merges and absorbed single-symbol Wegscheider ties) use
this exact string so the eq_hash canonicalizer normalizes consistently."""
const ANNOTATION_SUBSTITUTED = "  (substituted into v)"
```

- [ ] **Step 2: Rewrite `rate_equation_string(::M, ::ReducedMode)` for `EnzymeMechanism`.**

Around lines 772–778:

```julia
function rate_equation_string(::M, ::ReducedMode) where {M<:EnzymeMechanism}
    _, indep = _dependent_param_exprs(M)
    join(["(; $(join((indep..., :Keq, :E_total), ", "))) = params",
          "(; $(join(metabolites(M()), ", "))) = concs",
          _constraint_expr_strings(M)...,
          _rate_v_line(M)], "\n")
end
```

Replace with:

```julia
function rate_equation_string(::M, ::ReducedMode) where {M<:EnzymeMechanism}
    m = M()
    _, indep = _dependent_param_exprs(M)

    # User defined: equalities encoded in kinetic-group structure (type info)
    user_lines = String[]
    eq = equilibrium_steps(m)
    for g in kinetic_groups(m)
        idxs = steps_in_group(m, g)
        length(idxs) == 1 && continue
        rep = first(idxs)
        for idx in idxs
            idx == rep && continue
            if eq[idx]
                push!(user_lines, "K$idx = K$rep$ANNOTATION_SUBSTITUTED")
            else
                push!(user_lines, "k$(idx)f = k$(rep)f$ANNOTATION_SUBSTITUTED")
                push!(user_lines, "k$(idx)r = k$(rep)r$ANNOTATION_SUBSTITUTED")
            end
        end
    end

    # Wegscheider/Haldane: from raw dep_exprs (Pass-1 rename only).
    pass1_rename = Dict{Symbol, Symbol}()
    for g in kinetic_groups(m)
        idxs = steps_in_group(m, g)
        length(idxs) == 1 && continue
        rep = first(idxs)
        for idx in idxs
            idx == rep && continue
            eq[idx] ?
                (pass1_rename[Symbol("K$idx")] = Symbol("K$rep")) :
                (pass1_rename[Symbol("k$(idx)f")] = Symbol("k$(rep)f");
                 pass1_rename[Symbol("k$(idx)r")] = Symbol("k$(rep)r"))
        end
    end
    dep_raw, _ = _dependent_param_exprs_kernel(M, pass1_rename)

    weg_lines, hal_lines = String[], String[]
    keq_set = Set([:Keq])
    for (sym, expr) in sort(collect(dep_raw); by=p -> string(p[1]))
        is_haldane = _expr_references_any(expr, keq_set)
        is_single = expr isa Symbol
        suffix = is_single ? ANNOTATION_SUBSTITUTED : ""
        line = "$sym = $(string(expr))$suffix"
        push!(is_haldane ? hal_lines : weg_lines, line)
    end

    lines = ["(; $(join((indep..., :Keq, :E_total), ", "))) = params",
             "(; $(join(metabolites(m), ", "))) = concs"]
    isempty(user_lines) || (push!(lines, "# User defined constraints:"); append!(lines, user_lines))
    isempty(weg_lines)  || (push!(lines, "# Wegscheider constraints:");  append!(lines, weg_lines))
    isempty(hal_lines)  || (push!(lines, "# Haldane constraints:");      append!(lines, hal_lines))
    push!(lines, _rate_v_line(M))
    join(lines, "\n")
end
```

### Task 3.4: Update allosteric `rate_equation_string` for sections

**Files:**
- Modify: `src/rate_eq_derivation.jl:1647–1670`

- [ ] **Step 1: Apply the same partition logic.**

Update the allosteric overload to:
1. Compute `user_lines` from kinetic-group structure of `catalytic_mechanism(m)`.
2. Use `_dependent_param_exprs_kernel(CM, pass1_rename)` to get raw R-state deps.
3. Partition into Wegscheider vs Haldane via `_expr_references_any(expr, Set([:Keq]))`.
4. Build T-state equivalent lines via existing `_build_dep_assignments`.
5. Emit under the three section headers.

Skeleton:

```julia
function rate_equation_string(
    ::AllostericEnzymeMechanism{CM,CS,RS}, ::ReducedMode,
) where {CM,CS,RS}
    M = AllostericEnzymeMechanism{CM,CS,RS}
    m = M()
    cm = catalytic_mechanism(m)
    _, indep = _dependent_param_exprs(M)
    hw_params = (indep..., :Keq, :E_total)
    mets = metabolites(m)

    # User defined: from catalytic-mechanism kinetic-group structure
    user_lines = String[]
    eq = equilibrium_steps(cm)
    for g in kinetic_groups(cm)
        idxs = steps_in_group(cm, g)
        length(idxs) == 1 && continue
        rep = first(idxs)
        for idx in idxs
            idx == rep && continue
            if eq[idx]
                push!(user_lines, "K$idx = K$rep$ANNOTATION_SUBSTITUTED")
            else
                push!(user_lines, "k$(idx)f = k$(rep)f$ANNOTATION_SUBSTITUTED")
                push!(user_lines, "k$(idx)r = k$(rep)r$ANNOTATION_SUBSTITUTED")
            end
        end
    end

    # Wegscheider/Haldane R-state (raw): same partition logic as non-allosteric
    pass1_rename = ...  # build same as in non-allosteric overload
    dep_R_raw, _ = _dependent_param_exprs_kernel(typeof(cm), pass1_rename)

    weg_lines, hal_lines = String[], String[]
    keq_set = Set([:Keq])
    for (sym, expr) in sort(collect(dep_R_raw); by=p -> string(p[1]))
        is_haldane = _expr_references_any(expr, keq_set)
        is_single = expr isa Symbol
        suffix = is_single ? ANNOTATION_SUBSTITUTED : ""
        line = "$sym = $(string(expr))$suffix"
        push!(is_haldane ? hal_lines : weg_lines, line)
    end

    # T-state assignments from existing _build_dep_assignments logic;
    # partition into Wegscheider (no Keq) vs Haldane (has Keq) and append
    # to the corresponding section so the spec's three-section design
    # holds for allosteric mechanisms too.
    _, t_assignments_ = _build_dep_assignments(M)
    t_assignments = _t_state_dead(m) ? Expr[] : t_assignments_
    keq_set = Set([:Keq])
    for a in t_assignments
        sym = a.args[1]
        expr = a.args[2]
        is_haldane = _expr_references_any(expr, keq_set)
        line = "$sym = $(_expr_to_string(expr))"
        push!(is_haldane ? hal_lines : weg_lines, line)
    end

    full_num, full_den = _allosteric_num_den_exprs(M)
    v_line = "v = E_total * ($(_expr_to_string(full_num))) / ($(_expr_to_string(full_den)))"

    lines = [
        "(; $(join(hw_params, ", "))) = params",
        "(; $(join(mets, ", "))) = concs",
    ]
    isempty(user_lines) || (push!(lines, "# User defined constraints:"); append!(lines, user_lines))
    isempty(weg_lines)  || (push!(lines, "# Wegscheider constraints:");  append!(lines, weg_lines))
    isempty(hal_lines)  || (push!(lines, "# Haldane constraints:");      append!(lines, hal_lines))
    push!(lines, v_line)
    join(lines, "\n")
end
```

(T-state assignments fold into the Wegscheider/Haldane sections by Keq-reference predicate so allosteric output uses the same three-section structure as non-allosteric. Cosmetically, T-state lines may appear interleaved with R-state lines within each section — re-sort each section's lines lexicographically before emission if interleaving matters; otherwise rely on the canonicalizer's block-normalization for hashing.)

### Task 3.5: Update `_build_rate_body` to use post-absorption dep_exprs

**Files:**
- Modify: `src/thermodynamic_constr_for_rate_eq_derivation.jl:344–357`

- [ ] **Step 1: Confirm runtime body uses `_dependent_param_exprs` (post-absorption).**

The current `_build_rate_body(M, ::Type{ReducedMode})` calls `_dependent_param_exprs(M)`, which after Task 3.1+3.2 returns dep_exprs that exclude the absorbed single-symbol Wegscheider ties. Single-symbol absorbed ties are folded into the polynomial via the rename map and don't need runtime `K_a = K_b` assignment statements.

Verify the function at line 344:

```julia
function _build_rate_body(M, ::Type{ReducedMode})
    expr, _, conc_syms = _raw_rate_expr_and_symbols(M)
    dep_exprs, indep = _dependent_param_exprs(M)
    hw_params = (indep..., :Keq, :E_total)
    assignments = [Expr(:(=), sym, dep_exprs[sym])
                   for (sym, _) in sort(collect(dep_exprs); by=first)]
    Expr(:block,
        _destructuring_expr(hw_params, :params),
        _destructuring_expr(conc_syms, :concs),
        assignments...,
        expr)
end
```

This is unchanged structurally — just the dep_exprs now contain fewer entries. No edit needed here unless TODOs from Task 1.5 left any `_apply_kd_inversion` reference.

### Task 3.6: Add canonicalizer block normalization

**Files:**
- Modify: `src/identify_rate_equation.jl:151–191`

- [ ] **Step 1: Locate `_canonicalize_rate_eq_with_map`.**

The function around line 151 currently:
1. Strips destructure header lines.
2. Walks `parameters(m, Full)`, scans body for first-appearance positions, alpha-renames to `p_1, p_2, …`.

Add a new normalization step BEFORE the alpha-rename pass:
1. Strip section header lines (`# User defined constraints:`, `# Wegscheider constraints:`, `# Haldane constraints:`, `# T-state constraints:`).
2. Identify single-symbol equality lines (those with `(substituted into v)` annotation).
3. Collect them, sort lexicographically, re-emit as a single contiguous block at a fixed canonical position (immediately after the `concs` destructure line).

- [ ] **Step 2: Implement.**

Replace the body-prep section (around lines 156–164):

```julia
    # Strip ONLY the destructure header lines.
    body = join(
        filter(
            ln -> !occursin(
                r"^\s*\(; .* = (params|concs)$", ln),
            split(body, '\n')),
        '\n')
```

With:

```julia
    raw_lines = split(body, '\n')

    # 1. Identify the destructure header lines (kept verbatim) and
    #    section header lines (stripped — display-only).
    is_destructure(ln) = occursin(r"^\s*\(; .* = (params|concs)$", ln)
    is_section_header(ln) = occursin(r"^# .+ constraints:$", ln)

    # 2. Identify single-symbol equality lines (any section).
    is_single_eq(ln) = occursin(
        Regex("^\\s*\\w+\\s*=\\s*\\w+\\s*" *
              replace(ANNOTATION_SUBSTITUTED, "(" => "\\(", ")" => "\\)") *
              "\$"),
        ln)

    destructure_lines = String[ln for ln in raw_lines if is_destructure(ln)]
    single_eq_lines = sort!(String[ln for ln in raw_lines if is_single_eq(ln)])
    other_lines = String[
        ln for ln in raw_lines
        if !is_destructure(ln) && !is_section_header(ln) && !is_single_eq(ln)
    ]

    body = join(vcat(destructure_lines, single_eq_lines, other_lines), '\n')
    # Then strip destructure header lines for the alpha-rename phase
    body = join(
        filter(ln -> !is_destructure(ln), split(body, '\n')),
        '\n')
```

(Note: re-import `ANNOTATION_SUBSTITUTED` with `using EnzymeRates: ANNOTATION_SUBSTITUTED` at file top if needed; otherwise reference as `EnzymeRates.ANNOTATION_SUBSTITUTED`.)

### Task 3.7: Update existing fixtures with new section headers and annotations

**Files:**
- Modify: `test/mechanism_definitions_for_test_enzyme_derivation.jl`
- Modify: `test/test_rate_eq_derivation.jl:1142–1166`

- [ ] **Step 1: Identify fixtures with constraint lines.**

```bash
grep -n "k.r = \|K. = \|K.. = " test/mechanism_definitions_for_test_enzyme_derivation.jl test/test_rate_eq_derivation.jl
```

Each constraint in the fixtures needs the appropriate section header above it (User defined / Wegscheider / Haldane) and the `(substituted into v)` annotation if it's a single-symbol equality.

- [ ] **Step 2: For each fixture, structurally rewrite the expected string.**

Apply the new format. Example of a previous fixture string:

```julia
expected = """(; K1, k3f, Keq, E_total) = params
(; S, P) = concs
k3r = (1 / Keq) * K1 * (1 / K2) * k3f
v = E_total * (...)"""
```

Becomes:

```julia
expected = """(; K1, k3f, Keq, E_total) = params
(; S, P) = concs
# Haldane constraints:
k3r = (1 / Keq) * K1 * (1 / K2) * k3f
v = E_total * (...)"""
```

If there are user-defined merges or single-symbol Wegscheider ties, add a `# User defined constraints:` or `# Wegscheider constraints:` section with `(substituted into v)` annotations as appropriate.

The byte-identical allosteric fixture in `test/test_rate_eq_derivation.jl:1160–1164` also needs updating — possibly with a `# T-state constraints:` section if the implementation places T-state Haldane assignments there.

- [ ] **Step 3: Run existing fixture tests to confirm they FAIL under current code.**

```julia
include("test/test_rate_eq_derivation.jl")
```

Expected: byte-exact fixture tests fail. (Source-A and Source-B are still green; only the section-label additions are new failures.)

### Task 3.8: Run dedup tests; confirm Source C and section labels turn green

- [ ] **Step 1: Run.**

```julia
include("test/test_eq_hash_dedup.jl")
```

Expected:
- Sanity: green.
- Source A: green.
- Source B: green.
- Source C: **GREEN** ← this commit's target.
- Section labels: **GREEN** ← this commit's target.

If Source C is still red after Tasks 3.1–3.6: the canonicalizer's block-normalization regex isn't catching all single-symbol equalities. Diff a Source-C pair's canonical strings.

### Task 3.9: Create CSV-replay test

**Files:**
- Create: `test/test_dedup_csv_replay.jl`
- Modify: `test/runtests.jl`

- [ ] **Step 1: Write the test file.**

```julia
# ABOUTME: Replays params_estimate_{5,6,7,8}.csv through the new canonicalizer
# ABOUTME: and asserts eq_hash count matches distinct-loss count.

using CSV, DataFrames, Test
using EnzymeRates
using EnzymeRates: _canonical_rate_eq_hash

@testset "CSV dedup replay" begin
    for n in (5, 6, 7, 8)
        csv_path = joinpath(@__DIR__, "..", "dedup_investigation",
                            "params_estimate_$n.csv")
        if !isfile(csv_path)
            @info "skipping CSV replay for n=$n: $(csv_path) missing"
            continue
        end
        df = CSV.read(csv_path, DataFrame)

        new_hashes = map(eachrow(df)) do row
            m = eval(Meta.parse(row.mechanism_type))()
            _canonical_rate_eq_hash(m).short
        end
        df.new_hash = new_hashes

        # 1. Within-loss-group consistency.
        @testset "n=$n within-loss-group consistency" begin
            for g in groupby(df, :loss)
                @test length(unique(g.new_hash)) == 1
            end
        end

        # 2. Count match.
        @testset "n=$n hash count == loss count" begin
            n_loss = length(unique(round.(df.loss; sigdigits=10)))
            n_hash = length(unique(df.new_hash))
            @test n_loss == n_hash
        end
    end
end
```

- [ ] **Step 2: Add include to `test/runtests.jl`.**

Right after `include("test_eq_hash_dedup.jl")`:

```julia
    include("test_dedup_csv_replay.jl")
```

- [ ] **Step 3: Run the test in isolation.**

```julia
include("test/test_dedup_csv_replay.jl")
```

Expected:
- For each `n` in 5..8: within-loss-group consistency green; hash count equals loss count green.
- Predicted counts:
  - n=5: 11 distinct loss, 11 distinct new hash.
  - n=6: 24, 24.
  - n=7: 156, 156.
  - n=8: 467, 467.

If the count match fails for any n, dump the canonical strings for two rows with same loss but different new_hash and diff to find the residual quirk. Do not commit until classified.

### Task 3.10: Run full test suite

- [ ] **Step 1: Final full run.**

```bash
julia --project -e 'using Pkg; Pkg.test()'
```

Expected: ALL tests green. Aqua + JET checks still pass (no new stale deps from helper deletions; no new type-instability from the kernel extraction).

### Task 3.11: Commit 3

- [ ] **Step 1: Stage and commit.**

```bash
git add src/rate_eq_derivation.jl src/thermodynamic_constr_for_rate_eq_derivation.jl \
        src/identify_rate_equation.jl \
        test/mechanism_definitions_for_test_enzyme_derivation.jl \
        test/test_rate_eq_derivation.jl test/test_dedup_csv_replay.jl \
        test/runtests.jl
git commit -m "$(cat <<'EOF'
Source C: Wegscheider absorption + section labels + hash normalization

Source C fix: single-symbol Wegscheider RE ties are absorbed into the
kinetic-group rename map at polynomial level so v uses the
representative symbol only. _build_kinetic_rename_map gains Pass 2 that
calls _dependent_param_exprs_kernel with the user-only rename to
discover absorbable ties; the kernel is the existing Gaussian-elimination
logic extracted so it can be called twice with different rename maps.

rate_equation_string emits constraints under three section headers
("User defined constraints:", "Wegscheider constraints:",
"Haldane constraints:") with an explicit "(substituted into v)"
annotation on every single-symbol equality whose LHS got folded into
the rename map.

The eq_hash canonicalizer in identify_rate_equation strips section
headers (display-only) and normalizes the single-symbol equality lines
into a sorted block at a fixed canonical position, so Source C
duplicates collapse to one hash regardless of which section emitted
them.

Source-C and section-label dedup tests turn green. CSV-replay test
asserts eq_hash count equals distinct-loss count for n=5..8.

Co-Authored-By: Claude Opus 4.7 (1M context) <noreply@anthropic.com>
EOF
)"
```

---

## Final validation

### Task F.1: All tests green; final manual smoke check

- [ ] **Step 1: Full Pkg.test() one more time.**

```bash
julia --project -e 'using Pkg; Pkg.test()'
```

- [ ] **Step 2: Smoke check the CV results CSV against the new code.**

```julia
using CSV, DataFrames, EnzymeRates
using EnzymeRates: _canonical_rate_eq_hash

df = CSV.read("dedup_investigation/cv_results.csv", DataFrame)
df.new_hash = map(eachrow(df)) do row
    m = eval(Meta.parse(row.mechanism_type))()
    _canonical_rate_eq_hash(m).short
end

# Show the LDH n=7 cluster
cluster = df[in.(df.eq_hash, Ref(("831e36af", "9c7141ac", "89f33d51", "b362dd75"))), :]
println(unique(cluster.new_hash))
```

Expected: a single hash value across all four rows of the cluster.

- [ ] **Step 3: Verify the line count reduction.**

```bash
git diff --stat main src/ test/
```

Expected: net deletion of ~640 body lines (some + on tests for the dedup file, large − on src and the deleted test_sym_poly.jl).

### Task F.2: Verify CLAUDE.md is up to date

- [ ] **Step 1: Read the CLAUDE.md sections that mention rate-derivation files.**

Confirm the bullet for `src/rate_eq_derivation.jl` no longer mentions removed functions.

If any further references slipped through (e.g., to `_factor_poly` or the FactoredPoly types), update.

---

## Self-review notes

- **Spec coverage**: every requirement in `2026-05-09-rate-eq-dedup-and-simplification-design.md` has at least one task. Identifiability removal → Tasks 0.1–0.7. Kd-by-construction → Tasks 1.3–1.8. Always-expanded + helper deletion → Tasks 2.2–2.13. Section labels + Wegscheider absorption + canonicalizer → Tasks 3.1–3.6.
- **Doc updates**: CLAUDE.md (Task 0.7) and stale Phase G comments (Task 0.6).
- **TDD discipline**: failing tests come first (Tasks 0.8 dedup scaffold, 1.1+1.2 Source-B fixtures, 2.1 factored fixtures, 3.7 section fixtures). Each commit's success criterion is a test transitioning red→green.
- **Verification**: every commit ends with a Pkg.test() run; the final commit also runs CSV-replay against the saved 5/6/7/8 distributions.
- **Validation that the design's projected reductions hold**: Task F.1 Step 3.
