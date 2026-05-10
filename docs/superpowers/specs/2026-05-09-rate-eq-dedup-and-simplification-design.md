# Rate-equation deduplication and derivation-code simplification

## Background

`identify_rate_equation` produces many mechanisms whose rate equations are
algebraically identical but hash to different `eq_hash` values, so the LDH run
saved in `dedup_investigation/` ends up with `eq_hash` counts that are 3-7×
larger than the count of distinct loss values:

| n | distinct `eq_hash` | distinct loss | cross-hash dup % |
|---|------------------|----------------|--------------------|
| 5 |  13 |  11 | 15% |
| 6 |  62 |  24 | 61% |
| 7 | 493 | 156 | 68% |
| 8 |1390 | 467 | 66% |

`dedup_investigation/investigation_eq_hash_duplication.md` traces the
duplication to three independent textual sources:

- **Source A — factoring variants.** `_factor_poly` and
  `_try_algebraic_factor_sigma` produce different factorings of the same
  polynomial across mechanisms.
- **Source B — Ka↔Kd inversion artifact.** Polynomials are stored internally
  with K as the association constant Ka and inverted at print time via
  `_apply_kd_inversion` and `inverted_params` threading; the double-wrap
  produces `K = 1/(1/X)` lines for binding-K↔binding-K Haldane closures.
- **Source C — split-with-tie ≡ pre-merged.** When a mechanism's enumerator
  output has a kinetic-group split that Wegscheider then ties back together
  with a single-symbol equality (`K8 = K4`), the split form has the same
  rate equation as the pre-merged form but a different mechanism type.

Stage 2 fitting in `_beam_search` is already deduped by `eq_hash`
(`identify_rate_equation.jl:530-561` runs one fit per representative hash and
projects results to all members), so the cost of cross-hash duplication is
**redundant compile work plus inflated CSV row counts**, not redundant fits.

## Objectives

1. **Collapse all three sources at the `eq_hash` level** so distinct hashes
   match distinct loss values on the LDH validation CSVs.
2. **Reduce code volume** in the rate-equation derivation pipeline. Code
   simplification is a co-equal objective with the dedup fix; comments and
   docstrings don't count toward the line-reduction goal. Estimated net
   reduction is ~770 body lines across `rate_eq_derivation.jl`,
   `sym_poly_for_rate_eq_derivation.jl`,
   `thermodynamic_constr_for_rate_eq_derivation.jl`,
   `test/test_sym_poly.jl`,
   `test/mechanism_definitions_for_test_enzyme_derivation.jl`, and
   `test/test_rate_eq_derivation.jl`, plus ~15 lines added in
   `identify_rate_equation.jl` for the canonicalizer.
3. **Improve user-facing rate-equation transparency** by labeling derived
   constraints under three section headers: `# User defined constraints:`,
   `# Wegscheider constraints:`, `# Haldane constraints:`.
4. **Remove the parameter-identifiability code.** The current
   `structural_identifiability_deficit` implementation is not correct for
   the use cases that matter and will be redesigned later. Remove it now
   so that future work isn't constrained by the existing shape, and so
   the polynomial-expansion helpers it depends on can be removed too.

## Methodology

**TDD is required throughout this work, not optional.** For each of the
four commits, the order is:

1. Write or update tests that capture the post-change behavior (new dedup
   tests; updated rate-equation-string fixtures with the expanded /
   re-formatted expected strings).
2. Run the test suite and confirm the new and updated tests fail under
   current code.
3. Implement the change.
4. Run the test suite and confirm the targeted tests turn green and no
   existing tests regress.

This applies symmetrically to (a) the new dedup tests in
`test_eq_hash_dedup.jl`, (b) the CSV-replay test, and (c) the existing
rate-equation-string fixtures that need updating because the printed form
changes.

## Non-objectives

- Changing the optimizer pipeline, beam search, or model-selection logic.
- Within-hash multiplicity (multiple `mechanism_type` structs sharing one
  `eq_hash`). By design, multiple specs per hash get carried forward to the
  next expansion stage so split structures remain available for further
  refinement.
- The "extra parameter pinned at 1.0" pattern. A consequence of within-hash
  multiplicity at the current `n_params` level; resolves at higher `n_params`
  when the parameter actually enters the rate equation.
- Spec-level absorption / `_dedup_key` change. Source C duplicates remain
  distinct specs in the cache and pay the (cheap) compile cost; the fit cost
  is already deduped via `eq_hash`.

## Architecture

Four code changes, organized as four commits in one PR:

| # | Source | Change |
|---|--------|--------|
| 0 | — | Remove parameter-identifiability code (`structural_identifiability_deficit`, `_count_allosteric_rate_monomials`, the `expected_is_identifiable` test fixtures, the export, and the file-header references). Independent of the dedup fix; lands first because it unblocks the larger expansion-helper deletion in Commit 2. |
| 1 | **B** | Build polynomials directly with K_d. Eliminate the Ka↔Kd inversion layer. |
| 2 | **A** | Always emit the expanded polynomial. Drop algebraic factoring (allosteric retains the structurally-required MWC outer factoring; inner conformation polynomials expand). With identifiability code already gone, the polynomial-expansion helpers (`_expand_factored_*`, `_expand_to_poly`, `_estimate_expanded_term_count`) become unreachable and are deleted, along with the `FactoredPoly` / `FactoredSigma` / `DenomTerm` types and the entire `test/test_sym_poly.jl` file (299 lines, exclusively tests these types). |
| 3 | **C** | Polynomial-level absorption of single-symbol Wegscheider RE ties into the kinetic-group rename map. Section-labeled output. eq_hash canonicalizer normalizes single-symbol equality lines across sections. |

`EnzymeMechanism` and `AllostericEnzymeMechanism` types and their type
parameters are **unchanged**. Spec-level `_canonicalize!` and `_dedup_key`
in `mechanism_enumeration.jl` are **unchanged**.

## Step 0 — Remove parameter-identifiability code

`structural_identifiability_deficit` is being removed in full. The current
implementation over-counts identifiable degrees of freedom for factored
forms and is acknowledged in test code as not biophysically meaningful at
its current shape (`test/test_rate_eq_derivation.jl:479–481`). It will be
redesigned later; removing it now unblocks the polynomial-helper deletions
in Step 2 and avoids constraining the future redesign.

### What goes away

**Source code (`src/`):**

- `EnzymeRates.jl:16` — drop the `export structural_identifiability_deficit`
  line.
- `rate_eq_derivation.jl:780–800` — `structural_identifiability_deficit`
  for `EnzymeMechanism` (~20 lines including docstring).
- `rate_eq_derivation.jl:1672–1779` (file end) —
  `structural_identifiability_deficit` for `AllostericEnzymeMechanism` and
  `_count_allosteric_rate_monomials` (~110 lines combined; the latter
  function alone runs 1690–1778).
- `identify_rate_equation.jl:246` — reword the comment explaining the
  defensive `_proj` lookup to drop the "structurally-unidentifiable ghost
  param" framing. The defensive behavior itself (fallback when a fitted
  key isn't in the body) stays — it's still needed for `:NonequalRT`
  zeroed paths.

**Tests:**

- `test/test_rate_eq_derivation.jl:475–484` — `test_identifiability`
  helper function. Remove.
- `test/test_rate_eq_derivation.jl:751` — call-site to
  `test_identifiability(spec)`. Remove.
- `test/test_rate_eq_derivation.jl:2` — file-header comment lists
  "identifiability" as one of the validated aspects. Drop the word.
- `test/test_identify_rate_equation.jl:54` — comment reads "5
  identifiable params + 1 ghost". Reword to drop "identifiable".
- `test/mechanism_definitions_for_test_enzyme_derivation.jl:32` —
  `expected_is_identifiable::Bool` field on `MechanismTestSpec`. Remove
  the field.
- `test/mechanism_definitions_for_test_enzyme_derivation.jl` — 37
  occurrences of `expected_is_identifiable=true` /
  `expected_is_identifiable=false` across the fixture definitions.
  Remove each line. (`grep -c` reports 38 total when counting the
  field declaration as well.)

**CLAUDE.md:**

- `.claude/CLAUDE.md:291` — bullet describing
  `src/rate_eq_derivation.jl` mentions "identifiability checks" and
  `_count_allosteric_rate_monomials` and `structural_identifiability_deficit`.
  Update bullet to drop those references.

**Stale "Phase G" planning-phase doc references:**

- `src/identify_rate_equation.jl:147` — docstring for
  `_canonicalize_rate_eq_with_map` says "(Phase G.0)". Drop the parenthetical.
- `src/rate_eq_derivation.jl:22` — docstring for `parameters` mentions
  "Phase G's rate-equation canonicalizer". Reword to drop the phase label.

These reference an old planning phase that no longer exists; the
descriptions stay relevant once the phase tag is dropped.

### Why this lands first

The expansion helpers `_expand_factored_sigma`, `_expand_to_poly`,
`_expand_factored_poly`, and `_estimate_expanded_term_count` have callers
in three places after the rest of the cleanup:

1. `_kcat_components` and `_kcat_forward` (allosteric branch) — but Step 2
   makes `_raw_symbolic_rate_polys` return flat POLYs, eliminating these
   call sites.
2. `_allosteric_num_den_exprs` — same: Step 2 hands it flat POLYs.
3. `structural_identifiability_deficit` and `_count_allosteric_rate_monomials`
   — only removable if identifiability is removed.

Without Step 0, the helpers stay alive solely to support the broken
identifiability path, and Step 2 can't delete them. Putting Step 0 first
means Step 2 can eliminate the entire factored-poly type family cleanly.

## Step 1 — Kd-by-construction

### Polynomial construction

In `_compute_alpha` (`rate_eq_derivation.jl:149`), for each binding RE step
`E + S ⇌ ES`, the alpha relation flips so K appears in `alpha_den` (Kd
convention) instead of `alpha_num` (Ka convention):

```julia
# Before (Ka):
alpha_num[j_f] = poly_mul(poly_mul(alpha_num[cur], K), mp(m_l))
alpha_den[j_f] = poly_mul(alpha_den[cur], mp(m_r))

# After (Kd):
alpha_num[j_f] = poly_mul(alpha_num[cur], mp(m_l))
alpha_den[j_f] = poly_mul(poly_mul(alpha_den[cur], K), mp(m_r))
```

Non-binding RE steps (pure isomerization, both `m_l` and `m_r` empty) keep
K in `alpha_num` — Ka convention is correct for iso steps and unchanged.

### Haldane / Wegscheider elimination

`_dependent_param_exprs` in
`thermodynamic_constr_for_rate_eq_derivation.jl:198` builds the A-matrix
with a sign-flip on binding-K columns so the dep_exprs come out directly in
Kd form:

```julia
# Identify binding K's once (steps with metabolite on LHS).
binding_set = Set{Symbol}(...)

for i in 1:nc, j in 1:nsteps
    C[i, j] == 0 && continue
    if eq_steps[j]
        sym = Symbol("K$j")
        sym = get(rename, sym, sym)
        sign = sym in binding_set ? -1 : 1
        A[i, sym_col[sym]] += sign * C[i, j]
    else
        # SS step: kf gets +C[i,j], kr gets -C[i,j], no convention flip.
    end
end
```

Constraint lines render as plain `K8 = K2 * K3 / Keq` instead of double-wrapped
`K8 = 1/(1/K2 * 1/K3 / Keq)`.

### Code that goes away

- `_apply_kd_inversion` (`thermodynamic_constr_for_rate_eq_derivation.jl:295`)
  and all its callers.
- `_binding_K_symbols` for non-allosteric usage (still used internally by
  the allosteric MWC outer factoring; kept and reused).
- `inverted_params` parameter threaded through `_poly_to_expr`,
  `_factored_sigma_to_expr`, `_factored_poly_to_expr`,
  `_denom_terms_to_expr`, `to_rate_expr`. Remove the parameter; `_poly_to_expr`
  emits negative exponents naturally.
- `inv_fn` callbacks at line 312 (`K → 1/K`) and line 348 (`K → inv($K)`) of
  the constraint module.
- `inv_set = Set(_binding_K_symbols(M))` construction and threading in
  `_raw_rate_expr_and_symbols` and `_rate_v_line`.

### `_kcat_forward` cleanup

Both `EnzymeMechanism` and `AllostericEnzymeMechanism` overloads
(`rate_eq_derivation.jl:890` and `:929`) currently substitute `K → inv($K)` to
compensate for the Ka-stored / Kd-fitted mismatch. Step 1 deletes these
substitutions; `_kcat_components` and `_kcat_forward` operate on already-Kd
polynomials.

## Step 2 — Always-expanded emission

Drop algebraic factoring of the printed numerator and denominator for the
non-allosteric path. The `@generated rate_equation` body and
`rate_equation_string` both emit fully-expanded polynomials via
`_poly_to_expr`. The allosteric MWC outer factoring
(`E_total * (N_R + L*N_T) / (Q_R^N + L*Q_T^N)`) is structurally required and
stays. Inner conformation polynomials (N_R, Q_R, N_T, Q_T) expand to flat
POLYs.

### Code that goes away

- `_factor_poly` (`rate_eq_derivation.jl:405`).
- `_try_algebraic_factor_sigma` (`rate_eq_derivation.jl:210`).
- `_try_poly_power` (`rate_eq_derivation.jl:343`).
- `_haldane_equality_substitutions` (`rate_eq_derivation.jl:380`) — handled
  uniformly by Step 3's Pass 2 absorption.
- `_factored_sigma_to_expr`, `_factored_poly_to_expr`, `_denom_terms_to_expr`
  in `sym_poly_for_rate_eq_derivation.jl`.
- `FactoredPoly`, `FactoredSigma`, `DenomTerm` types and
  `unfactored_denom_term`. The `_rename_symbols` overloads for these types.
- `_estimate_expanded_term_count` and the `check_benefit` flag.
- `to_rate_expr`'s `Union{POLY, FactoredSigma}` branch.
- `_expand_factored_sigma`, `_expand_to_poly`, `_expand_factored_poly`
  in `sym_poly_for_rate_eq_derivation.jl`. Confirmed unreachable after
  Step 0 removes their identifiability callers and Step 2 makes
  `_raw_symbolic_rate_polys` return flat POLYs.
- `_poly_power` (`sym_poly_for_rate_eq_derivation.jl:48`). All four
  callers are gone after Step 0 + Step 2:
  `_count_allosteric_rate_monomials` (Step 0), `_try_poly_power` and
  `_expand_factored_poly` (Step 2). Confirmed unreachable.
- `_try_poly_exact_div` (`sym_poly_for_rate_eq_derivation.jl:69–116`).
  Only caller is `_try_algebraic_factor_sigma` (Step 2 deletion). ~50
  lines of polynomial exact-division logic become unreachable.
- `test/test_sym_poly.jl` — entire file (299 lines). Tests
  `FactoredPoly` / `FactoredSigma` / `DenomTerm` construction and the
  `_expand_*` helpers exclusively. With those types and helpers gone,
  there's nothing left to test.

### Code that stays

- `MAX_RATE_EQUATION_TERMS = 5000` cap. Already enforced by `sym_det`.
- `_poly_to_expr` and `_expr_to_string`. Former simplified (no
  `inverted_params` parameter).
- `_zero_symbols_in_poly` — still used by allosteric T-state polynomial
  construction on flat POLYs.

### `_rate_v_line` shrinks to ~4 lines

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

(versus 14 lines today at `rate_eq_derivation.jl:750`).

### `_raw_symbolic_rate_polys` shape change

Returns `(num::POLY, denom::POLY)` instead of
`(num_fs::FactoredSigma, denom_terms::Vector{DenomTerm})`.

### Allosteric path

`_allosteric_num_den_exprs` (`rate_eq_derivation.jl:1526`): N_R, Q_R, N_T,
Q_T come directly as flat POLYs from `_raw_symbolic_rate_polys(CM)` and
are converted to Expr via `_poly_to_expr`. No `_expand_*` calls. The MWC
outer assembly (`make_num_term`, `make_den_term`, the `(N_R + L*N_T)`
sum) is untouched.

### `_kcat_components` and allosteric `_kcat_forward`

Both consume the now-flat POLY return from `_raw_symbolic_rate_polys`
directly. The current `num = _expand_factored_sigma(num_fs)` /
`den = _expand_to_poly(denom_terms)` lines collapse to using `num` and
`den` from the destructured return. Allosteric `_kcat_forward`
(`rate_eq_derivation.jl:929`) similarly drops the expansion calls — the
T-state POLY construction via `_zero_symbols_in_poly` and `_rename_symbols`
applies directly to the flat POLY.

## Step 3 — Polynomial-level Wegscheider absorption + section labels

### Polynomial-level rename absorption (no spec-level change)

Extend `_build_kinetic_rename_map` (`rate_eq_derivation.jl:90`) to a
two-pass design:

```julia
function _build_kinetic_rename_map(M::Type{<:EnzymeMechanism})
    rename = Dict{Symbol, Symbol}()
    m = M()
    eq = equilibrium_steps(m)

    # Pass 1: user-defined kinetic-group merges (current logic).
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

    # Pass 2: single-symbol Wegscheider RE ties from raw dep_exprs.
    dep_exprs_raw, _ = _dependent_param_exprs_kernel(M, rename)
    for (lhs, rhs) in dep_exprs_raw
        rhs isa Symbol || continue
        is_re_K(lhs) && is_re_K(rhs) || continue
        rename[lhs] = get(rename, rhs, rhs)
    end

    rename
end
```

`_dependent_param_exprs_kernel` is the Gaussian-elimination kernel extracted
from `_dependent_param_exprs`; takes a rename-map argument so it can be
called with Pass-1-only or full rename. `_dependent_param_exprs` itself
calls the kernel with the full rename and returns dep_exprs that exclude
single-symbol Wegscheider ties (those are absorbed via the rename).

`EnzymeMechanism` type info is unchanged. The polynomial in `v` uses the
representative symbol everywhere.

### Section labels in `rate_equation_string`

```
(; K1, k4f, k4r, Keq, E_total) = params
(; NADH, Pyruvate, Lactate, NAD) = concs
# User defined constraints:
K7 = K1  (substituted into v)
# Wegscheider constraints:
K4 = K1  (substituted into v)
# Haldane constraints:
k10r = (1 / Keq) * k4f * (1 / k4r) * (1 / K1) * k10f
v = E_total * (... K1 ...) / (... K1 ...)
```

- **`# User defined constraints:`** — equalities derived from kinetic-group
  structure in the `EnzymeMechanism` type. Single-symbol RE ties from
  user-defined groups; SS group ties (both `k_f` and `k_r`).
- **`# Wegscheider constraints:`** — equalities from raw dep_exprs (via
  `_dependent_param_exprs_kernel` with Pass-1-only rename) where the RHS
  doesn't reference `Keq`. Includes both single-symbol absorbed ties (still
  printed for transparency) and multi-symbol expressions.
- **`# Haldane constraints:`** — equalities from raw dep_exprs where the
  RHS references `Keq`.

The annotation `(substituted into v)` is appended to every single-symbol
equality whose LHS got folded into the rename map (any such tie, whether it
appears in the User defined or Wegscheider section). Multi-symbol entries
get no annotation; their LHS appears in `v` as a runtime substitution.

A single constant prevents annotation drift:

```julia
const ANNOTATION_SUBSTITUTED = "  (substituted into v)"
```

Empty sections are suppressed (no header printed).

### Runtime body

`_build_rate_body` (used by `@generated rate_equation`) uses the
post-absorption `_dependent_param_exprs` output: only multi-symbol Wegscheider
and Haldane closures get runtime assignment statements. Single-symbol
absorbed ties don't need runtime assignments — the polynomial-rename already
substituted them in the Expr tree.

### eq_hash canonicalization

`_canonicalize_rate_eq_with_map` (`identify_rate_equation.jl:151`) gains one
normalization step. After stripping destructure header lines (current
behavior), it also:

1. Strips any line matching `^# .* constraints:$` (section headers — display
   only).
2. Collects every line whose body matches a single-symbol equality
   (`^\w+\s*=\s*\w+\s*\(substituted into v\)$`).
3. Sorts the collected lines lexicographically.
4. Re-emits them as a single contiguous block at a fixed canonical position
   (e.g., immediately after the `concs` destructure line).

Multi-symbol equality lines and the `v` line pass through unchanged. After
normalization, mechanisms α and β with the same canonical content (same set
of single-symbol ties + same multi-symbol closures + same v polynomial)
produce byte-identical text → matching `eq_hash`.

## Validation

### TDD tests (write before implementing each step)

Add `test/test_eq_hash_dedup.jl` with one test per source. Each test loads
specific mechanisms from `dedup_investigation/cv_results.csv` by parsing
the `mechanism_type` column with `Meta.parse` + `eval`.

```julia
using CSV, DataFrames, Test, EnzymeRates
using EnzymeRates: _canonical_rate_eq_hash

const _CV = CSV.read(joinpath(@__DIR__, "..", "dedup_investigation",
                              "cv_results.csv"), DataFrame)
_mech(row_idx) = eval(Meta.parse(_CV[row_idx, :mechanism_type]))()

@testset "eq_hash dedup" begin
    # Sanity: row indices match the investigation's eq_hashes.
    @test _CV[22, :eq_hash] == "831e36af"
    @test _CV[27, :eq_hash] == "9c7141ac"
    @test _CV[31, :eq_hash] == "89f33d51"
    @test _CV[36, :eq_hash] == "b362dd75"

    # Joint cluster collapse — does not isolate which source caused
    # the collapse. The CSV-replay test (Task 3.9) is the orthogonal
    # source-validation across all sources jointly.
    @testset "LDH n=7 cluster collapses to one hash" begin
        for j in (27, 31, 36)
            @test _canonical_rate_eq_hash(_mech(22)) ==
                  _canonical_rate_eq_hash(_mech(j))
        end
    end

    @testset "section labels render correctly" begin
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
```

Each test starts red under current code (current 4-mechanism cluster has 4
distinct hashes). Green-turning order matches commit order: Source-B tests
green after Commit 1; Source-A tests green after Commit 2; Source-C and
section-label tests green after Commit 3.

If any test stays red after its commit, that's the falsification signal — do
not merge until classified.

### CSV-replay test

Add `test/test_dedup_csv_replay.jl`:

```julia
@testset "CSV dedup replay" for n in (5, 6, 7, 8)
    csv_path = joinpath(@__DIR__, "..", "dedup_investigation",
                        "params_estimate_$n.csv")
    isfile(csv_path) || (@info "skipping n=$n: CSV missing"; continue)
    df = CSV.read(csv_path, DataFrame)

    new_hashes = map(eachrow(df)) do row
        m = eval(Meta.parse(row.mechanism_type))()
        _canonical_rate_eq_hash(m).short
    end
    df.new_hash = new_hashes

    # Within-loss-group consistency.
    for g in groupby(df, :loss)
        @test length(unique(g.new_hash)) == 1
    end

    # Count match.
    n_loss = length(unique(round.(df.loss; sigdigits=10)))
    n_hash = length(unique(df.new_hash))
    @test n_loss == n_hash
end
```

Expected pass after Commit 3:

| n | distinct loss | predicted `eq_hash` count |
|---|---|---|
| 5 | 11 | 11 |
| 6 | 24 | 24 |
| 7 | 156 | 156 |
| 8 | 467 | 467 |

The test runs in seconds (compile rate equations + hash, no fitting).

**Do not rerun `identify_rate_equation` end-to-end on LDH.** A full run takes
several hours and the saved CSVs already cover the validation surface.

### Existing rate-equation-string tests will break

Many tests in `test/test_rate_eq_derivation.jl` and
`test/mechanism_definitions_for_test_enzyme_derivation.jl` compare
rate-equation strings byte-for-byte:

- 20 `expected_factored_num` / `expected_factored_denom` fixtures in
  `mechanism_definitions_for_test_enzyme_derivation.jl` (lines 962–1707+).
  All break under Step 2 — the factored expectations need to be replaced
  with their expanded equivalents.
- One byte-identical allosteric fixture at
  `test/test_rate_eq_derivation.jl:1142–1166`. Breaks under all three steps:
  Step 1 removes `1/(1/K)` artifacts, Step 2 expands inner conformation
  polynomials, Step 3 adds section headers.

**TDD discipline for fixture updates.** Update expected strings *before*
implementing each step, watch them fail, then implement and watch them pass.

The mechanical transformations are:

- **Step 1**: replace every `1 / (1 / X)` substring with `X`. Rewrite Kd
  display: `K * met` → `met / K` for binding-K monomials. Mechanical text
  substitution; doable by hand on the 20 fixtures.
- **Step 2**: expand each factored denominator product to a sum of
  monomials. Tractable for bisubstrate (factored form has 2-4 factors of
  size 2-3). For larger mechanisms with many regulators, expansion may
  exceed manual reliability — for those fixtures, lock in numerical
  equivalence (`test_rate_equation_string` already asserts this) and
  drop the byte-exact comparison, OR derive the new expected by running
  the implementation and committing the output (acceptable here because
  numerical equivalence is the primary correctness guarantee).
- **Step 3**: add the three section-header lines + `(substituted into v)`
  annotations. Mechanical.

Mark in the spec which fixtures get hand-updated vs. which fall back to
numerical-equivalence-only. The byte-identical allosteric fixture at line
1142 is large enough to warrant the fallback approach; the 20 small-mechanism
fixtures should be hand-updated.

### Existing test suite

`Pkg.test()` must remain all-green: Aqua, JET, mechanism enumeration,
analytical kcat formulas, the existing `test_rate_equation_string` numerical
equivalence check (will pass — math is unchanged).

## Code-shrink summary

| File | Lines removed | Lines added | Net |
|---|---|---|---|
| `src/sym_poly_for_rate_eq_derivation.jl` | ~210 (FactoredPoly/Sigma/DenomTerm types and rename overloads, _factored_*_to_expr, _expand_factored_*, _expand_to_poly, _estimate_expanded_term_count, _poly_power, _try_poly_exact_div, inverted_params parameter) | ~5 | ≈ −205 |
| `src/rate_eq_derivation.jl` | ~200 (_factor_poly, _try_algebraic_factor_sigma, _try_poly_power, _haldane_equality_substitutions, structural_identifiability_deficit ×2, _count_allosteric_rate_monomials, _binding_K_symbols non-allosteric usage, inv_set threading, kcat inv($K) substitutions) | ~25 (Pass 2 of _build_kinetic_rename_map, _dependent_param_exprs_kernel extraction, sectioned rate_equation_string, ANNOTATION_SUBSTITUTED constant) | ≈ −175 |
| `src/thermodynamic_constr_for_rate_eq_derivation.jl` | ~25 (_apply_kd_inversion + callers, inv_fn plumbing) | ~5 (binding-K sign-flip in A-matrix) | ≈ −20 |
| `src/EnzymeRates.jl` | 1 (export of structural_identifiability_deficit) | 0 | ≈ −1 |
| `src/identify_rate_equation.jl` | 0 | ~15 (canonicalizer block normalization + section-header strip) | ≈ +15 |
| `test/test_sym_poly.jl` | 299 (entire file) | 0 | ≈ −299 |
| `test/test_rate_eq_derivation.jl` | ~12 (test_identifiability helper + caller) | ~30 (eq_hash dedup tests) | ≈ +18 |
| `test/mechanism_definitions_for_test_enzyme_derivation.jl` | ~45 (expected_is_identifiable field + 42 occurrences) | 0 | ≈ −45 |
| `test/test_eq_hash_dedup.jl` (new) | 0 | ~50 | ≈ +50 |
| `test/test_dedup_csv_replay.jl` (new) | 0 | ~25 | ≈ +25 |
| **Total** | **~922** | **~155** | **≈ −767** |

(Updated to include items added during 2nd-round review: inlined
`_kcat_components` and `_rate_v_line`, deleted `_constraint_expr_strings`,
deleted both `_binding_K_symbols` overloads, inlined `_classify_cycle`
into `_thermodynamic_constraints`, plus the existing 2-arg-overload
inlines.)

Comments and docstrings don't count toward the reduction. The number above
is body-line delta only, measured against current `git HEAD`. Per-file
estimates have ±20% uncertainty until the implementation lands; the
order-of-magnitude reduction is the binding commitment.

## Risks

1. **`MAX_RATE_EQUATION_TERMS = 5000` ceiling on always-expanded.**
   Pre-merge sanity check: pick the largest mechanisms in
   `params_estimate_8.csv` and confirm each compiles under the new path.

2. **Allosteric inner-polynomial expansion.** No allosteric mechanisms in
   the LDH CSVs, so the CSV-replay test won't exercise this path. Add one
   TDD test that an allosteric mechanism's `rate_equation_string` renders
   without errors and `_kcat_forward` returns a non-zero finite value.

3. **`_dependent_param_exprs_kernel` extraction.** Pass 2 of
   `_build_kinetic_rename_map` calls the elimination kernel with no
   rename-absorption; full `_dependent_param_exprs` calls it with the
   full rename. Refactor the kernel cleanly to avoid recursion. Risk is
   local.

4. **`(substituted into v)` annotation drift.** Single constant
   `ANNOTATION_SUBSTITUTED` used by both emitter sites prevents
   whitespace drift. The TDD section-label test verifies.

5. **`mechanism_type` column from CSVs needs `eval`-able strings.** The
   current `string(typeof(m))` output is parseable Julia syntax. Verified
   by the CSV-replay test itself — if any row's `mechanism_type` doesn't
   round-trip, that's caught immediately.

## PR plan

Four commits in dependency order, single PR:

1. **Commit 0 — Remove identifiability code.** Drop `structural_identifiability_deficit`
   (both methods), `_count_allosteric_rate_monomials`, the export, the
   `expected_is_identifiable` field and 42 fixture occurrences,
   `test_identifiability` helper and caller, the file-header reference,
   the `:54` ghost-param comment, and the CLAUDE.md bullet entry. The
   `_proj` defensive-lookup comment in `identify_rate_equation.jl` gets
   reworded to drop the "structurally-unidentifiable" framing. Pure
   deletion; tests must remain green after this commit (existing test
   suite minus the removed helper).

2. **Commit 1 — Kd-by-construction.** `_compute_alpha` flip; sign-flip in
   A-matrix; drop `_apply_kd_inversion`, `inverted_params` threading,
   `inv_set` plumbing, `inv` substitutions in `_kcat_forward`. Update
   Source-B fixture expectations first; Source-B TDD tests turn green.

3. **Commit 2 — Always-expanded + helper deletion.** Drop `_factor_poly`
   family; `_factored_*_to_expr`; `FactoredPoly`/`FactoredSigma`/`DenomTerm`;
   `_estimate_expanded_term_count`; `_haldane_equality_substitutions`;
   `_expand_factored_*` and `_expand_to_poly`; `_poly_power` and
   `_try_poly_exact_div` (now unreachable after the above); the entire
   `test/test_sym_poly.jl` file. Simplify `_rate_v_line`,
   `_raw_symbolic_rate_polys`. Update allosteric inner-polynomial emission
   to consume flat POLYs directly. Update the 20 factored-form fixtures
   (hand expansion) and the byte-identical allosteric fixture (regenerate
   from the new code after numerical-equivalence verification). Source-A
   TDD tests turn green.

4. **Commit 3 — Section labels + Wegscheider absorption + canonicalizer
   normalization.** Extract `_dependent_param_exprs_kernel`; Pass 2 of
   `_build_kinetic_rename_map`; sectioned `rate_equation_string`;
   `(substituted into v)` annotation; canonicalizer single-symbol-equality
   block normalization. Source-C TDD tests + section-label tests + CSV
   replay test pass.
