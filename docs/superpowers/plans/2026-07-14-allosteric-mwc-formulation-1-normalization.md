# Allosteric MWC formulation-1 normalization — Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Fix the allosteric MWC L-term leak by normalizing each conformation's King–Altman partition to its own free-enzyme segment weight `D[g_free]` before combining, so the derived rate equation matches the `n=1` two-conformation mass-action ground truth.

**Architecture:** Resurrect the reverted `D[g_free]` surfacing (`_raw_symbolic_rate_polys`/`_state_rate_polys` return a third value `d_free`), then apply a per-conformation rendering rule in `_allosteric_num_den_exprs` and `_kcat_forward`: raw combine when `D_A == D_I`, divide `Q/D` when `D` is a metabolite-free monomial, cross-weight by the other state's `D^n` otherwise. No discriminator on topology.

**Tech Stack:** Julia; the package's `POLY` Laurent-polynomial type (`src/sym_poly_for_rate_eq_derivation.jl`); `@generated` rate-equation derivation (`src/rate_eq_derivation.jl`); `Test` stdlib; the `n=1` ground-truth harness (`test/allosteric_ground_truth.jl`).

## Global Constraints

- 92-character line length, 4-space indentation. Match surrounding style.
- `rate_equation` MUST stay allocation-free and under 120 ns/call for every mechanism in `MECHANISM_TEST_SPECS` (`test_rate_equation_performance`, `test/test_rate_eq_derivation.jl`).
- All `Parameter → Symbol` rendering flows through the `name(p, m)` chokepoint. No stray `Symbol("K…")`/`Symbol("k…")` literals (guarded by `test/test_types.jl:1577`).
- Canonical Step Form is load-bearing; do not reorder steps/groups.
- Every ground truth self-validates (`L=0` → active-only rate; identical conformations → base rate, `L`-independent) before it gates the derivation.
- Every file starts with two `# ABOUTME:` lines. Do not add temporal/historical comments.
- Run tests before every commit. Commit frequently.

## Reference commits (adapt, do not cherry-pick blindly)

- `0e4c556` "Surface free-enzyme segment weight D[g_free] per allosteric state" — the exact D-surfacing threading. Task 1 reproduces it.
- `eae1c6c` "Cross-weight MWC free-enzyme normalization" — the cross-weight mechanics. Task 3/4 adapt it, **minus** its unconditional application, **minus** its ping-pong fail-loud guard, **plus** the `D_A==D_I` skip and the monomial-divide branch.

Both are on the current branch's history (reverted by `13a0020`). Read them with `git show <sha> -- src/rate_eq_derivation.jl`.

## File structure

- `src/rate_eq_derivation.jl` — the derivation. All source changes land here.
  - `_raw_symbolic_rate_polys` (2 methods), `_state_rate_polys` — return `d_free` (Task 1).
  - `_invert_monomial`, `_is_metabolite_free_monomial`, `_mwc_cross_weight` — helpers (Task 2).
  - `_allosteric_num_den_exprs` — the three-way rendering (Task 3).
  - `_kcat_forward` (allosteric method) — mirror the rendering at the poly level (Task 4).
- `test/allosteric_ground_truth.jl` — the acceptance gates: flip three `@test_broken`, switch the `:NonequalAI` reference, add metabolite-at-zero and ping-pong gates (Tasks 3, 5).
- `test/reference/allosteric_golden_reference.txt` — regenerate (Task 6).
- `docs/src/deriving/mwc_allostery.md` — rewrite the model description and formula (Task 7).

---

### Task 1: Surface `D[g_free]` per conformation

**Files:**
- Modify: `src/rate_eq_derivation.jl` — `_raw_symbolic_rate_polys(mech, …)` (~line 331), its 1-arg method (~line 389), `_raw_rate_expr_and_symbols` (~626), `_rate_v_line` (~742), `_kcat_forward` uni method (~889), `_kcat_forward` allo method (~942), `_dependent_param_exprs` (~1474), `_state_rate_polys` (~1248), `_allosteric_num_den_exprs` (~1602/1618).
- Test: `test/test_rate_eq_derivation.jl`

**Interfaces:**
- Produces: `_raw_symbolic_rate_polys(...) → (num::POLY, den::POLY, d_free::POLY)` and `_state_rate_polys(am, state) → (num::POLY, den::POLY, d_free::POLY)`. `d_free` is `D[form_to_group[i_free]]` renamed, where `i_free` is the form with empty `bound` and empty `residual`; `poly_one()` for a single-segment graph.

- [ ] **Step 1: Write the failing test** — append to `test/test_rate_eq_derivation.jl` (near the other `_state_rate_polys` tests):

```julia
@testset "D[g_free] surfaced per allosteric state" begin
    onlyA = @allosteric_mechanism begin
        substrates: S ; products: P ; catalytic_multiplicity: 1
        catalytic_steps: begin
            E + S ⇌ E(S)     :: OnlyA
            E(S) <--> E(P)   :: EqualAI
            E + P ⇌ E(P)     :: EqualAI
        end
    end
    am = EnzymeRates.AllostericMechanism(onlyA)
    _, _, dA = EnzymeRates._state_rate_polys(am, :A)
    _, _, dI = EnzymeRates._state_rate_polys(am, :I)
    @test dA == EnzymeRates.poly_one()                       # single active segment
    @test dI == EnzymeRates.POLY(EnzymeRates._mono(:k_ES_to_EP => 1) => 1)
end
```

- [ ] **Step 2: Run it, expect failure** — `julia --project -e 'using Pkg; Pkg.test()'` is slow; instead run just this file's relevant path. Expected: a `MethodError`/tuple-destructuring failure because `_state_rate_polys` returns a 2-tuple.

Run: `julia --project -e 'using EnzymeRates; const ER=EnzymeRates; m=@allosteric_mechanism begin substrates: S; products: P; catalytic_multiplicity: 1; catalytic_steps: begin E + S ⇌ E(S) :: OnlyA; E(S) <--> E(P) :: EqualAI; E + P ⇌ E(P) :: EqualAI end end; am=ER.AllostericMechanism(m); println(length(ER._state_rate_polys(am,:A)))'`
Expected: prints `2` (proves the third value is missing).

- [ ] **Step 3: Implement the D-surfacing** — reproduce `git show 0e4c556 -- src/rate_eq_derivation.jl` exactly. The core change in `_raw_symbolic_rate_polys` (after the `D = [...]` block, ~line 368):

```julia
    i_free = findfirst(f -> isempty(bound(f)) && isempty(residual(f)), enz_species)
    @assert i_free !== nothing "no free-enzyme form (empty bound, empty residual)"
    d_free = _rename_symbols(D[form_to_group[i_free]], rename_map)
```

Return `num, den, d_free` instead of `num, den`. Update the docstring to note the third return. Then thread the extra value through every caller: `_raw_symbolic_rate_polys(M::Type)` returns it; `_raw_rate_expr_and_symbols`, `_rate_v_line`, both `_kcat_forward` methods, and `_dependent_param_exprs` destructure with `num, den, _ = …`; `_state_rate_polys` returns `(num, den, d_free)` and updates its docstring; `_allosteric_num_den_exprs` destructures `num_A_poly, den_A_poly, _ = …` and `num_i_poly, den_i_poly, _ = …` (Task 3 will use the value).

- [ ] **Step 4: Run the test, expect pass** — rerun the Step-2 command; expected `3`. Then run the new testset:

Run: `julia --project -e 'using EnzymeRates; include("test/test_rate_eq_derivation.jl")'` (or the narrower testset runner the file uses). Expected: the `D[g_free] surfaced` testset passes.

- [ ] **Step 5: Guard against equation drift** — confirm no rate equation changed (this task is pure plumbing):

Run: `julia --project -e 'using EnzymeRates; include("test/test_allosteric_golden.jl")'`
Expected: golden testset passes unchanged (byte-identical).

- [ ] **Step 6: Commit**

```bash
git add src/rate_eq_derivation.jl test/test_rate_eq_derivation.jl
git commit -m "Surface free-enzyme segment weight D[g_free] per allosteric state"
```

---

### Task 2: Rendering helpers

**Files:**
- Modify: `src/rate_eq_derivation.jl` — add `_invert_monomial`, `_is_metabolite_free_monomial` near `_mwc_combine`/`_mwc_power_pair` (~line 1530); add `_mwc_cross_weight` (adapt from `eae1c6c`).
- Test: `test/test_rate_eq_derivation.jl`

**Interfaces:**
- Produces:
  - `_invert_monomial(p::POLY)::POLY` — negates exponents and inverts the coefficient of a single-term POLY; errors if `p` is not a monomial.
  - `_is_metabolite_free_monomial(p::POLY, mets::Set{Symbol})::Bool` — true iff `p` is one term whose monomial names no symbol in `mets`.
  - `_mwc_cross_weight(term, d_other_expr, n)` — `term` when `d_other_expr == 1`, else `:($(_power_expr(d_other_expr, n)) * $term)`.

- [ ] **Step 1: Write the failing test**

```julia
@testset "rendering helpers" begin
    ER = EnzymeRates
    k = ER.POLY(ER._mono(:k_ES_to_EP => 1) => 1)
    @test ER._invert_monomial(k) == ER.POLY(ER._mono(:k_ES_to_EP => -1) => 1)
    @test ER._invert_monomial(ER.poly_one()) == ER.poly_one()
    @test ER._is_metabolite_free_monomial(k, Set([:S, :P]))
    kb = ER.poly_mul(k, ER.POLY(ER._mono(:B => 1) => 1))          # k * B (has metabolite)
    @test !ER._is_metabolite_free_monomial(kb, Set([:B]))
    twoterm = ER.poly_add(k, ER.poly_one())                       # k + 1 (not a monomial)
    @test !ER._is_metabolite_free_monomial(twoterm, Set([:S]))
    @test ER._mwc_cross_weight(:foo, 1, 2) == :foo               # no-op when D==1
    @test ER._mwc_cross_weight(:foo, :D, 1) == :(D * foo)
end
```

- [ ] **Step 2: Run it, expect failure** (undefined functions).

Run: `julia --project -e 'using EnzymeRates; include("test/test_rate_eq_derivation.jl")'`
Expected: `UndefVarError: _invert_monomial`.

- [ ] **Step 3: Implement the helpers**

```julia
"""Inverse of a single-term (monomial) `POLY`: negate every exponent and invert
the coefficient. Errors unless `p` is exactly one term."""
function _invert_monomial(p::POLY)
    length(p) == 1 || error("_invert_monomial: not a monomial: $p")
    mono, coef = first(p)
    POLY(_mono((s => -e for (s, e) in mono)...) => inv(coef))
end

"""True when `p` is a single term whose monomial names no symbol in `mets`.
A metabolite-free monomial `D[g_free]` can be divided out of a `POLY` as a
Laurent factor; a metabolite-bearing or multi-term `D` cannot (division would
put a concentration or a rational in a denominator), so it is cross-weighted."""
function _is_metabolite_free_monomial(p::POLY, mets::Set{Symbol})
    length(p) == 1 || return false
    mono, _ = first(p)
    !any(s in mets for (s, _) in mono)
end
```

And, adapted from `eae1c6c` (place near `_mwc_power_pair`):

```julia
"""Cross-weight an MWC state term by the OTHER conformation's free-enzyme weight
`D_other^n` (`n = catalytic_multiplicity`). Restores a common free-enzyme basis
when the two conformations' `D[g_free]` differ and cannot be divided out (a
metabolite-bearing or multi-term `D`). A no-op when `d_other_expr == 1`."""
_mwc_cross_weight(term, d_other_expr, n) =
    d_other_expr == 1 ? term : :($(_power_expr(d_other_expr, n)) * $term)
```

- [ ] **Step 4: Run the test, expect pass.**

Run: `julia --project -e 'using EnzymeRates; include("test/test_rate_eq_derivation.jl")'`
Expected: the `rendering helpers` testset passes.

- [ ] **Step 5: Commit**

```bash
git add src/rate_eq_derivation.jl test/test_rate_eq_derivation.jl
git commit -m "Add MWC free-enzyme rendering helpers (monomial inverse, cross-weight)"
```

---

### Task 3: Three-way rendering in `_allosteric_num_den_exprs`, and flip the gates

**Files:**
- Modify: `src/rate_eq_derivation.jl` — `_allosteric_num_den_exprs` (~line 1595).
- Modify: `test/allosteric_ground_truth.jl` — flip three `@test_broken` → `@test`; switch the `:NonequalAI` gate's reference network to free-flip-only.

**Interfaces:**
- Consumes: `_state_rate_polys(am, state) → (num, den, d_free)` (Task 1); `_invert_monomial`, `_is_metabolite_free_monomial`, `_mwc_cross_weight` (Task 2).

- [ ] **Step 1: Add the free-flip-only `:NonequalAI` reference to the harness.** In `test/allosteric_ground_truth.jl`, add above the `:NonequalAI` gate testset:

```julia
# Formulation-1 reference: only the free enzyme flips conformation. Each
# conformation runs its own catalytic cycle, coupled only through the shared
# free-enzyme pool. Same forms/edges as `biuni_nonequalAI_flux` but with only
# the E_A<->E_I flip — the model this package derives (commit-when-free).
function biuni_nonequalAI_freeflip_flux(kon, koff, KB, KP; k_A, k_I, L, Keq, A, B, P, FAST=1e7)
    krA = k_A * kon * KP / (koff * KB * Keq); krI = k_I * kon * KP / (koff * KB * Keq)
    species = [:E_A, :EA_A, :EAB_A, :EP_A, :E_I, :EA_I, :EAB_I, :EP_I]
    edges = [
        (:E_A, :EA_A, kon*A), (:EA_A, :E_A, koff), (:E_I, :EA_I, kon*A), (:EA_I, :E_I, koff),
        (:EA_A, :EAB_A, FAST*B/KB), (:EAB_A, :EA_A, FAST),
        (:EA_I, :EAB_I, FAST*B/KB), (:EAB_I, :EA_I, FAST),
        (:E_A, :EP_A, FAST*P/KP), (:EP_A, :E_A, FAST),
        (:E_I, :EP_I, FAST*P/KP), (:EP_I, :E_I, FAST),
        (:E_A, :E_I, FAST*L), (:E_I, :E_A, FAST),                  # only free enzyme flips
        (:EAB_A, :EP_A, k_A), (:EP_A, :EAB_A, krA),
        (:EAB_I, :EP_I, k_I), (:EP_I, :EAB_I, krI),
    ]
    cat_edges = [(:EAB_A, :EP_A, k_A, krA), (:EAB_I, :EP_I, k_I, krI)]
    mwc_ground_truth_flux(species, edges, cat_edges, 1.0)
end
```

Add a self-validation testset immediately after it (mirrors the existing `:NonequalAI` harness self-validation, but against `biuni_nonequalAI_freeflip_flux`): assert `L=0` gives the base rate at `k_A`, and `k_I=k_A` gives the base rate independent of `L`. Use the existing `metab_dfree_base_flux` reference and the same random draws as the neighboring self-validation testset.

- [ ] **Step 2: Point the `:NonequalAI` derivation gate at the free-flip reference.** In the `":NonequalAI-catalysis MWC derivation matches mass-action ground truth"` testset (~line 453), change `v_gt = biuni_nonequalAI_flux(...)` to `v_gt = biuni_nonequalAI_freeflip_flux(...)` (same arguments). Leave it `@test` (it must stay green after the fix).

- [ ] **Step 3: Run the ground-truth file, expect the three `:OnlyA`/LDH gates broken and `:NonequalAI` now FAILING (red).**

Run: `julia --project -e 'using EnzymeRates; include("test/allosteric_ground_truth.jl")'`
Expected: uni-`:OnlyA`, multi-`:OnlyA`, metabolite-D still `@test_broken` (1 pass / 5 broken each); `:NonequalAI` now shows failures (raw derivation ≠ free-flip reference). This is the red state the fix turns green.

- [ ] **Step 4: Implement the three-way rendering.** Rewrite the body of `_allosteric_num_den_exprs`. Compute `cat_params`/`cat_mets` first, capture `d_free_A`/`d_free_I`, then branch. Replace lines ~1602–1657 with:

```julia
    num_A_poly, den_A_poly, d_free_A = _state_rate_polys(am, :A)
    cat_params = Set(_state_all_params(_state_mechanism(am, :A),
                                       _state_step_params(am, :A)))
    cat_mets = Set{Symbol}(metabolites(CM()))

    num_i_poly, den_i_poly, d_free_I = _state_rate_polys(am, :I)

    # Formulation-1 per-state free-enzyme normalization. Render the same value
    # three ways by how D[g_free] combines:
    #   D_A == D_I               → raw (identical conformations; the factor cancels)
    #   both metabolite-free monomials → divide Q/D (clean standard-MWC form)
    #   otherwise                → cross-weight by the other state's D^n (polynomial)
    D_A_expr = 1
    D_I_expr = 1
    if d_free_A == d_free_I
        # raw combine — leave the polynomials and D exprs as identities
    elseif _is_metabolite_free_monomial(d_free_A, cat_mets) &&
           _is_metabolite_free_monomial(d_free_I, cat_mets)
        inv_A = _invert_monomial(d_free_A)
        inv_I = _invert_monomial(d_free_I)
        num_A_poly = poly_mul(num_A_poly, inv_A); den_A_poly = poly_mul(den_A_poly, inv_A)
        num_i_poly = poly_mul(num_i_poly, inv_I); den_i_poly = poly_mul(den_i_poly, inv_I)
    else
        D_A_expr = _poly_to_expr(d_free_A, cat_params, cat_mets)
        D_I_expr = _poly_to_expr(d_free_I, cat_params, cat_mets)
    end

    N_A = _poly_to_expr(num_A_poly, cat_params, cat_mets)
    Q_A = _poly_to_expr(den_A_poly, cat_params, cat_mets)
    N_I = _poly_to_expr(num_i_poly, cat_params, cat_mets)
    Q_I = _poly_to_expr(den_i_poly, cat_params, cat_mets)

    reg_Q_A = Any[_reg_site_expr(am, i, false) for i in eachindex(RS)]
    reg_Q_I = Any[_reg_site_expr(am, i, true) for i in eachindex(RS)]
```

Keep `make_num_term`/`make_den_term` unchanged, then combine with the cross-weight (no-op for raw/divide since `D_*_expr == 1`):

```julia
    num_A = _mwc_cross_weight(make_num_term(N_A, Q_A, reg_Q_A), D_I_expr, CatN)
    den_A = _mwc_cross_weight(make_den_term(Q_A, reg_Q_A), D_I_expr, CatN)
    den_I = _mwc_cross_weight(make_den_term(Q_I, reg_Q_I), D_A_expr, CatN)
    full_den = _mwc_combine(den_A, den_I)

    if isempty(num_i_poly)
        num_A, full_den
    else
        num_I = _mwc_cross_weight(make_num_term(N_I, Q_I, reg_Q_I), D_A_expr, CatN)
        _mwc_combine(num_A, num_I), full_den
    end
```

Do NOT add any ping-pong guard. The single free form is found by `findfirst` in Task 1; ping-pong needs no special case.

- [ ] **Step 5: Flip the three broken gates.** In `test/allosteric_ground_truth.jl`, change `@test_broken isapprox(v_code, v_gt; rtol=1e-4)` to `@test isapprox(v_code, v_gt; rtol=1e-4)` in the uni-`:OnlyA` (~line 194), multi-`:OnlyA` (~line 227), and metabolite-D (~line 350) testsets. Remove the now-stale "KNOWN BUG" comments above each.

- [ ] **Step 6: Run the ground-truth file, expect all green.**

Run: `julia --project -e 'using EnzymeRates; include("test/allosteric_ground_truth.jl")'`
Expected: every testset passes — uni-`:OnlyA`, multi-`:OnlyA`, metabolite-D, and `:NonequalAI` (now against the free-flip reference). No `@test_broken` remain in the four derivation gates.

- [ ] **Step 7: Commit**

```bash
git add src/rate_eq_derivation.jl test/allosteric_ground_truth.jl
git commit -m "Normalize each MWC conformation by its free-enzyme weight (formulation 1)"
```

---

### Task 4: Mirror the rendering in `_kcat_forward`

**Files:**
- Modify: `src/rate_eq_derivation.jl` — `_kcat_forward(::AllostericEnzymeMechanism, …)` (~line 942).

**Interfaces:**
- Consumes: `_state_rate_polys` 3-tuple; `_invert_monomial`, `_is_metabolite_free_monomial`.

`_kcat_forward` must render the same value as `rate_equation`, so it applies the same per-state normalization — but at the POLY level (it works on polynomials, not Exprs). The kcat is a saturating-limit number, so raw/divide/cross-weight give the same value; matching `rate_equation`'s choice keeps the saturating-group extraction consistent.

- [ ] **Step 1: Write the failing test** — the kcat-rescaling round-trip must hold after the fix, on a fragmenting mechanism. Append to `test/test_rate_eq_derivation.jl`:

```julia
@testset "kcat consistent with rate_equation under normalization (uni-OnlyA)" begin
    ER = EnzymeRates
    m = @allosteric_mechanism begin
        substrates: S ; products: P ; catalytic_multiplicity: 1
        catalytic_steps: begin
            E + S ⇌ E(S) :: OnlyA ; E(S) <--> E(P) :: EqualAI ; E + P ⇌ E(P) :: EqualAI
        end
    end
    fp = ER.fitted_params(m)
    prm = NamedTuple{(fp..., :Keq, :E_total)}((0.9, 1.3, 2.1, 0.7, 3.0, 1.0))  # K_P_E,K_A_S_E,k,L,Keq,Etot
    rescaled = ER.rescale_parameter_values(m, prm, 5.0)   # ask for kcat = 5.0
    @test isapprox(ER._kcat_forward(m, rescaled), 5.0; rtol=1e-6)
end
```

- [ ] **Step 2: Run it, expect failure** (kcat computed on raw polys disagrees with the normalized `rate_equation`, so the rescale round-trip misses).

Run: `julia --project -e 'using EnzymeRates; include("test/test_rate_eq_derivation.jl")'`
Expected: the new testset fails (kcat off).

- [ ] **Step 3: Apply the poly-level normalization** in `_kcat_forward`. Capture `d_free_A`/`d_free_I`, and before building the `a_param_names`/`i_param_names` and the kcat groups, fold per the same rule:

```julia
    num_A_poly, den_A_poly, d_free_A = _state_rate_polys(am, :A)
    num_I_poly, den_I_poly, d_free_I = _state_rate_polys(am, :I)
    cat_mets = Set{Symbol}(metabolites(CM()))
    if d_free_A == d_free_I
        # raw
    elseif _is_metabolite_free_monomial(d_free_A, cat_mets) &&
           _is_metabolite_free_monomial(d_free_I, cat_mets)
        inv_A = _invert_monomial(d_free_A); inv_I = _invert_monomial(d_free_I)
        num_A_poly = poly_mul(num_A_poly, inv_A); den_A_poly = poly_mul(den_A_poly, inv_A)
        num_I_poly = poly_mul(num_I_poly, inv_I); den_I_poly = poly_mul(den_I_poly, inv_I)
    else
        num_A_poly = poly_mul(num_A_poly, d_free_I); den_A_poly = poly_mul(den_A_poly, d_free_I)
        num_I_poly = poly_mul(num_I_poly, d_free_A); den_I_poly = poly_mul(den_I_poly, d_free_A)
    end
```

Then extend `a_param_names` to include any non-metabolite symbol the folded A-polys now reference (mirrors `eae1c6c`):

```julia
    a_param_names = union(
        Set(_state_all_params(_state_mechanism(am, :A), _state_step_params(am, :A))),
        setdiff(union(_poly_param_syms(num_A_poly), _poly_param_syms(den_A_poly)), cat_mets))
```

Leave the existing `i_param_names = union(a_param_names, …)` line as-is (it already unions the I-poly params).

- [ ] **Step 4: Run the test, expect pass**, and confirm the existing kcat-rescaling suite stays green.

Run: `julia --project -e 'using EnzymeRates; include("test/test_rate_eq_derivation.jl")'`
Expected: the new kcat testset passes and no previously-passing kcat testset regressed.

- [ ] **Step 5: Commit**

```bash
git add src/rate_eq_derivation.jl test/test_rate_eq_derivation.jl
git commit -m "Apply free-enzyme normalization in _kcat_forward for rate-equation consistency"
```

---

### Task 5: Add the metabolite-at-zero and ping-pong gates

**Files:**
- Modify: `test/allosteric_ground_truth.jl`.

**Interfaces:**
- Consumes: `mwc_ground_truth_flux`, `metab_dfree_onlyA_flux`, `biuni_nonequalAI_freeflip_flux`, `metab_dfree_base_flux` (all in the harness).

- [ ] **Step 1: Metabolite-at-zero gate.** Add a testset that evaluates the derived `rate_equation` and the ground truth at `B = 0` for two mechanisms, asserting they agree:
  - metabolite-D `:OnlyA` (dead inactive): at `B=0` both give `0` (trap). Use `metabD` from the metabolite-D testset and `metab_dfree_onlyA_flux(...; B=0.0, ...)`.
  - `:NonequalAI` (productive): at `B=0`, `v ≠ 0` and equals `biuni_nonequalAI_freeflip_flux(...; B=0.0, ...)`.

```julia
@testset "metabolite in D[g_free] at zero concentration" begin
    ER = EnzymeRates
    metabD = @allosteric_mechanism begin
        substrates: S, B ; products: P ; catalytic_multiplicity: 1
        catalytic_steps: begin
            E + S <--> E(S) :: OnlyA ; E(S) + B ⇌ E(S, B) :: EqualAI
            E(S, B) <--> E(P) :: EqualAI ; E + P ⇌ E(P) :: EqualAI
        end
    end
    fp = ER.fitted_params(metabD)   # (:K_P_E,:kon_A_S_E,:koff_A_S_E,:k_EBS_to_EP,:K_B_ES,:L)
    kon,koff,KP,KB,k,L,Keq,S,P = 1.7,1.1,0.9,0.8,2.1,0.7,3.0,1.1,0.6
    prm = NamedTuple{(fp..., :Keq, :E_total)}((KP,kon,koff,k,KB,L,Keq,1.0))
    v0 = real(ER.rate_equation(metabD, (S=S,B=0.0,P=P), prm))
    @test isapprox(v0, metab_dfree_onlyA_flux(kon,koff,KB,KP,k; L=L,Keq=Keq,S=S,B=0.0,P=P); atol=1e-9)
    @test isapprox(v0, 0.0; atol=1e-9)                       # dead inactive traps at B=0

    allo = @allosteric_mechanism begin
        substrates: A, B ; products: P ; catalytic_multiplicity: 1
        catalytic_steps: begin
            E + A <--> E(A) :: EqualAI ; E(A) + B ⇌ E(A, B) :: EqualAI
            E(A, B) <--> E(P) :: NonequalAI ; E + P ⇌ E(P) :: EqualAI
        end
    end
    fpn = ER.fitted_params(allo)
    kon2,koff2,KP2,KB2,kA,kI,L2,Keq2,A2,P2 = 1.7,1.1,0.9,0.8,2.5,0.4,0.7,3.0,1.1,0.9
    d = Dict(:kon_A_E=>kon2,:koff_A_E=>koff2,:K_P_E=>KP2,:K_B_EA=>KB2,
             :k_A_EAB_to_EP=>kA,:k_I_EAB_to_EP=>kI,:L=>L2)
    prm2 = NamedTuple{(fpn..., :Keq, :E_total)}(((d[s] for s in fpn)..., Keq2, 1.0))
    vN = real(ER.rate_equation(allo, (A=A2,B=0.0,P=P2), prm2))
    @test isapprox(vN, biuni_nonequalAI_freeflip_flux(kon2,koff2,KB2,KP2;
        k_A=kA,k_I=kI,L=L2,Keq=Keq2,A=A2,B=0.0,P=P2); rtol=1e-4)
    @test abs(vN) > 1e-3                                     # reverse flux survives
end
```

- [ ] **Step 2: Ping-pong gate — build and SELF-VALIDATE the free-flip network before gating.** Add the eight-form free-flip ping-pong ground truth. Its two catalytic steps (`E(A)→E(res)+P` and `E(B;res)→E+Q`) are `:NonequalAI`; all bindings `:EqualAI`. Derive both reverse rate constants from a single overall `Keq` so both conformations share it. **The self-validation is mandatory: assert `L=0` gives the single-conformation (active-only) rate and `k_I=k_A` gives the base rate independent of `L`.** If it cannot be made to self-validate, STOP and leave a `@test_skip` with a comment — do not gate the derivation against an unvalidated network.

Build a single-conformation ping-pong mass-action reference first (`pingpong_base_flux`, four forms `E,EA,Eres,EBres`, one cycle) and check it against `rate_equation` of the non-allosteric ping-pong (validates the network construction and the parameter mapping). Then build the two-conformation free-flip version reusing those edges, add only the `E_A<->E_I` flip, and self-validate. Only then compare `rate_equation` of the `:NonequalAI` allosteric ping-pong to it.

Because the parameter-name mapping for ping-pong is intricate, extract the derived mechanism's rate-constant values through `_dependent_param_exprs` (as `test/allosteric_ground_truth.jl`'s other gates do for reverses) rather than hand-mapping names.

- [ ] **Step 3: Run the ground-truth file, expect all green** (including the new gates).

Run: `julia --project -e 'using EnzymeRates; include("test/allosteric_ground_truth.jl")'`
Expected: every testset passes.

- [ ] **Step 4: Commit**

```bash
git add test/allosteric_ground_truth.jl
git commit -m "Gate metabolite-at-zero and allosteric ping-pong against the ground truth"
```

---

### Task 6: Performance re-verification

**Files:**
- No source change expected. If the perf gate fails, STOP and report — do not weaken the bound.

- [ ] **Step 1: Run the performance gate** across `MECHANISM_TEST_SPECS`.

Run: `julia --project -e 'using EnzymeRates; include("test/test_rate_eq_derivation.jl")'` and inspect the `test_rate_equation_performance` testsets.
Expected: `allocs == 0` and `t < 120e-9` for every mechanism. The reverted cross-weight held this; the divide branch adds only rate-constant divisions computed once per call.

- [ ] **Step 2:** If green, no commit needed (no change). If red on a specific mechanism, capture the mechanism name and the measured allocs/time and STOP for review — the fix must not trade the perf contract.

---

### Task 7: Regenerate goldens and run the full suite

**Files:**
- Modify: `test/reference/allosteric_golden_reference.txt`.

- [ ] **Step 1: Regenerate the golden reference** — only fragmenting / `:NonequalAI` mechanisms change; all-`:EqualAI` and non-allosteric strings stay byte-identical.

Run:
```bash
julia --project -e 'using EnzymeRates; include("test/mechanism_definitions_for_test_enzyme_derivation.jl"); include("test/test_allosteric_golden.jl"); open("test/reference/allosteric_golden_reference.txt","w") do io; for l in _allosteric_golden_lines(); println(io, l); end; end'
```

- [ ] **Step 2: Review the golden diff** — confirm every changed line is an allosteric mechanism whose inactive fragments or whose catalysis is `:NonequalAI`, and that the new L-terms carry no bare catalytic rate constant beside a dimensionless `1` (the leak is gone). `git diff test/reference/allosteric_golden_reference.txt`.

- [ ] **Step 3: Run the full test suite.**

Run: `julia --project -e 'using Pkg; Pkg.test()'`
Expected: all green. Investigate any failure to root cause — do not delete or weaken tests.

- [ ] **Step 4: Commit**

```bash
git add test/reference/allosteric_golden_reference.txt
git commit -m "Regenerate allosteric golden reference under free-enzyme normalization"
```

---

### Task 8: Rewrite the MWC allostery docs

**Files:**
- Modify: `docs/src/deriving/mwc_allostery.md`.

The `@example` blocks execute the code, so the printed equations already reflect the fix once Tasks 1–7 land. Update the prose to match. Apply the elements-of-style skill.

- [ ] **Step 1: Rewrite the model paragraph** (~lines 84–96, "Derivation of the MWC rate equation"). State formulation 1 plainly: the enzyme flips between conformations only when free (free-enzyme ratio `L`) and commits to one conformation for a catalytic cycle, so the overall rate weights each conformation's catalytic cycle by its free-enzyme population. Note the contrast with a model where intermediates flip mid-cycle, which differs only for `:NonequalAI` catalysis. Keep the existing rapid-equilibrium statement for regulator binding.

- [ ] **Step 2: Replace the formula block** (~lines 98–114). Show the per-conformation normalization and the rendering:

```
v = E_total * num / den
den = (Q_A/D_A)^cat_n * W_A       + L * (Q_I/D_I)^cat_n       * W_I
num = (N_A/D_A)*(Q_A/D_A)^(cat_n-1)*W_A + L*(N_I/D_I)*(Q_I/D_I)^(cat_n-1)*W_I
```

Explain `D_A`, `D_I` as each conformation's free-enzyme segment weight; that `D = 1` for a single-segment graph (so nothing changes for the common case); and that the package renders `Q/D` directly when `D` is a single rate constant (standard-MWC form with a leading `1`) and cross-weights to a polynomial when `D` carries a metabolite. Keep the `N_I = 0` dead-inactive note.

- [ ] **Step 3: Build the docs to confirm the examples execute and doctests pass.**

Run: `julia --project=docs -e 'using Pkg; Pkg.develop(PackageSpec(path=".")); Pkg.instantiate(); include("docs/make.jl")'` (or the repo's documented docs-build command).
Expected: the build succeeds and the `@example` outputs render the normalized equations with no bare rate-constant leak.

- [ ] **Step 4: Commit**

```bash
git add docs/src/deriving/mwc_allostery.md
git commit -m "Document formulation-1 free-enzyme normalization in the MWC derivation"
```

---

## Acceptance

- `test/allosteric_ground_truth.jl` green: uni-`:OnlyA`, multi-`:OnlyA`, metabolite-D, `:NonequalAI` (free-flip reference), metabolite-at-zero, ping-pong. No `@test_broken` remain in the derivation gates.
- `test_rate_equation_performance`: `allocs == 0`, `t < 120e-9` for every mechanism.
- Full `Pkg.test()` green; golden reference regenerated and reviewed.
- Docs rewritten; `@example` outputs show the normalized equations.
