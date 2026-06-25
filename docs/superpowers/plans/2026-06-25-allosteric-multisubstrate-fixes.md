# Allosteric multi-substrate fixes Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Fix three defects in the allosteric rate-equation derivation that make the LDH `identify_rate_equation` run shed ~140k valid MWC candidates.

**Architecture:** All three live in `src/`. Fix A and Fix B are in `src/rate_eq_derivation.jl` (the `@generated` derivation). Fix C is in `src/identify_rate_equation.jl`. Each fix is independently testable; tests use hand-written `@allosteric_mechanism` fixtures verified to reproduce the bug pre-fix.

**Tech Stack:** Julia, `@generated` functions, `Test` stdlib. Run tests with `julia --project`.

**Spec:** `docs/superpowers/specs/2026-06-24-allosteric-multisubstrate-fixes-design.md`.

## Global Constraints

- 92-character line limit, 4-space indentation. Match surrounding code style.
- `rate_equation` MUST stay allocation-free and sub-100ns for every mechanism in `MECHANISM_TEST_SPECS` (enforced by `test_rate_equation_performance`). Fix A only changes dead-inactive-state allosteric bodies, which are not in that set — but Task 4 reruns the perf test to confirm.
- All `Parameter → Symbol` rendering flows through `name(p, m)`; do not introduce stray `Symbol("K…"/"k…"/"V…"/"L…")` literals (AST-walker guard at `test/test_types.jl`).
- Never weaken Canonical Step Form.
- TDD: write the failing test first, watch it fail, implement minimally, watch it pass, commit. Commit after every task.
- The branch is `allosteric-multisubstrate-fixes`. Commit there.

**kcat contract (Fix B):** `_kcat_forward` returns the *peak achievable forward turnover* — `max` over saturating patterns × regulator corners, at products = 0. This is correct, not a heuristic (verified: `max` equals the numerical grid-peak forward rate). Denominator-only substrate-inhibition regimes never match; product-containing matched monomials never win (lower net flux); effector inhibition is handled by the corner-`max`.

---

### Task 1: Fix C — preserve failing mechanism in CSV failure rows

**Files:**
- Modify: `src/identify_rate_equation.jl` (`_failure_row`, ~line 380)
- Test: `test/test_identify_rate_equation.jl`

**Interfaces:**
- Consumes: `EnzymeRates.FitFailure(mech, error::String)`, `EnzymeRates._failure_row(f)::NamedTuple`, `EnzymeRates.compile_mechanism(::Mechanism)`, `EnzymeRates.init_mechanisms(::EnzymeReaction)`.
- Produces: `_failure_row(f).mechanism_type` is the round-trippable parametric `EnzymeMechanism{Sig}` string when the mechanism compiles, else the bare concrete type name.

- [ ] **Step 1: Write the failing test**

Add inside the top-level `@testset "identify_rate_equation"` block in `test/test_identify_rate_equation.jl` (e.g. after the `_rows_to_dataframe with failure row` testset, ~line 192):

```julia
    @testset "failure row preserves round-trippable mechanism" begin
        rxn = @enzyme_reaction begin
            substrates: S[C]
            products:   P[C]
        end
        m = first(EnzymeRates.init_mechanisms(rxn))
        f = EnzymeRates.FitFailure(m, "boom")
        row = EnzymeRates._failure_row(f)
        # Round-trippable parametric Sig, not the bare concrete type name.
        @test row.mechanism_type == string(typeof(EnzymeRates.compile_mechanism(m)))
        @test row.mechanism_type != "EnzymeRates.Mechanism"
        T = Core.eval(EnzymeRates, Meta.parse(row.mechanism_type))
        @test EnzymeRates.Mechanism(T()) == m
        @test row.error == "boom"
    end
```

- [ ] **Step 2: Run test to verify it fails**

Run: `julia --project -e 'using Pkg; Pkg.test()' 2>&1 | grep -A3 "failure row preserves"`
Expected: FAIL — `row.mechanism_type` is `"EnzymeRates.Mechanism"`, so the first `@test` fails.

- [ ] **Step 3: Write minimal implementation**

In `src/identify_rate_equation.jl`, change the `mechanism_type` line in `_failure_row` (currently `mechanism_type = string(typeof(f.mech)),`) to:

```julia
     mechanism_type = try
         string(typeof(compile_mechanism(f.mech)))
     catch
         string(typeof(f.mech))
     end,
```

Update the explanatory comment above `_failure_row` to describe the new behavior: the round-trippable singleton type is recorded when the mechanism compiles; compilation failures fall back to the bare concrete type so the row still identifies the family.

- [ ] **Step 4: Run test to verify it passes**

Run: `julia --project -e 'using Pkg; Pkg.test()' 2>&1 | grep -A3 "failure row preserves"`
Expected: PASS.

- [ ] **Step 5: Commit**

```bash
git add src/identify_rate_equation.jl test/test_identify_rate_equation.jl
git commit -m "fix: record round-trippable mechanism in CSV failure rows"
```

---

### Task 2: Fix A — keep inactive-state assignments in the rate-equation body

**Files:**
- Modify: `src/rate_eq_derivation.jl` — `_build_allosteric_rate_body` (~line 1699) and `rate_equation_string` (~line 1758)
- Test: `test/test_rate_eq_derivation.jl`

**Interfaces:**
- Consumes: `@allosteric_mechanism`, `EnzymeRates.fitted_params(m)`, `rate_equation(m, concs, params, Reduced)`.
- Produces: `rate_equation` evaluates (no `UndefVarError`) for dead-inactive-state allosteric mechanisms.

- [ ] **Step 1: Write the failing test**

Add to `test/test_rate_eq_derivation.jl` (near the other allosteric tests):

```julia
@testset "Fix A: dead-inactive-state allosteric body defines all I-state symbols" begin
    # Random-order allosteric bi-bi with an :OnlyA catalytic step → dead inactive
    # state. Verified pre-fix to crash with `UndefVarError: koff_I_A_E`.
    m = @allosteric_mechanism begin
        substrates: A, B
        products: P, Q
        catalytic_multiplicity: 2
        catalytic_steps: begin
            E + A <--> E(A)        :: NonequalAI
            E + B <--> E(B)        :: NonequalAI
            E(A) + B <--> E(A, B)  :: NonequalAI
            E(B) + A <--> E(A, B)  :: NonequalAI
            E(A, B) <--> E(P, Q)   :: OnlyA
            E(P, Q) <--> E(P) + Q  :: NonequalAI
            E(P, Q) <--> E(Q) + P  :: NonequalAI
            E(P) <--> E + P        :: NonequalAI
            E(Q) <--> E + Q        :: NonequalAI
        end
    end
    pn = EnzymeRates.fitted_params(m)
    params = merge(NamedTuple{pn}(ntuple(_ -> 1.0, length(pn))),
                   (Keq = 1.0, E_total = 1.0))
    concs = (A = 1.0, B = 1.0, P = 1.0, Q = 1.0)
    @test isfinite(rate_equation(m, concs, params, Reduced))
end
```

- [ ] **Step 2: Run test to verify it fails**

Run: `julia --project -e 'using Pkg; Pkg.test()' 2>&1 | grep -A4 "Fix A: dead-inactive"`
Expected: FAIL — `UndefVarError: koff_I_A_E` (or another `*_I_*` symbol) when `rate_equation` runs.

- [ ] **Step 3: Write minimal implementation**

In `src/rate_eq_derivation.jl`, in `_build_allosteric_rate_body`, replace:

```julia
    a_assignments, i_assignments_ = _build_dep_assignments(M_type)
    # When the I-state cycle is broken, i_assignments (I-state Haldanes
    # and :EqualAI catalytic mirrors K_I = K) become dead code — they're
    # only referenced from the L*num_I branch, which is now elided.
    i_assignments = _i_state_dead(M_type()) ? Expr[] : i_assignments_
```

with:

```julia
    a_assignments, i_assignments_ = _build_dep_assignments(M_type)
    # Keep inactive-state assignments unconditionally: the retained Q_I
    # (`L * den_I`) references them. Deps whose RHS touches an :OnlyA symbol
    # are already zeroed in `_build_dep_assignments`, so nothing is undefined.
    i_assignments = i_assignments_
```

And in `rate_equation_string`, replace `i_assignments = _i_state_dead(m) ? Expr[] : i_assignments_` with `i_assignments = i_assignments_`.

- [ ] **Step 4: Run test to verify it passes**

Run: `julia --project -e 'using Pkg; Pkg.test()' 2>&1 | grep -A4 "Fix A: dead-inactive"`
Expected: PASS.

- [ ] **Step 5: Commit**

```bash
git add src/rate_eq_derivation.jl test/test_rate_eq_derivation.jl
git commit -m "fix: keep I-state assignments in dead-inactive-state allosteric body"
```

---

### Task 3: Fix B — `_kcat_forward` over multiple saturating patterns

**Files:**
- Modify: `src/rate_eq_derivation.jl` — allosteric `_kcat_forward` (`@generated`, ~lines 833–1026)
- Test: `test/test_rate_eq_derivation.jl`

**Interfaces:**
- Consumes: `@allosteric_mechanism`, `@enzyme_mechanism`, `EnzymeRates._kcat_forward(m, params)`, `rate_equation`.
- Produces: `_kcat_forward` returns a finite value (no assert) for allosteric mechanisms with `length(a_keys) > 1`, equal to the peak forward rate.

- [ ] **Step 1: Write the failing tests**

Add to `test/test_rate_eq_derivation.jl`:

```julia
@testset "Fix B: _kcat_forward handles multiple saturating patterns" begin
    # Random-order allosteric bi-bi, live I-state. Verified pre-fix to raise
    # "multiple saturating-substrate kcat components (9 found)".
    m = @allosteric_mechanism begin
        substrates: A, B
        products: P, Q
        catalytic_multiplicity: 2
        catalytic_steps: begin
            E + A <--> E(A)        :: NonequalAI
            E + B <--> E(B)        :: NonequalAI
            E(A) + B <--> E(A, B)  :: NonequalAI
            E(B) + A <--> E(A, B)  :: NonequalAI
            E(A, B) <--> E(P, Q)   :: NonequalAI
            E(P, Q) <--> E(P) + Q  :: NonequalAI
            E(P, Q) <--> E(Q) + P  :: NonequalAI
            E(P) <--> E + P        :: NonequalAI
            E(Q) <--> E + Q        :: NonequalAI
        end
    end
    rng = Random.MersenneTwister(1)
    pn = EnzymeRates.fitted_params(m)
    pv = NamedTuple{pn}(Tuple(0.2 + 2 * rand(rng) for _ in pn))
    kc = EnzymeRates._kcat_forward(m, merge(pv, (Keq = 1.0,)))
    @test isfinite(kc)
    # Peak-productive-turnover contract: equals the numerical peak forward rate.
    fp = merge(pv, (Keq = 1.0, E_total = 1.0))
    vmax = maximum(rate_equation(m, (A = x, B = y, P = 0.0, Q = 0.0), fp, Reduced)
                   for x in 10.0 .^ (0:1:9), y in 10.0 .^ (0:1:9))
    @test kc ≈ vmax rtol = 1e-3
end

@testset "Fix B: non-allosteric random-order bi-bi kcat = peak (contract guard)" begin
    # Already-correct non-allosteric path; pins the max contract Fix B mirrors.
    m = @enzyme_mechanism begin
        substrates: A, B
        products: P, Q
        steps: begin
            E + A <--> E(A)
            E + B <--> E(B)
            E(A) + B <--> E(A, B)
            E(B) + A <--> E(A, B)
            E(A, B) <--> E(P, Q)
            E(P, Q) <--> E(P) + Q
            E(P, Q) <--> E(Q) + P
            E(P) <--> E + P
            E(Q) <--> E + Q
        end
    end
    rng = Random.MersenneTwister(5)
    pn = EnzymeRates.fitted_params(m)
    pv = NamedTuple{pn}(Tuple(0.2 + 2 * rand(rng) for _ in pn))
    kc = EnzymeRates._kcat_forward(m, merge(pv, (Keq = 1.0,)))
    fp = merge(pv, (Keq = 1.0, E_total = 1.0))
    vmax = maximum(rate_equation(m, (A = x, B = y, P = 0.0, Q = 0.0), fp, Reduced)
                   for x in 10.0 .^ (0:1:9), y in 10.0 .^ (0:1:9))
    @test kc ≈ vmax rtol = 1e-3
end
```

- [ ] **Step 2: Run tests to verify they fail**

Run: `julia --project -e 'using Pkg; Pkg.test()' 2>&1 | grep -A4 "Fix B:"`
Expected: the allosteric test FAILS with `_kcat_forward: AllostericEnzymeMechanism with multiple saturating-substrate kcat components (9 found) is unsupported`. The non-allosteric guard test PASSES already (it characterizes existing behavior).

- [ ] **Step 3: Rewrite the allosteric `_kcat_forward` body**

In `src/rate_eq_derivation.jl`, replace the region from the `# Choose the saturating active-state met pattern` comment (~line 877) through the end of the regulator-corner result construction (~line 1025, the `return Expr(:block, ...)` just before the function's closing `end`) with the following. This (a) removes the `length(a_keys) == 1` assert, (b) keeps `i_assignments_` unconditionally (Fix A for this site), (c) hoists the regulator-corner setup out, and (d) loops the per-pattern A/B + corner construction over every `met_key`, accumulating one candidate per (pattern, corner) and returning their `max`:

```julia
    # kcat = peak forward turnover at saturation: max over saturating patterns
    # (met_key) and regulator corners. Denominator-only substrate-inhibition
    # patterns never appear in a_keys; product-containing matched patterns are
    # included but never win (product presence lowers net flux).
    a_keys = sort!([k for k in keys(num_A_groups) if haskey(den_A_groups, k)])
    isempty(a_keys) &&
        error("_kcat_forward: AllostericEnzymeMechanism produced no kcat " *
              "components — saturating-substrate pattern not found in numerator")
    empty_set = Set{Symbol}()

    a_assignments, i_assignments_ = _build_dep_assignments(M_type)
    # Keep inactive-state assignments unconditionally: B_I references them, and
    # deps touching an :OnlyA symbol are already zeroed in _build_dep_assignments.
    i_assignments = i_assignments_
    _, indep = _dependent_param_exprs(M_type)
    hw_params = (indep..., :Keq)
    i_state_dead = _i_state_dead(aem)

    # Regulator-corner setup (independent of the saturating pattern).
    all_ligs = AllostericRegulator[]
    for site in am.regulatory_sites
        for lig in site.ligands
            lig in all_ligs || push!(all_ligs, lig)
        end
    end
    n_ligs = length(all_ligs)
    lig_idx = Dict(lig => i - 1 for (i, lig) in enumerate(all_ligs))

    kcat_exprs = Any[]
    for met_key in a_keys
        num_k_A_expr = _poly_to_expr(num_A_groups[met_key], empty_set, empty_set)
        den_k_A_expr = _poly_to_expr(den_A_groups[met_key], empty_set, empty_set)
        num_I_p = get(num_I_groups, met_key, nothing)
        den_I_p = get(den_I_groups, met_key, nothing)
        num_k_I_expr = num_I_p === nothing ? 0 :
            _poly_to_expr(num_I_p, empty_set, empty_set)
        den_k_I_expr = den_I_p === nothing ? 0 :
            _poly_to_expr(den_I_p, empty_set, empty_set)
        i_pattern_dead = den_I_p === nothing

        A_A = CatN == 1 ? num_k_A_expr :
            :($(num_k_A_expr) * $(den_k_A_expr)^$(CatN - 1))
        B_A = :($(den_k_A_expr)^$(CatN))
        if i_pattern_dead
            A_I = 0
            B_I = 0
        else
            A_I = i_state_dead ? 0 :
                  CatN == 1 ? num_k_I_expr :
                              :($(num_k_I_expr) * $(den_k_I_expr)^$(CatN - 1))
            B_I = :($(den_k_I_expr)^$(CatN))
        end

        if isempty(RS)
            push!(kcat_exprs,
                  :($(CatN) * ($(A_A) + L * $(A_I)) / ($(B_A) + L * $(B_I))))
        else
            for mask in 0:(2^n_ligs - 1)
                W_A_factors = Any[]
                W_I_factors = Any[]
                for (site_idx, site) in enumerate(am.regulatory_sites)
                    n_reg = site.multiplicity
                    sat_terms_A = Any[]
                    sat_terms_I = Any[]
                    for (lig, tag) in zip(site.ligands, site.allo_states)
                        if (mask >> lig_idx[lig]) & 1 == 1
                            if tag !== :OnlyI
                                K_A_sym = name(Kreg(site, lig, :A), am)
                                push!(sat_terms_A, :(inv($K_A_sym)))
                            end
                            if tag !== :OnlyA
                                K_I_state = tag === :EqualAI ? :A : :I
                                K_I_sym = name(Kreg(site, lig, K_I_state), am)
                                push!(sat_terms_I, :(inv($K_I_sym)))
                            end
                        end
                    end
                    if !isempty(sat_terms_A)
                        q_A = length(sat_terms_A) == 1 ?
                            sat_terms_A[1] : _nest_binary(:+, sat_terms_A)
                        push!(W_A_factors, _power_expr(q_A, n_reg))
                    end
                    if !isempty(sat_terms_I)
                        q_I = length(sat_terms_I) == 1 ?
                            sat_terms_I[1] : _nest_binary(:+, sat_terms_I)
                        push!(W_I_factors, _power_expr(q_I, n_reg))
                    end
                end
                if isempty(W_A_factors) && isempty(W_I_factors)
                    kcat_expr = :($(CatN) * ($(A_A) + L * $(A_I)) /
                                  ($(B_A) + L * $(B_I)))
                else
                    W_A = isempty(W_A_factors) ? 1 :
                        length(W_A_factors) == 1 ? W_A_factors[1] :
                        _nest_binary(:*, W_A_factors)
                    W_I = isempty(W_I_factors) ? 1 :
                        length(W_I_factors) == 1 ? W_I_factors[1] :
                        _nest_binary(:*, W_I_factors)
                    kcat_expr = :($(CatN) *
                        ($(A_A) * $(W_A) + L * $(A_I) * $(W_I)) /
                        ($(B_A) * $(W_A) + L * $(B_I) * $(W_I)))
                end
                push!(kcat_exprs, kcat_expr)
            end
        end
    end

    result = length(kcat_exprs) == 1 ? kcat_exprs[1] :
        Expr(:call, :max, kcat_exprs...)
    return Expr(:block,
        _destructuring_expr(hw_params, :params),
        a_assignments...,
        i_assignments...,
        :(return $result))
```

Update the function's docstring (~lines 820–832) to state the kcat contract from the Global Constraints block above.

- [ ] **Step 4: Run tests to verify they pass**

Run: `julia --project -e 'using Pkg; Pkg.test()' 2>&1 | grep -A4 "Fix B:"`
Expected: both `Fix B:` testsets PASS.

- [ ] **Step 5: Commit**

```bash
git add src/rate_eq_derivation.jl test/test_rate_eq_derivation.jl
git commit -m "fix: _kcat_forward takes max over multiple saturating patterns"
```

---

### Task 4: Validation — performance contract + full suite

**Files:** none (verification only).

- [ ] **Step 1: Confirm the performance contract holds**

Run: `julia --project -e 'using Pkg; Pkg.test()' 2>&1 | grep -iE "Performance|alloc|Test Summary"`
Expected: the `Performance` testsets pass (0 allocations, < 100 ns), including the MWC dimer (spec 24B). Fix A only altered dead-inactive-state bodies; 24B has a live (`NonequalAI`) inactive state, so it is untouched.

- [ ] **Step 2: Run the full suite**

Run: `julia --project -e 'using Pkg; Pkg.test()'`
Expected: all tests pass, output pristine.

- [ ] **Step 3: Final commit (if any incidental fixes were needed)**

```bash
git status
# only if there are changes:
git add -A && git commit -m "test: allosteric multi-substrate fixes — full suite green"
```
