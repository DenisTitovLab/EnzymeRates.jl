# Rate-equation Numerator Reaction-Cut Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Fix `_compute_numerator` so the rate-equation numerator is the net flux across one complete steady-state reaction-cut (parallel routes summed, series routes counted once), instead of tracking one metabolite through SS steps.

**Architecture:** Replace the metabolite-tracking numerator with a reaction-cut. Build candidate cuts — *metabolite cuts* (steps that bind a substrate / release a product / iso-convert a substrate / iso-produce a product) and *central-species cuts* (steps that produce / consume an iso-step-endpoint form) — excluding dead-end steps that touch a substrate-product mixed complex. A candidate is usable iff all its steps are SS; the numerator is the oriented-flux sum over the fewest-step usable candidate (tie-break → a chemistry/iso cut). No usable candidate ⇒ a complete all-RE catalytic cycle exists (no finite rate) ⇒ raise. No graph traversal and no free-E reference (robust to multiple enzyme conformations). The same fix repairs the allosteric path and kcat (both consume this numerator). A separate construction guard forbids one reaction being both RE and SS.

**Tech Stack:** Julia, `@generated` King-Altman/Cha derivation in `src/rate_eq_derivation.jl`, POLY Laurent algebra, exact-steady-state ODE oracle already present in `test/test_rate_eq_derivation.jl`.

Full design + per-class verification: `docs/superpowers/specs/2026-06-16-numerator-conversion-step-design.md`.

---

## File Structure

- `src/rate_eq_derivation.jl` — replace `_compute_numerator` (lines ~414-493); add helpers `_reaction_step`, `_is_mixed_complex`. This is the whole numerator fix; the allosteric path and kcat call through it unchanged.
- `src/types.jl` — add a construction-time duplicate-step guard in the `Mechanism` (line ~483) and `AllostericMechanism` (line ~523) inner constructors.
- `test/test_rate_eq_derivation.jl` — add two top-level `@testset`s (mixed-RE/SS exact-rate cases; all-RE-cycle raises). Reuses the existing `ode_steady_state_flux`, `raw_to_ode_params`, `random_independent_params_concs` helpers in this file.
- `test/test_types.jl` — add a test for the duplicate-step construction guard.

---

## Background the implementer needs

- `_compute_numerator(mech, enz_name_to_form, step_params, alpha, form_to_group, D, subs_species, prods_species)` returns `(num::POLY, nu_ref::Int)`. Caller `_raw_symbolic_rate_polys` (same file) scales the denominator by `abs(nu_ref)` when `!= 1`; return `nu_ref = 1` (the cut is crossed once per turnover, so `num` is already `v·DEN/E_total`).
- Per-SS-step King-Altman flux is `rf·D[g1] − rr·D[g2]` in the **canonical** `from→to` direction, where `rf = _ss_contrib(poly_sym(name(step_params[idx][1], mech)), m_lhs, i_form, alpha)` and `rr = _ss_contrib(poly_sym(name(step_params[idx][2], mech)), m_rhs, j_form, alpha)`. This is exactly how the current loop builds it (lines 458-470).
- `_step_sides(s)` returns `(name(from_species(s)), name(to_species(s)), m_lhs, m_rhs)`.
- Accessors that exist: `from_species(s)`, `to_species(s)`, `bound_metabolite(s)`, `is_equilibrium(s)`, `bound(species)`, `has_residual(species)`, `name(species)`, `name(p, mech)`, `_flat_steps(mech)`, `_enumerate_species(mech)`. POLY ops: `poly_add`, `poly_sub`, `poly_mul`, `poly_neg`, `poly_zero`, `poly_sym`. Metabolite subtypes: `Substrate`, `Product`, `CompetitiveInhibitor`, `AllostericRegulator`.
- Forward-reaction orientation of a step: iso/chemistry and substrate-binding are forward in the canonical `from→to` direction; product-release is stored canonically as product *binding* (`E+P→EP`), so its forward (release) direction is `to→from`. Oriented flux: `+canon` for chem/bind, `poly_neg(canon)` for release.
- Run one test file fast (full suite is ~11 min): from the repo root,
  `julia --project=test -e 'using TestEnv; TestEnv.activate(); include("test/mechanism_definitions_for_test_enzyme_derivation.jl"); include("test/test_rate_eq_derivation.jl")'`
  — if `TestEnv` is unavailable, run the full suite `julia --project -e 'using Pkg; Pkg.test()'`.

---

### Task 1: Failing exact-rate tests for mixed-RE/SS mechanisms

**Files:**
- Test: `test/test_rate_eq_derivation.jl` (add a new top-level `@testset` near the other top-level ones, e.g. after the block ending at line 1255)

- [ ] **Step 1: Write the failing test**

Add to `test/test_rate_eq_derivation.jl`:

```julia
@testset "Numerator reaction-cut (mixed RE/SS)" begin
    # Random-order, chemistry SS: today undercounts (one binding branch, ~0.60x).
    m_chem_ss = @enzyme_mechanism begin
        substrates: S1, S2
        products: P
        steps: begin
            E + S1 ⇌ E(S1)
            E + S2 ⇌ E(S2)
            E(S1) + S2 <--> E(S1, S2)
            E(S2) + S1 <--> E(S1, S2)
            E(S1, S2) <--> E(P)
            E(P) ⇌ E + P
        end
    end
    # Catalytic isomerization RE, binding+release SS: today ~0.48x.
    m_re_chem = @enzyme_mechanism begin
        substrates: S1, S2
        products: P1, P2
        steps: begin
            E + S1 ⇌ E(S1)
            E + S2 ⇌ E(S2)
            E(S1) + S2 <--> E(S1, S2)
            E(S2) + S1 <--> E(S1, S2)
            E(S1, S2) ⇌ E(P1, P2)
            E(P1, P2) <--> E(P1) + P2
            E(P1, P2) <--> E(P2) + P1
            E(P1) ⇌ E + P1
            E(P2) ⇌ E + P2
        end
    end
    # Single-RE-segment redundant SS binding: today v ≡ 0; true rate is non-zero.
    m_degenerate = @enzyme_mechanism begin
        substrates: S1, S2
        products: P
        steps: begin
            E + S1 <--> E(S1)
            E + S2 ⇌ E(S2)
            E(S1) + S2 ⇌ E(S1, S2)
            E(S2) + S1 ⇌ E(S1, S2)
            E(S1, S2) <--> E(P)
            E(P) ⇌ E + P
        end
    end
    cases = [
        ("random chem-SS", m_chem_ss,   [:S1, :S2, :P]),
        ("RE-chemistry",   m_re_chem,   [:S1, :S2, :P1, :P2]),
        ("redundant bind", m_degenerate,[:S1, :S2, :P]),
    ]
    for (label, m, mets) in cases
        @testset "$label vs ODE steady state" begin
            rng = Random.MersenneTwister(2026)
            @test all(1:10) do _
                new_params, concs, all_params =
                    random_independent_params_concs(m, mets; rng=rng)
                ode_params = raw_to_ode_params(m, all_params)
                v_ode = ode_steady_state_flux(m, ode_params, concs)
                v_ka  = rate_equation(m, concs, new_params)
                isapprox(v_ode, v_ka; rtol=1e-3)
            end
        end
    end
end
```

- [ ] **Step 2: Run to verify it FAILS**

Run: `julia --project -e 'using Pkg; Pkg.test()'` (or the focused-run recipe above).
Expected: the new `@testset "Numerator reaction-cut (mixed RE/SS)"` FAILS all three cases (rate_equation disagrees with the ODE steady state).

- [ ] **Step 3: Commit the failing test**

```bash
git add test/test_rate_eq_derivation.jl
git commit -m "test: rate_equation must match exact steady state for mixed RE/SS mechanisms (red)"
```

---

### Task 2: Implement the reaction-cut numerator

**Files:**
- Modify: `src/rate_eq_derivation.jl` (replace `_compute_numerator`, lines ~414-493; add two helpers just above it)

- [ ] **Step 1: Add the orientation + mixed-complex helpers**

Insert immediately before `_compute_numerator` in `src/rate_eq_derivation.jl`:

```julia
"""
Forward-reaction endpoints and type of a step as `(from::Species, to::Species,
type::Symbol)`, or `nothing` if the step does not advance the reaction (binds an
inhibitor/regulator). `type` is `:chem` (iso substrate→product), `:bind`
(substrate from solution), or `:release` (product to solution). Product-release
steps are stored canonically as product *binding*, so their forward direction is
`to→from`. NB: a step binding a *substrate* to a regulator-bound enzyme is `:bind`
(its bound metabolite is the substrate) — regulator-bound parallel routes are kept.
"""
function _reaction_step(s::Step)
    bm = bound_metabolite(s)
    bm === nothing      && return (from_species(s), to_species(s), :chem)
    bm isa Substrate    && return (from_species(s), to_species(s), :bind)
    bm isa Product      && return (to_species(s), from_species(s), :release)
    return nothing
end

"""
A substrate-product *mixed complex*: a form carrying at least one bound `Substrate`
AND one bound `Product` simultaneously. Such forms are off the catalytic path
(reached only by product rebinding to a substrate complex, or vice versa) and carry
zero net flux; steps touching them are dead-ends and are excluded from cuts. A
regulator does not make a form mixed; ping-pong intermediates hold a `Residual`,
not bound substrate+product, so are not mixed.
"""
_is_mixed_complex(f::Species) =
    any(m -> m isa Substrate, bound(f)) && any(m -> m isa Product, bound(f))
```

- [ ] **Step 2: Replace `_compute_numerator`**

Replace the whole `_compute_numerator` function body with:

```julia
"""
Numerator = net flux across one complete steady-state reaction-cut. Each
per-turnover-conserved "event" is a candidate cut whose SS-step fluxes sum to v:
metabolite cuts (bind a substrate / release a product / iso-convert a substrate /
iso-produce a product) and central-species cuts (produce / consume an iso-step
endpoint form). Dead-end steps (touching a substrate-product mixed complex) are
excluded. A candidate is usable iff all its steps are SS. NUM = oriented-flux sum
over the fewest-step usable candidate (tie-break toward a chemistry/iso cut). No
usable candidate ⇒ a complete all-RE catalytic cycle ⇒ no finite rate ⇒ raise.
Returns `(num, 1)`.
"""
function _compute_numerator(
    mech::Mechanism, enz_name_to_form, step_params,
    alpha, form_to_group, D, subs_species, prods_species,
)
    flat = _flat_steps(mech)
    # Forward-oriented, non-dead-end reaction steps.
    rsteps = NamedTuple[]
    for (idx, (s, _)) in enumerate(flat)
        r = _reaction_step(s); r === nothing && continue
        ff, ft, typ = r
        (_is_mixed_complex(ff) || _is_mixed_complex(ft)) && continue
        push!(rsteps, (idx = idx, ff = ff, ft = ft, typ = typ, s = s))
    end
    subs_in(f)  = Set(name(m) for m in bound(f) if m isa Substrate)
    prods_in(f) = Set(name(m) for m in bound(f) if m isa Product)

    # Candidate cuts: each a Vector of indices into `rsteps`.
    cands = Vector{Int}[]
    add_cand!(pred) = (g = [k for k in eachindex(rsteps) if pred(rsteps[k])];
                       isempty(g) || push!(cands, g))
    for S in subs_species   # metabolite: bind S, or iso-convert S
        add_cand!(r -> r.typ === :bind && name(bound_metabolite(r.s)) == S)
        add_cand!(r -> r.typ === :chem && S in subs_in(r.ff) && !(S in subs_in(r.ft)))
    end
    for P in prods_species  # metabolite: release P, or iso-produce P
        add_cand!(r -> r.typ === :release && name(bound_metabolite(r.s)) == P)
        add_cand!(r -> r.typ === :chem && P in prods_in(r.ft) && !(P in prods_in(r.ff)))
    end
    central = Set{Species}()  # iso-step endpoint forms
    for r in rsteps; r.typ === :chem && (push!(central, r.ff); push!(central, r.ft)); end
    for X in central        # central-species: produce X (ft==X) / consume X (ff==X)
        add_cand!(r -> r.ft == X)
        add_cand!(r -> r.ff == X)
    end

    is_ss(r) = !is_equilibrium(r.s)
    usable = [c for c in cands if all(is_ss(rsteps[k]) for k in c)]
    isempty(usable) && error(
        "rate_equation: no rapid-equilibrium-consistent reaction cut — a complete " *
        "all-RE catalytic cycle exists, so the mechanism has no finite rate.")

    has_chem(c) = any(rsteps[k].typ === :chem for k in c)
    best = usable[argmin(i -> (length(usable[i]), has_chem(usable[i]) ? 0 : 1),
                         eachindex(usable))]
    num = poly_zero()
    for k in best
        r = rsteps[k]; idx = r.idx
        e_lhs, e_rhs, m_lhs, m_rhs = _step_sides(r.s)
        i_form = enz_name_to_form[e_lhs]; j_form = enz_name_to_form[e_rhs]
        g1, g2 = form_to_group[i_form], form_to_group[j_form]
        rf = _ss_contrib(poly_sym(name(step_params[idx][1], mech)), m_lhs, i_form, alpha)
        rr = _ss_contrib(poly_sym(name(step_params[idx][2], mech)), m_rhs, j_form, alpha)
        canon = poly_sub(poly_mul(rf, D[g1]), poly_mul(rr, D[g2]))
        num = poly_add(num, r.typ === :release ? poly_neg(canon) : canon)
    end
    num, 1
end
```

Note: `subs_species`/`prods_species` arrive as `Vector{Symbol}` (the caller wraps
them — see `_raw_symbolic_rate_polys(M::Type)`). Duplicate candidate sets (e.g.
`convert-S1` and `consume-E(S1,S2)` are often the same step set) are harmless. The
dead-end exclusion changes nothing on the simple bug-case mechanisms (no mixed
complexes) but protects the full enumeration's product-rebinding dead-ends.

- [ ] **Step 3: Run the Task 1 tests to verify they PASS**

Run: the focused-run recipe (or full suite).
Expected: `@testset "Numerator reaction-cut (mixed RE/SS)"` — all three cases PASS.

- [ ] **Step 4: Commit**

```bash
git add src/rate_eq_derivation.jl
git commit -m "fix: rate-equation numerator = net flux across an SS reaction-cut"
```

---

### Task 3: Full-suite regression + oracle string re-baseline

**Files:**
- Possibly modify: `test/test_rate_eq_derivation.jl` and/or fixture strings if flat-string/byte-identical oracle snapshots move (the numerator is now written via the chemistry step → standard Vmax form).

- [ ] **Step 1: Run the full test suite**

Run: `julia --project -e 'using Pkg; Pkg.test()'`
Expected: all numeric oracles (Segel/RE/ping-pong/inhibitor/allosteric analytical rate, ODE steady-state, kcat, Haldane-equilibrium, performance allocs==0/<100ns) PASS. JET/Aqua PASS. The only candidate failures are *string* snapshots whose numerator is now written via the chemistry-step constant instead of a binding/release constant — the numeric rate is unchanged.

- [ ] **Step 2: For each failed string snapshot, confirm the rate is unchanged, then re-baseline**

For any failing string/byte-identical test (e.g. `test_rate_equation_string`, the allosteric byte-identical fixture at line ~1619, or flat-string regressions): verify the same spec's numeric `test_analytical_rate`/`test_ode_steadystate` PASS, then update the expected string to the new (chemistry-form) numerator. Do NOT change a string snapshot whose numeric test also fails — that indicates a real regression to investigate, not a re-baseline.

- [ ] **Step 3: Re-run the full suite to confirm green**

Run: `julia --project -e 'using Pkg; Pkg.test()'`
Expected: entire suite PASS.

- [ ] **Step 4: Commit**

```bash
git add -u
git commit -m "test: re-baseline numerator string snapshots to the chemistry-cut form"
```

---

### Task 4: Raise on degenerate all-RE-catalytic-cycle mechanisms

**Files:**
- Test: `test/test_rate_eq_derivation.jl` (add a top-level `@testset`)

- [ ] **Step 1: Write the test**

Add to `test/test_rate_eq_derivation.jl`:

```julia
@testset "Numerator: all-RE catalytic cycle raises" begin
    # Binding stage mixed (S2 SS, S1 RE), release stage mixed (P1 SS, P2 RE),
    # chemistry RE ⇒ a complete all-RE catalytic cycle exists ⇒ no finite rate.
    m_allre = @enzyme_mechanism begin
        substrates: S1, S2
        products: P1, P2
        steps: begin
            E + S1 ⇌ E(S1)
            E + S2 ⇌ E(S2)
            E(S1) + S2 <--> E(S1, S2)
            E(S2) + S1 ⇌ E(S1, S2)
            E(S1, S2) ⇌ E(P1, P2)
            E(P1, P2) ⇌ E(P1) + P2
            E(P1, P2) <--> E(P2) + P1
            E(P1) ⇌ E + P1
            E(P2) ⇌ E + P2
        end
    end
    err = try
        rate_equation_string(m_allre); nothing
    catch e; e end
    @test err isa ErrorException
    @test occursin("all-RE catalytic cycle", err.msg)
end
```

- [ ] **Step 2: Run to verify it PASSES**

Run: the focused-run recipe (or full suite).
Expected: PASS — `_compute_numerator` (from Task 2) already raises with the all-RE-cycle message when no cut qualifies. (If it instead returns `v ≡ 0`, that is the old behavior — re-check Task 2.)

- [ ] **Step 3: Commit**

```bash
git add test/test_rate_eq_derivation.jl
git commit -m "test: degenerate all-RE catalytic-cycle mechanism raises an informative error"
```

---

### Task 5: Construction guard — forbid a reaction being both RE and SS

**Files:**
- Modify: `src/types.jl` — `Mechanism` inner ctor (line ~483) and `AllostericMechanism` inner ctor (line ~523)
- Test: `test/test_types.jl`

- [ ] **Step 1: Write the failing test**

Add to `test/test_types.jl`:

```julia
@testset "reject same reaction as both RE and SS" begin
    err = try
        @enzyme_mechanism begin
            substrates: S
            products: P
            steps: begin
                E + S <--> E(S)
                E + S ⇌ E(S)
                E(S) ⇌ E(P)
                E(P) ⇌ E + P
            end
        end
        nothing
    catch e; e end
    @test err isa ErrorException
    @test occursin("both rapid-equilibrium and steady-state", err.msg)
end
```

- [ ] **Step 2: Run to verify it FAILS**

Run: focused run of `test/test_types.jl` (or full suite).
Expected: FAIL — the mechanism currently constructs without error.

- [ ] **Step 3: Add the guard helper and call it from both constructors**

In `src/types.jl`, add a helper near `_canonical_group_order!` (line ~468):

```julia
"""
Reject a mechanism that contains the same physical reaction as both a
rapid-equilibrium and a steady-state step (e.g. `E + S <--> E(S)` AND
`E + S ⇌ E(S)`): a single reaction cannot be both fast and slow.
"""
function _assert_no_re_ss_duplicate(steps::Vector{Vector{Step}})
    seen = Dict{Tuple{Species, Species, Union{Metabolite, Nothing}}, Bool}()
    for group in steps, s in group
        k = (from_species(s), to_species(s), bound_metabolite(s))
        if haskey(seen, k) && seen[k] != is_equilibrium(s)
            error("Mechanism: reaction $(name(from_species(s))) → " *
                  "$(name(to_species(s))) appears as both rapid-equilibrium " *
                  "and steady-state; a reaction cannot be both.")
        end
        seen[k] = is_equilibrium(s)
    end
end
```

Then call `_assert_no_re_ss_duplicate(steps)` in the `Mechanism` inner ctor (after `permute!(steps, _canonical_group_order!(steps))`, line ~486) and `_assert_no_re_ss_duplicate(cat_steps)` in the `AllostericMechanism` inner ctor (after its `_canonical_group_order!(cat_steps)`, line ~547).

- [ ] **Step 4: Run to verify the test PASSES**

Run: focused run of `test/test_types.jl` (or full suite).
Expected: PASS.

- [ ] **Step 5: Run the full suite to confirm no mechanism legitimately relied on duplicates**

Run: `julia --project -e 'using Pkg; Pkg.test()'`
Expected: entire suite PASS (no existing spec contains an RE+SS duplicate).

- [ ] **Step 6: Commit**

```bash
git add src/types.jl test/test_types.jl
git commit -m "fix: reject a reaction declared both rapid-equilibrium and steady-state"
```

---

## Self-Review

**Spec coverage:**
- Numerator = net flux across an SS reaction-cut, parallel-summed / series-once → Task 2 (`_compute_numerator` candidate cuts: metabolite + central-species, dead-end excluded).
- Fewest-step cut, tie-break toward chemistry → Task 2 (`argmin` on `(length, rank, level)`).
- Raise on no all-SS cut (all-RE catalytic cycle) → Task 2 (`isempty(usable)`), tested in Task 4.
- Allosteric path + kcat fixed for free → covered by the full suite in Task 3 (they call `_compute_numerator`).
- Oracle preservation + string re-baseline → Task 3.
- Construction guard (RE+SS duplicate) → Task 5.
- Faithful-but-degenerate (uni-uni SS-binding) reproduced, not special-cased → no task needed; Task 2's algorithm yields `flux(bind)` for it (only all-SS cut), which equals the exact rate. (Pruning such mechanisms is #1, out of scope.)

**Placeholder scan:** No TBD/“handle edge cases”/“similar to” — every code step is complete.

**Type consistency:** `_reaction_step` returns `(Species, Species, Symbol)`; `_is_mixed_complex` takes a `Species`, returns `Bool`; `rsteps` is a `Vector{NamedTuple}` with fields `idx,ff,ft,typ,s`; `cands`/`usable`/`best` are `Vector{Int}` index lists into `rsteps`; `_compute_numerator` returns `(POLY, Int)` matching the existing caller; `subs_species`/`prods_species` are `Vector{Symbol}` and ARE read (per-metabolite candidates). `_assert_no_re_ss_duplicate` takes `Vector{Vector{Step}}` (the ctor's `steps`/`cat_steps`).

**Known limitation (documented, no current trigger):** completeness is not proven for a mechanism whose only SS bottleneck is neither a metabolite's full bind/release/convert set nor an iso-endpoint's full produce/consume set; such a mechanism would `raise` rather than mis-derive (none constructed). Dead-end steps touching a substrate-product mixed complex are excluded (zero net flux); regulator-bound forms are kept (parallel catalytic routes). Ping-pong covalent intermediates (Residual, not bound substrate+product) are handled by metabolite/iso-endpoint cuts.
