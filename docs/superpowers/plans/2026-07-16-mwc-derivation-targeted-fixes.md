# Targeted MWC derivation fixes — Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Fix the two ping-pong `:OnlyA` derivation defects (err1 NaN, err2 `_kcat_forward` crash) locally, complete the `:OnlyA` thermodynamic guard, and close the two largest oracle gaps — without rewriting the derivation.

**Architecture:** Five surgical changes. Task 1 tightens one predicate in `_reachable_from_free` (I-state pruning only). Task 2 fixes `_kcat_forward`'s group-key construction. Tasks 3-4 add mass-action oracle gates. Task 5 replaces `_onlya_haldane_violation`'s per-row sign test with exact ε-feasibility (constructor-only). Task 6 records the decision. The single-conformation King–Altman, the MWC combine, the normalization, and all `@generated` codegen are untouched.

**Tech Stack:** Julia 1.12, EnzymeRates.jl. Tests via `TestEnv` per-file drivers.

**Spec:** `docs/superpowers/specs/2026-07-16-mwc-derivation-targeted-fixes-design.md`. Read it before Task 1.

## Global Constraints

- 92-character lines, 4-space indentation. Match surrounding style.
- `rate_equation` MUST stay **allocation-free and < 120 ns** per call (`test_rate_equation_performance` in `test/test_rate_eq_derivation.jl`). Tasks 1-5 change **no codegen**, so this gate must stay green **unchanged**. If it moves, STOP — do not weaken the gate.
- `test/reference/allosteric_golden_reference.txt` must stay **byte-identical**. If any block moves, STOP and report — do not regenerate.
- All `Parameter → Symbol` rendering flows through `name(p, m)`. No `Symbol("K…")`/`Symbol("k…")`/`Symbol("V…")`/`Symbol("L…")` literals (AST-walker test at `test/test_types.jl:1577-1644`).
- All 12 `test/allosteric_ground_truth.jl` gates must stay green.
- Run Julia in the **FOREGROUND and WAIT**. Never background-and-yield. Use a 600000 ms timeout.
- Machine has ~5 GiB free RAM / 4 cores. Run **ONE** julia process at a time. Do NOT run full `Pkg.test()` during tasks — it is ~11 min and memory-heavy; it runs once, in Task 7.
- Never skip a pre-commit hook. Commit after each task.
- Do NOT edit a test to make it pass. A failing assertion is a finding for Denis.

## Per-file test driver (use this everywhere)

```bash
cd /home/denis.linux/.julia/dev/EnzymeRates
julia --project -e 'using TestEnv; TestEnv.activate(); include("test/allosteric_ground_truth.jl")' 2>&1 | tail -25
```

For `test_rate_eq_derivation.jl` the shared definitions include is required first:

```bash
cd /home/denis.linux/.julia/dev/EnzymeRates
julia --project -e 'using TestEnv; TestEnv.activate(); using Test, EnzymeRates, LinearAlgebra, Random; include("test/mechanism_definitions_for_test_enzyme_derivation.jl"); include("test/test_rate_eq_derivation.jl")' 2>&1 | tail -25
```

A `@test` failure does NOT abort an `include`. **Grep for a non-zero `Fail`/`Error` column**, not merely for the absence of `ERROR`.

---

## Background the implementer needs

**Formulation 1.** This package derives MWC formulation 1: *only the free enzyme flips* between the active (A) and inactive (I) conformations, ratio `L`; the enzyme commits to one conformation for a whole catalytic cycle. This is why the I-conformation is entered *only* via free `E`, and why a component free `E` cannot reach holds no inactive mass. Formulation 2 (every intermediate flips) differs by 0.1–3% for live `:NonequalAI` and is NOT the model. `test/allosteric_ground_truth.jl:437-499` carries both references and asserts they differ.

**`_reachable_from_free`** (`src/rate_eq_derivation.jl:1173-1193`) is used **only** by `_state_allo_mechanism`'s I-state pruning (`:1214-1234`). Nothing else calls it, so changing it cannot affect the active state.

**`:OnlyA` ⟹ dead inactive cycle.** The constructor guard forces an `:OnlyA` binding to be accompanied by an `:OnlyA` catalytic tag. Measured over 54 uni-uni tag combos: zero constructable `:OnlyA` mechanisms have a live inactive numerator.

---

### Task 1: err1 — tighten the I-state free-root seed

A ping-pong covalent intermediate has an empty `bound` but a non-empty `residual`, so it seeds itself as a free root. The `:OnlyA`-deleted I-graph then keeps a covalent island free `E` cannot reach, no spanning tree rooted at free `E` exists, and `D[g_free] = 0` → `rate_equation` is `NaN` at every concentration.

**Files:**
- Modify: `src/rate_eq_derivation.jl:1167-1193` (`_reachable_from_free` docstring + seed predicate)
- Test: `test/allosteric_ground_truth.jl` (append a new testset at the end)

**Interfaces:**
- Produces: `_reachable_from_free(groups)` seeds only from forms with an empty `bound` AND an empty `residual`. Signature unchanged.

- [ ] **Step 1: Write the failing test.** Append to the end of `test/allosteric_ground_truth.jl`:

```julia
# ── The gate: a ping-pong :OnlyA I-state must keep a reachable free-enzyme root ──
# A covalent intermediate carries no bound metabolite but does carry a residual.
# Seeding I-state reachability from it makes the pruned inactive graph retain a
# covalent island that free E cannot reach, so no spanning tree rooted at free E
# exists and D[g_free] = 0 — the normalization then divides by zero and
# `rate_equation` is NaN at every concentration. Under formulation 1 only free
# enzyme flips, so a component free E cannot reach holds no inactive mass and
# must be stranded. Both mechanisms below are accepted by the `:OnlyA`
# thermodynamic guard — they are valid, and the derivation must handle them.
@testset "ping-pong :OnlyA I-state keeps a reachable free-enzyme root" begin
    err1 = @allosteric_mechanism begin
        substrates: ATP, F6P
        products: ADP, F16BP
        catalytic_multiplicity: 1
        catalytic_steps: begin
            E + ATP ⇌ E(ATP)                                                       :: EqualAI
            E(ATP) <--> E(F16BP; residual = ATP - F16BP)                           :: OnlyA
            E(; residual = ATP - F16BP) + F16BP ⇌ E(F16BP; residual = ATP - F16BP) :: EqualAI
            E(; residual = ATP - F16BP) + F6P ⇌ E(F6P; residual = ATP - F16BP)     :: OnlyA
            E(F6P; residual = ATP - F16BP) ⇌ E(ADP)                                :: EqualAI
            E + ADP ⇌ E(ADP)                                                       :: EqualAI
        end
    end
    am = ER.AllostericMechanism(err1)
    @test ER._onlya_haldane_violation(ER.reaction(am), ER.steps(am),
                                      ER.cat_allo_states(am)) === nothing

    _, _, d_free_I = ER._state_rate_polys(am, :I)
    @test !isempty(d_free_I)

    fp = ER.fitted_params(err1)
    prm = NamedTuple{(fp..., :Keq, :E_total)}(((1.3 for _ in fp)..., 3.0, 1.0))
    concs = (ATP = 1.1, F6P = 0.7, ADP = 0.6, F16BP = 0.9)
    @test isfinite(real(ER.rate_equation(err1, concs, prm)))
    @test isfinite(ER._kcat_forward(err1, prm))
end
```

- [ ] **Step 2: Run it and watch it fail.**

Run the `allosteric_ground_truth.jl` driver above.
Expected: the new testset fails — `!isempty(d_free_I)` is `false`, `isfinite(real(...))` is `false` (NaN), and `_kcat_forward` errors with "produced no kcat components".

- [ ] **Step 3: Implement.** In `src/rate_eq_derivation.jl`, replace the `_reachable_from_free` docstring paragraph and the seed line.

Replace this docstring text (`:1168-1171`):

```
Names of enzyme forms in the connected component of a free (empty-bound) form
over ALL steps of `groups` (rapid-equilibrium and steady-state alike). Seeding
from every empty-bound form covers a ping-pong covalent intermediate, which
carries no bound metabolite and so is its own free root.
```

with:

```
Names of enzyme forms in the connected component of the free enzyme over ALL
steps of `groups` (rapid-equilibrium and steady-state alike). The free enzyme
is the form carrying neither a bound metabolite nor a residual — the only form
that interconverts between conformations under formulation 1, and so the only
root the inactive conformation can be entered through. A ping-pong covalent
intermediate carries a residual and is therefore not a root: a component the
free enzyme cannot reach holds no inactive mass, and leaving it in place would
strand the free-enzyme spanning tree (`D[g_free] = 0`).
```

Replace the seed line (`:1179`):

```julia
    reach = Set{Symbol}(name(f) for f in forms if isempty(bound(f)))
```

with:

```julia
    reach = Set{Symbol}(name(f) for f in forms
                        if isempty(bound(f)) && isempty(residual(f)))
```

- [ ] **Step 4: Run the test; confirm green.**

Run the `allosteric_ground_truth.jl` driver. Expected: **all 13 testsets pass, 0 Fail, 0 Error** (the 12 existing + the new one).

- [ ] **Step 5: Confirm no regression in the derivation suite.**

Run the `test_rate_eq_derivation.jl` driver. Expected: `Enzyme Derivation Tests | 1803 1803` (Pass == Total), **0 Fail / 0 Error**. This includes the golden reference and the perf gate. If the golden moved, STOP and report.

- [ ] **Step 6: Commit.**

```bash
git add src/rate_eq_derivation.jl test/allosteric_ground_truth.jl
git commit -m "Strand unreachable residual islands from the inactive MWC state"
```

---

### Task 2: err2 — `_kcat_forward` must not group on the normalization factor

`d_free_I` can carry a product. The cross-weight multiplies the A-numerator by `d_free_I^n`, so every saturating-substrate group key acquires a product factor; `_kcat_forward` evaluates at products = 0 and filters product-bearing keys, so `a_keys` empties and it throws. `rate_equation` is unaffected — the factor cancels between numerator and denominator (measured: err2's `rate_equation` = 0.0498, finite).

`_kcat_forward`'s own comment (`:960-965`) already states the normalization "is a common factor of the saturating-limit ratio, so it leaves kcat's value unchanged" — it is applied there only to match `rate_equation`'s branch.

**Files:**
- Modify: `src/rate_eq_derivation.jl:941-1100` (`_kcat_forward`)
- Test: `test/allosteric_ground_truth.jl` (append)

**Interfaces:**
- Consumes: Task 1's tightened `_reachable_from_free` (err2 is unaffected by it — its `d_free_I` is product-bearing, not zero).
- Produces: `_kcat_forward(::AllostericEnzymeMechanism, params)` returns a finite value for a product-bearing `d_free_I`. Signature unchanged.

- [ ] **Step 1: Write the failing test.** Append to `test/allosteric_ground_truth.jl`:

```julia
# ── The gate: a product-bearing d_free_I must not empty the kcat group keys ──
# The per-state free-enzyme normalization is a common factor of each
# conformation's saturating-limit ratio, so it cannot change kcat's value — but
# when `d_free_I` carries a product, cross-weighting the A-numerator by it puts
# that product into every saturating-substrate group key. kcat is evaluated at
# products = 0 and filters product-bearing keys, so the key set empties and the
# lookup throws. `rate_equation` is unaffected (the factor cancels), which is
# why this gate checks kcat specifically.
@testset "product-bearing d_free_I keeps kcat components" begin
    err2 = @allosteric_mechanism begin
        substrates: ATP, F6P
        products: ADP, F16BP
        catalytic_multiplicity: 1
        catalytic_steps: begin
            E + ATP ⇌ E(ATP)                                                       :: EqualAI
            E(ATP) <--> E(F16BP; residual = ATP - F16BP)                           :: EqualAI
            E(; residual = ATP - F16BP) + F16BP ⇌ E(F16BP; residual = ATP - F16BP) :: EqualAI
            E(; residual = ATP - F16BP) + F6P ⇌ E(F6P; residual = ATP - F16BP)     :: OnlyA
            E(F6P; residual = ATP - F16BP) ⇌ E(ADP)                                :: OnlyA
            E + ADP ⇌ E(ADP)                                                       :: EqualAI
        end
    end
    am = ER.AllostericMechanism(err2)
    @test ER._onlya_haldane_violation(ER.reaction(am), ER.steps(am),
                                      ER.cat_allo_states(am)) === nothing

    # The premise: d_free_I really does carry a metabolite here.
    _, _, d_free_I = ER._state_rate_polys(am, :I)
    mets = Set{Symbol}(ER.metabolites(typeof(ER.catalytic_mechanism(err2))()))
    @test any(any(s in mets for (s, _) in mono) for (mono, _) in d_free_I)

    fp = ER.fitted_params(err2)
    prm = NamedTuple{(fp..., :Keq, :E_total)}(((1.3 for _ in fp)..., 3.0, 1.0))
    @test isfinite(real(ER.rate_equation(err2, (ATP=1.1, F6P=0.7, ADP=0.6,
                                                F16BP=0.9), prm)))
    @test isfinite(ER._kcat_forward(err2, prm))
    @test ER._kcat_forward(err2, prm) > 0
end
```

- [ ] **Step 2: Run it and watch it fail.**

Run the `allosteric_ground_truth.jl` driver.
Expected: the `_kcat_forward` assertions error with "produced no kcat components". The `rate_equation` and `d_free_I` assertions should already PASS — confirming the defect is confined to kcat.

- [ ] **Step 3: Implement.** In `_kcat_forward` (`src/rate_eq_derivation.jl:941`), the normalization block at `:966-979` mutates `num_A_poly`/`den_A_poly`/`num_I_poly`/`den_I_poly` before `_kcat_groups_from_polys` builds the group keys. Because the factor is common to each state's numerator and denominator, it cancels in the ratio — so drop the cross-weight branch for kcat and group on the un-normalized polynomials.

Replace the `if d_free_A == d_free_I … else … end` block at `:966-979` with:

```julia
    # kcat is a per-state ratio of saturating limits, and the free-enzyme
    # normalization is a common factor of that state's numerator and
    # denominator, so it cancels and cannot change kcat's value. Applying it
    # here would only push `d_free`'s symbols — including a metabolite, when
    # `d_free` is metabolite-bearing — into every saturating-pattern group key,
    # which the products-are-zero filter below would then reject. Group on the
    # un-normalized polys; the metabolite-free monomial case is divided out
    # because it keeps the group keys identical while making the A/I patterns
    # directly comparable.
    if !_is_metabolite_free_monomial(d_free_A, cat_mets) ||
       !_is_metabolite_free_monomial(d_free_I, cat_mets)
        # leave the polys raw — a metabolite-bearing or multi-term `d_free`
        # cancels in each state's ratio and must not enter the group key
    elseif d_free_A != d_free_I
        inv_A = _invert_monomial(d_free_A); inv_I = _invert_monomial(d_free_I)
        num_A_poly = poly_mul(num_A_poly, inv_A); den_A_poly = poly_mul(den_A_poly, inv_A)
        num_I_poly = poly_mul(num_I_poly, inv_I); den_I_poly = poly_mul(den_I_poly, inv_I)
    end
```

- [ ] **Step 4: Run the test; confirm green.**

Run the `allosteric_ground_truth.jl` driver. Expected: all 14 testsets pass, 0 Fail, 0 Error.

- [ ] **Step 5: Confirm kcat's VALUE did not move.**

Run the `test_rate_eq_derivation.jl` driver. Expected: `1803 1803`, 0 Fail / 0 Error. The `kcat consistent with rate_equation` testset is the one that matters — it asserts kcat equals the numerical grid-peak forward rate, so it catches any value change this refactor causes.

**If that testset fails**, the common-factor argument does not hold for the cross-conformation pattern match (the normalization cancels *within* a conformation but not across the `L`-weighted A/I combine). Try the narrower fallback: keep the normalization for the value path, but strip metabolite symbols from the `met_key` only, in `_kcat_groups_from_polys`.

**If neither works, STOP and report to Denis.** Per the spec, this is the one finding that would legitimately reopen the solve-then-limit rewrite. Do NOT weaken the kcat gate and do NOT special-case the test mechanism.

- [ ] **Step 6: Commit.**

```bash
git add src/rate_eq_derivation.jl test/allosteric_ground_truth.jl
git commit -m "Group kcat saturating patterns on un-normalized polynomials"
```

---

### Task 3: n=2 / n=3 concerted-MWC oracle gate

n ≥ 2 currently has **no** mass-action oracle for any family — the single largest validation hole in the derivation. The `^n` cross term `L·N_I·D_I^{n-1}` is live only for `:NonequalAI` (every constructable `:OnlyA` mechanism has a dead inactive cycle), so this gate must use `:NonequalAI` catalysis.

**Files:**
- Modify: `test/allosteric_ground_truth.jl` (append; reuses `mwc_ground_truth_flux` at `:16` and `biuni_nonequalAI_freeflip_flux` at `:442`)

**Interfaces:**
- Consumes: `mwc_ground_truth_flux(species, edges, cat_edges, Etot)`, `biuni_nonequalAI_freeflip_flux(...)`, both already in the file.
- Produces: `biuni_mwc_oligomer_flux(nprot, kon, koff, KB, KP; k_A, k_I, L, Keq, A, B, P, FAST, freeflip)` — an n-protomer concerted-MWC oracle.

- [ ] **Step 1: Write the oracle and its self-validation.** Append to `test/allosteric_ground_truth.jl`:

```julia
# ── n-protomer concerted-MWC oracle (formulation 1) ─────────────────────────
# Concerted: all protomers share one conformation. Within a conformation the
# protomers are independent, so the joint occupancy state is a tuple. Only the
# FULLY-unliganded oligomer flips (`freeflip=true`) — the n-protomer extension
# of the free-flip-only model this package derives. `freeflip=false` flips every
# joint state (the classic per-form-flip MWC, formulation 2) and is kept ONLY as
# a discriminator: it must NOT match the derivation.
const OCC = (:E, :EA, :EAB, :EP)

function biuni_mwc_oligomer_flux(nprot, kon, koff, KB, KP; k_A, k_I, L, Keq,
                                 A, B, P, FAST=1e7, freeflip=true)
    krA = k_A * kon * KP / (koff * KB * Keq)
    krI = k_I * kon * KP / (koff * KB * Keq)
    prot_edges(kX, krX) = [
        (:E, :EA, kon * A), (:EA, :E, koff),
        (:EA, :EAB, FAST * B / KB), (:EAB, :EA, FAST),
        (:E, :EP, FAST * P / KP), (:EP, :E, FAST),
        (:EAB, :EP, kX), (:EP, :EAB, krX),
    ]
    tbl = Dict(:A => prot_edges(k_A, krA), :I => prot_edges(k_I, krI))
    catrate = Dict(:A => (k_A, krA), :I => (k_I, krI))

    occs = collect(Iterators.product(ntuple(_ -> OCC, nprot)...))
    sp(conf, o) = Symbol(conf, "_", join(o, "_"))
    species = Symbol[sp(conf, o) for conf in (:A, :I) for o in occs]
    setidx(o, i, v) = ntuple(j -> j == i ? v : o[j], length(o))

    edges = Tuple{Symbol,Symbol,Float64}[]
    cat_edges = Tuple{Symbol,Symbol,Float64,Float64}[]
    for conf in (:A, :I), o in occs, i in 1:nprot
        for (f, t, r) in tbl[conf]
            o[i] == f || continue
            push!(edges, (sp(conf, o), sp(conf, setidx(o, i, t)), r))
        end
        if o[i] == :EAB
            kf, kr = catrate[conf]
            push!(cat_edges, (sp(conf, o), sp(conf, setidx(o, i, :EP)), kf, kr))
        end
    end
    empty_o = ntuple(_ -> :E, nprot)
    if freeflip
        push!(edges, (sp(:A, empty_o), sp(:I, empty_o), FAST * L))
        push!(edges, (sp(:I, empty_o), sp(:A, empty_o), FAST))
    else
        for o in occs
            push!(edges, (sp(:A, o), sp(:I, o), FAST * L))
            push!(edges, (sp(:I, o), sp(:A, o), FAST))
        end
    end
    mwc_ground_truth_flux(species, edges, cat_edges, 1.0)
end

# Self-validation. Check (c) is the load-bearing one: it pins the oracle to
# formulation 1. Checks (a)/(b) pass for BOTH formulations and so cannot
# distinguish them on their own.
@testset "concerted-MWC oligomer oracle self-validation" begin
    rng = MersenneTwister(11)
    for nprot in (1, 2, 3), _ in 1:3
        kon = 0.5+2rand(rng); koff = 0.5+2rand(rng)
        KB = 0.5+2rand(rng); KP = 0.5+2rand(rng)
        kA = 0.5+2rand(rng); kI = 0.5+2rand(rng); Keq = 2.0+2rand(rng)
        A = 0.5+2rand(rng); B = 0.5+2rand(rng); P = 0.5+2rand(rng)
        L = 0.5+rand(rng)
        base = metab_dfree_base_flux(kon, koff, KB, KP, kA, Keq, A, B, P)

        # (a) L = 0 : inactive unpopulated -> nprot x the single-protomer rate.
        @test isapprox(biuni_mwc_oligomer_flux(nprot, kon, koff, KB, KP;
                k_A=kA, k_I=kI, L=0.0, Keq=Keq, A=A, B=B, P=P),
            nprot * base; rtol=1e-4)

        # (b) k_I = k_A : conformations identical -> L-independent.
        f1 = biuni_mwc_oligomer_flux(nprot, kon, koff, KB, KP;
                k_A=kA, k_I=kA, L=L, Keq=Keq, A=A, B=B, P=P)
        f5 = biuni_mwc_oligomer_flux(nprot, kon, koff, KB, KP;
                k_A=kA, k_I=kA, L=5.0, Keq=Keq, A=A, B=B, P=P)
        @test isapprox(f1, nprot * base; rtol=1e-4)
        @test isapprox(f1, f5; rtol=1e-4)

        # (d) v = 0 at the equilibrium metabolite ratio.
        Peq = Keq * A * B
        @test abs(biuni_mwc_oligomer_flux(nprot, kon, koff, KB, KP;
                k_A=kA, k_I=kI, L=L, Keq=Keq, A=A, B=B, P=Peq)) < 1e-6
    end

    # (c) THE formulation-1 pin: at nprot = 1 the oracle must reproduce the
    #     established free-flip reference, and must NOT equal the per-form-flip
    #     model. Without this, a formulation-2 oracle would pass (a), (b) and (d)
    #     and then disagree with the derivation by 0.1-3% for live :NonequalAI —
    #     a real number that is not a bug.
    args = (1.7, 1.1, 0.8, 0.9)
    kw = (k_A=2.5, k_I=0.4, L=0.7, Keq=3.0, A=1.1, B=0.5, P=0.6)
    @test isapprox(biuni_mwc_oligomer_flux(1, args...; kw...),
                   biuni_nonequalAI_freeflip_flux(args...; kw...); rtol=1e-4)
    @test !isapprox(biuni_mwc_oligomer_flux(1, args...; freeflip=false, kw...),
                    biuni_nonequalAI_freeflip_flux(args...; kw...); rtol=1e-4)
end
```

- [ ] **Step 2: Run it; the oracle must self-validate before it may gate anything.**

Run the `allosteric_ground_truth.jl` driver. Expected: `concerted-MWC oligomer oracle self-validation` passes. If check (c) fails, the oracle is not formulation 1 — fix the oracle, NOT the derivation.

- [ ] **Step 3: Write the derivation gate.** Append:

```julia
# ── The gate: the ^n combine with a LIVE inactive numerator ─────────────────
# `:OnlyA` always yields a dead inactive cycle (the guard forces an `:OnlyA`
# catalytic tag alongside an `:OnlyA` binding), so the numerator cross term
# `L*N_I*D_I^(n-1)` is live only for `:NonequalAI`. This is the only gate that
# exercises it, and the only mass-action gate at n >= 2 for any family.
@testset ":NonequalAI ^n cross term matches multi-protomer ground truth" begin
    rng = MersenneTwister(20260716)
    for nprot in (2, 3)
        allo = nprot == 2 ?
            @allosteric_mechanism(begin
                substrates: A, B ; products: P ; catalytic_multiplicity: 2
                catalytic_steps: begin
                    E + A <--> E(A)        :: EqualAI
                    E(A) + B ⇌ E(A, B)     :: EqualAI
                    E(A, B) <--> E(P)      :: NonequalAI
                    E + P ⇌ E(P)           :: EqualAI
                end
            end) :
            @allosteric_mechanism(begin
                substrates: A, B ; products: P ; catalytic_multiplicity: 3
                catalytic_steps: begin
                    E + A <--> E(A)        :: EqualAI
                    E(A) + B ⇌ E(A, B)     :: EqualAI
                    E(A, B) <--> E(P)      :: NonequalAI
                    E + P ⇌ E(P)           :: EqualAI
                end
            end)
        fp = ER.fitted_params(allo)
        @test fp == (:kon_A_E, :koff_A_E, :K_P_E, :K_B_EA,
                     :k_A_EAB_to_EP, :k_I_EAB_to_EP, :L)
        for _ in 1:6
            kon = 0.5+2rand(rng); koff = 0.5+2rand(rng)
            KP = 0.5+2rand(rng); KB = 0.5+2rand(rng)
            kA = 0.5+2rand(rng); kI = 0.5+2rand(rng)
            L = 0.5+rand(rng); Keq = 2.0+2rand(rng)
            A = 0.5+2rand(rng); B = 0.5+2rand(rng); P = 0.5+2rand(rng)
            d = Dict(:kon_A_E=>kon, :koff_A_E=>koff, :K_P_E=>KP, :K_B_EA=>KB,
                     :k_A_EAB_to_EP=>kA, :k_I_EAB_to_EP=>kI, :L=>L)
            prm = NamedTuple{(fp..., :Keq, :E_total)}(((d[s] for s in fp)..., Keq, 1.0))
            # `rate_equation` is per active site; the oracle is per oligomer.
            v_code = nprot * real(ER.rate_equation(allo, (A=A, B=B, P=P), prm))
            v_gt = biuni_mwc_oligomer_flux(nprot, kon, koff, KB, KP;
                k_A=kA, k_I=kI, L=L, Keq=Keq, A=A, B=B, P=P)
            @test isapprox(v_code, v_gt; rtol=1e-4)
        end
    end
end
```

- [ ] **Step 4: Run; confirm green.**

Run the `allosteric_ground_truth.jl` driver. Expected: green (measured max rel err 3.1e-7 at both n=2 and n=3 — the oracle's finite-FAST floor). The shipped `^n` combine is already correct; this gate pins it.

If it fails, check the param mapping and the per-active-site vs per-oligomer factor (`nprot *`) BEFORE concluding the derivation is wrong.

- [ ] **Step 5: Commit.**

```bash
git add test/allosteric_ground_truth.jl
git commit -m "Gate the MWC ^n cross term against a multi-protomer ground truth"
```

---

### Task 4: ping-pong value gate

`allosteric_ground_truth.jl:596-600` marks the ping-pong ground truth DEFERRED, and the testset at `:601` is only a self-consistency check. The shipped combine is **algebraically exact** for `:NonequalAI` ping-pong (measured 4.25e-16): the free-enzyme flip is the only cut between the two conformation subnetworks, so at steady state it carries zero net flux and sits at ratio `L` for any FAST — the fast-flip limit is exact here, not approximate. Gate at **1e-15**.

**Files:**
- Modify: `test/allosteric_ground_truth.jl:596-600` (delete the DEFERRED comment), append the gate.

**Interfaces:**
- Consumes: `mwc_ground_truth_flux`.
- Produces: `pingpong_nonequalAI_freeflip_flux(...)`.

- [ ] **Step 1: Read the DEFERRED comment** at `test/allosteric_ground_truth.jl:596-600` and the testset at `:601` so the replacement keeps whatever invariants it asserted.

- [ ] **Step 2: Write the oracle + gate.** Append:

```julia
# ── Two-conformation ping-pong bi-bi oracle (formulation 1) ────────────────
# Ping-pong has TWO empty-bound forms — free E and the covalent intermediate F
# — which is the case that broke the reverted cross-weighting fix. Only free E
# flips (formulation 1); F does not. The E_A<->E_I edge is the only cut between
# the two conformation subnetworks, so at steady state it carries zero net flux
# and sits exactly at ratio L for ANY FAST: the fast-flip limit is EXACT here,
# not approximate. That is why this gate runs at rtol 1e-15 rather than 1e-4.
function pingpong_nonequalAI_freeflip_flux(kon_A, koff_A, kon_B, koff_B;
        k_A, k_I, L, Keq, A, B, P, Q, FAST=1e7)
    krA = k_A / Keq
    krI = k_I / Keq
    species = [:E_A, :EA_A, :F_A, :FB_A, :E_I, :EA_I, :F_I, :FB_I]
    edges = Tuple{Symbol,Symbol,Float64}[]
    for (c, kf, kr) in ((:A, k_A, krA), (:I, k_I, krI))
        e, ea = Symbol(:E_, c), Symbol(:EA_, c)
        f, fb = Symbol(:F_, c), Symbol(:FB_, c)
        append!(edges, [
            (e, ea, kon_A * A), (ea, e, koff_A),
            (ea, f, kf), (f, ea, kr * P),
            (f, fb, kon_B * B), (fb, f, koff_B),
            (fb, e, kf), (e, fb, kr * Q),
        ])
    end
    push!(edges, (:E_A, :E_I, FAST * L), (:E_I, :E_A, FAST))
    cat_edges = [(:EA_A, :F_A, k_A, krA * P), (:EA_I, :F_I, k_I, krI * P)]
    mwc_ground_truth_flux(species, edges, cat_edges, 1.0)
end

@testset "allosteric ping-pong :NonequalAI matches mass-action ground truth" begin
    rng = MersenneTwister(20260716)
    # Self-validation: L = 0 -> the active-only ping-pong; k_I = k_A -> L-independent.
    for _ in 1:4
        ka = 0.5+2rand(rng); kb = 0.5+2rand(rng)
        oa = 0.5+2rand(rng); ob = 0.5+2rand(rng)
        kA = 0.5+2rand(rng); kI = 0.5+2rand(rng); Keq = 2.0+2rand(rng)
        A = 0.5+2rand(rng); B = 0.5+2rand(rng)
        P = 0.5+2rand(rng); Q = 0.5+2rand(rng); L = 0.5+rand(rng)
        base = pingpong_nonequalAI_freeflip_flux(ka, oa, kb, ob;
            k_A=kA, k_I=kA, L=0.0, Keq=Keq, A=A, B=B, P=P, Q=Q)
        f1 = pingpong_nonequalAI_freeflip_flux(ka, oa, kb, ob;
            k_A=kA, k_I=kA, L=L, Keq=Keq, A=A, B=B, P=P, Q=Q)
        @test isapprox(f1, base; rtol=1e-4)
        # v = 0 at the equilibrium metabolite ratio (P*Q / (A*B) = Keq^2).
        Qeq = Keq^2 * A * B / P
        @test abs(pingpong_nonequalAI_freeflip_flux(ka, oa, kb, ob;
            k_A=kA, k_I=kI, L=L, Keq=Keq, A=A, B=B, P=P, Q=Qeq)) < 1e-6
    end
end
```

- [ ] **Step 3: Run; confirm the oracle self-validates.**

Run the `allosteric_ground_truth.jl` driver. Expected: green.

**If `v = 0` at equilibrium fails**, the equilibrium condition or the Haldane wiring (`krA = k_A/Keq`) is wrong for this topology — fix the oracle before going further. A ping-pong's two half-reactions each carry their own Haldane; derive the ratio rather than guessing.

- [ ] **Step 4: Add the derivation gate against the shipped `@allosteric_mechanism`.**

Write the matching `@allosteric_mechanism` (the `E(; residual = A - P)` form at `:611` is the shape; make catalysis `:NonequalAI`), map `fitted_params` to the oracle's arguments by meaning, and assert `isapprox(v_code, v_gt; rtol=1e-15)`.

Document the param mapping in a comment, exactly as the `:NonequalAI` gate at `:530-534` does. If the mapping cannot be made to reproduce 1e-15, relax to 1e-4 **only** with a comment recording the measured value, and report the discrepancy — do NOT silently loosen it.

- [ ] **Step 5: Delete the DEFERRED comment** at `:596-600` now that the gap is closed.

- [ ] **Step 6: Run; confirm green, then commit.**

```bash
git add test/allosteric_ground_truth.jl
git commit -m "Gate allosteric ping-pong against a two-conformation ground truth"
```

---

### Task 5: complete the `:OnlyA` guard (Stiemke ε-feasibility)

The per-row sign test is sound but incomplete — its own docstring says so. Measured: it agrees with exact ε-feasibility on all 1042 assignments up to bi-bi, but admits 1074 of 16,384 ter-uni assignments that are provably infeasible; **136 of 13,005** mechanisms are enumeration-reachable at BFS depth 4, first at depth 3 with 6 groups. **Over-rejection is 0 over 17,814 assignments**, so completing it costs no valid mechanism.

Keep the sign test as a **fast sound pre-filter** and run ε-feasibility only when it passes. This preserves today's behaviour exactly for everything the sign test already catches.

**Files:**
- Modify: `src/thermodynamic_constr_for_rate_eq_derivation.jl:398-486` (`_onlya_haldane_violation` + two new helpers)
- Test: `test/test_types.jl` (append near the existing guard tests)

**Interfaces:**
- Produces: `_onlya_haldane_violation(rxn, cat_steps, cat_allo_states) → Union{Nothing, String}`. Signature unchanged.
- Internal: `_rational_nullspace(M) → Matrix{Rational{BigInt}}` (columns span `{x : M·x = 0}`); `_has_strict_positive_combination(N) → Bool` (is `{y : N·y > 0}` nonempty).

- [ ] **Step 1: Write the failing test.** Append to `test/test_types.jl`:

```julia
@testset ":OnlyA guard rejects a multi-cycle ter-substrate inconsistency" begin
    # A full random-order ter binding cube with seven :OnlyA edges. Every single
    # constraint row carries BOTH signs on its :OnlyA eps-exponents, so the
    # per-row sign test sees no violation — but the coupled system has no
    # strictly-positive solution: rows 2, 4 and 5 combine to force the
    # eps-exponent of K_B_EA to zero, i.e. K_I = K_A, contradicting its :OnlyA
    # tag. The inactive cube circulates flux around the E(A)->E(A,B)<-E(B)->
    # E(B,C)<-E(C)->E(A,C)<-E(A) hexagon at equilibrium — perpetual motion.
    @test_throws ErrorException @allosteric_mechanism begin
        substrates: A[C], B[N], C[O] ; products: P[C, N, O]
        catalytic_multiplicity: 1
        catalytic_steps: begin
            E + A ⇌ E(A)             :: OnlyA
            E + B ⇌ E(B)             :: OnlyA
            E + C ⇌ E(C)             :: OnlyA
            E(A) + B ⇌ E(A, B)       :: OnlyA
            E(A) + C ⇌ E(A, C)       :: EqualAI
            E(B) + A ⇌ E(A, B)       :: EqualAI
            E(B) + C ⇌ E(B, C)       :: EqualAI
            E(C) + A ⇌ E(A, C)       :: EqualAI
            E(C) + B ⇌ E(B, C)       :: EqualAI
            E(A, B) + C ⇌ E(A, B, C) :: OnlyA
            E(A, C) + B ⇌ E(A, B, C) :: OnlyA
            E(B, C) + A ⇌ E(A, B, C) :: OnlyA
            E(A, B, C) <--> E(P)     :: OnlyA
            E + P ⇌ E(P)             :: EqualAI
        end
    end
end
```

- [ ] **Step 2: Run it and watch it fail.**

```bash
cd /home/denis.linux/.julia/dev/EnzymeRates
julia --project -e 'using TestEnv; TestEnv.activate(); using Test, EnzymeRates; include("test/test_types.jl")' 2>&1 | tail -25
```

Expected: FAIL — the mechanism constructs fine today (the guard admits it).

- [ ] **Step 3: Implement the helpers.** Add above `_onlya_haldane_violation` in `src/thermodynamic_constr_for_rate_eq_derivation.jl`:

```julia
"""
Basis of `{x : M·x = 0}` as the COLUMNS of the returned matrix. Exact
Gauss-Jordan over the rationals; each non-pivot (free) column yields one basis
vector.
"""
function _rational_nullspace(M::AbstractMatrix{Rational{BigInt}})
    A = copy(M)
    m, n = size(A)
    pivots = Int[]
    r = 1
    for c in 1:n
        r > m && break
        p = findfirst(i -> A[i, c] != 0, r:m)
        p === nothing && continue
        p += r - 1
        A[[r, p], :] = A[[p, r], :]
        A[r, :] = A[r, :] ./ A[r, c]
        for i in 1:m
            (i == r || A[i, c] == 0) && continue
            A[i, :] = A[i, :] .- A[i, c] .* A[r, :]
        end
        push!(pivots, c)
        r += 1
    end
    free = setdiff(1:n, pivots)
    N = zeros(Rational{BigInt}, n, length(free))
    for (k, f) in enumerate(free)
        N[f, k] = 1
        for (i, c) in enumerate(pivots)
            N[c, k] = -A[i, f]
        end
    end
    N
end

"""
True when some `y` makes every row of `N·y` strictly positive — i.e. the open
cone `{y : N·y > 0}` is nonempty. Exact Fourier-Motzkin elimination: to drop
variable `v`, every (positive, negative) row pair is combined with positive
coefficients, which preserves strictness; a row that reduces to all-zero
encodes `0 > 0` and refutes the system. Returns `nothing` if elimination blows
up, so the caller can fall back to the sound per-row test rather than error on
the enumeration hot path.
"""
function _has_strict_positive_combination(N::AbstractMatrix{Rational{BigInt}})
    d = size(N, 2)
    d == 0 && return false
    rows = [collect(N[i, :]) for i in axes(N, 1)]
    for v in d:-1:1
        pos = [r for r in rows if r[v] > 0]
        neg = [r for r in rows if r[v] < 0]
        nxt = [r for r in rows if r[v] == 0]
        for p in pos, q in neg
            push!(nxt, p .* (-q[v]) .+ q .* p[v])
        end
        any(r -> all(iszero, r), nxt) && return false
        length(nxt) > 4000 && return nothing
        rows = nxt
    end
    true
end
```

- [ ] **Step 4: Wire the complete test into the guard.** In `_onlya_haldane_violation`, keep everything up to and including the `onlyA_cols` construction unchanged. Replace the per-row loop (`:470-485`) with:

```julia
    # Per-row sign test first: sound (an all-one-sign row forces a sum of
    # same-signed positive terms to vanish) and cheap, so it keeps the common
    # rejections on the fast path.
    for i in axes(A, 1)
        signs = Set{Int}()
        for (c, mult) in onlyA_cols
            A[i, c] == 0 || push!(signs, mult * A[i, c] > 0 ? 1 : -1)
        end
        isempty(signs) && continue
        length(signs) == 1 || continue
        offenders = sort!([string(columns[c]) for c in keys(onlyA_cols)
                           if A[i, c] != 0])
        return _onlya_violation_message(offenders)
    end

    # The complete condition: the `ε` exponents must admit a strictly-positive
    # solution of `M·w = 0` (Stiemke feasibility). The per-row test above only
    # inspects one row at a time, so from ter-substrate up a multi-cycle coupled
    # inconsistency passes it. `nothing` from the cone test means elimination
    # blew up; fall back to the per-row verdict (sound, just incomplete) rather
    # than reject a mechanism we could not decide.
    cols = sort!(collect(keys(onlyA_cols)))
    M = Rational{BigInt}[onlyA_cols[c] * A[i, c] for i in axes(A, 1), c in cols]
    N = _rational_nullspace(M)
    feasible = size(N, 2) == 0 ? false : _has_strict_positive_combination(N)
    feasible === nothing && return nothing
    feasible && return nothing
    return _onlya_violation_message(sort!([string(columns[c]) for c in cols]))
```

Add the shared message helper above `_onlya_haldane_violation` (the existing message text, extracted so both call sites share it):

```julia
_onlya_violation_message(offenders) =
    "an :OnlyA binding ($(join(offenders, ", "))) leaves a " *
    "thermodynamic (Haldane/Wegscheider) cycle unsatisfiable: the " *
    "inactive conformation cannot close that cycle at finite nonzero " *
    "affinity. Tag the cycle's chemical step :OnlyA, or tag an " *
    "opposing binding :OnlyA so the affinities diverge together."
```

Update the docstring's "The per-row sign test is a sufficient rejection condition, not a complete one… the gap is a checker-completeness contract issue" paragraph (`:421-429`) to state that the per-row test is now a fast pre-filter and the ε-feasibility test is the complete condition.

- [ ] **Step 5: Run the test; confirm the witness is now rejected.**

Run the `test_types.jl` driver. Expected: the new testset passes, **and every existing `test_types.jl` testset stays green**.

`test_types.jl:1779` asserts `uni(:OnlyA, :NonequalAI, :EqualAI)` is a violation. If the exact test now *admits* it, that is a **pre-existing over-rejection by the sign test** (a `:NonequalAI` iso group's free `k_I` ratio can absorb the imbalance, but the guard's `keep` filter drops only `:OnlyA` iso groups). **STOP and report to Denis** — do not edit that test.

- [ ] **Step 6: Confirm no enumeration regression.**

```bash
cd /home/denis.linux/.julia/dev/EnzymeRates
julia --project -e 'using TestEnv; TestEnv.activate(); using Test, EnzymeRates; include("test/test_mechanism_enumeration.jl")' 2>&1 | tail -25
```

Expected: green. Mechanism counts may drop (measured ~1% of allosteric mechanisms at ter). **Any count change must be explained by the guard** — report it in the commit message. If a count moves for a bi-bi-or-smaller reaction, STOP: the sign test and the exact test agree there (1042/1042 measured), so a change means the implementation is wrong.

- [ ] **Step 7: Run the derivation + ground-truth suites; confirm green, then commit.**

```bash
git add src/thermodynamic_constr_for_rate_eq_derivation.jl test/test_types.jl
git commit -m "Complete the :OnlyA guard with exact epsilon-feasibility"
```

---

### Task 6: record the decision

**Files:**
- Modify: `docs/superpowers/specs/2026-07-16-mwc-solve-then-limit-derivation-design.md` (status header)
- Modify: `docs/superpowers/plans/2026-07-16-mwc-solve-then-limit-derivation.md` (status header)
- Modify: `docs/superpowers/findings/2026-07-16-pingpong-onlya-kcat-bug.md` (Options section)

- [ ] **Step 1: Mark the declined spec.** Insert directly under its `# Solve-then-limit MWC allosteric derivation` title:

```markdown
> **Status: DECLINED (2026-07-16).** Superseded by
> `docs/superpowers/specs/2026-07-16-mwc-derivation-targeted-fixes-design.md`,
> which reproduces every claim below and records what measured true. Kept as a
> record. In particular: the normalization rewrite (Task 5) is an algebraic
> no-op; the `:OnlyA` limit equals graph deletion on every constructable
> mechanism; the guard gap is real but already documented in
> `_onlya_haldane_violation`'s docstring, and graph deletion is what makes it
> benign; the three-way `d_free` branch is required by `:NonequalAI` and does
> not go away; and the `^n` cross term is already correct at n=2 and n=3. The
> "12-43%" figure measured 0.08-86% (median 28%, n=400), and the dividing line
> is metabolite-in-`D`, not the cross-weight branch. Task 4's named witness
> (a lone `:OnlyA` edge) is provably impossible — the minimum is 7 of 12 cube
> edges — and the guard rule it proposes admits the witness class anyway,
> because the solve pins `K_I` finite rather than infinite.
```

- [ ] **Step 2: Mark the declined plan.** Insert the same style of banner under its title, pointing at `docs/superpowers/plans/2026-07-16-mwc-derivation-targeted-fixes.md`.

- [ ] **Step 3: Update the findings doc's Options section.** Replace options 1-3 with the measured outcome: option 3 (redesign) was declined; err1 is fixed by tightening `_reachable_from_free`'s seed to exclude residual-bearing forms; err2 is fixed by grouping kcat on un-normalized polynomials. Note that both root causes are confirmed deletion-induced (`d_free_I` is `1` on the undeleted graph for both) but that neither required the rewrite.

- [ ] **Step 4: Commit.**

```bash
git add docs/superpowers/
git commit -m "Record solve-then-limit as declined; point to the targeted fixes"
```

---

### Task 7: full-suite verification

- [ ] **Step 1: Kill any orphans.**

```bash
pgrep -af runtests.jl && pkill -9 -f runtests.jl; free -g | head -2
```

- [ ] **Step 2: Run the full suite ONCE, foreground.**

```bash
cd /home/denis.linux/.julia/dev/EnzymeRates
julia --project -e 'using Pkg; Pkg.test()' 2>&1 | tail -40
```

Expected: 0 failures, 0 errors. ~11 min. Do not run anything else concurrently.

- [ ] **Step 3: Confirm the two non-negotiables explicitly** in the output: `test_rate_equation_performance` (0 allocations, < 120 ns) green, and the allosteric golden byte-identical. If either moved, STOP and report.

- [ ] **Step 4: Report.** Summarize: which of err1/err2 are fixed, any enumeration count change from Task 5 and its explanation, and the LOC delta (`git diff --stat main..HEAD -- src/`).

## Self-Review

- **Spec coverage:** deliverable 1 → Task 1; deliverable 2 → Task 2; deliverable 3 → Task 5; deliverable 4 → Tasks 1, 2 (regression gates), 3 (n≥2 oracle), 4 (ping-pong); deliverable 5 → Task 6. Testing/risks → Tasks 5, 7. Covered.
- **Placeholders:** none. Task 4 Step 4 is the only step without literal code — the mechanism and mapping are discoverable from the two named in-file exemplars (`:611`, `:530-534`), and its acceptance criterion and its "do not silently loosen" rule are both explicit. Task 2 Step 5 names its fallback and its stop condition.
- **Type consistency:** `_rational_nullspace` / `_has_strict_positive_combination` / `_onlya_violation_message` are defined in Task 5 Step 3-4 and used only there. `biuni_mwc_oligomer_flux` is defined in Task 3 Step 1 and used in Task 3 Step 3. `pingpong_nonequalAI_freeflip_flux` is defined and used in Task 4. `_reachable_from_free`'s signature is unchanged.
- **Ordering:** proven fix first (Task 1), then the open one (Task 2) while context is fresh, then the two pure-addition gates (3, 4), then the behaviour-changing guard (5), docs (6), full suite (7).
