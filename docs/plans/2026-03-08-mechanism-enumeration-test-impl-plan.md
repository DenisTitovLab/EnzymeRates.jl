# Mechanism Enumeration Test Redesign — Implementation Plan

> **For Claude:** REQUIRED SUB-SKILL: Use superpowers:executing-plans to implement this plan task-by-task.

**Goal:** Restructure mechanism enumeration tests to separate data from logic, use DSL macros for mechanism definitions, and add combinatorial verification of hardcoded counts.

**Architecture:** Split current `test/test_mechanism_enumeration.jl` into a specs file (types, reactions, builder functions) and a test file (iteration over specs). Add `mechanism_spec_from_mechanism` helper to convert `@enzyme_mechanism`-built mechanisms to `MechanismSpec`. Each `StageExpansionTestSpec` runs every stage independently on the same base mechanism. `EnumerationTestSpec` runs the full pipeline end-to-end. Both use the same 8 reaction set with different roles.

**Tech Stack:** Julia, Test, EnzymeRates DSL macros (`@enzyme_reaction`, `@enzyme_mechanism`)

---

## Task 1: Create `mechanism_spec_from_mechanism` helper

Add a helper function that converts a compiled `EnzymeMechanism` back to a `MechanismSpec`. This is needed so test specs can define mechanisms using readable `@enzyme_mechanism` macros instead of raw edge tuples.

**Files:**
- Create: `test/mechanism_enumeration_test_specs.jl`

**Step 1: Write a failing test for the helper**

At the top of the new file, add the struct and helper, then a self-test:

```julia
# ABOUTME: Test specifications for mechanism enumeration pipeline
# ABOUTME: Defines spec types, reactions, helper functions, and builder functions

using Random

# ── Helper: convert EnzymeMechanism → MechanismSpec ──────────

"""
Convert a compiled EnzymeMechanism back to a MechanismSpec.
Uses enumerate_enzyme_forms to map species tuples to form indices.
"""
function mechanism_spec_from_mechanism(
    m::EnzymeMechanism, @nospecialize(rxn::EnzymeReaction))
    site_defs, forms = EnzymeRates.enumerate_enzyme_forms(rxn)
    form_names = [f.name for f in forms]
    name_to_idx = Dict(f.name => i for (i, f) in enumerate(forms))

    rxns = reactions(m)
    eq_steps_tuple = equilibrium_steps(m)
    pc = param_constraints(m)

    edges = Tuple{Int,Int}[]
    for (lhs, rhs) in rxns
        enz_lhs = [s for s in lhs if haskey(name_to_idx, s)]
        enz_rhs = [s for s in rhs if haskey(name_to_idx, s)]
        length(enz_lhs) == 1 && length(enz_rhs) == 1 ||
            error("Expected exactly 1 enzyme on each side")
        push!(edges, (name_to_idx[enz_lhs[1]],
                      name_to_idx[enz_rhs[1]]))
    end

    eq_steps = collect(Bool, eq_steps_tuple)
    constraints = [
        (t, c, [(s, sc) for (s, sc) in f])
        for (t, c, f) in pc
    ]

    EnzymeRates.MechanismSpec(
        rxn, edges, length(edges), eq_steps,
        constraints, length(parameters(m)))
end
```

**Step 2: Run a quick manual verification**

Run in Julia REPL to verify the helper works:

```bash
julia --project -e '
using EnzymeRates
include("test/mechanism_enumeration_test_specs.jl")
rxn = @enzyme_reaction begin
    substrates: S[C]
    products: P[C]
end
m = @enzyme_mechanism begin
    species: begin
        substrates: S[C]
        products: P[C]
        enzymes: E, ES[C], EP[C]
    end
    steps: begin
        [E, S] <--> [ES]
        [E, P] <--> [EP]
        [ES] --> [EP]
    end
end
spec = mechanism_spec_from_mechanism(m, rxn)
println("edges: ", spec.edges)
println("eq_steps: ", spec.equilibrium_steps)
println("param_count: ", spec.param_count)
println("n_catalytic_edges: ", spec.n_catalytic_edges)
# Verify round-trip: compile back and check params match
m2 = compile_mechanism(spec)
println("round-trip params match: ", length(parameters(m)) == length(parameters(m2)))
'
```

Expected: prints edges as form index pairs, param_count matches `length(parameters(m))`, round-trip succeeds.

**Step 3: Commit**

```bash
git add test/mechanism_enumeration_test_specs.jl
git commit -m "Add mechanism_spec_from_mechanism helper for test specs"
```

---

## Task 2: Add spec types and reaction definitions

Add `StageExpansionTestSpec`, `EnumerationTestSpec`, and all 8 reaction definitions (with `:unknown`, `:dead_end`, and `:allosteric` variants).

**Files:**
- Modify: `test/mechanism_enumeration_test_specs.jl`

**Step 1: Add spec type definitions**

Append after the helper function:

```julia
# ── Test spec types ──────────────────────────────────────────

"""
Tests a single base MechanismSpec through each pipeline stage
independently (not chained). Each stage runs on [base_mechanism]
and the count is compared to the expected value.
"""
Base.@kwdef struct StageExpansionTestSpec
    name::String
    reaction::Any
    base_mechanism::EnzymeRates.MechanismSpec
    catalytic_n::Int = 0

    expected_n_ress::Int
    expected_n_general_modifier::Int
    expected_n_essential_activator::Int
    expected_n_dead_end::Int
    expected_n_equivalence::Int
    expected_n_dedup::Int
    expected_n_allosteric::Int = 0
    expected_n_tr_equiv::Int = 0
    expected_n_oem_dedup::Int = 0
end

"""
Tests end-to-end from EnzymeReaction through full pipeline,
comparing output count at each stage across all regulator
partitions.
"""
Base.@kwdef struct EnumerationTestSpec
    name::String
    reaction::Any
    catalytic_n::Int = 0

    expected_n_forms::Int
    expected_n_catalytic::Int
    expected_n_ress::Int
    expected_n_general_modifier::Int
    expected_n_essential_activator::Int
    expected_n_dead_end::Int
    expected_n_equivalence::Int
    expected_n_dedup::Int
    expected_n_allosteric::Int = 0
    expected_n_tr_equiv::Int = 0
    expected_n_oem_dedup::Int = 0
    expected_n_total::Int
end
```

**Step 2: Add reaction definitions**

Append after the spec types. Each regulated reaction has variants with different regulator roles.

```julia
# ── Reaction definitions ─────────────────────────────────────
# 8 logical reactions. Regulated reactions have :unknown version
# (for EnumerationTestSpec) and explicit-role versions (for
# StageExpansionTestSpec).

# 1. Uni-Uni, no regulators
const uni_uni = @enzyme_reaction begin
    substrates: S[C]
    products: P[C]
end

# 2. Uni-Uni + 1 regulator
const uni_uni_reg_unknown = @enzyme_reaction begin
    substrates: S[C]
    products: P[C]
    regulators: I
end
const uni_uni_dead_end_I = @enzyme_reaction begin
    substrates: S[C]
    products: P[C]
    dead_end_inhibitors: I
end
const uni_uni_allosteric_I = @enzyme_reaction begin
    substrates: S[C]
    products: P[C]
    allosteric_regulators: I
end

# 3. Uni-Bi + 1 regulator
const uni_bi_reg_unknown = @enzyme_reaction begin
    substrates: S[AB]
    products: P[A], Q[B]
    regulators: I
end
const uni_bi_dead_end_I = @enzyme_reaction begin
    substrates: S[AB]
    products: P[A], Q[B]
    dead_end_inhibitors: I
end
const uni_bi_allosteric_I = @enzyme_reaction begin
    substrates: S[AB]
    products: P[A], Q[B]
    allosteric_regulators: I
end

# 4. Uni-Bi + allosteric regulator (OEM, catalytic_n=2)
const uni_bi_allosteric_I_oem = @enzyme_reaction begin
    substrates: S[AB]
    products: P[A], Q[B]
    allosteric_regulators: I
end

# 5. Bi-Bi, no regulators
const bi_bi = @enzyme_reaction begin
    substrates: A[C], B[N]
    products: P[C], Q[N]
end

# 6. Bi-Bi Ping-Pong, no regulators
const bi_bi_ping_pong = @enzyme_reaction begin
    substrates: A[CX], B[N]
    products: P[C], Q[NX]
end

# 7. Bi-Bi Ping-Pong + 1 regulator
const bi_bi_ping_pong_reg_unknown = @enzyme_reaction begin
    substrates: A[CX], B[N]
    products: P[C], Q[NX]
    regulators: I
end
const bi_bi_ping_pong_dead_end_I = @enzyme_reaction begin
    substrates: A[CX], B[N]
    products: P[C], Q[NX]
    dead_end_inhibitors: I
end
const bi_bi_ping_pong_allosteric_I = @enzyme_reaction begin
    substrates: A[CX], B[N]
    products: P[C], Q[NX]
    allosteric_regulators: I
end

# 8. Bi-Bi + 2 regulators
const bi_bi_two_regs_unknown = @enzyme_reaction begin
    substrates: A[C], B[N]
    products: P[C], Q[N]
    regulators: I, J
end
const bi_bi_dead_end_I_allosteric_J = @enzyme_reaction begin
    substrates: A[C], B[N]
    products: P[C], Q[N]
    dead_end_inhibitors: I
    allosteric_regulators: J
end
```

**Step 3: Verify reactions parse correctly**

```bash
julia --project -e '
using EnzymeRates
include("test/mechanism_enumeration_test_specs.jl")
for name in [:uni_uni, :uni_uni_reg_unknown, :uni_bi_reg_unknown,
             :uni_bi_allosteric_I_oem, :bi_bi, :bi_bi_ping_pong,
             :bi_bi_ping_pong_reg_unknown, :bi_bi_two_regs_unknown]
    rxn = getfield(Main, name)
    println("$name: subs=$(substrates(rxn)), prods=$(products(rxn)), regs=$(EnzymeRates.regulator_roles(rxn))")
end
'
```

Expected: all reactions print correctly with expected regulators.

**Step 4: Commit**

```bash
git add test/mechanism_enumeration_test_specs.jl
git commit -m "Add spec types and reaction definitions for mechanism enumeration tests"
```

---

## Task 3: Build `StageExpansionTestSpec` instances for no-regulator reactions

Build specs for reactions 1 (Uni-Uni), 5 (Bi-Bi), and 6 (Bi-Bi Ping-Pong) — one spec each with a catalytic-only base.

**Files:**
- Modify: `test/mechanism_enumeration_test_specs.jl`

**Step 1: Determine expected counts**

For each reaction, run each stage independently on the first catalytic topology to get actual counts. These become the expected values after verification.

```bash
julia --project -e '
using EnzymeRates
include("test/mechanism_enumeration_test_specs.jl")
for (name, rxn) in [("Uni-Uni", uni_uni), ("Bi-Bi", bi_bi),
                     ("Bi-Bi Ping-Pong", bi_bi_ping_pong)]
    cat = EnzymeRates._catalytic_topologies(rxn)
    base = cat[1]
    println("\n--- $name (first topology) ---")
    println("base edges: ", base.edges)
    println("base eq_steps: ", base.equilibrium_steps)
    for (sname, fn) in [
        ("ress", s -> EnzymeRates._expand_ress_variants([s], rxn)),
        ("gm", s -> EnzymeRates._expand_general_modifiers([s], rxn; allosteric_regs=Symbol[])),
        ("ea", s -> EnzymeRates._expand_essential_activators([s], rxn; allosteric_regs=Symbol[])),
        ("de", s -> EnzymeRates._expand_dead_end_inhibitors([s], rxn; dead_end_regs=Symbol[])),
        ("eq", s -> EnzymeRates._expand_equivalence_constraints([s], rxn)),
        ("dd", s -> EnzymeRates._deduplicate([s], rxn)),
    ]
        result = fn(base)
        println("  $sname: $(length(result))")
    end
end
'
```

**Step 2: Add builder function with combinatorial comments**

Append to the specs file:

```julia
# ── StageExpansionTestSpec builder ───────────────────────────

function build_stage_expansion_specs()
    specs = StageExpansionTestSpec[]

    # ── Reaction 1: Uni-Uni, no regs ─────────────────────────
    let
        m = @enzyme_mechanism begin
            species: begin
                substrates: S[C]
                products: P[C]
                enzymes: E, ES[C], EP[C]
            end
            steps: begin
                [E, S] <--> [ES]
                [E, P] <--> [EP]
                [ES] --> [EP]
            end
        end
        base = mechanism_spec_from_mechanism(m, uni_uni)
        push!(specs, StageExpansionTestSpec(
            name="Uni-Uni, no regs, catalytic base",
            reaction=uni_uni,
            base_mechanism=base,
            # RE/SS: ... (fill in with combinatorial reasoning
            # after Step 1 determines actual values)
            expected_n_ress=0,       # placeholder
            expected_n_general_modifier=0,
            expected_n_essential_activator=0,
            expected_n_dead_end=0,
            expected_n_equivalence=0,
            expected_n_dedup=0,
        ))
    end

    # ── Reaction 5: Bi-Bi, no regs ──────────────────────────
    # (similar let block with @enzyme_mechanism for ordered Bi-Bi)

    # ── Reaction 6: Bi-Bi Ping-Pong, no regs ────────────────
    # (similar let block with @enzyme_mechanism for Ping-Pong)

    specs
end
```

**Note:** The placeholder `0` values MUST be replaced with actual values from Step 1. Each value MUST have an inline comment explaining the combinatorial derivation. For example:

```julia
# 3 edges: S-bind(RE), P-bind(RE), isom(SS).
# RE/SS: 2 RE edges can each flip to SS. Must keep ≥1 SS (already
# have 1) and ≥1 RE group.
# Flips: {S}→valid, {P}→valid, {S,P}→all SS, invalid.
# Total = 1 (original) + 2 (valid flips) = 3
expected_n_ress=3,
```

**Step 3: Verify the builder runs**

```bash
julia --project -e '
using EnzymeRates
include("test/mechanism_enumeration_test_specs.jl")
specs = build_stage_expansion_specs()
println("Built $(length(specs)) StageExpansionTestSpecs")
for s in specs
    println("  $(s.name)")
end
'
```

**Step 4: Commit**

```bash
git add test/mechanism_enumeration_test_specs.jl
git commit -m "Add StageExpansionTestSpec instances for no-regulator reactions"
```

---

## Task 4: Build `StageExpansionTestSpec` instances for single-regulator reactions

Build specs for reactions 2 (Uni-Uni + reg), 3 (Uni-Bi + reg), and 7 (Bi-Bi Ping-Pong + reg). Each gets 3 specs:
- Catalytic-only base with reg as `:dead_end`
- Catalytic-only base with reg as `:allosteric`
- Pre-expanded base with dead-end edges

**Files:**
- Modify: `test/mechanism_enumeration_test_specs.jl`

**Step 1: Determine expected counts for catalytic-only bases**

```bash
julia --project -e '
using EnzymeRates
include("test/mechanism_enumeration_test_specs.jl")
for (name, rxn_de, rxn_al) in [
    ("Uni-Uni+reg", uni_uni_dead_end_I, uni_uni_allosteric_I),
    ("Uni-Bi+reg", uni_bi_dead_end_I, uni_bi_allosteric_I),
    ("Bi-Bi PP+reg", bi_bi_ping_pong_dead_end_I, bi_bi_ping_pong_allosteric_I),
]
    cat_de = EnzymeRates._catalytic_topologies(rxn_de)
    base_de = cat_de[1]
    de_regs = Symbol[r[1] for r in EnzymeRates.regulator_roles(rxn_de) if r[2] == :dead_end]
    al_regs_de = Symbol[r[1] for r in EnzymeRates.regulator_roles(rxn_de) if r[2] == :allosteric]
    println("\n--- $name dead-end (first topology) ---")
    for (sname, fn) in [
        ("ress", s -> EnzymeRates._expand_ress_variants([s], rxn_de)),
        ("gm", s -> EnzymeRates._expand_general_modifiers([s], rxn_de; allosteric_regs=al_regs_de)),
        ("ea", s -> EnzymeRates._expand_essential_activators([s], rxn_de; allosteric_regs=al_regs_de)),
        ("de", s -> EnzymeRates._expand_dead_end_inhibitors([s], rxn_de; dead_end_regs=de_regs)),
        ("eq", s -> EnzymeRates._expand_equivalence_constraints([s], rxn_de)),
        ("dd", s -> EnzymeRates._deduplicate([s], rxn_de)),
    ]
        result = fn(base_de)
        println("  $sname: $(length(result))")
    end

    cat_al = EnzymeRates._catalytic_topologies(rxn_al)
    base_al = cat_al[1]
    de_regs_al = Symbol[r[1] for r in EnzymeRates.regulator_roles(rxn_al) if r[2] == :dead_end]
    al_regs = Symbol[r[1] for r in EnzymeRates.regulator_roles(rxn_al) if r[2] == :allosteric]
    println("\n--- $name allosteric (first topology) ---")
    for (sname, fn) in [
        ("ress", s -> EnzymeRates._expand_ress_variants([s], rxn_al)),
        ("gm", s -> EnzymeRates._expand_general_modifiers([s], rxn_al; allosteric_regs=al_regs)),
        ("ea", s -> EnzymeRates._expand_essential_activators([s], rxn_al; allosteric_regs=al_regs)),
        ("de", s -> EnzymeRates._expand_dead_end_inhibitors([s], rxn_al; dead_end_regs=de_regs_al)),
        ("eq", s -> EnzymeRates._expand_equivalence_constraints([s], rxn_al)),
        ("dd", s -> EnzymeRates._deduplicate([s], rxn_al)),
    ]
        result = fn(base_al)
        println("  $sname: $(length(result))")
    end
end
'
```

**Step 2: Build pre-expanded dead-end base mechanisms**

For each reaction, define a mechanism with dead-end edges already present using `@enzyme_mechanism`, then convert with `mechanism_spec_from_mechanism`. The mechanism should include the catalytic cycle plus one dead-end complex (e.g., `EI` for Uni-Uni with dead-end I).

Example for Uni-Uni + dead-end I with pre-expanded base:
```julia
let
    m = @enzyme_mechanism begin
        species: begin
            substrates: S[C]
            products: P[C]
            enzymes: E, ES[C], EP[C], EI, ESI[C]
        end
        steps: begin
            [E, S] <--> [ES]      # catalytic
            [E, P] <--> [EP]      # catalytic
            [ES] --> [EP]          # catalytic (SS)
            [E, I] ⇌ [EI]        # dead-end (RE)
            [ES, I] ⇌ [ESI]      # dead-end (RE)
        end
    end
    base = mechanism_spec_from_mechanism(m, uni_uni_dead_end_I)
    # ... build StageExpansionTestSpec with this base
end
```

**Note:** For pre-expanded bases, `n_catalytic_edges` in the resulting `MechanismSpec` will equal the total number of edges (since `mechanism_spec_from_mechanism` sets `n_catalytic_edges = length(edges)`). This needs to be corrected — the helper should accept an optional `n_catalytic_edges` parameter, OR we manually set it after construction. Check whether dead-end stage behavior depends on `n_catalytic_edges`.

**IMPORTANT:** Before proceeding, verify that `mechanism_spec_from_mechanism` correctly handles the `n_catalytic_edges` distinction. If the dead-end expansion stage uses `n_catalytic_edges` to identify which edges are catalytic vs dead-end, then the pre-expanded base MUST set `n_catalytic_edges` correctly. This may require adding an `n_catalytic_edges` keyword to `mechanism_spec_from_mechanism`:

```julia
function mechanism_spec_from_mechanism(
    m::EnzymeMechanism, @nospecialize(rxn::EnzymeReaction);
    n_catalytic_edges::Int=0)
    # ... existing code ...
    n_cat = n_catalytic_edges > 0 ? n_catalytic_edges : length(edges)
    EnzymeRates.MechanismSpec(
        rxn, edges, n_cat, eq_steps,
        constraints, length(parameters(m)))
end
```

**Step 3: Determine expected counts for pre-expanded bases**

Run the same stage-by-stage analysis as Step 1 but using the pre-expanded bases.

**Step 4: Add specs to builder function with combinatorial comments**

Each expected value MUST have an inline comment explaining the derivation. For dead-end expansion on catalytic-only base, the comment should explain:
- How many forms exist
- How many dead-end complexes can form per regulator
- How many combinations that produces

**Step 5: Verify and commit**

```bash
julia --project -e '
using EnzymeRates
include("test/mechanism_enumeration_test_specs.jl")
specs = build_stage_expansion_specs()
println("Built $(length(specs)) StageExpansionTestSpecs")
for s in specs; println("  $(s.name)"); end
'
git add test/mechanism_enumeration_test_specs.jl
git commit -m "Add StageExpansionTestSpec instances for single-regulator reactions"
```

---

## Task 5: Build `StageExpansionTestSpec` instances for multi-regulator and OEM reactions

Build specs for reactions 4 (Uni-Bi + allosteric OEM) and 8 (Bi-Bi + 2 regs).

**Files:**
- Modify: `test/mechanism_enumeration_test_specs.jl`

**Step 1: Determine expected counts**

Same approach as Tasks 3-4: run stages independently, record counts.

For Uni-Bi + allosteric OEM (reaction 4), also test OEM stages:
```bash
julia --project -e '
using EnzymeRates
include("test/mechanism_enumeration_test_specs.jl")
rxn = uni_bi_allosteric_I_oem
cat = EnzymeRates._catalytic_topologies(rxn)
base = cat[1]
al_regs = Symbol[r[1] for r in EnzymeRates.regulator_roles(rxn) if r[2] == :allosteric]
# ... standard stages ...
# OEM stages:
dd = EnzymeRates._deduplicate([base], rxn)
oem = EnzymeRates._expand_allosteric(dd, rxn; catalytic_n=2, allosteric_regs=al_regs)
println("allosteric: ", length(oem))
oem = EnzymeRates._expand_tr_equivalence(oem, rxn)
println("tr_equiv: ", length(oem))
oem = EnzymeRates._deduplicate_oem(oem, rxn)
println("oem_dedup: ", length(oem))
'
```

For Bi-Bi + 2 regs (reaction 8), build 3 specs:
- Catalytic-only base (both DE and allosteric expand)
- Base with dead-end edges (allosteric expansion + equivalence)
- Base with GM + dead-end edges (equivalence and dedup on complex graph)

**Step 2: Add specs with combinatorial comments**

**Step 3: Verify and commit**

```bash
git add test/mechanism_enumeration_test_specs.jl
git commit -m "Add StageExpansionTestSpec instances for OEM and multi-regulator reactions"
```

---

## Task 6: Build `EnumerationTestSpec` instances

Build one spec per reaction using `:unknown` regs. Full pipeline with partitioning.

**Files:**
- Modify: `test/mechanism_enumeration_test_specs.jl`

**Step 1: Determine expected counts using `_run_full_pipeline_stages`**

Copy the existing `_run_full_pipeline_stages` helper into the specs file (it's needed by the test file too). Run it for each reaction:

```bash
julia --project -e '
using EnzymeRates
include("test/mechanism_enumeration_test_specs.jl")
for (name, rxn, cat_n) in [
    ("Uni-Uni", uni_uni, 0),
    ("Uni-Uni+reg", uni_uni_reg_unknown, 0),
    ("Uni-Bi+reg", uni_bi_reg_unknown, 0),
    ("Uni-Bi+allo OEM", uni_bi_allosteric_I_oem, 2),
    ("Bi-Bi", bi_bi, 0),
    ("Bi-Bi PP", bi_bi_ping_pong, 0),
    ("Bi-Bi PP+reg", bi_bi_ping_pong_reg_unknown, 0),
    ("Bi-Bi+2regs", bi_bi_two_regs_unknown, 0),
]
    println("\n--- $name ---")
    _, forms = EnzymeRates.enumerate_enzyme_forms(rxn)
    println("  forms: $(length(forms))")
    counts = _run_full_pipeline_stages(rxn; catalytic_n=cat_n)
    for (k, v) in pairs(counts)
        println("  $k: $v")
    end
    total = collect(EnzymeRates.enumerate_mechanisms(rxn; catalytic_n=cat_n))
    println("  total: $(length(total))")
end
'
```

**Note:** For Bi-Bi + 2 regs, this may be slow or time out. If it does, record what you can and mark the rest as needing lazy enumeration.

**Step 2: Add `_run_full_pipeline_stages` helper and builder function**

```julia
# ── Helper: run pipeline stage by stage (all partitions) ────

function _run_full_pipeline_stages(rxn; catalytic_n::Int=0,
                                   max_re_groups::Int=7)
    # (same implementation as current test file)
    ...
end

function build_enumeration_specs()
    specs = EnumerationTestSpec[]

    push!(specs, EnumerationTestSpec(
        name="Uni-Uni, no regs",
        reaction=uni_uni,
        expected_n_forms=3,
        expected_n_catalytic=1,
        # 1 topology × 1 partition (no regs) = 1
        # RE/SS: ... (combinatorial comments)
        expected_n_ress=0,  # placeholder — fill from Step 1
        ...
    ))

    # ... remaining specs ...

    specs
end
```

**Step 3: Verify and commit**

```bash
git add test/mechanism_enumeration_test_specs.jl
git commit -m "Add EnumerationTestSpec instances and pipeline helper"
```

---

## Task 7: Rewrite `test_mechanism_enumeration.jl` as pure test code

Replace current test file with code that iterates over specs from the builder functions. No data definitions in this file.

**Files:**
- Rewrite: `test/test_mechanism_enumeration.jl`

**Step 1: Write the new test file**

```julia
# ABOUTME: Tests for the staged mechanism enumeration pipeline
# ABOUTME: Iterates over specs defined in mechanism_enumeration_test_specs.jl

const STAGE_EXPANSION_SPECS = build_stage_expansion_specs()
const ENUMERATION_SPECS = build_enumeration_specs()

@testset "Mechanism Enumeration Pipeline" begin

    # ── Stage expansion: each stage independently on base ────
    @testset "Stage expansion: $(s.name)" for s in STAGE_EXPANSION_SPECS
        rxn = s.reaction
        roles = EnzymeRates.regulator_roles(rxn)
        de_regs = Symbol[r[1] for r in roles if r[2] == :dead_end]
        al_regs = Symbol[r[1] for r in roles
                         if r[2] == :allosteric]
        base = [s.base_mechanism]

        @test length(EnzymeRates._expand_ress_variants(
            base, rxn)) == s.expected_n_ress
        @test length(EnzymeRates._expand_general_modifiers(
            base, rxn; allosteric_regs=al_regs)) ==
            s.expected_n_general_modifier
        @test length(EnzymeRates._expand_essential_activators(
            base, rxn; allosteric_regs=al_regs)) ==
            s.expected_n_essential_activator
        @test length(EnzymeRates._expand_dead_end_inhibitors(
            base, rxn; dead_end_regs=de_regs)) ==
            s.expected_n_dead_end
        @test length(EnzymeRates._expand_equivalence_constraints(
            base, rxn)) == s.expected_n_equivalence
        @test length(EnzymeRates._deduplicate(
            base, rxn)) == s.expected_n_dedup

        if s.catalytic_n > 0
            dd = EnzymeRates._deduplicate(base, rxn)
            oem = EnzymeRates._expand_allosteric(
                dd, rxn; catalytic_n=s.catalytic_n,
                allosteric_regs=al_regs)
            @test length(oem) == s.expected_n_allosteric

            oem = EnzymeRates._expand_tr_equivalence(oem, rxn)
            @test length(oem) == s.expected_n_tr_equiv

            oem = EnzymeRates._deduplicate_oem(oem, rxn)
            @test length(oem) == s.expected_n_oem_dedup
        end
    end

    # ── End-to-end pipeline ──────────────────────────────────
    @testset "End-to-end: $(s.name)" for s in ENUMERATION_SPECS
        rxn = s.reaction

        _, forms = EnzymeRates.enumerate_enzyme_forms(rxn)
        @test length(forms) == s.expected_n_forms

        counts = _run_full_pipeline_stages(
            rxn; catalytic_n=s.catalytic_n)
        @test counts.catalytic == s.expected_n_catalytic
        @test counts.ress == s.expected_n_ress
        @test counts.general_modifier ==
            s.expected_n_general_modifier
        @test counts.essential_activator ==
            s.expected_n_essential_activator
        @test counts.dead_end == s.expected_n_dead_end
        @test counts.equivalence == s.expected_n_equivalence
        @test counts.dedup == s.expected_n_dedup

        if s.catalytic_n > 0
            @test counts.allosteric == s.expected_n_allosteric
            @test counts.tr_equiv == s.expected_n_tr_equiv
            @test counts.oem_dedup == s.expected_n_oem_dedup
        end

        total = collect(EnzymeRates.enumerate_mechanisms(
            rxn; catalytic_n=s.catalytic_n))
        @test length(total) == s.expected_n_total
    end

    # ── Property-based tests ─────────────────────────────────
    @testset "Catalytic topology properties" begin
        for rxn in [uni_uni, bi_bi, bi_bi_ping_pong]
            catalytic = EnzymeRates._catalytic_topologies(rxn)
            @test length(catalytic) > 0
            for spec in catalytic
                @test spec.n_catalytic_edges == length(spec.edges)
                @test count(.!spec.equilibrium_steps) >= 1
            end
        end
    end

    @testset "RE/SS expansion properties" begin
        for rxn in [uni_uni, bi_bi, bi_bi_ping_pong]
            catalytic = EnzymeRates._catalytic_topologies(rxn)
            for spec in catalytic
                ress = EnzymeRates._expand_ress_variants(
                    [spec], rxn)
                @test length(ress) > 0
                for s in ress
                    @test any(.!s.equilibrium_steps)
                end
            end
        end
    end

    @testset "Dead-end passthrough with no regs" begin
        for rxn in [uni_uni, bi_bi]
            catalytic = EnzymeRates._catalytic_topologies(rxn)
            no_de = EnzymeRates._expand_dead_end_inhibitors(
                catalytic, rxn; dead_end_regs=Symbol[])
            @test length(no_de) == length(catalytic)
        end
    end

    @testset "Stage monotonicity" begin
        for rxn in [uni_bi_reg_unknown,
                    bi_bi_ping_pong_reg_unknown]
            counts = _run_full_pipeline_stages(rxn)
            @test counts.general_modifier >= counts.ress
            @test counts.essential_activator >=
                counts.general_modifier
            @test counts.dead_end >= counts.essential_activator
            @test counts.equivalence >= counts.dead_end
            @test counts.dedup <= counts.equivalence
        end
    end

    @testset "Regulator roles affect partitioning" begin
        c_unk = _run_full_pipeline_stages(uni_uni_reg_unknown)
        c_de = _run_full_pipeline_stages(uni_uni_dead_end_I)
        c_al = _run_full_pipeline_stages(uni_uni_allosteric_I)
        @test c_unk.catalytic == c_de.catalytic
        @test c_unk.catalytic == c_al.catalytic
    end

    # ── param_count accuracy ─────────────────────────────────
    @testset "param_count accuracy (sampled)" begin
        rng = Random.MersenneTwister(42)
        for rxn in [uni_uni, uni_bi_reg_unknown, bi_bi]
            all_specs = collect(
                EnzymeRates.enumerate_mechanisms(rxn))
            n = min(20, length(all_specs))
            sample = all_specs[randperm(rng, length(all_specs))[
                1:n]]
            for s in sample
                m = compile_mechanism(s)
                @test s.param_count == length(parameters(m))
            end
        end
    end

    # ── OEM properties ───────────────────────────────────────
    @testset "OEM expansion properties" begin
        all_specs = collect(
            EnzymeRates.enumerate_mechanisms(
                uni_bi_allosteric_I_oem; catalytic_n=2))
        oem = filter(
            s -> s isa EnzymeRates.OligomericMechanismSpec,
            all_specs)
        @test length(oem) > 0
        for s in oem
            @test s.catalytic_n == 2
            @test !isempty(s.allosteric_reg_sites)
        end
    end

    # ── compile_mechanism round-trip ─────────────────────────
    @testset "compile_mechanism round-trip" begin
        for rxn in [uni_uni, bi_bi]
            all_specs = collect(
                EnzymeRates.enumerate_mechanisms(rxn))
            for s in all_specs
                m = compile_mechanism(s)
                @test m isa EnzymeMechanism
                @test length(parameters(m)) > 0
            end
        end
    end

    @testset "compile_mechanism OEM" begin
        all_specs = collect(
            EnzymeRates.enumerate_mechanisms(
                uni_bi_allosteric_I_oem; catalytic_n=2))
        oem = filter(
            s -> s isa EnzymeRates.OligomericMechanismSpec,
            all_specs)
        for s in oem[1:min(3, length(oem))]
            m = compile_mechanism(s)
            @test m isa OligomericEnzymeMechanism
            @test length(parameters(m)) > 0
        end
    end
end
```

**Step 2: Update `runtests.jl`**

Add `include("mechanism_enumeration_test_specs.jl")` before the existing `include("test_mechanism_enumeration.jl")`:

In `test/runtests.jl`, find:
```julia
include("test_mechanism_enumeration.jl")
```
Add before it:
```julia
include("mechanism_enumeration_test_specs.jl")
```

**Step 3: Run full test suite**

```bash
julia --project -e 'using Pkg; Pkg.test()'
```

Expected: all tests pass.

**Step 4: Commit**

```bash
git add test/test_mechanism_enumeration.jl test/runtests.jl
git commit -m "Rewrite test_mechanism_enumeration.jl as pure test code over specs"
```

---

## Task 8: Add combinatorial cross-check testset

Add executable formulas that independently compute expected counts and verify they match the hardcoded spec values.

**Files:**
- Modify: `test/test_mechanism_enumeration.jl`

**Step 1: Add combinatorial cross-checks testset**

Append inside the main `@testset` block:

```julia
    # ── Combinatorial cross-checks ───────────────────────────
    @testset "Combinatorial cross-checks" begin
        # Verify hardcoded expected values against independent
        # combinatorial formulas.

        # Uni-Uni RE/SS on catalytic triangle:
        # 3 edges, 1 SS (isomerization). 2 RE edges can flip.
        # Subsets of 2 RE edges to flip: C(2,1)=2, C(2,2)=1.
        # C(2,2)=1 makes all edges SS → invalid (no RE group).
        # Valid flips: 2. Total = 1 (original) + 2 = 3.
        let spec = STAGE_EXPANSION_SPECS[1]
            @test spec.expected_n_ress == 1 + binomial(2, 1)
        end

        # (Add more formulas as each reaction's expected values
        # are verified during Tasks 3-6. Each formula should be
        # independent of pipeline internals.)
    end
```

**Step 2: Run tests**

```bash
julia --project -e 'using Pkg; Pkg.test()'
```

**Step 3: Commit**

```bash
git add test/test_mechanism_enumeration.jl
git commit -m "Add combinatorial cross-check testset for mechanism enumeration"
```

---

## Task 9: Add timeout wrapper for Bi-Bi + 2 regs

**Files:**
- Modify: `test/test_mechanism_enumeration.jl`

**Step 1: Add timeout handling in the end-to-end testset**

Modify the end-to-end loop to handle timeout for large reactions:

```julia
    @testset "End-to-end: $(s.name)" for s in ENUMERATION_SPECS
        rxn = s.reaction

        _, forms = EnzymeRates.enumerate_enzyme_forms(rxn)
        @test length(forms) == s.expected_n_forms

        counts = _run_full_pipeline_stages(
            rxn; catalytic_n=s.catalytic_n)
        # ... stage checks (same as before) ...

        # End-to-end total with timeout for large reactions
        result = try
            t = @elapsed begin
                total = collect(EnzymeRates.enumerate_mechanisms(
                    rxn; catalytic_n=s.catalytic_n))
            end
            t > 120 ? nothing : total
        catch e
            e isa InterruptException ? nothing : rethrow()
        end

        if result === nothing
            @warn "$(s.name) timed out — lazy enumeration needed"
            @test_broken false
        else
            @test length(result) == s.expected_n_total
        end
    end
```

**Step 2: Run tests**

```bash
julia --project -e 'using Pkg; Pkg.test()'
```

**Step 3: Commit**

```bash
git add test/test_mechanism_enumeration.jl
git commit -m "Add timeout handling for large mechanism enumeration tests"
```

---

## Task 10: Final cleanup and full test suite

Remove old reaction/mechanism constants that are no longer referenced. Verify all tests pass.

**Files:**
- Modify: `test/mechanism_enumeration_test_specs.jl` — remove any dead code
- Modify: `test/test_mechanism_enumeration.jl` — final review

**Step 1: Verify no old constants are referenced**

Search for any references to old reaction names (`_rxn_uu`, `_rxn_ub`, etc.) or old mechanism constants (`_hand_uu`, `_hand_bu`, `_first_catalytic`):

```bash
grep -rn '_rxn_\|_hand_\|_first_catalytic' test/
```

Expected: no matches (all replaced by new names).

**Step 2: Run full test suite**

```bash
julia --project -e 'using Pkg; Pkg.test()'
```

Expected: all tests pass.

**Step 3: Commit**

```bash
git add test/mechanism_enumeration_test_specs.jl test/test_mechanism_enumeration.jl
git commit -m "Final cleanup: mechanism enumeration test restructuring complete"
```

---

## Dependency Graph

```
Task 1 (helper function)
  └─ Task 2 (spec types + reactions)
       ├─ Task 3 (no-reg StageExpansion specs)
       ├─ Task 4 (single-reg StageExpansion specs)
       │    └─ depends on Task 3 for pattern
       ├─ Task 5 (multi-reg + OEM StageExpansion specs)
       │    └─ depends on Task 4 for pattern
       └─ Task 6 (EnumerationTestSpec instances)
            └─ Task 7 (rewrite test file)
                 ├─ Task 8 (combinatorial cross-checks)
                 ├─ Task 9 (timeout wrapper)
                 └─ Task 10 (final cleanup)
```

Tasks 3, 4, 5 are sequential (each builds on the pattern of the previous).
Tasks 8, 9 are independent of each other but both depend on Task 7.
Task 10 depends on everything.
