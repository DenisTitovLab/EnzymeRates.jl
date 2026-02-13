# SPEC: Mechanism Enumeration from EnzymeFormSpecs

## Overview

Add `enumerate_mechanisms(reaction::EnzymeReaction)` that generates all valid `EnzymeMechanism` configurations for a given reaction, using the enzyme forms from `enumerate_enzyme_forms`. Returns a lazy iterator over lightweight `MechanismSpec` structs (not type-parameterized `EnzymeMechanism` instances) to avoid compiling millions of types. Users convert selected candidates to `EnzymeMechanism` on demand via `EnzymeMechanism(spec)`.

The enumeration is factored into independent combinatorial dimensions:
1. **Catalytic topology** — the graph structure of the core catalytic cycle (substrate/product binding, release, isomerization steps)
2. **Dead-end attachments** — which regulator-bound forms attach to each cycle form
3. **RE/SS assignments** — which steps are rapid-equilibrium vs steady-state
4. **Symmetry constraints** — constrained (K_i = K_j) and unconstrained variants for structurally symmetric steps

The iterator lazily computes the Cartesian product across dimensions without materializing all combinations.

---

## API

### Function Signature

```julia
enumerate_mechanisms(reaction::EnzymeReaction;
                     max_forms::Int = 2 * n_sites(reaction)
                    ) -> MechanismIterator
```

- `reaction`: the `EnzymeReaction` defining substrates, products, regulators with atoms and max binding sites.
- `max_forms`: maximum number of enzyme forms (including dead-end forms) in any single mechanism. Default = 2 × number of distinct binding sites in the `EnzymeFormSpec` (i.e., total sites from core substrates + core products + extra substrates + extra products + regulators). User can override.
- Returns a `MechanismIterator` that yields `MechanismSpec` structs.

### Helper

```julia
n_sites(reaction::EnzymeReaction)  # total binding sites (core + extra + regulator)
```

### MechanismSpec

Lightweight data struct (not type-parameterized):

```julia
struct MechanismSpec
    species::NTuple{4, Vector{Tuple{Symbol, Vector{Tuple{Symbol,Int}}}}}
    # (substrates, products, regulators, enzymes) — each a vector of (name, atoms)
    reactions::Vector{Tuple{NTuple{N,Symbol} where N, NTuple{M,Symbol} where M}}
    # Each reaction is (lhs_symbols, rhs_symbols)
    equilibrium_steps::Vector{Bool}
    param_constraints::Vector{Tuple{Symbol, Int, Vector{Tuple{Symbol,Int}}}}
end
```

### Conversion

```julia
EnzymeMechanism(spec::MechanismSpec) -> EnzymeMechanism{...}
```

Converts the lightweight spec to the type-parameterized singleton, triggering compilation of the rate equation machinery. Should be called only for mechanisms that pass filtering/selection.

### Iterator Protocol

```julia
Base.iterate(iter::MechanismIterator) -> (MechanismSpec, state) or nothing
Base.length(iter::MechanismIterator) -> Int  # total count across all dimensions
Base.eltype(::Type{MechanismIterator}) = MechanismSpec
```

`length` computes the total count from the factored dimensions without materializing:
```
total = sum over topologies of (dead_end_options × re_ss_options × constraint_options)
```

---

## Elementary Reaction Rules

Given two `EnzymeFormSpec` forms F1 and F2, a valid elementary reaction exists between them if exactly one of the following holds:

### 1. Binding

F1 and F2 differ in exactly one site: that site is unoccupied (`nothing`) in F1 and occupied (has atoms) in F2. All other sites are identical. The metabolite is the one associated with the changed site.

Reaction: `[F1_name, metabolite] → [F2_name]`

### 2. Release

Reverse of binding: F2 has one site unoccupied that is occupied in F1.

Reaction: `[F1_name] → [F2_name, metabolite]`

### 3. Isomerization (catalytic conversion)

F1 and F2 differ only in **core catalytic sites** (index-1 sites of substrates and products). All non-catalytic sites (extra substrate/product sites with index ≥ 2, and all regulator sites) must be identical between F1 and F2. The total atom content across all core catalytic sites must be conserved (equal on both sides). F1 and F2 must not be the same form.

Reaction: `[F1_name] → [F2_name]`

### Forbidden

- A metabolite cannot appear on both LHS and RHS of the same reaction (the existing constructor enforces this via the at-most-1-metabolite-per-side rule, plus the implicit constraint that the enzyme form changes).
- Self-loops (same enzyme form on both sides) are not valid elementary reactions.

---

## Catalytic Topology Enumeration

### Definitions

- **Core forms**: EnzymeFormSpecs with no regulator sites occupied.
- **Free enzyme**: the core form with all sites unoccupied (e.g., `E_0_0_0_0`).
- **Reaction graph**: nodes = core forms, edges = all valid elementary reactions between core forms (binding, release, isomerization per rules above).

### Valid Catalytic Topology

A catalytic topology is a connected subgraph of the core reaction graph such that:

1. **Contains the free enzyme form.**
2. **All forms are reachable** from the free enzyme via the included edges.
3. **Per-cycle stoichiometry = 1x**: Every simple cycle through the free enzyme has net stoichiometry equal to exactly 1× the reaction (each substrate consumed once, each product produced once), or 0× (futile exchange cycle in branched mechanisms).
4. **At least one 1x cycle exists**: The mechanism must have at least one cycle that catalyzes the net reaction.
5. **Form count ≤ max_forms** (before adding dead-ends).
6. **No isolated forms**: Every form participates in at least one reaction.

### Futile Cycles

Futile cycles with 0× stoichiometry (e.g., E→EA→EAB→EB→E in random-order Bi-Bi where net = bind A, bind B, release A, release B = 0) are allowed. They are natural consequences of branched mechanisms and represent thermodynamic equilibration between alternative binding orders.

### Enumeration Algorithm

The exact algorithm is implementation-defined but must:
- Start from the free enzyme and build outward
- Prune branches that cannot lead to valid stoichiometry
- Respect the `max_forms` limit
- Not create duplicate topologies

---

## Dead-End Regulator Attachments

### Definitions

- **Dead-end form**: an EnzymeFormSpec that has the same core substrate/product site states as some catalytic-cycle form, plus one or more regulator sites occupied.
- **Dead-end chain**: A sequence of dead-end forms where each successive form has one additional regulator site occupied. Chains are allowed when separate binding sites exist for the regulator (per `max_sites` from the EnzymeReaction). E.g., if G6P has 2 regulator sites, then E→E_G6P→E_G6P_G6P is a valid chain.

### Dead-End Options Per Cycle Form

For each catalytic-cycle form F, determine all EnzymeFormSpecs that:
- Match F's core substrate/product site states exactly
- Have one or more regulator sites additionally occupied
- Exist in the `enumerate_enzyme_forms` output

These dead-end extensions form a lattice ordered by set inclusion of occupied regulator sites.

A valid dead-end configuration for form F is any **downward-closed subset** of this lattice. "Downward-closed" means: if a form with regulators {R1, R2} is included, then forms with {R1} and {R2} must also be included (because binding is one-at-a-time elementary).

### Factoring

Dead-end options are independent across cycle forms. The total dead-end multiplier for a topology is the product of per-cycle-form option counts:

```
dead_end_options = product(n_dead_end_configs(F) for F in cycle_forms)
```

Each dead-end form also introduces binding/release reactions connecting it to its parent form.

---

## RE/SS Assignments

For a mechanism with N steps, each step can be rapid-equilibrium (RE, `true`) or steady-state (SS, `false`), with the constraint that **at least one step must be SS**.

Valid assignments: 2^N − 1 per topology (excluding all-RE).

This dimension is independent of topology and dead-ends.

---

## Symmetry Constraint Detection

Two steps in a mechanism are **structurally symmetric** if they represent the same metabolite binding to/releasing from equivalent sites. Specifically:

- Both steps involve the same metabolite
- Both are binding reactions (or both release reactions)
- The binding sites have the same metabolite species and differ only by site index (e.g., site 1 vs site 2 of the same metabolite's extra sites)

For each set of symmetric steps, generate two variants:
1. **Unconstrained**: all parameters independent
2. **Constrained**: equilibrium/rate constants are equal (e.g., `K_i = K_j` for RE steps, or `k_if = k_jf, k_ir = k_jr` for SS steps)

If there are M independent symmetry groups, this adds up to 2^M variants per mechanism.

---

## Size Limits

- **Default `max_forms`** = 2 × (number of distinct binding sites in EnzymeFormSpec). For the glucokinase example with 7 sites (Glu₁, ATP₁, G6P₁, ADP₁, Phosphate₁, G6P_reg₁, G6P_reg₂), this gives max_forms = 14.
- This cap applies to the total number of enzyme forms including dead-end forms.
- User can override with a different value.

---

## Compilation Considerations

### Key Principle

`enumerate_mechanisms` must NOT create `EnzymeMechanism` types during enumeration. All enumeration works with runtime data structures (`MechanismSpec`, `EnzymeFormSpec`, `SiteState`). Type-parameterized `EnzymeMechanism` instances are only created on explicit `EnzymeMechanism(spec)` conversion.

### Function Barriers

Following the existing pattern in `enumerate_enzyme_forms`:
- Use `@nospecialize` on reaction-specific arguments in hot enumeration paths
- Avoid `Iterators.product` splats that create 2^n specialized tuple types
- Use flat modular-arithmetic indexing for combinatorial enumeration

### Performance Targets

- Enumeration of the glucokinase reaction (7 sites, 3 regulators) should complete in < 30s cold, < 5s warm
- Individual `EnzymeMechanism(spec)` conversion should take < 1s
- Memory: the iterator state should be O(topologies), not O(total_mechanisms)

---

## Exported Symbols

Add to `EnzymeRates.jl`:
```julia
export enumerate_mechanisms, MechanismSpec
```

---

## Tests

### 1. Correctness: Uni-Uni

```julia
r = @enzyme_reaction begin
    substrates: S[C]
    products:   P[C]
end
mechs = collect(enumerate_mechanisms(r))
```

Expected mechanisms (catalytic topologies only, no regulators):
- **2-step**: E + S ⇌ ES, ES ⇌ E + P (ordered Uni-Uni, no isomerization)
- **3-step**: E + S ⇌ ES, ES ⇌ EP, EP ⇌ E + P (with isomerization intermediate)

Each topology × RE/SS assignments. Verify expected count.

### 2. Correctness: Bi-Bi

```julia
r = @enzyme_reaction begin
    substrates: A[C], B[N]
    products:   P[C], Q[N]
end
mechs = collect(enumerate_mechanisms(r))
```

Must include:
- Ordered Bi-Bi (sequential binding/release)
- Random-order Bi-Bi (branched binding)
- Verify that all enumerated mechanisms pass `EnzymeMechanism(spec)` constructor validation.

### 3. Correctness: With Regulators

```julia
r = @enzyme_reaction begin
    substrates: S[C]
    products:   P[C]
    regulators: I[N]
end
mechs = collect(enumerate_mechanisms(r))
```

Must include dead-end inhibitor complexes (EI, ESI, etc.) attached to cycle forms.

### 4. Validation: All Mechanisms Constructible

For every `MechanismSpec` from the Uni-Uni enumeration, verify:
```julia
m = EnzymeMechanism(spec)
# No errors thrown
# rate_equation(m) compiles and returns a function
```

### 5. Stoichiometry: Random-Order Is Valid

Verify that random-order Bi-Bi mechanisms satisfy per-cycle 1× stoichiometry even though total reaction sum > 1×.

### 6. Count Verification

For simple cases (Uni-Uni), manually count expected mechanisms and verify `length(iter)` matches.

### 7. Compilation Guard

```julia
rxn = @enzyme_reaction begin
    substrates: Glu[C6H12O6], ATP[C10H16N5O13P3]
    products: G6P[C6H13O9P], ADP[C10H15N5O10P2]
    regulators: Phosphate[PO4], G6P[C6H13O9P], G6P[C6H13O9P]
end
t1 = @elapsed begin
    iter = enumerate_mechanisms(rxn)
    count = length(iter)
end
t2 = @elapsed begin
    iter = enumerate_mechanisms(rxn)
    count = length(iter)
end
@test t1 < 30.0   # cold
@test t2 < 5.0    # warm
@test count < 10^6 # sanity check on mechanism count

# Conversion test: first mechanism should compile quickly
spec = first(iter)
t3 = @elapsed m = EnzymeMechanism(spec)
@test t3 < 2.0
```

### 8. Max Forms Limit

Verify that `max_forms` kwarg restricts mechanism size:
```julia
r = @enzyme_reaction begin
    substrates: S[C]
    products:   P[C]
end
small = collect(enumerate_mechanisms(r; max_forms=3))
large = collect(enumerate_mechanisms(r; max_forms=10))
@test length(small) <= length(large)
@test all(n_enzyme_forms(s) <= 3 for s in small)
```

### 9. Symmetry Constraints

For a reaction with equivalent extra binding sites:
```julia
r = @enzyme_reaction begin
    substrates: S[C, 2]
    products:   P[C]
end
mechs = collect(enumerate_mechanisms(r))
# Should include both constrained (K_bind_S1 = K_bind_S2)
# and unconstrained variants for mechanisms using both S sites
```

---

## File Changes

| File | Change |
|------|--------|
| `src/mechanism_enumeration.jl` | Add `MechanismSpec`, `MechanismIterator`, `enumerate_mechanisms`, reaction graph construction, topology enumeration, dead-end factoring, RE/SS enumeration, symmetry constraint detection |
| `src/types.jl` | Add `EnzymeMechanism(spec::MechanismSpec)` constructor overload |
| `src/EnzymeRates.jl` | Export `enumerate_mechanisms`, `MechanismSpec` |
| `test/test_mechanism_enumeration.jl` | Add all tests listed above |
| `.claude/CLAUDE.md` | Update source layout docs |
