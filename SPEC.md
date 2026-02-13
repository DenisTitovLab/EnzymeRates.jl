# SPEC: Mechanism Enumeration from EnzymeFormSpecs

## Overview

Add `enumerate_mechanisms(reaction::EnzymeReaction)` that generates all valid `EnzymeMechanism` configurations for a given reaction, using the enzyme forms from `enumerate_enzyme_forms`. Returns a lazy iterator over lightweight `MechanismSpec` structs (not type-parameterized `EnzymeMechanism` instances) to avoid compiling millions of types. Users convert selected candidates to `EnzymeMechanism` on demand via `EnzymeMechanism(spec)`.

The enumeration is factored into independent combinatorial dimensions:
1. **Catalytic topology** — the graph structure of the catalytic cycle, including substrate/product binding, release, isomerization steps, and optionally regulator binding/release (for activator mechanisms)
2. **Dead-end attachments** — which non-cycle forms (substrate, product, or regulator-bound) attach as dead-ends to cycle forms
3. **RE/SS assignments** — which steps are rapid-equilibrium vs steady-state
4. **Equivalent step constraints** — constrained (K_i = K_j) and unconstrained variants for equivalent steps

The iterator lazily computes the Cartesian product across dimensions without materializing all combinations.

---

## API

### Function Signature

```julia
enumerate_mechanisms(reaction::EnzymeReaction;
                     max_forms::Int = 3 * n_sites(reaction)
                    ) -> MechanismIterator
```

- `reaction`: the `EnzymeReaction` defining substrates, products, regulators with atoms and max binding sites.
- `max_forms`: maximum number of enzyme forms (including dead-end forms) in any single mechanism. Default = 3 × number of distinct binding sites in the `EnzymeFormSpec` (i.e., total sites from core substrates + core products + extra substrates + extra products + regulators). User can override.
- Returns a `MechanismIterator` that yields `MechanismSpec` structs.

### Helper

```julia
n_sites(reaction::EnzymeReaction)  # total binding sites (core + extra + regulator)
```

### MechanismSpec

Lightweight data struct (not type-parameterized):

```julia
struct MechanismSpec
    reaction::Any                    # source EnzymeReaction (untyped to avoid specialization)
    forms::Vector{Symbol}            # enzyme form names
    form_atoms::Vector{Vector{Pair{Symbol,Int}}}  # atoms for each enzyme form
    reactions::Vector{Tuple{Vector{Symbol}, Vector{Symbol}}}  # (lhs, rhs) per step
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

Reverse of binding. For standard binding/release (site goes empty↔occupied), the metabolite is the site's designated metabolite. For ping-pong partial release (site atoms decrease but site stays occupied), the metabolite is identified by the exact atom difference between the two site states, matched against known metabolite atoms.

Reaction: `[F1_name] → [F2_name, metabolite]`

### 3. Isomerization (catalytic conversion)

F1 must have ALL core substrate sites (index=1 for each substrate) occupied and ALL core product sites (index=1 for each product) unoccupied. F2 must have the reverse: ALL product sites occupied, ALL substrate sites unoccupied. (Or vice versa.) All non-core sites (extra substrate/product sites with index ≥ 2, and all regulator sites) must be identical between F1 and F2. The total atom content across all core catalytic sites must be conserved (equal on both sides). F1 and F2 must not be the same form.

Reaction: `[F1_name] → [F2_name]`

### Forbidden

- A metabolite cannot appear on both LHS and RHS of the same reaction (the existing constructor enforces this via the at-most-1-metabolite-per-side rule, plus the implicit constraint that the enzyme form changes).
- Self-loops (same enzyme form on both sides) are not valid elementary reactions.

---

## Catalytic Topology Enumeration

### Definitions

- **Free enzyme**: the EnzymeFormSpec with all sites unoccupied (e.g., `E_0_0_0_0`).
- **Reaction graph**: nodes = all EnzymeFormSpecs from `enumerate_enzyme_forms`, edges = all valid elementary reactions between forms (binding, release, isomerization per rules above). This includes forms with regulator sites occupied, enabling activator mechanisms where a regulator must be bound for catalysis.

### Valid Catalytic Topology

A catalytic topology is a connected subgraph of the reaction graph such that:

1. **Contains the free enzyme form.**
2. **All forms are reachable** from the free enzyme via the included edges.
3. **Per-cycle stoichiometry = 1x**: Every simple cycle through the free enzyme has net stoichiometry equal to exactly 1× the reaction (each substrate consumed once, each product produced once) with 0× for regulators (any regulator bound within the cycle must also be released), or 0× for all metabolites (futile exchange cycle in branched mechanisms).
4. **At least one 1x cycle exists**: The mechanism must have at least one cycle that catalyzes the net reaction.
5. **Form count ≤ max_forms** (before adding dead-ends).
6. **No isolated forms**: Every form participates in at least one reaction.

### Futile Cycles

Futile cycles with 0× stoichiometry (e.g., E→EA→EAB→EB→E in random-order Bi-Bi where net = bind A, bind B, release A, release B = 0) are allowed. They are natural consequences of branched mechanisms and represent thermodynamic equilibration between alternative binding orders. Cycles involving regulator binding/release (e.g., E→EA→EAS→EA→E for activator A) are valid as long as the regulator has 0 net stoichiometry and substrates/products have the correct 1× stoichiometry.

### Enumeration Algorithm

The exact algorithm is implementation-defined but must:
- Start from the free enzyme and build outward
- Prune branches that cannot lead to valid stoichiometry
- Respect the `max_forms` limit
- Not create duplicate topologies

---

## Dead-End Attachments

### Definitions

- **Dead-end form**: an EnzymeFormSpec that (a) is NOT part of the catalytic topology, (b) is connected to some cycle form via a valid elementary reaction (binding of one metabolite — substrate, product, or regulator), and (c) exists in the `enumerate_enzyme_forms` output.
- **Dead-end chain**: A sequence of dead-end forms where each successive form has one additional site occupied. Chains arise from multiple binding sites for the same metabolite (e.g., E→E_G6P→E_G6P_G6P for a regulator with 2 sites) or from successive binding of different metabolites (e.g., E→EP→EPQ for product inhibition).

Dead-end forms include:
- **Regulator-bound forms**: inhibitor or activator binding to cycle forms (e.g., EI, ESI)
- **Product-bound forms**: product inhibition (e.g., E + P ⇌ EP in a Bi-Bi mechanism where the cycle releases P from EQ, not E)
- **Abortive complexes**: forms with substrate/product combinations that cannot proceed catalytically (e.g., EAQ where substrate A and product Q are simultaneously bound)

Note: Forms that are part of the catalytic topology (e.g., EA in an essential activator mechanism E→EA→EAS→EA→E) are cycle forms, not dead-ends. Dead-ends are only those forms that branch off the cycle without participating in any cycle.

### Dead-End Options Per Cycle Form

For each form F in the catalytic topology, determine all EnzymeFormSpecs that:
- Differ from F by exactly one additional site occupied (any metabolite: substrate, product, or regulator)
- Are NOT already in the catalytic topology
- Exist in the `enumerate_enzyme_forms` output

These dead-end extensions form a lattice ordered by set inclusion of occupied sites.

A valid dead-end configuration for form F is any **downward-closed subset** of this lattice. "Downward-closed" means: if a form reachable via 2 binding steps is included, the intermediate form (reachable via 1 step) must also be included, because binding is one-at-a-time elementary.

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

## Equivalent Step Detection

Two steps in a mechanism are **equivalent** if they represent the same metabolite binding at the same site index but to different enzyme states. Specifically:

- Both steps are binding reactions
- Both bind the same metabolite
- Both bind at the same site index (e.g., both at site 1 of metabolite I)
- They differ only in the enzyme form they bind to (different enzyme state context)

For each set of equivalent steps, generate two variants:
1. **Unconstrained**: all parameters independent
2. **Constrained**: equilibrium/rate constants are equal (e.g., `K_i = K_j` for RE steps, or `k_if = k_jf, k_ir = k_jr` for SS steps)

If there are M independent equivalence groups, this adds up to 2^M variants per mechanism.

---

## Size Limits

- **Default `max_forms`** = 3 × (number of distinct binding sites in EnzymeFormSpec). For the glucokinase example with 7 sites (Glu₁, ATP₁, G6P₁, ADP₁, Phosphate₁, G6P_reg₁, G6P_reg₂), this gives max_forms = 21.
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
- **3-step**: E + S ⇌ ES, ES ⇌ EP, EP ⇌ E + P (with isomerization intermediate; simplest valid Uni-Uni)

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

### 3b. Correctness: Activator Mechanism

```julia
r = @enzyme_reaction begin
    substrates: S[C]
    products:   P[C]
    regulators: A[N]
end
mechs = collect(enumerate_mechanisms(r))
```

Must include activator mechanisms where A is part of the catalytic cycle, e.g.:
- E + A ⇌ EA, EA + S ⇌ EAS, EAS ⇌ EA + P (essential activator: A stays bound during catalysis, cycle = E→EA→EAS→EA→E)
- Verify that cycles through the free enzyme have 0× net stoichiometry for the regulator A (bound and released each cycle)

Must also include dead-end inhibitor mechanisms (same regulator as inhibitor: EI, ESI branching off cycle forms).

### 3c. Correctness: Product Inhibition / Abortive Complexes

```julia
r = @enzyme_reaction begin
    substrates: A[C], B[N]
    products:   P[C], Q[N]
end
mechs = collect(enumerate_mechanisms(r))
```

For an ordered Bi-Bi topology (E→EA→EAB→EPQ→EQ→E), must include mechanisms with dead-end forms such as:
- **Product inhibition**: EP (product P binding to free enzyme E, which normally only binds substrate A first)
- **Abortive complexes**: EAQ (substrate A + product Q simultaneously bound)
- Verify these dead-end forms are connected to cycle forms via valid elementary reactions

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
