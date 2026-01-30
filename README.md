# EnzymeRates.jl

A Julia package for building, validating, and analyzing enzyme kinetic mechanisms using the King-Altman method.

EnzymeRates lets you:

- **Define** enzyme mechanisms using a concise DSL or programmatic API
- **Validate** mechanisms via automatic atomic conservation checking
- **Analyze** mechanism structure: enzyme forms, metabolites, stoichiometry, connectivity
- **Count** independent kinetic parameters after Haldane/Wegscheider constraints
- **Derive** steady-state rate equations using King-Altman spanning-tree enumeration
- **Compile** rate equations into fast numeric functions
- **Enumerate** all valid mechanisms for a given reaction and parameter count

## Installation

```julia
using Pkg
Pkg.add(url="https://github.com/your-username/EnzymeRates.jl")
```

Or for local development:

```julia
] dev path/to/EnzymeRates
```

## Quick Start

```julia
using EnzymeRates

# Define a simple reversible Uni-Uni mechanism
m = @mechanism begin
    [E, S(C=1)] --> [ES]
    [ES] --> [E, P(C=1)]
end

# Validate atomic conservation
validate(m)  # true

# Get the compiled rate function
fn = rate_function(m)

# Evaluate: fn(rate_constants, concentrations) -> velocity
v = fn((k1f=3.2, k1r=0.8, k2f=2.5, k2r=1.1), (S=0.7, P=0.3))
# v ≈ 0.909
```

## Core Concepts

### Species

Every participant in a mechanism is a `Species` with a name, a role, and an optional atomic composition:

```julia
# Enzyme forms — no atoms (their atomic content is inferred from bound metabolites)
E  = Species(:E, enzyme)
ES = Species(:ES, enzyme)

# Metabolites — explicit atomic composition
S = Species(:S, metabolite, Dict(:C => 6, :H => 12, :O => 6))
P = Species(:P, metabolite, Dict(:C => 6, :H => 12, :O => 6))
```

The `role` field is a `SpeciesRole` enum with values `enzyme` and `metabolite`. Enzyme forms participate in the catalytic cycle and are conserved (their total concentration is constant). Metabolites are substrates, products, or regulators whose concentrations are independent variables.

### EnzymeMechanism

A mechanism is a list of elementary reversible steps. Each step is a `Pair` mapping a vector of reactant species to a vector of product species:

```julia
m = EnzymeMechanism([
    [E, S] => [ES],       # E + S ⇌ ES  (substrate binding)
    [ES]   => [E, P]      # ES ⇌ E + P  (catalysis + product release)
])
```

Every step must contain exactly one enzyme form on each side. Metabolites appear alongside the enzyme form when they bind to or release from the enzyme.

### Rate Constants

Each step automatically receives a forward and reverse rate constant named by its position:

| Step index | Forward | Reverse |
|:----------:|:-------:|:-------:|
| 1          | `k1f`   | `k1r`   |
| 2          | `k2f`   | `k2r`   |
| 3          | `k3f`   | `k3r`   |
| ...        | ...     | ...     |

Binding steps produce pseudo-first-order rate constants: the edge rate in the King-Altman graph is `k_if * [metabolite]`. Unimolecular steps (isomerization, product release) have edge rates equal to the rate constant alone.

## Defining Mechanisms

### The `@mechanism` Macro

The most direct way to define a mechanism. Each line is a step written as `[reactants] --> [products]`. Bare symbols are enzyme forms; symbols with parenthesized keyword arguments are metabolites with atomic compositions:

```julia
m = @mechanism begin
    [E, S(C=6, H=12, O=6)] --> [ES]
    [ES] --> [E, P(C=6, H=12, O=6)]
end
```

Atoms are specified as keyword arguments: `S(C=6, H=12, O=6)` means the metabolite S contains 6 carbon, 12 hydrogen, and 6 oxygen atoms. These are used by `validate()` to check atomic conservation at each step.

### The `@enzyme_reaction` Macro

Defines the overall reaction (substrates, products, regulators) without specifying a mechanism. This is used as input to `enumerate_mechanisms`:

```julia
spec = @enzyme_reaction begin
    substrates: S(C=6, H=12, O=6), ATP(C=10, H=16, N=5, O=13, P=3)
    products:   G6P(C=6, H=13, O=9, P=1), ADP(C=10, H=15, N=5, O=10, P=2)
    regulators: I(C=5, H=8, N=2)
end
```

Each line starts with a label (`substrates:`, `products:`, or `regulators:`) followed by one or more species with their atomic compositions. The `regulators` line is optional.

### Programmatic Construction

You can also build species and mechanisms directly without macros:

```julia
E   = Species(:E, enzyme)
ES  = Species(:ES, enzyme)
S   = Species(:S, metabolite, Dict(:C => 1))
P   = Species(:P, metabolite, Dict(:C => 1))

m = EnzymeMechanism([
    [E, S] => [ES],
    [ES]   => [E, P]
])
```

This is useful when generating mechanisms dynamically or from external data.

## Querying a Mechanism

### Structural Accessors

```julia
m = @mechanism begin
    [E, S(C=1, H=1)] --> [ES]
    [ES] --> [EPQ]
    [EPQ] --> [EQ, P(C=1)]
    [EQ] --> [E, Q(H=1)]
end

enzyme_forms(m)   # [E, ES, EPQ, EQ]  — distinct enzyme states
metabolites(m)    # [S, P, Q]         — distinct metabolites
n_states(m)       # 4                 — number of enzyme states
```

### Connectivity Graph

```julia
g, forms = graph(m)
# g is a Graphs.SimpleDiGraph where nodes correspond to enzyme forms.
# forms[i] gives the Species for node i.
# Edges represent reversible steps (both directions added).
```

### Stoichiometry Matrix

```julia
S = stoich_matrix(m)
# Rows = metabolites (same order as metabolites(m))
# Columns = steps (same order as m.steps)
# S[i,j] > 0 means metabolite i is produced in step j
# S[i,j] < 0 means metabolite i is consumed in step j
```

### Parameter Grouping

```julia
groups = param_groups(m)
# Vector of vectors of step indices.
# Steps involving the same metabolite are grouped together by default.
```

### Independent Parameters

```julia
n_independent_params(m)
# Number of independently adjustable kinetic parameters after accounting
# for Haldane relation and Wegscheider conditions.
#
# Formula: 2 * n_steps - n_constraints
# where n_constraints counts independent thermodynamic cycles.
```

For a simple 2-step Uni-Uni mechanism (4 raw rate constants, 1 Haldane constraint): 4 - 1 = 3 independent parameters.

### Validation

```julia
validate(m)  # returns true or false
```

Checks atomic conservation at every step. The algorithm:

1. Assigns atomic content to the free enzyme form (empty by convention).
2. Propagates via BFS: for each step, the destination enzyme form's atoms are computed from the source form's atoms plus consumed metabolite atoms minus produced metabolite atoms.
3. Verifies that total atoms on each side of every step match.

This catches errors like mismatched atom counts between substrates and products:

```julia
m_bad = @mechanism begin
    [E, S(C=2)] --> [ES]
    [ES] --> [E, P(C=1)]
end
validate(m_bad)  # false — carbon not conserved
```

## King-Altman Rate Equations

### Compiled Rate Function

```julia
fn = rate_function(m)
```

Returns a closure with signature `fn(params, concs) -> Float64`:

- `params` — a `NamedTuple` of rate constants: `(k1f=..., k1r=..., k2f=..., k2r=..., ...)`
- `concs` — a `NamedTuple` of metabolite concentrations: `(S=..., P=..., ...)`
- Optionally include `E_total` in `concs` to scale by total enzyme concentration (defaults to `1.0`)

The rate is computed using the King-Altman method:

```
v = E_total * (∏ forward_rates - ∏ reverse_rates) / Σ spanning_tree_products
```

where the denominator sums over all directed spanning trees of the enzyme-form graph, each weighted by the product of its edge rates.

### Rate Equation String

```julia
rate_equation_string(m)
```

Returns a human-readable string representation of the full rate equation:

```julia
m = @mechanism begin
    [E, S(C=1)] --> [ES]
    [ES] --> [E, P(C=1)]
end

println(rate_equation_string(m))
# (E_total * (k1f*S*k2f - k1r*k2r*P)) / ((k1r + k2f) + (k1f*S + k2r*P))
```

Denominator terms with combined edges (multiple rate constants contributing to the same directed edge) are shown as sums in parentheses, e.g. `(k1f*S + k2r*P)`.

## Mechanism Enumeration

```julia
mechanisms = enumerate_mechanisms(spec, n_params)
```

Automatically generates all valid mechanisms for a given overall reaction that have exactly `n_params` independent kinetic parameters.

### Algorithm

1. **Generate enzyme states**: free enzyme `E`, plus complexes formed by binding subsets of metabolites (`ES`, `EP`, `ESP`, etc.). For reactions with atom transfer between substrates and products, modified enzyme forms (e.g., `F` for ping-pong) are also generated.

2. **Generate candidate steps**: all chemically plausible elementary steps between pairs of enzyme forms — binding, unbinding, catalytic release (bound metabolite converted to a different product), catalytic exchange, and internal isomerization.

3. **Enumerate subsets**: try all combinations of candidate steps (up to 8 steps) and keep those that form a valid mechanism.

4. **Validation filters** (applied to each candidate):
   - Contains the free enzyme form `E`
   - Enzyme-form graph is connected
   - Net reaction matches the specification (substrates consumed, products formed)
   - Atomic conservation holds at every step
   - Independent parameter count matches target

### Example

```julia
spec = @enzyme_reaction begin
    substrates: S(C=1)
    products:   P(C=1)
end

mechanisms = enumerate_mechanisms(spec, 3)
# Returns all valid Uni-Uni mechanisms with 3 independent parameters.
# The simplest is [E,S]=>[ES], [ES]=>[E,P] (the classic Michaelis-Menten mechanism).

for m in mechanisms
    println(m.steps)
    println("  states: ", n_states(m))
    println("  valid:  ", validate(m))
    println()
end
```

## Tutorials

### Tutorial 1: Michaelis-Menten (Uni-Uni Reversible)

The simplest enzyme mechanism: one substrate, one product, two elementary steps.

```julia
using EnzymeRates

# Define the mechanism
m = @mechanism begin
    [E, S(C=1)] --> [ES]     # substrate binding
    [ES] --> [E, P(C=1)]     # catalysis + product release
end

# Inspect structure
println("Enzyme states: ", [s.name for s in enzyme_forms(m)])  # [:E, :ES]
println("Metabolites:   ", [s.name for s in metabolites(m)])   # [:S, :P]
println("States:        ", n_states(m))                        # 2
println("Indep params:  ", n_independent_params(m))            # 3

# Validate atomic conservation
@assert validate(m)

# Define kinetic parameters and concentrations
params = (k1f=3.2, k1r=0.8, k2f=2.5, k2r=1.1)
concs  = (S=0.7, P=0.3)

# Compute rate
fn = rate_function(m)
v = fn(params, concs)
println("Rate: ", round(v, digits=4))  # ≈ 0.9091

# Verify Haldane relation: at equilibrium, v = 0
# Keq = k1f*k2f / (k1r*k2r) = P_eq / S_eq
Keq = params.k1f * params.k2f / (params.k1r * params.k2r)
println("Keq = ", round(Keq, digits=4))  # 9.0909

S_eq = 1.0
P_eq = Keq * S_eq  # equilibrium product concentration
v_eq = fn(params, (S=S_eq, P=P_eq))
println("Rate at equilibrium: ", v_eq)  # ≈ 0 (within floating-point precision)

# Print the rate equation
println("\nRate equation:")
println(rate_equation_string(m))
```

### Tutorial 2: Ordered Bi-Bi Sequential Mechanism

Two substrates bind in order, catalysis occurs, then two products release in order.

```julia
using EnzymeRates

E   = Species(:E,   enzyme)
EA  = Species(:EA,  enzyme)
EAB = Species(:EAB, enzyme)
EPQ = Species(:EPQ, enzyme)
EQ  = Species(:EQ,  enzyme)

A = Species(:A, metabolite, Dict(:C => 2))
B = Species(:B, metabolite, Dict(:C => 3))
P = Species(:P, metabolite, Dict(:C => 2))
Q = Species(:Q, metabolite, Dict(:C => 3))

m = EnzymeMechanism([
    [E, A]   => [EA],     # first substrate binds
    [EA, B]  => [EAB],    # second substrate binds
    [EAB]    => [EPQ],    # catalytic interconversion
    [EPQ]    => [EQ, P],  # first product released
    [EQ]     => [E, Q]    # second product released
])

# Check the structure
println("States: ", n_states(m))   # 5
println("Valid:  ", validate(m))   # true
println("Independent parameters: ", n_independent_params(m))

# The stoichiometry matrix shows metabolite changes per step
S = stoich_matrix(m)
mets = metabolites(m)
println("\nStoichiometry matrix (rows = metabolites, cols = steps):")
for (i, met) in enumerate(mets)
    println("  ", met.name, ": ", S[i, :])
end

# Compute a rate
params = (k1f=2.0, k1r=0.5, k2f=3.0, k2r=0.4, k3f=1.5, k3r=0.3,
          k4f=4.0, k4r=0.6, k5f=2.5, k5r=0.7)
concs = (A=1.0, B=1.0, P=0.1, Q=0.1)
fn = rate_function(m)
v = fn(params, concs)
println("\nRate: ", round(v, digits=4))

# The rate equation string for a 5-state mechanism is large
println("\nRate equation:")
println(rate_equation_string(m))
```

### Tutorial 3: Ping-Pong Bi-Bi Mechanism

In a ping-pong mechanism, the first substrate binds and transfers a chemical group to the enzyme (forming a modified enzyme F), then the first product leaves. The second substrate binds to F, receives the group, and the second product leaves, regenerating free enzyme E.

```julia
using EnzymeRates

# A transfers an amino group (N) to enzyme; enzyme transfers it to B
m = @mechanism begin
    [E,  A(C=2, N=1)] --> [EA]   # A binds
    [EA]               --> [FP]   # amino group transferred to enzyme
    [FP]               --> [F, P(C=2)]  # first product (deaminated) released
    [F,  B(C=3)]       --> [FB]   # second substrate binds modified enzyme
    [FB]               --> [EQ]   # amino group transferred to product
    [EQ]               --> [E, Q(C=3, N=1)]  # second product (aminated) released
end

println("States: ", n_states(m))    # 6
println("Valid:  ", validate(m))    # true — nitrogen is conserved through F

# The modified enzyme form F carries the nitrogen atom.
# validate() infers this automatically: E carries nothing,
# after A(C=2,N=1) binds and P(C=2) leaves, the enzyme retains N=1.

# Enzyme-form connectivity
g, forms = graph(m)
println("\nEnzyme forms: ", [f.name for f in forms])

# Compute rate with total enzyme included in concentrations
params = (k1f=2.0, k1r=0.5, k2f=3.0, k2r=0.4, k3f=1.5, k3r=0.3,
          k4f=2.5, k4r=0.6, k5f=1.8, k5r=0.2, k6f=3.5, k6r=0.7)
concs = (A=1.0, P=0.1, B=1.0, Q=0.1, E_total=0.001)

fn = rate_function(m)
v = fn(params, concs)
println("\nRate (with E_total=0.001): ", v)
```

### Tutorial 4: Mechanism Enumeration

Automatically discover all valid mechanisms for a given overall reaction.

```julia
using EnzymeRates

# Define the overall reaction: S → P (uni-uni)
spec = @enzyme_reaction begin
    substrates: S(C=1)
    products:   P(C=1)
end

# Find all mechanisms with 3 independent kinetic parameters
mechanisms = enumerate_mechanisms(spec, 3)
println("Found $(length(mechanisms)) mechanisms with 3 independent parameters\n")

for (i, m) in enumerate(mechanisms)
    println("Mechanism $i:")
    for step in m.steps
        lhs_names = join([string(s.name) for s in step.first], " + ")
        rhs_names = join([string(s.name) for s in step.second], " + ")
        println("  $lhs_names ⇌ $rhs_names")
    end
    println("  States: ", n_states(m))
    println()
end

# Each discovered mechanism is guaranteed to:
# - include free enzyme E
# - have a connected enzyme-form graph
# - conserve atoms at every step
# - consume S and produce P (matching the reaction spec)
# - have exactly 3 independent kinetic parameters

# You can compute rate equations for any of them
fn = rate_function(mechanisms[1])
v = fn((k1f=1.0, k1r=0.1, k2f=5.0, k2r=0.5), (S=1.0, P=0.0))
println("Rate for mechanism 1: ", round(v, digits=4))
```

### Tutorial 5: Comparing Mechanisms

Use the rate function and Haldane relation to compare kinetic behavior across mechanisms.

```julia
using EnzymeRates

# Two-step Uni-Uni (Michaelis-Menten)
m1 = @mechanism begin
    [E, S(C=1)] --> [ES]
    [ES] --> [E, P(C=1)]
end

# Three-step Uni-Uni with an intermediate (E → ES → EP → E)
E  = Species(:E,  enzyme)
ES = Species(:ES, enzyme)
EP = Species(:EP, enzyme)
S  = Species(:S, metabolite, Dict(:C => 1))
P  = Species(:P, metabolite, Dict(:C => 1))

m2 = EnzymeMechanism([
    [E, S] => [ES],    # substrate binding
    [ES]   => [EP],    # isomerization
    [EP]   => [E, P]   # product release
])

println("Mechanism 1: ", n_states(m1), " states, ", n_independent_params(m1), " indep params")
println("Mechanism 2: ", n_states(m2), " states, ", n_independent_params(m2), " indep params")

# Both must satisfy the Haldane relation:
# Keq = [P]_eq / [S]_eq = ∏(k_forward) / ∏(k_reverse)

# Mechanism 1: Keq = k1f*k2f / (k1r*k2r)
p1 = (k1f=10.0, k1r=1.0, k2f=5.0, k2r=0.5)
Keq1 = p1.k1f * p1.k2f / (p1.k1r * p1.k2r)

# Mechanism 2: Keq = k1f*k2f*k3f / (k1r*k2r*k3r)
p2 = (k1f=10.0, k1r=1.0, k2f=8.0, k2r=0.8, k3f=5.0, k3r=0.5)
Keq2 = p2.k1f * p2.k2f * p2.k3f / (p2.k1r * p2.k2r * p2.k3r)

println("\nKeq (mechanism 1): ", Keq1)  # 100.0
println("Keq (mechanism 2): ", Keq2)    # 1000.0

# Compare rates at the same substrate/product concentrations
fn1 = rate_function(m1)
fn2 = rate_function(m2)

for S_val in [0.01, 0.1, 1.0, 10.0]
    v1 = fn1(p1, (S=S_val, P=0.0))
    v2 = fn2(p2, (S=S_val, P=0.0))
    println("S=$S_val:  v1=$(round(v1, digits=4)),  v2=$(round(v2, digits=4))")
end
```

## API Reference

### Types

| Type | Description |
|------|-------------|
| `Species(name, role[, atoms])` | A chemical species. `role` is `enzyme` or `metabolite`. `atoms` is a `Dict{Symbol,Int}`. |
| `ReactionSpec(substrates, products[, regulators])` | Overall reaction specification. |
| `EnzymeMechanism(steps)` | A mechanism: vector of `Pair{Vector{Species}, Vector{Species}}`. |
| `SpeciesRole` | Enum: `enzyme`, `metabolite`. |

### Macros

| Macro | Returns | Description |
|-------|---------|-------------|
| `@mechanism begin ... end` | `EnzymeMechanism` | Define a mechanism with `[lhs] --> [rhs]` steps. |
| `@enzyme_reaction begin ... end` | `ReactionSpec` | Define overall reaction with `substrates:`, `products:`, `regulators:` labels. |

### Functions

| Function | Returns | Description |
|----------|---------|-------------|
| `enzyme_forms(m)` | `Vector{Species}` | Distinct enzyme states in the mechanism. |
| `metabolites(m)` | `Vector{Species}` | Distinct metabolites in the mechanism. |
| `n_states(m)` | `Int` | Number of enzyme states. |
| `graph(m)` | `(SimpleDiGraph, Vector{Species})` | Enzyme-form connectivity graph and node labels. |
| `stoich_matrix(m)` | `Matrix{Int}` | Stoichiometry matrix (metabolites x steps). |
| `param_groups(m)` | `Vector{Vector{Int}}` | Default parameter grouping by metabolite identity. |
| `n_independent_params(m)` | `Int` | Independent kinetic parameters after thermodynamic constraints. |
| `validate(m)` | `Bool` | Check atomic conservation at every step. |
| `rate_function(m)` | `Function` | Compiled King-Altman rate equation: `(params, concs) -> Float64`. |
| `rate_equation_string(m)` | `String` | Human-readable rate equation formula. |
| `enumerate_mechanisms(spec, n)` | `Vector{EnzymeMechanism}` | All valid mechanisms with `n` independent parameters. |

## Dependencies

- [Graphs.jl](https://github.com/JuliaGraphs/Graphs.jl) — graph connectivity, spanning trees
- `LinearAlgebra` (stdlib) — used in King-Altman computations

## Running Tests

```julia
] test EnzymeRates
```

Tests cover four reference mechanisms (Uni-Uni, Uni-Bi ordered, Ping-Pong Bi-Bi, Sequential ordered Bi-Bi), each verified against an independent Laplacian cofactor reference implementation with random parameter sets, numeric spot checks, and Haldane equilibrium tests.
