# EnzymeRates.jl

A Julia package for deriving and evaluating steady-state enzyme rate equations via the quasi-steady-state approximation (QSSA).

## Features

- Define enzyme mechanisms using a concise DSL or programmatic API
- Validate mechanisms via atomic conservation checking
- Derive QSSA rate equations compiled into zero-allocation numeric functions
- Count independent kinetic parameters (Haldane/Wegscheider constraints)
- Enumerate all valid mechanisms for a given reaction and parameter count

## Installation

```julia
using Pkg
Pkg.add(url="https://github.com/your-username/EnzymeRates.jl")
```

## Quick Start

```julia
using EnzymeRates

# Define a reversible Uni-Uni mechanism
m = @mechanism begin
    [E, S(C=1)] --> [ES]
    [ES] --> [E, P(C=1)]
end

validate(m)  # true

# Compute rate: rate_equation(m, params, concs) -> Float64
v = rate_equation(m, (k1f=3.2, k1r=0.8, k2f=2.5, k2r=1.1), (S=0.7, P=0.3))

# Human-readable equation
rate_equation_string(m)
# (E_total * (k1f*S*k2f - k1r*k2r*P)) / ((k1r + k2f) + (k1f*S + k2r*P))
```

## Defining Mechanisms

### `@mechanism` macro

Each line is a reversible step `[reactants] --> [products]`. Bare symbols are enzyme forms; symbols with keyword arguments are metabolites with atomic compositions:

```julia
m = @mechanism begin
    [E, S(C=6, H=12, O=6)] --> [ES]
    [ES] --> [E, P(C=6, H=12, O=6)]
end
```

### `@enzyme_reaction` macro

Defines the overall reaction for use with `enumerate_mechanisms`:

```julia
spec = @enzyme_reaction begin
    substrates: S(C=6), ATP(C=10, N=5, P=3)
    products:   G6P(C=6, P=1), ADP(C=10, N=5, P=2)
    regulators: I(C=5)  # optional
end
```

### Programmatic construction

```julia
E  = Species(:E, enzyme)
ES = Species(:ES, enzyme)
S  = Species(:S, metabolite, Dict(:C => 1))
P  = Species(:P, metabolite, Dict(:C => 1))

m = EnzymeMechanism([[E, S] => [ES], [ES] => [E, P]])
```

## Rate Equations

`rate_equation(m, params, concs)` computes the steady-state rate, returning `Float64`:

- `params`: `NamedTuple` of rate constants (`k1f`, `k1r`, `k2f`, `k2r`, ...) and optionally `E_total` (defaults to `1.0`)
- `concs`: `NamedTuple` of metabolite concentrations (`S`, `P`, ...)

Rate constants are named by step index. Each step has a forward (`kNf`) and reverse (`kNr`) constant. Binding steps yield pseudo-first-order rates: `kNf * [metabolite]`.

## Querying a Mechanism

```julia
enzyme_forms(m)          # distinct enzyme states
metabolites(m)           # distinct metabolites
n_states(m)              # number of enzyme states
graph(m)                 # (SimpleDiGraph, Vector{Species})
stoich_matrix(m)         # metabolites x steps matrix
n_independent_params(m)  # independent params after thermodynamic constraints
validate(m)              # atomic conservation check
```

## Mechanism Enumeration

```julia
spec = @enzyme_reaction begin
    substrates: S(C=1)
    products:   P(C=1)
end

mechanisms = enumerate_mechanisms(spec, 3)  # all valid mechanisms with 3 indep. params
```

Each returned mechanism is guaranteed to contain free enzyme `E`, have a connected graph, conserve atoms, match the reaction spec, and have the requested parameter count.

## API Reference

### Types

| Type | Description |
|------|-------------|
| `Species(name, role[, atoms])` | Chemical species. `role` is `enzyme` or `metabolite`. |
| `ReactionSpec(substrates, products[, regulators])` | Overall reaction specification. |
| `EnzymeMechanism(steps)` | Mechanism: vector of `Pair{Vector{Species}, Vector{Species}}`. |

### Functions

| Function | Description |
|----------|-------------|
| `rate_equation(m, params, concs)` | Compiled QSSA rate equation → `Float64`. Zero allocations. |
| `rate_equation_string(m)` | Human-readable rate equation string. |
| `enzyme_forms(m)` | Distinct enzyme states. |
| `metabolites(m)` | Distinct metabolites. |
| `n_states(m)` | Number of enzyme states. |
| `graph(m)` | Enzyme-form connectivity graph. |
| `stoich_matrix(m)` | Stoichiometry matrix (metabolites x steps). |
| `n_independent_params(m)` | Independent parameters after thermodynamic constraints. |
| `validate(m)` | Check atomic conservation. |
| `enumerate_mechanisms(spec, n)` | All valid mechanisms with `n` independent parameters. |

## Running Tests

```julia
] test EnzymeRates
```
