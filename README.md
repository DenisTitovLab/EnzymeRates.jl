# EnzymeRates.jl

A Julia package for deriving and evaluating steady-state enzyme rate equations via the quasi-steady-state approximation (QSSA).

## Features

- Define enzyme mechanisms using a concise DSL or programmatic API
- Validate mechanisms via atomic conservation checking
- Derive QSSA rate equations compiled into zero-allocation numeric functions
- Automatic Haldane/Wegscheider constraint detection and dependent parameter elimination
- Human-readable rate equation output

## Installation

```julia
using Pkg
Pkg.add(url="https://github.com/DenisTitovLab/EnzymeRates.jl")
```

## Quick Start

```julia
using EnzymeRates

# Define a reversible Uni-Uni mechanism
m = @mechanism begin
    species: begin
        substrates: S[C]
        products:   P[C]
        enzymes:    E, ES[C]
    end
    steps: begin
        [E, S] <--> [ES]
        [ES] <--> [E, P]
    end
end

# Independent parameters + Keq + E_total (dependent k's are computed internally)
params = (k1f=3.2, k2f=2.5, k2r=1.1, Keq=5.0, E_total=1.0)
concs = (S=0.7, P=0.3)

# Compute rate: zero allocations, compiled at first call
v = rate_equation(m, params, concs)

# Human-readable equation
rate_equation_string(m)
```

## Defining Mechanisms

### `@mechanism` macro

`@mechanism` requires explicit species and steps blocks. Species include substrates, products, regulators, and enzyme forms. Atoms use chemical formula bracket syntax (`S[C6H12O6]`); bare symbols are allowed when all metabolites omit atoms.

Steps use the `<-->` arrow between bracketed sides:

```julia
m = @mechanism begin
    species: begin
        substrates: S[C6H12O6]
        products:   P[C6H12O6]
        enzymes:    E, ES[C6H12O6]
    end
    steps: begin
        [E, S] <--> [ES]
        [ES] <--> [E, P]
    end
end
```

### `@enzyme_reaction` macro

Define an overall reaction (substrates, products, optional regulators) without specifying the mechanism steps:

```julia
rxn = @enzyme_reaction begin
    substrates: A[C6H12O6], B[C10H16N5O13P3]
    products:   P[C6H12O6], Q[C10H16N5O13P3]
end
```

### Programmatic construction

```julia
species = (
    ((:S, ((:C, 1),)),),       # substrates
    ((:P, ((:C, 1),)),),       # products
    (),                        # regulators
    ((:E, ()), (:ES, ((:C, 1),))),  # enzyme forms
)
reactions = (
    ((:E, :S), (:ES,)),
    ((:ES,), (:E, :P)),
)

m = EnzymeMechanism(species, reactions)
```

## Rate Equations

`rate_equation(m, params, concs)` computes the steady-state rate (net consumption of the
first substrate, normalized by its stoichiometric coefficient):

- `params`: `NamedTuple` containing independent rate constants (`k1f`, `k1r`, ...), `Keq` (equilibrium constant), and `E_total` (total enzyme concentration). Dependent rate constants (determined by Haldane/Wegscheider constraints) are computed internally and should not be included.
- `concs`: `NamedTuple` of metabolite concentrations (`S`, `P`, ...)

Use `independent_parameters(m)` and `dependent_parameters(m)` to inspect which rate constants are independent vs. derived from constraints.

## Querying a Mechanism

Mechanisms are validated at construction (elementary-step structure, atomic
conservation, regulator balance). There is no separate `validate` API.

```julia
substrates(m)              # substrates with stoichiometric multiplicity
products(m)                # products with stoichiometric multiplicity
regulators(m)              # regulators
enzyme_forms(m)            # distinct enzyme states
metabolites(m)             # distinct metabolites
reactions(m)               # reaction steps as tuples of (lhs, rhs)
n_states(m)                # number of enzyme states
n_steps(m)                 # number of mechanism steps
graph(m)                   # (SimpleDiGraph, enzyme_forms)
stoich_matrix(m)           # metabolites × steps matrix
parameters(m)              # all rate constant names (k1f, k1r, ...)
independent_parameters(m)  # independent rate constant names
dependent_parameters(m)    # dependent params as (symbol, expression_string) pairs
```

## API Reference

### Types

| Type | Description |
|------|-------------|
| `EnzymeReaction{S,P,R}` | Overall reaction specification (substrates, products, regulators encoded in type parameters). |
| `EnzymeMechanism{Species,Reactions}` | Full mechanism with species and elementary steps encoded in type parameters. |

### Macros

| Macro | Description |
|-------|-------------|
| `@enzyme_reaction` | Create an `EnzymeReaction` from a DSL block. |
| `@mechanism` | Create an `EnzymeMechanism` from species + steps DSL blocks. |

### Functions

| Function | Description |
|----------|-------------|
| `rate_equation(m, params, concs)` | Compiled QSSA rate equation. Zero allocations. |
| `rate_equation_string(m)` | Human-readable rate equation string. |
| `substrates(m)` | Substrates (with stoichiometric multiplicity). |
| `products(m)` | Products (with stoichiometric multiplicity). |
| `regulators(m)` | Regulators. |
| `enzyme_forms(m)` | Distinct enzyme states. |
| `metabolites(m)` | Distinct metabolites. |
| `reactions(m)` | Reaction steps as `(lhs, rhs)` tuples. |
| `n_states(m)` | Number of enzyme states. |
| `n_steps(m)` | Number of mechanism steps. |
| `graph(m)` | Enzyme-form connectivity graph. |
| `stoich_matrix(m)` | Stoichiometry matrix (metabolites × steps). |
| `parameters(m)` | All rate constant names. |
| `all_parameters(m)` | Same as `parameters`. |
| `independent_parameters(m)` | Independent rate constant names (excludes dependent k's, Keq, E_total). |
| `dependent_parameters(m)` | Dependent parameters as `(symbol, expression_string)` pairs. |

## Running Tests

```julia
] test EnzymeRates
```
