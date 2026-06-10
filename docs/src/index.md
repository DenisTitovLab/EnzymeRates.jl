# EnzymeRates.jl

EnzymeRates.jl identifies the best enzyme rate equation from kinetic data.
Given a reaction definition and rate measurements at varying substrate and
product concentrations, the package enumerates the biochemically valid
mechanisms, fits each to the data, and selects the simplest equation that
describes the data.

## Installation

```julia
using Pkg
Pkg.add(url = "https://github.com/DenisTitovLab/EnzymeRates.jl")
```

## Where to go next

- [Getting Started](@ref) — an end-to-end example from reaction to identified
  rate equation.
- **Deriving rate equations** — how the package turns a mechanism into a
  symbolic rate law.
- **Fitting rate equations** — fitting a rate equation to kinetic data.
- **Identifying the best rate equation** — the model-selection search.
- [API Reference](@ref) — every exported name.
