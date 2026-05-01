# EnzymeRates.jl

Identify the best enzyme rate equation from kinetic data.

Given a reaction definition and experimental rate measurements at varying
substrate, product, and regulator concentrations, EnzymeRates enumerates
all biochemically valid mechanisms, fits each to the data, and selects
the simplest mechanism that adequately describes the data based on
leave-one-group-out cross-validation.

The package has first-class support for MWC allostery and for mechanisms
that mix steady-state and rapid-equilibrium elementary steps.
Thermodynamic constraints (Haldane, Wegscheider) are derived
automatically from the cycle structure of the mechanism, so users supply
only the independent rate constants plus a measured equilibrium
constant.

## Installation

```julia
# README-SKIP-IN-TEST
using Pkg
Pkg.add(url="https://github.com/DenisTitovLab/EnzymeRates.jl")
```

## Smoke test

```julia
using EnzymeRates
m = @enzyme_mechanism begin
    substrates: S
    products:   P
    steps: begin
        [E, S] <--> [ES]
        [ES] <--> [E, P]
    end
end
@assert :S in metabolites(m)
```
