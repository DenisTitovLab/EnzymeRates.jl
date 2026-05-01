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

## Define a mechanism, derive its rate equation, fit data

The running example is a uni-uni reaction `S ⇌ P` catalyzed by an MWC
homodimer in which substrate, product, the catalytic interconversion,
and an allosteric activator `A` all operate exclusively in the R
conformation. The T conformation is catalytically silent — a textbook
K-type allosteric activator.

```julia
using EnzymeRates

m = @allosteric_mechanism begin
    substrates: S
    products:   P
    allosteric_regulators: A::OnlyR

    site(:catalytic, 2): begin
        steps: begin
            [E, S] ⇌ [ES]      :: OnlyR
            [ES] <--> [EP]     :: OnlyR
            [EP] ⇌ [E, P]      :: OnlyR
        end
    end
end
```

The two `⇌` steps mark binding events at rapid equilibrium (one
binding constant `K` per step); the `<-->` step marks a steady-state
catalytic interconversion (independent forward and reverse rate
constants `kf`, `kr`). The `::OnlyR` annotation tells the framework
that those steps fire only in the R conformation, so the T-state
contribution to the rate equation is identically zero.

`parameters(m)` lists the names the framework needs at evaluation
time, `rate_equation_string(m)` prints the symbolic rate equation, and
`rate_equation(m, concs, params)` evaluates the rate numerically:

```julia
parameters(m)
rate_equation_string(m)
```
