# EnzymeRates.jl

[![Build Status](https://github.com/DenisTitovLab/EnzymeRates.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/DenisTitovLab/EnzymeRates.jl/actions/workflows/CI.yml?query=branch%3Amain)
[![Coverage](https://codecov.io/gh/DenisTitovLab/EnzymeRates.jl/branch/main/graph/badge.svg)](https://codecov.io/gh/DenisTitovLab/EnzymeRates.jl)
[![JET](https://img.shields.io/badge/%F0%9F%9B%A9%EF%B8%8F_tested_with-JET.jl-233f9a)](https://github.com/aviatesk/JET.jl)
[![Aqua](https://raw.githubusercontent.com/JuliaTesting/Aqua.jl/master/badge.svg)](https://github.com/JuliaTesting/Aqua.jl)

Identify the best enzyme rate equation from kinetic data. Given a reaction
definition and experimental rate measurements at varying substrate, product,
and regulator concentrations, EnzymeRates enumerates biochemically valid
mechanisms, fits each to the data, and selects the simplest mechanism that
generalizes by leave-one-group-out cross-validation. It has first-class
support for MWC allostery and for mechanisms that mix steady-state and
rapid-equilibrium elementary steps, and it derives Haldane/Wegscheider
thermodynamic constraints automatically from the mechanism's cycle structure.

## Installation

```julia
using Pkg
Pkg.add(url="https://github.com/DenisTitovLab/EnzymeRates.jl")
```

## Quickstart

Define a mechanism, derive its symbolic rate equation, and evaluate it:

```julia
using EnzymeRates

m = @enzyme_mechanism begin
    substrates: S
    products:   P
    steps: begin
        E + S ⇌ E(S)
        E(S) <--> E(P)
        E(P) ⇌ E + P
    end
end

parameters(m)                              # parameter names to supply
print(rate_equation_string(m))             # the symbolic rate equation
rate_equation(m, (S=1e-4, P=1e-5),        # evaluate numerically
    (K_P_E=1e-5, K_S_E=1e-4, k_ES_to_EP=100.0, Keq=2.0, E_total=1.0))
```

## Documentation

Full documentation — tutorials for deriving, fitting, and identifying rate
equations, plus the architecture and API reference — lives at the
[documentation site](https://DenisTitovLab.github.io/EnzymeRates.jl/stable/):

- [Getting Started](https://DenisTitovLab.github.io/EnzymeRates.jl/stable/getting_started/) — the end-to-end define → derive → fit → identify arc.
- [Deriving rate equations](https://DenisTitovLab.github.io/EnzymeRates.jl/stable/deriving/textbooks/) — RE vs steady state, the Cha/King–Altman algorithm, thermodynamic constraints, ping-pong, and MWC allostery.
- [Fitting rate equations](https://DenisTitovLab.github.io/EnzymeRates.jl/stable/fitting/tutorial/) — the data format, normalized vs absolute rate, and optimizer choice.
- [Identifying the best rate equation](https://DenisTitovLab.github.io/EnzymeRates.jl/stable/identify/tutorial/) — the enumeration engine, beam search, and cross-validation model selection.
- [API Reference](https://DenisTitovLab.github.io/EnzymeRates.jl/stable/api/) — every exported type, macro, and function.
