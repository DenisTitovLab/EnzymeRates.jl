# API Reference

```@meta
CurrentModule = EnzymeRates
```

This page lists every exported name. Maintainer-only internals are not
shown here — see the Developer / Architecture page for those.

## Index

```@index
Order = [:function, :macro, :type]
```

## Functions

```@docs
compile_rate_equation(::Mechanism)
dependent_parameter_values(::Mechanism, ::NamedTuple, ::Real)
fit_rate_equation(::FittingProblem, ::Any)
identify_rate_equation(::IdentifyRateEquationProblem)
metabolites(::Mechanism)
parameters(::Mechanism)
plot_rate_equation_on_data(::Mechanism, ::NamedTuple, ::Any, ::Real)
rate_equation_string(::Mechanism)
```

## Internal but usable

`test_rate_equation` is not exported; call it as `EnzymeRates.test_rate_equation`.
`RateEquationTooLargeError`, the error a rate law over the term cap throws, is not exported
either; catch it as `EnzymeRates.Derivation.RateEquationTooLargeError`.

```@docs
EnzymeRates.test_rate_equation(::Mechanism, ::NamedTuple, ::Any, ::Real)
EnzymeRates.Derivation.RateEquationTooLargeError
```

## Macros

```@docs
@enzyme_mechanism
@enzyme_reaction
```

## Types

```@docs
EnzymeReaction
FittingProblem
FittingProblem(::Mechanism, ::Any, ::Real)
IdentifyRateEquationProblem
IdentifyRateEquationResults
Mechanism
```
