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
fit_rate_equation(::FittingProblem, ::Any)
identify_rate_equation(::IdentifyRateEquationProblem)
metabolites(::M) where M<:Mechanisms.AbstractEnzymeMechanism
parameters(::Any)
dependent_parameter_values(::Union{Mechanisms.AbstractEnzymeMechanism, AllostericMechanism, Mechanism}, ::NamedTuple)
rate_equation(::M_type, ::NamedTuple, ::NamedTuple) where M_type<:Mechanisms.AbstractEnzymeMechanism
rate_equation_string(::Union{AllostericMechanism, Mechanism})
test_rate_equation(::Any, ::NamedTuple, ::Any)
```

## Macros

```@docs
@allosteric_mechanism
@enzyme_mechanism
@enzyme_reaction
```

## Types

```@docs
AllParameters
AllostericEnzymeMechanism
AllostericEnzymeMechanism(::AllostericMechanism)
EnzymeMechanism
EnzymeMechanism(::Mechanism)
EnzymeReaction
FittedParameters
FittingProblem
FittingProblem(::Mechanisms.AbstractEnzymeMechanism, ::Any)
FittingProblem(::Union{AllostericMechanism, Mechanism}, ::Any)
IdentifyRateEquationProblem
IdentifyRateEquationResults
IndependentParameters
```
