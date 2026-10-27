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
all_parameters(::Union{Mechanisms.AbstractEnzymeMechanism, AllostericMechanism, Mechanism}, ::NamedTuple)
fit_rate_equation(::FittingProblem, ::Any)
identify_rate_equation(::IdentifyRateEquationProblem)
metabolites(::M) where M<:Mechanisms.AbstractEnzymeMechanism
parameters(::M) where M<:Mechanisms.AbstractEnzymeMechanism
rate_equation(::M_type, ::NamedTuple, ::NamedTuple) where M_type<:Mechanisms.AbstractEnzymeMechanism
rate_equation_string(::Union{AllostericMechanism, Mechanism})
rescale_parameter_values(::Mechanisms.AbstractEnzymeMechanism, ::NamedTuple, ::AbstractVector{<:Pair})
```

## Macros

```@docs
@allosteric_mechanism
@enzyme_mechanism
@enzyme_reaction
```

## Types

```@docs
AllostericEnzymeMechanism
AllostericEnzymeMechanism(::AllostericMechanism)
EnzymeMechanism
EnzymeMechanism(::Mechanism)
EnzymeReaction
FittingProblem
FittingProblem(::Mechanisms.AbstractEnzymeMechanism, ::Any)
FittingProblem(::Union{AllostericMechanism, Mechanism}, ::Any)
IdentifyRateEquationProblem
IdentifyRateEquationResults
```
