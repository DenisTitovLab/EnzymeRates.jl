# ABOUTME: Top-level module; defines the public API and include order.
# ABOUTME: Exports the 19 public names and wires together src/ source files.
module EnzymeRates

# Types
export EnzymeReaction, EnzymeMechanism, AllostericEnzymeMechanism
export IndependentParameters, FittedParameters, AllParameters
export FittingProblem
export IdentifyRateEquationProblem, IdentifyRateEquationResults

# DSL
export @enzyme_reaction, @enzyme_mechanism, @allosteric_mechanism

# Core API
export rate_equation, rate_equation_string, parameters, metabolites
export rescale_parameter_values

# Fitting & model selection
export fit_rate_equation
export identify_rate_equation


using Dates
using Distributed

include("Mechanisms/Mechanisms.jl")
using .Mechanisms: Mechanisms, EnzymeReaction, EnzymeMechanism, AllostericEnzymeMechanism,
    metabolites, Mechanism, AllostericMechanism, compile_mechanism
include("Derivation/Derivation.jl")
using .Derivation: Derivation, rate_equation, rate_equation_string, parameters,
    IndependentParameters, FittedParameters, AllParameters
include("macros.jl")
include("Fitting/Fitting.jl")
using .Fitting: Fitting, FittingProblem, fit_rate_equation, rescale_parameter_values
include("Enumeration/Enumeration.jl")
include("identify_rate_equation.jl")

end
