module EnzymeRates

# Types
export EnzymeReaction, EnzymeMechanism
export FittingProblem
# export IdentifyRateEquationProblem, IdentifyRateEquationResults  # when implemented

# DSL
export @enzyme_reaction, @enzyme_mechanism

# Rate equation modes
export Full, Reduced

# Core API
export rate_equation, rate_equation_string, parameters, metabolites
export structural_identifiability_deficit

# Fitting & model selection
export fit_rate_equation
# export identify_rate_equation  # when implemented

# Mechanism enumeration (internal types, not public API)
# SiteState, EnzymeFormSpec, MechanismSpec, EnumerationStage subtypes
# enumerate_enzyme_forms, enumerate_mechanisms are accessible via EnzymeRates.*

using Tables
using Optimization

include("types.jl")
include("dsl.jl")
include("sym_poly_for_rate_eq_derivation.jl")
include("rate_eq_derivation.jl")
include("rate_eq_rewriting.jl")
include("fitting.jl")
include("mechanism_enumeration.jl")

end
