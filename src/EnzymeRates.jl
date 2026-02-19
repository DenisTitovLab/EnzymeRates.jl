module EnzymeRates

# Types
export AbstractEnzymeReaction, EnzymeReaction
export AbstractEnzymeMechanism, EnzymeMechanism

# DSL
export @enzyme_reaction, @enzyme_mechanism

# Core API
export rate_equation, rate_equation_string, parameters, metabolites

# Identifiability
export is_identifiable, structural_identifiability_deficit

# Fitting
export FittingProblem, fit_rate_equation
export fitted_params, metabolite_names

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
