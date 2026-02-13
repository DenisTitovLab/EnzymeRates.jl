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

# Mechanism enumeration
export SiteState, EnzymeFormSpec, enumerate_enzyme_forms
export enumerate_mechanisms, MechanismSpec, n_sites

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
