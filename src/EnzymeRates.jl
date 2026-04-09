module EnzymeRates

# Types
export EnzymeReaction, EnzymeMechanism, AllostericEnzymeMechanism
export FittingProblem
export IdentifyRateEquationProblem, IdentifyRateEquationResults

# DSL
export @enzyme_reaction, @enzyme_mechanism

# Rate equation modes
export Full, Reduced

# Core API
export rate_equation, rate_equation_string, parameters, metabolites
export structural_identifiability_deficit
export rescale_parameter_values

# Fitting & model selection
export fit_rate_equation
export identify_rate_equation


using Tables
using Optimization

include("types.jl")
include("dsl.jl")
include("sym_poly_for_rate_eq_derivation.jl")
include("rate_eq_derivation.jl")
include("thermodynamic_constr_for_rate_eq_derivation.jl")
include("fitting.jl")
include("mechanism_enumeration.jl")
include("identify_rate_equation.jl")

end
