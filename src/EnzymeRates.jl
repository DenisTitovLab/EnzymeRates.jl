# ABOUTME: Top-level module; defines the public API and include order.
# ABOUTME: Exports the 17 public names and wires together src/ source files.
module EnzymeRates

# Types
export EnzymeReaction, EnzymeMechanism, AllostericEnzymeMechanism
export FittingProblem
export IdentifyRateEquationProblem, IdentifyRateEquationResults

# DSL
export @enzyme_reaction, @enzyme_mechanism, @allosteric_mechanism

# Core API
export rate_equation, rate_equation_string, parameters, all_parameters, metabolites
export rescale_parameter_values

# Fitting & model selection
export fit_rate_equation
export identify_rate_equation


using Dates
using LinearAlgebra: det
using Tables
using Optimization
using Distributed
using FunctionWrappers: FunctionWrapper
import Unicode

include("types.jl")
include("names.jl")
include("dsl.jl")
include("sym_poly_for_rate_eq_derivation.jl")
include("rate_eq_derivation.jl")
include("thermo.jl")
include("poly.jl")
include("engine.jl")
include("admission.jl")
include("fitting.jl")
include("mechanism_enumeration.jl")
include("identify_rate_equation.jl")

# Every mechanism macro call runs the guard and the flip search runs `_faces_pass`, so both
# compile into the package image rather than in each session.
precompile(_guard, (Mechanism,))
precompile(_guard, (AllostericMechanism,))
precompile(_faces_pass, (Mechanism, BitVector))
precompile(_faces_pass, (AllostericMechanism, BitVector))

end
