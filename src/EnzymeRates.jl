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
using SHA: sha256

include("Mechanisms/Mechanisms.jl")
using .Mechanisms: Mechanisms, EnzymeReaction, EnzymeMechanism, AllostericEnzymeMechanism,
    metabolites, Mechanism, AllostericMechanism, compile_mechanism
include("macros.jl")
include("rate_eq_derivation.jl")
include("thermo.jl")
include("poly.jl")
include("engine.jl")
include("normalize.jl")
include("print.jl")
include("key.jl")
include("codegen.jl")
include("admission.jl")
include("fitting.jl")
include("mechanism_enumeration.jl")
include("identify_rate_equation.jl")

# Every mechanism macro call runs the guard, the flip search runs `_faces_pass` and the
# first `rate_equation` call of a mechanism runs `_rate_body`, so all three compile into
# the package image rather than in each session.
precompile(_guard, (Mechanisms.Mechanism,))
precompile(_guard, (Mechanisms.AllostericMechanism,))
precompile(_faces_pass, (Mechanisms.Mechanism, BitVector))
precompile(_faces_pass, (Mechanisms.AllostericMechanism, BitVector))
precompile(_rate_body, (Mechanisms.Mechanism,))
precompile(_rate_body, (Mechanisms.AllostericMechanism,))

end
