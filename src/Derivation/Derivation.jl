# ABOUTME: The Derivation module: the well-posedness guard, the joint thermodynamic solve
# ABOUTME: and each mechanism's rate law (engine, unit rule, printer, key, generated code).
module Derivation

using ..Mechanisms: Mechanisms
using LinearAlgebra: det
using SHA: sha256

const API = (:all_parameters, :check_valid_rate_equation, :fitted_params,
    :has_valid_rate_equation, :parameter_count, :parameters, :rate_equation,
    :rate_equation_complexity, :rate_equation_key, :rate_equation_string,
    :scale_rate_constants)

include("rate_eq_derivation.jl")
include("thermo.jl")
include("poly.jl")
include("engine.jl")
include("normalize.jl")
include("print.jl")
include("key.jl")
include("codegen.jl")
include("guard.jl")
include("limit_test.jl")
include("admission.jl")
include("messages.jl")
include("memo.jl")

# Every mechanism macro call and the flip, split and tagging searches run the guard, and
# the first `rate_equation` call of a mechanism runs `_rate_body`, so both compile into the
# package image rather than in each session.
precompile(_guard, (Mechanisms.Mechanism,))
precompile(_guard, (Mechanisms.AllostericMechanism,))
precompile(_rate_body, (Mechanisms.Mechanism,))
precompile(_rate_body, (Mechanisms.AllostericMechanism,))

end
