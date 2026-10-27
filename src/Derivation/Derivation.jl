# ABOUTME: The Derivation module: the well-posedness guard, the joint thermodynamic solve
# ABOUTME: and each mechanism's rate law (engine, unit rule, printer, key, generated code).
module Derivation

using ..Mechanisms: Mechanisms
using LinearAlgebra: det
using SHA: sha256

const API = (:_Canonical, :_NO_FAULT, :_admission, :_bottomless_segment, :_canonical,
    :_check_well_posed, :_eq_complexity, :_faces_pass, :_guard, :_has_turnover_cycle,
    :_independent_param_count, :_is_phantom, :_is_well_posed, :_law_text, :_rate_body,
    :_solve, :_ss_rate_constant_names, :_well_posed, :_with_law, :all_parameters,
    :fitted_params, :parameters, :rate_equation, :rate_equation_string)

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

end
