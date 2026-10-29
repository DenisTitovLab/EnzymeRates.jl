# ABOUTME: The Derivation module: the well-posedness guard, the joint thermodynamic solve
# ABOUTME: and each mechanism's rate law (engine, unit rule, printer, key, generated code).
module Derivation

using ..Mechanisms: Mechanisms
using LinearAlgebra: det
using SHA: sha256

const API = (:AllParameters, :FittedParameters, :IndependentParameters,
    :check_valid_rate_equation, :conformation_steps, :excludes, :has_turnover_cycle,
    :has_unbalanced_cycle, :parameter_count, :parameters, :rate_equation,
    :rate_equation_complexity, :rate_equation_dedup_key, :rate_equation_string,
    :re_segments, :scale_rate_constants, :segment_graph)

"""A plain or an allosteric mechanism."""
const _AnyMech = Union{Mechanisms.Mechanism, Mechanisms.AllostericMechanism}

"""A column of the joint solve (`constant_columns`): slot `slot` of step `edge`
of the catalytic graph, in the inactive conformation when `inactive`."""
const _Column = @NamedTuple{edge::Int, slot::Int, inactive::Bool}

"""The two conformations by the index the guard gives them: 1 active, 2 inactive."""
const _CONFORMATIONS = (:active, :inactive)

"""The number of copies of metabolite `name` in `held`, a reactant's modifier copies
included."""
_copies(held::Vector{Mechanisms.Met}, name::Symbol) = count(x -> x.name === name, held)

"""The rapid-equilibrium flag of each step of `m`."""
_re_flags(m::Mechanisms.Mechanism) = BitVector(s.rapid_equilibrium for s in m.steps)

include("conformations.jl")
include("segments.jl")
include("names.jl")
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
