# ABOUTME: The Derivation module: the well-posedness guard, the joint thermodynamic solve
# ABOUTME: and each mechanism's rate law (engine, unit rule, printer, key, generated code).
"""Derives a mechanism's rate equation and compiles it."""
module Derivation

using ..Mechanisms: Mechanisms
using LinearAlgebra: det
using SHA: sha256

public check_valid_rate_equation, compile_rate_equation, dependent_parameter_values,
    fitted_parameter_count, parameter_units, parameters, rate_equation_term_count,
    rate_equation_dedup_key, rate_equation_string

"""Whether `m` is allosteric: its kinetic groups carry allosteric states."""
_is_allosteric(m::Mechanisms.Mechanism) = !isempty(m.allosteric.allosteric_states)

"""The two conformations by the index the guard gives them: 1 active, 2 inactive."""
const _CONFORMATIONS = (:active, :inactive)

"""The number of copies of metabolite `name` in `held`, a reactant's modifier copies
included."""
_copies(held::Vector{Mechanisms.Metabolite}, name::Symbol) =
    count(x -> x.name === name, held)

"""The rapid-equilibrium flag of each step of `m`."""
_re_flags(m::Mechanisms.CatalyticScheme) = BitVector(s.rapid_equilibrium for s in m.steps)

"""Whether the equilibrium constant of step `e` at rapid equilibrium is named after the
reverse of its stored direction. A binding, a step that takes up or gives off one
metabolite and nothing else, has a dissociation constant, named after its release; so a
binding stored as its uptake names its constant in reverse. Every other step's constant is
named after its stored direction."""
_constant_named_in_reverse(e::Mechanisms.Step) =
    length(e.taken_up) == 1 && isempty(e.given_off)

"""The kinetic groups of `m`, numbered by their steps' `group`: each group's steps, in step
order."""
function _kinetic_groups(m::Mechanisms.CatalyticScheme)
    G = maximum(e -> e.group, m.steps; init = 0)
    [findall(e -> e.group == g, m.steps) for g in 1:G]
end

include("conformations.jl")
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
# the first call of a mechanism's compiled rate law decodes its Sig (`_derivation_input`)
# and runs `_rate_body`, so all of them compile into the package image rather than in
# each session.
precompile(_guard, (Mechanisms.Mechanism,))
precompile(_derivation_input, (Type{<:_CompiledRateEquation},))
precompile(Tuple{typeof(_mechanism_from_sig), Int, Vararg{Int}})
precompile(_rate_body, (Mechanisms.Mechanism,))

end
