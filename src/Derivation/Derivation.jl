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

"""The rapid-equilibrium flag of each step of `m`."""
_re_flags(m::Mechanisms.CatalyticScheme) = BitVector(s.rapid_equilibrium for s in m.steps)

"""Whether the equilibrium constant of step `e` at rapid equilibrium is named after the
reverse of its stored direction. A binding, a step that takes up or gives off one
metabolite and nothing else, has a dissociation constant, named after its release; so a
binding stored as its uptake names its constant in reverse. Every other step's constant is
named after its stored direction."""
_constant_named_in_reverse(e::Mechanisms.Step) =
    length(e.taken_up) == 1 && isempty(e.given_off)

"""The kinetic groups of `m`, numbered by their steps' `kinetic_group`: each group's steps,
in step order."""
function _kinetic_groups(m::Mechanisms.CatalyticScheme)
    G = maximum(e -> e.kinetic_group, m.steps; init = 0)
    [findall(e -> e.kinetic_group == g, m.steps) for g in 1:G]
end

# ─── Excluded steps and union-find, which several files use ─────────────

"""The steps conformation `c` of `m` (1 active, 2 inactive) excludes by their kinetic
group's allosteric state, one flag per catalytic step: the active conformation those of
`:OnlyI` groups, the inactive one those of `:OnlyA` groups. A plain mechanism excludes
none."""
function _excluded_steps(m::Mechanisms.Mechanism, c::Int)
    state, states = c == 1 ? :OnlyI : :OnlyA, m.allosteric.allosteric_states
    (e -> !isempty(states) && states[e.kinetic_group] === state).(m.catalytic.steps)
end

function _union_find_root(parent::Vector{Int}, i::Int)
    while parent[i] != i
        i = parent[i]
    end
    i
end

# ─── A reading: how each conformation switches off a step ───────────────

"""
How a reading switches off a catalytic step in one conformation (`_Reading`): the
conformation keeps the step (`_RUNS`); the step loses the rate of its stored direction
(`_LOSES_FORWARD`) or of its reverse (`_LOSES_REVERSE`), the direction that runs uphill,
into the step's higher side, while the other keeps its rate, so a rapid-equilibrium step
stays fast; it is stopped, both its rates vanishing (`_STOPPED`); or it cannot be switched
off (`_CANNOT_SWITCH_OFF`).
"""
@enum _SwitchOff::UInt8 _RUNS _LOSES_FORWARD _LOSES_REVERSE _STOPPED _CANNOT_SWITCH_OFF

"""How a mechanism's conformations switch off its catalytic steps (`_SwitchOff`):
`reading[k, c]` for step `k` in conformation `c`, 1 active and 2 inactive. A plain
mechanism's reading has its one conformation only."""
const _Reading = Matrix{_SwitchOff}

include("thermodynamic_constraints.jl")
include("rate_equation_derivation.jl")
include("rate_equation_dedup_key.jl")
include("rate_equation_correctness_guard.jl")
include("allosteric_rate_equation_correctness_guard.jl")
include("memo.jl")

# Every mechanism macro call and the flip, split and tagging searches run the guard, and
# the first call of a mechanism's compiled rate law decodes its Sig (`_derivation_input`)
# and runs `_rate_body`, so all of them compile into the package image rather than in
# each session.
precompile(_guard, (Mechanisms.Mechanism,))
precompile(_derivation_input, (Type{<:_CompiledRateEquation},))
precompile(_mechanism_from_sig, (Tuple,))
precompile(_rate_body, (Mechanisms.Mechanism,))

end
