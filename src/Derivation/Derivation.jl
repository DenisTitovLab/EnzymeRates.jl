# ABOUTME: The Derivation module: its public names, its files in load order, what the files
# ABOUTME: share (_SwitchOff, _Reading, _excluded_steps, helpers) and the derivation memo.
"""Derives a mechanism's rate equation and compiles it."""
module Derivation

using ..Mechanisms: Mechanisms
using LinearAlgebra: det
using SHA: sha256

public check_valid_rate_equation, compile_rate_equation, dependent_parameter_values,
    fitted_parameter_count, parameter_units, parameters, rate_equation_term_count,
    rate_equation_dedup_key, rate_equation_string, RateEquationTooLargeError

# ─── Helpers several files use ──────────────────────────────────────────

_is_allosteric(m::Mechanisms.Mechanism) = !isempty(m.allosteric.allosteric_states)

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
How a reading switches off a catalytic step in one conformation (`_Reading`):
- `_RUNS`: the conformation keeps the step;
- `_LOSES_FORWARD`, `_LOSES_REVERSE`: the step loses the rate of its stored direction, or
  of its reverse, the direction that runs uphill, into the step's higher side. The other
  direction keeps its rate, so a rapid-equilibrium step stays fast;
- `_STOPPED`: both its rates vanish;
- `_CANNOT_SWITCH_OFF`: the reading has no way to switch the step off, which makes the
  reading inconsistent. Only a built reading does this, to an excluded rapid-equilibrium
  step it cannot switch off (`_built_reading`). The ε-network (`_eps_network`) treats such
  a step as slow and vanishing both ways, as it treats `_STOPPED`. A searched reading
  (`_consistent_columns`) may switch off such a step when it lies between two islands, the
  sets of forms a conformation's kept steps join (`exclusion_components`).
The developer page, "Well-posedness guard", says which readings the guard accepts.
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

# ─── The derivation memo ────────────────────────────────────────────────

"""
The derivation memo of a task (`_memo`): the mechanism `m` last asked about, and its guard
verdict (`_guard`), joint solve (`_solve`) and law (`_law`), each filled on first need. It
holds one mechanism, replaced by the next one asked about (compared by `===`), and lives
in the task's local storage, so tasks share nothing. A compiled rate law's generated
methods (`_derivation_input`) never read it. The developer page, "The law, its code, its
text and its key", says why the memo exists.
"""
mutable struct _Memo
    m::Union{Nothing, Mechanisms.Mechanism}
    verdict::Union{Nothing, _Verdict}
    solve::Union{Nothing, _Solve}
    law::Union{Nothing, _Law, _AlloLaw}
end

"""The task's memo (`_Memo`), emptied and given to `m` unless it holds `m`."""
function _memo(m::Mechanisms.Mechanism)
    c = get!(() -> _Memo(nothing, nothing, nothing, nothing), task_local_storage(),
             _Memo)::_Memo
    c.m === m || (c.m = m; c.verdict = nothing; c.solve = nothing; c.law = nothing)
    c
end

function _memo_verdict(m::Mechanisms.Mechanism)
    v = _memo(m).verdict
    v === nothing || return v
    _memo(m).verdict = _guard(m)
end

function _memo_solve(m::Mechanisms.Mechanism)
    s = _memo(m).solve
    s === nothing || return s
    _memo(m).solve = _solve(m)
end

"""`m`'s law (`_law`) from its joint solve, through the memo: a `_Law` for a plain
mechanism, an `_AlloLaw` for an allosteric one."""
function _memo_law(m::Mechanisms.Mechanism)
    L = _memo(m).law
    L === nothing || return L::(_is_allosteric(m) ? _AlloLaw : _Law)
    _memo(m).law = _law(m, _memo_solve(m))
end

_with_memo_law(f, m::Mechanisms.Mechanism) = f(m, _memo_solve(m), _memo_law(m))

# Every mechanism macro call and the flip, split and tagging searches run the guard, and
# the first call of a mechanism's compiled rate law decodes its Sig (`_derivation_input`)
# and runs `_rate_body`, so all of them compile into the package image rather than in
# each session.
precompile(_guard, (Mechanisms.Mechanism,))
precompile(_derivation_input, (Type{<:_CompiledRateEquation},))
precompile(_mechanism_from_sig, (Tuple,))
precompile(_rate_body, (Mechanisms.Mechanism,))

end
