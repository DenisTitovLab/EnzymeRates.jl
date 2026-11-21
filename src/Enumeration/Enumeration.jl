# ABOUTME: The Enumeration module: the candidate mechanisms of a reaction, from its seed
# ABOUTME: mechanisms and the expansion moves that grow them, pruned by structural filters.
"""Generates candidate mechanisms for a reaction."""
module Enumeration

using ..Mechanisms: Mechanisms
using ..Derivation: Derivation
using Distributed: pmap

public CandidateMechanism, candidate_mechanisms

include("record.jl")

"""
    candidate_mechanisms(reaction; max_param_count = typemax(Int),
                         eq_complexity_filter = typemax(Int),
                         optional_allosteric_regulators = Symbol[],
                         optional_competitive_inhibitors = Symbol[])
        -> Vector{CandidateMechanism}
    candidate_mechanisms(reaction, parent::CandidateMechanism;
                         max_param_count = typemax(Int),
                         eq_complexity_filter = typemax(Int))
        -> Vector{CandidateMechanism}

The candidate mechanisms of the `EnzymeReaction` `reaction`. The first method returns the
mechanisms a search fits first: they bind every regulator `reaction` declares but those
the two optional lists name, and no regulator those lists name. The second returns the
children of `parent`, a candidate of `reaction`, by every expansion move: each changes
one thing, such as a step's speed, a kinetic group's split, a bound inhibitor or an
allosteric property.

Every candidate returned has a well-posed rate law (`Derivation.check_valid_rate_equation`)
with at most `eq_complexity_filter` terms (`Derivation.rate_equation_term_count`) and at
most `max_param_count` fitted constants (`Derivation.fitted_parameter_count`). The limits
drop only what a method returns: a mechanism over a limit is never returned, so it is
never expanded either.

Errors: the reaction's metabolite names spell one string two ways (such as `Ac`, `CoA`
and `AcCoA`), so that two forms would share a name (first method); the reaction has no
catalytic cycle; no mechanism binds every required regulator; an internal check of the
enumerator fails, such as atom conservation.
"""
function candidate_mechanisms end

"""The candidates of `cands` whose law has at most `eq_complexity_filter` terms and at
most `max_param_count` fitted constants, in order. A limit of `typemax(Int)` limits
nothing, so its count is not computed: the fitted count solves the mechanism."""
_candidates(cands::Vector{CandidateMechanism}, max_param_count::Int,
            eq_complexity_filter::Int) =
    filter(cands) do c
        (eq_complexity_filter == typemax(Int) ||
         Derivation.rate_equation_term_count(c.mechanism) <= eq_complexity_filter) &&
            (max_param_count == typemax(Int) ||
             Derivation.fitted_parameter_count(c.mechanism) <= max_param_count)
    end

"""Whether the enumerator keeps `m`: the guard (`Derivation.check_valid_rate_equation`)
accepts it, and no conformation keeps a form that only bindings it excludes enter
(`_kept_forms_keep_a_binding`). The test the enumerator's exits and its flip, split and
tagging searches apply to a candidate."""
_enumerator_keeps(m) =
    Derivation.check_valid_rate_equation(Bool, m) && _kept_forms_keep_a_binding(m)

"""The number of entries named `name` in the metabolite list `v`, a reactant's modifier
copies included."""
_copies(v::Vector{Mechanisms.Metabolite}, name::Symbol) = count(x -> x.name === name, v)

"""Whether `m` is allosteric: its kinetic groups carry allosteric states."""
_is_allosteric(m::Mechanisms.Mechanism) = !isempty(m.allosteric.allosteric_states)

"""The rapid-equilibrium flag of each step of `m`, in step order."""
_rapid_equilibrium_flags(m::Mechanisms.CatalyticScheme) =
    BitVector([e.rapid_equilibrium for e in m.steps])

"""Each kinetic group of graph `m`, its steps in step order."""
function _kinetic_groups(m::Mechanisms.CatalyticScheme)
    members = [Int[] for _ in 1:maximum(e -> e.group, m.steps; init = 0)]
    for (k, e) in enumerate(m.steps)
        push!(members[e.group], k)
    end
    members
end

"""How many more modifier copies of `name` the record of candidate `c` holds at the second
form of step `k` than at its first: positive when the step takes copies of `name` up,
negative when it gives them off, zero when it exchanges `name` only productively. A copy
binds a site of its own and takes no part in the chemistry, so it enters or leaves a form
only with the step that takes it up or gives it off."""
function _copy_change(c::CandidateMechanism, k::Int, name::Symbol)
    e = c.mechanism.catalytic.steps[k]
    count(==(name), c.forms[e.to].copies) - count(==(name), c.forms[e.from].copies)
end

"""Whether step `k` of candidate `c` is a plain binding: it takes up one metabolite and
gives off none (`_is_binding`), and the record of its second form is that of its first
with the metabolite added, among the copies when the step binds a copy (`_copy_change`),
the residual the same and the label free to change. Every other step is chemistry: an
isomerization, a fused binding (the last substrate binding straight into the
product-bound form), a Theorell–Chance step or an exchange."""
function _is_plain_binding(c::CandidateMechanism, k::Int)
    e = c.mechanism.catalytic.steps[k]
    _is_binding(e) || return false
    f, t, x = c.forms[e.from], c.forms[e.to], only(e.taken_up).name
    f.adds == t.adds && f.subtracts == t.subtracts || return false
    _copy_change(c, k, x) > 0 ?
        f.productive == t.productive && _one_more(t.copies, f.copies, x) :
        f.copies == t.copies && _one_more(t.productive, f.productive, x)
end

"""Whether the sorted list `t` is the sorted list `f` with one more `x`."""
function _one_more(t::Vector{Symbol}, f::Vector{Symbol}, x::Symbol)
    length(t) == length(f) + 1 || return false
    j, skipped = 1, false
    for y in t
        if j <= length(f) && y === f[j]
            j += 1
        elseif !skipped && y === x
            skipped = true
        else
            return false
        end
    end
    skipped
end

include("seeds.jl")
include("moves.jl")
include("filters.jl")
include("expand.jl")

end
