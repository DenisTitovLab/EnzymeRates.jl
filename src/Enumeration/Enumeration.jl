# ABOUTME: The Enumeration module: the candidate mechanisms of a reaction, from its seed
# ABOUTME: mechanisms and the expansion moves that grow them, pruned by structural filters.
"""Generates candidate mechanisms for a reaction."""
module Enumeration

using ..Mechanisms: Mechanisms
using ..Derivation: Derivation
using Distributed: pmap

public CandidateMechanism, candidate_mechanisms

"""
    CandidateMechanism

A mechanism the enumerator proposes for a reaction (`candidate_mechanisms`). Its field
`mechanism`, a `Mechanism` or an `AllostericMechanism`, is the mechanism to fit and
report. Two candidates are `==` and hash alike exactly when their mechanisms are.
"""
struct CandidateMechanism
    mechanism::Union{Mechanisms.Mechanism, Mechanisms.AllostericMechanism}
end

Base.:(==)(a::CandidateMechanism, b::CandidateMechanism) = a.mechanism == b.mechanism
Base.hash(c::CandidateMechanism, h::UInt) = hash(c.mechanism, h)

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

"""The candidates of `mechs` whose law has at most `eq_complexity_filter` terms and at
most `max_param_count` fitted constants, in order. A limit of `typemax(Int)` limits
nothing, so its count is not computed: the fitted count solves the mechanism."""
_candidates(mechs, max_param_count::Int, eq_complexity_filter::Int) =
    CandidateMechanism[CandidateMechanism(m) for m in mechs
                       if (eq_complexity_filter == typemax(Int) ||
                           Derivation.rate_equation_term_count(m) <= eq_complexity_filter) &&
                          (max_param_count == typemax(Int) ||
                           Derivation.fitted_parameter_count(m) <= max_param_count)]

"""Whether the enumerator keeps `m`: the guard (`Derivation.check_valid_rate_equation`)
accepts it, and no conformation keeps a form that only bindings it excludes enter
(`_kept_forms_keep_a_binding`). The test the enumerator's exits and its flip, split and
tagging searches apply to a candidate."""
_enumerator_keeps(m) =
    Derivation.check_valid_rate_equation(Bool, m) && _kept_forms_keep_a_binding(m)

"""The number of entries named `name` in the metabolite list `v`, a reactant's modifier
copies included."""
_copies(v::Vector{Mechanisms.Metabolite}, name::Symbol) = count(x -> x.name === name, v)

"""The rapid-equilibrium flag of each step of `m`, in step order."""
_rapid_equilibrium_flags(m::Mechanisms.Mechanism) =
    BitVector([e.rapid_equilibrium for e in m.steps])

include("seeds.jl")
include("moves.jl")
include("filters.jl")
include("expand.jl")

end
