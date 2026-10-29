# ABOUTME: The Enumeration module: the seed mechanisms of a reaction, the expansion moves
# ABOUTME: that grow them and the structural filters that prune the search.
module Enumeration

using ..Mechanisms: Mechanisms
using ..Derivation: Derivation
using Distributed: pmap

const API = (:base_mechanisms, :expand_mechanisms)

"""Whether the guard (`Derivation.check_valid_rate_equation`) accepts `m`: the test the
enumerator's exits and its flip, split and tagging searches apply to a candidate."""
_guard_accepts(m) = Derivation.check_valid_rate_equation(m; throw = false)

"""The number of entries named `name` in the metabolite list `v`, a reactant's modifier
copies included."""
_copies(v::Vector{Mechanisms.Met}, name::Symbol) = count(x -> x.name === name, v)

"""The rapid-equilibrium flag of each step of `m`, in step order."""
_rapid_equilibrium_flags(m::Mechanisms.Mechanism) =
    BitVector([e.rapid_equilibrium for e in m.steps])

include("seeds.jl")
include("moves.jl")
include("filters.jl")
include("expand.jl")

end
