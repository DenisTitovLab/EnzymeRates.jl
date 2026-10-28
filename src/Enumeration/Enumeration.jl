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

include("seeds.jl")
include("moves.jl")
include("filters.jl")
include("expand.jl")

end
