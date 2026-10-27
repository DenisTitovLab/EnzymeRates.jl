# ABOUTME: The Enumeration module: the seed mechanisms of a reaction, the expansion moves
# ABOUTME: that grow them and the structural filters that prune the search.
module Enumeration

using ..Mechanisms: Mechanisms
using ..Derivation: Derivation
using Distributed: pmap

const API = (:_assert_atom_conserving, :_degenerate, :_expand_re_to_ss, :expand_mechanisms,
    :init_mechanisms, :seed_mechanisms)

include("seeds.jl")
include("moves.jl")
include("filters.jl")
include("expand.jl")

end
