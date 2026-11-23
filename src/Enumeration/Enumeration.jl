# ABOUTME: The Enumeration module: the candidate mechanisms of a reaction, from its seed
# ABOUTME: mechanisms and the expansion moves that grow them, pruned by structural filters.
"""Generates candidate mechanisms for a reaction."""
module Enumeration

using ..Mechanisms: Mechanisms
using ..Derivation: Derivation
using Distributed: pmap

public CandidateMechanism, candidate_mechanisms

include("candidate_mechanism_type.jl")
include("expansion_moves.jl")
include("starting_mechanisms.jl")

end
