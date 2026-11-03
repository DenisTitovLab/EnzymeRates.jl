# ABOUTME: The DSL module: the macros @enzyme_reaction, @enzyme_mechanism and
# ABOUTME: @allosteric_mechanism, their reader, and parse_mechanism for mechanism text.
"""Turns macro-written reactions and mechanisms
into checked objects."""
module DSL

using ..Mechanisms: Mechanisms
using ..Derivation: Derivation

public parse_mechanism, @enzyme_reaction, @enzyme_mechanism, @allosteric_mechanism

include("read.jl")
include("macros.jl")

end
