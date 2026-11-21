# ABOUTME: The DSL module: the macros @enzyme_reaction and @enzyme_mechanism, their
# ABOUTME: reader, and parse_mechanism for mechanism text.
"""Turns macro-written reactions and mechanisms
into checked objects."""
module DSL

using ..Mechanisms: Mechanisms
using ..Derivation: Derivation

public parse_mechanism, @enzyme_reaction, @enzyme_mechanism

include("read.jl")
include("macros.jl")

end
