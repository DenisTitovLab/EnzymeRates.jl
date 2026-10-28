# ABOUTME: The DSL module: the macros @enzyme_reaction, @enzyme_mechanism and
# ABOUTME: @allosteric_mechanism, their reader, and parse_mechanism for mechanism text.
module DSL

using ..Mechanisms: Mechanisms
using ..Derivation: Derivation

const API = (:parse_mechanism, Symbol("@enzyme_reaction"), Symbol("@enzyme_mechanism"),
    Symbol("@allosteric_mechanism"))

include("read.jl")
include("macros.jl")

end
