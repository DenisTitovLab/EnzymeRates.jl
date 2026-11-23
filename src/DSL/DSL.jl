# ABOUTME: The DSL module: the macros @enzyme_reaction and @enzyme_mechanism, their
# ABOUTME: reader, and parse_mechanism for mechanism text.
"""Turns macro-written reactions and mechanisms
into checked objects."""
module DSL

using ..Mechanisms: Mechanisms
using ..Derivation: Derivation

public parse_mechanism, @enzyme_reaction, @enzyme_mechanism

# ─── Helpers both macros use ────────────────────────────────────────────

"""`v` as an `Int`; errors unless it is a positive integer literal, which `true` is not.
`what` names it."""
_positive_int(v, what) = v isa Integer && !(v isa Bool) && v >= 1 ? Int(v) :
    error("$what must be a positive Int, got $v.")

"""
Parse a labeled-line, returning `(label, values_vector)`. `values_vector` is the
list of args after the label, in source order. Values are returned as parsed
(bare Symbols, `name::Tag`, `Sym[atoms]`, `begin ... end` blocks, literals, ...);
each caller checks the shape it expects. Handles both Julia parse shapes:
  - `Expr(:call, :(:), label, value)` — single labeled value.
  - `Expr(:tuple, Expr(:call, :(:), label, first), rest...)` — multi-element labeled.
"""
function _parse_labeled_line(arg)
    head, rest = arg isa Expr && arg.head == :tuple && !isempty(arg.args) ?
                 (arg.args[1], arg.args[2:end]) : (arg, Any[])
    head isa Expr && head.head == :call && head.args[1] == :(:) ||
        error("Expected `label: value` or `label: v1, v2, ...`; got $arg")
    head.args[2], Any[head.args[3], rest...]
end

include("enzyme_reaction_macro.jl")
include("enzyme_mechanism_macro.jl")

end
