# ABOUTME: Form names, the two sides of a step, and forms and steps as DSL text.

# ─── Form names, step sides, and forms and steps as DSL text ───────────

"""Whether `x` is a reactant's modifier copy (`A::Mod`) in `m`; a modifier-only name has
no copy, so it is not one."""
_is_reactant_copy(m::Mechanism, x::Met) = x.modifier && _is_reactant(m, x.name)

"""`x` as names write it: its name, with `mod` on a reactant's modifier copy; a
modifier-only name stays bare."""
_entry_label(m::Mechanism, x::Met) =
    _is_reactant_copy(m, x) ? string(x.name, "mod") : String(x.name)

"""Node `f`'s name: its label, its bound metabolites' labels in sorted order, then any
residual as `_res`, `_plus<name>` per substrate it adds and `_minus<name>` per product it
subtracts (`E`, `EAB`, `EPmodQ`, `EB_res_plusA_minusP`)."""
function _form_name(m::Mechanism, f::Int)
    x = m.nodes[f]
    io = IOBuffer()
    print(io, x.conformation)
    foreach(y -> print(io, _entry_label(m, y)), x.bound)
    if _has_residual(x)
        print(io, "_res")
        foreach(n -> print(io, "_plus", n), x.added)
        foreach(n -> print(io, "_minus", n), x.subtracted)
    end
    String(take!(io))
end

"""`x` as a mechanism macro writes a metabolite of `m`: its name, with `::Mod` on a
reactant's modifier copy; a modifier-only name stays bare."""
_met_text(m::Mechanism, x::Met) =
    _is_reactant_copy(m, x) ? "$(x.name)::Mod" : String(x.name)

"""Node `f` as a mechanism macro writes it: its label alone when it binds nothing and
keeps no residual (`E`, `Estar`), else its label with its bound metabolites (`_met_text`)
and residual in parentheses: `E(NAD, P)`, `E(A::Mod)`, `E(B; residual = A - P)`,
`E(; residual = A - P)`."""
function form_text(m::Mechanism, f::Int)
    x = m.nodes[f]
    label = String(x.conformation)
    _is_free(x) && return label
    terms = join(x.added, " + ") * join(" - $n" for n in x.subtracted)
    res = _has_residual(x) ? "; residual = " * lstrip(terms) : ""
    label * "(" * join((_met_text(m, y) for y in x.bound), ", ") * res * ")"
end

"""A step side: node `f`'s name and the labels of the free metabolites `free`, joined by
`_` (`EA_B` for `E(A) + B`)."""
_side_name(m::Mechanism, f::Int, free::Vector{Met}) =
    join([_form_name(m, f); [_entry_label(m, x) for x in free]], "_")
"""A step's two sides in its stored direction; every step constant is named after
them."""
_sides(m::Mechanism, e::Step) =
    (_side_name(m, e.from, e.taken_up), _side_name(m, e.to, e.given_off))

"""The names of step `k`'s two sides in its stored direction, `(from, to)`; every step
constant is named after them."""
step_sides(m::Mechanism, k::Int) = _sides(m, m.steps[k])

"""Step `k` of `m` as a mechanism macro writes it, in its stored direction: each side is
its node (`form_text`) and the metabolites the step takes up or gives off there
(`_met_text`), joined by `⇌` (rapid equilibrium) or `<-->` (steady state):
`E + S ⇌ E(S)`, `E(A) + B <--> E(Q) + P`."""
function step_text(m::Mechanism, k::Int)
    e = m.steps[k]
    side(f, free) = join([form_text(m, f); [_met_text(m, x) for x in free]], " + ")
    side(e.from, e.taken_up) * (e.rapid_equilibrium ? " ⇌ " : " <--> ") *
        side(e.to, e.given_off)
end
