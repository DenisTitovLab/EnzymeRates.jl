# ABOUTME: The parameter-naming chokepoint name(p, m) and the renderers it builds on
# ABOUTME: (metabolite labels, form names, step sides), and forms and steps as DSL text.

# ─── Form names, step sides, and forms and steps as DSL text ───────────

"""Whether entry `i` of `m` is a reactant's modifier copy (`A::Mod`); a modifier-only
name has no copy, so its entry is not one."""
_is_reactant_copy(m::Mechanism, i::Int) =
    m.table[i].modifier && _is_reactant(m, m.table[i].name)

"""Entry `i` as names write it: its name, with `mod` on a reactant's modifier copy; a
modifier-only name stays bare."""
_entry_label(m::Mechanism, i::Int) =
    _is_reactant_copy(m, i) ? string(m.table[i].name, "mod") : String(m.table[i].name)

"""Form `f`'s name: its label, its bound entries' labels in table order, then any
residual as `_res`, `_plus<name>` per substrate it adds and `_minus<name>` per product it
subtracts (`E`, `EAB`, `EPmodQ`, `EB_res_plusA_minusP`)."""
function _form_name(m::Mechanism, f::Int)
    x = m.forms[f]
    io = IOBuffer()
    print(io, m.conformations[x.conformation])
    foreach(i -> print(io, _entry_label(m, i)), entries(x.bound))
    if _has_residual(x)
        print(io, "_res")
        foreach(i -> print(io, "_plus", m.table[i].name), entries(x.added))
        foreach(i -> print(io, "_minus", m.table[i].name), entries(x.subtracted))
    end
    String(take!(io))
end

"""Entry `i` of `m` as a mechanism macro writes a metabolite: its name, with `::Mod` on a
reactant's modifier copy; a modifier-only name stays bare."""
_met_text(m::Mechanism, i::Int) =
    _is_reactant_copy(m, i) ? "$(m.table[i].name)::Mod" : String(m.table[i].name)

"""Form `f` as a mechanism macro writes it: its label alone when it binds nothing and
keeps no residual (`E`, `Estar`), else its label with its bound metabolites (`_met_text`)
and residual in parentheses: `E(NAD, P)`, `E(A::Mod)`, `E(B; residual = A - P)`,
`E(; residual = A - P)`."""
function _form_text(m::Mechanism, f::Int)
    x = m.forms[f]
    label = String(m.conformations[x.conformation])
    _is_free(x) && return label
    terms = join(_entry_names(m, x.added), " + ") *
            join(" - $(m.table[i].name)" for i in entries(x.subtracted))
    res = _has_residual(x) ? "; residual = " * lstrip(terms) : ""
    label * "(" * join((_met_text(m, i) for i in entries(x.bound)), ", ") * res * ")"
end

"""A step side: form `f`'s name and the labels of the free entries `c`, joined by `_`
(`EA_B` for `E(A) + B`)."""
_side_name(m::Mechanism, f::Int, c::Counts) =
    join([_form_name(m, f); [_entry_label(m, i) for i in entries(c)]], "_")
"""An edge's two sides in its stored direction; every step constant is named after
them."""
_sides(m::Mechanism, e::Edge) =
    (_side_name(m, e.from, e.consumed), _side_name(m, e.to, e.released))
_sides(m::Mechanism, k::Int) = _sides(m, m.edges[k])

"""Edge `k` of `m` as a mechanism macro writes its step, in its stored direction: each
side is its form (`_form_text`) and the metabolites the step takes up or gives off there
(`_met_text`), joined by `⇌` (rapid equilibrium) or `<-->` (steady state):
`E + S ⇌ E(S)`, `E(A) + B <--> E(Q) + P`."""
function _step_text(m::Mechanism, k::Int)
    e = m.edges[k]
    side(f, c) = join([_form_text(m, f); [_met_text(m, i) for i in entries(c)]], " + ")
    side(e.from, e.consumed) * (e.rapid_equilibrium ? " ⇌ " : " <--> ") *
        side(e.to, e.released)
end

# ─── name(p::Parameter, m) chokepoint ─────────────────────────────────
#
# Single chokepoint for parameter Symbol production. Renders structural
# names (:K_ES_to_E_S, :k_ES_to_EP, :K_I_ES_to_E_S, :K_Ireg) via Parameter
# subtype dispatch. Routing all parameter-name production through one function
# keeps any name-scheme change a single-function edit.

"""
State token — placed right after the type prefix:
  :A       → "A_"  (allosteric active branch; distinguishes allosteric from non-)
  :I       → "I_"  (allosteric inactive branch: OnlyI or NonequalAI-inactive)
  :EqualAI → ""    (allosteric shared symbol; shared between A and I state)
  :None    → ""    (non-allosteric mechanism)
"""
_state_tag(state::Symbol) =
    state === :A ? "A_" : state === :I ? "I_" : state in (:EqualAI, :None) ? "" :
    error("_state_tag: unexpected Parameter.state $state " *
          "(must be one of :None, :EqualAI, :A, :I)")

"""
The constant of the reaction `a → b`: `prefix`, the state tag, then
`"<a>_to_<b>"` (e.g. `:k_ES_to_EP`, `:K_A_ES_to_E_S`).
"""
_render_reaction(prefix::String, (a, b)::Tuple{String, String}, state::Symbol) =
    Symbol(prefix, _state_tag(state), a, "_to_", b)

"""
The sides that name constant `slot` of edge `k` of `m`'s catalytic graph in `state`: those
of the kept member of its kinetic class (`kept`, `_kept_members(m)`). An `:EqualAI`
edge's constant is one class in both conformations; the inactive copy of an `:OnlyA`
edge is a class of its own, named after its own edge. A slot the edge's group has no
class for (the reverse rate of a rapid-equilibrium edge) takes the sides of its slot-1
class, so every parameter of every edge renders.
"""
function _constant_sides(m::_AnyMech, k::Int, slot::Int, state::Symbol,
                         kept::Dict{Tuple{Int, Int, Bool}, Int})
    cm = catalytic(m)
    g = cm.group[k]
    tag = m isa AllostericMechanism ? cat_allo_state(m, g) : :None
    state === :I && tag === :OnlyA && return _sides(cm, k)
    inactive = state === :I && tag === :NonequalAI
    _sides(cm, get(() -> kept[(g, 1, inactive)], kept, (g, slot, inactive)))
end

# Every edge constant is named after the reaction of its class's kept member: a rate
# constant (`k_`) after the direction it drives, an equilibrium constant (`K_`) after
# the direction whose products-over-reactants it is. A binding's K is named in the
# release direction, so it is a dissociation constant (`K_ES_to_E_S`); an
# isomerization's or other RE step's K in the stored direction (`K_ES_to_EP`). `kept` is
# `_kept_members(m)`, which a caller that names many constants of one mechanism computes
# once and passes.
name(p::Kfor, m::_AnyMech, kept = _kept_members(m)) =
    _render_reaction("k_", _constant_sides(m, p.edge, 1, p.state, kept), p.state)
name(p::Krev, m::_AnyMech, kept = _kept_members(m)) =
    _render_reaction("k_", reverse(_constant_sides(m, p.edge, 2, p.state, kept)),
                     p.state)
function name(p::Krapid, m::_AnyMech, kept = _kept_members(m))
    s = _constant_sides(m, p.edge, 1, p.state, kept)
    _render_reaction("K_", _is_binding(catalytic(m).edges[p.edge]) ? reverse(s) : s,
                     p.state)
end

"""
Regulator-site parameter: state tag + ligand name + "reg". No site index —
the AllostericMechanism constructor enforces that each ligand appears at most
once across all sites.
"""
name(p::Kreg, ::AllostericMechanism) =
    Symbol("K_", _state_tag(p.state), String(p.ligand), "reg")

"""
The Parameter(s) governing edge `k` of `m` in the given allosteric state, by its flag.
`_edge_names` applies it to every edge of one conformation, with that conformation's
tag, and `_column_parameter` picks a joint-solve column's.

Returns 1 element for an RE edge (`Krapid`) and 2 elements for an SS edge
(`Kfor`+`Krev`). `name` renders a binding's `Krapid`, plain or fused, as a
dissociation constant and every other edge's — an isomerization, a Theorell–Chance
step, a step with two metabolites on one side — as the equilibrium constant of its
stored direction.
"""
_edge_constants(m::Mechanism, k::Int, state::Symbol) =
    m.edges[k].rapid_equilibrium ? Parameter[Krapid(k, state)] :
                                   Parameter[Kfor(k, state), Krev(k, state)]
