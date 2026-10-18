# ABOUTME: The parameter-naming chokepoint name(p, m) and the renderers it builds on:
# ABOUTME: metabolite labels, enzyme-form names, step sides and step text.

# ─── Labels, form names and step text ─────────────────────────────────

"""Whether `x` is the competitive-inhibitor copy of a substrate or product of `ctx` (a
mechanism or a reaction; `_is_reactant`): a reactant bound as a modifier."""
_is_reactant_copy(x::Metabolite, ctx) =
    x isa CompetitiveInhibitor && _is_reactant(ctx, name(x))

"""
A metabolite as written in parameter names: its name, with suffix `mod` on a
competitive-inhibitor copy of a substrate or product of `ctx` (a mechanism or a
reaction; `_is_reactant_copy`), so the copy's forms and constants stay distinct from
those of the reactant itself. A competitive inhibitor that is no reactant binds only as
a modifier, so its name stands bare.
"""
_met_label(x::Metabolite, ctx) =
    _is_reactant_copy(x, ctx) ? String(name(x)) * "mod" : String(name(x))

"""
An enzyme form as parameter names write it for `ctx` (a mechanism or a reaction): its
conformation label, then the labels (`_met_label`) of its bound metabolites without a
separator, then, for a residual, `_res`, `_plus<name>` per substrate it adds and
`_minus<name>` per product it subtracts, each in name order: `E`, `EAB`, `EPmodQ`,
`EB_res_plusA_minusP`.
"""
function _form_label(sp::Species, ctx)
    head = String(conformation(sp)) * join(_met_label(x, ctx) for x in bound(sp))
    r = residual(sp)
    isempty(r) && return head
    head * "_res" * join("_plus" * String(name(a)) for a in added(r)) *
        join("_minus" * String(name(p)) for p in subtracted(r))
end

"""A step as text for `ctx` (a mechanism or a reaction): its entry side, an arrow (`⇌`
rapid equilibrium, `<-->` steady state) and its exit side. A side lists the enzyme form
as parameter names write it (`_form_label`), then the metabolites the step consumes or
releases."""
_step_text(s::Step, ctx) =
    join((_form_label(from_species(s), ctx), name.(consumed(s))...), " + ") *
    (is_equilibrium(s) ? " ⇌ " : " <--> ") *
    join((_form_label(to_species(s), ctx), name.(released(s))...), " + ")

# ─── Graph form names and step sides ─────────────────────────────────

"""Entry `i` as names write it: its name, with `mod` on a reactant's modifier copy; a
modifier-only name stays bare."""
_entry_label(m::Mechanism, i::Int) =
    m.table[i].modifier && _is_reactant(m, m.table[i].name) ?
        string(m.table[i].name, "mod") : String(m.table[i].name)

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

"""Form `f` as a mechanism macro writes it: `E(NAD, P)`, `E(A::Mod)`, `E(I)`,
`E(B; residual = A - P)`."""
function _form_text(m::Mechanism, f::Int)
    x = m.forms[f]
    mets = [_is_reactant(m, m.table[i].name) && m.table[i].modifier ?
            "$(m.table[i].name)::Mod" : "$(m.table[i].name)" for i in entries(x.bound)]
    terms = join(_entry_names(m, x.added), " + ") *
            join(" - $(m.table[i].name)" for i in entries(x.subtracted))
    res = _has_residual(x) ? "; residual = " * lstrip(terms) : ""
    "$(m.conformations[x.conformation])(" * join(mets, ", ") * res * ")"
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
tag, `_column_parameter` picks a joint-solve column's, and `_onlya_violation_message`
applies it to a group's first edge to label the group.

Returns 1 element for an RE edge (`Krapid`) and 2 elements for an SS edge
(`Kfor`+`Krev`). `name` renders a binding's `Krapid`, plain or fused, as a
dissociation constant and every other edge's — an isomerization, a Theorell–Chance
step, a step with two metabolites on one side — as the equilibrium constant of its
stored direction.
"""
_edge_constants(m::Mechanism, k::Int, state::Symbol) =
    m.edges[k].rapid_equilibrium ? Parameter[Krapid(k, state)] :
                                   Parameter[Kfor(k, state), Krev(k, state)]
