# ABOUTME: The parameter-naming chokepoint name(p, m) and the renderers it builds on:
# ABOUTME: metabolite labels, enzyme-form names, step sides and step text.

# ─── Labels, form names and step sides ────────────────────────────────

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

"""
A step side as parameter names write it for `ctx` (a mechanism or a reaction): its
enzyme form (`_form_label`) followed by its free metabolites (`_met_label`), joined by
`_`, e.g. `"EA_B"` for `E(A) + B`.
"""
_side_label(form::Species, mets::Vector{Metabolite}, ctx) =
    join([_form_label(form, ctx); [_met_label(x, ctx) for x in mets]], "_")

"""
The two sides of `s` in its stored direction, as parameter names write them for `ctx`
(a mechanism or a reaction): `from_species` with the consumed metabolites, then
`to_species` with the released ones. Every step constant is named after them.
"""
_forward_sides(s::Step, ctx) =
    (_side_label(from_species(s), consumed(s), ctx),
     _side_label(to_species(s), released(s), ctx))

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
    terms = join([m.table[i].name for i in entries(x.added)], " + ") *
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

"""What `name(p, m, sides)` reads to name `m`'s constants: `groups`, `m`'s kinetic groups as
`Step`s (`steps`), and `kept`, the `_forward_sides` of each kinetic class's kept member
(`_kept_members`), keyed by `(group, slot, inactive)`. `name(p, m)` computes it for each
call, so a caller that names many constants of one mechanism computes it once and passes
it."""
_class_sides(m::_AnyMech) = (groups = steps(m),
    kept = Dict{Tuple{Int, Int, Bool}, Tuple{String, String}}(
        k => _forward_sides(s, m) for (k, s) in _kept_members(m)))

"""
The sides that name constant `slot` of `step` in `state`: those of the kept member of
its kinetic class, from `sides` (`_class_sides(m)`). An `:EqualAI` step's constant is one
class in both conformations; the inactive copy of an `:OnlyA` step is a class of its own,
named after its own step. A slot the step's group has no class for (the reverse rate of
a rapid-equilibrium step) takes the sides of its slot-1 class, so every parameter of
every step renders.
"""
function _constant_sides(step::Step, slot::Int, state::Symbol, m::_AnyMech, sides)
    g = findfirst(group -> step in group, sides.groups)
    g === nothing && error("Step not found in mechanism: $step")
    tag = m isa AllostericMechanism ? cat_allo_state(m, g) : :None
    state === :I && tag === :OnlyA && return _forward_sides(step, m)
    inactive = state === :I && tag === :NonequalAI
    get(() -> sides.kept[(g, 1, inactive)], sides.kept, (g, slot, inactive))
end

# Every step constant is named after the reaction of its class's kept member: a rate
# constant (`k_`) after the direction it drives, an equilibrium constant (`K_`) after
# the direction whose products-over-reactants it is. A binding's K is named in the
# release direction, so it is a dissociation constant (`K_ES_to_E_S`); an
# isomerization's or other RE step's K in the stored direction (`K_ES_to_EP`). `sides`
# is `_class_sides(m)`.
name(p::Kfor, m::_AnyMech, sides = _class_sides(m)) =
    _render_reaction("k_", _constant_sides(p.step, 1, p.state, m, sides), p.state)
name(p::Krev, m::_AnyMech, sides = _class_sides(m)) =
    _render_reaction("k_", reverse(_constant_sides(p.step, 2, p.state, m, sides)),
                     p.state)
function name(p::Krapid, m::_AnyMech, sides = _class_sides(m))
    s = _constant_sides(p.step, 1, p.state, m, sides)
    _render_reaction("K_", is_binding(p.step) ? reverse(s) : s, p.state)
end

"""
Regulator-site parameter: state tag + ligand name + "reg". No site index —
the AllostericMechanism constructor enforces that each ligand appears at most
once across all sites.
"""
name(p::Kreg, ::AllostericMechanism) =
    Symbol("K_", _state_tag(p.state), String(p.ligand), "reg")

"""
The Parameter(s) governing step `s` with the given allosteric state, by
`is_equilibrium(s)`. `_step_names` applies it to every step of one conformation, with
that conformation's tag, `_column_parameter` picks a joint-solve column's, and
`_onlya_haldane_violation` applies it to a group's first step to label the group.

Returns 1 element for an RE step (`Krapid`) and 2 elements for an SS step
(`Kfor`+`Krev`). `name` renders a binding's `Krapid`, plain or fused, as a
dissociation constant and every other step's — an isomerization, a Theorell–Chance
step, a step with two metabolites on one side — as the equilibrium constant of its
stored direction.
"""
_step_constants(s::Step, state::Symbol) =
    is_equilibrium(s) ? Parameter[Krapid(s, state)] :
    Parameter[Kfor(s, state), Krev(s, state)]
