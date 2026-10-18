# ABOUTME: The parameter-naming chokepoint name(p, m) and the renderers it builds on:
# ABOUTME: metabolite labels, enzyme-form names, step sides and step text.

# ─── Labels, form names and step sides ────────────────────────────────

"""Whether `x` is the competitive-inhibitor copy of a substrate or product of `ctx` (a
mechanism, or the reaction of one under construction; `_is_reactant`): a reactant bound
as a modifier."""
_is_reactant_copy(x::Metabolite, ctx) =
    x isa CompetitiveInhibitor && _is_reactant(ctx, name(x))

"""
A metabolite as written in parameter names: its name, with suffix `mod` on a
competitive-inhibitor copy of a substrate or product of `ctx` (a mechanism, or the
reaction of one under construction; `_is_reactant_copy`), so the copy's forms and
constants stay distinct from those of the reactant itself. A competitive inhibitor that
is no reactant binds only as a modifier, so its name stands bare.
"""
_met_label(x::Metabolite, ctx) =
    _is_reactant_copy(x, ctx) ? String(name(x)) * "mod" : String(name(x))

"""
An enzyme form of `ctx` (a mechanism, or the reaction of one under construction) as
written in parameter names: its conformation label, then the labels (`_met_label`) of its
bound metabolites without a separator, then, for a residual, `_res`, `_plus<name>` per
substrate it adds and `_minus<name>` per product it subtracts, each in name order: `E`,
`EAB`, `EPmodQ`, `EB_res_plusA_minusP`.
"""
function _form_label(sp::Species, ctx)
    head = String(conformation(sp)) * join(_met_label(x, ctx) for x in bound(sp))
    r = residual(sp)
    isempty(r) && return head
    head * "_res" * join("_plus" * String(name(a)) for a in added(r)) *
        join("_minus" * String(name(p)) for p in subtracted(r))
end

"""
A step side of `ctx` (a mechanism, or the reaction of one under construction) as
written in parameter names: its enzyme form (`_form_label`) followed by its free
metabolites (`_met_label`), joined by `_`, e.g. `"EA_B"` for `E(A) + B`.
"""
_side_label(form::Species, mets::Vector{Metabolite}, ctx) =
    join([_form_label(form, ctx); [_met_label(x, ctx) for x in mets]], "_")

"""
The two sides of `s` in its stored direction, as parameter names write them for `ctx`
(a mechanism, or the reaction of one under construction): `from_species` with the
consumed metabolites, then `to_species` with the released ones. Every step constant is
named after them.
"""
_forward_sides(s::Step, ctx) =
    (_side_label(from_species(s), consumed(s), ctx),
     _side_label(to_species(s), released(s), ctx))

"""A step of `ctx` (a mechanism, or the reaction of one under construction) as text: its
entry side, an arrow (`⇌` rapid equilibrium, `<-->` steady state) and its exit side. A
side lists the enzyme form as parameter names write it (`_form_label`), then the
metabolites the step consumes or releases."""
_step_text(s::Step, ctx) =
    join((_form_label(from_species(s), ctx), name.(consumed(s))...), " + ") *
    (is_equilibrium(s) ? " ⇌ " : " <--> ") *
    join((_form_label(to_species(s), ctx), name.(released(s))...), " + ")

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

"""The `_forward_sides` of each kinetic class's kept member (`_kept_members`), keyed by
`(group, slot, inactive)`. `name(p, m)` computes them for each call, so a caller that
names many constants of one mechanism computes them once and passes them."""
_class_sides(m::_AnyMech) = Dict{Tuple{Int, Int, Bool}, Tuple{String, String}}(
    k => _forward_sides(s, m) for (k, s) in _kept_members(m))

"""
The sides that name constant `slot` of `step` in `state`: those of the kept member of
its kinetic class, from `sides` (`_class_sides(m)`). An `:EqualAI` step's constant is one
class in both conformations; the inactive copy of an `:OnlyA` step is a class of its own,
named after its own step. A slot the step's group has no class for (the reverse rate of
a rapid-equilibrium step) takes the sides of its slot-1 class, so every parameter of
every step renders.
"""
function _constant_sides(step::Step, slot::Int, state::Symbol, m::_AnyMech, sides)
    g = findfirst(group -> step in group, steps(m))
    g === nothing && error("Step not found in mechanism: $step")
    tag = m isa AllostericMechanism ? cat_allo_state(m, g) : :None
    state === :I && tag === :OnlyA && return _forward_sides(step, m)
    inactive = state === :I && tag === :NonequalAI
    get(() -> sides[(g, 1, inactive)], sides, (g, slot, inactive))
end

# Every step constant is named after the reaction of its class's kept member: a rate
# constant (`k_`) after the direction it drives, an equilibrium constant (`K_`) after
# the direction whose products-over-reactants it is. A binding's K is named in the
# release direction, so it is a dissociation constant (`K_ES_to_E_S`); an
# isomerization's or other RE step's K in the stored direction (`K_ES_to_EP`). `sides`
# holds the kept members' sides (`_class_sides`).
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
