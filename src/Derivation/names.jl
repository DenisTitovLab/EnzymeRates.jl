# ABOUTME: The parameter-naming chokepoint name(p, m), the rules it names by (kept class
# ABOUTME: members, regulator constants) and its renderers.

# Parameter family.
abstract type Parameter end

# Edge parameters: a rapid-equilibrium edge's equilibrium constant, and a steady-state
# edge's forward and reverse rate constants, `edge` indexing the edges of a mechanism's
# catalytic graph.
struct Krapid <: Parameter; edge::Int; state::Symbol end
struct Kfor   <: Parameter; edge::Int; state::Symbol end
struct Krev   <: Parameter; edge::Int; state::Symbol end

"""
Regulator-site parameter: the constant of the ligand named `ligand` in the A (active) or
I (inactive) branch of the polynomial (`state`). A ligand sits at one site, so its name
and state determine the constant.
"""
struct Kreg <: Parameter
    ligand::Symbol
    state::Symbol
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
edge's constant is one class in both conformations; the copy of an edge in the
conformation that excludes it (an `:OnlyA` edge's inactive copy, an `:OnlyI` edge's
active one) is a class of its own, named after its own edge. A slot the edge's group has no
class for (the reverse rate of a rapid-equilibrium edge) takes the sides of its slot-1
class, so every parameter of every edge renders.
"""
function _constant_sides(m::_AnyMech, k::Int, slot::Int, state::Symbol,
                         kept::Dict{Tuple{Int, Int, Bool}, Int})
    cm = Mechanisms.catalytic(m)
    g = cm.steps[k].group
    group_tag = m isa Mechanisms.AllostericMechanism ? m.mwc_tags[g] : :None
    Mechanisms.excludes(group_tag, state === :I ? :inactive : :active) &&
        return cm.step_sides[k]
    inactive = state === :I && group_tag in (:NonequalAI, :OnlyI)
    cm.step_sides[get(() -> kept[(g, 1, inactive)], kept, (g, slot, inactive))]
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
    binds = Mechanisms.catalytic(m).step_kinds[p.edge] in
            (:binding, :fused_binding)
    _render_reaction("K_", binds ? reverse(s) : s, p.state)
end

"""
Regulator-site parameter: state tag + ligand name + "reg". No site index —
the AllostericMechanism constructor enforces that each ligand appears at most
once across all sites.
"""
name(p::Kreg, ::Mechanisms.AllostericMechanism) =
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
_edge_constants(m::Mechanisms.Mechanism, k::Int, state::Symbol) =
    m.steps[k].rapid_equilibrium ? Parameter[Krapid(k, state)] :
                                   Parameter[Kfor(k, state), Krev(k, state)]

# ─── Kept members of kinetic classes ─────────────────────────────────

"""
The priority key of constant `slot` of edge `k` of `m` in conformation block `block` (2:
a copy outside its conformation's law, 1: an inactive copy in the law, 0: an active copy
in the law). The solve eliminates the larger key first, so a column with a larger key
becomes dependent sooner. In turn: the block; the edge's kind (3 an isomerization, 2
other chemistry, 1 a plain binding); whether both of its forms hold a metabolite; steady
state before rapid equilibrium; a steady-state constant that runs against the reaction
(for a product binding its uptake, for every other edge its reverse rate: the constructor
stores each edge in its physical-forward direction); then the edge's content as stored,
its rank `edge_rank[k]` (`_edge_base_rank`), the smaller rank kept.
"""
function _priority_key(m::Mechanisms.Mechanism, k::Int, slot::Int, block::Int,
                       edge_rank::Vector{Int})
    e = m.steps[k]
    sk = m.step_kinds[k]
    kind = sk === :isomerization ? 3 : sk === :binding ? 1 : 2
    occupied = !isempty(m.nodes[e.from].bound) && !isempty(m.nodes[e.to].bound)
    ss = !e.rapid_equilibrium
    against = ss && (slot == 2) != _binds_product(m, k)
    (block, kind, occupied, ss, against, edge_rank[k])
end

"""Each kinetic group's allosteric tag, `:None` for a plain mechanism."""
_solve_tags(m::Mechanisms.Mechanism) = fill(:None, length(m.groups))
_solve_tags(am::Mechanisms.AllostericMechanism) = am.mwc_tags

"""Whether the active column of edge `k` of `cm`, in a group tagged `tags`, enters the
law: its edge lies in the active conformation's law (`alaw`), or it is an `:EqualAI`
edge of the inactive conformation's law (`law`), whose one constant the inactive law
reads."""
_in_active_law(cm::Mechanisms.Mechanism, k::Int, tags::Vector{Symbol}, law::BitVector,
               alaw::BitVector) =
    alaw[k] || (tags[cm.steps[k].group] === :EqualAI && law[k])

"""Whether the steps of a group with allosteric tag `tag` take their own constants in the
inactive conformation: those of a `:NonequalAI`, `:OnlyA` or `:OnlyI` group do; an
`:EqualAI` group's read the active ones."""
_has_inactive_columns(tag::Symbol) = tag in (:NonequalAI, :OnlyA, :OnlyI)

"""
    _ordered_columns(m) -> (; columns, flat, tags, law, alaw)

The joint solve's columns of `m` in elimination order, which `_joint_system` and
`_kept_members` share: one column (`_Column`) per slot of every edge, and for an
allosteric mechanism one more per slot of every edge of a `:NonequalAI`, `:OnlyA` or
`:OnlyI` group, for its inactive copy, ordered by `_priority_key`, larger first. `flat`
holds the edges group by group (`_flat_edges`), `tags` each group's tag (`_solve_tags`),
`law[k]` whether edge `k` lies in the inactive conformation's law and `alaw[k]` whether it
lies in the active conformation's law (`active_steps`; `_in_active_law` adds the
`:EqualAI` constants the inactive law reads).
"""
function _ordered_columns(m::_AnyMech)
    cm = Mechanisms.catalytic(m)
    flat = _flat_edges(cm)
    tags = _solve_tags(m)
    allosteric = m isa Mechanisms.AllostericMechanism
    law = _conformation_steps(m, 2)
    alaw = _conformation_steps(m, 1)
    columns = _Column[]
    for k in flat, inactive in (false, true)
        inactive && !(allosteric && _has_inactive_columns(tags[cm.steps[k].group])) &&
            continue
        for slot in 1:(cm.steps[k].rapid_equilibrium ? 1 : 2)
            push!(columns, _Column((k, slot, inactive)))
        end
    end
    rank = _edge_base_rank(cm)
    block(c) = c.inactive ? (law[c.edge] ? 1 : 2) :
               _in_active_law(cm, c.edge, tags, law, alaw) ? 0 : 2
    priorities = [_priority_key(cm, c.edge, c.slot, block(c), rank) for c in columns]
    (; columns = columns[sortperm(priorities; rev = true)], flat, tags, law, alaw)
end

"""
The kept member of each kinetic class of `m`, keyed by `(group, slot, inactive)`: the
edge whose column comes last in elimination order (`cols`, `_ordered_columns(m)`), its
smallest `_priority_key`, so the solve keeps it whenever it keeps one of the class (E3).
A `:NonequalAI` or `:OnlyI` group has an inactive class, kept by a member in the inactive
law when it has one; an `:EqualAI` group's edges share their active class; the copy of an
edge in the conformation that excludes it (an `:OnlyA` edge's inactive copy, an `:OnlyI`
edge's active one) is a class of its own, so neither of these appears.
"""
function _kept_members(m::_AnyMech, cols = _ordered_columns(m))
    cm = Mechanisms.catalytic(m)
    kept = Dict{Tuple{Int, Int, Bool}, Int}()
    for c in cols.columns
        g = cm.steps[c.edge].group
        tag = cols.tags[g]
        (c.inactive ? !(tag in (:NonequalAI, :OnlyI)) : tag === :OnlyI) && continue
        kept[(g, c.slot, c.inactive)] = c.edge
    end
    kept
end

"""The parameter of column `c` of `m`: slot `c.slot` of its edge's `_edge_constants`, in
the state its name takes: `:None` in a plain mechanism; `:I` for an inactive copy;
`:EqualAI` or `:A` for an active constant of an allosteric mechanism."""
function _column_parameter(m::_AnyMech, c::_Column)
    cm = Mechanisms.catalytic(m)
    state = m isa Mechanisms.Mechanism ? :None : c.inactive ? :I :
            m.mwc_tags[cm.steps[c.edge].group] === :EqualAI ? :EqualAI : :A
    _edge_constants(cm, c.edge, state)[c.slot]
end

"""
    constant_columns(m) -> (; columns, names, in_law, flat, tags, law)

The constants of the steps of `m` as the joint solve's columns, in elimination order
(`_ordered_columns`), with the name each column takes. `columns[i]` is the column
`(; edge, slot, inactive)`: slot `slot` of step `edge` of the catalytic graph, in the
inactive conformation when `inactive`. Slot 1 is a rapid-equilibrium step's equilibrium
constant or a steady-state step's rate constant in its stored direction; slot 2 is a
steady-state step's rate constant in the reverse direction. A plain mechanism has one
column per slot of every step; an allosteric mechanism also has one per slot of every
step of a `:NonequalAI`, `:OnlyA` or `:OnlyI` group, for its inactive copy.

`names[i]` is the name of the column's kinetic class, that of its kept member
(`_kept_members`, through the `name(p, m)` chokepoint); `in_law[i]` whether the column
enters the law: an inactive copy when its step lies in the inactive conformation's law,
an active column when its step lies in the active one's or is an `:EqualAI` step of the
inactive one's, whose one constant the inactive law reads. `flat` holds the steps group
by group, `tags` each group's tag (`:None` for a plain mechanism) and `law` the steps of
the inactive conformation (`inactive_steps`).
"""
function constant_columns(m::_AnyMech)
    cols = _ordered_columns(m)
    kept = _kept_members(m, cols)
    cm = Mechanisms.catalytic(m)
    names = [name(_column_parameter(m, c), m, kept) for c in cols.columns]
    in_law = BitVector([c.inactive ? cols.law[c.edge] :
                        _in_active_law(cm, c.edge, cols.tags, cols.law, cols.alaw)
                        for c in cols.columns])
    (; cols.columns, names, in_law, cols.flat, cols.tags, cols.law)
end

# ─── Regulator constants ────────────────────────────────────

"""
    regulator_constant(am, lig, tag, inactive) -> Union{Symbol, Nothing}

The dissociation-constant name of the ligand named `lig` (allosteric tag `tag`) of `am`
in the active (`inactive = false`) or inactive conformation, or `nothing` when the
ligand is absent from that conformation (`:OnlyI` in the active state, `:OnlyA` in the
inactive one). An `:EqualAI` ligand has one bare name in both conformations
(`K_<ligand>reg`); an `:OnlyA`, `:OnlyI` or `:NonequalAI` ligand has its conformation's
`K_A_` or `K_I_` name. Renders the name via the `name(::Kreg, am)` chokepoint.
"""
function regulator_constant(am::Mechanisms.AllostericMechanism, lig::Symbol, tag::Symbol,
                            inactive::Bool)
    Mechanisms.excludes(tag, inactive ? :inactive : :active) && return nothing
    name(Kreg(lig, tag === :EqualAI ? :EqualAI : inactive ? :I : :A), am)
end

"""The steps group by group, each group in step order: the order the joint solve and the
derivation read them in."""
_flat_edges(m::Mechanisms.Mechanism) = reduce(vcat, m.groups; init = Int[])
"""Whether step `k` of `m` binds a product: it takes up one product, not as a modifier."""
_binds_product(m::Mechanisms.Mechanism, k::Int) =
    m.step_kinds[k] in (:binding, :fused_binding) &&
    (x = only(m.steps[k].taken_up); !x.modifier && insorted(x.name, m.products))
"""Each node's position in content order (`isless(::EnzymeForm)`: label, bound, added,
subtracted): its index, except that free E moves behind the nodes whose label sorts before
`E` (a free node labelled D)."""
_base_rank(m::Mechanisms.Mechanism) = invperm(sortperm(m.nodes))
"""Each step's position by the base ranks (`_base_rank`) of its two nodes, then its
taken-up and given-off metabolites as lists: the content order of the joint solve's ties
(`_priority_key`)."""
function _edge_base_rank(m::Mechanisms.Mechanism)
    rank = _base_rank(m)
    invperm(sortperm([Mechanisms.Step(rank[e.from], rank[e.to], e.taken_up, e.given_off,
                           e.rapid_equilibrium, e.group) for e in m.steps]))
end
