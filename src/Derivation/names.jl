# ABOUTME: The parameter types, the naming chokepoint _name(p, m), the rules it names by
# ABOUTME: (kept class members, regulator constants) and its renderers.

"""
A constant of a step, and a column of the joint solve (`_constant_columns`): slot `slot` of
step `step` of a mechanism's catalytic graph, in the inactive conformation when `inactive`.
Slot 1 is a rapid-equilibrium step's equilibrium constant or a steady-state step's rate
constant in its stored direction; slot 2 is a steady-state step's rate constant in the
reverse direction.
"""
struct _StepParameter
    step::Int
    slot::Int
    inactive::Bool
end

"""
Regulator-site parameter: the constant of the ligand named `ligand` in the active or
inactive conformation. A ligand sits at one site, so its name and conformation determine
the constant.
"""
struct _RegulatorParameter
    ligand::Symbol
    inactive::Bool
end

"""The allosteric constant `L` of an allosteric mechanism."""
struct _AllostericParameter end

"""A constant of a rate law. A value means nothing without its mechanism, which holds every
other fact about it (`_name(p, m)`). A union rather than an abstract type, so that a loop
over a list of them dispatches without a run-time lookup."""
const _Parameter = Union{_StepParameter, _RegulatorParameter, _AllostericParameter}

# ─── _name(p::_Parameter, m) chokepoint ─────────────────────────────────
#
# Single chokepoint for parameter Symbol production. Renders structural
# names (:K_ES_to_E_S, :k_ES_to_EP, :K_I_ES_to_E_S, :K_Ireg) by dispatch on
# the parameter's type. Routing all parameter-name production through one function
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
    error("_state_tag: unexpected state $state " *
          "(must be one of :None, :EqualAI, :A, :I)")

"""
The constant of the reaction `a → b`: `prefix`, the state tag, then
`"<a>_to_<b>"` (e.g. `:k_ES_to_EP`, `:K_A_ES_to_E_S`).
"""
_render_reaction(prefix::String, (a, b)::Tuple{String, String}, state::Symbol) =
    Symbol(prefix, _state_tag(state), a, "_to_", b)

"""The two sides of step `k` of `m` in its stored direction, which name its constants: each
side is its form's name (`enzyme_forms`), then the name of each metabolite the step takes
up or gives off there, in name order, joined by `_` (`EA_B` for `EA + B`)."""
function _sides(m::Mechanisms.CatalyticScheme, k::Int)
    e = m.steps[k]
    side(f, held) =
        join([String(m.enzyme_forms[f]); sort!(String.(m.metabolites[held]))], "_")
    (side(e.from, e.taken_up), side(e.to, e.given_off))
end

"""
The sides that name constant `slot` of edge `k` of `m`'s catalytic graph in `state`: those
of the kept member of its kinetic class (`kept`, `_kept_members(m)`). An `:EqualAI`
edge's constant is one class in both conformations; the copy of an edge in the
conformation that excludes it (an `:OnlyA` edge's inactive copy, an `:OnlyI` edge's
active one) is a class of its own, named after its own edge.
"""
function _constant_sides(m::Mechanisms.Mechanism, k::Int, slot::Int, state::Symbol,
                         kept::Dict{Tuple{Int, Int, Bool}, Int})
    cm = m.catalytic
    g = cm.steps[k].kinetic_group
    group_tag = _is_allosteric(m) ? m.allosteric.allosteric_states[g] : :None
    group_tag === (state === :I ? :OnlyA : :OnlyI) && return _sides(cm, k)
    inactive = state === :I && group_tag in (:NonequalAI, :OnlyI)
    _sides(cm, kept[(g, slot, inactive)])
end

"""The conformation tag of `p` in `m`, the one its name takes: `:None` in a plain
mechanism; `:I` for an inactive copy; `:EqualAI` or `:A` for an active constant of an
allosteric mechanism."""
_conformation_tag(m::Mechanisms.Mechanism, p::_StepParameter) =
    !_is_allosteric(m) ? :None : p.inactive ? :I :
    m.allosteric.allosteric_states[m.catalytic.steps[p.step].kinetic_group] === :EqualAI ?
    :EqualAI : :A

# Every edge constant is named after the reaction of its class's kept member: a rate
# constant (`k_`) after the direction it drives, an equilibrium constant (`K_`) after
# the direction whose products-over-reactants it is. A binding's K, plain or fused, is
# named in the release direction, so it is a dissociation constant (`K_ES_to_E_S`); every
# other RE step's K, an isomerization, a Theorell–Chance step or a step with two
# metabolites on one side, in the stored direction (`K_ES_to_EP`). `kept` is
# `_kept_members(m)`, which a caller that names many constants of one mechanism computes
# once and passes.
function _name(p::_StepParameter, m::Mechanisms.Mechanism, kept = _kept_members(m))
    cm = m.catalytic
    state = _conformation_tag(m, p)
    s = _constant_sides(m, p.step, p.slot, state, kept)
    e = cm.steps[p.step]
    if e.rapid_equilibrium
        _render_reaction("K_", _constant_named_in_reverse(e) ? reverse(s) : s, state)
    else
        _render_reaction("k_", p.slot == 1 ? s : reverse(s), state)
    end
end

"""
Regulator-site parameter: state tag + ligand name + "reg". No site index —
the Mechanism constructor enforces that each ligand appears at most
once across all sites. An `:EqualAI` ligand has one bare name in both conformations
(`K_<ligand>reg`); any other ligand takes its conformation's `K_A_` or `K_I_` name.
"""
function _name(p::_RegulatorParameter, am::Mechanisms.Mechanism)
    tag = only(t for site in am.allosteric.regulatory_sites
               for (lig, t) in zip(site.ligands, site.allosteric_states)
               if lig === p.ligand)
    state = tag === :EqualAI ? :EqualAI : p.inactive ? :I : :A
    Symbol("K_", _state_tag(state), String(p.ligand), "reg")
end

_name(::_AllostericParameter, ::Mechanisms.Mechanism) = :L

# ─── Kept members of kinetic classes ─────────────────────────────────

"""
The priority key of constant `slot` of step `k` of `m` in conformation block `block` (2:
a copy outside its conformation's law, 1: an inactive copy in the law, 0: an active copy
in the law), read from the graph alone. The solve eliminates the larger key first, so a
column with a larger key becomes dependent sooner. In turn: the block; the step's shape
(3 an isomerization, which takes up and gives off nothing; 2 a step that takes up or
gives off more than one metabolite; 1 a binding); whether neither of its forms is free E;
steady state before rapid equilibrium; a steady-state constant that runs against the
reaction (for a step stored against the reaction, the forward rate: as stored, it takes
up fewer substrates and gives off fewer products than it gives off substrates and takes
up products, as a product binding or a release of products, stored as its uptake, does;
for every other step its reverse rate, since the constructor stores it the way the
reaction runs, or the reaction gives it no direction); then its two forms' indices, the
smaller first, and its metabolites' indices, sorted; last its stored place `k`. The
indices follow the forms' and metabolites' names (`CatalyticScheme`), so the key does not
change with the order in which a writing gives them. The key decides which constants are
fitted, never the law.
"""
function _priority_key(m::Mechanisms.CatalyticScheme, k::Int, slot::Int, block::Int)
    e = m.steps[k]
    held = length(e.taken_up) + length(e.given_off)
    shape = held == 0 ? 3 : held == 1 ? 1 : 2
    off_free_E = e.from != 1 && e.to != 1
    ss = !e.rapid_equilibrium
    net(list) = count(in(m.substrates), list) - count(in(m.products), list)
    stored_against = net(e.taken_up) - net(e.given_off) < 0
    against = ss && (slot == 2) != stored_against
    # `k` decides only between two steps with the same forms and metabolites.
    first_form, last_form = minmax(e.from, e.to)
    (block, shape, off_free_E, ss, against, first_form, last_form,
     sort!([e.taken_up; e.given_off]), k)
end

"""Each kinetic group's allosteric state, `:None` for a plain mechanism, whose groups are
numbered 1 to the largest `group` of its steps."""
_solve_tags(m::Mechanisms.Mechanism) =
    _is_allosteric(m) ? m.allosteric.allosteric_states :
    fill(:None, maximum(e -> e.kinetic_group, m.catalytic.steps; init = 0))

"""Whether the active column of edge `k` of `cm`, in a group tagged `tags`, enters the
law: its edge lies in the active conformation's law (`alaw`), or it is an `:EqualAI`
edge of the inactive conformation's law (`law`), whose one constant the inactive law
reads."""
_in_active_law(cm::Mechanisms.CatalyticScheme, k::Int, tags::Vector{Symbol}, law::BitVector,
               alaw::BitVector) =
    alaw[k] || (tags[cm.steps[k].kinetic_group] === :EqualAI && law[k])

"""Whether the steps of a group with allosteric state `tag` take their own constants in the
inactive conformation: those of a `:NonequalAI`, `:OnlyA` or `:OnlyI` group do; an
`:EqualAI` group's read the active ones."""
_has_inactive_columns(tag::Symbol) = tag in (:NonequalAI, :OnlyA, :OnlyI)

"""
    _ordered_columns(m) -> (; columns, flat, tags, law, alaw)

The joint solve's columns of `m` in elimination order, which `_joint_system` and
`_kept_members` share: one column (`_StepParameter`) per slot of every edge, and for an
allosteric mechanism one more per slot of every edge of a `:NonequalAI`, `:OnlyA` or
`:OnlyI` group, for its inactive copy, ordered by `_priority_key`, larger first. `flat`
holds the edges group by group (`_flat_edges`), `tags` each group's tag (`_solve_tags`),
`law[k]` whether edge `k` lies in the inactive conformation's law and `alaw[k]` whether it
lies in the active conformation's law (`steps_in_active_conformation`; `_in_active_law`
adds the `:EqualAI` constants the inactive law reads). A plain mechanism's one
conformation runs every step, so both masks hold every edge.
"""
function _ordered_columns(m::Mechanisms.Mechanism)
    cm = m.catalytic
    flat = _flat_edges(cm)
    tags = _solve_tags(m)
    allosteric = _is_allosteric(m)
    law = Mechanisms.steps_in_inactive_conformation(m)
    alaw = Mechanisms.steps_in_active_conformation(m)
    columns = _StepParameter[]
    for k in flat, inactive in (false, true)
        inactive &&
            !(allosteric && _has_inactive_columns(tags[cm.steps[k].kinetic_group])) &&
            continue
        for slot in 1:(cm.steps[k].rapid_equilibrium ? 1 : 2)
            push!(columns, _StepParameter(k, slot, inactive))
        end
    end
    block(c) = c.inactive ? (law[c.step] ? 1 : 2) :
               _in_active_law(cm, c.step, tags, law, alaw) ? 0 : 2
    priorities = [_priority_key(cm, c.step, c.slot, block(c)) for c in columns]
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
function _kept_members(m::Mechanisms.Mechanism, cols = _ordered_columns(m))
    cm = m.catalytic
    kept = Dict{Tuple{Int, Int, Bool}, Int}()
    for c in cols.columns
        g = cm.steps[c.step].kinetic_group
        tag = cols.tags[g]
        (c.inactive ? !(tag in (:NonequalAI, :OnlyI)) : tag === :OnlyI) && continue
        kept[(g, c.slot, c.inactive)] = c.step
    end
    kept
end

"""
    _constant_columns(m) -> (; columns, names, in_law, kept_in_law, flat, tags, law)

The constants of the steps of `m` as the joint solve's columns, in elimination order
(`_ordered_columns`), with the name each column takes. `columns[i]` is a `_StepParameter`:
slot `slot` of step `step` of the catalytic graph, in the inactive conformation when
`inactive`. A plain mechanism has one column per slot of every step; an allosteric
mechanism also has one per slot of every step of a `:NonequalAI`, `:OnlyA` or `:OnlyI`
group, for its inactive copy.

`names[i]` is the name of the column's kinetic class, that of its kept member
(`_kept_members`, through the `_name(p, m)` chokepoint); `in_law[i]` whether the column
enters the law: an inactive copy when its step lies in the inactive conformation's law,
an active column when its step lies in the active one's or is an `:EqualAI` step of the
inactive one's, whose one constant the inactive law reads; `kept_in_law[i]` whether it
enters the law as its class's kept member, the last of its class (`_class`) in
elimination order, so that it stands for the class. `flat` holds the steps group by
group, `tags` each group's tag (`:None` for a plain mechanism) and `law` the steps of the
inactive conformation (`steps_in_inactive_conformation`; every step for a plain mechanism,
whose one conformation runs them all).
"""
function _constant_columns(m::Mechanisms.Mechanism)
    cols = _ordered_columns(m)
    kept = _kept_members(m, cols)
    cm = m.catalytic
    names = [_name(c, m, kept) for c in cols.columns]
    in_law = BitVector([c.inactive ? cols.law[c.step] :
                        _in_active_law(cm, c.step, cols.tags, cols.law, cols.alaw)
                        for c in cols.columns])
    excluded = (_excluded_steps(m, 1), _excluded_steps(m, 2))
    last_of = Dict(_class(m, excluded, c) => i for (i, c) in enumerate(cols.columns))
    kept_in_law = BitVector([in_law[i] && last_of[_class(m, excluded, c)] == i
                             for (i, c) in enumerate(cols.columns)])
    (; cols.columns, names, in_law, kept_in_law, cols.flat, cols.tags, cols.law)
end

"""
    _shared_constant_names(m) -> Vector{Symbol}

The names that two constants of the rate law of `m` share, regulator constants included,
each once; empty when each constant has a name of its own. The law would read two
constants with one name as one, so the guard rejects such a mechanism. Form names are
opaque, so a form's can imitate a side's: the side `ES + I` and a form named `ES_I` both
render `ES_I`, and the constructors accept a step between them.
"""
function _shared_constant_names(m::Mechanisms.Mechanism)
    cols = _constant_columns(m)
    names = [cols.names[cols.kept_in_law]; _regulator_names(m)]
    unique!([x for x in names if count(==(x), names) > 1])
end

# ─── Regulator constants ────────────────────────────────────

"""
    _regulator_constant(am, lig, tag, inactive) -> Union{Symbol, Nothing}

The dissociation-constant name of the ligand named `lig` (allosteric state `tag`) of `am`
in the active (`inactive = false`) or inactive conformation, or `nothing` when the
ligand is absent from that conformation (`:OnlyI` in the active state, `:OnlyA` in the
inactive one). An `:EqualAI` ligand has one bare name in both conformations
(`K_<ligand>reg`); an `:OnlyA`, `:OnlyI` or `:NonequalAI` ligand has its conformation's
`K_A_` or `K_I_` name. Renders the name via the `_name(::_RegulatorParameter, am)`
chokepoint.
"""
function _regulator_constant(am::Mechanisms.Mechanism, lig::Symbol, tag::Symbol,
                             inactive::Bool)
    tag === (inactive ? :OnlyA : :OnlyI) && return nothing
    _name(_RegulatorParameter(lig, inactive), am)
end

"""The steps group by group, each group in step order: the order the joint solve and the
derivation read them in."""
_flat_edges(m::Mechanisms.CatalyticScheme) = reduce(vcat, _kinetic_groups(m); init = Int[])
