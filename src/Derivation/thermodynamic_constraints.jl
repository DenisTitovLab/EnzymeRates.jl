# ABOUTME: A mechanism's constants under the joint thermodynamic solve: their types and
# ABOUTME: names, the solve, parameters, their units and the dependent constants' values.

# ─── The parameter types ────────────────────────────────────────────────

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
# Every parameter name (:K_ES_to_E_S, :k_ES_to_EP, :K_I_ES_to_E_S, :K_Ireg) comes from
# `_name`, by dispatch on the parameter's type, so a change of the naming scheme edits one
# function.

"""
The text a constant's name takes after its type prefix (`K_`, `k_`) for conformation tag
`state` (`_conformation_tag`): `A_` for `:A`, an allosteric mechanism's active constant;
`I_` for `:I`, an inactive one; none for `:EqualAI`, one constant of both conformations,
or for `:None`, a plain mechanism's.
"""
_state_tag(state::Symbol) = state === :A ? "A_" : state === :I ? "I_" : ""

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
The name of regulator constant `p` of `am`: `K_`, the conformation's tag, the ligand's
name and `reg`. An `:EqualAI` ligand has one name in both conformations (`K_<ligand>reg`);
any other ligand takes its conformation's `K_A_` or `K_I_` name. The name holds no site,
since a ligand sits at one site (`_RegulatorParameter`).
"""
function _name(p::_RegulatorParameter, am::Mechanisms.Mechanism)
    tag = only(t for site in am.allosteric.regulatory_sites
               for (lig, t) in zip(site.ligands, site.allosteric_states)
               if lig === p.ligand)
    state = tag === :EqualAI ? :EqualAI : p.inactive ? :I : :A
    Symbol("K_", _state_tag(state), String(p.ligand), "reg")
end

_name(::_AllostericParameter, ::Mechanisms.Mechanism) = :L

# ─── Kept members of kinetic classes ────────────────────────────────────

"""
The priority key of constant `slot` of step `k` of `m` in conformation block `block`, read
from the graph alone. The solve eliminates the larger key first, so a column with a larger
key becomes dependent sooner. The key decides which constants are fitted, never the law.
Its entries, in turn:
1. the block: 2 for a copy outside its conformation's law, 1 for an inactive copy in the
   law, 0 for an active copy in the law;
2. the step's shape: 3 for an isomerization, which takes up and gives off nothing; 2 for
   a step that takes up or gives off more than one metabolite; 1 for a binding;
3. whether neither of its forms is free E;
4. steady state before rapid equilibrium;
5. whether it is a steady-state constant that runs against the reaction: the forward rate
   of a step stored against the reaction, the reverse rate of any other step;
6. its two forms' indices, the smaller first, then its metabolites' indices, sorted;
7. its stored place `k`.
A step is stored against the reaction when, as stored, it takes up fewer substrates and
gives off fewer products than it gives off substrates and takes up products, as a product
binding or a release of products stored as its uptake does. The constructor stores every
other step the way the reaction runs, or the reaction gives it no direction. The indices
follow the forms' and metabolites' names (`CatalyticScheme`), so the key does not change
with the order in which a writing gives them.
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
lies in the active conformation's law (`steps_running_in_allosteric_conformations`;
`_in_active_law` adds the `:EqualAI` constants the inactive law reads). A plain
mechanism's one conformation runs every step, so both masks hold every edge.
"""
function _ordered_columns(m::Mechanisms.Mechanism)
    cm = m.catalytic
    flat = _flat_edges(cm)
    tags = _solve_tags(m)
    allosteric = _is_allosteric(m)
    runs = Mechanisms.steps_running_in_allosteric_conformations(m)
    law, alaw = runs.inactive, runs.active
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
edge whose column comes last in elimination order (`cols`, `_ordered_columns(m)`), the
one of smallest `_priority_key`. The solve keeps it whenever it keeps one of the class;
the tests call this rule E3. A `:NonequalAI` or `:OnlyI` group has an inactive class, kept
by a member in the inactive law when it has one; an `:EqualAI` group's edges share their
active class. The copy of an edge in the conformation that excludes it (an `:OnlyA` edge's
inactive copy, an `:OnlyI` edge's active one) is a class of its own, so neither of these
appears.
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
    _constant_columns(m) -> (; columns, names, in_law, kept_in_law, flat)

The constants of the steps of `m` as the joint solve's columns and `flat`, as
`_ordered_columns` gives them, with:
- `names[i]`, the name of column `i`'s kinetic class, that of its kept member
  (`_kept_members`, through the `_name(p, m)` chokepoint);
- `in_law[i]`, whether the column enters the law: an inactive copy when its step lies in
  the inactive conformation's law, an active column when `_in_active_law` holds;
- `kept_in_law[i]`, whether it enters the law as its class's kept member, the last of its
  class (`_class`) in elimination order, so that it stands for the class.
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
    (; cols.columns, names, in_law, kept_in_law, cols.flat)
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

# ─── Regulator constants ────────────────────────────────────────────────

"""
    _regulator_constant(am, lig, tag, inactive) -> Union{Symbol, Nothing}

The name of the dissociation constant (`_name`) of the ligand named `lig` (allosteric
state `tag`) of `am` in the active (`inactive = false`) or inactive conformation, or
`nothing` when the ligand is absent from that conformation (`:OnlyI` in the active state,
`:OnlyA` in the inactive one).
"""
function _regulator_constant(am::Mechanisms.Mechanism, lig::Symbol, tag::Symbol,
                             inactive::Bool)
    tag === (inactive ? :OnlyA : :OnlyI) && return nothing
    _name(_RegulatorParameter(lig, inactive), am)
end

"""The steps group by group, each group in step order: the order the joint solve and the
derivation read them in."""
_flat_edges(m::Mechanisms.CatalyticScheme) = reduce(vcat, _kinetic_groups(m); init = Int[])

# ─── Parameters API ─────────────────────────────────────────────────────

"""
    parameters(m) -> Tuple{Vararg{Symbol}}
    parameters(rate) -> Tuple{Vararg{Symbol}}

The names of the constants of the rate law of `m` that a fit estimates, as a tuple of
Symbols: the independent constants in name order, then an allosteric mechanism's regulator
constants and `L`. These are the names `params` holds in a call of
[`compile_rate_equation`](@ref)`(m)`, whose result `rate` gives the same names as a
compile-time constant. The equilibrium constant `Keq` is not among them, since it is
given, never fitted. The rate law computes every other constant from these and `Keq` (see
[`dependent_parameter_values`](@ref)). Errors when the rate law of `m` is not defined
(`check_valid_rate_equation`).
"""
function parameters(m::Mechanisms.Mechanism)
    check_valid_rate_equation(m)
    Tuple(_memo_solve(m).fitted)
end

"""
    parameter_units(m) -> (; rate::Vector{Int}, concentrations::Matrix{Int})

The units of the constants of [`parameters`](@ref)`(m)`, as integer powers: constant
`parameters(m)[i]` has the units of the rate per unit enzyme to the power `rate[i]` times
the concentration of each metabolite `metabolites(m)[j]` to the power
`concentrations[i, j]`. A steady-state step's rate constant has `rate` 1 and −1 for each
metabolite the step takes up in its direction. Every other constant has `rate` 0: a
rapid-equilibrium binding's dissociation constant +1 for the metabolite that binds, any
other rapid-equilibrium constant +1 for each metabolite given off and −1 for each taken
up, a regulator's constant +1 for its ligand, and `L` none. Errors when the rate law of
`m` is not defined (`check_valid_rate_equation`).
"""
function parameter_units(m::Mechanisms.Mechanism)
    check_valid_rate_equation(m)
    info = _memo_solve(m).info
    column = Dict(x => j for (j, x) in enumerate(Mechanisms.metabolites(m)))
    units = (; rate = zeros(Int, length(info)),
             concentrations = zeros(Int, length(info), length(column)))
    for (i, p) in enumerate(info)
        _write_units!(units, i, p, m, column)
    end
    units
end

"""Writes the units of `p`, the `i`-th fitted constant of `m`, into row `i` of `units`
(`parameter_units`), whose rows start at zero; `column` maps each metabolite of `m` to its
column."""
function _write_units!(units, i::Int, p::_StepParameter, m::Mechanisms.Mechanism,
                       column::Dict{Symbol, Int})
    c = m.catalytic
    e = c.steps[p.step]
    add!(xs, power) =
        foreach(x -> units.concentrations[i, column[c.metabolites[x]]] += power, xs)
    if e.rapid_equilibrium
        # An equilibrium constant is what the direction it is named after gives off over
        # what it takes up: a binding's dissociation constant has what binds.
        named = _constant_named_in_reverse(e) ? -1 : 1
        add!(e.given_off, named)
        add!(e.taken_up, -named)
    else
        units.rate[i] = 1
        add!(p.slot == 1 ? e.taken_up : e.given_off, -1)
    end
    nothing
end
_write_units!(units, i::Int, p::_RegulatorParameter, m::Mechanisms.Mechanism,
              column::Dict{Symbol, Int}) =
    (units.concentrations[i, column[p.ligand]] += 1; nothing)
_write_units!(units, i::Int, ::_AllostericParameter, m::Mechanisms.Mechanism,
              column::Dict{Symbol, Int}) = nothing

"""The dependent constants of the joint solve `s` in name order, less the kinetic-group
tie pins, which are not constants of the rate law."""
_dependents(s) =
    sort!(filter(d -> d.label !== :kinetic_group, s.dependents); by = d -> s.names[d.col])

"""
    dependent_parameter_values(m, params::NamedTuple, Keq::Real) -> NamedTuple

The values of the dependent constants of `m`'s rate law, by name in name order, each
computed from the fitted constants in `params` and the equilibrium constant `Keq` by its
Haldane or Wegscheider relation (the relations [`rate_equation_string`](@ref) prints).
`params` holds the values of [`parameters`](@ref)`(m)` by name; other names are ignored. An
error names every constant of `parameters(m)` that `params` lacks.
`merge(params, dependent_parameter_values(m, params, Keq))` is every constant.

```jldoctest
julia> using EnzymeRates

julia> m = @enzyme_mechanism begin
           substrates: S
           products: P
           steps: begin
               E + S ⇌ ES
               ES <--> EP
               EP ⇌ E + P
           end
       end;

julia> p = (K_EP_to_E_P = 1.0, K_ES_to_E_S = 2.0, k_ES_to_EP = 3.0);

julia> dependent_parameter_values(m, p, 4.0)
(k_EP_to_ES = 0.375,)
```
"""
function dependent_parameter_values(m::Mechanisms.Mechanism, params::NamedTuple, Keq::Real)
    check_valid_rate_equation(m)
    s = _memo_solve(m)
    absent = Symbol[k for k in s.fitted if !haskey(params, k)]
    isempty(absent) || error("dependent_parameter_values: params must hold every name " *
                             "of parameters(m); it lacks " * join(absent, ", "))
    # A negative integer exponent inverts the positive power; a fractional one is taken
    # in the value's own precision, so a BigFloat value keeps its full precision.
    pow(x, e) = !isone(denominator(e)) ?
                    x^(oftype(float(x), numerator(e)) / denominator(e)) :
                numerator(e) >= 0 ? x^numerator(e) : inv(x^-numerator(e))
    value(d) = prod((pow(params[s.names[f]], e) for (f, e) in d.factors);
                    init = pow(Keq, d.keq))
    deps = _dependents(s)
    NamedTuple{Tuple(s.names[d.col] for d in deps)}(Tuple(value(d) for d in deps))
end

# ─── The joint system ───────────────────────────────────────────────────

"""
The kinetic class of column `c` of `m`: its edge's group's constant in the column's slot
and conformation, except that the copy of an edge in the conformation that excludes it
(`excluded`, the steps each conformation excludes, `_excluded_steps`: an `:OnlyA` edge's
inactive copy, an `:OnlyI` edge's active one) is a class of its own.
"""
_class(m::Mechanisms.Mechanism, excluded::NTuple{2, BitVector}, c::_StepParameter) =
    excluded[c.inactive ? 2 : 1][c.step] ?
    (c.step, c.slot, c.inactive, true) :
    (m.catalytic.steps[c.step].kinetic_group, c.slot, c.inactive, false)

"""
    _fundamental_cycles(m, flat) -> (Z, multiples)

The fundamental cycles of a spanning forest of the graph of `m`'s edges `flat`, grown over
the edges in that order: row `i` of `Z` counts each edge of cycle `i` (column `j` is edge
`flat[j]`) with its sign in the edge's stored direction, and `multiples[i]` is the
multiple of the net reaction the cycle runs (0 for a Wegscheider cycle). Errors on a cycle
whose net change is no multiple of the net reaction, a modifier's included, which no
mechanism the guard admits holds (`_unbalanced_cycle`).

The net change reads what each edge takes up (-1) and gives off (+1) only. A form's name
holds no metabolites, so an isomerization, which takes up and gives off nothing,
contributes zero, and the binding and release steps count every metabolite once.
"""
function _fundamental_cycles(m::Mechanisms.CatalyticScheme, flat::Vector{Int})
    forms = eachindex(m.enzyme_forms)
    path = [zeros(Int, length(flat)) for _ in forms]   # signed tree edges from the root
    seen, tree = falses(length(forms)), falses(length(flat))
    for root in forms
        seen[root] && continue
        seen[root] = grown = true
        while grown
            grown = false
            for (j, k) in enumerate(flat)
                a, b = m.steps[k].from, m.steps[k].to
                seen[a] == seen[b] && continue
                from, to, sgn = seen[a] ? (a, b, 1) : (b, a, -1)
                path[to] .= path[from]
                path[to][j] += sgn
                seen[to] = tree[j] = grown = true
            end
        end
    end
    # Each metabolite's change in the net reaction, by its index; 0 for a modifier.
    net = zeros(Int, length(m.metabolites))
    net[m.substrates] .= -1
    net[m.products] .= 1
    x1 = first(m.substrates)
    rows, multiples = Vector{Int}[], Int[]
    for (j, k) in enumerate(flat)
        tree[j] && continue
        a, b = m.steps[k].from, m.steps[k].to
        z = path[a] - path[b]
        z[j] += 1
        change = zeros(Int, length(m.metabolites))
        for (i, kk) in enumerate(flat)
            z[i] == 0 && continue
            foreach(x -> change[x] += z[i], m.steps[kk].given_off)
            foreach(x -> change[x] -= z[i], m.steps[kk].taken_up)
        end
        c = change[x1] * net[x1]
        all(x -> change[x] == c * net[x], eachindex(net)) || error(
            "Cycle $(length(rows) + 1) produces metabolite change not proportional to " *
            "net reaction")
        push!(rows, z)
        push!(multiples, c)
    end
    Z = isempty(rows) ? zeros(Int, 0, length(flat)) : permutedims(reduce(hcat, rows))
    Z, multiples
end

"""
    _joint_system(m, cols = _constant_columns(m)) -> (; columns, A, rhs)

The thermodynamic system of `m`, `A · log(constants) = rhs · log(Keq)`, over one column
per bare constant in elimination order (`_constant_columns`).

Rows:
- every fundamental cycle (`_fundamental_cycles`) once per conformation, the inactive rows
  over the whole graph (islands included), an `:EqualAI` edge reading its active column.
  A cycle row sums, over its edges, the log of each edge's forward equilibrium ratio:
  `kf/kr` at steady state, `1/K` for a constant named after the reverse (such as a
  binding's dissociation constant), `K` otherwise. Its right side is the cycle's multiple
  of the net reaction;
- one row per further member of each kinetic class (`_class`), setting it equal to the
  class's first member, slot by slot.
"""
function _joint_system(m::Mechanisms.Mechanism, cols = _constant_columns(m))
    cm = m.catalytic
    (; columns, flat) = cols
    allosteric = _is_allosteric(m)
    index = Dict(c => i for (i, c) in enumerate(columns))
    # An inactive constant reads its own column when its group has one, else the active.
    column(k, slot, inactive) =
        get(() -> index[_StepParameter(k, slot, false)], index,
            _StepParameter(k, slot, inactive))
    Z, multiples = _fundamental_cycles(cm, flat)
    rows, rhs = Vector{Int}[], Int[]
    for inactive in (allosteric ? (false, true) : (false,)), i in eachindex(multiples)
        row = zeros(Int, length(columns))
        for (j, k) in enumerate(flat)
            z = Z[i, j]
            z == 0 && continue
            re = cm.steps[k].rapid_equilibrium
            reverse_named = re && _constant_named_in_reverse(cm.steps[k])
            row[column(k, 1, inactive)] += reverse_named ? -z : z
            re || (row[column(k, 2, inactive)] -= z)
        end
        push!(rows, row)
        push!(rhs, multiples[i])
    end
    first_of = Dict{Tuple{Int, Int, Bool, Bool}, Int}()
    excluded = (_excluded_steps(m, 1), _excluded_steps(m, 2))
    for (i, c) in enumerate(columns)
        f = get!(first_of, _class(m, excluded, c), i)
        f == i && continue
        row = zeros(Int, length(columns))
        row[i], row[f] = 1, -1
        push!(rows, row)
        push!(rhs, 0)
    end
    A = isempty(rows) ? zeros(Int, 0, length(columns)) : permutedims(reduce(hcat, rows))
    (; columns, A, rhs)
end

# ─── Exact reduction ────────────────────────────────────────────────────

"""
Bring `M` to reduced row echelon form in place over its first `n` columns and return the
pivot columns, in row order. `Rational{Int}` arithmetic throws `OverflowError` rather than
wrap.
"""
function _rref!(M::Matrix{<:Rational}, n::Int)
    pivots = Int[]
    for col in 1:n
        r = length(pivots) + 1
        r > size(M, 1) && break
        p = findfirst(i -> !iszero(M[i, col]), r:size(M, 1))
        p === nothing && continue
        p += r - 1
        for k in axes(M, 2)
            M[r, k], M[p, k] = M[p, k], M[r, k]
        end
        f = inv(M[r, col])
        for k in axes(M, 2)
            M[r, k] *= f
        end
        for i in axes(M, 1)
            (i == r || iszero(M[i, col])) && continue
            g = M[i, col]
            for k in axes(M, 2)
                M[i, k] -= g * M[r, k]
            end
        end
        push!(pivots, col)
    end
    pivots
end

"""
    _reduce_joint_system(A, rhs) -> (M, pivots)

The reduced row echelon form `M` of `[A rhs]` over `A`'s columns and its pivot columns,
in checked `Rational{Int}` arithmetic, redone in `Rational{BigInt}` on overflow. A pivot
column is a dependent constant: row `r` of `M` reads `x_p + Σ_f M[r, f]·x_f =
M[r, end]·log Keq` over the columns `f` that are not pivots. Errors when the rows combine
to `0 = log Keq`.
"""
function _reduce_joint_system(A::Matrix{Int}, rhs::Vector{Int})
    n = size(A, 2)
    M = Matrix{Rational{Int}}([A rhs])
    pivots = try
        _rref!(M, n)
    catch err
        err isa OverflowError || rethrow()
        M = Matrix{Rational{BigInt}}([A rhs])
        _rref!(M, n)
    end
    all(iszero, view(M, (length(pivots) + 1):size(M, 1), n + 1)) || error(
        "Thermodynamically contradictory mechanism: a combination of its " *
        "constraint rows reduces to 0 = log(Keq)")
    M, pivots
end

# ─── The solve ──────────────────────────────────────────────────────────

"""
A dependent constant of the joint solve: pivot column `col` equals `Keq^keq · Π x_f^e`
over the free columns `f => e` of `factors`, in column order. `label` says why it is
dependent: `:kinetic_group` for a column that is not its class's kept member, which it
equals; otherwise `:haldane` when the relation reads `Keq` and `:wegscheider` when it does
not.
"""
struct _Dependent
    col::Int
    keq::Rational{Int}
    factors::Vector{Pair{Int, Rational{Int}}}
    label::Symbol
end

"""
The joint solve of a mechanism (`_solve`). `columns` holds a column per bare constant, in
elimination order, and `names` each column's name, that of its class's kept member
(`_constant_columns`). `column` maps each bare constant to its column; an `:EqualAI`
edge's inactive constant reads its active column. `independents` lists the law's free
columns, sorted by name (`isless` on the strings, so `K_…` precedes `k_…`), and
`dependents` every pivot column of the law, in column order. `info` holds the fitted
parameters in `parameters` order: the independents' columns, then an allosteric
mechanism's regulator parameters (`_regulator_parameters`: those of the active
conformation, then the inactive conformation's own, site by site) and
`_AllostericParameter()`. `fitted` holds their names.
"""
struct _Solve
    columns::Vector{_StepParameter}
    names::Vector{Symbol}
    column::Dict{_StepParameter, Int}
    independents::Vector{Int}
    dependents::Vector{_Dependent}
    info::Vector{_Parameter}
    fitted::Vector{Symbol}
end

"""
    _solve(m) -> _Solve

The joint solve of `m` (`_joint_system`, `_reduce_joint_system`), its exact rationals
narrowed to `Rational{Int}`; a value outside `Int` throws. Every class of constants of the
law is independent or dependent under the name of its kept member (`_constant_columns`).
Columns outside the law (an `:OnlyA` edge's inactive copy, an `:OnlyI` edge's active copy,
a copy on a form its conformation prunes) are eliminated first, so no dependent of the
law reads them.
"""
function _solve(m::Mechanisms.Mechanism)
    cols = _constant_columns(m)
    sys = _joint_system(m, cols)
    M, pivots = _reduce_joint_system(sys.A, sys.rhs)
    n = length(sys.columns)
    names = cols.names
    free = setdiff(1:n, pivots)
    dependents = _Dependent[]
    for (r, p) in enumerate(pivots)
        cols.in_law[p] || continue
        keq, factors = _pivot_row_exponents(M, r, n + 1, free)
        label = !cols.kept_in_law[p] ? :kinetic_group :
                iszero(keq) ? :wegscheider : :haldane
        push!(dependents, _Dependent(p, keq, factors, label))
    end
    independents = sort!(filter(i -> cols.kept_in_law[i], free); by = i -> string(names[i]))
    column = Dict(c => i for (i, c) in enumerate(sys.columns))
    allosteric = _is_allosteric(m)
    if allosteric
        for (i, c) in enumerate(sys.columns)
            key = _StepParameter(c.step, c.slot, true)
            c.inactive || haskey(column, key) || (column[key] = i)
        end
    end
    info = _Parameter[sys.columns[independents]; _regulator_parameters(m)]
    fitted = [names[independents]; _regulator_names(m)]
    if allosteric
        push!(info, _AllostericParameter())
        push!(fitted, _name(_AllostericParameter(), m))
    end
    _Solve(sys.columns, names, column, independents, dependents, info, fitted)
end

"""The exponents row `r` of the reduced system `M` gives its pivot column (`_Dependent`):
`Keq`'s, from the column `keq_column`, and those of the columns `free` it reads, narrowed
to `Rational{Int}`."""
_pivot_row_exponents(M::Matrix{<:Rational}, r::Int, keq_column::Int, free::Vector{Int}) =
    (Rational{Int}(M[r, keq_column]),
     Pair{Int, Rational{Int}}[f => Rational{Int}(-M[r, f])
                              for f in free if !iszero(M[r, f])])

"""The regulator parameters of `m` in `parameters` order: those of the active
conformation, then the inactive conformation's own, site by site; an `:EqualAI` ligand's
one parameter is listed with the active ones, and a ligand has none in a conformation it
is absent from. A plain mechanism has none."""
function _regulator_parameters(am::Mechanisms.Mechanism)
    params(inactive) = _RegulatorParameter[
        _RegulatorParameter(lig, inactive) for site in am.allosteric.regulatory_sites
        for (lig, tag) in zip(site.ligands, site.allosteric_states)
        if !(inactive && tag === :EqualAI) && tag !== (inactive ? :OnlyA : :OnlyI)]
    [params(false); params(true)]
end

"""The parameter names of the regulator constants of `am` (`_regulator_parameters`), in
order."""
_regulator_names(am::Mechanisms.Mechanism) =
    Symbol[_name(p, am) for p in _regulator_parameters(am)]

"""Number of independent (fitted) constants of a mechanism, from the joint solve,
without compiling the mechanism."""
fitted_parameter_count(m::Mechanisms.Mechanism) = length(_solve(m).fitted)
