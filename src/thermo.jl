# ABOUTME: The joint thermodynamic solve: one exact row reduction over the constants of
# ABOUTME: every step in both conformations, with a fixed priority choosing the dependents.

"""
A bare constant of the joint solve: slot `slot` of flat step `step` (`_flat_steps`
order), in the inactive conformation when `inactive`. Slot 1 is a rapid-equilibrium
step's equilibrium constant or a steady-state step's rate constant in its stored
direction; slot 2 is a steady-state step's rate constant in the reverse direction.
"""
struct _Column
    step::Int
    slot::Int
    inactive::Bool
end

# ─── Content keys and orientation ──────────────────────────────────────────

"""A metabolite list as the content keys read it: each metabolite's name and whether it
binds as a modifier copy (a competitive inhibitor), in stored order."""
_metabolites_key(mets) = Tuple((String(name(x)), x isa CompetitiveInhibitor) for x in mets)

"""
An enzyme form's content key, built from base names: its conformation label, its bound
metabolites (`_metabolites_key`) and the names its residual adds and subtracts. Rendered
names never enter it.
"""
_form_base_key(sp::Species) =
    (String(conformation(sp)), _metabolites_key(bound(sp)),
     Tuple(String(name(x)) for x in added(residual(sp))),
     Tuple(String(name(x)) for x in subtracted(residual(sp))))

"""Whether `sp` is free E: conformation label `E`, nothing bound, no residual."""
_is_free_E(sp::Species) =
    conformation(sp) === :E && isempty(bound(sp)) && isempty(residual(sp))

"""
The canonical order of enzyme forms that orients the steps no chemistry orients (C9):
free E first, then by `_form_base_key`.
"""
_form_index_key(sp::Species) = (!_is_free_E(sp), _form_base_key(sp))

"""A step in a given direction as the canonical form order reads it: its two forms
(`_form_index_key`), then the metabolites it takes up and gives off."""
_step_index_key(s::Step) =
    (_form_index_key(from_species(s)), _form_index_key(to_species(s)),
     _metabolites_key(consumed(s)), _metabolites_key(released(s)))

"""
Whether each flat step of `groups`, as stored, runs in the joint solve's forward
direction. Chemistry (orientation tier 1) and graph context (tier 2) orient a step as
the `Mechanism` constructor does; a step neither orients (tier 3) runs from the form
earlier in the canonical form order (`_form_index_key`) to the later one, and the tied
steps of a kinetic group then turn with a lead ranked by `_step_index_key`
(`_orient_tied_steps`). The constructor breaks those last two ties by rendered name, so
for such a step the stored direction can differ from the solve's.
"""
function _solve_forward(groups::Vector{Vector{Step}})
    sides = _entry_sides(groups)
    fwd = Bool[]
    for group in groups
        decided = [_canonical_step_direction(s, sides) for s in group]
        tiers = [last(d) for d in decided]
        by_index = Step[t == 3 ? (_form_index_key(from_species(s)) <=
                                  _form_index_key(to_species(s)) ? s : _reversed(s)) : d
                        for (s, (d, t)) in zip(group, decided)]
        aligned = _orient_tied_steps(by_index, tiers; lead_key = _step_index_key)
        append!(fwd, [a == s for (a, s) in zip(aligned, group)])
    end
    fwd
end

"""
The priority key of constant `slot` of step `s`, whose stored direction is the solve's
forward one when `fwd`, in conformation block `block` (2: an inactive copy outside the
law, 1: another inactive copy, 0: active). The solve eliminates the larger key first, so
a column with a larger key becomes dependent sooner. In turn: the block; the step's kind
(3 an isomerization, 2 other chemistry, 1 a plain binding); whether both of its forms hold
a metabolite; steady state before rapid equilibrium; a steady-state constant that runs
against the reaction (the reverse of the solve's forward direction, except a product
binding's uptake); then the step's content in the solve's direction, its forms'
`_form_base_key`s and its metabolites, the smaller key kept.
"""
function _priority_key(s::Step, fwd::Bool, slot::Int, block::Int)
    kind = is_iso(s) ? 3 : _is_chemistry(s) ? 2 : 1
    occupied = !isempty(bound(from_species(s))) && !isempty(bound(to_species(s)))
    ss = !is_equilibrium(s)
    against = ss && ((slot == 2) == fwd) != (bound_metabolite(s) isa Product)
    a, b, c, r = fwd ? (from_species(s), to_species(s), consumed(s), released(s)) :
                       (to_species(s), from_species(s), released(s), consumed(s))
    (block, kind, occupied, ss, against,
     (_form_base_key(a), _form_base_key(b), _metabolites_key(c), _metabolites_key(r)))
end

# ─── The joint system ──────────────────────────────────────────────────────

"""Each kinetic group's allosteric tag, `:None` for a plain mechanism."""
_solve_tags(m::Mechanism) = fill(:None, length(steps(m)))
_solve_tags(am::AllostericMechanism) = cat_allo_states(am)

"""Whether each flat step lies in the inactive conformation's law (`_inactive_groups`);
false for every step of a plain mechanism."""
_inactive_law(m::Mechanism) = falses(length(_flat_steps(m)))
function _inactive_law(am::AllostericMechanism)
    kept = _inactive_groups(am)
    BitVector([s in kept[g] for (s, g) in _flat_steps(steps(am))])
end

"""Whether the steps of a group with allosteric tag `tag` take their own constants in the
inactive conformation: those of a `:NonequalAI` or `:OnlyA` group do; an `:EqualAI` group's
read the active ones."""
_has_inactive_columns(tag::Symbol) = tag in (:NonequalAI, :OnlyA)

"""
The kinetic class of column `c`, whose step lies in group `g` with tag `tag`: the group's
constant in the column's slot and conformation, except that the inactive copy of each
`:OnlyA` step is a class of its own.
"""
_class(c::_Column, g::Int, tag::Symbol) =
    c.inactive && tag === :OnlyA ? (c.step, c.slot, true, true) :
                                   (g, c.slot, c.inactive, false)

"""
    _fundamental_cycles(rxn, groups) -> (Z, multiples)

The fundamental cycles of a spanning forest of the graph whose nodes are the enzyme forms
of `groups` and whose edges are its flat steps: row `i` of `Z` counts each flat step of
cycle `i` with its sign in the step's stored direction, and `multiples[i]` is the multiple
of the net reaction the cycle runs (0 for a Wegscheider cycle). A metabolite counts by
name, so a reactant's modifier copy counts as the reactant. Errors on a cycle whose net
change is no multiple of the net reaction. Integer arithmetic throughout.
"""
function _fundamental_cycles(rxn::EnzymeReaction, groups::Vector{Vector{Step}})
    flat = _flat_steps(groups)
    forms = _forms(groups)
    idx = Dict(f => i for (i, f) in enumerate(forms))
    ends = [(idx[from_species(s)], idx[to_species(s)]) for (s, _) in flat]
    path = [zeros(Int, length(flat)) for _ in forms]   # signed tree steps from the root
    seen, tree = falses(length(forms)), falses(length(flat))
    for root in eachindex(forms)
        seen[root] && continue
        seen[root] = grown = true
        while grown
            grown = false
            for (k, (a, b)) in enumerate(ends)
                seen[a] == seen[b] && continue
                from, to, sgn = seen[a] ? (a, b, 1) : (b, a, -1)
                path[to] .= path[from]
                path[to][k] += sgn
                seen[to] = tree[k] = grown = true
            end
        end
    end
    net = Dict{Symbol, Int}()
    foreach(x -> net[name(x)] = -1, substrates(rxn))
    foreach(x -> net[name(x)] = 1, products(rxn))
    x1 = name(first(substrates(rxn)))
    rows, multiples = Vector{Int}[], Int[]
    for (k, (a, b)) in enumerate(ends)
        tree[k] && continue
        z = path[a] - path[b]
        z[k] += 1
        change = Dict{Symbol, Int}(x => 0 for x in keys(net))
        for (j, (s, _)) in enumerate(flat)
            z[j] == 0 && continue
            for x in released(s); haskey(net, name(x)) && (change[name(x)] += z[j]); end
            for x in consumed(s); haskey(net, name(x)) && (change[name(x)] -= z[j]); end
        end
        c = change[x1] * net[x1]
        all(x -> change[x] == c * net[x], keys(net)) || error(
            "Cycle $(length(rows) + 1) produces metabolite change not proportional to " *
            "net reaction")
        push!(rows, z)
        push!(multiples, c)
    end
    Z = isempty(rows) ? zeros(Int, 0, length(flat)) : permutedims(reduce(hcat, rows))
    Z, multiples
end

"""
    _joint_system(m) -> (; columns, A, rhs, law, tags, flat)

The thermodynamic system of `m`, `A · log(constants) = rhs · log(Keq)`, over one column
per bare constant (`_Column`): one per slot of every flat step, and for an allosteric
mechanism one more per slot of every step of a `:NonequalAI` or `:OnlyA` group, for its
inactive copy. `columns` lists them in elimination order: by `_priority_key`, larger
first, the steps oriented by `_solve_forward`.

Rows: every fundamental cycle (`_fundamental_cycles`) once per conformation, the inactive
rows over the whole graph (islands included) with an `:EqualAI` step reading its active
column; a cycle row sums, over its steps, the log of each step's forward equilibrium ratio
(`kf/kr` at steady state, `1/K` for a binding's dissociation constant, `K` otherwise), and
its right side is the cycle's multiple of the net reaction. Then one row per further
member of each kinetic class (`_class`), setting it equal to the class's first member,
slot by slot. `law[j]` tells whether flat step `j` lies in the inactive conformation's
law, `tags` holds each group's tag and `flat` the flat steps.
"""
function _joint_system(m::_AnyMech)
    groups = steps(m)
    flat = _flat_steps(groups)
    tags = _solve_tags(m)
    law = _inactive_law(m)
    fwd = _solve_forward(groups)
    allosteric = m isa AllostericMechanism
    columns = _Column[]
    for (j, (s, g)) in enumerate(flat), inactive in (false, true)
        inactive && !(allosteric && _has_inactive_columns(tags[g])) && continue
        for slot in 1:(is_equilibrium(s) ? 1 : 2)
            push!(columns, _Column(j, slot, inactive))
        end
    end
    block(c) = !c.inactive ? 0 : law[c.step] ? 1 : 2
    keys = [_priority_key(flat[c.step][1], fwd[c.step], c.slot, block(c)) for c in columns]
    columns = columns[sortperm(keys; rev = true)]
    index = Dict(c => i for (i, c) in enumerate(columns))
    column(j, slot, inactive) =
        index[_Column(j, slot, inactive && _has_inactive_columns(tags[flat[j][2]]))]
    Z, multiples = _fundamental_cycles(reaction(m), groups)
    rows, rhs = Vector{Int}[], Int[]
    for inactive in (allosteric ? (false, true) : (false,)), i in eachindex(multiples)
        row = zeros(Int, length(columns))
        for (j, (s, _)) in enumerate(flat)
            z = Z[i, j]
            z == 0 && continue
            row[column(j, 1, inactive)] += is_equilibrium(s) && is_binding(s) ? -z : z
            is_equilibrium(s) || (row[column(j, 2, inactive)] -= z)
        end
        push!(rows, row)
        push!(rhs, multiples[i])
    end
    first_of = Dict{NTuple{4, Any}, Int}()
    for (i, c) in enumerate(columns)
        g = flat[c.step][2]
        f = get!(first_of, _class(c, g, tags[g]), i)
        f == i && continue
        row = zeros(Int, length(columns))
        row[i], row[f] = 1, -1
        push!(rows, row)
        push!(rhs, 0)
    end
    A = isempty(rows) ? zeros(Int, 0, length(columns)) : permutedims(reduce(hcat, rows))
    (; columns, A, rhs, law, tags, flat)
end

# ─── Exact reduction ───────────────────────────────────────────────────────

"""
Bring `M` to reduced row echelon form in place over its first `n` columns and return the
pivot columns, in row order. `Rational{Int}` arithmetic throws `OverflowError` rather than
wrap, which `_reduce_joint_system` catches.
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

# ─── Class names ───────────────────────────────────────────────────────────

"""
The kept member of each kinetic class of `m`, keyed by `(group, slot, inactive)`: the
member with the smallest `_priority_key`, which comes last in the joint solve's
elimination order, so the solve keeps it whenever it keeps one of the class (E3). A
`:NonequalAI` group has an inactive class, kept by a member in the inactive law when it has
one; an `:EqualAI` group's steps share their active class, and each `:OnlyA` step's
inactive copy is a class of its own, so neither appears.
"""
function _kept_members(m::_AnyMech)
    flat = _flat_steps(steps(m))
    tags = _solve_tags(m)
    law = _inactive_law(m)
    fwd = _solve_forward(steps(m))
    best = Dict{Tuple{Int, Int, Bool}, Tuple{Any, Step}}()
    for (j, (s, g)) in enumerate(flat), inactive in (false, true)
        inactive && tags[g] !== :NonequalAI && continue
        for slot in 1:(is_equilibrium(s) ? 1 : 2)
            key = _priority_key(s, fwd[j], slot, !inactive ? 0 : law[j] ? 1 : 2)
            cur = get(best, (g, slot, inactive), nothing)
            (cur === nothing || key < first(cur)) && (best[(g, slot, inactive)] = (key, s))
        end
    end
    Dict{Tuple{Int, Int, Bool}, Step}(k => last(v) for (k, v) in best)
end
