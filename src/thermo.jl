# ABOUTME: The joint thermodynamic solve over the constants of both conformations, and the
# ABOUTME: cycle-space helpers of the :OnlyA guard and of the split move's counter.

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
_metabolites_key(mets) = Tuple{String, Bool}[(String(name(x)), x isa CompetitiveInhibitor)
                                                for x in mets]

"""
An enzyme form's content key, built from base names: its conformation label, its bound
metabolites (`_metabolites_key`) and the names its residual adds and subtracts. Rendered
names never enter it.
"""
_form_base_key(sp::Species) =
    (String(conformation(sp)), _metabolites_key(bound(sp)),
     String[String(name(x)) for x in added(residual(sp))],
     String[String(name(x)) for x in subtracted(residual(sp))])

"""Whether `sp` is free E: conformation label `E`, nothing bound, no residual."""
_is_free_E(sp::Species) =
    conformation(sp) === :E && isempty(bound(sp)) && isempty(residual(sp))

"""
The canonical order of enzyme forms that orients the steps neither chemistry nor graph
context orients (C9): free E first, then by `_form_base_key`.
"""
_form_index_key(sp::Species) = (!_is_free_E(sp), _form_base_key(sp))

"""A step in a given direction as the canonical form order reads it: its two forms
(`_form_index_key`), then the metabolites it takes up and gives off."""
_step_index_key(s::Step) =
    (_form_index_key(from_species(s)), _form_index_key(to_species(s)),
     _metabolites_key(consumed(s)), _metabolites_key(released(s)))

"""
The priority key of constant `slot` of step `s` in conformation block `block` (2: an
inactive copy outside the law, 1: another inactive copy, 0: active). The solve eliminates
the larger key first, so a column with a larger key becomes dependent sooner. In turn: the
block; the step's kind (3 an isomerization, 2 other chemistry, 1 a plain binding); whether
both of its forms hold a metabolite; steady state before rapid equilibrium; a steady-state
constant that runs against the reaction (for a product binding its uptake, for every other
step its reverse rate: the constructor stores each step in its physical-forward
direction); then the step's content as stored, its forms' `_form_base_key`s and its
metabolites, the smaller key kept.
"""
function _priority_key(s::Step, slot::Int, block::Int)
    kind = is_iso(s) ? 3 : _is_chemistry(s) ? 2 : 1
    occupied = !isempty(bound(from_species(s))) && !isempty(bound(to_species(s)))
    ss = !is_equilibrium(s)
    against = ss && (slot == 2) != (bound_metabolite(s) isa Product)
    (block, kind, occupied, ss, against,
     (_form_base_key(from_species(s)), _form_base_key(to_species(s)),
      _metabolites_key(consumed(s)), _metabolites_key(released(s))))
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
    _fundamental_cycles(m, groups) -> (Z, multiples)

The fundamental cycles of a spanning forest of the graph whose nodes are the enzyme forms
of `groups`, kinetic groups over `m`'s substrates and products, and whose edges are its
flat steps: row `i` of `Z` counts each flat step of cycle `i` with its sign in the step's
stored direction, and `multiples[i]` is the multiple of the net reaction the cycle runs
(0 for a Wegscheider cycle). A metabolite counts by name, so a reactant's modifier copy
counts as the reactant. Errors on a cycle whose net change is no multiple of the net
reaction. Integer arithmetic throughout.
"""
function _fundamental_cycles(m::_AnyMech, groups::Vector{Vector{Step}})
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
    foreach(x -> net[x] = -1, substrates(m))
    foreach(x -> net[x] = 1, products(m))
    x1 = first(substrates(m))
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
first.

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
    allosteric = m isa AllostericMechanism
    columns = _Column[]
    for (j, (s, g)) in enumerate(flat), inactive in (false, true)
        inactive && !(allosteric && _has_inactive_columns(tags[g])) && continue
        for slot in 1:(is_equilibrium(s) ? 1 : 2)
            push!(columns, _Column(j, slot, inactive))
        end
    end
    block(c) = !c.inactive ? 0 : law[c.step] ? 1 : 2
    priorities = [_priority_key(flat[c.step][1], c.slot, block(c)) for c in columns]
    columns = columns[sortperm(priorities; rev = true)]
    index = Dict(c => i for (i, c) in enumerate(columns))
    column(j, slot, inactive) =
        index[_Column(j, slot, inactive && _has_inactive_columns(tags[flat[j][2]]))]
    Z, multiples = _fundamental_cycles(m, groups)
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
    best = Dict{Tuple{Int, Int, Bool}, Tuple{Any, Step}}()
    for (j, (s, g)) in enumerate(flat), inactive in (false, true)
        inactive && tags[g] !== :NonequalAI && continue
        for slot in 1:(is_equilibrium(s) ? 1 : 2)
            key = _priority_key(s, slot, !inactive ? 0 : law[j] ? 1 : 2)
            cur = get(best, (g, slot, inactive), nothing)
            (cur === nothing || key < first(cur)) && (best[(g, slot, inactive)] = (key, s))
        end
    end
    Dict{Tuple{Int, Int, Bool}, Step}(k => last(v) for (k, v) in best)
end

"""The parameter of column `c` of `m`, whose flat steps are `flat`: slot `c.slot` of the
step's `_step_constants`, in the state its name takes: `:None` in a plain mechanism; `:I`
for an inactive copy; `:EqualAI` or `:A` for an active constant of an allosteric
mechanism."""
function _column_parameter(m::_AnyMech, flat, c::_Column)
    s, g = flat[c.step]
    state = m isa Mechanism ? :None : c.inactive ? :I :
            cat_allo_state(m, g) === :EqualAI ? :EqualAI : :A
    _step_constants(s, state)[c.slot]
end

# ─── Dependents ────────────────────────────────────────────────────────────

"""
    _dependent_param_exprs(m) -> (dep, indep)

The dependent constants of `m`'s law and its independent ones, from one joint solve
(`_joint_system`, `_reduce_joint_system`). Every class of constants of the law is
independent or dependent, under the name of its kept member (`_kept_members`):

- `indep` holds the classes whose kept column is not a pivot, sorted by name (`isless`
  on the strings, so `K_…` precedes `k_…`); an allosteric mechanism's list continues with
  its regulator constants (`_kreg_params`: active, then inactive, site by site) and `L`.
- `dep` maps each class whose kept column is a pivot to `Keq^q · Π independent^e` with
  rational exponents (`build_power_expr`); the empty product is `1.0`.

Columns outside the law (an `:OnlyA` step's inactive copy, an inactive copy on a pruned
form) are eliminated first, so no dependent of the law reads them; a pivot that is not
its class's kept member equals the kept one and adds no entry. The `Type` method decodes
a compiled mechanism and checks that its rate law is defined (`_derivation_input`).

Errors when two constants of the law, regulator constants included, render one name,
which the law would read as one constant: with a conformation labelled `A`, the uptake
constant of the `:EqualAI` step `A(Eb) <--> A + Eb` renders `k_A_Eb_to_AEb`, the name of
the active constant of a `:NonequalAI` step `E(b) <--> A(Eb)`.
"""
function _dependent_param_exprs(m::_AnyMech)
    sys = _joint_system(m)
    M, pivots = _reduce_joint_system(sys.A, sys.rhs)
    n = length(sys.columns)
    sides = _class_sides(m)
    names = [name(_column_parameter(m, sys.flat, c), m, sides) for c in sys.columns]
    class(c) = _class(c, sys.flat[c.step][2], sys.tags[sys.flat[c.step][2]])
    last_of = Dict(class(c) => i for (i, c) in enumerate(sys.columns))
    kept_in_law(i) = (c = sys.columns[i];
                      (!c.inactive || sys.law[c.step]) && last_of[class(c)] == i)
    reg = m isa Mechanism ? Symbol[] :
          Symbol[name(p, m) for p in [_kreg_params(m, :A); _kreg_params(m, :I)]]
    rendered = [names[filter(kept_in_law, 1:n)]; reg]
    allunique(rendered) || error(
        "Constants of this mechanism's rate equation share a name (" *
        join(unique(x for x in rendered if count(==(x), rendered) > 1), ", ") *
        "); rename a metabolite or a conformation label so that each constant has a " *
        "name of its own")
    free = setdiff(1:n, pivots)
    dep = Dict{Symbol, Union{Symbol, Expr, Float64}}()
    for (r, p) in enumerate(pivots)
        kept_in_law(p) || continue
        factors = [(names[f], -M[r, f]) for f in free if !iszero(M[r, f])]
        dep[names[p]] = build_power_expr(M[r, n + 1], factors)
    end
    indep = sort!([names[f] for f in free if kept_in_law(f)]; by = string)
    m isa Mechanism && return dep, Tuple(indep)
    dep, Tuple([indep; reg; :L])
end

_dependent_param_exprs(@nospecialize(M::Type{<:AbstractEnzymeMechanism})) =
    _dependent_param_exprs(_derivation_input(M))

"""Number of independent (fitted) constants of a concrete mechanism, from the joint solve,
without compiling the mechanism."""
_independent_param_count(m::_AnyMech) = length(_dependent_param_exprs(m)[2])

# ─── Cycle-space helpers of the :OnlyA guard and the split move's counter ─────

"""
Reduced row echelon form over `Rational{BigInt}`. Returns the pivot and
free column indices (pivot_cols in row-pivot order, so pivot_cols[i] is the
pivot at reduced-matrix row i) plus the reduced matrix R. Used by
`_rational_nullspace` (nullspace basis) and `_partition_independent_count` (rank).
"""
function _rref_partition(A::AbstractMatrix)
    m, n = size(A)
    R = Matrix{Rational{BigInt}}(A)
    pivot_cols = Int[]
    row = 1
    for col in 1:n
        piv = findfirst(r -> R[r, col] != 0, row:m)
        piv === nothing && continue
        piv += row - 1
        R[row, :], R[piv, :] = R[piv, :], R[row, :]
        R[row, :] ./= R[row, col]
        for r in 1:m
            r != row && R[r, col] != 0 && (R[r, :] .-= R[r, col] .* R[row, :])
        end
        push!(pivot_cols, col); row += 1
    end
    free_cols = setdiff(1:n, pivot_cols)
    return pivot_cols, free_cols, R
end

"""
Basis of `{x : M·x = 0}` as the COLUMNS of the returned matrix, exact over the
rationals: each free column of `M`'s reduced row echelon form yields one basis
vector, 1 at its free coordinate.
"""
function _rational_nullspace(M::AbstractMatrix)
    pivots, free, R = _rref_partition(M)
    N = zeros(Rational{BigInt}, size(M, 2), length(free))
    for (k, f) in enumerate(free)
        N[f, k] = 1
        for (r, c) in enumerate(pivots); N[c, k] = -R[r, f]; end
    end
    N
end

"""
`_rational_nullspace(A)` with each basis column scaled to coprime integers whose
first nonzero entry is positive.
"""
function _integer_nullspace(A::Matrix{Int})
    N = _rational_nullspace(A)
    result = zeros(Int, size(N))
    for k in axes(N, 2)
        col = Int.(N[:, k] .* lcm(denominator.(N[:, k])))
        col .÷= gcd(col)
        result[:, k] .= col[findfirst(!iszero, col)] < 0 ? -col : col
    end
    result
end

_thermodynamic_constraints(mech::Mechanism) = _thermodynamic_constraints(mech, steps(mech))

function _thermodynamic_constraints(m::_AnyMech, groups::Vector{Vector{Step}})
    flat = _flat_steps(groups)
    reactant_names = [substrates(m); products(m)]

    # Enzyme incidence matrix (rows = enzyme forms in the canonical form order, `_forms`)
    form = Dict(sp => i for (i, sp) in enumerate(_forms(groups)))
    B = zeros(Int, length(form), length(flat))
    for (j, (s, _)) in enumerate(flat)
        B[form[from_species(s)], j] -= 1
        B[form[to_species(s)], j] += 1
    end

    # Stoichiometry matrix (rows = metabolites, cols = steps), read from each
    # step's consumed (-1) and released (+1) lists.
    #
    # Iso steps carry no free-pool metabolite — their bound content is
    # encoded in the enzyme-form identity — so both lists are empty and they
    # contribute zero. Do NOT add a from_bound/to_bound diff for iso steps:
    # that double-counts metabolites already accounted for by the
    # binding/release steps and inflates the cycle's net change (e.g.
    # 1/Keq -> 1/Keq^2).
    met_idx = Dict(x => i for (i, x) in enumerate(reactant_names))
    stoich_mat = zeros(Int, length(reactant_names), length(flat))
    for (j, (s, _)) in enumerate(flat), (mets, sgn) in ((consumed(s), -1), (released(s), 1))
        for x in mets
            haskey(met_idx, name(x)) && (stoich_mat[met_idx[name(x)], j] += sgn)
        end
    end
    nu_net = [x in substrates(m) ? -1 : 1 for x in reactant_names]

    # Classify each null-space cycle as Haldane (proportional to the
    # net reaction → contributes log(Keq)) or Wegscheider (closed
    # cycle, zero net change). Errors on cycles that touch metabolites
    # but aren't proportional to the net reaction.
    C = _integer_nullspace(B)'
    j0 = findfirst(!iszero, nu_net)
    function classify(i)
        nu = stoich_mat * C[i, :]
        c = nu[j0] // nu_net[j0]
        nu == c .* nu_net && isinteger(c) ||
            error("Cycle $i produces metabolite change not proportional to net reaction")
        Int(c)
    end
    return C, Int[classify(i) for i in axes(C, 1)]
end

"""How a step's constants enter the constraint columns: `:ss` (a forward and a
reverse rate), `:binding_K` (a dissociation constant, whose column carries a sign
flip) or `:iso_K` (an equilibrium constant)."""
_count_kind(s::Step) = is_equilibrium(s) ? (is_binding(s) ? :binding_K : :iso_K) : :ss

"""
Add flat step `j`'s cycle incidence `C[:, j]` into the constraint matrix `A` at the
columns `cols` of its constants, by the step's constant kind (`_count_kind`). A
binding K enters negated, because it is a Kd in the polynomial while the cycle
product uses 1/Kd; an iso K enters as is; an SS step contributes `+kf` at `cols[1]`
and `-kr` at `cols[2]`.
"""
function _add_step_column!(A, C, j, kind::Symbol, cols)
    for i in axes(C, 1)
        c = C[i, j]
        c == 0 && continue
        if kind === :ss
            A[i, cols[1]] += c
            A[i, cols[2]] -= c
        else
            A[i, cols[1]] += kind === :binding_K ? -c : c
        end
    end
end

"""
    _partition_independent_count(parent::Mechanism) -> counter

Return `counter(group_of_step, kind_of_step = kinds of the parent's steps)`, the
independent-parameter count of the mechanism obtained by regrouping `parent`'s
flat steps (in `_flat_steps` order) into the groups labelled by `group_of_step`,
each step counted under the constant kind `kind_of_step` gives it (`_count_kind`).
Regrouping moves no edges and a flag changes none, so the cycle basis of the step
graph (`_thermodynamic_constraints`) is the same for every call and is computed
once here; each call only merges step columns by group and takes the rank. The
split move passes a step's own kind, or `:binding_K`/`:iso_K` for a steady-state
step it reverts to rapid equilibrium. Equals `_independent_param_count` of the
constructed child: the joint solve's independent set is its columns of the law less
its pivots, and a kinetic group's constants form one column here and one class there.
"""
function _partition_independent_count(parent::Mechanism)
    C, _ = _thermodynamic_constraints(parent)
    kinds = [_count_kind(s) for (s, _) in _flat_steps(parent)]
    function counter(group_of_step::AbstractVector{Int},
                     kind_of_step::AbstractVector{Symbol} = kinds)
        length(group_of_step) == length(kinds) ||
            error("group_of_step must label every flat step of the parent")
        length(kind_of_step) == length(kinds) ||
            error("kind_of_step must label every flat step of the parent")
        column = Dict{Tuple{Int, Int}, Int}()
        for (j, g) in enumerate(group_of_step)
            get!(column, (g, 1), length(column) + 1)
            kind_of_step[j] === :ss && get!(column, (g, 2), length(column) + 1)
        end
        A = zeros(Int, size(C, 1), length(column))
        for (j, g) in enumerate(group_of_step)
            kind = kind_of_step[j]
            cols = [column[(g, k)] for k in 1:(kind === :ss ? 2 : 1)]
            _add_step_column!(A, C, j, kind, cols)
        end
        length(column) - length(_rref_partition(A)[1])
    end
    counter
end

"""
True when some `y` makes every row of `N·y` strictly positive — i.e. the open
cone `{y : N·y > 0}` is nonempty. Exact Fourier-Motzkin elimination: to drop
variable `v`, every (positive, negative) row pair is combined with positive
coefficients, which preserves strictness; a row that reduces to all-zero
encodes `0 > 0` and refutes the system. Returns `nothing` if elimination blows
up, so the caller can fall back to the sound per-row test rather than error on
the enumeration hot path.
"""
function _has_strict_positive_combination(N::AbstractMatrix{Rational{BigInt}})
    d = size(N, 2)
    d == 0 && return false
    rows = [collect(N[i, :]) for i in axes(N, 1)]
    for v in d:-1:1
        pos = [r for r in rows if r[v] > 0]
        neg = [r for r in rows if r[v] < 0]
        nxt = [r for r in rows if r[v] == 0]
        length(nxt) + length(pos) * length(neg) > 4000 && return nothing
        for p in pos, q in neg
            push!(nxt, p .* (-q[v]) .+ q .* p[v])
        end
        any(r -> all(iszero, r), nxt) && return false
        rows = nxt
    end
    true
end

"""
The error text for the offending `:OnlyA` binding groups `groups` of `am`
(`_onlya_haldane_violation`), each named by its constant in `am`'s plain catalytic
mechanism: its first step's slot-1 constant without a state tag, which takes the name of
the group's active-class kept member.
"""
function _onlya_violation_message(am::AllostericMechanism, groups::Vector{Int})
    sides = _class_sides(am)
    labels = sort!([string(name(first(_step_constants(first(steps(am)[g]), :None)), am,
                                sides))
                    for g in groups])
    "an :OnlyA binding ($(join(labels, ", "))) leaves a " *
    "thermodynamic (Haldane/Wegscheider) cycle unsatisfiable: the " *
    "inactive conformation cannot close that cycle at finite nonzero " *
    "affinity. Tag the cycle's chemical step :OnlyA, or tag an " *
    "opposing binding :OnlyA so the affinities diverge together."
end

"""
    _onlya_haldane_violation(m, cat_steps, cat_allo_states) → Vector{Int}

The `:OnlyA` binding groups of `cat_steps`, kinetic groups over mechanism `m`'s
metabolites tagged `cat_allo_states`, that leave a catalytic thermodynamic cycle's
Haldane relation unsatisfiable under the `K_I = K_A/ε`, `ε → 0⁺` limit an `:OnlyA`
binding asserts; empty when every cycle stays satisfiable. The boundary check
(`_check_well_posed`) names them in its error (`_onlya_violation_message`); the moves read
only whether any group offends.

A cycle's Haldane carries `∏ε_p/∏ε_s` on the inactive side. The `ε` are
independent, so that monomial can be held finite only when its exponents
carry both signs. All-same-sign drives it to `0` or `∞`; the only thing that
absorbs the imbalance is a free inactive catalytic ratio `k_I_f/k_I_r`. An
`:OnlyA` catalytic tag supplies exactly that — its rate constants vanish from
the rate equation, so their ratio is free to satisfy the cycle's Haldane at any
affinity — which is why the inactive Haldane is present but never binding here.

The check graph drops `:OnlyA` chemistry groups (those holding a chemistry step,
`_is_chemistry`: an isomerization, a fused binding or a Theorell–Chance step), so
a cycle running through one never appears and never reports a violation: that
free `k_I` ratio is the escape. The `:OnlyA` bindings are the `:OnlyA` groups of
plain bindings alone; a fused binding tagged `:OnlyA` is chemistry and leaves the
graph. The check graph is read for its cycles only, so the kept groups enter as
stored, not as a mechanism: dropping a group that mixes a fused binding with
abortive bindings can leave abortive complexes in a rapid-equilibrium segment with
no bottom form, which the boundary check rejects (`_bottomless_re_segment`).
Bindings completing no cycle (competitive inhibitors, dead ends, regulator
sites) never enter a row and take no part. Both catalytic (Haldane) and
binding-only (Wegscheider, `rhs = 0`) cycle rows are inspected, so a one-sided
`:OnlyA` binding on a pure random-order binding square is caught.

The verdict comes from two stages. The per-row sign test runs first: it flags a
cycle whose single row's `:OnlyA` exponents are all one sign, which forces a sum
of same-signed positive terms to vanish and so is a sound rejection. It is cheap
and keeps the common rejections on the fast path, but it reads one row at a time
and a multi-cycle coupled inconsistency passes it. The complete condition
follows: writing `ε_p = δ^{w_p}`, every cycle's monomial is finite nonzero for
some `δ → 0⁺` exactly when the coupled `ε`-exponent system admits a
strictly-positive nullspace vector, `M·w = 0` with `w > 0` (Stiemke
feasibility). The two stages agree for every random-order mechanism up to bi-bi;
they part from ter-substrate up.

An `:OnlyA` group's `ε` exponent in a cycle is the sum of its steps' entries in
that row of the cycle basis (`_thermodynamic_constraints`), for RE and SS
bindings alike: an RE binding's `1/Kd` and an SS binding's `kon/koff` both enter
the cycle's product raised to the step's entry, so a cycle that mixes the two step
kinds reads its signs on one scale and a balanced pair is accepted.
"""
function _onlya_haldane_violation(m::_AnyMech, cat_steps::Vector{Vector{Step}},
                                  cat_allo_states::Vector{Symbol})
    onlya(g) = cat_allo_states[g] === :OnlyA
    keep = [g for g in eachindex(cat_steps)
            if !(onlya(g) && any(_is_chemistry, cat_steps[g]))]
    any(onlya, keep) || return Int[]
    kept = cat_steps[keep]
    C, _ = _thermodynamic_constraints(m, kept)
    flat = _flat_steps(kept)
    groups = [k for (k, g) in enumerate(keep) if onlya(g)]
    # Column k holds `:OnlyA` group k's `ε` exponents: its steps' cycle entries summed.
    M = Rational{BigInt}[sum((C[i, j] for (j, (_, g)) in enumerate(flat) if g == k);
                             init = 0) for i in axes(C, 1), k in groups]
    # Per-row sign test first: sound (an all-one-sign row forces a sum of
    # same-signed positive terms to vanish) and cheap, so it keeps the common
    # rejections on the fast path.
    for row in eachrow(M)
        nz = findall(!iszero, row)
        !isempty(nz) && allequal(sign.(row[nz])) && return keep[groups[nz]]
    end

    # The complete condition: the `ε` exponents must admit a strictly-positive
    # solution of `M·w = 0` (Stiemke feasibility). The per-row test above only
    # inspects one row at a time, so from ter-substrate up a multi-cycle coupled
    # inconsistency passes it. `nothing` from the cone test means elimination
    # blew up; fall back to the per-row verdict (sound, just incomplete) rather
    # than reject a mechanism we could not decide.
    _has_strict_positive_combination(_rational_nullspace(M)) === false || return Int[]
    return keep[groups]
end
