# ABOUTME: The joint thermodynamic solve over the constants of both conformations, and the
# ABOUTME: cycle-space helpers of class 3's limit test.

# ─── The joint system ──────────────────────────────────────────────────────

"""
The kinetic class of column `c`, whose edge lies in group `g` with tag `tag`: the group's
constant in the column's slot and conformation, except that the copy of an edge in the
conformation that excludes it (an `:OnlyA` edge's inactive copy, an `:OnlyI` edge's
active one) is a class of its own.
"""
_class(c::Mechanisms._Column, g::Int, tag::Symbol) =
    Mechanisms._excluded_by(tag, c.inactive) ? (c.edge, c.slot, c.inactive, true) :
                                    (g, c.slot, c.inactive, false)

"""
    _fundamental_cycles(m, flat) -> (Z, multiples)

The fundamental cycles of a spanning forest of the graph of graph `m`'s edges `flat`,
grown over the edges in that order: row `i` of `Z` counts each edge of cycle `i` (column
`j` is edge `flat[j]`) with its sign in the edge's stored direction, and `multiples[i]`
is the multiple of the net reaction the cycle runs (0 for a Wegscheider cycle). A
metabolite counts by name, so a reactant's modifier copy counts as the reactant. Errors
on a cycle whose net change is no multiple of the net reaction. Integer arithmetic
throughout.

The net change reads what each edge takes up (-1) and gives off (+1) only. An
isomerization carries no free-pool metabolite — its bound content is encoded in the
enzyme forms it joins — so it contributes zero. Do NOT add a from-bound/to-bound
difference for it: that double-counts metabolites already accounted for by the binding
and release steps and inflates the cycle's net change (e.g. 1/Keq -> 1/Keq^2).
"""
function _fundamental_cycles(m::Mechanisms.Mechanism, flat::Vector{Int})
    forms = Mechanisms.species(m)
    path = [zeros(Int, length(flat)) for _ in forms]   # signed tree edges from the root
    seen, tree = falses(length(forms)), falses(length(flat))
    for root in forms
        seen[root] && continue
        seen[root] = grown = true
        while grown
            grown = false
            for (j, k) in enumerate(flat)
                a, b = Mechanisms.ends(m, k)
                seen[a] == seen[b] && continue
                from, to, sgn = seen[a] ? (a, b, 1) : (b, a, -1)
                path[to] .= path[from]
                path[to][j] += sgn
                seen[to] = tree[j] = grown = true
            end
        end
    end
    net = Dict{Symbol, Int}()
    foreach(x -> net[x] = -1, Mechanisms.substrates(m))
    foreach(x -> net[x] = 1, Mechanisms.products(m))
    x1 = first(Mechanisms.substrates(m))
    rows, multiples = Vector{Int}[], Int[]
    for (j, k) in enumerate(flat)
        tree[j] && continue
        a, b = Mechanisms.ends(m, k)
        z = path[a] - path[b]
        z[j] += 1
        change = Dict{Symbol, Int}(x => 0 for x in keys(net))
        for (i, kk) in enumerate(flat)
            z[i] == 0 && continue
            up, off = Mechanisms.taken_up(m, kk), Mechanisms.given_off(m, kk)
            for x in keys(net)
                change[x] += z[i] * (Mechanisms.copies(off, x) - Mechanisms.copies(up, x))
            end
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
    _joint_system(m, cols = _ordered_columns(m)) -> (; columns, A, rhs, law, tags, flat)

The thermodynamic system of `m`, `A · log(constants) = rhs · log(Keq)`, over one column
per bare constant (`_Column`) in elimination order (`_ordered_columns`).

Rows: every fundamental cycle (`_fundamental_cycles`) once per conformation, the inactive
rows over the whole graph (islands included) with an `:EqualAI` edge reading its active
column; a cycle row sums, over its edges, the log of each edge's forward equilibrium ratio
(`kf/kr` at steady state, `1/K` for a binding's dissociation constant, `K` otherwise), and
its right side is the cycle's multiple of the net reaction. Then one row per further
member of each kinetic class (`_class`), setting it equal to the class's first member,
slot by slot. `law[k]` tells whether edge `k` lies in the inactive conformation's law,
`tags` holds each group's tag and `flat` the edges group by group.
"""
function _joint_system(m::Mechanisms._AnyMech, cols = Mechanisms._ordered_columns(m))
    cm = Mechanisms.catalytic(m)
    (; columns, flat, tags, law) = cols
    allosteric = m isa Mechanisms.AllostericMechanism
    index = Dict(c => i for (i, c) in enumerate(columns))
    column(k, slot, inactive) =
        index[Mechanisms._Column(k, slot, inactive &&
            Mechanisms._has_inactive_columns(tags[Mechanisms.group(cm, k)]))]
    Z, multiples = _fundamental_cycles(cm, flat)
    rows, rhs = Vector{Int}[], Int[]
    for inactive in (allosteric ? (false, true) : (false,)), i in eachindex(multiples)
        row = zeros(Int, length(columns))
        for (j, k) in enumerate(flat)
            z = Z[i, j]
            z == 0 && continue
            re = Mechanisms.is_rapid_equilibrium(cm, k)
            row[column(k, 1, inactive)] += re && Mechanisms.is_binding(cm, k) ? -z : z
            re || (row[column(k, 2, inactive)] -= z)
        end
        push!(rows, row)
        push!(rhs, multiples[i])
    end
    first_of = Dict{Tuple{Int, Int, Bool, Bool}, Int}()
    for (i, c) in enumerate(columns)
        g = Mechanisms.group(cm, c.edge)
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

"""The parameter of column `c` of `m`: slot `c.slot` of its edge's `_edge_constants`, in
the state its name takes: `:None` in a plain mechanism; `:I` for an inactive copy;
`:EqualAI` or `:A` for an active constant of an allosteric mechanism."""
function _column_parameter(m::Mechanisms._AnyMech, c::Mechanisms._Column)
    cm = Mechanisms.catalytic(m)
    state = m isa Mechanisms.Mechanism ? :None : c.inactive ? :I :
            Mechanisms.tag(m, Mechanisms.group(cm, c.edge)) === :EqualAI ? :EqualAI : :A
    Mechanisms._edge_constants(cm, c.edge, state)[c.slot]
end

# ─── The solve ─────────────────────────────────────────────────────────────

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
The joint solve of a mechanism (`_solve`). `columns` holds a column (`_Column`) per bare
constant, in elimination order, and `names` each column's name, that of its class's kept
member (`_kept_members`). `column` maps each bare constant to its column; an `:EqualAI`
edge's inactive constant reads its active column. `independents` lists the law's free
columns, sorted by name (`isless` on the strings, so `K_…` precedes `k_…`), and
`dependents` every pivot column of the law, in column order. `fitted` holds the fitted
constants in `fitted_params` order: the independents' names, then an allosteric
mechanism's regulator constants (`_kreg_params`: active, then inactive, site by site)
and `L`.
"""
struct _Solve
    columns::Vector{Mechanisms._Column}
    names::Vector{Symbol}
    column::Dict{Mechanisms._Column, Int}
    independents::Vector{Int}
    dependents::Vector{_Dependent}
    fitted::Vector{Symbol}
end

"""
    _solve(m) -> _Solve

The joint solve of `m` (`_joint_system`, `_reduce_joint_system`), its exact rationals
narrowed to `Rational{Int}`; a value outside `Int` throws. Every class of constants of the
law is independent or dependent under the name of its kept member (`_kept_members`).
Columns outside the law (an `:OnlyA` edge's inactive copy, an `:OnlyI` edge's active copy,
a copy on a form its conformation prunes) are eliminated first, so no dependent of the
law reads them. The `Type` method decodes a compiled mechanism and checks that its rate
law is defined (`_derivation_input`).

Errors when two constants of the law, regulator constants included, render one name,
which the law would read as one constant: with a conformation labelled `A`, the uptake
constant of the `:EqualAI` step `A(Eb) <--> A + Eb` renders `k_A_Eb_to_AEb`, the name of
the active constant of a `:NonequalAI` step `E(b) <--> A(Eb)`. The macros and the
constructors reject the labels `A` and `I` (`_check_conformation_label`). A metabolite
name can still imitate a residual's: with a modifier named `minusQ`, the side
`E(; residual = A - P) + minusQ` and the form `E(; residual = A - P - Q)` both render
`E_res_plusA_minusP_minusQ`, and the macros accept a step between them, so this check can
fire for a mechanism they build.
"""
function _solve(m::Mechanisms._AnyMech)
    cols = Mechanisms._ordered_columns(m)
    sys = _joint_system(m, cols)
    M, pivots = _reduce_joint_system(sys.A, sys.rhs)
    n = length(sys.columns)
    kept = Mechanisms._kept_members(m, cols)
    names = [Mechanisms.name(_column_parameter(m, c), m, kept) for c in sys.columns]
    cm = Mechanisms.catalytic(m)
    group = [Mechanisms.group(cm, k) for k in Mechanisms.steps(cm)]
    class(c) = _class(c, group[c.edge], sys.tags[group[c.edge]])
    last_of = Dict(class(c) => i for (i, c) in enumerate(sys.columns))
    in_law(c) = c.inactive ? sys.law[c.edge] :
                Mechanisms._in_active_law(c.edge, group, sys.tags, sys.law, cols.alaw)
    kept_in_law(i) = in_law(sys.columns[i]) && last_of[class(sys.columns[i])] == i
    reg = m isa Mechanisms.Mechanism ? Symbol[] :
          Symbol[Mechanisms.name(p, m)
                 for p in [Mechanisms._kreg_params(m, :A); Mechanisms._kreg_params(m, :I)]]
    rendered = [names[filter(kept_in_law, 1:n)]; reg]
    allunique(rendered) || error(
        "Constants of this mechanism's rate equation share a name (" *
        join(unique(x for x in rendered if count(==(x), rendered) > 1), ", ") *
        "); rename a metabolite or a conformation label so that each constant has a " *
        "name of its own")
    free = setdiff(1:n, pivots)
    dependents = _Dependent[]
    for (r, p) in enumerate(pivots)
        in_law(sys.columns[p]) || continue
        keq = Rational{Int}(M[r, n + 1])
        factors = Pair{Int, Rational{Int}}[f => Rational{Int}(-M[r, f])
                                           for f in free if !iszero(M[r, f])]
        label = !kept_in_law(p) ? :kinetic_group : iszero(keq) ? :wegscheider : :haldane
        push!(dependents, _Dependent(p, keq, factors, label))
    end
    independents = sort!(filter(kept_in_law, free); by = i -> string(names[i]))
    column = Dict(c => i for (i, c) in enumerate(sys.columns))
    if m isa Mechanisms.AllostericMechanism
        for (i, c) in enumerate(sys.columns)
            c.inactive || Mechanisms._has_inactive_columns(sys.tags[group[c.edge]]) ||
                (column[Mechanisms._Column(c.edge, c.slot, true)] = i)
        end
    end
    fitted = names[independents]
    m isa Mechanisms.AllostericMechanism && append!(fitted, reg, (:L,))
    _Solve(sys.columns, names, column, independents, dependents, fitted)
end

_solve(@nospecialize(M::Type{<:Mechanisms.AbstractEnzymeMechanism})) =
    _solve(_derivation_input(M))

"""Number of independent (fitted) constants of a concrete mechanism, from the joint solve,
without compiling the mechanism."""
parameter_count(m::Mechanisms._AnyMech) = length(_solve(m).fitted)

# ─── Cycle-space helpers of the limit test ─────────────────────────────────

"""
Reduced row echelon form over `Rational{BigInt}`. Returns the pivot and
free column indices (pivot_cols in row-pivot order, so pivot_cols[i] is the
pivot at reduced-matrix row i) plus the reduced matrix R. Used by
`_rational_nullspace` (nullspace basis).
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
