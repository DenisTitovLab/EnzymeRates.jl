# ABOUTME: The joint thermodynamic solve over the constants of both conformations, and the
# ABOUTME: cycle-space helpers of class 3's limit test.

# ─── The joint system ──────────────────────────────────────────────────────

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

The fundamental cycles of a spanning forest of the graph of graph `m`'s edges `flat`,
grown over the edges in that order: row `i` of `Z` counts each edge of cycle `i` (column
`j` is edge `flat[j]`) with its sign in the edge's stored direction, and `multiples[i]`
is the multiple of the net reaction the cycle runs (0 for a Wegscheider cycle). Errors
on a cycle whose net change is no multiple of the net reaction, a modifier's included,
which no mechanism the guard admits holds (`_unbalanced_cycle`). Integer arithmetic
throughout.

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
    _joint_system(m, cols = _constant_columns(m)) -> (; columns, A, rhs, law, tags, flat)

The thermodynamic system of `m`, `A · log(constants) = rhs · log(Keq)`, over one column
per bare constant in elimination order (`_constant_columns`).

Rows: every fundamental cycle (`_fundamental_cycles`) once per conformation, the inactive
rows over the whole graph (islands included) with an `:EqualAI` edge reading its active
column; a cycle row sums, over its edges, the log of each edge's forward equilibrium ratio
(`kf/kr` at steady state, `1/K` for a constant named after the reverse, such as a
binding's dissociation constant, `K` otherwise), and its right side is the cycle's
multiple of the net reaction. Then one row per further
member of each kinetic class (`_class`), setting it equal to the class's first member,
slot by slot. `law[k]` tells whether edge `k` lies in the inactive conformation's law,
`tags` holds each group's tag and `flat` the edges group by group.
"""
function _joint_system(m::Mechanisms.Mechanism, cols = _constant_columns(m))
    cm = m.catalytic
    (; columns, flat, tags, law) = cols
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

"""The names of the regulator parameters of `m` (`_regulator_parameters`), in their
order."""
_regulator_names(am::Mechanisms.Mechanism) =
    Symbol[_name(p, am) for p in _regulator_parameters(am)]

"""Number of independent (fitted) constants of a mechanism, from the joint solve,
without compiling the mechanism."""
fitted_parameter_count(m::Mechanisms.Mechanism) = length(_solve(m).fitted)

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
