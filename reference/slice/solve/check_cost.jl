# ABOUTME: Solve-stage check (f): the RREF's largest intermediate entry in checked Int64, the
# ABOUTME: system sizes, the cost of a solve and of the API, and RE-only cycle redundancy.
include(joinpath(@__DIR__, "lib.jl"))

# Only the converted mechanisms (no Pop records), so the timings run on a small heap.
records = [(set = String(l[1]), p = (id = c.id,), c = c)
           for l in [split(x, '\t') for x in readlines(joinpath(DATA, "MANIFEST.tsv"))[2:end]]
           for c in deserialize(joinpath(DATA, l[1] * ".jls"))
           if (try S.solve(c.mech); true catch; false end)]

"""`Slice.rref!` in checked `Rational{Int}` that also returns the largest numerator or
denominator of any entry it writes."""
function rref_track!(M::Matrix{Rational{Int}}, n::Int)
    piv, top = Int[], 1
    seen(x) = (top = max(top, abs(numerator(x)), denominator(x)); x)
    for col in 1:n
        r = length(piv) + 1
        r > size(M, 1) && break
        p = findfirst(i -> !iszero(M[i, col]), r:size(M, 1))
        p === nothing && continue
        p += r - 1
        for k in axes(M, 2)
            M[r, k], M[p, k] = M[p, k], M[r, k]
        end
        f = inv(M[r, col])
        foreach(k -> M[r, k] = seen(M[r, k] * f), axes(M, 2))
        for i in axes(M, 1)
            (i == r || iszero(M[i, col])) && continue
            g = M[i, col]
            foreach(k -> M[i, k] = seen(M[i, k] - g * M[r, k]), axes(M, 2))
        end
        push!(piv, col)
    end
    piv, top
end

"""Nanoseconds per call of `f`: the best of three 2 ms loops after one warm call."""
function cost(f)
    f()
    best = Inf
    for _ in 1:3
        n, t0 = 0, time_ns()
        while time_ns() - t0 < 2_000_000
            f()
            n += 1
        end
        best = min(best, (time_ns() - t0) / n)
    end
    best
end

"""Rows of the fundamental cycles of the rapid-equilibrium subgraph, in each state."""
function re_cycle_rows(m, cols)
    fwd = S.orient(m)
    index = Dict(c => i for (i, c) in enumerate(cols))
    seg, parent, order = S.re_forest(m)
    E = length(m.edges)
    path = [zeros(Int, E) for _ in m.forms]
    for i in order
        k = parent[i]
        k == 0 && continue
        e = m.edges[abs(k)]
        path[i] = copy(path[k > 0 ? e.a : e.b])
        path[i][abs(k)] += sign(k)
    end
    tree = Set(abs(k) for k in parent if k != 0)
    rows = Vector{Int}[]
    for (k, e) in enumerate(m.edges), inactive in (S.is_allosteric(m) ? (false, true) : (false,))
        e.re && !(k in tree) || continue
        z = path[e.a] - path[e.b]; z[k] += 1
        row = zeros(Int, length(cols))
        for (j, zj) in enumerate(z)
            zj == 0 && continue
            zf = fwd[j] ? zj : -zj
            o = S.oriented(m, j, fwd)
            c = index[S.Col(j, 1, inactive && m.tags[m.group[j]] !== :EqualAI)]
            row[c] += S.is_binding(o) ? -zf : zf
        end
        push!(rows, row)
    end
    rows
end

big_all, overflow, mism = 1, String[], String[]
rows_max, cols_max, re_rows, re_changed = 0, 0, 0, 0
byset = Dict{String, Vector{Any}}()
for r in records
    m = r.c.mech
    cols, A, b = S.constraint_system(m)
    n = length(cols)
    global rows_max = max(rows_max, size(A, 1))
    global cols_max = max(cols_max, n)
    M = Rational{Int}.([A b])
    mx = try
        piv, top = rref_track!(M, n)
        MB = Rational{BigInt}.([A b])
        S.rref!(MB, n) == piv && MB == M || push!(mism, "$(r.set) $(r.p.id)")
        top
    catch err
        err isa OverflowError || rethrow()
        push!(overflow, "$(r.set) $(r.p.id)")
        -1
    end
    global big_all = max(big_all, mx)
    extra = re_cycle_rows(m, cols)
    global re_rows += length(extra)
    if !isempty(extra)
        base = rank_exact([A b])
        R = permutedims(reduce(hcat, extra))
        rank_exact([A b; R zeros(Int, size(R, 1))]) == base || (global re_changed += 1)
    end
    t_solve = cost(() -> S.solve(m))
    t_sys = cost(() -> S.constraint_system(m))
    t_fit = cost(() -> S.fitted_params(m))
    pars = (; (x => 1.5 for x in S.parameters(m))...)
    t_all = cost(() -> S.all_parameters(m, pars))
    push!(get!(byset, r.set, Any[]), (mx, t_solve, t_sys, t_fit, t_all, size(A)...))
end
@printf("records: %d; system at most %d rows x %d columns\n", length(records), rows_max,
        cols_max)
@printf("RREF in checked Int64: overflow in %d; largest numerator or denominator written %d; \
        result differs from the Rational{BigInt} RREF in %d\n", length(overflow), big_all,
        length(mism))
@printf("RE-only cycle rows: %d added over all records; rank changed in %d\n", re_rows,
        re_changed)
q(v, p) = quantile(v, p)
println("per set: largest entry; solve, constraint_system, fitted_params, all_parameters in µs \
         (median / p90 / max)")
allt = [x for v in values(byset) for x in v]
for (set, v) in sort!(collect(byset); by = first)
    @printf("  %-17s %4d  big %d  solve %6.1f/%6.1f/%6.1f  system %6.1f  fitted %6.1f  \
            all %6.1f  rows<=%d cols<=%d\n", set, length(v), maximum(x[1] for x in v),
            q([x[2] for x in v], 0.5) / 1e3, q([x[2] for x in v], 0.9) / 1e3,
            maximum(x[2] for x in v) / 1e3, q([x[3] for x in v], 0.5) / 1e3,
            q([x[4] for x in v], 0.5) / 1e3, q([x[5] for x in v], 0.5) / 1e3,
            maximum(x[6] for x in v), maximum(x[7] for x in v))
end
@printf("all: solve median %.1f µs, p90 %.1f, max %.1f; fitted_params median %.1f µs; \
        all_parameters median %.1f µs\n", q([x[2] for x in allt], 0.5) / 1e3,
        q([x[2] for x in allt], 0.9) / 1e3, maximum(x[2] for x in allt) / 1e3,
        q([x[4] for x in allt], 0.5) / 1e3, q([x[5] for x in allt], 0.5) / 1e3)
foreach(x -> println("  overflow ", x), overflow)
foreach(x -> println("  mismatch ", x), mism)
