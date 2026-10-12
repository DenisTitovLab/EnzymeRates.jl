# ABOUTME: K1 classes per reaction pool: every merge verified through the oracles with the
# ABOUTME: canonical name map, compared with P8's classes; usage: check_classes.jl label pools…
include(joinpath(@__DIR__, "lib.jl"))

const LABEL, POOLS = ARGS[1], [split(p, "+") for p in ARGS[2:end]]
const P8DATA = joinpath(PROTO, "p8", "data")

"""Column values of `m`'s law from parameter values `x` (by name), as `random_values`."""
function column_values(m, sol, x)
    names = S.column_names(m, sol)
    cv = Dict{Int, BigFloat}(i => x[names[i]] for i in sol.indep)
    for d in sol.deps
        cv[d.col] = prod((S.pow(cv[f], e) for (f, e) in d.terms); init = S.pow(x[:Keq], d.keq))
    end
    cv
end

"""The largest relative difference, over `n` random points, between the oracle rates of the
representative `a` and the member `b`, `b`'s constants read off `a`'s through the canonical
positions (`oa[i]` is `ob[i]`) and its concentrations by position in `metabolites`."""
function merge_error(a, b, oa, ob, rng; n = 3)
    sa, sb = S.solve(a.c.mech), S.solve(b.c.mech)
    worst = 0.0
    ma, mb = S.metabolites(a.c.mech), S.metabolites(b.c.mech)
    for _ in 1:n
        x, _, conc = random_point(a, sa, rng)
        y = Dict{Symbol, BigFloat}(n => get(x, n, big(1)) for n in S.parameters(b.c.mech))
        foreach(((u, w),) -> y[w] = x[u], zip(oa, ob))
        cb = Dict(mb[i] => conc[ma[i]] for i in eachindex(ma))     # concentrations by position
        va = oracle_value(a, sa, x, column_values(a.c.mech, sa, x), conc)
        vb = oracle_value(b, sb, y, column_values(b.c.mech, sb, y), cb)
        worst = max(worst, relerr(vb, va))
    end
    worst
end

"""P8's partition of a population's ids by key `kk` of `variant` (`nothing` when P8 has none)."""
function p8_classes(pool, variant)
    f = joinpath(P8DATA, join(pool, "+") * ".jls")
    isfile(f) || return nothing
    recs = deserialize(f).recs
    Dict(r[:id] => r[(variant, :k1, :key)] for r in recs if haskey(r, (variant, :k1, :key)))
end

"""The pairs of ids that `p` puts in one class and `q` does not (both id => key)."""
function split_pairs(p, q)
    cls = Dict{Any, Vector{String}}()
    foreach(((id, k),) -> push!(get!(cls, k, String[]), id), p)
    [(a, b) for v in values(cls) for (i, a) in enumerate(v) for b in v[i + 1:end]
     if haskey(q, a) && haskey(q, b) && q[a] != q[b]]
end

rng = Xoshiro(20261012)
summary = NamedTuple[]
for pool in POOLS
    t0 = time()
    rows = [x for l in readdir(IDATA) if startswith(l, "pipe_")
            for x in deserialize(joinpath(IDATA, l)) if x.set in pool && x.ok]
    cls = Dict{String, Vector{Any}}()
    foreach(x -> push!(get!(cls, x.key, Any[]), x), sort(rows; by = x -> (x.set, x.index)))
    need = Set((x.set, x.index) for v in values(cls) if length(v) > 1 for x in v)
    byset = Dict(set => Dict(r.c.index => r for r in records_of(set)
                             if (set, r.c.index) in need) for set in pool)
    nmerge, worst, bad, cross, maxsize = 0, 0.0, String[], 0, 1
    ctrl, caught = 0, 0                 # negative control: two constants' images swapped
    for v in values(cls)
        length(v) < 2 && continue
        maxsize = max(maxsize, length(v))
        length(unique(x.set for x in v)) > 1 && (cross += 1)
        a = byset[v[1].set][v[1].index]
        for x in v[2:end]
            nmerge += 1
            e = merge_error(a, byset[x.set][x.index], v[1].order, x.order, rng)
            worst = max(worst, e)
            e > 1e-25 && push!(bad, "$(v[1].id) ~ $(x.id): $e")
            if length(x.order) >= 2 && ctrl < 200
                i, j = randperm(rng, length(x.order))[1:2]
                o = copy(x.order)
                o[i], o[j] = o[j], o[i]
                ctrl += 1
                caught += merge_error(a, byset[x.set][x.index], v[1].order, o, rng; n = 1) > 1e-6
            end
        end
    end
    @printf("%-40s %6d laws, %6d classes (largest %d, spanning sets %d); merges %d verified \
            through the oracles: worst %.1e, failures %d (%.0f s)\n", join(pool, "+"),
            length(rows), length(cls), maxsize, cross, nmerge, worst, length(bad), time() - t0)
    ctrl > 0 && println("    negative control (two constants swapped in the map): rates ",
                        "differ in $caught of $ctrl")
    foreach(b -> println("    FALSE MERGE? ", b), first(bad, 10))
    mine = Dict(x.id => x.key for x in rows)
    p8 = Dict{Symbol, Any}()
    for variant in (:redesign, :today)
        theirs = p8_classes(pool, variant)
        theirs === nothing && continue
        common = intersect(keys(mine), keys(theirs))
        a = Dict(i => mine[i] for i in common)
        b = Dict(i => theirs[i] for i in common)
        only_mine, only_p8 = split_pairs(a, b), split_pairs(b, a)
        @printf("    vs P8 %-8s on %d common laws: classes mine %d, P8 %d; pairs merged only \
                by mine %d, only by P8 %d\n", variant, length(common),
                length(unique(values(a))), length(unique(values(b))), length(only_mine),
                length(only_p8))
        foreach(p -> println("        only mine: ", p), first(only_mine, 6))
        foreach(p -> println("        only P8:   ", p), first(only_p8, 6))
        p8[variant] = (common = length(common), mine = length(unique(values(a))),
                       p8 = length(unique(values(b))), only_mine = length(only_mine),
                       only_p8 = length(only_p8))
    end
    push!(summary, (pool = join(pool, "+"), laws = length(rows), classes = length(cls),
                    merges = nmerge, worst, failures = length(bad), cross, p8))
    flush(stdout)
end
serialize(joinpath(IDATA, "classes_$LABEL.jls"), summary)
