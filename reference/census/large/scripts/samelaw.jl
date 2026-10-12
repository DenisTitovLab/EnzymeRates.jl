# ABOUTME: For a sample of rank-deficient fitted candidates, tests whether a candidate with fewer
# ABOUTME: fitted constants in the same run has the same law (Levenberg–Marquardt on coefficients).

include(joinpath(@__DIR__, "census_lib.jl"))
include(joinpath(@__DIR__, "summarize.jl"))

"""The law of slice mechanism `m` as `concentration monomial (text) => coefficient` for `N/H`
and for `D`, in lowest terms, with its variables and fitted names."""
function coefficient_law(m)
    sol = S.solve(m)
    L = S.rate_law(m, sol)
    vars, N, D, Nb = expand_law_checked(L)
    sh = S.conc_shift((N, D), vars)
    cs(P) = Dict(string(k) => v for (k, v) in S.split_conc(P * sh, vars))
    vars, cs(Nb), cs(D), slice_fitted(m, sol)
end

"""Log-coefficient ratios of a coefficient law at log values `lx`, against the reference
monomial `ref` of `D`, with their gradients in the log of each fitted variable (`fpos`)."""
function log_ratios(cN, cD, keysN, keysD, ref, lx, fpos, np)
    g, gr = zeros(np), zeros(np)
    ar = evalgrad!(gr, cD[ref], lx, fpos)
    vals, J = Float64[], zeros(0, np)
    rows = Vector{Float64}[]
    for (cset, ks) in ((cN, keysN), (cD, keysD)), k in ks
        cset === cD && k == ref && continue
        a = evalgrad!(g, cset[k], lx, fpos)
        push!(vals, log(a) - log(ar))
        push!(rows, g ./ a .- gr ./ ar)
    end
    vals, isempty(rows) ? zeros(0, np) : permutedims(reduce(hcat, rows))
end

"""
Whether law `C` (coefficient law of the larger candidate) reproduces law `S` at random points
of `S`: the same concentration skeleton, then Levenberg–Marquardt on `C`'s fitted logs to
match `S`'s log-coefficient ratios, from a few starts. Returns `(same, best_cost)`.
"""
function reproduces(Cl, Sl, rng; points = 2, starts = 6, iters = 300)
    varsC, NC, DC, fC = Cl
    varsS, NS, DS, fS = Sl
    (Set(keys(NC)) == Set(keys(NS)) && Set(keys(DC)) == Set(keys(DS))) || return false, Inf
    keysN, keysD = sort!(collect(keys(NC))), sort!(collect(keys(DC)))
    ref = keysD[1]
    cNC = Dict(k => cpoly(v) for (k, v) in NC); cDC = Dict(k => cpoly(v) for (k, v) in DC)
    cNS = Dict(k => cpoly(v) for (k, v) in NS); cDS = Dict(k => cpoly(v) for (k, v) in DS)
    fposC = zeros(Int, length(varsC.names))
    foreach(((j, n),) -> fposC[varsC.index[n]] = j, enumerate(fC))
    fposS = zeros(Int, length(varsS.names))
    keq = log(10.0) * (2 * rand(rng) - 1)
    worst = 0.0
    for _ in 1:points
        lxS = zeros(length(varsS.names))
        for (i, n) in enumerate(varsS.names)
            varsS.conc[i] && continue
            lxS[i] = n === :Keq ? keq : log(10.0) * (2 * rand(rng) - 1)
        end
        target, _ = log_ratios(cNS, cDS, keysN, keysD, ref, lxS, fposS, 0)
        best = Inf
        for _ in 1:starts
            lx = zeros(length(varsC.names))
            for (i, n) in enumerate(varsC.names)
                varsC.conc[i] && continue
                lx[i] = n === :Keq ? keq : log(10.0) * (4 * rand(rng) - 2)
            end
            λ = 1e-2
            cost(lx) = (r = log_ratios(cNC, cDC, keysN, keysD, ref, lx, fposC, length(fC));
                        (sum(abs2, r[1] .- target), r[1] .- target, r[2]))
            c, res, J = cost(lx)
            for _ in 1:iters
                c < 1e-24 && break
                A = J' * J
                step = (A + λ * Diagonal(diag(A) .+ 1e-12)) \ (J' * res)
                trial = copy(lx)
                for (j, n) in enumerate(fC)
                    trial[varsC.index[n]] -= step[j]
                end
                c2, res2, J2 = cost(trial)
                if c2 < c
                    lx, c, res, J = trial, c2, res2, J2
                    λ = max(λ / 3, 1e-12)
                else
                    λ *= 4
                    λ > 1e12 && break
                end
            end
            best = min(best, c)
            best < 1e-18 && break
        end
        worst = max(worst, best)
    end
    worst < 1e-14, worst
end

name = Symbol(ARGS[1])
K = parse(Int, ARGS[2])
data = deserialize(joinpath(@__DIR__, "..", "data", "$(name).jls"))
info, rounds = data.info, data.rounds
label!(rounds)
fit = [r for recs in rounds for r in recs if r[:fitted] && haskey(r, :smech) &&
       r[:smech] !== nothing && get(r, :rank, -1) >= 0]
byid = Dict(r[:id] => r for recs in rounds for r in recs)
rng = Xoshiro(99)
println("$name: same-law test on up to $K rank-deficient fitted candidates per round")
cache = Dict{Any, Any}()
claw(r) = get!(() -> coefficient_law(r[:smech]), cache, r[:id])
for (ri, recs) in enumerate(rounds)
    P = [r for r in recs if r[:fitted] && :P in get(r, :labels, Symbol[])]
    isempty(P) && continue
    pick = P[sort!(randperm(rng, length(P))[1:min(K, length(P))])]
    tally = Dict{String, Int}()
    t0 = time()
    for r in pick
        # Smaller candidates of the run with the same skeleton and the same rank; the
        # parent first.
        pool = [s for s in fit if s[:n_slice] < r[:n_slice] && s[:rank] == r[:rank] &&
                s[:skeleton] == r[:skeleton]]
        par = get(byid, r[:parent], nothing)
        par !== nothing && (filter!(s -> s !== par, pool); pushfirst!(pool, par))
        verdict = "no smaller candidate with the same skeleton and rank"
        tested = 0
        for s in pool[1:min(8, end)]
            (s[:n_slice] < r[:n_slice] && s[:rank] == r[:rank] &&
             s[:skeleton] == r[:skeleton]) || continue
            tested += 1
            same, c = try
                reproduces(claw(r), claw(s), rng)
            catch e
                println("    error: ", errmsg(e))
                false, Inf
            end
            if same
                verdict = s === par ? "same law as its parent" : "same law as another smaller candidate"
                break
            end
            verdict = "smaller candidates with the same skeleton and rank tested, none reproduced"
        end
        tally[verdict] = get(tally, verdict, 0) + 1
        sub = join(sort!(collect(keys(get(r, :explain, Dict{Symbol, Int}())))), ",")
        @printf("  round %d id %s (%s; slice count %d rank %d): %s [%d tested]\n", ri - 1,
                r[:id], sub, r[:n_slice], r[:rank], verdict, tested)
    end
    @printf("round %d: %d of %d P candidates sampled: %s (%.0f s)\n", ri - 1, length(pick),
            length(P), tally, time() - t0)
    flush(stdout)
end
