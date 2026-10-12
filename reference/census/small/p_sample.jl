# ABOUTME: For a sample of rank-deficient (P) candidates, tests whether a smaller candidate of the
# ABOUTME: same run has the same law (equal coefficient images); also spot-checks D key pairs.
# usage: p_sample.jl <reaction names, comma-separated> [per sub-class sample]
include(joinpath(@__DIR__, "lib.jl"))
const NAMES = split(ARGS[1], ",")
const PER = length(ARGS) >= 2 ? parse(Int, ARGS[2]) : 4

ppow(p, k) = k == 0 ? one(typeof(p)) : foldl(*, fill(p, k))

"""Numerator without the Haldane factor and denominator of a slice law, expanded, split by
concentration monomial (named), each with its coefficient polynomial in the constants."""
function coef_map(L)
    if L isa S.Law
        Np, Dp = L.Q, L.D
    else
        A = L
        QI, DI = S.inactive_terms(A)
        T = eltype(values(A.H.terms))
        Lp = S.term(T, S.mono(A.vars.index[:L] => 1))
        WA = foldl(*, (ppow(s[1], s[3]) for s in A.sites); init = one(S.Poly{T}))
        WI = foldl(*, (ppow(s[2], s[3]) for s in A.sites); init = one(S.Poly{T}))
        n = A.n
        Np = ppow(A.fI, n) * A.active.Q * ppow(A.active.D, n - 1) * WA +
             Lp * ppow(A.fA, n) * QI * ppow(DI, n - 1) * WI
        Dp = ppow(A.fI, n) * ppow(A.active.D, n) * WA + Lp * ppow(A.fA, n) * ppow(DI, n) * WI
    end
    vars = L.vars
    name(mo) = join(("$(vars.names[v])^$e" for (v, e) in mo), "*")
    q = Dict(name(mo) => p for (mo, p) in S.split_conc(Np, vars))
    d = Dict(name(mo) => p for (mo, p) in S.split_conc(Dp, vars))
    θv = [i for i in eachindex(vars.names) if !vars.conc[i] && vars.names[i] !== :Keq]
    (q = q, d = d, vars = vars, θv = θv,
     sig = join(sort!(collect(keys(q))), ",") * " / " * join(sort!(collect(keys(d))), ","))
end

"""Value and log-gradient numerators of polynomial `p` at `x` (Float64) for `θpos`."""
function pg64(p, x, θpos, nθ)
    F = 0.0
    G = zeros(nθ)
    for (mo, c) in p.terms
        t = Float64(c)
        for (v, e) in mo
            t *= x[v]^Float64(e)
        end
        F += t
        for (v, e) in mo
            j = get(θpos, v, 0)
            j == 0 || (G[j] += Float64(e) * t)
        end
    end
    F, G
end

"""The normalized log-coefficients of a coefficient map at constants `θ` (Keq fixed), with
their log-gradients: every coefficient over the reference one (D's constant term if any)."""
function logcoefs(cm, θ, keq)
    vars = cm.vars
    x = ones(length(vars.names))
    for (j, v) in enumerate(cm.θv)
        x[v] = θ[j]
    end
    haskey(vars.index, :Keq) && (x[vars.index[:Keq]] = keq)
    θpos = Dict(v => j for (j, v) in enumerate(cm.θv))
    nθ = length(cm.θv)
    ks = [("q", k) for k in sort!(collect(keys(cm.q)))]
    append!(ks, [("d", k) for k in sort!(collect(keys(cm.d)))])
    ref = haskey(cm.d, "") ? ("d", "") : ks[findfirst(k -> k[1] == "d", ks)]
    Fr, Gr = pg64((ref[1] == "q" ? cm.q : cm.d)[ref[2]], x, θpos, nθ)
    out = Tuple{Tuple{String, String}, Float64, Vector{Float64}}[]
    for k in ks
        k == ref && continue
        F, G = pg64((k[1] == "q" ? cm.q : cm.d)[k[2]], x, θpos, nθ)
        push!(out, (k, log(F) - log(Fr), G ./ F .- Gr ./ Fr))
    end
    out
end

"""Whether the law of `t` at constants `θt` lies in the image of `s`'s law (same support):
Levenberg–Marquardt on the log-coefficients from several starts; the best residual."""
function in_image(t, θt, s, keq, rng; starts = 10, iters = 400)
    target = Dict(k => y for (k, y, _) in logcoefs(t, θt, keq))
    ns = length(s.θv)
    best = Inf
    for _ in 1:starts
        u = randn(rng, ns) .* 1.5
        function rj(u)
            lc = logcoefs(s, exp.(clamp.(u, -60, 60)), keq)
            r = [y - target[k] for (k, y, _) in lc]
            J = reduce(vcat, (g' for (_, _, g) in lc))
            r, J
        end
        r, J = rj(u)
        f = sum(abs2, r)
        λ = 1e-2
        for _ in 1:iters
            A = J' * J
            δ = -((A + λ * Diagonal(diag(A) .+ 1e-12)) \ (J' * r))
            r2, J2 = rj(u .+ δ)
            f2 = sum(abs2, r2)
            if isfinite(f2) && f2 < f
                u .+= δ; r, J, f = r2, J2, f2; λ = max(λ / 3, 1e-12)
            else
                λ *= 4
            end
            (f < 1e-26 || λ > 1e12) && break
        end
        best = min(best, maximum(abs, r))
        best < 1e-9 && break
    end
    best
end

"""Same law, both ways (each image holds the other's law at 2 random points), or not."""
function same_law(c, s, rng)
    keq = 10.0^(2rand(rng) - 1)
    fwd = maximum(in_image(c, 10.0 .^ (2 .* rand(rng, length(c.θv)) .- 1), s, keq, rng) for _ in 1:2)
    bwd = maximum(in_image(s, 10.0 .^ (2 .* rand(rng, length(s.θv)) .- 1), c, keq, rng) for _ in 1:2)
    fwd, bwd
end

function main(name)
    d = deserialize(joinpath(CEN, "data", "$name.jls"))
    mechs = deserialize(joinpath(CEN, "data", "$name.mechs.jls"))
    recs = d.recs
    byg = Dict(r[:gid] => r for r in recs)
    fitted = [r for r in recs if r[:status] === :ok && r[:fitted] && haskey(r, :key)]
    println("\n== $name: $(length(fitted)) fitted")
    cmap = Dict{Int, Any}()
    t = @elapsed for r in fitted
        try
            sm = to_slice(mechs[r[:gid]])
            cmap[r[:gid]] = coef_map(S.rate_law(sm, S.solve(sm)))
        catch e
            println("  coef_map error #$(r[:gid]): ", first(errmsg(e), 100))
        end
    end
    @printf("  coefficient maps: %.1f s\n", t)
    bysig = Dict{String, Vector{Int}}()
    for (g, cm) in cmap
        push!(get!(bysig, cm.sig, Int[]), g)
    end
    # the P sample, stratified by primary sub-class
    prim(r) = (dp = get(r, :deficit_parts, Dict{Symbol, Int}());
               isempty(dp) ? :package_only :
               first(k for k in (CAUSES..., :other) if haskey(dp, k)))
    P = [r for r in fitted if get(r, :ranked, false) && r[:rank_pkg] < r[:np]]
    rng = Xoshiro(99)
    strata = Dict{Symbol, Vector{Int}}()
    for r in P
        push!(get!(strata, prim(r), Int[]), r[:gid])
    end
    tally = Dict{Tuple{Symbol, String}, Int}()
    for (k, gs) in sort!(collect(strata); by = first)
        pick = gs[sort!(randperm(rng, length(gs))[1:min(PER, length(gs))])]
        for g in pick
            r = byg[g]
            haskey(cmap, g) || continue
            c = cmap[g]
            cands = [h for h in bysig[c.sig] if h != g && byg[h][:np] < r[:np]]
            sort!(cands; by = h -> (byg[h][:np] != r[:rank_pkg], h != r[:parent], h))
            verdict, with = "no smaller candidate with the same coefficient support", 0
            for h in cands[1:min(end, 25)]
                f, b = same_law(c, cmap[h], rng)
                if f < 1e-8 && b < 1e-8
                    verdict, with = "same law as a smaller candidate", h
                    break
                elseif b < 1e-8 && verdict != "same law as a smaller candidate"
                    verdict, with = "contains a smaller candidate's law, not equal", h
                end
            end
            tally[(k, verdict)] = get(tally, (k, verdict), 0) + 1
            println("  [$k] #$g r$(r[:round]) $(r[:move]) np=$(r[:np]) rank=$(r[:rank_pkg]) ",
                    "($(length(cands)) smaller with its support) → $verdict",
                    with == 0 ? "" : " #$with (np=$(byg[with][:np]), $(with == r[:parent] ? "its parent" : "r$(byg[with][:round]) $(byg[with][:move])"))")
            println("      ", r[:text])
            with == 0 || println("      ", byg[with][:text])
            flush(stdout)
        end
    end
    println("  tally: ", join(("$(k[1]) / $(k[2]): $v" for (k, v) in sort!(collect(tally))), "; "))
    # D spot check: the two laws agree numerically under the canonical correspondence
    Dpairs = Tuple{Int, Int}[]
    first_of = Dict{String, Int}()
    for r in fitted
        haskey(first_of, r[:key]) ? push!(Dpairs, (r[:gid], first_of[r[:key]])) :
                                    (first_of[r[:key]] = r[:gid])
    end
    worst = 0.0
    nchk = 0
    for (a, b) in Dpairs[sort!(randperm(rng, length(Dpairs))[1:min(30, length(Dpairs))])]
        sa, sb = to_slice(mechs[a]), to_slice(mechs[b])
        La, Lb = S.rate_law(sa, S.solve(sa)), S.rate_law(sb, S.solve(sb))
        ca, cb = canonical_named(La), canonical_named(Lb)
        ca.key == cb.key || (println("  key mismatch on recompute #$a #$b"); continue)
        setprecision(BigFloat, 128) do
            vals = Dict(x => logu(rng) for x in cb.order)
            va = Dict(ca.order[i] => vals[cb.order[i]] for i in eachindex(ca.order))
            keq = logu(rng)
            Lv = logu(rng)
            for _ in 1:3
                conc = Dict{Symbol, BigFloat}()
                ev(L, pv) = begin
                    vars = L.vars
                    x = [vars.conc[i] ? get!(() -> logu(rng), conc, vars.names[i]) :
                         vars.names[i] === :Keq ? keq : vars.names[i] === :L ? Lv :
                         pv[vars.names[i]] for i in eachindex(vars.names)]
                    br = law_branches(L)
                    evalbranches(br.N, x) / evalbranches(br.D, x)
                end
                v1, v2 = ev(La, va), ev(Lb, vals)
                worst = max(worst, Float64(abs(v1 - v2) / max(abs(v1), abs(v2))))
            end
        end
        nchk += 1
    end
    println("  D spot check: $nchk pairs, worst relative difference under the key's correspondence $(worst)")
end

foreach(main, NAMES)
