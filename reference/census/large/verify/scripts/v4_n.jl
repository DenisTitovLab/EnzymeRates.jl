# ABOUTME: Class N by growth scans: exact leading exponents of every fitted law along many
# ABOUTME: growth directions, numeric scans of the law on a sample and of the oracle on fewer.
include(joinpath(@__DIR__, "vlib.jl"))
setprecision(BigFloat, 512)

const RX = (:ter_bi_rxn, :ter_ter_rxn, :pyruvate_carboxylase_rxn)

"""The numerator and denominator polynomials of law `L` (plain, or n=1 MWC with no site),
assembled from the law's parts."""
function law_polys(L)
    L isa S.Law && return L.H * L.Q, L.D
    (L.n == 1 && isempty(L.sites)) || error("only n = 1 without sites")
    T = typeof(L.H).parameters[1]
    QI, DI = S.inactive_terms(L)
    Lp = S.term(T, S.mono(L.vars.index[:L] => 1))
    L.H * (L.fI * L.active.Q + Lp * L.fA * QI), L.fI * L.active.D + Lp * L.fA * DI
end

"""Growth directions over `k` variables: each axis, the diagonal, each pair, and `nrand`
random directions in the positive orthant."""
function directions(k, nrand, rng)
    ds = [[i == j ? 1.0 : 0.0 for j in 1:k] for i in 1:k]
    push!(ds, ones(k))
    for i in 1:k, j in i + 1:k
        push!(ds, [(l == i || l == j) ? 1.0 : 0.0 for l in 1:k])
    end
    append!(ds, [rand(rng, k) for _ in 1:nrand])
    ds
end

"""Directions `w` (over the variables `gv`) along which `N/D` grows without bound with the
variables `zv` at zero: the leading exponent `max w·a` of `N`'s surviving terms exceeds
`D`'s."""
function unbounded_dirs(N, D, gv, zv, dirs)
    live(P) = [mo for mo in keys(P.terms) if all(z -> S.mexp(mo, z) == 0, zv)]
    nN, nD = live(N), live(D)
    isempty(nN) && return Vector{Float64}[]
    ex(mo) = [Float64(S.mexp(mo, g)) for g in gv]
    aN, aD = ex.(nN), ex.(nD)
    [w for w in dirs if maximum(a -> sum(w .* a), aN) >
                         maximum(b -> sum(w .* b), aD; init = -Inf) + 1e-9]
end

"""Numeric growth of `f` along `conc = base ∘ t^w` (the names `grow`, `zero` at zero) at
t = 1e10, 1e20, 1e30: the values, and whether they grow by more than half each decade
block."""
function scan(f, base, grow, zero, w)
    vs = map((big(10.0)^10, big(10.0)^20, big(10.0)^30)) do t
        c = copy(base)
        foreach(z -> c[z] = big(0.0), zero)
        for (x, e) in zip(grow, w)
            c[x] = base[x] * t^e
        end
        f(c)
    end
    (grows = abs(vs[3]) > 1.5 * abs(vs[2]) && abs(vs[2]) > 1.5 * abs(vs[1]), vs = vs)
end

rng = Xoshiro(20261011 + 4)
for name in RX
    data = load_run(name)
    t0 = time()
    tally = Dict{Tuple{Int, Symbol}, Int}()
    flagged = Any[]
    fit = [(ri - 1, r) for (ri, recs) in enumerate(data.rounds) for r in recs if r[:fitted]]
    numeric = Set(r[:id] for (_, r) in fit[sort!(randperm(rng, length(fit))[1:300])])
    orc = Set(r[:id] for (_, r) in fit[sort!(randperm(rng, length(fit))[1:12])])
    nnum, nnum_grow = 0, 0
    for (rnd, r) in fit
        m = r[:smech]
        L = S.rate_law(m)
        N, D = law_polys(L)
        idx(xs) = [L.vars.index[x] for x in xs]
        bad = String[]
        for (side, gv, zv) in (("fwd", idx(m.substrates), idx(m.products)),
                               ("rev", idx(m.products), idx(m.substrates)),
                               ("fwd, products>0", idx(m.substrates), Int[]))
            dirs = directions(length(gv), side == "fwd, products>0" ? 0 : 200, rng)
            side == "fwd, products>0" && (dirs = dirs[1:length(gv)])
            u = unbounded_dirs(N, D, gv, zv, dirs)
            isempty(u) || push!(bad, "$side $(length(u)) dirs, e.g. $(round.(u[1]; digits = 2))")
        end
        k = (rnd, isempty(bad) ? :bounded : :unbounded)
        tally[k] = get(tally, k, 0) + 1
        isempty(bad) || push!(flagged, (r[:id], bad, r[:text]))
        if r[:id] in numeric || r[:id] in orc
            vals = random_vals(m, rng)
            base = random_conc(m, rng)
            f(c) = law_value(L, lawx(L, vals, c))
            for (grow, zero) in ((m.substrates, m.products), (m.products, m.substrates))
                for w in directions(length(grow), 2, rng)
                    s = scan(f, base, grow, zero, w)
                    nnum += 1
                    s.grows && (nnum_grow += 1; println("   law grows: $(r[:id]) w=$w ", fmtb.(s.vs)))
                end
            end
            if r[:id] in orc
                cols, x, _, _ = all_columns(m, vals, rng)
                sp = random_speed(m, rng)
                g(c) = first(oracle(m, cols, x, vals, c; Λ = big(10.0)^80, speed = sp))
                out = String[]
                for (side, grow, zero) in (("fwd", m.substrates, m.products),
                                           ("rev", m.products, m.substrates))
                    for w in (ones(length(grow)), rand(rng, length(grow)))
                        s = scan(g, base, grow, zero, w)
                        push!(out, @sprintf("%s %s (%s, %s, %s)", side,
                                            s.grows ? "GROWS" : "bounded", fmtb(s.vs[1]),
                                            fmtb(s.vs[2]), fmtb(s.vs[3])))
                    end
                end
                println("   oracle scan $(r[:id]): ", join(out, "; "))
            end
        end
    end
    @printf("%s (%.0f s): exact leading-exponent test over 200 random + axis/pair/diagonal directions: (round, verdict) => %s\n",
            name, time() - t0, sort!(collect(tally)))
    println("   numeric law scans: $nnum scans on 300 candidates, $nnum_grow growing")
    for (id, bad, txt) in flagged[1:min(5, end)]
        println("   UNBOUNDED $id: ", join(bad, "; "), "\n      ", first(txt, 400))
    end
    flush(stdout)
end
