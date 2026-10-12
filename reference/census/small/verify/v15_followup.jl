# ABOUTME: Follow-ups: the :OnlyA family of #19519 with exponents from component heights, the rank of
# ABOUTME: uni_bi #119 in BigFloat, and the full-rank uni_bi pairs that share a law but not a key.
include(joinpath(@__DIR__, "vlib.jl"))
include(joinpath(@__DIR__, "moracle.jl"))
rng = Xoshiro(151515)

"""Exponents e_s = h(bound component) - h(unbound component) from longest-path heights on the
:OnlyA binding arcs between non-:OnlyA components (nothing when the arcs are not a DAG)."""
function height_expo(ms)
    nf = length(ms.forms)
    kept = [ms.tags[s.g] !== :OnlyA for s in ms.steps]
    p = collect(1:nf)
    f(i) = (while p[i] != i; i = p[i]; end; i)
    foreach(k -> kept[k] && (p[f(ms.steps[k].a)] = f(ms.steps[k].b)), eachindex(ms.steps))
    comp = [f(i) for i in 1:nf]
    arcs = [(comp[s.a], comp[s.b]) for (k, s) in enumerate(ms.steps) if !kept[k] && plain_binding(ms, s)]
    any(a -> a[1] == a[2], arcs) && return nothing
    h = Dict(c => 0 for c in unique(comp))
    for _ in 1:length(h) + 1
        changed = false
        for (a, b) in arcs
            h[b] < h[a] + 1 && (h[b] = h[a] + 1; changed = true)
        end
        changed || return [h[b] - h[a] for (a, b) in arcs]
    end
    nothing
end

d, mechs = loadrun("bi_bi_rxn.r2_100000")
ms = mstruct(mechs[19519])
ex = height_expo(ms)
println("#19519 height exponents: ", ex)
for k in 1:3
    A = aconsts(ms, rng)
    Ip = iconsts(ms, A, rng)
    xp = Ip.draw(rng)
    If = iconsts(ms, A, rng; family = true, expo = ex, pin = Ip.pins(xp))
    if If === nothing
        println("  draw $k: no family even with height exponents")
        continue
    end
    L = big(lu10(rng, -1, 1))
    σ = [big(lu10(rng, -1, 1)) for _ in 1:(length(ms.steps) + 1)]
    lp = Ip.lnk(xp); l1 = If.lnk(If.draw(rng)); l2 = If.lnk(If.draw(rng))
    for (lab, zs) in [("positive", Symbol[]); faces(ms)]
        c = Dict{Symbol, BigFloat}(x => big(lu10(rng, -1, 1)) for x in ms.mets)
        foreach(z -> c[z] = big(0), zs)
        vp = oracle_rate(ms, A, Ip, lp, L, c, big(10)^30, σ)
        v1 = oracle_rate(ms, A, If, l1, L, c, big(10)^30, σ)
        v2 = oracle_rate(ms, A, If, l2, L, c, big(10)^30, σ)
        @printf("  draw %d %-15s package model %.6g family %.6g / %.6g (spread %.1e, vs package %.1e)\n", k, lab,
                Float64(vp), Float64(v1), Float64(v2), Float64(rel(v1, v2)), Float64(max(rel(v1, vp), rel(v2, vp))))
    end
end

# the island samples and admitted samples again, with exponents from heights (the admitted ones
# where random exponents failed would show here)
SAMPLES = deserialize(joinpath(VER, "samples.jls"))
nok = 0; ntot = 0
for name in ["uni_bi_rxn", "bi_bi_rxn", "bi_bi_pp_rxn", "bi_bi_rxn.r2_100000", "bi_bi_pp_rxn.r2_1200"]
    dd, mm = name == "bi_bi_rxn.r2_100000" ? (d, mechs) : loadrun(name)
    for g in SAMPLES[name][:W_unflagged_allo]
        msg = mstruct(mm[g])
        e = height_expo(msg)
        A = aconsts(msg, rng); Ip = iconsts(msg, A, rng); xp = Ip.draw(rng)
        If = e === nothing ? nothing : iconsts(msg, A, rng; family = true, expo = e, pin = Ip.pins(xp))
        global ntot += 1
        If === nothing || (global nok += 1)
    end
end
println("admitted samples with a consistent family from height exponents: $nok of $ntot")

# uni_bi #119: rank in Complex{BigFloat}
setprecision(BigFloat, 256)
du, mu = loadrun("uni_bi_rxn")
pl = plaw(mu[119])
println("\nuni_bi #119 ", du.recs[119][:text], "  np $(length(pl.names)) census rank $(du.recs[119][:rank_pkg])")
for k in 1:3
    θ = [big(lu10(rng, -1, 1)) for _ in pl.names]
    keq = big(lu10(rng, -1, 1))
    n = length(θ)
    J = zeros(BigFloat, n + 12, n)
    h = big(10)^-60
    for r in 1:(n + 12)
        c = Dict{Symbol, BigFloat}(x => big(lu10(rng, -1, 1)) for x in pl.mets)
        v = real(prate(pl, θ, keq, c))
        for j in 1:n
            θc = complex.(θ); θc[j] += im * θ[j] * h
            J[r, j] = imag(prate(pl, θc, keq, c)) / h / v
        end
    end
    sv = svdvals(Float64.(J))
    println("  draw $k: Float64 singular values of the BigFloat log-Jacobian (relative): ",
            join((@sprintf("%.1e", s / sv[1]) for s in sv), " "))
    for (lo, hi) in ((-2.0, 2.0),)
        rk = cs_rank(pl, rng; lo, hi)
        println("  draw $k: Float64 complex-step at 10^[$lo,$hi]: rank $(rk.rank) kept $(rk.kept) dropped $(rk.dropped)")
    end
end

# the full-rank uni_bi pairs with one law and two keys
for (x, y) in ((126, 122), (127, 123), (134, 130), (136, 132))
    a, b = plaw(mu[x]), plaw(mu[y])
    keq = lu10(rng, -1, 1)
    ab = maximum(fit_into(a, [lu10(rng, -1, 1) for _ in a.names], b, keq, rng; starts = 30, iters = 600) for _ in 1:3)
    ba = maximum(fit_into(b, [lu10(rng, -1, 1) for _ in b.names], a, keq, rng; starts = 30, iters = 600) for _ in 1:3)
    @printf("\nuni_bi #%d (np %d rank %d) vs #%d (np %d rank %d): worst x-in-y %.1e, y-in-x %.1e; keys equal %s\n",
            x, du.recs[x][:np], du.recs[x][:rank_pkg], y, du.recs[y][:np], du.recs[y][:rank_pkg], ab, ba,
            du.recs[x][:key] == du.recs[y][:key])
    println("  #$x ", du.recs[x][:text], "\n  #$y ", du.recs[y][:text])
    for (g, m) in ((x, mu[x]), (y, mu[y]))
        s = ER.rate_equation_string(ER.compile_mechanism(m))
        for l in split(s, '\n')
            startswith(l, "(; S") && continue
            println("   #$g| ", first(l, 900))
        end
    end
end
