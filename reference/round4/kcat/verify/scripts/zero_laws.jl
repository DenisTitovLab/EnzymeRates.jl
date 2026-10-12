# ABOUTME: Rechecks the kcat-zero claims on the probe's populations with independent code: the
# ABOUTME: slice's class against numeric saturation, exact positivity, per-substrate gaps, traps.
include(joinpath(@__DIR__, "vlib.jl"))

"""The probe's populations: every core set, expand_terter, and its 3,000-law ter-ter sample."""
function populations()
    recs = NamedTuple[]
    for set in CORE
        foreach(r -> push!(recs, merge(r, (pop = "core",))), core_records(set))
    end
    foreach(r -> push!(recs, merge(r, (pop = "expand_terter",))), pop_records("expand_terter"))
    picks = deserialize(joinpath(KDATA, "terter_sample.jls"))
    for set in unique(first.(picks))
        idx = [i for (s, i, _) in picks if s == set]
        foreach(r -> push!(recs, merge(r, (pop = "terter_sample",))), pop_records(set, idx))
    end
    recs
end

"""Rates at products and modifiers zero along `t·base` of all substrates, t = 1e4, 1e8, 1e12."""
function scale_scan(nl, nc, subs, base)
    map((1e4, 1e8, 1e12)) do t
        c = zeros(nc)
        c[subs] .= t .* base
        rate(nl, c)
    end
end

"""Positivity at 20 random substrate points (products, modifiers zero), log-uniform 1e-4..1e4."""
function pos_scan(nl, nc, subs, rng)
    vs = map(1:20) do _
        c = zeros(nc)
        foreach(i -> c[i] = 10.0^(-4 + 8 * rand(rng)), subs)
        rate(nl, c)
    end
    minimum(vs), maximum(vs), count(v -> isfinite(v) && v > 0, vs)
end

"""
The trap of plain `m` as substrate index `xi` grows: oracle populations at X × 1e6 and × 1e9
from a random base (others log-uniform 1e-1..1e1, products and modifiers zero). Returns the
trap forms (population not falling, above 1e-9), their share, whether each holds X and a
product or residual, and the exits from the trap classified at zero products.
"""
function trap(m, kv, xi, rng)
    names, subs, prods, mods = roles(m)
    base = zeros(length(names))
    foreach(i -> base[i] = 10.0^(-1 + 2 * rand(rng)), subs)
    pops(g) = (c = big.(base); c[xi] *= g; oracle(m, kv, c)[2])
    p6, p9 = pops(1e6), pops(1e9)
    T = [i for i in eachindex(m.forms) if p9[i] > 1e-9 && p9[i] >= 0.5 * p6[i]]
    X = names[xi]
    zero_name(n) = n in m.products || !S.is_reactant(m, n)
    holdsX(i) = any(x -> x.name === X && !x.modifier, S.mets_of(m, m.forms[i].bound))
    holdsXany(i) = any(x -> x.name === X, S.mets_of(m, m.forms[i].bound))
    holdsP(i) = m.forms[i].res != 0 || any(x -> x.name in m.products, S.mets_of(m, m.forms[i].bound))
    exits = String[]
    kinds = Symbol[]
    for (k, e) in enumerate(m.edges)
        ina, inb = e.a in T, e.b in T
        ina == inb && continue
        cons, rel = ina ? (e.cons, e.rel) : (e.rel, e.cons)       # leaving the trap
        cn, rn = [x.name for x in S.mets_of(m, cons)], [x.name for x in S.mets_of(m, rel)]
        kind = any(zero_name, cn) ? :blocked_takes_up_product :
               isempty(cn) && !isempty(rn) && all(==(X), rn) ? :releases_X : :other
        push!(kinds, kind)
        push!(exits, String(kind) * ": " * S.step_text(m, e))
    end
    (X = X, forms = [S.form_text(m, i) for i in T], share = Float64(sum(p9[T]; init = big(0))),
     all_hold_X = all(holdsXany, T), all_hold_X_reactant = all(holdsX, T),
     all_hold_P = all(holdsP, T), kinds, exits)
end

function analyze(r, rng)
    m = r.m
    base = (pop = r.pop, set = r.set, index = r.index)
    ok, why, sol, L = derive_law(m)
    ok || return merge(base, (ok = false, why = why))
    allo = S.is_allosteric(m)
    cls, vars, terms = kcat_class(m, sol)
    names, subs, prods, mods = roles(m)
    @assert L.vars.names[1:length(names)] == names
    zf = zero_facts(m, L)
    zc = try
        zero_facts(m, allo ? S.allo_law(m, sol, S.kcat_mask(m)) :
                            S.law(m, sol, S.catalytic_core(m)))
    catch e
        (pos = false, num0 = false, den0 = false, gap = Any[nothing for _ in subs], ugap = nothing,
         err = first(sprint(showerror, e), 100))
    end
    nc = length(names)
    pts = map(1:2) do _
        pd = draw_params(L.vars, rng)
        nl = NL(L, pd)
        kc = cls in (:finite, :zero) ? kcat_value(vars, terms, pd) : NaN
        b1 = [10.0^(-1 + 2 * rand(rng)) for _ in subs]
        b2 = [10.0^(-1 + 2 * rand(rng)) for _ in subs]
        s1, s2 = scale_scan(nl, nc, subs, b1), scale_scan(nl, nc, subs, b2)
        pmin, pmax, npos = cls === :zero ? pos_scan(nl, nc, subs, rng) : (NaN, NaN, 0)
        (; kc, s1, s2, pmin, pmax, npos)
    end
    traps = Any[]
    orc = NaN
    if cls === :zero && !allo
        pd = draw_params(L.vars, rng)
        kv = edge_constants(m, sol, pd)
        nl = NL(L, pd)
        # law against oracle at two positive points with products zero
        orc = maximum(1:2) do _
            c = zeros(nc)
            foreach(i -> c[i] = 10.0^(-2 + 4 * rand(rng)), subs)
            vo, _ = oracle(m, kv, big.(c))
            Float64(abs(rate(nl, big.(c)) - vo) / abs(vo))
        end
        for (j, xi) in enumerate(subs)
            g, gc = zf.gap[j], zc.gap[j]
            (g !== nothing && g < 0) || (gc !== nothing && gc < 0) || continue
            push!(traps, merge(trap(m, kv, xi, rng), (core_culprit = gc !== nothing && gc < 0,
                                                      full_culprit = g !== nothing && g < 0)))
        end
    end
    merge(base, (ok = true, allo, n = m.multiplicity, nmods = length(mods), cls, zf, zc, pts,
                 traps, orc, lines = cls === :finite ? String[] : steps_text(m)))
end

t0 = time()
recs = populations()
@printf("%d records loaded in %.0f s\n", length(recs), time() - t0)
flush(stdout)
rng = Xoshiro(20261011)
rows = Any[]
for (i, r) in enumerate(recs)
    push!(rows, analyze(r, rng))
    i % 2000 == 0 && (@printf("  %d done at %.0f s\n", i, time() - t0); flush(stdout))
end
serialize(joinpath(VD, "data", "zero_rows.jls"), rows)
@printf("analysed %d records in %.0f s, maxrss %.0f MB\n", length(rows), time() - t0,
        Sys.maxrss() / 2^20)
