# ABOUTME: Per-law census of one population: kcat class (core and full law), zero-kcat analysis,
# ABOUTME: sup search at zero products and the reference-rate gauge; usage: analyze.jl <label>
include(joinpath(@__DIR__, "klib.jl"))

const LABEL = ARGS[1]

"""The records of a label: every core set, the ter-ter exports, or a ter-ter sample."""
function label_records(label)
    label == "core" && return [r for set in CORE_SETS for r in records_of(set)]
    label in CORE_SETS && return records_of(label)
    label == "expand_terter" && return records_of("expand_terter")
    picks = if label == "terter_sample"
        deserialize(joinpath(KDATA, "terter_sample.jls"))
    elseif label == "terter_nounique"
        rows = [x for k in 1:11 for x in deserialize(joinpath(IDATA, "pipe_terter$k.jls"))
                if x.ok && x.kcat isa String && occursin("depends on how", x.kcat)]
        [(x.set, x.index, x.id) for x in rows[randperm(Xoshiro(4), length(rows))[1:300]]]
    else
        error("unknown label $label")
    end
    out = []
    for set in unique(first.(picks))
        ps = deserialize(joinpath(IDATA, "pop_$set.jls"))
        for (s, i, _) in picks
            s == set && push!(out, (set = set, p = ps[i], c = S.Converted(ps[i], i)))
        end
        ps = nothing
        GC.gc()
    end
    out
end

"""kcat's class on the whole law (no dead-end branch removed): the plain `kcat_terms` with every
edge, or the MWC one on the full `allo_law`."""
function kcat_class_full(m, sol)
    vars, terms = try
        all_edges = trues(length(m.edges))
        S.is_allosteric(m) ? S.kcat_terms(m, S.allo_law(m, sol, all_edges)) :
                             S.kcat_terms(m, sol, all_edges)
    catch err
        msg = err isa ErrorException ? err.msg : sprint(showerror, err)
        cls = occursin("never saturates", msg) ? :infinite :
              occursin("depends on how", msg) ? :nounique : :other
        return (cls, first(msg, 160), nothing, nothing)
    end
    (all(t -> iszero(t[1]), terms) ? :zero : :finite, "", vars, terms)
end

"""Concentrations by name (BigFloat) for the oracle."""
conc_dict(m, c) = Dict{Symbol, BigFloat}(n => big(c[i]) for (i, n) in enumerate(S.metabolites(m)))

"""The forms holding the enzyme as substrate `X` (index `xi`) grows from 1e7 to 1e8 times its
base value (populations of the active network from the oracle): those whose population does
not fall; their share, and the exits from that set at zero products and modifiers."""
function trap_of(r, m, kA, cbase, xi)
    pf = pop_forms(r.c, r.p)
    X = S.metabolites(m)[xi]
    pops(g) = (c = copy(cbase); c[xi] *= g; Oracle.steady_state_populations(r.p, kA, conc_dict(m, c)))
    p7, p8 = pops(1e7), pops(1e8)
    T = [i for i in eachindex(m.forms) if p8[pf[i]] > 1e-30 && p8[pf[i]] > 0.5 * p7[pf[i]]]
    share = Float64(sum((p8[pf[i]] for i in T); init = big(0)))
    zeroname(nm) = nm in m.products || !S.is_reactant(m, nm)
    nx, blocked, other = 0, String[], String[]
    for (k, e) in enumerate(m.edges)
        ina, inb = e.a in T, e.b in T
        ina == inb && continue
        cons, rel = ina ? (e.cons, e.rel) : (e.rel, e.cons)
        cn = [y.name for y in S.mets_of(m, cons)]
        rn = [y.name for y in S.mets_of(m, rel)]
        if any(zeroname, cn)
            push!(blocked, S.step_text(m, e))
        elseif isempty(cn) && !isempty(rn) && all(==(X), rn)
            nx += 1
        else
            push!(other, S.step_text(m, e))
        end
    end
    (X = X, forms = [S.form_text(m, i) for i in T], holdsX = all(i -> holds(m, i, X) > 0, T),
     holdsP = all(i -> holds_product(m, i), T), share, nx, blocked, other)
end

"""State-resolved populations of an allosteric record at concentrations `c` (by name): the
active network (every step) and the inactive one (the inactive law's steps and the flip form),
each scaled by its conformation's share, as the two-conformation oracle (`v_flip`) weighs them.
Entries `(state, Pop species, population)`."""
function mwc_pops(r, sol, x, cv, c)
    m, p = r.c.mech, r.p
    mask = S.inactive_mask(m)
    kA, kI = step_constants(r, sol, cv, false), step_constants(r, sol, cv, true)
    actI = [mask[r.c.edge_of[j]] for j in eachindex(p.steps)]
    flip = pop_forms(r.c, p)[1]
    keepI = Set([flip; [f for (j, s) in enumerate(p.steps) if actI[j] for f in (s.from, s.to)]])
    rn = Dict(S.regulator_names(m))
    regA, regI = Dict{Tuple{Int, Symbol}, Any}(), Dict{Tuple{Int, Symbol}, Any}()
    for (si, site) in enumerate(p.sites), (lig, st) in zip(site.ligands, site.states)
        regA[(si, lig)] = st === :OnlyI ? nothing : x[rn[lig][1]]
        regI[(si, lig)] = st === :OnlyA ? nothing : x[rn[lig][end]]
    end
    re = [s.re for s in p.steps]
    setprecision(BigFloat, MWC.PREC) do
        function state(k, active, keep)
            forms = sort!(collect(keep))
            loc = Dict(f => i for (i, f) in enumerate(forms))
            length(forms) == 1 && return forms, [big(1)]
            xs, _ = MWC.solve_net(length(forms),
                                  MWC.cat_edges(p, k, re, active, c, Λ, f -> get(loc, f, nothing)))
            forms, xs
        end
        fA, xA = state(kA, trues(length(p.steps)), 1:length(p.species))
        fI, xI = state(kI, actI, keepI)
        n = p.cat_multiplicity
        wA = MWC.W_factor(p, regA, c) / xA[findfirst(==(flip), fA)]^n
        wI = BigFloat(x[:L]) * MWC.W_factor(p, regI, c) / xI[findfirst(==(flip), fI)]^n
        a = wA / (wA + wI)
        [[(:A, f, a * xA[i]) for (i, f) in enumerate(fA)];
         [(:I, f, (1 - a) * xI[i]) for (i, f) in enumerate(fI)]]
    end
end

"""The trap of an allosteric law as substrate `X` (index `xi`) grows from 1e7 to 1e8 times its
base value: the (conformation, form) entries whose population does not fall, their share, and
the exits from the trap inside each conformation's network at zero products and modifiers
(a trap holding the flip form can also leave by the conformational change)."""
function mwc_trap_of(r, sol, x, cv, cbase, xi)
    m = r.c.mech
    pf = pop_forms(r.c, r.p)
    X = S.metabolites(m)[xi]
    pops(g) = (c = copy(cbase); c[xi] *= g; mwc_pops(r, sol, x, cv, conc_dict(m, c)))
    p7, p8 = pops(1e7), pops(1e8)
    held = [e8[3] > 1e-30 && e8[3] > 0.5 * e7[3] for (e7, e8) in zip(p7, p8)]
    T = [(e[1], findfirst(==(e[2]), pf)) for (e, h) in zip(p8, held) if h]
    share = Float64(sum((e[3] for (e, h) in zip(p8, held) if h); init = big(0)))
    zeroname(nm) = nm in m.products || !S.is_reactant(m, nm)
    mask = S.inactive_mask(m)
    nx, blocked, other = 0, String[], String[]
    for st in (:A, :I)
        Ts = [i for (s, i) in T if s === st]
        for (k, e) in enumerate(m.edges)
            st === :I && !mask[k] && continue
            ina, inb = e.a in Ts, e.b in Ts
            ina == inb && continue
            cons, rel = ina ? (e.cons, e.rel) : (e.rel, e.cons)
            cn = [y.name for y in S.mets_of(m, cons)]
            rn = [y.name for y in S.mets_of(m, rel)]
            if any(zeroname, cn)
                push!(blocked, "$st: " * S.step_text(m, e))
            elseif isempty(cn) && !isempty(rn) && all(==(X), rn)
                nx += 1
            else
                push!(other, "$st: " * S.step_text(m, e))
            end
        end
        1 in Ts && push!(other, "$st: conformational change at the flip form")
    end
    (X = X, forms = ["$s:" * S.form_text(m, i) for (s, i) in T],
     holdsX = all(((s, i),) -> holds(m, i, X) > 0, T),
     holdsP = all(((s, i),) -> holds_product(m, i), T), share, nx, blocked, other)
end

"""Task 2 at one parameter point of a kcat-zero law: positivity at random substrate points,
growth scans, the oracle at the grown points, and the trap of each substrate that drives
`v` to zero."""
function zero_analysis(r, m, sol, nl, x, cv, rng)
    vpos = map(1:20) do _
        c = zeros(nl.nc)
        foreach(i -> c[i] = logu10(rng, -4, 4), nl.subs)
        rate(nl, c)
    end
    base, scans = growth_scans(nl, rng)
    slope(sc) = log10(sc[end] / sc[end - 1])
    to0(sc) = sc[end] < 1e-3 * maximum(sc) && slope(sc) < -0.5
    drives = [to0(sc) for sc in scans]                 # per substrate, then all together
    kA = step_constants(r, sol, cv, false)
    orc = Float64[]
    for (j, xi) in enumerate([nl.subs; 0]), g in (1.0, 1e8)
        c = copy(base)
        xi == 0 ? (c[nl.subs] .*= g) : (c[xi] *= g)
        o = oracle_value(r, sol, x, cv, conc_dict(m, c))
        push!(orc, Float64(abs(rate(nl, c) - o) / abs(o)))
    end
    traps = [nl.allo ? mwc_trap_of(r, sol, x, cv, base, xi) : trap_of(r, m, kA, base, xi)
             for (j, xi) in enumerate(nl.subs) if drives[j]]
    (vmin = minimum(vpos), vmax = maximum(vpos), npos = count(v -> v > 0 && isfinite(v), vpos),
     slopes = slope.(scans), drives, orc = maximum(orc), traps)
end

function analyze(r)
    m = r.c.mech
    base = (set = r.set, index = r.c.index, id = r.p.id)
    ok, why, sol, L = derivable(m)
    ok || return merge(base, (ok = false, why = why))
    rng = Xoshiro(hash((r.set, r.c.index)))
    allo = S.is_allosteric(m)
    cls, msg, vars, terms = kcat_class(m, sol)
    fcls, fmsg, fvars, fterms = kcat_class_full(m, sol)
    core = S.catalytic_core(m)
    feat = features(m, core)
    zp = allo ? nothing : zero_product_limits(m, L)
    zpc = allo ? nothing : zero_product_limits(m, S.law(m, sol, core))
    pts = map(1:3) do _
        x, cv = params_at(m, sol, rng)
        nl = NLaw(m, L, x)
        kc = cls in (:finite, :zero) ? kcat_value(vars, terms, x) : NaN
        kf = fcls in (:finite, :zero) ? kcat_value(fvars, fterms, x) : NaN
        g = gauge_rates(nl, rng)
        fin = filter(isfinite, g)
        sup, at = sup_search(nl, rng)
        sups, ats = isempty(nl.mods) ? (sup, at) : sup_search(nl, rng; mods = false)
        z2 = cls === :zero ? zero_analysis(r, m, sol, nl, x, cv, rng) : nothing
        (; kc, kf, sup, at, sups, ats, nnan = count(isnan, g), ninf = count(isinf, g),
         maxabs = maximum(abs, fin; init = 0.0), maxabs0 = maximum(abs, filter(isfinite, g[1:10]); init = 0.0),
         minabs = minimum(abs, fin; init = Inf), z2)
    end
    nD = L isa S.Law ? S.nterms(L.D) : S.nterms(L.active.D)
    ratio = cls === :finite ? maximum(p.sup / p.kc for p in pts) : NaN
    ratios = cls === :finite ? maximum(p.sups / p.kc for p in pts) : NaN
    keep_lines = cls !== :finite || fcls !== cls || ratio > 1.01 ||
                 any(p -> p.nnan > 0 || p.ninf > 0 || !(0 < p.maxabs < Inf), pts)
    merge(base, (ok = true, allo, n = m.multiplicity, nsubs = length(m.substrates),
                 mods = [n for n in S.metabolites(m) if !S.is_reactant(m, n)], nD, cls, msg,
                 fcls, fmsg, feat, ncore_removed = count(!, core), zp, zpc, pts, ratio, ratios,
                 lines = keep_lines ? step_lines(m; core) : String[]))
end

recs = label_records(LABEL)
println(LABEL, ": ", length(recs), " records")
flush(stdout)
rows = NamedTuple[]
t0 = time()
for (i, r) in enumerate(recs)
    push!(rows, analyze(r))
    i % 500 == 0 && (@printf("  %d done at %.0f s\n", i, time() - t0); flush(stdout))
end
serialize(joinpath(KDATA, "rows_$LABEL.jls"), rows)
ok = filter(x -> x.ok, rows)
cnt(f) = (d = Dict{Any, Int}(); foreach(x -> (k = f(x); d[k] = get(d, k, 0) + 1), ok); sort!(collect(d); by = string))
@printf("%s: %d records, %d derivable; %.0f s, maxrss %.0f MB\n", LABEL, length(rows), length(ok),
        time() - t0, Sys.maxrss() / 2^20)
println("kcat classes (core): ", cnt(x -> x.cls))
println("full-law classes:    ", cnt(x -> x.fcls))
println("sup > 1.01 kcat:     ", count(x -> x.cls === :finite && x.ratio > 1.01, ok), " of ",
        count(x -> x.cls === :finite, ok))
println("gauge failures:      ", count(x -> any(p -> p.nnan > 0 || p.ninf > 0 || !(0 < p.maxabs < Inf), x.pts), ok))
