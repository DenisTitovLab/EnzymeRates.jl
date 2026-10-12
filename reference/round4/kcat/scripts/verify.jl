# ABOUTME: Re-checks the sup > kcat laws at fresh parameter points against the BigFloat oracle, and
# ABOUTME: searches harder for overshoot in steady-state random-order mechanisms without abortive forms.
include(joinpath(@__DIR__, "klib.jl"))

conc_dict(m, c) = Dict{Symbol, BigFloat}(n => big(c[i]) for (i, n) in enumerate(S.metabolites(m)))

rows = [x for l in ("core", "expand_terter", "terter_sample")
        for x in deserialize(joinpath(KDATA, "rows_$l.jls")) if x.ok]
over = filter(x -> x.cls === :finite && x.ratio > 1.01, rows)
ssub = filter(x -> x.cls === :finite && !x.allo && x.feat.ss_nt_sub > 0 && x.feat.ss_nt_prod == 0, rows)
need = Set((x.set, x.index) for x in [over; ssub])
recs = Dict{Tuple{String, Int}, Any}()
for set in unique(first.(need))
    if set in CORE_SETS || set == "expand_terter"
        foreach(r -> (r.set, r.c.index) in need && (recs[(r.set, r.c.index)] = r), records_of(set))
    else
        ps = deserialize(joinpath(IDATA, "pop_$set.jls"))
        for (s, i) in need
            s == set && (recs[(s, i)] = (set = s, p = ps[i], c = S.Converted(ps[i], i)))
        end
        ps = nothing
        GC.gc()
    end
end
println(length(over), " sup > 1.01 kcat laws, ", length(ssub), " plain laws whose SS non-turnover cycles avoid product-holding forms")

# 1. fresh points: overshoot again, the oracle at the argmax, and v at saturating substrates
nver, worst, lims, nlaws = 0, 0.0, Float64[], 0
for x in over
    r = recs[(x.set, x.index)]
    m = r.c.mech
    _, _, sol, L = derivable(m)
    _, _, vars, terms = kcat_class(m, sol)
    rng = Xoshiro(hash((x.set, x.index, :verify)))
    hit = false
    for _ in 1:5
        xx, cv = params_at(m, sol, rng)
        nl = NLaw(m, L, xx)
        kc = kcat_value(vars, terms, xx)
        sup, at = sup_search(nl, rng)
        sup > 1.01 * kc || continue
        hit = true
        global nver += 1
        c = zeros(nl.nc)
        foreach(((j, v),) -> c[j] = v, at)
        o = oracle_value(r, sol, xx, cv, conc_dict(m, c))
        global worst = max(worst, Float64(abs(rate(nl, c) - o) / abs(o)))
        cl = zeros(nl.nc)
        foreach(j -> cl[j] = 1e12 * (1 + 9 * rand(rng)), nl.subs)
        push!(lims, rate(nl, cl) / kc)
    end
    global nlaws += hit
end
@printf("fresh points: overshoot again in %d of %d laws (%d points of %d); oracle vs law at the argmax, worst rel %.1e\n",
        nlaws, length(over), nver, 5 * length(over), worst)
@printf("v at substrates ~1e12 over kcat at those points: min %.4g, median %.4g, max %.4g\n",
        minimum(lims), median(lims), maximum(lims))

# 2. harder search in random-order SS mechanisms without SS cycles through product-holding forms
function wide_params(m, sol, rng)
    names = S.column_names(m, sol)
    x = Dict{Symbol, BigFloat}(n => big(10.0^(-3 + 6 * rand(rng))) for n in S.parameters(m))
    cv = Dict{Int, BigFloat}(i => x[names[i]] for i in sol.indep)
    for d in sol.deps
        cv[d.col] = prod((S.pow(cv[f], e) for (f, e) in d.terms); init = S.pow(x[:Keq], d.keq))
    end
    x, cv
end
function hard_search(m; npts = 30, rng)
    _, _, sol, L = derivable(m)
    _, _, vars, terms = kcat_class(m, sol)
    best = 0.0
    for _ in 1:npts
        xx, _ = wide_params(m, sol, rng)
        nl = NLaw(m, L, xx)
        sup, _ = sup_search(nl, rng; mods = false, lo = -6.0, hi = 6.0, nrand = 300 * length(nl.subs))
        best = max(best, sup / kcat_value(vars, terms, xx))
    end
    best
end
ratios = [hard_search(recs[(x.set, x.index)].c.mech; rng = Xoshiro(hash((x.set, x.index, :hard)))) for x in ssub]
@printf("random-order SS without product-holding SS cycles: %d laws, 30 points each (constants 1e-3..1e3, box 1e-6..1e6): sup > 1.01 kcat in %d; max ratio %.4f\n",
        length(ssub), count(>(1.01), ratios), maximum(ratios))
for x in ssub[1:3]
    println("  ", x.id, " (", x.set, "/", x.index, ")")
    foreach(l -> println("      ", l), step_lines(recs[(x.set, x.index)].c.mech))
end

# 3. the textbook random-order steady-state bi-bi with every step its own kinetic group
F(bound...) = S.RawForm(:E, [S.Met(b, false) for b in bound], Symbol[], Symbol[])
M(x) = [S.Met(x, false)]
none = S.Met[]
steps = [S.RawStep(F(), F(:A), M(:A), none, false), S.RawStep(F(), F(:B), M(:B), none, false),
         S.RawStep(F(:A), F(:A, :B), M(:B), none, false), S.RawStep(F(:B), F(:A, :B), M(:A), none, false),
         S.RawStep(F(:A, :B), F(:P, :Q), none, none, false), S.RawStep(F(:P, :Q), F(:Q), none, M(:P), true),
         S.RawStep(F(:Q), F(), none, M(:Q), true)]
m = S.Mechanism([:A, :B], [:P, :Q], steps, collect(1:7))
println("textbook random-order SS bi-bi (untied):")
foreach(l -> println("      ", l), step_lines(m))
r3 = hard_search(m; npts = 300, rng = Xoshiro(7))
@printf("  300 points: max sup/kcat %.4f\n", r3)
