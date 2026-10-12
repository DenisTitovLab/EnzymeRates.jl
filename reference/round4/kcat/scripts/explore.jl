# ABOUTME: First look: pipe-row verdicts per set, the kcat_zero ter-ter sample reproduced, and
# ABOUTME: the numeric law against the slice's BigFloat law_value on a few records.
include(joinpath(@__DIR__, "klib.jl"))

rows = [x for f in readdir(IDATA) if startswith(f, "pipe_") for x in deserialize(joinpath(IDATA, f))]
bys = Dict{String, Vector{Any}}()
foreach(x -> push!(get!(bys, population(x.set), Any[]), x), rows)
for (set, xs) in sort!(collect(bys); by = first)
    ok = filter(x -> x.ok, xs)
    kc = [x.kcat isa String ? (occursin("never", x.kcat) ? "never" : occursin("depends", x.kcat) ?
          "depends" : "other:" * first(x.kcat, 40)) : "candidates" for x in ok]
    cnt = Dict{String, Int}()
    foreach(k -> cnt[k] = get(cnt, k, 0) + 1, kc)
    @printf("%-18s records %6d derived %6d  %s\n", set, length(xs), length(ok),
            join(("$k=$v" for (k, v) in sort!(collect(cnt))), ", "))
end
println("pipe row fields: ", keys(first(rows)))

# the kcat_zero sample: rows in readdir order, ok and kcat not a String, randperm(Xoshiro(17))
xs = [x for f in readdir(IDATA) if startswith(f, "pipe_") for x in deserialize(joinpath(IDATA, f))
      if x.ok && !(x.kcat isa String) && startswith(x.set, "init_terter")]
pick = xs[randperm(Xoshiro(17), length(xs))[1:3000]]
println("init_terter non-error laws ", length(xs), "; sample ", length(pick), " first ",
        [(x.set, x.index) for x in pick[1:3]])
serialize(joinpath(KDATA, "terter_sample.jls"), [(x.set, x.index, x.id) for x in pick])

rng = Xoshiro(1)
for set in ("specs", "p5a_adv_cases", "init_bibi")
    recs = records_of(set)
    nerr, worst = 0, 0.0
    t = @elapsed for r in recs[1:min(end, 40)]
        m = r.c.mech
        ok, why, sol, L = derivable(m)
        ok || continue
        cls = kcat_class(m, sol)
        x, cv = params_at(m, sol, rng)
        nl = NLaw(m, L, x)
        for _ in 1:3
            conc = Dict(n => logu10(rng, -2, 2) for n in S.metabolites(m))
            c = [conc[n] for n in S.metabolites(m)]
            a = rate(nl, c)
            b = Float64(law_value(L, conc, x))
            worst = max(worst, abs(a - b) / abs(b))
        end
        if cls[1] == :zero
            println(set, " ", r.p.id, " kcat zero; ", L isa S.Law ? zero_product_limits(m, L) : "allo")
            foreach(l -> println("   ", l), step_lines(m))
        end
    end
    @printf("%s: numeric law vs law_value worst rel %.2e on %d records, %.1f s\n", set, worst,
            min(length(recs), 40), t)
end

# timing on a few ter-ter laws of the sample
ps = deserialize(joinpath(IDATA, "pop_init_terter.part2.jls"))
some = [x for x in pick if x.set == "init_terter.part2"][1:30]
t = @elapsed for x in some
    m = S.Converted(ps[x.index], x.index).mech
    ok, why, sol, L = derivable(m)
    cls = kcat_class(m, sol)
    xx, cv = params_at(m, sol, rng)
    nl = NLaw(m, L, xx)
    sup, at = sup_search(nl, rng)
    kc = cls[1] in (:finite, :zero) ? kcat_value(cls[3], cls[4], xx) : NaN
    f = features(m, S.catalytic_core(m))
    @printf("%-24s %-8s nD %5d sup %.4g kcat %.4g ratio %.4f %s\n", x.id, cls[1], x.nD, sup, kc,
            sup / kc, f)
end
@printf("30 ter-ter laws: %.1f s\n", t)
