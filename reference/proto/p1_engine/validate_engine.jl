# ABOUTME: Validates Engine.bare_rate_law on every plain mechanism: exact agreement with the package's raw
# ABOUTME: polynomials, the BigFloat oracle, detailed balance, speed identifiability, timings and sizes.
const C = joinpath(@__DIR__, "..", "common")
include(joinpath(C, "Oracle.jl")); include(joinpath(C, "Engine.jl"))
using Random, Printf, Statistics, Serialization
setprecision(BigFloat, Oracle.PREC)

const POPS = ["specs", "init_uniuni", "init_unibi", "init_bibi", "init_bibi_pp", "expand_bibi",
              "handwritten"]
const PR = Poly.P{Rational{Int}}
topoly(pp) = Dict{Vector{Pair{Symbol, Int}}, Rational{Int}}(mono => c for (mono, c) in pp)
"""A symbolic polynomial (`symbol => exponent` monomials) in the variables of `V`."""
fromsym(pp, V) = PR(Dict{Poly.Mono, Rational{Int}}(
    Poly.mono((Poly.var!(V, s) => e for (s, e) in mono)...) => c for (mono, c) in pp))
torat(p) = PR(Dict{Poly.Mono, Rational{Int}}(k => c for (k, c) in p.t))

"""Values of every variable of the bare law `L` at step constants `k` and concentrations `c`."""
function vals_of(L, m, k, c)
    vals = ones(BigFloat, length(L.vars.names))
    for (j, vs) in enumerate(L.stepvars), (i, v) in enumerate(vs)
        vals[v] = k[j][i]
    end
    for x in m.metabolites
        vals[L.vars.index[x]] = c[x]
    end
    vals
end
lawv(L, vals) = Poly.evaluate(L.num, vals) / Poly.evaluate(L.den, vals)
relerr(a, b) = Float64(abs(a - b) / abs(b))

"""SS steps whose kf and kr can be scaled together (×3) at constants `k` without changing v."""
function speed_free(L, m, k, c)
    v0 = lawv(L, vals_of(L, m, k, c))
    [j for (j, s) in enumerate(m.steps) if !s.re &&
     abs(lawv(L, vals_of(L, m, [i == j ? 3 .* q : q for (i, q) in enumerate(k)], c)) - v0) <=
     big(10)^-60 * abs(v0)]
end

mechs = [(p, m) for p in POPS for m in Pop.load_population(p) if !m.allosteric]
println("plain mechanisms: ", length(mechs))
for (_, m) in mechs[1:40]       # compile before timing
    try Engine.bare_rate_law(m); Engine.bare_rate_law(m; check = false) catch end
end

rows = Dict{Symbol, Any}[]
for (n, (pop, m)) in enumerate(mechs)
    r = Dict{Symbol, Any}(:pop => pop, :id => m.id, :source => m.source, :family => m.family,
                          :pkg_den => length(m.raw_den), :pkg_num => length(m.raw_num))
    push!(rows, r)
    L = try
        Engine.bare_rate_law(m)
    catch e
        r[:error] = sprint(showerror, e)
        continue
    end
    merge!(r, L.stats)
    r[:t_check] = minimum(@elapsed(Engine.bare_rate_law(m)) for _ in 1:3)
    r[:t_hill] = minimum(@elapsed(Engine.bare_rate_law(m; check = false)) for _ in 1:3)

    # Exact agreement with the package after renaming steps to their group names.
    num, den, dfree = Engine.grouped(L, m)
    V = L.vars
    r[:num_eq] = Poly.symbolic(num, V) == topoly(m.raw_num)
    r[:den_eq] = Poly.symbolic(den, V) == topoly(m.raw_den)
    r[:dfree_eq] = Poly.symbolic(dfree, V) == topoly(m.raw_dfree)
    r[:cross_eq] = torat(num) * fromsym(m.raw_den, V) == fromsym(m.raw_num, V) * torat(den)
    r[:grouped_den] = Poly.nterms(den)
    r[:grouped_num] = Poly.nterms(num)

    # Variables absent from the bare law.
    present = Set(v for p in (L.num, L.den) for mono in keys(p.t) for (v, _) in mono)
    r[:absent_steps] = [j for (j, vs) in enumerate(L.stepvars) if any(v -> v ∉ present, vs)]

    # Numeric: oracle at random bare constants (RE constants consistent on RE cycles), detailed
    # balance at consistent constants, and which SS speeds the law carries.
    rng = Random.Xoshiro(7919 * n)
    errs = Float64[]
    for _ in 1:2
        kc, _, _ = Oracle.consistent_constants(m, rng; Keq = Oracle.logu(rng))
        kb = [s.re ? kc[j] : [Oracle.logu(rng), Oracle.logu(rng)] for (j, s) in enumerate(m.steps)]
        c = Dict{Symbol, BigFloat}(x => Oracle.logu(rng) for x in m.metabolites)
        push!(errs, relerr(lawv(L, vals_of(L, m, kb, c)), Oracle.steady_state_v(m, kb, c)))
        r[:speed_free_bare] = speed_free(L, m, kb, c)
    end
    r[:oracle_err] = maximum(errs)
    Keq = Oracle.logu(rng)
    kc, _, _ = Oracle.consistent_constants(m, rng; Keq)
    c = Oracle.concs_at_ratio(m, rng, Keq)
    neq = Poly.evaluate(L.num, vals_of(L, m, kc, c))
    c[m.products[1]] /= 10
    r[:db] = Float64(abs(neq) / abs(Poly.evaluate(L.num, vals_of(L, m, kc, c))))
    c2 = Dict{Symbol, BigFloat}(x => Oracle.logu(rng) for x in m.metabolites)
    r[:speed_free_thermo] = speed_free(L, m, kc, c2)
end
serialize(joinpath(@__DIR__, "results.jls"), rows)

# ─── Report ───
ok(r) = !haskey(r, :error)
q(v, p) = isempty(v) ? NaN : quantile(v, p)
println("\n== Errors")
for r in rows
    ok(r) || println("  ", r[:id], "  ", first(r[:error], 200), "  ", r[:source])
end
println("\n== Exact agreement with the package's raw polynomials (step variables → group names)")
for pop in [POPS; "all"]
    rs = [r for r in rows if (pop == "all" || r[:pop] == pop) && ok(r)]
    @printf("  %-13s n=%3d  num %3d  den %3d  dfree %3d  num·den' = num'·den %3d\n", pop,
            length(rs), count(r -> r[:num_eq], rs), count(r -> r[:den_eq], rs),
            count(r -> r[:dfree_eq], rs), count(r -> r[:cross_eq], rs))
end
for r in rows
    ok(r) && !(r[:num_eq] && r[:den_eq] && r[:dfree_eq]) &&
        println("  MISMATCH ", r[:id], " num=", r[:num_eq], " den=", r[:den_eq], " dfree=",
                r[:dfree_eq], " cross=", r[:cross_eq], "  ", r[:source])
end
println("\n== Numeric (oracle at bare constants, 2 draws; detailed balance at Γ = Keq)")
for pop in [POPS; "all"]
    rs = [r for r in rows if (pop == "all" || r[:pop] == pop) && ok(r)]
    @printf("  %-13s n=%3d  oracle max rel %.1e (>1e-25: %d)  db max %.1e (>1e-30: %d)\n", pop,
            length(rs), maximum(r[:oracle_err] for r in rs), count(r -> r[:oracle_err] > 1e-25, rs),
            maximum(r[:db] for r in rs), count(r -> r[:db] > 1e-30, rs))
end
println("\n== Timing per derivation (µs; min of 3)")
pk = isfile(joinpath(@__DIR__, "pkg_times.jls")) ?
     deserialize(joinpath(@__DIR__, "pkg_times.jls")) : Dict{String, NamedTuple}()
for pop in [POPS; "all"]
    rs = [r for r in rows if (pop == "all" || r[:pop] == pop) && ok(r)]
    th = [r[:t_hill] for r in rs] .* 1e6; tc = [r[:t_check] for r in rs] .* 1e6
    tp = [pk[r[:id]].t_raw for r in rs if haskey(pk, r[:id])] ./ 1e3
    @printf("  %-13s engine (Hill only) median %7.1f max %8.1f | with edge-sum check median %7.1f max %8.1f | package raw median %7.1f max %8.1f\n",
            pop, median(th), maximum(th), median(tc), maximum(tc),
            isempty(tp) ? NaN : median(tp), isempty(tp) ? NaN : maximum(tp))
end
if !isempty(pk)
    rp = [r for r in rows if ok(r) && haskey(pk, r[:id])]
    ratio = [r[:t_hill] * 1e9 / pk[r[:id]].t_raw for r in rp]
    @printf("  engine/package time ratio: median %.2f, p90 %.2f, max %.2f; totals engine %.3f s, package %.3f s\n",
            median(ratio), quantile(ratio, 0.9), maximum(ratio), sum(r[:t_hill] for r in rp),
            sum(pk[r[:id]].t_raw for r in rp) / 1e9)
    worst = sort(rp; by = r -> -r[:t_hill])[1:5]
    for r in worst
        @printf("    slowest %-16s engine %.0f µs, package %.0f µs, den %d, G %d, cycles %d\n",
                r[:id], r[:t_hill] * 1e6, pk[r[:id]].t_raw / 1e3, r[:den_terms], r[:segments], r[:cycles])
    end
end
println("\n== Sizes")
for pop in [POPS; "all"]
    rs = [r for r in rows if (pop == "all" || r[:pop] == pop) && ok(r)]
    den = [r[:den_terms] for r in rs]; fac = [r[:factor_terms] + r[:sigma_terms] for r in rs]
    multi = count(r -> r[:blocks] >= 2, rs)
    @printf("  %-13s bare den terms median %4d max %4d | factored (block factors + σ) median %4d max %4d | ≥2 blocks %3d of %3d | grouped den max %4d\n",
            pop, median(den), maximum(den), median(fac), maximum(fac), multi, length(rs),
            maximum(r[:grouped_den] for r in rs))
end
allok = [r for r in rows if ok(r)]
ratios = [(r[:factor_terms] + r[:sigma_terms]) / r[:den_terms] for r in allok]
@printf("  factored/expanded size ratio: median %.2f, p10 %.2f, min %.3f; mechanisms where factored < expanded: %d of %d\n",
        median(ratios), quantile(ratios, 0.1), minimum(ratios), count(<(1), ratios), length(allok))
sig_only = [r[:sigma_terms] + r[:dsum_terms] for r in allok]
sig_blk = [r[:sigma_terms] + r[:factor_terms] for r in allok]
dens = [r[:den_terms] for r in allok]
@printf("  totals over %d mechanisms: expanded %d terms; σ factored only (Σ|σ_g| + Σ|D_g|) %d; σ and blocks %d\n",
        length(allok), sum(dens), sum(sig_only), sum(sig_blk))
@printf("  blocks shrink the σ-factored size in %d mechanisms (by %d terms in total); σ shrinks the expanded size in %d\n",
        count(i -> sig_blk[i] < sig_only[i], eachindex(allok)),
        sum(max(0, sig_only[i] - sig_blk[i]) for i in eachindex(allok)),
        count(i -> sig_only[i] < dens[i], eachindex(allok)))
hist = Dict{String, Int}()
for x in ratios
    k = x < 0.5 ? "<0.5" : x < 0.75 ? "0.5-0.75" : x < 1 ? "0.75-1" : x == 1 ? "=1" : ">1"
    hist[k] = get(hist, k, 0) + 1
end
println("  factored/expanded ratio histogram: ", join(["$k: $(get(hist, k, 0))" for k in ("<0.5", "0.5-0.75", "0.75-1", "=1", ">1")], ", "))
println("  max |n_κ| over all mechanisms: ", maximum(r[:max_n] for r in allok),
        "; mechanisms with |n_κ| = 2: ", count(r -> r[:max_n] == 2, allok))
nb = [r[:blocks] for r in allok]
println("  blocks per mechanism: ", join(["$b: $(count(==(b), nb))" for b in sort(unique(nb))], ", "),
        "; non-leaf multi-node blocks (size > 2) in mechanisms with ≥2 blocks: ",
        count(r -> r[:blocks] >= 2 && count(>(2), r[:block_sizes]) >= 2, allok))
bigs = [r for r in allok if r[:den_terms] >= 100]
for r in sort(bigs; by = r -> -r[:den_terms])[1:min(end, 8)]
    @printf("    %-16s den %4d factored %4d blocks %s σ-terms %d\n", r[:id], r[:den_terms],
            r[:factor_terms] + r[:sigma_terms], r[:block_sizes], r[:sigma_terms])
end
println("\n== Numerator: Hill vs edge-sum terms before cancellation")
for pop in [POPS; "all"]
    rs = [r for r in rows if (pop == "all" || r[:pop] == pop) && ok(r)]
    @printf("  %-13s Hill %7d  edge-sum %7d  after cancellation %6d  (Hill/edge-sum %.2f)\n", pop,
            sum(r[:hill_terms] for r in rs), sum(r[:omega_terms] for r in rs),
            sum(r[:num_terms] for r in rs),
            sum(r[:hill_terms] for r in rs) / sum(r[:omega_terms] for r in rs))
end
println("\n== Graph features")
for pop in [POPS; "all"]
    rs = [r for r in rows if (pop == "all" || r[:pop] == pop) && ok(r)]
    μ = [r[:ss] - r[:loops] - r[:segments] + 1 for r in rs]
    @printf("  %-13s segments ≤%2d  cycle rank max %d  cycles max %3d  self-loops %3d (silent %3d)  parallel pairs %3d  RE cycles %3d mech\n",
            pop, maximum(r[:segments] for r in rs), maximum(μ), maximum(r[:cycles] for r in rs),
            sum(r[:loops] for r in rs), sum(length(r[:silent_loops]) for r in rs),
            count(r -> r[:parallel] > 0, rs), count(r -> !isempty(r[:re_closing]), rs))
end
println("\n== Variables absent from the bare law vs the structural prediction (silent self-loops ∪ RE closing steps)")
bad = [r for r in allok if r[:absent_steps] != sort([r[:silent_loops]; r[:re_closing]])]
println("  mechanisms with absent step constants: ", count(r -> !isempty(r[:absent_steps]), allok),
        "; prediction wrong: ", length(bad))
for r in bad[1:min(end, 10)]
    println("    ", r[:id], " absent ", r[:absent_steps], " predicted ", sort([r[:silent_loops]; r[:re_closing]]))
end
println("\n== SS steps whose speed v does not carry (×3 on kf and kr), vs dead-end prediction")
for pop in [POPS; "all"]
    rs2 = [r for r in allok if pop == "all" || r[:pop] == pop]
    @printf("  %-13s thermo-consistent: %3d mech (%3d steps); bare: %3d mech (%3d steps); prediction = thermo in %3d of %3d\n",
            pop, count(r -> !isempty(r[:speed_free_thermo]), rs2), sum(length(r[:speed_free_thermo]) for r in rs2),
            count(r -> !isempty(r[:speed_free_bare]), rs2), sum(length(r[:speed_free_bare]) for r in rs2),
            count(r -> r[:dead_end] == r[:speed_free_thermo], rs2), length(rs2))
end
for r in [r for r in allok if r[:dead_end] != r[:speed_free_thermo]][1:min(end, 10)]
    println("    ", r[:id], " predicted ", r[:dead_end], " thermo ", r[:speed_free_thermo],
            " bare ", r[:speed_free_bare], "  ", r[:source])
end
