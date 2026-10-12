# ABOUTME: P4 census over every plain mechanism: constant-term criteria vs the polynomials, reference
# ABOUTME: monomials, units, the Hill factor after the joint solve, and an oracle check of each law.
include(joinpath(@__DIR__, "P4.jl"))
include(joinpath(@__DIR__, "..", "common", "Oracle.jl"))
using Random, Printf, Serialization
using .Poly: Mono, mono, mmul, mpow, mexp, P, nterms, evaluate

const PI = P{Int}
setprecision(BigFloat, 512)
rng = Xoshiro(41)
logu(r) = exp(log(1e-2) + rand(r) * (log(1e2) - log(1e-2)))

rows = Dict{Symbol, Any}[]
for pop in Pop.populations(), m in Pop.load_population(pop)
    m.allosteric && continue
    r = Dict{Symbol, Any}(:id => m.id, :source => m.source, :family => m.family)
    L = Engine.bare_rate_law(m; check = false)
    dpre, shift = P4.den_unshifted(L)
    cpre = P4.cmap(dpre, L.concvars); cpost = P4.cmap(L.den, L.concvars)
    r[:const_pre] = haskey(cpre, Mono()); r[:const_post] = haskey(cpost, Mono())
    r[:shift] = filter(p -> p.second != 0, shift)
    r[:lowdeg] = minimum(P4.degree, keys(cpost))
    cr = P4.criteria(m)
    r[:crit] = cr
    # IndependentParameters law
    NL = P4.law(m)
    D, N = P4.cmap(NL, NL.den), P4.cmap(NL, NL.num)
    S, Pr = P4.cmono(NL, m.substrates), P4.cmono(NL, m.products)
    r[:const_ind] = haskey(D, Mono())
    r[:piS_in_D] = haskey(D, S); r[:piP_in_D] = haskey(D, Pr)
    r[:piS_in_N] = haskey(N, S); r[:piP_in_N] = haskey(N, Pr)
    r[:c0_monomial] = haskey(D, Mono()) && nterms(D[Mono()]) == 1
    r[:cS_monomial] = haskey(D, S) && nterms(D[S]) == 1
    r[:c0_terms] = haskey(D, Mono()) ? nterms(D[Mono()]) : 0
    r[:cS_terms] = haskey(D, S) ? nterms(D[S]) : 0
    r[:multilinear_D] = all(e <= 1 for mo in keys(D) for (_, e) in mo)
    r[:multilinear_N] = all(e <= 1 for mo in keys(N) for (_, e) in mo)
    r[:zero_one_finite] = all(x -> any(mo -> mexp(mo, NL.conc[x]) == 0, keys(D)), m.metabolites)
    r[:homog_ind] = P4.homogeneous(NL)
    r[:homog_re] = P4.homogeneous(P4.law(m; cycles = :re))
    r[:homog_none] = P4.homogeneous(P4.law(m; cycles = :none))
    h = P4.hill(NL)
    r[:hill_factor] = h.Q !== nothing; r[:hill_poly] = h.poly
    r[:hill_positive] = h.positive; r[:hill_keqfree] = h.keq_free
    r[:hill_Q_terms] = h.Q === nothing ? 0 : nterms(h.Q); r[:num_terms] = nterms(NL.num)
    r[:haldane] = [(NL.sol.cols[d.col].name, P4.depclass(NL.sol, d)) for d in NL.sol.deps
                   if d.label === :haldane]
    r[:keq_deps] = count(d -> d.keq != 0, NL.sol.deps)
    r[:keq_in_den] = any(mexp(t, NL.keq) != 0 for t in keys(NL.den.t))
    # oracle at random independents, Keq and concentrations
    yv = Dict(v => logu(rng) for v in values(NL.colvar))
    keq = logu(rng)
    cols = Dict(i => BigFloat(yv[NL.colvar[i]])^NL.colsgn[i] for i in NL.sol.indep)
    k = Thermo.stored_constants(NL.sol, Thermo.all_values(NL.sol, cols, keq))
    concs = Dict(x => logu(rng) for x in m.metabolites)
    vals = zeros(BigFloat, length(NL.vars.names))
    for (v, y) in yv; vals[v] = y; end
    vals[NL.keq] = keq
    for (x, c) in concs; vals[NL.conc[x]] = c; end
    vlaw = evaluate(NL.num, vals) / evaluate(NL.den, vals)
    vor = Oracle.steady_state_v(m, k, concs)
    r[:oracle_rel] = Float64(abs(vlaw - vor) / abs(vor))
    push!(rows, r)
end
serialize(joinpath(@__DIR__, "census.jls"), rows)

n = length(rows)
println("plain mechanisms: ", n)
cnt(f) = count(f, rows)
@printf("oracle vs joint-solved law: max rel %.2e, failures (>1e-20) %d\n",
        maximum(r[:oracle_rel] for r in rows), cnt(r -> r[:oracle_rel] > 1e-20))
println("\n== Part 1: constant term ==")
println("constant term (shifted D): ", cnt(r -> r[:const_post]), "; unshifted D: ", cnt(r -> r[:const_pre]),
        "; IndependentParameters D: ", cnt(r -> r[:const_ind]))
println("nonzero lowest-terms shift: ", cnt(r -> !isempty(r[:shift])))
for (k, f) in ((:c0, r -> r[:crit].c0 == 1), (:c1, r -> r[:crit].c1 == 1), (:c2, r -> r[:crit].c2 == 1),
                (:c3, r -> r[:crit].c3 == 1))
    agree = cnt(r -> f(r) == r[:const_post])
    fp = [r[:id] for r in rows if f(r) && !r[:const_post]]
    fn = [r[:id] for r in rows if !f(r) && r[:const_post]]
    @printf("criterion %s: agrees %d/%d; predicts a constant term wrongly %d %s; misses one %d %s\n",
            k, agree, n, length(fp), first(fp, 8), length(fn), first(fn, 8))
end
allss = filter(r -> :all_ss in r[:family], rows)
println("all-SS mechanisms: ", length(allss), "; c0 agrees ", count(r -> (r[:crit].c0 == 1) == r[:const_post], allss))
println("terminal SCC counts (c3): ", sort(collect(Dict(k => count(r -> r[:crit].c3 == k, rows) for k in 0:4))))
println("single terminal SCC holds the free enzyme (c3): ", cnt(r -> r[:crit].c3 == 1 && r[:crit].free_terminal_c3),
        " of ", cnt(r -> r[:crit].c3 == 1))
println("mechanisms with a segment without exit form (c2): ", cnt(r -> r[:crit].exitless > 0),
        "; with a segment without a least-uptake form (c3): ", cnt(r -> r[:crit].no_minimal > 0),
        " ", [r[:id] for r in rows if r[:crit].no_minimal > 0])
println("segments whose least uptake is not the root's (nonzero minimum): ", cnt(r -> r[:crit].nonzero_min > 0))
println("no constant term, by c3 terminal SCC count: ",
        sort(collect(Dict(k => [r[:id] for r in rows if !r[:const_post] && r[:crit].c3 == k] for k in 0:4))))
bad = [r[:id] for r in rows if r[:crit].c3 >= 1 && r[:lowdeg] < r[:crit].c3 - 1]
println("lowest degree of D below (#terminal SCCs - 1): ", length(bad), " ", bad)
println("lowest degree == #tSCC - 1 when c3 >= 2: ",
        cnt(r -> r[:crit].c3 >= 2 && r[:lowdeg] == r[:crit].c3 - 1), " of ", cnt(r -> r[:crit].c3 >= 2))
println("\n== Part 2: references ==")
println("no constant term: ", cnt(r -> !r[:const_post]), "; of which ΠS in D: ", cnt(r -> !r[:const_post] && r[:piS_in_D]))
println("ΠS not in D: ", [(r[:id], r[:source]) for r in rows if !r[:piS_in_D]])
println("ΠS in D: ", cnt(r -> r[:piS_in_D]), "; ΠP in D: ", cnt(r -> r[:piP_in_D]),
        "; ΠS in N: ", cnt(r -> r[:piS_in_N]), "; ΠP in N: ", cnt(r -> r[:piP_in_N]))
println("constant coefficient a single monomial: ", cnt(r -> r[:c0_monomial]), " of ", cnt(r -> r[:const_ind]),
        "; ΠS coefficient a single monomial: ", cnt(r -> r[:cS_monomial]), " of ", cnt(r -> r[:piS_in_D]))
println("constant coefficient terms (median, max): ",
        let t = sort([r[:c0_terms] for r in rows if r[:const_ind]]); (t[cld(end, 2)], t[end]) end,
        "; ΠS coefficient terms (median, max): ",
        let t = sort([r[:cS_terms] for r in rows if r[:piS_in_D]]); (t[cld(end, 2)], t[end]) end)
println("D multilinear (Cleland form possible): ", cnt(r -> r[:multilinear_D]), "; N multilinear: ", cnt(r -> r[:multilinear_N]))
println("finite with any one metabolite at zero: ", cnt(r -> r[:zero_one_finite]))
println("units homogeneous: IndependentParameters ", cnt(r -> r[:homog_ind]), ", AllParameters+RE ",
        cnt(r -> r[:homog_re]), ", AllParameters groups only ", cnt(r -> r[:homog_none]))
println("\n== Part 3: Hill factor ==")
println("(ΠS − ΠP/Keq) divides N: ", cnt(r -> r[:hill_factor]), "; quotient polynomial in concentrations: ",
        cnt(r -> r[:hill_poly]), "; positive: ", cnt(r -> r[:hill_positive]), "; Keq-free: ", cnt(r -> r[:hill_keqfree]))
println("Keq appears in D: ", cnt(r -> r[:keq_in_den]))
@printf("numerator terms %d -> quotient terms %d\n", sum(r[:num_terms] for r in rows), sum(r[:hill_Q_terms] for r in rows))
cls = Dict{Any, Vector{String}}()
for r in rows, (nm, c) in r[:haldane]
    push!(get!(cls, (c, r[:hill_keqfree]), String[]), r[:id])
end
println("Haldane dependents by class (re, kind, against) and Keq-free quotient:")
for (k, v) in sort(collect(cls); by = x -> string(x[1]))
    println("  ", k, " => ", length(v), "  e.g. ", first(unique(v), 6))
end
println("Haldane dependents per mechanism: ", sort(collect(Dict(k => cnt(r -> length(r[:haldane]) == k) for k in 0:4))))
println("dependents carrying Keq per mechanism (Haldane plus tied members): ", sort(collect(Dict(k => cnt(r -> r[:keq_deps] == k) for k in 0:6))))
println("not Keq-free: ", [r[:id] for r in rows if r[:hill_factor] && !r[:hill_keqfree]][1:min(end, 30)])
