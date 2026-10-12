# ABOUTME: The skeptic's per-law checks: oracle, compiled code, printed forms, dependents and
# ABOUTME: kcat of one derived law against the independent route (indep.jl).
include(joinpath(@__DIR__, "indep.jl"))

"""`text` evaluated as Julia in a fresh module with `params` and `concs` bound; the value of
`v`. The coefficient form puts `v` first, so its `v` line moves to the end."""
function eval_text(text, params, concs; coefficient = false)
    mod = Module()
    Core.eval(mod, :(params = $params))
    Core.eval(mod, :(concs = $concs))
    lines = split(text, "\n")
    if coefficient
        iv = findfirst(l -> startswith(l, "v = "), lines)
        im = findfirst(l -> startswith(l, "# The coefficients"), lines)
        lines = [["(; $(join(keys(params), ", "))) = params",
                  "(; $(join(keys(concs), ", "))) = concs"];
                 lines[im + 1:end]; lines[iv + 1:im - 1]; lines[iv]]
    end
    Core.eval(mod, Meta.parseall(join(lines, "\n")))
    Core.eval(mod, :v)
end

"""The checks of one derived record `r` (law `d` from `derive`)."""
function check(r, d, rng)
    p, c, m = r.p, r.c, r.c.mech
    allo = S.is_allosteric(m)
    tags = allo ? p.group_tags : Symbol[]
    sys = my_system(p, tags)
    mine = my_solution(sys, rng)
    sol = d.solution
    count_ok = mine.nfree_law == length(sol.indep)
    lkeq = (4 * big(rand(rng)) - 2) * log(big(10))
    keq = exp(lkeq)
    vals = my_sample(sys, mine, rng, lkeq)
    kA = step_values(p, sys, vals, :A)
    kI = allo ? step_values(p, sys, vals, :I) : kA
    params = slice_params(c, sol, kA, kI, keq, rng)
    u = (0.3 + 1.7 * rand(rng)) * (rand(rng, Bool) ? 1 : -1)
    conc = Oracle.concs_at_ratio(p, rng, keq * big(10)^u)
    for x in S.metabolites(m)
        haskey(conc, x) || (conc[x] = Oracle.logu(rng))
    end
    concs = (; (x => conc[x] for x in S.metabolites(m))...)
    flip, keepI = allo ? my_inactive(p, tags) : (0, falses(0))
    function oracle(cc; Λ = big(10)^60)
        allo || return Oracle.steady_state_v(p, kA, cc; Λ)
        regvals = Dict(k => getfield(params, k) for k in keys(params))
        my_mwc(p, kA, kI, keepI, flip, params.L, my_W(m, regvals, cc, :A),
               my_W(m, regvals, cc, :I), m.multiplicity, cc; Λ)
    end
    vo = oracle(conc)
    em = S.compile(m)
    vB = S.rate_equation(em, concs, params)
    vF = S.rate_equation(em, map(Float64, concs), map(Float64, params))
    a = prod(BigFloat(conc[x]) for x in m.substrates)
    b = prod(BigFloat(conc[x]) for x in m.products) / keq
    cond = Float64((a + b) / abs(a - b))
    vP = eval_text(S.rate_equation_string(m), params, concs)
    vC = eval_text(S.coefficient_form(m), params, concs; coefficient = true)
    # every per-step constant of the law from all_parameters against the sample
    allp = S.all_parameters(m, params)
    fwd = S.orient(m)
    worst_dep = 0.0
    for (inactive, names) in ((false, S.constant_names(m)),
                              (true, allo ? S.constant_names(m; inactive = true) : nothing))
        names === nothing && continue
        mask = inactive ? S.inactive_mask(m) : trues(length(m.edges))
        for k in eachindex(m.edges), s in eachindex(names[k])
            mask[k] || continue
            g = m.group[k]
            nm = inactive && m.tags[g] === :EqualAI ? S.constant_names(m)[k][s] : names[k][s]
            haskey(allp, nm) || continue          # a constant no column of the law reads
            want = column_value(c, fwd, S.Col(k, s, inactive && m.tags[g] !== :EqualAI),
                                kA, kI)
            worst_dep = max(worst_dep, relerr(BigFloat(getfield(allp, nm)), want))
        end
    end
    # kcat against the oracle at saturation on the steps outside dead-end branches
    kc = "-"
    kerr = NaN
    core = my_core(p)
    if allo
        core .|= [s.re && Pop.is_free(p.species[s.from]) && Pop.is_free(p.species[s.to])
                  for s in p.steps]
    end
    function oracle_core(cc)
        Λ = big(10)^130
        if !allo
            x, J, _ = my_network(p, kA, core, sort!(unique([f for (j, s) in enumerate(p.steps)
                                  if core[j] for f in (s.from, s.to)])), cc; Λ)
            return J
        end
        regvals = Dict(k => getfield(params, k) for k in keys(params))
        keepIc = keepI .& core
        formsA = sort!(unique([f for (j, s) in enumerate(p.steps) if core[j]
                               for f in (s.from, s.to)]))
        xA, JA, locA = my_network(p, kA, core, formsA, cc; Λ)
        formsI = sort!(unique([flip; [f for (j, s) in enumerate(p.steps) if keepIc[j]
                                      for f in (s.from, s.to)]]))
        xI, JI, locI = my_network(p, kI, keepIc, formsI, cc; Λ)
        n = m.multiplicity
        wA = my_W(m, regvals, cc, :A) / xA[locA[flip]]^n
        wI = params.L * my_W(m, regvals, cc, :I) / xI[locI[flip]]^n
        (wA * JA + wI * JI) / (wA + wI)
    end
    kslice = try
        S.kcat(m, params)
    catch err
        err isa ErrorException || rethrow()
        err.msg
    end
    nmods = count(x -> !S.is_reactant(m, x), S.metabolites(m))
    if nmods > 8                              # 2^nmods corners: too many for the oracle
        kc = "skip(mods=$nmods)"
    elseif kslice isa String
        # never saturates: the oracle still grows tenfold up; depends on growth: two
        # substrate ratios at saturation disagree
        vs1 = oracle_corners(m, Xoshiro(1), oracle_core; X = big(10)^20, Y = big(10)^10)
        vs2 = oracle_corners(m, Xoshiro(1), oracle_core; X = big(10)^21, Y = big(10)^10)
        vs3 = oracle_corners(m, Xoshiro(2), oracle_core; X = big(10)^20, Y = big(10)^10)
        grows = maximum(Float64(abs(b2 / b1)) for (b1, b2) in zip(vs1, vs2))
        spread = maximum(relerr(b1, b3) for (b1, b3) in zip(vs1, vs3))
        kc = occursin("never saturates", kslice) ? (grows > 2 ? "grows-ok" : "grows-REFUTED") :
             occursin("depends on how", kslice) ? (spread > 1e-6 ? "depends-ok" :
                                                   "depends-REFUTED") : kslice
    else
        vk = maximum(oracle_corners(m, rng, oracle_core))
        kerr = relerr(BigFloat(kslice), vk)
        # both vanish at saturation (the oracle's rate decays like 1/X): compare to the rate
        both0 = abs(BigFloat(kslice)) + abs(vk) <= big(10)^-15 * abs(vo)
        kc = kerr < 1e-9 ? "kcat-ok" : both0 ? "kcat-ok(0)" : "kcat-DIFF"
    end
    (; count_ok, mine = mine.nfree_law, slice = length(sol.indep), cond,
     errB = relerr(vB, vo), errF = Float64(abs(vF - vo) / abs(vo)), errP = relerr(vP, vo),
     errC = relerr(vC, vo), worst_dep, kc, kerr)
end

