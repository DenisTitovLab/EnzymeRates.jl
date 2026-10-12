# ABOUTME: The skeptic's adversarial mechanisms: each built once, run through derive and checked
# ABOUTME: against the independent route; usage: adv.jl (log on stdout, rows in skeptic/data).
include(joinpath(@__DIR__, "advlib.jl"))

include(joinpath(@__DIR__, "adv_cases.jl"))

"""Denominators of the exponents in the law's polynomials."""
function denominators(L)
    ps = L isa S.Law ? [L.Q, L.D] :
         [L.active.Q, L.active.D, L.fA, L.fI, S.inactive_terms(L)...]
    sort!(unique(denominator(e) for P in ps for mo in keys(P.terms) for (_, e) in mo))
end

"""Whether `Keq` appears in `Q` (plain) or in either state's `Q`."""
keq_in_q(L) = (k = L.vars.index[:Keq];
    any(mo -> any(t -> t.first == k, mo),
        keys((L isa S.Law ? L : L.active).Q.terms)))

const ARTIFACTS = [r"\b1 \^", r"\^1\b(?![/0-9])", r"\* 1\b(?![.0-9])", r"\b1 \*", r"L \* 1\b"]

"""Per-call time and allocations of the compiled law through a function barrier."""
@noinline call(em, c, p) = S.rate_equation(em, c, p)
function timing(em, c, p)
    call(em, c, p)
    allocs = @allocated call(em, c, p)
    best = Inf
    acc = 0.0
    for _ in 1:5
        t = @elapsed for _ in 1:10_000
            acc += call(em, c, p)
        end
        best = min(best, t / 10_000)
    end
    isfinite(acc) || error("non-finite")
    allocs, 1e9 * best
end

function main()
    rng = Xoshiro(20261010)
    rows = NamedTuple[]
    for cs in CASES
        p = try
            cs.thunk()
        catch err
            println(rpad(cs.name, 30), "BUILD ERROR ", sprint(showerror, err))
            continue
        end
        r = try
            (set = "adv", p = p, c = S.Converted(p, 0))
        catch err
            msg = sprint(showerror, err)
            ok = cs.expect isa Pair && cs.expect.first === :error && occursin(cs.expect.second, msg)
            println(rpad(cs.name, 30), ok ? "OK " : "UNEXPECTED ", "constructor error: ",
                    first(msg, 200))
            push!(rows, (name = cs.name, outcome = :ctor_error, ok, msg))
            continue
        end
        m = r.c.mech
        d = try
            S.derive(m)
        catch err
            msg = sprint(showerror, err)
            kind = err isa ErrorException ? "error" : string(typeof(err))
            ok = cs.expect isa Pair && cs.expect.first === :reject &&
                 occursin(cs.expect.second, msg)
            println(rpad(cs.name, 30), ok ? "OK " : "UNEXPECTED ", "derive ", kind, ": ",
                    first(msg, 260))
            push!(rows, (name = cs.name, outcome = :rejected, ok, msg))
            continue
        end
        if cs.expect !== :derive
            println(rpad(cs.name, 30), "UNEXPECTED derived; expected ", cs.expect)
        end
        x = try
            check(r, d, rng)
        catch err
            println(rpad(cs.name, 30), "CHECK ERROR ", first(sprint(showerror, err), 400))
            continue
        end
        L = d.law
        sol = d.solution
        dens = denominators(L)
        nh = count(q -> q.label === :haldane, sol.deps)
        texts = [S.rate_equation_string(m), S.coefficient_form(m)]
        arts = [string(a) for a in ARTIFACTS if any(t -> occursin(a, t), texts)]
        concs = (; (y => 10.0^(4 * rand(rng) - 2) for y in S.metabolites(m))...)
        params = (; (y => (y === :E_total ? 1.0 : 10.0^(4 * rand(rng) - 2))
                     for y in S.parameters(m))...)
        em = S.compile(m)
        allocs, ns = timing(em, concs, params)
        # one consistent sample (independent route) for the zero and ε checks
        tags = p.group_tags
        sys = my_system(p, tags); mine = my_solution(sys, rng)
        keq = big(2.5)
        vals = my_sample(sys, mine, rng, log(keq))
        kA = step_values(p, sys, vals, :A)
        kI = S.is_allosteric(m) ? step_values(p, sys, vals, :I) : kA
        prm = slice_params(r.c, sol, kA, kI, keq, rng)
        cc = Oracle.concs_at_ratio(p, rng, keq / 7)
        for y in S.metabolites(m)
            haskey(cc, y) || (cc[y] = Oracle.logu(rng))
        end
        flip, keepI = S.is_allosteric(m) ? my_inactive(p, tags) : (0, falses(0))
        function orc(c2)
            S.is_allosteric(m) || return Oracle.steady_state_v(p, kA, c2; Λ = big(10)^60)
            rv = Dict(k => getfield(prm, k) for k in keys(prm))
            my_mwc(p, kA, kI, keepI, flip, prm.L, my_W(m, rv, c2, :A), my_W(m, rv, c2, :I),
                   m.multiplicity, c2)
        end
        # zero concentrations, one at a time: finite and equal to the oracle
        zero_bad = String[]
        for y in S.metabolites(m)
            c2 = copy(cc); c2[y] = BigFloat(10)^-30          # the oracle's "zero" (see indep.jl)
            cz = (; (z => Float64(c2[z]) for z in S.metabolites(m))...)
            v = S.rate_equation(em, cz, map(Float64, prm))
            vo = orc(c2)
            ok = isfinite(v) && (abs(vo) < big(10)^-40 ? abs(v) < 1e-30 :
                                 abs(v - vo) <= 1e-11 * abs(vo))
            ok || push!(zero_bad, @sprintf("%s(%.3g vs %.3g)", y, v, Float64(vo)))
        end
        # ε-limit of decision 6 (multiplicity 1, no sites)
        epsr = NaN
        if S.is_allosteric(m) && m.multiplicity == 1 && isempty(m.sites)
            ccn = (; (y => cc[y] for y in S.metabolites(m))...)
            vlaw = S.rate_equation(em, ccn, prm)
            errs = [relerr(eps_limit(p, tags, kA, kI, prm.L, cc, BigFloat(10)^(-e)), vlaw)
                    for e in (10, 20, 40)]
            epsr = errs[end]
            @printf("    ε-limit: rel. error at ε = 1e-10, 1e-20, 1e-40: %.1e %.1e %.1e\n", errs...)
        end
        key = d.key
        @printf("%-30s %s fit %d haldane %d dens %s KeqQ %s | B %.1e F %.1e P %.1e C %.1e dep %.1e | %s | %d allocs %.1f ns | zeros bad %s | art %s | %s\n",
                cs.name, x.count_ok ? "count-ok" : "COUNT-DIFF($(x.mine)/$(x.slice))",
                length(sol.indep), nh, dens, keq_in_q(L), x.errB, x.errF, x.errP, x.errC,
                x.worst_dep, x.kc, allocs, ns, isempty(zero_bad) ? "-" : join(zero_bad, ","),
                isempty(arts) ? "-" : join(arts, " "), first(key.key, 12))
        push!(rows, merge((name = cs.name, outcome = :derived, ok = true, dens, nh,
                           keqQ = keq_in_q(L), allocs, ns, zero_bad, arts, epsr,
                           key = key.key, fitted = S.fitted_params(m)), x))
        flush(stdout)
    end
    serialize(joinpath(SCR, "slice", "skeptic", "data", "adv.jls"), rows)
    println("cases: ", length(CASES), "; derived ", count(x -> x.outcome === :derived, rows))
end

main()
