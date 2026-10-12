# ABOUTME: Allosteric check (a): every allosteric record's MWC law against P5a's two-conformation
# ABOUTME: oracle (deletion semantics) in BigFloat and Float64, and the explicit network for small n.
include(joinpath(@__DIR__, "lib.jl"))

records = [r for r in load_records() if S.is_allosteric(r.c.mech)]
rng = Xoshiro(20261013)
stats = Dict{String, Vector{Any}}()   # set => [n, admitted, worst big (admitted), worst f64, worst rejected]
explicit = Float64[]
fails = String[]
t0 = time()
for r in records
    m, p = r.c.mech, r.p
    st = get!(() -> Any[0, 0, 0.0, 0.0, 0.0], stats, r.set)
    st[1] += 1
    ok = S.admit(m) === nothing
    st[2] += ok
    A, sol = try
        sol = S.solve(m)
        S.allo_law(m, sol), sol
    catch err
        push!(fails, "$(r.set) $(p.id): $(sprint(showerror, err))")
        continue
    end
    for _ in 1:2
        pars = random_params(m, rng)
        allp = S.all_parameters(m, pars)
        u = (0.3 + 1.7 * rand(rng)) * (rand(rng, Bool) ? 1 : -1)
        conc = Oracle.concs_at_ratio(p, rng, pars.Keq * big(10)^u)
        vo = oracle_v(r, allp, conc)
        eb = relerr(allo_value(A, conc, pars), vo)
        ef = relerr(allo_value(A, conc, pars, Float64), vo)
        if ok
            st[3], st[4] = max(st[3], eb), max(st[4], ef)
        else
            st[5] = max(st[5], eb)
        end
        (eb > 1e-25 || (ok && ef > 1e-12)) &&
            push!(fails, "$(r.set) $(p.id): big $eb f64 $ef admitted=$ok")
    end
    # the explicit 1- or 2-protomer network where its flip (the first free Pop species)
    # is the content rule's
    nreg = isempty(p.sites) ? 1 : length(MWC.reg_states(p))
    if ok && p.cat_multiplicity <= 2 && nreg * length(p.species)^p.cat_multiplicity <= 160 &&
       MWC.free_form(p) == pop_forms(r.c, p)[1] && length(explicit) < 400
        pars = random_params(m, rng)
        allp = S.all_parameters(m, pars)
        conc = Oracle.concs_at_ratio(p, rng, pars.Keq * 20)
        o = oracle_inputs(r, allp)
        ve = MWC.v_explicit(p, o.kA, o.kI, [s.re for s in p.steps], o.actI, o.keep, allp.L,
                            o.regA, o.regI, conc; nprot = p.cat_multiplicity, Λ)
        push!(explicit, relerr(allo_value(A, conc, pars), ve))
    end
end
@printf("%-18s %6s %8s %12s %12s %14s\n", "set", "n", "admitted", "worst big", "worst f64",
        "worst rejected")
for (set, st) in sort!(collect(stats))
    @printf("%-18s %6d %8d %12.1e %12.1e %14.1e\n", set, st...)
end
tot = reduce((a, b) -> [a[1] + b[1], a[2] + b[2], max(a[3], b[3]), max(a[4], b[4]),
                        max(a[5], b[5])], values(stats))
@printf("%-18s %6d %8d %12.1e %12.1e %14.1e\n", "all", tot...)
@printf("explicit network (n <= 2, small): %d mechanisms, worst %.1e\n", length(explicit),
        maximum(explicit; init = 0.0))
println("failures: ", length(fails))
foreach(println, first(fails, 30))
@printf("time %.0f s\n", time() - t0)
