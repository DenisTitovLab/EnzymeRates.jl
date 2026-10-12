# ABOUTME: Allosteric kcat check: each corner's analytic MWC limit against the law at large
# ABOUTME: concentrations, and the maximum against the package's sampled kcat on the same model.
include(joinpath(@__DIR__, "lib.jl"))

records = [r for r in load_records() if S.is_allosteric(r.c.mech) && S.admit(r.c.mech) === nothing]
rng = Xoshiro(20261016)
errs, corners, worst_num, cmp, errids, bad = Dict{String, Int}(), 0, 0.0, [], String[], []
zero_corners, zero_worst = 0, 0.0
t = @elapsed for r in records
    m, p = r.c.mech, r.p
    sol = S.solve(m)
    A, vars, terms = try
        A = S.allo_law(m, sol, S.kcat_mask(m))
        A, S.kcat_terms(m, A)...
    catch err
        k = replace(first(split(sprint(showerror, err), ";")), r"with [A-Za-z0-9_, ]+" => "with …")
        errs[k] = get(errs, k, 0) + 1
        push!(errids, "$(p.id): $k")
        continue
    end
    pars = random_params(m, rng)
    ev(x) = maximum(S.evaluate(nn, x) / S.evaluate(dd, x) for (nn, dd) in terms)
    others = [x for x in S.metabolites(m) if !S.is_reactant(m, x)]
    pure(sat) = all(s -> count(in(s.ligands), sat) <= 1, m.sites)
    for (sat, (nn, dd)) in zip(filter(pure, S.all_subsets(others)), terms)
        conc = Dict(x => x in m.substrates ? big(10)^100 : x in sat ? big(10)^10 : big(0)
                    for x in S.metabolites(m))
        x = [vars.conc[i] ? conc[n] : getfield(pars, n) for (i, n) in enumerate(vars.names)]
        k = S.evaluate(nn, x) / S.evaluate(dd, x)
        v = allo_value(A, conc, pars)
        if iszero(k)
            global zero_corners += 1
            global zero_worst = max(zero_worst, Float64(abs(v)))
        else
            e = Float64(abs(v - k) / abs(k))
            e > 1e-6 && push!(bad, (p.id, sat, Float64(k), Float64(v)))
            global worst_num = max(worst_num, e)
        end
        global corners += 1
    end
    # the package's kcat at its own samples, where its model is mine (equal fitted lengths)
    length(S.fitted_params(m)) == length(p.fitted_params) || continue
    names = S.column_names(m, sol)
    for s in p.samples
        isnan(s.kcat) && continue
        vals = pkg_values(p, s)
        kA, kI = pkg_step_constants(p, vals, :A), pkg_step_constants(p, vals, :I)
        kI = [kI[j] === nothing ? kA[j] : kI[j] for j in eachindex(kI)]
        xA, xI = S.edge_constants(r.c, kA), S.edge_constants(r.c, kI)
        col(i) = (q = sol.cols[i]; (q.inactive ? xI : xA)[q.edge][q.slot])
        pk(x) = haskey(vals, x) ? vals[x] : vals[Symbol(replace(String(x), "K_" => "K_A_"; count = 1))]
        mine = (; (names[i] => col(i) for i in sol.indep)...,
                (x => pk(x) for (_, ns) in S.regulator_names(m) for x in ns)...,
                L = vals[:L], Keq = vals[:Keq])
        x = [vars.conc[i] ? big(1) : getfield(mine, n) for (i, n) in enumerate(vars.names)]
        push!(cmp, (set = r.set, id = p.id, mine = Float64(ev(x)), pkg = s.kcat))
    end
end
println("admitted mechanisms: ", length(records), "; kcat errors: ", sum(values(errs); init = 0))
foreach(((k, v),) -> println("    ", v, " × ", k), sort!(collect(errs); by = last, rev = true))
foreach(x -> println("    ", x), errids)
@printf("corners checked against the law at substrates 1e100, corner names 1e10: %d; nonzero limits worst relative %.1e; %d zero limits, largest |v| there %.1e\n",
        corners, worst_num, zero_corners, zero_worst)
println("corners off by more than 1e-6: ", length(bad))
foreach(x -> println("    ", x), first(bad, 15))
ids = unique(x.id for x in cmp)
diff = [x for x in cmp if abs(x.mine - x.pkg) > 1e-9 * abs(x.pkg)]
@printf("package's sampled kcat on the same model: %d samples of %d mechanisms; differ (> 1e-9) in %d samples of %d mechanisms\n",
        length(cmp), length(ids), length(diff), length(unique(x.id for x in diff)))
foreach(x -> println("    ", x, "  ratio ", x.mine / x.pkg), first(diff, 25))
@printf("time %.0f s\n", t)
