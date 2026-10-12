# ABOUTME: Attacks the K1 key: renamings that keep or move positions, models that differ only in
# ABOUTME: which inhibitor or regulator binds, rewritten inputs, and the rates behind equal keys.
include(joinpath(@__DIR__, "advlib.jl"))

E(b...; kw...) = F(:E, b...; kw...)
key_of(m) = S.canonical(S.rate_law(m)).key

"""`m` with every metabolite name `x` renamed `f(x)`."""
function renamed(m, f)
    (subs, prods, steps, group), kw = S.raw(m)
    met(x) = S.Met(f(x.name), x.modifier)
    form(q) = S.RawForm(q.conf, met.(q.bound), f.(q.added), f.(q.subtracted))
    sites = [S.Site(f.(s.ligands), s.states, s.multiplicity) for s in kw.sites]
    S.Mechanism(f.(subs), f.(prods), [S.RawStep(form(s.from), form(s.to), met.(s.consumed),
                                                met.(s.released), s.re) for s in steps],
                group; kw.tags, kw.multiplicity, sites)
end

"""Rates of two mechanisms at the same concentrations by name and parameters by canonical
position (the key's map), at random values."""
function rates_by_position(m1, m2, rng)
    c1, c2 = S.canonical(S.rate_law(m1)), S.canonical(S.rate_law(m2))
    p1 = Dict{Symbol, Float64}(x => 10.0^(4 * rand(rng) - 2) for x in S.parameters(m1))
    p1[:E_total] = 1.0
    p2 = Dict{Symbol, Float64}(b => p1[a] for (a, b) in zip(c1.order, c2.order))
    p2[:Keq], p2[:E_total] = p1[:Keq], 1.0
    haskey(p1, :L) && (p2[:L] = p1[:L])
    names = union(S.metabolites(m1), S.metabolites(m2))
    conc = Dict(x => 10.0^(4 * rand(rng) - 2) for x in names)
    nt(d, ks) = (; (k => d[k] for k in ks)...)
    v1 = S.rate_equation(S.compile(m1), nt(conc, S.metabolites(m1)), nt(p1, S.parameters(m1)))
    v2 = S.rate_equation(S.compile(m2), nt(conc, S.metabolites(m2)), nt(p2, S.parameters(m2)))
    v1, v2, conc
end

function main()
    rng = Xoshiro(11)
    println("== 1. models that differ only in which inhibitor binds (same structure)")
    uni(inh) = S.Converted(build("uni_$inh", [:S], [:P], [
        st(E(), E(:S); c = [:S], re = true, g = 1), st(E(:S), E(); r = [:P], g = 2),
        st(E(), E(inh); c = [inh], re = true, g = 3)]), 0).mech
    mI, mJ = uni(:I), uni(:J)
    println("  metabolites ", S.metabolites(mI), " vs ", S.metabolites(mJ))
    println("  keys equal: ", key_of(mI) == key_of(mJ))
    v1, v2, conc = rates_by_position(mI, mJ, rng)
    @printf("  rates at one point (I = %.3g, J = %.3g): %.6g vs %.6g\n", conc[:I], conc[:J],
            v1, v2)
    println("== 2. models that differ only in which regulator binds")
    allo(reg) = S.Converted(build("allo_$reg", [:S], [:P], [
        st(E(), E(:S); c = [:S], re = true, g = 1), st(E(:S), E(); r = [:P], g = 2)];
        tags = [:NonequalAI, :NonequalAI], mult = 2,
        sites = [(ligands = [reg], states = [:OnlyA], multiplicity = 1)]), 0).mech
    mR, mT = allo(:R), allo(:T)
    println("  keys equal: ", key_of(mR) == key_of(mT))
    v1, v2, conc = rates_by_position(mR, mT, rng)
    @printf("  rates at one point (R = %.3g, T = %.3g): %.6g vs %.6g\n", conc[:R], conc[:T],
            v1, v2)
    println("== 3. inhibitor on E vs the same inhibitor on E(S) (a different model)")
    unc = S.Converted(build("unc", [:S], [:P], [
        st(E(), E(:S); c = [:S], re = true, g = 1), st(E(:S), E(); r = [:P], g = 2),
        st(E(:S), E(:S, :I); c = [:I], re = true, g = 3)]), 0).mech
    println("  keys equal: ", key_of(mI) == key_of(unc))
    println("== 4. renamings over the populations: substrates renamed to sort after every other name")
    for set in ("init_bibi", "init_bibi_pp", "expand_bibi", "allo_bibi", "p1_enum_inh",
                "allo_terbi")
        recs = skeptic_records(set)
        n = nkey = nsolve = nmeta = 0
        for r in recs[randperm(rng, length(recs))[1:min(150, length(recs))]]
            m = r.c.mech
            k0 = try key_of(m) catch; continue end
            f = x -> x in m.substrates ? Symbol("zz", x) : x
            m2 = try renamed(m, f) catch; continue end
            n += 1
            S.metabolites(m2) == f.(S.metabolites(m)) || (nmeta += 1)
            k2 = key_of(m2)
            nkey += k2 != k0
            s1, s2 = S.solve(m), S.solve(m2)
            nsolve += (s1.indep, [(q.col, q.terms, q.keq) for q in s1.deps]) !=
                      (s2.indep, [(q.col, q.terms, q.keq) for q in s2.deps])
        end
        @printf("  %-12s renamed %3d: positions moved %d, key changed %d, solve changed %d\n",
                set, n, nmeta, nkey, nsolve)
    end
    println("== 5. a rename that keeps metabolites(m) but moves the table (one record shown)")
    for set in ("init_bibi",), r in skeptic_records(set)
        m = r.c.mech
        f = x -> x in m.substrates ? Symbol("zz", x) : x
        m2 = renamed(m, f)
        key_of(m2) == key_of(m) && continue
        println("  ", r.p.id, ": fitted ", S.fitted_params(m), "\n   renamed: ",
                S.fitted_params(m2))
        break
    end
end

main()
