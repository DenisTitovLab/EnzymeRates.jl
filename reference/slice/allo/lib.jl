# ABOUTME: Helpers shared by the allosteric-stage check scripts: evaluating an AlloLaw, the map
# ABOUTME: from slice forms to Pop species, and P5a's two-conformation oracle at a chosen flip.
include(joinpath(@__DIR__, "..", "solve", "lib.jl"))
include(joinpath(PROTO, "p5a", "pkgvals.jl"))           # MWC oracle, package sample values
const Λ = MWC.big1024(10)^40

include(joinpath(@__DIR__, "value.jl"))

"""The Pop species index of each form of converted record `c`."""
function pop_forms(c, p)
    out = zeros(Int, length(c.mech.forms))
    for (j, s) in enumerate(p.steps)
        e = c.mech.edges[c.edge_of[j]]
        a, b = c.flipped[j] ? (e.b, e.a) : (e.a, e.b)
        out[a], out[b] = s.from, s.to
    end
    out
end

"""P5a's factorized formulation-1 oracle (`MWC.v_factorized`) with the conformation flipping
at Pop species `flip` instead of the first free species in Pop order."""
function v_flip(p, flip, kA, kI, reI, activeI, keepI, L, regA, regI, concs; n, lam = Λ)
    setprecision(BigFloat, MWC.PREC) do
        function state(k, re, active, keep)
            forms = sort!(collect(keep))
            loc = Dict(f => i for (i, f) in enumerate(forms))
            length(forms) == 1 && return big(1), big(0)
            x, J = MWC.solve_net(length(forms),
                                 MWC.cat_edges(p, k, re, active, concs, lam, f -> get(loc, f, nothing)))
            x[loc[flip]], J
        end
        πA, JA = state(kA, [s.re for s in p.steps], trues(length(p.steps)), 1:length(p.species))
        πI, JI = state(kI, reI, activeI, keepI)
        wA = MWC.W_factor(p, regA, concs) / πA^n
        wI = BigFloat(L) * MWC.W_factor(p, regI, concs) / πI^n
        (wA * JA + wI * JI) / (wA + wI)
    end
end

"""Oracle inputs of record `r` at all parameter values `allp` (`all_parameters`): the Pop
constants of both conformations (deletion: the inactive law's edges only), the inactive
steps and forms, the regulator constants and the flip species."""
function oracle_inputs(r, allp)
    m, c, p = r.c.mech, r.c, r.p
    mask = S.inactive_mask(m)
    kA = S.pop_constants(c, edge_values(m, allp))
    kI = S.pop_constants(c, edge_values(m, allp; inactive = true, mask))
    actI = [mask[c.edge_of[j]] for j in eachindex(p.steps)]
    flip = pop_forms(c, p)[1]
    keep = Set([flip; [f for (j, s) in enumerate(p.steps) if actI[j] for f in (s.from, s.to)]])
    rn = Dict(S.regulator_names(m))
    regA, regI = Dict{Tuple{Int, Symbol}, Any}(), Dict{Tuple{Int, Symbol}, Any}()
    for (si, site) in enumerate(p.sites), (lig, st) in zip(site.ligands, site.states)
        ns = [getfield(allp, x) for x in rn[lig]]
        regA[(si, lig)] = st === :OnlyI ? nothing : ns[1]
        regI[(si, lig)] = st === :OnlyA ? nothing : ns[end]
    end
    (; kA, kI, actI, keep, flip, regA, regI)
end

"""The oracle's `v/E_total` of record `r` at `allp` and `conc` under deletion semantics."""
function oracle_v(r, allp, conc)
    o = oracle_inputs(r, allp)
    v_flip(r.p, o.flip, o.kA, o.kI, [s.re for s in r.p.steps], o.actI, o.keep, allp.L,
           o.regA, o.regI, conc; n = r.p.cat_multiplicity)
end

relerr(a, b) = Float64(abs(a - b) / max(abs(b), big(1e-300)))
