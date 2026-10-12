# ABOUTME: Maps a package sample of an allosteric Pop.PMech to per-step constants of both
# ABOUTME: conformations (P5a), so the two-conformation oracle can evaluate the package's own model.
include(joinpath(@__DIR__, "MWC.jl"))

"""Every package value of a sample in BigFloat: the Reduced parameters and each dependent
evaluated from its monomial."""
function pkg_values(m, s)
    vals = Dict{Symbol, BigFloat}()
    for (k, v) in s.params; vals[k] = MWC.big1024(v); end
    setprecision(BigFloat, MWC.PREC) do
        for (k, mono) in m.dependent_monomials
            vals[k] = prod((vals[sym]^BigFloat(e) for (sym, e) in mono); init = big(1))
        end
    end
    vals
end

"""Per-step constants (stored orientation) of state `st` from the package's own names:
`step_names_A` for the active state, `step_names_I_pkg` (the pruned inactive graph, with
`inactive_reversed`) for the inactive one; a Pass-2 absorbed name takes its target's value.
A pruned inactive step gets `nothing`."""
function pkg_step_constants(m, vals, st)
    names = st === :A ? m.step_names_A : m.step_names_I_pkg
    ren = Dict(st === :A ? m.pass2_rename_A : m.pass2_rename_I)
    val(n) = vals[get(ren, n, n)]
    map(eachindex(m.steps)) do j
        ns = names[j]
        isempty(ns) && return nothing
        rev = st === :I && !isempty(m.inactive_reversed) && m.inactive_reversed[j]
        s = m.steps[j]
        if s.re
            K = val(ns[1])
            [rev && !Pop.is_binding(s) ? 1 / K : K]
        else
            rev ? [val(ns[2]), val(ns[1])] : [val(ns[1]), val(ns[2])]
        end
    end
end

"""Regulator dissociation constants per (site, ligand) for both states from the package
names `K_A_<lig>reg` / `K_I_<lig>reg` (`nothing` when the ligand is absent)."""
function pkg_regK(m, vals)
    A = Dict{Tuple{Int, Symbol}, Any}(); I = Dict{Tuple{Int, Symbol}, Any}()
    for (si, site) in enumerate(m.sites), (lig, t) in zip(site.ligands, site.states)
        A[(si, lig)] = t === :OnlyI ? nothing : vals[Symbol("K_A_", lig, "reg")]
        I[(si, lig)] = t === :OnlyA ? nothing :
                       t === :EqualAI ? vals[Symbol("K_A_", lig, "reg")] :
                       vals[Symbol("K_I_", lig, "reg")]
    end
    A, I
end

"""The package's own combine evaluated exactly from its exported per-state polynomials."""
function pkg_poly_v(m, vals, concs, regA, regI)
    setprecision(BigFloat, MWC.PREC) do
        allv = merge(vals, Dict(k => MWC.big1024(v) for (k, v) in concs))
        ev(p) = isempty(p) ? big(0) : Oracle.eval_poly(p, allv)
        NA, QA, DA = ev.(m.state_polys_A); NI, QI, DI = ev.(m.state_polys_I)
        n = m.cat_multiplicity
        WA = MWC.W_factor(m, regA, concs); WI = MWC.W_factor(m, regI, concs)
        L = vals[:L]
        num = NA / DA * (QA / DA)^(n - 1) * WA + L * NI / DI * (QI / DI)^(n - 1) * WI
        den = (QA / DA)^n * WA + L * (QI / DI)^n * WI
        num / den
    end
end

isdefined(Main, :Oracle) || include(joinpath(@__DIR__, "..", "common", "Oracle.jl"))
