# ABOUTME: Shared package-side helpers for the guard verification: package constants as
# ABOUTME: mass-action steps, law-vs-oracle comparison on faces, builder from plain-data specs.
using EnzymeRates, Serialization, Random, Printf, LinearAlgebra
const ER = EnzymeRates
include(joinpath(@__DIR__, "voracle.jl"))
const R = "/tmp/claude-501/-home-denis-linux--julia-dev-EnzymeRates/12f42995-f5ba-4fda-9510-dd9b87f7d069/scratchpad/round4/guard"

errline(e) = replace(first(sprint(showerror, e), 220), '\n' => ' ')

# ─── per-step mass-action constants of a package mechanism (plain) ──────────────────────
evalx(x::Symbol, P) = P[x]
evalx(x::Real, P) = big(x)
evalx(e::Expr, P) = getfield(Base, e.args[1])((evalx(a, P) for a in e.args[2:end])...)

"""Mass-action steps of `mech` with the constants the package's law uses for `params`
(fitted values, Keq): dependents from `_dependent_param_exprs`, absorbed Wegscheider
names through `_build_wegscheider_rename_map`, RE ratios as `_re_weight_ratio` reads them."""
function mass_action(mech::ER.Mechanism, params)
    dep, _ = ER._dependent_param_exprs(mech)
    rename = ER._build_wegscheider_rename_map(mech)
    P = Dict{Symbol, BigFloat}(k => big(v) for (k, v) in pairs(params))
    for (k, e) in dep
        P[k] = evalx(e, P)
    end
    look(k) = haskey(P, k) ? P[k] : P[rename[k]]
    map(ER._flat_steps(mech)) do (s, g)
        rep = first(ER.steps(mech)[g])
        f, t = string(ER.name(ER.from_species(s))), string(ER.name(ER.to_species(s)))
        c, r = Symbol[ER.name(x) for x in ER.consumed(s)], Symbol[ER.name(x) for x in ER.released(s)]
        if ER.is_equilibrium(s)
            K = look(ER.name(ER.Krapid(rep, :None), mech))
            kf, kr = ER.is_binding(s) ? (big(1), K) : (K, big(1))
        else
            kf, kr = look(ER.name(ER.Kfor(rep, :None), mech)), look(ER.name(ER.Krev(rep, :None), mech))
        end
        MStep(f, t, c, r, ER.is_equilibrium(s), kf, kr)
    end
end

"""Package law vs the oracle with the same constants: at a positive point (Λ → ∞), on the
faces `faces` (law value, interior limit, exact-face range over speed patterns)."""
function compare(label, em, faces, ref, rng; npar = 4)
    mech = ER.Mechanism(em)
    names = ER.fitted_params(em)
    mets = collect(ER.metabolites(em))
    println("\n=== $label  (fitted params: $(length(names)))")
    for ip in 1:npar
        params = (; (x => 10.0^(2rand(rng) - 1) for x in names)..., Keq = 10.0^(2rand(rng) - 1),
                  E_total = 1.0)
        steps = mass_action(mech, params)
        base = Dict{Symbol, BigFloat}(x => big(10.0^(2rand(rng) - 1)) for x in mets)
        nt(c) = (; (x => Float64(c[x]) for x in mets)...)
        vlaw = ER.rate_equation(em, nt(base), params)
        vor = rate_at(steps, base, big(10)^40, ones(BigFloat, length(steps)), ref)[1]
        @printf("  [param %d] positive point: law %.10g  oracle(Λ=1e40) %s  rel.diff %.2e\n", ip,
                vlaw, fmt(vor), abs(vlaw - Float64(vor)) / abs(vlaw))
        for zs in faces
            c0 = copy(base); foreach(z -> c0[z] = big(0), zs)
            lawface = ER.rate_equation(em, nt(c0), params)
            il = interior_limits(steps, base, zs, ref, rng)
            fl = face_limits(steps, base, zs, ref, rng)
            fv = fl.vals
            @printf("    face %-9s law %-12s interior limit %-12s exact face [%s, %s] (nterm %d) max|face-law|/|law| %s\n",
                    join(zs, ","), @sprintf("%.6g", lawface), join(unique(fmt.([r[end] for r in il])), "/"),
                    fmt(minimum(fv)), fmt(maximum(fv)), fl.nterm,
                    isnan(lawface) ? "law NaN" : fmt(maximum(abs.(fv .- lawface)) / abs(lawface)))
        end
    end
end

const RXN = @enzyme_reaction begin
    substrates: A[C], B[N]
    products: P[C], Q[N]
end

# ─── 1. the flagged enumerator children (probe's specs), through the public expand ───────
"""A package mechanism from the probe's plain-data spec, on reaction `rxn`."""
function mybuild(sp, rxn)
    mods = Set(x[1] for s in sp.steps for x in [s[3]; s[4]; s[1][2]; s[2][2]] if x[2])
    met(x) = x[2] ? ER.CompetitiveInhibitor(x[1]) :
             x[1] in sp.substrates ? ER.Substrate(x[1]) : ER.Product(x[1])
    form(f) = ER.Species(ER.Metabolite[met(x) for x in f[2]], f[1],
                         ER.Residual(ER.Substrate[ER.Substrate(x) for x in f[3]],
                                     ER.Product[ER.Product(x) for x in f[4]]))
    stp(s) = ER.Step(form(s[1]), form(s[2]), ER.Metabolite[met(x) for x in s[3]],
                     ER.Metabolite[met(x) for x in s[4]], s[5])
    groups = [ER.Step[stp(sp.steps[j]) for j in findall(==(g), sp.group)] for g in 1:maximum(sp.group)]
    isempty(sp.tags) && return ER.Mechanism(rxn, groups)
    sites = [ER.RegulatorySite(ER.AllostericRegulator[ER.AllostericRegulator(l) for l in s[1]],
                               s[3], collect(Symbol, s[2])) for s in sp.sites]
    ER.AllostericMechanism(rxn, groups, collect(Symbol, sp.tags), sp.multiplicity, sites)
end
dummy_rxn(sp) = begin
    mods = sort(unique([x[1] for s in sp.steps for x in [s[3]; s[4]; s[1][2]; s[2][2]] if x[2]]))
    ER.EnzymeReaction(
        [[ER.ReactantAtoms(ER.Substrate(s), [:C => length(sp.products)]) for s in sp.substrates];
         [ER.ReactantAtoms(ER.Product(p), [:C => length(sp.substrates)]) for p in sp.products]],
        [ER.RegulatorMults(ER.CompetitiveInhibitor(r), [1], :unspecified) for r in mods], [1])
end

