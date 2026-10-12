# ABOUTME: Slice-to-package conversion and the oracle face-gap test shared by the census scripts
# ABOUTME: (the limit from inside vs the exact-face RE limit on each one-reactant and full face).
"""The package mechanism of a slice mechanism (plain), on the DSL's placeholder reaction."""
function to_pkg(m::S.Mechanism)
    (subs, prods, steps, group), _ = S.raw(m)
    mods = sort(unique([x.name for s in steps for x in [s.consumed; s.released; s.from.bound; s.to.bound]
                        if x.modifier]))
    met(x) = x.modifier ? ER.CompetitiveInhibitor(x.name) :
             x.name in subs ? ER.Substrate(x.name) : ER.Product(x.name)
    form(f) = ER.Species(ER.Metabolite[met(x) for x in f.bound], f.conf,
                         ER.Residual(ER.Substrate[ER.Substrate(x) for x in f.added],
                                     ER.Product[ER.Product(x) for x in f.subtracted]))
    stp(s) = ER.Step(form(s.from), form(s.to), ER.Metabolite[met(x) for x in s.consumed],
                     ER.Metabolite[met(x) for x in s.released], s.re)
    rxn = ER.EnzymeReaction(
        [[ER.ReactantAtoms(ER.Substrate(s), [:C => length(prods)]) for s in subs];
         [ER.ReactantAtoms(ER.Product(p), [:C => length(subs)]) for p in prods]],
        [ER.RegulatorMults(ER.CompetitiveInhibitor(r), [1], :unspecified) for r in mods], [1])
    ER.Mechanism(rxn, [ER.Step[stp(steps[j]) for j in findall(==(g), group)] for g in 1:maximum(group)])
end

"""Per face (each reactant alone at zero, each side at zero): the largest relative gap
between the exact-face RE limit (3 speed patterns, Λ = 1e30) and the limit from inside
(Λ = 1e40 at 1e-15 of the face metabolites), or :nonunique."""
function face_gaps(mech, rng)
    rxn = ER.reaction(mech)
    subs = [ER.name(x) for x in ER.substrates(rxn)]
    prods = [ER.name(x) for x in ER.products(rxn)]
    _, indep = ER._dependent_param_exprs(mech)
    params = (; (x => 10.0^(2rand(rng) - 1) for x in indep)..., Keq = 10.0^(2rand(rng) - 1), E_total = 1.0)
    steps = mass_action(mech, params)
    allm = unique!([x for s in steps for x in [s.cons; s.rel]])
    base = Dict{Symbol, BigFloat}(x => big(10.0^(2rand(rng) - 1)) for x in allm)
    ref = subs[1]
    scale = abs(rate_at(steps, base, big(10)^30, ones(BigFloat, length(steps)), ref)[1])
    faces = unique!([[[p] for p in prods]; [prods]; [[s] for s in subs]; [subs]])
    out = Dict{Vector{Symbol}, Any}()
    for zs in faces
        c = copy(base); foreach(z -> c[z] = big(10)^-15 * lu10(rng, -1, 1), zs)
        vin = rate_at(steps, c, big(10)^40, ones(BigFloat, length(steps)), ref)[1]
        c0 = copy(base); foreach(z -> c0[z] = big(0), zs)
        gap = big(0); nonu = false
        for k in 1:3
            sp = k == 1 ? ones(BigFloat, length(steps)) : [lu10(rng, -3, 3) for _ in steps]
            v, _, _, nt = rate_at(steps, c0, big(10)^30, sp, ref)
            nt > 1 && (nonu = true; break)
            gap = max(gap, abs(v - vin) / max(abs(vin), big(1e-6) * scale, big(1e-300)))
        end
        out[zs] = nonu ? :nonunique : gap
    end
    out, subs, prods
end
