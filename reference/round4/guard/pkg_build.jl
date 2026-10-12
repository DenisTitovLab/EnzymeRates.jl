# ABOUTME: Builds a package mechanism from a plain-data spec, with the reaction the DSL would
# ABOUTME: build (dummy balanced atoms, modifier names as competitive inhibitors).
"""The package mechanism of a spec, with the reaction the DSL would build: substrates and
products with balanced dummy atoms, every name that binds as a modifier a competitive
inhibitor (`[1]`), catalytic multiplicities `[1]`."""
function build(sp)
    subs, prods = sp.substrates, sp.products
    all_mets = [x for s in sp.steps for x in [s[3]; s[4]; s[1][2]; s[2][2]]]
    mods = sort(unique([x[1] for x in all_mets if x[2]]))
    met(x) = x[2] ? ER.CompetitiveInhibitor(x[1]) :
             x[1] in subs ? ER.Substrate(x[1]) : ER.Product(x[1])
    form(f) = ER.Species(ER.Metabolite[met(x) for x in f[2]], f[1],
                         ER.Residual(ER.Substrate[ER.Substrate(x) for x in f[3]],
                                     ER.Product[ER.Product(x) for x in f[4]]))
    step(s) = ER.Step(form(s[1]), form(s[2]), ER.Metabolite[met(x) for x in s[3]],
                      ER.Metabolite[met(x) for x in s[4]], s[5])
    rxn = ER.EnzymeReaction(
        [[ER.ReactantAtoms(ER.Substrate(s), [:C => length(prods)]) for s in subs];
         [ER.ReactantAtoms(ER.Product(p), [:C => length(subs)]) for p in prods]],
        [ER.RegulatorMults(ER.CompetitiveInhibitor(r), [1], :unspecified) for r in mods],
        [1])
    G = maximum(sp.group)
    groups = [ER.Step[step(sp.steps[j]) for j in findall(==(g), sp.group)] for g in 1:G]
    isempty(sp.tags) && return ER.Mechanism(rxn, groups)
    sites = [ER.RegulatorySite(ER.AllostericRegulator[ER.AllostericRegulator(l) for l in s[1]],
                               s[3], collect(Symbol, s[2])) for s in sp.sites]
    ER.AllostericMechanism(rxn, groups, collect(Symbol, sp.tags), sp.multiplicity, sites)
end

