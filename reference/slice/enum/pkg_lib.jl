# ABOUTME: Shared package-side helpers of the enumeration stage: the reactions of the parent
# ABOUTME: sets and the conversions between package mechanisms and Pop.PMech structures.
using EnzymeRates
const ER = EnzymeRates

const bi_bi_pp_rxn = @enzyme_reaction begin
    substrates: A[CX], B[N]
    products: P[C], Q[NX]
end
const bi_bi_pp_inh = @enzyme_reaction begin
    substrates: A[CX], B[N]
    products: P[C], Q[NX]
    dead_end_inhibitors: I, A, Q
end
const bi_bi_allo_rxn = @enzyme_reaction begin
    substrates: A[C], B[N]
    products: P[C], Q[N]
    oligomeric_state: 2
end
const bi_bi_allo_inh = @enzyme_reaction begin
    substrates: A[C], B[N]
    products: P[C], Q[N]
    dead_end_inhibitors: I, A, Q
    oligomeric_state: 2
end
const uni_uni_with_reg = @enzyme_reaction begin
    substrates: S[C]
    products: P[C]
    dead_end_inhibitors: I
end
const bi_bi_with_inh = @enzyme_reaction begin
    substrates: A[C], B[N]
    products: P[C], Q[N]
    dead_end_inhibitors: I
end

_role(::ER.Substrate) = :S
_role(::ER.Product) = :P
_role(::ER.CompetitiveInhibitor) = :Inh
_role(::ER.AllostericRegulator) = :Reg
_pmets(xs) = Pop.PMet[Pop.PMet(ER.name(x), _role(x)) for x in xs]

"""The structure of a package mechanism as a Pop.PMech (no package outputs)."""
function to_pmech(m, id)
    groups = ER.steps(m)
    flat = [s for g in groups for s in g]
    gof = [g for (g, grp) in enumerate(groups) for _ in grp]
    species = ER._forms(groups)
    idx = Dict(sp => i for (i, sp) in enumerate(species))
    rxn = ER.reaction(m)
    allo = m isa ER.AllostericMechanism
    Pop.PMech(; id, source = id,
        substrates = Symbol[ER.name(x) for x in ER.substrates(rxn)],
        products = Symbol[ER.name(x) for x in ER.products(rxn)],
        species = [Pop.PSpecies(_pmets(ER.bound(sp)), ER.conformation(sp),
                                Symbol[ER.name(a) for a in ER.added(ER.residual(sp))],
                                Symbol[ER.name(p) for p in ER.subtracted(ER.residual(sp))],
                                ER.name(sp)) for sp in species],
        steps = [Pop.PStep(idx[ER.from_species(s)], idx[ER.to_species(s)],
                           _pmets(ER.consumed(s)), _pmets(ER.released(s)),
                           ER.is_equilibrium(s)) for s in flat],
        groups = [findall(==(g), gof) for g in eachindex(groups)],
        allosteric = allo,
        cat_multiplicity = allo ? ER.catalytic_multiplicity(m) : 1,
        group_tags = allo ? copy(ER.cat_allo_states(m)) : Symbol[],
        sites = allo ? [Pop.PSite(Symbol[ER.name(l) for l in ER.ligands(s)], ER.multiplicity(s),
                                  copy(ER.allo_states(s))) for s in ER.regulatory_sites(m)] :
                Pop.PSite[])
end

"""A package mechanism on reaction `rxn` from a Pop.PMech's structure."""
function to_er(p::Pop.PMech, rxn)
    met(x) = x.role === :S ? ER.Substrate(x.name) : x.role === :P ? ER.Product(x.name) :
             x.role === :Inh ? ER.CompetitiveInhibitor(x.name) : error("role $(x.role)")
    function sp(i)
        s = p.species[i]
        res = isempty(s.res_added) && isempty(s.res_subtracted) ? ER.Residual() :
              ER.Residual(ER.Substrate.(s.res_added), ER.Product.(s.res_subtracted))
        ER.Species(ER.Metabolite[met(x) for x in s.bound], s.conformation, res)
    end
    step(s) = ER.Step(sp(s.from), sp(s.to), ER.Metabolite[met(x) for x in s.consumed],
                      ER.Metabolite[met(x) for x in s.released], s.re)
    groups = [ER.Step[step(p.steps[j]) for j in js] for js in p.groups]
    p.allosteric || return ER.Mechanism(rxn, groups)
    sites = [ER.RegulatorySite(ER.AllostericRegulator.(s.ligands), s.multiplicity, s.states)
             for s in p.sites]
    ER.AllostericMechanism(rxn, groups, copy(p.group_tags), p.cat_multiplicity, sites)
end

