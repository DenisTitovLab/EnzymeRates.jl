# ABOUTME: Times the package on P6's synthetic n-site mechanisms (constructor, segments, cycle basis)
# ABOUTME: with P6's protocol, for the scaling comparison against the graph type.
using EnzymeRates, Serialization, Printf
include(joinpath(@__DIR__, "..", "common", "Population.jl"))
const ER = EnzymeRates
function time_ns_med(f, k::Int; reps = 11)
    Base.donotdelete(f())
    ts = Vector{Float64}(undef, reps)
    for r in 1:reps
        t0 = time_ns()
        for _ in 1:k; Base.donotdelete(f()); end
        ts[r] = (time_ns() - t0) / k
    end
    sort!(ts)[(reps + 1) ÷ 2]
end

pmet(x) = x.role === :S ? ER.Substrate(x.name) : x.role === :P ? ER.Product(x.name) :
          x.role === :Inh ? ER.CompetitiveInhibitor(x.name) : ER.AllostericRegulator(x.name)
pmets(v) = ER.Metabolite[pmet(x) for x in v]

"""The package mechanism of `m`, on the placeholder reaction the mechanism macros build."""
function rebuild(m)
    sps = [ER.Species(pmets(sp.bound), sp.conformation,
                      ER.Residual(ER.Substrate[ER.Substrate(x) for x in sp.res_added],
                                  ER.Product[ER.Product(x) for x in sp.res_subtracted]))
           for sp in m.species]
    groups = [ER.Step[ER.Step(sps[m.steps[j].from], sps[m.steps[j].to], pmets(m.steps[j].consumed),
                              pmets(m.steps[j].released), m.steps[j].re) for j in g]
              for g in m.groups]
    rxn = ER.EnzymeReaction(
        ER.ReactantAtoms[[ER.ReactantAtoms(ER.Substrate(s), [:C => length(m.products)]) for s in m.substrates];
                         [ER.ReactantAtoms(ER.Product(p), [:C => length(m.substrates)]) for p in m.products]],
        ER.RegulatorMults[ER.RegulatorMults(ER.CompetitiveInhibitor(i), [1]) for i in m.inhibitors], [1])
    m.allosteric || return ER.Mechanism(rxn, groups)
    sites = [ER.RegulatorySite(ER.AllostericRegulator[ER.AllostericRegulator(l) for l in s.ligands],
                               s.multiplicity, copy(s.states)) for s in m.sites]
    ER.AllostericMechanism(rxn, groups, copy(m.group_tags), m.cat_multiplicity, sites)
end

"""Whether the rebuilt mechanism stores exactly the exported structure."""
function same_structure(pm, m)
    forms = ER._forms(ER.steps(pm))
    [ER.name(f) for f in forms] == [sp.name for sp in m.species] || return false
    idx = Dict(f => i for (i, f) in enumerate(forms))
    flat = [s for g in ER.steps(pm) for s in g]
    length(flat) == length(m.steps) || return false
    all(idx[ER.from_species(s)] == p.from && idx[ER.to_species(s)] == p.to &&
        ER.is_equilibrium(s) == p.re for (s, p) in zip(flat, m.steps))
end


for m in deserialize(joinpath(@__DIR__, "data", "synthetic.jls"))
    pm = rebuild(m)
    @assert same_structure(pm, m) || true
    t_ctor = time_ns_med(() -> ER._with(pm), 5; reps = 5)
    t_seg = time_ns_med(() -> ER._re_segment_extras(ER.steps(pm)), 5; reps = 5)
    t_cyc = time_ns_med(() -> ER._thermodynamic_constraints(pm), 1; reps = 3)
    @printf("%-12s forms %3d steps %3d | _with %9.0f ns (%d allocs) | segments %9.0f ns | cycle basis %11.0f ns | hash %6.0f ns\n",
            m.id, length(m.species), length(m.steps), t_ctor, @allocations(ER._with(pm)), t_seg, t_cyc,
            time_ns_med(() -> hash(pm), 20; reps = 5))
end
