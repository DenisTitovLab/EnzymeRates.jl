# ABOUTME: Rebuilds every population's mechanisms with EnzymeRates and times its constructor, hash,
# ABOUTME: ==, analyses and inactive-state rebuild with P6's protocol; writes data/pkg_measure.jls.
using EnzymeRates, Serialization, Statistics, Printf
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

serbytes(x) = (io = IOBuffer(); serialize(io, x); position(io))
rows = NamedTuple[]
pops = ("specs", "handwritten", "init_uniuni", "init_unibi", "init_bibi", "init_bibi_pp",
        "expand_bibi", "allo_uniuni", "allo_bibi")
for pop in pops, m in Pop.load_population(pop)
    pm = rebuild(m)
    same = same_structure(pm, m)
    pm2 = ER._with(pm)
    @assert pm == pm2
    cm = m.allosteric ? ER.Mechanism(ER.reaction(pm), ER.steps(pm)) : pm
    t_ctor = time_ns_med(() -> ER._with(pm), 20)
    alloc_bytes = @allocated ER._with(pm)
    alloc_count = @allocations ER._with(pm)
    t_canon = time_ns_med(() -> ER._canonical_groups(ER.reaction(pm), ER.steps(pm)), 20)
    t_hash = time_ns_med(() -> hash(pm), 200)
    t_eq = time_ns_med(() -> pm == pm2, 200)
    t_segments = time_ns_med(() -> ER._re_segment_extras(ER.steps(pm)), 20)
    t_cycles = time_ns_med(() -> ER._thermodynamic_constraints(cm), 5)
    t_bottomless = time_ns_med(() -> ER._bottomless_re_segment(ER.reaction(pm), ER.steps(pm)), 20)
    t_inactive = NaN; t_state_parts = NaN; t_inactive_cycles = NaN
    if m.allosteric
        t_inactive = time_ns_med(() -> ER._inactive_groups(pm), 20)
        t_state_parts = time_ns_med(() -> ER._state_parts(pm, :I), 10)
        icm = first(ER._state_parts(pm, :I))
        t_inactive_cycles = time_ns_med(() -> (ER._re_segment_extras(ER.steps(icm));
                                               ER._thermodynamic_constraints(icm)), 5)
    end
    species, segments, extras = ER._re_segment_extras(ER.steps(pm))
    C, haldane = ER._thermodynamic_constraints(cm)
    push!(rows, (; id = m.id, allosteric = m.allosteric, same, t_ctor, alloc_bytes, alloc_count,
                 t_canon, t_hash, t_eq, t_segments, t_cycles, t_bottomless, t_inactive,
                 t_state_parts, t_inactive_cycles, size_fresh = Base.summarysize(pm2),
                 ser_fresh = serbytes(pm2),
                 segments = [sort([String(ER.name(species[i])) for i in s]) for s in segments],
                 extras = Dict(String(ER.name(species[i])) =>
                               Dict(String(k) => v for (k, v) in extras[i]) for i in eachindex(species)),
                 cycles = Matrix{Int}(C), haldane = collect(Int, haldane)))
end
serialize(joinpath(@__DIR__, "data", "pkg_measure.jls"), rows)

q(v) = (v = filter(!isnan, Float64.(v)); isempty(v) ? (NaN, NaN, NaN) : (median(v), quantile(v, 0.9), maximum(v)))
@printf("rebuilt mechanisms storing exactly the exported structure: %d of %d\n",
        count(r -> r.same, rows), length(rows))
for (label, sel) in (("plain", r -> !r.allosteric), ("allosteric", r -> r.allosteric), ("all", r -> true))
    rs = filter(sel, rows)
    println("== $label ($(length(rs)) mechanisms) — median / p90 / max")
    for f in (:t_ctor, :alloc_bytes, :alloc_count, :t_canon, :t_hash, :t_eq, :t_segments,
              :t_cycles, :t_bottomless, :t_inactive, :t_state_parts, :t_inactive_cycles,
              :size_fresh, :ser_fresh)
        @printf("  %-17s %10.1f %10.1f %10.1f\n", f, q(getfield.(rs, f))...)
    end
end
