# ABOUTME: Enumerates plain mechanisms of reactions with a dead-end inhibitor (init + two expand rounds) and
# ABOUTME: exports their structure only (no derivation) as Pop.PMech, to test the speed-free predicate on them.
using EnzymeRates, Random, Serialization, Printf
const C = joinpath(@__DIR__, "..", "common")
include(joinpath(C, "Population.jl"))
const ER = EnzymeRates

const rxns = [
    "uni_uni_with_reg" => @enzyme_reaction(begin
        substrates: S[C]
        products: P[C]
        dead_end_inhibitors: I
    end),
    "bi_bi_with_inh" => @enzyme_reaction(begin
        substrates: A[C], B[N]
        products: P[C], Q[N]
        dead_end_inhibitors: I
    end),
]

_role(::ER.Substrate) = :S
_role(::ER.Product) = :P
_role(::ER.CompetitiveInhibitor) = :Inh
_role(::ER.AllostericRegulator) = :Reg
_pmets(xs) = Pop.PMet[Pop.PMet(ER.name(x), _role(x)) for x in xs]

"""Structure of a plain mechanism as a Pop.PMech (no package outputs)."""
function to_struct(m, id, src)
    em = ER.compile_mechanism(m)
    lm = ER._concrete(em)
    groups = ER.steps(lm)
    flat = [s for g in groups for s in g]
    species = ER._forms(groups)
    idx = Dict(sp => i for (i, sp) in enumerate(species))
    rxn = ER.reaction(lm)
    gof = [g for (g, grp) in enumerate(groups) for _ in grp]
    Pop.PMech(; id, source = src,
        substrates = Symbol[ER.name(x) for x in ER.substrates(rxn)],
        products = Symbol[ER.name(x) for x in ER.products(rxn)],
        inhibitors = Symbol[ER.name(ER.regulator(r)) for r in ER.regulators(rxn)
                            if ER.regulator(r) isa ER.CompetitiveInhibitor],
        metabolites = collect(Symbol, ER.metabolites(em)),
        species = [Pop.PSpecies(_pmets(ER.bound(sp)), ER.conformation(sp),
                                Symbol[ER.name(a) for a in ER.added(ER.residual(sp))],
                                Symbol[ER.name(p) for p in ER.subtracted(ER.residual(sp))],
                                ER.name(sp)) for sp in species],
        steps = [Pop.PStep(idx[ER.from_species(s)], idx[ER.to_species(s)], _pmets(ER.consumed(s)),
                           _pmets(ER.released(s)), ER.is_equilibrium(s)) for s in flat],
        groups = [findall(==(g), gof) for g in eachindex(groups)],
        step_names = [Symbol[ER.name(p, lm) for p in ER._step_constants(s, :None)] for s in flat])
end

out = Pop.PMech[]
for (label, rxn) in rxns
    t0 = time()
    seeds = unique(ER.init_mechanisms(rxn))
    level = filter(m -> m isa ER.Mechanism, seeds)
    seen = Set{Any}(level)
    all_m = copy(level)
    for round in 1:2
        kids = Any[]
        for p in level
            try
                append!(kids, ER.expand_mechanisms([p], rxn))
            catch e
                println("expand failed: ", first(sprint(showerror, e), 200))
            end
        end
        new = [k for k in unique(kids) if k isa ER.Mechanism && !(k in seen)]
        union!(seen, new)
        cap = 400
        level = length(new) > cap ? new[sort!(randperm(Random.Xoshiro(round), length(new))[1:cap])] : new
        append!(all_m, level)
        @printf("%s round %d: %d new plain children (kept %d)\n", label, round, length(new), length(level))
    end
    for (i, m) in enumerate(all_m)
        push!(out, to_struct(m, "$label/$i", "$label[$i]"))
    end
    @printf("%s: %d plain mechanisms in %.1f s\n", label, length(all_m), time() - t0)
end
mkpath(joinpath(@__DIR__, "data"))
serialize(joinpath(@__DIR__, "data", "enum_inh.jls"), out)
println("wrote ", length(out))
