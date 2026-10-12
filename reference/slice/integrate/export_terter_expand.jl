# ABOUTME: Exports one round of expand_mechanisms children of 60 sampled ter-ter init mechanisms
# ABOUTME: with the package, structure only, as Pop.PMech records; usage: export_terter_expand.jl
using EnzymeRates, Random, Serialization, Printf, Statistics
const ER = EnzymeRates
const SCRATCH = "/tmp/claude-501/-home-denis-linux--julia-dev-EnzymeRates/12f42995-f5ba-4fda-9510-dd9b87f7d069/scratchpad"
include(joinpath(SCRATCH, "proto", "common", "Population.jl"))
include(joinpath(SCRATCH, "proto", "common", "export", "export_lib.jl"))
const OUT = joinpath(@__DIR__, "data")
const TESTDIR = joinpath(dirname(dirname(pathof(EnzymeRates))), "test")

"""Top-level `const` statements of a test file whose names are in `names`."""
function test_consts(file, names)
    ast = Meta.parseall(read(joinpath(TESTDIR, file), String))
    [ex for ex in ast.args if ex isa Expr && ex.head === :const &&
     ex.args[1] isa Expr && ex.args[1].args[1] in names]
end
foreach(ex -> Core.eval(Main, ex), test_consts("test_mechanism_enumeration.jl", (:ter_ter_rxn,)))

"""The structure of a package mechanism `m`, plain or allosteric, as a `Pop.PMech` (no
package outputs)."""
function structure(m, id, source)
    allo = m isa ER.AllostericMechanism
    groups = ER.steps(m)
    flat = [s for g in groups for s in g]
    species = ER._forms(groups)
    idx = Dict(sp => i for (i, sp) in enumerate(species))
    rxn = ER.reaction(m)
    subs = Symbol[ER.name(x) for x in ER.substrates(rxn)]
    prods = Symbol[ER.name(x) for x in ER.products(rxn)]
    off = cumsum([0; length.(groups)])
    tags = allo ? collect(Symbol, ER.cat_allo_states(m)) : Symbol[]
    sites = allo ? ER.regulatory_sites(m) : ER.RegulatorySite[]
    Pop.PMech(; id, source, substrates = subs, products = prods,
              metabolites = [subs; prods], allosteric = allo, group_tags = tags,
              cat_multiplicity = allo ? ER.catalytic_multiplicity(m) : 1,
              sites = Pop.PSite[Pop.PSite(Symbol[ER.name(l) for l in ER.ligands(s)],
                                          ER.multiplicity(s), collect(Symbol, ER.allo_states(s)))
                                for s in sites],
              species = Pop.PSpecies[_pspecies(sp) for sp in species],
              steps = Pop.PStep[Pop.PStep(idx[ER.from_species(s)], idx[ER.to_species(s)],
                                          _pmets(ER.consumed(s)), _pmets(ER.released(s)),
                                          ER.is_equilibrium(s)) for s in flat],
              groups = [collect(off[g] + 1:off[g + 1]) for g in eachindex(groups)],
              family = _family(groups, species, flat, tags, sites, subs, prods))
end

init = unique(ER.init_mechanisms(ter_ter_rxn))
pick = sort!(randperm(Random.Xoshiro(20261009), length(init))[1:60])
children, times, counts = Any[], Float64[], Int[]
for (k, i) in enumerate(pick)
    cs, err = [], ""
    t = @elapsed try
        cs = ER.expand_mechanisms([init[i]], ter_ter_rxn)
    catch e
        err = first(sprint(showerror, e), 300)
    end
    append!(children, cs); push!(times, t); push!(counts, length(cs))
    @printf("parent %2d init_terter/%d: %4d children in %.2f s%s\n", k, i, length(cs), t,
            isempty(err) ? "" : "; expand_mechanisms failed: " * err)
    flush(stdout)
end
kept = unique(children)
@printf("parents %s\n%d children, %d distinct, %d allosteric; expansion %.1f s (max %.1f s)\n",
        pick, length(children), length(kept), count(c -> c isa ER.AllostericMechanism, kept),
        sum(times), maximum(times))
recs = [structure(c, "expand_terter/$j", "expand_terter[$j]") for (j, c) in enumerate(kept)]
serialize(joinpath(OUT, "pop_expand_terter.jls"), recs)
@printf("wrote %d records; maxrss %.0f MB\n", length(recs), Sys.maxrss() / 2^20)
