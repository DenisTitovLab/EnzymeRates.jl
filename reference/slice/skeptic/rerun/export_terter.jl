# ABOUTME: Exports the plain init_mechanisms of ter_ter_rxn (test_mechanism_enumeration.jl) with
# ABOUTME: the package, structure only, as Pop.PMech records in chunks; usage: export_terter.jl
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

"""The structure of a plain package mechanism `m` as a `Pop.PMech` (no package outputs)."""
function structure(m, id, source)
    groups = ER.steps(m)
    flat = [s for g in groups for s in g]
    species = ER._forms(groups)
    idx = Dict(sp => i for (i, sp) in enumerate(species))
    rxn = ER.reaction(m)
    subs = Symbol[ER.name(x) for x in ER.substrates(rxn)]
    prods = Symbol[ER.name(x) for x in ER.products(rxn)]
    off = cumsum([0; length.(groups)])
    Pop.PMech(; id, source, substrates = subs, products = prods,
              metabolites = [subs; prods],
              species = Pop.PSpecies[_pspecies(sp) for sp in species],
              steps = Pop.PStep[Pop.PStep(idx[ER.from_species(s)], idx[ER.to_species(s)],
                                          _pmets(ER.consumed(s)), _pmets(ER.released(s)),
                                          ER.is_equilibrium(s)) for s in flat],
              groups = [collect(off[g] + 1:off[g + 1]) for g in eachindex(groups)],
              family = _family(groups, species, flat, Symbol[], ER.RegulatorySite[], subs,
                               prods))
end

t = @elapsed raw = ER.init_mechanisms(ter_ter_rxn)
ms = unique(raw)
@printf("init_mechanisms(ter_ter_rxn): %d returned, %d distinct, %.1f s; maxrss %.0f MB\n",
        length(raw), length(ms), t, Sys.maxrss() / 2^20)
raw = nothing
# The lift must be the mechanism itself, so exporting `m` exports what the derivation sees.
rng = Xoshiro(20261010)
lift_ok = count(i -> ER._concrete(ER.compile_mechanism(ms[i])) == ms[i],
                rand(rng, eachindex(ms), 200))
println("lift equals the mechanism in $lift_ok of 200 sampled")
const CHUNK = 25_000
nchunks = cld(length(ms), CHUNK)
nforms, nsteps, nre, fam = Int[], Int[], Int[], Dict{Symbol, Int}()
t0 = time()
for c in 1:nchunks
    part = Pop.PMech[]
    for i in ((c - 1) * CHUNK + 1):min(c * CHUNK, length(ms))
        p = structure(ms[i], "init_terter/$i", "init_terter[$i]")
        push!(part, p)
        push!(nforms, length(p.species)); push!(nsteps, length(p.steps))
        push!(nre, count(s -> s.re, p.steps))
        foreach(f -> fam[f] = get(fam, f, 0) + 1, p.family)
    end
    serialize(joinpath(OUT, "pop_init_terter.part$c.jls"), part)
    @printf("chunk %d: %d records, %.0f s\n", c, length(part), time() - t0)
    flush(stdout)
end
q(v) = @sprintf("median %d, max %d", median(v), maximum(v))
println("forms: ", q(nforms), "; steps: ", q(nsteps), "; RE steps: ", q(nre))
println("families: ", join(("$k=$v" for (k, v) in sort!(collect(fam); by = first)), ", "))
@printf("done in %.0f s; maxrss %.0f MB\n", time() - t0, Sys.maxrss() / 2^20)
