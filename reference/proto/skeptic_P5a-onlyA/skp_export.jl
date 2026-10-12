# ABOUTME: Skeptic P5a: exports allosteric :OnlyA populations the probe did not cover (ping-pong bi-bi,
# ABOUTME: uni-bi and ter-bi at oligomeric state 2, uni-uni with a competitive inhibitor and a regulator).
using EnzymeRates, Random, Serialization, Printf, Statistics
include(joinpath(@__DIR__, "..", "common", "Population.jl"))
include(joinpath(@__DIR__, "..", "common", "export", "export_lib.jl"))

const POP = ARGS[1]
const CAP = parse(Int, get(ARGS, 2, "400"))
subsample(xs, n, seed) = length(xs) <= n ? xs : xs[sort!(randperm(Random.Xoshiro(seed), length(xs))[1:n])]
function expand_each(parents, rxn)
    out = Union{ER.Mechanism, ER.AllostericMechanism}[]
    for p in parents
        try append!(out, ER.expand_mechanisms([p], rxn)) catch e; println("expand failed: ", _errtext(e)) end
    end
    out
end
rxn = if POP == "allo_bibi_pp"
    @enzyme_reaction begin
        substrates: A[CX], B[N]
        products: P[C], Q[NX]
        oligomeric_state: 2
    end
elseif POP == "allo_unibi"
    @enzyme_reaction begin
        substrates: S[AB]
        products: P[A], Q[B]
        oligomeric_state: 2
    end
elseif POP == "allo_uniuni_inh"
    @enzyme_reaction begin
        substrates: S[C]
        products: P[C]
        allosteric_regulators: R2
        competitive_inhibitors: R1
        oligomeric_state: 2
    end
elseif POP == "allo_terbi"
    @enzyme_reaction begin
        substrates: A[C], B[N], D[X]
        products: P[CN], Q[X]
        oligomeric_state: 2
    end
end
t0 = time()
seeds = unique(ER.init_mechanisms(rxn))
POP == "allo_terbi" && (seeds = subsample(seeds, 40, 7))
l1 = unique([c for m in seeds for c in ER._expand_to_allosteric(m, rxn)])
pool = filter(c -> c isa ER.AllostericMechanism && :OnlyA in ER.cat_allo_states(c), l1)
if POP in ("allo_uniuni_inh", "allo_unibi")
    l2 = unique(filter(c -> c isa ER.AllostericMechanism && !(c in l1), expand_each(l1, rxn)))
    append!(pool, filter(c -> :OnlyA in ER.cat_allo_states(c), l2))
    pool = unique(pool)
    POP == "allo_uniuni_inh" && begin
        l3 = unique(filter(c -> c isa ER.AllostericMechanism && !(c in l2) && !(c in l1),
                           expand_each(subsample(l2, 60, 9), rxn)))
        append!(pool, filter(c -> :OnlyA in ER.cat_allo_states(c), l3)); pool = unique(pool)
    end
end
println("$POP: $(length(seeds)) seeds, $(length(l1)) to_allosteric children, $(length(pool)) :OnlyA mechanisms in the pool (", round(time() - t0; digits = 1), " s)")
pool = subsample(pool, CAP, 11)
out = Pop.PMech[]
for (i, m) in enumerate(pool)
    pm = try
        to_pmech(m, "$POP/$i", "skp:$POP"; seed = 900_000 + i, do_bench = false, do_kcat = false)
    catch e
        println("FATAL $POP/$i: ", _errtext(e)); continue
    end
    push!(out, pm)
end
fam = Dict{Symbol, Int}()
for pm in out, f in pm.family; fam[f] = get(fam, f, 0) + 1; end
println("$POP: exported $(length(out)); not ok $(count(pm -> !pm.ok, out)); families ",
        join(("$k=$v" for (k, v) in sort!(collect(fam); by = first)), ", "))
for pm in out, (k, v) in pm.stage_errors
    println("  stage error ", pm.id, " ", k, ": ", first(replace(v, '\n' => ' '), 150))
end
serialize(joinpath(@__DIR__, "data", "$POP.jls"), out)
println("done in ", round(time() - t0; digits = 1), " s; maxrss ", round(Sys.maxrss() / 2^20), " MB")
