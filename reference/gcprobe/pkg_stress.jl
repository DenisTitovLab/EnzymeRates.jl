# ABOUTME: Mirrors the slice's F6 pattern with EnzymeRates objects: builds a bi-bi population
# ABOUTME: (init + one expand round), serializes it, reads it back, derives laws, frees, full GC.
using EnzymeRates, Serialization
const G = "/tmp/claude-501/-home-denis-linux--julia-dev-EnzymeRates/12f42995-f5ba-4fda-9510-dd9b87f7d069/scratchpad/phase3/gcprobe"
const ROUNDS = parse(Int, get(ENV, "ROUNDS", "4"))
const MAXN = parse(Int, get(ENV, "MAXN", "400"))

const rxn = @enzyme_reaction begin
    substrates: A[CX], B[N]
    products: P[C], Q[NX]
    allosteric_regulators: R
    oligomeric_state: 2
end

function build()
    seeds = EnzymeRates.init_mechanisms(rxn)
    kids = EnzymeRates.expand_mechanisms(seeds, rxn)
    pop = Union{EnzymeRates.Mechanism, EnzymeRates.AllostericMechanism}[seeds; kids]
    println("population: ", length(seeds), " seeds + ", length(kids), " children; ",
            count(m -> m isa EnzymeRates.AllostericMechanism, pop), " allosteric"); flush(stdout)
    serialize(joinpath(G, "pkg_pop.jls"), pop)
    length(pop)
end

function derive_all(path)
    pop = deserialize(path)
    n, chars, ndep = 0, 0, 0
    for m in pop[1:min(end, MAXN)]
        p = EnzymeRates.parameters(m)
        s = EnzymeRates.rate_equation_string(m)
        d = EnzymeRates._dependent_param_exprs(m)
        n += length(p); chars += length(s); ndep += length(d)
    end
    (n, chars, ndep)
end

function main()
    build()
    GC.gc(true)
    ref = nothing
    for r in 1:ROUNDS
        t = @elapsed (x = derive_all(joinpath(G, "pkg_pop.jls")))
        ref === nothing && (ref = x)
        x == ref || error("round $r differs: $x vs $ref")
        println("round ", r, " ", x, " ", round(t; digits = 1), "s"); flush(stdout)
        GC.gc(true)
        println("full gc ", r, " ok"); flush(stdout)
    end
    println("finished")
end
main()
