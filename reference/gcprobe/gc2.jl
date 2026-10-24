# ABOUTME: Builds the slice mechanisms of a few init_terter.part3 records, optionally drops the
# ABOUTME: deserialized population first, then derives their rate laws in a loop.
include("/tmp/claude-501/-home-denis-linux--julia-dev-EnzymeRates/12f42995-f5ba-4fda-9510-dd9b87f7d069/scratchpad/slice/skeptic/scripts/indep.jl")
const NM = parse(Int, get(ENV, "NM", "6"))
const ROUNDS = parse(Int, get(ENV, "ROUNDS", "4"))
const KEEP = get(ENV, "KEEP", "0") == "1"
const TODO = "/tmp/claude-501/-home-denis-linux--julia-dev-EnzymeRates/12f42995-f5ba-4fda-9510-dd9b87f7d069/scratchpad/phase3/gcprobe/todo.tsv"

function main()
    idx = [parse(Int, split(l, '\t')[3]) for l in readlines(TODO) if split(l, '\t')[2] == "init_terter.part3"][1:NM]
    ps = deserialize(joinpath(IDATA_IN, "pop_init_terter.part3.jls"))
    ms = [S.Converted(ps[i], i).mech for i in idx]
    KEEP || (ps = nothing)
    GC.gc()
    println("built ", length(ms)); flush(stdout)
    for r in 1:ROUNDS, (k, m) in enumerate(ms)
        L = S.rate_law(m)
        println(r, " ", k, " ", typeof(L), " ", L isa S.Law ? S.nterms(L.D) : S.nterms(L.active.D))
        flush(stdout)
    end
    KEEP && println(length(ps))
    println("finished")
end
main()
