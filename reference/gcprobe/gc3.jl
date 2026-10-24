# ABOUTME: The F6 workload with no deserialization at all: loads only the Slice module, builds
# ABOUTME: the init_terter.part3 mechanisms from source, and derives their rate laws in a loop.
include("/tmp/claude-501/-home-denis-linux--julia-dev-EnzymeRates/12f42995-f5ba-4fda-9510-dd9b87f7d069/scratchpad/slice/src/Slice.jl")
const S = Slice
include("/tmp/claude-501/-home-denis-linux--julia-dev-EnzymeRates/12f42995-f5ba-4fda-9510-dd9b87f7d069/scratchpad/phase3/gcprobe/mechs_src.jl")
const NM = parse(Int, get(ENV, "NM", "6"))
const ROUNDS = parse(Int, get(ENV, "ROUNDS", "4"))
const FN = get(ENV, "FN", "rate_law")

function main()
    ms = [S.Mechanism(a...; kw...) for (a, kw) in RAWS[1:NM]]
    println("built ", length(ms)); flush(stdout)
    for r in 1:ROUNDS, (k, m) in enumerate(ms)
        if FN == "rate_law"
            L = S.rate_law(m)
            println(r, " ", k, " ", typeof(L), " ", S.nterms(L.D))
        elseif FN == "solve"
            println(r, " ", k, " ", length(S.solve(m).deps))
        elseif FN == "bare_law"
            println(r, " ", k, " ", S.nterms(S.bare_law(m).den))
        end
        flush(stdout)
    end
    println("finished")
end
main()
