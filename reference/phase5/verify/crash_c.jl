# ABOUTME: GC crash test of a law-type variant (ARGS[1] = C, U, D or E): a caller loop derives a fresh law
# ABOUTME: per mechanism and merges a Poly field with a Vector field (L.H vs L.sites), then full GCs.
const SP = "/tmp/claude-501/-home-denis-linux--julia-dev-EnzymeRates/12f42995-f5ba-4fda-9510-dd9b87f7d069/scratchpad"
const VD = joinpath(SP, "phase45/research5/verify")
const VAR = ARGS[1]
include(VAR in ("D", "E") ? joinpath(VD, "variant$VAR", "Slice$VAR.jl") :
        joinpath(SP, "phase45/research5/law/variants", VAR, "Slice$VAR.jl"))
const S = getfield(Main, Symbol("Slice", VAR))
include(joinpath(VD, "sliceside.jl"))
pop = S.Mechanism[]
for l in eachline(joinpath(VD, "data", "scan_pop.tsv"))
    try
        push!(pop, mkmech(S, eval(Meta.parse(l))))
    catch
    end
end
println("julia ", VERSION, " variant ", VAR, ": ", length(pop), " mechanisms"); flush(stdout)

function churn(ms::Vector{S.Mechanism}, full::Bool)
    t = 0
    for m in ms
        L = S.rate_law(m)
        y = L isa S.Law ? L.H : L.sites
        t += y isa S.Poly ? S.nterms(y) : length(y)
        full && GC.gc(true)
    end
    t
end
for r in 1:3
    a = churn(pop, false)
    b = churn(pop[1:60], true)
    GC.gc(true)
    println("round $r: $a $b, full gc ok"); flush(stdout)
end
