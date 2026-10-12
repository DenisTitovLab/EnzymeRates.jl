# ABOUTME: Crash probe: deserialized pipeline rows partly alive, then the slice's law derivation
# ABOUTME: (COMPILE=1 adds the generated rate_equation call), then a full GC.
include(joinpath(@__DIR__, "..", "indep.jl"))
function main()
    rows = [x for k in 1:11 for x in deserialize(joinpath(IDATA_IN, "pipe_terter$k.jls")) if x.ok]
    keep = get(ENV, "KEEP", "some")
    sub = keep == "all" ? rows : keep == "some" ? rows[1:1000:end] : []
    rows = nothing
    GC.gc()
    recs = Dict(c.index => c for c in deserialize(joinpath(CORE, "expand_bibi.jls")))
    rng = Xoshiro(3)
    for index in (626, 642, 658)
        m = recs[index].mech
        L = S.rate_law(m)
        if get(ENV, "COMPILE", "0") == "1"
            em = S.compile(m)
            concs = (; (x => 10.0^(4 * rand(rng) - 2) for x in S.metabolites(m))...)
            params = (; (x => (x === :E_total ? 1.0 : 10.0^(4 * rand(rng) - 2))
                         for x in S.parameters(m))...)
            S.rate_equation(em, concs, params)
        end
        println(index, " ", S.nterms(L.D))
    end
    GC.gc()
    println("full GC ok; kept ", length(sub))
end
main()
