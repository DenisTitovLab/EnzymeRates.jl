# ABOUTME: Fresh-process stress closer to identify: reads the serialized population back, ships
# ABOUTME: it to Distributed workers (pmap), derives, compiles, evaluates, frees, and full-GCs.
using Distributed, Serialization
const NW = parse(Int, get(ENV, "NW", "2"))
const NEVAL = parse(Int, get(ENV, "NEVAL", "300"))
const ROUNDS = parse(Int, get(ENV, "ROUNDS", "3"))
NW > 0 && addprocs(NW; exeflags = ["--project=$(Base.active_project())", "--heap-size-hint=900M"])
@everywhere using EnzymeRates, Random
@everywhere function _testhelper_work(m, k)
    s = EnzymeRates.rate_equation_string(m)
    d = EnzymeRates._dependent_param_exprs(m)
    acc = 0.0
    if k <= $NEVAL
        em = EnzymeRates.compile_mechanism(m)
        fp = EnzymeRates.fitted_params(em)
        rng = Xoshiro(k)
        ps = NamedTuple{(fp..., :Keq, :E_total)}(Tuple(0.1 + 9.9rand(rng) for _ in 1:length(fp) + 2))
        ms = metabolites(em)
        cs = NamedTuple{Tuple(ms)}(Tuple(0.1 + 9.9rand(rng) for _ in ms))
        for _ in 1:50
            acc += rate_equation(em, cs, ps)
        end
    end
    k % 200 == 0 && GC.gc(true)
    (length(s), length(d), acc)
end
function main()
    pop = deserialize(ARGS[1])
    println("read ", length(pop), " mechanisms; workers ", nworkers()); flush(stdout)
    ref = nothing
    for r in 1:ROUNDS
        t = @elapsed (res = pmap(_testhelper_work, pop, eachindex(pop)))
        tot = (sum(x[1] for x in res), sum(x[2] for x in res), sum(x[3] for x in res))
        ref === nothing && (ref = tot)
        tot == ref || error("round $r differs: $tot vs $ref")
        println("round ", r, " ", tot, " ", round(t; digits = 1), "s"); flush(stdout)
        res = nothing
        GC.gc(true); NW > 0 && @everywhere GC.gc(true)
        println("full gc ", r, " ok"); flush(stdout)
    end
    pop = nothing
    GC.gc(true); NW > 0 && @everywhere GC.gc(true)
    println("finished")
end
main()
