# ABOUTME: The package's rate_equation on the slowest p1_enum_inh laws, rebuilt by the skeptic's
# ABOUTME: enumeration (bi-bi with a dead-end inhibitor, init + two expand rounds): T1/T1n/T2.
using EnzymeRates, Random, Serialization, Printf
const ER = EnzymeRates
include(joinpath(@__DIR__, "..", "..", "proto", "common", "Population.jl"))
const POP = deserialize(joinpath(@__DIR__, "..", "..", "proto", "skeptic_P1-engine", "data",
                                 "enum_inh.jls"))
rxn = @enzyme_reaction(begin
    substrates: A[C], B[N]
    products: P[C], Q[N]
    dead_end_inhibitors: I
end)
seeds = unique(ER.init_mechanisms(rxn))
level = filter(m -> m isa ER.Mechanism, seeds)
seen = Set{Any}(level)
all_m = copy(level)
for round in 1:2
    kids = Any[]
    for p in level
        try
            append!(kids, ER.expand_mechanisms([p], rxn))
        catch
        end
    end
    new = [k for k in unique(kids) if k isa ER.Mechanism && !(k in seen)]
    union!(seen, new)
    global level = length(new) > 400 ?
        new[sort!(randperm(Random.Xoshiro(round), length(new))[1:400])] : new
    append!(all_m, level)
end
println("rebuilt ", length(all_m), " bi_bi_with_inh mechanisms")

@noinline barrier(em, c, p) = ER.rate_equation(em, c, p)
function tloop(f::F, em, c, p) where {F}
    f(em, c, p)
    best = Inf
    for _ in 1:5
        acc = 0.0
        t = @elapsed for _ in 1:10_000
            acc += f(em, c, p)
        end
        best = min(best, t / 10_000)
        isfinite(acc) || error("non-finite")
    end
    1e9 * best
end
Base.cumulative_compile_timing(true)
for i in (1013, 959, 584, 977, 1015)
    id = "bi_bi_with_inh/$i"
    p = only(filter(x -> x.id == id, POP))
    m = all_m[i]
    em = ER.compile_mechanism(m)
    names = [ER.name(sp) for sp in ER._forms(ER.steps(ER._concrete(em)))]
    names == [sp.name for sp in p.species] || error("$id: rebuilt mechanism differs")
    rng = Xoshiro(i)
    c = NamedTuple{Tuple(ER.metabolites(em))}(Tuple(10.0^(4 * rand(rng) - 2)
                                                    for _ in ER.metabolites(em)))
    pr = NamedTuple{Tuple(ER.parameters(em))}(Tuple(x === :E_total ? 1.0 : 10.0^(4 * rand(rng) - 2)
                                                    for x in ER.parameters(em)))
    c0 = Base.cumulative_compile_time_ns()[1]
    j0 = Base.jit_total_bytes()
    t_first = @elapsed ER.rate_equation(em, c, pr)
    comp = (Base.cumulative_compile_time_ns()[1] - c0) / 1e9
    a = (ER.rate_equation(em, c, pr); @allocated ER.rate_equation(em, c, pr))
    @printf("%s package: first %.3f s comp %.3f s jit %d allocs %d T1 %.1f T1n %.1f ns\n", id,
            t_first, comp, Base.jit_total_bytes() - j0, a, tloop(ER.rate_equation, em, c, pr),
            tloop(barrier, em, c, pr))
end
