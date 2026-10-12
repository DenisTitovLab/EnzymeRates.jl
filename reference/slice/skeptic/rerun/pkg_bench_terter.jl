# ABOUTME: Times today's package on the ter-ter laws the slice found slowest (and on a sample of
# ABOUTME: others), with the same barrier protocol; usage: pkg_bench_terter.jl (er_oracle project)
using EnzymeRates, Random, Serialization, Printf, Statistics
const ER = EnzymeRates
const SCRATCH = "/tmp/claude-501/-home-denis-linux--julia-dev-EnzymeRates/12f42995-f5ba-4fda-9510-dd9b87f7d069/scratchpad"
const IDATA = joinpath(@__DIR__, "data")
const TESTDIR = joinpath(dirname(dirname(pathof(EnzymeRates))), "test")
ast = Meta.parseall(read(joinpath(TESTDIR, "test_mechanism_enumeration.jl"), String))
foreach(ex -> Core.eval(Main, ex), [ex for ex in ast.args if ex isa Expr && ex.head === :const &&
                                    ex.args[1] isa Expr && ex.args[1].args[1] === :ter_ter_rxn])

rows = [x for f in readdir(IDATA) if startswith(f, "bench_terter")
        for x in deserialize(joinpath(IDATA, f))]
slow = first(sort(rows; by = x -> -x.t1n_ns), 25)
pick = [slow; rows[randperm(Xoshiro(3), length(rows))[1:25]]]
init = unique(ER.init_mechanisms(ter_ter_rxn))
@noinline barrier(em, c, p) = ER.rate_equation(em, c, p)
function t1n(em, c, p)
    barrier(em, c, p)
    best = Inf
    for _ in 1:5
        acc = 0.0
        t = @elapsed for _ in 1:10_000
            acc += barrier(em, c, p)
        end
        best = min(best, t / 10_000)
    end
    best
end
rng = Xoshiro(7)
out = []
for x in pick
    i = parse(Int, split(x.id, "/")[2])
    em = ER.compile_mechanism(init[i])
    c = NamedTuple{Tuple(ER.metabolites(em))}(Tuple(10.0^(4rand(rng) - 2) for _ in ER.metabolites(em)))
    ps = (ER.fitted_params(em)..., :Keq, :E_total)
    p = NamedTuple{ps}(Tuple(n === :E_total ? 1.0 : 10.0^(4rand(rng) - 2) for n in ps))
    tc = @elapsed ER.rate_equation(em, c, p)
    a = (ER.rate_equation(em, c, p); @allocated ER.rate_equation(em, c, p))
    push!(out, (id = x.id, slice_ns = x.t1n_ns, pkg_ns = 1e9 * t1n(em, c, p), pkg_allocs = a,
                pkg_first_s = tc, slice_compile_s = x.compile, nD = x.nD))
    @printf("%-22s D %5d: slice %6.1f ns, package %7.1f ns (%d B); first call slice %.2f s, package %.2f s\n",
            x.id, x.nD, x.t1n_ns, last(out).pkg_ns, a, x.compile, tc)
    flush(stdout)
end
serialize(joinpath(IDATA, "pkg_bench_terter.jls"), out)
s, r = out[1:25], out[26:end]
@printf("slowest 25 by the slice: slice median %.1f ns, package median %.1f ns; over 120 ns: slice %d, package %d\n",
        median(x.slice_ns for x in s), median(x.pkg_ns for x in s), count(x -> x.slice_ns >= 120, s),
        count(x -> x.pkg_ns >= 120, s))
@printf("25 random compiled laws: slice median %.1f ns, package median %.1f ns; ratio median %.2f\n",
        median(x.slice_ns for x in r), median(x.pkg_ns for x in r),
        median(x.slice_ns / x.pkg_ns for x in r))
