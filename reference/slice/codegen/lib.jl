# ABOUTME: Helpers shared by the codegen-stage check scripts: records, random inputs as
# ABOUTME: NamedTuples, the BigFloat law value, and the package's timing protocol (T1).
include(joinpath(@__DIR__, "..", "engine", "lib.jl"))
using InteractiveUtils

const SHARED = ("specs", "handwritten", "init_uniuni", "init_unibi", "init_bibi",
                "init_bibi_pp", "expand_bibi")

"""Log-uniform Float64 in [1e-2, 1e2]."""
logu64(rng) = 10.0^(4 * rand(rng) - 2)

"""Random Float64 parameters for `parameters(m)` (E_total = 1)."""
rand_params(m, rng) = (; (n => (n === :E_total ? 1.0 : logu64(rng)) for n in S.parameters(m))...)

"""Random Float64 concentrations for `metabolites(m)`."""
rand_concs(m, rng) = (; (n => logu64(rng) for n in S.metabolites(m))...)

"""`v` of law `L` at Float64 NamedTuples, in BigFloat (`E_total` included)."""
law_big(L, concs, params) =
    BigFloat(params.E_total) * lawval(L, Dict(pairs(concs)), map(BigFloat, params))

"""The package's performance protocol: minimum over 5 batches of 10_000 calls (s/call)."""
function t1_of(f, m, c, p)
    f(m, c, p)
    best, acc = Inf, 0.0
    for _ in 1:5
        acc = 0.0
        t = @elapsed for _ in 1:10_000
            acc += f(m, c, p)
        end
        best = min(best, t / 10_000)
    end
    isfinite(acc) || error("non-finite result")
    best
end

"""Bytes one call allocates (after a warm-up call)."""
allocs_of(f, m, c, p) = (f(m, c, p); @allocated f(m, c, p))

Base.cumulative_compile_timing(true)
"""GC, then cumulative compile time (ns) and JIT code bytes."""
function snapshot()
    GC.gc(true)
    (compile = Int(Base.cumulative_compile_time_ns()[1]), jit = Int(Base.jit_total_bytes()))
end

"""Number of `:call` nodes in an expression (operation count of a generated body)."""
countops(b) = (n = Ref(0); w(x) = (x isa Expr && (x.head === :call && (n[] += 1);
                                                   foreach(w, x.args))); w(b); n[])
