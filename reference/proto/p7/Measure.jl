# ABOUTME: Shared measurement harness of the P7 probe: compile/JIT snapshots, the package's
# ABOUTME: performance-test protocol (T1), a 1000-row loop (T2), and sample NamedTuples.
using Random, Printf
Base.cumulative_compile_timing(true)

"""GC, then compile time (ns), JIT code bytes and pool live bytes."""
function snapshot()
    GC.gc(true); GC.gc(true)
    (compile = Int(Base.cumulative_compile_time_ns()[1]), jit = Int(Base.jit_total_bytes()),
     pool = Int(ccall(:jl_gc_pool_live_bytes, Int64, ())))
end

nt(ps) = NamedTuple{Tuple(first.(ps))}(Tuple(Float64.(last.(ps))))

"""Bytes one call allocates (after a warm-up call)."""
function allocs_of(f, m, c, p)
    f(m, c, p)
    @allocated f(m, c, p)
end

"""The package's `_testhelper_test_rate_equation_performance` protocol: minimum over 5 batches
of 10_000 calls of the per-call time."""
function t1_of(f, m, c, p)
    f(m, c, p)
    best = Inf
    acc = 0.0
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

"""Per-row time of a loop over `rows` (distinct concentrations, fixed parameters), the shape of
the fitter's loss: minimum over 5 passes."""
function t2_of(f, m, rows, p)
    best = Inf
    acc = 0.0
    for _ in 1:5
        acc = 0.0
        t = @elapsed for c in rows
            acc += f(m, c, p)
        end
        best = min(best, t / length(rows))
    end
    isfinite(acc) || error("non-finite result")
    best
end

"""1000 random concentration NamedTuples (log-uniform in [1e-2, 1e2]) for metabolites `names`."""
rows_for(names, seed) = (rng = Xoshiro(seed);
    [NamedTuple{Tuple(names)}(Tuple(10.0^(4 * rand(rng) - 2) for _ in names)) for _ in 1:1000])
