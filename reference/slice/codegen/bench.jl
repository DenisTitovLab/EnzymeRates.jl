# ABOUTME: Codegen checks (a) and (b): per law the generator time, first-call compile time and
# ABOUTME: native code, allocations and ns per call (T1, T1n, T2). Usage: bench.jl label sets…
include(joinpath(@__DIR__, "lib.jl"))

const LABEL, SETS = ARGS[1], ARGS[2:end]
"""A call the optimizer cannot hoist out of a loop: the per-call cost as executed."""
@noinline barrier(em, c, p) = S.rate_equation(em, c, p)

"""Per-row time of a loop over 1000 distinct concentration rows at fixed parameters (the
fitter's shape), minimum over 5 passes."""
function t2_of(em, rows, p)
    best, acc = Inf, 0.0
    for _ in 1:5
        acc = 0.0
        t = @elapsed for c in rows
            acc += S.rate_equation(em, c, p)
        end
        best = min(best, t / length(rows))
    end
    isfinite(acc) || error("non-finite result")
    best
end

function run_one(r, k; record = true)
    m = r.c.mech
    rng = Xoshiro(1000 + k)
    t_gen = minimum(@elapsed(S.body(S.program(m))) for _ in 1:3)
    b = S.body(S.program(m))
    em = S.compile(m)
    c, p = rand_concs(m, rng), rand_params(m, rng)
    record || return S.rate_equation(em, map(BigFloat, c), map(BigFloat, p))
    s0 = snapshot()
    t_first = @elapsed S.rate_equation(em, c, p)
    s1 = snapshot()
    a = allocs_of(S.rate_equation, em, c, p)
    t1 = t1_of(S.rate_equation, em, c, p)
    t1n = t1_of(barrier, em, c, p)
    t2 = t2_of(em, [rand_concs(m, rng) for _ in 1:1000], p)
    ci = code_typed(S.rate_equation, Tuple{typeof(em), typeof(c), typeof(p)})[1].first
    (set = r.set, id = r.p.id, t_gen, t_first, compile_s = (s1.compile - s0.compile) / 1e9,
     jit = s1.jit - s0.jit, allocs = a, t1_ns = 1e9 * t1, t1n_ns = 1e9 * t1n,
     t2_ns = 1e9 * t2, ops = countops(b), nstmts = length(ci.code),
     locals = length(S.program(m).defs), sig = length(S.sig(m)), allo = S.is_allosteric(m),
     n = m.multiplicity)
end

# Warm the generator path on BigFloat calls, so each measured Float64 call compiles only
# its own method instance.
foreach(((k, r),) -> run_one(r, -k; record = false), enumerate(load_records(("init_unibi",))))
rows = NamedTuple[]
for set in SETS, (k, r) in enumerate(load_records((set,)))
    x = try
        run_one(r, k)
    catch err
        println("$set $(r.p.id): ", sprint(showerror, err))
        continue
    end
    push!(rows, x)
    @printf("%-14s %-18s gen %6.2f ms first %6.3f s comp %6.3f s jit %6d allocs %d t1 %5.1f \
            t1n %5.1f t2 %5.1f ns ops %5d stmts %5d\n", set, x.id, 1e3 * x.t_gen, x.t_first,
            x.compile_s, x.jit, x.allocs, x.t1_ns, x.t1n_ns, x.t2_ns, x.ops, x.nstmts)
    flush(stdout)
end
serialize(joinpath(@__DIR__, "data", "bench_$LABEL.jls"), rows)
println("done: ", length(rows), " laws; maxrss ", round(Sys.maxrss() / 2^20), " MB")
