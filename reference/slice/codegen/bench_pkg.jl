# ABOUTME: The package's rate_equation on the shared populations, rebuilt as the exporter built
# ABOUTME: them: compile, native code, allocations, T1/T1n/T2. Usage: bench_pkg.jl label pops...
using EnzymeRates, Random, Serialization, Printf, InteractiveUtils
const EXPORT = "/tmp/claude-501/-home-denis-linux--julia-dev-EnzymeRates/12f42995-f5ba-4fda-9510-dd9b87f7d069/scratchpad/proto/common/export/export_pop.jl"
const LABEL, POPS = ARGS[1], ARGS[2:end]
empty!(ARGS)
push!(ARGS, "specs")                       # the exporter's preamble reads its population here
for ex in Meta.parseall(read(EXPORT, String); filename = EXPORT).args
    ex isa Expr && ex.head === :const && ex.args[1] isa Expr &&
        ex.args[1].args[1] === :MEMBERS && break
    Core.eval(Main, ex)                    # the exporter's definitions up to its export run
end

Base.cumulative_compile_timing(true)
snapshot() = (GC.gc(true); (compile = Int(Base.cumulative_compile_time_ns()[1]),
                            jit = Int(Base.jit_total_bytes())))
nt(ps) = NamedTuple{Tuple(first.(ps))}(Tuple(Float64.(last.(ps))))
@noinline barrier(em, c, p) = ER.rate_equation(em, c, p)
function tloop(f::F, em, c, p) where {F}
    f(em, c, p)
    best, acc = Inf, 0.0
    for _ in 1:5
        acc = 0.0
        t = @elapsed for _ in 1:10_000
            acc += f(em, c, p)
        end
        best = min(best, t / 10_000)
    end
    isfinite(acc) || error("non-finite result")
    1e9 * best
end
function t2(em, rows, p)
    best, acc = Inf, 0.0
    for _ in 1:5
        acc = 0.0
        t = @elapsed for c in rows
            acc += ER.rate_equation(em, c, p)
        end
        best = min(best, t / length(rows))
    end
    1e9 * best
end
countops(b) = (n = Ref(0); w(x) = (x isa Expr && (x.head === :call && (n[] += 1);
                                                   foreach(w, x.args))); w(b); n[])

function run_one(x, p, k; record = true)
    em = ER.compile_mechanism(x)
    record && ER.rate_equation_string(em) != p.rate_equation_string &&
        error("$(p.id): rebuilt mechanism differs from the exported one")
    s = p.samples[1]
    c, pr = nt(s.concs), nt(s.params)
    record || return ER.rate_equation(em, map(BigFloat, c), map(BigFloat, pr))
    s0 = snapshot()
    t_first = @elapsed ER.rate_equation(em, c, pr)
    s1 = snapshot()
    a = (ER.rate_equation(em, c, pr); @allocated ER.rate_equation(em, c, pr))
    rng = Xoshiro(1000 + k)
    rows = [NamedTuple{keys(c)}(Tuple(10.0^(4 * rand(rng) - 2) for _ in c)) for _ in 1:1000]
    M = typeof(em)
    ops = try
        dep, indep = ER._dependent_param_exprs(M)
        countops(ER._rate_body(M, (indep..., :Keq, :E_total), dep))
    catch
        -1
    end
    ci = code_typed(ER.rate_equation, Tuple{M, typeof(c), typeof(pr), typeof(ER.Reduced)})
    (id = p.id, t_first, compile_s = (s1.compile - s0.compile) / 1e9, jit = s1.jit - s0.jit,
     allocs = a, t1_ns = tloop(ER.rate_equation, em, c, pr), t1n_ns = tloop(barrier, em, c, pr),
     t2_ns = t2(em, rows, pr), ops, nstmts = length(ci[1].first.code), allo = p.allosteric,
     exported_ns = p.eval_ns)
end

# Warm the generator path on BigFloat calls, so each measured Float64 call compiles only
# its own method instance.
foreach(((k, (x, p)),) -> run_one(first(x), p, -k; record = false),
        enumerate(zip(members("init_unibi"), Pop.load_population("init_unibi"))))
rows = NamedTuple[]
for pop in POPS
    ms, recs = members(pop), Pop.load_population(pop)
    length(ms) == length(recs) || error("$pop: rebuilt $(length(ms)), exported $(length(recs))")
    for (k, ((x, _), p)) in enumerate(zip(ms, recs))
        p.ok || continue
        r = try
            run_one(x, p, k)
        catch err
            println("$(p.id): ", first(sprint(showerror, err), 300))
            continue
        end
        push!(rows, r)
        @printf("%-18s first %6.3f s comp %6.3f s jit %6d allocs %d t1 %6.1f t1n %6.1f t2 \
                %6.1f ns ops %5d stmts %5d (exported %6.1f)\n", r.id, r.t_first, r.compile_s,
                r.jit, r.allocs, r.t1_ns, r.t1n_ns, r.t2_ns, r.ops, r.nstmts, r.exported_ns)
        flush(stdout)
    end
end
serialize(joinpath(@__DIR__, "data", "bench_pkg_$LABEL.jls"), rows)
println("done: ", length(rows), " laws; maxrss ", round(Sys.maxrss() / 2^20), " MB")
