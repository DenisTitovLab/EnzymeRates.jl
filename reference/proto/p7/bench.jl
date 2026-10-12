# ABOUTME: P7 benchmark of the generated evaluator: per mechanism compile time, JIT bytes, allocations,
# ABOUTME: T1/T2 timing and op counts. Usage: bench.jl <flat|flatr|fact|facth|flath> <unique|shared>
include(joinpath(@__DIR__, "Encode.jl"))
include(joinpath(@__DIR__, "Select.jl"))
include(joinpath(@__DIR__, "Measure.jl"))
using Serialization, InteractiveUtils
const VARIANT, MODE = ARGS[1], ARGS[2]
const OPT = Dict("flat" => (false, false), "flatr" => (true, false), "fact" => (true, false),
                 "facth" => (true, true), "flath" => (true, true), "halh" => (true, true),
                 "hal" => (true, false), "segh" => (true, true))[VARIANT]
mklaw(m) = VARIANT in ("flat", "flatr", "flath") ? P7.flat_law(m) :
           P7.factored_law(m; expand_den = VARIANT in ("hal", "halh"), blocks = VARIANT != "segh")
countops(b) = (n = Ref(0); w(x) = (x isa Expr && (x.head === :call && (n[] += 1); foreach(w, x.args))); w(b); n[])

function run_one(m, k; record = true)
    t_law = @elapsed law = mklaw(m)
    code_ = code(law, (OPT[1], OPT[2], MODE == "unique" ? k : 0))
    L = Gen.RateLaw{Tuple(code_)}()
    pnames = [m.fitted_params; :Keq; :E_total]
    M = Gen.define!(Symbol(VARIANT, "_", MODE, "_", m.id), L, m.metabolites, pnames)
    s = m.samples[1]
    c, p = nt(s.concs), nt(s.params)
    s0 = snapshot()
    t_first = @elapsed v = Base.invokelatest(Gen.rate_equation, M, c, p)
    s1 = snapshot()
    err = maximum(abs(Base.invokelatest(Gen.rate_equation, M, nt(x.concs), nt(x.params)) - x.v) /
                  abs(x.v) for x in m.samples)
    a = Base.invokelatest(allocs_of, Gen.rate_equation, M, c, p)
    t1 = Base.invokelatest(t1_of, Gen.rate_equation, M, c, p)
    rows = rows_for(m.metabolites, 1000 + k)
    t2 = Base.invokelatest(t2_of, Gen.rate_equation, M, rows, p)
    ci = code_typed(Gen.evaluate, Tuple{typeof(L), NTuple{length(m.metabolites), Float64},
                                         NTuple{length(pnames), Float64}})[1].first
    record || return nothing
    (id = m.id, variant = VARIANT, mode = MODE, t_law = t_law, t_first = t_first,
     compile_s = (s1.compile - s0.compile) / 1e9, jit = s1.jit - s0.jit, pool = s1.pool - s0.pool,
     allocs = a, t1_ns = t1 * 1e9, t2_ns = t2 * 1e9, maxrelerr = err,
     nops = countops(Gen.body(Gen.decode(Tuple(code_)))), nstmts = length(ci.code),
     ints = length(code_), code_hash = hash(code_[1:7] ∪ code_[11:end]), code = code_)
end

# warm-up on mechanisms outside the measurement set
for (k, m) in enumerate(Pop.load_population("init_unibi"))
    run_one(m, -k; record = false)
end
rows = []
for (k, (pop, m)) in enumerate(selected())
    r = run_one(m, k)
    push!(rows, merge((pop = pop,), r))
    @printf("%-16s first %6.3fs comp %6.3fs jit %7d allocs %d t1 %6.1fns t2 %6.1fns ops %5d stmts %5d err %.1e\n",
            m.id, r.t_first, r.compile_s, r.jit, r.allocs, r.t1_ns, r.t2_ns, r.nops, r.nstmts, r.maxrelerr)
    flush(stdout)
end
serialize(joinpath(@__DIR__, "data", "bench_$(VARIANT)_$(MODE).jls"), rows)
println("done; maxrss ", Sys.maxrss() / 2^20, " MB")
