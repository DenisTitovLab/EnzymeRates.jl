# ABOUTME: P7 benchmark of the package's own rate_equation on the P7 measurement set: first-call
# ABOUTME: compile, JIT bytes, generator (derivation) time, allocations, T1/T2 timing, op counts.
using EnzymeRates, Random, Serialization, Printf, InteractiveUtils
const ER = EnzymeRates
include(joinpath(@__DIR__, "Select.jl"))
include(joinpath(@__DIR__, "Measure.jl"))
const TESTDIR = joinpath(dirname(dirname(pathof(EnzymeRates))), "test")
for ex in Meta.parseall(read(joinpath(TESTDIR, "test_mechanism_enumeration.jl"), String)).args
    ex isa Expr && ex.head === :const && ex.args[1] isa Expr &&
        ex.args[1].args[1] in (:uni_bi_rxn, :bi_bi_rxn, :bi_bi_pp_rxn) && Core.eval(Main, ex)
end
include(joinpath(TESTDIR, "mechanism_definitions_for_test_enzyme_derivation.jl"))

# The populations rebuilt exactly as common/export/export_pop.jl builds them.
init_bibi = unique(ER.init_mechanisms(bi_bi_rxn))
init_pp = unique(ER.init_mechanisms(bi_bi_pp_rxn))
pick = sort!(randperm(Random.Xoshiro(20261009), length(init_pp))[1:60])
children = unique(reduce(vcat, [ER.expand_mechanisms([p], bi_bi_pp_rxn) for p in init_pp[pick]]))
length(children) == 700 || error("expand_bibi rebuilt $(length(children)) children, not 700")
function package_mech(pop, m)
    i = parse(Int, split(m.id, "/")[2])
    x = pop == "specs" ? ER._concrete(MECHANISM_TEST_SPECS[i].mechanism) :
        pop == "init_bibi" ? init_bibi[i] : pop == "expand_bibi" ? children[i] :
        unique(ER.init_mechanisms(uni_bi_rxn))[i]
    em = ER.compile_mechanism(x)
    pop == "init_unibi" || ER.rate_equation_string(em) == m.rate_equation_string ||
        error("$(m.id): rebuilt mechanism differs from the exported one")
    em
end
countops(b) = (n = Ref(0); w(x) = (x isa Expr && (x.head === :call && (n[] += 1); foreach(w, x.args))); w(b); n[])

function run_one(pop, m, k; record = true)
    em = package_mech(pop, m)
    M = typeof(em)
    s = m.samples[1]
    c, p = nt(s.concs), nt(s.params)
    s0 = snapshot()
    t_first = @elapsed v = Base.invokelatest(ER.rate_equation, em, c, p)
    s1 = snapshot()
    err = maximum(abs(ER.rate_equation(em, nt(x.concs), nt(x.params)) - x.v) / abs(x.v) for x in m.samples)
    t_gen = minimum(@elapsed(begin
        dep, indep = ER._dependent_param_exprs(M)
        ER._rate_body(M, (indep..., :Keq, :E_total), dep)
    end) for _ in 1:3)
    dep, indep = ER._dependent_param_exprs(M)
    b = ER._rate_body(M, (indep..., :Keq, :E_total), dep)
    a = allocs_of(ER.rate_equation, em, c, p)
    t1 = t1_of(ER.rate_equation, em, c, p)
    rows = rows_for(m.metabolites, 1000 + k)
    t2 = t2_of(ER.rate_equation, em, rows, p)
    ci = code_typed(ER.rate_equation, Tuple{M, typeof(c), typeof(p), typeof(ER.Reduced)})[1].first
    record || return nothing
    (pop = pop, id = m.id, variant = "pkg", mode = "unique", t_first = t_first,
     compile_s = (s1.compile - s0.compile) / 1e9, jit = s1.jit - s0.jit, pool = s1.pool - s0.pool,
     t_gen = t_gen, allocs = a, t1_ns = t1 * 1e9, t2_ns = t2 * 1e9, maxrelerr = err,
     nops = countops(b), nstmts = length(ci.code), sig_bytes = Base.summarysize(M))
end

for (k, m) in enumerate(Pop.load_population("init_unibi"))
    run_one("init_unibi", m, -k; record = false)
end
rows = []
for (k, (pop, m)) in enumerate(selected())
    r = run_one(pop, m, k)
    push!(rows, r)
    @printf("%-16s first %6.3fs comp %6.3fs gen %6.4fs jit %7d allocs %d t1 %6.1fns t2 %6.1fns ops %5d stmts %5d err %.1e\n",
            m.id, r.t_first, r.compile_s, r.t_gen, r.jit, r.allocs, r.t1_ns, r.t2_ns, r.nops, r.nstmts, r.maxrelerr)
    flush(stdout)
end
serialize(joinpath(@__DIR__, "data", "bench_pkg_unique.jls"), rows)
println("done; maxrss ", Sys.maxrss() / 2^20, " MB")
