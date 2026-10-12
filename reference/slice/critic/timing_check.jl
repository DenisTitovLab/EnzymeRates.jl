# ABOUTME: Checks the speed protocol on the largest laws: fixed inputs against 64 varying inputs,
# ABOUTME: and the floating-point instructions LLVM keeps against the generated body's call count.
const SCR = "/tmp/claude-501/-home-denis-linux--julia-dev-EnzymeRates/12f42995-f5ba-4fda-9510-dd9b87f7d069/scratchpad"
include(joinpath(SCR, "proto", "common", "Population.jl"))
include(joinpath(SCR, "slice", "src", "Slice.jl"))
using Serialization, Random, InteractiveUtils, Printf
const SL = Slice

"""Number of `:call` expressions in a generated body (the integration stage's ops)."""
countops(b) = (n = Ref(0); w(x) = (x isa Expr && (x.head === :call && (n[] += 1);
                                                   foreach(w, x.args))); w(b); n[])

@noinline barrier(em, c, p) = SL.rate_equation(em, c, p)

"""Minimum over 5 batches of 10,000 calls (ns per call), inputs fixed or cycling 64 points."""
function per_call(em, pts, vary::Bool)
    best, acc = Inf, 0.0
    for _ in 1:5
        t = @elapsed for i in 1:10_000
            c, p = vary ? pts[(i & 63) + 1] : pts[1]
            acc += barrier(em, c, p)
        end
        best = min(best, t / 10_000)
    end
    isfinite(acc) || error("non-finite")
    1e9 * best
end

function check(label, m)
    rng = Xoshiro(7)
    cn, pn = Tuple(SL.metabolites(m)), SL.parameters(m)
    pts = [(NamedTuple{cn}(Tuple(10.0^(4 * rand(rng) - 2) for _ in cn)),
            NamedTuple{pn}(Tuple(x === :E_total ? 1.0 : 10.0^(4 * rand(rng) - 2) for x in pn)))
           for _ in 1:64]
    em = SL.compile(m)
    t_first = @elapsed barrier(em, pts[1]...)
    allocs = @allocated barrier(em, pts[1]...)
    ops = countops(SL.body(SL.program(m)))
    ir = sprint(io -> code_llvm(io, SL.rate_equation, typeof.((em, pts[1]...));
                                debuginfo = :none))
    fp = Dict(op => count(l -> occursin(Regex("= $op (fast )?double"), l), split(ir, '\n'))
              for op in ("fmul", "fadd", "fsub", "fdiv", "fneg"))
    calls = count(l -> occursin("@llvm.sqrt", l) || occursin("@llvm.pow", l) ||
                       occursin("@llvm.fma", l) || occursin("@llvm.fmuladd", l),
                  split(ir, '\n'))
    nfp = sum(values(fp))
    tf, tv = per_call(em, pts, false), per_call(em, pts, true)
    @printf("%-22s ops %5d | LLVM fp %5d (%s; sqrt/pow/fma calls %d) | fixed %.1f ns, varying %.1f ns | %.1f fp/ns | allocs %d | first call %.2f s\n",
            label, ops, nfp, join(("$k $(fp[k])" for k in sort!(collect(keys(fp)))), " "),
            calls, tf, tv, nfp / tv, allocs, t_first)
    flush(stdout)
end

conv = deserialize(joinpath(SCR, "slice", "core", "data", "expand_bibi.jls"))
for id in ("expand_bibi/642",)
    c = only(filter(x -> x.id == id, conv))
    check(id, c.mech)
end
inh = deserialize(joinpath(SCR, "slice", "core", "data", "p1_enum_inh.jls"))
check("bi_bi_with_inh/1013", only(filter(x -> x.id == "bi_bi_with_inh/1013", inh)).mech)
ps = deserialize(joinpath(SCR, "slice", "integrate", "data", "pop_init_terter.part3.jls"))
i = findfirst(p -> p.id == "init_terter/62921", ps)
check("init_terter/62921", SL.Converted(ps[i], i).mech)
println("done")
