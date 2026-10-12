# ABOUTME: Codegen check (c): method instances and native code that an identify-like flow
# ABOUTME: (PASS1, a fit and one LOOCV fold, as test_compile_budget.jl) compiles per mechanism.
const SCRATCH = "/tmp/claude-501/-home-denis-linux--julia-dev-EnzymeRates/12f42995-f5ba-4fda-9510-dd9b87f7d069/scratchpad"
include(joinpath(SCRATCH, "proto", "common", "Population.jl"))
include(joinpath(SCRATCH, "slice", "src", "Slice.jl"))
using Serialization, Random, Printf
const S = Slice
load(set) = [c.mech for c in deserialize(joinpath(SCRATCH, "slice", "core", "data", set * ".jls"))]

"""The package's FittingProblem, stdlib only: data columns, groups, Keq, log rates, buffer."""
struct FitProblem{M, D}
    mechanism::M
    data::D
    groups::Vector{Vector{Int}}
    Keq::Float64
    log_abs_rates::Vector{Float64}
    buf::Vector{Float64}
end
FitProblem(em, data, Keq) = FitProblem(em, data, [findall(==(g), data.group)
                                                  for g in unique(data.group)],
                                       Float64(Keq), log.(abs.(data.Rate)),
                                       zeros(length(data.Rate)))

"""The package's loss! (relative mode), line for line: per-group centered squared
log-ratios plus a flat penalty per sign mismatch."""
function loss!(x::Vector{Float64}, fp::FitProblem)
    buf = fp.buf
    ParamNames = S.fitted_params(fp.mechanism)
    MetNames = S.metabolites(fp.mechanism)
    N = length(ParamNames)
    K = length(MetNames)
    fitted = NamedTuple{ParamNames}(ntuple(i -> exp(x[i]), Val(N)))
    params = merge(fitted, (Keq = fp.Keq, E_total = 1.0))
    n_mismatch = 0
    @inbounds for i in eachindex(buf)
        concs = NamedTuple{MetNames}(ntuple(j -> getproperty(fp.data, MetNames[j])[i], Val(K)))
        pred = S.rate_equation(fp.mechanism, concs, params)
        if sign(pred) != sign(fp.data.Rate[i])
            buf[i] = 10.0
            n_mismatch += 1
        else
            buf[i] = log(abs(pred)) - fp.log_abs_rates[i]
        end
    end
    total = 0.0
    @inbounds for g in fp.groups
        mean_lr = 0.0
        for j in g
            mean_lr += buf[j]
        end
        mean_lr /= length(g)
        for j in g
            d = buf[j] - mean_lr
            total += d * d
        end
    end
    (total + 100.0 * n_mismatch) / length(buf)
end

"""A random-search optimizer compiled once: the objective is an untyped field read behind a
function barrier (the package hides it behind a FunctionWrapper for the same reason)."""
struct Objective
    f::Any
end
function minimize(obj::Objective, np::Int, iters::Int, rng)
    best = clamp.(2 .* randn(rng, np), -15, 15)
    fbest = obj.f(best)::Float64
    for _ in 1:iters, _ in 1:8
        x = clamp.(best .+ randn(rng, np), -15, 15)
        fx = obj.f(x)::Float64
        fx < fbest && ((best, fbest) = (x, fx))
    end
    best, fbest
end

fit(fp::FitProblem) = minimize(Objective(Base.Fix2(loss!, fp)),
                               length(S.fitted_params(fp.mechanism)), 20, Xoshiro(1))

instances(meths) = sum(m -> count(_ -> true, Base.specializations(m)), meths)
function all_instances()
    meths = Method[]
    Base.visit(m -> push!(meths, m), Core.methodtable)
    instances(meths)
end

"""Instances and native bytes PASS1, a fit and a LOOCV fold of `m` add, and compile time."""
function measure(m, rng)
    n0, b0 = all_instances(), Base.jit_total_bytes()
    t0 = Base.cumulative_compile_time_ns()[1]
    S.complexity(m)                                     # PASS1
    em = S.compile(m)
    S.fitted_params(em)
    S.rate_equation_string(m)
    names = S.metabolites(em)
    data = (; group = repeat(["G1", "G2"], inner = 4), Rate = 0.1 .+ rand(rng, 8),
            (x => 10.0 .^ (2 .* rand(rng, 8) .- 1) for x in names)...)
    fit(FitProblem(em, data, 2.0))                      # PASS2
    keep = data.group .!= "G1"                          # one LOOCV fold
    fit(FitProblem(em, map(c -> c[keep], data), 2.0))
    (all_instances() - n0, Base.jit_total_bytes() - b0,
     (Base.cumulative_compile_time_ns()[1] - t0) / 1e9)
end

Base.cumulative_compile_timing(true)
rng = Xoshiro(1)
ms = load("init_bibi")
n_sig, b_sig = all_instances(), Base.jit_total_bytes()
foreach(S.sig, ms)
@printf("Sig encoding of the %d init_bibi seeds: %d instances, %.1f KB\n", length(ms),
        all_instances() - n_sig, (Base.jit_total_bytes() - b_sig) / 1024)
# The first mechanism of each kind compiles what no mechanism owns (generator code paths a
# kind first takes); measured, then the fresh ones.
warm = [("init_bibi", 1), ("expand_bibi", 600), ("allo_bibi", 700), ("specs", 41)]
cases = [warm; [("init_bibi", i) for i in (2, 6, 11, 18)];
         [("expand_bibi", i) for i in (626, 642, 658)];
         [("allo_bibi", i) for i in (1, 500, 1000, 1400)]; [("specs", i) for i in (34, 36)]]
rows = []
for (set, i) in cases
    m = load(set)[i]
    n, b, t = measure(m, rng)
    push!(rows, (; set, i, n, b, t, fresh = !((set, i) in warm)))
    @printf("%-12s %4d: %4d method instances, %6.1f KB native code, %.3f s compile \
            (allosteric %s, %d locals)\n", set, i, n, b / 1024, t, S.is_allosteric(m),
            length(S.program(m).defs))
end
fr = filter(r -> r.fresh, rows)
@printf("per fresh mechanism: instances max %d (budget 250), native max %.1f KB (budget 300)\n",
        maximum(r -> r.n, fr), maximum(r -> r.b, fr) / 1024)
serialize(joinpath(@__DIR__, "data", "budget.jls"), rows)
