# ABOUTME: Does a kernel compiled for one mechanism serve another whose positional law is identical?
# ABOUTME: Per class of identical laws (expand_bibi, allo_bibi), first-member vs later-member compile cost.
include(joinpath(@__DIR__, "Encode.jl"))
include(joinpath(@__DIR__, "Measure.jl"))
using Statistics
const VARIANT = ARGS[1]
mk(m) = VARIANT == "flat" ? P7.flat_law(m) : P7.factored_law(m; expand_den = VARIANT == "halh")
opt = VARIANT == "flat" ? (false, false, 0) : (true, true, 0)
classes = Dict{Vector{Int}, Vector{Pop.PMech}}()
for pop in ("expand_bibi", "allo_bibi"), m in Pop.load_population(pop)
    push!(get!(classes, code(mk(m), opt), Pop.PMech[]), m)
end
multi = sort([v for v in values(classes) if length(v) > 1]; by = v -> v[1].id)
println("$VARIANT: $(length(multi)) classes with ≥ 2 members ($(sum(length, multi)) mechanisms)")
function first_call(m, L)
    M = Gen.define!(Symbol("share_", m.id), L, m.metabolites, [m.fitted_params; :Keq; :E_total])
    s = m.samples[1]; c, p = nt(s.concs), nt(s.params)
    s0 = snapshot()
    v = Base.invokelatest(Gen.rate_equation, M, c, p)
    s1 = snapshot()
    abs(v - s.v) <= 1e-12 * abs(s.v) || error("$(m.id): value")
    ((s1.compile - s0.compile) / 1e9, s1.jit - s0.jit)
end
# warm-up: one class outside the measured ones is enough to compile the harness
first_call(Pop.load_population("init_unibi")[1], encode(mk(Pop.load_population("init_unibi")[1]), opt))
firsts, laters = Tuple{Float64, Int}[], Tuple{Float64, Int}[]
for v in multi
    L = encode(mk(v[1]), opt)
    push!(firsts, first_call(v[1], L))
    for m in v[2:end]
        typeof(encode(mk(m), opt)) === typeof(L) || error("not identical")
        push!(laters, first_call(m, L))
    end
end
println("first member (kernel + wrapper): compile median $(round(median(first.(firsts)) * 1000, digits = 2)) ms, ",
        "JIT median $(median(last.(firsts))) B (n=$(length(firsts)))")
println("later members (wrapper only):   compile median $(round(median(first.(laters)) * 1000, digits = 2)) ms ",
        "(max $(round(maximum(first.(laters)) * 1000, digits = 2))), JIT median $(median(last.(laters))) B (n=$(length(laters)))")
