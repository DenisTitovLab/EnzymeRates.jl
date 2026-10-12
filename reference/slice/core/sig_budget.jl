# ABOUTME: Measures what encoding the Sigs of the 239 init_bibi mechanisms compiles, and what a
# ABOUTME: generator that decodes the Sig compiles per fresh mechanism (test_compile_budget way).
include(joinpath(@__DIR__, "..", "src", "Slice.jl"))
using Serialization, Printf
const DATA = joinpath(@__DIR__, "data")

instances(meths) = sum(m -> count(_ -> true, Base.specializations(m)), meths)
function all_instances()
    meths = Method[]
    Base.visit(m -> push!(meths, m), Core.methodtable)
    instances(meths)
end
specs(f) = sum(m -> count(_ -> true, Base.specializations(m)), methods(f))

ms = [c.mech for c in deserialize(joinpath(DATA, "init_bibi.jls"))]
println("init_bibi mechanisms: ", length(ms))

# Encoding every seed's Sig, once per process, as the package's gate measures `_sig_of`.
n0, b0 = all_instances(), Base.jit_total_bytes()
sigs = map(Slice.sig, ms)
@printf("SIG_INSTANCES:%d SIG_NATIVE:%d (%.1f KB); distinct Sig lengths %d\n",
        all_instances() - n0, Base.jit_total_bytes() - b0,
        (Base.jit_total_bytes() - b0) / 1024, length(unique(length, sigs)))

# Encoding more mechanisms of new Sig lengths afterwards compiles nothing more.
more = [c.mech for set in ("expand_bibi", "allo_bibi", "allo_terbi")
        for c in deserialize(joinpath(DATA, set * ".jls"))]
n0, b0 = all_instances(), Base.jit_total_bytes()
foreach(Slice.sig, more)
dn, db = all_instances() - n0, Base.jit_total_bytes() - b0
@printf("then %d more mechanisms (expand_bibi, allo_bibi, allo_terbi; %d new Sig lengths): \
         %d instances, %d B\n", length(more),
        length(setdiff(unique(length ∘ Slice.sig, more), unique(length, sigs))), dn, db)

# A generator that decodes its Sig by splatting, as `rate_equation` will.
@generated function nedges(::Slice.EnzymeMechanism{S}) where {S}
    length(Slice.decode(S...).edges)
end
Base.cumulative_compile_timing(true)
function measure(m)
    n0, b0 = all_instances(), Base.jit_total_bytes()
    t0 = Base.cumulative_compile_time_ns()[1]
    em = Slice.compile(m)
    n1, b1 = all_instances(), Base.jit_total_bytes()
    k = nedges(em)
    k == length(m.edges) || error("decoded the wrong mechanism")
    (n1 - n0, b1 - b0, all_instances() - n1, Base.jit_total_bytes() - b1,
     (Base.cumulative_compile_time_ns()[1] - t0) / 1e9)
end
measure(ms[1])
fresh = [2, 6, 11, 18]
for (k, i) in enumerate(fresh)
    a, b, c, d, t = measure(ms[i])
    @printf("fresh %d (init_bibi/%d, Sig length %d): compile(m) %d instances %d B; \
            generated call %d instances %d B; compile time %.4f s\n",
            k, i, length(sigs[i]), a, b, c, d, t)
end
rest = setdiff(eachindex(ms), [1; fresh])
tot = map(i -> measure(ms[i]), rest)
@printf("the other %d mechanisms: per mechanism compile(m) %.2f instances, generated call \
         %.2f instances, %.0f B native; max %d instances\n", length(rest),
        sum(x -> x[1], tot) / length(rest), sum(x -> x[3], tot) / length(rest),
        sum(x -> x[2] + x[4], tot) / length(rest), maximum(x -> x[1] + x[3], tot))
@printf("specializations: decode %d, sig %d, mechanism %d, compile %d\n", specs(Slice.decode),
        specs(Slice.sig), specs(Slice.mechanism), specs(Slice.compile))
gen = only(methods(nedges)).generator.gen
@printf("generator specializations %d for %d distinct Sig lengths among %d mechanisms\n",
        specs(gen), length(unique(length, sigs)), length(ms))
n0 = all_instances()
foreach(m -> nedges(Slice.compile(m)), more)
dn = all_instances() - n0
@printf("generated call on the %d further mechanisms: %d instances in all (%.2f each)\n",
        length(more), dn, dn / length(more))
mechs_back = map(Slice.mechanism ∘ Slice.compile, ms)
@printf("mechanism(compile(m)) == m for %d of %d\n", count(mechs_back .== ms), length(ms))
