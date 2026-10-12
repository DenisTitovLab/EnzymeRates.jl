# ABOUTME: The textbook random-order steady-state bi-bi (every step its own group): how often v(A) at
# ABOUTME: fixed B peaks above v(A -> infinity, same B), against how often sup v exceeds kcat.
include(joinpath(@__DIR__, "klib.jl"))

F(bound...) = S.RawForm(:E, [S.Met(b, false) for b in bound], Symbol[], Symbol[])
M(x) = [S.Met(x, false)]
none = S.Met[]
steps = [S.RawStep(F(), F(:A), M(:A), none, false), S.RawStep(F(), F(:B), M(:B), none, false),
         S.RawStep(F(:A), F(:A, :B), M(:B), none, false), S.RawStep(F(:B), F(:A, :B), M(:A), none, false),
         S.RawStep(F(:A, :B), F(:P, :Q), none, none, false), S.RawStep(F(:P, :Q), F(:Q), none, M(:P), true),
         S.RawStep(F(:Q), F(), none, M(:Q), true)]
m = S.Mechanism([:A, :B], [:P, :Q], steps, collect(1:7))
_, _, sol, L = derivable(m)
_, _, vars, terms = kcat_class(m, sol)
println("zero products: ", zero_product_limits(m, L))
rng = Xoshiro(11)
npts, fixedB_peak, over = 300, 0, 0
best = (0.0, nothing)
for _ in 1:npts
    names = S.column_names(m, sol)
    x = Dict{Symbol, BigFloat}(n => big(10.0^(-3 + 6 * rand(rng))) for n in S.parameters(m))
    nl = NLaw(m, L, x)
    kc = kcat_value(vars, terms, x)
    hit = false
    for B in 10.0 .^ (-3:1.0:3)
        vs = [rate(nl, [10.0^a, B, 0.0, 0.0]) for a in -6:0.05:12]
        q = maximum(vs) / vs[end]
        q > 1.01 && (hit = true)
        q > best[1] && (global best = (q, (B, x)))
    end
    global fixedB_peak += hit
    sup, _ = sup_search(nl, rng; lo = -6.0, hi = 6.0)
    global over += sup > 1.01 * kc
end
@printf("%d parameter points: v(A) at some fixed B peaks > 1%% above v(A -> 1e12): %d; sup over A, B > 1.01 kcat: %d; largest peak ratio %.3f\n",
        npts, fixedB_peak, over, best[1])
