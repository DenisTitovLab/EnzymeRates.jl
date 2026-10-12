# ABOUTME: P4 census 3: per-metabolite K-factorization of the constant-term-normalized denominator:
# ABOUTME: is every γ_m = c_m c_0^(d-1) / Π c_X^(m_X) a pure number (Schwartz-Zippel, mod 2^61-1)?
include(joinpath(@__DIR__, "P4.jl"))
using Random
using .Poly: Mono, mono, P, mexp, nterms
const PI = P{Int}
include(joinpath(@__DIR__, "Fp.jl"))
rng = Xoshiro(9)
ncon = 0; missingK = 0; factorizable = String[]; partial = 0
for pop in Pop.populations(), m in Pop.load_population(pop)
    m.allosteric && continue
    L = P4.law(m)
    D = P4.cmap(L, L.den)
    haskey(D, Mono()) || continue
    global ncon += 1
    lin = Dict(v => get(D, mono(v => 1), PI()) for v in values(L.conc))
    used = Set(v for mo in keys(D) for (v, _) in mo)
    if any(v -> iszero(lin[v]), used)
        global missingK += 1
        continue
    end
    pts = [[Fp(rand(rng, 1:Int(PR - 1))) for _ in L.vars.names] for _ in 1:2]
    γ(mo, x) = evalp(D[mo], x) * evalp(D[Mono()], x)^(P4.degree(mo) - 1) /
               foldl(*, (evalp(lin[v], x)^e for (v, e) in mo); init = Fp(1))
    numeric = all(mo -> γ(mo, pts[1]) == γ(mo, pts[2]), keys(D))
    numeric ? push!(factorizable, m.id) : (global partial += 1)
end
println("constant term: ", ncon, "; some metabolite in D has no linear term (K_X undefined): ", missingK)
println("every γ_m a pure number (fully K-factorizable): ", length(factorizable), "; not: ", partial)
println("fully K-factorizable by population: ", Dict(p => count(id -> startswith(id, p * "/"), factorizable)
        for p in unique(first.(split.(factorizable, "/")))))
println("examples: ", first(factorizable, 12))
