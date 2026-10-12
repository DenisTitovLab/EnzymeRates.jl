# ABOUTME: Checks which parameter-scaling directions leave the package's reduced A5 law unchanged.
# ABOUTME: Evaluates v from the reduced terms at θ and at θ with chosen parameters scaled by 3.
include(joinpath(@__DIR__, "SkLib.jl"))
using Random, Printf, Serialization
pop = deserialize(joinpath(@__DIR__, "data", "adv_pkg.jls"))
m = only(filter(x -> x.source == "A5_cyclic_dead_end_ss_square", pop))
θ = m.fitted_params
println(collect(enumerate(θ)))
println(m.dependent_monomials)
setprecision(BigFloat, 256) do
    N = Sk.reduced_terms(m.pkg_num, m, θ); D = Sk.reduced_terms(m.pkg_den, m, θ)
    rng = Random.Xoshiro(3)
    θv = [exp(big(log(10.0)) * (4 * big(rand(rng)) - 2)) for _ in θ]
    Keq = big(3)
    c = [exp(big(log(10.0)) * (6 * big(rand(rng)) - 3)) for _ in m.metabolites]
    val(ts, θv) = sum(BigFloat(numerator(t[1])) / denominator(t[1]) *
                      prod((θv[i]^t[2][i] for i in eachindex(θv) if t[2][i] != 0); init = big(1)) *
                      Keq^t[3] * prod((c[i]^t[4][i] for i in eachindex(c) if t[4][i] != 0); init = big(1))
                      for t in ts)
    v(θv) = val(N, θv) / val(D, θv)
    v0 = v(θv)
    for dir in ([8], [5, 9], [1, 3], [2, 4], [6], [7], [10], [1], [3])
        θ2 = copy(θv); θ2[dir] .*= 3
        @printf("scale %-8s (%s): rel change %.2e\n", string(dir), join(θ[dir], ","), Float64(abs(v(θ2) - v0) / abs(v0)))
    end
end
