# ABOUTME: Debugs the identifiable-rank helper on adversarial mechanisms with a known answer (A5, A14, A1).
# ABOUTME: Prints Float64 singular values of the column-normalized Jacobian next to the elimination rank.
include(joinpath(@__DIR__, "SkLib.jl"))
using Random, Printf, Serialization, LinearAlgebra
pop = deserialize(joinpath(@__DIR__, "data", "adv_pkg.jls"))
for id in ("A1_iso_parallel_ss", "A14_uniuni_two_intermediates_all_ss", "A5_cyclic_dead_end_ss_square",
           "A2_iso_parallel_re_iso")
    m = only(filter(x -> x.source == id, pop))
    θ = m.fitted_params
    nθ = length(θ)
    rng = Random.Xoshiro(5)
    setprecision(BigFloat, 256) do
        N = Sk.reduced_terms(m.pkg_num, m, θ); D = Sk.reduced_terms(m.pkg_den, m, θ)
        θv = [exp(big(log(10.0)) * (4 * big(rand(rng)) - 2)) for _ in 1:nθ]
        Keq = big(3)
        tv(t) = BigFloat(numerator(t[1])) / denominator(t[1]) *
                prod((θv[i]^t[2][i] for i in 1:nθ if t[2][i] != 0); init = big(1)) * Keq^t[3]
        Nθ = map(tv, N); Dθ = map(tv, D)
        np = 3 * nθ + 6
        J = zeros(BigFloat, np, nθ)
        for p in 1:np
            c = [exp(big(log(10.0)) * (6 * big(rand(rng)) - 3)) for _ in m.metabolites]
            cm(t) = prod((c[i]^t[4][i] for i in eachindex(c) if t[4][i] != 0); init = big(1))
            nv = [Nθ[k] * cm(N[k]) for k in eachindex(N)]
            dv = [Dθ[k] * cm(D[k]) for k in eachindex(D)]
            Nt, Dt = sum(nv), sum(dv)
            for j in 1:nθ
                Nj = sum(nv[k] * N[k][2][j] for k in eachindex(N))
                Dj = sum(dv[k] * D[k][2][j] for k in eachindex(D))
                J[p, j] = (Nj * Dt - Nt * Dj) / Dt^2
            end
        end
        Jn = copy(J)
        for j in 1:nθ; Jn[:, j] ./= sqrt(sum(abs2, Jn[:, j])); end
        sv = svdvals(Float64.(Jn))
        println(id, ": fitted ", nθ, ", elimination rank ", Sk.numeric_rank(J), ", Float64 singular values ",
                join([@sprintf("%.1e", s) for s in sv], " "))
        println("   fitted: ", θ)
        println("   num terms ", length(N), " den terms ", length(D), "; metabolites ", m.metabolites)
        # Rows: show the exponent table of θ in den terms (first 12)
        for t in D[1:min(end, 12)]
            println("     den term coef ", t[1], " θ-exps ", t[2], " Keq ", t[3], " conc ", t[4])
        end
    end
end
