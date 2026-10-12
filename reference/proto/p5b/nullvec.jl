# ABOUTME: Finds the null direction of dv/dlog(theta) for chosen allosteric mechanisms (probe P5b).
# ABOUTME: Usage: julia nullvec.jl <population> <id...>; prints the fitted names with the null vector.
using Random, LinearAlgebra
include(joinpath(@__DIR__, "Allo.jl"))
pop = Pop.load_population(ARGS[1])
for id in ARGS[2:end]
    m = only(filter(x -> x.id == id, pop))
    md = Allo.derive(m)
    fit = md.names.fitted
    rng = Xoshiro(7)
    setprecision(BigFloat, 256) do
        θ0 = [BigFloat(10)^(4 * big(rand(rng)) - 2) for _ in fit]; keq = big(3.7)
        f(θ, c) = Allo.rate(md, Allo.all_values(md, merge(Dict(fit .=> θ), Dict(:Keq => keq))), c)
        h = big(10)^-25; npts = 2length(fit) + 6
        J = zeros(BigFloat, npts, length(fit))
        for i in 1:npts
            c = Dict(x => BigFloat(10)^(4 * big(rand(rng)) - 2) for x in m.metabolites)
            for k in eachindex(fit)
                θp = copy(θ0); θp[k] *= exp(h); θm = copy(θ0); θm[k] *= exp(-h)
                J[i, k] = (f(θp, c) - f(θm, c)) / (2h)
            end
        end
        F = svd(Float64.(J ./ maximum(abs.(J))))
        println("== $id ", m.source, "  singular values: ", round.(F.S; sigdigits = 3))
        v = F.V[:, end]; v ./= v[argmax(abs.(v))]
        for (n, x) in zip(fit, v)
            abs(x) > 1e-8 && println("   ", rpad(n, 30), round(x; digits = 6))
        end
    end
end
