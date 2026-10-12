# ABOUTME: Checks whether my law's rank deficits are the one-subunit L phantom: the rank with L held
# ABOUTME: fixed equals the number of other fitted constants (probe P5b); also prints sample laws.
using Random, Serialization
include(joinpath(@__DIR__, "Allo.jl"))
include(joinpath(@__DIR__, "validate_lib.jl"))
function rank_gap(J; tol = big(10)^-35)
    A = copy(J); m, n = size(A); piv = BigFloat[]
    for k in 1:min(m, n)
        _, idx = findmax(abs.(A[k:m, k:n])); i, j = idx[1] + k - 1, idx[2] + k - 1
        A[[k, i], :] = A[[i, k], :]; A[:, [k, j]] = A[:, [j, k]]
        push!(piv, abs(A[k, k])); iszero(A[k, k]) && break
        for r in (k + 1):m; f = A[r, k] / A[k, k]; A[r, k:n] .-= f .* A[k, k:n]; end
    end
    count(p -> p > tol * piv[1], piv)
end
function jac(md, rng)
    setprecision(BigFloat, 256) do
        fit = md.names.fitted; θ0 = [BigFloat(10)^(4 * big(rand(rng)) - 2) for _ in fit]; keq = big(2.3)
        f(θ, c) = Allo.rate(md, Allo.all_values(md, merge(Dict(fit .=> θ), Dict(:Keq => keq))), c)
        h = big(10)^-25; J = zeros(BigFloat, 2length(fit) + 6, length(fit))
        for i in axes(J, 1)
            c = Dict(x => BigFloat(10)^(4 * big(rand(rng)) - 2) for x in md.m.metabolites)
            for k in eachindex(fit)
                θp = copy(θ0); θp[k] *= exp(h); θm = copy(θ0); θm[k] *= exp(-h)
                J[i, k] = (f(θp, c) - f(θm, c)) / (2h)
            end
        end
        J, fit
    end
end
rng = Xoshiro(17)
for name in ("expand_bibi", "allo_bibi", "allo_uniuni", "specs", "handwritten")
    rs = deserialize(joinpath(@__DIR__, "data", "results_$(name).jls"))
    P = Dict(m.id => m for m in Pop.load_population(name) if m.allosteric)
    tab = Dict{Any, Int}()
    for r in rs
        r[:count] - r[:rank_mine][1] > 0 || continue
        md = Allo.derive(P[r[:id]])
        J, fit = jac(md, rng)
        keep = findall(!=(:L), fit)
        k = (n = P[r[:id]].cat_multiplicity, deficit = r[:count] - r[:rank_mine][1],
             deficit_without_L = length(keep) - setprecision(() -> rank_gap(J[:, keep]), BigFloat, 256))
        tab[k] = get(tab, k, 0) + 1
    end
    println(name, ": ", sort!(collect(tab); by = string ∘ first))
end
