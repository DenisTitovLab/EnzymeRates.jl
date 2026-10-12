# ABOUTME: C5 follow-up: rank of the normalized-coefficient map of every plain law against its identifiable
# ABOUTME: rank (skeptic P1 census), coefficient counts for chosen specs, and Cleland ordered bi-bi checks.
include("/tmp/claude-501/-home-denis-linux--julia-dev-EnzymeRates/12f42995-f5ba-4fda-9510-dd9b87f7d069/scratchpad/proto/skeptic_P1-engine/SkLib.jl")
using Random, Printf, Serialization

const SKDATA = "/tmp/claude-501/-home-denis-linux--julia-dev-EnzymeRates/12f42995-f5ba-4fda-9510-dd9b87f7d069/scratchpad/proto/skeptic_P1-engine/data/rank_pop.jls"
const POPS = ["specs", "init_uniuni", "init_unibi", "init_bibi", "init_bibi_pp", "expand_bibi",
              "handwritten"]
logu(rng) = exp(big(log(10.0)) * (4 * big(rand(rng)) - 2))

"""Numerator and denominator coefficients of the Reduced law of `m` (dependents substituted),
one per concentration monomial: `(n_terms_N, n_terms_D, rank)`, where rank is the generic rank of
∂ log(c_m / c_ref)/∂ log θ over every coefficient c_m, with c_ref a lowest-degree D coefficient."""
function coef_rank(m, rng)
    θ = m.fitted_params
    nθ = length(θ)
    N = Sk.reduced_terms(m.pkg_num, m, θ); D = Sk.reduced_terms(m.pkg_den, m, θ)
    setprecision(BigFloat, 256) do
        θv = [logu(rng) for _ in 1:nθ]; Keq = logu(rng)
        coefs = Dict{Tuple{Bool, Vector{Int}}, Tuple{BigFloat, Vector{BigFloat}}}()
        for (isN, T) in ((true, N), (false, D)), t in T
            val = BigFloat(numerator(t[1])) / denominator(t[1]) *
                  prod((θv[i]^t[2][i] for i in 1:nθ if t[2][i] != 0); init = big(1)) * Keq^t[3]
            c, g = get(coefs, (isN, t[4]), (big(0), zeros(BigFloat, nθ)))
            coefs[(isN, t[4])] = (c + val, g .+ val .* BigFloat.(t[2]))
        end
        ks = collect(keys(coefs))
        ref = argmin(k -> (k[1], sum(k[2])), ks)          # a D coefficient of lowest degree
        rows = [coefs[k][2] ./ coefs[k][1] .- coefs[ref][2] ./ coefs[ref][1] for k in ks if k != ref]
        J = isempty(rows) ? zeros(BigFloat, 0, nθ) : permutedims(reduce(hcat, rows))
        count(k -> k[1], ks), count(k -> !k[1], ks), nθ == 0 ? 0 : Sk.numeric_rank(J)
    end
end

sk = Dict(r.id => r for r in deserialize(SKDATA))
println("== normalized-coefficient rank vs identifiable rank (skeptic P1), all plain mechanisms ==")
agree = 0; total = 0; mism = String[]; rows = []
for pop in POPS
    for (n, m) in enumerate(Pop.load_population(pop))
        (m.allosteric || !m.ok) && continue
        nN, nD, r1 = coef_rank(m, Random.Xoshiro(11 * n + length(pop)))
        _, _, r2 = coef_rank(m, Random.Xoshiro(53 * n + 7))
        rc = max(r1, r2)
        s = sk[m.id]
        global total += 1
        rc == s.rank ? (global agree += 1) : push!(mism, "$(m.id): coef $rc vs v $(s.rank)")
        push!(rows, (id = m.id, source = m.source, fitted = s.nθ, rank = s.rank, nN, nD, rc))
    end
end
println("coefficient-ratio rank == identifiable rank: $agree of $total; mismatches: ", mism)
println("deficient (fitted > rank): ", count(r -> r.fitted > r.rank, rows), " of ", length(rows),
        "; enumerated only: ", count(r -> r.fitted > r.rank && !startswith(r.id, "specs") &&
                                       !startswith(r.id, "handwritten"), rows), " of ",
        count(r -> !startswith(r.id, "specs") && !startswith(r.id, "handwritten"), rows))

println("\n== chosen mechanisms: fitted, rank, numerator/denominator coefficients, ratios ==")
for id in ["specs/1", "specs/2", "specs/5", "specs/7", "specs/11", "specs/13", "specs/20",
           "init_bibi_pp/61", "init_bibi_pp/62", "init_bibi_pp/261", "init_bibi_pp/262",
           "expand_bibi/303", "expand_bibi/304", "expand_bibi/305"]
    r = only(filter(x -> x.id == id, rows))
    @printf("  %-17s %-40s fitted %2d  rank %2d  N %d + D %2d coefficients = %2d ratios, ratio rank %2d\n",
            id, first(r.source, 40), r.fitted, r.rank, r.nN, r.nD, r.nN + r.nD - 1, r.rc)
end

# ── Cleland ordered bi-bi (Segel IX-87) from the 8 microscopic constants (P4's ratios) ──
println("\n== Cleland ordered bi-bi: names vs dimension ==")
function cleland_ordered(k)
    k1f, k1r, k2f, k2r, k3f, k3r, k4f, k4r = k
    V1 = k3f * k4f / (k3f + k4f); V2 = k2r * k1r / (k2r + k1r)
    (V1 = V1, V2 = V2, Ka = V1 / k1f, Kb = k4f * (k2r + k3f) / (k2f * (k3f + k4f)),
     Kp = k1r * (k2r + k3f) / (k3r * (k2r + k1r)), Kq = V2 / k4r, Kia = k1r / k1f,
     Kib = (k2r + k1r) / k2f, Kip = (k3f + k4f) / k3r, Kiq = k4f / k4r,
     Keq = k1f * k2f * k3f * k4f / (k1r * k2r * k3r * k4r))
end
function jac_rank(f, x)                                    # ∂ log f / ∂ log x by central differences
    h = big(10)^-30
    cols = [(log.(collect(f(x .* [i == j ? exp(h) : big(1) for i in eachindex(x)]))) .-
             log.(collect(f(x .* [i == j ? exp(-h) : big(1) for i in eachindex(x)])))) ./ (2h)
            for j in eachindex(x)]
    Sk.numeric_rank(reduce(hcat, cols); tol = big(10)^-20)
end
setprecision(BigFloat, 512) do
    rng = Random.Xoshiro(5)
    k = [logu(rng) for _ in 1:8]
    c = cleland_ordered(k)
    println("rank of (10 Cleland constants, Keq) in the 8 rate constants: ", jac_rank(cleland_ordered, k))
    H1 = c.V1 * c.Kp * c.Kiq / (c.V2 * c.Kia * c.Kb); H2 = (c.V1 / c.V2)^2 * c.Kip * c.Kq / (c.Ka * c.Kib)
    rel = c.Kb * c.Kia * c.V1 / (c.Kib * c.Ka) - (c.V1 + c.V2 - c.V1 * c.Kq / c.Kiq - c.V2 * c.Ka / c.Kia)
    @printf("H1/Keq - 1 = %.1e, H2/Keq - 1 = %.1e, non-monomial relation residual = %.1e\n",
            Float64(H1 / c.Keq - 1), Float64(H2 / c.Keq - 1), Float64(rel / c.Kb))
    # inequalities: k4f = V2*Kiq/Kq > V1 and k1r = V1*Kia/Ka > V2, at 10^5 random rate-constant draws
    viol = 0
    for _ in 1:100_000
        c = cleland_ordered([logu(rng) for _ in 1:8])
        (c.V2 * c.Kiq / c.Kq > c.V1 && c.V1 * c.Kia / c.Ka > c.V2) || (viol += 1)
    end
    println("draws violating V1 < V2*Kiq/Kq or V2 < V1*Kia/Ka: $viol of 100000")
    # The Cleland formula's 12 normalized coefficients are monomials in the constants; its dimension
    # with Keq fixed and H1 (V2 from Keq) or H1 and H2 (V2 and Kip from Keq) imposed.
    coefs(V1, V2, Ka, Kb, Kp, Kq, Kia, Kib, Kip, Kiq, Keq) =
        [V1 * V2, V1 * V2 / Keq, V2 * Kb, V2 * Ka, V1 * Kq / Keq, V1 * Kp / Keq, V2,
         V1 * Kq / (Kia * Keq), V2 * Ka / Kiq, V1 / Keq, V2 / Kip, V1 / (Kib * Keq)] ./
        (V2 * Kia * Kb)
    Keq = logu(rng)
    v2(V1, Ka, Kb, Kp, Kq, Kia, Kiq) = V1 * Kp * Kiq / (Keq * Kia * Kb)
    kip(V1, V2, Ka, Kq, Kib) = Keq * Ka * Kib / ((V1 / V2)^2 * Kq)
    f9(x) = coefs(x[1], v2(x[1], x[2], x[3], x[4], x[5], x[6], x[9]), x[2], x[3], x[4], x[5], x[6],
                  x[7], x[8], x[9], Keq)
    f8(x) = (V2 = v2(x[1], x[2], x[3], x[4], x[5], x[6], x[8]);
             coefs(x[1], V2, x[2], x[3], x[4], x[5], x[6], x[7], kip(x[1], V2, x[2], x[5], x[7]), x[8], Keq))
    println("Cleland formula, Keq fixed, V2 from H1 (9 free): dimension ", jac_rank(f9, [logu(rng) for _ in 1:9]))
    println("Cleland formula, Keq fixed, V2 from H1 and Kip from H2 (8 free): dimension ",
            jac_rank(f8, [logu(rng) for _ in 1:8]))
end
