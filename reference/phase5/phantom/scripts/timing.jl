# ABOUTME: Times the two exact L-phantom tests on sampled one-subunit variants: the rescaling
# ABOUTME: certificate and the modular rank of the variant and of its L → 0 law (warm, best of 3).
include(joinpath(@__DIR__, "plib.jl"))
include(joinpath(SP, "phase3", "census", "census_common.jl"))
const DATA = ARGS[1]
rows = read_rows(joinpath(DATA, "recs_bibi.tsv"))
facts = Dict{Int, String}()
for l in eachline(joinpath(DATA, "facts_bibi.tsv"))
    x = split(l, '\t')
    facts[parse(Int, x[1])] = l
end
pick = [r for r in rows[2:end] if r[5] == "1" && r[9] == "ok"]
pf = [r for r in pick if occursin("pf_bind=1", facts[parse(Int, r[1])])]
rng = Xoshiro(5)
sample = [pf[randperm(rng, length(pf))[1:200]]; pick[randperm(rng, length(pick))[1:200]]]
best(f) = minimum(@elapsed(f()) for _ in 1:3)
ts = Dict(k => Float64[] for k in ("law", "cert", "rank2", "rank1"))
for (i, r) in enumerate(sample)
    sm = parse_neutral(r[10], [:A, :B], [:P, :Q])
    pm = parse_neutral(r[10], [:A, :B], [:P, :Q]; plain = true)
    A = S.allo_law(sm)
    i == 1 && (certificate(A); rank_mod(blaw(A), colsof(blaw(A), collect(S.fitted_params(sm))), rng))
    push!(ts["law"], best(() -> S.allo_law(sm)))
    push!(ts["cert"], best(() -> certificate(A)))
    fv, fp = collect(S.fitted_params(sm)), collect(S.fitted_params(pm))
    LP = S.law(pm)
    push!(ts["rank2"], best(() -> (BV = blaw(A); BP = blaw(LP);
                                  rank_mod(BV, colsof(BV, fv), rng) == rank_mod(BP, colsof(BP, fp), rng))))
    push!(ts["rank1"], best(() -> (BV = blaw(A); BP = blaw(LP); Q1 = denlcm(BV); Q2 = denlcm(BP);
        rank_fp(jacobian(F1, BV, colsof(BV, fv), Q1, length(fv) + 6, rng)) ==
        rank_fp(jacobian(F1, BP, colsof(BP, fp), Q2, length(fp) + 6, rng)))))
end
q(xs, p) = sort(xs)[clamp(round(Int, p * length(xs)), 1, length(xs))]
for k in ("law", "cert", "rank1", "rank2")
    @printf("%-6s per variant (400 bi-bi variants: 200 pre-filter passes, 200 random): median %.3f ms, p95 %.3f ms, max %.3f ms\n",
            k, 1e3 * q(ts[k], 0.5), 1e3 * q(ts[k], 0.95), 1e3 * maximum(ts[k]))
end
println("law = slice allo_law; cert = certificate(A) on a built law; rank1 = one draw mod 2^61-1 of the variant and of the L→0 law; rank2 = both primes")
