# ABOUTME: Re-ranks the candidates whose complex-step rank was ambiguous or disagreed with the census,
# ABOUTME: in Complex{BigFloat} (256 bits) at three fresh draws with concentrations in 10^[-1,1].
include(joinpath(@__DIR__, "vlib.jl"))
setprecision(BigFloat, 256)

"""Rank of the log-Jacobian in Complex{BigFloat} complex step (h = 1e-60): Gaussian elimination
with complete pivoting, pivots above 1e-40 of the first; the smallest kept and largest dropped."""
function big_rank(pl, rng; extra = 12)
    n = length(pl.names)
    θ = [big(lu10(rng, -1, 1)) for _ in 1:n]
    keq = big(lu10(rng, -1, 1))
    h = big(10)^-60
    J = zeros(BigFloat, n + extra, n)
    for r in 1:(n + extra)
        c = Dict{Symbol, BigFloat}(x => big(lu10(rng, -1, 1)) for x in pl.mets)
        v = real(prate(pl, θ, keq, c))
        for j in 1:n
            θc = complex.(θ)
            θc[j] += im * θ[j] * h
            J[r, j] = imag(prate(pl, θc, keq, c)) / h / v
        end
    end
    A = copy(J); m = size(A, 1); piv = BigFloat[]
    for k in 1:min(m, n)
        _, idx = findmax(abs.(A[k:m, k:n]))
        i, j = idx[1] + k - 1, idx[2] + k - 1
        A[[k, i], :] = A[[i, k], :]; A[:, [k, j]] = A[:, [j, k]]
        push!(piv, abs(A[k, k]))
        iszero(A[k, k]) && break
        for r in (k + 1):m
            A[r, k:n] .-= (A[r, k] / A[k, k]) .* A[k, k:n]
        end
    end
    rk = count(q -> q > big(10)^-40 * piv[1], piv)
    (rank = rk, kept = Float64(piv[rk] / piv[1]), dropped = rk < length(piv) ? Float64(piv[rk + 1] / piv[1]) : 0.0)
end

const RECHECK = Dict(
    "uni_uni_allo_2reg" => Int[], "uni_bi_rxn" => Int[], "bi_bi_rxn" => Int[], "bi_bi_pp_rxn" => Int[],
    "bi_bi_rxn.r2_100000" => Int[], "bi_bi_pp_rxn.r2_1200" => Int[], "uni_uni_allo_reg" => Int[],
    "uni_uni_reg_and_inhibitor" => Int[], "uni_uni_with_reg" => Int[], "uni_uni_rxn" => Int[],
    "uni_uni_allo" => Int[])
# the ambiguous and disagreeing candidates of v4_P.log, by run
let cur = ""
    for ln in eachline(joinpath(VER, "logs", "v4_P.log"))
        startswith(ln, "== ") && (cur = ln[4:end]; continue)
        m = match(r"^\s+\[(P|unflagged)\] #(\d+) .*(DISAGREE|AMBIGUOUS)", ln)
        m === nothing || push!(RECHECK[cur], parse(Int, m.captures[2]))
    end
end
rng = Xoshiro(4242)
for (name, gs) in sort!(collect(RECHECK))
    isempty(gs) && continue
    d, mechs = loadrun(name)
    for g in unique(gs)
        r = d.recs[g]
        pl = plaw(mechs[g])
        rs = [big_rank(pl, rng) for _ in 1:3]
        @printf("%s #%d np %d census rank %d | BigFloat ranks %s | kept %s dropped %s %s\n", name, g,
                length(pl.names), r[:rank_pkg], join((x.rank for x in rs), ","),
                join((@sprintf("%.1e", x.kept) for x in rs), ","),
                join((@sprintf("%.1e", x.dropped) for x in rs), ","),
                all(x -> x.rank == r[:rank_pkg], rs) ? "agrees" : "DISAGREES")
        flush(stdout)
    end
end
