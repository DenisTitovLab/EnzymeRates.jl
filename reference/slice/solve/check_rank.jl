# ABOUTME: Solve-stage check (e), report only: the identifiable rank of the rate in the log
# ABOUTME: of the fitted constants, by central differences on the two-conformation oracle.
include(joinpath(@__DIR__, "lib.jl"))
include(joinpath(PROTO, "p5a", "pkgvals.jl"))

const Λ = MWC.big1024(10)^parse(Int, get(ENV, "LAMBDA_EXP", "40"))
const ONLY = split(get(ENV, "ONLY", ""), ",", keepempty = false)
const NRANDOM = parse(Int, get(ENV, "NRANDOM", "60"))
const MAXPICK = parse(Int, get(ENV, "MAXPICK", "100000"))

"""Rank of `J` by Gaussian elimination with complete pivoting (P5b's `rank_gap`): pivots
above `tol` times the first; the smallest kept and largest dropped pivot ratios."""
function rank_gap(J; tol = big(10)^-35)
    A = copy(J); m, n = size(A); piv = BigFloat[]
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
    rk = count(q -> q > tol * piv[1], piv)
    rk, Float64(piv[rk] / piv[1]), rk < length(piv) ? Float64(piv[rk + 1] / piv[1]) : 0.0
end

"""The two-conformation rate of record `r` at fitted values `θ` (in `fitted_params` order)."""
function rate(r, θ, keq, concs)
    m, p, c = r.c.mech, r.p, r.c
    pars = (; (S.fitted_params(m) .=> θ)..., Keq = keq, E_total = big(1))
    allp = S.all_parameters(m, pars)
    mask = S.inactive_mask(m)
    kA = S.pop_constants(c, edge_values(m, allp))
    kI = S.pop_constants(c, edge_values(m, allp; inactive = true, mask))
    act = [mask[c.edge_of[j]] for j in eachindex(p.steps)]
    keep = Set(f for (j, s) in enumerate(p.steps) if act[j] for f in (s.from, s.to))
    push!(keep, MWC.free_form(p))
    rn = Dict(S.regulator_names(m))
    regA, regI = Dict{Tuple{Int, Symbol}, Any}(), Dict{Tuple{Int, Symbol}, Any}()
    for (si, site) in enumerate(p.sites), (lig, state) in zip(site.ligands, site.states)
        ns = [getfield(allp, n) for n in rn[lig]]
        regA[(si, lig)] = state === :OnlyI ? nothing : ns[1]
        regI[(si, lig)] = state === :OnlyA ? nothing : ns[end]
    end
    MWC.v_factorized(p, kA, kI, [s.re for s in p.steps], act, keep, allp.L, regA, regI,
                     concs; n = p.cat_multiplicity, Λ = Λ)
end

records = [r for r in load_records(("specs", "handwritten", "allo_uniuni", "allo_bibi",
                                    "expand_bibi", "allo_uniuni_inh", "allo_unibi",
                                    "allo_bibi_pp", "allo_terbi"))
           if S.is_allosteric(r.c.mech)]
rng = Xoshiro(20261013)
pick = [r for r in records if r.set in ("specs", "handwritten", "allo_uniuni") ||
        length(S.fitted_params(r.c.mech)) != length(r.p.fitted_params)]
rest = [r for r in records if !(r in pick)]
append!(pick, rest[randperm(rng, length(rest))[1:NRANDOM]])
pick = pick[1:min(end, MAXPICK)]
isempty(ONLY) || (pick = [r for r in records if r.p.id in ONLY])
@printf("identifiable rank on %d of %d allosteric mechanisms\n", length(pick), length(records))

tally = Dict{Tuple{String, Int, Int}, Int}()     # (kind, n, count - rank) => mechanisms
gaps = Float64[]
t0 = time()
for r in pick
    θ0 = [Oracle.logu(rng) for _ in S.fitted_params(r.c.mech)]
    keq = Oracle.logu(rng)
    np = length(θ0)
    h = big(10)^-25
    J = zeros(BigFloat, 2np + 4, np)
    for i in axes(J, 1)
        concs = Dict(x => Oracle.logu(rng) for x in r.p.metabolites)
        for k in 1:np
            θp, θm = copy(θ0), copy(θ0)
            θp[k] *= exp(h); θm[k] *= exp(-h)
            J[i, k] = (rate(r, θp, keq, concs) - rate(r, θm, keq, concs)) / (2h)
        end
    end
    rk, kept, dropped = rank_gap(J)
    kind = length(S.fitted_params(r.c.mech)) != length(r.p.fitted_params) ?
           "count differs from the package" : r.set in ("specs", "handwritten",
           "allo_uniuni") ? "specs, handwritten, allo_uniuni" : "random sample"
    key = (kind, r.p.cat_multiplicity, np - rk)
    tally[key] = get(tally, key, 0) + 1
    push!(gaps, kept)
    kept < 1e-20 && @printf("  weak kept pivot: %s %s (n = %d, %d fitted; rank %d; kept %.1e, \
                             dropped %.1e)\n", r.set, r.p.id, r.p.cat_multiplicity, np, rk,
                            kept, dropped)
    np - rk > 0 && @printf("  deficit %d: %s %s (n = %d, %d fitted; kept %.1e, dropped %.1e)\n",
                           np - rk, r.set, r.p.id, r.p.cat_multiplicity, np, kept, dropped)
end
@printf("time %.0f s; smallest kept pivot ratio %.1e\n", time() - t0, minimum(gaps))
for k in sort!(collect(keys(tally)))
    @printf("  %-32s n = %d, fitted count - rank = %d: %d mechanisms\n", k..., tally[k])
end
