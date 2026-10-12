# ABOUTME: Recomputes the uni_bi round-2 census row end to end with independent code: identify's
# ABOUTME: pipeline through the public expand_mechanisms, then W, D, P, N, Z by this verifier's methods.
include(joinpath(@__DIR__, "vlib.jl"))
include(joinpath(@__DIR__, "moracle.jl"))

const RXN = uni_bi_rxn
rng = Xoshiro(717171)

# ── the pipeline (identify's order) ────────────────────────────────────────────────────
seen = Set{UInt64}()
memo = Dict{UInt64, Int}()
cands = NamedTuple[]          # every fresh candidate: mech, round, parent index, status, eqh, fitted
function batch!(items, round)
    nseen = 0
    first = length(cands) + 1
    for (m, par) in items
        h = hash(m)
        h in seen && (nseen += 1; continue)
        push!(seen, h)
        st, eqh, np = :ok, UInt64(0), -1
        try
            if ER._eq_complexity(m) > 337
                st = :cx
            else
                em = ER.compile_mechanism(m)
                np = length(ER.fitted_params(em))
                np > 20 ? (st = :np) : (eqh = ER._rate_eq_dedup_key(ER.rate_equation_string(em)))
            end
        catch e
            st = :err
        end
        push!(cands, (mech = m, round = round, parent = par, status = st, eqh = eqh, np = np,
                      fitted = Ref(false)))
    end
    for c in cands[first:end]
        c.status === :ok || continue
        haskey(memo, c.eqh) && continue
        memo[c.eqh] = length(memo) + 1
        c.fitted[] = true
    end
    b = cands[first:end]
    @printf("round %d: generated %d seen %d fresh %d cx %d np>20 %d err %d ok %d eqdup %d fitted %d\n",
            round, length(items), nseen, length(b), count(c -> c.status === :cx, b),
            count(c -> c.status === :np, b), count(c -> c.status === :err, b),
            count(c -> c.status === :ok, b), count(c -> c.status === :ok && !c.fitted[], b),
            count(c -> c.fitted[], b))
    first:length(cands)
end
function expand_round(parents, round)
    items = Tuple{Any, Int}[]
    firstp = Dict{Any, Int}()
    for i in parents
        for c in ER.expand_mechanisms(Union{ER.Mechanism, ER.AllostericMechanism}[cands[i].mech], RXN)
            haskey(firstp, c) && continue
            firstp[c] = i
            push!(items, (c, i))
        end
    end
    batch!(items, round)
end

seeds = unique!(collect(ER.init_mechanisms(RXN)))
base, _ = ER._base_tier(seeds, RXN)
r0 = batch!([(m, 0) for m in base], 0)
r1 = expand_round([i for i in r0 if cands[i].status === :ok], 1)
r2 = expand_round([i for i in r1 if cands[i].status === :ok], 2)

# ── match to the census's records ──────────────────────────────────────────────────────
d, mechs = loadrun("uni_bi_rxn")
cidx = Dict(m => g for (g, m) in enumerate(mechs))
mine2 = [i for i in r2 if cands[i].fitted[]]
cen2 = [r[:gid] for r in d.recs if r[:round] == 2 && r[:status] === :ok && r[:fitted]]
matched = [get(cidx, cands[i].mech, 0) for i in mine2]
println("round-2 fitted: mine $(length(mine2)), census $(length(cen2)); mine found in census $(count(>(0), matched)); ",
        "same order $(matched == cen2); same parent $(all(((i, g),) -> g > 0 && d.recs[g][:parent] == get(cidx, cands[cands[i].parent].mech, -1), zip(mine2, matched)))")

# ── classes, by this verifier's methods ────────────────────────────────────────────────
fitted_all = [i for i in eachindex(cands) if cands[i].fitted[]]
pls = Dict(i => plaw(cands[i].mech) for i in fitted_all)
mss = Dict(i => mstruct(cands[i].mech) for i in fitted_all)
rk = Dict(i => cs_rank(pls[i], rng) for i in fitted_all)

"""Integer asymptotic log-slopes of |v| in each metabolite (to ∞ and to 0) with every other
metabolite positive, with the products at zero and with the substrates at zero: equal for two
laws of one family at generic constants."""
function fingerprint(pl, ms)
    θ = [lu10(rng, -1, 1) for _ in pl.names]
    keq = lu10(rng, -1, 1)
    fp = Int[]
    for zero in (Symbol[], ms.prods, ms.subs)
        base = Dict{Symbol, Float64}(x => (x in zero ? 0.0 : lu10(rng, -1, 1)) for x in pl.mets)
        for x in pl.mets
            x in zero && continue
            v(t) = (c = copy(base); c[x] = t; abs(real(prate(pl, θ, keq, c))))
            hi = log10(v(1e10) / v(1e7)) / 3
            lo = log10(v(1e-10) / v(1e-7)) / -3
            push!(fp, isfinite(hi) ? round(Int, hi) : 99, isfinite(lo) ? round(Int, lo) : 99)
        end
    end
    fp
end
fps = Dict(i => (rk[i].rank, fingerprint(pls[i], mss[i])) for i in fitted_all)

"""Same family both ways, after a cheap one-way screen (one draw, three starts)."""
function same_law(a, b)
    keq = lu10(rng, -1, 1)
    fit_into(a, [lu10(rng, -1, 1) for _ in a.names], b, keq, rng; starts = 3, iters = 150) < 1e-6 ||
        return false
    ab, ba = same_family(a, b, rng)
    ab < 1e-8 && ba < 1e-8
end
D = Dict{Int, Int}()
ntest = 0
for (k, i) in enumerate(fitted_all)
    cands[i].round == 2 || continue
    for j in fitted_all[1:k-1]
        fps[j] == fps[i] || continue
        global ntest += 1
        if same_law(pls[i], pls[j])
            D[i] = j
            break
        end
    end
end
Z = Dict{Int, Bool}()
for i in mine2
    p = cands[i].parent
    pls_p = get(pls, p, nothing)
    pls_p === nothing && (pls_p = plaw(cands[p].mech))
    Z[i] = same_law(pls[i], pls_p)
end
W = Dict{Int, Vector{String}}()
for i in mine2
    ft = face_test(mss[i], rng)
    W[i] = [my_admit(mss[i]); ["face:" * f.face for f in ft.faces if f.fails]]
end
P = Dict(i => rk[i].rank < length(pls[i].names) for i in mine2)
N = Dict{Int, Bool}()
for i in mine2
    gr = growth(pls[i], mss[i], rng)
    N[i] = max(gr[:fwd], gr[:rev]) > 0.05
end
anyc = [i for i in mine2 if !isempty(W[i]) || haskey(D, i) || P[i] || N[i] || Z[i]]
n = length(mine2)
@printf("\nMY ROW uni_bi round 2: fitted %d | W %d | D(same law as an earlier fit) %d | P %d | N %d | Z %d | any %d | OK %d   (LM tests %d)\n",
        n, count(i -> !isempty(W[i]), mine2), count(i -> haskey(D, i), mine2), count(i -> P[i], mine2),
        count(i -> N[i], mine2), count(i -> Z[i], mine2), length(anyc), n - length(anyc), ntest)
# census flags of the same candidates
byg = Dict(r[:gid] => r for r in d.recs)
first_key = Dict{String, Int}()
cD = Set{Int}()
for r in d.recs
    (r[:status] === :ok && r[:fitted]) || continue
    haskey(first_key, r[:key]) ? push!(cD, r[:gid]) : (first_key[r[:key]] = r[:gid])
end
cP(g) = byg[g][:rank_pkg] < byg[g][:np]
cZ(g) = byg[g][:parent] != 0 && byg[byg[g][:parent]][:key] == byg[g][:key]
pairs_ = collect(zip(mine2, matched))
@printf("CENSUS ROW (same candidates): D %d | P %d | Z %d\n", count(((i, g),) -> g in cD, pairs_),
        count(((i, g),) -> cP(g), pairs_), count(((i, g),) -> cZ(g), pairs_))
@printf("agreement: D %d/%d, P %d/%d, Z %d/%d\n",
        count(((i, g),) -> haskey(D, i) == (g in cD), pairs_), n,
        count(((i, g),) -> P[i] == cP(g), pairs_), n, count(((i, g),) -> Z[i] == cZ(g), pairs_), n)
for (i, g) in pairs_
    (haskey(D, i) == (g in cD) && P[i] == cP(g) && Z[i] == cZ(g)) && continue
    j = get(D, i, 0)
    @printf("  #%d: mine D=%s P=%s Z=%s | census D=%s P=%s Z=%s | np %d rank mine %d census %d\n      %s\n",
            g, j == 0 ? "-" : "#$(get(cidx, cands[j].mech, 0))", P[i], Z[i], g in cD, cP(g), cZ(g),
            byg[g][:np], rk[i].rank, byg[g][:rank_pkg], byg[g][:text])
end
# same-law-as-earlier beyond D: how many are P (phantom) and which move
extra = [i for i in mine2 if haskey(D, i) && !(get(cidx, cands[i].mech, 0) in cD)]
println("same law as an earlier fit but not D: $(length(extra)); of which P: $(count(i -> P[i], extra))")
