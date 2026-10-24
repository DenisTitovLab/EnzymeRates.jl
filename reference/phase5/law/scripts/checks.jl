# ABOUTME: Engine questions on the slice's law (variant U): coefficient bounds against the cap (does
# ABOUTME: BigInt ever run), the cap's determinant in checked Int (Q10), RE turnover (Q9), class 5 (Q14).
include(joinpath(@__DIR__, "lib.jl"))
const S = load_variant("U")
using LinearAlgebra
const TINY = get(ENV, "TINY", "0") == "1"
t0 = time()

rows = NamedTuple[]
for r in census_rows()
    try
        push!(rows, (set = r.set, id = r.id, m = slice_mech(S, r.net)))
    catch
    end
end
TINY && (rows = rows[1:15:end])
"""The law masks of `m`: every edge, and for an allosteric mechanism the inactive one."""
masks(m) = S.is_allosteric(m) && any(S.inactive_mask(m)) ?
    [(trues(length(m.edges)), false), (S.inactive_mask(m), true)] :
    [(trues(length(m.edges)), false)]
println("julia ", VERSION, "; mechanisms: ", length(rows), " (",
        join(("$g $(count(r -> group(r.set) == g, rows))" for g in unique(group.(getfield.(rows, :set)))), ", "), ")")
flush(stdout)

# ── Q10: the cap's determinant in checked Int ───────────────────────────────────────────
"""Exact determinant of `A` by fraction-free (Bareiss) elimination in checked `Int`; an
`OverflowError` when an intermediate product leaves `Int`'s range."""
function det_bareiss(A::Matrix{Int})
    M, n = copy(A), size(A, 1)
    n == 0 && return 1
    sgn, prev = 1, 1
    for k in 1:(n - 1)
        if M[k, k] == 0
            p = findfirst(i -> M[i, k] != 0, (k + 1):n)
            p === nothing && return 0
            p += k
            for j in 1:n
                M[k, j], M[p, j] = M[p, j], M[k, j]
            end
            sgn = -sgn
        end
        for i in (k + 1):n, j in (k + 1):n
            M[i, j] = div(Base.checked_sub(Base.checked_mul(M[i, j], M[k, k]),
                                           Base.checked_mul(M[i, k], M[k, j])), prev)
        end
        prev = M[k, k]
    end
    sgn * M[n, n]
end

"""The segment graph's reduced Laplacian of `m`'s edges in `mask` and its node count."""
function laplacian(m, mask)
    node, G = S.segment_nodes(m, mask)
    lap = zeros(Int, G, G)
    for e in m.edges[mask]
        a, b = node[e.a], node[e.b]
        (e.re || a == b) && continue
        lap[a, a] += 1; lap[b, b] += 1; lap[a, b] -= 1; lap[b, a] -= 1
    end
    lap[2:end, 2:end], G
end

"""`complexity` in checked `Int`: the BigInt determinant only after an overflow, its value
saturated at `typemax(Int)`, which is over every cap."""
function complexity_int(m, mask)
    lap, G = laplacian(m, mask)
    G <= 1 && return G
    d = try
        det_bareiss(lap)
    catch err
        err isa OverflowError || rethrow()
        min(det(BigInt.(lap)), big(typemax(Int)))
    end
    Base.checked_mul(G, Int(d))
end

neq, ncalls, tB, tI = 0, 0, Float64[], Float64[]
for r in rows, (mask, inactive) in masks(r.m)
    cB = S.complexity(r.m, mask)
    cI = complexity_int(r.m, mask)
    global ncalls += 1
    global neq += cB == cI
    push!(tB, best(() -> S.complexity(r.m, mask)))
    push!(tI, best(() -> complexity_int(r.m, mask)))
end
println("Q10 det: checked-Int Bareiss equals the BigInt det in ", neq, " of ", ncalls,
        " (mechanism, conformation) pairs")
println("Q10 time per call µs: BigInt det ", stats(tB), " | Int Bareiss ", stats(tI))
for cap in (337, 5000)
    over = Dict{String, Int}()
    nallo = 0
    for r in rows
        S.is_allosteric(r.m) || continue
        nallo += 1
        ms = masks(r.m)
        cA = complexity_int(r.m, ms[1][1])
        cI = length(ms) == 2 ? complexity_int(r.m, ms[2][1]) : 0
        cA <= cap < Base.checked_add(cA, cI) && (over[group(r.set)] = get(over, group(r.set), 0) + 1)
    end
    println("Q10 cap $cap: allosteric mechanisms with active <= cap < active + inactive: ",
            isempty(over) ? "none" : over, " of ", nallo)
end
rng = Xoshiro(7)
nsyn, okay = 0, 0
for n in 2:14, trial in 1:30
    μ = rand(rng, 1:3)
    A = zeros(Int, n, n)
    for a in 1:n, b in (a + 1):n
        rand(rng) < 0.5 || continue
        k = rand(rng, 1:μ)
        A[a, a] += k; A[b, b] += k; A[a, b] -= k; A[b, a] -= k
    end
    L = A[2:end, 2:end]
    global nsyn += 1
    global okay += try det_bareiss(L) == det(BigInt.(L)) catch; false end
end
println("Q10 synthetic random multigraph Laplacians (2-14 nodes, multiplicity 1-3): Bareiss exact in ",
        okay, " of ", nsyn)
for (n, μ) in ((10, 1), (16, 2), (20, 3), (30, 3))
    K = fill(-μ, n, n)
    foreach(i -> K[i, i] = μ * (n - 1), 1:n)
    L = K[2:end, 2:end]
    exact = det(BigInt.(L))
    r = try
        string(det_bareiss(L))
    catch err
        string(typeof(err))
    end
    println("Q10 complete multigraph K_$n x$μ: trees ", exact, " (", exact > typemax(Int) ?
            "over" : "within", " Int64); checked Bareiss -> ", r)
end
flush(stdout)
function churn(f, n)
    k = 0
    while k < n
        for r in rows, (mask, _) in masks(r.m)
            f(r.m, mask)
            (k += 1) >= n && break
        end
    end
end
for (label, f) in (("Int Bareiss", complexity_int), ("BigInt det", S.complexity))
    churn(f, 200)
    GC.gc(true); GC.gc(true)
    lb0, g0 = Base.gc_live_bytes(), Base.gc_num()
    n = TINY ? 2_000 : 20_000
    t = @elapsed churn(f, n)
    g1 = Base.gc_num()
    GC.gc(true); GC.gc(true)
    lb1 = Base.gc_live_bytes()
    @printf("Q10 live bytes, %s: %d calls in %.2f s; after a full collection %+.2f KiB per call \
             (%+.1f MiB); %d collections, %d full\n", label, n, t, (lb1 - lb0) / 1024 / n,
            (lb1 - lb0) / 2^20, g1.pause - g0.pause, g1.full_sweep - g0.full_sweep)
end
flush(stdout)

# ── coefficient bounds: can checked Int overflow under the cap? ─────────────────────────
# Claim: |bare D|_1 = F·τ (F touched forms, τ multigraph spanning trees); |bare num|_1 <=
# 2τ·Σ|w|; substitution and the conc shift never raise |·|_1; dividing by H keeps every
# coefficient of Q at most |num|_1. So every coefficient <= max(F·τ, 2τ·Σ|w|).
nlaw, nbig, nid, nwithin, worst, worst_ratio, worst_id = 0, 0, 0, 0, big(0), 0.0, ""
maxF = 0
single = Tuple{BigInt, Int}[]                     # single-term state factors fI·DA, fA·DI: (coef, n)
for r in rows, (mask, inactive) in masks(r.m)
    m = r.m
    try
        S.complexity(m, mask) > S.MAX_TERMS && continue
        b = S.bare_law(m, mask)
        node, G = S.segment_nodes(m, mask)
        τ = G <= 1 ? 1 : Int(S.complexity(m, mask) ÷ G)
        F = count(>(0), node)
        den1 = sum((abs(BigInt(c)) for c in values(b.den.terms)); init = big(0))
        num1 = sum((abs(BigInt(c)) for c in values(b.num.terms)); init = big(0))
        global nid += den1 == F * τ
        L = S.law(m, S.solve(m), mask; inactive)
        global nbig += L isa S.Law{BigInt}
        mx = maximum(maxcoef, (b.num, b.den, L.Q, L.D))
        bound = max(big(F) * τ, num1)
        global nwithin += mx <= bound && num1 <= 2 * big(τ) * 30 * count(k -> mask[k] && !m.edges[k].re, eachindex(m.edges))
        global nlaw += 1
        global maxF = max(maxF, F)
        if mx > worst
            global worst, worst_id = mx, r.id
        end
        global worst_ratio = max(worst_ratio, Float64(mx) / Float64(bound))
    catch e
    end
end
println("bounds: ", nlaw, " laws (both conformations); BigInt needed in ", nbig,
        "; |bare den|_1 = F x trees in ", nid, " of ", nlaw,
        "; every coefficient <= max(F*trees, |bare num|_1) and |bare num|_1 <= 2*trees*30*SS-edges in ",
        nwithin, " of ", nlaw)
println("bounds: largest coefficient anywhere (bare num, bare den, Q, D) ", worst, " (", worst_id,
        "); largest ratio to the bound ", round(worst_ratio; digits = 4), "; most touched forms ", maxF,
        "; Int64 max ", typemax(Int))
for r in rows
    S.is_allosteric(r.m) || continue
    try
        A = S.rate_law(r.m)
        DI = A.inactive === nothing ? one(typeof(A.H)) : A.inactive.D
        for p in (A.fI * A.active.D, A.fA * DI)
            S.nterms(p) == 1 && push!(single, (abs(BigInt(only(values(p.terms)))), A.n))
        end
    catch
    end
end
println("codegen powers: single-term state factors (fI*DA or fA*DI) ", length(single),
        "; largest coefficient^n ", isempty(single) ? 0 : maximum(c^n for (c, n) in single))
flush(stdout)

# ── Q9 and Q14 on the census and on phase 4's hand cases ────────────────────────────────
hand = NamedTuple[]
for r in hand_rows()
    try
        push!(hand, (set = "hand", id = r.id, m = slice_mech(S, r.net)))
    catch e
        println("hand ", r.id, ": not built (", first(sprint(showerror, e), 80), ")")
    end
end
nre, t9, t5, fail5, nlaws = Dict{String, Int}(), Float64[], Float64[], Dict{String, Int}(), 0
tl = Float64[]
for r in [rows; hand]
    m = r.m
    flagged = [re_turnover_edge(S, m, mask) for (mask, _) in masks(m)]
    any(!=(0), flagged) && (nre[group(r.set)] = get(nre, group(r.set), 0) + 1)
    push!(t9, best(() -> re_turnover_edge(S, m, trues(length(m.edges)))))
    L = try
        sol = S.solve(m)
        push!(tl, best(() -> S.rate_law(m, sol); reps = 2))
        S.rate_law(m, sol)
    catch e
        r.set == "hand" && println(rpad(r.id, 44), " RE-turnover edge ", join(flagged, "/"),
                                   "; no law: ", first(sprint(showerror, e), 90))
        nothing
    end
    L === nothing && continue
    global nlaws += 1
    v = class5(m, L)
    push!(t5, best(() -> class5(m, L)))
    isempty(v) || (fail5[group(r.set)] = get(fail5, group(r.set), 0) + 1)
    if r.set == "hand" || startswith(r.set, "spec:")
        nq = hasproperty(L, :active) ? S.nterms(L.active.Q) : S.nterms(L.Q)
        (r.set == "hand" || any(!=(0), flagged) || !isempty(v)) &&
            println(rpad(r.set == "hand" ? r.id : r.set, 44), " RE-turnover edge ",
                    join(flagged, "/"), "; class 5 ", isempty(v) ? "pass" : "FAIL at " * v,
                    "; Q terms ", nq)
    end
end
println("Q9: mechanisms with an RE-turnover edge, by set: ", isempty(nre) ? "none" : nre)
println("Q9: check time µs ", stats(t9), " | law time µs ", stats(tl))
println("Q14: ", nlaws, " laws; class-5 failures by set: ", isempty(fail5) ? "none" : fail5)
println("Q14: class-5 time µs ", stats(t5))
println("finished in ", round(Int, time() - t0), " s")
