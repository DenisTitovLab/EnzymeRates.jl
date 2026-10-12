# ABOUTME: Card C4 extras: the ping-pong law printed after dividing by N = c(A)c(B)/c(AB), identities of
# ABOUTME: its γ's (Cleland ratios), sum-reference group examples, and a cross-check against census3.
include(joinpath(@__DIR__, "kform_lib.jl"))

rt(L, a, b) = b == one(PI) ? P4.fmt_factored(a, L.vars) :
              P4.fmt_factored(a, L.vars) * " / (" * P4.fmt_factored(b, L.vars) * ")"
cm(L, mo) = isempty(mo) ? "1" : P4.fmt_mono(mo, L.vars)

# 1. specs/7 printed under N = c(A)c(B)/c(AB): each term's coefficient d_m/N and n_m/N
m = byid("specs/7")
L = P4.law(m)
D, Nm = P4.cmap(L, L.den), P4.cmap(L, L.num)
c(xs...) = get(D, P4.cmono(L, xs), PI())
n(xs...) = get(Nm, P4.cmono(L, xs), PI())
Npoly = (c(:A) * c(:B), c(:A, :B))          # N as (numerator, denominator)
println("=== specs/7 divided by N = c(A)c(B)/c(AB)")
for (lab, M) in (("den", D), ("num", Nm))
    for mo in sort(collect(keys(M)); by = mo -> monokey(L, mo))
        a, b = P4.reduce_ratio(M[mo] * Npoly[2], Npoly[1])
        println("  ", lab, ": ", cm(L, mo), " × ", rt(L, a, b))
    end
end
# γ identities at random points (exact): γ(PQ) = V1/V2 with V1 = n(AB)/c(AB), V2 = -n(PQ)/c(PQ)
rng = Xoshiro(5)
for _ in 1:3
    x = randpoint(L, rng)
    e(p) = evalq(p, x)
    γPQ = e(c(:P, :Q)) * e(c(:A)) * e(c(:B)) / (e(c(:P)) * e(c(:Q)) * e(c(:A, :B)))
    V1 = e(n(:A, :B)) / e(c(:A, :B)); V2 = -e(n(:P, :Q)) / e(c(:P, :Q))
    kd_B = x[L.vars.index[:k_FB_to_F_B]] / x[L.vars.index[:k_F_B_to_FB]]
    Kb = e(c(:A)) / e(c(:A, :B))
    γBQ = e(c(:B, :Q)) * e(c(:A)) / (e(c(:Q)) * e(c(:A, :B)))
    println("  γ(PQ) == V1/V2: ", γPQ == V1 / V2, "; γ(BQ) == Kb/Kd(B) with Kd(B) = k_FB_to_F_B/k_F_B_to_FB: ", γBQ == Kb / kd_B)
end

# 2. groups of constant-term laws whose d_0 is a sum: smallest examples per shape of K_X = d_0/d_X
shape(p) = nterms(p) == 1 ? "mono" : "sum"
seen = Dict{Tuple{String, String}, Int}()
for pop in Pop.populations(), mm in Pop.load_population(pop)
    mm.allosteric && continue
    LL = P4.law(mm)
    DD = P4.cmap(LL, LL.den)
    haskey(DD, Mono()) || continue
    d0 = DD[Mono()]
    nterms(d0) > 1 || continue
    for mo in sort(collect(keys(DD)); by = mo -> monokey(LL, mo))
        P4.degree(mo) == 1 || continue
        a, b = P4.reduce_ratio(d0, DD[mo])      # K_X = a/b
        k = (shape(a), shape(b))
        get(seen, k, 0) >= 2 && continue
        nterms(a) + nterms(b) <= 8 || continue
        seen[k] = get(seen, k, 0) + 1
        println("K_", cm(LL, mo), " shape ", k, " [", mm.id, "]: ", rt(LL, a, b))
    end
end

# 3. census3 cross-check: laws with a constant term and a linear term for every metabolite of D
agree = 0; tot = 0; pure3 = 0
for pop in Pop.populations(), mm in Pop.load_population(pop)
    mm.allosteric && continue
    LL = P4.law(mm)
    DD = P4.cmap(LL, LL.den)
    haskey(DD, Mono()) || continue
    used = Set(v for mo in keys(DD) for (v, _) in mo)
    all(v -> haskey(DD, mono(v => 1)), used) || continue
    global tot += 1
    k = kform(LL, DD)
    pure = all(==(:one), k.cls)
    pure && (global pure3 += 1)
end
println("laws with a constant term and every linear term: ", tot, "; pure K-form (all γ = 1): ", pure3)
