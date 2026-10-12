# ABOUTME: Card C4 check of the printed ping-pong K-form of specs/7: v/E_t = V1(ab − pq)/(a + b + ab + p + q
# ABOUTME: + (Ka/KdA)ap + (Kb/KdB)bq + (V1/V2)pq), exact at random points, against the joint-solved law.
include(joinpath(@__DIR__, "kform_lib.jl"))

m = byid("specs/7")
L = P4.law(m)
D, Nm = P4.cmap(L, L.den), P4.cmap(L, L.num)
c(xs...) = D[P4.cmono(L, xs)]
n(xs...) = Nm[P4.cmono(L, xs)]
rng = Xoshiro(13)
for _ in 1:3
    x = randpoint(L, rng)
    e(p) = evalq(p, x)
    k(s) = x[L.vars.index[s]]
    N = e(c(:A)) * e(c(:B)) / e(c(:A, :B))
    Ka, Kb = e(c(:B)) / e(c(:A, :B)), e(c(:A)) / e(c(:A, :B))
    KP, KQ = N / e(c(:P)), N / e(c(:Q))
    V1, V2 = e(n(:A, :B)) / e(c(:A, :B)), -e(n(:P, :Q)) / e(c(:P, :Q))
    KdA, KdB = k(:k_EA_to_E_A) / k(:k_E_A_to_EA), k(:k_FB_to_F_B) / k(:k_F_B_to_FB)
    keq = x[L.keq]
    a, b, p, q = (x[L.conc[s]] for s in (:A, :B, :P, :Q))
    a, b, p, q = a / Ka, b / Kb, p / KP, q / KQ
    v_form = V1 * (a * b - p * q) / (a + b + a * b + p + q + (Ka / KdA) * a * p + (Kb / KdB) * b * q + (V1 / V2) * p * q)
    v_law = evalq(L.num, x) / evalq(L.den, x)
    kcat = k(:k_EA_to_F_P) * k(:k_FB_to_E_Q) / (k(:k_EA_to_F_P) + k(:k_FB_to_E_Q))
    println("v_form == v_law: ", v_form == v_law, "; V1 == k2k4/(k2+k4): ", V1 == kcat,
            "; KP*KQ == Keq*Ka*Kb: ", KP * KQ == keq * Ka * Kb)
end
