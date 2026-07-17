# ABOUTME: Prototype of Denis's "solve then limit" MWC allosteric derivation (n=1 uni-uni S<=>P).
# ABOUTME: Validates claims (A) correctness vs mass-action ground truth, (B) no cross-weighting, (C) collapse guard.
using LinearAlgebra, Random, Printf

# ── Ground-truth solver (copied from test/allosteric_ground_truth.jl) ─────────
function mwc_ground_truth_flux(species, edges, cat_edges, Etot)
    n = length(species); idx = Dict(s => i for (i, s) in enumerate(species))
    M = zeros(n, n)
    for (a, b, r) in edges
        M[idx[b], idx[a]] += r
        M[idx[a], idx[a]] -= r
    end
    A = copy(M); A[1, :] .= 1.0
    rhs = zeros(n); rhs[1] = Etot
    c = A \ rhs
    sum(kf * c[idx[r]] - kr * c[idx[p]] for (r, p, kf, kr) in cat_edges)
end

# ── Tag semantics ─────────────────────────────────────────────────────────────
# tags = (S, cat, P), each in (:EqualAI, :NonequalAI, :OnlyA).
# For a state we need K_I_S, k_I, K_I_P. EqualAI => equals the A value.
# We pass the full param set; helper resolves I-values.
function ivalues(tags, p)
    Stag, cattag, Ptag = tags
    K_I_S = Stag   == :EqualAI ? p.K_A_S : p.K_I_S
    k_I   = cattag == :EqualAI ? p.k_A   : p.k_I
    K_I_P = Ptag   == :EqualAI ? p.K_A_P : p.K_I_P
    (K_I_S=K_I_S, k_I=k_I, K_I_P=K_I_P)
end

# ── The GUARD (claim C): does the thermodynamic constraint solve FORCE an
#    OnlyA/NonequalAI param to collapse onto its A-counterpart? ────────────────
# The load-bearing physics: an :EqualAI CATALYTIC step shares BOTH microscopic
# rate constants k_f AND k_r (the catalytic step is conformation-independent).
# Each state obeys its own Haldane  k_r = k_f * K_P / (K_S * Keq).  Imposing
# k_I_r = k_A_r together with k_I_f = k_A_f gives the constraint (*)
#         K_I_P / K_I_S  =  K_A_P / K_A_S .
# (*) is active ONLY when cat==:EqualAI.  Given (*):
#   - P EqualAI (K_I_P=K_A_P)  => forces K_I_S = K_A_S.
#   - S EqualAI (K_I_S=K_A_S)  => forces K_I_P = K_A_P.
# A forced finite equality contradicts an :OnlyA (K_I->inf) or a free :NonequalAI
# intent -> collapse -> guard fires.
function guard_collapse(tags)
    Stag, cattag, Ptag = tags
    cattag == :EqualAI || return (fired=false, reason="")
    if Ptag == :EqualAI && Stag in (:OnlyA, :NonequalAI)
        return (fired=true,
            reason="EqualAI cat + EqualAI P forces K_I_S = K_A_S, contradicting $(Stag) S")
    end
    if Stag == :EqualAI && Ptag in (:OnlyA, :NonequalAI)
        return (fired=true,
            reason="EqualAI cat + EqualAI S forces K_I_P = K_A_P, contradicting $(Ptag) P")
    end
    (fired=false, reason="")
end

# ── The DERIVATION: solve-then-limit, n=1 combine  v = Et*(N_A + L*N_I)/(D_A + L*D_I) ─
# Each state's single-state King-Altman (RE S-binding, SS catalysis, RE P-binding):
#   N = k_f*S/K_S - k_r*P/K_P ,   D = 1 + S/K_S + P/K_P ,   k_r = k_f*K_P/(K_S*Keq).
# :OnlyA is applied as a LIMIT AFTER the (trivial, per-state) Haldane solve:
#   OnlyA S-binding: K_I_S->inf  => ES_I term drops (and k_I_r->0).
#   OnlyA P-binding: K_I_P->inf  => EP_I term drops (and reverse via EP_I ->0).
#   OnlyA catalysis: k_I->0      => I-state catalytic flux ->0.
# I-state catalysis exists only if the I cycle is complete (ES_I and EP_I present).
function v_new(tags, p, concs)
    g = guard_collapse(tags)
    g.fired && error("GUARD fired: $(g.reason)")
    Stag, cattag, Ptag = tags
    S, P = concs.S, concs.P
    iv = ivalues(tags, p)

    k_A_r = p.k_A * p.K_A_P / (p.K_A_S * p.Keq)
    N_A = p.k_A * S / p.K_A_S - k_A_r * P / p.K_A_P
    D_A = 1 + S / p.K_A_S + P / p.K_A_P

    ES_I = Stag != :OnlyA          # inactive S-binding present?
    EP_I = Ptag != :OnlyA          # inactive P-binding present?
    Icat = cattag != :OnlyA && ES_I && EP_I

    D_I = 1.0
    ES_I && (D_I += S / iv.K_I_S)
    EP_I && (D_I += P / iv.K_I_P)

    if Icat
        k_I_r = iv.k_I * iv.K_I_P / (iv.K_I_S * p.Keq)
        N_I = iv.k_I * S / iv.K_I_S - k_I_r * P / iv.K_I_P
    else
        N_I = 0.0
    end

    p.E_total * (N_A + p.L * N_I) / (D_A + p.L * D_I)
end

# ── General free-flip-only (formulation 1) ground-truth network for any tag combo ─
# Only the free enzyme flips (E_A<->E_I, ratio L); each conformation runs its own
# King-Altman cycle. Matches the n=1 combine exactly (see derivation above).
function gt_flux(tags, p, concs; FAST=1e7)
    Stag, cattag, Ptag = tags
    S, P = concs.S, concs.P
    iv = ivalues(tags, p)
    k_A_r = p.k_A * p.K_A_P / (p.K_A_S * p.Keq)

    ES_I = Stag != :OnlyA
    EP_I = Ptag != :OnlyA
    Icat = cattag != :OnlyA && ES_I && EP_I

    species = [:E_A, :ES_A, :EP_A, :E_I]
    ES_I && push!(species, :ES_I)
    EP_I && push!(species, :EP_I)

    edges = Tuple{Symbol,Symbol,Float64}[
        (:E_A, :ES_A, FAST * S / p.K_A_S), (:ES_A, :E_A, FAST),
        (:E_A, :EP_A, FAST * P / p.K_A_P), (:EP_A, :E_A, FAST),
        (:ES_A, :EP_A, p.k_A), (:EP_A, :ES_A, k_A_r),
        (:E_A, :E_I, FAST * p.L), (:E_I, :E_A, FAST),
    ]
    cat_edges = [(:ES_A, :EP_A, p.k_A, k_A_r)]
    if ES_I
        push!(edges, (:E_I, :ES_I, FAST * S / iv.K_I_S), (:ES_I, :E_I, FAST))
    end
    if EP_I
        push!(edges, (:E_I, :EP_I, FAST * P / iv.K_I_P), (:EP_I, :E_I, FAST))
    end
    if Icat
        k_I_r = iv.k_I * iv.K_I_P / (iv.K_I_S * p.Keq)
        push!(edges, (:ES_I, :EP_I, iv.k_I), (:EP_I, :ES_I, k_I_r))
        push!(cat_edges, (:ES_I, :EP_I, iv.k_I, k_I_r))
    end
    mwc_ground_truth_flux(species, edges, cat_edges, p.E_total)
end

# ── Random parameter draw (I-params present for all groups; used only when tag needs them) ─
function draw(rng)
    (K_A_S = 0.5 + 2rand(rng), k_A = 0.5 + 2rand(rng), K_A_P = 0.5 + 2rand(rng),
     K_I_S = 0.5 + 2rand(rng), k_I = 0.5 + 2rand(rng), K_I_P = 0.5 + 2rand(rng),
     L = 0.5 + rand(rng), Keq = 2.0 + 2rand(rng), E_total = 1.0)
end

# ── Test harness ──────────────────────────────────────────────────────────────
combos = [
    ("Family A  (S:OnlyA, cat:OnlyA, P:EqualAI)",      (:OnlyA, :OnlyA, :EqualAI)),
    ("TR        (S:EqualAI, cat:OnlyA, P:NonequalAI)", (:EqualAI, :OnlyA, :NonequalAI)),
    ("All-EqualAI",                                     (:EqualAI, :EqualAI, :EqualAI)),
    ("NonequalAI cat (S:EqualAI, cat:NonequalAI, P:EqualAI)", (:EqualAI, :NonequalAI, :EqualAI)),
    ("INCONSISTENT (S:OnlyA, cat:EqualAI, P:EqualAI)", (:OnlyA, :EqualAI, :EqualAI)),
]

println("="^92)
println("n=1 uni-uni  S<=>P   :  solve-then-limit  vs  mass-action ground truth")
println("="^92)
@printf("%-52s %-12s %-14s %-10s\n", "combo", "max relerr", "v@equilibrium", "guard")
println("-"^92)

results = []
for (name, tags) in combos
    g = guard_collapse(tags)
    if g.fired
        @printf("%-52s %-12s %-14s %-10s\n", name, "-", "-", "FIRED")
        push!(results, (name, tags, nothing, nothing, g))
        continue
    end
    rng = MersenneTwister(hash(tags))
    maxrel = 0.0
    for _ in 1:10
        p = draw(rng)
        concs = (S = 0.5 + 2rand(rng), P = 0.5 + 2rand(rng))
        vn = v_new(tags, p, concs)
        vg = gt_flux(tags, p, concs)
        rel = abs(vn - vg) / max(abs(vg), 1e-12)
        maxrel = max(maxrel, rel)
    end
    # equilibrium check: P/S = Keq  => v ~ 0
    p = draw(rng); Seq = 1.3
    ceq = (S = Seq, P = p.Keq * Seq)
    veq = v_new(tags, p, ceq)
    @printf("%-52s %-12.2e %-14.2e %-10s\n", name, maxrel, veq, "ok(no fire)")
    push!(results, (name, tags, maxrel, veq, g))
end
println("-"^92)

# ── Claim (B): Family A derived equation carries NO cross-weighting / d_free factor ─
println()
println("="^92)
println("CLAIM (B): Family A derived symbolic equation")
println("="^92)
println("""
Family A = (S:OnlyA, cat:OnlyA, P:EqualAI).  After the per-state Haldane solve and
the OnlyA limits (K_I_S -> inf, k_I -> 0):

  N_A = k_A*S/K_A_S - k_A_r*P/K_A_P        k_A_r = k_A*K_A_P/(K_A_S*Keq)
  D_A = 1 + S/K_A_S + P/K_A_P
  N_I = 0                                  (S:OnlyA drops ES_I, cat:OnlyA drops k_I)
  D_I = 1 + P/K_A_P                        (only EqualAI P survives in I)

  v = E_total * (N_A + L*0) / (D_A + L*D_I)
    = E_total * (k_A*S/K_A_S - k_A_r*P/K_A_P) / ( (1 + S/K_A_S + P/K_A_P) + L*(1 + P/K_A_P) )

The denominator is a PLAIN sum  D_A + L*D_I .  There is NO d_free_A/d_free_I
cross-weighting factor, no  D_I*Q_A + L*D_A*Q_I  product.  D_A and D_I are each
normalized to their own free-enzyme term (= 1), so the naive combine is exact.
""")

# programmatic confirmation: the denominator equals D_A + L*D_I with no cross product
let tags = (:OnlyA, :OnlyA, :EqualAI)
    rng = MersenneTwister(1)
    p = draw(rng); concs = (S=1.1, P=0.7)
    k_A_r = p.k_A*p.K_A_P/(p.K_A_S*p.Keq)
    N_A = p.k_A*concs.S/p.K_A_S - k_A_r*concs.P/p.K_A_P
    D_A = 1 + concs.S/p.K_A_S + concs.P/p.K_A_P
    D_I = 1 + concs.P/p.K_A_P
    v_manual_plain = p.E_total*N_A/(D_A + p.L*D_I)            # NO cross-weight
    println("  numeric check:  v_new = ", v_new(tags,p,concs),
            "   plain-combine = ", v_manual_plain,
            "   match=", isapprox(v_new(tags,p,concs), v_manual_plain; rtol=1e-12))
    # A hypothetical cross-weighted denominator would be D_I*D_A + L*D_A*D_I = D_A*D_I*(1+L)
    v_crossweight = p.E_total*N_A*D_I/(D_A*D_I + p.L*D_A*D_I)
    println("  cross-weighted form would give: ", v_crossweight,
            "   (differs from ground truth: ", !isapprox(v_crossweight, gt_flux(tags,p,concs); rtol=1e-4), ")")
    println("  ground truth = ", gt_flux(tags,p,concs))
end

# ── Claim (C): show the forced collapse explicitly for the inconsistent combo ─
println()
println("="^92)
println("CLAIM (C): forced collapse for INCONSISTENT combo (S:OnlyA, cat:EqualAI, P:EqualAI)")
println("="^92)
println("""
EqualAI catalysis shares k_f AND k_r.  Per-state Haldane:
  A:  k_A_r = k_A * K_A_P / (K_A_S * Keq)
  I:  k_I_r = k_I * K_I_P / (K_I_S * Keq)
EqualAI cat => k_I = k_A and (shared reverse) k_I_r = k_A_r.  Divide the two Haldanes:
  1 = (K_A_P/K_A_S) / (K_I_P/K_I_S)   =>   K_I_S = K_A_S * (K_I_P / K_A_P).
P is EqualAI => K_I_P = K_A_P  =>  K_I_S = K_A_S  (FORCED, finite).
But S is :OnlyA and wants K_I_S -> inf.  Contradiction -> collapse -> guard fires.
""")
let tags = (:OnlyA, :EqualAI, :EqualAI)
    g = guard_collapse(tags)
    println("  guard_collapse(", tags, ") => fired=", g.fired)
    println("  reason: ", g.reason)
    # demonstrate numerically: impose the shared-k_r constraint and watch K_I_S forced to K_A_S
    K_A_S, K_A_P, k_A, Keq = 1.3, 0.9, 2.1, 3.0
    K_I_P = K_A_P                          # P EqualAI
    k_A_r = k_A*K_A_P/(K_A_S*Keq)
    # solve k_I_r = k_A_r with k_I = k_A for K_I_S:
    K_I_S_forced = k_A*K_I_P/(k_A_r*Keq)   # from k_A_r = k_A*K_I_P/(K_I_S*Keq)
    println("  numeric: shared-k_r constraint forces K_I_S = ", K_I_S_forced,
            "  (== K_A_S = ", K_A_S, "?  ", isapprox(K_I_S_forced, K_A_S; rtol=1e-12), ")")
    try
        v_new(tags, (K_A_S=K_A_S,k_A=k_A,K_A_P=K_A_P,K_I_S=1e9,k_I=k_A,K_I_P=K_I_P,L=0.7,Keq=Keq,E_total=1.0), (S=1.1,P=0.7))
    catch e
        println("  v_new correctly refuses to evaluate: ", sprint(showerror, e))
    end
end

# sanity: my free-flip-only gt for Family A and all-EqualAI must match the
# canonical hand-written networks from the oracle (which add redundant EP/ES flips)
println()
println("="^92)
println("SANITY: free-flip-only gt matches oracle's hand-written networks")
println("="^92)
function uni_onlyA_flux(KA, KP, k; L, Keq, S, P, FAST=1e7)
    kr = k * KP / (Keq * KA)
    species = [:E_A, :ES_A, :EP_A, :E_I, :EP_I]
    edges = [
        (:E_A, :ES_A, FAST*S/KA), (:ES_A, :E_A, FAST),
        (:E_A, :EP_A, FAST*P/KP), (:EP_A, :E_A, FAST),
        (:E_I, :EP_I, FAST*P/KP), (:EP_I, :E_I, FAST),
        (:E_A, :E_I, FAST*L), (:E_I, :E_A, FAST),
        (:EP_A, :EP_I, FAST*L), (:EP_I, :EP_A, FAST),
        (:ES_A, :EP_A, k), (:EP_A, :ES_A, kr),
    ]
    mwc_ground_truth_flux(species, edges, [(:ES_A, :EP_A, k, kr)], 1.0)
end
let
    p = (K_A_S=1.3,k_A=2.1,K_A_P=0.9,K_I_S=9.9,k_I=9.9,K_I_P=0.9,L=0.7,Keq=3.0,E_total=1.0)
    concs = (S=1.1,P=0.6)
    mine = gt_flux((:OnlyA,:OnlyA,:EqualAI), p, concs)
    oracle = uni_onlyA_flux(1.3, 0.9, 2.1; L=0.7, Keq=3.0, S=1.1, P=0.6)
    println("  Family A:  my free-flip gt = ", mine, "   oracle uni_onlyA_flux = ", oracle,
            "   match=", isapprox(mine, oracle; rtol=1e-5))
end
