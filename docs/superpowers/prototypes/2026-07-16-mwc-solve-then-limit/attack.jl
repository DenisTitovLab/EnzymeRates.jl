# ABOUTME: Adversarial attack on "solve then limit". Tests the two regimes the build report
# ABOUTME: left unbuilt: metabolite-bearing-D (claim B) and the P:OnlyA-with-active-I-catalysis trap.
using LinearAlgebra, Random, Printf

# ── oracle solver (copied verbatim from test/allosteric_ground_truth.jl) ───────
function mwc_ground_truth_flux(species, edges, cat_edges, Etot)
    n = length(species); idx = Dict(s => i for (i, s) in enumerate(species))
    M = zeros(n, n)
    for (a, b, r) in edges
        M[idx[b], idx[a]] += r; M[idx[a], idx[a]] -= r
    end
    A = copy(M); A[1, :] .= 1.0
    rhs = zeros(n); rhs[1] = Etot
    c = A \ rhs
    (flux = sum(kf * c[idx[r]] - kr * c[idx[p]] for (r, p, kf, kr) in cat_edges),
     conc = Dict(s => c[idx[s]] for s in species))
end
flux(species, edges, cat_edges, Etot=1.0) = mwc_ground_truth_flux(species, edges, cat_edges, Etot).flux

# ─────────────────────────────────────────────────────────────────────────────
# ATTACK 1 (decisive claim B): metabolite-bearing-D ordered bi-uni, S:OnlyA.
# The oracle's free-flip-only ground truth for the LDH i-state case.
# Does the PLAIN combine (N_A + L*N_I)/(D_A + L*D_I) — with each state normalized
# to its own free-enzyme weight — match it, or is cross-weighting needed?
# ─────────────────────────────────────────────────────────────────────────────
# oracle: two-conformation free-flip-only, S+B<=>P, S:OnlyA(SS), B:EqualAI(RE),
# cat:OnlyA(SS), P:EqualAI(RE).  D_A carries metabolite B; D_I = 1 (inactive binds nothing).
function metab_dfree_onlyA_flux(kon, koff, KB, KP, k; L, Keq, S, B, P, FAST=1e7)
    kr = k * kon * KP / (koff * KB * Keq)
    species = [:E_A, :ES_A, :ESB_A, :EP_A, :E_I, :EP_I]
    edges = [
        (:E_A, :ES_A, kon * S), (:ES_A, :E_A, koff),
        (:ES_A, :ESB_A, FAST * B / KB), (:ESB_A, :ES_A, FAST),
        (:E_A, :EP_A, FAST * P / KP), (:EP_A, :E_A, FAST),
        (:E_I, :EP_I, FAST * P / KP), (:EP_I, :E_I, FAST),
        (:E_A, :E_I, FAST * L), (:E_I, :E_A, FAST),
        (:EP_A, :EP_I, FAST * L), (:EP_I, :EP_A, FAST),
        (:ESB_A, :EP_A, k), (:EP_A, :ESB_A, kr),
    ]
    flux(species, edges, [(:ESB_A, :EP_A, k, kr)], 1.0)
end

# single-conformation active-only network -> extract free-enzyme conc & flux,
# from which the free-flip combine needs N_A = v_A/[E_free], D_A = 1/[E_free].
function active_only_metabD(kon, koff, KB, KP, k, Keq, S, B, P; FAST=1e7)
    kr = k * kon * KP / (koff * KB * Keq)
    species = [:E, :ES, :ESB, :EP]
    edges = [
        (:E, :ES, kon * S), (:ES, :E, koff),
        (:ES, :ESB, FAST * B / KB), (:ESB, :ES, FAST),
        (:E, :EP, FAST * P / KP), (:EP, :E, FAST),
        (:ESB, :EP, k), (:EP, :ESB, kr),
    ]
    r = mwc_ground_truth_flux(species, edges, [(:ESB, :EP, k, kr)], 1.0)
    (Efree = r.conc[:E], vA = r.flux)          # E_total=1 so vA = N_A_abs/D_A_abs
end

println("="^92)
println("ATTACK 1 (claim B decisive): metabolite-bearing-D  plain combine vs cross-weight vs oracle")
println("="^92)
@printf("%-6s %-16s %-16s %-16s %-10s %-10s\n",
        "draw","oracle GT","plain combine","cross-weight","relerr_pl","relerr_cw")
rng = MersenneTwister(424242); maxpl = 0.0; maxcw = 0.0
for d in 1:10
    kon=0.5+2rand(rng); koff=0.5+2rand(rng); KB=0.5+2rand(rng); KP=0.5+2rand(rng)
    k=0.5+2rand(rng); L=0.5+rand(rng); Keq=2.0+2rand(rng)
    S=0.5+2rand(rng); B=0.5+2rand(rng); P=0.5+2rand(rng)

    vgt = metab_dfree_onlyA_flux(kon,koff,KB,KP,k; L=L,Keq=Keq,S=S,B=B,P=P)

    ao = active_only_metabD(kon,koff,KB,KP,k,Keq,S,B,P)
    D_A = 1/ao.Efree                     # active King-Altman den normalized to free E = 1
    N_A = ao.vA/ao.Efree                 # active numerator normalized likewise
    D_I = 1 + P/KP                       # inactive: {E_I, EP_I}, free E = 1, P:EqualAI RE
    N_I = 0.0                            # cat:OnlyA
    v_plain = (N_A + L*N_I)/(D_A + L*D_I)

    # what UNnormalized cross-weight would give if you (wrongly) used raw Q_A with free!=1
    # Q_A = D_A_abs = D_A*Efree_weight... here demonstrate the plain combine is the correct one.
    # cross-weight variant: den = D_I*Q_A + L*D_A_unit*Q_I  with Q_A un-normalized (free-weight w_E)
    wE = ao.Efree            # NOT the weight; placeholder to show a different normalization breaks
    v_cross = (N_A*D_I)/(D_A*D_I + L*D_A*D_I)   # = N_A/(D_A(1+L)) — a wrong normalization

    relpl = abs(v_plain-vgt)/max(abs(vgt),1e-12)
    relcw = abs(v_cross-vgt)/max(abs(vgt),1e-12)
    global maxpl = max(maxpl,relpl); global maxcw = max(maxcw,relcw)
    @printf("%-6d %-16.8g %-16.8g %-16.8g %-10.2e %-10.2e\n", d, vgt, v_plain, v_cross, relpl, relcw)
end
println("-"^92)
@printf("plain-combine max relerr = %.2e    cross-weight-variant max relerr = %.2e\n", maxpl, maxcw)
println()

# ─────────────────────────────────────────────────────────────────────────────
# ATTACK 2 (the trap + guard C): P:OnlyA with ACTIVE I-state catalysis.
# S:EqualAI, cat:NonequalAI, P:OnlyA.  Guard requires cat==:EqualAI so it stays
# SILENT here.  TRUE physics: ES_I is populated (S:EqualAI) and NonequalAI catalysis
# runs in I -> produces EP_I as a DEAD LEAF (P can't rebind, OnlyA). EP_I holds mass
# -> contributes to D_I.  The prototype's v_new (and its own gt_flux) DROP the EP_I
# term because they gate I-catalysis on EP_I-binding.  Build the true GT and compare.
# ─────────────────────────────────────────────────────────────────────────────
# prototype's derivation (copied logic from solve_then_limit.jl v_new, uni-uni)
function v_new(tags, p, concs)
    Stag, cattag, Ptag = tags
    S, P = concs.S, concs.P
    K_I_S = Stag=='E' ? p.K_A_S : p.K_I_S    # placeholder; overwritten below
    K_I_S = Stag == :EqualAI ? p.K_A_S : p.K_I_S
    k_I   = cattag == :EqualAI ? p.k_A : p.k_I
    K_I_P = Ptag == :EqualAI ? p.K_A_P : p.K_I_P
    k_A_r = p.k_A*p.K_A_P/(p.K_A_S*p.Keq)
    N_A = p.k_A*S/p.K_A_S - k_A_r*P/p.K_A_P
    D_A = 1 + S/p.K_A_S + P/p.K_A_P
    ES_I = Stag != :OnlyA
    EP_I = Ptag != :OnlyA
    Icat = cattag != :OnlyA && ES_I && EP_I
    D_I = 1.0
    ES_I && (D_I += S/K_I_S)
    EP_I && (D_I += P/K_I_P)
    if Icat
        k_I_r = k_I*K_I_P/(K_I_S*p.Keq)
        N_I = k_I*S/K_I_S - k_I_r*P/K_I_P
    else
        N_I = 0.0
    end
    p.E_total*(N_A + p.L*N_I)/(D_A + p.L*D_I)
end

# TRUE free-flip-only ground truth: I-state catalysis runs whenever ES_I exists and
# cat!=:OnlyA, PRODUCING EP_I even if P doesn't bind in I. EP_I is then a leaf whose
# only exits are reverse catalysis (and P release iff P binds in I).
function true_gt(tags, p, concs; FAST=1e7)
    Stag, cattag, Ptag = tags
    S, P = concs.S, concs.P
    K_I_S = Stag == :EqualAI ? p.K_A_S : p.K_I_S
    k_I   = cattag == :EqualAI ? p.k_A : p.k_I
    K_I_P = Ptag == :EqualAI ? p.K_A_P : p.K_I_P
    k_A_r = p.k_A*p.K_A_P/(p.K_A_S*p.Keq)
    ES_I = Stag != :OnlyA
    EP_I_binds = Ptag != :OnlyA
    Icat = cattag != :OnlyA && ES_I           # catalysis runs if substrate reachable in I
    EP_I_present = EP_I_binds || Icat          # EP_I exists if P binds OR catalysis makes it

    species = [:E_A, :ES_A, :EP_A, :E_I]
    ES_I && push!(species, :ES_I)
    EP_I_present && push!(species, :EP_I)
    edges = Tuple{Symbol,Symbol,Float64}[
        (:E_A,:ES_A, FAST*S/p.K_A_S), (:ES_A,:E_A, FAST),
        (:E_A,:EP_A, FAST*P/p.K_A_P), (:EP_A,:E_A, FAST),
        (:ES_A,:EP_A, p.k_A), (:EP_A,:ES_A, k_A_r),
        (:E_A,:E_I, FAST*p.L), (:E_I,:E_A, FAST),
    ]
    cat_edges = [(:ES_A,:EP_A,p.k_A,k_A_r)]
    ES_I && push!(edges, (:E_I,:ES_I, FAST*S/K_I_S), (:ES_I,:E_I, FAST))
    EP_I_binds && push!(edges, (:E_I,:EP_I, FAST*P/K_I_P), (:EP_I,:E_I, FAST))
    if Icat
        k_I_r = k_I*K_I_P/(K_I_S*p.Keq)
        push!(edges, (:ES_I,:EP_I, k_I), (:EP_I,:ES_I, k_I_r))
        push!(cat_edges, (:ES_I,:EP_I,k_I,k_I_r))
    end
    flux(species, edges, cat_edges, p.E_total)
end

function guard(tags)
    Stag, cattag, Ptag = tags
    cattag == :EqualAI || return false
    (Ptag==:EqualAI && Stag in (:OnlyA,:NonequalAI)) && return true
    (Stag==:EqualAI && Ptag in (:OnlyA,:NonequalAI)) && return true
    false
end

draw(rng) = (K_A_S=0.5+2rand(rng), k_A=0.5+2rand(rng), K_A_P=0.5+2rand(rng),
             K_I_S=0.5+2rand(rng), k_I=0.5+2rand(rng), K_I_P=0.5+2rand(rng),
             L=0.5+rand(rng), Keq=2.0+2rand(rng), E_total=1.0)

println("="^92)
println("ATTACK 2 (trap + guard C): P:OnlyA with active I-state catalysis")
println("="^92)
trap_combos = [
    (:EqualAI, :NonequalAI, :OnlyA),   # S binds in I, NonequalAI cat runs in I, P:OnlyA -> EP_I trap
    (:NonequalAI,:NonequalAI,:OnlyA),
    (:EqualAI, :EqualAI, :OnlyA),       # for reference (guard SHOULD fire here)
]
@printf("%-34s %-14s %-14s %-10s %-8s\n","(S,cat,P)","v_new","true GT","relerr","guard")
for tags in trap_combos
    if guard(tags)
        @printf("%-34s %-14s %-14s %-10s %-8s\n", string(tags), "-", "-", "-", "FIRED")
        continue
    end
    rng = MersenneTwister(hash(tags)); mr = 0.0; vn0=0.0; vg0=0.0
    for _ in 1:10
        p = draw(rng); concs = (S=0.5+2rand(rng), P=0.5+2rand(rng))
        vn = v_new(tags,p,concs); vg = true_gt(tags,p,concs)
        mr = max(mr, abs(vn-vg)/max(abs(vg),1e-12)); vn0=vn; vg0=vg
    end
    @printf("%-34s %-14.6g %-14.6g %-10.2e %-8s\n", string(tags), vn0, vg0, mr, "silent")
end
println()
println("If relerr is large and guard is 'silent', solve-then-limit returns a WRONG answer")
println("with no guard firing = a genuine break (the EP_I catalytic-trap term is dropped).")
