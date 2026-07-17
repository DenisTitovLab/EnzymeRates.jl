# ABOUTME: n=2 spot-check of the solve-then-limit ^n combine for Family A vs an explicit
# ABOUTME: 2-protomer concerted-MWC mass-action ground truth. Tests limit-commutes-with-^n.
using LinearAlgebra, Random, Printf

function mwc_ground_truth_flux(species, edges, cat_edges, Etot)
    n = length(species); idx = Dict(s => i for (i, s) in enumerate(species))
    M = zeros(n, n)
    for (a, b, r) in edges
        M[idx[b], idx[a]] += r; M[idx[a], idx[a]] -= r
    end
    A = copy(M); A[1, :] .= 1.0
    rhs = zeros(n); rhs[1] = Etot
    c = A \ rhs
    sum(kf * c[idx[r]] - kr * c[idx[p]] for (r, p, kf, kr) in cat_edges)
end

# ── Explicit 2-protomer concerted MWC, Family A (S:OnlyA, cat:OnlyA, P:EqualAI) ──
# Conformation R (both protomers active) or T (both inactive). Ordered protomers.
# R occupancy per protomer in {0=empty, S, P}; T occupancy in {0, P} (no S in T).
# Concerted flip R<->T allowed only when neither protomer bears S (ratio L; P is
# EqualAI so contributes ratio 1). Catalysis (protomer S->P) only in R.
function twoprotomer_familyA(; K_A_S, K_A_P, k_A, L, Keq, S, P, FAST=1e7)
    k_A_r = k_A * K_A_P / (K_A_S * Keq)
    Rocc = [(a,b) for a in (:0,:S,:P) for b in (:0,:S,:P)]
    Tocc = [(a,b) for a in (:0,:P) for b in (:0,:P)]
    Rname(o) = Symbol("R_", o[1], "_", o[2])
    Tname(o) = Symbol("T_", o[1], "_", o[2])
    species = vcat([Rname(o) for o in Rocc], [Tname(o) for o in Tocc])
    edges = Tuple{Symbol,Symbol,Float64}[]
    cat_edges = Tuple{Symbol,Symbol,Float64,Float64}[]
    setocc(o,i,v) = i==1 ? (v,o[2]) : (o[1],v)
    # R-state per-protomer bindings + catalysis
    for o in Rocc, i in 1:2
        cur = o[i]
        if cur == :0
            # empty -> S
            edges = push!(edges, (Rname(o), Rname(setocc(o,i,:S)), FAST*S/K_A_S))
            edges = push!(edges, (Rname(setocc(o,i,:S)), Rname(o), FAST))
            # empty -> P
            edges = push!(edges, (Rname(o), Rname(setocc(o,i,:P)), FAST*P/K_A_P))
            edges = push!(edges, (Rname(setocc(o,i,:P)), Rname(o), FAST))
        end
        if cur == :S
            # catalysis S -> P (this protomer)
            oP = setocc(o,i,:P)
            edges = push!(edges, (Rname(o), Rname(oP), k_A))
            edges = push!(edges, (Rname(oP), Rname(o), k_A_r))
            cat_edges = push!(cat_edges, (Rname(o), Rname(oP), k_A, k_A_r))
        end
    end
    # T-state per-protomer P binding
    for o in Tocc, i in 1:2
        if o[i] == :0
            edges = push!(edges, (Tname(o), Tname(setocc(o,i,:P)), FAST*P/K_A_P))
            edges = push!(edges, (Tname(setocc(o,i,:P)), Tname(o), FAST))
        end
    end
    # concerted flip R<->T for S-free occupancies, ratio L
    for o in Tocc
        edges = push!(edges, (Rname(o), Tname(o), FAST*L))
        edges = push!(edges, (Tname(o), Rname(o), FAST))
    end
    mwc_ground_truth_flux(species, edges, cat_edges, 1.0)
end

# ── solve-then-limit ^n combine, n=2, Family A (N_I=0) ─────────────────────────
# v_per_protomer = Et*(N_A*D_A^(n-1) + L*N_I*D_I^(n-1)) / (D_A^n + L*D_I^n)
# Family A: N_I = 0, D_I = 1 + P/K_A_P.  Molecule has n=2 protomers => x2.
function v_new_n2_familyA(; K_A_S, K_A_P, k_A, L, Keq, S, P, n=2)
    k_A_r = k_A * K_A_P / (K_A_S * Keq)
    N_A = k_A*S/K_A_S - k_A_r*P/K_A_P
    D_A = 1 + S/K_A_S + P/K_A_P
    D_I = 1 + P/K_A_P              # OnlyA S dropped, OnlyA cat -> N_I=0
    per_protomer = (N_A*D_A^(n-1)) / (D_A^n + L*D_I^n)
    n * per_protomer              # molecule flux = n protomers
end

println("="^92)
println("n=2 spot-check: Family A  solve-then-limit ^n combine  vs  2-protomer concerted MWC")
println("="^92)
@printf("%-10s %-18s %-18s %-12s\n", "draw", "combine(x n)", "2-protomer GT", "relerr")
rng = MersenneTwister(20260716)
maxrel = 0.0
for d in 1:8
    K_A_S = 0.5+2rand(rng); K_A_P = 0.5+2rand(rng); k_A = 0.5+2rand(rng)
    L = 0.5+rand(rng); Keq = 2.0+2rand(rng); S = 0.5+2rand(rng); P = 0.5+2rand(rng)
    vc = v_new_n2_familyA(; K_A_S,K_A_P,k_A,L,Keq,S,P)
    vg = twoprotomer_familyA(; K_A_S,K_A_P,k_A,L,Keq,S,P)
    rel = abs(vc-vg)/max(abs(vg),1e-12); global maxrel = max(maxrel, rel)
    @printf("%-10d %-18.6g %-18.6g %-12.2e\n", d, vc, vg, rel)
end
println("-"^92)
@printf("max relerr over 8 draws: %.2e\n", maxrel)

# equilibrium check
let K_A_S=1.3,K_A_P=0.9,k_A=2.1,L=0.7,Keq=3.0,S=1.4
    P = Keq*S
    println("v @ equilibrium (P/S=Keq): combine=", v_new_n2_familyA(;K_A_S,K_A_P,k_A,L,Keq,S,P),
            "   GT=", twoprotomer_familyA(;K_A_S,K_A_P,k_A,L,Keq,S,P))
end

# L=0 sanity: both must reduce to 2 * single-active-site rate
let K_A_S=1.3,K_A_P=0.9,k_A=2.1,Keq=3.0,S=1.1,P=0.6
    k_A_r=k_A*K_A_P/(K_A_S*Keq)
    site = (k_A*S/K_A_S - k_A_r*P/K_A_P)/(1+S/K_A_S+P/K_A_P)
    println("L=0: 2*single-site=", 2site, "  GT=", twoprotomer_familyA(;K_A_S,K_A_P,k_A,L=0.0,Keq,S,P),
            "  combine=", v_new_n2_familyA(;K_A_S,K_A_P,k_A,L=0.0,Keq,S,P))
end
