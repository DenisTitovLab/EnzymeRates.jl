# ABOUTME: Independent confirmation that the trap-family break is real physics, not a true_gt artifact.
# ABOUTME: Builds the FULL 6-form free-flip-1 network from scratch, cross-checks, shows trapped EP_I mass.
using LinearAlgebra, Random, Printf

# solver returns flux + full concentration dict
function solvenet(species, edges, cat_edges, Etot)
    n = length(species); idx = Dict(s => i for (i, s) in enumerate(species))
    M = zeros(n, n)
    for (a,b,r) in edges; M[idx[b],idx[a]] += r; M[idx[a],idx[a]] -= r; end
    A = copy(M); A[1,:] .= 1.0; rhs = zeros(n); rhs[1] = Etot
    c = A \ rhs
    (flux = sum(kf*c[idx[r]] - kr*c[idx[p]] for (r,p,kf,kr) in cat_edges),
     conc = Dict(s => c[idx[s]] for s in species))
end

# FULL 6-form network built from scratch with NO decision about which forms to drop.
# formulation 1: only E_A<->E_I flips. Every tag-permitted binding/catalysis edge added.
# A tag-permitted edge exists iff the ligand binds / catalysis runs in that conformation.
function full6(tags, p, concs; FAST=1e7)
    Stag, cattag, Ptag = tags; S, P = concs.S, concs.P
    K_I_S = Stag==:EqualAI ? p.K_A_S : p.K_I_S
    k_I   = cattag==:EqualAI ? p.k_A : p.k_I
    K_I_P = Ptag==:EqualAI ? p.K_A_P : p.K_I_P
    k_A_r = p.k_A*p.K_A_P/(p.K_A_S*p.Keq)
    k_I_r = k_I*K_I_P/(K_I_S*p.Keq)
    species = [:E_A,:ES_A,:EP_A,:E_I,:ES_I,:EP_I]      # ALL six forms, always
    edges = Tuple{Symbol,Symbol,Float64}[
        # active conformation: S binds, P binds, catalysis (all always present in A)
        (:E_A,:ES_A, FAST*S/p.K_A_S), (:ES_A,:E_A, FAST),
        (:E_A,:EP_A, FAST*P/p.K_A_P), (:EP_A,:E_A, FAST),
        (:ES_A,:EP_A, p.k_A), (:EP_A,:ES_A, k_A_r),
        # free-enzyme flip (formulation 1: ONLY this flip)
        (:E_A,:E_I, FAST*p.L), (:E_I,:E_A, FAST),
    ]
    cat_edges = Tuple{Symbol,Symbol,Float64,Float64}[(:ES_A,:EP_A,p.k_A,k_A_r)]
    # inactive conformation: add each edge iff its tag permits it
    Stag != :OnlyA && push!(edges, (:E_I,:ES_I, FAST*S/K_I_S), (:ES_I,:E_I, FAST))   # S binds in I
    Ptag != :OnlyA && push!(edges, (:E_I,:EP_I, FAST*P/K_I_P), (:EP_I,:E_I, FAST))   # P binds in I
    if cattag != :OnlyA                                                              # catalysis runs in I
        push!(edges, (:ES_I,:EP_I, k_I), (:EP_I,:ES_I, k_I_r))
        push!(cat_edges, (:ES_I,:EP_I,k_I,k_I_r))
    end
    r = solvenet(species, edges, cat_edges, p.E_total)
    (flux=r.flux, conc=r.conc, extras=(k_I=k_I, k_I_r=k_I_r, K_I_S=K_I_S))
end

# the prototype's derivation (verbatim logic)
function v_new(tags, p, concs)
    Stag, cattag, Ptag = tags; S, P = concs.S, concs.P
    K_I_S = Stag==:EqualAI ? p.K_A_S : p.K_I_S
    k_I   = cattag==:EqualAI ? p.k_A : p.k_I
    K_I_P = Ptag==:EqualAI ? p.K_A_P : p.K_I_P
    k_A_r = p.k_A*p.K_A_P/(p.K_A_S*p.Keq)
    N_A = p.k_A*S/p.K_A_S - k_A_r*P/p.K_A_P; D_A = 1 + S/p.K_A_S + P/p.K_A_P
    ES_I = Stag!=:OnlyA; EP_I = Ptag!=:OnlyA; Icat = cattag!=:OnlyA && ES_I && EP_I
    D_I = 1.0; ES_I && (D_I += S/K_I_S); EP_I && (D_I += P/K_I_P)
    if Icat; k_I_r = k_I*K_I_P/(K_I_S*p.Keq); N_I = k_I*S/K_I_S - k_I_r*P/K_I_P
    else; N_I = 0.0; end
    p.E_total*(N_A + p.L*N_I)/(D_A + p.L*D_I)
end

draw(rng) = (K_A_S=0.5+2rand(rng),k_A=0.5+2rand(rng),K_A_P=0.5+2rand(rng),
             K_I_S=0.5+2rand(rng),k_I=0.5+2rand(rng),K_I_P=0.5+2rand(rng),
             L=0.5+rand(rng),Keq=2.0+2rand(rng),E_total=1.0)

trap = (:EqualAI,:NonequalAI,:OnlyA)

println("="^92)
println("INDEPENDENT confirmation: full-6-form network built from scratch")
println("="^92)

# 1. Cross-check full6 vs v_new on the trap combo (should DIVERGE — the break)
rng = MersenneTwister(hash(trap)); mr=0.0
for _ in 1:10
    p=draw(rng); c=(S=0.5+2rand(rng),P=0.5+2rand(rng))
    global mr = max(mr, abs(v_new(trap,p,c)-full6(trap,p,c).flux)/abs(full6(trap,p,c).flux))
end
@printf("1. trap combo %-30s : v_new vs full6-from-scratch max relerr = %.3f  (%.0f%%)\n",
        string(trap), mr, 100mr)

# 2. Cross-check full6 vs v_new on a NO-TRAP combo (P binds in I) — should AGREE
notrap = (:EqualAI,:NonequalAI,:EqualAI)
rng = MersenneTwister(hash(notrap)); mr2=0.0
for _ in 1:10
    p=draw(rng); c=(S=0.5+2rand(rng),P=0.5+2rand(rng))
    global mr2 = max(mr2, abs(v_new(notrap,p,c)-full6(notrap,p,c).flux)/abs(full6(notrap,p,c).flux))
end
@printf("2. no-trap combo %-27s : v_new vs full6-from-scratch max relerr = %.2e  (agree)\n",
        string(notrap), mr2)

# 3. Show the trapped EP_I holds real mass, ratio [EP_I]/[ES_I] should equal k_I/k_I_r
let p=draw(MersenneTwister(9)), c=(S=1.3,P=0.6)
    r = full6(trap,p,c)
    ratio_num = r.conc[:EP_I]/r.conc[:ES_I]
    ratio_theory = r.extras.k_I/r.extras.k_I_r
    @printf("3. trapped mass: [EP_I]=%.4g  [ES_I]=%.4g   [EP_I]/[ES_I]=%.4g  vs  k_I/k_I_r=%.4g  match=%s\n",
            r.conc[:EP_I], r.conc[:ES_I], ratio_num, ratio_theory,
            isapprox(ratio_num,ratio_theory;rtol=1e-4))
    frac_trapped = r.conc[:EP_I]/sum(values(r.conc))
    @printf("   EP_I holds %.1f%% of total enzyme — mass v_new's D_I completely omits.\n", 100frac_trapped)
end

# 4. THERMODYNAMIC check: true physics must give v=0 at P/S=Keq (active cycle Haldane holds)
let p=draw(MersenneTwister(3)), Seq=1.3
    c=(S=Seq, P=p.Keq*Seq)
    @printf("4. equilibrium P/S=Keq: full6 flux = %.2e  (true physics still gives v=0 — reaction is thermodynamically consistent)\n",
            full6(trap,p,c).flux)
end

# 5. What D_I SHOULD be (incl. trapped EP_I) vs what v_new uses
let p=draw(MersenneTwister(9)), c=(S=1.3,P=0.6)
    r = full6(trap,p,c)
    # per-state-normalized D_I from the network: solve I-conformation alone (E_I,ES_I,EP_I), Efree=1 basis
    K_I_S=r.extras.K_I_S; k_I=r.extras.k_I; k_I_r=r.extras.k_I_r
    D_I_true = 1 + c.S/K_I_S + (c.S/K_I_S)*(k_I/k_I_r)   # 1 + [ES_I]/[E_I] + [EP_I]/[E_I]
    D_I_vnew = 1 + c.S/K_I_S                             # what v_new uses (EP_I dropped)
    @printf("5. D_I(true, incl. trapped EP_I)=%.4g   D_I(v_new, EP_I dropped)=%.4g   v_new's D_I is %.0f%% too small\n",
            D_I_true, D_I_vnew, 100*(D_I_true-D_I_vnew)/D_I_true)
end
