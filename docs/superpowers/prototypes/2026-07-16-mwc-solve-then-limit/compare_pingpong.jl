# ABOUTME: Compare shipped EnzymeRates.rate_equation for allosteric ping-pong :NonequalAI
# ABOUTME: against the independent E-flip mass-action oracle (BigFloat), with full param mapping.
using EnzymeRates, LinearAlgebra, Printf, Random
const ER = EnzymeRates
setprecision(BigFloat, 200)

# ---- Independent oracle (copied verbatim from redesign_pingpong.jl) ----------
# 8-state E-only free-flip ping-pong. `sel` picks which catalytic edges to sum:
#   :both -> P-edges + Q-edges (=2v), :P -> P production only (=v), :Q -> Q only (=v)
function oracle_flux(p, L, A, B, P, Q; sel=:both)
    species = [:E_A,:EA_A,:F_A,:FB_A,:E_I,:EA_I,:F_I,:FB_I]
    edges = [(:E_A,:EA_A,p.konA*A),(:EA_A,:E_A,p.koffA),(:EA_A,:F_A,p.kc1A),(:F_A,:EA_A,p.kc1rA*P),
       (:F_A,:FB_A,p.konB*B),(:FB_A,:F_A,p.koffB),(:FB_A,:E_A,p.kc2A),(:E_A,:FB_A,p.kc2rA*Q),
       (:E_I,:EA_I,p.konA*A),(:EA_I,:E_I,p.koffA),(:EA_I,:F_I,p.kc1I),(:F_I,:EA_I,p.kc1rI*P),
       (:F_I,:FB_I,p.konB*B),(:FB_I,:F_I,p.koffB),(:FB_I,:E_I,p.kc2I),(:E_I,:FB_I,p.kc2rI*Q),
       (:E_A,:E_I,big"1e7"*L),(:E_I,:E_A,big"1e7")]
    Pedges = [(:EA_A,:F_A,p.kc1A,p.kc1rA*P),(:EA_I,:F_I,p.kc1I,p.kc1rI*P)]
    Qedges = [(:FB_A,:E_A,p.kc2A,p.kc2rA*Q),(:FB_I,:E_I,p.kc2I,p.kc2rI*Q)]
    cat = sel === :P ? Pedges : sel === :Q ? Qedges : vcat(Pedges,Qedges)
    n=length(species); idx=Dict(s=>i for (i,s) in enumerate(species)); M=zeros(BigFloat,n,n)
    for (a,b,r) in edges; M[idx[b],idx[a]]+=r; M[idx[a],idx[a]]-=r; end
    Amat=copy(M); Amat[1,:].=1; rhs=zeros(BigFloat,n); rhs[1]=big(1.0); c=Amat\rhs
    sum(kf*c[idx[r]]-kr*c[idx[p]] for (r,p,kf,kr) in cat)
end

# ---- Shipped mechanism -------------------------------------------------------
m = @allosteric_mechanism begin
    substrates: A, B ; products: P, Q ; catalytic_multiplicity: 1
    catalytic_steps: begin
        E + A <--> E(A)                                       :: EqualAI
        E(A) <--> E(; residual = A - P) + P                   :: NonequalAI
        E(; residual = A - P) + B <--> E(B; residual = A - P) :: EqualAI
        E(B; residual = A - P) <--> E + Q                     :: NonequalAI
    end
end
fpn = ER.fitted_params(m)

# Symbols carrying `+`/`-` in their names (residual) must be built with Symbol().
S(x) = Symbol(x)
kAQ = S("kon_A_Q_EB_res_+A_-P");  kAQr = S("koff_A_Q_EB_res_+A_-P")
kBon = S("kon_B_E_res_+A_-P");    kBoff = S("koff_B_E_res_+A_-P")
kIQ  = S("kon_I_Q_EB_res_+A_-P"); kIQr  = S("koff_I_Q_EB_res_+A_-P")

# Parameter mapping (fitted -> oracle rate constants), by MEANING:
#   kon_A_E      = konA (step1 fwd, A binding)          koff_A_E = koffA  [HALDANE-derived]
#   kon_A_P_EA   = kc1A (step2 fwd, EA->F+P active)     koff_A_P_EA = kc1rA
#   kon_A_Q_EB.. = kc2A (step4 fwd, FB->E+Q active)     koff_A_Q_EB.. = kc2rA
#   kon_B_E_res..= konB (step3 fwd, B binding)          koff_B_E_res..= koffB
#   kon_I_P_EA   = kc1I (step2 fwd inactive)            koff_I_P_EA = kc1rI [WEGSCHEIDER-derived]
#   kon_I_Q_EB.. = kc2I (step4 fwd inactive)            koff_I_Q_EB.. = kc2rI
bf(x)=big(x)

function run_draw(rng)
    # draw the 10 free rate constants + L + Keq (all O(1), positive)
    konA=0.5+2rand(rng); kc1A=0.5+2rand(rng); kc1rA=0.3+rand(rng)
    kc2A=0.5+2rand(rng); kc2rA=0.3+rand(rng); konB=0.5+2rand(rng); koffB=0.5+2rand(rng)
    kc1I=0.5+2rand(rng); kc2I=0.5+2rand(rng); kc2rI=0.3+rand(rng)
    L=0.5+rand(rng); Keq=2.0+2rand(rng)
    # Derived dependent constants, EXACTLY as the shipped reduced eq states them:
    # koff_A_E = (1/Keq)(1/kc1rA)(1/kc2rA)(1/koffB) konA kc1A kc2A konB
    koffA = (1/Keq)*(1/kc1rA)*(1/kc2rA)*(1/koffB)*konA*kc1A*kc2A*konB
    # koff_I_P_EA(=kc1rI) = kc1rA*kc2rA*(1/kc2rI)*(1/kc1A)*(1/kc2A)*kc1I*kc2I
    kc1rI = kc1rA*kc2rA*(1/kc2rI)*(1/kc1A)*(1/kc2A)*kc1I*kc2I

    # shipped param NamedTuple (fitted names + Keq + E_total)
    d = Dict(:kon_A_E=>konA, :kon_A_P_EA=>kc1A, :koff_A_P_EA=>kc1rA,
             kAQ=>kc2A, kAQr=>kc2rA, kBon=>konB, kBoff=>koffB,
             :kon_I_P_EA=>kc1I, kIQ=>kc2I, kIQr=>kc2rI, :L=>L)
    prm = NamedTuple{(fpn..., :Keq, :E_total)}(((d[s] for s in fpn)..., Keq, 1.0))

    # oracle param NamedTuple (BigFloat)
    p = (konA=bf(konA),koffA=bf(koffA),konB=bf(konB),koffB=bf(koffB),
         kc1A=bf(kc1A),kc1rA=bf(kc1rA),kc2A=bf(kc2A),kc2rA=bf(kc2rA),
         kc1I=bf(kc1I),kc1rI=bf(kc1rI),kc2I=bf(kc2I),kc2rI=bf(kc2rI))
    (prm=prm, p=p, L=bf(L), Keq=Keq)
end

println("="^100)
println("Allosteric ping-pong :NonequalAI  —  SHIPPED rate_equation vs E-flip oracle")
println("="^100)
@printf "%-4s %-14s %-14s %-14s %-11s %-11s %-9s\n" "draw" "shipped" "oracle_P(=v)" "oracle_both" "relerr_P" "relerr_½both" "shp/orclB"
rng = MersenneTwister(20260716)
maxrel = 0.0
for t in 1:12
    dr = run_draw(rng)
    A=0.5+2rand(rng); B=0.5+2rand(rng); P=0.3+rand(rng); Q=0.3+rand(rng)
    shp = real(ER.rate_equation(m, (A=A,B=B,P=P,Q=Q), dr.prm))
    oP = oracle_flux(dr.p, dr.L, bf(A),bf(B),bf(P),bf(Q); sel=:P)
    oQ = oracle_flux(dr.p, dr.L, bf(A),bf(B),bf(P),bf(Q); sel=:Q)
    oB = oracle_flux(dr.p, dr.L, bf(A),bf(B),bf(P),bf(Q); sel=:both)
    relP = Float64(abs(shp-oP)/abs(oP))
    relH = Float64(abs(shp-oB/2)/abs(oB/2))
    global maxrel = max(maxrel, relP)
    @printf "%-4d %-14.8f %-14.8f %-14.8f %-11.2e %-11.2e %-9.5f\n" t shp Float64(oP) Float64(oB) relP relH Float64(shp/oB)
    # sanity: P-edge flux == Q-edge flux (proves both==2v)
    @assert Float64(abs(oP-oQ)/abs(oP)) < 1e-30 "P/Q flux mismatch draw $t"
end
@printf "\nMAX relative error (shipped vs oracle_P=v): %.3e\n" maxrel

# ---- Equilibrium check: A*B/(P*Q) = 1/Keq  =>  both v = 0 --------------------
println("\nEquilibrium check (Q set so P*Q/(A*B) = Keq): shipped and oracle should be ~0")
rng2 = MersenneTwister(999)
for t in 1:4
    dr = run_draw(rng2)
    A=0.5+2rand(rng2); B=0.5+2rand(rng2); P=0.3+rand(rng2)
    Q = dr.Keq*A*B/P                       # forces P*Q/(A*B) = Keq
    shp = real(ER.rate_equation(m, (A=A,B=B,P=P,Q=Q), dr.prm))
    oB  = oracle_flux(dr.p, dr.L, bf(A),bf(B),bf(P),bf(Q); sel=:both)
    @printf "  draw %d  shipped=% .3e   oracle=% .3e\n" t shp Float64(oB)
end
