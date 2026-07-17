# ABOUTME: (1) attack.jl's cat:OnlyA mechanism: does SKIP break too? (2) ping-pong v=0 at eq.
using LinearAlgebra, Printf
setprecision(BigFloat, 200)
function flux(species, edges, cat_edges, Etot)
    n=length(species); idx=Dict(s=>i for (i,s) in enumerate(species)); M=zeros(BigFloat,n,n)
    for (a,b,r) in edges; M[idx[b],idx[a]]+=r; M[idx[a],idx[a]]-=r; end
    A=copy(M); A[1,:].=1; rhs=zeros(BigFloat,n); rhs[1]=Etot; c=A\rhs
    sum(kf*c[idx[r]]-kr*c[idx[p]] for (r,p,kf,kr) in cat_edges)
end
function ND(forms, edges, cat_edges, free)
    n=length(forms); idx=Dict(f=>i for (i,f) in enumerate(forms)); K=zeros(BigFloat,n,n)
    for (a,b,r) in edges; K[idx[a],idx[a]]+=r; K[idx[b],idx[a]]-=r; end
    D=zeros(BigFloat,n); for i in 1:n; keep=[j for j in 1:n if j!=i]; D[i]= isempty(keep) ? big(1) : det(K[keep,keep]); end
    N= isempty(cat_edges) ? big(0) : sum(kf*D[idx[r]]-kr*D[idx[p]] for (r,p,kf,kr) in cat_edges)
    (N=N,D=sum(D),dfree=D[idx[free]])
end
vnorm(A,I,L)=(A.N/A.dfree+L*I.N/I.dfree)/(A.D/A.dfree+L*I.D/I.dfree)
vskip(A,I,L)=(A.N+L*I.N)/(A.D+L*I.D)
const FAST=big"1e7"

# ---- attack.jl mechanism: S:OnlyA, cat:OnlyA, P:EqualAI. inactive = {E_I,EP_I} ----
A_full(kon,koff,KB,KP,k,kr,S,B,P)=ND([:E,:ES,:ESB,:EP],
    [(:E,:ES,kon*S),(:ES,:E,koff),(:ES,:ESB,FAST*B/KB),(:ESB,:ES,FAST),
     (:E,:EP,FAST*P/KP),(:EP,:E,FAST),(:ESB,:EP,k),(:EP,:ESB,kr)],[(:ESB,:EP,k,kr)],:E)
I_EP(KP,P)=ND([:E,:EP],[(:E,:EP,FAST*P/KP),(:EP,:E,FAST)],Tuple{Symbol,Symbol,BigFloat,BigFloat}[],:E)
oracle_catOnlyA(kon,koff,KB,KP,k,kr,L,S,B,P)=flux(
    [:E_A,:ES_A,:ESB_A,:EP_A,:E_I,:EP_I],
    [(:E_A,:ES_A,kon*S),(:ES_A,:E_A,koff),(:ES_A,:ESB_A,FAST*B/KB),(:ESB_A,:ES_A,FAST),
     (:E_A,:EP_A,FAST*P/KP),(:EP_A,:E_A,FAST),(:ESB_A,:EP_A,k),(:EP_A,:ESB_A,kr),
     (:E_I,:EP_I,FAST*P/KP),(:EP_I,:E_I,FAST),
     (:E_A,:E_I,FAST*L),(:E_I,:E_A,FAST),(:EP_A,:EP_I,FAST*L),(:EP_I,:EP_A,FAST)],
    [(:ESB_A,:EP_A,k,kr)],big(1))

println("attack.jl mechanism (S:OnlyA, cat:OnlyA, P:EqualAI) -- sweep B: dfA/dfI, NORM, SKIP vs oracle")
kon,koff,KB,KP,k,L,Keq,S,P=big(1.7),big(1.1),big(0.8),big(0.9),big(2.1),big(0.7),big(3.0),big(1.1),big(0.6)
kr=k*kon*KP/(koff*KB*Keq)
@printf "%-7s %-16s %-12s %-11s %-11s\n" "B" "dfA/dfI" "reqdL'=L*dfI/dfA" "NORMerr" "SKIPerr"
rats=BigFloat[]
for B in (big(0.1),big(0.3),big(1.0),big(3.0),big(8.0))
    A=A_full(kon,koff,KB,KP,k,kr,S,B,P); I=I_EP(KP,P)
    o=oracle_catOnlyA(kon,koff,KB,KP,k,kr,L,S,B,P); r=A.dfree/I.dfree; push!(rats,r)
    @printf "%-7.2f %-16.4e %-12.4e %-11.2e %-11.2e\n" Float64(B) Float64(r) Float64(L*I.dfree/A.dfree) Float64(abs(vnorm(A,I,L)-o)/abs(o)) Float64(abs(vskip(A,I,L)-o)/abs(o))
end
@printf "reqd-L' concentration span: %.4e .. %.4e  (%.3fx across B)\n" Float64(L*minimum(rats)/1) Float64(L*maximum(rats)/1) Float64(maximum(rats)/minimum(rats))

# ---- ping-pong: is the free-flip oracle thermodynamically consistent (v=0 at eq)? ----
println("\nPing-pong :NonequalAI free-flip oracle: v at equilibrium (Q/P = Keq_pp) for L=0.5,2.0")
function pp(species,edges,cat,Et); flux(species,edges,cat,Et); end
function pp_oracle(p,L,A,B,P,Q)
    flux([:E_A,:EA_A,:F_A,:FB_A,:E_I,:EA_I,:F_I,:FB_I],
      [(:E_A,:EA_A,p.konA*A),(:EA_A,:E_A,p.koffA),(:EA_A,:F_A,p.kc1A),(:F_A,:EA_A,p.kc1rA*P),
       (:F_A,:FB_A,p.konB*B),(:FB_A,:F_A,p.koffB),(:FB_A,:E_A,p.kc2A),(:E_A,:FB_A,p.kc2rA*Q),
       (:E_I,:EA_I,p.konA*A),(:EA_I,:E_I,p.koffA),(:EA_I,:F_I,p.kc1I),(:F_I,:EA_I,p.kc1rI*P),
       (:F_I,:FB_I,p.konB*B),(:FB_I,:F_I,p.koffB),(:FB_I,:E_I,p.kc2I),(:E_I,:FB_I,p.kc2rI*Q),
       (:E_A,:E_I,FAST*L),(:E_I,:E_A,FAST)],
      [(:EA_A,:F_A,p.kc1A,p.kc1rA*P),(:FB_A,:E_A,p.kc2A,p.kc2rA*Q),
       (:EA_I,:F_I,p.kc1I,p.kc1rI*P),(:FB_I,:E_I,p.kc2I,p.kc2rI*Q)],big(1))
end
# Build a Haldane-consistent :NonequalAI param set: each conformation's own Keq must match.
# Reaction A+B<->P+Q. For conf X: Keq_X = (kc1*kc2)/(kc1r*kc2r) * (1/(KA_bind..)) -- enforce by
# choosing reverse cats so both conformations share Keq_pp := kc1*kc2*konA*konB/(kc1r*kc2r*koffA*koffB).
bf(x)=big(x)
p0=(konA=bf(1.3),koffA=bf(1.1),konB=bf(0.9),koffB=bf(1.2),
    kc1A=bf(2.1),kc2A=bf(1.7),kc1I=bf(0.7),kc2I=bf(0.6))
# pick Keq_pp, then set reverse catalytic constants so BOTH confs obey it:
Keqpp=bf(4.0)
# per conf: kc1*kc2*konA*konB/(kc1r*kc2r*koffA*koffB) = Keqpp  => kc1r*kc2r = kc1*kc2*konA*konB/(Keqpp*koffA*koffB)
# split evenly: kc1r=kc2r=sqrt(that)
prodA=p0.kc1A*p0.kc2A*p0.konA*p0.konB/(Keqpp*p0.koffA*p0.koffB); cA=sqrt(prodA)
prodI=p0.kc1I*p0.kc2I*p0.konA*p0.konB/(Keqpp*p0.koffA*p0.koffB); cI=sqrt(prodI)
p=(p0...,kc1rA=cA,kc2rA=cA,kc1rI=cI,kc2rI=cI)
A0,P0,Q0=bf(1.1),bf(0.6),bf(0.5); B0 = Keqpp*P0*Q0/(A0)   # A*B/(P*Q)=Keqpp => B=Keqpp*P*Q/A
@printf "  A*B/(P*Q)=%.4f (Keq_pp=%.1f)\n" Float64(A0*B0/(P0*Q0)) Float64(Keqpp)
for L in (bf(0.5),bf(2.0))
    @printf "  L=%.1f : oracle v(eq)=%.3e\n" Float64(L) Float64(pp_oracle(p,L,A0,B0,P0,Q0))
end
