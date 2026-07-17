# ABOUTME: Ping-pong with Keq-CONSISTENT params: v=0 at true eq (A*B/(P*Q)=1/Keq_pp), and
# ABOUTME: NORM correct / SKIP wrong off-equilibrium. Closes the ping-pong consistency loop.
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
    D=zeros(BigFloat,n); for i in 1:n; keep=[j for j in 1:n if j!=i]; D[i]=det(K[keep,keep]); end
    N=sum(kf*D[idx[r]]-kr*D[idx[p]] for (r,p,kf,kr) in cat_edges); (N=N,D=sum(D),dfree=D[idx[free]])
end
vnorm(A,I,L)=(A.N/A.dfree+L*I.N/I.dfree)/(A.D/A.dfree+L*I.D/I.dfree)
vskip(A,I,L)=(A.N+L*I.N)/(A.D+L*I.D)
const FAST=big"1e7"; bf(x)=big(x)
conf(p,kc1,kc1r,kc2,kc2r,A,B,P,Q)=ND([:E,:EA,:F,:FB],
   [(:E,:EA,p.konA*A),(:EA,:E,p.koffA),(:EA,:F,kc1),(:F,:EA,kc1r*P),
    (:F,:FB,p.konB*B),(:FB,:F,p.koffB),(:FB,:E,kc2),(:E,:FB,kc2r*Q)],
   [(:EA,:F,kc1,kc1r*P),(:FB,:E,kc2,kc2r*Q)],:E)
function oracle(p,L,A,B,P,Q)
    flux([:E_A,:EA_A,:F_A,:FB_A,:E_I,:EA_I,:F_I,:FB_I],
      [(:E_A,:EA_A,p.konA*A),(:EA_A,:E_A,p.koffA),(:EA_A,:F_A,p.kc1A),(:F_A,:EA_A,p.kc1rA*P),
       (:F_A,:FB_A,p.konB*B),(:FB_A,:F_A,p.koffB),(:FB_A,:E_A,p.kc2A),(:E_A,:FB_A,p.kc2rA*Q),
       (:E_I,:EA_I,p.konA*A),(:EA_I,:E_I,p.koffA),(:EA_I,:F_I,p.kc1I),(:F_I,:EA_I,p.kc1rI*P),
       (:F_I,:FB_I,p.konB*B),(:FB_I,:F_I,p.koffB),(:FB_I,:E_I,p.kc2I),(:E_I,:FB_I,p.kc2rI*Q),
       (:E_A,:E_I,FAST*L),(:E_I,:E_A,FAST)],
      [(:EA_A,:F_A,p.kc1A,p.kc1rA*P),(:FB_A,:E_A,p.kc2A,p.kc2rA*Q),
       (:EA_I,:F_I,p.kc1I,p.kc1rI*P),(:FB_I,:E_I,p.kc2I,p.kc2rI*Q)],big(1))
end
# Keq-consistent params: enforce SAME Keq_pp in both conformations.
Keqpp=bf(4.0)
p0=(konA=bf(1.3),koffA=bf(1.1),konB=bf(0.9),koffB=bf(1.2),kc1A=bf(2.1),kc2A=bf(1.7),kc1I=bf(0.7),kc2I=bf(0.6))
cA=sqrt(p0.kc1A*p0.kc2A*p0.konA*p0.konB/(Keqpp*p0.koffA*p0.koffB))
cI=sqrt(p0.kc1I*p0.kc2I*p0.konA*p0.konB/(Keqpp*p0.koffA*p0.koffB))
p=(p0...,kc1rA=cA,kc2rA=cA,kc1rI=cI,kc2rI=cI)
A0,P0,Q0=bf(1.1),bf(0.6),bf(0.5); Beq=P0*Q0/(Keqpp*A0)   # A*B/(P*Q)=1/Keq_pp
println("Ping-pong, Keq-consistent params. True equilibrium at A*B/(P*Q)=1/Keq_pp=", Float64(1/Keqpp))
@printf "  A*B/(P*Q)=%.5f  v(eq): L=0.5 -> %.2e ,  L=2.0 -> %.2e  (both ~0 => oracle CONSISTENT)\n" Float64(A0*Beq/(P0*Q0)) Float64(oracle(p,bf(0.5),A0,Beq,P0,Q0)) Float64(oracle(p,bf(2.0),A0,Beq,P0,Q0))
println("\nOff-equilibrium (Keq-consistent params): NORM vs SKIP over B")
@printf "%-6s %-13s %-11s %-11s %-9s\n" "B" "oracle" "NORMerr" "SKIPerr" "dfA/dfI"
for B in (bf(0.3),bf(1.0),bf(3.0),bf(8.0))
    Ap=conf(p,p.kc1A,p.kc1rA,p.kc2A,p.kc2rA,A0,B,P0,Q0); Ip=conf(p,p.kc1I,p.kc1rI,p.kc2I,p.kc2rI,A0,B,P0,Q0)
    o=oracle(p,bf(0.8),A0,B,P0,Q0)
    @printf "%-6.2f %-13.7f %-11.2e %-11.2e %-9.4f\n" Float64(B) Float64(o) Float64(abs(vnorm(Ap,Ip,bf(0.8))-o)/abs(o)) Float64(abs(vskip(Ap,Ip,bf(0.8))-o)/abs(o)) Float64(Ap.dfree/Ip.dfree)
end
