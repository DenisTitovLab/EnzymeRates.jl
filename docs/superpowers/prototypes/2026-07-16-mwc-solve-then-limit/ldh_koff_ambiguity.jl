# ABOUTME: Is the LDH :OnlyA SS-binding limit well-defined? Sweep the inactive S-release
# ABOUTME: rate koff_I (kon_I->0 fixed). If oracle flux depends on koff_I, the limit is ambiguous.
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
    N=sum(kf*D[idx[r]]-kr*D[idx[p]] for (r,p,kf,kr) in cat_edges)
    (N=N,D=sum(D),dfree=D[idx[free]])
end
vnorm(A,I,L,Et)=Et*(A.N/A.dfree+L*I.N/I.dfree)/(A.D/A.dfree+L*I.D/I.dfree)
vskip(A,I,L,Et)=Et*(A.N+L*I.N)/(A.D+L*I.D)
const FAST=big"1e7"; const KONI=big"1e-10"   # inactive S-assoc ~ 0 (K_I -> inf)

# I-state pieces with a chosen inactive S-release rate koffI (kon_I fixed ~0)
Ipieces(koff,KB,KP,k,kr,S,B,P,koffI)=ND([:E,:ES,:ESB,:EP],
    [(:E,:ES,KONI*S),(:ES,:E,koffI),(:ES,:ESB,FAST*B/KB),(:ESB,:ES,FAST),
     (:E,:EP,FAST*P/KP),(:EP,:E,FAST),(:ESB,:EP,k),(:EP,:ESB,kr)],[(:ESB,:EP,k,kr)],:E)
Apieces(kon,koff,KB,KP,k,kr,S,B,P)=ND([:E,:ES,:ESB,:EP],
    [(:E,:ES,kon*S),(:ES,:E,koff),(:ES,:ESB,FAST*B/KB),(:ESB,:ES,FAST),
     (:E,:EP,FAST*P/KP),(:EP,:E,FAST),(:ESB,:EP,k),(:EP,:ESB,kr)],[(:ESB,:EP,k,kr)],:E)
function oracle(kon,koff,KB,KP,k,kr,L,S,B,P,koffI)
    flux([:E_A,:ES_A,:ESB_A,:EP_A,:E_I,:ES_I,:ESB_I,:EP_I],
      [(:E_A,:ES_A,kon*S),(:ES_A,:E_A,koff),(:ES_A,:ESB_A,FAST*B/KB),(:ESB_A,:ES_A,FAST),
       (:E_A,:EP_A,FAST*P/KP),(:EP_A,:E_A,FAST),(:ESB_A,:EP_A,k),(:EP_A,:ESB_A,kr),
       (:E_I,:ES_I,KONI*S),(:ES_I,:E_I,koffI),(:ES_I,:ESB_I,FAST*B/KB),(:ESB_I,:ES_I,FAST),
       (:E_I,:EP_I,FAST*P/KP),(:EP_I,:E_I,FAST),(:ESB_I,:EP_I,k),(:EP_I,:ESB_I,kr),
       (:E_A,:E_I,FAST*L),(:E_I,:E_A,FAST),(:EP_A,:EP_I,FAST*L),(:EP_I,:EP_A,FAST)],
      [(:ESB_A,:EP_A,k,kr),(:ESB_I,:EP_I,k,kr)],big(1.0))
end

kon,koff,KB,KP,k,L,Keq,S,B,P=big(1.7),big(1.1),big(0.8),big(0.9),big(2.1),big(0.7),big(3.0),big(1.1),big(0.9),big(0.6)
kr=k*kon*KP/(koff*KB*Keq)
Ap=Apieces(kon,koff,KB,KP,k,kr,S,B,P)
println("Does the LDH :OnlyA limit depend on the inactive S-release rate koff_I?  (kon_I~0 fixed)")
println("koff_I convention        oracle_flux    N_I(full)   dfA/dfI    NORMerr   SKIPerr")
for (lbl,koffI) in [("->0 (step vanishes=DELETE)",big"1e-10"),("0.5*koff",big(0.55)),("=koff (my earlier test)",big(1.1)),("10*koff",big(11.0)),("->inf (fast release)",big"1e10")]
    Ip=Ipieces(koff,KB,KP,k,kr,S,B,P,koffI)
    o=oracle(kon,koff,KB,KP,k,kr,L,S,B,P,koffI)
    @printf "%-24s %-14.8f %-11.5f %-10.5f %-9.1e %-9.1e\n" lbl Float64(o) Float64(Ip.N) Float64(Ap.dfree/Ip.dfree) Float64(abs(vnorm(Ap,Ip,L,big(1))-o)/abs(o)) Float64(abs(vskip(Ap,Ip,L,big(1))-o)/abs(o))
end
println("\nInterpretation:")
println("  If oracle_flux varies across rows -> the :OnlyA SS limit is AMBIGUOUS (koff_I is a real")
println("  modeling choice: does inactive enzyme release :OnlyA substrate reached via reverse cat?).")
println("  NORMerr shows whether full-graph NORM tracks the oracle for each convention;")
println("  SKIPerr shows whether the un-normalized combine does.")
