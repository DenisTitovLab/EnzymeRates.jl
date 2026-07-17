# ABOUTME: Which koff_I convention is THERMODYNAMICALLY CONSISTENT? A legal network gives
# ABOUTME: v=0 at equilibrium for ALL L. Find equilibrium (active v=0), test 2-conf v vs L.
using LinearAlgebra, Printf
setprecision(BigFloat, 200)
function flux(species, edges, cat_edges, Etot)
    n=length(species); idx=Dict(s=>i for (i,s) in enumerate(species)); M=zeros(BigFloat,n,n)
    for (a,b,r) in edges; M[idx[b],idx[a]]+=r; M[idx[a],idx[a]]-=r; end
    A=copy(M); A[1,:].=1; rhs=zeros(BigFloat,n); rhs[1]=Etot; c=A\rhs
    sum(kf*c[idx[r]]-kr*c[idx[p]] for (r,p,kf,kr) in cat_edges)
end
const FAST=big"1e7"; const KONI=big"1e-10"
active(kon,koff,KB,KP,k,kr,S,B,P)=flux([:E,:ES,:ESB,:EP],
    [(:E,:ES,kon*S),(:ES,:E,koff),(:ES,:ESB,FAST*B/KB),(:ESB,:ES,FAST),
     (:E,:EP,FAST*P/KP),(:EP,:E,FAST),(:ESB,:EP,k),(:EP,:ESB,kr)],[(:ESB,:EP,k,kr)],big(1))
function twoconf(kon,koff,KB,KP,k,kr,L,S,B,P,koffI)
    flux([:E_A,:ES_A,:ESB_A,:EP_A,:E_I,:ES_I,:ESB_I,:EP_I],
      [(:E_A,:ES_A,kon*S),(:ES_A,:E_A,koff),(:ES_A,:ESB_A,FAST*B/KB),(:ESB_A,:ES_A,FAST),
       (:E_A,:EP_A,FAST*P/KP),(:EP_A,:E_A,FAST),(:ESB_A,:EP_A,k),(:EP_A,:ESB_A,kr),
       (:E_I,:ES_I,KONI*S),(:ES_I,:E_I,koffI),(:ES_I,:ESB_I,FAST*B/KB),(:ESB_I,:ES_I,FAST),
       (:E_I,:EP_I,FAST*P/KP),(:EP_I,:E_I,FAST),(:ESB_I,:EP_I,k),(:EP_I,:ESB_I,kr),
       (:E_A,:E_I,FAST*L),(:E_I,:E_A,FAST),(:EP_A,:EP_I,FAST*L),(:EP_I,:EP_A,FAST)],
      [(:ESB_A,:EP_A,k,kr),(:ESB_I,:EP_I,k,kr)],big(1))
end
kon,koff,KB,KP,k,L,Keq=big(1.7),big(1.1),big(0.8),big(0.9),big(2.1),big(0.7),big(3.0)
kr=k*kon*KP/(koff*KB*Keq)
# find equilibrium B* (active flux=0) at fixed S,P by bisection
S,P=big(1.1),big(0.6)
function bisect_eq(kon,koff,KB,KP,k,kr,S,P)
    lo,hi=big(1e-6),big(1e3)
    for _ in 1:200; mid=sqrt(lo*hi); active(kon,koff,KB,KP,k,kr,S,mid,P)>0 ? (hi=mid) : (lo=mid); end
    sqrt(lo*hi)
end
Bstar=bisect_eq(kon,koff,KB,KP,k,kr,S,P)
@printf "Equilibrium B* (active flux=0): B*=%.8f   check active(B*)=%.2e   S*B*/P=%.5f  (1/Keq=%.5f)\n" Float64(Bstar) Float64(active(kon,koff,KB,KP,k,kr,S,Bstar,P)) Float64(S*Bstar/P) Float64(1/Keq)
println("\nAt equilibrium B*, a THERMODYNAMICALLY CONSISTENT network gives v=0 for ALL L.")
println("Two-conformation v(B*) for each koff_I convention, L=0.7 and L=3.0:\n")
@printf "%-26s %-15s %-15s %-10s\n" "koff_I convention" "v(B*,L=0.7)" "v(B*,L=3.0)" "verdict"
for (lbl,koffI) in [("->0 (DELETE)",big"1e-12"),("=koff_A",big(1.1)),("->inf (fast release)",big"1e12")]
    v1=twoconf(kon,koff,KB,KP,k,kr,big(0.7),S,Bstar,P,koffI)
    v2=twoconf(kon,koff,KB,KP,k,kr,big(3.0),S,Bstar,P,koffI)
    verdict = (abs(v1)<1e-6 && abs(v2)<1e-6) ? "CONSISTENT" : "VIOLATES 2nd law"
    @printf "%-26s %-15.2e %-15.2e %-10s\n" lbl Float64(v1) Float64(v2) verdict
end
println("\n(A nonzero v at equilibrium = the inactive cycle pumps flux with no thermodynamic drive.)")
