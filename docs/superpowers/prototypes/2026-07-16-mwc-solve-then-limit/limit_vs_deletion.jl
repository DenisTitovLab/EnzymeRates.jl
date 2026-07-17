# ABOUTME: Tests whether :OnlyA-as-LIMIT (task spec) == :OnlyA-as-DELETION (shipped harness).
# ABOUTME: If they agree, the normalization verdict rests on the right oracle.
using LinearAlgebra, Random, Printf

function flux(species, edges, cat_edges, Etot)
    n=length(species); idx=Dict(s=>i for (i,s) in enumerate(species)); M=zeros(n,n)
    for (a,b,r) in edges; M[idx[b],idx[a]]+=r; M[idx[a],idx[a]]-=r; end
    A=copy(M); A[1,:].=1.0; rhs=zeros(n); rhs[1]=Etot; c=A\rhs
    sum(kf*c[idx[r]]-kr*c[idx[p]] for (r,p,kf,kr) in cat_edges)
end
function cofactors(forms, edges)
    n=length(forms); idx=Dict(f=>i for (i,f) in enumerate(forms)); K=zeros(n,n)
    for (a,b,r) in edges; K[idx[a],idx[a]]+=r; K[idx[b],idx[a]]-=r; end
    D=zeros(n); for i in 1:n; keep=[j for j in 1:n if j!=i]; D[i]=isempty(keep) ? 1.0 : det(K[keep,keep]); end
    D, idx
end
function ND(forms, edges, cat_edges, free)
    D,idx=cofactors(forms,edges); N=isempty(cat_edges) ? 0.0 : sum(kf*D[idx[r]]-kr*D[idx[p]] for (r,p,kf,kr) in cat_edges)
    (N=N,D=sum(D),dfree=D[idx[free]])
end
v_norm(A,I,L)=(A.N/A.dfree+L*I.N/I.dfree)/(A.D/A.dfree+L*I.D/I.dfree)
v_unorm(A,I,L)=(A.N+L*I.N)/(A.D+L*I.D)
const FAST=1e7; const BIG=1e10

# -------- Uni-uni :OnlyA (S binds via RAPID-EQUILIBRIUM step) --------
# Deletion: no inactive S-binding edge. Limit: forward FAST*S/BIG (~0), reverse FAST.
function uni_oracle(KA,KP,k,kr;L,S,P,mode)
    sp=[:E_A,:ES_A,:EP_A,:E_I,:EP_I,:ES_I]
    ed=[(:E_A,:ES_A,FAST*S/KA),(:ES_A,:E_A,FAST),(:E_A,:EP_A,FAST*P/KP),(:EP_A,:E_A,FAST),
        (:E_I,:EP_I,FAST*P/KP),(:EP_I,:E_I,FAST),(:E_A,:E_I,FAST*L),(:E_I,:E_A,FAST),
        (:EP_A,:EP_I,FAST*L),(:EP_I,:EP_A,FAST),(:ES_A,:EP_A,k),(:EP_A,:ES_A,kr),
        (:ES_I,:EP_I,k),(:EP_I,:ES_I,kr)]
    if mode==:limit; append!(ed,[(:E_I,:ES_I,FAST*S/BIG),(:ES_I,:E_I,FAST)]); end
    flux(sp,ed,[(:ES_A,:EP_A,k,kr),(:ES_I,:EP_I,k,kr)],1.0)
end

# -------- LDH i-state (S binds via STEADY-STATE step kon*S/koff) --------
# Deletion: no inactive S-binding edge. Limit: forward kon*S/BIG (~0), reverse koff.
function ldh_oracle(kon,koff,KB,KP,k,kr;L,S,B,P,mode)
    sp=[:E_A,:ES_A,:ESB_A,:EP_A,:E_I,:EP_I,:ESB_I,:ES_I]
    ed=[(:E_A,:ES_A,kon*S),(:ES_A,:E_A,koff),(:ES_A,:ESB_A,FAST*B/KB),(:ESB_A,:ES_A,FAST),
        (:ES_I,:ESB_I,FAST*B/KB),(:ESB_I,:ES_I,FAST),(:E_A,:EP_A,FAST*P/KP),(:EP_A,:E_A,FAST),
        (:E_I,:EP_I,FAST*P/KP),(:EP_I,:E_I,FAST),(:E_A,:E_I,FAST*L),(:E_I,:E_A,FAST),
        (:EP_A,:EP_I,FAST*L),(:EP_I,:EP_A,FAST),(:ESB_A,:EP_A,k),(:EP_A,:ESB_A,kr),
        (:ESB_I,:EP_I,k),(:EP_I,:ESB_I,kr)]
    if mode==:limit; append!(ed,[(:E_I,:ES_I,kon*S/BIG),(:ES_I,:E_I,koff)]); end
    flux(sp,ed,[(:ESB_A,:EP_A,k,kr),(:ESB_I,:EP_I,k,kr)],1.0)
end
# package pieces (deletion, as derived)
ldh_A(kon,koff,KB,KP,k,kr,S,B,P)=ND([:E,:ES,:ESB,:EP],
    [(:E,:ES,kon*S),(:ES,:E,koff),(:ES,:ESB,FAST*B/KB),(:ESB,:ES,FAST),
     (:E,:EP,FAST*P/KP),(:EP,:E,FAST),(:ESB,:EP,k),(:EP,:ESB,kr)],[(:ESB,:EP,k,kr)],:E)
ldh_I(kon,koff,KB,KP,k,kr,S,B,P)=ND([:E,:EP,:ESB,:ES],
    [(:ES,:ESB,FAST*B/KB),(:ESB,:ES,FAST),(:E,:EP,FAST*P/KP),(:EP,:E,FAST),
     (:ESB,:EP,k),(:EP,:ESB,kr)],[(:ESB,:EP,k,kr)],:E)

println("="^76)
println("Uni-uni :OnlyA (RE S-binding): does LIMIT oracle == DELETION oracle?")
println("="^76)
@printf "%-4s %-14s %-14s %-12s\n" "d" "oracle_delete" "oracle_limit" "rel-diff"
rng=MersenneTwister(3)
for t in 1:6
    KA=.5+2rand(rng);KP=.5+2rand(rng);k=.5+2rand(rng);L=.5+rand(rng);Keq=2+2rand(rng);S=.5+2rand(rng);P=.5+2rand(rng)
    kr=k*KP/(Keq*KA)
    od=uni_oracle(KA,KP,k,kr;L=L,S=S,P=P,mode=:delete)
    ol=uni_oracle(KA,KP,k,kr;L=L,S=S,P=P,mode=:limit)
    @printf "%-4d %-14.6g %-14.6g %-12.3e\n" t od ol abs(od-ol)/abs(od)
end

println("\n"*"="^76)
println("LDH i-state (SS S-binding): LIMIT vs DELETION oracle, and v_norm vs BOTH")
println("="^76)
@printf "%-4s %-13s %-13s %-11s %-11s %-11s\n" "d" "orcl_delete" "orcl_limit" "del-lim" "vnorm-del" "vnorm-lim"
rng=MersenneTwister(1)
for t in 1:6
    kon=.5+2rand(rng);koff=.5+2rand(rng);KB=.5+2rand(rng);KP=.5+2rand(rng)
    k=.5+2rand(rng);L=.5+rand(rng);Keq=2+2rand(rng);S=.5+2rand(rng);B=.5+2rand(rng);P=.5+2rand(rng)
    kr=k*kon*KP/(koff*KB*Keq)
    od=ldh_oracle(kon,koff,KB,KP,k,kr;L=L,S=S,B=B,P=P,mode=:delete)
    ol=ldh_oracle(kon,koff,KB,KP,k,kr;L=L,S=S,B=B,P=P,mode=:limit)
    A=ldh_A(kon,koff,KB,KP,k,kr,S,B,P); I=ldh_I(kon,koff,KB,KP,k,kr,S,B,P)
    vn=v_norm(A,I,L)
    @printf "%-4d %-13.6g %-13.6g %-11.2e %-11.2e %-11.2e\n" t od ol abs(od-ol)/abs(od) abs(vn-od)/abs(od) abs(vn-ol)/abs(ol)
end

# ---- Keq-convention check: where is oracle == 0 ? ----
println("\n"*"="^76); println("Keq-convention / v=0 discrimination check (deletion oracle)"); println("="^76)
kon,koff,KB,KP,k,L,Keq,S,P=1.7,1.1,0.8,0.9,2.1,0.7,3.0,1.1,0.6; kr=k*kon*KP/(koff*KB*Keq)
for (lbl,Bval) in [("S*B/P=Keq", Keq*P/S),("S*B/P=1/Keq",(1/Keq)*P/S)]
    A=ldh_A(kon,koff,KB,KP,k,kr,S,Bval,P); I=ldh_I(kon,koff,KB,KP,k,kr,S,Bval,P)
    o=ldh_oracle(kon,koff,KB,KP,k,kr;L=L,S=S,B=Bval,P=P,mode=:delete)
    @printf "%-14s B=%-8.4f oracle=%-11.3e vnorm=%-11.3e vunorm=%-11.3e\n" lbl Bval o v_norm(A,I,L) v_unorm(A,I,L)
end
