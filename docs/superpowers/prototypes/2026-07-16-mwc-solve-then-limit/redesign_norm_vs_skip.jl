# ABOUTME: The redesign's ACTUAL open question: full-graph pieces + :OnlyA LIMIT oracle.
# ABOUTME: Does the un-normalized (SKIP) combine match the limit oracle, or only NORM?
using LinearAlgebra, Printf, Random
setprecision(BigFloat, 200)

# steady-state flux of an explicit network (BigFloat to survive the K_I->inf stiffness)
function flux(species, edges, cat_edges, Etot)
    n=length(species); idx=Dict(s=>i for (i,s) in enumerate(species)); M=zeros(BigFloat,n,n)
    for (a,b,r) in edges; M[idx[b],idx[a]]+=r; M[idx[a],idx[a]]-=r; end
    A=copy(M); A[1,:].=1; rhs=zeros(BigFloat,n); rhs[1]=Etot; c=A\rhs
    sum(kf*c[idx[r]]-kr*c[idx[p]] for (r,p,kf,kr) in cat_edges)
end
# single-conformation King-Altman pieces (N, D, dfree) via matrix-tree cofactors
function ND(forms, edges, cat_edges, free)
    n=length(forms); idx=Dict(f=>i for (i,f) in enumerate(forms)); K=zeros(BigFloat,n,n)
    for (a,b,r) in edges; K[idx[a],idx[a]]+=r; K[idx[b],idx[a]]-=r; end
    D=zeros(BigFloat,n); for i in 1:n; keep=[j for j in 1:n if j!=i]; D[i]=det(K[keep,keep]); end
    N=sum(kf*D[idx[r]]-kr*D[idx[p]] for (r,p,kf,kr) in cat_edges)
    (N=N,D=sum(D),dfree=D[idx[free]])
end
vnorm(A,I,L,Et)=Et*(A.N/A.dfree+L*I.N/I.dfree)/(A.D/A.dfree+L*I.D/I.dfree)
vskip(A,I,L,Et)=Et*(A.N+L*I.N)/(A.D+L*I.D)

const FAST=big"1e6"

# ---- LDH i-state, REDESIGN frame: I-state is the FULL graph with S-binding at
#      near-:OnlyA (kon_I = kon/BIG, koff_I = koff).  B,P,cat are EqualAI (=A). ----
function pieces(kon,koff,KB,KP,k,kr,S,B,P,BIG)
    A=ND([:E,:ES,:ESB,:EP],
        [(:E,:ES,kon*S),(:ES,:E,koff),(:ES,:ESB,FAST*B/KB),(:ESB,:ES,FAST),
         (:E,:EP,FAST*P/KP),(:EP,:E,FAST),(:ESB,:EP,k),(:EP,:ESB,kr)],[(:ESB,:EP,k,kr)],:E)
    # FULL inactive graph (undeleted): S-binding present but forward rate ~0 (K_I->inf)
    I=ND([:E,:ES,:ESB,:EP],
        [(:E,:ES,kon*S/BIG),(:ES,:E,koff),(:ES,:ESB,FAST*B/KB),(:ESB,:ES,FAST),
         (:E,:EP,FAST*P/KP),(:EP,:E,FAST),(:ESB,:EP,k),(:EP,:ESB,kr)],[(:ESB,:EP,k,kr)],:E)
    A,I
end
# LIMIT oracle: two-conformation free-flip (E and EP segment flip, ratio L); I-state
# S-binding near-:OnlyA (kon/BIG forward, koff reverse -- the koff drain SURVIVES).
function oracle(kon,koff,KB,KP,k,kr,L,S,B,P,BIG)
    flux([:E_A,:ES_A,:ESB_A,:EP_A,:E_I,:ES_I,:ESB_I,:EP_I],
        [(:E_A,:ES_A,kon*S),(:ES_A,:E_A,koff),(:ES_A,:ESB_A,FAST*B/KB),(:ESB_A,:ES_A,FAST),
         (:E_A,:EP_A,FAST*P/KP),(:EP_A,:E_A,FAST),(:ESB_A,:EP_A,k),(:EP_A,:ESB_A,kr),
         (:E_I,:ES_I,kon*S/BIG),(:ES_I,:E_I,koff),(:ES_I,:ESB_I,FAST*B/KB),(:ESB_I,:ES_I,FAST),
         (:E_I,:EP_I,FAST*P/KP),(:EP_I,:E_I,FAST),(:ESB_I,:EP_I,k),(:EP_I,:ESB_I,kr),
         (:E_A,:E_I,FAST*L),(:E_I,:E_A,FAST),(:EP_A,:EP_I,FAST*L),(:EP_I,:EP_A,FAST)],
        [(:ESB_A,:EP_A,k,kr),(:ESB_I,:EP_I,k,kr)],big(1.0))
end

println("="^88)
println("LDH i-state, REDESIGN frame (FULL I-graph + :OnlyA LIMIT oracle, koff drain survives)")
println("BIG-convergence check on draw 1 (NORM & SKIP rel-err vs limit oracle):")
println("="^88)
kon,koff,KB,KP,k,L,Keq,S,B,P = big(1.7),big(1.1),big(0.8),big(0.9),big(2.1),big(0.7),big(3.0),big(1.1),big(0.9),big(0.6)
kr=k*kon*KP/(koff*KB*Keq)
for BIG in (big"1e4",big"1e6",big"1e8",big"1e10",big"1e12")
    A,I=pieces(kon,koff,KB,KP,k,kr,S,B,P,BIG); o=oracle(kon,koff,KB,KP,k,kr,L,S,B,P,BIG)
    @printf "  BIG=%-7.0e  oracle=%-12.7f  NORMerr=%-10.2e  SKIPerr=%-10.2e\n" Float64(BIG) Float64(o) Float64(abs(vnorm(A,I,L,big(1))-o)/abs(o)) Float64(abs(vskip(A,I,L,big(1))-o)/abs(o))
end

BIG=big"1e10"
println("\nRandom draws (BIG=1e10):")
@printf "%-5s %-13s %-13s %-13s %-11s %-11s\n" "draw" "limit_oracle" "NORM" "SKIP" "NORMerr" "SKIPerr"
rng=MersenneTwister(1)
for t in 1:6
    kon=big(0.5+2rand(rng));koff=big(0.5+2rand(rng));KB=big(0.5+2rand(rng));KP=big(0.5+2rand(rng))
    k=big(0.5+2rand(rng));L=big(0.5+rand(rng));Keq=big(2+2rand(rng));S=big(0.5+2rand(rng));B=big(0.5+2rand(rng));P=big(0.5+2rand(rng))
    kr=k*kon*KP/(koff*KB*Keq)
    A,I=pieces(kon,koff,KB,KP,k,kr,S,B,P,BIG); o=oracle(kon,koff,KB,KP,k,kr,L,S,B,P,BIG)
    vn=vnorm(A,I,L,big(1)); vs=vskip(A,I,L,big(1))
    @printf "%-5d %-13.7f %-13.7f %-13.7f %-11.2e %-11.2e\n" t Float64(o) Float64(vn) Float64(vs) Float64(abs(vn-o)/abs(o)) Float64(abs(vs-o)/abs(o))
end

println("\nd_free_A, d_free_I (FULL near-limit graph) and ratio vs B  (is it metabolite-bearing?):")
kon,koff,KB,KP,k,L,Keq,S,P = big(1.7),big(1.1),big(0.8),big(0.9),big(2.1),big(0.7),big(3.0),big(1.1),big(0.6); kr=k*kon*KP/(koff*KB*Keq)
@printf "%-7s %-16s %-16s %-12s %-12s %-12s\n" "B" "dfree_A" "dfree_I" "dfA/dfI" "NORMerr" "SKIPerr"
rats=BigFloat[]
for Bf in (big(0.1),big(0.3),big(0.6),big(1.0),big(2.0),big(4.0),big(8.0))
    A,I=pieces(kon,koff,KB,KP,k,kr,S,Bf,P,BIG); o=oracle(kon,koff,KB,KP,k,kr,L,S,Bf,P,BIG)
    r=A.dfree/I.dfree; push!(rats,r)
    @printf "%-7.2f %-16.6g %-16.6g %-12.5f %-12.2e %-12.2e\n" Float64(Bf) Float64(A.dfree) Float64(I.dfree) Float64(r) Float64(abs(vnorm(A,I,L,big(1))-o)/abs(o)) Float64(abs(vskip(A,I,L,big(1))-o)/abs(o))
end
@printf "dfA/dfI span: %.5f .. %.5f  (%.3fx)\n" Float64(minimum(rats)) Float64(maximum(rats)) Float64(maximum(rats)/minimum(rats))

# sanity: L=0 must give active-only rate for BOTH combines and oracle
println("\nSanity L=0 (inactive unpopulated): oracle vs NORM vs SKIP must all equal active-only rate:")
kon,koff,KB,KP,k,Keq,S,B,P=big(1.7),big(1.1),big(0.8),big(0.9),big(2.1),big(3.0),big(1.1),big(0.9),big(0.6);kr=k*kon*KP/(koff*KB*Keq)
A,I=pieces(kon,koff,KB,KP,k,kr,S,B,P,BIG); o0=oracle(kon,koff,KB,KP,k,kr,big(0),S,B,P,BIG)
@printf "  oracle(L=0)=%.8f  NORM(L=0)=%.8f  SKIP(L=0)=%.8f  active-only(A.N/A.D)=%.8f\n" Float64(o0) Float64(vnorm(A,I,big(0),big(1))) Float64(vskip(A,I,big(0),big(1))) Float64(A.N/A.D)
