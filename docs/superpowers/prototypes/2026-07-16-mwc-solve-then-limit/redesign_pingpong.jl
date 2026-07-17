# ABOUTME: Ping-pong bi-bi (:NonequalAI cat) in the redesign frame: full-graph pieces,
# ABOUTME: free-flip oracle. d_free_A vs d_free_I differ here (cat rates enter D[:E]). NORM vs SKIP?
using LinearAlgebra, Printf, Random
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

# one conformation's King-Altman pieces (free ref = E)
function conf(p, kc1,kc1r,kc2,kc2r, A,B,P,Q)
    ND([:E,:EA,:F,:FB],
       [(:E,:EA,p.konA*A),(:EA,:E,p.koffA),(:EA,:F,kc1),(:F,:EA,kc1r*P),
        (:F,:FB,p.konB*B),(:FB,:F,p.koffB),(:FB,:E,kc2),(:E,:FB,kc2r*Q)],
       [(:EA,:F,kc1,kc1r*P),(:FB,:E,kc2,kc2r*Q)],:E)
end
# free-flip oracle: only E flips (E_A<->E_I ratio L); F does NOT flip (formulation 1)
function oracle_Eflip(p, L,A,B,P,Q)
    flux([:E_A,:EA_A,:F_A,:FB_A,:E_I,:EA_I,:F_I,:FB_I],
      [(:E_A,:EA_A,p.konA*A),(:EA_A,:E_A,p.koffA),(:EA_A,:F_A,p.kc1A),(:F_A,:EA_A,p.kc1rA*P),
       (:F_A,:FB_A,p.konB*B),(:FB_A,:F_A,p.koffB),(:FB_A,:E_A,p.kc2A),(:E_A,:FB_A,p.kc2rA*Q),
       (:E_I,:EA_I,p.konA*A),(:EA_I,:E_I,p.koffA),(:EA_I,:F_I,p.kc1I),(:F_I,:EA_I,p.kc1rI*P),
       (:F_I,:FB_I,p.konB*B),(:FB_I,:F_I,p.koffB),(:FB_I,:E_I,p.kc2I),(:E_I,:FB_I,p.kc2rI*Q),
       (:E_A,:E_I,big"1e7"*L),(:E_I,:E_A,big"1e7")],
      [(:EA_A,:F_A,p.kc1A,p.kc1rA*P),(:FB_A,:E_A,p.kc2A,p.kc2rA*Q),
       (:EA_I,:F_I,p.kc1I,p.kc1rI*P),(:FB_I,:E_I,p.kc2I,p.kc2rI*Q)],big(1.0))
end

bf(x)=big(x)
println("="^90)
println("PING-PONG bi-bi :NonequalAI, REDESIGN frame (full-graph pieces, E-only free-flip oracle)")
println("="^90)
@printf "%-5s %-13s %-13s %-13s %-11s %-11s %-9s\n" "draw" "oracle" "NORM" "SKIP" "NORMerr" "SKIPerr" "dfA/dfI"
rng=MersenneTwister(7)
for t in 1:8
    p=(konA=bf(.5+2rand(rng)),koffA=bf(.5+2rand(rng)),konB=bf(.5+2rand(rng)),koffB=bf(.5+2rand(rng)),
       kc1A=bf(.5+2rand(rng)),kc1rA=bf(.3+rand(rng)),kc2A=bf(.5+2rand(rng)),kc2rA=bf(.3+rand(rng)),
       kc1I=bf(.5+2rand(rng)),kc1rI=bf(.3+rand(rng)),kc2I=bf(.5+2rand(rng)),kc2rI=bf(.3+rand(rng)))
    L=bf(.5+rand(rng));A=bf(.5+2rand(rng));B=bf(.5+2rand(rng));P=bf(.3+rand(rng));Q=bf(.3+rand(rng))
    Ap=conf(p,p.kc1A,p.kc1rA,p.kc2A,p.kc2rA,A,B,P,Q)
    Ip=conf(p,p.kc1I,p.kc1rI,p.kc2I,p.kc2rI,A,B,P,Q)
    o=oracle_Eflip(p,L,A,B,P,Q); vn=vnorm(Ap,Ip,L,big(1)); vs=vskip(Ap,Ip,L,big(1))
    @printf "%-5d %-13.7f %-13.7f %-13.7f %-11.2e %-11.2e %-9.4f\n" t Float64(o) Float64(vn) Float64(vs) Float64(abs(vn-o)/abs(o)) Float64(abs(vs-o)/abs(o)) Float64(Ap.dfree/Ip.dfree)
end

# Sweep B: is dfA/dfI concentration-dependent? does SKIP error track it?
println("\nSweep B (fixed rates): d_free ratio and SKIP error vs B")
p=(konA=bf(1.3),koffA=bf(1.1),konB=bf(0.9),koffB=bf(1.2),kc1A=bf(2.1),kc1rA=bf(0.6),kc2A=bf(1.7),kc2rA=bf(0.5),
   kc1I=bf(0.7),kc1rI=bf(0.4),kc2I=bf(0.6),kc2rI=bf(0.3))
L=bf(0.8);A=bf(1.1);P=bf(0.6);Q=bf(0.5)
@printf "%-7s %-14s %-14s %-11s %-11s\n" "B" "dfA/dfI" "reqd L'=L*dfI/dfA" "NORMerr" "SKIPerr"
rats=BigFloat[]
for B in (bf(0.1),bf(0.3),bf(0.6),bf(1.0),bf(2.0),bf(4.0),bf(8.0))
    Ap=conf(p,p.kc1A,p.kc1rA,p.kc2A,p.kc2rA,A,B,P,Q); Ip=conf(p,p.kc1I,p.kc1rI,p.kc2I,p.kc2rI,A,B,P,Q)
    o=oracle_Eflip(p,L,A,B,P,Q); r=Ap.dfree/Ip.dfree; push!(rats,r)
    @printf "%-7.2f %-14.6f %-14.6f %-11.2e %-11.2e\n" Float64(B) Float64(r) Float64(L*Ip.dfree/Ap.dfree) Float64(abs(vnorm(Ap,Ip,L,big(1))-o)/abs(o)) Float64(abs(vskip(Ap,Ip,L,big(1))-o)/abs(o))
end
@printf "dfA/dfI span: %.5f .. %.5f  (%.3fx)\n" Float64(minimum(rats)) Float64(maximum(rats)) Float64(maximum(rats)/minimum(rats))
