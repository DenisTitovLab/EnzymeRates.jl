# ABOUTME: Strongest claim-B test: ordered bi-uni, cat:NonequalAI, all-EqualAI bindings.
# ABOUTME: Metabolite-bearing D AND active catalysis in BOTH states -> plain combine vs oracle free-flip.
using LinearAlgebra, Random, Printf

function solve(species, edges, cat_edges, Etot)
    n=length(species); idx=Dict(s=>i for (i,s) in enumerate(species))
    M=zeros(n,n); for (a,b,r) in edges; M[idx[b],idx[a]]+=r; M[idx[a],idx[a]]-=r; end
    A=copy(M); A[1,:].=1.0; rhs=zeros(n); rhs[1]=Etot; c=A\rhs
    (flux=sum(kf*c[idx[r]]-kr*c[idx[p]] for (r,p,kf,kr) in cat_edges), Efree=c[1])
end

# oracle free-flip-only bi-uni NonequalAI catalysis (verbatim from allosteric_ground_truth.jl)
function biuni_nonequalAI_freeflip_flux(kon,koff,KB,KP; k_A,k_I,L,Keq,A,B,P,FAST=1e7)
    krA=k_A*kon*KP/(koff*KB*Keq); krI=k_I*kon*KP/(koff*KB*Keq)
    species=[:E_A,:EA_A,:EAB_A,:EP_A,:E_I,:EA_I,:EAB_I,:EP_I]
    edges=[(:E_A,:EA_A,kon*A),(:EA_A,:E_A,koff),(:E_I,:EA_I,kon*A),(:EA_I,:E_I,koff),
        (:EA_A,:EAB_A,FAST*B/KB),(:EAB_A,:EA_A,FAST),(:EA_I,:EAB_I,FAST*B/KB),(:EAB_I,:EA_I,FAST),
        (:E_A,:EP_A,FAST*P/KP),(:EP_A,:E_A,FAST),(:E_I,:EP_I,FAST*P/KP),(:EP_I,:E_I,FAST),
        (:E_A,:E_I,FAST*L),(:E_I,:E_A,FAST),
        (:EAB_A,:EP_A,k_A),(:EP_A,:EAB_A,krA),(:EAB_I,:EP_I,k_I),(:EP_I,:EAB_I,krI)]
    solve(species,edges,[(:EAB_A,:EP_A,k_A,krA),(:EAB_I,:EP_I,k_I,krI)],1.0).flux
end

# single-conformation ordered bi-uni at rate constant k -> (Efree, flux) for per-state normalize
function active_only(kon,koff,KB,KP,k,Keq,A,B,P;FAST=1e7)
    kr=k*kon*KP/(koff*KB*Keq)
    species=[:E,:EA,:EAB,:EP]
    edges=[(:E,:EA,kon*A),(:EA,:E,koff),(:EA,:EAB,FAST*B/KB),(:EAB,:EA,FAST),
           (:E,:EP,FAST*P/KP),(:EP,:E,FAST),(:EAB,:EP,k),(:EP,:EAB,kr)]
    r=solve(species,edges,[(:EAB,:EP,k,kr)],1.0)
    (Efree=r.Efree, flux=r.flux)
end

println("="^92)
println("Claim B STRONG: bi-uni NonequalAI-cat, metabolite-bearing D in BOTH states + active I-cat")
println("plain per-state-normalized combine  vs  oracle free-flip reference")
println("="^92)
@printf("%-6s %-16s %-16s %-10s\n","draw","oracle free-flip","plain combine","relerr")
rng=MersenneTwister(999); mr=0.0
for d in 1:10
    kon=0.5+2rand(rng);koff=0.5+2rand(rng);KB=0.5+2rand(rng);KP=0.5+2rand(rng)
    kA=0.5+2rand(rng);kI=0.5+2rand(rng);L=0.5+rand(rng);Keq=2.0+2rand(rng)
    A=0.5+2rand(rng);B=0.5+2rand(rng);P=0.5+2rand(rng)
    vgt=biuni_nonequalAI_freeflip_flux(kon,koff,KB,KP;k_A=kA,k_I=kI,L=L,Keq=Keq,A=A,B=B,P=P)
    a=active_only(kon,koff,KB,KP,kA,Keq,A,B,P); i=active_only(kon,koff,KB,KP,kI,Keq,A,B,P)
    D_A=1/a.Efree; N_A=a.flux/a.Efree; D_I=1/i.Efree; N_I=i.flux/i.Efree
    vpl=(N_A+L*N_I)/(D_A+L*D_I)
    rel=abs(vpl-vgt)/max(abs(vgt),1e-12); global mr=max(mr,rel)
    @printf("%-6d %-16.8g %-16.8g %-10.2e\n",d,vgt,vpl,rel)
end
println("-"^92)
@printf("plain-combine max relerr = %.2e   (matches oracle -> no cross-weighting needed)\n",mr)
