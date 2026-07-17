# ABOUTME: Confirms solve-then-limit (v_new) matches the independent full-6-form GT on EVERY
# ABOUTME: constructable uni-uni combo, and diverges only on the constructor-rejected ones.
using LinearAlgebra, Random, Printf

function solvenet(species, edges, cat_edges, Etot)
    n=length(species); idx=Dict(s=>i for (i,s) in enumerate(species)); M=zeros(n,n)
    for (a,b,r) in edges; M[idx[b],idx[a]]+=r; M[idx[a],idx[a]]-=r; end
    A=copy(M); A[1,:].=1.0; rhs=zeros(n); rhs[1]=Etot; c=A\rhs
    sum(kf*c[idx[r]]-kr*c[idx[p]] for (r,p,kf,kr) in cat_edges)
end

# independent full-6-form free-flip-1 network (no dropping of any form)
function full6(tags, p, concs; FAST=1e7)
    Stag,cattag,Ptag=tags; S,P=concs.S,concs.P
    K_I_S = Stag==:EqualAI ? p.K_A_S : p.K_I_S
    k_I   = cattag==:EqualAI ? p.k_A : p.k_I
    K_I_P = Ptag==:EqualAI ? p.K_A_P : p.K_I_P
    k_A_r = p.k_A*p.K_A_P/(p.K_A_S*p.Keq); k_I_r = k_I*K_I_P/(K_I_S*p.Keq)
    edges=Tuple{Symbol,Symbol,Float64}[
        (:E_A,:ES_A,FAST*S/p.K_A_S),(:ES_A,:E_A,FAST),
        (:E_A,:EP_A,FAST*P/p.K_A_P),(:EP_A,:E_A,FAST),
        (:ES_A,:EP_A,p.k_A),(:EP_A,:ES_A,k_A_r),
        (:E_A,:E_I,FAST*p.L),(:E_I,:E_A,FAST)]
    cat=Tuple{Symbol,Symbol,Float64,Float64}[(:ES_A,:EP_A,p.k_A,k_A_r)]
    Stag!=:OnlyA && push!(edges,(:E_I,:ES_I,FAST*S/K_I_S),(:ES_I,:E_I,FAST))
    Ptag!=:OnlyA && push!(edges,(:E_I,:EP_I,FAST*P/K_I_P),(:EP_I,:E_I,FAST))
    # I-catalysis runs only if ES_I or EP_I is reachable from the free pool (S or P binds in I);
    # if BOTH S and P are :OnlyA the I-catalytic forms are unreachable (inert inactive), zero population.
    if cattag!=:OnlyA && (Stag!=:OnlyA || Ptag!=:OnlyA)
        push!(edges,(:ES_I,:EP_I,k_I),(:EP_I,:ES_I,k_I_r)); push!(cat,(:ES_I,:EP_I,k_I,k_I_r))
    end
    # species = every node that appears in an edge, E_A first (carries conservation row)
    species=[:E_A]; for e in edges, s in (e[1],e[2]); s in species || push!(species,s); end
    solvenet(species,edges,cat,p.E_total)
end

function v_new(tags,p,concs)
    Stag,cattag,Ptag=tags; S,P=concs.S,concs.P
    K_I_S=Stag==:EqualAI ? p.K_A_S : p.K_I_S; k_I=cattag==:EqualAI ? p.k_A : p.k_I
    K_I_P=Ptag==:EqualAI ? p.K_A_P : p.K_I_P
    k_A_r=p.k_A*p.K_A_P/(p.K_A_S*p.Keq)
    N_A=p.k_A*S/p.K_A_S-k_A_r*P/p.K_A_P; D_A=1+S/p.K_A_S+P/p.K_A_P
    ES_I=Stag!=:OnlyA; EP_I=Ptag!=:OnlyA; Icat=cattag!=:OnlyA && ES_I && EP_I
    D_I=1.0; ES_I && (D_I+=S/K_I_S); EP_I && (D_I+=P/K_I_P)
    if Icat; k_I_r=k_I*K_I_P/(K_I_S*p.Keq); N_I=k_I*S/K_I_S-k_I_r*P/K_I_P; else; N_I=0.0; end
    p.E_total*(N_A+p.L*N_I)/(D_A+p.L*D_I)
end

draw(rng)=(K_A_S=0.5+2rand(rng),k_A=0.5+2rand(rng),K_A_P=0.5+2rand(rng),
           K_I_S=0.5+2rand(rng),k_I=0.5+2rand(rng),K_I_P=0.5+2rand(rng),
           L=0.5+rand(rng),Keq=2.0+2rand(rng),E_total=1.0)

# CONSTRUCTABLE set (confirmed via the package constructor in brief_enum_check.jl)
constructable = [
    (:OnlyA,:OnlyA,:EqualAI),      # Family A
    (:EqualAI,:OnlyA,:NonequalAI), # TR
    (:EqualAI,:EqualAI,:EqualAI),  # all-EqualAI
    (:EqualAI,:NonequalAI,:EqualAI),# NonequalAI cat / no-trap
    (:EqualAI,:OnlyA,:OnlyA),      # repair-1 (I-cat off)
    (:OnlyA,:EqualAI,:OnlyA),      # repair-2 (inert I)
]
# constructor-REJECTED set (trap + K-collapse)
rejected = [
    (:EqualAI,:NonequalAI,:OnlyA),
    (:NonequalAI,:NonequalAI,:OnlyA),
    (:EqualAI,:EqualAI,:OnlyA),
    (:OnlyA,:EqualAI,:EqualAI),
]

println("="^92)
println("solve-then-limit vs independent full-6-form GT, over all CONSTRUCTABLE uni-uni combos")
println("="^92)
worst = 0.0
for tags in constructable
    rng=MersenneTwister(hash(tags)); mr=0.0
    for _ in 1:20
        p=draw(rng); c=(S=0.5+2rand(rng),P=0.5+2rand(rng))
        mr=max(mr, abs(v_new(tags,p,c)-full6(tags,p,c))/max(abs(full6(tags,p,c)),1e-12))
    end
    global worst=max(worst,mr)
    @printf("  %-40s  max relerr = %.2e   %s\n", string(tags), mr, mr<1e-4 ? "MATCH" : "BREAK")
end
@printf("worst relerr over all constructable combos: %.2e\n\n", worst)

println("="^92)
println("same comparison over constructor-REJECTED combos (shown for contrast — these never reach derivation)")
println("="^92)
for tags in rejected
    rng=MersenneTwister(hash(tags)); mr=0.0
    for _ in 1:20
        p=draw(rng); c=(S=0.5+2rand(rng),P=0.5+2rand(rng))
        mr=max(mr, abs(v_new(tags,p,c)-full6(tags,p,c))/max(abs(full6(tags,p,c)),1e-12))
    end
    @printf("  %-40s  max relerr = %.2e   %s\n", string(tags), mr, mr<1e-4 ? "match" : "BREAK (constructor-rejected)")
end
