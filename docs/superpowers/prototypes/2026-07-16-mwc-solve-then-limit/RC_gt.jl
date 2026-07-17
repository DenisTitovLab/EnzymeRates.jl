# ABOUTME: n=1 two-conformation ground-truth for the two PFK ping-pong :OnlyA mechs.
# ABOUTME: Formulation 1 (only free E flips). Tests detailed balance at equilibrium.
using LinearAlgebra

# Generic steady-state solver: species names, directed edges (from,to,rate),
# conservation on species[1]. Species in a graph component DISCONNECTED from
# species[1] (the conserved enzyme pool) hold zero enzyme physically, so we
# solve only on species[1]'s connected component and return 0 for the rest.
# A disconnected reversible island is the structural signature we care about.
function solve(species, edges, Etot)
    # undirected adjacency for connectivity from species[1]
    adj = Dict(s=>Set{Symbol}() for s in species)
    for (a,b,_) in edges; push!(adj[a],b); push!(adj[b],a); end
    reach = Set([species[1]]); stack=[species[1]]
    while !isempty(stack)
        u=pop!(stack)
        for w in adj[u]; w in reach || (push!(reach,w); push!(stack,w)); end
    end
    isolated = setdiff(Set(species), reach)
    keep = [s for s in species if s in reach]
    kedges = [(a,b,r) for (a,b,r) in edges if a in reach && b in reach]
    n=length(keep); idx=Dict(s=>i for (i,s) in enumerate(keep))
    M=zeros(n,n)
    for (a,b,r) in kedges
        M[idx[b],idx[a]] += r
        M[idx[a],idx[a]] -= r
    end
    A=copy(M); A[1,:] .= 1.0
    rhs=zeros(n); rhs[1]=Etot
    c=A\rhs
    d=Dict(s=>c[idx[s]] for s in keep)
    for s in isolated; d[s]=0.0; end
    d[:__isolated__count__] = length(isolated)
    d
end

# Full two-conformation ping-pong network.
#   drop_in_I :: Set of group symbols dropped from the inactive conformation
#   (dropped = that group is :OnlyA, so no inactive-conformation edge).
function build(drop_in_I; K_ATP,K_ADP,K_F16BP,K_F6P,q4,k3r,L,Keq,
                          ATP,ADP,F16BP,F6P, FAST=1e7)
    k3f = k3r * (Keq*K_ATP*K_F6P/(K_F16BP*K_ADP)) / q4   # active-cycle Haldane

    EA,ATPa,F16a,COVa,F6a,ADPa = :EA,:E_ATP_A,:E_F16_A,:E_COV_A,:E_F6_A,:E_ADP_A
    EI,ATPi,F16i,COVi,F6i,ADPi = :EI,:E_ATP_I,:E_F16_I,:E_COV_I,:E_F6_I,:E_ADP_I

    edges = Tuple{Symbol,Symbol,Float64}[
        (EA,ATPa, FAST*ATP/K_ATP), (ATPa,EA, FAST),          # g2 bind ATP
        (ATPa,F16a, k3f), (F16a,ATPa, k3r),                  # g3 SS iso chem1
        (COVa,F16a, FAST*F16BP/K_F16BP), (F16a,COVa, FAST),  # g5 bind F16BP
        (COVa,F6a, FAST*F6P/K_F6P), (F6a,COVa, FAST),        # g6 bind F6P
        (F6a,ADPa, FAST*q4), (ADPa,F6a, FAST),               # g4 RE iso chem2
        (EA,ADPa, FAST*ADP/K_ADP), (ADPa,EA, FAST),          # g1 bind ADP
    ]
    species = Symbol[EA,ATPa,F16a,COVa,F6a,ADPa]

    Iforms = Set{Symbol}()
    function iedge!(g,a,b,rab,rba)
        g in drop_in_I && return
        push!(edges,(a,b,rab)); push!(edges,(b,a,rba))
        push!(Iforms,a); push!(Iforms,b)
    end
    iedge!(:g2, EI,ATPi, FAST*ATP/K_ATP, FAST)
    iedge!(:g3, ATPi,F16i, k3f, k3r)
    iedge!(:g5, COVi,F16i, FAST*F16BP/K_F16BP, FAST)
    iedge!(:g6, COVi,F6i, FAST*F6P/K_F6P, FAST)
    iedge!(:g4, F6i,ADPi, FAST*q4, FAST)
    iedge!(:g1, EI,ADPi, FAST*ADP/K_ADP, FAST)

    push!(edges,(EA,EI, FAST*L)); push!(edges,(EI,EA, FAST))   # only free E flips
    push!(species, EI)
    for f in Iforms; f==EI || push!(species,f); end

    c = solve(species, edges, 1.0)
    v_g3  = k3f*c[ATPa] - k3r*c[F16a]                 # turnover across active g3
    Jg1A  = FAST*ADP/K_ADP*c[EA] - FAST*c[ADPa]
    Jg1I  = haskey(c,ADPi) ? (FAST*ADP/K_ADP*c[EI] - FAST*c[ADPi]) : 0.0
    v_adp = -(Jg1A + Jg1I)                            # net ADP production
    (v_g3=v_g3, v_adp=v_adp, c=c)
end

function report(label, drop_in_I)
    K_ATP,K_ADP,K_F16BP,K_F6P = 0.7,1.3,0.9,1.1
    q4,k3r,L,Keq = 1.4,0.8,2.5,3.0
    ADP,F16BP,ATP = 0.6,0.5,1.1
    F6P_eq = ADP*F16BP/(Keq*ATP)                       # ATP*F6P/(ADP*F16BP)=1/Keq
    r_eq = build(drop_in_I; K_ATP,K_ADP,K_F16BP,K_F6P,q4,k3r,L,Keq,ATP,ADP,F16BP,F6P=F6P_eq)
    r_ne = build(drop_in_I; K_ATP,K_ADP,K_F16BP,K_F6P,q4,k3r,L,Keq,ATP,ADP,F16BP,F6P=5.0*F6P_eq)
    r_rv = build(drop_in_I; K_ATP,K_ADP,K_F16BP,K_F6P,q4,k3r,L,Keq,ATP,ADP,F16BP,F6P=0.1*F6P_eq)
    println("=== $label ===")
    println("  drop_in_I = $(sort(collect(drop_in_I)))")
    println("  at EQUILIBRIUM ratio:  v_g3=$(round(r_eq.v_g3,sigdigits=4))   v_adp=$(round(r_eq.v_adp,sigdigits=4))")
    println("  fwd-driven (F6P x5):   v_g3=$(round(r_ne.v_g3,sigdigits=4))   v_adp=$(round(r_ne.v_adp,sigdigits=4))")
    println("  rev-driven (F6P x0.1): v_g3=$(round(r_rv.v_g3,sigdigits=4))   v_adp=$(round(r_rv.v_adp,sigdigits=4))")
    ci=get(r_eq.c,:E_COV_I,NaN); f16i=get(r_eq.c,:E_F16_I,NaN); f6i=get(r_eq.c,:E_F6_I,NaN)
    println("  I-covalent masses @eq: E_COV_I=$(round(ci,sigdigits=3))  E_F16_I=$(round(f16i,sigdigits=3))  E_F6_I=$(round(f6i,sigdigits=3))")
    println("  # isolated (disconnected-from-E) forms in graph: $(Int(get(r_eq.c,:__isolated__count__,0)))")
end

report("err1  chem1(g3)=OnlyA, chem2(g4)=EqualAI, F6P(g6)=OnlyA", Set([:g3,:g6]))
report("err2  chem1(g3)=EqualAI, chem2(g4)=OnlyA, F6P(g6)=OnlyA", Set([:g4,:g6]))
