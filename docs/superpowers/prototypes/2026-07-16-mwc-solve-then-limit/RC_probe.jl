# ABOUTME: Cross-check hand graph vs package I-state pruning; reproduce NaN + kcat error.
using EnzymeRates, Random
const ER = EnzymeRates
const SP = "/tmp/claude-501/-home-denis-linux--julia-dev-EnzymeRates/fb96a560-d2b8-4694-a474-ff2f47968c52/scratchpad"

function load(path)
    m = Core.eval(ER, Meta.parse(strip(read(path,String))))()
    m
end

function dump_state(tag, m)
    am = ER.AllostericMechanism(m)
    println("\n########## $tag ##########")
    println("guard _onlya_haldane_violation => ",
        ER._onlya_haldane_violation(ER.reaction(am),
            collect(ER.steps(am)), collect(ER.cat_allo_states(am))))
    println("--- ACTIVE catalytic groups (tag) ---")
    for (g,grp) in enumerate(ER.steps(am))
        st = ER.cat_allo_state(am, g)
        for s in grp
            bd = ER.is_binding(s) ? "  bind $(ER.name(ER.bound_metabolite(s)))" : "  ISO"
            println("  g$g [$st]  ", ER.name(ER.from_species(s)), " -> ",
                ER.name(ER.to_species(s)), bd, "  ", ER.is_equilibrium(s) ? "RE" : "SS")
        end
    end
    # package's own I-state graph after :OnlyA drop + stranding prune
    sam = ER._state_allo_mechanism(am, :I)
    println("--- INACTIVE (:I) graph the DERIVATION builds ---")
    Iforms = Set{Symbol}()
    for (g,grp) in enumerate(ER.steps(sam)), s in grp
        push!(Iforms, ER.name(ER.from_species(s))); push!(Iforms, ER.name(ER.to_species(s)))
        println("  Ig$g [$(ER.cat_allo_state(sam,g))]  ", ER.name(ER.from_species(s)), " -> ",
            ER.name(ER.to_species(s)), ER.is_binding(s) ? "  bind $(ER.name(ER.bound_metabolite(s)))" : "  ISO")
    end
    # connectivity of the I-graph: components (undirected)
    edges = Tuple{Symbol,Symbol}[]
    for grp in ER.steps(sam), s in grp
        push!(edges, (ER.name(ER.from_species(s)), ER.name(ER.to_species(s))))
    end
    comps = components(collect(Iforms), edges)
    println("  I-graph forms: ", sort(collect(Iforms)))
    println("  I-graph connected components: ", length(comps))
    for (i,c) in enumerate(comps); println("    comp$i: ", sort(collect(c))); end
    am
end

function components(nodes, edges)
    adj = Dict(n=>Set{Symbol}() for n in nodes)
    for (a,b) in edges; push!(adj[a],b); push!(adj[b],a); end
    seen=Set{Symbol}(); comps=Vector{Set{Symbol}}()
    for n in nodes
        n in seen && continue
        c=Set{Symbol}(); st=[n]; push!(seen,n)
        while !isempty(st)
            u=pop!(st); push!(c,u)
            for w in adj[u]; w in seen||(push!(seen,w);push!(st,w)); end
        end
        push!(comps,c)
    end
    comps
end

function repro_symptoms(tag, m)
    println("\n--- $tag: derivation symptoms ---")
    rng = MersenneTwister(1)
    fp = ER.fitted_params(m)
    println("  fitted_params = ", fp)
    Keq = 3.0
    K_ATP,K_ADP,K_F16,K_F6P = 0.7,1.3,0.9,1.1
    ADP,F16BP,ATP = 0.6,0.5,1.1
    F6P = ADP*F16BP/(Keq*ATP)
    concs = (ATP=ATP, F6P=F6P, ADP=ADP, F16BP=F16BP,
             Citrate=0.0, F26BP=0.0, Phosphate=0.0)
    for trial in 1:3
        vals = Dict(s => 0.5+2rand(rng) for s in fp)
        prm = NamedTuple{(fp..., :Keq, :E_total)}(((vals[s] for s in fp)..., Keq, 1.0))
        v = try real(ER.rate_equation(m, concs, ER.Reduced, prm)) catch e; "ERR:$e"; end
        println("  trial$trial rate_equation@eq = ", v)
    end
    kc = try ER._kcat_forward(m, NamedTuple{(fp..., :Keq, :E_total)}(
            ((0.5+2rand(rng) for s in fp)..., Keq, 1.0)))
         catch e; "ERR: $(sprint(showerror,e))"; end
    println("  _kcat_forward = ", kc)
end

m1 = load("$SP/pfkp_err1.txt")
m2 = load("$SP/pfkp_err2.txt")
dump_state("ERR1  (g3 chem1 :OnlyA)", m1)
dump_state("ERR2  (g4 chem2 :OnlyA)", m2)
repro_symptoms("ERR1", m1)
repro_symptoms("ERR2", m2)
