# ABOUTME: Scaling check for P6: rapid-equilibrium random-order mechanisms with n independent
# ABOUTME: substrate/product sites (3^n forms), written as Pop.PMech and timed as GMech.
include(joinpath(@__DIR__, "harness.jl"))
using Printf, Serialization

"""Sites `i = 1:n`, each empty or holding substrate `S_i` or product `P_i`; every binding is
RE and shares its metabolite's group; one SS chemistry step joins all-substrate to all-product."""
function sitemech(n)
    subs = [Symbol("S", i) for i in 1:n]; prods = [Symbol("P", i) for i in 1:n]
    states = vec(collect(Iterators.product(fill(0:2, n)...)))
    content(st) = sort([st[i] == 1 ? Pop.PMet(subs[i], :S) : Pop.PMet(prods[i], :P)
                        for i in 1:n if st[i] > 0]; by = x -> x.name)
    species = [Pop.PSpecies(content(st), :E, Symbol[], Symbol[],
                            Symbol("E", join(x.name for x in content(st)))) for st in states]
    idx = Dict(st => k for (k, st) in enumerate(states))
    steps = Pop.PStep[]; groups = [Int[] for _ in 1:2n + 1]
    for st in states, i in 1:n, v in 1:2
        st[i] == 0 || continue
        to = ntuple(j -> j == i ? v : st[j], n)
        push!(steps, Pop.PStep(idx[st], idx[to], [v == 1 ? Pop.PMet(subs[i], :S) : Pop.PMet(prods[i], :P)],
                               Pop.PMet[], true))
        push!(groups[2(i - 1) + v], length(steps))
    end
    push!(steps, Pop.PStep(idx[ntuple(_ -> 1, n)], idx[ntuple(_ -> 2, n)], Pop.PMet[], Pop.PMet[], false))
    push!(groups[end], length(steps))
    Pop.PMech(; id = "synthetic/$n", source = "synthetic: $n independent sites", substrates = subs,
              products = prods, metabolites = [subs; prods], species, steps, groups)
end

mechs = [sitemech(n) for n in 2:4]
serialize(joinpath(@__DIR__, "data", "synthetic.jls"), mechs)
for m in mechs
    raw = _testhelper_toraw(m)
    build() = GM.canonicalize(raw[1], raw[2], raw[3], raw[4]; raw[5]...)
    g = first(build())
    gv = _testhelper_gm(_testhelper_toraw(m; rng = Xoshiro(3), flip = bitrand(Xoshiro(4), length(m.steps))))
    a = GM.analyze(g)
    @printf("%-12s forms %3d edges %3d cycles %3d | ctor %7.0f ns (%d allocs) | analyze %7.0f ns | hash %5.0f ns | variant equal %s\n",
            m.id, length(g.forms), length(g.edges), length(a.cycles), _testhelper_time_ns(build, 10),
            @allocations(build()), _testhelper_time_ns(() -> GM.analyze(g), 10),
            _testhelper_time_ns(() -> hash(g), 100), gv == g && hash(gv) == hash(g))
end
