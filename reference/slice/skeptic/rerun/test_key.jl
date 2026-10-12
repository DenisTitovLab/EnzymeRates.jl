# ABOUTME: Toy tests of the K1 key (key.jl): renaming and reordering invariance, lowest terms,
# ABOUTME: distinct functions, symmetric constants that need search, rational exponents.
include(joinpath(@__DIR__, "lib.jl"))
using Test

"""A toy plain law over concentrations `A`, `P`: `H = A − P/Keq`, `Q` and `D` given as
`[(name => exponent, …) => coefficient]`, the fitted constants declared in `params` order."""
function _testhelper_toy(params, num, den)
    vars = S.Vars()
    foreach(x -> S.var!(vars, x; conc = true), (:A, :P))
    foreach(x -> S.var!(vars, x), params)
    S.var!(vars, :Keq)
    poly(ts) = foldl((p, (ms, c)) -> S.addto!(p, S.mono((vars.index[s] => e for (s, e) in ms)...),
                                              c), ts; init = S.Poly{Int}())
    D = poly(den)
    S.Law{Int}(vars, poly([((:A => 1,), 1), ((:P => 1, :Keq => -1), -1)]), poly(num), D,
               [D], [one(S.Poly{Int})], S.Mono())
end

_testhelper_mm(a, b, c, d; order = [a, b, c, d]) =
    _testhelper_toy(order, [((a => 1, b => 1), 1)],
                    [((c => 1,), 1), ((b => 1,), 1), ((a => 1, :A => 1), 1),
                     ((d => 1, :P => 1), 1)])

"""Whether `L2` is `L1` with each constant `o1[i]` renamed `o2[i]`: equal values at random
points."""
function _testhelper_same_under(L1, o1, L2, o2, rng)
    all(1:3) do _
        pars1 = Dict{Symbol, BigFloat}(n => Oracle.logu(rng) for n in L1.vars.names)
        pars2 = Dict{Symbol, BigFloat}(n => get(pars1, n, big(1)) for n in L2.vars.names)
        foreach(((x, y),) -> pars2[y] = pars1[x], zip(o1, o2))
        conc = Dict(:A => pars1[:A], :P => pars1[:P])
        v1 = law_value(L1, conc, (; (n => pars1[n] for n in L1.vars.names)...))
        v2 = law_value(L2, conc, (; (n => pars2[n] for n in L2.vars.names)...))
        abs(v1 - v2) <= 1e-60 * abs(v1)
    end
end

rng = Xoshiro(7)
@testset "renaming and reordering the fitted constants keep the key" begin
    a = S.canonical(_testhelper_mm(:a, :b, :c, :d))
    b = S.canonical(_testhelper_mm(:w, :x, :y, :z; order = [:z, :y, :x, :w]))
    c = S.canonical(_testhelper_mm(:z, :y, :x, :w))
    @test a.key == b.key == c.key
    @test Dict(zip(a.order, c.order)) == Dict(:a => :z, :b => :y, :c => :x, :d => :w)
end

@testset "lowest terms: a common monomial and an integer factor" begin
    L = _testhelper_mm(:a, :b, :c, :d)
    s = S.mono(L.vars.index[:a] => -2, L.vars.index[:A] => 1)
    L2 = S.Law{Int}(L.vars, L.H, S.scale(L.Q * s, 3), S.scale(L.D * s, 3), L.sigma, L.trees,
                    L.shift)
    @test S.canonical(L2).key == S.canonical(L).key
end

@testset "distinct functions get distinct keys" begin
    other = _testhelper_toy([:a, :b, :c, :d], [((:a => 1, :b => 1), 1)],
                            [((:c => 1,), 1), ((:b => 1,), 1), ((:b => 1, :A => 1), 1),
                             ((:d => 1, :P => 1), 1)])
    @test S.canonical(other).key != S.canonical(_testhelper_mm(:a, :b, :c, :d)).key
    absent = _testhelper_mm(:a, :b, :c, :d; order = [:a, :b, :c, :d, :e])
    ca = S.canonical(absent)
    @test ca.key != S.canonical(_testhelper_mm(:a, :b, :c, :d)).key && ca.order[end] == :e
end

@testset "symmetric constants need search" begin
    sym(a, b, c, d) = _testhelper_toy([a, b, c, d], [((a => 1, b => 1), 1)],
                                      [((), 1), ((a => 1, b => 1, :A => 1), 1),
                                       ((c => 1, d => 1, :P => 1), 1)])
    L1, L2 = sym(:a, :b, :c, :d), sym(:d2, :a2, :c2, :b2)
    c1, c2 = S.canonical(L1), S.canonical(L2)
    @test c1.key == c2.key
    @test _testhelper_same_under(L1, c1.order, L2, c2.order, rng)
    # every constant has degree 2 in both: colour refinement alone cannot tell a 6-cycle of
    # constants from two triangles
    ps = [:p1, :p2, :p3, :p4, :p5, :p6]
    graph(es, order = ps) = _testhelper_toy(order, [((), 1)],
                                            [((x => 1, y => 1), 1) for (x, y) in es])
    six = [(:p1, :p2), (:p2, :p3), (:p3, :p4), (:p4, :p5), (:p5, :p6), (:p6, :p1)]
    two = [(:p1, :p2), (:p2, :p3), (:p3, :p1), (:p4, :p5), (:p5, :p6), (:p6, :p4)]
    @test S.canonical(graph(six)).key != S.canonical(graph(two)).key
    @test S.canonical(graph(six, shuffle(rng, ps))).key == S.canonical(graph(six)).key
end

@testset "rational exponents" begin
    half(a, b, e) = _testhelper_toy([a, b], [((a => e, :Keq => -1 // 2), 1)],
                                    [((), 1), ((b => 1, :A => 1), 1), ((a => e, :P => 1), 1)])
    c1, c2 = S.canonical(half(:a, :b, 1 // 2)), S.canonical(half(:y, :x, 1 // 2))
    @test c1.key == c2.key
    @test _testhelper_same_under(half(:a, :b, 1 // 2), c1.order, half(:y, :x, 1 // 2),
                                 c2.order, rng)
    @test S.canonical(half(:a, :b, 1 // 1)).key != c1.key
    @test occursin("1//2", c1.text)
end

@testset "a derived law: renamed metabolites in their order keep the key" begin
    r = only(filter(r -> r.p.id == "init_bibi/55", records_of("init_bibi")))
    m = r.c.mech
    (subs, prods, steps, group), kw = S.raw(m)
    ren(x) = Symbol("z", x)
    met(x) = S.Met(ren(x.name), x.modifier)
    form(f) = S.RawForm(f.conf, met.(f.bound), ren.(f.added), ren.(f.subtracted))
    m2 = S.Mechanism(ren.(subs), ren.(prods),
                     [S.RawStep(form(s.from), form(s.to), met.(s.consumed), met.(s.released),
                                s.re) for s in steps], group; kw...)
    d, d2 = S.derive(m), S.derive(m2)
    @test d.key.key == d2.key.key
    pars = random_params(m, rng)
    pars2 = (; (y => getfield(pars, x) for (x, y) in zip(d.key.order, d2.key.order))...,
             Keq = pars.Keq, E_total = pars.E_total)
    conc = Dict(x => Oracle.logu(rng) for x in S.metabolites(m))
    v, v2 = law_value(d.law, conc, pars), law_value(d2.law, Dict(ren(x) => c for (x, c) in conc),
                                                    pars2)
    @test abs(v - v2) <= 1e-60 * abs(v)
end
