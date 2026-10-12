# ABOUTME: Toy checks of the P8 key (Canon): renaming invariance, normalization, symmetric
# ABOUTME: parameters that need search, distinct functions, and the common-factor diagnostic.
include(joinpath(@__DIR__, "Canon.jl"))
using Test, Random

P(ts...) = Canon.Poly(sort!([s => e for (s, e) in m]; by = first) => big(c) // 1 for (m, c) in ts)
law(num, den; params) = Canon.Law([num, den], Set([:A, :B, :P, :Keq]), [[1, 2]], params, "plain")

# Michaelis–Menten, reversible: v = (k1 k2 A − k3 k4 P / ...)… written over named constants.
mm(k1, k2, k3, k4) = law(P([k1 => 1, k2 => 1, :A => 1] => 1, [k3 => 1, k4 => 1, :P => 1] => -1),
                         P([k3 => 1] => 1, [k2 => 1] => 1, [k1 => 1, :A => 1] => 1,
                           [k4 => 1, :P => 1] => 1); params = [k1, k2, k3, k4])

@testset "renaming invariance" begin
    a = Canon.canonical(mm(:a, :b, :c, :d))
    b = Canon.canonical(mm(:w, :x, :y, :z))
    c = Canon.canonical(mm(:z, :y, :x, :w))
    @test a.key == b.key == c.key
    mp = Canon.name_map(a, c)
    @test mp == Dict(:a => :z, :b => :y, :c => :x, :d => :w)
    @test Canon.same_under(mp, mm(:a, :b, :c, :d), mm(:z, :y, :x, :w))
    @test a.leaves == 1
end

@testset "normalization: common monomial and constant factors" begin
    base = mm(:a, :b, :c, :d)
    scaled = Canon.Law([Canon.Poly(sort!([m; :a => -2]; by = first) => 3c for (m, c) in p)
                        for p in base.comps], base.fixed, base.units, base.params, base.header)
    # merging a's exponents: build the product monomial properly
    fix(p) = Canon.Poly(begin
        d = Dict{Symbol, Int}()
        for (s, e) in m; d[s] = get(d, s, 0) + e; end
        sort!([s => e for (s, e) in d if e != 0]; by = first)
    end => c for (m, c) in p)
    scaled = Canon.Law(fix.(scaled.comps), base.fixed, base.units, base.params, base.header)
    @test Canon.canonical(scaled).key == Canon.canonical(base).key
    neg = Canon.Law([Canon.Poly(m => -c for (m, c) in p) for p in base.comps], base.fixed,
                    base.units, base.params, base.header)
    @test Canon.canonical(neg).key == Canon.canonical(base).key
end

@testset "distinct functions get distinct keys" begin
    # swapping two parameters' roles in the denominator only changes the function
    other = law(P([:a => 1, :b => 1, :A => 1] => 1, [:c => 1, :d => 1, :P => 1] => -1),
                P([:c => 1] => 1, [:b => 1] => 1, [:b => 1, :A => 1] => 1, [:d => 1, :P => 1] => 1);
                params = [:a, :b, :c, :d])
    @test Canon.canonical(other).key != Canon.canonical(mm(:a, :b, :c, :d)).key
    # a different declared parameter count changes the key
    l5 = mm(:a, :b, :c, :d)
    l5 = Canon.Law(l5.comps, l5.fixed, l5.units, [l5.params; :e], l5.header)
    @test Canon.canonical(l5).key != Canon.canonical(mm(:a, :b, :c, :d)).key
end

@testset "symmetric parameters need search" begin
    # a and b appear only as the product a·b, c and d only as c·d: two symmetric pairs
    s = law(P([:a => 1, :b => 1, :A => 1] => 1, [:c => 1, :d => 1, :P => 1] => -1),
            P([] => 1, [:a => 1, :b => 1, :A => 1] => 1, [:c => 1, :d => 1, :P => 1] => 1);
            params = [:a, :b, :c, :d])
    cs = Canon.canonical(s)
    @test cs.leaves == 4 && cs.hits == 4
    s2 = Canon.rename(s, Dict(:a => :d2, :b => :a2, :c => :c2, :d => :b2))
    c2 = Canon.canonical(s2)
    @test c2.key == cs.key
    @test Canon.same_under(Canon.name_map(cs, c2), s, s2)
    # regular structure that 1-WL cannot split: a 6-cycle of parameters vs two triangles
    # (each term is one edge x·y; every parameter has degree 2)
    edges6 = [(:p1, :p2), (:p2, :p3), (:p3, :p4), (:p4, :p5), (:p5, :p6), (:p6, :p1)]
    edges33 = [(:p1, :p2), (:p2, :p3), (:p3, :p1), (:p4, :p5), (:p5, :p6), (:p6, :p4)]
    g(es) = law(P([:A => 1] => 1), P(([x => 1, y => 1] => 1 for (x, y) in es)...);
                params = [:p1, :p2, :p3, :p4, :p5, :p6])
    k6, k33 = Canon.canonical(g(edges6)), Canon.canonical(g(edges33))
    @test k6.key != k33.key
    @test k6.leaves > 1 && k33.leaves > 1
    perm = Dict(zip([:p1, :p2, :p3, :p4, :p5, :p6], shuffle(Xoshiro(3), [:p1, :p2, :p3, :p4, :p5, :p6])))
    @test Canon.canonical(Canon.rename(g(edges6), perm)).key == k6.key
end

@testset "common factor diagnostic" begin
    rng = Xoshiro(1)
    @test isempty(Canon.common_factor_vars(mm(:a, :b, :c, :d), rng))
    # multiply num and den by (a + A): a common factor in a and A
    base = mm(:a, :b, :c, :d)
    f = P([:a => 1] => 1, [:A => 1] => 1)
    mul(p, q) = (out = Canon.Poly();
                 for (m1, c1) in p, (m2, c2) in q
                     d = Dict{Symbol, Int}()
                     for (s, e) in [m1; m2]; d[s] = get(d, s, 0) + e; end
                     k = sort!([s => e for (s, e) in d if e != 0]; by = first)
                     out[k] = get(out, k, big(0) // 1) + c1 * c2
                 end; filter!(kv -> !iszero(kv.second), out))
    l2 = Canon.Law([mul(base.comps[1], f), mul(base.comps[2], f)], base.fixed, base.units,
                   base.params, base.header)
    @test sort(Canon.common_factor_vars(l2, rng)) == [:A, :a]
end

include(joinpath(@__DIR__, "Factor.jl"))
@testset "factor cancellation" begin
    base = mm(:a, :b, :c, :d)
    mul(p, q) = (out = Canon.Poly();
                 for (m1, c1) in p, (m2, c2) in q
                     d = Dict{Symbol, Int}()
                     for (s, e) in [m1; m2]; d[s] = get(d, s, 0) + e; end
                     k = sort!([s => e for (s, e) in d if e != 0]; by = first)
                     out[k] = get(out, k, big(0) // 1) + c1 * c2
                 end; filter!(kv -> !iszero(kv.second), out))
    for f in (P([:a => 1] => 1, [:A => 1] => 1),                                   # a + A
              P([:e => 1, :g => 1] => 1, [:A => 1, :g => 1] => 1, [:P => 1, :e => 1] => 1), # eg + Ag + Pe
              P([:e => 2] => 1, [:B => 1, :e => 1] => 2, [:B => 2] => 1))          # (e + B)^2
        l2 = Canon.Law([mul(base.comps[1], f), mul(base.comps[2], f)], base.fixed, base.units,
                       sort!(unique([base.params; :e; :g])), base.header)
        lc, removed, unres = Factor.cancel(l2)
        @test isempty(unres) && !isempty(removed)
        @test Canon.normalize(lc) == Canon.normalize(Canon.Law(base.comps, base.fixed, base.units,
                                                              l2.params, base.header))
        @test isempty(Canon.common_factor_vars(lc, Xoshiro(2)))
    end
end
