# ABOUTME: Tests of the dedup key (canonical labelling, lowest terms, the leaf cap, which
# ABOUTME: metabolite binds, five pinned spec keys) and the L-phantom rescaling certificate.

# ── Toy laws over named constants ───────────────────────────────────────────

"""The polynomial Σ c·Π x^e of `terms`, each `c => [x => e, …]`, its variables looked up
by name in `vars`."""
_testhelper_poly(vars::EnzymeRates._Vars, terms) =
    foldl(+, (EnzymeRates._term(EnzymeRates._monomial(
                  (vars.index[x] => e for (x, e) in xs)...), c) for (c, xs) in terms);
          init = zero(EnzymeRates._Poly))

"""A plain law over the concentrations `concs`, the constants `constants` and `:Keq`,
with `H`, `Q` and `D` written as `_testhelper_poly` terms; the key reads no other field."""
function _testhelper_toy_law(concs, constants, H, Q, D)
    vars = EnzymeRates._Vars()
    foreach(x -> EnzymeRates._var!(vars, x; conc = true), concs)
    foreach(x -> EnzymeRates._var!(vars, x), constants)
    EnzymeRates._var!(vars, :Keq)
    p(ts) = _testhelper_poly(vars, ts)
    EnzymeRates._Law(vars, p(H), p(Q), p(D), [p(D)], [one(EnzymeRates._Poly)],
                     EnzymeRates._Mono())
end

"""Reversible Michaelis–Menten over the constants `k1`, `k2`, `k3`, `k4` (and `extra`,
which no term holds): `H = A − P/Keq`, `Q = k1·k2`, `D = k3 + k2 + k1·A + k4·P`."""
_testhelper_mm_law(k1, k2, k3, k4; extra = Symbol[]) = _testhelper_toy_law(
    [:A, :P], [k1, k2, k3, k4, extra...],
    [1 => [:A => 1], -1 => [:P => 1, :Keq => -1]], [1 => [k1 => 1, k2 => 1]],
    [1 => [k3 => 1], 1 => [k2 => 1], 1 => [k1 => 1, :A => 1], 1 => [k4 => 1, :P => 1]])

"""A law whose constants come in two symmetric pairs: `a`, `b` only as `a·b` and `c`, `d`
only as `c·d`."""
_testhelper_pairs_law(a, b, c, d) = _testhelper_toy_law(
    [:A, :P], [a, b, c, d], [1 => [:A => 1], -1 => [:P => 1, :Keq => -1]],
    [1 => Pair{Symbol, Int}[]],
    [1 => Pair{Symbol, Int}[], 1 => [a => 1, b => 1, :A => 1],
     1 => [c => 1, d => 1, :P => 1]])

"""A law whose denominator holds one term `x·y` per edge `(x, y)` of a graph on the six
constants p1…p6: every constant has degree 2, so refinement alone cannot split them."""
_testhelper_graph_law(edges) = _testhelper_toy_law(
    [:A, :P], [:p1, :p2, :p3, :p4, :p5, :p6],
    [1 => [:A => 1], -1 => [:P => 1, :Keq => -1]], [1 => Pair{Symbol, Int}[]],
    [1 => [x => 1, y => 1] for (x, y) in edges])

"""A law over the concentration S and nine constants a1…a9 that multiply S alike in
`D = 1 + Σ S·a_i`: 9! labellings tie."""
_testhelper_nine_law() = _testhelper_toy_law(
    [:S], [Symbol(:a, i) for i in 1:9], [1 => [:S => 1]], [1 => Pair{Symbol, Int}[]],
    [1 => Pair{Symbol, Int}[]; [1 => [Symbol(:a, i) => 1, :S => 1] for i in 1:9]])

const _testhelper_cycle6 = [(:p1, :p2), (:p2, :p3), (:p3, :p4), (:p4, :p5), (:p5, :p6),
                            (:p6, :p1)]

@testset "renaming fitted constants keeps the key and gives the induced order" begin
    a = EnzymeRates._canonical(_testhelper_mm_law(:a, :b, :c, :d), "toy")
    b = EnzymeRates._canonical(_testhelper_mm_law(:w, :x, :y, :z), "toy")
    c = EnzymeRates._canonical(_testhelper_mm_law(:z, :y, :x, :w), "toy")
    @test a.key == b.key == c.key
    @test length(a.key) == 64 && all(in("0123456789abcdef"), a.key)
    @test Dict(zip(a.order, c.order)) == Dict(:a => :z, :b => :y, :c => :x, :d => :w)
    @test a.leaves == 1 && !a.capped
    # Random renamings of three toys, each a permutation of the constant names.
    rng = Random.MersenneTwister(94)
    for L in (_testhelper_mm_law(:a, :b, :c, :d), _testhelper_pairs_law(:a, :b, :c, :d),
              _testhelper_graph_law(_testhelper_cycle6))
        names = filter(!=(:Keq), L.vars.names[.!L.vars.conc])
        cL = EnzymeRates._canonical(L, "toy")
        for _ in 1:3
            R = _testhelper_renamed_law(L, Dict(zip(names, shuffle(rng, names))))
            cR = EnzymeRates._canonical(R, "toy")
            @test cR.key == cL.key
            @test _testhelper_same_law_under_map(L, R, cL.order, cR.order; rng)
        end
    end
    # Swapping two constants' roles in the denominator changes the function and the key.
    swapped = _testhelper_toy_law([:A, :P], [:a, :b, :c, :d],
        [1 => [:A => 1], -1 => [:P => 1, :Keq => -1]], [1 => [:a => 1, :b => 1]],
        [1 => [:c => 1], 1 => [:b => 1], 1 => [:b => 1, :A => 1], 1 => [:d => 1, :P => 1]])
    @test EnzymeRates._canonical(swapped, "toy").key != a.key
end

@testset "the key brings each unit to lowest terms" begin
    # Q and D times 3·a⁻²·A, and H times 2: each unit's lowest-terms form, and so the key,
    # stay the same.
    base = _testhelper_mm_law(:a, :b, :c, :d)
    v = base.vars
    m = EnzymeRates._monomial(v.index[:a] => -2, v.index[:A] => 1)
    scaled = EnzymeRates._Law(v, EnzymeRates._scale(base.H, 2),
                              EnzymeRates._scale(base.Q * m, 3),
                              EnzymeRates._scale(base.D * m, 3), base.sigma, base.trees,
                              base.shift)
    @test EnzymeRates._canonical(scaled, "toy").key ==
          EnzymeRates._canonical(base, "toy").key
end

@testset "the key brings each MWC unit to lowest terms" begin
    # Each unit times its own integer and monomial: H by 2, (Q_A, D_A, fA) by 3·k⁻²·c,
    # (Q_I, D_I, fI) by 5·k and the two sites' (W_A, W_I) by 7·c² and 11·c³, with c a
    # concentration and k a constant. The rate doubles and the key stays the same. Scaling
    # fA alone, fI alone or one site's W_A alone changes the rate law and the key.
    rng = Random.MersenneTwister(219)
    for (name, nsites) in ("MWC Dimer [AllostericEnzymeMechanism]" => 0,
                           "MWC Tetramer Random Bi-Bi RE + Two Allosteric Sites" => 2)
        @testset "$name" begin
            A = EnzymeRates._law(_testhelper_lift(
                only(s for s in MECHANISM_TEST_SPECS if s.name == name).mechanism))
            v = A.vars
            c = findfirst(v.conc)
            k = findfirst(i -> !v.conc[i] && v.names[i] ∉ (:Keq, :L), eachindex(v.names))
            ma = EnzymeRates._monomial(k => -2, c => 1)
            mi = EnzymeRates._monomial(k => 1)
            ms = [EnzymeRates._monomial(c => i + 1) for i in 1:2]
            ns = (7, 11)
            sc(p, mo, n) = EnzymeRates._scale(p * mo, n)
            key(L) = EnzymeRates._canonical(L, "x").key
            function _testhelper_scaled_conformation(C::EnzymeRates._Law, mo, n)
                EnzymeRates._Law(C.vars, C.H, sc(C.Q, mo, n), sc(C.D, mo, n), C.sigma,
                                 C.trees, C.shift)
            end
            function _testhelper_mwc(; H = A.H, active = A.active, inactive = A.inactive,
                                     fA = A.fA, fI = A.fI, sites = A.sites)
                EnzymeRates._AlloLaw(v, H, active, inactive, fA, fI, sites, A.n)
            end
            scaled = _testhelper_mwc(;
                H = EnzymeRates._scale(A.H, 2),
                active = _testhelper_scaled_conformation(A.active, ma, 3),
                inactive = _testhelper_scaled_conformation(A.inactive, mi, 5),
                fA = sc(A.fA, ma, 3), fI = sc(A.fI, mi, 5),
                sites = [(sc(WA, ms[i], ns[i]), sc(WI, ms[i], ns[i]), n)
                         for (i, (WA, WI, n)) in enumerate(A.sites)])
            env = Dict{Symbol, BigFloat}(x => exp(2 * big(rand(rng)) - 1) for x in v.names)
            rate(L) = _testhelper_evaluate_law(L, env)
            @test abs(rate(scaled) / rate(A) - 2) <= big(10)^-20
            @test key(scaled) == key(A)
            @test key(_testhelper_mwc(; fA = sc(A.fA, ma, 3))) != key(A)
            @test key(_testhelper_mwc(; fI = sc(A.fI, mi, 5))) != key(A)
            @test length(A.sites) == nsites
            for i in eachindex(A.sites)
                WA, WI, n = A.sites[i]
                sites = [A.sites[1:(i - 1)]; (sc(WA, ms[i], ns[i]), WI, n);
                         A.sites[(i + 1):end]]
                @test key(_testhelper_mwc(; sites)) != key(A)
            end
        end
    end
end

@testset "symmetric constants take 4 leaves and one key" begin
    # a and b tie, and c and d: two branches for each pair, four labellings in all.
    s = _testhelper_pairs_law(:a, :b, :c, :d)
    cs = EnzymeRates._canonical(s, "toy")
    @test cs.leaves == 4 && !cs.capped
    s2 = _testhelper_renamed_law(s, Dict(:a => :d, :b => :a, :c => :c, :d => :b))
    c2 = EnzymeRates._canonical(s2, "toy")
    @test c2.key == cs.key
    @test _testhelper_same_law_under_map(s, s2, cs.order, c2.order;
                                         rng = Random.MersenneTwister(95))
end

@testset "a 6-cycle and two triangles get different keys" begin
    # Both graphs are 2-regular, so refinement alone ties every constant; individualizing
    # tells them apart, and a relabelled 6-cycle keeps its key.
    cycle = _testhelper_graph_law(_testhelper_cycle6)
    triangles = _testhelper_graph_law([(:p1, :p2), (:p2, :p3), (:p3, :p1), (:p4, :p5),
                                       (:p5, :p6), (:p6, :p4)])
    k6 = EnzymeRates._canonical(cycle, "toy")
    k33 = EnzymeRates._canonical(triangles, "toy")
    @test k6.key != k33.key
    @test k6.leaves > 1 && k33.leaves > 1 && !k6.capped && !k33.capped
    names = [:p1, :p2, :p3, :p4, :p5, :p6]
    rename = Dict(zip(names, shuffle(Random.MersenneTwister(96), names)))
    @test EnzymeRates._canonical(_testhelper_renamed_law(cycle, rename), "toy").key ==
          k6.key
end

@testset "an extra declared constant changes the key" begin
    plain = EnzymeRates._canonical(_testhelper_mm_law(:a, :b, :c, :d), "toy")
    extra = EnzymeRates._canonical(_testhelper_mm_law(:a, :b, :c, :d; extra = [:e]), "toy")
    @test extra.key != plain.key
    @test extra.order[end] === :e
    @test occursin("|P=5|absent=1", extra.text) && occursin("|P=4|absent=0", plain.text)
end

@testset "which metabolite binds changes the key" begin
    # Concentrations are labelled by name, so models that differ only in which inhibitor
    # or regulator binds, or where an inhibitor binds, get different keys.
    key(m) = EnzymeRates._canonical(EnzymeRates._law(_testhelper_lift(m)), "attack").key
    uni_I = @enzyme_mechanism begin
        substrates: S
        products: P
        modifiers: I
        steps: begin
            E + S ⇌ E(S)
            E(S) <--> E + P
            E + I ⇌ E(I)
        end
    end
    uni_J = @enzyme_mechanism begin
        substrates: S
        products: P
        modifiers: J
        steps: begin
            E + S ⇌ E(S)
            E(S) <--> E + P
            E + J ⇌ E(J)
        end
    end
    unc_I = @enzyme_mechanism begin
        substrates: S
        products: P
        modifiers: I
        steps: begin
            E + S ⇌ E(S)
            E(S) <--> E + P
            E(S) + I ⇌ E(S, I)
        end
    end
    allo_R = @allosteric_mechanism begin
        substrates: S
        products: P
        allosteric_regulators: R
        catalytic_multiplicity: 2
        steps: begin
            E + S ⇌ E(S)       :: NonequalAI
            E(S) <--> E + P    :: NonequalAI
        end
        regulatory_site(1): R::OnlyA
    end
    allo_T = @allosteric_mechanism begin
        substrates: S
        products: P
        allosteric_regulators: T
        catalytic_multiplicity: 2
        steps: begin
            E + S ⇌ E(S)       :: NonequalAI
            E(S) <--> E + P    :: NonequalAI
        end
        regulatory_site(1): T::OnlyA
    end
    @test key(uni_I) != key(uni_J)
    @test key(uni_I) != key(unc_I)
    @test key(allo_R) != key(allo_T)
end

@testset "at 100 leaves a law gets its mechanism's key" begin
    # Nine constants that tie: the search stops at the cap, and the law takes the key of
    # its mechanism's id, so twins of it are fitted apart.
    id = "1234567e89012345"
    c = EnzymeRates._canonical(_testhelper_nine_law(), id)
    @test c.capped && c.leaves == EnzymeRates._KEY_LEAF_CAP == 100
    @test c.text == "mechanism|" * id
    @test c.key == bytes2hex(EnzymeRates.sha256("mechanism|" * id))
    @test c.order == [Symbol(:a, i) for i in 1:9]
end

@testset "five spec keys are pinned" begin
    # A key is the same in every process and Julia version (spec section 2, criterion 3).
    # These digests were read once and checked on Julia 1.12.5 and 1.13.1.
    pinned = ("Uni-Uni" =>
                  "e654c0e40f4ecd1af22ade6f28d2307f9bb15b90e6a409efbdf2ef8fb3ec0fb5",
              "Segel Ordered Bi Bi" =>
                  "8709d323651f7c90d4139f4b0f459c3275d245cfc6e0157510a794a82e3c0a39",
              "Random-order Bi-Bi" =>
                  "373b0fe0beb850947b7a3401aca91854b7618ba93f9af8518d346379092dac39",
              "Segel Ping Pong Bi Bi" =>
                  "d981644bde31a60073ec02634d1eb0ac94b81e88e228ff61d97d318c3d965fac",
              "MWC Dimer [AllostericEnzymeMechanism]" =>
                  "8b54530459132227ea8fe854cce576346e7ae872ac59f137eab1f6a5d913cbbf")
    for (name, key) in pinned
        spec = only(s for s in MECHANISM_TEST_SPECS if s.name == name)
        m = _testhelper_lift(spec.mechanism)
        @test EnzymeRates._canonical(EnzymeRates._law(m), "pin").key == key
    end
end

# ── One-subunit L phantoms ─────────────────────────────────────────────────

@testset "a certificate maps the phantom onto its L→0 law" begin
    # The merged uni-uni K-type at one catalytic subunit: E + S <--> E(P) and
    # E + P ⇌ E(P), both :OnlyA. Its inactive conformation is free E alone, so the law is
    # v = N/(Y + L·X) with X a term of Y. The certificate scales P's dissociation constant
    # by (1 + L) and the uptake rate constant by 1/(1 + L): at every L the variant is its
    # law at L = 0, the law of its catalytic mechanism, at those constants (spec
    # section 13).
    am = EnzymeRates.AllostericMechanism(@allosteric_mechanism begin
        substrates: S
        products: P
        steps: begin
            E + S <--> E(P)    :: OnlyA
            E + P ⇌ E(P)       :: OnlyA
        end
    end)
    A = EnzymeRates._law(am)
    θ = [A.vars.names[v] for v in eachindex(A.vars.names)
         if !A.vars.conc[v] && !(A.vars.names[v] in (:Keq, :L))]
    @test θ == [:K_A_EP_to_E_P, :k_A_E_S_to_EP]
    @test EnzymeRates._phantom_weights(A) == [1 // 1, -1 // 1]
    @test EnzymeRates._is_phantom(am, A)
    @test _testhelper_phantom_map_holds(am; rng = Random.MersenneTwister(21))
end

@testset "no phantom without the pre-filter's conditions" begin
    # The spec's pre-filter stays, and the certificate runs inside it. The merged
    # uni-uni K-type at two subunits, with P's binding :NonequalAI (the inactive
    # conformation binds P), and with an :OnlyI site ligand fails the pre-filter and has
    # no certificate. With an :OnlyA site ligand it passes the pre-filter, but the
    # regulator reveals L and no certificate exists. The Theorell–Chance K-type only_a
    # has a certificate, yet its inactive conformation binds Q at free E, so it is fitted.
    two_subunits = EnzymeRates.AllostericMechanism(@allosteric_mechanism begin
        substrates: S
        products: P
        catalytic_multiplicity: 2
        steps: begin
            E + S <--> E(P)    :: OnlyA
            E + P ⇌ E(P)       :: OnlyA
        end
    end)
    p_relaxed = EnzymeRates.AllostericMechanism(@allosteric_mechanism begin
        substrates: S
        products: P
        catalytic_multiplicity: 1
        steps: begin
            E + S <--> E(P)    :: OnlyA
            E + P ⇌ E(P)       :: NonequalAI
        end
    end)
    onlyi_site = EnzymeRates.AllostericMechanism(@allosteric_mechanism begin
        substrates: S
        products: P
        allosteric_regulators: R
        steps: begin
            E + S <--> E(P)    :: OnlyA
            E + P ⇌ E(P)       :: OnlyA
        end
        regulatory_site(1): R::OnlyI
    end)
    onlya_site = EnzymeRates.AllostericMechanism(@allosteric_mechanism begin
        substrates: S
        products: P
        allosteric_regulators: R
        steps: begin
            E + S <--> E(P)    :: OnlyA
            E + P ⇌ E(P)       :: OnlyA
        end
        regulatory_site(1): R::OnlyA
    end)
    only_a = EnzymeRates.AllostericMechanism(@allosteric_mechanism begin
        substrates: A, B ; products: P, Q ; catalytic_multiplicity: 1
        steps: begin
            E + A ⇌ E(A)                :: OnlyA
            E(A) + B <--> E(Q) + P      :: OnlyA
            E + Q ⇌ E(Q)                :: EqualAI
        end
    end)
    for am in (two_subunits, p_relaxed, onlyi_site, onlya_site, only_a)
        @test !EnzymeRates._is_phantom(am, EnzymeRates._law(am))
    end
    for am in (two_subunits, p_relaxed, onlyi_site, onlya_site)
        @test isempty(EnzymeRates._phantom_weights(EnzymeRates._law(am)))
    end
    @test !any(EnzymeRates._inactive_mask(onlya_site))
    @test !isempty(EnzymeRates._phantom_weights(EnzymeRates._law(only_a)))
    @test any(EnzymeRates._inactive_mask(only_a))
end

@testset "a certificate's rows overflow to no certificate" begin
    # A one-subunit toy law v = A·k1^p / (1 + k1^p·A + k1·k2^q·A + L), its inactive
    # conformation free E alone. Its certificate weighs k1 by −1/p and k2 by
    # (1 − p)/(p·q): −1/2 and −1/6 for p = 2, q = 3. For p = 2^40 + 15 and q = 2^40 + 3 the
    # second weight's denominator leaves Int64, so the row reduction overflows and the
    # variant gets no certificate: identify fits it.
    function _testhelper_exponent_toy_law(p, q)
        vars = EnzymeRates._Vars()
        a = EnzymeRates._var!(vars, :A; conc = true)
        k1, k2 = EnzymeRates._var!(vars, :k_one), EnzymeRates._var!(vars, :k_two)
        EnzymeRates._var!(vars, :Keq)
        EnzymeRates._var!(vars, :L)
        t(ps...) = EnzymeRates._term(EnzymeRates._monomial(ps...))
        D = one(EnzymeRates._Poly) + t(k1 => p, a => 1) + t(k1 => 1, k2 => q, a => 1)
        LA = EnzymeRates._Law(vars, t(a => 1), t(k1 => p), D, [D],
                              [one(EnzymeRates._Poly)], EnzymeRates._Mono())
        EnzymeRates._AlloLaw(vars, t(a => 1), LA, EnzymeRates._flip_only_law(LA),
                             one(EnzymeRates._Poly), one(EnzymeRates._Poly),
                             Tuple{EnzymeRates._Poly, EnzymeRates._Poly, Int}[], 1)
    end
    @test EnzymeRates._phantom_weights(_testhelper_exponent_toy_law(2, 3)) ==
          [-1 // 2, -1 // 6]
    p, q = 2^40 + 15, 2^40 + 3
    @test denominator(big(1 - p) // (big(p) * q)) > typemax(Int)
    @test isempty(EnzymeRates._phantom_weights(_testhelper_exponent_toy_law(p, q)))
end
