# ABOUTME: Tests of the Sig, the compiled rate law's only type parameter: it holds the whole
# ABOUTME: mechanism, decodes back to it, and its decoder compiles once for every mechanism.

const ER = EnzymeRates

"""Whether every leaf of the nested tuple `x`, a Sig or a part of one, is a Symbol, an Int
or a Bool."""
_testhelper_sig_leaves(x) =
    x isa Tuple ? all(_testhelper_sig_leaves, x) : x isa Union{Symbol, Int, Bool}

@testset "a mechanism's Sig decodes to the mechanism" begin
    @testset "$(spec.name)" for spec in MECHANISM_TEST_SPECS
        x = _testhelper_build(spec.raw)
        sig = ER.Derivation._sig_of(x)
        @test sig isa Tuple && _testhelper_sig_leaves(sig)
        @test _testhelper_decode(sig) == x
        rate = compile_rate_equation(x)
        @test typeof(rate).parameters[1] === sig
        @test _testhelper_decode(typeof(rate).parameters[1]) == x
    end
    m = @enzyme_mechanism begin
        substrates: S
        products: P
        steps: begin
            E + P ⇌ EP
            E + S ⇌ ES
            ES <--> EP
        end
    end
    @test _testhelper_decode(ER.Derivation._sig_of(m)) == m
    am = ER.Mechanism(
        m.catalytic, ER.Mechanisms.Allostery([:EqualAI, :NonequalAI, :OnlyA], 4,
        [ER.Mechanisms.RegulatorySite([:X, :R], 2, [:OnlyA, :OnlyI])]))
    @test _testhelper_decode(ER.Derivation._sig_of(am)) == am
    # A name outside ASCII, sorting last among 31 substrates.
    names = [[Symbol("S", i) for i in 1:30]; :αKG; :P]
    wide = _testhelper_plain_mechanism(
        names, [:E; [Symbol("E", x) for x in names]], collect(1:31), [32],
        [ER.Mechanisms.Step(1, i + 1, [i], Int[], true, i) for i in 1:32])
    @test wide.catalytic.metabolites[31] === :αKG
    @test _testhelper_decode(ER.Derivation._sig_of(wide)) == wide
end

"""Uni-uni with `n` modifiers, `I1` to `In`, each binding free E alone at rapid
equilibrium: a Sig whose type no other test decodes, its length set by `n`."""
function _testhelper_n_modifiers(n)
    mods = [Symbol("I", i) for i in 1:n]
    forms = [:E, :EP, :ES, [Symbol("E", x) for x in mods]...]
    steps = [ER.Mechanisms.Step(1, 3, [1], Int[], true, 1),
             ER.Mechanisms.Step(3, 2, Int[], Int[], false, 2),
             ER.Mechanisms.Step(1, 2, [2], Int[], true, 3),
             [ER.Mechanisms.Step(1, 3 + i, [2 + i], Int[], true, 3 + i) for i in 1:n]...]
    _testhelper_plain_mechanism([:S; :P; mods], forms, [1], [2], steps)
end

@testset "the Sig decoder compiles once, whatever the mechanism" begin
    # A compiled law's generated methods decode its Sig without specializing on its type,
    # so the decoder compiles once: after one decode, decoding the Sigs of three other
    # mechanisms, each of a Sig type of its own, adds no method instance.
    instances(f) = count(_ -> true, Base.specializations(only(methods(f))))
    _testhelper_allosteric_n_modifiers(n) =
        ER.Mechanism(_testhelper_n_modifiers(n).catalytic,
                     ER.Mechanisms.Allostery(fill(:EqualAI, n + 3), 2,
                                             ER.Mechanisms.RegulatorySite[]))
    decoder = ER.Derivation._mechanism_from_sig
    for build in (_testhelper_n_modifiers, _testhelper_allosteric_n_modifiers)
        counts = 41:44
        sigs = [typeof(compile_rate_equation(build(n))).parameters[1] for n in counts]
        @test allunique(typeof.(sigs))
        @test _testhelper_decode(sigs[1]) == build(first(counts))
        before = instances(decoder)
        foreach(_testhelper_decode, sigs[2:4])
        @test instances(decoder) == before
    end
end

@testset "_sig_of / _mechanism_from_sig roundtrip" begin
    # One fixture built from the walker's value: two substrates, two products, a
    # modifier-only name bound by a step, a covalent residual, a Theorell–Chance step,
    # and steps at rapid equilibrium and at steady state. It holds no reactant's
    # modifier copy.
    m = ER.DSL._construct(ER.DSL.raw_mechanism(:(begin
        substrates: A, B
        products: P, Q
        modifiers: I
        steps: begin
            E + A ⇌ EA
            EA <--> E_res_plusA_minusP + P
            E_res_plusA_minusP + B <--> EQ
            EA + B <--> EQ + P
            E + Q ⇌ EQ
            E + I ⇌ EI
        end
    end)))

    sig = ER.Derivation._sig_of(m)
    @test sig isa Tuple && _testhelper_sig_leaves(sig)
    # Its leaves are Symbols for the names, Ints for indices and the oligomeric state, and
    # Bools for the rapid-equilibrium flags.
    @test sig[1] === (:A, :B, :P, :Q, :I) && sig[7] === 1
    @test all(s -> s[5] isa Bool, sig[5])
    @test ER.Derivation._mechanism_from_sig(sig) == m   # roundtrip

    # The Sig is the compiled rate law's only parameter.
    rate = compile_rate_equation(m)
    @test typeof(rate) === ER.Derivation._CompiledRateEquation{sig}

    # Roundtrip through the type-parameter form preserves everything.
    @test ER.Derivation._mechanism_from_sig(typeof(rate).parameters[1]) == m

    # A regulator that no step binds is no part of a mechanism, so the init mechanism
    # of a reaction that declares one compiles to the law of the reaction's without
    # it.
    without_u = @enzyme_reaction begin
        substrates: S[C]
        products: P[C]
    end
    with_u = @enzyme_reaction begin
        substrates: S[C]
        products: P[C]
        competitive_inhibitors: U
    end
    @test compile_rate_equation(
        only(_testhelper_mechanisms(ER.Enumeration.init_mechanisms(with_u)))) ===
          compile_rate_equation(
              only(_testhelper_mechanisms(ER.Enumeration.init_mechanisms(without_u))))
end
