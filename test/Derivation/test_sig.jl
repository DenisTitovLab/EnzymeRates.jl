# ABOUTME: Tests of the Sig, the compiled rate law's only type parameter: it holds the whole
# ABOUTME: mechanism, decodes back to it, and its decoder compiles once for every mechanism.

const ER = EnzymeRates

@testset "a mechanism's Sig decodes to the mechanism" begin
    @testset "$(spec.name)" for spec in MECHANISM_TEST_SPECS
        x = _testhelper_build(spec.raw)
        sig = ER.Derivation._sig_of(x)
        @test sig isa Tuple
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

"""Uni-uni whose substrate's name holds `n` characters: a Sig no other test decodes."""
_testhelper_long_named(n) = let S = Symbol("S", "x"^(n - 1)), ES = Symbol("ES", "x"^(n - 1))
    _testhelper_build(EnzymeRates.DSL.raw_mechanism(:(begin
        substrates: $S
        products: P
        steps: begin
            E + $S ⇌ $ES
            $ES <--> EP
            EP ⇌ E + P
        end
    end)))
end

@testset "the Sig decoder compiles once, whatever the mechanism" begin
    # A compiled law's generated methods decode its Sig without specializing on its type,
    # so the decoder compiles once: after one decode, decoding the Sigs of three other
    # mechanisms adds no method instance.
    instances(f) = count(_ -> true, Base.specializations(only(methods(f))))
    _testhelper_allosteric_long_named(n) =
        ER.Mechanism(_testhelper_long_named(n).catalytic,
                     ER.Mechanisms.Allostery([:EqualAI, :EqualAI, :OnlyA], 2,
                                             EnzymeRates.Mechanisms.RegulatorySite[]))
    decoder = ER.Derivation._mechanism_from_sig
    for build in (_testhelper_long_named, _testhelper_allosteric_long_named)
        lengths = (301, 401, 501, 601)
        sigs = [typeof(compile_rate_equation(build(n))).parameters[1] for n in lengths]
        @test _testhelper_decode(sigs[1]) == build(lengths[1])
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
    @test sig isa Tuple
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
