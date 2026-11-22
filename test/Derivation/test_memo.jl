# ABOUTME: Tests of src/Derivation/memo.jl: each answer about a mechanism equals the
# ABOUTME: answer on a fresh copy of it, in any order.

@testset "per-mechanism questions answer as on a fresh mechanism, in any order" begin
    # Identify asks one concrete mechanism several questions in a row. Each answer (a
    # value, or the error the call throws) equals the answer on a fresh copy of the
    # mechanism asked alone: in order, in reverse, and with two mechanisms alternated.
    # A mechanism the guard rejects (a bottomless rapid-equilibrium segment) is among
    # them.
    D = EnzymeRates.Derivation
    rejected = EnzymeRates.DSL._construct(EnzymeRates.DSL.raw_mechanism(:(begin
        substrates: S
        products: P, Q
        steps: begin
            E + S ⇌ ES
            ES <--> EPQ
            EPQ ⇌ EQ + P
            EPQ ⇌ EP + Q
            EP <--> E + P
            EQ <--> E + Q
        end
    end)))
    mechs = [[s.mechanism for s in MECHANISM_TEST_SPECS]; rejected]
    questions = (m -> D.check_valid_rate_equation(Bool, m),
                 D.check_valid_rate_equation,
                 m -> parameters(m), rate_equation_string,
                 m -> rate_equation_string(m; coefficients = true),
                 D.rate_equation_dedup_key)
    answer(q, m) = try q(m) catch e; e end
    expected = [[answer(q, deepcopy(m)) for q in questions] for m in mechs]
    @test expected[end][1] === false && expected[end][2] isa ErrorException
    for (m, want) in zip(mechs, expected)
        @test [answer(q, m) for q in questions] == want
        again = deepcopy(m)
        @test reverse([answer(q, again) for q in reverse(questions)]) == want
    end
    for i in 1:(length(mechs) - 1)
        a, b = deepcopy(mechs[i]), deepcopy(mechs[i + 1])
        got = [(answer(q, a), answer(q, b)) for q in questions]
        @test first.(got) == expected[i] && last.(got) == expected[i + 1]
    end
end
