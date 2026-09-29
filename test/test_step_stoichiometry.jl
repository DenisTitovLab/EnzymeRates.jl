# ABOUTME: Tests for Steps that carry explicit consumed/released metabolite lists:
# ABOUTME: construction, canonical orientation, kinds, names, orientation-free derivation.

const ER = EnzymeRates
_sp(bound, conf = :E) = ER.Species(ER.Metabolite[bound...], conf)
_sp(bound, conf, res) = ER.Species(ER.Metabolite[bound...], conf, res)

@testset "Step: explicit consumed/released lists" begin
    A, B, P, Q = ER.Substrate(:A), ER.Substrate(:B), ER.Product(:P), ER.Product(:Q)
    E, EA, EQ, EAB = _sp([]), _sp([A]), _sp([Q]), _sp([A, B])

    @testset "pure binding keeps its written orientation" begin
        s = ER.Step(E, EA, [A], ER.Metabolite[], true)
        @test ER.from_species(s) == E && ER.to_species(s) == EA
        @test ER.consumed(s) == ER.Metabolite[A] && isempty(ER.released(s))
        @test ER.ligand(s) == A && ER.is_binding(s) && !ER.is_iso(s)
    end

    @testset "a pure release is stored as the binding it reverses" begin
        s = ER.Step(EA, E, ER.Metabolite[], [A], false)
        @test s == ER.Step(E, EA, [A], ER.Metabolite[], false)
        # conformation change allowed: E*(A) → E + A is the binding E + A → E*(A)
        Estar_A = _sp([A], :Estar)
        r = ER.Step(Estar_A, E, ER.Metabolite[], [A], true)
        @test ER.from_species(r) == E && ER.to_species(r) == Estar_A
        @test ER.ligand(r) == A
    end

    @testset "isomerization and transformations" begin
        iso = ER.Step(EAB, _sp([P, Q]), ER.Metabolite[], ER.Metabolite[], false)
        @test ER.is_iso(iso) && ER.ligand(iso) === nothing && !ER.is_binding(iso)
        fused = ER.Step(EAB, EQ, ER.Metabolite[], [P], false)        # chemistry + release
        @test !ER.is_iso(fused) && ER.ligand(fused) === nothing
        @test ER.released(fused) == ER.Metabolite[P]
        tc = ER.Step(EA, EQ, [B], [P], false)                          # Theorell–Chance
        @test ER.consumed(tc) == ER.Metabolite[B] && ER.released(tc) == ER.Metabolite[P]
        @test ER.ligand(tc) === nothing
        two = ER.Step(E, EAB, [B, A], ER.Metabolite[], true)           # lists are sorted
        @test ER.consumed(two) == ER.Metabolite[A, B] && ER.ligand(two) === nothing
    end

    @testset "covalent residual: binding onto a residual form vs chemistry" begin
        res = ER.Residual([A], [P])
        F, FB = _sp([], :E, res), _sp([B], :E, res)
        @test ER.ligand(ER.Step(F, FB, [B], ER.Metabolite[], true)) == B
        chem = ER.Step(EA, F, ER.Metabolite[], [P], false)            # E(A) → F + P
        @test ER.ligand(chem) === nothing
        m = ER.Mechanism(
            @enzyme_reaction(begin
                substrates: A[CX], B[N]
                products: P[C], Q[NX]
            end),
            [[ER.Step(E, EA, [A], ER.Metabolite[], false)], [chem],
             [ER.Step(F, FB, [B], ER.Metabolite[], false)],
             [ER.Step(FB, E, ER.Metabolite[], [Q], false)]])
        @test_throws ErrorException ER._assert_chemistry_is_iso(m)
    end

    @testset "inhibitor copy stays distinct from the substrate" begin
        Ai = ER.CompetitiveInhibitor(:A)
        s_sub = ER.Step(E, EA, [A], ER.Metabolite[], true)
        s_inh = ER.Step(E, _sp([Ai]), [Ai], ER.Metabolite[], true)
        @test s_sub != s_inh && hash(s_sub) != hash(s_inh)
        @test ER.ligand(s_inh) == Ai
    end

    @testset "rejections" begin
        @test_throws ErrorException ER.Step(E, E, [A], ER.Metabolite[], true)
        @test_throws ErrorException ER.Step(EA, EQ, [A], [A], false)
    end

    @testset "sort key reproduces today's order" begin
        s = ER.Step(E, EA, [A], ER.Metabolite[], true)
        @test ER._step_canonical_key(s) == ("E", "EA", "A", "", true)
        iso = ER.Step(EAB, _sp([P, Q]), ER.Metabolite[], ER.Metabolite[], false)
        @test ER._step_canonical_key(iso) == ("EAB", "EPQ", "", "", false)
    end

    @testset "signature round trip" begin
        tc = ER.Step(EA, EQ, [B], [P], false)
        @test ER._step_from_sig(ER._to_sig(tc)) == tc
    end
end

@testset "transformation steps are named by form pair" begin
    m = @enzyme_mechanism begin
        substrates: A, B
        products: P, Q
        steps: begin
            E + A <--> E(A)
            E(A) + B <--> E(A, B)
            E(A, B) <--> E(Q) + P
            E(Q) <--> E + Q
        end
    end
    names = Set(ER.parameters(m, ER.Full))
    @test :k_EAB_to_EQ in names && :k_EQ_to_EAB in names
    @test !(:kon_P_EAB in names) && !(:koff_P_EAB in names)
    @test :kon_A_E in names && :koff_Q_E in names
end
