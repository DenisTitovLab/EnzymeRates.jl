@testset "DSL" begin
    @testset "DSL: decomposed-Species notation parses to Mechanism" begin
        # Decomposed call-form on every step → parser emits a structured
        # Mechanism (not a synthesized opaque Symbol). The to_species of
        # `E + S => E(S)` is a Species with conformation :E and bound
        # metabolite [Substrate(:S)].
        m = @enzyme_mechanism begin
            substrates: S
            products:   P
            steps: begin
                E + S ⇌ E(S)
                E(S) ⇌ E(P)
                E(P) ⇌ E + P
            end
        end
        mech = EnzymeRates.Mechanism(m)
        @test length(mech.steps) == 3
        es_step = first(mech.steps)[1]                  # E + S ⇌ E(S)
        @test EnzymeRates.conformation(es_step.to_species) == :E
        @test EnzymeRates.bound(es_step.to_species) ==
              EnzymeRates.Metabolite[EnzymeRates.Substrate(:S)]
        # The E + S side has conformation :E with no bound metabolite.
        @test EnzymeRates.conformation(es_step.from_species) == :E
        @test isempty(EnzymeRates.bound(es_step.from_species))

        # Second step: E(S) ⇌ E(P) — RE iso. The Mechanism constructor
        # canonicalizes iso direction physical-forward (substrate-bound
        # `from`, product-bound `to`) via `_canonical_iso_direction`, so
        # `E_S` is `from_species` and `E_P` is `to_species`.
        iso_step = mech.steps[2][1]
        @test EnzymeRates.bound(iso_step.from_species) ==
              EnzymeRates.Metabolite[EnzymeRates.Substrate(:S)]
        @test EnzymeRates.bound(iso_step.to_species) ==
              EnzymeRates.Metabolite[EnzymeRates.Product(:P)]

        # Multi-bound: E(A, B) → bound [Substrate(:A), Substrate(:B)]
        # (sorted by name). Mixed substrate/product/regulator roles
        # resolved from the declared metabolite blocks.
        m_multi = @enzyme_mechanism begin
            substrates: A, B
            products:   P
            regulators: I
            steps: begin
                E + A <--> E(A)
                E(A) + B <--> E(A, B)
                E(A, B) <--> E(P)
                E(P) <--> E + P
                E + I <--> E(I)
            end
        end
        mech_multi = EnzymeRates.Mechanism(m_multi)
        eab = mech_multi.steps[2][1].to_species         # E(A, B)
        @test EnzymeRates.conformation(eab) == :E
        @test EnzymeRates.bound(eab) ==
              EnzymeRates.Metabolite[EnzymeRates.Substrate(:A),
                                     EnzymeRates.Substrate(:B)]
        # Dead-end inhibitor lookup picks the correct Metabolite subtype.
        ei = mech_multi.steps[5][1].to_species          # E(I)
        @test EnzymeRates.bound(ei) ==
              EnzymeRates.Metabolite[
                  EnzymeRates.CompetitiveInhibitor(:I)]
    end

    @testset "@enzyme_mechanism: + step-side syntax" begin
        # New form: + separator, no brackets.
        m = @enzyme_mechanism begin
            substrates: S
            products:   P
            steps: begin
                E + S <--> E(S)
                E(S) <--> E + P
            end
        end
        @test m isa EnzymeMechanism
        @test EnzymeRates.n_steps(m) == 2
        @test Set(EnzymeRates.enzyme_forms(m)) == Set([:E, :ES])
    end

    @testset "@enzyme_mechanism multi-product balances placeholders" begin
        # Uni-bi (1 substrate, 2 products) — asymmetric reactant counts.
        # The placeholder atoms emitted by the macro must balance across
        # the substrate/product sides so the EnzymeReaction atom-balance
        # check passes.
        m = @enzyme_mechanism begin
            substrates: A
            products:   P, Q
            steps: begin
                E + A <--> E(A)
                E(A) <--> E(P, Q)
                E(P, Q) <--> E(Q) + P
                E(Q) <--> E + Q
            end
        end
        @test m isa EnzymeRates.EnzymeMechanism
    end

    @testset "@enzyme_mechanism decomposed-Species grammar" begin
        m = @enzyme_mechanism begin
            substrates: S
            products:   P
            regulators: I

            steps: begin
                (E + S ⇌ E(S), E(P) + S ⇌ E(P, S))
                E(S) + I ⇌ E(S, I)
                E(S) <--> E(P)
                E(P) ⇌ E + P
            end
        end
        @test EnzymeRates.substrates(m) == (:S,)
        @test EnzymeRates.products(m) == (:P,)
        @test EnzymeRates.regulators(m) == (:I,)
        @test EnzymeRates.kinetic_group(m, 1) == EnzymeRates.kinetic_group(m, 2)
        @test EnzymeRates.kinetic_group(m, 3) != EnzymeRates.kinetic_group(m, 4)

        # Function-call species notation: E(S) ≡ species with conformation :E
        # and bound metabolite :S. Synthesized form name is :ES
        # (matching `name(::Species)` from src/types.jl).
        m_call = @enzyme_mechanism begin
            substrates: S
            products:   P
            steps: begin
                E + S <--> E(S)
                E(S) <--> E(P)
                E(P) <--> E + P
            end
        end
        @test m_call isa EnzymeMechanism
        @test Set(EnzymeRates.enzyme_forms(m_call)) == Set([:E, :ES, :EP])
        @test EnzymeRates.n_steps(m_call) == 3

        # Multi-bound species: E(S, P) → :EPS (sorted alphabetically).
        m_multi = @enzyme_mechanism begin
            substrates: A, B
            products:   P, Q
            steps: begin
                E + A <--> E(A)
                E(A) + B <--> E(A, B)
                E(A, B) <--> E(P, Q)
                E(P, Q) <--> E(Q) + P
                E(Q) <--> E + Q
            end
        end
        @test m_multi isa EnzymeMechanism
        @test :EAB in EnzymeRates.enzyme_forms(m_multi)
        @test :EPQ in EnzymeRates.enzyme_forms(m_multi)

        # Residual notation: Estar(; residual = A - P).
        m_res = @enzyme_mechanism begin
            substrates: A, B
            products:   P, Q
            steps: begin
                E + A <--> E(A)
                E(A) <--> Estar(; residual = A - P)
                Estar(; residual = A - P) <--> Estar(Q; residual = A - P)
                Estar(Q; residual = A - P) <--> Estar(; residual = A - P) + Q
                Estar(; residual = A - P) + B <--> Estar(B; residual = A - P)
                Estar(B; residual = A - P) <--> E + P
            end
        end
        @test m_res isa EnzymeMechanism
        @test EnzymeRates.n_steps(m_res) == 6

        # Reject atom bracket syntax in substrates:
        @test_throws Exception eval(:(@enzyme_mechanism begin
            substrates: S[C]
            products:   P
            steps: begin
                E + S ⇌ E(S)
                E(S) <--> E(P)
                E(P) ⇌ E + P
            end
        end))

        # Reject allosteric-only syntax (regulatory_site(...))
        @test_throws Exception eval(:(@enzyme_mechanism begin
            substrates: S
            products:   P
            regulatory_site(multiplicity = 2): begin
                ligands: A
            end
            steps: begin
                E + S ⇌ E(S)
                E(S) <--> E(P)
                E(P) ⇌ E + P
            end
        end))
    end

    @testset "@allosteric_mechanism (parsing & validation)" begin
        m = @allosteric_mechanism begin
            substrates: F6P
            products:   F16BP
            catalytic_multiplicity: 2
            allosteric_regulators: I::OnlyI

            catalytic_steps: begin
                E + F6P ⇌ E(F6P)         :: EqualAI
                E(F6P) <--> E(F16BP)     :: EqualAI
                E(F16BP) ⇌ E + F16BP     :: EqualAI
            end
        end
        @test m isa EnzymeRates.AllostericEnzymeMechanism
        @test EnzymeRates.allosteric_regulators(m) ⊇ ((:I, :OnlyI),)

        # Reject untagged catalytic step
        @test_throws Exception eval(:(@allosteric_mechanism begin
            substrates: F6P
            products:   F16BP
            catalytic_multiplicity: 2
            catalytic_steps: begin
                E + F6P ⇌ E(F6P) :: EqualAI
                E(F6P) <--> E(F16BP)
                E(F16BP) ⇌ E + F16BP :: EqualAI
            end
        end))

        # Reject :OnlyI on a catalytic step (V-type allostery not supported)
        @test_throws Exception eval(:(@allosteric_mechanism begin
            substrates: F6P
            products:   F16BP
            catalytic_multiplicity: 2
            catalytic_steps: begin
                E + F6P ⇌ E(F6P) :: EqualAI
                E(F6P) <--> E(F16BP) :: OnlyI
                E(F16BP) ⇌ E + F16BP :: EqualAI
            end
        end))

        # Reject untagged allosteric regulator
        @test_throws Exception eval(:(@allosteric_mechanism begin
            substrates: F6P
            products:   F16BP
            catalytic_multiplicity: 2
            allosteric_regulators: I, J::OnlyI
            catalytic_steps: begin
                E + F6P ⇌ E(F6P) :: EqualAI
                E(F6P) <--> E(F16BP) :: EqualAI
                E(F16BP) ⇌ E + F16BP :: EqualAI
            end
        end))

        # Reject parenthesized step group without ::AlloState
        @test_throws Exception eval(:(@allosteric_mechanism begin
            substrates: S, A
            products:   P
            catalytic_multiplicity: 2
            catalytic_steps: begin
                (E + S ⇌ E(S), E(A) + S ⇌ E(A, S))
                E(A, S) <--> E(P)   :: EqualAI
                E(P) ⇌ E + P        :: EqualAI
            end
        end))
    end

    @testset "@enzyme_reaction" begin
        spec = @enzyme_reaction begin
            substrates: S[C]
            products:   P[C]
        end
        @test spec isa EnzymeReaction
        @test EnzymeRates.substrates(spec) == [EnzymeRates.Substrate(:S)]
        @test EnzymeRates.products(spec) == [EnzymeRates.Product(:P)]
        @test EnzymeRates.reactants(spec)[1] == EnzymeRates.ReactantAtoms(
            EnzymeRates.Product(:P), [:C => 1])
        @test EnzymeRates.reactants(spec)[2] == EnzymeRates.ReactantAtoms(
            EnzymeRates.Substrate(:S), [:C => 1])
        @test EnzymeRates.regulators(spec) == EnzymeRates.RegulatorMults[]

        spec2 = @enzyme_reaction begin
            substrates: S[C6H12O6], ATP[C10H16N5O13P3]
            products:   G6P[C6H13O9P], ADP[C10H15N5O10P2]
            competitive_inhibitors: I
        end
        @test length(EnzymeRates.substrates(spec2)) == 2
        @test length(EnzymeRates.products(spec2)) == 2
        @test length(EnzymeRates.regulators(spec2)) == 1
        @test EnzymeRates.regulator(EnzymeRates.regulators(spec2)[1]) ==
            EnzymeRates.CompetitiveInhibitor(:I)
    end

    @testset "multi-atom metabolites" begin
        rxn = @enzyme_reaction begin
            substrates: A[C2H3], B[N,P]
            products: P[C2,N], Q[H3,P]
        end
        subs = EnzymeRates.substrates(rxn)
        @test subs[1] == EnzymeRates.Substrate(:A)
        @test subs[2] == EnzymeRates.Substrate(:B)
        ra = EnzymeRates.reactants(rxn)
        ra_map = Dict(EnzymeRates.name(EnzymeRates.metabolite(r)) =>
                          EnzymeRates.atoms(r) for r in ra)
        @test ra_map[:A] == [:C => 2, :H => 3]
        @test ra_map[:B] == [:N => 1, :P => 1]
        prods = EnzymeRates.products(rxn)
        @test prods[1] == EnzymeRates.Product(:P)
        @test prods[2] == EnzymeRates.Product(:Q)
        @test ra_map[:P] == [:C => 2, :N => 1]
        @test ra_map[:Q] == [:H => 3, :P => 1]
    end

    @testset "@enzyme_reaction regulator kinds" begin
        # dead_end_inhibitors: and competitive_inhibitors: both emit
        # CompetitiveInhibitor entries. allosteric_regulators: emits
        # AllostericRegulator and requires per-name multiplicities.
        spec_kinds = @enzyme_reaction begin
            substrates: S[C]
            products: P[C]
            dead_end_inhibitors: I
            allosteric_regulators: A(1)
            competitive_inhibitors: R
        end
        @test spec_kinds isa EnzymeReaction
        regs = EnzymeRates.regulators(spec_kinds)
        reg_by_name = Dict(EnzymeRates.name(EnzymeRates.regulator(rm)) => rm
                           for rm in regs)
        @test Set(keys(reg_by_name)) == Set([:I, :A, :R])
        @test EnzymeRates.regulator(reg_by_name[:I]) ==
            EnzymeRates.CompetitiveInhibitor(:I)
        @test EnzymeRates.regulator(reg_by_name[:A]) ==
            EnzymeRates.AllostericRegulator(:A)
        @test EnzymeRates.regulator(reg_by_name[:R]) ==
            EnzymeRates.CompetitiveInhibitor(:R)
    end

    @testset "@enzyme_reaction rejects bare `regulators:` label" begin
        # The @enzyme_reaction grammar requires `competitive_inhibitors:`,
        # `dead_end_inhibitors:`, or `allosteric_regulators:`. A bare
        # `regulators:` label must be reported as unknown.
        @test_throws Exception eval(:(@enzyme_reaction begin
            substrates: S[C]
            products: P[C]
            regulators: R1, R2
        end))
    end

    @testset "opaque bound-form names are rejected" begin
        @test_throws "opaque bound-form name" eval(:(@enzyme_mechanism begin
            substrates: S
            products: P
            steps: begin
                E + S <--> ES
                ES <--> E + P
            end
        end))
    end

    @testset "@allosteric_mechanism rejects opaque bound-form names" begin
        @test_throws "opaque bound-form name" eval(:(@allosteric_mechanism begin
            substrates: S
            products: P
            allosteric_regulators: I::OnlyI
            catalytic_steps: begin
                E + S <--> ES :: EqualAI
                ES <--> E + P :: EqualAI
            end
        end))
    end

    @testset "@allosteric_mechanism opaque rejection names itself" begin
        err = try
            eval(:(@allosteric_mechanism begin
                substrates: S
                products: P
                allosteric_regulators: I::OnlyI
                catalytic_steps: begin
                    E + S <--> ES :: EqualAI
                    ES <--> E + P :: EqualAI
                end
            end))
            nothing
        catch e
            e
        end
        @test err !== nothing
        msg = err isa LoadError ? sprint(showerror, err.error) :
              sprint(showerror, err)
        @test occursin("@allosteric_mechanism", msg)
        @test !occursin("@enzyme_mechanism", msg)
    end

    @testset "@enzyme_reaction with oligomeric_state" begin
        rxn = @enzyme_reaction begin
            substrates: S[C]
            products: P[C]
            oligomeric_state: 4
        end
        @test EnzymeRates.allowed_catalytic_multiplicities(rxn) == [4]

        # Without oligomeric_state defaults to 1
        rxn2 = @enzyme_reaction begin
            substrates: S[C]
            products: P[C]
        end
        @test EnzymeRates.allowed_catalytic_multiplicities(rxn2) == [1]

        # Explicit allowed_catalytic_multiplicities tuple
        rxn3 = @enzyme_reaction begin
            substrates: S[C]
            products: P[C]
            allowed_catalytic_multiplicities: (1, 2, 4)
        end
        @test EnzymeRates.allowed_catalytic_multiplicities(rxn3) == [1, 2, 4]
    end

    @testset "@enzyme_mechanism" begin
        m = @enzyme_mechanism begin
            substrates: S
            products:   P
            steps: begin
                E + S <--> E(S)
                E(S) <--> E + P
            end
        end
        @test m isa EnzymeMechanism
        @test EnzymeRates.n_steps(m) == 2
        @test EnzymeRates.n_states(m) == 2
        @test Set(EnzymeRates.enzyme_forms(m)) == Set([:E, :ES])
        @test Set(metabolites(m)) == Set([:S, :P])

        # Numeric check: same as Uni-Uni spot check
        Keq = 3.2 * 2.5 / (0.8 * 1.1)
        params = (kon_S_E=3.2, kon_P_ES=2.5, koff_P_ES=1.1, Keq=Keq, E_total=1.0)
        concs = (S=0.7, P=0.3)
        @test rate_equation(m, concs, params) ≈ 0.9091 atol=0.001

        # Multi-step mechanism
        m2 = @enzyme_mechanism begin
            substrates: A, B
            products:   P, Q
            steps: begin
                E + A <--> E(A)
                E(A) <--> F(P)
                F(P) <--> F + P
                F + B <--> F(B)
                F(B) <--> E(Q)
                E(Q) <--> E + Q
            end
        end
        @test EnzymeRates.n_states(m2) == 6
    end

    @testset "Elementary steps" begin
        # The decomposed grammar's parser validates ≤1 metabolite/side
        # at macro-expansion time. Wrap in `eval(:(...))` so the macro
        # expansion happens at runtime where `@test_throws` can catch
        # the LoadError that wraps the parser's exception (same pattern
        # as lines 149, 160 above).
        @test_throws Exception eval(:(@enzyme_mechanism begin
            substrates: S
            products:   P
            steps: begin
                E + S + P <--> E(S, P)
            end
        end))

        spec = @enzyme_reaction begin
            substrates: S[C]
            products:   P[C]
            competitive_inhibitors: I
        end
        @test spec isa EnzymeReaction

        # Dead-end inhibitor: valid mechanism (competitive inhibition)
        m = @enzyme_mechanism begin
            substrates: S
            products:   P
            regulators: I
            steps: begin
                E + S <--> E(S)
                E(S) <--> E + P
                E + I <--> E(I)
            end
        end
        @test m isa EnzymeMechanism
    end

    @testset "No-atom species" begin
        # All metabolites without atoms
        m = @enzyme_mechanism begin
            substrates: S
            products:   P
            steps: begin
                E + S <--> E(S)
                E(S) <--> E + P
            end
        end
        @test m isa EnzymeMechanism
        @test EnzymeRates.n_steps(m) == 2
    end

    @testset "Constraint DSL parsing" begin
        # Bi-bi random with two K_A binding steps in shared kinetic group
        m = @enzyme_mechanism begin
            substrates: A, B
            products:   P, Q
            steps: begin
                (E + A ⇌ E(A), E(B) + A ⇌ E(A, B))
                E + B ⇌ E(B)
                E(A) + B ⇌ E(A, B)
                E(A, B) <--> E(P, Q)
                E(P, Q) ⇌ E(Q) + P
                E(Q) ⇌ E + Q
            end
        end
        @test m isa EnzymeMechanism
        @test EnzymeRates.kinetic_group(m, 1) == EnzymeRates.kinetic_group(m, 2)
        @test EnzymeRates.kinetic_group(m, 3) != EnzymeRates.kinetic_group(m, 1)
    end

    @testset "::Inh role tag: product that also competitively inhibits" begin
        m = @enzyme_mechanism begin
            substrates: S
            products:   P
            steps: begin
                E + S <--> E(S)
                E(S) <--> E(P)
                E(P) <--> E + P
                E + P::Inh <--> E(P::Inh)
            end
        end
        forms = EnzymeRates.enzyme_forms(m)
        @test :EP in forms        # product-bound form (catalytic release)
        @test :EPinh in forms     # inhibitor-bound form — DISTINCT from :EP
        @test :EP != :EPinh
        # rate_equation must reference concs.P for the inhibitor term, not concs.Pinh:
        re_str = rate_equation_string(m)
        @test occursin("P", re_str)             # real metabolite name appears
        @test occursin("(; S, P) = concs", re_str)  # concentration is P, never Pinh
        @test !occursin("Pinh = concs", re_str)     # the inh marker never names a concentration
    end

    @testset "fused catalytic release: metabolite in neither bound list dissociates" begin
        # E(A) <--> E(Q) + P : P is produced by the step (in neither E(A) nor
        # E(Q) bound list, and both sides are 1-bound). It must reconstruct as
        # a release (P on the rhs), not a binding.
        m = @enzyme_mechanism begin
            substrates: A
            products: P, Q
            steps: begin
                E + A <--> E(A)
                E(A) <--> E(Q) + P
                E(Q) <--> E + Q
            end
        end
        rxns = EnzymeRates.reactions(m)
        mid = rxns[2]            # (lhs, rhs, is_eq, g) for E(A) <--> E(Q) + P
        @test mid[1] == (:EA,)  # lhs: only the enzyme form
        @test :P in mid[2]       # P released (rhs)
        @test :P ∉ mid[1]        # P not consumed (lhs)
    end
end
