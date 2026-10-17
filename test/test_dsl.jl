# ABOUTME: Tests for the @enzyme_reaction / @enzyme_mechanism DSL parsing,
# ABOUTME: including decomposed-Species notation and opaque-form rejection.
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
        # Steps are canonicalized; pick each by content, not position.
        es_step = only(s for g in mech.steps for s in g           # E + S ⇌ E(S)
            if EnzymeRates.bound_metabolite(s) == EnzymeRates.Substrate(:S))
        @test EnzymeRates.conformation(es_step.to_species) == :E
        @test EnzymeRates.bound(es_step.to_species) ==
              EnzymeRates.Metabolite[EnzymeRates.Substrate(:S)]
        # The E + S side has conformation :E with no bound metabolite.
        @test EnzymeRates.conformation(es_step.from_species) == :E
        @test isempty(EnzymeRates.bound(es_step.from_species))

        # Second step: E(S) ⇌ E(P) — RE iso. The Mechanism constructor
        # canonicalizes iso direction physical-forward (substrate-bound
        # `from`, product-bound `to`) via `_canonical_step_direction`, so
        # `E_S` is `from_species` and `E_P` is `to_species`.
        iso_step = only(s for g in mech.steps for s in g
                        if EnzymeRates.bound_metabolite(s) === nothing)
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
            modifiers: I
            steps: begin
                E + A <--> E(A)
                E(A) + B <--> E(A, B)
                E(A, B) <--> E(P)
                E(P) <--> E + P
                E + I <--> E(I)
            end
        end
        mech_multi = EnzymeRates.Mechanism(m_multi)
        eab = only(EnzymeRates.to_species(s)            # E(A, B)
            for g in mech_multi.steps for s in g
            if Set(EnzymeRates.name(b)
                   for b in EnzymeRates.bound(EnzymeRates.to_species(s))) ==
               Set([:A, :B]))
        @test EnzymeRates.conformation(eab) == :E
        @test EnzymeRates.bound(eab) ==
              EnzymeRates.Metabolite[EnzymeRates.Substrate(:A),
                                     EnzymeRates.Substrate(:B)]
        # Dead-end inhibitor lookup picks the correct Metabolite subtype.
        ei = only(EnzymeRates.to_species(s)             # E(I)
            for g in mech_multi.steps for s in g
            if EnzymeRates.bound_metabolite(s) ==
               EnzymeRates.CompetitiveInhibitor(:I))
        @test EnzymeRates.bound(ei) ==
              EnzymeRates.Metabolite[
                  EnzymeRates.CompetitiveInhibitor(:I)]
    end

    @testset "@enzyme_mechanism decomposed-Species grammar" begin
        m = @enzyme_mechanism begin
            substrates: S
            products:   P
            modifiers: I

            steps: begin
                (E + S ⇌ E(S), E(P) + S ⇌ E(P, S))
                E(S) + I ⇌ E(S, I)
                E(S) <--> E(P)
                E(P) ⇌ E + P
            end
        end
        @test _testhelper_substrates(m) == [:S]
        @test _testhelper_products(m) == [:P]
        @test _testhelper_regulators(m) == [:I]
        # One shared 2-step kinetic group (the S-binding pair) + 3 singletons →
        # 5 steps, 4 groups (order-independent: canonicalization reorders steps).
        @test length(EnzymeRates.steps(EnzymeRates.Mechanism(m))) == 4

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
        @test sum(length, EnzymeRates.steps(EnzymeRates.Mechanism(m_res))) == 6

        # Reject atom bracket syntax in substrates:
        @test_throws "expects bare names" eval(:(@enzyme_mechanism begin
            substrates: S[C]
            products:   P
            steps: begin
                E + S ⇌ E(S)
                E(S) <--> E(P)
                E(P) ⇌ E + P
            end
        end))

        # Reject allosteric-only syntax (regulatory_site(...))
        @test_throws "belong in @allosteric_mechanism" eval(:(@enzyme_mechanism begin
            substrates: S
            products:   P
            regulatory_site(2): A
            steps: begin
                E + S ⇌ E(S)
                E(S) <--> E(P)
                E(P) ⇌ E + P
            end
        end))
    end

    @testset "@allosteric_mechanism (parsing & validation)" begin
        # Reject untagged catalytic step
        @test_throws Exception eval(:(@allosteric_mechanism begin
            substrates: F6P
            products:   F16BP
            catalytic_multiplicity: 2
            steps: begin
                E + F6P ⇌ E(F6P) :: EqualAI
                E(F6P) <--> E(F16BP)
                E(F16BP) ⇌ E + F16BP :: EqualAI
            end
        end))

        # Reject :OnlyI on a catalytic step (V-type allostery not supported)
        @test_throws ":OnlyI is rejected for catalytic groups" eval(
            :(@allosteric_mechanism begin
                substrates: F6P
                products:   F16BP
                catalytic_multiplicity: 2
                steps: begin
                    E + F6P ⇌ E(F6P) :: EqualAI
                    E(F6P) <--> E(F16BP) :: OnlyI
                    E(F16BP) ⇌ E + F16BP :: EqualAI
                end
            end))

        # Reject an allosteric regulator without a tag: it sits at no regulatory site
        @test_throws("@allosteric_mechanism: allosteric regulator `I` sits at no " *
                     "regulatory site; list it on one `regulatory_site(n):` line, as in " *
                     "`regulatory_site(2): I::OnlyA`.", eval(:(@allosteric_mechanism begin
            substrates: F6P
            products:   F16BP
            catalytic_multiplicity: 2
            allosteric_regulators: I, J
            steps: begin
                E + F6P ⇌ E(F6P) :: EqualAI
                E(F6P) <--> E(F16BP) :: EqualAI
                E(F16BP) ⇌ E + F16BP :: EqualAI
            end
            regulatory_site(2): J::OnlyI
        end)))

        # Reject a parenthesized step group without a tag.
        @test_throws("@allosteric_mechanism: parenthesized step group " *
                     "`(E + S ⇌ E(S), E(A) + S ⇌ E(A, S))` is missing its tag, " *
                     "`:: EqualAI`, `:: NonequalAI` or `:: OnlyA`; add it after the " *
                     "closing paren.", eval(:(@allosteric_mechanism begin
            substrates: S, A
            products:   P
            catalytic_multiplicity: 2
            steps: begin
                (E + S ⇌ E(S), E(A) + S ⇌ E(A, S))
                E(A, S) <--> E(P)   :: EqualAI
                E(P) ⇌ E + P        :: EqualAI
            end
        end)))
    end

    @testset "@allosteric_mechanism checks names and sites; constructors check states" begin
        # A name listed twice in `allosteric_regulators:` is rejected, with one site line
        # for it or two.
        @test_throws "lists `A` more than once" eval(:(@allosteric_mechanism begin
            substrates: S
            products:   P
            allosteric_regulators: A, A
            steps: begin
                E + S ⇌ E(S)      :: EqualAI
                E(S) <--> E(P)    :: EqualAI
                E(P) ⇌ E + P      :: EqualAI
            end
            regulatory_site(2): A::OnlyA
        end))
        @test_throws "lists `A` more than once" eval(:(@allosteric_mechanism begin
            substrates: S
            products:   P
            allosteric_regulators: A, A
            steps: begin
                E + S ⇌ E(S)      :: EqualAI
                E(S) <--> E(P)    :: EqualAI
                E(P) ⇌ E + P      :: EqualAI
            end
            regulatory_site(1): A::OnlyA
            regulatory_site(1): A::OnlyA
        end))

        # An unknown regulator state is reported by RegulatorySite when the expansion
        # runs.
        @test_throws "allo state Bogus must be one of" eval(:(@allosteric_mechanism begin
            substrates: S
            products:   P
            allosteric_regulators: X
            steps: begin
                E + S ⇌ E(S)      :: EqualAI
                E(S) <--> E(P)    :: EqualAI
                E(P) ⇌ E + P      :: EqualAI
            end
            regulatory_site(1): X::Bogus
        end))

        # Unknown and :OnlyI catalytic states are reported by AllostericMechanism.
        for state in (:Bogus, :OnlyI)
            @test_throws "catalytic group" eval(:(@allosteric_mechanism begin
                substrates: S
                products:   P
                steps: begin
                    E + S ⇌ E(S)      :: EqualAI
                    E(S) <--> E(P)    :: EqualAI
                    E(P) ⇌ E + P     :: $state
                end
            end))
        end

        # A ligand on two regulatory sites is reported by the macro, which reads the
        # site lines.
        @test_throws("@allosteric_mechanism: allosteric regulator `A` is listed 2 times " *
                     "on `regulatory_site(n):` lines; an allosteric regulator sits at " *
                     "exactly one regulatory site.", eval(
            :(@allosteric_mechanism begin
                substrates: S
                products:   P
                allosteric_regulators: A
                steps: begin
                    E + S ⇌ E(S)      :: EqualAI
                    E(S) <--> E(P)    :: EqualAI
                    E(P) ⇌ E + P      :: EqualAI
                end
                regulatory_site(2): A::OnlyA
                regulatory_site(4): A::OnlyA
            end)))
    end

    @testset "a mechanism macro's single-valued labels are written once" begin
        # A repeated label errors and names the label; it never keeps the last line.
        @test_throws "`steps:` given more than once" eval(:(@enzyme_mechanism begin
            substrates: S
            products:   P
            steps: begin
                E + S ⇌ E(S)
                E(S) <--> E(P)
                E(P) ⇌ E + P
            end
            steps: begin
                E + S ⇌ E(S)
                E(S) <--> E(P)
                E(P) ⇌ E + P
            end
        end))
        @test_throws "`catalytic_multiplicity:` given more than once" eval(
            :(@allosteric_mechanism begin
                substrates: S
                products:   P
                catalytic_multiplicity: 2
                catalytic_multiplicity: 4
                steps: begin
                    E + S ⇌ E(S)      :: EqualAI
                    E(S) <--> E(P)    :: EqualAI
                    E(P) ⇌ E + P      :: EqualAI
                end
            end))
    end

    @testset "a mechanism macro's list labels take bare names" begin
        # An atom bracket and a role tag each get a hint where they go; another value
        # gets none.
        @test_throws("`substrates:` expects bare names; got S[C]. Atom brackets belong " *
                     "in @enzyme_reaction.", eval(:(@enzyme_mechanism begin
            substrates: S[C]
            products:   P
            steps: begin
                E + S ⇌ E(S)
                E(S) <--> E(P)
                E(P) ⇌ E + P
            end
        end)))
        @test_throws("`modifiers:` expects bare names; got P::Mod. List `P` bare; a role " *
                     "tag goes on the metabolite in a step.", eval(:(@enzyme_mechanism begin
            substrates: S
            products:   P
            modifiers: P::Mod
            steps: begin
                E + S ⇌ E(S)
                E(S) <--> E(P)
                E(P) ⇌ E + P
                E + P::Mod ⇌ E(P::Mod)
            end
        end)))
        @test_throws("`allosteric_regulators:` expects bare names; got X(2).",
                     eval(:(@allosteric_mechanism begin
            substrates: S
            products:   P
            allosteric_regulators: X(2)
            steps: begin
                E + S ⇌ E(S)      :: EqualAI
                E(S) <--> E(P)    :: EqualAI
                E(P) ⇌ E + P      :: EqualAI
            end
            regulatory_site(2): X::OnlyA
        end)))
    end

    @testset "@enzyme_mechanism rejects each allosteric-only label" begin
        for line in (:(allosteric_regulators: A),
                     :(catalytic_multiplicity: 2),
                     :(regulatory_site(2): A::OnlyA))
            @test_throws "belong in @allosteric_mechanism" eval(:(@enzyme_mechanism begin
                substrates: S
                products:   P
                $line
                steps: begin
                    E + S ⇌ E(S)
                    E(S) <--> E(P)
                    E(P) ⇌ E + P
                end
            end))
        end
    end

    @testset "mechanism macros name themselves in body errors" begin
        # A macro argument that is not a `begin ... end` block.
        @test_throws "@enzyme_mechanism: expected a `begin ... end` block" eval(
            :(@enzyme_mechanism S))
        @test_throws "@allosteric_mechanism: expected a `begin ... end` block" eval(
            :(@allosteric_mechanism S))

        # A missing required label, in either macro.
        @test_throws "@enzyme_mechanism: `steps:` not specified" eval(
            :(@enzyme_mechanism begin
                substrates: S
                products:   P
            end))
        @test_throws "@enzyme_mechanism: `substrates:` not specified" eval(
            :(@enzyme_mechanism begin
                products: P
            end))
        @test_throws "@allosteric_mechanism: `steps:` not specified" eval(
            :(@allosteric_mechanism begin
                substrates: S
                products:   P
            end))

        # A repeated `steps:` line names the label.
        @test_throws "@allosteric_mechanism: `steps:` given more than once" eval(
            :(@allosteric_mechanism begin
                substrates: S
                products:   P
                steps: begin
                    E + S ⇌ E(S)      :: EqualAI
                    E(S) <--> E(P)    :: EqualAI
                    E(P) ⇌ E + P      :: EqualAI
                end
                steps: begin
                    E + S ⇌ E(S)      :: EqualAI
                    E(S) <--> E(P)    :: EqualAI
                    E(P) ⇌ E + P      :: EqualAI
                end
            end))
    end

    @testset "regulatory_site takes exactly one multiplicity" begin
        # The multiplicity is one positional Int: keyword, repeated, missing, misnamed and
        # second arguments are all errors.
        for site in (:(regulatory_site(multiplicity = 2)),
                     :(regulatory_site(multiplicity = 2, multiplicity = 4)),
                     :(regulatory_site()),
                     :(regulatory_site(size = 2)),
                     :(regulatory_site(2, 4)))
            @test_throws("@allosteric_mechanism: `regulatory_site(n):` takes one " *
                         "positive Int multiplicity, as in `regulatory_site(2): " *
                         "A::OnlyA, I::OnlyI`; got `$site`.",
                         eval(:(@allosteric_mechanism begin
                substrates: S
                products:   P
                allosteric_regulators: A
                steps: begin
                    E + S ⇌ E(S)      :: EqualAI
                    E(S) <--> E(P)    :: EqualAI
                    E(P) ⇌ E + P      :: EqualAI
                end
                $site: A::OnlyA
            end)))
        end
    end

    @testset "a tagged tuple of one step is a one-step group" begin
        # A parenthesized group is a tuple of steps, `(step,)` for one: tagged, it is the
        # same mechanism as the bare tagged step.
        one_tuple = EnzymeRates.AllostericMechanism(@allosteric_mechanism begin
            substrates: F6P
            products:   F16BP
            catalytic_multiplicity: 2
            allosteric_regulators: A, I

            steps: begin
                E + F6P ⇌ E(F6P)        :: EqualAI
                E(F6P) <--> E(F16BP)    :: EqualAI
                (E(F16BP) ⇌ E + F16BP,) :: EqualAI
            end

            regulatory_site(4): A::OnlyA
            regulatory_site(4): I::OnlyI
        end)
        bare = EnzymeRates.AllostericMechanism(@allosteric_mechanism begin
            substrates: F6P
            products:   F16BP
            catalytic_multiplicity: 2
            allosteric_regulators: A, I

            steps: begin
                E + F6P ⇌ E(F6P)        :: EqualAI
                E(F6P) <--> E(F16BP)    :: EqualAI
                E(F16BP) ⇌ E + F16BP    :: EqualAI
            end

            regulatory_site(4): A::OnlyA
            regulatory_site(4): I::OnlyI
        end)
        @test one_tuple == bare
        @test length(EnzymeRates.steps(one_tuple)) == 3

        # A plain mechanism has no group tags, parenthesized or not.
        @test_throws("@enzyme_mechanism: tag annotation `(E(P) ⇌ E + P)::EqualAI` is not " *
                     "allowed", eval(:(@enzyme_mechanism begin
            substrates: S
            products:   P
            steps: begin
                E + S ⇌ E(S)
                E(S) <--> E(P)
                (E(P) ⇌ E + P) :: EqualAI
            end
        end)))
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
    end

    @testset "@enzyme_reaction rejects a name listed as substrate and product" begin
        # Concentrations and constants are keyed by name, so a name may hold only
        # one reactant role; either label order is rejected.
        both_roles = ErrorException(
            "EnzymeReaction: S named as both a substrate and a product; " *
            "concentrations and constants are keyed by name")
        @test_throws both_roles @enzyme_reaction begin
            products:   S[C]
            substrates: S[C]
        end
        @test_throws both_roles @enzyme_reaction begin
            substrates: S[C]
            products:   S[C]
        end
    end

    @testset "multi-atom metabolites" begin
        rxn = @enzyme_reaction begin
            substrates: A[C2H3], B[NP]
            products: P[C2N], Q[H3P]
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

        # A reactant takes one formula in brackets: several formulas, an empty bracket and
        # a bare name are errors.
        for entry in (:(B[N, P]), :(B[]), :B)
            @test_throws("@enzyme_reaction: `substrates:` lists each reactant with one " *
                         "chemical formula in brackets, such as `S[C6H12O6]`; got " *
                         "`$entry`.", eval(:(@enzyme_reaction begin
                substrates: A[C2H3], $entry
                products: P[C2N], Q[H3P]
            end)))
        end
    end

    @testset "chemical formulas and multiplicities are validated" begin
        # A formula is a run of elements, each an uppercase letter with optional
        # lowercase letters and an optional count.
        for formula in (:c6, :C6x, 6)
            @test_throws "Invalid chemical formula: \"$formula\"" eval(
                :(@enzyme_reaction begin
                    substrates: S[$formula]
                    products:   P[C6]
                end))
        end

        # Every multiplicity is a positive Int.
        @test_throws "regulator A multiplicity must be a positive Int, got 0." eval(
            :(@enzyme_reaction begin
                substrates: S[C]
                products:   P[C]
                allosteric_regulators: A(0)
            end))
        @test_throws "`oligomeric_state:` entry must be a positive Int, got 0." eval(
            :(@enzyme_reaction begin
                substrates: S[C]
                products:   P[C]
                oligomeric_state: 0
            end))
        @test_throws "`catalytic_multiplicity:` must be a positive Int, got 2.0." eval(
            :(@allosteric_mechanism begin
                substrates: S
                products:   P
                catalytic_multiplicity: 2.0
                steps: begin
                    E + S ⇌ E(S)      :: EqualAI
                    E(S) <--> E(P)    :: EqualAI
                    E(P) ⇌ E + P      :: EqualAI
                end
            end))
        @test_throws "`regulatory_site` multiplicity must be a positive Int, got 0." eval(
            :(@allosteric_mechanism begin
                substrates: S
                products:   P
                allosteric_regulators: A
                steps: begin
                    E + S ⇌ E(S)      :: EqualAI
                    E(S) <--> E(P)    :: EqualAI
                    E(P) ⇌ E + P      :: EqualAI
                end
                regulatory_site(0): A::OnlyA
            end))
    end

    @testset "@enzyme_reaction regulator kinds" begin
        # competitive_inhibitors: emits CompetitiveInhibitor entries, and its lines
        # accumulate. allosteric_regulators: emits AllostericRegulator; a bare entry
        # takes the catalytic multiplicities.
        spec_kinds = @enzyme_reaction begin
            substrates: S[C]
            products: P[C]
            competitive_inhibitors: I
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

    @testset "@enzyme_reaction dual-role regulator (allosteric + competitive)" begin
        # One metabolite MAY be declared as BOTH an allosteric regulator and a
        # competitive inhibitor; the two roles render to distinct parameter
        # names, so both RegulatorMults are carried.
        rxn = @enzyme_reaction begin
            substrates: S[C]
            products: P[C]
            competitive_inhibitors: ATP
            allosteric_regulators: ATP(1)
        end
        @test rxn isa EnzymeReaction
        atp_regs = [rm for rm in EnzymeRates.regulators(rxn)
                    if EnzymeRates.name(EnzymeRates.regulator(rm)) == :ATP]
        @test length(atp_regs) == 2
        @test Set(typeof(EnzymeRates.regulator(rm)) for rm in atp_regs) ==
              Set([EnzymeRates.AllostericRegulator,
                   EnzymeRates.CompetitiveInhibitor])
        # The order of the two labels does not matter: both orders build one reaction.
        allosteric_first = @enzyme_reaction begin
            substrates: S[C]
            products: P[C]
            allosteric_regulators: ATP(1)
            competitive_inhibitors: ATP
        end
        @test allosteric_first == rxn
        @test hash(allosteric_first) == hash(rxn)
        # Two allosteric ATPs are rejected.
        @test_throws Exception eval(:(@enzyme_reaction begin
            substrates: S[C]
            products: P[C]
            allosteric_regulators: ATP(1), ATP(2)
        end))
    end

    @testset "@enzyme_reaction regulator types" begin
        rxn = @enzyme_reaction begin
            substrates: A[C6H12O6]
            products: B[C6H12O6]
            allosteric_regulators: X::Activator, Y::Inhibitor, Z
            oligomeric_state: 2
        end
        reg_types = Dict(EnzymeRates.name(EnzymeRates.regulator(rm)) =>
                             EnzymeRates.reg_type(rm)
                         for rm in EnzymeRates.regulators(rxn))
        @test reg_types[:X] == :activator
        @test reg_types[:Y] == :inhibitor
        @test reg_types[:Z] == :unspecified

        # Call-form entry with explicit multiplicities plus a type tag:
        # `X(1,2)::Inhibitor` — the type peels, then the inner `X(1,2)`
        # parses through the call-form branch for its multiplicities.
        rxn_call = @enzyme_reaction begin
            substrates: A[C6H12O6]
            products: B[C6H12O6]
            allosteric_regulators: X(1, 2)::Inhibitor
            oligomeric_state: 2
        end
        rm_call = only(EnzymeRates.regulators(rxn_call))
        @test EnzymeRates.allowed_multiplicities(rm_call) == [1, 2]
        @test EnzymeRates.reg_type(rm_call) == :inhibitor

        @test_throws Exception eval(:(@enzyme_reaction begin
            substrates: A[C6H12O6]
            products: B[C6H12O6]
            allosteric_regulators: X::Bogus
            oligomeric_state: 2
        end))
        @test_throws Exception eval(:(@enzyme_reaction begin
            substrates: A[C6H12O6]
            products: B[C6H12O6]
            competitive_inhibitors: X::Activator
        end))
        @test_throws("@enzyme_reaction: regulator type ::Inhibitor on X::Inhibitor is " *
                     "only valid on allosteric_regulators: entries.",
                     eval(:(@enzyme_reaction begin
            substrates: A[C6H12O6]
            products: B[C6H12O6]
            competitive_inhibitors: X::Inhibitor
        end)))
    end

    @testset "@enzyme_reaction names its labels for any other" begin
        # Regulators go under `competitive_inhibitors:` or `allosteric_regulators:`, the
        # catalytic multiplicities under `oligomeric_state:`; any other label is reported
        # as unknown, with the valid ones.
        for label in (:regulators, :dead_end_inhibitors, :allowed_catalytic_multiplicities)
            @test_throws("@enzyme_reaction: unknown label `$label:`. Valid labels: " *
                         "substrates, products, competitive_inhibitors, " *
                         "allosteric_regulators, oligomeric_state, shared_catalytic_site.",
                         eval(:(@enzyme_reaction begin
                substrates: S[C]
                products: P[C]
                $label: R1, R2
            end)))
        end
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

    @testset "every conformation label is conformation-shaped" begin
        # A call head is a conformation label like a bare name: `ES(P)` is an opaque
        # bound-form name written as a call.
        @test_throws "`ES` looks like an opaque bound-form name" eval(
            :(@enzyme_mechanism begin
                substrates: S
                products: P
                steps: begin
                    E + S <--> E(S)
                    E(S) <--> ES(P)
                    ES(P) <--> E + P
                end
            end))
        # A bare `ER` is no more acceptable for also heading a call elsewhere.
        @test_throws "`ER` looks like an opaque bound-form name" eval(
            :(@enzyme_mechanism begin
                substrates: S
                products: P
                steps: begin
                    E + S <--> ER(S)
                    ER(S) <--> ER
                    ER <--> E + P
                end
            end))
    end

    @testset "an undeclared metabolite is reported ahead of the opaque-name check" begin
        # `ATP` is undeclared: it is a misspelled metabolite, not a bound-form name.
        @test_throws "bound metabolite `ATP` in species `E(ATP)` is not declared" eval(
            :(@enzyme_mechanism begin
                substrates: S
                products: P
                steps: begin
                    E + ATP ⇌ E(ATP)
                    E(S) <--> E(P)
                    E(P) ⇌ E + P
                end
            end))
        @test_throws "bound metabolite `ATP` in species `E(ATP)` is not declared" eval(
            :(@allosteric_mechanism begin
                substrates: S
                products: P
                steps: begin
                    E + ATP ⇌ E(ATP)    :: EqualAI
                    E + S ⇌ E(S)        :: EqualAI
                    E(S) <--> E(P)      :: EqualAI
                    E(P) ⇌ E + P        :: EqualAI
                end
            end))
    end

    @testset "@allosteric_mechanism opaque rejection names itself" begin
        err = try
            eval(:(@allosteric_mechanism begin
                substrates: S
                products: P
                allosteric_regulators: I
                steps: begin
                    E + S <--> ES :: EqualAI
                    ES <--> E + P :: EqualAI
                end
                regulatory_site(1): I::OnlyI
            end))
            nothing
        catch e
            e
        end
        @test err !== nothing
        msg = err isa LoadError ? sprint(showerror, err.error) :
              sprint(showerror, err)
        @test occursin("opaque bound-form name", msg)
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

        # A tuple of multiplicities
        rxn3 = @enzyme_reaction begin
            substrates: S[C]
            products: P[C]
            oligomeric_state: (1, 2, 4)
        end
        @test EnzymeRates.allowed_catalytic_multiplicities(rxn3) == [1, 2, 4]
    end

    @testset "@enzyme_reaction sets the catalytic multiplicities once" begin
        # A second `oligomeric_state:` line errors; it never replaces the first.
        @test_throws "@enzyme_reaction: `oligomeric_state:` given more than once." eval(
            :(@enzyme_reaction begin
                substrates: S[C]
                products:   P[C]
                oligomeric_state: 2
                oligomeric_state: (1, 2)
            end))

        # `oligomeric_state:` takes one Int or one tuple of them, so a search can span
        # several multiplicities.
        one_tuple = eval(:(@enzyme_reaction begin
            substrates: S[C]
            products:   P[C]
            oligomeric_state: (2,)
        end))
        @test EnzymeRates.allowed_catalytic_multiplicities(one_tuple) == [2]
        several = @enzyme_reaction begin
            substrates: S[C]
            products:   P[C]
            oligomeric_state: (4, 1, 2)
        end
        @test EnzymeRates.allowed_catalytic_multiplicities(several) == [1, 2, 4]
        # An empty tuple leaves the search no multiplicity.
        @test_throws("@enzyme_reaction: `oligomeric_state:` takes one Int or a tuple of " *
                     "one or more Ints, such as `(1, 2, 4)`; got ().",
                     eval(:(@enzyme_reaction begin
            substrates: S[C]
            products:   P[C]
            oligomeric_state: ()
        end)))
        # Two values after the label, not one tuple.
        @test_throws("@enzyme_reaction: `oligomeric_state:` takes one Int or a tuple of " *
                     "one or more Ints, such as `(1, 2, 4)`; got 1, 2.", eval(Meta.parse("""
            @enzyme_reaction begin
                substrates: S[C]
                products:   P[C]
                oligomeric_state: 1, 2
            end""")))
    end

    @testset "@enzyme_reaction shared_catalytic_site" begin
        rxn = @enzyme_reaction begin
            substrates: A[C], B[C]
            products:   P[C], Q[C]
            shared_catalytic_site: (A, P), (B, Q)
        end
        @test EnzymeRates.shared_catalytic_site(rxn) == [(:A, :P), (:B, :Q)]

        # Single pair.
        rxn1 = @enzyme_reaction begin
            substrates: A[C], B[C]
            products:   P[C], Q[C]
            shared_catalytic_site: (P, A)
        end
        @test EnzymeRates.shared_catalytic_site(rxn1) == [(:A, :P)]

        # Malformed: bare (unparenthesized) names error.
        @test_throws Exception eval(:(@enzyme_reaction begin
            substrates: A[C]
            products:   P[C]
            shared_catalytic_site: A, P
        end))

        # Malformed: three-name tuple errors.
        @test_throws Exception eval(:(@enzyme_reaction begin
            substrates: A[C], B[C]
            products:   P[C], Q[C]
            shared_catalytic_site: (A, P, Q)
        end))

        # Malformed: non-Symbol element in a pair errors.
        @test_throws Exception eval(:(@enzyme_reaction begin
            substrates: A[C]
            products:   P[C]
            shared_catalytic_site: (A, 5)
        end))
    end

    @testset "a reaction prints as the @enzyme_reaction text that builds it" begin
        # The reaction shapes the tests use: one formula per reactant, competitive
        # inhibitors with a reactant among them, allosteric regulators with types and
        # multiplicities, a name in both regulator roles, declared in either order, one
        # or several catalytic multiplicities and shared catalytic sites. The printed
        # text builds an equal reaction.
        reactions = [
            @enzyme_reaction(begin
                substrates: S[C]
                products:   P[C]
            end),
            @enzyme_reaction(begin
                substrates: S[C6H12O6], ATP[C10H16N5O13P3]
                products:   P[C6H13O9P], ADP[C10H15N5O10P2]
                competitive_inhibitors: I, P
                allosteric_regulators: X::Activator, Y(1, 2)::Inhibitor, Z, W(4)
                oligomeric_state: 4
            end),
            @enzyme_reaction(begin
                substrates: A[C2H3], B[NP]
                products:   P[C2N], Q[H3P]
                competitive_inhibitors: ATP, I(1, 2)
                allosteric_regulators: ATP(1)
                oligomeric_state: (4, 1, 2)
                shared_catalytic_site: (A, P), (Q, B)
            end),
            @enzyme_reaction(begin
                substrates: S[C]
                products:   P[C]
                allosteric_regulators: ATP(1)
                competitive_inhibitors: ATP
                oligomeric_state: 2
            end)]
        for rxn in reactions
            @test eval(Meta.parse(sprint(show, rxn))) == rxn
        end
        # A regulator whose multiplicities are the catalytic ones prints bare (W), any
        # other with its multiplicities (Y).
        @test sprint(show, reactions[2]) == """
            @enzyme_reaction begin
                substrates: ATP[C10H16N5O13P3], S[C6H12O6]
                products: ADP[C10H15N5O10P2], P[C6H13O9P]
                competitive_inhibitors: I, P
                allosteric_regulators: W, X::Activator, Y(1, 2)::Inhibitor, Z
                oligomeric_state: 4
            end"""
        @test sprint(show, reactions[3]) == """
            @enzyme_reaction begin
                substrates: A[C2H3], B[NP]
                products: P[C2N], Q[H3P]
                competitive_inhibitors: ATP, I(1, 2)
                allosteric_regulators: ATP(1)
                oligomeric_state: (1, 2, 4)
                shared_catalytic_site: (A, P), (B, Q)
            end"""
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
        @test sum(length, EnzymeRates.steps(EnzymeRates.Mechanism(m))) == 2
        @test Set(_testhelper_enzyme_forms(_testhelper_flat_steps(m))) == Set([:E, :ES])
        @test Set(metabolites(m)) == Set([:S, :P])
    end

    @testset "Elementary steps" begin
        # The decomposed grammar's parser requires an enzyme form on each
        # step side at macro-expansion time. Wrap in `eval(:(...))` so the
        # macro expansion happens at runtime where `@test_throws` can catch
        # the LoadError that wraps the parser's exception.
        @test_throws "no enzyme-form term" eval(:(@enzyme_mechanism begin
            substrates: S
            products:   P
            steps: begin
                S <--> E(S)
            end
        end))
    end

    @testset "a call with fewer than three arguments is not a step" begin
        # `f(x)` and `-E` parse as two-argument calls; the parser reports the line
        # instead of reading a right-hand side the call lacks.
        for line in (:(f(x)), :(-E))
            @test_throws "@enzyme_mechanism: expected step or step-group; got $line" eval(
                :(@enzyme_mechanism begin
                    substrates: S
                    products:   P
                    steps: begin
                        E + S ⇌ E(S)
                        $line
                    end
                end))
        end
    end

    @testset "::Mod role tag: product that also competitively inhibits" begin
        m = @enzyme_mechanism begin
            substrates: S
            products:   P
            modifiers: P
            steps: begin
                E + S <--> E(S)
                E(S) <--> E(P)
                E(P) <--> E + P
                E + P::Mod <--> E(P::Mod)
            end
        end
        forms = _testhelper_enzyme_forms(_testhelper_flat_steps(m))
        @test :EP in forms        # product-bound form (catalytic release)
        @test :EPinh in forms     # inhibitor-bound form — DISTINCT from :EP
        @test :EP != :EPinh
        # rate_equation must reference concs.P for the inhibitor term, not concs.Pmod:
        re_str = rate_equation_string(m)
        @test occursin("P", re_str)             # real metabolite name appears
        @test occursin("(; S, P) = concs", re_str)  # concentration is P, never Pmod
        # The mod marker never names a concentration.
        @test !occursin("Pmod = concs", re_str)
        _testhelper_check_law_against_oracle(m)
    end

    @testset "dual-role names bind by their catalytic role unless tagged ::Mod" begin
        _testhelper_bound_types(m) = Dict(
            EnzymeRates.name(EnzymeRates.to_species(s)) =>
                typeof(EnzymeRates.bound_metabolite(s))
            for g in EnzymeRates.steps(m) for s in g if EnzymeRates.is_binding(s))
        # A substrate also declared a modifier binds as the substrate in a bare
        # catalytic step; `E(A::Mod)` writes its modifier form.
        m = EnzymeRates.Mechanism(@enzyme_mechanism begin
            substrates: A
            products:   P
            modifiers: A
            steps: begin
                E + A ⇌ E(A)
                E(A) <--> E(P)
                E(P) ⇌ E + P
                E + A::Mod ⇌ E(A::Mod)
            end
        end)
        types = _testhelper_bound_types(m)
        @test types[:EA] === EnzymeRates.Substrate
        @test types[:EAinh] === EnzymeRates.CompetitiveInhibitor
        @test types[:EP] === EnzymeRates.Product

        # A substrate and a product also declared as allosteric regulators bind
        # as the substrate and the product in catalytic steps; as regulators
        # they bind only at their regulatory site.
        am = EnzymeRates.AllostericMechanism(@allosteric_mechanism begin
            substrates: A
            products:   P
            allosteric_regulators: A, P
            catalytic_multiplicity: 2
            steps: begin
                E + A ⇌ E(A)          :: EqualAI
                E(A) <--> E(P)        :: NonequalAI
                E(P) ⇌ E + P          :: EqualAI
            end
            regulatory_site(2): A::OnlyI, P::OnlyA
        end)
        types = _testhelper_bound_types(am)
        @test types[:EA] === EnzymeRates.Substrate
        @test types[:EP] === EnzymeRates.Product
        site = only(EnzymeRates.regulatory_sites(am))
        @test Set(EnzymeRates.name.(EnzymeRates.ligands(site))) == Set([:A, :P])
        @test all(l -> l isa EnzymeRates.AllostericRegulator, EnzymeRates.ligands(site))
    end

    @testset "an allosteric regulator binds only at its regulatory site" begin
        # R is declared only in `allosteric_regulators:`: a free term, a bound
        # metabolite or a residual entry of a catalytic step names it in error.
        @test_throws "`R` is an allosteric regulator" eval(:(@allosteric_mechanism begin
            substrates: S
            products:   P
            allosteric_regulators: R
            steps: begin
                E + S ⇌ E(S)      :: EqualAI
                E(S) <--> E(P)    :: EqualAI
                E(P) ⇌ E + P      :: EqualAI
                E + R ⇌ E(R)      :: EqualAI
            end
            regulatory_site(1): R::OnlyA
        end))
        @test_throws "`R` is an allosteric regulator" eval(:(@allosteric_mechanism begin
            substrates: S
            products:   P
            allosteric_regulators: R
            steps: begin
                E + S ⇌ E(S)        :: EqualAI
                E(S) <--> E(P, R)   :: EqualAI
                E(P, R) ⇌ E + P     :: EqualAI
            end
            regulatory_site(1): R::OnlyA
        end))
        @test_throws "got `R`" eval(:(@allosteric_mechanism begin
            substrates: S
            products:   P
            allosteric_regulators: R
            steps: begin
                E + S ⇌ E(S)                       :: EqualAI
                E(S) <--> Estar(; residual = S - R) :: EqualAI
                Estar(; residual = S - R) ⇌ E + P   :: EqualAI
            end
            regulatory_site(1): R::OnlyA
        end))
        # Tagged `R::Mod`, it is neither a substrate nor a product to bind as a modifier.
        # The macro raises the error while it expands, so `eval` wraps it in a LoadError:
        # match the text.
        @test_throws(
            "@allosteric_mechanism: `R::Mod` marks a substrate or product bound as a " *
            "modifier, and `R` is neither; write a name under `modifiers:` bare.",
            eval(:(@allosteric_mechanism begin
            substrates: S
            products:   P
            allosteric_regulators: R
            steps: begin
                E + S ⇌ E(S)              :: EqualAI
                E(S) <--> E(P)            :: EqualAI
                E(P) ⇌ E + P              :: EqualAI
                E + R::Mod ⇌ E(R::Mod)    :: EqualAI
            end
            regulatory_site(1): R::OnlyA
        end)))
        # Listed in both, R binds its regulatory site and, as a modifier, the catalytic
        # site.
        dual = EnzymeRates.AllostericMechanism(@allosteric_mechanism begin
            substrates: S
            products:   P
            modifiers: R
            allosteric_regulators: R
            steps: begin
                E + S ⇌ E(S)              :: EqualAI
                E(S) <--> E(P)            :: EqualAI
                E(P) ⇌ E + P              :: EqualAI
                E + R ⇌ E(R)    :: EqualAI
            end
            regulatory_site(1): R::OnlyA
        end)
        @test EnzymeRates.CompetitiveInhibitor(:R) in
              [EnzymeRates.bound_metabolite(g[1]) for g in EnzymeRates.steps(dual)]
        @test EnzymeRates.AllostericRegulator(:R) in
              EnzymeRates.ligands(only(EnzymeRates.regulatory_sites(dual)))
    end

    @testset "X::Mod names a declared modifier" begin
        # A substrate or a product written `X::Mod` must also be listed under
        # `modifiers:`; the reaction then lists it as a competitive inhibitor.
        for x in (:S, :P)
            @test_throws(
                "@enzyme_mechanism: `$x::Mod` names `$x`, which `modifiers:` does " *
                "not declare; add `$x` to `modifiers:` to let it bind the catalytic " *
                "site as a modifier.", eval(:(@enzyme_mechanism begin
                substrates: S
                products:   P
                steps: begin
                    E + S ⇌ E(S)
                    E(S) <--> E(P)
                    E(P) ⇌ E + P
                    E + $x::Mod ⇌ E($x::Mod)
                end
            end)))
        end
        declared = EnzymeRates.Mechanism(@enzyme_mechanism begin
            substrates: S
            products:   P
            modifiers: P
            steps: begin
                E + S ⇌ E(S)
                E(S) <--> E(P)
                E(P) ⇌ E + P
                E + P::Mod ⇌ E(P::Mod)
            end
        end)
        @test EnzymeRates.regulators(EnzymeRates.reaction(declared)) ==
              [EnzymeRates.RegulatorMults(EnzymeRates.CompetitiveInhibitor(:P), [1])]
    end

    @testset "a step side holds one enzyme form" begin
        # A side with two enzyme forms is rejected.
        @test_throws "more than one enzyme-form term" eval(:(@enzyme_mechanism begin
            substrates: A, B
            products:   P
            steps: begin
                E(A) + E(B) <--> E(A, B)
                E(A, B) <--> E + P
            end
        end))
    end

    @testset "a residual adds substrates and subtracts products" begin
        # `P - S` adds a product and subtracts a substrate; `P` and `-S` each do one of
        # the two; `S - I` subtracts a modifier.
        for residual in (:(P - S), :P, :(-S), :(S - I))
            @test_throws "adds substrates and subtracts products" eval(
                :(@enzyme_mechanism begin
                    substrates: S
                    products:   P
                    modifiers: I
                    steps: begin
                        E + S ⇌ E(S)
                        E(S) <--> Estar(; residual = $residual)
                        Estar(; residual = $residual) ⇌ E + P
                    end
                end))
        end
    end

    @testset "::Mod is never a step tag" begin
        # `E(S::Mod) ⇌ E + S::Mod :: EqualAI` parses as `(S::Mod)::EqualAI`: the state
        # after `::Mod` is the step's tag, and the step reads like its binding direction.
        bind_inh = EnzymeRates.AllostericMechanism(@allosteric_mechanism begin
            substrates: S
            products:   P
            modifiers: S
            steps: begin
                E + S ⇌ E(S)             :: EqualAI
                E(S) <--> E(P)           :: EqualAI
                E(P) ⇌ E + P             :: EqualAI
                E + S::Mod ⇌ E(S::Mod)   :: EqualAI
            end
        end)
        release_inh = EnzymeRates.AllostericMechanism(@allosteric_mechanism begin
            substrates: S
            products:   P
            modifiers: S
            steps: begin
                E + S ⇌ E(S)             :: EqualAI
                E(S) <--> E(P)           :: EqualAI
                E(P) ⇌ E + P             :: EqualAI
                E(S::Mod) ⇌ E + S::Mod   :: EqualAI
            end
        end)
        @test release_inh == bind_inh

        # Without a state the step is missing its tag: the role tag that ends it is no
        # state.
        @test_throws("@allosteric_mechanism: step `E(S::Mod) ⇌ E + S::Mod` is missing " *
                     "its tag, `:: EqualAI`, `:: NonequalAI` or `:: OnlyA`; add it after " *
                     "the step expression.", eval(:(@allosteric_mechanism begin
            substrates: S
            products:   P
            modifiers: S
            steps: begin
                E + S ⇌ E(S)           :: EqualAI
                E(S) <--> E(P)         :: EqualAI
                E(P) ⇌ E + P           :: EqualAI
                E(S::Mod) ⇌ E + S::Mod
            end
        end)))

        # A plain mechanism may release a modifier copy as well as bind it: the role tag
        # that ends the release is no step tag.
        bind_plain = EnzymeRates.Mechanism(@enzyme_mechanism begin
            substrates: S
            products:   P
            modifiers: S
            steps: begin
                E + S ⇌ E(S)
                E(S) <--> E(P)
                E(P) ⇌ E + P
                E + S::Mod ⇌ E(S::Mod)
            end
        end)
        release_plain = EnzymeRates.Mechanism(@enzyme_mechanism begin
            substrates: S
            products:   P
            modifiers: S
            steps: begin
                E + S ⇌ E(S)
                E(S) <--> E(P)
                E(P) ⇌ E + P
                E(S::Mod) ⇌ E + S::Mod
            end
        end)
        @test release_plain == bind_plain

        # A step tag that is not a name is rejected.
        @test_throws "must be a Symbol" eval(:(@allosteric_mechanism begin
            substrates: S
            products:   P
            steps: begin
                E + S ⇌ E(S)    :: EqualAI
                E(S) <--> E(P)  :: EqualAI
                E(P) ⇌ E + P    :: 3
            end
        end))
    end
end

@testset "DSL: modifiers, ::Mod and regulatory_site(n)" begin
    @testset "a modifier binds bare, a reactant's modifier copy as X::Mod" begin
        # I binds only as a modifier; P is a product that also binds as a modifier.
        m = @enzyme_mechanism begin
            substrates: S
            products:   P
            modifiers:  I, P
            steps: begin
                E + S ⇌ E(S)
                E(S) <--> E(P)
                E(P) ⇌ E + P
                E(S) + I ⇌ E(S, I)
                E + P::Mod ⇌ E(P::Mod)
            end
        end
        @test m isa EnzymeMechanism
        @test EnzymeRates.fitted_params(m) == (:K_EIS_to_ES_I, :K_EP_to_E_P,
            :K_EPmod_to_E_Pmod, :K_ES_to_E_S, :k_ES_to_EP)
        @test metabolites(m) == (:S, :P, :I)
        _testhelper_check_law_against_oracle(m)
    end

    @testset "the allosteric grammar: one tag per statement, ligands at their sites" begin
        am = @allosteric_mechanism begin
            substrates: A, B
            products:   P, Q
            modifiers:  I, P
            allosteric_regulators: X, Y
            catalytic_multiplicity: 2
            steps: begin
                E + A ⇌ E(A)                              :: EqualAI
                E(A) + B ⇌ E(A, B)                        :: NonequalAI
                E(A, B) <--> E(P, Q)                      :: OnlyA
                E(P, Q) ⇌ E(Q) + P                        :: EqualAI
                E(Q) ⇌ E + Q                              :: EqualAI
                (E + I ⇌ E(I), E(A) + I ⇌ E(A, I))        :: EqualAI
                E(Q) + P::Mod ⇌ E(P::Mod, Q)              :: EqualAI
            end
            regulatory_site(2): X::OnlyA, Y::OnlyI
        end
        @test am isa AllostericEnzymeMechanism
        @test EnzymeRates.catalytic_multiplicity(am) == 2
        @test EnzymeRates.fitted_params(am) == (:K_A_EAB_to_EA_B, :K_EA_to_E_A,
            :K_EI_to_E_I, :K_EPQ_to_EQ_P, :K_EPmodQ_to_EQ_Pmod, :K_EQ_to_E_Q,
            :K_I_EAB_to_EA_B, :k_A_EAB_to_EPQ, :K_A_Xreg, :K_I_Yreg, :L)
        @test metabolites(am) == (:A, :B, :P, :Q, :I, :X, :Y)
        _testhelper_check_law_against_oracle(am)
    end

    @testset "catalytic_multiplicity: defaults to 1" begin
        am = @allosteric_mechanism begin
            substrates: S
            products:   P
            allosteric_regulators: X
            steps: begin
                E + S ⇌ E(S)      :: EqualAI
                E(S) <--> E(P)    :: OnlyA
                E(P) ⇌ E + P      :: EqualAI
            end
            regulatory_site(1): X::OnlyI
        end
        @test EnzymeRates.catalytic_multiplicity(am) == 1
        @test EnzymeRates.fitted_params(am) ==
              (:K_EP_to_E_P, :K_ES_to_E_S, :k_A_ES_to_EP, :K_I_Xreg, :L)
        _testhelper_check_law_against_oracle(am)
    end

    @testset "a name can be a modifier and an allosteric regulator" begin
        # R binds the catalytic site as a modifier and its regulatory site as a ligand.
        dual = @allosteric_mechanism begin
            substrates: S
            products:   P
            modifiers:  R
            allosteric_regulators: R
            steps: begin
                E + S ⇌ E(S)      :: EqualAI
                E(S) <--> E(P)    :: EqualAI
                E(P) ⇌ E + P      :: EqualAI
                E + R ⇌ E(R)      :: EqualAI
            end
            regulatory_site(1): R::OnlyA
        end
        @test EnzymeRates.fitted_params(dual) ==
              (:K_EP_to_E_P, :K_ER_to_E_R, :K_ES_to_E_S, :k_ES_to_EP, :K_A_Rreg, :L)
        @test metabolites(dual) == (:S, :P, :R)
        _testhelper_check_law_against_oracle(dual)
    end

    @testset "the walker returns each block as written" begin
        # Steps keep their written order and direction, each statement is one kinetic
        # group label, and a metabolite is a modifier when written bare under
        # `modifiers:` or as a reactant's X::Mod.
        M, F = EnzymeRates.Met, EnzymeRates.RawForm
        s, p, i, smod = M(:S, false), M(:P, false), M(:I, true), M(:S, true)
        raw = EnzymeRates._raw_mechanism(quote
            substrates: S
            products:   P
            modifiers:  I, S
            steps: begin
                E(S) ⇌ E + S
                (E + I ⇌ E(I), E(S) + I ⇌ E(I, S))
                E(S) <--> Estar(; residual = S - P) + P
                E + S::Mod ⇌ E(S::Mod)
            end
        end, false)
        @test raw == EnzymeRates.RawMechanism([:S], [:P], [:I, :S], [
            EnzymeRates.RawStep(F(:E, [s]), F(:E, M[]), M[], [s], true),
            EnzymeRates.RawStep(F(:E, M[]), F(:E, [i]), [i], M[], true),
            EnzymeRates.RawStep(F(:E, [s]), F(:E, [i, s]), [i], M[], true),
            EnzymeRates.RawStep(F(:E, [s]), F(:Estar, M[], [:S], [:P]), M[], [p], false),
            EnzymeRates.RawStep(F(:E, M[]), F(:E, [smod]), [smod], M[], true)],
            [1, 2, 2, 3, 4])
        # An allosteric block adds a tag per group label, the catalytic multiplicity
        # (1 when omitted) and its sites, each ligand list as written.
        allo = EnzymeRates._raw_mechanism(quote
            substrates: S
            products:   P
            allosteric_regulators: X, Y
            steps: begin
                E + S ⇌ E(S)              :: EqualAI
                E(S) <--> E(P)            :: OnlyA
                (E(P) ⇌ E + P,)           :: NonequalAI
            end
            regulatory_site(2): Y::OnlyI, X::OnlyA
        end, true)
        @test allo == EnzymeRates.RawAllostericMechanism(
            EnzymeRates.RawMechanism([:S], [:P], Symbol[], [
                EnzymeRates.RawStep(F(:E, M[]), F(:E, [s]), [s], M[], true),
                EnzymeRates.RawStep(F(:E, [s]), F(:E, [p]), M[], M[], false),
                EnzymeRates.RawStep(F(:E, [p]), F(:E, M[]), M[], [p], true)],
                [1, 2, 3]),
            [:EqualAI, :OnlyA, :NonequalAI], 1,
            [EnzymeRates.RawSite([:Y, :X], 2, [:OnlyI, :OnlyA])])
    end

    @testset "::Mod marks a substrate or product under modifiers: only" begin
        # A modifier-only name binds bare, so ::Mod on it is an error.
        @test_throws("@enzyme_mechanism: `I::Mod` marks a substrate or product bound " *
                     "as a modifier, and `I` is neither; write a name under " *
                     "`modifiers:` bare.", eval(:(@enzyme_mechanism begin
            substrates: S
            products:   P
            modifiers:  I
            steps: begin
                E + S ⇌ E(S)
                E(S) <--> E(P)
                E(P) ⇌ E + P
                E + I::Mod ⇌ E(I::Mod)
            end
        end)))
        # A reactant written X::Mod must be listed under modifiers:.
        @test_throws("@enzyme_mechanism: `P::Mod` names `P`, which `modifiers:` does " *
                     "not declare; add `P` to `modifiers:` to let it bind the catalytic " *
                     "site as a modifier.", eval(:(@enzyme_mechanism begin
            substrates: S
            products:   P
            steps: begin
                E + S ⇌ E(S)
                E(S) <--> E(P)
                E(P) ⇌ E + P
                E + P::Mod ⇌ E(P::Mod)
            end
        end)))
        # An allosteric regulator binds only at its site, bare or tagged.
        @test_throws "`R` is an allosteric regulator" eval(:(@allosteric_mechanism begin
            substrates: S
            products:   P
            allosteric_regulators: R
            steps: begin
                E + S ⇌ E(S)      :: EqualAI
                E(S) <--> E(P)    :: EqualAI
                E(P) ⇌ E + P      :: EqualAI
                E + R ⇌ E(R)      :: EqualAI
            end
            regulatory_site(1): R::OnlyA
        end))
        # ::Mod is the one role tag, also where a role tag ends a step.
        for line in (:(E + P::Inh ⇌ E(P::Inh)), :(E(P) ⇌ E + P::Inh))
            @test_throws("@enzyme_mechanism: unknown role tag ::Inh on P; the one role " *
                         "tag is ::Mod, which marks a substrate or product bound as a " *
                         "modifier.", eval(:(@enzyme_mechanism begin
                substrates: S
                products:   P
                modifiers:  P
                steps: begin
                    E + S ⇌ E(S)
                    E(S) <--> E(P)
                    E(P) ⇌ E + P
                    $line
                end
            end)))
        end
    end

    @testset "every allosteric regulator sits at exactly one regulatory site" begin
        # The mechanism `@allosteric_mechanism` builds from the block below with the site
        # lines `sites` appended.
        declarations = quote
            substrates: S
            products:   P
            allosteric_regulators: X, Y
            steps: begin
                E + S ⇌ E(S)      :: EqualAI
                E(S) <--> E(P)    :: OnlyA
                E(P) ⇌ E + P      :: EqualAI
            end
        end
        _testhelper_sites(sites...) = eval(:(@allosteric_mechanism $(
            Expr(:block, declarations.args..., sites...))))
        @test_throws("@allosteric_mechanism: allosteric regulator `Y` sits at no " *
                     "regulatory site; list it on one `regulatory_site(n):` line, as " *
                     "in `regulatory_site(2): Y::OnlyA`.",
                     _testhelper_sites(:(regulatory_site(2): X::OnlyA)))
        @test_throws("@allosteric_mechanism: allosteric regulator `X` is listed 2 times " *
                     "on `regulatory_site(n):` lines; an allosteric regulator sits at " *
                     "exactly one regulatory site.",
                     _testhelper_sites(:(regulatory_site(2): X::OnlyA, Y::OnlyI),
                                       :(regulatory_site(4): X::OnlyA)))
        @test_throws("@allosteric_mechanism: allosteric regulator `X` is listed 2 times",
                     _testhelper_sites(:(regulatory_site(2): X::OnlyA, Y::OnlyI,
                                                             X::OnlyI)))
        @test_throws("@allosteric_mechanism: ligand `Z` on a `regulatory_site(n):` line " *
                     "is not declared in `allosteric_regulators:`.",
                     _testhelper_sites(:(regulatory_site(2): X::OnlyA, Y::OnlyI,
                                                             Z::OnlyI)))
        @test_throws("@allosteric_mechanism: `allosteric_regulators:` lists `X` more " *
                     "than once.", eval(:(@allosteric_mechanism begin
            substrates: S
            products:   P
            allosteric_regulators: X, X
            steps: begin
                E + S ⇌ E(S)      :: EqualAI
            end
            regulatory_site(2): X::OnlyA
        end)))
        # A site lists each ligand with its tag after a positive multiplicity.
        @test_throws("@allosteric_mechanism: `regulatory_site(n):` lists each ligand " *
                     "with its tag, as in `regulatory_site(2): A::OnlyA, I::OnlyI`; " *
                     "got `X`.", _testhelper_sites(:(regulatory_site(2): X, Y::OnlyI)))
        @test_throws("`regulatory_site` multiplicity must be a positive Int, got 0.",
                     _testhelper_sites(:(regulatory_site(0): X::OnlyA, Y::OnlyI)))
        @test_throws("@allosteric_mechanism: `regulatory_site(n):` lists each ligand " *
                     "with its tag, as in `regulatory_site(2): A::OnlyA, I::OnlyI`; " *
                     "got a `begin ... end` block.",
                     _testhelper_sites(:(regulatory_site(2): begin
                                             ligands: X, Y
                                         end)))
        # The allosteric regulators list names; each tag sits on a site line. A tagged
        # regulator is rejected whether or not a site line lists it.
        tagged = "@allosteric_mechanism: `allosteric_regulators:` lists names; got " *
                 "`X::OnlyA`. A ligand's tag goes on its `regulatory_site(n):` line, as " *
                 "in `regulatory_site(2): X::OnlyA`."
        @test_throws tagged eval(:(@allosteric_mechanism begin
            substrates: S
            products:   P
            allosteric_regulators: X::OnlyA
            steps: begin
                E + S ⇌ E(S)      :: EqualAI
            end
            regulatory_site(2): X::OnlyA
        end))
        @test_throws tagged eval(:(@allosteric_mechanism begin
            substrates: S
            products:   P
            allosteric_regulators: X::OnlyA
            steps: begin
                E + S ⇌ E(S)      :: EqualAI
            end
        end))
    end

    @testset "allosteric steps carry one tag per statement, plain steps none" begin
        @test_throws("@allosteric_mechanism: step `E(S) <--> E(P)` is missing its tag, " *
                     "`:: EqualAI`, `:: NonequalAI` or `:: OnlyA`; add it after the step " *
                     "expression.", eval(:(@allosteric_mechanism begin
            substrates: S
            products:   P
            steps: begin
                E + S ⇌ E(S)      :: EqualAI
                E(S) <--> E(P)
                E(P) ⇌ E + P      :: EqualAI
            end
        end)))
        @test_throws("@enzyme_mechanism: tag annotation on `E(S) <--> E(P)::EqualAI` is " *
                     "not allowed", eval(:(@enzyme_mechanism begin
            substrates: S
            products:   P
            steps: begin
                E + S ⇌ E(S)
                E(S) <--> E(P)    :: EqualAI
                E(P) ⇌ E + P
            end
        end)))
        # A single step carries its tag after its right-hand side; a parenthesized group
        # is a tuple of steps, `(step,)` for one.
        @test_throws("@allosteric_mechanism: a single step carries its tag after its " *
                     "right-hand side, as in `E + S ⇌ E(S) :: EqualAI`. Got " *
                     "`(E(P) ⇌ E + P)::EqualAI`.", eval(:(@allosteric_mechanism begin
            substrates: S
            products:   P
            steps: begin
                E + S ⇌ E(S)              :: EqualAI
                E(S) <--> E(P)            :: EqualAI
                (E(P) ⇌ E + P)            :: EqualAI
            end
        end)))
        # A tag after a parenthesized right-hand side tags the step, as after a bare one.
        M, F = EnzymeRates.Met, EnzymeRates.RawForm
        s, p = M(:S, false), M(:P, false)
        @test EnzymeRates._raw_mechanism(quote
            substrates: S
            products:   P
            steps: begin
                E + S ⇌ E(S)        :: EqualAI
                E(S) <--> E(P)      :: OnlyA
                E(P) ⇌ (E + P)      :: EqualAI
            end
        end, true) == EnzymeRates.RawAllostericMechanism(
            EnzymeRates.RawMechanism([:S], [:P], Symbol[], [
                EnzymeRates.RawStep(F(:E, M[]), F(:E, [s]), [s], M[], true),
                EnzymeRates.RawStep(F(:E, [s]), F(:E, [p]), M[], M[], false),
                EnzymeRates.RawStep(F(:E, [p]), F(:E, M[]), M[], [p], true)],
                [1, 2, 3]),
            [:EqualAI, :OnlyA, :EqualAI], 1, EnzymeRates.RawSite[])
        @test_throws "@enzyme_mechanism: a parenthesized group holds no step" eval(
            :(@enzyme_mechanism begin
                substrates: S
                products:   P
                steps: begin
                    E + S ⇌ E(S)
                    ()
                end
            end))
    end

    @testset "labels: the valid ones, each single-valued label once with one value" begin
        plain_labels = "Valid labels: substrates, products, modifiers, steps. " *
                       "Allosteric declarations (`allosteric_regulators:`, " *
                       "`catalytic_multiplicity:`, `regulatory_site(n):`) belong in " *
                       "@allosteric_mechanism."
        for label in (:regulators, :catalytic_steps, :catalytic_inhibitors)
            @test_throws("@enzyme_mechanism: unknown label `$label:`. " * plain_labels,
                         eval(:(@enzyme_mechanism begin
                substrates: S
                products:   P
                $label: I
                steps: begin
                    E + S ⇌ E(S)
                end
            end)))
        end
        allo_labels = "Valid labels: substrates, products, modifiers, " *
                      "allosteric_regulators, catalytic_multiplicity, steps, " *
                      "regulatory_site(n)."
        for label in (:catalytic_steps, :catalytic_inhibitors, :regulators)
            @test_throws("@allosteric_mechanism: unknown label `$label:`. " * allo_labels,
                         eval(:(@allosteric_mechanism begin
                substrates: S
                products:   P
                $label: I
                steps: begin
                    E + S ⇌ E(S)      :: EqualAI
                end
            end)))
        end
        @test_throws("@allosteric_mechanism: `catalytic_multiplicity:` takes one " *
                     "positive Int; got 2, 4.", eval(Meta.parse("""
            @allosteric_mechanism begin
                substrates: S
                products:   P
                catalytic_multiplicity: 2, 4
                steps: begin
                    E + S ⇌ E(S)      :: EqualAI
                end
            end""")))
        @test_throws("@enzyme_mechanism: `steps:` takes a `begin ... end` block of " *
                     "steps; got S.", eval(:(@enzyme_mechanism begin
            substrates: S
            products:   P
            steps: S
        end)))
        @test_throws("@enzyme_mechanism: `steps:` takes one `begin ... end` block; got " *
                     "2 values.", eval(Meta.parse("""
            @enzyme_mechanism begin
                substrates: S
                products:   P
                steps: begin
                    E + S ⇌ E(S)
                end, S
            end""")))
        @test_throws "`steps:` given more than once" eval(:(@enzyme_mechanism begin
            substrates: S
            products:   P
            steps: begin
                E + S ⇌ E(S)
            end
            steps: begin
                E(S) <--> E + P
            end
        end))
        @test_throws "`catalytic_multiplicity:` given more than once" eval(
            :(@allosteric_mechanism begin
                substrates: S
                products:   P
                catalytic_multiplicity: 2
                catalytic_multiplicity: 2
                steps: begin
                    E + S ⇌ E(S)      :: EqualAI
                end
            end))
    end

    @testset "a conformation label: a capital, then lowercase letters and digits" begin
        @test_throws("@enzyme_mechanism: `E_c` is not a conformation label: a " *
                     "conformation label is a capital letter followed by lowercase " *
                     "letters and digits, such as `E`, `Estar` or `F2`.",
                     eval(:(@enzyme_mechanism begin
            substrates: S
            products:   P
            steps: begin
                E_c + S ⇌ E_c(S)
                E_c(S) <--> E_c + P
            end
        end)))
        m = @enzyme_mechanism begin
            substrates: S
            products:   P
            steps: begin
                E + S ⇌ E(S)
                E(S) <--> Estar2(P)
                Estar2(P) ⇌ Estar2 + P
                Estar2 <--> E
            end
        end
        @test EnzymeRates.fitted_params(m) == (:K_ES_to_E_S, :K_Estar2P_to_Estar2_P,
            :k_ES_to_Estar2P, :k_E_to_Estar2, :k_Estar2_to_E)
        _testhelper_check_law_against_oracle(m)
    end

    @testset "constructor errors throw when the expanded code runs" begin
        # The walker reads the block at expansion; the constructors check the mechanism
        # when the macro's code runs, so these throw outside `eval`.
        @test_throws ":OnlyI is rejected for catalytic groups" @allosteric_mechanism begin
            substrates: S
            products:   P
            steps: begin
                E + S ⇌ E(S)      :: EqualAI
                E(S) <--> E(P)    :: OnlyI
                E(P) ⇌ E + P      :: EqualAI
            end
        end
        @test_throws "allo state Bogus must be one of" @allosteric_mechanism begin
            substrates: S
            products:   P
            allosteric_regulators: X
            steps: begin
                E + S ⇌ E(S)      :: EqualAI
                E(S) <--> E(P)    :: EqualAI
                E(P) ⇌ E + P      :: EqualAI
            end
            regulatory_site(2): X::Bogus
        end
        @test_throws "S named as both a substrate and a product" @enzyme_mechanism begin
            substrates: S
            products:   S
            steps: begin
                E + S ⇌ E(S)
            end
        end
    end
end
