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
            regulators: I

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

    @testset "@allosteric_mechanism leaves states and ligands to the constructors" begin
        # A name listed twice in `allosteric_regulators:` is rejected, even when an
        # explicit site would let the last tag silently replace the first.
        @test_throws "lists `A` more than once" eval(:(@allosteric_mechanism begin
            substrates: S
            products:   P
            allosteric_regulators: A::Bogus, A::OnlyA
            catalytic_steps: begin
                E + S ⇌ E(S)      :: EqualAI
                E(S) <--> E(P)    :: EqualAI
                E(P) ⇌ E + P      :: EqualAI
            end
            regulatory_site(multiplicity = 2): begin
                ligands: A
            end
        end))
        @test_throws "lists `A` more than once" eval(:(@allosteric_mechanism begin
            substrates: S
            products:   P
            allosteric_regulators: A::OnlyA, A::OnlyA
            catalytic_steps: begin
                E + S ⇌ E(S)      :: EqualAI
                E(S) <--> E(P)    :: EqualAI
                E(P) ⇌ E + P      :: EqualAI
            end
        end))

        # An unknown regulator state is reported by RegulatorySite when the expansion
        # runs.
        @test_throws "allo state Bogus must be one of" eval(:(@allosteric_mechanism begin
            substrates: S
            products:   P
            allosteric_regulators: X::Bogus
            catalytic_steps: begin
                E + S ⇌ E(S)      :: EqualAI
                E(S) <--> E(P)    :: EqualAI
                E(P) ⇌ E + P      :: EqualAI
            end
        end))

        # Unknown and :OnlyI catalytic states are reported by AllostericMechanism.
        for state in (:Bogus, :OnlyI)
            @test_throws "catalytic group" eval(:(@allosteric_mechanism begin
                substrates: S
                products:   P
                catalytic_steps: begin
                    E + S ⇌ E(S)      :: EqualAI
                    E(S) <--> E(P)    :: EqualAI
                    E(P) ⇌ E + P     :: $state
                end
            end))
        end

        # A ligand on two regulatory sites is reported by AllostericMechanism.
        @test_throws "appears in two distinct regulatory sites" eval(
            :(@allosteric_mechanism begin
                substrates: S
                products:   P
                allosteric_regulators: A::OnlyA
                catalytic_steps: begin
                    E + S ⇌ E(S)      :: EqualAI
                    E(S) <--> E(P)    :: EqualAI
                    E(P) ⇌ E + P      :: EqualAI
                end
                regulatory_site(multiplicity = 2): begin
                    ligands: A
                end
                regulatory_site(multiplicity = 4): begin
                    ligands: A
                end
            end))
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
                catalytic_steps: begin
                    E + S ⇌ E(S)      :: EqualAI
                    E(S) <--> E(P)    :: EqualAI
                    E(P) ⇌ E + P      :: EqualAI
                end
            end))
    end

    @testset "@enzyme_mechanism rejects each allosteric-only label" begin
        for line in (:(allosteric_regulators: A::OnlyA),
                     :(catalytic_inhibitors: I),
                     :(catalytic_multiplicity: 2),
                     :(catalytic_steps: begin
                           E + S ⇌ E(S)
                       end),
                     :(regulatory_site(multiplicity = 2): begin
                           ligands: A
                       end))
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
        @test_throws "@allosteric_mechanism: `catalytic_steps:` not specified" eval(
            :(@allosteric_mechanism begin
                substrates: S
                products:   P
            end))

        # A repeated `catalytic_steps:` line names the label.
        @test_throws "`catalytic_steps:` given more than once" eval(
            :(@allosteric_mechanism begin
                substrates: S
                products:   P
                catalytic_steps: begin
                    E + S ⇌ E(S)      :: EqualAI
                    E(S) <--> E(P)    :: EqualAI
                    E(P) ⇌ E + P      :: EqualAI
                end
                catalytic_steps: begin
                    E + S ⇌ E(S)      :: EqualAI
                    E(S) <--> E(P)    :: EqualAI
                    E(P) ⇌ E + P      :: EqualAI
                end
            end))
    end

    @testset "regulatory_site takes exactly one multiplicity" begin
        # Repeated, missing, misnamed and positional arguments are all errors.
        for site in (:(regulatory_site(multiplicity = 2, multiplicity = 4)),
                     :(regulatory_site()),
                     :(regulatory_site(size = 2)),
                     :(regulatory_site(2)))
            @test_throws "exactly one `multiplicity = N`" eval(
                :(@allosteric_mechanism begin
                    substrates: S
                    products:   P
                    allosteric_regulators: A::OnlyA
                    catalytic_steps: begin
                        E + S ⇌ E(S)      :: EqualAI
                        E(S) <--> E(P)    :: EqualAI
                        E(P) ⇌ E + P      :: EqualAI
                    end
                    $site: begin
                        ligands: A
                    end
                end))
        end
    end

    @testset "a parenthesized single step carries a group tag" begin
        # The @allosteric_mechanism docstring's example: its last step is a
        # parenthesized one-step group, the same mechanism as the bare tagged step.
        parenthesized = EnzymeRates.AllostericMechanism(@allosteric_mechanism begin
            substrates: F6P
            products:   F16BP
            catalytic_multiplicity: 2
            allosteric_regulators: A::OnlyA, I::OnlyI

            catalytic_steps: begin
                E + F6P ⇌ E(F6P)        :: EqualAI
                E(F6P) <--> E(F16BP)    :: EqualAI
                (E(F16BP) ⇌ E + F16BP)  :: EqualAI
            end

            regulatory_site(multiplicity = 4): begin
                ligands: A
            end
            regulatory_site(multiplicity = 4): begin
                ligands: I
            end
        end)
        bare = EnzymeRates.AllostericMechanism(@allosteric_mechanism begin
            substrates: F6P
            products:   F16BP
            catalytic_multiplicity: 2
            allosteric_regulators: A::OnlyA, I::OnlyI

            catalytic_steps: begin
                E + F6P ⇌ E(F6P)        :: EqualAI
                E(F6P) <--> E(F16BP)    :: EqualAI
                E(F16BP) ⇌ E + F16BP    :: EqualAI
            end

            regulatory_site(multiplicity = 4): begin
                ligands: A
            end
            regulatory_site(multiplicity = 4): begin
                ligands: I
            end
        end)
        @test parenthesized == bare
        @test length(EnzymeRates.steps(parenthesized)) == 3

        # A plain mechanism has no group tags, parenthesized or not.
        @test_throws "tag annotation" eval(:(@enzyme_mechanism begin
            substrates: S
            products:   P
            steps: begin
                E + S ⇌ E(S)
                E(S) <--> E(P)
                (E(P) ⇌ E + P) :: EqualAI
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
    end

    @testset "@enzyme_reaction lists substrates before products" begin
        # A substrate and a product may share a name; the constructor's name sort keeps
        # ties in input order, so label order must not decide it.
        products_first = @enzyme_reaction begin
            products:   S[C]
            substrates: S[C]
        end
        substrates_first = @enzyme_reaction begin
            substrates: S[C]
            products:   S[C]
        end
        @test products_first == substrates_first
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
        # Two allosteric ATPs still rejected.
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
        @test_throws Exception eval(:(@enzyme_reaction begin
            substrates: A[C6H12O6]
            products: B[C6H12O6]
            dead_end_inhibitors: X::Inhibitor
        end))
    end

    @testset "@enzyme_reaction rejects bare `regulators:` label" begin
        # The @enzyme_reaction grammar requires `competitive_inhibitors:`,
        # `dead_end_inhibitors:`, or `allosteric_regulators:`. A bare
        # `regulators:` label must be reported as unknown.
        @test_throws "unknown label `regulators:`" eval(:(@enzyme_reaction begin
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
                catalytic_steps: begin
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

        # Explicit allowed_catalytic_multiplicities tuple
        rxn3 = @enzyme_reaction begin
            substrates: S[C]
            products: P[C]
            allowed_catalytic_multiplicities: (1, 2, 4)
        end
        @test EnzymeRates.allowed_catalytic_multiplicities(rxn3) == [1, 2, 4]
    end

    @testset "@enzyme_reaction sets the catalytic multiplicities once" begin
        # `oligomeric_state:` and `allowed_catalytic_multiplicities:` set the same
        # value, so a second line of either label errors, in either order.
        @test_throws "`allowed_catalytic_multiplicities:` sets the catalytic" eval(
            :(@enzyme_reaction begin
                substrates: S[C]
                products:   P[C]
                oligomeric_state: 2
                allowed_catalytic_multiplicities: (1, 2)
            end))
        @test_throws "`oligomeric_state:` sets the catalytic" eval(
            :(@enzyme_reaction begin
                substrates: S[C]
                products:   P[C]
                allowed_catalytic_multiplicities: (1, 2)
                oligomeric_state: 2
            end))
        @test_throws "`allowed_catalytic_multiplicities:` sets the catalytic" eval(
            :(@enzyme_reaction begin
                substrates: S[C]
                products:   P[C]
                allowed_catalytic_multiplicities: (1, 2)
                allowed_catalytic_multiplicities: (4,)
            end))

        # `oligomeric_state:` takes one multiplicity, bare or as a one-element tuple.
        one_tuple = eval(:(@enzyme_reaction begin
            substrates: S[C]
            products:   P[C]
            oligomeric_state: (2,)
        end))
        @test EnzymeRates.allowed_catalytic_multiplicities(one_tuple) == [2]
        @test_throws "`oligomeric_state:` takes a single Int" eval(
            :(@enzyme_reaction begin
                substrates: S[C]
                products:   P[C]
                oligomeric_state: (1, 2)
            end))
        # Two values after the label, not one tuple.
        @test_throws "`oligomeric_state:` takes a single Int" eval(Meta.parse("""
            @enzyme_reaction begin
                substrates: S[C]
                products:   P[C]
                oligomeric_state: 1, 2
            end"""))
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
        forms = _testhelper_enzyme_forms(_testhelper_flat_steps(m))
        @test :EP in forms        # product-bound form (catalytic release)
        @test :EPinh in forms     # inhibitor-bound form — DISTINCT from :EP
        @test :EP != :EPinh
        # rate_equation must reference concs.P for the inhibitor term, not concs.Pinh:
        re_str = rate_equation_string(m)
        @test occursin("P", re_str)             # real metabolite name appears
        @test occursin("(; S, P) = concs", re_str)  # concentration is P, never Pinh
        @test !occursin("Pinh = concs", re_str)     # the inh marker never names a concentration
    end

    @testset "dual-role names bind by their catalytic role unless tagged ::Inh" begin
        _testhelper_bound_types(m) = Dict(
            EnzymeRates.name(EnzymeRates.to_species(s)) =>
                typeof(EnzymeRates.bound_metabolite(s))
            for g in EnzymeRates.steps(m) for s in g if EnzymeRates.is_binding(s))
        # A substrate also declared as a competitive inhibitor binds as the
        # substrate in a bare catalytic step; `E(A::Inh)` writes its inhibitor
        # form.
        m = EnzymeRates.Mechanism(@enzyme_mechanism begin
            substrates: A
            products:   P
            regulators: A
            steps: begin
                E + A ⇌ E(A)
                E(A) <--> E(P)
                E(P) ⇌ E + P
                E + A::Inh ⇌ E(A::Inh)
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
            allosteric_regulators: A::OnlyI, P::OnlyA
            catalytic_multiplicity: 2
            catalytic_steps: begin
                E + A ⇌ E(A)          :: EqualAI
                E(A) <--> E(P)        :: NonequalAI
                E(P) ⇌ E + P          :: EqualAI
            end
            regulatory_site(multiplicity = 2): begin
                ligands: A, P
            end
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
            allosteric_regulators: R::OnlyA
            catalytic_steps: begin
                E + S ⇌ E(S)      :: EqualAI
                E(S) <--> E(P)    :: EqualAI
                E(P) ⇌ E + P      :: EqualAI
                E + R ⇌ E(R)      :: EqualAI
            end
        end))
        @test_throws "`R` is an allosteric regulator" eval(:(@allosteric_mechanism begin
            substrates: S
            products:   P
            allosteric_regulators: R::OnlyA
            catalytic_steps: begin
                E + S ⇌ E(S)        :: EqualAI
                E(S) <--> E(P, R)   :: EqualAI
                E(P, R) ⇌ E + P     :: EqualAI
            end
        end))
        @test_throws "got `R`" eval(:(@allosteric_mechanism begin
            substrates: S
            products:   P
            allosteric_regulators: R::OnlyA
            catalytic_steps: begin
                E + S ⇌ E(S)                       :: EqualAI
                E(S) <--> Estar(; residual = S - R) :: EqualAI
                Estar(; residual = S - R) ⇌ E + P   :: EqualAI
            end
        end))
    end

    @testset "several metabolites on a step side" begin
        # A side with two enzyme forms is still rejected.
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
        # the two; `S - I` subtracts a competitive inhibitor.
        for residual in (:(P - S), :P, :(-S), :(S - I))
            @test_throws "adds substrates and subtracts products" eval(
                :(@enzyme_mechanism begin
                    substrates: S
                    products:   P
                    regulators: I
                    steps: begin
                        E + S ⇌ E(S)
                        E(S) <--> Estar(; residual = $residual)
                        Estar(; residual = $residual) ⇌ E + P
                    end
                end))
        end
    end

    @testset "::Inh is never a step tag" begin
        # `E(S::Inh) ⇌ E + S::Inh :: EqualAI` parses as `(S::Inh)::EqualAI`: the state
        # after `::Inh` is the step's tag, and the step reads like its binding direction.
        bind_inh = EnzymeRates.AllostericMechanism(@allosteric_mechanism begin
            substrates: S
            products:   P
            catalytic_inhibitors: S
            catalytic_steps: begin
                E + S ⇌ E(S)             :: EqualAI
                E(S) <--> E(P)           :: EqualAI
                E(P) ⇌ E + P             :: EqualAI
                E + S::Inh ⇌ E(S::Inh)   :: EqualAI
            end
        end)
        release_inh = EnzymeRates.AllostericMechanism(@allosteric_mechanism begin
            substrates: S
            products:   P
            catalytic_inhibitors: S
            catalytic_steps: begin
                E + S ⇌ E(S)             :: EqualAI
                E(S) <--> E(P)           :: EqualAI
                E(P) ⇌ E + P             :: EqualAI
                E(S::Inh) ⇌ E + S::Inh   :: EqualAI
            end
        end)
        @test release_inh == bind_inh

        # Without a state the step is missing its annotation.
        @test_throws "is missing" eval(:(@allosteric_mechanism begin
            substrates: S
            products:   P
            catalytic_inhibitors: I
            catalytic_steps: begin
                E + S ⇌ E(S)           :: EqualAI
                E(S) <--> E(P)         :: EqualAI
                E(P) ⇌ E + P           :: EqualAI
                E(I::Inh) ⇌ E + I::Inh
            end
        end))

        # A plain mechanism may release an inhibitor as well as bind it.
        bind_plain = EnzymeRates.Mechanism(@enzyme_mechanism begin
            substrates: S
            products:   P
            regulators: I
            steps: begin
                E + S ⇌ E(S)
                E(S) <--> E(P)
                E(P) ⇌ E + P
                E + I::Inh ⇌ E(I::Inh)
            end
        end)
        release_plain = EnzymeRates.Mechanism(@enzyme_mechanism begin
            substrates: S
            products:   P
            regulators: I
            steps: begin
                E + S ⇌ E(S)
                E(S) <--> E(P)
                E(P) ⇌ E + P
                E(I::Inh) ⇌ E + I::Inh
            end
        end)
        @test release_plain == bind_plain

        # A step tag that is not a name is rejected.
        @test_throws "must be a Symbol" eval(:(@allosteric_mechanism begin
            substrates: S
            products:   P
            catalytic_steps: begin
                E + S ⇌ E(S)    :: EqualAI
                E(S) <--> E(P)  :: EqualAI
                E(P) ⇌ E + P    :: 3
            end
        end))
    end
end
