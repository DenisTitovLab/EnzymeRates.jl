# ABOUTME: Tests of src/DSL/enzyme_reaction_macro.jl: the reactions @enzyme_reaction
# ABOUTME: builds and the blocks it rejects.

@testset "@enzyme_reaction" begin
    spec = @enzyme_reaction begin
        substrates: S[C]
        products:   P[C]
    end
    @test spec isa EnzymeReaction
    @test spec.substrates == [EnzymeRates.Mechanisms.Reactant(:S, [:C => 1])]
    @test spec.products == [EnzymeRates.Mechanisms.Reactant(:P, [:C => 1])]
    @test spec.regulators == EnzymeRates.Mechanisms.Regulator[]
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
    @test [s.name for s in rxn.substrates] == [:A, :B]
    ra_map = Dict(r.name => r.atoms for r in [rxn.substrates; rxn.products])
    @test ra_map[:A] == [:C => 2, :H => 3]
    @test ra_map[:B] == [:N => 1, :P => 1]
    @test [p.name for p in rxn.products] == [:P, :Q]
    @test ra_map[:P] == [:C => 2, :N => 1]
    @test ra_map[:Q] == [:H => 3, :P => 1]

    # A formula that repeats an element means the summed count: S[CHC] is S[C2H], and
    # both print as S[C2H].
    repeated = @enzyme_reaction begin
        substrates: S[CHC]
        products: P[C2H]
    end
    @test repeated == @enzyme_reaction(begin
        substrates: S[C2H]
        products: P[C2H]
    end)
    @test sprint(show, repeated) == """
        @enzyme_reaction begin
            substrates: S[C2H]
            products: P[C2H]
        end"""

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
    @test_throws "`oligomeric_state:` entry must be a positive Int, got 0." eval(
        :(@enzyme_reaction begin
            substrates: S[C]
            products:   P[C]
            oligomeric_state: 0
        end))
    @test_throws "`oligomeric_state:` must be a positive Int, got 2.0." eval(
        :(@enzyme_mechanism begin
            substrates: S
            products:   P
            oligomeric_state: 2.0
            steps: begin
                E + S ⇌ ES      :: EqualAI
                ES <--> EP    :: EqualAI
                EP ⇌ E + P      :: EqualAI
            end
        end))
    @test_throws "`regulatory_site` multiplicity must be a positive Int, got 0." eval(
        :(@enzyme_mechanism begin
            substrates: S
            products:   P
            steps: begin
                E + S ⇌ ES      :: EqualAI
                ES <--> EP    :: EqualAI
                EP ⇌ E + P      :: EqualAI
            end
            regulatory_site(0): A::OnlyA
        end))
    # `true` is no multiplicity, though Julia counts a Bool an integer.
    @test_throws "`oligomeric_state:` entry must be a positive Int, got true." eval(
        :(@enzyme_reaction begin
            substrates: S[C]
            products:   P[C]
            oligomeric_state: true
        end))
    @test_throws "`oligomeric_state:` must be a positive Int, got true." eval(
        :(@enzyme_mechanism begin
            substrates: S
            products:   P
            oligomeric_state: true
            steps: begin
                E + S ⇌ ES      :: EqualAI
                ES <--> EP    :: EqualAI
                EP ⇌ E + P      :: EqualAI
            end
        end))
    @test_throws("`regulatory_site` multiplicity must be a positive Int, got true.",
                 eval(:(@enzyme_mechanism begin
            substrates: S
            products:   P
            steps: begin
                E + S ⇌ ES      :: EqualAI
                ES <--> EP    :: EqualAI
                EP ⇌ E + P      :: EqualAI
            end
            regulatory_site(true): A::OnlyA
        end)))
end

@testset "a regulator multiplicity list is a parse error" begin
    # A regulator takes the oligomeric states that `oligomeric_state:` lists; a list
    # after its name, tagged or not, is rejected and names the regulator.
    @test_throws("@enzyme_reaction: regulator I lists multiplicities, got `I(1, 2)`; " *
                 "a regulator takes the oligomeric states that `oligomeric_state:` " *
                 "lists, so list I bare.", eval(:(@enzyme_reaction begin
        substrates: S[C]
        products:   P[C]
        competitive_inhibitors: I(1, 2)
    end)))
    @test_throws("@enzyme_reaction: regulator A lists multiplicities, got " *
                 "`A(4)::Inhibitor`; a regulator takes the oligomeric states that " *
                 "`oligomeric_state:` lists, so list A bare.",
                 eval(:(@enzyme_reaction begin
        substrates: S[C]
        products:   P[C]
        allosteric_regulators: A(4)::Inhibitor
        oligomeric_state: 4
    end)))
end

@testset "@enzyme_reaction regulator kinds" begin
    # competitive_inhibitors: emits :competitive regulators, and its lines
    # accumulate. allosteric_regulators: emits :allosteric regulators.
    spec_kinds = @enzyme_reaction begin
        substrates: S[C]
        products: P[C]
        competitive_inhibitors: I
        allosteric_regulators: A
        competitive_inhibitors: R
    end
    @test spec_kinds isa EnzymeReaction
    reg_by_name = Dict(reg.name => reg for reg in spec_kinds.regulators)
    @test Set(keys(reg_by_name)) == Set([:I, :A, :R])
    @test reg_by_name[:I] == EnzymeRates.Mechanisms.Regulator(:I, :competitive)
    @test reg_by_name[:A] == EnzymeRates.Mechanisms.Regulator(:A, :allosteric)
    @test reg_by_name[:R] == EnzymeRates.Mechanisms.Regulator(:R, :competitive)
end

@testset "@enzyme_reaction dual-role regulator (allosteric + competitive)" begin
    # One metabolite MAY be declared as BOTH an allosteric regulator and a
    # competitive inhibitor; the two roles render to distinct parameter
    # names, so both Regulators are carried.
    rxn = @enzyme_reaction begin
        substrates: S[C]
        products: P[C]
        competitive_inhibitors: ATP
        allosteric_regulators: ATP
    end
    @test rxn isa EnzymeReaction
    atp_regs = [reg for reg in rxn.regulators if reg.name == :ATP]
    @test length(atp_regs) == 2
    @test Set(reg.kind for reg in atp_regs) == Set([:allosteric, :competitive])
    # The order of the two labels does not matter: both orders build one reaction.
    allosteric_first = @enzyme_reaction begin
        substrates: S[C]
        products: P[C]
        allosteric_regulators: ATP
        competitive_inhibitors: ATP
    end
    @test allosteric_first == rxn
    @test hash(allosteric_first) == hash(rxn)
    # Two allosteric ATPs are rejected.
    @test_throws Exception eval(:(@enzyme_reaction begin
        substrates: S[C]
        products: P[C]
        allosteric_regulators: ATP, ATP
    end))
end

@testset "@enzyme_reaction regulator types" begin
    rxn = @enzyme_reaction begin
        substrates: A[C6H12O6]
        products: B[C6H12O6]
        allosteric_regulators: X::Activator, Y::Inhibitor, Z
        oligomeric_state: 2
    end
    reg_types = Dict(reg.name => reg.effect for reg in rxn.regulators)
    @test reg_types[:X] == :activator
    @test reg_types[:Y] == :inhibitor
    @test reg_types[:Z] == :unspecified

    @test_throws Exception eval(:(@enzyme_reaction begin
        substrates: A[C6H12O6]
        products: B[C6H12O6]
        allosteric_regulators: X::Bogus
        oligomeric_state: 2
    end))
    # A competitive inhibitor takes no type tag; the reaction's regulator rejects it.
    for (tag, t) in ((:Activator, :activator), (:Inhibitor, :inhibitor))
        @test_throws(ErrorException("Regulator: X is a competitive " *
                                    "inhibitor, which takes no regulator type, got " *
                                    ":$t; only an allosteric regulator is an " *
                                    "activator (::Activator) or an inhibitor " *
                                    "(::Inhibitor)"),
                     eval(:(@enzyme_reaction begin
            substrates: A[C6H12O6]
            products: B[C6H12O6]
            competitive_inhibitors: X::$tag
        end)))
    end
end

@testset "@enzyme_reaction names its labels for any other" begin
    # Regulators go under `competitive_inhibitors:` or `allosteric_regulators:`, the
    # oligomeric states under `oligomeric_state:`; any other label is reported
    # as unknown, with the valid ones.
    for label in (:regulators, :dead_end_inhibitors, :allowed_oligomeric_states)
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

@testset "@enzyme_reaction with oligomeric_state" begin
    rxn = @enzyme_reaction begin
        substrates: S[C]
        products: P[C]
        oligomeric_state: 4
    end
    @test rxn.allowed_oligomeric_states == [4]

    # Without oligomeric_state defaults to 1
    rxn2 = @enzyme_reaction begin
        substrates: S[C]
        products: P[C]
    end
    @test rxn2.allowed_oligomeric_states == [1]

    # A tuple of multiplicities
    rxn3 = @enzyme_reaction begin
        substrates: S[C]
        products: P[C]
        oligomeric_state: (1, 2, 4)
    end
    @test rxn3.allowed_oligomeric_states == [1, 2, 4]
end

@testset "@enzyme_reaction sets the oligomeric states once" begin
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
    @test one_tuple.allowed_oligomeric_states == [2]
    several = @enzyme_reaction begin
        substrates: S[C]
        products:   P[C]
        oligomeric_state: (4, 1, 2)
    end
    @test several.allowed_oligomeric_states == [1, 2, 4]
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
    # Each multiplicity once.
    @test_throws("EnzymeReaction: the allowed oligomeric states list a multiplicity " *
                 "more than once; got [2, 2].", eval(:(@enzyme_reaction begin
        substrates: S[C]
        products:   P[C]
        allosteric_regulators: A
        oligomeric_state: (2, 2)
    end)))
end

@testset "@enzyme_reaction needs a substrate and a product" begin
    @test_throws(ErrorException("EnzymeReaction: no substrate is given; a reaction " *
                                "has at least one substrate."),
                 @enzyme_reaction begin
                     products: P[C]
                 end)
    @test_throws(ErrorException("EnzymeReaction: no product is given; a reaction " *
                                "has at least one product."),
                 @enzyme_reaction begin
                     substrates: S[C]
                 end)
end

@testset "@enzyme_reaction shared_catalytic_site" begin
    rxn = @enzyme_reaction begin
        substrates: A[C], B[C]
        products:   P[C], Q[C]
        shared_catalytic_site: (A, P), (B, Q)
    end
    @test rxn.shared_catalytic_site == [(:A, :P), (:B, :Q)]

    # Single pair.
    rxn1 = @enzyme_reaction begin
        substrates: A[C], B[C]
        products:   P[C], Q[C]
        shared_catalytic_site: (P, A)
    end
    @test rxn1.shared_catalytic_site == [(:A, :P)]

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
