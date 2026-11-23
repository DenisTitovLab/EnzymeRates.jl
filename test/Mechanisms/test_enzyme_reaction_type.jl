# ABOUTME: Tests of src/Mechanisms/enzyme_reaction_type.jl: EnzymeReaction, Reactant,
# ABOUTME: Regulator, the names a metabolite may take and a reaction's text.

# ── The reaction's types and their checks ─────────────────────────────────────────────

@testset "reaction fields are its declaration" begin
    rxn = @enzyme_reaction begin
        substrates: B[N], A[C2H]
        products: P[C2HN]
        competitive_inhibitors: I
        allosteric_regulators: X::Activator, Y::Inhibitor, Z
        oligomeric_state: (2, 1)
        shared_catalytic_site: (A, P)
    end
    @test [s.name for s in rxn.substrates] == [:A, :B]
    @test rxn.substrates[1].atoms == [:C => 2, :H => 1]
    @test rxn.products == [ER.Mechanisms.Reactant(:P, [:N => 1, :H => 1, :C => 2])]
    @test rxn.regulators == [ER.Mechanisms.Regulator(:I, :competitive, :unspecified),
                             ER.Mechanisms.Regulator(:X, :allosteric, :activator),
                             ER.Mechanisms.Regulator(:Y, :allosteric, :inhibitor),
                             ER.Mechanisms.Regulator(:Z, :allosteric, :unspecified)]
    @test rxn.allowed_oligomeric_states == [1, 2]
    @test rxn.shared_catalytic_site == [(:A, :P)]
    @test sprint(show, rxn) ==
          "@enzyme_reaction begin\n    substrates: A[C2H], B[N]\n" *
          "    products: P[C2HN]\n    competitive_inhibitors: I\n" *
          "    allosteric_regulators: X::Activator, Y::Inhibitor, Z\n" *
          "    oligomeric_state: (1, 2)\n    shared_catalytic_site: (A, P)\nend"
end

@testset "Reactant and Regulator compare by content" begin
    @test ER.Mechanisms.Reactant(:X, [:C => 1]) == ER.Mechanisms.Reactant(:X, [:C => 1])
    @test ER.Mechanisms.Reactant(:X, [:C => 1]) != ER.Mechanisms.Reactant(:Y, [:C => 1])
    @test hash(ER.Mechanisms.Reactant(:X, [:C => 1])) ==
          hash(ER.Mechanisms.Reactant(:X, [:C => 1]))
    a = ER.Mechanisms.Regulator(:A, :allosteric)
    @test a == ER.Mechanisms.Regulator(:A, :allosteric, :unspecified)
    @test hash(a) == hash(ER.Mechanisms.Regulator(:A, :allosteric))
    @test a != ER.Mechanisms.Regulator(:A, :competitive)
    @test a != ER.Mechanisms.Regulator(:A, :allosteric, :activator)
end

@testset "Reactant canonicalizes atom ordering" begin
    ra1 = ER.Mechanisms.Reactant(:ATP, [:C => 10, :H => 16, :N => 5])
    ra2 = ER.Mechanisms.Reactant(:ATP, [:N => 5, :H => 16, :C => 10])
    @test ra1 == ra2
    @test hash(ra1) == hash(ra2)
    @test ra1.name === :ATP
    @test ra1.atoms == [:C => 10, :H => 16, :N => 5]

    # An element listed twice counts once with the summed count, so every writing of
    # one formula gives one value.
    ra3 = ER.Mechanisms.Reactant(:ATP, [:C => 4, :N => 5, :C => 6, :H => 16])
    @test ra3 == ra1 && hash(ra3) == hash(ra1)
    @test ra3.atoms == [:C => 10, :H => 16, :N => 5]
end

@testset "Reactant validation" begin
    # Mandatory atoms: empty atom list is rejected.
    @test_throws ErrorException ER.Mechanisms.Reactant(:S, Pair{Symbol,Int}[])
    # Positive counts: zero / negative rejected.
    @test_throws ErrorException ER.Mechanisms.Reactant(:S, [:C => 0])
    @test_throws ErrorException ER.Mechanisms.Reactant(:S, [:C => -1])
    # Bool is not a valid count. Use untyped vector so Bool survives to
    # the constructor (Pair{Symbol,Int}[:C => true] would convert true→1).
    @test_throws ErrorException ER.Mechanisms.Reactant(:S, [:C => true])
    # An element is what a formula in `@enzyme_reaction` writes: a capital letter and
    # any lowercase letters.
    for el in (:c, :CL, Symbol("C1"))
        @test_throws(
            ErrorException("Reactant: S has the element $el; an " *
                           "element is a capital letter and any lowercase " *
                           "letters, such as C or Cl"),
            ER.Mechanisms.Reactant(:S, [el => 1]))
    end
    # Valid construction still works.
    @test ER.Mechanisms.Reactant(:S, [:C => 6, :H => 12]).atoms == [:C => 6, :H => 12]
end

@testset "EnzymeReaction validation" begin
    M = ER.Mechanisms
    S = M.Reactant(:S, [:C => 1])
    P = M.Reactant(:P, [:C => 1])
    Punbal = M.Reactant(:P, [:C => 2])
    noregs = M.Regulator[]
    regS = M.Regulator(:S, :competitive)
    regA = M.Regulator(:A, :competitive)
    regA2 = M.Regulator(:A, :competitive)

    # Empty substrate set rejected.
    @test_throws ErrorException ER.EnzymeReaction(M.Reactant[], [P], noregs, Int[1])
    # Empty product set rejected.
    @test_throws ErrorException ER.EnzymeReaction([S], M.Reactant[], noregs, Int[1])
    # Duplicate product names rejected.
    @test_throws ErrorException ER.EnzymeReaction([S], [P, P], noregs, Int[1])
    # Atom imbalance rejected (S has 1 C, Punbal has 2 C).
    @test_throws ErrorException ER.EnzymeReaction([S], [Punbal], noregs, Int[1])
    # Two regulators of the same kind with the same name rejected.
    @test_throws ErrorException ER.EnzymeReaction([S], [P], [regA, regA2], Int[1])
    # A regulator MAY share a substrate/product name (`::Mod` role tag:
    # one metabolite binds as a competitive inhibitor under its real name).
    rxn_inh = ER.EnzymeReaction([S], [P], [regS], Int[1])
    @test :S in [s.name for s in rxn_inh.substrates]
    # Valid balanced reaction with a distinct-named regulator still constructs.
    rxn = ER.EnzymeReaction([S], [P], [regA], Int[1])
    @test [s.name for s in rxn.substrates] == [:S]
    @test [p.name for p in rxn.products] == [:P]
    # A name listed as both a substrate and a product is rejected, whether it
    # is the only reactant name (X → X) or one of several (A + B → A + C).
    Xs = M.Reactant(:X, [:C => 2])
    @test_throws ErrorException(
        "EnzymeReaction: X named as both a substrate and a product; " *
        "concentrations and constants are keyed by name") ER.EnzymeReaction(
        [Xs], [Xs], noregs, Int[1])
    As = M.Reactant(:A, [:C => 1])
    Bs = M.Reactant(:B, [:N => 1])
    Cp = M.Reactant(:C, [:N => 1])
    @test_throws ErrorException(
        "EnzymeReaction: A named as both a substrate and a product; " *
        "concentrations and constants are keyed by name") ER.EnzymeReaction(
        [As, Bs], [As, Cp], noregs, Int[1])
    # A reactant or a regulator takes no name a rate equation keeps for its constants
    # and the local names of its code, and no name holding `_`, which constant names
    # use to join names.
    _testhelper_reserved(n) = ErrorException(
        "EnzymeReaction: $n is a reserved name; Keq, L, N, params, concs " *
        "and v are kept for a rate equation's constants and local names")
    for n in (:Keq, :L, :N, :params, :concs, :v)
        @test_throws _testhelper_reserved(n) ER.EnzymeReaction(
            [M.Reactant(n, [:C => 1])], [P], noregs, Int[1])
        @test_throws _testhelper_reserved(n) ER.EnzymeReaction(
            [S], [P], [M.Regulator(n, :allosteric)], Int[1])
    end
    @test_throws _testhelper_reserved(:Keq) @enzyme_reaction begin
        substrates: S[C]
        products: P[C]
        allosteric_regulators: Keq
    end
    _testhelper_underscore(n) = ErrorException(
        "EnzymeReaction: the metabolite name $n holds `_`; constant names join " *
        "metabolite names with `_`, so a metabolite name holds none")
    @test_throws _testhelper_underscore(:L_Ala) ER.EnzymeReaction(
        [M.Reactant(:L_Ala, [:C => 1])], [P], noregs, Int[1])
    @test_throws _testhelper_underscore(:I_X) @enzyme_reaction begin
        substrates: S[C]
        products: P[C]
        competitive_inhibitors: I_X
    end
    # Nor is a metabolite named E, the free enzyme form's label, which would make the
    # enumerator's forms read as the metabolite; nor does it take a name that no DSL
    # text can write.
    free_E = ErrorException("EnzymeReaction: E cannot name a metabolite; it labels " *
                            "the free enzyme form")
    @test_throws free_E ER.EnzymeReaction(
        [M.Reactant(:E, [:C => 1])], [P], noregs, Int[1])
    @test_throws free_E @enzyme_reaction begin
        substrates: S[C]
        products: P[C]
        allosteric_regulators: E
    end
    for n in (Symbol("1A"), Symbol("A B"), :end, :in)
        @test_throws(ErrorException("EnzymeReaction: \"$n\" cannot name a " *
                                    "metabolite; a metabolite name is an identifier " *
                                    "that is no keyword"),
                     ER.EnzymeReaction([S], [P], [M.Regulator(n, :competitive)],
                                       Int[1]))
    end
    # A contextual keyword is a keyword only in its own construct (`abstract type`,
    # `for outer`, `var"..."`), which no DSL text writes, so it names a metabolite and
    # the reaction's text reads back.
    for n in (:abstract, :mutable, :outer, :primitive, :public, :type, :var)
        r = ER.EnzymeReaction([M.Reactant(n, [:C => 1])], [P], noregs, Int[1])
        @test eval(Meta.parse(sprint(show, r))) == r
    end
    # Julia reads some identifiers as others: a micro sign as μ, a decomposed é as one
    # character. Such a name prints as text that builds another reaction, so it is
    # rejected; a name already in the form Julia reads, Greek or with a subscript,
    # reads back. The rejected strings are written with escapes: Julia would read the
    # characters themselves as the other names.
    for (n, read) in (("\u00b5M", "\u03bcM"), ("e\u0301", "\u00e9"))
        @test_throws(ErrorException("EnzymeReaction: Julia reads the name \"$n\" as " *
                                    "\"$read\", so no DSL text can write it; use " *
                                    "\"$read\""),
                     ER.EnzymeReaction([M.Reactant(Symbol(n), [:C => 1])], [P],
                                       noregs, Int[1]))
    end
    for n in (:α, :Δ, :μM, :S₁, :é)
        r = ER.EnzymeReaction([M.Reactant(n, [:C => 1])], [P], noregs, Int[1])
        @test eval(Meta.parse(sprint(show, r))) == r
    end
    # A single name MAY be declared in BOTH regulator roles: one allosteric
    # regulator and one competitive inhibitor. The two roles render to distinct
    # parameter names, so both are kept.
    regAllo = M.Regulator(:A, :allosteric)
    rxn_dual = ER.EnzymeReaction([S], [P], [regA, regAllo], Int[1])
    @test length(rxn_dual.regulators) == 2
    @test Set(reg.kind for reg in rxn_dual.regulators) ==
          Set([:allosteric, :competitive])
    @test all(reg.name == :A for reg in rxn_dual.regulators)
    # Two allosteric regulators of the same name still rejected.
    @test_throws ErrorException ER.EnzymeReaction(
        [S], [P], [regAllo, M.Regulator(:A, :allosteric)], Int[1])
end

@testset "EnzymeReaction shared_catalytic_site" begin
    M = ER.Mechanisms
    A = M.Reactant(:A, [:C => 1])
    B = M.Reactant(:B, [:C => 1])
    P = M.Reactant(:P, [:C => 1])
    Q = M.Reactant(:Q, [:C => 1])
    noregs = M.Regulator[]

    # Accepts a valid (substrate, product) pair; stored normalized + sorted.
    r = ER.EnzymeReaction([A, B], [P, Q], noregs, Int[1];
        shared_catalytic_site = [(:A, :P)])
    @test r.shared_catalytic_site == [(:A, :P)]

    # Normalizes product-first input to (substrate, product).
    r2 = ER.EnzymeReaction([A, B], [P, Q], noregs, Int[1];
        shared_catalytic_site = [(:P, :A)])
    @test r2.shared_catalytic_site == [(:A, :P)]

    # Sorts multiple pairs.
    r3 = ER.EnzymeReaction([A, B], [P, Q], noregs, Int[1];
        shared_catalytic_site = [(:B, :Q), (:A, :P)])
    @test r3.shared_catalytic_site == [(:A, :P), (:B, :Q)]

    # Default is empty.
    @test ER.EnzymeReaction([A], [P], noregs, Int[1]).shared_catalytic_site ==
        Tuple{Symbol,Symbol}[]

    # Rejections.
    @test_throws ErrorException ER.EnzymeReaction(  # unknown name
        [A, B], [P, Q], noregs, Int[1]; shared_catalytic_site = [(:A, :Z)])
    @test_throws ErrorException ER.EnzymeReaction(  # two substrates
        [A, B], [P, Q], noregs, Int[1]; shared_catalytic_site = [(:A, :B)])
    @test_throws ErrorException ER.EnzymeReaction(  # two products
        [A, B], [P, Q], noregs, Int[1]; shared_catalytic_site = [(:P, :Q)])
    @test_throws ErrorException ER.EnzymeReaction(  # duplicate pair
        [A, B], [P, Q], noregs, Int[1]; shared_catalytic_site = [(:A, :P), (:P, :A)])

    # A metabolite that is BOTH a substrate and a competitive inhibitor may
    # still be the substrate side of a shared pair (validation keys on the
    # substrate/product role, not the inhibitor role).
    regA = M.Regulator(:A, :competitive)
    r_dual = ER.EnzymeReaction([A, B], [P, Q], [regA], Int[1];
        shared_catalytic_site = [(:A, :P)])
    @test r_dual.shared_catalytic_site == [(:A, :P)]

    # == / hash reflect the field.
    @test r != ER.EnzymeReaction([A, B], [P, Q], noregs, Int[1])
    @test hash(r) != hash(ER.EnzymeReaction([A, B], [P, Q], noregs, Int[1]))
    @test r == ER.EnzymeReaction([A, B], [P, Q], noregs, Int[1];
        shared_catalytic_site = [(:A, :P)])
end

@testset "Regulator kind and effect" begin
    M = ER.Mechanisms
    @test M.Regulator(:X, :allosteric).effect == :unspecified      # default
    @test M.Regulator(:X, :allosteric, :activator).effect == :activator
    @test M.Regulator(:X, :allosteric) != M.Regulator(:X, :allosteric, :activator)
    @test_throws(ErrorException("Regulator: kind must be :allosteric or " *
                                ":competitive, got bogus"),
                 M.Regulator(:X, :bogus))
    @test_throws(ErrorException("Regulator: effect must be :activator, :inhibitor, " *
                                "or :unspecified, got bogus"),
                 M.Regulator(:X, :allosteric, :bogus))
    # A competitive inhibitor takes no type: only an allosteric regulator activates or
    # inhibits through the conformation it favors.
    for t in (:activator, :inhibitor)
        @test_throws(
            ErrorException("Regulator: X is a competitive " *
                           "inhibitor, which takes no regulator type, got " *
                           ":$t; only an allosteric regulator is an " *
                           "activator (::Activator) or an inhibitor " *
                           "(::Inhibitor)"),
            M.Regulator(:X, :competitive, t))
    end
end

@testset "EnzymeReaction canonicalizes reactant + regulator ordering" begin
    M = ER.Mechanisms
    r1 = EnzymeReaction(
        [M.Reactant(:B, [:C => 1]), M.Reactant(:A, [:C => 1])],
        [M.Reactant(:P, [:C => 2])],
        [M.Regulator(:Y, :allosteric), M.Regulator(:X, :allosteric)],
        [3, 1, 2],
    )
    r2 = EnzymeReaction(
        [M.Reactant(:A, [:C => 1]), M.Reactant(:B, [:C => 1])],
        [M.Reactant(:P, [:C => 2])],
        [M.Regulator(:X, :allosteric), M.Regulator(:Y, :allosteric)],
        [1, 2, 3],
    )
    @test r1 == r2
    @test hash(r1) == hash(r2)
end

@testset "EnzymeReaction takes one or more oligomeric states, each ≥ 1" begin
    _testhelper_reaction(mults) =
        EnzymeReaction([ER.Mechanisms.Reactant(:S, [:C => 1])],
                       [ER.Mechanisms.Reactant(:P, [:C => 1])],
                       ER.Mechanisms.Regulator[], mults)
    @test_throws(ErrorException("EnzymeReaction: the allowed oligomeric states " *
                                "must be a nonempty list of Ints ≥ 1, got [0, 1]"),
                 _testhelper_reaction([0, 1]))
    # An empty list would leave the search no multiplicity to assign.
    @test_throws(ErrorException("EnzymeReaction: the allowed oligomeric states " *
                                "must be a nonempty list of Ints ≥ 1, got Int64[]"),
                 _testhelper_reaction(Int[]))
    # A multiplicity listed twice is no list a DSL text writes.
    @test_throws(ErrorException("EnzymeReaction: the allowed oligomeric states list " *
                                "a multiplicity more than once; got [2, 2]."),
                 _testhelper_reaction([2, 2]))
    # A Bool is no multiplicity: the list is a Vector{Int}, and none is converted.
    @test_throws MethodError _testhelper_reaction(Bool[true])
end

@testset "a regulatory-site ligand takes no reserved name and no `_`" begin
    # An allosteric regulator reaches the constructor only as a site ligand, so the
    # Mechanism constructor checks the ligands' names.
    @test_throws(ErrorException(
        "Mechanism: L is a reserved name; Keq, L, N, params, concs " *
        "and v are kept for a rate equation's constants and local names"),
        @enzyme_mechanism begin
            substrates: S
            products:   P
            steps: begin
                E + S ⇌ ES      :: EqualAI
                ES <--> EP    :: OnlyA
                EP ⇌ E + P      :: EqualAI
            end
            regulatory_site(2): L::OnlyI
        end)
    @test_throws(ErrorException(
        "Mechanism: the metabolite name A_R holds `_`; constant names join " *
        "metabolite names with `_`, so a metabolite name holds none"),
        @enzyme_mechanism begin
            substrates: S
            products:   P
            steps: begin
                E + S ⇌ ES      :: EqualAI
                ES <--> EP    :: OnlyA
                EP ⇌ E + P      :: EqualAI
            end
            regulatory_site(2): A_R::EqualAI
        end)
end

# ── The names a metabolite may not take ───────────────────────────────────────────────

# The rate equation binds its constants and the concentrations as local names of one
# scope, so no metabolite takes a constant's name: a reserved name, or one holding `_`,
# is rejected when a mechanism macro's code runs.
@testset "a metabolite named after a constant is rejected" begin
    # An allosteric regulator named L would take the conformational constant's name.
    @test_throws(ErrorException(
        "Mechanism: L is a reserved name; Keq, L, N, params, concs " *
        "and v are kept for a rate equation's constants and local names"),
        @enzyme_mechanism begin
            substrates: S
            products: P
            oligomeric_state: 2
            steps: begin
                E + S ⇌ ES     :: EqualAI
                ES <--> EP   :: OnlyA
                EP ⇌ E + P     :: EqualAI
            end
            regulatory_site(2): L::OnlyI
        end)
    # An inhibitor named after the reverse chemistry constant holds `_`.
    @test_throws(ErrorException(
        "Mechanism: the metabolite name k_EP_to_ES holds `_`; constant names join " *
        "metabolite names with `_`, so a metabolite name holds none"),
        @enzyme_mechanism begin
            substrates: S
            products: P
            modifiers: k_EP_to_ES
            steps: begin
                E + S ⇌ ES
                ES <--> EP
                EP ⇌ E + P
                E + k_EP_to_ES ⇌ Ek_EP_to_ES
            end
        end)
end

# Every constant of a rate equation needs a name of its own. An :EqualAI ligand's
# constant is bare, K_<ligand>reg, while an :OnlyA or :NonequalAI ligand's names its
# conformation, K_A_<ligand>reg or K_I_<ligand>reg, so a ligand named A_R or I_X would
# take the constant of a ligand R or X. A ligand name holds no `_`.
@testset "a ligand name holding `_` is rejected" begin
    @test_throws(ErrorException(
        "Mechanism: the metabolite name A_R holds `_`; constant names join " *
        "metabolite names with `_`, so a metabolite name holds none"),
        @enzyme_mechanism begin
            substrates: S
            products: P
            oligomeric_state: 2
            steps: begin
                E + S ⇌ ES      :: EqualAI
                ES <--> EP    :: NonequalAI
                EP ⇌ E + P      :: EqualAI
            end
            regulatory_site(2): R::OnlyA
            regulatory_site(2): A_R::EqualAI
        end)
    @test_throws(ErrorException(
        "Mechanism: the metabolite name I_X holds `_`; constant names join " *
        "metabolite names with `_`, so a metabolite name holds none"),
        @enzyme_mechanism begin
            substrates: S
            products: P
            oligomeric_state: 2
            steps: begin
                E + S ⇌ ES      :: EqualAI
                ES <--> EP    :: NonequalAI
                EP ⇌ E + P      :: EqualAI
            end
            regulatory_site(2): X::NonequalAI, I_X::EqualAI
        end)

    # Renamed, the :EqualAI ligand has a constant of its own.
    renamed = @enzyme_mechanism begin
        substrates: S
        products: P
        oligomeric_state: 2
        steps: begin
            E + S ⇌ ES      :: EqualAI
            ES <--> EP    :: NonequalAI
            EP ⇌ E + P      :: EqualAI
        end
        regulatory_site(2): R::OnlyA
        regulatory_site(2): T::EqualAI
    end
    @test count(==(:K_A_Rreg), parameters(renamed)) == 1
    @test :K_Treg in parameters(renamed)
end

# ── A reaction's text ─────────────────────────────────────────────────────────────────

@testset "a reaction prints as the @enzyme_reaction text that builds it" begin
    # The reaction shapes the tests use: one formula per reactant, competitive
    # inhibitors with a reactant among them, allosteric regulators with and without
    # types, a name in both regulator roles, declared in either order, one
    # or several oligomeric states and shared catalytic sites. The printed
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
            allosteric_regulators: X::Activator, Y::Inhibitor, Z, W
            oligomeric_state: 4
        end),
        @enzyme_reaction(begin
            substrates: A[C2H3], B[NP]
            products:   P[C2N], Q[H3P]
            competitive_inhibitors: ATP, I
            allosteric_regulators: ATP
            oligomeric_state: (4, 1, 2)
            shared_catalytic_site: (A, P), (Q, B)
        end),
        @enzyme_reaction(begin
            substrates: S[C]
            products:   P[C]
            allosteric_regulators: ATP
            competitive_inhibitors: ATP
            oligomeric_state: 2
        end)]
    for rxn in reactions
        @test eval(Meta.parse(sprint(show, rxn))) == rxn
    end
    @test sprint(show, reactions[2]) == """
        @enzyme_reaction begin
            substrates: ATP[C10H16N5O13P3], S[C6H12O6]
            products: ADP[C10H15N5O10P2], P[C6H13O9P]
            competitive_inhibitors: I, P
            allosteric_regulators: W, X::Activator, Y::Inhibitor, Z
            oligomeric_state: 4
        end"""
    @test sprint(show, reactions[3]) == """
        @enzyme_reaction begin
            substrates: A[C2H3], B[NP]
            products: P[C2N], Q[H3P]
            competitive_inhibitors: ATP, I
            allosteric_regulators: ATP
            oligomeric_state: (1, 2, 4)
            shared_catalytic_site: (A, P), (B, Q)
        end"""
end
