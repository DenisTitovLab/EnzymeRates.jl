# ABOUTME: Tests for the DSL module: the two macros, the walker's raw values and
# ABOUTME: parse_mechanism, which reads a mechanism back from its text.

using CSV
using DataFrames

@testset "DSL module entry points" begin
    @testset "the macros work when imported alone" begin
        importer = Module(:_testhelper_MacroOnly)
        Core.eval(importer, :(using EnzymeRates: @enzyme_mechanism, @enzyme_reaction))
        m = Core.eval(importer, :(@enzyme_mechanism begin
            substrates: S
            products: P
            steps: begin
                E + S ⇌ ES
                ES <--> EP
                EP ⇌ E + P
            end
        end))
        @test !_testhelper_is_allosteric(m)
        rxn = Core.eval(importer, :(@enzyme_reaction begin
            substrates: S[C]
            products: P[C]
            competitive_inhibitors: I
        end))
        @test rxn isa EnzymeRates.EnzymeReaction
    end

    @testset "parse_mechanism reads show back" begin
        _testhelper_text_reads_back(m) =
            sprint(show, EnzymeRates.DSL.parse_mechanism(sprint(show, m))) ==
            sprint(show, m)
        @testset "$(spec.name)" for spec in MECHANISM_TEST_SPECS
            @test _testhelper_text_reads_back(spec.mechanism)
        end
        bi_bi = @enzyme_reaction begin
            substrates: A[C], B[N]
            products: P[C], Q[N]
        end
        init = _testhelper_mechanisms(EnzymeRates.Enumeration.init_mechanisms(bi_bi))
        @test length(init) > 0
        @test all(_testhelper_text_reads_back, init)
    end
end

@testset "one macro" begin
    M = EnzymeRates.Mechanisms
    # The statements carry allosteric states, so the mechanism is allosteric. Its first
    # statement, the chemistry, is OnlyA.
    allosteric = @enzyme_mechanism begin
        substrates: S
        products: P
        oligomeric_state: 2
        steps: begin
            ES <--> EP      :: OnlyA
            E + S ⇌ ES        :: EqualAI
            EP ⇌ E + P        :: EqualAI
        end
        regulatory_site(2): X::OnlyA
    end
    none = Int[]
    # The steps in written order, each labelled by its statement, as the macro reads them,
    # between the forms E, ES and EP over the metabolites S and P.
    scheme = M.CatalyticScheme([:S, :P], [:E, :ES, :EP], [1], [2],
                               [M.Step(2, 3, none, none, false, 1),
                                M.Step(1, 2, [1], none, true, 2),
                                M.Step(3, 1, none, [2], true, 3)])
    sites = [M.RegulatorySite([:X], 2, [:OnlyA])]
    built = EnzymeRates.Mechanism(scheme,
                                  M.Allostery([:OnlyA, :EqualAI, :EqualAI], 2, sites))
    @test allosteric == built && hash(allosteric) == hash(built)
    plain = @enzyme_mechanism begin
        substrates: S
        products: P
        steps: begin
            ES <--> EP
            E + S ⇌ ES
            EP ⇌ E + P
        end
    end
    @test plain ==
          EnzymeRates.Mechanism(scheme, M.Allostery(Symbol[], 1, M.RegulatorySite[]))
    @test isempty(plain.allosteric.allosteric_states)
    @test plain.allosteric.oligomeric_state == 1
    # One catalytic scheme, plain and allosteric: two mechanisms.
    @test plain.catalytic == allosteric.catalytic
    @test plain != allosteric && hash(plain) != hash(allosteric)
    for m in (allosteric, plain)
        @test EnzymeRates.DSL.parse_mechanism(sprint(show, m)) == m
    end
    # `oligomeric_state:` and `regulatory_site(n):` belong to an allosteric mechanism, and
    # every statement carries a state or none does: each break is an ArgumentError, from
    # the macro and from the parser, that names the keyword or the statement.
    uniuni = """
        steps: begin
            E + S ⇌ ES
            ES <--> EP
            EP ⇌ E + P
        end
    """
    allosteric_only(keyword) =
        "@enzyme_mechanism: `$keyword` is read only for an allosteric mechanism, whose " *
        "statements each carry an allosteric state (`:: EqualAI`, `:: NonequalAI`, " *
        "`:: OnlyA` or `:: OnlyI`); no statement carries one."
    mixed = "@enzyme_mechanism: the statement `EP ⇌ E + P` carries no allosteric " *
            "state, while `E + S ⇌ ES::EqualAI` carries one; a mechanism is allosteric " *
            "when every statement ends in its state (`:: EqualAI`, `:: NonequalAI`, " *
            "`:: OnlyA` or `:: OnlyI`), and plain when none does."
    for (text, message) in (
            ("oligomeric_state: 2\n" * uniuni, allosteric_only("oligomeric_state:")),
            (uniuni * "regulatory_site(2): X::OnlyA\n",
             allosteric_only("regulatory_site(n):")),
            ("""
             steps: begin
                 E + S ⇌ ES        :: EqualAI
                 ES <--> EP      :: OnlyA
                 EP ⇌ E + P
             end
             """, mixed))
        block = "@enzyme_mechanism begin\nsubstrates: S\nproducts: P\n" * text * "end"
        @test_throws ArgumentError(message) macroexpand(@__MODULE__, Meta.parse(block))
        @test_throws ArgumentError(message) EnzymeRates.DSL.parse_mechanism(block)
    end
end

"""The step of mechanism `x` that binds the metabolite named `name`: it takes up that
metabolite and gives off nothing."""
_testhelper_binding_of(x, name) =
    only(e for e in x.catalytic.steps
         if isempty(e.given_off) && x.catalytic.metabolites[e.taken_up] == [name])

"""The names of the two forms step `e` of mechanism `x` joins, in its stored direction."""
_testhelper_ends(x, e) = x.catalytic.enzyme_forms[[e.from, e.to]]

@testset "DSL" begin
    @testset "DSL: a name on a step side is a declared metabolite or a form" begin
        # The parser reads each name on a step side that is no declared metabolite as an
        # enzyme form, used as written: in `E + S ⇌ ES`, E and ES are forms and S the
        # substrate.
        mech = @_testhelper_unguarded @enzyme_mechanism begin
            substrates: S
            products:   P
            steps: begin
                E + S ⇌ ES
                ES ⇌ EP
                EP ⇌ E + P
            end
        end
        # Every step at rapid equilibrium closes a fast cycle that turns the reaction
        # over (class 2); the parse is the subject, so the constructor builds it
        # unchecked.
        @test EnzymeRates.Derivation._guard(mech).fault ===
              EnzymeRates.Derivation._RE_TURNOVER
        @test length(_testhelper_group_members(mech.catalytic)) == 3
        @test mech.catalytic.metabolites == [:S, :P]
        @test mech.catalytic.enzyme_forms == [:E, :EP, :ES]
        # Steps are canonicalized; pick each by content, not position.
        es_step = _testhelper_binding_of(mech, :S)
        @test _testhelper_ends(mech, es_step) == [:E, :ES]

        # ES ⇌ EP, an RE isomerization. The Mechanism constructor stores an isomerization
        # the way the reaction turns over through it (`_orientation`), so ES is its `from`
        # form and EP its `to` form.
        iso_step = only(filter(EnzymeRates.Mechanisms._is_iso, mech.catalytic.steps))
        @test _testhelper_ends(mech, iso_step) == [:ES, :EP]

        # Substrates, products and modifiers by their declared roles.
        m_multi = @enzyme_mechanism begin
            substrates: A, B
            products:   P
            modifiers: I
            steps: begin
                E + A <--> EA
                EA + B <--> EAB
                EAB <--> EP
                EP <--> E + P
                E + I <--> EI
            end
        end
        @test m_multi.catalytic.metabolites == [:A, :B, :P, :I]
        @test (m_multi.catalytic.substrates, m_multi.catalytic.products) == ([1, 2], [3])
        @test m_multi.catalytic.enzyme_forms == [:E, :EA, :EAB, :EI, :EP]
        # The modifier I binds E.
        @test _testhelper_ends(m_multi, _testhelper_binding_of(m_multi, :I)) == [:E, :EI]
    end

    @testset "@enzyme_mechanism step grammar: groups, forms by name, modifiers" begin
        m = @enzyme_mechanism begin
            substrates: S
            products:   P
            modifiers: I

            steps: begin
                (E + S ⇌ ES, EP + S ⇌ EPS)
                ES + I ⇌ EIS
                ES <--> EP
                EP ⇌ E + P
            end
        end
        @test _testhelper_substrates(m) == [:S]
        @test _testhelper_products(m) == [:P]
        @test _testhelper_regulators(m) == [:I]
        # One shared 2-step kinetic group (the S-binding pair) + 3 singletons →
        # 5 steps, 4 groups (order-independent: canonicalization reorders steps).
        @test length(_testhelper_group_members(m.catalytic)) == 4

        # A covalent form, named for its residual: Estar_res_plusA_minusP.
        m_res = @enzyme_mechanism begin
            substrates: A, B
            products:   P, Q
            steps: begin
                E + A <--> EA
                EA <--> Estar_res_plusA_minusP
                Estar_res_plusA_minusP <--> EstarQ_res_plusA_minusP
                EstarQ_res_plusA_minusP <--> Estar_res_plusA_minusP + Q
                Estar_res_plusA_minusP + B <--> EstarB_res_plusA_minusP
                EstarB_res_plusA_minusP <--> E + P
            end
        end
        @test length(m_res.catalytic.steps) == 6

        # Reject atom bracket syntax in substrates:
        @test_throws "expects bare names" eval(:(@enzyme_mechanism begin
            substrates: S[C]
            products:   P
            steps: begin
                E + S ⇌ ES
                ES <--> EP
                EP ⇌ E + P
            end
        end))

        # Reject allosteric-only syntax (regulatory_site(...)) on a plain mechanism
        @test_throws "`regulatory_site(n):` is read only for an allosteric mechanism" eval(
            :(@enzyme_mechanism begin
                substrates: S
                products:   P
                regulatory_site(2): A::OnlyA
                steps: begin
                    E + S ⇌ ES
                    ES <--> EP
                    EP ⇌ E + P
                end
            end))
    end

    @testset "@enzyme_mechanism (parsing & validation)" begin
        # Reject untagged catalytic step
        @test_throws Exception eval(:(@enzyme_mechanism begin
            substrates: F6P
            products:   F16BP
            oligomeric_state: 2
            steps: begin
                E + F6P ⇌ EF6P :: EqualAI
                EF6P <--> EF16BP
                EF16BP ⇌ E + F16BP :: EqualAI
            end
        end))

        # Reject :OnlyI on the only catalytic step: the active conformation must catalyze.
        @test_throws ("Mechanism: without its :OnlyI steps EF6P <--> " *
                      "EF16BP the active conformation turns no substrate into " *
                      "product") eval(
            :(@enzyme_mechanism begin
                substrates: F6P
                products:   F16BP
                oligomeric_state: 2
                steps: begin
                    E + F6P ⇌ EF6P :: EqualAI
                    EF6P <--> EF16BP :: OnlyI
                    EF16BP ⇌ E + F16BP :: EqualAI
                end
            end))

        # Reject a parenthesized step group without a tag beside tagged statements.
        @test_throws("@enzyme_mechanism: the statement " *
                     "`(E + S ⇌ ES, EA + S ⇌ EAS)` carries no allosteric state, " *
                     "while `EAS <--> EP::EqualAI` carries one; a mechanism is " *
                     "allosteric when every statement ends in its state (`:: EqualAI`, " *
                     "`:: NonequalAI`, `:: OnlyA` or `:: OnlyI`), and plain when none " *
                     "does.", eval(:(@enzyme_mechanism begin
            substrates: S, A
            products:   P
            oligomeric_state: 2
            steps: begin
                (E + S ⇌ ES, EA + S ⇌ EAS)
                EAS <--> EP   :: EqualAI
                EP ⇌ E + P        :: EqualAI
            end
        end)))
    end

    @testset "@enzyme_mechanism checks names and sites; constructors check states" begin
        # An unknown regulator state is reported by RegulatorySite when the expansion
        # runs.
        @test_throws "allosteric state Bogus must be one of" eval(:(@enzyme_mechanism begin
            substrates: S
            products:   P
            steps: begin
                E + S ⇌ ES      :: EqualAI
                ES <--> EP    :: EqualAI
                EP ⇌ E + P      :: EqualAI
            end
            regulatory_site(1): X::Bogus
        end))

        # An unknown catalytic state is reported by the Mechanism constructor.
        for state in (:Bogus,)
            @test_throws "catalytic group" eval(:(@enzyme_mechanism begin
                substrates: S
                products:   P
                steps: begin
                    E + S ⇌ ES      :: EqualAI
                    ES <--> EP    :: EqualAI
                    EP ⇌ E + P     :: $state
                end
            end))
        end
        # The guard rejects a rapid-equilibrium :OnlyI binding whose two sides the active
        # conformation's other steps still join.
        joined = "Mechanism: in the active conformation the rapid-equilibrium " *
                 "step E + P ⇌ EP is tagged :OnlyI, but the steps that conformation " *
                 "keeps still join its two sides."
        @test_throws joined eval(:(@enzyme_mechanism begin
            substrates: S
            products:   P
            steps: begin
                E + S ⇌ ES      :: EqualAI
                ES <--> EP    :: EqualAI
                EP ⇌ E + P      :: OnlyI
            end
        end))

        # A ligand on two regulatory sites is reported by the macro, which reads the
        # site lines.
        @test_throws("@enzyme_mechanism: allosteric regulator `A` is listed 2 times " *
                     "on `regulatory_site(n):` lines; an allosteric regulator sits at " *
                     "exactly one regulatory site.", eval(
            :(@enzyme_mechanism begin
                substrates: S
                products:   P
                steps: begin
                    E + S ⇌ ES      :: EqualAI
                    ES <--> EP    :: EqualAI
                    EP ⇌ E + P      :: EqualAI
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
                E + S ⇌ ES
                ES <--> EP
                EP ⇌ E + P
            end
            steps: begin
                E + S ⇌ ES
                ES <--> EP
                EP ⇌ E + P
            end
        end))
        @test_throws "`oligomeric_state:` given more than once" eval(
            :(@enzyme_mechanism begin
                substrates: S
                products:   P
                oligomeric_state: 2
                oligomeric_state: 4
                steps: begin
                    E + S ⇌ ES      :: EqualAI
                    ES <--> EP    :: EqualAI
                    EP ⇌ E + P      :: EqualAI
                end
            end))
    end

    @testset "list labels may repeat; their names accumulate" begin
        # Each macro reads a list label written on two lines as one line listing both
        # lines' entries.
        @test @enzyme_mechanism(begin
            substrates: A
            products:   P
            modifiers: I
            substrates: B
            products:   Q
            modifiers: J
            steps: begin
                E + A ⇌ EA
                EA + B ⇌ EAB
                EAB <--> EPQ
                EPQ ⇌ EQ + P
                EQ ⇌ E + Q
                E + I ⇌ EI
                E + J ⇌ EJ
            end
        end) == @enzyme_mechanism(begin
            substrates: A, B
            products:   P, Q
            modifiers: I, J
            steps: begin
                E + A ⇌ EA
                EA + B ⇌ EAB
                EAB <--> EPQ
                EPQ ⇌ EQ + P
                EQ ⇌ E + Q
                E + I ⇌ EI
                E + J ⇌ EJ
            end
        end)
        @test @enzyme_mechanism(begin
            substrates: S
            products:   P
            oligomeric_state: 2
            steps: begin
                E + S ⇌ ES      :: EqualAI
                ES <--> EP    :: EqualAI
                EP ⇌ E + P      :: EqualAI
            end
            regulatory_site(2): X::OnlyA, Y::OnlyI
        end) == @enzyme_mechanism(begin
            substrates: S
            products:   P
            oligomeric_state: 2
            steps: begin
                E + S ⇌ ES      :: EqualAI
                ES <--> EP    :: EqualAI
                EP ⇌ E + P      :: EqualAI
            end
            regulatory_site(2): X::OnlyA, Y::OnlyI
        end)
        @test @enzyme_reaction(begin
            substrates: A[C]
            products:   P[C]
            allosteric_regulators: X
            shared_catalytic_site: (A, P)
            substrates: B[N]
            products:   Q[N]
            allosteric_regulators: Y::Inhibitor
            shared_catalytic_site: (B, Q)
            oligomeric_state: 2
        end) == @enzyme_reaction(begin
            substrates: A[C], B[N]
            products:   P[C], Q[N]
            allosteric_regulators: X, Y::Inhibitor
            shared_catalytic_site: (A, P), (B, Q)
            oligomeric_state: 2
        end)
    end

    @testset "a mechanism macro's list labels take bare names" begin
        # An atom bracket gets a hint where it goes, and a tagged name the hint to list it
        # bare; another value gets none.
        @test_throws("`substrates:` expects bare names; got S[C]. Atom brackets belong " *
                     "in @enzyme_reaction.", eval(:(@enzyme_mechanism begin
            substrates: S[C]
            products:   P
            steps: begin
                E + S ⇌ ES
                ES <--> EP
                EP ⇌ E + P
            end
        end)))
        @test_throws("`modifiers:` expects bare names; got P::Mod. List `P` bare.",
                     eval(:(@enzyme_mechanism begin
            substrates: S
            products:   P
            modifiers: P::Mod
            steps: begin
                E + S ⇌ ES
                ES <--> EP
                EP ⇌ E + P
                E + P ⇌ EPmod
            end
        end)))
        # Under nested tags the hint still names the metabolite alone.
        @test_throws("`modifiers:` expects bare names; got (P::Mod)::Mod. List `P` bare.",
                     eval(:(@enzyme_mechanism begin
            substrates: S
            products:   P
            modifiers: P::Mod::Mod
            steps: begin
                E + S ⇌ ES
                ES <--> EP
                EP ⇌ E + P
                E + P ⇌ EPmod
            end
        end)))
        # A multiplicity call (`P(2)`, `X(1, 2)`) under the tag leaves the metabolite's
        # name in the hint; a value that is no name gets no hint (the whole message is
        # pinned).
        @test_throws(
            ErrorException("`modifiers:` expects bare names; got P(2)::Mod. " *
                           "List `P` bare."),
            EnzymeRates.DSL._bare_symbols_from_values(Any[:(P(2)::Mod)], :modifiers))
        @test_throws(ErrorException("`modifiers:` expects bare names; got X(1, 2)::Mod. " *
                                    "List `X` bare."),
                     EnzymeRates.DSL._bare_symbols_from_values(Any[:(X(1, 2)::Mod)],
                                                                      :modifiers))
        @test_throws(
            ErrorException("`modifiers:` expects bare names; got 1::Mod."),
            EnzymeRates.DSL._bare_symbols_from_values(Any[:(1::Mod)], :modifiers))
        # A tag with nothing before it names no metabolite, so it gets no hint.
        @test EnzymeRates.DSL._bare_name(Expr(:(::), :OnlyA)) === nothing
        @test_throws(ErrorException("`modifiers:` expects bare names; got ::Mod."),
                     EnzymeRates.DSL._bare_symbols_from_values(
                         Any[Expr(:(::), :Mod)], :modifiers))
        # An operator call is no multiplicity call: its head names no metabolite, so it
        # gets no hint. An atom bracket under the tag gets the atom-bracket hint.
        @test_throws(ErrorException("`modifiers:` expects bare names; got (A + B)::Mod."),
                     EnzymeRates.DSL._bare_symbols_from_values(Any[:((A + B)::Mod)],
                                                                      :modifiers))
        @test_throws(
            ErrorException("`modifiers:` expects bare names; got S[C]::Mod. " *
                           "Atom brackets belong in @enzyme_reaction."),
            EnzymeRates.DSL._bare_symbols_from_values(Any[:(S[C]::Mod)], :modifiers))
        # `Base.isoperator` does not flag the word operators `in` and `isa`, so `in(2)`
        # counts as a multiplicity call; the infix `A in B` and `A isa B` do not (their
        # arguments are no Integers), nor does an operator called with one Integer
        # (`(-(2))`, which prints as `(-2)`): none gets a hint.
        @test_throws(ErrorException("`modifiers:` expects bare names; got (A in B)::Mod."),
                     EnzymeRates.DSL._bare_symbols_from_values(Any[:((A in B)::Mod)],
                                                                      :modifiers))
        @test_throws(ErrorException("`modifiers:` expects bare names; got (A isa B)::Mod."),
                     EnzymeRates.DSL._bare_symbols_from_values(
                         Any[:((A isa B)::Mod)], :modifiers))
        @test_throws(ErrorException("`modifiers:` expects bare names; got (in)(2)::Mod. " *
                                    "List `in` bare."),
                     EnzymeRates.DSL._bare_symbols_from_values(
                         Any[:(in(2)::Mod)], :modifiers))
        @test_throws(ErrorException("`modifiers:` expects bare names; got (-2)::Mod."),
                     EnzymeRates.DSL._bare_symbols_from_values(Any[:((-(2))::Mod)],
                                                                      :modifiers))
    end

    @testset "a plain mechanism rejects each allosteric-only label" begin
        allosteric_only = "is read only for an allosteric mechanism"
        for (line, message) in (
                (:(allosteric_regulators: A), "unknown label `allosteric_regulators:`"),
                (:(oligomeric_state: 2), "`oligomeric_state:` $allosteric_only"),
                (:(regulatory_site(2): A::OnlyA), "`regulatory_site(n):` $allosteric_only"))
            @test_throws message eval(:(@enzyme_mechanism begin
                substrates: S
                products:   P
                $line
                steps: begin
                    E + S ⇌ ES
                    ES <--> EP
                    EP ⇌ E + P
                end
            end))
        end
    end

    @testset "mechanism macros name themselves in body errors" begin
        # A macro argument that is not a `begin ... end` block.
        @test_throws "@enzyme_mechanism: expected a `begin ... end` block" eval(
            :(@enzyme_mechanism S))

        # A missing required label.
        @test_throws "@enzyme_mechanism: `steps:` not specified" eval(
            :(@enzyme_mechanism begin
                substrates: S
                products:   P
            end))
        @test_throws "@enzyme_mechanism: `substrates:` not specified" eval(
            :(@enzyme_mechanism begin
                products: P
            end))

        # A repeated `steps:` line names the label.
        @test_throws "@enzyme_mechanism: `steps:` given more than once" eval(
            :(@enzyme_mechanism begin
                substrates: S
                products:   P
                steps: begin
                    E + S ⇌ ES      :: EqualAI
                    ES <--> EP    :: EqualAI
                    EP ⇌ E + P      :: EqualAI
                end
                steps: begin
                    E + S ⇌ ES      :: EqualAI
                    ES <--> EP    :: EqualAI
                    EP ⇌ E + P      :: EqualAI
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
            @test_throws("@enzyme_mechanism: `regulatory_site(n):` takes one " *
                         "positive Int multiplicity, as in `regulatory_site(2): " *
                         "A::OnlyA, I::OnlyI`; got `$site`.",
                         eval(:(@enzyme_mechanism begin
                substrates: S
                products:   P
                steps: begin
                    E + S ⇌ ES      :: EqualAI
                    ES <--> EP    :: EqualAI
                    EP ⇌ E + P      :: EqualAI
                end
                $site: A::OnlyA
            end)))
        end
    end

    @testset "a tagged tuple of one step is a one-step group" begin
        # A parenthesized group is a tuple of steps, `(step,)` for one: tagged, it is the
        # same mechanism as the bare tagged step.
        one_tuple = @enzyme_mechanism begin
            substrates: F6P
            products:   F16BP
            oligomeric_state: 2

            steps: begin
                E + F6P ⇌ EF6P        :: EqualAI
                EF6P <--> EF16BP    :: EqualAI
                (EF16BP ⇌ E + F16BP,) :: EqualAI
            end

            regulatory_site(4): A::OnlyA
            regulatory_site(4): I::OnlyI
        end
        bare = @enzyme_mechanism begin
            substrates: F6P
            products:   F16BP
            oligomeric_state: 2

            steps: begin
                E + F6P ⇌ EF6P        :: EqualAI
                EF6P <--> EF16BP    :: EqualAI
                EF16BP ⇌ E + F16BP    :: EqualAI
            end

            regulatory_site(4): A::OnlyA
            regulatory_site(4): I::OnlyI
        end
        @test one_tuple == bare
        @test length(_testhelper_group_members(one_tuple.catalytic)) == 3

        # A parenthesized step is no tuple: its tag goes after its right-hand side.
        @test_throws("@enzyme_mechanism: a single step carries its tag after its " *
                     "right-hand side, as in `E + S ⇌ ES :: EqualAI`. Got " *
                     "`(EP ⇌ E + P)::EqualAI`.", eval(:(@enzyme_mechanism begin
            substrates: S
            products:   P
            steps: begin
                E + S ⇌ ES
                ES <--> EP
                (EP ⇌ E + P) :: EqualAI
            end
        end)))
    end

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

    @testset "the macros work in a module that names the package otherwise" begin
        # A module that imports the package under another name, or imports only the
        # macros, calls each macro and gets the value it builds here: an expansion names
        # nothing from the caller's scope.
        calls = [
            :(@enzyme_reaction begin
                substrates: S[C]
                products: P[C]
            end),
            :(@enzyme_mechanism begin
                substrates: S
                products: P
                steps: begin
                    E + S ⇌ ES
                    ES <--> EP
                    EP ⇌ E + P
                end
            end),
            :(@enzyme_mechanism begin
                substrates: S
                products: P
                steps: begin
                    E + S ⇌ ES     :: EqualAI
                    ES <--> EP   :: OnlyA
                    EP ⇌ E + P     :: EqualAI
                end
                regulatory_site(2): X::OnlyI
            end)]
        values = [eval(call) for call in calls]
        renamed, macros_only = Module(), Module()
        Core.eval(renamed, :(import EnzymeRates as ER2))
        Core.eval(macros_only, :(using EnzymeRates: @enzyme_reaction, @enzyme_mechanism))
        for (call, value) in zip(calls, values)
            qualified = Expr(:macrocall, Expr(:., :ER2, QuoteNode(call.args[1])),
                             call.args[2:end]...)
            @test Core.eval(renamed, qualified) == value
            @test Core.eval(macros_only, call) == value
        end
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

    @testset "an undeclared metabolite on a step side reads as a second form" begin
        # A name on a step side that is no declared metabolite names an enzyme form, so an
        # undeclared `ATP` beside the form E is a second form on its side.
        two_forms = "@enzyme_mechanism: step side has more than one enzyme-form term " *
                    "(E, ATP); each elementary step has exactly one enzyme form per " *
                    "side, and a name that is no declared metabolite names an enzyme form."
        @test_throws two_forms eval(
            :(@enzyme_mechanism begin
                substrates: S
                products: P
                steps: begin
                    E + ATP ⇌ EATP
                    ES <--> EP
                    EP ⇌ E + P
                end
            end))
        @test_throws two_forms eval(
            :(@enzyme_mechanism begin
                substrates: S
                products: P
                steps: begin
                    E + ATP ⇌ EATP    :: EqualAI
                    E + S ⇌ ES        :: EqualAI
                    ES <--> EP      :: EqualAI
                    EP ⇌ E + P        :: EqualAI
                end
            end))
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
        @test_throws("@enzyme_reaction: `oligomeric_state:` lists a multiplicity more " *
                     "than once; got (2, 2).", eval(:(@enzyme_reaction begin
            substrates: S[C]
            products:   P[C]
            allosteric_regulators: A
            oligomeric_state: (2, 2)
        end)))
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

    @testset "@enzyme_mechanism" begin
        m = @enzyme_mechanism begin
            substrates: S
            products:   P
            steps: begin
                E + S <--> ES
                ES <--> E + P
            end
        end
        @test length(m.catalytic.steps) == 2
        @test Set(_testhelper_graph(m).forms) == Set(["E", "ES"])
        @test Set(metabolites(m)) == Set([:S, :P])
    end

    @testset "Elementary steps" begin
        # The parser requires an enzyme form on each step side at macro-expansion time.
        # Wrap in `eval(:(...))` so the macro expansion happens at runtime where
        # `@test_throws` can catch the LoadError that wraps the parser's exception.
        @test_throws "no enzyme-form term" eval(:(@enzyme_mechanism begin
            substrates: S
            products:   P
            steps: begin
                S <--> ES
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
                        E + S ⇌ ES
                        $line
                    end
                end))
        end
    end

    @testset "a product bound as a modifier names its form apart" begin
        m = @enzyme_mechanism begin
            substrates: S
            products:   P
            steps: begin
                E + S <--> ES
                ES <--> EP
                EP <--> E + P
                E + P <--> EPmod
            end
        end
        forms = _testhelper_graph(m).forms
        @test "EP" in forms       # product-bound form (catalytic release)
        @test "EPmod" in forms    # inhibitor-bound form — DISTINCT from EP
        @test "EP" != "EPmod"
        # The rate equation reads concs.P for the inhibitor term, not concs.Pmod:
        re_str = rate_equation_string(m)
        @test occursin("P", re_str)             # real metabolite name appears
        @test occursin("(; S, P) = concs", re_str)  # concentration is P, never Pmod
        # The mod marker never names a concentration.
        @test !occursin("Pmod = concs", re_str)
        _testhelper_check_law_against_oracle(m)
    end

    @testset "a dual-role name binds in its catalytic role" begin
        # The role of the metabolite each binding of `m` takes up, by the name of the form
        # it binds into.
        _testhelper_bound_roles(m) = (cm = m.catalytic;
            Dict(String(cm.enzyme_forms[e.to]) =>
                     (only(e.taken_up) in cm.substrates ? :substrate :
                      only(e.taken_up) in cm.products ? :product : :modifier)
                 for e in cm.steps if EnzymeRates.Enumeration._is_binding(e)))
        # A substrate that also binds as a modifier is taken up as the substrate; the
        # form it makes, EAmod, tells the two bindings apart.
        m = @enzyme_mechanism begin
            substrates: A
            products:   P
            steps: begin
                E + A ⇌ EA
                EA <--> EP
                EP ⇌ E + P
                E + A ⇌ EAmod
            end
        end
        roles = _testhelper_bound_roles(m)
        @test roles["EA"] === :substrate
        @test roles["EAmod"] === :substrate
        @test roles["EP"] === :product

        # A substrate and a product also declared as allosteric regulators bind
        # as the substrate and the product in catalytic steps; as regulators
        # they bind only at their regulatory site.
        am = @enzyme_mechanism begin
            substrates: A
            products:   P
            oligomeric_state: 2
            steps: begin
                E + A ⇌ EA          :: EqualAI
                EA <--> EP        :: NonequalAI
                EP ⇌ E + P          :: EqualAI
            end
            regulatory_site(2): A::OnlyI, P::OnlyA
        end
        roles = _testhelper_bound_roles(am)
        @test roles["EA"] === :substrate
        @test roles["EP"] === :product
        site = only(am.allosteric.regulatory_sites)
        @test Set(site.ligands) == Set([:A, :P])
    end

    @testset "an allosteric regulator binds only at its regulatory site" begin
        # R is only a ligand of a regulatory site: a free term of a catalytic step names
        # it in error.
        @test_throws "`R` is an allosteric regulator" eval(:(@enzyme_mechanism begin
            substrates: S
            products:   P
            steps: begin
                E + S ⇌ ES      :: EqualAI
                ES <--> EP    :: EqualAI
                EP ⇌ E + P      :: EqualAI
                E + R ⇌ ER      :: EqualAI
            end
            regulatory_site(1): R::OnlyA
        end))
        # Listed in both, R binds its regulatory site and, as a modifier, the catalytic
        # site.
        dual = @enzyme_mechanism begin
            substrates: S
            products:   P
            modifiers: R
            steps: begin
                E + S ⇌ ES              :: EqualAI
                ES <--> EP            :: EqualAI
                EP ⇌ E + P              :: EqualAI
                E + R ⇌ ER    :: EqualAI
            end
            regulatory_site(1): R::OnlyA
        end
        r = _testhelper_metabolite(dual.catalytic, :R)
        @test count(e -> e.taken_up == [r] && isempty(e.given_off),
                    dual.catalytic.steps) == 1
        @test :R in only(dual.allosteric.regulatory_sites).ligands
    end

    @testset "X::Mod is an error that says to take X up" begin
        # A step that binds a reactant, or a regulator, as a modifier takes the metabolite
        # itself up; a tag `::Mod` asks to write the bare name.
        for x in (:S, :P, :R)
            @test_throws(
                "@enzyme_mechanism: `$x::Mod` writes a modifier copy, which a mechanism " *
                "does not hold: a step that binds $x as a modifier takes $x up; write " *
                "`$x`.", eval(:(@enzyme_mechanism begin
                substrates: S
                products:   P
                steps: begin
                    E + S ⇌ ES
                    ES <--> EP
                    EP ⇌ E + P
                    E + $x::Mod ⇌ $(Symbol(:E, x, :mod))
                end
            end)))
        end
        declared = @enzyme_mechanism begin
            substrates: S
            products:   P
            steps: begin
                E + S ⇌ ES
                ES <--> EP
                EP ⇌ E + P
                E + P ⇌ EPmod
            end
        end
        @test declared.catalytic.metabolites == [:S, :P]
    end

    @testset "a step side holds one enzyme form" begin
        # A side with two enzyme forms is rejected.
        @test_throws("@enzyme_mechanism: step side has more than one enzyme-form term " *
                     "(EA, EB); each elementary step has exactly one enzyme form per " *
                     "side, and a name that is no declared metabolite names an enzyme " *
                     "form.", eval(:(@enzyme_mechanism begin
            substrates: A, B
            products:   P
            steps: begin
                EA + EB <--> EAB
                EAB <--> E + P
            end
        end)))
    end

    @testset "a binding written as its release is one step, its state at the end" begin
        # `ESmod ⇌ E + S :: EqualAI` parses as `S :: EqualAI` on its last term: the state
        # is the step's tag, and the step reads like its binding direction.
        bind_mod = @enzyme_mechanism begin
            substrates: S
            products:   P
            steps: begin
                E + S ⇌ ES             :: EqualAI
                ES <--> EP           :: EqualAI
                EP ⇌ E + P             :: EqualAI
                E + S ⇌ ESmod   :: EqualAI
            end
        end
        release_mod = @enzyme_mechanism begin
            substrates: S
            products:   P
            steps: begin
                E + S ⇌ ES             :: EqualAI
                ES <--> EP           :: EqualAI
                EP ⇌ E + P             :: EqualAI
                ESmod ⇌ E + S   :: EqualAI
            end
        end
        @test release_mod == bind_mod

        # Without a state the step is missing its tag.
        @test_throws("@enzyme_mechanism: the statement `ESmod ⇌ E + S` " *
                     "carries no allosteric state, while `E + S ⇌ ES::EqualAI` carries " *
                     "one; a mechanism is allosteric when every statement ends in its " *
                     "state (`:: EqualAI`, `:: NonequalAI`, `:: OnlyA` or `:: OnlyI`), " *
                     "and plain when none does.", eval(:(@enzyme_mechanism begin
            substrates: S
            products:   P
            steps: begin
                E + S ⇌ ES           :: EqualAI
                ES <--> EP         :: EqualAI
                EP ⇌ E + P           :: EqualAI
                ESmod ⇌ E + S
            end
        end)))

        # A plain mechanism may write the second binding of S as its release.
        bind_plain = @enzyme_mechanism begin
            substrates: S
            products:   P
            steps: begin
                E + S ⇌ ES
                ES <--> EP
                EP ⇌ E + P
                E + S ⇌ ESmod
            end
        end
        release_plain = @enzyme_mechanism begin
            substrates: S
            products:   P
            steps: begin
                E + S ⇌ ES
                ES <--> EP
                EP ⇌ E + P
                ESmod ⇌ E + S
            end
        end
        @test release_plain == bind_plain

        # A step tag that is not a name is rejected.
        @test_throws "must be a Symbol" eval(:(@enzyme_mechanism begin
            substrates: S
            products:   P
            steps: begin
                E + S ⇌ ES    :: EqualAI
                ES <--> EP  :: EqualAI
                EP ⇌ E + P    :: 3
            end
        end))
    end

    @testset "a malformed step is rejected with the walker's message" begin
        # Each statement replaces the S binding of a uni-uni block with modifier I.
        _testhelper_block_with(step) = :(@enzyme_mechanism begin
            substrates: S
            products:   P
            modifiers: I
            steps: begin
                $step
                ES <--> EP
                EP ⇌ E + P
                E + I ⇌ EI
            end
        end)
        for (step, text) in (
                (:(E + S → ES), "expected ⇌ or <--> step operator; got →"),
                (:(E + 2 ⇌ ES),
                 "expected a metabolite or an enzyme form's name on a step side; got `2`."),
                (:(E + 2::Mod ⇌ ES),
                 "expected a metabolite or an enzyme form's name on a step side; got " *
                 "`2::Mod`."))
            # A syntax error throws at expansion, so `eval` wraps it in a LoadError.
            @test_throws "@enzyme_mechanism: " * text eval(_testhelper_block_with(step))
        end

        # Call notation is an error. Where every argument is a bare name the error names
        # the form, its head and then the arguments sorted; otherwise it asks for a name.
        for (step, name) in ((:(E + S ⇌ E(S)), "ES"), (:(E + S ⇌ E(B, A)), "EAB"))
            @test_throws("@enzyme_mechanism: `$(step.args[3])` is call notation; an " *
                         "enzyme form is a name, so write `$name`.",
                         eval(_testhelper_block_with(step)))
        end
        for step in (:(E + S ⇌ E(S, f(I))), :(E + S ⇌ E(S; residual = 2S)),
                     :(E + S ⇌ E(S; foo = P)), :(E + S ⇌ E(S; foo)),
                     :(E + S ⇌ E((I + S)::Mod)), :(E + S ⇌ E(S::Mod)),
                     :(E + S ⇌ E(B; residual = A - P)))
            @test_throws("@enzyme_mechanism: `$(step.args[3])` is call notation; an " *
                         "enzyme form is a name, so write one, such as `EAB` for " *
                         "`E(A, B)`.",
                         eval(_testhelper_block_with(step)))
        end
        # A count is an error: a step names a metabolite once for each copy.
        @test_throws("@enzyme_mechanism: `2S` writes a count; a step names a metabolite " *
                     "once for each copy it takes up or gives off, as in " *
                     "`E + S + S ⇌ ESS`.",
                     eval(_testhelper_block_with(:(E + 2S ⇌ ES))))

        # A step group's tag is a name.
        @test_throws("@enzyme_mechanism: step-group tag must be a Symbol; got 1",
                     eval(:(@enzyme_mechanism begin
            substrates: S
            products:   P
            steps: begin
                (E + S ⇌ ES, ES <--> EP) :: 1
                EP ⇌ E + P                   :: EqualAI
            end
        end)))
        # A step inside a group carries no state of its own: the group's state follows its
        # closing parenthesis.
        @test_throws("@enzyme_mechanism: the step `ES + I ⇌ EIS::EqualAI` carries a " *
                     "state inside a parenthesized group; a group's state is written " *
                     "after its closing parenthesis, as in " *
                     "`(E + I ⇌ EI, ES + I ⇌ EIS) :: EqualAI`.",
                     eval(:(@enzyme_mechanism begin
            substrates: S
            products:   P
            modifiers: I
            steps: begin
                E + S ⇌ ES                           :: EqualAI
                ES <--> EP                           :: OnlyA
                EP ⇌ E + P                           :: EqualAI
                (E + I ⇌ EI, ES + I ⇌ EIS :: EqualAI)
            end
        end)))
    end

    @testset "the guard judges :OnlyI on a step that does more than bind a modifier" begin
        # Three shapes, each :OnlyI: a reactant binding, an isomerization of an inhibitor
        # complex, and a Theorell–Chance step that takes up a modifier and gives off a
        # product. The active conformation cannot switch off the rapid-equilibrium S
        # binding, whose two sides the product release still joins. The guard's built
        # readings stop the isomerization, which leaves FI an island; the search lets
        # it lose its rate into FI instead and admits the mechanism. The Theorell–Chance
        # step is stopped. Each admitted shape has a law equal to the limit on every face,
        # and each group tagged :NonequalAI passes the guard.
        D = EnzymeRates.Derivation
        rng = Random.MersenneTwister(31)
        _testhelper_block(steps; guarded = true) = Meta.parse(
            (guarded ? "" : "@_testhelper_unguarded ") * "@enzyme_mechanism begin\n" *
            "substrates: S\nproducts: P\nmodifiers: I\nsteps: begin\n$steps\nend\nend")
        rest = "ES <--> E + P :: EqualAI\nE + I ⇌ EI :: EqualAI"
        binding = "E + S ⇌ ES :: OnlyI\n$rest"
        isomerization = "E + S ⇌ ES :: EqualAI\n$rest\nEI <--> FI :: OnlyI"
        exchange = "E + S ⇌ ES :: EqualAI\n$rest\nES + I <--> EI + P :: OnlyI"
        am = eval(_testhelper_block(binding; guarded = false))
        v = D._guard(am)
        @test (v.fault, v.conformation) == (D._FAST_STEP_PINNED, 1)
        @test EnzymeRates.Mechanisms.step_text(am.catalytic, v.edge) == "E + S ⇌ ES"
        @test_throws(ErrorException(D._guard_message(am)),
                     eval(_testhelper_block(binding)))
        am = eval(_testhelper_block(isomerization; guarded = false))
        v = D._admission_interior(am, D._built_reading(am))
        @test (v.fault, v.conformation) == (D._ISLAND, 1)
        @test EnzymeRates.Mechanisms.step_text(am.catalytic, v.edge) == "EI <--> FI"
        @test eval(_testhelper_block(isomerization)) == am
        @test all(_testhelper_admission_agreement(am; rng = Random.MersenneTwister(32),
                                                  corners = true))
        _testhelper_check_law_against_limit(am)
        admitted = eval(_testhelper_block(exchange))
        @test all(_testhelper_admission_agreement(admitted; rng))
        @test all(==(:agree), _testhelper_respeed_verdicts(admitted; rng))
        for text in (binding, isomerization, exchange)
            @test _testhelper_is_allosteric(
                eval(_testhelper_block(replace(text, ":: OnlyI" => ":: NonequalAI"))))
        end
    end
end

@testset "DSL: modifiers and regulatory_site(n)" begin
    @testset "a modifier binds bare, and a reactant binds as a modifier by its name" begin
        # I binds only as a modifier; P is a product that also binds as a modifier, into
        # the form EPmod.
        m = @enzyme_mechanism begin
            substrates: S
            products:   P
            modifiers:  I
            steps: begin
                E + S ⇌ ES
                ES <--> EP
                EP ⇌ E + P
                ES + I ⇌ EIS
                E + P ⇌ EPmod
            end
        end
        @test parameters(m) == (:K_EIS_to_ES_I, :K_EP_to_E_P,
            :K_EPmod_to_E_P, :K_ES_to_E_S, :k_ES_to_EP)
        @test metabolites(m) == [:S, :P, :I]
        _testhelper_check_law_against_oracle(m)
    end

    @testset "the allosteric grammar: one tag per statement, ligands at their sites" begin
        am = @enzyme_mechanism begin
            substrates: A, B
            products:   P, Q
            modifiers:  I
            oligomeric_state: 2
            steps: begin
                E + A ⇌ EA                              :: EqualAI
                EA + B ⇌ EAB                        :: NonequalAI
                EAB <--> EPQ                      :: OnlyA
                EPQ ⇌ EQ + P                        :: EqualAI
                EQ ⇌ E + Q                              :: EqualAI
                (E + I ⇌ EI, EA + I ⇌ EAI)        :: EqualAI
                EQ + P ⇌ EPmodQ              :: EqualAI
            end
            regulatory_site(2): X::OnlyA, Y::OnlyI
        end
        @test _testhelper_is_allosteric(am)
        @test am.allosteric.oligomeric_state == 2
        @test parameters(am) == (:K_A_EAB_to_EA_B, :K_EA_to_E_A,
            :K_EI_to_E_I, :K_EPQ_to_EQ_P, :K_EPmodQ_to_EQ_P, :K_EQ_to_E_Q,
            :K_I_EAB_to_EA_B, :k_A_EAB_to_EPQ, :K_A_Xreg, :K_I_Yreg, :L)
        @test metabolites(am) == [:A, :B, :P, :Q, :I, :X, :Y]
        _testhelper_check_law_against_oracle(am)
    end

    @testset "oligomeric_state: defaults to 1" begin
        am = @enzyme_mechanism begin
            substrates: S
            products:   P
            steps: begin
                E + S ⇌ ES      :: EqualAI
                ES <--> EP    :: OnlyA
                EP ⇌ E + P      :: EqualAI
            end
            regulatory_site(1): X::OnlyI
        end
        @test am.allosteric.oligomeric_state == 1
        @test parameters(am) ==
              (:K_EP_to_E_P, :K_ES_to_E_S, :k_A_ES_to_EP, :K_I_Xreg, :L)
        _testhelper_check_law_against_oracle(am)
    end

    @testset "a name can be a modifier and an allosteric regulator" begin
        # R binds the catalytic site as a modifier and its regulatory site as a ligand.
        dual = @enzyme_mechanism begin
            substrates: S
            products:   P
            modifiers:  R
            steps: begin
                E + S ⇌ ES      :: EqualAI
                ES <--> EP    :: EqualAI
                EP ⇌ E + P      :: EqualAI
                E + R ⇌ ER      :: EqualAI
            end
            regulatory_site(1): R::OnlyA
        end
        @test parameters(dual) ==
              (:K_EP_to_E_P, :K_ER_to_E_R, :K_ES_to_E_S, :k_ES_to_EP, :K_A_Rreg, :L)
        @test metabolites(dual) == [:S, :P, :R]
        _testhelper_check_law_against_oracle(dual)
    end

    @testset "the walker returns each block as written" begin
        # Steps keep their written order and direction, each joining its two written
        # forms, named in the order the sides first write them, and taking up and giving
        # off metabolites by their place among the declared names; each statement is one
        # kinetic group label.
        St = EnzymeRates.Mechanisms.Step
        none = Int[]
        raw = EnzymeRates.DSL.raw_mechanism(quote
            substrates: S
            products:   P
            modifiers:  I
            steps: begin
                ES ⇌ E + S
                (E + I ⇌ EI, ES + I ⇌ EIS)
                ES <--> Estar_res_plusS_minusP + P
                E + S ⇌ ESmod
            end
        end)
        @test raw == (; substrates = [:S], products = [:P], modifiers = [:I],
            forms = [:ES, :E, :EI, :EIS, :Estar_res_plusS_minusP, :ESmod],
            steps = [St(1, 2, none, [1], true, 1), St(2, 3, [3], none, true, 2),
                     St(1, 4, [3], none, true, 2), St(1, 5, none, [2], false, 3),
                     St(2, 6, [1], none, true, 4)],
            allosteric_states = Symbol[], oligomeric_state = 1,
            regulatory_sites = EnzymeRates.DSL._WrittenSite[])
        # An allosteric block adds a tag per group label, the oligomeric state
        # (1 when omitted) and its sites, each ligand list as written.
        allo = EnzymeRates.DSL.raw_mechanism(quote
            substrates: S
            products:   P
            steps: begin
                E + S ⇌ ES              :: EqualAI
                ES <--> EP            :: OnlyA
                (EP ⇌ E + P,)           :: NonequalAI
            end
            regulatory_site(2): Y::OnlyI, X::OnlyA
        end)
        @test allo == (; substrates = [:S], products = [:P], modifiers = Symbol[],
            forms = [:E, :ES, :EP],
            steps = [St(1, 2, [1], none, true, 1), St(2, 3, none, none, false, 2),
                     St(3, 1, none, [2], true, 3)],
            allosteric_states = [:EqualAI, :OnlyA, :NonequalAI], oligomeric_state = 1,
            regulatory_sites = [(; ligands = [:Y, :X], multiplicity = 2,
                                 allosteric_states = [:OnlyI, :OnlyA])])
    end

    @testset "a tag on a metabolite is an error that says what to write" begin
        # A modifier, a reactant: `::Mod` asks to write the bare name.
        for (x, declared) in ((:I, :(modifiers: I)), (:P, :(modifiers: J)))
            @test_throws("@enzyme_mechanism: `$x::Mod` writes a modifier copy, which a " *
                         "mechanism does not hold: a step that binds $x as a modifier " *
                         "takes $x up; write `$x`.", eval(:(@enzyme_mechanism begin
                substrates: S
                products:   P
                $declared
                steps: begin
                    E + S ⇌ ES
                    ES <--> EP
                    EP ⇌ E + P
                    E + $x::Mod ⇌ $(Symbol(:E, x, :mod))
                end
            end)))
        end
        # An allosteric regulator binds only at its site.
        @test_throws "`R` is an allosteric regulator" eval(:(@enzyme_mechanism begin
            substrates: S
            products:   P
            steps: begin
                E + S ⇌ ES      :: EqualAI
                ES <--> EP    :: EqualAI
                EP ⇌ E + P      :: EqualAI
                E + R ⇌ ER      :: EqualAI
            end
            regulatory_site(1): R::OnlyA
        end))
        # Another tag on a left-hand metabolite is unknown; on the last right-hand term it
        # is the step's state, so one statement carries a state and the others none.
        @test_throws("@enzyme_mechanism: unknown tag ::Inh on P; a step ends in its " *
                     "allosteric state, written after its right-hand side.",
                     eval(:(@enzyme_mechanism begin
            substrates: S
            products:   P
            steps: begin
                E + S ⇌ ES
                ES <--> EP
                EP ⇌ E + P
                E + P::Inh ⇌ EPInh
            end
        end)))
        @test_throws("carries no allosteric state, while `EP ⇌ E + P::Inh` carries one",
                     eval(:(@enzyme_mechanism begin
            substrates: S
            products:   P
            steps: begin
                E + S ⇌ ES
                ES <--> EP
                EP ⇌ E + P::Inh
            end
        end)))
    end

    @testset "every allosteric regulator sits at exactly one regulatory site" begin
        # The mechanism `@enzyme_mechanism` builds from the block below with the site
        # lines `sites` appended.
        declarations = quote
            substrates: S
            products:   P
            steps: begin
                E + S ⇌ ES      :: EqualAI
                ES <--> EP    :: OnlyA
                EP ⇌ E + P      :: EqualAI
            end
        end
        _testhelper_sites(sites...) = eval(:(@enzyme_mechanism $(
            Expr(:block, declarations.args..., sites...))))
        @test_throws("@enzyme_mechanism: allosteric regulator `X` is listed 2 times " *
                     "on `regulatory_site(n):` lines; an allosteric regulator sits at " *
                     "exactly one regulatory site.",
                     _testhelper_sites(:(regulatory_site(2): X::OnlyA, Y::OnlyI),
                                       :(regulatory_site(4): X::OnlyA)))
        @test_throws("@enzyme_mechanism: allosteric regulator `X` is listed 2 times",
                     _testhelper_sites(:(regulatory_site(2): X::OnlyA, Y::OnlyI,
                                                             X::OnlyI)))
        # A site lists each ligand with its tag after a positive multiplicity.
        @test_throws("@enzyme_mechanism: `regulatory_site(n):` lists each ligand " *
                     "with its tag, as in `regulatory_site(2): A::OnlyA, I::OnlyI`; " *
                     "got `X`.", _testhelper_sites(:(regulatory_site(2): X, Y::OnlyI)))
        @test_throws("`regulatory_site` multiplicity must be a positive Int, got 0.",
                     _testhelper_sites(:(regulatory_site(0): X::OnlyA, Y::OnlyI)))
        @test_throws("@enzyme_mechanism: `regulatory_site(n):` lists each ligand " *
                     "with its tag, as in `regulatory_site(2): A::OnlyA, I::OnlyI`; " *
                     "got a `begin ... end` block.",
                     _testhelper_sites(:(regulatory_site(2): begin
                                             ligands: X, Y
                                         end)))
    end

    @testset "every statement carries an allosteric state, or none does" begin
        mixed(without, with) =
            "@enzyme_mechanism: the statement `$without` carries no allosteric state, " *
            "while `$with` carries one; a mechanism is allosteric when every statement " *
            "ends in its state (`:: EqualAI`, `:: NonequalAI`, `:: OnlyA` or " *
            "`:: OnlyI`), and plain when none does."
        @test_throws(mixed("ES <--> EP", "E + S ⇌ ES::EqualAI"),
                     eval(:(@enzyme_mechanism begin
            substrates: S
            products:   P
            steps: begin
                E + S ⇌ ES      :: EqualAI
                ES <--> EP
                EP ⇌ E + P      :: EqualAI
            end
        end)))
        @test_throws(mixed("E + S ⇌ ES", "ES <--> EP::EqualAI"),
                     eval(:(@enzyme_mechanism begin
            substrates: S
            products:   P
            steps: begin
                E + S ⇌ ES
                ES <--> EP    :: EqualAI
                EP ⇌ E + P
            end
        end)))
        # A single step carries its tag after its right-hand side; a parenthesized group
        # is a tuple of steps, `(step,)` for one.
        @test_throws("@enzyme_mechanism: a single step carries its tag after its " *
                     "right-hand side, as in `E + S ⇌ ES :: EqualAI`. Got " *
                     "`(EP ⇌ E + P)::EqualAI`.", eval(:(@enzyme_mechanism begin
            substrates: S
            products:   P
            steps: begin
                E + S ⇌ ES              :: EqualAI
                ES <--> EP            :: EqualAI
                (EP ⇌ E + P)            :: EqualAI
            end
        end)))
        # A tag after a parenthesized right-hand side tags the step, as after a bare one.
        St = EnzymeRates.Mechanisms.Step
        none = Int[]
        @test EnzymeRates.DSL.raw_mechanism(quote
            substrates: S
            products:   P
            steps: begin
                E + S ⇌ ES        :: EqualAI
                ES <--> EP      :: OnlyA
                EP ⇌ (E + P)      :: EqualAI
            end
        end) == (; substrates = [:S], products = [:P], modifiers = Symbol[],
            forms = [:E, :ES, :EP],
            steps = [St(1, 2, [1], none, true, 1), St(2, 3, none, none, false, 2),
                     St(3, 1, none, [2], true, 3)],
            allosteric_states = [:EqualAI, :OnlyA, :EqualAI], oligomeric_state = 1,
            regulatory_sites = EnzymeRates.DSL._WrittenSite[])
        @test_throws "@enzyme_mechanism: a parenthesized group holds no step" eval(
            :(@enzyme_mechanism begin
                substrates: S
                products:   P
                steps: begin
                    E + S ⇌ ES
                    ()
                end
            end))
    end

    @testset "labels: the valid ones, each single-valued label once with one value" begin
        valid = "Valid labels: substrates, products, modifiers, oligomeric_state, steps, " *
                "regulatory_site(n)."
        for label in (:regulators, :catalytic_steps, :catalytic_inhibitors)
            @test_throws("@enzyme_mechanism: unknown label `$label:`. " * valid,
                         eval(:(@enzyme_mechanism begin
                substrates: S
                products:   P
                $label: I
                steps: begin
                    E + S ⇌ ES
                end
            end)))
        end
        for label in (:catalytic_steps, :catalytic_inhibitors, :regulators)
            @test_throws("@enzyme_mechanism: unknown label `$label:`. " * valid,
                         eval(:(@enzyme_mechanism begin
                substrates: S
                products:   P
                $label: I
                steps: begin
                    E + S ⇌ ES      :: EqualAI
                end
            end)))
        end
        @test_throws("@enzyme_mechanism: `oligomeric_state:` takes one " *
                     "positive Int; got 2, 4.", eval(Meta.parse("""
            @enzyme_mechanism begin
                substrates: S
                products:   P
                oligomeric_state: 2, 4
                steps: begin
                    E + S ⇌ ES      :: EqualAI
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
                    E + S ⇌ ES
                end, S
            end""")))
        @test_throws "`steps:` given more than once" eval(:(@enzyme_mechanism begin
            substrates: S
            products:   P
            steps: begin
                E + S ⇌ ES
            end
            steps: begin
                ES <--> E + P
            end
        end))
        @test_throws "`oligomeric_state:` given more than once" eval(
            :(@enzyme_mechanism begin
                substrates: S
                products:   P
                oligomeric_state: 2
                oligomeric_state: 2
                steps: begin
                    E + S ⇌ ES      :: EqualAI
                end
            end))
    end

    @testset "a form's name is opaque: capitals, lowercase letters and digits" begin
        m = @enzyme_mechanism begin
            substrates: S
            products:   P
            steps: begin
                E + S ⇌ ES
                ES <--> Estar2P
                Estar2P ⇌ Estar2 + P
                Estar2 <--> E
            end
        end
        @test parameters(m) == (:K_ES_to_E_S, :K_Estar2P_to_Estar2_P,
            :k_ES_to_Estar2P, :k_E_to_Estar2, :k_Estar2_to_E)
        _testhelper_check_law_against_oracle(m)
    end

    @testset "a name on a step side is a metabolite before it is a form" begin
        # A declared metabolite's name on a side is the metabolite, so a side of two
        # metabolites names no form; a ligand's name on a side is the ligand, which binds
        # only at its regulatory site.
        @test_throws("@enzyme_mechanism: step side has no enzyme-form term (terms: S)",
                     eval(:(@enzyme_mechanism begin
            substrates: S
            products:   P
            steps: begin
                E + S ⇌ ES
                ES <--> S
                ES ⇌ E + P
            end
        end)))
        @test_throws("@enzyme_mechanism: `FP` is an allosteric regulator; it binds only " *
                     "at its regulatory site, unless `modifiers:` lists it too.",
                     eval(:(@enzyme_mechanism begin
            substrates: S
            products:   P
            steps: begin
                E + S ⇌ ES      :: EqualAI
                ES <--> FP    :: EqualAI
                FP ⇌ EP       :: EqualAI
                EP ⇌ E + P      :: EqualAI
            end
            regulatory_site(1): FP::OnlyI
        end)))
    end

    @testset "constructor errors throw when the expanded code runs" begin
        # The walker reads the block at expansion; the constructors and the guard check
        # the mechanism when the macro's code runs, so these throw outside `eval`.
        no_turnover = "Mechanism: without its :OnlyI steps ES <--> EP the " *
                      "active conformation turns no substrate into product"
        @test_throws no_turnover @enzyme_mechanism begin
            substrates: S
            products:   P
            steps: begin
                E + S ⇌ ES      :: EqualAI
                ES <--> EP    :: OnlyI
                EP ⇌ E + P      :: EqualAI
            end
        end
        @test_throws "allosteric state Bogus must be one of" @enzyme_mechanism begin
            substrates: S
            products:   P
            steps: begin
                E + S ⇌ ES      :: EqualAI
                ES <--> EP    :: EqualAI
                EP ⇌ E + P      :: EqualAI
            end
            regulatory_site(2): X::Bogus
        end
        @test_throws "S named as both a substrate and a product" @enzyme_mechanism begin
            substrates: S
            products:   S
            steps: begin
                E + S ⇌ ES
            end
        end
    end
end

@testset "DSL text: parse_mechanism reads a mechanism back from its text" begin
    @testset "every mechanism the macros accept reads back" begin
        # Names a block may use, each through the macro, then printed and read back:
        # contextual keywords, Greek letters and a subscript as names; form names with
        # digits, of one letter, or that once named tags (`ModP`, `Inh`); forms named
        # after a residual, with and without bound metabolites, one that only subtracts;
        # a reactant bound twice as a modifier; a name that is both a modifier and an
        # allosteric regulator, a substrate that is also an allosteric regulator; an
        # allosteric mechanism without sites; multiplicities of two digits.
        mechanisms = [
            @enzyme_mechanism(begin
                substrates: abstract, mutable
                products: outer, primitive
                modifiers: public, type, var
                steps: begin
                    E + abstract ⇌ Eabstract
                    Eabstract + mutable ⇌ Eabstractmutable
                    Eabstractmutable <--> Eouterprimitive
                    Eouterprimitive ⇌ Eprimitive + outer
                    Eprimitive ⇌ E + primitive
                    E + public ⇌ Epublic
                    E + type ⇌ Etype
                    E + var ⇌ Evar
                end
            end),
            @enzyme_mechanism(begin
                substrates: α
                products: β₂
                modifiers: Δ
                steps: begin
                    E + α ⇌ Eα
                    Eα <--> Eβ₂
                    Eβ₂ ⇌ E + β₂
                    E + Δ ⇌ EΔ
                end
            end),
            @enzyme_mechanism(begin
                substrates: S
                products: P
                steps: begin
                    E + S <--> ES
                    ES <--> Z9S
                    Z9S <--> ModP
                    ModP <--> Inh + P
                    Inh <--> F
                    F <--> E
                end
            end),
            @enzyme_mechanism(begin
                substrates: A, B
                products: P, Q
                steps: begin
                    E + A <--> EA
                    EA <--> F_res_plusA_minusP + P
                    F_res_plusA_minusP + B <--> FB_res_plusA_minusP
                    FB_res_plusA_minusP <--> EQ
                    EQ <--> E + Q
                    EQ <--> Estar_res_minusQ
                end
            end),
            @enzyme_mechanism(begin
                substrates: S
                products: P
                steps: begin
                    E + S ⇌ ES
                    ES <--> EP
                    EP ⇌ E + P
                    ES + S ⇌ ESSmod
                    ESSmod + S ⇌ ESSmodSmod
                end
            end),
            @enzyme_mechanism(begin
                substrates: S
                products: P
                modifiers: I
                oligomeric_state: 12
                steps: begin
                    E + S ⇌ ES                         :: EqualAI
                    ES <--> EP                       :: OnlyA
                    EP ⇌ E + P                         :: EqualAI
                    (E + I ⇌ EI, ES + I ⇌ EIS)   :: NonequalAI
                end
                regulatory_site(16): I::OnlyI, S::EqualAI
            end),
            @enzyme_mechanism(begin
                substrates: S
                products: P
                steps: begin
                    E + S ⇌ ES     :: EqualAI
                    ES <--> EP   :: OnlyA
                    EP ⇌ E + P     :: EqualAI
                end
            end)]
        for m in mechanisms
            @test _testhelper_reads_back(m)
            @test EnzymeRates.DSL.parse_mechanism(sprint(show, m)) == m
        end
        # A form prints by its name, as written.
        @test contains(sprint(show, mechanisms[4]), "EQ <--> Estar_res_minusQ")
    end

    @testset "the text is one mechanism block, read without evaluation" begin
        text = """
            @enzyme_mechanism begin
                substrates: S
                products: P
                steps: begin
                    E + P ⇌ EP
                    E + S ⇌ ES
                    ES <--> EP
                end
            end"""
        @test EnzymeRates.DSL.parse_mechanism(text) ==
            @enzyme_mechanism begin
                substrates: S
                products: P
                steps: begin
                    E + S ⇌ ES
                    ES <--> EP
                    EP ⇌ E + P
                end
            end
        # Blank lines around the block, as a text copied from a CSV cell may hold, are
        # whitespace, not text.
        @test EnzymeRates.DSL.parse_mechanism("\n  " * text * "\n\n") ==
              EnzymeRates.DSL.parse_mechanism(text)
        # The text as a user reads it back from a saved results file.
        buf = IOBuffer()
        CSV.write(buf, DataFrame(mechanism = [text]))
        seekstart(buf)
        @test EnzymeRates.DSL.parse_mechanism(only(CSV.read(buf, DataFrame).mechanism)) ==
              EnzymeRates.DSL.parse_mechanism(text)
        # Anything but one unqualified mechanism macro call is no mechanism text.
        for other in ("1 + 2", "EnzymeRates." * text,
                      "@enzyme_reaction begin\n    substrates: S[C]\n" *
                      "    products: P[C]\nend")
            ex = Meta.parse(other)
            @test_throws(ErrorException("parse_mechanism: expected the text of one " *
                                        "@enzyme_mechanism block; got $ex"),
                         EnzymeRates.DSL.parse_mechanism(other))
        end
        @test_throws Meta.ParseError EnzymeRates.DSL.parse_mechanism(text * "\nE")
        # An interpolation stays unevaluated: the walker meets an expression where a
        # name belongs.
        unevaluated = replace(text, "substrates: S" =>
                                    "substrates: \$(error(\"evaluated\"))")
        @test_throws("`substrates:` expects bare names; got \$(error(\"evaluated\")).",
                     EnzymeRates.DSL.parse_mechanism(unevaluated))
    end

    @testset "a text the macros reject is rejected with their message" begin
        # A syntax error, constructor errors and the boundary check, each as the macro
        # reports it, thrown at once rather than wrapped in a LoadError.
        @test_throws(ErrorException("@enzyme_mechanism: `steps:` not specified."),
                     EnzymeRates.DSL.parse_mechanism("""
            @enzyme_mechanism begin
                substrates: S
                products: P
            end"""))
        @test_throws(ArgumentError("Mechanism: no enzyme form is named `E`; the free " *
                                   "enzyme is the form named `E`"),
                     EnzymeRates.DSL.parse_mechanism("""
            @enzyme_mechanism begin
                substrates: S
                products: P
                steps: begin
                    F + S ⇌ FS
                    FS <--> FP
                    FP ⇌ F + P
                end
            end"""))
        bottomless = """
            @enzyme_mechanism begin
                substrates: S
                products: P, Q
                steps: begin
                    E + S ⇌ ES
                    E + P <--> EP
                    E + Q <--> EQ
                    EP + Q ⇌ EPQ
                    EQ + P ⇌ EPQ
                    ES <--> EPQ
                end
            end"""
        macro_error = try
            eval(Meta.parse(bottomless))
            nothing
        catch e
            e
        end
        @test macro_error isa ErrorException &&
              occursin("rapid-equilibrium segment", macro_error.msg)
        @test_throws macro_error EnzymeRates.DSL.parse_mechanism(bottomless)
    end
end

@testset "DSL entry points: reaction_expr and raw_mechanism" begin
    M, DSL = EnzymeRates.Mechanisms, EnzymeRates.DSL
    # The `EnzymeReaction` expression a block writes builds the reaction the macro builds;
    # its names resolve in `DSL`.
    @test Core.eval(DSL, DSL.reaction_expr(:(begin
        substrates: S[C]
        products: P[C]
        competitive_inhibitors: I
    end))) == @enzyme_reaction begin
        substrates: S[C]
        products: P[C]
        competitive_inhibitors: I
    end

    # A mechanism block as written: each step in written order and direction between its
    # two written forms, one group label per statement.
    none = Int[]
    plain = (; substrates = [:S], products = [:P], modifiers = Symbol[],
             forms = [:E, :ES, :EP],
             steps = [M.Step(1, 2, [1], none, true, 1), M.Step(2, 3, none, none, false, 2),
                      M.Step(3, 1, none, [2], true, 3)],
             allosteric_states = Symbol[], oligomeric_state = 1,
             regulatory_sites = DSL._WrittenSite[])
    @test DSL.raw_mechanism(:(begin
        substrates: S
        products: P
        steps: begin
            E + S ⇌ ES
            ES <--> EP
            EP ⇌ E + P
        end
    end)) == plain
    allosteric = DSL.raw_mechanism(:(begin
        substrates: S
        products: P
        oligomeric_state: 2
        steps: begin
            E + S ⇌ ES      :: NonequalAI
            ES <--> EP    :: OnlyA
            EP ⇌ E + P      :: EqualAI
        end
        regulatory_site(2): R::OnlyI
    end))
    @test allosteric == (; plain..., allosteric_states = [:NonequalAI, :OnlyA, :EqualAI],
                         oligomeric_state = 2,
                         regulatory_sites = [(; ligands = [:R], multiplicity = 2,
                                              allosteric_states = [:OnlyI])])


    # The constructors build the value a block writes (`_construct`) from its declared
    # names, forms and steps, and check the names: each is listed once, and some step
    # takes up or gives off each.
    @test DSL._construct(plain) == @enzyme_mechanism begin
        substrates: S
        products: P
        steps: begin
            E + S ⇌ ES
            ES <--> EP
            EP ⇌ E + P
        end
    end
    @test_throws(ErrorException("Mechanism: duplicate modifier names"),
                 @enzyme_mechanism begin
        substrates: S
        products: P
        modifiers: I, I
        steps: begin
            E + S ⇌ ES
            ES <--> EP
            EP ⇌ E + P
            E + I ⇌ EI
        end
    end)
    @test_throws(ErrorException("Mechanism: modifier I is declared, but no step binds it"),
                 DSL._construct(DSL.raw_mechanism(:(begin
        substrates: S
        products: P
        modifiers: I
        steps: begin
            E + S ⇌ ES
            ES <--> EP
            EP ⇌ E + P
            E + P ⇌ EPmod
        end
    end))))
    # A reactant declared again under `modifiers:`: a step that binds it as a modifier
    # takes the reactant itself up.
    @test_throws(ErrorException("Mechanism: P is declared both as a reactant and as a " *
                                "modifier; a step that binds a reactant as a modifier " *
                                "simply takes it up, as `ES + S ⇌ ESS` does"),
                 @enzyme_mechanism begin
        substrates: S
        products: P
        modifiers: P
        steps: begin
            E + S ⇌ ES
            ES <--> EP
            EP ⇌ E + P
        end
    end)
end
