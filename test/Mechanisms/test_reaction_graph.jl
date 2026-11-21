# ABOUTME: Tests of the reaction-graph mechanism types: identity, canonical form and
# ABOUTME: checks, its sub-graph reads and the speed of its constructor.

using Statistics: median

const ER = EnzymeRates

"""Uni-uni with rapid-equilibrium bindings, written field by field in canonical form:
metabolites S and P, forms E, EP, ES; the bindings of P and S, then ES → EP."""
function _testhelper_graph_uniuni()
    M = ER.Mechanisms
    none = Int[]
    _testhelper_plain_mechanism(
        [:S, :P], [:E, :EP, :ES], [1], [2],
        [M.Step(1, 2, [2], none, true, 1), M.Step(1, 3, [1], none, true, 2),
         M.Step(3, 2, none, none, false, 3)])
end

"""Whether a value of type `T` hashes the same in every process: a symbol, an integer or
a flag, a vector of such values, or one of the types whose hash folds its name and its
fields (`_testhelper_hash_folds`), whose fields are such values in turn."""
function _testhelper_hash_is_content(T)
    T in (Symbol, Bool, Int) && return true
    T <: Vector && return _testhelper_hash_is_content(eltype(T))
    T in (EnzymeRates.Mechanisms.RegulatorySite, ER.Mechanisms.Step,
          ER.Mechanisms.CatalyticScheme, ER.Mechanisms.Allostery, ER.Mechanism) ||
        return false
    all(f -> _testhelper_hash_is_content(fieldtype(T, f)), _testhelper_identity_fields(T))
end

"""Whether `x` hashes its type's name, then the fields it is compared by
(`_testhelper_identity_fields`) in declaration order."""
_testhelper_hash_folds(x, h) =
    hash(x, h) == foldl((acc, f) -> hash(getfield(x, f), acc),
                        _testhelper_identity_fields(typeof(x));
                        init = hash(nameof(typeof(x)), h))

@testset "a graph hashes its content, the same in every process" begin
    # Julia hashes symbols, integers, flags and vectors of them alike in every process,
    # but a struct's default hash mixes in its type's identity, which changes between
    # precompile sessions. So each graph type hashes its name, then its fields, and its
    # fields hold only such values.
    m = _testhelper_graph_uniuni()
    c = m.catalytic
    site = EnzymeRates.Mechanisms.RegulatorySite([:R], 2, [:OnlyI])
    allostery(tags, n, sites) = ER.Mechanisms.Allostery(tags, n, sites)
    am = ER.Mechanism(c, allostery([:EqualAI, :EqualAI, :OnlyA], 2, [site]))
    h = UInt(17)
    for x in (c.steps[3], c, m.allosteric, m, am.allosteric, am, site)
        @test _testhelper_hash_folds(x, h)
    end
    @test all(_testhelper_hash_is_content,
              (ER.Mechanisms.Step, ER.Mechanisms.CatalyticScheme, ER.Mechanisms.Allostery,
               ER.Mechanism))
    # Each field takes part in `==` and `hash`.
    @test m == _testhelper_graph_uniuni() && hash(m) == hash(_testhelper_graph_uniuni())
    fields = _testhelper_identity_fields(ER.Mechanisms.CatalyticScheme)
    for f in fields
        values = Any[getfield(c, g) for g in fields]
        i = findfirst(==(f), fields)
        values[i] = f === :substrates ? [2] : f === :products ? [1] :
                    reverse(getfield(c, f))
        other = ER.Mechanisms._unchecked_scheme(values...)
        @test other != c && hash(other) != hash(c)
    end
    # The steps' groups take part too.
    regrouped = ER.Mechanisms._unchecked_scheme(c.metabolites, c.enzyme_forms, c.substrates,
        c.products,
        [ER.Mechanisms.Step(e.from, e.to, e.taken_up, e.given_off, e.rapid_equilibrium, g)
         for (e, g) in zip(c.steps, [1, 3, 2])])
    @test regrouped != c && hash(regrouped) != hash(c)
    @test am == ER.Mechanism(c, allostery([:EqualAI, :EqualAI, :OnlyA], 2, [site]))
    @test am != ER.Mechanism(c, allostery([:EqualAI, :OnlyA, :OnlyA], 2, [site]))
    @test am != ER.Mechanism(c, allostery([:EqualAI, :EqualAI, :OnlyA], 4, [site]))
    @test am != ER.Mechanism(c, allostery([:EqualAI, :EqualAI, :OnlyA], 2,
                                          EnzymeRates.Mechanisms.RegulatorySite[]))
end

@testset "a mechanism's fields are its graph" begin
    # No reaction, no cache and no fact derived from the graph: a mechanism is its
    # catalytic scheme and its allostery. The scheme holds the names of the metabolites
    # and the forms, the reactants by index and the steps with their groups; the
    # allostery the states, the oligomeric state and the sites.
    @test fieldnames(ER.Mechanism) == (:catalytic, :allosteric)
    @test fieldnames(ER.Mechanisms.CatalyticScheme) ==
          (:metabolites, :enzyme_forms, :substrates, :products, :steps)
    @test fieldnames(ER.Mechanisms.Step) ==
          (:from, :to, :taken_up, :given_off, :rapid_equilibrium, :kinetic_group)
    @test fieldnames(ER.Mechanisms.Allostery) ==
          (:allosteric_states, :oligomeric_state, :regulatory_sites)
end

@testset "construction is canonical" begin
    # Michaelis–Menten from its fields: (a) in canonical form, (b) its metabolites and
    # forms in another order, (c) every step reversed, (d) its group labels renumbered.
    # Every writing builds one mechanism, with one hash and one printed text.
    M = ER.Mechanisms
    none = Int[]
    steps = [M.Step(1, 2, [2], none, true, 1), M.Step(1, 3, [1], none, true, 2),
             M.Step(3, 2, none, none, false, 3)]
    canonical = _testhelper_plain_mechanism([:S, :P], [:E, :EP, :ES], [1], [2], steps)
    @test (canonical.catalytic.enzyme_forms, canonical.catalytic.steps) ==
          ([:E, :EP, :ES], steps)
    @test canonical == @enzyme_mechanism begin
        substrates: S
        products: P
        steps: begin
            E + S ⇌ ES
            ES <--> EP
            EP ⇌ E + P
        end
    end
    permuted = _testhelper_plain_mechanism([:P, :S], [:ES, :E, :EP], [2], [1],
                                           [M.Step(2, 3, [1], none, true, 1),
                                            M.Step(2, 1, [2], none, true, 2),
                                            M.Step(1, 3, none, none, false, 3)])
    reversed = _testhelper_plain_mechanism([:S, :P], [:E, :EP, :ES], [1], [2],
                                           [M.Step(2, 1, none, [2], true, 1),
                                            M.Step(3, 1, none, [1], true, 2),
                                            M.Step(2, 3, none, none, false, 3)])
    relabelled_steps = [M.Step(1, 2, [2], none, true, 3), M.Step(1, 3, [1], none, true, 1),
                        M.Step(3, 2, none, none, false, 2)]
    relabelled = _testhelper_plain_mechanism([:S, :P], [:E, :EP, :ES], [1], [2],
                                             relabelled_steps)
    for m in (permuted, reversed, relabelled)
        @test m == canonical && hash(m) == hash(canonical)
        @test sprint(show, m) == sprint(show, canonical)
    end
    # A step's lists in any order: A and B taken up in one step, written as [B, A].
    _testhelper_bi_uni(up) = _testhelper_plain_mechanism(
        [:A, :B, :Q], [:E, :EAB], [1, 2], [3],
        [M.Step(1, 2, up, none, true, 1), M.Step(2, 1, none, [3], false, 2)])
    unsorted = _testhelper_bi_uni([2, 1])
    @test unsorted == _testhelper_bi_uni([1, 2])
    @test [e.taken_up for e in unsorted.catalytic.steps] == [[1, 2], [3]]
    # An allosteric mechanism from its written parts takes each tag by written group label:
    # label 1 (S's binding) is canonical group 2, label 2 (the chemistry) group 3 and
    # label 3 (P's binding) group 1.
    site = M.RegulatorySite([:R], 2, [:OnlyI])
    am = ER.Mechanism(
        ER.Mechanisms.CatalyticScheme([:S, :P], [:E, :EP, :ES], [1], [2], relabelled_steps),
        ER.Mechanisms.Allostery([:EqualAI, :OnlyA, :NonequalAI], 2, [site]))
    expected = ER.Mechanism(
        canonical.catalytic,
        ER.Mechanisms.Allostery([:NonequalAI, :EqualAI, :OnlyA], 2, [site]))
    @test am == expected && hash(am) == hash(expected)
    @test sprint(show, am) == sprint(show, expected)
end

@testset "an allosteric graph sorts its sites and checks its tags and ligands" begin
    m = _testhelper_graph_uniuni()
    tags = [:EqualAI, :EqualAI, :OnlyA]
    sites = [EnzymeRates.Mechanisms.RegulatorySite([:X], 2, [:OnlyI]),
             EnzymeRates.Mechanisms.RegulatorySite([:B, :A], 2, [:OnlyA, :EqualAI])]
    am = ER.Mechanism(m.catalytic, ER.Mechanisms.Allostery(tags, 2, sites))
    @test [s.ligands for s in am.allosteric.regulatory_sites] == [[:A, :B], [:X]]
    @test am.allosteric.regulatory_sites[1].allosteric_states == [:EqualAI, :OnlyA]
    tags[1] = :OnlyA
    @test am.allosteric.allosteric_states == [:EqualAI, :EqualAI, :OnlyA]
    @test [am.allosteric.allosteric_states[g] for g in 1:3] == [:EqualAI, :EqualAI, :OnlyA]
    @test_throws(ErrorException("Mechanism: allosteric_states length 2 must " *
                                "match the number of kinetic groups 3"),
                 ER.Mechanism(m.catalytic,
                              ER.Mechanisms.Allostery([:EqualAI, :EqualAI], 2, sites)))
    @test_throws(ErrorException("Mechanism: oligomeric_state must be ≥ " *
                                "1, got 0"),
                 ER.Mechanism(m.catalytic, ER.Mechanisms.Allostery(
                     [:EqualAI, :EqualAI, :OnlyA], 0, sites)))
    # Any group may be :OnlyI, and the guard judges the mechanism: the active
    # conformation cannot switch off the rapid-equilibrium S binding (step 2), whose two
    # sides the chemistry and the P binding still join. Neither fix the message offers
    # passes: each leaves the active conformation no turnover. Untagging the binding
    # passes, but is not one of them.
    onlyi = ER.Mechanism(m.catalytic,
                         ER.Mechanisms.Allostery([:EqualAI, :OnlyI, :OnlyA], 2, sites))
    v = ER.Derivation._guard(onlyi)
    @test (v.fault, v.conformation, v.edge) == (ER.Derivation._FAST_STEP_PINNED, 1, 2)
    @test ER.Derivation._guard_message(onlyi) == _testhelper_pinned_cause(onlyi, 2, 1)
    @test all(t -> ER.Derivation._retag_passes(onlyi, [2], t), (:EqualAI, :NonequalAI))
    @test_throws(ErrorException(
        "Mechanism: ligand X appears in two distinct regulatory sites; " *
        "rendered regulator parameter names would collide. " *
        "Same-ligand-two-sites is not enumerated."),
                 ER.Mechanism(
                     m.catalytic, ER.Mechanisms.Allostery([:EqualAI, :EqualAI, :OnlyA], 2,
                     [sites; EnzymeRates.Mechanisms.RegulatorySite([:X], 4, [:OnlyA])])))
    @test_throws(ErrorException(
        "Mechanism: L is a reserved name; Keq, L, N, params, concs " *
        "and v are kept for a rate equation's constants and local names"),
        ER.Mechanism(m.catalytic, ER.Mechanisms.Allostery([:EqualAI, :EqualAI, :OnlyA], 2,
                               [EnzymeRates.Mechanisms.RegulatorySite([:L], 2, [:OnlyI])])))
    # A ligand named like a form would read as that form in the mechanism's text.
    @test_throws(ErrorException(
        "Mechanism: ligand ES is also the name of an enzyme form; a text could not tell " *
        "them apart"),
        ER.Mechanism(m.catalytic,
                     ER.Mechanisms.Allostery([:EqualAI, :EqualAI, :OnlyA], 2,
                                             [ER.Mechanisms.RegulatorySite([:ES], 2,
                                                                           [:OnlyI])])))
end

"""`raw`, a spec's raw value (`raw_mechanism`), written another way, the same mechanism:
with `order`, its forms and its steps in reverse order; with `direction`, each step from
its other end; with `labels`, its group labels numbered from the last. Its declared names
come in reverse order, the steps' lists following them. Its tags follow their labels, and
its sites and each site's ligands come in reverse order."""
function _testhelper_rewrite(raw::NamedTuple; order = true, direction = true,
                             labels = true)
    M = EnzymeRates.Mechanisms
    forms, steps, tags = raw.forms, raw.steps, raw.allosteric_states
    G, n = maximum(e -> e.kinetic_group, steps; init = 0), length(forms)
    declared = [raw.substrates; raw.products; raw.modifiers]
    redeclared = [reverse(raw.substrates); reverse(raw.products); reverse(raw.modifiers)]
    at = [findfirst(==(x), redeclared) for x in declared]
    steps = [M.Step(e.from, e.to, at[e.taken_up], at[e.given_off], e.rapid_equilibrium,
                    e.kinetic_group) for e in steps]
    direction && (steps = [M.Step(e.to, e.from, e.given_off, e.taken_up,
                                  e.rapid_equilibrium, e.kinetic_group) for e in steps])
    labels && ((steps, tags) = ([M.Step(e.from, e.to, e.taken_up, e.given_off,
                                        e.rapid_equilibrium, G + 1 - e.kinetic_group)
                                 for e in steps], reverse(tags)))
    order && ((forms, steps) = (reverse(forms),
                                reverse([M.Step(n + 1 - e.from, n + 1 - e.to, e.taken_up,
                                                e.given_off, e.rapid_equilibrium,
                                                e.kinetic_group)
                                         for e in steps])))
    (; substrates = reverse(raw.substrates), products = reverse(raw.products),
     modifiers = reverse(raw.modifiers), forms, steps, allosteric_states = tags,
     oligomeric_state = raw.oligomeric_state,
     regulatory_sites = [(; ligands = reverse(s.ligands), s.multiplicity,
                          allosteric_states = reverse(s.allosteric_states))
                         for s in reverse(raw.regulatory_sites)])
end

@testset "every writing of a mechanism gives one canonical graph" begin
    @testset "$(spec.name)" for spec in MECHANISM_TEST_SPECS
        m = _testhelper_build(spec.raw)
        for (order, direction, labels) in ((true, false, false), (false, true, false),
                                           (false, false, true), (true, true, true))
            other = _testhelper_build(
                _testhelper_rewrite(spec.raw; order, direction, labels))
            @test other == m && hash(other) == hash(m)
            @test compile_rate_equation(other) === compile_rate_equation(m)
        end
    end
end

"""Whether the constructor stores each step of the mechanism block `block` turned from the
way it is written: the stored step joining its two forms with its metabolites runs from
its written second form."""
function _testhelper_flipped(block)
    raw = EnzymeRates.DSL.raw_mechanism(block)
    c = _testhelper_scheme(raw)
    declared = [raw.substrates; raw.products; raw.modifiers]
    side(f, list) = (f, sort(list))
    stored = [(side(c.enzyme_forms[s.from], c.metabolites[s.taken_up]),
               side(c.enzyme_forms[s.to], c.metabolites[s.given_off])) for s in c.steps]
    map(raw.steps) do e
        a = side(raw.forms[e.from], declared[e.taken_up])
        b = side(raw.forms[e.to], declared[e.given_off])
        (a, b) in stored && return false
        (b, a) in stored && return true
        error("no stored step joins $a and $b")
    end
end

@testset "canonical form stores each step the way the graph gives it, a group as one" begin
    # A step with one metabolite takes it up, so the fused release FP → F + P is stored
    # as F + P → FP. The isomerizations lie on the one cycle that turns the reaction
    # over, E → ES → FS → FP → F → E, and run its way.
    @test _testhelper_flipped(:(begin
        substrates: S
        products: P
        steps: begin
            E + S <--> ES
            FS <--> ES
            FS <--> FP
            FP <--> F + P
            F <--> E
        end
    end)) == [false, true, false, true, false]
    @test _testhelper_flipped(:(begin
        substrates: S
        products: P
        steps: begin
            ES <--> E + S
            ES <--> FS
            FP <--> FS
            F + P <--> FP
            E <--> F
        end
    end)) == [true, false, true, false, true]
    # F → E runs the way the reaction turns over and leads its kinetic group, nearer free
    # E than its mirror at the I-bound forms, which keeps the direction it is written in
    # relative to F → E. Both written the other way round, the group is stored the same
    # way; written against each other, they write another mechanism. The release of P is
    # stored as the binding it reverses.
    _testhelper_mirrored(free, bound) = :(begin
        substrates: A
        products: P
        modifiers: I
        steps: begin
            E + A <--> EA
            EA <--> EP
            EP <--> F + P
            ($free, $bound)
            (E + I ⇌ EI, F + I ⇌ FI)
        end
    end)
    c = _testhelper_mirrored(:(F <--> E), :(FI <--> EI))
    d = _testhelper_mirrored(:(E <--> F), :(EI <--> FI))
    against = _testhelper_mirrored(:(F <--> E), :(EI <--> FI))
    built(block) = _testhelper_build(EnzymeRates.DSL.raw_mechanism(block))
    @test _testhelper_flipped(c) == [false, false, true, false, false, false, false]
    @test _testhelper_flipped(d) == [false, false, true, true, true, false, false]
    @test built(c) == built(d)
    @test _testhelper_flipped(against) == [false, false, true, false, false, false, false]
    @test built(against) != built(c)
end

@testset "a graph holds any number of metabolites, and a step any number of one" begin
    M = ER.Mechanisms
    names = [[Symbol("S", i) for i in 1:32]; :P]
    # Each metabolite bound to free E into a form of its own, in a group of its own.
    _testhelper_bindings(n) = _testhelper_plain_mechanism(
        names[[1:n; 33]], [:E; [Symbol("E", x) for x in names[[1:n; 33]]]],
        collect(1:n), [n + 1],
        [M.Step(1, i + 1, [i], Int[], true, i) for i in 1:(n + 1)])
    @test length(_testhelper_bindings(31).catalytic.metabolites) == 32
    @test length(_testhelper_bindings(32).catalytic.metabolites) == 33
    for n in (15, 16)
        c = M.CatalyticScheme([:A, :P], [:E, :F], [1], [2],
                              [M.Step(1, 2, fill(1, n), [2], false, 1)])
        @test count(==(1), only(c.steps).taken_up) == n
    end
end

"""The raw value of `block`, a mechanism macro's block, read by the DSL walker."""
_testhelper_walk(block) = EnzymeRates.DSL.raw_mechanism(block)

@testset "the constructors build one value from every writing of a mechanism" begin
    @testset "$(spec.name)" for spec in MECHANISM_TEST_SPECS
        x = _testhelper_build(spec.raw)
        # The mechanism holds the canonical scheme, its groups numbered by first step.
        c = _testhelper_scheme(spec.raw)
        @test (x.catalytic.metabolites, x.catalytic.enzyme_forms, x.catalytic.substrates,
               x.catalytic.products) ==
              (c.metabolites, c.enzyme_forms, c.substrates, c.products)
        unlabelled(steps) = [(e.from, e.to, e.taken_up, e.given_off, e.rapid_equilibrium)
                             for e in steps]
        @test unlabelled(x.catalytic.steps) == unlabelled(c.steps)
        # The same partition of the steps into kinetic groups, numbered by first step.
        groups(s) = Set(Set(ks) for ks in _testhelper_group_members(s))
        @test groups(x.catalytic) == groups(c)
        @test unique(e.kinetic_group for e in x.catalytic.steps) == 1:length(groups(c))
        other = _testhelper_build(_testhelper_rewrite(spec.raw))
        @test other == x && hash(other) == hash(x)
    end
    # Each tag goes to the canonical group its label became.
    am = _testhelper_build(_testhelper_walk(:(begin
        substrates: S
        products: P
        oligomeric_state: 2
        steps: begin
            E + S ⇌ ES :: EqualAI
            ES <--> EP :: OnlyA
            EP ⇌ E + P :: NonequalAI
        end
    end)))
    @test am.allosteric.allosteric_states == [:NonequalAI, :EqualAI, :OnlyA]
end

"""Whether `v` takes the memory of its copy, which `copy` stores at its length."""
_testhelper_at_length(v) = Base.summarysize(v) == Base.summarysize(copy(v))

"""Whether `m`, a catalytic scheme or a mechanism, stores each of its vectors at its
length: `m` takes the memory of the value rebuilt from a copy of each vector, its steps'
lists and its sites included, which `copy` stores at its length."""
_testhelper_stored_at_length(m::ER.Mechanisms.CatalyticScheme) =
    Base.summarysize(m) == Base.summarysize(
        ER.Mechanisms._unchecked_scheme(copy(m.metabolites), copy(m.enzyme_forms),
            copy(m.substrates), copy(m.products),
            [ER.Mechanisms.Step(e.from, e.to, copy(e.taken_up), copy(e.given_off),
                                e.rapid_equilibrium, e.kinetic_group) for e in m.steps]))
_testhelper_stored_at_length(am::ER.Mechanism) =
    _testhelper_stored_at_length(am.catalytic) && Base.summarysize(am) ==
    Base.summarysize(ER.Mechanism(am.catalytic,
        ER.Mechanisms.Allostery(copy(am.allosteric.allosteric_states),
                                am.allosteric.oligomeric_state,
                                [ER.Mechanisms.RegulatorySite(copy(s.ligands),
                                                              s.multiplicity,
                                                              copy(s.allosteric_states))
                                 for s in am.allosteric.regulatory_sites])))

@testset "a mechanism stores each vector at its length" begin
    # A search keeps every mechanism it emits, so memory a vector holds past its length
    # is paid once per mechanism: the ter-ter init set holds 250,855 of them.
    @testset "$(spec.name)" for spec in MECHANISM_TEST_SPECS
        @test _testhelper_stored_at_length(_testhelper_build(spec.raw))
    end
    rxn = @enzyme_reaction begin
        substrates: A[C], B[N]
        products: P[C], Q[N]
        competitive_inhibitors: I
    end
    init = _testhelper_mechanisms(EnzymeRates.Enumeration.init_mechanisms(rxn))
    children = _testhelper_children(rxn, init[1:10])
    @test any(c -> _testhelper_is_allosteric(c), children)
    @test all(_testhelper_stored_at_length, init)
    @test all(_testhelper_stored_at_length, children)
end

@testset "the constructors reject what no graph represents, each with its message" begin
    # Uni-uni over S and P, E + S ⇌ ES, ES <--> EP, EP ⇌ E + P, forms E, ES and EP (1 to
    # 3), with the forms `extra` (4, 5, …) and `steps` added.
    M = ER.Mechanisms
    none = Int[]
    uniuni = [M.Step(1, 2, [1], none, true, 1), M.Step(2, 3, none, none, false, 2),
              M.Step(3, 1, none, [2], true, 3)]
    _testhelper_uniuni(extra = Symbol[], steps...; metabolites = [:S, :P],
                       substrates = [1], products = [2], forms = [:E, :ES, :EP],
                       group = [1, 2, 3]) =
        _testhelper_plain_mechanism(metabolites, [forms; extra], substrates, products,
                     [[M.Step(e.from, e.to, e.taken_up, e.given_off, e.rapid_equilibrium, g)
                       for (e, g) in zip(uniuni, group)]; steps...])
    # Uni-uni with the modifier `x`, metabolite 3, bound to free E by a step of its own.
    _testhelper_binding(x; kw...) =
        _testhelper_uniuni([Symbol(:E, x)], M.Step(1, 4, [3], none, true, 4);
                           metabolites = [:S, :P, x], kw...)
    _testhelper_allosteric(tags, n, sites) =
        ER.Mechanism(ER.Mechanisms.CatalyticScheme([:S, :P], [:E, :ES, :EP], [1], [2],
                                                   uniuni),
                     ER.Mechanisms.Allostery(tags, n, sites))
    tags = [:EqualAI, :OnlyA, :EqualAI]
    site = M.RegulatorySite
    @test _testhelper_uniuni() == _testhelper_graph_uniuni()
    # Regulatory sites.
    @test_throws(ErrorException("RegulatorySite: needs at least one ligand"),
                 _testhelper_allosteric(tags, 2, [site(Symbol[], 2, Symbol[])]))
    @test_throws(ErrorException("RegulatorySite: length(ligands)=1 must equal " *
                                "length(allosteric_states)=0"),
                 _testhelper_allosteric(tags, 2, [site([:R], 2, Symbol[])]))
    @test_throws(ErrorException("RegulatorySite: multiplicity must be ≥ 1, got 0"),
                 _testhelper_allosteric(tags, 2, [site([:R], 0, [:OnlyA])]))
    @test_throws(ErrorException("RegulatorySite: allosteric state Bogus must be one of " *
                                "(:OnlyA, :OnlyI, :EqualAI, :NonequalAI)"),
                 _testhelper_allosteric(tags, 2, [site([:R], 2, [:Bogus])]))
    @test_throws(ErrorException("RegulatorySite: ligand R appears twice in one " *
                                "regulatory site"),
                 _testhelper_allosteric(tags, 2, [site([:R, :X, :R], 2,
                                                       [:OnlyA, :OnlyI, :OnlyI])]))
    @test_throws(ErrorException(
        "Mechanism: ligand R appears in two distinct regulatory sites; " *
        "rendered regulator parameter names would collide. " *
        "Same-ligand-two-sites is not enumerated."),
        _testhelper_allosteric(tags, 2, [site([:R], 2, [:OnlyA]), site([:R], 4, [:OnlyI])]))
    # A ligand takes no name a metabolite cannot take, and no form's name.
    @test_throws(ErrorException("Mechanism: \"1R\" cannot name a metabolite; a " *
                                "metabolite name is an identifier that is no keyword"),
                 _testhelper_allosteric(tags, 2, [site([Symbol("1R")], 2, [:OnlyA])]))
    @test_throws(ErrorException("Mechanism: E cannot name a metabolite; it " *
                                "labels the free enzyme form"),
                 _testhelper_allosteric(tags, 2, [site([:E], 2, [:OnlyA])]))
    @test_throws(ErrorException("Mechanism: ligand ES is also the name of an enzyme " *
                                "form; a text could not tell them apart"),
                 _testhelper_allosteric(tags, 2, [site([:ES], 2, [:OnlyA])]))
    # Tags, by written group label: label 2 is canonical group 3.
    @test_throws(ErrorException("Mechanism: allosteric_states length 2 must " *
                                "match the number of kinetic groups 3"),
                 _testhelper_allosteric(tags[1:2], 2, M.RegulatorySite[]))
    @test_throws(ErrorException("Mechanism: oligomeric_state must be ≥ " *
                                "1, got 0"),
                 _testhelper_allosteric(tags, 0, M.RegulatorySite[]))
    # A Bool is no multiplicity: each constructor takes an Int and converts no Bool.
    @test_throws MethodError site([:R], true, [:OnlyA])
    @test_throws MethodError _testhelper_allosteric(tags, true, M.RegulatorySite[])
    @test_throws MethodError ER.Mechanism(_testhelper_uniuni().catalytic,
                                          ER.Mechanisms.Allostery(tags, true,
                                                                  M.RegulatorySite[]))
    # Any group may be :OnlyI: the constructor builds the chemistry :OnlyI, and the guard
    # rejects it, since the active conformation then turns nothing over.
    onlyi = _testhelper_allosteric([:EqualAI, :OnlyI, :EqualAI], 2, M.RegulatorySite[])
    v = ER.Derivation._guard(onlyi)
    @test (v.fault, v.conformation) == (ER.Derivation._NO_ACTIVE_TURNOVER, 1)
    @test startswith(ER.Derivation._guard_message(onlyi),
                     "Mechanism: without its :OnlyI steps ES <--> EP the " *
                     "active conformation turns no substrate into product")
    # Declarations and group labels.
    @test_throws(ErrorException("Mechanism: substrates must not be empty"),
                 _testhelper_uniuni(substrates = Int[]))
    @test_throws(ErrorException("Mechanism: products must not be empty"),
                 _testhelper_uniuni(products = Int[]))
    @test_throws(ErrorException("Mechanism: duplicate substrate names"),
                 _testhelper_uniuni(metabolites = [:S, :P, :S], substrates = [1, 3]))
    @test_throws(ErrorException("Mechanism: duplicate product names"),
                 _testhelper_uniuni(metabolites = [:S, :P, :P], products = [2, 3]))
    @test_throws(ErrorException("Mechanism: S named as both a substrate and a product; " *
                                "concentrations and constants are keyed by name"),
                 _testhelper_uniuni(products = [2, 1]))
    @test_throws(ErrorException("Mechanism: duplicate modifier names"),
                 _testhelper_uniuni(metabolites = [:S, :P, :I, :I]))
    # A reactant that binds as a modifier is taken up as itself, so it is not declared
    # again as a modifier.
    @test_throws(ErrorException("Mechanism: S is declared both as a reactant and as a " *
                                "modifier; a step that binds a reactant as a modifier " *
                                "simply takes it up, as `ES + S ⇌ ESS` does"),
                 _testhelper_uniuni(metabolites = [:S, :P, :S]))
    @test_throws(ErrorException("Mechanism: each step takes one kinetic group label " *
                                "from 1 to the number of groups"),
                 _testhelper_uniuni(group = [0, 1, 2]))
    @test_throws(ErrorException("Mechanism: kinetic group 2 is empty; a kinetic group " *
                                "holds at least one step"),
                 _testhelper_uniuni(group = [1, 3, 3]))
    # A step's ends index the forms, and its lists the metabolites.
    for (from, to) in ((2, 4), (0, 1))
        @test_throws(ErrorException("Mechanism: step 4 joins forms $from and $to, but " *
                                    "there are 3 forms; a step's ends index the enzyme " *
                                    "forms"),
                     _testhelper_uniuni(Symbol[], M.Step(from, to, none, none, false, 4)))
    end
    @test_throws(ErrorException("Mechanism: step 4 takes up [3] and gives off Int64[], " *
                                "but there are 2 metabolites; a step's lists index the " *
                                "metabolites"),
                 _testhelper_uniuni([:EI], M.Step(1, 4, [3], none, true, 4)))
    # Names a metabolite cannot take: no identifier, a keyword, E, a reserved name, a
    # name holding `_`.
    _testhelper_unwritable(x) = ErrorException(
        "Mechanism: \"$x\" cannot name a metabolite; a metabolite name is an identifier " *
        "that is no keyword")
    @test_throws _testhelper_unwritable("1I") _testhelper_binding(Symbol("1I"))
    @test_throws(_testhelper_unwritable("end"),
                 _testhelper_uniuni(metabolites = [:S, :P, :end], products = [2, 3]))
    @test_throws(ErrorException("Mechanism: E cannot name a metabolite; it labels the " *
                                "free enzyme form"),
                 _testhelper_uniuni(metabolites = [:S, :P, :E], substrates = [1, 3]))
    @test_throws(ErrorException("Mechanism: Keq is a reserved name; Keq, L, N, " *
                                "params, concs and v are kept for a rate equation's " *
                                "constants and local names"),
                 _testhelper_uniuni(metabolites = [:S, :P, :Keq], substrates = [1, 3]))
    @test_throws(ErrorException("Mechanism: the metabolite name S_1 holds `_`; constant " *
                                "names join metabolite names with `_`, so a metabolite " *
                                "name holds none"),
                 _testhelper_binding(:S_1))
    # A name Julia reads as another (a micro sign as μ, a decomposed é as one character)
    # prints as text that builds another mechanism, so no DSL text can write it. The
    # strings are written with escapes: Julia would read the characters themselves as
    # the other names.
    _testhelper_read_as(who, x, read) = ErrorException(
        "$who: Julia reads the name \"$x\" as \"$read\", so no DSL text can write it; " *
        "use \"$read\"")
    @test_throws(_testhelper_read_as("Mechanism", "µM", "μM"),
                 _testhelper_binding(Symbol("µM")))
    @test_throws(_testhelper_read_as("Mechanism", "é", "é"),
                 _testhelper_uniuni(metabolites = [:S, :P, Symbol("é")],
                                    substrates = [1, 3]))
    @test_throws(_testhelper_read_as("Mechanism", "ɛ", "ε"),
                 _testhelper_allosteric(tags, 2, [site([Symbol("ɛ")], 2, [:OnlyA])]))
    # Form names: no identifier, read as another, given twice, a metabolite's name.
    @test_throws(ErrorException("Mechanism: \"1E\" cannot name an enzyme form; an " *
                                "enzyme form name is an identifier that is no keyword"),
                 _testhelper_uniuni([Symbol("1E")], M.Step(1, 4, none, none, false, 4)))
    @test_throws(_testhelper_read_as("Mechanism", "Eµ", "Eμ"),
                 _testhelper_uniuni([Symbol("Eµ")], M.Step(1, 4, none, none, false, 4)))
    @test_throws(ErrorException("Mechanism: two enzyme forms are named ES"),
                 _testhelper_uniuni([:ES], M.Step(1, 4, none, none, false, 4)))
    @test_throws(ErrorException("Mechanism: S names both an enzyme form and a " *
                                "metabolite; a text could not tell them apart"),
                 _testhelper_uniuni([:S], M.Step(1, 4, none, none, false, 4)))
    # Steps: a self-loop, a metabolite taken up and given off.
    @test_throws(ErrorException("Step: both ends are the form ES; a step joins two forms"),
                 _testhelper_uniuni(Symbol[], M.Step(2, 2, none, none, false, 4)))
    @test_throws(ErrorException("Step: both ends are the form EA; a step joins two forms"),
                 _testhelper_plain_mechanism(
                     [:A, :P], [:E, :EA], [1], [2],
                     [M.Step(1, 2, [1], none, true, 1), M.Step(2, 1, none, [2], false, 2),
                      M.Step(2, 2, none, none, false, 3)]))
    @test_throws(ErrorException("Step EA → EQ: B both consumed and released"),
                 _testhelper_build(_testhelper_walk(:(begin
                     substrates: A, B
                     products: Q
                     steps: begin
                         E + A ⇌ EA
                         EA + B <--> EQ + B
                         EQ ⇌ E + Q
                     end
                 end))))
    # One copy of B taken up and two given off: the sides share B though their counts
    # differ. The step gives off more substrate than it takes up, so it is stored the
    # other way round.
    @test_throws(ErrorException("Step EQ → EA: B both consumed and released"),
                 _testhelper_build(_testhelper_walk(:(begin
                     substrates: A, B
                     products: Q
                     steps: begin
                         E + A ⇌ EA
                         EA + B <--> EQ + B + B
                         EQ ⇌ E + Q
                     end
                 end))))
    # A kinetic group of steps of different kinds or flags.
    for (steps, text) in (
            (:(begin
                 E + A <--> EA
                 (EA + B <--> EAB, EA + B <--> EQ + P)
                 EQ <--> E + Q
             end), "EA_B → EAB (SS) and EA_B → EQ_P (SS)"),
            (:(begin
                 (E + A ⇌ EA, EQ + A <--> EAQ)
                 EA <--> EQ + P
                 EQ <--> E + Q
             end), "E_A → EA (RE) and EQ_A → EAQ (SS)"))
        @test_throws(ErrorException(
            "Mechanism: a kinetic group holds $text, which differ in kind or in RE/SS; " *
            "the steps of a kinetic group share their constants, so they must be the " *
            "same kind of step with the same flag"),
            _testhelper_build(_testhelper_walk(:(begin
                substrates: A, B
                products: P, Q
                steps: $steps
            end))))
    end
    # One reaction twice: in two groups, with different flags, stored once each way, and in
    # one group.
    _testhelper_twice(steps) = _testhelper_build(_testhelper_walk(:(begin
        substrates: S
        products: P
        steps: $steps
    end)))
    @test_throws(ErrorException(
        "Mechanism: kinetic groups 2 and 3 both hold the reaction E_S ⇌ ES (SS in group " *
        "2, RE in group 3); a reaction belongs to one kinetic group, whose constants are " *
        "named after it"),
        _testhelper_twice(:(begin
            E + S <--> ES
            E + S ⇌ ES
            ES ⇌ EP
            EP ⇌ E + P
        end)))
    @test_throws(ErrorException(
        "Mechanism: kinetic groups 3 and 4 both hold the reaction ES ⇌ EP; a reaction " *
        "belongs to one kinetic group, whose constants are named after it"),
        _testhelper_twice(:(begin
            E + S <--> ES
            ES <--> EP
            EP <--> ES
            EP <--> E + P
        end)))
    # The tied pair and the lone step both store F → E, which runs the way the reaction
    # turns over, in groups 5 and 6.
    @test_throws(ErrorException(
        "Mechanism: kinetic groups 5 and 6 both hold the reaction F ⇌ E; a reaction " *
        "belongs to one kinetic group, whose constants are named after it"),
        _testhelper_build(_testhelper_walk(:(begin
            substrates: A
            products: P
            modifiers: I
            steps: begin
                E + A <--> EA
                EA <--> EP
                EP <--> F + P
                F <--> E
                (E <--> F, EI <--> FI)
                E + I ⇌ EI
                F + I ⇌ FI
                EI + P <--> EIP
                FI + A <--> FAI
            end
        end))))
    @test_throws(ErrorException(
        "Mechanism: kinetic group 2 holds the reaction E_S ⇌ ES twice; a mechanism holds " *
        "each reaction once"),
        _testhelper_twice(:(begin
            (E + S <--> ES, ES <--> E + S)
            ES <--> EP
            EP <--> E + P
        end)))
    # No form named E: no step at all, a free form named F, forms whose names start with
    # E but are not E, or a ping-pong whose free form is F and whose covalent forms'
    # names start with E.
    no_E = ArgumentError("Mechanism: no enzyme form is named `E`; the free enzyme is the " *
                         "form named `E`")
    @test_throws no_E _testhelper_plain_mechanism([:S, :P], Symbol[], [1], [2], M.Step[])
    @test_throws no_E _testhelper_build(_testhelper_walk(:(begin
        substrates: S
        products: P
        steps: begin
            F + S ⇌ FS
            FS <--> FP
            FP ⇌ F + P
        end
    end)))
    @test_throws no_E _testhelper_build(_testhelper_walk(:(begin
        substrates: S
        products: P
        steps: begin
            F + S ⇌ ES
            ES <--> EP
            EP ⇌ F + P
        end
    end)))
    @test_throws no_E _testhelper_build(_testhelper_walk(:(begin
        substrates: A, B
        products: P, Q
        steps: begin
            F + A ⇌ FA
            FA <--> E_res_plusA_minusP + P
            E_res_plusA_minusP + B ⇌ EB_res_plusA_minusP
            EB_res_plusA_minusP <--> FQ
            FQ ⇌ F + Q
        end
    end)))
    # Forms that no steps join to free E: an isomerization between unreached forms.
    @test_throws(ErrorException("Mechanism: the steps do not form one connected graph; " *
                                "no steps lead from free E to EstarP, EstarS"),
                 _testhelper_build(_testhelper_walk(:(begin
        substrates: S
        products: P
        steps: begin
            E + S ⇌ ES
            ES <--> EP
            EP ⇌ E + P
            EstarS <--> EstarP
        end
    end))))
    # Declared names no step exchanges: a substrate, a product, a modifier.
    @test_throws(ErrorException("Mechanism: substrate T is declared, but no step uses it"),
                 _testhelper_uniuni(metabolites = [:S, :P, :T], substrates = [1, 3]))
    @test_throws(ErrorException("Mechanism: product Q is declared, but no step uses it"),
                 _testhelper_uniuni(metabolites = [:S, :P, :Q], products = [2, 3]))
    @test_throws(ErrorException("Mechanism: modifier I is declared, but no step binds it"),
                 _testhelper_uniuni(metabolites = [:S, :P, :I]))
end

@testset "a graph written back as its forms and steps builds the same graph" begin
    @testset "$(spec.name)" for spec in MECHANISM_TEST_SPECS
        x = _testhelper_build(spec.raw)
        raw = _testhelper_written(x)
        @test _testhelper_build(raw) == x
        # Its fields are stored as given: in canonical order, each step in its stored
        # direction, labelled by group.
        c = _testhelper_scheme(raw)
        @test c == x.catalytic
    end
end

# ── The canonical form pinned as text, and the sub-graph reads against oracles ──

"""The regression census's six reactions: uni-uni, bi-bi, ping-pong, two sites, ter-bi and
the reaction whose names hold one another as prefixes."""
function _testhelper_census_reactions()
    [
        @enzyme_reaction(begin
            substrates: S[C]
            products: P[C]
            competitive_inhibitors: I
            allosteric_regulators: R
            oligomeric_state: 2
        end),
        @enzyme_reaction(begin
            substrates: A[C], B[N]
            products: P[C], Q[N]
            competitive_inhibitors: A, I
        end),
        @enzyme_reaction(begin
            substrates: A[CX], B[N]
            products: P[C], Q[NX]
            allosteric_regulators: R
            oligomeric_state: 2
            shared_catalytic_site: (A, P)
        end),
        @enzyme_reaction(begin
            substrates: A[C6H12O6]
            products: B[C6H12O6]
            allosteric_regulators: X::Activator, Y::Inhibitor
            oligomeric_state: 2
        end),
        @enzyme_reaction(begin
            substrates: A[C], B[N], D[X]
            products: P[CN], Q[X]
        end),
        @enzyme_reaction(begin
            substrates: G[CX], GT[N]
            products: Gl[C], U[NX]
            competitive_inhibitors: GTI
            shared_catalytic_site: (G, Gl)
        end)]
end

"""The mechanisms the regression census enumerates first: the init populations of its six
reactions, then the children of every init mechanism of its bi-bi reaction."""
function _testhelper_census_population()
    reactions = _testhelper_census_reactions()
    population = EnzymeRates.Mechanism[]
    foreach(r -> append!(population, _testhelper_mechanisms(
        EnzymeRates.Enumeration.init_mechanisms(r))), reactions)
    bibi = reactions[2]
    append!(population,
            _testhelper_children(bibi, _testhelper_mechanisms(
                EnzymeRates.Enumeration.init_mechanisms(bibi))))
end

"""The census's first allosteric mechanisms: the allosteric children of the init
mechanisms of its three reactions that declare an allosteric regulator."""
function _testhelper_census_allosteric()
    children = [_testhelper_children(
        r, _testhelper_mechanisms(EnzymeRates.Enumeration.init_mechanisms(r)))
                for r in _testhelper_census_reactions()[[1, 3, 4]]]
    filter(x -> _testhelper_is_allosteric(x), reduce(vcat, children))
end

"""Graph `x`'s canonical form as text lines: its names, its forms in order, each edge in
its stored direction (`_testhelper_edge_text`) with its group (and the group's tag), then
any oligomeric state and sites."""
function _testhelper_canonical_text(x)
    m, allo = x.catalytic, _testhelper_is_allosteric(x)
    g = _testhelper_graph(x)
    lines = ["substrates: $(join(g.substrates, ", ")) | products: " *
             "$(join(g.products, ", "))" *
             (isempty(g.modifiers) ? "" : " | modifiers: $(join(g.modifiers, ", "))"),
             "forms: " * join(m.enzyme_forms, ", ")]
    for k in eachindex(m.steps)
        g = m.steps[k].kinetic_group
        push!(lines, "$k: $(_testhelper_edge_text(m, k)) | group $g" *
                     (allo ? " :: $(x.allosteric.allosteric_states[g])" : ""))
    end
    allo && push!(lines, "oligomeric_state: $(x.allosteric.oligomeric_state)")
    for s in x.allosteric.regulatory_sites
        push!(lines, "regulatory_site($(s.multiplicity)): " *
                     join(["$l::$t" for (l, t) in zip(s.ligands, s.allosteric_states)], ", "))
    end
    lines
end

"""
The inactive edges of the allosteric mechanism whose test reader (`_testhelper_graph`) is
`g`, read off `g`'s steps alone, one flag per step: no step of an `:OnlyA` group, and every
other step whose two forms the steps outside `:OnlyA` groups reach from free E.
"""
function _testhelper_inactive_reach(g)
    outside = [g.tags[s.group] !== :OnlyA for s in g.steps]
    reach = falses(length(g.forms))
    reach[1] = true
    grew = true
    while grew
        grew = false
        for (j, s) in enumerate(g.steps)
            outside[j] && reach[s.from] != reach[s.to] || continue
            reach[s.from] = reach[s.to] = true
            grew = true
        end
    end
    BitVector([outside[j] && reach[s.from] && reach[s.to] for (j, s) in enumerate(g.steps)])
end

"""
Whether the sub-graph reads of `x` (`rapid_equilibrium_segments`,
`steps_in_inactive_conformation`) give what the oracles read off its test reader's steps
give (`_testhelper_re_segments`, `_testhelper_inactive_reach`): the catalytic graph's
segments, form by form and breadth-first, with each form's positive offsets, and, for an
allosteric mechanism, the inactive conformation's edges and their segments.
"""
function _testhelper_subgraphs_agree(x)
    m, g = x.catalytic, _testhelper_graph(x)
    # Step `j` of the reader is edge `flat[j]`.
    flat = ER.Derivation._flat_edges(m)
    function _testhelper_segments_agree(mask, kept)
        s = ER.Mechanisms.rapid_equilibrium_segments(m, mask)
        segments, extras = _testhelper_re_segments(g, kept)
        present = findall(!=(0), s.segment_of)
        names = m.metabolites
        _testhelper_offsets(f) = Dict(names[n] => s.offsets[n, f]
                                      for n in eachindex(names) if s.offsets[n, f] > 0)
        s.segments == segments && [_testhelper_offsets(f) for f in present] == extras
    end
    _testhelper_segments_agree(trues(length(m.steps)), trues(length(flat))) || return false
    _testhelper_is_allosteric(x) || return true
    inactive = ER.Mechanisms.steps_in_inactive_conformation(x)
    kept = _testhelper_inactive_reach(g)
    inactive[flat] == kept && _testhelper_segments_agree(inactive, kept)
end

# The census's first populations (`_testhelper_census_population`,
# `_testhelper_census_allosteric`), enumerated once for the three testsets that read them
# and released after them.
let population = _testhelper_census_population(),
    allosteric = _testhelper_census_allosteric()
@testset "a graph's canonical form, pinned as text, read back from its Sig and raw" begin
    @testset "$(spec.name)" for spec in MECHANISM_TEST_SPECS
        @test spec.mechanism == _testhelper_build(spec.raw)
    end
    # The specs' canonical text, a blank line between specs, pinned as a fixture.
    text = reduce((a, b) -> [a; ""; b],
                  [["# $(spec.name)";
                    _testhelper_canonical_text(_testhelper_build(spec.raw))]
                   for spec in MECHANISM_TEST_SPECS])
    fixture = joinpath(@__DIR__, "..", "fixtures", "canonical_text_golden.txt")
    if !isfile(fixture)
        mkpath(dirname(fixture)); write(fixture, join(text, "\n") * "\n")
        @warn "bootstrapped canonical text golden fixture; commit it" fixture
    end
    @test text == readlines(fixture)
    # Free E comes first, also where another label sorts before E.
    raw = _testhelper_walk(:(begin
        substrates: A
        products: P
        steps: begin
            E + A <--> EA
            EA <--> EP
            EP <--> D + P
            D <--> E
        end
    end))
    x = _testhelper_build(raw)
    @test x.catalytic.enzyme_forms[1:2] == [:E, :D]
    @test _testhelper_graph(x).forms[1:2] == ["E", "D"]
    # The census's first populations: each graph decodes from its Sig and builds again
    # from the forms and steps it writes back. The count: the six reactions' init
    # mechanisms, then the bi-bi reaction's children.
    @test length(population) == 5508 + 5471
    @test all(_testhelper_decode(ER.Derivation._sig_of(x)) == x for x in population)
    @test all(_testhelper_build(_testhelper_written(x)) == x for x in population)
end

@testset "a sub-graph's segments and inactive edges are what its steps give" begin
    @testset "$(spec.name)" for spec in MECHANISM_TEST_SPECS
        @test _testhelper_subgraphs_agree(spec.mechanism)
    end
    @test isempty([i for (i, x) in enumerate(population)
                   if !_testhelper_subgraphs_agree(x)])
end

@testset "a graph's DSL text reads back as the graph" begin
    # `show` prints the macro text that builds a mechanism and `DSL.parse_mechanism`
    # reads it back, for the specs, plain and allosteric, and the census's first
    # populations, plain and allosteric.
    @testset "$(spec.name)" for spec in MECHANISM_TEST_SPECS
        @test _testhelper_reads_back(spec.mechanism)
        @test ER.DSL.parse_mechanism(sprint(show, spec.mechanism)) == spec.mechanism
    end
    # The allosteric children of the uni-uni, ping-pong and two-site reactions.
    @test length(allosteric) == 5 + 885 + 5
    @test isempty([i for (i, x) in enumerate(population) if !_testhelper_reads_back(x)])
    @test isempty([i for (i, x) in enumerate(allosteric) if !_testhelper_reads_back(x)])
end
end

@testset "cycles, distances and a group's lead orient steps alike in every writing" begin
    # Isomerizations a turnover cycle orients, exchanges that run away from free E, and an
    # isomerization that keeps its direction relative to its group's lead, each case in two
    # writings, against each case's pinned canonical text. In the tied fixture neither
    # exchange lies on a cycle: each runs from its form nearer free E, EI and FI,
    # taking up J.
    ss_text = ["substrates: S | products: P",
               "forms: E, ES, F, FP, FS",
               "1: E + S <--> ES | group 1",
               "2: ES <--> FS | group 2",
               "3: F <--> E | group 3",
               "4: F + P <--> FP | group 4",
               "5: FS <--> FP | group 5"]
    tied_text = ["substrates: A | products: P | modifiers: I, J",
                 "forms: E, EA, EAI, EI, EJ, EJP, EP, F, FI, FJ",
                 "1: E + A <--> EA | group 1",
                 "2: E + I ⇌ EI | group 2",
                 "3: EA <--> EP | group 3",
                 "4: EI + A ⇌ EAI | group 4",
                 "5: EI + J <--> EJ + I | group 5",
                 "6: EJ + P ⇌ EJP | group 6",
                 "7: F <--> E | group 7",
                 "8: F + P <--> EP | group 8",
                 "9: F + I ⇌ FI | group 2",
                 "10: FI <--> EI | group 7",
                 "11: FI + J <--> FJ + I | group 5"]
    allosteric_text = ["substrates: A | products: P | modifiers: I, J",
                       "forms: E, EA, EI, EJ, EP, F, FI, FJ",
                       "1: E + A <--> EA | group 1 :: EqualAI",
                       "2: E + I ⇌ EI | group 2 :: EqualAI",
                       "3: EA <--> EP | group 3 :: EqualAI",
                       "4: EI + J <--> EJ + I | group 4 :: NonequalAI",
                       "5: F ⇌ E | group 5 :: EqualAI",
                       "6: F + P <--> EP | group 6 :: EqualAI",
                       "7: F + I ⇌ FI | group 2 :: EqualAI",
                       "8: FI + J <--> FJ + I | group 4 :: NonequalAI",
                       "oligomeric_state: 1"]
    _testhelper_tied_steps(free, bound, exchange) = :(begin
        E + A <--> EA
        EA <--> EP
        EP <--> F + P
        ($free, $bound)
        (E + I ⇌ EI, F + I ⇌ FI)
        EI + A ⇌ EAI
        EJ + P ⇌ EJP
        $exchange
    end)
    blocks = [
        :(begin
            substrates: S
            products: P
            steps: begin
                E + S <--> ES
                FS <--> ES
                FS <--> FP
                FP <--> F + P
                F <--> E
            end
        end),
        :(begin
            substrates: S
            products: P
            steps: begin
                ES <--> FS
                E + S <--> ES
                F + P <--> FP
                FP <--> FS
                E <--> F
            end
        end)]
    for (free, bound, exchange) in (
            (:(F <--> E), :(FI <--> EI),
             :((EI + J <--> EJ + I, FI + J <--> FJ + I))),
            (:(E <--> F), :(EI <--> FI),
             :((FJ + I <--> FI + J, EJ + I <--> EI + J))))
        push!(blocks, :(begin
            substrates: A
            products: P
            modifiers: I, J
            steps: $(_testhelper_tied_steps(free, bound, exchange))
        end))
    end
    texts = [ss_text, ss_text, tied_text, tied_text]
    for (raw, text) in zip(_testhelper_walk.(blocks), texts)
        @test _testhelper_canonical_text(_testhelper_build(raw)) == text
    end
    for exchange in (:((EI + J <--> EJ + I, FI + J <--> FJ + I)),
                     :((FJ + I <--> FI + J, EJ + I <--> EI + J)))
        raw = _testhelper_walk(:(begin
            substrates: A
            products: P
            modifiers: I, J
            steps: begin
                E + A <--> EA                  :: EqualAI
                EA <--> EP                   :: EqualAI
                EP <--> F + P                  :: EqualAI
                F ⇌ E                            :: EqualAI
                (E + I ⇌ EI, F + I ⇌ FI)     :: EqualAI
                $exchange                        :: NonequalAI
            end
        end))
        @test _testhelper_canonical_text(_testhelper_build(raw)) == allosteric_text
    end
end

"""Whether `x` survives its round trips: its printed text reads back as `x` and prints
again as that text, two reads of the text hash alike, and `x` decodes from the Sig of its
compiled law as itself."""
function _testhelper_round_trips(x)
    text = sprint(show, x)
    a, b = ER.DSL.parse_mechanism(text), ER.DSL.parse_mechanism(text)
    a == x && sprint(show, a) == text && hash(a) == hash(b) == hash(x) &&
        _testhelper_decode(ER.Derivation._sig_of(x)) == x
end

@testset "storage round trips" begin
    # The specs, plain and allosteric, and the bi-bi init set, the enumeration golden set
    # (`phase2_init_golden.txt`).
    rxn = @enzyme_reaction begin
        substrates: A[C], B[N]
        products: P[C], Q[N]
    end
    @testset "$(spec.name)" for spec in MECHANISM_TEST_SPECS
        @test _testhelper_round_trips(spec.mechanism)
    end
    init = _testhelper_mechanisms(EnzymeRates.Enumeration.init_mechanisms(rxn))
    @test length(init) == 239
    @test isempty([i for (i, x) in enumerate(init) if !_testhelper_round_trips(x)])
end

"""Whether the steps each conformation of `am` runs (`steps_in_active_conformation`,
`steps_in_inactive_conformation`) are those its tags do not exclude
(`_testhelper_excludes`), on the forms they join to free E (`_testhelper_exclusion`). True
for a plain mechanism, whose one conformation runs every step."""
function _testhelper_conformations_hold(am::ER.Mechanism)
    cm = am.catalytic
    M = ER.Mechanisms
    _testhelper_is_allosteric(am) ||
        return M.steps_in_active_conformation(am) == trues(length(cm.steps))
    all(((1, M.steps_in_active_conformation(am)),
         (2, M.steps_in_inactive_conformation(am)))) do (c, runs)
        excluded = BitVector([_testhelper_excludes(am, k, c) for k in eachindex(cm.steps)])
        component = _testhelper_exclusion(am, c).component
        runs == BitVector([!excluded[k] && component[e.from] == 1
                           for (k, e) in enumerate(cm.steps)])
    end
end

@testset "the steps a conformation runs are what its states give" begin
    # Every construction path: the constructors (the specs' writings, the bi-bi seeds)
    # and the Sig decoder (the specs' compiled laws).
    @testset "$(spec.name)" for spec in MECHANISM_TEST_SPECS
        @test _testhelper_conformations_hold(_testhelper_build(spec.raw))
        @test _testhelper_conformations_hold(
            _testhelper_decode(ER.Derivation._sig_of(spec.mechanism)))
    end
    rxn = @enzyme_reaction begin
        substrates: A[C], B[N]
        products: P[C], Q[N]
    end
    @test all(_testhelper_conformations_hold,
              _testhelper_mechanisms(EnzymeRates.Enumeration.init_mechanisms(rxn)))
end

@testset "forms and steps are the storage" begin
    m = _testhelper_build(_testhelper_walk(:(begin
        substrates: S
        products: P
        steps: begin
            E + S ⇌ ES
            ES <--> EP
            EP ⇌ E + P
        end
    end)))
    c = m.catalytic
    @test (c.metabolites, c.enzyme_forms, c.substrates, c.products) ==
          ([:S, :P], [:E, :EP, :ES], [1], [2])
    # The release of P is stored as the binding it reverses, and sorts first.
    @test [e.taken_up for e in c.steps[1:2]] == [[2], [1]]
    @test (c.steps[2].from, c.steps[2].to, c.steps[2].given_off) == (1, 3, Int[])
    @test c.steps[1].rapid_equilibrium && c.steps[2].rapid_equilibrium &&
          !c.steps[3].rapid_equilibrium
    @test [s.kinetic_group for s in c.steps] == [1, 2, 3]
end

@testset "the constructor builds a spec's graph in under 30 µs, the median over specs" begin
    # Each spec's time per call is the best of five batches; the gate is on their median.
    # The spec's goal is about 2 µs and the median is about 5 µs locally; shared CI
    # runners take about 5× as long (23 µs on Ubuntu), so the gate is 30 µs.
    per_call = map(MECHANISM_TEST_SPECS) do spec
        _testhelper_build(spec.raw)
        minimum(@elapsed(for _ in 1:200
                             _testhelper_build(spec.raw)
                         end) / 200 for _ in 1:5)
    end
    @test median(per_call) < 30e-6
end

"""The bytes one construction of `x` allocates, from its own canonical fields, after a
warm-up call: `Mechanism(CatalyticScheme(metabolites, enzyme_forms, substrates, products,
steps), allostery)` with the allostery of `x`."""
function _testhelper_construction_bytes(x)
    cm, a = x.catalytic, x.allosteric
    _testhelper_construct() =
        ER.Mechanism(ER.Mechanisms.CatalyticScheme(cm.metabolites, cm.enzyme_forms,
                                                   cm.substrates, cm.products, cm.steps),
                     ER.Mechanisms.Allostery(a.allosteric_states, a.oligomeric_state,
                                             a.regulatory_sites))
    _testhelper_construct()
    @allocated _testhelper_construct()
end

@testset "constructor allocation budget" begin
    # Each construction allocates at most 1.3 times the bytes it allocated at bf49d6e1 on
    # Julia 1.13 (6,672, 22,264 and 25,912 bytes): a mechanism stores the facts it computes
    # once, and they may not make every construction in a search meaningfully dearer.
    budget(bytes) = ceil(Int, 1.3 * bytes)
    michaelis_menten = @enzyme_mechanism begin
        substrates: S
        products: P
        steps: begin
            E + S <--> ES
            ES <--> E + P
        end
    end
    @test _testhelper_construction_bytes(michaelis_menten) <= budget(6_672)
    # The largest seed of the bi-bi init set (A, B → P, Q): 9 forms, 13 steps.
    random_bi_bi = @enzyme_mechanism begin
        substrates: A, B
        products: P, Q
        steps: begin
            (E + A ⇌ EA, EB + A ⇌ EAB, EP + A ⇌ EAP)
            (E + B ⇌ EB, EA + B ⇌ EAB, EQ + B ⇌ EBQ)
            (E + P ⇌ EP, EA + P ⇌ EAP, EQ + P ⇌ EPQ)
            (E + Q ⇌ EQ, EB + Q ⇌ EBQ, EP + Q ⇌ EPQ)
            EAB <--> EPQ
        end
    end
    @test _testhelper_construction_bytes(random_bi_bi) <= budget(22_264)
    # Hexokinase: 14 steps, two :OnlyA groups and a regulatory site.
    hk = only(s for s in MECHANISM_TEST_SPECS if s.name == "HK").mechanism
    @test _testhelper_construction_bytes(hk) <= budget(25_912)
end
