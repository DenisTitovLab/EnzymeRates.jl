# ABOUTME: Tests of the reaction-graph mechanism types: content order and identity,
# ABOUTME: canonical form and checks, its sub-graph reads, the Sig and speed.

using Statistics: median

const ER = EnzymeRates

"""Uni-uni with rapid-equilibrium bindings, written field by field in canonical form:
nodes E, E(P), E(S); the bindings of P and S, then E(S) → E(P)."""
function _testhelper_graph_uniuni()
    M = ER.Mechanisms
    P, S, none = [M.Metabolite(:P, false)], [M.Metabolite(:S, false)], M.Metabolite[]
    _testhelper_plain_mechanism(
        [:S], [:P],
        [M.EnzymeForm(:E, none, Symbol[], Symbol[]),
         M.EnzymeForm(:E, P, Symbol[], Symbol[]),
         M.EnzymeForm(:E, S, Symbol[], Symbol[])],
        [M.Step(1, 2, P, none, true, 1), M.Step(1, 3, S, none, true, 2),
         M.Step(3, 2, none, none, false, 3)])
end

"""The metabolites `A`, `B`, `C`, `D` (`Metabolite`s, none a modifier), indexed 1 to 4."""
const _testhelper_METS = [ER.Mechanisms.Metabolite(x, false) for x in (:A, :B, :C, :D)]

"""The sorted metabolite list holding `_testhelper_METS[i]` once for each `i` of
`list`."""
_testhelper_mets(list) = sort!([_testhelper_METS[i] for i in list])

@testset "metabolite lists order lexicographically, a list before the lists it begins" begin
    # Every list over four metabolites with at most two copies of each: nodes that differ
    # in their bound lists alone order as Julia orders the lists of metabolite indices.
    lists = [reduce(vcat, [fill(i, n[i]) for i in 1:4]; init = Int[])
             for n in Iterators.product(0:2, 0:2, 0:2, 0:2)]
    node(list) = ER.Mechanisms.EnzymeForm(:E, _testhelper_mets(list), Symbol[], Symbol[])
    @test all(isless(node(a), node(b)) == isless(a, b) for a in lists, b in lists)
end

@testset "nodes order by label, bound, added, subtracted; steps by ends, lists, flag" begin
    M = ER.Mechanisms
    A, B, AB, none = _testhelper_mets([1]), _testhelper_mets([2]), _testhelper_mets([1, 2]),
                     M.Metabolite[]
    node(c, bound, added = Symbol[], subtracted = Symbol[]) =
        M.EnzymeForm(c, bound, added, subtracted)
    @test node(:E, B) < node(:F, none)
    @test node(:E, AB) < node(:E, B)
    @test node(:E, B) < node(:E, B, [:A]) < node(:E, B, [:A, :B])
    @test node(:E, B, [:A]) < node(:E, B, [:A], [:B])
    @test !(node(:E, B, [:A], [:B]) < node(:E, B, [:A], [:B]))
    @test M.Step(1, 3, none, none, false, 1) < M.Step(2, 1, none, none, false, 1)
    @test M.Step(1, 2, none, none, false, 1) < M.Step(1, 3, none, none, false, 1)
    @test M.Step(1, 2, AB, none, true, 1) < M.Step(1, 2, B, none, false, 1)
    @test M.Step(1, 2, A, AB, true, 1) < M.Step(1, 2, A, B, false, 1)
    @test M.Step(1, 2, A, none, false, 1) < M.Step(1, 2, A, none, true, 1)
    # The group takes no part in the order.
    @test !(M.Step(1, 2, A, none, false, 1) < M.Step(1, 2, A, none, false, 2)) &&
          !(M.Step(1, 2, A, none, false, 2) < M.Step(1, 2, A, none, false, 1))
    @test M._reversed(M.Step(1, 2, A, B, true, 3)) == M.Step(2, 1, B, A, true, 3)
end

"""Whether a value of type `T` hashes the same in every process: a symbol, an integer or
a flag, a vector of such values, or one of the types whose hash folds its name and its
fields (`_testhelper_hash_folds`), whose fields are such values in turn."""
function _testhelper_hash_is_content(T)
    T in (Symbol, Bool, Int) && return true
    T <: Vector && return _testhelper_hash_is_content(eltype(T))
    T in (EnzymeRates.Mechanisms.Metabolite, EnzymeRates.Mechanisms.RegulatorySite,
          ER.Mechanisms.EnzymeForm, ER.Mechanisms.Step, ER.Mechanisms.CatalyticScheme,
          ER.Mechanisms.Allostery, ER.Mechanism) || return false
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
    for x in (c.nodes[2], c.steps[3], c, m.allosteric, m, am.allosteric, am,
              c.nodes[2].bound[1], site)
        @test _testhelper_hash_folds(x, h)
    end
    @test all(
        _testhelper_hash_is_content,
        (ER.Mechanisms.EnzymeForm, ER.Mechanisms.Step, ER.Mechanisms.CatalyticScheme,
         ER.Mechanisms.Allostery, ER.Mechanism))
    # Each canonical field takes part in `==` and `hash`.
    @test m == _testhelper_graph_uniuni() && hash(m) == hash(_testhelper_graph_uniuni())
    canonical = _testhelper_identity_fields(ER.Mechanisms.CatalyticScheme)
    for f in canonical
        fields = Any[getfield(c, g) for g in canonical]
        i = findfirst(==(f), canonical)
        fields[i] = f === :substrates ? [:A] : f === :products ? [:Q] :
                    reverse(getfield(c, f))
        other = ER.Mechanisms._unchecked_scheme(fields...)
        @test other != c && hash(other) != hash(c)
    end
    # The steps' groups take part too.
    regrouped = ER.Mechanisms._unchecked_scheme(c.substrates, c.products, c.nodes,
        [ER.Mechanisms.Step(e.from, e.to, e.taken_up, e.given_off, e.rapid_equilibrium, g)
         for (e, g) in zip(c.steps, [1, 3, 2])])
    @test regrouped != c && hash(regrouped) != hash(c)
    @test am == ER.Mechanism(c, allostery([:EqualAI, :EqualAI, :OnlyA], 2, [site]))
    @test am != ER.Mechanism(c, allostery([:EqualAI, :OnlyA, :OnlyA], 2, [site]))
    @test am != ER.Mechanism(c, allostery([:EqualAI, :EqualAI, :OnlyA], 4, [site]))
    @test am != ER.Mechanism(c, allostery([:EqualAI, :EqualAI, :OnlyA], 2,
                                          EnzymeRates.Mechanisms.RegulatorySite[]))
end

@testset "a mechanism's fields are its graph and the facts it gives" begin
    # No reaction and no cache: a mechanism is its catalytic scheme and its allostery. The
    # scheme holds the names, the nodes and the steps with their groups, then the facts
    # computed from them once; the allostery the tags, multiplicity and sites.
    @test fieldnames(ER.Mechanism) == (:catalytic, :allosteric)
    @test fieldnames(ER.Mechanisms.CatalyticScheme) ==
          (:substrates, :products, :nodes, :steps, :enzyme_forms, :modifiers, :groups,
           :step_kinds, :toward_products, :step_sides)
    @test fieldnames(ER.Mechanisms.EnzymeForm) ==
          (:conformation, :bound, :added, :subtracted)
    @test fieldnames(ER.Mechanisms.Step) ==
          (:from, :to, :taken_up, :given_off, :rapid_equilibrium, :group)
    @test fieldnames(ER.Mechanisms.Allostery) ==
          (:allosteric_states, :oligomeric_state, :regulatory_sites)
end

@testset "construction is canonical" begin
    # Michaelis–Menten from nodes and steps: (a) in canonical form, (b) its nodes permuted
    # and one written twice, (c) every step reversed, (d) its group labels renumbered.
    # Every writing builds one mechanism, with one hash and one printed text.
    M = ER.Mechanisms
    S, P, none = [M.Metabolite(:S, false)], [M.Metabolite(:P, false)], M.Metabolite[]
    E, EP, ES = M.EnzymeForm(:E, none, Symbol[], Symbol[]),
                M.EnzymeForm(:E, P, Symbol[], Symbol[]),
                M.EnzymeForm(:E, S, Symbol[], Symbol[])
    steps = [M.Step(1, 2, P, none, true, 1), M.Step(1, 3, S, none, true, 2),
             M.Step(3, 2, none, none, false, 3)]
    canonical = _testhelper_plain_mechanism([:S], [:P], [E, EP, ES], steps)
    @test (canonical.catalytic.nodes, canonical.catalytic.steps) == ([E, EP, ES], steps)
    @test canonical == @enzyme_mechanism begin
        substrates: S
        products: P
        steps: begin
            E + S ⇌ E(S)
            E(S) <--> E(P)
            E(P) ⇌ E + P
        end
    end
    permuted = _testhelper_plain_mechanism([:S], [:P], [ES, E, EP, E],
                            [M.Step(4, 3, P, none, true, 1), M.Step(2, 1, S, none, true, 2),
                             M.Step(1, 3, none, none, false, 3)])
    reversed = _testhelper_plain_mechanism([:S], [:P], [E, EP, ES],
                            [M.Step(2, 1, none, P, true, 1), M.Step(3, 1, none, S, true, 2),
                             M.Step(2, 3, none, none, false, 3)])
    relabelled_steps = [M.Step(1, 2, P, none, true, 3), M.Step(1, 3, S, none, true, 1),
                        M.Step(3, 2, none, none, false, 2)]
    relabelled = _testhelper_plain_mechanism([:S], [:P], [E, EP, ES], relabelled_steps)
    for m in (permuted, reversed, relabelled)
        @test m == canonical && hash(m) == hash(canonical)
        @test sprint(show, m) == sprint(show, canonical)
    end
    # Each list of a node or a step in any order: ordered bi-uni with E(A, B) written as
    # E(B, A).
    A, B, Q = M.Metabolite(:A, false), M.Metabolite(:B, false), M.Metabolite(:Q, false)
    EA, EAB = M.EnzymeForm(:E, [A], Symbol[], Symbol[]),
              M.EnzymeForm(:E, [A, B], Symbol[], Symbol[])
    _testhelper_bi_uni(eab) = _testhelper_plain_mechanism(
        [:A, :B], [:Q], [M.EnzymeForm(:E, none, Symbol[], Symbol[]), EA, eab],
        [M.Step(1, 2, [A], none, true, 1), M.Step(2, 3, [B], none, true, 2),
         M.Step(3, 1, none, [Q], false, 3)])
    unsorted = _testhelper_bi_uni(M.EnzymeForm(:E, [B, A], Symbol[], Symbol[]))
    @test unsorted == _testhelper_bi_uni(EAB) && EAB in unsorted.catalytic.nodes
    # An allosteric mechanism from its written parts takes each tag by written group label:
    # label 1 (S's binding) is canonical group 2, label 2 (the chemistry) group 3 and
    # label 3 (P's binding) group 1.
    site = M.RegulatorySite([:R], 2, [:OnlyI])
    am = ER.Mechanism(
        ER.Mechanisms.CatalyticScheme([:S], [:P], [E, EP, ES], relabelled_steps),
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
end

"""`raw`, a spec's raw value (`raw_mechanism`), written another way, the same mechanism:
with `order`, its nodes and its steps in reverse order; with `direction`, each step from
its other end; with `labels`, its group labels numbered from the last. Its declared
names come in reverse order. Its tags follow their labels, and its sites and each site's
ligands come in reverse order."""
function _testhelper_rewrite(raw::NamedTuple; order = true, direction = true,
                             labels = true)
    M = EnzymeRates.Mechanisms
    nodes, steps, tags = raw.nodes, raw.steps, raw.allosteric_states
    G, n = maximum(e -> e.group, steps; init = 0), length(nodes)
    direction && (steps = [M.Step(e.to, e.from, e.given_off, e.taken_up,
                                  e.rapid_equilibrium, e.group) for e in steps])
    labels && ((steps, tags) = ([M.Step(e.from, e.to, e.taken_up, e.given_off,
                                        e.rapid_equilibrium, G + 1 - e.group)
                                 for e in steps], reverse(tags)))
    order && ((nodes, steps) = (reverse(nodes),
                                reverse([M.Step(n + 1 - e.from, n + 1 - e.to, e.taken_up,
                                                e.given_off, e.rapid_equilibrium, e.group)
                                         for e in steps])))
    (; substrates = reverse(raw.substrates), products = reverse(raw.products),
     modifiers = reverse(raw.modifiers), nodes, steps, allosteric_states = tags,
     oligomeric_state = raw.oligomeric_state,
     regulatory_sites = [(; ligands = reverse(s.ligands), s.multiplicity,
                          allosteric_states = reverse(s.allosteric_states))
                         for s in reverse(raw.regulatory_sites)])
end

@testset "every writing of a mechanism gives one canonical graph" begin
    @testset "$(spec.name)" for spec in MECHANISM_TEST_SPECS
        m = _testhelper_canonicalize(spec.raw).scheme
        for (order, direction, labels) in ((true, false, false), (false, true, false),
                                           (false, false, true), (true, true, true))
            other = _testhelper_canonicalize(
                _testhelper_rewrite(spec.raw; order, direction, labels)).scheme
            @test other == m && hash(other) == hash(m)
        end
    end
end

@testset "canonical form maps each written step to its edge, group and forms" begin
    c = _testhelper_canonicalize(EnzymeRates.DSL.raw_mechanism(:(begin
        substrates: S
        products: P
        steps: begin
            E + S ⇌ E(S)
            E(S) <--> E(P)
            E(P) ⇌ E + P
        end
    end)))
    # Forms E, E(P), E(S); the release of P is stored as the binding it reverses.
    @test c.scheme == _testhelper_graph_uniuni().catalytic
    @test c.edge_of == [2, 3, 1]
    @test c.flipped == [false, false, true]
    @test c.group_of == [2, 3, 1]
    @test c.node_of == [1, 3, 3, 2, 2, 1]
end

@testset "canonical form stores each step the way the graph gives it, a group as one" begin
    _testhelper_flipped(block) = _testhelper_canonicalize(
        EnzymeRates.DSL.raw_mechanism(block)).flipped
    # A step with one metabolite takes it up, so the fused release F(P) → F + P is stored
    # as F + P → F(P). The isomerizations lie on the one cycle that turns the reaction
    # over, E → E(S) → F(S) → F(P) → F → E, and run its way.
    @test _testhelper_flipped(:(begin
        substrates: S
        products: P
        steps: begin
            E + S <--> E(S)
            F(S) <--> E(S)
            F(S) <--> F(P)
            F(P) <--> F + P
            F <--> E
        end
    end)) == [false, true, false, true, false]
    @test _testhelper_flipped(:(begin
        substrates: S
        products: P
        steps: begin
            E(S) <--> E + S
            E(S) <--> F(S)
            F(P) <--> F(S)
            F + P <--> F(P)
            E <--> F
        end
    end)) == [true, false, true, false, true]
    # F → E runs the way the reaction turns over and leads its kinetic group, nearer free
    # E than its mirror at the I-bound forms, which keeps the direction it is written in
    # relative to F → E. Both written the other way round, the group is stored the same
    # way; written against each other, they write another mechanism. The release of P is
    # stored as the binding it reverses.
    _testhelper_mirrored(free, bound) = EnzymeRates.DSL.raw_mechanism(:(begin
        substrates: A
        products: P
        modifiers: I
        steps: begin
            E + A <--> E(A)
            E(A) <--> E(P)
            E(P) <--> F + P
            ($free, $bound)
            (E + I ⇌ E(I), F + I ⇌ F(I))
        end
    end))
    c = _testhelper_canonicalize(_testhelper_mirrored(:(F <--> E), :(F(I) <--> E(I))))
    d = _testhelper_canonicalize(_testhelper_mirrored(:(E <--> F), :(E(I) <--> F(I))))
    against = _testhelper_canonicalize(
        _testhelper_mirrored(:(F <--> E), :(E(I) <--> F(I))))
    @test c.flipped == [false, false, true, false, false, false, false]
    @test d.flipped == [false, false, true, true, true, false, false]
    @test c.scheme == d.scheme
    @test against.flipped == [false, false, true, false, false, false, false]
    @test against.scheme != c.scheme
end

"""Bindings of each metabolite of `mets` to free E, node 1, each into a node of its own and
in a group of its own, at rapid equilibrium: the nodes and the steps."""
function _testhelper_bindings(mets::Vector{ER.Mechanisms.Metabolite})
    M = ER.Mechanisms
    nodes = [M.EnzymeForm(:E, M.Metabolite[], Symbol[], Symbol[]);
             [M.EnzymeForm(:E, [x], Symbol[], Symbol[]) for x in mets]]
    nodes, [M.Step(1, i + 1, [x], M.Metabolite[], true, i) for (i, x) in enumerate(mets)]
end

@testset "a graph holds any number of metabolites and of copies of one" begin
    M = ER.Mechanisms
    names = [Symbol("S", i) for i in 1:32]
    plain(xs) = [M.Metabolite(x, false) for x in xs]
    _testhelper_entries(substrates, mets) = length(M._table(
        M._canonicalize(substrates, [:P], _testhelper_bindings(mets)...).scheme))
    @test _testhelper_entries(names[1:31], plain([names[1:31]; :P])) == 32
    @test _testhelper_entries(names, plain([names; :P])) == 33
    # A reactant and its modifier copy are two entries.
    @test _testhelper_entries(
        names[1:16],
        [plain([names[1:16]; :P]); [M.Metabolite(x, true) for x in names[1:16]]]) == 33
    _testhelper_holding(n) = M._canonicalize(
        [:A], [:P],
        [M.EnzymeForm(c, fill(M.Metabolite(:A, false), n), Symbol[], Symbol[])
         for c in (:E, :F)],
        [M.Step(1, 2, M.Metabolite[], M.Metabolite[], false, 1)])
    for n in (15, 16)
        m = _testhelper_holding(n).scheme
        @test count(x -> x.name === :A, m.nodes[1].bound) == n
    end
end

"""The raw value of `block`, a mechanism macro's block, read by the DSL walker."""
_testhelper_walk(block) = EnzymeRates.DSL.raw_mechanism(block)

@testset "the constructors build one value from every writing of a mechanism" begin
    @testset "$(spec.name)" for spec in MECHANISM_TEST_SPECS
        x = _testhelper_build(spec.raw)
        @test x.catalytic == _testhelper_canonicalize(spec.raw).scheme
        other = _testhelper_build(_testhelper_rewrite(spec.raw))
        @test other == x && hash(other) == hash(x)
    end
    # Each tag goes to the canonical group its label became.
    am = _testhelper_build(_testhelper_walk(:(begin
        substrates: S
        products: P
        oligomeric_state: 2
        steps: begin
            E + S ⇌ E(S) :: EqualAI
            E(S) <--> E(P) :: OnlyA
            E(P) ⇌ E + P :: NonequalAI
        end
    end)))
    @test am.allosteric.allosteric_states == [:NonequalAI, :EqualAI, :OnlyA]
end

"""Whether `v` takes the memory of its copy, which `copy` stores at its length."""
_testhelper_at_length(v) = Base.summarysize(v) == Base.summarysize(copy(v))

"""Whether `m`, a catalytic scheme or a mechanism, stores each of its vectors at its
length: `m` takes the memory of the value rebuilt from a copy of each vector, its nodes'
and steps' lists and its sites included, which `copy` stores at its length, and each fact
a scheme stores takes the memory of its copy."""
_testhelper_stored_at_length(m::ER.Mechanisms.CatalyticScheme) =
    Base.summarysize(m) == Base.summarysize(
        ER.Mechanisms._unchecked_scheme(copy(m.substrates), copy(m.products),
            [ER.Mechanisms.EnzymeForm(n.conformation, copy(n.bound), copy(n.added),
                                      copy(n.subtracted)) for n in m.nodes],
            [ER.Mechanisms.Step(e.from, e.to, copy(e.taken_up), copy(e.given_off),
                                e.rapid_equilibrium, e.group) for e in m.steps])) &&
    all(_testhelper_at_length, (m.modifiers, m.groups, m.step_kinds, m.toward_products,
                                m.step_sides)) && all(_testhelper_at_length, m.groups)
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
    # Uni-uni, E + S ⇌ E(S), E(S) <--> E(P), E(P) ⇌ E + P, nodes 1 to 3, with the nodes
    # `extra` (4, 5, …) and `steps` added; every name a node or a step binds as a modifier
    # is a modifier.
    M = ER.Mechanisms
    S, P, I = (M.Metabolite(x, false) for x in (:S, :P, :I))
    Imod, none = M.Metabolite(:I, true), M.Metabolite[]
    node(c, bound...) = M.EnzymeForm(c, collect(M.Metabolite, bound), Symbol[], Symbol[])
    nodes = [node(:E), node(:E, S), node(:E, P)]
    uniuni = [M.Step(1, 2, [S], none, true, 1), M.Step(2, 3, none, none, false, 2),
              M.Step(3, 1, none, [P], true, 3)]
    _testhelper_uniuni(extra = M.EnzymeForm[], steps...; substrates = [:S], products = [:P],
                       group = [1, 2, 3]) =
        _testhelper_plain_mechanism(substrates, products, [nodes; extra],
                     [[M.Step(e.from, e.to, e.taken_up, e.given_off, e.rapid_equilibrium, g)
                       for (e, g) in zip(uniuni, group)]; steps...])
    # Uni-uni with `x` bound to free E by a step of its own.
    _testhelper_binding(x; kw...) =
        _testhelper_uniuni([node(:E, x)], M.Step(1, 4, [x], none, true, 4); kw...)
    _testhelper_allosteric(tags, n, sites) =
        ER.Mechanism(ER.Mechanisms.CatalyticScheme([:S], [:P], nodes, uniuni),
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
    # A ligand takes no name a metabolite cannot take.
    @test_throws(ErrorException("Mechanism: \"1R\" cannot name a metabolite; a " *
                                "metabolite name is an identifier that is no keyword"),
                 _testhelper_allosteric(tags, 2, [site([Symbol("1R")], 2, [:OnlyA])]))
    @test_throws(ErrorException("Mechanism: E cannot name a metabolite; it " *
                                "labels the free enzyme form"),
                 _testhelper_allosteric(tags, 2, [site([:E], 2, [:OnlyA])]))
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
                     "Mechanism: without its :OnlyI steps E(S) <--> E(P) the " *
                     "active conformation turns no substrate into product")
    # Declarations and group labels.
    @test_throws(ErrorException("Mechanism: substrates must not be empty"),
                 _testhelper_uniuni(substrates = Symbol[]))
    @test_throws(ErrorException("Mechanism: products must not be empty"),
                 _testhelper_uniuni(products = Symbol[]))
    @test_throws(ErrorException("Mechanism: duplicate substrate names"),
                 _testhelper_uniuni(substrates = [:S, :S]))
    @test_throws(ErrorException("Mechanism: duplicate product names"),
                 _testhelper_uniuni(products = [:P, :P]))
    @test_throws(ErrorException("Mechanism: S named as both a substrate and a product; " *
                                "concentrations and constants are keyed by name"),
                 _testhelper_uniuni(products = [:P, :S]))
    @test_throws(ErrorException("Mechanism: each step takes one kinetic group label " *
                                "from 1 to the number of groups"),
                 _testhelper_uniuni(group = [0, 1, 2]))
    @test_throws(ErrorException("Mechanism: kinetic group 2 is empty; a kinetic group " *
                                "holds at least one step"),
                 _testhelper_uniuni(group = [1, 3, 3]))
    # A step's ends index the nodes.
    for (from, to) in ((2, 4), (0, 1))
        @test_throws(ErrorException("Step 4: its ends are nodes $from and $to, but there " *
                                    "are 3 nodes; a step's ends index the nodes"),
                     _testhelper_uniuni(M.EnzymeForm[],
                                        M.Step(from, to, none, none, false, 4)))
    end
    # Declared and bound names a metabolite cannot take: no identifier, a keyword, E, a
    # reserved name, a name holding `_`.
    _testhelper_unwritable(x) = ErrorException(
        "Mechanism: \"$x\" cannot name a metabolite; a metabolite name is an identifier " *
        "that is no keyword")
    @test_throws(_testhelper_unwritable("1I"),
                 _testhelper_binding(M.Metabolite(Symbol("1I"), true)))
    @test_throws(_testhelper_unwritable("end"), _testhelper_uniuni(products = [:P, :end]))
    @test_throws(ErrorException("Mechanism: E cannot name a metabolite; it labels the " *
                                "free enzyme form"),
                 _testhelper_uniuni(substrates = [:S, :E]))
    @test_throws(ErrorException("Mechanism: Keq is a reserved name; Keq, L, N, " *
                                "params, concs and v are kept for a rate equation's " *
                                "constants and local names"),
                 _testhelper_uniuni(substrates = [:S, :Keq]))
    @test_throws(ErrorException("Mechanism: the metabolite name S_1 holds `_`; constant " *
                                "names join metabolite names with `_`, so a metabolite " *
                                "name holds none"),
                 _testhelper_binding(M.Metabolite(:S_1, true)))
    # A name Julia reads as another (a micro sign as μ, a decomposed é as one character)
    # prints as text that builds another mechanism, so no DSL text can write it. The
    # strings are written with escapes: Julia would read the characters themselves as
    # the other names.
    _testhelper_read_as(who, x, read) = ErrorException(
        "$who: Julia reads the name \"$x\" as \"$read\", so no DSL text can write it; " *
        "use \"$read\"")
    @test_throws(_testhelper_read_as("Mechanism", "\u00b5M", "\u03bcM"),
                 _testhelper_binding(M.Metabolite(Symbol("\u00b5M"), true)))
    @test_throws(_testhelper_read_as("Mechanism", "e\u0301", "\u00e9"),
                 _testhelper_uniuni(substrates = [:S, Symbol("e\u0301")]))
    @test_throws(_testhelper_read_as("Mechanism", "\u025b", "\u03b5"),
                 _testhelper_allosteric(tags, 2, [site([Symbol("\u025b")], 2, [:OnlyA])]))
    # Names: undeclared; a reactant's modifier copy; a plain entry of a name bound as a
    # modifier; residuals by role.
    @test_throws(ErrorException("Mechanism: P is not declared; declare it under " *
                                "substrates:, products: or modifiers:"),
                 _testhelper_uniuni(products = [:Q]))
    @test !_testhelper_is_allosteric(_testhelper_binding(M.Metabolite(:P, true)))
    @test !_testhelper_is_allosteric(_testhelper_binding(Imod))
    @test_throws(ErrorException("Mechanism: I is a modifier, so it binds only as a " *
                                "modifier"),
                 _testhelper_uniuni([node(:E, Imod), node(:E, I)],
                                    M.Step(1, 4, [Imod], none, true, 4),
                                    M.Step(1, 5, [I], none, true, 5)))
    for (added, subtracted, verb, x) in (([:P], Symbol[], "adds", :P),
                                         (Symbol[], [:S], "subtracts", :S))
        @test_throws(
            ErrorException("Mechanism: a residual adds substrates and subtracts " *
                           "products; one $verb $x"),
            _testhelper_uniuni([M.EnzymeForm(:E, none, added, subtracted)],
                               M.Step(2, 4, none, none, false, 4)))
    end
    # Steps: a self-loop, a metabolite taken up and given off, a role change.
    @test_throws(ErrorException("Step: both ends are the form EA; a step joins two forms"),
                 _testhelper_plain_mechanism(
                     [:A], [:P], [node(:E, M.Metabolite(:A, false))],
                     [M.Step(1, 1, none, none, false, 1)]))
    @test_throws(ErrorException("Step: both ends are the form ES; a step joins two forms"),
                 _testhelper_uniuni(M.EnzymeForm[], M.Step(2, 2, none, none, false, 4)))
    @test_throws(ErrorException("Step EA → EQ: B both consumed and released"),
                 _testhelper_build(_testhelper_walk(:(begin
                     substrates: A, B
                     products: Q
                     steps: begin
                         E + A ⇌ E(A)
                         E(A) + B <--> E(Q) + B
                         E(Q) ⇌ E + Q
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
                         E + A ⇌ E(A)
                         E(A) + B <--> E(Q) + B + B
                         E(Q) ⇌ E + Q
                     end
                 end))))
    _testhelper_role_change(from, to) = ErrorException(
        "Step $from → $to changes a metabolite's role: $to holds $from's metabolites " *
        "plus P by name, not by role; write a modifier copy as X::Mod on both sides")
    for (line, from, to) in ((:(E + P::Mod ⇌ E(P)), :E, :EP),
                             (:(E(P) ⇌ E + P::Mod), :E, :EP),
                             (:(E + P ⇌ E(P::Mod)), :E, :EPmod),
                             (:(E(P) + P::Mod ⇌ E(P, P)), :EP, :EPP))
        block = :(begin
            substrates: S
            products: P
            modifiers: P
            steps: begin
                E + S ⇌ E(S)
                E(S) <--> E(P)
                E(P) ⇌ E + P
                $line
            end
        end)
        @test_throws(_testhelper_role_change(from, to),
                     _testhelper_build(_testhelper_walk(block)))
    end
    # Two forms that render one name: bound names that concatenate, a label that holds a
    # bound name, a reactant's modifier copy beside a modifier named Xmod.
    _testhelper_one_name(a, b, name) = ErrorException(
        "Mechanism: the enzyme forms $a and $b both render the name $name; every " *
        "constant is named after the forms its step joins, so rename a metabolite to " *
        "tell the two forms apart")
    @test_throws _testhelper_one_name("E(NAD, P)", "E(NADP)", "ENADP") _testhelper_build(
        _testhelper_walk(:(begin
            substrates: NAD, P
            products: NADP
            steps: begin
                E + NAD ⇌ E(NAD)
                E(NAD) + P ⇌ E(NAD, P)
                E(NAD, P) <--> E(NADP)
                E(NADP) ⇌ E + NADP
            end
        end)))
    @test_throws _testhelper_one_name("E(star)", "Estar", "Estar") _testhelper_build(
        _testhelper_walk(:(begin
            substrates: star
            products: P
            steps: begin
                E + star ⇌ E(star)
                E(star) <--> E(P)
                E(P) ⇌ E + P
                E <--> Estar
            end
        end)))
    @test_throws _testhelper_one_name("E(P::Mod)", "E(Pmod)", "EPmod") _testhelper_build(
        _testhelper_walk(:(begin
            substrates: S
            products: P
            modifiers: P, Pmod
            steps: begin
                E + S ⇌ E(S)
                E(S) <--> E(P)
                E(P) ⇌ E + P
                E + P::Mod ⇌ E(P::Mod)
                E + Pmod ⇌ E(Pmod)
            end
        end)))
    # A kinetic group of steps of different kinds or flags.
    for (steps, text) in (
            (:(begin
                 E + A <--> E(A)
                 (E(A) + B <--> E(A, B), E(A) + B <--> E(Q) + P)
                 E(Q) <--> E + Q
             end), "EA_B → EAB (SS) and EA_B → EQ_P (SS)"),
            (:(begin
                 (E + A ⇌ E(A), E(Q) + A <--> E(A, Q))
                 E(A) <--> E(Q) + P
                 E(Q) <--> E + Q
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
            E + S <--> E(S)
            E + S ⇌ E(S)
            E(S) ⇌ E(P)
            E(P) ⇌ E + P
        end)))
    @test_throws(ErrorException(
        "Mechanism: kinetic groups 3 and 4 both hold the reaction ES ⇌ EP; a reaction " *
        "belongs to one kinetic group, whose constants are named after it"),
        _testhelper_twice(:(begin
            E + S <--> E(S)
            E(S) <--> E(P)
            E(P) <--> E(S)
            E(P) <--> E + P
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
                E + A <--> E(A)
                E(A) <--> E(P)
                E(P) <--> F + P
                F <--> E
                (E <--> F, E(I) <--> F(I))
                E + I ⇌ E(I)
                F + I ⇌ F(I)
                E(I) + P <--> E(I, P)
                F(I) + A <--> F(A, I)
            end
        end))))
    @test_throws(ErrorException(
        "Mechanism: kinetic group 2 holds the reaction E_S ⇌ ES twice; a mechanism holds " *
        "each reaction once"),
        _testhelper_twice(:(begin
            (E + S <--> E(S), E(S) <--> E + S)
            E(S) <--> E(P)
            E(P) <--> E + P
        end)))
    # No free E: no step at all, a free form labelled F, E only with S or P bound, or E
    # with nothing bound only where it holds a residual.
    no_free_E = ErrorException("Mechanism: no free enzyme form E; every mechanism holds " *
                               "the free form E, with nothing bound and no residual")
    @test_throws no_free_E _testhelper_plain_mechanism([:S], [:P], M.EnzymeForm[], M.Step[])
    @test_throws no_free_E _testhelper_build(_testhelper_walk(:(begin
        substrates: S
        products: P
        steps: begin
            F + S ⇌ F(S)
            F(S) <--> F(P)
            F(P) ⇌ F + P
        end
    end)))
    @test_throws no_free_E _testhelper_build(_testhelper_walk(:(begin
        substrates: S
        products: P
        steps: begin
            F + S ⇌ E(S)
            E(S) <--> E(P)
            E(P) ⇌ F + P
        end
    end)))
    # A ping-pong enzyme whose free form is F: every declared name takes part, and E
    # with nothing bound appears only after A's chemistry has left its residual.
    @test_throws no_free_E _testhelper_build(_testhelper_walk(:(begin
        substrates: A, B
        products: P, Q
        steps: begin
            F + A ⇌ F(A)
            F(A) <--> E(; residual = A - P) + P
            E(; residual = A - P) + B ⇌ E(B; residual = A - P)
            E(B; residual = A - P) <--> F(Q)
            F(Q) ⇌ F + Q
        end
    end)))
    # Forms that no steps join to free E: an isomerization between unreached forms.
    @test_throws(ErrorException("Mechanism: the steps do not form one connected graph; " *
                                "no steps lead from free E to EstarP"),
                 _testhelper_build(_testhelper_walk(:(begin
        substrates: S
        products: P
        steps: begin
            E + S ⇌ E(S)
            E(S) <--> E(P)
            E(P) ⇌ E + P
            Estar(S) <--> Estar(P)
        end
    end))))
    # Declared names no step exchanges: a substrate, a product; a substrate exchanged
    # only as its modifier copy; a modifier, and a reactant's modifier copy, bound only
    # in a node.
    _testhelper_unused(kind, x) = ErrorException(kind == "modifier" ?
        "Mechanism: modifier $x is bound at a node, but no step binds it" :
        "Mechanism: $kind $x is declared, but no step uses it")
    @test_throws(_testhelper_unused("substrate", :T),
                 _testhelper_uniuni(substrates = [:S, :T]))
    @test_throws(_testhelper_unused("product", :Q), _testhelper_uniuni(products = [:P, :Q]))
    @test_throws(_testhelper_unused("substrate", :T),
                 _testhelper_binding(M.Metabolite(:T, true); substrates = [:S, :T]))
    for x in (Imod, M.Metabolite(:P, true))
        @test_throws(_testhelper_unused("modifier", x.name),
                     _testhelper_uniuni([node(:E, S, x)],
                                        M.Step(2, 4, none, none, false, 4)))
    end
    # Conformation labels: a label shaped otherwise than E, Estar or F2; an opaque
    # bound-form name; A or I, which constant names use as conformation tags; a
    # metabolite's name, of the catalytic graph or of a ligand.
    _testhelper_turning_to(label) =
        _testhelper_uniuni([node(label)], M.Step(1, 4, none, none, false, 4))
    @test_throws(ErrorException("Mechanism: `Ex_1` is not a conformation label: a " *
                                "conformation label is a capital letter followed by " *
                                "lowercase letters and digits, such as `E`, `Estar` or " *
                                "`F2`."),
                 _testhelper_turning_to(:Ex_1))
    @test_throws(ErrorException("Mechanism: `ES` looks like an opaque bound-form name; " *
                                "write it as decomposed call notation, e.g. `E(S)` or " *
                                "`E(A, B)`."),
                 _testhelper_turning_to(:ES))
    for label in (:A, :I)
        @test_throws(ErrorException("Mechanism: `$label` cannot be a conformation label: " *
                                    "constant names tag a conformation with A_ or I_ " *
                                    "(`k_A_ES_to_EP`), so the labels A and I could give " *
                                    "two constants one name; choose another label."),
                     _testhelper_turning_to(label))
    end
    @test_throws(ErrorException("Mechanism: conformation label `S` collides with " *
                                "declared metabolite `S`; choose a different " *
                                "conformation label."),
                 _testhelper_turning_to(:S))
    @test_throws(ErrorException("Mechanism: conformation label `F` collides " *
                                "with ligand `F`; choose a different conformation " *
                                "label."),
                 ER.Mechanism(
                     ER.Mechanisms.CatalyticScheme(
                         [:S], [:P], [nodes; node(:F)],
                         [uniuni; M.Step(1, 4, none, none, false, 4)]),
                     ER.Mechanisms.Allostery([tags; :EqualAI], 2,
                                             [site([:F], 2, [:OnlyA])])))
end

@testset "a form's name and text read as constants and the macros write them" begin
    m = _testhelper_build(_testhelper_walk(:(begin
        substrates: A, B
        products: P, Q
        modifiers: I, Q
        steps: begin
            E + A <--> E(A)
            E(A) <--> E(; residual = A - P) + P
            E(; residual = A - P) + B <--> E(B; residual = A - P)
            E(B; residual = A - P) <--> E(Q)
            E(Q) <--> E + Q
            E(Q) + Q::Mod ⇌ E(Q, Q::Mod)
            E + I ⇌ E(I)
        end
    end)))
    # Forms in content order: a residual form binds nothing, so it follows free E.
    @test [ER.Mechanisms._form_name(m.catalytic, f)
           for f in eachindex(m.catalytic.nodes)] ==
          ["E", "E_res_plusA_minusP", "EA", "EB_res_plusA_minusP", "EI", "EQ", "EQQmod"]
    @test [ER.Mechanisms.form_text(m.catalytic, f) for f in eachindex(m.catalytic.nodes)] ==
          ["E", "E(; residual = A - P)", "E(A)", "E(B; residual = A - P)", "E(I)",
           "E(Q)", "E(Q, Q::Mod)"]
    # Each edge's sides in its stored direction: releases stored as bindings.
    @test [join(m.catalytic.step_sides[k], " → ") for k in eachindex(m.catalytic.steps)] ==
          ["E_A → EA", "E_I → EI", "E_Q → EQ", "E_res_plusA_minusP_P → EA",
           "E_res_plusA_minusP_B → EB_res_plusA_minusP", "EB_res_plusA_minusP → EQ",
           "EQ_Qmod → EQQmod"]
end

"""The graph whose Sig (`Derivation._sig_of`) is `sig`. `sig` is splatted with its length
left open, as a compiled law's generated methods splat it: inference does not specialize
on it, so the decoder compiles once, whatever the Sig length."""
Base.@nospecializeinfer _testhelper_decode(@nospecialize(sig)) =
    ER.Derivation._mechanism_from_sig(sig...)

@testset "a graph's Sig decodes to the graph, and equal graphs share one law" begin
    @testset "$(spec.name)" for spec in MECHANISM_TEST_SPECS
        x = _testhelper_build(spec.raw)
        sig = ER.Derivation._sig_of(x)
        @test sig isa Tuple{Vararg{Int}}
        @test _testhelper_decode(sig) == x
        rate = compile_rate_equation(x)
        @test typeof(rate).parameters[1] === sig
        @test _testhelper_decode(typeof(rate).parameters[1]) == x
        @test compile_rate_equation(_testhelper_build(_testhelper_rewrite(spec.raw))) ===
              rate
    end
    m = _testhelper_graph_uniuni()
    @test _testhelper_decode(ER.Derivation._sig_of(m)) == m
    am = ER.Mechanism(
        m.catalytic, ER.Mechanisms.Allostery([:EqualAI, :NonequalAI, :OnlyA], 4,
        [EnzymeRates.Mechanisms.RegulatorySite([:X, :R], 2, [:OnlyA, :OnlyI])]))
    @test _testhelper_decode(ER.Derivation._sig_of(am)) == am
    # A name outside ASCII, sorting last among 32 metabolite entries.
    names = [[Symbol("S", i) for i in 1:30]; :αKG]
    wide = _testhelper_plain_mechanism(names, [:P], _testhelper_bindings(
        [EnzymeRates.Mechanisms.Metabolite(x, false) for x in [names; :P]])...)
    @test ER.Mechanisms._table(wide.catalytic)[32] ==
          EnzymeRates.Mechanisms.Metabolite(:αKG, false)
    @test _testhelper_decode(ER.Derivation._sig_of(wide)) == wide
end

"""Uni-uni whose substrate's name holds `n` characters, so its Sig is longer than any
other test's: a Sig length no other test decodes."""
_testhelper_long_named(n) = let S = Symbol("S", "x"^(n - 1))
    _testhelper_build(_testhelper_walk(:(begin
        substrates: $S
        products: P
        steps: begin
            E + $S ⇌ E($S)
            E($S) <--> E(P)
            E(P) ⇌ E + P
        end
    end)))
end

@testset "the Sig decoder compiles once, whatever a Sig's length" begin
    # A compiled law's generated methods splat its Sig with its length left open, so the
    # varargs decoder compiles once, whatever the Sig length: after one decode, decoding
    # the Sigs of three lengths no other test decodes adds no method instance. A decoder
    # that took its Sig as an `NTuple{N, Int}` would compile once per length.
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

@testset "a graph written back as its nodes and steps builds the same graph" begin
    @testset "$(spec.name)" for spec in MECHANISM_TEST_SPECS
        x = _testhelper_build(spec.raw)
        raw = _testhelper_written(x)
        @test _testhelper_build(raw) == x
        # Its steps are the edges in stored order and direction, labelled by group.
        c = _testhelper_canonicalize(raw)
        @test c.edge_of == eachindex(c.edge_of) && !any(c.flipped)
        @test c.group_of == 1:length(x.catalytic.groups)
        @test c.node_of == eachindex(c.node_of)
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
    mods = m.modifiers
    lines = ["substrates: $(join(m.substrates, ", ")) | products: " *
             "$(join(m.products, ", "))" *
             (isempty(mods) ? "" : " | modifiers: $(join(mods, ", "))"),
             "forms: " *
             join([ER.Mechanisms._form_name(m, f) for f in eachindex(m.nodes)], ", ")]
    for k in eachindex(m.steps)
        g = m.steps[k].group
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
        names = ER.Mechanisms._distinct_names(m)
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
            E + A <--> E(A)
            E(A) <--> E(P)
            E(P) <--> D + P
            D <--> E
        end
    end))
    x = _testhelper_build(raw)
    none = EnzymeRates.Mechanisms.Metabolite[]
    @test x.catalytic.nodes[1:2] == [ER.Mechanisms.EnzymeForm(:E, none, Symbol[], Symbol[]),
                           ER.Mechanisms.EnzymeForm(:D, none, Symbol[], Symbol[])]
    @test _testhelper_graph(x).forms[1:2] == ["E", "D"]
    # The census's first populations: each graph decodes from its Sig and builds again
    # from the nodes and steps it writes back. The count: the six reactions' init
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
    # exchange lies on a cycle: each runs from its form nearer free E, E(I) and F(I),
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
        E + A <--> E(A)
        E(A) <--> E(P)
        E(P) <--> F + P
        ($free, $bound)
        (E + I ⇌ E(I), F + I ⇌ F(I))
        E(I) + A ⇌ E(A, I)
        E(J) + P ⇌ E(J, P)
        $exchange
    end)
    blocks = [
        :(begin
            substrates: S
            products: P
            steps: begin
                E + S <--> E(S)
                F(S) <--> E(S)
                F(S) <--> F(P)
                F(P) <--> F + P
                F <--> E
            end
        end),
        :(begin
            substrates: S
            products: P
            steps: begin
                E(S) <--> F(S)
                E + S <--> E(S)
                F + P <--> F(P)
                F(P) <--> F(S)
                E <--> F
            end
        end)]
    for (free, bound, exchange) in (
            (:(F <--> E), :(F(I) <--> E(I)),
             :((E(I) + J <--> E(J) + I, F(I) + J <--> F(J) + I))),
            (:(E <--> F), :(E(I) <--> F(I)),
             :((F(J) + I <--> F(I) + J, E(J) + I <--> E(I) + J))))
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
    for exchange in (:((E(I) + J <--> E(J) + I, F(I) + J <--> F(J) + I)),
                     :((F(J) + I <--> F(I) + J, E(J) + I <--> E(I) + J)))
        raw = _testhelper_walk(:(begin
            substrates: A
            products: P
            modifiers: I, J
            steps: begin
                E + A <--> E(A)                  :: EqualAI
                E(A) <--> E(P)                   :: EqualAI
                E(P) <--> F + P                  :: EqualAI
                F ⇌ E                            :: EqualAI
                (E + I ⇌ E(I), F + I ⇌ F(I))     :: EqualAI
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

"""Whether the facts the catalytic scheme `m` stores are those its canonical fields give,
and the stored fields hold them: its form names, modifiers and groups, read off its nodes
and steps here, and each step's kind, direction toward products and side names, as the
graph reads compute them from the canonical fields."""
function _testhelper_facts_hold(m::ER.Mechanisms.CatalyticScheme)
    M = ER.Mechanisms
    lists = [[n.bound for n in m.nodes]; [e.taken_up for e in m.steps];
             [e.given_off for e in m.steps]]
    modifiers = sort!(unique!(Symbol[x.name for list in lists for x in list if x.modifier]))
    groups = [findall(e -> e.group == g, m.steps)
              for g in 1:maximum(e -> e.group, m.steps)]
    m.enzyme_forms == [Symbol(M._form_name(m, f)) for f in eachindex(m.nodes)] &&
        m.modifiers == modifiers && m.groups == groups &&
        m.step_kinds == [M._step_kind(m, e) for e in m.steps] &&
        m.toward_products == [M._runs_toward_products(m, e) for e in m.steps] &&
        m.step_sides == [M._sides(m, e) for e in m.steps]
end

"""Whether the facts of the mechanism `am` are those its fields give: its catalytic
graph's (`_testhelper_facts_hold`), and for an allosteric mechanism per conformation the
steps it runs (`steps_in_active_conformation`, `steps_in_inactive_conformation`): those
its tags do not exclude (`_testhelper_excludes`), on the forms they join to free E
(`_testhelper_exclusion`)."""
function _testhelper_facts_hold(am::ER.Mechanism)
    cm = am.catalytic
    _testhelper_facts_hold(cm) || return false
    _testhelper_is_allosteric(am) || return true
    M = ER.Mechanisms
    all(((1, M.steps_in_active_conformation(am)),
         (2, M.steps_in_inactive_conformation(am)))) do (c, runs)
        excluded = BitVector([_testhelper_excludes(am, k, c) for k in eachindex(cm.steps)])
        component = _testhelper_exclusion(am, c).component
        runs == BitVector([!excluded[k] && component[e.from] == 1
                           for (k, e) in enumerate(cm.steps)])
    end
end

"""Whether `a` and `b` hold equal values in every field, a catalytic scheme's and an
allostery's included."""
_testhelper_all_fields_equal(a, b) = all(fieldnames(typeof(a))) do f
    x, y = getfield(a, f), getfield(b, f)
    x isa Union{ER.Mechanisms.CatalyticScheme, ER.Mechanisms.Allostery} ?
        _testhelper_all_fields_equal(x, y) : x == y
end

@testset "stored facts equal what the canonical fields give" begin
    # Every construction path: the constructors (the specs' writings, the bi-bi seeds)
    # and the Sig decoder (the specs' compiled laws).
    @testset "$(spec.name)" for spec in MECHANISM_TEST_SPECS
        @test _testhelper_facts_hold(_testhelper_build(spec.raw))
        @test _testhelper_facts_hold(
            _testhelper_decode(ER.Derivation._sig_of(spec.mechanism)))
    end
    rxn = @enzyme_reaction begin
        substrates: A[C], B[N]
        products: P[C], Q[N]
    end
    @test all(_testhelper_facts_hold,
              _testhelper_mechanisms(EnzymeRates.Enumeration.init_mechanisms(rxn)))
end

@testset "stored facts survive the Sig round trip" begin
    @testset "$(spec.name)" for spec in MECHANISM_TEST_SPECS
        x = _testhelper_build(spec.raw)
        @test _testhelper_all_fields_equal(_testhelper_decode(ER.Derivation._sig_of(x)), x)
    end
    rxn = @enzyme_reaction begin
        substrates: A[C], B[N]
        products: P[C], Q[N]
    end
    @test all(x -> _testhelper_all_fields_equal(
                  _testhelper_decode(ER.Derivation._sig_of(x)), x),
              _testhelper_mechanisms(EnzymeRates.Enumeration.init_mechanisms(rxn)))
end

@testset "derived fields stay out of == and hash" begin
    # A mechanism's facts are functions of its canonical fields, so two mechanisms with
    # equal canonical fields are equal and hash alike whatever their facts hold: here a
    # rebuilt copy whose every fact is overwritten.
    m = @enzyme_mechanism begin
        substrates: S
        products: P
        steps: begin
            E + S ⇌ E(S)
            E(S) <--> E(P)
            E(P) ⇌ E + P
        end
    end
    c = m.catalytic
    x = ER.Mechanisms._unchecked_scheme(c.substrates, c.products, c.nodes, c.steps)
    x.enzyme_forms[2] = :F
    push!(x.modifiers, :S)
    x.groups[1] = [2, 3]
    x.step_kinds[1] = :chemistry
    x.toward_products .= .!x.toward_products
    x.step_sides[1] = ("E", "E")
    @test x == c && hash(x) == hash(c)
    @test length(unique!([c, x])) == 1
    # A mechanism stores no fact of its own: one around the overwritten copy equals one
    # around the original.
    am = ER.Mechanism(c, ER.Mechanisms.Allostery([:EqualAI, :OnlyA, :NonequalAI], 2,
                                                 ER.Mechanisms.RegulatorySite[]))
    y = ER.Mechanism(x, ER.Mechanisms.Allostery([:EqualAI, :OnlyA, :NonequalAI], 2,
                                                ER.Mechanisms.RegulatorySite[]))
    @test y == am && hash(y) == hash(am)
    @test length(unique!([am, y])) == 1
end

@testset "nodes and steps are the storage" begin
    m = _testhelper_build(_testhelper_walk(:(begin
        substrates: S
        products: P
        steps: begin
            E + S ⇌ E(S)
            E(S) <--> E(P)
            E(P) ⇌ E + P
        end
    end)))
    @test fieldnames(ER.Mechanisms.CatalyticScheme)[1:4] ==
          (:substrates, :products, :nodes, :steps)
    @test length(m.catalytic.nodes) == 3
    @test m.catalytic.nodes[1] == ER.Mechanisms.EnzymeForm(
        :E, EnzymeRates.Mechanisms.Metabolite[], Symbol[], Symbol[])
    # Nodes E, E(P), E(S); the release of P is stored as the binding it reverses, and
    # sorts first.
    @test [n.bound for n in m.catalytic.nodes[2:3]] ==
          [[EnzymeRates.Mechanisms.Metabolite(:P, false)],
           [EnzymeRates.Mechanisms.Metabolite(:S, false)]]
    @test m.catalytic.steps[1].taken_up == [EnzymeRates.Mechanisms.Metabolite(:P, false)]
    @test m.catalytic.steps[2].taken_up == [EnzymeRates.Mechanisms.Metabolite(:S, false)]
    @test (m.catalytic.steps[2].from, m.catalytic.steps[2].to,
           m.catalytic.steps[2].given_off) == (1, 3, EnzymeRates.Mechanisms.Metabolite[])
    @test m.catalytic.steps[1].rapid_equilibrium &&
          m.catalytic.steps[2].rapid_equilibrium && !m.catalytic.steps[3].rapid_equilibrium
    @test [s.group for s in m.catalytic.steps] == [1, 2, 3]
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
warm-up call: `Mechanism(CatalyticScheme(substrates, products, nodes, steps), allostery)`
with the allostery of `x`."""
function _testhelper_construction_bytes(x)
    cm, a = x.catalytic, x.allosteric
    _testhelper_construct() =
        ER.Mechanism(ER.Mechanisms.CatalyticScheme(cm.substrates, cm.products, cm.nodes,
                                                   cm.steps),
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
            E + S <--> E(S)
            E(S) <--> E + P
        end
    end
    @test _testhelper_construction_bytes(michaelis_menten) <= budget(6_672)
    # The largest seed of the bi-bi init set (A, B → P, Q): 9 forms, 13 steps.
    random_bi_bi = @enzyme_mechanism begin
        substrates: A, B
        products: P, Q
        steps: begin
            (E + A ⇌ E(A), E(B) + A ⇌ E(A, B), E(P) + A ⇌ E(A, P))
            (E + B ⇌ E(B), E(A) + B ⇌ E(A, B), E(Q) + B ⇌ E(B, Q))
            (E + P ⇌ E(P), E(A) + P ⇌ E(A, P), E(Q) + P ⇌ E(P, Q))
            (E + Q ⇌ E(Q), E(B) + Q ⇌ E(B, Q), E(P) + Q ⇌ E(P, Q))
            E(A, B) <--> E(P, Q)
        end
    end
    @test _testhelper_construction_bytes(random_bi_bi) <= budget(22_264)
    # Hexokinase: 14 steps, two :OnlyA groups and a regulatory site.
    hk = only(s for s in MECHANISM_TEST_SPECS if s.name == "HK").mechanism
    @test _testhelper_construction_bytes(hk) <= budget(25_912)
end
