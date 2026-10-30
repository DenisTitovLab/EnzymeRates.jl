# ABOUTME: Tests of the reaction-graph mechanism types: content order and identity,
# ABOUTME: canonical form and checks, its sub-graph reads, the Sig and speed.

using Statistics: median

const ER = EnzymeRates

"""Uni-uni with rapid-equilibrium bindings, written field by field in canonical form:
nodes E, E(P), E(S); the bindings of P and S, then E(S) → E(P)."""
function _testhelper_graph_uniuni()
    M = ER.Mechanisms
    P, S, none = [M.Met(:P, false)], [M.Met(:S, false)], M.Met[]
    ER.Mechanism(
        [:S], [:P],
        [M.Node(:E, none, Symbol[], Symbol[]), M.Node(:E, P, Symbol[], Symbol[]),
         M.Node(:E, S, Symbol[], Symbol[])],
        [M.Step(1, 2, P, none, true, 1), M.Step(1, 3, S, none, true, 2),
         M.Step(3, 2, none, none, false, 3)])
end

"""The metabolites `A`, `B`, `C`, `D` (`Met`s, none a modifier), indexed 1 to 4."""
const _testhelper_METS = [ER.Mechanisms.Met(x, false) for x in (:A, :B, :C, :D)]

"""The sorted metabolite list holding `_testhelper_METS[i]` once for each `i` of
`list`."""
_testhelper_mets(list) = sort!([_testhelper_METS[i] for i in list])

@testset "metabolite lists order lexicographically, a list before the lists it begins" begin
    # Every list over four metabolites with at most two copies of each: nodes that differ
    # in their bound lists alone order as Julia orders the lists of metabolite indices.
    lists = [reduce(vcat, [fill(i, n[i]) for i in 1:4]; init = Int[])
             for n in Iterators.product(0:2, 0:2, 0:2, 0:2)]
    node(list) = ER.Mechanisms.Node(:E, _testhelper_mets(list), Symbol[], Symbol[])
    @test all(isless(node(a), node(b)) == isless(a, b) for a in lists, b in lists)
end

@testset "nodes order by label, bound, added, subtracted; steps by ends, lists, flag" begin
    M = ER.Mechanisms
    A, B, AB, none = _testhelper_mets([1]), _testhelper_mets([2]), _testhelper_mets([1, 2]),
                     M.Met[]
    node(c, bound, added = Symbol[], subtracted = Symbol[]) =
        M.Node(c, bound, added, subtracted)
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
    T in (EnzymeRates.Mechanisms.Met, EnzymeRates.Mechanisms.RegulatorySite,
          ER.Mechanisms.Node, ER.Mechanisms.Step, ER.Mechanism,
          ER.AllostericMechanism) || return false
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
    site = EnzymeRates.Mechanisms.RegulatorySite([:R], 2, [:OnlyI])
    am = ER.AllostericMechanism(m, [:EqualAI, :EqualAI, :OnlyA], 2, [site])
    h = UInt(17)
    for x in (m.nodes[2], m.steps[3], m, am, m.nodes[2].bound[1], site)
        @test _testhelper_hash_folds(x, h)
    end
    @test all(
        _testhelper_hash_is_content,
        (ER.Mechanisms.Node, ER.Mechanisms.Step, ER.Mechanism, ER.AllostericMechanism))
    # Each canonical field takes part in `==` and `hash`.
    @test m == _testhelper_graph_uniuni() && hash(m) == hash(_testhelper_graph_uniuni())
    canonical = _testhelper_identity_fields(ER.Mechanism)
    for f in canonical
        fields = Any[getfield(m, g) for g in canonical]
        i = findfirst(==(f), canonical)
        fields[i] = f === :substrates ? [:A] : f === :products ? [:Q] :
                    reverse(getfield(m, f))
        other = ER.Mechanisms._unchecked_mechanism(fields...)
        @test other != m && hash(other) != hash(m)
    end
    # The steps' groups take part too.
    regrouped = ER.Mechanisms._unchecked_mechanism(m.substrates, m.products, m.nodes,
        [ER.Mechanisms.Step(e.from, e.to, e.taken_up, e.given_off, e.rapid_equilibrium, g)
         for (e, g) in zip(m.steps, [1, 3, 2])])
    @test regrouped != m && hash(regrouped) != hash(m)
    @test am == ER.AllostericMechanism(m, [:EqualAI, :EqualAI, :OnlyA], 2, [site])
    @test am != ER.AllostericMechanism(m, [:EqualAI, :OnlyA, :OnlyA], 2, [site])
    @test am != ER.AllostericMechanism(m, [:EqualAI, :EqualAI, :OnlyA], 4, [site])
    @test am != ER.AllostericMechanism(m, [:EqualAI, :EqualAI, :OnlyA], 2,
                                       EnzymeRates.Mechanisms.RegulatorySite[])
end

@testset "a mechanism's fields are its graph and the facts it gives" begin
    # No reaction and no cache: the names, the nodes and the steps with their groups, then
    # the facts computed from them once; for an allosteric mechanism its tags,
    # multiplicity and sites, then each conformation's steps and the steps it excludes.
    @test fieldnames(ER.Mechanism) == (:substrates, :products, :nodes, :steps, :modifiers,
                                       :groups, :step_kinds, :toward_products, :step_sides)
    @test fieldnames(ER.Mechanisms.Node) == (:conformation, :bound, :added, :subtracted)
    @test fieldnames(ER.Mechanisms.Step) ==
          (:from, :to, :taken_up, :given_off, :rapid_equilibrium, :group)
    @test fieldnames(ER.AllostericMechanism) ==
          (:catalytic, :mwc_tags, :catalytic_multiplicity, :regulatory_sites, :active_steps,
           :inactive_steps, :excluded_active, :excluded_inactive)
end

@testset "construction is canonical" begin
    # Michaelis–Menten from nodes and steps: (a) in canonical form, (b) its nodes permuted
    # and one written twice, (c) every step reversed, (d) its group labels renumbered.
    # Every writing builds one mechanism, with one hash and one printed text.
    M = ER.Mechanisms
    S, P, none = [M.Met(:S, false)], [M.Met(:P, false)], M.Met[]
    E, EP, ES = M.Node(:E, none, Symbol[], Symbol[]), M.Node(:E, P, Symbol[], Symbol[]),
                M.Node(:E, S, Symbol[], Symbol[])
    steps = [M.Step(1, 2, P, none, true, 1), M.Step(1, 3, S, none, true, 2),
             M.Step(3, 2, none, none, false, 3)]
    canonical = ER.Mechanism([:S], [:P], [E, EP, ES], steps)
    @test (canonical.nodes, canonical.steps) == ([E, EP, ES], steps)
    @test canonical == ER.Mechanism(@enzyme_mechanism begin
        substrates: S
        products: P
        steps: begin
            E + S ⇌ E(S)
            E(S) <--> E(P)
            E(P) ⇌ E + P
        end
    end)
    permuted = ER.Mechanism([:S], [:P], [ES, E, EP, E],
                            [M.Step(4, 3, P, none, true, 1), M.Step(2, 1, S, none, true, 2),
                             M.Step(1, 3, none, none, false, 3)])
    reversed = ER.Mechanism([:S], [:P], [E, EP, ES],
                            [M.Step(2, 1, none, P, true, 1), M.Step(3, 1, none, S, true, 2),
                             M.Step(2, 3, none, none, false, 3)])
    relabelled_steps = [M.Step(1, 2, P, none, true, 3), M.Step(1, 3, S, none, true, 1),
                        M.Step(3, 2, none, none, false, 2)]
    relabelled = ER.Mechanism([:S], [:P], [E, EP, ES], relabelled_steps)
    for m in (permuted, reversed, relabelled)
        @test m == canonical && hash(m) == hash(canonical)
        @test sprint(show, m) == sprint(show, canonical)
    end
    # Each list of a node or a step in any order: ordered bi-uni with E(A, B) written as
    # E(B, A).
    A, B, Q = M.Met(:A, false), M.Met(:B, false), M.Met(:Q, false)
    EA, EAB = M.Node(:E, [A], Symbol[], Symbol[]), M.Node(:E, [A, B], Symbol[], Symbol[])
    _testhelper_bi_uni(eab) = ER.Mechanism(
        [:A, :B], [:Q], [M.Node(:E, none, Symbol[], Symbol[]), EA, eab],
        [M.Step(1, 2, [A], none, true, 1), M.Step(2, 3, [B], none, true, 2),
         M.Step(3, 1, none, [Q], false, 3)])
    unsorted = _testhelper_bi_uni(M.Node(:E, [B, A], Symbol[], Symbol[]))
    @test unsorted == _testhelper_bi_uni(EAB) && EAB in unsorted.nodes
    # An allosteric mechanism from its written parts takes each tag by written group label:
    # label 1 (S's binding) is canonical group 2, label 2 (the chemistry) group 3 and
    # label 3 (P's binding) group 1.
    site = M.RegulatorySite([:R], 2, [:OnlyI])
    am = ER.AllostericMechanism([:S], [:P], [E, EP, ES], relabelled_steps,
                                [:EqualAI, :OnlyA, :NonequalAI], 2, [site])
    expected = ER.AllostericMechanism(canonical, [:NonequalAI, :EqualAI, :OnlyA], 2, [site])
    @test am == expected && hash(am) == hash(expected)
    @test sprint(show, am) == sprint(show, expected)
end

@testset "an allosteric graph sorts its sites and checks its tags and ligands" begin
    m = _testhelper_graph_uniuni()
    tags = [:EqualAI, :EqualAI, :OnlyA]
    sites = [EnzymeRates.Mechanisms.RegulatorySite([:X], 2, [:OnlyI]),
             EnzymeRates.Mechanisms.RegulatorySite([:B, :A], 2, [:OnlyA, :EqualAI])]
    am = ER.AllostericMechanism(m, tags, 2, sites)
    @test [s.ligands for s in am.regulatory_sites] == [[:A, :B], [:X]]
    @test am.regulatory_sites[1].mwc_tags == [:EqualAI, :OnlyA]
    tags[1] = :OnlyA
    @test am.mwc_tags == [:EqualAI, :EqualAI, :OnlyA]
    @test [am.mwc_tags[g] for g in 1:3] == [:EqualAI, :EqualAI, :OnlyA]
    @test_throws(ErrorException("AllostericMechanism: mwc_tags length 2 must " *
                                "match the number of kinetic groups 3"),
                 ER.AllostericMechanism(m, [:EqualAI, :EqualAI], 2, sites))
    @test_throws(ErrorException("AllostericMechanism: catalytic_multiplicity must be ≥ " *
                                "1, got 0"),
                 ER.AllostericMechanism(m, [:EqualAI, :EqualAI, :OnlyA], 0, sites))
    @test_throws(ErrorException(
        "AllostericMechanism: E + S ⇌ E(S) is tagged :OnlyI but does more than bind a " *
        "modifier. Only the active conformation catalyzes, so :OnlyI fits only a group " *
        "whose steps take up or give off modifiers alone. Tag its group :NonequalAI."),
        ER.AllostericMechanism(m, [:EqualAI, :OnlyI, :OnlyA], 2, sites))
    @test_throws(ErrorException(
        "AllostericMechanism: ligand X appears in two distinct regulatory sites; " *
        "rendered Kreg names would collide. Same-ligand-two-sites is not enumerated."),
                 ER.AllostericMechanism(
                     m, [:EqualAI, :EqualAI, :OnlyA], 2,
                     [sites; EnzymeRates.Mechanisms.RegulatorySite([:X], 4, [:OnlyA])]))
    @test_throws(ErrorException(
        "AllostericMechanism: L is a reserved name; E_total, Keq, L, N, params, concs " *
        "and v are kept for a rate equation's constants and local names"),
        ER.AllostericMechanism(m, [:EqualAI, :EqualAI, :OnlyA], 2,
                               [EnzymeRates.Mechanisms.RegulatorySite([:L], 2, [:OnlyI])]))
end

"""`raw`, a spec's raw value (`raw_mechanism`), written another way, the same mechanism:
with `order`, its nodes and its steps in reverse order; with `direction`, each step from
its other end; with `labels`, its group labels numbered from the last. Its declared
names come in reverse order. An allosteric value's tags follow their labels, and its
sites and each site's ligands come in reverse order."""
function _testhelper_rewrite(raw::NamedTuple; order = true, direction = true,
                             labels = true)
    M = EnzymeRates.Mechanisms
    nodes, steps, tags = raw.nodes, raw.steps, get(raw, :mwc_tags, Symbol[])
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
    plain = (; substrates = reverse(raw.substrates), products = reverse(raw.products),
             modifiers = reverse(raw.modifiers), nodes, steps)
    haskey(raw, :mwc_tags) || return plain
    (; plain..., mwc_tags = tags, catalytic_multiplicity = raw.catalytic_multiplicity,
     regulatory_sites = [(; ligands = reverse(s.ligands), s.multiplicity,
                          mwc_tags = reverse(s.mwc_tags))
                         for s in reverse(raw.regulatory_sites)])
end

@testset "every writing of a mechanism gives one canonical graph" begin
    @testset "$(spec.name)" for spec in MECHANISM_TEST_SPECS
        m = _testhelper_canonicalize(spec.raw).mechanism
        for (order, direction, labels) in ((true, false, false), (false, true, false),
                                           (false, false, true), (true, true, true))
            other = _testhelper_canonicalize(
                _testhelper_rewrite(spec.raw; order, direction, labels)).mechanism
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
    end), false))
    # Forms E, E(P), E(S); the release of P is stored as the binding it reverses.
    @test c.mechanism == _testhelper_graph_uniuni()
    @test c.edge_of == [2, 3, 1]
    @test c.flipped == [false, false, true]
    @test c.group_of == [2, 3, 1]
    @test c.node_of == [1, 3, 3, 2, 2, 1]
    @test isempty(c.misaligned)
end

@testset "canonical form runs each step physical-forward, tied steps as one" begin
    _testhelper_flipped(block) = _testhelper_canonicalize(
        EnzymeRates.DSL.raw_mechanism(block, false)).flipped
    # A step with one metabolite takes it up; a fused release F(P) → F + P too. Then the
    # side with more substrates comes first (E(S) → E(P)); then the form where only
    # products enter or leave (F → E); then the earlier form, F(S) → E(S) running E(S)
    # → F(S), E(S) and F(S) seeing no metabolite enter or leave.
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
    # The mirror of F → E at the I-bound forms would run E(I) → F(I) on its own (the
    # earlier form); in one kinetic group with F → E it runs F(I) → E(I). The release
    # of P is stored as the binding it reverses.
    _testhelper_mirrored(iso) = EnzymeRates.DSL.raw_mechanism(:(begin
        substrates: A
        products: P
        modifiers: I
        steps: begin
            E + A <--> E(A)
            E(A) <--> E(P)
            E(P) <--> F + P
            (F <--> E, $iso)
            (E + I ⇌ E(I), F + I ⇌ F(I))
        end
    end), false)
    c, d = _testhelper_canonicalize(_testhelper_mirrored(:(F(I) <--> E(I)))),
           _testhelper_canonicalize(_testhelper_mirrored(:(E(I) <--> F(I))))
    @test c.flipped == [false, false, true, false, false, false, false]
    @test d.flipped == [false, false, true, false, true, false, false]
    @test c.mechanism == d.mechanism
    # Tied steps that change different pairs of conformations: no orientation aligns
    # them, and canonical form returns the lead and the first that misses it.
    c = _testhelper_canonicalize(EnzymeRates.DSL.raw_mechanism(:(begin
        substrates: A
        products: P
        steps: begin
            E + A <--> E(A)
            E(A) <--> E(P)
            E(P) <--> E + P
            (E ⇌ F, F ⇌ G, G ⇌ E)
        end
    end), false))
    none = EnzymeRates.Mechanisms.Met[]
    @test c.misaligned == [ER.Mechanisms.Step(1, 4, none, none, true, 4),
                           ER.Mechanisms.Step(4, 5, none, none, true, 4)]
end

"""Bindings of each metabolite of `mets` to free E, node 1, each into a node of its own and
in a group of its own, at rapid equilibrium: the nodes and the steps."""
function _testhelper_bindings(mets::Vector{ER.Mechanisms.Met})
    M = ER.Mechanisms
    nodes = [M.Node(:E, M.Met[], Symbol[], Symbol[]);
             [M.Node(:E, [x], Symbol[], Symbol[]) for x in mets]]
    nodes, [M.Step(1, i + 1, [x], M.Met[], true, i) for (i, x) in enumerate(mets)]
end

@testset "a graph holds any number of metabolites and of copies of one" begin
    M = ER.Mechanisms
    names = [Symbol("S", i) for i in 1:32]
    plain(xs) = [M.Met(x, false) for x in xs]
    _testhelper_entries(substrates, mets) = length(M._table(
        M._canonicalize(substrates, [:P], _testhelper_bindings(mets)...).mechanism))
    @test _testhelper_entries(names[1:31], plain([names[1:31]; :P])) == 32
    @test _testhelper_entries(names, plain([names; :P])) == 33
    # A reactant and its modifier copy are two entries.
    @test _testhelper_entries(names[1:16], [plain([names[1:16]; :P]);
                                            [M.Met(x, true) for x in names[1:16]]]) == 33
    _testhelper_holding(n) = M._canonicalize(
        [:A], [:P],
        [M.Node(c, fill(M.Met(:A, false), n), Symbol[], Symbol[]) for c in (:E, :F)],
        [M.Step(1, 2, M.Met[], M.Met[], false, 1)])
    for n in (15, 16)
        m = _testhelper_holding(n).mechanism
        @test count(x -> x.name === :A, m.nodes[1].bound) == n
    end
end

"""The raw value of `block`, a mechanism macro's block, read by the DSL walker."""
_testhelper_walk(block; allosteric = false) =
    EnzymeRates.DSL.raw_mechanism(block, allosteric)

@testset "the constructors build one value from every writing of a mechanism" begin
    @testset "$(spec.name)" for spec in MECHANISM_TEST_SPECS
        x = _testhelper_build(spec.raw)
        @test ER.Mechanisms.catalytic(x) == _testhelper_canonicalize(spec.raw).mechanism
        other = _testhelper_build(_testhelper_rewrite(spec.raw))
        @test other == x && hash(other) == hash(x)
    end
    # Each tag goes to the canonical group its label became.
    am = _testhelper_build(_testhelper_walk(:(begin
        substrates: S
        products: P
        catalytic_multiplicity: 2
        steps: begin
            E + S ⇌ E(S) :: EqualAI
            E(S) <--> E(P) :: OnlyA
            E(P) ⇌ E + P :: NonequalAI
        end
    end); allosteric = true))
    @test am.mwc_tags == [:NonequalAI, :EqualAI, :OnlyA]
end

"""Whether `v` takes the memory of its copy, which `copy` stores at its length."""
_testhelper_at_length(v) = Base.summarysize(v) == Base.summarysize(copy(v))

"""Whether `m` stores each of its vectors at its length: `m` takes the memory of the
value rebuilt from a copy of each vector, its nodes' and steps' lists included, which
`copy` stores at its length, and each fact it stores takes the memory of its copy."""
_testhelper_stored_at_length(m::ER.Mechanism) = Base.summarysize(m) == Base.summarysize(
    ER.Mechanisms._unchecked_mechanism(copy(m.substrates), copy(m.products),
                 [ER.Mechanisms.Node(n.conformation, copy(n.bound), copy(n.added),
                                     copy(n.subtracted)) for n in m.nodes],
                 [ER.Mechanisms.Step(e.from, e.to, copy(e.taken_up), copy(e.given_off),
                                     e.rapid_equilibrium, e.group) for e in m.steps])) &&
    all(_testhelper_at_length, (m.modifiers, m.groups, m.step_kinds, m.toward_products,
                                m.step_sides)) && all(_testhelper_at_length, m.groups)
_testhelper_stored_at_length(am::ER.AllostericMechanism) =
    _testhelper_stored_at_length(am.catalytic) && Base.summarysize(am) ==
    Base.summarysize(ER.AllostericMechanism(am.catalytic, copy(am.mwc_tags),
        am.catalytic_multiplicity,
        [ER.Mechanisms.RegulatorySite(copy(s.ligands), s.multiplicity, copy(s.mwc_tags))
         for s in am.regulatory_sites])) &&
    all(_testhelper_at_length, (am.active_steps, am.inactive_steps, am.excluded_active,
                                am.excluded_inactive))

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
    init = EnzymeRates.Enumeration.init_mechanisms(rxn)
    children = EnzymeRates.Enumeration.expand_mechanisms(init[1:10], rxn)
    @test any(c -> c isa ER.AllostericMechanism, children)
    @test all(_testhelper_stored_at_length, init)
    @test all(_testhelper_stored_at_length, children)
end

@testset "the constructors reject what no graph represents, each with its message" begin
    # Uni-uni, E + S ⇌ E(S), E(S) <--> E(P), E(P) ⇌ E + P, nodes 1 to 3, with the nodes
    # `extra` (4, 5, …) and `steps` added; every name a node or a step binds as a modifier
    # is a modifier.
    M = ER.Mechanisms
    S, P, I, Imod, none = M.Met(:S, false), M.Met(:P, false), M.Met(:I, false),
                          M.Met(:I, true), M.Met[]
    node(c, bound...) = M.Node(c, collect(M.Met, bound), Symbol[], Symbol[])
    nodes = [node(:E), node(:E, S), node(:E, P)]
    uniuni = [M.Step(1, 2, [S], none, true, 1), M.Step(2, 3, none, none, false, 2),
              M.Step(3, 1, none, [P], true, 3)]
    _testhelper_uniuni(extra = M.Node[], steps...; substrates = [:S], products = [:P],
                       group = [1, 2, 3]) =
        ER.Mechanism(substrates, products, [nodes; extra],
                     [[M.Step(e.from, e.to, e.taken_up, e.given_off, e.rapid_equilibrium, g)
                       for (e, g) in zip(uniuni, group)]; steps...])
    # Uni-uni with `x` bound to free E by a step of its own.
    _testhelper_binding(x; kw...) =
        _testhelper_uniuni([node(:E, x)], M.Step(1, 4, [x], none, true, 4); kw...)
    _testhelper_allosteric(tags, n, sites) =
        ER.AllostericMechanism([:S], [:P], nodes, uniuni, tags, n, sites)
    tags = [:EqualAI, :OnlyA, :EqualAI]
    site = M.RegulatorySite
    @test _testhelper_uniuni() == _testhelper_graph_uniuni()
    # Regulatory sites.
    @test_throws(ErrorException("RegulatorySite: needs at least one ligand"),
                 _testhelper_allosteric(tags, 2, [site(Symbol[], 2, Symbol[])]))
    @test_throws(ErrorException("RegulatorySite: length(ligands)=1 must equal " *
                                "length(mwc_tags)=0"),
                 _testhelper_allosteric(tags, 2, [site([:R], 2, Symbol[])]))
    @test_throws(ErrorException("RegulatorySite: multiplicity must be ≥ 1, got 0"),
                 _testhelper_allosteric(tags, 2, [site([:R], 0, [:OnlyA])]))
    @test_throws(ErrorException("RegulatorySite: MWC tag Bogus must be one of " *
                                "(:OnlyA, :OnlyI, :EqualAI, :NonequalAI)"),
                 _testhelper_allosteric(tags, 2, [site([:R], 2, [:Bogus])]))
    @test_throws(ErrorException("RegulatorySite: ligand R appears twice in one " *
                                "regulatory site"),
                 _testhelper_allosteric(tags, 2, [site([:R, :X, :R], 2,
                                                       [:OnlyA, :OnlyI, :OnlyI])]))
    @test_throws(ErrorException(
        "AllostericMechanism: ligand R appears in two distinct regulatory sites; " *
        "rendered Kreg names would collide. Same-ligand-two-sites is not enumerated."),
        _testhelper_allosteric(tags, 2, [site([:R], 2, [:OnlyA]), site([:R], 4, [:OnlyI])]))
    # A ligand takes no name a metabolite cannot take.
    @test_throws(ErrorException("AllostericMechanism: \"1R\" cannot name a metabolite; a " *
                                "metabolite name is an identifier that is no keyword"),
                 _testhelper_allosteric(tags, 2, [site([Symbol("1R")], 2, [:OnlyA])]))
    @test_throws(ErrorException("AllostericMechanism: E cannot name a metabolite; it " *
                                "labels the free enzyme form"),
                 _testhelper_allosteric(tags, 2, [site([:E], 2, [:OnlyA])]))
    # Tags, by written group label: label 2 is canonical group 3.
    @test_throws(ErrorException("AllostericMechanism: mwc_tags length 2 must " *
                                "match the number of kinetic groups 3"),
                 _testhelper_allosteric(tags[1:2], 2, M.RegulatorySite[]))
    @test_throws(ErrorException("AllostericMechanism: catalytic_multiplicity must be ≥ " *
                                "1, got 0"),
                 _testhelper_allosteric(tags, 0, M.RegulatorySite[]))
    # A Bool is no multiplicity: each constructor takes an Int and converts no Bool.
    @test_throws MethodError site([:R], true, [:OnlyA])
    @test_throws MethodError _testhelper_allosteric(tags, true, M.RegulatorySite[])
    @test_throws MethodError ER.AllostericMechanism(_testhelper_uniuni(), tags, true,
                                                    M.RegulatorySite[])
    @test_throws(ErrorException(
        "AllostericMechanism: E(S) <--> E(P) is tagged :OnlyI but does more than bind a " *
        "modifier. Only the active conformation catalyzes, so :OnlyI fits only a group " *
        "whose steps take up or give off modifiers alone. Tag its group :NonequalAI."),
                 _testhelper_allosteric([:EqualAI, :OnlyI, :EqualAI], 2,
                                        M.RegulatorySite[]))
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
                     _testhelper_uniuni(M.Node[], M.Step(from, to, none, none, false, 4)))
    end
    # Declared and bound names a metabolite cannot take: no identifier, a keyword, E, a
    # reserved name, a name holding `_`.
    _testhelper_unwritable(x) = ErrorException(
        "Mechanism: \"$x\" cannot name a metabolite; a metabolite name is an identifier " *
        "that is no keyword")
    @test_throws(_testhelper_unwritable("1I"),
                 _testhelper_binding(M.Met(Symbol("1I"), true)))
    @test_throws(_testhelper_unwritable("end"), _testhelper_uniuni(products = [:P, :end]))
    @test_throws(ErrorException("Mechanism: E cannot name a metabolite; it labels the " *
                                "free enzyme form"),
                 _testhelper_uniuni(substrates = [:S, :E]))
    @test_throws(ErrorException("Mechanism: Keq is a reserved name; E_total, Keq, L, N, " *
                                "params, concs and v are kept for a rate equation's " *
                                "constants and local names"),
                 _testhelper_uniuni(substrates = [:S, :Keq]))
    @test_throws(ErrorException("Mechanism: the metabolite name S_1 holds `_`; constant " *
                                "names join metabolite names with `_`, so a metabolite " *
                                "name holds none"),
                 _testhelper_binding(M.Met(:S_1, true)))
    # A name Julia reads as another (a micro sign as μ, a decomposed é as one character)
    # prints as text that builds another mechanism, so no DSL text can write it. The
    # strings are written with escapes: Julia would read the characters themselves as
    # the other names.
    _testhelper_read_as(who, x, read) = ErrorException(
        "$who: Julia reads the name \"$x\" as \"$read\", so no DSL text can write it; " *
        "use \"$read\"")
    @test_throws(_testhelper_read_as("Mechanism", "\u00b5M", "\u03bcM"),
                 _testhelper_binding(M.Met(Symbol("\u00b5M"), true)))
    @test_throws(_testhelper_read_as("Mechanism", "e\u0301", "\u00e9"),
                 _testhelper_uniuni(substrates = [:S, Symbol("e\u0301")]))
    @test_throws(_testhelper_read_as("AllostericMechanism", "\u025b", "\u03b5"),
                 _testhelper_allosteric(tags, 2, [site([Symbol("\u025b")], 2, [:OnlyA])]))
    # Names: undeclared; a reactant's modifier copy; a plain entry of a name bound as a
    # modifier; residuals by role.
    @test_throws(ErrorException("Mechanism: P is not declared; declare it under " *
                                "substrates:, products: or modifiers:"),
                 _testhelper_uniuni(products = [:Q]))
    @test _testhelper_binding(M.Met(:P, true)) isa ER.Mechanism
    @test _testhelper_binding(Imod) isa ER.Mechanism
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
            _testhelper_uniuni([M.Node(:E, none, added, subtracted)],
                               M.Step(2, 4, none, none, false, 4)))
    end
    # Steps: a self-loop, a metabolite taken up and given off, a role change.
    @test_throws(ErrorException("Step: both ends are the form EA; a step joins two forms"),
                 ER.Mechanism([:A], [:P], [node(:E, M.Met(:A, false))],
                              [M.Step(1, 1, none, none, false, 1)]))
    @test_throws(ErrorException("Step: both ends are the form ES; a step joins two forms"),
                 _testhelper_uniuni(M.Node[], M.Step(2, 2, none, none, false, 4)))
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
    # differ.
    @test_throws(ErrorException("Step EA → EQ: B both consumed and released"),
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
    # Tied steps that make different conformational changes.
    @test_throws(ErrorException(
        "Mechanism: a kinetic group holds E → F and F → G, which make different " *
        "conformational changes; the steps of a kinetic group share their constants, so " *
        "they must make the same conformational change"),
        _testhelper_build(_testhelper_walk(:(begin
            substrates: A
            products: P
            steps: begin
                E + A <--> E(A)
                E(A) <--> E(P)
                E(P) <--> E + P
                (E ⇌ F, F ⇌ G, G ⇌ E)
            end
        end))))
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
    # The tied pair stores E → F in group 3, the lone step F → E in group 6.
    @test_throws(ErrorException(
        "Mechanism: kinetic groups 3 and 6 both hold the reaction E ⇌ F; a reaction " *
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
    @test_throws no_free_E ER.Mechanism([:S], [:P], M.Node[], M.Step[])
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
                 _testhelper_binding(M.Met(:T, true); substrates = [:S, :T]))
    for x in (Imod, M.Met(:P, true))
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
    @test_throws(ErrorException("AllostericMechanism: conformation label `F` collides " *
                                "with ligand `F`; choose a different conformation " *
                                "label."),
                 ER.AllostericMechanism([:S], [:P], [nodes; node(:F)],
                                        [uniuni; M.Step(1, 4, none, none, false, 4)],
                                        [tags; :EqualAI], 2, [site([:F], 2, [:OnlyA])]))
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
    @test [ER.Mechanisms._form_name(m, f) for f in eachindex(m.nodes)] ==
          ["E", "E_res_plusA_minusP", "EA", "EB_res_plusA_minusP", "EI", "EQ", "EQQmod"]
    @test [ER.Mechanisms.form_text(m, f) for f in eachindex(m.nodes)] ==
          ["E", "E(; residual = A - P)", "E(A)", "E(B; residual = A - P)", "E(I)",
           "E(Q)", "E(Q, Q::Mod)"]
    # Each edge's sides in its stored direction: releases stored as bindings.
    @test [join(m.step_sides[k], " → ") for k in eachindex(m.steps)] ==
          ["E_A → EA", "E_I → EI", "E_Q → EQ", "E_res_plusA_minusP_P → EA",
           "E_res_plusA_minusP_B → EB_res_plusA_minusP", "EB_res_plusA_minusP → EQ",
           "EQ_Qmod → EQQmod"]
end

"""The graph whose Sig is `sig`, `x`'s kind of graph. `sig` is splatted with its length
left open, as the lifts splat it, so the decoders compile once, whatever the Sig length."""
_testhelper_decode(::ER.Mechanism, @nospecialize(sig)) =
    ER.Mechanisms._mechanism_from_sig(sig...)
_testhelper_decode(::ER.AllostericMechanism, @nospecialize(sig)) =
    ER.Mechanisms._allosteric_mechanism_from_sig(sig...)

@testset "a graph's Sig decodes to the graph, and equal graphs share one type" begin
    @testset "$(spec.name)" for spec in MECHANISM_TEST_SPECS
        x = _testhelper_build(spec.raw)
        sig = ER.Mechanisms._sig_of(x)
        @test sig isa Tuple{Vararg{Int}}
        @test _testhelper_decode(x, sig) == x
        cx = ER.compile_mechanism(x)
        @test typeof(cx).parameters[1] === sig
        @test ER.Mechanisms._concrete(cx) == x
        @test ER.compile_mechanism(_testhelper_build(_testhelper_rewrite(spec.raw))) === cx
    end
    m = _testhelper_graph_uniuni()
    @test ER.Mechanism(ER.EnzymeMechanism(m)) == m
    am = ER.AllostericMechanism(
        m, [:EqualAI, :NonequalAI, :OnlyA], 4,
        [EnzymeRates.Mechanisms.RegulatorySite([:X, :R], 2, [:OnlyA, :OnlyI])])
    aem = ER.AllostericEnzymeMechanism(am)
    @test ER.AllostericMechanism(aem) == am
    # A name outside ASCII, sorting last among 32 metabolite entries.
    names = [[Symbol("S", i) for i in 1:30]; :αKG]
    wide = ER.Mechanism(names, [:P], _testhelper_bindings(
        [EnzymeRates.Mechanisms.Met(x, false) for x in [names; :P]])...)
    @test ER.Mechanisms._table(wide)[32] == EnzymeRates.Mechanisms.Met(:αKG, false)
    @test _testhelper_decode(wide, ER.Mechanisms._sig_of(wide)) == wide
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

@testset "the Sig decoders compile once, whatever a Sig's length" begin
    # The lifts splat a Sig whose length their code leaves open, so each varargs decoder
    # compiles once, whatever the Sig length: after one lift, lifting mechanisms of three
    # Sig lengths no other test decodes adds no method instance. A decoder that took its
    # Sig as an `NTuple{N, Int}` would compile once per length.
    instances(f) = count(_ -> true, Base.specializations(only(methods(f))))
    _testhelper_allosteric_long_named(n) =
        ER.AllostericMechanism(_testhelper_long_named(n), [:EqualAI, :EqualAI, :OnlyA], 2,
                               EnzymeRates.Mechanisms.RegulatorySite[])
    for (decoder, build) in ((ER.Mechanisms._mechanism_from_sig, _testhelper_long_named),
                             (ER.Mechanisms._allosteric_mechanism_from_sig,
                              _testhelper_allosteric_long_named))
        lengths = (301, 401, 501, 601)
        compiled = [ER.compile_mechanism(build(n)) for n in lengths]
        @test ER.Mechanisms._concrete(compiled[1]) == build(lengths[1])
        before = instances(decoder)
        foreach(ER.Mechanisms._concrete, compiled[2:4])
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
        @test c.group_of == 1:length(ER.Mechanisms.catalytic(x).groups)
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
    population = Union{EnzymeRates.Mechanism, EnzymeRates.AllostericMechanism}[]
    foreach(r -> append!(population, EnzymeRates.Enumeration.init_mechanisms(r)), reactions)
    bibi = reactions[2]
    append!(population,
            EnzymeRates.Enumeration.expand_mechanisms(
                EnzymeRates.Enumeration.init_mechanisms(bibi), bibi))
end

"""The census's first allosteric mechanisms: the allosteric children of the init
mechanisms of its three reactions that declare an allosteric regulator."""
function _testhelper_census_allosteric()
    children = [EnzymeRates.Enumeration.expand_mechanisms(
        EnzymeRates.Enumeration.init_mechanisms(r), r)
                for r in _testhelper_census_reactions()[[1, 3, 4]]]
    filter(x -> x isa EnzymeRates.AllostericMechanism, reduce(vcat, children))
end

"""Graph `x`'s canonical form as text lines: its names, its forms in order, each edge in
its stored direction (`_testhelper_edge_text`) with its group (and the group's tag), then
any catalytic multiplicity and sites."""
function _testhelper_canonical_text(x)
    m, allo = ER.Mechanisms.catalytic(x), x isa ER.AllostericMechanism
    mods = m.modifiers
    lines = ["substrates: $(join(m.substrates, ", ")) | products: " *
             "$(join(m.products, ", "))" *
             (isempty(mods) ? "" : " | modifiers: $(join(mods, ", "))"),
             "forms: " *
             join([ER.Mechanisms._form_name(m, f) for f in eachindex(m.nodes)], ", ")]
    for k in eachindex(m.steps)
        g = m.steps[k].group
        push!(lines, "$k: $(_testhelper_edge_text(m, k)) | group $g" *
                     (allo ? " :: $(x.mwc_tags[g])" : ""))
    end
    allo && push!(lines, "catalytic_multiplicity: $(x.catalytic_multiplicity)")
    for s in (allo ? x.regulatory_sites : EnzymeRates.Mechanisms.RegulatorySite[])
        push!(lines, "regulatory_site($(s.multiplicity)): " *
                     join(["$l::$t" for (l, t) in zip(s.ligands, s.mwc_tags)], ", "))
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
Whether the sub-graph reads of `x` (`re_segments`, `conformation_steps`) give what the
oracles read off its test reader's steps give (`_testhelper_re_segments`,
`_testhelper_inactive_reach`): the catalytic graph's segments, form by form and
breadth-first, with each form's positive offsets, and, for an allosteric mechanism, the
inactive conformation's edges and their segments.
"""
function _testhelper_subgraphs_agree(x)
    m, g = ER.Mechanisms.catalytic(x), _testhelper_graph(x)
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
    x isa ER.AllostericMechanism || return true
    inactive = x.inactive_steps
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
        @test ER.Mechanisms._concrete(spec.mechanism) == _testhelper_build(spec.raw)
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
    none = EnzymeRates.Mechanisms.Met[]
    @test x.nodes[1:2] == [ER.Mechanisms.Node(:E, none, Symbol[], Symbol[]),
                           ER.Mechanisms.Node(:D, none, Symbol[], Symbol[])]
    @test _testhelper_graph(x).forms[1:2] == ["E", "D"]
    # The census's first populations: each graph decodes from its Sig and builds again
    # from the nodes and steps it writes back.
    @test length(population) == 5508 + 5505
    @test all(_testhelper_decode(x, ER.Mechanisms._sig_of(x)) == x for x in population)
    @test all(_testhelper_build(_testhelper_written(x)) == x for x in population)
end

@testset "a sub-graph's segments and inactive edges are what its steps give" begin
    @testset "$(spec.name)" for spec in MECHANISM_TEST_SPECS
        @test _testhelper_subgraphs_agree(_testhelper_lift(spec.mechanism))
    end
    @test isempty([i for (i, x) in enumerate(population)
                   if !_testhelper_subgraphs_agree(x)])
end

@testset "a graph's DSL text reads back as the graph" begin
    # `show` prints the macro text that builds a mechanism and `DSL.parse_mechanism`
    # reads it back, for the specs, plain and allosteric, and the census's first
    # populations, plain and allosteric.
    @testset "$(spec.name)" for spec in MECHANISM_TEST_SPECS
        @test _testhelper_reads_back(_testhelper_lift(spec.mechanism))
        @test ER.compile_mechanism(
            ER.DSL.parse_mechanism(sprint(show, spec.mechanism))) === spec.mechanism
    end
    @test length(allosteric) == 5 + 890 + 5
    @test isempty([i for (i, x) in enumerate(population) if !_testhelper_reads_back(x)])
    @test isempty([i for (i, x) in enumerate(allosteric) if !_testhelper_reads_back(x)])
end
end

@testset "tier 3 and a tied group's lead orient steps alike in every writing" begin
    # No spec and no census population holds a step that only form order orients or a
    # tied group whose alignment turns a step, so these cases are checked here, each in
    # two writings, against each case's pinned canonical text.
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
                 "5: EJ + I <--> EI + J | group 5",
                 "6: EJ + P ⇌ EJP | group 6",
                 "7: F <--> E | group 7",
                 "8: F + P <--> EP | group 8",
                 "9: F + I ⇌ FI | group 2",
                 "10: FI <--> EI | group 7",
                 "11: FJ + I <--> FI + J | group 5"]
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
                       "catalytic_multiplicity: 1"]
    _testhelper_tied_steps(iso, exchange) = :(begin
        E + A <--> E(A)
        E(A) <--> E(P)
        E(P) <--> F + P
        ($iso, F(I) <--> E(I))
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
    for (iso, exchange) in (
            (:(F <--> E), :((E(I) + J <--> E(J) + I, F(I) + J <--> F(J) + I))),
            (:(E <--> F), :((F(J) + I <--> F(I) + J, E(J) + I <--> E(I) + J))))
        push!(blocks, :(begin
            substrates: A
            products: P
            modifiers: I, J
            steps: $(_testhelper_tied_steps(iso, exchange))
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
        end); allosteric = true)
        @test _testhelper_canonical_text(_testhelper_build(raw)) == allosteric_text
    end
end

"""Whether `x` survives its round trips: its printed text reads back as `x` and prints
again as that text, two reads of the text hash alike, and `x` decodes from its compiled
type as itself."""
function _testhelper_round_trips(x)
    text = sprint(show, x)
    a, b = ER.DSL.parse_mechanism(text), ER.DSL.parse_mechanism(text)
    a == x && sprint(show, a) == text && hash(a) == hash(b) == hash(x) &&
        _testhelper_lift(ER.compile_mechanism(x)) == x
end

@testset "storage round trips" begin
    # The specs, plain and allosteric, and the bi-bi init set, the enumeration golden set
    # (`phase2_init_golden.txt`).
    rxn = @enzyme_reaction begin
        substrates: A[C], B[N]
        products: P[C], Q[N]
    end
    @testset "$(spec.name)" for spec in MECHANISM_TEST_SPECS
        @test _testhelper_round_trips(_testhelper_lift(spec.mechanism))
    end
    init = EnzymeRates.Enumeration.init_mechanisms(rxn)
    @test length(init) == 239
    @test isempty([i for (i, x) in enumerate(init) if !_testhelper_round_trips(x)])
end

"""Whether the facts `m` stores are those its canonical fields give, and the functions that
read them return them: its modifiers and groups, read off its nodes and steps here, and
each step's kind, direction toward products and side names, as the graph reads compute
them from the canonical fields."""
function _testhelper_facts_hold(m::ER.Mechanism)
    M = ER.Mechanisms
    lists = [[n.bound for n in m.nodes]; [e.taken_up for e in m.steps];
             [e.given_off for e in m.steps]]
    modifiers = sort!(unique!(Symbol[x.name for list in lists for x in list if x.modifier]))
    groups = [findall(e -> e.group == g, m.steps)
              for g in 1:maximum(e -> e.group, m.steps)]
    m.modifiers == modifiers && m.groups == groups &&
        m.step_kinds == [M._step_kind(m, e) for e in m.steps] &&
        m.toward_products == [M._runs_toward_products(m, e) for e in m.steps] &&
        m.step_sides == [M._sides(m, e) for e in m.steps]
end

"""Whether the facts `am` stores are those its fields give: its catalytic graph's
(`_testhelper_facts_hold`), and per conformation the steps its tags exclude
(`_testhelper_excludes`) and the steps it keeps on the forms they join to free E
(`_testhelper_exclusion`)."""
function _testhelper_facts_hold(am::ER.AllostericMechanism)
    cm = am.catalytic
    _testhelper_facts_hold(cm) || return false
    all(((1, am.excluded_active, am.active_steps),
         (2, am.excluded_inactive, am.inactive_steps))) do (c, out, kept)
        excluded = BitVector([_testhelper_excludes(am, k, c) for k in eachindex(cm.steps)])
        component = _testhelper_exclusion(am, c).component
        steps = BitVector([!excluded[k] && component[e.from] == 1
                           for (k, e) in enumerate(cm.steps)])
        (out, kept) == (excluded, steps)
    end
end

"""Whether `a` and `b` hold equal values in every field, a catalytic graph's included."""
_testhelper_all_fields_equal(a, b) = all(fieldnames(typeof(a))) do f
    x, y = getfield(a, f), getfield(b, f)
    x isa ER.Mechanism ? _testhelper_all_fields_equal(x, y) : x == y
end

@testset "stored facts equal what the canonical fields give" begin
    # Every construction path: the constructors (the specs' writings, the bi-bi seeds)
    # and the Sig decoder (the specs' compiled forms).
    @testset "$(spec.name)" for spec in MECHANISM_TEST_SPECS
        @test _testhelper_facts_hold(_testhelper_build(spec.raw))
        @test _testhelper_facts_hold(_testhelper_lift(spec.mechanism))
    end
    rxn = @enzyme_reaction begin
        substrates: A[C], B[N]
        products: P[C], Q[N]
    end
    @test all(_testhelper_facts_hold, EnzymeRates.Enumeration.init_mechanisms(rxn))
end

@testset "stored facts survive the singleton round trip" begin
    @testset "$(spec.name)" for spec in MECHANISM_TEST_SPECS
        x = _testhelper_build(spec.raw)
        @test _testhelper_all_fields_equal(_testhelper_lift(ER.compile_mechanism(x)), x)
    end
    rxn = @enzyme_reaction begin
        substrates: A[C], B[N]
        products: P[C], Q[N]
    end
    @test all(x -> _testhelper_all_fields_equal(_testhelper_lift(ER.compile_mechanism(x)),
                                                x),
              EnzymeRates.Enumeration.init_mechanisms(rxn))
end

@testset "derived fields stay out of == and hash" begin
    # A mechanism's facts are functions of its canonical fields, so two mechanisms with
    # equal canonical fields are equal and hash alike whatever their facts hold: here a
    # rebuilt copy whose every fact is overwritten.
    m = ER.Mechanism(@enzyme_mechanism begin
        substrates: S
        products: P
        steps: begin
            E + S ⇌ E(S)
            E(S) <--> E(P)
            E(P) ⇌ E + P
        end
    end)
    x = ER.Mechanisms._unchecked_mechanism(m.substrates, m.products, m.nodes, m.steps)
    push!(x.modifiers, :S)
    x.groups[1] = [2, 3]
    x.step_kinds[1] = :chemistry
    x.toward_products .= .!x.toward_products
    x.step_sides[1] = ("E", "E")
    @test x == m && hash(x) == hash(m)
    @test length(unique!([m, x])) == 1
    am = ER.AllostericMechanism(m, [:EqualAI, :OnlyA, :NonequalAI], 2,
                                ER.Mechanisms.RegulatorySite[])
    y = ER.AllostericMechanism(m, [:EqualAI, :OnlyA, :NonequalAI], 2,
                               ER.Mechanisms.RegulatorySite[])
    for mask in (y.active_steps, y.inactive_steps, y.excluded_active, y.excluded_inactive)
        mask .= .!mask
    end
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
    @test fieldnames(ER.Mechanism)[1:4] == (:substrates, :products, :nodes, :steps)
    @test length(m.nodes) == 3
    @test m.nodes[1] == ER.Mechanisms.Node(:E, EnzymeRates.Mechanisms.Met[], Symbol[],
                                           Symbol[])
    # Nodes E, E(P), E(S); the release of P is stored as the binding it reverses, and
    # sorts first.
    @test [n.bound for n in m.nodes[2:3]] ==
          [[EnzymeRates.Mechanisms.Met(:P, false)], [EnzymeRates.Mechanisms.Met(:S, false)]]
    @test m.steps[1].taken_up == [EnzymeRates.Mechanisms.Met(:P, false)]
    @test m.steps[2].taken_up == [EnzymeRates.Mechanisms.Met(:S, false)]
    @test (m.steps[2].from, m.steps[2].to, m.steps[2].given_off) ==
          (1, 3, EnzymeRates.Mechanisms.Met[])
    @test m.steps[1].rapid_equilibrium && m.steps[2].rapid_equilibrium &&
          !m.steps[3].rapid_equilibrium
    @test [s.group for s in m.steps] == [1, 2, 3]
end

@testset "the constructor builds a spec's graph in under 20 µs, the median over specs" begin
    # Each spec's time per call is the best of five batches; the gate is on their median.
    # The spec's goal is about 2 µs and the median is about 3.5 µs locally; shared CI
    # runners take 3-4× as long, so the gate is 20 µs.
    per_call = map(MECHANISM_TEST_SPECS) do spec
        _testhelper_build(spec.raw)
        minimum(@elapsed(for _ in 1:200
                             _testhelper_build(spec.raw)
                         end) / 200 for _ in 1:5)
    end
    @test median(per_call) < 20e-6
end

"""The bytes one construction of `x` allocates, from its own canonical fields, after a
warm-up call: `Mechanism(substrates, products, nodes, steps)`, or for an allosteric
mechanism the constructor that also takes its tags, multiplicity and sites."""
function _testhelper_construction_bytes(x)
    cm = ER.Mechanisms.catalytic(x)
    _testhelper_construct() = x isa ER.Mechanism ?
        ER.Mechanism(cm.substrates, cm.products, cm.nodes, cm.steps) :
        ER.AllostericMechanism(cm.substrates, cm.products, cm.nodes, cm.steps, x.mwc_tags,
                               x.catalytic_multiplicity, x.regulatory_sites)
    _testhelper_construct()
    @allocated _testhelper_construct()
end

@testset "constructor allocation budget" begin
    # Each construction allocates at most 1.3 times the bytes it allocated at bf49d6e1 on
    # Julia 1.13 (6,672, 22,264 and 25,912 bytes): a mechanism stores the facts it computes
    # once, and they may not make every construction in a search meaningfully dearer.
    budget(bytes) = ceil(Int, 1.3 * bytes)
    michaelis_menten = ER.Mechanism(@enzyme_mechanism begin
        substrates: S
        products: P
        steps: begin
            E + S <--> E(S)
            E(S) <--> E + P
        end
    end)
    @test _testhelper_construction_bytes(michaelis_menten) <= budget(6_672)
    # The largest seed of the bi-bi init set (A, B → P, Q): 9 forms, 13 steps.
    random_bi_bi = ER.Mechanism(@enzyme_mechanism begin
        substrates: A, B
        products: P, Q
        steps: begin
            (E + A ⇌ E(A), E(B) + A ⇌ E(A, B), E(P) + A ⇌ E(A, P))
            (E + B ⇌ E(B), E(A) + B ⇌ E(A, B), E(Q) + B ⇌ E(B, Q))
            (E + P ⇌ E(P), E(A) + P ⇌ E(A, P), E(Q) + P ⇌ E(P, Q))
            (E + Q ⇌ E(Q), E(B) + Q ⇌ E(B, Q), E(P) + Q ⇌ E(P, Q))
            E(A, B) <--> E(P, Q)
        end
    end)
    @test _testhelper_construction_bytes(random_bi_bi) <= budget(22_264)
    # Hexokinase: 14 steps, two :OnlyA groups and a regulatory site.
    hk = _testhelper_lift(only(s for s in MECHANISM_TEST_SPECS if s.name == "HK").mechanism)
    @test _testhelper_construction_bytes(hk) <= budget(25_912)
end
