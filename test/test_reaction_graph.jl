# ABOUTME: Tests of the reaction-graph mechanism types: content order and identity,
# ABOUTME: canonical form and checks, its sub-graph reads, the Sig and speed.

using Statistics: median

const ER = EnzymeRates

"""The counts of `list`, table entries listed with repeats."""
_testhelper_pack(list) =
    foldl((c, i) -> c + (ER.Counts(1) << (4 * (i - 1))), list; init = ER.Counts(0))

"""Uni-uni with rapid-equilibrium bindings, written field by field in canonical form:
forms E, E(P), E(S); the bindings of P and S, then E(S) → E(P)."""
function _testhelper_graph_uniuni()
    P, S = _testhelper_pack([1]), _testhelper_pack([2])
    ER.Mechanism([:S], [:P], [EnzymeRates.Met(:P, false), EnzymeRates.Met(:S, false)],
                 [:E], [ER.Form(1, 0, 0, 0), ER.Form(1, P, 0, 0), ER.Form(1, S, 0, 0)],
                 [ER.Edge(1, 2, P, 0, true), ER.Edge(1, 3, S, 0, true),
                  ER.Edge(3, 2, 0, 0, false)],
                 [1, 2, 3])
end

@testset "packed counts list the table entries of a form or a step side" begin
    c = _testhelper_pack([1, 1, 3])
    @test (ER.copies(c, 1), ER.copies(c, 2), ER.copies(c, 3)) == (2, 0, 1)
    @test ER.total(c) == 3
    @test ER.entries(c) == [1, 1, 3]
    @test ER.entries(ER.Counts(0)) == Int[] && ER.total(ER.Counts(0)) == 0
    # The last of 32 entries, 15 copies: the top four bits.
    top = ER.Counts(15) << 124
    @test ER.copies(top, 32) == 15 && ER.total(top) == 15 && ER.entries(top) == fill(32, 15)
end

@testset "entry lists order lexicographically, a list before every list it begins" begin
    # Every list over four entries with at most two copies of each, against Julia's order
    # of the listed entries; the packed numbers order differently ([1, 2] < [2], 17 > 16).
    lists = [reduce(vcat, [fill(i, n[i]) for i in 1:4]; init = Int[])
             for n in Iterators.product(0:2, 0:2, 0:2, 0:2)]
    @test all(ER._seq_less(_testhelper_pack(a), _testhelper_pack(b)) == isless(a, b)
              for a in lists, b in lists)
    @test _testhelper_pack([1, 2]) > _testhelper_pack([2])
end

@testset "forms order by label, bound, added, subtracted; edges by ends, lists, flag" begin
    A, B, AB = _testhelper_pack([1]), _testhelper_pack([2]), _testhelper_pack([1, 2])
    @test ER.Form(1, B, 0, 0) < ER.Form(2, 0, 0, 0)
    @test ER.Form(1, AB, 0, 0) < ER.Form(1, B, 0, 0)
    @test ER.Form(1, B, 0, 0) < ER.Form(1, B, A, 0) < ER.Form(1, B, AB, 0)
    @test ER.Form(1, B, A, 0) < ER.Form(1, B, A, B)
    @test !(ER.Form(1, B, A, B) < ER.Form(1, B, A, B))
    @test ER.Edge(1, 3, 0, 0, false) < ER.Edge(2, 1, 0, 0, false)
    @test ER.Edge(1, 2, 0, 0, false) < ER.Edge(1, 3, 0, 0, false)
    @test ER.Edge(1, 2, AB, 0, true) < ER.Edge(1, 2, B, 0, false)
    @test ER.Edge(1, 2, A, AB, true) < ER.Edge(1, 2, A, B, false)
    @test ER.Edge(1, 2, A, 0, false) < ER.Edge(1, 2, A, 0, true)
    @test ER._reversed(ER.Edge(1, 2, A, B, true)) == ER.Edge(2, 1, B, A, true)
end

"""Whether a value of type `T` hashes the same in every process: a symbol, an integer or
a flag, a vector of such values, or one of the types whose hash folds its name and its
fields (`_testhelper_hash_folds`), whose fields are such values in turn."""
function _testhelper_hash_is_content(T)
    T in (Symbol, Bool, Int, UInt128) && return true
    T <: Vector && return _testhelper_hash_is_content(eltype(T))
    T in (EnzymeRates.Met, EnzymeRates.RegulatorySite, ER.Form, ER.Edge, ER.Mechanism,
          ER.AllostericMechanism) || return false
    all(_testhelper_hash_is_content, fieldtypes(T))
end

"""Whether `x` hashes its type's name, then its fields in declaration order."""
_testhelper_hash_folds(x, h) =
    hash(x, h) == foldl((acc, f) -> hash(getfield(x, f), acc), fieldnames(typeof(x));
                        init = hash(nameof(typeof(x)), h))

@testset "a graph hashes its content, the same in every process" begin
    # Julia hashes symbols, integers, flags and vectors of them alike in every process,
    # but a struct's default hash mixes in its type's identity, which changes between
    # precompile sessions. So each graph type hashes its name, then its fields, and its
    # fields hold only such values.
    m = _testhelper_graph_uniuni()
    site = EnzymeRates.RegulatorySite([:R], 2, [:OnlyI])
    am = ER.AllostericMechanism(m, [:EqualAI, :EqualAI, :OnlyA], 2, [site])
    h = UInt(17)
    for x in (m.forms[2], m.edges[3], m, am, m.table[1], site)
        @test _testhelper_hash_folds(x, h)
    end
    @test all(_testhelper_hash_is_content,
              (ER.Form, ER.Edge, ER.Mechanism, ER.AllostericMechanism))
    # Each field takes part in `==` and `hash`.
    @test m == _testhelper_graph_uniuni() && hash(m) == hash(_testhelper_graph_uniuni())
    for f in fieldnames(ER.Mechanism)
        fields = Any[getfield(m, g) for g in fieldnames(ER.Mechanism)]
        i = findfirst(==(f), fieldnames(ER.Mechanism))
        fields[i] = f === :group ? [1, 3, 2] : f === :substrates ? [:A] :
                    f === :products ? [:Q] : f === :conformations ? [:F] :
                    f === :table ? reverse(m.table) : reverse(getfield(m, f))
        other = ER.Mechanism(fields...)
        @test other != m && hash(other) != hash(m)
    end
    @test am == ER.AllostericMechanism(m, [:EqualAI, :EqualAI, :OnlyA], 2, [site])
    @test am != ER.AllostericMechanism(m, [:EqualAI, :OnlyA, :OnlyA], 2, [site])
    @test am != ER.AllostericMechanism(m, [:EqualAI, :EqualAI, :OnlyA], 4, [site])
    @test am != ER.AllostericMechanism(m, [:EqualAI, :EqualAI, :OnlyA], 2,
                                       EnzymeRates.RegulatorySite[])
end

@testset "a mechanism's fields are its graph and nothing else" begin
    # No reaction and no cache: the names, the table, the forms, the edges and their
    # groups, and for an allosteric mechanism its tags, multiplicity and sites.
    @test fieldnames(ER.Mechanism) ==
          (:substrates, :products, :table, :conformations, :forms, :edges, :group)
    @test fieldnames(ER.Form) == (:conformation, :bound, :added, :subtracted)
    @test fieldnames(ER.Edge) == (:from, :to, :consumed, :released, :rapid_equilibrium)
    @test fieldnames(ER.AllostericMechanism) ==
          (:catalytic, :cat_allo_states, :catalytic_multiplicity, :regulatory_sites)
end

@testset "an allosteric graph sorts its sites and checks its tags and ligands" begin
    m = _testhelper_graph_uniuni()
    tags = [:EqualAI, :EqualAI, :OnlyA]
    sites = [EnzymeRates.RegulatorySite([:X], 2, [:OnlyI]),
             EnzymeRates.RegulatorySite([:B, :A], 2, [:OnlyA, :EqualAI])]
    am = ER.AllostericMechanism(m, tags, 2, sites)
    @test [s.ligands for s in am.regulatory_sites] == [[:A, :B], [:X]]
    @test am.regulatory_sites[1].allo_states == [:EqualAI, :OnlyA]
    tags[1] = :OnlyA
    @test am.cat_allo_states == [:EqualAI, :EqualAI, :OnlyA]
    @test [ER.cat_allo_state(am, g) for g in 1:3] == [:EqualAI, :EqualAI, :OnlyA]
    @test_throws(ErrorException("AllostericMechanism: cat_allo_states length 2 must " *
                                "match the number of kinetic groups 3"),
                 ER.AllostericMechanism(m, [:EqualAI, :EqualAI], 2, sites))
    @test_throws(ErrorException("AllostericMechanism: catalytic_multiplicity must be ≥ " *
                                "1, got 0"),
                 ER.AllostericMechanism(m, [:EqualAI, :EqualAI, :OnlyA], 0, sites))
    @test_throws(ErrorException(
        "AllostericMechanism: catalytic group 2 has invalid allo state OnlyI (must be " *
        "one of (:OnlyA, :EqualAI, :NonequalAI)); :OnlyI is rejected for catalytic " *
        "groups (active-state-active convention)"),
        ER.AllostericMechanism(m, [:EqualAI, :OnlyI, :OnlyA], 2, sites))
    @test_throws(ErrorException(
        "AllostericMechanism: ligand X appears in two distinct regulatory sites; " *
        "rendered Kreg names would collide. Same-ligand-two-sites is not enumerated."),
        ER.AllostericMechanism(m, [:EqualAI, :EqualAI, :OnlyA], 2,
                               [sites; EnzymeRates.RegulatorySite([:X], 4, [:OnlyA])]))
    @test_throws(ErrorException(
        "AllostericMechanism: L is a reserved name; E_total, Keq, L, N, params, concs " *
        "and v are kept for a rate equation's constants and local names"),
        ER.AllostericMechanism(m, [:EqualAI, :EqualAI, :OnlyA], 2,
                               [EnzymeRates.RegulatorySite([:L], 2, [:OnlyI])]))
end

"""`raw` written another way, the same mechanism: with `order`, its steps in reverse
order; with `direction`, each from its other end; with `labels`, its group labels
numbered from the last. Its declared names come in reverse order."""
function _testhelper_rewrite(raw::EnzymeRates.RawMechanism;
                             order = true, direction = true, labels = true)
    steps, group = raw.steps, raw.group
    direction && (steps = [EnzymeRates.RawStep(s.to, s.from, s.released, s.consumed,
                                               s.rapid_equilibrium) for s in steps])
    labels && (group = maximum(group; init = 0) + 1 .- group)
    order && ((steps, group) = (reverse(steps), reverse(group)))
    EnzymeRates.RawMechanism(reverse(raw.substrates), reverse(raw.products),
                             reverse(raw.modifiers), steps, group)
end

"""The catalytic part of `raw`, a spec's raw value."""
_testhelper_plain(raw::EnzymeRates.RawMechanism) = raw
_testhelper_plain(raw::EnzymeRates.RawAllostericMechanism) = raw.catalytic

@testset "every writing of a mechanism gives one canonical graph" begin
    @testset "$(spec.name)" for spec in MECHANISM_TEST_SPECS
        raw = _testhelper_plain(spec.raw)
        m = ER._canonicalize(raw).mechanism
        for (order, direction, labels) in ((true, false, false), (false, true, false),
                                           (false, false, true), (true, true, true))
            other = ER._canonicalize(
                _testhelper_rewrite(raw; order, direction, labels)).mechanism
            @test other == m && hash(other) == hash(m)
        end
    end
end

@testset "canonical form maps each written step to its edge, group and forms" begin
    c = ER._canonicalize(EnzymeRates._raw_mechanism(:(begin
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
    @test c.form_of == [1, 3, 3, 2, 2, 1]
    @test isempty(c.misaligned)
end

@testset "canonical form runs each step physical-forward, tied steps as one" begin
    _testhelper_flipped(block) =
        ER._canonicalize(EnzymeRates._raw_mechanism(block, false)).flipped
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
    _testhelper_mirrored(iso) = EnzymeRates._raw_mechanism(:(begin
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
    c, d = ER._canonicalize(_testhelper_mirrored(:(F(I) <--> E(I)))),
           ER._canonicalize(_testhelper_mirrored(:(E(I) <--> F(I))))
    @test c.flipped == [false, false, true, false, false, false, false]
    @test d.flipped == [false, false, true, false, true, false, false]
    @test c.mechanism == d.mechanism
    # Tied steps that change different pairs of conformations: no orientation aligns
    # them, and canonical form returns the lead and the first that misses it.
    c = ER._canonicalize(EnzymeRates._raw_mechanism(:(begin
        substrates: A
        products: P
        steps: begin
            E + A <--> E(A)
            E(A) <--> E(P)
            E(P) <--> E + P
            (E ⇌ F, F ⇌ G, G ⇌ E)
        end
    end), false))
    @test c.misaligned == [ER.Edge(1, 4, 0, 0, true), ER.Edge(4, 5, 0, 0, true)]
end

"""Bindings of each name of `names` to free E, at rapid equilibrium."""
_testhelper_bindings(names) =
    [EnzymeRates.RawStep(EnzymeRates.RawForm(:E, EnzymeRates.Met[]),
                         EnzymeRates.RawForm(:E, [EnzymeRates.Met(x, false)]),
                         [EnzymeRates.Met(x, false)], EnzymeRates.Met[], true)
     for x in names]

@testset "a graph holds at most 32 table entries and 15 copies of one entry" begin
    names = [Symbol("S", i) for i in 1:32]
    _testhelper_raw(n) = EnzymeRates.RawMechanism(names[1:(n - 1)], [:P], Symbol[],
                                      _testhelper_bindings([names[1:(n - 1)]; :P]),
                                      collect(1:n))
    @test length(ER._canonicalize(_testhelper_raw(32)).mechanism.table) == 32
    @test_throws(ErrorException("Mechanism: a mechanism holds at most 32 metabolite " *
                                "entries (a name and its modifier copy count apart); " *
                                "this one has 33"),
                 ER._canonicalize(_testhelper_raw(33)))
    # A reactant and its modifier copy are two entries.
    _testhelper_copy_of(x) = EnzymeRates.RawStep(EnzymeRates.RawForm(:E, EnzymeRates.Met[]),
                                     EnzymeRates.RawForm(:E, [EnzymeRates.Met(x, true)]),
                                     [EnzymeRates.Met(x, true)], EnzymeRates.Met[], true)
    both = EnzymeRates.RawMechanism(names[1:16], [:P], names[1:16],
                                    [_testhelper_bindings([names[1:16]; :P]);
                                     _testhelper_copy_of.(names[1:16])], collect(1:33))
    @test_throws(ErrorException("Mechanism: a mechanism holds at most 32 metabolite " *
                                "entries (a name and its modifier copy count apart); " *
                                "this one has 33"),
                 ER._canonicalize(both))
    _testhelper_holding(n) = EnzymeRates.RawMechanism([:A], [:P], Symbol[],
        [EnzymeRates.RawStep(EnzymeRates.RawForm(:E, fill(EnzymeRates.Met(:A, false), n)),
                             EnzymeRates.RawForm(:F, fill(EnzymeRates.Met(:A, false), n)),
                             EnzymeRates.Met[], EnzymeRates.Met[], false)], [1])
    @test ER.copies(
        ER._canonicalize(_testhelper_holding(15)).mechanism.forms[1].bound, 1) == 15
    @test_throws(ErrorException("Mechanism: a form or a step side holds A more than 15 " *
                                "times; a mechanism holds at most 15 copies of one " *
                                "metabolite there"),
                 ER._canonicalize(_testhelper_holding(16)))
end

"""An allosteric `raw` written another way, the same mechanism: its catalytic steps
rewritten (`_testhelper_rewrite`), its tags following their labels, and its sites and
each site's ligands in reverse order."""
_testhelper_rewrite(raw::EnzymeRates.RawAllostericMechanism) =
    EnzymeRates.RawAllostericMechanism(
        _testhelper_rewrite(raw.catalytic), reverse(raw.cat_allo_states),
        raw.catalytic_multiplicity,
        [EnzymeRates.RawSite(reverse(s.ligands), s.multiplicity, reverse(s.allo_states))
         for s in reverse(raw.sites)])

"""The graph a spec's raw value `raw` writes, built and checked."""
_testhelper_build(raw::EnzymeRates.RawMechanism) = ER.Mechanism(raw)
_testhelper_build(raw::EnzymeRates.RawAllostericMechanism) = ER.AllostericMechanism(raw)

"""The raw value of `block`, a mechanism macro's block, read by the DSL walker."""
_testhelper_walk(block; allosteric = false) =
    EnzymeRates._raw_mechanism(block, allosteric)

@testset "the constructors build one value from every writing of a mechanism" begin
    @testset "$(spec.name)" for spec in MECHANISM_TEST_SPECS
        x = _testhelper_build(spec.raw)
        @test ER.catalytic(x) == ER._canonicalize(_testhelper_plain(spec.raw)).mechanism
        other = _testhelper_build(_testhelper_rewrite(spec.raw))
        @test other == x && hash(other) == hash(x)
    end
    # Each tag goes to the canonical group its label became.
    am = ER.AllostericMechanism(_testhelper_walk(:(begin
        substrates: S
        products: P
        catalytic_multiplicity: 2
        steps: begin
            E + S ⇌ E(S) :: EqualAI
            E(S) <--> E(P) :: OnlyA
            E(P) ⇌ E + P :: NonequalAI
        end
    end); allosteric = true))
    @test am.cat_allo_states == [:NonequalAI, :EqualAI, :OnlyA]
end

"""Whether `m` stores each of its vectors at its length: `m` takes the memory of the
value rebuilt from a copy of each vector, which `copy` stores at its length."""
_testhelper_stored_at_length(m::ER.Mechanism) = Base.summarysize(m) == Base.summarysize(
    ER.Mechanism((copy(getfield(m, f)) for f in fieldnames(ER.Mechanism))...))
_testhelper_stored_at_length(am::ER.AllostericMechanism) =
    _testhelper_stored_at_length(am.catalytic) && Base.summarysize(am) ==
    Base.summarysize(ER.AllostericMechanism(am.catalytic, copy(am.cat_allo_states),
        am.catalytic_multiplicity,
        [ER.RegulatorySite(copy(s.ligands), s.multiplicity, copy(s.allo_states))
         for s in am.regulatory_sites]))

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
    init = EnzymeRates.init_mechanisms(rxn)
    children = EnzymeRates.expand_mechanisms(init[1:10], rxn)
    @test any(c -> c isa ER.AllostericMechanism, children)
    @test all(_testhelper_stored_at_length, init)
    @test all(_testhelper_stored_at_length, children)
end

@testset "the constructors reject what no graph represents, each with its message" begin
    _testhelper_uniuni_raw(; substrates = [:S], products = [:P], modifiers = Symbol[],
                           group = [1, 2, 3]) =
        EnzymeRates.RawMechanism(substrates, products, modifiers,
                                 _testhelper_walk(:(begin
                                     substrates: S
                                     products: P
                                     steps: begin
                                         E + S ⇌ E(S)
                                         E(S) <--> E(P)
                                         E(P) ⇌ E + P
                                     end
                                 end)).steps, group)
    _testhelper_allosteric(tags, n, sites) =
        ER.AllostericMechanism(EnzymeRates.RawAllostericMechanism(
            _testhelper_uniuni_raw(), tags, n, sites))
    tags = [:EqualAI, :OnlyA, :EqualAI]
    site(ligands, n, states) = EnzymeRates.RawSite(ligands, n, states)
    # Regulatory sites.
    @test_throws(ErrorException("RegulatorySite: needs at least one ligand"),
                 _testhelper_allosteric(tags, 2, [site(Symbol[], 2, Symbol[])]))
    @test_throws(ErrorException("RegulatorySite: length(ligands)=1 must equal " *
                                "length(allo_states)=0"),
                 _testhelper_allosteric(tags, 2, [site([:R], 2, Symbol[])]))
    @test_throws(ErrorException("RegulatorySite: multiplicity must be ≥ 1, got 0"),
                 _testhelper_allosteric(tags, 2, [site([:R], 0, [:OnlyA])]))
    @test_throws(ErrorException("RegulatorySite: allo state Bogus must be one of " *
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
    @test_throws(ErrorException("AllostericMechanism: cat_allo_states length 2 must " *
                                "match the number of kinetic groups 3"),
                 _testhelper_allosteric(tags[1:2], 2, EnzymeRates.RawSite[]))
    @test_throws(ErrorException("AllostericMechanism: catalytic_multiplicity must be ≥ " *
                                "1, got 0"),
                 _testhelper_allosteric(tags, 0, EnzymeRates.RawSite[]))
    # A Bool is no multiplicity: each constructor takes an Int and converts no Bool.
    @test_throws MethodError site([:R], true, [:OnlyA])
    @test_throws MethodError EnzymeRates.RawAllostericMechanism(
        _testhelper_uniuni_raw(), tags, true, EnzymeRates.RawSite[])
    @test_throws MethodError EnzymeRates.RegulatorySite([:R], true, [:OnlyA])
    @test_throws MethodError ER.AllostericMechanism(
        ER.Mechanism(_testhelper_uniuni_raw()), tags, true, EnzymeRates.RegulatorySite[])
    @test_throws(ErrorException(
        "AllostericMechanism: catalytic group 2 has invalid allo state OnlyI (must be " *
        "one of (:OnlyA, :EqualAI, :NonequalAI)); :OnlyI is rejected for catalytic " *
        "groups (active-state-active convention)"),
        _testhelper_allosteric([:EqualAI, :OnlyI, :EqualAI], 2, EnzymeRates.RawSite[]))
    # Declarations and group labels.
    @test_throws(ErrorException("Mechanism: substrates must not be empty"),
                 ER.Mechanism(_testhelper_uniuni_raw(substrates = Symbol[])))
    @test_throws(ErrorException("Mechanism: products must not be empty"),
                 ER.Mechanism(_testhelper_uniuni_raw(products = Symbol[])))
    @test_throws(ErrorException("Mechanism: duplicate substrate names"),
                 ER.Mechanism(_testhelper_uniuni_raw(substrates = [:S, :S])))
    @test_throws(ErrorException("Mechanism: duplicate product names"),
                 ER.Mechanism(_testhelper_uniuni_raw(products = [:P, :P])))
    @test_throws(ErrorException("Mechanism: duplicate modifier names"),
                 ER.Mechanism(_testhelper_uniuni_raw(modifiers = [:I, :I])))
    @test_throws(ErrorException("Mechanism: S named as both a substrate and a product; " *
                                "concentrations and constants are keyed by name"),
                 ER.Mechanism(_testhelper_uniuni_raw(products = [:P, :S])))
    for group in ([1, 2], [1, 2, 3, 3], [0, 1, 2])
        @test_throws(ErrorException("Mechanism: each step takes one kinetic group label " *
                                    "from 1 to the number of groups"),
                     ER.Mechanism(_testhelper_uniuni_raw(; group)))
    end
    @test_throws(ErrorException("Mechanism: kinetic group 2 is empty; a kinetic group " *
                                "holds at least one step"),
                 ER.Mechanism(_testhelper_uniuni_raw(group = [1, 3, 3])))
    # Declared names a metabolite cannot take: no identifier, a keyword, E, a reserved
    # name, a name holding `_`.
    _testhelper_unwritable(x) = ErrorException(
        "Mechanism: \"$x\" cannot name a metabolite; a metabolite name is an identifier " *
        "that is no keyword")
    @test_throws(_testhelper_unwritable("1I"),
                 ER.Mechanism(_testhelper_uniuni_raw(modifiers = [Symbol("1I")])))
    @test_throws(_testhelper_unwritable("end"),
                 ER.Mechanism(_testhelper_uniuni_raw(products = [:P, :end])))
    @test_throws(ErrorException("Mechanism: E cannot name a metabolite; it labels the " *
                                "free enzyme form"),
                 ER.Mechanism(_testhelper_uniuni_raw(substrates = [:S, :E])))
    @test_throws(ErrorException("Mechanism: Keq is a reserved name; E_total, Keq, L, N, " *
                                "params, concs and v are kept for a rate equation's " *
                                "constants and local names"),
                 ER.Mechanism(_testhelper_uniuni_raw(substrates = [:S, :Keq])))
    @test_throws(ErrorException("Mechanism: the metabolite name S_1 holds `_`; constant " *
                                "names join metabolite names with `_`, so a metabolite " *
                                "name holds none"),
                 ER.Mechanism(_testhelper_uniuni_raw(modifiers = [:S_1])))
    # A name Julia reads as another (a micro sign as μ, a decomposed é as one character)
    # prints as text that builds another mechanism, so no DSL text can write it. The
    # strings are written with escapes: Julia would read the characters themselves as
    # the other names.
    _testhelper_read_as(who, x, read) = ErrorException(
        "$who: Julia reads the name \"$x\" as \"$read\", so no DSL text can write it; " *
        "use \"$read\"")
    @test_throws(_testhelper_read_as("Mechanism", "\u00b5M", "\u03bcM"),
                 ER.Mechanism(_testhelper_uniuni_raw(modifiers = [Symbol("\u00b5M")])))
    @test_throws(_testhelper_read_as("Mechanism", "e\u0301", "\u00e9"),
                 ER.Mechanism(_testhelper_uniuni_raw(substrates = [:S, Symbol("e\u0301")])))
    @test_throws(_testhelper_read_as("AllostericMechanism", "\u025b", "\u03b5"),
                 _testhelper_allosteric(tags, 2, [site([Symbol("\u025b")], 2, [:OnlyA])]))
    # Names: undeclared; a modifier entry of a name modifiers: does not list; a plain
    # entry of a modifier-only name; residuals by role.
    @test_throws(ErrorException("Mechanism: P is not declared; declare it under " *
                                "substrates:, products: or modifiers:"),
                 ER.Mechanism(_testhelper_uniuni_raw(products = [:Q])))
    mod_P = EnzymeRates.Met(:P, true)
    free = EnzymeRates.RawForm(:E, EnzymeRates.Met[])
    copy_binding = EnzymeRates.RawStep(free, EnzymeRates.RawForm(:E, [mod_P]), [mod_P],
                                       EnzymeRates.Met[], true)
    _testhelper_adding(raw, steps...) = EnzymeRates.RawMechanism(
        raw.substrates, raw.products, raw.modifiers, [raw.steps; steps...],
        [raw.group; length(raw.steps) .+ (1:length(steps))])
    @test_throws(ErrorException("Mechanism: P binds as a modifier, but modifiers: does " *
                                "not declare it"),
                 ER.Mechanism(_testhelper_adding(_testhelper_uniuni_raw(), copy_binding)))
    @test ER.Mechanism(_testhelper_adding(_testhelper_uniuni_raw(modifiers = [:P]),
                                          copy_binding)) isa ER.Mechanism
    I = EnzymeRates.Met(:I, false)
    @test_throws(ErrorException("Mechanism: I is a modifier, so it binds only as a " *
                                "modifier"),
                 ER.Mechanism(_testhelper_adding(
                     _testhelper_uniuni_raw(modifiers = [:I]),
                     EnzymeRates.RawStep(free, EnzymeRates.RawForm(:E, [I]), [I],
                                         EnzymeRates.Met[], true))))
    S = EnzymeRates.Met(:S, false)
    for (added, subtracted, verb, x) in (([:P], Symbol[], "adds", :P),
                                         (Symbol[], [:S], "subtracts", :S))
        residual = EnzymeRates.RawForm(:E, EnzymeRates.Met[], added, subtracted)
        @test_throws(ErrorException("Mechanism: a residual adds substrates and subtracts " *
                                    "products; one $verb $x"),
                     ER.Mechanism(_testhelper_adding(_testhelper_uniuni_raw(),
                         EnzymeRates.RawStep(EnzymeRates.RawForm(:E, [S]), residual,
                                             EnzymeRates.Met[], EnzymeRates.Met[],
                                             false))))
    end
    # Steps: a self-loop, a metabolite taken up and given off, a role change.
    EA = EnzymeRates.RawForm(:E, [EnzymeRates.Met(:A, false)])
    @test_throws(ErrorException("Step: both ends are the form EA; a step joins two forms"),
                 ER.Mechanism(EnzymeRates.RawMechanism([:A], [:P], Symbol[],
                     [EnzymeRates.RawStep(EA, EA, EnzymeRates.Met[], EnzymeRates.Met[],
                                          false)], [1])))
    @test_throws(ErrorException("Step EA → EQ: B both consumed and released"),
                 ER.Mechanism(_testhelper_walk(:(begin
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
                 ER.Mechanism(_testhelper_walk(:(begin
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
        @test_throws _testhelper_role_change(from, to) ER.Mechanism(_testhelper_walk(:(begin
            substrates: S
            products: P
            modifiers: P
            steps: begin
                E + S ⇌ E(S)
                E(S) <--> E(P)
                E(P) ⇌ E + P
                $line
            end
        end)))
    end
    # Tied steps that make different conformational changes.
    @test_throws(ErrorException(
        "Mechanism: a kinetic group holds E → F and F → G, which make different " *
        "conformational changes; the steps of a kinetic group share their constants, so " *
        "they must make the same conformational change"),
        ER.Mechanism(_testhelper_walk(:(begin
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
    @test_throws _testhelper_one_name("E(NAD, P)", "E(NADP)", "ENADP") ER.Mechanism(
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
    @test_throws _testhelper_one_name("E(star)", "Estar", "Estar") ER.Mechanism(
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
    @test_throws _testhelper_one_name("E(P::Mod)", "E(Pmod)", "EPmod") ER.Mechanism(
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
            ER.Mechanism(_testhelper_walk(:(begin
                substrates: A, B
                products: P, Q
                steps: $steps
            end))))
    end
    # One reaction twice: in two groups, with different flags, stored once each way, and in
    # one group.
    _testhelper_twice(steps) = ER.Mechanism(_testhelper_walk(:(begin
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
        ER.Mechanism(_testhelper_walk(:(begin
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
    @test_throws no_free_E ER.Mechanism(EnzymeRates.RawMechanism(
        [:S], [:P], Symbol[], EnzymeRates.RawStep[], Int[]))
    @test_throws no_free_E ER.Mechanism(_testhelper_walk(:(begin
        substrates: S
        products: P
        steps: begin
            F + S ⇌ F(S)
            F(S) <--> F(P)
            F(P) ⇌ F + P
        end
    end)))
    @test_throws no_free_E ER.Mechanism(_testhelper_walk(:(begin
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
    @test_throws no_free_E ER.Mechanism(_testhelper_walk(:(begin
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
                 ER.Mechanism(_testhelper_walk(:(begin
        substrates: S
        products: P
        steps: begin
            E + S ⇌ E(S)
            E(S) <--> E(P)
            E(P) ⇌ E + P
            Estar(S) <--> Estar(P)
        end
    end))))
    # Declared names no step exchanges: a substrate, a product, a modifier; a substrate
    # exchanged only as its modifier copy; a reactant declared a modifier whose copy no
    # step binds.
    _testhelper_unused(kind, x) = ErrorException(
        "Mechanism: $kind $x is declared, but no step " *
        (kind == "modifier" ? "binds" : "uses") * " it")
    @test_throws(_testhelper_unused("substrate", :T),
                 ER.Mechanism(_testhelper_uniuni_raw(substrates = [:S, :T])))
    @test_throws(_testhelper_unused("product", :Q),
                 ER.Mechanism(_testhelper_uniuni_raw(products = [:P, :Q])))
    @test_throws(_testhelper_unused("modifier", :I),
                 ER.Mechanism(_testhelper_uniuni_raw(modifiers = [:I])))
    mod_T = EnzymeRates.Met(:T, true)
    @test_throws(_testhelper_unused("substrate", :T),
                 ER.Mechanism(_testhelper_adding(
                     _testhelper_uniuni_raw(substrates = [:S, :T], modifiers = [:T]),
                     EnzymeRates.RawStep(free, EnzymeRates.RawForm(:E, [mod_T]), [mod_T],
                                         EnzymeRates.Met[], true))))
    @test_throws(_testhelper_unused("modifier", :P),
                 ER.Mechanism(_testhelper_uniuni_raw(modifiers = [:P])))
    # Conformation labels: a label shaped otherwise than E, Estar or F2; an opaque
    # bound-form name; A or I, which constant names use as conformation tags; a
    # metabolite's name, of the catalytic graph or of a ligand.
    _testhelper_turning_to(label) = _testhelper_adding(_testhelper_uniuni_raw(),
        EnzymeRates.RawStep(free, EnzymeRates.RawForm(label, EnzymeRates.Met[]),
                            EnzymeRates.Met[], EnzymeRates.Met[], false))
    @test_throws(ErrorException("Mechanism: `Ex_1` is not a conformation label: a " *
                                "conformation label is a capital letter followed by " *
                                "lowercase letters and digits, such as `E`, `Estar` or " *
                                "`F2`."),
                 ER.Mechanism(_testhelper_turning_to(:Ex_1)))
    @test_throws(ErrorException("Mechanism: `ES` looks like an opaque bound-form name; " *
                                "write it as decomposed call notation, e.g. `E(S)` or " *
                                "`E(A, B)`."),
                 ER.Mechanism(_testhelper_turning_to(:ES)))
    for label in (:A, :I)
        @test_throws(ErrorException("Mechanism: `$label` cannot be a conformation label: " *
                                    "constant names tag a conformation with A_ or I_ " *
                                    "(`k_A_ES_to_EP`), so the labels A and I could give " *
                                    "two constants one name; choose another label."),
                     ER.Mechanism(_testhelper_turning_to(label)))
    end
    @test_throws(ErrorException("Mechanism: conformation label `S` collides with " *
                                "declared metabolite `S`; choose a different " *
                                "conformation label."),
                 ER.Mechanism(_testhelper_turning_to(:S)))
    @test_throws(ErrorException("AllostericMechanism: conformation label `F` collides " *
                                "with ligand `F`; choose a different conformation " *
                                "label."),
                 ER.AllostericMechanism(EnzymeRates.RawAllostericMechanism(
                     _testhelper_turning_to(:F), [tags; :EqualAI], 2,
                     [site([:F], 2, [:OnlyA])])))
end

@testset "a form's name and text read as constants and the macros write them" begin
    m = ER.Mechanism(_testhelper_walk(:(begin
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
    @test [ER._form_name(m, f) for f in eachindex(m.forms)] ==
          ["E", "E_res_plusA_minusP", "EA", "EB_res_plusA_minusP", "EI", "EQ", "EQQmod"]
    @test [ER._form_text(m, f) for f in eachindex(m.forms)] ==
          ["E", "E(; residual = A - P)", "E(A)", "E(B; residual = A - P)", "E(I)",
           "E(Q)", "E(Q, Q::Mod)"]
    # Each edge's sides in its stored direction: releases stored as bindings.
    @test [join(ER._sides(m, k), " → ") for k in eachindex(m.edges)] ==
          ["E_A → EA", "E_I → EI", "E_Q → EQ", "E_res_plusA_minusP_P → EA",
           "E_res_plusA_minusP_B → EB_res_plusA_minusP", "EB_res_plusA_minusP → EQ",
           "EQ_Qmod → EQQmod"]
end

@testset "a sub-graph reads the forms, segments and bottom of the edges a mask keeps" begin
    # Random release of P and Q at rapid equilibrium between steady-state releases: the
    # segment {E(P), E(P, Q), E(Q)} has every form holding a product.
    m = ER.Mechanism(_testhelper_walk(:(begin
        substrates: S
        products: P, Q
        steps: begin
            E + S ⇌ E(S)
            E(S) <--> E(P, Q)
            E(P, Q) ⇌ E(Q) + P
            E(P, Q) ⇌ E(P) + Q
            E(P) <--> E + P
            E(Q) <--> E + Q
        end
    end)))
    names(forms) = [ER._form_name(m, f) for f in forms]
    @test names(eachindex(m.forms)) == ["E", "EP", "EPQ", "EQ", "ES"]
    # One edge per group, so the edges group by group are the edges in order.
    @test [_testhelper_edge_text(m, k) for k in ER._flat_edges(m)] ==
          ["E + P <--> EP", "E + Q <--> EQ", "E + S ⇌ ES", "EP + Q ⇌ EPQ", "EQ + P ⇌ EPQ",
           "ES <--> EPQ"]
    @test ER._flags(m) == BitVector([0, 0, 1, 1, 1, 0])
    @test [ER._is_free(f) for f in m.forms] == [true, false, false, false, false]
    every = trues(length(m.edges))
    @test ER._present_forms(m, every) == trues(5)
    # Each segment breadth-first from its first form; offsets per metabolite name (rows
    # P, Q, S) above the segment's lowest.
    s = ER._segments(m, every, ER._flags(m))
    @test [names(seg) for seg in s.segments] == [["E", "ES"], ["EP", "EPQ", "EQ"]]
    @test s.segment_of == [1, 2, 2, 2, 1]
    @test s.offsets == [0 1 1 0 0; 0 0 1 1 0; 0 0 0 0 1]
    @test names(ER._bottomless_segment(m, every, ER._flags(m))) == ["EP", "EPQ", "EQ"]
    # Without the release of P from E(P, Q), E(Q) stands alone and E(P) is the bottom of
    # {E(P), E(P, Q)}.
    no_p = BitVector([1, 1, 1, 1, 0, 1])
    @test ER._present_forms(m, no_p) == trues(5)
    s = ER._segments(m, no_p, ER._flags(m))
    @test [names(seg) for seg in s.segments] == [["E", "ES"], ["EP", "EPQ"], ["EQ"]]
    @test isempty(ER._bottomless_segment(m, no_p, ER._flags(m)))
    # A mask keeps only the forms its edges touch; a form it leaves is in no segment.
    bindings = BitVector([1, 1, 0, 0, 0, 0])
    @test ER._present_forms(m, bindings) == BitVector([1, 1, 0, 1, 0])
    s = ER._segments(m, bindings, ER._flags(m))
    @test [names(seg) for seg in s.segments] == [["E"], ["EP"], ["EQ"]]
    @test s.segment_of == [1, 2, 0, 3, 0]
    # `re` picks the edges that join segments: with none, each present form stands alone.
    @test length(ER._segments(m, every, falses(6)).segments) == 5

    # The inactive conformation keeps the edges outside :OnlyA groups between forms that
    # such edges reach from a free form: here E(A) and E(A, I) lie behind the :OnlyA
    # binding of A.
    am = ER.AllostericMechanism(_testhelper_walk(:(begin
        substrates: A
        products: P
        modifiers: I
        catalytic_multiplicity: 1
        steps: begin
            E + A <--> E(A)                                   :: OnlyA
            E(A) <--> E(P)                                    :: OnlyA
            E + P ⇌ E(P)                                      :: EqualAI
            (E(A) + I ⇌ E(A, I), E(P) + I ⇌ E(I, P))          :: NonequalAI
        end
    end); allosteric = true))
    inactive = ER._inactive_mask(am)
    @test [_testhelper_edge_text(am, k) for k in findall(inactive)] ==
          ["E + P ⇌ EP", "EP + I ⇌ EIP"]
    @test [ER._form_name(am.catalytic, f)
           for f in findall(ER._present_forms(am.catalytic, inactive))] ==
          ["E", "EIP", "EP"]
end

@testset "base ranks order forms and edges by content, free E among the forms" begin
    # Free E is form 1, but by content it sorts after the free form D.
    m = ER.Mechanism(_testhelper_walk(:(begin
        substrates: A
        products: P
        steps: begin
            E + A <--> E(A)
            E(A) <--> E(P)
            E(P) <--> D + P
            D <--> E
        end
    end)))
    @test [ER._form_name(m, f) for f in eachindex(m.forms)] == ["E", "D", "EA", "EP"]
    @test ER._base_rank(m) == [2, 1, 3, 4]
    # P enters at D and A at E, so the isomerization runs from D (orientation tier 2).
    texts = [_testhelper_edge_text(m, k) for k in eachindex(m.edges)]
    @test texts == ["E + A <--> EA", "D <--> E", "D + P <--> EP", "EA <--> EP"]
    # By base rank the edges from D come first.
    @test ER._edge_base_rank(m) == [3, 1, 2, 4]
    @test [ER._is_iso(e) for e in m.edges] == [false, true, false, true]
    @test [ER._binds_product(m, e) for e in m.edges] == [false, false, true, false]
end

"""The graph whose Sig is `sig`, `x`'s kind of graph. `sig` is splatted with its length
left open, as the lifts splat it, so the decoders compile once, whatever the Sig length."""
_testhelper_decode(::ER.Mechanism, @nospecialize(sig)) = ER._mechanism_from_sig(sig...)
_testhelper_decode(::ER.AllostericMechanism, @nospecialize(sig)) =
    ER._allosteric_mechanism_from_sig(sig...)

@testset "a graph's Sig decodes to the graph, and equal graphs share one type" begin
    @testset "$(spec.name)" for spec in MECHANISM_TEST_SPECS
        x = _testhelper_build(spec.raw)
        sig = ER._sig_of(x)
        @test sig isa Tuple{Vararg{Int}}
        @test _testhelper_decode(x, sig) == x
        cx = ER.compile_mechanism(x)
        @test typeof(cx).parameters[1] === sig
        @test ER._concrete(cx) == x
        @test ER.compile_mechanism(_testhelper_build(_testhelper_rewrite(spec.raw))) === cx
    end
    m = _testhelper_graph_uniuni()
    @test ER.Mechanism(ER.EnzymeMechanism(m)) == m
    am = ER.AllostericMechanism(m, [:EqualAI, :NonequalAI, :OnlyA], 4,
                                [EnzymeRates.RegulatorySite([:X, :R], 2, [:OnlyA, :OnlyI])])
    aem = ER.AllostericEnzymeMechanism(am)
    @test ER.AllostericMechanism(aem) == am
    @test ER.catalytic_mechanism(aem) === ER.EnzymeMechanism(m)
    @test ER.catalytic_multiplicity(aem) == 4
    # A name outside ASCII, and table entries past the sixteenth, whose counts fill the
    # high half of a packed count: αKG sorts last, entry 32, the top four bits.
    names = [[Symbol("S", i) for i in 1:30]; :αKG]
    wide = ER.Mechanism(EnzymeRates.RawMechanism(
        names, [:P], Symbol[], _testhelper_bindings([names; :P]), collect(1:32)))
    @test wide.table[32] == EnzymeRates.Met(:αKG, false)
    @test _testhelper_decode(wide, ER._sig_of(wide)) == wide
end

"""Uni-uni whose substrate's name holds `n` characters, so its Sig is longer than any
other test's: a Sig length no other test decodes."""
_testhelper_long_named(n) = let S = Symbol("S", "x"^(n - 1))
    ER.Mechanism(_testhelper_walk(:(begin
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
                               EnzymeRates.RegulatorySite[])
    for (decoder, build) in ((ER._mechanism_from_sig, _testhelper_long_named),
                             (ER._allosteric_mechanism_from_sig,
                              _testhelper_allosteric_long_named))
        lengths = (301, 401, 501, 601)
        compiled = [ER.compile_mechanism(build(n)) for n in lengths]
        @test ER._concrete(compiled[1]) == build(lengths[1])
        before = instances(decoder)
        foreach(ER._concrete, compiled[2:4])
        @test instances(decoder) == before
    end
end

@testset "a graph written back as raw input builds the same graph" begin
    @testset "$(spec.name)" for spec in MECHANISM_TEST_SPECS
        x = _testhelper_build(spec.raw)
        raw = ER._raw(x)
        @test _testhelper_build(raw) == x
        # Its steps are the edges in stored order and direction, labelled by group.
        c = ER._canonicalize(_testhelper_plain(raw))
        @test c.edge_of == eachindex(c.edge_of) && !any(c.flipped)
        @test c.group_of == 1:ER.ngroups(ER.catalytic(x))
        @test _testhelper_plain(raw).modifiers == ER.modifiers(ER.catalytic(x))
    end
    # Every vector is fresh: editing the raw value leaves the graph as it was.
    am = ER.AllostericMechanism(_testhelper_graph_uniuni(), [:EqualAI, :OnlyA, :OnlyA], 2,
                                [EnzymeRates.RegulatorySite([:R], 2, [:OnlyI])])
    before = deepcopy(am)
    raw = ER._raw(am)
    for v in (raw.catalytic.substrates, raw.catalytic.products, raw.catalytic.modifiers,
              raw.sites[1].ligands)
        push!(v, :X)
    end
    push!(raw.catalytic.group, 4)
    push!(raw.cat_allo_states, :EqualAI)
    raw.sites[1].allo_states[1] = :OnlyA
    empty!(raw.catalytic.steps)
    @test am == before
end

# ── The canonical form pinned as text, and the sub-graph reads against oracles ──

"""The mechanisms the regression census enumerates first: the init populations of its six
reactions, then the children of every init mechanism of its bi-bi reaction."""
function _testhelper_census_population()
    reactions = [
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
    population = Union{EnzymeRates.Mechanism, EnzymeRates.AllostericMechanism}[]
    foreach(r -> append!(population, EnzymeRates.init_mechanisms(r)), reactions)
    bibi = reactions[2]
    append!(population,
            EnzymeRates.expand_mechanisms(EnzymeRates.init_mechanisms(bibi), bibi))
end

"""A step side of graph `m` as text: form `f`'s name, then its free entries `c`."""
_testhelper_side_text(m::ER.Mechanism, f::Int, c::ER.Counts) =
    join([ER._form_name(m, f); [ER._entry_label(m, i) for i in ER.entries(c)]], " + ")

"""Graph `x`'s canonical form as text lines: its names, its forms in order, each edge in
its stored direction with its group (and the group's tag), then any catalytic
multiplicity and sites."""
function _testhelper_canonical_text(x)
    m, allo = ER.catalytic(x), x isa ER.AllostericMechanism
    mods = ER.modifiers(m)
    lines = ["substrates: $(join(m.substrates, ", ")) | products: " *
             "$(join(m.products, ", "))" *
             (isempty(mods) ? "" : " | modifiers: $(join(mods, ", "))"),
             "forms: " * join([ER._form_name(m, f) for f in eachindex(m.forms)], ", ")]
    for (k, e) in enumerate(m.edges)
        g = m.group[k]
        push!(lines, "$k: $(_testhelper_side_text(m, e.from, e.consumed)) " *
                     "$(e.rapid_equilibrium ? "⇌" : "<-->") " *
                     "$(_testhelper_side_text(m, e.to, e.released)) " *
                     "| group $g" * (allo ? " :: $(x.cat_allo_states[g])" : ""))
    end
    allo && push!(lines, "catalytic_multiplicity: $(x.catalytic_multiplicity)")
    for s in (allo ? x.regulatory_sites : EnzymeRates.RegulatorySite[])
        push!(lines, "regulatory_site($(s.multiplicity)): " *
                     join(["$l::$t" for (l, t) in zip(s.ligands, s.allo_states)], ", "))
    end
    lines
end

"""
The inactive edges of the allosteric mechanism whose test reader (`_testhelper_graph`) is
`g`, read off `g`'s steps alone, one flag per step: no step of an `:OnlyA` group, and every
other step whose two forms the steps outside `:OnlyA` groups reach from a free form (no
metabolite, no residual).
"""
function _testhelper_inactive_reach(g)
    outside = [g.tags[s.group] !== :OnlyA for s in g.steps]
    reach = falses(length(g.forms))
    for (j, s) in enumerate(g.steps)
        outside[j] || continue
        g.free[s.from] && (reach[s.from] = true)
        g.free[s.to] && (reach[s.to] = true)
    end
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
Whether the sub-graph reads of `x` (`_segments`, `_inactive_mask`) give what the oracles
read off its test reader's steps give (`_testhelper_re_segments`,
`_testhelper_inactive_reach`): the catalytic graph's segments, form by form and
breadth-first, with each form's positive offsets, and, for an allosteric mechanism, the
inactive conformation's edges and their segments.
"""
function _testhelper_subgraphs_agree(x)
    m, g = ER.catalytic(x), _testhelper_graph(x)
    # Step `j` of the reader is edge `flat[j]`.
    flat = ER._flat_edges(m)
    function _testhelper_segments_agree(mask, kept)
        s = ER._segments(m, mask, ER._flags(m))
        segments, extras = _testhelper_re_segments(g, kept)
        present = findall(!=(0), s.segment_of)
        names = unique!([e.name for e in m.table])
        _testhelper_offsets(f) = Dict(names[n] => s.offsets[n, f]
                                      for n in eachindex(names) if s.offsets[n, f] > 0)
        s.segments == segments && [_testhelper_offsets(f) for f in present] == extras
    end
    _testhelper_segments_agree(trues(length(m.edges)), trues(length(flat))) || return false
    x isa ER.AllostericMechanism || return true
    inactive = ER._inactive_mask(x)
    kept = _testhelper_inactive_reach(g)
    inactive[flat] == kept && _testhelper_segments_agree(inactive, kept)
end

# The census's first populations (`_testhelper_census_population`), enumerated once for
# the two testsets that read them and released after them.
let population = _testhelper_census_population()
@testset "a graph's canonical form, pinned as text, read back from its Sig and raw" begin
    @testset "$(spec.name)" for spec in MECHANISM_TEST_SPECS
        @test ER._concrete(spec.mechanism) == _testhelper_build(spec.raw)
    end
    # The specs' canonical text, a blank line between specs, pinned as a fixture.
    text = reduce((a, b) -> [a; ""; b],
                  [["# $(spec.name)";
                    _testhelper_canonical_text(_testhelper_build(spec.raw))]
                   for spec in MECHANISM_TEST_SPECS])
    fixture = joinpath(@__DIR__, "fixtures", "phase3_canonical_golden.txt")
    if !isfile(fixture)
        mkpath(dirname(fixture)); write(fixture, join(text, "\n") * "\n")
        @warn "bootstrapped phase3 canonical golden fixture; commit it" fixture
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
    x = ER.Mechanism(raw)
    @test x.conformations == [:D, :E] && x.forms[1:2] == [ER.Form(2, 0, 0, 0),
                                                          ER.Form(1, 0, 0, 0)]
    @test _testhelper_graph(x).forms[1:2] == ["E", "D"]
    # The census's first populations: each graph decodes from its Sig and builds again
    # from the raw input it writes back.
    @test length(population) == 5508 + 5505
    @test all(_testhelper_decode(x, ER._sig_of(x)) == x for x in population)
    @test all(_testhelper_build(ER._raw(x)) == x for x in population)
end

@testset "a sub-graph's segments and inactive edges are what its steps give" begin
    @testset "$(spec.name)" for spec in MECHANISM_TEST_SPECS
        @test _testhelper_subgraphs_agree(_testhelper_lift(spec.mechanism))
    end
    @test isempty([i for (i, x) in enumerate(population)
                   if !_testhelper_subgraphs_agree(x)])
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
        @test _testhelper_canonical_text(ER.Mechanism(raw)) == text
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
        @test _testhelper_canonical_text(ER.AllostericMechanism(raw)) == allosteric_text
    end
end

@testset "the constructor builds a spec's graph in under 10 µs, the median over specs" begin
    # Each spec's time per call is the best of five batches; the gate is on their median.
    # The spec's goal is about 2 µs; the gate leaves room for shared CI runners.
    per_call = map(MECHANISM_TEST_SPECS) do spec
        _testhelper_build(spec.raw)
        minimum(@elapsed(for _ in 1:200
                             _testhelper_build(spec.raw)
                         end) / 200 for _ in 1:5)
    end
    @test median(per_call) < 10e-6
end
