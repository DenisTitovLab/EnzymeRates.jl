# ABOUTME: Tests of the reaction-graph mechanism types: content order and identity,
# ABOUTME: canonical form and checks, agreement with today's types, the Sig and speed.

const RG = EnzymeRates.ReactionGraph

"""The counts of `list`, table entries listed with repeats."""
_testhelper_pack(list) =
    foldl((c, i) -> c + (RG.Counts(1) << (4 * (i - 1))), list; init = RG.Counts(0))

"""Uni-uni with rapid-equilibrium bindings, written field by field in canonical form:
forms E, E(P), E(S); the bindings of P and S, then E(S) → E(P)."""
function _testhelper_graph_uniuni()
    P, S = _testhelper_pack([1]), _testhelper_pack([2])
    RG.Mechanism([:S], [:P], [EnzymeRates.Met(:P, false), EnzymeRates.Met(:S, false)],
                 [:E], [RG.Form(1, 0, 0, 0), RG.Form(1, P, 0, 0), RG.Form(1, S, 0, 0)],
                 [RG.Edge(1, 2, P, 0, true), RG.Edge(1, 3, S, 0, true),
                  RG.Edge(3, 2, 0, 0, false)],
                 [1, 2, 3])
end

@testset "packed counts list the table entries of a form or a step side" begin
    c = _testhelper_pack([1, 1, 3])
    @test (RG.copies(c, 1), RG.copies(c, 2), RG.copies(c, 3)) == (2, 0, 1)
    @test RG.total(c) == 3
    @test RG.entries(c) == [1, 1, 3]
    @test RG.entries(RG.Counts(0)) == Int[] && RG.total(RG.Counts(0)) == 0
    # The last of 32 entries, 15 copies: the top four bits.
    top = RG.Counts(15) << 124
    @test RG.copies(top, 32) == 15 && RG.total(top) == 15 && RG.entries(top) == fill(32, 15)
end

@testset "entry lists order lexicographically, a list before every list it begins" begin
    # Every list over four entries with at most two copies of each, against Julia's order
    # of the listed entries; the packed numbers order differently ([1, 2] < [2], 17 > 16).
    lists = [reduce(vcat, [fill(i, n[i]) for i in 1:4]; init = Int[])
             for n in Iterators.product(0:2, 0:2, 0:2, 0:2)]
    @test all(RG._seq_less(_testhelper_pack(a), _testhelper_pack(b)) == isless(a, b)
              for a in lists, b in lists)
    @test _testhelper_pack([1, 2]) > _testhelper_pack([2])
end

@testset "forms order by label, bound, added, subtracted; edges by ends, lists, flag" begin
    A, B, AB = _testhelper_pack([1]), _testhelper_pack([2]), _testhelper_pack([1, 2])
    @test RG.Form(1, B, 0, 0) < RG.Form(2, 0, 0, 0)
    @test RG.Form(1, AB, 0, 0) < RG.Form(1, B, 0, 0)
    @test RG.Form(1, B, 0, 0) < RG.Form(1, B, A, 0) < RG.Form(1, B, AB, 0)
    @test RG.Form(1, B, A, 0) < RG.Form(1, B, A, B)
    @test !(RG.Form(1, B, A, B) < RG.Form(1, B, A, B))
    @test RG.Edge(1, 3, 0, 0, false) < RG.Edge(2, 1, 0, 0, false)
    @test RG.Edge(1, 2, 0, 0, false) < RG.Edge(1, 3, 0, 0, false)
    @test RG.Edge(1, 2, AB, 0, true) < RG.Edge(1, 2, B, 0, false)
    @test RG.Edge(1, 2, A, AB, true) < RG.Edge(1, 2, A, B, false)
    @test RG.Edge(1, 2, A, 0, false) < RG.Edge(1, 2, A, 0, true)
    @test RG._reversed(RG.Edge(1, 2, A, B, true)) == RG.Edge(2, 1, B, A, true)
end

"""Whether a value of type `T` hashes the same in every process: a symbol, an integer or
a flag, a vector of such values, or one of the types whose hash folds its name and its
fields (`_testhelper_hash_folds`), whose fields are such values in turn."""
function _testhelper_hash_is_content(T)
    T in (Symbol, Bool, Int, UInt128) && return true
    T <: Vector && return _testhelper_hash_is_content(eltype(T))
    T in (EnzymeRates.Met, EnzymeRates.RegulatorySite, RG.Form, RG.Edge, RG.Mechanism,
          RG.AllostericMechanism) || return false
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
    am = RG.AllostericMechanism(m, [:EqualAI, :EqualAI, :OnlyA], 2, [site])
    h = UInt(17)
    for x in (m.forms[2], m.edges[3], m, am, m.table[1], site)
        @test _testhelper_hash_folds(x, h)
    end
    @test all(_testhelper_hash_is_content,
              (RG.Form, RG.Edge, RG.Mechanism, RG.AllostericMechanism))
    # Each field takes part in `==` and `hash`.
    @test m == _testhelper_graph_uniuni() && hash(m) == hash(_testhelper_graph_uniuni())
    for f in fieldnames(RG.Mechanism)
        fields = Any[getfield(m, g) for g in fieldnames(RG.Mechanism)]
        i = findfirst(==(f), fieldnames(RG.Mechanism))
        fields[i] = f === :group ? [1, 3, 2] : f === :substrates ? [:A] :
                    f === :products ? [:Q] : f === :conformations ? [:F] :
                    f === :table ? reverse(m.table) : reverse(getfield(m, f))
        other = RG.Mechanism(fields...)
        @test other != m && hash(other) != hash(m)
    end
    @test am == RG.AllostericMechanism(m, [:EqualAI, :EqualAI, :OnlyA], 2, [site])
    @test am != RG.AllostericMechanism(m, [:EqualAI, :OnlyA, :OnlyA], 2, [site])
    @test am != RG.AllostericMechanism(m, [:EqualAI, :EqualAI, :OnlyA], 4, [site])
    @test am != RG.AllostericMechanism(m, [:EqualAI, :EqualAI, :OnlyA], 2,
                                       EnzymeRates.RegulatorySite[])
end

@testset "an allosteric graph sorts its sites and checks its tags and ligands" begin
    m = _testhelper_graph_uniuni()
    tags = [:EqualAI, :EqualAI, :OnlyA]
    sites = [EnzymeRates.RegulatorySite([:X], 2, [:OnlyI]),
             EnzymeRates.RegulatorySite([:B, :A], 2, [:OnlyA, :EqualAI])]
    am = RG.AllostericMechanism(m, tags, 2, sites)
    @test [s.ligands for s in am.regulatory_sites] == [[:A, :B], [:X]]
    @test am.regulatory_sites[1].allo_states == [:EqualAI, :OnlyA]
    tags[1] = :OnlyA
    @test am.cat_allo_states == [:EqualAI, :EqualAI, :OnlyA]
    @test_throws(ErrorException("AllostericMechanism: cat_allo_states length 2 must " *
                                "match the number of kinetic groups 3"),
                 RG.AllostericMechanism(m, [:EqualAI, :EqualAI], 2, sites))
    @test_throws(ErrorException("AllostericMechanism: catalytic_multiplicity must be ≥ " *
                                "1, got 0"),
                 RG.AllostericMechanism(m, [:EqualAI, :EqualAI, :OnlyA], 0, sites))
    @test_throws(ErrorException(
        "AllostericMechanism: catalytic group 2 has invalid allo state OnlyI (must be " *
        "one of (:OnlyA, :EqualAI, :NonequalAI)); :OnlyI is rejected for catalytic " *
        "groups (active-state-active convention)"),
        RG.AllostericMechanism(m, [:EqualAI, :OnlyI, :OnlyA], 2, sites))
    @test_throws(ErrorException(
        "AllostericMechanism: ligand X appears in two distinct regulatory sites; " *
        "rendered Kreg names would collide. Same-ligand-two-sites is not enumerated."),
        RG.AllostericMechanism(m, [:EqualAI, :EqualAI, :OnlyA], 2,
                               [sites; EnzymeRates.RegulatorySite([:X], 4, [:OnlyA])]))
    @test_throws(ErrorException(
        "AllostericMechanism: L is a reserved name; E_total, Keq, L, N, params, concs " *
        "and v are kept for a rate equation's constants and local names"),
        RG.AllostericMechanism(m, [:EqualAI, :EqualAI, :OnlyA], 2,
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
        m = RG._canonicalize(raw).mechanism
        for (order, direction, labels) in ((true, false, false), (false, true, false),
                                           (false, false, true), (true, true, true))
            other = RG._canonicalize(
                _testhelper_rewrite(raw; order, direction, labels)).mechanism
            @test other == m && hash(other) == hash(m)
        end
    end
end

@testset "canonical form maps each written step to its edge, group and forms" begin
    c = RG._canonicalize(EnzymeRates._raw_mechanism(:(begin
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
        RG._canonicalize(EnzymeRates._raw_mechanism(block, false)).flipped
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
    c, d = RG._canonicalize(_testhelper_mirrored(:(F(I) <--> E(I)))),
           RG._canonicalize(_testhelper_mirrored(:(E(I) <--> F(I))))
    @test c.flipped == [false, false, true, false, false, false, false]
    @test d.flipped == [false, false, true, false, true, false, false]
    @test c.mechanism == d.mechanism
    # Tied steps that change different pairs of conformations: no orientation aligns
    # them, and canonical form returns the lead and the first that misses it.
    c = RG._canonicalize(EnzymeRates._raw_mechanism(:(begin
        substrates: A
        products: P
        steps: begin
            E + A <--> E(A)
            E(A) <--> E(P)
            E(P) <--> E + P
            (E ⇌ F, F ⇌ G, G ⇌ E)
        end
    end), false))
    @test c.misaligned == [RG.Edge(1, 4, 0, 0, true), RG.Edge(4, 5, 0, 0, true)]
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
    @test length(RG._canonicalize(_testhelper_raw(32)).mechanism.table) == 32
    @test_throws(ErrorException("Mechanism: a mechanism holds at most 32 metabolite " *
                                "entries (a name and its modifier copy count apart); " *
                                "this one has 33"),
                 RG._canonicalize(_testhelper_raw(33)))
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
                 RG._canonicalize(both))
    _testhelper_holding(n) = EnzymeRates.RawMechanism([:A], [:P], Symbol[],
        [EnzymeRates.RawStep(EnzymeRates.RawForm(:E, fill(EnzymeRates.Met(:A, false), n)),
                             EnzymeRates.RawForm(:F, fill(EnzymeRates.Met(:A, false), n)),
                             EnzymeRates.Met[], EnzymeRates.Met[], false)], [1])
    @test RG.copies(
        RG._canonicalize(_testhelper_holding(15)).mechanism.forms[1].bound, 1) == 15
    @test_throws(ErrorException("Mechanism: a form or a step side holds A more than 15 " *
                                "times; a mechanism holds at most 15 copies of one " *
                                "metabolite there"),
                 RG._canonicalize(_testhelper_holding(16)))
end
