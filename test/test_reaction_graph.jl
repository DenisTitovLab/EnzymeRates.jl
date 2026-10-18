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
