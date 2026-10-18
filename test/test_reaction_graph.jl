# ABOUTME: Tests of the reaction-graph mechanism types: content order and identity,
# ABOUTME: canonical form and checks, agreement with today's types, the Sig and speed.

using Statistics: median

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
_testhelper_build(raw::EnzymeRates.RawMechanism) = RG.Mechanism(raw)
_testhelper_build(raw::EnzymeRates.RawAllostericMechanism) = RG.AllostericMechanism(raw)

"""The raw value of `block`, a mechanism macro's block, read by the DSL walker."""
_testhelper_walk(block; allosteric = false) =
    EnzymeRates._raw_mechanism(block, allosteric)

@testset "the constructors build one value from every writing of a mechanism" begin
    @testset "$(spec.name)" for spec in MECHANISM_TEST_SPECS
        x = _testhelper_build(spec.raw)
        @test RG.catalytic(x) == RG._canonicalize(_testhelper_plain(spec.raw)).mechanism
        other = _testhelper_build(_testhelper_rewrite(spec.raw))
        @test other == x && hash(other) == hash(x)
    end
    # Each tag goes to the canonical group its label became.
    am = RG.AllostericMechanism(_testhelper_walk(:(begin
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

@testset "the constructors reject what no graph represents, with today's messages" begin
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
        RG.AllostericMechanism(EnzymeRates.RawAllostericMechanism(
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
    # Tags, by written group label: label 2 is canonical group 3.
    @test_throws(ErrorException("AllostericMechanism: cat_allo_states length 2 must " *
                                "match the number of kinetic groups 3"),
                 _testhelper_allosteric(tags[1:2], 2, EnzymeRates.RawSite[]))
    @test_throws(ErrorException("AllostericMechanism: catalytic_multiplicity must be ≥ " *
                                "1, got 0"),
                 _testhelper_allosteric(tags, 0, EnzymeRates.RawSite[]))
    @test_throws(ErrorException(
        "AllostericMechanism: catalytic group 2 has invalid allo state OnlyI (must be " *
        "one of (:OnlyA, :EqualAI, :NonequalAI)); :OnlyI is rejected for catalytic " *
        "groups (active-state-active convention)"),
        _testhelper_allosteric([:EqualAI, :OnlyI, :EqualAI], 2, EnzymeRates.RawSite[]))
    # Declarations and group labels.
    @test_throws(ErrorException("Mechanism: substrates must not be empty"),
                 RG.Mechanism(_testhelper_uniuni_raw(substrates = Symbol[])))
    @test_throws(ErrorException("Mechanism: products must not be empty"),
                 RG.Mechanism(_testhelper_uniuni_raw(products = Symbol[])))
    @test_throws(ErrorException("Mechanism: duplicate substrate names"),
                 RG.Mechanism(_testhelper_uniuni_raw(substrates = [:S, :S])))
    @test_throws(ErrorException("Mechanism: duplicate product names"),
                 RG.Mechanism(_testhelper_uniuni_raw(products = [:P, :P])))
    @test_throws(ErrorException("Mechanism: duplicate modifier names"),
                 RG.Mechanism(_testhelper_uniuni_raw(modifiers = [:I, :I])))
    @test_throws(ErrorException("Mechanism: S named as both a substrate and a product; " *
                                "concentrations and constants are keyed by name"),
                 RG.Mechanism(_testhelper_uniuni_raw(products = [:P, :S])))
    for group in ([1, 2], [1, 2, 3, 3], [0, 1, 2])
        @test_throws(ErrorException("Mechanism: each step takes one kinetic group label " *
                                    "from 1 to the number of groups"),
                     RG.Mechanism(_testhelper_uniuni_raw(; group)))
    end
    @test_throws(ErrorException("Mechanism: kinetic group 2 is empty; a kinetic group " *
                                "holds at least one step"),
                 RG.Mechanism(_testhelper_uniuni_raw(group = [1, 3, 3])))
    # Names: undeclared; a modifier entry of a name modifiers: does not list; a plain
    # entry of a modifier-only name; residuals by role.
    @test_throws(ErrorException("Mechanism: P is not declared; declare it under " *
                                "substrates:, products: or modifiers:"),
                 RG.Mechanism(_testhelper_uniuni_raw(products = [:Q])))
    mod_P = EnzymeRates.Met(:P, true)
    free = EnzymeRates.RawForm(:E, EnzymeRates.Met[])
    copy_binding = EnzymeRates.RawStep(free, EnzymeRates.RawForm(:E, [mod_P]), [mod_P],
                                       EnzymeRates.Met[], true)
    _testhelper_adding(raw, steps...) = EnzymeRates.RawMechanism(
        raw.substrates, raw.products, raw.modifiers, [raw.steps; steps...],
        [raw.group; length(raw.steps) .+ (1:length(steps))])
    @test_throws(ErrorException("Mechanism: P binds as a modifier, but modifiers: does " *
                                "not declare it"),
                 RG.Mechanism(_testhelper_adding(_testhelper_uniuni_raw(), copy_binding)))
    @test RG.Mechanism(_testhelper_adding(_testhelper_uniuni_raw(modifiers = [:P]),
                                          copy_binding)) isa RG.Mechanism
    I = EnzymeRates.Met(:I, false)
    @test_throws(ErrorException("Mechanism: I is a modifier, so it binds only as a " *
                                "modifier"),
                 RG.Mechanism(_testhelper_adding(
                     _testhelper_uniuni_raw(modifiers = [:I]),
                     EnzymeRates.RawStep(free, EnzymeRates.RawForm(:E, [I]), [I],
                                         EnzymeRates.Met[], true))))
    S = EnzymeRates.Met(:S, false)
    for (added, subtracted, verb, x) in (([:P], Symbol[], "adds", :P),
                                         (Symbol[], [:S], "subtracts", :S))
        residual = EnzymeRates.RawForm(:E, EnzymeRates.Met[], added, subtracted)
        @test_throws(ErrorException("Mechanism: a residual adds substrates and subtracts " *
                                    "products; one $verb $x"),
                     RG.Mechanism(_testhelper_adding(_testhelper_uniuni_raw(),
                         EnzymeRates.RawStep(EnzymeRates.RawForm(:E, [S]), residual,
                                             EnzymeRates.Met[], EnzymeRates.Met[],
                                             false))))
    end
    # Steps: a self-loop, a metabolite taken up and given off, a role change.
    EA = EnzymeRates.RawForm(:E, [EnzymeRates.Met(:A, false)])
    @test_throws(ErrorException("Step: both ends are the form EA; a step joins two forms"),
                 RG.Mechanism(EnzymeRates.RawMechanism([:A], [:P], Symbol[],
                     [EnzymeRates.RawStep(EA, EA, EnzymeRates.Met[], EnzymeRates.Met[],
                                          false)], [1])))
    @test_throws(ErrorException("Step EA → EQ: B both consumed and released"),
                 RG.Mechanism(_testhelper_walk(:(begin
                     substrates: A, B
                     products: Q
                     steps: begin
                         E + A ⇌ E(A)
                         E(A) + B <--> E(Q) + B
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
        @test_throws _testhelper_role_change(from, to) RG.Mechanism(_testhelper_walk(:(begin
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
        RG.Mechanism(_testhelper_walk(:(begin
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
    @test_throws _testhelper_one_name("E(NAD, P)", "E(NADP)", "ENADP") RG.Mechanism(
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
    @test_throws _testhelper_one_name("E(star)", "Estar()", "Estar") RG.Mechanism(
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
    @test_throws _testhelper_one_name("E(P::Mod)", "E(Pmod)", "EPmod") RG.Mechanism(
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
            RG.Mechanism(_testhelper_walk(:(begin
                substrates: A, B
                products: P, Q
                steps: $steps
            end))))
    end
    # One reaction twice: in two groups, with different flags, stored once each way, and in
    # one group.
    _testhelper_twice(steps) = RG.Mechanism(_testhelper_walk(:(begin
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
        RG.Mechanism(_testhelper_walk(:(begin
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
end

@testset "a form's name and text read as constants and the macros write them" begin
    m = RG.Mechanism(_testhelper_walk(:(begin
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
    @test [RG._form_name(m, f) for f in eachindex(m.forms)] ==
          ["E", "E_res_plusA_minusP", "EA", "EB_res_plusA_minusP", "EI", "EQ", "EQQmod"]
    @test [RG._form_text(m, f) for f in eachindex(m.forms)] ==
          ["E()", "E(; residual = A - P)", "E(A)", "E(B; residual = A - P)", "E(I)",
           "E(Q)", "E(Q, Q::Mod)"]
    # Each edge's sides in its stored direction: releases stored as bindings.
    @test [join(RG._sides(m, k), " → ") for k in eachindex(m.edges)] ==
          ["E_A → EA", "E_I → EI", "E_Q → EQ", "E_res_plusA_minusP_P → EA",
           "E_res_plusA_minusP_B → EB_res_plusA_minusP", "EB_res_plusA_minusP → EQ",
           "EQ_Qmod → EQQmod"]
end

"""The graph whose Sig is `sig`, `x`'s kind of graph. `sig` is splatted with its length
left open, as the lifts splat it, so the decoders compile once for every length."""
_testhelper_decode(::RG.Mechanism, @nospecialize(sig)) = RG._mechanism_from_sig(sig...)
_testhelper_decode(::RG.AllostericMechanism, @nospecialize(sig)) =
    RG._allosteric_mechanism_from_sig(sig...)

@testset "a graph's Sig decodes to the graph, and equal graphs share one type" begin
    @testset "$(spec.name)" for spec in MECHANISM_TEST_SPECS
        x = _testhelper_build(spec.raw)
        sig = RG._sig_of(x)
        @test sig isa Tuple{Vararg{Int}}
        @test _testhelper_decode(x, sig) == x
        cx = RG.compile_mechanism(x)
        @test typeof(cx).parameters[1] === sig
        @test RG._concrete(cx) == x
        @test RG.compile_mechanism(_testhelper_build(_testhelper_rewrite(spec.raw))) === cx
    end
    m = _testhelper_graph_uniuni()
    @test RG.Mechanism(RG.EnzymeMechanism(m)) == m
    am = RG.AllostericMechanism(m, [:EqualAI, :NonequalAI, :OnlyA], 4,
                                [EnzymeRates.RegulatorySite([:X, :R], 2, [:OnlyA, :OnlyI])])
    aem = RG.AllostericEnzymeMechanism(am)
    @test RG.AllostericMechanism(aem) == am
    @test RG.catalytic_mechanism(aem) === RG.EnzymeMechanism(m)
    @test RG.catalytic_multiplicity(aem) == 4
    # A name outside ASCII, and table entries past the sixteenth, whose counts fill the
    # high half of a packed count: αKG sorts last, entry 32, the top four bits.
    names = [[Symbol("S", i) for i in 1:30]; :αKG]
    wide = RG.Mechanism(EnzymeRates.RawMechanism(
        names, [:P], Symbol[], _testhelper_bindings([names; :P]), collect(1:32)))
    @test wide.table[32] == EnzymeRates.Met(:αKG, false)
    @test _testhelper_decode(wide, RG._sig_of(wide)) == wide
end

"""Uni-uni whose substrate's name holds `n` characters, so its Sig is longer than any
other test's: a Sig length no other test decodes."""
_testhelper_long_named(n) = let S = Symbol("S", "x"^(n - 1))
    RG.Mechanism(_testhelper_walk(:(begin
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
    # compiles once for every length: after one lift, lifting mechanisms of three Sig
    # lengths no other test decodes adds no method instance. A decoder that took its Sig
    # as an `NTuple{N, Int}` would compile once per length.
    instances(f) = count(_ -> true, Base.specializations(only(methods(f))))
    _testhelper_allosteric_long_named(n) =
        RG.AllostericMechanism(_testhelper_long_named(n), [:EqualAI, :EqualAI, :OnlyA], 2,
                               EnzymeRates.RegulatorySite[])
    for (decoder, build) in ((RG._mechanism_from_sig, _testhelper_long_named),
                             (RG._allosteric_mechanism_from_sig,
                              _testhelper_allosteric_long_named))
        lengths = (301, 401, 501, 601)
        compiled = [RG.compile_mechanism(build(n)) for n in lengths]
        @test RG._concrete(compiled[1]) == build(lengths[1])
        before = instances(decoder)
        foreach(RG._concrete, compiled[2:4])
        @test instances(decoder) == before
    end
end

@testset "a graph written back as raw input builds the same graph" begin
    @testset "$(spec.name)" for spec in MECHANISM_TEST_SPECS
        x = _testhelper_build(spec.raw)
        raw = RG._raw(x)
        @test _testhelper_build(raw) == x
        # Its steps are the edges in stored order and direction, labelled by group.
        c = RG._canonicalize(_testhelper_plain(raw))
        @test c.edge_of == eachindex(c.edge_of) && !any(c.flipped)
        @test c.group_of == 1:RG.ngroups(RG.catalytic(x))
        @test _testhelper_plain(raw).modifiers == RG.modifiers(RG.catalytic(x))
    end
    # Every vector is fresh: editing the raw value leaves the graph as it was.
    am = RG.AllostericMechanism(_testhelper_graph_uniuni(), [:EqualAI, :OnlyA, :OnlyA], 2,
                                [EnzymeRates.RegulatorySite([:R], 2, [:OnlyI])])
    before = deepcopy(am)
    raw = RG._raw(am)
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

# ── The graph against today's types ─────────────────────────────────────────

"""Metabolite `x` of today's types as raw input: its name, and whether it binds as a
modifier (a competitive-inhibitor copy)."""
_testhelper_met(x) = EnzymeRates.Met(EnzymeRates.name(x),
                                     x isa EnzymeRates.CompetitiveInhibitor)

"""Form `sp` of today's types as raw input."""
_testhelper_raw_form(sp) = EnzymeRates.RawForm(EnzymeRates.conformation(sp),
    EnzymeRates.Met[_testhelper_met(x) for x in EnzymeRates.bound(sp)],
    Symbol[EnzymeRates.name(x) for x in EnzymeRates.added(EnzymeRates.residual(sp))],
    Symbol[EnzymeRates.name(x) for x in EnzymeRates.subtracted(EnzymeRates.residual(sp))])

"""The steps of today's mechanism `m` as raw input, written against the way `m` stores
them: in reverse order, each from its other end, each labelled with its kinetic group."""
function _testhelper_raw_steps(m)
    flat = reverse(EnzymeRates._flat_steps(EnzymeRates.steps(m)))
    EnzymeRates.RawMechanism(EnzymeRates.substrates(m), EnzymeRates.products(m),
        EnzymeRates.modifiers(m),
        [EnzymeRates.RawStep(_testhelper_raw_form(EnzymeRates.to_species(s)),
                             _testhelper_raw_form(EnzymeRates.from_species(s)),
                             _testhelper_met.(EnzymeRates.released(s)),
                             _testhelper_met.(EnzymeRates.consumed(s)),
                             EnzymeRates.is_equilibrium(s)) for (s, _) in flat],
        [g for (_, g) in flat])
end

"""Today's mechanism `m` as raw input (`_testhelper_raw_steps`); an allosteric one with
its tags, catalytic multiplicity and sites."""
_testhelper_raw_of(m::EnzymeRates.Mechanism) = _testhelper_raw_steps(m)
_testhelper_raw_of(am::EnzymeRates.AllostericMechanism) =
    EnzymeRates.RawAllostericMechanism(_testhelper_raw_steps(am),
        copy(EnzymeRates.cat_allo_states(am)), EnzymeRates.catalytic_multiplicity(am),
        [EnzymeRates.RawSite(copy(EnzymeRates.ligands(s)), EnzymeRates.multiplicity(s),
                             copy(EnzymeRates.allo_states(s)))
         for s in EnzymeRates.regulatory_sites(am)])

"""Today's mechanism `m` (concrete or compiled) as `_testhelper_graph` reads it, without
the constant names, which the graph does not render yet, and with its concentration
names."""
function _testhelper_today_view(m)
    g = _testhelper_graph(m)
    (forms = g.forms, free = g.free,
     steps = [(; s.from, s.to, s.consumed, s.released, s.consumed_modifier,
               s.released_modifier, s.rapid_equilibrium, s.group) for s in g.steps],
     g.tags, g.catalytic_multiplicity, g.sites, g.substrates, g.products, g.modifiers,
     metabolites = EnzymeRates.metabolites(_testhelper_lift(m)))
end

"""Edge `k` of graph `m` as `_testhelper_graph` reads a step."""
function _testhelper_edge_view(m::RG.Mechanism, k::Int)
    e = m.edges[k]
    c, r = RG.entries(e.consumed), RG.entries(e.released)
    (from = e.from, to = e.to, consumed = Symbol[m.table[i].name for i in c],
     released = Symbol[m.table[i].name for i in r],
     consumed_modifier = Bool[m.table[i].modifier for i in c],
     released_modifier = Bool[m.table[i].modifier for i in r],
     rapid_equilibrium = e.rapid_equilibrium, group = m.group[k])
end

"""Graph `x` read as `_testhelper_today_view` reads today's mechanism: its steps group
by group, each group in edge order."""
function _testhelper_graph_view(x)
    m, allo = RG.catalytic(x), x isa RG.AllostericMechanism
    (forms = [RG._form_name(m, f) for f in eachindex(m.forms)],
     free = Bool[f.bound == 0 && f.added == 0 && f.subtracted == 0 for f in m.forms],
     steps = [_testhelper_edge_view(m, k) for members in RG._group_members(m)
              for k in members],
     tags = allo ? RG.cat_allo_states(x) : Symbol[],
     catalytic_multiplicity = allo ? RG.catalytic_multiplicity(x) : nothing,
     sites = [(multiplicity = s.multiplicity, ligands = s.ligands,
               allo_states = s.allo_states)
              for s in (allo ? RG.regulatory_sites(x) : EnzymeRates.RegulatorySite[])],
     substrates = RG.substrates(m), products = RG.products(m),
     modifiers = RG.modifiers(m), metabolites = RG.metabolites(x))
end

"""Today's mechanisms the regression census enumerates first: the init populations of
its six reactions, then the children of every init mechanism of its bi-bi reaction."""
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
_testhelper_side_text(m::RG.Mechanism, f::Int, c::RG.Counts) =
    join([RG._form_name(m, f); [RG._entry_label(m, i) for i in RG.entries(c)]], " + ")

"""Graph `x`'s canonical form as text lines: its names, its forms in order, each edge in
its stored direction with its group (and the group's tag), then any catalytic
multiplicity and sites."""
function _testhelper_canonical_text(x)
    m, allo = RG.catalytic(x), x isa RG.AllostericMechanism
    mods = RG.modifiers(m)
    lines = ["substrates: $(join(m.substrates, ", ")) | products: " *
             "$(join(m.products, ", "))" *
             (isempty(mods) ? "" : " | modifiers: $(join(mods, ", "))"),
             "forms: " * join([RG._form_name(m, f) for f in eachindex(m.forms)], ", ")]
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

@testset "the graph agrees with today's canonical form, in order" begin
    @testset "$(spec.name)" for spec in MECHANISM_TEST_SPECS
        @test _testhelper_graph_view(_testhelper_build(spec.raw)) ==
              _testhelper_today_view(spec.mechanism)
    end
    # The specs' canonical text, a blank line between specs, kept so the comparison
    # outlives today's types.
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
    x = RG.Mechanism(raw)
    @test x.conformations == [:D, :E] && x.forms[1:2] == [RG.Form(2, 0, 0, 0),
                                                          RG.Form(1, 0, 0, 0)]
    @test _testhelper_graph_view(x) == _testhelper_today_view(EnzymeRates.Mechanism(raw))
    # The census's first populations, each written against its stored form; in
    # agreement, the graph also decodes from its Sig and from its raw form.
    population = _testhelper_census_population()
    graphs = [_testhelper_build(_testhelper_raw_of(m)) for m in population]
    disagree = [i for (i, (m, x)) in enumerate(zip(population, graphs))
                if _testhelper_graph_view(x) != _testhelper_today_view(m)]
    @test isempty(disagree)
    @test length(population) == 5508 + 5505
    @test length(unique(graphs)) == length(unique(population))
    @test all(_testhelper_decode(x, RG._sig_of(x)) == x for x in graphs)
    @test all(_testhelper_build(RG._raw(x)) == x for x in graphs)
end

"""Today's mechanism that `raw` writes, built by today's constructors."""
_testhelper_today(raw::EnzymeRates.RawMechanism) = EnzymeRates.Mechanism(raw)
_testhelper_today(raw::EnzymeRates.RawAllostericMechanism) =
    EnzymeRates.AllostericMechanism(raw)

@testset "the graph orients as today where tier 3 and tied groups decide" begin
    # No spec and no census population holds a step that only form order orients or a
    # tied group whose alignment turns a step, so these cases are compared here, each in
    # two writings.
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
    for raw in _testhelper_walk.(blocks)
        @test _testhelper_graph_view(RG.Mechanism(raw)) ==
              _testhelper_today_view(_testhelper_today(raw))
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
        @test _testhelper_graph_view(RG.AllostericMechanism(raw)) ==
              _testhelper_today_view(_testhelper_today(raw))
    end
end

@testset "the constructor builds a spec's graph in under 10 µs, the median over specs" begin
    # Each spec's time per call is the best of five batches; the gate is on their median
    # (R17: the design aims at about 2 µs and leaves room for shared CI runners).
    per_call = map(MECHANISM_TEST_SPECS) do spec
        _testhelper_build(spec.raw)
        minimum(@elapsed(for _ in 1:200
                             _testhelper_build(spec.raw)
                         end) / 200 for _ in 1:5)
    end
    @test median(per_call) < 10e-6
end
