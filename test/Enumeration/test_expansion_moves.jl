# ABOUTME: Tests of src/Enumeration/expansion_moves.jl: the seven moves, their searches,
# ABOUTME: the emission rules and candidate_mechanisms(reaction, parent).

# The names of the forms of graph `m` where `name` (its modifier copy when `modifier`
# holds) is free on some edge (`_forms_where_free`).
_testhelper_free_at(m, name, modifier = false) =
    Set(String(m.catalytic.enzyme_forms[f])
        for f in findall(EnzymeRates.Enumeration._forms_where_free(_testhelper_candidate(m),
                                                                   name, modifier)))

# The number of rapid-equilibrium segments of the catalytic graph of `m`.
_testhelper_segment_count(m) =
    (cm = m.catalytic;
     EnzymeRates.Enumeration._re_segment_count(
         cm, _testhelper_flags(cm)))

# Sorted fitted-parameter count of each child minus its parent's.
function _testhelper_param_deltas(parent, kids)
    base = _testhelper_fitted(parent)
    sort([_testhelper_fitted(k) - base for k in kids])
end

# The names that the forms of the enumerated mechanism `m` bind as modifier copies, a pure
# inhibitor among them, read off the records of its forms (`_testhelper_candidate`).
_testhelper_modifier_names(m) =
    Set(x for f in _testhelper_candidate(m).forms for x in f.copies)

# Whether `child` keeps the reaction of `parent`: the same metabolites as substrates and
# products, and every name `parent` binds as a modifier copy still bound so.
_testhelper_keeps_reaction(child, parent) =
    child.catalytic.metabolites[child.catalytic.substrates] ==
        parent.catalytic.metabolites[parent.catalytic.substrates] &&
    child.catalytic.metabolites[child.catalytic.products] ==
        parent.catalytic.metabolites[parent.catalytic.products] &&
    _testhelper_modifier_names(parent) ⊆ _testhelper_modifier_names(child)

# The steps of mechanism or catalytic scheme `m` apart from their kinetic groups: what a
# regrouping move keeps.
_testhelper_ungrouped(m::EnzymeRates.Mechanism) = _testhelper_ungrouped(m.catalytic)
_testhelper_ungrouped(c::EnzymeRates.Mechanisms.CatalyticScheme) =
    [(e.from, e.to, e.taken_up, e.given_off, e.rapid_equilibrium) for e in c.steps]

# The kinetic groups of `m` whose steps are isomerizations.
_testhelper_iso_groups(m) =
    unique(s.group for s in _testhelper_graph(m).steps if _testhelper_is_iso(s))

# The kinetic-group index of every step of `m` that binds the metabolite named `x`.
_testhelper_groups_binding(m, x) = [s.group for s in _testhelper_graph(m).steps
                                    if _testhelper_is_binding(s) && only(s.consumed) == x]

# A kinetic group is catalytic iff it holds a chemistry step (anything but a plain
# binding: an isomerization, a fused binding or a Theorell–Chance step); otherwise it is a
# binding group. Mirrors `_chemistry_groups`, which `_expand_to_allosteric` reads to choose
# the sets of `:OnlyA` groups it pairs with a regulator.
_testhelper_is_catalytic_group(m, g) =
    (c = _testhelper_candidate(m);
     any(e -> e.kinetic_group == g &&
              (!EnzymeRates.Enumeration._is_binding(e) || _testhelper_is_fused(c, e)),
         m.catalytic.steps))

const uni_uni_allo = @enzyme_reaction begin
    substrates: S[C]
    products: P[C]
    oligomeric_state: 2
end

const uni_uni_allo_2reg = @enzyme_reaction begin
    substrates: S[C]
    products: P[C]
    allosteric_regulators: R1, R2
    oligomeric_state: 2
end

# One allosteric regulator (R2) and one competitive inhibitor (R1): the two
# regulator kinds must stay in their own expansion moves — R2 goes to a
# regulatory site, R1 to a dead-end binding.
const uni_uni_reg_and_inhibitor = @enzyme_reaction begin
    substrates: S[C]
    products: P[C]
    allosteric_regulators: R2
    competitive_inhibitors: R1
    oligomeric_state: 2
end

"""
    _testhelper_enumerate_all_mechanism(rxn::EnzymeReaction; max_params::Int=typemax(Int))
        -> Dict{Int, Vector{Mechanism}}

Enumerate all mechanisms reachable by init → expand, bucketed by ACTUAL
fitted-parameter count. `max_params` caps the search: mechanisms whose
actual fitted count exceeds it are dropped. Uses the same advancing-target
sweep as `_beam_search` (but expands every swept mechanism — no beam
selection) so Δ=0 expansion children (same param count as the parent) are
not lost.
"""
function _testhelper_enumerate_all_mechanism(rxn; max_params::Int=typemax(Int))
    M = EnzymeRates.Mechanism
    actual(m) = _testhelper_fitted(m)
    frontier = Dict{Int, Vector{M}}()
    function _testhelper_add!(m)
        pc = actual(m)
        pc <= max_params && push!(get!(frontier, pc, M[]), m)
    end
    for m in unique!(
        collect(_testhelper_mechanisms(EnzymeRates.Enumeration.init_mechanisms(rxn))))
        _testhelper_add!(m)
    end
    results = Dict{Int, Vector{M}}()
    isempty(frontier) && return results
    target = minimum(keys(frontier))
    while !isempty(frontier)
        swept = M[]
        for c in collect(keys(frontier))
            c <= target && append!(swept, pop!(frontier, c))
        end
        swept = unique!(swept)
        for m in swept
            push!(get!(results, actual(m), M[]), m)
        end
        for child in _testhelper_children(rxn, swept)
            _testhelper_add!(child)
        end
        isempty(frontier) && break
        target = max(target + 1, minimum(keys(frontier)))
    end
    results
end

# ─── _inhibitor_competition_patterns ────────────────────────────────────
@testset "_inhibitor_competition_patterns" begin
# Uni-uni, no existing inhibitors
pats = EnzymeRates.Enumeration._inhibitor_competition_patterns(
    Set([:S]), Set([:P]), Symbol[])
@test length(pats) == 1
@test pats[1] == (Set([:S]), Set([:P]), Set{Symbol}())

# Each existing inhibitor doubles the inhibitor-competition combinations: k existing
# inhibitors give 2^k, combined with the base patterns (bi-bi 3×3 = 9, ter-ter
# 7×7 = 49).
for (S, P, existing, n) in [
    ([:A, :B], [:P, :Q], Symbol[], 9),
    ([:A, :B, :C], [:P, :Q, :R], Symbol[], 49),
    ([:A, :B], [:P, :Q], [:I1__reg], 18),                        # 9 × 2
    ([:A, :B], [:P, :Q], [:I1__reg, :I2__reg], 36),              # 9 × 4
    ([:A, :B], [:P, :Q], [:I1__reg, :I2__reg, :I3__reg], 72),    # 9 × 8
]
    pats = EnzymeRates.Enumeration._inhibitor_competition_patterns(Set(S), Set(P), existing)
    @test length(pats) == n
end
end

# ─── _forms_where_free ──────────────────────────────────────────────────
@testset "_forms_where_free" begin
@testset "_forms_where_free" begin
# Uni-uni: S binds to E, P binds to E
m_uu = @enzyme_mechanism begin
    substrates: S
    products: P
    steps: begin
        E + P ⇌ EP
        E + S ⇌ ES
        ES <--> EP
    end
end
mech_uu = m_uu
@test _testhelper_free_at(mech_uu, :S) == Set(["E"])
@test _testhelper_free_at(mech_uu, :P) == Set(["E"])

# Bi-bi random: B binds to E and E_A
m_bb = @enzyme_mechanism begin
    substrates: A, B
    products: P, Q
    steps: begin
        E + A ⇌ EA
        EB + A ⇌ EAB
        E + B ⇌ EB
        EA + B ⇌ EAB
        E + P ⇌ EP
        EP + Q ⇌ EPQ
        E + Q ⇌ EQ
        EQ + P ⇌ EPQ
        EAB <--> EPQ
    end
end
mech_bb = m_bb
@test _testhelper_free_at(mech_bb, :B) == Set(["E", "EA"])
@test _testhelper_free_at(mech_bb, :A) == Set(["E", "EB"])
@test _testhelper_free_at(mech_bb, :P) == Set(["E", "EQ"])
@test _testhelper_free_at(mech_bb, :Q) == Set(["E", "EP"])
end

@testset "a copy binds at its own sites, apart from its reactant's" begin
    # Ordered bi-bi with A as its own inhibitor at {E, EQ}: A binds productively
    # at E, its copy at E and EQ; Q binds at E and at the copy's form E(A*).
    m = @enzyme_mechanism begin
        substrates: A, B
        products: P, Q
        steps: begin
            E + A ⇌ EA
            EA + B ⇌ EAB
            EAB <--> EPQ
            EQ + P ⇌ EPQ
            (E + Q ⇌ EQ, EAmod + Q ⇌ EAmodQ)
            (E + A ⇌ EAmod, EQ + A ⇌ EAmodQ)
        end
    end
    @test _testhelper_free_at(m, :A) == Set(["E"])
    @test _testhelper_free_at(m, :A, true) == Set(["E", "EQ"])
    @test _testhelper_free_at(m, :Q) == Set(["E", "EAmod"])
    @test isempty(_testhelper_free_at(m, :Q, true))
end

@testset "a Theorell–Chance step frees B at EA and P at EQ" begin
    # EA + B → EQ + P takes B up at EA and gives P off at EQ; A and Q
    # bind at E.
    m = @enzyme_mechanism begin
        substrates: A, B
        products: P, Q
        steps: begin
            E + A <--> EA
            EA + B <--> EQ + P
            E + Q <--> EQ
        end
    end
    @test _testhelper_free_at(m, :A) == Set(["E"])
    @test _testhelper_free_at(m, :B) == Set(["EA"])
    @test _testhelper_free_at(m, :P) == Set(["EQ"])
    @test _testhelper_free_at(m, :Q) == Set(["E"])
end
end

# ─── _expand_re_to_ss ──────────────────────────────────────────────────
@testset "_expand_re_to_ss" begin

@testset "Mechanism — flip past a bottomless RE segment" begin
    # PARENT: uni-bi with random product release, all binding RE. The product
    # ring E–EP–EPQ–EQ–E needs two cuts to raise the segment count, so
    # every pair of ring edges gains. Cutting the two edges at E ({P@E, Q@E})
    # would leave {EP, EQ, EPQ} with no bottom form, which class 3 rejects at
    # zero products: the flip move's guard (`check_valid_rate_equation`) fails that pair,
    # so it counts as failed and is not emitted.
    # Its supersets are reached by one more flip from the other pairs' children.
    # The two pairs at EP and at EQ isolate one form joined to the rest by two
    # edges of stored weights −1 and +1, whose 2-cycle sums to 0 (a balanced block):
    # neither flipped group carries flux, so they fail the same way and keep their
    # parent's rank. Children: the S flip plus the three other pairs.
    parent = @enzyme_mechanism begin
        substrates: S
        products: P, Q
        steps: begin
            E + S ⇌ ES
            E + P ⇌ EP
            E + Q ⇌ EQ
            EP + Q ⇌ EPQ
            EQ + P ⇌ EPQ
            ES <--> EPQ
        end
    end
    children = _testhelper_move(EnzymeRates.Enumeration._expand_re_to_ss, parent)
    expected = [
        (@enzyme_mechanism begin      # S alone: ES becomes its own segment
            substrates: S
            products: P, Q
            steps: begin
                E + S <--> ES
                E + P ⇌ EP
                E + Q ⇌ EQ
                EP + Q ⇌ EPQ
                EQ + P ⇌ EPQ
                ES <--> EPQ
            end
        end),
        (@enzyme_mechanism begin      # the two edges at EPQ
            substrates: S
            products: P, Q
            steps: begin
                E + S ⇌ ES
                E + P ⇌ EP
                E + Q ⇌ EQ
                EP + Q <--> EPQ
                EQ + P <--> EPQ
                ES <--> EPQ
            end
        end),
        (@enzyme_mechanism begin      # P@E and P@EQ
            substrates: S
            products: P, Q
            steps: begin
                E + S ⇌ ES
                E + P <--> EP
                E + Q ⇌ EQ
                EP + Q ⇌ EPQ
                EQ + P <--> EPQ
                ES <--> EPQ
            end
        end),
        (@enzyme_mechanism begin      # Q@E and Q@EP
            substrates: S
            products: P, Q
            steps: begin
                E + S ⇌ ES
                E + P ⇌ EP
                E + Q <--> EQ
                EP + Q <--> EPQ
                EQ + P ⇌ EPQ
                ES <--> EPQ
            end
        end)
    ]
    @test length(children) == 4
    @test Set(children) == Set(expected)
    for (i, c) in enumerate(children)
        @test all(==(:agree),
                  _testhelper_respeed_verdicts(c; rng = Random.MersenneTwister(1302 + i)))
    end
    # The pair {P@E, Q@E} fails class 3 at zero products.
    pq = [_testhelper_edge(parent, "E + P ⇌ EP"), _testhelper_edge(parent, "E + Q ⇌ EQ")]
    pair = _testhelper_flags(parent)
    pair[pq] .= false
    v = EnzymeRates.Derivation._guard(
        EnzymeRates.Enumeration._with_flags(_testhelper_candidate(parent), pair).mechanism)
    @test v.fault === EnzymeRates.Derivation._BOTTOMLESS
    # The two zero-flux pairs are absent, and each has exactly its parent's rank.
    r0 = _testhelper_identifiable_rank(parent)
    for absent in ([
        (@enzyme_mechanism begin      # P@E and Q@EP
            substrates: S
            products: P, Q
            steps: begin
                E + S ⇌ ES
                E + P <--> EP
                E + Q ⇌ EQ
                EP + Q <--> EPQ
                EQ + P ⇌ EPQ
                ES <--> EPQ
            end
        end),
        (@enzyme_mechanism begin      # Q@E and P@EQ
            substrates: S
            products: P, Q
            steps: begin
                E + S ⇌ ES
                E + P ⇌ EP
                E + Q <--> EQ
                EP + Q ⇌ EPQ
                EQ + P <--> EPQ
                ES <--> EPQ
            end
        end)
    ])
        @test !(absent in children)
        @test _testhelper_segment_count(absent) > _testhelper_segment_count(parent)
        @test _testhelper_identifiable_rank(absent) == r0
    end
end

@testset "Mechanism — the flip skips the flanks of a qualifying chain" begin
    # Ordered RE seed: E + A ⇌ EA, EA + B ⇌ EAB, EAB → EPQ (SS),
    # EQ + P ⇌ EPQ, E + Q ⇌ EQ. EAB and EPQ each have one other step, a
    # binding into them, and every group holds one step: a qualifying chain (RSR). Its
    # flanks, the B and P groups, never flip: a flipped flank gives its parent's family
    # plus a phantom (chain lemma). The A and Q flips each cut a segment and carry flux.
    seed = @enzyme_mechanism begin
        substrates: A, B; products: P, Q
        steps: begin
            E + A ⇌ EA
            EA + B ⇌ EAB
            EAB <--> EPQ
            EQ + P ⇌ EPQ
            E + Q ⇌ EQ
        end
    end
    a_flip = @enzyme_mechanism begin
        substrates: A, B; products: P, Q
        steps: begin
            E + A <--> EA
            EA + B ⇌ EAB
            EAB <--> EPQ
            EQ + P ⇌ EPQ
            E + Q ⇌ EQ
        end
    end
    q_flip = @enzyme_mechanism begin
        substrates: A, B; products: P, Q
        steps: begin
            E + A ⇌ EA
            EA + B ⇌ EAB
            EAB <--> EPQ
            EQ + P ⇌ EPQ
            E + Q <--> EQ
        end
    end
    kids = _testhelper_move(EnzymeRates.Enumeration._expand_re_to_ss, seed)
    @test length(kids) == 2
    @test Set(kids) == Set([a_flip, q_flip])
    # The seed has 5 fitted constants, all identifiable; each child adds one and
    # gains one. Each flank flip, absent, adds a constant and no rank.
    @test _testhelper_fitted(seed) == _testhelper_identifiable_rank(seed) == 5
    for k in kids
        @test _testhelper_fitted(k) == _testhelper_identifiable_rank(k) == 6
    end
    for absent in ([
        (@enzyme_mechanism begin      # the B flank
            substrates: A, B; products: P, Q
            steps: begin
                E + A ⇌ EA
                EA + B <--> EAB
                EAB <--> EPQ
                EQ + P ⇌ EPQ
                E + Q ⇌ EQ
            end
        end),
        (@enzyme_mechanism begin      # the P flank
            substrates: A, B; products: P, Q
            steps: begin
                E + A ⇌ EA
                EA + B ⇌ EAB
                EAB <--> EPQ
                EQ + P <--> EPQ
                E + Q ⇌ EQ
            end
        end)
    ])
        @test !(absent in kids)
        @test _testhelper_fitted(absent) == 6
        @test _testhelper_identifiable_rank(absent) == 5
    end
end

@testset "Mechanism — a flank whose group is shared still flips" begin
    # The ordered seed with B's abortive binding at EQ in B's group. One RE segment
    # holds every form; the iso is steady state. EAB's binding shares its group
    # with EQ + B ⇌ EBQ, so the chain EAB → EPQ does not qualify and
    # no group is a flank. Each flip cuts the segment alone:
    #   A: {EA, EAB} | {E, EQ, EBQ, EPQ}             → 2 segments
    #   B: {EAB} | {EBQ} | {E, EA, EQ, EPQ}           → 3 segments
    #   P: {EPQ} | {E, EA, EAB, EQ, EBQ}              → 2 segments
    #   Q: {E, EA, EAB} | {EQ, EBQ, EPQ}              → 2 segments
    # Every segment keeps a bottom form, and each flipped group holds a step of the
    # cycle E → EA → EAB → EPQ → EQ → E, so it carries flux: four children,
    # each one constant and one rank above the seed's 5.
    seed = @enzyme_mechanism begin
        substrates: A, B; products: P, Q
        steps: begin
            E + A ⇌ EA
            (EA + B ⇌ EAB, EQ + B ⇌ EBQ)
            EAB <--> EPQ
            EQ + P ⇌ EPQ
            E + Q ⇌ EQ
        end
    end
    a_flip = @enzyme_mechanism begin
        substrates: A, B; products: P, Q
        steps: begin
            E + A <--> EA
            (EA + B ⇌ EAB, EQ + B ⇌ EBQ)
            EAB <--> EPQ
            EQ + P ⇌ EPQ
            E + Q ⇌ EQ
        end
    end
    b_flip = @enzyme_mechanism begin
        substrates: A, B; products: P, Q
        steps: begin
            E + A ⇌ EA
            (EA + B <--> EAB, EQ + B <--> EBQ)
            EAB <--> EPQ
            EQ + P ⇌ EPQ
            E + Q ⇌ EQ
        end
    end
    p_flip = @enzyme_mechanism begin
        substrates: A, B; products: P, Q
        steps: begin
            E + A ⇌ EA
            (EA + B ⇌ EAB, EQ + B ⇌ EBQ)
            EAB <--> EPQ
            EQ + P <--> EPQ
            E + Q ⇌ EQ
        end
    end
    q_flip = @enzyme_mechanism begin
        substrates: A, B; products: P, Q
        steps: begin
            E + A ⇌ EA
            (EA + B ⇌ EAB, EQ + B ⇌ EBQ)
            EAB <--> EPQ
            EQ + P ⇌ EPQ
            E + Q <--> EQ
        end
    end
    @test isempty(EnzymeRates.Enumeration._chain_flank_groups(seed))
    kids = _testhelper_move(EnzymeRates.Enumeration._expand_re_to_ss, seed)
    @test length(kids) == 4
    @test Set(kids) == Set([a_flip, b_flip, p_flip, q_flip])
    @test _testhelper_fitted(seed) == _testhelper_identifiable_rank(seed) == 5
    for k in kids
        @test _testhelper_fitted(k) == _testhelper_identifiable_rank(k) == 6
    end
end

@testset "Allosteric mechanism — a chain flank does not flip at multiplicity 2" begin
    # A variant of the ordered seed over two catalytic subunits: the iso and the A
    # binding :OnlyA, the other bindings :EqualAI. The B and P groups are the flanks of
    # the iso, so only the A and Q flips are children. The parent needs a hyperbolic
    # scheme (`_hyperbolic_catalysis`). With two segments an arborescence toward a root
    # holds one edge, which leaves the other segment, so the scheme fails only when a
    # metabolite scores 2 on one edge or one segment, or scores in one segment and on
    # an edge leaving the other. Each flip passes:
    #   A: {E, EQ, EPQ} | {EA, EAB}. A scores on E → EA only; B in the
    #      second segment and on EAB → EPQ, which leaves it; P and Q in the
    #      first and on EPQ → EAB, which leaves it.
    #   Q: the mirror of A.
    # Two children, tags kept. The parent is identifiable at 6; each child adds one
    # constant and gains one.
    seed = @enzyme_mechanism begin
        substrates: A, B; products: P, Q
        steps: begin
            E + A ⇌ EA
            EA + B ⇌ EAB
            EAB <--> EPQ
            EQ + P ⇌ EPQ
            E + Q ⇌ EQ
        end
    end
    rxn2 = @enzyme_reaction begin
        substrates: A[C], B[N]
        products: P[C], Q[N]
        oligomeric_state: 2
    end
    am = @enzyme_mechanism begin
        substrates: A, B ; products: P, Q ; oligomeric_state: 2
        steps: begin
            E + A ⇌ EA                :: OnlyA
            EA + B ⇌ EAB          :: EqualAI
            EAB <--> EPQ        :: OnlyA
            EQ + P ⇌ EPQ          :: EqualAI
            E + Q ⇌ EQ                :: EqualAI
        end
    end
    @test am in _testhelper_move(EnzymeRates.Enumeration._expand_to_allosteric, seed, rxn2)
    a_flip = @enzyme_mechanism begin
        substrates: A, B ; products: P, Q ; oligomeric_state: 2
        steps: begin
            E + A <--> EA             :: OnlyA
            EA + B ⇌ EAB          :: EqualAI
            EAB <--> EPQ        :: OnlyA
            EQ + P ⇌ EPQ          :: EqualAI
            E + Q ⇌ EQ                :: EqualAI
        end
    end
    q_flip = @enzyme_mechanism begin
        substrates: A, B ; products: P, Q ; oligomeric_state: 2
        steps: begin
            E + A ⇌ EA                :: OnlyA
            EA + B ⇌ EAB          :: EqualAI
            EAB <--> EPQ        :: OnlyA
            EQ + P ⇌ EPQ          :: EqualAI
            E + Q <--> EQ             :: EqualAI
        end
    end
    kids = _testhelper_move(EnzymeRates.Enumeration._expand_re_to_ss, am)
    @test length(kids) == 2
    @test Set(kids) == Set([a_flip, q_flip])
    # candidate_mechanisms returns every flip, and every child it returns at two subunits
    # is hyperbolic.
    returned = _testhelper_children(rxn2, [am])
    @test issubset(kids, returned)
    @test all(EnzymeRates.Enumeration._hyperbolic_catalysis,
              _testhelper_candidate.(returned))
    @test _testhelper_fitted(am) == _testhelper_identifiable_rank(am) == 6
    for k in kids
        @test _testhelper_fitted(k) == _testhelper_identifiable_rank(k) == 7
    end
end

@testset "an allosteric mechanism's chain flanks are not flipped" begin
    # The iso EAB <--> EPQ is alone in its group at steady state, with one binding on
    # each side: B into EAB and P into EPQ. Those two groups are its flanks, whatever
    # the conformation states of the groups. Only the A and Q bindings can flip.
    for n in 1:2
        am = n == 1 ? (@enzyme_mechanism begin
            substrates: A, B ; products: P, Q ; oligomeric_state: 1
            steps: begin
                E + A ⇌ EA            :: OnlyA
                EA + B ⇌ EAB          :: EqualAI
                EAB <--> EPQ          :: OnlyA
                EQ + P ⇌ EPQ          :: EqualAI
                E + Q ⇌ EQ            :: EqualAI
            end
        end) : (@enzyme_mechanism begin
            substrates: A, B ; products: P, Q ; oligomeric_state: 2
            steps: begin
                E + A ⇌ EA            :: OnlyA
                EA + B ⇌ EAB          :: EqualAI
                EAB <--> EPQ          :: OnlyA
                EQ + P ⇌ EPQ          :: EqualAI
                E + Q ⇌ EQ            :: EqualAI
            end
        end)
        a_flip = n == 1 ? (@enzyme_mechanism begin
            substrates: A, B ; products: P, Q ; oligomeric_state: 1
            steps: begin
                E + A <--> EA         :: OnlyA
                EA + B ⇌ EAB          :: EqualAI
                EAB <--> EPQ          :: OnlyA
                EQ + P ⇌ EPQ          :: EqualAI
                E + Q ⇌ EQ            :: EqualAI
            end
        end) : (@enzyme_mechanism begin
            substrates: A, B ; products: P, Q ; oligomeric_state: 2
            steps: begin
                E + A <--> EA         :: OnlyA
                EA + B ⇌ EAB          :: EqualAI
                EAB <--> EPQ          :: OnlyA
                EQ + P ⇌ EPQ          :: EqualAI
                E + Q ⇌ EQ            :: EqualAI
            end
        end)
        q_flip = n == 1 ? (@enzyme_mechanism begin
            substrates: A, B ; products: P, Q ; oligomeric_state: 1
            steps: begin
                E + A ⇌ EA            :: OnlyA
                EA + B ⇌ EAB          :: EqualAI
                EAB <--> EPQ          :: OnlyA
                EQ + P ⇌ EPQ          :: EqualAI
                E + Q <--> EQ         :: EqualAI
            end
        end) : (@enzyme_mechanism begin
            substrates: A, B ; products: P, Q ; oligomeric_state: 2
            steps: begin
                E + A ⇌ EA            :: OnlyA
                EA + B ⇌ EAB          :: EqualAI
                EAB <--> EPQ          :: OnlyA
                EQ + P ⇌ EPQ          :: EqualAI
                E + Q <--> EQ         :: EqualAI
            end
        end)
        members = _testhelper_group_members(am.catalytic)
        flanks = Set(EnzymeRates.Mechanisms.step_text(am.catalytic, only(members[g]))
                     for g in EnzymeRates.Enumeration._chain_flank_groups(am))
        @test flanks == Set(["EA + B ⇌ EAB", "EQ + P ⇌ EPQ"])
        kids = _testhelper_move(EnzymeRates.Enumeration._expand_re_to_ss, am)
        @test length(kids) == 2
        @test Set(kids) == Set([a_flip, q_flip])
    end
end

@testset "_chain_flank_groups: the flanks of steady-state chains only" begin
    # Ping-pong seed. EA → EstarP is steady state, alone in its group, and its two
    # ends each have one other step, a binding into that end, alone in its group: the
    # A and P groups are its flanks. EstarB ⇌ EQ is a chain of the same shape
    # whose isomerization is at rapid equilibrium, so the B and Q groups are no
    # flanks.
    m = @enzyme_mechanism begin
        substrates: A, B
        products: P, Q
        steps: begin
            E + A ⇌ EA
            Estar + B ⇌ EstarB
            E + Q ⇌ EQ
            Estar + P ⇌ EstarP
            EA <--> EstarP
            EstarB ⇌ EQ
        end
    end
    members = _testhelper_group_members(m.catalytic)
    flanks = Set(members[g] for g in EnzymeRates.Enumeration._chain_flank_groups(m))
    @test length(flanks) == 2
    @test all(ks -> length(ks) == 1 &&
              EnzymeRates.Enumeration._is_binding(m.catalytic.steps[only(ks)]),
              flanks)
    @test Set(m.catalytic.metabolites[only(m.catalytic.steps[only(ks)].taken_up)]
              for ks in flanks) == Set([:A, :P])
end

@testset "_chain_flank_groups: a flank binds into its end of the chain" begin
    # Iso uni-uni with P also binding ES abortively. Both isomerizations are steady
    # state and alone in their groups, and no chain qualifies:
    #   ES → EstarP: ES has two other steps, the S binding into it and the P
    #     binding out of it (ES + P ⇌ EPS), so it is no chain end.
    #   Estar → E: each end has exactly one other step, a binding alone in its group,
    #     but out of the end, not into it (Estar + P ⇌ EstarP, E + S ⇌ ES), so
    #     neither is a flank.
    # Without the orientation condition the second isomerization would name the S and P
    # groups, the flanks the first would name without the abortive step.
    m = @enzyme_mechanism begin
        substrates: S
        products: P
        steps: begin
            E + S ⇌ ES
            ES + P ⇌ EPS
            ES <--> EstarP
            Estar + P ⇌ EstarP
            Estar <--> E
        end
    end
    @test isempty(EnzymeRates.Enumeration._chain_flank_groups(m))
end

@testset "Mechanism — bi-bi ping-pong: 5 RE groups → 3 variants" begin
    # SEED: bi-bi ping-pong with Estar (residual) form. 5 singleton RE
    # groups + 1 SS iso group. EA → EstarP is a qualifying chain with a
    # steady-state isomerization: EA and EstarP each have one other step,
    # a binding into that form alone in its group, so its flanks, the A and P
    # groups, never flip (`_chain_flank_groups`). EstarB ⇌ EQ has the same
    # shape at rapid equilibrium, so the B and Q groups stay units. The B, Q
    # and RE iso groups each flip alone → 3.
    m = @enzyme_mechanism begin
        substrates: A, B
        products: P, Q
        steps: begin
            E + A ⇌ EA
            Estar + B ⇌ EstarB
            E + Q ⇌ EQ
            Estar + P ⇌ EstarP
            EA <--> EstarP
            EstarB ⇌ EQ
        end
    end
    _testhelper_assert_mechanism_invariants(m)

    result = _testhelper_move(EnzymeRates.Enumeration._expand_re_to_ss, m)

    # 1. count: 5 all-RE groups, of which A and P are chain flanks → 3 variants.
    @test length(result) == 3
    for r in result
        @test !_testhelper_is_allosteric(r)
        _testhelper_assert_mechanism_invariants(r)
    end

    # 2. property-style: each variant flips exactly one initial RE
    # group to all-SS. The flipped group is distinct across variants.
    flipped_groups = Int[]
    for r in result
        for (gi, (old_grp, new_grp)) in
                enumerate(zip(_testhelper_group_flags(m), _testhelper_group_flags(r)))
            all(old_grp) && !any(new_grp) && push!(flipped_groups, gi)
        end
    end
    @test length(unique(flipped_groups)) == 3

    # 3. preservation
    for r in result
        @test _testhelper_keeps_reaction(r, m)
    end

    # 4. The A and P flank flips, absent, keep the seed's rank (chain lemma).
    r0 = _testhelper_identifiable_rank(m)
    for absent in ([
        (@enzyme_mechanism begin      # the A flank
            substrates: A, B
            products: P, Q
            steps: begin
                E + A <--> EA
                Estar + B ⇌ EstarB
                E + Q ⇌ EQ
                Estar + P ⇌ EstarP
                EA <--> EstarP
                EstarB ⇌ EQ
            end
        end),
        (@enzyme_mechanism begin      # the P flank
            substrates: A, B
            products: P, Q
            steps: begin
                E + A ⇌ EA
                Estar + B ⇌ EstarB
                E + Q ⇌ EQ
                Estar + P <--> EstarP
                EA <--> EstarP
                EstarB ⇌ EQ
            end
        end)
    ])
        @test !(absent in result)
        @test _testhelper_identifiable_rank(absent) == r0
    end
end

@testset "Mechanism — ter-ter sequential" begin
    # SEED: ter-ter sequential ordered. 6 singleton RE binding groups
    # + 1 SS iso. The D and P groups are the flanks of the qualifying chain
    # EABD → EPQR (`_chain_flank_groups`) and never flip; the other
    # four each cut the one segment alone → 4 variants.
    m = @enzyme_mechanism begin
        substrates: A, B, D
        products: P, Q, R
        steps: begin
            E + A ⇌ EA
            EA + B ⇌ EAB
            EAB + D ⇌ EABD
            E + R ⇌ ER
            ER + Q ⇌ EQR
            EQR + P ⇌ EPQR
            EABD <--> EPQR
        end
    end
    _testhelper_assert_mechanism_invariants(m)

    result = _testhelper_move(EnzymeRates.Enumeration._expand_re_to_ss, m)
    @test length(result) == 4
    for r in result
        @test !_testhelper_is_allosteric(r)
        _testhelper_assert_mechanism_invariants(r)
    end
    # The D and P flank flips, absent, keep the seed's rank (chain lemma).
    r0 = _testhelper_identifiable_rank(m)
    for absent in ([
        (@enzyme_mechanism begin      # the D flank
            substrates: A, B, D
            products: P, Q, R
            steps: begin
                E + A ⇌ EA
                EA + B ⇌ EAB
                EAB + D <--> EABD
                E + R ⇌ ER
                ER + Q ⇌ EQR
                EQR + P ⇌ EPQR
                EABD <--> EPQR
            end
        end),
        (@enzyme_mechanism begin      # the P flank
            substrates: A, B, D
            products: P, Q, R
            steps: begin
                E + A ⇌ EA
                EA + B ⇌ EAB
                EAB + D ⇌ EABD
                E + R ⇌ ER
                ER + Q ⇌ EQR
                EQR + P <--> EPQR
                EABD <--> EPQR
            end
        end)
    ])
        @test !(absent in result)
        @test _testhelper_identifiable_rank(absent) == r0
    end
end

@testset "Allosteric mechanism — uni-uni has no flip: both bindings are chain flanks" begin
    # The binding of S into ES and the binding of P into EP are alone in their groups and
    # flank the steady-state isomerization ES <--> EP, whatever the tags of the three
    # groups, so no group is a unit.
    for am in [
        (@enzyme_mechanism begin
            substrates: S
            products: P
            oligomeric_state: 2
            steps: begin
                E + P ⇌ EP      :: EqualAI
                E + S ⇌ ES      :: EqualAI
                ES <--> EP    :: EqualAI
            end
        end),
        (@enzyme_mechanism begin
            substrates: S
            products: P
            oligomeric_state: 2
            steps: begin
                E + P ⇌ EP      :: OnlyA
                E + S ⇌ ES      :: OnlyA
                ES <--> EP    :: EqualAI
            end
        end),
        (@enzyme_mechanism begin
            substrates: S
            products: P
            oligomeric_state: 2
            steps: begin
                E + P ⇌ EP      :: EqualAI
                E + S ⇌ ES      :: NonequalAI
                ES <--> EP    :: EqualAI
            end
        end),
        (@enzyme_mechanism begin
            substrates: S
            products: P
            oligomeric_state: 2
            steps: begin
                E + P ⇌ EP      :: NonequalAI
                E + S ⇌ ES      :: NonequalAI
                ES <--> EP    :: NonequalAI
            end
        end),
    ]
        _testhelper_assert_mechanism_invariants(am)
        @test length(EnzymeRates.Enumeration._chain_flank_groups(am)) == 2
        @test isempty(_testhelper_move(EnzymeRates.Enumeration._expand_re_to_ss, am))
    end
end

@testset "Mechanism — a substrate's modifier copy never flips; catalytic groups do" begin
    # Ordered bi-bi with A also bound as a competitive-inhibitor copy at E and
    # EQ, the Q binding mirrored onto the copy forms in the Q group, as the
    # dead-end move builds it. The copy group binds a regulator and is never a
    # unit. A, B, P and Q each divide the one segment alone (the Q flip takes its
    # mirror with it and isolates {EQ, EPQ, EAmodQ}), and each flipped
    # group shares a block with the chemistry edge, so all four carry flux. The copy
    # binds at E and EQ, so EAB and EPQ keep one step each beside the
    # isomerization, a binding into that form alone in its group: B and P are the
    # flanks of a qualifying chain and never flip (`_chain_flank_groups`). Their
    # flips, absent, keep the parent's rank. Children: the A and Q flips.
    m = @enzyme_mechanism begin
        substrates: A, B
        products: P, Q
        steps: begin
            E + A ⇌ EA
            EA + B ⇌ EAB
            EAB <--> EPQ
            EQ + P ⇌ EPQ
            (E + Q ⇌ EQ, EAmod + Q ⇌ EAmodQ)
            (E + A ⇌ EAmod, EQ + A ⇌ EAmodQ)
        end
    end
    expected = [
        (@enzyme_mechanism begin
            substrates: A, B
            products: P, Q
            steps: begin
                E + A <--> EA
                EA + B ⇌ EAB
                EAB <--> EPQ
                EQ + P ⇌ EPQ
                (E + Q ⇌ EQ, EAmod + Q ⇌ EAmodQ)
                (E + A ⇌ EAmod, EQ + A ⇌ EAmodQ)
            end
        end),
        (@enzyme_mechanism begin
            substrates: A, B
            products: P, Q
            steps: begin
                E + A ⇌ EA
                EA + B ⇌ EAB
                EAB <--> EPQ
                EQ + P ⇌ EPQ
                (E + Q <--> EQ, EAmod + Q <--> EAmodQ)
                (E + A ⇌ EAmod, EQ + A ⇌ EAmodQ)
            end
        end),
    ]
    kids = _testhelper_move(EnzymeRates.Enumeration._expand_re_to_ss, m)
    @test length(kids) == 2
    @test Set(kids) == Set(expected)
    r0 = _testhelper_identifiable_rank(m)
    for absent in [
        (@enzyme_mechanism begin      # the B flank
            substrates: A, B
            products: P, Q
            steps: begin
                E + A ⇌ EA
                EA + B <--> EAB
                EAB <--> EPQ
                EQ + P ⇌ EPQ
                (E + Q ⇌ EQ, EAmod + Q ⇌ EAmodQ)
                (E + A ⇌ EAmod, EQ + A ⇌ EAmodQ)
            end
        end),
        @enzyme_mechanism begin      # the P flank
            substrates: A, B
            products: P, Q
            steps: begin
                E + A ⇌ EA
                EA + B ⇌ EAB
                EAB <--> EPQ
                EQ + P <--> EPQ
                (E + Q ⇌ EQ, EAmod + Q ⇌ EAmodQ)
                (E + A ⇌ EAmod, EQ + A ⇌ EAmodQ)
            end
        end
    ]
        @test !(absent in kids)
        @test _testhelper_identifiable_rank(absent) == r0
    end
    for r in kids
        _testhelper_assert_mechanism_invariants(r)
    end
    for r in kids
        c = _testhelper_candidate(r)
        for k in eachindex(r.catalytic.steps)
            _testhelper_binds_copy(c, k) && @test r.catalytic.steps[k].rapid_equilibrium
        end
    end
end

@testset "Allosteric mechanism — the modifier copy never flips; tags are preserved" begin
    # The same mechanism as an allosteric one with a dead inactive conformation
    # (the chemistry `:OnlyA`, every binding `:EqualAI`) and one catalytic subunit.
    # The copy binds at E and EQ, so EAB and EPQ keep one step each beside the
    # isomerization: B and P are the flanks of a chain and never flip
    # (`_chain_flank_groups`). The A and Q groups flip. Every child keeps one `:OnlyA`
    # group, the chemistry, and the multiplicity and (empty) regulatory sites.
    am = @enzyme_mechanism begin
        substrates: A, B
        products: P, Q
        oligomeric_state: 1
        steps: begin
            E + A ⇌ EA                                            :: EqualAI
            EA + B ⇌ EAB                                      :: EqualAI
            EAB <--> EPQ                                    :: OnlyA
            EQ + P ⇌ EPQ                                      :: EqualAI
            (E + Q ⇌ EQ, EAmod + Q ⇌ EAmodQ)            :: EqualAI
            (E + A ⇌ EAmod, EQ + A ⇌ EAmodQ)  :: EqualAI
        end
    end
    expected = [
        (@enzyme_mechanism begin
            substrates: A, B
            products: P, Q
            oligomeric_state: 1
            steps: begin
                E + A <--> EA                                         :: EqualAI
                EA + B ⇌ EAB                                      :: EqualAI
                EAB <--> EPQ                                    :: OnlyA
                EQ + P ⇌ EPQ                                      :: EqualAI
                (E + Q ⇌ EQ, EAmod + Q ⇌ EAmodQ)            :: EqualAI
                (E + A ⇌ EAmod, EQ + A ⇌ EAmodQ)  :: EqualAI
            end
        end),
        (@enzyme_mechanism begin
            substrates: A, B
            products: P, Q
            oligomeric_state: 1
            steps: begin
                E + A ⇌ EA                                            :: EqualAI
                EA + B ⇌ EAB                                      :: EqualAI
                EAB <--> EPQ                                    :: OnlyA
                EQ + P ⇌ EPQ                                      :: EqualAI
                (E + Q <--> EQ, EAmod + Q <--> EAmodQ)      :: EqualAI
                (E + A ⇌ EAmod, EQ + A ⇌ EAmodQ)  :: EqualAI
            end
        end),
    ]
    kids = _testhelper_move(EnzymeRates.Enumeration._expand_re_to_ss, am)
    @test length(kids) == 2
    @test Set(kids) == Set(expected)
    for r in kids
        _testhelper_assert_mechanism_invariants(r)
    end
    for r in kids
        @test count(==(:OnlyA), r.allosteric.allosteric_states) == 1
        onlya = findfirst(==(:OnlyA), r.allosteric.allosteric_states)
        cm = r.catalytic
        @test EnzymeRates.Mechanisms._is_iso(
            cm.steps[findfirst(e -> e.kinetic_group == onlya, cm.steps)])
        @test r.allosteric.oligomeric_state == am.allosteric.oligomeric_state
        @test r.allosteric.regulatory_sites == am.allosteric.regulatory_sites
        @test _testhelper_keeps_reaction(r, am)
    end
end

@testset "_expand_re_to_ss flips inhibitor-bound mirrors together" begin
    # A catalytic binding and its inhibitor-bound mirror (the same binding
    # with every inhibitor stripped) sit in separate all-RE kinetic groups
    # on a closed inhibitor square: A binding to E, and A binding to the
    # inhibitor-bound form EI, with I binding both E and EA. Neither
    # cuts a rapid-equilibrium segment alone, because the other supplies an
    # RE route around the square, so the minimal-set search emits them only
    # together: they flip to SS as a pair, never one without the other.
    # Built via the macro; a separate catalytic inhibitor I keeps the
    # substrate/inhibitor roles unambiguous. Binding is ordered (A then B)
    # so the pair flip stays hyperbolic; a conformational parent drops a
    # child whose catalytic scheme carries a concentration power, and the
    # random-order pair flip would carry B². The B and P bindings are the flanks of the
    # isomerization EAB <--> EP (`_chain_flank_groups`) and never flip, so the one
    # child flips the A binding with its mirror.
    m = EnzymeRates.@enzyme_mechanism begin
        substrates: A, B
        products: P
        modifiers: I
        steps: begin
            E + A ⇌ EA          :: EqualAI
            EA + B ⇌ EAB    :: EqualAI
            EAB <--> EP      :: EqualAI
            E + P ⇌ EP          :: EqualAI
            E + I ⇌ EI          :: EqualAI
            EA + I ⇌ EAI    :: EqualAI
            EI + A ⇌ EAI    :: EqualAI
        end
    end
    kids = _testhelper_move(EnzymeRates.Enumeration._expand_re_to_ss, m)
    expected = [
        @enzyme_mechanism(begin
            substrates: A, B
            products: P
            modifiers: I
            steps: begin
                E + A <--> EA        :: EqualAI
                EA + B ⇌ EAB    :: EqualAI
                EAB <--> EP      :: EqualAI
                E + P ⇌ EP          :: EqualAI
                E + I ⇌ EI          :: EqualAI
                EA + I ⇌ EAI    :: EqualAI
                EI + A <--> EAI  :: EqualAI
            end
        end)]
    @test length(kids) == 1
    @test Set(kids) == Set(expected)
    # Mirror-lock: in every child the A binding at E and its mirror at EI carry one
    # flag.
    _testhelper_flag(c, from, to) = only(e.rapid_equilibrium for e in c.catalytic.steps
        if (String(c.catalytic.enzyme_forms[e.from]),
            String(c.catalytic.enzyme_forms[e.to])) == (from, to))
    @test all(c -> _testhelper_flag(c, "E", "EA") == _testhelper_flag(c, "EI", "EAI"), kids)
end

@testset "Mechanism — all-SS seed: empty (negative)" begin
    # If every catalytic group is already SS, the move has no RE
    # group to fire on → empty.
    m_all_ss = @enzyme_mechanism begin
        substrates: S
        products: P
        steps: begin
            E + P <--> EP
            E + S <--> ES
            ES <--> EP
        end
    end
    @test isempty(_testhelper_move(EnzymeRates.Enumeration._expand_re_to_ss, m_all_ss))
end

end

# ─── _expand_split_kinetic_group ───────────────────────────────────────
@testset "_expand_split_kinetic_group" begin

@testset "Mechanism — mixed RE/SS multi-step groups: every child gains" begin
    # SEED: bi-bi random where the A-binding kinetic group is SS
    # (size-2, both SS) and the B-binding kinetic group is RE
    # (size-2, both RE). The remaining P-binding (size-2 RE),
    # Q-binding (size-2 RE), and iso (singleton SS) groups are
    # unchanged. Total: 9 steps, 5 kinetic groups. Each multi-step
    # group divides by binding context, and a division is emitted only
    # when it raises the independent-parameter count.
    m_seed = @enzyme_mechanism begin
        substrates: A, B
        products: P, Q
        steps: begin
            (E + A <--> EA, EB + A <--> EAB)
            (E + B ⇌ EB, EA + B ⇌ EAB)
            (E + P ⇌ EP, EQ + P ⇌ EPQ)
            (E + Q ⇌ EQ, EP + Q ⇌ EPQ)
            EAB <--> EPQ
        end
    end
    _testhelper_assert_mechanism_invariants(m_seed)

    result = _testhelper_move(EnzymeRates.Enumeration._expand_split_kinetic_group, m_seed)
    @test !isempty(result)
    base = EnzymeRates.Derivation.fitted_parameter_count(m_seed)
    for r in result
        @test EnzymeRates.Derivation.fitted_parameter_count(r) > base
        @test !_testhelper_is_allosteric(r)
        _testhelper_assert_mechanism_invariants(r)
        @test length(r.catalytic.steps) == length(m_seed.catalytic.steps)
        @test _testhelper_keeps_reaction(r, m_seed)
    end
    # The SS A-binding group has no equilibrium constant to tie, so its
    # context split by B gains alone and is emitted as a single split.
    _testhelper_lone_ss_a(r, ks) =
        length(ks) == 1 && !r.catalytic.steps[only(ks)].rapid_equilibrium &&
        r.catalytic.metabolites[only(r.catalytic.steps[only(ks)].taken_up)] == :A
    a_split =
        [r for r in result
         if length(_testhelper_group_members(r.catalytic)) ==
                length(_testhelper_group_members(m_seed.catalytic)) + 1 &&
             any(ks -> _testhelper_lone_ss_a(r, ks),
                 _testhelper_group_members(r.catalytic))]
    @test length(a_split) == 1
end

@testset "Allosteric mechanism — SS multi-step :NonequalAI split" begin
    # SEED: bi-bi allosteric where one multi-step group is BOTH SS AND
    # :NonequalAI. Splitting this group costs more parameters than
    # splitting a :EqualAI RE group (factor 2 for SS pair × factor 2
    # for R/T-state pair).
    am = @enzyme_mechanism begin
        substrates: A, B
        products: P, Q
        oligomeric_state: 2
        steps: begin
            (E + A <--> EA, EB + A <--> EAB)    :: NonequalAI
            (E + B ⇌ EB, EA + B ⇌ EAB)          :: EqualAI
            (E + P ⇌ EP, EQ + P ⇌ EPQ)          :: EqualAI
            (E + Q ⇌ EQ, EP + Q ⇌ EPQ)          :: EqualAI
            EAB <--> EPQ                        :: EqualAI
        end
    end
    _testhelper_assert_mechanism_invariants(am)

    result = _testhelper_move(EnzymeRates.Enumeration._expand_split_kinetic_group, am)
    # The A group divides by B alone, a steady-state group with no constant to tie;
    # the P and Q groups divide only together, each by the other product. Each part
    # keeps its group's tag.
    expected = [
        @enzyme_mechanism(begin
            substrates: A, B
            products: P, Q
            oligomeric_state: 2
            steps: begin
                E + A <--> EA                     :: NonequalAI
                EB + A <--> EAB               :: NonequalAI
                (E + B ⇌ EB, EA + B ⇌ EAB)  :: EqualAI
                (E + P ⇌ EP, EQ + P ⇌ EPQ)  :: EqualAI
                (E + Q ⇌ EQ, EP + Q ⇌ EPQ)  :: EqualAI
                EAB <--> EPQ                :: EqualAI
            end
        end),
        @enzyme_mechanism(begin
            substrates: A, B
            products: P, Q
            oligomeric_state: 2
            steps: begin
                (E + A <--> EA, EB + A <--> EAB)  :: NonequalAI
                (E + B ⇌ EB, EA + B ⇌ EAB)        :: EqualAI
                E + P ⇌ EP                              :: EqualAI
                EQ + P ⇌ EPQ                        :: EqualAI
                E + Q ⇌ EQ                              :: EqualAI
                EP + Q ⇌ EPQ                        :: EqualAI
                EAB <--> EPQ                      :: EqualAI
            end
        end)]
    @test length(result) == 2
    @test Set(result) == Set(expected)

    # 1. every emitted child raises the independent-parameter count.
    @test !isempty(result)
    @test all(r -> EnzymeRates.Derivation.fitted_parameter_count(r) >
                   EnzymeRates.Derivation.fitted_parameter_count(am), result)

    # 2. the :NonequalAI A-binding group is SS, so it has no equilibrium
    # constant to tie and its context split gains on its own: a child
    # that splits that group and nothing else is emitted.
    cm = am.catalytic
    a_g = cm.steps[_testhelper_edge(cm, "E + A <--> EA")].kinetic_group
    a_group = findall(e -> e.kinetic_group == a_g, cm.steps)
    @test any(r -> length(_testhelper_group_members(r.catalytic)) ==
              length(_testhelper_group_members(cm)) + 1 &&
              _testhelper_ungrouped(r.catalytic) == _testhelper_ungrouped(cm) &&
              length(unique(r.catalytic.steps[k].kinetic_group for k in a_group)) == 2,
              result)

    # 3. kind and invariants
    for r in result
        @test _testhelper_is_allosteric(r)
        _testhelper_assert_mechanism_invariants(r)
    end

    # 4. tag inheritance: each result group inherits the tag of the parent
    # group whose steps contain it. A child subdivides one or more groups;
    # both halves of each carry that group's tag, and every other group
    # keeps its own. Group ORDER is canonical (not source-preserved), so match each
    # result group to its parent by step content rather than by position.
    for r in result
        @test length(r.allosteric.allosteric_states) ==
              length(_testhelper_group_members(r.catalytic))
        @test _testhelper_ungrouped(r.catalytic) == _testhelper_ungrouped(cm)
        for (k, g) in enumerate(e.kinetic_group for e in r.catalytic.steps)
            @test r.allosteric.allosteric_states[g] ==
                  am.allosteric.allosteric_states[cm.steps[k].kinetic_group]
        end
    end

    # 5. preservation
    for r in result
        @test r.allosteric.oligomeric_state == am.allosteric.oligomeric_state
        @test r.allosteric.regulatory_sites == am.allosteric.regulatory_sites
        @test _testhelper_keeps_reaction(r, am)
    end
end

@testset "Mechanism — all singleton groups: empty (negative)" begin
    # If every group is a singleton, no split is possible.
    m_seed = @enzyme_mechanism begin
        substrates: S
        products: P
        steps: begin
            E + P ⇌ EP
            E + S ⇌ ES
            ES <--> EP
        end
    end
    @test isempty(_testhelper_move(EnzymeRates.Enumeration._expand_split_kinetic_group,
                                   m_seed))
end

@testset "Allosteric mechanism — bi-bi RE groups: every emitted split gains" begin
    # SEED: bi-bi allosteric with mixed tags (:NonequalAI, :EqualAI),
    # but both multi-step groups (A, B) are RE bindings, each closing
    # its own per-conformer thermodynamic cycle, so a division of one
    # group alone may be tied back by that cycle. Tag inheritance is
    # covered by "Allosteric mechanism — SS multi-step :NonequalAI
    # split" above.
    m_seed = @enzyme_mechanism begin
        substrates: A, B
        products: P, Q
        oligomeric_state: 2
        steps: begin
            (E + A ⇌ EA, EB + A ⇌ EAB)    :: NonequalAI
            (E + B ⇌ EB, EA + B ⇌ EAB)    :: EqualAI
            E + P ⇌ EP             :: EqualAI
            EP + Q ⇌ EPQ       :: EqualAI
            E + Q ⇌ EQ             :: EqualAI
            EQ + P ⇌ EPQ       :: EqualAI
            EAB <--> EPQ     :: EqualAI
        end
    end
    result = _testhelper_move(EnzymeRates.Enumeration._expand_split_kinetic_group, m_seed)
    @test !isempty(result)
    @test all(r -> EnzymeRates.Derivation.fitted_parameter_count(r) >
                   EnzymeRates.Derivation.fitted_parameter_count(m_seed), result)
end

@testset "Mechanism — a split part without flux is emitted at rapid equilibrium" begin
    # Ordered bi-bi whose steady-state B group holds the catalytic binding
    # EA + B and the abortive binding EQ + B. Splitting by context leaves
    # EQ + B → EBQ alone: EBQ is a dead end, so the step is a pendant
    # bridge of the segment graph and carries no flux. Its two constants would
    # enter the rate only as their ratio, so the part is emitted at rapid
    # equilibrium, one constant fewer than the raw split and the same family.
    # EBQ lies on no cycle, so the new constant raises the independent count: one
    # child.
    m = @enzyme_mechanism begin
        substrates: A, B
        products: P, Q
        steps: begin
            E + A ⇌ EA
            (EA + B <--> EAB, EQ + B <--> EBQ)
            EAB <--> EPQ
            EQ + P ⇌ EPQ
            E + Q ⇌ EQ
        end
    end
    reverted = @enzyme_mechanism begin
        substrates: A, B
        products: P, Q
        steps: begin
            E + A ⇌ EA
            EA + B <--> EAB
            EQ + B ⇌ EBQ
            EAB <--> EPQ
            EQ + P ⇌ EPQ
            E + Q ⇌ EQ
        end
    end
    raw = @enzyme_mechanism begin
        substrates: A, B
        products: P, Q
        steps: begin
            E + A ⇌ EA
            EA + B <--> EAB
            EQ + B <--> EBQ
            EAB <--> EPQ
            EQ + P ⇌ EPQ
            E + Q ⇌ EQ
        end
    end
    kids = _testhelper_move(EnzymeRates.Enumeration._expand_split_kinetic_group, m)
    @test length(kids) == 1
    @test Set(kids) == Set([reverted])
    @test !(raw in kids)
    @test EnzymeRates.Derivation.fitted_parameter_count(reverted) ==
        EnzymeRates.Derivation.fitted_parameter_count(m) + 1
    r_rev = _testhelper_identifiable_rank(reverted)
    r_raw = _testhelper_identifiable_rank(raw)
    @test r_rev == r_raw                                          # the same family
    @test _testhelper_fitted(reverted) == _testhelper_fitted(raw) - 1  # one constant fewer
    # Both children keep one phantom of another class: once the abortive step no
    # longer pins kf_B/kr_B, EA + B → EAB is a steady-state binding into a
    # form with one exit, and its three constants enter the law through two
    # combinations. The revert removes exactly the zero-flux phantom.
    @test _testhelper_fitted(raw) - r_raw == _testhelper_fitted(reverted) - r_rev + 1
end

@testset "Mechanism — a split leaving a redundant copy part is not emitted" begin
    # Ordered bi-bi with B also bound as a competitive-inhibitor copy at EA and
    # EQ. EABmod has the composition of EAB; EBmodQ is new, so
    # the placement stands. Splitting the copy group by context would leave
    # EA + B alone, a group that binds only where it duplicates EAB, and
    # EAB is formed by the B group alone and left by the chemistry alone, so the
    # gauge exists: the child would be its parent plus a phantom. Every other group
    # holds one step, so the move emits nothing.
    m = @enzyme_mechanism begin
        substrates: A, B
        products: P, Q
        steps: begin
            E + A ⇌ EA
            EA + B ⇌ EAB
            EAB <--> EPQ
            EQ + P ⇌ EPQ
            E + Q ⇌ EQ
            (EA + B ⇌ EABmod, EQ + B ⇌ EBmodQ)
        end
    end
    @test isempty(_testhelper_move(EnzymeRates.Enumeration._expand_split_kinetic_group, m))
    absent = @enzyme_mechanism begin
        substrates: A, B
        products: P, Q
        steps: begin
            E + A ⇌ EA
            EA + B ⇌ EAB
            EAB <--> EPQ
            EQ + P ⇌ EPQ
            E + Q ⇌ EQ
            EA + B ⇌ EABmod
            EQ + B ⇌ EBmodQ
        end
    end
    @test EnzymeRates.Derivation.fitted_parameter_count(absent) ==
        EnzymeRates.Derivation.fitted_parameter_count(m) + 1
    @test _testhelper_identifiable_rank(absent) == _testhelper_identifiable_rank(m)
end

@testset "Mechanism — a split that makes a copy group redundant is not emitted" begin
    # The shared A group with abortive EAQ and A as its own inhibitor at E: the
    # copy duplicates EA, but the A group also binds at EQ, so the gauge fails and
    # the mechanism is valid. Splitting the A group alone or the Q group alone gains
    # nothing (the cycle E, EA, EAQ, EQ ties the new constant back); splitting
    # both frees one. In that child the A binding at E sits alone in its group and the
    # Q binding at EA in its own: every group holds one ratio, the gauge exists, and
    # the copy is redundant (7 fitted, rank 6). Its family is that of the split parent
    # without the copy, so the move emits nothing.
    m = @enzyme_mechanism begin
        substrates: A, B
        products: P, Q
        steps: begin
            (E + A ⇌ EA, EQ + A ⇌ EAQ)
            EA + B ⇌ EAB
            EAB <--> EPQ
            EQ + P ⇌ EPQ
            (E + Q ⇌ EQ, EA + Q ⇌ EAQ)
            E + A ⇌ EAmod
        end
    end
    absent = @enzyme_mechanism begin
        substrates: A, B
        products: P, Q
        steps: begin
            E + A ⇌ EA
            EQ + A ⇌ EAQ
            EA + B ⇌ EAB
            EAB <--> EPQ
            EQ + P ⇌ EPQ
            E + Q ⇌ EQ
            EA + Q ⇌ EAQ
            E + A ⇌ EAmod
        end
    end
    @test EnzymeRates.Enumeration._assert_emission_rules(
        _testhelper_candidate(m)) === nothing
    @test isempty(_testhelper_move(EnzymeRates.Enumeration._expand_split_kinetic_group, m))
    @test !isempty(
        EnzymeRates.Enumeration._redundant_copy_groups(_testhelper_candidate(absent)))
    @test EnzymeRates.Derivation.fitted_parameter_count(absent) ==
        EnzymeRates.Derivation.fitted_parameter_count(m) + 1
    @test _testhelper_fitted(absent) == 7 && _testhelper_identifiable_rank(absent) == 6
end

@testset "Mechanism — two copies split down to complexes of one composition" begin
    # Ordered bi-bi with Q bound as a competitive-inhibitor copy at EA and at EQ in two
    # groups, and A as a copy at EA and EQ in one group. EAmodQ has the
    # composition {A, Q} of EAQmod, a copy's complex, and no productive form has
    # it, so neither part of the A copy's context bipartition is twin-only; every other
    # group holds one step, and the move emits the pair. The child carries one phantom:
    # with EQ + A alone in its group, its complex and EAQmod enter the law
    # through one coefficient, K_A·K_Q* + K_Q·K_A*, so 9 fitted constants have rank 8.
    # The rule tolerates it: rejecting it would also reject identifiable splits whose
    # other copy is pinned.
    m = @enzyme_mechanism begin
        substrates: A, B
        products: P, Q
        steps: begin
            E + A ⇌ EA
            EA + B ⇌ EAB
            EAB <--> EPQ
            EQ + P ⇌ EPQ
            E + Q ⇌ EQ
            EA + Q ⇌ EAQmod
            EQ + Q ⇌ EQQmod
            (EA + A ⇌ EAAmod, EQ + A ⇌ EAmodQ)
        end
    end
    emitted = @enzyme_mechanism begin
        substrates: A, B
        products: P, Q
        steps: begin
            E + A ⇌ EA
            EA + B ⇌ EAB
            EAB <--> EPQ
            EQ + P ⇌ EPQ
            E + Q ⇌ EQ
            EA + Q ⇌ EAQmod
            EQ + Q ⇌ EQQmod
            EA + A ⇌ EAAmod
            EQ + A ⇌ EAmodQ
        end
    end
    kids = _testhelper_move(EnzymeRates.Enumeration._expand_split_kinetic_group, m)
    @test length(kids) == 1
    @test Set(kids) == Set([emitted])
    @test _testhelper_fitted(emitted) == _testhelper_identifiable_rank(emitted) + 1
    @test _testhelper_fitted(m) == _testhelper_identifiable_rank(m)
end

@testset "Allosteric mechanism — the split count of a law over the term cap" begin
    # Merged random-order ter-bi with the abortive complexes EBP and EABP. Every
    # steady-state binding sits alone in its group except the A bindings at EB and
    # EBP, and the active state has more denominator terms than
    # MAX_RATE_EQUATION_TERMS. Splitting the P group by A leaves EAB + P alone, and the
    # square EB, EAB, EABP, EBP, closed by the shared A group, ties its K
    # to the other part's: the joint solve counts that child, which derives no rate
    # polynomial, at the parent's count, so the split gains nothing and is not emitted.
    # The move emits the two splits that gain, the P group by B and the A group by P, and
    # candidate_mechanisms keeps them beside the other moves' children.
    am = @enzyme_mechanism begin
        substrates: A, B, D
        products: P, Q
        steps: begin
            E + A <--> EA                                         :: OnlyA
            (EB + A <--> EAB, EBP + A <--> EABP)    :: OnlyA
            ED + A <--> EAD                                   :: OnlyA
            EBD + A <--> EPQ                                :: OnlyA
            E + B <--> EB                                         :: OnlyA
            EA + B <--> EAB                                   :: OnlyA
            ED + B <--> EBD                                   :: OnlyA
            EAD + B <--> EPQ                                :: OnlyA
            E + D <--> ED                                         :: OnlyA
            EA + D <--> EAD                                   :: OnlyA
            EB + D <--> EBD                                   :: OnlyA
            EAB + D <--> EPQ                                :: OnlyA
            (EAB + P ⇌ EABP, EB + P ⇌ EBP,
             EQ + P ⇌ EPQ)                                    :: OnlyA
            E + Q <--> EQ                                         :: EqualAI
        end
    end
    p_by_a = @enzyme_mechanism begin
        substrates: A, B, D
        products: P, Q
        steps: begin
            E + A <--> EA                                         :: OnlyA
            (EB + A <--> EAB, EBP + A <--> EABP)    :: OnlyA
            ED + A <--> EAD                                   :: OnlyA
            EBD + A <--> EPQ                                :: OnlyA
            E + B <--> EB                                         :: OnlyA
            EA + B <--> EAB                                   :: OnlyA
            ED + B <--> EBD                                   :: OnlyA
            EAD + B <--> EPQ                                :: OnlyA
            E + D <--> ED                                         :: OnlyA
            EA + D <--> EAD                                   :: OnlyA
            EB + D <--> EBD                                   :: OnlyA
            EAB + D <--> EPQ                                :: OnlyA
            EAB + P ⇌ EABP                                :: OnlyA
            (EB + P ⇌ EBP, EQ + P ⇌ EPQ)                :: OnlyA
            E + Q <--> EQ                                         :: EqualAI
        end
    end
    p_by_b = @enzyme_mechanism begin
        substrates: A, B, D
        products: P, Q
        steps: begin
            E + A <--> EA                                         :: OnlyA
            (EB + A <--> EAB, EBP + A <--> EABP)    :: OnlyA
            ED + A <--> EAD                                   :: OnlyA
            EBD + A <--> EPQ                                :: OnlyA
            E + B <--> EB                                         :: OnlyA
            EA + B <--> EAB                                   :: OnlyA
            ED + B <--> EBD                                   :: OnlyA
            EAD + B <--> EPQ                                :: OnlyA
            E + D <--> ED                                         :: OnlyA
            EA + D <--> EAD                                   :: OnlyA
            EB + D <--> EBD                                   :: OnlyA
            EAB + D <--> EPQ                                :: OnlyA
            (EAB + P ⇌ EABP, EB + P ⇌ EBP)          :: OnlyA
            EQ + P ⇌ EPQ                                      :: OnlyA
            E + Q <--> EQ                                         :: EqualAI
        end
    end
    a_by_p = @enzyme_mechanism begin
        substrates: A, B, D
        products: P, Q
        steps: begin
            E + A <--> EA                                         :: OnlyA
            EB + A <--> EAB                                   :: OnlyA
            EBP + A <--> EABP                             :: OnlyA
            ED + A <--> EAD                                   :: OnlyA
            EBD + A <--> EPQ                                :: OnlyA
            E + B <--> EB                                         :: OnlyA
            EA + B <--> EAB                                   :: OnlyA
            ED + B <--> EBD                                   :: OnlyA
            EAD + B <--> EPQ                                :: OnlyA
            E + D <--> ED                                         :: OnlyA
            EA + D <--> EAD                                   :: OnlyA
            EB + D <--> EBD                                   :: OnlyA
            EAB + D <--> EPQ                                :: OnlyA
            (EAB + P ⇌ EABP, EB + P ⇌ EBP,
             EQ + P ⇌ EPQ)                                    :: OnlyA
            E + Q <--> EQ                                         :: EqualAI
        end
    end
    @test EnzymeRates.Derivation.rate_equation_term_count(am) >
          EnzymeRates.Derivation.MAX_RATE_EQUATION_TERMS
    # Both conformations count: each inactive conformation keeps only E + Q <--> EQ, two
    # segments joined by one tree, so the sum is the active law's count plus 2.
    for x in (am, p_by_a, p_by_b, a_by_p)
        runs = EnzymeRates.Mechanisms.steps_running_in_allosteric_conformations(x)
        @test _testhelper_laplacian_tree_count(x.catalytic, runs.inactive) == 2
        @test EnzymeRates.Derivation.rate_equation_term_count(x) ==
              _testhelper_laplacian_tree_count(x.catalytic, runs.active) + 2
    end
    @test EnzymeRates.Derivation.fitted_parameter_count(p_by_a) ==
          EnzymeRates.Derivation.fitted_parameter_count(am)
    # The count derives no law; deriving the law itself still refuses one this large.
    @test_throws(EnzymeRates.Derivation.RateEquationTooLargeError,
                 rate_equation_string(p_by_a))
    @test_throws "polynomial terms" rate_equation_string(p_by_a)
    kids = _testhelper_move(EnzymeRates.Enumeration._expand_split_kinetic_group, am)
    @test length(kids) == 2
    @test Set(kids) == Set([p_by_b, a_by_p])
    others = [_testhelper_move(EnzymeRates.Enumeration._expand_re_to_ss, am);
              _testhelper_move(EnzymeRates.Enumeration._expand_change_allo_state, am)]
    @test Set(_testhelper_children(ter_bi_rxn, [am])) ==
          Set([others; p_by_b; a_by_p])
end

@testset "Allosteric mechanism — a parent over the term cap splits by its count" begin
    # Merged random-order ter-bi over MAX_RATE_EQUATION_TERMS whose P group is split by A:
    # EAB + P sits alone, its K tied to the other part's by the square EB, EAB,
    # EABP, EBP that the shared A group closes. The joint solve counts the parent
    # without deriving its law, so the move weighs every split against that count: the A
    # group split by P and the remaining P group split by B gain, and candidate_mechanisms
    # keeps them beside the other moves' children.
    am = @enzyme_mechanism begin
        substrates: A, B, D
        products: P, Q
        steps: begin
            E + A <--> EA                                         :: OnlyA
            (EB + A <--> EAB, EBP + A <--> EABP)    :: OnlyA
            ED + A <--> EAD                                   :: OnlyA
            EBD + A <--> EPQ                                :: OnlyA
            E + B <--> EB                                         :: OnlyA
            EA + B <--> EAB                                   :: OnlyA
            ED + B <--> EBD                                   :: OnlyA
            EAD + B <--> EPQ                                :: OnlyA
            E + D <--> ED                                         :: OnlyA
            EA + D <--> EAD                                   :: OnlyA
            EB + D <--> EBD                                   :: OnlyA
            EAB + D <--> EPQ                                :: OnlyA
            EAB + P ⇌ EABP                                :: OnlyA
            (EB + P ⇌ EBP, EQ + P ⇌ EPQ)                :: OnlyA
            E + Q <--> EQ                                         :: EqualAI
        end
    end
    a_by_p = @enzyme_mechanism begin
        substrates: A, B, D
        products: P, Q
        steps: begin
            E + A <--> EA                                         :: OnlyA
            EB + A <--> EAB                                   :: OnlyA
            EBP + A <--> EABP                             :: OnlyA
            ED + A <--> EAD                                   :: OnlyA
            EBD + A <--> EPQ                                :: OnlyA
            E + B <--> EB                                         :: OnlyA
            EA + B <--> EAB                                   :: OnlyA
            ED + B <--> EBD                                   :: OnlyA
            EAD + B <--> EPQ                                :: OnlyA
            E + D <--> ED                                         :: OnlyA
            EA + D <--> EAD                                   :: OnlyA
            EB + D <--> EBD                                   :: OnlyA
            EAB + D <--> EPQ                                :: OnlyA
            EAB + P ⇌ EABP                                :: OnlyA
            (EB + P ⇌ EBP, EQ + P ⇌ EPQ)                :: OnlyA
            E + Q <--> EQ                                         :: EqualAI
        end
    end
    p_by_b = @enzyme_mechanism begin
        substrates: A, B, D
        products: P, Q
        steps: begin
            E + A <--> EA                                         :: OnlyA
            (EB + A <--> EAB, EBP + A <--> EABP)    :: OnlyA
            ED + A <--> EAD                                   :: OnlyA
            EBD + A <--> EPQ                                :: OnlyA
            E + B <--> EB                                         :: OnlyA
            EA + B <--> EAB                                   :: OnlyA
            ED + B <--> EBD                                   :: OnlyA
            EAD + B <--> EPQ                                :: OnlyA
            E + D <--> ED                                         :: OnlyA
            EA + D <--> EAD                                   :: OnlyA
            EB + D <--> EBD                                   :: OnlyA
            EAB + D <--> EPQ                                :: OnlyA
            EAB + P ⇌ EABP                                :: OnlyA
            EB + P ⇌ EBP                                      :: OnlyA
            EQ + P ⇌ EPQ                                      :: OnlyA
            E + Q <--> EQ                                         :: EqualAI
        end
    end
    @test EnzymeRates.Derivation.rate_equation_term_count(am) >
          EnzymeRates.Derivation.MAX_RATE_EQUATION_TERMS
    # Both conformations count: the inactive conformation keeps only E + Q <--> EQ.
    for x in (am, a_by_p, p_by_b)
        runs = EnzymeRates.Mechanisms.steps_running_in_allosteric_conformations(x)
        @test _testhelper_laplacian_tree_count(x.catalytic, runs.inactive) == 2
        @test EnzymeRates.Derivation.rate_equation_term_count(x) ==
              _testhelper_laplacian_tree_count(x.catalytic, runs.active) + 2
    end
    @test EnzymeRates.Derivation.fitted_parameter_count(am) == 22
    @test EnzymeRates.Derivation.fitted_parameter_count(a_by_p) == 24
    @test EnzymeRates.Derivation.fitted_parameter_count(p_by_b) == 23
    # The test's own count of each model's dimension, which derives no law either.
    for x in (am, a_by_p, p_by_b)
        @test EnzymeRates.Derivation.fitted_parameter_count(x) ==
              _testhelper_declared_dimension(x)
    end
    kids = _testhelper_move(EnzymeRates.Enumeration._expand_split_kinetic_group, am)
    @test length(kids) == 2
    @test Set(kids) == Set([a_by_p, p_by_b])
    others = [_testhelper_move(EnzymeRates.Enumeration._expand_re_to_ss, am);
              _testhelper_move(EnzymeRates.Enumeration._expand_change_allo_state, am)]
    @test Set(_testhelper_children(ter_bi_rxn, [am])) ==
          Set([others; a_by_p; p_by_b])
end

@testset "Allosteric mechanism — a stopped :OnlyA binding strands EA" begin
    # Merged random-order bi-bi with the abortive complex EAQ. The steady-state A group
    # binds A at E and, fused, at EB into EPQ, all :OnlyA. In the inactive
    # conformation the :EqualAI steps E + Q <--> EQ, EQ + A <--> EAQ and
    # EA + Q <--> EAQ still join E to EA, the two sides of the :OnlyA binding
    # E + A <--> EA, so that conformation stops it. At zero P and Q, EA then leaves
    # only by taking up Q: it holds enzyme apart from free E, so the steady state there is
    # not unique (`_NONUNIQUE`), as both oracles find. A split's parts inherit their
    # group's tag, so no split has a limit either. The A group tagged :NonequalAI restores
    # the binding in the inactive conformation, and the guard accepts that child. Every
    # other group tagged alone is refused: B's group, Q's group and Q's abortive binding of
    # A keep EA stranded; P's group joins the two sides of the rapid-equilibrium binding
    # EA + B ⇌ EPQ, which then cannot be switched off (`_FAST_STEP_PINNED`). The search
    # finds every set built from these four refused too, so the A group's child is the
    # move's only one. The flip move emits nothing. Four split sets raise the count, but
    # the split search asks the guard on each child, so they count as failed, and no
    # superset passes: the move emits no split.
    am = @_testhelper_unguarded @enzyme_mechanism begin
        substrates: A, B
        products: P, Q
        steps: begin
            (E + A <--> EA, EB + A <--> EPQ)                         :: OnlyA
            (E + B ⇌ EB, EA + B ⇌ EPQ)                               :: OnlyA
            (E + P <--> EP, EQ + P <--> EPQ)                         :: OnlyA
            (E + Q <--> EQ, EA + Q <--> EAQ, EP + Q <--> EPQ)  :: EqualAI
            EQ + A <--> EAQ                                            :: EqualAI
        end
    end
    D = EnzymeRates.Derivation
    v = D._guard(am)
    @test v.fault === D._NONUNIQUE
    @test EnzymeRates.Mechanisms.step_text(am.catalytic, v.edge) == "EA + Q <--> EAQ"
    @test D._face_names(am.catalytic, D._faces(am.catalytic)[v.face]) == [:P, :Q]
    @test all(_testhelper_admission_agreement(am; rng = Random.MersenneTwister(2922)))
    msg = D._guard_message(am)
    @test_throws ErrorException(msg) @enzyme_mechanism begin
        substrates: A, B
        products: P, Q
        steps: begin
            (E + A <--> EA, EB + A <--> EPQ)                         :: OnlyA
            (E + B ⇌ EB, EA + B ⇌ EPQ)                               :: OnlyA
            (E + P <--> EP, EQ + P <--> EPQ)                         :: OnlyA
            (E + Q <--> EQ, EA + Q <--> EAQ, EP + Q <--> EPQ)  :: EqualAI
            EQ + A <--> EAQ                                            :: EqualAI
        end
    end
    p_relaxed = @_testhelper_unguarded @enzyme_mechanism begin
        substrates: A, B ; products: P, Q ; oligomeric_state: 1
        steps: begin
            (E + A <--> EA, EB + A <--> EPQ)                        :: OnlyA
            (E + B ⇌ EB, EA + B ⇌ EPQ)                              :: OnlyA
            (E + P <--> EP, EQ + P <--> EPQ)                        :: NonequalAI
            (E + Q <--> EQ, EA + Q <--> EAQ, EP + Q <--> EPQ) :: EqualAI
            EQ + A <--> EAQ                                           :: EqualAI
        end
    end
    q_relaxed = @_testhelper_unguarded @enzyme_mechanism begin
        substrates: A, B ; products: P, Q ; oligomeric_state: 1
        steps: begin
            (E + A <--> EA, EB + A <--> EPQ)                        :: OnlyA
            (E + B ⇌ EB, EA + B ⇌ EPQ)                              :: OnlyA
            (E + P <--> EP, EQ + P <--> EPQ)                        :: OnlyA
            (E + Q <--> EQ, EA + Q <--> EAQ, EP + Q <--> EPQ) :: NonequalAI
            EQ + A <--> EAQ                                           :: EqualAI
        end
    end
    qa_relaxed = @_testhelper_unguarded @enzyme_mechanism begin
        substrates: A, B ; products: P, Q ; oligomeric_state: 1
        steps: begin
            (E + A <--> EA, EB + A <--> EPQ)                        :: OnlyA
            (E + B ⇌ EB, EA + B ⇌ EPQ)                              :: OnlyA
            (E + P <--> EP, EQ + P <--> EPQ)                        :: OnlyA
            (E + Q <--> EQ, EA + Q <--> EAQ, EP + Q <--> EPQ) :: EqualAI
            EQ + A <--> EAQ                                           :: NonequalAI
        end
    end
    b_relaxed = @_testhelper_unguarded @enzyme_mechanism begin
        substrates: A, B ; products: P, Q ; oligomeric_state: 1
        steps: begin
            (E + A <--> EA, EB + A <--> EPQ)                        :: OnlyA
            (E + B ⇌ EB, EA + B ⇌ EPQ)                              :: NonequalAI
            (E + P <--> EP, EQ + P <--> EPQ)                        :: OnlyA
            (E + Q <--> EQ, EA + Q <--> EAQ, EP + Q <--> EPQ) :: EqualAI
            EQ + A <--> EAQ                                           :: EqualAI
        end
    end
    a_relaxed = @enzyme_mechanism begin
        substrates: A, B ; products: P, Q ; oligomeric_state: 1
        steps: begin
            (E + A <--> EA, EB + A <--> EPQ)                        :: NonequalAI
            (E + B ⇌ EB, EA + B ⇌ EPQ)                              :: OnlyA
            (E + P <--> EP, EQ + P <--> EPQ)                        :: OnlyA
            (E + Q <--> EQ, EA + Q <--> EAQ, EP + Q <--> EPQ) :: EqualAI
            EQ + A <--> EAQ                                           :: EqualAI
        end
    end
    relaxed = _testhelper_move(EnzymeRates.Enumeration._expand_change_allo_state, am)
    @test length(relaxed) == 1
    @test Set(relaxed) == Set([a_relaxed])
    split = _testhelper_move(EnzymeRates.Enumeration._expand_split_kinetic_group, am)
    @test length(split) == 0
    @test Set(split) == Set()
    flips = _testhelper_move(EnzymeRates.Enumeration._expand_re_to_ss, am)
    @test length(flips) == 0
    @test Set(flips) == Set()
    refused = [b_relaxed, p_relaxed, q_relaxed, qa_relaxed]
    @test Dict(c => D._guard(c).fault for c in refused) ==
          Dict(b_relaxed => D._NONUNIQUE, p_relaxed => D._FAST_STEP_PINNED,
               q_relaxed => D._NONUNIQUE, qa_relaxed => D._NONUNIQUE)
    kids = _testhelper_children(bi_bi_rxn, [am])
    @test length(kids) == 1
    @test Set(kids) == Set([a_relaxed])
end

@testset "Allosteric mechanism — beside an inhibitor, stopped :OnlyA A strands EA" begin
    # Merged random-order bi-bi with the competitive inhibitor I and the abortive
    # complexes EAP and EBQ; the A, B, I and Q groups :OnlyA. In the inactive
    # conformation the :EqualAI steps E + P <--> EP, EP + A <--> EAP and
    # EA + P <--> EAP join E to EA, the two sides of the :OnlyA binding
    # E + A <--> EA, so that conformation stops it. At zero P and Q, EA then leaves
    # only by taking up P: it holds enzyme apart from free E, so the steady state there is
    # not unique (`_NONUNIQUE`), as both oracles find. A split's parts inherit their
    # group's tag, so no split has a limit either. The A group tagged :NonequalAI restores
    # the binding in the inactive conformation, and the guard accepts that child. Every
    # other group tagged alone is refused: B's group, I's group, P's group and P's abortive
    # binding of A keep EA stranded; Q's group joins the two sides of the
    # rapid-equilibrium binding EA + B ⇌ EPQ, which then cannot be switched off
    # (`_FAST_STEP_PINNED`). The search finds every set built from these five refused too,
    # so the A group's child is the move's only one. The flip move emits nothing. Six split
    # sets raise the count, but the split search asks the guard on each child, so they
    # count as failed, and no superset passes: the move emits no split.
    am = @_testhelper_unguarded @enzyme_mechanism begin
        substrates: A, B
        products: P, Q
        modifiers: I
        steps: begin
            (E + A <--> EA, EB + A <--> EPQ, EI + A <--> EAI)   :: OnlyA
            (E + B ⇌ EB, EA + B ⇌ EPQ, EQ + B ⇌ EBQ)            :: OnlyA
            (E + I ⇌ EI, EA + I ⇌ EAI)                                :: OnlyA
            (E + P <--> EP, EA + P <--> EAP, EQ + P <--> EPQ)   :: EqualAI
            (E + Q <--> EQ, EB + Q <--> EBQ, EP + Q <--> EPQ)   :: OnlyA
            EP + A <--> EAP                                             :: EqualAI
        end
    end
    D = EnzymeRates.Derivation
    v = D._guard(am)
    @test v.fault === D._NONUNIQUE
    @test EnzymeRates.Mechanisms.step_text(am.catalytic, v.edge) == "EA + P <--> EAP"
    @test D._face_names(am.catalytic, D._faces(am.catalytic)[v.face]) == [:P, :Q]
    @test all(_testhelper_admission_agreement(am; rng = Random.MersenneTwister(95)))
    msg = D._guard_message(am)
    @test_throws ErrorException(msg) @enzyme_mechanism begin
        substrates: A, B
        products: P, Q
        modifiers: I
        steps: begin
            (E + A <--> EA, EB + A <--> EPQ, EI + A <--> EAI)   :: OnlyA
            (E + B ⇌ EB, EA + B ⇌ EPQ, EQ + B ⇌ EBQ)            :: OnlyA
            (E + I ⇌ EI, EA + I ⇌ EAI)                                :: OnlyA
            (E + P <--> EP, EA + P <--> EAP, EQ + P <--> EPQ)   :: EqualAI
            (E + Q <--> EQ, EB + Q <--> EBQ, EP + Q <--> EPQ)   :: OnlyA
            EP + A <--> EAP                                             :: EqualAI
        end
    end
    rxn = @enzyme_reaction begin
        substrates: A[C], B[N]
        products: P[C], Q[N]
        competitive_inhibitors: I
    end
    i_relaxed = @_testhelper_unguarded @enzyme_mechanism begin
        substrates: A, B ; products: P, Q ; modifiers: I ; oligomeric_state: 1
        steps: begin
            (E + A <--> EA, EB + A <--> EPQ, EI + A <--> EAI) :: OnlyA
            (E + B ⇌ EB, EA + B ⇌ EPQ, EQ + B ⇌ EBQ)          :: OnlyA
            (E + I ⇌ EI, EA + I ⇌ EAI)                              :: NonequalAI
            (E + P <--> EP, EA + P <--> EAP, EQ + P <--> EPQ) :: EqualAI
            (E + Q <--> EQ, EB + Q <--> EBQ, EP + Q <--> EPQ) :: OnlyA
            EP + A <--> EAP                                           :: EqualAI
        end
    end
    p_relaxed = @_testhelper_unguarded @enzyme_mechanism begin
        substrates: A, B ; products: P, Q ; modifiers: I ; oligomeric_state: 1
        steps: begin
            (E + A <--> EA, EB + A <--> EPQ, EI + A <--> EAI) :: OnlyA
            (E + B ⇌ EB, EA + B ⇌ EPQ, EQ + B ⇌ EBQ)          :: OnlyA
            (E + I ⇌ EI, EA + I ⇌ EAI)                              :: OnlyA
            (E + P <--> EP, EA + P <--> EAP, EQ + P <--> EPQ) :: NonequalAI
            (E + Q <--> EQ, EB + Q <--> EBQ, EP + Q <--> EPQ) :: OnlyA
            EP + A <--> EAP                                           :: EqualAI
        end
    end
    q_relaxed = @_testhelper_unguarded @enzyme_mechanism begin
        substrates: A, B ; products: P, Q ; modifiers: I ; oligomeric_state: 1
        steps: begin
            (E + A <--> EA, EB + A <--> EPQ, EI + A <--> EAI) :: OnlyA
            (E + B ⇌ EB, EA + B ⇌ EPQ, EQ + B ⇌ EBQ)          :: OnlyA
            (E + I ⇌ EI, EA + I ⇌ EAI)                              :: OnlyA
            (E + P <--> EP, EA + P <--> EAP, EQ + P <--> EPQ) :: EqualAI
            (E + Q <--> EQ, EB + Q <--> EBQ, EP + Q <--> EPQ) :: NonequalAI
            EP + A <--> EAP                                           :: EqualAI
        end
    end
    pa_relaxed = @_testhelper_unguarded @enzyme_mechanism begin
        substrates: A, B ; products: P, Q ; modifiers: I ; oligomeric_state: 1
        steps: begin
            (E + A <--> EA, EB + A <--> EPQ, EI + A <--> EAI) :: OnlyA
            (E + B ⇌ EB, EA + B ⇌ EPQ, EQ + B ⇌ EBQ)          :: OnlyA
            (E + I ⇌ EI, EA + I ⇌ EAI)                              :: OnlyA
            (E + P <--> EP, EA + P <--> EAP, EQ + P <--> EPQ) :: EqualAI
            (E + Q <--> EQ, EB + Q <--> EBQ, EP + Q <--> EPQ) :: OnlyA
            EP + A <--> EAP                                           :: NonequalAI
        end
    end
    b_relaxed = @_testhelper_unguarded @enzyme_mechanism begin
        substrates: A, B ; products: P, Q ; modifiers: I ; oligomeric_state: 1
        steps: begin
            (E + A <--> EA, EB + A <--> EPQ, EI + A <--> EAI) :: OnlyA
            (E + B ⇌ EB, EA + B ⇌ EPQ, EQ + B ⇌ EBQ)          :: NonequalAI
            (E + I ⇌ EI, EA + I ⇌ EAI)                              :: OnlyA
            (E + P <--> EP, EA + P <--> EAP, EQ + P <--> EPQ) :: EqualAI
            (E + Q <--> EQ, EB + Q <--> EBQ, EP + Q <--> EPQ) :: OnlyA
            EP + A <--> EAP                                           :: EqualAI
        end
    end
    a_relaxed = @enzyme_mechanism begin
        substrates: A, B ; products: P, Q ; modifiers: I ; oligomeric_state: 1
        steps: begin
            (E + A <--> EA, EB + A <--> EPQ, EI + A <--> EAI) :: NonequalAI
            (E + B ⇌ EB, EA + B ⇌ EPQ, EQ + B ⇌ EBQ)          :: OnlyA
            (E + I ⇌ EI, EA + I ⇌ EAI)                              :: OnlyA
            (E + P <--> EP, EA + P <--> EAP, EQ + P <--> EPQ) :: EqualAI
            (E + Q <--> EQ, EB + Q <--> EBQ, EP + Q <--> EPQ) :: OnlyA
            EP + A <--> EAP                                           :: EqualAI
        end
    end
    relaxed = _testhelper_move(EnzymeRates.Enumeration._expand_change_allo_state, am)
    @test length(relaxed) == 1
    @test Set(relaxed) == Set([a_relaxed])
    split = _testhelper_move(EnzymeRates.Enumeration._expand_split_kinetic_group, am)
    @test length(split) == 0
    @test Set(split) == Set()
    flips = _testhelper_move(EnzymeRates.Enumeration._expand_re_to_ss, am)
    @test length(flips) == 0
    @test Set(flips) == Set()
    refused = [b_relaxed, i_relaxed, p_relaxed, q_relaxed, pa_relaxed]
    @test Dict(c => D._guard(c).fault for c in refused) ==
          Dict(b_relaxed => D._NONUNIQUE, i_relaxed => D._NONUNIQUE,
               p_relaxed => D._NONUNIQUE, q_relaxed => D._FAST_STEP_PINNED,
               pa_relaxed => D._NONUNIQUE)
    kids = _testhelper_children(rxn, [am])
    @test length(kids) == 1
    @test Set(kids) == Set([a_relaxed])
end

@testset "Allosteric mechanism — a split keeps its parent's interior verdict" begin
    # An aggregate check over the bi-bi seed set (the exception to writing every fixture
    # inline): each split child of a seed tagged as below, one parent in
    # ten. Both parts of a split group inherit its tag, and admission's interior reads
    # each edge's tag, so every split the move builds (zero-flux parts at rapid
    # equilibrium) has its parent's interior verdict: the split search's guard rejects
    # every split of a parent the interior refuses.
    # Each split child's every-chemistry tag set with each nonempty set of its
    # binding groups :OnlyA, children of the tagging move, among them the parents
    # admission refuses, which the reject direction needs.
    function _testhelper_unfiltered_tags(c)
        chem = EnzymeRates.Enumeration._chemistry_groups(_testhelper_candidate(c))
        vt = [x ? :OnlyA : :EqualAI for x in chem]
        [(t = copy(vt); t[s] .= :OnlyA; t)
         for s in EnzymeRates.Enumeration._subsets(findall(!, chem))[2:end]]
    end
    parents = unique(EnzymeRates.Mechanism[
        EnzymeRates.Mechanism(c.catalytic, EnzymeRates.Mechanisms.Allostery(
            t, 1, EnzymeRates.Mechanisms.RegulatorySite[]))
        for m in _testhelper_mechanisms(EnzymeRates.Enumeration.init_mechanisms(bi_bi_rxn))
        for c in _testhelper_move(EnzymeRates.Enumeration._expand_split_kinetic_group, m)
        for t in _testhelper_unfiltered_tags(c)])[1:10:end]
    accepted = rejected = 0
    for am in parents
        cm = am.catalytic
        verdict = EnzymeRates.Derivation._admission_interior(
            am, EnzymeRates.Derivation._built_reading(am)).fault
        flux = EnzymeRates.Enumeration._flux_carrying_edges(
            cm, _testhelper_flags(cm))
        units = [(g, bp) for g in 1:length(_testhelper_group_members(cm))
                 for bp in EnzymeRates.Enumeration._context_bipartitions(
                     _testhelper_candidate(cm), g)]
        reverted = [EnzymeRates.Enumeration._revert_zero_flux_parts(cm, bp, flux)
                    for (_, bp) in units]
        for sel in [[[u] for u in eachindex(units)];
                    [[u, v] for u in eachindex(units) for v in u+1:length(units)
                     if units[u][1] != units[v][1]]]
            re = _testhelper_flags(cm)
            foreach(u -> re[reverted[u]] .= true, sel)
            child = EnzymeRates.Enumeration._apply_bipartitions(_testhelper_candidate(am),
                                                                units[sel], re).mechanism
            @test EnzymeRates.Derivation._admission_interior(
                child, EnzymeRates.Derivation._built_reading(child)).fault === verdict
        end
        verdict === EnzymeRates.Derivation._NO_FAULT ? (accepted += 1) : (rejected += 1)
    end
    # Aggregate pin over the sampled parents: both verdicts occur, so the agreement
    # means something.
    @test (accepted, rejected) == (761, 227)
    @test all(
        am -> _testhelper_move(EnzymeRates.Enumeration._expand_split_kinetic_group, am) isa
              Vector, parents)
end

@testset "Allosteric mechanism — an :EqualAI group counts once across conformations" begin
    # Bi-bi with A bound before B, P and Q released in either order, and the abortive
    # complexes EAP and EBQ; the A group and the chemistry are :OnlyA. The
    # inactive conformation binds no A, so of the :EqualAI B group it reaches only
    # EQ + B, yet the B constant is one class in both conformations: the parent fits six
    # constants, its dimension. Splitting the B group gains one, since EBQ lies on no
    # cycle. The square E, EA, EAP, EP, closed by the shared A and P groups, ties
    # back the split of either group by the other, and the square E, EP, EPQ, EQ
    # does the same for the P and Q groups, so splitting the P group by A or by Q gains
    # nothing alone; with the A or the Q group split too, it gains one. The move emits
    # those three children, each fitting as many constants as its dimension.
    am = @enzyme_mechanism begin
        substrates: A, B
        products: P, Q
        oligomeric_state: 2
        steps: begin
            (E + A ⇌ EA, EP + A ⇌ EAP)                          :: OnlyA
            (E + P ⇌ EP, EA + P ⇌ EAP, EQ + P ⇌ EPQ)      :: EqualAI
            (E + Q ⇌ EQ, EP + Q ⇌ EPQ)                          :: EqualAI
            (EA + B ⇌ EAB, EQ + B ⇌ EBQ)                    :: EqualAI
            EAB <--> EPQ                                        :: OnlyA
        end
    end
    b_by_q = @enzyme_mechanism begin
        substrates: A, B
        products: P, Q
        oligomeric_state: 2
        steps: begin
            (E + A ⇌ EA, EP + A ⇌ EAP)                          :: OnlyA
            (E + P ⇌ EP, EA + P ⇌ EAP, EQ + P ⇌ EPQ)      :: EqualAI
            (E + Q ⇌ EQ, EP + Q ⇌ EPQ)                          :: EqualAI
            EA + B ⇌ EAB                                          :: EqualAI
            EAB <--> EPQ                                        :: OnlyA
            EQ + B ⇌ EBQ                                          :: EqualAI
        end
    end
    a_by_p_and_p_by_a = @enzyme_mechanism begin
        substrates: A, B
        products: P, Q
        oligomeric_state: 2
        steps: begin
            E + A ⇌ EA                                                :: OnlyA
            (E + P ⇌ EP, EQ + P ⇌ EPQ)                          :: EqualAI
            (E + Q ⇌ EQ, EP + Q ⇌ EPQ)                          :: EqualAI
            (EA + B ⇌ EAB, EQ + B ⇌ EBQ)                    :: EqualAI
            EA + P ⇌ EAP                                          :: EqualAI
            EAB <--> EPQ                                        :: OnlyA
            EP + A ⇌ EAP                                          :: OnlyA
        end
    end
    q_by_p_and_p_by_q = @enzyme_mechanism begin
        substrates: A, B
        products: P, Q
        oligomeric_state: 2
        steps: begin
            (E + A ⇌ EA, EP + A ⇌ EAP)                          :: OnlyA
            (E + P ⇌ EP, EA + P ⇌ EAP)                          :: EqualAI
            E + Q ⇌ EQ                                                :: EqualAI
            (EA + B ⇌ EAB, EQ + B ⇌ EBQ)                    :: EqualAI
            EAB <--> EPQ                                        :: OnlyA
            EP + Q ⇌ EPQ                                          :: EqualAI
            EQ + P ⇌ EPQ                                          :: EqualAI
        end
    end
    p_by_a = @enzyme_mechanism begin
        substrates: A, B
        products: P, Q
        oligomeric_state: 2
        steps: begin
            (E + A ⇌ EA, EP + A ⇌ EAP)                          :: OnlyA
            (E + P ⇌ EP, EQ + P ⇌ EPQ)                          :: EqualAI
            (E + Q ⇌ EQ, EP + Q ⇌ EPQ)                          :: EqualAI
            (EA + B ⇌ EAB, EQ + B ⇌ EBQ)                    :: EqualAI
            EA + P ⇌ EAP                                          :: EqualAI
            EAB <--> EPQ                                        :: OnlyA
        end
    end
    p_by_q = @enzyme_mechanism begin
        substrates: A, B
        products: P, Q
        oligomeric_state: 2
        steps: begin
            (E + A ⇌ EA, EP + A ⇌ EAP)                          :: OnlyA
            (E + P ⇌ EP, EA + P ⇌ EAP)                          :: EqualAI
            (E + Q ⇌ EQ, EP + Q ⇌ EPQ)                          :: EqualAI
            (EA + B ⇌ EAB, EQ + B ⇌ EBQ)                    :: EqualAI
            EAB <--> EPQ                                        :: OnlyA
            EQ + P ⇌ EPQ                                          :: EqualAI
        end
    end
    kids = _testhelper_move(EnzymeRates.Enumeration._expand_split_kinetic_group, am)
    @test length(kids) == 3
    @test Set(kids) == Set([b_by_q, a_by_p_and_p_by_a, q_by_p_and_p_by_q])
    @test !(p_by_a in kids)
    @test !(p_by_q in kids)
    # The test's own count of each model's dimension, which derives no law.
    for (x, n) in ((am, 6), (b_by_q, 7), (a_by_p_and_p_by_a, 7), (q_by_p_and_p_by_q, 7),
                   (p_by_a, 6), (p_by_q, 6))
        @test EnzymeRates.Derivation.fitted_parameter_count(x) ==
              _testhelper_declared_dimension(x) == n
    end
end

end

# ─── _expand_add_dead_end_regulator ────────────────────────────────────
@testset "_expand_add_dead_end_regulator" begin

@testset "Mechanism — bi-bi + I: variants" begin
    # Each child is the seed with the foreign inhibitor I bound at the forms of one
    # competition pattern's placement, plus a mirror of every seed step that joins two
    # of those forms in that step's group. I has no productive twin, so every placement
    # is kept.
    for (m, rxn, n_variants, multi_form, expected) in [
        # SEED: bi-bi sequential. The expansion should produce 4 form sets after dedup:
        # {E, EQ}, {E}, {EA, EQ} and {E, EA}.
        (@enzyme_mechanism(begin
            substrates: A, B
            products: P, Q
            steps: begin
                E + A ⇌ EA
                EA + B ⇌ EAB
                E + Q ⇌ EQ
                EQ + P ⇌ EPQ
                EAB <--> EPQ
            end
        end), @enzyme_reaction(begin
            substrates: A[C], B[N]
            products: P[C], Q[N]
            competitive_inhibitors: I
        end), 4, false, [
            @enzyme_mechanism(begin
                substrates: A, B; products: P, Q; modifiers: I
                steps: begin
                    E + A ⇌ EA
                    EA + B ⇌ EAB
                    (E + Q ⇌ EQ, EI + Q ⇌ EIQ)
                    EQ + P ⇌ EPQ
                    EAB <--> EPQ
                    (E + I ⇌ EI, EQ + I ⇌ EIQ)
                end
            end),
            @enzyme_mechanism(begin
                substrates: A, B; products: P, Q; modifiers: I
                steps: begin
                    E + A ⇌ EA
                    EA + B ⇌ EAB
                    E + Q ⇌ EQ
                    EQ + P ⇌ EPQ
                    EAB <--> EPQ
                    E + I ⇌ EI
                end
            end),
            @enzyme_mechanism(begin
                substrates: A, B; products: P, Q; modifiers: I
                steps: begin
                    E + A ⇌ EA
                    EA + B ⇌ EAB
                    E + Q ⇌ EQ
                    EQ + P ⇌ EPQ
                    EAB <--> EPQ
                    (EA + I ⇌ EAI, EQ + I ⇌ EIQ)
                end
            end),
            @enzyme_mechanism(begin
                substrates: A, B; products: P, Q; modifiers: I
                steps: begin
                    (E + A ⇌ EA, EI + A ⇌ EAI)
                    EA + B ⇌ EAB
                    E + Q ⇌ EQ
                    EQ + P ⇌ EPQ
                    EAB <--> EPQ
                    (E + I ⇌ EI, EA + I ⇌ EAI)
                end
            end)]),
        # SEED: bi-bi random. It has 5 eligible forms (E, E_A, E_B, E_P, E_Q);
        # competition patterns × dedup → 9 variants. Some variants bind I at several
        # forms: {E, EB, EQ}, {E, EB, EP}, {E, EB}, {E, EA, EQ},
        # {E, EA, EP}, {E, EA}, {E, EQ}, {E, EP} and {E}.
        (@enzyme_mechanism(begin
            substrates: A, B
            products: P, Q
            steps: begin
                (E + A ⇌ EA, EB + A ⇌ EAB)
                (E + B ⇌ EB, EA + B ⇌ EAB)
                (E + P ⇌ EP, EQ + P ⇌ EPQ)
                (E + Q ⇌ EQ, EP + Q ⇌ EPQ)
                EAB <--> EPQ
            end
        end), @enzyme_reaction(begin
            substrates: A[C], B[N]
            products: P[C], Q[N]
            competitive_inhibitors: I
        end), 9, true, [
            @enzyme_mechanism(begin
                substrates: A, B; products: P, Q; modifiers: I
                steps: begin
                    (E + A ⇌ EA, EB + A ⇌ EAB)
                    (E + B ⇌ EB, EA + B ⇌ EAB, EI + B ⇌ EBI)
                    (E + P ⇌ EP, EQ + P ⇌ EPQ)
                    (E + Q ⇌ EQ, EP + Q ⇌ EPQ, EI + Q ⇌ EIQ)
                    EAB <--> EPQ
                    (E + I ⇌ EI, EB + I ⇌ EBI, EQ + I ⇌ EIQ)
                end
            end),
            @enzyme_mechanism(begin
                substrates: A, B; products: P, Q; modifiers: I
                steps: begin
                    (E + A ⇌ EA, EB + A ⇌ EAB)
                    (E + B ⇌ EB, EA + B ⇌ EAB, EI + B ⇌ EBI)
                    (E + P ⇌ EP, EQ + P ⇌ EPQ, EI + P ⇌ EIP)
                    (E + Q ⇌ EQ, EP + Q ⇌ EPQ)
                    EAB <--> EPQ
                    (E + I ⇌ EI, EB + I ⇌ EBI, EP + I ⇌ EIP)
                end
            end),
            @enzyme_mechanism(begin
                substrates: A, B; products: P, Q; modifiers: I
                steps: begin
                    (E + A ⇌ EA, EB + A ⇌ EAB)
                    (E + B ⇌ EB, EA + B ⇌ EAB, EI + B ⇌ EBI)
                    (E + P ⇌ EP, EQ + P ⇌ EPQ)
                    (E + Q ⇌ EQ, EP + Q ⇌ EPQ)
                    EAB <--> EPQ
                    (E + I ⇌ EI, EB + I ⇌ EBI)
                end
            end),
            @enzyme_mechanism(begin
                substrates: A, B; products: P, Q; modifiers: I
                steps: begin
                    (E + A ⇌ EA, EB + A ⇌ EAB, EI + A ⇌ EAI)
                    (E + B ⇌ EB, EA + B ⇌ EAB)
                    (E + P ⇌ EP, EQ + P ⇌ EPQ)
                    (E + Q ⇌ EQ, EP + Q ⇌ EPQ, EI + Q ⇌ EIQ)
                    EAB <--> EPQ
                    (E + I ⇌ EI, EA + I ⇌ EAI, EQ + I ⇌ EIQ)
                end
            end),
            @enzyme_mechanism(begin
                substrates: A, B; products: P, Q; modifiers: I
                steps: begin
                    (E + A ⇌ EA, EB + A ⇌ EAB, EI + A ⇌ EAI)
                    (E + B ⇌ EB, EA + B ⇌ EAB)
                    (E + P ⇌ EP, EQ + P ⇌ EPQ, EI + P ⇌ EIP)
                    (E + Q ⇌ EQ, EP + Q ⇌ EPQ)
                    EAB <--> EPQ
                    (E + I ⇌ EI, EA + I ⇌ EAI, EP + I ⇌ EIP)
                end
            end),
            @enzyme_mechanism(begin
                substrates: A, B; products: P, Q; modifiers: I
                steps: begin
                    (E + A ⇌ EA, EB + A ⇌ EAB, EI + A ⇌ EAI)
                    (E + B ⇌ EB, EA + B ⇌ EAB)
                    (E + P ⇌ EP, EQ + P ⇌ EPQ)
                    (E + Q ⇌ EQ, EP + Q ⇌ EPQ)
                    EAB <--> EPQ
                    (E + I ⇌ EI, EA + I ⇌ EAI)
                end
            end),
            @enzyme_mechanism(begin
                substrates: A, B; products: P, Q; modifiers: I
                steps: begin
                    (E + A ⇌ EA, EB + A ⇌ EAB)
                    (E + B ⇌ EB, EA + B ⇌ EAB)
                    (E + P ⇌ EP, EQ + P ⇌ EPQ)
                    (E + Q ⇌ EQ, EP + Q ⇌ EPQ, EI + Q ⇌ EIQ)
                    EAB <--> EPQ
                    (E + I ⇌ EI, EQ + I ⇌ EIQ)
                end
            end),
            @enzyme_mechanism(begin
                substrates: A, B; products: P, Q; modifiers: I
                steps: begin
                    (E + A ⇌ EA, EB + A ⇌ EAB)
                    (E + B ⇌ EB, EA + B ⇌ EAB)
                    (E + P ⇌ EP, EQ + P ⇌ EPQ, EI + P ⇌ EIP)
                    (E + Q ⇌ EQ, EP + Q ⇌ EPQ)
                    EAB <--> EPQ
                    (E + I ⇌ EI, EP + I ⇌ EIP)
                end
            end),
            @enzyme_mechanism(begin
                substrates: A, B; products: P, Q; modifiers: I
                steps: begin
                    (E + A ⇌ EA, EB + A ⇌ EAB)
                    (E + B ⇌ EB, EA + B ⇌ EAB)
                    (E + P ⇌ EP, EQ + P ⇌ EPQ)
                    (E + Q ⇌ EQ, EP + Q ⇌ EPQ)
                    EAB <--> EPQ
                    E + I ⇌ EI
                end
            end)]),
        # SEED: bi-bi ping-pong with Estar. It has 4 eligible forms (E, E_A, Estar,
        # E_Q); competition patterns × dedup → 3 variants: {E, Estar}, {E} and {Estar}.
        (@enzyme_mechanism(begin
            substrates: A, B
            products: P, Q
            steps: begin
                E + A ⇌ EA
                Estar + B ⇌ EstarB
                E + Q ⇌ EQ
                Estar + P ⇌ EstarP
                EA <--> EstarP
                EstarB ⇌ EQ
            end
        end), @enzyme_reaction(begin
            substrates: A[CX], B[N]
            products: P[C], Q[NX]
            competitive_inhibitors: I
        end), 3, false, [
            @enzyme_mechanism(begin
                substrates: A, B; products: P, Q; modifiers: I
                steps: begin
                    E + A ⇌ EA
                    Estar + B ⇌ EstarB
                    E + Q ⇌ EQ
                    Estar + P ⇌ EstarP
                    EA <--> EstarP
                    EstarB ⇌ EQ
                    (E + I ⇌ EI, Estar + I ⇌ EstarI)
                end
            end),
            @enzyme_mechanism(begin
                substrates: A, B; products: P, Q; modifiers: I
                steps: begin
                    E + A ⇌ EA
                    Estar + B ⇌ EstarB
                    E + Q ⇌ EQ
                    Estar + P ⇌ EstarP
                    EA <--> EstarP
                    EstarB ⇌ EQ
                    E + I ⇌ EI
                end
            end),
            @enzyme_mechanism(begin
                substrates: A, B; products: P, Q; modifiers: I
                steps: begin
                    E + A ⇌ EA
                    Estar + B ⇌ EstarB
                    E + Q ⇌ EQ
                    Estar + P ⇌ EstarP
                    EA <--> EstarP
                    EstarB ⇌ EQ
                    Estar + I ⇌ EstarI
                end
            end)]),
    ]
        _testhelper_assert_mechanism_invariants(m)

        result = _testhelper_move(EnzymeRates.Enumeration._expand_add_dead_end_regulator, m,
                                  rxn)

        # 1. count and the exact children
        @test length(result) == n_variants
        @test Set(result) == Set(expected)

        # 2. Δ params: +1 each (one new K_I parameter), measured against ground-truth
        # `parameters(...)` — the canonical source
        # for exact parameter counts.
        @test _testhelper_param_deltas(m, result) == fill(1, n_variants)

        # 3. invariants
        for r in result
            _testhelper_assert_mechanism_invariants(r)
        end

        # 4. property: each variant has ≥1 I-binding step, and all I-binding steps in a
        # single variant live in the same kinetic group (outer-vector index) — one K_I,
        # not one per form.
        for r in result
            i_groups = _testhelper_groups_binding(r, :I)
            @test !isempty(i_groups)
            @test length(unique(i_groups)) == 1
        end

        # 5. some variant binds I at multiple forms: ≥2 I-binding steps, in one group.
        multi_form && @test any(
            r -> length(_testhelper_groups_binding(r, :I)) >= 2, result)
    end
end

@testset "Mechanism — dead-end sites on a Theorell–Chance step" begin
    # Ordered Theorell–Chance, all SS: E + A → EA, EA + B → EQ + P, E + Q → EQ.
    # Sites: A is free at E, B at EA (taken up by the TC step), P at EQ (given off
    # there), Q at E. A foreign inhibitor I competes with one or both substrates and one
    # or both products; a site is skipped when it already holds a competing reactant.
    # {A}×{P}: E, EQ. {A}×{Q}: E. {A}×{P,Q}: E (EQ holds Q). {B}×{P}: EA, EQ.
    # {B}×{Q}: E, EA. {B}×{P,Q}: E, EA (EQ holds Q). {A,B}×{P}: E, EQ (EA
    # holds A). {A,B}×{Q}: E. {A,B}×{P,Q}: E. Distinct placements: {E}, {E, EQ},
    # {EA, EQ}, {E, EA}. Mirrors join two sites: E–EQ by Q's binding, E–EA by
    # A's, EA–EQ by the TC step, each in its parent's group, at its parent's flag.
    # A foreign inhibitor has no twins, so the copy rule keeps every placement.
    rxn = @enzyme_reaction begin
        substrates: A[C], B[N]
        products: P[C], Q[N]
        competitive_inhibitors: I
    end
    m = @enzyme_mechanism begin
        substrates: A, B
        products: P, Q
        steps: begin
            E + A <--> EA
            EA + B <--> EQ + P
            E + Q <--> EQ
        end
    end
    at_E = @enzyme_mechanism begin
        substrates: A, B; products: P, Q; modifiers: I
        steps: begin
            E + A <--> EA
            EA + B <--> EQ + P
            E + Q <--> EQ
            E + I ⇌ EI
        end
    end
    at_E_EQ = @enzyme_mechanism begin
        substrates: A, B; products: P, Q; modifiers: I
        steps: begin
            E + A <--> EA
            EA + B <--> EQ + P
            (E + Q <--> EQ, EI + Q <--> EIQ)
            (E + I ⇌ EI, EQ + I ⇌ EIQ)
        end
    end
    at_EA_EQ = @enzyme_mechanism begin
        substrates: A, B; products: P, Q; modifiers: I
        steps: begin
            E + A <--> EA
            (EA + B <--> EQ + P, EAI + B <--> EIQ + P)
            E + Q <--> EQ
            (EA + I ⇌ EAI, EQ + I ⇌ EIQ)
        end
    end
    at_E_EA = @enzyme_mechanism begin
        substrates: A, B; products: P, Q; modifiers: I
        steps: begin
            (E + A <--> EA, EI + A <--> EAI)
            EA + B <--> EQ + P
            E + Q <--> EQ
            (E + I ⇌ EI, EA + I ⇌ EAI)
        end
    end
    kids = _testhelper_move(EnzymeRates.Enumeration._expand_add_dead_end_regulator, m, rxn)
    @test length(kids) == 4
    @test Set(kids) == Set([at_E, at_E_EQ, at_EA_EQ, at_E_EA])
    # The parent's six rate constants less one Haldane relation leave 5 fitted, all
    # identifiable. Each child adds K_I, and the law sees it: at {E} the free enzyme's
    # terms gain the factor 1 + I/K_I; at two sites each site and its complex weigh
    # 1 + I/K_I times the site, the mirror carries the step between the sites
    # unchanged, and each step from a site to the third form is divided by that
    # factor, so I enters squared. 6 fitted, rank 6.
    @test _testhelper_fitted(m) == _testhelper_identifiable_rank(m) == 5
    for k in kids
        @test _testhelper_fitted(k) == _testhelper_identifiable_rank(k) == 6
    end
end

@testset "Mechanism — a copy whose complexes sort among the parent's forms" begin
    # A as its own competitive inhibitor: the complexes where A binds as a copy, named
    # with Amod, sort in among the parent's forms by name, so later forms of the parent
    # move in the child. The copy rule reads each candidate's own twins: A at E alone
    # duplicates EA and is dropped; at EQ, with E + Q mirrored onto EAmod, and at EA and
    # EQ it is new.
    rxn = @enzyme_reaction begin
        substrates: A[C], B[N]
        products: P[C], Q[N]
        competitive_inhibitors: A
    end
    parent = @enzyme_mechanism begin
        substrates: A, B
        products: P, Q
        steps: begin
            E + A ⇌ EA
            EA + B ⇌ EAB
            EAB <--> EPQ
            EQ + P ⇌ EPQ
            E + Q ⇌ EQ
        end
    end
    kids = _testhelper_move(EnzymeRates.Enumeration._expand_add_dead_end_regulator, parent,
                            rxn)
    expected = [
        @enzyme_mechanism(begin
            substrates: A, B
            products: P, Q
            steps: begin
                E + A ⇌ EA
                (E + A ⇌ EAmod, EQ + A ⇌ EAmodQ)
                (E + Q ⇌ EQ, EAmod + Q ⇌ EAmodQ)
                EA + B ⇌ EAB
                EAB <--> EPQ
                EQ + P ⇌ EPQ
            end
        end),
        @enzyme_mechanism(begin
            substrates: A, B
            products: P, Q
            steps: begin
                E + A ⇌ EA
                E + Q ⇌ EQ
                (EA + A ⇌ EAAmod, EQ + A ⇌ EAmodQ)
                EA + B ⇌ EAB
                EAB <--> EPQ
                EQ + P ⇌ EPQ
            end
        end),
        @enzyme_mechanism(begin
            substrates: A, B
            products: P, Q
            steps: begin
                (E + A ⇌ EA, EAmod + A ⇌ EAAmod)
                (E + A ⇌ EAmod, EA + A ⇌ EAAmod)
                E + Q ⇌ EQ
                EA + B ⇌ EAB
                EAB <--> EPQ
                EQ + P ⇌ EPQ
            end
        end)]
    @test length(kids) == 3
    @test Set(kids) == Set(expected)
    # The dropped copy: EAmod at E alone, whose twin EA the A binding forms.
    dropped = @enzyme_mechanism begin
        substrates: A, B
        products: P, Q
        steps: begin
            E + A ⇌ EA
            E + A ⇌ EAmod
            E + Q ⇌ EQ
            EA + B ⇌ EAB
            EAB <--> EPQ
            EQ + P ⇌ EPQ
        end
    end
    @test EnzymeRates.Enumeration._redundant_copy_groups(_testhelper_candidate(dropped)) ==
          [dropped.catalytic.steps[
               _testhelper_edge(dropped, "E + A ⇌ EAmod")].kinetic_group]
    # A bound where it does not react is A taken up: the children hold the parent's
    # metabolites.
    @test parent.catalytic.metabolites == [:A, :B, :P, :Q]
    @test all(k -> k.catalytic.metabolites == [:A, :B, :P, :Q], kids)
end

@testset "Mechanism — Two regulators chain: J added after I preserves I" begin
    # After step A (add I) and step B (add J), J-binding steps must exist
    # and I-binding steps must remain.
    m = @enzyme_mechanism begin
        substrates: S
        products: P
        steps: begin
            E + P ⇌ EP
            E + S ⇌ ES
            ES <--> EP
        end
    end
    _testhelper_assert_mechanism_invariants(m)
    rxn = @enzyme_reaction begin
        substrates: S[C]
        products: P[C]
        competitive_inhibitors: I, J
    end

    # Step A: 2 eligible regs (I, J), 1 form each → 2 variants total.
    i_or_j_ms = _testhelper_move(EnzymeRates.Enumeration._expand_add_dead_end_regulator, m,
                                 rxn)
    @test length(i_or_j_ms) == 2
    with_i = first(filter(r -> !isempty(_testhelper_groups_binding(r, :I)), i_or_j_ms))

    # Step B: add J on top of the I-bound variant.
    j_ms = _testhelper_move(EnzymeRates.Enumeration._expand_add_dead_end_regulator, with_i,
                            rxn)
    @test !isempty(j_ms)

    # Property: each result has ≥1 J-binding step AND ≥1 I-binding
    # step (adding J must not remove I-binding steps).
    for r in j_ms
        @test !isempty(_testhelper_groups_binding(r, :J))
        @test !isempty(_testhelper_groups_binding(r, :I))
    end
end

@testset "Mechanism — Two regulators competition: 17 variants" begin
    # Bi-bi random with two dead-end inhibitors. Step A: add both regs
    # (9 + 9 = 18 variants). Step B: pick variant with I1 at ≥2 forms
    # and add I2 on top → 17 variants (one less than 18 after dedup
    # accounts for I2's active-form set intersection with existing I1).
    m = @enzyme_mechanism begin
        substrates: A, B
        products: P, Q
        steps: begin
            (E + A ⇌ EA, EB + A ⇌ EAB)
            (E + B ⇌ EB, EA + B ⇌ EAB)
            (E + P ⇌ EP, EQ + P ⇌ EPQ)
            (E + Q ⇌ EQ, EP + Q ⇌ EPQ)
            EAB <--> EPQ
        end
    end
    _testhelper_assert_mechanism_invariants(m)
    rxn = @enzyme_reaction begin
        substrates: A[C], B[N]
        products: P[C], Q[N]
        competitive_inhibitors: I1, I2
    end

    # Step A
    result1 = _testhelper_move(EnzymeRates.Enumeration._expand_add_dead_end_regulator, m,
                               rxn)
    @test length(result1) == 18

    # Pick variant where I1 binds at multiple forms.
    multi = filter(result1) do r
        i1_forms = Set(e.to for e in r.catalytic.steps
                       if EnzymeRates.Enumeration._is_binding(e) &&
                          r.catalytic.metabolites[only(e.taken_up)] === :I1)
        length(i1_forms) >= 2
    end
    @test !isempty(multi)
    m_i1 = first(multi)

    # Step B
    result2 = _testhelper_move(EnzymeRates.Enumeration._expand_add_dead_end_regulator, m_i1,
                               rxn)
    @test length(result2) == 17

    # Property: ≥1 variant has I1 + I2 coexisting on the same enzyme
    # form (non-competing); ≥1 variant has I2 forms that never
    # coexist with I1 (fully-competing). Each form's bound names, read
    # off the variant's graph.
    _testhelper_bound_names(r) =
        [Set(_testhelper_bound(f)) for f in _testhelper_candidate(r).forms]
    has_coexist = any(result2) do r
        any(b -> :I1 in b && :I2 in b, _testhelper_bound_names(r))
    end
    @test has_coexist
    has_compete = any(result2) do r
        forms = _testhelper_bound_names(r)
        has_i2 = any(b -> :I2 in b, forms)
        no_coexist = !any(b -> :I1 in b && :I2 in b, forms)
        has_i2 && no_coexist
    end
    @test has_compete
end

@testset "Mechanism — uni-uni: a reactant as its own inhibitor has no placement" begin
    # The only site that carries neither S nor P is free E, and E·S* has the
    # composition of ES (E·P* that of EP), formed by the S group alone: the
    # gauge exists, and the copy's constant would enter the rate only as
    # K_S + K_S,inh with kcat rescaled. No child is emitted.
    rxn_S = @enzyme_reaction begin
        substrates: S[C]
        products: P[C]
        competitive_inhibitors: S
    end
    rxn_P = @enzyme_reaction begin
        substrates: S[C]
        products: P[C]
        competitive_inhibitors: P
    end
    for rxn in (rxn_S, rxn_P)
        m = @enzyme_mechanism begin
            substrates: S
            products: P
            steps: begin
                E + S ⇌ ES
                ES <--> EP
                E + P ⇌ EP
            end
        end
        @test isempty(_testhelper_move(
            EnzymeRates.Enumeration._expand_add_dead_end_regulator, m, rxn))
    end
end

@testset "Mechanism — ordered bi-bi with A as its own inhibitor: three placements" begin
    # A binds productively at E only, so A* competing with P puts it at {E, EQ},
    # with Q puts it at {E} alone, and competing with B puts it at {EA, EQ} or
    # {E, EA}. E·A* has the composition of EA; EQ·A* and EA·A* are new. The
    # placement at E alone is all-twin, and its gauge exists: EA is formed by the
    # A group alone, so the A group's ratio is 1/ρ and the B group's ρ, each over
    # one step. It is redundant and is not emitted. Mirrors: E + Q ⇌ EQ onto
    # {E, EQ} in the Q group; E + A ⇌ EA onto {E, EA} in the A group.
    rxn = @enzyme_reaction begin
        substrates: A[C], B[N]
        products: P[C], Q[N]
        competitive_inhibitors: A
    end
    m = @enzyme_mechanism begin
        substrates: A, B
        products: P, Q
        steps: begin
            E + A ⇌ EA
            EA + B ⇌ EAB
            EAB <--> EPQ
            EQ + P ⇌ EPQ
            E + Q ⇌ EQ
        end
    end
    at_E_EQ = @enzyme_mechanism begin
        substrates: A, B
        products: P, Q
        steps: begin
            E + A ⇌ EA
            EA + B ⇌ EAB
            EAB <--> EPQ
            EQ + P ⇌ EPQ
            (E + Q ⇌ EQ, EAmod + Q ⇌ EAmodQ)
            (E + A ⇌ EAmod, EQ + A ⇌ EAmodQ)
        end
    end
    at_EA_EQ = @enzyme_mechanism begin
        substrates: A, B
        products: P, Q
        steps: begin
            E + A ⇌ EA
            EA + B ⇌ EAB
            EAB <--> EPQ
            EQ + P ⇌ EPQ
            E + Q ⇌ EQ
            (EA + A ⇌ EAAmod, EQ + A ⇌ EAmodQ)
        end
    end
    at_E_EA = @enzyme_mechanism begin
        substrates: A, B
        products: P, Q
        steps: begin
            (E + A ⇌ EA, EAmod + A ⇌ EAAmod)
            EA + B ⇌ EAB
            EAB <--> EPQ
            EQ + P ⇌ EPQ
            E + Q ⇌ EQ
            (E + A ⇌ EAmod, EA + A ⇌ EAAmod)
        end
    end
    at_E = @enzyme_mechanism begin
        substrates: A, B
        products: P, Q
        steps: begin
            E + A ⇌ EA
            EA + B ⇌ EAB
            EAB <--> EPQ
            EQ + P ⇌ EPQ
            E + Q ⇌ EQ
            E + A ⇌ EAmod
        end
    end
    kids = _testhelper_move(EnzymeRates.Enumeration._expand_add_dead_end_regulator, m, rxn)
    @test length(kids) == 3
    @test Set(kids) == Set([at_E_EQ, at_EA_EQ, at_E_EA])
    # The absent placement adds a constant the data cannot see.
    @test !(at_E in kids)
    # Each child binds A as a copy, in one new kinetic group, and adds exactly one
    # fitted constant; every child has full rank.
    base_fitted = _testhelper_fitted(m)
    for r in kids
        c = _testhelper_candidate(r)
        @test any(k -> _testhelper_binds_copy(c, k) &&
                       r.catalytic.metabolites[only(r.catalytic.steps[k].taken_up)] === :A,
                  first.(_testhelper_group_members(r.catalytic)))
        @test length(_testhelper_group_members(r.catalytic)) ==
              length(_testhelper_group_members(m.catalytic)) + 1
        fitted = _testhelper_fitted(r)
        @test fitted == base_fitted + 1
        @test fitted == _testhelper_identifiable_rank(r)
    end
end

@testset "Mechanism — shared A group with abortive EAQ: the gauge decides" begin
    # A binds at E and, in the same kinetic group, at EQ to form the abortive
    # complex EAQ. A copy of A at E duplicates EA and at EQ duplicates
    # EAQ, so the placements {E} and {E, EQ} are all-twin, and the dwell gauge
    # decides. At {E, EQ}, with E + Q ⇌ EQ mirrored onto the copy forms in the Q
    # group, both twins come from the one A group and share ρ: every step of the A
    # group has ratio 1/ρ and every step of the Q group ratio 1, so the gauge exists
    # and the copy is redundant (not emitted; 6 fitted, rank 5). At {E} the A group
    # holds E + A ⇌ EA, ratio 1/ρ, beside EQ + A ⇌ EAQ, ratio 1: no gauge,
    # since the shared group pins K_A at EQ, where no copy binds. The child is
    # emitted and has full rank (6 fitted, rank 6). {EA, EQ} and {E, EA} create
    # E(A, A*) and are emitted; the second mirrors E + A ⇌ EA in the A group.
    rxn = @enzyme_reaction begin
        substrates: A[C], B[N]
        products: P[C], Q[N]
        competitive_inhibitors: A
    end
    m = @enzyme_mechanism begin
        substrates: A, B
        products: P, Q
        steps: begin
            (E + A ⇌ EA, EQ + A ⇌ EAQ)
            EA + B ⇌ EAB
            EAB <--> EPQ
            EQ + P ⇌ EPQ
            (E + Q ⇌ EQ, EA + Q ⇌ EAQ)
        end
    end
    at_EA_EQ = @enzyme_mechanism begin
        substrates: A, B
        products: P, Q
        steps: begin
            (E + A ⇌ EA, EQ + A ⇌ EAQ)
            EA + B ⇌ EAB
            EAB <--> EPQ
            EQ + P ⇌ EPQ
            (E + Q ⇌ EQ, EA + Q ⇌ EAQ)
            (EA + A ⇌ EAAmod, EQ + A ⇌ EAmodQ)
        end
    end
    at_E_EA = @enzyme_mechanism begin
        substrates: A, B
        products: P, Q
        steps: begin
            (E + A ⇌ EA, EQ + A ⇌ EAQ, EAmod + A ⇌ EAAmod)
            EA + B ⇌ EAB
            EAB <--> EPQ
            EQ + P ⇌ EPQ
            (E + Q ⇌ EQ, EA + Q ⇌ EAQ)
            (E + A ⇌ EAmod, EA + A ⇌ EAAmod)
        end
    end
    at_E = @enzyme_mechanism begin
        substrates: A, B
        products: P, Q
        steps: begin
            (E + A ⇌ EA, EQ + A ⇌ EAQ)
            EA + B ⇌ EAB
            EAB <--> EPQ
            EQ + P ⇌ EPQ
            (E + Q ⇌ EQ, EA + Q ⇌ EAQ)
            E + A ⇌ EAmod
        end
    end
    at_E_EQ = @enzyme_mechanism begin
        substrates: A, B
        products: P, Q
        steps: begin
            (E + A ⇌ EA, EQ + A ⇌ EAQ)
            EA + B ⇌ EAB
            EAB <--> EPQ
            EQ + P ⇌ EPQ
            (E + Q ⇌ EQ, EA + Q ⇌ EAQ, EAmod + Q ⇌ EAmodQ)
            (E + A ⇌ EAmod, EQ + A ⇌ EAmodQ)
        end
    end
    kids = _testhelper_move(EnzymeRates.Enumeration._expand_add_dead_end_regulator, m, rxn)
    @test length(kids) == 3
    @test Set(kids) == Set([at_EA_EQ, at_E_EA, at_E])
    @test !(at_E_EQ in kids)
    @test _testhelper_identifiable_rank(m) == _testhelper_fitted(m)
end

@testset "Mechanism — a second copy may share an older copy's composition" begin
    # Ordered bi-bi with Q as its own inhibitor at {E, EA}, the A binding
    # mirrored onto the copy forms. Competition is decided per site: Q* binds a
    # dead-end site of its own, so E(Q*) carries no productive Q. A binds
    # productively at E and E(Q*), B at EA, P at EQ, Q at E; Q* binds at E and
    # EA. A copy of A goes to {EA, EQ} competing with B and P, to {E, EA}
    # competing with B and Q, to {E, EQ, E(Q*)} competing with A and P, to
    # {E, E(Q*)} competing with A and Q, to {E, EQ} competing with A, P and Q*, to
    # {E, EA, EQ} competing with B, P and Q*, and to {E} competing with A, Q and
    # Q*. E·A* duplicates EA, so {E} is all-twin; but the A group that forms EA
    # also binds at E(Q*), where no A* binds, so it holds a step of ratio 1/ρ beside
    # one of ratio 1, the gauge fails, and {E} is emitted. EA·A* is new, and
    # EQ·A* and E(Q*)·A* have the composition of E(A, Q*), a copy's complex, which
    # is no twin source. Where A* binds at EA, E(A, A*) pins K_A* and the child has
    # full rank. In the child at {EA, EQ} the older Q* group binds at E, a twin of
    # EQ, and at EA, whose complex E(A, Q*) has the composition of E(A*, Q); its
    # constant still shows, since the A·Q term is A·Q/(K_A·K_Q*) + A·Q/(K_Q·K_A*).
    # Where A* binds at E and otherwise, if anywhere, only where its complex has the
    # composition of E(A, Q*), in the children at {E, EQ, E(Q*)}, {E, E(Q*)},
    # {E, EQ} and {E}, no site pins either copy: K_A, K_A*, K_Q and K_Q* enter the
    # law only through the A, Q and A·Q coefficients, so the child carries one
    # phantom. The rule tolerates it: it rejects a copy only on a proof.
    rxn = @enzyme_reaction begin
        substrates: A[C], B[N]
        products: P[C], Q[N]
        competitive_inhibitors: A, Q
    end
    m = @enzyme_mechanism begin
        substrates: A, B
        products: P, Q
        steps: begin
            (E + A ⇌ EA, EQmod + A ⇌ EAQmod)
            EA + B ⇌ EAB
            EAB <--> EPQ
            EQ + P ⇌ EPQ
            E + Q ⇌ EQ
            (E + Q ⇌ EQmod, EA + Q ⇌ EAQmod)
        end
    end
    at_EA_EQ = @enzyme_mechanism begin
        substrates: A, B
        products: P, Q
        steps: begin
            (E + A ⇌ EA, EQmod + A ⇌ EAQmod)
            EA + B ⇌ EAB
            EAB <--> EPQ
            EQ + P ⇌ EPQ
            E + Q ⇌ EQ
            (E + Q ⇌ EQmod, EA + Q ⇌ EAQmod)
            (EA + A ⇌ EAAmod, EQ + A ⇌ EAmodQ)
        end
    end
    at_E_EA = @enzyme_mechanism begin
        substrates: A, B
        products: P, Q
        steps: begin
            (E + A ⇌ EA, EQmod + A ⇌ EAQmod, EAmod + A ⇌ EAAmod)
            EA + B ⇌ EAB
            EAB <--> EPQ
            EQ + P ⇌ EPQ
            E + Q ⇌ EQ
            (E + Q ⇌ EQmod, EA + Q ⇌ EAQmod)
            (E + A ⇌ EAmod, EA + A ⇌ EAAmod)
        end
    end
    at_E_EQ_EQinh = @enzyme_mechanism begin
        substrates: A, B
        products: P, Q
        steps: begin
            (E + A ⇌ EA, EQmod + A ⇌ EAQmod)
            EA + B ⇌ EAB
            EAB <--> EPQ
            EQ + P ⇌ EPQ
            (E + Q ⇌ EQ, EAmod + Q ⇌ EAmodQ)
            (E + Q ⇌ EQmod, EA + Q ⇌ EAQmod,
             EAmod + Q ⇌ EAmodQmod)
            (E + A ⇌ EAmod, EQ + A ⇌ EAmodQ,
             EQmod + A ⇌ EAmodQmod)
        end
    end
    at_E_EQinh = @enzyme_mechanism begin
        substrates: A, B
        products: P, Q
        steps: begin
            (E + A ⇌ EA, EQmod + A ⇌ EAQmod)
            EA + B ⇌ EAB
            EAB <--> EPQ
            EQ + P ⇌ EPQ
            E + Q ⇌ EQ
            (E + Q ⇌ EQmod, EA + Q ⇌ EAQmod,
             EAmod + Q ⇌ EAmodQmod)
            (E + A ⇌ EAmod, EQmod + A ⇌ EAmodQmod)
        end
    end
    at_E_EQ = @enzyme_mechanism begin
        substrates: A, B
        products: P, Q
        steps: begin
            (E + A ⇌ EA, EQmod + A ⇌ EAQmod)
            EA + B ⇌ EAB
            EAB <--> EPQ
            EQ + P ⇌ EPQ
            (E + Q ⇌ EQ, EAmod + Q ⇌ EAmodQ)
            (E + Q ⇌ EQmod, EA + Q ⇌ EAQmod)
            (E + A ⇌ EAmod, EQ + A ⇌ EAmodQ)
        end
    end
    at_E_EA_EQ = @enzyme_mechanism begin
        substrates: A, B
        products: P, Q
        steps: begin
            (E + A ⇌ EA, EQmod + A ⇌ EAQmod, EAmod + A ⇌ EAAmod)
            EA + B ⇌ EAB
            EAB <--> EPQ
            EQ + P ⇌ EPQ
            (E + Q ⇌ EQ, EAmod + Q ⇌ EAmodQ)
            (E + Q ⇌ EQmod, EA + Q ⇌ EAQmod)
            (E + A ⇌ EAmod, EA + A ⇌ EAAmod,
             EQ + A ⇌ EAmodQ)
        end
    end
    at_E = @enzyme_mechanism begin
        substrates: A, B
        products: P, Q
        steps: begin
            (E + A ⇌ EA, EQmod + A ⇌ EAQmod)
            EA + B ⇌ EAB
            EAB <--> EPQ
            EQ + P ⇌ EPQ
            E + Q ⇌ EQ
            (E + Q ⇌ EQmod, EA + Q ⇌ EAQmod)
            E + A ⇌ EAmod
        end
    end
    kids = _testhelper_move(EnzymeRates.Enumeration._expand_add_dead_end_regulator, m, rxn)
    @test length(kids) == 7
    @test Set(kids) == Set([at_EA_EQ, at_E_EA, at_E_EQ_EQinh, at_E_EQinh, at_E_EQ,
                            at_E_EA_EQ, at_E])
    for r in kids
        @test isempty(
            EnzymeRates.Enumeration._redundant_copy_groups(_testhelper_candidate(r)))
        @test EnzymeRates.Enumeration._assert_emission_rules(
            _testhelper_candidate(r)) === nothing
    end
    for k in (at_EA_EQ, at_E_EA, at_E_EA_EQ)
        @test _testhelper_identifiable_rank(k) == _testhelper_fitted(k)
    end
    for k in (at_E_EQ_EQinh, at_E_EQinh, at_E_EQ, at_E)
        @test _testhelper_identifiable_rank(k) == _testhelper_fitted(k) - 1
    end
end

@testset "Mechanism — a copy's form hosts a copy that competes with its reactant" begin
    # Ordered bi-bi with A as its own inhibitor at {E, EQ}, placing a copy of Q.
    # Competition is decided per site: A* binds a dead-end site of its own, so
    # E(A*) carries no productive A and is a site for a Q* that competes with A. A
    # binds productively at E, B at EA, P at EQ, Q at E and E(A*); A* binds at E
    # and EQ. EAB and EPQ hold every substrate or every product and are no
    # site. Q* goes to {E, EQ} competing with A and P, to {E, E(A*)} competing
    # with A and Q, to {EA, EQ} competing with B and P, to {E, EA, E(A*)}
    # competing with B and Q, to {E, EA, EQ} competing with B, P and A*, to
    # {E, EA} competing with B, Q and A*, and to {E} competing with A, Q and A*.
    # E·Q* duplicates EQ, so {E} is all-twin; but the Q group that forms EQ also
    # binds at E(A*), where no Q* binds, so it holds a step of ratio 1/ρ beside one of
    # ratio 1, the gauge fails, and {E} is emitted. EQ·Q* is new, and so are
    # EA·Q* and E(A*)·Q*: E(A, Q*) and E(A*, Q*) have the composition of E(A*, Q),
    # a copy's complex, which is no twin source. Mirrors: E + A ⇌ EA onto E(Q*) in
    # A's group wherever EA is a site, E + Q ⇌ EQ in Q's group wherever EQ is,
    # and E + A* ⇌ E(A*) in A*'s group wherever E(A*) is. Where Q* binds at EQ,
    # E(Q, Q*) pins K_Q*, and {EA, EQ} has no twin site; those children have full
    # rank. Where Q* binds at E and otherwise, if anywhere, only at EA or E(A*), no
    # site pins K_Q* apart from K_Q, so the child carries one phantom, which the rule
    # tolerates: it rejects a copy only on a proof.
    rxn = @enzyme_reaction begin
        substrates: A[C], B[N]
        products: P[C], Q[N]
        competitive_inhibitors: A, Q
    end
    m = @enzyme_mechanism begin
        substrates: A, B
        products: P, Q
        steps: begin
            E + A ⇌ EA
            EA + B ⇌ EAB
            EAB <--> EPQ
            EQ + P ⇌ EPQ
            (E + Q ⇌ EQ, EAmod + Q ⇌ EAmodQ)
            (E + A ⇌ EAmod, EQ + A ⇌ EAmodQ)
        end
    end
    at_E_EQ = @enzyme_mechanism begin
        substrates: A, B
        products: P, Q
        steps: begin
            E + A ⇌ EA
            EA + B ⇌ EAB
            EAB <--> EPQ
            EQ + P ⇌ EPQ
            (E + Q ⇌ EQ, EAmod + Q ⇌ EAmodQ, EQmod + Q ⇌ EQQmod)
            (E + A ⇌ EAmod, EQ + A ⇌ EAmodQ)
            (E + Q ⇌ EQmod, EQ + Q ⇌ EQQmod)
        end
    end
    at_E_EAinh = @enzyme_mechanism begin
        substrates: A, B
        products: P, Q
        steps: begin
            E + A ⇌ EA
            EA + B ⇌ EAB
            EAB <--> EPQ
            EQ + P ⇌ EPQ
            (E + Q ⇌ EQ, EAmod + Q ⇌ EAmodQ)
            (E + A ⇌ EAmod, EQ + A ⇌ EAmodQ,
             EQmod + A ⇌ EAmodQmod)
            (E + Q ⇌ EQmod, EAmod + Q ⇌ EAmodQmod)
        end
    end
    at_EA_EQ = @enzyme_mechanism begin
        substrates: A, B
        products: P, Q
        steps: begin
            E + A ⇌ EA
            EA + B ⇌ EAB
            EAB <--> EPQ
            EQ + P ⇌ EPQ
            (E + Q ⇌ EQ, EAmod + Q ⇌ EAmodQ)
            (E + A ⇌ EAmod, EQ + A ⇌ EAmodQ)
            (EA + Q ⇌ EAQmod, EQ + Q ⇌ EQQmod)
        end
    end
    at_E_EA_EAinh = @enzyme_mechanism begin
        substrates: A, B
        products: P, Q
        steps: begin
            (E + A ⇌ EA, EQmod + A ⇌ EAQmod)
            EA + B ⇌ EAB
            EAB <--> EPQ
            EQ + P ⇌ EPQ
            (E + Q ⇌ EQ, EAmod + Q ⇌ EAmodQ)
            (E + A ⇌ EAmod, EQ + A ⇌ EAmodQ,
             EQmod + A ⇌ EAmodQmod)
            (E + Q ⇌ EQmod, EA + Q ⇌ EAQmod,
             EAmod + Q ⇌ EAmodQmod)
        end
    end
    at_E_EA_EQ = @enzyme_mechanism begin
        substrates: A, B
        products: P, Q
        steps: begin
            (E + A ⇌ EA, EQmod + A ⇌ EAQmod)
            EA + B ⇌ EAB
            EAB <--> EPQ
            EQ + P ⇌ EPQ
            (E + Q ⇌ EQ, EAmod + Q ⇌ EAmodQ, EQmod + Q ⇌ EQQmod)
            (E + A ⇌ EAmod, EQ + A ⇌ EAmodQ)
            (E + Q ⇌ EQmod, EA + Q ⇌ EAQmod,
             EQ + Q ⇌ EQQmod)
        end
    end
    at_E_EA = @enzyme_mechanism begin
        substrates: A, B
        products: P, Q
        steps: begin
            (E + A ⇌ EA, EQmod + A ⇌ EAQmod)
            EA + B ⇌ EAB
            EAB <--> EPQ
            EQ + P ⇌ EPQ
            (E + Q ⇌ EQ, EAmod + Q ⇌ EAmodQ)
            (E + A ⇌ EAmod, EQ + A ⇌ EAmodQ)
            (E + Q ⇌ EQmod, EA + Q ⇌ EAQmod)
        end
    end
    at_E = @enzyme_mechanism begin
        substrates: A, B
        products: P, Q
        steps: begin
            E + A ⇌ EA
            EA + B ⇌ EAB
            EAB <--> EPQ
            EQ + P ⇌ EPQ
            (E + Q ⇌ EQ, EAmod + Q ⇌ EAmodQ)
            (E + A ⇌ EAmod, EQ + A ⇌ EAmodQ)
            E + Q ⇌ EQmod
        end
    end
    kids = _testhelper_move(EnzymeRates.Enumeration._expand_add_dead_end_regulator, m, rxn)
    @test length(kids) == 7
    @test Set(kids) == Set([at_E_EQ, at_E_EAinh, at_EA_EQ, at_E_EA_EAinh, at_E_EA_EQ,
                            at_E_EA, at_E])
    @test at_E_EAinh in kids
    for r in kids
        @test EnzymeRates.Enumeration._assert_emission_rules(
            _testhelper_candidate(r)) === nothing
    end
    for k in (at_E_EQ, at_EA_EQ, at_E_EA_EQ)
        @test _testhelper_identifiable_rank(k) == _testhelper_fitted(k)
    end
    for k in (at_E_EAinh, at_E_EA_EAinh, at_E_EA, at_E)
        @test _testhelper_identifiable_rank(k) == _testhelper_fitted(k) - 1
    end
end

@testset "Mechanism — an inhibitor that competes with A binds the copy's form E(A*)" begin
    # Ordered bi-bi with A as its own inhibitor at {E, EQ}, placing a foreign
    # inhibitor I. A* binds a dead-end site of its own, so E(A*) carries no
    # productive A and is a site for an I that competes with A. A binds productively
    # at E, B at EA, P at EQ, Q at E and E(A*); A* binds at E and EQ. I goes
    # to {E, EQ} competing with A and P, to {E, E(A*)} competing with A and Q, to
    # {EA, EQ} competing with B and P, to {E, EA, E(A*)} competing with B and
    # Q, to {E, EA, EQ} competing with B, P and A*, to {E, EA} competing with
    # B, Q and A*, and to {E} competing with A, Q and A*. I copies no reactant, so
    # every complex it forms is new: all seven children are emitted, each with full
    # rank. Mirrors: E + A ⇌ EA onto EI in A's group wherever EA is a site,
    # E + Q ⇌ EQ in Q's group wherever EQ is, and E + A* ⇌ E(A*) in A*'s group
    # wherever E(A*) is.
    rxn = @enzyme_reaction begin
        substrates: A[C], B[N]
        products: P[C], Q[N]
        competitive_inhibitors: A, I
    end
    m = @enzyme_mechanism begin
        substrates: A, B
        products: P, Q
        steps: begin
            E + A ⇌ EA
            EA + B ⇌ EAB
            EAB <--> EPQ
            EQ + P ⇌ EPQ
            (E + Q ⇌ EQ, EAmod + Q ⇌ EAmodQ)
            (E + A ⇌ EAmod, EQ + A ⇌ EAmodQ)
        end
    end
    at_E_EQ = @enzyme_mechanism begin
        substrates: A, B
        products: P, Q
        modifiers: I
        steps: begin
            E + A ⇌ EA
            EA + B ⇌ EAB
            EAB <--> EPQ
            EQ + P ⇌ EPQ
            (E + Q ⇌ EQ, EAmod + Q ⇌ EAmodQ, EI + Q ⇌ EIQ)
            (E + A ⇌ EAmod, EQ + A ⇌ EAmodQ)
            (E + I ⇌ EI, EQ + I ⇌ EIQ)
        end
    end
    at_E_EAinh = @enzyme_mechanism begin
        substrates: A, B
        products: P, Q
        modifiers: I
        steps: begin
            E + A ⇌ EA
            EA + B ⇌ EAB
            EAB <--> EPQ
            EQ + P ⇌ EPQ
            (E + Q ⇌ EQ, EAmod + Q ⇌ EAmodQ)
            (E + A ⇌ EAmod, EQ + A ⇌ EAmodQ,
             EI + A ⇌ EAmodI)
            (E + I ⇌ EI, EAmod + I ⇌ EAmodI)
        end
    end
    at_EA_EQ = @enzyme_mechanism begin
        substrates: A, B
        products: P, Q
        modifiers: I
        steps: begin
            E + A ⇌ EA
            EA + B ⇌ EAB
            EAB <--> EPQ
            EQ + P ⇌ EPQ
            (E + Q ⇌ EQ, EAmod + Q ⇌ EAmodQ)
            (E + A ⇌ EAmod, EQ + A ⇌ EAmodQ)
            (EA + I ⇌ EAI, EQ + I ⇌ EIQ)
        end
    end
    at_E_EA_EAinh = @enzyme_mechanism begin
        substrates: A, B
        products: P, Q
        modifiers: I
        steps: begin
            (E + A ⇌ EA, EI + A ⇌ EAI)
            EA + B ⇌ EAB
            EAB <--> EPQ
            EQ + P ⇌ EPQ
            (E + Q ⇌ EQ, EAmod + Q ⇌ EAmodQ)
            (E + A ⇌ EAmod, EQ + A ⇌ EAmodQ,
             EI + A ⇌ EAmodI)
            (E + I ⇌ EI, EA + I ⇌ EAI, EAmod + I ⇌ EAmodI)
        end
    end
    at_E_EA_EQ = @enzyme_mechanism begin
        substrates: A, B
        products: P, Q
        modifiers: I
        steps: begin
            (E + A ⇌ EA, EI + A ⇌ EAI)
            EA + B ⇌ EAB
            EAB <--> EPQ
            EQ + P ⇌ EPQ
            (E + Q ⇌ EQ, EAmod + Q ⇌ EAmodQ, EI + Q ⇌ EIQ)
            (E + A ⇌ EAmod, EQ + A ⇌ EAmodQ)
            (E + I ⇌ EI, EA + I ⇌ EAI, EQ + I ⇌ EIQ)
        end
    end
    at_E_EA = @enzyme_mechanism begin
        substrates: A, B
        products: P, Q
        modifiers: I
        steps: begin
            (E + A ⇌ EA, EI + A ⇌ EAI)
            EA + B ⇌ EAB
            EAB <--> EPQ
            EQ + P ⇌ EPQ
            (E + Q ⇌ EQ, EAmod + Q ⇌ EAmodQ)
            (E + A ⇌ EAmod, EQ + A ⇌ EAmodQ)
            (E + I ⇌ EI, EA + I ⇌ EAI)
        end
    end
    at_E = @enzyme_mechanism begin
        substrates: A, B
        products: P, Q
        modifiers: I
        steps: begin
            E + A ⇌ EA
            EA + B ⇌ EAB
            EAB <--> EPQ
            EQ + P ⇌ EPQ
            (E + Q ⇌ EQ, EAmod + Q ⇌ EAmodQ)
            (E + A ⇌ EAmod, EQ + A ⇌ EAmodQ)
            E + I ⇌ EI
        end
    end
    kids = _testhelper_move(EnzymeRates.Enumeration._expand_add_dead_end_regulator, m, rxn)
    @test length(kids) == 7
    @test Set(kids) == Set([at_E_EQ, at_E_EAinh, at_EA_EQ, at_E_EA_EAinh, at_E_EA_EQ,
                            at_E_EA, at_E])
    for k in kids
        @test EnzymeRates.Enumeration._assert_emission_rules(
            _testhelper_candidate(k)) === nothing
        @test _testhelper_fitted(k) == _testhelper_fitted(m) + 1
        @test _testhelper_identifiable_rank(k) == _testhelper_fitted(k)
    end
end

@testset "Mechanism — a form that holds only a copy has room for an inhibitor" begin
    # Uni-bi with S as its own inhibitor at {E, EQ}. The capacity test counts
    # productive bindings only: E(S*) holds S's copy at a dead-end site of its own
    # and no substrate, so it is a site although S is the only substrate. S binds
    # productively at E, P at EQ, Q at E and E(S*); S* binds at E and EQ. A
    # foreign inhibitor I goes to {E, EQ} competing with S and P, to {E, E(S*)}
    # competing with S and Q, and to {E} competing with S, Q and S*. The child at
    # {E, EQ} mirrors E + Q ⇌ EQ onto EI in Q's group, the one at {E, E(S*)}
    # mirrors E + S* ⇌ E(S*) onto EI in S*'s group. Each child has full rank.
    rxn = @enzyme_reaction begin
        substrates: S[CN]
        products: P[C], Q[N]
        competitive_inhibitors: S, I
    end
    m = @enzyme_mechanism begin
        substrates: S
        products: P, Q
        steps: begin
            E + S ⇌ ES
            ES <--> EPQ
            EQ + P ⇌ EPQ
            (E + Q ⇌ EQ, ESmod + Q ⇌ EQSmod)
            (E + S ⇌ ESmod, EQ + S ⇌ EQSmod)
        end
    end
    at_E_EQ = @enzyme_mechanism begin
        substrates: S
        products: P, Q
        modifiers: I
        steps: begin
            E + S ⇌ ES
            ES <--> EPQ
            EQ + P ⇌ EPQ
            (E + Q ⇌ EQ, ESmod + Q ⇌ EQSmod, EI + Q ⇌ EIQ)
            (E + S ⇌ ESmod, EQ + S ⇌ EQSmod)
            (E + I ⇌ EI, EQ + I ⇌ EIQ)
        end
    end
    at_E_ESinh = @enzyme_mechanism begin
        substrates: S
        products: P, Q
        modifiers: I
        steps: begin
            E + S ⇌ ES
            ES <--> EPQ
            EQ + P ⇌ EPQ
            (E + Q ⇌ EQ, ESmod + Q ⇌ EQSmod)
            (E + S ⇌ ESmod, EQ + S ⇌ EQSmod,
             EI + S ⇌ EISmod)
            (E + I ⇌ EI, ESmod + I ⇌ EISmod)
        end
    end
    at_E = @enzyme_mechanism begin
        substrates: S
        products: P, Q
        modifiers: I
        steps: begin
            E + S ⇌ ES
            ES <--> EPQ
            EQ + P ⇌ EPQ
            (E + Q ⇌ EQ, ESmod + Q ⇌ EQSmod)
            (E + S ⇌ ESmod, EQ + S ⇌ EQSmod)
            E + I ⇌ EI
        end
    end
    kids = _testhelper_move(EnzymeRates.Enumeration._expand_add_dead_end_regulator, m, rxn)
    @test length(kids) == 3
    @test Set(kids) == Set([at_E_EQ, at_E_ESinh, at_E])
    for k in kids
        @test EnzymeRates.Enumeration._assert_emission_rules(
            _testhelper_candidate(k)) === nothing
        @test _testhelper_identifiable_rank(k) == _testhelper_fitted(k)
    end
end

@testset "Mechanism — competing with B targets where B binds, not where B* binds" begin
    # Ordered bi-bi with B as its own inhibitor at E (placed competing with A and Q),
    # placing a foreign inhibitor I. Competition with a reactant targets the forms
    # where it binds productively: A at E, B at EA, P at EQ, Q at E; the copy B*
    # binds at E, a target only for competition with B* itself. Capacity leaves out
    # EAB and EPQ; E(B*) is never a target, since nothing binds there. I goes
    # to {E, EQ} competing with A and P, to {E} competing with A and Q, to
    # {EA, EQ} competing with B and P, to {E, EA} competing with B and Q, and to
    # {E, EA, EQ} competing with B, P and B*. In the pattern B and P, free E is a
    # site of B's copy, holds nothing that competes and is not a site of P: targeting
    # by name, B and its copy counting as one competitor, would add it and turn
    # {EA, EQ} into {E, EA, EQ}, so {EA, EQ} comes only from per-site
    # targeting. Mirrors: E + Q ⇌ EQ onto EI in Q's group wherever E and EQ are
    # sites, E + A ⇌ EA in A's group wherever E and EA are. I is foreign, so each
    # child adds one constant and has full rank (7 fitted, rank 7).
    rxn = @enzyme_reaction begin
        substrates: A[C], B[N]
        products: P[C], Q[N]
        competitive_inhibitors: B, I
    end
    m = @enzyme_mechanism begin
        substrates: A, B
        products: P, Q
        steps: begin
            E + A ⇌ EA
            EA + B ⇌ EAB
            EAB <--> EPQ
            EQ + P ⇌ EPQ
            E + Q ⇌ EQ
            E + B ⇌ EBmod
        end
    end
    at_E_EQ = @enzyme_mechanism begin
        substrates: A, B
        products: P, Q
        modifiers: I
        steps: begin
            E + A ⇌ EA
            EA + B ⇌ EAB
            EAB <--> EPQ
            EQ + P ⇌ EPQ
            (E + Q ⇌ EQ, EI + Q ⇌ EIQ)
            E + B ⇌ EBmod
            (E + I ⇌ EI, EQ + I ⇌ EIQ)
        end
    end
    at_E = @enzyme_mechanism begin
        substrates: A, B
        products: P, Q
        modifiers: I
        steps: begin
            E + A ⇌ EA
            EA + B ⇌ EAB
            EAB <--> EPQ
            EQ + P ⇌ EPQ
            E + Q ⇌ EQ
            E + B ⇌ EBmod
            E + I ⇌ EI
        end
    end
    at_EA_EQ = @enzyme_mechanism begin
        substrates: A, B
        products: P, Q
        modifiers: I
        steps: begin
            E + A ⇌ EA
            EA + B ⇌ EAB
            EAB <--> EPQ
            EQ + P ⇌ EPQ
            E + Q ⇌ EQ
            E + B ⇌ EBmod
            (EA + I ⇌ EAI, EQ + I ⇌ EIQ)
        end
    end
    at_E_EA = @enzyme_mechanism begin
        substrates: A, B
        products: P, Q
        modifiers: I
        steps: begin
            (E + A ⇌ EA, EI + A ⇌ EAI)
            EA + B ⇌ EAB
            EAB <--> EPQ
            EQ + P ⇌ EPQ
            E + Q ⇌ EQ
            E + B ⇌ EBmod
            (E + I ⇌ EI, EA + I ⇌ EAI)
        end
    end
    at_E_EA_EQ = @enzyme_mechanism begin
        substrates: A, B
        products: P, Q
        modifiers: I
        steps: begin
            (E + A ⇌ EA, EI + A ⇌ EAI)
            EA + B ⇌ EAB
            EAB <--> EPQ
            EQ + P ⇌ EPQ
            (E + Q ⇌ EQ, EI + Q ⇌ EIQ)
            E + B ⇌ EBmod
            (E + I ⇌ EI, EA + I ⇌ EAI, EQ + I ⇌ EIQ)
        end
    end
    kids = _testhelper_move(EnzymeRates.Enumeration._expand_add_dead_end_regulator, m, rxn)
    @test length(kids) == 5
    @test Set(kids) == Set([at_E_EQ, at_E, at_EA_EQ, at_E_EA, at_E_EA_EQ])
    @test _testhelper_identifiable_rank(m) == _testhelper_fitted(m) == 6
    for k in kids
        @test EnzymeRates.Enumeration._assert_emission_rules(
            _testhelper_candidate(k)) === nothing
        @test _testhelper_identifiable_rank(k) == _testhelper_fitted(k) == 7
    end
end

@testset "Allosteric mechanism — a copy new only in the inactive state is kept" begin
    # Uni-uni whose S binding and chemistry are `:OnlyA`: the inactive conformation
    # binds P but not S. An `:EqualAI` copy of S at E duplicates ES in the active
    # state and is the only S-bound form in the inactive one, where its constant is
    # visible; the placement stands and adds an identifiable constant. So does an
    # `:OnlyI` copy, which binds the inactive state alone. An `:OnlyA` copy binds the
    # active state alone, where it duplicates ES, so it is redundant and dropped. With S
    # binding `:NonequalAI` the copy duplicates ES in each state where it binds, the S
    # group alone forms ES in each, so the gauge exists in each, and all three are dropped.
    rxn = @enzyme_reaction begin
        substrates: S[C]
        products: P[C]
        competitive_inhibitors: S
        oligomeric_state: 2
    end
    onlya = @enzyme_mechanism begin
        substrates: S
        products: P
        oligomeric_state: 2
        steps: begin
            E + S ⇌ ES      :: OnlyA
            ES <--> EP    :: OnlyA
            E + P ⇌ EP      :: EqualAI
        end
    end
    kept = @enzyme_mechanism begin
        substrates: S
        products: P
        oligomeric_state: 2
        steps: begin
            E + S ⇌ ES                :: OnlyA
            ES <--> EP              :: OnlyA
            E + P ⇌ EP                :: EqualAI
            E + S ⇌ ESmod      :: EqualAI
        end
    end
    only_i = @enzyme_mechanism begin
        substrates: S
        products: P
        oligomeric_state: 2
        steps: begin
            E + S ⇌ ES                :: OnlyA
            ES <--> EP              :: OnlyA
            E + P ⇌ EP                :: EqualAI
            E + S ⇌ ESmod      :: OnlyI
        end
    end
    only_a = @enzyme_mechanism begin
        substrates: S
        products: P
        oligomeric_state: 2
        steps: begin
            E + S ⇌ ES                :: OnlyA
            ES <--> EP              :: OnlyA
            E + P ⇌ EP                :: EqualAI
            E + S ⇌ ESmod      :: OnlyA
        end
    end
    kids = _testhelper_move(EnzymeRates.Enumeration._expand_add_dead_end_regulator, onlya,
                            rxn)
    @test length(kids) == 2
    @test Set(kids) == Set([kept, only_i])
    for k in (kept, only_i)
        @test _testhelper_identifiable_rank(k) == _testhelper_identifiable_rank(onlya) + 1
        @test isempty(
            EnzymeRates.Enumeration._redundant_copy_groups(_testhelper_candidate(k)))
    end
    @test !isempty(
        EnzymeRates.Enumeration._redundant_copy_groups(_testhelper_candidate(only_a)))
    @test _testhelper_identifiable_rank(only_a) == _testhelper_identifiable_rank(onlya)
    nonequal = @enzyme_mechanism begin
        substrates: S
        products: P
        oligomeric_state: 2
        steps: begin
            E + S ⇌ ES      :: NonequalAI
            ES <--> EP    :: OnlyA
            E + P ⇌ EP      :: EqualAI
        end
    end
    @test isempty(_testhelper_move(EnzymeRates.Enumeration._expand_add_dead_end_regulator,
                                   nonequal, rxn))
end

@testset "Allosteric mechanism — free enzyme binds a copy in every conformation" begin
    # Uni-uni whose every catalytic group is `:OnlyA`: the inactive conformation has no
    # step, yet it holds its free enzyme. An `:EqualAI` copy of S at E duplicates ES in
    # the active state and is the only S-bound form in the inactive one, where its
    # constant is visible; the placement stands and adds an identifiable constant. So does
    # an `:OnlyI` copy, which binds the inactive state alone. An `:OnlyA` copy binds the
    # active state alone, where it duplicates ES, so it is redundant and dropped.
    rxn = @enzyme_reaction begin
        substrates: S[C]
        products: P[C]
        competitive_inhibitors: S
        oligomeric_state: 2
    end
    dead_inactive = @enzyme_mechanism begin
        substrates: S
        products: P
        oligomeric_state: 2
        steps: begin
            E + S ⇌ ES      :: OnlyA
            ES <--> EP    :: OnlyA
            E + P ⇌ EP      :: OnlyA
        end
    end
    kept = @enzyme_mechanism begin
        substrates: S
        products: P
        oligomeric_state: 2
        steps: begin
            E + S ⇌ ES                :: OnlyA
            ES <--> EP              :: OnlyA
            E + P ⇌ EP                :: OnlyA
            E + S ⇌ ESmod      :: EqualAI
        end
    end
    only_i = @enzyme_mechanism begin
        substrates: S
        products: P
        oligomeric_state: 2
        steps: begin
            E + S ⇌ ES                :: OnlyA
            ES <--> EP              :: OnlyA
            E + P ⇌ EP                :: OnlyA
            E + S ⇌ ESmod      :: OnlyI
        end
    end
    only_a = @enzyme_mechanism begin
        substrates: S
        products: P
        oligomeric_state: 2
        steps: begin
            E + S ⇌ ES                :: OnlyA
            ES <--> EP              :: OnlyA
            E + P ⇌ EP                :: OnlyA
            E + S ⇌ ESmod      :: OnlyA
        end
    end
    kids = _testhelper_move(EnzymeRates.Enumeration._expand_add_dead_end_regulator,
                            dead_inactive, rxn)
    @test length(kids) == 2
    @test Set(kids) == Set([kept, only_i])
    for k in (kept, only_i)
        @test isempty(
            EnzymeRates.Enumeration._redundant_copy_groups(_testhelper_candidate(k)))
        @test _testhelper_identifiable_rank(k) ==
            _testhelper_identifiable_rank(dead_inactive) + 1
    end
    @test !isempty(
        EnzymeRates.Enumeration._redundant_copy_groups(_testhelper_candidate(only_a)))
    @test _testhelper_identifiable_rank(only_a) ==
          _testhelper_identifiable_rank(dead_inactive)
end

@testset "Allosteric mechanism — no relaxation makes a kept copy redundant" begin
    # The copy of S at E is new only in the inactive state, where S does not bind.
    # Relaxing the S binding to `:NonequalAI` brings ES into the inactive state:
    # the copy then duplicates ES in both states, formed in each by the S group
    # alone, so the gauge exists in each, its constant shows only beside K_S, and
    # the child is not emitted (one fitted constant above its rank).
    # The chemistry tagged :NonequalAI while S binds the active state only is refused:
    # the chemistry and the P binding would join E to ES in the inactive state, leaving
    # S's affinity no path to zero (`_FAST_STEP_PINNED`). Tagged with the S binding, it
    # leaves the copy redundant again. The P binding and the copy relax, each alone, and
    # the search stops there. The relaxed copy, with one
    # constant per state, duplicates ES in the active state and is new in the
    # inactive one, so it keeps one phantom, its active-state constant; the copy rule
    # keeps such a `:NonequalAI` copy.
    kept = @enzyme_mechanism begin
        substrates: S
        products: P
        oligomeric_state: 2
        steps: begin
            E + S ⇌ ES                :: OnlyA
            ES <--> EP              :: OnlyA
            E + P ⇌ EP                :: EqualAI
            E + S ⇌ ESmod      :: EqualAI
        end
    end
    relaxed_p = @enzyme_mechanism begin
        substrates: S
        products: P
        oligomeric_state: 2
        steps: begin
            E + S ⇌ ES                :: OnlyA
            ES <--> EP              :: OnlyA
            E + P ⇌ EP                :: NonequalAI
            E + S ⇌ ESmod      :: EqualAI
        end
    end
    relaxed_copy = @enzyme_mechanism begin
        substrates: S
        products: P
        oligomeric_state: 2
        steps: begin
            E + S ⇌ ES                :: OnlyA
            ES <--> EP              :: OnlyA
            E + P ⇌ EP                :: EqualAI
            E + S ⇌ ESmod      :: NonequalAI
        end
    end
    relaxed_s = @enzyme_mechanism begin
        substrates: S
        products: P
        oligomeric_state: 2
        steps: begin
            E + S ⇌ ES                :: NonequalAI
            ES <--> EP              :: OnlyA
            E + P ⇌ EP                :: EqualAI
            E + S ⇌ ESmod      :: EqualAI
        end
    end
    kids = _testhelper_move(EnzymeRates.Enumeration._expand_change_allo_state, kept)
    @test length(kids) == 2
    @test Set(kids) == Set([relaxed_p, relaxed_copy])
    @test !isempty(
        EnzymeRates.Enumeration._redundant_copy_groups(_testhelper_candidate(relaxed_s)))
    @test (_testhelper_fitted(relaxed_s),
           _testhelper_identifiable_rank(relaxed_s)) == (6, 5)
    @test _testhelper_fitted(relaxed_copy) ==
          _testhelper_identifiable_rank(relaxed_copy) + 1
end

@testset "Mechanism — a copy that only matches a form through a conformational isomer" begin
    # E ⇌ E* at rapid equilibrium; S binds E* only, P binds E only. S* competing
    # with P is placed at E (where P binds) and E* (where S binds), with the
    # isomerization mirrored between the copy forms. E*·S* has E*(S)'s composition.
    # E·S* has a composition no form has, but E and E* share offsets, so E·S* and
    # E*(S) have proportional weights: the pattern is all-twin. Both complexes have
    # E*(S) as their twin, the mirrored isomerization joins the two complexes, both at
    # the copy's one factor, so its ratio is 1 like the original's, and only the S
    # binding and the chemistry touch E*(S), each alone in its group: the gauge exists,
    # and the pattern is not emitted. The would-be child has its parent's rank.
    rxn = @enzyme_reaction begin
        substrates: S[C]
        products: P[C]
        competitive_inhibitors: S
    end
    m = @enzyme_mechanism begin
        substrates: S
        products: P
        steps: begin
            E ⇌ Estar
            Estar + S ⇌ EstarS
            EstarS <--> EP
            E + P ⇌ EP
        end
    end
    @test isempty(
        _testhelper_move(EnzymeRates.Enumeration._expand_add_dead_end_regulator, m, rxn))
    absent = @enzyme_mechanism begin
        substrates: S
        products: P
        steps: begin
            (E ⇌ Estar, ESmod ⇌ EstarSmod)
            Estar + S ⇌ EstarS
            EstarS <--> EP
            E + P ⇌ EP
            (E + S ⇌ ESmod, Estar + S ⇌ EstarSmod)
        end
    end
    @test _testhelper_identifiable_rank(absent) == _testhelper_identifiable_rank(m)
end

@testset "Allosteric mechanism — dead-end binding tagged :EqualAI, :OnlyA and :OnlyI" begin
    # Uni-uni allosteric (catalytic_n=2) with mixed regs:
    # :I dead-end, :R allosteric. The move must (a) exclude :R from
    # eligible regs (allosteric ligand), (b) tag the new dead-end group
    # :EqualAI, :OnlyA and :OnlyI, one child each.
    am = @enzyme_mechanism begin
        substrates: S
        products: P
        oligomeric_state: 2
        steps: begin
            E + P ⇌ EP    :: EqualAI
            E + S ⇌ ES    :: EqualAI
            ES <--> EP  :: EqualAI
        end
    end
    _testhelper_assert_mechanism_invariants(am)
    rxn = @enzyme_reaction begin
        substrates: S[C]
        products: P[C]
        competitive_inhibitors: I
        allosteric_regulators: R
        oligomeric_state: 2
    end

    expected = [
        (@enzyme_mechanism begin
             substrates: S ; products: P ; modifiers: I ; oligomeric_state: 2
             steps: begin
                 E + P ⇌ EP   :: EqualAI
                 E + S ⇌ ES   :: EqualAI
                 ES <--> EP :: EqualAI
                 E + I ⇌ EI   :: EqualAI
             end
         end),
        (@enzyme_mechanism begin
             substrates: S ; products: P ; modifiers: I ; oligomeric_state: 2
             steps: begin
                 E + P ⇌ EP   :: EqualAI
                 E + S ⇌ ES   :: EqualAI
                 ES <--> EP :: EqualAI
                 E + I ⇌ EI   :: OnlyA
             end
         end),
        (@enzyme_mechanism begin
             substrates: S ; products: P ; modifiers: I ; oligomeric_state: 2
             steps: begin
                 E + P ⇌ EP   :: EqualAI
                 E + S ⇌ ES   :: EqualAI
                 ES <--> EP :: EqualAI
                 E + I ⇌ EI   :: OnlyI
             end
         end)]

    result = _testhelper_move(EnzymeRates.Enumeration._expand_add_dead_end_regulator, am,
                              rxn)

    # 1. count: 3, one per state of I's group at its one placement, free E, which both
    # conformations hold, so none is dropped (:R is excluded as allosteric).
    @test length(result) == 3
    @test Set(result) == Set(expected)

    # 2. Δ params: +1 (the new group → one K, whatever its state), measured by
    # `parameters`, the exact count of the constants a fit estimates.
    base_fitted = _testhelper_fitted(am)
    for r in result
        r_fitted = _testhelper_fitted(r)
        @test r_fitted == base_fitted + 1
    end

    # 3. kind and invariants
    for r in result
        @test _testhelper_is_allosteric(r)
        _testhelper_assert_mechanism_invariants(r)
    end

    # 4. property: exactly one new catalytic group, with a state of its own.
    for r in result
        @test length(_testhelper_group_members(r.catalytic)) ==
              length(_testhelper_group_members(am.catalytic)) + 1
        @test length(r.allosteric.allosteric_states) ==
              length(am.allosteric.allosteric_states) + 1
    end

    # 5. preservation: oligomeric_state and regulatory_sites
    # unchanged (dead-end add does not touch the regulatory-site
    # vector).
    for r in result
        @test r.allosteric.oligomeric_state == am.allosteric.oligomeric_state
        @test r.allosteric.regulatory_sites == am.allosteric.regulatory_sites
    end
end

@testset "Allosteric mechanism — allosteric-only regulator → empty" begin
    # Rxn declares only :R as an allosteric regulator
    # (no dead-end inhibitors); all declared regulators are
    # allosteric ligands → eligible_regs is empty → result is empty.
    am = @enzyme_mechanism begin
        substrates: S
        products: P
        oligomeric_state: 2
        steps: begin
            E + P ⇌ EP    :: EqualAI
            E + S ⇌ ES    :: EqualAI
            ES <--> EP  :: EqualAI
        end
    end
    @test isempty(_testhelper_move(EnzymeRates.Enumeration._expand_add_dead_end_regulator,
                                   am, uni_uni_allo_reg))
end

@testset "the dead-end move gives its inhibitor three states on an allosteric parent" begin
    # The uni-uni mechanism at two subunits with every group :OnlyA: the inactive
    # conformation holds free E alone. I competes with S and P, which are free only at E,
    # so it has one placement, E. Its new group is :EqualAI, :OnlyA and :OnlyI in turn:
    # each binds E, which both conformations hold, and EI has no productive twin, so none
    # is dropped. Each fits one constant more than the parent.
    rxn = @enzyme_reaction begin
        substrates: S[C]
        products: P[C]
        competitive_inhibitors: I
        oligomeric_state: 2
    end
    parent = @enzyme_mechanism begin
        substrates: S ; products: P ; oligomeric_state: 2
        steps: begin
            E + S ⇌ ES :: OnlyA
            ES <--> EP :: OnlyA
            E + P ⇌ EP :: OnlyA
        end
    end
    expected = [
        (@enzyme_mechanism begin
             substrates: S ; products: P ; modifiers: I ; oligomeric_state: 2
             steps: begin
                 E + S ⇌ ES :: OnlyA
                 ES <--> EP :: OnlyA
                 E + P ⇌ EP :: OnlyA
                 E + I ⇌ EI :: EqualAI
             end
         end),
        (@enzyme_mechanism begin
             substrates: S ; products: P ; modifiers: I ; oligomeric_state: 2
             steps: begin
                 E + S ⇌ ES :: OnlyA
                 ES <--> EP :: OnlyA
                 E + P ⇌ EP :: OnlyA
                 E + I ⇌ EI :: OnlyA
             end
         end),
        (@enzyme_mechanism begin
             substrates: S ; products: P ; modifiers: I ; oligomeric_state: 2
             steps: begin
                 E + S ⇌ ES :: OnlyA
                 ES <--> EP :: OnlyA
                 E + P ⇌ EP :: OnlyA
                 E + I ⇌ EI :: OnlyI
             end
         end)]
    kids = _testhelper_move(EnzymeRates.Enumeration._expand_add_dead_end_regulator, parent,
                            rxn)
    @test length(kids) == 3
    @test Set(kids) == Set(expected)
    @test kids == expected
    @test all(k -> _testhelper_fitted(k) == _testhelper_fitted(parent) + 1, expected)
end

@testset "the dead-end move drops an :OnlyI inhibitor the inactive state cannot bind" begin
    # The three-step ping-pong with both half-reactions :OnlyA: the inactive conformation
    # runs Q's binding alone and never reaches the covalent form E_res_plusA_minusP. I has
    # three placements: free E, where A and Q are free; the covalent form, where B and P
    # are free; and both, the first half-reaction mirrored between EI and
    # EI_res_plusA_minusP in its own group. At the covalent form alone an :OnlyI group
    # binds in neither conformation and adds no constant, so the move keeps the
    # :EqualAI and the :OnlyA child there; the other two placements give all three. Every
    # child it keeps fits one constant more than the parent.
    rxn = @enzyme_reaction begin
        substrates: A[CX], B[N]
        products: P[C], Q[NX]
        competitive_inhibitors: I
    end
    parent = @enzyme_mechanism begin
        substrates: A, B ; products: P, Q ; oligomeric_state: 1
        steps: begin
            E + Q <--> EQ                     :: EqualAI
            E + A <--> E_res_plusA_minusP + P :: OnlyA
            E_res_plusA_minusP + B <--> EQ    :: OnlyA
        end
    end
    # Each child's comment names its placement and the state of I's group.
    # Untyped: inside `T[...]` the parser reads `begin` as an index keyword.
    expected = [
        # E, :EqualAI
        (@enzyme_mechanism begin
             substrates: A, B ; products: P, Q ; modifiers: I ; oligomeric_state: 1
             steps: begin
                 E + Q <--> EQ                     :: EqualAI
                 E + A <--> E_res_plusA_minusP + P :: OnlyA
                 E_res_plusA_minusP + B <--> EQ    :: OnlyA
                 E + I ⇌ EI                        :: EqualAI
             end
         end),
        # E, :OnlyA
        (@enzyme_mechanism begin
             substrates: A, B ; products: P, Q ; modifiers: I ; oligomeric_state: 1
             steps: begin
                 E + Q <--> EQ                     :: EqualAI
                 E + A <--> E_res_plusA_minusP + P :: OnlyA
                 E_res_plusA_minusP + B <--> EQ    :: OnlyA
                 E + I ⇌ EI                        :: OnlyA
             end
         end),
        # E, :OnlyI
        (@enzyme_mechanism begin
             substrates: A, B ; products: P, Q ; modifiers: I ; oligomeric_state: 1
             steps: begin
                 E + Q <--> EQ                     :: EqualAI
                 E + A <--> E_res_plusA_minusP + P :: OnlyA
                 E_res_plusA_minusP + B <--> EQ    :: OnlyA
                 E + I ⇌ EI                        :: OnlyI
             end
         end),
        # the covalent form, :EqualAI
        (@enzyme_mechanism begin
             substrates: A, B ; products: P, Q ; modifiers: I ; oligomeric_state: 1
             steps: begin
                 E + Q <--> EQ                                :: EqualAI
                 E + A <--> E_res_plusA_minusP + P            :: OnlyA
                 E_res_plusA_minusP + B <--> EQ               :: OnlyA
                 E_res_plusA_minusP + I ⇌ EI_res_plusA_minusP :: EqualAI
             end
         end),
        # the covalent form, :OnlyA
        (@enzyme_mechanism begin
             substrates: A, B ; products: P, Q ; modifiers: I ; oligomeric_state: 1
             steps: begin
                 E + Q <--> EQ                                :: EqualAI
                 E + A <--> E_res_plusA_minusP + P            :: OnlyA
                 E_res_plusA_minusP + B <--> EQ               :: OnlyA
                 E_res_plusA_minusP + I ⇌ EI_res_plusA_minusP :: OnlyA
             end
         end),
        # both, :EqualAI
        (@enzyme_mechanism begin
             substrates: A, B ; products: P, Q ; modifiers: I ; oligomeric_state: 1
             steps: begin
                 E + Q <--> EQ                                         :: EqualAI
                 (E + A <--> E_res_plusA_minusP + P,
                  EI + A <--> EI_res_plusA_minusP + P)                 :: OnlyA
                 E_res_plusA_minusP + B <--> EQ                        :: OnlyA
                 (E + I ⇌ EI, E_res_plusA_minusP + I ⇌ EI_res_plusA_minusP) :: EqualAI
             end
         end),
        # both, :OnlyA
        (@enzyme_mechanism begin
             substrates: A, B ; products: P, Q ; modifiers: I ; oligomeric_state: 1
             steps: begin
                 E + Q <--> EQ                                         :: EqualAI
                 (E + A <--> E_res_plusA_minusP + P,
                  EI + A <--> EI_res_plusA_minusP + P)                 :: OnlyA
                 E_res_plusA_minusP + B <--> EQ                        :: OnlyA
                 (E + I ⇌ EI, E_res_plusA_minusP + I ⇌ EI_res_plusA_minusP) :: OnlyA
             end
         end),
        # both, :OnlyI
        (@enzyme_mechanism begin
             substrates: A, B ; products: P, Q ; modifiers: I ; oligomeric_state: 1
             steps: begin
                 E + Q <--> EQ                                         :: EqualAI
                 (E + A <--> E_res_plusA_minusP + P,
                  EI + A <--> EI_res_plusA_minusP + P)                 :: OnlyA
                 E_res_plusA_minusP + B <--> EQ                        :: OnlyA
                 (E + I ⇌ EI, E_res_plusA_minusP + I ⇌ EI_res_plusA_minusP) :: OnlyI
             end
         end)]
    nowhere = @enzyme_mechanism begin
        substrates: A, B ; products: P, Q ; modifiers: I ; oligomeric_state: 1
        steps: begin
            E + Q <--> EQ                                :: EqualAI
            E + A <--> E_res_plusA_minusP + P            :: OnlyA
            E_res_plusA_minusP + B <--> EQ               :: OnlyA
            E_res_plusA_minusP + I ⇌ EI_res_plusA_minusP :: OnlyI
        end
    end
    kids = _testhelper_move(EnzymeRates.Enumeration._expand_add_dead_end_regulator, parent,
                            rxn)
    @test length(kids) == 8
    @test Set(kids) == Set(expected)
    @test all(k -> _testhelper_fitted(k) == _testhelper_fitted(parent) + 1, expected)
    @test !EnzymeRates.Enumeration._inhibitors_bind(_testhelper_candidate(nowhere))
    @test _testhelper_fitted(nowhere) == _testhelper_fitted(parent)
end

@testset "an inhibitor added before or after to_allosteric gives the same mechanisms" begin
    # The uni-uni mechanism, I bound at free E, its one placement, and the chemistry alone
    # tagged :OnlyA among the mechanism's groups. Binding I and then tagging gives three
    # children with that :OnlyA set, I's group :EqualAI, :OnlyA and :OnlyI; tagging and
    # then binding I gives the same three.
    rxn = @enzyme_reaction begin
        substrates: S[C]
        products: P[C]
        competitive_inhibitors: I
    end
    plain = @enzyme_mechanism begin
        substrates: S ; products: P
        steps: begin
            E + S ⇌ ES
            ES <--> EP
            E + P ⇌ EP
        end
    end
    with_i = @enzyme_mechanism begin
        substrates: S ; products: P ; modifiers: I
        steps: begin
            E + S ⇌ ES
            ES <--> EP
            E + P ⇌ EP
            E + I ⇌ EI
        end
    end
    chemistry_only_a = @enzyme_mechanism begin
        substrates: S ; products: P ; oligomeric_state: 1
        steps: begin
            E + S ⇌ ES :: EqualAI
            ES <--> EP :: OnlyA
            E + P ⇌ EP :: EqualAI
        end
    end
    expected = [
        (@enzyme_mechanism begin
             substrates: S ; products: P ; modifiers: I ; oligomeric_state: 1
             steps: begin
                 E + S ⇌ ES :: EqualAI
                 ES <--> EP :: OnlyA
                 E + P ⇌ EP :: EqualAI
                 E + I ⇌ EI :: EqualAI
             end
         end),
        (@enzyme_mechanism begin
             substrates: S ; products: P ; modifiers: I ; oligomeric_state: 1
             steps: begin
                 E + S ⇌ ES :: EqualAI
                 ES <--> EP :: OnlyA
                 E + P ⇌ EP :: EqualAI
                 E + I ⇌ EI :: OnlyA
             end
         end),
        (@enzyme_mechanism begin
             substrates: S ; products: P ; modifiers: I ; oligomeric_state: 1
             steps: begin
                 E + S ⇌ ES :: EqualAI
                 ES <--> EP :: OnlyA
                 E + P ⇌ EP :: EqualAI
                 E + I ⇌ EI :: OnlyI
             end
         end)]
    dead_end = EnzymeRates.Enumeration._expand_add_dead_end_regulator
    to_allosteric = EnzymeRates.Enumeration._expand_to_allosteric
    @test _testhelper_move(dead_end, plain, rxn) == [with_i]
    # The children that tag the chemistry alone :OnlyA among S's binding, P's and the
    # chemistry; each keeps the catalytic scheme of `with_i`.
    group(text) = with_i.catalytic.steps[_testhelper_edge(with_i, text)].kinetic_group
    catalytic = group.(["E + S ⇌ ES", "ES <--> EP", "E + P ⇌ EP"])
    bound_then_tagged = filter(_testhelper_move(to_allosteric, with_i, rxn)) do c
        c.allosteric.allosteric_states[catalytic] == [:EqualAI, :OnlyA, :EqualAI]
    end
    @test length(bound_then_tagged) == 3
    @test Set(bound_then_tagged) == Set(expected)
    @test chemistry_only_a in _testhelper_move(to_allosteric, plain, rxn)
    tagged_then_bound = _testhelper_move(dead_end, chemistry_only_a, rxn)
    @test length(tagged_then_bound) == 3
    @test Set(tagged_then_bound) == Set(expected)
end

@testset "the dead-end move adds a second inhibitor beside an :OnlyI one" begin
    # The uni-uni mechanism with I bound at free E and tagged :OnlyI, every other group
    # :EqualAI: the active conformation lacks EI. J competes with S and P, and with I or
    # not; each is free only at E, so J has one placement, E. Its group takes the three
    # states, each binds E in both conformations, and EJ has no productive twin, so none is
    # dropped. Each fits one constant more than the parent.
    rxn = @enzyme_reaction begin
        substrates: S[C]
        products: P[C]
        competitive_inhibitors: I, J
    end
    parent = @enzyme_mechanism begin
        substrates: S ; products: P ; modifiers: I ; oligomeric_state: 1
        steps: begin
            E + I ⇌ EI :: OnlyI
            E + P ⇌ EP :: EqualAI
            E + S ⇌ ES :: EqualAI
            ES <--> EP :: EqualAI
        end
    end
    expected = [
        (@enzyme_mechanism begin
             substrates: S ; products: P ; modifiers: I, J ; oligomeric_state: 1
             steps: begin
                 E + I ⇌ EI :: OnlyI
                 E + P ⇌ EP :: EqualAI
                 E + S ⇌ ES :: EqualAI
                 ES <--> EP :: EqualAI
                 E + J ⇌ EJ :: EqualAI
             end
         end),
        (@enzyme_mechanism begin
             substrates: S ; products: P ; modifiers: I, J ; oligomeric_state: 1
             steps: begin
                 E + I ⇌ EI :: OnlyI
                 E + P ⇌ EP :: EqualAI
                 E + S ⇌ ES :: EqualAI
                 ES <--> EP :: EqualAI
                 E + J ⇌ EJ :: OnlyA
             end
         end),
        (@enzyme_mechanism begin
             substrates: S ; products: P ; modifiers: I, J ; oligomeric_state: 1
             steps: begin
                 E + I ⇌ EI :: OnlyI
                 E + P ⇌ EP :: EqualAI
                 E + S ⇌ ES :: EqualAI
                 ES <--> EP :: EqualAI
                 E + J ⇌ EJ :: OnlyI
             end
         end)]
    kids = _testhelper_move(EnzymeRates.Enumeration._expand_add_dead_end_regulator, parent,
                            rxn)
    @test length(kids) == 3
    @test Set(kids) == Set(expected)
    @test all(k -> _testhelper_fitted(k) == _testhelper_fitted(parent) + 1, expected)
end

@testset "Mechanism — no regulators: empty (negative)" begin
    m = @enzyme_mechanism begin
        substrates: S
        products: P
        steps: begin
            E + S ⇌ ES
            ES <--> EP
            EP ⇌ E + P
        end
    end
    rxn = @enzyme_reaction begin
        substrates: S[C]
        products: P[C]
    end
    @test isempty(
        _testhelper_move(EnzymeRates.Enumeration._expand_add_dead_end_regulator, m, rxn))
end

end

# The tagging move would drop an allosteric parent's states.
@testset "the tagging move refuses an allosteric mechanism" begin
    rxn = @enzyme_reaction begin
        substrates: S[C]
        products: P[C]
        allosteric_regulators: R
    end
    allosteric = @enzyme_mechanism begin
        substrates: S
        products: P
        oligomeric_state: 2
        steps: begin
            E + S ⇌ ES     :: EqualAI
            ES <--> EP   :: OnlyA
            E + P ⇌ EP     :: EqualAI
        end
        regulatory_site(2): R::OnlyI
    end
    @test_throws(ArgumentError(
        "_expand_to_allosteric: the mechanism is allosteric; the move makes a plain " *
        "one allosteric"),
        EnzymeRates.Enumeration._expand_to_allosteric(_testhelper_candidate(allosteric),
                                                      rxn))
end

# ─── _expand_to_allosteric ─────────────────────────────────────────────
@testset "_expand_to_allosteric" begin

@testset "Mechanism — bi-bi sequential: every :OnlyA set, and two pairs :NonequalAI" begin
    # Five kinetic groups: the bindings of A, B, Q and P and the chemistry EAB <--> EPQ.
    # At two subunits with no regulator the move emits the 2^5 − 1 = 31 nonempty sets of
    # groups tagged :OnlyA, each adding L alone. Tagged :NonequalAI, seven minimal sets
    # gain a constant beyond L: the chemistry alone and each pair of bindings. Tagged
    # :OnlyA, two of the pairs are refused: A with B leaves an inactive conformation that
    # binds no substrate and keeps the chemistry, and Q with P one that binds no product
    # (`_BRIDGE`). The move emits those two tagged :NonequalAI, each adding two
    # constants.
    m = @enzyme_mechanism begin
        substrates: A, B
        products: P, Q
        steps: begin
            E + A ⇌ EA
            EA + B ⇌ EAB
            E + Q ⇌ EQ
            EQ + P ⇌ EPQ
            EAB <--> EPQ
        end
    end
    _testhelper_assert_mechanism_invariants(m)
    bi_bi_allo_rxn = @enzyme_reaction begin
        substrates: A[C], B[N]
        products: P[C], Q[N]
        oligomeric_state: 2
    end
    a_b_nonequal = @enzyme_mechanism begin
        substrates: A, B ; products: P, Q ; oligomeric_state: 2
        steps: begin
            E + A ⇌ EA   :: NonequalAI
            EA + B ⇌ EAB :: NonequalAI
            E + Q ⇌ EQ   :: EqualAI
            EQ + P ⇌ EPQ :: EqualAI
            EAB <--> EPQ :: EqualAI
        end
    end
    q_p_nonequal = @enzyme_mechanism begin
        substrates: A, B ; products: P, Q ; oligomeric_state: 2
        steps: begin
            E + A ⇌ EA   :: EqualAI
            EA + B ⇌ EAB :: EqualAI
            E + Q ⇌ EQ   :: NonequalAI
            EQ + P ⇌ EPQ :: NonequalAI
            EAB <--> EPQ :: EqualAI
        end
    end
    a_b_only_a = @_testhelper_unguarded @enzyme_mechanism begin
        substrates: A, B ; products: P, Q ; oligomeric_state: 2
        steps: begin
            E + A ⇌ EA   :: OnlyA
            EA + B ⇌ EAB :: OnlyA
            E + Q ⇌ EQ   :: EqualAI
            EQ + P ⇌ EPQ :: EqualAI
            EAB <--> EPQ :: EqualAI
        end
    end
    q_p_only_a = @_testhelper_unguarded @enzyme_mechanism begin
        substrates: A, B ; products: P, Q ; oligomeric_state: 2
        steps: begin
            E + A ⇌ EA   :: EqualAI
            EA + B ⇌ EAB :: EqualAI
            E + Q ⇌ EQ   :: OnlyA
            EQ + P ⇌ EPQ :: OnlyA
            EAB <--> EPQ :: EqualAI
        end
    end
    result = _testhelper_move(EnzymeRates.Enumeration._expand_to_allosteric, m,
                              bi_bi_allo_rxn)
    @test length(result) == 31 + 2
    @test allunique(result)
    nonequal = filter(r -> :NonequalAI in r.allosteric.allosteric_states, result)
    @test Set(nonequal) == Set([a_b_nonequal, q_p_nonequal])
    only_a = setdiff(result, nonequal)
    @test all(r -> issubset(r.allosteric.allosteric_states, (:OnlyA, :EqualAI)) &&
                   isempty(r.allosteric.regulatory_sites), only_a)
    for (refused, step) in ((a_b_only_a, "EA + B ⇌ EAB"), (q_p_only_a, "EQ + P ⇌ EPQ"))
        @test refused in only_a
        v = EnzymeRates.Derivation._guard(refused)
        @test (v.fault, EnzymeRates.Mechanisms.step_text(m.catalytic, v.edge)) ==
              (EnzymeRates.Derivation._BRIDGE, step)
    end
    fitted = EnzymeRates.Derivation.fitted_parameter_count
    @test all(r -> fitted(r) == fitted(m) + 1, only_a)
    @test all(r -> fitted(r) == fitted(m) + 2, nonequal)
    for r in result
        _testhelper_assert_mechanism_invariants(r)
        @test !all(==(:EqualAI), r.allosteric.allosteric_states)
    end
end

@testset "Mechanism — bi-bi ping-pong at two subunits: every :OnlyA set" begin
    # SEED: bi-bi ping-pong, 6 kinetic groups: 4 bindings (A and Q at free E, B and P at
    # the covalent form E_res_plusA_minusP) and the two half-reactions EA <-->
    # EP_res_plusA_minusP and EB_res_plusA_minusP ⇌ EQ. At two subunits with no
    # regulator the move emits the 63 nonempty sets of groups tagged :OnlyA. Tagged
    # :NonequalAI, six minimal sets gain a constant beyond L; four of them are refused
    # tagged :OnlyA (Q with the second half, A with B, Q with B, B with the second half),
    # and the move emits those four tagged :NonequalAI. candidate_mechanisms returns 14 of
    # the :OnlyA sets and the four :NonequalAI children.
    m = @enzyme_mechanism begin
        substrates: A, B
        products: P, Q
        steps: begin
            E + A ⇌ EA
            E_res_plusA_minusP + B ⇌ EB_res_plusA_minusP
            E + Q ⇌ EQ
            E_res_plusA_minusP + P ⇌ EP_res_plusA_minusP
            EA <--> EP_res_plusA_minusP
            EB_res_plusA_minusP ⇌ EQ
        end
    end
    _testhelper_assert_mechanism_invariants(m)
    bi_bi_pp_allo_rxn = @enzyme_reaction begin
        substrates: A[CX], B[N]
        products: P[C], Q[NX]
        oligomeric_state: 2
    end
    q_half2_nonequal = @enzyme_mechanism begin
        substrates: A, B ; products: P, Q ; oligomeric_state: 2
        steps: begin
            E + A ⇌ EA                                   :: EqualAI
            E_res_plusA_minusP + B ⇌ EB_res_plusA_minusP :: EqualAI
            E + Q ⇌ EQ                                   :: NonequalAI
            E_res_plusA_minusP + P ⇌ EP_res_plusA_minusP :: EqualAI
            EA <--> EP_res_plusA_minusP                  :: EqualAI
            EB_res_plusA_minusP ⇌ EQ                     :: NonequalAI
        end
    end
    a_b_nonequal = @enzyme_mechanism begin
        substrates: A, B ; products: P, Q ; oligomeric_state: 2
        steps: begin
            E + A ⇌ EA                                   :: NonequalAI
            E_res_plusA_minusP + B ⇌ EB_res_plusA_minusP :: NonequalAI
            E + Q ⇌ EQ                                   :: EqualAI
            E_res_plusA_minusP + P ⇌ EP_res_plusA_minusP :: EqualAI
            EA <--> EP_res_plusA_minusP                  :: EqualAI
            EB_res_plusA_minusP ⇌ EQ                     :: EqualAI
        end
    end
    q_b_nonequal = @enzyme_mechanism begin
        substrates: A, B ; products: P, Q ; oligomeric_state: 2
        steps: begin
            E + A ⇌ EA                                   :: EqualAI
            E_res_plusA_minusP + B ⇌ EB_res_plusA_minusP :: NonequalAI
            E + Q ⇌ EQ                                   :: NonequalAI
            E_res_plusA_minusP + P ⇌ EP_res_plusA_minusP :: EqualAI
            EA <--> EP_res_plusA_minusP                  :: EqualAI
            EB_res_plusA_minusP ⇌ EQ                     :: EqualAI
        end
    end
    b_half2_nonequal = @enzyme_mechanism begin
        substrates: A, B ; products: P, Q ; oligomeric_state: 2
        steps: begin
            E + A ⇌ EA                                   :: EqualAI
            E_res_plusA_minusP + B ⇌ EB_res_plusA_minusP :: NonequalAI
            E + Q ⇌ EQ                                   :: EqualAI
            E_res_plusA_minusP + P ⇌ EP_res_plusA_minusP :: EqualAI
            EA <--> EP_res_plusA_minusP                  :: EqualAI
            EB_res_plusA_minusP ⇌ EQ                     :: NonequalAI
        end
    end
    # The first half-reaction :OnlyA, alone or with A's or P's binding: class 3 rejects
    # each on the faces with every substrate zero and with B zero, a split whose limit
    # the speeds do not move but which differs from the law (limit test `differ`).
    a_half1 = @_testhelper_unguarded @enzyme_mechanism begin
        substrates: A, B ; products: P, Q ; oligomeric_state: 2
        steps: begin
            E + A ⇌ EA                                   :: OnlyA
            E_res_plusA_minusP + B ⇌ EB_res_plusA_minusP :: EqualAI
            E + Q ⇌ EQ                                   :: EqualAI
            E_res_plusA_minusP + P ⇌ EP_res_plusA_minusP :: EqualAI
            EA <--> EP_res_plusA_minusP                  :: OnlyA
            EB_res_plusA_minusP ⇌ EQ                     :: EqualAI
        end
    end
    p_half1 = @_testhelper_unguarded @enzyme_mechanism begin
        substrates: A, B ; products: P, Q ; oligomeric_state: 2
        steps: begin
            E + A ⇌ EA                                   :: EqualAI
            E_res_plusA_minusP + B ⇌ EB_res_plusA_minusP :: EqualAI
            E + Q ⇌ EQ                                   :: EqualAI
            E_res_plusA_minusP + P ⇌ EP_res_plusA_minusP :: OnlyA
            EA <--> EP_res_plusA_minusP                  :: OnlyA
            EB_res_plusA_minusP ⇌ EQ                     :: EqualAI
        end
    end
    ap_half1 = @_testhelper_unguarded @enzyme_mechanism begin
        substrates: A, B ; products: P, Q ; oligomeric_state: 2
        steps: begin
            E + A ⇌ EA                                   :: OnlyA
            E_res_plusA_minusP + B ⇌ EB_res_plusA_minusP :: EqualAI
            E + Q ⇌ EQ                                   :: EqualAI
            E_res_plusA_minusP + P ⇌ EP_res_plusA_minusP :: OnlyA
            EA <--> EP_res_plusA_minusP                  :: OnlyA
            EB_res_plusA_minusP ⇌ EQ                     :: EqualAI
        end
    end
    # Both half-reactions :OnlyA beside the bindings of A, of Q, of A and Q, of A and P,
    # and of A, Q and P: each is returned.
    a_both = @enzyme_mechanism begin
        substrates: A, B ; products: P, Q ; oligomeric_state: 2
        steps: begin
            E + A ⇌ EA                                   :: OnlyA
            E_res_plusA_minusP + B ⇌ EB_res_plusA_minusP :: EqualAI
            E + Q ⇌ EQ                                   :: EqualAI
            E_res_plusA_minusP + P ⇌ EP_res_plusA_minusP :: EqualAI
            EA <--> EP_res_plusA_minusP                  :: OnlyA
            EB_res_plusA_minusP ⇌ EQ                     :: OnlyA
        end
    end
    q_both = @enzyme_mechanism begin
        substrates: A, B ; products: P, Q ; oligomeric_state: 2
        steps: begin
            E + A ⇌ EA                                   :: EqualAI
            E_res_plusA_minusP + B ⇌ EB_res_plusA_minusP :: EqualAI
            E + Q ⇌ EQ                                   :: OnlyA
            E_res_plusA_minusP + P ⇌ EP_res_plusA_minusP :: EqualAI
            EA <--> EP_res_plusA_minusP                  :: OnlyA
            EB_res_plusA_minusP ⇌ EQ                     :: OnlyA
        end
    end
    aq_both = @enzyme_mechanism begin
        substrates: A, B ; products: P, Q ; oligomeric_state: 2
        steps: begin
            E + A ⇌ EA                                   :: OnlyA
            E_res_plusA_minusP + B ⇌ EB_res_plusA_minusP :: EqualAI
            E + Q ⇌ EQ                                   :: OnlyA
            E_res_plusA_minusP + P ⇌ EP_res_plusA_minusP :: EqualAI
            EA <--> EP_res_plusA_minusP                  :: OnlyA
            EB_res_plusA_minusP ⇌ EQ                     :: OnlyA
        end
    end
    ap_both = @enzyme_mechanism begin
        substrates: A, B ; products: P, Q ; oligomeric_state: 2
        steps: begin
            E + A ⇌ EA                                   :: OnlyA
            E_res_plusA_minusP + B ⇌ EB_res_plusA_minusP :: EqualAI
            E + Q ⇌ EQ                                   :: EqualAI
            E_res_plusA_minusP + P ⇌ EP_res_plusA_minusP :: OnlyA
            EA <--> EP_res_plusA_minusP                  :: OnlyA
            EB_res_plusA_minusP ⇌ EQ                     :: OnlyA
        end
    end
    apq_both = @enzyme_mechanism begin
        substrates: A, B ; products: P, Q ; oligomeric_state: 2
        steps: begin
            E + A ⇌ EA                                   :: OnlyA
            E_res_plusA_minusP + B ⇌ EB_res_plusA_minusP :: EqualAI
            E + Q ⇌ EQ                                   :: OnlyA
            E_res_plusA_minusP + P ⇌ EP_res_plusA_minusP :: OnlyA
            EA <--> EP_res_plusA_minusP                  :: OnlyA
            EB_res_plusA_minusP ⇌ EQ                     :: OnlyA
        end
    end
    kept = [a_both, q_both, aq_both, ap_both, apq_both]
    result = _testhelper_move(EnzymeRates.Enumeration._expand_to_allosteric, m,
                              bi_bi_pp_allo_rxn)
    @test length(result) == 63 + 4
    @test allunique(result)
    nonequal = filter(r -> :NonequalAI in r.allosteric.allosteric_states, result)
    @test Set(nonequal) == Set([q_half2_nonequal, a_b_nonequal, q_b_nonequal,
                                b_half2_nonequal])
    # The other 63, distinct and each with only :OnlyA and :EqualAI groups, are the 63
    # nonempty sets tagged :OnlyA.
    @test all(setdiff(result, nonequal)) do r
        r.catalytic == m.catalytic && r.allosteric.oligomeric_state == 2 &&
            isempty(r.allosteric.regulatory_sites) &&
            :OnlyA in r.allosteric.allosteric_states &&
            issubset(r.allosteric.allosteric_states, (:OnlyA, :EqualAI))
    end
    @test issubset([a_half1, p_half1, ap_half1, kept...], result)
    for c in (a_half1, p_half1, ap_half1)
        v = EnzymeRates.Derivation._guard(c)
        names =
            EnzymeRates.Derivation._face_names(
                m.catalytic, EnzymeRates.Derivation._faces(m.catalytic)[v.face])
        @test (v.fault, v.limit, names) == (EnzymeRates.Derivation._SPLIT,
                                            EnzymeRates.Derivation._LIMIT_DIFFER, [:A, :B])
    end
    returned = [c for c in _testhelper_children(bi_bi_pp_allo_rxn, [m])
                if _testhelper_is_allosteric(c)]
    @test length(returned) == 14 + 4
    @test issubset([kept; nonequal], returned)
    @test isdisjoint([a_half1, p_half1, ap_half1], returned)
    # Δ params: each kept :OnlyA child adds L only.
    @test _testhelper_param_deltas(m, kept) == fill(1, 5)
    foreach(_testhelper_assert_mechanism_invariants, result)
end

@testset "Mechanism — merged abortive ping-pong: hyperbolic laws above one subunit" begin
    # A ping-pong whose A binding is fused into the product-bound covalent form
    # EP_res_plusA_minusP and whose B binding is fused into EQ, with the abortive complex
    # EAQ: A binds EQ in A's group, and Q binds EP_res_plusA_minusP in Q's group. Its
    # scheme is hyperbolic, so two subunits stay open. The move tags each nonempty set of
    # the four groups :OnlyA: fifteen children, none :NonequalAI (each minimal set that
    # gains is returned tagged :OnlyA). At two subunits candidate_mechanisms returns a
    # child only when both conformations' schemes are hyperbolic (`_enumerator_returns`).
    # A's group alone :OnlyA leaves the inactive conformation Q's abortive binding, and
    # its scheme fails `_hyperbolic_catalysis`; likewise Q's group alone, which leaves
    # A's. The guard rejects both as well (`_FAST_STEP_PINNED`), so the testset pins that
    # the rule refuses them and that candidate_mechanisms does not return them.
    m = @enzyme_mechanism begin
        substrates: A, B
        products: P, Q
        steps: begin
            (E + A ⇌ EP_res_plusA_minusP, EQ + A ⇌ EAQ)
            (E + Q ⇌ EQ, EP_res_plusA_minusP + Q ⇌ EAQ)
            E_res_plusA_minusP + P <--> EP_res_plusA_minusP
            E_res_plusA_minusP + B <--> EQ
        end
    end
    rxn = @enzyme_reaction begin
        substrates: A[CX], B[N]
        products: P[C], Q[NX]
        oligomeric_state: 2
    end
    a = @_testhelper_unguarded @enzyme_mechanism begin
        substrates: A, B ; products: P, Q ; oligomeric_state: 2
        steps: begin
            (E + A ⇌ EP_res_plusA_minusP, EQ + A ⇌ EAQ)     :: OnlyA
            (E + Q ⇌ EQ, EP_res_plusA_minusP + Q ⇌ EAQ)     :: EqualAI
            E_res_plusA_minusP + P <--> EP_res_plusA_minusP :: EqualAI
            E_res_plusA_minusP + B <--> EQ                  :: EqualAI
        end
    end
    q = @_testhelper_unguarded @enzyme_mechanism begin
        substrates: A, B ; products: P, Q ; oligomeric_state: 2
        steps: begin
            (E + A ⇌ EP_res_plusA_minusP, EQ + A ⇌ EAQ)     :: EqualAI
            (E + Q ⇌ EQ, EP_res_plusA_minusP + Q ⇌ EAQ)     :: OnlyA
            E_res_plusA_minusP + P <--> EP_res_plusA_minusP :: EqualAI
            E_res_plusA_minusP + B <--> EQ                  :: EqualAI
        end
    end
    a_q = @_testhelper_unguarded @enzyme_mechanism begin
        substrates: A, B ; products: P, Q ; oligomeric_state: 2
        steps: begin
            (E + A ⇌ EP_res_plusA_minusP, EQ + A ⇌ EAQ)     :: OnlyA
            (E + Q ⇌ EQ, EP_res_plusA_minusP + Q ⇌ EAQ)     :: OnlyA
            E_res_plusA_minusP + P <--> EP_res_plusA_minusP :: EqualAI
            E_res_plusA_minusP + B <--> EQ                  :: EqualAI
        end
    end
    p = @_testhelper_unguarded @enzyme_mechanism begin
        substrates: A, B ; products: P, Q ; oligomeric_state: 2
        steps: begin
            (E + A ⇌ EP_res_plusA_minusP, EQ + A ⇌ EAQ)     :: EqualAI
            (E + Q ⇌ EQ, EP_res_plusA_minusP + Q ⇌ EAQ)     :: EqualAI
            E_res_plusA_minusP + P <--> EP_res_plusA_minusP :: OnlyA
            E_res_plusA_minusP + B <--> EQ                  :: EqualAI
        end
    end
    a_p = @_testhelper_unguarded @enzyme_mechanism begin
        substrates: A, B ; products: P, Q ; oligomeric_state: 2
        steps: begin
            (E + A ⇌ EP_res_plusA_minusP, EQ + A ⇌ EAQ)     :: OnlyA
            (E + Q ⇌ EQ, EP_res_plusA_minusP + Q ⇌ EAQ)     :: EqualAI
            E_res_plusA_minusP + P <--> EP_res_plusA_minusP :: OnlyA
            E_res_plusA_minusP + B <--> EQ                  :: EqualAI
        end
    end
    q_p = @_testhelper_unguarded @enzyme_mechanism begin
        substrates: A, B ; products: P, Q ; oligomeric_state: 2
        steps: begin
            (E + A ⇌ EP_res_plusA_minusP, EQ + A ⇌ EAQ)     :: EqualAI
            (E + Q ⇌ EQ, EP_res_plusA_minusP + Q ⇌ EAQ)     :: OnlyA
            E_res_plusA_minusP + P <--> EP_res_plusA_minusP :: OnlyA
            E_res_plusA_minusP + B <--> EQ                  :: EqualAI
        end
    end
    a_q_p = @_testhelper_unguarded @enzyme_mechanism begin
        substrates: A, B ; products: P, Q ; oligomeric_state: 2
        steps: begin
            (E + A ⇌ EP_res_plusA_minusP, EQ + A ⇌ EAQ)     :: OnlyA
            (E + Q ⇌ EQ, EP_res_plusA_minusP + Q ⇌ EAQ)     :: OnlyA
            E_res_plusA_minusP + P <--> EP_res_plusA_minusP :: OnlyA
            E_res_plusA_minusP + B <--> EQ                  :: EqualAI
        end
    end
    b = @_testhelper_unguarded @enzyme_mechanism begin
        substrates: A, B ; products: P, Q ; oligomeric_state: 2
        steps: begin
            (E + A ⇌ EP_res_plusA_minusP, EQ + A ⇌ EAQ)     :: EqualAI
            (E + Q ⇌ EQ, EP_res_plusA_minusP + Q ⇌ EAQ)     :: EqualAI
            E_res_plusA_minusP + P <--> EP_res_plusA_minusP :: EqualAI
            E_res_plusA_minusP + B <--> EQ                  :: OnlyA
        end
    end
    a_b = @_testhelper_unguarded @enzyme_mechanism begin
        substrates: A, B ; products: P, Q ; oligomeric_state: 2
        steps: begin
            (E + A ⇌ EP_res_plusA_minusP, EQ + A ⇌ EAQ)     :: OnlyA
            (E + Q ⇌ EQ, EP_res_plusA_minusP + Q ⇌ EAQ)     :: EqualAI
            E_res_plusA_minusP + P <--> EP_res_plusA_minusP :: EqualAI
            E_res_plusA_minusP + B <--> EQ                  :: OnlyA
        end
    end
    q_b = @_testhelper_unguarded @enzyme_mechanism begin
        substrates: A, B ; products: P, Q ; oligomeric_state: 2
        steps: begin
            (E + A ⇌ EP_res_plusA_minusP, EQ + A ⇌ EAQ)     :: EqualAI
            (E + Q ⇌ EQ, EP_res_plusA_minusP + Q ⇌ EAQ)     :: OnlyA
            E_res_plusA_minusP + P <--> EP_res_plusA_minusP :: EqualAI
            E_res_plusA_minusP + B <--> EQ                  :: OnlyA
        end
    end
    a_q_b = @_testhelper_unguarded @enzyme_mechanism begin
        substrates: A, B ; products: P, Q ; oligomeric_state: 2
        steps: begin
            (E + A ⇌ EP_res_plusA_minusP, EQ + A ⇌ EAQ)     :: OnlyA
            (E + Q ⇌ EQ, EP_res_plusA_minusP + Q ⇌ EAQ)     :: OnlyA
            E_res_plusA_minusP + P <--> EP_res_plusA_minusP :: EqualAI
            E_res_plusA_minusP + B <--> EQ                  :: OnlyA
        end
    end
    p_b = @_testhelper_unguarded @enzyme_mechanism begin
        substrates: A, B ; products: P, Q ; oligomeric_state: 2
        steps: begin
            (E + A ⇌ EP_res_plusA_minusP, EQ + A ⇌ EAQ)     :: EqualAI
            (E + Q ⇌ EQ, EP_res_plusA_minusP + Q ⇌ EAQ)     :: EqualAI
            E_res_plusA_minusP + P <--> EP_res_plusA_minusP :: OnlyA
            E_res_plusA_minusP + B <--> EQ                  :: OnlyA
        end
    end
    a_p_b = @_testhelper_unguarded @enzyme_mechanism begin
        substrates: A, B ; products: P, Q ; oligomeric_state: 2
        steps: begin
            (E + A ⇌ EP_res_plusA_minusP, EQ + A ⇌ EAQ)     :: OnlyA
            (E + Q ⇌ EQ, EP_res_plusA_minusP + Q ⇌ EAQ)     :: EqualAI
            E_res_plusA_minusP + P <--> EP_res_plusA_minusP :: OnlyA
            E_res_plusA_minusP + B <--> EQ                  :: OnlyA
        end
    end
    q_p_b = @_testhelper_unguarded @enzyme_mechanism begin
        substrates: A, B ; products: P, Q ; oligomeric_state: 2
        steps: begin
            (E + A ⇌ EP_res_plusA_minusP, EQ + A ⇌ EAQ)     :: EqualAI
            (E + Q ⇌ EQ, EP_res_plusA_minusP + Q ⇌ EAQ)     :: OnlyA
            E_res_plusA_minusP + P <--> EP_res_plusA_minusP :: OnlyA
            E_res_plusA_minusP + B <--> EQ                  :: OnlyA
        end
    end
    a_q_p_b = @_testhelper_unguarded @enzyme_mechanism begin
        substrates: A, B ; products: P, Q ; oligomeric_state: 2
        steps: begin
            (E + A ⇌ EP_res_plusA_minusP, EQ + A ⇌ EAQ)     :: OnlyA
            (E + Q ⇌ EQ, EP_res_plusA_minusP + Q ⇌ EAQ)     :: OnlyA
            E_res_plusA_minusP + P <--> EP_res_plusA_minusP :: OnlyA
            E_res_plusA_minusP + B <--> EQ                  :: OnlyA
        end
    end
    @test EnzymeRates.Enumeration._hyperbolic_catalysis(_testhelper_candidate(m))
    kids = _testhelper_move(EnzymeRates.Enumeration._expand_to_allosteric, m, rxn)
    @test length(kids) == 15
    @test Set(kids) == Set([a, q, a_q, p, a_p, q_p, a_q_p, b, a_b, q_b, a_q_b, p_b, a_p_b,
                            q_p_b, a_q_p_b])
    returned = _testhelper_children(rxn, [m])
    for (dropped, step) in ((a, "E + A ⇌ EP_res_plusA_minusP"), (q, "E + Q ⇌ EQ"))
        @test !EnzymeRates.Enumeration._hyperbolic_where_required(
            _testhelper_candidate(dropped))
        runs = EnzymeRates.Mechanisms.steps_running_in_allosteric_conformations(dropped)
        @test EnzymeRates.Enumeration._hyperbolic_catalysis(
            _testhelper_candidate(m), runs.active)
        @test !EnzymeRates.Enumeration._hyperbolic_catalysis(
            _testhelper_candidate(m), runs.inactive)
        v = EnzymeRates.Derivation._guard(dropped)
        @test (v.fault, EnzymeRates.Mechanisms.step_text(m.catalytic, v.edge)) ==
              (EnzymeRates.Derivation._FAST_STEP_PINNED, step)
        @test dropped ∉ returned
    end
    hyperbolic = setdiff(kids, [a, q])
    @test all(EnzymeRates.Enumeration._hyperbolic_where_required,
              _testhelper_candidate.(hyperbolic))
    @test Set(c for c in returned if _testhelper_is_allosteric(c)) == Set(hyperbolic)
end

@testset "Mechanism — a tagging child whose law is not hyperbolic is not returned" begin
    # A uni-bi parent: S binds at rapid equilibrium, the chemistry ES → EPQ is at steady
    # state, and the products leave in either order at rapid equilibrium, each product's
    # two bindings in one group. Its scheme is hyperbolic, so two subunits stay open, and
    # at two subunits candidate_mechanisms returns a child only when both conformations'
    # schemes are hyperbolic (`_enumerator_returns`). The move tags each nonempty set of
    # the four groups :OnlyA, and tags P's and Q's groups :NonequalAI together, a minimal
    # gaining set refused as :OnlyA (an inactive conformation that binds no product and
    # keeps the chemistry, `_BRIDGE`): sixteen children. Two fail the rule: with P's group
    # alone :OnlyA the inactive conformation keeps Q's two bindings, and its scheme fails
    # `_hyperbolic_catalysis`; likewise Q's group alone. The guard rejects both as well,
    # so the testset pins that the rule refuses them and that candidate_mechanisms does
    # not return them.
    m = @enzyme_mechanism begin
        substrates: S
        products: P, Q
        steps: begin
            E + S ⇌ ES
            ES <--> EPQ
            (E + P ⇌ EP, EQ + P ⇌ EPQ)
            (E + Q ⇌ EQ, EP + Q ⇌ EPQ)
        end
    end
    rxn = @enzyme_reaction begin
        substrates: S[CN]
        products: P[C], Q[N]
        oligomeric_state: 2
    end
    s = @_testhelper_unguarded @enzyme_mechanism begin
        substrates: S ; products: P, Q ; oligomeric_state: 2
        steps: begin
            E + S ⇌ ES                 :: OnlyA
            ES <--> EPQ                :: EqualAI
            (E + P ⇌ EP, EQ + P ⇌ EPQ) :: EqualAI
            (E + Q ⇌ EQ, EP + Q ⇌ EPQ) :: EqualAI
        end
    end
    chem = @_testhelper_unguarded @enzyme_mechanism begin
        substrates: S ; products: P, Q ; oligomeric_state: 2
        steps: begin
            E + S ⇌ ES                 :: EqualAI
            ES <--> EPQ                :: OnlyA
            (E + P ⇌ EP, EQ + P ⇌ EPQ) :: EqualAI
            (E + Q ⇌ EQ, EP + Q ⇌ EPQ) :: EqualAI
        end
    end
    s_chem = @_testhelper_unguarded @enzyme_mechanism begin
        substrates: S ; products: P, Q ; oligomeric_state: 2
        steps: begin
            E + S ⇌ ES                 :: OnlyA
            ES <--> EPQ                :: OnlyA
            (E + P ⇌ EP, EQ + P ⇌ EPQ) :: EqualAI
            (E + Q ⇌ EQ, EP + Q ⇌ EPQ) :: EqualAI
        end
    end
    p = @_testhelper_unguarded @enzyme_mechanism begin
        substrates: S ; products: P, Q ; oligomeric_state: 2
        steps: begin
            E + S ⇌ ES                 :: EqualAI
            ES <--> EPQ                :: EqualAI
            (E + P ⇌ EP, EQ + P ⇌ EPQ) :: OnlyA
            (E + Q ⇌ EQ, EP + Q ⇌ EPQ) :: EqualAI
        end
    end
    s_p = @_testhelper_unguarded @enzyme_mechanism begin
        substrates: S ; products: P, Q ; oligomeric_state: 2
        steps: begin
            E + S ⇌ ES                 :: OnlyA
            ES <--> EPQ                :: EqualAI
            (E + P ⇌ EP, EQ + P ⇌ EPQ) :: OnlyA
            (E + Q ⇌ EQ, EP + Q ⇌ EPQ) :: EqualAI
        end
    end
    chem_p = @_testhelper_unguarded @enzyme_mechanism begin
        substrates: S ; products: P, Q ; oligomeric_state: 2
        steps: begin
            E + S ⇌ ES                 :: EqualAI
            ES <--> EPQ                :: OnlyA
            (E + P ⇌ EP, EQ + P ⇌ EPQ) :: OnlyA
            (E + Q ⇌ EQ, EP + Q ⇌ EPQ) :: EqualAI
        end
    end
    s_chem_p = @_testhelper_unguarded @enzyme_mechanism begin
        substrates: S ; products: P, Q ; oligomeric_state: 2
        steps: begin
            E + S ⇌ ES                 :: OnlyA
            ES <--> EPQ                :: OnlyA
            (E + P ⇌ EP, EQ + P ⇌ EPQ) :: OnlyA
            (E + Q ⇌ EQ, EP + Q ⇌ EPQ) :: EqualAI
        end
    end
    q = @_testhelper_unguarded @enzyme_mechanism begin
        substrates: S ; products: P, Q ; oligomeric_state: 2
        steps: begin
            E + S ⇌ ES                 :: EqualAI
            ES <--> EPQ                :: EqualAI
            (E + P ⇌ EP, EQ + P ⇌ EPQ) :: EqualAI
            (E + Q ⇌ EQ, EP + Q ⇌ EPQ) :: OnlyA
        end
    end
    s_q = @_testhelper_unguarded @enzyme_mechanism begin
        substrates: S ; products: P, Q ; oligomeric_state: 2
        steps: begin
            E + S ⇌ ES                 :: OnlyA
            ES <--> EPQ                :: EqualAI
            (E + P ⇌ EP, EQ + P ⇌ EPQ) :: EqualAI
            (E + Q ⇌ EQ, EP + Q ⇌ EPQ) :: OnlyA
        end
    end
    chem_q = @_testhelper_unguarded @enzyme_mechanism begin
        substrates: S ; products: P, Q ; oligomeric_state: 2
        steps: begin
            E + S ⇌ ES                 :: EqualAI
            ES <--> EPQ                :: OnlyA
            (E + P ⇌ EP, EQ + P ⇌ EPQ) :: EqualAI
            (E + Q ⇌ EQ, EP + Q ⇌ EPQ) :: OnlyA
        end
    end
    s_chem_q = @_testhelper_unguarded @enzyme_mechanism begin
        substrates: S ; products: P, Q ; oligomeric_state: 2
        steps: begin
            E + S ⇌ ES                 :: OnlyA
            ES <--> EPQ                :: OnlyA
            (E + P ⇌ EP, EQ + P ⇌ EPQ) :: EqualAI
            (E + Q ⇌ EQ, EP + Q ⇌ EPQ) :: OnlyA
        end
    end
    p_q = @_testhelper_unguarded @enzyme_mechanism begin
        substrates: S ; products: P, Q ; oligomeric_state: 2
        steps: begin
            E + S ⇌ ES                 :: EqualAI
            ES <--> EPQ                :: EqualAI
            (E + P ⇌ EP, EQ + P ⇌ EPQ) :: OnlyA
            (E + Q ⇌ EQ, EP + Q ⇌ EPQ) :: OnlyA
        end
    end
    s_p_q = @_testhelper_unguarded @enzyme_mechanism begin
        substrates: S ; products: P, Q ; oligomeric_state: 2
        steps: begin
            E + S ⇌ ES                 :: OnlyA
            ES <--> EPQ                :: EqualAI
            (E + P ⇌ EP, EQ + P ⇌ EPQ) :: OnlyA
            (E + Q ⇌ EQ, EP + Q ⇌ EPQ) :: OnlyA
        end
    end
    chem_p_q = @_testhelper_unguarded @enzyme_mechanism begin
        substrates: S ; products: P, Q ; oligomeric_state: 2
        steps: begin
            E + S ⇌ ES                 :: EqualAI
            ES <--> EPQ                :: OnlyA
            (E + P ⇌ EP, EQ + P ⇌ EPQ) :: OnlyA
            (E + Q ⇌ EQ, EP + Q ⇌ EPQ) :: OnlyA
        end
    end
    s_chem_p_q = @_testhelper_unguarded @enzyme_mechanism begin
        substrates: S ; products: P, Q ; oligomeric_state: 2
        steps: begin
            E + S ⇌ ES                 :: OnlyA
            ES <--> EPQ                :: OnlyA
            (E + P ⇌ EP, EQ + P ⇌ EPQ) :: OnlyA
            (E + Q ⇌ EQ, EP + Q ⇌ EPQ) :: OnlyA
        end
    end
    p_q_nonequal = @enzyme_mechanism begin
        substrates: S ; products: P, Q ; oligomeric_state: 2
        steps: begin
            E + S ⇌ ES                 :: EqualAI
            ES <--> EPQ                :: EqualAI
            (E + P ⇌ EP, EQ + P ⇌ EPQ) :: NonequalAI
            (E + Q ⇌ EQ, EP + Q ⇌ EPQ) :: NonequalAI
        end
    end
    @test EnzymeRates.Enumeration._hyperbolic_catalysis(_testhelper_candidate(m))
    kids = _testhelper_move(EnzymeRates.Enumeration._expand_to_allosteric, m, rxn)
    @test length(kids) == 16
    @test Set(kids) == Set([s, chem, s_chem, p, s_p, chem_p, s_chem_p, q, s_q, chem_q,
                            s_chem_q, p_q, s_p_q, chem_p_q, s_chem_p_q, p_q_nonequal])
    v = EnzymeRates.Derivation._guard(p_q)
    @test (v.fault, EnzymeRates.Mechanisms.step_text(m.catalytic, v.edge)) ==
          (EnzymeRates.Derivation._BRIDGE, "EQ + P ⇌ EPQ")
    returned = [c for c in _testhelper_children(rxn, [m]) if _testhelper_is_allosteric(c)]
    for (dropped, step) in ((p, "E + P ⇌ EP"), (q, "E + Q ⇌ EQ"))
        @test !EnzymeRates.Enumeration._hyperbolic_where_required(
            _testhelper_candidate(dropped))
        @test dropped ∉ returned
        runs = EnzymeRates.Mechanisms.steps_running_in_allosteric_conformations(dropped)
        @test EnzymeRates.Enumeration._hyperbolic_catalysis(
            _testhelper_candidate(m), runs.active)
        @test !EnzymeRates.Enumeration._hyperbolic_catalysis(
            _testhelper_candidate(m), runs.inactive)
        v = EnzymeRates.Derivation._guard(dropped)
        @test (v.fault, EnzymeRates.Mechanisms.step_text(m.catalytic, v.edge)) ==
              (EnzymeRates.Derivation._FAST_STEP_PINNED, step)
    end
    # S's binding alone :OnlyA and both products :OnlyA are refused by the guard alone;
    # every other child is returned.
    @test Set(returned) == Set(setdiff(kids, [p, q, s, p_q]))
    @test all(EnzymeRates.Enumeration._hyperbolic_where_required,
              _testhelper_candidate.(returned))
end

@testset "to_allosteric emits every :OnlyA set of a uni-uni mechanism" begin
    # Three kinetic groups: S's binding, the chemistry and P's binding. At one subunit
    # with no regulator the move tags each nonempty set of them :OnlyA and the rest
    # :EqualAI: seven children. Tagged :NonequalAI, one binding alone gains no constant
    # beyond L: the inactive conformation runs the whole cycle, so it must give the same
    # Keq, and thermodynamics gives it the active conformation's constant. The two
    # bindings together gain one, and so does the chemistry alone. Both sets are returned
    # tagged :OnlyA, so no child is tagged :NonequalAI.
    rxn = @enzyme_reaction begin
        substrates: S[C]
        products: P[C]
    end
    m = @enzyme_mechanism begin
        substrates: S
        products: P
        steps: begin
            E + S ⇌ ES
            ES <--> EP
            E + P ⇌ EP
        end
    end
    only_s = @_testhelper_unguarded @enzyme_mechanism begin
        substrates: S ; products: P ; oligomeric_state: 1
        steps: begin
            E + S ⇌ ES :: OnlyA
            ES <--> EP :: EqualAI
            E + P ⇌ EP :: EqualAI
        end
    end
    only_chem = @_testhelper_unguarded @enzyme_mechanism begin
        substrates: S ; products: P ; oligomeric_state: 1
        steps: begin
            E + S ⇌ ES :: EqualAI
            ES <--> EP :: OnlyA
            E + P ⇌ EP :: EqualAI
        end
    end
    s_chem = @_testhelper_unguarded @enzyme_mechanism begin
        substrates: S ; products: P ; oligomeric_state: 1
        steps: begin
            E + S ⇌ ES :: OnlyA
            ES <--> EP :: OnlyA
            E + P ⇌ EP :: EqualAI
        end
    end
    only_p = @_testhelper_unguarded @enzyme_mechanism begin
        substrates: S ; products: P ; oligomeric_state: 1
        steps: begin
            E + S ⇌ ES :: EqualAI
            ES <--> EP :: EqualAI
            E + P ⇌ EP :: OnlyA
        end
    end
    s_p = @_testhelper_unguarded @enzyme_mechanism begin
        substrates: S ; products: P ; oligomeric_state: 1
        steps: begin
            E + S ⇌ ES :: OnlyA
            ES <--> EP :: EqualAI
            E + P ⇌ EP :: OnlyA
        end
    end
    chem_p = @_testhelper_unguarded @enzyme_mechanism begin
        substrates: S ; products: P ; oligomeric_state: 1
        steps: begin
            E + S ⇌ ES :: EqualAI
            ES <--> EP :: OnlyA
            E + P ⇌ EP :: OnlyA
        end
    end
    every = @_testhelper_unguarded @enzyme_mechanism begin
        substrates: S ; products: P ; oligomeric_state: 1
        steps: begin
            E + S ⇌ ES :: OnlyA
            ES <--> EP :: OnlyA
            E + P ⇌ EP :: OnlyA
        end
    end
    s_nonequal = @enzyme_mechanism begin
        substrates: S ; products: P ; oligomeric_state: 1
        steps: begin
            E + S ⇌ ES :: NonequalAI
            ES <--> EP :: EqualAI
            E + P ⇌ EP :: EqualAI
        end
    end
    s_p_nonequal = @enzyme_mechanism begin
        substrates: S ; products: P ; oligomeric_state: 1
        steps: begin
            E + S ⇌ ES :: NonequalAI
            ES <--> EP :: EqualAI
            E + P ⇌ EP :: NonequalAI
        end
    end
    chem_nonequal = @enzyme_mechanism begin
        substrates: S ; products: P ; oligomeric_state: 1
        steps: begin
            E + S ⇌ ES :: EqualAI
            ES <--> EP :: NonequalAI
            E + P ⇌ EP :: EqualAI
        end
    end
    kids = _testhelper_move(EnzymeRates.Enumeration._expand_to_allosteric, m, rxn)
    @test length(kids) == 7
    @test Set(kids) == Set([only_s, only_chem, s_chem, only_p, s_p, chem_p, every])
    fitted = EnzymeRates.Derivation.fitted_parameter_count
    @test fitted(s_nonequal) == fitted(m) + 1
    @test fitted(s_p_nonequal) == fitted(chem_nonequal) == fitted(m) + 2
    @test all(EnzymeRates.Enumeration._enumerator_returns,
              _testhelper_candidate.([s_p, only_chem]))
    foreach(_testhelper_assert_mechanism_invariants, kids)
end

@testset "to_allosteric on a three-step ping-pong" begin
    # The three-step ping-pong of census pingpong seed 8fe71fded09a27ec: free E takes up
    # A and gives off P in one Theorell–Chance step, leaving the covalent form E(;
    # residual = A - P), which takes up B into EQ (fused); Q's binding is the one binding
    # group. At one subunit with no regulator the move tags each nonempty set of the
    # three groups :OnlyA: seven children. Tagged :NonequalAI, each group alone gains a
    # constant beyond L, and each is returned tagged :OnlyA, so no child is tagged
    # :NonequalAI. candidate_mechanisms returns all seven.
    m = @enzyme_mechanism begin
        substrates: A, B
        products: P, Q
        steps: begin
            E + Q <--> EQ
            E + A <--> E_res_plusA_minusP + P
            E_res_plusA_minusP + B <--> EQ
        end
    end
    rxn = @enzyme_reaction begin
        substrates: A[CX], B[N]
        products: P[C], Q[NX]
    end
    only_q = @_testhelper_unguarded @enzyme_mechanism begin
        substrates: A, B ; products: P, Q ; oligomeric_state: 1
        steps: begin
            E + Q <--> EQ                     :: OnlyA
            E + A <--> E_res_plusA_minusP + P :: EqualAI
            E_res_plusA_minusP + B <--> EQ    :: EqualAI
        end
    end
    only_half1 = @_testhelper_unguarded @enzyme_mechanism begin
        substrates: A, B ; products: P, Q ; oligomeric_state: 1
        steps: begin
            E + Q <--> EQ                     :: EqualAI
            E + A <--> E_res_plusA_minusP + P :: OnlyA
            E_res_plusA_minusP + B <--> EQ    :: EqualAI
        end
    end
    q_half1 = @_testhelper_unguarded @enzyme_mechanism begin
        substrates: A, B ; products: P, Q ; oligomeric_state: 1
        steps: begin
            E + Q <--> EQ                     :: OnlyA
            E + A <--> E_res_plusA_minusP + P :: OnlyA
            E_res_plusA_minusP + B <--> EQ    :: EqualAI
        end
    end
    only_half2 = @_testhelper_unguarded @enzyme_mechanism begin
        substrates: A, B ; products: P, Q ; oligomeric_state: 1
        steps: begin
            E + Q <--> EQ                     :: EqualAI
            E + A <--> E_res_plusA_minusP + P :: EqualAI
            E_res_plusA_minusP + B <--> EQ    :: OnlyA
        end
    end
    q_half2 = @_testhelper_unguarded @enzyme_mechanism begin
        substrates: A, B ; products: P, Q ; oligomeric_state: 1
        steps: begin
            E + Q <--> EQ                     :: OnlyA
            E + A <--> E_res_plusA_minusP + P :: EqualAI
            E_res_plusA_minusP + B <--> EQ    :: OnlyA
        end
    end
    both_halves = @_testhelper_unguarded @enzyme_mechanism begin
        substrates: A, B ; products: P, Q ; oligomeric_state: 1
        steps: begin
            E + Q <--> EQ                     :: EqualAI
            E + A <--> E_res_plusA_minusP + P :: OnlyA
            E_res_plusA_minusP + B <--> EQ    :: OnlyA
        end
    end
    every = @_testhelper_unguarded @enzyme_mechanism begin
        substrates: A, B ; products: P, Q ; oligomeric_state: 1
        steps: begin
            E + Q <--> EQ                     :: OnlyA
            E + A <--> E_res_plusA_minusP + P :: OnlyA
            E_res_plusA_minusP + B <--> EQ    :: OnlyA
        end
    end
    expected = [only_q, only_half1, q_half1, only_half2, q_half2, both_halves, every]
    kids = _testhelper_move(EnzymeRates.Enumeration._expand_to_allosteric, m, rxn)
    @test length(kids) == 7
    @test Set(kids) == Set(expected)
    @test Set(c for c in _testhelper_children(rxn, [m])
              if _testhelper_is_allosteric(c)) == Set(expected)
    # Every group :OnlyA cuts EQ and the covalent form off in the inactive conformation;
    # both oracles agree with the law on every face.
    @test all(_testhelper_admission_agreement(
        every; rng = Random.MersenneTwister(15209), corners = true))
    _testhelper_check_law_against_limit(every)
end

@testset "to_allosteric pairs a chemistry-only :OnlyA set with each regulator" begin
    # The uni-uni mechanism at two subunits, with the allosteric regulator R declared.
    # Beside the seven :OnlyA sets, the set that holds the chemistry alone comes twice
    # more, with R at a new site of multiplicity 2 tagged :OnlyA and :OnlyI: a version
    # that stops only chemistry shows nothing in a substrate or product titration that
    # the plain mechanism lacks. A set with a binding group gets no regulator, and no
    # child is tagged :NonequalAI (see the uni-uni testset above).
    m = @enzyme_mechanism begin
        substrates: S
        products: P
        steps: begin
            E + S ⇌ ES
            ES <--> EP
            E + P ⇌ EP
        end
    end
    only_s = @_testhelper_unguarded @enzyme_mechanism begin
        substrates: S ; products: P ; oligomeric_state: 2
        steps: begin
            E + S ⇌ ES :: OnlyA
            ES <--> EP :: EqualAI
            E + P ⇌ EP :: EqualAI
        end
    end
    only_chem = @_testhelper_unguarded @enzyme_mechanism begin
        substrates: S ; products: P ; oligomeric_state: 2
        steps: begin
            E + S ⇌ ES :: EqualAI
            ES <--> EP :: OnlyA
            E + P ⇌ EP :: EqualAI
        end
    end
    s_chem = @_testhelper_unguarded @enzyme_mechanism begin
        substrates: S ; products: P ; oligomeric_state: 2
        steps: begin
            E + S ⇌ ES :: OnlyA
            ES <--> EP :: OnlyA
            E + P ⇌ EP :: EqualAI
        end
    end
    only_p = @_testhelper_unguarded @enzyme_mechanism begin
        substrates: S ; products: P ; oligomeric_state: 2
        steps: begin
            E + S ⇌ ES :: EqualAI
            ES <--> EP :: EqualAI
            E + P ⇌ EP :: OnlyA
        end
    end
    s_p = @_testhelper_unguarded @enzyme_mechanism begin
        substrates: S ; products: P ; oligomeric_state: 2
        steps: begin
            E + S ⇌ ES :: OnlyA
            ES <--> EP :: EqualAI
            E + P ⇌ EP :: OnlyA
        end
    end
    chem_p = @_testhelper_unguarded @enzyme_mechanism begin
        substrates: S ; products: P ; oligomeric_state: 2
        steps: begin
            E + S ⇌ ES :: EqualAI
            ES <--> EP :: OnlyA
            E + P ⇌ EP :: OnlyA
        end
    end
    every = @_testhelper_unguarded @enzyme_mechanism begin
        substrates: S ; products: P ; oligomeric_state: 2
        steps: begin
            E + S ⇌ ES :: OnlyA
            ES <--> EP :: OnlyA
            E + P ⇌ EP :: OnlyA
        end
    end
    chem_r_a = @_testhelper_unguarded @enzyme_mechanism begin
        substrates: S ; products: P ; oligomeric_state: 2
        steps: begin
            E + S ⇌ ES :: EqualAI
            ES <--> EP :: OnlyA
            E + P ⇌ EP :: EqualAI
        end
        regulatory_site(2): R::OnlyA
    end
    chem_r_i = @_testhelper_unguarded @enzyme_mechanism begin
        substrates: S ; products: P ; oligomeric_state: 2
        steps: begin
            E + S ⇌ ES :: EqualAI
            ES <--> EP :: OnlyA
            E + P ⇌ EP :: EqualAI
        end
        regulatory_site(2): R::OnlyI
    end
    kids = _testhelper_move(EnzymeRates.Enumeration._expand_to_allosteric, m,
                            uni_uni_allo_reg)
    @test length(kids) == 9
    @test Set(kids) == Set([only_s, only_chem, s_chem, only_p, s_p, chem_p, every,
                            chem_r_a, chem_r_i])
    @test count(k -> k.allosteric.allosteric_states ==
                     only_chem.allosteric.allosteric_states, kids) == 3
    # Whether `k` tags a binding group :OnlyA.
    _testhelper_binding_only_a(k) =
        any(g -> k.allosteric.allosteric_states[g] === :OnlyA &&
                 !_testhelper_is_catalytic_group(m, g),
            eachindex(k.allosteric.allosteric_states))
    @test all(k -> isempty(k.allosteric.regulatory_sites),
              filter(_testhelper_binding_only_a, kids))
    @test !any(k -> !isempty(k.allosteric.regulatory_sites) &&
                    :NonequalAI in k.allosteric.allosteric_states, kids)
end

@testset "to_allosteric emits a :NonequalAI set only where its :OnlyA set is refused" begin
    # The six-step ping-pong with both half-reactions at steady state: A and Q bind free
    # E, B and P bind the covalent form E_res_plusA_minusP. At one subunit with no
    # regulator the move emits the 63 :OnlyA sets. Tagged :NonequalAI, eight minimal sets
    # gain a constant beyond L: each half-reaction alone, and the pairs of bindings A with
    # Q, A with B, Q with B, A with P, Q with P and B with P; one binding alone gains
    # nothing. Tagged :OnlyA, three of the pairs are refused, A with B, Q with P and B with
    # P, and the move emits those three tagged :NonequalAI. The other five are returned
    # tagged :OnlyA, and `_expand_change_allo_state` reaches their :NonequalAI versions
    # from those children.
    m = @enzyme_mechanism begin
        substrates: A, B
        products: P, Q
        steps: begin
            E + A ⇌ EA
            E_res_plusA_minusP + B ⇌ EB_res_plusA_minusP
            E + Q ⇌ EQ
            E_res_plusA_minusP + P ⇌ EP_res_plusA_minusP
            EA <--> EP_res_plusA_minusP
            EB_res_plusA_minusP <--> EQ
        end
    end
    rxn = @enzyme_reaction begin
        substrates: A[CX], B[N]
        products: P[C], Q[NX]
    end
    a_b_nonequal = @enzyme_mechanism begin
        substrates: A, B ; products: P, Q ; oligomeric_state: 1
        steps: begin
            E + A ⇌ EA                                   :: NonequalAI
            E_res_plusA_minusP + B ⇌ EB_res_plusA_minusP :: NonequalAI
            E + Q ⇌ EQ                                   :: EqualAI
            E_res_plusA_minusP + P ⇌ EP_res_plusA_minusP :: EqualAI
            EA <--> EP_res_plusA_minusP                  :: EqualAI
            EB_res_plusA_minusP <--> EQ                  :: EqualAI
        end
    end
    q_p_nonequal = @enzyme_mechanism begin
        substrates: A, B ; products: P, Q ; oligomeric_state: 1
        steps: begin
            E + A ⇌ EA                                   :: EqualAI
            E_res_plusA_minusP + B ⇌ EB_res_plusA_minusP :: EqualAI
            E + Q ⇌ EQ                                   :: NonequalAI
            E_res_plusA_minusP + P ⇌ EP_res_plusA_minusP :: NonequalAI
            EA <--> EP_res_plusA_minusP                  :: EqualAI
            EB_res_plusA_minusP <--> EQ                  :: EqualAI
        end
    end
    b_p_nonequal = @enzyme_mechanism begin
        substrates: A, B ; products: P, Q ; oligomeric_state: 1
        steps: begin
            E + A ⇌ EA                                   :: EqualAI
            E_res_plusA_minusP + B ⇌ EB_res_plusA_minusP :: NonequalAI
            E + Q ⇌ EQ                                   :: EqualAI
            E_res_plusA_minusP + P ⇌ EP_res_plusA_minusP :: NonequalAI
            EA <--> EP_res_plusA_minusP                  :: EqualAI
            EB_res_plusA_minusP <--> EQ                  :: EqualAI
        end
    end
    a_b_only_a = @_testhelper_unguarded @enzyme_mechanism begin
        substrates: A, B ; products: P, Q ; oligomeric_state: 1
        steps: begin
            E + A ⇌ EA                                   :: OnlyA
            E_res_plusA_minusP + B ⇌ EB_res_plusA_minusP :: OnlyA
            E + Q ⇌ EQ                                   :: EqualAI
            E_res_plusA_minusP + P ⇌ EP_res_plusA_minusP :: EqualAI
            EA <--> EP_res_plusA_minusP                  :: EqualAI
            EB_res_plusA_minusP <--> EQ                  :: EqualAI
        end
    end
    q_p_only_a = @_testhelper_unguarded @enzyme_mechanism begin
        substrates: A, B ; products: P, Q ; oligomeric_state: 1
        steps: begin
            E + A ⇌ EA                                   :: EqualAI
            E_res_plusA_minusP + B ⇌ EB_res_plusA_minusP :: EqualAI
            E + Q ⇌ EQ                                   :: OnlyA
            E_res_plusA_minusP + P ⇌ EP_res_plusA_minusP :: OnlyA
            EA <--> EP_res_plusA_minusP                  :: EqualAI
            EB_res_plusA_minusP <--> EQ                  :: EqualAI
        end
    end
    b_p_only_a = @_testhelper_unguarded @enzyme_mechanism begin
        substrates: A, B ; products: P, Q ; oligomeric_state: 1
        steps: begin
            E + A ⇌ EA                                   :: EqualAI
            E_res_plusA_minusP + B ⇌ EB_res_plusA_minusP :: OnlyA
            E + Q ⇌ EQ                                   :: EqualAI
            E_res_plusA_minusP + P ⇌ EP_res_plusA_minusP :: OnlyA
            EA <--> EP_res_plusA_minusP                  :: EqualAI
            EB_res_plusA_minusP <--> EQ                  :: EqualAI
        end
    end
    a_q_only_a = @_testhelper_unguarded @enzyme_mechanism begin
        substrates: A, B ; products: P, Q ; oligomeric_state: 1
        steps: begin
            E + A ⇌ EA                                   :: OnlyA
            E_res_plusA_minusP + B ⇌ EB_res_plusA_minusP :: EqualAI
            E + Q ⇌ EQ                                   :: OnlyA
            E_res_plusA_minusP + P ⇌ EP_res_plusA_minusP :: EqualAI
            EA <--> EP_res_plusA_minusP                  :: EqualAI
            EB_res_plusA_minusP <--> EQ                  :: EqualAI
        end
    end
    a_p_only_a = @_testhelper_unguarded @enzyme_mechanism begin
        substrates: A, B ; products: P, Q ; oligomeric_state: 1
        steps: begin
            E + A ⇌ EA                                   :: OnlyA
            E_res_plusA_minusP + B ⇌ EB_res_plusA_minusP :: EqualAI
            E + Q ⇌ EQ                                   :: EqualAI
            E_res_plusA_minusP + P ⇌ EP_res_plusA_minusP :: OnlyA
            EA <--> EP_res_plusA_minusP                  :: EqualAI
            EB_res_plusA_minusP <--> EQ                  :: EqualAI
        end
    end
    q_b_only_a = @_testhelper_unguarded @enzyme_mechanism begin
        substrates: A, B ; products: P, Q ; oligomeric_state: 1
        steps: begin
            E + A ⇌ EA                                   :: EqualAI
            E_res_plusA_minusP + B ⇌ EB_res_plusA_minusP :: OnlyA
            E + Q ⇌ EQ                                   :: OnlyA
            E_res_plusA_minusP + P ⇌ EP_res_plusA_minusP :: EqualAI
            EA <--> EP_res_plusA_minusP                  :: EqualAI
            EB_res_plusA_minusP <--> EQ                  :: EqualAI
        end
    end
    half1_only_a = @_testhelper_unguarded @enzyme_mechanism begin
        substrates: A, B ; products: P, Q ; oligomeric_state: 1
        steps: begin
            E + A ⇌ EA                                   :: EqualAI
            E_res_plusA_minusP + B ⇌ EB_res_plusA_minusP :: EqualAI
            E + Q ⇌ EQ                                   :: EqualAI
            E_res_plusA_minusP + P ⇌ EP_res_plusA_minusP :: EqualAI
            EA <--> EP_res_plusA_minusP                  :: OnlyA
            EB_res_plusA_minusP <--> EQ                  :: EqualAI
        end
    end
    half2_only_a = @_testhelper_unguarded @enzyme_mechanism begin
        substrates: A, B ; products: P, Q ; oligomeric_state: 1
        steps: begin
            E + A ⇌ EA                                   :: EqualAI
            E_res_plusA_minusP + B ⇌ EB_res_plusA_minusP :: EqualAI
            E + Q ⇌ EQ                                   :: EqualAI
            E_res_plusA_minusP + P ⇌ EP_res_plusA_minusP :: EqualAI
            EA <--> EP_res_plusA_minusP                  :: EqualAI
            EB_res_plusA_minusP <--> EQ                  :: OnlyA
        end
    end
    kids = _testhelper_move(EnzymeRates.Enumeration._expand_to_allosteric, m, rxn)
    @test length(kids) == 63 + 3
    @test allunique(kids)
    nonequal = filter(k -> :NonequalAI in k.allosteric.allosteric_states, kids)
    @test length(nonequal) == 3
    @test Set(nonequal) == Set([a_b_nonequal, q_p_nonequal, b_p_nonequal])
    # The other 63, distinct and each with only :OnlyA and :EqualAI groups, are the 63
    # nonempty sets tagged :OnlyA.
    @test all(setdiff(kids, nonequal)) do k
        k.catalytic == m.catalytic && k.allosteric.oligomeric_state == 1 &&
            isempty(k.allosteric.regulatory_sites) &&
            :OnlyA in k.allosteric.allosteric_states &&
            issubset(k.allosteric.allosteric_states, (:OnlyA, :EqualAI))
    end
    returns(x) = EnzymeRates.Enumeration._enumerator_returns(_testhelper_candidate(x))
    @test !any(returns, [a_b_only_a, q_p_only_a, b_p_only_a])
    @test all(returns, [a_q_only_a, a_p_only_a, q_b_only_a, half1_only_a, half2_only_a])
end

@testset "to_allosteric on a random-order ping-pong: the oracles on two children" begin
    # A ping-pong of census pingpong: A and Q bind free E at random, P and B bind the
    # covalent form E_res_plusA_minusP, and two chemistry groups, EA to EP_res_plusA_minusP
    # and EB_res_plusA_minusP to EQ, make the two half-reactions. At one subunit the move
    # emits the 63 sets of the six groups tagged :OnlyA. Tagged :NonequalAI, six minimal
    # sets gain a constant beyond L, each a pair of groups: A's with Q's, A's with the
    # second half-reaction, A's with B's, Q's with B's, Q's with the second half-reaction,
    # and B's with the second half-reaction. Tagged :OnlyA the last four are refused, and
    # the move emits them tagged :NonequalAI. The testset writes those four out, and two of
    # the 63 :OnlyA children. With both half-reactions :OnlyA (`both_halves`) the inactive
    # conformation never reaches the covalent forms: both oracles find the law to be the
    # limit on every face, and the enumerator returns it. With the first half-reaction
    # alone :OnlyA (`half1`) admission's interior accepts it, but a face holds a
    # rapid-equilibrium segment whose lowest forms split (`_SPLIT`), so the enumerator does
    # not return it.
    m = @enzyme_mechanism begin
        substrates: A, B
        products: P, Q
        steps: begin
            (E + A ⇌ EA, EQ + A ⇌ EAQ)
            (E + Q ⇌ EQ, EA + Q ⇌ EAQ)
            E_res_plusA_minusP + B ⇌ EB_res_plusA_minusP
            E_res_plusA_minusP + P ⇌ EP_res_plusA_minusP
            EA <--> EP_res_plusA_minusP
            EB_res_plusA_minusP ⇌ EQ
        end
    end
    rxn = @enzyme_reaction begin
        substrates: A[CX], B[N]
        products: P[C], Q[NX]
    end
    half1 = @_testhelper_unguarded @enzyme_mechanism begin
        substrates: A, B ; products: P, Q ; oligomeric_state: 1
        steps: begin
            (E + A ⇌ EA, EQ + A ⇌ EAQ)                 :: EqualAI
            (E + Q ⇌ EQ, EA + Q ⇌ EAQ)                 :: EqualAI
            E_res_plusA_minusP + B ⇌ EB_res_plusA_minusP :: EqualAI
            E_res_plusA_minusP + P ⇌ EP_res_plusA_minusP :: EqualAI
            EA <--> EP_res_plusA_minusP                   :: OnlyA
            EB_res_plusA_minusP ⇌ EQ                      :: EqualAI
        end
    end
    both_halves = @_testhelper_unguarded @enzyme_mechanism begin
        substrates: A, B ; products: P, Q ; oligomeric_state: 1
        steps: begin
            (E + A ⇌ EA, EQ + A ⇌ EAQ)                 :: EqualAI
            (E + Q ⇌ EQ, EA + Q ⇌ EAQ)                 :: EqualAI
            E_res_plusA_minusP + B ⇌ EB_res_plusA_minusP :: EqualAI
            E_res_plusA_minusP + P ⇌ EP_res_plusA_minusP :: EqualAI
            EA <--> EP_res_plusA_minusP                   :: OnlyA
            EB_res_plusA_minusP ⇌ EQ                      :: OnlyA
        end
    end
    a_b_nonequal = @enzyme_mechanism begin
        substrates: A, B ; products: P, Q ; oligomeric_state: 1
        steps: begin
            (E + A ⇌ EA, EQ + A ⇌ EAQ)                 :: NonequalAI
            (E + Q ⇌ EQ, EA + Q ⇌ EAQ)                 :: EqualAI
            E_res_plusA_minusP + B ⇌ EB_res_plusA_minusP :: NonequalAI
            E_res_plusA_minusP + P ⇌ EP_res_plusA_minusP :: EqualAI
            EA <--> EP_res_plusA_minusP                   :: EqualAI
            EB_res_plusA_minusP ⇌ EQ                      :: EqualAI
        end
    end
    q_b_nonequal = @enzyme_mechanism begin
        substrates: A, B ; products: P, Q ; oligomeric_state: 1
        steps: begin
            (E + A ⇌ EA, EQ + A ⇌ EAQ)                 :: EqualAI
            (E + Q ⇌ EQ, EA + Q ⇌ EAQ)                 :: NonequalAI
            E_res_plusA_minusP + B ⇌ EB_res_plusA_minusP :: NonequalAI
            E_res_plusA_minusP + P ⇌ EP_res_plusA_minusP :: EqualAI
            EA <--> EP_res_plusA_minusP                   :: EqualAI
            EB_res_plusA_minusP ⇌ EQ                      :: EqualAI
        end
    end
    q_half2_nonequal = @enzyme_mechanism begin
        substrates: A, B ; products: P, Q ; oligomeric_state: 1
        steps: begin
            (E + A ⇌ EA, EQ + A ⇌ EAQ)                 :: EqualAI
            (E + Q ⇌ EQ, EA + Q ⇌ EAQ)                 :: NonequalAI
            E_res_plusA_minusP + B ⇌ EB_res_plusA_minusP :: EqualAI
            E_res_plusA_minusP + P ⇌ EP_res_plusA_minusP :: EqualAI
            EA <--> EP_res_plusA_minusP                   :: EqualAI
            EB_res_plusA_minusP ⇌ EQ                      :: NonequalAI
        end
    end
    b_half2_nonequal = @enzyme_mechanism begin
        substrates: A, B ; products: P, Q ; oligomeric_state: 1
        steps: begin
            (E + A ⇌ EA, EQ + A ⇌ EAQ)                 :: EqualAI
            (E + Q ⇌ EQ, EA + Q ⇌ EAQ)                 :: EqualAI
            E_res_plusA_minusP + B ⇌ EB_res_plusA_minusP :: NonequalAI
            E_res_plusA_minusP + P ⇌ EP_res_plusA_minusP :: EqualAI
            EA <--> EP_res_plusA_minusP                   :: EqualAI
            EB_res_plusA_minusP ⇌ EQ                      :: NonequalAI
        end
    end
    kids = _testhelper_move(EnzymeRates.Enumeration._expand_to_allosteric, m, rxn)
    @test length(kids) == 63 + 4
    @test allunique(kids)
    nonequal = filter(k -> :NonequalAI in k.allosteric.allosteric_states, kids)
    @test Set(nonequal) ==
          Set([a_b_nonequal, q_b_nonequal, q_half2_nonequal, b_half2_nonequal])
    # The other 63, distinct and each with only :OnlyA and :EqualAI groups, are the 63
    # nonempty sets tagged :OnlyA.
    @test all(setdiff(kids, nonequal)) do k
        k.catalytic == m.catalytic && k.allosteric.oligomeric_state == 1 &&
            isempty(k.allosteric.regulatory_sites) &&
            :OnlyA in k.allosteric.allosteric_states &&
            issubset(k.allosteric.allosteric_states, (:OnlyA, :EqualAI))
    end
    @test half1 in kids && both_halves in kids
    D = EnzymeRates.Derivation
    returns(x) = EnzymeRates.Enumeration._enumerator_returns(_testhelper_candidate(x))
    @test D._admission_interior(half1, D._built_reading(half1)) === D._WELL_POSED
    @test D._guard(half1).fault === D._SPLIT
    @test !returns(half1)
    @test all(_testhelper_admission_agreement(both_halves;
                                              rng = Random.MersenneTwister(15327),
                                              corners = true))
    _testhelper_check_law_against_limit(both_halves)
    @test returns(both_halves)
end

@testset "to_allosteric on a ter-bi ping-pong: the oracles on every group :OnlyA" begin
    # Census ter-bi parent 0dc399e7dab6e4f5 (census/base6): its covalent form E(; residual
    # = B + D - Q) takes up A and gives off P. Four groups hold chemistry, those of B, D, A
    # and Q, each mixing a plain binding with a fused one; P's binding at E is the only
    # binding group. At one subunit the move emits the 31 sets of the five groups tagged
    # :OnlyA. Tagged :NonequalAI, four minimal sets gain a constant beyond L: the groups of
    # B, of P and of Q alone, and those of D and A together. Tagged :OnlyA, B's group and
    # the pair of D's and A's are refused, and the move emits those two tagged
    # :NonequalAI. The testset writes those two out, and three of the 31 :OnlyA children.
    # The four chemistry groups :OnlyA (`chemistry`) fail the guard: at zero A, B and D the
    # fused A binding, vanishing one way into the covalent form, fills it faster than it
    # empties (`_FILLED`). B's group alone (`b_alone`) passes the guard but keeps EB,
    # which only the binding it excludes enters (`_kept_forms_keep_a_binding`). With every
    # group :OnlyA (`every`) both oracles find the law to be the limit on every face, and
    # the enumerator returns it, but neither of the other two.
    m = @enzyme_mechanism begin
        substrates: A, B, D
        products: P, Q
        steps: begin
            (E + B <--> EB, ED + B <--> EQ_res_plusB_plusD_minusQ)
            (E + D ⇌ ED, EB + D ⇌ EQ_res_plusB_plusD_minusQ)
            E + P <--> EP
            (E_res_plusB_plusD_minusQ + A ⇌ EP,
             EQ_res_plusB_plusD_minusQ + A ⇌ EAQ_res_plusB_plusD_minusQ)
            (E_res_plusB_plusD_minusQ + Q <--> EQ_res_plusB_plusD_minusQ,
             EP + Q <--> EAQ_res_plusB_plusD_minusQ)
        end
    end
    rxn = @enzyme_reaction begin
        substrates: A[C], B[N], D[X]
        products: P[CN], Q[X]
    end
    chemistry = @_testhelper_unguarded @enzyme_mechanism begin
        substrates: A, B, D ; products: P, Q ; oligomeric_state: 1
        steps: begin
            (E + B <--> EB, ED + B <--> EQ_res_plusB_plusD_minusQ)      :: OnlyA
            (E + D ⇌ ED, EB + D ⇌ EQ_res_plusB_plusD_minusQ)            :: OnlyA
            E + P <--> EP                                                  :: EqualAI
            (E_res_plusB_plusD_minusQ + A ⇌ EP,
             EQ_res_plusB_plusD_minusQ + A ⇌ EAQ_res_plusB_plusD_minusQ) :: OnlyA
            (E_res_plusB_plusD_minusQ + Q <--> EQ_res_plusB_plusD_minusQ,
             EP + Q <--> EAQ_res_plusB_plusD_minusQ)                    :: OnlyA
        end
    end
    b_alone = @enzyme_mechanism begin
        substrates: A, B, D ; products: P, Q ; oligomeric_state: 1
        steps: begin
            (E + B <--> EB, ED + B <--> EQ_res_plusB_plusD_minusQ)      :: OnlyA
            (E + D ⇌ ED, EB + D ⇌ EQ_res_plusB_plusD_minusQ)            :: EqualAI
            E + P <--> EP                                                  :: EqualAI
            (E_res_plusB_plusD_minusQ + A ⇌ EP,
             EQ_res_plusB_plusD_minusQ + A ⇌ EAQ_res_plusB_plusD_minusQ) :: EqualAI
            (E_res_plusB_plusD_minusQ + Q <--> EQ_res_plusB_plusD_minusQ,
             EP + Q <--> EAQ_res_plusB_plusD_minusQ)                    :: EqualAI
        end
    end
    every = @enzyme_mechanism begin
        substrates: A, B, D ; products: P, Q ; oligomeric_state: 1
        steps: begin
            (E + B <--> EB, ED + B <--> EQ_res_plusB_plusD_minusQ)      :: OnlyA
            (E + D ⇌ ED, EB + D ⇌ EQ_res_plusB_plusD_minusQ)            :: OnlyA
            E + P <--> EP                                                  :: OnlyA
            (E_res_plusB_plusD_minusQ + A ⇌ EP,
             EQ_res_plusB_plusD_minusQ + A ⇌ EAQ_res_plusB_plusD_minusQ) :: OnlyA
            (E_res_plusB_plusD_minusQ + Q <--> EQ_res_plusB_plusD_minusQ,
             EP + Q <--> EAQ_res_plusB_plusD_minusQ)                    :: OnlyA
        end
    end
    b_nonequal = @enzyme_mechanism begin
        substrates: A, B, D ; products: P, Q ; oligomeric_state: 1
        steps: begin
            (E + B <--> EB, ED + B <--> EQ_res_plusB_plusD_minusQ)      :: NonequalAI
            (E + D ⇌ ED, EB + D ⇌ EQ_res_plusB_plusD_minusQ)            :: EqualAI
            E + P <--> EP                                                  :: EqualAI
            (E_res_plusB_plusD_minusQ + A ⇌ EP,
             EQ_res_plusB_plusD_minusQ + A ⇌ EAQ_res_plusB_plusD_minusQ) :: EqualAI
            (E_res_plusB_plusD_minusQ + Q <--> EQ_res_plusB_plusD_minusQ,
             EP + Q <--> EAQ_res_plusB_plusD_minusQ)                    :: EqualAI
        end
    end
    d_a_nonequal = @enzyme_mechanism begin
        substrates: A, B, D ; products: P, Q ; oligomeric_state: 1
        steps: begin
            (E + B <--> EB, ED + B <--> EQ_res_plusB_plusD_minusQ)      :: EqualAI
            (E + D ⇌ ED, EB + D ⇌ EQ_res_plusB_plusD_minusQ)            :: NonequalAI
            E + P <--> EP                                                  :: EqualAI
            (E_res_plusB_plusD_minusQ + A ⇌ EP,
             EQ_res_plusB_plusD_minusQ + A ⇌ EAQ_res_plusB_plusD_minusQ) :: NonequalAI
            (E_res_plusB_plusD_minusQ + Q <--> EQ_res_plusB_plusD_minusQ,
             EP + Q <--> EAQ_res_plusB_plusD_minusQ)                    :: EqualAI
        end
    end
    kids = _testhelper_move(EnzymeRates.Enumeration._expand_to_allosteric, m, rxn)
    @test length(kids) == 31 + 2
    @test allunique(kids)
    nonequal = filter(k -> :NonequalAI in k.allosteric.allosteric_states, kids)
    @test Set(nonequal) == Set([b_nonequal, d_a_nonequal])
    # The other 31, distinct and each with only :OnlyA and :EqualAI groups, are the 31
    # nonempty sets tagged :OnlyA.
    @test all(setdiff(kids, nonequal)) do k
        k.catalytic == m.catalytic && k.allosteric.oligomeric_state == 1 &&
            isempty(k.allosteric.regulatory_sites) &&
            :OnlyA in k.allosteric.allosteric_states &&
            issubset(k.allosteric.allosteric_states, (:OnlyA, :EqualAI))
    end
    @test chemistry in kids && b_alone in kids && every in kids
    D = EnzymeRates.Derivation
    v = D._guard(chemistry)
    @test (v.fault, EnzymeRates.Mechanisms.step_text(chemistry.catalytic, v.edge)) ==
          (D._FILLED, "E_res_plusB_plusD_minusQ + A ⇌ EP")
    @test _testhelper_is_valid(b_alone) &&
          !EnzymeRates.Enumeration._kept_forms_keep_a_binding(b_alone)
    @test all(_testhelper_admission_agreement(every; rng = Random.MersenneTwister(13964)))
    _testhelper_check_law_against_limit(every)
    returns(x) = EnzymeRates.Enumeration._enumerator_returns(_testhelper_candidate(x))
    @test returns(every) && !returns(chemistry) && !returns(b_alone)
end

@testset "a ter-bi ping-pong's chemistry groups, each :OnlyA alone, are refused" begin
    # A ter-bi ping-pong of the census ter-bi reaction: B and D bind free E in either
    # order, the second fused into the product-bound covalent form E(Q; residual = B + D -
    # Q) (B's and D's groups, each chemistry); Q's release from it leaves E(; residual =
    # B + D - Q), which takes up A into EP, fused (A's group), and Q's group holds Q's
    # abortive binding at EP beside its release. P's binding at E is the one binding
    # group. At one subunit the move emits the 31 sets of the five groups tagged :OnlyA.
    # Tagged :NonequalAI, four minimal sets gain a constant beyond L: the groups of B, of D
    # and of P alone, and those of Q and A together. Tagged :OnlyA, all but P's group are
    # refused, and the move emits those three tagged :NonequalAI. The testset writes those
    # three out, and the four :OnlyA children that hold one chemistry group each, none of
    # which the enumerator returns. Q's group alone keeps the two sides of its
    # rapid-equilibrium release joined in the inactive conformation (`_FAST_STEP_PINNED`),
    # and A's alone those of its fused binding; the guard names the group's first step in
    # step order. B's group alone passes the guard but keeps EB, which only the binding it
    # excludes enters (`_kept_forms_keep_a_binding`), and D's alone keeps ED so.
    m = @enzyme_mechanism begin
        substrates: A, B, D
        products: P, Q
        steps: begin
            (E + B <--> EB, ED + B <--> EQ_res_plusB_plusD_minusQ)
            (E + D <--> ED, EB + D <--> EQ_res_plusB_plusD_minusQ)
            E + P <--> EP
            (E_res_plusB_plusD_minusQ + A ⇌ EP,
             EQ_res_plusB_plusD_minusQ + A ⇌ EAQ_res_plusB_plusD_minusQ)
            (E_res_plusB_plusD_minusQ + Q ⇌ EQ_res_plusB_plusD_minusQ,
             EP + Q ⇌ EAQ_res_plusB_plusD_minusQ)
        end
    end
    rxn = @enzyme_reaction begin
        substrates: A[C], B[N], D[X]
        products: P[CN], Q[X]
    end
    q_alone = @_testhelper_unguarded @enzyme_mechanism begin
        substrates: A, B, D ; products: P, Q ; oligomeric_state: 1
        steps: begin
            (E + B <--> EB, ED + B <--> EQ_res_plusB_plusD_minusQ)      :: EqualAI
            (E + D <--> ED, EB + D <--> EQ_res_plusB_plusD_minusQ)      :: EqualAI
            E + P <--> EP                                                  :: EqualAI
            (E_res_plusB_plusD_minusQ + A ⇌ EP,
             EQ_res_plusB_plusD_minusQ + A ⇌ EAQ_res_plusB_plusD_minusQ) :: EqualAI
            (E_res_plusB_plusD_minusQ + Q ⇌ EQ_res_plusB_plusD_minusQ,
             EP + Q ⇌ EAQ_res_plusB_plusD_minusQ)                       :: OnlyA
        end
    end
    a_alone = @_testhelper_unguarded @enzyme_mechanism begin
        substrates: A, B, D ; products: P, Q ; oligomeric_state: 1
        steps: begin
            (E + B <--> EB, ED + B <--> EQ_res_plusB_plusD_minusQ)      :: EqualAI
            (E + D <--> ED, EB + D <--> EQ_res_plusB_plusD_minusQ)      :: EqualAI
            E + P <--> EP                                                  :: EqualAI
            (E_res_plusB_plusD_minusQ + A ⇌ EP,
             EQ_res_plusB_plusD_minusQ + A ⇌ EAQ_res_plusB_plusD_minusQ) :: OnlyA
            (E_res_plusB_plusD_minusQ + Q ⇌ EQ_res_plusB_plusD_minusQ,
             EP + Q ⇌ EAQ_res_plusB_plusD_minusQ)                       :: EqualAI
        end
    end
    b_alone = @enzyme_mechanism begin
        substrates: A, B, D ; products: P, Q ; oligomeric_state: 1
        steps: begin
            (E + B <--> EB, ED + B <--> EQ_res_plusB_plusD_minusQ)      :: OnlyA
            (E + D <--> ED, EB + D <--> EQ_res_plusB_plusD_minusQ)      :: EqualAI
            E + P <--> EP                                                  :: EqualAI
            (E_res_plusB_plusD_minusQ + A ⇌ EP,
             EQ_res_plusB_plusD_minusQ + A ⇌ EAQ_res_plusB_plusD_minusQ) :: EqualAI
            (E_res_plusB_plusD_minusQ + Q ⇌ EQ_res_plusB_plusD_minusQ,
             EP + Q ⇌ EAQ_res_plusB_plusD_minusQ)                       :: EqualAI
        end
    end
    d_alone = @enzyme_mechanism begin
        substrates: A, B, D ; products: P, Q ; oligomeric_state: 1
        steps: begin
            (E + B <--> EB, ED + B <--> EQ_res_plusB_plusD_minusQ)      :: EqualAI
            (E + D <--> ED, EB + D <--> EQ_res_plusB_plusD_minusQ)      :: OnlyA
            E + P <--> EP                                                  :: EqualAI
            (E_res_plusB_plusD_minusQ + A ⇌ EP,
             EQ_res_plusB_plusD_minusQ + A ⇌ EAQ_res_plusB_plusD_minusQ) :: EqualAI
            (E_res_plusB_plusD_minusQ + Q ⇌ EQ_res_plusB_plusD_minusQ,
             EP + Q ⇌ EAQ_res_plusB_plusD_minusQ)                       :: EqualAI
        end
    end
    b_nonequal = @enzyme_mechanism begin
        substrates: A, B, D ; products: P, Q ; oligomeric_state: 1
        steps: begin
            (E + B <--> EB, ED + B <--> EQ_res_plusB_plusD_minusQ)      :: NonequalAI
            (E + D <--> ED, EB + D <--> EQ_res_plusB_plusD_minusQ)      :: EqualAI
            E + P <--> EP                                                  :: EqualAI
            (E_res_plusB_plusD_minusQ + A ⇌ EP,
             EQ_res_plusB_plusD_minusQ + A ⇌ EAQ_res_plusB_plusD_minusQ) :: EqualAI
            (E_res_plusB_plusD_minusQ + Q ⇌ EQ_res_plusB_plusD_minusQ,
             EP + Q ⇌ EAQ_res_plusB_plusD_minusQ)                       :: EqualAI
        end
    end
    d_nonequal = @enzyme_mechanism begin
        substrates: A, B, D ; products: P, Q ; oligomeric_state: 1
        steps: begin
            (E + B <--> EB, ED + B <--> EQ_res_plusB_plusD_minusQ)      :: EqualAI
            (E + D <--> ED, EB + D <--> EQ_res_plusB_plusD_minusQ)      :: NonequalAI
            E + P <--> EP                                                  :: EqualAI
            (E_res_plusB_plusD_minusQ + A ⇌ EP,
             EQ_res_plusB_plusD_minusQ + A ⇌ EAQ_res_plusB_plusD_minusQ) :: EqualAI
            (E_res_plusB_plusD_minusQ + Q ⇌ EQ_res_plusB_plusD_minusQ,
             EP + Q ⇌ EAQ_res_plusB_plusD_minusQ)                       :: EqualAI
        end
    end
    q_a_nonequal = @enzyme_mechanism begin
        substrates: A, B, D ; products: P, Q ; oligomeric_state: 1
        steps: begin
            (E + B <--> EB, ED + B <--> EQ_res_plusB_plusD_minusQ)      :: EqualAI
            (E + D <--> ED, EB + D <--> EQ_res_plusB_plusD_minusQ)      :: EqualAI
            E + P <--> EP                                                  :: EqualAI
            (E_res_plusB_plusD_minusQ + A ⇌ EP,
             EQ_res_plusB_plusD_minusQ + A ⇌ EAQ_res_plusB_plusD_minusQ) :: NonequalAI
            (E_res_plusB_plusD_minusQ + Q ⇌ EQ_res_plusB_plusD_minusQ,
             EP + Q ⇌ EAQ_res_plusB_plusD_minusQ)                       :: NonequalAI
        end
    end
    kids = _testhelper_move(EnzymeRates.Enumeration._expand_to_allosteric, m, rxn)
    @test length(kids) == 31 + 3
    @test allunique(kids)
    nonequal = filter(k -> :NonequalAI in k.allosteric.allosteric_states, kids)
    @test Set(nonequal) == Set([b_nonequal, d_nonequal, q_a_nonequal])
    # The other 31, distinct and each with only :OnlyA and :EqualAI groups, are the 31
    # nonempty sets tagged :OnlyA.
    @test all(setdiff(kids, nonequal)) do k
        k.catalytic == m.catalytic && k.allosteric.oligomeric_state == 1 &&
            isempty(k.allosteric.regulatory_sites) &&
            :OnlyA in k.allosteric.allosteric_states &&
            issubset(k.allosteric.allosteric_states, (:OnlyA, :EqualAI))
    end
    @test all(in(kids), [q_alone, a_alone, b_alone, d_alone])
    D = EnzymeRates.Derivation
    for (x, step) in (
            (q_alone, "EP + Q ⇌ EAQ_res_plusB_plusD_minusQ"),
            (a_alone, "EQ_res_plusB_plusD_minusQ + A ⇌ EAQ_res_plusB_plusD_minusQ"))
        v = D._guard(x)
        @test (v.fault, EnzymeRates.Mechanisms.step_text(x.catalytic, v.edge)) ==
              (D._FAST_STEP_PINNED, step)
    end
    for x in (b_alone, d_alone)
        @test _testhelper_is_valid(x)
        @test !EnzymeRates.Enumeration._kept_forms_keep_a_binding(x)
    end
    @test !any(x -> EnzymeRates.Enumeration._enumerator_returns(_testhelper_candidate(x)),
               [q_alone, a_alone, b_alone, d_alone])
end

@testset "Mechanism — enumerates all allowed multiplicities" begin
    # Multi-valued allowed_oligomeric_states → the move's children are emitted once per
    # multiplicity: each nonempty set of the three groups :OnlyA. No set tagged
    # :NonequalAI is refused as :OnlyA, so there is no :NonequalAI child.
    rxn = @enzyme_reaction begin
        substrates: S[C]
        products: P[C]
        oligomeric_state: (2, 4)
    end
    m = @enzyme_mechanism begin
        substrates: S
        products: P
        steps: begin
            E + S ⇌ ES
            ES <--> EP
            EP ⇌ E + P
        end
    end
    allo = _testhelper_move(EnzymeRates.Enumeration._expand_to_allosteric, m, rxn)
    # Each child's comment names its :OnlyA groups; the first seven are at multiplicity
    # 2, the last seven at 4.
    # Untyped: inside `T[...]` the parser reads `begin` as an index keyword.
    expected = [
        # S
        (@_testhelper_unguarded @enzyme_mechanism begin
             substrates: S ; products: P ; oligomeric_state: 2
             steps: begin
                 E + S ⇌ ES :: OnlyA
                 ES <--> EP :: EqualAI
                 EP ⇌ E + P :: EqualAI
             end
         end),
        # the chemistry
        (@_testhelper_unguarded @enzyme_mechanism begin
             substrates: S ; products: P ; oligomeric_state: 2
             steps: begin
                 E + S ⇌ ES :: EqualAI
                 ES <--> EP :: OnlyA
                 EP ⇌ E + P :: EqualAI
             end
         end),
        # S, the chemistry
        (@_testhelper_unguarded @enzyme_mechanism begin
             substrates: S ; products: P ; oligomeric_state: 2
             steps: begin
                 E + S ⇌ ES :: OnlyA
                 ES <--> EP :: OnlyA
                 EP ⇌ E + P :: EqualAI
             end
         end),
        # P
        (@_testhelper_unguarded @enzyme_mechanism begin
             substrates: S ; products: P ; oligomeric_state: 2
             steps: begin
                 E + S ⇌ ES :: EqualAI
                 ES <--> EP :: EqualAI
                 EP ⇌ E + P :: OnlyA
             end
         end),
        # S, P
        (@_testhelper_unguarded @enzyme_mechanism begin
             substrates: S ; products: P ; oligomeric_state: 2
             steps: begin
                 E + S ⇌ ES :: OnlyA
                 ES <--> EP :: EqualAI
                 EP ⇌ E + P :: OnlyA
             end
         end),
        # the chemistry, P
        (@_testhelper_unguarded @enzyme_mechanism begin
             substrates: S ; products: P ; oligomeric_state: 2
             steps: begin
                 E + S ⇌ ES :: EqualAI
                 ES <--> EP :: OnlyA
                 EP ⇌ E + P :: OnlyA
             end
         end),
        # S, the chemistry, P
        (@_testhelper_unguarded @enzyme_mechanism begin
             substrates: S ; products: P ; oligomeric_state: 2
             steps: begin
                 E + S ⇌ ES :: OnlyA
                 ES <--> EP :: OnlyA
                 EP ⇌ E + P :: OnlyA
             end
         end),
        # S
        (@_testhelper_unguarded @enzyme_mechanism begin
             substrates: S ; products: P ; oligomeric_state: 4
             steps: begin
                 E + S ⇌ ES :: OnlyA
                 ES <--> EP :: EqualAI
                 EP ⇌ E + P :: EqualAI
             end
         end),
        # the chemistry
        (@_testhelper_unguarded @enzyme_mechanism begin
             substrates: S ; products: P ; oligomeric_state: 4
             steps: begin
                 E + S ⇌ ES :: EqualAI
                 ES <--> EP :: OnlyA
                 EP ⇌ E + P :: EqualAI
             end
         end),
        # S, the chemistry
        (@_testhelper_unguarded @enzyme_mechanism begin
             substrates: S ; products: P ; oligomeric_state: 4
             steps: begin
                 E + S ⇌ ES :: OnlyA
                 ES <--> EP :: OnlyA
                 EP ⇌ E + P :: EqualAI
             end
         end),
        # P
        (@_testhelper_unguarded @enzyme_mechanism begin
             substrates: S ; products: P ; oligomeric_state: 4
             steps: begin
                 E + S ⇌ ES :: EqualAI
                 ES <--> EP :: EqualAI
                 EP ⇌ E + P :: OnlyA
             end
         end),
        # S, P
        (@_testhelper_unguarded @enzyme_mechanism begin
             substrates: S ; products: P ; oligomeric_state: 4
             steps: begin
                 E + S ⇌ ES :: OnlyA
                 ES <--> EP :: EqualAI
                 EP ⇌ E + P :: OnlyA
             end
         end),
        # the chemistry, P
        (@_testhelper_unguarded @enzyme_mechanism begin
             substrates: S ; products: P ; oligomeric_state: 4
             steps: begin
                 E + S ⇌ ES :: EqualAI
                 ES <--> EP :: OnlyA
                 EP ⇌ E + P :: OnlyA
             end
         end),
        # S, the chemistry, P
        (@_testhelper_unguarded @enzyme_mechanism begin
             substrates: S ; products: P ; oligomeric_state: 4
             steps: begin
                 E + S ⇌ ES :: OnlyA
                 ES <--> EP :: OnlyA
                 EP ⇌ E + P :: OnlyA
             end
         end)]
    @test length(allo) == 14
    @test Set(allo) == Set(expected)
    mults = Set(am.allosteric.oligomeric_state for am in allo)
    @test mults == Set([2, 4])
    # Each multiplicity carries the same seven state sets: no all-:EqualAI one, and
    # three, three and one with one, two and three groups :OnlyA.
    for cn in (2, 4)
        cn_variants = filter(am -> am.allosteric.oligomeric_state == cn, allo)
        @test count(am -> all(==(:EqualAI), am.allosteric.allosteric_states),
                    cn_variants) == 0
        @test [count(am -> count(==(:OnlyA), am.allosteric.allosteric_states) == k,
                     cn_variants) for k in 1:3] == [3, 3, 1]
    end
end

@testset "Regulator children: parameter counts, no competitive inhibitor at a site" begin
    # uni_uni_allo_2reg declares two allosteric regulators (R1, R2). From a uni-uni seed
    # (RE bindings, SS chemistry), _expand_to_allosteric emits:
    #   * the 7 nonempty sets of the three groups tagged :OnlyA, no regulatory site, at
    #     +1 param each, the conformational constant L;
    #   * the set that holds the chemistry alone again with one regulator at a site, at
    #     +2 params each, L plus the regulator's K: one per (regulator, tag), regulator ∈
    #     {R1, R2}, tag ∈ {:OnlyA, :OnlyI}.
    # No set tagged :NonequalAI is refused as :OnlyA, so no child is :NonequalAI. Each
    # child's comment names its :OnlyA groups.
    seed = @enzyme_mechanism begin
        substrates: S
        products: P
        steps: begin
            E + S ⇌ ES
            ES <--> EP
            E + P ⇌ EP
        end
    end
    fitted = EnzymeRates.Derivation.fitted_parameter_count
    Δ(am) = fitted(am) - fitted(seed)
    outs = _testhelper_move(EnzymeRates.Enumeration._expand_to_allosteric, seed,
                            uni_uni_allo_2reg)
    # Untyped: inside `T[...]` the parser reads `begin` as an index keyword.
    only_a_expected = [
        # P
        (@_testhelper_unguarded @enzyme_mechanism begin
             substrates: S ; products: P ; oligomeric_state: 2
             steps: begin
                 E + P ⇌ EP :: OnlyA
                 E + S ⇌ ES :: EqualAI
                 ES <--> EP :: EqualAI
             end
         end),
        # S
        (@_testhelper_unguarded @enzyme_mechanism begin
             substrates: S ; products: P ; oligomeric_state: 2
             steps: begin
                 E + P ⇌ EP :: EqualAI
                 E + S ⇌ ES :: OnlyA
                 ES <--> EP :: EqualAI
             end
         end),
        # P, S
        (@_testhelper_unguarded @enzyme_mechanism begin
             substrates: S ; products: P ; oligomeric_state: 2
             steps: begin
                 E + P ⇌ EP :: OnlyA
                 E + S ⇌ ES :: OnlyA
                 ES <--> EP :: EqualAI
             end
         end),
        # the chemistry
        (@_testhelper_unguarded @enzyme_mechanism begin
             substrates: S ; products: P ; oligomeric_state: 2
             steps: begin
                 E + P ⇌ EP :: EqualAI
                 E + S ⇌ ES :: EqualAI
                 ES <--> EP :: OnlyA
             end
         end),
        # P, the chemistry
        (@_testhelper_unguarded @enzyme_mechanism begin
             substrates: S ; products: P ; oligomeric_state: 2
             steps: begin
                 E + P ⇌ EP :: OnlyA
                 E + S ⇌ ES :: EqualAI
                 ES <--> EP :: OnlyA
             end
         end),
        # S, the chemistry
        (@_testhelper_unguarded @enzyme_mechanism begin
             substrates: S ; products: P ; oligomeric_state: 2
             steps: begin
                 E + P ⇌ EP :: EqualAI
                 E + S ⇌ ES :: OnlyA
                 ES <--> EP :: OnlyA
             end
         end),
        # P, S, the chemistry
        (@_testhelper_unguarded @enzyme_mechanism begin
             substrates: S ; products: P ; oligomeric_state: 2
             steps: begin
                 E + P ⇌ EP :: OnlyA
                 E + S ⇌ ES :: OnlyA
                 ES <--> EP :: OnlyA
             end
         end)]
    regulated_expected = [
        # the chemistry, R1::OnlyA
        (@_testhelper_unguarded @enzyme_mechanism begin
             substrates: S ; products: P ; oligomeric_state: 2
             steps: begin
                 E + P ⇌ EP :: EqualAI
                 E + S ⇌ ES :: EqualAI
                 ES <--> EP :: OnlyA
             end
             regulatory_site(2): R1::OnlyA
         end),
        # the chemistry, R1::OnlyI
        (@_testhelper_unguarded @enzyme_mechanism begin
             substrates: S ; products: P ; oligomeric_state: 2
             steps: begin
                 E + P ⇌ EP :: EqualAI
                 E + S ⇌ ES :: EqualAI
                 ES <--> EP :: OnlyA
             end
             regulatory_site(2): R1::OnlyI
         end),
        # the chemistry, R2::OnlyA
        (@_testhelper_unguarded @enzyme_mechanism begin
             substrates: S ; products: P ; oligomeric_state: 2
             steps: begin
                 E + P ⇌ EP :: EqualAI
                 E + S ⇌ ES :: EqualAI
                 ES <--> EP :: OnlyA
             end
             regulatory_site(2): R2::OnlyA
         end),
        # the chemistry, R2::OnlyI
        (@_testhelper_unguarded @enzyme_mechanism begin
             substrates: S ; products: P ; oligomeric_state: 2
             steps: begin
                 E + P ⇌ EP :: EqualAI
                 E + S ⇌ ES :: EqualAI
                 ES <--> EP :: OnlyA
             end
             regulatory_site(2): R2::OnlyI
         end)]
    @test length(outs) == 11
    @test Set(outs) == Set([only_a_expected; regulated_expected])
    @test all(am -> Δ(am) == 2, regulated_expected)
    @test all(am -> Δ(am) == 1, only_a_expected)

    # A regulator already bound as a competitive inhibitor in the seed is not used at a
    # site; only the other (free) allosteric regulator is. In uni_uni_reg_and_inhibitor,
    # R1 is the competitive inhibitor and R2 the free allosteric regulator; the seed binds
    # R1 at free E. Four groups: the 15 nonempty sets :OnlyA, R1's binding among them,
    # then R1's group :OnlyI with each of the 8 sets of the other three groups :OnlyA,
    # the empty set among them, then the chemistry alone with R2 at a site.
    seed_ci = @enzyme_mechanism begin
        substrates: S; products: P
        modifiers: R1
        steps: begin
            E + P ⇌ EP
            E + R1 ⇌ ER1
            E + S ⇌ ES
            ES <--> EP
        end
    end
    @test seed_ci in _testhelper_children(
        uni_uni_reg_and_inhibitor,
        _testhelper_mechanisms(
            EnzymeRates.Enumeration.init_mechanisms(uni_uni_reg_and_inhibitor)))
    kids_ci =
        _testhelper_move(EnzymeRates.Enumeration._expand_to_allosteric, seed_ci,
                         uni_uni_reg_and_inhibitor)
    expected_ci = [
        # P
        (@_testhelper_unguarded @enzyme_mechanism begin
             substrates: S ; products: P ; modifiers: R1 ; oligomeric_state: 2
             steps: begin
                 E + P ⇌ EP   :: OnlyA
                 E + R1 ⇌ ER1 :: EqualAI
                 E + S ⇌ ES   :: EqualAI
                 ES <--> EP   :: EqualAI
             end
         end),
        # R1
        (@_testhelper_unguarded @enzyme_mechanism begin
             substrates: S ; products: P ; modifiers: R1 ; oligomeric_state: 2
             steps: begin
                 E + P ⇌ EP   :: EqualAI
                 E + R1 ⇌ ER1 :: OnlyA
                 E + S ⇌ ES   :: EqualAI
                 ES <--> EP   :: EqualAI
             end
         end),
        # P, R1
        (@_testhelper_unguarded @enzyme_mechanism begin
             substrates: S ; products: P ; modifiers: R1 ; oligomeric_state: 2
             steps: begin
                 E + P ⇌ EP   :: OnlyA
                 E + R1 ⇌ ER1 :: OnlyA
                 E + S ⇌ ES   :: EqualAI
                 ES <--> EP   :: EqualAI
             end
         end),
        # S
        (@_testhelper_unguarded @enzyme_mechanism begin
             substrates: S ; products: P ; modifiers: R1 ; oligomeric_state: 2
             steps: begin
                 E + P ⇌ EP   :: EqualAI
                 E + R1 ⇌ ER1 :: EqualAI
                 E + S ⇌ ES   :: OnlyA
                 ES <--> EP   :: EqualAI
             end
         end),
        # P, S
        (@_testhelper_unguarded @enzyme_mechanism begin
             substrates: S ; products: P ; modifiers: R1 ; oligomeric_state: 2
             steps: begin
                 E + P ⇌ EP   :: OnlyA
                 E + R1 ⇌ ER1 :: EqualAI
                 E + S ⇌ ES   :: OnlyA
                 ES <--> EP   :: EqualAI
             end
         end),
        # R1, S
        (@_testhelper_unguarded @enzyme_mechanism begin
             substrates: S ; products: P ; modifiers: R1 ; oligomeric_state: 2
             steps: begin
                 E + P ⇌ EP   :: EqualAI
                 E + R1 ⇌ ER1 :: OnlyA
                 E + S ⇌ ES   :: OnlyA
                 ES <--> EP   :: EqualAI
             end
         end),
        # P, R1, S
        (@_testhelper_unguarded @enzyme_mechanism begin
             substrates: S ; products: P ; modifiers: R1 ; oligomeric_state: 2
             steps: begin
                 E + P ⇌ EP   :: OnlyA
                 E + R1 ⇌ ER1 :: OnlyA
                 E + S ⇌ ES   :: OnlyA
                 ES <--> EP   :: EqualAI
             end
         end),
        # the chemistry
        (@_testhelper_unguarded @enzyme_mechanism begin
             substrates: S ; products: P ; modifiers: R1 ; oligomeric_state: 2
             steps: begin
                 E + P ⇌ EP   :: EqualAI
                 E + R1 ⇌ ER1 :: EqualAI
                 E + S ⇌ ES   :: EqualAI
                 ES <--> EP   :: OnlyA
             end
         end),
        # P, the chemistry
        (@_testhelper_unguarded @enzyme_mechanism begin
             substrates: S ; products: P ; modifiers: R1 ; oligomeric_state: 2
             steps: begin
                 E + P ⇌ EP   :: OnlyA
                 E + R1 ⇌ ER1 :: EqualAI
                 E + S ⇌ ES   :: EqualAI
                 ES <--> EP   :: OnlyA
             end
         end),
        # R1, the chemistry
        (@_testhelper_unguarded @enzyme_mechanism begin
             substrates: S ; products: P ; modifiers: R1 ; oligomeric_state: 2
             steps: begin
                 E + P ⇌ EP   :: EqualAI
                 E + R1 ⇌ ER1 :: OnlyA
                 E + S ⇌ ES   :: EqualAI
                 ES <--> EP   :: OnlyA
             end
         end),
        # P, R1, the chemistry
        (@_testhelper_unguarded @enzyme_mechanism begin
             substrates: S ; products: P ; modifiers: R1 ; oligomeric_state: 2
             steps: begin
                 E + P ⇌ EP   :: OnlyA
                 E + R1 ⇌ ER1 :: OnlyA
                 E + S ⇌ ES   :: EqualAI
                 ES <--> EP   :: OnlyA
             end
         end),
        # S, the chemistry
        (@_testhelper_unguarded @enzyme_mechanism begin
             substrates: S ; products: P ; modifiers: R1 ; oligomeric_state: 2
             steps: begin
                 E + P ⇌ EP   :: EqualAI
                 E + R1 ⇌ ER1 :: EqualAI
                 E + S ⇌ ES   :: OnlyA
                 ES <--> EP   :: OnlyA
             end
         end),
        # P, S, the chemistry
        (@_testhelper_unguarded @enzyme_mechanism begin
             substrates: S ; products: P ; modifiers: R1 ; oligomeric_state: 2
             steps: begin
                 E + P ⇌ EP   :: OnlyA
                 E + R1 ⇌ ER1 :: EqualAI
                 E + S ⇌ ES   :: OnlyA
                 ES <--> EP   :: OnlyA
             end
         end),
        # R1, S, the chemistry
        (@_testhelper_unguarded @enzyme_mechanism begin
             substrates: S ; products: P ; modifiers: R1 ; oligomeric_state: 2
             steps: begin
                 E + P ⇌ EP   :: EqualAI
                 E + R1 ⇌ ER1 :: OnlyA
                 E + S ⇌ ES   :: OnlyA
                 ES <--> EP   :: OnlyA
             end
         end),
        # P, R1, S, the chemistry
        (@_testhelper_unguarded @enzyme_mechanism begin
             substrates: S ; products: P ; modifiers: R1 ; oligomeric_state: 2
             steps: begin
                 E + P ⇌ EP   :: OnlyA
                 E + R1 ⇌ ER1 :: OnlyA
                 E + S ⇌ ES   :: OnlyA
                 ES <--> EP   :: OnlyA
             end
         end),
        # R1 :OnlyI
        (@_testhelper_unguarded @enzyme_mechanism begin
             substrates: S ; products: P ; modifiers: R1 ; oligomeric_state: 2
             steps: begin
                 E + P ⇌ EP   :: EqualAI
                 E + R1 ⇌ ER1 :: OnlyI
                 E + S ⇌ ES   :: EqualAI
                 ES <--> EP   :: EqualAI
             end
         end),
        # P; R1 :OnlyI
        (@_testhelper_unguarded @enzyme_mechanism begin
             substrates: S ; products: P ; modifiers: R1 ; oligomeric_state: 2
             steps: begin
                 E + P ⇌ EP   :: OnlyA
                 E + R1 ⇌ ER1 :: OnlyI
                 E + S ⇌ ES   :: EqualAI
                 ES <--> EP   :: EqualAI
             end
         end),
        # S; R1 :OnlyI
        (@_testhelper_unguarded @enzyme_mechanism begin
             substrates: S ; products: P ; modifiers: R1 ; oligomeric_state: 2
             steps: begin
                 E + P ⇌ EP   :: EqualAI
                 E + R1 ⇌ ER1 :: OnlyI
                 E + S ⇌ ES   :: OnlyA
                 ES <--> EP   :: EqualAI
             end
         end),
        # P, S; R1 :OnlyI
        (@_testhelper_unguarded @enzyme_mechanism begin
             substrates: S ; products: P ; modifiers: R1 ; oligomeric_state: 2
             steps: begin
                 E + P ⇌ EP   :: OnlyA
                 E + R1 ⇌ ER1 :: OnlyI
                 E + S ⇌ ES   :: OnlyA
                 ES <--> EP   :: EqualAI
             end
         end),
        # the chemistry; R1 :OnlyI
        (@_testhelper_unguarded @enzyme_mechanism begin
             substrates: S ; products: P ; modifiers: R1 ; oligomeric_state: 2
             steps: begin
                 E + P ⇌ EP   :: EqualAI
                 E + R1 ⇌ ER1 :: OnlyI
                 E + S ⇌ ES   :: EqualAI
                 ES <--> EP   :: OnlyA
             end
         end),
        # P, the chemistry; R1 :OnlyI
        (@_testhelper_unguarded @enzyme_mechanism begin
             substrates: S ; products: P ; modifiers: R1 ; oligomeric_state: 2
             steps: begin
                 E + P ⇌ EP   :: OnlyA
                 E + R1 ⇌ ER1 :: OnlyI
                 E + S ⇌ ES   :: EqualAI
                 ES <--> EP   :: OnlyA
             end
         end),
        # S, the chemistry; R1 :OnlyI
        (@_testhelper_unguarded @enzyme_mechanism begin
             substrates: S ; products: P ; modifiers: R1 ; oligomeric_state: 2
             steps: begin
                 E + P ⇌ EP   :: EqualAI
                 E + R1 ⇌ ER1 :: OnlyI
                 E + S ⇌ ES   :: OnlyA
                 ES <--> EP   :: OnlyA
             end
         end),
        # P, S, the chemistry; R1 :OnlyI
        (@_testhelper_unguarded @enzyme_mechanism begin
             substrates: S ; products: P ; modifiers: R1 ; oligomeric_state: 2
             steps: begin
                 E + P ⇌ EP   :: OnlyA
                 E + R1 ⇌ ER1 :: OnlyI
                 E + S ⇌ ES   :: OnlyA
                 ES <--> EP   :: OnlyA
             end
         end),
        # the chemistry, R2::OnlyA
        (@_testhelper_unguarded @enzyme_mechanism begin
             substrates: S ; products: P ; modifiers: R1 ; oligomeric_state: 2
             steps: begin
                 E + P ⇌ EP   :: EqualAI
                 E + R1 ⇌ ER1 :: EqualAI
                 E + S ⇌ ES   :: EqualAI
                 ES <--> EP   :: OnlyA
             end
             regulatory_site(2): R2::OnlyA
         end),
        # the chemistry, R2::OnlyI
        (@_testhelper_unguarded @enzyme_mechanism begin
             substrates: S ; products: P ; modifiers: R1 ; oligomeric_state: 2
             steps: begin
                 E + P ⇌ EP   :: EqualAI
                 E + R1 ⇌ ER1 :: EqualAI
                 E + S ⇌ ES   :: EqualAI
                 ES <--> EP   :: OnlyA
             end
             regulatory_site(2): R2::OnlyI
         end)]
    @test length(kids_ci) == 25
    @test Set(kids_ci) == Set(expected_ci)
    site_regs = Set(l for am in kids_ci
                    for s in am.allosteric.regulatory_sites
                    for l in s.ligands)
    @test site_regs == Set([:R2])
end

@testset "Mechanism — non-hyperbolic catalytic scheme: no children" begin
    # Uni-bi with random SS product release carries P² and Q² in its own
    # denominator; a conformational equilibrium on top would stack a second
    # source of concentration powers, so the promotion emits nothing.
    m = @enzyme_mechanism begin
        substrates: S
        products: P, Q
        steps: begin
            E + S ⇌ ES
            ES <--> EPQ
            (E + P <--> EP, EQ + P <--> EPQ)
            (E + Q <--> EQ, EP + Q <--> EPQ)
        end
    end
    rxn = @enzyme_reaction begin
        substrates: S[AB]
        products: P[A], Q[B]
        oligomeric_state: 2
    end
    @test isempty(_testhelper_move(EnzymeRates.Enumeration._expand_to_allosteric, m, rxn))
end

@testset "Mechanism — non-hyperbolic catalytic scheme: multiplicity 1 only" begin
    # The same uni-bi with random SS product release, with multiplicities 1 and 2
    # allowed. One catalytic subunit adds no concentration power of its own, so the
    # multiplicity-1 children are emitted exactly as for a reaction allowing only
    # multiplicity 1; multiplicity 2 would stack a second source of powers and emits
    # nothing. Four groups: the 15 nonempty sets :OnlyA, and P's group and Q's group each
    # tagged :NonequalAI alone, two minimal gaining sets refused as :OnlyA (each comment
    # names the child's :OnlyA groups, or its :NonequalAI one).
    m = @enzyme_mechanism begin
        substrates: S
        products: P, Q
        steps: begin
            E + S ⇌ ES
            ES <--> EPQ
            (E + P <--> EP, EQ + P <--> EPQ)
            (E + Q <--> EQ, EP + Q <--> EPQ)
        end
    end
    rxn_1_2 = @enzyme_reaction begin
        substrates: S[AB]
        products: P[A], Q[B]
        oligomeric_state: (1, 2)
    end
    rxn_1 = @enzyme_reaction begin
        substrates: S[AB]
        products: P[A], Q[B]
        oligomeric_state: (1,)
    end
    # Untyped: inside `T[...]` the parser reads `begin` as an index keyword.
    expected = [
        # S
        (@_testhelper_unguarded @enzyme_mechanism begin
             substrates: S ; products: P, Q ; oligomeric_state: 1
             steps: begin
                 E + S ⇌ ES                       :: OnlyA
                 ES <--> EPQ                      :: EqualAI
                 (E + P <--> EP, EQ + P <--> EPQ) :: EqualAI
                 (E + Q <--> EQ, EP + Q <--> EPQ) :: EqualAI
             end
         end),
        # the chemistry
        (@_testhelper_unguarded @enzyme_mechanism begin
             substrates: S ; products: P, Q ; oligomeric_state: 1
             steps: begin
                 E + S ⇌ ES                       :: EqualAI
                 ES <--> EPQ                      :: OnlyA
                 (E + P <--> EP, EQ + P <--> EPQ) :: EqualAI
                 (E + Q <--> EQ, EP + Q <--> EPQ) :: EqualAI
             end
         end),
        # S, the chemistry
        (@_testhelper_unguarded @enzyme_mechanism begin
             substrates: S ; products: P, Q ; oligomeric_state: 1
             steps: begin
                 E + S ⇌ ES                       :: OnlyA
                 ES <--> EPQ                      :: OnlyA
                 (E + P <--> EP, EQ + P <--> EPQ) :: EqualAI
                 (E + Q <--> EQ, EP + Q <--> EPQ) :: EqualAI
             end
         end),
        # P
        (@_testhelper_unguarded @enzyme_mechanism begin
             substrates: S ; products: P, Q ; oligomeric_state: 1
             steps: begin
                 E + S ⇌ ES                       :: EqualAI
                 ES <--> EPQ                      :: EqualAI
                 (E + P <--> EP, EQ + P <--> EPQ) :: OnlyA
                 (E + Q <--> EQ, EP + Q <--> EPQ) :: EqualAI
             end
         end),
        # S, P
        (@_testhelper_unguarded @enzyme_mechanism begin
             substrates: S ; products: P, Q ; oligomeric_state: 1
             steps: begin
                 E + S ⇌ ES                       :: OnlyA
                 ES <--> EPQ                      :: EqualAI
                 (E + P <--> EP, EQ + P <--> EPQ) :: OnlyA
                 (E + Q <--> EQ, EP + Q <--> EPQ) :: EqualAI
             end
         end),
        # the chemistry, P
        (@_testhelper_unguarded @enzyme_mechanism begin
             substrates: S ; products: P, Q ; oligomeric_state: 1
             steps: begin
                 E + S ⇌ ES                       :: EqualAI
                 ES <--> EPQ                      :: OnlyA
                 (E + P <--> EP, EQ + P <--> EPQ) :: OnlyA
                 (E + Q <--> EQ, EP + Q <--> EPQ) :: EqualAI
             end
         end),
        # S, the chemistry, P
        (@_testhelper_unguarded @enzyme_mechanism begin
             substrates: S ; products: P, Q ; oligomeric_state: 1
             steps: begin
                 E + S ⇌ ES                       :: OnlyA
                 ES <--> EPQ                      :: OnlyA
                 (E + P <--> EP, EQ + P <--> EPQ) :: OnlyA
                 (E + Q <--> EQ, EP + Q <--> EPQ) :: EqualAI
             end
         end),
        # Q
        (@_testhelper_unguarded @enzyme_mechanism begin
             substrates: S ; products: P, Q ; oligomeric_state: 1
             steps: begin
                 E + S ⇌ ES                       :: EqualAI
                 ES <--> EPQ                      :: EqualAI
                 (E + P <--> EP, EQ + P <--> EPQ) :: EqualAI
                 (E + Q <--> EQ, EP + Q <--> EPQ) :: OnlyA
             end
         end),
        # S, Q
        (@_testhelper_unguarded @enzyme_mechanism begin
             substrates: S ; products: P, Q ; oligomeric_state: 1
             steps: begin
                 E + S ⇌ ES                       :: OnlyA
                 ES <--> EPQ                      :: EqualAI
                 (E + P <--> EP, EQ + P <--> EPQ) :: EqualAI
                 (E + Q <--> EQ, EP + Q <--> EPQ) :: OnlyA
             end
         end),
        # the chemistry, Q
        (@_testhelper_unguarded @enzyme_mechanism begin
             substrates: S ; products: P, Q ; oligomeric_state: 1
             steps: begin
                 E + S ⇌ ES                       :: EqualAI
                 ES <--> EPQ                      :: OnlyA
                 (E + P <--> EP, EQ + P <--> EPQ) :: EqualAI
                 (E + Q <--> EQ, EP + Q <--> EPQ) :: OnlyA
             end
         end),
        # S, the chemistry, Q
        (@_testhelper_unguarded @enzyme_mechanism begin
             substrates: S ; products: P, Q ; oligomeric_state: 1
             steps: begin
                 E + S ⇌ ES                       :: OnlyA
                 ES <--> EPQ                      :: OnlyA
                 (E + P <--> EP, EQ + P <--> EPQ) :: EqualAI
                 (E + Q <--> EQ, EP + Q <--> EPQ) :: OnlyA
             end
         end),
        # P, Q
        (@_testhelper_unguarded @enzyme_mechanism begin
             substrates: S ; products: P, Q ; oligomeric_state: 1
             steps: begin
                 E + S ⇌ ES                       :: EqualAI
                 ES <--> EPQ                      :: EqualAI
                 (E + P <--> EP, EQ + P <--> EPQ) :: OnlyA
                 (E + Q <--> EQ, EP + Q <--> EPQ) :: OnlyA
             end
         end),
        # S, P, Q
        (@_testhelper_unguarded @enzyme_mechanism begin
             substrates: S ; products: P, Q ; oligomeric_state: 1
             steps: begin
                 E + S ⇌ ES                       :: OnlyA
                 ES <--> EPQ                      :: EqualAI
                 (E + P <--> EP, EQ + P <--> EPQ) :: OnlyA
                 (E + Q <--> EQ, EP + Q <--> EPQ) :: OnlyA
             end
         end),
        # the chemistry, P, Q
        (@_testhelper_unguarded @enzyme_mechanism begin
             substrates: S ; products: P, Q ; oligomeric_state: 1
             steps: begin
                 E + S ⇌ ES                       :: EqualAI
                 ES <--> EPQ                      :: OnlyA
                 (E + P <--> EP, EQ + P <--> EPQ) :: OnlyA
                 (E + Q <--> EQ, EP + Q <--> EPQ) :: OnlyA
             end
         end),
        # S, the chemistry, P, Q
        (@_testhelper_unguarded @enzyme_mechanism begin
             substrates: S ; products: P, Q ; oligomeric_state: 1
             steps: begin
                 E + S ⇌ ES                       :: OnlyA
                 ES <--> EPQ                      :: OnlyA
                 (E + P <--> EP, EQ + P <--> EPQ) :: OnlyA
                 (E + Q <--> EQ, EP + Q <--> EPQ) :: OnlyA
             end
         end),
        # P's group :NonequalAI
        (@_testhelper_unguarded @enzyme_mechanism begin
             substrates: S ; products: P, Q ; oligomeric_state: 1
             steps: begin
                 E + S ⇌ ES                       :: EqualAI
                 ES <--> EPQ                      :: EqualAI
                 (E + P <--> EP, EQ + P <--> EPQ) :: NonequalAI
                 (E + Q <--> EQ, EP + Q <--> EPQ) :: EqualAI
             end
         end),
        # Q's group :NonequalAI
        (@_testhelper_unguarded @enzyme_mechanism begin
             substrates: S ; products: P, Q ; oligomeric_state: 1
             steps: begin
                 E + S ⇌ ES                       :: EqualAI
                 ES <--> EPQ                      :: EqualAI
                 (E + P <--> EP, EQ + P <--> EPQ) :: EqualAI
                 (E + Q <--> EQ, EP + Q <--> EPQ) :: NonequalAI
             end
         end)]
    children = _testhelper_move(EnzymeRates.Enumeration._expand_to_allosteric, m, rxn_1_2)
    @test length(children) == 17
    @test Set(children) == Set(expected)
    @test all(am -> am.allosteric.oligomeric_state == 1, children)
    @test Set(children) == Set(
        _testhelper_move(EnzymeRates.Enumeration._expand_to_allosteric, m, rxn_1))
end

@testset "Mechanism — catalytic-site abortive complex: no children" begin
    # Ordered SS bi-bi with A binding EQ as an abortive complex at its
    # catalytic site: the term rooted at {EQ, EAQ} holds E → EA, so the
    # scheme carries A², and a conformational equilibrium on top would stack a
    # second source of concentration powers. The promotion emits nothing.
    m = @enzyme_mechanism begin
        substrates: A, B
        products: P, Q
        steps: begin
            E + A <--> EA
            EA + B <--> EAB
            E + Q <--> EQ
            EQ + P <--> EPQ
            EAB <--> EPQ
            EQ + A ⇌ EAQ
        end
    end
    rxn = @enzyme_reaction begin
        substrates: A[C], B[N]
        products: P[C], Q[N]
        oligomeric_state: 2
    end
    @test isempty(_testhelper_move(EnzymeRates.Enumeration._expand_to_allosteric, m, rxn))
end

@testset "Mechanism — substrate declared as an inhibitor: promotion proceeds" begin
    # Ordered SS bi-bi with A also declared a competitive inhibitor: EAmod, where A binds
    # without reacting, is a form apart from the catalytic A-bound forms. The inhibitor's
    # own group binds E and EQ; Q's catalytic binding is mirrored onto the
    # inhibitor-bound form EAmod, sharing its kinetic group, exactly as
    # `_expand_add_dead_end_regulator` builds it.
    parent = @enzyme_mechanism begin
        substrates: A, B
        products: P, Q
        steps: begin
            E + A <--> EA
            EA + B <--> EAB
            (E + Q <--> EQ, EAmod + Q <--> EAmodQ)
            EQ + P <--> EPQ
            EAB <--> EPQ
            (E + A ⇌ EAmod, EQ + A ⇌ EAmodQ)
        end
    end
    rxn = @enzyme_reaction begin
        substrates: A[C], B[N]
        products: P[C], Q[N]
        competitive_inhibitors: A
        oligomeric_state: 2
    end
    # Six groups, the inhibitor's own among them, which takes :EqualAI, :OnlyA or :OnlyI
    # while the other five take :EqualAI or :OnlyA: 2^5 · 3 − 1 = 95 assignments, none
    # dropped, since the copy binds free E, which both conformations hold, and no complex
    # of it has a productive twin. Four steady-state groups each gain alone tagged
    # :NonequalAI, A's catalytic binding, B's, Q's and P's, whose :OnlyA versions are
    # refused. The reaction declares no allosteric regulator, so no child has a site. The
    # derived denominator carries A², but from the inhibitor site, so the catalytic scheme
    # is hyperbolic and the promotion proceeds.
    a_nonequal = @enzyme_mechanism begin
        substrates: A, B ; products: P, Q ; oligomeric_state: 2
        steps: begin
            E + A <--> EA                          :: NonequalAI
            EA + B <--> EAB                        :: EqualAI
            (E + Q <--> EQ, EAmod + Q <--> EAmodQ) :: EqualAI
            EQ + P <--> EPQ                        :: EqualAI
            EAB <--> EPQ                           :: EqualAI
            (E + A ⇌ EAmod, EQ + A ⇌ EAmodQ)       :: EqualAI
        end
    end
    b_nonequal = @enzyme_mechanism begin
        substrates: A, B ; products: P, Q ; oligomeric_state: 2
        steps: begin
            E + A <--> EA                          :: EqualAI
            EA + B <--> EAB                        :: NonequalAI
            (E + Q <--> EQ, EAmod + Q <--> EAmodQ) :: EqualAI
            EQ + P <--> EPQ                        :: EqualAI
            EAB <--> EPQ                           :: EqualAI
            (E + A ⇌ EAmod, EQ + A ⇌ EAmodQ)       :: EqualAI
        end
    end
    q_nonequal = @enzyme_mechanism begin
        substrates: A, B ; products: P, Q ; oligomeric_state: 2
        steps: begin
            E + A <--> EA                          :: EqualAI
            EA + B <--> EAB                        :: EqualAI
            (E + Q <--> EQ, EAmod + Q <--> EAmodQ) :: NonequalAI
            EQ + P <--> EPQ                        :: EqualAI
            EAB <--> EPQ                           :: EqualAI
            (E + A ⇌ EAmod, EQ + A ⇌ EAmodQ)       :: EqualAI
        end
    end
    p_nonequal = @enzyme_mechanism begin
        substrates: A, B ; products: P, Q ; oligomeric_state: 2
        steps: begin
            E + A <--> EA                          :: EqualAI
            EA + B <--> EAB                        :: EqualAI
            (E + Q <--> EQ, EAmod + Q <--> EAmodQ) :: EqualAI
            EQ + P <--> EPQ                        :: NonequalAI
            EAB <--> EPQ                           :: EqualAI
            (E + A ⇌ EAmod, EQ + A ⇌ EAmodQ)       :: EqualAI
        end
    end
    children = _testhelper_move(EnzymeRates.Enumeration._expand_to_allosteric, parent, rxn)
    @test length(children) == 95 + 4
    @test allunique(children)
    nonequal = filter(c -> :NonequalAI in c.allosteric.allosteric_states, children)
    @test Set(nonequal) == Set([a_nonequal, b_nonequal, q_nonequal, p_nonequal])
    # The other 95, distinct, are the 95 assignments: no site, some group not :EqualAI,
    # and :OnlyI on the inhibitor's group alone.
    copy_group =
        parent.catalytic.steps[_testhelper_edge(parent, "E + A ⇌ EAmod")].kinetic_group
    @test all(setdiff(children, nonequal)) do c
        states = c.allosteric.allosteric_states
        c.catalytic == parent.catalytic && c.allosteric.oligomeric_state == 2 &&
            isempty(c.allosteric.regulatory_sites) && any(!=(:EqualAI), states) &&
            issubset(states, (:OnlyA, :EqualAI, :OnlyI)) &&
            all(g -> states[g] !== :OnlyI || g == copy_group, eachindex(states))
    end
    @test all(EnzymeRates.Enumeration._hyperbolic_catalysis,
              _testhelper_candidate.(children))
end
end

@testset "to_allosteric emits its children in the order its docstring states" begin
    # A uni-uni with both bindings at steady state, at two subunits with R declared, so
    # that every kind of child appears. The groups are numbered P's binding 1, S's 2, the
    # chemistry 3. The move emits the nonempty sets of groups :OnlyA in `_subsets` order
    # of those numbers, then the chemistry alone with R :OnlyA and then :OnlyI, then the
    # minimal gaining sets that are refused tagged :OnlyA, tagged :NonequalAI by size and
    # then lexicographically: P's binding, then S's. Each steady-state binding alone gains
    # its rate constant in the inactive conformation and is refused tagged :OnlyA. A move
    # that emitted its children in hash order would break the vector.
    m = @enzyme_mechanism begin
        substrates: S
        products: P
        steps: begin
            E + P <--> EP
            E + S <--> ES
            ES <--> EP
        end
    end
    @test [EnzymeRates.Mechanisms.step_text(m.catalytic, only(ks))
           for ks in _testhelper_group_members(m.catalytic)] ==
          ["E + P <--> EP", "E + S <--> ES", "ES <--> EP"]
    # Each child's comment names its :OnlyA groups, or its :NonequalAI one.
    # Untyped: inside `T[...]` the parser reads `begin` as an index keyword.
    expected = [
        # P
        (@_testhelper_unguarded @enzyme_mechanism begin
             substrates: S ; products: P ; oligomeric_state: 2
             steps: begin
                 E + P <--> EP :: OnlyA
                 E + S <--> ES :: EqualAI
                 ES <--> EP    :: EqualAI
             end
         end),
        # S
        (@_testhelper_unguarded @enzyme_mechanism begin
             substrates: S ; products: P ; oligomeric_state: 2
             steps: begin
                 E + P <--> EP :: EqualAI
                 E + S <--> ES :: OnlyA
                 ES <--> EP    :: EqualAI
             end
         end),
        # P, S
        (@_testhelper_unguarded @enzyme_mechanism begin
             substrates: S ; products: P ; oligomeric_state: 2
             steps: begin
                 E + P <--> EP :: OnlyA
                 E + S <--> ES :: OnlyA
                 ES <--> EP    :: EqualAI
             end
         end),
        # the chemistry
        (@_testhelper_unguarded @enzyme_mechanism begin
             substrates: S ; products: P ; oligomeric_state: 2
             steps: begin
                 E + P <--> EP :: EqualAI
                 E + S <--> ES :: EqualAI
                 ES <--> EP    :: OnlyA
             end
         end),
        # P, the chemistry
        (@_testhelper_unguarded @enzyme_mechanism begin
             substrates: S ; products: P ; oligomeric_state: 2
             steps: begin
                 E + P <--> EP :: OnlyA
                 E + S <--> ES :: EqualAI
                 ES <--> EP    :: OnlyA
             end
         end),
        # S, the chemistry
        (@_testhelper_unguarded @enzyme_mechanism begin
             substrates: S ; products: P ; oligomeric_state: 2
             steps: begin
                 E + P <--> EP :: EqualAI
                 E + S <--> ES :: OnlyA
                 ES <--> EP    :: OnlyA
             end
         end),
        # P, S, the chemistry
        (@_testhelper_unguarded @enzyme_mechanism begin
             substrates: S ; products: P ; oligomeric_state: 2
             steps: begin
                 E + P <--> EP :: OnlyA
                 E + S <--> ES :: OnlyA
                 ES <--> EP    :: OnlyA
             end
         end),
        # the chemistry, R::OnlyA
        (@_testhelper_unguarded @enzyme_mechanism begin
             substrates: S ; products: P ; oligomeric_state: 2
             steps: begin
                 E + P <--> EP :: EqualAI
                 E + S <--> ES :: EqualAI
                 ES <--> EP    :: OnlyA
             end
             regulatory_site(2): R::OnlyA
         end),
        # the chemistry, R::OnlyI
        (@_testhelper_unguarded @enzyme_mechanism begin
             substrates: S ; products: P ; oligomeric_state: 2
             steps: begin
                 E + P <--> EP :: EqualAI
                 E + S <--> ES :: EqualAI
                 ES <--> EP    :: OnlyA
             end
             regulatory_site(2): R::OnlyI
         end),
        # P :NonequalAI
        (@_testhelper_unguarded @enzyme_mechanism begin
             substrates: S ; products: P ; oligomeric_state: 2
             steps: begin
                 E + P <--> EP :: NonequalAI
                 E + S <--> ES :: EqualAI
                 ES <--> EP    :: EqualAI
             end
         end),
        # S :NonequalAI
        (@_testhelper_unguarded @enzyme_mechanism begin
             substrates: S ; products: P ; oligomeric_state: 2
             steps: begin
                 E + P <--> EP :: EqualAI
                 E + S <--> ES :: NonequalAI
                 ES <--> EP    :: EqualAI
             end
         end)]
    kids = _testhelper_move(EnzymeRates.Enumeration._expand_to_allosteric, m,
                            uni_uni_allo_reg)
    @test kids == expected
end

@testset "to_allosteric gives an inhibitor's group three states" begin
    # The uni-uni mechanism with the competitive inhibitor I bound at free E: four groups,
    # numbered I's binding 1, P's 2, S's 3 and the chemistry 4. At one subunit with no
    # regulator I's group takes :EqualAI, :OnlyA or :OnlyI and every other group
    # :EqualAI or :OnlyA: 2^3 · 3 − 1 = 23 children, in counting order with I's group
    # changing fastest. None is dropped: I binds free E, which both conformations hold,
    # and its complex EI has no productive twin. Tagged :NonequalAI, three minimal sets
    # gain a constant beyond L, I's group, the chemistry, and P's binding with S's; each
    # is returned tagged :OnlyA, so no child is tagged :NonequalAI. The child with only
    # I's group :OnlyI is returned and fits one constant more than the plain mechanism.
    rxn = @enzyme_reaction begin
        substrates: S[C]
        products: P[C]
        competitive_inhibitors: I
    end
    m = @enzyme_mechanism begin
        substrates: S
        products: P
        modifiers: I
        steps: begin
            E + I ⇌ EI
            E + P ⇌ EP
            E + S ⇌ ES
            ES <--> EP
        end
    end
    @test [EnzymeRates.Mechanisms.step_text(m.catalytic, only(ks))
           for ks in _testhelper_group_members(m.catalytic)] ==
          ["E + I ⇌ EI", "E + P ⇌ EP", "E + S ⇌ ES", "ES <--> EP"]
    # Each child's comment names its groups that are not :EqualAI.
    # Untyped: inside `T[...]` the parser reads `begin` as an index keyword.
    expected = [
        # I :OnlyA
        (@_testhelper_unguarded @enzyme_mechanism begin
             substrates: S ; products: P ; modifiers: I ; oligomeric_state: 1
             steps: begin
                 E + I ⇌ EI :: OnlyA
                 E + P ⇌ EP :: EqualAI
                 E + S ⇌ ES :: EqualAI
                 ES <--> EP :: EqualAI
             end
         end),
        # I :OnlyI
        (@_testhelper_unguarded @enzyme_mechanism begin
             substrates: S ; products: P ; modifiers: I ; oligomeric_state: 1
             steps: begin
                 E + I ⇌ EI :: OnlyI
                 E + P ⇌ EP :: EqualAI
                 E + S ⇌ ES :: EqualAI
                 ES <--> EP :: EqualAI
             end
         end),
        # P :OnlyA
        (@_testhelper_unguarded @enzyme_mechanism begin
             substrates: S ; products: P ; modifiers: I ; oligomeric_state: 1
             steps: begin
                 E + I ⇌ EI :: EqualAI
                 E + P ⇌ EP :: OnlyA
                 E + S ⇌ ES :: EqualAI
                 ES <--> EP :: EqualAI
             end
         end),
        # I, P :OnlyA
        (@_testhelper_unguarded @enzyme_mechanism begin
             substrates: S ; products: P ; modifiers: I ; oligomeric_state: 1
             steps: begin
                 E + I ⇌ EI :: OnlyA
                 E + P ⇌ EP :: OnlyA
                 E + S ⇌ ES :: EqualAI
                 ES <--> EP :: EqualAI
             end
         end),
        # P :OnlyA; I :OnlyI
        (@_testhelper_unguarded @enzyme_mechanism begin
             substrates: S ; products: P ; modifiers: I ; oligomeric_state: 1
             steps: begin
                 E + I ⇌ EI :: OnlyI
                 E + P ⇌ EP :: OnlyA
                 E + S ⇌ ES :: EqualAI
                 ES <--> EP :: EqualAI
             end
         end),
        # S :OnlyA
        (@_testhelper_unguarded @enzyme_mechanism begin
             substrates: S ; products: P ; modifiers: I ; oligomeric_state: 1
             steps: begin
                 E + I ⇌ EI :: EqualAI
                 E + P ⇌ EP :: EqualAI
                 E + S ⇌ ES :: OnlyA
                 ES <--> EP :: EqualAI
             end
         end),
        # I, S :OnlyA
        (@_testhelper_unguarded @enzyme_mechanism begin
             substrates: S ; products: P ; modifiers: I ; oligomeric_state: 1
             steps: begin
                 E + I ⇌ EI :: OnlyA
                 E + P ⇌ EP :: EqualAI
                 E + S ⇌ ES :: OnlyA
                 ES <--> EP :: EqualAI
             end
         end),
        # S :OnlyA; I :OnlyI
        (@_testhelper_unguarded @enzyme_mechanism begin
             substrates: S ; products: P ; modifiers: I ; oligomeric_state: 1
             steps: begin
                 E + I ⇌ EI :: OnlyI
                 E + P ⇌ EP :: EqualAI
                 E + S ⇌ ES :: OnlyA
                 ES <--> EP :: EqualAI
             end
         end),
        # P, S :OnlyA
        (@_testhelper_unguarded @enzyme_mechanism begin
             substrates: S ; products: P ; modifiers: I ; oligomeric_state: 1
             steps: begin
                 E + I ⇌ EI :: EqualAI
                 E + P ⇌ EP :: OnlyA
                 E + S ⇌ ES :: OnlyA
                 ES <--> EP :: EqualAI
             end
         end),
        # I, P, S :OnlyA
        (@_testhelper_unguarded @enzyme_mechanism begin
             substrates: S ; products: P ; modifiers: I ; oligomeric_state: 1
             steps: begin
                 E + I ⇌ EI :: OnlyA
                 E + P ⇌ EP :: OnlyA
                 E + S ⇌ ES :: OnlyA
                 ES <--> EP :: EqualAI
             end
         end),
        # P, S :OnlyA; I :OnlyI
        (@_testhelper_unguarded @enzyme_mechanism begin
             substrates: S ; products: P ; modifiers: I ; oligomeric_state: 1
             steps: begin
                 E + I ⇌ EI :: OnlyI
                 E + P ⇌ EP :: OnlyA
                 E + S ⇌ ES :: OnlyA
                 ES <--> EP :: EqualAI
             end
         end),
        # the chemistry :OnlyA
        (@_testhelper_unguarded @enzyme_mechanism begin
             substrates: S ; products: P ; modifiers: I ; oligomeric_state: 1
             steps: begin
                 E + I ⇌ EI :: EqualAI
                 E + P ⇌ EP :: EqualAI
                 E + S ⇌ ES :: EqualAI
                 ES <--> EP :: OnlyA
             end
         end),
        # I, the chemistry :OnlyA
        (@_testhelper_unguarded @enzyme_mechanism begin
             substrates: S ; products: P ; modifiers: I ; oligomeric_state: 1
             steps: begin
                 E + I ⇌ EI :: OnlyA
                 E + P ⇌ EP :: EqualAI
                 E + S ⇌ ES :: EqualAI
                 ES <--> EP :: OnlyA
             end
         end),
        # the chemistry :OnlyA; I :OnlyI
        (@_testhelper_unguarded @enzyme_mechanism begin
             substrates: S ; products: P ; modifiers: I ; oligomeric_state: 1
             steps: begin
                 E + I ⇌ EI :: OnlyI
                 E + P ⇌ EP :: EqualAI
                 E + S ⇌ ES :: EqualAI
                 ES <--> EP :: OnlyA
             end
         end),
        # P, the chemistry :OnlyA
        (@_testhelper_unguarded @enzyme_mechanism begin
             substrates: S ; products: P ; modifiers: I ; oligomeric_state: 1
             steps: begin
                 E + I ⇌ EI :: EqualAI
                 E + P ⇌ EP :: OnlyA
                 E + S ⇌ ES :: EqualAI
                 ES <--> EP :: OnlyA
             end
         end),
        # I, P, the chemistry :OnlyA
        (@_testhelper_unguarded @enzyme_mechanism begin
             substrates: S ; products: P ; modifiers: I ; oligomeric_state: 1
             steps: begin
                 E + I ⇌ EI :: OnlyA
                 E + P ⇌ EP :: OnlyA
                 E + S ⇌ ES :: EqualAI
                 ES <--> EP :: OnlyA
             end
         end),
        # P, the chemistry :OnlyA; I :OnlyI
        (@_testhelper_unguarded @enzyme_mechanism begin
             substrates: S ; products: P ; modifiers: I ; oligomeric_state: 1
             steps: begin
                 E + I ⇌ EI :: OnlyI
                 E + P ⇌ EP :: OnlyA
                 E + S ⇌ ES :: EqualAI
                 ES <--> EP :: OnlyA
             end
         end),
        # S, the chemistry :OnlyA
        (@_testhelper_unguarded @enzyme_mechanism begin
             substrates: S ; products: P ; modifiers: I ; oligomeric_state: 1
             steps: begin
                 E + I ⇌ EI :: EqualAI
                 E + P ⇌ EP :: EqualAI
                 E + S ⇌ ES :: OnlyA
                 ES <--> EP :: OnlyA
             end
         end),
        # I, S, the chemistry :OnlyA
        (@_testhelper_unguarded @enzyme_mechanism begin
             substrates: S ; products: P ; modifiers: I ; oligomeric_state: 1
             steps: begin
                 E + I ⇌ EI :: OnlyA
                 E + P ⇌ EP :: EqualAI
                 E + S ⇌ ES :: OnlyA
                 ES <--> EP :: OnlyA
             end
         end),
        # S, the chemistry :OnlyA; I :OnlyI
        (@_testhelper_unguarded @enzyme_mechanism begin
             substrates: S ; products: P ; modifiers: I ; oligomeric_state: 1
             steps: begin
                 E + I ⇌ EI :: OnlyI
                 E + P ⇌ EP :: EqualAI
                 E + S ⇌ ES :: OnlyA
                 ES <--> EP :: OnlyA
             end
         end),
        # P, S, the chemistry :OnlyA
        (@_testhelper_unguarded @enzyme_mechanism begin
             substrates: S ; products: P ; modifiers: I ; oligomeric_state: 1
             steps: begin
                 E + I ⇌ EI :: EqualAI
                 E + P ⇌ EP :: OnlyA
                 E + S ⇌ ES :: OnlyA
                 ES <--> EP :: OnlyA
             end
         end),
        # I, P, S, the chemistry :OnlyA
        (@_testhelper_unguarded @enzyme_mechanism begin
             substrates: S ; products: P ; modifiers: I ; oligomeric_state: 1
             steps: begin
                 E + I ⇌ EI :: OnlyA
                 E + P ⇌ EP :: OnlyA
                 E + S ⇌ ES :: OnlyA
                 ES <--> EP :: OnlyA
             end
         end),
        # P, S, the chemistry :OnlyA; I :OnlyI
        (@_testhelper_unguarded @enzyme_mechanism begin
             substrates: S ; products: P ; modifiers: I ; oligomeric_state: 1
             steps: begin
                 E + I ⇌ EI :: OnlyI
                 E + P ⇌ EP :: OnlyA
                 E + S ⇌ ES :: OnlyA
                 ES <--> EP :: OnlyA
             end
         end)]
    kids = _testhelper_move(EnzymeRates.Enumeration._expand_to_allosteric, m, rxn)
    @test length(kids) == 23
    @test Set(kids) == Set(expected)
    @test kids == expected
    only_i = expected[2]
    @test only_i in _testhelper_children(rxn, [m])
    fitted = EnzymeRates.Derivation.fitted_parameter_count
    @test fitted(only_i) == fitted(m) + 1
    # The three minimal gaining sets, each returned tagged :OnlyA.
    for (only_a, nonequal) in (
            (expected[1],
             @_testhelper_unguarded(@enzyme_mechanism begin
                 substrates: S ; products: P ; modifiers: I ; oligomeric_state: 1
                 steps: begin
                     E + I ⇌ EI :: NonequalAI
                     E + P ⇌ EP :: EqualAI
                     E + S ⇌ ES :: EqualAI
                     ES <--> EP :: EqualAI
                 end
             end)),
            (expected[12],
             @_testhelper_unguarded(@enzyme_mechanism begin
                 substrates: S ; products: P ; modifiers: I ; oligomeric_state: 1
                 steps: begin
                     E + I ⇌ EI :: EqualAI
                     E + P ⇌ EP :: EqualAI
                     E + S ⇌ ES :: EqualAI
                     ES <--> EP :: NonequalAI
                 end
             end)),
            (expected[9],
             @_testhelper_unguarded(@enzyme_mechanism begin
                 substrates: S ; products: P ; modifiers: I ; oligomeric_state: 1
                 steps: begin
                     E + I ⇌ EI :: EqualAI
                     E + P ⇌ EP :: NonequalAI
                     E + S ⇌ ES :: NonequalAI
                     ES <--> EP :: EqualAI
                 end
             end)))
        @test fitted(nonequal) == fitted(m) + 2
        @test EnzymeRates.Enumeration._enumerator_returns(_testhelper_candidate(only_a))
    end
end

@testset "to_allosteric drops an :OnlyI inhibitor that binds nowhere" begin
    # The three-step ping-pong with the competitive inhibitor I bound at the covalent form
    # E_res_plusA_minusP alone, where B and P are free. Four groups, numbered Q's binding
    # 1, the first half-reaction 2, I's binding 3 and the second half-reaction 4. The
    # inactive conformation reaches the covalent form from free E through the first
    # half-reaction, or through Q's binding and the second half-reaction. So with the
    # first half-reaction :OnlyA and Q's binding or the second half-reaction :OnlyA, I
    # :OnlyI binds in neither conformation and adds no constant: the move drops those
    # three of the 2^3 · 3 − 1 = 23 assignments and keeps the same :OnlyA sets with I
    # :EqualAI and :OnlyA. Tagged :NonequalAI each group alone gains a constant beyond L
    # and is returned tagged :OnlyA, so no child is tagged :NonequalAI.
    rxn = @enzyme_reaction begin
        substrates: A[CX], B[N]
        products: P[C], Q[NX]
        competitive_inhibitors: I
    end
    m = @enzyme_mechanism begin
        substrates: A, B
        products: P, Q
        modifiers: I
        steps: begin
            E + Q <--> EQ
            E + A <--> E_res_plusA_minusP + P
            E_res_plusA_minusP + I ⇌ EI_res_plusA_minusP
            E_res_plusA_minusP + B <--> EQ
        end
    end
    @test [EnzymeRates.Mechanisms.step_text(m.catalytic, only(ks))
           for ks in _testhelper_group_members(m.catalytic)] ==
          ["E + Q <--> EQ", "E + A <--> E_res_plusA_minusP + P",
           "E_res_plusA_minusP + I ⇌ EI_res_plusA_minusP",
           "E_res_plusA_minusP + B <--> EQ"]
    # Each child's comment names its groups that are not :EqualAI, in the move's order.
    # Untyped: inside `T[...]` the parser reads `begin` as an index keyword.
    expected = [
        # Q :OnlyA
        (@_testhelper_unguarded @enzyme_mechanism begin
             substrates: A, B ; products: P, Q ; modifiers: I ; oligomeric_state: 1
             steps: begin
                 E + Q <--> EQ                                :: OnlyA
                 E + A <--> E_res_plusA_minusP + P            :: EqualAI
                 E_res_plusA_minusP + I ⇌ EI_res_plusA_minusP :: EqualAI
                 E_res_plusA_minusP + B <--> EQ               :: EqualAI
             end
         end),
        # the first half-reaction :OnlyA
        (@_testhelper_unguarded @enzyme_mechanism begin
             substrates: A, B ; products: P, Q ; modifiers: I ; oligomeric_state: 1
             steps: begin
                 E + Q <--> EQ                                :: EqualAI
                 E + A <--> E_res_plusA_minusP + P            :: OnlyA
                 E_res_plusA_minusP + I ⇌ EI_res_plusA_minusP :: EqualAI
                 E_res_plusA_minusP + B <--> EQ               :: EqualAI
             end
         end),
        # Q, the first half-reaction :OnlyA
        (@_testhelper_unguarded @enzyme_mechanism begin
             substrates: A, B ; products: P, Q ; modifiers: I ; oligomeric_state: 1
             steps: begin
                 E + Q <--> EQ                                :: OnlyA
                 E + A <--> E_res_plusA_minusP + P            :: OnlyA
                 E_res_plusA_minusP + I ⇌ EI_res_plusA_minusP :: EqualAI
                 E_res_plusA_minusP + B <--> EQ               :: EqualAI
             end
         end),
        # I :OnlyA
        (@_testhelper_unguarded @enzyme_mechanism begin
             substrates: A, B ; products: P, Q ; modifiers: I ; oligomeric_state: 1
             steps: begin
                 E + Q <--> EQ                                :: EqualAI
                 E + A <--> E_res_plusA_minusP + P            :: EqualAI
                 E_res_plusA_minusP + I ⇌ EI_res_plusA_minusP :: OnlyA
                 E_res_plusA_minusP + B <--> EQ               :: EqualAI
             end
         end),
        # Q, I :OnlyA
        (@_testhelper_unguarded @enzyme_mechanism begin
             substrates: A, B ; products: P, Q ; modifiers: I ; oligomeric_state: 1
             steps: begin
                 E + Q <--> EQ                                :: OnlyA
                 E + A <--> E_res_plusA_minusP + P            :: EqualAI
                 E_res_plusA_minusP + I ⇌ EI_res_plusA_minusP :: OnlyA
                 E_res_plusA_minusP + B <--> EQ               :: EqualAI
             end
         end),
        # the first half-reaction, I :OnlyA
        (@_testhelper_unguarded @enzyme_mechanism begin
             substrates: A, B ; products: P, Q ; modifiers: I ; oligomeric_state: 1
             steps: begin
                 E + Q <--> EQ                                :: EqualAI
                 E + A <--> E_res_plusA_minusP + P            :: OnlyA
                 E_res_plusA_minusP + I ⇌ EI_res_plusA_minusP :: OnlyA
                 E_res_plusA_minusP + B <--> EQ               :: EqualAI
             end
         end),
        # Q, the first half-reaction, I :OnlyA
        (@_testhelper_unguarded @enzyme_mechanism begin
             substrates: A, B ; products: P, Q ; modifiers: I ; oligomeric_state: 1
             steps: begin
                 E + Q <--> EQ                                :: OnlyA
                 E + A <--> E_res_plusA_minusP + P            :: OnlyA
                 E_res_plusA_minusP + I ⇌ EI_res_plusA_minusP :: OnlyA
                 E_res_plusA_minusP + B <--> EQ               :: EqualAI
             end
         end),
        # I :OnlyI
        (@_testhelper_unguarded @enzyme_mechanism begin
             substrates: A, B ; products: P, Q ; modifiers: I ; oligomeric_state: 1
             steps: begin
                 E + Q <--> EQ                                :: EqualAI
                 E + A <--> E_res_plusA_minusP + P            :: EqualAI
                 E_res_plusA_minusP + I ⇌ EI_res_plusA_minusP :: OnlyI
                 E_res_plusA_minusP + B <--> EQ               :: EqualAI
             end
         end),
        # Q :OnlyA; I :OnlyI
        (@_testhelper_unguarded @enzyme_mechanism begin
             substrates: A, B ; products: P, Q ; modifiers: I ; oligomeric_state: 1
             steps: begin
                 E + Q <--> EQ                                :: OnlyA
                 E + A <--> E_res_plusA_minusP + P            :: EqualAI
                 E_res_plusA_minusP + I ⇌ EI_res_plusA_minusP :: OnlyI
                 E_res_plusA_minusP + B <--> EQ               :: EqualAI
             end
         end),
        # the first half-reaction :OnlyA; I :OnlyI
        (@_testhelper_unguarded @enzyme_mechanism begin
             substrates: A, B ; products: P, Q ; modifiers: I ; oligomeric_state: 1
             steps: begin
                 E + Q <--> EQ                                :: EqualAI
                 E + A <--> E_res_plusA_minusP + P            :: OnlyA
                 E_res_plusA_minusP + I ⇌ EI_res_plusA_minusP :: OnlyI
                 E_res_plusA_minusP + B <--> EQ               :: EqualAI
             end
         end),
        # the second half-reaction :OnlyA
        (@_testhelper_unguarded @enzyme_mechanism begin
             substrates: A, B ; products: P, Q ; modifiers: I ; oligomeric_state: 1
             steps: begin
                 E + Q <--> EQ                                :: EqualAI
                 E + A <--> E_res_plusA_minusP + P            :: EqualAI
                 E_res_plusA_minusP + I ⇌ EI_res_plusA_minusP :: EqualAI
                 E_res_plusA_minusP + B <--> EQ               :: OnlyA
             end
         end),
        # Q, the second half-reaction :OnlyA
        (@_testhelper_unguarded @enzyme_mechanism begin
             substrates: A, B ; products: P, Q ; modifiers: I ; oligomeric_state: 1
             steps: begin
                 E + Q <--> EQ                                :: OnlyA
                 E + A <--> E_res_plusA_minusP + P            :: EqualAI
                 E_res_plusA_minusP + I ⇌ EI_res_plusA_minusP :: EqualAI
                 E_res_plusA_minusP + B <--> EQ               :: OnlyA
             end
         end),
        # the first half-reaction, the second half-reaction :OnlyA
        (@_testhelper_unguarded @enzyme_mechanism begin
             substrates: A, B ; products: P, Q ; modifiers: I ; oligomeric_state: 1
             steps: begin
                 E + Q <--> EQ                                :: EqualAI
                 E + A <--> E_res_plusA_minusP + P            :: OnlyA
                 E_res_plusA_minusP + I ⇌ EI_res_plusA_minusP :: EqualAI
                 E_res_plusA_minusP + B <--> EQ               :: OnlyA
             end
         end),
        # Q, the first half-reaction, the second half-reaction :OnlyA
        (@_testhelper_unguarded @enzyme_mechanism begin
             substrates: A, B ; products: P, Q ; modifiers: I ; oligomeric_state: 1
             steps: begin
                 E + Q <--> EQ                                :: OnlyA
                 E + A <--> E_res_plusA_minusP + P            :: OnlyA
                 E_res_plusA_minusP + I ⇌ EI_res_plusA_minusP :: EqualAI
                 E_res_plusA_minusP + B <--> EQ               :: OnlyA
             end
         end),
        # I, the second half-reaction :OnlyA
        (@_testhelper_unguarded @enzyme_mechanism begin
             substrates: A, B ; products: P, Q ; modifiers: I ; oligomeric_state: 1
             steps: begin
                 E + Q <--> EQ                                :: EqualAI
                 E + A <--> E_res_plusA_minusP + P            :: EqualAI
                 E_res_plusA_minusP + I ⇌ EI_res_plusA_minusP :: OnlyA
                 E_res_plusA_minusP + B <--> EQ               :: OnlyA
             end
         end),
        # Q, I, the second half-reaction :OnlyA
        (@_testhelper_unguarded @enzyme_mechanism begin
             substrates: A, B ; products: P, Q ; modifiers: I ; oligomeric_state: 1
             steps: begin
                 E + Q <--> EQ                                :: OnlyA
                 E + A <--> E_res_plusA_minusP + P            :: EqualAI
                 E_res_plusA_minusP + I ⇌ EI_res_plusA_minusP :: OnlyA
                 E_res_plusA_minusP + B <--> EQ               :: OnlyA
             end
         end),
        # the first half-reaction, I, the second half-reaction :OnlyA
        (@_testhelper_unguarded @enzyme_mechanism begin
             substrates: A, B ; products: P, Q ; modifiers: I ; oligomeric_state: 1
             steps: begin
                 E + Q <--> EQ                                :: EqualAI
                 E + A <--> E_res_plusA_minusP + P            :: OnlyA
                 E_res_plusA_minusP + I ⇌ EI_res_plusA_minusP :: OnlyA
                 E_res_plusA_minusP + B <--> EQ               :: OnlyA
             end
         end),
        # Q, the first half-reaction, I, the second half-reaction :OnlyA
        (@_testhelper_unguarded @enzyme_mechanism begin
             substrates: A, B ; products: P, Q ; modifiers: I ; oligomeric_state: 1
             steps: begin
                 E + Q <--> EQ                                :: OnlyA
                 E + A <--> E_res_plusA_minusP + P            :: OnlyA
                 E_res_plusA_minusP + I ⇌ EI_res_plusA_minusP :: OnlyA
                 E_res_plusA_minusP + B <--> EQ               :: OnlyA
             end
         end),
        # the second half-reaction :OnlyA; I :OnlyI
        (@_testhelper_unguarded @enzyme_mechanism begin
             substrates: A, B ; products: P, Q ; modifiers: I ; oligomeric_state: 1
             steps: begin
                 E + Q <--> EQ                                :: EqualAI
                 E + A <--> E_res_plusA_minusP + P            :: EqualAI
                 E_res_plusA_minusP + I ⇌ EI_res_plusA_minusP :: OnlyI
                 E_res_plusA_minusP + B <--> EQ               :: OnlyA
             end
         end),
        # Q, the second half-reaction :OnlyA; I :OnlyI
        (@_testhelper_unguarded @enzyme_mechanism begin
             substrates: A, B ; products: P, Q ; modifiers: I ; oligomeric_state: 1
             steps: begin
                 E + Q <--> EQ                                :: OnlyA
                 E + A <--> E_res_plusA_minusP + P            :: EqualAI
                 E_res_plusA_minusP + I ⇌ EI_res_plusA_minusP :: OnlyI
                 E_res_plusA_minusP + B <--> EQ               :: OnlyA
             end
         end)]
    nowhere = [
        # Q, the first half-reaction :OnlyA; I :OnlyI
        (@_testhelper_unguarded @enzyme_mechanism begin
             substrates: A, B ; products: P, Q ; modifiers: I ; oligomeric_state: 1
             steps: begin
                 E + Q <--> EQ                                :: OnlyA
                 E + A <--> E_res_plusA_minusP + P            :: OnlyA
                 E_res_plusA_minusP + I ⇌ EI_res_plusA_minusP :: OnlyI
                 E_res_plusA_minusP + B <--> EQ               :: EqualAI
             end
         end),
        # the first half-reaction, the second half-reaction :OnlyA; I :OnlyI
        (@_testhelper_unguarded @enzyme_mechanism begin
             substrates: A, B ; products: P, Q ; modifiers: I ; oligomeric_state: 1
             steps: begin
                 E + Q <--> EQ                                :: EqualAI
                 E + A <--> E_res_plusA_minusP + P            :: OnlyA
                 E_res_plusA_minusP + I ⇌ EI_res_plusA_minusP :: OnlyI
                 E_res_plusA_minusP + B <--> EQ               :: OnlyA
             end
         end),
        # Q, the first half-reaction, the second half-reaction :OnlyA; I :OnlyI
        (@_testhelper_unguarded @enzyme_mechanism begin
             substrates: A, B ; products: P, Q ; modifiers: I ; oligomeric_state: 1
             steps: begin
                 E + Q <--> EQ                                :: OnlyA
                 E + A <--> E_res_plusA_minusP + P            :: OnlyA
                 E_res_plusA_minusP + I ⇌ EI_res_plusA_minusP :: OnlyI
                 E_res_plusA_minusP + B <--> EQ               :: OnlyA
             end
         end)]
    kids = _testhelper_move(EnzymeRates.Enumeration._expand_to_allosteric, m, rxn)
    @test length(kids) == 20
    @test Set(kids) == Set(expected)
    @test kids == expected
    fitted = EnzymeRates.Derivation.fitted_parameter_count
    i_binding = _testhelper_edge(m, "E_res_plusA_minusP + I ⇌ EI_res_plusA_minusP")
    for am in nowhere
        runs = EnzymeRates.Mechanisms.steps_running_in_allosteric_conformations(am)
        @test !runs.active[i_binding] && !runs.inactive[i_binding]
        @test !EnzymeRates.Enumeration._inhibitors_bind(_testhelper_candidate(am))
        @test fitted(am) == fitted(m)
        # The enumerator would return it, so only the move keeps it from the search.
        @test EnzymeRates.Enumeration._enumerator_returns(_testhelper_candidate(am))
    end
    @test isempty(intersect(nowhere, _testhelper_children(rxn, [m])))
end

@testset "to_allosteric drops a reactant copy an :OnlyI inhibitor makes redundant" begin
    # An ordered bi-bi whose B binding forms EPQ, with A bound a second time as a copy at
    # free E (EAmod) and the inhibitor I bound at EA and EQ. Six groups, numbered A's
    # binding 1, A's copy 2, Q's binding 3, I's binding 4, B's 5 and P's 6. The copy's
    # complex EAmod has the productive twin EA. Tag A's binding, the copy and B's binding
    # :OnlyA and I's group :OnlyI: the active conformation, the only one the copy binds,
    # then no longer binds I at EA, and the copy is redundant (`_redundant_copy_groups`),
    # so the move drops that child. With I's group :EqualAI or :OnlyA instead, the active
    # conformation binds I at EA, the copy is not redundant, and both children are
    # emitted. The two groups that bind an inhibitor take three states, the other four
    # two: 2^4 · 3^2 − 1 = 143 assignments. 12 are dropped because I binds in neither
    # conformation and 6 because the copy is redundant, which leaves 125, and Q's binding
    # tagged :NonequalAI is the one minimal gaining set refused tagged :OnlyA: 126
    # children. They are not all written out: the smallest mechanism that shows this drop
    # has these six groups and 126 children, and a search over the bi-bi, ping-pong and
    # uni-bi starting mechanisms with one or two inhibitors bound found none smaller.
    rxn = @enzyme_reaction begin
        substrates: A[C], B[N]
        products: P[C], Q[N]
        competitive_inhibitors: A, I
    end
    m = @enzyme_mechanism begin
        substrates: A, B
        products: P, Q
        modifiers: I
        steps: begin
            E + A ⇌ EA
            E + A ⇌ EAmod
            E + Q <--> EQ
            (EA + I ⇌ EAI, EQ + I ⇌ EIQ)
            EA + B <--> EPQ
            EQ + P ⇌ EPQ
        end
    end
    @test [EnzymeRates.Mechanisms.step_text(m.catalytic, first(ks))
           for ks in _testhelper_group_members(m.catalytic)] ==
          ["E + A ⇌ EA", "E + A ⇌ EAmod", "E + Q <--> EQ", "EA + I ⇌ EAI",
           "EA + B <--> EPQ", "EQ + P ⇌ EPQ"]
    redundant = @_testhelper_unguarded @enzyme_mechanism begin
        substrates: A, B ; products: P, Q ; modifiers: I ; oligomeric_state: 1
        steps: begin
            E + A ⇌ EA                   :: OnlyA
            E + A ⇌ EAmod                :: OnlyA
            E + Q <--> EQ                :: EqualAI
            (EA + I ⇌ EAI, EQ + I ⇌ EIQ) :: OnlyI
            EA + B <--> EPQ              :: OnlyA
            EQ + P ⇌ EPQ                 :: EqualAI
        end
    end
    i_equal = @_testhelper_unguarded @enzyme_mechanism begin
        substrates: A, B ; products: P, Q ; modifiers: I ; oligomeric_state: 1
        steps: begin
            E + A ⇌ EA                   :: OnlyA
            E + A ⇌ EAmod                :: OnlyA
            E + Q <--> EQ                :: EqualAI
            (EA + I ⇌ EAI, EQ + I ⇌ EIQ) :: EqualAI
            EA + B <--> EPQ              :: OnlyA
            EQ + P ⇌ EPQ                 :: EqualAI
        end
    end
    i_only_a = @_testhelper_unguarded @enzyme_mechanism begin
        substrates: A, B ; products: P, Q ; modifiers: I ; oligomeric_state: 1
        steps: begin
            E + A ⇌ EA                   :: OnlyA
            E + A ⇌ EAmod                :: OnlyA
            E + Q <--> EQ                :: EqualAI
            (EA + I ⇌ EAI, EQ + I ⇌ EIQ) :: OnlyA
            EA + B <--> EPQ              :: OnlyA
            EQ + P ⇌ EPQ                 :: EqualAI
        end
    end
    q_nonequal = @enzyme_mechanism begin
        substrates: A, B ; products: P, Q ; modifiers: I ; oligomeric_state: 1
        steps: begin
            E + A ⇌ EA                   :: EqualAI
            E + A ⇌ EAmod                :: EqualAI
            E + Q <--> EQ                :: NonequalAI
            (EA + I ⇌ EAI, EQ + I ⇌ EIQ) :: EqualAI
            EA + B <--> EPQ              :: EqualAI
            EQ + P ⇌ EPQ                 :: EqualAI
        end
    end
    redundant_groups = EnzymeRates.Enumeration._redundant_copy_groups
    @test redundant_groups(_testhelper_candidate(redundant)) == [2]
    @test isempty(redundant_groups(_testhelper_candidate(i_equal)))
    @test isempty(redundant_groups(_testhelper_candidate(i_only_a)))
    # I binds the inactive conformation, so the copy alone drops the child; and the
    # enumerator would return it, so only the move keeps it from the search.
    @test EnzymeRates.Enumeration._inhibitors_bind(_testhelper_candidate(redundant))
    @test EnzymeRates.Enumeration._enumerator_returns(_testhelper_candidate(redundant))
    children = _testhelper_move(EnzymeRates.Enumeration._expand_to_allosteric, m, rxn)
    @test length(children) == 126
    @test allunique(children)
    @test all(c -> isempty(redundant_groups(_testhelper_candidate(c))), children)
    @test redundant ∉ children
    @test i_equal in children
    @test i_only_a in children
    @test filter(c -> :NonequalAI in c.allosteric.allosteric_states, children) ==
          [q_nonequal]
end

@testset "no child holds :OnlyI without a competitive inhibitor" begin
    # The uni-uni mechanism binds no inhibitor, so each of its three groups takes
    # :EqualAI or :OnlyA alone: the seven sets of groups :OnlyA, at one subunit with no
    # regulator, whether or not the reaction declares the inhibitor I. No child holds an
    # :OnlyI group.
    m = @enzyme_mechanism begin
        substrates: S
        products: P
        steps: begin
            E + S ⇌ ES
            ES <--> EP
            E + P ⇌ EP
        end
    end
    without_i = @enzyme_reaction begin
        substrates: S[C]
        products: P[C]
    end
    with_i = @enzyme_reaction begin
        substrates: S[C]
        products: P[C]
        competitive_inhibitors: I
    end
    # Untyped: inside `T[...]` the parser reads `begin` as an index keyword.
    expected = [
        (@_testhelper_unguarded @enzyme_mechanism begin
             substrates: S ; products: P ; oligomeric_state: 1
             steps: begin
                 E + S ⇌ ES :: OnlyA
                 ES <--> EP :: EqualAI
                 E + P ⇌ EP :: EqualAI
             end
         end),
        (@_testhelper_unguarded @enzyme_mechanism begin
             substrates: S ; products: P ; oligomeric_state: 1
             steps: begin
                 E + S ⇌ ES :: EqualAI
                 ES <--> EP :: OnlyA
                 E + P ⇌ EP :: EqualAI
             end
         end),
        (@_testhelper_unguarded @enzyme_mechanism begin
             substrates: S ; products: P ; oligomeric_state: 1
             steps: begin
                 E + S ⇌ ES :: OnlyA
                 ES <--> EP :: OnlyA
                 E + P ⇌ EP :: EqualAI
             end
         end),
        (@_testhelper_unguarded @enzyme_mechanism begin
             substrates: S ; products: P ; oligomeric_state: 1
             steps: begin
                 E + S ⇌ ES :: EqualAI
                 ES <--> EP :: EqualAI
                 E + P ⇌ EP :: OnlyA
             end
         end),
        (@_testhelper_unguarded @enzyme_mechanism begin
             substrates: S ; products: P ; oligomeric_state: 1
             steps: begin
                 E + S ⇌ ES :: OnlyA
                 ES <--> EP :: EqualAI
                 E + P ⇌ EP :: OnlyA
             end
         end),
        (@_testhelper_unguarded @enzyme_mechanism begin
             substrates: S ; products: P ; oligomeric_state: 1
             steps: begin
                 E + S ⇌ ES :: EqualAI
                 ES <--> EP :: OnlyA
                 E + P ⇌ EP :: OnlyA
             end
         end),
        (@_testhelper_unguarded @enzyme_mechanism begin
             substrates: S ; products: P ; oligomeric_state: 1
             steps: begin
                 E + S ⇌ ES :: OnlyA
                 ES <--> EP :: OnlyA
                 E + P ⇌ EP :: OnlyA
             end
         end)]
    for rxn in (without_i, with_i)
        kids = _testhelper_move(EnzymeRates.Enumeration._expand_to_allosteric, m, rxn)
        @test length(kids) == 7
        @test Set(kids) == Set(expected)
        @test all(k -> :OnlyI ∉ k.allosteric.allosteric_states, kids)
    end
end

# ─── _expand_add_allosteric_regulator ──────────────────────────────────
@testset "_expand_add_allosteric_regulator" begin

@testset "Allosteric mechanism — a first allo regulator: 3 variants" begin
    # SEED: uni-uni allosteric with all groups :EqualAI and no allosteric regulator
    # added yet. A declared allosteric regulator x is the only un-added one, so it gets
    # one new site (0 existing reg sites) with each of the 3 non-:EqualAI tags. The
    # :EqualAI branch is gated to "existing site with ≥1 non-:EqualAI ligand", so it
    # does not apply.
    am = @enzyme_mechanism begin
        substrates: S; products: P
        oligomeric_state: 2
        steps: begin
            E + P ⇌ EP    :: EqualAI
            E + S ⇌ ES    :: EqualAI
            ES <--> EP  :: EqualAI
        end
    end
    _testhelper_assert_mechanism_invariants(am)

    # Each reaction with its expected children: x on a site of its own at each of the
    # three non-:EqualAI tags.
    for (rxn, x, expected) in [
        # R is declared only as an allosteric regulator.
        (@enzyme_reaction(begin
            substrates: S[C]
            products: P[C]
            allosteric_regulators: R
            oligomeric_state: 2
        end), :R, [
            (@enzyme_mechanism begin
                substrates: S; products: P
                oligomeric_state: 2
                steps: begin
                    E + P ⇌ EP    :: EqualAI
                    E + S ⇌ ES    :: EqualAI
                    ES <--> EP  :: EqualAI
                end
                regulatory_site(2): R::OnlyA
            end),
            (@enzyme_mechanism begin
                substrates: S; products: P
                oligomeric_state: 2
                steps: begin
                    E + P ⇌ EP    :: EqualAI
                    E + S ⇌ ES    :: EqualAI
                    ES <--> EP  :: EqualAI
                end
                regulatory_site(2): R::OnlyI
            end),
            (@enzyme_mechanism begin
                substrates: S; products: P
                oligomeric_state: 2
                steps: begin
                    E + P ⇌ EP    :: EqualAI
                    E + S ⇌ ES    :: EqualAI
                    ES <--> EP  :: EqualAI
                end
                regulatory_site(2): R::NonequalAI
            end),
        ]),
        # S is both the substrate and an allosteric regulator; it still plays its
        # catalytic substrate role in the base mechanism.
        (@enzyme_reaction(begin
            substrates: S[C]
            products: P[C]
            allosteric_regulators: S
            oligomeric_state: 2
        end), :S, [
            (@enzyme_mechanism begin
                substrates: S; products: P
                oligomeric_state: 2
                steps: begin
                    E + P ⇌ EP    :: EqualAI
                    E + S ⇌ ES    :: EqualAI
                    ES <--> EP  :: EqualAI
                end
                regulatory_site(2): S::OnlyA
            end),
            (@enzyme_mechanism begin
                substrates: S; products: P
                oligomeric_state: 2
                steps: begin
                    E + P ⇌ EP    :: EqualAI
                    E + S ⇌ ES    :: EqualAI
                    ES <--> EP  :: EqualAI
                end
                regulatory_site(2): S::OnlyI
            end),
            (@enzyme_mechanism begin
                substrates: S; products: P
                oligomeric_state: 2
                steps: begin
                    E + P ⇌ EP    :: EqualAI
                    E + S ⇌ ES    :: EqualAI
                    ES <--> EP  :: EqualAI
                end
                regulatory_site(2): S::NonequalAI
            end),
        ]),
        # P is both the product and an allosteric regulator.
        (@enzyme_reaction(begin
            substrates: S[C]
            products: P[C]
            allosteric_regulators: P
            oligomeric_state: 2
        end), :P, [
            (@enzyme_mechanism begin
                substrates: S; products: P
                oligomeric_state: 2
                steps: begin
                    E + P ⇌ EP    :: EqualAI
                    E + S ⇌ ES    :: EqualAI
                    ES <--> EP  :: EqualAI
                end
                regulatory_site(2): P::OnlyA
            end),
            (@enzyme_mechanism begin
                substrates: S; products: P
                oligomeric_state: 2
                steps: begin
                    E + P ⇌ EP    :: EqualAI
                    E + S ⇌ ES    :: EqualAI
                    ES <--> EP  :: EqualAI
                end
                regulatory_site(2): P::OnlyI
            end),
            (@enzyme_mechanism begin
                substrates: S; products: P
                oligomeric_state: 2
                steps: begin
                    E + P ⇌ EP    :: EqualAI
                    E + S ⇌ ES    :: EqualAI
                    ES <--> EP  :: EqualAI
                end
                regulatory_site(2): P::NonequalAI
            end),
        ]),
    ]
        result = _testhelper_move(EnzymeRates.Enumeration._expand_add_allosteric_regulator,
                                  am, rxn)

        # 1. count: 3 non-:EqualAI tags × 1 site option (new site only) = 3.
        @test length(result) == 3
        @test Set(result) == Set(expected)

        # 2. Δ params: :OnlyA/:OnlyI add one K_x each (+1); :NonequalAI adds K_x and
        # K_x_T (+2). Sorted [1, 1, 2]. Measured against the actual compiled fitted
        # count.
        @test _testhelper_param_deltas(am, result) == [1, 1, 2]

        # 3. each variant has exactly one regulatory site holding the single ligand x,
        # and the tags across the 3 variants are exactly {:OnlyA, :OnlyI, :NonequalAI}.
        # The added site's multiplicity is am.allosteric.oligomeric_state.
        @test all(r -> length(r.allosteric.regulatory_sites) == 1, result)
        @test Set((site.ligands,
                   collect(site.allosteric_states),
                   site.multiplicity)
                  for r in result for site in r.allosteric.regulatory_sites) == Set([
            ([x], [:OnlyA],      am.allosteric.oligomeric_state),
            ([x], [:OnlyI],      am.allosteric.oligomeric_state),
            ([x], [:NonequalAI], am.allosteric.oligomeric_state),
        ])

        # 4. invariants and preservation: catalytic side and allosteric_states untouched.
        for r in result
            _testhelper_assert_mechanism_invariants(r)
            @test r.allosteric.oligomeric_state == am.allosteric.oligomeric_state
            @test r.allosteric.allosteric_states == am.allosteric.allosteric_states
            @test _testhelper_keeps_reaction(r, am)
        end
    end
end

@testset "Allosteric mechanism — competitive-only name: no allo variant" begin
    # Build a uni-uni allosteric mechanism that already has :I bound as a
    # dead-end (added via init→dead-end→allosteric on the Mechanism path).
    # `rxn` declares :I only as a competitive inhibitor, not an allosteric
    # regulator, so `_expand_add_allosteric_regulator(am, rxn)` has no
    # allosteric regulator to add and returns empty.
    rxn = @enzyme_reaction begin
        substrates: S[C]
        products: P[C]
        competitive_inhibitors: I
        oligomeric_state: 2
    end
    am = @enzyme_mechanism begin
        substrates: S
        products: P
        modifiers: I
        oligomeric_state: 2
        steps: begin
            E + S ⇌ ES    :: EqualAI
            ES <--> EP  :: OnlyA
            EP ⇌ E + P    :: EqualAI
            E + I ⇌ EI    :: OnlyA
        end
    end
    _testhelper_assert_mechanism_invariants(am)
    @test isempty(
        _testhelper_move(EnzymeRates.Enumeration._expand_add_allosteric_regulator, am, rxn))
end

@testset "Allosteric mechanism — Two regulators with site options: count = 6" begin
    # SEED: allosteric uni-uni with R1 already added as :OnlyA
    # (existing site has one
    # non-:EqualAI ligand). R2 is un-added; the :EqualAI-at-existing
    # branch fires because R1 qualifies.
    am = @enzyme_mechanism begin
        substrates: S; products: P
        oligomeric_state: 2
        steps: begin
            E + P ⇌ EP    :: EqualAI
            E + S ⇌ ES    :: EqualAI
            ES <--> EP  :: EqualAI
        end
        regulatory_site(2): R1::OnlyA
    end
    _testhelper_assert_mechanism_invariants(am)

    result = _testhelper_move(EnzymeRates.Enumeration._expand_add_allosteric_regulator, am,
                              uni_uni_allo_2reg)

    # The exact children: R2 on a site of its own as :OnlyA, :OnlyI or :NonequalAI, and
    # on R1's :OnlyA site as :OnlyA, :NonequalAI or :EqualAI.
    expected = [
        (@enzyme_mechanism begin
            substrates: S; products: P
            oligomeric_state: 2
            steps: begin
                E + P ⇌ EP    :: EqualAI
                E + S ⇌ ES    :: EqualAI
                ES <--> EP  :: EqualAI
            end
            regulatory_site(2): R1::OnlyA
            regulatory_site(2): R2::OnlyA
        end),
        (@enzyme_mechanism begin
            substrates: S; products: P
            oligomeric_state: 2
            steps: begin
                E + P ⇌ EP    :: EqualAI
                E + S ⇌ ES    :: EqualAI
                ES <--> EP  :: EqualAI
            end
            regulatory_site(2): R1::OnlyA
            regulatory_site(2): R2::OnlyI
        end),
        (@enzyme_mechanism begin
            substrates: S; products: P
            oligomeric_state: 2
            steps: begin
                E + P ⇌ EP    :: EqualAI
                E + S ⇌ ES    :: EqualAI
                ES <--> EP  :: EqualAI
            end
            regulatory_site(2): R1::OnlyA
            regulatory_site(2): R2::NonequalAI
        end),
        (@enzyme_mechanism begin
            substrates: S; products: P
            oligomeric_state: 2
            steps: begin
                E + P ⇌ EP    :: EqualAI
                E + S ⇌ ES    :: EqualAI
                ES <--> EP  :: EqualAI
            end
            regulatory_site(2): R1::OnlyA, R2::OnlyA
        end),
        (@enzyme_mechanism begin
            substrates: S; products: P
            oligomeric_state: 2
            steps: begin
                E + P ⇌ EP    :: EqualAI
                E + S ⇌ ES    :: EqualAI
                ES <--> EP  :: EqualAI
            end
            regulatory_site(2): R1::OnlyA, R2::NonequalAI
        end),
        (@enzyme_mechanism begin
            substrates: S; products: P
            oligomeric_state: 2
            steps: begin
                E + P ⇌ EP    :: EqualAI
                E + S ⇌ ES    :: EqualAI
                ES <--> EP  :: EqualAI
            end
            regulatory_site(2): R1::OnlyA, R2::EqualAI
        end)]

    # 1. count: 3 non-:EqualAI tag flavors × 2 site options (new + R1's
    # existing) = 6, minus the redundant :OnlyI-onto-R1's-:OnlyA-site
    # append (disjoint conformations) = 5. Plus 1 :EqualAI-at-existing
    # variant (gated on R1 being non-:EqualAI). → 6.
    @test length(result) == 6
    @test Set(result) == Set(expected)

    # 2. Δ params: four variants add one parameter (:OnlyA new/existing,
    # :OnlyI new, and the :EqualAI-at-existing one), two :NonequalAI
    # variants add two (K_R2 + K_R2_T). Sorted [1, 1, 1, 1, 2, 2].
    # Measured against the actual compiled fitted count.
    @test _testhelper_param_deltas(am, result) == [1, 1, 1, 1, 2, 2]

    # 3. structural: every result has :R2 in some regulatory site.
    for r in result
        has_r2 = any(r.allosteric.regulatory_sites) do site
            any(l -> l === :R2,
                site.ligands)
        end
        @test has_r2
    end

    # 4. structural: at least one result has :R2 :EqualAI at site 1 (the same site as
    # R1). This is the :EqualAI-at-existing branch, gated on R1 (the existing ligand)
    # being non-:EqualAI.
    target = findfirst(result) do r
        length(r.allosteric.regulatory_sites) == 1 || return false
        site = r.allosteric.regulatory_sites[1]
        ligs = site.ligands
        states = site.allosteric_states
        idx = findfirst(l -> l === :R2, ligs)
        idx === nothing && return false
        states[idx] == :EqualAI
    end
    @test target !== nothing

    # 5. property-style: separate new-site (#sites grows to 2) and existing-site
    # (#sites stays at 1) placements.
    new_site_variants = filter(r -> length(r.allosteric.regulatory_sites) == 2, result)
    existing_site_variants = filter(r -> length(r.allosteric.regulatory_sites) == 1, result)
    @test length(new_site_variants) == 3   # 3 non-:EqualAI × new site
    # :OnlyA + :NonequalAI appends onto R1's :OnlyA site + 1 :EqualAI; the :OnlyI
    # append is skipped as redundant (disjoint conformations).
    @test length(existing_site_variants) == 3

    # 6. invariants and preservation: catalytic side and allosteric_states untouched.
    for r in result
        _testhelper_assert_mechanism_invariants(r)
        @test r.allosteric.oligomeric_state == am.allosteric.oligomeric_state
        @test r.allosteric.allosteric_states == am.allosteric.allosteric_states
        @test _testhelper_keeps_reaction(r, am)
    end
end

@testset "Allosteric mechanism — Adding :EqualAI R2 at site with :OnlyI R1" begin
    # SEED: R1::OnlyI at site 1. Adding R2 enumerates non-:EqualAI
    # tags × 2 site options
    # plus :EqualAI-at-existing (gated on R1 being non-:EqualAI). 7 total.
    am = @enzyme_mechanism begin
        substrates: S
        products: P
        oligomeric_state: 2
        steps: begin
            E + P ⇌ EP      :: EqualAI
            E + S ⇌ ES      :: EqualAI
            ES <--> EP    :: EqualAI
        end
        regulatory_site(2): R1::OnlyI
    end
    _testhelper_assert_mechanism_invariants(am)

    result = _testhelper_move(EnzymeRates.Enumeration._expand_add_allosteric_regulator, am,
                              uni_uni_allo_2reg)

    # 1. count: 3 non-:EqualAI × 2 sites + 1 :EqualAI-at-existing = 7,
    # minus the redundant :OnlyA-onto-R1's-:OnlyI-site append (disjoint
    # conformations) = 6.
    @test length(result) == 6

    # 2. Δ params: same multiset as the :OnlyA seed — four +1 variants
    # and two :NonequalAI +2 variants. Sorted [1, 1, 1, 1, 2, 2].
    # Measured against the actual compiled fitted count.
    @test _testhelper_param_deltas(am, result) == [1, 1, 1, 1, 2, 2]

    # 3. invariants; the deltas above compile every result.
    foreach(_testhelper_assert_mechanism_invariants, result)

    # 4. property: at least one variant has :R2 :EqualAI at site 1
    # (the :EqualAI-at-existing branch, gated on R1 being non-:EqualAI).
    has_eq_at_site1 = any(result) do r
        length(r.allosteric.regulatory_sites) == 1 || return false
        site = r.allosteric.regulatory_sites[1]
        ligs = site.ligands
        states = site.allosteric_states
        idx = findfirst(l -> l === :R2, ligs)
        idx === nothing && return false
        states[idx] == :EqualAI
    end
    @test has_eq_at_site1
end

@testset "Allosteric mechanism — all declared regs already present → empty" begin
    # SEED: allosteric uni-uni with :R already bound; rxn declares only :R → eligible_regs
    # is empty → result is empty.
    rxn = @enzyme_reaction begin
        substrates: S[C]
        products: P[C]
        allosteric_regulators: R
        oligomeric_state: 2
    end
    am = @enzyme_mechanism begin
        substrates: S
        products: P
        oligomeric_state: 2
        steps: begin
            E + P ⇌ EP      :: EqualAI
            E + S ⇌ ES      :: EqualAI
            ES <--> EP    :: EqualAI
        end
        regulatory_site(2): R::OnlyA
    end
    _testhelper_assert_mechanism_invariants(am)
    @test isempty(
        _testhelper_move(EnzymeRates.Enumeration._expand_add_allosteric_regulator, am, rxn))
end

@testset "Allosteric mechanism — dual-role name gains an allosteric site" begin
    # A reaction may declare one name (ATP) in BOTH roles. Starting from an
    # allosteric mechanism that already binds ATP as a competitive dead-end
    # inhibitor, `_expand_add_allosteric_regulator` also adds ATP at a
    # regulatory site — the two roles carry distinct parameter names.
    rxn = @enzyme_reaction begin
        substrates: S[C]
        products: P[C]
        competitive_inhibitors: ATP
        allosteric_regulators: ATP
    end
    am = @enzyme_mechanism begin
        substrates: S
        products: P
        modifiers: ATP
        oligomeric_state: 1
        steps: begin
            E + S ⇌ ES      :: EqualAI
            ES <--> EP    :: OnlyA
            EP ⇌ E + P      :: EqualAI
            E + ATP ⇌ EATP  :: OnlyA
        end
    end
    # am already binds ATP as a competitive dead-end inhibitor step.
    @test any(k -> _testhelper_binds_copy(_testhelper_candidate(am), k) &&
                   am.catalytic.metabolites[only(am.catalytic.steps[k].taken_up)] === :ATP,
              eachindex(am.catalytic.steps))

    result = _testhelper_move(EnzymeRates.Enumeration._expand_add_allosteric_regulator, am,
                              rxn)
    # ATP reaches a regulatory site despite its dead-end role: one site of
    # multiplicity 1 per tag the site can carry.
    expected = [(@enzyme_mechanism begin
                    substrates: S
                    products: P
                    modifiers: ATP
                    oligomeric_state: 1
                    steps: begin
                        E + S ⇌ ES      :: EqualAI
                        ES <--> EP    :: OnlyA
                        EP ⇌ E + P      :: EqualAI
                        E + ATP ⇌ EATP  :: OnlyA
                    end
                    regulatory_site(1): ATP::OnlyA
                end),
                (@enzyme_mechanism begin
                    substrates: S
                    products: P
                    modifiers: ATP
                    oligomeric_state: 1
                    steps: begin
                        E + S ⇌ ES      :: EqualAI
                        ES <--> EP    :: OnlyA
                        EP ⇌ E + P      :: EqualAI
                        E + ATP ⇌ EATP  :: OnlyA
                    end
                    regulatory_site(1): ATP::OnlyI
                end),
                @enzyme_mechanism begin
                    substrates: S
                    products: P
                    modifiers: ATP
                    oligomeric_state: 1
                    steps: begin
                        E + S ⇌ ES      :: EqualAI
                        ES <--> EP    :: OnlyA
                        EP ⇌ E + P      :: EqualAI
                        E + ATP ⇌ EATP  :: OnlyA
                    end
                    regulatory_site(1): ATP::NonequalAI
                end]
    @test length(result) == 3
    @test Set(result) == Set(expected)
    for r in result
        @test any(r.allosteric.regulatory_sites) do site
            any(l -> l === :ATP,
                site.ligands)
        end
        _testhelper_assert_mechanism_invariants(r)
    end
end

@testset "Allosteric mechanism — dual-role name derives distinct params" begin
    # ATP binds BOTH at a regulatory site AND as a competitive dead-end step, so it
    # takes two distinct constants (K_A_EATP_to_E_ATP at the catalytic site, K_A_ATPreg
    # and K_I_ATPreg at the regulatory one), never one name for both, and the mechanism
    # compiles.
    am = @enzyme_mechanism begin
        substrates: S
        products: P
        modifiers: ATP
        oligomeric_state: 1
        steps: begin
            E + S ⇌ ES            :: EqualAI
            ES <--> EP          :: OnlyA
            E + P ⇌ EP            :: EqualAI
            E + ATP ⇌ EATP        :: OnlyA
        end
        regulatory_site(1): ATP::NonequalAI
    end
    _testhelper_assert_mechanism_invariants(am)
    @test parameters(am) ==
          (:K_A_EATP_to_E_ATP, :K_EP_to_E_P, :K_ES_to_E_S, :k_A_ES_to_EP, :K_A_ATPreg,
           :K_I_ATPreg, :L)
    _testhelper_check_law_against_oracle(am)
end

@testset "Allosteric mechanism — dead-end move skips a regulatory-site ligand" begin
    # The reaction declares ATP both as a competitive inhibitor and as an allosteric
    # regulator, and the parent binds ATP at a regulatory site. The dead-end move leaves
    # a site ligand on its site, so its children bind I, at free E, its one placement,
    # with I's group :EqualAI, :OnlyA and :OnlyI.
    rxn = @enzyme_reaction begin
        substrates: S[C]
        products: P[C]
        competitive_inhibitors: ATP, I
        allosteric_regulators: ATP
        oligomeric_state: 2
    end
    parent = @enzyme_mechanism begin
        substrates: S; products: P
        oligomeric_state: 2
        steps: begin
            E + P ⇌ EP    :: EqualAI
            E + S ⇌ ES    :: EqualAI
            ES <--> EP  :: EqualAI
        end
        regulatory_site(2): ATP::OnlyI
    end
    children = _testhelper_move(EnzymeRates.Enumeration._expand_add_dead_end_regulator,
                                parent, rxn)
    expected = [
        (@enzyme_mechanism begin
            substrates: S; products: P
            modifiers: I
            oligomeric_state: 2
            steps: begin
                E + P ⇌ EP            :: EqualAI
                E + S ⇌ ES            :: EqualAI
                ES <--> EP          :: EqualAI
                E + I ⇌ EI  :: EqualAI
            end
            regulatory_site(2): ATP::OnlyI
        end),
        (@enzyme_mechanism begin
            substrates: S; products: P
            modifiers: I
            oligomeric_state: 2
            steps: begin
                E + P ⇌ EP            :: EqualAI
                E + S ⇌ ES            :: EqualAI
                ES <--> EP          :: EqualAI
                E + I ⇌ EI  :: OnlyA
            end
            regulatory_site(2): ATP::OnlyI
        end),
        (@enzyme_mechanism begin
            substrates: S; products: P
            modifiers: I
            oligomeric_state: 2
            steps: begin
                E + P ⇌ EP            :: EqualAI
                E + S ⇌ ES            :: EqualAI
                ES <--> EP          :: EqualAI
                E + I ⇌ EI  :: OnlyI
            end
            regulatory_site(2): ATP::OnlyI
        end)]
    @test length(children) == 3
    @test Set(children) == Set(expected)
    @test all(c -> :ATP ∉ _testhelper_modifier_names(c), children)
end

@testset "skips redundant OnlyI-onto-OnlyA-site append (disjoint states)" begin
    # Appending an :OnlyI regulator onto a site holding only an :OnlyA
    # regulator makes an [OnlyA, OnlyI] single site — disjoint conformations,
    # so it derives to the same rate as adding it at a new site. Skip that
    # append; keep the new-site form and the non-disjoint appends.
    rxn = @enzyme_reaction begin
        substrates: S[C]
        products: P[C]
        allosteric_regulators: A::Activator, I::Inhibitor
        oligomeric_state: 4
    end
    am = @enzyme_mechanism begin
        substrates: S; products: P
        oligomeric_state: 4
        steps: begin
            E + P ⇌ EP    :: OnlyA
            E + S ⇌ ES    :: OnlyA
            ES <--> EP  :: OnlyA
        end
        regulatory_site(4): A::OnlyA
    end
    kids = _testhelper_move(EnzymeRates.Enumeration._expand_add_allosteric_regulator, am,
                            rxn)
    # I on a site of its own at each pure or unequal tag, and on A's site as :OnlyA,
    # :NonequalAI or :EqualAI; not on A's site as :OnlyI.
    expected = [
        (@enzyme_mechanism begin
            substrates: S; products: P
            oligomeric_state: 4
            steps: begin
                E + P ⇌ EP    :: OnlyA
                E + S ⇌ ES    :: OnlyA
                ES <--> EP  :: OnlyA
            end
            regulatory_site(4): A::OnlyA
            regulatory_site(4): I::OnlyA
        end),
        (@enzyme_mechanism begin
            substrates: S; products: P
            oligomeric_state: 4
            steps: begin
                E + P ⇌ EP    :: OnlyA
                E + S ⇌ ES    :: OnlyA
                ES <--> EP  :: OnlyA
            end
            regulatory_site(4): A::OnlyA
            regulatory_site(4): I::OnlyI
        end),
        (@enzyme_mechanism begin
            substrates: S; products: P
            oligomeric_state: 4
            steps: begin
                E + P ⇌ EP    :: OnlyA
                E + S ⇌ ES    :: OnlyA
                ES <--> EP  :: OnlyA
            end
            regulatory_site(4): A::OnlyA
            regulatory_site(4): I::NonequalAI
        end),
        (@enzyme_mechanism begin
            substrates: S; products: P
            oligomeric_state: 4
            steps: begin
                E + P ⇌ EP    :: OnlyA
                E + S ⇌ ES    :: OnlyA
                ES <--> EP  :: OnlyA
            end
            regulatory_site(4): A::OnlyA, I::OnlyA
        end),
        (@enzyme_mechanism begin
            substrates: S; products: P
            oligomeric_state: 4
            steps: begin
                E + P ⇌ EP    :: OnlyA
                E + S ⇌ ES    :: OnlyA
                ES <--> EP  :: OnlyA
            end
            regulatory_site(4): A::OnlyA, I::NonequalAI
        end),
        @enzyme_mechanism begin
            substrates: S; products: P
            oligomeric_state: 4
            steps: begin
                E + P ⇌ EP    :: OnlyA
                E + S ⇌ ES    :: OnlyA
                ES <--> EP  :: OnlyA
            end
            regulatory_site(4): A::OnlyA, I::EqualAI
        end]
    @test length(kids) == 6
    @test Set(kids) == Set(expected)
    _testhelper_st(child, sidx, lig) = begin
        s = child.allosteric.regulatory_sites[sidx]
        i = findfirst(l -> l == lig, s.ligands)
        i === nothing ? nothing : s.allosteric_states[i]
    end
    _testhelper_one_AI_site(c) = length(c.allosteric.regulatory_sites) == 1 &&
        Set(only(c.allosteric.regulatory_sites).ligands) == Set([:A, :I])
    # redundant: one site with A:OnlyA + I:OnlyI — must be skipped
    @test !any(c -> _testhelper_one_AI_site(c) && _testhelper_st(c,1,:A)==:OnlyA &&
                    _testhelper_st(c,1,:I)==:OnlyI, kids)
    # kept: two-site new-site form carrying I:OnlyI on its own site
    @test any(kids) do c
        length(c.allosteric.regulatory_sites) == 2 &&
            any(s -> s.ligands == [:I] &&
                     s.allosteric_states == [:OnlyI],
                c.allosteric.regulatory_sites)
    end
    # kept: both-OnlyA append (A:OnlyA + I:OnlyA on one site — not disjoint)
    @test any(c -> _testhelper_one_AI_site(c) && _testhelper_st(c,1,:A)==:OnlyA &&
                   _testhelper_st(c,1,:I)==:OnlyA, kids)
end

@testset "adding Y at X's site and X at Y's site give one mechanism" begin
    # A site holding X::OnlyA and Y::NonequalAI is reached from a parent with X alone
    # (Y appended to X's site) and from a parent with Y alone (X appended to Y's
    # site). Both children are the same mechanism, so the beam dedups them.
    rxn = @enzyme_reaction begin
        substrates: S[C]
        products: P[C]
        allosteric_regulators: X, Y
        oligomeric_state: 2
    end
    am_x = @enzyme_mechanism begin
        substrates: S; products: P
        oligomeric_state: 2
        steps: begin
            E + P ⇌ EP    :: EqualAI
            E + S ⇌ ES    :: EqualAI
            ES <--> EP  :: EqualAI
        end
        regulatory_site(2): X::OnlyA
    end
    am_y = @enzyme_mechanism begin
        substrates: S; products: P
        oligomeric_state: 2
        steps: begin
            E + P ⇌ EP    :: EqualAI
            E + S ⇌ ES    :: EqualAI
            ES <--> EP  :: EqualAI
        end
        regulatory_site(2): Y::NonequalAI
    end
    # The site both routes reach: X::OnlyA and Y::NonequalAI share one site.
    shared = @enzyme_mechanism begin
        substrates: S; products: P
        oligomeric_state: 2
        steps: begin
            E + P ⇌ EP    :: EqualAI
            E + S ⇌ ES    :: EqualAI
            ES <--> EP  :: EqualAI
        end
        regulatory_site(2): X::OnlyA, Y::NonequalAI
    end
    # X::OnlyA and Y::NonequalAI on sites of their own: also reached from either parent.
    x_y_apart = @enzyme_mechanism begin
        substrates: S; products: P
        oligomeric_state: 2
        steps: begin
            E + P ⇌ EP    :: EqualAI
            E + S ⇌ ES    :: EqualAI
            ES <--> EP  :: EqualAI
        end
        regulatory_site(2): X::OnlyA
        regulatory_site(2): Y::NonequalAI
    end

    # From am_x, Y joins as :OnlyA or :NonequalAI on a new site or X's site, as :OnlyI on a
    # new site only (an :OnlyI ligand on X's all-:OnlyA site acts on a disjoint
    # conformation), and as :EqualAI on X's site.
    y_onlya_new = @enzyme_mechanism begin
        substrates: S; products: P
        oligomeric_state: 2
        steps: begin
            E + P ⇌ EP    :: EqualAI
            E + S ⇌ ES    :: EqualAI
            ES <--> EP  :: EqualAI
        end
        regulatory_site(2): X::OnlyA
        regulatory_site(2): Y::OnlyA
    end
    y_onlya_at_x = @enzyme_mechanism begin
        substrates: S; products: P
        oligomeric_state: 2
        steps: begin
            E + P ⇌ EP    :: EqualAI
            E + S ⇌ ES    :: EqualAI
            ES <--> EP  :: EqualAI
        end
        regulatory_site(2): X::OnlyA, Y::OnlyA
    end
    y_onlyi_new = @enzyme_mechanism begin
        substrates: S; products: P
        oligomeric_state: 2
        steps: begin
            E + P ⇌ EP    :: EqualAI
            E + S ⇌ ES    :: EqualAI
            ES <--> EP  :: EqualAI
        end
        regulatory_site(2): X::OnlyA
        regulatory_site(2): Y::OnlyI
    end
    y_equalai_at_x = @enzyme_mechanism begin
        substrates: S; products: P
        oligomeric_state: 2
        steps: begin
            E + P ⇌ EP    :: EqualAI
            E + S ⇌ ES    :: EqualAI
            ES <--> EP  :: EqualAI
        end
        regulatory_site(2): X::OnlyA, Y::EqualAI
    end
    from_x = _testhelper_move(EnzymeRates.Enumeration._expand_add_allosteric_regulator,
                              am_x, rxn)
    @test length(from_x) == 6
    @test Set(from_x) == Set([y_onlya_new, y_onlya_at_x, y_onlyi_new, x_y_apart, shared,
                              y_equalai_at_x])

    # From am_y, X joins with each non-:EqualAI tag on a new site or Y's site (Y's
    # :NonequalAI site acts on both conformations), and as :EqualAI on Y's site.
    x_onlyi_new = @enzyme_mechanism begin
        substrates: S; products: P
        oligomeric_state: 2
        steps: begin
            E + P ⇌ EP    :: EqualAI
            E + S ⇌ ES    :: EqualAI
            ES <--> EP  :: EqualAI
        end
        regulatory_site(2): X::OnlyI
        regulatory_site(2): Y::NonequalAI
    end
    x_onlyi_at_y = @enzyme_mechanism begin
        substrates: S; products: P
        oligomeric_state: 2
        steps: begin
            E + P ⇌ EP    :: EqualAI
            E + S ⇌ ES    :: EqualAI
            ES <--> EP  :: EqualAI
        end
        regulatory_site(2): X::OnlyI, Y::NonequalAI
    end
    x_nonequalai_new = @enzyme_mechanism begin
        substrates: S; products: P
        oligomeric_state: 2
        steps: begin
            E + P ⇌ EP    :: EqualAI
            E + S ⇌ ES    :: EqualAI
            ES <--> EP  :: EqualAI
        end
        regulatory_site(2): X::NonequalAI
        regulatory_site(2): Y::NonequalAI
    end
    x_nonequalai_at_y = @enzyme_mechanism begin
        substrates: S; products: P
        oligomeric_state: 2
        steps: begin
            E + P ⇌ EP    :: EqualAI
            E + S ⇌ ES    :: EqualAI
            ES <--> EP  :: EqualAI
        end
        regulatory_site(2): X::NonequalAI, Y::NonequalAI
    end
    x_equalai_at_y = @enzyme_mechanism begin
        substrates: S; products: P
        oligomeric_state: 2
        steps: begin
            E + P ⇌ EP    :: EqualAI
            E + S ⇌ ES    :: EqualAI
            ES <--> EP  :: EqualAI
        end
        regulatory_site(2): X::EqualAI, Y::NonequalAI
    end
    from_y = _testhelper_move(EnzymeRates.Enumeration._expand_add_allosteric_regulator,
                              am_y, rxn)
    @test length(from_y) == 7
    @test Set(from_y) == Set([x_y_apart, shared, x_onlyi_new, x_onlyi_at_y,
                              x_nonequalai_new, x_nonequalai_at_y, x_equalai_at_y])

    # The shared site comes out of each route once, as one mechanism under == and hash.
    @test count(==(shared), from_x) == 1
    @test count(==(shared), from_y) == 1
    via_x = EnzymeRates.Enumeration._make_am_with_added_reg(_testhelper_candidate(am_x), :Y,
                                                            :NonequalAI, 1).mechanism
    via_y = EnzymeRates.Enumeration._make_am_with_added_reg(_testhelper_candidate(am_y), :X,
                                                            :OnlyA, 1).mechanism
    @test via_x == via_y == shared
    @test hash(via_x) == hash(via_y)
end

end

# ─── regulator-kind routing ────────────────────────────────────────────
@testset "regulator kinds route to their own moves" begin
# uni_uni_reg_and_inhibitor declares R2 as an allosteric regulator and R1
# as a competitive inhibitor. The two kinds must not cross moves: the
# allosteric-regulator move places only R2 (at a regulatory site), and the
# dead-end move binds only R1.
@testset "allosteric-regulator addition places only R2" begin
    am = @enzyme_mechanism begin
        substrates: S
        products: P
        oligomeric_state: 2
        steps: begin
            E + S ⇌ ES    :: EqualAI
            ES <--> EP  :: OnlyA
            EP ⇌ E + P    :: OnlyA
        end
    end
    added = _testhelper_move(EnzymeRates.Enumeration._expand_add_allosteric_regulator, am,
                             uni_uni_reg_and_inhibitor)
    # R2 at a new site, once per tag the site can carry.
    expected = [
        (@enzyme_mechanism begin
            substrates: S
            products: P
            oligomeric_state: 2
            steps: begin
                E + S ⇌ ES    :: EqualAI
                ES <--> EP  :: OnlyA
                EP ⇌ E + P    :: OnlyA
            end
            regulatory_site(2): R2::OnlyA
        end),
        (@enzyme_mechanism begin
            substrates: S
            products: P
            oligomeric_state: 2
            steps: begin
                E + S ⇌ ES    :: EqualAI
                ES <--> EP  :: OnlyA
                EP ⇌ E + P    :: OnlyA
            end
            regulatory_site(2): R2::OnlyI
        end),
        (@enzyme_mechanism begin
            substrates: S
            products: P
            oligomeric_state: 2
            steps: begin
                E + S ⇌ ES    :: EqualAI
                ES <--> EP  :: OnlyA
                EP ⇌ E + P    :: OnlyA
            end
            regulatory_site(2): R2::NonequalAI
        end),
    ]
    @test length(added) == 3
    @test Set(added) == Set(expected)
    site_regs = Set(l for m in added
                    for s in m.allosteric.regulatory_sites
                    for l in s.ligands)
    @test site_regs == Set([:R2])
end

@testset "dead-end expansion binds only R1" begin
    seed = @enzyme_mechanism begin
        substrates: S
        products: P
        steps: begin
            E + S ⇌ ES
            ES <--> EP
            EP ⇌ E + P
        end
    end
    added = _testhelper_move(EnzymeRates.Enumeration._expand_add_dead_end_regulator, seed,
                             uni_uni_reg_and_inhibitor)
    # The one child binds R1 at free E.
    @test added == [@enzyme_mechanism begin
        substrates: S
        products: P
        modifiers: R1
        steps: begin
            E + S ⇌ ES
            ES <--> EP
            EP ⇌ E + P
            E + R1 ⇌ ER1
        end
    end]
    de_regs = Set(x for m in added for x in _testhelper_modifier_names(m))
    @test de_regs == Set([:R1])
end
end

# ─── _expand_change_allo_state ─────────────────────────────────────────
@testset "_expand_change_allo_state" begin

@testset "change_allo_state tags an :OnlyA half-reaction :NonequalAI alone" begin
    # The six-step ping-pong with both half-reactions at steady state, at one subunit, its
    # first half-reaction :OnlyA. The inactive conformation then turns nothing over, so
    # each group tagged :NonequalAI alone gains: each binding its one constant, the second
    # half-reaction both of its rate constants, since no cycle of the inactive
    # conformation ties either, and the stopped half-reaction one, since running in the
    # inactive conformation it closes the cycle whose Keq ties the other. Six children.
    parent = @enzyme_mechanism begin
        substrates: A, B ; products: P, Q ; oligomeric_state: 1
        steps: begin
            E + A ⇌ EA                                       :: EqualAI
            E_res_plusA_minusP + B ⇌ EB_res_plusA_minusP :: EqualAI
            E + Q ⇌ EQ                                       :: EqualAI
            E_res_plusA_minusP + P ⇌ EP_res_plusA_minusP :: EqualAI
            EA <--> EP_res_plusA_minusP                   :: OnlyA
            EB_res_plusA_minusP <--> EQ                   :: EqualAI
        end
    end
    first_half = @enzyme_mechanism begin
        substrates: A, B ; products: P, Q ; oligomeric_state: 1
        steps: begin
            E + A ⇌ EA                                       :: EqualAI
            E_res_plusA_minusP + B ⇌ EB_res_plusA_minusP :: EqualAI
            E + Q ⇌ EQ                                       :: EqualAI
            E_res_plusA_minusP + P ⇌ EP_res_plusA_minusP :: EqualAI
            EA <--> EP_res_plusA_minusP                   :: NonequalAI
            EB_res_plusA_minusP <--> EQ                   :: EqualAI
        end
    end
    second_half = @enzyme_mechanism begin
        substrates: A, B ; products: P, Q ; oligomeric_state: 1
        steps: begin
            E + A ⇌ EA                                       :: EqualAI
            E_res_plusA_minusP + B ⇌ EB_res_plusA_minusP :: EqualAI
            E + Q ⇌ EQ                                       :: EqualAI
            E_res_plusA_minusP + P ⇌ EP_res_plusA_minusP :: EqualAI
            EA <--> EP_res_plusA_minusP                   :: OnlyA
            EB_res_plusA_minusP <--> EQ                   :: NonequalAI
        end
    end
    a = @enzyme_mechanism begin
        substrates: A, B ; products: P, Q ; oligomeric_state: 1
        steps: begin
            E + A ⇌ EA                                       :: NonequalAI
            E_res_plusA_minusP + B ⇌ EB_res_plusA_minusP :: EqualAI
            E + Q ⇌ EQ                                       :: EqualAI
            E_res_plusA_minusP + P ⇌ EP_res_plusA_minusP :: EqualAI
            EA <--> EP_res_plusA_minusP                   :: OnlyA
            EB_res_plusA_minusP <--> EQ                   :: EqualAI
        end
    end
    b = @enzyme_mechanism begin
        substrates: A, B ; products: P, Q ; oligomeric_state: 1
        steps: begin
            E + A ⇌ EA                                       :: EqualAI
            E_res_plusA_minusP + B ⇌ EB_res_plusA_minusP :: NonequalAI
            E + Q ⇌ EQ                                       :: EqualAI
            E_res_plusA_minusP + P ⇌ EP_res_plusA_minusP :: EqualAI
            EA <--> EP_res_plusA_minusP                   :: OnlyA
            EB_res_plusA_minusP <--> EQ                   :: EqualAI
        end
    end
    q = @enzyme_mechanism begin
        substrates: A, B ; products: P, Q ; oligomeric_state: 1
        steps: begin
            E + A ⇌ EA                                       :: EqualAI
            E_res_plusA_minusP + B ⇌ EB_res_plusA_minusP :: EqualAI
            E + Q ⇌ EQ                                       :: NonequalAI
            E_res_plusA_minusP + P ⇌ EP_res_plusA_minusP :: EqualAI
            EA <--> EP_res_plusA_minusP                   :: OnlyA
            EB_res_plusA_minusP <--> EQ                   :: EqualAI
        end
    end
    p = @enzyme_mechanism begin
        substrates: A, B ; products: P, Q ; oligomeric_state: 1
        steps: begin
            E + A ⇌ EA                                       :: EqualAI
            E_res_plusA_minusP + B ⇌ EB_res_plusA_minusP :: EqualAI
            E + Q ⇌ EQ                                       :: EqualAI
            E_res_plusA_minusP + P ⇌ EP_res_plusA_minusP :: NonequalAI
            EA <--> EP_res_plusA_minusP                   :: OnlyA
            EB_res_plusA_minusP <--> EQ                   :: EqualAI
        end
    end
    kids = _testhelper_move(EnzymeRates.Enumeration._expand_change_allo_state, parent)
    @test length(kids) == 6
    @test Set(kids) == Set([first_half, second_half, a, b, q, p])
    # A steady-state group adds two constants where a cycle of the inactive conformation
    # does not tie one of them, a binding at rapid equilibrium one.
    @test _testhelper_param_deltas(parent, kids) == [1, 1, 1, 1, 1, 2]
end

@testset "change_allo_state tags a pair of bindings together" begin
    # The six-step ping-pong with both half-reactions at steady state, at one subunit,
    # A's binding and P's :OnlyA. A tagged :NonequalAI beside P :OnlyA, or P beside A, is
    # refused: the inactive conformation then joins the two sides of the other binding,
    # whose affinity has no path to zero (`_FAST_STEP_PINNED`). Tagged together they
    # gain. The first half-reaction tagged :NonequalAI adds nothing, since the inactive
    # conformation, binding neither A nor P, never reaches its forms; it is not tagged.
    # Q's binding, B's and the second half-reaction gain alone. Four children.
    parent = @enzyme_mechanism begin
        substrates: A, B ; products: P, Q ; oligomeric_state: 1
        steps: begin
            E + A ⇌ EA                                       :: OnlyA
            E_res_plusA_minusP + B ⇌ EB_res_plusA_minusP :: EqualAI
            E + Q ⇌ EQ                                       :: EqualAI
            E_res_plusA_minusP + P ⇌ EP_res_plusA_minusP :: OnlyA
            EA <--> EP_res_plusA_minusP                   :: EqualAI
            EB_res_plusA_minusP <--> EQ                   :: EqualAI
        end
    end
    a_and_p = @enzyme_mechanism begin
        substrates: A, B ; products: P, Q ; oligomeric_state: 1
        steps: begin
            E + A ⇌ EA                                       :: NonequalAI
            E_res_plusA_minusP + B ⇌ EB_res_plusA_minusP :: EqualAI
            E + Q ⇌ EQ                                       :: EqualAI
            E_res_plusA_minusP + P ⇌ EP_res_plusA_minusP :: NonequalAI
            EA <--> EP_res_plusA_minusP                   :: EqualAI
            EB_res_plusA_minusP <--> EQ                   :: EqualAI
        end
    end
    b = @enzyme_mechanism begin
        substrates: A, B ; products: P, Q ; oligomeric_state: 1
        steps: begin
            E + A ⇌ EA                                       :: OnlyA
            E_res_plusA_minusP + B ⇌ EB_res_plusA_minusP :: NonequalAI
            E + Q ⇌ EQ                                       :: EqualAI
            E_res_plusA_minusP + P ⇌ EP_res_plusA_minusP :: OnlyA
            EA <--> EP_res_plusA_minusP                   :: EqualAI
            EB_res_plusA_minusP <--> EQ                   :: EqualAI
        end
    end
    q = @enzyme_mechanism begin
        substrates: A, B ; products: P, Q ; oligomeric_state: 1
        steps: begin
            E + A ⇌ EA                                       :: OnlyA
            E_res_plusA_minusP + B ⇌ EB_res_plusA_minusP :: EqualAI
            E + Q ⇌ EQ                                       :: NonequalAI
            E_res_plusA_minusP + P ⇌ EP_res_plusA_minusP :: OnlyA
            EA <--> EP_res_plusA_minusP                   :: EqualAI
            EB_res_plusA_minusP <--> EQ                   :: EqualAI
        end
    end
    second_half = @enzyme_mechanism begin
        substrates: A, B ; products: P, Q ; oligomeric_state: 1
        steps: begin
            E + A ⇌ EA                                       :: OnlyA
            E_res_plusA_minusP + B ⇌ EB_res_plusA_minusP :: EqualAI
            E + Q ⇌ EQ                                       :: EqualAI
            E_res_plusA_minusP + P ⇌ EP_res_plusA_minusP :: OnlyA
            EA <--> EP_res_plusA_minusP                   :: EqualAI
            EB_res_plusA_minusP <--> EQ                   :: NonequalAI
        end
    end
    a_alone = @_testhelper_unguarded @enzyme_mechanism begin
        substrates: A, B ; products: P, Q ; oligomeric_state: 1
        steps: begin
            E + A ⇌ EA                                       :: NonequalAI
            E_res_plusA_minusP + B ⇌ EB_res_plusA_minusP :: EqualAI
            E + Q ⇌ EQ                                       :: EqualAI
            E_res_plusA_minusP + P ⇌ EP_res_plusA_minusP :: OnlyA
            EA <--> EP_res_plusA_minusP                   :: EqualAI
            EB_res_plusA_minusP <--> EQ                   :: EqualAI
        end
    end
    first_half = @enzyme_mechanism begin
        substrates: A, B ; products: P, Q ; oligomeric_state: 1
        steps: begin
            E + A ⇌ EA                                       :: OnlyA
            E_res_plusA_minusP + B ⇌ EB_res_plusA_minusP :: EqualAI
            E + Q ⇌ EQ                                       :: EqualAI
            E_res_plusA_minusP + P ⇌ EP_res_plusA_minusP :: OnlyA
            EA <--> EP_res_plusA_minusP                   :: NonequalAI
            EB_res_plusA_minusP <--> EQ                   :: EqualAI
        end
    end
    kids = _testhelper_move(EnzymeRates.Enumeration._expand_change_allo_state, parent)
    @test length(kids) == 4
    @test Set(kids) == Set([a_and_p, b, q, second_half])
    @test EnzymeRates.Derivation._guard(a_alone).fault ===
          EnzymeRates.Derivation._FAST_STEP_PINNED
    @test !EnzymeRates.Enumeration._enumerator_returns(_testhelper_candidate(a_alone))
    @test _testhelper_fitted(first_half) == _testhelper_fitted(parent)
end

@testset "change_allo_state tags no binding of an all-:EqualAI ping-pong alone" begin
    # The six-step ping-pong with both half-reactions at steady state, at one subunit,
    # every group :EqualAI. The inactive conformation runs the whole cycle, so it must give
    # the same Keq: one binding tagged :NonequalAI gets K_I = K_A from thermodynamics and
    # adds nothing, and two bindings gain one constant, the second tied to the first.
    # Either half-reaction alone gains, its inactive Haldane relation fixing one of its
    # two rate constants. The eight minimal gaining sets: each half-reaction, and the six
    # pairs of bindings.
    parent = @enzyme_mechanism begin
        substrates: A, B ; products: P, Q ; oligomeric_state: 1
        steps: begin
            E + A ⇌ EA                                       :: EqualAI
            E_res_plusA_minusP + B ⇌ EB_res_plusA_minusP :: EqualAI
            E + Q ⇌ EQ                                       :: EqualAI
            E_res_plusA_minusP + P ⇌ EP_res_plusA_minusP :: EqualAI
            EA <--> EP_res_plusA_minusP                   :: EqualAI
            EB_res_plusA_minusP <--> EQ                   :: EqualAI
        end
    end
    first_half = @enzyme_mechanism begin
        substrates: A, B ; products: P, Q ; oligomeric_state: 1
        steps: begin
            E + A ⇌ EA                                       :: EqualAI
            E_res_plusA_minusP + B ⇌ EB_res_plusA_minusP :: EqualAI
            E + Q ⇌ EQ                                       :: EqualAI
            E_res_plusA_minusP + P ⇌ EP_res_plusA_minusP :: EqualAI
            EA <--> EP_res_plusA_minusP                   :: NonequalAI
            EB_res_plusA_minusP <--> EQ                   :: EqualAI
        end
    end
    second_half = @enzyme_mechanism begin
        substrates: A, B ; products: P, Q ; oligomeric_state: 1
        steps: begin
            E + A ⇌ EA                                       :: EqualAI
            E_res_plusA_minusP + B ⇌ EB_res_plusA_minusP :: EqualAI
            E + Q ⇌ EQ                                       :: EqualAI
            E_res_plusA_minusP + P ⇌ EP_res_plusA_minusP :: EqualAI
            EA <--> EP_res_plusA_minusP                   :: EqualAI
            EB_res_plusA_minusP <--> EQ                   :: NonequalAI
        end
    end
    a_and_b = @enzyme_mechanism begin
        substrates: A, B ; products: P, Q ; oligomeric_state: 1
        steps: begin
            E + A ⇌ EA                                       :: NonequalAI
            E_res_plusA_minusP + B ⇌ EB_res_plusA_minusP :: NonequalAI
            E + Q ⇌ EQ                                       :: EqualAI
            E_res_plusA_minusP + P ⇌ EP_res_plusA_minusP :: EqualAI
            EA <--> EP_res_plusA_minusP                   :: EqualAI
            EB_res_plusA_minusP <--> EQ                   :: EqualAI
        end
    end
    a_and_q = @enzyme_mechanism begin
        substrates: A, B ; products: P, Q ; oligomeric_state: 1
        steps: begin
            E + A ⇌ EA                                       :: NonequalAI
            E_res_plusA_minusP + B ⇌ EB_res_plusA_minusP :: EqualAI
            E + Q ⇌ EQ                                       :: NonequalAI
            E_res_plusA_minusP + P ⇌ EP_res_plusA_minusP :: EqualAI
            EA <--> EP_res_plusA_minusP                   :: EqualAI
            EB_res_plusA_minusP <--> EQ                   :: EqualAI
        end
    end
    a_and_p = @enzyme_mechanism begin
        substrates: A, B ; products: P, Q ; oligomeric_state: 1
        steps: begin
            E + A ⇌ EA                                       :: NonequalAI
            E_res_plusA_minusP + B ⇌ EB_res_plusA_minusP :: EqualAI
            E + Q ⇌ EQ                                       :: EqualAI
            E_res_plusA_minusP + P ⇌ EP_res_plusA_minusP :: NonequalAI
            EA <--> EP_res_plusA_minusP                   :: EqualAI
            EB_res_plusA_minusP <--> EQ                   :: EqualAI
        end
    end
    b_and_q = @enzyme_mechanism begin
        substrates: A, B ; products: P, Q ; oligomeric_state: 1
        steps: begin
            E + A ⇌ EA                                       :: EqualAI
            E_res_plusA_minusP + B ⇌ EB_res_plusA_minusP :: NonequalAI
            E + Q ⇌ EQ                                       :: NonequalAI
            E_res_plusA_minusP + P ⇌ EP_res_plusA_minusP :: EqualAI
            EA <--> EP_res_plusA_minusP                   :: EqualAI
            EB_res_plusA_minusP <--> EQ                   :: EqualAI
        end
    end
    b_and_p = @enzyme_mechanism begin
        substrates: A, B ; products: P, Q ; oligomeric_state: 1
        steps: begin
            E + A ⇌ EA                                       :: EqualAI
            E_res_plusA_minusP + B ⇌ EB_res_plusA_minusP :: NonequalAI
            E + Q ⇌ EQ                                       :: EqualAI
            E_res_plusA_minusP + P ⇌ EP_res_plusA_minusP :: NonequalAI
            EA <--> EP_res_plusA_minusP                   :: EqualAI
            EB_res_plusA_minusP <--> EQ                   :: EqualAI
        end
    end
    q_and_p = @enzyme_mechanism begin
        substrates: A, B ; products: P, Q ; oligomeric_state: 1
        steps: begin
            E + A ⇌ EA                                       :: EqualAI
            E_res_plusA_minusP + B ⇌ EB_res_plusA_minusP :: EqualAI
            E + Q ⇌ EQ                                       :: NonequalAI
            E_res_plusA_minusP + P ⇌ EP_res_plusA_minusP :: NonequalAI
            EA <--> EP_res_plusA_minusP                   :: EqualAI
            EB_res_plusA_minusP <--> EQ                   :: EqualAI
        end
    end
    a_alone = @enzyme_mechanism begin
        substrates: A, B ; products: P, Q ; oligomeric_state: 1
        steps: begin
            E + A ⇌ EA                                       :: NonequalAI
            E_res_plusA_minusP + B ⇌ EB_res_plusA_minusP :: EqualAI
            E + Q ⇌ EQ                                       :: EqualAI
            E_res_plusA_minusP + P ⇌ EP_res_plusA_minusP :: EqualAI
            EA <--> EP_res_plusA_minusP                   :: EqualAI
            EB_res_plusA_minusP <--> EQ                   :: EqualAI
        end
    end
    kids = _testhelper_move(EnzymeRates.Enumeration._expand_change_allo_state, parent)
    @test length(kids) == 8
    @test Set(kids) == Set([first_half, second_half, a_and_b, a_and_q, a_and_p, b_and_q,
                            b_and_p, q_and_p])
    @test _testhelper_param_deltas(parent, kids) == fill(1, 8)
    @test _testhelper_fitted(a_alone) == _testhelper_fitted(parent)
end

@testset "change_allo_state tags each group alone once two bindings are :NonequalAI" begin
    # The six-step ping-pong with both half-reactions at steady state, at one subunit, A's
    # binding and B's :NonequalAI and every other group :EqualAI: a child of the pair A
    # with B, the second round of the move. What a group gains never falls as other groups
    # become :NonequalAI: its constants are new unknowns, the thermodynamic equations that
    # bind them are fixed, and they bind less as more constants are free. The pair has
    # already taken the one constant the inactive cycle ties, so Q's binding and P's each
    # gain one alone, and each half-reaction gains both of its rate constants, where from
    # the all-:EqualAI parent it gains one. Four children.
    parent = @enzyme_mechanism begin
        substrates: A, B ; products: P, Q ; oligomeric_state: 1
        steps: begin
            E + A ⇌ EA                                       :: NonequalAI
            E_res_plusA_minusP + B ⇌ EB_res_plusA_minusP :: NonequalAI
            E + Q ⇌ EQ                                       :: EqualAI
            E_res_plusA_minusP + P ⇌ EP_res_plusA_minusP :: EqualAI
            EA <--> EP_res_plusA_minusP                   :: EqualAI
            EB_res_plusA_minusP <--> EQ                   :: EqualAI
        end
    end
    q = @enzyme_mechanism begin
        substrates: A, B ; products: P, Q ; oligomeric_state: 1
        steps: begin
            E + A ⇌ EA                                       :: NonequalAI
            E_res_plusA_minusP + B ⇌ EB_res_plusA_minusP :: NonequalAI
            E + Q ⇌ EQ                                       :: NonequalAI
            E_res_plusA_minusP + P ⇌ EP_res_plusA_minusP :: EqualAI
            EA <--> EP_res_plusA_minusP                   :: EqualAI
            EB_res_plusA_minusP <--> EQ                   :: EqualAI
        end
    end
    p = @enzyme_mechanism begin
        substrates: A, B ; products: P, Q ; oligomeric_state: 1
        steps: begin
            E + A ⇌ EA                                       :: NonequalAI
            E_res_plusA_minusP + B ⇌ EB_res_plusA_minusP :: NonequalAI
            E + Q ⇌ EQ                                       :: EqualAI
            E_res_plusA_minusP + P ⇌ EP_res_plusA_minusP :: NonequalAI
            EA <--> EP_res_plusA_minusP                   :: EqualAI
            EB_res_plusA_minusP <--> EQ                   :: EqualAI
        end
    end
    first_half = @enzyme_mechanism begin
        substrates: A, B ; products: P, Q ; oligomeric_state: 1
        steps: begin
            E + A ⇌ EA                                       :: NonequalAI
            E_res_plusA_minusP + B ⇌ EB_res_plusA_minusP :: NonequalAI
            E + Q ⇌ EQ                                       :: EqualAI
            E_res_plusA_minusP + P ⇌ EP_res_plusA_minusP :: EqualAI
            EA <--> EP_res_plusA_minusP                   :: NonequalAI
            EB_res_plusA_minusP <--> EQ                   :: EqualAI
        end
    end
    second_half = @enzyme_mechanism begin
        substrates: A, B ; products: P, Q ; oligomeric_state: 1
        steps: begin
            E + A ⇌ EA                                       :: NonequalAI
            E_res_plusA_minusP + B ⇌ EB_res_plusA_minusP :: NonequalAI
            E + Q ⇌ EQ                                       :: EqualAI
            E_res_plusA_minusP + P ⇌ EP_res_plusA_minusP :: EqualAI
            EA <--> EP_res_plusA_minusP                   :: EqualAI
            EB_res_plusA_minusP <--> EQ                   :: NonequalAI
        end
    end
    kids = _testhelper_move(EnzymeRates.Enumeration._expand_change_allo_state, parent)
    @test length(kids) == 4
    @test Set(kids) == Set([q, p, first_half, second_half])
    @test _testhelper_param_deltas(parent, kids) == [1, 1, 2, 2]
end

@testset "a ping-pong whose fast second half is :NonequalAI with the first is refused" begin
    # The six-step ping-pong with the first half-reaction at steady state and the second
    # at rapid equilibrium, at one subunit, R at an :OnlyA site, both half-reactions
    # :NonequalAI and every binding :EqualAI. The conformations then run different
    # chemistry on one rapid-equilibrium segment, and with every substrate zero, or B
    # zero, the covalent forms hold the enzyme in both, joined only through free E: a
    # split whose limit depends on the speeds (`_SPLIT`). So the enumerator never returns
    # it, whichever move builds it.
    both = @_testhelper_unguarded @enzyme_mechanism begin
        substrates: A, B ; products: P, Q ; oligomeric_state: 1
        steps: begin
            E + A ⇌ EA                                       :: EqualAI
            E_res_plusA_minusP + B ⇌ EB_res_plusA_minusP :: EqualAI
            E + Q ⇌ EQ                                       :: EqualAI
            E_res_plusA_minusP + P ⇌ EP_res_plusA_minusP :: EqualAI
            EA <--> EP_res_plusA_minusP                   :: NonequalAI
            EB_res_plusA_minusP ⇌ EQ                      :: NonequalAI
        end
        regulatory_site(1): R::OnlyA
    end
    D = EnzymeRates.Derivation
    v = D._guard(both)
    names = D._face_names(both.catalytic, D._faces(both.catalytic)[v.face])
    @test (v.fault, v.limit, names) == (D._SPLIT, D._LIMIT_SPEED, [:A, :B])
    @test !EnzymeRates.Enumeration._enumerator_returns(_testhelper_candidate(both))
end

@testset "change_allo_state on a parent with every group :NonequalAI" begin
    # No group is left to tag, and no regulatory site: no child.
    parent = @enzyme_mechanism begin
        substrates: A, B ; products: P, Q ; oligomeric_state: 1
        steps: begin
            E + A ⇌ EA                                       :: NonequalAI
            E_res_plusA_minusP + B ⇌ EB_res_plusA_minusP :: NonequalAI
            E + Q ⇌ EQ                                       :: NonequalAI
            E_res_plusA_minusP + P ⇌ EP_res_plusA_minusP :: NonequalAI
            EA <--> EP_res_plusA_minusP                   :: NonequalAI
            EB_res_plusA_minusP <--> EQ                   :: NonequalAI
        end
    end
    kids = _testhelper_move(EnzymeRates.Enumeration._expand_change_allo_state, parent)
    @test length(kids) == 0
    @test Set(kids) == Set()
end

@testset "change_allo_state tags an :OnlyI inhibitor :NonequalAI" begin
    # The uni-uni mechanism with the competitive inhibitor I bound at free E, I's group
    # :OnlyI and every other group :EqualAI. Three minimal sets gain a constant: I's
    # group, which then binds the active conformation too with a constant of its own;
    # the chemistry; and P's binding with S's, since one binding alone gets back the
    # active conformation's constant from thermodynamics.
    parent = @enzyme_mechanism begin
        substrates: S ; products: P ; modifiers: I ; oligomeric_state: 1
        steps: begin
            E + I ⇌ EI :: OnlyI
            E + P ⇌ EP :: EqualAI
            E + S ⇌ ES :: EqualAI
            ES <--> EP :: EqualAI
        end
    end
    i_nonequal = @enzyme_mechanism begin
        substrates: S ; products: P ; modifiers: I ; oligomeric_state: 1
        steps: begin
            E + I ⇌ EI :: NonequalAI
            E + P ⇌ EP :: EqualAI
            E + S ⇌ ES :: EqualAI
            ES <--> EP :: EqualAI
        end
    end
    chem_nonequal = @enzyme_mechanism begin
        substrates: S ; products: P ; modifiers: I ; oligomeric_state: 1
        steps: begin
            E + I ⇌ EI :: OnlyI
            E + P ⇌ EP :: EqualAI
            E + S ⇌ ES :: EqualAI
            ES <--> EP :: NonequalAI
        end
    end
    p_s_nonequal = @enzyme_mechanism begin
        substrates: S ; products: P ; modifiers: I ; oligomeric_state: 1
        steps: begin
            E + I ⇌ EI :: OnlyI
            E + P ⇌ EP :: NonequalAI
            E + S ⇌ ES :: NonequalAI
            ES <--> EP :: EqualAI
        end
    end
    p_nonequal = @enzyme_mechanism begin
        substrates: S ; products: P ; modifiers: I ; oligomeric_state: 1
        steps: begin
            E + I ⇌ EI :: OnlyI
            E + P ⇌ EP :: NonequalAI
            E + S ⇌ ES :: EqualAI
            ES <--> EP :: EqualAI
        end
    end
    kids = _testhelper_move(EnzymeRates.Enumeration._expand_change_allo_state, parent)
    @test length(kids) == 3
    @test Set(kids) == Set([i_nonequal, chem_nonequal, p_s_nonequal])
    fitted = EnzymeRates.Derivation.fitted_parameter_count
    @test fitted(i_nonequal) == fitted(chem_nonequal) == fitted(p_s_nonequal) ==
          fitted(parent) + 1
    @test fitted(p_nonequal) == fitted(parent)
end

@testset "change_allo_state keeps a regulatory ligand's child that adds a constant" begin
    # The uni-uni mechanism at two subunits with every kinetic group :NonequalAI and R
    # :OnlyA at a site. R tagged :NonequalAI binds the inactive conformation with a
    # constant of its own, one more fitted constant: a site's constants enter no
    # thermodynamic cycle. The only child.
    parent = @enzyme_mechanism begin
        substrates: S ; products: P ; oligomeric_state: 2
        steps: begin
            E + S ⇌ ES      :: NonequalAI
            ES <--> EP    :: NonequalAI
            E + P ⇌ EP      :: NonequalAI
        end
        regulatory_site(2): R::OnlyA
    end
    r = @enzyme_mechanism begin
        substrates: S ; products: P ; oligomeric_state: 2
        steps: begin
            E + S ⇌ ES      :: NonequalAI
            ES <--> EP    :: NonequalAI
            E + P ⇌ EP      :: NonequalAI
        end
        regulatory_site(2): R::NonequalAI
    end
    kids = _testhelper_move(EnzymeRates.Enumeration._expand_change_allo_state, parent)
    @test length(kids) == 1
    @test Set(kids) == Set([r])
    @test _testhelper_param_deltas(parent, kids) == [1]
end

@testset "change_allo_state returns in time on a ten-group parent" begin
    # An ordered ter-ter with three abortive complexes, ten kinetic groups, at one
    # subunit, R's binding at free E and the chemistry :OnlyA. The inactive conformation
    # binds A, B and C and R at EAB, and never reaches EPQR, EQR, ER or EAR: A, B, C and R
    # at either form each gain alone, and the chemistry, P, Q and A at EQR or at ER add
    # nothing or are refused in every set the search builds from them. The search tries
    # 36 sets. Five children.
    parent = @enzyme_mechanism begin
        substrates: A, B, C ; products: P, Q, R ; oligomeric_state: 1
        steps: begin
            E + A ⇌ EA            :: EqualAI
            EA + B ⇌ EAB          :: EqualAI
            EAB + C ⇌ EABC        :: EqualAI
            EABC <--> EPQR        :: OnlyA
            EQR + P ⇌ EPQR        :: EqualAI
            ER + Q ⇌ EQR          :: EqualAI
            E + R ⇌ ER            :: OnlyA
            ER + A ⇌ EAR          :: EqualAI
            EQR + A ⇌ EAQR        :: EqualAI
            EAB + R ⇌ EABR        :: EqualAI
        end
    end
    a = @enzyme_mechanism begin
        substrates: A, B, C ; products: P, Q, R ; oligomeric_state: 1
        steps: begin
            E + A ⇌ EA            :: NonequalAI
            EA + B ⇌ EAB          :: EqualAI
            EAB + C ⇌ EABC        :: EqualAI
            EABC <--> EPQR        :: OnlyA
            EQR + P ⇌ EPQR        :: EqualAI
            ER + Q ⇌ EQR          :: EqualAI
            E + R ⇌ ER            :: OnlyA
            ER + A ⇌ EAR          :: EqualAI
            EQR + A ⇌ EAQR        :: EqualAI
            EAB + R ⇌ EABR        :: EqualAI
        end
    end
    b = @enzyme_mechanism begin
        substrates: A, B, C ; products: P, Q, R ; oligomeric_state: 1
        steps: begin
            E + A ⇌ EA            :: EqualAI
            EA + B ⇌ EAB          :: NonequalAI
            EAB + C ⇌ EABC        :: EqualAI
            EABC <--> EPQR        :: OnlyA
            EQR + P ⇌ EPQR        :: EqualAI
            ER + Q ⇌ EQR          :: EqualAI
            E + R ⇌ ER            :: OnlyA
            ER + A ⇌ EAR          :: EqualAI
            EQR + A ⇌ EAQR        :: EqualAI
            EAB + R ⇌ EABR        :: EqualAI
        end
    end
    c = @enzyme_mechanism begin
        substrates: A, B, C ; products: P, Q, R ; oligomeric_state: 1
        steps: begin
            E + A ⇌ EA            :: EqualAI
            EA + B ⇌ EAB          :: EqualAI
            EAB + C ⇌ EABC        :: NonequalAI
            EABC <--> EPQR        :: OnlyA
            EQR + P ⇌ EPQR        :: EqualAI
            ER + Q ⇌ EQR          :: EqualAI
            E + R ⇌ ER            :: OnlyA
            ER + A ⇌ EAR          :: EqualAI
            EQR + A ⇌ EAQR        :: EqualAI
            EAB + R ⇌ EABR        :: EqualAI
        end
    end
    r = @enzyme_mechanism begin
        substrates: A, B, C ; products: P, Q, R ; oligomeric_state: 1
        steps: begin
            E + A ⇌ EA            :: EqualAI
            EA + B ⇌ EAB          :: EqualAI
            EAB + C ⇌ EABC        :: EqualAI
            EABC <--> EPQR        :: OnlyA
            EQR + P ⇌ EPQR        :: EqualAI
            ER + Q ⇌ EQR          :: EqualAI
            E + R ⇌ ER            :: NonequalAI
            ER + A ⇌ EAR          :: EqualAI
            EQR + A ⇌ EAQR        :: EqualAI
            EAB + R ⇌ EABR        :: EqualAI
        end
    end
    r_at_eab = @enzyme_mechanism begin
        substrates: A, B, C ; products: P, Q, R ; oligomeric_state: 1
        steps: begin
            E + A ⇌ EA            :: EqualAI
            EA + B ⇌ EAB          :: EqualAI
            EAB + C ⇌ EABC        :: EqualAI
            EABC <--> EPQR        :: OnlyA
            EQR + P ⇌ EPQR        :: EqualAI
            ER + Q ⇌ EQR          :: EqualAI
            E + R ⇌ ER            :: OnlyA
            ER + A ⇌ EAR          :: EqualAI
            EQR + A ⇌ EAQR        :: EqualAI
            EAB + R ⇌ EABR        :: NonequalAI
        end
    end
    move = EnzymeRates.Enumeration._expand_change_allo_state
    _testhelper_move(move, parent)
    seconds = @elapsed kids = _testhelper_move(move, parent)
    @test seconds < 10
    @test length(kids) == 5
    @test Set(kids) == Set([a, b, c, r, r_at_eab])
end

@testset "Allosteric mechanism — regulator tag removal delta" begin
    # SEED: uni-uni allosteric with one regulator R tagged :OnlyA;
    # all 3 catalytic groups :EqualAI. The inactive conformation runs the whole cycle, so
    # one binding tagged :NonequalAI alone gets K_I = K_A from thermodynamics and adds
    # nothing: the chemistry alone, the two bindings together and R each add one.
    am = @enzyme_mechanism begin
        substrates: S; products: P
        oligomeric_state: 2
        steps: begin
            E + P ⇌ EP    :: EqualAI
            E + S ⇌ ES    :: EqualAI
            ES <--> EP  :: EqualAI
        end
        regulatory_site(2): R::OnlyA
    end
    _testhelper_assert_mechanism_invariants(am)
    chemistry = @enzyme_mechanism begin
        substrates: S; products: P
        oligomeric_state: 2
        steps: begin
            E + P ⇌ EP    :: EqualAI
            E + S ⇌ ES    :: EqualAI
            ES <--> EP  :: NonequalAI
        end
        regulatory_site(2): R::OnlyA
    end
    bindings = @enzyme_mechanism begin
        substrates: S; products: P
        oligomeric_state: 2
        steps: begin
            E + P ⇌ EP    :: NonequalAI
            E + S ⇌ ES    :: NonequalAI
            ES <--> EP  :: EqualAI
        end
        regulatory_site(2): R::OnlyA
    end
    ligand = @enzyme_mechanism begin
        substrates: S; products: P
        oligomeric_state: 2
        steps: begin
            E + P ⇌ EP    :: EqualAI
            E + S ⇌ ES    :: EqualAI
            ES <--> EP  :: EqualAI
        end
        regulatory_site(2): R::NonequalAI
    end

    result = _testhelper_move(EnzymeRates.Enumeration._expand_change_allo_state, am)

    # 1. the chemistry, the two bindings together, and the R ligand.
    @test length(result) == 3
    @test Set(result) == Set([chemistry, bindings, ligand])

    # 2. exactly one variant has the R ligand flipped to :NonequalAI
    # (allosteric_states preserved); locate it structurally.
    r_removal = filter(result) do r
        r.allosteric.allosteric_states == am.allosteric.allosteric_states &&
            any(s -> any(t -> t == :NonequalAI,
                         s.allosteric_states),
                r.allosteric.regulatory_sites)
    end
    @test length(r_removal) == 1

    # 3. ground-truth Δ params via `parameters`: the
    # R-ligand :OnlyA → :NonequalAI flip adds exactly one new
    # independent parameter (Δ=+1).
    seed_truth = _testhelper_fitted(am)
    @test _testhelper_fitted(only(r_removal)) == seed_truth + 1

    # 4. kind and invariants of each variant.
    for r in result
        @test _testhelper_is_allosteric(r)
        _testhelper_assert_mechanism_invariants(r)
    end
end

@testset "Allosteric mechanism — :OnlyI regulator-ligand relaxation" begin
    # SEED: uni-uni allosteric with R::OnlyI and 3 :EqualAI cat
    # groups: the chemistry alone, the two bindings together, and the R ligand.
    am = @enzyme_mechanism begin
        substrates: S
        products: P
        oligomeric_state: 2
        steps: begin
            E + P ⇌ EP      :: EqualAI
            E + S ⇌ ES      :: EqualAI
            ES <--> EP    :: EqualAI
        end
        regulatory_site(2): R::OnlyI
    end
    _testhelper_assert_mechanism_invariants(am)
    chemistry = @enzyme_mechanism begin
        substrates: S
        products: P
        oligomeric_state: 2
        steps: begin
            E + P ⇌ EP      :: EqualAI
            E + S ⇌ ES      :: EqualAI
            ES <--> EP    :: NonequalAI
        end
        regulatory_site(2): R::OnlyI
    end
    bindings = @enzyme_mechanism begin
        substrates: S
        products: P
        oligomeric_state: 2
        steps: begin
            E + P ⇌ EP      :: NonequalAI
            E + S ⇌ ES      :: NonequalAI
            ES <--> EP    :: EqualAI
        end
        regulatory_site(2): R::OnlyI
    end
    ligand = @enzyme_mechanism begin
        substrates: S
        products: P
        oligomeric_state: 2
        steps: begin
            E + P ⇌ EP      :: EqualAI
            E + S ⇌ ES      :: EqualAI
            ES <--> EP    :: EqualAI
        end
        regulatory_site(2): R::NonequalAI
    end

    result = _testhelper_move(EnzymeRates.Enumeration._expand_change_allo_state, am)

    # 1. the chemistry, the two bindings together, and the R ligand.
    @test length(result) == 3
    @test Set(result) == Set([chemistry, bindings, ligand])

    # 2. ground-truth Δ multiset via `parameters`. Under strict
    # `:EqualAI`, a single binding group tagged :NonequalAI while its
    # partners stay :EqualAI is degenerate: a thermodynamic cycle forbids
    # that affinity split, so it collapses (K_I = K_A) and adds no fitted
    # parameter, and the move tags none alone. The two bindings together, the
    # chemistry (its derived reverse absorbs the cycle) and the reg-ligand
    # relaxation each add one: `[1, 1, 1]`.
    @test _testhelper_param_deltas(am, result) == [1, 1, 1]

    # 3. invariants; the deltas above compile every result.
    foreach(_testhelper_assert_mechanism_invariants, result)

    # 4. exactly one ligand-relaxation variant: the R ligand's tag
    # flipped from :OnlyI to :NonequalAI (located structurally).
    n_r_relaxed = count(result) do r
        r.allosteric.allosteric_states == am.allosteric.allosteric_states &&
            any(s -> any(t -> t == :NonequalAI,
                         s.allosteric_states),
                r.allosteric.regulatory_sites)
    end
    @test n_r_relaxed == 1
end

@testset "Allosteric mechanism — multiple regulator ligands at independent tags" begin
    # SEED: allosteric uni-uni with R1::OnlyA and R2::OnlyI at two
    # regulatory sites, all 3 cat groups :EqualAI: the chemistry alone, the two
    # bindings together, and each ligand.
    am = @enzyme_mechanism begin
        substrates: S
        products: P
        oligomeric_state: 2
        steps: begin
            E + P ⇌ EP      :: EqualAI
            E + S ⇌ ES      :: EqualAI
            ES <--> EP    :: EqualAI
        end
        regulatory_site(2): R1::OnlyA
        regulatory_site(2): R2::OnlyI
    end
    _testhelper_assert_mechanism_invariants(am)
    chemistry = @enzyme_mechanism begin
        substrates: S
        products: P
        oligomeric_state: 2
        steps: begin
            E + P ⇌ EP      :: EqualAI
            E + S ⇌ ES      :: EqualAI
            ES <--> EP    :: NonequalAI
        end
        regulatory_site(2): R1::OnlyA
        regulatory_site(2): R2::OnlyI
    end
    bindings = @enzyme_mechanism begin
        substrates: S
        products: P
        oligomeric_state: 2
        steps: begin
            E + P ⇌ EP      :: NonequalAI
            E + S ⇌ ES      :: NonequalAI
            ES <--> EP    :: EqualAI
        end
        regulatory_site(2): R1::OnlyA
        regulatory_site(2): R2::OnlyI
    end
    r1 = @enzyme_mechanism begin
        substrates: S
        products: P
        oligomeric_state: 2
        steps: begin
            E + P ⇌ EP      :: EqualAI
            E + S ⇌ ES      :: EqualAI
            ES <--> EP    :: EqualAI
        end
        regulatory_site(2): R1::NonequalAI
        regulatory_site(2): R2::OnlyI
    end
    r2 = @enzyme_mechanism begin
        substrates: S
        products: P
        oligomeric_state: 2
        steps: begin
            E + P ⇌ EP      :: EqualAI
            E + S ⇌ ES      :: EqualAI
            ES <--> EP    :: EqualAI
        end
        regulatory_site(2): R1::OnlyA
        regulatory_site(2): R2::NonequalAI
    end

    result = _testhelper_move(EnzymeRates.Enumeration._expand_change_allo_state, am)

    # 1. the chemistry, the two bindings together, and each ligand.
    @test length(result) == 4
    @test Set(result) == Set([chemistry, bindings, r1, r2])

    # 2. ground-truth Δ multiset via `parameters`. Under strict
    # `:EqualAI`, a single binding group tagged :NonequalAI while its
    # partners stay :EqualAI collapses that affinity (K_I = K_A) and adds no
    # fitted parameter (a thermodynamic cycle forbids the split), and the move
    # tags none alone. The two bindings together, the chemistry and the two
    # reg-ligand relaxations each add one: `[1, 1, 1, 1]`.
    @test _testhelper_param_deltas(am, result) == [1, 1, 1, 1]

    # 3. invariants; the deltas above compile every result.
    foreach(_testhelper_assert_mechanism_invariants, result)

    # 4. property: each ligand has exactly one variant where ONLY it
    # is relaxed to :NonequalAI (cat states preserved, the other
    # ligand keeps its seed tag). Locate by per-ligand state.
    function _testhelper_site_lig_state(r, lig_name)
        for site in r.allosteric.regulatory_sites
            for (l, st) in zip(site.ligands,
                               site.allosteric_states)
                l === lig_name && return st
            end
        end
        error("ligand $lig_name not found")
    end
    n_r1_relaxed_only = count(result) do r
        r.allosteric.allosteric_states == am.allosteric.allosteric_states &&
            _testhelper_site_lig_state(r, :R1) == :NonequalAI &&
            _testhelper_site_lig_state(r, :R2) == :OnlyI
    end
    n_r2_relaxed_only = count(result) do r
        r.allosteric.allosteric_states == am.allosteric.allosteric_states &&
            _testhelper_site_lig_state(r, :R1) == :OnlyA &&
            _testhelper_site_lig_state(r, :R2) == :NonequalAI
    end
    @test n_r1_relaxed_only == 1
    @test n_r2_relaxed_only == 1
end

@testset "Allosteric mechanism — non-default site multiplicity (cat=4, reg=2)" begin
    # SEED: catalytic 4-mer with R::OnlyA at a multiplicity-2 reg site
    # (less than oligomeric_state).
    am = @enzyme_mechanism begin
        substrates: S
        products: P
        oligomeric_state: 4
        steps: begin
            E + P ⇌ EP      :: EqualAI
            E + S ⇌ ES      :: EqualAI
            ES <--> EP    :: EqualAI
        end
        regulatory_site(2): R::OnlyA
    end
    _testhelper_assert_mechanism_invariants(am)
    @test am.allosteric.oligomeric_state == 4
    @test [s.multiplicity for s in am.allosteric.regulatory_sites] == [2]

    # _expand_change_allo_state must preserve both
    # oligomeric_state and per-site multiplicity independently.
    result = _testhelper_move(EnzymeRates.Enumeration._expand_change_allo_state, am)
    @test !isempty(result)
    for r in result
        _testhelper_assert_mechanism_invariants(r)
        @test r.allosteric.oligomeric_state == 4
        @test [s.multiplicity for s in r.allosteric.regulatory_sites] == [2]
    end
end

@testset "Allosteric mechanism — uni-uni all-:EqualAI: bindings tag in pairs" begin
    # SEED: uni-uni allosteric with all 3 catalytic groups tagged
    # :EqualAI and no regulatory sites. The inactive conformation runs the whole cycle, so
    # either binding tagged :NonequalAI alone gets K_I = K_A and adds nothing; the
    # chemistry alone and the two bindings together each add one constant.
    rxn = @enzyme_reaction begin
        substrates: S[C]
        products: P[C]
        oligomeric_state: 2
    end
    am = @enzyme_mechanism begin
        substrates: S; products: P
        oligomeric_state: 2
        steps: begin
            E + P ⇌ EP    :: EqualAI
            E + S ⇌ ES    :: EqualAI
            ES <--> EP  :: EqualAI
        end
    end
    @test all(t -> t == :EqualAI, am.allosteric.allosteric_states)

    result = _testhelper_move(EnzymeRates.Enumeration._expand_change_allo_state, am)

    # 1. the chemistry, and the two bindings together; no regulatory ligand.
    expected = [
        (@enzyme_mechanism begin
            substrates: S; products: P
            oligomeric_state: 2
            steps: begin
                E + P ⇌ EP    :: NonequalAI
                E + S ⇌ ES    :: NonequalAI
                ES <--> EP  :: EqualAI
            end
        end),
        (@enzyme_mechanism begin
            substrates: S; products: P
            oligomeric_state: 2
            steps: begin
                E + P ⇌ EP    :: EqualAI
                E + S ⇌ ES    :: EqualAI
                ES <--> EP  :: NonequalAI
            end
        end),
    ]
    @test length(result) == 2
    @test Set(result) == Set(expected)

    # 2. a binding tagged alone has the seed's fitted count.
    p_alone = @enzyme_mechanism begin
        substrates: S; products: P
        oligomeric_state: 2
        steps: begin
            E + P ⇌ EP    :: NonequalAI
            E + S ⇌ ES    :: EqualAI
            ES <--> EP  :: EqualAI
        end
    end
    @test _testhelper_fitted(p_alone) == _testhelper_fitted(am)
    @test _testhelper_param_deltas(am, result) == [1, 1]

    # 3. preservation: reaction, multiplicity, sites.
    for r in result
        @test _testhelper_keeps_reaction(r, am)
        @test r.allosteric.oligomeric_state == am.allosteric.oligomeric_state
        @test r.allosteric.regulatory_sites == am.allosteric.regulatory_sites
    end
end

@testset "Allosteric mechanism — already-:NonequalAI: empty (negative)" begin
    # If every tag is :NonequalAI, no relaxation is possible.
    am_all_neq = @enzyme_mechanism begin
        substrates: S
        products: P
        oligomeric_state: 2
        steps: begin
            E + S ⇌ ES      :: NonequalAI
            ES <--> EP    :: NonequalAI
            EP ⇌ E + P      :: NonequalAI
        end
    end
    @test isempty(
        _testhelper_move(EnzymeRates.Enumeration._expand_change_allo_state, am_all_neq))
end

@testset "Allosteric mechanism — a relaxation that is not hyperbolic is not returned" begin
    # The split abortive ping-pong at two subunits (testset "_hyperbolic_catalysis"): the
    # second half-reaction and both abortive bindings, A's at EP_res_plusA_minusP and
    # P's at EA, :OnlyA. Each other binding tagged :NonequalAI keeps both laws
    # hyperbolic. An abortive binding tagged :NonequalAI brings its complex into the
    # inactive conformation, whose law then reaches A² or P² (`_hyperbolic_catalysis`), so
    # `_enumerator_returns` refuses both and candidate_mechanisms returns neither. The
    # predicate reads the graph alone: the guard rejects this parent and each of the six,
    # since the inactive conformation keeps the two sides of the rapid-equilibrium step
    # EB_res_plusA_minusP ⇌ EQ joined (`_FAST_STEP_PINNED`). The guard refuses every set
    # of groups the move's search builds, so the move has no child.
    am = @_testhelper_unguarded @enzyme_mechanism begin
        substrates: A, B ; products: P, Q ; oligomeric_state: 2
        steps: begin
            E + A ⇌ EA                                           :: EqualAI
            EP_res_plusA_minusP + A ⇌ EAP_res_plusA_minusP :: OnlyA
            E + Q ⇌ EQ                                           :: EqualAI
            E_res_plusA_minusP + B ⇌ EB_res_plusA_minusP     :: EqualAI
            E_res_plusA_minusP + P ⇌ EP_res_plusA_minusP     :: EqualAI
            EA + P ⇌ EAP                                     :: OnlyA
            EA <--> EP_res_plusA_minusP                       :: EqualAI
            EB_res_plusA_minusP ⇌ EQ                          :: OnlyA
        end
    end
    a_relaxed = @_testhelper_unguarded @enzyme_mechanism begin
        substrates: A, B ; products: P, Q ; oligomeric_state: 2
        steps: begin
            E + A ⇌ EA                                           :: NonequalAI
            EP_res_plusA_minusP + A ⇌ EAP_res_plusA_minusP :: OnlyA
            E + Q ⇌ EQ                                           :: EqualAI
            E_res_plusA_minusP + B ⇌ EB_res_plusA_minusP     :: EqualAI
            E_res_plusA_minusP + P ⇌ EP_res_plusA_minusP     :: EqualAI
            EA + P ⇌ EAP                                     :: OnlyA
            EA <--> EP_res_plusA_minusP                       :: EqualAI
            EB_res_plusA_minusP ⇌ EQ                          :: OnlyA
        end
    end
    q_relaxed = @_testhelper_unguarded @enzyme_mechanism begin
        substrates: A, B ; products: P, Q ; oligomeric_state: 2
        steps: begin
            E + A ⇌ EA                                           :: EqualAI
            EP_res_plusA_minusP + A ⇌ EAP_res_plusA_minusP :: OnlyA
            E + Q ⇌ EQ                                           :: NonequalAI
            E_res_plusA_minusP + B ⇌ EB_res_plusA_minusP     :: EqualAI
            E_res_plusA_minusP + P ⇌ EP_res_plusA_minusP     :: EqualAI
            EA + P ⇌ EAP                                     :: OnlyA
            EA <--> EP_res_plusA_minusP                       :: EqualAI
            EB_res_plusA_minusP ⇌ EQ                          :: OnlyA
        end
    end
    b_relaxed = @_testhelper_unguarded @enzyme_mechanism begin
        substrates: A, B ; products: P, Q ; oligomeric_state: 2
        steps: begin
            E + A ⇌ EA                                           :: EqualAI
            EP_res_plusA_minusP + A ⇌ EAP_res_plusA_minusP :: OnlyA
            E + Q ⇌ EQ                                           :: EqualAI
            E_res_plusA_minusP + B ⇌ EB_res_plusA_minusP     :: NonequalAI
            E_res_plusA_minusP + P ⇌ EP_res_plusA_minusP     :: EqualAI
            EA + P ⇌ EAP                                     :: OnlyA
            EA <--> EP_res_plusA_minusP                       :: EqualAI
            EB_res_plusA_minusP ⇌ EQ                          :: OnlyA
        end
    end
    p_relaxed = @_testhelper_unguarded @enzyme_mechanism begin
        substrates: A, B ; products: P, Q ; oligomeric_state: 2
        steps: begin
            E + A ⇌ EA                                           :: EqualAI
            EP_res_plusA_minusP + A ⇌ EAP_res_plusA_minusP :: OnlyA
            E + Q ⇌ EQ                                           :: EqualAI
            E_res_plusA_minusP + B ⇌ EB_res_plusA_minusP     :: EqualAI
            E_res_plusA_minusP + P ⇌ EP_res_plusA_minusP     :: NonequalAI
            EA + P ⇌ EAP                                     :: OnlyA
            EA <--> EP_res_plusA_minusP                       :: EqualAI
            EB_res_plusA_minusP ⇌ EQ                          :: OnlyA
        end
    end
    abortive_a_relaxed = @_testhelper_unguarded @enzyme_mechanism begin
        substrates: A, B ; products: P, Q ; oligomeric_state: 2
        steps: begin
            E + A ⇌ EA                                           :: EqualAI
            EP_res_plusA_minusP + A ⇌ EAP_res_plusA_minusP :: NonequalAI
            E + Q ⇌ EQ                                           :: EqualAI
            E_res_plusA_minusP + B ⇌ EB_res_plusA_minusP     :: EqualAI
            E_res_plusA_minusP + P ⇌ EP_res_plusA_minusP     :: EqualAI
            EA + P ⇌ EAP                                     :: OnlyA
            EA <--> EP_res_plusA_minusP                       :: EqualAI
            EB_res_plusA_minusP ⇌ EQ                          :: OnlyA
        end
    end
    abortive_p_relaxed = @_testhelper_unguarded @enzyme_mechanism begin
        substrates: A, B ; products: P, Q ; oligomeric_state: 2
        steps: begin
            E + A ⇌ EA                                           :: EqualAI
            EP_res_plusA_minusP + A ⇌ EAP_res_plusA_minusP :: OnlyA
            E + Q ⇌ EQ                                           :: EqualAI
            E_res_plusA_minusP + B ⇌ EB_res_plusA_minusP     :: EqualAI
            E_res_plusA_minusP + P ⇌ EP_res_plusA_minusP     :: EqualAI
            EA + P ⇌ EAP                                     :: NonequalAI
            EA <--> EP_res_plusA_minusP                       :: EqualAI
            EB_res_plusA_minusP ⇌ EQ                          :: OnlyA
        end
    end
    rxn = @enzyme_reaction begin
        substrates: A[CX], B[N]
        products: P[C], Q[NX]
        oligomeric_state: 2
    end
    hyperbolic = [a_relaxed, q_relaxed, b_relaxed, p_relaxed]
    kids = _testhelper_move(EnzymeRates.Enumeration._expand_change_allo_state, am)
    @test length(kids) == 0
    @test Set(kids) == Set()
    @test all(EnzymeRates.Enumeration._hyperbolic_catalysis,
              _testhelper_candidate.(hyperbolic))
    returned = _testhelper_children(rxn, [am])
    for d in (abortive_a_relaxed, abortive_p_relaxed)
        @test !EnzymeRates.Enumeration._hyperbolic_where_required(_testhelper_candidate(d))
        @test d ∉ returned
    end
    @test all(c -> EnzymeRates.Derivation._guard(c).fault ===
                   EnzymeRates.Derivation._FAST_STEP_PINNED,
              [am; hyperbolic; abortive_a_relaxed; abortive_p_relaxed])
end

end

# ─── Reachability of the allosteric versions ────────────────────────────

"""
    _testhelper_irreducible_versions(c, oligomeric_state)
        -> (; irreducible, keeps_reads_nonequal)

Every version of the plain candidate `c` at `oligomeric_state` with no regulatory site,
tried one by one: each kinetic group `:EqualAI`, `:OnlyA` or `:NonequalAI`, or also
`:OnlyI` when it binds a competitive inhibitor (`_bound_copy`), but not every group
`:EqualAI`. A mechanism shows when the enumerator returns it (`_enumerator_returns`), every
group that binds an inhibitor runs in a conformation (`_inhibitors_bind`) and none binds it
redundantly (`_redundant_copy_groups`). `irreducible` holds the versions that show and in
which making any one `:NonequalAI` group `:EqualAI` gives a mechanism with fewer fitted
constants (`Derivation.fitted_parameter_count`) or one that does not show. The mechanism
with every group `:EqualAI` is judged like the others: it fits `c`'s constants and `L`.
`keeps_reads_nonequal` holds the versions whose `_enumerator_keeps` answer changes when
every `:NonequalAI` group is made `:EqualAI`.
"""
function _testhelper_irreducible_versions(c, oligomeric_state)
    EN = EnzymeRates.Enumeration
    choices = [EN._bound_copy(c, first(ks)) === nothing ?
                   [:EqualAI, :OnlyA, :NonequalAI] :
                   [:EqualAI, :OnlyA, :NonequalAI, :OnlyI]
               for ks in EN._kinetic_groups(c.mechanism.catalytic)]
    no_site = EnzymeRates.Mechanisms.RegulatorySite[]
    mechanism = Dict{Vector{Symbol}, EnzymeRates.Mechanism}()
    keeps = Dict{Vector{Symbol}, Bool}()
    shows = Dict{Vector{Symbol}, Bool}()
    # One pass asks the guard of each version once: it remembers the last mechanism only.
    for states in Iterators.product(choices...)
        s = collect(states)
        v = EN._candidate(c, s, oligomeric_state, no_site)
        mechanism[s] = v.mechanism
        keeps[s] = EN._enumerator_keeps(v.mechanism)
        shows[s] = keeps[s] && EN._enumerator_returns(v) && EN._inhibitors_bind(v) &&
                   isempty(EN._redundant_copy_groups(v))
    end
    fitted = Dict{Vector{Symbol}, Int}()
    fitted_of(s) =
        get!(() -> EnzymeRates.Derivation.fitted_parameter_count(mechanism[s]), fitted, s)
    versions = filter(s -> any(!=(:EqualAI), s), collect(keys(mechanism)))
    keeps_reads_nonequal =
        filter(s -> keeps[s] != keeps[replace(s, :NonequalAI => :EqualAI)], versions)
    irreducible = Set{EnzymeRates.Mechanism}()
    for s in versions
        shows[s] || continue
        lowers(g) = (r = copy(s); r[g] = :EqualAI; fitted_of(r) < fitted_of(s) || !shows[r])
        all(lowers, findall(==(:NonequalAI), s)) && push!(irreducible, mechanism[s])
    end
    (; irreducible, keeps_reads_nonequal)
end

"""
    _testhelper_reached_versions(c, rxn) -> (; reached, no_gain, not_returned)

The candidates that the two allosteric moves reach from the plain candidate `c` of `rxn`,
keyed by their mechanisms: the children of `_expand_to_allosteric` that the enumerator
returns (`_enumerator_returns`), then, round by round, the returned children of
`_expand_change_allo_state` of each candidate reached. `no_gain` holds each pair
`parent => child` of a child, returned or not, whose fitted count is not above its
parent's, and `not_returned` each such pair of a child of `_expand_change_allo_state` that
the enumerator does not return.
"""
function _testhelper_reached_versions(c, rxn)
    EN = EnzymeRates.Enumeration
    fitted = Dict{EnzymeRates.Mechanism, Int}()
    fitted_of(m) = get!(() -> EnzymeRates.Derivation.fitted_parameter_count(m), fitted, m)
    no_gain = Pair{EnzymeRates.Mechanism, EnzymeRates.Mechanism}[]
    not_returned = Pair{EnzymeRates.Mechanism, EnzymeRates.Mechanism}[]
    reached = Dict{EnzymeRates.Mechanism, EN.CandidateMechanism}()
    level = [c]
    while !isempty(level)
        next = EN.CandidateMechanism[]
        for p in level
            children = p === c ? EN._expand_to_allosteric(c, rxn) :
                                 EN._expand_change_allo_state(p)
            for k in children
                fitted_of(k.mechanism) > fitted_of(p.mechanism) ||
                    push!(no_gain, p.mechanism => k.mechanism)
                if !EN._enumerator_returns(k)
                    p === c || push!(not_returned, p.mechanism => k.mechanism)
                    continue
                end
                haskey(reached, k.mechanism) && continue
                reached[k.mechanism] = k
                push!(next, k)
            end
        end
        level = next
    end
    (; reached, no_gain, not_returned)
end

@testset "to_allosteric and change_allo_state reach every irreducible version" begin
    # An aggregate check over samples of three populations of starting mechanisms, the
    # exception to writing every mechanism inline. For each sampled plain mechanism, at one
    # subunit with no regulator, every version is tried: each irreducible one, which the
    # enumerator would return and each of whose :NonequalAI groups adds a fitted constant,
    # is reached by `_expand_to_allosteric` and then `_expand_change_allo_state`; every
    # version reached binds each inhibitor in a conformation and none redundantly; every
    # child of the two moves, returned or not, has more fitted constants than its parent,
    # and the enumerator returns every child of `_expand_change_allo_state`; and
    # `_enumerator_keeps` reads the :OnlyA and :OnlyI groups alone, which the argument that
    # the moves reach every irreducible version assumes. The samples: every 7th of the 278
    # base-tier mechanisms of the bi-bi reaction (40), every 150th of the 1,548 that bind
    # a competitive inhibitor I (11), and every 3,200th of the 15,852 of at most seven
    # groups that bind I and a copy of A (5). The base-tier sample runs again at two
    # subunits, on its mechanisms that pass `_hyperbolic_catalysis` (10): the hyperbolic
    # rule (`_hyperbolic_where_required`), which the two searches for minimal gaining sets
    # read, can refuse a version only above one subunit. A group takes 3 states, an
    # inhibitor's 4, so seven groups, two of them inhibitors', give 3,887 versions; the
    # samples keep the testset under two minutes.
    EN = EnzymeRates.Enumeration
    plain_rxn = @enzyme_reaction begin
        substrates: A[CX], B[N]
        products: P[C], Q[NX]
    end
    two_subunit_rxn = @enzyme_reaction begin
        substrates: A[CX], B[N]
        products: P[C], Q[NX]
        oligomeric_state: 2
    end
    inhibitor_rxn = @enzyme_reaction begin
        substrates: A[CX], B[N]
        products: P[C], Q[NX]
        competitive_inhibitors: I
    end
    copy_rxn = @enzyme_reaction begin
        substrates: A[CX], B[N]
        products: P[C], Q[NX]
        competitive_inhibitors: I, A
    end
    at_most_seven_groups(c) = length(EN._kinetic_groups(c.mechanism.catalytic)) <= 7
    populations = [
        (plain_rxn, EN.candidate_mechanisms(plain_rxn), 7),
        (inhibitor_rxn, EN.candidate_mechanisms(inhibitor_rxn), 150),
        (copy_rxn, filter(at_most_seven_groups, EN.candidate_mechanisms(copy_rxn)), 3200)]
    @test [length(population) for (_, population, _) in populations] == [278, 1548, 15852]
    samples = [(rxn, population[1:step:end]) for (rxn, population, step) in populations]
    @test [length(parents) for (_, parents) in samples] == [40, 11, 5]
    (_, base_sample) = first(samples)
    hyperbolic = filter(EN._hyperbolic_catalysis, base_sample)
    @test length(hyperbolic) == 10
    runs = [[(rxn, 1, parents) for (rxn, parents) in samples];
            (two_subunit_rxn, 2, hyperbolic)]
    for (rxn, oligomeric_state, parents) in runs, c in parents
        versions = _testhelper_irreducible_versions(c, oligomeric_state)
        reach = _testhelper_reached_versions(c, rxn)
        @test !isempty(versions.irreducible)
        @test isempty(setdiff(versions.irreducible, keys(reach.reached)))
        @test all(EN._inhibitors_bind, values(reach.reached))
        @test all(k -> isempty(EN._redundant_copy_groups(k)), values(reach.reached))
        @test isempty(reach.no_gain)
        @test isempty(reach.not_returned)
        @test isempty(versions.keeps_reads_nonequal)
    end
end

@testset "change_allo_state gains a constant on mechanisms with two regulators" begin
    # An aggregate check over a sample of starting mechanisms, the exception to writing
    # every mechanism inline, and the only testset that runs the gain check on regulatory
    # ligands: every 4,000th of the 19,456 starting mechanisms that bind R and S at sites
    # of their own (5), and the children of `_expand_merge_regulatory_sites` that put the
    # two at one site and that the enumerator returns. Every child of
    # `_expand_change_allo_state`, a kinetic group's or a ligand's, has more fitted
    # constants than its parent. On this sample every ligand not yet :NonequalAI gives a
    # child, which adds one constant, so the check drops no ligand's child.
    EN = EnzymeRates.Enumeration
    fitted(c) = EnzymeRates.Derivation.fitted_parameter_count(c.mechanism)
    rxn = @enzyme_reaction begin
        substrates: A[CX], B[N]
        products: P[C], Q[NX]
        allosteric_regulators: R, S
    end
    starting = EN.candidate_mechanisms(rxn)
    @test length(starting) == 19456
    separate = starting[1:4000:end]
    shared = filter(EN._enumerator_returns,
                    [k for p in separate for k in EN._expand_merge_regulatory_sites(p)])
    @test length(separate) == 5
    @test all(p -> length(p.mechanism.allosteric.regulatory_sites) == 2, separate)
    @test !isempty(shared)
    @test all(p -> length(p.mechanism.allosteric.regulatory_sites) == 1, shared)
    parents = [separate; shared]
    children = [EN._expand_change_allo_state(p) for p in parents]
    no_gain = [p => k for (p, ks) in zip(parents, children) for k in ks
               if fitted(k) <= fitted(p)]
    @test isempty(no_gain)
    # A ligand's child keeps its parent's kinetic states; a kinetic group's changes them.
    ligand_children = [p => k for (p, ks) in zip(parents, children) for k in ks
                       if k.mechanism.allosteric.allosteric_states ==
                          p.mechanism.allosteric.allosteric_states]
    untagged = sum(count(!=(:NonequalAI), site.allosteric_states)
                   for p in parents for site in p.mechanism.allosteric.regulatory_sites)
    @test !isempty(ligand_children)
    @test length(ligand_children) == untagged
    @test all(((p, k),) -> fitted(k) == fitted(p) + 1, ligand_children)
end

# ─── _site_active_states ─────────────────────────────────────────────────
@testset "_site_active_states" begin
_testhelper_mk(states) = EnzymeRates.Mechanisms.RegulatorySite(
    [Symbol("R", i) for i in eachindex(states)],
    4, collect(Symbol, states))
@test EnzymeRates.Enumeration._site_active_states(_testhelper_mk([:OnlyA])) ==
      Set([:active])
@test EnzymeRates.Enumeration._site_active_states(_testhelper_mk([:OnlyI])) ==
      Set([:inactive])
@test EnzymeRates.Enumeration._site_active_states(_testhelper_mk([:EqualAI])) ==
      Set([:active, :inactive])
@test EnzymeRates.Enumeration._site_active_states(_testhelper_mk([:NonequalAI])) ==
      Set([:active, :inactive])
@test EnzymeRates.Enumeration._site_active_states(_testhelper_mk([:OnlyA, :OnlyA])) ==
      Set([:active])
@test EnzymeRates.Enumeration._site_active_states(_testhelper_mk([:OnlyA, :OnlyI])) ==
      Set([:active, :inactive])
end

# ─── _expand_merge_regulatory_sites ─────────────────────────────────────
@testset "_expand_merge_regulatory_sites" begin
# A four-subunit reaction with a designated activator and inhibitor.
merge_rxn = @enzyme_reaction begin
    substrates: S[C]
    products: P[C]
    allosteric_regulators: A::Activator, I::Inhibitor
    oligomeric_state: 4
end
parent = @enzyme_mechanism begin
    substrates: S; products: P
    oligomeric_state: 4
    steps: begin
        E + P ⇌ EP    :: OnlyA
        E + S ⇌ ES    :: OnlyA
        ES <--> EP  :: OnlyA
    end
    regulatory_site(4): A::OnlyA
    regulatory_site(4): I::OnlyI
end

children = _testhelper_move(EnzymeRates.Enumeration._expand_merge_regulatory_sites, parent)

# The allosteric state of ligand `lig` in a child's single merged site.
_testhelper_merged_state(child, lig) = begin
    site = only(child.allosteric.regulatory_sites)
    idx = findfirst(l -> l == lig,
                    site.ligands)
    site.allosteric_states[idx]
end
_testhelper_single_site_names(child) = Set(
    only(child.allosteric.regulatory_sites).ligands)

@testset "disjoint OnlyA/OnlyI co-binding skipped; antagonists kept" begin
    # A (:OnlyA) acts only on the active state, I (:OnlyI) only on the
    # inactive one, so co-binding them on one site derives to the same
    # equation as separate sites — that all-keep merge is redundant and
    # skipped. The two antagonist retags stay (each :EqualAI ligand now
    # acts on both states, a genuinely distinct mechanism).
    @test all(c -> _testhelper_single_site_names(c) == Set([:A, :I]), children)
    states = Set((_testhelper_merged_state(c, :A), _testhelper_merged_state(c, :I))
                 for c in children)
    @test !((:OnlyA, :OnlyI) in states)  # redundant co-binding skipped
    @test (:EqualAI, :OnlyI) in states   # activator → antagonist
    @test (:OnlyA, :EqualAI) in states   # inhibitor → antagonist
    @test !((:EqualAI, :EqualAI) in states)  # all-EqualAI dropped
    expected = [
        (@enzyme_mechanism begin
            substrates: S; products: P
            oligomeric_state: 4
            steps: begin
                E + P ⇌ EP    :: OnlyA
                E + S ⇌ ES    :: OnlyA
                ES <--> EP  :: OnlyA
            end
            regulatory_site(4): A::EqualAI, I::OnlyI
        end),
        @enzyme_mechanism begin
            substrates: S; products: P
            oligomeric_state: 4
            steps: begin
                E + P ⇌ EP    :: OnlyA
                E + S ⇌ ES    :: OnlyA
                ES <--> EP  :: OnlyA
            end
            regulatory_site(4): A::OnlyA, I::EqualAI
        end]
    @test length(children) == 2
    @test Set(children) == Set(expected)
end

@testset "every child is Δ0 (same fitted-param count as parent)" begin
    np_parent = _testhelper_fitted(parent)
    for c in children
        @test _testhelper_fitted(c) == np_parent
    end
end

@testset "parent and children are rate-equation-distinct" begin
    eqs = [EnzymeRates.rate_equation_string(parent);
           [EnzymeRates.rate_equation_string(c) for c in children]]
    @test length(unique(eqs)) == length(eqs)
end

@testset "reg_type filter keeps all when activator/inhibitor differ in type" begin
    # A::Activator + I::Inhibitor share a site with opposite types: no
    # ligand violates its type, so nothing is dropped — every merge child
    # survives.
    kept = filter(c -> EnzymeRates.Enumeration._respects_reg_type(c, merge_rxn), children)
    @test Set(kept) == Set(children)
    @test length(kept) == 2
end

@testset "reg_type filter drops same-type antagonist retags" begin
    # Two designated activators sharing one site: co-binding {OnlyA, OnlyA}
    # survives, but each antagonist retag places an :EqualAI beside a
    # same-type activator sibling, so both antagonist forms are dropped.
    rxn_aa = @enzyme_reaction begin
        substrates: S[C]
        products: P[C]
        allosteric_regulators: A1::Activator, A2::Activator
        oligomeric_state: 4
    end
    p_aa = @enzyme_mechanism begin
        substrates: S; products: P
        oligomeric_state: 4
        steps: begin
            E + P ⇌ EP    :: OnlyA
            E + S ⇌ ES    :: OnlyA
            ES <--> EP  :: OnlyA
        end
        regulatory_site(4): A1::OnlyA
        regulatory_site(4): A2::OnlyA
    end
    kids = _testhelper_move(EnzymeRates.Enumeration._expand_merge_regulatory_sites, p_aa)
    expected = [
        (@enzyme_mechanism begin
            substrates: S; products: P
            oligomeric_state: 4
            steps: begin
                E + P ⇌ EP    :: OnlyA
                E + S ⇌ ES    :: OnlyA
                ES <--> EP  :: OnlyA
            end
            regulatory_site(4): A1::OnlyA, A2::OnlyA
        end),
        (@enzyme_mechanism begin
            substrates: S; products: P
            oligomeric_state: 4
            steps: begin
                E + P ⇌ EP    :: OnlyA
                E + S ⇌ ES    :: OnlyA
                ES <--> EP  :: OnlyA
            end
            regulatory_site(4): A1::EqualAI, A2::OnlyA
        end),
        @enzyme_mechanism begin
            substrates: S; products: P
            oligomeric_state: 4
            steps: begin
                E + P ⇌ EP    :: OnlyA
                E + S ⇌ ES    :: OnlyA
                ES <--> EP  :: OnlyA
            end
            regulatory_site(4): A1::OnlyA, A2::EqualAI
        end]
    @test length(kids) == 3
    @test Set(kids) == Set(expected)
    kept = filter(c -> EnzymeRates.Enumeration._respects_reg_type(c, rxn_aa), kids)
    @test length(kept) == 1
    surviving = only(kept)
    @test Set(only(surviving.allosteric.regulatory_sites).allosteric_states) ==
          Set([:OnlyA])
end

@testset "different merge routes to one 3-way site dedup" begin
    rxn3 = @enzyme_reaction begin
        substrates: S[C]
        products: P[C]
        allosteric_regulators: A, B, C
        oligomeric_state: 4
    end
    # Route 1: {A,B} co-site + {C} single site → merge → {A,B,C}.
    p_ab_c = @enzyme_mechanism begin
        substrates: S; products: P
        oligomeric_state: 4
        steps: begin
            E + P ⇌ EP    :: OnlyA
            E + S ⇌ ES    :: OnlyA
            ES <--> EP  :: OnlyA
        end
        regulatory_site(4): A::OnlyA, B::OnlyA
        regulatory_site(4): C::OnlyA
    end
    # Route 2: {A,C} co-site + {B} single site → merge → {A,C,B}.
    p_ac_b = @enzyme_mechanism begin
        substrates: S; products: P
        oligomeric_state: 4
        steps: begin
            E + P ⇌ EP    :: OnlyA
            E + S ⇌ ES    :: OnlyA
            ES <--> EP  :: OnlyA
        end
        regulatory_site(4): A::OnlyA, C::OnlyA
        regulatory_site(4): B::OnlyA
    end
    # Both routes reach one site of all three ligands: co-binding, or one ligand
    # retagged :EqualAI.
    expected = [
        (@enzyme_mechanism begin
            substrates: S; products: P
            oligomeric_state: 4
            steps: begin
                E + P ⇌ EP    :: OnlyA
                E + S ⇌ ES    :: OnlyA
                ES <--> EP  :: OnlyA
            end
            regulatory_site(4): A::OnlyA, B::OnlyA, C::OnlyA
        end),
        (@enzyme_mechanism begin
            substrates: S; products: P
            oligomeric_state: 4
            steps: begin
                E + P ⇌ EP    :: OnlyA
                E + S ⇌ ES    :: OnlyA
                ES <--> EP  :: OnlyA
            end
            regulatory_site(4): A::EqualAI, B::OnlyA, C::OnlyA
        end),
        (@enzyme_mechanism begin
            substrates: S; products: P
            oligomeric_state: 4
            steps: begin
                E + P ⇌ EP    :: OnlyA
                E + S ⇌ ES    :: OnlyA
                ES <--> EP  :: OnlyA
            end
            regulatory_site(4): A::OnlyA, B::EqualAI, C::OnlyA
        end),
        @enzyme_mechanism begin
            substrates: S; products: P
            oligomeric_state: 4
            steps: begin
                E + P ⇌ EP    :: OnlyA
                E + S ⇌ ES    :: OnlyA
                ES <--> EP  :: OnlyA
            end
            regulatory_site(4): A::OnlyA, B::OnlyA, C::EqualAI
        end]
    kids_ab_c = _testhelper_move(EnzymeRates.Enumeration._expand_merge_regulatory_sites,
                                 p_ab_c)
    kids_ac_b = _testhelper_move(EnzymeRates.Enumeration._expand_merge_regulatory_sites,
                                 p_ac_b)
    @test length(kids_ab_c) == length(kids_ac_b) == 4
    @test Set(kids_ab_c) == Set(kids_ac_b) == Set(expected)
    _testhelper_threeway(kids) = only(filter(kids) do c
        sites = c.allosteric.regulatory_sites
        length(sites) == 1 &&
            length(only(sites).ligands) == 3 &&
            all(==(:OnlyA), only(sites).allosteric_states)
    end)
    c1 = _testhelper_threeway(kids_ab_c)
    c2 = _testhelper_threeway(kids_ac_b)
    @test c1 == c2
    @test hash(c1) == hash(c2)
end

@testset "ligand↔state pairing survives the name-sort (non-identity perm)" begin
    # Merge a 2-ligand site {A::OnlyA, C::OnlyI} with a 1-ligand site
    # {B::OnlyA}. The merged vcat order is ligands [A, C, B] / states
    # [OnlyA, OnlyI, OnlyA]; the RegulatorySite constructor's name-sort reorders
    # ligands to [A, B, C] (perm [1, 3, 2]). C's :OnlyI must ride along to the new
    # C slot — this only holds if the states are permuted in lockstep with the
    # ligands.
    rxn3 = @enzyme_reaction begin
        substrates: S[C]
        products: P[C]
        allosteric_regulators: A, B, C
        oligomeric_state: 4
    end
    p = @enzyme_mechanism begin
        substrates: S; products: P
        oligomeric_state: 4
        steps: begin
            E + P ⇌ EP    :: OnlyA
            E + S ⇌ ES    :: OnlyA
            ES <--> EP  :: OnlyA
        end
        regulatory_site(4): A::OnlyA, C::OnlyI
        regulatory_site(4): B::OnlyA
    end
    kids = _testhelper_move(EnzymeRates.Enumeration._expand_merge_regulatory_sites, p)
    expected = [
        (@enzyme_mechanism begin
            substrates: S; products: P
            oligomeric_state: 4
            steps: begin
                E + P ⇌ EP    :: OnlyA
                E + S ⇌ ES    :: OnlyA
                ES <--> EP  :: OnlyA
            end
            regulatory_site(4): A::OnlyA, B::OnlyA, C::OnlyI
        end),
        (@enzyme_mechanism begin
            substrates: S; products: P
            oligomeric_state: 4
            steps: begin
                E + P ⇌ EP    :: OnlyA
                E + S ⇌ ES    :: OnlyA
                ES <--> EP  :: OnlyA
            end
            regulatory_site(4): A::EqualAI, B::OnlyA, C::OnlyI
        end),
        (@enzyme_mechanism begin
            substrates: S; products: P
            oligomeric_state: 4
            steps: begin
                E + P ⇌ EP    :: OnlyA
                E + S ⇌ ES    :: OnlyA
                ES <--> EP  :: OnlyA
            end
            regulatory_site(4): A::OnlyA, B::EqualAI, C::OnlyI
        end),
        @enzyme_mechanism begin
            substrates: S; products: P
            oligomeric_state: 4
            steps: begin
                E + P ⇌ EP    :: OnlyA
                E + S ⇌ ES    :: OnlyA
                ES <--> EP  :: OnlyA
            end
            regulatory_site(4): A::OnlyA, B::OnlyA, C::EqualAI
        end]
    @test length(kids) == 4
    @test Set(kids) == Set(expected)
    # The co-binding child keeps every state (no :EqualAI).
    cobind = only(filter(kids) do c
        !(:EqualAI in only(c.allosteric.regulatory_sites).allosteric_states)
    end)
    @test _testhelper_merged_state(cobind, :A) == :OnlyA
    @test _testhelper_merged_state(cobind, :B) == :OnlyA
    @test _testhelper_merged_state(cobind, :C) == :OnlyI
end

@testset "no-op: single-site allosteric mechanism" begin
    single = @enzyme_mechanism begin
        substrates: S; products: P
        oligomeric_state: 4
        steps: begin
            E + P ⇌ EP    :: OnlyA
            E + S ⇌ ES    :: OnlyA
            ES <--> EP  :: OnlyA
        end
        regulatory_site(4): A::OnlyA
    end
    @test isempty(
        _testhelper_move(EnzymeRates.Enumeration._expand_merge_regulatory_sites, single))
end

@testset "OnlyA+OnlyI merge is derivation-redundant (premise guard)" begin
# The reason the all-keep OnlyA/OnlyI merge is skipped: merged onto one
# site it evaluates to the same rate as on separate sites. Confirm
# numerically over random parameters and concentrations.
separate = @enzyme_mechanism begin
    substrates: S; products: P
    oligomeric_state: 4
    steps: begin
        E + P ⇌ EP    :: OnlyA
        E + S ⇌ ES    :: OnlyA
        ES <--> EP  :: OnlyA
    end
    regulatory_site(4): A::OnlyA
    regulatory_site(4): I::OnlyI
end
merged = @enzyme_mechanism begin
    substrates: S; products: P
    oligomeric_state: 4
    steps: begin
        E + P ⇌ EP    :: OnlyA
        E + S ⇌ ES    :: OnlyA
        ES <--> EP  :: OnlyA
    end
    regulatory_site(4): A::OnlyA, I::OnlyI
end
fps = parameters(separate)
fpm = parameters(merged)
@test Set(fps) == Set(fpm)
separate_rate = EnzymeRates.compile_rate_equation(separate)
merged_rate = EnzymeRates.compile_rate_equation(merged)
Random.seed!(42)
for _ in 1:25
    vals = Dict(p => exp(randn()) for p in union(fps, fpm))
    ps = (; (p => vals[p] for p in fps)...)
    pm = (; (p => vals[p] for p in fpm)...)
    c = (; S=exp(randn()), P=exp(randn()), A=exp(randn()), I=exp(randn()))
    vs = real(separate_rate(c, ps, 100.0))
    vm = real(merged_rate(c, pm, 100.0))
    @test isapprox(vs, vm; rtol=1e-9)
end
end
end

@testset "candidate_mechanisms: a parent's children" begin
    # E + S → EP fused (SS), E + P ⇌ EP (RE). Flip: the P group cuts E from EP,
    # carries flux, no bottomless segment: one child, both steps SS. Split: single-step
    # groups, nothing. Dead end: no regulator. To-allosteric at multiplicity 1: each
    # nonempty set of the two groups :OnlyA. P's binding alone is refused: the inactive
    # conformation's fused step still joins its two sides (`_FAST_STEP_PINNED`). The
    # chemistry alone and both groups are returned; no set tagged :NonequalAI is refused
    # as :OnlyA. No regulator, so no site.
    parent = @enzyme_mechanism begin
        substrates: S
        products: P
        steps: begin
            E + S <--> EP
            E + P ⇌ EP
        end
    end
    flipped = @enzyme_mechanism begin
        substrates: S
        products: P
        steps: begin
            E + S <--> EP
            E + P <--> EP
        end
    end
    chemistry_only_a = @enzyme_mechanism begin
        substrates: S
        products: P
        steps: begin
            E + S <--> EP    :: OnlyA
            E + P ⇌ EP       :: EqualAI
        end
    end
    both_only_a = @enzyme_mechanism begin
        substrates: S
        products: P
        steps: begin
            E + S <--> EP    :: OnlyA
            E + P ⇌ EP       :: OnlyA
        end
    end
    children = EnzymeRates.Enumeration.candidate_mechanisms(
        uni_uni_rxn, _testhelper_candidate(parent))
    @test length(children) == 3
    @test Set(c.mechanism for c in children) ==
          Set([flipped, chemistry_only_a, both_only_a])
    # The limits drop children. Each child fits three constants. The flipped child's law
    # has 4 terms, two steady-state forms times the two steps that join them. With both
    # groups :OnlyA the active law has 1, one rapid-equilibrium segment, and the inactive
    # conformation runs no step; with the chemistry alone :OnlyA the inactive conformation
    # runs P's binding and adds a second term.
    @test isempty(EnzymeRates.Enumeration.candidate_mechanisms(
        uni_uni_rxn, _testhelper_candidate(parent); max_param_count = 2))
    @test [c.mechanism for c in EnzymeRates.Enumeration.candidate_mechanisms(
        uni_uni_rxn, _testhelper_candidate(parent); eq_complexity_filter = 2)] ==
          [chemistry_only_a, both_only_a]
end

# ─── a parent's children (candidate_mechanisms) ────────────────────────
@testset "candidate_mechanisms of a parent" begin

@testset "Mechanism — Returns candidates" begin
    # SEED: uni-uni RE-only, 3 singleton kinetic groups.
    m = @enzyme_mechanism begin
        substrates: S
        products: P
        steps: begin
            E + P ⇌ EP
            E + S ⇌ ES
            ES <--> EP
        end
    end
    _testhelper_assert_mechanism_invariants(m)
    result = EnzymeRates.Enumeration.candidate_mechanisms(uni_uni_rxn,
                                                          _testhelper_candidate(m))
    @test result isa Vector{EnzymeRates.Enumeration.CandidateMechanism}
    @test !isempty(result)
end

@testset "Mechanism — Allosteric expansion included" begin
    # SEED: uni-uni RE-only attached to an oligomeric reaction.
    # candidate_mechanisms must include allosteric mechanism variants in
    # its output.
    m = @enzyme_mechanism begin
        substrates: S
        products: P
        steps: begin
            E + P ⇌ EP
            E + S ⇌ ES
            ES <--> EP
        end
    end
    _testhelper_assert_mechanism_invariants(m)
    result = _testhelper_children(uni_uni_allo, [m])
    allo_count = count(s -> _testhelper_is_allosteric(s), result)
    # _expand_to_allosteric on this uni-uni seed (3 kinetic groups: 2
    # binding + 1 catalytic, no regulator declared) tags each nonempty set
    # of the groups :OnlyA. The exits keep 5 of the 7: one binding alone :OnlyA
    # is refused, since the inactive conformation still joins its two sides
    # (`_FAST_STEP_PINNED`). No set tagged :NonequalAI is refused as :OnlyA, so
    # there is no :NonequalAI child. Other moves do not produce allosteric
    # mechanism output from a plain mechanism, so exactly 5 exist.
    @test allo_count == 5
end

@testset "Mechanism — expansion never reduces param count" begin
    # SEED: uni-uni RE-only, base actual fitted count = 3. Expansion
    # moves never REDUCE the fitted-param count, so every child has
    # actual >= base. (The old estimate-based test asserted strictly
    # `> base`; that is FALSE for actual counts in general — some moves,
    # e.g. splitting a kinetic group whose peeled-off parameter is
    # already pinned by a Wegscheider cycle, leave the count UNCHANGED
    # at Δ=0. This particular uni-uni seed happens to add one parameter
    # per child, but the invariant we pin is the monotone `>=`.)
    m = @enzyme_mechanism begin
        substrates: S
        products: P
        steps: begin
            E + P ⇌ EP
            E + S ⇌ ES
            ES <--> EP
        end
    end
    _testhelper_assert_mechanism_invariants(m)
    base_fitted = _testhelper_fitted(m)
    result = _testhelper_children(uni_uni_rxn, [m])
    for child in result
        @test _testhelper_fitted(child) >= base_fitted
    end
end

@testset "Allosteric mechanism — rewrap preserves structure" begin
    # SEED: all-:EqualAI allosteric mechanism uni-uni. Passing this to
    # candidate_mechanisms must produce allosteric mechanism expansions
    # (RE→SS rewrapped as allosteric, etc.).
    am = @enzyme_mechanism begin
        substrates: S; products: P
        oligomeric_state: 2
        steps: begin
            E + P ⇌ EP    :: EqualAI
            E + S ⇌ ES    :: EqualAI
            ES <--> EP  :: EqualAI
        end
    end
    result = _testhelper_children(uni_uni_allo, [am])
    allo_results = filter(s -> _testhelper_is_allosteric(s),
                          result)
    @test !isempty(allo_results)
    # Every rewrapped allosteric result must preserve the input's
    # oligomeric_state and reaction — the catalytic steps may differ (a
    # base move may have changed them) but the allosteric-side
    # metadata is preserved.
    for r in allo_results
        @test r.allosteric.oligomeric_state == am.allosteric.oligomeric_state
        @test _testhelper_keeps_reaction(r, am)
    end
end

@testset "Allosteric mechanism — Dead-end excludes allosteric regs" begin
    # SEED: allosteric mechanism uni-uni with R already added as an
    # allosteric regulator (:OnlyA). candidate_mechanisms must never add R
    # as a dead-end inhibitor — no step in any expansion may bind R
    # (allosteric regulators live in regulatory_sites, not in the
    # catalytic steps).
    am = @enzyme_mechanism begin
        substrates: S; products: P
        oligomeric_state: 2
        steps: begin
            E + P ⇌ EP    :: EqualAI
            E + S ⇌ ES    :: EqualAI
            ES <--> EP  :: EqualAI
        end
        regulatory_site(2): R::OnlyA
    end
    result = _testhelper_children(uni_uni_allo_reg, [am])
    for r in result
        @test !any(s -> _testhelper_is_binding(s) && only(s.consumed) === :R,
                   _testhelper_graph(r).steps)
    end
end

@testset "Reg_type filter drops violating children; no-op when typeless" begin
    # True iff `m` places ligand `reg` in allosteric state `st` at any site.
    function _testhelper_has_reg_state(m, reg, st)
        _testhelper_is_allosteric(m) || return false
        for site in m.allosteric.regulatory_sites
            for (lig, s) in zip(site.ligands,
                                site.allosteric_states)
                lig == reg && s == st && return true
            end
        end
        false
    end

    rxn_act = @enzyme_reaction begin
        substrates: S[C]
        products: P[C]
        allosteric_regulators: R::Activator
        oligomeric_state: 2
    end
    rxn_plain = @enzyme_reaction begin
        substrates: S[C]
        products: P[C]
        allosteric_regulators: R
        oligomeric_state: 2
    end
    seed = @enzyme_mechanism begin
        substrates: S
        products: P
        steps: begin
            E + S ⇌ ES
            ES <--> EP
            EP ⇌ E + P
        end
    end

    # Unfiltered children reproduce candidate_mechanisms' moves on a Mechanism without its
    # filters: the regulator-type filter and `_enumerator_returns`.
    raw_act = _testhelper_plain_moves(seed, rxn_act)
    raw_plain = _testhelper_plain_moves(seed, rxn_plain)

    children_act = _testhelper_children(rxn_act, [seed])
    children_plain = _testhelper_children(rxn_plain, [seed])
    # The five sets of groups :OnlyA that the exits keep (one binding alone :OnlyA is
    # refused), and the chemistry alone :OnlyA with R at an :OnlyA site; the typeless
    # reaction adds the chemistry alone :OnlyA with R at an :OnlyI site.
    returned_act = [
        (@enzyme_mechanism begin
            substrates: S
            products: P
            oligomeric_state: 2
            steps: begin
                E + S ⇌ ES    :: OnlyA
                ES <--> EP  :: EqualAI
                EP ⇌ E + P    :: OnlyA
            end
        end),
        (@enzyme_mechanism begin
            substrates: S
            products: P
            oligomeric_state: 2
            steps: begin
                E + S ⇌ ES    :: EqualAI
                ES <--> EP  :: OnlyA
                EP ⇌ E + P    :: EqualAI
            end
        end),
        (@enzyme_mechanism begin
            substrates: S
            products: P
            oligomeric_state: 2
            steps: begin
                E + S ⇌ ES    :: EqualAI
                ES <--> EP  :: OnlyA
                EP ⇌ E + P    :: OnlyA
            end
        end),
        (@enzyme_mechanism begin
            substrates: S
            products: P
            oligomeric_state: 2
            steps: begin
                E + S ⇌ ES    :: OnlyA
                ES <--> EP  :: OnlyA
                EP ⇌ E + P    :: EqualAI
            end
        end),
        (@enzyme_mechanism begin
            substrates: S
            products: P
            oligomeric_state: 2
            steps: begin
                E + S ⇌ ES    :: OnlyA
                ES <--> EP  :: OnlyA
                EP ⇌ E + P    :: OnlyA
            end
        end),
        (@enzyme_mechanism begin
            substrates: S
            products: P
            oligomeric_state: 2
            steps: begin
                E + S ⇌ ES    :: EqualAI
                ES <--> EP  :: OnlyA
                EP ⇌ E + P    :: EqualAI
            end
            regulatory_site(2): R::OnlyA
        end),
    ]
    only_i = @enzyme_mechanism begin
        substrates: S
        products: P
        oligomeric_state: 2
        steps: begin
            E + S ⇌ ES    :: EqualAI
            ES <--> EP  :: OnlyA
            EP ⇌ E + P    :: EqualAI
        end
        regulatory_site(2): R::OnlyI
    end
    @test length(children_act) == 6
    @test Set(children_act) == Set(returned_act)
    @test length(children_plain) == 7
    @test Set(children_plain) == Set([returned_act; only_i])
    returned(m) = EnzymeRates.Enumeration._enumerator_returns(_testhelper_candidate(m))

    # A designated activator's raw expansion carries an :OnlyI R child; the reg_type
    # filter drops it, leaving strictly fewer.
    @test any(m -> _testhelper_has_reg_state(m, :R, :OnlyI), raw_act)
    @test !any(m -> _testhelper_has_reg_state(m, :R, :OnlyI), children_act)
    @test children_act ==
          filter(c -> EnzymeRates.Enumeration._respects_reg_type(c, rxn_act) &&
                      returned(c), raw_act)
    @test length(children_act) < length(raw_act)

    # A typeless reaction declares no type, so the filter is a no-op:
    # candidate_mechanisms returns exactly the raw children the enumerator returns,
    # :OnlyI included.
    @test any(m -> _testhelper_has_reg_state(m, :R, :OnlyI), children_plain)
    @test children_plain == filter(returned, raw_plain)
end

@testset "Merge move wired for two-site allosteric; no-op for Mechanism" begin
    # A two-regulatory-site Mechanism: A::Activator on one site,
    # I::Inhibitor on the other. candidate_mechanisms yields the Δ0
    # one-site merged children (both ligands sharing a single site).
    merge_rxn = @enzyme_reaction begin
        substrates: S[C]
        products: P[C]
        allosteric_regulators: A::Activator, I::Inhibitor
        oligomeric_state: 4
    end
    parent = @enzyme_mechanism begin
        substrates: S; products: P
        oligomeric_state: 4
        steps: begin
            E + P ⇌ EP    :: OnlyA
            E + S ⇌ ES    :: OnlyA
            ES <--> EP  :: OnlyA
        end
        regulatory_site(4): A::OnlyA
        regulatory_site(4): I::OnlyI
    end
    merged = _testhelper_move(EnzymeRates.Enumeration._expand_merge_regulatory_sites,
                              parent)
    expected = [
        (@enzyme_mechanism begin
            substrates: S; products: P
            oligomeric_state: 4
            steps: begin
                E + P ⇌ EP    :: OnlyA
                E + S ⇌ ES    :: OnlyA
                ES <--> EP  :: OnlyA
            end
            regulatory_site(4): A::EqualAI, I::OnlyI
        end),
        @enzyme_mechanism begin
            substrates: S; products: P
            oligomeric_state: 4
            steps: begin
                E + P ⇌ EP    :: OnlyA
                E + S ⇌ ES    :: OnlyA
                ES <--> EP  :: OnlyA
            end
            regulatory_site(4): A::OnlyA, I::EqualAI
        end]
    @test length(merged) == 2
    @test Set(merged) == Set(expected)

    children = _testhelper_children(merge_rxn, [parent])
    _testhelper_is_merged(c) = _testhelper_is_allosteric(c) &&
        length(c.allosteric.regulatory_sites) == 1 &&
        Set(only(c.allosteric.regulatory_sites).ligands) == Set([:A, :I])
    # The two Δ0 antagonist merge children survive the reg_type filter and
    # appear in the output; the redundant OnlyA/OnlyI co-binding is skipped.
    @test count(_testhelper_is_merged, children) == 2
    @test issubset(
        Set(_testhelper_move(EnzymeRates.Enumeration._expand_merge_regulatory_sites,
                             parent)),
        Set(children))

    # A plain Mechanism has no regulatory sites, so the merge move
    # contributes nothing: candidate_mechanisms equals the output of the moves that apply
    # to a Mechanism, filtered by regulator type and by `_enumerator_returns`.
    m = @enzyme_mechanism begin
        substrates: S; products: P
        steps: begin
            E + P ⇌ EP
            E + S ⇌ ES
            ES <--> EP
        end
    end
    raw = _testhelper_plain_moves(m, uni_uni_rxn)
    baseline = filter(raw) do c
        EnzymeRates.Enumeration._respects_reg_type(c, uni_uni_rxn) &&
            EnzymeRates.Enumeration._enumerator_returns(_testhelper_candidate(c))
    end
    @test _testhelper_children(uni_uni_rxn, [m]) == baseline
end

@testset "candidate_mechanisms reaches ≥2 distinct-metabolite catalytic :OnlyA" begin
    rxn = @enzyme_reaction begin
        substrates: A[C], B[N]
        products:   P[C1N1]
        allosteric_regulators: R
        oligomeric_state: 2
    end

    # #distinct metabolites bound by an :OnlyA catalytic group (iso → skip)
    _testhelper_onlya_mets(am) = (cm = am.catalytic;
        Set(cm.metabolites[only(cm.steps[first(ks)].taken_up)]
            for (g, ks) in enumerate(_testhelper_group_members(cm))
            if am.allosteric.allosteric_states[g] === :OnlyA &&
               !EnzymeRates.Mechanisms._is_iso(cm.steps[first(ks)])))

    seen = Set{UInt64}()
    frontier = EnzymeRates.Mechanism[]
    for m in _testhelper_mechanisms(EnzymeRates.Enumeration.init_mechanisms(rxn))
        h = hash(m); h in seen || (push!(seen, h); push!(frontier, m))
    end
    maxdistinct = 0
    gen = 0
    # The first generation reaches the target; three also run the regulator moves.
    while !isempty(frontier) && gen < 3
        nextf = eltype(frontier)[]
        for c in _testhelper_children(rxn, frontier)
            h = hash(c); h in seen || (push!(seen, h); push!(nextf, c))
            _testhelper_is_allosteric(c) &&
                (maxdistinct = max(maxdistinct, length(_testhelper_onlya_mets(c))))
        end
        frontier = nextf; gen += 1
    end

    @test maxdistinct >= 2
end
end

# ─── enumerate_all ─────────────────────────────────────────────────────
@testset "Integration" begin

uni = _testhelper_enumerate_all_mechanism(uni_uni_rxn; max_params=8)

@testset "Mechanism — Uni-uni full enumeration" begin
    # _testhelper_enumerate_all_mechanism buckets by ACTUAL fitted-parameter count. At
    # least 2 param-count buckets, and consecutive buckets separated by at
    # most 4 (max single-move delta).
    @test !isempty(uni)
    pcs = sort(collect(keys(uni)))
    @test length(pcs) >= 2
    @test all(pcs[i+1] - pcs[i] <= 4 for i in 1:length(pcs)-1)
end

@testset "Mechanism — Bi-bi init-tier (actual-count buckets)" begin
    # Full multi-tier bi-bi enumeration would have to compile every
    # reachable mechanism to bucket it by actual fitted count (hundreds
    # of @generated derivations) — too slow for the suite. The init tier
    # suffices to verify bi-bi enumeration produces mechanisms that
    # compile and that their actual fitted-param counts fall in the
    # expected {5,6,7} band. (Multi-tier actual-count enumeration is
    # exercised by the uni-uni / dead-end / allosteric callers below.)
    init = unique!(
        collect(_testhelper_mechanisms(EnzymeRates.Enumeration.init_mechanisms(bi_bi_rxn))))
    @test !isempty(init)
    _testhelper_holds_tc(m) = any(_testhelper_is_tc, m.catalytic.steps)
    n_ss(m) = count(fs -> !first(fs), _testhelper_group_flags(m))
    _testhelper_kind(m) = _testhelper_holds_iso(m) ? :seed :
                          _testhelper_holds_tc(m) ? :theorell_chance :
                          n_ss(m) == 2 ? :merged_two_groups : :merged_three_groups
    # The groups whose every step binds into a form that has no other step.
    function _testhelper_dead_end_groups(m)
        degree = zeros(Int, length(m.catalytic.enzyme_forms))
        for e in m.catalytic.steps
            degree[e.from] += 1; degree[e.to] += 1
        end
        count(ks -> all(k -> degree[m.catalytic.steps[k].to] == 1, ks),
              _testhelper_group_members(m.catalytic))
    end
    counts = Set{Int}()
    tally = Dict{Tuple{Symbol, Int}, Int}()
    for m in init
        n = _testhelper_fitted(m)
        push!(counts, n)
        tally[(_testhelper_kind(m), n)] = get(tally, (_testhelper_kind(m), n), 0) + 1
        _testhelper_kind(m) == :theorell_chance &&
            @test n == 5 + _testhelper_dead_end_groups(m)
    end
    # {5,6,7}: the seeds and their two-group merged variants fit 5, the
    # three-group merged variants 6, and a Theorell–Chance variant 5 plus one
    # for each group the eliminated steps leave holding only dead-end steps.
    @test counts == Set([5, 6, 7])
    @test tally == Dict((:seed, 5) => 55, (:merged_two_groups, 5) => 108,
                        (:merged_three_groups, 6) => 56, (:theorell_chance, 5) => 8,
                        (:theorell_chance, 6) => 8, (:theorell_chance, 7) => 4)
    @test all(m -> _testhelper_kind(m) == :theorell_chance ||
                   _testhelper_dead_end_groups(m) == 0, init)
end

@testset "Mechanism — With allosteric regulators" begin
    # A uni-uni allosteric reaction enumerates allosteric mechanisms.
    rxn = @enzyme_reaction begin
        substrates: S[C]
        products: P[C]
        allosteric_regulators: R
        oligomeric_state: 2
    end
    results = _testhelper_enumerate_all_mechanism(rxn; max_params=8)
    has_allo = any(
        any(_testhelper_is_allosteric(s) for s in mechs)
        for (_, mechs) in results)
    @test has_allo
end

@testset "Mechanism — With dead-end regulator" begin
    # Mechanism-form parallel. Total population with a dead-end
    # regulator strictly exceeds the plain uni-uni population at the
    # same cap, since the regulator opens additional expansion moves.
    rxn = @enzyme_reaction begin
        substrates: S[C]
        products: P[C]
        competitive_inhibitors: I
    end
    results = _testhelper_enumerate_all_mechanism(rxn; max_params=8)
    @test !isempty(results)
    total_with_reg = sum(length(v) for v in values(results))
    total_plain = sum(length(v) for v in values(uni))
    @test total_with_reg > total_plain
end
end

# ─── Expansion-move helpers ──────────────────────────────────────────────
# Helpers and building blocks of the RE→SS flip move and the kinetic-group
# split move. Every fixture is written out with the mechanism macro so the
# mechanism under test is visible in the testset.

@testset "_unbalanced_blocks" begin
    # Two parallel edges between the same two vertices, of different weight, close a
    # cycle of nonzero weight: one block, unbalanced.
    @test EnzymeRates.Enumeration._unbalanced_blocks(2, [(1, 2), (1, 2)], [1, 2]) ==
          ([1, 1], Set([1]))
    # A triangle whose third edge is stored against the walk 1 → 2 → 3: along its stored
    # direction 3 → 1 it weighs -3, so the cycle weighs 1 + 2 - 3 = 0 and is balanced.
    @test EnzymeRates.Enumeration._unbalanced_blocks(
        3, [(1, 2), (2, 3), (3, 1)], [1, 2, -3]) == ([1, 1, 1], Set{Int}())
    # The search reaches vertex 2 from vertex 1 over the edge stored as (2, 1), against
    # its direction, so the potential falls by its weight: (0, -1, 0). The cycle
    # 1 → 2 → 3 → 1 weighs -1 + 1 + 0 = 0, balanced.
    @test EnzymeRates.Enumeration._unbalanced_blocks(
        3, [(2, 1), (2, 3), (3, 1)], [1, 1, 0]) == ([1, 1, 1], Set{Int}())
    # Two components. Vertices 1 and 2 carry an edge and its reverse of opposite weight
    # (balanced). Vertices 3, 4 and 5 form a triangle whose paths 3 → 4 → 5 and 3 → 5
    # weigh 2 and 1 (unbalanced), with a bridge 5 → 6 hanging off it. Blocks are
    # numbered in the order the search closes them: the pair, the bridge, the triangle.
    @test EnzymeRates.Enumeration._unbalanced_blocks(
        6, [(1, 2), (2, 1), (3, 4), (4, 5), (3, 5), (5, 6)], [1, -1, 1, 1, 1, 4]) ==
          ([1, 1, 3, 3, 3, 2], Set([3]))
end

@testset "_flux_carrying_groups" begin
    # Every edge of `m` steady state.
    _testhelper_all_ss(m) = falses(length(_testhelper_flags(m)))
    # Ordered uni-uni: every step is on the catalytic cycle.
    m = @enzyme_mechanism begin
        substrates: S
        products: P
        steps: begin
            E + S ⇌ ES
            ES <--> EP
            E + P ⇌ EP
        end
    end
    @test all(EnzymeRates.Enumeration._flux_carrying_groups(m.catalytic,
                                                            _testhelper_all_ss(m)))

    # A dead-end leaf hanging off the cycle is a bridge: not flux-carrying.
    leaf = @enzyme_mechanism begin
        substrates: S
        products: P
        steps: begin
            E + S ⇌ ES
            ES <--> EP
            E + P ⇌ EP
            EP + S ⇌ EPS
        end
    end
    fc = EnzymeRates.Enumeration._flux_carrying_groups(leaf.catalytic,
                                                       _testhelper_all_ss(leaf))
    # The leaf group is the one whose step forms the doubly-bound EPS.
    @test !fc[leaf.catalytic.steps[_testhelper_edge(leaf, "EP + S ⇌ EPS")].kinetic_group]
    @test count(!, fc) == 1

    # Ping-pong: the second chemistry step starts at rapid equilibrium and lies
    # on the catalytic cycle, so every group is flux-carrying.
    pingpong = @enzyme_mechanism begin
        substrates: A, B
        products: P, Q
        steps: begin
            E + A ⇌ EA
            Estar + B ⇌ EstarB
            E + Q ⇌ EQ
            Estar + P ⇌ EstarP
            EA <--> EstarP
            EstarB ⇌ EQ
        end
    end
    @test all(EnzymeRates.Enumeration._flux_carrying_groups(pingpong.catalytic,
                                                            _testhelper_all_ss(pingpong)))

    # An inhibitor bound to both E and ES, with the mirror EI+S ⇌ EIS: the
    # inhibitor square shares the edge E→ES with the catalytic cycle, so the
    # cycle E→EI→EIS→ES→EP→E passes through chemistry and every step,
    # inhibitor binding included, carries net flux. Contrast the pendant square
    # below, whose inhibitors bind only free E.
    mirror = EnzymeRates.@enzyme_mechanism begin
        substrates: S
        products: P
        modifiers: I
        steps: begin
            E + S ⇌ ES          :: EqualAI
            ES <--> EP        :: EqualAI
            E + P ⇌ EP          :: EqualAI
            E + I ⇌ EI          :: EqualAI
            EI + S ⇌ EIS    :: EqualAI
            ES + I ⇌ EIS    :: EqualAI
        end
    end
    @test all(EnzymeRates.Enumeration._flux_carrying_groups(
        mirror.catalytic, _testhelper_all_ss(mirror.catalytic)))

    # Two inhibitors binding only free E form a binding-only square joined to the
    # cycle at the single form E: no cycle through it contains chemistry, so its
    # four groups are not flux-carrying.
    pendant = EnzymeRates.@enzyme_mechanism begin
        substrates: S
        products: P
        modifiers: I, J
        steps: begin
            E + S ⇌ ES          :: EqualAI
            ES <--> EP        :: EqualAI
            E + P ⇌ EP          :: EqualAI
            E + I ⇌ EI          :: EqualAI
            E + J ⇌ EJ          :: EqualAI
            EI + J ⇌ EIJ    :: EqualAI
            EJ + I ⇌ EIJ    :: EqualAI
        end
    end
    cp = pendant.catalytic
    pfc = EnzymeRates.Enumeration._flux_carrying_groups(cp, _testhelper_all_ss(cp))
    inhibitor_groups =
        unique(e.kinetic_group for e in cp.steps
               if EnzymeRates.Enumeration._is_binding(e) &&
                   cp.metabolites[only(e.taken_up)] in (:I, :J))
    @test length(inhibitor_groups) == 4
    @test !any(pfc[g] for g in inhibitor_groups)
    @test all(pfc[g] for g in eachindex(pfc) if !(g in inhibitor_groups))

    # Ping-pong: every group lies on the single catalytic cycle, so every
    # group is flux-carrying.
    pp = @enzyme_mechanism begin
        substrates: A, B
        products: P, Q
        steps: begin
            E + A ⇌ EA
            EA <--> EP_res_plusA_minusP
            EP_res_plusA_minusP ⇌ E_res_plusA_minusP + P
            E_res_plusA_minusP + B ⇌ EB_res_plusA_minusP
            EB_res_plusA_minusP <--> EQ
            EQ ⇌ E + Q
        end
    end
    @test all(EnzymeRates.Enumeration._flux_carrying_groups(pp.catalytic,
                                                            _testhelper_all_ss(pp)))
    # A Theorell–Chance step consumes B and releases P in one edge, so the
    # mechanism has no isomerization step. The cycle E → EA → EQ → E has
    # weight 1 + 2 + 1 = 4 = one turnover times the four reactants, so it is
    # unbalanced and every step carries flux.
    tc = @enzyme_mechanism begin
        substrates: A, B
        products: P, Q
        steps: begin
            E + A ⇌ EA
            EA + B <--> EQ + P
            E + Q ⇌ EQ
        end
    end
    @test all(EnzymeRates.Enumeration._flux_carrying_groups(tc.catalytic,
                                                            _testhelper_all_ss(tc)))
    @test EnzymeRates.Enumeration._flux_carrying_groups(
        tc.catalytic, _testhelper_flags(tc)) ==
          [!first(fs) for fs in _testhelper_group_flags(tc)]

    # A merged central complex X with two fused steps into it and no
    # isomerization: E + A ⇌ EA, EA + B → X, EQ + P → X, E + Q ⇌ EQ.
    # Every step lies on the unbalanced cycle E → EA → X → EQ → E.
    merged = @enzyme_mechanism begin
        substrates: A, B
        products: P, Q
        steps: begin
            E + A ⇌ EA
            EA + B <--> EPQ
            EQ + P <--> EPQ
            E + Q ⇌ EQ
        end
    end
    @test all(
        EnzymeRates.Enumeration._flux_carrying_groups(merged.catalytic,
                                                      _testhelper_all_ss(merged)))

    # An RE shunt. E + A and E + Q are steady-state while the abortive route
    # EQ + A ⇌ EAQ ⇌ EA + Q stays at rapid equilibrium. The segment
    # {E} joins the rest by two parallel edges of weight 0 (A uptake +1 against
    # EA's offset +1; Q uptake −1 against EQ's offset −1): a balanced block,
    # so neither steady-state binding carries flux. The turnover runs inside the
    # big segment, where the chemistry step is an unbalanced self-loop of
    # weight 4 (EAB sits at A + B, EPQ at −P − Q). On the all-steady-state
    # graph every group carries flux.
    shunt = @enzyme_mechanism begin
        substrates: A, B
        products: P, Q
        steps: begin
            E + A <--> EA
            EQ + A ⇌ EAQ
            E + Q <--> EQ
            EA + Q ⇌ EAQ
            EA + B ⇌ EAB
            EAB <--> EPQ
            EQ + P ⇌ EPQ
        end
    end
    # Each group's first edge.
    _testhelper_firsts(m) =
        [m.catalytic.steps[first(ks)] for ks in _testhelper_group_members(m.catalytic)]
    ss_binding(e) = EnzymeRates.Enumeration._is_binding(e) && !e.rapid_equilibrium
    flags = EnzymeRates.Enumeration._flux_carrying_groups(
        shunt.catalytic, _testhelper_flags(shunt))
    for (g, e) in enumerate(_testhelper_firsts(shunt))
        if ss_binding(e)
            @test String(shunt.catalytic.enzyme_forms[e.from]) == "E" &&
                  !flags[g]
        elseif EnzymeRates.Mechanisms._is_iso(e)
            @test flags[g]
        else
            @test e.rapid_equilibrium && !flags[g]
        end
    end
    @test count(ss_binding, _testhelper_firsts(shunt)) == 2
    @test all(
        EnzymeRates.Enumeration._flux_carrying_groups(shunt.catalytic,
                                                      _testhelper_all_ss(shunt)))

    # A balanced self-loop: E + A is steady-state, but E and EA already share
    # a segment through E ⇌ EQ ⇌ EAQ ⇌ EA, so the step is a self-loop
    # of weight 1 + 0 − 1 = 0 and carries no flux. The chemistry self-loop in
    # the same segment has weight 4 and does.
    loop = @_testhelper_unguarded @enzyme_mechanism begin
        substrates: A, B
        products: P, Q
        steps: begin
            E + A <--> EA
            EQ + A ⇌ EAQ
            (E + Q ⇌ EQ, EA + Q ⇌ EAQ)
            EA + B ⇌ EAB
            EAB <--> EPQ
            EQ + P ⇌ EPQ
        end
    end
    # At zero products its segment's steady state depends on how fast the
    # rapid-equilibrium steps run (class 3); the flux test reads the graph alone.
    v = EnzymeRates.Derivation._guard(loop)
    @test (v.fault, v.limit, v.face) ===
          (EnzymeRates.Derivation._SPLIT, EnzymeRates.Derivation._LIMIT_SPEED, 2)
    @test EnzymeRates.Mechanisms.step_text(loop.catalytic, v.edge) == "E + Q ⇌ EQ"
    @test _testhelper_respeed_verdicts(loop; rng = Random.MersenneTwister(92),
                                       law = :limit)[2] === :pattern
    lflags = EnzymeRates.Enumeration._flux_carrying_groups(
        loop.catalytic, _testhelper_flags(loop))
    for (g, e) in enumerate(_testhelper_firsts(loop))
        if ss_binding(e)
            @test !lflags[g]
        elseif EnzymeRates.Mechanisms._is_iso(e)
            @test lflags[g]
        end
    end
    @test count(ss_binding, _testhelper_firsts(loop)) == 1
    @test EnzymeRates.Enumeration._re_segment_count(
        loop.catalytic, _testhelper_flags(loop)) == 1
end

@testset "_re_segment_count" begin
    m = @enzyme_mechanism begin
        substrates: S
        products: P
        steps: begin
            E + S ⇌ ES
            ES <--> EP
            E + P ⇌ EP
        end
    end
    # Flipping E + S ⇌ ES cuts the one segment in two.
    @test EnzymeRates.Enumeration._re_segment_count(
        m.catalytic, _testhelper_flags(m)) == 1
    flipped = _testhelper_flags(m)
    flipped[_testhelper_edge(m, "E + S ⇌ ES")] = false
    @test EnzymeRates.Enumeration._re_segment_count(m.catalytic, flipped) == 2

    # Allosteric: measured on its catalytic groups.
    am = EnzymeRates.@enzyme_mechanism begin
        substrates: S
        products: P
        oligomeric_state: 2
        steps: begin
            E + S ⇌ ES   :: EqualAI
            ES <--> EP :: EqualAI
            E + P ⇌ EP   :: EqualAI
        end
    end
    @test _testhelper_segment_count(am) == 1
    am_flipped = _testhelper_flags(am.catalytic)
    am_flipped[_testhelper_edge(am, "E + S ⇌ ES")] = false
    @test EnzymeRates.Enumeration._re_segment_count(am.catalytic, am_flipped) == 2

    # Aggregate regression pin over the bi-bi seed set: on every single RE→SS group flip
    # the count agrees with the segmentation the test reader's steps of the flipped graph
    # give (`_testhelper_re_segments`).
    flips = [(seed, g)
             for seed in _testhelper_mechanisms(
                 EnzymeRates.Enumeration.init_mechanisms(bi_bi_rxn))
             for (g, ks) in enumerate(_testhelper_group_members(seed.catalytic))
             if all(k -> seed.catalytic.steps[k].rapid_equilibrium, ks)]
    @test length(flips) == 508
    @test all(flips) do (seed, g)
        re = _testhelper_flags(seed)
        re[[e.kinetic_group == g for e in seed.catalytic.steps]] .= false
        cm = seed.catalytic
        flipped = _testhelper_plain_of(EnzymeRates.Mechanisms._canonical(
            cm.metabolites, cm.enzyme_forms, cm.substrates, cm.products,
            [EnzymeRates.Mechanisms.Step(e.from, e.to, e.taken_up, e.given_off, re[k],
                                         e.kinetic_group)
             for (k, e) in enumerate(cm.steps)]))
        EnzymeRates.Enumeration._re_segment_count(seed.catalytic, re) ==
            length(first(_testhelper_re_segments(_testhelper_graph(flipped))))
    end
end

@testset "_minimal_gaining_sets" begin
    # Units 1..4. Sets gain iff they contain {1,2} or contain 3.
    gains(set) = (1 in set && 2 in set) || 3 in set
    sets = EnzymeRates.Enumeration._minimal_gaining_sets(gains, _ -> 1:4)
    @test sets == [[3], [1, 2]]
    # Partner pruning: unit 2 may never join unit 1, so {1,2} is unreachable.
    partners(set) = (1 in set || 2 in set) ? [3, 4] : 1:4
    @test EnzymeRates.Enumeration._minimal_gaining_sets(gains, partners) == [[3]]
    # Nothing gains: empty result, and the search terminates.
    @test isempty(EnzymeRates.Enumeration._minimal_gaining_sets(_ -> false, _ -> 1:3))
    # No units at all.
    @test isempty(EnzymeRates.Enumeration._minimal_gaining_sets(_ -> true, _ -> 1:0))
end

"Whole-group flip of groups `gs` of `m` to steady state (test helper)."
function _testhelper_flip_groups(m, gs)
    re = _testhelper_flags(m)
    re[[e.kinetic_group in gs for e in m.catalytic.steps]] .= false
    EnzymeRates.Enumeration._with_flags(_testhelper_candidate(m), re).mechanism
end

@testset "_context_bipartitions" begin
    # A binds E, EB, and EP: contexts B and P give two bipartitions.
    m = @enzyme_mechanism begin
        substrates: A, B
        products: P, Q
        steps: begin
            (E + A ⇌ EA, EB + A ⇌ EAB, EP + A ⇌ EAP)
            (E + B ⇌ EB, EA + B ⇌ EAB)
            (E + P ⇌ EP, EQ + P ⇌ EPQ, EA + P ⇌ EAP)
            (E + Q ⇌ EQ, EP + Q ⇌ EPQ)
            EAB <--> EPQ
        end
    end
    # The names of the forms the edges of `part`, a part of graph `ctx`, start from.
    _testhelper_src_forms(part, ctx) =
        Set(ctx.catalytic.enzyme_forms[ctx.catalytic.steps[k].from]
            for k in part)
    # The kinetic group of the edge of `ctx` whose text is `text`.
    _testhelper_group_of(ctx, text) =
        ctx.catalytic.steps[_testhelper_edge(ctx, text)].kinetic_group
    ga = _testhelper_group_of(m, "E + A ⇌ EA")
    a_group = findall(e -> e.kinetic_group == ga, m.catalytic.steps)
    @test length(a_group) == 3
    bps = EnzymeRates.Enumeration._context_bipartitions(_testhelper_candidate(m), ga)
    @test length(bps) == 2
    for (with, without) in bps
        @test !isempty(with) && !isempty(without)
        @test length(with) + length(without) == 3
        @test first(a_group) in with
        @test isempty(intersect(with, without))
    end
    # Contexts are ordered by ligand role (product before substrate) then
    # name: bps[1] is the P-context division, bps[2] the B-context division.
    @test bps[1][2] == [_testhelper_edge(m, "EP + A ⇌ EAP")]
    @test bps[2][2] == [_testhelper_edge(m, "EB + A ⇌ EAB")]
    # By P: the P-free forms E, EB against the P-bound form EP.
    @test _testhelper_src_forms(bps[1][1], m) == Set([:E, :EB])
    @test _testhelper_src_forms(bps[1][2], m) == Set([:EP])
    # By B: the B-free forms E, EP against the B-bound form EB.
    @test _testhelper_src_forms(bps[2][1], m) == Set([:E, :EP])
    @test _testhelper_src_forms(bps[2][2], m) == Set([:EB])
    # _context_form: a binding is stored taking its metabolite up, so the context form
    # is its source.
    @test EnzymeRates.Enumeration._context_form(m.catalytic, first(a_group)) ==
          m.catalytic.steps[first(a_group)].from
    # _context_form: a steady-state step that gives off one metabolite found in neither
    # end's bound list (the Segel ping-pong step shape) is stored as the binding it
    # reverses, F + P → EA, so its context form is F, the form P binds to.
    segel = @enzyme_mechanism begin
        substrates: A
        products: P
        steps: begin
            E + A ⇌ EA
            EA <--> F + P
            F <--> E
        end
    end
    segel_k = _testhelper_edge(segel, "F + P <--> EA")
    @test EnzymeRates.Enumeration._context_form(segel.catalytic, segel_k) ==
          segel.catalytic.steps[segel_k].from ==
          only(f for f in eachindex(segel.catalytic.enzyme_forms)
               if String(segel.catalytic.enzyme_forms[f]) == "F")
    # A two-step group with one context has one bipartition; a group whose
    # source forms carry no other ligand has none.
    @test length(EnzymeRates.Enumeration._context_bipartitions(_testhelper_candidate(m),
                                                               _testhelper_group_of(m,
                                                                     "E + B ⇌ EB"))) == 1
    @test isempty(EnzymeRates.Enumeration._context_bipartitions(_testhelper_candidate(m),
                                                                _testhelper_group_of(m,
                                                                     "EAB <--> EPQ")))

    # Ter-ter substrate cube: each of the two other substrates divides the A
    # group's four steps two and two.
    cube = @enzyme_mechanism begin
        substrates: A, B, C
        products: P, Q, R
        steps: begin
            (E + A ⇌ EA, EB + A ⇌ EAB, EC + A ⇌ EAC,
             EBC + A ⇌ EABC)
            (E + B ⇌ EB, EA + B ⇌ EAB, EC + B ⇌ EBC,
             EAC + B ⇌ EABC)
            (E + C ⇌ EC, EA + C ⇌ EAC, EB + C ⇌ EBC,
             EAB + C ⇌ EABC)
            EABC <--> EPQR
            EQR + P ⇌ EPQR
            ER + Q ⇌ EQR
            E + R ⇌ ER
        end
    end
    cube_bps = EnzymeRates.Enumeration._context_bipartitions(
        _testhelper_candidate(cube), _testhelper_group_of(cube, "E + A ⇌ EA"))
    @test length(cube_bps) == 2
    # Both contexts are substrates, so they are ordered by name: B then C.
    # By B: the B-free forms E, EC against the B-bound forms EB, EBC.
    @test _testhelper_src_forms(cube_bps[1][1], cube) == Set([:E, :EC])
    @test _testhelper_src_forms(cube_bps[1][2], cube) == Set([:EB, :EBC])
    # By C: the C-free forms E, EB against the C-bound forms EC, EBC.
    @test _testhelper_src_forms(cube_bps[2][1], cube) == Set([:E, :EB])
    @test _testhelper_src_forms(cube_bps[2][2], cube) == Set([:EC, :EBC])

    # Residual as a context: Q binds free E, the modified enzyme, and the
    # modified enzyme with B bound. Context B divides {E, E(;res)} from
    # {E(B;res)}; the residual divides {E} from the two modified forms.
    res = @enzyme_mechanism begin
        substrates: A, B
        products: P, Q
        steps: begin
            (E + Q ⇌ EQ, E_res_plusA_minusP + Q ⇌ EQ_res_plusA_minusP,
             EB_res_plusA_minusP + Q ⇌ EBQ_res_plusA_minusP)
            E + A ⇌ EA
            EA <--> EP_res_plusA_minusP
            EP_res_plusA_minusP ⇌ E_res_plusA_minusP + P
            E_res_plusA_minusP + B ⇌ EB_res_plusA_minusP
            EB_res_plusA_minusP <--> EQ
        end
    end
    rbps = EnzymeRates.Enumeration._context_bipartitions(
        _testhelper_candidate(res), _testhelper_group_of(res, "E + Q ⇌ EQ"))
    @test length(rbps) == 2
    Er, EBr = :E_res_plusA_minusP, :EB_res_plusA_minusP
    # By B, then by the residual.
    @test _testhelper_src_forms(rbps[1][1], res) == Set([:E, Er]) &&
          _testhelper_src_forms(rbps[1][2], res) == Set([EBr])
    @test _testhelper_src_forms(rbps[2][1], res) == Set([:E]) &&
          _testhelper_src_forms(rbps[2][2], res) == Set([Er, EBr])

    # Residual contexts come in content order: by the names a residual adds, then by
    # those it subtracts. Q binds at E, at E_res_plusA_minusP and at
    # E_res_plusA_plusB_minusP, so each residual parts one binding from the other two:
    # no residual first, then A - P, then A + B - P.
    residues = @enzyme_mechanism begin
        substrates: A, B
        products: P, Q
        steps: begin
            (E + Q ⇌ EQ, E_res_plusA_minusP + Q ⇌ EQ_res_plusA_minusP,
             E_res_plusA_plusB_minusP + Q ⇌ EQ_res_plusA_plusB_minusP)
            E + A ⇌ EA
            EA <--> EP_res_plusA_minusP
            EP_res_plusA_minusP ⇌ E_res_plusA_minusP + P
            E_res_plusA_minusP + B ⇌ EB_res_plusA_minusP
            EB_res_plusA_minusP <--> E_res_plusA_plusB_minusP
            E_res_plusA_plusB_minusP <--> EQ
        end
    end
    none, a_p, ab_p = [_testhelper_edge(residues, t) for t in (
        "E + Q ⇌ EQ", "E_res_plusA_minusP + Q ⇌ EQ_res_plusA_minusP",
        "E_res_plusA_plusB_minusP + Q ⇌ EQ_res_plusA_plusB_minusP")]
    @test EnzymeRates.Enumeration._context_bipartitions(
              _testhelper_candidate(residues),
              residues.catalytic.steps[none].kinetic_group) ==
          [([none], [a_p, ab_p]), ([none, ab_p], [a_p]), ([none, a_p], [ab_p])]

    # Conformation as a context: A binds E, Estar, and EstarB, which B's binding to
    # Estar forms.
    conf = @enzyme_mechanism begin
        substrates: A, B
        products: P, Q
        steps: begin
            (E + A ⇌ EA, Estar + A ⇌ EstarA, EstarB + A ⇌ EstarAB)
            Estar + B ⇌ EstarB
            EA <--> EstarA
            EstarAB <--> EstarPQ
            EstarPQ ⇌ EstarP + Q
            EstarP ⇌ E + P
        end
    end
    cbps = EnzymeRates.Enumeration._context_bipartitions(_testhelper_candidate(conf),
                                                         _testhelper_group_of(conf,
                                                                      "E + A ⇌ EA"))
    @test length(cbps) == 2
    # By B, then by the conformation.
    @test _testhelper_src_forms(cbps[1][1], conf) == Set([:E, :Estar])
    @test _testhelper_src_forms(cbps[1][2], conf) == Set([:EstarB])
    @test _testhelper_src_forms(cbps[2][1], conf) == Set([:E])
    @test _testhelper_src_forms(cbps[2][2], conf) == Set([:Estar, :EstarB])
end

@testset "_context_bipartitions separates an inhibitor-bound mirror" begin
    am = EnzymeRates.@enzyme_mechanism begin
        substrates: S
        products: P
        modifiers: I
        steps: begin
            (E + S ⇌ ES, EI + S ⇌ EIS)   :: EqualAI
            ES <--> EP                        :: EqualAI
            E + P ⇌ EP                          :: EqualAI
            E + I ⇌ EI                          :: EqualAI
        end
    end
    cm = am.catalytic
    bps = EnzymeRates.Enumeration._context_bipartitions(
        _testhelper_candidate(cm),
        cm.steps[_testhelper_edge(cm, "E + S ⇌ ES")].kinetic_group)
    @test bps == [([_testhelper_edge(cm, "E + S ⇌ ES")],
                   [_testhelper_edge(cm, "EI + S ⇌ EIS")])]
end

@testset "_with_flags" begin
    m = @enzyme_mechanism begin
        substrates: S
        products: P
        steps: begin
            E + S ⇌ ES
            ES <--> EP
            EP ⇌ E + P
        end
    end
    # Canonical form stores the steps E + P → EP, E + S → ES, ES → EP.
    @test EnzymeRates.Enumeration._with_flags(_testhelper_candidate(m),
                                              BitVector([1, 0, 0])).mechanism ==
          @enzyme_mechanism begin
              substrates: S
              products: P
              steps: begin
                  E + S <--> ES
                  ES <--> EP
                  EP ⇌ E + P
              end
          end
    @test EnzymeRates.Enumeration._with_flags(_testhelper_candidate(m),
                                              BitVector([1, 1, 0])).mechanism == m

    am = @enzyme_mechanism begin
        substrates: S
        products: P
        steps: begin
            E + S ⇌ ES      :: NonequalAI
            ES <--> EP    :: OnlyA
            EP ⇌ E + P      :: EqualAI
        end
    end
    @test EnzymeRates.Enumeration._with_flags(_testhelper_candidate(am),
                                              BitVector([1, 0, 0])).mechanism ==
          @enzyme_mechanism begin
              substrates: S
              products: P
              steps: begin
                  E + S <--> ES   :: NonequalAI
                  ES <--> EP    :: OnlyA
                  EP ⇌ E + P      :: EqualAI
              end
          end
end

@testset "_regrouped" begin
    flags = _testhelper_flags
    m = @enzyme_mechanism begin
        substrates: A, B
        products: P, Q
        steps: begin
            (E + A ⇌ EA, EB + A ⇌ EAB)
            (E + B ⇌ EB, EA + B ⇌ EAB)
            EAB <--> EPQ
            EPQ <--> EQ + P
            EQ <--> E + Q
        end
    end
    # Canonical form stores the steps E + A → EA, E + B → EB, E + Q → EQ,
    # EA + B → EAB, EAB → EPQ, EB + A → EAB, EQ + P → EPQ, in
    # the groups [1, 2, 3, 2, 4, 1, 5]. The groups below put EB + A in a group of its own.
    split = [1, 2, 3, 2, 4, 6, 5]
    @test EnzymeRates.Enumeration._regrouped(_testhelper_candidate(m), split,
                                             flags(m)).mechanism ==
          @enzyme_mechanism begin
              substrates: A, B
              products: P, Q
              steps: begin
                  E + A ⇌ EA
                  EB + A ⇌ EAB
                  (E + B ⇌ EB, EA + B ⇌ EAB)
                  EAB <--> EPQ
                  EPQ <--> EQ + P
                  EQ <--> E + Q
              end
          end
    # New groups and new flags at once: EB + A also turns steady state.
    @test EnzymeRates.Enumeration._regrouped(_testhelper_candidate(m), split,
                                             BitVector([1, 1, 0, 1, 0, 0, 0])).mechanism ==
          @enzyme_mechanism begin
              substrates: A, B
              products: P, Q
              steps: begin
                  E + A ⇌ EA
                  EB + A <--> EAB
                  (E + B ⇌ EB, EA + B ⇌ EAB)
                  EAB <--> EPQ
                  EPQ <--> EQ + P
                  EQ <--> E + Q
              end
          end
    @test EnzymeRates.Enumeration._regrouped(
        _testhelper_candidate(m), [e.kinetic_group for e in m.catalytic.steps],
        flags(m)).mechanism == m

    # Each new group of an allosteric mechanism takes the tag of its first step's group.
    am = @enzyme_mechanism begin
        substrates: A, B
        products: P, Q
        steps: begin
            (E + A ⇌ EA, EB + A ⇌ EAB)  :: NonequalAI
            (E + B ⇌ EB, EA + B ⇌ EAB)  :: EqualAI
            EAB <--> EPQ                :: OnlyA
            EPQ <--> EQ + P               :: EqualAI
            EQ <--> E + Q                     :: EqualAI
        end
    end
    @test EnzymeRates.Enumeration._regrouped(_testhelper_candidate(am), split,
                                             flags(am.catalytic)).mechanism ==
          @enzyme_mechanism begin
              substrates: A, B
              products: P, Q
              steps: begin
                  E + A ⇌ EA                        :: NonequalAI
                  EB + A ⇌ EAB                  :: NonequalAI
                  (E + B ⇌ EB, EA + B ⇌ EAB)  :: EqualAI
                  EAB <--> EPQ                :: OnlyA
                  EPQ <--> EQ + P               :: EqualAI
                  EQ <--> E + Q                     :: EqualAI
              end
          end
end

@testset "_apply_bipartitions" begin
    m = @enzyme_mechanism begin
        substrates: A, B
        products: P, Q
        steps: begin
            (E + A ⇌ EA, EB + A ⇌ EAB)
            (E + B ⇌ EB, EA + B ⇌ EAB)
            (E + P ⇌ EP, EQ + P ⇌ EPQ)
            (E + Q ⇌ EQ, EP + Q ⇌ EPQ)
            EAB <--> EPQ
        end
    end
    g = m.catalytic.steps[_testhelper_edge(m, "E + A ⇌ EA")].kinetic_group
    bp = only(EnzymeRates.Enumeration._context_bipartitions(_testhelper_candidate(m), g))
    child = EnzymeRates.Enumeration._apply_bipartitions(_testhelper_candidate(m),
                                                        [(g, bp)]).mechanism
    # The steps stay; only the group of the second part is new.
    @test length(_testhelper_group_members(child.catalytic)) ==
          length(_testhelper_group_members(m.catalytic)) + 1
    @test child.catalytic.enzyme_forms == m.catalytic.enzyme_forms &&
          _testhelper_ungrouped(child) == _testhelper_ungrouped(m)
    _testhelper_members(c, k) =
        findall(e -> e.kinetic_group == c.catalytic.steps[k].kinetic_group,
                c.catalytic.steps)
    @test _testhelper_members(child, first(bp[1])) == bp[1]
    @test _testhelper_members(child, first(bp[2])) == bp[2]

    am = @enzyme_mechanism begin
        substrates: A, B
        products: P, Q
        oligomeric_state: 2
        steps: begin
            (E + A ⇌ EA, EB + A ⇌ EAB)    :: NonequalAI
            (E + B ⇌ EB, EA + B ⇌ EAB)    :: EqualAI
            E + P ⇌ EP             :: EqualAI
            EP + Q ⇌ EPQ       :: EqualAI
            E + Q ⇌ EQ             :: EqualAI
            EQ + P ⇌ EPQ       :: EqualAI
            EAB <--> EPQ     :: EqualAI
        end
    end
    cm = am.catalytic
    ga = cm.steps[_testhelper_edge(cm, "E + A ⇌ EA")].kinetic_group
    bpa = only(EnzymeRates.Enumeration._context_bipartitions(_testhelper_candidate(cm), ga))
    achild = EnzymeRates.Enumeration._apply_bipartitions(_testhelper_candidate(am),
                                                         [(ga, bpa)]).mechanism
    @test _testhelper_ungrouped(achild.catalytic) == _testhelper_ungrouped(cm)
    @test [achild.allosteric.allosteric_states[
               achild.catalytic.steps[first(part)].kinetic_group]
           for part in bpa] == [:NonequalAI, :NonequalAI]
    @test achild.allosteric.oligomeric_state == 2
    @test achild.allosteric.regulatory_sites ==
          am.allosteric.regulatory_sites
end

@testset "_revert_zero_flux_parts" begin
    # Random-order product release with B also bound abortively at EQ and at EP,
    # the two abortive bindings one steady-state group. Neither abortive step carries
    # flux, so with both flags false both parts of the group's context bipartition
    # revert to rapid equilibrium; with both flags true neither does.
    m = @enzyme_mechanism begin
        substrates: A, B
        products: P, Q
        steps: begin
            E + A ⇌ EA
            EA + B ⇌ EAB
            EAB <--> EPQ
            EQ + P ⇌ EPQ
            EP + Q ⇌ EPQ
            E + P ⇌ EP
            E + Q ⇌ EQ
            (EQ + B <--> EBQ, EP + B <--> EBP)
        end
    end
    abortive = [_testhelper_edge(m, "EQ + B <--> EBQ"),
                _testhelper_edge(m, "EP + B <--> EBP")]
    g = only(unique(m.catalytic.steps[k].kinetic_group for k in abortive))
    @test sort(findall(e -> e.kinetic_group == g, m.catalytic.steps)) == sort(abortive)
    @test !any(k -> m.catalytic.steps[k].rapid_equilibrium, abortive)
    @test EnzymeRates.Enumeration._flux_carrying_edges(
        m.catalytic, _testhelper_flags(m))[abortive] == falses(2)
    bp = only(EnzymeRates.Enumeration._context_bipartitions(_testhelper_candidate(m), g))
    @test all(part -> length(part) == 1, bp)
    reverted =
        EnzymeRates.Enumeration._revert_zero_flux_parts(m.catalytic, bp,
                                                        falses(length(m.catalytic.steps)))
    @test sort(reverted) == sort(abortive)
    @test isempty(
        EnzymeRates.Enumeration._revert_zero_flux_parts(m.catalytic, bp,
                                                        trues(length(m.catalytic.steps))))
    # The child built with those flags holds both parts at rapid equilibrium, each its
    # forms and lists unchanged.
    re = _testhelper_flags(m)
    re[reverted] .= true
    child = EnzymeRates.Enumeration._apply_bipartitions(_testhelper_candidate(m), [(g, bp)],
                                                        re).mechanism
    @test [(e.from, e.to, e.taken_up, e.given_off) for e in child.catalytic.steps] ==
          [(e.from, e.to, e.taken_up, e.given_off) for e in m.catalytic.steps]
    @test all(k -> child.catalytic.steps[k].rapid_equilibrium, abortive)
end

"Closure of `seed` under `gen`, by structural identity."
function _testhelper_closure(seed, gen; maxn = 5_000)
    seen = Dict{UInt64, Any}(hash(seed) => seed); queue = Any[seed]
    while !isempty(queue)
        m = popfirst!(queue)
        for c in gen(m)
            h = hash(c)
            haskey(seen, h) && continue
            seen[h] = c; push!(queue, c)
            length(seen) > maxn && error("closure exceeded $maxn")
        end
    end
    collect(values(seen))
end

@testset "_expand_re_to_ss (group-set flips)" begin
    @testset "_expand_re_to_ss: seed child count (204 over bi-bi seeds)" begin
        # Aggregate regression pin over the whole bi-bi init set: the 55 seeds, each
        # holding its isomerization, and their 184 merged and Theorell–Chance variants.
        # 8 of the 55 seeds hold a qualifying chain whose isomerization is steady state;
        # its two flanks never flip (`_chain_flank_groups`), and with them the seeds
        # would have 220 flip children. The variants hold no isomerization, so no chain,
        # and each flip child of a variant flips one group: each of the 108 two-group
        # merged variants holds two flip units, each of the 56 three-group ones one, and
        # the 20 Theorell–Chance ones none, so the variants have 2 × 108 + 56 = 272.
        init = _testhelper_mechanisms(EnzymeRates.Enumeration.init_mechanisms(bi_bi_rxn))
        seeds = filter(_testhelper_holds_iso, init)
        variants = filter(!_testhelper_holds_iso, init)
        @test (length(seeds), length(variants)) == (55, 184)
        @test count(
            m -> !isempty(EnzymeRates.Enumeration._chain_flank_groups(m)), init) == 8
        @test sum(length(_testhelper_move(EnzymeRates.Enumeration._expand_re_to_ss, m))
                  for m in seeds) == 204
        @test sum(length(_testhelper_move(EnzymeRates.Enumeration._expand_re_to_ss, m))
                  for m in variants) == 272
    end

    @testset "_expand_re_to_ss: uni-uni emits no flip, both bindings are chain flanks" begin
        # Uni-uni steady-state and rapid-equilibrium laws have the same form: ES →
        # EP is a qualifying chain whose flanks are the S and P groups, so neither
        # flips (`_chain_flank_groups`). The S flip, absent, keeps the parent's rank; the
        # P flip is its mirror under reversing the reaction.
        m = @enzyme_mechanism begin
            substrates: S
            products: P
            steps: begin
                E + S ⇌ ES
                ES <--> EP
                E + P ⇌ EP
            end
        end
        absent = @enzyme_mechanism begin
            substrates: S
            products: P
            steps: begin
                E + S <--> ES
                ES <--> EP
                E + P ⇌ EP
            end
        end
        kids = _testhelper_move(EnzymeRates.Enumeration._expand_re_to_ss, m)
        @test EnzymeRates.Enumeration._chain_flank_groups(m) ==
              Set(g for (g, ks) in enumerate(_testhelper_group_members(m.catalytic))
                  if !EnzymeRates.Mechanisms._is_iso(m.catalytic.steps[only(ks)]))
        @test isempty(kids)
        @test _testhelper_identifiable_rank(absent) == _testhelper_identifiable_rank(m)
    end

    @testset "_expand_re_to_ss: random-order bi-bi emits one flip per metabolite" begin
        # Every metabolite's two binding steps share one group, so its single group cuts a
        # segment alone, raising the segment count on its own; each metabolite is emitted
        # alone and no pair is minimal.
        m = @enzyme_mechanism begin
            substrates: A, B
            products: P, Q
            steps: begin
                (E + A ⇌ EA, EB + A ⇌ EAB)
                (E + B ⇌ EB, EA + B ⇌ EAB)
                (E + P ⇌ EP, EQ + P ⇌ EPQ)
                (E + Q ⇌ EQ, EP + Q ⇌ EPQ)
                EAB <--> EPQ
            end
        end
        flipA = @enzyme_mechanism begin
            substrates: A, B
            products: P, Q
            steps: begin
                (E + A <--> EA, EB + A <--> EAB)
                (E + B ⇌ EB, EA + B ⇌ EAB)
                (E + P ⇌ EP, EQ + P ⇌ EPQ)
                (E + Q ⇌ EQ, EP + Q ⇌ EPQ)
                EAB <--> EPQ
            end
        end
        flipB = @enzyme_mechanism begin
            substrates: A, B
            products: P, Q
            steps: begin
                (E + A ⇌ EA, EB + A ⇌ EAB)
                (E + B <--> EB, EA + B <--> EAB)
                (E + P ⇌ EP, EQ + P ⇌ EPQ)
                (E + Q ⇌ EQ, EP + Q ⇌ EPQ)
                EAB <--> EPQ
            end
        end
        flipP = @enzyme_mechanism begin
            substrates: A, B
            products: P, Q
            steps: begin
                (E + A ⇌ EA, EB + A ⇌ EAB)
                (E + B ⇌ EB, EA + B ⇌ EAB)
                (E + P <--> EP, EQ + P <--> EPQ)
                (E + Q ⇌ EQ, EP + Q ⇌ EPQ)
                EAB <--> EPQ
            end
        end
        flipQ = @enzyme_mechanism begin
            substrates: A, B
            products: P, Q
            steps: begin
                (E + A ⇌ EA, EB + A ⇌ EAB)
                (E + B ⇌ EB, EA + B ⇌ EAB)
                (E + P ⇌ EP, EQ + P ⇌ EPQ)
                (E + Q <--> EQ, EP + Q <--> EPQ)
                EAB <--> EPQ
            end
        end
        kids = _testhelper_move(EnzymeRates.Enumeration._expand_re_to_ss, m)
        @test length(kids) == 4
        @test Set(kids) == Set([flipA, flipB, flipP, flipQ])
    end

    @testset "_expand_re_to_ss: a dead-end leaf group never flips" begin
        # EP + S ⇌ EPS is a bridge: no cycle through it contains chemistry,
        # so at steady state it carries no net flux and flipping it would add a
        # parameter the data cannot see. Only the two catalytic bindings flip.
        m = @enzyme_mechanism begin
            substrates: S
            products: P
            steps: begin
                E + S ⇌ ES
                ES <--> EP
                E + P ⇌ EP
                EP + S ⇌ EPS
            end
        end
        expected = [
            (@enzyme_mechanism begin
                substrates: S
                products: P
                steps: begin
                    E + S <--> ES
                    ES <--> EP
                    E + P ⇌ EP
                    EP + S ⇌ EPS
                end
            end),
            (@enzyme_mechanism begin
                substrates: S
                products: P
                steps: begin
                    E + S ⇌ ES
                    ES <--> EP
                    E + P <--> EP
                    EP + S ⇌ EPS
                end
            end),
        ]
        kids = _testhelper_move(EnzymeRates.Enumeration._expand_re_to_ss, m)
        @test length(kids) == 2
        @test Set(kids) == Set(expected)
    end

    @testset "_expand_re_to_ss: a split metabolite flips as a pair, never singly" begin
        # A's two binding steps sit in separate groups. Flipping either alone
        # leaves E and EA joined through the other A step (E–EB–EAB–EA),
        # so the segment count does not rise and no such child exists. The pair
        # cuts the square, so it is emitted; B, P, and Q each flip alone, since each
        # cuts a segment on its own.
        m = @enzyme_mechanism begin
            substrates: A, B
            products: P, Q
            steps: begin
                E + A ⇌ EA
                EB + A ⇌ EAB
                (E + B ⇌ EB, EA + B ⇌ EAB)
                (E + P ⇌ EP, EQ + P ⇌ EPQ)
                (E + Q ⇌ EQ, EP + Q ⇌ EPQ)
                EAB <--> EPQ
            end
        end
        flipA = @enzyme_mechanism begin
            substrates: A, B
            products: P, Q
            steps: begin
                E + A <--> EA
                EB + A <--> EAB
                (E + B ⇌ EB, EA + B ⇌ EAB)
                (E + P ⇌ EP, EQ + P ⇌ EPQ)
                (E + Q ⇌ EQ, EP + Q ⇌ EPQ)
                EAB <--> EPQ
            end
        end
        flipB = @enzyme_mechanism begin
            substrates: A, B
            products: P, Q
            steps: begin
                E + A ⇌ EA
                EB + A ⇌ EAB
                (E + B <--> EB, EA + B <--> EAB)
                (E + P ⇌ EP, EQ + P ⇌ EPQ)
                (E + Q ⇌ EQ, EP + Q ⇌ EPQ)
                EAB <--> EPQ
            end
        end
        flipP = @enzyme_mechanism begin
            substrates: A, B
            products: P, Q
            steps: begin
                E + A ⇌ EA
                EB + A ⇌ EAB
                (E + B ⇌ EB, EA + B ⇌ EAB)
                (E + P <--> EP, EQ + P <--> EPQ)
                (E + Q ⇌ EQ, EP + Q ⇌ EPQ)
                EAB <--> EPQ
            end
        end
        flipQ = @enzyme_mechanism begin
            substrates: A, B
            products: P, Q
            steps: begin
                E + A ⇌ EA
                EB + A ⇌ EAB
                (E + B ⇌ EB, EA + B ⇌ EAB)
                (E + P ⇌ EP, EQ + P ⇌ EPQ)
                (E + Q <--> EQ, EP + Q <--> EPQ)
                EAB <--> EPQ
            end
        end
        kids = _testhelper_move(EnzymeRates.Enumeration._expand_re_to_ss, m)
        @test length(kids) == 4
        @test Set(kids) == Set([flipA, flipB, flipP, flipQ])
        # The two single-A flips are absent: each is segment-flat (E and EA stay joined
        # through the other A step), so the move does not emit it.
        # The first also fails class 3 at zero substrates: its segment's steady state
        # depends on how fast the rapid-equilibrium steps run.
        single_A = @_testhelper_unguarded @enzyme_mechanism begin
            substrates: A, B
            products: P, Q
            steps: begin
                E + A <--> EA
                EB + A ⇌ EAB
                (E + B ⇌ EB, EA + B ⇌ EAB)
                (E + P ⇌ EP, EQ + P ⇌ EPQ)
                (E + Q ⇌ EQ, EP + Q ⇌ EPQ)
                EAB <--> EPQ
            end
        end
        v = EnzymeRates.Derivation._guard(single_A)
        @test (v.fault, v.limit, v.face) ===
              (EnzymeRates.Derivation._TRAP, EnzymeRates.Derivation._LIMIT_SPEED, 3)
        @test EnzymeRates.Mechanisms.step_text(single_A.catalytic, v.edge) ==
              "EA + B ⇌ EAB"
        @test _testhelper_respeed_verdicts(single_A; rng = Random.MersenneTwister(93),
                                           law = :limit)[3] === :pattern
        for single in (
            single_A,
            @enzyme_mechanism begin
                substrates: A, B
                products: P, Q
                steps: begin
                    E + A ⇌ EA
                    EB + A <--> EAB
                    (E + B ⇌ EB, EA + B ⇌ EAB)
                    (E + P ⇌ EP, EQ + P ⇌ EPQ)
                    (E + Q ⇌ EQ, EP + Q ⇌ EPQ)
                    EAB <--> EPQ
                end
            end)
            @test !(single in kids)
            @test _testhelper_segment_count(single) == _testhelper_segment_count(m)
        end
    end

    @testset "_expand_re_to_ss: ter-ter substrate cube emits one flip per group" begin
        # Random substrate addition, ordered product release. Each of the three
        # substrate groups carries all four cube edges for its metabolite, so
        # flipping one cuts every route between the forms it joins; each product
        # release is a single step on the ordered path. Six groups, six children.
        m = @enzyme_mechanism begin
            substrates: A, B, C
            products: P, Q, R
            steps: begin
                (E + A ⇌ EA, EB + A ⇌ EAB, EC + A ⇌ EAC,
                 EBC + A ⇌ EABC)
                (E + B ⇌ EB, EA + B ⇌ EAB, EC + B ⇌ EBC,
                 EAC + B ⇌ EABC)
                (E + C ⇌ EC, EA + C ⇌ EAC, EB + C ⇌ EBC,
                 EAB + C ⇌ EABC)
                EABC <--> EPQR
                EQR + P ⇌ EPQR
                ER + Q ⇌ EQR
                E + R ⇌ ER
            end
        end
        flipA = @enzyme_mechanism begin
            substrates: A, B, C
            products: P, Q, R
            steps: begin
                (E + A <--> EA, EB + A <--> EAB,
                 EC + A <--> EAC, EBC + A <--> EABC)
                (E + B ⇌ EB, EA + B ⇌ EAB, EC + B ⇌ EBC,
                 EAC + B ⇌ EABC)
                (E + C ⇌ EC, EA + C ⇌ EAC, EB + C ⇌ EBC,
                 EAB + C ⇌ EABC)
                EABC <--> EPQR
                EQR + P ⇌ EPQR
                ER + Q ⇌ EQR
                E + R ⇌ ER
            end
        end
        flipB = @enzyme_mechanism begin
            substrates: A, B, C
            products: P, Q, R
            steps: begin
                (E + A ⇌ EA, EB + A ⇌ EAB, EC + A ⇌ EAC,
                 EBC + A ⇌ EABC)
                (E + B <--> EB, EA + B <--> EAB,
                 EC + B <--> EBC, EAC + B <--> EABC)
                (E + C ⇌ EC, EA + C ⇌ EAC, EB + C ⇌ EBC,
                 EAB + C ⇌ EABC)
                EABC <--> EPQR
                EQR + P ⇌ EPQR
                ER + Q ⇌ EQR
                E + R ⇌ ER
            end
        end
        flipC = @enzyme_mechanism begin
            substrates: A, B, C
            products: P, Q, R
            steps: begin
                (E + A ⇌ EA, EB + A ⇌ EAB, EC + A ⇌ EAC,
                 EBC + A ⇌ EABC)
                (E + B ⇌ EB, EA + B ⇌ EAB, EC + B ⇌ EBC,
                 EAC + B ⇌ EABC)
                (E + C <--> EC, EA + C <--> EAC,
                 EB + C <--> EBC, EAB + C <--> EABC)
                EABC <--> EPQR
                EQR + P ⇌ EPQR
                ER + Q ⇌ EQR
                E + R ⇌ ER
            end
        end
        flipP = @enzyme_mechanism begin
            substrates: A, B, C
            products: P, Q, R
            steps: begin
                (E + A ⇌ EA, EB + A ⇌ EAB, EC + A ⇌ EAC,
                 EBC + A ⇌ EABC)
                (E + B ⇌ EB, EA + B ⇌ EAB, EC + B ⇌ EBC,
                 EAC + B ⇌ EABC)
                (E + C ⇌ EC, EA + C ⇌ EAC, EB + C ⇌ EBC,
                 EAB + C ⇌ EABC)
                EABC <--> EPQR
                EQR + P <--> EPQR
                ER + Q ⇌ EQR
                E + R ⇌ ER
            end
        end
        flipQ = @enzyme_mechanism begin
            substrates: A, B, C
            products: P, Q, R
            steps: begin
                (E + A ⇌ EA, EB + A ⇌ EAB, EC + A ⇌ EAC,
                 EBC + A ⇌ EABC)
                (E + B ⇌ EB, EA + B ⇌ EAB, EC + B ⇌ EBC,
                 EAC + B ⇌ EABC)
                (E + C ⇌ EC, EA + C ⇌ EAC, EB + C ⇌ EBC,
                 EAB + C ⇌ EABC)
                EABC <--> EPQR
                EQR + P ⇌ EPQR
                ER + Q <--> EQR
                E + R ⇌ ER
            end
        end
        flipR = @enzyme_mechanism begin
            substrates: A, B, C
            products: P, Q, R
            steps: begin
                (E + A ⇌ EA, EB + A ⇌ EAB, EC + A ⇌ EAC,
                 EBC + A ⇌ EABC)
                (E + B ⇌ EB, EA + B ⇌ EAB, EC + B ⇌ EBC,
                 EAC + B ⇌ EABC)
                (E + C ⇌ EC, EA + C ⇌ EAC, EB + C ⇌ EBC,
                 EAB + C ⇌ EABC)
                EABC <--> EPQR
                EQR + P ⇌ EPQR
                ER + Q ⇌ EQR
                E + R <--> ER
            end
        end
        kids = _testhelper_move(EnzymeRates.Enumeration._expand_re_to_ss, m)
        @test length(kids) == 6
        @test Set(kids) == Set([flipA, flipB, flipC, flipP, flipQ, flipR])
    end

    @testset "_expand_re_to_ss: split ter-ter pairs the A and B parts" begin
        # Ter-ter after one context split: A is split by whether B is bound and B
        # by whether A is bound. Each of the four parts alone leaves E and EA
        # joined through its sibling, so nothing flips one part by itself; C, P, Q
        # and R still flip alone. Of the six pairs, {A1, B1} cuts both binding edges
        # at E and EC and leaves {EA, EB, EAB, EAC, EBC, EABC}
        # as a segment with no substrate-free form, which the boundary check rejects.
        # {A1, B2} isolates {EA, EAC}: its four edges to the rest all have
        # weight 1 inward (A uptake) and −1 outward (B uptake against the +2 offset
        # of EAB), so every cycle through them is balanced and the pair carries
        # no flux; the turnover runs through EB inside the big segment. {A2, B1}
        # isolates {EB, EBC} the same way. Both count as failed and have no
        # extension, since every three-set containing them holds a gaining pair.
        # The other three pairs cut a block that also holds the chemistry edge
        # (weight 4 or 5 against flipped edges of weight 2 or 1; every cycle through
        # it sums to 6) and are emitted.
        m = @enzyme_mechanism begin
            substrates: A, B, C
            products: P, Q, R
            steps: begin
                (E + A ⇌ EA, EC + A ⇌ EAC)
                (EB + A ⇌ EAB, EBC + A ⇌ EABC)
                (E + B ⇌ EB, EC + B ⇌ EBC)
                (EA + B ⇌ EAB, EAC + B ⇌ EABC)
                (E + C ⇌ EC, EA + C ⇌ EAC, EB + C ⇌ EBC,
                 EAB + C ⇌ EABC)
                EABC <--> EPQR
                EQR + P ⇌ EPQR
                ER + Q ⇌ EQR
                E + R ⇌ ER
            end
        end
        # A group of graph `x`, edges `ks`, by the metabolite its first step binds (`:iso`
        # for none) and the forms (`_form_name`) its steps bind to.
        _testhelper_group_key(x, ks) =
            (e = x.catalytic.steps[first(ks)];
             (EnzymeRates.Enumeration._is_binding(e) ?
              x.catalytic.metabolites[only(e.taken_up)] : :iso,
              Set(x.catalytic.enzyme_forms[x.catalytic.steps[k].from]
                for k in ks)))
        # Each group is pinned by the metabolite it binds and the forms it binds to.
        A1, A2 = (:A, Set([:E, :EC])), (:A, Set([:EB, :EBC]))
        B1, B2 = (:B, Set([:E, :EC])), (:B, Set([:EA, :EAC]))
        C = (:C, Set([:E, :EA, :EB, :EAB]))
        P, Q, R = (:P, Set([:EQR])), (:Q, Set([:ER])), (:R, Set([:E]))
        _testhelper_group_of(key) =
            only(g for (g, ks) in enumerate(_testhelper_group_members(m.catalytic))
                 if _testhelper_group_key(m, ks) == key)
        flip(keys...) = _testhelper_flip_groups(m, [_testhelper_group_of(k) for k in keys])
        expected = [
            flip(C),        # {C}
            flip(P),        # {P}
            flip(Q),        # {Q}
            flip(R),        # {R}
            flip(A1, A2),   # {A1, A2}
            flip(B1, B2),   # {B1, B2}
            flip(A2, B2),   # {A2, B2}
        ]
        kids = _testhelper_move(EnzymeRates.Enumeration._expand_re_to_ss, m)
        @test length(kids) == 7
        @test Set(kids) == Set(expected)
        for c in kids
            ss = count((A1, A2, B1, B2)) do key
                ks = only(ks for ks in _testhelper_group_members(c.catalytic)
                          if _testhelper_group_key(c, ks) == key)
                !any(k -> c.catalytic.steps[k].rapid_equilibrium, ks)
            end
            @test ss != 1
        end
        # The two zero-flux pairs divide a segment yet keep the parent's rank;
        # flip(A2, B1) is flip(A1, B2) with A and B swapped.
        for absent in (flip(A1, B2), flip(A2, B1))
            @test !(absent in kids)
            @test _testhelper_segment_count(absent) > _testhelper_segment_count(m)
        end
        @test _testhelper_identifiable_rank(flip(A1, B2)) ==
            _testhelper_identifiable_rank(m)
    end

    @testset "_expand_re_to_ss: ping-pong" begin
        # The rapid-equilibrium subgraph is two trees, {E, EA, EQ} and
        # {E(P; res), E(; res), E(B; res)}, so every RE group is a bridge and
        # flipping any one of them alone raises the segment count. Both
        # isomerizations are steady state, and each of their ends has one other
        # step, a binding into that end alone in its group: EA → E(P; res) and
        # E(B; res) → EQ are qualifying chains whose flanks are all four binding
        # groups (`_chain_flank_groups`). No group is a unit, so no child. The A flip,
        # absent, keeps the parent's rank; the B, P and Q flips are its images under
        # swapping the two half-reactions and reversing the reaction.
        m = @enzyme_mechanism begin
            substrates: A, B
            products: P, Q
            steps: begin
                E + A ⇌ EA
                EA <--> EP_res_plusA_minusP
                EP_res_plusA_minusP ⇌ E_res_plusA_minusP + P
                E_res_plusA_minusP + B ⇌ EB_res_plusA_minusP
                EB_res_plusA_minusP <--> EQ
                EQ ⇌ E + Q
            end
        end
        flipA = @enzyme_mechanism begin
            substrates: A, B
            products: P, Q
            steps: begin
                E + A <--> EA
                EA <--> EP_res_plusA_minusP
                EP_res_plusA_minusP ⇌ E_res_plusA_minusP + P
                E_res_plusA_minusP + B ⇌ EB_res_plusA_minusP
                EB_res_plusA_minusP <--> EQ
                EQ ⇌ E + Q
            end
        end
        kids = _testhelper_move(EnzymeRates.Enumeration._expand_re_to_ss, m)
        @test isempty(kids)
        @test _testhelper_identifiable_rank(flipA) == _testhelper_identifiable_rank(m)
    end

    @testset "_expand_re_to_ss: an allosteric parent keeps a hyperbolic scheme" begin
        # Uni-bi with random product release. Flipping the S group alone splits
        # off {ES} and the equation stays degree 1. Flipping the P group splits
        # off {EP, EPQ}, whose term carries Q from EPQ beyond its
        # bottom EP and again from the edge EQ → EPQ, so the equation
        # carries Q²; flipping Q mirrors it with P². The move emits all three, for a
        # plain parent and for an allosteric one over two subunits. For the allosteric
        # parent candidate_mechanisms returns only the S flip (`_enumerator_returns`):
        # the guard keeps the P and Q flips, and the hyperbolic rule refuses them.
        plain = @enzyme_mechanism begin
            substrates: S
            products: P, Q
            steps: begin
                E + S ⇌ ES
                ES <--> EPQ
                (E + P ⇌ EP, EQ + P ⇌ EPQ)
                (E + Q ⇌ EQ, EP + Q ⇌ EPQ)
            end
        end
        flipS = @enzyme_mechanism begin
            substrates: S
            products: P, Q
            steps: begin
                E + S <--> ES
                ES <--> EPQ
                (E + P ⇌ EP, EQ + P ⇌ EPQ)
                (E + Q ⇌ EQ, EP + Q ⇌ EPQ)
            end
        end
        flipP = @enzyme_mechanism begin
            substrates: S
            products: P, Q
            steps: begin
                E + S ⇌ ES
                ES <--> EPQ
                (E + P <--> EP, EQ + P <--> EPQ)
                (E + Q ⇌ EQ, EP + Q ⇌ EPQ)
            end
        end
        flipQ = @enzyme_mechanism begin
            substrates: S
            products: P, Q
            steps: begin
                E + S ⇌ ES
                ES <--> EPQ
                (E + P ⇌ EP, EQ + P ⇌ EPQ)
                (E + Q <--> EQ, EP + Q <--> EPQ)
            end
        end
        kids = _testhelper_move(EnzymeRates.Enumeration._expand_re_to_ss, plain)
        @test length(kids) == 3
        @test Set(kids) == Set([flipS, flipP, flipQ])

        allo = @enzyme_mechanism begin
            substrates: S
            products: P, Q
            oligomeric_state: 2
            steps: begin
                E + S ⇌ ES                                :: NonequalAI
                ES <--> EPQ                           :: NonequalAI
                (E + P ⇌ EP, EQ + P ⇌ EPQ)          :: NonequalAI
                (E + Q ⇌ EQ, EP + Q ⇌ EPQ)          :: NonequalAI
            end
        end
        allo_flipS = @enzyme_mechanism begin
            substrates: S
            products: P, Q
            oligomeric_state: 2
            steps: begin
                E + S <--> ES                             :: NonequalAI
                ES <--> EPQ                           :: NonequalAI
                (E + P ⇌ EP, EQ + P ⇌ EPQ)          :: NonequalAI
                (E + Q ⇌ EQ, EP + Q ⇌ EPQ)          :: NonequalAI
            end
        end
        allo_flipP = @enzyme_mechanism begin
            substrates: S
            products: P, Q
            oligomeric_state: 2
            steps: begin
                E + S ⇌ ES                                :: NonequalAI
                ES <--> EPQ                           :: NonequalAI
                (E + P <--> EP, EQ + P <--> EPQ)    :: NonequalAI
                (E + Q ⇌ EQ, EP + Q ⇌ EPQ)          :: NonequalAI
            end
        end
        allo_flipQ = @enzyme_mechanism begin
            substrates: S
            products: P, Q
            oligomeric_state: 2
            steps: begin
                E + S ⇌ ES                                :: NonequalAI
                ES <--> EPQ                           :: NonequalAI
                (E + P ⇌ EP, EQ + P ⇌ EPQ)          :: NonequalAI
                (E + Q <--> EQ, EP + Q <--> EPQ)    :: NonequalAI
            end
        end
        rxn = @enzyme_reaction begin
            substrates: S[CN]
            products: P[C], Q[N]
            oligomeric_state: 2
        end
        allo_move = _testhelper_move(EnzymeRates.Enumeration._expand_re_to_ss, allo)
        @test length(allo_move) == 3
        @test Set(allo_move) == Set([allo_flipS, allo_flipP, allo_flipQ])
        # The flips are the children with the parent's allostery: every group is
        # :NonequalAI already, so no group relaxes, and a split adds a group.
        allo_kids = [c for c in _testhelper_children(rxn, [allo])
                     if c.allosteric == allo.allosteric]
        @test length(allo_kids) == 1
        @test Set(allo_kids) == Set([allo_flipS])
        for d in (allo_flipP, allo_flipQ)
            @test EnzymeRates.Enumeration._enumerator_keeps(d)
            @test !EnzymeRates.Enumeration._hyperbolic_where_required(
                _testhelper_candidate(d))
        end
    end

    @testset "_expand_re_to_ss: a one-subunit allosteric parent keeps every flip" begin
        # The same uni-bi as an allosteric mechanism with one catalytic subunit.
        # The conformational equilibrium adds no concentration power, so
        # candidate_mechanisms returns the P and Q flips it refuses for two subunits,
        # and the same three flips as for the plain mechanism.
        allo_one = @enzyme_mechanism begin
            substrates: S
            products: P, Q
            oligomeric_state: 1
            steps: begin
                E + S ⇌ ES                                :: NonequalAI
                ES <--> EPQ                           :: NonequalAI
                (E + P ⇌ EP, EQ + P ⇌ EPQ)          :: NonequalAI
                (E + Q ⇌ EQ, EP + Q ⇌ EPQ)          :: NonequalAI
            end
        end
        flipS_one = @enzyme_mechanism begin
            substrates: S
            products: P, Q
            oligomeric_state: 1
            steps: begin
                E + S <--> ES                             :: NonequalAI
                ES <--> EPQ                           :: NonequalAI
                (E + P ⇌ EP, EQ + P ⇌ EPQ)          :: NonequalAI
                (E + Q ⇌ EQ, EP + Q ⇌ EPQ)          :: NonequalAI
            end
        end
        flipP_one = @enzyme_mechanism begin
            substrates: S
            products: P, Q
            oligomeric_state: 1
            steps: begin
                E + S ⇌ ES                                :: NonequalAI
                ES <--> EPQ                           :: NonequalAI
                (E + P <--> EP, EQ + P <--> EPQ)    :: NonequalAI
                (E + Q ⇌ EQ, EP + Q ⇌ EPQ)          :: NonequalAI
            end
        end
        flipQ_one = @enzyme_mechanism begin
            substrates: S
            products: P, Q
            oligomeric_state: 1
            steps: begin
                E + S ⇌ ES                                :: NonequalAI
                ES <--> EPQ                           :: NonequalAI
                (E + P ⇌ EP, EQ + P ⇌ EPQ)          :: NonequalAI
                (E + Q <--> EQ, EP + Q <--> EPQ)    :: NonequalAI
            end
        end
        rxn = @enzyme_reaction begin
            substrates: S[CN]
            products: P[C], Q[N]
        end
        # The flips are the children with the parent's allostery, as above.
        kids_one = [c for c in _testhelper_children(rxn, [allo_one])
                    if c.allosteric == allo_one.allosteric]
        @test length(kids_one) == 3
        @test Set(kids_one) == Set([flipS_one, flipP_one, flipQ_one])
    end

    @testset "_expand_re_to_ss: a zero-flux pair is absent, its supersets reachable" begin
        # Ordered bi-bi with the abortive complex EAQ, after both the A group and
        # the Q group were split by context: six rapid-equilibrium binding groups,
        # each holding one step. Singles: EA + B and EQ + P are bridges of the RE
        # graph (EAB and EPQ are leaves there), so each would flip alone, but
        # they are the flanks of the qualifying chain EAB → EPQ and never flip
        # (`_chain_flank_groups`); the four bindings on the square E–EA–EAQ–EQ
        # are no bridges. Pairs of square edges all divide the segment.
        # {E + A, E + Q} isolates {E}, joined to the rest by
        # two parallel edges of weight 0 (A uptake +1 against EA's offset +1; Q
        # uptake −1 against EQ's offset −1): a balanced block, so neither flipped
        # group carries flux and the child would be its parent plus two phantoms.
        # {EA + Q, EQ + A} isolates {EAQ} the same way (weights −1 + 1 and
        # 1 − 1). The other four pairs each cut a block that also holds the chemistry
        # edge, whose weight 3 or 4 against a binding edge of weight 1 makes the block
        # unbalanced. Each of the two zero-flux pairs has no extension: every
        # three-set containing it also contains a gaining pair, so its supersets are
        # reached by a later flip of those children.
        m = @enzyme_mechanism begin
            substrates: A, B
            products: P, Q
            steps: begin
                E + A ⇌ EA
                EQ + A ⇌ EAQ
                E + Q ⇌ EQ
                EA + Q ⇌ EAQ
                EA + B ⇌ EAB
                EAB <--> EPQ
                EQ + P ⇌ EPQ
            end
        end
        expected = [
            (@enzyme_mechanism begin      # {E + A, EA + Q}
                substrates: A, B
                products: P, Q
                steps: begin
                    E + A <--> EA
                    EQ + A ⇌ EAQ
                    E + Q ⇌ EQ
                    EA + Q <--> EAQ
                    EA + B ⇌ EAB
                    EAB <--> EPQ
                    EQ + P ⇌ EPQ
                end
            end),
            (@enzyme_mechanism begin      # {E + A, EQ + A}
                substrates: A, B
                products: P, Q
                steps: begin
                    E + A <--> EA
                    EQ + A <--> EAQ
                    E + Q ⇌ EQ
                    EA + Q ⇌ EAQ
                    EA + B ⇌ EAB
                    EAB <--> EPQ
                    EQ + P ⇌ EPQ
                end
            end),
            (@enzyme_mechanism begin      # {E + Q, EA + Q}
                substrates: A, B
                products: P, Q
                steps: begin
                    E + A ⇌ EA
                    EQ + A ⇌ EAQ
                    E + Q <--> EQ
                    EA + Q <--> EAQ
                    EA + B ⇌ EAB
                    EAB <--> EPQ
                    EQ + P ⇌ EPQ
                end
            end),
            (@enzyme_mechanism begin      # {E + Q, EQ + A}
                substrates: A, B
                products: P, Q
                steps: begin
                    E + A ⇌ EA
                    EQ + A <--> EAQ
                    E + Q <--> EQ
                    EA + Q ⇌ EAQ
                    EA + B ⇌ EAB
                    EAB <--> EPQ
                    EQ + P ⇌ EPQ
                end
            end),
        ]
        kids = _testhelper_move(EnzymeRates.Enumeration._expand_re_to_ss, m)
        @test length(kids) == 4
        @test Set(kids) == Set(expected)
        # The flank flip of B and the two zero-flux pairs are absent; each divides a
        # segment, and each has exactly its parent's rank. The flank flip of P is the
        # flip of B with the reaction reversed.
        r0 = _testhelper_identifiable_rank(m)
        for absent in (
            (@enzyme_mechanism begin  # {EA + B}
                substrates: A, B
                products: P, Q
                steps: begin
                    E + A ⇌ EA
                    EQ + A ⇌ EAQ
                    E + Q ⇌ EQ
                    EA + Q ⇌ EAQ
                    EA + B <--> EAB
                    EAB <--> EPQ
                    EQ + P ⇌ EPQ
                end
            end),
            (@enzyme_mechanism begin  # {E + A, E + Q}
                substrates: A, B
                products: P, Q
                steps: begin
                    E + A <--> EA
                    EQ + A ⇌ EAQ
                    E + Q <--> EQ
                    EA + Q ⇌ EAQ
                    EA + B ⇌ EAB
                    EAB <--> EPQ
                    EQ + P ⇌ EPQ
                end
            end),
            @enzyme_mechanism begin  # {EA + Q, EQ + A}
                substrates: A, B
                products: P, Q
                steps: begin
                    E + A ⇌ EA
                    EQ + A <--> EAQ
                    E + Q ⇌ EQ
                    EA + Q <--> EAQ
                    EA + B ⇌ EAB
                    EAB <--> EPQ
                    EQ + P ⇌ EPQ
                end
            end)
            @test !(absent in kids)
            @test _testhelper_segment_count(absent) > _testhelper_segment_count(m)
            @test _testhelper_identifiable_rank(absent) == r0
        end
    end

    @testset "_expand_re_to_ss: a Theorell–Chance parent flips each binding" begin
        # No isomerization step. Every step lies on the unbalanced cycle
        # E → EA → EQ → E (weights 1, 2, 1), so both rapid-equilibrium bindings
        # are units; each divides the one segment alone.
        m = @enzyme_mechanism begin
            substrates: A, B
            products: P, Q
            steps: begin
                E + A ⇌ EA
                EA + B <--> EQ + P
                E + Q ⇌ EQ
            end
        end
        expected = [
            (@enzyme_mechanism begin
                substrates: A, B
                products: P, Q
                steps: begin
                    E + A <--> EA
                    EA + B <--> EQ + P
                    E + Q ⇌ EQ
                end
            end),
            (@enzyme_mechanism begin
                substrates: A, B
                products: P, Q
                steps: begin
                    E + A ⇌ EA
                    EA + B <--> EQ + P
                    E + Q <--> EQ
                end
            end),
        ]
        kids = _testhelper_move(EnzymeRates.Enumeration._expand_re_to_ss, m)
        @test length(kids) == 2
        @test Set(kids) == Set(expected)
    end
end

@testset "_expand_split_kinetic_group (context bipartitions)" begin
    @testset "_expand_split_kinetic_group: random-order bi-bi frees indep params" begin
        # The split closure of this seed, which has 5 independent parameters, holds 16
        # structures and reaches the form with every step in its own group and 9
        # independent parameters.
        m = @enzyme_mechanism begin
            substrates: A, B
            products: P, Q
            steps: begin
                (E + A ⇌ EA, EB + A ⇌ EAB, EP + A ⇌ EAP)
                (E + B ⇌ EB, EA + B ⇌ EAB, EQ + B ⇌ EBQ)
                (E + P ⇌ EP, EQ + P ⇌ EPQ, EA + P ⇌ EAP)
                (E + Q ⇌ EQ, EP + Q ⇌ EPQ, EB + Q ⇌ EBQ)
                EAB <--> EPQ
            end
        end
        @test EnzymeRates.Derivation.fitted_parameter_count(m) == 5
        cl = _testhelper_closure(
            m, p -> _testhelper_move(EnzymeRates.Enumeration._expand_split_kinetic_group,
                                     p))
        @test maximum(length(_testhelper_group_members(m.catalytic)) for m in cl) == 13
        @test maximum(EnzymeRates.Derivation.fitted_parameter_count(m) for m in cl) == 9
        @test length(cl) == 16
    end

    @testset "_expand_split_kinetic_group: bi-bi without dead ends splits a square" begin
        # A single context split (A by B alone) is tied straight back by the A/B
        # square, so no child splits one group by itself. Each child frees one
        # interaction factor by separating both groups of one square.
        m = @enzyme_mechanism begin
            substrates: A, B
            products: P, Q
            steps: begin
                (E + A ⇌ EA, EB + A ⇌ EAB)
                (E + B ⇌ EB, EA + B ⇌ EAB)
                (E + P ⇌ EP, EQ + P ⇌ EPQ)
                (E + Q ⇌ EQ, EP + Q ⇌ EPQ)
                EAB <--> EPQ
            end
        end
        splitAB = @enzyme_mechanism begin
            substrates: A, B
            products: P, Q
            steps: begin
                E + A ⇌ EA
                EB + A ⇌ EAB
                E + B ⇌ EB
                EA + B ⇌ EAB
                (E + P ⇌ EP, EQ + P ⇌ EPQ)
                (E + Q ⇌ EQ, EP + Q ⇌ EPQ)
                EAB <--> EPQ
            end
        end
        splitPQ = @enzyme_mechanism begin
            substrates: A, B
            products: P, Q
            steps: begin
                (E + A ⇌ EA, EB + A ⇌ EAB)
                (E + B ⇌ EB, EA + B ⇌ EAB)
                E + P ⇌ EP
                EQ + P ⇌ EPQ
                E + Q ⇌ EQ
                EP + Q ⇌ EPQ
                EAB <--> EPQ
            end
        end
        @test EnzymeRates.Derivation.fitted_parameter_count(m) == 5
        kids = _testhelper_move(EnzymeRates.Enumeration._expand_split_kinetic_group, m)
        @test length(kids) == 2
        @test Set(kids) == Set([splitAB, splitPQ])
        for c in kids
            @test EnzymeRates.Derivation.fitted_parameter_count(c) == 6
        end
    end

    @testset "_expand_split_kinetic_group: bi-bi with dead ends splits each square" begin
        # Four four-cycles run through this graph — the A/B and P/Q catalytic
        # squares and the A/P and B/Q dead-end squares — and each gives one
        # child that frees the interaction factor of its doubly-bound form.
        m = @enzyme_mechanism begin
            substrates: A, B
            products: P, Q
            steps: begin
                (E + A ⇌ EA, EB + A ⇌ EAB, EP + A ⇌ EAP)
                (E + B ⇌ EB, EA + B ⇌ EAB, EQ + B ⇌ EBQ)
                (E + P ⇌ EP, EQ + P ⇌ EPQ, EA + P ⇌ EAP)
                (E + Q ⇌ EQ, EP + Q ⇌ EPQ, EB + Q ⇌ EBQ)
                EAB <--> EPQ
            end
        end
        splitAP = @enzyme_mechanism begin
            substrates: A, B
            products: P, Q
            steps: begin
                (E + A ⇌ EA, EB + A ⇌ EAB)
                EP + A ⇌ EAP
                (E + B ⇌ EB, EA + B ⇌ EAB, EQ + B ⇌ EBQ)
                (E + P ⇌ EP, EQ + P ⇌ EPQ)
                EA + P ⇌ EAP
                (E + Q ⇌ EQ, EP + Q ⇌ EPQ, EB + Q ⇌ EBQ)
                EAB <--> EPQ
            end
        end
        splitAB = @enzyme_mechanism begin
            substrates: A, B
            products: P, Q
            steps: begin
                (E + A ⇌ EA, EP + A ⇌ EAP)
                EB + A ⇌ EAB
                (E + B ⇌ EB, EQ + B ⇌ EBQ)
                EA + B ⇌ EAB
                (E + P ⇌ EP, EQ + P ⇌ EPQ, EA + P ⇌ EAP)
                (E + Q ⇌ EQ, EP + Q ⇌ EPQ, EB + Q ⇌ EBQ)
                EAB <--> EPQ
            end
        end
        splitBQ = @enzyme_mechanism begin
            substrates: A, B
            products: P, Q
            steps: begin
                (E + A ⇌ EA, EB + A ⇌ EAB, EP + A ⇌ EAP)
                (E + B ⇌ EB, EA + B ⇌ EAB)
                EQ + B ⇌ EBQ
                (E + P ⇌ EP, EQ + P ⇌ EPQ, EA + P ⇌ EAP)
                (E + Q ⇌ EQ, EP + Q ⇌ EPQ)
                EB + Q ⇌ EBQ
                EAB <--> EPQ
            end
        end
        splitPQ = @enzyme_mechanism begin
            substrates: A, B
            products: P, Q
            steps: begin
                (E + A ⇌ EA, EB + A ⇌ EAB, EP + A ⇌ EAP)
                (E + B ⇌ EB, EA + B ⇌ EAB, EQ + B ⇌ EBQ)
                (E + P ⇌ EP, EA + P ⇌ EAP)
                EQ + P ⇌ EPQ
                (E + Q ⇌ EQ, EB + Q ⇌ EBQ)
                EP + Q ⇌ EPQ
                EAB <--> EPQ
            end
        end
        @test EnzymeRates.Derivation.fitted_parameter_count(m) == 5
        kids = _testhelper_move(EnzymeRates.Enumeration._expand_split_kinetic_group, m)
        @test length(kids) == 4
        @test Set(kids) == Set([splitAP, splitAB, splitBQ, splitPQ])
        for c in kids
            @test EnzymeRates.Derivation.fitted_parameter_count(c) == 6
        end
    end

    @testset "_expand_split_kinetic_group: ter-ter cube splits a substrate pair" begin
        # Each child says "the affinity for one substrate depends on whether the
        # other is bound", a 2 + 2 division of two four-step groups that no carve
        # of a single step could produce. One child per substrate pair.
        m = @enzyme_mechanism begin
            substrates: A, B, C
            products: P, Q, R
            steps: begin
                (E + A ⇌ EA, EB + A ⇌ EAB, EC + A ⇌ EAC,
                 EBC + A ⇌ EABC)
                (E + B ⇌ EB, EA + B ⇌ EAB, EC + B ⇌ EBC,
                 EAC + B ⇌ EABC)
                (E + C ⇌ EC, EA + C ⇌ EAC, EB + C ⇌ EBC,
                 EAB + C ⇌ EABC)
                EABC <--> EPQR
                EQR + P ⇌ EPQR
                ER + Q ⇌ EQR
                E + R ⇌ ER
            end
        end
        splitAB = @enzyme_mechanism begin
            substrates: A, B, C
            products: P, Q, R
            steps: begin
                (E + A ⇌ EA, EC + A ⇌ EAC)
                (EB + A ⇌ EAB, EBC + A ⇌ EABC)
                (E + B ⇌ EB, EC + B ⇌ EBC)
                (EA + B ⇌ EAB, EAC + B ⇌ EABC)
                (E + C ⇌ EC, EA + C ⇌ EAC, EB + C ⇌ EBC,
                 EAB + C ⇌ EABC)
                EABC <--> EPQR
                EQR + P ⇌ EPQR
                ER + Q ⇌ EQR
                E + R ⇌ ER
            end
        end
        splitAC = @enzyme_mechanism begin
            substrates: A, B, C
            products: P, Q, R
            steps: begin
                (E + A ⇌ EA, EB + A ⇌ EAB)
                (EC + A ⇌ EAC, EBC + A ⇌ EABC)
                (E + B ⇌ EB, EA + B ⇌ EAB, EC + B ⇌ EBC,
                 EAC + B ⇌ EABC)
                (E + C ⇌ EC, EB + C ⇌ EBC)
                (EA + C ⇌ EAC, EAB + C ⇌ EABC)
                EABC <--> EPQR
                EQR + P ⇌ EPQR
                ER + Q ⇌ EQR
                E + R ⇌ ER
            end
        end
        splitBC = @enzyme_mechanism begin
            substrates: A, B, C
            products: P, Q, R
            steps: begin
                (E + A ⇌ EA, EB + A ⇌ EAB, EC + A ⇌ EAC,
                 EBC + A ⇌ EABC)
                (E + B ⇌ EB, EA + B ⇌ EAB)
                (EC + B ⇌ EBC, EAC + B ⇌ EABC)
                (E + C ⇌ EC, EA + C ⇌ EAC)
                (EB + C ⇌ EBC, EAB + C ⇌ EABC)
                EABC <--> EPQR
                EQR + P ⇌ EPQR
                ER + Q ⇌ EQR
                E + R ⇌ ER
            end
        end
        @test EnzymeRates.Derivation.fitted_parameter_count(m) == 7
        kids = _testhelper_move(EnzymeRates.Enumeration._expand_split_kinetic_group, m)
        @test length(kids) == 3
        @test Set(kids) == Set([splitAB, splitAC, splitBC])
        for c in kids
            @test EnzymeRates.Derivation.fitted_parameter_count(c) == 8
        end
    end

    @testset "_expand_split_kinetic_group: bi-bi seeds emit 102 children" begin
        # Aggregate regression pin over the whole bi-bi init set: 102 children from the
        # 55 seeds and 420 from their 184 merged and Theorell–Chance variants, which are
        # 176 from the 108 two-group merged, 236 from the 56 three-group merged and 8
        # from the 20 Theorell–Chance ones.
        init = _testhelper_mechanisms(EnzymeRates.Enumeration.init_mechanisms(bi_bi_rxn))
        n_children(ms) =
            sum(length(_testhelper_move(
                    EnzymeRates.Enumeration._expand_split_kinetic_group, m)) for m in ms)
        @test n_children(filter(_testhelper_holds_iso, init)) == 102
        @test n_children(filter(!_testhelper_holds_iso, init)) == 420
    end

    @testset "_expand_split_kinetic_group: a rejected split is the parent's model" begin
        # Level-1 candidates the count test rejects have the parent's identifiable
        # rank; the rendered equation may differ only in which tied name survives. The
        # B, P and Q groups' candidates are images of the A group's under relabeling
        # (A and B swapped with P and Q swapped) and reversing the reaction, so the
        # rank is checked on the A group's.
        m = @enzyme_mechanism begin
            substrates: A, B
            products: P, Q
            steps: begin
                (E + A ⇌ EA, EB + A ⇌ EAB, EP + A ⇌ EAP)
                (E + B ⇌ EB, EA + B ⇌ EAB, EQ + B ⇌ EBQ)
                (E + P ⇌ EP, EQ + P ⇌ EPQ, EA + P ⇌ EAP)
                (E + Q ⇌ EQ, EP + Q ⇌ EPQ, EB + Q ⇌ EBQ)
                EAB <--> EPQ
            end
        end
        base = EnzymeRates.Derivation.fitted_parameter_count(m)
        r0 = _testhelper_identifiable_rank(m)
        @test r0 == base
        rejected = 0
        g = m.catalytic.steps[_testhelper_edge(m, "E + A ⇌ EA")].kinetic_group
        for bp in EnzymeRates.Enumeration._context_bipartitions(_testhelper_candidate(m), g)
            child = EnzymeRates.Enumeration._apply_bipartitions(_testhelper_candidate(m),
                                                                [(g, bp)]).mechanism
            EnzymeRates.Derivation.fitted_parameter_count(child) > base && continue
            rejected += 1
            @test _testhelper_identifiable_rank(child) == r0
        end
        @test rejected > 0
    end

    @testset "_expand_split_kinetic_group: a four-step group splits two and two" begin
        m = @enzyme_mechanism begin
            substrates: A, B
            products: P, Q
            steps: begin
                (E + A ⇌ EA, EB + A ⇌ EAB, EP + A ⇌ EAP,
                 EBP + A ⇌ EABP)
                (E + B ⇌ EB, EA + B ⇌ EAB, EP + B ⇌ EBP,
                 EAP + B ⇌ EABP)
                (E + P ⇌ EP, EQ + P ⇌ EPQ, EA + P ⇌ EAP,
                 EB + P ⇌ EBP, EAB + P ⇌ EABP)
                (E + Q ⇌ EQ, EP + Q ⇌ EPQ)
                EAB <--> EPQ
            end
        end
        gi = m.catalytic.steps[_testhelper_edge(m, "E + A ⇌ EA")].kinetic_group
        a_group = findall(e -> e.kinetic_group == gi, m.catalytic.steps)
        @test length(a_group) == 4
        bps = EnzymeRates.Enumeration._context_bipartitions(_testhelper_candidate(m), gi)
        # Context B and context P each divide the four A steps two and two, a
        # division no carve of a single step can produce.
        @test count(bp -> length(bp[1]) == 2 && length(bp[2]) == 2, bps) == 2
        even_bp = first(bp for bp in bps if length(bp[1]) == 2)
        child = EnzymeRates.Enumeration._apply_bipartitions(_testhelper_candidate(m),
                                                            [(gi, even_bp)]).mechanism
        @test Set(findall(e -> e.kinetic_group == child.catalytic.steps[k].kinetic_group,
                          child.catalytic.steps)
                  for k in a_group) ==
              Set(even_bp)
    end
end

@testset "_expand_split_kinetic_group: a split member of a tied group turns alone" begin
    # E <--> D and ES <--> DS share their constants. The reaction turns over through
    # them in opposite directions, E → D and DS → ES, so the group turns as one with
    # its lead E → D, nearer free E, and stores ES → DS as written relative to it.
    # Split by S, ES <--> DS is a group of its own and runs the way the reaction
    # turns over through it, from DS. The child goes through the checked constructor,
    # which turns it; a child built from the parent's stored directions would keep
    # ES → DS.
    m = @enzyme_mechanism begin
        substrates: S
        products: P
        steps: begin
            E + S <--> ES
            ES <--> EP
            EP <--> E + P
            D + S <--> DS
            (D <--> E, DS <--> ES)
        end
    end
    @test m.catalytic.steps[_testhelper_edge(m, "ES <--> DS")].kinetic_group ==
          m.catalytic.steps[_testhelper_edge(m, "E <--> D")].kinetic_group
    kids = _testhelper_move(EnzymeRates.Enumeration._expand_split_kinetic_group, m)
    expected = [@enzyme_mechanism begin
        substrates: S
        products: P
        steps: begin
            E + S <--> ES
            ES <--> EP
            EP <--> E + P
            D + S <--> DS
            D <--> E
            DS <--> ES
        end
    end]
    @test length(kids) == 1
    @test Set(kids) == Set(expected)
    child = only(kids)
    @test _testhelper_edge(child, "DS <--> ES") isa Int
    # Written back as stored, the child turns no step.
    @test _testhelper_scheme(_testhelper_written(child)).steps == child.catalytic.steps
end

@testset "_expand_split_kinetic_group: by conformation" begin
    # S binds both conformations in one kinetic group; EstarS is a dead end
    # whose dissociation constant no cycle ties, so dividing the group by
    # conformation frees one parameter. The other groups are singletons.
    m = @enzyme_mechanism begin
        substrates: S
        products: P
        steps: begin
            (E + S ⇌ ES, Estar + S ⇌ EstarS)
            ES <--> EstarP
            EstarP ⇌ Estar + P
            Estar <--> E
        end
    end
    by_conformation = @enzyme_mechanism begin
        substrates: S
        products: P
        steps: begin
            E + S ⇌ ES
            Estar + S ⇌ EstarS
            ES <--> EstarP
            EstarP ⇌ Estar + P
            Estar <--> E
        end
    end
    kids = _testhelper_move(EnzymeRates.Enumeration._expand_split_kinetic_group, m)
    @test length(kids) == 1
    @test Set(kids) == Set([by_conformation])
end

@testset "_expand_split_kinetic_group: by residual" begin
    # Q binds free E on the cycle and the two modified forms as dead ends. The
    # Q group divides by B ({E, E(; res)} | {E(B; res)}) and by residual
    # ({E} | {E(; res), E(B; res)}); each division frees a dead-end constant
    # on its own, and only one bipartition per group is ever applied.
    m = @enzyme_mechanism begin
        substrates: A, B
        products: P, Q
        steps: begin
            (E + Q ⇌ EQ, E_res_plusA_minusP + Q ⇌ EQ_res_plusA_minusP,
             EB_res_plusA_minusP + Q ⇌ EBQ_res_plusA_minusP)
            E + A ⇌ EA
            EA <--> EP_res_plusA_minusP
            EP_res_plusA_minusP ⇌ E_res_plusA_minusP + P
            E_res_plusA_minusP + B ⇌ EB_res_plusA_minusP
            EB_res_plusA_minusP <--> EQ
        end
    end
    by_B = @enzyme_mechanism begin
        substrates: A, B
        products: P, Q
        steps: begin
            (E + Q ⇌ EQ, E_res_plusA_minusP + Q ⇌ EQ_res_plusA_minusP)
            EB_res_plusA_minusP + Q ⇌ EBQ_res_plusA_minusP
            E + A ⇌ EA
            EA <--> EP_res_plusA_minusP
            EP_res_plusA_minusP ⇌ E_res_plusA_minusP + P
            E_res_plusA_minusP + B ⇌ EB_res_plusA_minusP
            EB_res_plusA_minusP <--> EQ
        end
    end
    by_residual = @enzyme_mechanism begin
        substrates: A, B
        products: P, Q
        steps: begin
            E + Q ⇌ EQ
            (E_res_plusA_minusP + Q ⇌ EQ_res_plusA_minusP,
             EB_res_plusA_minusP + Q ⇌ EBQ_res_plusA_minusP)
            E + A ⇌ EA
            EA <--> EP_res_plusA_minusP
            EP_res_plusA_minusP ⇌ E_res_plusA_minusP + P
            E_res_plusA_minusP + B ⇌ EB_res_plusA_minusP
            EB_res_plusA_minusP <--> EQ
        end
    end
    kids = _testhelper_move(EnzymeRates.Enumeration._expand_split_kinetic_group, m)
    @test length(kids) == 2
    @test Set(kids) == Set([by_B, by_residual])
end

@testset "expansion moves: ter-ter random-order seed within budget" begin
    # The seed of `init_mechanisms(terter)` with the most steps (55) is the enumeration's
    # worst case: random order, each substrate and product binding at nine forms in one
    # rapid-equilibrium group. Outside garbage collection, the split alone must take under
    # 90 s and the four moves that apply to a plain mechanism (split, flip, dead end,
    # to-allosteric) under 120 s together; they return 151 children: 12 splits, 6 flips,
    # the 2^7 − 1 = 127 nonempty sets of the seven groups tagged :OnlyA, and six pairs
    # tagged :NonequalAI (the reaction declares no regulator). Tagged :NonequalAI the
    # chemistry alone and each pair of bindings gain a constant beyond L; the pairs of
    # substrates and the pairs of products are refused tagged :OnlyA. The split takes 4 s
    # in a full local suite and 7-11 s on a CI runner, but a slow macOS runner reached
    # 60 s.
    terter = @enzyme_reaction begin
        substrates: A[C], B[N], C[O]
        products: P[C], Q[N], R[O]
    end
    seed = @enzyme_mechanism begin
        substrates: A, B, C
        products: P, Q, R
        steps: begin
            (E + A ⇌ EA, EB + A ⇌ EAB, EBC + A ⇌ EABC,
             EBP + A ⇌ EABP, EC + A ⇌ EAC, ECQ + A ⇌ EACQ,
             EP + A ⇌ EAP, EPQ + A ⇌ EAPQ, EQ + A ⇌ EAQ)
            (E + B ⇌ EB, EA + B ⇌ EAB, EAC + B ⇌ EABC,
             EAP + B ⇌ EABP, EC + B ⇌ EBC, ECR + B ⇌ EBCR,
             EP + B ⇌ EBP, EPR + B ⇌ EBPR, ER + B ⇌ EBR)
            (E + C ⇌ EC, EA + C ⇌ EAC, EAB + C ⇌ EABC,
             EAQ + C ⇌ EACQ, EB + C ⇌ EBC, EBR + C ⇌ EBCR,
             EQ + C ⇌ ECQ, EQR + C ⇌ ECQR, ER + C ⇌ ECR)
            (E + P ⇌ EP, EA + P ⇌ EAP, EAB + P ⇌ EABP,
             EAQ + P ⇌ EAPQ, EB + P ⇌ EBP, EBR + P ⇌ EBPR,
             EQ + P ⇌ EPQ, EQR + P ⇌ EPQR, ER + P ⇌ EPR)
            (E + Q ⇌ EQ, EA + Q ⇌ EAQ, EAC + Q ⇌ EACQ,
             EAP + Q ⇌ EAPQ, EC + Q ⇌ ECQ, ECR + Q ⇌ ECQR,
             EP + Q ⇌ EPQ, EPR + Q ⇌ EPQR, ER + Q ⇌ EQR)
            (E + R ⇌ ER, EB + R ⇌ EBR, EBC + R ⇌ EBCR,
             EBP + R ⇌ EBPR, EC + R ⇌ ECR, ECQ + R ⇌ ECQR,
             EP + R ⇌ EPR, EPQ + R ⇌ EPQR, EQ + R ⇌ EQR)
            EABC <--> EPQR
        end
    end
    worst = seed
    @test length(worst.catalytic.steps) == 55
    split = @timed _testhelper_move(EnzymeRates.Enumeration._expand_split_kinetic_group,
                                    worst)
    rest = @timed vcat(
        _testhelper_move(EnzymeRates.Enumeration._expand_re_to_ss, worst),
        _testhelper_move(EnzymeRates.Enumeration._expand_add_dead_end_regulator, worst,
                         terter),
        _testhelper_move(EnzymeRates.Enumeration._expand_to_allosteric, worst, terter))
    t_split = split.time - split.gctime
    t_rest = rest.time - rest.gctime
    @info "ter-ter seed's moves: split $(round(split.time; digits = 1)) s, " *
          "$(round(split.gctime; digits = 1)) s of it in garbage collection; " *
          "the other moves $(round(rest.time; digits = 1)) s, " *
          "$(round(rest.gctime; digits = 1)) s of it in garbage collection"
    @test length(split.value) == 12
    @test t_split < 90
    @test length(rest.value) == 6 + 127 + 6
    @test t_split + t_rest < 120
end

@testset "catalytic moves on bi-bi to depth 2: counts and both rules on every child" begin
    # Aggregate regression pin over the whole enumerated population of the bi-bi
    # reaction whose atoms admit ping-pong, every reactant also a competitive
    # inhibitor: `init_mechanisms` plus two levels of the flip, split and dead-end
    # moves, deduplicated across levels. Every mechanism satisfies both emission rules.
    rxn = @enzyme_reaction begin
        substrates: A[CX], B[N]
        products: P[C], Q[NX]
        competitive_inhibitors: A, B, P, Q
    end
    _testhelper_moves(m, rxn) =
        vcat(_testhelper_move(EnzymeRates.Enumeration._expand_re_to_ss, m),
             _testhelper_move(EnzymeRates.Enumeration._expand_split_kinetic_group, m),
             _testhelper_move(EnzymeRates.Enumeration._expand_add_dead_end_regulator, m,
                              rxn))
    iso_kept = true
    function _testhelper_levels(rxn)
        level = unique!(
            _testhelper_mechanisms(EnzymeRates.Enumeration.init_mechanisms(rxn)))
        seen = Set(level)
        out = [level]
        for _ in 1:2
            next = EnzymeRates.Mechanism[]
            for m in level, c in _testhelper_moves(m, rxn)
                iso_kept &= _testhelper_holds_iso(c) == _testhelper_holds_iso(m)
                c in seen && continue
                push!(seen, c); push!(next, c)
            end
            level = next
            push!(out, level)
        end
        out
    end
    obeys_rules(m) = EnzymeRates.Enumeration._assert_emission_rules(
        _testhelper_candidate(m)) === nothing
    free(m) = isempty(
        EnzymeRates.Enumeration._bound_comp_inhibitors(_testhelper_candidate(m)))
    copies = _testhelper_levels(rxn)
    @test iso_kept
    @test all(obeys_rules, Iterators.flatten(copies))
    @test any(!free, copies[2])

    # The counts that follow are of the levels grown from the 62 seeds that hold an
    # isomerization; the moves neither add nor remove one. Without the copy rule
    # level 1 would hold 1,749: the 120 seed-level placements whose every complex has
    # a productive twin and whose dwell gauge is consistent are not emitted, and the
    # 240 whose gauge fails, the shared-group family of Case 3, are. 8,694 level-2
    # mechanisms hold a copy group whose every complex has a productive twin and
    # whose gauge fails. Without the flip rule the levels would hold 62, 1,649 and
    # 31,730: it leaves out 20 seed children at level 1 and 348 mechanisms at level 2,
    # each holding a qualifying chain whose isomerization and at least one flank are
    # steady state. The 202 merged and Theorell–Chance variants, which hold none, grow
    # apart from the seeds: 5,109 mechanisms at level 1 and 97,143 at level 2. Grown
    # from all 264 init mechanisms the levels hold the sums, 264, 6,738 and 128,525.
    @test [count(_testhelper_holds_iso, l) for l in copies] == [62, 1629, 31382]
    @test length.(copies) == [264, 6738, 128525]

    # The moves never remove a step, so a mechanism that binds no inhibitor descends
    # only from inhibitor-free ones and sits at the level it has in the same reaction
    # declaring no inhibitor: the plain bi-bi reaction. Without the flux and
    # new-complex rules its levels would hold 62, 369 and 1,388; the flux rule leaves
    # out 148 zero-flux flip children at level 2, and the 40 zero-flux split children,
    # emitted at rapid equilibrium, duplicate level-2 flip children, so 1,200 would
    # remain; it affects no seed or level-1 child. The flip rule
    # (`_chain_flank_groups`) then leaves out 20 seed children at level 1 and 66
    # mechanisms at level 2, each holding a qualifying chain whose isomerization and
    # at least one flank are steady state, and each with its twin whose flanks are
    # both at rapid equilibrium in the population: 62, 349 and 1,134 remain from the
    # seeds that hold an isomerization. The variants hold none: 669 mechanisms at
    # level 1 and 1,269 at level 2. Class 3 inside the flip move's guard drops
    # 16 level-2 variant flips (`_SPLIT`) and adds 48 superset flips, all binding no
    # inhibitor. Grown from all 264 init mechanisms the levels hold the sums, 264, 1,018
    # and 2,403.
    @test [count(m -> free(m) && _testhelper_holds_iso(m), l) for l in copies] ==
          [62, 349, 1134]
    @test [count(m -> free(m) && !_testhelper_holds_iso(m), l) for l in copies] ==
          [202, 669, 1269]
    @test [count(free, l) for l in copies] == [264, 1018, 2403]
end

@testset "every bi-bi init mechanism and child is what its nodes and steps build" begin
    # An aggregate regression pin over the bi-bi seed set and its first expansion level
    # (the exception to writing every fixture inline). Every mechanism the moves emit
    # goes through the checked constructor, so writing it back as its nodes and steps
    # (`_testhelper_written`) and building that again gives the same mechanism: no child
    # keeps a direction, an order or a group number the constructor would change.
    init = _testhelper_mechanisms(EnzymeRates.Enumeration.init_mechanisms(bi_bi_rxn))
    kids = _testhelper_children(bi_bi_rxn, init)
    @test length(init) == 239
    @test length(kids) == 4779
    @test count(k -> _testhelper_is_allosteric(k), kids) == 3781
    @test all(c -> _testhelper_build(_testhelper_written(c)) == c, [init; kids])
end

@testset "candidate_mechanisms on a merged uni-uni" begin
    # E + S → EP fused (SS), E + P ⇌ EP (RE). Flip: the P group cuts E from EP,
    # carries flux, no bottomless segment: one child, both steps SS. Split: single-step
    # groups, nothing. Dead end: no regulator. To-allosteric at multiplicity 1: the
    # chemistry alone :OnlyA and both groups :OnlyA; P's binding alone is refused
    # (`_FAST_STEP_PINNED`), and no set tagged :NonequalAI is refused as :OnlyA.
    rxn = @enzyme_reaction begin
        substrates: S[C]
        products: P[C]
    end
    m = @enzyme_mechanism begin
        substrates: S
        products: P
        steps: begin
            E + S <--> EP
            E + P ⇌ EP
        end
    end
    flipped = @enzyme_mechanism begin
        substrates: S
        products: P
        steps: begin
            E + S <--> EP
            E + P <--> EP
        end
    end
    chemistry_only_a = @enzyme_mechanism begin
        substrates: S
        products: P
        steps: begin
            E + S <--> EP    :: OnlyA
            E + P ⇌ EP       :: EqualAI
        end
    end
    both_only_a = @enzyme_mechanism begin
        substrates: S
        products: P
        steps: begin
            E + S <--> EP    :: OnlyA
            E + P ⇌ EP       :: OnlyA
        end
    end
    kids = _testhelper_children(rxn, [m])
    @test length(kids) == 3
    @test Set(kids) == Set([flipped, chemistry_only_a, both_only_a])
    # The parent's rate is k·(S − P/Keq)/(1 + P/Kp), the Haldane relation fixing the
    # fused step's reverse constant: 2 fitted, rank 2. The flipped child's two forms at
    # steady state give V·(S − P/Keq)/(1 + a·S + b·P) with V, a and b free: 3 fitted,
    # rank 3. With both groups :OnlyA the inactive conformation binds nothing and runs no
    # chemistry, so L enters only as 1 + L beside the free enzyme and folds into k and
    # Kp: 3 fitted, rank 2.
    @test _testhelper_fitted(m) == _testhelper_identifiable_rank(m) == 2
    @test _testhelper_fitted(flipped) == _testhelper_identifiable_rank(flipped) == 3
    @test _testhelper_fitted(both_only_a) == 3 &&
          _testhelper_identifiable_rank(both_only_a) == 2
    # That child, whose inactive conformation is free E alone, is an L phantom: its
    # rescaling certificate maps it onto the law of its catalytic mechanism.
    @test EnzymeRates.Derivation._is_phantom(EnzymeRates.Derivation._memo_law(both_only_a))
    @test _testhelper_phantom_map_holds(both_only_a; rng = Random.MersenneTwister(1))
end

@testset "_expand_change_allo_state on a merged uni-uni with both groups :OnlyA" begin
    # The fused step E + S → EP is the chemistry group. Relaxing P to :NonequalAI
    # keeps the chemistry :OnlyA: the inactive conformation binds P with a constant of
    # its own and runs no chemistry, a valid child. The chemistry tagged :NonequalAI
    # alone is refused: restored, it would join E to EP in the inactive conformation
    # beside the :OnlyA P binding, which then cannot be switched off
    # (`_FAST_STEP_PINNED`). P's set gains, so the search builds no larger one. One child.
    both_only_a = @enzyme_mechanism begin
        substrates: S
        products: P
        oligomeric_state: 1
        steps: begin
            E + S <--> EP    :: OnlyA
            E + P ⇌ EP       :: OnlyA
        end
    end
    p_relaxed = @enzyme_mechanism begin
        substrates: S
        products: P
        oligomeric_state: 1
        steps: begin
            E + S <--> EP    :: OnlyA
            E + P ⇌ EP       :: NonequalAI
        end
    end
    kids = _testhelper_move(EnzymeRates.Enumeration._expand_change_allo_state, both_only_a)
    @test length(kids) == 1
    @test Set(kids) == Set([p_relaxed])
    # The child's rate is k·(S − P/Keq)/(1 + L + P·(1/Kp_A + L/Kp_I)): dividing by
    # 1 + L leaves k/(1 + L) and one P coefficient, so k, Kp_A, Kp_I and L give
    # 4 fitted, rank 2.
    @test _testhelper_fitted(p_relaxed) == 4 &&
          _testhelper_identifiable_rank(p_relaxed) == 2
end

@testset "_expand_to_allosteric on an ordered Theorell–Chance mechanism" begin
    # EA + B → EQ + P takes up B and gives off P in one step: the chemistry group. At one
    # subunit with no regulator the move tags each nonempty set of the three groups
    # :OnlyA: seven children, none :NonequalAI, since each minimal set that gains is
    # returned tagged :OnlyA. With the chemistry :OnlyA the inactive conformation keeps
    # the A and Q bindings it is not tagged on as the tree EA – E – EQ, so every
    # binding points away from E and admission refuses none.
    m = @enzyme_mechanism begin
        substrates: A, B
        products: P, Q
        steps: begin
            E + A ⇌ EA
            EA + B <--> EQ + P
            E + Q ⇌ EQ
        end
    end
    rxn = @enzyme_reaction begin
        substrates: A[C], B[N]
        products: P[C], Q[N]
    end
    a = @_testhelper_unguarded @enzyme_mechanism begin
        substrates: A, B ; products: P, Q ; oligomeric_state: 1
        steps: begin
            E + A ⇌ EA         :: OnlyA
            EA + B <--> EQ + P :: EqualAI
            E + Q ⇌ EQ         :: EqualAI
        end
    end
    chem = @_testhelper_unguarded @enzyme_mechanism begin
        substrates: A, B ; products: P, Q ; oligomeric_state: 1
        steps: begin
            E + A ⇌ EA         :: EqualAI
            EA + B <--> EQ + P :: OnlyA
            E + Q ⇌ EQ         :: EqualAI
        end
    end
    a_chem = @_testhelper_unguarded @enzyme_mechanism begin
        substrates: A, B ; products: P, Q ; oligomeric_state: 1
        steps: begin
            E + A ⇌ EA         :: OnlyA
            EA + B <--> EQ + P :: OnlyA
            E + Q ⇌ EQ         :: EqualAI
        end
    end
    q = @_testhelper_unguarded @enzyme_mechanism begin
        substrates: A, B ; products: P, Q ; oligomeric_state: 1
        steps: begin
            E + A ⇌ EA         :: EqualAI
            EA + B <--> EQ + P :: EqualAI
            E + Q ⇌ EQ         :: OnlyA
        end
    end
    a_q = @_testhelper_unguarded @enzyme_mechanism begin
        substrates: A, B ; products: P, Q ; oligomeric_state: 1
        steps: begin
            E + A ⇌ EA         :: OnlyA
            EA + B <--> EQ + P :: EqualAI
            E + Q ⇌ EQ         :: OnlyA
        end
    end
    chem_q = @_testhelper_unguarded @enzyme_mechanism begin
        substrates: A, B ; products: P, Q ; oligomeric_state: 1
        steps: begin
            E + A ⇌ EA         :: EqualAI
            EA + B <--> EQ + P :: OnlyA
            E + Q ⇌ EQ         :: OnlyA
        end
    end
    a_chem_q = @_testhelper_unguarded @enzyme_mechanism begin
        substrates: A, B ; products: P, Q ; oligomeric_state: 1
        steps: begin
            E + A ⇌ EA         :: OnlyA
            EA + B <--> EQ + P :: OnlyA
            E + Q ⇌ EQ         :: OnlyA
        end
    end
    kids = _testhelper_move(EnzymeRates.Enumeration._expand_to_allosteric, m, rxn)
    @test length(kids) == 7
    @test Set(kids) == Set([a, chem, a_chem, q, a_q, chem_q, a_chem_q])
    # The parent's rate is (k/Ka)·(A·B − P·Q/Keq)/(1 + A/Ka + Q/Kq): 3 fitted, rank 3.
    # Each child with the chemistry :OnlyA has a dead inactive conformation that adds
    # (1 + L) times the bindings it keeps to the bindings it loses, so dividing by 1 + L
    # rescales the lost bindings' constants and the turnover and L disappears: 4 fitted,
    # rank 3. The rank is checked on a_chem and a_chem_q; chem_q is a_chem with the
    # reaction reversed.
    @test _testhelper_fitted(m) == _testhelper_identifiable_rank(m) == 3
    for k in (a_chem, chem_q, a_chem_q)
        @test _testhelper_fitted(k) == 4
    end
    for k in (a_chem, a_chem_q)
        @test _testhelper_identifiable_rank(k) == 3
    end
    # All three are L phantoms, each mapped onto the parent's law by its certificate,
    # a_chem and chem_q though their inactive conformation keeps a binding at free E.
    for k in (a_chem, chem_q, a_chem_q)
        @test EnzymeRates.Derivation._is_phantom(EnzymeRates.Derivation._memo_law(k))
        @test _testhelper_phantom_map_holds(k; rng = Random.MersenneTwister(2))
    end
end

@testset "_expand_to_allosteric on a merged ordered bi-bi at multiplicity 2" begin
    # The ordered bi-bi seed with EAB merged into EPQ: B's binding is fused with the
    # chemistry (EA + B ⇌ EPQ), and the A and P bindings are steady state. The scheme is
    # hyperbolic, so multiplicity 2 stays open: the RE segments are {E, EQ}, where EQ
    # carries one Q, and {EA, EPQ}, where EPQ carries one B; the steady-state A and P
    # bindings join them, and each metabolite scores on one edge only, an edge that leaves
    # the segment carrying that metabolite (`_hyperbolic_catalysis`). No regulator: the
    # move tags each nonempty set of the four groups :OnlyA, 15 children, and A's
    # steady-state binding alone :NonequalAI, a minimal gaining set whose :OnlyA version
    # the enumerator does not return: 16 children. The guard refuses A's and B's groups
    # :OnlyA together: EA is cut off in the inactive conformation, and at zero A and B the
    # fused binding, vanishing one way into EA, fills it faster than A's steady-state
    # binding empties it (`_FILLED`).
    m2 = @enzyme_mechanism begin
        substrates: A, B
        products: P, Q
        steps: begin
            E + A <--> EA
            EA + B ⇌ EPQ
            EQ + P <--> EPQ
            E + Q ⇌ EQ
        end
    end
    rxn2 = @enzyme_reaction begin
        substrates: A[C], B[N]
        products: P[C], Q[N]
        oligomeric_state: 2
    end
    a = @_testhelper_unguarded @enzyme_mechanism begin
        substrates: A, B ; products: P, Q ; oligomeric_state: 2
        steps: begin
            E + A <--> EA   :: OnlyA
            EA + B ⇌ EPQ    :: EqualAI
            EQ + P <--> EPQ :: EqualAI
            E + Q ⇌ EQ      :: EqualAI
        end
    end
    b = @_testhelper_unguarded @enzyme_mechanism begin
        substrates: A, B ; products: P, Q ; oligomeric_state: 2
        steps: begin
            E + A <--> EA   :: EqualAI
            EA + B ⇌ EPQ    :: OnlyA
            EQ + P <--> EPQ :: EqualAI
            E + Q ⇌ EQ      :: EqualAI
        end
    end
    a_b = @_testhelper_unguarded @enzyme_mechanism begin
        substrates: A, B ; products: P, Q ; oligomeric_state: 2
        steps: begin
            E + A <--> EA   :: OnlyA
            EA + B ⇌ EPQ    :: OnlyA
            EQ + P <--> EPQ :: EqualAI
            E + Q ⇌ EQ      :: EqualAI
        end
    end
    p = @_testhelper_unguarded @enzyme_mechanism begin
        substrates: A, B ; products: P, Q ; oligomeric_state: 2
        steps: begin
            E + A <--> EA   :: EqualAI
            EA + B ⇌ EPQ    :: EqualAI
            EQ + P <--> EPQ :: OnlyA
            E + Q ⇌ EQ      :: EqualAI
        end
    end
    a_p = @_testhelper_unguarded @enzyme_mechanism begin
        substrates: A, B ; products: P, Q ; oligomeric_state: 2
        steps: begin
            E + A <--> EA   :: OnlyA
            EA + B ⇌ EPQ    :: EqualAI
            EQ + P <--> EPQ :: OnlyA
            E + Q ⇌ EQ      :: EqualAI
        end
    end
    b_p = @_testhelper_unguarded @enzyme_mechanism begin
        substrates: A, B ; products: P, Q ; oligomeric_state: 2
        steps: begin
            E + A <--> EA   :: EqualAI
            EA + B ⇌ EPQ    :: OnlyA
            EQ + P <--> EPQ :: OnlyA
            E + Q ⇌ EQ      :: EqualAI
        end
    end
    a_b_p = @_testhelper_unguarded @enzyme_mechanism begin
        substrates: A, B ; products: P, Q ; oligomeric_state: 2
        steps: begin
            E + A <--> EA   :: OnlyA
            EA + B ⇌ EPQ    :: OnlyA
            EQ + P <--> EPQ :: OnlyA
            E + Q ⇌ EQ      :: EqualAI
        end
    end
    q = @_testhelper_unguarded @enzyme_mechanism begin
        substrates: A, B ; products: P, Q ; oligomeric_state: 2
        steps: begin
            E + A <--> EA   :: EqualAI
            EA + B ⇌ EPQ    :: EqualAI
            EQ + P <--> EPQ :: EqualAI
            E + Q ⇌ EQ      :: OnlyA
        end
    end
    a_q = @_testhelper_unguarded @enzyme_mechanism begin
        substrates: A, B ; products: P, Q ; oligomeric_state: 2
        steps: begin
            E + A <--> EA   :: OnlyA
            EA + B ⇌ EPQ    :: EqualAI
            EQ + P <--> EPQ :: EqualAI
            E + Q ⇌ EQ      :: OnlyA
        end
    end
    b_q = @_testhelper_unguarded @enzyme_mechanism begin
        substrates: A, B ; products: P, Q ; oligomeric_state: 2
        steps: begin
            E + A <--> EA   :: EqualAI
            EA + B ⇌ EPQ    :: OnlyA
            EQ + P <--> EPQ :: EqualAI
            E + Q ⇌ EQ      :: OnlyA
        end
    end
    a_b_q = @_testhelper_unguarded @enzyme_mechanism begin
        substrates: A, B ; products: P, Q ; oligomeric_state: 2
        steps: begin
            E + A <--> EA   :: OnlyA
            EA + B ⇌ EPQ    :: OnlyA
            EQ + P <--> EPQ :: EqualAI
            E + Q ⇌ EQ      :: OnlyA
        end
    end
    p_q = @_testhelper_unguarded @enzyme_mechanism begin
        substrates: A, B ; products: P, Q ; oligomeric_state: 2
        steps: begin
            E + A <--> EA   :: EqualAI
            EA + B ⇌ EPQ    :: EqualAI
            EQ + P <--> EPQ :: OnlyA
            E + Q ⇌ EQ      :: OnlyA
        end
    end
    a_p_q = @_testhelper_unguarded @enzyme_mechanism begin
        substrates: A, B ; products: P, Q ; oligomeric_state: 2
        steps: begin
            E + A <--> EA   :: OnlyA
            EA + B ⇌ EPQ    :: EqualAI
            EQ + P <--> EPQ :: OnlyA
            E + Q ⇌ EQ      :: OnlyA
        end
    end
    b_p_q = @_testhelper_unguarded @enzyme_mechanism begin
        substrates: A, B ; products: P, Q ; oligomeric_state: 2
        steps: begin
            E + A <--> EA   :: EqualAI
            EA + B ⇌ EPQ    :: OnlyA
            EQ + P <--> EPQ :: OnlyA
            E + Q ⇌ EQ      :: OnlyA
        end
    end
    a_b_p_q = @_testhelper_unguarded @enzyme_mechanism begin
        substrates: A, B ; products: P, Q ; oligomeric_state: 2
        steps: begin
            E + A <--> EA   :: OnlyA
            EA + B ⇌ EPQ    :: OnlyA
            EQ + P <--> EPQ :: OnlyA
            E + Q ⇌ EQ      :: OnlyA
        end
    end
    a_nonequal = @enzyme_mechanism begin
        substrates: A, B ; products: P, Q ; oligomeric_state: 2
        steps: begin
            E + A <--> EA   :: NonequalAI
            EA + B ⇌ EPQ    :: EqualAI
            EQ + P <--> EPQ :: EqualAI
            E + Q ⇌ EQ      :: EqualAI
        end
    end
    @test EnzymeRates.Enumeration._hyperbolic_catalysis(_testhelper_candidate(m2))
    kids = _testhelper_move(EnzymeRates.Enumeration._expand_to_allosteric, m2, rxn2)
    @test all(c -> _testhelper_is_allosteric(c), kids)
    @test length(kids) == 16
    @test Set(kids) == Set([a, b, a_b, p, a_p, b_p, a_b_p, q, a_q, b_q, a_b_q, p_q, a_p_q,
                            b_p_q, a_b_p_q, a_nonequal])
    @test !EnzymeRates.Enumeration._enumerator_returns(_testhelper_candidate(a))
    v = EnzymeRates.Derivation._guard(a_b)
    @test (v.fault, EnzymeRates.Mechanisms.step_text(a_b.catalytic, v.edge)) ==
          (EnzymeRates.Derivation._FILLED, "EA + B ⇌ EPQ")
    # The parent has 5 fitted constants, all identifiable. Each child with B's group
    # :OnlyA that the guard admits adds L; at multiplicity 2 each conformation's binding
    # polynomial enters squared, so dividing by 1 + L no longer folds L into the other
    # constants: 6 fitted, rank 6.
    @test _testhelper_fitted(m2) == _testhelper_identifiable_rank(m2) == 5
    for k in (b_p, b_q, a_b_p, a_b_q, b_p_q, a_b_p_q)
        @test _testhelper_fitted(k) == _testhelper_identifiable_rank(k) == 6
    end
    # At two subunits no child is an L phantom.
    @test !any(k -> EnzymeRates.Derivation._is_phantom(EnzymeRates.Derivation._memo_law(k)),
               (a_b, b_p, b_q, a_b_p, a_b_q, b_p_q, a_b_p_q))
end

@testset "_expand_to_allosteric on a merged bi-bi whose B group mixes fused and plain" begin
    # The merged ordered bi-bi with B's abortive binding at EQ in B's own group:
    # EA + B ⇌ EPQ is fused, EQ + B ⇌ EBQ is plain, one RE kinetic group. A group
    # holding a chemistry step is chemistry, so the whole B group, its plain step
    # included, is the one set the move pairs with the declared regulator R. The move
    # emits the 15 nonempty sets of the four groups :OnlyA, B's group alone :OnlyA again
    # with R :OnlyA and :OnlyI, and A's steady-state binding alone :NonequalAI, a minimal
    # gaining set whose :OnlyA version the enumerator does not return: 18 children. The
    # guard refuses A's and B's groups :OnlyA together: EA is cut off in the inactive
    # conformation, and at zero A and B the fused binding, vanishing one way into EA,
    # fills it faster than A's steady-state binding empties it (`_FILLED`).
    m = @enzyme_mechanism begin
        substrates: A, B
        products: P, Q
        steps: begin
            E + A <--> EA
            (EA + B ⇌ EPQ, EQ + B ⇌ EBQ)
            EQ + P <--> EPQ
            E + Q ⇌ EQ
        end
    end
    rxn = @enzyme_reaction begin
        substrates: A[C], B[N]
        products: P[C], Q[N]
        allosteric_regulators: R
    end
    a = @_testhelper_unguarded @enzyme_mechanism begin
        substrates: A, B ; products: P, Q ; oligomeric_state: 1
        steps: begin
            E + A <--> EA                :: OnlyA
            (EA + B ⇌ EPQ, EQ + B ⇌ EBQ) :: EqualAI
            EQ + P <--> EPQ              :: EqualAI
            E + Q ⇌ EQ                   :: EqualAI
        end
    end
    b = @_testhelper_unguarded @enzyme_mechanism begin
        substrates: A, B ; products: P, Q ; oligomeric_state: 1
        steps: begin
            E + A <--> EA                :: EqualAI
            (EA + B ⇌ EPQ, EQ + B ⇌ EBQ) :: OnlyA
            EQ + P <--> EPQ              :: EqualAI
            E + Q ⇌ EQ                   :: EqualAI
        end
    end
    a_b = @_testhelper_unguarded @enzyme_mechanism begin
        substrates: A, B ; products: P, Q ; oligomeric_state: 1
        steps: begin
            E + A <--> EA                :: OnlyA
            (EA + B ⇌ EPQ, EQ + B ⇌ EBQ) :: OnlyA
            EQ + P <--> EPQ              :: EqualAI
            E + Q ⇌ EQ                   :: EqualAI
        end
    end
    p = @_testhelper_unguarded @enzyme_mechanism begin
        substrates: A, B ; products: P, Q ; oligomeric_state: 1
        steps: begin
            E + A <--> EA                :: EqualAI
            (EA + B ⇌ EPQ, EQ + B ⇌ EBQ) :: EqualAI
            EQ + P <--> EPQ              :: OnlyA
            E + Q ⇌ EQ                   :: EqualAI
        end
    end
    a_p = @_testhelper_unguarded @enzyme_mechanism begin
        substrates: A, B ; products: P, Q ; oligomeric_state: 1
        steps: begin
            E + A <--> EA                :: OnlyA
            (EA + B ⇌ EPQ, EQ + B ⇌ EBQ) :: EqualAI
            EQ + P <--> EPQ              :: OnlyA
            E + Q ⇌ EQ                   :: EqualAI
        end
    end
    b_p = @_testhelper_unguarded @enzyme_mechanism begin
        substrates: A, B ; products: P, Q ; oligomeric_state: 1
        steps: begin
            E + A <--> EA                :: EqualAI
            (EA + B ⇌ EPQ, EQ + B ⇌ EBQ) :: OnlyA
            EQ + P <--> EPQ              :: OnlyA
            E + Q ⇌ EQ                   :: EqualAI
        end
    end
    a_b_p = @_testhelper_unguarded @enzyme_mechanism begin
        substrates: A, B ; products: P, Q ; oligomeric_state: 1
        steps: begin
            E + A <--> EA                :: OnlyA
            (EA + B ⇌ EPQ, EQ + B ⇌ EBQ) :: OnlyA
            EQ + P <--> EPQ              :: OnlyA
            E + Q ⇌ EQ                   :: EqualAI
        end
    end
    q = @_testhelper_unguarded @enzyme_mechanism begin
        substrates: A, B ; products: P, Q ; oligomeric_state: 1
        steps: begin
            E + A <--> EA                :: EqualAI
            (EA + B ⇌ EPQ, EQ + B ⇌ EBQ) :: EqualAI
            EQ + P <--> EPQ              :: EqualAI
            E + Q ⇌ EQ                   :: OnlyA
        end
    end
    a_q = @_testhelper_unguarded @enzyme_mechanism begin
        substrates: A, B ; products: P, Q ; oligomeric_state: 1
        steps: begin
            E + A <--> EA                :: OnlyA
            (EA + B ⇌ EPQ, EQ + B ⇌ EBQ) :: EqualAI
            EQ + P <--> EPQ              :: EqualAI
            E + Q ⇌ EQ                   :: OnlyA
        end
    end
    b_q = @_testhelper_unguarded @enzyme_mechanism begin
        substrates: A, B ; products: P, Q ; oligomeric_state: 1
        steps: begin
            E + A <--> EA                :: EqualAI
            (EA + B ⇌ EPQ, EQ + B ⇌ EBQ) :: OnlyA
            EQ + P <--> EPQ              :: EqualAI
            E + Q ⇌ EQ                   :: OnlyA
        end
    end
    a_b_q = @_testhelper_unguarded @enzyme_mechanism begin
        substrates: A, B ; products: P, Q ; oligomeric_state: 1
        steps: begin
            E + A <--> EA                :: OnlyA
            (EA + B ⇌ EPQ, EQ + B ⇌ EBQ) :: OnlyA
            EQ + P <--> EPQ              :: EqualAI
            E + Q ⇌ EQ                   :: OnlyA
        end
    end
    p_q = @_testhelper_unguarded @enzyme_mechanism begin
        substrates: A, B ; products: P, Q ; oligomeric_state: 1
        steps: begin
            E + A <--> EA                :: EqualAI
            (EA + B ⇌ EPQ, EQ + B ⇌ EBQ) :: EqualAI
            EQ + P <--> EPQ              :: OnlyA
            E + Q ⇌ EQ                   :: OnlyA
        end
    end
    a_p_q = @_testhelper_unguarded @enzyme_mechanism begin
        substrates: A, B ; products: P, Q ; oligomeric_state: 1
        steps: begin
            E + A <--> EA                :: OnlyA
            (EA + B ⇌ EPQ, EQ + B ⇌ EBQ) :: EqualAI
            EQ + P <--> EPQ              :: OnlyA
            E + Q ⇌ EQ                   :: OnlyA
        end
    end
    b_p_q = @_testhelper_unguarded @enzyme_mechanism begin
        substrates: A, B ; products: P, Q ; oligomeric_state: 1
        steps: begin
            E + A <--> EA                :: EqualAI
            (EA + B ⇌ EPQ, EQ + B ⇌ EBQ) :: OnlyA
            EQ + P <--> EPQ              :: OnlyA
            E + Q ⇌ EQ                   :: OnlyA
        end
    end
    a_b_p_q = @_testhelper_unguarded @enzyme_mechanism begin
        substrates: A, B ; products: P, Q ; oligomeric_state: 1
        steps: begin
            E + A <--> EA                :: OnlyA
            (EA + B ⇌ EPQ, EQ + B ⇌ EBQ) :: OnlyA
            EQ + P <--> EPQ              :: OnlyA
            E + Q ⇌ EQ                   :: OnlyA
        end
    end
    b_r_only_a = @_testhelper_unguarded @enzyme_mechanism begin
        substrates: A, B ; products: P, Q ; oligomeric_state: 1
        steps: begin
            E + A <--> EA                :: EqualAI
            (EA + B ⇌ EPQ, EQ + B ⇌ EBQ) :: OnlyA
            EQ + P <--> EPQ              :: EqualAI
            E + Q ⇌ EQ                   :: EqualAI
        end
        regulatory_site(1): R::OnlyA
    end
    b_r_only_i = @_testhelper_unguarded @enzyme_mechanism begin
        substrates: A, B ; products: P, Q ; oligomeric_state: 1
        steps: begin
            E + A <--> EA                :: EqualAI
            (EA + B ⇌ EPQ, EQ + B ⇌ EBQ) :: OnlyA
            EQ + P <--> EPQ              :: EqualAI
            E + Q ⇌ EQ                   :: EqualAI
        end
        regulatory_site(1): R::OnlyI
    end
    a_nonequal = @enzyme_mechanism begin
        substrates: A, B ; products: P, Q ; oligomeric_state: 1
        steps: begin
            E + A <--> EA                :: NonequalAI
            (EA + B ⇌ EPQ, EQ + B ⇌ EBQ) :: EqualAI
            EQ + P <--> EPQ              :: EqualAI
            E + Q ⇌ EQ                   :: EqualAI
        end
    end
    kids = _testhelper_move(EnzymeRates.Enumeration._expand_to_allosteric, m, rxn)
    @test length(kids) == 18
    @test Set(kids) == Set([a, b, a_b, p, a_p, b_p, a_b_p, q, a_q, b_q, a_b_q, p_q, a_p_q,
                            b_p_q, a_b_p_q, b_r_only_a, b_r_only_i, a_nonequal])
    @test !EnzymeRates.Enumeration._enumerator_returns(_testhelper_candidate(a))
    v = EnzymeRates.Derivation._guard(a_b)
    @test (v.fault, EnzymeRates.Mechanisms.step_text(a_b.catalytic, v.edge)) ==
          (EnzymeRates.Derivation._FILLED, "EA + B ⇌ EPQ")
    # The parent has 5 fitted constants, all identifiable. Each child with B's group
    # :OnlyA that the guard admits adds L, and its dead inactive conformation keeps the
    # bindings it reaches from the free enzyme through :EqualAI groups: A's where A is
    # :EqualAI; Q's where Q is, then P's at EQ where P is too. Where it keeps one, L
    # multiplies a binding polynomial that carries A or Q, which no rescaling of the
    # parent's constants absorbs: 6 fitted, rank 6. Where A and Q are both :OnlyA it binds
    # nothing, so L multiplies the free enzyme's terms alone, and 1 + L folds into A's
    # binding rate and Q's dissociation constant: 6 fitted, rank 5.
    @test _testhelper_fitted(m) == _testhelper_identifiable_rank(m) == 5
    for k in (b_p, b_q, a_b_p, b_p_q)
        @test _testhelper_fitted(k) == _testhelper_identifiable_rank(k) == 6
    end
    for k in (a_b_q, a_b_p_q)
        @test _testhelper_fitted(k) == 6 && _testhelper_identifiable_rank(k) == 5
    end
    # The two whose inactive conformation binds nothing are L phantoms, each mapped onto
    # the parent's law by its certificate; the five that keep a binding there have no
    # certificate and are fitted.
    for k in (a_b_q, a_b_p_q)
        @test EnzymeRates.Derivation._is_phantom(EnzymeRates.Derivation._memo_law(k))
        @test _testhelper_phantom_map_holds(k; rng = Random.MersenneTwister(3))
    end
    for k in (a_b, b_p, b_q, a_b_p, b_p_q)
        @test !EnzymeRates.Derivation._is_phantom(EnzymeRates.Derivation._memo_law(k))
    end
end

@testset "_expand_to_allosteric on a random-order merged bi-bi led by plain bindings" begin
    # The random-order merged bi-bi: A and B bind free E in either order, and the second
    # binding is fused with the chemistry. A's steady-state group holds E + A → EA and
    # the fused EB + A → EPQ; B's rapid-equilibrium group holds E + B ⇌ EB and
    # the fused EA + B ⇌ EPQ. Canonical order sorts a group's steps by source form,
    # E before EA and EB, so each group opens with its plain binding. A group holding
    # a chemistry step is chemistry whatever its first step, so the sets the move pairs
    # with the declared regulator R are A's group, B's and both. The reaction allows one
    # catalytic subunit. The move emits the 15 nonempty sets of the four groups :OnlyA,
    # the three chemistry sets again with R :OnlyA and :OnlyI, and A's steady-state group
    # alone :NonequalAI, a minimal gaining set whose :OnlyA version the enumerator does
    # not return: 22 children. Read by its first step, neither group would be chemistry,
    # and no set would be paired with R.
    m = @enzyme_mechanism begin
        substrates: A, B
        products: P, Q
        steps: begin
            (E + A <--> EA, EB + A <--> EPQ)
            (E + B ⇌ EB, EA + B ⇌ EPQ)
            E + Q ⇌ EQ
            EQ + P <--> EPQ
        end
    end
    rxn = @enzyme_reaction begin
        substrates: A[C], B[N]
        products: P[C], Q[N]
        allosteric_regulators: R
    end
    c = _testhelper_candidate(m)
    _testhelper_plain(k) = (e = m.catalytic.steps[k];
                            EnzymeRates.Enumeration._is_binding(e) &&
                                !_testhelper_is_fused(c, e))
    chemistry = [ks for ks in _testhelper_group_members(m.catalytic)
                 if !all(_testhelper_plain, ks)]
    @test length(chemistry) == 2
    @test all(ks -> _testhelper_plain(first(ks)), chemistry)
    a = @_testhelper_unguarded @enzyme_mechanism begin
        substrates: A, B ; products: P, Q ; oligomeric_state: 1
        steps: begin
            (E + A <--> EA, EB + A <--> EPQ) :: OnlyA
            (E + B ⇌ EB, EA + B ⇌ EPQ)       :: EqualAI
            E + Q ⇌ EQ                       :: EqualAI
            EQ + P <--> EPQ                  :: EqualAI
        end
    end
    b = @_testhelper_unguarded @enzyme_mechanism begin
        substrates: A, B ; products: P, Q ; oligomeric_state: 1
        steps: begin
            (E + A <--> EA, EB + A <--> EPQ) :: EqualAI
            (E + B ⇌ EB, EA + B ⇌ EPQ)       :: OnlyA
            E + Q ⇌ EQ                       :: EqualAI
            EQ + P <--> EPQ                  :: EqualAI
        end
    end
    a_b = @_testhelper_unguarded @enzyme_mechanism begin
        substrates: A, B ; products: P, Q ; oligomeric_state: 1
        steps: begin
            (E + A <--> EA, EB + A <--> EPQ) :: OnlyA
            (E + B ⇌ EB, EA + B ⇌ EPQ)       :: OnlyA
            E + Q ⇌ EQ                       :: EqualAI
            EQ + P <--> EPQ                  :: EqualAI
        end
    end
    q = @_testhelper_unguarded @enzyme_mechanism begin
        substrates: A, B ; products: P, Q ; oligomeric_state: 1
        steps: begin
            (E + A <--> EA, EB + A <--> EPQ) :: EqualAI
            (E + B ⇌ EB, EA + B ⇌ EPQ)       :: EqualAI
            E + Q ⇌ EQ                       :: OnlyA
            EQ + P <--> EPQ                  :: EqualAI
        end
    end
    a_q = @_testhelper_unguarded @enzyme_mechanism begin
        substrates: A, B ; products: P, Q ; oligomeric_state: 1
        steps: begin
            (E + A <--> EA, EB + A <--> EPQ) :: OnlyA
            (E + B ⇌ EB, EA + B ⇌ EPQ)       :: EqualAI
            E + Q ⇌ EQ                       :: OnlyA
            EQ + P <--> EPQ                  :: EqualAI
        end
    end
    b_q = @_testhelper_unguarded @enzyme_mechanism begin
        substrates: A, B ; products: P, Q ; oligomeric_state: 1
        steps: begin
            (E + A <--> EA, EB + A <--> EPQ) :: EqualAI
            (E + B ⇌ EB, EA + B ⇌ EPQ)       :: OnlyA
            E + Q ⇌ EQ                       :: OnlyA
            EQ + P <--> EPQ                  :: EqualAI
        end
    end
    a_b_q = @_testhelper_unguarded @enzyme_mechanism begin
        substrates: A, B ; products: P, Q ; oligomeric_state: 1
        steps: begin
            (E + A <--> EA, EB + A <--> EPQ) :: OnlyA
            (E + B ⇌ EB, EA + B ⇌ EPQ)       :: OnlyA
            E + Q ⇌ EQ                       :: OnlyA
            EQ + P <--> EPQ                  :: EqualAI
        end
    end
    p = @_testhelper_unguarded @enzyme_mechanism begin
        substrates: A, B ; products: P, Q ; oligomeric_state: 1
        steps: begin
            (E + A <--> EA, EB + A <--> EPQ) :: EqualAI
            (E + B ⇌ EB, EA + B ⇌ EPQ)       :: EqualAI
            E + Q ⇌ EQ                       :: EqualAI
            EQ + P <--> EPQ                  :: OnlyA
        end
    end
    a_p = @_testhelper_unguarded @enzyme_mechanism begin
        substrates: A, B ; products: P, Q ; oligomeric_state: 1
        steps: begin
            (E + A <--> EA, EB + A <--> EPQ) :: OnlyA
            (E + B ⇌ EB, EA + B ⇌ EPQ)       :: EqualAI
            E + Q ⇌ EQ                       :: EqualAI
            EQ + P <--> EPQ                  :: OnlyA
        end
    end
    b_p = @_testhelper_unguarded @enzyme_mechanism begin
        substrates: A, B ; products: P, Q ; oligomeric_state: 1
        steps: begin
            (E + A <--> EA, EB + A <--> EPQ) :: EqualAI
            (E + B ⇌ EB, EA + B ⇌ EPQ)       :: OnlyA
            E + Q ⇌ EQ                       :: EqualAI
            EQ + P <--> EPQ                  :: OnlyA
        end
    end
    a_b_p = @_testhelper_unguarded @enzyme_mechanism begin
        substrates: A, B ; products: P, Q ; oligomeric_state: 1
        steps: begin
            (E + A <--> EA, EB + A <--> EPQ) :: OnlyA
            (E + B ⇌ EB, EA + B ⇌ EPQ)       :: OnlyA
            E + Q ⇌ EQ                       :: EqualAI
            EQ + P <--> EPQ                  :: OnlyA
        end
    end
    q_p = @_testhelper_unguarded @enzyme_mechanism begin
        substrates: A, B ; products: P, Q ; oligomeric_state: 1
        steps: begin
            (E + A <--> EA, EB + A <--> EPQ) :: EqualAI
            (E + B ⇌ EB, EA + B ⇌ EPQ)       :: EqualAI
            E + Q ⇌ EQ                       :: OnlyA
            EQ + P <--> EPQ                  :: OnlyA
        end
    end
    a_q_p = @_testhelper_unguarded @enzyme_mechanism begin
        substrates: A, B ; products: P, Q ; oligomeric_state: 1
        steps: begin
            (E + A <--> EA, EB + A <--> EPQ) :: OnlyA
            (E + B ⇌ EB, EA + B ⇌ EPQ)       :: EqualAI
            E + Q ⇌ EQ                       :: OnlyA
            EQ + P <--> EPQ                  :: OnlyA
        end
    end
    b_q_p = @_testhelper_unguarded @enzyme_mechanism begin
        substrates: A, B ; products: P, Q ; oligomeric_state: 1
        steps: begin
            (E + A <--> EA, EB + A <--> EPQ) :: EqualAI
            (E + B ⇌ EB, EA + B ⇌ EPQ)       :: OnlyA
            E + Q ⇌ EQ                       :: OnlyA
            EQ + P <--> EPQ                  :: OnlyA
        end
    end
    a_b_q_p = @_testhelper_unguarded @enzyme_mechanism begin
        substrates: A, B ; products: P, Q ; oligomeric_state: 1
        steps: begin
            (E + A <--> EA, EB + A <--> EPQ) :: OnlyA
            (E + B ⇌ EB, EA + B ⇌ EPQ)       :: OnlyA
            E + Q ⇌ EQ                       :: OnlyA
            EQ + P <--> EPQ                  :: OnlyA
        end
    end
    a_r_only_a = @_testhelper_unguarded @enzyme_mechanism begin
        substrates: A, B ; products: P, Q ; oligomeric_state: 1
        steps: begin
            (E + A <--> EA, EB + A <--> EPQ) :: OnlyA
            (E + B ⇌ EB, EA + B ⇌ EPQ)       :: EqualAI
            E + Q ⇌ EQ                       :: EqualAI
            EQ + P <--> EPQ                  :: EqualAI
        end
        regulatory_site(1): R::OnlyA
    end
    a_r_only_i = @_testhelper_unguarded @enzyme_mechanism begin
        substrates: A, B ; products: P, Q ; oligomeric_state: 1
        steps: begin
            (E + A <--> EA, EB + A <--> EPQ) :: OnlyA
            (E + B ⇌ EB, EA + B ⇌ EPQ)       :: EqualAI
            E + Q ⇌ EQ                       :: EqualAI
            EQ + P <--> EPQ                  :: EqualAI
        end
        regulatory_site(1): R::OnlyI
    end
    b_r_only_a = @_testhelper_unguarded @enzyme_mechanism begin
        substrates: A, B ; products: P, Q ; oligomeric_state: 1
        steps: begin
            (E + A <--> EA, EB + A <--> EPQ) :: EqualAI
            (E + B ⇌ EB, EA + B ⇌ EPQ)       :: OnlyA
            E + Q ⇌ EQ                       :: EqualAI
            EQ + P <--> EPQ                  :: EqualAI
        end
        regulatory_site(1): R::OnlyA
    end
    b_r_only_i = @_testhelper_unguarded @enzyme_mechanism begin
        substrates: A, B ; products: P, Q ; oligomeric_state: 1
        steps: begin
            (E + A <--> EA, EB + A <--> EPQ) :: EqualAI
            (E + B ⇌ EB, EA + B ⇌ EPQ)       :: OnlyA
            E + Q ⇌ EQ                       :: EqualAI
            EQ + P <--> EPQ                  :: EqualAI
        end
        regulatory_site(1): R::OnlyI
    end
    a_b_r_only_a = @_testhelper_unguarded @enzyme_mechanism begin
        substrates: A, B ; products: P, Q ; oligomeric_state: 1
        steps: begin
            (E + A <--> EA, EB + A <--> EPQ) :: OnlyA
            (E + B ⇌ EB, EA + B ⇌ EPQ)       :: OnlyA
            E + Q ⇌ EQ                       :: EqualAI
            EQ + P <--> EPQ                  :: EqualAI
        end
        regulatory_site(1): R::OnlyA
    end
    a_b_r_only_i = @_testhelper_unguarded @enzyme_mechanism begin
        substrates: A, B ; products: P, Q ; oligomeric_state: 1
        steps: begin
            (E + A <--> EA, EB + A <--> EPQ) :: OnlyA
            (E + B ⇌ EB, EA + B ⇌ EPQ)       :: OnlyA
            E + Q ⇌ EQ                       :: EqualAI
            EQ + P <--> EPQ                  :: EqualAI
        end
        regulatory_site(1): R::OnlyI
    end
    a_nonequal = @enzyme_mechanism begin
        substrates: A, B ; products: P, Q ; oligomeric_state: 1
        steps: begin
            (E + A <--> EA, EB + A <--> EPQ) :: NonequalAI
            (E + B ⇌ EB, EA + B ⇌ EPQ)       :: EqualAI
            E + Q ⇌ EQ                       :: EqualAI
            EQ + P <--> EPQ                  :: EqualAI
        end
    end
    kids = _testhelper_move(EnzymeRates.Enumeration._expand_to_allosteric, m, rxn)
    @test length(kids) == 22
    @test Set(kids) == Set([a, b, a_b, q, a_q, b_q, a_b_q, p, a_p, b_p, a_b_p, q_p, a_q_p,
                            b_q_p, a_b_q_p, a_r_only_a, a_r_only_i, b_r_only_a, b_r_only_i,
                            a_b_r_only_a, a_b_r_only_i, a_nonequal])
    @test !EnzymeRates.Enumeration._enumerator_returns(_testhelper_candidate(a))
    # The parent has 5 fitted constants, all identifiable. Each child with A's and B's
    # groups :OnlyA adds L. Its dead inactive conformation reaches no EB, since B's group
    # is :OnlyA, so it adds L·E to E's rapid-equilibrium segment {E, EB, EQ}, times
    # 1 + Q/K_Q where Q stays :EqualAI. A's group takes A up at E and at EB under one rate
    # constant, so the denominator's A·B term is 2·B/K_B times its A term whatever L is:
    # K_B is pinned and cannot also absorb the 1 + L at E. L shows in every such child,
    # the two whose inactive conformation binds nothing included: 6 fitted, rank 6.
    @test _testhelper_fitted(m) == _testhelper_identifiable_rank(m) == 5
    for k in (a_b_q, a_b_p, a_b_q_p)
        @test _testhelper_fitted(k) == _testhelper_identifiable_rank(k) == 6
    end
    # None is an L phantom: L raises every such child's rank, so no certificate exists,
    # though a_b_q and a_b_q_p bind nothing in the inactive conformation.
    for k in (a_b_q, a_b_q_p)
        @test !any(
            EnzymeRates.Mechanisms.steps_running_in_allosteric_conformations(k).inactive)
    end
    @test !any(k -> EnzymeRates.Derivation._is_phantom(EnzymeRates.Derivation._memo_law(k)),
               (a_b_q, a_b_p, a_b_q_p))
end

@testset "candidate_mechanisms expands a parent whose binding changes conformation" begin
    rxn = @enzyme_reaction begin
        substrates: A[CX], B[N]
        products: P[C], Q[NX]
    end
    # A binding step may change conformation; only the residual is chemistry.
    conf = @enzyme_mechanism begin
        substrates: A, B
        products: P, Q
        steps: begin
            E + A ⇌ EstarA
            EstarA + B ⇌ EstarAB
            EstarAB <--> EstarPQ
            EstarPQ ⇌ EstarP + Q
            EstarP ⇌ E + P
        end
    end
    @test _testhelper_children(rxn, [conf]) isa Vector
end

@testset "candidate_mechanisms rejects a parent with a zero-flux steady-state group" begin
    # E + A and E + Q at steady state while EQ + A ⇌ EAQ ⇌ EA + Q stays at
    # rapid equilibrium: {E} joins the rest by two edges of weight 0, so neither
    # group carries flux and their constants enter the rate only as ratios. The
    # flip never emits this child; a hand-written parent is refused.
    rxn = @enzyme_reaction begin
        substrates: A[C], B[N]
        products: P[C], Q[N]
    end
    shunt = @enzyme_mechanism begin
        substrates: A, B
        products: P, Q
        steps: begin
            E + A <--> EA
            EQ + A ⇌ EAQ
            E + Q <--> EQ
            EA + Q ⇌ EAQ
            EA + B ⇌ EAB
            EAB <--> EPQ
            EQ + P ⇌ EPQ
        end
    end
    err = try
        _testhelper_children(rxn, [shunt]); nothing
    catch e
        e
    end
    @test err isa ErrorException
    @test occursin("carries net flux", err.msg)
    @test occursin("{E + A <--> EA}", err.msg)
end

@testset "candidate_mechanisms rejects a parent whose copy group is redundant" begin
    # A bound as its own competitive inhibitor at E alone: EAmod has the
    # composition of EA, formed by the A group alone, so the gauge exists and the
    # copy's constant enters the rate only added to K_A's. The dead-end move never
    # emits this child; a hand-written parent is refused.
    rxn = @enzyme_reaction begin
        substrates: A[C], B[N]
        products: P[C], Q[N]
        competitive_inhibitors: A
    end
    twin_only = @enzyme_mechanism begin
        substrates: A, B
        products: P, Q
        steps: begin
            E + A ⇌ EA
            EA + B ⇌ EAB
            EAB <--> EPQ
            EQ + P ⇌ EPQ
            E + Q ⇌ EQ
            E + A ⇌ EAmod
        end
    end
    err = try
        _testhelper_children(rxn, [twin_only]); nothing
    catch e
        e
    end
    @test err isa ErrorException
    @test occursin("duplicates a productive form", err.msg)
    @test occursin("dwell gauge", err.msg)
    # The same copy placed where it also creates a new complex is a valid parent.
    kept = @enzyme_mechanism begin
        substrates: A, B
        products: P, Q
        steps: begin
            E + A ⇌ EA
            EA + B ⇌ EAB
            EAB <--> EPQ
            EQ + P ⇌ EPQ
            (E + Q ⇌ EQ, EAmod + Q ⇌ EAmodQ)
            (E + A ⇌ EAmod, EQ + A ⇌ EAmodQ)
        end
    end
    @test !isempty(_testhelper_children(rxn, [kept]))
    @test EnzymeRates.Enumeration._assert_emission_rules(
        _testhelper_candidate(kept)) === nothing
    # A copy at E alone whose A group also binds at EQ, forming EAQ: every
    # complex duplicates a productive form, but the gauge fails, so it is a valid
    # parent.
    no_gauge = @enzyme_mechanism begin
        substrates: A, B
        products: P, Q
        steps: begin
            (E + A ⇌ EA, EQ + A ⇌ EAQ)
            EA + B ⇌ EAB
            EAB <--> EPQ
            EQ + P ⇌ EPQ
            (E + Q ⇌ EQ, EA + Q ⇌ EAQ)
            E + A ⇌ EAmod
        end
    end
    @test EnzymeRates.Enumeration._assert_emission_rules(
        _testhelper_candidate(no_gauge)) === nothing
    @test !isempty(_testhelper_children(rxn, [no_gauge]))
end

@testset "candidate_mechanisms checks that every child's rate law is defined" begin
    # Random product release at rapid equilibrium between steady-state releases leaves
    # the segment {EP, EPQ, EQ} without a bottom form. The macro rejects such a
    # mechanism; the constructor, which checks the representation only, builds it from
    # the walker's raw value. Its dead-end children keep the segment, the first of them
    # with the copy's complex EIQ in it, and `candidate_mechanisms` drops them.
    rxn = @enzyme_reaction begin
        substrates: S[C2]
        products: P[C], Q[C]
        competitive_inhibitors: I
    end
    parent = EnzymeRates.DSL._construct(EnzymeRates.DSL.raw_mechanism(:(begin
        substrates: S
        products: P, Q
        steps: begin
            E + S ⇌ ES
            ES <--> EPQ
            EPQ ⇌ EQ + P
            EPQ ⇌ EP + Q
            EP <--> E + P
            EQ <--> E + Q
        end
    end)))
    @test !isempty(_testhelper_move(EnzymeRates.Enumeration._expand_add_dead_end_regulator,
                                    parent, rxn))
    # The exit keeps the two flips that cut the segment, one product release steady state.
    expected = [
        (@enzyme_mechanism begin
            substrates: S
            products: P, Q
            steps: begin
                E + S ⇌ ES
                ES <--> EPQ
                EPQ ⇌ EQ + P
                EPQ <--> EP + Q
                EP <--> E + P
                EQ <--> E + Q
            end
        end),
        (@enzyme_mechanism begin
            substrates: S
            products: P, Q
            steps: begin
                E + S ⇌ ES
                ES <--> EPQ
                EPQ <--> EQ + P
                EPQ ⇌ EP + Q
                EP <--> E + P
                EQ <--> E + Q
            end
        end)]
    kids = _testhelper_children(rxn, [parent])
    @test length(kids) == 2
    @test Set(kids) == Set(expected)
end

@testset "candidate_mechanisms binds an inhibitor I beside a substrate Iinh" begin
    # The forms EI and EIinh render apart (EI, EIinh), so the dead-end child holds
    # both. The other five children tag the seed's groups :OnlyA (no oligomeric state is
    # declared, so they keep one catalytic subunit): the chemistry alone or with the
    # substrate binding, the product binding or both, and the two bindings without the
    # chemistry. One binding alone :OnlyA is refused (`_FAST_STEP_PINNED`).
    rxn = @enzyme_reaction begin
        substrates: Iinh[C]
        products: P[C]
        competitive_inhibitors: I
    end
    seed = @enzyme_mechanism begin
        substrates: Iinh
        products: P
        steps: begin
            E + Iinh ⇌ EIinh
            EIinh <--> EP
            EP ⇌ E + P
        end
    end
    @test _testhelper_mechanisms(EnzymeRates.Enumeration.init_mechanisms(rxn)) == [seed]
    children = _testhelper_children(rxn, [seed])
    dead_end = @enzyme_mechanism begin
        substrates: Iinh
        products: P
        modifiers: I
        steps: begin
            E + Iinh ⇌ EIinh
            EIinh <--> EP
            EP ⇌ E + P
            E + I ⇌ EI
        end
    end
    iinh_only_a = @enzyme_mechanism begin
        substrates: Iinh
        products: P
        oligomeric_state: 1
        steps: begin
            E + Iinh ⇌ EIinh      :: OnlyA
            EIinh <--> EP       :: OnlyA
            EP ⇌ E + P            :: EqualAI
        end
    end
    p_only_a = @enzyme_mechanism begin
        substrates: Iinh
        products: P
        oligomeric_state: 1
        steps: begin
            E + Iinh ⇌ EIinh      :: EqualAI
            EIinh <--> EP       :: OnlyA
            EP ⇌ E + P            :: OnlyA
        end
    end
    both_only_a = @enzyme_mechanism begin
        substrates: Iinh
        products: P
        oligomeric_state: 1
        steps: begin
            E + Iinh ⇌ EIinh      :: OnlyA
            EIinh <--> EP       :: OnlyA
            EP ⇌ E + P            :: OnlyA
        end
    end
    chemistry_only_a = @enzyme_mechanism begin
        substrates: Iinh
        products: P
        oligomeric_state: 1
        steps: begin
            E + Iinh ⇌ EIinh      :: EqualAI
            EIinh <--> EP       :: OnlyA
            EP ⇌ E + P            :: EqualAI
        end
    end
    bindings_only_a = @enzyme_mechanism begin
        substrates: Iinh
        products: P
        oligomeric_state: 1
        steps: begin
            E + Iinh ⇌ EIinh      :: OnlyA
            EIinh <--> EP       :: EqualAI
            EP ⇌ E + P            :: OnlyA
        end
    end
    @test length(children) == 6
    @test Set(children) == Set([dead_end, iinh_only_a, p_only_a, both_only_a,
                                chemistry_only_a, bindings_only_a])
end

@testset "_expand_add_dead_end_regulator: inhibitor mirrors one half-reaction" begin
    # The inhibitor binds the forms that bind a competing ligand, never a form
    # already carrying one, so a competing substrate's binding step is never
    # mirrored and no inhibitor-bound branch completes the net reaction. With A
    # binding EQ and P binding E(B; res) as dead ends, the pattern "compete
    # with A and P" puts I on E, EQ, E(; res) and E(B; res): the second
    # half-reaction (B binds, chemistry, Q leaves) runs with I bound, and A
    # cannot bind EI. Six patterns give distinct targets.
    rxn = @enzyme_reaction begin
        substrates: A[CX], B[N]
        products: P[C], Q[NX]
        competitive_inhibitors: I
    end
    m = @enzyme_mechanism begin
        substrates: A, B
        products: P, Q
        steps: begin
            E + A ⇌ EA
            EA <--> EP_res_plusA_minusP
            EP_res_plusA_minusP ⇌ E_res_plusA_minusP + P
            E_res_plusA_minusP + B ⇌ EB_res_plusA_minusP
            EB_res_plusA_minusP <--> EQ
            EQ ⇌ E + Q
            EQ + A ⇌ EAQ
            EB_res_plusA_minusP + P ⇌ EBP_res_plusA_minusP
        end
    end
    second_half = @enzyme_mechanism begin      # I competes with A and P
        substrates: A, B
        products: P, Q
        modifiers: I
        steps: begin
            E + A ⇌ EA
            EA <--> EP_res_plusA_minusP
            EP_res_plusA_minusP ⇌ E_res_plusA_minusP + P
            (E_res_plusA_minusP + B ⇌ EB_res_plusA_minusP,
             EI_res_plusA_minusP + B ⇌ EBI_res_plusA_minusP)
            (EB_res_plusA_minusP <--> EQ,
             EBI_res_plusA_minusP <--> EIQ)
            (EQ ⇌ E + Q, EIQ ⇌ EI + Q)
            EQ + A ⇌ EAQ
            EB_res_plusA_minusP + P ⇌ EBP_res_plusA_minusP
            (E + I ⇌ EI, EQ + I ⇌ EIQ,
             E_res_plusA_minusP + I ⇌ EI_res_plusA_minusP,
             EB_res_plusA_minusP + I ⇌ EBI_res_plusA_minusP)
        end
    end
    only_E = @enzyme_mechanism begin           # I competes with A and Q
        substrates: A, B
        products: P, Q
        modifiers: I
        steps: begin
            E + A ⇌ EA
            EA <--> EP_res_plusA_minusP
            EP_res_plusA_minusP ⇌ E_res_plusA_minusP + P
            E_res_plusA_minusP + B ⇌ EB_res_plusA_minusP
            EB_res_plusA_minusP <--> EQ
            EQ ⇌ E + Q
            EQ + A ⇌ EAQ
            EB_res_plusA_minusP + P ⇌ EBP_res_plusA_minusP
            E + I ⇌ EI
        end
    end
    B_binding = @enzyme_mechanism begin        # I competes with A, P and Q
        substrates: A, B
        products: P, Q
        modifiers: I
        steps: begin
            E + A ⇌ EA
            EA <--> EP_res_plusA_minusP
            EP_res_plusA_minusP ⇌ E_res_plusA_minusP + P
            (E_res_plusA_minusP + B ⇌ EB_res_plusA_minusP,
             EI_res_plusA_minusP + B ⇌ EBI_res_plusA_minusP)
            EB_res_plusA_minusP <--> EQ
            EQ ⇌ E + Q
            EQ + A ⇌ EAQ
            EB_res_plusA_minusP + P ⇌ EBP_res_plusA_minusP
            (E + I ⇌ EI,
             E_res_plusA_minusP + I ⇌ EI_res_plusA_minusP,
             EB_res_plusA_minusP + I ⇌ EBI_res_plusA_minusP)
        end
    end
    only_Eres = @enzyme_mechanism begin        # I competes with B and P
        substrates: A, B
        products: P, Q
        modifiers: I
        steps: begin
            E + A ⇌ EA
            EA <--> EP_res_plusA_minusP
            EP_res_plusA_minusP ⇌ E_res_plusA_minusP + P
            E_res_plusA_minusP + B ⇌ EB_res_plusA_minusP
            EB_res_plusA_minusP <--> EQ
            EQ ⇌ E + Q
            EQ + A ⇌ EAQ
            EB_res_plusA_minusP + P ⇌ EBP_res_plusA_minusP
            E_res_plusA_minusP + I ⇌ EI_res_plusA_minusP
        end
    end
    E_and_Eres = @enzyme_mechanism begin       # I competes with B and Q
        substrates: A, B
        products: P, Q
        modifiers: I
        steps: begin
            E + A ⇌ EA
            EA <--> EP_res_plusA_minusP
            EP_res_plusA_minusP ⇌ E_res_plusA_minusP + P
            E_res_plusA_minusP + B ⇌ EB_res_plusA_minusP
            EB_res_plusA_minusP <--> EQ
            EQ ⇌ E + Q
            EQ + A ⇌ EAQ
            EB_res_plusA_minusP + P ⇌ EBP_res_plusA_minusP
            (E + I ⇌ EI, E_res_plusA_minusP + I ⇌ EI_res_plusA_minusP)
        end
    end
    Q_release = @enzyme_mechanism begin        # I competes with A, B and P
        substrates: A, B
        products: P, Q
        modifiers: I
        steps: begin
            E + A ⇌ EA
            EA <--> EP_res_plusA_minusP
            EP_res_plusA_minusP ⇌ E_res_plusA_minusP + P
            E_res_plusA_minusP + B ⇌ EB_res_plusA_minusP
            EB_res_plusA_minusP <--> EQ
            (EQ ⇌ E + Q, EIQ ⇌ EI + Q)
            EQ + A ⇌ EAQ
            EB_res_plusA_minusP + P ⇌ EBP_res_plusA_minusP
            (E + I ⇌ EI, EQ + I ⇌ EIQ,
             E_res_plusA_minusP + I ⇌ EI_res_plusA_minusP)
        end
    end
    expected = [second_half, only_E, B_binding, only_Eres, E_and_Eres, Q_release]
    kids = _testhelper_move(EnzymeRates.Enumeration._expand_add_dead_end_regulator, m, rxn)
    @test length(kids) == 6
    @test Set(kids) == Set(expected)
    @test all(c -> _testhelper_keeps_reaction(c, m) &&
              :I in _testhelper_modifier_names(c),
              kids)
end

@testset "_expand_add_dead_end_regulator: inhibitor never runs the net reaction" begin
    # Aggregate pin over the ping-pong seed set: in every regulator child some
    # substrate has no binding step inside the inhibitor-bound branch (the steps
    # whose two forms both carry I), so no cycle of inhibitor-bound forms
    # completes the reaction. Follows from competition with at least one
    # substrate: a form carrying a competing ligand never receives I, so that
    # ligand's binding step is never mirrored.
    rxn = @enzyme_reaction begin
        substrates: A[CX], B[N]
        products: P[C], Q[NX]
        competitive_inhibitors: I
    end
    # Whether form `f` of graph `c` carries a modifier, read off its record.
    _testhelper_inhibitor_bound(c, f) = !isempty(_testhelper_candidate(c).forms[f].copies)
    n_children = 0; n_half = 0
    seeds = _testhelper_mechanisms(EnzymeRates.Enumeration.init_mechanisms(rxn))
    for m in seeds, c in _testhelper_move(
        EnzymeRates.Enumeration._expand_add_dead_end_regulator, m, rxn)
        branch = [e for e in c.catalytic.steps if _testhelper_inhibitor_bound(c, e.from) &&
                                        _testhelper_inhibitor_bound(c, e.to)]
        bound_in_branch = Set(c.catalytic.metabolites[only(e.taken_up)]
                              for e in branch if EnzymeRates.Enumeration._is_binding(e))
        @test !(Set([:A, :B]) ⊆ bound_in_branch)
        n_children += 1
        any(EnzymeRates.Mechanisms._is_iso, branch) && (n_half += 1)
    end
    @test n_children > 0
    @test n_half > 0          # the half-reaction mirror really occurs on the seeds
end

@testset "to_allosteric shares :EqualAI downstream constants across conformations" begin
    # The inactive conformation drops each `:OnlyA` binding step, and its law holds only
    # the forms its untagged steps join to free E: that ligand's downstream complex holds
    # inactive enzyme only if other untagged steps join it to free E.
    # A lumped `:EqualAI` group is one class in both conformations, named after the
    # kept member of the full mechanism, so the inactive conformation cannot give it a
    # second name and the joint solve counts it once.
    #
    # `_expand_to_allosteric` tags every nonempty set of groups `:OnlyA`, so the case
    # this guards is the child that sets the free-enzyme bindings (a substrate/product
    # pair on free E) and every chemical step `:OnlyA`, a dead inactive conformation,
    # leaving the lumped downstream group's binding `:EqualAI`, outside the inactive law,
    # which holds only the forms its untagged steps join to free E; the constant stays
    # merged, so the child adds only `L` — parent + 1 fitted params. A child that instead
    # sets a downstream binding `:OnlyA` orphans a bound form of the lumped group; its
    # multi-`:OnlyA` derivation is a separate, documented concern and is not asserted
    # here.
    rxn = @enzyme_reaction begin
        substrates: NADH[C21H29N7O14P2], Pyruvate[C3H4O3]
        products:   Lactate[C3H6O3], NAD[C21H27N7O14P2]
        oligomeric_state: 4
    end
    # A mechanism exercises the path only when a downstream kinetic group lumps
    # binding to ≥2 bound complexes (no free-E member) — the split the `:OnlyA`
    # orphaning would otherwise break.
    # A form's record holds nothing bound.
    _testhelper_empty(m, f) = isempty(_testhelper_bound(_testhelper_candidate(m).forms[f]))
    _testhelper_has_lumped_bound_group(m) =
        any(_testhelper_group_members(m.catalytic)) do ks
        length(ks) ≥ 2 && all(k -> !_testhelper_empty(m, m.catalytic.steps[k].from), ks)
    end
    # Whether the first step of a group, edges `ks` of `m`, binds into a form that holds
    # nothing bound.
    _testhelper_is_free_e_binding(m, ks) = (e = m.catalytic.steps[first(ks)];
        _testhelper_empty(m, e.from) && EnzymeRates.Enumeration._is_binding(e))
    n_reproducers = 0
    for m in _testhelper_mechanisms(EnzymeRates.Enumeration.init_mechanisms(rxn))
        _testhelper_has_lumped_bound_group(m) || continue
        pn = _testhelper_fitted(m)
        free_e = Set(g for (g, ks) in enumerate(_testhelper_group_members(m.catalytic))
                     if _testhelper_is_free_e_binding(m, ks))
        iso = Set(_testhelper_iso_groups(m))
        children = _testhelper_move(EnzymeRates.Enumeration._expand_to_allosteric, m, rxn)
        @test all(_testhelper_is_allosteric, children)
        # The free-enzyme-normalization child (the free-E bindings :OnlyA, plus
        # all chemical steps :OnlyA — a dead inactive conformation) leaves the
        # lumped downstream group's binding :EqualAI, so its shared constant
        # stays merged → Δ = 1.
        clean = filter(children) do c
            Set(g for g in 1:length(c.allosteric.allosteric_states)
                if c.allosteric.allosteric_states[g] == :OnlyA) ==
                union(free_e, iso)
        end
        @test length(clean) == 1
        @test _testhelper_fitted(only(clean)) == pn + 1
        n_reproducers += 1
        n_reproducers ≥ 3 && break
    end
    @test n_reproducers ≥ 1
end

# Tag of each catalytic kinetic group, keyed by the metabolite it binds
# (`:chem` for the chemical step). Group index is never a stable key: the
# Mechanism constructor canonicalizes group order.
_testhelper_tags_by_bound_metabolite(x) = Dict(
    (s = first(steps); _testhelper_is_binding(s) ? only(s.consumed) : :chem) =>
        x.allosteric.allosteric_states[g]
    for (g, steps) in enumerate(_testhelper_groups(_testhelper_graph(x))))

# ─── _expand_change_allo_state admission filter ────────────────────────
@testset "_expand_change_allo_state drops the chemical step's relaxation" begin
    ER = EnzymeRates

    # Uni-uni S ⇌ P. The S binding's :OnlyA is one-sided; only the :OnlyA
    # chemical step (`k_I = 0`) makes the parent legal. Relaxing that
    # chemical step to :NonequalAI restores a finite `k_I` and strands the
    # binding, so that child has no thermodynamic reading and is dropped.
    parent = @enzyme_mechanism begin
        substrates: S
        products: P
        oligomeric_state: 2
        steps: begin
            E + S ⇌ ES      :: OnlyA
            ES <--> EP    :: OnlyA
            E + P ⇌ EP      :: EqualAI
        end
    end

    result = _testhelper_move(ER.Enumeration._expand_change_allo_state, parent)
    foreach(_testhelper_assert_mechanism_invariants, result)
    # Three tags are relaxable; the chemical step's relaxation is dropped.
    @test length(result) == 2
    @test !any(r -> _testhelper_tags_by_bound_metabolite(r)[:chem] == :NonequalAI, result)
    @test Set(_testhelper_tags_by_bound_metabolite(r) for r in result) == Set([
        Dict(:S => :OnlyA, :chem => :OnlyA, :P => :NonequalAI),
        Dict(:S => :NonequalAI, :chem => :OnlyA, :P => :EqualAI)])

    # Balanced parent: both bindings and the chemical step :OnlyA — the inactive
    # conformation binds nothing (a fully-inert T-state). Relaxing either binding
    # is retained (it leaves all chemical steps :OnlyA). The chemical step tagged
    # :NonequalAI adds no fitted constant: the inactive conformation, binding nothing,
    # never reaches ES or EP, so k_I enters nowhere, and the move drops that child.
    balanced = @enzyme_mechanism begin
        substrates: S
        products: P
        oligomeric_state: 2
        steps: begin
            E + S ⇌ ES      :: OnlyA
            ES <--> EP    :: OnlyA
            E + P ⇌ EP      :: OnlyA
        end
    end
    kids = _testhelper_move(ER.Enumeration._expand_change_allo_state, balanced)
    @test length(kids) == 2           # the chemical-step relaxation is dropped
    @test !any(_testhelper_tags_by_bound_metabolite(k)[:chem] == :NonequalAI for k in kids)
    @test Set(_testhelper_tags_by_bound_metabolite(k) for k in kids) == Set([
        Dict(:S => :NonequalAI, :chem => :OnlyA, :P => :OnlyA),
        Dict(:S => :OnlyA, :chem => :OnlyA, :P => :NonequalAI)])

    # A regulatory ligand's tag does not enter admission —
    # a regulator site completes no catalytic cycle — so relaxing one can
    # never strand an :OnlyA binding; it adds one constant, and the move keeps it.
    reg_parent = @enzyme_mechanism begin
        substrates: S
        products: P
        oligomeric_state: 2
        steps: begin
            E + S ⇌ ES      :: OnlyA
            ES <--> EP    :: OnlyA
            E + P ⇌ EP      :: EqualAI
        end
        regulatory_site(2): R::OnlyA
    end
    reg_kids = _testhelper_move(ER.Enumeration._expand_change_allo_state, reg_parent)
    # 2 surviving cat relaxations + the R-ligand relaxation, which survives.
    @test length(reg_kids) == 3
    @test count(k -> only(k.allosteric.regulatory_sites).allosteric_states ==
                     [:NonequalAI], reg_kids) == 1
end

@testset "change_allo_state tags each group of a three-step ping-pong alone" begin
    # The three-step ping-pong (census pingpong seed 8fe71fded09a27ec) with the second
    # half-reaction and Q's binding :OnlyA, the first half-reaction running in both
    # conformations. The inactive conformation turns nothing over, so each group tagged
    # :NonequalAI alone adds constants, and the guard accepts each child: Q's binding,
    # the first half-reaction (a Theorell–Chance step) and the second (B's fused binding
    # into EQ). Three children.
    parent = @_testhelper_unguarded @enzyme_mechanism begin
        substrates: A, B ; products: P, Q ; oligomeric_state: 1
        steps: begin
            E + A <--> E_res_plusA_minusP + P :: EqualAI
            E + Q <--> EQ                      :: OnlyA
            E_res_plusA_minusP + B <--> EQ  :: OnlyA
        end
    end
    q = @enzyme_mechanism begin
        substrates: A, B ; products: P, Q ; oligomeric_state: 1
        steps: begin
            E + A <--> E_res_plusA_minusP + P :: EqualAI
            E + Q <--> EQ                      :: NonequalAI
            E_res_plusA_minusP + B <--> EQ  :: OnlyA
        end
    end
    first_half = @enzyme_mechanism begin
        substrates: A, B ; products: P, Q ; oligomeric_state: 1
        steps: begin
            E + A <--> E_res_plusA_minusP + P :: NonequalAI
            E + Q <--> EQ                      :: OnlyA
            E_res_plusA_minusP + B <--> EQ  :: OnlyA
        end
    end
    second_half = @enzyme_mechanism begin
        substrates: A, B ; products: P, Q ; oligomeric_state: 1
        steps: begin
            E + A <--> E_res_plusA_minusP + P :: EqualAI
            E + Q <--> EQ                      :: OnlyA
            E_res_plusA_minusP + B <--> EQ  :: NonequalAI
        end
    end
    kids = _testhelper_move(EnzymeRates.Enumeration._expand_change_allo_state, parent)
    @test length(kids) == 3
    @test Set(kids) == Set([q, first_half, second_half])
    @test _testhelper_param_deltas(parent, kids) == [2, 2, 2]
    @test _testhelper_is_valid(parent)
end

@testset "enumeration emits only children the guard admits" begin
    # A plain ping-pong (census 3fdbe6e1dc001bb7) whose A binding is fused into the
    # product-bound covalent form and whose B binding is fused into EQ: A's and B's groups
    # are the chemistry groups. Its one rapid-equilibrium group flips to steady state:
    # one plain child. The tagging move emits the 15 nonempty sets of the four groups
    # :OnlyA, and A's group, B's and both again with R :OnlyA and :OnlyI: 21 children,
    # none :NonequalAI, since each minimal set that gains is returned tagged :OnlyA. The
    # guard rejects A's fused binding alone :OnlyA, with R or without: the inactive
    # conformation still joins its two sides (`_FAST_STEP_PINNED`). The exit keeps the
    # other 18 children and the flip. A second census ping-pong, whose A binding is fused
    # at steady state and whose B binding is fused at rapid equilibrium, shows more drops.
    # Its move emits 22 children: the same 21 sets, and P's and B's bindings tagged
    # :NonequalAI together, a minimal gaining set refused tagged :OnlyA (`_BRIDGE`). The
    # exit keeps 11 of them and the two flips. Among those it drops is A's, P's and B's
    # groups :OnlyA, where at zero A and B the fused B binding fills the covalent form
    # E_res_plusA_minusP faster than it empties (`_FILLED`).
    rxn = @enzyme_reaction begin
        substrates: A[CX], B[N]
        products: P[C], Q[NX]
        allosteric_regulators: R
        oligomeric_state: 2
        shared_catalytic_site: (A, P)
    end
    parent = @enzyme_mechanism begin
        substrates: A, B
        products: P, Q
        steps: begin
            E + A ⇌ EP_res_plusA_minusP
            E_res_plusA_minusP + P <--> EP_res_plusA_minusP
            E_res_plusA_minusP + B <--> EQ
            E + Q <--> EQ
        end
    end
    flip = @enzyme_mechanism begin
        substrates: A, B
        products: P, Q
        steps: begin
            E + A <--> EP_res_plusA_minusP
            E_res_plusA_minusP + P <--> EP_res_plusA_minusP
            E_res_plusA_minusP + B <--> EQ
            E + Q <--> EQ
        end
    end
    a = @_testhelper_unguarded @enzyme_mechanism begin
        substrates: A, B ; products: P, Q ; oligomeric_state: 2
        steps: begin
            E + A ⇌ EP_res_plusA_minusP                     :: OnlyA
            E_res_plusA_minusP + P <--> EP_res_plusA_minusP :: EqualAI
            E_res_plusA_minusP + B <--> EQ                  :: EqualAI
            E + Q <--> EQ                                   :: EqualAI
        end
    end
    p = @_testhelper_unguarded @enzyme_mechanism begin
        substrates: A, B ; products: P, Q ; oligomeric_state: 2
        steps: begin
            E + A ⇌ EP_res_plusA_minusP                     :: EqualAI
            E_res_plusA_minusP + P <--> EP_res_plusA_minusP :: OnlyA
            E_res_plusA_minusP + B <--> EQ                  :: EqualAI
            E + Q <--> EQ                                   :: EqualAI
        end
    end
    a_p = @_testhelper_unguarded @enzyme_mechanism begin
        substrates: A, B ; products: P, Q ; oligomeric_state: 2
        steps: begin
            E + A ⇌ EP_res_plusA_minusP                     :: OnlyA
            E_res_plusA_minusP + P <--> EP_res_plusA_minusP :: OnlyA
            E_res_plusA_minusP + B <--> EQ                  :: EqualAI
            E + Q <--> EQ                                   :: EqualAI
        end
    end
    b = @_testhelper_unguarded @enzyme_mechanism begin
        substrates: A, B ; products: P, Q ; oligomeric_state: 2
        steps: begin
            E + A ⇌ EP_res_plusA_minusP                     :: EqualAI
            E_res_plusA_minusP + P <--> EP_res_plusA_minusP :: EqualAI
            E_res_plusA_minusP + B <--> EQ                  :: OnlyA
            E + Q <--> EQ                                   :: EqualAI
        end
    end
    a_b = @_testhelper_unguarded @enzyme_mechanism begin
        substrates: A, B ; products: P, Q ; oligomeric_state: 2
        steps: begin
            E + A ⇌ EP_res_plusA_minusP                     :: OnlyA
            E_res_plusA_minusP + P <--> EP_res_plusA_minusP :: EqualAI
            E_res_plusA_minusP + B <--> EQ                  :: OnlyA
            E + Q <--> EQ                                   :: EqualAI
        end
    end
    p_b = @_testhelper_unguarded @enzyme_mechanism begin
        substrates: A, B ; products: P, Q ; oligomeric_state: 2
        steps: begin
            E + A ⇌ EP_res_plusA_minusP                     :: EqualAI
            E_res_plusA_minusP + P <--> EP_res_plusA_minusP :: OnlyA
            E_res_plusA_minusP + B <--> EQ                  :: OnlyA
            E + Q <--> EQ                                   :: EqualAI
        end
    end
    a_p_b = @_testhelper_unguarded @enzyme_mechanism begin
        substrates: A, B ; products: P, Q ; oligomeric_state: 2
        steps: begin
            E + A ⇌ EP_res_plusA_minusP                     :: OnlyA
            E_res_plusA_minusP + P <--> EP_res_plusA_minusP :: OnlyA
            E_res_plusA_minusP + B <--> EQ                  :: OnlyA
            E + Q <--> EQ                                   :: EqualAI
        end
    end
    q = @_testhelper_unguarded @enzyme_mechanism begin
        substrates: A, B ; products: P, Q ; oligomeric_state: 2
        steps: begin
            E + A ⇌ EP_res_plusA_minusP                     :: EqualAI
            E_res_plusA_minusP + P <--> EP_res_plusA_minusP :: EqualAI
            E_res_plusA_minusP + B <--> EQ                  :: EqualAI
            E + Q <--> EQ                                   :: OnlyA
        end
    end
    a_q = @_testhelper_unguarded @enzyme_mechanism begin
        substrates: A, B ; products: P, Q ; oligomeric_state: 2
        steps: begin
            E + A ⇌ EP_res_plusA_minusP                     :: OnlyA
            E_res_plusA_minusP + P <--> EP_res_plusA_minusP :: EqualAI
            E_res_plusA_minusP + B <--> EQ                  :: EqualAI
            E + Q <--> EQ                                   :: OnlyA
        end
    end
    p_q = @_testhelper_unguarded @enzyme_mechanism begin
        substrates: A, B ; products: P, Q ; oligomeric_state: 2
        steps: begin
            E + A ⇌ EP_res_plusA_minusP                     :: EqualAI
            E_res_plusA_minusP + P <--> EP_res_plusA_minusP :: OnlyA
            E_res_plusA_minusP + B <--> EQ                  :: EqualAI
            E + Q <--> EQ                                   :: OnlyA
        end
    end
    a_p_q = @_testhelper_unguarded @enzyme_mechanism begin
        substrates: A, B ; products: P, Q ; oligomeric_state: 2
        steps: begin
            E + A ⇌ EP_res_plusA_minusP                     :: OnlyA
            E_res_plusA_minusP + P <--> EP_res_plusA_minusP :: OnlyA
            E_res_plusA_minusP + B <--> EQ                  :: EqualAI
            E + Q <--> EQ                                   :: OnlyA
        end
    end
    b_q = @_testhelper_unguarded @enzyme_mechanism begin
        substrates: A, B ; products: P, Q ; oligomeric_state: 2
        steps: begin
            E + A ⇌ EP_res_plusA_minusP                     :: EqualAI
            E_res_plusA_minusP + P <--> EP_res_plusA_minusP :: EqualAI
            E_res_plusA_minusP + B <--> EQ                  :: OnlyA
            E + Q <--> EQ                                   :: OnlyA
        end
    end
    a_b_q = @_testhelper_unguarded @enzyme_mechanism begin
        substrates: A, B ; products: P, Q ; oligomeric_state: 2
        steps: begin
            E + A ⇌ EP_res_plusA_minusP                     :: OnlyA
            E_res_plusA_minusP + P <--> EP_res_plusA_minusP :: EqualAI
            E_res_plusA_minusP + B <--> EQ                  :: OnlyA
            E + Q <--> EQ                                   :: OnlyA
        end
    end
    p_b_q = @_testhelper_unguarded @enzyme_mechanism begin
        substrates: A, B ; products: P, Q ; oligomeric_state: 2
        steps: begin
            E + A ⇌ EP_res_plusA_minusP                     :: EqualAI
            E_res_plusA_minusP + P <--> EP_res_plusA_minusP :: OnlyA
            E_res_plusA_minusP + B <--> EQ                  :: OnlyA
            E + Q <--> EQ                                   :: OnlyA
        end
    end
    a_p_b_q = @_testhelper_unguarded @enzyme_mechanism begin
        substrates: A, B ; products: P, Q ; oligomeric_state: 2
        steps: begin
            E + A ⇌ EP_res_plusA_minusP                     :: OnlyA
            E_res_plusA_minusP + P <--> EP_res_plusA_minusP :: OnlyA
            E_res_plusA_minusP + B <--> EQ                  :: OnlyA
            E + Q <--> EQ                                   :: OnlyA
        end
    end
    a_r_only_a = @_testhelper_unguarded @enzyme_mechanism begin
        substrates: A, B ; products: P, Q ; oligomeric_state: 2
        steps: begin
            E + A ⇌ EP_res_plusA_minusP                     :: OnlyA
            E_res_plusA_minusP + P <--> EP_res_plusA_minusP :: EqualAI
            E_res_plusA_minusP + B <--> EQ                  :: EqualAI
            E + Q <--> EQ                                   :: EqualAI
        end
        regulatory_site(2): R::OnlyA
    end
    a_r_only_i = @_testhelper_unguarded @enzyme_mechanism begin
        substrates: A, B ; products: P, Q ; oligomeric_state: 2
        steps: begin
            E + A ⇌ EP_res_plusA_minusP                     :: OnlyA
            E_res_plusA_minusP + P <--> EP_res_plusA_minusP :: EqualAI
            E_res_plusA_minusP + B <--> EQ                  :: EqualAI
            E + Q <--> EQ                                   :: EqualAI
        end
        regulatory_site(2): R::OnlyI
    end
    b_r_only_a = @_testhelper_unguarded @enzyme_mechanism begin
        substrates: A, B ; products: P, Q ; oligomeric_state: 2
        steps: begin
            E + A ⇌ EP_res_plusA_minusP                     :: EqualAI
            E_res_plusA_minusP + P <--> EP_res_plusA_minusP :: EqualAI
            E_res_plusA_minusP + B <--> EQ                  :: OnlyA
            E + Q <--> EQ                                   :: EqualAI
        end
        regulatory_site(2): R::OnlyA
    end
    b_r_only_i = @_testhelper_unguarded @enzyme_mechanism begin
        substrates: A, B ; products: P, Q ; oligomeric_state: 2
        steps: begin
            E + A ⇌ EP_res_plusA_minusP                     :: EqualAI
            E_res_plusA_minusP + P <--> EP_res_plusA_minusP :: EqualAI
            E_res_plusA_minusP + B <--> EQ                  :: OnlyA
            E + Q <--> EQ                                   :: EqualAI
        end
        regulatory_site(2): R::OnlyI
    end
    a_b_r_only_a = @_testhelper_unguarded @enzyme_mechanism begin
        substrates: A, B ; products: P, Q ; oligomeric_state: 2
        steps: begin
            E + A ⇌ EP_res_plusA_minusP                     :: OnlyA
            E_res_plusA_minusP + P <--> EP_res_plusA_minusP :: EqualAI
            E_res_plusA_minusP + B <--> EQ                  :: OnlyA
            E + Q <--> EQ                                   :: EqualAI
        end
        regulatory_site(2): R::OnlyA
    end
    a_b_r_only_i = @_testhelper_unguarded @enzyme_mechanism begin
        substrates: A, B ; products: P, Q ; oligomeric_state: 2
        steps: begin
            E + A ⇌ EP_res_plusA_minusP                     :: OnlyA
            E_res_plusA_minusP + P <--> EP_res_plusA_minusP :: EqualAI
            E_res_plusA_minusP + B <--> EQ                  :: OnlyA
            E + Q <--> EQ                                   :: EqualAI
        end
        regulatory_site(2): R::OnlyI
    end
    moved = _testhelper_move(EnzymeRates.Enumeration._expand_to_allosteric, parent, rxn)
    @test length(moved) == 21
    @test Set(moved) == Set([a, p, a_p, b, a_b, p_b, a_p_b, q, a_q, p_q, a_p_q, b_q, a_b_q,
                            p_b_q, a_p_b_q, a_r_only_a, a_r_only_i, b_r_only_a, b_r_only_i,
                            a_b_r_only_a, a_b_r_only_i])
    v = EnzymeRates.Derivation._guard(a)
    @test (v.fault, EnzymeRates.Mechanisms.step_text(parent.catalytic, v.edge)) ==
          (EnzymeRates.Derivation._FAST_STEP_PINNED, "E + A ⇌ EP_res_plusA_minusP")
    # Both half-reactions' fused bindings :OnlyA, beside Q's binding, P's or both, or with
    # R: both oracles agree with the guard on every face.
    tagged = [a_b_q, a_p_b, a_p_b_q, a_b_r_only_a, a_b_r_only_i]
    @test all(_testhelper_is_valid, tagged)
    @test all(t -> all(_testhelper_admission_agreement(
                           t; rng = Random.MersenneTwister(81), corners = true)), tagged)
    @test _testhelper_move(EnzymeRates.Enumeration._expand_re_to_ss, parent) == [flip]
    kids = _testhelper_children(rxn, [parent])
    @test length(kids) == 19
    @test Set(kids) == Set([flip, p, a_p, b, a_b, p_b, a_p_b, q, a_q, p_q, a_p_q, b_q,
                            a_b_q, p_b_q, a_p_b_q, b_r_only_a, b_r_only_i, a_b_r_only_a,
                            a_b_r_only_i])
    second = @enzyme_mechanism begin
        substrates: A, B
        products: P, Q
        steps: begin
            E + A <--> EP_res_plusA_minusP
            E_res_plusA_minusP + P ⇌ EP_res_plusA_minusP
            E_res_plusA_minusP + B ⇌ EQ
            E + Q <--> EQ
        end
    end
    second_flips = [
        (@enzyme_mechanism begin
            substrates: A, B
            products: P, Q
            steps: begin
                E + A <--> EP_res_plusA_minusP
                E_res_plusA_minusP + P <--> EP_res_plusA_minusP
                E_res_plusA_minusP + B ⇌ EQ
                E + Q <--> EQ
            end
        end),
        @enzyme_mechanism begin
            substrates: A, B
            products: P, Q
            steps: begin
                E + A <--> EP_res_plusA_minusP
                E_res_plusA_minusP + P ⇌ EP_res_plusA_minusP
                E_res_plusA_minusP + B <--> EQ
                E + Q <--> EQ
            end
        end]
    second_a = @_testhelper_unguarded @enzyme_mechanism begin
        substrates: A, B ; products: P, Q ; oligomeric_state: 2
        steps: begin
            E + A <--> EP_res_plusA_minusP               :: OnlyA
            E_res_plusA_minusP + P ⇌ EP_res_plusA_minusP :: EqualAI
            E_res_plusA_minusP + B ⇌ EQ                  :: EqualAI
            E + Q <--> EQ                                :: EqualAI
        end
    end
    second_p = @_testhelper_unguarded @enzyme_mechanism begin
        substrates: A, B ; products: P, Q ; oligomeric_state: 2
        steps: begin
            E + A <--> EP_res_plusA_minusP               :: EqualAI
            E_res_plusA_minusP + P ⇌ EP_res_plusA_minusP :: OnlyA
            E_res_plusA_minusP + B ⇌ EQ                  :: EqualAI
            E + Q <--> EQ                                :: EqualAI
        end
    end
    second_a_p = @_testhelper_unguarded @enzyme_mechanism begin
        substrates: A, B ; products: P, Q ; oligomeric_state: 2
        steps: begin
            E + A <--> EP_res_plusA_minusP               :: OnlyA
            E_res_plusA_minusP + P ⇌ EP_res_plusA_minusP :: OnlyA
            E_res_plusA_minusP + B ⇌ EQ                  :: EqualAI
            E + Q <--> EQ                                :: EqualAI
        end
    end
    second_b = @_testhelper_unguarded @enzyme_mechanism begin
        substrates: A, B ; products: P, Q ; oligomeric_state: 2
        steps: begin
            E + A <--> EP_res_plusA_minusP               :: EqualAI
            E_res_plusA_minusP + P ⇌ EP_res_plusA_minusP :: EqualAI
            E_res_plusA_minusP + B ⇌ EQ                  :: OnlyA
            E + Q <--> EQ                                :: EqualAI
        end
    end
    second_a_b = @_testhelper_unguarded @enzyme_mechanism begin
        substrates: A, B ; products: P, Q ; oligomeric_state: 2
        steps: begin
            E + A <--> EP_res_plusA_minusP               :: OnlyA
            E_res_plusA_minusP + P ⇌ EP_res_plusA_minusP :: EqualAI
            E_res_plusA_minusP + B ⇌ EQ                  :: OnlyA
            E + Q <--> EQ                                :: EqualAI
        end
    end
    second_p_b = @_testhelper_unguarded @enzyme_mechanism begin
        substrates: A, B ; products: P, Q ; oligomeric_state: 2
        steps: begin
            E + A <--> EP_res_plusA_minusP               :: EqualAI
            E_res_plusA_minusP + P ⇌ EP_res_plusA_minusP :: OnlyA
            E_res_plusA_minusP + B ⇌ EQ                  :: OnlyA
            E + Q <--> EQ                                :: EqualAI
        end
    end
    second_a_p_b = @_testhelper_unguarded @enzyme_mechanism begin
        substrates: A, B ; products: P, Q ; oligomeric_state: 2
        steps: begin
            E + A <--> EP_res_plusA_minusP               :: OnlyA
            E_res_plusA_minusP + P ⇌ EP_res_plusA_minusP :: OnlyA
            E_res_plusA_minusP + B ⇌ EQ                  :: OnlyA
            E + Q <--> EQ                                :: EqualAI
        end
    end
    second_q = @_testhelper_unguarded @enzyme_mechanism begin
        substrates: A, B ; products: P, Q ; oligomeric_state: 2
        steps: begin
            E + A <--> EP_res_plusA_minusP               :: EqualAI
            E_res_plusA_minusP + P ⇌ EP_res_plusA_minusP :: EqualAI
            E_res_plusA_minusP + B ⇌ EQ                  :: EqualAI
            E + Q <--> EQ                                :: OnlyA
        end
    end
    second_a_q = @_testhelper_unguarded @enzyme_mechanism begin
        substrates: A, B ; products: P, Q ; oligomeric_state: 2
        steps: begin
            E + A <--> EP_res_plusA_minusP               :: OnlyA
            E_res_plusA_minusP + P ⇌ EP_res_plusA_minusP :: EqualAI
            E_res_plusA_minusP + B ⇌ EQ                  :: EqualAI
            E + Q <--> EQ                                :: OnlyA
        end
    end
    second_p_q = @_testhelper_unguarded @enzyme_mechanism begin
        substrates: A, B ; products: P, Q ; oligomeric_state: 2
        steps: begin
            E + A <--> EP_res_plusA_minusP               :: EqualAI
            E_res_plusA_minusP + P ⇌ EP_res_plusA_minusP :: OnlyA
            E_res_plusA_minusP + B ⇌ EQ                  :: EqualAI
            E + Q <--> EQ                                :: OnlyA
        end
    end
    second_a_p_q = @_testhelper_unguarded @enzyme_mechanism begin
        substrates: A, B ; products: P, Q ; oligomeric_state: 2
        steps: begin
            E + A <--> EP_res_plusA_minusP               :: OnlyA
            E_res_plusA_minusP + P ⇌ EP_res_plusA_minusP :: OnlyA
            E_res_plusA_minusP + B ⇌ EQ                  :: EqualAI
            E + Q <--> EQ                                :: OnlyA
        end
    end
    second_b_q = @_testhelper_unguarded @enzyme_mechanism begin
        substrates: A, B ; products: P, Q ; oligomeric_state: 2
        steps: begin
            E + A <--> EP_res_plusA_minusP               :: EqualAI
            E_res_plusA_minusP + P ⇌ EP_res_plusA_minusP :: EqualAI
            E_res_plusA_minusP + B ⇌ EQ                  :: OnlyA
            E + Q <--> EQ                                :: OnlyA
        end
    end
    second_a_b_q = @_testhelper_unguarded @enzyme_mechanism begin
        substrates: A, B ; products: P, Q ; oligomeric_state: 2
        steps: begin
            E + A <--> EP_res_plusA_minusP               :: OnlyA
            E_res_plusA_minusP + P ⇌ EP_res_plusA_minusP :: EqualAI
            E_res_plusA_minusP + B ⇌ EQ                  :: OnlyA
            E + Q <--> EQ                                :: OnlyA
        end
    end
    second_p_b_q = @_testhelper_unguarded @enzyme_mechanism begin
        substrates: A, B ; products: P, Q ; oligomeric_state: 2
        steps: begin
            E + A <--> EP_res_plusA_minusP               :: EqualAI
            E_res_plusA_minusP + P ⇌ EP_res_plusA_minusP :: OnlyA
            E_res_plusA_minusP + B ⇌ EQ                  :: OnlyA
            E + Q <--> EQ                                :: OnlyA
        end
    end
    second_a_p_b_q = @_testhelper_unguarded @enzyme_mechanism begin
        substrates: A, B ; products: P, Q ; oligomeric_state: 2
        steps: begin
            E + A <--> EP_res_plusA_minusP               :: OnlyA
            E_res_plusA_minusP + P ⇌ EP_res_plusA_minusP :: OnlyA
            E_res_plusA_minusP + B ⇌ EQ                  :: OnlyA
            E + Q <--> EQ                                :: OnlyA
        end
    end
    second_a_r_only_a = @_testhelper_unguarded @enzyme_mechanism begin
        substrates: A, B ; products: P, Q ; oligomeric_state: 2
        steps: begin
            E + A <--> EP_res_plusA_minusP               :: OnlyA
            E_res_plusA_minusP + P ⇌ EP_res_plusA_minusP :: EqualAI
            E_res_plusA_minusP + B ⇌ EQ                  :: EqualAI
            E + Q <--> EQ                                :: EqualAI
        end
        regulatory_site(2): R::OnlyA
    end
    second_a_r_only_i = @_testhelper_unguarded @enzyme_mechanism begin
        substrates: A, B ; products: P, Q ; oligomeric_state: 2
        steps: begin
            E + A <--> EP_res_plusA_minusP               :: OnlyA
            E_res_plusA_minusP + P ⇌ EP_res_plusA_minusP :: EqualAI
            E_res_plusA_minusP + B ⇌ EQ                  :: EqualAI
            E + Q <--> EQ                                :: EqualAI
        end
        regulatory_site(2): R::OnlyI
    end
    second_b_r_only_a = @_testhelper_unguarded @enzyme_mechanism begin
        substrates: A, B ; products: P, Q ; oligomeric_state: 2
        steps: begin
            E + A <--> EP_res_plusA_minusP               :: EqualAI
            E_res_plusA_minusP + P ⇌ EP_res_plusA_minusP :: EqualAI
            E_res_plusA_minusP + B ⇌ EQ                  :: OnlyA
            E + Q <--> EQ                                :: EqualAI
        end
        regulatory_site(2): R::OnlyA
    end
    second_b_r_only_i = @_testhelper_unguarded @enzyme_mechanism begin
        substrates: A, B ; products: P, Q ; oligomeric_state: 2
        steps: begin
            E + A <--> EP_res_plusA_minusP               :: EqualAI
            E_res_plusA_minusP + P ⇌ EP_res_plusA_minusP :: EqualAI
            E_res_plusA_minusP + B ⇌ EQ                  :: OnlyA
            E + Q <--> EQ                                :: EqualAI
        end
        regulatory_site(2): R::OnlyI
    end
    second_a_b_r_only_a = @_testhelper_unguarded @enzyme_mechanism begin
        substrates: A, B ; products: P, Q ; oligomeric_state: 2
        steps: begin
            E + A <--> EP_res_plusA_minusP               :: OnlyA
            E_res_plusA_minusP + P ⇌ EP_res_plusA_minusP :: EqualAI
            E_res_plusA_minusP + B ⇌ EQ                  :: OnlyA
            E + Q <--> EQ                                :: EqualAI
        end
        regulatory_site(2): R::OnlyA
    end
    second_a_b_r_only_i = @_testhelper_unguarded @enzyme_mechanism begin
        substrates: A, B ; products: P, Q ; oligomeric_state: 2
        steps: begin
            E + A <--> EP_res_plusA_minusP               :: OnlyA
            E_res_plusA_minusP + P ⇌ EP_res_plusA_minusP :: EqualAI
            E_res_plusA_minusP + B ⇌ EQ                  :: OnlyA
            E + Q <--> EQ                                :: EqualAI
        end
        regulatory_site(2): R::OnlyI
    end
    second_p_b_nonequal = @enzyme_mechanism begin
        substrates: A, B ; products: P, Q ; oligomeric_state: 2
        steps: begin
            E + A <--> EP_res_plusA_minusP               :: EqualAI
            E_res_plusA_minusP + P ⇌ EP_res_plusA_minusP :: NonequalAI
            E_res_plusA_minusP + B ⇌ EQ                  :: NonequalAI
            E + Q <--> EQ                                :: EqualAI
        end
    end
    second_moved = _testhelper_move(EnzymeRates.Enumeration._expand_to_allosteric, second,
                                    rxn)
    @test length(second_moved) == 22
    @test Set(second_moved) == Set([second_a, second_p, second_a_p, second_b, second_a_b,
                                   second_p_b, second_a_p_b, second_q, second_a_q,
                                   second_p_q, second_a_p_q, second_b_q, second_a_b_q,
                                   second_p_b_q, second_a_p_b_q, second_a_r_only_a,
                                   second_a_r_only_i, second_b_r_only_a, second_b_r_only_i,
                                   second_a_b_r_only_a, second_a_b_r_only_i,
                                   second_p_b_nonequal])
    second_returned = [second_a, second_q, second_a_q, second_a_p, second_a_p_q, second_b_q,
                       second_a_b_q, second_a_p_b_q, second_a_r_only_a, second_a_r_only_i,
                       second_p_b_nonequal]
    @test all(_testhelper_is_valid, second_returned)
    v = EnzymeRates.Derivation._guard(second_a_p_b)
    @test (v.fault, EnzymeRates.Mechanisms.step_text(second.catalytic, v.edge)) ==
          (EnzymeRates.Derivation._FILLED, "E_res_plusA_minusP + B ⇌ EQ")
    @test all(_testhelper_admission_agreement(second_a_p_b;
                                              rng = Random.MersenneTwister(82),
                                              corners = true))
    v = EnzymeRates.Derivation._guard(second_p_b)
    @test (v.fault, EnzymeRates.Mechanisms.step_text(second.catalytic, v.edge)) ==
          (EnzymeRates.Derivation._BRIDGE, "E_res_plusA_minusP + B ⇌ EQ")
    second_kids = _testhelper_children(rxn, [second])
    @test length(second_kids) == 13
    @test Set(second_kids) == Set([second_flips; second_returned])
end

@testset "class 3 inside the flip move emits the superset flips" begin
    # Random-order merged bi-bi with the abortive complex EAP (census
    # 0988ab459f638946). Flipping the A group alone gains a segment but leaves a face
    # whose steady state moves with the relative speeds of the rapid-equilibrium steps, a
    # class-3 rejection (`_LIMIT_SPEED`) the RE-speed oracle reads as a pattern. The flip
    # move's guard (`check_valid_rate_equation`) fails that set, so the move extends it:
    # the A group with each of the other three rapid-equilibrium groups. The other two
    # minimal sets pass and stay: five children; the minimal sets alone, without class 3
    # in the predicate, would give two.
    parent = @enzyme_mechanism begin
        substrates: A, B
        products: P, Q
        steps: begin
            (E + A ⇌ EA, EB + A ⇌ EPQ)
            (E + B <--> EB, EA + B <--> EPQ)
            E + P ⇌ EP
            EA + P ⇌ EAP
            EP + A ⇌ EAP
            EP + Q <--> EPQ
        end
    end
    a_alone = @_testhelper_unguarded @enzyme_mechanism begin
        substrates: A, B
        products: P, Q
        steps: begin
            (E + A <--> EA, EB + A <--> EPQ)
            (E + B <--> EB, EA + B <--> EPQ)
            E + P ⇌ EP
            EA + P ⇌ EAP
            EP + A ⇌ EAP
            EP + Q <--> EPQ
        end
    end
    v = EnzymeRates.Derivation._guard(a_alone)
    @test (v.fault, v.limit, v.face) ===
          (EnzymeRates.Derivation._SPLIT, EnzymeRates.Derivation._LIMIT_SPEED, 2)
    @test EnzymeRates.Mechanisms.step_text(a_alone.catalytic, v.edge) == "E + P ⇌ EP"
    @test _testhelper_respeed_verdicts(a_alone; rng = Random.MersenneTwister(91),
                                       law = :limit)[v.face] === :pattern
    flags = _testhelper_flags(parent)
    flags[[_testhelper_edge(parent, "E + A ⇌ EA"),
           _testhelper_edge(parent, "EB + A ⇌ EPQ")]] .= false
    @test EnzymeRates.Enumeration._with_flags(_testhelper_candidate(parent),
                                              flags).mechanism ==
          a_alone
    expected = [
        (@enzyme_mechanism begin      # A and P@E
            substrates: A, B
            products: P, Q
            steps: begin
                (E + A <--> EA, EB + A <--> EPQ)
                (E + B <--> EB, EA + B <--> EPQ)
                E + P <--> EP
                EA + P ⇌ EAP
                EP + A ⇌ EAP
                EP + Q <--> EPQ
            end
        end),
        (@enzyme_mechanism begin      # A and P@EA
            substrates: A, B
            products: P, Q
            steps: begin
                (E + A <--> EA, EB + A <--> EPQ)
                (E + B <--> EB, EA + B <--> EPQ)
                E + P ⇌ EP
                EA + P <--> EAP
                EP + A ⇌ EAP
                EP + Q <--> EPQ
            end
        end),
        (@enzyme_mechanism begin      # A and A@EP
            substrates: A, B
            products: P, Q
            steps: begin
                (E + A <--> EA, EB + A <--> EPQ)
                (E + B <--> EB, EA + B <--> EPQ)
                E + P ⇌ EP
                EA + P ⇌ EAP
                EP + A <--> EAP
                EP + Q <--> EPQ
            end
        end),
        (@enzyme_mechanism begin      # P@E and P@EA
            substrates: A, B
            products: P, Q
            steps: begin
                (E + A ⇌ EA, EB + A ⇌ EPQ)
                (E + B <--> EB, EA + B <--> EPQ)
                E + P <--> EP
                EA + P <--> EAP
                EP + A ⇌ EAP
                EP + Q <--> EPQ
            end
        end),
        (@enzyme_mechanism begin      # P@E and A@EP
            substrates: A, B
            products: P, Q
            steps: begin
                (E + A ⇌ EA, EB + A ⇌ EPQ)
                (E + B <--> EB, EA + B <--> EPQ)
                E + P <--> EP
                EA + P ⇌ EAP
                EP + A <--> EAP
                EP + Q <--> EPQ
            end
        end)]
    children = _testhelper_move(EnzymeRates.Enumeration._expand_re_to_ss, parent)
    @test length(children) == 5
    @test Set(children) == Set(expected)
    @test a_alone ∉ children
end

@testset "the flip move's guard reads both conformations of an allosteric child" begin
    # A ping-pong whose abortive complex EAP forms in the active conformation only.
    # Flipping the A group leaves, at zero A and B, the rapid-equilibrium segment of
    # F = E_res_plusA_minusP and FP split between the two conformations, which only
    # B's uptake joins. The split moves with the relative speeds of the rapid-equilibrium
    # steps (`_LIMIT_SPEED`), a pattern to the RE-speed oracle. The catalytic graph alone,
    # one conformation, passes. The move's guard reads the child's two conformations, so
    # the A flip fails, and none of its supersets is minimal: the children are the four
    # other single flips.
    parent = @enzyme_mechanism begin
        substrates: A, B
        products: P, Q
        oligomeric_state: 1
        steps: begin
            (E + A ⇌ EA,
             EP_res_plusA_minusP + A ⇌ EAP_res_plusA_minusP)  :: EqualAI
            E + Q ⇌ EQ                                             :: EqualAI
            E_res_plusA_minusP + B ⇌ EB_res_plusA_minusP       :: EqualAI
            E_res_plusA_minusP + P ⇌ EP_res_plusA_minusP       :: EqualAI
            EA + P ⇌ EAP                                       :: OnlyA
            EA <--> EP_res_plusA_minusP                         :: EqualAI
            EB_res_plusA_minusP ⇌ EQ                            :: EqualAI
        end
    end
    a_flip = @_testhelper_unguarded @enzyme_mechanism begin
        substrates: A, B
        products: P, Q
        oligomeric_state: 1
        steps: begin
            (E + A <--> EA,
             EP_res_plusA_minusP + A <--> EAP_res_plusA_minusP)  :: EqualAI
            E + Q ⇌ EQ                                               :: EqualAI
            E_res_plusA_minusP + B ⇌ EB_res_plusA_minusP         :: EqualAI
            E_res_plusA_minusP + P ⇌ EP_res_plusA_minusP         :: EqualAI
            EA + P ⇌ EAP                                         :: OnlyA
            EA <--> EP_res_plusA_minusP                           :: EqualAI
            EB_res_plusA_minusP ⇌ EQ                              :: EqualAI
        end
    end
    v = EnzymeRates.Derivation._guard(a_flip)
    @test (v.fault, v.limit, v.face) ===
          (EnzymeRates.Derivation._SPLIT, EnzymeRates.Derivation._LIMIT_SPEED, 3)
    @test _testhelper_respeed_verdicts(a_flip; rng = Random.MersenneTwister(29),
                                       law = :limit)[v.face] === :pattern
    cm = parent.catalytic
    flags = _testhelper_flags(cm)
    a = cm.steps[_testhelper_edge(parent, "E + A ⇌ EA")].kinetic_group
    flags[[e.kinetic_group == a for e in cm.steps]] .= false
    @test EnzymeRates.Enumeration._with_flags(_testhelper_candidate(parent),
                                              flags).mechanism ==
          a_flip
    @test _testhelper_is_valid(
        EnzymeRates.Enumeration._with_flags(_testhelper_candidate(cm), flags).mechanism)
    expected = [
        (@enzyme_mechanism begin      # Q@E
            substrates: A, B
            products: P, Q
            oligomeric_state: 1
            steps: begin
                (E + A ⇌ EA,
                 EP_res_plusA_minusP + A ⇌ EAP_res_plusA_minusP)  :: EqualAI
                E + Q <--> EQ                                          :: EqualAI
                E_res_plusA_minusP + B ⇌ EB_res_plusA_minusP       :: EqualAI
                E_res_plusA_minusP + P ⇌ EP_res_plusA_minusP       :: EqualAI
                EA + P ⇌ EAP                                       :: OnlyA
                EA <--> EP_res_plusA_minusP                         :: EqualAI
                EB_res_plusA_minusP ⇌ EQ                            :: EqualAI
            end
        end),
        (@enzyme_mechanism begin      # B@F
            substrates: A, B
            products: P, Q
            oligomeric_state: 1
            steps: begin
                (E + A ⇌ EA,
                 EP_res_plusA_minusP + A ⇌ EAP_res_plusA_minusP)  :: EqualAI
                E + Q ⇌ EQ                                             :: EqualAI
                E_res_plusA_minusP + B <--> EB_res_plusA_minusP    :: EqualAI
                E_res_plusA_minusP + P ⇌ EP_res_plusA_minusP       :: EqualAI
                EA + P ⇌ EAP                                       :: OnlyA
                EA <--> EP_res_plusA_minusP                         :: EqualAI
                EB_res_plusA_minusP ⇌ EQ                            :: EqualAI
            end
        end),
        (@enzyme_mechanism begin      # P@F
            substrates: A, B
            products: P, Q
            oligomeric_state: 1
            steps: begin
                (E + A ⇌ EA,
                 EP_res_plusA_minusP + A ⇌ EAP_res_plusA_minusP)  :: EqualAI
                E + Q ⇌ EQ                                             :: EqualAI
                E_res_plusA_minusP + B ⇌ EB_res_plusA_minusP       :: EqualAI
                E_res_plusA_minusP + P <--> EP_res_plusA_minusP    :: EqualAI
                EA + P ⇌ EAP                                       :: OnlyA
                EA <--> EP_res_plusA_minusP                         :: EqualAI
                EB_res_plusA_minusP ⇌ EQ                            :: EqualAI
            end
        end),
        (@enzyme_mechanism begin      # the second half's chemistry
            substrates: A, B
            products: P, Q
            oligomeric_state: 1
            steps: begin
                (E + A ⇌ EA,
                 EP_res_plusA_minusP + A ⇌ EAP_res_plusA_minusP)  :: EqualAI
                E + Q ⇌ EQ                                             :: EqualAI
                E_res_plusA_minusP + B ⇌ EB_res_plusA_minusP       :: EqualAI
                E_res_plusA_minusP + P ⇌ EP_res_plusA_minusP       :: EqualAI
                EA + P ⇌ EAP                                       :: OnlyA
                EA <--> EP_res_plusA_minusP                         :: EqualAI
                EB_res_plusA_minusP <--> EQ                         :: EqualAI
            end
        end)]
    children = _testhelper_move(EnzymeRates.Enumeration._expand_re_to_ss, parent)
    @test length(children) == 4
    @test Set(children) == Set(expected)
    @test a_flip ∉ children
end

@testset "a seed that binds no inhibitor has no :OnlyI group a round later" begin
    # Aggregate check over the whole uni-uni seed set with an allosteric regulator and a
    # competitive inhibitor, and its first expansion round: no seed or child tags a
    # kinetic group :OnlyI, since a move gives the tag to a competitive inhibitor's group
    # alone and the seed binds no inhibitor (a hand-written mechanism may tag any group),
    # while the tagging move still gives regulatory sites the tag. The round's seven
    # allosteric children are the five sets of the seed's groups :OnlyA that the exits
    # keep and the chemistry alone :OnlyA with R at a site, tagged :OnlyA and :OnlyI.
    rxn = @enzyme_reaction begin
        substrates: S[C]
        products: P[C]
        allosteric_regulators: R
        competitive_inhibitors: I
        oligomeric_state: 2
    end
    seeds = _testhelper_mechanisms(EnzymeRates.Enumeration.init_mechanisms(rxn))
    children = _testhelper_children(rxn, seeds)
    allo = [m for m in [seeds; children] if _testhelper_is_allosteric(m)]
    @test length(seeds) == 1
    @test length(allo) == 7
    @test all(m -> :OnlyI ∉ m.allosteric.allosteric_states, allo)
    @test any(m -> any(s -> :OnlyI in s.allosteric_states,
                       m.allosteric.regulatory_sites), allo)
end

@testset "the flip move keeps the flips that leave a tagged step between two groups" begin
    # The guard admits the seed. Its chemistry ES ⇌ EP sits in rapid equilibrium with
    # both sides, so the seed is degenerate and expanded. Six flip sets gain a segment and
    # pass the guard. Three of them (`between`: the S binding, the chemistry, or the pair
    # EP ⇌ FP and F + P ⇌ EP made steady state) leave the :OnlyA step E <--> F
    # between two groups of forms of the inactive conformation; the guard finds a reading
    # of it that gives a limit, and both oracles agree on every face. Two of the other
    # three keep the chemistry in rapid equilibrium with both sides, which makes them
    # degenerate. The base tier keeps the four others.
    ER = EnzymeRates
    rxn = @enzyme_reaction begin
        substrates: S[C]
        products: P[C]
    end
    seed = @enzyme_mechanism begin
        substrates: S
        products: P
        steps: begin
            E + S ⇌ ES        :: OnlyA
            ES ⇌ EP         :: EqualAI
            F + P ⇌ EP        :: EqualAI
            F <--> E            :: OnlyA
            F + S ⇌ FS        :: EqualAI
            FS <--> FP      :: EqualAI
            F + P ⇌ FP        :: EqualAI
            FP ⇌ EP         :: EqualAI
        end
    end
    # Untyped: inside `T[...]` the parser reads `begin` as an index keyword.
    between = [
        @_testhelper_unguarded(@enzyme_mechanism begin
            substrates: S
            products: P
            steps: begin
                E + S <--> ES     :: OnlyA
                ES ⇌ EP         :: EqualAI
                F + P ⇌ EP        :: EqualAI
                F <--> E            :: OnlyA
                F + S ⇌ FS        :: EqualAI
                FS <--> FP      :: EqualAI
                F + P ⇌ FP        :: EqualAI
                FP ⇌ EP         :: EqualAI
            end
        end),
        @_testhelper_unguarded(@enzyme_mechanism begin
            substrates: S
            products: P
            steps: begin
                E + S ⇌ ES        :: OnlyA
                ES <--> EP      :: EqualAI
                F + P ⇌ EP        :: EqualAI
                F <--> E            :: OnlyA
                F + S ⇌ FS        :: EqualAI
                FS <--> FP      :: EqualAI
                F + P ⇌ FP        :: EqualAI
                FP ⇌ EP         :: EqualAI
            end
        end),
        @_testhelper_unguarded(@enzyme_mechanism begin
            substrates: S
            products: P
            steps: begin
                E + S ⇌ ES        :: OnlyA
                ES ⇌ EP         :: EqualAI
                F + P <--> EP     :: EqualAI
                F <--> E            :: OnlyA
                F + S ⇌ FS        :: EqualAI
                FS <--> FP      :: EqualAI
                F + P ⇌ FP        :: EqualAI
                FP <--> EP      :: EqualAI
            end
        end)]
    kept = @enzyme_mechanism begin
        substrates: S
        products: P
        steps: begin
            E + S ⇌ ES        :: OnlyA
            ES ⇌ EP         :: EqualAI
            F + P <--> EP     :: EqualAI
            F <--> E            :: OnlyA
            F + S ⇌ FS        :: EqualAI
            FS <--> FP      :: EqualAI
            F + P <--> FP     :: EqualAI
            FP ⇌ EP         :: EqualAI
        end
    end
    expected = [
        (@enzyme_mechanism begin    # F + S
            substrates: S
            products: P
            steps: begin
                E + S ⇌ ES        :: OnlyA
                ES ⇌ EP         :: EqualAI
                F + P ⇌ EP        :: EqualAI
                F <--> E            :: OnlyA
                F + S <--> FS     :: EqualAI
                FS <--> FP      :: EqualAI
                F + P ⇌ FP        :: EqualAI
                FP ⇌ EP         :: EqualAI
            end
        end),
        (@enzyme_mechanism begin    # F + P ⇌ FP, FP ⇌ EP
            substrates: S
            products: P
            steps: begin
                E + S ⇌ ES        :: OnlyA
                ES ⇌ EP         :: EqualAI
                F + P ⇌ EP        :: EqualAI
                F <--> E            :: OnlyA
                F + S ⇌ FS        :: EqualAI
                FS <--> FP      :: EqualAI
                F + P <--> FP     :: EqualAI
                FP <--> EP      :: EqualAI
            end
        end),
        kept]
    @test _testhelper_is_valid(seed)
    @test ER.Enumeration._degenerate(_testhelper_candidate(seed))
    kids = _testhelper_move(ER.Enumeration._expand_re_to_ss, seed)
    @test length(kids) == 6
    @test Set(kids) == Set([expected; between])
    for (i, k) in enumerate(between)
        @test all(_testhelper_admission_agreement(
            k; rng = Random.MersenneTwister(16000 + i), corners = true))
        _testhelper_check_law_against_limit(k)
    end
    @test !any(ER.Enumeration._degenerate, _testhelper_candidate.(between))
    @test count(ER.Enumeration._degenerate, _testhelper_candidate.(kids)) == 2
    base = _testhelper_mechanisms(
        ER.Enumeration._base_tier(_testhelper_candidate.([seed]), rxn))
    @test length(base) == 4
    @test Set(base) == Set([kept; between])
end

@testset "the flip search extends a flip set whose child the guard rejects" begin
    # An allosteric ter-bi ping-pong among the census ter-bi children: B and D bind in
    # either order, EBD turns at rapid equilibrium into the covalent form E(Q; residual
    # = B + D - Q), which releases Q, and the covalent form E_res_plusB_plusD_minusQ takes
    # up A and gives off P through the steady-state EA_res_plusB_plusD_minusQ <--> EP.
    # B's group, P's binding and both chemistry steps are :OnlyA. Every step but the
    # second chemistry step is at rapid equilibrium, so the seed is degenerate and
    # expanded. Flipping D's group alone, or Q's group alone, gains a segment, but at zero
    # P and Q the inactive conformation then feeds the covalent forms, which leave only
    # through the vanishing chemistry step (`_PUMPED`). That is admission's face
    # condition, part of the guard the flip search asks of each candidate, so the two sets
    # count as failed and are extended: D's and Q's groups together pass. A's and P's
    # bindings are the flanks of the steady-state EA_res_plusB_plusD_minusQ <--> EP and
    # never flip (`_chain_flank_groups`). The flips of B's group and the first chemistry
    # step pass alone: three children. One keeps the seed's degeneracy; the base tier
    # keeps the other two. The extension rests on Q's flip being refused: under the
    # guard's first built reading the ε-path oracle too finds no limit at zero P and Q. A
    # guard that admitted Q's flip would emit it alone and leave D's flip with no failed
    # partner to extend with.
    rxn = @enzyme_reaction begin
        substrates: A[C], B[N], D[X]
        products: P[CN], Q[X]
    end
    seed = @enzyme_mechanism begin
        substrates: A, B, D ; products: P, Q ; oligomeric_state: 1
        steps: begin
            (E + B ⇌ EB, ED + B ⇌ EBD, EDQ + B ⇌ EBDQ,
             EQ_res_plusB_plusD_minusQ + B ⇌ EBQ_res_plusB_plusD_minusQ)  :: OnlyA
            (E + D ⇌ ED, EB + D ⇌ EBD, EBQ + D ⇌ EBDQ,
             EQ_res_plusB_plusD_minusQ + D ⇌ EDQ_res_plusB_plusD_minusQ)  :: EqualAI
            E + P ⇌ EP                                                      :: OnlyA
            E_res_plusB_plusD_minusQ + A ⇌ EA_res_plusB_plusD_minusQ        :: EqualAI
            (E_res_plusB_plusD_minusQ + Q ⇌ EQ_res_plusB_plusD_minusQ,
             EB + Q ⇌ EBQ, EBD + Q ⇌ EBDQ,
             ED + Q ⇌ EDQ)                                              :: EqualAI
            EA_res_plusB_plusD_minusQ <--> EP                              :: OnlyA
            EBD ⇌ EQ_res_plusB_plusD_minusQ                              :: OnlyA
        end
    end
    b_flip = @enzyme_mechanism begin
        substrates: A, B, D ; products: P, Q ; oligomeric_state: 1
        steps: begin
            (E + B <--> EB, ED + B <--> EBD,
             EDQ + B <--> EBDQ,
             EQ_res_plusB_plusD_minusQ + B <--> EBQ_res_plusB_plusD_minusQ):: OnlyA
            (E + D ⇌ ED, EB + D ⇌ EBD, EBQ + D ⇌ EBDQ,
             EQ_res_plusB_plusD_minusQ + D ⇌ EDQ_res_plusB_plusD_minusQ)  :: EqualAI
            E + P ⇌ EP                                                      :: OnlyA
            E_res_plusB_plusD_minusQ + A ⇌ EA_res_plusB_plusD_minusQ        :: EqualAI
            (E_res_plusB_plusD_minusQ + Q ⇌ EQ_res_plusB_plusD_minusQ,
             EB + Q ⇌ EBQ, EBD + Q ⇌ EBDQ,
             ED + Q ⇌ EDQ)                                              :: EqualAI
            EA_res_plusB_plusD_minusQ <--> EP                              :: OnlyA
            EBD ⇌ EQ_res_plusB_plusD_minusQ                              :: OnlyA
        end
    end
    chemistry_flip = @enzyme_mechanism begin
        substrates: A, B, D ; products: P, Q ; oligomeric_state: 1
        steps: begin
            (E + B ⇌ EB, ED + B ⇌ EBD, EDQ + B ⇌ EBDQ,
             EQ_res_plusB_plusD_minusQ + B ⇌ EBQ_res_plusB_plusD_minusQ)  :: OnlyA
            (E + D ⇌ ED, EB + D ⇌ EBD, EBQ + D ⇌ EBDQ,
             EQ_res_plusB_plusD_minusQ + D ⇌ EDQ_res_plusB_plusD_minusQ)  :: EqualAI
            E + P ⇌ EP                                                      :: OnlyA
            E_res_plusB_plusD_minusQ + A ⇌ EA_res_plusB_plusD_minusQ        :: EqualAI
            (E_res_plusB_plusD_minusQ + Q ⇌ EQ_res_plusB_plusD_minusQ,
             EB + Q ⇌ EBQ, EBD + Q ⇌ EBDQ,
             ED + Q ⇌ EDQ)                                              :: EqualAI
            EA_res_plusB_plusD_minusQ <--> EP                              :: OnlyA
            EBD <--> EQ_res_plusB_plusD_minusQ                           :: OnlyA
        end
    end
    dq_flip = @enzyme_mechanism begin
        substrates: A, B, D ; products: P, Q ; oligomeric_state: 1
        steps: begin
            (E + B ⇌ EB, ED + B ⇌ EBD, EDQ + B ⇌ EBDQ,
             EQ_res_plusB_plusD_minusQ + B ⇌ EBQ_res_plusB_plusD_minusQ)  :: OnlyA
            (E + D <--> ED, EB + D <--> EBD,
             EBQ + D <--> EBDQ,
             EQ_res_plusB_plusD_minusQ + D <--> EDQ_res_plusB_plusD_minusQ):: EqualAI
            E + P ⇌ EP                                                      :: OnlyA
            E_res_plusB_plusD_minusQ + A ⇌ EA_res_plusB_plusD_minusQ        :: EqualAI
            (E_res_plusB_plusD_minusQ + Q <--> EQ_res_plusB_plusD_minusQ,
             EB + Q <--> EBQ, EBD + Q <--> EBDQ,
             ED + Q <--> EDQ)                                           :: EqualAI
            EA_res_plusB_plusD_minusQ <--> EP                              :: OnlyA
            EBD ⇌ EQ_res_plusB_plusD_minusQ                              :: OnlyA
        end
    end
    d_flip = @_testhelper_unguarded @enzyme_mechanism begin
        substrates: A, B, D ; products: P, Q ; oligomeric_state: 1
        steps: begin
            (E + B ⇌ EB, ED + B ⇌ EBD, EDQ + B ⇌ EBDQ,
             EQ_res_plusB_plusD_minusQ + B ⇌ EBQ_res_plusB_plusD_minusQ)  :: OnlyA
            (E + D <--> ED, EB + D <--> EBD,
             EBQ + D <--> EBDQ,
             EQ_res_plusB_plusD_minusQ + D <--> EDQ_res_plusB_plusD_minusQ):: EqualAI
            E + P ⇌ EP                                                      :: OnlyA
            E_res_plusB_plusD_minusQ + A ⇌ EA_res_plusB_plusD_minusQ        :: EqualAI
            (E_res_plusB_plusD_minusQ + Q ⇌ EQ_res_plusB_plusD_minusQ,
             EB + Q ⇌ EBQ, EBD + Q ⇌ EBDQ,
             ED + Q ⇌ EDQ)                                              :: EqualAI
            EA_res_plusB_plusD_minusQ <--> EP                              :: OnlyA
            EBD ⇌ EQ_res_plusB_plusD_minusQ                              :: OnlyA
        end
    end
    q_flip = @_testhelper_unguarded @enzyme_mechanism begin
        substrates: A, B, D ; products: P, Q ; oligomeric_state: 1
        steps: begin
            (E + B ⇌ EB, ED + B ⇌ EBD, EDQ + B ⇌ EBDQ,
             EQ_res_plusB_plusD_minusQ + B ⇌ EBQ_res_plusB_plusD_minusQ)  :: OnlyA
            (E + D ⇌ ED, EB + D ⇌ EBD, EBQ + D ⇌ EBDQ,
             EQ_res_plusB_plusD_minusQ + D ⇌ EDQ_res_plusB_plusD_minusQ)  :: EqualAI
            E + P ⇌ EP                                                      :: OnlyA
            E_res_plusB_plusD_minusQ + A ⇌ EA_res_plusB_plusD_minusQ        :: EqualAI
            (E_res_plusB_plusD_minusQ + Q <--> EQ_res_plusB_plusD_minusQ,
             EB + Q <--> EBQ, EBD + Q <--> EBDQ,
             ED + Q <--> EDQ)                                           :: EqualAI
            EA_res_plusB_plusD_minusQ <--> EP                              :: OnlyA
            EBD ⇌ EQ_res_plusB_plusD_minusQ                              :: OnlyA
        end
    end

    D = EnzymeRates.Derivation
    @test _testhelper_is_valid(seed)
    @test EnzymeRates.Enumeration._degenerate(_testhelper_candidate(seed))
    for c in (d_flip, q_flip)
        v = D._guard(c)
        @test (v.fault, D._face_names(c.catalytic, D._faces(c.catalytic)[v.face])) ==
              (D._PUMPED, [:P, :Q])
    end
    # The message names the feed whose reverse is dead, the release of Q.
    @test D._guard_message(d_flip) ==
          "Mechanism: at zero P and Q, E_res_plusB_plusD_minusQ + Q ⇌ " *
          "EQ_res_plusB_plusD_minusQ feeds the forms EA_res_plusB_plusD_minusQ, " *
          "E_res_plusB_plusD_minusQ of the inactive conformation, which leave only " *
          "through the vanishing step EA_res_plusB_plusD_minusQ <--> EP while " *
          "nothing runs back the way it came, so admission cannot show that their " *
          "weight does not depend on how fast that step vanishes."
    kids = _testhelper_move(EnzymeRates.Enumeration._expand_re_to_ss, seed)
    @test length(kids) == 3
    @test Set(kids) == Set([b_flip, chemistry_flip, dq_flip])
    @test Set(filter(m -> EnzymeRates.Enumeration._degenerate(_testhelper_candidate(m)),
                     kids)) ==
          Set([b_flip])
    base = _testhelper_mechanisms(
        EnzymeRates.Enumeration._base_tier(_testhelper_candidate.([seed]), rxn))
    @test length(base) == 2
    @test Set(base) == Set([chemistry_flip, dq_flip])
end
