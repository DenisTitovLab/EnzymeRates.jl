# ABOUTME: Tests for mechanism enumeration pipeline
# ABOUTME: Unit tests per move + integration tests per reaction

# Fixture rules for this file are in CLAUDE.md, "Enumeration-engine tests":
# every mechanism inline via the macro, every move test asserts the exact
# child set.

# The names a form's record (`Enumeration._FormRecord`) binds, productively or as copies.
_testhelper_bound(f) = [f.productive; f.copies]

# Orients the isomerization `s` of the written graph `t` (a topology) substrate-rich side
# → product-rich side, counting the names in `subs` as substrates (ties by content), and
# returns the records of its (source, destination) forms. A topology step keeps the
# direction it was written in; the mechanism constructors orient isomerizations
# (`_orientation`).
function _testhelper_iso_orient(t, s::EnzymeRates.Mechanisms.Step, subs)
    from, to = t.forms[s.from], t.forms[s.to]
    nf, nt = (count(in(subs), _testhelper_bound(f)) for f in (from, to))
    key(f) = (f.label, _testhelper_bound(f), f.adds, f.subtracts)
    nf > nt || (nf == nt && !isless(key(to), key(from))) ? (from, to) : (to, from)
end

# Whether the record of a form shows a covalent residual.
_testhelper_has_residual(f::EnzymeRates.Enumeration._FormRecord) =
    !isempty(f.adds) || !isempty(f.subtracts)

# Build a Mechanism over `rxn`'s names from a topology, a written graph (`_trimmed`) whose
# steps index the reactants in the reaction's order, each step its own kinetic group.
function _testhelper_topo_mech(rxn, t)
    subs, prods = [s.name for s in rxn.substrates], [p.name for p in rxn.products]
    EnzymeRates.Enumeration._candidate(
        [subs; prods], t.forms, collect(eachindex(subs)),
        collect(length(subs) .+ eachindex(prods)),
        [EnzymeRates.Mechanisms.Step(e.from, e.to, e.taken_up, e.given_off,
                                     e.rapid_equilibrium, k)
         for (k, e) in enumerate(t.steps)]).mechanism
end

# The names of the forms of the written graph `t` of a topology or a seed, over `rxn`'s
# names.
_testhelper_form_names(rxn, t) = Set(_testhelper_topo_mech(rxn, t).catalytic.enzyme_forms)

# The names of the forms of graph `m` where `name` (its modifier copy when `modifier`
# holds) is free on some edge (`_forms_where_free`).
_testhelper_free_at(m, name, modifier = false) =
    Set(String(m.catalytic.enzyme_forms[f])
        for f in findall(EnzymeRates.Enumeration._forms_where_free(_testhelper_candidate(m),
                                                                   name, modifier)))

# A mechanism as a written graph of the records of its forms (`_testhelper_candidate`) and
# its steps in their stored direction, each in group 0: the shape a topology of
# `_catalytic_topologies` takes, every binding taking its metabolite up. Its steps index
# the metabolites in stored order, the reaction's order for a reaction whose substrates
# and products are each written in name order.
_testhelper_flat_topo(m) =
    (; forms = _testhelper_candidate(m).forms,
     steps = [EnzymeRates.Mechanisms.Step(e.from, e.to, e.taken_up, e.given_off,
                                          e.rapid_equilibrium, 0)
              for e in m.catalytic.steps])

# Connectivity invariant: two enzyme forms identical in conformation+residual
# whose bound-metabolite sets differ by exactly one metabolite MUST be joined by
# a binding step. Returns the list of (formA, formB) pairs that violate it among the
# forms (their records) and steps of `t`: a topology or a seed (a written graph), or a
# mechanism.
_testhelper_connectivity_violations(m::EnzymeRates.Mechanism) =
    _testhelper_connectivity_violations(
        (; forms = _testhelper_candidate(m).forms, steps = m.catalytic.steps))
function _testhelper_connectivity_violations(t)
    _bset(f) = Set(_testhelper_bound(f))
    _testhelper_key(f) = (f.label, f.adds, f.subtracts, Tuple(sort(collect(_bset(f)))))
    forms = Dict{Any,Any}()
    edges = Set{Tuple{Any,Any}}()
    for s in t.steps
        a, b = t.forms[s.from], t.forms[s.to]
        forms[_testhelper_key(a)] = a; forms[_testhelper_key(b)] = b
        push!(edges, (_testhelper_key(a), _testhelper_key(b)))
        push!(edges, (_testhelper_key(b), _testhelper_key(a)))
    end
    viol = Tuple{EnzymeRates.Enumeration._FormRecord, EnzymeRates.Enumeration._FormRecord}[]
    fv = collect(values(forms))
    for s1 in fv, s2 in fv
        (s1.label == s2.label && (s1.adds, s1.subtracts) == (s2.adds, s2.subtracts)) ||
            continue
        b1, b2 = _bset(s1), _bset(s2)
        if length(b2) == length(b1) + 1 && issubset(b1, b2) &&
           !((_testhelper_key(s1), _testhelper_key(s2)) in edges)
            push!(viol, (s1, s2))
        end
    end
    viol
end

"""Numeric degeneracy probe (test oracle only): true when the rate at zero products does not
need some substrate (ratio of the rate at 10⁻⁸ of it to the rate at 1 stays above 10⁻³), or
grows without bound as the substrates scale (10⁹ against 10⁶), or the mirror for
products."""
function _testhelper_degenerate(m; ndraws = 2)
    fp = collect(parameters(m))
    g = _testhelper_graph(m)
    subs, prods = g.substrates, g.products
    mets = sort(collect(EnzymeRates.metabolites(m)))
    rate = compile_rate_equation(m)
    rng = MersenneTwister(7)
    for _ in 1:ndraws
        θ = exp.(randn(rng, length(fp)))
        p = NamedTuple{Tuple(fp)}(Tuple(θ))
        keq = exp(randn(rng))
        _testhelper_v(d) =
            rate(NamedTuple{Tuple(mets)}(Tuple(get(d, x, 0.0) for x in mets)), p, keq)
        for side in (subs, prods)
            base = _testhelper_v(Dict(x => 1.0 for x in side))
            for x in side
                r = _testhelper_v(Dict(y => (y == x ? 1e-8 : 1.0) for y in side)) / base
                (isfinite(r) && abs(r) < 1e-3) || return true
            end
            big = _testhelper_v(Dict(x => 1e9 for x in side)) /
                  _testhelper_v(Dict(x => 1e6 for x in side))
            big > 10 && return true
        end
    end
    false
end

# The number of rapid-equilibrium segments of the catalytic graph of `m`.
_testhelper_segment_count(m) =
    (cm = m.catalytic;
     EnzymeRates.Enumeration._re_segment_count(
         cm, _testhelper_flags(cm)))

# Ground truth for the parameter counts the enumeration moves report: the number of
# fitted parameters of the compiled rate equation.
_testhelper_fitted(x) = length(parameters(x))

# Sorted fitted-parameter count of each child minus its parent's.
function _testhelper_param_deltas(parent, kids)
    base = _testhelper_fitted(parent)
    sort([_testhelper_fitted(k) - base for k in kids])
end

# The names that the forms of the enumerated mechanism `m` bind as modifier copies, a pure
# inhibitor among them, read off the records of its forms (`_testhelper_candidate`).
_testhelper_modifier_names(m) =
    Set(x for f in _testhelper_candidate(m).forms for x in f.copies)

# Whether step `k` of candidate `c` binds a modifier copy: a binding into a form whose
# record holds one copy more than its source's.
_testhelper_binds_copy(c, k) =
    (e = c.mechanism.catalytic.steps[k];
     EnzymeRates.Enumeration._is_binding(e) &&
         length(c.forms[e.to].copies) == length(c.forms[e.from].copies) + 1)

# The kinetic groups of graph `m` whose bindings take up `x` as a modifier copy.
_testhelper_groups_binding_copy(m, x) =
    (c = _testhelper_candidate(m);
     unique(e.kinetic_group for (k, e) in enumerate(m.catalytic.steps)
            if _testhelper_binds_copy(c, k) &&
               m.catalytic.metabolites[only(e.taken_up)] === x))

# Whether step `e` of candidate `c` is a fused binding: it takes up one metabolite and
# gives off none, but its target's record is not its source's with that metabolite bound
# and the residual unchanged.
function _testhelper_is_fused(c, e)
    EnzymeRates.Enumeration._is_binding(e) || return false
    f, t = c.forms[e.from], c.forms[e.to]
    x = c.mechanism.catalytic.metabolites[only(e.taken_up)]
    !((f.adds, f.subtracts) == (t.adds, t.subtracts) &&
      sort!([_testhelper_bound(f); x]) == sort!(_testhelper_bound(t)))
end

# Whether `child` keeps the reaction of `parent`: the same metabolites as substrates and
# products, and every name `parent` binds as a modifier copy still bound so.
_testhelper_keeps_reaction(child, parent) =
    child.catalytic.metabolites[child.catalytic.substrates] ==
        parent.catalytic.metabolites[parent.catalytic.substrates] &&
    child.catalytic.metabolites[child.catalytic.products] ==
        parent.catalytic.metabolites[parent.catalytic.products] &&
    _testhelper_modifier_names(parent) ⊆ _testhelper_modifier_names(child)

# Whether `m` holds an isomerization step. It reads the graph's steps, not
# `_testhelper_graph`: the population tests call it on every child of two expansion
# levels, where building the test reader each time would dominate their run time.
_testhelper_holds_iso(m) =
    any(EnzymeRates.Mechanisms._is_iso, m.catalytic.steps)

# The rapid-equilibrium flags of each kinetic group's edges of `m`, group by group.
_testhelper_group_flags(m) = (cm = m.catalytic;
                              [[cm.steps[k].rapid_equilibrium for k in ks]
                               for ks in _testhelper_group_members(cm)])

# The steps of mechanism or catalytic scheme `m` apart from their kinetic groups: what a
# regrouping move keeps.
_testhelper_ungrouped(m::EnzymeRates.Mechanism) = _testhelper_ungrouped(m.catalytic)
_testhelper_ungrouped(c::EnzymeRates.Mechanisms.CatalyticScheme) =
    [(e.from, e.to, e.taken_up, e.given_off, e.rapid_equilibrium) for e in c.steps]

# The kinetic groups of `m` whose steps are isomerizations.
_testhelper_iso_groups(m) =
    unique(s.group for s in _testhelper_graph(m).steps if _testhelper_is_iso(s))

# Whether a step of graph `m` binds a metabolite named in `xs`.
_testhelper_binds(m, xs...) = e -> EnzymeRates.Enumeration._is_binding(e) &&
    m.catalytic.metabolites[only(e.taken_up)] in xs

# The kinetic-group index of every step of `m` that binds the metabolite named `x`.
_testhelper_groups_binding(m, x) = [s.group for s in _testhelper_graph(m).steps
                                    if _testhelper_is_binding(s) && only(s.consumed) == x]

# Whether a step is a Theorell–Chance step: it takes up one metabolite and gives off
# another.
_testhelper_is_tc(e) = !isempty(e.taken_up) && !isempty(e.given_off)

# The flags of graph `m` with every kinetic group that holds an edge satisfying one of
# `preds` turned steady state.
function _testhelper_flip_flags(m, preds...)
    re = _testhelper_flags(m)
    hit = EnzymeRates.Enumeration._groups_with(
        k -> any(p -> p(m.catalytic.steps[k]), preds), m.catalytic)
    re[hit[[e.kinetic_group for e in m.catalytic.steps]]] .= false
    re
end

# Graph `m` with every kinetic group that holds an edge satisfying one of `preds` turned
# steady state (`_testhelper_flip_flags`), built through the checked constructor.
_testhelper_flip_matching(m, preds...) =
    EnzymeRates.Enumeration._with_flags(_testhelper_candidate(m),
                                        _testhelper_flip_flags(m, preds...)).mechanism

# The graph that `steps`, between the forms of the seed `m`, write, built through the
# checked constructor as a seed variant builds its base: the forms no step joins left out
# (`_trimmed`, on the records of `m`'s forms), the group labels numbered
# (`_numbered_groups`).
function _testhelper_base(m, steps)
    g = EnzymeRates.Enumeration._trimmed(_testhelper_candidate(m).forms,
                                         EnzymeRates.Enumeration._numbered_groups(steps))
    cm = m.catalytic
    EnzymeRates.Enumeration._candidate(cm.metabolites, g.forms, cm.substrates, cm.products,
                                       g.steps).mechanism
end

# The bases of the seed `m` (`_seed_bases`), each built through the checked constructor as
# a seed variant builds it.
_testhelper_bases(m) =
    [EnzymeRates.Enumeration._candidate(m.catalytic.metabolites, b.forms,
                                        m.catalytic.substrates, m.catalytic.products,
                                        b.steps).mechanism
     for b in EnzymeRates.Enumeration._seed_bases(_testhelper_candidate(m))]

# The productive twins of graph `m` over its whole graph (`_productive_twin`).
_testhelper_twin(m) = EnzymeRates.Enumeration._productive_twin(
    _testhelper_candidate(m), trues(length(m.catalytic.steps)))

# A kinetic group is catalytic iff it holds a chemistry step (anything but a plain
# binding: an isomerization, a fused binding or a Theorell–Chance step); otherwise it is a
# binding group. Mirrors the rule `_expand_to_allosteric` uses to decide whether a group's
# `:OnlyA` flip needs a paired regulator to be distinguishable.
_testhelper_is_catalytic_group(m, g) =
    (c = _testhelper_candidate(m);
     any(e -> e.kinetic_group == g &&
              (!EnzymeRates.Enumeration._is_binding(e) || _testhelper_is_fused(c, e)),
         m.catalytic.steps))

const uni_uni_rxn = @enzyme_reaction begin
    substrates: S[C]
    products: P[C]
end

const uni_bi_rxn = @enzyme_reaction begin
    substrates: S[AB]
    products: P[A], Q[B]
end

const bi_bi_rxn = @enzyme_reaction begin
    substrates: A[C], B[N]
    products: P[C], Q[N]
end

const bi_bi_pp_rxn = @enzyme_reaction begin
    substrates: A[CX], B[N]
    products: P[C], Q[NX]
end

const uni_uni_with_reg = @enzyme_reaction begin
    substrates: S[C]
    products: P[C]
    competitive_inhibitors: I
end

const uni_uni_allo = @enzyme_reaction begin
    substrates: S[C]
    products: P[C]
    oligomeric_state: 2
end

const uni_uni_allo_reg = @enzyme_reaction begin
    substrates: S[C]
    products: P[C]
    allosteric_regulators: R
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
# regulatory site (and a V-type), R1 to a dead-end binding.
const uni_uni_reg_and_inhibitor = @enzyme_reaction begin
    substrates: S[C]
    products: P[C]
    allosteric_regulators: R2
    competitive_inhibitors: R1
    oligomeric_state: 2
end

const ter_ter_rxn = @enzyme_reaction begin
    substrates: A[C], B[N], D[X]
    products: P[C], Q[N], R[X]
end

const ter_bi_rxn = @enzyme_reaction begin
    substrates: A[C], B[N], D[X]
    products: P[CN], Q[X]
end

# Pyruvate carboxylase: Pyr + HCO3 + ATP = OAA + ADP + Pi
# Mechanism: ATP+HCO3 → ADP+Pi+CO2_residual,
#            then Pyr+CO2 → OAA
const pyruvate_carboxylase_rxn = @enzyme_reaction begin
    substrates: Pyr[C3H3O3], HCO3[HCO3], ATP[C10H16N5O13P3]
    products: OAA[C4H3O5], ADP[C10H15N5O10P2], Pi[H2PO4]
end

# Pyruvate dehydrogenase: Pyr + NAD + CoA = AcCoA + NADH + CO2
# Mechanism: Pyr → CO2+residual, CoA+residual → AcCoA+residual,
#            NAD+residual → NADH
const pyruvate_dehydrogenase_rxn = @enzyme_reaction begin
    substrates: Pyr[C3H3O3], NAD[C21H28N7O14P2], CoA[C21H36N7O16P3S]
    products: AcCoA[C23H38N7O17P3S], NADH[C21H29N7O14P2], CO2[CO2]
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

@testset "_candidate merges a record given twice and refuses two records of one name" begin
    R, Step = EnzymeRates.Enumeration._FormRecord, EnzymeRates.Mechanisms.Step
    none = Symbol[]
    E, ES, EP = R(:E, none, none, none, none), R(:E, [:S], none, none, none),
                R(:E, [:P], none, none, none)
    steps = [Step(1, 2, [1], Int[], true, 1), Step(2, 3, Int[], Int[], false, 2),
             Step(4, 3, [2], Int[], true, 3)]
    # E is given twice, as forms 1 and 4: one form.
    c = EnzymeRates.Enumeration._candidate([:S, :P], [E, ES, EP, E], [1], [2], steps)
    @test c.mechanism.catalytic.enzyme_forms == [:E, :EP, :ES]
    @test c.forms == [E, EP, ES]
    # The conformation Estar with nothing bound and E with a metabolite star bound both
    # take the name Estar.
    Estar, Ewithstar = R(:Estar, none, none, none, none), R(:E, none, [:star], none, none)
    @test_throws(ErrorException("_candidate: two different form records are both named " *
                                "Estar"),
                 EnzymeRates.Enumeration._candidate(
                     [:S, :P, :star], [E, ES, EP, Estar, Ewithstar], [1], [2],
                     [steps; Step(1, 4, Int[], Int[], false, 4);
                      Step(1, 5, [3], Int[], true, 5)]))
end

@testset "the test reader's records of a seed equal the enumerator's" begin
    # An aggregate check over whole seed sets (CLAUDE.md enumeration rule 1's exception):
    # the records `_testhelper_candidate` reads from each seed's form names equal the
    # records the enumerator built the seed from, the copies of A and the inhibitor I
    # among them.
    rxn = @enzyme_reaction begin
        substrates: A[C], B[N]
        products: P[C], Q[N]
        competitive_inhibitors: A, I
    end
    seeds = [EnzymeRates.Enumeration.init_mechanisms(rxn);
             EnzymeRates.Enumeration.seed_mechanisms(rxn, Set{Symbol}(), Set([:A, :I]))]
    @test any(c -> any(f -> :A in f.copies, c.forms), seeds) &&
          any(c -> any(f -> :I in f.copies, c.forms), seeds)
    @test all(c -> _testhelper_candidate(c.mechanism).forms == c.forms, seeds)
end

@testset "Canonical-by-construction representation independence" begin
    # One kinetic group holding two bindings of the same metabolite S (context-
    # shared: S binds free E, and S binds A-bound E, with the same K), A's binding, the
    # chemistry and P's release. Written two ways: the statements in reverse order AND
    # the group's two steps swapped. Canonical-by-construction must collapse them to one
    # struct.
    m_orderA = @enzyme_mechanism begin
        substrates: S, A
        products: P
        steps: begin
            (E + S ⇌ ES, EA + S ⇌ EAS)
            E + A ⇌ EA
            EAS <--> EP
            EP ⇌ E + P
        end
    end
    m_orderB = @enzyme_mechanism begin
        substrates: S, A
        products: P
        steps: begin
            EP ⇌ E + P
            EAS <--> EP
            E + A ⇌ EA
            (EA + S ⇌ EAS, E + S ⇌ ES)
        end
    end
    @test m_orderA == m_orderB
    @test hash(m_orderA) == hash(m_orderB)

    # unique! collapses the duplicate orderings; the mechanisms are not mutated.
    # m_split groups the two binding steps separately, so it is structurally
    # distinct and survives alongside the collapsed orderA/orderB.
    m_split = @enzyme_mechanism begin
        substrates: S, A
        products: P
        steps: begin
            E + S ⇌ ES
            EA + S ⇌ EAS
            E + A ⇌ EA
            EAS <--> EP
            EP ⇌ E + P
        end
    end
    mechs = [m_orderA, m_split, m_orderB]
    snapshot = deepcopy(mechs)
    result = unique!(mechs)
    @test length(result) == 2
    @test all(r -> any(==(r), snapshot), result)
    @test m_orderA == snapshot[1]
    @test m_split  == snapshot[2]
    @test m_orderB == snapshot[3]
    @test (m_orderA.catalytic.enzyme_forms, m_orderA.catalytic.steps) ==
          (snapshot[1].catalytic.enzyme_forms, snapshot[1].catalytic.steps)
end

@testset "_testhelper_assert_mechanism_invariants: coverage of every reactant" begin
    # POSITIVE: the init mechanisms of a reaction with an inhibitor bind no inhibitor
    # yet; the coverage check reads substrates and products only.
    rxn_inh = @enzyme_reaction begin
        substrates: S[C]
        products:   P[C]
        competitive_inhibitors: R
    end
    for m in _testhelper_mechanisms(EnzymeRates.Enumeration.init_mechanisms(rxn_inh))
        @test _testhelper_assert_mechanism_invariants(m) === nothing
    end

    # NEGATIVE: a declared substrate T that no step takes up. The macro rejects such a
    # mechanism; the canonical graph of the walker's raw value, which no constructor
    # check has read, holds it, and the helper must fail on it. This half checks that the
    # oracle can fail.
    @test_throws(ErrorException("Mechanism: substrate T is declared, but no step uses it"),
                 @enzyme_mechanism begin
        substrates: S, T
        products: P
        steps: begin
            E + S ⇌ ES
            ES <--> EP
            EP ⇌ E + P
        end
    end)
    raw_unused = EnzymeRates.DSL.raw_mechanism(:(begin
            substrates: S, T
            products: P
            steps: begin
                E + S ⇌ ES
                ES <--> EP
                EP ⇌ E + P
            end
        end))
    m_unused = _testhelper_plain_of(EnzymeRates.Mechanisms._canonical(
        [raw_unused.substrates; raw_unused.products], raw_unused.forms, [1, 2], [3],
        raw_unused.steps))
    @test_throws(ErrorException("declared substrate/product T appears in no step"),
                 _testhelper_assert_mechanism_invariants(m_unused))
end

@testset "_assert_atom_conserving" begin
    # Lactate dehydrogenase: NADH + Pyr ⇌ Lac + NAD. The substrate and
    # product atom inventories differ per metabolite (NADH ≠ Lac), so an iso
    # that transmutes bound NADH directly into bound Lac with no residual is
    # atom-non-conserving and MUST be rejected.
    ldh_rxn = @enzyme_reaction begin
        substrates: NADH[C21H29N7O14P2], Pyr[C3H3O3]
        products:   Lac[C3H5O3], NAD[C21H27N7O14P2]
    end

    # NEGATIVE: bound NADH → bound Lac, no residual (atoms don't balance). Pyr and NAD
    # bind free E in dead ends, which conserve atoms, so that step is the one offender.
    # Its cycle turns NADH into Lac without Pyr and NAD, which the guard rejects.
    bad_m = @_testhelper_unguarded @enzyme_mechanism begin
        substrates: NADH, Pyr
        products: Lac, NAD
        steps: begin
            E + NADH ⇌ ENADH
            ENADH <--> ELac
            ELac ⇌ E + Lac
            E + Pyr ⇌ EPyr
            E + NAD ⇌ ENAD
        end
    end
    @test_throws ErrorException EnzymeRates.Enumeration._assert_atom_conserving(
        _testhelper_candidate(bad_m), ldh_rxn)
    @test_throws("atom-non-conserving step ENADH → ELac",
                 EnzymeRates.Enumeration._assert_atom_conserving(
                     _testhelper_candidate(bad_m), ldh_rxn))
    # The label names the difference the message prints: Lac has 18 fewer C than NADH.
    label = "atoms(to) + atoms(released) − atoms(from) − atoms(consumed) = "
    @test_throws label EnzymeRates.Enumeration._assert_atom_conserving(
        _testhelper_candidate(bad_m), ldh_rxn)
    @test_throws ":C => -18" EnzymeRates.Enumeration._assert_atom_conserving(
        _testhelper_candidate(bad_m), ldh_rxn)

    # POSITIVE: the bi-bi ping-pong worked example carries a real covalent
    # residual (+A −P) on conformation :E; every step conserves atoms.
    good_m = @enzyme_mechanism begin
        substrates: A, B
        products: P, Q
        steps: begin
            E + A ⇌ EA
            EA <--> EP_res_plusA_minusP
            EP_res_plusA_minusP ⇌ E_res_plusA_minusP + P
            E_res_plusA_minusP + B ⇌ EB_res_plusA_minusP
            EB_res_plusA_minusP ⇌ EQ
            EQ ⇌ E + Q
        end
    end
    @test EnzymeRates.Enumeration._assert_atom_conserving(_testhelper_candidate(good_m),
                                                          bi_bi_pp_rxn) === nothing
end


# ═══════════════════════════════════════════════════════════════════════
# 1. Support functions
# ═══════════════════════════════════════════════════════════════════════

# ─── _catalytic_topologies ──────────────────────────────────────────────
@testset "_catalytic_topologies" begin

ter_ter_topos = EnzymeRates.Enumeration._catalytic_topologies(ter_ter_rxn)
uni_bi_topos = EnzymeRates.Enumeration._catalytic_topologies(uni_bi_rxn)

@testset "topology counts" begin
    for (topos, n) in [
        (EnzymeRates.Enumeration._catalytic_topologies(uni_uni_rxn), 1),
        (uni_bi_topos, 3),
        # 9 sequential topologies. A[C]→P[C] and B[N]→Q[N] each leave no
        # covalent residue, so the only ping-pong is degenerate
        # (empty-residue) and is rejected by the admissible-residual rule.
        (EnzymeRates.Enumeration._catalytic_topologies(bi_bi_rxn), 9),
        (EnzymeRates.Enumeration._catalytic_topologies(bi_bi_pp_rxn), 10),
        (ter_ter_topos, 223),
        # 45 = 39 sequential + 6 genuine-residual ping-pong. P[CN] combines
        # A+B, so ping-pong covalent residuals persist on :E; the degenerate
        # empty-residue variants are rejected by the admissible-residual rule.
        (EnzymeRates.Enumeration._catalytic_topologies(ter_bi_rxn), 45),
    ]
        @test length(topos) == n
        @test all(isempty(_testhelper_connectivity_violations(t)) for t in topos)
        # Every topology has exactly one SS step.
        @test all(count(!s.rapid_equilibrium for s in t.steps) == 1 for t in topos)
    end
    @test !any(_testhelper_is_allosteric(_testhelper_topo_mech(uni_bi_rxn, t))
               for t in uni_bi_topos)
end

@testset "a topology writes each release as the binding it reverses" begin
    ER = EnzymeRates
    form(xs...) =
        ER.Enumeration._FormRecord(:E, Symbol[xs...], Symbol[], Symbol[], Symbol[])
    # Uni-uni: the route binds S at E, isomerizes ES to EP and gives P off from EP.
    # The topology holds the forms in the order the route meets them and writes that
    # release as the binding E + P → EP, which takes P up; it lists its bindings first,
    # by their forms in content order, then its steady-state isomerization.
    # The steps index the reactants S and P, 1 and 2.
    none, S, P = Int[], [1], [2]
    @test ER.Enumeration._catalytic_topologies(uni_uni_rxn) ==
          [(; forms = [form(), form(:S), form(:P)],
            steps = [ER.Mechanisms.Step(1, 3, P, none, true, 0),
                     ER.Mechanisms.Step(1, 2, S, none, true, 0),
                     ER.Mechanisms.Step(2, 3, none, none, false, 0)])]
end

@testset "isomerization patterns come in content order" begin
    # Ter-bi's ping-pongs run their first half-reaction on A and D or on B and D. The
    # patterns sort by their isomerizations' forms in content order, so the one with
    # EA_res_plusB_plusD_minusQ comes before the one with EAD: that residual form binds
    # fewer metabolites, though its name sorts after EAD.
    # Each topology's first isomerization in that order is its steady-state step.
    function _testhelper_iso_texts(t)
        m = _testhelper_topo_mech(ter_bi_rxn, t)
        sort([EnzymeRates.Mechanisms.step_text(m.catalytic, k)
              for (k, e) in enumerate(m.catalytic.steps)
              if EnzymeRates.Mechanisms._is_iso(e)])
    end
    @test unique(_testhelper_iso_texts.(
        EnzymeRates.Enumeration._catalytic_topologies(ter_bi_rxn))) ==
          [["EABD <--> EPQ"],
           ["EA_res_plusB_plusD_minusQ <--> EP", "EBD ⇌ EQ_res_plusB_plusD_minusQ"],
           ["EAD <--> EQ_res_plusA_plusD_minusQ", "EB_res_plusA_plusD_minusQ ⇌ EP"]]
end

@testset "admissible-residual ping-pong" begin
    # The admissible-residual rule rejects degenerate empty-residue
    # ping-pong (which would return the enzyme to apo E mid-cycle,
    # splitting the reaction into disconnected half-cycles). A genuine
    # ping-pong intermediate carries a non-empty covalent residual,
    # always on conformation :E (never a separate conformation). For
    # ter-ter, residues form by combining substrates (e.g. bind A+B,
    # release P[C] leaving an N residue), so ping-pong topologies survive.
    topos = ter_ter_topos
    # Every enzyme form lives on conformation :E.
    @test all(f.label === :E for t in topos for f in t.forms)
    # Surviving ping-pong topologies each carry a genuine (non-empty)
    # covalent intermediate — no empty-residue ping-pong remains.
    pingpong = filter(t -> count(EnzymeRates.Enumeration._is_iso, t.steps) >= 2, topos)
    @test !isempty(pingpong)
    @test all(any(_testhelper_has_residual, t.forms) for t in pingpong)
end

@testset "isomerization constraints" begin
    topos = ter_ter_topos
    sub_names_set = Set([:A, :B, :D])

    # C5: at most max(n_subs, n_prods) = 3 metabolites bound on any form.
    @test all(length(_testhelper_bound(f)) <= 3 for t in topos for f in t.forms)

    # C7: every iso source must contain at
    # least one substrate. Use bound list directly (name-parsing is
    # ambiguous with the concat form naming convention).
    @test all(any(b -> b ∈ sub_names_set,
                  _testhelper_bound(_testhelper_iso_orient(t, s, sub_names_set)[1]))
              for t in topos for s in t.steps if EnzymeRates.Enumeration._is_iso(s))

    # C8: an iso's product-side (destination) form is built with only
    # products bound, never substrates.
    @test all(b ∉ sub_names_set
              for t in topos for s in t.steps if EnzymeRates.Enumeration._is_iso(s)
              for b in _testhelper_bound(_testhelper_iso_orient(t, s, sub_names_set)[2]))
end

_testhelper_bound_names(f) = Set(_testhelper_bound(f))
# The record of an enzyme form binding exactly the metabolites `names`, with or without a
# covalent residual.
_testhelper_form(names, residual) = f ->
    _testhelper_bound_names(f) == Set(names) && _testhelper_has_residual(f) == residual
# Whether the written graph `spec` holds an isomerization with one end satisfying `p` and
# the other `q`, whichever way the step is oriented.
_testhelper_iso_between(spec, p, q) = any(spec.steps) do s
    EnzymeRates.Enumeration._is_iso(s) || return false
    f, t = spec.forms[s.from], spec.forms[s.to]
    (p(f) && q(t)) || (p(t) && q(f))
end

@testset "pyruvate carboxylase mechanism" begin
    topos = EnzymeRates.Enumeration._catalytic_topologies(
        pyruvate_carboxylase_rxn)
    @test all(isempty(_testhelper_connectivity_violations(t)) for t in topos)

    # Known mechanism: ATP+HCO3 → ADP+Pi leaving a CO2 covalent
    # residual on :E, then Pyr+CO2 → OAA. The carboxylation iso converts
    # the {ATP,HCO3}-bound form into a residual-bearing form; the
    # carboxyl-transfer iso converts a residual-bearing Pyr form into the
    # bare EOAA.
    @test any(topos) do spec
        _testhelper_iso_between(spec, _testhelper_form([:ATP, :HCO3], false),
                                _testhelper_has_residual) &&
        _testhelper_iso_between(spec, _testhelper_form([:OAA], false),
            t -> :Pyr in _testhelper_bound_names(t) && _testhelper_has_residual(t))
    end

    # 312 = 169 seq + 143 pp, classified by iso-step count: sequential
    # topologies have one iso step, ping-pong ≥2. Every topology (seq or
    # pp) has exactly one SS step — for ping-pong only one of the iso
    # steps is steady-state, the rest are rapid-equilibrium.
    @test length(topos) == 312
    seq_count = count(t -> count(EnzymeRates.Enumeration._is_iso, t.steps) == 1, topos)
    pp_count = length(topos) - seq_count
    @test seq_count == 169
    @test pp_count == 143
end

@testset "pyruvate dehydrogenase mechanism" begin
    topos = EnzymeRates.Enumeration._catalytic_topologies(
        pyruvate_dehydrogenase_rxn)
    @test all(isempty(_testhelper_connectivity_violations(t)) for t in topos)

    # Known mechanism, with covalent residuals on :E:
    # Pyr→CO2 (leaves an acetyl residual),
    # CoA+acetyl→AcCoA (leaves a hydride residual),
    # NAD+hydride→NADH (residual cancels → bare ENADH).
    @test any(topos) do spec
        _testhelper_iso_between(spec, _testhelper_form([:Pyr], false),
                                _testhelper_form([:CO2], true)) &&
        _testhelper_iso_between(spec, _testhelper_form([:CoA], true),
                                _testhelper_form([:AcCoA], true)) &&
        _testhelper_iso_between(spec, _testhelper_form([:NAD], true),
                                _testhelper_form([:NADH], false))
    end

    # 334 = 169 seq + 165 pp, classified by iso-step count: sequential
    # topologies have one iso step, ping-pong ≥2. Every topology (seq or
    # pp) has exactly one SS step — for ping-pong only one of the iso
    # steps is steady-state, the rest are rapid-equilibrium.
    @test length(topos) == 334
    seq_count = count(t -> count(EnzymeRates.Enumeration._is_iso, t.steps) == 1, topos)
    pp_count = length(topos) - seq_count
    @test seq_count == 169
    @test pp_count == 165
end

@testset "quad-quad: C6 forces ping-pong" begin
    # Quad-quad reaction: 4 subs, 4 prods
    # With C6 (iso ≤ 3×3), 4→4 sequential iso is blocked
    # All topologies must use ping-pong (at least 2 iso steps)
    quad_rxn = @enzyme_reaction begin
        substrates: A[C], B[N], D[X], F[Y]
        products: P[C], Q[N], R[X], S[Y]
    end
    topos = EnzymeRates.Enumeration._catalytic_topologies(quad_rxn)
    @test all(isempty(_testhelper_connectivity_violations(t)) for t in topos)
    @test length(topos) > 0
    # Every topology must have ≥ 2 iso steps (no 4→4)
    @test all(count(EnzymeRates.Enumeration._is_iso, t.steps) >= 2 for t in topos)
end

@testset "no admissible catalytic cycle throws" begin
    # Uni-quad: one isomerization would have to give off all four products (C6), and
    # one that gives off fewer leaves a covalent residue after the last substrate has
    # bound.
    uni_quad_rxn = @enzyme_reaction begin
        substrates: S[CNOX]
        products: P[C], Q[N], R[O], T[X]
    end
    @test_throws(ErrorException(
        "no catalytic cycle for S ⇌ P + Q + R + T: every route needs an isomerization " *
        "converting more than three substrates, or more than three products counting a " *
        "covalent residue as one"), _testhelper_mechanisms(
            EnzymeRates.Enumeration.init_mechanisms(uni_quad_rxn)))
end

end

# ─── _competition_patterns ──────────────────────────────────────────────
@testset "_competition_patterns" begin
# Uni-uni: 1×1, only 1 pattern (single edge)
pats_11 = EnzymeRates.Enumeration._competition_patterns(
    Set([:S]), Set([:P]))
@test length(pats_11) == 1
@test pats_11[1] == Set([(:S, :P)])

# Uni-bi: 1×2, S competes with both P and Q
pats_12 = EnzymeRates.Enumeration._competition_patterns(
    Set([:S]), Set([:P, :Q]))
@test length(pats_12) == 1
@test pats_12[1] == Set([(:S, :P), (:S, :Q)])

# Bi-uni: symmetric
pats_21 = EnzymeRates.Enumeration._competition_patterns(
    Set([:A, :B]), Set([:P]))
@test length(pats_21) == 1
@test pats_21[1] == Set([(:A, :P), (:B, :P)])

# Invalid: {A↔P, B↔P} leaves Q uncovered
@test Set([(:A, :P), (:B, :P)]) ∉ EnzymeRates.Enumeration._competition_patterns(
    Set([:A, :B]), Set([:P, :Q]))

# Every pattern covers all vertices
_testhelper_covers(pat, S, P) = all(any(p -> (s, p) in pat, P) for s in S) &&
                                all(any(s -> (s, p) in pat, S) for p in P)
for (S, P, n) in [
    ([:A, :B], [:P, :Q], 7),              # bi-bi
    ([:A, :B, :C], [:P, :Q, :R], 265),    # ter-ter
    # 2 substrates × 3 products. Bipartite-cover count on K(2,3)
    # by inclusion-exclusion = 25; by symmetry 3 × 2 is also 25.
    ([:A, :B], [:P, :Q, :R], 25),
    ([:A, :B, :D], [:P, :Q], 25),
]
    pats = EnzymeRates.Enumeration._competition_patterns(Set(S), Set(P))
    @test length(pats) == n
    @test all(pat -> _testhelper_covers(pat, S, P), pats)
end
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

@testset "_is_plain_binding and _copy_change read a step's two form records" begin
    # Ordered bi-bi whose last substrate binds straight into the product-bound form, with
    # A's modifier copy at E and EQ, and Q binding at the copy's form too: the mirror of
    # E + Q ⇌ EQ, in its group.
    fused = @enzyme_mechanism begin
        substrates: A, B
        products: P, Q
        steps: begin
            E + A ⇌ EA
            EA + B <--> EPQ
            EQ + P ⇌ EPQ
            (E + Q ⇌ EQ, EAmod + Q ⇌ EAmodQ)
            (E + A ⇌ EAmod, EQ + A ⇌ EAmodQ)
        end
    end
    # A Theorell–Chance step takes up B and gives off P in one step.
    tc = @enzyme_mechanism begin
        substrates: A, B
        products: P, Q
        steps: begin
            E + A ⇌ EA
            EA + B <--> EQ + P
            E + Q ⇌ EQ
        end
    end
    # An isomerization takes up and gives off nothing.
    iso = @enzyme_mechanism begin
        substrates: S
        products: P
        steps: begin
            E + S ⇌ ES
            ES <--> EP
            E + P ⇌ EP
        end
    end
    plain = EnzymeRates.Enumeration._is_plain_binding
    copies = EnzymeRates.Enumeration._copy_change
    # Candidate `m` and its step written `text`.
    at(m, text) = (_testhelper_candidate(m),
                   only(k for k in eachindex(m.catalytic.steps)
                        if EnzymeRates.Mechanisms.step_text(m.catalytic, k) == text))
    # A plain binding of a reactant.
    c, k = at(fused, "E + A ⇌ EA")
    @test (plain(c, k), copies(c, k, :A)) == (true, 0)
    # The fused binding: B binds into EPQ, not into EAB.
    c, k = at(fused, "EA + B <--> EPQ")
    @test (plain(c, k), copies(c, k, :B)) == (false, 0)
    # A binding of the modifier copy: one more copy of A, the productive A unchanged.
    c, k = at(fused, "E + A ⇌ EAmod")
    @test (plain(c, k), copies(c, k, :A)) == (true, 1)
    c, k = at(fused, "EQ + A ⇌ EAmodQ")
    @test (plain(c, k), copies(c, k, :A), copies(c, k, :Q)) == (true, 1, 0)
    # The mirror: Q binds productively beside the copy, which stays.
    c, k = at(fused, "EAmod + Q ⇌ EAmodQ")
    @test (plain(c, k), copies(c, k, :Q), copies(c, k, :A)) == (true, 0, 0)
    c, k = at(tc, "EA + B <--> EQ + P")
    @test (plain(c, k), copies(c, k, :B), copies(c, k, :P)) == (false, 0, 0)
    c, k = at(iso, "ES <--> EP")
    @test (plain(c, k), copies(c, k, :S), copies(c, k, :P)) == (false, 0, 0)
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

# ─── _expand_substrate_product_dead_ends ────────────────────────────────
@testset "_expand_substrate_product_dead_ends" begin
@testset "_expand_substrate_product_dead_ends" begin

@testset "Uni-Uni: no dead-end forms" begin
    # 3 forms: E, E_S[C], E_P[C]. E_S has all subs,
    # E_P has all prods. No mixed dead-end possible.
    # → 0 dead-end forms, 1 variant (bare topology)
    m = @enzyme_mechanism begin
        substrates: S
        products: P
        steps: begin
            E + P ⇌ EP
            E + S ⇌ ES
            ES <--> EP
        end
    end
    topo = _testhelper_flat_topo(m)
    result =
        EnzymeRates.Enumeration._expand_substrate_product_dead_ends(
            [topo],uni_uni_rxn)
    @test all(isempty(_testhelper_connectivity_violations(seed)) for seed in result)
    @test length(result) == 1
end

@testset "Bi-Bi random: 4 dead-end forms" begin
    # 7 forms: E, E_A, E_B, E_A_B, E_P, E_Q, E_P_Q
    # Eligible dead-end forms (mixed sub+prod binding):
    #   E_A: +P→E_A_P(mixed✓), +Q→E_A_Q(mixed✓)
    #   E_B: +P→E_B_P(mixed✓), +Q→E_B_Q(mixed✓)
    #   E_P: +A→E_A_P(same), +B→E_B_P(same)
    #   E_Q: +A→E_A_Q(same), +B→E_B_Q(same)
    # 4 unique mixed-substrate-product forms across competition patterns.
    # Competition patterns for bi-bi (2 subs × 2 prods): 7 patterns
    # (the count from _competition_patterns(2, 2)). Each pattern produces
    # a distinct dead-end-form set:
    #   {A↔P, B↔Q}: forbids E_A_P, E_B_Q → emits {E_A_Q, E_B_P}
    #   {A↔Q, B↔P}: forbids E_A_Q, E_B_P → emits {E_A_P, E_B_Q}
    #   ... (one set per pattern, all distinct)
    #   {A↔P, A↔Q, B↔P, B↔Q}: forbids all → emits {} (bare topology)
    # All 7 sets are distinct → 7 variants after dedup.
    m = @enzyme_mechanism begin
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
    topo = _testhelper_flat_topo(m)
    result =
        EnzymeRates.Enumeration._expand_substrate_product_dead_ends(
            [topo],bi_bi_rxn)
    @test all(isempty(_testhelper_connectivity_violations(seed)) for seed in result)
    # 4 unique dead-end forms, 7 competition patterns,
    # all 7 produce distinct dead-end sets → 7 variants
    @test length(result) == 7

    # Each variant adds exactly the dead-end forms whose (substrate, product)
    # pair its competition pattern leaves unforbidden.
    seed_forms = _testhelper_form_names(bi_bi_rxn, topo)
    @test Set(setdiff(_testhelper_form_names(bi_bi_rxn, r), seed_forms)
              for r in result) == Set([
        Set([:EAQ, :EBP]),     # forbids A↔P, B↔Q
        Set([:EAP, :EBQ]),     # forbids A↔Q, B↔P
        Set([:EBQ]),           # forbids all but B↔Q
        Set([:EBP]),           # forbids all but B↔P
        Set([:EAQ]),           # forbids all but A↔Q
        Set([:EAP]),           # forbids all but A↔P
        Set{Symbol}(),         # forbids all four pairs
    ])
end

@testset "Uni-Bi ordered: no dead-end forms" begin
    # 4 forms: E, E_S, E_P_Q, E_Q
    # E+P→E_P: single-product → rejected (need mixed)
    # E_Q+S→E_S_Q: has all subs → rejected
    # → 0 dead-end forms, 1 variant
    m = @enzyme_mechanism begin
        substrates: S
        products: P, Q
        steps: begin
            E + Q ⇌ EQ
            EQ + P ⇌ EPQ
            E + S ⇌ ES
            ES <--> EPQ
        end
    end
    topo = _testhelper_flat_topo(m)
    result =
        EnzymeRates.Enumeration._expand_substrate_product_dead_ends(
            [topo],uni_bi_rxn)
    @test all(isempty(_testhelper_connectivity_violations(seed)) for seed in result)
    @test length(result) == 1
end

@testset "Bi-Bi Ping-Pong: 6 dead-end forms → 7 variants" begin
    # Forms: E, E_A, Estar, Estar_P, Estar_B, E_Q
    # 6 dead-end forms total, each binding one substrate and one product
    # (E-side: E_A_P, E_A_Q, E_B_Q; Estar-side: Estar_A_P, Estar_B_P,
    # Estar_B_Q). 7 competition patterns; every (substrate, product) pair
    # names at least one dead-end form, so each pattern yields a distinct
    # dead-end-form set after dedup → 7 variants.
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
    topo = _testhelper_flat_topo(m)
    result =
        EnzymeRates.Enumeration._expand_substrate_product_dead_ends(
            [topo],bi_bi_pp_rxn)
    @test all(isempty(_testhelper_connectivity_violations(seed)) for seed in result)
    # 6 dead-end forms (E_A_P, E_A_Q, E_B_Q from
    # E-side + Estar_A_P, Estar_B_P, Estar_B_Q from
    # Estar-side), competition-filtered
    @test length(result) == 7

    # Each variant adds exactly the dead-end forms whose (substrate, product)
    # pair its competition pattern leaves unforbidden.
    seed_forms = _testhelper_form_names(bi_bi_pp_rxn, topo)
    @test Set(setdiff(_testhelper_form_names(bi_bi_pp_rxn, r), seed_forms)
              for r in result) == Set([
        Set([:EAQ, :EstarBP]),                  # forbids A↔P, B↔Q
        Set([:EAP, :EBQ, :EstarAP, :EstarBQ]),  # forbids A↔Q, B↔P
        Set([:EBQ, :EstarBQ]),                  # forbids all but B↔Q
        Set([:EstarBP]),                        # forbids all but B↔P
        Set([:EAQ]),                            # forbids all but A↔Q
        Set([:EAP, :EstarAP]),                  # forbids all but A↔P
        Set{Symbol}(),                          # forbids all four pairs
    ])

    # Assert that some result variants contain Estar-prefixed dead-end
    # forms (proving dead-end forms inherit the base form's Estar
    # conformation).
    new_estar_forms = Set{Symbol}()
    for r in result
        new_forms = setdiff(_testhelper_form_names(bi_bi_pp_rxn, r), seed_forms)
        for f in new_forms
            startswith(string(f), "Estar") && string(f) != "Estar" &&
                push!(new_estar_forms, f)
        end
    end
    @test !isempty(new_estar_forms)
end

@testset "Ter-ter per-topology (OOM on full init)" begin
    # Test that competition filtering works
    # on representative ter-ter topologies.
    topos = EnzymeRates.Enumeration._catalytic_topologies(
        ter_ter_rxn)
    @test length(topos) == 223
    # Test first (random, most forms) and last topology
    for topo in [topos[1], topos[end]]
        result =
            EnzymeRates.Enumeration._expand_substrate_product_dead_ends(
                [topo], ter_ter_rxn)
        @test all(isempty(_testhelper_connectivity_violations(seed)) for seed in result)
        # Competition patterns reduce 2^27 to
        # ≤265 variants per topology
        @test length(result) > 0
        @test length(result) <= 265
    end
end

@testset "Ter-ter random: 27 dead-end forms, 12 on the diagonal pattern" begin
    # The random ter-ter topology (the one with most forms) has 27 dead-end forms.
    # The diagonal competition pattern {A↔P, B↔Q, D↔R} keeps 12 of them:
    #   1S+1P: 6 (the diagonal pairs A-P, B-Q, D-R are forbidden)
    #   2S+1P: 3 (EABR, EADQ, EBDP)
    #   1S+2P: 3 (EAQR, EBPR, EDPQ)
    topos = EnzymeRates.Enumeration._catalytic_topologies(ter_ter_rxn)
    _, idx = findmax(length(_testhelper_form_names(ter_ter_rxn, t)) for t in topos)
    result = EnzymeRates.Enumeration._expand_substrate_product_dead_ends([topos[idx]],
                                                                         ter_ter_rxn)
    added = [setdiff(_testhelper_form_names(ter_ter_rxn, seed),
                     _testhelper_form_names(ter_ter_rxn, topos[idx]))
             for seed in result]
    @test length(union(added...)) == 27
    @test Set([:EAQ, :EAR, :EBP, :EBR, :EDP, :EDQ, :EABR, :EADQ, :EBDP, :EAQR, :EBPR,
               :EDPQ]) in added
end

end

@testset "shared_catalytic_site prunes catalytic-site co-occupancy" begin
rxn = @enzyme_reaction begin
    substrates: A[C], B[C]
    products:   P[C], Q[C]
    shared_catalytic_site: (A, P)
end
# Whether the form whose record is `f` holds both A and P as reactants, not as modifier
# copies.
_testhelper_binds_both(f) = :A in f.productive && :P in f.productive
@test !any(_testhelper_binds_both(f)
           for m in _testhelper_mechanisms(EnzymeRates.Enumeration.init_mechanisms(rxn))
           for f in _testhelper_candidate(m).forms)
end

@testset "shared_catalytic_site strictly reduces mechanism count" begin
base = @enzyme_reaction begin
    substrates: A[C], B[C]
    products:   P[C], Q[C]
end
constrained = @enzyme_reaction begin
    substrates: A[C], B[C]
    products:   P[C], Q[C]
    shared_catalytic_site: (A, P)
end
n_base = length(
    unique!(collect(_testhelper_mechanisms(EnzymeRates.Enumeration.init_mechanisms(base)))))
n_con  = length(unique!(
    collect(_testhelper_mechanisms(EnzymeRates.Enumeration.init_mechanisms(constrained)))))
@test n_con < n_base
end
end

# ═══════════════════════════════════════════════════════════════════════
# Testsets covering non-enumeration features (allosteric mechanism
# accessor identity)
# ═══════════════════════════════════════════════════════════════════════

@testset "Allosteric mechanism TR equivalence" begin
# Three distinct tags on three groups: S binding (RE, :EqualAI),
# P binding (RE, :NonequalAI), iso (SS, :OnlyA). The :OnlyA sits on
# the chemical step (dropped from the Haldane graph), so the mechanism
# is Haldane-valid with no one-sided :OnlyA binding.
m_compiled = @enzyme_mechanism begin
    substrates: S
    products: P
    oligomeric_state: 2
    steps: begin
        E + S ⇌ ES             :: EqualAI
        E + P ⇌ EP             :: NonequalAI
        ES <--> EP           :: OnlyA
    end
end
# Group order is canonical; allosteric tags stay bound to their steps.
g = _testhelper_graph(m_compiled)
# The tag of the one group whose first step satisfies `pred`.
_testhelper_state_of(pred) = only(g.tags[s.group] for s in _testhelper_first_steps(g)
                                  if pred(s))
_testhelper_binds_in(names) = s -> _testhelper_is_binding(s) && only(s.consumed) in names
@test _testhelper_state_of(_testhelper_binds_in(g.substrates)) == :EqualAI
@test _testhelper_state_of(_testhelper_binds_in(g.products)) == :NonequalAI
@test _testhelper_state_of(!_testhelper_is_binding) == :OnlyA
end

# ═══════════════════════════════════════════════════════════════════════
# 2. Initialization (init_mechanisms)
# ═══════════════════════════════════════════════════════════════════════

# ─── init_mechanisms ───────────────────────────────────────────────────
@testset "init_mechanisms" begin

# The init mechanisms of uni_uni_rxn, uni_bi_rxn, bi_bi_rxn and bi_bi_pp_rxn.
uu, ub, bb, pp = (_testhelper_mechanisms(EnzymeRates.Enumeration.init_mechanisms(rxn))
                  for rxn in [uni_uni_rxn, uni_bi_rxn, bi_bi_rxn, bi_bi_pp_rxn])

@testset "init mechanisms are connected, with one SS step per seed" begin
    # Uni-uni topology: 1 catalytic topology × 1 dead-end variant
    # (none possible — see test_expand_substrate_product_dead_ends
    # uni-uni case). Hence init produces exactly 1 mechanism.
    @test length(uu) == 1
    # A seed holds its chemistry as an isomerization, and exactly one step,
    # an isomerization, is SS by construction. Subsequent RE→SS expansions add
    # more SS steps. The merged and Theorell–Chance variants that follow the
    # seeds hold no isomerization and two or three SS groups.
    for specs in (uu, ub, bb, pp)
        @test all(isempty(_testhelper_connectivity_violations(m)) for m in specs)
        for s in specs
            if _testhelper_holds_iso(s)
                @test count(e -> !e.rapid_equilibrium, s.catalytic.steps) == 1
            else
                @test count(fs -> !first(fs), _testhelper_group_flags(s)) in (2, 3)
            end
        end
    end
end

@testset "_seed_groups labels one metabolite's bindings at one flag as one group" begin
    ER = EnzymeRates
    # A step from form `f` to form `t` of E, EA, EB, EAB, EPQ, EP and EQ,
    # numbered 1 to 7, taking up `xs`.
    # The metabolites A and B are 1 and 2.
    index = Dict(:A => 1, :B => 2)
    step(f, t, xs, r) = ER.Mechanisms.Step(f, t, Int[index[x] for x in xs], Int[], r, 0)
    # The A bindings at rapid equilibrium share label 1; B's rapid-equilibrium and
    # steady-state bindings take labels 2 and 3; each isomerization has a label of its
    # own though both take up and give off nothing. Labels run in order of first use.
    steps = [step(1, 2, [:A], true), step(1, 3, [:B], true), step(3, 4, [:A], true),
             step(2, 4, [:B], false), step(4, 5, [], false), step(6, 7, [], false)]
    @test ER.Enumeration._seed_groups(steps) == [1, 2, 1, 3, 4, 5]
end

@testset "Same-metabolite RE bindings share kinetic_group" begin
    # _seed_groups collapses all RE binding steps for
    # the same metabolite into one kinetic group (one shared K).
    # For bi-bi, metabolites like :B appear in multiple binding steps
    # (e.g. E+B⇌E_B and E_A+B⇌E_A_B) — these must share one kinetic_group.
    @test !isempty(bb)
    n_assertions_fired = 0
    for spec in bb
        # Map each metabolite to the kinetic-group index (inner-vector
        # position) of every RE binding step that binds it. Same-group
        # sharing is structural: same inner vector == same kinetic group.
        by_metabolite = Dict{Symbol, Vector{Int}}()
        for step in _testhelper_graph(spec).steps
            step.rapid_equilibrium && _testhelper_is_binding(step) || continue
            push!(get!(by_metabolite, only(step.consumed), Int[]), step.group)
        end
        for (_met, gis) in by_metabolite
            length(gis) >= 2 || continue
            @test length(Set(gis)) == 1
            n_assertions_fired += 1
        end
    end
    @test n_assertions_fired >= 1   # at least one multi-binding case existed
end

@testset "init_mechanisms: the seeds, then their merged and Theorell–Chance variants" begin
    ER = EnzymeRates
    function _testhelper_kind(m)
        tc = any(_testhelper_is_tc, m.catalytic.steps)
        c = _testhelper_candidate(m)
        fused = any(e -> _testhelper_is_fused(c, e), m.catalytic.steps)
        (tc ? (fused ? :half_tc : :tc) : :merged,
         count(fs -> !first(fs), _testhelper_group_flags(m)))
    end
    _testhelper_tally(ms) = Dict(k => count(m -> _testhelper_kind(m) == k, ms)
                                 for k in unique(_testhelper_kind.(ms)))
    # Aggregate pins over the whole seed sets. bi_bi_pp_rxn: 62 seeds (55 sequential, 7
    # ping-pong), each holding its isomerization, then 202 variants, none holding one.
    # Two-group merged: 52 from the 20 ordered/ordered seeds (3 each, less the 8 lumping
    # twins), 56 from the 28 ordered/random and random/ordered seeds (2 each), 14 from the
    # 7 ping-pong seeds (PP-AQ and PP-BP each). Three-group merged: 28 from the
    # ordered/random and random/ordered seeds (1 each) and 28 from the 7 random/random
    # seeds (4 each). Theorell–Chance: 20, one per ordered/ordered seed. Half
    # Theorell–Chance: 4, both halves of the undecorated ping-pong seed and one half of
    # each of two decorated ones. The 7 ping-pong seeds are degenerate and no variant is.
    @test length(pp) == 264 && allunique(pp)
    @test all(_testhelper_holds_iso, pp[1:62]) && !any(_testhelper_holds_iso, pp[63:end])
    @test _testhelper_tally(pp[63:end]) == Dict((:merged, 2) => 122, (:merged, 3) => 56,
                                                (:tc, 3) => 20, (:half_tc, 3) => 4)
    @test count(ER.Enumeration._degenerate, _testhelper_candidate.(pp[1:62])) == 7 &&
          !any(ER.Enumeration._degenerate, _testhelper_candidate.(pp[63:end]))
    @test all(m -> ER.Enumeration._assert_emission_rules(
        _testhelper_candidate(m)) === nothing, pp)
    # bi_bi_rxn's atoms admit no ping-pong: the 55 sequential seeds and their 184 variants.
    @test length(bb) == 239 && allunique(bb)
    @test all(_testhelper_holds_iso, bb[1:55]) && !any(_testhelper_holds_iso, bb[56:end])
    @test _testhelper_tally(bb[56:end]) ==
          Dict((:merged, 2) => 108, (:merged, 3) => 56, (:tc, 3) => 20)
end

@testset "init_mechanisms on uni-bi: one merged variant per ordered seed" begin
    ER = EnzymeRates.Enumeration
    # Aggregate pin over the uni-bi seed set: the two ordered seeds (P or Q released first)
    # and the random-release seed. An ordered seed's merged base holds S (E + S → EPQ),
    # the inner product's group (into EPQ) and the outer product's group. Singletons
    # fail V; {S, inner} is the lumping twin; {inner, outer} has no reverse maximal rate;
    # {S, outer} is the one variant (n + m − 2 = 1). Its Theorell–Chance base gives none:
    # each single group leaves the TC step RE or closes the cycle, and {TC, outer} turned
    # RE for products closes it as well. The random-release base holds S, Pˣ and Qˣ, each
    # product group with a step into EPQ; {S, Pˣ} and {S, Qˣ} pass and are no twins
    # (EPQ has three steps), {Pˣ, Qˣ} has no reverse maximal rate, and ELIM refuses the
    # three-step complex. Each variant fits 2 + 1 + 2 − 1 = 4 (2 + 2 + 1 − 1 = 4 for the
    # random seed's), the seeds' count.
    ordered_p = @enzyme_mechanism begin
        substrates: S; products: P, Q
        steps: begin
            E + S <--> EPQ
            EQ + P ⇌ EPQ
            E + Q <--> EQ
        end
    end
    ordered_q = @enzyme_mechanism begin
        substrates: S; products: P, Q
        steps: begin
            E + S <--> EPQ
            EP + Q ⇌ EPQ
            E + P <--> EP
        end
    end
    random_p = @enzyme_mechanism begin
        substrates: S; products: P, Q
        steps: begin
            E + S <--> EPQ
            (EQ + P <--> EPQ, E + P <--> EP)
            (EP + Q ⇌ EPQ, E + Q ⇌ EQ)
        end
    end
    random_q = @enzyme_mechanism begin
        substrates: S; products: P, Q
        steps: begin
            E + S <--> EPQ
            (EQ + P ⇌ EPQ, E + P ⇌ EP)
            (EP + Q <--> EPQ, E + Q <--> EQ)
        end
    end
    @test length(ub) == 7
    @test all(_testhelper_holds_iso, ub[1:3])
    @test Set(ub[4:7]) == Set([ordered_p, ordered_q, random_p, random_q])
    for seed in ub[1:3]
        @test length(_testhelper_move(ER._seed_variants, seed, uni_bi_rxn)) ==
              (length(seed.catalytic.steps) == 4 ? 1 : 2)
    end
    for v in ub[4:7]
        fitted = _testhelper_fitted(v)
        @test (fitted, _testhelper_identifiable_rank(v)) == (4, 4)
    end
end

"""The value of one call of `f()` and the bytes the call allocates, after a warm-up call
that compiles it."""
function _testhelper_value_and_bytes(f)
    f()
    run = @timed f()
    run.value, run.bytes
end

@testset "the enumeration holds and allocates a bounded number of bytes" begin
    # Aggregate check over the whole bi-bi starting set. On the macOS CI runner the ter-ter
    # call spent 13.9-116.1 s in garbage collection in three runs after mechanisms stopped
    # storing facts derived from their steps, against 171.0 s in one run before, so no
    # gain there is established (ROADMAP). Time and peak memory depend on the machine, so
    # these bounds gate what does not: the same bytes in every process on a Julia version,
    # measured on Julia 1.11.9 and 1.13.4.
    rxn = @enzyme_reaction begin
        substrates: A[C], B[N]
        products: P[C], Q[N]
    end
    parent = @enzyme_mechanism begin
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
    start, start_bytes = _testhelper_value_and_bytes(
        () -> EnzymeRates.Enumeration.candidate_mechanisms(rxn))
    held = Base.summarysize(start) / length(start)
    candidate = _testhelper_candidate(parent)
    children, expand_bytes = _testhelper_value_and_bytes(
        () -> EnzymeRates.Enumeration.candidate_mechanisms(rxn, candidate))
    @info "bi-bi enumeration: $(round(Int, held)) bytes held per candidate; " *
          "$start_bytes bytes allocated by the starting set, $expand_bytes by one expansion"
    # The bounds hold for a known amount of work: fewer candidates or children would
    # allocate less and pass.
    @test length(start) == 239
    @test length(children) == 21
    # 1,593.6 bytes, the same in every process on Julia 1.11.9 and 1.13.4, so the margin
    # of 3.5% only allows for an untested Julia version. A field derived from the steps
    # and stored on `Mechanism`, `Step` or the candidate's record of a form raises it: one
    # `Int` vector on `Mechanism` takes it to 1,700.
    @test held < 1_650
    # Each allocation bound sits 25% above the larger of its two values; more work per
    # candidate raises it. 42,545,112 bytes on 1.11, 38,788,376 on 1.13:
    @test start_bytes < 53_000_000
    # 13,750,248 bytes on 1.13, 12,683,288 on 1.11:
    @test expand_bytes < 17_200_000
end

@testset "init_mechanisms on ter-ter within its time and allocation budgets" begin
    # Aggregate pin over the whole ter-ter seed set: 35,665 seeds and their 215,442 merged
    # and Theorell–Chance variants that the guard admits. Garbage collection varies with
    # the heap, from about 23% of the call on CI to about 70% under a local heap-size hint,
    # so the budget excludes it. Without it the call takes 38 s in a full local suite and
    # 59-109 s on a CI runner; a macOS runner that took 430 s in all would be near 330 s.
    # The bytes the call allocates do not depend on garbage collection and were the same
    # in every process on a Julia version; peak resident memory depends on how the
    # collector sizes the heap, so it is logged and not bounded.
    run = @timed _testhelper_mechanisms(
        EnzymeRates.Enumeration.init_mechanisms(ter_ter_rxn))
    @info "init_mechanisms on ter-ter: $(round(run.time; digits = 1)) s, " *
          "$(round(run.gctime; digits = 1)) s of it in garbage collection; " *
          "$(round(run.bytes / 2^30; digits = 2)) GiB allocated; " *
          "the process's peak resident memory $(round(Sys.maxrss() / 2^30; digits = 2)) GiB"
    @test length(run.value) == 251107
    @test run.time - run.gctime < 400
    # 89.3 GiB on Julia 1.11.9, 81.2 GiB on 1.13.4; more work per candidate raises it.
    @test run.bytes < 112 * 2^30
end

@testset "Drops unbound regulators from init Mechanism" begin
    # init_mechanisms produces Mechanisms without dead-end regulators
    # bound. In a built Mechanism, the regulator must NOT
    # appear among the modifiers its steps bind (`_testhelper_regulators`)
    # — only the catalytic mechanism is built. After an expansion
    # adds the dead-end regulator, it should appear.
    init_mechs = _testhelper_mechanisms(
        EnzymeRates.Enumeration.init_mechanisms(uni_uni_with_reg))
    @test all(isempty(_testhelper_connectivity_violations(m)) for m in init_mechs)
    @test !isempty(init_mechs)
    for m in init_mechs
        @test :I ∉ _testhelper_regulators(m)
    end

    expanded = _testhelper_children(uni_uni_with_reg, init_mechs)
    found_with_reg = false
    for mm in expanded
        if :I in _testhelper_regulators(mm)
            found_with_reg = true
            break
        end
    end
    @test found_with_reg
end

@testset "Substrate-as-product overlap (racemase shape)" begin
    # A substrate racemase has differently-named substrate/product
    # (e.g., L-Ala → D-Ala) but the same atomic composition. Init
    # mechanisms must compile correctly.
    rxn = @enzyme_reaction begin
        substrates: LAla[CHN]
        products: DAla[CHN]
    end
    specs = _testhelper_mechanisms(EnzymeRates.Enumeration.init_mechanisms(rxn))
    @test all(isempty(_testhelper_connectivity_violations(m)) for m in specs)
    @test !isempty(specs)
    @test !any(_testhelper_is_allosteric, specs)
end

end

# ═══════════════════════════════════════════════════════════════════════
# 3. Expansion moves
# ═══════════════════════════════════════════════════════════════════════

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

@testset "Allosteric mechanism — a chain flank still flips at multiplicity 2" begin
    # A K-type variant of the ordered seed over two catalytic subunits: the iso and the
    # A binding :OnlyA, the other bindings :EqualAI. Its flank flips stay children: two
    # conformations over two subunits can see a flank's flip. The A-state graph is the
    # plain seed's, so each single flip cuts the one segment in two, keeps a bottom
    # form and carries flux. The parent needs a hyperbolic scheme
    # (`_hyperbolic_catalysis`). With two segments an arborescence toward a root holds
    # one edge, which leaves the other segment, so the scheme fails only when a
    # metabolite scores 2 on one edge or one segment, or scores in one segment and on
    # an edge leaving the other. Each flip passes:
    #   A: {E, EQ, EPQ} | {EA, EAB}. A scores on E → EA only; B in the
    #      second segment and on EAB → EPQ, which leaves it; P and Q in the
    #      first and on EPQ → EAB, which leaves it.
    #   Q: the mirror of A.
    #   B: {E, EA, EQ, EPQ} | {EAB}. A, P and Q score in the first segment
    #      and on edges leaving it (EA → EAB, EPQ → EAB); B scores 1 on
    #      EA → EAB and in no segment.
    #   P: the mirror of B.
    # Four children, tags kept. The parent is identifiable at 6; every child adds one
    # constant and gains one, the B and P flank flips included, where the plain
    # seed's flank flips gain none.
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
    @test isempty(EnzymeRates.Enumeration._chain_flank_groups(am))
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
    b_flip = @enzyme_mechanism begin
        substrates: A, B ; products: P, Q ; oligomeric_state: 2
        steps: begin
            E + A ⇌ EA                :: OnlyA
            EA + B <--> EAB       :: EqualAI
            EAB <--> EPQ        :: OnlyA
            EQ + P ⇌ EPQ          :: EqualAI
            E + Q ⇌ EQ                :: EqualAI
        end
    end
    p_flip = @enzyme_mechanism begin
        substrates: A, B ; products: P, Q ; oligomeric_state: 2
        steps: begin
            E + A ⇌ EA                :: OnlyA
            EA + B ⇌ EAB          :: EqualAI
            EAB <--> EPQ        :: OnlyA
            EQ + P <--> EPQ       :: EqualAI
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
    @test length(kids) == 4
    @test Set(kids) == Set([a_flip, b_flip, p_flip, q_flip])
    @test all(EnzymeRates.Enumeration._hyperbolic_catalysis, _testhelper_candidate.(kids))
    @test _testhelper_fitted(am) == _testhelper_identifiable_rank(am) == 6
    for k in kids
        @test _testhelper_fitted(k) == _testhelper_identifiable_rank(k) == 7
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

@testset "Allosteric mechanism — uni-uni RE→SS flips one group and preserves the tags" begin
    # Each seed is a uni-uni allosteric mechanism whose two binding groups are all-RE and
    # whose isomerization is SS: _expand_re_to_ss fires per all-RE group → 2 variants,
    # one with the P binding steady state and one with the S binding steady state, each
    # keeping every tag.
    for (am, deltas, expected) in [
        # :EqualAI groups. Each converted group stays :EqualAI (RE K → SS (kf, kr),
        # shared R/T). Cheap-tag RE→SS adds 1 fitted param per variant (no separate
        # T-state pair).
        ((@enzyme_mechanism begin
            substrates: S
            products: P
            oligomeric_state: 2
            steps: begin
                E + P ⇌ EP      :: EqualAI
                E + S ⇌ ES      :: EqualAI
                ES <--> EP    :: EqualAI
            end
        end), [1, 1], [
            (@enzyme_mechanism begin
                substrates: S
                products: P
                oligomeric_state: 2
                steps: begin
                    E + P <--> EP   :: EqualAI
                    E + S ⇌ ES      :: EqualAI
                    ES <--> EP    :: EqualAI
                end
            end),
            (@enzyme_mechanism begin
                substrates: S
                products: P
                oligomeric_state: 2
                steps: begin
                    E + P ⇌ EP      :: EqualAI
                    E + S <--> ES   :: EqualAI
                    ES <--> EP    :: EqualAI
                end
            end)]),
        # Both binding groups :OnlyA (a balanced substrate/product pair, so the
        # mechanism is Haldane-valid) and :EqualAI catalysis. :OnlyA groups live in the
        # R-state only; the T-state contributes no kf_T/kr_T after RE→SS. The move MUST
        # NOT change R/T-state semantics, so the converted group stays :OnlyA.
        # Connected-component pruning populates the inactive state from free E over ALL
        # surviving steps, so because catalysis is :EqualAI the base mechanism's Q_I
        # ALREADY carries ES and EP through the inactive catalytic step — even
        # though both bindings are :OnlyA and cannot bind directly in the T-state.
        # Those forms are therefore present in the inactive state of the base AND of
        # every variant, so flipping a binding group RE→SS adds only that group's own SS
        # rate constant (+1) — it does not introduce a new inactive-state form. Both
        # variants → Δ=1.
        ((@enzyme_mechanism begin
            substrates: S
            products: P
            oligomeric_state: 2
            steps: begin
                E + P ⇌ EP      :: OnlyA
                E + S ⇌ ES      :: OnlyA
                ES <--> EP    :: EqualAI
            end
        end), [1, 1], [
            (@enzyme_mechanism begin
                substrates: S
                products: P
                oligomeric_state: 2
                steps: begin
                    E + P <--> EP   :: OnlyA
                    E + S ⇌ ES      :: OnlyA
                    ES <--> EP    :: EqualAI
                end
            end),
            (@enzyme_mechanism begin
                substrates: S
                products: P
                oligomeric_state: 2
                steps: begin
                    E + P ⇌ EP      :: OnlyA
                    E + S <--> ES   :: OnlyA
                    ES <--> EP    :: EqualAI
                end
            end)]),
        # One :NonequalAI group (S-binding), the others :EqualAI. When RE→SS converts a
        # :NonequalAI group, BOTH the R-state K and the T-state K_T must split into
        # (kf, kr) and (kf_T, kr_T): P-binding :EqualAI → +1; S-binding :NonequalAI → +2.
        # SS :NonequalAI adds kf_T, kr_T on top of the R-state pair.
        ((@enzyme_mechanism begin
            substrates: S
            products: P
            oligomeric_state: 2
            steps: begin
                E + P ⇌ EP      :: EqualAI
                E + S ⇌ ES      :: NonequalAI
                ES <--> EP    :: EqualAI
            end
        end), [1, 2], [
            (@enzyme_mechanism begin
                substrates: S
                products: P
                oligomeric_state: 2
                steps: begin
                    E + P <--> EP   :: EqualAI
                    E + S ⇌ ES      :: NonequalAI
                    ES <--> EP    :: EqualAI
                end
            end),
            (@enzyme_mechanism begin
                substrates: S
                products: P
                oligomeric_state: 2
                steps: begin
                    E + P ⇌ EP      :: EqualAI
                    E + S <--> ES   :: NonequalAI
                    ES <--> EP    :: EqualAI
                end
            end)]),
        # Every group :NonequalAI — a distinguishable mechanism (an all-:EqualAI seed, the
        # first entry above, is a model-space no-op that is never enumerated; the move is
        # still tested on it). :NonequalAI RE→SS adds 2 fitted params per variant (the A-
        # and I-state binding K each split into (kf, kr)).
        ((@enzyme_mechanism begin
            substrates: S
            products: P
            oligomeric_state: 2
            steps: begin
                E + S ⇌ ES    :: NonequalAI
                E + P ⇌ EP    :: NonequalAI
                ES <--> EP  :: NonequalAI
            end
        end), [2, 2], [
            (@enzyme_mechanism begin
                substrates: S
                products: P
                oligomeric_state: 2
                steps: begin
                    E + P <--> EP   :: NonequalAI
                    E + S ⇌ ES      :: NonequalAI
                    ES <--> EP    :: NonequalAI
                end
            end),
            (@enzyme_mechanism begin
                substrates: S
                products: P
                oligomeric_state: 2
                steps: begin
                    E + P ⇌ EP      :: NonequalAI
                    E + S <--> ES   :: NonequalAI
                    ES <--> EP    :: NonequalAI
                end
            end)]),
    ]
        _testhelper_assert_mechanism_invariants(am)
        result = _testhelper_move(EnzymeRates.Enumeration._expand_re_to_ss, am)
        @test length(result) == 2
        @test Set(result) == Set(expected)
        @test _testhelper_param_deltas(am, result) == deltas
        for r in result
            _testhelper_assert_mechanism_invariants(r)
            # Exactly one group newly all-SS; every tag, the multiplicity, the
            # regulatory sites and the reaction preserved.
            n_newly_ss = count(((old, new),) -> all(old) && !any(new),
                               zip(_testhelper_group_flags(am), _testhelper_group_flags(r)))
            @test n_newly_ss == 1
            @test r.allosteric.allosteric_states == am.allosteric.allosteric_states
            @test r.allosteric.oligomeric_state == am.allosteric.oligomeric_state
            @test r.allosteric.regulatory_sites == am.allosteric.regulatory_sites
            @test _testhelper_keeps_reaction(r, am)
        end
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
    # All four groups A, B, P and Q flip: an allosteric mechanism has no chain
    # flanks (`_chain_flank_groups`). Every child keeps one `:OnlyA` group, the
    # chemistry, and the multiplicity and (empty) regulatory sites.
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
                EA + B <--> EAB                                   :: EqualAI
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
                EQ + P <--> EPQ                                   :: EqualAI
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
    @test length(kids) == 4
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
    # random-order pair flip would carry B². The three children flip the P
    # binding, the B binding, or the A binding with its mirror.
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
                E + A ⇌ EA          :: EqualAI
                EA + B ⇌ EAB    :: EqualAI
                EAB <--> EP      :: EqualAI
                E + P <--> EP       :: EqualAI
                E + I ⇌ EI          :: EqualAI
                EA + I ⇌ EAI    :: EqualAI
                EI + A ⇌ EAI    :: EqualAI
            end
        end),
        @enzyme_mechanism(begin
            substrates: A, B
            products: P
            modifiers: I
            steps: begin
                E + A ⇌ EA          :: EqualAI
                EA + B <--> EAB  :: EqualAI
                EAB <--> EP      :: EqualAI
                E + P ⇌ EP          :: EqualAI
                E + I ⇌ EI          :: EqualAI
                EA + I ⇌ EAI    :: EqualAI
                EI + A ⇌ EAI    :: EqualAI
            end
        end),
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
    @test length(kids) == 3
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

@testset "Allosteric mechanism — a stopped :OnlyA binding strands EA; no child" begin
    # Merged random-order bi-bi with the abortive complex EAQ. The steady-state A group
    # binds A at E and, fused, at EB into EPQ, all :OnlyA. In the inactive
    # conformation the :EqualAI steps E + Q <--> EQ, EQ + A <--> EAQ and
    # EA + Q <--> EAQ still join E to EA, the two sides of the :OnlyA binding
    # E + A <--> EA, so that conformation stops it. At zero P and Q, EA then leaves
    # only by taking up Q: it holds enzyme apart from free E, so the steady state there is
    # not unique (`_NONUNIQUE`), as both oracles find. Every move keeps that binding :OnlyA
    # and those steps untagged: a split's parts inherit their group's tag, and the
    # change-tag move relaxes binding groups only, while the A group holds chemistry. So no
    # child has a limit either. Of change-tag's three relaxations, Q's group and Q's
    # abortive binding of A keep EA stranded; P's group joins the two sides of the
    # rapid-equilibrium binding EA + B ⇌ EPQ, which then cannot be switched off
    # (`_FAST_STEP_PINNED`). The flip move emits nothing, and the exit drops the three
    # relaxations. Four split sets raise the count, but the split search asks the guard on
    # each child, so they count as failed, and no superset passes: the move emits no
    # split.
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
    relaxed = _testhelper_move(EnzymeRates.Enumeration._expand_change_allo_state, am)
    @test length(relaxed) == 3
    @test Set(relaxed) == Set([p_relaxed, q_relaxed, qa_relaxed])
    split = _testhelper_move(EnzymeRates.Enumeration._expand_split_kinetic_group, am)
    @test length(split) == 0
    @test Set(split) == Set()
    moved = [_testhelper_move(EnzymeRates.Enumeration._expand_re_to_ss, am); relaxed]
    @test Dict(c => D._guard(c).fault for c in moved) ==
          Dict(p_relaxed => D._FAST_STEP_PINNED, q_relaxed => D._NONUNIQUE,
               qa_relaxed => D._NONUNIQUE)
    kids = _testhelper_children(bi_bi_rxn, [am])
    @test length(kids) == 0
    @test Set(kids) == Set()
end

@testset "Allosteric mechanism — beside an inhibitor, stopped :OnlyA A strands EA" begin
    # Merged random-order bi-bi with the competitive inhibitor I and the abortive
    # complexes EAP and EBQ; the A, B, I and Q groups :OnlyA. In the inactive
    # conformation the :EqualAI steps E + P <--> EP, EP + A <--> EAP and
    # EA + P <--> EAP join E to EA, the two sides of the :OnlyA binding
    # E + A <--> EA, so that conformation stops it. At zero P and Q, EA then leaves
    # only by taking up P: it holds enzyme apart from free E, so the steady state there is
    # not unique (`_NONUNIQUE`), as both oracles find. No child has a limit either. Of
    # change-tag's four relaxations, I's group, P's group and P's abortive binding of A
    # keep EA stranded; Q's group joins the two sides of the rapid-equilibrium binding
    # EA + B ⇌ EPQ, which then cannot be switched off (`_FAST_STEP_PINNED`). The flip
    # move emits nothing, and the exit drops the four relaxations. Six split sets raise the
    # count, but the split search asks the guard on each child, so they count as failed,
    # and no superset passes: the move emits no split.
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
    relaxed = _testhelper_move(EnzymeRates.Enumeration._expand_change_allo_state, am)
    @test length(relaxed) == 4
    @test Set(relaxed) == Set([i_relaxed, p_relaxed, q_relaxed, pa_relaxed])
    split = _testhelper_move(EnzymeRates.Enumeration._expand_split_kinetic_group, am)
    @test length(split) == 0
    @test Set(split) == Set()
    moved = [_testhelper_move(EnzymeRates.Enumeration._expand_re_to_ss, am); relaxed]
    @test Dict(c => D._guard(c).fault for c in moved) ==
          Dict(i_relaxed => D._NONUNIQUE, p_relaxed => D._NONUNIQUE,
               q_relaxed => D._FAST_STEP_PINNED, pa_relaxed => D._NONUNIQUE)
    kids = _testhelper_children(rxn, [am])
    @test length(kids) == 0
    @test Set(kids) == Set()
end

@testset "Allosteric mechanism — a split keeps its parent's interior verdict" begin
    # An aggregate check over the bi-bi seed set (the exception to writing every fixture
    # inline): each split child of a seed tagged as below, one parent in
    # ten. Both parts of a split group inherit its tag, and admission's interior reads
    # each edge's tag, so every split the move builds (zero-flux parts at rapid
    # equilibrium) has its parent's interior verdict: the split search's guard rejects
    # every split of a parent the interior refuses.
    # Each split child's every-chemistry tag set with each nonempty set of its
    # binding groups :OnlyA: the tagging move's K-types on that set, among them the
    # parents admission refuses, which the reject direction needs.
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
    # visible; the placement stands and adds an identifiable constant. With S
    # binding `:NonequalAI` the copy duplicates ES in both states, the S group
    # alone forms ES in each, so the gauge exists in each, and it is skipped.
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
    kids = _testhelper_move(EnzymeRates.Enumeration._expand_add_dead_end_regulator, onlya,
                            rxn)
    @test length(kids) == 1
    @test Set(kids) == Set([kept])
    @test _testhelper_identifiable_rank(kept) == _testhelper_identifiable_rank(onlya) + 1
    @test isempty(
        EnzymeRates.Enumeration._redundant_copy_groups(_testhelper_candidate(kept)))
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
    # constant is visible; the placement stands and adds an identifiable constant.
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
    kids = _testhelper_move(EnzymeRates.Enumeration._expand_add_dead_end_regulator,
                            dead_inactive, rxn)
    @test length(kids) == 1
    @test Set(kids) == Set([kept])
    @test isempty(
        EnzymeRates.Enumeration._redundant_copy_groups(_testhelper_candidate(kept)))
    @test _testhelper_identifiable_rank(kept) ==
        _testhelper_identifiable_rank(dead_inactive) + 1
end

@testset "Allosteric mechanism — no relaxation makes a kept copy redundant" begin
    # The copy of S at E is new only in the inactive state, where S does not bind.
    # Relaxing the S binding to `:NonequalAI` brings ES into the inactive state:
    # the copy then duplicates ES in both states, formed in each by the S group
    # alone, so the gauge exists in each, its constant shows only beside K_S, and
    # the child is not emitted (one fitted constant above its rank).
    # The move never relaxes the chemistry while S binds the active state only, and
    # admission would refuse it: the chemistry and the P binding would join E to ES
    # in the inactive state, leaving S's affinity no path to zero. The P binding and
    # the copy relax. The relaxed copy, with one
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

@testset "Allosteric mechanism — dead-end binding tagged :EqualAI" begin
    # Uni-uni allosteric (catalytic_n=2) with mixed regs:
    # :I dead-end, :R allosteric. The move must (a) exclude :R from
    # eligible regs (allosteric ligand), (b) append the new dead-end
    # group's allosteric_states tag as :EqualAI.
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

    result = _testhelper_move(EnzymeRates.Enumeration._expand_add_dead_end_regulator, am,
                              rxn)

    # 1. count: 1 (same as plain uni-uni + I; :R is excluded as
    # allosteric).
    @test length(result) == 1

    # 2. Δ params: +1 (:EqualAI new group → one shared K), measured by `parameters`,
    # the exact count of the constants a fit estimates.
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

    # 4. property: exactly one new catalytic group; its tag is
    # :EqualAI. Pre-existing tags are unchanged.
    for r in result
        @test length(_testhelper_group_members(r.catalytic)) ==
              length(_testhelper_group_members(am.catalytic)) + 1
        @test length(r.allosteric.allosteric_states) ==
              length(am.allosteric.allosteric_states) + 1
        # The new group's tag is :EqualAI (appended at the end by the dead-end
        # move).
        @test r.allosteric.allosteric_states[end] == :EqualAI
        # Pre-existing tags preserved positionally.
        @test r.allosteric.allosteric_states[1:length(am.allosteric.allosteric_states)] ==
            am.allosteric.allosteric_states
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

# ═══════════════════════════════════════════════════════════════════════
# 4. Allosteric expansion moves
# ═══════════════════════════════════════════════════════════════════════

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

@testset "Mechanism — Bi-bi sequential: binding-group :OnlyA variants only" begin
    # 5 kinetic groups: 4 bindings (A, B substrate-side; P, Q
    # product-side) + 1 chemical step (EAB <--> EPQ). The
    # catalytic group is never a bare primary :OnlyA (that V-type needs
    # a regulator, none declared). Instead each nonempty set of binding
    # groups is tagged :OnlyA on top of the base whose chemistry is :OnlyA,
    # and the move emits each; the exits keep those the guard accepts
    # (`check_valid_rate_equation`).
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

    result = _testhelper_move(EnzymeRates.Enumeration._expand_to_allosteric, m,
                              bi_bi_allo_rxn)

    # 1. count: every non-empty subset of the 4 binding groups is set
    # :OnlyA, each with the chemical step also :OnlyA (a catalytically-dead
    # inactive conformation). 2^4 - 1 = 15 subsets, all admitted.
    n_groups = length(_testhelper_group_members(m.catalytic))
    n_cat = count(g -> _testhelper_is_catalytic_group(m, g), 1:n_groups)
    @test n_cat == 1
    @test length(result) == 15

    # 2. Δ params: every variant is +1 — the conformational constant L.
    # Each :OnlyA binding is K-type (the inactive conformation cannot bind
    # that ligand, adding no inactive binding constant), and the dead
    # inactive conformation runs no chemistry, so it carries no free
    # catalytic constant. All 15 variants are Δ=1.
    @test _testhelper_param_deltas(m, result) == fill(1, 15)

    # 3. invariants (the deltas above compile every result), and no result is the
    # all-:EqualAI baseline.
    for r in result
        _testhelper_assert_mechanism_invariants(r)
        @test !all(==(:EqualAI), r.allosteric.allosteric_states)
    end
end

@testset "Mechanism — Bi-bi ping-pong: binding-group :OnlyA variants only" begin
    # SEED: bi-bi ping-pong, 6 kinetic groups: 4 bindings (A and Q at free E, B and P at
    # the covalent form E_res_plusA_minusP) and the two half-reactions EA <-->
    # EP_res_plusA_minusP and EB_res_plusA_minusP ⇌ EQ. Tagging both halves
    # :OnlyA cuts the covalent forms off in the inactive conformation, and the guard
    # admits it, so that set is the one base. No regulator is declared, so no V-type: the
    # children are K-types, each nonempty set of the four bindings :OnlyA on the base, all
    # 15 emitted. The exits keep five, the sets without B's binding that tag P's binding
    # only beside A's: a set with B's binding pumps a basin at zero products, and one with
    # P's binding but not A's pumps one at zero substrates.
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
    a_half1 = @_testhelper_unguarded @enzyme_mechanism begin
        substrates: A, B ; products: P, Q ; oligomeric_state: 2
        steps: begin
            E + A ⇌ EA                                          :: OnlyA
            E_res_plusA_minusP + B ⇌ EB_res_plusA_minusP    :: EqualAI
            E + Q ⇌ EQ                                          :: EqualAI
            E_res_plusA_minusP + P ⇌ EP_res_plusA_minusP    :: EqualAI
            EA <--> EP_res_plusA_minusP                      :: OnlyA
            EB_res_plusA_minusP ⇌ EQ                         :: EqualAI
        end
    end
    p_half1 = @_testhelper_unguarded @enzyme_mechanism begin
        substrates: A, B ; products: P, Q ; oligomeric_state: 2
        steps: begin
            E + A ⇌ EA                                          :: EqualAI
            E_res_plusA_minusP + B ⇌ EB_res_plusA_minusP    :: EqualAI
            E + Q ⇌ EQ                                          :: EqualAI
            E_res_plusA_minusP + P ⇌ EP_res_plusA_minusP    :: OnlyA
            EA <--> EP_res_plusA_minusP                      :: OnlyA
            EB_res_plusA_minusP ⇌ EQ                         :: EqualAI
        end
    end
    ap_half1 = @_testhelper_unguarded @enzyme_mechanism begin
        substrates: A, B ; products: P, Q ; oligomeric_state: 2
        steps: begin
            E + A ⇌ EA                                          :: OnlyA
            E_res_plusA_minusP + B ⇌ EB_res_plusA_minusP    :: EqualAI
            E + Q ⇌ EQ                                          :: EqualAI
            E_res_plusA_minusP + P ⇌ EP_res_plusA_minusP    :: OnlyA
            EA <--> EP_res_plusA_minusP                      :: OnlyA
            EB_res_plusA_minusP ⇌ EQ                         :: EqualAI
        end
    end
    a_every = @_testhelper_unguarded @enzyme_mechanism begin
        substrates: A, B ; products: P, Q ; oligomeric_state: 2
        steps: begin
            E + A ⇌ EA                                       :: OnlyA
            E + Q ⇌ EQ                                       :: EqualAI
            E_res_plusA_minusP + B ⇌ EB_res_plusA_minusP :: EqualAI
            E_res_plusA_minusP + P ⇌ EP_res_plusA_minusP :: EqualAI
            EA <--> EP_res_plusA_minusP                   :: OnlyA
            EB_res_plusA_minusP ⇌ EQ                      :: OnlyA
        end
    end
    q_every = @_testhelper_unguarded @enzyme_mechanism begin
        substrates: A, B ; products: P, Q ; oligomeric_state: 2
        steps: begin
            E + A ⇌ EA                                       :: EqualAI
            E + Q ⇌ EQ                                       :: OnlyA
            E_res_plusA_minusP + B ⇌ EB_res_plusA_minusP :: EqualAI
            E_res_plusA_minusP + P ⇌ EP_res_plusA_minusP :: EqualAI
            EA <--> EP_res_plusA_minusP                   :: OnlyA
            EB_res_plusA_minusP ⇌ EQ                      :: OnlyA
        end
    end
    aq_every = @_testhelper_unguarded @enzyme_mechanism begin
        substrates: A, B ; products: P, Q ; oligomeric_state: 2
        steps: begin
            E + A ⇌ EA                                       :: OnlyA
            E + Q ⇌ EQ                                       :: OnlyA
            E_res_plusA_minusP + B ⇌ EB_res_plusA_minusP :: EqualAI
            E_res_plusA_minusP + P ⇌ EP_res_plusA_minusP :: EqualAI
            EA <--> EP_res_plusA_minusP                   :: OnlyA
            EB_res_plusA_minusP ⇌ EQ                      :: OnlyA
        end
    end
    b_every = @_testhelper_unguarded @enzyme_mechanism begin
        substrates: A, B ; products: P, Q ; oligomeric_state: 2
        steps: begin
            E + A ⇌ EA                                       :: EqualAI
            E + Q ⇌ EQ                                       :: EqualAI
            E_res_plusA_minusP + B ⇌ EB_res_plusA_minusP :: OnlyA
            E_res_plusA_minusP + P ⇌ EP_res_plusA_minusP :: EqualAI
            EA <--> EP_res_plusA_minusP                   :: OnlyA
            EB_res_plusA_minusP ⇌ EQ                      :: OnlyA
        end
    end
    ab_every = @_testhelper_unguarded @enzyme_mechanism begin
        substrates: A, B ; products: P, Q ; oligomeric_state: 2
        steps: begin
            E + A ⇌ EA                                       :: OnlyA
            E + Q ⇌ EQ                                       :: EqualAI
            E_res_plusA_minusP + B ⇌ EB_res_plusA_minusP :: OnlyA
            E_res_plusA_minusP + P ⇌ EP_res_plusA_minusP :: EqualAI
            EA <--> EP_res_plusA_minusP                   :: OnlyA
            EB_res_plusA_minusP ⇌ EQ                      :: OnlyA
        end
    end
    bq_every = @_testhelper_unguarded @enzyme_mechanism begin
        substrates: A, B ; products: P, Q ; oligomeric_state: 2
        steps: begin
            E + A ⇌ EA                                       :: EqualAI
            E + Q ⇌ EQ                                       :: OnlyA
            E_res_plusA_minusP + B ⇌ EB_res_plusA_minusP :: OnlyA
            E_res_plusA_minusP + P ⇌ EP_res_plusA_minusP :: EqualAI
            EA <--> EP_res_plusA_minusP                   :: OnlyA
            EB_res_plusA_minusP ⇌ EQ                      :: OnlyA
        end
    end
    abq_every = @_testhelper_unguarded @enzyme_mechanism begin
        substrates: A, B ; products: P, Q ; oligomeric_state: 2
        steps: begin
            E + A ⇌ EA                                       :: OnlyA
            E + Q ⇌ EQ                                       :: OnlyA
            E_res_plusA_minusP + B ⇌ EB_res_plusA_minusP :: OnlyA
            E_res_plusA_minusP + P ⇌ EP_res_plusA_minusP :: EqualAI
            EA <--> EP_res_plusA_minusP                   :: OnlyA
            EB_res_plusA_minusP ⇌ EQ                      :: OnlyA
        end
    end
    p_every = @_testhelper_unguarded @enzyme_mechanism begin
        substrates: A, B ; products: P, Q ; oligomeric_state: 2
        steps: begin
            E + A ⇌ EA                                       :: EqualAI
            E + Q ⇌ EQ                                       :: EqualAI
            E_res_plusA_minusP + B ⇌ EB_res_plusA_minusP :: EqualAI
            E_res_plusA_minusP + P ⇌ EP_res_plusA_minusP :: OnlyA
            EA <--> EP_res_plusA_minusP                   :: OnlyA
            EB_res_plusA_minusP ⇌ EQ                      :: OnlyA
        end
    end
    ap_every = @_testhelper_unguarded @enzyme_mechanism begin
        substrates: A, B ; products: P, Q ; oligomeric_state: 2
        steps: begin
            E + A ⇌ EA                                       :: OnlyA
            E + Q ⇌ EQ                                       :: EqualAI
            E_res_plusA_minusP + B ⇌ EB_res_plusA_minusP :: EqualAI
            E_res_plusA_minusP + P ⇌ EP_res_plusA_minusP :: OnlyA
            EA <--> EP_res_plusA_minusP                   :: OnlyA
            EB_res_plusA_minusP ⇌ EQ                      :: OnlyA
        end
    end
    pq_every = @_testhelper_unguarded @enzyme_mechanism begin
        substrates: A, B ; products: P, Q ; oligomeric_state: 2
        steps: begin
            E + A ⇌ EA                                       :: EqualAI
            E + Q ⇌ EQ                                       :: OnlyA
            E_res_plusA_minusP + B ⇌ EB_res_plusA_minusP :: EqualAI
            E_res_plusA_minusP + P ⇌ EP_res_plusA_minusP :: OnlyA
            EA <--> EP_res_plusA_minusP                   :: OnlyA
            EB_res_plusA_minusP ⇌ EQ                      :: OnlyA
        end
    end
    apq_every = @_testhelper_unguarded @enzyme_mechanism begin
        substrates: A, B ; products: P, Q ; oligomeric_state: 2
        steps: begin
            E + A ⇌ EA                                       :: OnlyA
            E + Q ⇌ EQ                                       :: OnlyA
            E_res_plusA_minusP + B ⇌ EB_res_plusA_minusP :: EqualAI
            E_res_plusA_minusP + P ⇌ EP_res_plusA_minusP :: OnlyA
            EA <--> EP_res_plusA_minusP                   :: OnlyA
            EB_res_plusA_minusP ⇌ EQ                      :: OnlyA
        end
    end
    bp_every = @_testhelper_unguarded @enzyme_mechanism begin
        substrates: A, B ; products: P, Q ; oligomeric_state: 2
        steps: begin
            E + A ⇌ EA                                       :: EqualAI
            E + Q ⇌ EQ                                       :: EqualAI
            E_res_plusA_minusP + B ⇌ EB_res_plusA_minusP :: OnlyA
            E_res_plusA_minusP + P ⇌ EP_res_plusA_minusP :: OnlyA
            EA <--> EP_res_plusA_minusP                   :: OnlyA
            EB_res_plusA_minusP ⇌ EQ                      :: OnlyA
        end
    end
    abp_every = @_testhelper_unguarded @enzyme_mechanism begin
        substrates: A, B ; products: P, Q ; oligomeric_state: 2
        steps: begin
            E + A ⇌ EA                                       :: OnlyA
            E + Q ⇌ EQ                                       :: EqualAI
            E_res_plusA_minusP + B ⇌ EB_res_plusA_minusP :: OnlyA
            E_res_plusA_minusP + P ⇌ EP_res_plusA_minusP :: OnlyA
            EA <--> EP_res_plusA_minusP                   :: OnlyA
            EB_res_plusA_minusP ⇌ EQ                      :: OnlyA
        end
    end
    bpq_every = @_testhelper_unguarded @enzyme_mechanism begin
        substrates: A, B ; products: P, Q ; oligomeric_state: 2
        steps: begin
            E + A ⇌ EA                                       :: EqualAI
            E + Q ⇌ EQ                                       :: OnlyA
            E_res_plusA_minusP + B ⇌ EB_res_plusA_minusP :: OnlyA
            E_res_plusA_minusP + P ⇌ EP_res_plusA_minusP :: OnlyA
            EA <--> EP_res_plusA_minusP                   :: OnlyA
            EB_res_plusA_minusP ⇌ EQ                      :: OnlyA
        end
    end
    abpq_every = @_testhelper_unguarded @enzyme_mechanism begin
        substrates: A, B ; products: P, Q ; oligomeric_state: 2
        steps: begin
            E + A ⇌ EA                                       :: OnlyA
            E + Q ⇌ EQ                                       :: OnlyA
            E_res_plusA_minusP + B ⇌ EB_res_plusA_minusP :: OnlyA
            E_res_plusA_minusP + P ⇌ EP_res_plusA_minusP :: OnlyA
            EA <--> EP_res_plusA_minusP                   :: OnlyA
            EB_res_plusA_minusP ⇌ EQ                      :: OnlyA
        end
    end
    result = _testhelper_move(EnzymeRates.Enumeration._expand_to_allosteric, m,
                              bi_bi_pp_allo_rxn)
    @test length(result) == 15
    @test Set(result) == Set([a_every, q_every, aq_every, b_every, ab_every, bq_every,
                              abq_every, p_every, ap_every, pq_every, apq_every, bp_every,
                              abp_every, bpq_every, abpq_every])
    # The first half alone :OnlyA is no base: class 3 rejects its K-types on the faces
    # with every substrate zero and with B zero, a split whose limit the speeds do not
    # move but which differs from the law (limit test `differ`).
    for c in (a_half1, p_half1, ap_half1)
        v = EnzymeRates.Derivation._guard(c)
        names =
            EnzymeRates.Derivation._face_names(
                m.catalytic, EnzymeRates.Derivation._faces(m.catalytic)[v.face])
        @test (v.fault, v.limit, names) == (EnzymeRates.Derivation._SPLIT,
                                            EnzymeRates.Derivation._LIMIT_DIFFER, [:A, :B])
    end
    kept = [a_every, q_every, aq_every, ap_every, apq_every]
    @test all(_testhelper_is_valid, kept)
    @test !any(_testhelper_is_valid, setdiff(result, kept))
    @test Set(c for c in _testhelper_children(bi_bi_pp_allo_rxn, [m])
              if _testhelper_is_allosteric(c)) == Set(kept)
    # Δ params: each kept child adds L only; every tag is :EqualAI or :OnlyA.
    @test _testhelper_param_deltas(m, kept) == fill(1, 5)
    foreach(_testhelper_assert_mechanism_invariants, kept)
end

@testset "Mechanism — merged abortive ping-pong: hyperbolic laws above one subunit" begin
    # SEED: the bi-bi ping-pong above with A's abortive binding at EP_res_plusA_minusP
    # merged into A's group and P's at EA into P's group. At two subunits the move keeps
    # a child only when both conformations' laws are hyperbolic (`_hyperbolic_catalysis`).
    # Both half-reactions :OnlyA is the one base: the inactive conformation runs no
    # chemistry, so its abortive complexes put no power into its law, and the move emits
    # all 15 K-types. The first half alone :OnlyA is no base: admission accepts it, class
    # 3 rejects it. Its K-types (A's group, P's group or both :OnlyA) would have
    # hyperbolic laws, but class 3 rejects all three. The exits keep five K-types, as for
    # the seed without abortive bindings: a set with B's binding pumps a basin at zero
    # products, and one with P's group but not A's pumps one at zero substrates.
    m = @enzyme_mechanism begin
        substrates: A, B
        products: P, Q
        steps: begin
            (E + A ⇌ EA, EP_res_plusA_minusP + A ⇌ EAP_res_plusA_minusP)
            E + Q ⇌ EQ
            E_res_plusA_minusP + B ⇌ EB_res_plusA_minusP
            (E_res_plusA_minusP + P ⇌ EP_res_plusA_minusP, EA + P ⇌ EAP)
            EA <--> EP_res_plusA_minusP
            EB_res_plusA_minusP ⇌ EQ
        end
    end
    rxn = @enzyme_reaction begin
        substrates: A[CX], B[N]
        products: P[C], Q[NX]
        oligomeric_state: 2
    end
    a_half1 = @_testhelper_unguarded @enzyme_mechanism begin
        substrates: A, B ; products: P, Q ; oligomeric_state: 2
        steps: begin
            (E + A ⇌ EA,
             EP_res_plusA_minusP + A ⇌ EAP_res_plusA_minusP)   :: OnlyA
            E + Q ⇌ EQ                                               :: EqualAI
            E_res_plusA_minusP + B ⇌ EB_res_plusA_minusP         :: EqualAI
            (E_res_plusA_minusP + P ⇌ EP_res_plusA_minusP,
             EA + P ⇌ EAP)                                       :: EqualAI
            EA <--> EP_res_plusA_minusP                           :: OnlyA
            EB_res_plusA_minusP ⇌ EQ                              :: EqualAI
        end
    end
    p_half1 = @_testhelper_unguarded @enzyme_mechanism begin
        substrates: A, B ; products: P, Q ; oligomeric_state: 2
        steps: begin
            (E + A ⇌ EA,
             EP_res_plusA_minusP + A ⇌ EAP_res_plusA_minusP)   :: EqualAI
            E + Q ⇌ EQ                                               :: EqualAI
            E_res_plusA_minusP + B ⇌ EB_res_plusA_minusP         :: EqualAI
            (E_res_plusA_minusP + P ⇌ EP_res_plusA_minusP,
             EA + P ⇌ EAP)                                       :: OnlyA
            EA <--> EP_res_plusA_minusP                           :: OnlyA
            EB_res_plusA_minusP ⇌ EQ                              :: EqualAI
        end
    end
    ap_half1 = @_testhelper_unguarded @enzyme_mechanism begin
        substrates: A, B ; products: P, Q ; oligomeric_state: 2
        steps: begin
            (E + A ⇌ EA,
             EP_res_plusA_minusP + A ⇌ EAP_res_plusA_minusP)   :: OnlyA
            E + Q ⇌ EQ                                               :: EqualAI
            E_res_plusA_minusP + B ⇌ EB_res_plusA_minusP         :: EqualAI
            (E_res_plusA_minusP + P ⇌ EP_res_plusA_minusP,
             EA + P ⇌ EAP)                                       :: OnlyA
            EA <--> EP_res_plusA_minusP                           :: OnlyA
            EB_res_plusA_minusP ⇌ EQ                              :: EqualAI
        end
    end
    a_every = @_testhelper_unguarded @enzyme_mechanism begin
        substrates: A, B ; products: P, Q ; oligomeric_state: 2
        steps: begin
            (E + A ⇌ EA,
             EP_res_plusA_minusP + A ⇌ EAP_res_plusA_minusP) :: OnlyA
            E + Q ⇌ EQ                                             :: EqualAI
            E_res_plusA_minusP + B ⇌ EB_res_plusA_minusP       :: EqualAI
            (E_res_plusA_minusP + P ⇌ EP_res_plusA_minusP,
             EA + P ⇌ EAP)                                     :: EqualAI
            EA <--> EP_res_plusA_minusP                         :: OnlyA
            EB_res_plusA_minusP ⇌ EQ                            :: OnlyA
        end
    end
    q_every = @_testhelper_unguarded @enzyme_mechanism begin
        substrates: A, B ; products: P, Q ; oligomeric_state: 2
        steps: begin
            (E + A ⇌ EA,
             EP_res_plusA_minusP + A ⇌ EAP_res_plusA_minusP) :: EqualAI
            E + Q ⇌ EQ                                             :: OnlyA
            E_res_plusA_minusP + B ⇌ EB_res_plusA_minusP       :: EqualAI
            (E_res_plusA_minusP + P ⇌ EP_res_plusA_minusP,
             EA + P ⇌ EAP)                                     :: EqualAI
            EA <--> EP_res_plusA_minusP                         :: OnlyA
            EB_res_plusA_minusP ⇌ EQ                            :: OnlyA
        end
    end
    aq_every = @_testhelper_unguarded @enzyme_mechanism begin
        substrates: A, B ; products: P, Q ; oligomeric_state: 2
        steps: begin
            (E + A ⇌ EA,
             EP_res_plusA_minusP + A ⇌ EAP_res_plusA_minusP) :: OnlyA
            E + Q ⇌ EQ                                             :: OnlyA
            E_res_plusA_minusP + B ⇌ EB_res_plusA_minusP       :: EqualAI
            (E_res_plusA_minusP + P ⇌ EP_res_plusA_minusP,
             EA + P ⇌ EAP)                                     :: EqualAI
            EA <--> EP_res_plusA_minusP                         :: OnlyA
            EB_res_plusA_minusP ⇌ EQ                            :: OnlyA
        end
    end
    b_every = @_testhelper_unguarded @enzyme_mechanism begin
        substrates: A, B ; products: P, Q ; oligomeric_state: 2
        steps: begin
            (E + A ⇌ EA,
             EP_res_plusA_minusP + A ⇌ EAP_res_plusA_minusP) :: EqualAI
            E + Q ⇌ EQ                                             :: EqualAI
            E_res_plusA_minusP + B ⇌ EB_res_plusA_minusP       :: OnlyA
            (E_res_plusA_minusP + P ⇌ EP_res_plusA_minusP,
             EA + P ⇌ EAP)                                     :: EqualAI
            EA <--> EP_res_plusA_minusP                         :: OnlyA
            EB_res_plusA_minusP ⇌ EQ                            :: OnlyA
        end
    end
    ab_every = @_testhelper_unguarded @enzyme_mechanism begin
        substrates: A, B ; products: P, Q ; oligomeric_state: 2
        steps: begin
            (E + A ⇌ EA,
             EP_res_plusA_minusP + A ⇌ EAP_res_plusA_minusP) :: OnlyA
            E + Q ⇌ EQ                                             :: EqualAI
            E_res_plusA_minusP + B ⇌ EB_res_plusA_minusP       :: OnlyA
            (E_res_plusA_minusP + P ⇌ EP_res_plusA_minusP,
             EA + P ⇌ EAP)                                     :: EqualAI
            EA <--> EP_res_plusA_minusP                         :: OnlyA
            EB_res_plusA_minusP ⇌ EQ                            :: OnlyA
        end
    end
    bq_every = @_testhelper_unguarded @enzyme_mechanism begin
        substrates: A, B ; products: P, Q ; oligomeric_state: 2
        steps: begin
            (E + A ⇌ EA,
             EP_res_plusA_minusP + A ⇌ EAP_res_plusA_minusP) :: EqualAI
            E + Q ⇌ EQ                                             :: OnlyA
            E_res_plusA_minusP + B ⇌ EB_res_plusA_minusP       :: OnlyA
            (E_res_plusA_minusP + P ⇌ EP_res_plusA_minusP,
             EA + P ⇌ EAP)                                     :: EqualAI
            EA <--> EP_res_plusA_minusP                         :: OnlyA
            EB_res_plusA_minusP ⇌ EQ                            :: OnlyA
        end
    end
    abq_every = @_testhelper_unguarded @enzyme_mechanism begin
        substrates: A, B ; products: P, Q ; oligomeric_state: 2
        steps: begin
            (E + A ⇌ EA,
             EP_res_plusA_minusP + A ⇌ EAP_res_plusA_minusP) :: OnlyA
            E + Q ⇌ EQ                                             :: OnlyA
            E_res_plusA_minusP + B ⇌ EB_res_plusA_minusP       :: OnlyA
            (E_res_plusA_minusP + P ⇌ EP_res_plusA_minusP,
             EA + P ⇌ EAP)                                     :: EqualAI
            EA <--> EP_res_plusA_minusP                         :: OnlyA
            EB_res_plusA_minusP ⇌ EQ                            :: OnlyA
        end
    end
    p_every = @_testhelper_unguarded @enzyme_mechanism begin
        substrates: A, B ; products: P, Q ; oligomeric_state: 2
        steps: begin
            (E + A ⇌ EA,
             EP_res_plusA_minusP + A ⇌ EAP_res_plusA_minusP) :: EqualAI
            E + Q ⇌ EQ                                             :: EqualAI
            E_res_plusA_minusP + B ⇌ EB_res_plusA_minusP       :: EqualAI
            (E_res_plusA_minusP + P ⇌ EP_res_plusA_minusP,
             EA + P ⇌ EAP)                                     :: OnlyA
            EA <--> EP_res_plusA_minusP                         :: OnlyA
            EB_res_plusA_minusP ⇌ EQ                            :: OnlyA
        end
    end
    ap_every = @_testhelper_unguarded @enzyme_mechanism begin
        substrates: A, B ; products: P, Q ; oligomeric_state: 2
        steps: begin
            (E + A ⇌ EA,
             EP_res_plusA_minusP + A ⇌ EAP_res_plusA_minusP) :: OnlyA
            E + Q ⇌ EQ                                             :: EqualAI
            E_res_plusA_minusP + B ⇌ EB_res_plusA_minusP       :: EqualAI
            (E_res_plusA_minusP + P ⇌ EP_res_plusA_minusP,
             EA + P ⇌ EAP)                                     :: OnlyA
            EA <--> EP_res_plusA_minusP                         :: OnlyA
            EB_res_plusA_minusP ⇌ EQ                            :: OnlyA
        end
    end
    pq_every = @_testhelper_unguarded @enzyme_mechanism begin
        substrates: A, B ; products: P, Q ; oligomeric_state: 2
        steps: begin
            (E + A ⇌ EA,
             EP_res_plusA_minusP + A ⇌ EAP_res_plusA_minusP) :: EqualAI
            E + Q ⇌ EQ                                             :: OnlyA
            E_res_plusA_minusP + B ⇌ EB_res_plusA_minusP       :: EqualAI
            (E_res_plusA_minusP + P ⇌ EP_res_plusA_minusP,
             EA + P ⇌ EAP)                                     :: OnlyA
            EA <--> EP_res_plusA_minusP                         :: OnlyA
            EB_res_plusA_minusP ⇌ EQ                            :: OnlyA
        end
    end
    apq_every = @_testhelper_unguarded @enzyme_mechanism begin
        substrates: A, B ; products: P, Q ; oligomeric_state: 2
        steps: begin
            (E + A ⇌ EA,
             EP_res_plusA_minusP + A ⇌ EAP_res_plusA_minusP) :: OnlyA
            E + Q ⇌ EQ                                             :: OnlyA
            E_res_plusA_minusP + B ⇌ EB_res_plusA_minusP       :: EqualAI
            (E_res_plusA_minusP + P ⇌ EP_res_plusA_minusP,
             EA + P ⇌ EAP)                                     :: OnlyA
            EA <--> EP_res_plusA_minusP                         :: OnlyA
            EB_res_plusA_minusP ⇌ EQ                            :: OnlyA
        end
    end
    bp_every = @_testhelper_unguarded @enzyme_mechanism begin
        substrates: A, B ; products: P, Q ; oligomeric_state: 2
        steps: begin
            (E + A ⇌ EA,
             EP_res_plusA_minusP + A ⇌ EAP_res_plusA_minusP) :: EqualAI
            E + Q ⇌ EQ                                             :: EqualAI
            E_res_plusA_minusP + B ⇌ EB_res_plusA_minusP       :: OnlyA
            (E_res_plusA_minusP + P ⇌ EP_res_plusA_minusP,
             EA + P ⇌ EAP)                                     :: OnlyA
            EA <--> EP_res_plusA_minusP                         :: OnlyA
            EB_res_plusA_minusP ⇌ EQ                            :: OnlyA
        end
    end
    abp_every = @_testhelper_unguarded @enzyme_mechanism begin
        substrates: A, B ; products: P, Q ; oligomeric_state: 2
        steps: begin
            (E + A ⇌ EA,
             EP_res_plusA_minusP + A ⇌ EAP_res_plusA_minusP) :: OnlyA
            E + Q ⇌ EQ                                             :: EqualAI
            E_res_plusA_minusP + B ⇌ EB_res_plusA_minusP       :: OnlyA
            (E_res_plusA_minusP + P ⇌ EP_res_plusA_minusP,
             EA + P ⇌ EAP)                                     :: OnlyA
            EA <--> EP_res_plusA_minusP                         :: OnlyA
            EB_res_plusA_minusP ⇌ EQ                            :: OnlyA
        end
    end
    bpq_every = @_testhelper_unguarded @enzyme_mechanism begin
        substrates: A, B ; products: P, Q ; oligomeric_state: 2
        steps: begin
            (E + A ⇌ EA,
             EP_res_plusA_minusP + A ⇌ EAP_res_plusA_minusP) :: EqualAI
            E + Q ⇌ EQ                                             :: OnlyA
            E_res_plusA_minusP + B ⇌ EB_res_plusA_minusP       :: OnlyA
            (E_res_plusA_minusP + P ⇌ EP_res_plusA_minusP,
             EA + P ⇌ EAP)                                     :: OnlyA
            EA <--> EP_res_plusA_minusP                         :: OnlyA
            EB_res_plusA_minusP ⇌ EQ                            :: OnlyA
        end
    end
    abpq_every = @_testhelper_unguarded @enzyme_mechanism begin
        substrates: A, B ; products: P, Q ; oligomeric_state: 2
        steps: begin
            (E + A ⇌ EA, EP_res_plusA_minusP + A ⇌ EAP_res_plusA_minusP) :: OnlyA
            E + Q ⇌ EQ                                                           :: OnlyA
            E_res_plusA_minusP + B ⇌ EB_res_plusA_minusP                     :: OnlyA
            (E_res_plusA_minusP + P ⇌ EP_res_plusA_minusP,
             EA + P ⇌ EAP)                                                   :: OnlyA
            EA <--> EP_res_plusA_minusP                                       :: OnlyA
            EB_res_plusA_minusP ⇌ EQ                                          :: OnlyA
        end
    end
    kids = _testhelper_move(EnzymeRates.Enumeration._expand_to_allosteric, m, rxn)
    @test length(kids) == 15
    @test Set(kids) == Set([a_every, q_every, aq_every, b_every, ab_every, bq_every,
                            abq_every, p_every, ap_every, pq_every, apq_every, bp_every,
                            abp_every, bpq_every, abpq_every])
    @test all(EnzymeRates.Enumeration._hyperbolic_where_required,
              _testhelper_candidate.(kids))
    kept = [a_every, q_every, aq_every, ap_every, apq_every]
    @test all(_testhelper_is_valid, kept)
    @test !any(_testhelper_is_valid, setdiff(kids, kept))
    @test all(EnzymeRates.Enumeration._hyperbolic_where_required,
              _testhelper_candidate.([a_half1, p_half1, ap_half1]))
    @test !any(_testhelper_is_valid, [a_half1, p_half1, ap_half1])
    @test Set(c for c in _testhelper_children(rxn, [m])
              if _testhelper_is_allosteric(c)) == Set(kept)
end

@testset "Mechanism — the tagging move drops a K-type whose law is not hyperbolic" begin
    # A uni-bi parent: S binds fused into EPQ at steady state, and the products leave
    # in either order, each at steady state from EPQ and at rapid equilibrium from
    # EP or EQ. Its scheme is hyperbolic, so two subunits stay open, and at two
    # subunits the move keeps a child only when both conformations' schemes are hyperbolic
    # (`_hyperbolic_where_required`). The one base tags the fused S binding :OnlyA, and the
    # move builds a K-type for each nonempty set of the four binding groups :OnlyA on it.
    # It drops two: with P's binding at E alone :OnlyA too, the inactive conformation
    # reaches EP only through steady-state steps and its scheme fails
    # `_hyperbolic_catalysis`; likewise Q's binding at E alone. The guard rejects both as
    # well (the inactive conformation still joins the two sides of the excluded binding),
    # so the exits would drop them anyway; this testset pins the move's own drop.
    m = @enzyme_mechanism begin
        substrates: S
        products: P, Q
        steps: begin
            E + S <--> EPQ
            EP + Q <--> EPQ
            EQ + P <--> EPQ
            E + P ⇌ EP
            E + Q ⇌ EQ
        end
    end
    rxn = @enzyme_reaction begin
        substrates: S[CN]
        products: P[C], Q[N]
        oligomeric_state: 2
    end
    k_pe = @_testhelper_unguarded @enzyme_mechanism begin
        substrates: S ; products: P, Q ; oligomeric_state: 2
        steps: begin
            E + S <--> EPQ       :: OnlyA
            EP + Q <--> EPQ    :: EqualAI
            EQ + P <--> EPQ    :: EqualAI
            E + P ⇌ EP             :: OnlyA
            E + Q ⇌ EQ             :: EqualAI
        end
    end
    k_qe = @_testhelper_unguarded @enzyme_mechanism begin
        substrates: S ; products: P, Q ; oligomeric_state: 2
        steps: begin
            E + S <--> EPQ       :: OnlyA
            EP + Q <--> EPQ    :: EqualAI
            EQ + P <--> EPQ    :: EqualAI
            E + P ⇌ EP             :: EqualAI
            E + Q ⇌ EQ             :: OnlyA
        end
    end
    k_q_ep = @enzyme_mechanism begin
        substrates: S ; products: P, Q ; oligomeric_state: 2
        steps: begin
            E + S <--> EPQ       :: OnlyA
            EP + Q <--> EPQ    :: OnlyA
            EQ + P <--> EPQ    :: EqualAI
            E + P ⇌ EP             :: EqualAI
            E + Q ⇌ EQ             :: EqualAI
        end
    end
    k_p_eq = @enzyme_mechanism begin
        substrates: S ; products: P, Q ; oligomeric_state: 2
        steps: begin
            E + S <--> EPQ       :: OnlyA
            EP + Q <--> EPQ    :: EqualAI
            EQ + P <--> EPQ    :: OnlyA
            E + P ⇌ EP             :: EqualAI
            E + Q ⇌ EQ             :: EqualAI
        end
    end
    k_pe_qe = @enzyme_mechanism begin
        substrates: S ; products: P, Q ; oligomeric_state: 2
        steps: begin
            E + S <--> EPQ       :: OnlyA
            EP + Q <--> EPQ    :: EqualAI
            EQ + P <--> EPQ    :: EqualAI
            E + P ⇌ EP             :: OnlyA
            E + Q ⇌ EQ             :: OnlyA
        end
    end
    k_pe_q_ep = @enzyme_mechanism begin
        substrates: S ; products: P, Q ; oligomeric_state: 2
        steps: begin
            E + S <--> EPQ       :: OnlyA
            EP + Q <--> EPQ    :: OnlyA
            EQ + P <--> EPQ    :: EqualAI
            E + P ⇌ EP             :: OnlyA
            E + Q ⇌ EQ             :: EqualAI
        end
    end
    k_pe_p_eq = @enzyme_mechanism begin
        substrates: S ; products: P, Q ; oligomeric_state: 2
        steps: begin
            E + S <--> EPQ       :: OnlyA
            EP + Q <--> EPQ    :: EqualAI
            EQ + P <--> EPQ    :: OnlyA
            E + P ⇌ EP             :: OnlyA
            E + Q ⇌ EQ             :: EqualAI
        end
    end
    k_qe_q_ep = @enzyme_mechanism begin
        substrates: S ; products: P, Q ; oligomeric_state: 2
        steps: begin
            E + S <--> EPQ       :: OnlyA
            EP + Q <--> EPQ    :: OnlyA
            EQ + P <--> EPQ    :: EqualAI
            E + P ⇌ EP             :: EqualAI
            E + Q ⇌ EQ             :: OnlyA
        end
    end
    k_qe_p_eq = @enzyme_mechanism begin
        substrates: S ; products: P, Q ; oligomeric_state: 2
        steps: begin
            E + S <--> EPQ       :: OnlyA
            EP + Q <--> EPQ    :: EqualAI
            EQ + P <--> EPQ    :: OnlyA
            E + P ⇌ EP             :: EqualAI
            E + Q ⇌ EQ             :: OnlyA
        end
    end
    k_q_ep_p_eq = @enzyme_mechanism begin
        substrates: S ; products: P, Q ; oligomeric_state: 2
        steps: begin
            E + S <--> EPQ       :: OnlyA
            EP + Q <--> EPQ    :: OnlyA
            EQ + P <--> EPQ    :: OnlyA
            E + P ⇌ EP             :: EqualAI
            E + Q ⇌ EQ             :: EqualAI
        end
    end
    k_pe_qe_q_ep = @enzyme_mechanism begin
        substrates: S ; products: P, Q ; oligomeric_state: 2
        steps: begin
            E + S <--> EPQ       :: OnlyA
            EP + Q <--> EPQ    :: OnlyA
            EQ + P <--> EPQ    :: EqualAI
            E + P ⇌ EP             :: OnlyA
            E + Q ⇌ EQ             :: OnlyA
        end
    end
    k_pe_qe_p_eq = @enzyme_mechanism begin
        substrates: S ; products: P, Q ; oligomeric_state: 2
        steps: begin
            E + S <--> EPQ       :: OnlyA
            EP + Q <--> EPQ    :: EqualAI
            EQ + P <--> EPQ    :: OnlyA
            E + P ⇌ EP             :: OnlyA
            E + Q ⇌ EQ             :: OnlyA
        end
    end
    k_pe_q_ep_p_eq = @enzyme_mechanism begin
        substrates: S ; products: P, Q ; oligomeric_state: 2
        steps: begin
            E + S <--> EPQ       :: OnlyA
            EP + Q <--> EPQ    :: OnlyA
            EQ + P <--> EPQ    :: OnlyA
            E + P ⇌ EP             :: OnlyA
            E + Q ⇌ EQ             :: EqualAI
        end
    end
    k_qe_q_ep_p_eq = @enzyme_mechanism begin
        substrates: S ; products: P, Q ; oligomeric_state: 2
        steps: begin
            E + S <--> EPQ       :: OnlyA
            EP + Q <--> EPQ    :: OnlyA
            EQ + P <--> EPQ    :: OnlyA
            E + P ⇌ EP             :: EqualAI
            E + Q ⇌ EQ             :: OnlyA
        end
    end
    k_pe_qe_q_ep_p_eq = @enzyme_mechanism begin
        substrates: S ; products: P, Q ; oligomeric_state: 2
        steps: begin
            E + S <--> EPQ       :: OnlyA
            EP + Q <--> EPQ    :: OnlyA
            EQ + P <--> EPQ    :: OnlyA
            E + P ⇌ EP             :: OnlyA
            E + Q ⇌ EQ             :: OnlyA
        end
    end
    @test EnzymeRates.Enumeration._hyperbolic_catalysis(_testhelper_candidate(m))
    kids = _testhelper_move(EnzymeRates.Enumeration._expand_to_allosteric, m, rxn)
    @test length(kids) == 13
    @test Set(kids) == Set([k_q_ep, k_p_eq, k_pe_qe, k_pe_q_ep, k_pe_p_eq, k_qe_q_ep,
                            k_qe_p_eq, k_q_ep_p_eq, k_pe_qe_q_ep, k_pe_qe_p_eq,
                            k_pe_q_ep_p_eq, k_qe_q_ep_p_eq, k_pe_qe_q_ep_p_eq])
    @test all(EnzymeRates.Enumeration._hyperbolic_where_required,
              _testhelper_candidate.(kids))
    for (dropped, step) in ((k_pe, "E + P ⇌ EP"), (k_qe, "E + Q ⇌ EQ"))
        @test dropped ∉ kids
        runs = EnzymeRates.Mechanisms.steps_running_in_allosteric_conformations(dropped)
        @test EnzymeRates.Enumeration._hyperbolic_catalysis(
            _testhelper_candidate(m), runs.active)
        @test !EnzymeRates.Enumeration._hyperbolic_catalysis(
            _testhelper_candidate(m), runs.inactive)
        v = EnzymeRates.Derivation._guard(dropped)
        @test (v.fault, EnzymeRates.Mechanisms.step_text(m.catalytic, v.edge)) ==
              (EnzymeRates.Derivation._FAST_STEP_PINNED, step)
    end
    @test Set(c for c in _testhelper_children(rxn, [m])
              if _testhelper_is_allosteric(c)) == Set(kids)
end

@testset "Mechanism — uni-uni: binding-group :OnlyA variants only" begin
    # SEED: uni-uni, no declared regulator. 3 kinetic groups: 2 bindings (S, P) + 1
    # catalytic. The catalytic group is never a bare primary :OnlyA (that V-type needs a
    # regulator, none declared), so every child is a K-type: every non-empty subset of
    # the 2 bindings (S, P) is set :OnlyA on top of the base whose chemical step is
    # :OnlyA (dead inactive), and the guard keeps all 3 (`check_valid_rate_equation`).
    m = @enzyme_mechanism begin
        substrates: S
        products: P
        steps: begin
            E + S ⇌ ES
            ES <--> EP
            E + P ⇌ EP
        end
    end
    only_s = @enzyme_mechanism begin
        substrates: S
        products: P
        oligomeric_state: 2
        steps: begin
            E + S ⇌ ES      :: OnlyA
            ES <--> EP    :: OnlyA
            E + P ⇌ EP      :: EqualAI
        end
    end
    only_p = @enzyme_mechanism begin
        substrates: S
        products: P
        oligomeric_state: 2
        steps: begin
            E + S ⇌ ES      :: EqualAI
            ES <--> EP    :: OnlyA
            E + P ⇌ EP      :: OnlyA
        end
    end
    both = @enzyme_mechanism begin
        substrates: S
        products: P
        oligomeric_state: 2
        steps: begin
            E + S ⇌ ES      :: OnlyA
            ES <--> EP    :: OnlyA
            E + P ⇌ EP      :: OnlyA
        end
    end

    result = _testhelper_move(EnzymeRates.Enumeration._expand_to_allosteric, m,
                              uni_uni_allo)
    @test length(result) == 3
    @test Set(result) == Set([only_s, only_p, both])

    for r in result
        _testhelper_assert_mechanism_invariants(r)
    end
end

@testset "Mechanism — enumerates all allowed multiplicities" begin
    # Multi-valued allowed_oligomeric_states → the binding-group
    # :OnlyA variant set is emitted once per multiplicity.
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
    # The S binding, the P binding, or both :OnlyA with the chemistry, at each
    # multiplicity.
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
            oligomeric_state: 4
            steps: begin
                E + S ⇌ ES    :: EqualAI
                ES <--> EP  :: OnlyA
                EP ⇌ E + P    :: OnlyA
            end
        end),
        (@enzyme_mechanism begin
            substrates: S
            products: P
            oligomeric_state: 4
            steps: begin
                E + S ⇌ ES    :: OnlyA
                ES <--> EP  :: OnlyA
                EP ⇌ E + P    :: EqualAI
            end
        end),
        (@enzyme_mechanism begin
            substrates: S
            products: P
            oligomeric_state: 4
            steps: begin
                E + S ⇌ ES    :: OnlyA
                ES <--> EP  :: OnlyA
                EP ⇌ E + P    :: OnlyA
            end
        end),
    ]
    @test length(allo) == 6
    @test Set(allo) == Set(expected)
    mults = Set(am.allosteric.oligomeric_state for am in allo)
    @test mults == Set([2, 4])

    # No regulator declared → the binding-completion variant set (see the
    # uni-uni case: 3 variants) is emitted once per multiplicity.
    n_groups = length(_testhelper_group_members(m.catalytic))
    n_cat = count(g -> _testhelper_is_catalytic_group(m, g), 1:n_groups)
    @test n_cat == 1
    @test length(allo) == 2 * 3

    # Each multiplicity carries the same dead-inactive allosteric state set: no
    # all-:EqualAI baseline, and 3 variants ({S}, {P}, {S, P} bindings
    # :OnlyA, each with the chemical step :OnlyA) — :OnlyA counts 2, 2, 3.
    for cn in (2, 4)
        cn_variants = filter(
            am -> am.allosteric.oligomeric_state == cn, allo)
        @test count(
            am -> all(==(:EqualAI), am.allosteric.allosteric_states),
            cn_variants) == 0
        @test count(
            am -> count(==(:OnlyA), am.allosteric.allosteric_states) == 2,
            cn_variants) == 2
        @test count(
            am -> count(==(:OnlyA), am.allosteric.allosteric_states) == 3,
            cn_variants) == 1
    end
end

@testset "V-type param counts + competitive-inhibitor exclusion" begin
    # uni_uni_allo_2reg declares two allosteric regulators (R1, R2). From a
    # uni-uni seed (RE binding, SS catalytic), _expand_to_allosteric emits:
    #   * 3 K-type variants (binding-group :OnlyA, no regulatory site) at
    #     +1 param each — the conformational constant L; each tags a nonempty
    #     set of the bindings :OnlyA on top of the base whose chemical step is
    #     :OnlyA, and the guard keeps all 3 (`check_valid_rate_equation`);
    #   * 4 V-type variants (the catalytic group :OnlyA paired with one
    #     regulator at a site) at +2 params each — L plus the regulator's K
    #     — one per (regulator, tag), regulator ∈ {R1, R2}, tag ∈
    #     {:OnlyA, :OnlyI}.
    seed = @enzyme_mechanism begin
        substrates: S
        products: P
        steps: begin
            E + S ⇌ ES
            ES <--> EP
            E + P ⇌ EP
        end
    end
    base = _testhelper_fitted(seed)
    Δ(am) = _testhelper_fitted(am) - base
    outs = _testhelper_move(EnzymeRates.Enumeration._expand_to_allosteric, seed,
                            uni_uni_allo_2reg)
    # Tags: the P binding, the S binding, the chemistry.
    ktypes_expected = [
        @enzyme_mechanism(begin
            substrates: S; products: P
            oligomeric_state: 2
            steps: begin
                E + P ⇌ EP    :: OnlyA
                E + S ⇌ ES    :: EqualAI
                ES <--> EP  :: OnlyA
            end
        end),
        @enzyme_mechanism(begin
            substrates: S; products: P
            oligomeric_state: 2
            steps: begin
                E + P ⇌ EP    :: EqualAI
                E + S ⇌ ES    :: OnlyA
                ES <--> EP  :: OnlyA
            end
        end),
        @enzyme_mechanism(begin
            substrates: S; products: P
            oligomeric_state: 2
            steps: begin
                E + P ⇌ EP    :: OnlyA
                E + S ⇌ ES    :: OnlyA
                ES <--> EP  :: OnlyA
            end
        end)]
    vtypes_expected = [
        @enzyme_mechanism(begin
            substrates: S; products: P
            oligomeric_state: 2
            steps: begin
                E + P ⇌ EP    :: EqualAI
                E + S ⇌ ES    :: EqualAI
                ES <--> EP  :: OnlyA
            end
            regulatory_site(2): R1::OnlyA
        end),
        @enzyme_mechanism(begin
            substrates: S; products: P
            oligomeric_state: 2
            steps: begin
                E + P ⇌ EP    :: EqualAI
                E + S ⇌ ES    :: EqualAI
                ES <--> EP  :: OnlyA
            end
            regulatory_site(2): R1::OnlyI
        end),
        @enzyme_mechanism(begin
            substrates: S; products: P
            oligomeric_state: 2
            steps: begin
                E + P ⇌ EP    :: EqualAI
                E + S ⇌ ES    :: EqualAI
                ES <--> EP  :: OnlyA
            end
            regulatory_site(2): R2::OnlyA
        end),
        @enzyme_mechanism(begin
            substrates: S; products: P
            oligomeric_state: 2
            steps: begin
                E + P ⇌ EP    :: EqualAI
                E + S ⇌ ES    :: EqualAI
                ES <--> EP  :: OnlyA
            end
            regulatory_site(2): R2::OnlyI
        end)]
    @test length(outs) == 7
    @test Set(outs) == Set([ktypes_expected; vtypes_expected])
    @test all(am -> Δ(am) == 2, vtypes_expected)
    @test all(am -> Δ(am) == 1, ktypes_expected)
    # Every V-type has its catalytic group :OnlyA.
    @test all(vtypes_expected) do am
        cat_group = only(g for g in 1:length(_testhelper_group_members(am.catalytic))
                         if _testhelper_is_catalytic_group(am, g))
        am.allosteric.allosteric_states[cat_group] == :OnlyA
    end

    # A regulator already bound as a competitive inhibitor in the seed is
    # NOT re-used for a V-type — only the other (free) allosteric regulator
    # is. In uni_uni_reg_and_inhibitor, R1 is the competitive inhibitor and
    # R2 the free allosteric regulator; the seed binds R1 at free E.
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
    # Tags: the P binding, the R1 binding, the S binding, the chemistry. Every binding
    # subset :OnlyA, R1's included, then R2 alone on a site.
    expected_ci = [
        @enzyme_mechanism(begin
            substrates: S; products: P
            modifiers: R1
            oligomeric_state: 2
            steps: begin
                E + P ⇌ EP    :: OnlyA
                E + R1 ⇌ ER1  :: EqualAI
                E + S ⇌ ES    :: EqualAI
                ES <--> EP  :: OnlyA
            end
        end),
        @enzyme_mechanism(begin
            substrates: S; products: P
            modifiers: R1
            oligomeric_state: 2
            steps: begin
                E + P ⇌ EP    :: EqualAI
                E + R1 ⇌ ER1  :: OnlyA
                E + S ⇌ ES    :: EqualAI
                ES <--> EP  :: OnlyA
            end
        end),
        @enzyme_mechanism(begin
            substrates: S; products: P
            modifiers: R1
            oligomeric_state: 2
            steps: begin
                E + P ⇌ EP    :: OnlyA
                E + R1 ⇌ ER1  :: OnlyA
                E + S ⇌ ES    :: EqualAI
                ES <--> EP  :: OnlyA
            end
        end),
        @enzyme_mechanism(begin
            substrates: S; products: P
            modifiers: R1
            oligomeric_state: 2
            steps: begin
                E + P ⇌ EP    :: EqualAI
                E + R1 ⇌ ER1  :: EqualAI
                E + S ⇌ ES    :: OnlyA
                ES <--> EP  :: OnlyA
            end
        end),
        @enzyme_mechanism(begin
            substrates: S; products: P
            modifiers: R1
            oligomeric_state: 2
            steps: begin
                E + P ⇌ EP    :: OnlyA
                E + R1 ⇌ ER1  :: EqualAI
                E + S ⇌ ES    :: OnlyA
                ES <--> EP  :: OnlyA
            end
        end),
        @enzyme_mechanism(begin
            substrates: S; products: P
            modifiers: R1
            oligomeric_state: 2
            steps: begin
                E + P ⇌ EP    :: EqualAI
                E + R1 ⇌ ER1  :: OnlyA
                E + S ⇌ ES    :: OnlyA
                ES <--> EP  :: OnlyA
            end
        end),
        @enzyme_mechanism(begin
            substrates: S; products: P
            modifiers: R1
            oligomeric_state: 2
            steps: begin
                E + P ⇌ EP    :: OnlyA
                E + R1 ⇌ ER1  :: OnlyA
                E + S ⇌ ES    :: OnlyA
                ES <--> EP  :: OnlyA
            end
        end),
        @enzyme_mechanism(begin
            substrates: S; products: P
            modifiers: R1
            oligomeric_state: 2
            steps: begin
                E + P ⇌ EP    :: EqualAI
                E + R1 ⇌ ER1  :: EqualAI
                E + S ⇌ ES    :: EqualAI
                ES <--> EP  :: OnlyA
            end
            regulatory_site(2): R2::OnlyA
        end),
        @enzyme_mechanism(begin
            substrates: S; products: P
            modifiers: R1
            oligomeric_state: 2
            steps: begin
                E + P ⇌ EP    :: EqualAI
                E + R1 ⇌ ER1  :: EqualAI
                E + S ⇌ ES    :: EqualAI
                ES <--> EP  :: OnlyA
            end
            regulatory_site(2): R2::OnlyI
        end)]
    @test length(kids_ci) == 9
    @test Set(kids_ci) == Set(expected_ci)
    site_regs = Set(l for am in kids_ci
                    for s in am.allosteric.regulatory_sites
                    for l in s.ligands)
    @test site_regs == Set([:R2])
end

@testset "distinguishability invariant on every emitted mechanism" begin
    # Every mechanism _expand_to_allosteric returns is distinguishable
    # from a simpler mechanism: not all-:EqualAI, and either a binding
    # group carries a non-:EqualAI tag or a regulatory site is present.
    base = @enzyme_mechanism begin
        substrates: S
        products: P
        steps: begin
            E + S ⇌ ES
            ES <--> EP
            EP ⇌ E + P
        end
    end
    # The K-types: the S binding, the P binding, or both :OnlyA with the chemistry.
    k_types = [
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
    ]
    # With R declared, also the V-types: the chemistry :OnlyA with R at a site.
    v_types = [
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
        (@enzyme_mechanism begin
            substrates: S
            products: P
            oligomeric_state: 2
            steps: begin
                E + S ⇌ ES    :: EqualAI
                ES <--> EP  :: OnlyA
                EP ⇌ E + P    :: EqualAI
            end
            regulatory_site(2): R::OnlyI
        end),
    ]
    scenarios = [(uni_uni_allo, k_types), (uni_uni_allo_reg, [k_types; v_types])]
    for (rxn, expected) in scenarios
        children = _testhelper_move(EnzymeRates.Enumeration._expand_to_allosteric, base,
                                    rxn)
        @test length(children) == length(expected)
        @test Set(children) == Set(expected)
        for am in children
            states = am.allosteric.allosteric_states
            @test !all(==(:EqualAI), states)
            binding_gs = [g for g in 1:length(states)
                          if !_testhelper_is_catalytic_group(base, g)]
            binding_non_equal = any(g -> states[g] != :EqualAI, binding_gs)
            @test binding_non_equal ||
                !isempty(am.allosteric.regulatory_sites)
        end
    end
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
    # allowed. One catalytic subunit adds no concentration power of its own, so
    # the multiplicity-1 variants are emitted exactly as for a reaction allowing
    # only multiplicity 1; multiplicity 2 would stack a second source of powers
    # and emits nothing. Three binding groups (S, P, Q), each subset :OnlyA with
    # the chemistry :OnlyA: 2^3 - 1 = 7 K-type children, all valid
    # (`check_valid_rate_equation`).
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
    children = _testhelper_move(EnzymeRates.Enumeration._expand_to_allosteric, m, rxn_1_2)
    @test length(children) == 7
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

    # Five binding groups (A, B, Q with its inhibitor-bound mirror, P, and the
    # inhibitor's own group), each subset :OnlyA with the chemistry :OnlyA:
    # 2^5 - 1 = 31 K-type children, all emitted. The reaction declares no
    # allosteric regulator, so no V-type variant is emitted. N = 31. The
    # derived denominator carries A²,
    # but from the inhibitor site, so the catalytic scheme is hyperbolic and
    # the promotion proceeds.
    N = 31
    children = _testhelper_move(EnzymeRates.Enumeration._expand_to_allosteric, parent, rxn)
    @test length(children) == N
    @test all(EnzymeRates.Enumeration._hyperbolic_catalysis,
              _testhelper_candidate.(children))
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
    # a site ligand on its site, so its only child binds I.
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
        @enzyme_mechanism begin
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
        end]
    @test length(children) == 1
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

@testset "Allosteric mechanism — regulator tag removal delta" begin
    # SEED: uni-uni allosteric with one regulator R tagged :OnlyA;
    # all 3 catalytic groups :EqualAI.
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

    result = _testhelper_move(EnzymeRates.Enumeration._expand_change_allo_state, am)

    # 1. count: 3 cat-group relaxations + 1 reg-ligand relaxation = 4.
    @test length(result) == 4

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
    # groups. Each non-:NonequalAI entry contributes one variant:
    # 3 cat-group + 1 reg-ligand = 4.
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

    result = _testhelper_move(EnzymeRates.Enumeration._expand_change_allo_state, am)

    # 1. count: 3 cat-group relaxations + 1 reg-ligand relaxation.
    @test length(result) == 4

    # 2. ground-truth Δ multiset via `parameters`. Under strict
    # `:EqualAI`, relaxing a single binding group to :NonequalAI while its
    # partners stay :EqualAI is degenerate: a thermodynamic cycle forbids
    # that affinity split, so it collapses (K_I = K_A) and adds NO fitted
    # parameter. Only the catalytic relaxation (its derived reverse absorbs
    # the cycle) and the reg-ligand relaxation each add one. So the two
    # binding relaxations contribute 0: `[0, 0, 1, 1]`. (The enumerator
    # will skip such degenerate configs in a follow-up PR.)
    @test _testhelper_param_deltas(am, result) == [0, 0, 1, 1]

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
    # SEED: allosteric uni-uni with R1::OnlyA + R2::OnlyI at the same
    # regulatory site, all 3 cat groups :EqualAI. Each non-:NonequalAI
    # entry contributes one variant: 3 cat + 2 reg-ligand = 5.
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

    result = _testhelper_move(EnzymeRates.Enumeration._expand_change_allo_state, am)

    # 1. count: 3 cat-group relaxations + 2 reg-ligand relaxations = 5.
    @test length(result) == 5

    # 2. ground-truth Δ multiset via `parameters`. Under strict
    # `:EqualAI`, relaxing a single binding group to :NonequalAI while its
    # partners stay :EqualAI collapses that affinity (K_I = K_A) and adds NO
    # fitted parameter (a thermodynamic cycle forbids the split). The two
    # binding relaxations contribute 0; the catalytic relaxation and the two
    # reg-ligand relaxations each add one: `[0, 0, 1, 1, 1]`. (The enumerator
    # will skip such degenerate configs in a follow-up PR.)
    @test _testhelper_param_deltas(am, result) == [0, 0, 1, 1, 1]

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

@testset "Allosteric mechanism — uni-uni all-:EqualAI: 3 cat relaxations" begin
    # SEED: uni-uni allosteric with all 3 catalytic groups tagged
    # :EqualAI and no regulatory sites. Each non-:NonequalAI cat-group
    # tag contributes one variant (flip to :NonequalAI).
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

    # 1. count: 3 cat-group relaxations + 0 reg-ligand relaxations, each group in turn.
    expected = [
        (@enzyme_mechanism begin
            substrates: S; products: P
            oligomeric_state: 2
            steps: begin
                E + P ⇌ EP    :: NonequalAI
                E + S ⇌ ES    :: EqualAI
                ES <--> EP  :: EqualAI
            end
        end),
        (@enzyme_mechanism begin
            substrates: S; products: P
            oligomeric_state: 2
            steps: begin
                E + P ⇌ EP    :: EqualAI
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
    @test length(result) == 3
    @test Set(result) == Set(expected)

    # 2. each variant has exactly one :NonequalAI entry in
    # allosteric_states; the rest match the seed.
    for r in result
        @test count(t -> t == :NonequalAI, r.allosteric.allosteric_states) == 1
        for i in 1:length(am.allosteric.allosteric_states)
            @test r.allosteric.allosteric_states[i] == :NonequalAI ||
                  r.allosteric.allosteric_states[i] == am.allosteric.allosteric_states[i]
        end
    end

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

@testset "Allosteric mechanism — drops a relaxation whose inactive law has a power" begin
    # The split abortive ping-pong at two subunits (testset "_hyperbolic_catalysis"): the
    # second half-reaction and both abortive bindings, A's at EP_res_plusA_minusP and
    # P's at EA, :OnlyA. Each other binding relaxes to :NonequalAI and keeps both laws
    # hyperbolic. Relaxing an abortive binding brings its complex into the inactive
    # conformation, whose law then reaches A² or P² (`_hyperbolic_catalysis`), so the move
    # drops both relaxations. The chemistry does not relax while a binding is :OnlyA. Four
    # children. The predicate reads the graph alone: the guard rejects this parent and each
    # relaxation, since the inactive conformation keeps the two sides of the
    # rapid-equilibrium step EB_res_plusA_minusP ⇌ EQ joined (`_FAST_STEP_PINNED`).
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
    kids = _testhelper_move(EnzymeRates.Enumeration._expand_change_allo_state, am)
    @test length(kids) == 4
    @test Set(kids) == Set([a_relaxed, q_relaxed, b_relaxed, p_relaxed])
    @test all(EnzymeRates.Enumeration._hyperbolic_catalysis, _testhelper_candidate.(kids))
    for d in (abortive_a_relaxed, abortive_p_relaxed)
        @test !EnzymeRates.Enumeration._hyperbolic_catalysis(_testhelper_candidate(d))
    end
    @test all(c -> EnzymeRates.Derivation._guard(c).fault ===
                   EnzymeRates.Derivation._FAST_STEP_PINNED,
              [am; kids; abortive_a_relaxed; abortive_p_relaxed])
end

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

# ─── _respects_reg_type ──────────────────────────────────────────────────
@testset "_respects_reg_type" begin

@testset "_respects_reg_type — one- and two-ligand sites" begin
    rxn = @enzyme_reaction begin
        substrates: S[C]
        products: P[C]
        allosteric_regulators: A::Activator, A2::Activator, I::Inhibitor, I2::Inhibitor, U
        oligomeric_state: 2
    end
    core = @enzyme_mechanism begin
        substrates: S
        products: P
        oligomeric_state: 2
        steps: begin
            E + S ⇌ ES   :: EqualAI
            E + P ⇌ EP   :: EqualAI
            ES <--> EP :: EqualAI
        end
    end
    # Whether `core` with one site holding `ligs` at `states` respects the types `r`
    # declares.
    _testhelper_respects(ligs, states, r = rxn) =
        EnzymeRates.Enumeration._respects_reg_type(EnzymeRates.Mechanism(
            core.catalytic, EnzymeRates.Mechanisms.Allostery(
            core.allosteric.allosteric_states, core.allosteric.oligomeric_state,
            [EnzymeRates.Mechanisms.RegulatorySite(
                collect(Symbol, ligs), 2, collect(states))])), r)

    # An :unspecified ligand passes in every state, alone or beside another
    # :unspecified ligand: U is declared untyped and Nope is not declared at all.
    for state in (:OnlyA, :OnlyI, :EqualAI, :NonequalAI)
        @test _testhelper_respects((:U,), (state,))
        @test _testhelper_respects((:Nope,), (state,))
        @test _testhelper_respects((:U, :Nope), (state, :EqualAI))
    end
    # A name declared only as a competitive inhibitor, which takes no type, is
    # :unspecified at a regulatory site.
    rxn_ci = @enzyme_reaction begin
        substrates: S[C]
        products: P[C]
        competitive_inhibitors: X
    end
    @test _testhelper_respects((:X,), (:OnlyI,), rxn_ci)
    # A designated effector is never the opposite pure state; the matching pure state
    # and :NonequalAI always pass. :EqualAI (an antagonist) passes alone and beside an
    # opposite-type or :unspecified sibling, and fails beside a same-type designated
    # sibling.
    for (ligs, states, ok) in [((:A,), (:OnlyI,), false),
                               ((:A,), (:OnlyA,), true),
                               ((:A,), (:NonequalAI,), true),
                               ((:A,), (:EqualAI,), true),
                               ((:A, :I), (:EqualAI, :OnlyI), true),
                               ((:A, :U), (:EqualAI, :OnlyA), true),
                               ((:A, :A2), (:EqualAI, :OnlyA), false),
                               ((:I,), (:OnlyA,), false),
                               ((:I,), (:OnlyI,), true),
                               ((:I,), (:NonequalAI,), true),
                               ((:I,), (:EqualAI,), true),
                               ((:I, :A), (:EqualAI, :OnlyA), true),
                               ((:I, :I2), (:EqualAI, :OnlyI), false)]
        @test _testhelper_respects(ligs, states) == ok
    end
end

@testset "_respects_reg_type as a filter" begin
    rxn = @enzyme_reaction begin
        substrates: S[C]
        products: P[C]
        allosteric_regulators: R::Activator
        oligomeric_state: 2
    end
    m_seed = @enzyme_mechanism begin
        substrates: S; products: P
        steps: begin
            E + P ⇌ EP
            E + S ⇌ ES
            ES <--> EP
        end
    end
    # The activator R as :OnlyI, against its type, and as :OnlyA.
    am_bad = @enzyme_mechanism begin
        substrates: S; products: P
        oligomeric_state: 2
        steps: begin
            E + P ⇌ EP    :: EqualAI
            E + S ⇌ ES    :: EqualAI
            ES <--> EP  :: EqualAI
        end
        regulatory_site(2): R::OnlyI
    end
    am_good = @enzyme_mechanism begin
        substrates: S; products: P
        oligomeric_state: 2
        steps: begin
            E + P ⇌ EP    :: EqualAI
            E + S ⇌ ES    :: EqualAI
            ES <--> EP  :: EqualAI
        end
        regulatory_site(2): R::OnlyA
    end

    kept = filter(c -> EnzymeRates.Enumeration._respects_reg_type(c, rxn),
        EnzymeRates.Mechanism[
            m_seed, am_bad, am_good])
    @test m_seed in kept   # Mechanism (no sites) passes trivially
    @test am_good in kept
    @test !(am_bad in kept)
    @test length(kept) == 2
end

@testset "_respects_reg_type — two-ligand site siblings" begin
    # A single site holding two ligands: each ligand's :EqualAI is judged against
    # the declared type of the other ligand at the same site.

    # DROP: two same-type activators, both :EqualAI. Each is the other's
    # same-type EqualAI sibling, so both are rejected.
    rxn_aa = @enzyme_reaction begin
        substrates: S[C]
        products: P[C]
        allosteric_regulators: A1::Activator, A2::Activator
        oligomeric_state: 2
    end
    am_aa = @enzyme_mechanism begin
        substrates: S; products: P
        oligomeric_state: 2
        steps: begin
            E + P ⇌ EP    :: EqualAI
            E + S ⇌ ES    :: EqualAI
            ES <--> EP  :: EqualAI
        end
        regulatory_site(2): A1::EqualAI, A2::EqualAI
    end
    @test !EnzymeRates.Enumeration._respects_reg_type(am_aa, rxn_aa)

    # KEEP: an :EqualAI activator beside an :OnlyI inhibitor. The
    # activator's only sibling is opposite-type, so EqualAI passes; the
    # inhibitor's :OnlyI matches its type.
    rxn_ai = @enzyme_reaction begin
        substrates: S[C]
        products: P[C]
        allosteric_regulators: A::Activator, I::Inhibitor
        oligomeric_state: 2
    end
    am_ai = @enzyme_mechanism begin
        substrates: S; products: P
        oligomeric_state: 2
        steps: begin
            E + P ⇌ EP    :: EqualAI
            E + S ⇌ ES    :: EqualAI
            ES <--> EP  :: EqualAI
        end
        regulatory_site(2): A::EqualAI, I::OnlyI
    end
    @test EnzymeRates.Enumeration._respects_reg_type(am_ai, rxn_ai)
end
end

# ═══════════════════════════════════════════════════════════════════════
# 5. Composition (dedup, candidate_mechanisms)
# ═══════════════════════════════════════════════════════════════════════

# ─── canonical by construction ─────────────────────────────────────────

@testset "Allosteric mechanism — canonical by construction" begin
# Site permutation with DISTINCT multiplicities. The constructor must permute
# allosteric_states alongside the catalytic groups and produce the same
# regulatory_sites ordering (regulatory side) regardless of input site order.
am_ab = @enzyme_mechanism begin
    substrates: S; products: P
    oligomeric_state: 2
    steps: begin
        E + P ⇌ EP    :: EqualAI
        E + S ⇌ ES    :: EqualAI
        ES <--> EP  :: EqualAI
    end
    regulatory_site(2): A::OnlyA
    regulatory_site(4): B::OnlyI
end
am_ba = @enzyme_mechanism begin   # sites swapped
    substrates: S; products: P
    oligomeric_state: 2
    steps: begin
        E + P ⇌ EP    :: EqualAI
        E + S ⇌ ES    :: EqualAI
        ES <--> EP  :: EqualAI
    end
    regulatory_site(4): B::OnlyI
    regulatory_site(2): A::OnlyA
end
_testhelper_assert_mechanism_invariants(am_ab)
_testhelper_assert_mechanism_invariants(am_ba)
@test am_ab == am_ba
end

# ─── _dedup_key ────────────────────────────────────────────────────────
@testset "_dedup_key" begin

# Mechanism dedup keys: struct equality. Mechanisms are canonical at
# construction, so dedup is just `unique!`, which relies on
# `Base.==` / `Base.hash` on the struct itself. This testset locks in the
# struct-equality contract that powers that dedup.
@testset "Allosteric mechanism — site multiplicity is part of the dedup key" begin
# Mechanism: differing site multiplicities → unequal because
# multiplicities are part of the structural identity.
am_m2 = @enzyme_mechanism begin
    substrates: S; products: P
    oligomeric_state: 2
    steps: begin
        E + P ⇌ EP    :: EqualAI
        E + S ⇌ ES    :: EqualAI
        ES <--> EP  :: EqualAI
    end
    regulatory_site(2): A::OnlyA
end
am_m4 = @enzyme_mechanism begin   # different multiplicity
    substrates: S; products: P
    oligomeric_state: 2
    steps: begin
        E + P ⇌ EP    :: EqualAI
        E + S ⇌ ES    :: EqualAI
        ES <--> EP  :: EqualAI
    end
    regulatory_site(4): A::OnlyA
end
_testhelper_assert_mechanism_invariants(am_m2)
_testhelper_assert_mechanism_invariants(am_m4)
@test am_m2 != am_m4
end
end

# ─── mechanism dedup (unique!) ──────────────────────────────────────────
@testset "Dedup" begin

@testset "Mechanism — expansion-path overlap: dedup actually fires" begin
    # Run two rounds of expansion on the uni-uni init seeds
    # (Mechanism path), then unique!. Assert that the flat vector shrinks,
    # proving that two different expansion paths — two moves, or one move
    # reached through two parents — produced equivalent Mechanisms.
    init_mechs = collect(
        _testhelper_mechanisms(EnzymeRates.Enumeration.init_mechanisms(uni_uni_rxn)))
    pool = unique!(_testhelper_children(uni_uni_rxn, init_mechs))
    expanded = _testhelper_children(uni_uni_rxn, pool)
    pre = length(expanded)
    unique!(expanded)
    # dedup fired: two different expansion paths produced equivalent
    # Mechanisms, so the flat vector shrank.
    @test length(expanded) < pre
end

@testset "Mechanism — permuted groups collapse via canonicalization" begin
    # Two Mechanisms with the same physics but with their outer
    # kinetic-group order arbitrarily rearranged should collapse to one
    # after unique!. This exercises the constructor's outer-group
    # sort with a non-trivial permutation, confirming that any
    # permutation of the outer Vector canonicalizes back to the same
    # struct.
    m_seed = @enzyme_mechanism begin
        substrates: A, B
        products: P, Q
        steps: begin
            E + A ⇌ EA
            E + Q ⇌ EQ
            (EA + B ⇌ EAB, EQ + B ⇌ EBQ)
            (EA + P ⇌ EAP, EQ + P ⇌ EPQ)
            EAB <--> EPQ
        end
    end
    _testhelper_assert_mechanism_invariants(m_seed)
    n_groups = length(_testhelper_group_members(m_seed.catalytic))
    @assert n_groups >= 3 "bi-bi init seed must have ≥3 kinetic groups " *
        "to exercise a non-trivial permutation"
    # Cyclic rotation by one, [g1, g2, …, gN] → [g2, …, gN, g1]: group g written as group
    # g - 1 and group 1 last, the steps written group by group in that order.
    label = [mod1(e.kinetic_group - 1, n_groups) for e in m_seed.catalytic.steps]
    order = sortperm(label)
    m_rotated = _testhelper_plain_mechanism(
        m_seed.catalytic.metabolites, m_seed.catalytic.enzyme_forms,
        m_seed.catalytic.substrates, m_seed.catalytic.products,
        [EnzymeRates.Mechanisms.Step(e.from, e.to, e.taken_up, e.given_off,
                                     e.rapid_equilibrium, label[k])
         for (k, e) in enumerate(m_seed.catalytic.steps)][order])
    v = EnzymeRates.Mechanism[m_seed, m_rotated]
    unique!(v)
    @test length(v) == 1
end
end

# ─── candidate_mechanisms ──────────────────────────────────────────────
@testset "candidate_mechanisms: the starting set" begin
    # The uni-uni reaction has one seed, which is not degenerate, so the starting set is
    # that seed: rapid-equilibrium bindings around a steady-state chemistry.
    seed = @enzyme_mechanism begin
        substrates: S
        products: P
        steps: begin
            E + S ⇌ ES
            ES <--> EP
            EP ⇌ E + P
        end
    end
    flipped = @enzyme_mechanism begin
        substrates: S
        products: P
        steps: begin
            E + S <--> ES
            ES <--> EP
            EP ⇌ E + P
        end
    end
    start = EnzymeRates.Enumeration.candidate_mechanisms(uni_uni_rxn)
    @test start isa Vector{EnzymeRates.Enumeration.CandidateMechanism}
    @test length(start) == 1
    @test Set(c.mechanism for c in start) == Set([seed])
    # A candidate is its mechanism to `==` and `hash`.
    @test only(start) == _testhelper_candidate(seed) != _testhelper_candidate(flipped)
    @test hash(only(start)) == hash(seed)
end

@testset "candidate_mechanisms: the size limits" begin
    # An aggregate check over the bi-bi starting set (the exception to writing every
    # fixture inline): each limit keeps, in order, exactly the candidates within it.
    D, En = EnzymeRates.Derivation, EnzymeRates.Enumeration
    every = En.candidate_mechanisms(bi_bi_rxn)
    few_constants = En.candidate_mechanisms(bi_bi_rxn; max_param_count = 6)
    @test few_constants == filter(c -> D.fitted_parameter_count(c.mechanism) <= 6, every)
    @test 0 < length(few_constants) < length(every)
    few_terms = En.candidate_mechanisms(bi_bi_rxn; eq_complexity_filter = 40)
    @test few_terms == filter(c -> D.rate_equation_term_count(c.mechanism) <= 40, every)
    @test 0 < length(few_terms) < length(every)
end

@testset "eq_complexity_filter: random bi-bi passes, ter-ter dropped (default 337)" begin
    random_bibi = @enzyme_mechanism begin        # V×τ = 336 — the ceiling; must pass
        substrates: S1, S2
        products:   P1, P2
        steps: begin
            E + S1 <--> ES1
            E + S2 <--> ES2
            ES1 + S2 <--> ES1S2
            ES2 + S1 <--> ES1S2
            ES1S2 <--> EP1P2
            EP1P2 <--> EP1 + P2
            EP1P2 <--> EP2 + P1
            EP1 <--> E + P1
            EP2 <--> E + P2
        end
    end
    # The enumerator returns a mechanism only within identify's limits
    # (`Enumeration._within_size_limits`). At the default 337 a random-order bi-bi
    # (V×τ = 336) passes, while a random-order ter-ter (V×τ ≈ 5.9M) is dropped by its term
    # count, before any derivation. A tighter limit drops the bi-bi too.
    batch = [random_bibi,
             _testhelper_random_terter()]
    within(limit) = [c.mechanism
                     for c in EnzymeRates.Enumeration._within_size_limits(
                         _testhelper_candidate.(batch), typemax(Int), limit)]
    @test within(337) == [random_bibi]
    @test isempty(within(100))
end

@testset "eq_complexity_filter counts both conformations" begin
    # The 336-term random-order bi-bi at one catalytic subunit with its chemistry and its
    # product bindings at free E :OnlyA, a K-type the tagging move emits: the inactive
    # conformation keeps the substrate square E, ES1, ES1S2, ES2, four segments and
    # four spanning trees, so the law counts 336 + 16 = 352. The default filter drops it
    # (`Enumeration._within_size_limits`, the enumerator's limits), where the active count
    # alone would let it through; the parameter limit reads its 16 fitted constants.
    ktype = @enzyme_mechanism begin
        substrates: S1, S2
        products: P1, P2
        oligomeric_state: 1
        steps: begin
            E + P1 <--> EP1 :: OnlyA
            E + P2 <--> EP2 :: OnlyA
            E + S1 <--> ES1 :: EqualAI
            E + S2 <--> ES2 :: EqualAI
            EP1 + P2 <--> EP1P2 :: EqualAI
            EP2 + P1 <--> EP1P2 :: EqualAI
            ES1 + S2 <--> ES1S2 :: EqualAI
            ES1S2 <--> EP1P2 :: OnlyA
            ES2 + S1 <--> ES1S2 :: EqualAI
        end
    end
    cm = ktype.catalytic
    runs = EnzymeRates.Mechanisms.steps_running_in_allosteric_conformations(ktype)
    @test _testhelper_laplacian_tree_count(cm, runs.active) == 336
    @test _testhelper_laplacian_tree_count(cm, runs.inactive) == 16
    @test EnzymeRates.Derivation.rate_equation_term_count(ktype) == 352
    @test EnzymeRates.Derivation.fitted_parameter_count(ktype) == 16
    within(max_param_count, eq_complexity_filter) =
        [c.mechanism for c in EnzymeRates.Enumeration._within_size_limits(
             _testhelper_candidate.([ktype]), max_param_count, eq_complexity_filter)]
    @test isempty(within(typemax(Int), 337))
    @test within(typemax(Int), 352) == [ktype]
    @test isempty(within(15, typemax(Int)))
    @test within(16, typemax(Int)) == [ktype]
end

@testset "candidate_mechanisms: a parent's children" begin
    # E + S → EP fused (SS), E + P ⇌ EP (RE). Flip: the P group cuts E from EP,
    # carries flux, no bottomless segment: one child, both steps SS. Split: single-step
    # groups, nothing. Dead end: no regulator. To-allosteric at multiplicity 1: the
    # chemistry group (the fused step) is :OnlyA, the binding subsets range over {P}:
    # one K-type child; no regulator, so no V-type.
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
    k_type = @enzyme_mechanism begin
        substrates: S
        products: P
        steps: begin
            E + S <--> EP    :: OnlyA
            E + P ⇌ EP       :: OnlyA
        end
    end
    children = EnzymeRates.Enumeration.candidate_mechanisms(
        uni_uni_rxn, _testhelper_candidate(parent))
    @test length(children) == 2
    @test Set(c.mechanism for c in children) == Set([flipped, k_type])
    # The limits drop children. Each child fits three constants. The flipped child's law
    # has 4 terms, two steady-state forms times the two steps that join them; the
    # K-type's active law has 1, one rapid-equilibrium segment, and its inactive
    # conformation runs no step.
    @test isempty(EnzymeRates.Enumeration.candidate_mechanisms(
        uni_uni_rxn, _testhelper_candidate(parent); max_param_count = 2))
    @test [c.mechanism for c in EnzymeRates.Enumeration.candidate_mechanisms(
        uni_uni_rxn, _testhelper_candidate(parent); eq_complexity_filter = 2)] == [k_type]
end

@testset "candidate_mechanisms: an undecodable reaction" begin
    # E with Ac and CoA bound and E with AcCoA bound would both be named EAcCoA, and the
    # seeds of Ac + CoA ⇌ AcCoA hold both. The reaction is a valid reaction; the
    # enumerator, which names the forms it builds, refuses it.
    rxn = @enzyme_reaction begin
        substrates: Ac[C], CoA[S]
        products: AcCoA[CS]
    end
    @test rxn isa EnzymeReaction
    @test_throws(ErrorException(
        "candidate_mechanisms: Ac + CoA and AcCoA both spell AcCoA: enzyme forms and " *
        "their constants join the names of the metabolites they bind without a " *
        "separator, so two forms would share a name; rename a metabolite"),
        EnzymeRates.Enumeration.candidate_mechanisms(rxn))
end

@testset "_required_regulators" begin
    rxn = @enzyme_reaction begin
        substrates: S[C]
        products: P[C]
        allosteric_regulators: A, B
        competitive_inhibitors: I, J
        oligomeric_state: 2
    end
    # Empty optional lists: every declared regulator is required, split by kind.
    ra, rc = EnzymeRates.Enumeration._required_regulators(rxn, Symbol[], Symbol[])
    @test ra == Set([:A, :B])
    @test rc == Set([:I, :J])
    # All regulators optional: both required sets empty (so the search starts from the
    # unregulated seeds).
    ra_opt, rc_opt =
        EnzymeRates.Enumeration._required_regulators(rxn, [:A, :B], [:I, :J])
    @test isempty(ra_opt)
    @test isempty(rc_opt)
    # Mixed: one of each kind optional, the other required.
    ra_mix, rc_mix = EnzymeRates.Enumeration._required_regulators(rxn, [:A], [:J])
    @test ra_mix == Set([:B])
    @test rc_mix == Set([:I])

    # A dual-role name (ATP declared as BOTH allosteric regulator and
    # competitive inhibitor) is required/optional per role, independently:
    # marking ATP optional as a competitive inhibitor leaves it required as an
    # allosteric regulator.
    dual_rxn = @enzyme_reaction begin
        substrates: S[C]
        products: P[C]
        allosteric_regulators: ATP
        competitive_inhibitors: ATP
        oligomeric_state: 2
    end
    ra_dual, rc_dual =
        EnzymeRates.Enumeration._required_regulators(dual_rxn, Symbol[], [:ATP])
    @test :ATP in ra_dual        # still required as an allosteric regulator
    @test !(:ATP in rc_dual)     # optional as a competitive inhibitor
end

@testset "_check_decodable_names: a reaction's names spell each name one way" begin
    # The enumerator names a form by joining the names of the metabolites it binds, so
    # `candidate_mechanisms` refuses a reaction whose names spell one string two ways; a
    # reactant also declared a competitive inhibitor binds as its copy, named Xmod.
    _testhelper_two_ways(a, b, x) = ErrorException(
        "candidate_mechanisms: $a and $b both spell $x: enzyme forms and their " *
        "constants join the names of the metabolites they bind without a separator, so " *
        "two forms would share a name; rename a metabolite")
    check(rxn) = EnzymeRates.Enumeration._check_decodable_names(rxn)
    @test_throws _testhelper_two_ways("Ac + CoA", "AcCoA", "AcCoA") check(
        @enzyme_reaction begin
            substrates: Ac[C], CoA[S]
            products:   AcCoA[CS]
        end)
    @test_throws _testhelper_two_ways("A + I", "AI", "AI") check(@enzyme_reaction begin
        substrates: A[C], AI[N]
        products:   P[CN]
        competitive_inhibitors: I
    end)
    # The split can take several steps to find: XY + Z and X + YZ both spell XYZ.
    @test_throws _testhelper_two_ways("XY + Z", "X + YZ", "XYZ") check(
        @enzyme_reaction begin
            substrates: X[C], YZ[N]
            products:   XY[C], Z[N]
        end)
    # A reactant's modifier copy is named Xmod, and the error says which name is the copy.
    @test_throws _testhelper_two_ways("Pmod (P's modifier copy)", "Pmod", "Pmod") check(
        @enzyme_reaction begin
            substrates: S[C]
            products:   P[C]
            competitive_inhibitors: P, Pmod
        end)
    @test_throws _testhelper_two_ways("X + mod", "Xmod (X's modifier copy)", "Xmod") check(
        @enzyme_reaction begin
            substrates: X[C]
            products:   mod[C]
            competitive_inhibitors: X
        end)
    # NAD and NADH spell no name two ways; with H, NAD + H spells NADH.
    @test check(@enzyme_reaction begin
        substrates: NAD[C]
        products:   NADH[C]
    end) === nothing
    @test_throws _testhelper_two_ways("NAD + H", "NADH", "NADH") check(
        @enzyme_reaction begin
            substrates: NAD[C], H[N]
            products:   NADH[CN]
        end)
    # An allosteric regulator binds no enzyme form, so its name takes no part.
    @test check(@enzyme_reaction begin
        substrates: A[C]
        products:   P[C]
        competitive_inhibitors: I
        allosteric_regulators: AI
    end) === nothing
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
    # of the bindings :OnlyA on top of the base whose chemical step is
    # :OnlyA, 3 distinct allosteric mechanism variants that the exit guard
    # (`check_valid_rate_equation`) all keeps (the catalytic group is never a bare
    # primary :OnlyA — that needs a regulator). Other moves do not
    # produce allosteric mechanism output from a plain mechanism, so
    # exactly 3 exist.
    @test allo_count == 3
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

    # Unfiltered children reproduce candidate_mechanisms' moves on a Mechanism without the
    # filter.
    raw_act = _testhelper_plain_moves(seed, rxn_act)
    raw_plain = _testhelper_plain_moves(seed, rxn_plain)

    children_act = _testhelper_children(rxn_act, [seed])
    children_plain = _testhelper_children(rxn_plain, [seed])
    # The three K-types and the V-type with R at an :OnlyA site; the typeless reaction
    # adds the V-type with R at an :OnlyI site.
    k_and_v = [
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
    @test length(children_act) == 4
    @test Set(children_act) == Set(k_and_v)
    @test length(children_plain) == 5
    @test Set(children_plain) == Set([k_and_v; only_i])

    # A designated activator's raw expansion carries an :OnlyI R child (a
    # V-type variant); the reg_type filter drops it, leaving strictly fewer.
    @test any(m -> _testhelper_has_reg_state(m, :R, :OnlyI), raw_act)
    @test !any(m -> _testhelper_has_reg_state(m, :R, :OnlyI), children_act)
    @test children_act ==
          filter(c -> EnzymeRates.Enumeration._respects_reg_type(c, rxn_act), raw_act)
    @test length(children_act) < length(raw_act)

    # A typeless reaction declares no type, so the filter is a no-op:
    # candidate_mechanisms returns exactly the raw children, :OnlyI included.
    @test any(m -> _testhelper_has_reg_state(m, :R, :OnlyI), children_plain)
    @test children_plain == raw_plain
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
    # contributes nothing: candidate_mechanisms equals the reg-type-filtered output of the
    # moves that apply to a Mechanism.
    m = @enzyme_mechanism begin
        substrates: S; products: P
        steps: begin
            E + P ⇌ EP
            E + S ⇌ ES
            ES <--> EP
        end
    end
    raw = _testhelper_plain_moves(m, uni_uni_rxn)
    baseline = filter(c -> EnzymeRates.Enumeration._respects_reg_type(c, uni_uni_rxn), raw)
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

# ═══════════════════════════════════════════════════════════════════════
# 6. Integration (enumerate_all)
# ═══════════════════════════════════════════════════════════════════════

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

# ═══════════════════════════════════════════════════════════════════════
# Testsets covering downstream concerns (canonicalization parameter-naming;
# move-on-allosteric polymorphism). Adjacent to enumeration but tests
# rate-equation-derivation and allosteric mechanism integration.
# ═══════════════════════════════════════════════════════════════════════

@testset "Tagged groups exclude inactive-state params" begin
# A K-type uni-uni: the P binding and the chemistry :OnlyA, the S binding :EqualAI.
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
# The kinetic group of `am` that holds the edge written `text`.
_testhelper_group_of(text) = am.catalytic.steps[_testhelper_edge(am, text)].kinetic_group

# The :OnlyA group g of am has no inactive-state constant (`K_I_…`/`k_I_…`, rendered
# through the naming chokepoint) among the parameters. With every :OnlyA group of am
# tagged :NonequalAI those names do appear, so the first assertion can fail.
function _testhelper_test_no_inactive_constants(am, g)
    cm = am.catalytic
    i_names = EnzymeRates.Derivation._name.(
        _testhelper_step_parameters(am, first(_testhelper_group_members(cm)[g]), true),
        Ref(am))
    @test isempty(intersect(i_names, parameters(am)))
    nonequal = EnzymeRates.Mechanism(
        cm, EnzymeRates.Mechanisms.Allostery(
                replace(am.allosteric.allosteric_states, :OnlyA => :NonequalAI),
                am.allosteric.oligomeric_state, am.allosteric.regulatory_sites))
    @test !isempty(
        intersect(i_names, parameters(nonequal)))
end

@testset ":OnlyA binding group: no K_I_ param" begin
    # Must NOT be an iso-only group (iso `:OnlyA` is just a relabel — the test wants a
    # binding group whose K param disappears in the inactive state).
    _testhelper_test_no_inactive_constants(am, _testhelper_group_of("E + P ⇌ EP"))
end

@testset ":OnlyA iso group: no k_I_ params" begin
    _testhelper_test_no_inactive_constants(am, _testhelper_group_of("ES <--> EP"))
end

@testset "zero inactive-state numerator keeps the :NonequalAI binding constants" begin
    # K-type allosteric uni-uni: catalytic step is :OnlyA (so the I-state
    # numerator is zero), but binding steps are :NonequalAI.
    # When the I-state numerator is zero, the binding partition function
    # for :NonequalAI groups must still emit their I-state constants in the
    # I-state denominator so they appear in the rate-equation body and in
    # parameters(m).
    m = @enzyme_mechanism begin
        substrates: S
        products: P
        oligomeric_state: 2
        steps: begin
            E + S ⇌ ES      :: NonequalAI
            E + P ⇌ EP      :: NonequalAI
            ES <--> EP    :: OnlyA
        end
    end
    @test iszero(EnzymeRates.Derivation._memo_law(m).inactive.Q)
    params = parameters(m)
    # K_I_ES_to_E_S and K_I_EP_to_E_P are referenced in the I-state
    # denominator of the body (the binding partition function for :NonequalAI
    # groups is built whether or not the I-state numerator is zero, since the
    # I-state denominator always appears in the rate's denominator).
    @test :K_I_ES_to_E_S in params
    @test :K_I_EP_to_E_P in params
    _testhelper_check_law_against_oracle(m)
end
end


# ─── Expansion-move helpers ──────────────────────────────────────────────
# Helpers and building blocks of the RE→SS flip move and the kinetic-group
# split move. Every fixture is written out with the mechanism macro so the
# mechanism under test is visible in the testset.

@testset "fitted_parameter_count matches the fitted parameters on the test specs" begin
    for spec in MECHANISM_TEST_SPECS
        m = spec.mechanism
        @test EnzymeRates.Derivation.fitted_parameter_count(m) ==
              _testhelper_fitted(m)
    end
end

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

@testset "MERGE, ELIM and the seed tests on the ordered bi-bi seed" begin
    ER = EnzymeRates
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
    # Each name's uptake sign: the names are A, B, P, Q.
    @test ER.Enumeration._name_signs(seed.catalytic) == [1, 1, -1, -1]
    iso = only(filter(EnzymeRates.Enumeration._is_iso, seed.catalytic.steps))
    merged = ER.Enumeration._merge_isomerization(seed.catalytic.steps, iso)
    # The isomerization's product side EPQ takes over EAB's steps.
    base = _testhelper_base(seed, merged)
    # The expected groups turn the reaction over at rapid equilibrium, which the guard
    # rejects; the constructor builds them unchecked.
    @test base == @_testhelper_unguarded(@enzyme_mechanism begin
        substrates: A, B; products: P, Q
        steps: begin
            E + A ⇌ EA
            EA + B ⇌ EPQ
            EQ + P ⇌ EPQ
            E + Q ⇌ EQ
        end
    end)
    # The merged cycle E → EA → EPQ → EQ → E takes up A and B and gives off P and
    # Q, weight 1 + 1 + 1 + 1 = 4, all at rapid equilibrium. The seed's cycle runs
    # through its steady-state isomerization, and its RE steps form a tree.
    # every step RE: infinite rate
    @test ER.Mechanisms.has_turnover_cycle(base.catalytic,
                                           _testhelper_flags(base))
    @test !ER.Mechanisms.has_turnover_cycle(seed.catalytic,
                                            _testhelper_flags(seed))
    _testhelper_variant(ss) = .!BitVector([_testhelper_binds(base, ss...)(e)
                                           for e in base.catalytic.steps])
    _testhelper_clean(m, re) = !ER.Mechanisms.has_turnover_cycle(m.catalytic, re) &&
        ER.Enumeration._has_vmax(m.catalytic, re, m.catalytic.substrates) &&
        ER.Enumeration._has_vmax(m.catalytic, re, m.catalytic.products) &&
        !ER.Enumeration._chemistry_equilibrates_both_sides(_testhelper_candidate(m), re)
    # Track 4 table 4.3: {A, Pˣ}, {Bˣ, Q}, {Bˣ, Pˣ} clean; {A, Q} fails C (B is not
    # needed); {A, Bˣ} has no forward Vmax; {Pˣ, Q} no reverse Vmax; singletons none.
    # The merged groups are A (E + A ⇌ EA), Bˣ (EA + B → EPQ), Pˣ
    # (EQ + P → EPQ) and Q (E + Q ⇌ EQ); ˣ marks a step at the merged complex
    # EPQ, the one chemistry node. V for substrates turns A and Bˣ RE and fails
    # exactly when Pˣ and Q are RE too, closing the cycle: the steady-state set lies in
    # {A, Bˣ}. V for products fails exactly when it lies in {Pˣ, Q}. Every singleton
    # lies in one of the two. C fails exactly when Bˣ (giving B off from EPQ) and Pˣ
    # (giving P off) both stay RE: among the pairs, {A, Q} only.
    @test _testhelper_clean(base, _testhelper_variant([:A, :P])) &&
          _testhelper_clean(base, _testhelper_variant([:B, :Q])) &&
          _testhelper_clean(base, _testhelper_variant([:B, :P]))
    @test ER.Enumeration._chemistry_equilibrates_both_sides(_testhelper_candidate(base),
                                                            _testhelper_variant([:A, :Q]))
    @test !ER.Enumeration._has_vmax(base.catalytic, _testhelper_variant([:A, :B]),
                                    base.catalytic.substrates)
    @test !ER.Enumeration._has_vmax(base.catalytic, _testhelper_variant([:P, :Q]),
                                    base.catalytic.products)
    @test all(x -> !_testhelper_clean(base, _testhelper_variant([x])), (:A, :B, :P, :Q))
    # ELIM of the merged complex: one step EA + B → EQ + P, written from the binding
    # of a substrate, in a group of its own after the others; all three SS is clean.
    # The seed's form holding `xs`, by its record; the metabolite named `x`.
    at(xs...) = only(findall(f -> _testhelper_bound(f) == [xs...],
                             _testhelper_candidate(seed).forms))
    metabolite(x) = _testhelper_metabolite(seed.catalytic, x)
    tc_steps = ER.Enumeration._eliminate_form(merged, iso.to, seed.catalytic.substrates)
    # Labelled 6, after the merged groups 1, 2, 3 and 5 (the isomerization's group 4
    # went with it).
    @test last(tc_steps) ==
          ER.Mechanisms.Step(at(:A), at(:Q), [metabolite(:B)], [metabolite(:P)], true, 6)
    @test [e.kinetic_group for e in tc_steps] == [1, 2, 6]
    tc = _testhelper_base(seed, tc_steps)
    @test tc == @_testhelper_unguarded(@enzyme_mechanism begin
        substrates: A, B; products: P, Q
        steps: begin
            E + A ⇌ EA
            EA + B ⇌ EQ + P
            E + Q ⇌ EQ
        end
    end)
    group_of(pred) =
        tc.catalytic.steps[only(findall(pred, tc.catalytic.steps))].kinetic_group
    gA, gQ = group_of(_testhelper_binds(tc, :A)), group_of(_testhelper_binds(tc, :Q))
    gTC = group_of(_testhelper_is_tc)
    # E → EA → EQ → E, all RE
    @test ER.Mechanisms.has_turnover_cycle(tc.catalytic, _testhelper_flags(tc))
    _testhelper_ss(which) =
        BitVector([!(e.kinetic_group in which) for e in tc.catalytic.steps])
    # The Theorell–Chance step takes up B and gives off P, so V turns it RE for either
    # side. {A, Q, TC}: V for substrates leaves Q steady state, V for products leaves A,
    # and no RE step is left for C. {A, Q}: the RE Theorell–Chance step gives B off one
    # way and P the other, so C fails. {A, TC}: V for substrates closes the cycle with
    # the RE Q step, no forward Vmax; {Q, TC} is the mirror.
    @test _testhelper_clean(tc, _testhelper_ss([gA, gQ, gTC]))
    @test !_testhelper_clean(tc, _testhelper_ss([gA, gQ])) &&
          !_testhelper_clean(tc, _testhelper_ss([gA, gTC])) &&
          !_testhelper_clean(tc, _testhelper_ss([gQ, gTC]))
    @test ER.Enumeration._chemistry_equilibrates_both_sides(_testhelper_candidate(tc),
                                                            _testhelper_ss([gA, gQ]))
    @test !ER.Enumeration._has_vmax(tc.catalytic, _testhelper_ss([gA, gTC]),
                                    tc.catalytic.substrates)
    @test !ER.Enumeration._has_vmax(tc.catalytic, _testhelper_ss([gQ, gTC]),
                                    tc.catalytic.products)
    # ELIM needs two bindings into the form: EA's second step binds B out of it.
    @test ER.Enumeration._eliminate_form(merged, at(:A),
                                         seed.catalytic.substrates) === nothing
    # The seed itself is not degenerate; the probe agrees on every verdict above.
    @test !ER.Enumeration._degenerate(_testhelper_candidate(seed))
    @test !_testhelper_degenerate(seed)
    pairs = [(:A, :B), (:A, :P), (:A, :Q), (:B, :P), (:B, :Q), (:P, :Q)]
    candidates = vcat(
        [(string(x), base, _testhelper_variant([x])) for x in (:A, :B, :P, :Q)],
        [(string(x, y), base, _testhelper_variant([x, y])) for (x, y) in pairs],
        [("TC: A Q TC", tc, _testhelper_ss([gA, gQ, gTC])),
         ("TC: A Q", tc, _testhelper_ss([gA, gQ])),
         ("TC: A TC", tc, _testhelper_ss([gA, gTC])),
         ("TC: Q TC", tc, _testhelper_ss([gQ, gTC]))])
    for (label, g, re) in candidates
        m = ER.Enumeration._with_flags(_testhelper_candidate(g), re).mechanism
        @test (label, ER.Enumeration._degenerate(_testhelper_candidate(m))) == (
            label, !_testhelper_clean(g, re))
        @test (label, _testhelper_degenerate(m)) == (label, !_testhelper_clean(g, re))
    end
end

@testset "a reactant's modifier copy counts as its reactant in the seed tests" begin
    ER = EnzymeRates
    # Ordered bi-bi with A also bound as a modifier copy at EQ, where it does not react.
    # That binding takes A up, so it takes A's uptake sign and exchanges a substrate: the
    # metabolites are A, B, P, Q.
    m = @enzyme_mechanism begin
        substrates: A, B; products: P, Q
        steps: begin
            E + A ⇌ EA
            EA + B ⇌ EAB
            EAB <--> EPQ
            EQ + P ⇌ EPQ
            E + Q ⇌ EQ
            EQ + A ⇌ EAmodQ
        end
    end
    @test m.catalytic.metabolites == [:A, :B, :P, :Q]
    @test ER.Enumeration._name_signs(m.catalytic) == [1, 1, -1, -1]
    copy_binding = _testhelper_edge(m, "EQ + A ⇌ EAmodQ")
    @test ER.Enumeration._uptake_weight(m.catalytic, copy_binding) == 1
    @test ER.Enumeration._exchanges(m.catalytic, copy_binding, m.catalytic.substrates)
    @test !ER.Enumeration._exchanges(m.catalytic, copy_binding, m.catalytic.products)
end

@testset "ELIM refuses the uni-uni merged complex" begin
    ER = EnzymeRates
    seed = @enzyme_mechanism begin
        substrates: S; products: P
        steps: begin
            E + S ⇌ ES
            ES <--> EP
            E + P ⇌ EP
        end
    end
    iso = only(filter(EnzymeRates.Enumeration._is_iso, seed.catalytic.steps))
    merged = ER.Enumeration._merge_isomerization(seed.catalytic.steps, iso)
    # S binds E into the merged complex EP, and P binds E into it as well.
    # The expected groups turn the reaction over at rapid equilibrium, which the guard
    # rejects; the constructor builds them unchecked.
    @test _testhelper_base(seed, merged) == @_testhelper_unguarded(@enzyme_mechanism begin
        substrates: S; products: P
        steps: begin
            E + S ⇌ EP
            E + P ⇌ EP
        end
    end)
    # Both of EP's steps start at E, so a Theorell–Chance step E + S → E + P would join
    # E to itself.
    @test ER.Enumeration._eliminate_form(merged, iso.to,
                                         seed.catalytic.substrates) === nothing
end

@testset "has_turnover_cycle: a balanced rapid-equilibrium cycle runs no turnover" begin
    ER = EnzymeRates
    # Random-order binding of A and B at rapid equilibrium closes the cycle
    # E → EA → EAB ← EB ← E of weight 1 + 1 − 1 − 1 = 0: it runs no net reaction.
    # With the isomerization RE as well, E → EA → EAB → EPQ → EQ → E has
    # weight 4, one turnover, in the same block.
    m = @enzyme_mechanism begin
        substrates: A, B; products: P, Q
        steps: begin
            E + A ⇌ EA
            E + B ⇌ EB
            EA + B ⇌ EAB
            EB + A ⇌ EAB
            EAB <--> EPQ
            EQ + P ⇌ EPQ
            E + Q ⇌ EQ
        end
    end
    @test !ER.Mechanisms.has_turnover_cycle(m.catalytic, _testhelper_flags(m))
    @test ER.Mechanisms.has_turnover_cycle(m.catalytic, trues(length(m.catalytic.steps)))
end

@testset "C on the ping-pong seed and its single flips" begin
    ER = EnzymeRates
    seed = @enzyme_mechanism begin
        substrates: A, B; products: P, Q
        steps: begin
            E + A ⇌ EA
            EA <--> EP_res_plusA_minusP
            E_res_plusA_minusP + P ⇌ EP_res_plusA_minusP
            E_res_plusA_minusP + B ⇌ EB_res_plusA_minusP
            EB_res_plusA_minusP ⇌ EQ
            E + Q ⇌ EQ
        end
    end
    # The steady-state isomerization EA → E(P; res) cuts the one catalytic cycle, and
    # no step it would turn RE touches a substrate or a product, so both maximal rates
    # exist. The RE isomerization E(B; res) ⇌ EQ makes {E(B; res), EQ} a chemistry
    # node, left by the RE B binding (traversed away from E(B; res) it gives B off) and
    # by the RE Q binding (away from EQ it gives Q off): C fails. At zero products the
    # forms past the steady-state step drain to E and EA, so B is not needed.
    re = _testhelper_flags(seed)
    @test ER.Enumeration._has_vmax(seed.catalytic, re, seed.catalytic.substrates) &&
          ER.Enumeration._has_vmax(seed.catalytic, re, seed.catalytic.products)
    @test ER.Enumeration._chemistry_equilibrates_both_sides(_testhelper_candidate(seed), re)
    @test ER.Enumeration._degenerate(_testhelper_candidate(seed))
    @test _testhelper_degenerate(seed)
    # Flipping B or Q leaves the node one RE exit, Q's or B's: C holds (the B flip is the
    # steady-state half of the SRR ping-pong). Flipping A or P leaves both exits RE: C
    # still fails. V never turns the steady-state isomerization RE (it touches no
    # metabolite), so it still cuts the cycle and every flip keeps both maximal rates.
    for x in (:A, :B, :P, :Q)
        re = _testhelper_flip_flags(seed, _testhelper_binds(seed, x))
        m = ER.Enumeration._with_flags(_testhelper_candidate(seed), re).mechanism
        fails_c = x in (:A, :P)
        @test ER.Enumeration._has_vmax(seed.catalytic, re, seed.catalytic.substrates) &&
              ER.Enumeration._has_vmax(seed.catalytic, re, seed.catalytic.products)
        @test (x, ER.Enumeration._chemistry_equilibrates_both_sides(
            _testhelper_candidate(seed), re)) ==
              (x, fails_c)
        @test (x, ER.Enumeration._degenerate(_testhelper_candidate(m))) == (x, fails_c)
        @test (x, _testhelper_degenerate(m)) == (x, fails_c)
    end
    # A conformational variant is read on its active-state graph, the seed's: still
    # degenerate. The variant tags both half-reactions :OnlyA, the tagging move's one base
    # for this seed. Tagging the second half alone would leave the two sides of its
    # rapid-equilibrium isomerization joined in the inactive conformation, which the guard
    # rejects (`_FAST_STEP_PINNED`).
    am = @enzyme_mechanism begin
        substrates: A, B ; products: P, Q ; oligomeric_state: 1
        steps: begin
            E + A ⇌ EA                                          :: EqualAI
            EA <--> EP_res_plusA_minusP                      :: OnlyA
            E_res_plusA_minusP + P ⇌ EP_res_plusA_minusP    :: EqualAI
            E_res_plusA_minusP + B ⇌ EB_res_plusA_minusP    :: EqualAI
            EB_res_plusA_minusP ⇌ EQ                         :: OnlyA
            E + Q ⇌ EQ                                          :: EqualAI
        end
    end
    @test all(_testhelper_admission_agreement(am; rng = Random.MersenneTwister(9624),
                                              corners = true))
    _testhelper_check_law_against_limit(am)
    @test ER.Enumeration._degenerate(_testhelper_candidate(am))
    # The probe reads the derived allosteric law, whose inactive conformation runs no
    # chemistry and keeps only A's and Q's bindings at free E: it calls the law degenerate,
    # as the seed's.
    @test _testhelper_degenerate(am) == ER.Enumeration._degenerate(
        _testhelper_candidate(am))
end

@testset "_seed_variants: the uni-uni seed has none" begin
    ER = EnzymeRates
    seed = @enzyme_mechanism begin
        substrates: S; products: P
        steps: begin
            E + S ⇌ ES
            ES <--> EP
            E + P ⇌ EP
        end
    end
    # The merged base holds S (E + S → EP, fused) and P (E + P ⇌ EP), both units.
    # Either group alone at steady state has no maximal rate one way: V turns it rapid
    # equilibrium and closes the turnover cycle E → EP → E with the other. {S, P} passes
    # every test but is the lumping twin of the seed (both steps at EP steady state, each
    # alone in its group), so it is skipped. Both steps at EP start at E, so ELIM gives no
    # Theorell–Chance base.
    @test isempty(_testhelper_move(ER.Enumeration._seed_variants, seed, uni_uni_rxn))
    only_s = @enzyme_mechanism begin
        substrates: S; products: P
        steps: begin
            E + S <--> EP
            E + P ⇌ EP
        end
    end
    only_p = @enzyme_mechanism begin
        substrates: S; products: P
        steps: begin
            E + S ⇌ EP
            E + P <--> EP
        end
    end
    twin = @enzyme_mechanism begin
        substrates: S; products: P
        steps: begin
            E + S <--> EP
            E + P <--> EP
        end
    end
    @test !ER.Enumeration._has_vmax(only_s.catalytic, _testhelper_flags(only_s),
                                    only_s.catalytic.substrates)
    @test !ER.Enumeration._has_vmax(only_p.catalytic, _testhelper_flags(only_p),
                                    only_p.catalytic.products)
    @test _testhelper_degenerate(only_s) && _testhelper_degenerate(only_p)
    # The twin is skipped as the seed's twin, not as degenerate.
    @test !ER.Enumeration._degenerate(
        _testhelper_candidate(twin)) && !_testhelper_degenerate(twin)
end

@testset "_seed_variants: the undecorated ordered seed" begin
    ER = EnzymeRates.Enumeration
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
    # The merged base holds A (E + A ⇌ EA), Bˣ (EA + B → EPQ), Pˣ
    # (EQ + P → EPQ) and Q (E + Q ⇌ EQ), every group a unit. Every singleton has no
    # maximal rate one way; among the pairs {A, Q} fails C, {A, Bˣ} has no forward and
    # {Pˣ, Q} no reverse maximal rate (verdicts and probes pinned in "MERGE, ELIM and the
    # seed tests on the ordered bi-bi seed"). {A, Pˣ} and {Bˣ, Q} are variants; {Bˣ, Pˣ} is
    # the seed's lumping twin and is skipped; every triple holds a passing pair. The
    # Theorell–Chance base (A, Q, TC = EA + B → EQ + P) fails every set of one or two
    # groups and gives {A, Q, TC}. Each variant fits two SS groups and two RE groups
    # (2 + 2 + 1 + 1) or three SS groups (2 + 2 + 2), less the Haldane relation: 5.
    a_p = @enzyme_mechanism begin
        substrates: A, B; products: P, Q
        steps: begin
            E + A <--> EA
            EA + B ⇌ EPQ
            EQ + P <--> EPQ
            E + Q ⇌ EQ
        end
    end
    b_q = @enzyme_mechanism begin
        substrates: A, B; products: P, Q
        steps: begin
            E + A ⇌ EA
            EA + B <--> EPQ
            EQ + P ⇌ EPQ
            E + Q <--> EQ
        end
    end
    tc = @enzyme_mechanism begin
        substrates: A, B; products: P, Q
        steps: begin
            E + A <--> EA
            EA + B <--> EQ + P
            E + Q <--> EQ
        end
    end
    variants = _testhelper_move(ER._seed_variants, seed, bi_bi_rxn)
    @test length(variants) == 3
    @test Set(variants) == Set([a_p, b_q, tc])
    for v in variants
        fitted = _testhelper_fitted(v)
        @test (fitted, _testhelper_identifiable_rank(v)) == (5, 5)
        @test !_testhelper_degenerate(v)
    end
end

@testset "_seed_variants: an ordered seed whose B group is shared" begin
    ER = EnzymeRates.Enumeration
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
    # The merged base holds A, Bˣ (EA + B → EPQ with the abortive EQ + B ⇌ EBQ),
    # Pˣ and Q; the abortive step is a dead end, and its group carries flux through the
    # fused step, so every group is a unit. The verdicts are the undecorated seed's: the
    # abortive step lies on no cycle, so it changes no turnover test, and it leaves no node
    # (EPQ is the one node, with exits Bˣ and Pˣ). Singletons fail V; {A, Q} fails C;
    # {A, Bˣ} and {Pˣ, Q} fail V. {A, Pˣ}, {Bˣ, Q} and {Bˣ, Pˣ} are variants: {Bˣ, Pˣ} is no
    # lumping twin, because the Bˣ group holds two steps. Each fits 2 + 2 + 1 + 1 − 1 = 5.
    # ELIM of EPQ leaves the abortive step alone in its group, which carries no flux
    # in the all-steady-state base, so it is no unit and stays RE. As in the undecorated
    # base, every set of one or two of A, Q and TC fails, and {A, Q, TC} fits
    # 2 + 2 + 2 + 1 − 1 = 6.
    b_p = @enzyme_mechanism begin
        substrates: A, B; products: P, Q
        steps: begin
            E + A ⇌ EA
            (EA + B <--> EPQ, EQ + B <--> EBQ)
            EQ + P <--> EPQ
            E + Q ⇌ EQ
        end
    end
    a_p = @enzyme_mechanism begin
        substrates: A, B; products: P, Q
        steps: begin
            E + A <--> EA
            (EA + B ⇌ EPQ, EQ + B ⇌ EBQ)
            EQ + P <--> EPQ
            E + Q ⇌ EQ
        end
    end
    b_q = @enzyme_mechanism begin
        substrates: A, B; products: P, Q
        steps: begin
            E + A ⇌ EA
            (EA + B <--> EPQ, EQ + B <--> EBQ)
            EQ + P ⇌ EPQ
            E + Q <--> EQ
        end
    end
    tc = @enzyme_mechanism begin
        substrates: A, B; products: P, Q
        steps: begin
            E + A <--> EA
            EQ + B ⇌ EBQ
            EA + B <--> EQ + P
            E + Q <--> EQ
        end
    end
    variants = _testhelper_move(ER._seed_variants, seed, bi_bi_rxn)
    @test length(variants) == 4
    @test Set(variants) == Set([b_p, a_p, b_q, tc])
    for (v, n) in ((b_p, 5), (a_p, 5), (b_q, 5), (tc, 6))
        fitted = _testhelper_fitted(v)
        @test (fitted, _testhelper_identifiable_rank(v)) == (n, n)
        @test !_testhelper_degenerate(v)
    end
    # The rejected candidates, each degenerate by the probe as well.
    base, tc_base = _testhelper_bases(seed)
    bA, bQ = _testhelper_binds(tc_base, :A), _testhelper_binds(tc_base, :Q)
    rejected = vcat(
        [_testhelper_flip_matching(base, _testhelper_binds(base, x))
         for x in (:A, :B, :P, :Q)],
        [_testhelper_flip_matching(base, _testhelper_binds(base, xs...))
         for xs in ([:A, :Q], [:A, :B], [:P, :Q])],
        [_testhelper_flip_matching(tc_base, ps...)
         for ps in ([bA], [bQ], [_testhelper_is_tc], [bA, bQ],
                    [bA, _testhelper_is_tc], [bQ, _testhelper_is_tc])])
    for m in rejected
        @test ER._degenerate(_testhelper_candidate(m)) && _testhelper_degenerate(m)
    end
end

@testset "_seed_variants: the undecorated ping-pong seed" begin
    ER = EnzymeRates
    seed = @enzyme_mechanism begin
        substrates: A, B; products: P, Q
        steps: begin
            E + A ⇌ EA
            EA <--> EP_res_plusA_minusP
            E_res_plusA_minusP + P ⇌ EP_res_plusA_minusP
            E_res_plusA_minusP + B ⇌ EB_res_plusA_minusP
            EB_res_plusA_minusP ⇌ EQ
            E + Q ⇌ EQ
        end
    end
    # Merging both isomerizations gives A (E + A → E(P; res), fused), P, B
    # (E(; res) + B → EQ, fused) and Q, all units, with the merged complexes E(P; res) and
    # EQ as C's nodes. Every singleton fails V. {A, Q} (PP-AQ) and {B, P} (PP-BP) pass:
    # each node keeps one RE exit (P from E(P; res) and B from EQ, or A and Q), and each
    # side's steps turned RE leave the other side's steady-state group on the cycle.
    # {A, P} fails C at EQ (RE exits B and Q), {B, Q} at E(P; res) (A and P); {A, B} has
    # no forward and {P, Q} no reverse maximal rate. Every triple holds PP-AQ or PP-BP. No
    # merged complex has both its steps steady state in a variant, so neither is a twin.
    # ELIM of E(P; res) gives TC1 = E + A → E(; res) + P beside B and Q; of EQ,
    # TC2 = E(; res) + B → E + Q beside A and P. In each Theorell–Chance base a set of one
    # or two groups leaves the TC step RE (C fails) or turns both of one side's groups RE
    # (V fails); all three steady state pass. Each variant fits 2 + 2 + 1 + 1 − 1 or
    # 2 + 2 + 2 − 1 = 5.
    pp_aq = @enzyme_mechanism begin
        substrates: A, B; products: P, Q
        steps: begin
            E + A <--> EP_res_plusA_minusP
            E_res_plusA_minusP + P ⇌ EP_res_plusA_minusP
            E_res_plusA_minusP + B ⇌ EQ
            E + Q <--> EQ
        end
    end
    pp_bp = @enzyme_mechanism begin
        substrates: A, B; products: P, Q
        steps: begin
            E + A ⇌ EP_res_plusA_minusP
            E_res_plusA_minusP + P <--> EP_res_plusA_minusP
            E_res_plusA_minusP + B <--> EQ
            E + Q ⇌ EQ
        end
    end
    half_tc1 = @enzyme_mechanism begin
        substrates: A, B; products: P, Q
        steps: begin
            E + A <--> E_res_plusA_minusP + P
            E_res_plusA_minusP + B <--> EQ
            E + Q <--> EQ
        end
    end
    half_tc2 = @enzyme_mechanism begin
        substrates: A, B; products: P, Q
        steps: begin
            E + A <--> EP_res_plusA_minusP
            E_res_plusA_minusP + P <--> EP_res_plusA_minusP
            E_res_plusA_minusP + B <--> E + Q
        end
    end
    variants = _testhelper_move(ER.Enumeration._seed_variants, seed, bi_bi_pp_rxn)
    @test length(variants) == 4
    @test Set(variants) == Set([pp_aq, pp_bp, half_tc1, half_tc2])
    for v in variants
        fitted = _testhelper_fitted(v)
        @test (fitted, _testhelper_identifiable_rank(v)) == (5, 5)
        @test !_testhelper_degenerate(v)
    end
    # The rejected candidates, each degenerate by the probe as well.
    base, tc1, tc2 = _testhelper_bases(seed)
    # A Theorell–Chance base holds its TC step and the bindings of the other half's two
    # metabolites, x and y.
    function _testhelper_tc_rejected(g)
        x, y = [g.catalytic.metabolites[only(e.taken_up)] for e in g.catalytic.steps
                if EnzymeRates.Enumeration._is_binding(e)]
        bx, by = _testhelper_binds(g, x), _testhelper_binds(g, y)
        [_testhelper_flip_matching(g, ps...)
         for ps in ([_testhelper_is_tc], [bx], [by], [_testhelper_is_tc, bx],
                    [_testhelper_is_tc, by], [bx, by])]
    end
    rejected = vcat(
        [_testhelper_flip_matching(base, _testhelper_binds(base, x))
         for x in (:A, :B, :P, :Q)],
        [_testhelper_flip_matching(base, _testhelper_binds(base, xs...))
         for xs in ([:A, :P], [:B, :Q], [:A, :B], [:P, :Q])],
        _testhelper_tc_rejected(tc1), _testhelper_tc_rejected(tc2))
    for m in rejected
        @test ER.Enumeration._degenerate(
            _testhelper_candidate(m)) && _testhelper_degenerate(m)
    end
end

@testset "_seed_variants: the ordered/random seed" begin
    ER = EnzymeRates
    seed = @enzyme_mechanism begin
        substrates: A, B; products: P, Q
        steps: begin
            E + A ⇌ EA
            EA + B ⇌ EAB
            EAB <--> EPQ
            (EQ + P ⇌ EPQ, E + P ⇌ EP)
            (EP + Q ⇌ EPQ, E + Q ⇌ EQ)
        end
    end
    # The merged base holds A, Bˣ (EA + B → EPQ), Pˣ (EQ + P → EPQ,
    # E + P ⇌ EP) and Qˣ (EP + Q → EPQ, E + Q ⇌ EQ), all units. The node EPQ
    # has three exits: Bˣ gives B off, Pˣ's step gives P off and Qˣ's gives Q off. Every
    # singleton fails V.
    # {A, Pˣ} keeps the RE exits Bˣ and Qˣ, {A, Qˣ} keeps Bˣ and Pˣ: both fail C. {A, Bˣ}
    # turned RE with Pˣ and Qˣ closes a turnover cycle (no forward maximal rate), and
    # {Pˣ, Qˣ} likewise backward. {Bˣ, Pˣ} passes: V for products turns Pˣ and Qˣ RE, which
    # closes only the balanced square E → EP → EPQ ← EQ ← E, and the node keeps one
    # RE exit, Q; {Bˣ, Qˣ} is its mirror. EPQ has three steps, so neither is a twin.
    # Among the triples only {A, Pˣ, Qˣ} has three failing pairs, and it passes: the node
    # keeps the one RE exit Bˣ. Fitted: 1 + 2 + 2 + 1 − 1 = 5 for the pairs, and
    # 2 + 1 + 2 + 2 − 1 = 6 for the triple; the shared P and Q groups satisfy the square's
    # Wegscheider relation by themselves.
    b_p = @enzyme_mechanism begin
        substrates: A, B; products: P, Q
        steps: begin
            E + A ⇌ EA
            EA + B <--> EPQ
            (EQ + P <--> EPQ, E + P <--> EP)
            (EP + Q ⇌ EPQ, E + Q ⇌ EQ)
        end
    end
    b_q = @enzyme_mechanism begin
        substrates: A, B; products: P, Q
        steps: begin
            E + A ⇌ EA
            EA + B <--> EPQ
            (EQ + P ⇌ EPQ, E + P ⇌ EP)
            (EP + Q <--> EPQ, E + Q <--> EQ)
        end
    end
    a_p_q = @enzyme_mechanism begin
        substrates: A, B; products: P, Q
        steps: begin
            E + A <--> EA
            EA + B ⇌ EPQ
            (EQ + P <--> EPQ, E + P <--> EP)
            (EP + Q <--> EPQ, E + Q <--> EQ)
        end
    end
    variants = _testhelper_move(ER.Enumeration._seed_variants, seed, bi_bi_rxn)
    @test length(variants) == 3
    @test Set(variants) == Set([b_p, b_q, a_p_q])
    for (v, n) in ((b_p, 5), (b_q, 5), (a_p_q, 6))
        fitted = _testhelper_fitted(v)
        @test (fitted, _testhelper_identifiable_rank(v)) == (n, n)
        @test !_testhelper_degenerate(v)
    end
    # ELIM refuses the merged complex: it has three steps.
    iso = only(filter(EnzymeRates.Enumeration._is_iso, seed.catalytic.steps))
    merged = ER.Enumeration._merge_isomerization(seed.catalytic.steps, iso)
    @test count(s -> s.to == iso.to, merged) == 3
    @test ER.Enumeration._eliminate_form(merged, iso.to,
                                         seed.catalytic.substrates) === nothing
    # The merged base, the only one: the isomerization merged onto EPQ, every step RE.
    base = only(_testhelper_bases(seed))
    # The expected groups turn the reaction over at rapid equilibrium, which the guard
    # rejects; the constructor builds them unchecked.
    @test base == @_testhelper_unguarded(@enzyme_mechanism begin
        substrates: A, B; products: P, Q
        steps: begin
            E + A ⇌ EA
            EA + B ⇌ EPQ
            (EQ + P ⇌ EPQ, E + P ⇌ EP)
            (EP + Q ⇌ EPQ, E + Q ⇌ EQ)
        end
    end)
    # The rejected candidates, each degenerate by the probe as well.
    for xs in ([:A], [:B], [:P], [:Q], [:A, :P], [:A, :Q], [:A, :B], [:P, :Q])
        m = _testhelper_flip_matching(base, _testhelper_binds(base, xs...))
        @test (xs, ER.Enumeration._degenerate(_testhelper_candidate(m)),
               _testhelper_degenerate(m)) ==
              (xs, true, true)
    end
end

@testset "_seed_variants: the variants follow each base's canonical group order" begin
    ER = EnzymeRates.Enumeration
    # A ter-ter seed whose merged base moves the D binding into EPQR. That
    # binding's group then sorts after the P group, which binds P at EAB first, so in
    # the canonical group order of the base the P group comes before the D group. The
    # minimal flip sets, by level and then by their groups' canonical places, give the
    # variants in this order: the merged base flipped at {B, P}, {R, D}, {P, D}, {P, A},
    # {D, Q}, then the Theorell–Chance base (TC = EAB + D → EQR + P, which
    # replaces EPQR) flipped at {B, R, TC}, {B, TC, Q}, {R, TC, A}, {TC, A, Q}.
    seed = @enzyme_mechanism begin
        substrates: A, B, D; products: P, Q, R
        steps: begin
            E + B ⇌ EB
            E + R ⇌ ER
            EAB + D ⇌ EABD
            (EAB + P ⇌ EABP, EB + P ⇌ EBP, EQR + P ⇌ EPQR)
            EABD <--> EPQR
            (EB + A ⇌ EAB, EBP + A ⇌ EABP)
            ER + Q ⇌ EQR
        end
    end
    expected = [
        (@enzyme_mechanism begin
            substrates: A, B, D; products: P, Q, R
            steps: begin
                E + B <--> EB
                E + R ⇌ ER
                EAB + D ⇌ EPQR
                (EAB + P <--> EABP, EB + P <--> EBP,
                 EQR + P <--> EPQR)
                (EB + A ⇌ EAB, EBP + A ⇌ EABP)
                ER + Q ⇌ EQR
            end
        end),
        (@enzyme_mechanism begin
            substrates: A, B, D; products: P, Q, R
            steps: begin
                E + B ⇌ EB
                E + R <--> ER
                EAB + D <--> EPQR
                (EAB + P ⇌ EABP, EB + P ⇌ EBP, EQR + P ⇌ EPQR)
                (EB + A ⇌ EAB, EBP + A ⇌ EABP)
                ER + Q ⇌ EQR
            end
        end),
        (@enzyme_mechanism begin
            substrates: A, B, D; products: P, Q, R
            steps: begin
                E + B ⇌ EB
                E + R ⇌ ER
                EAB + D <--> EPQR
                (EAB + P <--> EABP, EB + P <--> EBP,
                 EQR + P <--> EPQR)
                (EB + A ⇌ EAB, EBP + A ⇌ EABP)
                ER + Q ⇌ EQR
            end
        end),
        (@enzyme_mechanism begin
            substrates: A, B, D; products: P, Q, R
            steps: begin
                E + B ⇌ EB
                E + R ⇌ ER
                EAB + D ⇌ EPQR
                (EAB + P <--> EABP, EB + P <--> EBP,
                 EQR + P <--> EPQR)
                (EB + A <--> EAB, EBP + A <--> EABP)
                ER + Q ⇌ EQR
            end
        end),
        (@enzyme_mechanism begin
            substrates: A, B, D; products: P, Q, R
            steps: begin
                E + B ⇌ EB
                E + R ⇌ ER
                EAB + D <--> EPQR
                (EAB + P ⇌ EABP, EB + P ⇌ EBP, EQR + P ⇌ EPQR)
                (EB + A ⇌ EAB, EBP + A ⇌ EABP)
                ER + Q <--> EQR
            end
        end),
        (@enzyme_mechanism begin
            substrates: A, B, D; products: P, Q, R
            steps: begin
                E + B <--> EB
                E + R <--> ER
                (EAB + P ⇌ EABP, EB + P ⇌ EBP)
                (EB + A ⇌ EAB, EBP + A ⇌ EABP)
                ER + Q ⇌ EQR
                EAB + D <--> EQR + P
            end
        end),
        (@enzyme_mechanism begin
            substrates: A, B, D; products: P, Q, R
            steps: begin
                E + B <--> EB
                E + R ⇌ ER
                (EAB + P ⇌ EABP, EB + P ⇌ EBP)
                (EB + A ⇌ EAB, EBP + A ⇌ EABP)
                ER + Q <--> EQR
                EAB + D <--> EQR + P
            end
        end),
        (@enzyme_mechanism begin
            substrates: A, B, D; products: P, Q, R
            steps: begin
                E + B ⇌ EB
                E + R <--> ER
                (EAB + P ⇌ EABP, EB + P ⇌ EBP)
                (EB + A <--> EAB, EBP + A <--> EABP)
                ER + Q ⇌ EQR
                EAB + D <--> EQR + P
            end
        end),
        @enzyme_mechanism begin
            substrates: A, B, D; products: P, Q, R
            steps: begin
                E + B ⇌ EB
                E + R ⇌ ER
                (EAB + P ⇌ EABP, EB + P ⇌ EBP)
                (EB + A <--> EAB, EBP + A <--> EABP)
                ER + Q <--> EQR
                EAB + D <--> EQR + P
            end
        end
    ]
    @test _testhelper_move(ER._seed_variants, seed, ter_ter_rxn) == expected
end

@testset "_seed_candidate_screen agrees with the group predicates on every candidate" begin
    ER = EnzymeRates
    # The screen reads the turnover, maximal-rate and chemistry tests of a candidate from
    # index arrays. On every subset of the groups of eight bases flipped to steady state it
    # must give the verdict of `has_turnover_cycle`, `_has_vmax` both ways and
    # `_chemistry_equilibrates_both_sides` on the flipped flags.
    ordered = @enzyme_mechanism begin
        substrates: A, B; products: P, Q
        steps: begin
            E + A ⇌ EA
            EA + B ⇌ EAB
            EAB <--> EPQ
            EQ + P ⇌ EPQ
            E + Q ⇌ EQ
        end
    end
    shared_b = @enzyme_mechanism begin
        substrates: A, B; products: P, Q
        steps: begin
            E + A ⇌ EA
            (EA + B ⇌ EAB, EQ + B ⇌ EBQ)
            EAB <--> EPQ
            EQ + P ⇌ EPQ
            E + Q ⇌ EQ
        end
    end
    ordered_random = @enzyme_mechanism begin
        substrates: A, B; products: P, Q
        steps: begin
            E + A ⇌ EA
            EA + B ⇌ EAB
            EAB <--> EPQ
            (EQ + P ⇌ EPQ, E + P ⇌ EP)
            (EP + Q ⇌ EPQ, E + Q ⇌ EQ)
        end
    end
    ping_pong = @enzyme_mechanism begin
        substrates: A, B; products: P, Q
        steps: begin
            E + A ⇌ EA
            EA <--> EP_res_plusA_minusP
            E_res_plusA_minusP + P ⇌ EP_res_plusA_minusP
            E_res_plusA_minusP + B ⇌ EB_res_plusA_minusP
            EB_res_plusA_minusP ⇌ EQ
            E + Q ⇌ EQ
        end
    end
    # Bases: merged and Theorell–Chance for the two ordered seeds, merged only for the
    # ordered/random seed, merged and both halves for the ping-pong seed.
    all_bases = [b for seed in (ordered, shared_b, ordered_random, ping_pong)
                 for b in _testhelper_bases(seed)]
    @test [length(_testhelper_group_members(b.catalytic)) for b in all_bases] ==
          [4, 3, 4, 4, 4, 4, 3, 3]
    # The screen's verdict and the predicates' verdict on every subset of the groups of the
    # base `m`, every edge of which is at rapid equilibrium.
    function _testhelper_sweep(m)
        screen = ER.Enumeration._seed_candidate_screen(_testhelper_candidate(m))
        map(0:(2^length(_testhelper_group_members(m.catalytic)) - 1)) do bits
            ss = BitVector([isodd(bits >> (g - 1))
                            for g in 1:length(_testhelper_group_members(m.catalytic))])
            re = BitVector([!ss[e.kinetic_group] for e in m.catalytic.steps])
            expected = !ER.Mechanisms.has_turnover_cycle(m.catalytic, re) &&
                ER.Enumeration._has_vmax(m.catalytic, re, m.catalytic.substrates) &&
                ER.Enumeration._has_vmax(m.catalytic, re, m.catalytic.products) &&
                !ER.Enumeration._chemistry_equilibrates_both_sides(_testhelper_candidate(m),
                                                                   re)
            (ss, screen(re), expected)
        end
    end
    verdicts = Bool[]
    for base in all_bases, (ss, screened, expected) in _testhelper_sweep(base)
        @test (ss, screened) == (ss, expected)
        push!(verdicts, expected)
    end
    @test length(verdicts) == 16 + 8 + 16 + 16 + 16 + 16 + 8 + 8
    @test 0 < count(verdicts) < length(verdicts)

    # Aggregate pin over the whole bi_bi_pp_rxn seed set: the 86 bases of its 62 seeds
    # (merged, and Theorell–Chance where ELIM applies) and every subset of their groups.
    # 8 bases of decorated ping-pong seeds hold a fused binding of a product, which makes
    # no chemistry node: the merge turns a dead-end binding of P at the substrate-side
    # form into one, and leaves that merged complex with three steps. The 35 merged bases
    # of the seeds with a random-order side hold a merged complex of three or four steps.
    pp_seeds = filter(_testhelper_holds_iso,
                      _testhelper_mechanisms(ER.Enumeration.init_mechanisms(bi_bi_pp_rxn)))
    pp_bases = [b for seed in pp_seeds for b in _testhelper_bases(seed)]
    # The fused bindings of base `b` that bind one of the metabolites `ms`, not as a
    # modifier copy.
    _testhelper_fused(b, ms) =
        (c = _testhelper_candidate(b);
         [e for (k, e) in enumerate(b.catalytic.steps)
          if _testhelper_is_fused(c, e) && !_testhelper_binds_copy(c, k) &&
             only(e.taken_up) in ms])
    # Whether a fused binding of a substrate enters a form of three or more edges.
    _testhelper_three_step_complex(b) =
        any(_testhelper_fused(b, b.catalytic.substrates)) do e
        count(d -> e.to in (d.from, d.to), b.catalytic.steps) >= 3
    end
    @test (length(pp_seeds), length(pp_bases)) == (62, 86)
    @test count(b -> !isempty(_testhelper_fused(b, b.catalytic.products)), pp_bases) == 8
    @test count(b -> !isempty(_testhelper_fused(b, b.catalytic.products)) &&
                     _testhelper_three_step_complex(b), pp_bases) == 8
    @test count(_testhelper_three_step_complex, pp_bases) == 35 + 8
    outcomes = [o for b in pp_bases for o in _testhelper_sweep(b)]
    @test length(outcomes) == 1344
    @test isempty([(ss, screened) for (ss, screened, expected) in outcomes
                   if screened != expected])
    @test count(o -> o[3], outcomes) == 484
end

@testset "_seed_candidate_screen refuses a graph that holds an isomerization" begin
    # A base has every isomerization of its seed merged, and the screen takes each
    # chemistry node for a merged complex, so it refuses the seed itself.
    seed = @enzyme_mechanism begin
        substrates: S
        products: P
        steps: begin
            E + S ⇌ ES
            ES <--> EP
            EP ⇌ E + P
        end
    end
    @test_throws(ErrorException("_seed_candidate_screen: the seed base holds an " *
                                "isomerization; a base has every isomerization merged"),
                 EnzymeRates.Enumeration._seed_candidate_screen(
                     _testhelper_candidate(seed)))
end

@testset "_productive_twin" begin
    ER = EnzymeRates
    # The form of graph `m` whose name (`_form_name`) is `nm`.
    _testhelper_form_called(m, nm) =
        only(f for f in eachindex(m.catalytic.enzyme_forms)
             if String(m.catalytic.enzyme_forms[f]) == "$nm")
    # A twin is asked for a modifier copy by its name, and 0 means none.

    # Ordered bi-bi. A copy of A at E has the composition of EA, its twin; at EQ
    # it is new; a copy of B at E is new; a foreign inhibitor is never a twin.
    ordered = @enzyme_mechanism begin
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
    twin = _testhelper_twin(ordered)
    @test twin(_testhelper_form_called(ordered, :E), :A) ==
          _testhelper_form_called(ordered, :EA)
    @test twin(_testhelper_form_called(ordered, :EQ), :A) == 0
    @test twin(_testhelper_form_called(ordered, :E), :B) == 0
    @test twin(_testhelper_form_called(ordered, :E), :I) == 0 &&
          twin(_testhelper_form_called(ordered, :EQ), :I) == 0

    # With the abortive complex EAQ present, the copy at EQ is a twin too.
    abortive = @enzyme_mechanism begin
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
    @test _testhelper_twin(abortive)(_testhelper_form_called(abortive, :EQ), :A) ==
          _testhelper_form_called(abortive, :EAQ)

    # A copy whose complex already exists is judged against the productive forms only.
    with_copy = @enzyme_mechanism begin
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
    twin_c = _testhelper_twin(with_copy)
    @test twin_c(_testhelper_form_called(with_copy, :E), :A) ==
          _testhelper_form_called(with_copy, :EA)
    @test twin_c(_testhelper_form_called(with_copy, :EQ), :A) == 0

    # The offsets key. E ⇌ E* is a rapid-equilibrium isomerization, so E and E*
    # have equal offsets; a copy of S at E has a composition no form has, yet its
    # weight is proportional to E*(S)'s. With the isomerization at steady state the
    # two forms sit in different segments and the copy is not a twin.
    iso = @enzyme_mechanism begin
        substrates: S
        products: P
        steps: begin
            E ⇌ Estar
            Estar + S ⇌ EstarS
            EstarS <--> EP
            E + P ⇌ EP
        end
    end
    estar_s = _testhelper_form_called(iso, :EstarS)
    @test _testhelper_twin(iso)(_testhelper_form_called(iso, :E), :S) == estar_s
    # Twins are found by offsets; composition only breaks a tie between forms with the
    # complex's offsets. This twin differs from the complex in conformation or residual,
    # so no form has the complex's composition.
    records = _testhelper_candidate(iso).forms
    _testhelper_label_residual(f) =
        (records[f].label, records[f].adds, records[f].subtracts)
    @test _testhelper_label_residual(estar_s) !=
          _testhelper_label_residual(_testhelper_form_called(iso, :E))
    iso_ss = @enzyme_mechanism begin
        substrates: S
        products: P
        steps: begin
            E <--> Estar
            Estar + S ⇌ EstarS
            EstarS <--> EP
            E + P ⇌ EP
        end
    end
    @test _testhelper_twin(iso_ss)(_testhelper_form_called(iso_ss, :E), :S) == 0

    # Several forms can tie as twins. E ⇌ F and E ⇌ F2 put FS and F2S at the offsets of
    # a copy of S at E, and no form holds the copy's composition, S at E. The twin is then
    # the first of them in form order, F2S, by name.
    tied = @enzyme_mechanism begin
        substrates: S
        products: P
        steps: begin
            E ⇌ F
            E ⇌ F2
            F + S ⇌ FS
            F2 + S ⇌ F2S
            FS <--> EP
            E + P ⇌ EP
        end
    end
    @test _testhelper_twin(tied)(_testhelper_form_called(tied, :E), :S) ==
          _testhelper_form_called(tied, :F2S)

    # A form of the complex's composition outside the site's segment is no twin. With
    # A bound at steady state, EA has the composition of E·A* but sits in another
    # segment, its weight set by the steady state rather than by K·A·[E]: the two
    # weights are not proportional.
    ss_bound = @enzyme_mechanism begin
        substrates: A, B
        products: P, Q
        steps: begin
            E + A <--> EA
            EA + B ⇌ EAB
            EAB <--> EPQ
            EQ + P ⇌ EPQ
            E + Q ⇌ EQ
        end
    end
    @test _testhelper_twin(ss_bound)(_testhelper_form_called(ss_bound, :E), :A) == 0

    # Ping-pong with the second chemistry step at rapid equilibrium, an abortive
    # E(B, P; residual) and a copy of P at E. A copy of Q at EPmod has the new
    # composition {P, Q}, but E(B; residual) ⇌ EQ puts E(B, P; residual) at the
    # offsets P + Q of E's segment, the same as EPmod·Q: a twin. At
    # E(; residual) the copy is not.
    pp = @enzyme_mechanism begin
        substrates: A, B
        products: P, Q
        steps: begin
            E + A ⇌ EA
            EA <--> EP_res_plusA_minusP
            EP_res_plusA_minusP ⇌ E_res_plusA_minusP + P
            E_res_plusA_minusP + B ⇌ EB_res_plusA_minusP
            EB_res_plusA_minusP ⇌ EQ
            EQ ⇌ E + Q
            EB_res_plusA_minusP + P ⇌ EBP_res_plusA_minusP
            E + P ⇌ EPmod
        end
    end
    twin_pp = _testhelper_twin(pp)
    @test twin_pp(_testhelper_form_called(pp, :EPmod), :Q) ==
          _testhelper_form_called(pp, :EBP_res_plusA_minusP)
    @test twin_pp(_testhelper_form_called(pp, :E_res_plusA_minusP), :Q) == 0

    # Groups that bind a copy redundantly: only at twin sites, with a consistent gauge.
    @test isempty(ER.Enumeration._redundant_copy_groups(_testhelper_candidate(with_copy)))

    # Two copies. E(A, Q*) has the composition of E(A*, Q), but a form bound to a
    # competitive inhibitor is no twin source: only productive forms are. So the Q*
    # site at EA is not a twin, the Q* site at E duplicates EQ, and both copy
    # groups hold a site with no twin.
    two_copies = @enzyme_mechanism begin
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
    twin_two = _testhelper_twin(two_copies)
    @test twin_two(_testhelper_form_called(two_copies, :EA), :Q) == 0
    @test twin_two(_testhelper_form_called(two_copies, :E), :Q) ==
          _testhelper_form_called(two_copies, :EQ)
    @test twin_two(_testhelper_form_called(two_copies, :EQ), :A) == 0
    @test isempty(ER.Enumeration._redundant_copy_groups(_testhelper_candidate(two_copies)))

    # Conformational states. A copy binds the active state always and the inactive state
    # unless it is `:OnlyA`. It is redundant when every complex has a twin in every state
    # where it binds and the gauge holds over the states together. In each fixture below,
    # the copy E + S* ⇌ E(S*) and every other kinetic group are single steps, and a single
    # step has one rescaling, so the gauge holds in each state where E(S*) has a twin. The
    # states agree too. The only `:EqualAI` groups are E + P ⇌ EP, which touches neither
    # ES nor E(S*) and takes factor 1 at both ends in each state, and an `:EqualAI` copy,
    # which takes 1 at E and its one factor s at E(S*) in each. So a copy here is
    # redundant exactly when E(S*) has a twin in every state where it binds. In the active
    # state that twin is ES, formed from E by the RE S binding.
    # S binding `:OnlyA`: the inactive state holds E, EP and an `:EqualAI` copy's
    # E(S*), and no productive form bound to S. The copy is new there and kept.
    onlya_with_equalai_copy = @enzyme_mechanism begin
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
    @test isempty(ER.Enumeration._redundant_copy_groups(
        _testhelper_candidate(onlya_with_equalai_copy)))
    @test _testhelper_fitted(onlya_with_equalai_copy) ==
        _testhelper_identifiable_rank(onlya_with_equalai_copy)
    # An `:OnlyA` copy binds the active state only, where ES is its twin: redundant,
    # its constant one fitted parameter above the rank.
    onlya_with_onlya_copy = @enzyme_mechanism begin
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
    @test !isempty(
        ER.Enumeration._redundant_copy_groups(_testhelper_candidate(onlya_with_onlya_copy)))
    @test _testhelper_fitted(onlya_with_onlya_copy) ==
        _testhelper_identifiable_rank(onlya_with_onlya_copy) + 1
    # S binding `:NonequalAI`: ES exists in both states and is the twin in each, so an
    # `:EqualAI` copy is redundant.
    nonequal_with_equalai_copy = @enzyme_mechanism begin
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
    @test !isempty(ER.Enumeration._redundant_copy_groups(
        _testhelper_candidate(nonequal_with_equalai_copy)))
    @test _testhelper_fitted(nonequal_with_equalai_copy) ==
        _testhelper_identifiable_rank(nonequal_with_equalai_copy) + 1
    # Every catalytic group `:OnlyA`: the inactive state has no catalytic step, yet it
    # holds its free enzyme, so an `:EqualAI` copy binds E there, where E(S*) has no
    # twin: kept. An `:OnlyA` copy binds the active state only and is redundant.
    all_onlya_with_equalai_copy = @enzyme_mechanism begin
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
    @test isempty(ER.Enumeration._redundant_copy_groups(
        _testhelper_candidate(all_onlya_with_equalai_copy)))
    @test _testhelper_fitted(all_onlya_with_equalai_copy) ==
        _testhelper_identifiable_rank(all_onlya_with_equalai_copy)
    all_onlya_with_onlya_copy = @enzyme_mechanism begin
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
    @test !isempty(ER.Enumeration._redundant_copy_groups(
        _testhelper_candidate(all_onlya_with_onlya_copy)))
    @test _testhelper_fitted(all_onlya_with_onlya_copy) ==
        _testhelper_identifiable_rank(all_onlya_with_onlya_copy) + 1
end

@testset "_redundant_copy_groups: an :OnlyI copy binds the inactive state alone" begin
    ER = EnzymeRates
    # The copy binds E in the inactive conformation, where ES is its twin and the
    # :NonequalAI S binding's inactive constant absorbs it. Its binding at EJS lies
    # in neither law: the active law drops the :OnlyI group, and the inactive law prunes
    # EJS behind the :OnlyA J binding. Read over every step, EJSSmod would
    # have no twin and the copy would be kept.
    copy_beyond = @enzyme_mechanism begin
        substrates: S
        products: P
        modifiers: J
        oligomeric_state: 2
        steps: begin
            E + S ⇌ ES                :: NonequalAI
            ES <--> EP              :: OnlyA
            E + P ⇌ EP                :: EqualAI
            ES + J ⇌ EJS          :: OnlyA
            (E + S ⇌ ESmod, EJS + S ⇌ EJSSmod) :: OnlyI
        end
    end
    g = only(_testhelper_groups_binding_copy(copy_beyond, :S))
    @test ER.Enumeration._redundant_copy_groups(_testhelper_candidate(copy_beyond)) == [g]
    @test _testhelper_fitted(copy_beyond) == _testhelper_identifiable_rank(copy_beyond) + 1
    # Tagged :EqualAI, the copy also binds the active conformation, at E and at EJS,
    # where EJSSmod has no twin: kept.
    equalai_copy = @enzyme_mechanism begin
        substrates: S
        products: P
        modifiers: J
        oligomeric_state: 2
        steps: begin
            E + S ⇌ ES                :: NonequalAI
            ES <--> EP              :: OnlyA
            E + P ⇌ EP                :: EqualAI
            ES + J ⇌ EJS          :: OnlyA
            (E + S ⇌ ESmod, EJS + S ⇌ EJSSmod) :: EqualAI
        end
    end
    @test isempty(
        ER.Enumeration._redundant_copy_groups(_testhelper_candidate(equalai_copy)))
    @test _testhelper_fitted(equalai_copy) == _testhelper_identifiable_rank(equalai_copy)
    # Every catalytic group is :OnlyA, so the inactive law holds free E alone, and the
    # active law drops the :OnlyI I binding: the copy binds in neither state, and its
    # constant is not fitted.
    neither = @enzyme_mechanism begin
        substrates: S
        products: P
        modifiers: I
        oligomeric_state: 2
        steps: begin
            E + S ⇌ ES                :: OnlyA
            ES <--> EP              :: OnlyA
            E + P ⇌ EP                :: OnlyA
            ES + I ⇌ EIS          :: OnlyI
        end
    end
    @test isempty(ER.Enumeration._redundant_copy_groups(_testhelper_candidate(neither)))
    @test _testhelper_fitted(neither) == _testhelper_identifiable_rank(neither) == 4
end

@testset "_redundant_copy_groups: every complex a twin and a consistent dwell gauge" begin
    # A copy group is redundant when every complex duplicates a productive form and the
    # dwell gauge exists: scale factors σ, ρ_class on
    # each twin, one factor s on every complex of the copy and 1 elsewhere, under which
    # every kinetic group, the copy's own included, rescales its constants alike: an RE
    # group needs one ratio σ(from)/σ(to), an SS group one σ(from) and one σ(to). Twins
    # formed by RE bindings of the copied ligand in one kinetic group share a class
    # (ρ = K/(K + K*)); every other twin is its own. The copy's group, 1 at each site,
    # holds all its complexes at one factor, and a complex's weight moves opposite to
    # its twin's, so no complex takes its twin's factor.
    ER = EnzymeRates
    # The one group of `m` that binds a modifier copy.
    _testhelper_copy_group(m) =
        (c = _testhelper_candidate(m);
         only(unique(e.kinetic_group for (k, e) in enumerate(m.catalytic.steps)
                     if _testhelper_binds_copy(c, k))))

    # B1: ordered bi-bi with A* at E. The twin EA is formed
    # by E + A ⇌ EA, alone in its group, so the A group's ratio is 1/ρ, the B group's
    # ρ and the chemistry's forward factor ρ, each over one step: the gauge exists.
    # 6 fitted, rank 5.
    b1 = @enzyme_mechanism begin
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
    @test ER.Enumeration._redundant_copy_groups(
        _testhelper_candidate(b1)) == [_testhelper_copy_group(b1)]
    @test _testhelper_fitted(b1) == 6 && _testhelper_identifiable_rank(b1) == 5

    # Case 3: the A group also binds at EQ, forming the abortive EAQ. The copy at
    # E duplicates EA, but the A group then holds a step with ratio 1/ρ (at E) and
    # one with ratio 1 (at EQ, no copy): no gauge, and the copy's constant is
    # identifiable. 6 fitted, rank 6.
    case3 = @enzyme_mechanism begin
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
    g3 = _testhelper_copy_group(case3)
    @test ER.Enumeration._all_twin(_testhelper_candidate(case3),
                                   _testhelper_group_members(case3.catalytic)[g3],
                                   _testhelper_twin(case3))
    @test isempty(ER.Enumeration._redundant_copy_groups(_testhelper_candidate(case3)))
    @test _testhelper_fitted(case3) == 6 && _testhelper_identifiable_rank(case3) == 6

    # B14: the same parent with A* at E and EQ, the Q binding E + Q ⇌ EQ mirrored
    # onto the copy forms. Both twins, EA and EAQ, are formed by the one A
    # group, so they share ρ: the A group's two steps both have ratio 1/ρ, and the Q
    # group's three steps (E → EQ, EA → EAQ and the mirror between the two
    # complexes, both at s) all have ratio 1. The gauge exists. 6 fitted, rank 5.
    b14 = @enzyme_mechanism begin
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
    @test ER.Enumeration._redundant_copy_groups(
        _testhelper_candidate(b14)) == [_testhelper_copy_group(b14)]
    @test _testhelper_fitted(b14) == 6 && _testhelper_identifiable_rank(b14) == 5

    # H1: A binds E in one group and EB in another, so the twins EA and EAB
    # have different ρ. The B group shared by EA + B ⇌ EAB (ratio ρ₁/ρ₂) and
    # EQ + B ⇌ EBQ (ratio 1) cannot absorb both: no gauge. (The B group at E
    # and its mirror between the complexes have ratio 1.) 7 fitted, rank 7.
    h1 = @enzyme_mechanism begin
        substrates: A, B
        products: P, Q
        steps: begin
            E + A ⇌ EA
            EB + A ⇌ EAB
            (E + B ⇌ EB, EAmod + B ⇌ EAmodB)
            (EA + B ⇌ EAB, EQ + B ⇌ EBQ)
            EAB <--> EPQ
            EQ + P ⇌ EPQ
            E + Q ⇌ EQ
            (E + A ⇌ EAmod, EB + A ⇌ EAmodB)
        end
    end
    gh = _testhelper_copy_group(h1)
    @test ER.Enumeration._all_twin(_testhelper_candidate(h1),
                                   _testhelper_group_members(h1.catalytic)[gh],
                                   _testhelper_twin(h1))
    @test isempty(ER.Enumeration._redundant_copy_groups(_testhelper_candidate(h1)))
    @test _testhelper_fitted(h1) == 7 && _testhelper_identifiable_rank(h1) == 7

    # A step that leaves a complex takes the complex's factor s, never its twin's. B14
    # plus P* at the complex E(A*) and at the twin EAQ: the P* group holds
    # E(A*) → E(A*, P*) at s beside EAQ → E(A, P*, Q) at ρ, two ratios, so no gauge,
    # and the copy's constant shows (7 fitted, rank 7, as without the copy). The
    # mechanism is a valid parent.
    p_at_twin = @enzyme_mechanism begin
        substrates: A, B
        products: P, Q
        steps: begin
            (E + A ⇌ EA, EQ + A ⇌ EAQ)
            EA + B ⇌ EAB
            EAB <--> EPQ
            EQ + P ⇌ EPQ
            (E + Q ⇌ EQ, EA + Q ⇌ EAQ, EAmod + Q ⇌ EAmodQ)
            (E + A ⇌ EAmod, EQ + A ⇌ EAmodQ)
            (EAmod + P ⇌ EAmodPmod, EAQ + P ⇌ EAPmodQ)
        end
    end
    @test isempty(ER.Enumeration._redundant_copy_groups(_testhelper_candidate(p_at_twin)))
    @test ER.Enumeration._assert_emission_rules(
        _testhelper_candidate(p_at_twin)) === nothing
    @test _testhelper_fitted(p_at_twin) == 7 &&
          _testhelper_identifiable_rank(p_at_twin) == 7
    # P* at E(A*) and at its own twin EA instead: E(A*, P*) pairs with E(A, P*), and a
    # phantom remains (7 fitted, rank 6). The gauge does not reach it: it reads
    # E(A*) → E(A*, P*) at s beside EA → E(A, P*) at ρ, and the copy is kept.
    p_at_own_twin = @enzyme_mechanism begin
        substrates: A, B
        products: P, Q
        steps: begin
            (E + A ⇌ EA, EQ + A ⇌ EAQ)
            EA + B ⇌ EAB
            EAB <--> EPQ
            EQ + P ⇌ EPQ
            (E + Q ⇌ EQ, EA + Q ⇌ EAQ, EAmod + Q ⇌ EAmodQ)
            (E + A ⇌ EAmod, EQ + A ⇌ EAmodQ)
            (EAmod + P ⇌ EAmodPmod, EA + P ⇌ EAPmod)
        end
    end
    @test isempty(
        ER.Enumeration._redundant_copy_groups(_testhelper_candidate(p_at_own_twin)))
    @test _testhelper_fitted(p_at_own_twin) == 7 &&
          _testhelper_identifiable_rank(p_at_own_twin) == 6

    # An RE mirror joins two complexes, both at s, so its ratio is 1, as is its parent
    # step's between two sites, whatever the twins' classes. Random-order A in two
    # groups with A* at {E, EB}: twins EA and EAB of two classes, and
    # E + B ⇌ EB mirrored onto the complexes. Redundant: 7 fitted, rank 6.
    two_classes = @enzyme_mechanism begin
        substrates: A, B
        products: P, Q
        steps: begin
            E + A ⇌ EA
            EB + A ⇌ EAB
            (E + B ⇌ EB, EAmod + B ⇌ EAmodB)
            EA + B ⇌ EAB
            EAB <--> EPQ
            EQ + P ⇌ EPQ
            E + Q ⇌ EQ
            (E + A ⇌ EAmod, EB + A ⇌ EAmodB)
        end
    end
    @test ER.Enumeration._redundant_copy_groups(_testhelper_candidate(two_classes)) ==
          [_testhelper_copy_group(two_classes)]
    @test _testhelper_fitted(two_classes) == 7 &&
          _testhelper_identifiable_rank(two_classes) == 6
    # B14 with the A bindings and EA + Q ⇌ EAQ each in a group of their own: the
    # twins EA and EAQ have two classes, and the Q mirror joins their complexes.
    # Redundant: 7 fitted, rank 6. With the Q group at steady state the mirror carries
    # flux between the complexes, its forward constant needs s beside its parent step's
    # 1, and the gauge fails.
    split_b14 = @enzyme_mechanism begin
        substrates: A, B
        products: P, Q
        steps: begin
            E + A ⇌ EA
            EQ + A ⇌ EAQ
            EA + B ⇌ EAB
            EAB <--> EPQ
            EQ + P ⇌ EPQ
            (E + Q ⇌ EQ, EAmod + Q ⇌ EAmodQ)
            EA + Q ⇌ EAQ
            (E + A ⇌ EAmod, EQ + A ⇌ EAmodQ)
        end
    end
    ss_mirror = @_testhelper_unguarded @enzyme_mechanism begin
        substrates: A, B
        products: P, Q
        steps: begin
            E + A ⇌ EA
            EQ + A ⇌ EAQ
            EA + B ⇌ EAB
            EAB <--> EPQ
            EQ + P ⇌ EPQ
            (E + Q <--> EQ, EAmod + Q <--> EAmodQ)
            EA + Q ⇌ EAQ
            (E + A ⇌ EAmod, EQ + A ⇌ EAmodQ)
        end
    end
    @test ER.Enumeration._redundant_copy_groups(_testhelper_candidate(split_b14)) ==
          [_testhelper_copy_group(split_b14)]
    @test _testhelper_fitted(split_b14) == 7 &&
          _testhelper_identifiable_rank(split_b14) == 6
    @test isempty(ER.Enumeration._redundant_copy_groups(_testhelper_candidate(ss_mirror)))
    # The guard rejects ss_mirror at zero substrates, where its segment's steady state
    # depends on how fast its rapid-equilibrium steps run; its law is still counted
    # and ranked, unchecked.
    v = ER.Derivation._guard(ss_mirror)
    @test (v.fault, v.limit, v.face) ===
          (ER.Derivation._SPLIT, ER.Derivation._LIMIT_SPEED, 3)
    @test ER.Mechanisms.step_text(ss_mirror.catalytic, v.edge) == "E + A ⇌ EA"
    @test _testhelper_respeed_verdicts(ss_mirror; rng = Random.MersenneTwister(94),
                                       law = :limit)[3] === :pattern
    @test ER.Derivation.fitted_parameter_count(ss_mirror) == 8 &&
          _testhelper_identifiable_rank(ss_mirror) == 6
    # One factor for all of a copy's complexes, whatever their twins' classes: the
    # previous mechanism plus P* at the two complexes E(A*) and E(A*, Q), the Q binding
    # mirrored between the P* complexes. Both P* steps leave a complex at s, one ratio,
    # and the gauge exists: 8 fitted, rank 7.
    p_at_complexes = @enzyme_mechanism begin
        substrates: A, B
        products: P, Q
        steps: begin
            E + A ⇌ EA
            EQ + A ⇌ EAQ
            EA + B ⇌ EAB
            EAB <--> EPQ
            EQ + P ⇌ EPQ
            (E + Q ⇌ EQ, EAmod + Q ⇌ EAmodQ,
             EAmodPmod + Q ⇌ EAmodPmodQ)
            EA + Q ⇌ EAQ
            (E + A ⇌ EAmod, EQ + A ⇌ EAmodQ)
            (EAmod + P ⇌ EAmodPmod,
             EAmodQ + P ⇌ EAmodPmodQ)
        end
    end
    @test ER.Enumeration._redundant_copy_groups(_testhelper_candidate(p_at_complexes)) ==
          _testhelper_groups_binding_copy(p_at_complexes, :A)
    @test _testhelper_fitted(p_at_complexes) == 8 &&
          _testhelper_identifiable_rank(p_at_complexes) == 7

    # Case 1, uni-uni with S* at E: the twin ES is formed by the S group alone, and
    # the chemistry leaves it alone. 4 fitted, rank 3.
    uni = @enzyme_mechanism begin
        substrates: S
        products: P
        steps: begin
            E + S ⇌ ES
            ES <--> EP
            E + P ⇌ EP
            E + S ⇌ ESmod
        end
    end
    @test ER.Enumeration._redundant_copy_groups(
        _testhelper_candidate(uni)) == [_testhelper_copy_group(uni)]
    @test _testhelper_fitted(uni) == 4 && _testhelper_identifiable_rank(uni) == 3

    # A fused RE binding forms a twin as a plain one does. Merged uni-uni with A* at E:
    # E + A ⇌ EP is fused and rapid equilibrium, E + P → EP steady state, so EP
    # shares E's segment at the offsets {A: 1}, those of E·A*, and is its twin though
    # no form has the complex's composition. The fused group forms that twin alone, so
    # the twin's class is that group: its ratio is 1/ρ, P's binding takes ρ at EP,
    # and the copy s at its complex; the gauge exists. The law
    # E_t·k·(Keq·A − P)/(1 + A/K_A + A/K*), k the P binding's rate constant, sees K_A
    # and K* only through 1/K_A + 1/K*: 3 fitted, rank 2, against 2 and 2 without A*.
    fused_re = @enzyme_mechanism begin
        substrates: A
        products: P
        steps: begin
            E + A ⇌ EP
            E + P <--> EP
            E + A ⇌ EAmod
        end
    end
    fused_parent = @enzyme_mechanism begin
        substrates: A
        products: P
        steps: begin
            E + A ⇌ EP
            E + P <--> EP
        end
    end
    members = _testhelper_group_members(fused_re.catalytic)
    copy_g = _testhelper_copy_group(fused_re)
    fused_k = _testhelper_edge(fused_re, "E + A ⇌ EP")
    fused_g = fused_re.catalytic.steps[fused_k].kinetic_group
    p_g =
        fused_re.catalytic.steps[_testhelper_edge(fused_re, "E + P <--> EP")].kinetic_group
    twin = _testhelper_twin(fused_re)
    @test twin(fused_re.catalytic.steps[only(members[copy_g])].from, :A) ==
          fused_re.catalytic.steps[fused_k].to
    C(g) = ER.Enumeration._GaugeClass(g, 0, 1)
    @test ER.Enumeration._gauge_rescaling(_testhelper_candidate(fused_re), members, copy_g,
                                          twin, ER.Enumeration._one_state_class) ==
          Dict(copy_g => (ER.Enumeration._NO_RESCALING, C(copy_g)),
               fused_g => (ER.Enumeration._NO_RESCALING, C(fused_g)),
               p_g => (ER.Enumeration._NO_RESCALING, C(fused_g)))
    @test ER.Enumeration._redundant_copy_groups(_testhelper_candidate(fused_re)) == [copy_g]
    @test _testhelper_fitted(fused_re) == 3 && _testhelper_identifiable_rank(fused_re) == 2
    @test _testhelper_fitted(fused_parent) ==
          _testhelper_identifiable_rank(fused_parent) == 2

    # Ping-pong with the second chemistry step at rapid equilibrium, an abortive
    # E(B, P; residual) and P* at E, then Q* at E(P*). E(P*, Q*) has no form's
    # composition, but it shares E's segment and the offsets P + Q with
    # E(B, P; residual), an offsets twin, its own class. Only E(B; residual) + P ⇌
    # E(B, P; residual) touches the twin, alone in its group: the gauge exists, and
    # Q* adds one phantom to the parent's one (9 fitted, rank 7, against 8 and 7).
    pp = @enzyme_mechanism begin
        substrates: A, B
        products: P, Q
        steps: begin
            E + A ⇌ EA
            EA <--> EP_res_plusA_minusP
            EP_res_plusA_minusP ⇌ E_res_plusA_minusP + P
            E_res_plusA_minusP + B ⇌ EB_res_plusA_minusP
            EB_res_plusA_minusP ⇌ EQ
            EQ ⇌ E + Q
            EB_res_plusA_minusP + P ⇌ EBP_res_plusA_minusP
            E + P ⇌ EPmod
        end
    end
    pp_q = @enzyme_mechanism begin
        substrates: A, B
        products: P, Q
        steps: begin
            E + A ⇌ EA
            EA <--> EP_res_plusA_minusP
            EP_res_plusA_minusP ⇌ E_res_plusA_minusP + P
            E_res_plusA_minusP + B ⇌ EB_res_plusA_minusP
            EB_res_plusA_minusP ⇌ EQ
            EQ ⇌ E + Q
            EB_res_plusA_minusP + P ⇌ EBP_res_plusA_minusP
            E + P ⇌ EPmod
            EPmod + Q ⇌ EPmodQmod
        end
    end
    @test isempty(ER.Enumeration._redundant_copy_groups(_testhelper_candidate(pp)))
    @test ER.Enumeration._redundant_copy_groups(_testhelper_candidate(pp_q)) ==
          _testhelper_groups_binding_copy(pp_q, :Q)
    @test _testhelper_fitted(pp) == 8 && _testhelper_identifiable_rank(pp) == 7
    @test _testhelper_fitted(pp_q) == 9 && _testhelper_identifiable_rank(pp_q) == 7

    # A twin must share its site's segment. A binds E at steady state and E·A*
    # duplicates EA by composition only: EA lies in another segment, so merging
    # E·A* into it is no reparameterization, and the copy's constant shows (7 fitted,
    # rank 7, one above the mechanism without the copy).
    ss_twin = @enzyme_mechanism begin
        substrates: A, B
        products: P, Q
        steps: begin
            E + A <--> EA
            (E + P ⇌ EP, EQ + P ⇌ EPQ)
            (E + Q ⇌ EQ, EP + Q ⇌ EPQ)
            EA + B ⇌ EAB
            EAB <--> EPQ
            E + A ⇌ EAmod
        end
    end
    ss_parent = @enzyme_mechanism begin
        substrates: A, B
        products: P, Q
        steps: begin
            E + A <--> EA
            (E + P ⇌ EP, EQ + P ⇌ EPQ)
            (E + Q ⇌ EQ, EP + Q ⇌ EPQ)
            EA + B ⇌ EAB
            EAB <--> EPQ
        end
    end
    @test isempty(ER.Enumeration._redundant_copy_groups(_testhelper_candidate(ss_twin)))
    @test _testhelper_fitted(ss_twin) == 7 && _testhelper_identifiable_rank(ss_twin) == 7
    @test _testhelper_fitted(ss_parent) == 6 &&
          _testhelper_identifiable_rank(ss_parent) == 6

    # An `:EqualAI` group has one constant in both states, so it must take one rescaling.
    # A binds `:NonequalAI`, so the twin EA has ρ = K_A/(K_A + K*) with a different K_A
    # in each state, while the `:EqualAI` B binding leaves EA in both: the active state
    # asks K_B·ρ_A, the inactive K_B·ρ_I. Each state alone has a gauge; the two together
    # have none, and the copy adds an identifiable constant (rank one above the
    # mechanism without it).
    tied = @enzyme_mechanism begin
        substrates: A, B
        products: P, Q
        oligomeric_state: 2
        steps: begin
            E + A ⇌ EA                :: NonequalAI
            EA + B ⇌ EAB          :: EqualAI
            EAB <--> EPQ        :: NonequalAI
            EQ + P ⇌ EPQ          :: EqualAI
            E + Q ⇌ EQ                :: EqualAI
            E + A ⇌ EAmod      :: EqualAI
        end
    end
    tied_parent = @enzyme_mechanism begin
        substrates: A, B
        products: P, Q
        oligomeric_state: 2
        steps: begin
            E + A ⇌ EA                :: NonequalAI
            EA + B ⇌ EAB          :: EqualAI
            EAB <--> EPQ        :: NonequalAI
            EQ + P ⇌ EPQ          :: EqualAI
            E + Q ⇌ EQ                :: EqualAI
        end
    end
    @test isempty(ER.Enumeration._redundant_copy_groups(_testhelper_candidate(tied)))
    @test _testhelper_identifiable_rank(tied) ==
        _testhelper_identifiable_rank(tied_parent) + 1
end

@testset "_gauge_rescaling: group 2 and form 2 are two classes" begin
    # A binds E and Estar, which E ⇌ Estar keeps at one weight ratio, so both
    # complexes have the twin EA, form 2, in a class of its own; the copy's group is
    # group 2. P binds the complex EAmod and the twin EA in one group, which
    # then needs the complexes' factor at one step and the twin's at the other: the
    # A gauge fails. A class keyed by one index would make group 2 and form 2 one
    # number, that group would rescale alike, and the gauge would hold.
    ER = EnzymeRates
    m = @enzyme_mechanism begin
        substrates: A
        products: P
        steps: begin
            E + A ⇌ EA
            EA <--> EP
            E + P ⇌ EP
            E ⇌ Estar
            (E + A ⇌ EAmod, Estar + A ⇌ EstarAmod)
            (EAmod + P ⇌ EAmodPmod, EA + P ⇌ EAPmod)
        end
    end
    members = _testhelper_group_members(m.catalytic)
    copy_g = m.catalytic.steps[_testhelper_edge(m, "E + A ⇌ EAmod")].kinetic_group
    twin = _testhelper_twin(m)
    @test copy_g == 2
    @test String(m.catalytic.enzyme_forms[2]) == "EA"
    @test [twin(m.catalytic.steps[k].from, :A) for k in members[copy_g]] == [2, 2]
    @test ER.Enumeration._GaugeClass(2, 0, 1) != ER.Enumeration._GaugeClass(0, 2, 1)
    @test ER.Enumeration._gauge_rescaling(_testhelper_candidate(m), members, copy_g, twin,
                                          ER.Enumeration._one_state_class) === nothing
    @test ER.Enumeration._gauge_rescaling(
        _testhelper_candidate(m), members, copy_g, twin,
        (h, t) -> ER.Enumeration._GaugeClass(h + t, 0, 1)) !== nothing
    @test isempty(ER.Enumeration._redundant_copy_groups(_testhelper_candidate(m)))
end

@testset "_hyperbolic_catalysis" begin
    # Ordered SS bi-bi: every metabolite binds on one edge and every segment is
    # a single form, so no denominator term carries a concentration twice.
    ordered = @enzyme_mechanism begin
        substrates: A, B
        products: P, Q
        steps: begin
            E + A <--> EA
            EA + B <--> EAB
            E + Q <--> EQ
            EQ + P <--> EPQ
            EAB <--> EPQ
        end
    end
    @test EnzymeRates.Enumeration._hyperbolic_catalysis(_testhelper_candidate(ordered))

    # Random SS bi-bi: A binds on E → EA and on EB → EAB; both edges lie
    # in one tree toward EAB, so the denominator carries A².
    random_ss = @enzyme_mechanism begin
        substrates: A, B
        products: P, Q
        steps: begin
            (E + A <--> EA, EB + A <--> EAB)
            (E + B <--> EB, EA + B <--> EAB)
            (E + P <--> EP, EQ + P <--> EPQ)
            (E + Q <--> EQ, EP + Q <--> EPQ)
            EAB <--> EPQ
        end
    end
    @test !EnzymeRates.Enumeration._hyperbolic_catalysis(_testhelper_candidate(random_ss))

    # Random RE bi-bi with the chemistry step as the only SS step: one segment,
    # so the equation is the rapid-equilibrium law, degree 1 throughout.
    random_re = @enzyme_mechanism begin
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
    @test EnzymeRates.Enumeration._hyperbolic_catalysis(_testhelper_candidate(random_re))

    # Random RE bi-bi with the A group flipped: the segment {EA, EAB} has
    # EAB carrying B beyond its bottom EA, and the edge EB → EAB into
    # it leaves a form carrying B, so the term rooted there carries B².
    random_flip_a = @enzyme_mechanism begin
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
    @test !EnzymeRates.Enumeration._hyperbolic_catalysis(
        _testhelper_candidate(random_flip_a))

    # Ordered SS bi-bi with substrate A binding EQ as an abortive complex at
    # rapid equilibrium: EAQ sits one A above EQ in its segment, and the
    # term rooted there holds E → EA, so the equation carries A². A binds its
    # catalytic site, so the power counts even though EAQ is a dead end.
    dead_end = @enzyme_mechanism begin
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
    @test !EnzymeRates.Enumeration._hyperbolic_catalysis(_testhelper_candidate(dead_end))

    # Uni-bi with random SS product release: P binds on E → EP and on
    # EQ → EPQ, both toward EPQ, so the denominator carries P².
    unibi_ss = @enzyme_mechanism begin
        substrates: S
        products: P, Q
        steps: begin
            E + S ⇌ ES
            ES <--> EPQ
            (E + P <--> EP, EQ + P <--> EPQ)
            (E + Q <--> EQ, EP + Q <--> EPQ)
        end
    end
    @test !EnzymeRates.Enumeration._hyperbolic_catalysis(_testhelper_candidate(unibi_ss))

    # Ping-pong with one SS chemistry step per half-reaction: two segments,
    # each metabolite carried once per term.
    pingpong = @enzyme_mechanism begin
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
    @test EnzymeRates.Enumeration._hyperbolic_catalysis(_testhelper_candidate(pingpong))

    # An allosteric mechanism is scored on its catalytic steps: uni-bi with the
    # P group at steady state carries Q² (the segment {EP, EPQ} plus the
    # edge EQ → EPQ).
    allo_unibi_flip_p = @enzyme_mechanism begin
        substrates: S
        products: P, Q
        oligomeric_state: 2
        steps: begin
            E + S ⇌ ES                                        :: NonequalAI
            ES <--> EPQ                                   :: NonequalAI
            (E + P <--> EP, EQ + P <--> EPQ)            :: NonequalAI
            (E + Q ⇌ EQ, EP + Q ⇌ EPQ)                  :: NonequalAI
        end
    end
    @test !EnzymeRates.Enumeration._hyperbolic_catalysis(
        _testhelper_candidate(allo_unibi_flip_p))

    # Ordered SS bi-bi whose A group also binds A as an abortive complex on
    # EQ, as the enumerator groups it. The abortive step went to steady state
    # with its group; the tree toward EAQ holds E → EA and EQ → EAQ,
    # both binding A, so the equation carries A².
    grouped_dead_end = @enzyme_mechanism begin
        substrates: A, B
        products: P, Q
        steps: begin
            (E + A <--> EA, EQ + A <--> EAQ)
            EA + B <--> EAB
            E + Q <--> EQ
            EQ + P <--> EPQ
            EAB <--> EPQ
        end
    end
    @test !EnzymeRates.Enumeration._hyperbolic_catalysis(
        _testhelper_candidate(grouped_dead_end))

    # Ping-pong whose second chemistry step is at rapid equilibrium, with B also
    # bound as an abortive complex on EQ in the same kinetic group as its
    # catalytic binding. One segment spans both halves: E sits one B above
    # E(; residual) and EBQ two above, so the segment weight carries B².
    pingpong_abortive = @enzyme_mechanism begin
        substrates: A, B
        products: P, Q
        steps: begin
            E + A ⇌ EA
            EA <--> EP_res_plusA_minusP
            EP_res_plusA_minusP ⇌ E_res_plusA_minusP + P
            (E_res_plusA_minusP + B ⇌ EB_res_plusA_minusP, EQ + B ⇌ EBQ)
            EB_res_plusA_minusP ⇌ EQ
            EQ ⇌ E + Q
        end
    end
    @test !EnzymeRates.Enumeration._hyperbolic_catalysis(
        _testhelper_candidate(pingpong_abortive))

    # Ordered SS bi-bi with the abortive complex EAQ reachable from both
    # EA + Q and EQ + A, each step grouped with its metabolite's catalytic
    # binding as the enumerator writes it. The tree toward EAQ holds
    # E → EA and EQ → EAQ, both binding A, so the equation carries A².
    two_sided_abortive = @enzyme_mechanism begin
        substrates: A, B
        products: P, Q
        steps: begin
            (E + A <--> EA, EQ + A <--> EAQ)
            EA + B <--> EAB
            (E + Q <--> EQ, EA + Q <--> EAQ)
            EQ + P <--> EPQ
            EAB <--> EPQ
        end
    end
    @test !EnzymeRates.Enumeration._hyperbolic_catalysis(
        _testhelper_candidate(two_sided_abortive))

    # A competitive-inhibitor square (I on E and on EA, with A binding the
    # inhibitor-bound form too) leaves the verdict unchanged: an inhibitor binds
    # a site of its own, so every step touching a form that carries it is left
    # out.
    rxn_with_i = @enzyme_reaction begin
        substrates: A[C], B[N]
        products: P[C], Q[N]
        competitive_inhibitors: I
    end
    ordered_ss = @enzyme_mechanism begin
        substrates: A, B
        products: P, Q
        steps: begin
            E + A <--> EA
            EA + B <--> EAB
            E + Q <--> EQ
            EQ + P <--> EPQ
            EAB <--> EPQ
        end
    end
    with_i = _testhelper_move(EnzymeRates.Enumeration._expand_add_dead_end_regulator,
                              ordered_ss, rxn_with_i)
    @test !isempty(with_i)
    @test all(EnzymeRates.Enumeration._hyperbolic_catalysis, _testhelper_candidate.(with_i))

    # A substrate declared as a dead-end inhibitor binds its inhibitor site under
    # the substrate's own name. Its bindings put A² and A³ into the derived
    # denominator, but those powers come from the inhibitor site, which the
    # predicate leaves out: the four placements on the ordered scheme pass, with
    # degrees 1, 2, 2 and 3. A binds at steady state, so EA sits in a segment of
    # its own and E·A* is no twin of it: the copy at E alone is emitted, and it adds
    # a phantom (10 fitted, rank 7, over the parent's 9 and 7).
    rxn_a_inhibits = @enzyme_reaction begin
        substrates: A[C], B[N]
        products: P[C], Q[N]
        competitive_inhibitors: A
    end
    # The forms (`_form_name`) where a step binds the copy A.
    _testhelper_copy_sites(k) =
        (c = _testhelper_candidate(k);
         Set(k.catalytic.enzyme_forms[e.from] for (j, e) in enumerate(k.catalytic.steps)
             if _testhelper_binds_copy(c, j)))
    with_a =
        _testhelper_move(EnzymeRates.Enumeration._expand_add_dead_end_regulator, ordered_ss,
                         rxn_a_inhibits)
    _testhelper_den_a_degree(k) = (L = EnzymeRates.Derivation._memo_law(k);
        maximum(numerator(EnzymeRates.Derivation._mexp(mono, L.vars.index[:A]))
                for mono in keys(L.D.terms)))
    @test length(with_a) == 4
    @test Set(_testhelper_copy_sites.(with_a)) ==
        Set([Set([:E]), Set([:E, :EQ]), Set([:EA, :EQ]), Set([:E, :EA])])
    @test all(EnzymeRates.Enumeration._hyperbolic_catalysis, _testhelper_candidate.(with_a))
    @test sort(_testhelper_den_a_degree.(with_a)) == [1, 2, 2, 3]
    at_e(ks) = only(k for k in ks if _testhelper_copy_sites(k) == Set([:E]))
    @test _testhelper_identifiable_rank(ordered_ss) == 7 &&
          _testhelper_fitted(ordered_ss) == 9
    @test _testhelper_identifiable_rank(at_e(with_a)) == 7 &&
          _testhelper_fitted(at_e(with_a)) == 10

    # Inhibitor-bound forms leave, but catalytic-site powers stay: adding the
    # inhibitor role of A to the abortive-complex scheme keeps it non-hyperbolic.
    # EA is formed at steady state, so E·A* has no twin and no placement is
    # redundant: all four are emitted. At {E, EQ}, where EQ·A* duplicates the
    # abortive EAQ, the copy adds an identifiable constant (11 fitted, rank 9, over
    # the parent's 10 and 8); at E alone it adds a phantom (11 fitted, rank 8).
    dead_end_with_a =
        _testhelper_move(EnzymeRates.Enumeration._expand_add_dead_end_regulator, dead_end,
                         rxn_a_inhibits)
    @test length(dead_end_with_a) == 4
    @test Set(_testhelper_copy_sites.(dead_end_with_a)) ==
        Set([Set([:E]), Set([:E, :EQ]), Set([:EA, :EQ]), Set([:E, :EA])])
    @test !any(EnzymeRates.Enumeration._hyperbolic_catalysis,
               _testhelper_candidate.(dead_end_with_a))
    at_e_eq = only(k for k in dead_end_with_a
                   if _testhelper_copy_sites(k) == Set([:E, :EQ]))
    @test _testhelper_identifiable_rank(dead_end) == 8 && _testhelper_fitted(dead_end) == 10
    @test _testhelper_identifiable_rank(at_e_eq) == 9 && _testhelper_fitted(at_e_eq) == 11
    @test _testhelper_identifiable_rank(at_e(dead_end_with_a)) == 8

    # An allosteric mechanism passes when each conformation's law does. The merged
    # abortive ping-pong at two subunits, Q's binding and the second half-reaction
    # :OnlyA: its catalytic graph is hyperbolic, but the inactive conformation, without
    # the second half, has two rapid-equilibrium segments joined only by the steady-state
    # first half, and each holds an abortive complex, EAP or E(A, P; residual = A -
    # P). A term rooted at one segment carries that complex's A or P, and the first half
    # leaves the other segment from a form that carries the same metabolite, so the
    # inactive law reaches A² and P². The predicate reads the graph alone: the guard
    # rejects this mechanism, since in the inactive conformation the :OnlyA step
    # EB_res_plusA_minusP ⇌ EQ bridges forms that the other steps keep apart.
    merged_abortive = @_testhelper_unguarded @enzyme_mechanism begin
        substrates: A, B ; products: P, Q ; oligomeric_state: 2
        steps: begin
            (E + A ⇌ EA,
             EP_res_plusA_minusP + A ⇌ EAP_res_plusA_minusP)   :: EqualAI
            E + Q ⇌ EQ                                               :: OnlyA
            E_res_plusA_minusP + B ⇌ EB_res_plusA_minusP         :: EqualAI
            (E_res_plusA_minusP + P ⇌ EP_res_plusA_minusP,
             EA + P ⇌ EAP)                                       :: EqualAI
            EA <--> EP_res_plusA_minusP                           :: EqualAI
            EB_res_plusA_minusP ⇌ EQ                              :: OnlyA
        end
    end
    @test EnzymeRates.Derivation._guard(merged_abortive).fault ===
          EnzymeRates.Derivation._BRIDGE
    _testhelper_inactive_degree(am, x) = (C = EnzymeRates.Derivation._memo_law(am).inactive;
        maximum(numerator(EnzymeRates.Derivation._mexp(mono, C.vars.index[x]))
                for mono in keys(C.D.terms)))
    @test EnzymeRates.Enumeration._hyperbolic_catalysis(
        _testhelper_candidate(merged_abortive.catalytic))
    @test !EnzymeRates.Enumeration._hyperbolic_catalysis(
        _testhelper_candidate(merged_abortive))
    @test _testhelper_inactive_degree(merged_abortive, :A) == 2
    @test _testhelper_inactive_degree(merged_abortive, :P) == 2
    # Split apart, the two abortive bindings each take a group of their own, and tagged
    # :OnlyA beside the second half-reaction they leave the inactive conformation no
    # abortive complex: both laws are hyperbolic. The guard rejects this one too: the
    # inactive conformation keeps the two sides of the rapid-equilibrium step
    # EB_res_plusA_minusP ⇌ EQ joined.
    split_abortive = @_testhelper_unguarded @enzyme_mechanism begin
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
    @test EnzymeRates.Derivation._guard(split_abortive).fault ===
          EnzymeRates.Derivation._FAST_STEP_PINNED
    @test EnzymeRates.Enumeration._hyperbolic_catalysis(
        _testhelper_candidate(split_abortive))
    @test all(x -> _testhelper_inactive_degree(split_abortive, x) <= 1, (:A, :B, :P, :Q))

    @test !EnzymeRates.Enumeration._requires_hyperbolic_catalysis(ordered)
    @test EnzymeRates.Enumeration._requires_hyperbolic_catalysis(allo_unibi_flip_p)
    # One catalytic subunit: the conformational equilibrium reweights each form
    # without raising any catalytic-site binding to a power, so the rule does
    # not apply.
    allo_unibi_flip_p_one_subunit =
        @enzyme_mechanism begin
            substrates: S
            products: P, Q
            oligomeric_state: 1
            steps: begin
                E + S ⇌ ES                                        :: NonequalAI
                ES <--> EPQ                                   :: NonequalAI
                (E + P <--> EP, EQ + P <--> EPQ)            :: NonequalAI
                (E + Q ⇌ EQ, EP + Q ⇌ EPQ)                  :: NonequalAI
            end
        end
    @test !EnzymeRates.Enumeration._requires_hyperbolic_catalysis(
        allo_unibi_flip_p_one_subunit)
end

@testset "_all_reach" begin
    # Segment graph 1 ⇄ 2 ⇄ 3 with edge tuples (source, target, form, mets).
    edges = [(1, 2, 0, Symbol[]), (2, 1, 0, Symbol[]),
             (2, 3, 0, Symbol[]), (3, 2, 0, Symbol[])]
    # Unpinned, every segment reaches every root.
    @test EnzymeRates.Enumeration._all_reach(3, edges, 3, ())
    # Pinning 1 → 2 keeps the path 1 → 2 → 3.
    @test EnzymeRates.Enumeration._all_reach(3, edges, 3, (edges[1],))
    # Pinning 2 → 1 strands 2 (and 1) away from root 3.
    @test !EnzymeRates.Enumeration._all_reach(3, edges, 3, (edges[2],))
    # A pinned 2-cycle (1 → 2 and 2 → 1) reaches no root outside it.
    @test !EnzymeRates.Enumeration._all_reach(3, edges, 3, (edges[1], edges[2]))
    # A root that is a pinned edge's target is fine.
    @test EnzymeRates.Enumeration._all_reach(3, edges, 2, (edges[1], edges[4]))
end

@testset "_hyperbolic_catalysis matches the derived denominator" begin
    # The structural predicate against the exponents of the derived denominator,
    # over every mechanism reachable from the seeds in two expansion levels. The
    # ping-pong-capable reaction's seeds include every bi-bi seed (pinned below),
    # and the moves, the predicate and the derivation read metabolite names, never atoms
    # (atoms enter enumeration only in seed construction and in the atom-conservation
    # assertion, which throws rather than filters), so its population holds every
    # bi-bi mechanism and the allosteric children of the bi-bi allosteric mechanisms.
    # The reactions declare no inhibitors, so the predicate sees every step and is
    # compared with each mechanism's own derived denominator. For an allosteric
    # mechanism the predicate is compared with both states, and the I-state is checked
    # to be hyperbolic whenever the A-state is. Allosteric mechanisms with the same
    # state graph share one derivation per state: their tags only rename parameters,
    # which changes no metabolite exponent.
    _testhelper_hyperbolic(C, mets) = all(e <= 1 for mono in keys(C.D.terms)
                                          for (v, e) in mono if C.vars.names[v] in mets)
    unibi = @enzyme_reaction begin
        substrates: S[AB]
        products: P[A], Q[B]
        oligomeric_state: 2
    end
    pingpong = @enzyme_reaction begin
        substrates: A[CX], B[N]
        products: P[C], Q[NX]
        oligomeric_state: 2
    end
    # An aggregate pin over two whole seed sets (the exception to writing every
    # fixture inline): every bi-bi seed is a seed of the ping-pong-capable reaction.
    bibi = @enzyme_reaction begin
        substrates: A[C], B[N]
        products: P[C], Q[N]
        oligomeric_state: 2
    end
    @test issubset(
        Set(_testhelper_mechanisms(EnzymeRates.Enumeration.init_mechanisms(bibi))),
        Set(_testhelper_mechanisms(EnzymeRates.Enumeration.init_mechanisms(pingpong))))
    n_checked = 0
    n_nonhyperbolic = 0
    for rxn in (unibi, pingpong)
        mets = Symbol[x.name for x in vcat(rxn.substrates, rxn.products)]
        level = EnzymeRates.Mechanism[
            m for m in _testhelper_mechanisms(EnzymeRates.Enumeration.init_mechanisms(rxn))]
        mechs = copy(level)
        for _ in 1:2
            level = unique!(_testhelper_children(rxn, level))
            append!(mechs, level)
        end
        unique!(mechs)
        derived = Dict{Tuple{Symbol, EnzymeRates.Mechanisms.CatalyticScheme, BitVector},
                       Bool}()
        function _testhelper_state_hyperbolic(m, state)
            runs = EnzymeRates.Mechanisms.steps_running_in_allosteric_conformations(m)
            get!(derived, (state, m.catalytic,
                           state === :A ? runs.active : runs.inactive)) do
                L = EnzymeRates.Derivation._memo_law(m)
                _testhelper_hyperbolic(state === :A ? L.active : L.inactive, mets)
            end
        end
        for m in mechs
            EnzymeRates.Derivation.rate_equation_term_count(m) <= 337 || continue
            structural = EnzymeRates.Enumeration._hyperbolic_catalysis(
                _testhelper_candidate(m))
            if !_testhelper_is_allosteric(m)
                @test structural ==
                      _testhelper_hyperbolic(EnzymeRates.Derivation._memo_law(m), mets)
            else
                hyp_a = _testhelper_state_hyperbolic(m, :A)
                @test structural == (hyp_a && _testhelper_state_hyperbolic(m, :I))
                @test _testhelper_state_hyperbolic(m, :I) || !hyp_a
            end
            n_checked += 1
            structural || (n_nonhyperbolic += 1)
        end
    end
    # Both classes must be exercised for the comparison to mean anything.
    @test n_checked > 100
    @test n_nonhyperbolic > 0
    @test n_nonhyperbolic < n_checked
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
        # carries Q²; flipping Q mirrors it with P². A plain Mechanism emits all
        # three; an allosteric parent emits only the S flip.
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
        allo_kids = _testhelper_move(EnzymeRates.Enumeration._expand_re_to_ss, allo)
        @test length(allo_kids) == 1
        @test Set(allo_kids) == Set([allo_flipS])
    end

    @testset "_expand_re_to_ss: a one-subunit allosteric parent keeps every flip" begin
        # The same uni-bi as an allosteric mechanism with one catalytic subunit.
        # The conformational equilibrium adds no concentration power, so the
        # flip move keeps the P and Q flips it drops for two subunits, and emits
        # the same three flips as the plain mechanism.
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
        kids_one = _testhelper_move(EnzymeRates.Enumeration._expand_re_to_ss, allo_one)
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
    # to-allosteric) under 120 s together; they return 81 children: 12 splits, 6 flips and
    # 63 K-types (the reaction declares no regulator). The split takes 4 s in a full local
    # suite and 7-11 s on a CI runner, but a slow macOS runner reached 60 s.
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
    @test length(rest.value) == 69
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
    @test length(kids) == 2797
    @test count(k -> _testhelper_is_allosteric(k), kids) == 1799
    @test all(c -> _testhelper_build(_testhelper_written(c)) == c, [init; kids])
end

@testset "candidate_mechanisms on a merged uni-uni" begin
    # E + S → EP fused (SS), E + P ⇌ EP (RE). Flip: the P group cuts E from EP,
    # carries flux, no bottomless segment: one child, both steps SS. Split: single-step
    # groups, nothing. Dead end: no regulator. To-allosteric at multiplicity 1: the
    # chemistry group (the fused step) is :OnlyA, the binding subsets range over {P}:
    # one K-type child; no regulator, so no V-type.
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
    k_type = @enzyme_mechanism begin
        substrates: S
        products: P
        steps: begin
            E + S <--> EP    :: OnlyA
            E + P ⇌ EP       :: OnlyA
        end
    end
    kids = _testhelper_children(rxn, [m])
    @test length(kids) == 2
    @test Set(kids) == Set([flipped, k_type])
    # The parent's rate is k·(S − P/Keq)/(1 + P/Kp), the Haldane relation fixing the
    # fused step's reverse constant: 2 fitted, rank 2. The flipped child's two forms at
    # steady state give V·(S − P/Keq)/(1 + a·S + b·P) with V, a and b free: 3 fitted,
    # rank 3. In the K-type child the inactive conformation binds nothing and runs no
    # chemistry, so L enters only as 1 + L beside the free enzyme and folds into k and
    # Kp: 3 fitted, rank 2.
    @test _testhelper_fitted(m) == _testhelper_identifiable_rank(m) == 2
    @test _testhelper_fitted(flipped) == _testhelper_identifiable_rank(flipped) == 3
    @test _testhelper_fitted(k_type) == 3 && _testhelper_identifiable_rank(k_type) == 2
    # The K-type child, whose inactive conformation is free E alone, is an L phantom
    # (spec section 13): its rescaling certificate maps it onto the law of its catalytic
    # mechanism.
    @test EnzymeRates.Derivation._is_phantom(EnzymeRates.Derivation._memo_law(k_type))
    @test _testhelper_phantom_map_holds(k_type; rng = Random.MersenneTwister(1))
end

@testset "_expand_change_allo_state on a merged uni-uni K-type" begin
    # The fused step E + S → EP is the chemistry group. Relaxing P to :NonequalAI
    # keeps the chemistry :OnlyA: the inactive conformation binds P with a constant of
    # its own and runs no chemistry, a valid child. The chemistry is not relaxed while
    # P's binding is :OnlyA: restored, it would join E to EP in the inactive
    # conformation beside the :OnlyA P binding, which then has no path to its limit
    # (admission's existence test). One child.
    k_type = @enzyme_mechanism begin
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
    kids = _testhelper_move(EnzymeRates.Enumeration._expand_change_allo_state, k_type)
    @test length(kids) == 1
    @test Set(kids) == Set([p_relaxed])
    # The child's rate is k·(S − P/Keq)/(1 + L + P·(1/Kp_A + L/Kp_I)): dividing by
    # 1 + L leaves k/(1 + L) and one P coefficient, so k, Kp_A, Kp_I and L give
    # 4 fitted, rank 2.
    @test _testhelper_fitted(p_relaxed) == 4 &&
          _testhelper_identifiable_rank(p_relaxed) == 2
end

@testset "_expand_to_allosteric on an ordered Theorell–Chance mechanism" begin
    # EA + B → EQ + P takes up B and gives off P in one step: the chemistry group.
    # Each K-type child sets it :OnlyA and ranges the binding subsets over {A, Q}:
    # {A}, {Q} and {A, Q}. In the inactive conformation the :OnlyA chemistry leaves
    # the A and Q bindings as the tree EA – E – EQ, so every subset's bindings
    # point away from E and admission refuses none. No regulator, so no V-type:
    # three children.
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
    only_a = @enzyme_mechanism begin
        substrates: A, B ; products: P, Q ; oligomeric_state: 1
        steps: begin
            E + A ⇌ EA                :: OnlyA
            EA + B <--> EQ + P      :: OnlyA
            E + Q ⇌ EQ                :: EqualAI
        end
    end
    only_q = @enzyme_mechanism begin
        substrates: A, B ; products: P, Q ; oligomeric_state: 1
        steps: begin
            E + A ⇌ EA                :: EqualAI
            EA + B <--> EQ + P      :: OnlyA
            E + Q ⇌ EQ                :: OnlyA
        end
    end
    a_and_q = @enzyme_mechanism begin
        substrates: A, B ; products: P, Q ; oligomeric_state: 1
        steps: begin
            E + A ⇌ EA                :: OnlyA
            EA + B <--> EQ + P      :: OnlyA
            E + Q ⇌ EQ                :: OnlyA
        end
    end
    kids = _testhelper_move(EnzymeRates.Enumeration._expand_to_allosteric, m, rxn)
    @test length(kids) == 3
    @test Set(kids) == Set([only_a, only_q, a_and_q])
    # The parent's rate is (k/Ka)·(A·B − P·Q/Keq)/(1 + A/Ka + Q/Kq): 3 fitted, rank 3.
    # Each child's dead inactive conformation adds (1 + L) times the bindings it keeps
    # to the bindings it loses, so dividing by 1 + L rescales the lost bindings'
    # constants and the turnover and L disappears: 4 fitted, rank 3. The rank is
    # checked on only_a and a_and_q; only_q is only_a with the reaction reversed.
    @test _testhelper_fitted(m) == _testhelper_identifiable_rank(m) == 3
    for k in kids
        @test _testhelper_fitted(k) == 4
    end
    for k in (only_a, a_and_q)
        @test _testhelper_identifiable_rank(k) == 3
    end
    # All three are L phantoms, each mapped onto the parent's law by its certificate,
    # only_a and only_q though their inactive conformation keeps a binding at free E.
    for k in (only_a, only_q, a_and_q)
        @test EnzymeRates.Derivation._is_phantom(EnzymeRates.Derivation._memo_law(k))
        @test _testhelper_phantom_map_holds(k; rng = Random.MersenneTwister(2))
    end
end

@testset "_expand_to_allosteric on a merged ordered bi-bi at multiplicity 2" begin
    # The ordered bi-bi seed with EAB merged into EPQ: B's binding is fused
    # with the chemistry (EA + B ⇌ EPQ), and the A and P bindings are steady
    # state. The fused B group is the chemistry group, :OnlyA in every K-type child,
    # and the subsets range over the A, P and Q groups: 2^3 − 1 = 7. Without the
    # :OnlyA chemistry group the bindings form the tree EA – E – EQ – EPQ. The
    # guard refuses one child: with A's binding :OnlyA too, EA is cut off in the
    # inactive conformation, and at zero A and B the fused binding, vanishing one way into
    # EA, fills it faster than A's steady-state binding empties it (`_FILLED`). The
    # scheme is hyperbolic, so
    # multiplicity 2 stays open: the RE segments are {E, EQ}, where EQ carries one
    # Q, and {EA, EPQ}, where EPQ carries one B; the steady-state A and P
    # bindings join them, and each metabolite scores on one edge only, an edge that
    # leaves the segment carrying that metabolite (`_hyperbolic_catalysis`). No
    # regulator, so no V-type: seven children.
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
    only_a = @_testhelper_unguarded @enzyme_mechanism begin
        substrates: A, B ; products: P, Q ; oligomeric_state: 2
        steps: begin
            E + A <--> EA             :: OnlyA
            EA + B ⇌ EPQ          :: OnlyA
            EQ + P <--> EPQ       :: EqualAI
            E + Q ⇌ EQ                :: EqualAI
        end
    end
    only_p = @enzyme_mechanism begin
        substrates: A, B ; products: P, Q ; oligomeric_state: 2
        steps: begin
            E + A <--> EA             :: EqualAI
            EA + B ⇌ EPQ          :: OnlyA
            EQ + P <--> EPQ       :: OnlyA
            E + Q ⇌ EQ                :: EqualAI
        end
    end
    only_q = @enzyme_mechanism begin
        substrates: A, B ; products: P, Q ; oligomeric_state: 2
        steps: begin
            E + A <--> EA             :: EqualAI
            EA + B ⇌ EPQ          :: OnlyA
            EQ + P <--> EPQ       :: EqualAI
            E + Q ⇌ EQ                :: OnlyA
        end
    end
    a_and_p = @enzyme_mechanism begin
        substrates: A, B ; products: P, Q ; oligomeric_state: 2
        steps: begin
            E + A <--> EA             :: OnlyA
            EA + B ⇌ EPQ          :: OnlyA
            EQ + P <--> EPQ       :: OnlyA
            E + Q ⇌ EQ                :: EqualAI
        end
    end
    a_and_q = @enzyme_mechanism begin
        substrates: A, B ; products: P, Q ; oligomeric_state: 2
        steps: begin
            E + A <--> EA             :: OnlyA
            EA + B ⇌ EPQ          :: OnlyA
            EQ + P <--> EPQ       :: EqualAI
            E + Q ⇌ EQ                :: OnlyA
        end
    end
    p_and_q = @enzyme_mechanism begin
        substrates: A, B ; products: P, Q ; oligomeric_state: 2
        steps: begin
            E + A <--> EA             :: EqualAI
            EA + B ⇌ EPQ          :: OnlyA
            EQ + P <--> EPQ       :: OnlyA
            E + Q ⇌ EQ                :: OnlyA
        end
    end
    all_three = @enzyme_mechanism begin
        substrates: A, B ; products: P, Q ; oligomeric_state: 2
        steps: begin
            E + A <--> EA             :: OnlyA
            EA + B ⇌ EPQ          :: OnlyA
            EQ + P <--> EPQ       :: OnlyA
            E + Q ⇌ EQ                :: OnlyA
        end
    end
    @test EnzymeRates.Enumeration._hyperbolic_catalysis(_testhelper_candidate(m2))
    kids = _testhelper_move(EnzymeRates.Enumeration._expand_to_allosteric, m2, rxn2)
    @test all(c -> _testhelper_is_allosteric(c), kids)
    @test length(kids) == 7
    @test Set(kids) == Set([only_a, only_p, only_q, a_and_p, a_and_q, p_and_q,
                            all_three])
    v = EnzymeRates.Derivation._guard(only_a)
    @test (v.fault, EnzymeRates.Mechanisms.step_text(only_a.catalytic, v.edge)) ==
          (EnzymeRates.Derivation._FILLED, "EA + B ⇌ EPQ")
    # The parent has 5 fitted constants, all identifiable. Each child the guard admits
    # adds L; at multiplicity 2 each conformation's binding polynomial enters squared, so
    # dividing by 1 + L no longer folds L into the other constants: 6 fitted, rank 6.
    @test _testhelper_fitted(m2) == _testhelper_identifiable_rank(m2) == 5
    for k in setdiff(kids, [only_a])
        @test _testhelper_fitted(k) == _testhelper_identifiable_rank(k) == 6
    end
    # At two subunits no variant is an L phantom.
    @test !any(
        k -> EnzymeRates.Derivation._is_phantom(EnzymeRates.Derivation._memo_law(k)), kids)
end

@testset "_expand_to_allosteric on a merged bi-bi whose B group mixes fused and plain" begin
    # The merged ordered bi-bi with B's abortive binding at EQ in B's own group:
    # EA + B ⇌ EPQ is fused, EQ + B ⇌ EBQ is plain, one RE kinetic group. A
    # group holding a chemistry step is chemistry, so the whole B group, its plain step
    # included, is :OnlyA in every K-type child, and the subsets range over the A, P
    # and Q groups: 2^3 − 1 = 7. Without the :OnlyA B group the inactive conformation
    # keeps the tree EA – E – EQ – EPQ, and EBQ lies behind B's plain
    # binding at EQ. The guard refuses one child: with A's binding :OnlyA too, EA is
    # cut off in the inactive conformation, and at zero A and B the fused binding,
    # vanishing one way into EA, fills it faster than A's steady-state binding empties
    # it (`_FILLED`). No regulator, so no V-type: seven children.
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
    end
    only_a = @_testhelper_unguarded @enzyme_mechanism begin
        substrates: A, B ; products: P, Q ; oligomeric_state: 1
        steps: begin
            E + A <--> EA                              :: OnlyA
            (EA + B ⇌ EPQ, EQ + B ⇌ EBQ)     :: OnlyA
            EQ + P <--> EPQ                        :: EqualAI
            E + Q ⇌ EQ                                 :: EqualAI
        end
    end
    only_p = @enzyme_mechanism begin
        substrates: A, B ; products: P, Q ; oligomeric_state: 1
        steps: begin
            E + A <--> EA                              :: EqualAI
            (EA + B ⇌ EPQ, EQ + B ⇌ EBQ)     :: OnlyA
            EQ + P <--> EPQ                        :: OnlyA
            E + Q ⇌ EQ                                 :: EqualAI
        end
    end
    only_q = @enzyme_mechanism begin
        substrates: A, B ; products: P, Q ; oligomeric_state: 1
        steps: begin
            E + A <--> EA                              :: EqualAI
            (EA + B ⇌ EPQ, EQ + B ⇌ EBQ)     :: OnlyA
            EQ + P <--> EPQ                        :: EqualAI
            E + Q ⇌ EQ                                 :: OnlyA
        end
    end
    a_and_p = @enzyme_mechanism begin
        substrates: A, B ; products: P, Q ; oligomeric_state: 1
        steps: begin
            E + A <--> EA                              :: OnlyA
            (EA + B ⇌ EPQ, EQ + B ⇌ EBQ)     :: OnlyA
            EQ + P <--> EPQ                        :: OnlyA
            E + Q ⇌ EQ                                 :: EqualAI
        end
    end
    a_and_q = @enzyme_mechanism begin
        substrates: A, B ; products: P, Q ; oligomeric_state: 1
        steps: begin
            E + A <--> EA                              :: OnlyA
            (EA + B ⇌ EPQ, EQ + B ⇌ EBQ)     :: OnlyA
            EQ + P <--> EPQ                        :: EqualAI
            E + Q ⇌ EQ                                 :: OnlyA
        end
    end
    p_and_q = @enzyme_mechanism begin
        substrates: A, B ; products: P, Q ; oligomeric_state: 1
        steps: begin
            E + A <--> EA                              :: EqualAI
            (EA + B ⇌ EPQ, EQ + B ⇌ EBQ)     :: OnlyA
            EQ + P <--> EPQ                        :: OnlyA
            E + Q ⇌ EQ                                 :: OnlyA
        end
    end
    all_three = @enzyme_mechanism begin
        substrates: A, B ; products: P, Q ; oligomeric_state: 1
        steps: begin
            E + A <--> EA                              :: OnlyA
            (EA + B ⇌ EPQ, EQ + B ⇌ EBQ)     :: OnlyA
            EQ + P <--> EPQ                        :: OnlyA
            E + Q ⇌ EQ                                 :: OnlyA
        end
    end
    kids = _testhelper_move(EnzymeRates.Enumeration._expand_to_allosteric, m, rxn)
    @test length(kids) == 7
    @test Set(kids) == Set([only_a, only_p, only_q, a_and_p, a_and_q, p_and_q,
                            all_three])
    v = EnzymeRates.Derivation._guard(only_a)
    @test (v.fault, EnzymeRates.Mechanisms.step_text(only_a.catalytic, v.edge)) ==
          (EnzymeRates.Derivation._FILLED, "EA + B ⇌ EPQ")
    # The parent has 5 fitted constants, all identifiable. Each child the guard admits
    # adds L, and its dead inactive conformation keeps the bindings it reaches from the
    # free enzyme through :EqualAI groups: A's where A is :EqualAI; Q's where Q is, then
    # P's at EQ where P is too. Where it keeps one, L multiplies a binding
    # polynomial that carries A or Q, which no rescaling of the parent's constants
    # absorbs: 6 fitted, rank 6. Where A and Q are both :OnlyA it binds nothing, so L
    # multiplies the free enzyme's terms alone, and 1 + L folds into A's binding rate and
    # Q's dissociation constant: 6 fitted, rank 5.
    @test _testhelper_fitted(m) == _testhelper_identifiable_rank(m) == 5
    for k in (only_p, only_q, a_and_p, p_and_q)
        @test _testhelper_fitted(k) == _testhelper_identifiable_rank(k) == 6
    end
    for k in (a_and_q, all_three)
        @test _testhelper_fitted(k) == 6 && _testhelper_identifiable_rank(k) == 5
    end
    # The two whose inactive conformation binds nothing are L phantoms, each mapped onto
    # the parent's law by its certificate; the five that keep a binding there have no
    # certificate and are fitted.
    for k in (a_and_q, all_three)
        @test EnzymeRates.Derivation._is_phantom(EnzymeRates.Derivation._memo_law(k))
        @test _testhelper_phantom_map_holds(k; rng = Random.MersenneTwister(3))
    end
    for k in (only_a, only_p, only_q, a_and_p, p_and_q)
        @test !EnzymeRates.Derivation._is_phantom(EnzymeRates.Derivation._memo_law(k))
    end
end

@testset "_expand_to_allosteric on a random-order merged bi-bi led by plain bindings" begin
    # The random-order merged bi-bi: A and B bind free E in either order, and the second
    # binding is fused with the chemistry. A's steady-state group holds E + A → EA and
    # the fused EB + A → EPQ; B's rapid-equilibrium group holds E + B ⇌ EB and
    # the fused EA + B ⇌ EPQ. Canonical order sorts a group's steps by source form,
    # E before EA and EB, so each group opens with its plain binding. A group holding
    # a chemistry step is chemistry whatever its first step, so both are :OnlyA in every
    # K-type child, and the subsets range over the Q and P groups: 2^2 − 1 = 3. Without
    # the :OnlyA A and B groups the inactive conformation keeps the path E – EQ –
    # EPQ, and EA and EB lie behind A's and B's bindings at E, so admission
    # refuses none. The reaction allows one catalytic subunit
    # and declares no regulator, so no V-type: three children. Read by its first step,
    # neither group would be chemistry, and the subsets would range over all four groups.
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
    end
    c = _testhelper_candidate(m)
    _testhelper_plain(k) = (e = m.catalytic.steps[k];
                            EnzymeRates.Enumeration._is_binding(e) &&
                                !_testhelper_is_fused(c, e))
    chemistry = [ks for ks in _testhelper_group_members(m.catalytic)
                 if !all(_testhelper_plain, ks)]
    @test length(chemistry) == 2
    @test all(ks -> _testhelper_plain(first(ks)), chemistry)
    only_q = @enzyme_mechanism begin
        substrates: A, B ; products: P, Q ; oligomeric_state: 1
        steps: begin
            (E + A <--> EA, EB + A <--> EPQ)     :: OnlyA
            (E + B ⇌ EB, EA + B ⇌ EPQ)           :: OnlyA
            E + Q ⇌ EQ                                 :: OnlyA
            EQ + P <--> EPQ                        :: EqualAI
        end
    end
    only_p = @enzyme_mechanism begin
        substrates: A, B ; products: P, Q ; oligomeric_state: 1
        steps: begin
            (E + A <--> EA, EB + A <--> EPQ)     :: OnlyA
            (E + B ⇌ EB, EA + B ⇌ EPQ)           :: OnlyA
            E + Q ⇌ EQ                                 :: EqualAI
            EQ + P <--> EPQ                        :: OnlyA
        end
    end
    p_and_q = @enzyme_mechanism begin
        substrates: A, B ; products: P, Q ; oligomeric_state: 1
        steps: begin
            (E + A <--> EA, EB + A <--> EPQ)     :: OnlyA
            (E + B ⇌ EB, EA + B ⇌ EPQ)           :: OnlyA
            E + Q ⇌ EQ                                 :: OnlyA
            EQ + P <--> EPQ                        :: OnlyA
        end
    end
    kids = _testhelper_move(EnzymeRates.Enumeration._expand_to_allosteric, m, rxn)
    @test length(kids) == 3
    @test Set(kids) == Set([only_q, only_p, p_and_q])
    # The parent has 5 fitted constants, all identifiable. Each child adds L. Its dead
    # inactive conformation reaches no EB, since B's group is :OnlyA, so it adds L·E to
    # E's rapid-equilibrium segment {E, EB, EQ}, times 1 + Q/K_Q where Q stays
    # :EqualAI. A's group takes A up at E and at EB under one rate constant, so the
    # denominator's A·B term is 2·B/K_B times its A term whatever L is: K_B is pinned and
    # cannot also absorb the 1 + L at E. L shows in every child, the two whose inactive
    # conformation binds nothing included: 6 fitted, rank 6.
    @test _testhelper_fitted(m) == _testhelper_identifiable_rank(m) == 5
    for k in (only_q, only_p, p_and_q)
        @test _testhelper_fitted(k) == _testhelper_identifiable_rank(k) == 6
    end
    # None is an L phantom: L raises every child's rank, so no certificate exists, though
    # only_q and p_and_q bind nothing in the inactive conformation.
    for k in (only_q, p_and_q)
        @test !any(
            EnzymeRates.Mechanisms.steps_running_in_allosteric_conformations(k).inactive)
    end
    @test !any(k -> EnzymeRates.Derivation._is_phantom(EnzymeRates.Derivation._memo_law(k)),
               (only_q, only_p, p_and_q))
end

@testset "K-types build on the every-chemistry set when the guard rejects it" begin
    # Census bi-bi parent 1291abb91f537c6b (census/base6), one of the 12 bi-bi parents
    # whose V-type set fails the guard. A's group mixes its plain bindings at E
    # and EP with its fused binding into EPQ, and B's mixes its plain binding at E
    # with its fused one, so both are chemistry. Tagged :OnlyA with every binding group
    # :EqualAI, the rapid-equilibrium fused binding EB + A ⇌ EPQ keeps both sides
    # joined, through Q's bindings, P's and B's binding at EQ, so it cannot be switched
    # off (`_FAST_STEP_PINNED`). No smaller set of chemistry groups stops the inactive
    # conformation's catalysis, so there is no base, and the K-types build on that set
    # itself: a tagged Q binding or B binding at EQ cuts the join. The move emits all
    # seven candidates; the guard accepts five. P's bindings alone leave B's binding at E
    # joined, so the inactive conformation stops it, and at zero P and Q EB holds enzyme
    # apart from free E (`_NONUNIQUE`). B's binding at EQ alone leaves the fused A binding
    # filling EB at zero A and B faster than EB empties (`_FILLED`).
    m = @enzyme_mechanism begin
        substrates: A, B
        products: P, Q
        steps: begin
            (E + A ⇌ EA, EP + A ⇌ EAP, EB + A ⇌ EPQ)
            (E + B <--> EB, EA + B <--> EPQ)
            (E + P <--> EP, EA + P <--> EAP, EQ + P <--> EPQ)
            (E + Q <--> EQ, EB + Q <--> EBQ, EP + Q <--> EPQ)
            EQ + B <--> EBQ
        end
    end
    rxn = @enzyme_reaction begin
        substrates: A[C], B[N]
        products: P[C], Q[N]
    end
    q_tagged = @_testhelper_unguarded @enzyme_mechanism begin
        substrates: A, B ; products: P, Q ; oligomeric_state: 1
        steps: begin
            (E + A ⇌ EA, EP + A ⇌ EAP, EB + A ⇌ EPQ)          :: OnlyA
            (E + B <--> EB, EA + B <--> EPQ)                        :: OnlyA
            (E + P <--> EP, EA + P <--> EAP, EQ + P <--> EPQ) :: EqualAI
            (E + Q <--> EQ, EB + Q <--> EBQ, EP + Q <--> EPQ) :: OnlyA
            EQ + B <--> EBQ                                           :: EqualAI
        end
    end
    b_at_eq_tagged = @_testhelper_unguarded @enzyme_mechanism begin
        substrates: A, B ; products: P, Q ; oligomeric_state: 1
        steps: begin
            (E + A ⇌ EA, EP + A ⇌ EAP, EB + A ⇌ EPQ)          :: OnlyA
            (E + B <--> EB, EA + B <--> EPQ)                        :: OnlyA
            (E + P <--> EP, EA + P <--> EAP, EQ + P <--> EPQ) :: EqualAI
            (E + Q <--> EQ, EB + Q <--> EBQ, EP + Q <--> EPQ) :: EqualAI
            EQ + B <--> EBQ                                           :: OnlyA
        end
    end
    p_q_tagged = @_testhelper_unguarded @enzyme_mechanism begin
        substrates: A, B ; products: P, Q ; oligomeric_state: 1
        steps: begin
            (E + A ⇌ EA, EP + A ⇌ EAP, EB + A ⇌ EPQ)          :: OnlyA
            (E + B <--> EB, EA + B <--> EPQ)                        :: OnlyA
            (E + P <--> EP, EA + P <--> EAP, EQ + P <--> EPQ) :: OnlyA
            (E + Q <--> EQ, EB + Q <--> EBQ, EP + Q <--> EPQ) :: OnlyA
            EQ + B <--> EBQ                                           :: EqualAI
        end
    end
    p_b_at_eq_tagged = @_testhelper_unguarded @enzyme_mechanism begin
        substrates: A, B ; products: P, Q ; oligomeric_state: 1
        steps: begin
            (E + A ⇌ EA, EP + A ⇌ EAP, EB + A ⇌ EPQ)          :: OnlyA
            (E + B <--> EB, EA + B <--> EPQ)                        :: OnlyA
            (E + P <--> EP, EA + P <--> EAP, EQ + P <--> EPQ) :: OnlyA
            (E + Q <--> EQ, EB + Q <--> EBQ, EP + Q <--> EPQ) :: EqualAI
            EQ + B <--> EBQ                                           :: OnlyA
        end
    end
    q_b_at_eq_tagged = @_testhelper_unguarded @enzyme_mechanism begin
        substrates: A, B ; products: P, Q ; oligomeric_state: 1
        steps: begin
            (E + A ⇌ EA, EP + A ⇌ EAP, EB + A ⇌ EPQ)          :: OnlyA
            (E + B <--> EB, EA + B <--> EPQ)                        :: OnlyA
            (E + P <--> EP, EA + P <--> EAP, EQ + P <--> EPQ) :: EqualAI
            (E + Q <--> EQ, EB + Q <--> EBQ, EP + Q <--> EPQ) :: OnlyA
            EQ + B <--> EBQ                                           :: OnlyA
        end
    end
    all_tagged = @_testhelper_unguarded @enzyme_mechanism begin
        substrates: A, B ; products: P, Q ; oligomeric_state: 1
        steps: begin
            (E + A ⇌ EA, EP + A ⇌ EAP, EB + A ⇌ EPQ)          :: OnlyA
            (E + B <--> EB, EA + B <--> EPQ)                        :: OnlyA
            (E + P <--> EP, EA + P <--> EAP, EQ + P <--> EPQ) :: OnlyA
            (E + Q <--> EQ, EB + Q <--> EBQ, EP + Q <--> EPQ) :: OnlyA
            EQ + B <--> EBQ                                           :: OnlyA
        end
    end
    p_tagged = @_testhelper_unguarded @enzyme_mechanism begin
        substrates: A, B ; products: P, Q ; oligomeric_state: 1
        steps: begin
            (E + A ⇌ EA, EP + A ⇌ EAP, EB + A ⇌ EPQ)          :: OnlyA
            (E + B <--> EB, EA + B <--> EPQ)                        :: OnlyA
            (E + P <--> EP, EA + P <--> EAP, EQ + P <--> EPQ) :: OnlyA
            (E + Q <--> EQ, EB + Q <--> EBQ, EP + Q <--> EPQ) :: EqualAI
            EQ + B <--> EBQ                                           :: EqualAI
        end
    end
    function _testhelper_admission_fault(t)
        am = EnzymeRates.Mechanism(m.catalytic, EnzymeRates.Mechanisms.Allostery(t, 1,
                                             EnzymeRates.Mechanisms.RegulatorySite[]))
        EnzymeRates.Derivation._verdict_under(
            am, EnzymeRates.Derivation._built_reading(am)).fault
    end
    vt = [c ? :OnlyA : :EqualAI for c in EnzymeRates.Enumeration._chemistry_groups(
        _testhelper_candidate(m))]
    @test _testhelper_admission_fault(vt) === EnzymeRates.Derivation._FAST_STEP_PINNED
    @test isempty(EnzymeRates.Enumeration._tagging_bases(_testhelper_candidate(m)))
    kids = _testhelper_move(EnzymeRates.Enumeration._expand_to_allosteric, m, rxn)
    @test length(kids) == 7
    @test Set(kids) == Set([p_tagged, q_tagged, b_at_eq_tagged, p_q_tagged,
                            p_b_at_eq_tagged, q_b_at_eq_tagged, all_tagged])
    @test all(_testhelper_is_valid, setdiff(kids, [p_tagged, b_at_eq_tagged]))
    @test EnzymeRates.Derivation._guard(p_tagged).fault ===
          EnzymeRates.Derivation._NONUNIQUE
    @test EnzymeRates.Derivation._guard(b_at_eq_tagged).fault ===
          EnzymeRates.Derivation._FILLED
end

@testset "a ter-bi ping-pong's one base is the smallest stopping set the guard admits" begin
    # Census ter-bi parent 0dc399e7dab6e4f5 (census/base6): its covalent form E(; residual
    # = B + D - Q) takes up A and gives off P. Four groups hold chemistry, those of B, D, A
    # and Q, each mixing a plain binding with a fused one; P's binding at E is the only
    # binding group. Tagging every chemistry group :OnlyA fails the guard: at zero A, B and
    # D the fused A binding, vanishing one way into the covalent form, fills it faster than
    # it empties (`_FILLED`). So the bases are the smallest stopping sets that pass. Q's
    # group alone passes and is the one base. B's alone passes the guard but keeps EB,
    # which only the binding it excludes enters (`_kept_forms_keep_a_binding`); D's alone
    # and A's alone each leave the sides of a rapid-equilibrium binding joined. Each K-type
    # adds P's binding to the base or to the every-chemistry set, and the guard admits
    # both. The reaction declares no regulator, so there is no V-type: the exits keep the
    # two K-types.
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
    base_q = @enzyme_mechanism begin
        substrates: A, B, D ; products: P, Q ; oligomeric_state: 1
        steps: begin
            (E + B <--> EB, ED + B <--> EQ_res_plusB_plusD_minusQ)      :: EqualAI
            (E + D ⇌ ED, EB + D ⇌ EQ_res_plusB_plusD_minusQ)            :: EqualAI
            E + P <--> EP                                                  :: EqualAI
            (E_res_plusB_plusD_minusQ + A ⇌ EP,
             EQ_res_plusB_plusD_minusQ + A ⇌ EAQ_res_plusB_plusD_minusQ) :: EqualAI
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
    p_on_q = @enzyme_mechanism begin
        substrates: A, B, D ; products: P, Q ; oligomeric_state: 1
        steps: begin
            (E + B <--> EB, ED + B <--> EQ_res_plusB_plusD_minusQ)      :: EqualAI
            (E + D ⇌ ED, EB + D ⇌ EQ_res_plusB_plusD_minusQ)            :: EqualAI
            E + P <--> EP                                                  :: OnlyA
            (E_res_plusB_plusD_minusQ + A ⇌ EP,
             EQ_res_plusB_plusD_minusQ + A ⇌ EAQ_res_plusB_plusD_minusQ) :: EqualAI
            (E_res_plusB_plusD_minusQ + Q <--> EQ_res_plusB_plusD_minusQ,
             EP + Q <--> EAQ_res_plusB_plusD_minusQ)                    :: OnlyA
        end
    end
    p_on_every = @enzyme_mechanism begin
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
    D = EnzymeRates.Derivation
    tags = EnzymeRates.Enumeration._tagging_bases(_testhelper_candidate(m))
    @test length(tags) == 1
    @test EnzymeRates.Mechanism(
        m.catalytic, EnzymeRates.Mechanisms.Allostery(
            only(tags), 1, EnzymeRates.Mechanisms.RegulatorySite[])) == base_q
    vt = [c ? :OnlyA : :EqualAI for c in EnzymeRates.Enumeration._chemistry_groups(
        _testhelper_candidate(m))]
    v = D._guard(EnzymeRates.Mechanism(
        m.catalytic,
        EnzymeRates.Mechanisms.Allostery(vt, 1, EnzymeRates.Mechanisms.RegulatorySite[])))
    @test (v.fault, EnzymeRates.Mechanisms.step_text(m.catalytic, v.edge)) ==
          (D._FILLED, "E_res_plusB_plusD_minusQ + A ⇌ EP")
    @test _testhelper_is_valid(b_alone) &&
          !EnzymeRates.Enumeration._kept_forms_keep_a_binding(b_alone)
    kids = _testhelper_move(EnzymeRates.Enumeration._expand_to_allosteric, m, rxn)
    @test length(kids) == 2
    @test Set(kids) == Set([p_on_q, p_on_every])
    # P's binding :OnlyA on the every-chemistry set gives a law both oracles find to be
    # the limit on every face.
    @test all(_testhelper_admission_agreement(p_on_every;
                                              rng = Random.MersenneTwister(13964)))
    _testhelper_check_law_against_limit(p_on_every)
    @test Set(c for c in _testhelper_children(rxn, [m])
              if _testhelper_is_allosteric(c)) == Set(kids)
end

@testset "change_allo_state emits each relaxation; the exit drops one that is pinned" begin
    # A K-type of census bi-bi parent 1291abb91f537c6b (testset "K-types build on the
    # every-chemistry set when the guard rejects it"): A's and B's mixed groups and Q's
    # bindings :OnlyA. Relaxing Q's bindings to :NonequalAI gives back the every-chemistry
    # set, whose rapid-equilibrium fused binding EB + A ⇌ EPQ has its two sides
    # joined in the inactive conformation, through Q's bindings, P's and B's binding at
    # EQ, so it cannot be switched off (`_FAST_STEP_PINNED`), and the exit drops it.
    # Relaxing P's bindings or B's binding at EQ leaves Q's cut, and the guard accepts
    # both. The chemistry does not relax while a binding is :OnlyA. Three children.
    q_tagged = @_testhelper_unguarded @enzyme_mechanism begin
        substrates: A, B ; products: P, Q ; oligomeric_state: 1
        steps: begin
            (E + A ⇌ EA, EP + A ⇌ EAP, EB + A ⇌ EPQ)          :: OnlyA
            (E + B <--> EB, EA + B <--> EPQ)                        :: OnlyA
            (E + P <--> EP, EA + P <--> EAP, EQ + P <--> EPQ) :: EqualAI
            (E + Q <--> EQ, EB + Q <--> EBQ, EP + Q <--> EPQ) :: OnlyA
            EQ + B <--> EBQ                                           :: EqualAI
        end
    end
    p_relaxed = @_testhelper_unguarded @enzyme_mechanism begin
        substrates: A, B ; products: P, Q ; oligomeric_state: 1
        steps: begin
            (E + A ⇌ EA, EP + A ⇌ EAP, EB + A ⇌ EPQ)          :: OnlyA
            (E + B <--> EB, EA + B <--> EPQ)                        :: OnlyA
            (E + P <--> EP, EA + P <--> EAP, EQ + P <--> EPQ) :: NonequalAI
            (E + Q <--> EQ, EB + Q <--> EBQ, EP + Q <--> EPQ) :: OnlyA
            EQ + B <--> EBQ                                           :: EqualAI
        end
    end
    b_at_eq_relaxed = @_testhelper_unguarded @enzyme_mechanism begin
        substrates: A, B ; products: P, Q ; oligomeric_state: 1
        steps: begin
            (E + A ⇌ EA, EP + A ⇌ EAP, EB + A ⇌ EPQ)          :: OnlyA
            (E + B <--> EB, EA + B <--> EPQ)                        :: OnlyA
            (E + P <--> EP, EA + P <--> EAP, EQ + P <--> EPQ) :: EqualAI
            (E + Q <--> EQ, EB + Q <--> EBQ, EP + Q <--> EPQ) :: OnlyA
            EQ + B <--> EBQ                                           :: NonequalAI
        end
    end
    q_relaxed = @_testhelper_unguarded @enzyme_mechanism begin
        substrates: A, B ; products: P, Q ; oligomeric_state: 1
        steps: begin
            (E + A ⇌ EA, EP + A ⇌ EAP, EB + A ⇌ EPQ)          :: OnlyA
            (E + B <--> EB, EA + B <--> EPQ)                        :: OnlyA
            (E + P <--> EP, EA + P <--> EAP, EQ + P <--> EPQ) :: EqualAI
            (E + Q <--> EQ, EB + Q <--> EBQ, EP + Q <--> EPQ) :: NonequalAI
            EQ + B <--> EBQ                                           :: EqualAI
        end
    end
    kids = _testhelper_move(EnzymeRates.Enumeration._expand_change_allo_state, q_tagged)
    @test length(kids) == 3
    @test Set(kids) == Set([p_relaxed, b_at_eq_relaxed, q_relaxed])
    @test EnzymeRates.Derivation._guard(q_relaxed).fault ===
          EnzymeRates.Derivation._FAST_STEP_PINNED
    @test all(_testhelper_is_valid, [p_relaxed, b_at_eq_relaxed])
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
    # both. The other three children tag the seed's groups (no oligomeric state is
    # declared, so they keep one catalytic subunit): the substrate binding, the product
    # binding, or both, :OnlyA with the chemistry.
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
    @test length(children) == 4
    @test Set(children) == Set([dead_end, iinh_only_a, p_only_a, both_only_a])
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

# ═══════════════════════════════════════════════════════════════════════
# Rate-equation dedup key
# ═══════════════════════════════════════════════════════════════════════

@testset "Rate-equation dedup key" begin
    @testset "Distinct mechanisms produce distinct dedup keys" begin
        # Ordered binding (A-first only) vs random binding (both A-first
        # and B-first paths). Same substrate set, same product set, but
        # the random mechanism has a different enzyme-form graph. The
        # dedup key must differ.
        m_ordered = @enzyme_mechanism begin
            substrates: A, B
            products: P
            steps: begin
                E + A ⇌ EA
                EA + B ⇌ EAB
                EAB <--> EP
                E + P ⇌ EP
            end
        end
        m_random = @enzyme_mechanism begin
            substrates: A, B
            products: P
            steps: begin
                E + A ⇌ EA
                E + B ⇌ EB
                EA + B ⇌ EAB
                EB + A ⇌ EAB
                EAB <--> EP
                E + P ⇌ EP
            end
        end
        key(m) = EnzymeRates.Derivation._canonical(
            EnzymeRates.Derivation._memo_law(m), m).key
        @test key(m_ordered) != key(m_random)
    end

    @testset "LDH Pattern-A: graph-distinct mechanisms of equal v share a dedup key" begin
        # LDH Pattern-A pair (11 steps each). m_a and m_b differ in step
        # ordering AND in which intermediate forms appear (m_a has
        # Lactate-binding via E_NADH only; m_b adds a Lactate-binding via
        # E_NAD path). After Pass-1 absorption their v polynomials are
        # equivalent. 11 steps is the minimal known case for this property
        # — smaller mechanisms produce graph-equivalent topologies
        # (already collapsed by step sorting) rather than graph-distinct
        # yet v-equivalent ones.
        ldh_m_a = @enzyme_mechanism begin
            substrates: NADH, Pyruvate
            products: Lactate, NAD
            steps: begin
                (E + Lactate ⇌ ELactate, ENADH + Lactate ⇌ ELactateNADH)
                (E + NAD ⇌ ENAD, ELactate + NAD ⇌ ELactateNAD,
                 EPyruvate + NAD ⇌ ENADPyruvate)
                (E + NADH ⇌ ENADH, ELactate + NADH ⇌ ELactateNADH)
                (E + Pyruvate ⇌ EPyruvate, ENAD + Pyruvate ⇌ ENADPyruvate,
                 ENADH + Pyruvate ⇌ ENADHPyruvate)
                ENADHPyruvate <--> ELactateNAD
            end
        end

        ldh_m_b = @enzyme_mechanism begin
            substrates: NADH, Pyruvate
            products: Lactate, NAD
            steps: begin
                (E + Lactate ⇌ ELactate, ENAD + Lactate ⇌ ELactateNAD,
                 ENADH + Lactate ⇌ ELactateNADH)
                (E + NAD ⇌ ENAD, EPyruvate + NAD ⇌ ENADPyruvate)
                (E + NADH ⇌ ENADH, ELactate + NADH ⇌ ELactateNADH)
                (E + Pyruvate ⇌ EPyruvate, ENAD + Pyruvate ⇌ ENADPyruvate,
                 ENADH + Pyruvate ⇌ ENADHPyruvate)
                ENADHPyruvate <--> ELactateNAD
            end
        end

        law_a = EnzymeRates.Derivation._memo_law(ldh_m_a)
        law_b = EnzymeRates.Derivation._memo_law(ldh_m_b)
        canon_a = EnzymeRates.Derivation._canonical(law_a, ldh_m_a)
        canon_b = EnzymeRates.Derivation._canonical(law_b, ldh_m_b)
        key_a = canon_a.key
        @test key_a == canon_b.key
        @test _testhelper_same_law_under_map(law_a, law_b, canon_a.order, canon_b.order;
                                             rng = Random.MersenneTwister(13))

        # Negative control: ldh_m_b with the E_NAD+Lactate step's
        # kinetic_group changed from 1 to 12 breaks the Lactate-
        # binding-via-NAD path's sharing with the direct E+Lactate
        # binding step. The resulting v polynomial differs.
        ldh_m_c = @enzyme_mechanism begin
            substrates: NADH, Pyruvate
            products: Lactate, NAD
            steps: begin
                (E + Lactate ⇌ ELactate, ENADH + Lactate ⇌ ELactateNADH)
                (E + NAD ⇌ ENAD, EPyruvate + NAD ⇌ ENADPyruvate)
                (E + NADH ⇌ ENADH, ELactate + NADH ⇌ ELactateNADH)
                (E + Pyruvate ⇌ EPyruvate, ENAD + Pyruvate ⇌ ENADPyruvate,
                 ENADH + Pyruvate ⇌ ENADHPyruvate)
                ENAD + Lactate ⇌ ELactateNAD
                ENADHPyruvate <--> ELactateNAD
            end
        end
        @test key_a != EnzymeRates.Derivation._canonical(
            EnzymeRates.Derivation._memo_law(ldh_m_c), ldh_m_c).key
    end

    @testset "rate_equation_string emits section labels" begin
        # Minimal uni-uni 3-step mechanism.
        uni_uni_3step = @enzyme_mechanism begin
            substrates: S
            products: P
            steps: begin
                E + S ⇌ ES
                ES <--> EP
                E + P ⇌ EP
            end
        end

        # Random bi-uni with substrate-side mirror sharing: A-binding
        # steps share kinetic_group=1, B-binding steps share kg=2.
        biuni_mirror = @enzyme_mechanism begin
            substrates: A, B
            products: P
            steps: begin
                (E + A ⇌ EA, EB + A ⇌ EAB)
                (E + B ⇌ EB, EA + B ⇌ EAB)
                EAB <--> EP
                E + P ⇌ EP
            end
        end

        # With structural names, shared kinetic_group members collapse
        # via the value-context chokepoint — no separate user-defined
        # section is emitted.
        s_user = rate_equation_string(biuni_mirror)
        @test !occursin("# User defined constraints:", s_user)

        # Haldane section: any RE binding mechanism with Keq has it.
        # The Wegscheider section is not emitted on minimal mechanisms;
        # the LDH Pattern-A test above is the indirect regression for
        # that section's stripping.
        s_hal = rate_equation_string(uni_uni_3step)
        @test occursin("# Haldane constraints:", s_hal)

        # Wegscheider section: emitted when the thermodynamic constraint
        # system produces cycle equalities. Random bi-bi with all-
        # singleton kinetic groups (no Pass-1 absorption) preserves the
        # Wegscheider cycle relations and renders them as multi-symbol
        # RHSes in this section.
        m_weg = @enzyme_mechanism begin
            substrates: A, B
            products: P, Q
            steps: begin
                E + A ⇌ EA
                E + B ⇌ EB
                EA + B ⇌ EAB
                EB + A ⇌ EAB
                EAB <--> EPQ
                EP + Q ⇌ EPQ
                EQ + P ⇌ EPQ
                E + P ⇌ EP
                E + Q ⇌ EQ
            end
        end
        s_weg = rate_equation_string(m_weg)
        @test occursin("# Wegscheider constraints:", s_weg)
    end
end

@testset "bi-bi init_mechanisms structural golden" begin
    # Canonical, derivation-free key for one Mechanism: per kinetic group,
    # the sorted set of (from form, to form, bound metabolite, RE/SS) tuples;
    # groups themselves sorted so the key is order-independent.
    function _testhelper_mech_struct_key(m::EnzymeRates.Mechanism)
        g = _testhelper_graph(m)
        grpkeys = [join(sort([
                       string((Symbol(g.forms[s.from]), Symbol(g.forms[s.to]),
                               _testhelper_is_binding(s) ? only(s.consumed) : :iso,
                               s.rapid_equilibrium))
                       for s in g.steps if s.group == k]), "|")
                   for k in unique(s.group for s in g.steps)]
        join(sort(grpkeys), " ;; ")
    end

    init = _testhelper_mechanisms(EnzymeRates.Enumeration.init_mechanisms(bi_bi_rxn))
    mech_keys = sort([_testhelper_mech_struct_key(m) for m in init])
    @test length(unique(mech_keys)) == length(mech_keys)   # no structural dups post-dedup

    fixture = joinpath(@__DIR__, "..", "fixtures", "phase2_init_golden.txt")
    if !isfile(fixture)
        mkpath(dirname(fixture)); write(fixture, join(mech_keys, "\n"))
        @warn "bootstrapped phase2 golden fixture; commit it" fixture
    end
    golden = readlines(fixture)
    @test mech_keys == golden          # permanent structural regression gate
end

@testset "init division-freeness (bi_bi_pp)" begin
    # Every enumerated init mechanism's derived rate equation must stay finite
    # when any single metabolite concentration is zero (real data has zeros).
    # `isfinite` suffices here: a residual coupling that survived the
    # concentration-GCD would put the same 1/conc factor in BOTH numerator and
    # denominator, so a zeroed metabolite yields Inf/Inf = NaN, which isfinite
    # catches (no separate != 0 check needed).
    mets = [:A, :B, :P, :Q]
    for m in unique!(collect(
        _testhelper_mechanisms(EnzymeRates.Enumeration.init_mechanisms(bi_bi_pp_rxn))))
        params, Keq =
            _testhelper_random_independent_params(m; rng = Random.MersenneTwister(1))
        rate = compile_rate_equation(m)
        for zeroed in mets
            cvals = Tuple(n == zeroed ? 0.0 : 1.0 for n in mets)
            concs = NamedTuple{Tuple(mets)}(cvals)
            v = rate(concs, params, Keq)
            @test isfinite(v)
        end
    end
end

@testset "to_allosteric shares :EqualAI downstream constants across conformations" begin
    # The inactive conformation drops each `:OnlyA` binding step, and its law holds only
    # the forms its untagged steps join to free E: that ligand's downstream complex holds
    # inactive enzyme only if other untagged steps join it to free E.
    # A lumped `:EqualAI` group is one class in both conformations, named after the
    # kept member of the full mechanism, so the inactive conformation cannot give it a
    # second name and the joint solve counts it once.
    #
    # On these seeds, whose every-chemistry set admission accepts, the K-types of
    # `_expand_to_allosteric` are dead-inactive `:OnlyA`-binding combos (all chemical
    # steps `:OnlyA`), so the case this guards is the child that sets the free-enzyme
    # bindings (a substrate/product pair on free E) `:OnlyA`, leaving the lumped
    # downstream group's binding `:EqualAI`, outside the inactive law, which holds only
    # the forms its untagged steps join to free E; the constant stays merged, so the child
    # adds only `L` — parent + 1 fitted params. A child that instead sets a downstream
    # binding `:OnlyA` orphans a bound form of the lumped group; its multi-`:OnlyA`
    # derivation is a separate, documented concern and is not asserted here.
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

@testset "seed_mechanisms" begin
    # Tag-stripped skeleton: normalize every regulator state to :OnlyA so the
    # 2^n one-ligand-site state assignments of a lineage collapse to one key.
    # The catalytic backbone and allosteric_states are preserved.
    _testhelper_skeleton_key(m) = hash(EnzymeRates.Mechanism(
        m.catalytic, EnzymeRates.Mechanisms.Allostery(m.allosteric.allosteric_states,
        m.allosteric.oligomeric_state,
        [EnzymeRates.Mechanisms.RegulatorySite(
            s.ligands, s.multiplicity,
            fill(:OnlyA, length(s.allosteric_states)))
                 for s in m.allosteric.regulatory_sites])))

    # The allosteric state of ligand `reg` wherever it sits, or :none.
    function _testhelper_state_of(m, reg)
        for site in m.allosteric.regulatory_sites
            for (lig, st) in zip(site.ligands,
                                 site.allosteric_states)
                lig == reg && return st
            end
        end
        :none
    end

    # The seeds of a reaction with two undesignated allosteric regulators, both required.
    # Their names differ from the `rxn` and `seeds` that other subtests assign.
    rxn_xy = @enzyme_reaction begin
        substrates: A[C6H12O6]
        products:   B[C6H12O6]
        allosteric_regulators: X, Y
        oligomeric_state: 2
    end
    seeds_xy = _testhelper_mechanisms(
        EnzymeRates.Enumeration.seed_mechanisms(rxn_xy, Set([:X, :Y]), Set{Symbol}()))

    @testset "undesignated: skeletons × 2²" begin
        @test !isempty(seeds_xy)
        # Every seed is a fully-regulated, cheap-state allosteric mechanism with X
        # and Y each at their own single-ligand site.
        for s in seeds_xy
            @test _testhelper_is_allosteric(s)
            @test EnzymeRates.Enumeration._bound_allo_regs(s) == Set([:X, :Y])
            sites = s.allosteric.regulatory_sites
            @test length(sites) == 2
            @test all(length(si.ligands) == 1 for si in sites)
            @test :NonequalAI ∉ [
                s.allosteric.allosteric_states;
                [st for si in sites for st in si.allosteric_states]]
        end
        # Each lineage contributes exactly 2² = 4 state assignments.
        skels = Dict{UInt64, Int}()
        for s in seeds_xy
            k = _testhelper_skeleton_key(s)
            skels[k] = get(skels, k, 0) + 1
        end
        @test all(==(4), values(skels))
        @test length(seeds_xy) == 4 * length(skels)
    end

    @testset "designated types collapse to ×1" begin
        rxn_d = @enzyme_reaction begin
            substrates: A[C6H12O6]
            products:   B[C6H12O6]
            allosteric_regulators: X::Activator, Y::Inhibitor
            oligomeric_state: 2
        end
        seeds_d =
            _testhelper_mechanisms(EnzymeRates.Enumeration.seed_mechanisms(
                rxn_d, Set([:X, :Y]), Set{Symbol}()))
        @test !isempty(seeds_d)
        # One state assignment per lineage: seed count equals the skeleton count.
        skels_d = Set(_testhelper_skeleton_key(s) for s in seeds_d)
        @test length(seeds_d) == length(skels_d)
        # The same catalytic-allostery skeletons as the undesignated build.
        @test length(seeds_xy) == 4 * length(seeds_d)
        # An activator seeds as :OnlyA, an inhibitor as :OnlyI.
        for s in seeds_d
            @test _testhelper_state_of(s, :X) == :OnlyA
            @test _testhelper_state_of(s, :Y) == :OnlyI
        end
    end

    @testset "required vs optional" begin
        rxn = @enzyme_reaction begin
            substrates: A[C6H12O6]
            products:   B[C6H12O6]
            allosteric_regulators: X, Y, Z
            oligomeric_state: 2
        end
        seeds = _testhelper_mechanisms(
            EnzymeRates.Enumeration.seed_mechanisms(rxn, Set([:X, :Y]), Set{Symbol}()))
        @test !isempty(seeds)
        for s in seeds
            @test issubset(Set([:X, :Y]), EnzymeRates.Enumeration._bound_allo_regs(s))
            @test !(:Z in EnzymeRates.Enumeration._bound_allo_regs(s))
        end
    end

    @testset "per-lineage floor invariant" begin
        n_required = 2
        # A required allosteric regulator adds one dissociation constant; the two
        # required regulators additionally lift the mechanism to two
        # conformations, so a seed carries its base catalytic count + n_required
        # + 1 (L). Compile is slow; sample a few seeds.
        for s in seeds_xy[1:min(3, length(seeds_xy))]
            base = _testhelper_plain_of(s)
            base_count = _testhelper_fitted(base)
            seed_count = _testhelper_fitted(s)
            @test seed_count == base_count + n_required + 1
        end
    end

    @testset "competitive-only: non-allosteric seeds, no L" begin
        rxn = @enzyme_reaction begin
            substrates: A[C6H12O6]
            products:   B[C6H12O6]
            competitive_inhibitors: I
        end
        seeds = _testhelper_mechanisms(
            EnzymeRates.Enumeration.seed_mechanisms(rxn, Set{Symbol}(), Set([:I])))
        @test !isempty(seeds)
        # No allosteric-lifting move fires: every seed is a plain mechanism that
        # binds I as a dead-end competitive binding, and carries no L.
        for s in seeds
            @test !_testhelper_is_allosteric(s)
            @test :I in EnzymeRates.Enumeration._bound_comp_inhibitors(
                _testhelper_candidate(s))
        end
        # Floor is base + n_required_comp (no L). Uni-uni has one catalytic
        # backbone, so every seed shares one base lineage.
        inits = _testhelper_mechanisms(EnzymeRates.Enumeration.init_mechanisms(rxn))
        @test length(inits) == 1
        base_count = _testhelper_fitted(inits[1])
        for s in seeds[1:min(2, length(seeds))]
            fp = parameters(s)
            @test length(fp) == base_count + 1
            @test !(:L in fp)
        end
    end
end

@testset "seed_mechanisms errors when no mechanism binds every required regulator" begin
    # Uni-uni with S as its own competitive inhibitor: the only site that carries
    # neither S nor P is free E, where the copy duplicates ES. No seed exists,
    # and the beam must say so rather than return nothing.
    rxn = @enzyme_reaction begin
        substrates: S[C]
        products: P[C]
        competitive_inhibitors: S
    end
    err = try
        _testhelper_mechanisms(
            EnzymeRates.Enumeration.seed_mechanisms(rxn, Set{Symbol}(), Set([:S]))); nothing
    catch e
        e
    end
    @test err isa ErrorException
    # The user calls `candidate_mechanisms`, so the message names it.
    @test startswith(err.msg, "candidate_mechanisms: no mechanism binds every required " *
                              "regulator")
    @test occursin("required regulator", err.msg)
    @test occursin("competitive inhibitors: S", err.msg)
    @test occursin("allosteric regulators: none", err.msg)
    @test occursin("a uni-uni seed has no such site", err.msg)
    @test occursin("optional_competitive_inhibitors", err.msg)
    # Bi-bi with A as its own inhibitor has placements that create a new complex.
    bibi = @enzyme_reaction begin
        substrates: A[C], B[N]
        products: P[C], Q[N]
        competitive_inhibitors: A
    end
    seeds = _testhelper_mechanisms(
        EnzymeRates.Enumeration.seed_mechanisms(bibi, Set{Symbol}(), Set([:A])))
    @test !isempty(seeds)
    @test all(m -> EnzymeRates.Enumeration._assert_emission_rules(
        _testhelper_candidate(m)) === nothing, seeds)
    @test all(m -> :A in EnzymeRates.Enumeration._bound_comp_inhibitors(
        _testhelper_candidate(m)), seeds)
end

@testset "seed_mechanisms grows the merged and Theorell–Chance seeds" begin
    ER = EnzymeRates
    # A required dead-end inhibitor is placed on every init mechanism, the merged and
    # Theorell–Chance variants included, so some seed carries a fused binding.
    rxn = @enzyme_reaction begin
        substrates: A[C], B[N]
        products: P[C], Q[N]
        competitive_inhibitors: I
    end
    seeds = _testhelper_mechanisms(
        ER.Enumeration.seed_mechanisms(rxn, Set{Symbol}(), Set([:I])))
    _testhelper_fused(m) =
        (c = _testhelper_candidate(m);
         any(e -> _testhelper_is_fused(c, e), m.catalytic.steps))
    @test all(m -> :I in ER.Enumeration._bound_comp_inhibitors(_testhelper_candidate(m)),
              seeds)
    @test any(_testhelper_fused, seeds)
end

@testset "seed_mechanisms wave-parallel equivalence" begin
    # Inline serial FIFO BFS = the reference the parallel version must match.
    _testhelper_binds_all_required(m, req_allo, req_comp) =
        issubset(req_allo, EnzymeRates.Enumeration._bound_allo_regs(m)) &&
        issubset(req_comp,
                 EnzymeRates.Enumeration._bound_comp_inhibitors(_testhelper_candidate(m)))
    function _testhelper_serial_seed_reference(rxn, req_allo, req_comp)
        visited = Set{UInt64}()
        queue = EnzymeRates.Mechanism[]
        seeds = EnzymeRates.Mechanism[]
        _testhelper_enq(m) = begin
            h = hash(m)
            h in visited && return
            push!(visited, h); push!(queue, m)
            _testhelper_binds_all_required(m, req_allo, req_comp) && push!(seeds, m)
        end
        for m in _testhelper_mechanisms(
            EnzymeRates.Enumeration.init_mechanisms(rxn)); _testhelper_enq(m); end
        while !isempty(queue)
            m = popfirst!(queue)
            lifted = isempty(req_allo) ? [] : !_testhelper_is_allosteric(m) ?
                _testhelper_move(EnzymeRates.Enumeration._expand_to_allosteric, m, rxn) :
                _testhelper_move(EnzymeRates.Enumeration._expand_add_allosteric_regulator,
                                 m, rxn)
            for c in [lifted;
                      _testhelper_move(
                          EnzymeRates.Enumeration._expand_add_dead_end_regulator, m, rxn)]
                EnzymeRates.Enumeration._is_seed_node(_testhelper_candidate(c), rxn,
                                                      req_allo, req_comp) &&
                    _testhelper_enq(c)
            end
        end
        seeds
    end

    req = Set([:R])
    empty = Set{Symbol}()
    got = _testhelper_mechanisms(
        EnzymeRates.Enumeration.seed_mechanisms(uni_uni_allo_reg, req, empty))
    ref = _testhelper_serial_seed_reference(uni_uni_allo_reg, req, empty)

    @test got == ref                                   # same seeds, same order
    @test !isempty(got)                                # the case is non-trivial
    @test allunique(hash.(got))                        # no duplicate structures
    @test all(m -> _testhelper_binds_all_required(m, req, empty), got)
    # deterministic
    @test _testhelper_mechanisms(
        EnzymeRates.Enumeration.seed_mechanisms(uni_uni_allo_reg, req, empty)) == got
end

@testset "seed_mechanisms drops a seed the guard rejects" begin
    ER = EnzymeRates
    rxn = @enzyme_reaction begin
        substrates: A[CX], B[N]
        products: P[C], Q[NX]
        allosteric_regulators: R
        oligomeric_state: 2
        shared_catalytic_site: (A, P)
    end
    seeds = _testhelper_mechanisms(
        ER.Enumeration.seed_mechanisms(rxn, Set([:R]), Set{Symbol}()))
    # Aggregate regression pin over the whole seed set: the search reaches 1,728 seeds and
    # the guard drops 90 of them, 80 `_PUMPED` and 10 `_FILLED`.
    @test length(seeds) == 1638
    @test all(_testhelper_is_valid, seeds)
    # The six-step ping-pong with its first half-reaction tagged :OnlyA fails class 3 at
    # zero A and B, so it is dropped. With its second half tagged instead, the inactive
    # conformation keeps the two sides of the rapid-equilibrium step
    # EB_res_plusA_minusP ⇌ EQ joined (`_FAST_STEP_PINNED`), so it is no seed either.
    # With both halves tagged it stays.
    half1 = @_testhelper_unguarded @enzyme_mechanism begin
        substrates: A, B
        products: P, Q
        oligomeric_state: 2
        steps: begin
            E + A ⇌ EA                                        :: EqualAI
            EA <--> EP_res_plusA_minusP                    :: OnlyA
            EP_res_plusA_minusP ⇌ E_res_plusA_minusP + P  :: EqualAI
            E_res_plusA_minusP + B ⇌ EB_res_plusA_minusP  :: EqualAI
            EB_res_plusA_minusP ⇌ EQ                       :: EqualAI
            EQ ⇌ E + Q                                        :: EqualAI
        end
        regulatory_site(2): R::OnlyA
    end
    v = ER.Derivation._guard(half1)
    @test v.fault === ER.Derivation._SPLIT
    @test ER.Derivation._face_names(half1.catalytic,
                                    ER.Derivation._faces(half1.catalytic)[v.face]) ==
          [:A, :B]
    @test half1 ∉ seeds
    half2 = @_testhelper_unguarded @enzyme_mechanism begin
        substrates: A, B
        products: P, Q
        oligomeric_state: 2
        steps: begin
            E + A ⇌ EA                                        :: EqualAI
            EA <--> EP_res_plusA_minusP                    :: EqualAI
            EP_res_plusA_minusP ⇌ E_res_plusA_minusP + P  :: EqualAI
            E_res_plusA_minusP + B ⇌ EB_res_plusA_minusP  :: EqualAI
            EB_res_plusA_minusP ⇌ EQ                       :: OnlyA
            EQ ⇌ E + Q                                        :: EqualAI
        end
        regulatory_site(2): R::OnlyA
    end
    @test ER.Derivation._guard(half2).fault === ER.Derivation._FAST_STEP_PINNED
    @test half2 ∉ seeds
    both = @enzyme_mechanism begin
        substrates: A, B
        products: P, Q
        oligomeric_state: 2
        steps: begin
            E + A ⇌ EA                                        :: EqualAI
            EA <--> EP_res_plusA_minusP                    :: OnlyA
            EP_res_plusA_minusP ⇌ E_res_plusA_minusP + P  :: EqualAI
            E_res_plusA_minusP + B ⇌ EB_res_plusA_minusP  :: EqualAI
            EB_res_plusA_minusP ⇌ EQ                       :: OnlyA
            EQ ⇌ E + Q                                        :: EqualAI
        end
        regulatory_site(2): R::OnlyA
    end
    @test both ∈ seeds
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
    # is retained (it leaves all chemical steps :OnlyA); the chemical step is not
    # relaxed while a binding is :OnlyA, because an :OnlyA binding requires a
    # catalytically-dead inactive conformation. That skipped
    # :NonequalAI-chemistry variant is anyway rate-equivalent to this fully-dead
    # form (the inactive binds nothing, so k_I is unobservable), so no hypothesis
    # is lost.
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
    # never strand an :OnlyA binding, and that branch needs no filter.
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

@testset "the tagging move builds a three-step ping-pong on both half-reactions" begin
    # The three-step ping-pong of census pingpong seed 8fe71fded09a27ec: free E takes up
    # A and gives off P in one Theorell–Chance step, leaving the covalent form E(;
    # residual = A - P), which takes up B into EQ (fused); Q's release is the one
    # binding group. Tagging both half-reactions :OnlyA cuts the covalent form off in the
    # inactive conformation, and the guard admits it, so that set is the one base, though
    # the guard admits each half alone too. The move emits Q's binding :OnlyA on the base
    # and the two V-types that pair it with the declared regulator R; the exits keep all
    # three.
    m = @enzyme_mechanism begin
        substrates: A, B
        products: P, Q
        steps: begin
            E + A <--> E_res_plusA_minusP + P
            E + Q <--> EQ
            E_res_plusA_minusP + B <--> EQ
        end
    end
    rxn = @enzyme_reaction begin
        substrates: A[CX], B[N]
        products: P[C], Q[NX]
        allosteric_regulators: R
    end
    half1 = @enzyme_mechanism begin
        substrates: A, B ; products: P, Q ; oligomeric_state: 1
        steps: begin
            E + A <--> E_res_plusA_minusP + P :: OnlyA
            E + Q <--> EQ                      :: EqualAI
            E_res_plusA_minusP + B <--> EQ  :: EqualAI
        end
    end
    half2 = @enzyme_mechanism begin
        substrates: A, B ; products: P, Q ; oligomeric_state: 1
        steps: begin
            E + A <--> E_res_plusA_minusP + P :: EqualAI
            E + Q <--> EQ                      :: EqualAI
            E_res_plusA_minusP + B <--> EQ  :: OnlyA
        end
    end
    every = @enzyme_mechanism begin
        substrates: A, B ; products: P, Q ; oligomeric_state: 1
        steps: begin
            E + A <--> E_res_plusA_minusP + P :: OnlyA
            E + Q <--> EQ                      :: EqualAI
            E_res_plusA_minusP + B <--> EQ  :: OnlyA
        end
    end
    k_q_every = @enzyme_mechanism begin
        substrates: A, B ; products: P, Q ; oligomeric_state: 1
        steps: begin
            E + A <--> E_res_plusA_minusP + P :: OnlyA
            E + Q <--> EQ                      :: OnlyA
            E_res_plusA_minusP + B <--> EQ  :: OnlyA
        end
    end
    v_every_a = @enzyme_mechanism begin
        substrates: A, B ; products: P, Q ; oligomeric_state: 1
        steps: begin
            E + A <--> E_res_plusA_minusP + P :: OnlyA
            E + Q <--> EQ                      :: EqualAI
            E_res_plusA_minusP + B <--> EQ  :: OnlyA
        end
        regulatory_site(1): R::OnlyA
    end
    v_every_i = @enzyme_mechanism begin
        substrates: A, B ; products: P, Q ; oligomeric_state: 1
        steps: begin
            E + A <--> E_res_plusA_minusP + P :: OnlyA
            E + Q <--> EQ                      :: EqualAI
            E_res_plusA_minusP + B <--> EQ  :: OnlyA
        end
        regulatory_site(1): R::OnlyI
    end
    @test all(_testhelper_is_valid, [half1, half2])
    @test Set(
        EnzymeRates.Mechanism(m.catalytic, EnzymeRates.Mechanisms.Allostery(
            t, 1, EnzymeRates.Mechanisms.RegulatorySite[]))
        for t in EnzymeRates.Enumeration._tagging_bases(
            _testhelper_candidate(m))) == Set([every])
    kids = _testhelper_move(EnzymeRates.Enumeration._expand_to_allosteric, m, rxn)
    @test length(kids) == 3
    @test Set(kids) == Set([k_q_every, v_every_a, v_every_i])
    # Q's binding :OnlyA on the base cuts EQ off as well; both oracles agree with the
    # law on every face.
    @test all(_testhelper_admission_agreement(
        k_q_every; rng = Random.MersenneTwister(15209), corners = true))
    _testhelper_check_law_against_limit(k_q_every)
    @test Set(c for c in _testhelper_children(rxn, [m])
              if _testhelper_is_allosteric(c)) == Set(kids)
end

@testset "the tagging move builds a random-order ping-pong on both half-reactions" begin
    # A ping-pong of census pingpong (a plain parent of the census): A and Q bind free E
    # at random, P and B bind the covalent form E_res_plusA_minusP, and two chemistry
    # groups, EA to EP_res_plusA_minusP and EB_res_plusA_minusP to EQ, make the
    # two half-reactions. Tagging both chemistry groups :OnlyA (`every`) cuts the covalent
    # forms off in the inactive conformation, and the guard admits it, so that set is the
    # one base. Admission's interior accepts the first half's chemistry :OnlyA alone
    # (`half1`), but the guard does not: a face holds a rapid-equilibrium segment whose
    # lowest forms split (`_SPLIT`). The move emits two V-types with R on the base and a
    # K-type for each nonempty set of the four bindings :OnlyA on it (15). The exits keep
    # seven: the V-types and the K-types without B's binding that tag P's binding only
    # beside A's.
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
        allosteric_regulators: R
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
    every = @_testhelper_unguarded @enzyme_mechanism begin
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
    D = EnzymeRates.Derivation
    @test D._admission_interior(half1, D._built_reading(half1)) === D._WELL_POSED
    @test D._guard(half1).fault === D._SPLIT
    @test all(_testhelper_admission_agreement(every; rng = Random.MersenneTwister(15327),
                                              corners = true))
    _testhelper_check_law_against_limit(every)
    bases = EnzymeRates.Enumeration._tagging_bases(_testhelper_candidate(m))
    @test length(bases) == 1
    @test EnzymeRates.Mechanism(m.catalytic, EnzymeRates.Mechanisms.Allostery(
              only(bases), 1, EnzymeRates.Mechanisms.RegulatorySite[])) == every
    a_every = @_testhelper_unguarded @enzyme_mechanism begin
        substrates: A, B ; products: P, Q ; oligomeric_state: 1
        steps: begin
            (E + A ⇌ EA, EQ + A ⇌ EAQ)                 :: OnlyA
            (E + Q ⇌ EQ, EA + Q ⇌ EAQ)                 :: EqualAI
            E_res_plusA_minusP + B ⇌ EB_res_plusA_minusP :: EqualAI
            E_res_plusA_minusP + P ⇌ EP_res_plusA_minusP :: EqualAI
            EA <--> EP_res_plusA_minusP                   :: OnlyA
            EB_res_plusA_minusP ⇌ EQ                      :: OnlyA
        end
    end
    q_every = @_testhelper_unguarded @enzyme_mechanism begin
        substrates: A, B ; products: P, Q ; oligomeric_state: 1
        steps: begin
            (E + A ⇌ EA, EQ + A ⇌ EAQ)                 :: EqualAI
            (E + Q ⇌ EQ, EA + Q ⇌ EAQ)                 :: OnlyA
            E_res_plusA_minusP + B ⇌ EB_res_plusA_minusP :: EqualAI
            E_res_plusA_minusP + P ⇌ EP_res_plusA_minusP :: EqualAI
            EA <--> EP_res_plusA_minusP                   :: OnlyA
            EB_res_plusA_minusP ⇌ EQ                      :: OnlyA
        end
    end
    aq_every = @_testhelper_unguarded @enzyme_mechanism begin
        substrates: A, B ; products: P, Q ; oligomeric_state: 1
        steps: begin
            (E + A ⇌ EA, EQ + A ⇌ EAQ)                 :: OnlyA
            (E + Q ⇌ EQ, EA + Q ⇌ EAQ)                 :: OnlyA
            E_res_plusA_minusP + B ⇌ EB_res_plusA_minusP :: EqualAI
            E_res_plusA_minusP + P ⇌ EP_res_plusA_minusP :: EqualAI
            EA <--> EP_res_plusA_minusP                   :: OnlyA
            EB_res_plusA_minusP ⇌ EQ                      :: OnlyA
        end
    end
    b_every = @_testhelper_unguarded @enzyme_mechanism begin
        substrates: A, B ; products: P, Q ; oligomeric_state: 1
        steps: begin
            (E + A ⇌ EA, EQ + A ⇌ EAQ)                 :: EqualAI
            (E + Q ⇌ EQ, EA + Q ⇌ EAQ)                 :: EqualAI
            E_res_plusA_minusP + B ⇌ EB_res_plusA_minusP :: OnlyA
            E_res_plusA_minusP + P ⇌ EP_res_plusA_minusP :: EqualAI
            EA <--> EP_res_plusA_minusP                   :: OnlyA
            EB_res_plusA_minusP ⇌ EQ                      :: OnlyA
        end
    end
    ab_every = @_testhelper_unguarded @enzyme_mechanism begin
        substrates: A, B ; products: P, Q ; oligomeric_state: 1
        steps: begin
            (E + A ⇌ EA, EQ + A ⇌ EAQ)                 :: OnlyA
            (E + Q ⇌ EQ, EA + Q ⇌ EAQ)                 :: EqualAI
            E_res_plusA_minusP + B ⇌ EB_res_plusA_minusP :: OnlyA
            E_res_plusA_minusP + P ⇌ EP_res_plusA_minusP :: EqualAI
            EA <--> EP_res_plusA_minusP                   :: OnlyA
            EB_res_plusA_minusP ⇌ EQ                      :: OnlyA
        end
    end
    bq_every = @_testhelper_unguarded @enzyme_mechanism begin
        substrates: A, B ; products: P, Q ; oligomeric_state: 1
        steps: begin
            (E + A ⇌ EA, EQ + A ⇌ EAQ)                 :: EqualAI
            (E + Q ⇌ EQ, EA + Q ⇌ EAQ)                 :: OnlyA
            E_res_plusA_minusP + B ⇌ EB_res_plusA_minusP :: OnlyA
            E_res_plusA_minusP + P ⇌ EP_res_plusA_minusP :: EqualAI
            EA <--> EP_res_plusA_minusP                   :: OnlyA
            EB_res_plusA_minusP ⇌ EQ                      :: OnlyA
        end
    end
    abq_every = @_testhelper_unguarded @enzyme_mechanism begin
        substrates: A, B ; products: P, Q ; oligomeric_state: 1
        steps: begin
            (E + A ⇌ EA, EQ + A ⇌ EAQ)                 :: OnlyA
            (E + Q ⇌ EQ, EA + Q ⇌ EAQ)                 :: OnlyA
            E_res_plusA_minusP + B ⇌ EB_res_plusA_minusP :: OnlyA
            E_res_plusA_minusP + P ⇌ EP_res_plusA_minusP :: EqualAI
            EA <--> EP_res_plusA_minusP                   :: OnlyA
            EB_res_plusA_minusP ⇌ EQ                      :: OnlyA
        end
    end
    p_every = @_testhelper_unguarded @enzyme_mechanism begin
        substrates: A, B ; products: P, Q ; oligomeric_state: 1
        steps: begin
            (E + A ⇌ EA, EQ + A ⇌ EAQ)                 :: EqualAI
            (E + Q ⇌ EQ, EA + Q ⇌ EAQ)                 :: EqualAI
            E_res_plusA_minusP + B ⇌ EB_res_plusA_minusP :: EqualAI
            E_res_plusA_minusP + P ⇌ EP_res_plusA_minusP :: OnlyA
            EA <--> EP_res_plusA_minusP                   :: OnlyA
            EB_res_plusA_minusP ⇌ EQ                      :: OnlyA
        end
    end
    ap_every = @_testhelper_unguarded @enzyme_mechanism begin
        substrates: A, B ; products: P, Q ; oligomeric_state: 1
        steps: begin
            (E + A ⇌ EA, EQ + A ⇌ EAQ)                 :: OnlyA
            (E + Q ⇌ EQ, EA + Q ⇌ EAQ)                 :: EqualAI
            E_res_plusA_minusP + B ⇌ EB_res_plusA_minusP :: EqualAI
            E_res_plusA_minusP + P ⇌ EP_res_plusA_minusP :: OnlyA
            EA <--> EP_res_plusA_minusP                   :: OnlyA
            EB_res_plusA_minusP ⇌ EQ                      :: OnlyA
        end
    end
    pq_every = @_testhelper_unguarded @enzyme_mechanism begin
        substrates: A, B ; products: P, Q ; oligomeric_state: 1
        steps: begin
            (E + A ⇌ EA, EQ + A ⇌ EAQ)                 :: EqualAI
            (E + Q ⇌ EQ, EA + Q ⇌ EAQ)                 :: OnlyA
            E_res_plusA_minusP + B ⇌ EB_res_plusA_minusP :: EqualAI
            E_res_plusA_minusP + P ⇌ EP_res_plusA_minusP :: OnlyA
            EA <--> EP_res_plusA_minusP                   :: OnlyA
            EB_res_plusA_minusP ⇌ EQ                      :: OnlyA
        end
    end
    apq_every = @_testhelper_unguarded @enzyme_mechanism begin
        substrates: A, B ; products: P, Q ; oligomeric_state: 1
        steps: begin
            (E + A ⇌ EA, EQ + A ⇌ EAQ)                 :: OnlyA
            (E + Q ⇌ EQ, EA + Q ⇌ EAQ)                 :: OnlyA
            E_res_plusA_minusP + B ⇌ EB_res_plusA_minusP :: EqualAI
            E_res_plusA_minusP + P ⇌ EP_res_plusA_minusP :: OnlyA
            EA <--> EP_res_plusA_minusP                   :: OnlyA
            EB_res_plusA_minusP ⇌ EQ                      :: OnlyA
        end
    end
    bp_every = @_testhelper_unguarded @enzyme_mechanism begin
        substrates: A, B ; products: P, Q ; oligomeric_state: 1
        steps: begin
            (E + A ⇌ EA, EQ + A ⇌ EAQ)                 :: EqualAI
            (E + Q ⇌ EQ, EA + Q ⇌ EAQ)                 :: EqualAI
            E_res_plusA_minusP + B ⇌ EB_res_plusA_minusP :: OnlyA
            E_res_plusA_minusP + P ⇌ EP_res_plusA_minusP :: OnlyA
            EA <--> EP_res_plusA_minusP                   :: OnlyA
            EB_res_plusA_minusP ⇌ EQ                      :: OnlyA
        end
    end
    abp_every = @_testhelper_unguarded @enzyme_mechanism begin
        substrates: A, B ; products: P, Q ; oligomeric_state: 1
        steps: begin
            (E + A ⇌ EA, EQ + A ⇌ EAQ)                 :: OnlyA
            (E + Q ⇌ EQ, EA + Q ⇌ EAQ)                 :: EqualAI
            E_res_plusA_minusP + B ⇌ EB_res_plusA_minusP :: OnlyA
            E_res_plusA_minusP + P ⇌ EP_res_plusA_minusP :: OnlyA
            EA <--> EP_res_plusA_minusP                   :: OnlyA
            EB_res_plusA_minusP ⇌ EQ                      :: OnlyA
        end
    end
    bpq_every = @_testhelper_unguarded @enzyme_mechanism begin
        substrates: A, B ; products: P, Q ; oligomeric_state: 1
        steps: begin
            (E + A ⇌ EA, EQ + A ⇌ EAQ)                 :: EqualAI
            (E + Q ⇌ EQ, EA + Q ⇌ EAQ)                 :: OnlyA
            E_res_plusA_minusP + B ⇌ EB_res_plusA_minusP :: OnlyA
            E_res_plusA_minusP + P ⇌ EP_res_plusA_minusP :: OnlyA
            EA <--> EP_res_plusA_minusP                   :: OnlyA
            EB_res_plusA_minusP ⇌ EQ                      :: OnlyA
        end
    end
    abpq_every = @_testhelper_unguarded @enzyme_mechanism begin
        substrates: A, B ; products: P, Q ; oligomeric_state: 1
        steps: begin
            (E + A ⇌ EA, EQ + A ⇌ EAQ)                 :: OnlyA
            (E + Q ⇌ EQ, EA + Q ⇌ EAQ)                 :: OnlyA
            E_res_plusA_minusP + B ⇌ EB_res_plusA_minusP :: OnlyA
            E_res_plusA_minusP + P ⇌ EP_res_plusA_minusP :: OnlyA
            EA <--> EP_res_plusA_minusP                   :: OnlyA
            EB_res_plusA_minusP ⇌ EQ                      :: OnlyA
        end
    end
    v_a = @_testhelper_unguarded @enzyme_mechanism begin
        substrates: A, B ; products: P, Q ; oligomeric_state: 1
        steps: begin
            (E + A ⇌ EA, EQ + A ⇌ EAQ)                 :: EqualAI
            (E + Q ⇌ EQ, EA + Q ⇌ EAQ)                 :: EqualAI
            E_res_plusA_minusP + B ⇌ EB_res_plusA_minusP :: EqualAI
            E_res_plusA_minusP + P ⇌ EP_res_plusA_minusP :: EqualAI
            EA <--> EP_res_plusA_minusP                   :: OnlyA
            EB_res_plusA_minusP ⇌ EQ                      :: OnlyA
        end
        regulatory_site(1): R::OnlyA
    end
    v_i = @_testhelper_unguarded @enzyme_mechanism begin
        substrates: A, B ; products: P, Q ; oligomeric_state: 1
        steps: begin
            (E + A ⇌ EA, EQ + A ⇌ EAQ)                 :: EqualAI
            (E + Q ⇌ EQ, EA + Q ⇌ EAQ)                 :: EqualAI
            E_res_plusA_minusP + B ⇌ EB_res_plusA_minusP :: EqualAI
            E_res_plusA_minusP + P ⇌ EP_res_plusA_minusP :: EqualAI
            EA <--> EP_res_plusA_minusP                   :: OnlyA
            EB_res_plusA_minusP ⇌ EQ                      :: OnlyA
        end
        regulatory_site(1): R::OnlyI
    end
    kids = _testhelper_move(EnzymeRates.Enumeration._expand_to_allosteric, m, rxn)
    @test length(kids) == 17
    @test Set(kids) == Set([a_every, q_every, aq_every, b_every, ab_every, bq_every,
                            abq_every, p_every, ap_every, pq_every, apq_every, bp_every,
                            abp_every, bpq_every, abpq_every, v_a, v_i])
    kept = [a_every, q_every, aq_every, ap_every, apq_every, v_a, v_i]
    @test all(_testhelper_is_valid, kept)
    @test !any(_testhelper_is_valid, setdiff(kids, kept))
    @test Set(c for c in _testhelper_children(rxn, [m])
              if _testhelper_is_allosteric(c)) == Set(kept)
end

@testset "the tagging search extends a tag set whose child the guard rejects" begin
    # A ter-bi ping-pong of the census ter-bi reaction: B and D bind free E in either
    # order, the second fused into the product-bound covalent form E(Q; residual = B + D -
    # Q) (B's and D's groups, each chemistry); Q's release from it leaves E(; residual =
    # B + D - Q), which takes up A into EP, fused (A's group), and Q's group holds Q's
    # abortive binding at EP beside its release. P's binding at E is the one binding
    # group. Every chemistry group :OnlyA fails the guard: at zero A, B and D the fused A
    # binding, vanishing one way into the covalent form, fills it faster than it empties
    # (`_FILLED`). So the bases are the smallest sets of chemistry groups that stop the
    # inactive conformation's catalysis and pass, and a set that fails leaves its
    # supersets in the search. Q's group alone keeps the two sides of its
    # rapid-equilibrium release joined in the inactive conformation (`_FAST_STEP_PINNED`),
    # and the search extends it to B's and Q's groups and to D's and Q's, which pass. B's
    # group alone passes the guard but keeps EB, which only the binding it excludes
    # enters (`_kept_forms_keep_a_binding`), so it fails too and extends to B's and D's
    # groups; D's alone likewise. A's group alone keeps the sides of its fused binding
    # joined, and each pair that holds it fails the guard too. Three bases: the move emits
    # a K-type with P's binding :OnlyA on each and on the every-chemistry set, and two
    # V-types with R on each base. The exits keep all ten.
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
        allosteric_regulators: R
    end
    base_bd = @enzyme_mechanism begin
        substrates: A, B, D ; products: P, Q ; oligomeric_state: 1
        steps: begin
            (E + B <--> EB, ED + B <--> EQ_res_plusB_plusD_minusQ)      :: OnlyA
            (E + D <--> ED, EB + D <--> EQ_res_plusB_plusD_minusQ)      :: OnlyA
            E + P <--> EP                                                  :: EqualAI
            (E_res_plusB_plusD_minusQ + A ⇌ EP,
             EQ_res_plusB_plusD_minusQ + A ⇌ EAQ_res_plusB_plusD_minusQ) :: EqualAI
            (E_res_plusB_plusD_minusQ + Q ⇌ EQ_res_plusB_plusD_minusQ,
             EP + Q ⇌ EAQ_res_plusB_plusD_minusQ)                       :: EqualAI
        end
    end
    base_bq = @enzyme_mechanism begin
        substrates: A, B, D ; products: P, Q ; oligomeric_state: 1
        steps: begin
            (E + B <--> EB, ED + B <--> EQ_res_plusB_plusD_minusQ)      :: OnlyA
            (E + D <--> ED, EB + D <--> EQ_res_plusB_plusD_minusQ)      :: EqualAI
            E + P <--> EP                                                  :: EqualAI
            (E_res_plusB_plusD_minusQ + A ⇌ EP,
             EQ_res_plusB_plusD_minusQ + A ⇌ EAQ_res_plusB_plusD_minusQ) :: EqualAI
            (E_res_plusB_plusD_minusQ + Q ⇌ EQ_res_plusB_plusD_minusQ,
             EP + Q ⇌ EAQ_res_plusB_plusD_minusQ)                       :: OnlyA
        end
    end
    base_dq = @enzyme_mechanism begin
        substrates: A, B, D ; products: P, Q ; oligomeric_state: 1
        steps: begin
            (E + B <--> EB, ED + B <--> EQ_res_plusB_plusD_minusQ)      :: EqualAI
            (E + D <--> ED, EB + D <--> EQ_res_plusB_plusD_minusQ)      :: OnlyA
            E + P <--> EP                                                  :: EqualAI
            (E_res_plusB_plusD_minusQ + A ⇌ EP,
             EQ_res_plusB_plusD_minusQ + A ⇌ EAQ_res_plusB_plusD_minusQ) :: EqualAI
            (E_res_plusB_plusD_minusQ + Q ⇌ EQ_res_plusB_plusD_minusQ,
             EP + Q ⇌ EAQ_res_plusB_plusD_minusQ)                       :: OnlyA
        end
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
    p_on_bd = @enzyme_mechanism begin
        substrates: A, B, D ; products: P, Q ; oligomeric_state: 1
        steps: begin
            (E + B <--> EB, ED + B <--> EQ_res_plusB_plusD_minusQ)      :: OnlyA
            (E + D <--> ED, EB + D <--> EQ_res_plusB_plusD_minusQ)      :: OnlyA
            E + P <--> EP                                                  :: OnlyA
            (E_res_plusB_plusD_minusQ + A ⇌ EP,
             EQ_res_plusB_plusD_minusQ + A ⇌ EAQ_res_plusB_plusD_minusQ) :: EqualAI
            (E_res_plusB_plusD_minusQ + Q ⇌ EQ_res_plusB_plusD_minusQ,
             EP + Q ⇌ EAQ_res_plusB_plusD_minusQ)                       :: EqualAI
        end
    end
    p_on_bq = @enzyme_mechanism begin
        substrates: A, B, D ; products: P, Q ; oligomeric_state: 1
        steps: begin
            (E + B <--> EB, ED + B <--> EQ_res_plusB_plusD_minusQ)      :: OnlyA
            (E + D <--> ED, EB + D <--> EQ_res_plusB_plusD_minusQ)      :: EqualAI
            E + P <--> EP                                                  :: OnlyA
            (E_res_plusB_plusD_minusQ + A ⇌ EP,
             EQ_res_plusB_plusD_minusQ + A ⇌ EAQ_res_plusB_plusD_minusQ) :: EqualAI
            (E_res_plusB_plusD_minusQ + Q ⇌ EQ_res_plusB_plusD_minusQ,
             EP + Q ⇌ EAQ_res_plusB_plusD_minusQ)                       :: OnlyA
        end
    end
    p_on_dq = @enzyme_mechanism begin
        substrates: A, B, D ; products: P, Q ; oligomeric_state: 1
        steps: begin
            (E + B <--> EB, ED + B <--> EQ_res_plusB_plusD_minusQ)      :: EqualAI
            (E + D <--> ED, EB + D <--> EQ_res_plusB_plusD_minusQ)      :: OnlyA
            E + P <--> EP                                                  :: OnlyA
            (E_res_plusB_plusD_minusQ + A ⇌ EP,
             EQ_res_plusB_plusD_minusQ + A ⇌ EAQ_res_plusB_plusD_minusQ) :: EqualAI
            (E_res_plusB_plusD_minusQ + Q ⇌ EQ_res_plusB_plusD_minusQ,
             EP + Q ⇌ EAQ_res_plusB_plusD_minusQ)                       :: OnlyA
        end
    end
    p_on_every = @enzyme_mechanism begin
        substrates: A, B, D ; products: P, Q ; oligomeric_state: 1
        steps: begin
            (E + B <--> EB, ED + B <--> EQ_res_plusB_plusD_minusQ)      :: OnlyA
            (E + D <--> ED, EB + D <--> EQ_res_plusB_plusD_minusQ)      :: OnlyA
            E + P <--> EP                                                  :: OnlyA
            (E_res_plusB_plusD_minusQ + A ⇌ EP,
             EQ_res_plusB_plusD_minusQ + A ⇌ EAQ_res_plusB_plusD_minusQ) :: OnlyA
            (E_res_plusB_plusD_minusQ + Q ⇌ EQ_res_plusB_plusD_minusQ,
             EP + Q ⇌ EAQ_res_plusB_plusD_minusQ)                       :: OnlyA
        end
    end
    r_on_bd_a = @enzyme_mechanism begin
        substrates: A, B, D ; products: P, Q ; oligomeric_state: 1
        steps: begin
            (E + B <--> EB, ED + B <--> EQ_res_plusB_plusD_minusQ)      :: OnlyA
            (E + D <--> ED, EB + D <--> EQ_res_plusB_plusD_minusQ)      :: OnlyA
            E + P <--> EP                                                  :: EqualAI
            (E_res_plusB_plusD_minusQ + A ⇌ EP,
             EQ_res_plusB_plusD_minusQ + A ⇌ EAQ_res_plusB_plusD_minusQ) :: EqualAI
            (E_res_plusB_plusD_minusQ + Q ⇌ EQ_res_plusB_plusD_minusQ,
             EP + Q ⇌ EAQ_res_plusB_plusD_minusQ)                       :: EqualAI
        end
        regulatory_site(1): R::OnlyA
    end
    r_on_bd_i = @enzyme_mechanism begin
        substrates: A, B, D ; products: P, Q ; oligomeric_state: 1
        steps: begin
            (E + B <--> EB, ED + B <--> EQ_res_plusB_plusD_minusQ)      :: OnlyA
            (E + D <--> ED, EB + D <--> EQ_res_plusB_plusD_minusQ)      :: OnlyA
            E + P <--> EP                                                  :: EqualAI
            (E_res_plusB_plusD_minusQ + A ⇌ EP,
             EQ_res_plusB_plusD_minusQ + A ⇌ EAQ_res_plusB_plusD_minusQ) :: EqualAI
            (E_res_plusB_plusD_minusQ + Q ⇌ EQ_res_plusB_plusD_minusQ,
             EP + Q ⇌ EAQ_res_plusB_plusD_minusQ)                       :: EqualAI
        end
        regulatory_site(1): R::OnlyI
    end
    r_on_bq_a = @enzyme_mechanism begin
        substrates: A, B, D ; products: P, Q ; oligomeric_state: 1
        steps: begin
            (E + B <--> EB, ED + B <--> EQ_res_plusB_plusD_minusQ)      :: OnlyA
            (E + D <--> ED, EB + D <--> EQ_res_plusB_plusD_minusQ)      :: EqualAI
            E + P <--> EP                                                  :: EqualAI
            (E_res_plusB_plusD_minusQ + A ⇌ EP,
             EQ_res_plusB_plusD_minusQ + A ⇌ EAQ_res_plusB_plusD_minusQ) :: EqualAI
            (E_res_plusB_plusD_minusQ + Q ⇌ EQ_res_plusB_plusD_minusQ,
             EP + Q ⇌ EAQ_res_plusB_plusD_minusQ)                       :: OnlyA
        end
        regulatory_site(1): R::OnlyA
    end
    r_on_bq_i = @enzyme_mechanism begin
        substrates: A, B, D ; products: P, Q ; oligomeric_state: 1
        steps: begin
            (E + B <--> EB, ED + B <--> EQ_res_plusB_plusD_minusQ)      :: OnlyA
            (E + D <--> ED, EB + D <--> EQ_res_plusB_plusD_minusQ)      :: EqualAI
            E + P <--> EP                                                  :: EqualAI
            (E_res_plusB_plusD_minusQ + A ⇌ EP,
             EQ_res_plusB_plusD_minusQ + A ⇌ EAQ_res_plusB_plusD_minusQ) :: EqualAI
            (E_res_plusB_plusD_minusQ + Q ⇌ EQ_res_plusB_plusD_minusQ,
             EP + Q ⇌ EAQ_res_plusB_plusD_minusQ)                       :: OnlyA
        end
        regulatory_site(1): R::OnlyI
    end
    r_on_dq_a = @enzyme_mechanism begin
        substrates: A, B, D ; products: P, Q ; oligomeric_state: 1
        steps: begin
            (E + B <--> EB, ED + B <--> EQ_res_plusB_plusD_minusQ)      :: EqualAI
            (E + D <--> ED, EB + D <--> EQ_res_plusB_plusD_minusQ)      :: OnlyA
            E + P <--> EP                                                  :: EqualAI
            (E_res_plusB_plusD_minusQ + A ⇌ EP,
             EQ_res_plusB_plusD_minusQ + A ⇌ EAQ_res_plusB_plusD_minusQ) :: EqualAI
            (E_res_plusB_plusD_minusQ + Q ⇌ EQ_res_plusB_plusD_minusQ,
             EP + Q ⇌ EAQ_res_plusB_plusD_minusQ)                       :: OnlyA
        end
        regulatory_site(1): R::OnlyA
    end
    r_on_dq_i = @enzyme_mechanism begin
        substrates: A, B, D ; products: P, Q ; oligomeric_state: 1
        steps: begin
            (E + B <--> EB, ED + B <--> EQ_res_plusB_plusD_minusQ)      :: EqualAI
            (E + D <--> ED, EB + D <--> EQ_res_plusB_plusD_minusQ)      :: OnlyA
            E + P <--> EP                                                  :: EqualAI
            (E_res_plusB_plusD_minusQ + A ⇌ EP,
             EQ_res_plusB_plusD_minusQ + A ⇌ EAQ_res_plusB_plusD_minusQ) :: EqualAI
            (E_res_plusB_plusD_minusQ + Q ⇌ EQ_res_plusB_plusD_minusQ,
             EP + Q ⇌ EAQ_res_plusB_plusD_minusQ)                       :: OnlyA
        end
        regulatory_site(1): R::OnlyI
    end

    D = EnzymeRates.Derivation
    # The guard names the group's first step in step order.
    v = D._guard(q_alone)
    @test (v.fault, EnzymeRates.Mechanisms.step_text(m.catalytic, v.edge)) ==
          (D._FAST_STEP_PINNED, "EP + Q ⇌ EAQ_res_plusB_plusD_minusQ")
    @test _testhelper_is_valid(b_alone)
    @test !EnzymeRates.Enumeration._kept_forms_keep_a_binding(b_alone)
    @test _testhelper_is_valid(d_alone)
    @test !EnzymeRates.Enumeration._kept_forms_keep_a_binding(d_alone)
    v = D._guard(a_alone)
    @test (v.fault, EnzymeRates.Mechanisms.step_text(m.catalytic, v.edge)) ==
          (D._FAST_STEP_PINNED,
           "EQ_res_plusB_plusD_minusQ + A ⇌ EAQ_res_plusB_plusD_minusQ")
    # Each single group stops the inactive conformation's catalysis, so each is a
    # candidate the search tries before the pairs.
    @test all(c -> c.catalytic == m.catalytic &&
                   EnzymeRates.Enumeration._stops_inactive_turnover(
                       _testhelper_candidate(m), c.allosteric.allosteric_states),
              [q_alone, b_alone, d_alone, a_alone])
    bases = EnzymeRates.Enumeration._tagging_bases(_testhelper_candidate(m))
    @test length(bases) == 3
    @test Set(EnzymeRates.Mechanism(m.catalytic, EnzymeRates.Mechanisms.Allostery(t, 1,
                                              EnzymeRates.Mechanisms.RegulatorySite[]))
              for t in bases) == Set([base_bd, base_bq, base_dq])
    kids = _testhelper_move(EnzymeRates.Enumeration._expand_to_allosteric, m, rxn)
    @test length(kids) == 10
    @test Set(kids) == Set([p_on_bd, p_on_bq, p_on_dq, p_on_every, r_on_bd_a, r_on_bd_i,
                            r_on_bq_a, r_on_bq_i, r_on_dq_a, r_on_dq_i])
    @test Set(c for c in _testhelper_children(rxn, [m])
              if _testhelper_is_allosteric(c)) == Set(kids)
end

@testset "change_allo_state relaxes a ping-pong K-type's binding, not its chemistry" begin
    # A hand-built K-type of the three-step ping-pong (census pingpong seed
    # 8fe71fded09a27ec): the second half-reaction and Q's binding :OnlyA, the first half
    # running in both conformations. The tagging move builds this seed's K-type on both
    # half-reactions instead (testset "the tagging move builds a three-step ping-pong on
    # both half-reactions"). Q's binding relaxes, and admission accepts it; the chemistry
    # does not relax while a binding is :OnlyA. One child.
    k_type = @_testhelper_unguarded @enzyme_mechanism begin
        substrates: A, B ; products: P, Q ; oligomeric_state: 1
        steps: begin
            E + A <--> E_res_plusA_minusP + P :: EqualAI
            E + Q <--> EQ                      :: OnlyA
            E_res_plusA_minusP + B <--> EQ  :: OnlyA
        end
    end
    q_relaxed = @_testhelper_unguarded @enzyme_mechanism begin
        substrates: A, B ; products: P, Q ; oligomeric_state: 1
        steps: begin
            E + A <--> E_res_plusA_minusP + P :: EqualAI
            E + Q <--> EQ                      :: NonequalAI
            E_res_plusA_minusP + B <--> EQ  :: OnlyA
        end
    end
    kids = _testhelper_move(EnzymeRates.Enumeration._expand_change_allo_state, k_type)
    @test length(kids) == 1
    @test Set(kids) == Set([q_relaxed])
    @test _testhelper_is_valid(k_type) &&
          _testhelper_is_valid(q_relaxed)
end

@testset "change_allo_state relaxes a ping-pong V-type, leaving class 3 to the exit" begin
    # The tagging move's V-type on the six-step bi-bi ping-pong (testset "Mechanism —
    # Bi-bi ping-pong: …" above), at one subunit: both half-reactions :OnlyA, every binding
    # :EqualAI, R at an :OnlyA site. Each binding relaxes alone, the two half-reactions
    # relax together, and the site ligand relaxes; admission accepts every catalytic
    # relaxation. With both halves :NonequalAI the conformations run different
    # chemistry on one rapid-equilibrium segment, and with every substrate zero or B zero
    # the covalent forms hold the enzyme in both, joined only through free E: a split whose
    # limit depends on the speeds. Class 3 rejects it, so candidate_mechanisms drops that
    # child.
    rxn = @enzyme_reaction begin
        substrates: A[CX], B[N]
        products: P[C], Q[NX]
        allosteric_regulators: R
    end
    v_type = @_testhelper_unguarded @enzyme_mechanism begin
        substrates: A, B ; products: P, Q ; oligomeric_state: 1
        steps: begin
            E + A ⇌ EA                                       :: EqualAI
            E_res_plusA_minusP + B ⇌ EB_res_plusA_minusP :: EqualAI
            E + Q ⇌ EQ                                       :: EqualAI
            E_res_plusA_minusP + P ⇌ EP_res_plusA_minusP :: EqualAI
            EA <--> EP_res_plusA_minusP                   :: OnlyA
            EB_res_plusA_minusP ⇌ EQ                      :: OnlyA
        end
        regulatory_site(1): R::OnlyA
    end
    a_relaxed = @_testhelper_unguarded @enzyme_mechanism begin
        substrates: A, B ; products: P, Q ; oligomeric_state: 1
        steps: begin
            E + A ⇌ EA                                       :: NonequalAI
            E_res_plusA_minusP + B ⇌ EB_res_plusA_minusP :: EqualAI
            E + Q ⇌ EQ                                       :: EqualAI
            E_res_plusA_minusP + P ⇌ EP_res_plusA_minusP :: EqualAI
            EA <--> EP_res_plusA_minusP                   :: OnlyA
            EB_res_plusA_minusP ⇌ EQ                      :: OnlyA
        end
        regulatory_site(1): R::OnlyA
    end
    b_relaxed = @_testhelper_unguarded @enzyme_mechanism begin
        substrates: A, B ; products: P, Q ; oligomeric_state: 1
        steps: begin
            E + A ⇌ EA                                       :: EqualAI
            E_res_plusA_minusP + B ⇌ EB_res_plusA_minusP :: NonequalAI
            E + Q ⇌ EQ                                       :: EqualAI
            E_res_plusA_minusP + P ⇌ EP_res_plusA_minusP :: EqualAI
            EA <--> EP_res_plusA_minusP                   :: OnlyA
            EB_res_plusA_minusP ⇌ EQ                      :: OnlyA
        end
        regulatory_site(1): R::OnlyA
    end
    q_relaxed = @_testhelper_unguarded @enzyme_mechanism begin
        substrates: A, B ; products: P, Q ; oligomeric_state: 1
        steps: begin
            E + A ⇌ EA                                       :: EqualAI
            E_res_plusA_minusP + B ⇌ EB_res_plusA_minusP :: EqualAI
            E + Q ⇌ EQ                                       :: NonequalAI
            E_res_plusA_minusP + P ⇌ EP_res_plusA_minusP :: EqualAI
            EA <--> EP_res_plusA_minusP                   :: OnlyA
            EB_res_plusA_minusP ⇌ EQ                      :: OnlyA
        end
        regulatory_site(1): R::OnlyA
    end
    p_relaxed = @_testhelper_unguarded @enzyme_mechanism begin
        substrates: A, B ; products: P, Q ; oligomeric_state: 1
        steps: begin
            E + A ⇌ EA                                       :: EqualAI
            E_res_plusA_minusP + B ⇌ EB_res_plusA_minusP :: EqualAI
            E + Q ⇌ EQ                                       :: EqualAI
            E_res_plusA_minusP + P ⇌ EP_res_plusA_minusP :: NonequalAI
            EA <--> EP_res_plusA_minusP                   :: OnlyA
            EB_res_plusA_minusP ⇌ EQ                      :: OnlyA
        end
        regulatory_site(1): R::OnlyA
    end
    chemistry_relaxed = @_testhelper_unguarded @enzyme_mechanism begin
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
    site_relaxed = @_testhelper_unguarded @enzyme_mechanism begin
        substrates: A, B ; products: P, Q ; oligomeric_state: 1
        steps: begin
            E + A ⇌ EA                                       :: EqualAI
            E_res_plusA_minusP + B ⇌ EB_res_plusA_minusP :: EqualAI
            E + Q ⇌ EQ                                       :: EqualAI
            E_res_plusA_minusP + P ⇌ EP_res_plusA_minusP :: EqualAI
            EA <--> EP_res_plusA_minusP                   :: OnlyA
            EB_res_plusA_minusP ⇌ EQ                      :: OnlyA
        end
        regulatory_site(1): R::NonequalAI
    end
    kids = _testhelper_move(EnzymeRates.Enumeration._expand_change_allo_state, v_type)
    @test length(kids) == 6
    @test Set(kids) == Set([a_relaxed, b_relaxed, q_relaxed, p_relaxed,
                            chemistry_relaxed, site_relaxed])
    v = EnzymeRates.Derivation._guard(chemistry_relaxed)
    chem_graph = chemistry_relaxed.catalytic
    names = EnzymeRates.Derivation._face_names(
        chem_graph, EnzymeRates.Derivation._faces(chem_graph)[v.face])
    @test (v.fault, v.limit, names) ==
          (EnzymeRates.Derivation._SPLIT, EnzymeRates.Derivation._LIMIT_SPEED, [:A, :B])
    kept = [a_relaxed, b_relaxed, q_relaxed, p_relaxed, site_relaxed]
    @test all(_testhelper_is_valid, kept)
    exits = _testhelper_children(rxn, [v_type])
    @test chemistry_relaxed ∉ exits && all(in(exits), kept)
    # Every kept child turns the reaction over, at zero products too: its inactive
    # conformation runs no chemistry, so no form holds the enzyme there.
    for k in kept
        fp = parameters(k)
        prm = NamedTuple{fp}(Tuple(1.3 for _ in fp))
        rate = EnzymeRates.compile_rate_equation(k)
        @test rate((A = 1.3, B = 1.3, P = 0.0, Q = 0.0, R = 1.0), prm, 3.0) > 0
        @test rate((A = 1.3, B = 1.3, P = 0.1, Q = 0.1, R = 1.0), prm, 3.0) > 0
    end
end

@testset "change_allo_state relaxes the bindings of a hand-built partial parent" begin
    # A partial parent: every binding :OnlyA beside a live chemistry step. The
    # enumeration builds such a parent only as a K-type on a ping-pong's half-reaction
    # base; this one is written by hand. Each binding relaxation leaves three :OnlyA
    # bindings whose arrows still point away from E and close no cycle, so admission
    # accepts all four. The chemistry relaxation is not tried while a
    # binding is :OnlyA.
    partial = @enzyme_mechanism begin
        substrates: A, B ; products: P, Q ; oligomeric_state: 1
        steps: begin
            E + A ⇌ EA                :: OnlyA
            EA + B ⇌ EAB          :: OnlyA
            EAB <--> EPQ        :: EqualAI
            EQ + P ⇌ EPQ          :: OnlyA
            E + Q ⇌ EQ                :: OnlyA
        end
    end
    a_relaxed = @enzyme_mechanism begin
        substrates: A, B ; products: P, Q ; oligomeric_state: 1
        steps: begin
            E + A ⇌ EA                :: NonequalAI
            EA + B ⇌ EAB          :: OnlyA
            EAB <--> EPQ        :: EqualAI
            EQ + P ⇌ EPQ          :: OnlyA
            E + Q ⇌ EQ                :: OnlyA
        end
    end
    b_relaxed = @enzyme_mechanism begin
        substrates: A, B ; products: P, Q ; oligomeric_state: 1
        steps: begin
            E + A ⇌ EA                :: OnlyA
            EA + B ⇌ EAB          :: NonequalAI
            EAB <--> EPQ        :: EqualAI
            EQ + P ⇌ EPQ          :: OnlyA
            E + Q ⇌ EQ                :: OnlyA
        end
    end
    p_relaxed = @enzyme_mechanism begin
        substrates: A, B ; products: P, Q ; oligomeric_state: 1
        steps: begin
            E + A ⇌ EA                :: OnlyA
            EA + B ⇌ EAB          :: OnlyA
            EAB <--> EPQ        :: EqualAI
            EQ + P ⇌ EPQ          :: NonequalAI
            E + Q ⇌ EQ                :: OnlyA
        end
    end
    q_relaxed = @enzyme_mechanism begin
        substrates: A, B ; products: P, Q ; oligomeric_state: 1
        steps: begin
            E + A ⇌ EA                :: OnlyA
            EA + B ⇌ EAB          :: OnlyA
            EAB <--> EPQ        :: EqualAI
            EQ + P ⇌ EPQ          :: OnlyA
            E + Q ⇌ EQ                :: NonequalAI
        end
    end
    kids = _testhelper_move(EnzymeRates.Enumeration._expand_change_allo_state, partial)
    @test length(kids) == 4
    @test Set(kids) == Set([a_relaxed, b_relaxed, p_relaxed, q_relaxed])
    # candidate_mechanisms passes them on.
    @test issubset(kids,
                   _testhelper_children(bi_bi_rxn, [partial]))
end

@testset "enumeration emits only children the guard admits" begin
    # A plain ping-pong (census 3fdbe6e1dc001bb7) whose A binding is fused into the
    # product-bound covalent form. Its one rapid-equilibrium group flips to steady state:
    # one plain child. Tagging both half-reactions' fused bindings :OnlyA cuts the
    # covalent forms off in the inactive conformation, and the guard admits it, so that
    # set is the one base: the move emits three K-types and two V-types on it (`tagged`),
    # the guard admits all five, and the exit keeps them and the flip. A second census
    # ping-pong, whose A binding is fused at steady state and whose B binding is fused at
    # rapid equilibrium, shows the exit's drop. Its one base tags A's fused binding alone.
    # Of the six K-types on the base and on the every-chemistry set the guard rejects
    # one: P's binding :OnlyA on the every-chemistry set, where at zero A and B the fused
    # B binding fills the covalent form E_res_plusA_minusP faster than it empties
    # (`_FILLED`). The exit keeps the other seven children and the two flips.
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
    # Untyped: inside `T[...]` the parser reads `begin` as an index keyword.
    tagged = [
        @_testhelper_unguarded(@enzyme_mechanism begin
            substrates: A, B
            products: P, Q
            oligomeric_state: 2
            steps: begin
                E + A ⇌ EP_res_plusA_minusP                        :: OnlyA
                E_res_plusA_minusP + P <--> EP_res_plusA_minusP :: EqualAI
                E_res_plusA_minusP + B <--> EQ                   :: OnlyA
                E + Q <--> EQ                                       :: OnlyA
            end
        end),
        @_testhelper_unguarded(@enzyme_mechanism begin
            substrates: A, B
            products: P, Q
            oligomeric_state: 2
            steps: begin
                E + A ⇌ EP_res_plusA_minusP                        :: OnlyA
                E_res_plusA_minusP + P <--> EP_res_plusA_minusP :: OnlyA
                E_res_plusA_minusP + B <--> EQ                   :: OnlyA
                E + Q <--> EQ                                       :: EqualAI
            end
        end),
        @_testhelper_unguarded(@enzyme_mechanism begin
            substrates: A, B
            products: P, Q
            oligomeric_state: 2
            steps: begin
                E + A ⇌ EP_res_plusA_minusP                        :: OnlyA
                E_res_plusA_minusP + P <--> EP_res_plusA_minusP :: OnlyA
                E_res_plusA_minusP + B <--> EQ                   :: OnlyA
                E + Q <--> EQ                                       :: OnlyA
            end
        end),
        @_testhelper_unguarded(@enzyme_mechanism begin
            substrates: A, B
            products: P, Q
            oligomeric_state: 2
            steps: begin
                E + A ⇌ EP_res_plusA_minusP                        :: OnlyA
                E_res_plusA_minusP + P <--> EP_res_plusA_minusP :: EqualAI
                E_res_plusA_minusP + B <--> EQ                   :: OnlyA
                E + Q <--> EQ                                       :: EqualAI
            end
            regulatory_site(2): R::OnlyA
        end),
        @_testhelper_unguarded(@enzyme_mechanism begin
            substrates: A, B
            products: P, Q
            oligomeric_state: 2
            steps: begin
                E + A ⇌ EP_res_plusA_minusP                        :: OnlyA
                E_res_plusA_minusP + P <--> EP_res_plusA_minusP :: EqualAI
                E_res_plusA_minusP + B <--> EQ                   :: OnlyA
                E + Q <--> EQ                                       :: EqualAI
            end
            regulatory_site(2): R::OnlyI
        end)]
    moved = _testhelper_move(EnzymeRates.Enumeration._expand_to_allosteric, parent, rxn)
    @test length(moved) == 5
    @test Set(moved) == Set(tagged)
    @test all(_testhelper_is_valid, tagged)
    @test all(t -> all(_testhelper_admission_agreement(
                           t; rng = Random.MersenneTwister(81), corners = true)), tagged)
    @test _testhelper_move(EnzymeRates.Enumeration._expand_re_to_ss, parent) == [flip]
    kids = _testhelper_children(rxn, [parent])
    @test length(kids) == 6
    @test Set(kids) == Set([flip; tagged])
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
    # Untyped: inside `T[...]` the parser reads `begin` as an index keyword.
    second_admitted = [
        @_testhelper_unguarded(@enzyme_mechanism begin
            substrates: A, B ; products: P, Q ; oligomeric_state: 2
            steps: begin
                E + A <--> EP_res_plusA_minusP                     :: OnlyA
                E_res_plusA_minusP + P ⇌ EP_res_plusA_minusP    :: EqualAI
                E_res_plusA_minusP + B ⇌ EQ                      :: EqualAI
                E + Q <--> EQ                                       :: OnlyA
            end
        end),
        @_testhelper_unguarded(@enzyme_mechanism begin
            substrates: A, B ; products: P, Q ; oligomeric_state: 2
            steps: begin
                E + A <--> EP_res_plusA_minusP                     :: OnlyA
                E_res_plusA_minusP + P ⇌ EP_res_plusA_minusP    :: OnlyA
                E_res_plusA_minusP + B ⇌ EQ                      :: EqualAI
                E + Q <--> EQ                                       :: EqualAI
            end
        end),
        @_testhelper_unguarded(@enzyme_mechanism begin
            substrates: A, B ; products: P, Q ; oligomeric_state: 2
            steps: begin
                E + A <--> EP_res_plusA_minusP                     :: OnlyA
                E_res_plusA_minusP + P ⇌ EP_res_plusA_minusP    :: OnlyA
                E_res_plusA_minusP + B ⇌ EQ                      :: EqualAI
                E + Q <--> EQ                                       :: OnlyA
            end
        end),
        @_testhelper_unguarded(@enzyme_mechanism begin
            substrates: A, B ; products: P, Q ; oligomeric_state: 2
            steps: begin
                E + A <--> EP_res_plusA_minusP                     :: OnlyA
                E_res_plusA_minusP + P ⇌ EP_res_plusA_minusP    :: EqualAI
                E_res_plusA_minusP + B ⇌ EQ                      :: OnlyA
                E + Q <--> EQ                                       :: OnlyA
            end
        end),
        @_testhelper_unguarded(@enzyme_mechanism begin
            substrates: A, B ; products: P, Q ; oligomeric_state: 2
            steps: begin
                E + A <--> EP_res_plusA_minusP                     :: OnlyA
                E_res_plusA_minusP + P ⇌ EP_res_plusA_minusP    :: OnlyA
                E_res_plusA_minusP + B ⇌ EQ                      :: OnlyA
                E + Q <--> EQ                                       :: OnlyA
            end
        end),
        @_testhelper_unguarded(@enzyme_mechanism begin
            substrates: A, B ; products: P, Q ; oligomeric_state: 2
            steps: begin
                E + A <--> EP_res_plusA_minusP                     :: OnlyA
                E_res_plusA_minusP + P ⇌ EP_res_plusA_minusP    :: EqualAI
                E_res_plusA_minusP + B ⇌ EQ                      :: EqualAI
                E + Q <--> EQ                                       :: EqualAI
            end
            regulatory_site(2): R::OnlyA
        end),
        @_testhelper_unguarded(@enzyme_mechanism begin
            substrates: A, B ; products: P, Q ; oligomeric_state: 2
            steps: begin
                E + A <--> EP_res_plusA_minusP                     :: OnlyA
                E_res_plusA_minusP + P ⇌ EP_res_plusA_minusP    :: EqualAI
                E_res_plusA_minusP + B ⇌ EQ                      :: EqualAI
                E + Q <--> EQ                                       :: EqualAI
            end
            regulatory_site(2): R::OnlyI
        end)]
    second_rejected = @_testhelper_unguarded @enzyme_mechanism begin
        substrates: A, B ; products: P, Q ; oligomeric_state: 2
        steps: begin
            E + A <--> EP_res_plusA_minusP                     :: OnlyA
            E_res_plusA_minusP + P ⇌ EP_res_plusA_minusP    :: OnlyA
            E_res_plusA_minusP + B ⇌ EQ                      :: OnlyA
            E + Q <--> EQ                                       :: EqualAI
        end
    end
    second_moved = _testhelper_move(EnzymeRates.Enumeration._expand_to_allosteric, second,
                                    rxn)
    @test length(second_moved) == 8
    @test Set(second_moved) == Set([second_admitted; second_rejected])
    @test all(_testhelper_is_valid, second_admitted)
    v = EnzymeRates.Derivation._guard(second_rejected)
    @test (v.fault, EnzymeRates.Mechanisms.step_text(second.catalytic, v.edge)) ==
          (EnzymeRates.Derivation._FILLED, "E_res_plusA_minusP + B ⇌ EQ")
    @test all(_testhelper_admission_agreement(second_rejected;
                                              rng = Random.MersenneTwister(82),
                                              corners = true))
    second_kids = _testhelper_children(rxn, [second])
    @test length(second_kids) == 9
    @test Set(second_kids) == Set([second_flips; second_admitted])
end

@testset "the enumerator does not build a kept form that no running binding enters" begin
    En = EnzymeRates.Enumeration
    # S binding steady-state :OnlyA, the rest :EqualAI: the inactive conformation keeps
    # its two sides joined and stops it, while ES stays in its law, reached from EP by
    # the chemistry, and the guard admits the law. No binding that runs there enters
    # ES, so `_kept_forms_keep_a_binding` is false and the enumerator does not build it.
    shut_door = @enzyme_mechanism begin
        substrates: S
        products: P
        steps: begin
            E + S <--> ES     :: OnlyA
            ES <--> EP      :: EqualAI
            E + P ⇌ EP        :: EqualAI
        end
    end
    @test EnzymeRates.Derivation.check_valid_rate_equation(Bool, shut_door)
    @test !En._kept_forms_keep_a_binding(shut_door)
    @test !En._enumerator_keeps(shut_door)
    # The chemistry :OnlyA too: ES is cut off.
    k_type = @enzyme_mechanism begin
        substrates: S
        products: P
        steps: begin
            E + S <--> ES     :: OnlyA
            ES <--> EP      :: OnlyA
            E + P ⇌ EP        :: EqualAI
        end
    end
    @test En._kept_forms_keep_a_binding(k_type) && En._enumerator_keeps(k_type)
    # A merged complex still made by a product binding stays.
    merged_v_type = @enzyme_mechanism begin
        substrates: A, B
        products: P, Q
        steps: begin
            E + A ⇌ EA               :: EqualAI
            EA + B <--> EPQ      :: OnlyA
            EQ + P ⇌ EPQ         :: EqualAI
            E + Q ⇌ EQ               :: EqualAI
        end
    end
    @test En._kept_forms_keep_a_binding(merged_v_type)
    # The active conformation is read alike: S's binding :OnlyI keeps ES, reached from
    # EP, in the active law with no binding into it that runs there. (The guard rejects
    # this one too: the active conformation turns nothing over.)
    shut_door_active = @_testhelper_unguarded @enzyme_mechanism begin
        substrates: S
        products: P
        steps: begin
            E + S <--> ES     :: OnlyI
            ES <--> EP      :: EqualAI
            E + P ⇌ EP        :: EqualAI
        end
    end
    @test !En._kept_forms_keep_a_binding(shut_door_active)
    plain_uni_uni = @enzyme_mechanism begin
        substrates: S
        products: P
        steps: begin
            E + S ⇌ ES
            ES <--> EP
            E + P ⇌ EP
        end
    end
    @test En._kept_forms_keep_a_binding(plain_uni_uni)
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

@testset "no move tags a catalytic group :OnlyI" begin
    # Aggregate check over the whole uni-uni seed set with an allosteric regulator and a
    # competitive inhibitor, and its first expansion round: no seed or child tags a
    # catalytic group :OnlyI (only a hand-written mechanism may), while the tagging move
    # still gives regulatory sites the tag.
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
    @test length(allo) == 5
    @test all(m -> :OnlyI ∉ m.allosteric.allosteric_states, allo)
    @test any(m -> any(s -> :OnlyI in s.allosteric_states,
                       m.allosteric.regulatory_sites), allo)
end

@testset "candidate_mechanisms expands degenerate seeds instead of fitting them" begin
    rxn = @enzyme_reaction begin
        substrates: A[CX], B[N]
        products: P[C], Q[NX]
    end
    seeds = unique!(
        collect(_testhelper_mechanisms(EnzymeRates.Enumeration.init_mechanisms(rxn))))
    degenerate = filter(m -> EnzymeRates.Enumeration._degenerate(_testhelper_candidate(m)),
                        seeds)
    @test length(degenerate) == 7                  # the seven ping-pong seeds
    base = [c.mechanism for c in EnzymeRates.Enumeration.candidate_mechanisms(rxn)]
    @test !any(EnzymeRates.Enumeration._degenerate, _testhelper_candidate.(base))
    cures =
        unique!([c for m in degenerate
                 for c in _testhelper_move(EnzymeRates.Enumeration._expand_re_to_ss, m)
                 if !EnzymeRates.Enumeration._degenerate(_testhelper_candidate(c))])
    @test length(cures) == 21                      # three per ping-pong seed
    # Only a flip changes which steps are at rapid equilibrium, so every other move leaves
    # a degenerate seed's child degenerate.
    kids = _testhelper_children(
        rxn, EnzymeRates.Mechanism[degenerate...])
    @test Set(filter(m -> !EnzymeRates.Enumeration._degenerate(_testhelper_candidate(m)),
                     kids)) == Set(cures)
    @test Set(base) ==
          Set(vcat(filter(m -> !EnzymeRates.Enumeration._degenerate(
                               _testhelper_candidate(m)), seeds), cures))
    @test length(base) == 278
end

@testset "_base_tier replaces a degenerate seed by the flips that cure it" begin
    # The undecorated ping-pong seed. Its steady-state chemistry EA → E(P; res) is a
    # qualifying chain, so the flip rule keeps its flanks, the A and P bindings, at rapid
    # equilibrium. The rapid-equilibrium chemistry E(B; res) ⇌ EQ is left by the B and
    # the Q binding, both rapid equilibrium, so the seed is degenerate. Each of the other
    # three groups alone raises the segment count, and each flip cures the seed: it turns
    # that chemistry steady state, or one of its two exits.
    ER = EnzymeRates
    rxn = @enzyme_reaction begin
        substrates: A[CX], B[N]
        products: P[C], Q[NX]
    end
    seed = @enzyme_mechanism begin
        substrates: A, B; products: P, Q
        steps: begin
            E + A ⇌ EA
            EA <--> EP_res_plusA_minusP
            E_res_plusA_minusP + P ⇌ EP_res_plusA_minusP
            E_res_plusA_minusP + B ⇌ EB_res_plusA_minusP
            EB_res_plusA_minusP ⇌ EQ
            E + Q ⇌ EQ
        end
    end
    chemistry_flipped = @enzyme_mechanism begin
        substrates: A, B; products: P, Q
        steps: begin
            E + A ⇌ EA
            EA <--> EP_res_plusA_minusP
            E_res_plusA_minusP + P ⇌ EP_res_plusA_minusP
            E_res_plusA_minusP + B ⇌ EB_res_plusA_minusP
            EB_res_plusA_minusP <--> EQ
            E + Q ⇌ EQ
        end
    end
    b_flipped = @enzyme_mechanism begin
        substrates: A, B; products: P, Q
        steps: begin
            E + A ⇌ EA
            EA <--> EP_res_plusA_minusP
            E_res_plusA_minusP + P ⇌ EP_res_plusA_minusP
            E_res_plusA_minusP + B <--> EB_res_plusA_minusP
            EB_res_plusA_minusP ⇌ EQ
            E + Q ⇌ EQ
        end
    end
    q_flipped = @enzyme_mechanism begin
        substrates: A, B; products: P, Q
        steps: begin
            E + A ⇌ EA
            EA <--> EP_res_plusA_minusP
            E_res_plusA_minusP + P ⇌ EP_res_plusA_minusP
            E_res_plusA_minusP + B ⇌ EB_res_plusA_minusP
            EB_res_plusA_minusP ⇌ EQ
            E + Q <--> EQ
        end
    end
    @test ER.Enumeration._degenerate(_testhelper_candidate(seed))
    base = _testhelper_mechanisms(
        ER.Enumeration._base_tier(_testhelper_candidate.([seed]), rxn))
    @test length(base) == 3
    @test Set(base) == Set([chemistry_flipped, b_flipped, q_flipped])
end

@testset "_base_tier throws a degenerate seed's expansion error" begin
    # A gives its O to B. The isomerization EA ⇌ EP turns bound A into bound P with no
    # covalent residual to hold the O, so it does not conserve atoms (EPB ⇌ EPQ makes the
    # O from nothing in turn), and the atom-conservation assertion raises on the seed's
    # flip children. The cycle takes up A and B and gives off P and Q, so the guard admits
    # the seed and its flip children. With the Q binding steady state the seed has no
    # maximal rate in products (turning that binding rapid equilibrium closes an all-RE
    # turnover cycle), so it is degenerate and expanded.
    rxn_bad = @enzyme_reaction begin
        substrates: A[CO], B[N]
        products:   P[C], Q[NO]
    end
    m_bad = @enzyme_mechanism begin
        substrates: A, B
        products: P, Q
        steps: begin
            E + A ⇌ EA
            EA ⇌ EP
            EP + B ⇌ EPB
            EPB ⇌ EPQ
            EPQ ⇌ EQ + P
            E + Q <--> EQ
        end
    end
    @test EnzymeRates.Enumeration._degenerate(_testhelper_candidate(m_bad))
    err = try
        _testhelper_mechanisms(EnzymeRates.Enumeration._base_tier(
            _testhelper_candidate.([m_bad]), rxn_bad)); nothing
    catch e
        e
    end
    @test err isa ErrorException
    # The error is the assertion's: the isomerization EA → EP loses A's O.
    @test occursin("atom-non-conserving step EA → EP", err.msg)
    @test occursin("= Dict(:O => -1)", err.msg)
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
    # count as failed and are extended: D's and Q's groups together pass. The flips of
    # B's group, P's binding, A's binding and the first chemistry step pass alone: five
    # children. Three keep the seed's degeneracy; the base tier keeps the other two. The
    # extension rests on Q's flip being refused: under the guard's first built reading the
    # ε-path oracle too finds no limit at zero P and Q. A guard that admitted Q's flip
    # would emit it alone and leave D's flip with no failed partner to extend with.
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
    p_flip = @enzyme_mechanism begin
        substrates: A, B, D ; products: P, Q ; oligomeric_state: 1
        steps: begin
            (E + B ⇌ EB, ED + B ⇌ EBD, EDQ + B ⇌ EBDQ,
             EQ_res_plusB_plusD_minusQ + B ⇌ EBQ_res_plusB_plusD_minusQ)  :: OnlyA
            (E + D ⇌ ED, EB + D ⇌ EBD, EBQ + D ⇌ EBDQ,
             EQ_res_plusB_plusD_minusQ + D ⇌ EDQ_res_plusB_plusD_minusQ)  :: EqualAI
            E + P <--> EP                                                   :: OnlyA
            E_res_plusB_plusD_minusQ + A ⇌ EA_res_plusB_plusD_minusQ        :: EqualAI
            (E_res_plusB_plusD_minusQ + Q ⇌ EQ_res_plusB_plusD_minusQ,
             EB + Q ⇌ EBQ, EBD + Q ⇌ EBDQ,
             ED + Q ⇌ EDQ)                                              :: EqualAI
            EA_res_plusB_plusD_minusQ <--> EP                              :: OnlyA
            EBD ⇌ EQ_res_plusB_plusD_minusQ                              :: OnlyA
        end
    end
    a_flip = @enzyme_mechanism begin
        substrates: A, B, D ; products: P, Q ; oligomeric_state: 1
        steps: begin
            (E + B ⇌ EB, ED + B ⇌ EBD, EDQ + B ⇌ EBDQ,
             EQ_res_plusB_plusD_minusQ + B ⇌ EBQ_res_plusB_plusD_minusQ)  :: OnlyA
            (E + D ⇌ ED, EB + D ⇌ EBD, EBQ + D ⇌ EBDQ,
             EQ_res_plusB_plusD_minusQ + D ⇌ EDQ_res_plusB_plusD_minusQ)  :: EqualAI
            E + P ⇌ EP                                                      :: OnlyA
            E_res_plusB_plusD_minusQ + A <--> EA_res_plusB_plusD_minusQ     :: EqualAI
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
    @test length(kids) == 5
    @test Set(kids) == Set([b_flip, p_flip, a_flip, chemistry_flip, dq_flip])
    @test Set(filter(m -> EnzymeRates.Enumeration._degenerate(_testhelper_candidate(m)),
                     kids)) ==
          Set([b_flip, p_flip, a_flip])
    base = _testhelper_mechanisms(
        EnzymeRates.Enumeration._base_tier(_testhelper_candidate.([seed]), rxn))
    @test length(base) == 2
    @test Set(base) == Set([chemistry_flip, dq_flip])
end

@testset "_base_tier keeps every flip child the guard admits" begin
    # Every step is rapid equilibrium, so the catalytic cycle E, ES, EP, F turns over
    # at rapid equilibrium: the seed is degenerate and expanded. The six flips of two steps
    # of the ring cut it, and each passes every face of the guard: the conformations flip
    # at E, so F needs no rapid-equilibrium link to it. No failure is recorded.
    ER = EnzymeRates
    rxn = @enzyme_reaction begin
        substrates: S[C]
        products: P[C]
    end
    seed = @_testhelper_unguarded @enzyme_mechanism begin
        substrates: S
        products: P
        steps: begin
            E + S ⇌ ES      :: EqualAI
            ES ⇌ EP       :: EqualAI
            EP ⇌ F + P      :: EqualAI
            F ⇌ E             :: EqualAI
        end
    end
    # Untyped: inside `T[...]` the parser reads `begin` as an index keyword.
    kept = [
        (@enzyme_mechanism begin
            substrates: S
            products: P
            steps: begin
                E + S <--> ES   :: EqualAI
                ES <--> EP    :: EqualAI
                EP ⇌ F + P      :: EqualAI
                F ⇌ E             :: EqualAI
            end
        end),
        (@enzyme_mechanism begin
            substrates: S
            products: P
            steps: begin
                E + S <--> ES   :: EqualAI
                ES ⇌ EP       :: EqualAI
                EP <--> F + P   :: EqualAI
                F ⇌ E             :: EqualAI
            end
        end),
        (@enzyme_mechanism begin
            substrates: S
            products: P
            steps: begin
                E + S ⇌ ES      :: EqualAI
                ES <--> EP    :: EqualAI
                EP <--> F + P   :: EqualAI
                F ⇌ E             :: EqualAI
            end
        end),
        (@enzyme_mechanism begin
            substrates: S
            products: P
            steps: begin
                E + S <--> ES   :: EqualAI
                ES ⇌ EP       :: EqualAI
                EP ⇌ F + P      :: EqualAI
                F <--> E          :: EqualAI
            end
        end),
        (@enzyme_mechanism begin
            substrates: S
            products: P
            steps: begin
                E + S ⇌ ES      :: EqualAI
                ES <--> EP    :: EqualAI
                EP ⇌ F + P      :: EqualAI
                F <--> E          :: EqualAI
            end
        end),
        @enzyme_mechanism begin
            substrates: S
            products: P
            steps: begin
                E + S ⇌ ES      :: EqualAI
                ES ⇌ EP       :: EqualAI
                EP <--> F + P   :: EqualAI
                F <--> E          :: EqualAI
            end
        end]
    @test ER.Enumeration._degenerate(_testhelper_candidate(seed))
    @test length(_testhelper_move(ER.Enumeration._expand_re_to_ss, seed)) == 6
    @test all(_testhelper_is_valid, _testhelper_move(ER.Enumeration._expand_re_to_ss, seed))
    base = _testhelper_mechanisms(
        ER.Enumeration._base_tier(_testhelper_candidate.([seed]), rxn))
    @test length(base) == 6
    @test Set(base) == Set(kept)
end

@testset "candidate_mechanisms starts from seeds that bind every required regulator" begin
    # Aggregate regression pin over a whole seed set. A required competitive inhibitor
    # makes the beam start from `seed_mechanisms`; the unregulated `init_mechanisms` seed
    # is not among them.
    ER = EnzymeRates
    rxn = @enzyme_reaction begin
        substrates: S[C]
        products: P[C]
        competitive_inhibitors: I
    end
    base = [c.mechanism for c in ER.Enumeration.candidate_mechanisms(rxn)]
    @test length(base) == 1
    @test all(m -> :I in ER.Enumeration._bound_comp_inhibitors(_testhelper_candidate(m)),
              base)
    @test Set(base) == Set(_testhelper_mechanisms(
        ER.Enumeration.seed_mechanisms(rxn, Set{Symbol}(), Set([:I]))))
end
