# ABOUTME: Tests of src/Enumeration/starting_mechanisms.jl: topologies, decorations, seed
# ABOUTME: variants, init_mechanisms, seed_mechanisms and candidate_mechanisms(reaction).

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

const uni_bi_rxn = @enzyme_reaction begin
    substrates: S[AB]
    products: P[A], Q[B]
end

const uni_uni_with_reg = @enzyme_reaction begin
    substrates: S[C]
    products: P[C]
    competitive_inhibitors: I
end

const ter_ter_rxn = @enzyme_reaction begin
    substrates: A[C], B[N], D[X]
    products: P[C], Q[N], R[X]
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
    # Aggregate check over the whole bi-bi starting set. It bounds the bytes a candidate
    # holds and the bytes the enumeration allocates, measured on Julia 1.11.9 and 1.13.4
    # (one expansion on 1.13.4 alone).
    # On one Julia version the bytes held are the same in every process. In the last CI
    # run, on all three runners, each candidate held 1,594 bytes and the bytes allocated
    # agreed within 2.5%. Time and peak memory depend on the machine. On the macOS CI
    # runner the next testset's ter-ter `init_mechanisms` call spent 13.9-116.1 s in garbage
    # collection in three runs of nearly the same code; ROADMAP's known issue on the slow
    # macOS tests says their cause is not known.
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
    # allocate less and pass. The 33 children are the 6 of the flip and split moves, the
    # 25 sets of the parent's five groups tagged :OnlyA that the enumerator returns, and
    # the pairs of substrate bindings and of product bindings tagged :NonequalAI, which it
    # refuses tagged :OnlyA.
    @test length(start) == 239
    @test length(children) == 33
    # 1,593.6 bytes, the same in every process on Julia 1.11.9 and 1.13.4, so the margin
    # of 3.5% only allows for an untested Julia version. A field derived from the steps
    # and stored on `Mechanism`, `Step` or the candidate's record of a form raises it: one
    # `Int` vector on `Mechanism` takes it to 1,700.
    @test held < 1_650
    # More work per candidate raises the bytes allocated. The starting set's bound sits 25%
    # above the larger of its two values, 42,545,112 bytes on 1.11 and 38,788,376 on 1.13:
    @test start_bytes < 53_000_000
    # The expansion's bound sits 3.5% above its 30,804,736 bytes on 1.13:
    @test expand_bytes < 31_883_000
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

    fixture = joinpath(@__DIR__, "..", "fixtures", "bi_bi_init_golden.txt")
    if !isfile(fixture)
        mkpath(dirname(fixture)); write(fixture, join(mech_keys, "\n"))
        @warn "bootstrapped the bi-bi starting set's golden fixture; commit it" fixture
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
    # Inline serial FIFO BFS = the reference the parallel version must match. It expands
    # every seed node and returns those the enumerator returns (`_enumerator_returns`).
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
        for m in _testhelper_mechanisms(EnzymeRates.Enumeration._base_tier(
                unique(EnzymeRates.Enumeration.init_mechanisms(rxn)), rxn))
            _testhelper_enq(m)
        end
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
        filter(m -> EnzymeRates.Enumeration._enumerator_returns(_testhelper_candidate(m)),
               seeds)
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
    # Aggregate regression pin over the whole seed set: the search reaches 4,432 seeds and
    # the enumerator drops 996 of them. 260 fail the hyperbolic rule, 60 keep a form that
    # only bindings they exclude enter (`_kept_forms_keep_a_binding`), and the guard
    # rejects 676: 236 `_BRIDGE`, 204 `_FAST_STEP_PINNED`, 188 `_FILLED` and 48 `_PUMPED`.
    @test length(seeds) == 3436
    @test all(_testhelper_is_valid, seeds)
    # The six-step ping-pong with every binding and the second half-reaction at rapid
    # equilibrium is degenerate, so the search binds R to its flip children instead
    # (`_base_tier`). Three versions of two of them, each with R::OnlyA. With both
    # half-reactions at steady state and both tagged :OnlyA the guard admits the
    # mechanism, and it is a seed. With Q's release at steady state instead, the same
    # tags fail admission's face condition (`_BRIDGE`), and the second half-reaction
    # tagged alone leaves the inactive conformation the two sides of the
    # rapid-equilibrium step EB_res_plusA_minusP ⇌ EQ joined (`_FAST_STEP_PINNED`):
    # neither is a seed. On every such flip child the guard admits the first
    # half-reaction alone :OnlyA, which it refused (`_SPLIT`) on the degenerate mechanism.
    both = @enzyme_mechanism begin
        substrates: A, B
        products: P, Q
        oligomeric_state: 2
        steps: begin
            E + A ⇌ EA                                        :: EqualAI
            EA <--> EP_res_plusA_minusP                    :: OnlyA
            EP_res_plusA_minusP ⇌ E_res_plusA_minusP + P  :: EqualAI
            E_res_plusA_minusP + B ⇌ EB_res_plusA_minusP  :: EqualAI
            EB_res_plusA_minusP <--> EQ                    :: OnlyA
            EQ ⇌ E + Q                                        :: EqualAI
        end
        regulatory_site(2): R::OnlyA
    end
    @test both ∈ seeds
    both_q_slow = @_testhelper_unguarded @enzyme_mechanism begin
        substrates: A, B
        products: P, Q
        oligomeric_state: 2
        steps: begin
            E + A ⇌ EA                                        :: EqualAI
            EA <--> EP_res_plusA_minusP                    :: OnlyA
            EP_res_plusA_minusP ⇌ E_res_plusA_minusP + P  :: EqualAI
            E_res_plusA_minusP + B ⇌ EB_res_plusA_minusP  :: EqualAI
            EB_res_plusA_minusP ⇌ EQ                       :: OnlyA
            EQ <--> E + Q                                     :: EqualAI
        end
        regulatory_site(2): R::OnlyA
    end
    @test ER.Derivation._guard(both_q_slow).fault === ER.Derivation._BRIDGE
    @test both_q_slow ∉ seeds
    half2_q_slow = @_testhelper_unguarded @enzyme_mechanism begin
        substrates: A, B
        products: P, Q
        oligomeric_state: 2
        steps: begin
            E + A ⇌ EA                                        :: EqualAI
            EA <--> EP_res_plusA_minusP                    :: EqualAI
            EP_res_plusA_minusP ⇌ E_res_plusA_minusP + P  :: EqualAI
            E_res_plusA_minusP + B ⇌ EB_res_plusA_minusP  :: EqualAI
            EB_res_plusA_minusP ⇌ EQ                       :: OnlyA
            EQ <--> E + Q                                     :: EqualAI
        end
        regulatory_site(2): R::OnlyA
    end
    v = ER.Derivation._guard(half2_q_slow)
    @test (v.fault, ER.Mechanisms.step_text(half2_q_slow.catalytic, v.edge)) ==
          (ER.Derivation._FAST_STEP_PINNED, "EB_res_plusA_minusP ⇌ EQ")
    @test half2_q_slow ∉ seeds
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

@testset "seed_mechanisms starts from the base tier" begin
    # Aggregate regression pins over two whole starting sets. The search binds R to the 278
    # mechanisms of the base tier, none of them degenerate, and the moves that bind it keep
    # every rapid-equilibrium flag, so no starting mechanism is degenerate.
    ER = EnzymeRates
    rxn = @enzyme_reaction begin
        substrates: A[CX], B[N]
        products: P[C], Q[NX]
        allosteric_regulators: R
        oligomeric_state: (1, 2)
    end
    starting = ER.Enumeration.candidate_mechanisms(rxn)
    @test length(starting) == 14078
    @test !any(ER.Enumeration._degenerate, starting)
    rxn_inhibited = @enzyme_reaction begin
        substrates: A[CX], B[N]
        products: P[C], Q[NX]
        competitive_inhibitors: I
        oligomeric_state: (1, 2)
    end
    @test length(ER.Enumeration.candidate_mechanisms(rxn_inhibited)) == 1548
end

@testset "starting mechanisms hold one half-reaction :OnlyA and a :NonequalAI pair" begin
    # With R required the starting mechanisms are grown from the base tier by the tagging
    # and regulator moves (`seed_mechanisms`), never by `_expand_change_allo_state`. On the
    # six-step ping-pong with both half-reactions at steady state they hold the version
    # whose inactive conformation stops the first half-reaction alone, and the version
    # whose inactive conformation binds A and B with constants of its own: A's and B's
    # bindings tagged :OnlyA are refused, so the tagging move tags them :NonequalAI. Each
    # binds R.
    rxn = @enzyme_reaction begin
        substrates: A[CX], B[N]
        products: P[C], Q[NX]
        allosteric_regulators: R
    end
    half1 = @enzyme_mechanism begin
        substrates: A, B ; products: P, Q ; oligomeric_state: 1
        steps: begin
            E + A ⇌ EA                                   :: EqualAI
            E_res_plusA_minusP + B ⇌ EB_res_plusA_minusP :: EqualAI
            E + Q ⇌ EQ                                   :: EqualAI
            E_res_plusA_minusP + P ⇌ EP_res_plusA_minusP :: EqualAI
            EA <--> EP_res_plusA_minusP                  :: OnlyA
            EB_res_plusA_minusP <--> EQ                  :: EqualAI
        end
        regulatory_site(1): R::OnlyA
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
        regulatory_site(1): R::OnlyA
    end
    starting = _testhelper_mechanisms(EnzymeRates.Enumeration.candidate_mechanisms(rxn))
    @test half1 in starting
    @test a_b_nonequal in starting
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
