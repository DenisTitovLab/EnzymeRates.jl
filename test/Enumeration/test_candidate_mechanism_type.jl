# ABOUTME: Tests of src/Enumeration/candidate_mechanism_type.jl: records, form names,
# ABOUTME: structural tests, size limits, _enumerator_keeps and _enumerator_returns.

# ── Records, form names and the structural tests ──────────────────────────────────────

# The kinetic groups of graph `m` whose bindings take up `x` as a modifier copy.
_testhelper_groups_binding_copy(m, x) =
    (c = _testhelper_candidate(m);
     unique(e.kinetic_group for (k, e) in enumerate(m.catalytic.steps)
            if _testhelper_binds_copy(c, k) &&
               m.catalytic.metabolites[only(e.taken_up)] === x))

# The productive twins of graph `m` over its whole graph (`_productive_twin`).
_testhelper_twin(m) = EnzymeRates.Enumeration._productive_twin(
    _testhelper_candidate(m), trues(length(m.catalytic.steps)))

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
    one_state(h, t) = ER.Enumeration._GaugeClass(h, t, 1)
    @test ER.Enumeration._gauge_rescaling(_testhelper_candidate(fused_re), members, copy_g,
                                          twin, one_state) ==
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
                                          (h, t) -> ER.Enumeration._GaugeClass(h, t, 1)) ===
          nothing
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

@testset "_enumerator_returns" begin
    En = EnzymeRates.Enumeration
    # A uni-bi whose S binds fused into EPQ at steady state and whose products leave in
    # either order. With S's binding and P's binding at E :OnlyA, the inactive conformation
    # reaches EP only through steady-state steps, and its scheme fails
    # `_hyperbolic_catalysis` (test_expansion_moves.jl, "Mechanism — a tagging child whose
    # law is not hyperbolic is not returned"). The guard rejects it too.
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
    k_pe_one = @_testhelper_unguarded @enzyme_mechanism begin
        substrates: S ; products: P, Q ; oligomeric_state: 1
        steps: begin
            E + S <--> EPQ       :: OnlyA
            EP + Q <--> EPQ    :: EqualAI
            EQ + P <--> EPQ    :: EqualAI
            E + P ⇌ EP             :: OnlyA
            E + Q ⇌ EQ             :: EqualAI
        end
    end
    @test !En._hyperbolic_where_required(_testhelper_candidate(k_pe))
    @test !En._enumerator_returns(_testhelper_candidate(k_pe))
    @test En._enumerator_returns(_testhelper_candidate(k_pe_one)) ==
          En._enumerator_keeps(k_pe_one)
    @test En._enumerator_returns(_testhelper_candidate(m)) == En._enumerator_keeps(m)
    # A uni-bi with random product release and P's group at steady state, every group
    # :NonequalAI: its scheme carries Q² (test_expansion_moves.jl, "_expand_re_to_ss: an
    # allosteric parent keeps a hyperbolic scheme"). The guard keeps it at either
    # subunit count, so at two subunits only the hyperbolic rule refuses it.
    flip_p = @enzyme_mechanism begin
        substrates: S ; products: P, Q ; oligomeric_state: 2
        steps: begin
            E + S ⇌ ES                                :: NonequalAI
            ES <--> EPQ                           :: NonequalAI
            (E + P <--> EP, EQ + P <--> EPQ)    :: NonequalAI
            (E + Q ⇌ EQ, EP + Q ⇌ EPQ)          :: NonequalAI
        end
    end
    flip_p_one = @enzyme_mechanism begin
        substrates: S ; products: P, Q ; oligomeric_state: 1
        steps: begin
            E + S ⇌ ES                                :: NonequalAI
            ES <--> EPQ                           :: NonequalAI
            (E + P <--> EP, EQ + P <--> EPQ)    :: NonequalAI
            (E + Q ⇌ EQ, EP + Q ⇌ EPQ)          :: NonequalAI
        end
    end
    @test En._enumerator_keeps(flip_p)
    @test !En._enumerator_returns(_testhelper_candidate(flip_p))
    @test En._enumerator_keeps(flip_p_one)
    @test En._enumerator_returns(_testhelper_candidate(flip_p_one))
end

# ── Redundant copies per conformation ─────────────────────────────────────────────────

@testset "_redundant_copy_groups builds each conformation from its own law" begin
    # An aggregate regression pin over the thirteen :OnlyI fixtures: the copy rule reads
    # the active conformation without its :OnlyI groups and the inactive one without its
    # :OnlyA groups, and it flags none of them.
    for (label, am) in pairs(_testhelper_onlyi_fixtures())
        @test (label, EnzymeRates.Enumeration._redundant_copy_groups(
            _testhelper_candidate(am))) == (label, Int[])
    end
end

# ── Degeneracy of an allosteric law ───────────────────────────────────────────────────

@testset "_degenerate reads the active conformation's law" begin
    # An activator X that binds only the inactive conformation opens a route whose
    # chemistry only the inactive law holds. The active law, three steady-state steps,
    # has a maximal rate both ways; read over every step, the route's product release
    # turns rapid equilibrium and closes a turnover cycle with the X complex's substrate
    # binding, so the whole graph has no reverse maximal rate. The guard admits the
    # mechanism.
    route = @enzyme_mechanism begin
        substrates: S
        products: P
        modifiers: X
        steps: begin
            E + S <--> ES         :: EqualAI
            ES <--> EP          :: EqualAI
            EP <--> E + P         :: EqualAI
            E + X ⇌ EX            :: OnlyI
            EX + S ⇌ ESX      :: EqualAI
            ESX <--> EX + P   :: EqualAI
        end
    end
    @test !EnzymeRates.Enumeration._degenerate(_testhelper_candidate(route))
    @test EnzymeRates.Enumeration._degenerate(
        _testhelper_candidate(@enzyme_mechanism begin
        substrates: S
        products: P
        modifiers: X
        steps: begin
            E + S <--> ES
            ES <--> EP
            EP <--> E + P
            E + X ⇌ EX
            EX + S ⇌ ESX
            ESX <--> EX + P
        end
    end))
    # The maximal-rate test turns rapid equilibrium every step that takes up or gives
    # off a substrate, so masking only the flags would bring the excluded steps back:
    # here the X route's two Theorell–Chance steps, each taking up S and giving off P,
    # would close a turnover cycle. Read through the edge mask, the active law keeps
    # only its three steady-state steps, which have a maximal rate both ways.
    tc_route = @_testhelper_unguarded @enzyme_mechanism begin
        substrates: S
        products: P
        modifiers: X
        steps: begin
            E + S <--> ES           :: EqualAI
            ES <--> EP            :: EqualAI
            EP <--> E + P           :: EqualAI
            E + X ⇌ EX              :: OnlyI
            EX + S <--> FX + P    :: EqualAI
            FX + S <--> EX + P    :: EqualAI
        end
    end
    cm = tc_route.catalytic
    on = EnzymeRates.Mechanisms.steps_running_in_allosteric_conformations(
        tc_route).active
    re = _testhelper_flags(cm)
    @test !EnzymeRates.Enumeration._has_vmax(cm, re .& on, cm.substrates)
    @test EnzymeRates.Enumeration._has_vmax(cm, re, cm.substrates, on)
    @test !EnzymeRates.Enumeration._degenerate(_testhelper_candidate(tc_route))
    # Read over every step, the X route's merged complex EPX sits in rapid
    # equilibrium with both sides (condition C), while the maximal rate exists both
    # ways: the steady-state FX <--> EX keeps the route from closing a
    # rapid-equilibrium turnover cycle. Read through the edge mask, the route lies
    # outside the active law.
    merged_route = @_testhelper_unguarded @enzyme_mechanism begin
        substrates: S
        products: P
        modifiers: X
        steps: begin
            E + S <--> ES           :: EqualAI
            ES <--> EP            :: EqualAI
            EP <--> E + P           :: EqualAI
            E + X ⇌ EX              :: OnlyI
            EX + S ⇌ EPX        :: EqualAI
            EPX ⇌ FX + P        :: EqualAI
            FX <--> EX            :: EqualAI
        end
    end
    cm = merged_route.catalytic
    on = EnzymeRates.Mechanisms.steps_running_in_allosteric_conformations(
        merged_route).active
    re = _testhelper_flags(cm)
    subs = cm.substrates
    @test cm.enzyme_forms[EnzymeRates.Enumeration._merged_complexes(
        _testhelper_candidate(cm), subs)] == [:EPX]
    @test isempty(
        EnzymeRates.Enumeration._merged_complexes(_testhelper_candidate(cm), subs, on))
    @test EnzymeRates.Enumeration._chemistry_equilibrates_both_sides(
        _testhelper_candidate(cm), re)
    @test !EnzymeRates.Enumeration._chemistry_equilibrates_both_sides(
        _testhelper_candidate(cm), re, on)
    @test EnzymeRates.Enumeration._has_vmax(cm, re, cm.substrates) &&
          EnzymeRates.Enumeration._has_vmax(cm, re, cm.products)
    @test !EnzymeRates.Enumeration._degenerate(_testhelper_candidate(merged_route))
    # The same route through a rapid-equilibrium Theorell–Chance step. Read over every
    # step, that step takes up S and gives off P in rapid equilibrium (condition C), and
    # the maximal rate again exists both ways; read through the edge mask, the step lies
    # outside the active law.
    re_tc_route = @_testhelper_unguarded @enzyme_mechanism begin
        substrates: S
        products: P
        modifiers: X
        steps: begin
            E + S <--> ES           :: EqualAI
            ES <--> EP            :: EqualAI
            EP <--> E + P           :: EqualAI
            E + X ⇌ EX              :: OnlyI
            EX + S ⇌ FX + P       :: EqualAI
            FX <--> EX            :: EqualAI
        end
    end
    cm = re_tc_route.catalytic
    on = EnzymeRates.Mechanisms.steps_running_in_allosteric_conformations(
        re_tc_route).active
    re = _testhelper_flags(cm)
    @test EnzymeRates.Enumeration._chemistry_equilibrates_both_sides(
        _testhelper_candidate(cm), re)
    @test !EnzymeRates.Enumeration._chemistry_equilibrates_both_sides(
        _testhelper_candidate(cm), re, on)
    @test EnzymeRates.Enumeration._has_vmax(cm, re, cm.substrates) &&
          EnzymeRates.Enumeration._has_vmax(cm, re, cm.products)
    @test !EnzymeRates.Enumeration._degenerate(_testhelper_candidate(re_tc_route))
end
