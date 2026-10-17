# ABOUTME: Tests of the joint thermodynamic solve: the priority key and orientation, the
# ABOUTME: cycles, the exact reduction, kept members, and its constants against the oracle.
using Random

const _testhelper_JS = EnzymeRates

"""
A parameter point of `m` built from its joint solve: every column that is not a pivot
takes a random log value, `log Keq` one more, and each pivot column the value its row of
the reduced system gives. Returns the steps of each conformation as
`_testhelper_steady_state_rate` reads them, with `L` and the regulator constants drawn at
random.
"""
function _testhelper_solve_point(m, rng)
    mech = _testhelper_lift(m)
    sys = _testhelper_JS._joint_system(mech)
    M, pivots = _testhelper_JS._reduce_joint_system(sys.A, sys.rhs)
    n = length(sys.columns)
    x = zeros(BigFloat, n)
    free = setdiff(1:n, pivots)
    foreach(f -> x[f] = 4 * big(rand(rng)) - 2, free)
    log_keq = 4 * big(rand(rng)) - 2
    for (r, p) in enumerate(pivots)
        x[p] = M[r, n + 1] * log_keq -
               sum((M[r, f] * x[f] for f in free); init = big(0))
    end
    value = Dict(c => exp(x[i]) for (i, c) in enumerate(sys.columns))
    allo = mech isa _testhelper_JS.AllostericMechanism
    function _testhelper_conformation(inactive)
        out = []
        for (j, (s, g)) in enumerate(sys.flat)
            inactive && !sys.law[j] && continue
            own = inactive && _testhelper_JS._has_inactive_columns(sys.tags[g])
            col(slot) = value[_testhelper_JS._Column(j, slot, own)]
            if _testhelper_JS.is_equilibrium(s)
                K = col(1)
                push!(out, (s, _testhelper_JS.is_binding(s) ? 1 / K : K, nothing))
            else
                push!(out, (s, col(1) / col(2), col(1)))
            end
        end
        out
    end
    vals = Dict{Symbol, BigFloat}()
    L = big(0)
    if allo
        L = exp(4 * big(rand(rng)) - 2)
        for p in [_testhelper_JS._kreg_params(mech, :A);
                  _testhelper_JS._kreg_params(mech, :I)]
            vals[_testhelper_JS.name(p, mech)] = exp(4 * big(rand(rng)) - 2)
        end
    end
    (; vals, Keq = exp(log_keq), active = _testhelper_conformation(false),
     inactive = allo ? _testhelper_conformation(true) : [], L)
end

# Reversible Michaelis–Menten with every step at steady state.
const _testhelper_mm_ss = @enzyme_mechanism begin
    substrates: S
    products: P
    steps: begin
        E + S <--> E(S)
        E(S) <--> E(P)
        E(P) <--> E + P
    end
end

@testset "Joint thermodynamic solve" begin
    @testset "the priority key orders the columns" begin
        mech = EnzymeRates.Mechanism(_testhelper_mm_ss)
        sys = EnzymeRates._joint_system(mech)
        describe(c) = (EnzymeRates._step_text(sys.flat[c.step][1]), c.slot)
        # The isomerization goes first, its reverse rate (against the reaction) before
        # its forward one; then the two bindings' constants that run against the
        # reaction, the release of S and the uptake of P (a product binding's uptake runs
        # against it), the larger content key first (S sorts after P); then the others.
        @test describe.(sys.columns) == [("ES <--> EP", 2), ("ES <--> EP", 1),
                                         ("E + S <--> ES", 2), ("E + P <--> EP", 1),
                                         ("E + S <--> ES", 1), ("E + P <--> EP", 2)]
        # One Haldane row: the first column, the isomerization's reverse rate, is the
        # dependent one.
        M, pivots = EnzymeRates._reduce_joint_system(sys.A, sys.rhs)
        @test pivots == [1]
        @test M[1, end] == -1
    end

    @testset "kind, occupancy, flag and direction rank in that order" begin
        E = EnzymeRates.Species(EnzymeRates.Metabolite[], :E)
        A, B = EnzymeRates.Substrate(:A), EnzymeRates.Substrate(:B)
        EA = EnzymeRates.Species(EnzymeRates.Metabolite[A], :E)
        EAB = EnzymeRates.Species(EnzymeRates.Metabolite[A, B], :E)
        plain = EnzymeRates.Step(E, EA, [A], EnzymeRates.Metabolite[], false)
        occupied = EnzymeRates.Step(EA, EAB, [B], EnzymeRates.Metabolite[], false)
        occupied_re = EnzymeRates.Step(EA, EAB, [B], EnzymeRates.Metabolite[], true)
        key = EnzymeRates._priority_key
        # block first: an inactive copy outside the law before any other column
        @test key(plain, true, 1, 2) > key(occupied, true, 2, 1) > key(occupied, true, 2, 0)
        # kind next, ahead of every later column: an isomerization, then other chemistry (a
        # binding fused with chemistry), then a plain binding. In the first chain the later
        # columns tie or sort the other way, so kind alone orders it.
        P, Q = EnzymeRates.Product(:P), EnzymeRates.Product(:Q)
        EB = EnzymeRates.Species(EnzymeRates.Metabolite[B], :E)
        EP = EnzymeRates.Species(EnzymeRates.Metabolite[P], :E)
        EPQ = EnzymeRates.Species(EnzymeRates.Metabolite[P, Q], :E)
        F = EnzymeRates.Species(EnzymeRates.Metabolite[], :F)
        no_mets = EnzymeRates.Metabolite[]
        iso = EnzymeRates.Step(EA, EP, no_mets, no_mets, false)
        fused = EnzymeRates.Step(EA, EPQ, [B], no_mets, false)
        binding = EnzymeRates.Step(EB, EAB, [A], no_mets, false)
        free_iso = EnzymeRates.Step(E, F, no_mets, no_mets, false)
        @test key(iso, true, 2, 0) > key(fused, true, 2, 0) > key(binding, true, 2, 0)
        # kind before the constant against the reaction, and before occupancy
        @test key(iso, true, 1, 0) > key(fused, true, 2, 0)
        @test key(free_iso, true, 1, 0) > key(occupied, true, 1, 0)
        # a binding between two occupied forms before one that touches free E
        @test key(occupied, true, 1, 0) > key(plain, true, 2, 0)
        # steady state before rapid equilibrium, then the constant against the reaction
        @test key(occupied, true, 1, 0) > key(occupied_re, true, 1, 0)
        @test key(occupied, true, 2, 0) > key(occupied, true, 1, 0)
    end

    @testset "a modifier copy keys after the reactant of its name" begin
        A = EnzymeRates.Substrate(:A)
        A_copy = EnzymeRates.CompetitiveInhibitor(:A)
        @test EnzymeRates._metabolites_key([A]) < EnzymeRates._metabolites_key([A_copy])
    end

    @testset "a step neither chemistry nor graph context orients runs from the \
              earlier form (C9)" begin
        # D <--> E is decided by neither atom progression nor graph context, so it runs
        # from free E, the first form of the canonical form order, however it is
        # written: its forward rate constant is k_E_to_D.
        m = EnzymeRates.Mechanism(@enzyme_mechanism begin
            substrates: S
            products: P
            steps: begin
                E + S <--> E(S)
                E(S) <--> E(P)
                E(P) <--> E + P
                D <--> E
            end
        end)
        iso = only(s for (s, _) in EnzymeRates._flat_steps(m)
                   if EnzymeRates.conformation(EnzymeRates.to_species(s)) === :D ||
                      EnzymeRates.conformation(EnzymeRates.from_species(s)) === :D)
        @test EnzymeRates.conformation(EnzymeRates.from_species(iso)) === :E
        kf, kr = EnzymeRates._step_constants(iso, :None)
        @test EnzymeRates.name(kf, m) === :k_E_to_D
        @test EnzymeRates.name(kr, m) === :k_D_to_E
        # At rapid equilibrium the step's one constant is the equilibrium constant of that
        # direction, [D]/[E].
        re = @enzyme_mechanism begin
            substrates: S
            products: P
            steps: begin
                E + S <--> E(S)
                E(S) <--> E(P)
                E(P) <--> E + P
                D ⇌ E
            end
        end
        @test EnzymeRates.fitted_params(re) == (:K_E_to_D, :k_EP_to_E_P, :k_ES_to_EP,
                                                :k_ES_to_E_S, :k_E_P_to_EP, :k_E_S_to_ES)
        _testhelper_check_law_against_oracle(re)
    end

    @testset "a tied group turns with the lead the form order ranks first (C9)" begin
        # E <--> D and E(S) <--> D(S) share their constants, and neither atom progression
        # nor graph context orients them. On its own, E(S) <--> D(S) would run from D(S),
        # which precedes E(S) in the canonical form order, but the group turns as one with
        # its lead, the tied step ranked first by `_step_index_key`: E → D, which starts
        # at free E. So both steps run from the E side to the D side.
        tied_mm = @enzyme_mechanism begin
            substrates: S
            products: P
            steps: begin
                E + S <--> E(S)
                E(S) <--> E(P)
                E(P) <--> E + P
                D + S <--> D(S)
                (D <--> E, D(S) <--> E(S))
            end
        end
        m = EnzymeRates.Mechanism(tied_mm)
        @test EnzymeRates._step_text.(first.(EnzymeRates._flat_steps(m))) ==
              ["D + S <--> DS", "E <--> D", "ES <--> DS", "E + P <--> EP", "E + S <--> ES",
               "ES <--> EP"]
        tied = only(g for g in EnzymeRates.steps(m) if length(g) == 2)
        @test all(s -> EnzymeRates.conformation(EnzymeRates.from_species(s)) === :E &&
                       EnzymeRates.conformation(EnzymeRates.to_species(s)) === :D, tied)
        _testhelper_check_law_against_oracle(tied_mm)
    end

    @testset "fundamental cycles carry their multiple of the net reaction" begin
        m = EnzymeRates.Mechanism(@enzyme_mechanism begin
            substrates: A, B
            products: P
            steps: begin
                E + A ⇌ E(A)
                E + B ⇌ E(B)
                E(A) + B ⇌ E(A, B)
                E(B) + A ⇌ E(A, B)
                E(A, B) <--> E(P)
                E(P) ⇌ E + P
            end
        end)
        Z, multiples = EnzymeRates._fundamental_cycles(EnzymeRates.reaction(m),
                                                       EnzymeRates.steps(m))
        @test size(Z) == (2, 6)
        @test sort(abs.(multiples)) == [0, 1]
        # A cycle that turns A into P without B and Q is no multiple of the net reaction.
        slip = EnzymeRates.Mechanism(@enzyme_mechanism begin
            substrates: A, B
            products: P, Q
            steps: begin
                E + A <--> E(A)
                E(A) + B <--> E(A, B)
                E(A, B) <--> E(P, Q)
                E(P, Q) <--> E(Q) + P
                E(Q) <--> E + Q
                E(A) <--> E(P)
                E(P) <--> E + P
            end
        end)
        err = @test_throws ErrorException EnzymeRates._fundamental_cycles(
            EnzymeRates.reaction(slip), EnzymeRates.steps(slip))
        @test occursin("produces metabolite change not proportional to net reaction",
                       err.value.msg)
    end

    @testset "the reduction runs in Int64 and retries in BigInt on overflow" begin
        M, pivots = EnzymeRates._reduce_joint_system([1 -1; 1 1], [1, 0])
        @test eltype(M) == Rational{Int} && pivots == [1, 2]
        A = [2^40 1; 3 2^40]
        M, pivots = EnzymeRates._reduce_joint_system(A, [1, 0])
        @test eltype(M) == Rational{BigInt} && pivots == [1, 2]
        @test Rational{BigInt}.(A) * M[:, 3] == [1, 0]
        err = @test_throws ErrorException EnzymeRates._reduce_joint_system(
            reshape([1, 1], 2, 1), [0, 1])
        @test occursin("Thermodynamically contradictory mechanism", err.value.msg)
    end

    @testset "a class keeps its last column, a member of the inactive law first" begin
        # I binds E(A) and E(P) under one :NonequalAI group. In the inactive
        # conformation E(A) and E(A, I) lie behind the :OnlyA binding of A, so only the
        # E(P) + I member is in that conformation's law.
        am = EnzymeRates.AllostericMechanism(@allosteric_mechanism begin
            substrates: A
            products: P
            catalytic_inhibitors: I
            catalytic_multiplicity: 1
            catalytic_steps: begin
                E + A <--> E(A)                                   :: OnlyA
                E(A) <--> E(P)                                    :: OnlyA
                E + P ⇌ E(P)                                      :: EqualAI
                (E(A) + I ⇌ E(A, I), E(P) + I ⇌ E(I, P))          :: NonequalAI
            end
        end)
        kept = EnzymeRates._kept_members(am)
        g = only(g for (g, group) in enumerate(EnzymeRates.steps(am)) if length(group) == 2)
        # Active: E(A) + I sorts before E(P) + I (A before P), so it is kept; inactive:
        # the member in the law is kept.
        @test EnzymeRates._step_text(kept[(g, 1, false)]) == "EA + I ⇌ EAIinh"
        @test EnzymeRates._step_text(kept[(g, 1, true)]) == "EP + I ⇌ EIinhP"
        sys = EnzymeRates._joint_system(am)
        for (key, step) in kept
            last_col = findlast(c -> c.slot == key[2] && c.inactive == key[3] &&
                                     sys.flat[c.step][2] == key[1], sys.columns)
            @test sys.flat[sys.columns[last_col].step][1] == step
        end
        _testhelper_check_law_against_oracle(EnzymeRates.compile_mechanism(am))
    end

    @testset "the inactive cycle row runs over the whole graph, islands included" begin
        # E(A) and E(A, I) lie behind the :OnlyA binding of A, so the inactive law holds
        # only E, E(P) and E(I, P). The cycle through E(A) still gets an inactive row,
        # which reads the inactive copies of the two :OnlyA steps, outside the law.
        am = EnzymeRates.AllostericMechanism(@allosteric_mechanism begin
            substrates: A
            products: P
            catalytic_inhibitors: I
            catalytic_multiplicity: 1
            catalytic_steps: begin
                E + A <--> E(A)                                   :: OnlyA
                E(A) <--> E(P)                                    :: OnlyA
                E + P ⇌ E(P)                                      :: EqualAI
                (E(A) + I ⇌ E(A, I), E(P) + I ⇌ E(I, P))          :: NonequalAI
            end
        end)
        sys = EnzymeRates._joint_system(am)
        # The cycle once per conformation, then the second member of the :NonequalAI
        # group's active class and of its inactive class.
        @test size(sys.A) == (4, 13)
        # Row 1 is the active cycle, row 2 the inactive one. The :EqualAI binding of P has
        # no inactive constant, so the inactive row reads its active one.
        reads = Set((EnzymeRates._step_text(sys.flat[c.step][1]), c.slot, c.inactive)
                    for c in sys.columns[findall(!iszero, sys.A[2, :])])
        @test reads == Set([("E + A <--> EA", 1, true), ("E + A <--> EA", 2, true),
                            ("EA <--> EP", 1, true), ("EA <--> EP", 2, true),
                            ("E + P ⇌ EP", 1, false)])
    end

    @testset "every spec's solve counts its declared dimension" begin
        for spec in MECHANISM_TEST_SPECS
            @testset "$(spec.name)" begin
                mech = _testhelper_lift(spec.mechanism)
                sys = EnzymeRates._joint_system(mech)
                _, pivots = EnzymeRates._reduce_joint_system(sys.A, sys.rhs)
                class(c) = (g = sys.flat[c.step][2]; EnzymeRates._class(c, g, sys.tags[g]))
                last_of = Dict(class(c) => i for (i, c) in enumerate(sys.columns))
                n_law = count(i -> (c = sys.columns[i]; (!c.inactive || sys.law[c.step]) &&
                                                        last_of[class(c)] == i &&
                                                        i ∉ pivots), eachindex(sys.columns))
                n_reg = mech isa EnzymeRates.AllostericMechanism ?
                        sum((t === :NonequalAI ? 2 : 1
                             for site in EnzymeRates.regulatory_sites(mech)
                             for t in EnzymeRates.allo_states(site)); init = 0) + 1 : 0
                @test n_law + n_reg == _testhelper_declared_dimension(spec.mechanism)
            end
        end
    end

    @testset "every spec keeps the last column of each class (E3)" begin
        for spec in MECHANISM_TEST_SPECS
            @testset "$(spec.name)" begin
                mech = _testhelper_lift(spec.mechanism)
                sys = EnzymeRates._joint_system(mech)
                # In elimination order the last column of a class overwrites the earlier
                # ones. The inactive copy of an :OnlyA step is a class of its own, which
                # `_kept_members` leaves out.
                last_of_class = Dict{Tuple{Int, Int, Bool}, EnzymeRates.Step}()
                for c in sys.columns
                    g = sys.flat[c.step][2]
                    c.inactive && sys.tags[g] === :OnlyA && continue
                    last_of_class[(g, c.slot, c.inactive)] = sys.flat[c.step][1]
                end
                @test EnzymeRates._kept_members(mech) == last_of_class
            end
        end
    end

    @testset "every spec's solve satisfies detailed balance in the oracle" begin
        for (i, spec) in enumerate(MECHANISM_TEST_SPECS)
            @testset "$(spec.name)" begin
                rng = Random.MersenneTwister(i)
                mech = _testhelper_lift(spec.mechanism)
                point = setprecision(BigFloat, 256) do
                    _testhelper_solve_point(mech, rng)
                end
                rxn = EnzymeRates.reaction(mech)
                subs = EnzymeRates.name.(EnzymeRates.substrates(rxn))
                prods = EnzymeRates.name.(EnzymeRates.products(rxn))
                concs = Dict{Symbol, BigFloat}(
                    x => exp(2 * big(rand(rng)) - 1)
                    for x in EnzymeRates.metabolites(spec.mechanism))
                p1 = first(prods)
                concs[p1] = point.Keq * prod(concs[x] for x in subs) /
                            prod((concs[x] for x in prods if x != p1); init = big(1))
                off = copy(concs)
                off[p1] /= 10
                v_eq = _testhelper_steady_state_rate(mech, point, concs)
                v_off = _testhelper_steady_state_rate(mech, point, off)
                @test abs(v_eq) <= big(10)^-30 * abs(v_off)
            end
        end
    end

    @testset "the content key reads base names, never rendered names" begin
        # The form holding A and C has the smaller base-name key, since its first bound
        # name, A, sorts before AB, though its rendered name, EAC, sorts after EAB.
        m = EnzymeRates.Mechanism(@enzyme_mechanism begin
            substrates: A, C
            products: AB
            steps: begin
                E + A ⇌ E(A)
                E(A) + C ⇌ E(A, C)
                E(A, C) <--> E(AB)
                E(AB) ⇌ E + AB
            end
        end)
        forms = EnzymeRates._forms(EnzymeRates.steps(m))
        eac = only(f for f in forms if EnzymeRates.name.(EnzymeRates.bound(f)) == [:A, :C])
        eab = only(f for f in forms if EnzymeRates.name.(EnzymeRates.bound(f)) == [:AB])
        @test EnzymeRates._form_base_key(eac) < EnzymeRates._form_base_key(eab)
        @test String(EnzymeRates.name(eac)) > String(EnzymeRates.name(eab))
    end

    @testset "content keys have one concrete type whatever a form holds" begin
        # A key type that follows a form's content makes the solve's sort and comparisons
        # compile again for every new shape a search meets.
        pingpong = @enzyme_mechanism begin
            substrates: A, B
            products: P, Q
            steps: begin
                E + A <--> E(A)
                E(A) <--> E(P; residual = A - P)
                E(P; residual = A - P) ⇌ E(; residual = A - P) + P
                E(; residual = A - P) + B <--> E(Q)
                E(Q) ⇌ E + Q
            end
        end
        mechanisms = [[spec.mechanism for spec in MECHANISM_TEST_SPECS]; pingpong]
        forms = unique(f for m in mechanisms
                       for f in EnzymeRates._forms(EnzymeRates.steps(_testhelper_lift(m))))
        @test length(unique(typeof.(EnzymeRates._form_base_key.(forms)))) == 1
        @test isconcretetype(only(Base.return_types(EnzymeRates._form_base_key,
                                                    (EnzymeRates.Species,))))
        @test length(unique(typeof.(EnzymeRates._metabolites_key.(
            EnzymeRates.bound.(forms))))) == 1
    end
end
