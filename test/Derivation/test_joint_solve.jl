# ABOUTME: Tests of the joint thermodynamic solve: the priority key and orientation, the
# ABOUTME: cycles, the exact reduction, kept members, and its constants against the oracle.
using Random

const _testhelper_JS = EnzymeRates

"""
A parameter point of `m` built from its joint solve: every column that is not a pivot
takes a random log value, `log Keq` one more, and each pivot column the value its row of
the reduced system gives. Returns the steps of each conformation as
`_testhelper_steady_state_rate` reads them (the steps of `_testhelper_graph`, whose order
is the solve's flat order, `_flat_edges`), with `L` and the regulator constants drawn at
random.
"""
function _testhelper_solve_point(m, rng)
    mech = _testhelper_lift(m)
    g = _testhelper_graph(mech)
    sys = _testhelper_JS.Derivation._joint_system(mech)
    edges = _testhelper_JS.Mechanisms._flat_edges(_testhelper_JS.Mechanisms.catalytic(mech))
    M, pivots = _testhelper_JS.Derivation._reduce_joint_system(sys.A, sys.rhs)
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
    allo = g.catalytic_multiplicity !== nothing
    function _testhelper_conformation(inactive)
        out = []
        for (j, s) in enumerate(g.steps)
            k = edges[j]
            inactive && !sys.law[k] && continue
            own = inactive &&
                _testhelper_JS.Mechanisms._has_inactive_columns(sys.tags[s.group])
            col(slot) = value[(edge = k, slot = slot, inactive = own)]
            if s.rapid_equilibrium
                K = col(1)
                push!(out, (s, _testhelper_is_binding(s) ? 1 / K : K, nothing))
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
        for K in [_testhelper_regulator_constants(mech, :A);
                  _testhelper_regulator_constants(mech, :I)]
            vals[K] = exp(4 * big(rand(rng)) - 2)
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
        sys = EnzymeRates.Derivation._joint_system(mech)
        describe(c) = (_testhelper_edge_text(mech, c.edge), c.slot)
        # The isomerization goes first, its reverse rate (against the reaction) before
        # its forward one; then the two bindings' constants that run against the
        # reaction, the release of S and the uptake of P (a product binding's uptake runs
        # against it), the larger content key first (S sorts after P); then the others.
        @test describe.(sys.columns) == [("ES <--> EP", 2), ("ES <--> EP", 1),
                                         ("E + S <--> ES", 2), ("E + P <--> EP", 1),
                                         ("E + S <--> ES", 1), ("E + P <--> EP", 2)]
        # One Haldane row: the first column, the isomerization's reverse rate, is the
        # dependent one.
        M, pivots = EnzymeRates.Derivation._reduce_joint_system(sys.A, sys.rhs)
        @test pivots == [1]
        @test M[1, end] == -1
    end

    @testset "kind, occupancy, flag and direction rank in that order" begin
        # Every step at steady state: a plain binding at free E, a binding between two
        # occupied forms, an isomerization, a binding fused with chemistry, a binding of A
        # at E(B) and an isomerization of free E.
        m = EnzymeRates.Mechanism(@enzyme_mechanism begin
            substrates: A, B
            products: P, Q
            steps: begin
                E + A <--> E(A)
                E(A) + B <--> E(A, B)
                E(A) <--> E(P)
                E(A) + B <--> E(P, Q)
                E + B <--> E(B)
                E(B) + A <--> E(A, B)
                E <--> F
                E(P) <--> E + P
                E(P, Q) <--> E(P) + Q
            end
        end)
        # The binding between occupied forms at rapid equilibrium.
        m_re = EnzymeRates.Mechanism(@enzyme_mechanism begin
            substrates: A, B
            products: P
            steps: begin
                E + A <--> E(A)
                E(A) + B ⇌ E(A, B)
                E(A, B) <--> E(P)
                E(P) <--> E + P
            end
        end)
        _testhelper_key(x, text, slot, block) =
            EnzymeRates.Mechanisms._priority_key(x, _testhelper_edge(x, text), slot, block,
                                                 EnzymeRates.Mechanisms._edge_base_rank(x))
        key(text, slot, block) = _testhelper_key(m, text, slot, block)
        plain, occupied = "E + A <--> EA", "EA + B <--> EAB"
        # block first: an inactive copy outside the law before any other column
        @test key(plain, 1, 2) > key(occupied, 2, 1) > key(occupied, 2, 0)
        # kind next, ahead of every later column: an isomerization, then other chemistry (a
        # binding fused with chemistry), then a plain binding. In the first chain the later
        # columns tie or sort the other way, so kind alone orders it.
        iso, fused, binding = "EA <--> EP", "EA + B <--> EPQ", "EB + A <--> EAB"
        @test key(iso, 2, 0) > key(fused, 2, 0) > key(binding, 2, 0)
        # kind before the constant against the reaction, and before occupancy
        @test key(iso, 1, 0) > key(fused, 2, 0)
        @test key("E <--> F", 1, 0) > key(occupied, 1, 0)
        # a binding between two occupied forms before one that touches free E
        @test key(occupied, 1, 0) > key(plain, 2, 0)
        # steady state before rapid equilibrium, then the constant against the reaction
        @test key(occupied, 1, 0) > _testhelper_key(m_re, "EA + B ⇌ EAB", 1, 0)
        @test key(occupied, 2, 0) > key(occupied, 1, 0)
    end

    @testset "a modifier copy keys after the reactant of its name" begin
        @test isless(EnzymeRates.Mechanisms.Met(:A, false),
                     EnzymeRates.Mechanisms.Met(:A, true))
        m = EnzymeRates.Mechanism(@enzyme_mechanism begin
            substrates: A
            products: P
            modifiers: A
            steps: begin
                E + A ⇌ E(A)
                E(A) <--> E(P)
                E(P) ⇌ E + P
                E + A::Mod ⇌ E(A::Mod)
            end
        end)
        @test m.table ==
              [EnzymeRates.Mechanisms.Met(:A, false), EnzymeRates.Mechanisms.Met(:A, true),
               EnzymeRates.Mechanisms.Met(:P, false)]
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
        iso = only(k for k in eachindex(m.edges)
                   if occursin("D", _testhelper_edge_text(m, k)))
        @test _testhelper_edge_text(m, iso) == "E <--> D"
        kf, kr = EnzymeRates.Mechanisms._edge_constants(m, iso, :None)
        @test EnzymeRates.Mechanisms.name(kf, m) === :k_E_to_D
        @test EnzymeRates.Mechanisms.name(kr, m) === :k_D_to_E
        # Written E <--> D, the step is stored the same way.
        @test EnzymeRates.Mechanism(@enzyme_mechanism begin
            substrates: S
            products: P
            steps: begin
                E + S <--> E(S)
                E(S) <--> E(P)
                E(P) <--> E + P
                E <--> D
            end
        end) == m
        _testhelper_check_stored_orientation(m)
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
        @test parameters(re, FittedParameters()) == (:K_E_to_D, :k_EP_to_E_P, :k_ES_to_EP,
                                                :k_ES_to_E_S, :k_E_P_to_EP, :k_E_S_to_ES)
        _testhelper_check_law_against_oracle(re)
        _testhelper_check_stored_orientation(re)
    end

    @testset "a tied group turns with the lead the form order ranks first (C9)" begin
        # E <--> D and E(S) <--> D(S) share their constants, and neither atom progression
        # nor graph context orients them. On its own, E(S) <--> D(S) would run from D(S),
        # which precedes E(S) in the canonical form order, but the group turns as one with
        # its lead, the tied step `_align_tied_groups!` ranks first: one that changes
        # conformation, then the lower tier, then the earlier oriented edge. Both steps
        # change conformation at tier 3, and E → D, which starts at free E, is the earlier
        # edge. So both steps run from the E side to the D side.
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
        # Groups and steps in canonical order: the steps at free E first, then by form.
        @test [_testhelper_edge_text(m, k)
               for k in EnzymeRates.Mechanisms._flat_edges(m)] ==
              ["E <--> D", "ES <--> DS", "E + P <--> EP", "E + S <--> ES", "D + S <--> DS",
               "ES <--> EP"]
        label(f) = m.conformations[m.forms[f].conformation]
        tied =
            only(ks for ks in EnzymeRates.Mechanisms._group_members(m) if length(ks) == 2)
        @test all(k -> label(m.edges[k].from) === :E && label(m.edges[k].to) === :D, tied)
        # Written E <--> D and E(S) <--> D(S), the group is stored the same way.
        @test EnzymeRates.Mechanism(@enzyme_mechanism begin
            substrates: S
            products: P
            steps: begin
                E + S <--> E(S)
                E(S) <--> E(P)
                E(P) <--> E + P
                D + S <--> D(S)
                (E <--> D, E(S) <--> D(S))
            end
        end) == m
        _testhelper_check_stored_orientation(m)
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
        Z, multiples = EnzymeRates.Derivation._fundamental_cycles(
            m, EnzymeRates.Mechanisms._flat_edges(m))
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
        err = @test_throws ErrorException EnzymeRates.Derivation._fundamental_cycles(
            slip, EnzymeRates.Mechanisms._flat_edges(slip))
        @test occursin("produces metabolite change not proportional to net reaction",
                       err.value.msg)
    end

    @testset "the reduction runs in Int64 and retries in BigInt on overflow" begin
        M, pivots = EnzymeRates.Derivation._reduce_joint_system([1 -1; 1 1], [1, 0])
        @test eltype(M) == Rational{Int} && pivots == [1, 2]
        A = [2^40 1; 3 2^40]
        M, pivots = EnzymeRates.Derivation._reduce_joint_system(A, [1, 0])
        @test eltype(M) == Rational{BigInt} && pivots == [1, 2]
        @test Rational{BigInt}.(A) * M[:, 3] == [1, 0]
        err = @test_throws ErrorException EnzymeRates.Derivation._reduce_joint_system(
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
            modifiers: I
            catalytic_multiplicity: 1
            steps: begin
                E + A <--> E(A)                                   :: OnlyA
                E(A) <--> E(P)                                    :: OnlyA
                E + P ⇌ E(P)                                      :: EqualAI
                (E(A) + I ⇌ E(A, I), E(P) + I ⇌ E(I, P))          :: NonequalAI
            end
        end)
        kept = EnzymeRates.Mechanisms._kept_members(am)
        group = am.catalytic.group
        g = only(g for (g, members) in enumerate(
            EnzymeRates.Mechanisms._group_members(am.catalytic)) if length(members) == 2)
        # Active: E(A) + I sorts before E(P) + I (A before P), so it is kept; inactive:
        # the member in the law is kept.
        @test _testhelper_edge_text(am, kept[(g, 1, false)]) == "EA + I ⇌ EAI"
        @test _testhelper_edge_text(am, kept[(g, 1, true)]) == "EP + I ⇌ EIP"
        sys = EnzymeRates.Derivation._joint_system(am)
        for (key, edge) in kept
            last_col = findlast(c -> c.slot == key[2] && c.inactive == key[3] &&
                                     group[c.edge] == key[1], sys.columns)
            @test sys.columns[last_col].edge == edge
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
            modifiers: I
            catalytic_multiplicity: 1
            steps: begin
                E + A <--> E(A)                                   :: OnlyA
                E(A) <--> E(P)                                    :: OnlyA
                E + P ⇌ E(P)                                      :: EqualAI
                (E(A) + I ⇌ E(A, I), E(P) + I ⇌ E(I, P))          :: NonequalAI
            end
        end)
        sys = EnzymeRates.Derivation._joint_system(am)
        # The cycle once per conformation, then the second member of the :NonequalAI
        # group's active class and of its inactive class.
        @test size(sys.A) == (4, 13)
        # Row 1 is the active cycle, row 2 the inactive one. The :EqualAI binding of P has
        # no inactive constant, so the inactive row reads its active one.
        reads = Set((_testhelper_edge_text(am, c.edge), c.slot, c.inactive)
                    for c in sys.columns[findall(!iszero, sys.A[2, :])])
        @test reads == Set([("E + A <--> EA", 1, true), ("E + A <--> EA", 2, true),
                            ("EA <--> EP", 1, true), ("EA <--> EP", 2, true),
                            ("E + P ⇌ EP", 1, false)])
    end

    @testset "every spec's solve counts its declared dimension" begin
        for spec in MECHANISM_TEST_SPECS
            @testset "$(spec.name)" begin
                mech = _testhelper_lift(spec.mechanism)
                sys = EnzymeRates.Derivation._joint_system(mech)
                _, pivots = EnzymeRates.Derivation._reduce_joint_system(sys.A, sys.rhs)
                class(c) = EnzymeRates.Derivation._class(mech, c)
                last_of = Dict(class(c) => i for (i, c) in enumerate(sys.columns))
                n_law = count(i -> (c = sys.columns[i]; (!c.inactive || sys.law[c.edge]) &&
                                                        last_of[class(c)] == i &&
                                                        i ∉ pivots), eachindex(sys.columns))
                n_reg = mech isa EnzymeRates.AllostericMechanism ?
                    sum((t === :NonequalAI ? 2 : 1
                         for site in EnzymeRates.Mechanisms.regulatory_sites(mech)
                         for t in EnzymeRates.Mechanisms.mwc_tags(site));
                        init = 0) + 1 : 0
                @test n_law + n_reg == _testhelper_declared_dimension(spec.mechanism)
            end
        end
    end

    @testset "every spec keeps the last column of each class (E3)" begin
        for spec in MECHANISM_TEST_SPECS
            @testset "$(spec.name)" begin
                mech = _testhelper_lift(spec.mechanism)
                sys = EnzymeRates.Derivation._joint_system(mech)
                # In elimination order the last column of a class overwrites the earlier
                # ones. The inactive copy of an :OnlyA step is a class of its own, which
                # `_kept_members` leaves out.
                group = EnzymeRates.Mechanisms.catalytic(mech).group
                last_of_class = Dict{Tuple{Int, Int, Bool}, Int}()
                for c in sys.columns
                    g = group[c.edge]
                    c.inactive && sys.tags[g] === :OnlyA && continue
                    last_of_class[(g, c.slot, c.inactive)] = c.edge
                end
                @test EnzymeRates.Mechanisms._kept_members(mech) == last_of_class
            end
        end
    end

    @testset "every spec stores its steps in the direction the orientation tiers give" begin
        # The solve reads each stored step as the physical-forward direction of its
        # constants: the rate constant against the reaction is the reverse one, and the
        # content key reads the step as stored. The constructor orients every step by
        # the tiers and turns each tied group with its lead, so orienting the stored
        # steps again leaves each one as stored, and writing every step the other way
        # round, in the reverse order, builds the same mechanism.
        # The C9 testsets above and the exchange testsets of test_types.jl check the same
        # on steps that Tier 3 orients and on tied groups.
        for spec in MECHANISM_TEST_SPECS
            @testset "$(spec.name)" begin
                _testhelper_check_stored_orientation(spec.mechanism)
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
                subs, prods = _testhelper_substrates(mech), _testhelper_products(mech)
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

    @testset "the content order reads base names, never rendered names" begin
        # The form holding A and C comes first in the content order, since its first
        # bound name, A, sorts before AB, though its rendered name, EAC, sorts after EAB.
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
        names = [EnzymeRates.Mechanisms._form_name(m, f) for f in eachindex(m.forms)]
        eac, eab = findfirst(==("EAC"), names), findfirst(==("EAB"), names)
        @test isless(m.forms[eac], m.forms[eab]) && eac < eab
        @test names[eac] > names[eab]
        # The solve's content tie ranks the graph's forms the same way (`_base_rank`).
        rank = EnzymeRates.Mechanisms._base_rank(m)
        @test rank[eac] < rank[eab]
    end

    @testset "the content order's ranks have one concrete type whatever a graph holds" begin
        # A rank type that follows a graph's content makes the sorts and comparisons that
        # read it compile again for every new shape a search meets: the solve breaks its
        # ties by integer ranks of forms and edges, which compare by content.
        @test isconcretetype(only(Base.return_types(EnzymeRates.Mechanisms._base_rank,
                                                    (EnzymeRates.Mechanism,))))
        @test isconcretetype(only(Base.return_types(EnzymeRates.Mechanisms._edge_base_rank,
                                                    (EnzymeRates.Mechanism,))))
        @test isbitstype(EnzymeRates.Mechanisms.Form) &&
              isbitstype(EnzymeRates.Mechanisms.Edge)
    end
end

@testset "_solve's half-integer dependent has keq -1//2 and exact exponents" begin
    # The conformational square ties two isomerizations of one group, so the Haldane
    # relation fixes k_F_to_E through square roots: Keq = (k_E_to_F / k_F_to_E)² ·
    # (k_F_S_to_ES / k_ES_to_F_S) · (k_FS_to_E_P / k_E_P_to_FS).
    half = EnzymeRates.Mechanism(@enzyme_mechanism begin
        substrates: S
        products: P
        steps: begin
            (E(S) <--> F(S), E <--> F)
            F(S) <--> E + P
            F + S <--> E(S)
        end
    end)
    s = EnzymeRates.Derivation._solve(half)
    @test issorted([d.col for d in s.dependents])
    @test sort([d.label for d in s.dependents]) ==
          [:haldane, :kinetic_group, :kinetic_group]
    d = only(d for d in s.dependents if d.label === :haldane)
    @test s.names[d.col] === :k_F_to_E && d.keq == -1//2
    @test Dict(s.names[f] => e for (f, e) in d.factors) ==
          Dict(:k_E_to_F => 1//1, :k_E_P_to_FS => -1//2, :k_FS_to_E_P => 1//2,
               :k_ES_to_F_S => -1//2, :k_F_S_to_ES => 1//2)
    # The rendered exponents are exact `Rational{Int}` literals (`_testhelper_power_expr`).
    @test _testhelper_dependent_param_exprs(s)[1] ==
          Dict(:k_F_to_E => :((Keq ^ $(-1//2)) * (k_ES_to_F_S ^ $(-1//2)) *
                              (k_E_P_to_FS ^ $(-1//2)) * k_E_to_F *
                              (k_FS_to_E_P ^ $(1//2)) * (k_F_S_to_ES ^ $(1//2))))
end

@testset "_solve labels kinetic-group, Haldane and Wegscheider dependents" begin
    # A binds E and E(B) under one group (a kinetic group: its occupied member's
    # constant equals the kept one); B's two bindings close a binding square with A's,
    # which ties them (Wegscheider); the cycle through the catalytic step fixes its
    # reverse rate (Haldane).
    folded = EnzymeRates.Mechanism(@enzyme_mechanism begin
        substrates: A, B
        products: P
        steps: begin
            (E + A ⇌ E(A), E(B) + A ⇌ E(A, B))
            E + B ⇌ E(B)
            E(A) + B ⇌ E(A, B)
            E(A, B) <--> E(P)
            E(P) ⇌ E + P
        end
    end)
    s = EnzymeRates.Derivation._solve(folded)
    described = [(s.names[d.col], d.label, d.keq, [s.names[f] => e for (f, e) in d.factors])
                 for d in s.dependents]
    @test described == [
        (:k_EP_to_EAB, :haldane, -1//1, [:k_EAB_to_EP => 1//1, :K_EP_to_E_P => 1//1,
                                         :K_EB_to_E_B => -1//1, :K_EA_to_E_A => -1//1]),
        (:K_EA_to_E_A, :kinetic_group, 0//1, [:K_EA_to_E_A => 1//1]),
        (:K_EAB_to_EA_B, :wegscheider, 0//1, [:K_EB_to_E_B => 1//1])]
    @test s.fitted == [:K_EA_to_E_A, :K_EB_to_E_B, :K_EP_to_E_P, :k_EAB_to_EP]
    @test s.names[s.independents] == s.fitted
    @test _testhelper_dependent_param_exprs(s) ==
          (Dict(:K_EAB_to_EA_B => :K_EB_to_E_B,
                :k_EP_to_EAB => :((1 / Keq) * (1 / K_EA_to_E_A) * (1 / K_EB_to_E_B) *
                                  K_EP_to_E_P * k_EAB_to_EP)),
           (:K_EA_to_E_A, :K_EB_to_E_B, :K_EP_to_E_P, :k_EAB_to_EP))
end

@testset "_solve maps an :EqualAI edge's inactive constant to its active column" begin
    am = EnzymeRates.AllostericMechanism(@allosteric_mechanism begin
        substrates: S
        products: P
        steps: begin
            E + S ⇌ E(S)        :: EqualAI
            E(S) <--> E + P     :: NonequalAI
        end
    end)
    s = EnzymeRates.Derivation._solve(am)
    C(edge, slot, inactive) = (; edge, slot, inactive)
    binding = _testhelper_edge(am, "E + S ⇌ ES")
    release = _testhelper_edge(am, "E + P <--> ES")
    @test length(s.columns) == 5
    @test s.column[C(binding, 1, true)] == s.column[C(binding, 1, false)]
    @test s.columns[s.column[C(binding, 1, false)]] == C(binding, 1, false)
    for slot in 1:2, inactive in (false, true)
        @test s.columns[s.column[C(release, slot, inactive)]] == C(release, slot, inactive)
    end
    @test s.fitted == [:K_ES_to_E_S, :k_A_ES_to_E_P, :k_I_ES_to_E_P, :L]
end

@testset "_solve lists the fitted names in parameter order" begin
    for spec in MECHANISM_TEST_SPECS
        @test Tuple(
            EnzymeRates.Derivation._solve(_testhelper_lift(spec.mechanism)).fitted) ==
              parameters(spec.mechanism, FittedParameters())
    end
end
