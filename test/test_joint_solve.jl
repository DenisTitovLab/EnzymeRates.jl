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
        # a binding between two occupied forms before one that touches free E
        @test key(occupied, true, 1, 0) > key(plain, true, 2, 0)
        # steady state before rapid equilibrium, then the constant against the reaction
        @test key(occupied, true, 1, 0) > key(occupied_re, true, 1, 0)
        @test key(occupied, true, 2, 0) > key(occupied, true, 1, 0)
    end

    @testset "steps no chemistry orients run from the earlier form (C9)" begin
        # E ⇌ D is decided by neither atom progression nor graph context. The
        # constructor stores it by rendered name, D before E; the solve orients it from
        # free E, the first form of the canonical form order.
        m = EnzymeRates.Mechanism(@enzyme_mechanism begin
            substrates: S
            products: P
            steps: begin
                E + S <--> E(S)
                E(S) <--> E(P)
                E(P) <--> E + P
                E <--> D
            end
        end)
        fwd = EnzymeRates._solve_forward(EnzymeRates.steps(m))
        flat = EnzymeRates._flat_steps(m)
        j = only(j for (j, (s, _)) in enumerate(flat)
                 if EnzymeRates.name(EnzymeRates.from_species(s)) === :D)
        @test fwd[j] === false
        @test all(fwd[k] for k in eachindex(fwd) if k != j)
        # The stored forward rate D → E runs against the solve's direction.
        @test EnzymeRates._priority_key(first(flat[j]), false, 1, 0)[5] === true
        @test EnzymeRates._form_index_key(EnzymeRates.to_species(first(flat[j]))) <
              EnzymeRates._form_index_key(EnzymeRates.from_species(first(flat[j])))
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
    end

    @testset "every spec's solve counts its declared dimension" begin
        for spec in MECHANISM_TEST_SPECS
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

    @testset "every spec's solve satisfies detailed balance in the oracle" begin
        rng = Random.MersenneTwister(11)
        for spec in MECHANISM_TEST_SPECS
            mech = _testhelper_lift(spec.mechanism)
            point = setprecision(() -> _testhelper_solve_point(mech, rng), BigFloat, 256)
            rxn = EnzymeRates.reaction(mech)
            subs = EnzymeRates.name.(EnzymeRates.substrates(rxn))
            prods = EnzymeRates.name.(EnzymeRates.products(rxn))
            concs = Dict{Symbol, BigFloat}(x => exp(2 * big(rand(rng)) - 1)
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
