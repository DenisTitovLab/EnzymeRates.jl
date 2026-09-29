# ABOUTME: Tests for Steps that carry explicit consumed/released metabolite lists:
# ABOUTME: construction, canonical orientation, kinds, names, orientation-free derivation.

const ER = EnzymeRates
_testhelper_sp(bound, conf = :E) = ER.Species(ER.Metabolite[bound...], conf)
_testhelper_sp(bound, conf, res) = ER.Species(ER.Metabolite[bound...], conf, res)

@testset "Step: explicit consumed/released lists" begin
    A, B, P, Q = ER.Substrate(:A), ER.Substrate(:B), ER.Product(:P), ER.Product(:Q)
    E, EA = _testhelper_sp([]), _testhelper_sp([A])
    EQ, EAB = _testhelper_sp([Q]), _testhelper_sp([A, B])

    @testset "pure binding keeps its written orientation" begin
        s = ER.Step(E, EA, [A], ER.Metabolite[], true)
        @test ER.from_species(s) == E && ER.to_species(s) == EA
        @test ER.consumed(s) == ER.Metabolite[A] && isempty(ER.released(s))
        @test ER.ligand(s) == A && ER.is_binding(s) && !ER.is_iso(s)
    end

    @testset "a pure release is stored as the binding it reverses" begin
        s = ER.Step(EA, E, ER.Metabolite[], [A], false)
        @test s == ER.Step(E, EA, [A], ER.Metabolite[], false)
        # conformation change allowed: E*(A) → E + A is the binding E + A → E*(A)
        Estar_A = _testhelper_sp([A], :Estar)
        r = ER.Step(Estar_A, E, ER.Metabolite[], [A], true)
        @test ER.from_species(r) == E && ER.to_species(r) == Estar_A
        @test ER.ligand(r) == A
    end

    @testset "isomerization and transformations" begin
        iso = ER.Step(EAB, _testhelper_sp([P, Q]), ER.Metabolite[], ER.Metabolite[], false)
        @test ER.is_iso(iso) && ER.ligand(iso) === nothing && !ER.is_binding(iso)
        fused = ER.Step(EAB, EQ, ER.Metabolite[], [P], false)        # chemistry + release
        @test !ER.is_iso(fused) && ER.ligand(fused) === nothing
        @test ER.released(fused) == ER.Metabolite[P]
        tc = ER.Step(EA, EQ, [B], [P], false)                          # Theorell–Chance
        @test ER.consumed(tc) == ER.Metabolite[B] && ER.released(tc) == ER.Metabolite[P]
        @test ER.ligand(tc) === nothing
        two = ER.Step(E, EAB, [B, A], ER.Metabolite[], true)           # lists are sorted
        @test ER.consumed(two) == ER.Metabolite[A, B] && ER.ligand(two) === nothing
    end

    @testset "covalent residual: binding onto a residual form vs chemistry" begin
        res = ER.Residual([A], [P])
        F, FB = _testhelper_sp([], :E, res), _testhelper_sp([B], :E, res)
        @test ER.ligand(ER.Step(F, FB, [B], ER.Metabolite[], true)) == B
        chem = ER.Step(EA, F, ER.Metabolite[], [P], false)            # E(A) → F + P
        @test ER.ligand(chem) === nothing
        m = ER.Mechanism(
            @enzyme_reaction(begin
                substrates: A[CX], B[N]
                products: P[C], Q[NX]
            end),
            [[ER.Step(E, EA, [A], ER.Metabolite[], false)], [chem],
             [ER.Step(F, FB, [B], ER.Metabolite[], false)],
             [ER.Step(FB, E, ER.Metabolite[], [Q], false)]])
        @test_throws ErrorException ER._assert_chemistry_is_iso(m)
    end

    @testset "inhibitor copy stays distinct from the substrate" begin
        Ai = ER.CompetitiveInhibitor(:A)
        s_sub = ER.Step(E, EA, [A], ER.Metabolite[], true)
        s_inh = ER.Step(E, _testhelper_sp([Ai]), [Ai], ER.Metabolite[], true)
        @test s_sub != s_inh && hash(s_sub) != hash(s_inh)
        @test ER.ligand(s_inh) == Ai
    end

    @testset "rejections" begin
        @test_throws ErrorException ER.Step(E, E, [A], ER.Metabolite[], true)
        @test_throws ErrorException ER.Step(EA, EQ, [A], [A], false)
    end

    @testset "sort key reproduces today's order" begin
        s = ER.Step(E, EA, [A], ER.Metabolite[], true)
        @test ER._step_canonical_key(s) == ("E", "EA", "A", "", true)
        iso = ER.Step(EAB, _testhelper_sp([P, Q]), ER.Metabolite[], ER.Metabolite[], false)
        @test ER._step_canonical_key(iso) == ("EAB", "EPQ", "", "", false)
    end

    @testset "signature round trip" begin
        tc = ER.Step(EA, EQ, [B], [P], false)
        @test ER._step_from_sig(ER._to_sig(tc)) == tc
    end
end

@testset "transformation steps are named by form pair" begin
    m = @enzyme_mechanism begin
        substrates: A, B
        products: P, Q
        steps: begin
            E + A <--> E(A)
            E(A) + B <--> E(A, B)
            E(A, B) <--> E(Q) + P
            E(Q) <--> E + Q
        end
    end
    names = Set(ER.parameters(m, ER.Full))
    @test :k_EAB_to_EQ in names && :k_EQ_to_EAB in names
    @test !(:kon_P_EAB in names) && !(:koff_P_EAB in names)
    @test :kon_A_E in names && :koff_Q_E in names
end

# ── Brute-force mass-action oracle (no Cha segments) ────────────────────────
# A consistent parameter point: free energies g for every enzyme form and
# metabolite make each step's association constant
#   Ka = exp(g_from + Σ g_consumed − g_to − Σ g_released),
# so every cycle multiplies to Keq^n by construction. Each step is its own
# kinetic group in the mechanisms used here.
function _testhelper_consistent_point(em, rng)
    mech = ER.Mechanism(em)
    flat = ER._flat_steps(mech)
    forms = unique(vcat([[ER.from_species(s), ER.to_species(s)] for (s, _) in flat]...))
    g = Dict{Any, BigFloat}(f => 4 * big(rand(rng)) - 2 for f in forms)
    rxn = ER.reaction(mech)
    mets = vcat([ER.name(x) for x in ER.substrates(rxn)],
                [ER.name(x) for x in ER.products(rxn)],
                [ER.name(ER.regulator(r)) for r in ER.regulators(rxn)])
    for x in mets; g[x] = 4 * big(rand(rng)) - 2; end
    Keq = exp(sum(g[ER.name(x)] for x in ER.substrates(rxn)) -
              sum(g[ER.name(x)] for x in ER.products(rxn)))
    steps = []
    vals = Dict{Symbol, BigFloat}()
    sp = ER._step_parameters(mech)
    for (idx, (s, _)) in enumerate(flat)
        Ka = exp(g[ER.from_species(s)] +
                 sum((g[ER.name(m)] for m in ER.consumed(s)); init = big(0)) -
                 g[ER.to_species(s)] -
                 sum((g[ER.name(m)] for m in ER.released(s)); init = big(0)))
        if ER.is_equilibrium(s)
            p = sp[idx][1]
            vals[ER.name(p, mech)] = p isa ER.Kd ? 1 / Ka : Ka
            push!(steps, (s, Ka, nothing))
        else
            kf = exp(4 * big(rand(rng)) - 2)
            vals[ER.name(sp[idx][1], mech)] = kf
            vals[ER.name(sp[idx][2], mech)] = kf / Ka
            push!(steps, (s, Ka, kf))
        end
    end
    (; steps, vals, Keq)
end

# v per unit enzyme from the full linear steady state over every form; RE steps
# run at rate scale Λ. v = net consumption of the first substrate.
function _testhelper_mass_action_rate(em, point, concs::Dict{Symbol, BigFloat};
                                      Λ = big(10)^60)
    setprecision(BigFloat, 512) do
        mech = ER.Mechanism(em)
        forms = unique(vcat([[ER.from_species(s), ER.to_species(s)]
                             for (s, _, _) in point.steps]...))
        idx = Dict(f => i for (i, f) in enumerate(forms))
        n = length(forms)
        M = zeros(BigFloat, n, n)
        rates = []
        for (s, Ka, kf) in point.steps
            f_ = kf === nothing ? Λ : kf
            r_ = kf === nothing ? Λ / Ka : kf / Ka
            rf = f_ * prod((concs[ER.name(m)] for m in ER.consumed(s)); init = big(1))
            rr = r_ * prod((concs[ER.name(m)] for m in ER.released(s)); init = big(1))
            a, b = idx[ER.from_species(s)], idx[ER.to_species(s)]
            M[b, a] += rf; M[a, a] -= rf; M[a, b] += rr; M[b, b] -= rr
            push!(rates, (s, a, b, rf, rr))
        end
        M[n, :] .= 1
        rhs = zeros(BigFloat, n); rhs[n] = 1
        x = M \ rhs
        first_sub = ER.name(first(ER.substrates(ER.reaction(mech))))
        sum(((count(m -> ER.name(m) == first_sub, ER.consumed(s)) -
              count(m -> ER.name(m) == first_sub, ER.released(s))) * (rf * x[a] - rr * x[b])
             for (s, a, b, rf, rr) in rates); init = big(0))
    end
end

function _testhelper_check_against_mass_action(em; n_points = 3, seed = 1)
    rng = Random.MersenneTwister(seed)
    fitted = ER.fitted_params(em)
    mets = ER.metabolites(em)
    for _ in 1:n_points
        point = _testhelper_consistent_point(em, rng)
        concs = Dict{Symbol, BigFloat}(x => exp(3 * big(rand(rng)) - 1.5) for x in mets)
        params = NamedTuple{(fitted..., :Keq, :E_total)}(
            (Float64.(getindex.(Ref(point.vals), fitted))..., Float64(point.Keq), 1.0))
        concs_nt = NamedTuple{Tuple(mets)}(Tuple(Float64(concs[x]) for x in mets))
        v_pkg = ER.rate_equation(em, concs_nt, params)
        v_ref = _testhelper_mass_action_rate(em, point, concs)
        @test isapprox(v_pkg, Float64(v_ref); rtol = 1e-8)
    end
end

@testset "fused steps derive the mass-action rate" begin
    cases = [
        # uni-uni: fused release / fused binding, both orientations, both flags
        @enzyme_mechanism(begin substrates: S; products: P; steps: begin
            E + S <--> E(S); E(S) <--> E + P end end),
        @enzyme_mechanism(begin substrates: S; products: P; steps: begin
            E + S ⇌ E(S); E(S) <--> E + P end end),
        @enzyme_mechanism(begin substrates: S; products: P; steps: begin
            E + S <--> E(S); E(S) ⇌ E + P end end),
        @enzyme_mechanism(begin substrates: S; products: P; steps: begin
            E + S <--> E(S); E + P <--> E(S) end end),
        @enzyme_mechanism(begin substrates: S; products: P; steps: begin
            E + S <--> E(P); E(P) <--> E + P end end),
        @enzyme_mechanism(begin substrates: S; products: P; steps: begin
            E(P) <--> E + S; E + P <--> E(P) end end),
        @enzyme_mechanism(begin substrates: S; products: P; steps: begin
            E(P) ⇌ E + S; E(P) <--> E + P end end),
        @enzyme_mechanism(begin substrates: S; products: P; steps: begin
            E(P) <--> E + S; E(P) ⇌ E + P end end),
        # ordered bi-bi, one central complex (Segel IX-87) and mixes
        @enzyme_mechanism(begin substrates: A, B; products: P, Q; steps: begin
            E + A <--> E(A); E(A) + B <--> E(A, B); E(A, B) <--> E(Q) + P
            E(Q) <--> E + Q end end),
        @enzyme_mechanism(begin substrates: A, B; products: P, Q; steps: begin
            E + A <--> E(A); E(A) + B <--> E(A, B); E(A, B) ⇌ E(Q) + P
            E(Q) <--> E + Q end end),
        @enzyme_mechanism(begin substrates: A, B; products: P, Q; steps: begin
            E + A ⇌ E(A); E(P, Q) <--> E(A) + B; E(P, Q) <--> E(Q) + P
            E(Q) ⇌ E + Q end end),
        @enzyme_mechanism(begin substrates: A, B; products: P, Q; steps: begin
            E + A ⇌ E(A); E(P, Q) <--> E(A) + B; E(P, Q) ⇌ E(Q) + P
            E(Q) ⇌ E + Q end end),
        # ping-pong with an RE fused release
        @enzyme_mechanism(begin substrates: A, B; products: P, Q; steps: begin
            E + A <--> E(A); E(A) ⇌ E(; residual = A - P) + P
            E(; residual = A - P) + B <--> E(B; residual = A - P)
            E(B; residual = A - P) <--> E + Q end end),
        # merged random bi-bi: RE bindings and releases (former missing cut),
        # and all SS with a reversed fused binding (former mixed cut)
        @enzyme_mechanism(begin substrates: A, B; products: P, Q; steps: begin
            E + A ⇌ E(A); E + B ⇌ E(B); E(A) + B ⇌ E(A, B); E(B) + A ⇌ E(A, B)
            E(A, B) <--> E(Q) + P; E(A, B) <--> E(P) + Q; E(Q) ⇌ E + Q; E(P) ⇌ E + P
        end end),
        @enzyme_mechanism(begin substrates: A, B; products: P, Q; steps: begin
            E + A <--> E(A); E + B <--> E(B); E(A) + B <--> E(P, Q)
            E(P, Q) <--> E(B) + A; E(P, Q) <--> E(Q) + P; E(P, Q) <--> E(P) + Q
            E(Q) <--> E + Q; E(P) <--> E + P end end),
        # dead ends on a merged complex
        @enzyme_mechanism(begin substrates: A, B; products: P, Q; regulators: I
            steps: begin
            E + A ⇌ E(A); E(A) + B ⇌ E(A, B); E(A, B) <--> E(Q) + P; E(Q) ⇌ E + Q
            E + I ⇌ E(I); E(Q) + I ⇌ E(I, Q) end end),
        # ping-pong: parallel product-release route via a mixed substrate/product complex
        @enzyme_mechanism(begin substrates: A, B; products: P, Q; regulators: I
            steps: begin
            E + A ⇌ E(A); E(I) + A ⇌ E(A, I); E(Q) + A ⇌ E(A, Q)
            E + I ⇌ E(I); E(A) + I ⇌ E(A, I)
            E(P; residual = A - P) + I ⇌ E(I, P; residual = A - P)
            E(; residual = A - P) + I ⇌ E(I; residual = A - P)
            E + Q <--> E(Q); E(A) + Q <--> E(A, Q)
            E(A) <--> E(P; residual = A - P)
            E(A, I) <--> E(I, P; residual = A - P)
            E(B; residual = A - P) + P ⇌ E(B, P; residual = A - P)
            E(I; residual = A - P) + P ⇌ E(I, P; residual = A - P)
            E(; residual = A - P) + P ⇌ E(P; residual = A - P)
            E(B; residual = A - P) ⇌ E(Q)
            E(P; residual = A - P) + B ⇌ E(B, P; residual = A - P)
            E(; residual = A - P) + B ⇌ E(B; residual = A - P) end end),
    ]
    for em in cases
        _testhelper_check_against_mass_action(em)
    end
end
