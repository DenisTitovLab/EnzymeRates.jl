# ABOUTME: Tests for mechanism enumeration pipeline
# ABOUTME: Unit tests per move + integration tests per reaction

# Test-local adapters presenting a `Step`'s form-name view: a binding
# step reads as reactants `[from_form, metabolite]` → products `[to_form]`;
# an iso step as `[source_form]` → `[dest_form]`, oriented substrate-rich
# side → product-rich side (ties by lex on form name). This matches the
# shape the enumeration assertions below were authored against; the Step
# itself canonicalizes iso direction by lex only.
function _iso_orient(s::EnzymeRates.Step)
    from, to = EnzymeRates.from_species(s), EnzymeRates.to_species(s)
    nf = count(b -> b isa EnzymeRates.Substrate, EnzymeRates.bound(from))
    nt = count(b -> b isa EnzymeRates.Substrate, EnzymeRates.bound(to))
    forward = nf > nt || (nf == nt &&
        string(EnzymeRates.name(from)) <= string(EnzymeRates.name(to)))
    forward ? (from, to) : (to, from)
end
_t_reactants(s::EnzymeRates.Step) =
    EnzymeRates.bound_metabolite(s) === nothing ?
        [EnzymeRates.name(_iso_orient(s)[1])] :
        [EnzymeRates.name(EnzymeRates.from_species(s)),
         EnzymeRates.name(EnzymeRates.bound_metabolite(s))]
_t_products(s::EnzymeRates.Step) =
    EnzymeRates.bound_metabolite(s) === nothing ?
        [EnzymeRates.name(_iso_orient(s)[2])] :
        [EnzymeRates.name(EnzymeRates.to_species(s))]

# Build a Mechanism from a flat topology Step list (each step its own
# kinetic group, in source order) — mirrors how init_mechanisms groups.
_topo_mech(rxn, t::Vector{EnzymeRates.Step}) =
    EnzymeRates.Mechanism(
        rxn, EnzymeRates._to_group_list(t, collect(1:length(t))))

# Form-name set and form→bound-metabolite-name map for a flat topology
# Step list, derived from the decomposed Species of each Step.
_form_names(t::Vector{EnzymeRates.Step}) = Set{Symbol}(
    EnzymeRates.name(sp)
    for s in t for sp in (EnzymeRates.from_species(s),
                          EnzymeRates.to_species(s)))
_boundmap(t::Vector{EnzymeRates.Step}) = Dict{Symbol, Set{Symbol}}(
    EnzymeRates.name(sp) =>
        Set(EnzymeRates.name(b) for b in EnzymeRates.bound(sp))
    for s in t for sp in (EnzymeRates.from_species(s),
                          EnzymeRates.to_species(s)))
_form_species(t::Vector{EnzymeRates.Step}) = Dict{Symbol, EnzymeRates.Species}(
    EnzymeRates.name(sp) => sp
    for s in t for sp in (EnzymeRates.from_species(s),
                          EnzymeRates.to_species(s)))

# Flat topology (Vector{Step}) from a compiled mechanism, matching the
# shape _catalytic_topologies now returns (each step its own kinetic group).
_flat_topo(m) = EnzymeRates.Step[
    s for g in EnzymeRates.Mechanism(m).steps for s in g]

# Connectivity invariant: two enzyme forms identical in conformation+residual
# whose bound-metabolite sets differ by exactly one metabolite MUST be joined by
# a binding step. Returns the list of (formA, formB) pairs that violate it.
# Accepts a flat `Vector{Step}` (a topology) or a `Vector{Vector{Step}}`
# (a Mechanism's kinetic groups).
function _connectivity_violations(steps)
    flat = eltype(steps) <: AbstractVector ?
           collect(Iterators.flatten(steps)) : steps
    _bset(sp) = Set(EnzymeRates.name(m) for m in EnzymeRates.bound(sp))
    _key(sp) = (EnzymeRates.conformation(sp), EnzymeRates.residual(sp),
                Tuple(sort(collect(_bset(sp)))))
    forms = Dict{Any,Any}()
    edges = Set{Tuple{Any,Any}}()
    for s in flat
        a, b = EnzymeRates.from_species(s), EnzymeRates.to_species(s)
        forms[_key(a)] = a; forms[_key(b)] = b
        push!(edges, (_key(a), _key(b))); push!(edges, (_key(b), _key(a)))
    end
    viol = Tuple{Symbol,Symbol}[]
    fv = collect(values(forms))
    for s1 in fv, s2 in fv
        (EnzymeRates.conformation(s1) == EnzymeRates.conformation(s2) &&
         EnzymeRates.residual(s1) == EnzymeRates.residual(s2)) || continue
        b1, b2 = _bset(s1), _bset(s2)
        if length(b2) == length(b1) + 1 && issubset(b1, b2) &&
           !((_key(s1), _key(s2)) in edges)
            push!(viol, (EnzymeRates.name(s1), EnzymeRates.name(s2)))
        end
    end
    viol
end

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
    dead_end_inhibitors: I
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
    enumerate_all_mechanism(rxn::EnzymeReaction; max_params::Int=typemax(Int))
        -> Dict{Int, Vector{Union{Mechanism, AllostericMechanism}}}

Enumerate all mechanisms reachable by init → expand, bucketed by ACTUAL
fitted-parameter count. `max_params` caps the search: mechanisms whose
actual fitted count exceeds it are dropped. Uses the same advancing-target
sweep as `_beam_search` (but expands every swept mechanism — no beam
selection) so Δ=0 expansion children (same param count as the parent) are
not lost.
"""
function enumerate_all_mechanism(rxn; max_params::Int=typemax(Int))
    M = Union{EnzymeRates.Mechanism, EnzymeRates.AllostericMechanism}
    actual(m) = length(
        EnzymeRates.fitted_params(EnzymeRates.compile_mechanism(m)))
    frontier = Dict{Int, Vector{M}}()
    function add!(m)
        pc = actual(m)
        pc <= max_params && push!(get!(frontier, pc, M[]), m)
    end
    for m in unique!(collect(EnzymeRates.init_mechanisms(rxn)))
        add!(m)
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
        for child in EnzymeRates.expand_mechanisms(swept, rxn)
            add!(child)
        end
        isempty(frontier) && break
        target = max(target + 1, minimum(keys(frontier)))
    end
    results
end

@testset "Canonical-by-construction representation independence" begin
    rxn = @enzyme_reaction begin
        substrates: S[C], A[N]
        products:   P[CN]
    end
    # One kinetic group holding two binding steps + a second group with an
    # iso step. Built two ways: reversed outer group order AND swapped inner
    # step order. Canonical-by-construction must collapse them to one struct.
    bind_S = EnzymeRates.Step(
        EnzymeRates.Species([EnzymeRates.Substrate(:S)], :E),
        EnzymeRates.Species([EnzymeRates.Substrate(:S)], :E_S),
        EnzymeRates.Substrate(:S), true)
    bind_A = EnzymeRates.Step(
        EnzymeRates.Species([EnzymeRates.Substrate(:A)], :E),
        EnzymeRates.Species([EnzymeRates.Substrate(:A)], :E_A),
        EnzymeRates.Substrate(:A), true)
    iso = EnzymeRates.Step(
        EnzymeRates.Species([EnzymeRates.Substrate(:S)], :E_S),
        EnzymeRates.Species([EnzymeRates.Product(:P)], :E_P),
        nothing, false)

    m_orderA = EnzymeRates.Mechanism(rxn, [[bind_S, bind_A], [iso]])
    m_orderB = EnzymeRates.Mechanism(rxn, [[iso], [bind_A, bind_S]])
    @test m_orderA == m_orderB
    @test hash(m_orderA) == hash(m_orderB)

    # unique! collapses the duplicate orderings; the mechanisms are not mutated.
    # m_split groups the two binding steps separately, so it is structurally
    # distinct and survives alongside the collapsed orderA/orderB.
    m_split = EnzymeRates.Mechanism(rxn, [[bind_S], [bind_A], [iso]])
    mechs = [m_orderA, m_split, m_orderB]
    snapshot = deepcopy(mechs)
    result = unique!(mechs)
    @test length(result) == 2
    @test all(r -> any(==(r), snapshot), result)
    @test m_orderA == snapshot[1]
    @test m_split  == snapshot[2]
    @test m_orderB == snapshot[3]
    @test EnzymeRates.steps(m_orderA) == EnzymeRates.steps(snapshot[1])
end

@testset "_assert_mechanism_invariants on uni-uni init" begin
    rxn = @enzyme_reaction begin
        substrates: S[C]
        products:   P[C]
    end
    m = first(EnzymeRates.init_mechanisms(rxn))
    @test EnzymeRates._assert_mechanism_invariants(m) === nothing
end

@testset "_assert_mechanism_invariants: ported coverage + group composition" begin
    # POSITIVE: an init mechanism with an unbound declared inhibitor must NOT
    # error — regulators are intentionally excluded from the coverage check
    # (init_mechanisms declares dead-end inhibitors that no step binds yet).
    rxn_inh = @enzyme_reaction begin
        substrates: S[C]
        products:   P[C]
        competitive_inhibitors: R
    end
    for m in EnzymeRates.init_mechanisms(rxn_inh)
        @test EnzymeRates._assert_mechanism_invariants(m) === nothing
    end

    # NEGATIVE 1: a declared SUBSTRATE that no step binds → error.
    rxn_unused = @enzyme_reaction begin
        substrates: S[C], T[C]
        products:   P[C2]
    end
    s1 = EnzymeRates.Step(EnzymeRates.Species([EnzymeRates.Substrate(:S)], :E),
                          EnzymeRates.Species([EnzymeRates.Substrate(:S)], :E_S),
                          EnzymeRates.Substrate(:S), true)
    s2 = EnzymeRates.Step(EnzymeRates.Species([EnzymeRates.Substrate(:S)], :E_S),
                          EnzymeRates.Species([EnzymeRates.Product(:P)], :E_P),
                          nothing, false)
    s3 = EnzymeRates.Step(EnzymeRates.Species([EnzymeRates.Product(:P)], :E_P),
                          EnzymeRates.Species(EnzymeRates.Metabolite[], :E),
                          EnzymeRates.Product(:P), true)
    m_unused = EnzymeRates.Mechanism(rxn_unused, [[s1], [s2], [s3]])
    @test_throws ErrorException EnzymeRates._assert_mechanism_invariants(m_unused)

    # NEGATIVE 2: a kinetic group binding two different metabolites → error.
    rxn2 = @enzyme_reaction begin
        substrates: S[C], A[N]
        products:   P[CN]
    end
    g1a = EnzymeRates.Step(EnzymeRates.Species([EnzymeRates.Substrate(:S)], :E),
                           EnzymeRates.Species([EnzymeRates.Substrate(:S)], :E_S),
                           EnzymeRates.Substrate(:S), true)
    g1b = EnzymeRates.Step(EnzymeRates.Species([EnzymeRates.Substrate(:A)], :E),
                           EnzymeRates.Species([EnzymeRates.Substrate(:A)], :E_A),
                           EnzymeRates.Substrate(:A), true)
    g2  = EnzymeRates.Step(EnzymeRates.Species([EnzymeRates.Substrate(:S)], :E_S),
                           EnzymeRates.Species([EnzymeRates.Product(:P)], :E_P),
                           nothing, false)
    m_mixed = EnzymeRates.Mechanism(rxn2, [[g1a, g1b], [g2]])
    @test_throws ErrorException EnzymeRates._assert_mechanism_invariants(m_mixed)
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

    # NEGATIVE: bound NADH → bound Lac, no residual (atoms don't balance).
    bad_iso = EnzymeRates.Step(
        EnzymeRates.Species([EnzymeRates.Substrate(:NADH)], :E),
        EnzymeRates.Species([EnzymeRates.Product(:Lac)], :E),
        nothing, false)
    bad_m = EnzymeRates.Mechanism(ldh_rxn, [[bad_iso]])
    @test_throws ErrorException EnzymeRates._assert_atom_conserving(bad_m)

    # POSITIVE: the bi-bi ping-pong worked example carries a real covalent
    # residual (+A −P) on conformation :E; every step conserves atoms.
    A = EnzymeRates.Substrate(:A); B = EnzymeRates.Substrate(:B)
    P = EnzymeRates.Product(:P);   Q = EnzymeRates.Product(:Q)
    res_AP = EnzymeRates.Residual([A], [P])
    E       = EnzymeRates.Species(EnzymeRates.Metabolite[], :E)
    E_A     = EnzymeRates.Species([A], :E)
    E_P_res = EnzymeRates.Species(EnzymeRates.Metabolite[P], :E, res_AP)
    F       = EnzymeRates.Species(EnzymeRates.Metabolite[], :E, res_AP)
    E_B_res = EnzymeRates.Species(EnzymeRates.Metabolite[B], :E, res_AP)
    E_Q     = EnzymeRates.Species([Q], :E)
    good_steps = [
        [EnzymeRates.Step(E, E_A, A, true)],
        [EnzymeRates.Step(E_A, E_P_res, nothing, false)],
        [EnzymeRates.Step(E_P_res, F, P, true)],
        [EnzymeRates.Step(F, E_B_res, B, true)],
        [EnzymeRates.Step(E_B_res, E_Q, nothing, true)],
        [EnzymeRates.Step(E_Q, E, Q, true)],
    ]
    good_m = EnzymeRates.Mechanism(bi_bi_pp_rxn, good_steps)
    @test EnzymeRates._assert_atom_conserving(good_m) === nothing
end

@testset "Mechanism Enumeration" begin

# ═══════════════════════════════════════════════════════════════════════
# 1. Support functions
# ═══════════════════════════════════════════════════════════════════════

# ─── _catalytic_topologies ──────────────────────────────────────────────
@testset "_catalytic_topologies" begin

    @testset "Uni-Uni" begin
        topos = EnzymeRates._catalytic_topologies(uni_uni_rxn)
        @test all(isempty(_connectivity_violations(t)) for t in topos)
        @test length(topos) == 1
        # Every topology has exactly one SS step (the iso).
        for t in topos
            @test count(
                !s.is_equilibrium for s in t) == 1
        end
    end

    @testset "Uni-Bi" begin
        topos = EnzymeRates._catalytic_topologies(uni_bi_rxn)
        @test all(isempty(_connectivity_violations(t)) for t in topos)
        @test length(topos) == 3
        for t in topos
            @test count(
                !s.is_equilibrium for s in t) == 1
            m = EnzymeMechanism(_topo_mech(uni_bi_rxn, t))
            @test m isa EnzymeMechanism
        end
    end

    @testset "Bi-Bi" begin
        topos = EnzymeRates._catalytic_topologies(bi_bi_rxn)
        @test all(isempty(_connectivity_violations(t)) for t in topos)
        # 9 sequential topologies. A[C]→P[C] and B[N]→Q[N] each leave no
        # covalent residue, so the only ping-pong is degenerate
        # (empty-residue) and is rejected by the admissible-residual rule.
        @test length(topos) == 9
        for t in topos
            @test count(
                !s.is_equilibrium for s in t) == 1
        end
    end

    @testset "Bi-Bi Ping-Pong" begin
        topos = EnzymeRates._catalytic_topologies(
            bi_bi_pp_rxn)
        @test all(isempty(_connectivity_violations(t)) for t in topos)
        @test length(topos) == 10
        for t in topos
            @test count(
                !s.is_equilibrium for s in t) == 1
        end
    end

    @testset "Ter-Ter" begin
        topos = EnzymeRates._catalytic_topologies(
            ter_ter_rxn)
        @test all(isempty(_connectivity_violations(t)) for t in topos)
        @test length(topos) == 223
        for t in topos
            @test count(
                !s.is_equilibrium for s in t) == 1
        end
    end

    @testset "Ter-Bi" begin
        topos = EnzymeRates._catalytic_topologies(
            ter_bi_rxn)
        @test all(isempty(_connectivity_violations(t)) for t in topos)
        # 45 = 39 sequential + 6 genuine-residual ping-pong. P[CN] combines
        # A+B, so ping-pong covalent residuals persist on :E; the degenerate
        # empty-residue variants are rejected by the admissible-residual rule.
        @test length(topos) == 45
        for t in topos
            @test count(
                !s.is_equilibrium for s in t) == 1
        end
    end

    @testset "admissible-residual ping-pong" begin
        # The admissible-residual rule rejects degenerate empty-residue
        # ping-pong (which would return the enzyme to apo E mid-cycle,
        # splitting the reaction into disconnected half-cycles). A genuine
        # ping-pong intermediate carries a non-empty covalent residual,
        # always on conformation :E (never a separate conformation). For
        # ter-ter, residues form by combining substrates (e.g. bind A+B,
        # release P[C] leaving an N residue), so ping-pong topologies survive.
        ter_ter = @enzyme_reaction begin
            substrates: A[C], B[N], D[X]
            products: P[C], Q[N], R[X]
        end
        topos = EnzymeRates._catalytic_topologies(ter_ter)
        @test all(isempty(_connectivity_violations(t)) for t in topos)
        # Every enzyme form lives on conformation :E.
        for t in topos, s in t,
            sp in (EnzymeRates.from_species(s), EnzymeRates.to_species(s))
            @test EnzymeRates.conformation(sp) === :E
        end
        # Surviving ping-pong topologies each carry a genuine (non-empty)
        # covalent intermediate — no empty-residue ping-pong remains.
        pingpong = filter(t -> count(EnzymeRates.is_iso, t) >= 2, topos)
        @test !isempty(pingpong)
        for t in pingpong
            @test any(
                EnzymeRates.has_residual(sp)
                for s in t
                for sp in (EnzymeRates.from_species(s),
                           EnzymeRates.to_species(s)))
        end
    end

    @testset "weak-ordering combining" begin
        # For bi-bi: 9 sequential (the degenerate empty-residue ping-pong
        # is rejected by the admissible-residual rule).
        bi_bi_rxn_test = @enzyme_reaction begin
            substrates: A[C], B[N]
            products: P[C], Q[N]
        end
        topos = EnzymeRates._catalytic_topologies(
            bi_bi_rxn_test)
        @test all(isempty(_connectivity_violations(t)) for t in topos)
        @test length(topos) == 9

        topos_tt = EnzymeRates._catalytic_topologies(
            ter_ter_rxn)
        @test all(isempty(_connectivity_violations(t)) for t in topos_tt)
        @test length(topos_tt) == 223
    end

    @testset "isomerization constraints" begin
        topos = EnzymeRates._catalytic_topologies(
            ter_ter_rxn)
        @test all(isempty(_connectivity_violations(t)) for t in topos)

        sub_names_set = Set([:A, :B, :D])

        # C5: at most max(n_subs, n_prods) = 3 metabolites bound on any form.
        for spec in topos
            for s in spec
                for sp in (EnzymeRates.from_species(s),
                           EnzymeRates.to_species(s))
                    @test length(EnzymeRates.bound(sp)) <= 3
                end
            end
        end

        # C7: every iso source must contain at
        # least one substrate
        for spec in topos
            for s in spec
                if length(_t_reactants(s)) == 1 &&
                        length(_t_products(s)) == 1
                    # Use bound list directly (name-parsing is ambiguous
                    # with the concat form naming convention).
                    src_sp, _ = _iso_orient(s)
                    has_sub = any(
                        b -> EnzymeRates.name(b) ∈ sub_names_set,
                        EnzymeRates.bound(src_sp))
                    @test has_sub
                end
            end
        end

        # C8: an iso's product-side (destination) form is built with only
        # products bound, never substrates.
        for spec in topos
            for s in spec
                EnzymeRates.is_iso(s) || continue
                dst = _iso_orient(s)[2]
                for b in EnzymeRates.bound(dst)
                    @test !(EnzymeRates.name(b) ∈ sub_names_set)
                end
            end
        end
    end

    @testset "pyruvate carboxylase mechanism" begin
        topos = EnzymeRates._catalytic_topologies(
            pyruvate_carboxylase_rxn)
        @test all(isempty(_connectivity_violations(t)) for t in topos)

        # Known mechanism: ATP+HCO3 → ADP+Pi leaving a CO2 covalent
        # residual on :E, then Pyr+CO2 → OAA. The carboxylation iso converts
        # the {ATP,HCO3}-bound form into a residual-bearing form; the
        # carboxyl-transfer iso converts a residual-bearing Pyr form into the
        # bare E(OAA).
        _bset(sp) = Set(EnzymeRates.name(b) for b in EnzymeRates.bound(sp))
        found = false
        for spec in topos
            has_carboxylation = any(spec) do s
                EnzymeRates.is_iso(s) || return false
                f, t = EnzymeRates.from_species(s), EnzymeRates.to_species(s)
                (_bset(f) == Set([:ATP, :HCO3]) &&
                    !EnzymeRates.has_residual(f) &&
                    EnzymeRates.has_residual(t)) ||
                (_bset(t) == Set([:ATP, :HCO3]) &&
                    !EnzymeRates.has_residual(t) &&
                    EnzymeRates.has_residual(f))
            end
            has_carboxyl_transfer = any(spec) do s
                EnzymeRates.is_iso(s) || return false
                f, t = EnzymeRates.from_species(s), EnzymeRates.to_species(s)
                (_bset(f) == Set([:OAA]) && !EnzymeRates.has_residual(f) &&
                    :Pyr in _bset(t) && EnzymeRates.has_residual(t)) ||
                (_bset(t) == Set([:OAA]) && !EnzymeRates.has_residual(t) &&
                    :Pyr in _bset(f) && EnzymeRates.has_residual(f))
            end
            if has_carboxylation && has_carboxyl_transfer
                found = true
                break
            end
        end
        @test found

        # 312 = 169 seq + 143 pp, classified by iso-step count: sequential
        # topologies have one iso step, ping-pong ≥2. Every topology (seq or
        # pp) has exactly one SS step — for ping-pong only one of the iso
        # steps is steady-state, the rest are rapid-equilibrium.
        @test length(topos) == 312
        seq_count = count(t -> count(EnzymeRates.is_iso, t) == 1, topos)
        pp_count = length(topos) - seq_count
        @test seq_count == 169
        @test pp_count == 143
    end

    @testset "pyruvate dehydrogenase mechanism" begin
        topos = EnzymeRates._catalytic_topologies(
            pyruvate_dehydrogenase_rxn)
        @test all(isempty(_connectivity_violations(t)) for t in topos)

        # Known mechanism, with covalent residuals on :E:
        # Pyr→CO2 (leaves an acetyl residual),
        # CoA+acetyl→AcCoA (leaves a hydride residual),
        # NAD+hydride→NADH (residual cancels → bare E(NADH)).
        _bset(sp) = Set(EnzymeRates.name(b) for b in EnzymeRates.bound(sp))
        found = false
        for spec in topos
            has_pyr = any(spec) do s
                EnzymeRates.is_iso(s) || return false
                f, t = EnzymeRates.from_species(s), EnzymeRates.to_species(s)
                (_bset(f) == Set([:Pyr]) && !EnzymeRates.has_residual(f) &&
                    _bset(t) == Set([:CO2]) && EnzymeRates.has_residual(t)) ||
                (_bset(t) == Set([:Pyr]) && !EnzymeRates.has_residual(t) &&
                    _bset(f) == Set([:CO2]) && EnzymeRates.has_residual(f))
            end
            has_coa = any(spec) do s
                EnzymeRates.is_iso(s) || return false
                f, t = EnzymeRates.from_species(s), EnzymeRates.to_species(s)
                (_bset(f) == Set([:CoA]) && EnzymeRates.has_residual(f) &&
                    _bset(t) == Set([:AcCoA]) && EnzymeRates.has_residual(t)) ||
                (_bset(t) == Set([:CoA]) && EnzymeRates.has_residual(t) &&
                    _bset(f) == Set([:AcCoA]) && EnzymeRates.has_residual(f))
            end
            has_nad = any(spec) do s
                EnzymeRates.is_iso(s) || return false
                f, t = EnzymeRates.from_species(s), EnzymeRates.to_species(s)
                (_bset(f) == Set([:NAD]) && EnzymeRates.has_residual(f) &&
                    _bset(t) == Set([:NADH]) && !EnzymeRates.has_residual(t)) ||
                (_bset(t) == Set([:NAD]) && EnzymeRates.has_residual(t) &&
                    _bset(f) == Set([:NADH]) && !EnzymeRates.has_residual(f))
            end
            if has_pyr && has_coa && has_nad
                found = true
                break
            end
        end
        @test found

        # 334 = 169 seq + 165 pp, classified by iso-step count: sequential
        # topologies have one iso step, ping-pong ≥2. Every topology (seq or
        # pp) has exactly one SS step — for ping-pong only one of the iso
        # steps is steady-state, the rest are rapid-equilibrium.
        @test length(topos) == 334
        seq_count = count(t -> count(EnzymeRates.is_iso, t) == 1, topos)
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
        topos = EnzymeRates._catalytic_topologies(quad_rxn)
        @test all(isempty(_connectivity_violations(t)) for t in topos)
        @test length(topos) > 0
        # Every topology must have ≥ 2 iso steps (no 4→4)
        for spec in topos
            n_iso = count(spec) do s
                length(_t_reactants(s)) == 1 &&
                    length(_t_products(s)) == 1
            end
            @test n_iso >= 2
        end
    end

end

# ─── _competition_patterns ──────────────────────────────────────────────
@testset "_competition_patterns" begin
    # Uni-uni: 1×1, only 1 pattern (single edge)
    pats_11 = EnzymeRates._competition_patterns(
        Set([:S]), Set([:P]))
    @test length(pats_11) == 1
    @test pats_11[1] == Set([(:S, :P)])

    # Uni-bi: 1×2, S competes with both P and Q
    pats_12 = EnzymeRates._competition_patterns(
        Set([:S]), Set([:P, :Q]))
    @test length(pats_12) == 1
    @test pats_12[1] == Set([(:S, :P), (:S, :Q)])

    # Bi-uni: symmetric
    pats_21 = EnzymeRates._competition_patterns(
        Set([:A, :B]), Set([:P]))
    @test length(pats_21) == 1
    @test pats_21[1] == Set([(:A, :P), (:B, :P)])

    # Bi-bi: 7 patterns
    pats_22 = EnzymeRates._competition_patterns(
        Set([:A, :B]), Set([:P, :Q]))
    @test length(pats_22) == 7
    # Every pattern covers all vertices
    for pat in pats_22
        for s in [:A, :B]
            @test any(
                p -> (s, p) in pat, [:P, :Q])
        end
        for p in [:P, :Q]
            @test any(
                s -> (s, p) in pat, [:A, :B])
        end
    end
    # Invalid: {A↔P, B↔P} leaves Q uncovered
    @test Set([(:A, :P), (:B, :P)]) ∉ pats_22

    # Ter-ter: 265 patterns
    pats_33 = EnzymeRates._competition_patterns(
        Set([:A, :B, :C]),
        Set([:P, :Q, :R]))
    @test length(pats_33) == 265
    for pat in pats_33
        for s in [:A, :B, :C]
            @test any(
                p -> (s, p) in pat,
                [:P, :Q, :R])
        end
        for p in [:P, :Q, :R]
            @test any(
                s -> (s, p) in pat,
                [:A, :B, :C])
        end
    end

    @testset "Asymmetric: 2 × 3" begin
        # 2 substrates × 3 products. Bipartite-cover count on K(2,3)
        # by inclusion-exclusion = 25.
        pats = EnzymeRates._competition_patterns(Set([:A, :B]), Set([:P, :Q, :R]))
        @test length(pats) == 25
        for pat in pats
            for s in [:A, :B]
                @test any((s, p) in pat for p in [:P, :Q, :R])
            end
            for p in [:P, :Q, :R]
                @test any((s, p) in pat for s in [:A, :B])
            end
        end
    end

    @testset "Asymmetric: 3 × 2" begin
        # By symmetry with 2 × 3, count is also 25.
        pats = EnzymeRates._competition_patterns(Set([:A, :B, :D]), Set([:P, :Q]))
        @test length(pats) == 25
        for pat in pats
            for s in [:A, :B, :D]
                @test any((s, p) in pat for p in [:P, :Q])
            end
            for p in [:P, :Q]
                @test any((s, p) in pat for s in [:A, :B, :D])
            end
        end
    end
end

# ─── _inhibitor_competition_patterns ────────────────────────────────────
@testset "_inhibitor_competition_patterns" begin
    # Uni-uni, no existing inhibitors
    pats = EnzymeRates._inhibitor_competition_patterns(
        Set([:S]), Set([:P]), Symbol[])
    @test length(pats) == 1
    @test pats[1] == (Set([:S]), Set([:P]), Set{Symbol}())

    # Bi-bi, no existing inhibitors: 3×3 = 9
    pats_bb = EnzymeRates._inhibitor_competition_patterns(
        Set([:A, :B]), Set([:P, :Q]), Symbol[])
    @test length(pats_bb) == 9

    # Ter-ter, no existing inhibitors: 7×7 = 49
    pats_tt = EnzymeRates._inhibitor_competition_patterns(
        Set([:A, :B, :C]), Set([:P, :Q, :R]), Symbol[])
    @test length(pats_tt) == 49

    # Bi-bi, 1 existing inhibitor: 9 × 2 = 18
    pats_1i = EnzymeRates._inhibitor_competition_patterns(
        Set([:A, :B]), Set([:P, :Q]), [:I1__reg])
    @test length(pats_1i) == 18

    # Bi-bi, 2 existing inhibitors: 9 × 4 = 36
    pats_2i = EnzymeRates._inhibitor_competition_patterns(
        Set([:A, :B]), Set([:P, :Q]),
        [:I1__reg, :I2__reg])
    @test length(pats_2i) == 36

    @testset "Bi-bi with 3 existing inhibitors: 9 × 8 = 72" begin
        # 3 existing inhibitors → 2^3 = 8 inhibitor-competition combinations.
        # Combined with 9 base patterns → 9 × 8 = 72 variants.
        pats = EnzymeRates._inhibitor_competition_patterns(
            Set([:A, :B]), Set([:P, :Q]),
            [:I1__reg, :I2__reg, :I3__reg])
        @test length(pats) == 72
    end
end

# ─── _forms_with_binding_step ───────────────────────────────────────────
@testset "_forms_with_binding_step" begin
    # Uni-uni: S binds to E, P binds to E
    m_uu = @enzyme_mechanism begin
        substrates: S
        products: P
        steps: begin
            E + P ⇌ E(P)
            E + S ⇌ E(S)
            E(S) <--> E(P)
        end
    end
    mech_uu = EnzymeRates.Mechanism(m_uu)
    @test EnzymeRates._forms_with_binding_step_native(
        mech_uu, :S) == Set([:E])
    @test EnzymeRates._forms_with_binding_step_native(
        mech_uu, :P) == Set([:E])

    # Bi-bi random: B binds to E and E_A
    m_bb = @enzyme_mechanism begin
        substrates: A, B
        products: P, Q
        steps: begin
            E + A ⇌ E(A)
            E(B) + A ⇌ E(A, B)
            E + B ⇌ E(B)
            E(A) + B ⇌ E(A, B)
            E + P ⇌ E(P)
            E(P) + Q ⇌ E(P, Q)
            E + Q ⇌ E(Q)
            E(Q) + P ⇌ E(P, Q)
            E(A, B) <--> E(P, Q)
        end
    end
    mech_bb = EnzymeRates.Mechanism(m_bb)
    @test EnzymeRates._forms_with_binding_step_native(
        mech_bb, :B) == Set([:E, :EA])
    @test EnzymeRates._forms_with_binding_step_native(
        mech_bb, :A) == Set([:E, :EB])
    @test EnzymeRates._forms_with_binding_step_native(
        mech_bb, :P) == Set([:E, :EQ])
    @test EnzymeRates._forms_with_binding_step_native(
        mech_bb, :Q) == Set([:E, :EP])
end

# ─── _substrate_product_dead_end_opportunities ──────────────────────────
@testset "_substrate_product_dead_end_opportunities" begin
    # Random ter-ter topology has 27 possible
    # dead-end forms. With diagonal competition
    # {A↔P, B↔Q, D↔R}:
    #   1S+1P: 6 allowed, 3 forbidden (the 3
    #     diagonal pairs A-P, B-Q, D-R)
    #   2S+1P: 3 allowed (EABR, EADQ,
    #     EBDP), 6 forbidden
    #   1S+2P: 3 allowed (EAQR, EBPR,
    #     EDPQ), 6 forbidden
    #   Total: 12 allowed out of 27
    topos = EnzymeRates._catalytic_topologies(
        ter_ter_rxn)
    # Pick a sequential topo with most forms
    _, idx = findmax(
        length(_form_names(t))
        for t in topos)
    random_topo = topos[idx]
    bound = _boundmap(random_topo)
    form_sp = _form_species(random_topo)
    sub_names = Set([:A, :B, :D])
    prod_names = Set([:P, :Q, :R])
    cat_forms = _form_names(random_topo)
    _role(m) = m in sub_names ? EnzymeRates.Substrate(m) :
                                EnzymeRates.Product(m)
    _add(sp, m) = EnzymeRates.Species(
        EnzymeRates.Metabolite[EnzymeRates.bound(sp)..., _role(m)],
        EnzymeRates.conformation(sp), EnzymeRates.residual(sp))
    _sp_de_opps =
        EnzymeRates._substrate_product_dead_end_opportunities
    de_opps = _sp_de_opps(
        form_sp, bound, cat_forms,
        sub_names, prod_names, _add)
    # Group dead-end forms
    de_forms = Dict{Symbol,
        Vector{Tuple{Symbol, Symbol}}}()
    for (f, m) in de_opps
        de_name = EnzymeRates.name(_add(form_sp[f], m))
        push!(get!(de_forms, de_name,
            Tuple{Symbol, Symbol}[]), (f, m))
    end
    de_form_names =
        sort(collect(keys(de_forms)))
    @test length(de_form_names) == 27

    # Build de_bound mapping
    de_bound = Dict{Symbol, Set{Symbol}}()
    for de_name in de_form_names
        f, m = first(de_forms[de_name])
        de_bound[de_name] = union(
            bound[f], Set([m]))
    end

    # Apply diagonal competition filter
    diagonal =
        Set([(:A, :P), (:B, :Q), (:D, :R)])
    allowed = Symbol[]
    for de_name in de_form_names
        mets = de_bound[de_name]
        de_subs = intersect(mets, sub_names)
        de_prods =
            intersect(mets, prod_names)
        has_conflict = any(
            (s, p) in diagonal
            for s in de_subs
            for p in de_prods)
        has_conflict ||
            push!(allowed, de_name)
    end
    @test length(allowed) == 12

    # Verify specific allowed forms
    @test :EAQ in allowed   # 1S+1P
    @test :EAR in allowed
    @test :EBP in allowed
    @test :EBR in allowed
    @test :EDP in allowed
    @test :EDQ in allowed
    @test :EABR in allowed # 2S+1P
    @test :EADQ in allowed
    @test :EBDP in allowed
    @test :EAQR in allowed # 1S+2P
    @test :EBPR in allowed
    @test :EDPQ in allowed

    # Verify specific forbidden forms
    @test :EAP ∉ allowed    # A↔P diagonal
    @test :EBQ ∉ allowed    # B↔Q diagonal
    @test :EDR ∉ allowed    # D↔R diagonal
end

# ─── _expand_substrate_product_dead_ends ────────────────────────────────
@testset "_expand_substrate_product_dead_ends" begin

    @testset "Uni-Uni: no dead-end forms" begin
        # 3 forms: E, E_S[C], E_P[C]. E_S has all subs,
        # E_P has all prods. No mixed dead-end possible.
        # → 0 dead-end forms, 1 variant (bare topology)
        m = @enzyme_mechanism begin
            substrates: S
            products: P
            steps: begin
                E + P ⇌ E(P)
                E + S ⇌ E(S)
                E(S) <--> E(P)
            end
        end
        topo = _flat_topo(m)
        result =
            EnzymeRates._expand_substrate_product_dead_ends(
                [topo],uni_uni_rxn)
        @test all(isempty(_connectivity_violations(steps))
                  for (steps, _groups) in result)
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
                E + A ⇌ E(A)
                E(B) + A ⇌ E(A, B)
                E + B ⇌ E(B)
                E(A) + B ⇌ E(A, B)
                E + P ⇌ E(P)
                E(P) + Q ⇌ E(P, Q)
                E + Q ⇌ E(Q)
                E(Q) + P ⇌ E(P, Q)
                E(A, B) <--> E(P, Q)
            end
        end
        topo = _flat_topo(m)
        result =
            EnzymeRates._expand_substrate_product_dead_ends(
                [topo],bi_bi_rxn)
        @test all(isempty(_connectivity_violations(steps))
                  for (steps, _groups) in result)
        # 4 unique dead-end forms, 7 competition patterns,
        # all 7 produce distinct dead-end sets → 7 variants
        @test length(result) == 7
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
                E + Q ⇌ E(Q)
                E(Q) + P ⇌ E(P, Q)
                E + S ⇌ E(S)
                E(S) <--> E(P, Q)
            end
        end
        topo = _flat_topo(m)
        result =
            EnzymeRates._expand_substrate_product_dead_ends(
                [topo],uni_bi_rxn)
        @test all(isempty(_connectivity_violations(steps))
                  for (steps, _groups) in result)
        @test length(result) == 1
    end

    @testset "Bi-Bi Ping-Pong: 5 dead-end forms → 7 variants" begin
        # Forms: E, E_A, Estar, Estar_A_P, Estar_B, E_Q
        # 5 dead-end forms total (E-side: E_A_P, E_A_Q, E_B_Q; Estar-side:
        # Estar_B_P, Estar_B_Q). 7 competition patterns; each yields a
        # distinct dead-end-form set after dedup → 7 variants.
        m = @enzyme_mechanism begin
            substrates: A, B
            products: P, Q
            steps: begin
                E + A ⇌ E(A)
                Estar + B ⇌ Estar(B)
                E + Q ⇌ E(Q)
                Estar + P ⇌ Estar(A, P)
                E(A) <--> Estar(A, P)
                Estar(B) ⇌ E(Q)
            end
        end
        topo = _flat_topo(m)
        result =
            EnzymeRates._expand_substrate_product_dead_ends(
                [topo],bi_bi_pp_rxn)
        @test all(isempty(_connectivity_violations(steps))
                  for (steps, _groups) in result)
        # 5 dead-end forms (E_A_P, E_A_Q, E_B_Q from
        # E-side + Estar_B_P, Estar_B_Q from
        # Estar-side), competition-filtered
        @test length(result) == 7

        # Assert that some result variants contain Estar-prefixed dead-end
        # forms (proving dead-end forms inherit the base form's Estar
        # conformation).
        seed_forms = _form_names(topo)
        new_estar_forms = Set{Symbol}()
        for r in result
            new_forms = setdiff(_form_names(r[1]), seed_forms)
            for f in new_forms
                startswith(string(f), "Estar") && string(f) != "Estar" && push!(new_estar_forms, f)
            end
        end
        @test !isempty(new_estar_forms)
    end

    @testset "Dead-end filtering by competition" begin

        # Shared bi-bi random mechanism for multiple tests
        m_bb = @enzyme_mechanism begin
            substrates: A, B
            products: P, Q
            steps: begin
                E + A ⇌ E(A)
                E(B) + A ⇌ E(A, B)
                E + B ⇌ E(B)
                E(A) + B ⇌ E(A, B)
                E + P ⇌ E(P)
                E(P) + Q ⇌ E(P, Q)
                E + Q ⇌ E(Q)
                E(Q) + P ⇌ E(P, Q)
                E(A, B) <--> E(P, Q)
            end
        end
        topo_bb = _flat_topo(m_bb)

        # Note: the variant count assertion (length == 7) is covered by
        # the "Bi-Bi random: 4 dead-end forms" sub-testset higher in this
        # same parent testset. The sub-testsets below probe the SHAPE of
        # those 7 variants — which forms appear in each, and which
        # patterns produce empty/full dead-end sets.

        @testset "Bi-bi random: complete competition → bare topology" begin
            result =
                EnzymeRates._expand_substrate_product_dead_ends(
                    [topo_bb],bi_bi_rxn)
            @test all(isempty(_connectivity_violations(steps))
                      for (steps, _groups) in result)
            # Complete pattern {A↔P,A↔Q,B↔P,B↔Q} forbids
            # all dead-end forms → 1 variant has no dead-end
            # steps (same step count as original)
            bare = filter(
                r -> length(r[1]) == length(topo_bb),
                result)
            @test length(bare) == 1
        end

        @testset "Bi-bi random: diagonal has exactly 2 dead-end forms" begin
            result =
                EnzymeRates._expand_substrate_product_dead_ends(
                    [topo_bb],bi_bi_rxn)
            @test all(isempty(_connectivity_violations(steps))
                      for (steps, _groups) in result)
            # Diagonal patterns {A↔P,B↔Q} and {A↔Q,B↔P}
            # each allow exactly 2 dead-end forms.
            two_de = filter(result) do r
                de_forms = setdiff(
                    _form_names(r[1]),
                    _form_names(topo_bb))
                length(de_forms) == 2
            end
            @test length(two_de) == 2  # diagonal + anti-diagonal
        end

        @testset "Ter-ter per-topology (OOM on full init)" begin
            # Test that competition filtering works
            # on representative ter-ter topologies.
            topos = EnzymeRates._catalytic_topologies(
                ter_ter_rxn)
            @test length(topos) == 223
            # Test first (random, most forms) and last topology
            for topo in [topos[1], topos[end]]
                result =
                    EnzymeRates._expand_substrate_product_dead_ends(
                        [topo], ter_ter_rxn)
                @test all(isempty(_connectivity_violations(steps))
                          for (steps, _groups) in result)
                # Competition patterns reduce 2^27 to
                # ≤265 variants per topology
                @test length(result) > 0
                @test length(result) <= 265
            end
        end
    end

end

# ═══════════════════════════════════════════════════════════════════════
# Testsets covering non-enumeration features (atom balance from
# @enzyme_reaction; AllostericEnzymeMechanism accessor identity)
# ═══════════════════════════════════════════════════════════════════════

@testset "test reaction atom balance" begin
    for rxn in [pyruvate_carboxylase_rxn,
                pyruvate_dehydrogenase_rxn]
        sub_atoms = Dict{Symbol,Int}()
        prod_atoms = Dict{Symbol,Int}()
        for ra in EnzymeRates.reactants(rxn)
            target = EnzymeRates.metabolite(ra) isa EnzymeRates.Substrate ?
                     sub_atoms : prod_atoms
            for (a, c) in EnzymeRates.atoms(ra)
                target[a] = get(target, a, 0) + c
            end
        end
        @test sub_atoms == prod_atoms
    end
end

@testset "AllostericEnzymeMechanism TR equivalence" begin
    # Steps with kinetic_group 1 = S binding (RE, :EqualAI),
    # 2 = P binding (RE, :OnlyA), 3 = iso (SS, :NonequalAI).
    m_compiled = @allosteric_mechanism begin
        substrates: S
        products: P
        catalytic_multiplicity: 2
        catalytic_steps: begin
            E + S ⇌ E(S)             :: EqualAI
            E + P ⇌ E(P)             :: OnlyA
            E(S) <--> E(P)           :: NonequalAI
        end
    end
    # Group order is canonical; allosteric tags stay bound to their steps.
    am = EnzymeRates.AllostericMechanism(m_compiled)
    state_of(pred) = EnzymeRates.cat_allo_state(am,
        only(g for g in EnzymeRates.kinetic_groups(am)
             if pred(EnzymeRates.bound_metabolite(EnzymeRates.rep_step(am, g)))))
    @test state_of(bm -> bm isa EnzymeRates.Substrate) == :EqualAI
    @test state_of(bm -> bm isa EnzymeRates.Product) == :OnlyA
    @test state_of(bm -> bm === nothing) == :NonequalAI
end

# ═══════════════════════════════════════════════════════════════════════
# 2. Initialization (compile_mechanism + init_mechanisms)
# ═══════════════════════════════════════════════════════════════════════

# ─── compile_mechanism dispatch ────────────────────────────────────────
@testset "compile_mechanism dispatch" begin
    # `compile_mechanism` dispatches a `Mechanism` to the
    # `EnzymeMechanism` constructor and an `AllostericMechanism` to the
    # `AllostericEnzymeMechanism` constructor. Verify both legs.
    m = first(EnzymeRates.init_mechanisms(uni_uni_rxn))
    @test EnzymeRates.compile_mechanism(m) === EnzymeMechanism(m)

    am_seed = @allosteric_mechanism begin
        substrates: S
        products: P
        catalytic_multiplicity: 2
        catalytic_steps: begin
            E + P ⇌ E(P)   :: EqualAI
            E + S ⇌ E(S)   :: EqualAI
            E(S) <--> E(P) :: EqualAI
        end
    end
    am = EnzymeRates.AllostericMechanism(am_seed)
    @test EnzymeRates.compile_mechanism(am) === AllostericEnzymeMechanism(am)
end

# ─── init_mechanisms ───────────────────────────────────────────────────
@testset "init_mechanisms" begin

    @testset "init mechanisms are connected" begin
        for (rxn, n_s, n_p) in [
            (uni_uni_rxn, 1, 1),
            (uni_bi_rxn, 1, 2),
            (bi_bi_rxn, 2, 2),
            (bi_bi_pp_rxn, 2, 2),
        ]
            specs = EnzymeRates.init_mechanisms(rxn)
            @test all(isempty(_connectivity_violations(
                EnzymeRates.steps(m))) for m in specs)
        end
    end

    @testset "exactly 1 SS step per init mechanism" begin
        # init_mechanisms produces minimum-parameter mechanisms — exactly
        # one isomerization step, which is SS by construction. Subsequent
        # RE→SS expansions add more SS steps; init never does.
        for rxn in [uni_uni_rxn, uni_bi_rxn,
                    bi_bi_rxn, bi_bi_pp_rxn]
            specs = EnzymeRates.init_mechanisms(rxn)
            @test all(isempty(_connectivity_violations(
                EnzymeRates.steps(m))) for m in specs)
            for s in specs
                @test count(st -> !st.is_equilibrium,
                            Iterators.flatten(s.steps)) == 1
            end
        end
    end

    @testset "Same-metabolite RE bindings share kinetic_group" begin
        # _apply_equivalence_grouping collapses all RE binding steps for
        # the same metabolite into one kinetic group (one shared K).
        # For bi-bi, metabolites like :B appear in multiple binding steps
        # (e.g. E+B⇌E_B and E_A+B⇌E_A_B) — these must share one kinetic_group.
        specs = EnzymeRates.init_mechanisms(bi_bi_rxn)
        @test all(isempty(_connectivity_violations(
            EnzymeRates.steps(m))) for m in specs)
        @test !isempty(specs)
        n_assertions_fired = 0
        for spec in specs
            # Map each metabolite to the kinetic-group index (inner-vector
            # position) of every RE binding step that binds it. Same-group
            # sharing is structural: same inner vector == same kinetic group.
            by_metabolite = Dict{Symbol, Vector{Int}}()
            for (gi, group) in enumerate(spec.steps)
                for step in group
                    EnzymeRates.is_equilibrium(step) || continue
                    bm = EnzymeRates.bound_metabolite(step)
                    bm === nothing && continue
                    push!(get!(by_metabolite, EnzymeRates.name(bm),
                               Int[]), gi)
                end
            end
            for (_met, gis) in by_metabolite
                length(gis) >= 2 || continue
                @test length(Set(gis)) == 1
                n_assertions_fired += 1
            end
        end
        @test n_assertions_fired >= 1   # at least one multi-binding case existed
    end

    @testset "Uni-uni: exactly 1 init mechanism" begin
        # Uni-uni topology: 1 catalytic topology × 1 dead-end variant
        # (none possible — see test_expand_substrate_product_dead_ends
        # uni-uni case). Hence init produces exactly 1 mechanism.
        specs = EnzymeRates.init_mechanisms(uni_uni_rxn)
        @test all(isempty(_connectivity_violations(
            EnzymeRates.steps(m))) for m in specs)
        @test length(specs) == 1
    end

    @testset "Init compiles for all small reactions" begin
        # Every init mechanism must compile to a valid EnzymeMechanism.
        # Tests first 5 mechanisms per reaction to cap @generated cost.
        for rxn in [uni_uni_rxn, bi_bi_rxn, bi_bi_pp_rxn]
            specs = EnzymeRates.init_mechanisms(rxn)
            @test all(isempty(_connectivity_violations(
                EnzymeRates.steps(m))) for m in specs)
            for spec in first(specs, 5)
                m = EnzymeMechanism(spec)
                @test m isa EnzymeMechanism
            end
        end
    end

    @testset "bi-bi exit gate: init mechanisms derive (subset)" begin
        mechs = EnzymeRates.init_mechanisms(bi_bi_rxn)
        @test all(isempty(_connectivity_violations(
            EnzymeRates.steps(m))) for m in mechs)
        @test length(mechs) == 55
        # Derive a small subset only — full derivation is slow. Pick the 5
        # smallest by step count (cheapest to compile).
        by_size = sort(mechs; by = m -> EnzymeRates.n_steps(m))
        for m in by_size[1:5]
            s = EnzymeRates.rate_equation_string(EnzymeRates.compile_mechanism(m))
            @test s isa AbstractString && !isempty(s)
        end
    end

    @testset "Drops unbound regulators from init Mechanism" begin
        # init_mechanisms produces Mechanisms without dead-end regulators
        # bound. When compiled to EnzymeMechanism, the regulator must NOT
        # appear in the regulators tuple — only the catalytic mechanism is
        # built. After expand_mechanisms adds the dead-end regulator, it
        # should appear.
        init_mechs = EnzymeRates.init_mechanisms(uni_uni_with_reg)
        @test all(isempty(_connectivity_violations(
            EnzymeRates.steps(m))) for m in init_mechs)
        @test !isempty(init_mechs)
        for m in init_mechs
            em = EnzymeRates.compile_mechanism(m)
            @test :I ∉ EnzymeRates.regulators(em)
        end

        expanded = EnzymeRates.expand_mechanisms(init_mechs, uni_uni_with_reg)
        found_with_reg = false
        for mm in expanded
            em = EnzymeRates.compile_mechanism(mm)
            if :I in EnzymeRates.regulators(em)
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
            substrates: L_Ala[CHN]
            products: D_Ala[CHN]
        end
        specs = EnzymeRates.init_mechanisms(rxn)
        @test all(isempty(_connectivity_violations(
            EnzymeRates.steps(m))) for m in specs)
        @test !isempty(specs)
        for spec in first(specs, min(3, length(specs)))
            m = EnzymeMechanism(spec)
            @test m isa EnzymeMechanism
        end
    end

end

# ═══════════════════════════════════════════════════════════════════════
# 3. Expansion moves
# ═══════════════════════════════════════════════════════════════════════

# ─── _expand_re_to_ss ──────────────────────────────────────────────────
@testset "_expand_re_to_ss" begin

    @testset "Mechanism — bi-bi sequential: 4 RE binding groups → 4 variants" begin
        # SEED: bi-bi sequential ordered, 4 singleton RE binding groups + 1
        # SS iso. _expand_re_to_ss fires per RE group → 4 variants.
        em_seed = @enzyme_mechanism begin
            substrates: A, B
            products: P, Q
            steps: begin
                E + A ⇌ E(A)
                E(A) + B ⇌ E(A, B)
                E + Q ⇌ E(Q)
                E(Q) + P ⇌ E(P, Q)
                E(A, B) <--> E(P, Q)
            end
        end
        m = EnzymeRates.Mechanism(em_seed)
        EnzymeRates._assert_mechanism_invariants(m)

        result = EnzymeRates._expand_re_to_ss(m)

        # 1. count: 4 all-RE singleton groups (A, B, Q, P bindings). Iso SS. → 4.
        @test length(result) == 4
        for r in result
            @test r isa EnzymeRates.Mechanism
            EnzymeRates._assert_mechanism_invariants(r)
            @test EnzymeRates.reaction(r) == EnzymeRates.reaction(m)
            @test length(r.steps) == length(m.steps)
        end

        # 2. property-style: each variant flips exactly one initial RE
        # group to all-SS; all other groups unchanged.
        for r in result
            n_groups_newly_ss = count(zip(m.steps, r.steps)) do (old_grp, new_grp)
                length(old_grp) == length(new_grp) &&
                    all(EnzymeRates.is_equilibrium, old_grp) &&
                    !any(EnzymeRates.is_equilibrium, new_grp)
            end
            @test n_groups_newly_ss == 1
        end

        # 3. distinct flipped group across variants (one per RE group).
        flipped_groups = Int[]
        for r in result
            for (gi, (old_grp, new_grp)) in enumerate(zip(m.steps, r.steps))
                if all(EnzymeRates.is_equilibrium, old_grp) &&
                   !any(EnzymeRates.is_equilibrium, new_grp)
                    push!(flipped_groups, gi)
                end
            end
        end
        @test length(unique(flipped_groups)) == 4
    end

    @testset "Mechanism — bi-bi multi-step kinetic group: atomic conversion" begin
        # SEED: bi-bi random with 4 multi-step RE groups (A, B, P, Q each
        # binding at two forms, shared via parens) + 1 SS iso. RE→SS fires
        # per group atomically — both steps in a multi-step group flip together.
        em_seed = @enzyme_mechanism begin
            substrates: A, B
            products: P, Q
            steps: begin
                (E + A ⇌ E(A), E(B) + A ⇌ E(A, B))
                (E + B ⇌ E(B), E(A) + B ⇌ E(A, B))
                (E + P ⇌ E(P), E(Q) + P ⇌ E(P, Q))
                (E + Q ⇌ E(Q), E(P) + Q ⇌ E(P, Q))
                E(A, B) <--> E(P, Q)
            end
        end
        m = EnzymeRates.Mechanism(em_seed)
        EnzymeRates._assert_mechanism_invariants(m)
        # Sanity: the seed has 4 multi-step groups + 1 singleton iso.
        @test count(g -> length(g) >= 2, m.steps) == 4
        @test count(g -> length(g) == 1, m.steps) == 1

        result = EnzymeRates._expand_re_to_ss(m)

        # 1. count: 4 all-RE multi-step groups → 4 variants. Iso SS excluded.
        @test length(result) == 4
        for r in result
            @test r isa EnzymeRates.Mechanism
            EnzymeRates._assert_mechanism_invariants(r)
        end

        # 2. property-style: in each variant, exactly one initial RE
        # multi-step group has ALL its steps now SS (atomic conversion).
        for r in result
            n_all_ss_multi = count(zip(m.steps, r.steps)) do (old_grp, new_grp)
                length(old_grp) >= 2 &&
                    all(EnzymeRates.is_equilibrium, old_grp) &&
                    !any(EnzymeRates.is_equilibrium, new_grp)
            end
            @test n_all_ss_multi == 1
        end

        # 3. preservation: reaction unchanged on every variant.
        for r in result
            @test EnzymeRates.reaction(r) == EnzymeRates.reaction(m)
        end
    end

    @testset "Mechanism — bi-bi ping-pong: 5 RE groups → 5 variants" begin
        # SEED: bi-bi ping-pong with Estar (residual) form. 5 singleton RE
        # groups + 1 SS iso group. _expand_re_to_ss fires per RE group → 5.
        em_seed = @enzyme_mechanism begin
            substrates: A, B
            products: P, Q
            steps: begin
                E + A ⇌ E(A)
                Estar + B ⇌ Estar(B)
                E + Q ⇌ E(Q)
                Estar + P ⇌ Estar(A, P)
                E(A) <--> Estar(A, P)
                Estar(B) ⇌ E(Q)
            end
        end
        m = EnzymeRates.Mechanism(em_seed)
        EnzymeRates._assert_mechanism_invariants(m)

        result = EnzymeRates._expand_re_to_ss(m)

        # 1. count: 5 all-RE groups → 5 variants.
        @test length(result) == 5
        for r in result
            @test r isa EnzymeRates.Mechanism
            EnzymeRates._assert_mechanism_invariants(r)
        end

        # 2. property-style: each variant flips exactly one initial RE
        # group to all-SS. The flipped group is distinct across variants.
        flipped_groups = Int[]
        for r in result
            for (gi, (old_grp, new_grp)) in enumerate(zip(m.steps, r.steps))
                if all(EnzymeRates.is_equilibrium, old_grp) &&
                   !any(EnzymeRates.is_equilibrium, new_grp)
                    push!(flipped_groups, gi)
                end
            end
        end
        @test length(unique(flipped_groups)) == 5

        # 3. preservation
        for r in result
            @test EnzymeRates.reaction(r) == EnzymeRates.reaction(m)
        end
    end

    @testset "Mechanism — ter-ter sequential" begin
        # SEED: ter-ter sequential ordered. 6 singleton RE binding groups
        # + 1 SS iso. _expand_re_to_ss fires per RE group → 6 variants.
        em_seed = @enzyme_mechanism begin
            substrates: A, B, D
            products: P, Q, R
            steps: begin
                E + A ⇌ E(A)
                E(A) + B ⇌ E(A, B)
                E(A, B) + D ⇌ E(A, B, D)
                E + R ⇌ E(R)
                E(R) + Q ⇌ E(Q, R)
                E(Q, R) + P ⇌ E(P, Q, R)
                E(A, B, D) <--> E(P, Q, R)
            end
        end
        m = EnzymeRates.Mechanism(em_seed)
        EnzymeRates._assert_mechanism_invariants(m)

        result = EnzymeRates._expand_re_to_ss(m)
        @test length(result) == 6
        for r in result
            @test r isa EnzymeRates.Mechanism
            EnzymeRates._assert_mechanism_invariants(r)
            @test EnzymeRates.compile_mechanism(r) isa EnzymeMechanism
        end
    end

    @testset "Mechanism — :EqualAI group: RE→SS preserves tag" begin
        # SEED: uni-uni allosteric with all catalytic groups :EqualAI.
        # _expand_re_to_ss fires per all-RE group → 2 variants. Each
        # converted group stays :EqualAI (RE K → SS (kf, kr), shared R/T).
        em_seed = @allosteric_mechanism begin
            substrates: S
            products: P
            catalytic_multiplicity: 2
            catalytic_steps: begin
                E + P ⇌ E(P)      :: EqualAI
                E + S ⇌ E(S)      :: EqualAI
                E(S) <--> E(P)    :: EqualAI
            end
        end
        am = EnzymeRates.AllostericMechanism(em_seed)
        EnzymeRates._assert_mechanism_invariants(am)

        result = EnzymeRates._expand_re_to_ss(am)

        # 1. count: 2 all-RE groups; iso SS. → 2.
        @test length(result) == 2

        # 2. Δ params: cheap-tag RE→SS adds 1 fitted param per variant
        # (RE K → SS (kf, kr), shared across R/T — no separate T-state pair).
        # Measured against the actual compiled fitted-param count.
        base_fitted = length(EnzymeRates.fitted_params(
            EnzymeRates.compile_mechanism(am)))
        deltas = sort([length(EnzymeRates.fitted_params(
            EnzymeRates.compile_mechanism(r))) - base_fitted for r in result])
        @test deltas == [1, 1]

        # 3. compilability — must produce AllostericEnzymeMechanism.
        for r in result
            @test r isa EnzymeRates.AllostericMechanism
            EnzymeRates._assert_mechanism_invariants(r)
            @test EnzymeRates.compile_mechanism(r) isa AllostericEnzymeMechanism
        end

        # 4. exactly one group newly all-SS; all cat_allo_states preserved.
        for r in result
            n_newly_ss = count(zip(am.cat_steps, r.cat_steps)) do (old, new)
                all(EnzymeRates.is_equilibrium, old) &&
                    !any(EnzymeRates.is_equilibrium, new)
            end
            @test n_newly_ss == 1
            @test r.cat_allo_states == am.cat_allo_states
        end

        # 5. preservation: multiplicity, reg sites, reaction untouched.
        for r in result
            @test r.catalytic_multiplicity == am.catalytic_multiplicity
            @test r.regulatory_sites == am.regulatory_sites
            @test EnzymeRates.reaction(r) == EnzymeRates.reaction(am)
        end
    end

    @testset "Mechanism — :OnlyA group: RE→SS preserves tag" begin
        # SEED: uni-uni allosteric with one :OnlyA group (S-binding,
        # group 2), others :EqualAI. :OnlyA groups live in the R-state
        # only; T-state contributes no kf_T/kr_T after RE→SS.
        em_seed = @allosteric_mechanism begin
            substrates: S
            products: P
            catalytic_multiplicity: 2
            catalytic_steps: begin
                E + P ⇌ E(P)      :: EqualAI
                E + S ⇌ E(S)      :: OnlyA
                E(S) <--> E(P)    :: EqualAI
            end
        end
        am = EnzymeRates.AllostericMechanism(em_seed)
        EnzymeRates._assert_mechanism_invariants(am)

        result = EnzymeRates._expand_re_to_ss(am)

        # 1. count: 2 all-RE groups (P-binding :EqualAI, S-binding :OnlyA).
        # Iso SS. → 2 variants.
        @test length(result) == 2

        # 2. Δ params. Connected-component pruning populates the inactive state
        # from free E over ALL surviving steps, so because catalysis is :EqualAI
        # the base mechanism's Q_I ALREADY carries E(S) via REVERSE catalysis
        # E(P)→E(S) — even though S-binding is :OnlyA and cannot bind directly in
        # the T-state. E(S) is therefore present in the inactive state of the
        # base AND of every variant, so flipping a binding group RE→SS adds only
        # that group's own SS rate constant (+1) — it does not introduce a new
        # inactive-state form. Both variants → Δ=1. Equilibrium flux stays 0 for
        # both (verified at Keq mass-action ratio, |v_eq|/|v| ~ 1e-16).
        base_fitted = length(EnzymeRates.fitted_params(
            EnzymeRates.compile_mechanism(am)))
        deltas = sort([length(EnzymeRates.fitted_params(
            EnzymeRates.compile_mechanism(r))) - base_fitted for r in result])
        @test deltas == [1, 1]

        # 3. compilability
        for r in result
            @test r isa EnzymeRates.AllostericMechanism
            EnzymeRates._assert_mechanism_invariants(r)
            @test EnzymeRates.compile_mechanism(r) isa AllostericEnzymeMechanism
        end

        # 4. exactly one group flipped to SS; ALL cat_allo_states preserved
        # (including the converted group's :OnlyA — move MUST NOT change
        # R/T-state semantics).
        for r in result
            n_newly_ss = count(zip(am.cat_steps, r.cat_steps)) do (old, new)
                all(EnzymeRates.is_equilibrium, old) &&
                    !any(EnzymeRates.is_equilibrium, new)
            end
            @test n_newly_ss == 1
            @test r.cat_allo_states == am.cat_allo_states
        end

        # 5. preservation
        for r in result
            @test r.catalytic_multiplicity == am.catalytic_multiplicity
            @test r.regulatory_sites == am.regulatory_sites
            @test EnzymeRates.reaction(r) == EnzymeRates.reaction(am)
        end
    end

    @testset "Mechanism — :NonequalAI group: RE→SS adds 2 params" begin
        # SEED: uni-uni allosteric with one :NonequalAI group (S-binding),
        # others :EqualAI. When RE→SS converts a :NonequalAI group, BOTH
        # the R-state K and the T-state K_T must split into (kf, kr) and
        # (kf_T, kr_T). Δ for :NonequalAI = +2; cheap-tag conversion = +1.
        em_seed = @allosteric_mechanism begin
            substrates: S
            products: P
            catalytic_multiplicity: 2
            catalytic_steps: begin
                E + P ⇌ E(P)      :: EqualAI
                E + S ⇌ E(S)      :: NonequalAI
                E(S) <--> E(P)    :: EqualAI
            end
        end
        am = EnzymeRates.AllostericMechanism(em_seed)
        EnzymeRates._assert_mechanism_invariants(am)

        result = EnzymeRates._expand_re_to_ss(am)

        # 1. count: 2 RE groups (P-binding :EqualAI, S-binding :NonequalAI).
        # → 2 variants.
        @test length(result) == 2

        # 2. Δ params: P-binding :EqualAI → +1; S-binding :NonequalAI → +2.
        # Measured against the actual compiled fitted-param count (ground
        # truth): SS :NonequalAI adds kf_T, kr_T on top of the R-state pair.
        base_fitted = length(EnzymeRates.fitted_params(
            EnzymeRates.compile_mechanism(am)))
        deltas = sort([length(EnzymeRates.fitted_params(
            EnzymeRates.compile_mechanism(r))) - base_fitted for r in result])
        @test deltas == [1, 2]

        # 3. compilability
        for r in result
            @test r isa EnzymeRates.AllostericMechanism
            EnzymeRates._assert_mechanism_invariants(r)
            @test EnzymeRates.compile_mechanism(r) isa AllostericEnzymeMechanism
        end

        # 4. exactly one group flipped to SS; cat_allo_states preserved.
        for r in result
            n_newly_ss = count(zip(am.cat_steps, r.cat_steps)) do (old, new)
                all(EnzymeRates.is_equilibrium, old) &&
                    !any(EnzymeRates.is_equilibrium, new)
            end
            @test n_newly_ss == 1
            @test r.cat_allo_states == am.cat_allo_states
        end

        # 5. preservation
        for r in result
            @test r.catalytic_multiplicity == am.catalytic_multiplicity
            @test r.regulatory_sites == am.regulatory_sites
            @test EnzymeRates.reaction(r) == EnzymeRates.reaction(am)
        end
    end

    @testset "Mechanism — Substrate-as-dead-end-inhibitor overlap" begin
        # SEED: uni-uni where S is BOTH substrate AND dead-end inhibitor.
        # Build via init + dead-end expansion to get the S/__reg overlap.
        # _expand_re_to_ss should treat the substrate-S and inhibitor-S
        # kinetic groups as independent.
        rxn = @enzyme_reaction begin
            substrates: S[C]
            products: P[C]
            dead_end_inhibitors: S
        end
        init_ms = EnzymeRates.init_mechanisms(rxn)
        @test length(init_ms) == 1   # uni-uni: 1 catalytic topology
        seed = first(init_ms)
        de_ms = EnzymeRates._expand_add_dead_end_regulator(seed, rxn)
        @test !isempty(de_ms)
        m = first(de_ms)
        EnzymeRates._assert_mechanism_invariants(m)

        result = EnzymeRates._expand_re_to_ss(m)

        # 1. count: after dead-end expansion, groups are substrate-binding (RE),
        # product-binding (RE), iso (SS), dead-end-S__reg-binding (RE).
        # → 3 RE groups → 3 variants.
        @test length(result) == 3
        for r in result
            @test r isa EnzymeRates.Mechanism
            EnzymeRates._assert_mechanism_invariants(r)
            @test EnzymeRates.compile_mechanism(r) isa EnzymeMechanism
        end

        # 2. property-style: each variant flips exactly one initial RE
        # group to all-SS, and the flipped group covers 3 distinct values
        # across the 3 variants. This proves the move treats the substrate
        # and inhibitor kinetic groups as independent even when they share
        # the metabolite name :S.
        flipped_groups = Int[]
        for r in result
            for (gi, (old_grp, new_grp)) in enumerate(zip(m.steps, r.steps))
                if all(EnzymeRates.is_equilibrium, old_grp) &&
                   !any(EnzymeRates.is_equilibrium, new_grp)
                    push!(flipped_groups, gi)
                end
            end
        end
        @test length(unique(flipped_groups)) == 3

        # 3. preservation: reaction unchanged.
        for r in result
            @test EnzymeRates.reaction(r) == EnzymeRates.reaction(m)
        end
    end

    @testset "Mechanism — Allosteric substrate-as-dead-end-I overlap" begin
        # AllostericMechanism counterpart of the substrate-as-I overlap.
        # Same count derivation: 3 RE groups → 3 variants. Each flip should
        # preserve cat_allo_states tags exactly.
        rxn = @enzyme_reaction begin
            substrates: S[C]
            products: P[C]
            dead_end_inhibitors: S
            oligomeric_state: 2
        end
        init_ms = EnzymeRates.init_mechanisms(rxn)
        seed = first(init_ms)
        de_ms = EnzymeRates._expand_add_dead_end_regulator(seed, rxn)
        @test !isempty(de_ms)
        plain = first(de_ms)
        # Convert to allosteric
        allo_ms = EnzymeRates._expand_to_allosteric(plain, rxn)
        @test !isempty(allo_ms)
        am = first(allo_ms)
        EnzymeRates._assert_mechanism_invariants(am)

        result = EnzymeRates._expand_re_to_ss(am)

        # 1. count: 3 RE groups → 3 variants.
        @test length(result) == 3
        for r in result
            @test r isa EnzymeRates.AllostericMechanism
            EnzymeRates._assert_mechanism_invariants(r)
            @test EnzymeRates.compile_mechanism(r) isa AllostericEnzymeMechanism
        end

        # 2. cat_allo_states preserved on every variant.
        for r in result
            @test r.cat_allo_states == am.cat_allo_states
            @test r.catalytic_multiplicity == am.catalytic_multiplicity
            @test r.regulatory_sites == am.regulatory_sites
            @test EnzymeRates.reaction(r) == EnzymeRates.reaction(am)
        end
    end

    @testset "Mechanism — uni-uni: 2 RE binding groups → 2 variants" begin
        # SEED: uni-uni init mechanism. 3 kinetic groups: 2 RE binding
        # (S-binding, P-binding) and 1 SS iso. _expand_re_to_ss fires per
        # all-RE group → 2 variants.
        m = first(EnzymeRates.init_mechanisms(uni_uni_rxn))
        @test m isa EnzymeRates.Mechanism

        result = EnzymeRates._expand_re_to_ss(m)

        # 1. count: 2 all-RE groups (P-binding, S-binding). Iso already SS. → 2.
        @test length(result) == 2
        for r in result
            @test r isa EnzymeRates.Mechanism
            @test EnzymeRates.reaction(r) == EnzymeRates.reaction(m)
            @test length(r.steps) == length(m.steps)
        end

        # 2. property-style: in each variant, exactly one initial RE
        # group has all its steps newly SS; all other groups unchanged.
        for r in result
            n_groups_newly_ss = count(zip(m.steps, r.steps)) do (old_grp, new_grp)
                length(old_grp) == length(new_grp) &&
                    all(EnzymeRates.is_equilibrium, old_grp) &&
                    !any(EnzymeRates.is_equilibrium, new_grp)
            end
            @test n_groups_newly_ss == 1
        end

        # 3. step structure preserved: each Step's chemistry (from/to
        # species + bound metabolite) in the result matches the
        # corresponding Step in the seed, position-for-position.
        for r in result
            for (old_grp, new_grp) in zip(m.steps, r.steps)
                @test [EnzymeRates.from_species(s) for s in old_grp] ==
                      [EnzymeRates.from_species(s) for s in new_grp]
                @test [EnzymeRates.to_species(s) for s in old_grp] ==
                      [EnzymeRates.to_species(s) for s in new_grp]
                @test [EnzymeRates.bound_metabolite(s) for s in old_grp] ==
                      [EnzymeRates.bound_metabolite(s) for s in new_grp]
            end
        end
    end

    @testset "Mechanism — all-SS seed: empty (negative)" begin
        # If every catalytic group is already SS, the move has no RE
        # group to fire on → empty.
        m_seed_with_all_ss = first(EnzymeRates.init_mechanisms(uni_uni_rxn))
        # Convert every group to SS by applying _flip_group_to_ss
        # repeatedly. (Use the same helper the production code uses.)
        all_ss_groups = m_seed_with_all_ss.steps
        for g in 1:length(all_ss_groups)
            all_ss_groups = EnzymeRates._flip_group_to_ss(all_ss_groups, g)
        end
        m_all_ss = EnzymeRates.Mechanism(
            EnzymeRates.reaction(m_seed_with_all_ss), all_ss_groups)
        @test isempty(EnzymeRates._expand_re_to_ss(m_all_ss))
    end

    @testset "AllostericMechanism — :EqualAI: 2 variants, tags preserved" begin
        # SEED: uni-uni allosteric init via Mechanism path → promote to
        # AllostericMechanism. Default tag is :EqualAI post init (since
        # init mechanisms are plain → :EqualAI after _expand_to_allosteric).
        init_mechs = EnzymeRates.init_mechanisms(uni_uni_allo)
        m_seed = first(init_mechs)
        allo_mechs = EnzymeRates._expand_to_allosteric(m_seed, uni_uni_allo)
        @test !isempty(allo_mechs)
        am = first(allo_mechs)

        result = EnzymeRates._expand_re_to_ss(am)

        # 1. count: 2 all-RE binding groups (P-, S-binding). Iso is SS. → 2.
        @test length(result) == 2

        # 2. Δ params: :EqualAI RE→SS adds 1 fitted param per variant
        # (RE K → SS (kf, kr) with K/Wegscheider absorbing one). Measured
        # against the actual compiled count.
        base_fitted = length(EnzymeRates.fitted_params(
            EnzymeRates.compile_mechanism(am)))
        deltas = sort([length(EnzymeRates.fitted_params(
            EnzymeRates.compile_mechanism(r))) - base_fitted for r in result])
        @test deltas == [1, 1]

        # 3. allosteric state preserved on every variant.
        for r in result
            @test r isa EnzymeRates.AllostericMechanism
            @test r.cat_allo_states == am.cat_allo_states
            @test r.catalytic_multiplicity == am.catalytic_multiplicity
            @test r.regulatory_sites == am.regulatory_sites
            @test EnzymeRates.reaction(r) == EnzymeRates.reaction(am)
        end

        # 4. exactly one group newly all-SS.
        for r in result
            n_newly_ss = count(zip(am.cat_steps, r.cat_steps)) do (old, new)
                all(EnzymeRates.is_equilibrium, old) &&
                    !any(EnzymeRates.is_equilibrium, new)
            end
            @test n_newly_ss == 1
        end
    end

end

# ─── _expand_split_kinetic_group ───────────────────────────────────────
@testset "_expand_split_kinetic_group" begin

    @testset "Mechanism — mixed RE/SS multi-step groups split per member" begin
        # SEED: bi-bi random where the A-binding kinetic group is SS
        # (size-2, both SS) and the B-binding kinetic group is RE
        # (size-2, both RE). The remaining P-binding (size-2 RE),
        # Q-binding (size-2 RE), and iso (singleton SS) groups are
        # unchanged. Total: 9 steps, 5 kinetic groups. Each multi-step
        # group splits per member, peeling one step into a new group.
        m_seed = @enzyme_mechanism begin
            substrates: A, B
            products: P, Q
            steps: begin
                (E + A <--> E(A), E(B) + A <--> E(A, B))
                (E + B ⇌ E(B), E(A) + B ⇌ E(A, B))
                (E + P ⇌ E(P), E(Q) + P ⇌ E(P, Q))
                (E + Q ⇌ E(Q), E(P) + Q ⇌ E(P, Q))
                E(A, B) <--> E(P, Q)
            end
        end
        m = EnzymeRates.Mechanism(m_seed)
        EnzymeRates._assert_mechanism_invariants(m)

        result = EnzymeRates._expand_split_kinetic_group(m)

        # 1. count: 4 multi-step groups (A SS×2, B RE×2, P RE×2, Q RE×2),
        # 4 × 2 = 8 candidates. Each candidate is canonicalized and dropped
        # if it collapses back to the parent (a Wegscheider-tied binding-K
        # rename — a model-space no-op). The A and B splits survive; the P
        # and Q splits are self-loops and are dropped, leaving 4 variants.
        @test length(result) == 4

        # 2. Δ params measured against the actual compiled fitted count.
        # All 4 surviving splits add one parameter (Δ=1); the dropped P/Q
        # splits were the Δ=0 self-loops that canonicalization removes.
        base_fitted = length(EnzymeRates.fitted_params(
            EnzymeRates.compile_mechanism(m)))
        deltas = sort([length(EnzymeRates.fitted_params(
            EnzymeRates.compile_mechanism(r))) - base_fitted for r in result])
        @test deltas == [1, 1, 1, 1]

        # 3. compilability
        for r in result
            @test r isa EnzymeRates.Mechanism
            EnzymeRates._assert_mechanism_invariants(r)
            @test EnzymeRates.compile_mechanism(r) isa EnzymeMechanism
        end

        # 4. property-style: each result introduces exactly one new
        # trailing kinetic group with exactly one step in it.
        for r in result
            @test length(r.steps) == length(m.steps) + 1
            @test length(last(r.steps)) == 1
        end

        # 5. total step count preserved across the split.
        for r in result
            @test EnzymeRates.n_steps(r) == EnzymeRates.n_steps(m)
        end

        # 6. preservation
        for r in result
            @test EnzymeRates.reaction(r) == EnzymeRates.reaction(m)
        end
    end

    @testset "AllostericMechanism — SS multi-step :NonequalAI split" begin
        # SEED: bi-bi allosteric where one multi-step group is BOTH SS AND
        # :NonequalAI. Splitting this group costs more parameters than
        # splitting a :EqualAI RE group (factor 2 for SS pair × factor 2
        # for R/T-state pair).
        em_seed = @allosteric_mechanism begin
            substrates: A, B
            products: P, Q
            catalytic_multiplicity: 2
            catalytic_steps: begin
                (E + A <--> E(A), E(B) + A <--> E(A, B))    :: NonequalAI
                (E + B ⇌ E(B), E(A) + B ⇌ E(A, B))          :: EqualAI
                (E + P ⇌ E(P), E(Q) + P ⇌ E(P, Q))          :: EqualAI
                (E + Q ⇌ E(Q), E(P) + Q ⇌ E(P, Q))          :: EqualAI
                E(A, B) <--> E(P, Q)                        :: EqualAI
            end
        end
        am = EnzymeRates.AllostericMechanism(em_seed)
        EnzymeRates._assert_mechanism_invariants(am)

        result = EnzymeRates._expand_split_kinetic_group(am)

        # 1. count: 4 multi-step groups (A-binding SS×2 :NonequalAI,
        # B-binding RE×2 :EqualAI, P-binding RE×2 :EqualAI,
        # Q-binding RE×2 :EqualAI), 4 × 2 members = 8 candidates. Each
        # candidate is canonicalized and dropped if it collapses back to
        # the parent (a Wegscheider-tied self-loop). 4 of the 6 RE splits
        # are such self-loops; the 2 SS splits never tie (no equilibrium
        # constant to absorb them), so 4 variants survive.
        @test length(result) == 4

        # 2. Δ params: 4 surviving variants, deltas measured against the
        # actual compiled fitted-param count (ground truth — true count
        # after thermo-cycle bookkeeping). 2 add 1 (the surviving RE
        # splits), 2 add 2 (the :NonequalAI SS splits, doubled by the
        # R/T-state pair).
        base_fitted = length(EnzymeRates.fitted_params(
            EnzymeRates.compile_mechanism(am)))
        deltas = sort([length(EnzymeRates.fitted_params(
            EnzymeRates.compile_mechanism(r))) - base_fitted for r in result])
        @test deltas == [1, 1, 2, 2]

        # 3. compilability
        for r in result
            @test r isa EnzymeRates.AllostericMechanism
            EnzymeRates._assert_mechanism_invariants(r)
            @test EnzymeRates.compile_mechanism(r) isa AllostericEnzymeMechanism
        end

        # 4. tag inheritance: each result group inherits the tag of the parent
        # group whose steps contain it. A split subdivides one group into two;
        # both halves carry that group's tag, and every other group keeps its
        # own. Group ORDER is canonical (not source-preserved), so match each
        # result group to its parent by step content rather than by position.
        for r in result
            @test length(r.cat_allo_states) == length(am.cat_allo_states) + 1
            for (g, grp) in enumerate(r.cat_steps)
                sset = Set(grp)
                parent = only(ag for ag in 1:length(am.cat_steps)
                              if sset ⊆ Set(am.cat_steps[ag]))
                @test r.cat_allo_states[g] == am.cat_allo_states[parent]
            end
        end

        # 5. preservation
        for r in result
            @test r.catalytic_multiplicity == am.catalytic_multiplicity
            @test r.regulatory_sites == am.regulatory_sites
            @test EnzymeRates.reaction(r) == EnzymeRates.reaction(am)
        end
    end

    @testset "Mechanism — bi-bi random: all splits Wegscheider-tied (negative)" begin
        # SEED: bi-bi random with 4 multi-step kinetic groups (A, B, P, Q)
        # forming a single closed thermodynamic cycle. Splitting off any one
        # member of any group produces a binding-K that is a single-symbol
        # Wegscheider rename of an existing parameter, so every one of the
        # 4 × 2 = 8 candidates canonicalizes back to the parent and is
        # dropped as a self-loop.
        m_seed = @enzyme_mechanism begin
            substrates: A, B
            products: P, Q
            steps: begin
                (E + A ⇌ E(A), E(B) + A ⇌ E(A, B))
                (E + B ⇌ E(B), E(A) + B ⇌ E(A, B))
                (E + P ⇌ E(P), E(Q) + P ⇌ E(P, Q))
                (E + Q ⇌ E(Q), E(P) + Q ⇌ E(P, Q))
                E(A, B) <--> E(P, Q)
            end
        end
        m = EnzymeRates.Mechanism(m_seed)
        @test isempty(EnzymeRates._expand_split_kinetic_group(m))
    end

    @testset "Mechanism — all singleton groups: empty (negative)" begin
        # If every group is a singleton, no split is possible.
        m_seed = @enzyme_mechanism begin
            substrates: S
            products: P
            steps: begin
                E + P ⇌ E(P)
                E + S ⇌ E(S)
                E(S) <--> E(P)
            end
        end
        m = EnzymeRates.Mechanism(m_seed)
        @test isempty(EnzymeRates._expand_split_kinetic_group(m))
    end

    @testset "AllostericMechanism — bi-bi RE groups all Wegscheider-tied (negative)" begin
        # SEED: bi-bi allosteric with mixed tags (:NonequalAI, :EqualAI),
        # but both multi-step groups (A, B) are RE bindings, each closing
        # its own per-conformer thermodynamic cycle. Every one of the
        # 2 × 2 = 4 candidates renames to an existing parameter and
        # canonicalizes back to the parent. Tag inheritance on a surviving
        # (non-self-loop) split is covered by "AllostericMechanism — SS
        # multi-step :NonequalAI split" above.
        m_seed = @allosteric_mechanism begin
            substrates: A, B
            products: P, Q
            catalytic_multiplicity: 2
            catalytic_steps: begin
                (E + A ⇌ E(A), E(B) + A ⇌ E(A, B))    :: NonequalAI
                (E + B ⇌ E(B), E(A) + B ⇌ E(A, B))    :: EqualAI
                E + P ⇌ E(P)             :: EqualAI
                E(P) + Q ⇌ E(P, Q)       :: EqualAI
                E + Q ⇌ E(Q)             :: EqualAI
                E(Q) + P ⇌ E(P, Q)       :: EqualAI
                E(A, B) <--> E(P, Q)     :: EqualAI
            end
        end
        am = EnzymeRates.AllostericMechanism(m_seed)
        @test isempty(EnzymeRates._expand_split_kinetic_group(am))
    end
end

# ─── _expand_add_dead_end_regulator ────────────────────────────────────
@testset "_expand_add_dead_end_regulator" begin

    @testset "Mechanism — Sequential bi-bi + I: 4 variants" begin
        # SEED: bi-bi sequential. The expansion should produce 4 form sets
        # after dedup.
        em_seed = @enzyme_mechanism begin
            substrates: A, B
            products: P, Q
            steps: begin
                E + A ⇌ E(A)
                E(A) + B ⇌ E(A, B)
                E + Q ⇌ E(Q)
                E(Q) + P ⇌ E(P, Q)
                E(A, B) <--> E(P, Q)
            end
        end
        m = EnzymeRates.Mechanism(em_seed)
        EnzymeRates._assert_mechanism_invariants(m)
        rxn = @enzyme_reaction begin
            substrates: A[C], B[N]
            products: P[C], Q[N]
            dead_end_inhibitors: I
        end

        result = EnzymeRates._expand_add_dead_end_regulator(m, rxn)

        # 1. count
        @test length(result) == 4

        # 2. Δ params: +1 each (one new K_I parameter), measured against
        # ground-truth `fitted_params(compile_mechanism(...))` — the
        # canonical source for exact parameter counts.
        base_fitted = length(EnzymeRates.fitted_params(
            EnzymeRates.compile_mechanism(m)))
        for r in result
            r_fitted = length(EnzymeRates.fitted_params(
                EnzymeRates.compile_mechanism(r)))
            @test r_fitted == base_fitted + 1
        end

        # 3. compilability
        for r in result
            @test r isa EnzymeRates.Mechanism
            EnzymeRates._assert_mechanism_invariants(r)
            @test EnzymeRates.EnzymeMechanism(r) isa EnzymeMechanism
        end

        # 4. property: each variant has ≥1 I-binding step, and all
        # I-binding steps in a single variant live in the same kinetic
        # group (outer-vector index) — one K_I, not multiple.
        for r in result
            i_groups = Int[]
            for (gi, group) in enumerate(r.steps), s in group
                bm = EnzymeRates.bound_metabolite(s)
                bm !== nothing && EnzymeRates.name(bm) === :I &&
                    push!(i_groups, gi)
            end
            @test !isempty(i_groups)
            @test length(unique(i_groups)) == 1
        end
    end

    @testset "Mechanism — Bi-bi random + I: 9 variants" begin
        # Bi-bi random has 5 eligible forms (E, E_A, E_B, E_P, E_Q);
        # competition patterns × dedup → 9 variants.
        em_seed = @enzyme_mechanism begin
            substrates: A, B
            products: P, Q
            steps: begin
                (E + A ⇌ E(A), E(B) + A ⇌ E(A, B))
                (E + B ⇌ E(B), E(A) + B ⇌ E(A, B))
                (E + P ⇌ E(P), E(Q) + P ⇌ E(P, Q))
                (E + Q ⇌ E(Q), E(P) + Q ⇌ E(P, Q))
                E(A, B) <--> E(P, Q)
            end
        end
        m = EnzymeRates.Mechanism(em_seed)
        EnzymeRates._assert_mechanism_invariants(m)
        rxn = @enzyme_reaction begin
            substrates: A[C], B[N]
            products: P[C], Q[N]
            dead_end_inhibitors: I
        end

        result = EnzymeRates._expand_add_dead_end_regulator(m, rxn)

        # 1. count
        @test length(result) == 9

        # 2. Δ params: +1 each (one new K_I parameter), measured against
        # ground-truth `fitted_params(compile_mechanism(...))` — exact
        # parameter counts come from the compiled mechanism.
        base_fitted = length(EnzymeRates.fitted_params(
            EnzymeRates.compile_mechanism(m)))
        for r in result
            r_fitted = length(EnzymeRates.fitted_params(
                EnzymeRates.compile_mechanism(r)))
            @test r_fitted == base_fitted + 1
        end

        # 3. compilability
        for r in result
            @test r isa EnzymeRates.Mechanism
            EnzymeRates._assert_mechanism_invariants(r)
            @test EnzymeRates.EnzymeMechanism(r) isa EnzymeMechanism
        end

        # 4. property: every variant has ≥1 I-binding step; all I-binding
        # steps in a single variant share the same kinetic group.
        for r in result
            i_groups = Int[]
            for (gi, group) in enumerate(r.steps), s in group
                bm = EnzymeRates.bound_metabolite(s)
                bm !== nothing && EnzymeRates.name(bm) === :I &&
                    push!(i_groups, gi)
            end
            @test !isempty(i_groups)
            @test length(unique(i_groups)) == 1
        end
    end

    @testset "Mechanism — Bi-bi PP + I: 3 variants" begin
        # Ping-pong with Estar; 4 eligible forms (E, E_A, Estar, E_Q);
        # competition patterns × dedup → 3 variants.
        em_seed = @enzyme_mechanism begin
            substrates: A, B
            products: P, Q
            steps: begin
                E + A ⇌ E(A)
                Estar + B ⇌ Estar(B)
                E + Q ⇌ E(Q)
                Estar + P ⇌ Estar(A, P)
                E(A) <--> Estar(A, P)
                Estar(B) ⇌ E(Q)
            end
        end
        m = EnzymeRates.Mechanism(em_seed)
        EnzymeRates._assert_mechanism_invariants(m)
        rxn = @enzyme_reaction begin
            substrates: A[CX], B[N]
            products: P[C], Q[NX]
            dead_end_inhibitors: I
        end

        result = EnzymeRates._expand_add_dead_end_regulator(m, rxn)

        # 1. count
        @test length(result) == 3

        # 2. Δ params: +1 each, measured against ground-truth
        # `fitted_params(compile_mechanism(...))`.
        base_fitted = length(EnzymeRates.fitted_params(
            EnzymeRates.compile_mechanism(m)))
        for r in result
            r_fitted = length(EnzymeRates.fitted_params(
                EnzymeRates.compile_mechanism(r)))
            @test r_fitted == base_fitted + 1
        end

        # 3. compilability
        for r in result
            @test r isa EnzymeRates.Mechanism
            EnzymeRates._assert_mechanism_invariants(r)
            @test EnzymeRates.EnzymeMechanism(r) isa EnzymeMechanism
        end

        # 4. property: every variant has ≥1 I-binding step; all I-binding
        # steps in a single variant share the same kinetic group.
        for r in result
            i_groups = Int[]
            for (gi, group) in enumerate(r.steps), s in group
                bm = EnzymeRates.bound_metabolite(s)
                bm !== nothing && EnzymeRates.name(bm) === :I &&
                    push!(i_groups, gi)
            end
            @test !isempty(i_groups)
            @test length(unique(i_groups)) == 1
        end
    end

    @testset "Mechanism — Two regulators chain: J added after I preserves I" begin
        # After step A (add I) and step B (add J), J-binding steps must exist
        # and I-binding steps must remain.
        em_seed = @enzyme_mechanism begin
            substrates: S
            products: P
            steps: begin
                E + P ⇌ E(P)
                E + S ⇌ E(S)
                E(S) <--> E(P)
            end
        end
        m = EnzymeRates.Mechanism(em_seed)
        EnzymeRates._assert_mechanism_invariants(m)
        rxn = @enzyme_reaction begin
            substrates: S[C]
            products: P[C]
            dead_end_inhibitors: I, J
        end

        # Step A: 2 eligible regs (I, J), 1 form each → 2 variants total.
        i_or_j_ms = EnzymeRates._expand_add_dead_end_regulator(m, rxn)
        @test length(i_or_j_ms) == 2
        with_i = first(filter(i_or_j_ms) do r
            any(r.steps) do group
                any(group) do s
                    bm = EnzymeRates.bound_metabolite(s)
                    bm !== nothing && EnzymeRates.name(bm) === :I
                end
            end
        end)

        # Step B: add J on top of the I-bound variant.
        j_ms = EnzymeRates._expand_add_dead_end_regulator(with_i, rxn)
        @test !isempty(j_ms)

        # Property: each result has ≥1 J-binding step AND ≥1 I-binding
        # step (adding J must not remove I-binding steps).
        for r in j_ms
            j_present = false
            i_present = false
            for group in r.steps, s in group
                bm = EnzymeRates.bound_metabolite(s)
                bm === nothing && continue
                EnzymeRates.name(bm) === :J && (j_present = true)
                EnzymeRates.name(bm) === :I && (i_present = true)
            end
            @test j_present
            @test i_present
        end
    end

    @testset "Mechanism — Two regulators competition: 17 variants" begin
        # Bi-bi random with two dead-end inhibitors. Step A: add both regs
        # (9 + 9 = 18 variants). Step B: pick variant with I1 at ≥2 forms
        # and add I2 on top → 17 variants (one less than 18 after dedup
        # accounts for I2's active-form set intersection with existing I1).
        em_seed = @enzyme_mechanism begin
            substrates: A, B
            products: P, Q
            steps: begin
                (E + A ⇌ E(A), E(B) + A ⇌ E(A, B))
                (E + B ⇌ E(B), E(A) + B ⇌ E(A, B))
                (E + P ⇌ E(P), E(Q) + P ⇌ E(P, Q))
                (E + Q ⇌ E(Q), E(P) + Q ⇌ E(P, Q))
                E(A, B) <--> E(P, Q)
            end
        end
        m = EnzymeRates.Mechanism(em_seed)
        EnzymeRates._assert_mechanism_invariants(m)
        rxn = @enzyme_reaction begin
            substrates: A[C], B[N]
            products: P[C], Q[N]
            dead_end_inhibitors: I1, I2
        end

        # Step A
        result1 = EnzymeRates._expand_add_dead_end_regulator(m, rxn)
        @test length(result1) == 18

        # Pick variant where I1 binds at multiple forms.
        multi = filter(result1) do r
            i1_forms = Set{Symbol}()
            for group in r.steps, s in group
                bm = EnzymeRates.bound_metabolite(s)
                if bm !== nothing && EnzymeRates.name(bm) === :I1
                    push!(i1_forms, EnzymeRates.name(
                        EnzymeRates.to_species(s)))
                end
            end
            length(i1_forms) >= 2
        end
        @test !isempty(multi)
        m_i1 = first(multi)

        # Step B
        result2 = EnzymeRates._expand_add_dead_end_regulator(m_i1, rxn)
        @test length(result2) == 17

        # Property: ≥1 variant has I1 + I2 coexisting on the same enzyme
        # form (non-competing); ≥1 variant has I2 forms that never
        # coexist with I1 (fully-competing). Dead-end species carry the
        # inhibitor as a `CompetitiveInhibitor` in `bound`, rendered with
        # an `inh` marker in the form name (e.g., `:E_I1inh_I2inh`). The
        # substring checks below match the bare inhibitor names within
        # those rendered form names. Collect ALL form names per variant
        # from both from_species and to_species across every step.
        function _all_forms(r)
            forms = String[]
            for group in r.steps, s in group
                push!(forms, string(EnzymeRates.name(
                    EnzymeRates.from_species(s))))
                push!(forms, string(EnzymeRates.name(
                    EnzymeRates.to_species(s))))
            end
            forms
        end
        has_coexist = any(result2) do r
            any(f -> contains(f, "I1") && contains(f, "I2"),
                _all_forms(r))
        end
        @test has_coexist
        has_compete = any(result2) do r
            forms = _all_forms(r)
            has_i2 = any(f -> contains(f, "I2"), forms)
            no_coexist = !any(f ->
                contains(f, "I1") && contains(f, "I2"), forms)
            has_i2 && no_coexist
        end
        @test has_compete
    end

    @testset "Mechanism — Substrate-as-dead-end-inhibitor overlap" begin
        # :S declared as BOTH substrate and dead-end inhibitor. The move
        # must treat the substrate-:S and the inhibitor-:S binding kinetic
        # groups as independent — Mechanism stores the inhibitor as a
        # CompetitiveInhibitor Metabolite (separate type from the
        # Substrate), so there is no name collision in bound_metabolite.
        rxn_overlap = @enzyme_reaction begin
            substrates: S[C]
            products: P[C]
            dead_end_inhibitors: S
        end
        init_ms = EnzymeRates.init_mechanisms(rxn_overlap)
        @test length(init_ms) == 1
        m = first(init_ms)
        EnzymeRates._assert_mechanism_invariants(m)

        result = EnzymeRates._expand_add_dead_end_regulator(m, rxn_overlap)

        # 1. count: 1 variant (uni-uni; eligible form = E only).
        @test length(result) == 1

        # 2. Δ params: +1 (one new K_I parameter), measured against
        # ground-truth `fitted_params(compile_mechanism(...))`. The
        # dead-end inhibitor binds as a `CompetitiveInhibitor`, so its
        # form renders `:ESinh` — distinct from the substrate-bound
        # `:ES` (`Species([Substrate(:S)], :E)`). Compiled `fitted_params`
        # is the canonical source for exact counts.
        base_fitted = length(EnzymeRates.fitted_params(
            EnzymeRates.compile_mechanism(m)))
        for r in result
            r_fitted = length(EnzymeRates.fitted_params(
                EnzymeRates.compile_mechanism(r)))
            @test r_fitted == base_fitted + 1
        end

        # 3. compilability
        for r in result
            @test r isa EnzymeRates.Mechanism
            EnzymeRates._assert_mechanism_invariants(r)
            @test EnzymeRates.EnzymeMechanism(r) isa EnzymeMechanism
        end

        # 4. property: a new step has bound_metabolite::CompetitiveInhibitor
        # named :S (proving the substrate-:S vs inhibitor-:S distinction is
        # preserved). Exactly one new outer-vector kinetic group was added
        # for the inhibitor binding.
        for r in result
            has_inh_s = any(r.steps) do group
                any(group) do s
                    bm = EnzymeRates.bound_metabolite(s)
                    bm isa EnzymeRates.CompetitiveInhibitor &&
                        EnzymeRates.name(bm) === :S
                end
            end
            @test has_inh_s
            @test length(r.steps) == length(m.steps) + 1
        end
    end

    @testset "AllostericMechanism — dead-end binding tagged :EqualAI" begin
        # Uni-uni allosteric (catalytic_n=2) with mixed regs:
        # :I dead-end, :R allosteric. The move must (a) exclude :R from
        # eligible regs (allosteric ligand), (b) append the new dead-end
        # group's cat_allo_states tag as :EqualAI.
        em_seed = @allosteric_mechanism begin
            substrates: S
            products: P
            catalytic_multiplicity: 2
            catalytic_steps: begin
                E + P ⇌ E(P)    :: EqualAI
                E + S ⇌ E(S)    :: EqualAI
                E(S) <--> E(P)  :: EqualAI
            end
        end
        am = EnzymeRates.AllostericMechanism(em_seed)
        EnzymeRates._assert_mechanism_invariants(am)
        rxn = @enzyme_reaction begin
            substrates: S[C]
            products: P[C]
            dead_end_inhibitors: I
            allosteric_regulators: R
            oligomeric_state: 2
        end

        result = EnzymeRates._expand_add_dead_end_regulator(am, rxn)

        # 1. count: 1 (same as plain uni-uni + I; :R is excluded as
        # allosteric).
        @test length(result) == 1

        # 2. Δ params: +1 (:EqualAI new group → one shared K), measured
        # against ground-truth `fitted_params(compile_mechanism(...))` —
        # exact parameter counts come from the compiled mechanism.
        base_fitted = length(EnzymeRates.fitted_params(
            EnzymeRates.compile_mechanism(am)))
        for r in result
            r_fitted = length(EnzymeRates.fitted_params(
                EnzymeRates.compile_mechanism(r)))
            @test r_fitted == base_fitted + 1
        end

        # 3. compilability
        for r in result
            @test r isa EnzymeRates.AllostericMechanism
            EnzymeRates._assert_mechanism_invariants(r)
            @test EnzymeRates.AllostericEnzymeMechanism(r) isa
                AllostericEnzymeMechanism
        end

        # 4. property: exactly one new cat_steps group; its tag is
        # :EqualAI. Pre-existing tags are unchanged.
        for r in result
            @test length(r.cat_steps) == length(am.cat_steps) + 1
            @test length(r.cat_allo_states) == length(am.cat_allo_states) + 1
            # The new group's tag is :EqualAI (appended at the end by the
            # AllostericMechanism wrap kernel).
            @test r.cat_allo_states[end] == :EqualAI
            # Pre-existing tags preserved positionally.
            @test r.cat_allo_states[1:length(am.cat_allo_states)] ==
                am.cat_allo_states
        end

        # 5. preservation: catalytic_multiplicity and regulatory_sites
        # unchanged (dead-end add does not touch the regulatory-site
        # vector).
        for r in result
            @test r.catalytic_multiplicity == am.catalytic_multiplicity
            @test r.regulatory_sites == am.regulatory_sites
        end
    end

    @testset "AllostericMechanism — allosteric-only regulator → empty" begin
        # Rxn declares only :R as an allosteric regulator
        # (no dead-end inhibitors); all declared regulators are
        # allosteric ligands → eligible_regs is empty → result is empty.
        em_seed = @allosteric_mechanism begin
            substrates: S
            products: P
            catalytic_multiplicity: 2
            catalytic_steps: begin
                E + P ⇌ E(P)    :: EqualAI
                E + S ⇌ E(S)    :: EqualAI
                E(S) <--> E(P)  :: EqualAI
            end
        end
        am = EnzymeRates.AllostericMechanism(em_seed)
        @test isempty(EnzymeRates._expand_add_dead_end_regulator(
            am, uni_uni_allo_reg))
    end

    @testset "Mechanism — I-binding steps share one kinetic group" begin
        # Variants where :I binds at multiple forms MUST keep all
        # I-binding steps in a single outer-vector kinetic group (one K_I
        # parameter, not one per form). This is invariant across all
        # variants that have multi-form I binding.
        em_seed = @enzyme_mechanism begin
            substrates: A, B
            products: P, Q
            steps: begin
                (E + A ⇌ E(A), E(B) + A ⇌ E(A, B))
                (E + B ⇌ E(B), E(A) + B ⇌ E(A, B))
                (E + P ⇌ E(P), E(Q) + P ⇌ E(P, Q))
                (E + Q ⇌ E(Q), E(P) + Q ⇌ E(P, Q))
                E(A, B) <--> E(P, Q)
            end
        end
        m = EnzymeRates.Mechanism(em_seed)
        rxn = @enzyme_reaction begin
            substrates: A[C], B[N]
            products: P[C], Q[N]
            dead_end_inhibitors: I
        end
        result = EnzymeRates._expand_add_dead_end_regulator(m, rxn)

        # Pick variants with ≥2 I-binding steps (multi-form inhibitor binding).
        multi = filter(result) do r
            n = 0
            for group in r.steps, s in group
                bm = EnzymeRates.bound_metabolite(s)
                bm !== nothing && EnzymeRates.name(bm) === :I && (n += 1)
            end
            n >= 2
        end
        @test !isempty(multi)
        for r in multi
            i_groups = Int[]
            for (gi, group) in enumerate(r.steps), s in group
                bm = EnzymeRates.bound_metabolite(s)
                bm !== nothing && EnzymeRates.name(bm) === :I &&
                    push!(i_groups, gi)
            end
            @test length(unique(i_groups)) == 1
        end
    end

    @testset "Mechanism — uni-uni + I: 1 variant" begin
        # SEED: uni-uni init from rxn that declares :I as a dead-end inhibitor.
        # The Mechanism overload requires the caller to pass the declared
        # `rxn` separately because the Mechanism's own .reaction only carries
        # regulators that are already bound by its steps (`I` isn't bound yet).
        rxn = @enzyme_reaction begin
            substrates: S[C]
            products: P[C]
            dead_end_inhibitors: I
        end
        m = first(EnzymeRates.init_mechanisms(rxn))
        @test m isa EnzymeRates.Mechanism

        result = EnzymeRates._expand_add_dead_end_regulator(m, rxn)

        # 1. count: 1 variant.
        @test length(result) == 1
        for r in result
            @test r isa EnzymeRates.Mechanism
        end

        # 2. compilability: the result compiles via EnzymeMechanism.
        for r in result
            @test EnzymeRates.EnzymeMechanism(r) isa
                EnzymeRates.EnzymeMechanism
        end

        # 3. structural: a new step exists with :I bound to E.
        r1 = first(result)
        has_i_step = any(r1.steps) do group
            any(group) do s
                EnzymeRates.bound_metabolite(s) !== nothing &&
                    EnzymeRates.name(
                        EnzymeRates.bound_metabolite(s)) === :I
            end
        end
        @test has_i_step
    end

    @testset "Mechanism — exclude_regs suppresses regulator addition" begin
        rxn_ij = @enzyme_reaction begin
            substrates: S[C]
            products: P[C]
            dead_end_inhibitors: I, J
        end
        m = first(EnzymeRates.init_mechanisms(rxn_ij))

        baseline = EnzymeRates._expand_add_dead_end_regulator(m, rxn_ij)
        @test length(baseline) == 2

        excluded = EnzymeRates._expand_add_dead_end_regulator(
            m, rxn_ij; exclude_regs=Set([:I]))
        @test length(excluded) == 1
        # The remaining variant must bind :J, not :I.
        has_i = any(excluded) do r
            any(r.steps) do group
                any(group) do s
                    EnzymeRates.bound_metabolite(s) !== nothing &&
                        EnzymeRates.name(
                            EnzymeRates.bound_metabolite(s)) === :I
                end
            end
        end
        @test !has_i

        @test isempty(EnzymeRates._expand_add_dead_end_regulator(
            m, rxn_ij; exclude_regs=Set([:I, :J])))
    end

    @testset "Mechanism — no regulators: empty (negative)" begin
        m = first(EnzymeRates.init_mechanisms(uni_uni_rxn))
        rxn = @enzyme_reaction begin
            substrates: S[C]
            products: P[C]
        end
        @test isempty(EnzymeRates._expand_add_dead_end_regulator(m, rxn))
    end

    @testset "AllostericMechanism — dead-end on allosteric base" begin
        # Build an AllostericMechanism with a dead-end-eligible regulator
        # declared in the reaction (but not yet bound).
        rxn = @enzyme_reaction begin
            substrates: S[C]
            products: P[C]
            dead_end_inhibitors: I
            oligomeric_state: 2
        end
        init_mechs = EnzymeRates.init_mechanisms(rxn)
        allo_mechs = EnzymeRates._expand_to_allosteric(first(init_mechs), rxn)
        am = first(allo_mechs)

        result = EnzymeRates._expand_add_dead_end_regulator(am, rxn)
        @test !isempty(result)
        for r in result
            @test r isa EnzymeRates.AllostericMechanism
            @test r.catalytic_multiplicity == am.catalytic_multiplicity
        end
    end
end

# ═══════════════════════════════════════════════════════════════════════
# 4. Allosteric expansion moves
# ═══════════════════════════════════════════════════════════════════════

# ─── _expand_to_allosteric ─────────────────────────────────────────────
@testset "_expand_to_allosteric" begin

    @testset "Mechanism — oligomeric_state from reaction" begin
        # The catalytic_multiplicity of the resulting AllostericMechanism
        # is taken from rxn4's oligomeric_state, not hardcoded to 2.
        rxn4 = @enzyme_reaction begin
            substrates: S[C]
            products: P[C]
            oligomeric_state: 4
        end
        em_seed = @enzyme_mechanism begin
            substrates: S
            products: P
            steps: begin
                E + P ⇌ E(P)
                E + S ⇌ E(S)
                E(S) <--> E(P)
            end
        end
        m = EnzymeRates.Mechanism(em_seed)
        EnzymeRates._assert_mechanism_invariants(m)

        result = EnzymeRates._expand_to_allosteric(m, rxn4)
        @test !isempty(result)
        for r in result
            @test r isa EnzymeRates.AllostericMechanism
            EnzymeRates._assert_mechanism_invariants(r)
            @test r.catalytic_multiplicity == 4
        end
    end

    @testset "Mechanism — Bi-bi sequential: 5 groups → 6 variants" begin
        # _expand_to_allosteric emits 1 baseline + 1 :OnlyA per group.
        # 5 kinetic groups → 6 variants.
        em_seed = @enzyme_mechanism begin
            substrates: A, B
            products: P, Q
            steps: begin
                E + A ⇌ E(A)
                E(A) + B ⇌ E(A, B)
                E + Q ⇌ E(Q)
                E(Q) + P ⇌ E(P, Q)
                E(A, B) <--> E(P, Q)
            end
        end
        m = EnzymeRates.Mechanism(em_seed)
        EnzymeRates._assert_mechanism_invariants(m)
        bi_bi_allo_rxn = @enzyme_reaction begin
            substrates: A[C], B[N]
            products: P[C], Q[N]
            oligomeric_state: 2
        end

        result = EnzymeRates._expand_to_allosteric(m, bi_bi_allo_rxn)

        # 1. count: 5 kinetic groups → 1 + 5 = 6 variants.
        @test length(result) == 6

        # 2. Δ params: +1 per variant for the allosteric L, and nothing more.
        # Making a group :OnlyA can leave its bound form reverse-catalysis-
        # populated in the inactive state (connected-component pruning derives
        # Q_I from free E over ALL surviving steps, so with :EqualAI catalysis a
        # form whose direct binding is now forbidden is still reached through
        # catalysis), but the weight on that form is the reverse catalytic rate —
        # Haldane-DEPENDENT, derived from the forward rate and Keq, not fit. So it
        # is NOT in fitted_params. All 6 variants (baseline + 5 per-group :OnlyA)
        # are therefore Δ=1 (just L). Equilibrium flux is 0 for all six (verified
        # at the Keq mass-action ratio, |v_eq|/|v| ~ 1e-16).
        base_fitted = length(EnzymeRates.fitted_params(
            EnzymeRates.compile_mechanism(m)))
        deltas = sort([length(EnzymeRates.fitted_params(
            EnzymeRates.compile_mechanism(r))) - base_fitted for r in result])
        @test deltas == [1, 1, 1, 1, 1, 1]

        # 3. compilability — must produce AllostericEnzymeMechanism.
        for r in result
            @test r isa EnzymeRates.AllostericMechanism
            EnzymeRates._assert_mechanism_invariants(r)
            @test EnzymeRates.compile_mechanism(r) isa AllostericEnzymeMechanism
        end
    end

    @testset "Mechanism — Bi-bi ping-pong: n_groups + 1 variants" begin
        # SEED: bi-bi ping-pong topology. Move emits 1 baseline + 1
        # :OnlyA variant per kinetic group.
        em_seed = @enzyme_mechanism begin
            substrates: A, B
            products: P, Q
            steps: begin
                E + A ⇌ E(A)
                Estar + B ⇌ Estar(B)
                E + Q ⇌ E(Q)
                Estar + P ⇌ Estar(A, P)
                E(A) <--> Estar(A, P)
                Estar(B) ⇌ E(Q)
            end
        end
        m = EnzymeRates.Mechanism(em_seed)
        EnzymeRates._assert_mechanism_invariants(m)
        bi_bi_pp_allo_rxn = @enzyme_reaction begin
            substrates: A[CX], B[N]
            products: P[C], Q[NX]
            oligomeric_state: 2
        end

        result = EnzymeRates._expand_to_allosteric(m, bi_bi_pp_allo_rxn)

        # 1. count: n_groups + 1 variants. Verify the actual group count
        # from the seed before accepting this number.
        n_groups = length(m.steps)
        @test length(result) == n_groups + 1

        # 2. Δ params: +1 per variant for the allosteric L, and nothing more.
        # In this ping-pong topology (:EqualAI catalysis) making a binding group
        # :OnlyA forbids its direct binding, but reverse catalysis still populates
        # the bound form in the inactive state, so Q_I keeps it. The weight on
        # that form is the reverse catalytic rate — Haldane-DEPENDENT (derived
        # from the forward rate and Keq), not fit — so it is NOT in fitted_params.
        # All 7 variants (baseline + 6 per-group :OnlyA) are therefore Δ=1 (just
        # L). Equilibrium flux is 0 for all seven (verified at the Keq mass-action
        # ratio, |v_eq|/|v| ~ 1e-16).
        base_fitted = length(EnzymeRates.fitted_params(
            EnzymeRates.compile_mechanism(m)))
        deltas = sort([length(EnzymeRates.fitted_params(
            EnzymeRates.compile_mechanism(r))) - base_fitted for r in result])
        @test deltas == [1, 1, 1, 1, 1, 1, 1]

        # 3. compilability
        for r in result
            @test r isa EnzymeRates.AllostericMechanism
            EnzymeRates._assert_mechanism_invariants(r)
            @test EnzymeRates.compile_mechanism(r) isa AllostericEnzymeMechanism
        end
    end

    @testset "Mechanism — uni-uni: n_groups + 1 variants" begin
        # SEED: uni-uni init Mechanism. _expand_to_allosteric emits
        # baseline + one per group with that group :OnlyA.
        rxn = @enzyme_reaction begin
            substrates: S[C]
            products: P[C]
            oligomeric_state: 2
        end
        m = first(EnzymeRates.init_mechanisms(rxn))

        result = EnzymeRates._expand_to_allosteric(m, rxn)

        # 1. count: baseline + per-group :OnlyA = n_groups + 1.
        n_groups = length(m.steps)
        @test length(result) == n_groups + 1

        # 2. each result is an AllostericMechanism with correct
        # multiplicity and empty regulatory_sites.
        for r in result
            @test r isa EnzymeRates.AllostericMechanism
            @test r.catalytic_multiplicity == 2
            @test isempty(r.regulatory_sites)
            @test EnzymeRates.reaction(r) == EnzymeRates.reaction(m)
            @test length(r.cat_allo_states) == n_groups
        end

        # 3. baseline (first) is all-:EqualAI.
        @test all(t -> t == :EqualAI, result[1].cat_allo_states)

        # 4. each subsequent variant has exactly one :OnlyA tag.
        for i in 2:length(result)
            @test count(t -> t == :OnlyA,
                        result[i].cat_allo_states) == 1
            @test count(t -> t == :EqualAI,
                        result[i].cat_allo_states) == n_groups - 1
        end

        # 5. compilability: each AllostericMechanism compiles to an
        # AllostericEnzymeMechanism.
        for r in result
            @test EnzymeRates.compile_mechanism(r) isa
                EnzymeRates.AllostericEnzymeMechanism
        end
    end

    @testset "Mechanism — enumerates all allowed multiplicities" begin
        # Multi-valued allowed_catalytic_multiplicities → the variant set
        # (baseline + per-group :OnlyA) is emitted once per multiplicity.
        S = EnzymeRates.ReactantAtoms(EnzymeRates.Substrate(:S), [:C => 1])
        P = EnzymeRates.ReactantAtoms(EnzymeRates.Product(:P), [:C => 1])
        rxn = EnzymeRates.EnzymeReaction(
            [S, P], EnzymeRates.RegulatorMults[], Int[2, 4])
        m = first(EnzymeRates.init_mechanisms(rxn))
        allo = EnzymeRates._expand_to_allosteric(m, rxn)
        mults = Set(EnzymeRates.catalytic_multiplicity(am) for am in allo)
        @test mults == Set([2, 4])

        # Each multiplicity gets the full per-group variant set.
        n_groups = length(EnzymeRates.steps(m))
        @test length(allo) == 2 * (n_groups + 1)

        # Each multiplicity carries the same allo-state set: one all-:EqualAI
        # baseline plus one variant per group with exactly one :OnlyA.
        for cn in (2, 4)
            cn_variants = filter(
                am -> EnzymeRates.catalytic_multiplicity(am) == cn, allo)
            @test count(
                am -> all(==(:EqualAI), EnzymeRates.cat_allo_states(am)),
                cn_variants) == 1
            @test count(
                am -> count(==(:OnlyA), EnzymeRates.cat_allo_states(am)) == 1,
                cn_variants) == n_groups
        end

        # Single-valued case unchanged (regression guard).
        rxn1 = EnzymeRates.EnzymeReaction(
            [S, P], EnzymeRates.RegulatorMults[], Int[2])
        allo1 = EnzymeRates._expand_to_allosteric(
            first(EnzymeRates.init_mechanisms(rxn1)), rxn1)
        @test all(
            EnzymeRates.catalytic_multiplicity(am) == 2 for am in allo1)
    end

    @testset "AllostericMechanism — no-op (negative)" begin
        rxn = @enzyme_reaction begin
            substrates: S[C]
            products: P[C]
            oligomeric_state: 2
        end
        m = first(EnzymeRates.init_mechanisms(rxn))
        allo_variants = EnzymeRates._expand_to_allosteric(m, rxn)
        am = first(allo_variants)
        @test isempty(EnzymeRates._expand_to_allosteric(am, rxn))
    end
end

# ─── _expand_add_allosteric_regulator ──────────────────────────────────
@testset "_expand_add_allosteric_regulator" begin

    @testset "AllostericMechanism — uni-uni + first allo regulator R: 3 variants" begin
        # SEED: uni-uni allosteric with all groups :EqualAI and no
        # allosteric regulator added yet; rxn declares :R as the only
        # allo regulator.
        em_seed = @allosteric_mechanism begin
            substrates: S; products: P
            catalytic_multiplicity: 2
            catalytic_steps: begin
                E + P ⇌ E(P)    :: EqualAI
                E + S ⇌ E(S)    :: EqualAI
                E(S) <--> E(P)  :: EqualAI
            end
        end
        am = EnzymeRates.AllostericMechanism(em_seed)
        EnzymeRates._assert_mechanism_invariants(am)

        result = EnzymeRates._expand_add_allosteric_regulator(
            am, uni_uni_allo_reg)

        # 1. count: R is the only un-added allosteric regulator. 0 existing
        # reg sites → 3 non-:EqualAI tags × 1 site option (new site only) = 3.
        # :EqualAI branch is gated to "existing site with ≥1 non-:EqualAI
        # ligand" → not applicable here.
        @test length(result) == 3

        # 2. Δ params: :OnlyA/:OnlyI add one K_R each (+1); :NonequalAI adds
        # K_R and K_R_T (+2). Sorted [1, 1, 2]. Measured against the actual
        # compiled fitted count.
        base_fitted = length(EnzymeRates.fitted_params(
            EnzymeRates.compile_mechanism(am)))
        deltas = sort([length(EnzymeRates.fitted_params(
            EnzymeRates.compile_mechanism(r))) - base_fitted for r in result])
        @test deltas == [1, 1, 2]

        # 3. equivalence-style structural: each variant has exactly one
        # regulatory site holding the single ligand :R; tags across the 3
        # variants are exactly {:OnlyA, :OnlyI, :NonequalAI}. This captures
        # the same contract without depending on type-parameter encoding;
        # compiled type identity is stricter than semantic equivalence.
        triples = Set{Tuple{Vector{Symbol}, Vector{Symbol}, Int}}()
        for r in result
            @test length(r.regulatory_sites) == 1
            site = r.regulatory_sites[1]
            ligs = Symbol[EnzymeRates.name(l) for l in EnzymeRates.ligands(site)]
            states = collect(EnzymeRates.allo_states(site))
            push!(triples,
                  (ligs, states, EnzymeRates.multiplicity(site)))
        end
        @test triples == Set([
            ([:R], [:OnlyA],      am.catalytic_multiplicity),
            ([:R], [:OnlyI],      am.catalytic_multiplicity),
            ([:R], [:NonequalAI], am.catalytic_multiplicity),
        ])

        # 4. compilability + invariants on each variant.
        for r in result
            @test r isa EnzymeRates.AllostericMechanism
            EnzymeRates._assert_mechanism_invariants(r)
            @test EnzymeRates.compile_mechanism(r) isa AllostericEnzymeMechanism
        end

        # 5. preservation: catalytic side and cat_allo_states untouched;
        # new-site multiplicity inherits am.catalytic_multiplicity.
        for r in result
            @test r.catalytic_multiplicity == am.catalytic_multiplicity
            @test r.cat_allo_states == am.cat_allo_states
            @test EnzymeRates.reaction(r) == EnzymeRates.reaction(am)
            # The added site is appended last; its multiplicity matches
            # am.catalytic_multiplicity (new-site multiplicity contract).
            new_site = r.regulatory_sites[end]
            @test EnzymeRates.multiplicity(new_site) == am.catalytic_multiplicity
        end
    end

    @testset "AllostericMechanism — existing_de exclusion prevents adding bound dead-end" begin
        # Build a uni-uni AllostericMechanism that already has :I bound as a
        # dead-end (added via init→dead-end→allosteric on the Mechanism path).
        # Then `_expand_add_allosteric_regulator(am, rxn)` must exclude :I
        # (it's in existing_de). With :I the only declared regulator in
        # rxn, result is empty.
        rxn = @enzyme_reaction begin
            substrates: S[C]
            products: P[C]
            competitive_inhibitors: I
            oligomeric_state: 2
        end
        init_mechs = EnzymeRates.init_mechanisms(rxn)
        de_mechs = EnzymeRates._expand_add_dead_end_regulator(first(init_mechs), rxn)
        @test !isempty(de_mechs)
        plain_with_i = first(de_mechs)
        allo_mechs = EnzymeRates._expand_to_allosteric(plain_with_i, rxn)
        @test !isempty(allo_mechs)
        am = first(allo_mechs)
        EnzymeRates._assert_mechanism_invariants(am)
        @test isempty(EnzymeRates._expand_add_allosteric_regulator(am, rxn))
    end

    @testset "AllostericMechanism — Two regulators with site options: count = 7" begin
        # SEED: allosteric uni-uni with R1 already added as :OnlyA
        # (existing site has one
        # non-:EqualAI ligand). R2 is un-added; the :EqualAI-at-existing
        # branch fires because R1 qualifies.
        em_seed = @allosteric_mechanism begin
            substrates: S; products: P
            allosteric_regulators: R1::OnlyA
            catalytic_multiplicity: 2
            catalytic_steps: begin
                E + P ⇌ E(P)    :: EqualAI
                E + S ⇌ E(S)    :: EqualAI
                E(S) <--> E(P)  :: EqualAI
            end
        end
        am = EnzymeRates.AllostericMechanism(em_seed)
        EnzymeRates._assert_mechanism_invariants(am)

        result = EnzymeRates._expand_add_allosteric_regulator(
            am, uni_uni_allo_2reg)

        # 1. count: 3 non-:EqualAI tag flavors × 2 site options
        # (new + R1's existing) = 6. Plus 1 :EqualAI-at-existing variant
        # (gated on R1 being non-:EqualAI). → 7.
        @test length(result) == 7

        # 2. Δ params: five variants add one parameter (:OnlyA/:OnlyI plus
        # the :EqualAI-at-existing one), two :NonequalAI variants add two
        # (K_R2 + K_R2_T). Sorted [1, 1, 1, 1, 1, 2, 2]. Measured against
        # the actual compiled fitted count.
        base_fitted = length(EnzymeRates.fitted_params(
            EnzymeRates.compile_mechanism(am)))
        deltas = sort([length(EnzymeRates.fitted_params(
            EnzymeRates.compile_mechanism(r))) - base_fitted for r in result])
        @test deltas == [1, 1, 1, 1, 1, 2, 2]

        # 3. structural: every result has :R2 in some regulatory site.
        for r in result
            has_r2 = any(r.regulatory_sites) do site
                any(l -> EnzymeRates.name(l) === :R2,
                    EnzymeRates.ligands(site))
            end
            @test has_r2
        end

        # 4. compilability + invariants.
        for r in result
            @test r isa EnzymeRates.AllostericMechanism
            EnzymeRates._assert_mechanism_invariants(r)
            @test EnzymeRates.compile_mechanism(r) isa AllostericEnzymeMechanism
        end

        # 5. preservation: catalytic side and cat_allo_states untouched.
        for r in result
            @test r.catalytic_multiplicity == am.catalytic_multiplicity
            @test r.cat_allo_states == am.cat_allo_states
            @test EnzymeRates.reaction(r) == EnzymeRates.reaction(am)
        end
    end

    @testset "AllostericMechanism — EqualAI ligand reachable at existing reg site" begin
        # Same SEED as the 7-variant case (R1::OnlyA). Adding R2 must
        # produce at least one variant where R2 is :EqualAI at site 1
        # (the same site as R1).
        em_seed = @allosteric_mechanism begin
            substrates: S; products: P
            allosteric_regulators: R1::OnlyA
            catalytic_multiplicity: 2
            catalytic_steps: begin
                E + P ⇌ E(P)    :: EqualAI
                E + S ⇌ E(S)    :: EqualAI
                E(S) <--> E(P)  :: EqualAI
            end
        end
        am = EnzymeRates.AllostericMechanism(em_seed)
        EnzymeRates._assert_mechanism_invariants(am)

        result = EnzymeRates._expand_add_allosteric_regulator(
            am, uni_uni_allo_2reg)

        # 1. count: 7 (same derivation as the 7-variant case).
        @test length(result) == 7

        # 4. structural: at least one result has :R2 :EqualAI at site 1
        # (the same site as R1). This is the :EqualAI-at-existing branch,
        # gated on R1 (the existing ligand) being non-:EqualAI.
        target = findfirst(result) do r
            length(r.regulatory_sites) == 1 || return false
            site = r.regulatory_sites[1]
            ligs = EnzymeRates.ligands(site)
            states = EnzymeRates.allo_states(site)
            idx = findfirst(l -> EnzymeRates.name(l) === :R2, ligs)
            idx === nothing && return false
            states[idx] == :EqualAI
        end
        @test target !== nothing

        # 5. preservation
        for r in result
            @test r.catalytic_multiplicity == am.catalytic_multiplicity
            @test r.cat_allo_states == am.cat_allo_states
            @test EnzymeRates.reaction(r) == EnzymeRates.reaction(am)
        end
    end

    @testset "AllostericMechanism — Adding :EqualAI R2 at site with :OnlyI R1" begin
        # SEED: R1::OnlyI at site 1. Adding R2 enumerates non-:EqualAI
        # tags × 2 site options
        # plus :EqualAI-at-existing (gated on R1 being non-:EqualAI). 7 total.
        em_seed = @allosteric_mechanism begin
            substrates: S
            products: P
            allosteric_regulators: R1::OnlyI
            catalytic_multiplicity: 2
            catalytic_steps: begin
                E + P ⇌ E(P)      :: EqualAI
                E + S ⇌ E(S)      :: EqualAI
                E(S) <--> E(P)    :: EqualAI
            end
        end
        am = EnzymeRates.AllostericMechanism(em_seed)
        EnzymeRates._assert_mechanism_invariants(am)

        result = EnzymeRates._expand_add_allosteric_regulator(
            am, uni_uni_allo_2reg)

        # 1. count: 3 non-:EqualAI × 2 sites + 1 :EqualAI-at-existing = 7.
        @test length(result) == 7

        # 2. Δ params: same multiset as the :OnlyA seed — five +1 variants
        # and two :NonequalAI +2 variants. Sorted [1, 1, 1, 1, 1, 2, 2].
        # Measured against the actual compiled fitted count.
        base_fitted = length(EnzymeRates.fitted_params(
            EnzymeRates.compile_mechanism(am)))
        deltas = sort([length(EnzymeRates.fitted_params(
            EnzymeRates.compile_mechanism(r))) - base_fitted for r in result])
        @test deltas == [1, 1, 1, 1, 1, 2, 2]

        # 3. compilability
        for r in result
            @test r isa EnzymeRates.AllostericMechanism
            EnzymeRates._assert_mechanism_invariants(r)
            @test EnzymeRates.compile_mechanism(r) isa AllostericEnzymeMechanism
        end

        # 4. property: at least one variant has :R2 :EqualAI at site 1
        # (the :EqualAI-at-existing branch, gated on R1 being non-:EqualAI).
        has_eq_at_site1 = any(result) do r
            length(r.regulatory_sites) == 1 || return false
            site = r.regulatory_sites[1]
            ligs = EnzymeRates.ligands(site)
            states = EnzymeRates.allo_states(site)
            idx = findfirst(l -> EnzymeRates.name(l) === :R2, ligs)
            idx === nothing && return false
            states[idx] == :EqualAI
        end
        @test has_eq_at_site1
    end

    @testset "AllostericMechanism — Substrate-as-allosteric-regulator overlap" begin
        # SEED: allosteric uni-uni; rxn declares :S as both substrate and
        # allosteric regulator.
        rxn_allo_overlap = @enzyme_reaction begin
            substrates: S[C]
            products: P[C]
            allosteric_regulators: S
            oligomeric_state: 2
        end
        em_seed = @allosteric_mechanism begin
            substrates: S; products: P
            catalytic_multiplicity: 2
            catalytic_steps: begin
                E + P ⇌ E(P)    :: EqualAI
                E + S ⇌ E(S)    :: EqualAI
                E(S) <--> E(P)  :: EqualAI
            end
        end
        am = EnzymeRates.AllostericMechanism(em_seed)
        EnzymeRates._assert_mechanism_invariants(am)

        result = EnzymeRates._expand_add_allosteric_regulator(
            am, rxn_allo_overlap)

        # 1. count: :S is the only un-added allosteric regulator. 0 existing
        # reg sites → 3 non-:EqualAI tags × 1 site = 3. → 3.
        @test length(result) == 3

        # 2. Δ params: :OnlyA/:OnlyI add one K_S each (+1); :NonequalAI adds
        # K_S and K_S_T (+2). Sorted [1, 1, 2]. Measured against the actual
        # compiled fitted count.
        base_fitted = length(EnzymeRates.fitted_params(
            EnzymeRates.compile_mechanism(am)))
        deltas = sort([length(EnzymeRates.fitted_params(
            EnzymeRates.compile_mechanism(r))) - base_fitted for r in result])
        @test deltas == [1, 1, 2]

        # 3. structural: :S appears in some regulatory site of every result.
        # :S still plays its catalytic substrate role in the base mechanism.
        for r in result
            has_s = any(r.regulatory_sites) do site
                any(l -> EnzymeRates.name(l) === :S,
                    EnzymeRates.ligands(site))
            end
            @test has_s
        end

        # 4. compilability + invariants (explicit since dual-role is unusual).
        for r in result
            @test r isa EnzymeRates.AllostericMechanism
            EnzymeRates._assert_mechanism_invariants(r)
            @test EnzymeRates.compile_mechanism(r) isa AllostericEnzymeMechanism
        end

        # 5. preservation
        for r in result
            @test r.catalytic_multiplicity == am.catalytic_multiplicity
            @test r.cat_allo_states == am.cat_allo_states
            @test EnzymeRates.reaction(r) == EnzymeRates.reaction(am)
        end
    end

    @testset "AllostericMechanism — Two regulators at different sites" begin
        # SEED: allosteric uni-uni with R1 already at site 1; add R2 with
        # rxn declaring both. Verifies the new-site vs existing-site
        # placement split.
        rxn = @enzyme_reaction begin
            substrates: S[C]
            products: P[C]
            allosteric_regulators: R1, R2
            oligomeric_state: 2
        end
        em_seed = @allosteric_mechanism begin
            substrates: S
            products: P
            allosteric_regulators: R1::OnlyA
            catalytic_multiplicity: 2
            catalytic_steps: begin
                E + P ⇌ E(P)      :: EqualAI
                E + S ⇌ E(S)      :: EqualAI
                E(S) <--> E(P)    :: EqualAI
            end
        end
        am = EnzymeRates.AllostericMechanism(em_seed)
        EnzymeRates._assert_mechanism_invariants(am)

        result = EnzymeRates._expand_add_allosteric_regulator(am, rxn)

        # 1. count: 3 non-:EqualAI × 2 sites + 1 :EqualAI-at-existing = 7.
        @test length(result) == 7

        # 4. property-style: separate new-site (#sites grows to 2) and
        # existing-site (#sites stays at 1) placements.
        new_site_variants = filter(r -> length(r.regulatory_sites) == 2, result)
        existing_site_variants =
            filter(r -> length(r.regulatory_sites) == 1, result)
        @test length(new_site_variants) == 3   # 3 non-:EqualAI × new site
        @test length(existing_site_variants) == 4  # 3 non-:EqualAI + 1 :EqualAI
    end

    @testset "AllostericMechanism — Product-as-allosteric-regulator overlap" begin
        # SEED: uni-uni allosteric where product :P is ALSO declared as an allosteric
        # regulator. Adding :P as allo regulator → 3 tag variants × 1 site = 3.
        rxn = @enzyme_reaction begin
            substrates: S[C]
            products: P[C]
            allosteric_regulators: P
            oligomeric_state: 2
        end
        em_seed = @allosteric_mechanism begin
            substrates: S
            products: P
            catalytic_multiplicity: 2
            catalytic_steps: begin
                E + P ⇌ E(P)      :: EqualAI
                E + S ⇌ E(S)      :: EqualAI
                E(S) <--> E(P)    :: EqualAI
            end
        end
        am = EnzymeRates.AllostericMechanism(em_seed)
        EnzymeRates._assert_mechanism_invariants(am)

        result = EnzymeRates._expand_add_allosteric_regulator(am, rxn)

        # 1. count: 3 non-:EqualAI tags × 1 new site = 3 variants.
        @test length(result) == 3

        # 2. Δ params: :OnlyA/:OnlyI add one K_P each (+1); :NonequalAI adds
        # K_P and K_P_T (+2). Sorted [1, 1, 2]. Measured against the actual
        # compiled fitted count.
        base_fitted = length(EnzymeRates.fitted_params(
            EnzymeRates.compile_mechanism(am)))
        deltas = sort([length(EnzymeRates.fitted_params(
            EnzymeRates.compile_mechanism(r))) - base_fitted for r in result])
        @test deltas == [1, 1, 2]

        # 3. compilability
        for r in result
            @test r isa EnzymeRates.AllostericMechanism
            EnzymeRates._assert_mechanism_invariants(r)
            @test EnzymeRates.compile_mechanism(r) isa AllostericEnzymeMechanism
        end

        # 4. property: :P appears in some regulatory site of every result.
        for r in result
            @test any(r.regulatory_sites) do site
                any(l -> EnzymeRates.name(l) === :P,
                    EnzymeRates.ligands(site))
            end
        end
    end

    @testset "AllostericMechanism — all declared regs already present → empty" begin
        # SEED: allosteric uni-uni with :R already bound; rxn declares only :R → eligible_regs
        # is empty → result is empty.
        rxn = @enzyme_reaction begin
            substrates: S[C]
            products: P[C]
            allosteric_regulators: R
            oligomeric_state: 2
        end
        em_seed = @allosteric_mechanism begin
            substrates: S
            products: P
            allosteric_regulators: R::OnlyA
            catalytic_multiplicity: 2
            catalytic_steps: begin
                E + P ⇌ E(P)      :: EqualAI
                E + S ⇌ E(S)      :: EqualAI
                E(S) <--> E(P)    :: EqualAI
            end
        end
        am = EnzymeRates.AllostericMechanism(em_seed)
        EnzymeRates._assert_mechanism_invariants(am)
        @test isempty(EnzymeRates._expand_add_allosteric_regulator(am, rxn))
    end

    @testset "AllostericMechanism — uni-uni + R: enumerate variants" begin
        # SEED: allosteric uni-uni with NO allosteric regulator bound yet,
        # but :R declared in the reaction. The Mechanism overload requires
        # passing the declared rxn because the AllostericMechanism's own
        # .reaction strips not-yet-bound regulators.
        rxn = @enzyme_reaction begin
            substrates: S[C]
            products: P[C]
            allosteric_regulators: R
            oligomeric_state: 2
        end
        # Build a baseline AllostericMechanism from the init → allo move.
        init_mechs = EnzymeRates.init_mechanisms(rxn)
        allo_mechs = EnzymeRates._expand_to_allosteric(first(init_mechs), rxn)
        am = first(allo_mechs)

        result = EnzymeRates._expand_add_allosteric_regulator(am, rxn)

        # 1. non-empty: at least one variant adds :R.
        @test !isempty(result)

        # 2. each result is an AllostericMechanism preserving multiplicity.
        for r in result
            @test r isa EnzymeRates.AllostericMechanism
            @test r.catalytic_multiplicity == am.catalytic_multiplicity
        end

        # 3. :R appears in regulatory_sites for every variant.
        for r in result
            has_r = any(r.regulatory_sites) do site
                any(l -> EnzymeRates.name(l) === :R,
                    EnzymeRates.ligands(site))
            end
            @test has_r
        end
    end

    @testset "Mechanism — no-op (negative)" begin
        rxn = @enzyme_reaction begin
            substrates: S[C]
            products: P[C]
            allosteric_regulators: R
            oligomeric_state: 2
        end
        m = first(EnzymeRates.init_mechanisms(rxn))
        @test isempty(EnzymeRates._expand_add_allosteric_regulator(m, rxn))
    end

end

# ─── _expand_change_allo_state ─────────────────────────────────────────
@testset "_expand_change_allo_state" begin

    @testset "AllostericMechanism — regulator tag removal delta" begin
        # SEED: uni-uni allosteric with one regulator R tagged :OnlyA;
        # all 3 catalytic groups :EqualAI.
        em_seed = @allosteric_mechanism begin
            substrates: S; products: P
            allosteric_regulators: R::OnlyA
            catalytic_multiplicity: 2
            catalytic_steps: begin
                E + P ⇌ E(P)    :: EqualAI
                E + S ⇌ E(S)    :: EqualAI
                E(S) <--> E(P)  :: EqualAI
            end
        end
        am = EnzymeRates.AllostericMechanism(em_seed)
        EnzymeRates._assert_mechanism_invariants(am)

        result = EnzymeRates._expand_change_allo_state(am)

        # 1. count: 3 cat-group relaxations + 1 reg-ligand relaxation = 4.
        @test length(result) == 4

        # 2. exactly one variant has the R ligand flipped to :NonequalAI
        # (cat_allo_states preserved); locate it structurally.
        r_removal = filter(result) do r
            r.cat_allo_states == am.cat_allo_states &&
                any(s -> any(t -> t == :NonequalAI,
                             EnzymeRates.allo_states(s)),
                    r.regulatory_sites)
        end
        @test length(r_removal) == 1

        # 3. ground-truth Δ params via compiled `fitted_params`: the
        # R-ligand :OnlyA → :NonequalAI flip adds exactly one new
        # independent parameter (Δ=+1).
        seed_truth = length(EnzymeRates.fitted_params(
            EnzymeRates.compile_mechanism(am)))
        @test length(EnzymeRates.fitted_params(
            EnzymeRates.compile_mechanism(only(r_removal)))) ==
            seed_truth + 1

        # 4. compilability + invariants on each variant.
        for r in result
            @test r isa EnzymeRates.AllostericMechanism
            EnzymeRates._assert_mechanism_invariants(r)
            @test EnzymeRates.compile_mechanism(r) isa AllostericEnzymeMechanism
        end
    end

    @testset "AllostericMechanism — :OnlyI regulator-ligand relaxation" begin
        # SEED: uni-uni allosteric with R::OnlyI and 3 :EqualAI cat
        # groups. Each non-:NonequalAI entry contributes one variant:
        # 3 cat-group + 1 reg-ligand = 4.
        em_seed = @allosteric_mechanism begin
            substrates: S
            products: P
            allosteric_regulators: R::OnlyI
            catalytic_multiplicity: 2
            catalytic_steps: begin
                E + P ⇌ E(P)      :: EqualAI
                E + S ⇌ E(S)      :: EqualAI
                E(S) <--> E(P)    :: EqualAI
            end
        end
        am = EnzymeRates.AllostericMechanism(em_seed)
        EnzymeRates._assert_mechanism_invariants(am)

        result = EnzymeRates._expand_change_allo_state(am)

        # 1. count: 3 cat-group relaxations + 1 reg-ligand relaxation.
        @test length(result) == 4

        # 2. ground-truth Δ multiset via `fitted_params`. Under strict
        # `:EqualAI`, relaxing a single binding group to :NonequalAI while its
        # partners stay :EqualAI is degenerate: a thermodynamic cycle forbids
        # that affinity split, so it collapses (K_I = K_A) and adds NO fitted
        # parameter. Only the catalytic relaxation (its derived reverse absorbs
        # the cycle) and the reg-ligand relaxation each add one. So the two
        # binding relaxations contribute 0: `[0, 0, 1, 1]`. (The enumerator
        # will skip such degenerate configs in a follow-up PR.)
        seed_truth = length(EnzymeRates.fitted_params(
            EnzymeRates.compile_mechanism(am)))
        truth_deltas = sort([
            length(EnzymeRates.fitted_params(
                EnzymeRates.compile_mechanism(r))) - seed_truth
            for r in result
        ])
        @test truth_deltas == [0, 0, 1, 1]

        # 3. compilability + invariants.
        for r in result
            @test r isa EnzymeRates.AllostericMechanism
            EnzymeRates._assert_mechanism_invariants(r)
            @test EnzymeRates.compile_mechanism(r) isa AllostericEnzymeMechanism
        end

        # 4. exactly one ligand-relaxation variant: the R ligand's tag
        # flipped from :OnlyI to :NonequalAI (located structurally).
        n_r_relaxed = count(result) do r
            r.cat_allo_states == am.cat_allo_states &&
                any(s -> any(t -> t == :NonequalAI,
                             EnzymeRates.allo_states(s)),
                    r.regulatory_sites)
        end
        @test n_r_relaxed == 1
    end

    @testset "AllostericMechanism — multiple regulator ligands at independent tags" begin
        # SEED: allosteric uni-uni with R1::OnlyA + R2::OnlyI at the same
        # regulatory site, all 3 cat groups :EqualAI. Each non-:NonequalAI
        # entry contributes one variant: 3 cat + 2 reg-ligand = 5.
        em_seed = @allosteric_mechanism begin
            substrates: S
            products: P
            allosteric_regulators: R1::OnlyA, R2::OnlyI
            catalytic_multiplicity: 2
            catalytic_steps: begin
                E + P ⇌ E(P)      :: EqualAI
                E + S ⇌ E(S)      :: EqualAI
                E(S) <--> E(P)    :: EqualAI
            end
        end
        am = EnzymeRates.AllostericMechanism(em_seed)
        EnzymeRates._assert_mechanism_invariants(am)

        result = EnzymeRates._expand_change_allo_state(am)

        # 1. count: 3 cat-group relaxations + 2 reg-ligand relaxations = 5.
        @test length(result) == 5

        # 2. ground-truth Δ multiset via `fitted_params`. Under strict
        # `:EqualAI`, relaxing a single binding group to :NonequalAI while its
        # partners stay :EqualAI collapses that affinity (K_I = K_A) and adds NO
        # fitted parameter (a thermodynamic cycle forbids the split). The two
        # binding relaxations contribute 0; the catalytic relaxation and the two
        # reg-ligand relaxations each add one: `[0, 0, 1, 1, 1]`. (The enumerator
        # will skip such degenerate configs in a follow-up PR.)
        seed_truth = length(EnzymeRates.fitted_params(
            EnzymeRates.compile_mechanism(am)))
        truth_deltas = sort([
            length(EnzymeRates.fitted_params(
                EnzymeRates.compile_mechanism(r))) - seed_truth
            for r in result
        ])
        @test truth_deltas == [0, 0, 1, 1, 1]

        # 3. compilability + invariants.
        for r in result
            @test r isa EnzymeRates.AllostericMechanism
            EnzymeRates._assert_mechanism_invariants(r)
            @test EnzymeRates.compile_mechanism(r) isa AllostericEnzymeMechanism
        end

        # 4. property: each ligand has exactly one variant where ONLY it
        # is relaxed to :NonequalAI (cat states preserved, the other
        # ligand keeps its seed tag). Locate by per-ligand state.
        function _site_lig_state(r, lig_name)
            for site in r.regulatory_sites
                for (l, st) in zip(EnzymeRates.ligands(site),
                                   EnzymeRates.allo_states(site))
                    EnzymeRates.name(l) === lig_name && return st
                end
            end
            error("ligand $lig_name not found")
        end
        n_r1_relaxed_only = count(result) do r
            r.cat_allo_states == am.cat_allo_states &&
                _site_lig_state(r, :R1) == :NonequalAI &&
                _site_lig_state(r, :R2) == :OnlyI
        end
        n_r2_relaxed_only = count(result) do r
            r.cat_allo_states == am.cat_allo_states &&
                _site_lig_state(r, :R1) == :OnlyA &&
                _site_lig_state(r, :R2) == :NonequalAI
        end
        @test n_r1_relaxed_only == 1
        @test n_r2_relaxed_only == 1
    end

    @testset "AllostericMechanism — non-default site multiplicity (cat=4, reg=2)" begin
        # SEED: catalytic 4-mer with R::OnlyA at a multiplicity-2 reg site
        # (less than catalytic_multiplicity).
        em_seed = @allosteric_mechanism begin
            substrates: S
            products: P
            allosteric_regulators: R::OnlyA
            catalytic_multiplicity: 4
            catalytic_steps: begin
                E + P ⇌ E(P)      :: EqualAI
                E + S ⇌ E(S)      :: EqualAI
                E(S) <--> E(P)    :: EqualAI
            end
            regulatory_site(multiplicity = 2): begin
                ligands: R
            end
        end
        am = EnzymeRates.AllostericMechanism(em_seed)
        EnzymeRates._assert_mechanism_invariants(am)
        @test am.catalytic_multiplicity == 4
        @test [EnzymeRates.multiplicity(s) for s in am.regulatory_sites] == [2]

        # _expand_change_allo_state must preserve both
        # catalytic_multiplicity and per-site multiplicity independently.
        result = EnzymeRates._expand_change_allo_state(am)
        @test !isempty(result)
        for r in result
            @test r isa EnzymeRates.AllostericMechanism
            EnzymeRates._assert_mechanism_invariants(r)
            @test r.catalytic_multiplicity == 4
            @test [EnzymeRates.multiplicity(s) for s in r.regulatory_sites] == [2]
        end
    end

    @testset "AllostericMechanism — uni-uni all-:EqualAI: 3 cat relaxations" begin
        # SEED: uni-uni allosteric with all 3 catalytic groups tagged
        # :EqualAI and no regulatory sites. Each non-:NonequalAI cat-group
        # tag contributes one variant (flip to :NonequalAI).
        rxn = @enzyme_reaction begin
            substrates: S[C]
            products: P[C]
            oligomeric_state: 2
        end
        init_mechs = EnzymeRates.init_mechanisms(rxn)
        allo_mechs = EnzymeRates._expand_to_allosteric(first(init_mechs), rxn)
        # Pick the all-:EqualAI baseline (first variant per impl).
        am = first(allo_mechs)
        @test all(t -> t == :EqualAI, am.cat_allo_states)

        result = EnzymeRates._expand_change_allo_state(am)

        # 1. count: 3 cat-group relaxations + 0 reg-ligand relaxations.
        @test length(result) == length(am.cat_allo_states)

        # 2. each variant has exactly one :NonequalAI entry in
        # cat_allo_states; the rest match the seed.
        for r in result
            @test r isa EnzymeRates.AllostericMechanism
            @test count(t -> t == :NonequalAI, r.cat_allo_states) == 1
            for i in 1:length(am.cat_allo_states)
                @test r.cat_allo_states[i] == :NonequalAI ||
                      r.cat_allo_states[i] == am.cat_allo_states[i]
            end
        end

        # 3. preservation: reaction, multiplicity, sites.
        for r in result
            @test EnzymeRates.reaction(r) == EnzymeRates.reaction(am)
            @test r.catalytic_multiplicity == am.catalytic_multiplicity
            @test r.regulatory_sites == am.regulatory_sites
        end
    end

    @testset "AllostericMechanism — already-:NonequalAI: empty (negative)" begin
        # If every tag is :NonequalAI, no relaxation is possible.
        rxn = @enzyme_reaction begin
            substrates: S[C]
            products: P[C]
            oligomeric_state: 2
        end
        init_mechs = EnzymeRates.init_mechanisms(rxn)
        allo_mechs = EnzymeRates._expand_to_allosteric(first(init_mechs), rxn)
        am_seed = first(allo_mechs)
        # Manually build an all-:NonequalAI version.
        am_all_neq = EnzymeRates.AllostericMechanism(
            EnzymeRates.reaction(am_seed),
            copy(am_seed.cat_steps),
            fill(:NonequalAI, length(am_seed.cat_allo_states)),
            am_seed.catalytic_multiplicity,
            copy(am_seed.regulatory_sites))
        @test isempty(EnzymeRates._expand_change_allo_state(am_all_neq))
    end

    @testset "Mechanism — no-op (negative)" begin
        m = first(EnzymeRates.init_mechanisms(uni_uni_rxn))
        @test isempty(EnzymeRates._expand_change_allo_state(m))
    end

end

# ═══════════════════════════════════════════════════════════════════════
# 5. Composition (dedup, expand_mechanisms)
# ═══════════════════════════════════════════════════════════════════════

# ─── canonical by construction ─────────────────────────────────────────

@testset "Mechanism — canonical by construction" begin
    # Test 1: outer kinetic-group order does not matter. Building from the
    # same steps with the outer groups in reversed order yields a
    # struct-equal Mechanism (the basis for dedup via `unique!`).
    m_seed = first(EnzymeRates.init_mechanisms(uni_uni_rxn))
    EnzymeRates._assert_mechanism_invariants(m_seed)
    m_perm = EnzymeRates.Mechanism(
        EnzymeRates.reaction(m_seed), reverse(m_seed.steps))
    @test m_seed == m_perm

    # Test 2: AllostericMechanism — site permutation with DISTINCT
    # multiplicities. The constructor must permute cat_allo_states alongside
    # cat_steps (catalytic side) and produce the same regulatory_sites
    # ordering (regulatory side) regardless of input site order.
    base = first(EnzymeRates.init_mechanisms(uni_uni_allo))
    cat_states = [:EqualAI for _ in base.steps]
    site_a = EnzymeRates.RegulatorySite(
        [EnzymeRates.AllostericRegulator(:A)], 2, [:OnlyA])
    site_b = EnzymeRates.RegulatorySite(
        [EnzymeRates.AllostericRegulator(:B)], 4, [:OnlyI])
    am_ab = EnzymeRates.AllostericMechanism(
        EnzymeRates.reaction(base),
        [copy(g) for g in base.steps], cat_states, 2, [site_a, site_b])
    am_ba = EnzymeRates.AllostericMechanism(   # sites swapped
        EnzymeRates.reaction(base),
        [copy(g) for g in base.steps], cat_states, 2, [site_b, site_a])
    EnzymeRates._assert_mechanism_invariants(am_ab)
    EnzymeRates._assert_mechanism_invariants(am_ba)
    @test am_ab == am_ba
end

# ─── _dedup_key ────────────────────────────────────────────────────────

# Mechanism dedup keys: struct equality. Mechanisms are canonical at
# construction, so dedup is just `unique!`, which relies on
# `Base.==` / `Base.hash` on the struct itself. This testset locks in the
# struct-equality contract that powers that dedup.
@testset "Mechanism — dedup key via struct equality" begin
    # Same content → equal (and equal hashes).
    m_seed = first(EnzymeRates.init_mechanisms(uni_uni_rxn))
    EnzymeRates._assert_mechanism_invariants(m_seed)
    m_copy = EnzymeRates.Mechanism(
        EnzymeRates.reaction(m_seed),
        Vector{EnzymeRates.Step}[copy(g) for g in m_seed.steps])
    @test m_seed == m_copy
    @test hash(m_seed) == hash(m_copy)

    # AllostericMechanism: differing site multiplicities → unequal because
    # multiplicities are part of the structural identity.
    base = first(EnzymeRates.init_mechanisms(uni_uni_allo))
    cat_states = [:EqualAI for _ in base.steps]
    am_m2 = EnzymeRates.AllostericMechanism(
        EnzymeRates.reaction(base),
        [copy(g) for g in base.steps], cat_states, 2,
        [EnzymeRates.RegulatorySite(
            [EnzymeRates.AllostericRegulator(:A)], 2, [:OnlyA])])
    am_m4 = EnzymeRates.AllostericMechanism(   # different multiplicity
        EnzymeRates.reaction(base),
        [copy(g) for g in base.steps], cat_states, 2,
        [EnzymeRates.RegulatorySite(
            [EnzymeRates.AllostericRegulator(:A)], 4, [:OnlyA])])
    EnzymeRates._assert_mechanism_invariants(am_m2)
    EnzymeRates._assert_mechanism_invariants(am_m4)
    @test am_m2 != am_m4
end

# ─── mechanism dedup (unique!) ──────────────────────────────────────────
@testset "Dedup" begin

    @testset "Mechanism — same physics, different group order" begin
        # Two Mechanisms representing the same physics but with their
        # outer kinetic-group order swapped should collapse to one.
        m_seed = first(EnzymeRates.init_mechanisms(uni_uni_rxn))
        @test m_seed isa EnzymeRates.Mechanism
        # Build a permuted copy by reversing the group order.
        permuted_steps = reverse(m_seed.steps)
        m_perm = EnzymeRates.Mechanism(
            EnzymeRates.reaction(m_seed), permuted_steps)
        v = EnzymeRates.Mechanism[m_seed, m_perm]
        unique!(v)
        @test length(v) == 1
    end

    @testset "Mechanism — different mechanisms preserved" begin
        # Surviving Mechanisms must be pairwise distinct under
        # compile-time equality (EnzymeMechanism singleton type).
        mechs = collect(EnzymeRates.init_mechanisms(bi_bi_rxn))
        unique!(mechs)
        compiled = Set(EnzymeRates.EnzymeMechanism(m) for m in mechs)
        @test length(mechs) == length(compiled)
        @test length(mechs) >= 2
    end

    @testset "Mechanism — idempotent" begin
        mechs = collect(EnzymeRates.init_mechanisms(bi_bi_rxn))
        unique!(mechs)
        n1 = length(mechs)
        unique!(mechs)
        @test length(mechs) == n1
    end

    @testset "Mechanism — bi-bi init: dedup leaves canonical seeds intact" begin
        # init_mechanisms produces mechanisms that are already in canonical
        # form (no two are presentation-variants of each other). unique!
        # is therefore a no-op on the count.
        mechs = collect(EnzymeRates.init_mechanisms(bi_bi_rxn))
        n = length(mechs)
        unique!(mechs)
        @test length(mechs) == n
    end

    @testset "AllostericMechanism — same physics, site permutation" begin
        # Build two AllostericMechanisms representing the same physics with
        # sites in different order.
        base = first(EnzymeRates.init_mechanisms(uni_uni_allo))
        cat_states = [:NonequalAI for _ in base.steps]
        site_a = EnzymeRates.RegulatorySite(
            [EnzymeRates.AllostericRegulator(:A)], 2, [:NonequalAI])
        site_b = EnzymeRates.RegulatorySite(
            [EnzymeRates.AllostericRegulator(:B)], 2, [:NonequalAI])
        am_ab = EnzymeRates.AllostericMechanism(
            EnzymeRates.reaction(base),
            [copy(g) for g in base.steps], cat_states, 2, [site_a, site_b])
        am_ba = EnzymeRates.AllostericMechanism(
            EnzymeRates.reaction(base),
            [copy(g) for g in base.steps], cat_states, 2, [site_b, site_a])
        v = EnzymeRates.AllostericMechanism[am_ab, am_ba]
        unique!(v)
        @test length(v) == 1
    end

    @testset "Mechanism — inter-move overlap: dedup actually fires" begin
        # Run expand_mechanisms on a bi-bi init seed (Mechanism path), then
        # unique!. Assert that the flat vector shrinks, proving that two
        # different expansion paths produced equivalent Mechanisms.
        init_mechs = collect(EnzymeRates.init_mechanisms(bi_bi_rxn))
        expanded = EnzymeRates.expand_mechanisms(init_mechs, bi_bi_rxn)
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
        m_seed = first(EnzymeRates.init_mechanisms(bi_bi_rxn))
        EnzymeRates._assert_mechanism_invariants(m_seed)
        n_groups = length(m_seed.steps)
        @assert n_groups >= 3 "bi-bi init seed must have ≥3 kinetic groups " *
            "to exercise a non-trivial permutation"
        # Cyclic-rotate-by-1 permutation: [g1, g2, …, gN] → [g2, …, gN, g1].
        perm = vcat(2:n_groups, 1)
        m_rotated = EnzymeRates.Mechanism(
            EnzymeRates.reaction(m_seed), m_seed.steps[perm])
        v = EnzymeRates.Mechanism[m_seed, m_rotated]
        unique!(v)
        @test length(v) == 1
    end
end

# ─── expand_mechanisms ─────────────────────────────────────────────────
@testset "expand_mechanisms" begin

    @testset "Mechanism — Empty input" begin
        # Mechanism-form expand_mechanisms with empty input returns empty Vector.
        empty_in = Union{EnzymeRates.Mechanism,
                         EnzymeRates.AllostericMechanism}[]
        @test isempty(EnzymeRates.expand_mechanisms(empty_in, uni_uni_rxn))
    end

    @testset "Mechanism — Returns flat vector" begin
        # SEED: uni-uni RE-only, 3 singleton kinetic groups.
        em_seed = @enzyme_mechanism begin
            substrates: S
            products: P
            steps: begin
                E + P ⇌ E(P)
                E + S ⇌ E(S)
                E(S) <--> E(P)
            end
        end
        m = EnzymeRates.Mechanism(em_seed)
        EnzymeRates._assert_mechanism_invariants(m)
        result = EnzymeRates.expand_mechanisms([m], uni_uni_rxn)
        @test result isa Vector{Union{
            EnzymeRates.Mechanism, EnzymeRates.AllostericMechanism}}
        @test !isempty(result)
    end

    @testset "Mechanism — Allosteric expansion included" begin
        # SEED: uni-uni RE-only attached to an oligomeric reaction.
        # expand_mechanisms must include AllostericMechanism variants in
        # its output.
        em_seed = @enzyme_mechanism begin
            substrates: S
            products: P
            steps: begin
                E + P ⇌ E(P)
                E + S ⇌ E(S)
                E(S) <--> E(P)
            end
        end
        m = EnzymeRates.Mechanism(em_seed)
        EnzymeRates._assert_mechanism_invariants(m)
        result = EnzymeRates.expand_mechanisms([m], uni_uni_allo)
        allo_count = count(s -> s isa EnzymeRates.AllostericMechanism, result)
        # _expand_to_allosteric on a uni-uni seed with 3 kinetic groups
        # produces n_groups+1=4 AllostericMechanism variants (one per
        # group + one for L-only). Other moves do not produce
        # AllostericMechanism output from a plain Mechanism, so at least 4
        # exist.
        @test allo_count >= 4
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
        em_seed = @enzyme_mechanism begin
            substrates: S
            products: P
            steps: begin
                E + P ⇌ E(P)
                E + S ⇌ E(S)
                E(S) <--> E(P)
            end
        end
        m = EnzymeRates.Mechanism(em_seed)
        EnzymeRates._assert_mechanism_invariants(m)
        base_fitted = length(EnzymeRates.fitted_params(
            EnzymeRates.compile_mechanism(m)))
        result = EnzymeRates.expand_mechanisms([m], uni_uni_rxn)
        for child in result
            @test length(EnzymeRates.fitted_params(
                EnzymeRates.compile_mechanism(child))) >= base_fitted
        end
    end

    @testset "AllostericMechanism — rewrap preserves structure" begin
        # SEED: all-:EqualAI AllostericMechanism uni-uni. Passing this to
        # expand_mechanisms must produce AllostericMechanism expansions
        # (RE→SS rewrapped as allosteric, etc.).
        aem = @allosteric_mechanism begin
            substrates: S; products: P
            catalytic_multiplicity: 2
            catalytic_steps: begin
                E + P ⇌ E(P)    :: EqualAI
                E + S ⇌ E(S)    :: EqualAI
                E(S) <--> E(P)  :: EqualAI
            end
        end
        am = EnzymeRates.AllostericMechanism(aem)
        result = EnzymeRates.expand_mechanisms([am], uni_uni_allo)
        allo_results = filter(s -> s isa EnzymeRates.AllostericMechanism,
                              result)
        @test !isempty(allo_results)
        # Every rewrapped allosteric result must preserve the input's
        # catalytic_multiplicity and reaction — cat_steps may differ (a
        # base move may have changed them) but the allosteric-side
        # metadata is preserved.
        for r in allo_results
            @test r.catalytic_multiplicity == am.catalytic_multiplicity
            @test EnzymeRates.reaction(r) == EnzymeRates.reaction(am)
        end
    end

    @testset "AllostericMechanism — Dead-end excludes allosteric regs" begin
        # SEED: AllostericMechanism uni-uni with R already added as an
        # allosteric regulator (:OnlyA). expand_mechanisms must never add R
        # as a dead-end inhibitor — no Step in any expansion may have
        # bound_metabolite named :R (allosteric regulators live in
        # regulatory_sites, not in cat_steps).
        aem = @allosteric_mechanism begin
            substrates: S; products: P
            allosteric_regulators: R::OnlyA
            catalytic_multiplicity: 2
            catalytic_steps: begin
                E + P ⇌ E(P)    :: EqualAI
                E + S ⇌ E(S)    :: EqualAI
                E(S) <--> E(P)  :: EqualAI
            end
        end
        am = EnzymeRates.AllostericMechanism(aem)
        result = EnzymeRates.expand_mechanisms([am], uni_uni_allo_reg)
        for r in result
            groups = r isa EnzymeRates.AllostericMechanism ?
                r.cat_steps : r.steps
            for group in groups, st in group
                bm = EnzymeRates.bound_metabolite(st)
                bm === nothing && continue
                @test EnzymeRates.name(bm) !== :R
            end
        end
    end
end

# ═══════════════════════════════════════════════════════════════════════
# 6. Integration (enumerate_all)
# ═══════════════════════════════════════════════════════════════════════

# ─── enumerate_all ─────────────────────────────────────────────────────
@testset "Integration" begin

    @testset "Mechanism — Uni-uni full enumeration" begin
        # enumerate_all_mechanism buckets by ACTUAL fitted-parameter count,
        # so each bucket key `pc` must equal every member's fitted count.
        results = enumerate_all_mechanism(uni_uni_rxn; max_params=8)
        @test !isempty(results)
        pcs = sort(collect(keys(results)))
        @test issorted(pcs)
        for (pc, mechs) in results
            for m in mechs
                @test length(EnzymeRates.fitted_params(
                    EnzymeRates.compile_mechanism(m))) == pc
            end
        end
    end

    @testset "Mechanism — Bi-bi init-tier (actual-count buckets)" begin
        # Full multi-tier bi-bi enumeration would have to compile every
        # reachable mechanism to bucket it by actual fitted count (hundreds
        # of @generated derivations) — too slow for the suite. The init tier
        # suffices to verify bi-bi enumeration produces mechanisms that
        # compile and that their actual fitted-param counts fall in the
        # expected {5,6} band. (Multi-tier actual-count enumeration is
        # exercised by the uni-uni / dead-end / allosteric callers below.)
        init = unique!(
            collect(EnzymeRates.init_mechanisms(bi_bi_rxn)))
        @test !isempty(init)
        counts = Set{Int}()
        for m in init
            em = EnzymeRates.compile_mechanism(m)
            push!(counts, length(EnzymeRates.fitted_params(em)))
        end
        # {5,6}: ordered binding identifies one more thermodynamic
        # constraint than random binding, so it fits one fewer parameter.
        @test issubset(counts, Set([5, 6]))
        @test 5 in counts
    end

    @testset "Mechanism — With allosteric regulators" begin
        # Sample-based per bucket on a uni-uni allosteric reaction. The
        # bucket key `pc` IS the actual fitted-parameter count.
        rxn = @enzyme_reaction begin
            substrates: S[C]
            products: P[C]
            allosteric_regulators: R
            oligomeric_state: 2
        end
        results = enumerate_all_mechanism(rxn; max_params=8)
        has_allo = any(
            any(s isa EnzymeRates.AllostericMechanism for s in mechs)
            for (_, mechs) in results)
        @test has_allo
        for (pc, mechs) in results
            for m in first(mechs, 5)
                @test length(EnzymeRates.fitted_params(
                    EnzymeRates.compile_mechanism(m))) == pc
            end
        end
    end

    @testset "Mechanism — With dead-end regulator" begin
        # Mechanism-form parallel. Total population with a dead-end
        # regulator strictly exceeds the plain uni-uni population at the
        # same cap, since the regulator opens additional expansion moves.
        rxn = @enzyme_reaction begin
            substrates: S[C]
            products: P[C]
            dead_end_inhibitors: I
        end
        results = enumerate_all_mechanism(rxn; max_params=8)
        @test !isempty(results)
        plain = enumerate_all_mechanism(uni_uni_rxn; max_params=8)
        total_with_reg = sum(length(v) for v in values(results))
        total_plain = sum(length(v) for v in values(plain))
        @test total_with_reg > total_plain
    end

    @testset "Mechanism — Multiple levels populated" begin
        # Mechanism-form parallel. At least 2 param-count buckets, and
        # consecutive buckets separated by at most 4 (max single-move
        # delta).
        results = enumerate_all_mechanism(uni_uni_rxn; max_params=8)
        @test length(results) >= 2
        pcs = sort(collect(keys(results)))
        @test all(pcs[i+1] - pcs[i] <= 4 for i in 1:length(pcs)-1)
    end
end

# ═══════════════════════════════════════════════════════════════════════
# Testsets covering downstream concerns (canonicalization parameter-naming;
# move-on-allosteric polymorphism). Adjacent to enumeration but tests
# rate-equation-derivation and AllostericEnzymeMechanism integration.
# ═══════════════════════════════════════════════════════════════════════

@testset "Tagged groups exclude T-state params" begin
    init_mechs = EnzymeRates.init_mechanisms(uni_uni_allo)
    m_seed = first(init_mechs)
    allo_mechs = EnzymeRates._expand_to_allosteric(m_seed, uni_uni_allo)

    @testset ":OnlyA binding group: no K_T param" begin
        only_r = first(filter(allo_mechs) do am
            any(EnzymeRates.kinetic_groups(am)) do g
                EnzymeRates.cat_allo_state(am, g) === :OnlyA || return false
                # Must NOT be an iso-only group (iso `:OnlyA` is just a relabel
                # — the test wants a binding group whose K param disappears in T).
                group_steps = am.cat_steps[g]
                any(s -> EnzymeRates.bound_metabolite(s) !== nothing,
                    group_steps)
            end
        end)
        m = EnzymeRates.compile_mechanism(only_r)
        params = parameters(m)
        t_params = filter(
            p -> endswith(string(p), "_T"), params)
        @test isempty(t_params)
    end

    @testset ":OnlyA iso group: no kf_T/kr_T param" begin
        only_r_iso = first(filter(allo_mechs) do am
            any(EnzymeRates.kinetic_groups(am)) do g
                EnzymeRates.cat_allo_state(am, g) === :OnlyA || return false
                group_steps = am.cat_steps[g]
                all(s -> !EnzymeRates.is_equilibrium(s) &&
                         EnzymeRates.bound_metabolite(s) === nothing,
                    group_steps)
            end
        end)
        m = EnzymeRates.compile_mechanism(only_r_iso)
        params = parameters(m)
        t_k_params = filter(
            p -> contains(string(p), "f_T") ||
                 contains(string(p), "r_T"), params)
        @test isempty(t_k_params)
    end

    @testset "t_state_dead with :NonequalAI: K_T in body must be in parameters(Full)" begin
        # K-type allosteric uni-uni: catalytic step is :OnlyA (so
        # `_i_state_num_zero == true`), but binding steps are :NonequalAI.
        # When `_i_state_num_zero == true`, the binding partition function
        # for :NonequalAI groups must still emit K1_T / K2_T in `den_T`
        # so they appear in the rate-equation body and in parameters(Full).
        m = @allosteric_mechanism begin
            substrates: S
            products: P
            catalytic_multiplicity: 2
            catalytic_steps: begin
                E_c + S ⇌ E_c(S)      :: NonequalAI
                E_c + P ⇌ E_c(P)      :: NonequalAI
                E_c(S) <--> E_c(P)    :: OnlyA
            end
        end
        @test EnzymeRates._i_state_num_zero(EnzymeRates.AllostericMechanism(m))
        params_full = parameters(m, Full)
        # K_I_S_E_c and K_I_P_E_c are referenced in `den_T` of the body
        # (the binding partition function for :NonequalAI groups
        # is built regardless of `t_state_dead` since `den_T`
        # always appears in the denominator).
        @test :K_I_S_E_c in params_full
        @test :K_I_P_E_c in params_full
    end
end

end # top-level testset

# ═══════════════════════════════════════════════════════════════════════
# Rate-equation dedup key
# ═══════════════════════════════════════════════════════════════════════

@testset "Rate-equation dedup key" begin
    let
        # Shared exemplars used by multiple testsets below.

        # Minimal uni-uni 3-step mechanism.
        uni_uni_3step = @enzyme_mechanism begin
            substrates: S
            products: P
            steps: begin
                E + S ⇌ E(S)
                E(S) <--> E(P)
                E + P ⇌ E(P)
            end
        end

        # Random bi-uni with substrate-side mirror sharing: A-binding
        # steps share kinetic_group=1, B-binding steps share kg=2.
        biuni_mirror = @enzyme_mechanism begin
            substrates: A, B
            products: P
            steps: begin
                (E + A ⇌ E(A), E(B) + A ⇌ E(A, B))
                (E + B ⇌ E(B), E(A) + B ⇌ E(A, B))
                E(A, B) <--> E(P)
                E + P ⇌ E(P)
            end
        end

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
                (E + Lactate ⇌ E(Lactate), E(NADH) + Lactate ⇌ E(Lactate, NADH))
                (E + NAD ⇌ E(NAD), E(Lactate) + NAD ⇌ E(Lactate, NAD), E(Pyruvate) + NAD ⇌ E(NAD, Pyruvate))
                (E + NADH ⇌ E(NADH), E(Lactate) + NADH ⇌ E(Lactate, NADH))
                (E + Pyruvate ⇌ E(Pyruvate), E(NAD) + Pyruvate ⇌ E(NAD, Pyruvate), E(NADH) + Pyruvate ⇌ E(NADH, Pyruvate))
                E(NADH, Pyruvate) <--> E(Lactate, NAD)
            end
        end

        ldh_m_b = @enzyme_mechanism begin
            substrates: NADH, Pyruvate
            products: Lactate, NAD
            steps: begin
                (E + Lactate ⇌ E(Lactate), E(NAD) + Lactate ⇌ E(Lactate, NAD), E(NADH) + Lactate ⇌ E(Lactate, NADH))
                (E + NAD ⇌ E(NAD), E(Pyruvate) + NAD ⇌ E(NAD, Pyruvate))
                (E + NADH ⇌ E(NADH), E(Lactate) + NADH ⇌ E(Lactate, NADH))
                (E + Pyruvate ⇌ E(Pyruvate), E(NAD) + Pyruvate ⇌ E(NAD, Pyruvate), E(NADH) + Pyruvate ⇌ E(NADH, Pyruvate))
                E(NADH, Pyruvate) <--> E(Lactate, NAD)
            end
        end

        @testset "Distinct mechanisms produce distinct dedup keys" begin
            # Ordered binding (A-first only) vs random binding (both A-first
            # and B-first paths). Same substrate set, same product set, but
            # the random mechanism has a different enzyme-form graph. The
            # dedup key must differ.
            m_ordered = @enzyme_mechanism begin
                substrates: A, B
                products: P
                steps: begin
                    E + A ⇌ E(A)
                    E(A) + B ⇌ E(A, B)
                    E(A, B) <--> E(P)
                    E + P ⇌ E(P)
                end
            end
            m_random = @enzyme_mechanism begin
                substrates: A, B
                products: P
                steps: begin
                    E + A ⇌ E(A)
                    E + B ⇌ E(B)
                    E(A) + B ⇌ E(A, B)
                    E(B) + A ⇌ E(A, B)
                    E(A, B) <--> E(P)
                    E + P ⇌ E(P)
                end
            end
            @test EnzymeRates._rate_eq_dedup_key(rate_equation_string(m_ordered)) !=
                  EnzymeRates._rate_eq_dedup_key(rate_equation_string(m_random))
        end

        @testset "LDH Pattern-A: graph-distinct mechanisms with equivalent v share a dedup key" begin
            @test EnzymeRates._rate_eq_dedup_key(rate_equation_string(ldh_m_a)) ==
                  EnzymeRates._rate_eq_dedup_key(rate_equation_string(ldh_m_b))

            # Negative control: ldh_m_b with the E_NAD+Lactate step's
            # kinetic_group changed from 1 to 12 breaks the Lactate-
            # binding-via-NAD path's sharing with the direct E+Lactate
            # binding step. The resulting v polynomial differs.
            ldh_m_c = @enzyme_mechanism begin
                substrates: NADH, Pyruvate
                products: Lactate, NAD
                steps: begin
                    (E + Lactate ⇌ E(Lactate), E(NADH) + Lactate ⇌ E(Lactate, NADH))
                    (E + NAD ⇌ E(NAD), E(Pyruvate) + NAD ⇌ E(NAD, Pyruvate))
                    (E + NADH ⇌ E(NADH), E(Lactate) + NADH ⇌ E(Lactate, NADH))
                    (E + Pyruvate ⇌ E(Pyruvate), E(NAD) + Pyruvate ⇌ E(NAD, Pyruvate), E(NADH) + Pyruvate ⇌ E(NADH, Pyruvate))
                    E(NAD) + Lactate ⇌ E(Lactate, NAD)
                    E(NADH, Pyruvate) <--> E(Lactate, NAD)
                end
            end
            @test EnzymeRates._rate_eq_dedup_key(rate_equation_string(ldh_m_a)) !=
                  EnzymeRates._rate_eq_dedup_key(rate_equation_string(ldh_m_c))
        end

        @testset "rate_equation_string emits section labels" begin
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
            # singleton kinetic_groups (no Pass-1 absorption) preserves the
            # Wegscheider cycle relations and renders them as multi-symbol
            # RHSes in this section.
            m_weg = @enzyme_mechanism begin
                substrates: A, B
                products: P, Q
                steps: begin
                    E + A ⇌ E(A)
                    E + B ⇌ E(B)
                    E(A) + B ⇌ E(A, B)
                    E(B) + A ⇌ E(A, B)
                    E(A, B) <--> E(P, Q)
                    E(P) + Q ⇌ E(P, Q)
                    E(Q) + P ⇌ E(P, Q)
                    E + P ⇌ E(P)
                    E + Q ⇌ E(Q)
                end
            end
            s_weg = rate_equation_string(m_weg)
            @test occursin("# Wegscheider constraints:", s_weg)
        end

    end # let
end

@testset "bi-bi init_mechanisms structural golden" begin
    # Canonical, derivation-free key for one Mechanism: per kinetic group,
    # the sorted set of (from form, to form, bound metabolite, RE/SS) tuples;
    # groups themselves sorted so the key is order-independent.
    function _mech_struct_key(m::EnzymeRates.Mechanism)
        grpkeys = String[]
        for grp in m.steps
            stepkeys = sort([
                string((EnzymeRates.name(EnzymeRates.from_species(s)),
                        EnzymeRates.name(EnzymeRates.to_species(s)),
                        EnzymeRates.bound_metabolite(s) === nothing ? :iso :
                            EnzymeRates.name(EnzymeRates.bound_metabolite(s)),
                        EnzymeRates.is_equilibrium(s)))
                for s in grp])
            push!(grpkeys, join(stepkeys, "|"))
        end
        join(sort(grpkeys), " ;; ")
    end

    init = EnzymeRates.init_mechanisms(bi_bi_rxn)
    mech_keys = sort([_mech_struct_key(m) for m in init])
    @test length(unique(mech_keys)) == length(mech_keys)   # no structural dups post-dedup

    fixture = joinpath(@__DIR__, "fixtures", "phase2_init_golden.txt")
    if !isfile(fixture)
        mkpath(dirname(fixture)); write(fixture, join(mech_keys, "\n"))
        @warn "bootstrapped phase2 golden fixture; commit it" fixture
    end
    golden = readlines(fixture)
    @test mech_keys == golden          # permanent structural regression gate
    @test length(init) == length(golden)   # init_mechanisms count invariant
end

@testset "mechanism dedup via unique!" begin
    rxn = @enzyme_reaction begin
        substrates:S[C]
        products:P[C]
    end
    ms = collect(EnzymeRates.init_mechanisms(rxn))
    dup = vcat(ms, deepcopy(ms))          # every mechanism twice
    out = unique!(dup)
    @test length(out) == length(unique!(collect(ms)))
    @test length(out) <= length(dup)
    @test unique!(Union{EnzymeRates.Mechanism,
        EnzymeRates.AllostericMechanism}[]) == []
end

@testset "expand_mechanisms output is canonical" begin
    rxn = @enzyme_reaction begin
        substrates: NADH[C21H29N7O14P2], Pyruvate[C3H4O3]
        products: Lactate[C3H6O3], NAD[C21H27N7O14P2]
        oligomeric_state: 4
    end
    MECH = Union{EnzymeRates.Mechanism, EnzymeRates.AllostericMechanism}
    base = collect(EnzymeRates.init_mechanisms(rxn))
    children = EnzymeRates.expand_mechanisms(MECH[base...], rxn)
    noncanon = [c for c in children if EnzymeRates._canonical_mechanism(c) != c]
    @test isempty(noncanon)
end

@testset "init division-freeness (bi_bi_pp)" begin
    # Every enumerated init mechanism's derived rate equation must stay finite
    # when any single metabolite concentration is zero (real data has zeros).
    # `isfinite` suffices here: a residual coupling that survived the
    # concentration-GCD would put the same 1/conc factor in BOTH numerator and
    # denominator, so a zeroed metabolite yields Inf/Inf = NaN, which isfinite
    # catches (no separate != 0 check needed).
    mets = [:A, :B, :P, :Q]
    for m in unique!(collect(EnzymeRates.init_mechanisms(bi_bi_pp_rxn)))
        cm = EnzymeRates.compile_mechanism(m)
        params = random_reduced_params(cm; rng = Random.MersenneTwister(1))
        for zeroed in mets
            cvals = Tuple(n == zeroed ? 0.0 : 1.0 for n in mets)
            concs = NamedTuple{Tuple(mets)}(cvals)
            v = rate_equation(cm, concs, params)
            @test isfinite(v)
        end
    end
end
