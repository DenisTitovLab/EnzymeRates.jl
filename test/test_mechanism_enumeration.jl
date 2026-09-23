# ABOUTME: Tests for mechanism enumeration pipeline
# ABOUTME: Unit tests per move + integration tests per reaction

# Fixture rules for this file are in CLAUDE.md, "Enumeration-engine tests":
# every mechanism inline via the macro, every move test asserts the exact
# child set.

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

# A kinetic group is catalytic iff its representative step binds no
# metabolite (an isomerization/conversion step); otherwise it is a binding
# group. Mirrors the rule `_expand_to_allosteric` uses to decide whether a
# group's `:OnlyA` flip needs a paired regulator to be distinguishable.
_is_catalytic_group(m, g) =
    EnzymeRates.bound_metabolite(EnzymeRates.rep_step(m, g)) === nothing

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

@testset "shared_catalytic_site prunes catalytic-site co-occupancy" begin
rxn = @enzyme_reaction begin
    substrates: A[C], B[C]
    products:   P[C], Q[C]
    shared_catalytic_site: (A, P)
end
_reactant_names(sp) = Set(EnzymeRates.name(b)
    for b in EnzymeRates.bound(sp) if b isa EnzymeRates.Reactant)
binds_both(sp) = :A in _reactant_names(sp) && :P in _reactant_names(sp)
@test !any(binds_both(sp)
    for m in EnzymeRates.init_mechanisms(rxn)
    for g in EnzymeRates.steps(m) for s in g
    for sp in (EnzymeRates.from_species(s), EnzymeRates.to_species(s)))
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
n_base = length(unique!(collect(EnzymeRates.init_mechanisms(base))))
n_con  = length(unique!(collect(EnzymeRates.init_mechanisms(constrained))))
@test n_con < n_base
end

@testset "_add_competitive_inhibitor preserves shared_catalytic_site" begin
rxn = @enzyme_reaction begin
    substrates: A[C], B[C]
    products:   P[C], Q[C]
    shared_catalytic_site: (A, P)
end
result = EnzymeRates._add_competitive_inhibitor(rxn, :I)
@test EnzymeRates.shared_catalytic_site(result) == [(:A, :P)]
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
# Three distinct tags on three groups: S binding (RE, :EqualAI),
# P binding (RE, :NonequalAI), iso (SS, :OnlyA). The :OnlyA sits on
# the chemical step (dropped from the Haldane graph), so the mechanism
# is Haldane-valid with no one-sided :OnlyA binding.
m_compiled = @allosteric_mechanism begin
    substrates: S
    products: P
    catalytic_multiplicity: 2
    catalytic_steps: begin
        E + S ⇌ E(S)             :: EqualAI
        E + P ⇌ E(P)             :: NonequalAI
        E(S) <--> E(P)           :: OnlyA
    end
end
# Group order is canonical; allosteric tags stay bound to their steps.
am = EnzymeRates.AllostericMechanism(m_compiled)
state_of(pred) = EnzymeRates.cat_allo_state(am,
    only(g for g in EnzymeRates.kinetic_groups(am)
         if pred(EnzymeRates.bound_metabolite(EnzymeRates.rep_step(am, g)))))
@test state_of(bm -> bm isa EnzymeRates.Substrate) == :EqualAI
@test state_of(bm -> bm isa EnzymeRates.Product) == :NonequalAI
@test state_of(bm -> bm === nothing) == :OnlyA
end

# ═══════════════════════════════════════════════════════════════════════
# 2. Initialization (compile_mechanism + init_mechanisms)
# ═══════════════════════════════════════════════════════════════════════

# ─── compile_mechanism dispatch ────────────────────────────────────────
@testset "compile_mechanism dispatch" begin
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

@testset "Mechanism — flip past a bottomless RE segment" begin
    # PARENT: uni-bi with random product release, all binding RE. The product
    # ring E–E(P)–E(P, Q)–E(Q)–E needs two cuts to raise the segment count, so
    # every pair of ring edges gains. Cutting the two edges at E ({P@E, Q@E})
    # would leave {E(P), E(Q), E(P, Q)} with no bottom form — a mechanism the
    # constructor rejects — so that pair counts as failed and is not emitted.
    # Its supersets are reached by one more flip from the other pairs' children.
    # Children: the S flip plus the five other pairs.
    parent = EnzymeRates.Mechanism(@enzyme_mechanism begin
        substrates: S
        products: P, Q
        steps: begin
            E + S ⇌ E(S)
            E + P ⇌ E(P)
            E + Q ⇌ E(Q)
            E(P) + Q ⇌ E(P, Q)
            E(Q) + P ⇌ E(P, Q)
            E(S) <--> E(P, Q)
        end
    end)
    children = EnzymeRates._expand_re_to_ss(parent)
    expected = EnzymeRates.Mechanism.([
        (@enzyme_mechanism begin      # S alone: E(S) becomes its own segment
            substrates: S
            products: P, Q
            steps: begin
                E + S <--> E(S)
                E + P ⇌ E(P)
                E + Q ⇌ E(Q)
                E(P) + Q ⇌ E(P, Q)
                E(Q) + P ⇌ E(P, Q)
                E(S) <--> E(P, Q)
            end
        end),
        (@enzyme_mechanism begin      # the two edges at E(P)
            substrates: S
            products: P, Q
            steps: begin
                E + S ⇌ E(S)
                E + P <--> E(P)
                E + Q ⇌ E(Q)
                E(P) + Q <--> E(P, Q)
                E(Q) + P ⇌ E(P, Q)
                E(S) <--> E(P, Q)
            end
        end),
        (@enzyme_mechanism begin      # the two edges at E(Q)
            substrates: S
            products: P, Q
            steps: begin
                E + S ⇌ E(S)
                E + P ⇌ E(P)
                E + Q <--> E(Q)
                E(P) + Q ⇌ E(P, Q)
                E(Q) + P <--> E(P, Q)
                E(S) <--> E(P, Q)
            end
        end),
        (@enzyme_mechanism begin      # the two edges at E(P, Q)
            substrates: S
            products: P, Q
            steps: begin
                E + S ⇌ E(S)
                E + P ⇌ E(P)
                E + Q ⇌ E(Q)
                E(P) + Q <--> E(P, Q)
                E(Q) + P <--> E(P, Q)
                E(S) <--> E(P, Q)
            end
        end),
        (@enzyme_mechanism begin      # P@E and P@E(Q)
            substrates: S
            products: P, Q
            steps: begin
                E + S ⇌ E(S)
                E + P <--> E(P)
                E + Q ⇌ E(Q)
                E(P) + Q ⇌ E(P, Q)
                E(Q) + P <--> E(P, Q)
                E(S) <--> E(P, Q)
            end
        end),
        (@enzyme_mechanism begin      # Q@E and Q@E(P)
            substrates: S
            products: P, Q
            steps: begin
                E + S ⇌ E(S)
                E + P ⇌ E(P)
                E + Q <--> E(Q)
                E(P) + Q <--> E(P, Q)
                E(Q) + P ⇌ E(P, Q)
                E(S) <--> E(P, Q)
            end
        end)
    ])
    @test length(children) == 6
    @test Set(children) == Set(expected)
end

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
    # SEED: uni-uni allosteric with both binding groups :OnlyA (a
    # balanced substrate/product pair, so the mechanism is Haldane-valid)
    # and :EqualAI catalysis. :OnlyA groups live in the R-state only; the
    # T-state contributes no kf_T/kr_T after RE→SS.
    em_seed = @allosteric_mechanism begin
        substrates: S
        products: P
        catalytic_multiplicity: 2
        catalytic_steps: begin
            E + P ⇌ E(P)      :: OnlyA
            E + S ⇌ E(S)      :: OnlyA
            E(S) <--> E(P)    :: EqualAI
        end
    end
    am = EnzymeRates.AllostericMechanism(em_seed)
    EnzymeRates._assert_mechanism_invariants(am)

    result = EnzymeRates._expand_re_to_ss(am)

    # 1. count: 2 all-RE binding groups (both :OnlyA). Iso SS. → 2 variants.
    @test length(result) == 2

    # 2. Δ params. Connected-component pruning populates the inactive state
    # from free E over ALL surviving steps, so because catalysis is :EqualAI
    # the base mechanism's Q_I ALREADY carries E(S) and E(P) through the
    # inactive catalytic step — even though both bindings are :OnlyA and
    # cannot bind directly in the T-state. Those forms are therefore present
    # in the inactive state of the base AND of every variant, so flipping a
    # binding group RE→SS adds only that group's own SS rate constant (+1) —
    # it does not introduce a new inactive-state form. Both variants → Δ=1.
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
    # _expand_re_to_ss keeps inhibitor bindings RE-only, so the
    # dead-end-inhibitor-S kinetic group never flips, even though the
    # substrate-S and inhibitor-S groups are otherwise independent.
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
    # product-binding (RE), iso (SS), dead-end-S__reg-binding (RE, but
    # binds a Regulator so it is skipped). → 2 RE groups → 2 variants.
    @test length(result) == 2
    for r in result
        @test r isa EnzymeRates.Mechanism
        EnzymeRates._assert_mechanism_invariants(r)
        @test EnzymeRates.compile_mechanism(r) isa EnzymeMechanism
    end

    # 2. property-style: each variant flips exactly one initial RE
    # group to all-SS, and the flipped group covers 2 distinct values
    # (substrate-binding, product-binding) across the 2 variants; the
    # dead-end-inhibitor-S group never appears among them.
    flipped_groups = Int[]
    for r in result
        for (gi, (old_grp, new_grp)) in enumerate(zip(m.steps, r.steps))
            if all(EnzymeRates.is_equilibrium, old_grp) &&
               !any(EnzymeRates.is_equilibrium, new_grp)
                push!(flipped_groups, gi)
            end
        end
    end
    @test length(unique(flipped_groups)) == 2

    # 3. preservation: reaction unchanged.
    for r in result
        @test EnzymeRates.reaction(r) == EnzymeRates.reaction(m)
    end
end

@testset "Mechanism — Allosteric substrate-as-dead-end-I overlap" begin
    # AllostericMechanism counterpart of the substrate-as-I overlap.
    # Same count derivation: the dead-end-inhibitor-S group binds a
    # Regulator and stays RE, leaving 2 RE groups → 2 variants. Each
    # flip should preserve cat_allo_states tags exactly.
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

    # 1. count: 2 RE groups (dead-end-inhibitor-S excluded). → 2 variants.
    @test length(result) == 2
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

@testset "_expand_re_to_ss keeps inhibitor bindings RE" begin
    rxn = @enzyme_reaction begin
        substrates: S[C]
        products: P[C]
        dead_end_inhibitors: S
    end
    m = first(EnzymeRates._expand_add_dead_end_regulator(
              first(EnzymeRates.init_mechanisms(rxn)), rxn))
    for r in EnzymeRates._expand_re_to_ss(m)
        for grp in EnzymeRates.steps(r), s in grp
            if EnzymeRates.bound_metabolite(s) isa EnzymeRates.Regulator
                @test EnzymeRates.is_equilibrium(s)   # inhibitor binding stays RE
            end
        end
    end
end

@testset "_expand_re_to_ss flips inhibitor-bound mirrors together" begin
    # A catalytic binding and its inhibitor-bound mirror (the same binding
    # with every inhibitor stripped) sit in separate all-RE kinetic groups
    # on a closed inhibitor square: A binding to E, and A binding to the
    # inhibitor-bound form E(I), with I binding both E and E(A). Neither
    # cuts a rapid-equilibrium segment alone, because the other supplies an
    # RE route around the square, so the minimal-set search emits them only
    # together: they flip to SS as a pair, never one without the other.
    # Built via the macro; a separate catalytic inhibitor I keeps the
    # substrate/inhibitor roles unambiguous. Binding is ordered (A then B)
    # so the pair flip stays hyperbolic; a conformational parent drops a
    # child whose catalytic scheme carries a concentration power, and the
    # random-order pair flip would carry B².
    m = EnzymeRates.AllostericMechanism(EnzymeRates.@allosteric_mechanism begin
        substrates: A, B
        products: P
        catalytic_inhibitors: I
        catalytic_steps: begin
            E + A ⇌ E(A)          :: EqualAI
            E(A) + B ⇌ E(A, B)    :: EqualAI
            E(A, B) <--> E(P)      :: EqualAI
            E + P ⇌ E(P)          :: EqualAI
            E + I ⇌ E(I)          :: EqualAI
            E(A) + I ⇌ E(A, I)    :: EqualAI
            E(I) + A ⇌ E(A, I)    :: EqualAI
        end
    end)
    # Non-vacuity: the A-binding group and its inhibitor-bound mirror are
    # separate all-RE groups, each bridged by the other, so they can only
    # flip as a pair — and some child does flip the pair.
    kids = EnzymeRates._expand_re_to_ss(m)
    @test !isempty(kids)
    core(s) = begin
        strip(sp) = EnzymeRates.Species(
            EnzymeRates.Metabolite[b for b in EnzymeRates.bound(sp)
                                   if !(b isa EnzymeRates.Regulator)],
            EnzymeRates.conformation(sp), EnzymeRates.residual(sp))
        (strip(EnzymeRates.from_species(s)), strip(EnzymeRates.to_species(s)),
         EnzymeRates.bound_metabolite(s))
    end
    # The base E + A ⇌ E(A) and its mirror E(I) + A ⇌ E(A, I) share a core.
    base = only(s for grp in EnzymeRates.steps(m) for s in grp
                if EnzymeRates.name(EnzymeRates.from_species(s)) == :E &&
                   EnzymeRates.name(EnzymeRates.to_species(s)) == :EA)
    pair_ss = [EnzymeRates.Step(EnzymeRates.from_species(s),
                                EnzymeRates.to_species(s),
                                EnzymeRates.bound_metabolite(s), false)
               for grp in EnzymeRates.steps(m) for s in grp
               if core(s) == core(base)]
    @test length(pair_ss) == 2
    ss_steps(c) = Set(s for grp in EnzymeRates.steps(c) for s in grp
                      if !EnzymeRates.is_equilibrium(s))
    @test any(c -> all(in(ss_steps(c)), pair_ss), kids)
    # Mirror-lock: no re_to_ss variant leaves a mirror RE while its base is SS.
    for r in kids
        status = Dict{Any, Bool}()
        for grp in EnzymeRates.steps(r)
            allss = !any(EnzymeRates.is_equilibrium, grp)
            for s in grp
                c = core(s)
                haskey(status, c) ? (@test status[c] == allss) :
                                    (status[c] = allss)
            end
        end
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

@testset "AllostericMechanism — :NonequalAI: 2 variants, tags preserved" begin
    # SEED: uni-uni allosteric with every group :NonequalAI — a
    # distinguishable mechanism (an all-:EqualAI seed is a model-space
    # no-op that is never enumerated). Built via @allosteric_mechanism,
    # then lifted to the decomposed form the expansion moves operate on.
    am = EnzymeRates.AllostericMechanism(@allosteric_mechanism begin
        substrates: S
        products: P
        catalytic_multiplicity: 2
        catalytic_steps: begin
            E + S ⇌ E(S)    :: NonequalAI
            E + P ⇌ E(P)    :: NonequalAI
            E(S) <--> E(P)  :: NonequalAI
        end
    end)

    result = EnzymeRates._expand_re_to_ss(am)

    # 1. count: 2 all-RE binding groups (P-, S-binding). Iso is SS. → 2.
    @test length(result) == 2

    # 2. Δ params: :NonequalAI RE→SS adds 2 fitted params per variant (the
    # A- and I-state binding K each split into (kf, kr)). Measured against
    # the actual compiled count.
    base_fitted = length(EnzymeRates.fitted_params(
        EnzymeRates.compile_mechanism(am)))
    deltas = sort([length(EnzymeRates.fitted_params(
        EnzymeRates.compile_mechanism(r))) - base_fitted for r in result])
    @test deltas == [2, 2]

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
    @test !isempty(result)
    base = EnzymeRates._independent_param_count(m)
    for r in result
        @test EnzymeRates._independent_param_count(r) > base
        @test r isa EnzymeRates.Mechanism
        EnzymeRates._assert_mechanism_invariants(r)
        @test EnzymeRates.compile_mechanism(r) isa EnzymeMechanism
        @test EnzymeRates.n_steps(r) == EnzymeRates.n_steps(m)
        @test EnzymeRates.reaction(r) == EnzymeRates.reaction(m)
    end
    # The SS A-binding group has no equilibrium constant to tie, so its
    # context split by B gains alone and is emitted as a single split.
    a_split = [r for r in result
               if length(EnzymeRates.steps(r)) == length(m.steps) + 1 &&
               any(grp -> length(grp) == 1 && !EnzymeRates.is_equilibrium(only(grp)) &&
                   EnzymeRates.name(EnzymeRates.bound_metabolite(only(grp))) == :A,
                   EnzymeRates.steps(r))]
    @test length(a_split) == 1
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

    # 1. every emitted child raises the independent-parameter count.
    @test !isempty(result)
    @test all(r -> EnzymeRates._independent_param_count(r) >
                   EnzymeRates._independent_param_count(am), result)

    # 2. the :NonequalAI A-binding group is SS, so it has no equilibrium
    # constant to tie and its context split gains on its own: a child
    # that splits that group and nothing else is emitted.
    a_group = only(grp for grp in am.cat_steps
                   if length(grp) == 2 && EnzymeRates.name(
                       EnzymeRates.bound_metabolite(first(grp))) == :A)
    @test any(r -> length(EnzymeRates.steps(r)) == length(am.cat_steps) + 1 &&
                   count(grp -> Set(grp) ⊆ Set(a_group),
                         EnzymeRates.steps(r)) == 2, result)

    # 3. compilability
    for r in result
        @test r isa EnzymeRates.AllostericMechanism
        EnzymeRates._assert_mechanism_invariants(r)
        @test EnzymeRates.compile_mechanism(r) isa AllostericEnzymeMechanism
    end

    # 4. tag inheritance: each result group inherits the tag of the parent
    # group whose steps contain it. A child subdivides one or more groups;
    # both halves of each carry that group's tag, and every other group
    # keeps its own. Group ORDER is canonical (not source-preserved), so match each
    # result group to its parent by step content rather than by position.
    for r in result
        @test length(r.cat_allo_states) == length(r.cat_steps)
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

@testset "Mechanism — bi-bi random: single splits tied, a pair frees a constant" begin
    # SEED: bi-bi random with 4 multi-step kinetic groups (A, B, P, Q)
    # forming a single closed thermodynamic cycle. Dividing any one group
    # produces a binding-K that the cycle ties straight back to an existing
    # parameter, so no single division gains; the cycle is broken only by
    # dividing a second group on it as well.
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
    result = EnzymeRates._expand_split_kinetic_group(m)
    @test !isempty(result)
    base = EnzymeRates._independent_param_count(m)
    for r in result
        @test EnzymeRates._independent_param_count(r) > base
        # Every child splits at least two groups: no single split gains here.
        @test length(EnzymeRates.steps(r)) >= length(m.steps) + 2
    end
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

@testset "AllostericMechanism — bi-bi RE groups: every emitted split gains" begin
    # SEED: bi-bi allosteric with mixed tags (:NonequalAI, :EqualAI),
    # but both multi-step groups (A, B) are RE bindings, each closing
    # its own per-conformer thermodynamic cycle, so a division of one
    # group alone may be tied back by that cycle. Tag inheritance is
    # covered by "AllostericMechanism — SS multi-step :NonequalAI
    # split" above.
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
    result = EnzymeRates._expand_split_kinetic_group(am)
    @test !isempty(result)
    @test all(r -> EnzymeRates._independent_param_count(r) >
                   EnzymeRates._independent_param_count(am), result)
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

@testset "Mechanism — Bi-bi sequential: binding-group :OnlyA variants only" begin
    # 5 kinetic groups: 4 bindings (A, B substrate-side; P, Q
    # product-side) + 1 chemical step (E(A,B) <--> E(P,Q)). The
    # catalytic group is never a bare primary :OnlyA (that V-type needs
    # a regulator, none declared). Each binding's one-sided :OnlyA is
    # instead closed over its minimal Haldane completions.
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

    # 1. count: every non-empty subset of the 4 binding groups is set
    # :OnlyA, each with the chemical step also :OnlyA (a catalytically-dead
    # inactive conformation). 2^4 - 1 = 15 subsets, all Wegscheider-valid.
    n_groups = length(m.steps)
    n_cat = count(g -> _is_catalytic_group(m, g), 1:n_groups)
    @test n_cat == 1
    @test length(result) == 15

    # 2. Δ params: every variant is +1 — the conformational constant L.
    # Each :OnlyA binding is K-type (the inactive conformation cannot bind
    # that ligand, adding no inactive binding constant), and the dead
    # inactive conformation runs no chemistry, so it carries no free
    # catalytic constant. All 15 variants are Δ=1.
    base_fitted = length(EnzymeRates.fitted_params(
        EnzymeRates.compile_mechanism(m)))
    deltas = sort([length(EnzymeRates.fitted_params(
        EnzymeRates.compile_mechanism(r))) - base_fitted for r in result])
    @test deltas == fill(1, 15)

    # 3. compilability — must produce AllostericEnzymeMechanism, and no
    # result is the all-:EqualAI baseline.
    for r in result
        @test r isa EnzymeRates.AllostericMechanism
        EnzymeRates._assert_mechanism_invariants(r)
        @test EnzymeRates.compile_mechanism(r) isa AllostericEnzymeMechanism
        @test !all(==(:EqualAI), EnzymeRates.cat_allo_states(r))
    end
end

@testset "Mechanism — Bi-bi ping-pong: binding-group :OnlyA variants only" begin
    # SEED: bi-bi ping-pong topology, 6 kinetic groups: 4 bindings
    # (A, B substrate-side; P, Q product-side) + 2 chemical steps (the
    # two half-reaction interconversions E(A)<-->Estar(A,P) and
    # Estar(B)⇌E(Q)). Neither catalytic group is a bare primary :OnlyA
    # (that V-type needs a regulator, none declared); each binding's
    # one-sided :OnlyA is closed over its minimal Haldane completions.
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

    # 1. count: every non-empty subset of the 4 binding groups is set
    # :OnlyA, each with BOTH chemical steps also :OnlyA (a catalytically-dead
    # inactive conformation). 2^4 - 1 = 15 subsets, all Wegscheider-valid.
    n_groups = length(m.steps)
    n_cat = count(g -> _is_catalytic_group(m, g), 1:n_groups)
    @test n_cat == 2
    @test length(result) == 15

    # 2. Δ params: every variant is +1 — the conformational constant L.
    # Each :OnlyA binding is K-type (no inactive binding constant), and the
    # dead inactive conformation runs no chemistry, so it carries no free
    # catalytic constant. All 15 variants are Δ=1.
    base_fitted = length(EnzymeRates.fitted_params(
        EnzymeRates.compile_mechanism(m)))
    deltas = sort([length(EnzymeRates.fitted_params(
        EnzymeRates.compile_mechanism(r))) - base_fitted for r in result])
    @test deltas == fill(1, 15)

    # 3. compilability, and no result is the all-:EqualAI baseline.
    for r in result
        @test r isa EnzymeRates.AllostericMechanism
        EnzymeRates._assert_mechanism_invariants(r)
        @test EnzymeRates.compile_mechanism(r) isa AllostericEnzymeMechanism
        @test !all(==(:EqualAI), EnzymeRates.cat_allo_states(r))
    end
end

@testset "Mechanism — uni-uni: binding-group :OnlyA variants only" begin
    # SEED: uni-uni init Mechanism, no declared regulator. 3 kinetic
    # groups: 2 bindings (S, P) + 1 catalytic. The catalytic group is
    # never a bare primary :OnlyA (that V-type needs a regulator, none
    # declared); each binding's one-sided :OnlyA is closed over its
    # minimal Haldane completions.
    rxn = @enzyme_reaction begin
        substrates: S[C]
        products: P[C]
        oligomeric_state: 2
    end
    m = first(EnzymeRates.init_mechanisms(rxn))

    result = EnzymeRates._expand_to_allosteric(m, rxn)

    # 1. count: every non-empty subset of the 2 binding groups (S, P) is
    # set :OnlyA, each with the chemical step also :OnlyA (dead inactive):
    # {S}, {P}, {S, P} → 3 variants.
    n_groups = length(m.steps)
    n_cat = count(g -> _is_catalytic_group(m, g), 1:n_groups)
    @test n_cat == 1
    @test length(result) == 3

    # 2. each result is an AllostericMechanism with correct
    # multiplicity and empty regulatory_sites (K-type: bare).
    for r in result
        @test r isa EnzymeRates.AllostericMechanism
        @test r.catalytic_multiplicity == 2
        @test isempty(r.regulatory_sites)
        @test EnzymeRates.reaction(r) == EnzymeRates.reaction(m)
        @test length(r.cat_allo_states) == n_groups
    end

    # 3. no result is the all-:EqualAI baseline.
    @test all(r -> !all(==(:EqualAI), r.cat_allo_states), result)

    # 4. every variant is dead-inactive: the chemical (iso) group is
    # :OnlyA, plus a non-empty subset of the binding groups — so the
    # :OnlyA count is 2 ({one binding} + chem) or 3 ({S, P} + chem), and
    # at least one binding group is :OnlyA.
    for r in result
        iso_gs = [g for g in 1:length(r.cat_allo_states)
                  if _is_catalytic_group(m, g)]
        @test all(r.cat_allo_states[g] == :OnlyA for g in iso_gs)
        @test count(==(:OnlyA), r.cat_allo_states) in (2, 3)
        onlya_gs = [g for g in 1:length(r.cat_allo_states)
                    if r.cat_allo_states[g] == :OnlyA]
        @test any(g -> !_is_catalytic_group(m, g), onlya_gs)
    end

    # 5. compilability: each AllostericMechanism compiles to an
    # AllostericEnzymeMechanism.
    for r in result
        @test EnzymeRates.compile_mechanism(r) isa
            EnzymeRates.AllostericEnzymeMechanism
    end
end

@testset "Mechanism — enumerates all allowed multiplicities" begin
    # Multi-valued allowed_catalytic_multiplicities → the binding-group
    # :OnlyA variant set is emitted once per multiplicity.
    rxn = @enzyme_reaction begin
        substrates: S[C]
        products: P[C]
        allowed_catalytic_multiplicities: (2, 4)
    end
    m = first(EnzymeRates.init_mechanisms(rxn))
    allo = EnzymeRates._expand_to_allosteric(m, rxn)
    mults = Set(EnzymeRates.catalytic_multiplicity(am) for am in allo)
    @test mults == Set([2, 4])

    # No regulator declared → the binding-completion variant set (see the
    # uni-uni case: 3 variants) is emitted once per multiplicity.
    n_groups = length(EnzymeRates.steps(m))
    n_cat = count(g -> _is_catalytic_group(m, g), 1:n_groups)
    @test n_cat == 1
    @test length(allo) == 2 * 3

    # Each multiplicity carries the same dead-inactive allo-state set: no
    # all-:EqualAI baseline, and 3 variants ({S}, {P}, {S, P} bindings
    # :OnlyA, each with the chemical step :OnlyA) — :OnlyA counts 2, 2, 3.
    for cn in (2, 4)
        cn_variants = filter(
            am -> EnzymeRates.catalytic_multiplicity(am) == cn, allo)
        @test count(
            am -> all(==(:EqualAI), EnzymeRates.cat_allo_states(am)),
            cn_variants) == 0
        @test count(
            am -> count(==(:OnlyA), EnzymeRates.cat_allo_states(am)) == 2,
            cn_variants) == 2
        @test count(
            am -> count(==(:OnlyA), EnzymeRates.cat_allo_states(am)) == 3,
            cn_variants) == 1
    end

    # Single-valued case unchanged (regression guard).
    rxn1 = @enzyme_reaction begin
        substrates: S[C]
        products: P[C]
        oligomeric_state: 2
    end
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

@testset "no all-:EqualAI baseline is ever emitted" begin
    # Structural guard: across a regulator-free and a
    # regulator-declaring reaction, no emitted mechanism has every
    # group tagged :EqualAI — the conformational constant L would
    # cancel and be unobservable.
    rxn_no_reg = @enzyme_reaction begin
        substrates: S[C]
        products: P[C]
        oligomeric_state: 4
    end
    for rxn in (rxn_no_reg, uni_uni_allo_reg)
        base = first(EnzymeRates.init_mechanisms(rxn))
        allo = EnzymeRates._expand_to_allosteric(base, rxn)
        @test !isempty(allo)
        @test all(
            am -> !all(==(:EqualAI), EnzymeRates.cat_allo_states(am)), allo)
    end
end

@testset "catalysis-:OnlyA is V-type only (needs a regulator)" begin
    is_bare_cat_onlyA(base, am) = begin
        g = findfirst(==(:OnlyA), EnzymeRates.cat_allo_states(am))
        g !== nothing && _is_catalytic_group(base, g) &&
            isempty(EnzymeRates.regulatory_sites(am))
    end

    # Regulator-free reaction: no bare catalysis-:OnlyA is ever emitted.
    rxn0 = @enzyme_reaction begin
        substrates: S[C]
        products: P[C]
        oligomeric_state: 2
    end
    base0 = first(EnzymeRates.init_mechanisms(rxn0))
    allo0 = EnzymeRates._expand_to_allosteric(base0, rxn0)
    @test !any(am -> is_bare_cat_onlyA(base0, am), allo0)

    # Regulator-declaring reaction: a V-type (catalysis-:OnlyA paired
    # with the declared regulator) IS reachable, and still no bare one
    # is emitted.
    base1 = first(EnzymeRates.init_mechanisms(uni_uni_allo_reg))
    allo1 = EnzymeRates._expand_to_allosteric(base1, uni_uni_allo_reg)
    @test !any(am -> is_bare_cat_onlyA(base1, am), allo1)
    @test any(am -> begin
        g = findfirst(==(:OnlyA), EnzymeRates.cat_allo_states(am))
        g !== nothing && _is_catalytic_group(base1, g) &&
            !isempty(EnzymeRates.regulatory_sites(am))
    end, allo1)
end

@testset "V-type param counts + competitive-inhibitor exclusion" begin
    # uni_uni_allo_2reg declares two allosteric regulators (R1, R2). From a
    # uni-uni seed (RE binding, SS catalytic), _expand_to_allosteric emits:
    #   * 3 K-type variants (binding-group :OnlyA, no regulatory site) at
    #     +1 param each — the conformational constant L; each is a binding
    #     promotion closed with its Haldane completion (a second :OnlyA
    #     binding, or the dropped chemical step), which adds no free param;
    #   * 4 V-type variants (the catalytic group :OnlyA paired with one
    #     regulator at a site) at +2 params each — L plus the regulator's K
    #     — one per (regulator, tag), regulator ∈ {R1, R2}, tag ∈
    #     {:OnlyA, :OnlyI}.
    # The uni-uni seed (RE binding, SS catalysis) carrying
    # uni_uni_allo_2reg's two declared regulators, built explicitly —
    # first(init_mechanisms(...)) is not guaranteed to stay this mechanism.
    seed = EnzymeRates.Mechanism(uni_uni_allo_2reg, EnzymeRates.steps(
        EnzymeRates.Mechanism(@enzyme_mechanism begin
            substrates: S
            products: P
            steps: begin
                E + S ⇌ E(S)
                E(S) <--> E(P)
                E + P ⇌ E(P)
            end
        end)))
    base = length(EnzymeRates.fitted_params(
        EnzymeRates.compile_mechanism(seed)))
    Δ(am) = length(EnzymeRates.fitted_params(
        EnzymeRates.compile_mechanism(am))) - base
    outs = EnzymeRates._expand_to_allosteric(seed, uni_uni_allo_2reg)
    vtypes = filter(am -> !isempty(EnzymeRates.regulatory_sites(am)), outs)
    ktypes = filter(am -> isempty(EnzymeRates.regulatory_sites(am)), outs)

    @test length(vtypes) == 4
    @test all(am -> Δ(am) == 2, vtypes)
    @test length(ktypes) == 3
    @test all(am -> Δ(am) == 1, ktypes)

    # The four V-types cover R1/R2 × :OnlyA/:OnlyI, each a single regulator
    # at a single site.
    combos = Set((only([EnzymeRates.name(l)
                        for s in EnzymeRates.regulatory_sites(am)
                        for l in EnzymeRates.ligands(s)]),
                  only([t for s in EnzymeRates.regulatory_sites(am)
                        for t in EnzymeRates.allo_states(s)]))
                 for am in vtypes)
    @test combos == Set([(:R1, :OnlyA), (:R1, :OnlyI),
                         (:R2, :OnlyA), (:R2, :OnlyI)])

    # A regulator already bound as a competitive inhibitor in the seed is
    # NOT re-used for a V-type — only the other (free) allosteric regulator
    # is. In uni_uni_reg_and_inhibitor, R1 is the competitive inhibitor and
    # R2 the free allosteric regulator.
    pool = unique!(EnzymeRates.expand_mechanisms(
        EnzymeRates.Mechanism[EnzymeRates.init_mechanisms(
            uni_uni_reg_and_inhibitor)...], uni_uni_reg_and_inhibitor))
    seed_ci = first(filter(pool) do m
        m isa EnzymeRates.Mechanism && any(
            (bm = EnzymeRates.bound_metabolite(s);
             bm isa EnzymeRates.Regulator && EnzymeRates.name(bm) == :R1)
            for g in EnzymeRates.steps(m) for s in g)
    end)
    vt_ci = filter(am -> !isempty(EnzymeRates.regulatory_sites(am)),
                   EnzymeRates._expand_to_allosteric(
                       seed_ci, uni_uni_reg_and_inhibitor))
    site_regs = Set(EnzymeRates.name(l) for am in vt_ci
                    for s in EnzymeRates.regulatory_sites(am)
                    for l in EnzymeRates.ligands(s))
    @test !isempty(vt_ci)
    @test site_regs == Set([:R2])
end

@testset "distinguishability invariant on every emitted mechanism" begin
    # Every mechanism _expand_to_allosteric returns is distinguishable
    # from a simpler mechanism: not all-:EqualAI, and either a binding
    # group carries a non-:EqualAI tag or a regulatory site is present.
    scenarios = [
        (first(EnzymeRates.init_mechanisms(uni_uni_allo)), uni_uni_allo),
        (first(EnzymeRates.init_mechanisms(uni_uni_allo_reg)), uni_uni_allo_reg),
    ]
    for (base, rxn) in scenarios
        for am in EnzymeRates._expand_to_allosteric(base, rxn)
            states = EnzymeRates.cat_allo_states(am)
            @test !all(==(:EqualAI), states)
            binding_gs = [g for g in 1:length(states)
                          if !_is_catalytic_group(base, g)]
            binding_non_equal = any(g -> states[g] != :EqualAI, binding_gs)
            @test binding_non_equal ||
                !isempty(EnzymeRates.regulatory_sites(am))
        end
    end
end

@testset "Mechanism — non-hyperbolic catalytic scheme: no children" begin
    # Uni-bi with random SS product release carries P² and Q² in its own
    # denominator; a conformational equilibrium on top would stack a second
    # source of concentration powers, so the promotion emits nothing.
    m = EnzymeRates.Mechanism(@enzyme_mechanism begin
        substrates: S
        products: P, Q
        steps: begin
            E + S ⇌ E(S)
            E(S) <--> E(P, Q)
            (E + P <--> E(P), E(Q) + P <--> E(P, Q))
            (E + Q <--> E(Q), E(P) + Q <--> E(P, Q))
        end
    end)
    rxn = @enzyme_reaction begin
        substrates: S[AB]
        products: P[A], Q[B]
        oligomeric_state: 2
    end
    @test isempty(EnzymeRates._expand_to_allosteric(m, rxn))
end

@testset "Mechanism — substrate dead-end keeps the scheme hyperbolic" begin
    # Ordered SS bi-bi with A as a dead end on E(Q): the dead-end's A² term is
    # substrate inhibition, not random-order steady state, so the promotion
    # proceeds. Five binding groups (A, B, Q, P, the dead end), each subset
    # :OnlyA with the chemistry :OnlyA: 2^5 - 1 = 31 K-type children.
    m = EnzymeRates.Mechanism(@enzyme_mechanism begin
        substrates: A, B
        products: P, Q
        steps: begin
            E + A <--> E(A)
            E(A) + B <--> E(A, B)
            E + Q <--> E(Q)
            E(Q) + P <--> E(P, Q)
            E(A, B) <--> E(P, Q)
            E(Q) + A ⇌ E(A, Q)
        end
    end)
    rxn = @enzyme_reaction begin
        substrates: A[C], B[N]
        products: P[C], Q[N]
        oligomeric_state: 2
    end
    children = EnzymeRates._expand_to_allosteric(m, rxn)
    @test length(children) == 31
    @test all(EnzymeRates._hyperbolic_catalysis, children)
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

@testset "AllostericMechanism — competitive-only name: no allo variant" begin
    # Build a uni-uni AllostericMechanism that already has :I bound as a
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

@testset "AllostericMechanism — Two regulators with site options: count = 6" begin
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

    # 1. count: 3 non-:EqualAI tag flavors × 2 site options (new + R1's
    # existing) = 6, minus the redundant :OnlyI-onto-R1's-:OnlyA-site
    # append (disjoint conformations) = 5. Plus 1 :EqualAI-at-existing
    # variant (gated on R1 being non-:EqualAI). → 6.
    @test length(result) == 6

    # 2. Δ params: four variants add one parameter (:OnlyA new/existing,
    # :OnlyI new, and the :EqualAI-at-existing one), two :NonequalAI
    # variants add two (K_R2 + K_R2_T). Sorted [1, 1, 1, 1, 2, 2].
    # Measured against the actual compiled fitted count.
    base_fitted = length(EnzymeRates.fitted_params(
        EnzymeRates.compile_mechanism(am)))
    deltas = sort([length(EnzymeRates.fitted_params(
        EnzymeRates.compile_mechanism(r))) - base_fitted for r in result])
    @test deltas == [1, 1, 1, 1, 2, 2]

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

    # 1. count: 6 (same derivation as the 6-variant case — the redundant
    # :OnlyI-onto-R1's-:OnlyA-site append is skipped).
    @test length(result) == 6

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

    # 1. count: 3 non-:EqualAI × 2 sites + 1 :EqualAI-at-existing = 7,
    # minus the redundant :OnlyA-onto-R1's-:OnlyI-site append (disjoint
    # conformations) = 6.
    @test length(result) == 6

    # 2. Δ params: same multiset as the :OnlyA seed — four +1 variants
    # and two :NonequalAI +2 variants. Sorted [1, 1, 1, 1, 2, 2].
    # Measured against the actual compiled fitted count.
    base_fitted = length(EnzymeRates.fitted_params(
        EnzymeRates.compile_mechanism(am)))
    deltas = sort([length(EnzymeRates.fitted_params(
        EnzymeRates.compile_mechanism(r))) - base_fitted for r in result])
    @test deltas == [1, 1, 1, 1, 2, 2]

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

    # 1. count: 3 non-:EqualAI × 2 sites + 1 :EqualAI-at-existing = 7,
    # minus the redundant :OnlyI-onto-R1's-:OnlyA-site append = 6.
    @test length(result) == 6

    # 4. property-style: separate new-site (#sites grows to 2) and
    # existing-site (#sites stays at 1) placements.
    new_site_variants = filter(r -> length(r.regulatory_sites) == 2, result)
    existing_site_variants =
        filter(r -> length(r.regulatory_sites) == 1, result)
    @test length(new_site_variants) == 3   # 3 non-:EqualAI × new site
    # :OnlyA + :NonequalAI appends onto R1's :OnlyA site + 1 :EqualAI;
    # the :OnlyI append is skipped as redundant (disjoint conformations).
    @test length(existing_site_variants) == 3
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

@testset "AllostericMechanism — dual-role name gains an allosteric site" begin
    # A reaction may declare one name (ATP) in BOTH roles. Starting from an
    # allosteric mechanism that already binds ATP as a competitive dead-end
    # inhibitor, `_expand_add_allosteric_regulator` now also adds ATP at a
    # regulatory site — the two roles carry distinct parameter names.
    rxn = @enzyme_reaction begin
        substrates: S[C]
        products: P[C]
        competitive_inhibitors: ATP
        allosteric_regulators: ATP(1)
    end
    base = first(EnzymeRates.init_mechanisms(rxn))
    with_de = EnzymeRates._expand_add_dead_end_regulator(base, rxn)
    @test !isempty(with_de)
    allo_mechs = EnzymeRates._expand_to_allosteric(first(with_de), rxn)
    @test !isempty(allo_mechs)
    am = first(allo_mechs)
    # am already binds ATP as a competitive dead-end inhibitor step.
    @test any(EnzymeRates.steps(am)) do g
        any(g) do s
            bm = EnzymeRates.bound_metabolite(s)
            bm isa EnzymeRates.CompetitiveInhibitor &&
                EnzymeRates.name(bm) == :ATP
        end
    end

    result = EnzymeRates._expand_add_allosteric_regulator(am, rxn)
    # ATP is now reachable at a regulatory site despite its dead-end role.
    @test !isempty(result)
    for r in result
        @test any(r.regulatory_sites) do site
            any(l -> EnzymeRates.name(l) === :ATP,
                EnzymeRates.ligands(site))
        end
        @test r isa EnzymeRates.AllostericMechanism
        EnzymeRates._assert_mechanism_invariants(r)
    end
end

@testset "AllostericMechanism — dual-role name derives distinct params" begin
    # ATP bound BOTH at a regulatory site AND as a competitive dead-end step
    # must yield two distinct constants (…ATPreg vs …ATPinh…), never a name
    # collision, and the mechanism must compile.
    rxn = @enzyme_reaction begin
        substrates: S[C]
        products: P[C]
        competitive_inhibitors: ATP
        allosteric_regulators: ATP(1)
    end
    base = first(EnzymeRates.init_mechanisms(rxn))
    with_de = first(EnzymeRates._expand_add_dead_end_regulator(base, rxn))
    am0 = first(EnzymeRates._expand_to_allosteric(with_de, rxn))
    # Force an ATP regulatory site on top of the ATP dead-end step.
    am = EnzymeRates._make_am_with_added_reg(am0, :ATP, :NonequalAI, 0)
    EnzymeRates._assert_mechanism_invariants(am)

    fp = collect(EnzymeRates.fitted_params(EnzymeRates.compile_mechanism(am)))
    @test length(fp) == length(unique(fp))          # no duplicate symbol
    atp_params = filter(s -> occursin("ATP", String(s)), fp)
    @test any(s -> occursin("ATPreg", String(s)), atp_params)
    @test any(s -> occursin("ATPinh", String(s)), atp_params)
    @test EnzymeRates.compile_mechanism(am) isa AllostericEnzymeMechanism
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
    base = first(EnzymeRates.init_mechanisms(rxn))
    cat = Symbol[:OnlyA for _ in 1:length(EnzymeRates.steps(base))]
    site_a = EnzymeRates.RegulatorySite(
        [EnzymeRates.AllostericRegulator(:A)], 4, [:OnlyA])
    am = EnzymeRates.AllostericMechanism(
        rxn, copy(EnzymeRates.steps(base)), cat, 4, [site_a])
    kids = EnzymeRates._expand_add_allosteric_regulator(am, rxn)
    st(child, sidx, lig) = begin
        s = EnzymeRates.regulatory_sites(child)[sidx]
        i = findfirst(l -> EnzymeRates.name(l) == lig, EnzymeRates.ligands(s))
        i === nothing ? nothing : EnzymeRates.allo_states(s)[i]
    end
    one_AI_site(c) = length(EnzymeRates.regulatory_sites(c)) == 1 &&
        Set(EnzymeRates.name(l) for l in EnzymeRates.ligands(
            only(EnzymeRates.regulatory_sites(c)))) == Set([:A, :I])
    # redundant: one site with A:OnlyA + I:OnlyI — must be skipped
    @test !any(c -> one_AI_site(c) && st(c,1,:A)==:OnlyA && st(c,1,:I)==:OnlyI, kids)
    # kept: two-site new-site form carrying I:OnlyI on its own site
    @test any(kids) do c
        length(EnzymeRates.regulatory_sites(c)) == 2 &&
            any(s -> [EnzymeRates.name(l) for l in EnzymeRates.ligands(s)] == [:I] &&
                     EnzymeRates.allo_states(s) == [:OnlyI],
                EnzymeRates.regulatory_sites(c))
    end
    # kept: both-OnlyA append (A:OnlyA + I:OnlyA on one site — not disjoint)
    @test any(c -> one_AI_site(c) && st(c,1,:A)==:OnlyA && st(c,1,:I)==:OnlyA, kids)
end

end

# ─── regulator-kind routing ────────────────────────────────────────────
@testset "regulator kinds route to their own moves" begin
# uni_uni_reg_and_inhibitor declares R2 as an allosteric regulator and R1
# as a competitive inhibitor. The two kinds must not cross moves: the
# allosteric-regulator move places only R2 (at a regulatory site), and the
# dead-end move binds only R1.
seed = first(EnzymeRates.init_mechanisms(uni_uni_reg_and_inhibitor))

@testset "allosteric-regulator addition places only R2" begin
    am = first(EnzymeRates._expand_to_allosteric(
        seed, uni_uni_reg_and_inhibitor))
    added = EnzymeRates._expand_add_allosteric_regulator(
        am, uni_uni_reg_and_inhibitor)
    site_regs = Set(EnzymeRates.name(l) for m in added
                    for s in EnzymeRates.regulatory_sites(m)
                    for l in EnzymeRates.ligands(s))
    @test !isempty(added)
    @test site_regs == Set([:R2])
end

@testset "dead-end expansion binds only R1" begin
    added = EnzymeRates._expand_add_dead_end_regulator(
        seed, uni_uni_reg_and_inhibitor)
    de_regs = Set(EnzymeRates.name(bm) for m in added
                  for g in EnzymeRates.steps(m) for s in g
                  for bm in (EnzymeRates.bound_metabolite(s),)
                  if bm isa EnzymeRates.Regulator)
    @test !isempty(added)
    @test de_regs == Set([:R1])
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
    m_seed = first(init_mechs)
    n_groups = length(EnzymeRates.steps(m_seed))
    am = EnzymeRates.AllostericMechanism(
        EnzymeRates.reaction(m_seed), copy(EnzymeRates.steps(m_seed)),
        Symbol[:EqualAI for _ in 1:n_groups], 2, EnzymeRates.RegulatorySite[])
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

# ─── _site_active_states ─────────────────────────────────────────────────
@testset "_site_active_states" begin
@testset "_site_active_states" begin
mk(states) = EnzymeRates.RegulatorySite(
    [EnzymeRates.AllostericRegulator(Symbol("R", i)) for i in eachindex(states)],
    4, collect(Symbol, states))
@test EnzymeRates._site_active_states(mk([:OnlyA])) == Set([:active])
@test EnzymeRates._site_active_states(mk([:OnlyI])) == Set([:inactive])
@test EnzymeRates._site_active_states(mk([:EqualAI])) == Set([:active, :inactive])
@test EnzymeRates._site_active_states(mk([:NonequalAI])) == Set([:active, :inactive])
@test EnzymeRates._site_active_states(mk([:OnlyA, :OnlyA])) == Set([:active])
@test EnzymeRates._site_active_states(mk([:OnlyA, :OnlyI])) ==
      Set([:active, :inactive])
end

@testset "_merged_site_state_assignments drop_all_keep" begin
base = [:OnlyA, :OnlyI]
keep = EnzymeRates._merged_site_state_assignments(base)
@test [:OnlyA, :OnlyI] in keep          # all-keep present by default
@test [:EqualAI, :OnlyI] in keep
@test [:OnlyA, :EqualAI] in keep
@test length(keep) == 3
dropped = EnzymeRates._merged_site_state_assignments(base; drop_all_keep=true)
@test !([:OnlyA, :OnlyI] in dropped)    # all-keep omitted
@test [:EqualAI, :OnlyI] in dropped     # antagonist retags retained
@test [:OnlyA, :EqualAI] in dropped
@test length(dropped) == 2
end
end

# ─── _expand_merge_regulatory_sites ─────────────────────────────────────
@testset "_expand_merge_regulatory_sites" begin
@testset "_expand_merge_regulatory_sites" begin
# A four-subunit reaction with a designated activator and inhibitor.
merge_rxn = @enzyme_reaction begin
    substrates: S[C]
    products: P[C]
    allosteric_regulators: A::Activator, I::Inhibitor
    oligomeric_state: 4
end
base = first(EnzymeRates.init_mechanisms(merge_rxn))
cat = Symbol[:OnlyA for _ in 1:length(EnzymeRates.steps(base))]
site_a = EnzymeRates.RegulatorySite(
    [EnzymeRates.AllostericRegulator(:A)], 4, [:OnlyA])
site_i = EnzymeRates.RegulatorySite(
    [EnzymeRates.AllostericRegulator(:I)], 4, [:OnlyI])
parent = EnzymeRates.AllostericMechanism(
    merge_rxn, copy(EnzymeRates.steps(base)), cat, 4, [site_a, site_i])

children = EnzymeRates._expand_merge_regulatory_sites(parent)

# The allo state of ligand `lig` in a child's single merged site.
merged_state(child, lig) = begin
    site = only(EnzymeRates.regulatory_sites(child))
    idx = findfirst(l -> EnzymeRates.name(l) == lig,
                    EnzymeRates.ligands(site))
    EnzymeRates.allo_states(site)[idx]
end
single_site_names(child) =
    Set(EnzymeRates.name(l)
        for l in EnzymeRates.ligands(only(EnzymeRates.regulatory_sites(child))))

@testset "disjoint OnlyA/OnlyI co-binding skipped; antagonists kept" begin
    # A (:OnlyA) acts only on the active state, I (:OnlyI) only on the
    # inactive one, so co-binding them on one site derives to the same
    # equation as separate sites — that all-keep merge is redundant and
    # skipped. The two antagonist retags stay (each :EqualAI ligand now
    # acts on both states, a genuinely distinct mechanism).
    @test all(c -> single_site_names(c) == Set([:A, :I]), children)
    states = Set((merged_state(c, :A), merged_state(c, :I)) for c in children)
    @test !((:OnlyA, :OnlyI) in states)  # redundant co-binding skipped
    @test (:EqualAI, :OnlyI) in states   # activator → antagonist
    @test (:OnlyA, :EqualAI) in states   # inhibitor → antagonist
    @test !((:EqualAI, :EqualAI) in states)  # all-EqualAI dropped
    @test length(children) == 2
end

@testset "every child is Δ0 (same fitted-param count as parent)" begin
    np_parent = length(
        EnzymeRates.fitted_params(EnzymeRates.compile_mechanism(parent)))
    for c in children
        @test length(
            EnzymeRates.fitted_params(EnzymeRates.compile_mechanism(c))) ==
              np_parent
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
    kept = EnzymeRates._filter_by_reg_type(children, merge_rxn)
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
    b_aa = first(EnzymeRates.init_mechanisms(rxn_aa))
    cat_aa = Symbol[:OnlyA for _ in 1:length(EnzymeRates.steps(b_aa))]
    s1 = EnzymeRates.RegulatorySite(
        [EnzymeRates.AllostericRegulator(:A1)], 4, [:OnlyA])
    s2 = EnzymeRates.RegulatorySite(
        [EnzymeRates.AllostericRegulator(:A2)], 4, [:OnlyA])
    p_aa = EnzymeRates.AllostericMechanism(
        rxn_aa, copy(EnzymeRates.steps(b_aa)), cat_aa, 4, [s1, s2])
    kids = EnzymeRates._expand_merge_regulatory_sites(p_aa)
    @test length(kids) == 3
    kept = EnzymeRates._filter_by_reg_type(kids, rxn_aa)
    @test length(kept) == 1
    surviving = only(kept)
    @test Set(EnzymeRates.allo_states(
        only(EnzymeRates.regulatory_sites(surviving)))) == Set([:OnlyA])
end

@testset "different merge routes to one 3-way site dedup" begin
    rxn3 = @enzyme_reaction begin
        substrates: S[C]
        products: P[C]
        allosteric_regulators: A, B, C
        oligomeric_state: 4
    end
    b3 = first(EnzymeRates.init_mechanisms(rxn3))
    cat3 = Symbol[:OnlyA for _ in 1:length(EnzymeRates.steps(b3))]
    mk1(n) = EnzymeRates.RegulatorySite(
        [EnzymeRates.AllostericRegulator(n)], 4, [:OnlyA])
    mk2(n1, n2) = EnzymeRates.RegulatorySite(
        [EnzymeRates.AllostericRegulator(n1),
         EnzymeRates.AllostericRegulator(n2)], 4, [:OnlyA, :OnlyA])
    # Route 1: {A,B} co-site + {C} single site → merge → {A,B,C}.
    p_ab_c = EnzymeRates.AllostericMechanism(
        rxn3, copy(EnzymeRates.steps(b3)), cat3, 4, [mk2(:A, :B), mk1(:C)])
    # Route 2: {A,C} co-site + {B} single site → merge → {A,C,B}.
    p_ac_b = EnzymeRates.AllostericMechanism(
        rxn3, copy(EnzymeRates.steps(b3)), cat3, 4, [mk2(:A, :C), mk1(:B)])
    threeway(kids) = only(filter(kids) do c
        sites = EnzymeRates.regulatory_sites(c)
        length(sites) == 1 &&
            length(EnzymeRates.ligands(only(sites))) == 3 &&
            all(==(:OnlyA), EnzymeRates.allo_states(only(sites)))
    end)
    c1 = threeway(EnzymeRates._expand_merge_regulatory_sites(p_ab_c))
    c2 = threeway(EnzymeRates._expand_merge_regulatory_sites(p_ac_b))
    @test c1 == c2
    @test hash(c1) == hash(c2)
end

@testset "ligand↔state pairing survives the name-sort (non-identity perm)" begin
    # Merge a 2-ligand site {A::OnlyA, C::OnlyI} with a 1-ligand site
    # {B::OnlyA}. The merged vcat order is ligands [A, C, B] / states
    # [OnlyA, OnlyI, OnlyA]; the name-sort reorders ligands to [A, B, C]
    # (perm [1, 3, 2]). C's :OnlyI must ride along to the new C slot — this
    # only holds if base_states is permuted in lockstep with the ligands.
    rxn3 = @enzyme_reaction begin
        substrates: S[C]
        products: P[C]
        allosteric_regulators: A, B, C
        oligomeric_state: 4
    end
    b3 = first(EnzymeRates.init_mechanisms(rxn3))
    cat3 = Symbol[:OnlyA for _ in 1:length(EnzymeRates.steps(b3))]
    site_ac = EnzymeRates.RegulatorySite(
        [EnzymeRates.AllostericRegulator(:A),
         EnzymeRates.AllostericRegulator(:C)], 4, [:OnlyA, :OnlyI])
    site_b = EnzymeRates.RegulatorySite(
        [EnzymeRates.AllostericRegulator(:B)], 4, [:OnlyA])
    p = EnzymeRates.AllostericMechanism(
        rxn3, copy(EnzymeRates.steps(b3)), cat3, 4, [site_ac, site_b])
    kids = EnzymeRates._expand_merge_regulatory_sites(p)
    # The co-binding child keeps every state (no :EqualAI).
    cobind = only(filter(kids) do c
        !(:EqualAI in EnzymeRates.allo_states(
            only(EnzymeRates.regulatory_sites(c))))
    end)
    @test merged_state(cobind, :A) == :OnlyA
    @test merged_state(cobind, :B) == :OnlyA
    @test merged_state(cobind, :C) == :OnlyI
end

@testset "no-op: Mechanism and single-site AllostericMechanism" begin
    m = first(EnzymeRates.init_mechanisms(uni_uni_rxn))
    @test EnzymeRates._expand_merge_regulatory_sites(m) ==
          EnzymeRates.AllostericMechanism[]
    single = EnzymeRates.AllostericMechanism(
        merge_rxn, copy(EnzymeRates.steps(base)), cat, 4, [site_a])
    @test isempty(EnzymeRates._expand_merge_regulatory_sites(single))
end
end

@testset "OnlyA+OnlyI merge is derivation-redundant (premise guard)" begin
# The reason the all-keep OnlyA/OnlyI merge is skipped: merged onto one
# site it evaluates to the same rate as on separate sites. Confirm
# numerically over random parameters and concentrations.
rxn = @enzyme_reaction begin
    substrates: S[C]
    products: P[C]
    allosteric_regulators: A::Activator, I::Inhibitor
    oligomeric_state: 4
end
base = first(EnzymeRates.init_mechanisms(rxn))
cat = Symbol[:OnlyA for _ in 1:length(EnzymeRates.steps(base))]
mkm(sites) = EnzymeRates.AllostericMechanism(
    rxn, copy(EnzymeRates.steps(base)), cat, 4, sites)
A = EnzymeRates.AllostericRegulator(:A)
I = EnzymeRates.AllostericRegulator(:I)
separate = mkm([EnzymeRates.RegulatorySite([A], 4, [:OnlyA]),
                EnzymeRates.RegulatorySite([I], 4, [:OnlyI])])
merged = mkm([EnzymeRates.RegulatorySite([A, I], 4, [:OnlyA, :OnlyI])])
fps = EnzymeRates.fitted_params(EnzymeRates.compile_mechanism(separate))
fpm = EnzymeRates.fitted_params(EnzymeRates.compile_mechanism(merged))
@test Set(fps) == Set(fpm)
Random.seed!(42)
for _ in 1:25
    vals = Dict(p => exp(randn()) for p in union(fps, fpm))
    ps = (; (p => vals[p] for p in fps)..., Keq=100.0, E_total=1.0)
    pm = (; (p => vals[p] for p in fpm)..., Keq=100.0, E_total=1.0)
    c = (; S=exp(randn()), P=exp(randn()), A=exp(randn()), I=exp(randn()))
    vs = real(EnzymeRates.rate_equation(separate, c, ps))
    vm = real(EnzymeRates.rate_equation(merged, c, pm))
    @test isapprox(vs, vm; rtol=1e-9)
end
end
end

# ─── _declared_reg_type / _state_respects_reg_type / _filter_by_reg_type ──────────
@testset "_declared_reg_type / _state_respects_reg_type / _filter_by_reg_type" begin

@testset "_declared_reg_type" begin
    rxn = @enzyme_reaction begin
        substrates: S[C]
        products: P[C]
        allosteric_regulators: A::Activator, I::Inhibitor, U
        oligomeric_state: 2
    end
    @test EnzymeRates._declared_reg_type(:A, rxn) == :activator
    @test EnzymeRates._declared_reg_type(:I, rxn) == :inhibitor
    @test EnzymeRates._declared_reg_type(:U, rxn) == :unspecified
    # No matching regulator at all → :unspecified.
    @test EnzymeRates._declared_reg_type(:Nope, rxn) == :unspecified

    # A CompetitiveInhibitor carrying a reg_type (only reachable by direct
    # construction; the DSL restricts type tags to allosteric_regulators:)
    # is not an AllostericRegulator, so its reg_type is ignored.
    rxn_ci = EnzymeRates.EnzymeReaction(
        [EnzymeRates.ReactantAtoms(EnzymeRates.Substrate(:S), [:C => 1]),
         EnzymeRates.ReactantAtoms(EnzymeRates.Product(:P), [:C => 1])],
        [EnzymeRates.RegulatorMults(
            EnzymeRates.CompetitiveInhibitor(:X), [1], :activator)],
        [1])
    @test EnzymeRates._declared_reg_type(:X, rxn_ci) == :unspecified
end

@testset "_state_respects_reg_type — :unspecified always passes" begin
    rxn = @enzyme_reaction begin
        substrates: S[C]
        products: P[C]
        allosteric_regulators: U
        oligomeric_state: 2
    end
    for state in (:OnlyA, :OnlyI, :EqualAI, :NonequalAI)
        @test EnzymeRates._state_respects_reg_type(:U, state, Symbol[], rxn)
        @test EnzymeRates._state_respects_reg_type(:U, state, [:Other], rxn)
    end
end

@testset "_state_respects_reg_type — activator" begin
    rxn = @enzyme_reaction begin
        substrates: S[C]
        products: P[C]
        allosteric_regulators: A::Activator, A2::Activator, I::Inhibitor
        oligomeric_state: 2
    end
    # never the opposite pure state
    @test !EnzymeRates._state_respects_reg_type(:A, :OnlyI, Symbol[], rxn)
    # matching pure state / NonequalAI always allowed
    @test EnzymeRates._state_respects_reg_type(:A, :OnlyA, Symbol[], rxn)
    @test EnzymeRates._state_respects_reg_type(:A, :NonequalAI, Symbol[], rxn)
    # EqualAI: no siblings → allowed
    @test EnzymeRates._state_respects_reg_type(:A, :EqualAI, Symbol[], rxn)
    # EqualAI: opposite-type sibling → allowed
    @test EnzymeRates._state_respects_reg_type(:A, :EqualAI, [:I], rxn)
    # EqualAI: unspecified sibling → allowed
    rxn_u = @enzyme_reaction begin
        substrates: S[C]
        products: P[C]
        allosteric_regulators: A::Activator, U
        oligomeric_state: 2
    end
    @test EnzymeRates._state_respects_reg_type(:A, :EqualAI, [:U], rxn_u)
    # EqualAI: same-type sibling → rejected
    @test !EnzymeRates._state_respects_reg_type(:A, :EqualAI, [:A2], rxn)
end

@testset "_state_respects_reg_type — inhibitor (symmetric)" begin
    rxn = @enzyme_reaction begin
        substrates: S[C]
        products: P[C]
        allosteric_regulators: I::Inhibitor, I2::Inhibitor, A::Activator
        oligomeric_state: 2
    end
    @test !EnzymeRates._state_respects_reg_type(:I, :OnlyA, Symbol[], rxn)
    @test EnzymeRates._state_respects_reg_type(:I, :OnlyI, Symbol[], rxn)
    @test EnzymeRates._state_respects_reg_type(:I, :NonequalAI, Symbol[], rxn)
    @test EnzymeRates._state_respects_reg_type(:I, :EqualAI, Symbol[], rxn)
    @test EnzymeRates._state_respects_reg_type(:I, :EqualAI, [:A], rxn)
    @test !EnzymeRates._state_respects_reg_type(:I, :EqualAI, [:I2], rxn)
end

@testset "_filter_by_reg_type" begin
    rxn = @enzyme_reaction begin
        substrates: S[C]
        products: P[C]
        allosteric_regulators: R::Activator
        oligomeric_state: 2
    end
    m_seed = first(EnzymeRates.init_mechanisms(rxn))
    n_groups = length(EnzymeRates.steps(m_seed))
    cat_states = Symbol[:EqualAI for _ in 1:n_groups]

    site_bad = EnzymeRates.RegulatorySite(
        [EnzymeRates.AllostericRegulator(:R)], 2, [:OnlyI])
    am_bad = EnzymeRates.AllostericMechanism(
        rxn, copy(EnzymeRates.steps(m_seed)), cat_states, 2, [site_bad])

    site_good = EnzymeRates.RegulatorySite(
        [EnzymeRates.AllostericRegulator(:R)], 2, [:OnlyA])
    am_good = EnzymeRates.AllostericMechanism(
        rxn, copy(EnzymeRates.steps(m_seed)), cat_states, 2, [site_good])

    kept = EnzymeRates._filter_by_reg_type(
        Union{EnzymeRates.Mechanism, EnzymeRates.AllostericMechanism}[
            m_seed, am_bad, am_good], rxn)
    @test m_seed in kept   # Mechanism (no sites) passes trivially
    @test am_good in kept
    @test !(am_bad in kept)
    @test length(kept) == 2
end

@testset "_filter_by_reg_type — real two-ligand site (sibling extraction)" begin
    # A single site holding two ligands exercises the sibling-extraction
    # path (`j != i`) that single-ligand fixtures never run: each ligand's
    # sibling_names come from the site itself, not a hand-fed vector.

    # DROP: two same-type activators, both :EqualAI. Each is the other's
    # same-type EqualAI sibling, so both are rejected.
    rxn_aa = @enzyme_reaction begin
        substrates: S[C]
        products: P[C]
        allosteric_regulators: A1::Activator, A2::Activator
        oligomeric_state: 2
    end
    m_aa = first(EnzymeRates.init_mechanisms(rxn_aa))
    cat_aa = Symbol[:EqualAI for _ in 1:length(EnzymeRates.steps(m_aa))]
    site_aa = EnzymeRates.RegulatorySite(
        [EnzymeRates.AllostericRegulator(:A1),
         EnzymeRates.AllostericRegulator(:A2)], 2, [:EqualAI, :EqualAI])
    am_aa = EnzymeRates.AllostericMechanism(
        rxn_aa, copy(EnzymeRates.steps(m_aa)), cat_aa, 2, [site_aa])
    @test isempty(EnzymeRates._filter_by_reg_type(
        EnzymeRates.AllostericMechanism[am_aa], rxn_aa))

    # KEEP: an :EqualAI activator beside an :OnlyI inhibitor. The
    # activator's only sibling is opposite-type, so EqualAI passes; the
    # inhibitor's :OnlyI matches its type.
    rxn_ai = @enzyme_reaction begin
        substrates: S[C]
        products: P[C]
        allosteric_regulators: A::Activator, I::Inhibitor
        oligomeric_state: 2
    end
    m_ai = first(EnzymeRates.init_mechanisms(rxn_ai))
    cat_ai = Symbol[:EqualAI for _ in 1:length(EnzymeRates.steps(m_ai))]
    site_ai = EnzymeRates.RegulatorySite(
        [EnzymeRates.AllostericRegulator(:A),
         EnzymeRates.AllostericRegulator(:I)], 2, [:EqualAI, :OnlyI])
    am_ai = EnzymeRates.AllostericMechanism(
        rxn_ai, copy(EnzymeRates.steps(m_ai)), cat_ai, 2, [site_ai])
    kept_ai = EnzymeRates._filter_by_reg_type(
        EnzymeRates.AllostericMechanism[am_ai], rxn_ai)
    @test kept_ai == [am_ai]
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
@testset "_dedup_key" begin

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

@testset "Mechanism — expansion-path overlap: dedup actually fires" begin
    # Run two rounds of expand_mechanisms on the uni-uni init seeds
    # (Mechanism path), then unique!. Assert that the flat vector shrinks,
    # proving that two different expansion paths — two moves, or one move
    # reached through two parents — produced equivalent Mechanisms.
    init_mechs = collect(EnzymeRates.init_mechanisms(uni_uni_rxn))
    pool = unique!(EnzymeRates.expand_mechanisms(init_mechs, uni_uni_rxn))
    expanded = EnzymeRates.expand_mechanisms(pool, uni_uni_rxn)
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
    # _expand_to_allosteric on this uni-uni seed (3 kinetic groups: 2
    # binding + 1 catalytic, no regulator declared) closes each binding
    # group's one-sided :OnlyA over its Haldane completion (the dropped
    # chemical step, or the opposing binding), yielding 3 distinct
    # AllostericMechanism variants (the catalytic group is never a bare
    # primary :OnlyA — that needs a regulator). Other moves do not
    # produce AllostericMechanism output from a plain Mechanism, so
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

@testset "Reg_type filter drops violating children; no-op when typeless" begin
    # True iff `m` places ligand `reg` in allo state `st` at any site.
    function has_reg_state(m, reg, st)
        m isa EnzymeRates.AllostericMechanism || return false
        for site in EnzymeRates.regulatory_sites(m)
            for (lig, s) in zip(EnzymeRates.ligands(site),
                                EnzymeRates.allo_states(site))
                EnzymeRates.name(lig) == reg && s == st && return true
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
    m_act = first(EnzymeRates.init_mechanisms(rxn_act))
    m_plain = first(EnzymeRates.init_mechanisms(rxn_plain))

    # Raw children reproduce expand_mechanisms' moves without the filter.
    raw_act = Union{EnzymeRates.Mechanism,
                    EnzymeRates.AllostericMechanism}[]
    EnzymeRates._add_expansions_mech!(raw_act, m_act, rxn_act)
    raw_plain = Union{EnzymeRates.Mechanism,
                      EnzymeRates.AllostericMechanism}[]
    EnzymeRates._add_expansions_mech!(raw_plain, m_plain, rxn_plain)

    children_act = EnzymeRates.expand_mechanisms([m_act], rxn_act)
    children_plain = EnzymeRates.expand_mechanisms([m_plain], rxn_plain)

    # A designated activator's raw expansion carries an :OnlyI R child (a
    # V-type variant); the reg_type filter drops it, leaving strictly fewer.
    @test any(m -> has_reg_state(m, :R, :OnlyI), raw_act)
    @test !any(m -> has_reg_state(m, :R, :OnlyI), children_act)
    @test children_act == EnzymeRates._filter_by_reg_type(raw_act, rxn_act)
    @test length(children_act) < length(raw_act)

    # A typeless reaction declares no type, so the filter is a no-op:
    # expand_mechanisms returns exactly the raw children, :OnlyI included.
    @test any(m -> has_reg_state(m, :R, :OnlyI), children_plain)
    @test children_plain == raw_plain
end

@testset "Merge move wired for two-site allosteric; no-op for Mechanism" begin
    # A two-regulatory-site AllostericMechanism: A::Activator on one site,
    # I::Inhibitor on the other. expand_mechanisms must now yield the Δ0
    # one-site merged children (both ligands sharing a single site).
    merge_rxn = @enzyme_reaction begin
        substrates: S[C]
        products: P[C]
        allosteric_regulators: A::Activator, I::Inhibitor
        oligomeric_state: 4
    end
    base = first(EnzymeRates.init_mechanisms(merge_rxn))
    cat = Symbol[:OnlyA for _ in 1:length(EnzymeRates.steps(base))]
    site_a = EnzymeRates.RegulatorySite(
        [EnzymeRates.AllostericRegulator(:A)], 4, [:OnlyA])
    site_i = EnzymeRates.RegulatorySite(
        [EnzymeRates.AllostericRegulator(:I)], 4, [:OnlyI])
    parent = EnzymeRates.AllostericMechanism(
        merge_rxn, copy(EnzymeRates.steps(base)), cat, 4, [site_a, site_i])

    children = EnzymeRates.expand_mechanisms([parent], merge_rxn)
    is_merged(c) =
        c isa EnzymeRates.AllostericMechanism &&
        length(EnzymeRates.regulatory_sites(c)) == 1 &&
        Set(EnzymeRates.name(l) for l in EnzymeRates.ligands(
            only(EnzymeRates.regulatory_sites(c)))) == Set([:A, :I])
    # The two Δ0 antagonist merge children survive the reg_type filter and
    # appear in the output; the redundant OnlyA/OnlyI co-binding is skipped.
    @test count(is_merged, children) == 2
    @test issubset(
        Set(EnzymeRates._expand_merge_regulatory_sites(parent)),
        Set(children))

    # A plain Mechanism has no regulatory sites, so the merge move
    # contributes nothing: expand_mechanisms equals the six non-merge
    # moves' reg-type-filtered output, unchanged by the wiring.
    m = first(EnzymeRates.init_mechanisms(uni_uni_rxn))
    raw6 = Union{EnzymeRates.Mechanism,
                 EnzymeRates.AllostericMechanism}[]
    append!(raw6, EnzymeRates._expand_re_to_ss(m))
    append!(raw6, EnzymeRates._expand_split_kinetic_group(m))
    append!(raw6, EnzymeRates._expand_add_dead_end_regulator(m, uni_uni_rxn))
    append!(raw6, EnzymeRates._expand_to_allosteric(m, uni_uni_rxn))
    append!(raw6, EnzymeRates._expand_add_allosteric_regulator(m, uni_uni_rxn))
    append!(raw6, EnzymeRates._expand_change_allo_state(m))
    baseline = EnzymeRates._filter_by_reg_type(raw6, uni_uni_rxn)
    @test EnzymeRates.expand_mechanisms([m], uni_uni_rxn) == baseline
end

@testset "expand_mechanisms reaches ≥2 distinct-metabolite catalytic :OnlyA" begin
    rxn = @enzyme_reaction begin
        substrates: A[C], B[N]
        products:   P[C1N1]
        allosteric_regulators: R
        oligomeric_state: 2
    end

    # #distinct metabolites bound by an :OnlyA catalytic group (iso → skip)
    onlya_mets(am) = Set(EnzymeRates.name(EnzymeRates.bound_metabolite(
                            EnzymeRates.rep_step(am, g)))
        for g in EnzymeRates.kinetic_groups(am)
        if EnzymeRates.cat_allo_states(am)[g] === :OnlyA &&
           !EnzymeRates.is_iso(EnzymeRates.rep_step(am, g)))

    seen = Set{UInt64}()
    frontier = Union{EnzymeRates.Mechanism, EnzymeRates.AllostericMechanism}[]
    for m in EnzymeRates.init_mechanisms(rxn)
        h = hash(m); h in seen || (push!(seen, h); push!(frontier, m))
    end
    maxdistinct = 0
    gen = 0
    while !isempty(frontier) && gen < 14
        nextf = eltype(frontier)[]
        for c in EnzymeRates.expand_mechanisms(frontier, rxn)
            h = hash(c); h in seen || (push!(seen, h); push!(nextf, c))
            c isa EnzymeRates.AllostericMechanism &&
                (maxdistinct = max(maxdistinct, length(onlya_mets(c))))
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
# uni_uni_allo_reg (not uni_uni_allo): the iso group's :OnlyA is only
# emitted paired with a regulator (V-type), so a declared regulator is
# required to reach the ":OnlyA iso group" case below.
init_mechs = EnzymeRates.init_mechanisms(uni_uni_allo_reg)
m_seed = first(init_mechs)
allo_mechs = EnzymeRates._expand_to_allosteric(m_seed, uni_uni_allo_reg)

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


# ─── Expansion-move helpers ──────────────────────────────────────────────
# Helpers and building blocks of the RE→SS flip move and the kinetic-group
# split move. Every fixture is written out with the mechanism macro so the
# mechanism under test is visible in the testset.
const _testhelper_bibi_rxn = @enzyme_reaction begin
    substrates: A[C], B[N]
    products: P[C], Q[N]
end

@testset "_independent_param_count matches fitted_params on the test specs" begin
    for spec in MECHANISM_TEST_SPECS
        m = spec.mechanism isa EnzymeRates.AllostericEnzymeMechanism ?
            EnzymeRates.AllostericMechanism(spec.mechanism) :
            EnzymeRates.Mechanism(spec.mechanism)
        @test EnzymeRates._independent_param_count(m) ==
              length(EnzymeRates.fitted_params(EnzymeRates.compile_mechanism(m)))
    end
end

@testset "_flux_carrying_groups" begin
    # Ordered uni-uni: every step is on the catalytic cycle.
    m = EnzymeRates.Mechanism(@enzyme_mechanism begin
        substrates: S
        products: P
        steps: begin
            E + S ⇌ E(S)
            E(S) <--> E(P)
            E + P ⇌ E(P)
        end
    end)
    @test all(EnzymeRates._flux_carrying_groups(m))

    # A dead-end leaf hanging off the cycle is a bridge: not flux-carrying.
    leaf = EnzymeRates.Mechanism(@enzyme_mechanism begin
        substrates: S
        products: P
        steps: begin
            E + S ⇌ E(S)
            E(S) <--> E(P)
            E + P ⇌ E(P)
            E(P) + S ⇌ E(P, S)
        end
    end)
    fc = EnzymeRates._flux_carrying_groups(leaf)
    # The leaf group is the one whose step forms the doubly-bound E(P, S).
    leaf_group = only(g for (g, grp) in enumerate(EnzymeRates.steps(leaf))
                      if any(s -> length(EnzymeRates.bound(
                                       EnzymeRates.to_species(s))) == 2, grp))
    @test !fc[leaf_group]
    @test count(!, fc) == 1

    # Ping-pong: the second chemistry step starts at rapid equilibrium and lies
    # on the catalytic cycle, so every group is flux-carrying.
    pingpong = EnzymeRates.Mechanism(@enzyme_mechanism begin
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
    end)
    @test all(EnzymeRates._flux_carrying_groups(pingpong))

    # An inhibitor bound to both E and E(S), with the mirror E(I)+S ⇌ E(I,S): the
    # inhibitor square shares the edge E→E(S) with the catalytic cycle, so the
    # cycle E→E(I)→E(I,S)→E(S)→E(P)→E passes through chemistry and every step,
    # inhibitor binding included, carries net flux. Contrast the pendant square
    # below, whose inhibitors bind only free E.
    mirror = EnzymeRates.AllostericMechanism(EnzymeRates.@allosteric_mechanism begin
        substrates: S
        products: P
        catalytic_inhibitors: I
        catalytic_steps: begin
            E + S ⇌ E(S)          :: EqualAI
            E(S) <--> E(P)        :: EqualAI
            E + P ⇌ E(P)          :: EqualAI
            E + I ⇌ E(I)          :: EqualAI
            E(I) + S ⇌ E(I, S)    :: EqualAI
            E(S) + I ⇌ E(I, S)    :: EqualAI
        end
    end)
    @test all(EnzymeRates._flux_carrying_groups(mirror))

    # Two inhibitors binding only free E form a binding-only square joined to the
    # cycle at the single form E: no cycle through it contains chemistry, so its
    # four groups are not flux-carrying.
    pendant = EnzymeRates.AllostericMechanism(EnzymeRates.@allosteric_mechanism begin
        substrates: S
        products: P
        catalytic_inhibitors: I, J
        catalytic_steps: begin
            E + S ⇌ E(S)          :: EqualAI
            E(S) <--> E(P)        :: EqualAI
            E + P ⇌ E(P)          :: EqualAI
            E + I ⇌ E(I)          :: EqualAI
            E + J ⇌ E(J)          :: EqualAI
            E(I) + J ⇌ E(I, J)    :: EqualAI
            E(J) + I ⇌ E(I, J)    :: EqualAI
        end
    end)
    pfc = EnzymeRates._flux_carrying_groups(pendant)
    inhibits(grp) = (bm = EnzymeRates.bound_metabolite(first(grp));
                     bm !== nothing && EnzymeRates.name(bm) in (:I, :J))
    inhibitor_groups = [g for (g, grp) in enumerate(EnzymeRates.steps(pendant))
                        if inhibits(grp)]
    @test length(inhibitor_groups) == 4
    @test !any(pfc[g] for g in inhibitor_groups)
    @test all(pfc[g] for g in eachindex(pfc) if !(g in inhibitor_groups))

    # Ping-pong: every group lies on the single catalytic cycle, so every
    # group is flux-carrying.
    pp = EnzymeRates.Mechanism(@enzyme_mechanism begin
        substrates: A, B
        products: P, Q
        steps: begin
            E + A ⇌ E(A)
            E(A) <--> E(P; residual = A - P)
            E(P; residual = A - P) ⇌ E(; residual = A - P) + P
            E(; residual = A - P) + B ⇌ E(B; residual = A - P)
            E(B; residual = A - P) <--> E(Q)
            E(Q) ⇌ E + Q
        end
    end)
    @test all(EnzymeRates._flux_carrying_groups(pp))
end

@testset "_hyperbolic_catalysis" begin
    # Ordered SS bi-bi: every metabolite binds on one edge and every segment is
    # a single form, so no denominator term carries a concentration twice.
    ordered = EnzymeRates.Mechanism(@enzyme_mechanism begin
        substrates: A, B
        products: P, Q
        steps: begin
            E + A <--> E(A)
            E(A) + B <--> E(A, B)
            E + Q <--> E(Q)
            E(Q) + P <--> E(P, Q)
            E(A, B) <--> E(P, Q)
        end
    end)
    @test EnzymeRates._hyperbolic_catalysis(ordered)

    # Random SS bi-bi: A binds on E → E(A) and on E(B) → E(A, B); both edges lie
    # in one tree toward E(A, B), so the denominator carries A².
    random_ss = EnzymeRates.Mechanism(@enzyme_mechanism begin
        substrates: A, B
        products: P, Q
        steps: begin
            (E + A <--> E(A), E(B) + A <--> E(A, B))
            (E + B <--> E(B), E(A) + B <--> E(A, B))
            (E + P <--> E(P), E(Q) + P <--> E(P, Q))
            (E + Q <--> E(Q), E(P) + Q <--> E(P, Q))
            E(A, B) <--> E(P, Q)
        end
    end)
    @test !EnzymeRates._hyperbolic_catalysis(random_ss)

    # Random RE bi-bi with the chemistry step as the only SS step: one segment,
    # so the equation is the rapid-equilibrium law, degree 1 throughout.
    random_re = EnzymeRates.Mechanism(@enzyme_mechanism begin
        substrates: A, B
        products: P, Q
        steps: begin
            (E + A ⇌ E(A), E(B) + A ⇌ E(A, B))
            (E + B ⇌ E(B), E(A) + B ⇌ E(A, B))
            (E + P ⇌ E(P), E(Q) + P ⇌ E(P, Q))
            (E + Q ⇌ E(Q), E(P) + Q ⇌ E(P, Q))
            E(A, B) <--> E(P, Q)
        end
    end)
    @test EnzymeRates._hyperbolic_catalysis(random_re)

    # Random RE bi-bi with the A group flipped: the segment {E(A), E(A, B)} has
    # E(A, B) carrying B beyond its bottom E(A), and the edge E(B) → E(A, B) into
    # it leaves a form carrying B, so the term rooted there carries B².
    random_flip_a = EnzymeRates.Mechanism(@enzyme_mechanism begin
        substrates: A, B
        products: P, Q
        steps: begin
            (E + A <--> E(A), E(B) + A <--> E(A, B))
            (E + B ⇌ E(B), E(A) + B ⇌ E(A, B))
            (E + P ⇌ E(P), E(Q) + P ⇌ E(P, Q))
            (E + Q ⇌ E(Q), E(P) + Q ⇌ E(P, Q))
            E(A, B) <--> E(P, Q)
        end
    end)
    @test !EnzymeRates._hyperbolic_catalysis(random_flip_a)

    # Ordered SS bi-bi with substrate A as a dead end on E(Q): the derived
    # denominator carries A² (substrate inhibition), but the dead-end group is
    # not flux-carrying and the predicate ignores it.
    dead_end = EnzymeRates.Mechanism(@enzyme_mechanism begin
        substrates: A, B
        products: P, Q
        steps: begin
            E + A <--> E(A)
            E(A) + B <--> E(A, B)
            E + Q <--> E(Q)
            E(Q) + P <--> E(P, Q)
            E(A, B) <--> E(P, Q)
            E(Q) + A ⇌ E(A, Q)
        end
    end)
    @test EnzymeRates._hyperbolic_catalysis(dead_end)

    # Uni-bi with random SS product release: P binds on E → E(P) and on
    # E(Q) → E(P, Q), both toward E(P, Q), so the denominator carries P².
    unibi_ss = EnzymeRates.Mechanism(@enzyme_mechanism begin
        substrates: S
        products: P, Q
        steps: begin
            E + S ⇌ E(S)
            E(S) <--> E(P, Q)
            (E + P <--> E(P), E(Q) + P <--> E(P, Q))
            (E + Q <--> E(Q), E(P) + Q <--> E(P, Q))
        end
    end)
    @test !EnzymeRates._hyperbolic_catalysis(unibi_ss)

    # Ping-pong with one SS chemistry step per half-reaction: two segments,
    # each metabolite carried once per term.
    pingpong = EnzymeRates.Mechanism(@enzyme_mechanism begin
        substrates: A, B
        products: P, Q
        steps: begin
            E + A ⇌ E(A)
            E(A) <--> E(P; residual = A - P)
            E(P; residual = A - P) ⇌ E(; residual = A - P) + P
            E(; residual = A - P) + B ⇌ E(B; residual = A - P)
            E(B; residual = A - P) <--> E(Q)
            E(Q) ⇌ E + Q
        end
    end)
    @test EnzymeRates._hyperbolic_catalysis(pingpong)

    # An allosteric mechanism is scored on its catalytic steps: uni-bi with the
    # P group at steady state carries Q² (the segment {E(P), E(P, Q)} plus the
    # edge E(Q) → E(P, Q)).
    allo_unibi_flip_p = EnzymeRates.AllostericMechanism(@allosteric_mechanism begin
        substrates: S
        products: P, Q
        catalytic_multiplicity: 2
        catalytic_steps: begin
            E + S ⇌ E(S)                                        :: NonequalAI
            E(S) <--> E(P, Q)                                   :: NonequalAI
            (E + P <--> E(P), E(Q) + P <--> E(P, Q))            :: NonequalAI
            (E + Q ⇌ E(Q), E(P) + Q ⇌ E(P, Q))                  :: NonequalAI
        end
    end)
    @test !EnzymeRates._hyperbolic_catalysis(allo_unibi_flip_p)

    # Ordered SS bi-bi whose A group also binds A as an abortive complex on
    # E(Q). The abortive step went to steady state with its group, but it is a
    # dead end and is ignored; counted, its edge into E(A, Q) together with
    # E → E(A) would carry A twice in one tree.
    grouped_dead_end = EnzymeRates.Mechanism(@enzyme_mechanism begin
        substrates: A, B
        products: P, Q
        steps: begin
            (E + A <--> E(A), E(Q) + A <--> E(A, Q))
            E(A) + B <--> E(A, B)
            E + Q <--> E(Q)
            E(Q) + P <--> E(P, Q)
            E(A, B) <--> E(P, Q)
        end
    end)
    @test EnzymeRates._hyperbolic_catalysis(grouped_dead_end)

    flags = vcat(EnzymeRates._flux_carrying_steps(grouped_dead_end)...)
    all_steps = vcat(EnzymeRates.steps(grouped_dead_end)...)
    @test count(!, flags) == 1
    abortive = all_steps[findfirst(!, flags)]
    abortive_bound = Symbol[EnzymeRates.name(x)
                             for x in EnzymeRates.bound(EnzymeRates.to_species(abortive))]
    @test length(abortive_bound) == 2 && :Q in abortive_bound

    # Ping-pong whose second chemistry step is at rapid equilibrium, with B also
    # bound as an abortive complex on E(Q) in the same kinetic group as its
    # catalytic binding. One segment spans both halves: E sits one B above
    # E(; residual) and E(B, Q) two above, so with the abortive form counted the
    # segment weight would carry B². The abortive step is a dead end and is
    # ignored, leaving degree 1.
    pingpong_abortive = EnzymeRates.Mechanism(@enzyme_mechanism begin
        substrates: A, B
        products: P, Q
        steps: begin
            E + A ⇌ E(A)
            E(A) <--> E(P; residual = A - P)
            E(P; residual = A - P) ⇌ E(; residual = A - P) + P
            (E(; residual = A - P) + B ⇌ E(B; residual = A - P), E(Q) + B ⇌ E(B, Q))
            E(B; residual = A - P) ⇌ E(Q)
            E(Q) ⇌ E + Q
        end
    end)
    @test EnzymeRates._hyperbolic_catalysis(pingpong_abortive)

    @test !EnzymeRates._requires_hyperbolic_catalysis(ordered)
    @test EnzymeRates._requires_hyperbolic_catalysis(allo_unibi_flip_p)
end

@testset "_hyperbolic_catalysis matches the derived denominator" begin
    # The structural predicate against the exponents of the derived denominator,
    # over every mechanism reachable from the seeds in a bounded number of
    # expansion levels. Dead-end steps are stripped before deriving, since the
    # predicate ignores them; the reactions declare no regulators, so stripping
    # never leaves the reaction naming a metabolite no step binds. For an
    # allosteric mechanism the predicate is compared with the A-state, and the
    # I-state is checked to be hyperbolic whenever the A-state is.
    _testhelper_poly_hyperbolic(p, mets) =
        all(e <= 1 for mono in keys(p) for (s, e) in mono if s in mets)
    _testhelper_mets(rxn) = vcat(
        Symbol[EnzymeRates.name(s) for s in EnzymeRates.substrates(rxn)],
        Symbol[EnzymeRates.name(p) for p in EnzymeRates.products(rxn)])
    function _testhelper_den_hyperbolic(m::EnzymeRates.Mechanism)
        rxn = EnzymeRates.reaction(m)
        subs = Symbol[EnzymeRates.name(s) for s in EnzymeRates.substrates(rxn)]
        prods = Symbol[EnzymeRates.name(p) for p in EnzymeRates.products(rxn)]
        _, den, _ = EnzymeRates._raw_symbolic_rate_polys(
            m, EnzymeRates._step_parameters(m),
            EnzymeRates._build_wegscheider_rename_map(m), subs, prods)
        _testhelper_poly_hyperbolic(den, _testhelper_mets(rxn))
    end
    function _testhelper_flux_only(m::EnzymeRates.Mechanism)
        flux = EnzymeRates._flux_carrying_steps(m)
        EnzymeRates.Mechanism(EnzymeRates.reaction(m),
            [group[f] for (group, f) in zip(EnzymeRates.steps(m), flux) if any(f)])
    end
    function _testhelper_flux_only(am::EnzymeRates.AllostericMechanism)
        flux = EnzymeRates._flux_carrying_steps(am)
        keep = [g for (g, f) in enumerate(flux) if any(f)]
        EnzymeRates._with_steps_and_cat_states(
            am, [EnzymeRates.steps(am)[g][flux[g]] for g in keep],
            EnzymeRates.cat_allo_states(am)[keep])
    end
    function _testhelper_levels(rxn, depth)
        M = Union{EnzymeRates.Mechanism, EnzymeRates.AllostericMechanism}
        level = M[m for m in EnzymeRates.init_mechanisms(rxn)]
        mechs = unique(level)
        for _ in 1:depth
            level = unique!(EnzymeRates.expand_mechanisms(level, rxn))
            append!(mechs, level)
        end
        unique!(mechs)
    end

    unibi = @enzyme_reaction begin
        substrates: S[AB]
        products: P[A], Q[B]
        oligomeric_state: 2
    end
    bibi = @enzyme_reaction begin
        substrates: A[C], B[N]
        products: P[C], Q[N]
        oligomeric_state: 2
    end
    pingpong = @enzyme_reaction begin
        substrates: A[CX], B[N]
        products: P[C], Q[NX]
        oligomeric_state: 2
    end
    n_checked = 0
    n_nonhyperbolic = 0
    # stripped mechanism => (A hyperbolic, I hyperbolic)
    derived = Dict{Any, Tuple{Bool, Bool}}()
    for (rxn, depth) in ((unibi, 2), (bibi, 1), (pingpong, 2))
        mets = _testhelper_mets(rxn)
        for m in _testhelper_levels(rxn, depth)
            EnzymeRates._eq_complexity(m) <= 337 || continue
            stripped = _testhelper_flux_only(m)
            structural = EnzymeRates._hyperbolic_catalysis(m)
            hyp_a, hyp_i = get!(derived, stripped) do
                if stripped isa EnzymeRates.Mechanism
                    h = _testhelper_den_hyperbolic(stripped)
                    (h, h)
                else
                    _, den_a, _ = EnzymeRates._state_rate_polys(stripped, :A)
                    _, den_i, _ = EnzymeRates._state_rate_polys(stripped, :I)
                    (_testhelper_poly_hyperbolic(den_a, mets),
                     _testhelper_poly_hyperbolic(den_i, mets))
                end
            end
            @test structural == hyp_a
            @test hyp_i || !hyp_a
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
    m = EnzymeRates.Mechanism(@enzyme_mechanism begin
        substrates: S
        products: P
        steps: begin
            E + S ⇌ E(S)
            E(S) <--> E(P)
            E + P ⇌ E(P)
        end
    end)
    @test EnzymeRates._re_segment_count(m) == 1
    groups = EnzymeRates.steps(m)
    g = findfirst(grp -> all(EnzymeRates.is_equilibrium, grp), groups)
    flipped = EnzymeRates.Mechanism(EnzymeRates.reaction(m),
                                    EnzymeRates._flip_group_to_ss(groups, g))
    @test EnzymeRates._re_segment_count(flipped) == 2

    # Allosteric: measured on the A-state projection.
    am = EnzymeRates.AllostericMechanism(EnzymeRates.@allosteric_mechanism begin
        substrates: S
        products: P
        catalytic_multiplicity: 2
        catalytic_steps: begin
            E + S ⇌ E(S)   :: EqualAI
            E(S) <--> E(P) :: EqualAI
            E + P ⇌ E(P)   :: EqualAI
        end
    end)
    @test EnzymeRates._re_segment_count(am) == 1
    am_groups = EnzymeRates.steps(am)
    am_g = findfirst(grp -> all(EnzymeRates.is_equilibrium, grp), am_groups)
    am_flipped = EnzymeRates._with_steps(am,
                                        EnzymeRates._flip_group_to_ss(am_groups, am_g))
    @test EnzymeRates._re_segment_count(am_flipped) == 2
end

@testset "_minimal_gaining_sets" begin
    # Units 1..4. Sets gain iff they contain {1,2} or contain 3.
    gains(set) = (1 in set && 2 in set) || 3 in set
    sets = EnzymeRates._minimal_gaining_sets(gains, _ -> 1:4)
    @test sets == [[3], [1, 2]]
    # Partner pruning: unit 2 may never join unit 1, so {1,2} is unreachable.
    partners(set) = (1 in set || 2 in set) ? [3, 4] : 1:4
    @test EnzymeRates._minimal_gaining_sets(gains, partners) == [[3]]
    # Nothing gains: empty result, and the search terminates.
    @test isempty(EnzymeRates._minimal_gaining_sets(_ -> false, _ -> 1:3))
    # No units at all.
    @test isempty(EnzymeRates._minimal_gaining_sets(_ -> true, _ -> 1:0))
end

"Whole-group flip of groups `gs` (test helper; production uses _flip_group_to_ss)."
function _testhelper_flip_groups(m, gs)
    groups = EnzymeRates.steps(m)
    for g in gs
        groups = EnzymeRates._flip_group_to_ss(groups, g)
    end
    EnzymeRates._with_steps(m, groups)
end

@testset "_re_segment_count_after_flip agrees with the built child" begin
    function check(m)
        groups = EnzymeRates.steps(m)
        eligible = [g for g in eachindex(groups)
                    if all(EnzymeRates.is_equilibrium, groups[g])]
        for g in eligible
            @test EnzymeRates._re_segment_count_after_flip(m, Set([g])) ==
                  EnzymeRates._re_segment_count(_testhelper_flip_groups(m, [g]))
        end
        for (i, g) in enumerate(eligible), h in eligible[(i + 1):end]
            @test EnzymeRates._re_segment_count_after_flip(m, Set([g, h])) ==
                  EnzymeRates._re_segment_count(_testhelper_flip_groups(m, [g, h]))
        end
    end
    # Aggregate pin over the seed set; the inline fixtures below cover the shapes.
    for m in EnzymeRates.init_mechanisms(_testhelper_bibi_rxn)
        check(m)
    end
    am = EnzymeRates.AllostericMechanism(EnzymeRates.@allosteric_mechanism begin
        substrates: S
        products: P
        catalytic_multiplicity: 2
        catalytic_steps: begin
            E + S ⇌ E(S)   :: EqualAI
            E(S) <--> E(P) :: EqualAI
            E + P ⇌ E(P)   :: EqualAI
        end
    end)
    check(am)
end

@testset "_context_bipartitions" begin
    # A binds E, E(B), and E(P): contexts B and P give two bipartitions.
    m = EnzymeRates.Mechanism(@enzyme_mechanism begin
        substrates: A, B
        products: P, Q
        steps: begin
            (E + A ⇌ E(A), E(B) + A ⇌ E(A, B), E(P) + A ⇌ E(A, P))
            (E + B ⇌ E(B), E(A) + B ⇌ E(A, B))
            (E + P ⇌ E(P), E(Q) + P ⇌ E(P, Q), E(A) + P ⇌ E(A, P))
            (E + Q ⇌ E(Q), E(P) + Q ⇌ E(P, Q))
            E(A, B) <--> E(P, Q)
        end
    end)
    src_forms(part) = Set(EnzymeRates.name(EnzymeRates.from_species(s)) for s in part)
    a_group = only(grp for grp in EnzymeRates.steps(m)
                   if length(grp) == 3 &&
                      EnzymeRates.name(EnzymeRates.bound_metabolite(first(grp))) == :A)
    bps = EnzymeRates._context_bipartitions(a_group)
    @test length(bps) == 2
    for (with, without) in bps
        @test !isempty(with) && !isempty(without)
        @test length(with) + length(without) == 3
        @test first(a_group) in with
        @test isempty(intersect(with, without))
    end
    # Contexts are ordered by ligand role (Product before Substrate) then
    # name: bps[1] is the P-context division, bps[2] the B-context division.
    p_step = only(s for s in a_group
                  if any(b -> b isa EnzymeRates.Product,
                         EnzymeRates.bound(EnzymeRates.from_species(s))))
    b_step = only(s for s in a_group
                  if any(b -> b isa EnzymeRates.Substrate,
                         EnzymeRates.bound(EnzymeRates.from_species(s))))
    @test bps[1][2] == [p_step]
    @test bps[2][2] == [b_step]
    # By P: the P-free forms E, E(B) against the P-bound form E(P).
    @test src_forms(bps[1][1]) == Set([:E, :EB])
    @test src_forms(bps[1][2]) == Set([:EP])
    # By B: the B-free forms E, E(P) against the B-bound form E(B).
    @test src_forms(bps[2][1]) == Set([:E, :EP])
    @test src_forms(bps[2][2]) == Set([:EB])
    # _context_form: canonical RE binding puts the metabolite on to_species,
    # so the context form is from_species.
    @test EnzymeRates._context_form(first(a_group)) ==
          EnzymeRates.from_species(first(a_group))
    # _context_form: an SS dissociation step whose bound metabolite is in
    # neither endpoint's bound list (the Segel ping-pong step shape) puts
    # the metabolite on to_species.
    ping_pong_step = EnzymeRates.Step(
        EnzymeRates.Species([EnzymeRates.Substrate(:A)], :E),
        EnzymeRates.Species(EnzymeRates.Metabolite[], :F),
        EnzymeRates.Product(:P), false)
    @test EnzymeRates._context_form(ping_pong_step) ==
          EnzymeRates.to_species(ping_pong_step)
    # A two-step group with one context has one bipartition; a group whose
    # source forms carry no other ligand has none.
    b_group = only(grp for grp in EnzymeRates.steps(m)
                   if length(grp) == 2 &&
                      EnzymeRates.name(EnzymeRates.bound_metabolite(first(grp))) == :B)
    @test length(EnzymeRates._context_bipartitions(b_group)) == 1
    iso_group = only(grp for grp in EnzymeRates.steps(m) if EnzymeRates.is_iso(first(grp)))
    @test isempty(EnzymeRates._context_bipartitions(iso_group))

    # Ter-ter substrate cube: each of the two other substrates divides the A
    # group's four steps two and two.
    cube = EnzymeRates.Mechanism(@enzyme_mechanism begin
        substrates: A, B, C
        products: P, Q, R
        steps: begin
            (E + A ⇌ E(A), E(B) + A ⇌ E(A, B), E(C) + A ⇌ E(A, C),
             E(B, C) + A ⇌ E(A, B, C))
            (E + B ⇌ E(B), E(A) + B ⇌ E(A, B), E(C) + B ⇌ E(B, C),
             E(A, C) + B ⇌ E(A, B, C))
            (E + C ⇌ E(C), E(A) + C ⇌ E(A, C), E(B) + C ⇌ E(B, C),
             E(A, B) + C ⇌ E(A, B, C))
            E(A, B, C) <--> E(P, Q, R)
            E(Q, R) + P ⇌ E(P, Q, R)
            E(R) + Q ⇌ E(Q, R)
            E + R ⇌ E(R)
        end
    end)
    cube_a = only(grp for grp in EnzymeRates.steps(cube)
                  if length(grp) == 4 &&
                     EnzymeRates.name(EnzymeRates.bound_metabolite(first(grp))) == :A)
    cube_bps = EnzymeRates._context_bipartitions(cube_a)
    @test length(cube_bps) == 2
    # Both contexts are substrates, so they are ordered by name: B then C.
    # By B: the B-free forms E, E(C) against the B-bound forms E(B), E(B,C).
    @test src_forms(cube_bps[1][1]) == Set([:E, :EC])
    @test src_forms(cube_bps[1][2]) == Set([:EB, :EBC])
    # By C: the C-free forms E, E(B) against the C-bound forms E(C), E(B,C).
    @test src_forms(cube_bps[2][1]) == Set([:E, :EB])
    @test src_forms(cube_bps[2][2]) == Set([:EC, :EBC])

    # Residual as a context: Q binds free E, the modified enzyme, and the
    # modified enzyme with B bound. Context B divides {E, E(;res)} from
    # {E(B;res)}; the residual divides {E} from the two modified forms.
    res = EnzymeRates.Mechanism(@enzyme_mechanism begin
        substrates: A, B
        products: P, Q
        steps: begin
            (E + Q ⇌ E(Q), E(; residual = A - P) + Q ⇌ E(Q; residual = A - P),
             E(B; residual = A - P) + Q ⇌ E(B, Q; residual = A - P))
            E + A ⇌ E(A)
            E(A) <--> E(P; residual = A - P)
            E(P; residual = A - P) ⇌ E(; residual = A - P) + P
            E(; residual = A - P) + B ⇌ E(B; residual = A - P)
            E(B; residual = A - P) <--> E(Q)
        end
    end)
    q_group = only(grp for grp in EnzymeRates.steps(res) if length(grp) == 3)
    rbps = EnzymeRates._context_bipartitions(q_group)
    @test length(rbps) == 2
    Er, EBr = Symbol("E_res_+A_-P"), Symbol("EB_res_+A_-P")
    # By B, then by the residual.
    @test src_forms(rbps[1][1]) == Set([:E, Er]) && src_forms(rbps[1][2]) == Set([EBr])
    @test src_forms(rbps[2][1]) == Set([:E]) && src_forms(rbps[2][2]) == Set([Er, EBr])

    # Conformation as a context: A binds E, Estar, and Estar(B).
    conf = EnzymeRates.Mechanism(@enzyme_mechanism begin
        substrates: A, B
        products: P, Q
        steps: begin
            (E + A ⇌ E(A), Estar + A ⇌ Estar(A), Estar(B) + A ⇌ Estar(A, B))
            E(A) <--> Estar(A)
            Estar(A, B) <--> Estar(P, Q)
            Estar(P, Q) ⇌ Estar(P) + Q
            Estar(P) ⇌ E + P
        end
    end)
    a_conf_group = only(grp for grp in EnzymeRates.steps(conf) if length(grp) == 3)
    cbps = EnzymeRates._context_bipartitions(a_conf_group)
    @test length(cbps) == 2
    # By B, then by the conformation.
    @test src_forms(cbps[1][1]) == Set([:E, :Estar])
    @test src_forms(cbps[1][2]) == Set([:EstarB])
    @test src_forms(cbps[2][1]) == Set([:E])
    @test src_forms(cbps[2][2]) == Set([:Estar, :EstarB])
end

@testset "_context_bipartitions separates an inhibitor-bound mirror" begin
    am = EnzymeRates.AllostericMechanism(EnzymeRates.@allosteric_mechanism begin
        substrates: S
        products: P
        catalytic_inhibitors: I
        catalytic_steps: begin
            (E + S ⇌ E(S), E(I) + S ⇌ E(I, S))   :: EqualAI
            E(S) <--> E(P)                        :: EqualAI
            E + P ⇌ E(P)                          :: EqualAI
            E + I ⇌ E(I)                          :: EqualAI
        end
    end)
    s_group = only(grp for grp in EnzymeRates.steps(am) if length(grp) == 2)
    bps = EnzymeRates._context_bipartitions(s_group)
    @test length(bps) == 1
    @test all(part -> length(part) == 1, bps[1])
end

@testset "_apply_bipartitions" begin
    m = EnzymeRates.Mechanism(@enzyme_mechanism begin
        substrates: A, B
        products: P, Q
        steps: begin
            (E + A ⇌ E(A), E(B) + A ⇌ E(A, B))
            (E + B ⇌ E(B), E(A) + B ⇌ E(A, B))
            (E + P ⇌ E(P), E(Q) + P ⇌ E(P, Q))
            (E + Q ⇌ E(Q), E(P) + Q ⇌ E(P, Q))
            E(A, B) <--> E(P, Q)
        end
    end)
    groups = EnzymeRates.steps(m)
    g = findfirst(grp -> length(grp) == 2, groups)
    bp = only(EnzymeRates._context_bipartitions(groups[g]))
    child = EnzymeRates._apply_bipartitions(m, [(g, bp)])
    @test length(EnzymeRates.steps(child)) == length(groups) + 1
    @test EnzymeRates.n_steps(child) == EnzymeRates.n_steps(m)
    @test Set(s for grp in EnzymeRates.steps(child) for s in grp) ==
          Set(s for grp in groups for s in grp)
    @test any(grp -> Set(grp) == Set(bp[1]), EnzymeRates.steps(child))
    @test any(grp -> Set(grp) == Set(bp[2]), EnzymeRates.steps(child))

    am = EnzymeRates.AllostericMechanism(@allosteric_mechanism begin
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
    end)
    ga = findfirst(grp -> length(grp) == 2 &&
                   EnzymeRates.name(EnzymeRates.bound_metabolite(first(grp))) == :A,
                   EnzymeRates.steps(am))
    bpa = only(EnzymeRates._context_bipartitions(EnzymeRates.steps(am)[ga]))
    achild = EnzymeRates._apply_bipartitions(am, [(ga, bpa)])
    @test length(EnzymeRates.cat_allo_states(achild)) == length(EnzymeRates.steps(achild))
    for (gi, grp) in enumerate(EnzymeRates.steps(achild))
        Set(grp) ⊆ Set(bpa[1]) || Set(grp) ⊆ Set(bpa[2]) || continue
        @test EnzymeRates.cat_allo_state(achild, gi) == :NonequalAI
    end
    @test EnzymeRates.catalytic_multiplicity(achild) == 2
    @test EnzymeRates.regulatory_sites(achild) == EnzymeRates.regulatory_sites(am)
end

@testset "_partition_independent_count agrees with _independent_param_count" begin
    seeds = EnzymeRates.init_mechanisms(_testhelper_bibi_rxn)
    checked = 0
    for m in seeds
        counter = EnzymeRates._partition_independent_count(m)
        flat = EnzymeRates._flat_steps(m)
        parent_ids = [g for (_, g) in flat]
        @test counter(parent_ids) == EnzymeRates._independent_param_count(m)
        pos = Dict(s => j for (j, (s, _)) in enumerate(flat))
        groups = EnzymeRates.steps(m)
        for g in eachindex(groups), bp in EnzymeRates._context_bipartitions(groups[g])
            ids = copy(parent_ids)
            for s in bp[2]
                ids[pos[s]] = length(groups) + 1
            end
            child = EnzymeRates._apply_bipartitions(m, [(g, bp)])
            @test counter(ids) == EnzymeRates._independent_param_count(child)
            checked += 1
        end
    end
    @test checked > 100
end

"Finite-difference rank of ∂v/∂log θ over the fitted parameters (test oracle only)."
function _testhelper_identifiable_rank(m; npts = 60, ndraws = 3, h = 1e-5)
    em = EnzymeRates.compile_mechanism(m)
    fp = collect(EnzymeRates.fitted_params(em))
    cm = m isa EnzymeRates.Mechanism ? m : EnzymeRates._state_mechanism(m, :A)
    mets = sort!(collect(EnzymeRates._concentration_symbols(cm)))
    # Seeded from the rendered equation: a struct hash mixes in objectid and
    # would not reproduce a failure in a later session.
    rng = MersenneTwister(hash(rate_equation_string(em)) % 2^31)
    best = 0
    for _ in 1:ndraws
        θ = exp.(randn(rng, length(fp)))
        keq = exp(randn(rng))
        concs = [NamedTuple{Tuple(mets)}(Tuple(exp.(2 .* randn(rng, length(mets)))))
                 for _ in 1:npts]
        J = zeros(npts, length(fp))
        for j in eachindex(fp), sgn in (1, -1)
            θp = copy(θ); θp[j] *= exp(sgn * h)
            p = NamedTuple{(fp..., :Keq, :E_total)}((θp..., keq, 1.0))
            for (i, c) in enumerate(concs)
                J[i, j] += sgn * rate_equation(em, c, p) / (2h)
            end
        end
        all(isfinite, J) || continue
        sv = svdvals(J)
        best = max(best, count(>(1e-7 * sv[1]), sv))
    end
    best
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
    @testset "_expand_re_to_ss: seed child count is unchanged (220 over bi-bi seeds)" begin
        seeds = EnzymeRates.init_mechanisms(_testhelper_bibi_rxn)
        @test length(seeds) == 55
        @test sum(length(EnzymeRates._expand_re_to_ss(m)) for m in seeds) == 220
    end

    @testset "_expand_re_to_ss: uni-uni emits both single-group flips exactly" begin
        # Uni-uni steady-state and rapid-equilibrium laws have the same form, so
        # these two children are documented no-ops; the move still emits them
        # because their rejection has no structural proof.
        m = EnzymeRates.Mechanism(@enzyme_mechanism begin
            substrates: S
            products: P
            steps: begin
                E + S ⇌ E(S)
                E(S) <--> E(P)
                E + P ⇌ E(P)
            end
        end)
        expected = [
            EnzymeRates.Mechanism(@enzyme_mechanism begin
                substrates: S
                products: P
                steps: begin
                    E + S <--> E(S)
                    E(S) <--> E(P)
                    E + P ⇌ E(P)
                end
            end),
            EnzymeRates.Mechanism(@enzyme_mechanism begin
                substrates: S
                products: P
                steps: begin
                    E + S ⇌ E(S)
                    E(S) <--> E(P)
                    E + P <--> E(P)
                end
            end),
        ]
        kids = EnzymeRates._expand_re_to_ss(m)
        @test length(kids) == 2
        @test Set(kids) == Set(expected)
    end

    @testset "_expand_re_to_ss: random-order bi-bi emits one flip per metabolite" begin
        # Every metabolite's two binding steps share one group, so its single group cuts a
        # segment alone, raising the segment count on its own; each metabolite is emitted
        # alone and no pair is minimal.
        m = EnzymeRates.Mechanism(@enzyme_mechanism begin
            substrates: A, B
            products: P, Q
            steps: begin
                (E + A ⇌ E(A), E(B) + A ⇌ E(A, B))
                (E + B ⇌ E(B), E(A) + B ⇌ E(A, B))
                (E + P ⇌ E(P), E(Q) + P ⇌ E(P, Q))
                (E + Q ⇌ E(Q), E(P) + Q ⇌ E(P, Q))
                E(A, B) <--> E(P, Q)
            end
        end)
        flipA = EnzymeRates.Mechanism(@enzyme_mechanism begin
            substrates: A, B
            products: P, Q
            steps: begin
                (E + A <--> E(A), E(B) + A <--> E(A, B))
                (E + B ⇌ E(B), E(A) + B ⇌ E(A, B))
                (E + P ⇌ E(P), E(Q) + P ⇌ E(P, Q))
                (E + Q ⇌ E(Q), E(P) + Q ⇌ E(P, Q))
                E(A, B) <--> E(P, Q)
            end
        end)
        flipB = EnzymeRates.Mechanism(@enzyme_mechanism begin
            substrates: A, B
            products: P, Q
            steps: begin
                (E + A ⇌ E(A), E(B) + A ⇌ E(A, B))
                (E + B <--> E(B), E(A) + B <--> E(A, B))
                (E + P ⇌ E(P), E(Q) + P ⇌ E(P, Q))
                (E + Q ⇌ E(Q), E(P) + Q ⇌ E(P, Q))
                E(A, B) <--> E(P, Q)
            end
        end)
        flipP = EnzymeRates.Mechanism(@enzyme_mechanism begin
            substrates: A, B
            products: P, Q
            steps: begin
                (E + A ⇌ E(A), E(B) + A ⇌ E(A, B))
                (E + B ⇌ E(B), E(A) + B ⇌ E(A, B))
                (E + P <--> E(P), E(Q) + P <--> E(P, Q))
                (E + Q ⇌ E(Q), E(P) + Q ⇌ E(P, Q))
                E(A, B) <--> E(P, Q)
            end
        end)
        flipQ = EnzymeRates.Mechanism(@enzyme_mechanism begin
            substrates: A, B
            products: P, Q
            steps: begin
                (E + A ⇌ E(A), E(B) + A ⇌ E(A, B))
                (E + B ⇌ E(B), E(A) + B ⇌ E(A, B))
                (E + P ⇌ E(P), E(Q) + P ⇌ E(P, Q))
                (E + Q <--> E(Q), E(P) + Q <--> E(P, Q))
                E(A, B) <--> E(P, Q)
            end
        end)
        kids = EnzymeRates._expand_re_to_ss(m)
        @test length(kids) == 4
        @test Set(kids) == Set([flipA, flipB, flipP, flipQ])
    end

    @testset "_expand_re_to_ss: a dead-end leaf group never flips" begin
        # E(P) + S ⇌ E(P, S) is a bridge: no cycle through it contains chemistry,
        # so at steady state it carries no net flux and flipping it would add a
        # parameter the data cannot see. Only the two catalytic bindings flip.
        m = EnzymeRates.Mechanism(@enzyme_mechanism begin
            substrates: S
            products: P
            steps: begin
                E + S ⇌ E(S)
                E(S) <--> E(P)
                E + P ⇌ E(P)
                E(P) + S ⇌ E(P, S)
            end
        end)
        expected = [
            EnzymeRates.Mechanism(@enzyme_mechanism begin
                substrates: S
                products: P
                steps: begin
                    E + S <--> E(S)
                    E(S) <--> E(P)
                    E + P ⇌ E(P)
                    E(P) + S ⇌ E(P, S)
                end
            end),
            EnzymeRates.Mechanism(@enzyme_mechanism begin
                substrates: S
                products: P
                steps: begin
                    E + S ⇌ E(S)
                    E(S) <--> E(P)
                    E + P <--> E(P)
                    E(P) + S ⇌ E(P, S)
                end
            end),
        ]
        kids = EnzymeRates._expand_re_to_ss(m)
        @test length(kids) == 2
        @test Set(kids) == Set(expected)
    end

    @testset "_expand_re_to_ss: a split metabolite flips as a pair, never singly" begin
        # A's two binding steps sit in separate groups. Flipping either alone
        # leaves E and E(A) joined through the other A step (E–E(B)–E(A,B)–E(A)),
        # so the segment count does not rise and no such child exists. The pair
        # cuts the square, so it is emitted; B, P, and Q each flip alone, since each
        # cuts a segment on its own.
        m = EnzymeRates.Mechanism(@enzyme_mechanism begin
            substrates: A, B
            products: P, Q
            steps: begin
                E + A ⇌ E(A)
                E(B) + A ⇌ E(A, B)
                (E + B ⇌ E(B), E(A) + B ⇌ E(A, B))
                (E + P ⇌ E(P), E(Q) + P ⇌ E(P, Q))
                (E + Q ⇌ E(Q), E(P) + Q ⇌ E(P, Q))
                E(A, B) <--> E(P, Q)
            end
        end)
        flipA = EnzymeRates.Mechanism(@enzyme_mechanism begin
            substrates: A, B
            products: P, Q
            steps: begin
                E + A <--> E(A)
                E(B) + A <--> E(A, B)
                (E + B ⇌ E(B), E(A) + B ⇌ E(A, B))
                (E + P ⇌ E(P), E(Q) + P ⇌ E(P, Q))
                (E + Q ⇌ E(Q), E(P) + Q ⇌ E(P, Q))
                E(A, B) <--> E(P, Q)
            end
        end)
        flipB = EnzymeRates.Mechanism(@enzyme_mechanism begin
            substrates: A, B
            products: P, Q
            steps: begin
                E + A ⇌ E(A)
                E(B) + A ⇌ E(A, B)
                (E + B <--> E(B), E(A) + B <--> E(A, B))
                (E + P ⇌ E(P), E(Q) + P ⇌ E(P, Q))
                (E + Q ⇌ E(Q), E(P) + Q ⇌ E(P, Q))
                E(A, B) <--> E(P, Q)
            end
        end)
        flipP = EnzymeRates.Mechanism(@enzyme_mechanism begin
            substrates: A, B
            products: P, Q
            steps: begin
                E + A ⇌ E(A)
                E(B) + A ⇌ E(A, B)
                (E + B ⇌ E(B), E(A) + B ⇌ E(A, B))
                (E + P <--> E(P), E(Q) + P <--> E(P, Q))
                (E + Q ⇌ E(Q), E(P) + Q ⇌ E(P, Q))
                E(A, B) <--> E(P, Q)
            end
        end)
        flipQ = EnzymeRates.Mechanism(@enzyme_mechanism begin
            substrates: A, B
            products: P, Q
            steps: begin
                E + A ⇌ E(A)
                E(B) + A ⇌ E(A, B)
                (E + B ⇌ E(B), E(A) + B ⇌ E(A, B))
                (E + P ⇌ E(P), E(Q) + P ⇌ E(P, Q))
                (E + Q <--> E(Q), E(P) + Q <--> E(P, Q))
                E(A, B) <--> E(P, Q)
            end
        end)
        kids = EnzymeRates._expand_re_to_ss(m)
        @test length(kids) == 4
        @test Set(kids) == Set([flipA, flipB, flipP, flipQ])
        # The two single-A flips are absent: each is segment-flat (E and E(A) stay joined
        # through the other A step), so the move does not emit it.
        for single in (
            EnzymeRates.Mechanism(@enzyme_mechanism begin
                substrates: A, B
                products: P, Q
                steps: begin
                    E + A <--> E(A)
                    E(B) + A ⇌ E(A, B)
                    (E + B ⇌ E(B), E(A) + B ⇌ E(A, B))
                    (E + P ⇌ E(P), E(Q) + P ⇌ E(P, Q))
                    (E + Q ⇌ E(Q), E(P) + Q ⇌ E(P, Q))
                    E(A, B) <--> E(P, Q)
                end
            end),
            EnzymeRates.Mechanism(@enzyme_mechanism begin
                substrates: A, B
                products: P, Q
                steps: begin
                    E + A ⇌ E(A)
                    E(B) + A <--> E(A, B)
                    (E + B ⇌ E(B), E(A) + B ⇌ E(A, B))
                    (E + P ⇌ E(P), E(Q) + P ⇌ E(P, Q))
                    (E + Q ⇌ E(Q), E(P) + Q ⇌ E(P, Q))
                    E(A, B) <--> E(P, Q)
                end
            end))
            @test !(single in kids)
            @test EnzymeRates._re_segment_count(single) == EnzymeRates._re_segment_count(m)
        end
    end

    @testset "_expand_re_to_ss: ter-ter substrate cube emits one flip per group" begin
        # Random substrate addition, ordered product release. Each of the three
        # substrate groups carries all four cube edges for its metabolite, so
        # flipping one cuts every route between the forms it joins; each product
        # release is a single step on the ordered path. Six groups, six children.
        m = EnzymeRates.Mechanism(@enzyme_mechanism begin
            substrates: A, B, C
            products: P, Q, R
            steps: begin
                (E + A ⇌ E(A), E(B) + A ⇌ E(A, B), E(C) + A ⇌ E(A, C),
                 E(B, C) + A ⇌ E(A, B, C))
                (E + B ⇌ E(B), E(A) + B ⇌ E(A, B), E(C) + B ⇌ E(B, C),
                 E(A, C) + B ⇌ E(A, B, C))
                (E + C ⇌ E(C), E(A) + C ⇌ E(A, C), E(B) + C ⇌ E(B, C),
                 E(A, B) + C ⇌ E(A, B, C))
                E(A, B, C) <--> E(P, Q, R)
                E(Q, R) + P ⇌ E(P, Q, R)
                E(R) + Q ⇌ E(Q, R)
                E + R ⇌ E(R)
            end
        end)
        flipA = EnzymeRates.Mechanism(@enzyme_mechanism begin
            substrates: A, B, C
            products: P, Q, R
            steps: begin
                (E + A <--> E(A), E(B) + A <--> E(A, B),
                 E(C) + A <--> E(A, C), E(B, C) + A <--> E(A, B, C))
                (E + B ⇌ E(B), E(A) + B ⇌ E(A, B), E(C) + B ⇌ E(B, C),
                 E(A, C) + B ⇌ E(A, B, C))
                (E + C ⇌ E(C), E(A) + C ⇌ E(A, C), E(B) + C ⇌ E(B, C),
                 E(A, B) + C ⇌ E(A, B, C))
                E(A, B, C) <--> E(P, Q, R)
                E(Q, R) + P ⇌ E(P, Q, R)
                E(R) + Q ⇌ E(Q, R)
                E + R ⇌ E(R)
            end
        end)
        flipB = EnzymeRates.Mechanism(@enzyme_mechanism begin
            substrates: A, B, C
            products: P, Q, R
            steps: begin
                (E + A ⇌ E(A), E(B) + A ⇌ E(A, B), E(C) + A ⇌ E(A, C),
                 E(B, C) + A ⇌ E(A, B, C))
                (E + B <--> E(B), E(A) + B <--> E(A, B),
                 E(C) + B <--> E(B, C), E(A, C) + B <--> E(A, B, C))
                (E + C ⇌ E(C), E(A) + C ⇌ E(A, C), E(B) + C ⇌ E(B, C),
                 E(A, B) + C ⇌ E(A, B, C))
                E(A, B, C) <--> E(P, Q, R)
                E(Q, R) + P ⇌ E(P, Q, R)
                E(R) + Q ⇌ E(Q, R)
                E + R ⇌ E(R)
            end
        end)
        flipC = EnzymeRates.Mechanism(@enzyme_mechanism begin
            substrates: A, B, C
            products: P, Q, R
            steps: begin
                (E + A ⇌ E(A), E(B) + A ⇌ E(A, B), E(C) + A ⇌ E(A, C),
                 E(B, C) + A ⇌ E(A, B, C))
                (E + B ⇌ E(B), E(A) + B ⇌ E(A, B), E(C) + B ⇌ E(B, C),
                 E(A, C) + B ⇌ E(A, B, C))
                (E + C <--> E(C), E(A) + C <--> E(A, C),
                 E(B) + C <--> E(B, C), E(A, B) + C <--> E(A, B, C))
                E(A, B, C) <--> E(P, Q, R)
                E(Q, R) + P ⇌ E(P, Q, R)
                E(R) + Q ⇌ E(Q, R)
                E + R ⇌ E(R)
            end
        end)
        flipP = EnzymeRates.Mechanism(@enzyme_mechanism begin
            substrates: A, B, C
            products: P, Q, R
            steps: begin
                (E + A ⇌ E(A), E(B) + A ⇌ E(A, B), E(C) + A ⇌ E(A, C),
                 E(B, C) + A ⇌ E(A, B, C))
                (E + B ⇌ E(B), E(A) + B ⇌ E(A, B), E(C) + B ⇌ E(B, C),
                 E(A, C) + B ⇌ E(A, B, C))
                (E + C ⇌ E(C), E(A) + C ⇌ E(A, C), E(B) + C ⇌ E(B, C),
                 E(A, B) + C ⇌ E(A, B, C))
                E(A, B, C) <--> E(P, Q, R)
                E(Q, R) + P <--> E(P, Q, R)
                E(R) + Q ⇌ E(Q, R)
                E + R ⇌ E(R)
            end
        end)
        flipQ = EnzymeRates.Mechanism(@enzyme_mechanism begin
            substrates: A, B, C
            products: P, Q, R
            steps: begin
                (E + A ⇌ E(A), E(B) + A ⇌ E(A, B), E(C) + A ⇌ E(A, C),
                 E(B, C) + A ⇌ E(A, B, C))
                (E + B ⇌ E(B), E(A) + B ⇌ E(A, B), E(C) + B ⇌ E(B, C),
                 E(A, C) + B ⇌ E(A, B, C))
                (E + C ⇌ E(C), E(A) + C ⇌ E(A, C), E(B) + C ⇌ E(B, C),
                 E(A, B) + C ⇌ E(A, B, C))
                E(A, B, C) <--> E(P, Q, R)
                E(Q, R) + P ⇌ E(P, Q, R)
                E(R) + Q <--> E(Q, R)
                E + R ⇌ E(R)
            end
        end)
        flipR = EnzymeRates.Mechanism(@enzyme_mechanism begin
            substrates: A, B, C
            products: P, Q, R
            steps: begin
                (E + A ⇌ E(A), E(B) + A ⇌ E(A, B), E(C) + A ⇌ E(A, C),
                 E(B, C) + A ⇌ E(A, B, C))
                (E + B ⇌ E(B), E(A) + B ⇌ E(A, B), E(C) + B ⇌ E(B, C),
                 E(A, C) + B ⇌ E(A, B, C))
                (E + C ⇌ E(C), E(A) + C ⇌ E(A, C), E(B) + C ⇌ E(B, C),
                 E(A, B) + C ⇌ E(A, B, C))
                E(A, B, C) <--> E(P, Q, R)
                E(Q, R) + P ⇌ E(P, Q, R)
                E(R) + Q ⇌ E(Q, R)
                E + R <--> E(R)
            end
        end)
        kids = EnzymeRates._expand_re_to_ss(m)
        @test length(kids) == 6
        @test Set(kids) == Set([flipA, flipB, flipC, flipP, flipQ, flipR])
    end

    @testset "_expand_re_to_ss: split ter-ter pairs the A and B parts" begin
        # Ter-ter after one context split: A is split by whether B is bound and B
        # by whether A is bound. Each of the four parts alone leaves E and E(A)
        # joined through its sibling, so nothing flips one part by itself; five
        # of the six pairs that cut a segment are emitted, and C, P, Q, R still
        # flip alone. The sixth pair, {A1, B1}, cuts both binding edges at E and
        # E(C) and leaves {E(A), E(B), E(A, B), E(A, C), E(B, C), E(A, B, C)} as a
        # segment with no substrate-free form, which the constructor rejects;
        # its supersets are reached from the {A1, A2} and {B1, B2} children.
        m = EnzymeRates.Mechanism(@enzyme_mechanism begin
            substrates: A, B, C
            products: P, Q, R
            steps: begin
                (E + A ⇌ E(A), E(C) + A ⇌ E(A, C))
                (E(B) + A ⇌ E(A, B), E(B, C) + A ⇌ E(A, B, C))
                (E + B ⇌ E(B), E(C) + B ⇌ E(B, C))
                (E(A) + B ⇌ E(A, B), E(A, C) + B ⇌ E(A, B, C))
                (E + C ⇌ E(C), E(A) + C ⇌ E(A, C), E(B) + C ⇌ E(B, C),
                 E(A, B) + C ⇌ E(A, B, C))
                E(A, B, C) <--> E(P, Q, R)
                E(Q, R) + P ⇌ E(P, Q, R)
                E(R) + Q ⇌ E(Q, R)
                E + R ⇌ E(R)
            end
        end)
        binds(grp) = (bm = EnzymeRates.bound_metabolite(first(grp));
                      bm === nothing ? :iso : EnzymeRates.name(bm))
        sources(grp) = Set(EnzymeRates.name(EnzymeRates.from_species(s)) for s in grp)
        # Each group is pinned by the metabolite it binds and the forms it binds to.
        A1, A2 = (:A, Set([:E, :EC])), (:A, Set([:EB, :EBC]))
        B1, B2 = (:B, Set([:E, :EC])), (:B, Set([:EA, :EAC]))
        C = (:C, Set([:E, :EA, :EB, :EAB]))
        P, Q, R = (:P, Set([:EQR])), (:Q, Set([:ER])), (:R, Set([:E]))
        group_of(key) = only(g for (g, grp) in enumerate(EnzymeRates.steps(m))
                             if (binds(grp), sources(grp)) == key)
        flip(keys...) = _testhelper_flip_groups(m, [group_of(k) for k in keys])
        expected = [
            flip(C),        # {C}
            flip(P),        # {P}
            flip(Q),        # {Q}
            flip(R),        # {R}
            flip(A1, A2),   # {A1, A2}
            flip(B1, B2),   # {B1, B2}
            flip(A1, B2),   # {A1, B2}
            flip(A2, B1),   # {A2, B1}
            flip(A2, B2),   # {A2, B2}
        ]
        kids = EnzymeRates._expand_re_to_ss(m)
        @test length(kids) == 9
        @test Set(kids) == Set(expected)
        for c in kids
            ss = count((A1, A2, B1, B2)) do key
                grp = only(grp for grp in EnzymeRates.steps(c)
                           if (binds(grp), sources(grp)) == key)
                !any(EnzymeRates.is_equilibrium, grp)
            end
            @test ss != 1
        end
    end

    @testset "_expand_re_to_ss: ping-pong" begin
        # The rapid-equilibrium subgraph is two trees, {E, E(A), E(Q)} and
        # {E(P; res), E(; res), E(B; res)}, so every RE group is a bridge and
        # flipping any one of them alone raises the segment count; no pair is
        # minimal.
        m = EnzymeRates.Mechanism(@enzyme_mechanism begin
            substrates: A, B
            products: P, Q
            steps: begin
                E + A ⇌ E(A)
                E(A) <--> E(P; residual = A - P)
                E(P; residual = A - P) ⇌ E(; residual = A - P) + P
                E(; residual = A - P) + B ⇌ E(B; residual = A - P)
                E(B; residual = A - P) <--> E(Q)
                E(Q) ⇌ E + Q
            end
        end)
        flipA = EnzymeRates.Mechanism(@enzyme_mechanism begin
            substrates: A, B
            products: P, Q
            steps: begin
                E + A <--> E(A)
                E(A) <--> E(P; residual = A - P)
                E(P; residual = A - P) ⇌ E(; residual = A - P) + P
                E(; residual = A - P) + B ⇌ E(B; residual = A - P)
                E(B; residual = A - P) <--> E(Q)
                E(Q) ⇌ E + Q
            end
        end)
        flipP = EnzymeRates.Mechanism(@enzyme_mechanism begin
            substrates: A, B
            products: P, Q
            steps: begin
                E + A ⇌ E(A)
                E(A) <--> E(P; residual = A - P)
                E(P; residual = A - P) <--> E(; residual = A - P) + P
                E(; residual = A - P) + B ⇌ E(B; residual = A - P)
                E(B; residual = A - P) <--> E(Q)
                E(Q) ⇌ E + Q
            end
        end)
        flipB = EnzymeRates.Mechanism(@enzyme_mechanism begin
            substrates: A, B
            products: P, Q
            steps: begin
                E + A ⇌ E(A)
                E(A) <--> E(P; residual = A - P)
                E(P; residual = A - P) ⇌ E(; residual = A - P) + P
                E(; residual = A - P) + B <--> E(B; residual = A - P)
                E(B; residual = A - P) <--> E(Q)
                E(Q) ⇌ E + Q
            end
        end)
        flipQ = EnzymeRates.Mechanism(@enzyme_mechanism begin
            substrates: A, B
            products: P, Q
            steps: begin
                E + A ⇌ E(A)
                E(A) <--> E(P; residual = A - P)
                E(P; residual = A - P) ⇌ E(; residual = A - P) + P
                E(; residual = A - P) + B ⇌ E(B; residual = A - P)
                E(B; residual = A - P) <--> E(Q)
                E(Q) <--> E + Q
            end
        end)
        kids = EnzymeRates._expand_re_to_ss(m)
        @test length(kids) == 4
        @test Set(kids) == Set([flipA, flipP, flipB, flipQ])
    end

    @testset "_expand_re_to_ss: an allosteric parent keeps a hyperbolic scheme" begin
        # Uni-bi with random product release. Flipping the S group alone splits
        # off {E(S)} and the equation stays degree 1. Flipping the P group splits
        # off {E(P), E(P, Q)}, whose term carries Q from E(P, Q) beyond its
        # bottom E(P) and again from the edge E(Q) → E(P, Q), so the equation
        # carries Q²; flipping Q mirrors it with P². A plain Mechanism emits all
        # three; an allosteric parent emits only the S flip.
        plain = EnzymeRates.Mechanism(@enzyme_mechanism begin
            substrates: S
            products: P, Q
            steps: begin
                E + S ⇌ E(S)
                E(S) <--> E(P, Q)
                (E + P ⇌ E(P), E(Q) + P ⇌ E(P, Q))
                (E + Q ⇌ E(Q), E(P) + Q ⇌ E(P, Q))
            end
        end)
        flipS = EnzymeRates.Mechanism(@enzyme_mechanism begin
            substrates: S
            products: P, Q
            steps: begin
                E + S <--> E(S)
                E(S) <--> E(P, Q)
                (E + P ⇌ E(P), E(Q) + P ⇌ E(P, Q))
                (E + Q ⇌ E(Q), E(P) + Q ⇌ E(P, Q))
            end
        end)
        flipP = EnzymeRates.Mechanism(@enzyme_mechanism begin
            substrates: S
            products: P, Q
            steps: begin
                E + S ⇌ E(S)
                E(S) <--> E(P, Q)
                (E + P <--> E(P), E(Q) + P <--> E(P, Q))
                (E + Q ⇌ E(Q), E(P) + Q ⇌ E(P, Q))
            end
        end)
        flipQ = EnzymeRates.Mechanism(@enzyme_mechanism begin
            substrates: S
            products: P, Q
            steps: begin
                E + S ⇌ E(S)
                E(S) <--> E(P, Q)
                (E + P ⇌ E(P), E(Q) + P ⇌ E(P, Q))
                (E + Q <--> E(Q), E(P) + Q <--> E(P, Q))
            end
        end)
        kids = EnzymeRates._expand_re_to_ss(plain)
        @test length(kids) == 3
        @test Set(kids) == Set([flipS, flipP, flipQ])

        allo = EnzymeRates.AllostericMechanism(@allosteric_mechanism begin
            substrates: S
            products: P, Q
            catalytic_multiplicity: 2
            catalytic_steps: begin
                E + S ⇌ E(S)                                :: NonequalAI
                E(S) <--> E(P, Q)                           :: NonequalAI
                (E + P ⇌ E(P), E(Q) + P ⇌ E(P, Q))          :: NonequalAI
                (E + Q ⇌ E(Q), E(P) + Q ⇌ E(P, Q))          :: NonequalAI
            end
        end)
        allo_flipS = EnzymeRates.AllostericMechanism(@allosteric_mechanism begin
            substrates: S
            products: P, Q
            catalytic_multiplicity: 2
            catalytic_steps: begin
                E + S <--> E(S)                             :: NonequalAI
                E(S) <--> E(P, Q)                           :: NonequalAI
                (E + P ⇌ E(P), E(Q) + P ⇌ E(P, Q))          :: NonequalAI
                (E + Q ⇌ E(Q), E(P) + Q ⇌ E(P, Q))          :: NonequalAI
            end
        end)
        allo_kids = EnzymeRates._expand_re_to_ss(allo)
        @test length(allo_kids) == 1
        @test Set(allo_kids) == Set([allo_flipS])
    end
end

@testset "_expand_split_kinetic_group (context bipartitions)" begin
    @testset "_expand_split_kinetic_group: random-order bi-bi frees indep params" begin
        # Today's canonicalization drops every split of this seed. The split closure
        # must reach the form with every step in its own group and 9 independent
        # parameters (measured: 16 structures).
        m = EnzymeRates.Mechanism(@enzyme_mechanism begin
            substrates: A, B
            products: P, Q
            steps: begin
                (E + A ⇌ E(A), E(B) + A ⇌ E(A, B), E(P) + A ⇌ E(A, P))
                (E + B ⇌ E(B), E(A) + B ⇌ E(A, B), E(Q) + B ⇌ E(B, Q))
                (E + P ⇌ E(P), E(Q) + P ⇌ E(P, Q), E(A) + P ⇌ E(A, P))
                (E + Q ⇌ E(Q), E(P) + Q ⇌ E(P, Q), E(B) + Q ⇌ E(B, Q))
                E(A, B) <--> E(P, Q)
            end
        end)
        @test EnzymeRates._independent_param_count(m) == 5
        cl = _testhelper_closure(m, EnzymeRates._expand_split_kinetic_group)
        @test maximum(length(EnzymeRates.steps(m)) for m in cl) == 13
        @test maximum(EnzymeRates._independent_param_count(m) for m in cl) == 9
        @test length(cl) == 16
    end

    @testset "_expand_split_kinetic_group: bi-bi without dead ends splits a square" begin
        # A single context split (A by B alone) is tied straight back by the A/B
        # square, so no child splits one group by itself. Each child frees one
        # interaction factor by separating both groups of one square.
        m = EnzymeRates.Mechanism(@enzyme_mechanism begin
            substrates: A, B
            products: P, Q
            steps: begin
                (E + A ⇌ E(A), E(B) + A ⇌ E(A, B))
                (E + B ⇌ E(B), E(A) + B ⇌ E(A, B))
                (E + P ⇌ E(P), E(Q) + P ⇌ E(P, Q))
                (E + Q ⇌ E(Q), E(P) + Q ⇌ E(P, Q))
                E(A, B) <--> E(P, Q)
            end
        end)
        splitAB = EnzymeRates.Mechanism(@enzyme_mechanism begin
            substrates: A, B
            products: P, Q
            steps: begin
                E + A ⇌ E(A)
                E(B) + A ⇌ E(A, B)
                E + B ⇌ E(B)
                E(A) + B ⇌ E(A, B)
                (E + P ⇌ E(P), E(Q) + P ⇌ E(P, Q))
                (E + Q ⇌ E(Q), E(P) + Q ⇌ E(P, Q))
                E(A, B) <--> E(P, Q)
            end
        end)
        splitPQ = EnzymeRates.Mechanism(@enzyme_mechanism begin
            substrates: A, B
            products: P, Q
            steps: begin
                (E + A ⇌ E(A), E(B) + A ⇌ E(A, B))
                (E + B ⇌ E(B), E(A) + B ⇌ E(A, B))
                E + P ⇌ E(P)
                E(Q) + P ⇌ E(P, Q)
                E + Q ⇌ E(Q)
                E(P) + Q ⇌ E(P, Q)
                E(A, B) <--> E(P, Q)
            end
        end)
        @test EnzymeRates._independent_param_count(m) == 5
        kids = EnzymeRates._expand_split_kinetic_group(m)
        @test length(kids) == 2
        @test Set(kids) == Set([splitAB, splitPQ])
        for c in kids
            @test EnzymeRates._independent_param_count(c) == 6
        end
    end

    @testset "_expand_split_kinetic_group: bi-bi with dead ends splits each square" begin
        # Four four-cycles run through this graph — the A/B and P/Q catalytic
        # squares and the A/P and B/Q dead-end squares — and each gives one
        # child that frees the interaction factor of its doubly-bound form.
        m = EnzymeRates.Mechanism(@enzyme_mechanism begin
            substrates: A, B
            products: P, Q
            steps: begin
                (E + A ⇌ E(A), E(B) + A ⇌ E(A, B), E(P) + A ⇌ E(A, P))
                (E + B ⇌ E(B), E(A) + B ⇌ E(A, B), E(Q) + B ⇌ E(B, Q))
                (E + P ⇌ E(P), E(Q) + P ⇌ E(P, Q), E(A) + P ⇌ E(A, P))
                (E + Q ⇌ E(Q), E(P) + Q ⇌ E(P, Q), E(B) + Q ⇌ E(B, Q))
                E(A, B) <--> E(P, Q)
            end
        end)
        splitAP = EnzymeRates.Mechanism(@enzyme_mechanism begin
            substrates: A, B
            products: P, Q
            steps: begin
                (E + A ⇌ E(A), E(B) + A ⇌ E(A, B))
                E(P) + A ⇌ E(A, P)
                (E + B ⇌ E(B), E(A) + B ⇌ E(A, B), E(Q) + B ⇌ E(B, Q))
                (E + P ⇌ E(P), E(Q) + P ⇌ E(P, Q))
                E(A) + P ⇌ E(A, P)
                (E + Q ⇌ E(Q), E(P) + Q ⇌ E(P, Q), E(B) + Q ⇌ E(B, Q))
                E(A, B) <--> E(P, Q)
            end
        end)
        splitAB = EnzymeRates.Mechanism(@enzyme_mechanism begin
            substrates: A, B
            products: P, Q
            steps: begin
                (E + A ⇌ E(A), E(P) + A ⇌ E(A, P))
                E(B) + A ⇌ E(A, B)
                (E + B ⇌ E(B), E(Q) + B ⇌ E(B, Q))
                E(A) + B ⇌ E(A, B)
                (E + P ⇌ E(P), E(Q) + P ⇌ E(P, Q), E(A) + P ⇌ E(A, P))
                (E + Q ⇌ E(Q), E(P) + Q ⇌ E(P, Q), E(B) + Q ⇌ E(B, Q))
                E(A, B) <--> E(P, Q)
            end
        end)
        splitBQ = EnzymeRates.Mechanism(@enzyme_mechanism begin
            substrates: A, B
            products: P, Q
            steps: begin
                (E + A ⇌ E(A), E(B) + A ⇌ E(A, B), E(P) + A ⇌ E(A, P))
                (E + B ⇌ E(B), E(A) + B ⇌ E(A, B))
                E(Q) + B ⇌ E(B, Q)
                (E + P ⇌ E(P), E(Q) + P ⇌ E(P, Q), E(A) + P ⇌ E(A, P))
                (E + Q ⇌ E(Q), E(P) + Q ⇌ E(P, Q))
                E(B) + Q ⇌ E(B, Q)
                E(A, B) <--> E(P, Q)
            end
        end)
        splitPQ = EnzymeRates.Mechanism(@enzyme_mechanism begin
            substrates: A, B
            products: P, Q
            steps: begin
                (E + A ⇌ E(A), E(B) + A ⇌ E(A, B), E(P) + A ⇌ E(A, P))
                (E + B ⇌ E(B), E(A) + B ⇌ E(A, B), E(Q) + B ⇌ E(B, Q))
                (E + P ⇌ E(P), E(A) + P ⇌ E(A, P))
                E(Q) + P ⇌ E(P, Q)
                (E + Q ⇌ E(Q), E(B) + Q ⇌ E(B, Q))
                E(P) + Q ⇌ E(P, Q)
                E(A, B) <--> E(P, Q)
            end
        end)
        @test EnzymeRates._independent_param_count(m) == 5
        kids = EnzymeRates._expand_split_kinetic_group(m)
        @test length(kids) == 4
        @test Set(kids) == Set([splitAP, splitAB, splitBQ, splitPQ])
        for c in kids
            @test EnzymeRates._independent_param_count(c) == 6
        end
    end

    @testset "_expand_split_kinetic_group: ter-ter cube splits a substrate pair" begin
        # Each child says "the affinity for one substrate depends on whether the
        # other is bound", a 2 + 2 division of two four-step groups that no carve
        # of a single step could produce. One child per substrate pair.
        m = EnzymeRates.Mechanism(@enzyme_mechanism begin
            substrates: A, B, C
            products: P, Q, R
            steps: begin
                (E + A ⇌ E(A), E(B) + A ⇌ E(A, B), E(C) + A ⇌ E(A, C),
                 E(B, C) + A ⇌ E(A, B, C))
                (E + B ⇌ E(B), E(A) + B ⇌ E(A, B), E(C) + B ⇌ E(B, C),
                 E(A, C) + B ⇌ E(A, B, C))
                (E + C ⇌ E(C), E(A) + C ⇌ E(A, C), E(B) + C ⇌ E(B, C),
                 E(A, B) + C ⇌ E(A, B, C))
                E(A, B, C) <--> E(P, Q, R)
                E(Q, R) + P ⇌ E(P, Q, R)
                E(R) + Q ⇌ E(Q, R)
                E + R ⇌ E(R)
            end
        end)
        splitAB = EnzymeRates.Mechanism(@enzyme_mechanism begin
            substrates: A, B, C
            products: P, Q, R
            steps: begin
                (E + A ⇌ E(A), E(C) + A ⇌ E(A, C))
                (E(B) + A ⇌ E(A, B), E(B, C) + A ⇌ E(A, B, C))
                (E + B ⇌ E(B), E(C) + B ⇌ E(B, C))
                (E(A) + B ⇌ E(A, B), E(A, C) + B ⇌ E(A, B, C))
                (E + C ⇌ E(C), E(A) + C ⇌ E(A, C), E(B) + C ⇌ E(B, C),
                 E(A, B) + C ⇌ E(A, B, C))
                E(A, B, C) <--> E(P, Q, R)
                E(Q, R) + P ⇌ E(P, Q, R)
                E(R) + Q ⇌ E(Q, R)
                E + R ⇌ E(R)
            end
        end)
        splitAC = EnzymeRates.Mechanism(@enzyme_mechanism begin
            substrates: A, B, C
            products: P, Q, R
            steps: begin
                (E + A ⇌ E(A), E(B) + A ⇌ E(A, B))
                (E(C) + A ⇌ E(A, C), E(B, C) + A ⇌ E(A, B, C))
                (E + B ⇌ E(B), E(A) + B ⇌ E(A, B), E(C) + B ⇌ E(B, C),
                 E(A, C) + B ⇌ E(A, B, C))
                (E + C ⇌ E(C), E(B) + C ⇌ E(B, C))
                (E(A) + C ⇌ E(A, C), E(A, B) + C ⇌ E(A, B, C))
                E(A, B, C) <--> E(P, Q, R)
                E(Q, R) + P ⇌ E(P, Q, R)
                E(R) + Q ⇌ E(Q, R)
                E + R ⇌ E(R)
            end
        end)
        splitBC = EnzymeRates.Mechanism(@enzyme_mechanism begin
            substrates: A, B, C
            products: P, Q, R
            steps: begin
                (E + A ⇌ E(A), E(B) + A ⇌ E(A, B), E(C) + A ⇌ E(A, C),
                 E(B, C) + A ⇌ E(A, B, C))
                (E + B ⇌ E(B), E(A) + B ⇌ E(A, B))
                (E(C) + B ⇌ E(B, C), E(A, C) + B ⇌ E(A, B, C))
                (E + C ⇌ E(C), E(A) + C ⇌ E(A, C))
                (E(B) + C ⇌ E(B, C), E(A, B) + C ⇌ E(A, B, C))
                E(A, B, C) <--> E(P, Q, R)
                E(Q, R) + P ⇌ E(P, Q, R)
                E(R) + Q ⇌ E(Q, R)
                E + R ⇌ E(R)
            end
        end)
        @test EnzymeRates._independent_param_count(m) == 7
        kids = EnzymeRates._expand_split_kinetic_group(m)
        @test length(kids) == 3
        @test Set(kids) == Set([splitAB, splitAC, splitBC])
        for c in kids
            @test EnzymeRates._independent_param_count(c) == 8
        end
    end

    @testset "_expand_split_kinetic_group: bi-bi seeds emit 102 children" begin
        seeds = EnzymeRates.init_mechanisms(_testhelper_bibi_rxn)
        @test sum(length(EnzymeRates._expand_split_kinetic_group(m)) for m in seeds) == 102
    end

    @testset "_expand_split_kinetic_group: a rejected split is the parent's model" begin
        # Level-1 candidates the count test rejects have the parent's identifiable
        # rank; the rendered equation may differ only in which tied name survives.
        m = EnzymeRates.Mechanism(@enzyme_mechanism begin
            substrates: A, B
            products: P, Q
            steps: begin
                (E + A ⇌ E(A), E(B) + A ⇌ E(A, B), E(P) + A ⇌ E(A, P))
                (E + B ⇌ E(B), E(A) + B ⇌ E(A, B), E(Q) + B ⇌ E(B, Q))
                (E + P ⇌ E(P), E(Q) + P ⇌ E(P, Q), E(A) + P ⇌ E(A, P))
                (E + Q ⇌ E(Q), E(P) + Q ⇌ E(P, Q), E(B) + Q ⇌ E(B, Q))
                E(A, B) <--> E(P, Q)
            end
        end)
        counter = EnzymeRates._partition_independent_count(m)
        flat = EnzymeRates._flat_steps(m)
        pos = Dict(s => j for (j, (s, _)) in enumerate(flat))
        ids0 = [g for (_, g) in flat]
        base = counter(ids0)
        r0 = _testhelper_identifiable_rank(m)
        @test r0 == base
        rejected = 0
        for (g, grp) in enumerate(EnzymeRates.steps(m)),
            bp in EnzymeRates._context_bipartitions(grp)
            ids = copy(ids0)
            for s in bp[2]; ids[pos[s]] = length(EnzymeRates.steps(m)) + 1; end
            counter(ids) > base && continue
            rejected += 1
            child = EnzymeRates._apply_bipartitions(m, [(g, bp)])
            @test _testhelper_identifiable_rank(child) == r0
        end
        @test rejected > 0
    end

    @testset "_expand_split_kinetic_group: a four-step group splits two and two" begin
        m = EnzymeRates.Mechanism(@enzyme_mechanism begin
            substrates: A, B
            products: P, Q
            steps: begin
                (E + A ⇌ E(A), E(B) + A ⇌ E(A, B), E(P) + A ⇌ E(A, P),
                 E(B, P) + A ⇌ E(A, B, P))
                (E + B ⇌ E(B), E(A) + B ⇌ E(A, B), E(P) + B ⇌ E(B, P),
                 E(A, P) + B ⇌ E(A, B, P))
                (E + P ⇌ E(P), E(Q) + P ⇌ E(P, Q), E(A) + P ⇌ E(A, P),
                 E(B) + P ⇌ E(B, P), E(A, B) + P ⇌ E(A, B, P))
                (E + Q ⇌ E(Q), E(P) + Q ⇌ E(P, Q))
                E(A, B) <--> E(P, Q)
            end
        end)
        binder(grp) = EnzymeRates.name(EnzymeRates.bound_metabolite(first(grp)))
        a_group = only(grp for grp in EnzymeRates.steps(m)
                       if length(grp) == 4 && binder(grp) == :A)
        bps = EnzymeRates._context_bipartitions(a_group)
        # Context B and context P each divide the four A steps two and two, a
        # division no carve of a single step can produce.
        @test count(bp -> length(bp[1]) == 2 && length(bp[2]) == 2, bps) == 2
        gi = findfirst(==(a_group), EnzymeRates.steps(m))
        even_bp = first(bp for bp in bps if length(bp[1]) == 2)
        child = EnzymeRates._apply_bipartitions(m, [(gi, even_bp)])
        @test count(grp -> length(grp) == 2 && Set(grp) ⊆ Set(a_group),
                    EnzymeRates.steps(child)) == 2
    end

    @testset "_expand_split_kinetic_group: allosteric parent" begin
        am = EnzymeRates.AllostericMechanism(@allosteric_mechanism begin
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
        end)
        base = EnzymeRates._independent_param_count(am)
        kids = EnzymeRates._expand_split_kinetic_group(am)
        @test !isempty(kids)
        for c in kids
            @test c isa EnzymeRates.AllostericMechanism
            @test EnzymeRates._independent_param_count(c) > base
            @test length(EnzymeRates.cat_allo_states(c)) == length(EnzymeRates.steps(c))
        end
    end

    @testset "_expand_split_kinetic_group: ter-ter random-order seed within budget" begin
        terter = @enzyme_reaction begin
            substrates: A[C], B[N], C[O]
            products: P[C], Q[N], R[O]
        end
        seeds = EnzymeRates.init_mechanisms(terter)
        nst = [EnzymeRates.n_steps(m) for m in seeds]
        worst = seeds[argmax(nst)]
        @test EnzymeRates.n_steps(worst) == 55
        t = @elapsed kids = EnzymeRates._expand_split_kinetic_group(worst)
        @test length(kids) == 12
        @test t < 60
    end
end

@testset "_expand_split_kinetic_group: by conformation" begin
    # S binds both conformations in one kinetic group; Estar(S) is a dead end
    # whose dissociation constant no cycle ties, so dividing the group by
    # conformation frees one parameter. The other groups are singletons.
    m = EnzymeRates.Mechanism(@enzyme_mechanism begin
        substrates: S
        products: P
        steps: begin
            (E + S ⇌ E(S), Estar + S ⇌ Estar(S))
            E(S) <--> Estar(P)
            Estar(P) ⇌ Estar + P
            Estar <--> E
        end
    end)
    by_conformation = EnzymeRates.Mechanism(@enzyme_mechanism begin
        substrates: S
        products: P
        steps: begin
            E + S ⇌ E(S)
            Estar + S ⇌ Estar(S)
            E(S) <--> Estar(P)
            Estar(P) ⇌ Estar + P
            Estar <--> E
        end
    end)
    kids = EnzymeRates._expand_split_kinetic_group(m)
    @test length(kids) == 1
    @test Set(kids) == Set([by_conformation])
end

@testset "_expand_split_kinetic_group: by residual" begin
    # Q binds free E on the cycle and the two modified forms as dead ends. The
    # Q group divides by B ({E, E(; res)} | {E(B; res)}) and by residual
    # ({E} | {E(; res), E(B; res)}); each division frees a dead-end constant
    # on its own, and only one bipartition per group is ever applied.
    m = EnzymeRates.Mechanism(@enzyme_mechanism begin
        substrates: A, B
        products: P, Q
        steps: begin
            (E + Q ⇌ E(Q), E(; residual = A - P) + Q ⇌ E(Q; residual = A - P),
             E(B; residual = A - P) + Q ⇌ E(B, Q; residual = A - P))
            E + A ⇌ E(A)
            E(A) <--> E(P; residual = A - P)
            E(P; residual = A - P) ⇌ E(; residual = A - P) + P
            E(; residual = A - P) + B ⇌ E(B; residual = A - P)
            E(B; residual = A - P) <--> E(Q)
        end
    end)
    by_B = EnzymeRates.Mechanism(@enzyme_mechanism begin
        substrates: A, B
        products: P, Q
        steps: begin
            (E + Q ⇌ E(Q), E(; residual = A - P) + Q ⇌ E(Q; residual = A - P))
            E(B; residual = A - P) + Q ⇌ E(B, Q; residual = A - P)
            E + A ⇌ E(A)
            E(A) <--> E(P; residual = A - P)
            E(P; residual = A - P) ⇌ E(; residual = A - P) + P
            E(; residual = A - P) + B ⇌ E(B; residual = A - P)
            E(B; residual = A - P) <--> E(Q)
        end
    end)
    by_residual = EnzymeRates.Mechanism(@enzyme_mechanism begin
        substrates: A, B
        products: P, Q
        steps: begin
            E + Q ⇌ E(Q)
            (E(; residual = A - P) + Q ⇌ E(Q; residual = A - P),
             E(B; residual = A - P) + Q ⇌ E(B, Q; residual = A - P))
            E + A ⇌ E(A)
            E(A) <--> E(P; residual = A - P)
            E(P; residual = A - P) ⇌ E(; residual = A - P) + P
            E(; residual = A - P) + B ⇌ E(B; residual = A - P)
            E(B; residual = A - P) <--> E(Q)
        end
    end)
    kids = EnzymeRates._expand_split_kinetic_group(m)
    @test length(kids) == 2
    @test Set(kids) == Set([by_B, by_residual])
end

@testset "expand_mechanisms: ter-ter random-order seed within budget" begin
    # Aggregate pin over the seed set: the seed with the most steps (55) is the
    # enumeration's worst case; measured 26 s for all seven moves in a cold
    # focused run, JIT included.
    terter = @enzyme_reaction begin
        substrates: A[C], B[N], C[O]
        products: P[C], Q[N], R[O]
    end
    seeds = EnzymeRates.init_mechanisms(terter)
    worst = seeds[argmax([EnzymeRates.n_steps(m) for m in seeds])]
    @test EnzymeRates.n_steps(worst) == 55
    t = @elapsed kids = EnzymeRates.expand_mechanisms([worst], terter)
    @test length(kids) == 81
    @test t < 120
end

@testset "expand_mechanisms rejects chemistry folded into a release step" begin
    # The moves take the isomerization step as the chemistry step, which is
    # how the enumerator writes every mechanism. A mechanism written for the
    # derivation with chemistry folded into a release (the ping-pong docs
    # page) is not a valid parent.
    rxn = @enzyme_reaction begin
        substrates: A[CX], B[N]
        products: P[C], Q[NX]
    end
    folded = EnzymeRates.Mechanism(@enzyme_mechanism begin
        substrates: A, B
        products: P, Q
        steps: begin
            E + A ⇌ E(A)
            E(A) <--> E(; residual = A - P) + P
            E(; residual = A - P) + B ⇌ E(B; residual = A - P)
            E(B; residual = A - P) <--> E + Q
        end
    end)
    err = try
        EnzymeRates.expand_mechanisms([folded], rxn); nothing
    catch e
        e
    end
    @test err isa ErrorException &&
          occursin("changes the covalent residual", sprint(showerror, err))
    # Aggregate pin over the ping-pong seed set: the enumerator itself never
    # writes a binding step that changes the residual.
    @test all(EnzymeRates._assert_chemistry_is_iso(m) === nothing
              for m in EnzymeRates.init_mechanisms(rxn))
    # A binding step may change conformation; only the residual is chemistry.
    conf = EnzymeRates.Mechanism(@enzyme_mechanism begin
        substrates: A, B
        products: P, Q
        steps: begin
            E + A ⇌ Estar(A)
            Estar(A) + B ⇌ Estar(A, B)
            Estar(A, B) <--> Estar(P, Q)
            Estar(P, Q) ⇌ Estar(P) + Q
            Estar(P) ⇌ E + P
        end
    end)
    @test EnzymeRates.expand_mechanisms([conf], rxn) isa Vector
end

@testset "_expand_add_dead_end_regulator: inhibitor mirrors one half-reaction" begin
    # The inhibitor binds the forms that bind a competing ligand, never a form
    # already carrying one, so a competing substrate's binding step is never
    # mirrored and no inhibitor-bound branch completes the net reaction. With A
    # binding E(Q) and P binding E(B; res) as dead ends, the pattern "compete
    # with A and P" puts I on E, E(Q), E(; res) and E(B; res): the second
    # half-reaction (B binds, chemistry, Q leaves) runs with I bound, and A
    # cannot bind E(I). Six patterns give distinct targets.
    rxn = @enzyme_reaction begin
        substrates: A[CX], B[N]
        products: P[C], Q[NX]
        competitive_inhibitors: I
    end
    m = EnzymeRates.Mechanism(@enzyme_mechanism begin
        substrates: A, B
        products: P, Q
        steps: begin
            E + A ⇌ E(A)
            E(A) <--> E(P; residual = A - P)
            E(P; residual = A - P) ⇌ E(; residual = A - P) + P
            E(; residual = A - P) + B ⇌ E(B; residual = A - P)
            E(B; residual = A - P) <--> E(Q)
            E(Q) ⇌ E + Q
            E(Q) + A ⇌ E(A, Q)
            E(B; residual = A - P) + P ⇌ E(B, P; residual = A - P)
        end
    end)
    mech(block) = EnzymeRates.Mechanism(block)
    second_half = mech(@enzyme_mechanism begin      # I competes with A and P
        substrates: A, B
        products: P, Q
        regulators: I
        steps: begin
            E + A ⇌ E(A)
            E(A) <--> E(P; residual = A - P)
            E(P; residual = A - P) ⇌ E(; residual = A - P) + P
            (E(; residual = A - P) + B ⇌ E(B; residual = A - P),
             E(I::Inh; residual = A - P) + B ⇌ E(B, I::Inh; residual = A - P))
            (E(B; residual = A - P) <--> E(Q),
             E(B, I::Inh; residual = A - P) <--> E(I::Inh, Q))
            (E(Q) ⇌ E + Q, E(I::Inh, Q) ⇌ E(I::Inh) + Q)
            E(Q) + A ⇌ E(A, Q)
            E(B; residual = A - P) + P ⇌ E(B, P; residual = A - P)
            (E + I ⇌ E(I::Inh), E(Q) + I ⇌ E(I::Inh, Q),
             E(; residual = A - P) + I ⇌ E(I::Inh; residual = A - P),
             E(B; residual = A - P) + I ⇌ E(B, I::Inh; residual = A - P))
        end
    end)
    only_E = mech(@enzyme_mechanism begin           # I competes with A and Q
        substrates: A, B
        products: P, Q
        regulators: I
        steps: begin
            E + A ⇌ E(A)
            E(A) <--> E(P; residual = A - P)
            E(P; residual = A - P) ⇌ E(; residual = A - P) + P
            E(; residual = A - P) + B ⇌ E(B; residual = A - P)
            E(B; residual = A - P) <--> E(Q)
            E(Q) ⇌ E + Q
            E(Q) + A ⇌ E(A, Q)
            E(B; residual = A - P) + P ⇌ E(B, P; residual = A - P)
            E + I ⇌ E(I::Inh)
        end
    end)
    B_binding = mech(@enzyme_mechanism begin        # I competes with A, P and Q
        substrates: A, B
        products: P, Q
        regulators: I
        steps: begin
            E + A ⇌ E(A)
            E(A) <--> E(P; residual = A - P)
            E(P; residual = A - P) ⇌ E(; residual = A - P) + P
            (E(; residual = A - P) + B ⇌ E(B; residual = A - P),
             E(I::Inh; residual = A - P) + B ⇌ E(B, I::Inh; residual = A - P))
            E(B; residual = A - P) <--> E(Q)
            E(Q) ⇌ E + Q
            E(Q) + A ⇌ E(A, Q)
            E(B; residual = A - P) + P ⇌ E(B, P; residual = A - P)
            (E + I ⇌ E(I::Inh),
             E(; residual = A - P) + I ⇌ E(I::Inh; residual = A - P),
             E(B; residual = A - P) + I ⇌ E(B, I::Inh; residual = A - P))
        end
    end)
    only_Eres = mech(@enzyme_mechanism begin        # I competes with B and P
        substrates: A, B
        products: P, Q
        regulators: I
        steps: begin
            E + A ⇌ E(A)
            E(A) <--> E(P; residual = A - P)
            E(P; residual = A - P) ⇌ E(; residual = A - P) + P
            E(; residual = A - P) + B ⇌ E(B; residual = A - P)
            E(B; residual = A - P) <--> E(Q)
            E(Q) ⇌ E + Q
            E(Q) + A ⇌ E(A, Q)
            E(B; residual = A - P) + P ⇌ E(B, P; residual = A - P)
            E(; residual = A - P) + I ⇌ E(I::Inh; residual = A - P)
        end
    end)
    E_and_Eres = mech(@enzyme_mechanism begin       # I competes with B and Q
        substrates: A, B
        products: P, Q
        regulators: I
        steps: begin
            E + A ⇌ E(A)
            E(A) <--> E(P; residual = A - P)
            E(P; residual = A - P) ⇌ E(; residual = A - P) + P
            E(; residual = A - P) + B ⇌ E(B; residual = A - P)
            E(B; residual = A - P) <--> E(Q)
            E(Q) ⇌ E + Q
            E(Q) + A ⇌ E(A, Q)
            E(B; residual = A - P) + P ⇌ E(B, P; residual = A - P)
            (E + I ⇌ E(I::Inh), E(; residual = A - P) + I ⇌ E(I::Inh; residual = A - P))
        end
    end)
    Q_release = mech(@enzyme_mechanism begin        # I competes with A, B and P
        substrates: A, B
        products: P, Q
        regulators: I
        steps: begin
            E + A ⇌ E(A)
            E(A) <--> E(P; residual = A - P)
            E(P; residual = A - P) ⇌ E(; residual = A - P) + P
            E(; residual = A - P) + B ⇌ E(B; residual = A - P)
            E(B; residual = A - P) <--> E(Q)
            (E(Q) ⇌ E + Q, E(I::Inh, Q) ⇌ E(I::Inh) + Q)
            E(Q) + A ⇌ E(A, Q)
            E(B; residual = A - P) + P ⇌ E(B, P; residual = A - P)
            (E + I ⇌ E(I::Inh), E(Q) + I ⇌ E(I::Inh, Q),
             E(; residual = A - P) + I ⇌ E(I::Inh; residual = A - P))
        end
    end)
    expected = [second_half, only_E, B_binding, only_Eres, E_and_Eres, Q_release]
    kids = EnzymeRates._expand_add_dead_end_regulator(m, rxn)
    @test length(kids) == 6
    @test Set(EnzymeRates.steps.(kids)) == Set(EnzymeRates.steps.(expected))
    @test all(c -> EnzymeRates.reaction(c) ==
                   EnzymeRates._add_competitive_inhibitor(rxn, :I), kids)
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
    inhibitor_bound(sp) =
        any(b -> b isa EnzymeRates.CompetitiveInhibitor, EnzymeRates.bound(sp))
    in_branch(s) = inhibitor_bound(EnzymeRates.from_species(s)) &&
                   inhibitor_bound(EnzymeRates.to_species(s))
    n_children = 0; n_half = 0
    seeds = EnzymeRates.init_mechanisms(rxn)
    for m in seeds, c in EnzymeRates._expand_add_dead_end_regulator(m, rxn)
        branch = [s for grp in EnzymeRates.steps(c) for s in grp if in_branch(s)]
        bound_in_branch = Set(EnzymeRates.name(EnzymeRates.bound_metabolite(s))
                              for s in branch
                              if EnzymeRates.bound_metabolite(s) !== nothing)
        @test !(Set([:A, :B]) ⊆ bound_in_branch)
        n_children += 1
        any(EnzymeRates.is_iso, branch) && (n_half += 1)
    end
    @test n_children > 0
    @test n_half > 0          # the half-reaction mirror really occurs on the seeds
end

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

@testset "to_allosteric shares :EqualAI downstream constants across conformations" begin
    # An inactive-conformation graph (`_state_mechanism(am, :I)`) drops each
    # `:OnlyA` binding step, so that ligand's downstream complex loses its
    # binding-in edge and is reached only by conformational flip / reverse
    # catalysis. It is still a bound form, never free enzyme, so `_free_enz_set`
    # excludes it — keeping a lumped `:EqualAI` group's naming rep identical in
    # both conformations. Without that exclusion the inactive-state rep flips
    # (e.g. `K_Lactate_ENADH` vs the active-state `K_Lactate_ENAD`), un-lumping a
    # shared constant so the combined solve fails to merge it and over-counts.
    #
    # `_expand_to_allosteric` emits dead-inactive `:OnlyA`-binding combos (all
    # chemical steps `:OnlyA`), so the case this guards is the child that sets
    # the free-enzyme bindings (a substrate/product pair on free E) `:OnlyA`,
    # leaving the lumped downstream group's binding `:EqualAI` — reached in the
    # inactive conformation only by flip. That group's shared constant stays
    # merged, so the child adds only `L` — parent + 1 fitted params. A child
    # that instead sets a downstream binding `:OnlyA` orphans a bound form of
    # the lumped group; its multi-`:OnlyA` derivation is a separate, documented
    # concern and is not asserted here.
    rxn = @enzyme_reaction begin
        substrates: NADH[C21H29N7O14P2], Pyruvate[C3H4O3]
        products:   Lactate[C3H6O3], NAD[C21H27N7O14P2]
        oligomeric_state: 4
    end
    # A mechanism exercises the path only when a downstream kinetic group lumps
    # binding to ≥2 bound complexes (no free-E member) — the split the `:OnlyA`
    # orphaning would otherwise break.
    has_lumped_bound_group(m) = any(EnzymeRates.steps(m)) do g
        length(g) ≥ 2 &&
            all(s -> !isempty(EnzymeRates.bound(EnzymeRates.from_species(s))), g)
    end
    is_free_e_binding(m, g) = begin
        rs = EnzymeRates.rep_step(m, g)
        isempty(EnzymeRates.bound(EnzymeRates.from_species(rs))) &&
            EnzymeRates.bound_metabolite(rs) !== nothing
    end
    n_reproducers = 0
    for m in EnzymeRates.init_mechanisms(rxn)
        has_lumped_bound_group(m) || continue
        pn = length(EnzymeRates.fitted_params(EnzymeRates.compile_mechanism(m)))
        free_e = Set(g for g in EnzymeRates.kinetic_groups(m)
                     if is_free_e_binding(m, g))
        iso = Set(g for g in 1:length(EnzymeRates.steps(m))
                  if EnzymeRates.is_iso(EnzymeRates.steps(m)[g][1]))
        children = EnzymeRates._expand_to_allosteric(m, rxn)
        # Every emitted child is Haldane-valid (Task-7-safe) and compiles.
        for c in children
            @test EnzymeRates._onlya_haldane_violation(
                EnzymeRates.reaction(c), EnzymeRates.steps(c),
                EnzymeRates.cat_allo_states(c)) === nothing
            @test EnzymeRates.compile_mechanism(c) isa AllostericEnzymeMechanism
        end
        # The free-enzyme-normalization child (the free-E bindings :OnlyA, plus
        # all chemical steps :OnlyA — a dead inactive conformation) leaves the
        # lumped downstream group's binding :EqualAI, so its shared constant
        # stays merged → Δ = 1.
        clean = filter(children) do c
            Set(g for g in 1:length(EnzymeRates.cat_allo_states(c))
                if EnzymeRates.cat_allo_states(c)[g] == :OnlyA) == union(free_e, iso)
        end
        @test length(clean) == 1
        @test length(EnzymeRates.fitted_params(
            EnzymeRates.compile_mechanism(only(clean)))) == pn + 1
        n_reproducers += 1
        n_reproducers ≥ 3 && break
    end
    @test n_reproducers ≥ 1
end

@testset "seed_mechanisms" begin
    # Tag-stripped skeleton: normalize every regulator state to :OnlyA so the
    # 2^n one-ligand-site state assignments of a lineage collapse to one key.
    # The catalytic backbone and cat_allo_states are preserved.
    skeleton_key(m) = hash(EnzymeRates.AllostericMechanism(
        EnzymeRates.reaction(m), EnzymeRates.steps(m),
        EnzymeRates.cat_allo_states(m), EnzymeRates.catalytic_multiplicity(m),
        [EnzymeRates.RegulatorySite(EnzymeRates.ligands(s),
             EnzymeRates.multiplicity(s),
             fill(:OnlyA, length(EnzymeRates.allo_states(s))))
         for s in EnzymeRates.regulatory_sites(m)]))

    # The allosteric state of ligand `reg` wherever it sits, or :none.
    function state_of(m, reg)
        for site in EnzymeRates.regulatory_sites(m)
            for (lig, st) in zip(EnzymeRates.ligands(site),
                                 EnzymeRates.allo_states(site))
                EnzymeRates.name(lig) == reg && return st
            end
        end
        :none
    end

    @testset "undesignated: skeletons × 2²" begin
        rxn = @enzyme_reaction begin
            substrates: A[C6H12O6]
            products:   B[C6H12O6]
            allosteric_regulators: X, Y
            oligomeric_state: 2
        end
        seeds = EnzymeRates.seed_mechanisms(rxn, Set([:X, :Y]), Set{Symbol}())
        @test !isempty(seeds)
        # Every seed is a fully-regulated, cheap-state AllostericMechanism with X
        # and Y each at their own single-ligand site.
        for s in seeds
            @test s isa EnzymeRates.AllostericMechanism
            @test EnzymeRates._bound_allo_regs(s) == Set([:X, :Y])
            sites = EnzymeRates.regulatory_sites(s)
            @test length(sites) == 2
            @test all(length(EnzymeRates.ligands(si)) == 1 for si in sites)
            @test !EnzymeRates._has_nonequalai(s)
        end
        # Each lineage contributes exactly 2² = 4 state assignments.
        skels = Dict{UInt64, Int}()
        for s in seeds
            k = skeleton_key(s)
            skels[k] = get(skels, k, 0) + 1
        end
        @test all(==(4), values(skels))
        @test length(seeds) == 4 * length(skels)
    end

    @testset "designated types collapse to ×1" begin
        rxn = @enzyme_reaction begin
            substrates: A[C6H12O6]
            products:   B[C6H12O6]
            allosteric_regulators: X, Y
            oligomeric_state: 2
        end
        rxn_d = @enzyme_reaction begin
            substrates: A[C6H12O6]
            products:   B[C6H12O6]
            allosteric_regulators: X::Activator, Y::Inhibitor
            oligomeric_state: 2
        end
        seeds = EnzymeRates.seed_mechanisms(rxn, Set([:X, :Y]), Set{Symbol}())
        seeds_d = EnzymeRates.seed_mechanisms(rxn_d, Set([:X, :Y]), Set{Symbol}())
        @test !isempty(seeds_d)
        # One state assignment per lineage: seed count equals the skeleton count.
        skels_d = Set(skeleton_key(s) for s in seeds_d)
        @test length(seeds_d) == length(skels_d)
        # The same catalytic-allostery skeletons as the undesignated build.
        @test length(seeds) == 4 * length(seeds_d)
        # An activator seeds as :OnlyA, an inhibitor as :OnlyI.
        for s in seeds_d
            @test state_of(s, :X) == :OnlyA
            @test state_of(s, :Y) == :OnlyI
        end
    end

    @testset "required vs optional" begin
        rxn = @enzyme_reaction begin
            substrates: A[C6H12O6]
            products:   B[C6H12O6]
            allosteric_regulators: X, Y, Z
            oligomeric_state: 2
        end
        seeds = EnzymeRates.seed_mechanisms(rxn, Set([:X, :Y]), Set{Symbol}())
        @test !isempty(seeds)
        for s in seeds
            @test issubset(Set([:X, :Y]), EnzymeRates._bound_allo_regs(s))
            @test !(:Z in EnzymeRates._bound_allo_regs(s))
        end
    end

    @testset "per-lineage floor invariant" begin
        rxn = @enzyme_reaction begin
            substrates: A[C6H12O6]
            products:   B[C6H12O6]
            allosteric_regulators: X, Y
            oligomeric_state: 2
        end
        seeds = EnzymeRates.seed_mechanisms(rxn, Set([:X, :Y]), Set{Symbol}())
        n_required = 2
        # A required allosteric regulator adds one dissociation constant; the two
        # required regulators additionally lift the mechanism to two
        # conformations, so a seed carries its base catalytic count + n_required
        # + 1 (L). Compile is slow; sample a few seeds.
        for s in seeds[1:min(3, length(seeds))]
            base = EnzymeRates.Mechanism(
                EnzymeRates.reaction(s), copy(EnzymeRates.steps(s)))
            base_count = length(
                EnzymeRates.fitted_params(EnzymeRates.compile_mechanism(base)))
            seed_count = length(
                EnzymeRates.fitted_params(EnzymeRates.compile_mechanism(s)))
            @test seed_count == base_count + n_required + 1
        end
    end

    @testset "competitive-only: non-allosteric seeds, no L" begin
        rxn = @enzyme_reaction begin
            substrates: A[C6H12O6]
            products:   B[C6H12O6]
            competitive_inhibitors: I
        end
        seeds = EnzymeRates.seed_mechanisms(rxn, Set{Symbol}(), Set([:I]))
        @test !isempty(seeds)
        # No allosteric-lifting move fires: every seed is a plain Mechanism that
        # binds I as a dead-end competitive binding, and carries no L.
        for s in seeds
            @test s isa EnzymeRates.Mechanism
            @test :I in EnzymeRates._bound_comp_inhibitors(s)
        end
        # Floor is base + n_required_comp (no L). Uni-uni has one catalytic
        # backbone, so every seed shares one base lineage.
        inits = EnzymeRates.init_mechanisms(rxn)
        @test length(inits) == 1
        base_count = length(
            EnzymeRates.fitted_params(EnzymeRates.compile_mechanism(inits[1])))
        for s in seeds[1:min(2, length(seeds))]
            fp = EnzymeRates.fitted_params(EnzymeRates.compile_mechanism(s))
            @test length(fp) == base_count + 1
            @test !(:L in fp)
        end
    end
end

@testset "seed_mechanisms wave-parallel equivalence" begin
    # Inline serial FIFO BFS = the reference the parallel version must match.
    function serial_seed_reference(rxn, req_allo, req_comp)
        visited = Set{UInt64}()
        queue = Union{EnzymeRates.Mechanism, EnzymeRates.AllostericMechanism}[]
        seeds = Union{EnzymeRates.Mechanism, EnzymeRates.AllostericMechanism}[]
        enq(m) = begin
            h = hash(m)
            h in visited && return
            push!(visited, h); push!(queue, m)
            EnzymeRates._binds_all_required(m, req_allo, req_comp) && push!(seeds, m)
        end
        for m in EnzymeRates.init_mechanisms(rxn); enq(m); end
        while !isempty(queue)
            m = popfirst!(queue)
            for c in EnzymeRates._seed_children(m, rxn, req_allo)
                EnzymeRates._is_seed_node(c, rxn, req_allo, req_comp) && enq(c)
            end
        end
        seeds
    end

    req = Set([:R])
    empty = Set{Symbol}()
    got = EnzymeRates.seed_mechanisms(uni_uni_allo_reg, req, empty)
    ref = serial_seed_reference(uni_uni_allo_reg, req, empty)

    @test got == ref                                   # same seeds, same order
    @test !isempty(got)                                # the case is non-trivial
    @test allunique(hash.(got))                        # no duplicate structures
    @test all(m -> EnzymeRates._binds_all_required(m, req, empty), got)
    @test EnzymeRates.seed_mechanisms(uni_uni_allo_reg, req, empty) == got  # deterministic
end

# ─── _expand_to_allosteric Haldane closure ─────────────────────────────
@testset "_expand_to_allosteric emits only Haldane-valid mechanisms" begin
    ER = EnzymeRates
    # Uni-uni S ⇌ P: 3 kinetic groups (S binding, chemical step, P binding).
    # No allosteric regulator is declared, so the chemical group's own bare
    # :OnlyA variant is not emitted; every child therefore starts from a
    # binding promotion, which alone leaves the Haldane unsatisfiable.
    rxn = @enzyme_reaction begin
        substrates: S[C]
        products: P[C]
        oligomeric_state: 2
    end
    m = ER.Mechanism(@enzyme_mechanism begin
        substrates: S
        products: P
        steps: begin
            E + S ⇌ E(S)
            E(S) <--> E(P)
            E + P ⇌ E(P)
        end
    end)
    result = ER._expand_to_allosteric(m, rxn)

    @test !isempty(result)
    for r in result
        @test ER._onlya_haldane_violation(
            ER.reaction(r), ER.steps(r), ER.cat_allo_states(r)) === nothing
        # the move's purpose: every child still carries an :OnlyA
        @test :OnlyA in ER.cat_allo_states(r)
        ER._assert_mechanism_invariants(r)
        @test ER.compile_mechanism(r) isa AllostericEnzymeMechanism
    end
    @test length(unique(result)) == length(result)

    # Read tags by what each group binds, never by group index: the
    # AllostericMechanism constructor canonicalizes group order.
    function shape(x)
        tags = ER.cat_allo_states(x)
        by_ligand = Dict{Symbol, Symbol}()
        for (g, grp) in enumerate(ER.steps(x))
            bm = ER.bound_metabolite(grp[1])
            by_ligand[bm === nothing ? :chem : ER.name(bm)] = tags[g]
        end
        (by_ligand[:S], by_ligand[:chem], by_ligand[:P])
    end
    shapes = Set(shape(r) for r in result)

    # Every non-empty binding subset is emitted dead-inactive (the chemical
    # step :OnlyA): {S}, {P}, {S, P}. Three shapes, all with the chemical step
    # :OnlyA.
    @test (:OnlyA, :OnlyA, :EqualAI) in shapes    # {S} binding :OnlyA
    @test (:EqualAI, :OnlyA, :OnlyA) in shapes    # {P} binding :OnlyA
    @test (:OnlyA, :OnlyA, :OnlyA) in shapes      # {S, P} bindings :OnlyA
    @test length(shapes) == 3
end

# Tag of each catalytic kinetic group, keyed by the metabolite it binds
# (`:chem` for the chemical step). Group index is never a stable key: the
# AllostericMechanism constructor canonicalizes group order.
tags_by_ligand(x) = Dict(
    (bm = EnzymeRates.bound_metabolite(grp[1]);
     bm === nothing ? :chem : EnzymeRates.name(bm)) =>
        EnzymeRates.cat_allo_states(x)[g]
    for (g, grp) in enumerate(EnzymeRates.steps(x)))

haldane_ok(x) = EnzymeRates._onlya_haldane_violation(
    EnzymeRates.reaction(x), EnzymeRates.steps(x),
    EnzymeRates.cat_allo_states(x)) === nothing

# ─── _expand_change_allo_state Haldane filter ──────────────────────────
@testset "_expand_change_allo_state emits only Haldane-valid mechanisms" begin
    ER = EnzymeRates

    # Uni-uni S ⇌ P. The S binding's :OnlyA is one-sided; only the :OnlyA
    # chemical step (`k_I = 0`) makes the parent legal. Relaxing that
    # chemical step to :NonequalAI restores a finite `k_I` and strands the
    # binding, so that child has no thermodynamic reading and is dropped.
    parent = ER.AllostericMechanism(@allosteric_mechanism begin
        substrates: S
        products: P
        catalytic_multiplicity: 2
        catalytic_steps: begin
            E + S ⇌ E(S)      :: OnlyA
            E(S) <--> E(P)    :: OnlyA
            E + P ⇌ E(P)      :: EqualAI
        end
    end)
    @test haldane_ok(parent)

    result = ER._expand_change_allo_state(parent)
    for r in result
        @test haldane_ok(r)
        ER._assert_mechanism_invariants(r)
    end
    # Three tags are relaxable; the chemical step's relaxation is dropped.
    @test length(result) == 2
    @test !any(r -> tags_by_ligand(r)[:chem] == :NonequalAI, result)
    @test Set(tags_by_ligand(r) for r in result) == Set([
        Dict(:S => :OnlyA, :chem => :OnlyA, :P => :NonequalAI),
        Dict(:S => :NonequalAI, :chem => :OnlyA, :P => :EqualAI)])

    # Balanced parent: both bindings and the chemical step :OnlyA — the inactive
    # conformation binds nothing (a fully-inert T-state). Relaxing either binding
    # is retained (it leaves all chemical steps :OnlyA); relaxing the chemical
    # step is dropped (`_partial_onlya_catalysis`), because an :OnlyA binding
    # requires a catalytically-dead inactive conformation. That dropped
    # :NonequalAI-chemistry variant is anyway rate-equivalent to this fully-dead
    # form (the inactive binds nothing, so k_I is unobservable), so no hypothesis
    # is lost.
    balanced = ER.AllostericMechanism(@allosteric_mechanism begin
        substrates: S
        products: P
        catalytic_multiplicity: 2
        catalytic_steps: begin
            E + S ⇌ E(S)      :: OnlyA
            E(S) <--> E(P)    :: OnlyA
            E + P ⇌ E(P)      :: OnlyA
        end
    end)
    @test haldane_ok(balanced)
    kids = ER._expand_change_allo_state(balanced)
    @test length(kids) == 2           # the chemical-step relaxation is dropped
    @test !any(tags_by_ligand(k)[:chem] == :NonequalAI for k in kids)
    @test Set(tags_by_ligand(k) for k in kids) == Set([
        Dict(:S => :NonequalAI, :chem => :OnlyA, :P => :OnlyA),
        Dict(:S => :OnlyA, :chem => :OnlyA, :P => :NonequalAI)])

    # A regulatory ligand's tag is not an argument to the Haldane check —
    # a regulator site completes no catalytic cycle — so relaxing one can
    # never strand an :OnlyA binding, and that branch needs no filter.
    reg_parent = ER.AllostericMechanism(@allosteric_mechanism begin
        substrates: S
        products: P
        allosteric_regulators: R::OnlyA
        catalytic_multiplicity: 2
        catalytic_steps: begin
            E + S ⇌ E(S)      :: OnlyA
            E(S) <--> E(P)    :: OnlyA
            E + P ⇌ E(P)      :: EqualAI
        end
    end)
    @test haldane_ok(reg_parent)
    reg_kids = ER._expand_change_allo_state(reg_parent)
    # 2 surviving cat relaxations + the R-ligand relaxation, which survives.
    @test length(reg_kids) == 3
    @test count(k -> ER.allo_states(only(ER.regulatory_sites(k))) ==
                     [:NonequalAI], reg_kids) == 1
    for k in reg_kids
        @test haldane_ok(k)
    end
end

@testset ":OnlyA to_allosteric emits dead-inactive combos only" begin
    # 6-step iso-bearing ping-pong: 4 binding groups (ATP, F16BP, F6P, ADP) + 2
    # iso steps. Allostericise it and require every :OnlyA-tagged child to have
    # BOTH iso steps :OnlyA (no partial), and the bare (regulator-free) :OnlyA
    # child count to equal the valid binding subsets (2^4 - 1 = 15, all valid).
    m = EnzymeRates.Mechanism(@enzyme_mechanism begin
        substrates: ATP, F6P
        products: ADP, F16BP
        steps: begin
            E + ATP ⇌ E(ATP)
            E(ATP) <--> E(F16BP; residual = ATP - F16BP)
            E(; residual = ATP - F16BP) + F16BP ⇌ E(F16BP; residual = ATP - F16BP)
            E(; residual = ATP - F16BP) + F6P ⇌ E(F6P; residual = ATP - F16BP)
            E(F6P; residual = ATP - F16BP) ⇌ E(ADP)
            E + ADP ⇌ E(ADP)
        end
    end)
    rxn = EnzymeRates.reaction(m)
    kids = EnzymeRates._expand_to_allosteric(m, rxn)
    onlya = [k for k in kids if any(==(:OnlyA), EnzymeRates.cat_allo_states(k))]
    isochem(k) = [g for g in eachindex(EnzymeRates.steps(k))
                  if EnzymeRates.is_iso(EnzymeRates.steps(k)[g][1])]
    # every :OnlyA child is dead-inactive: all iso steps :OnlyA
    for k in onlya
        @test all(EnzymeRates.cat_allo_states(k)[g] === :OnlyA for g in isochem(k))
    end
    bare = [k for k in onlya if isempty(EnzymeRates.regulatory_sites(k))]
    @test length(bare) == 15
end

@testset "change_allo_state drops partial-catalysis relaxations" begin
    # Dead-inactive 6-step ping-pong: F6P binding + both iso steps :OnlyA. Relaxing
    # an iso step back would leave an :OnlyA binding beside a live chemical step (the
    # sink); the move must drop those while keeping the :OnlyA-binding relaxation.
    dead = @allosteric_mechanism begin
        substrates: ATP, F6P ; products: ADP, F16BP ; catalytic_multiplicity: 1
        catalytic_steps: begin
            E + ATP ⇌ E(ATP)                                                       :: EqualAI
            E(ATP) <--> E(F16BP; residual = ATP - F16BP)                           :: OnlyA
            E(; residual = ATP - F16BP) + F16BP ⇌ E(F16BP; residual = ATP - F16BP) :: EqualAI
            E(; residual = ATP - F16BP) + F6P ⇌ E(F6P; residual = ATP - F16BP)     :: OnlyA
            E(F6P; residual = ATP - F16BP) ⇌ E(ADP)                                :: OnlyA
            E + ADP ⇌ E(ADP)                                                       :: EqualAI
        end
    end
    am = EnzymeRates.AllostericMechanism(dead)
    kids = EnzymeRates._expand_change_allo_state(am)
    isochem(k) = [g for g in eachindex(EnzymeRates.steps(k))
                  if EnzymeRates.is_iso(EnzymeRates.steps(k)[g][1])]
    hasonlyabind(k) = any(EnzymeRates.cat_allo_states(k)[g] === :OnlyA &&
                          EnzymeRates.is_binding(EnzymeRates.steps(k)[g][1])
                          for g in eachindex(EnzymeRates.steps(k)))
    for k in kids
        @test !(hasonlyabind(k) &&
                !all(EnzymeRates.cat_allo_states(k)[g] === :OnlyA for g in isochem(k)))
    end
    @test !isempty(kids)   # the :OnlyA-binding relaxation is retained (non-vacuous)
end

@testset "change_allo_state drops mixed-iso V-type relaxations (no :OnlyA binding)" begin
    # A dead-inactive V-type shape: both chemical (iso) steps :OnlyA, all
    # bindings :EqualAI, no :OnlyA binding. Relaxing ONE iso step to
    # :NonequalAI leaves a partial-catalysis inactive conformation (one iso
    # live, one dead) with NO :OnlyA binding — the mixed-iso kinetic sink. The
    # filter must drop it even though there is no :OnlyA binding to key on.
    vtype = @allosteric_mechanism begin
        substrates: ATP, F6P ; products: ADP, F16BP ; catalytic_multiplicity: 1
        catalytic_steps: begin
            E + ATP ⇌ E(ATP)                                                       :: EqualAI
            E(ATP) <--> E(F16BP; residual = ATP - F16BP)                           :: OnlyA
            E(; residual = ATP - F16BP) + F16BP ⇌ E(F16BP; residual = ATP - F16BP) :: EqualAI
            E(; residual = ATP - F16BP) + F6P ⇌ E(F6P; residual = ATP - F16BP)     :: EqualAI
            E(F6P; residual = ATP - F16BP) ⇌ E(ADP)                                :: OnlyA
            E + ADP ⇌ E(ADP)                                                       :: EqualAI
        end
    end
    am = EnzymeRates.AllostericMechanism(vtype)
    isochem(k) = [g for g in eachindex(EnzymeRates.steps(k))
                  if EnzymeRates.is_iso(EnzymeRates.steps(k)[g][1])]
    kids = EnzymeRates._expand_change_allo_state(am)
    for k in kids
        tags = EnzymeRates.cat_allo_states(k)
        onlya_iso = any(tags[g] === :OnlyA for g in isochem(k))
        live_iso = any(tags[g] !== :OnlyA for g in isochem(k))
        @test !(onlya_iso && live_iso)   # never a mixed-iso partial
        # end-to-end: every surviving child's saturating turnover is finite —
        # no relaxation reintroduces the covalent sink that crashes kcat.
        cm = EnzymeRates.compile_mechanism(k)
        fp = EnzymeRates.fitted_params(cm)
        prm = NamedTuple{(fp..., :Keq, :E_total)}(((1.3 for _ in fp)..., 3.0, 1.0))
        @test isfinite(EnzymeRates._kcat_forward(cm, prm))
    end
    # The chemical steps relax together, so the fully-productive inactive form
    # (every chemical step :NonequalAI) is reachable in one move — even though
    # the one-iso-at-a-time intermediate is a filtered partial.
    @test any(all(EnzymeRates.cat_allo_states(k)[g] === :NonequalAI
                  for g in isochem(k)) for k in kids)
end
