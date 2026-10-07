# ABOUTME: Tests for mechanism enumeration pipeline
# ABOUTME: Unit tests per move + integration tests per reaction

# Fixture rules for this file are in CLAUDE.md, "Enumeration-engine tests":
# every mechanism inline via the macro, every move test asserts the exact
# child set.

# Orients an isomerization Step substrate-rich side → product-rich side (ties
# by lex on form name) and returns its (source, destination) Species. A raw topology
# step keeps the direction it was built in; the mechanism constructors orient
# isomerizations (`_canonicalize_step_directions`).
function _iso_orient(s::EnzymeRates.Step)
    from, to = EnzymeRates.from_species(s), EnzymeRates.to_species(s)
    nf = count(b -> b isa EnzymeRates.Substrate, EnzymeRates.bound(from))
    nt = count(b -> b isa EnzymeRates.Substrate, EnzymeRates.bound(to))
    forward = nf > nt || (nf == nt &&
        string(EnzymeRates.name(from)) <= string(EnzymeRates.name(to)))
    forward ? (from, to) : (to, from)
end

# Build a Mechanism from a flat topology Step list, each step its own
# kinetic group.
_topo_mech(rxn, t::Vector{EnzymeRates.Step}) =
    EnzymeRates.Mechanism(rxn, [[s] for s in t])

# Form-name set for a flat topology Step list, derived from the decomposed
# Species of each Step.
_form_names(t::Vector{EnzymeRates.Step}) = Set{Symbol}(
    EnzymeRates.name(sp)
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

"Finite-difference rank of ∂v/∂log θ over the fitted parameters (test oracle only)."
function _testhelper_identifiable_rank(m; npts = 60, ndraws = 3, h = 1e-5)
    em = EnzymeRates.compile_mechanism(m)
    fp = collect(EnzymeRates.fitted_params(em))
    cm = m isa EnzymeRates.Mechanism ? m : EnzymeRates._state_mechanism(m, :A)
    mets = sort!(collect(EnzymeRates._concentration_symbols(cm)))
    # A fixed seed reproduces a failure in a later session; a struct hash mixes in
    # objectid and would not.
    rng = MersenneTwister(1)
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

"""Numeric degeneracy probe (test oracle only): true when the rate at zero products does not
need some substrate (ratio of the rate at 10⁻⁸ of it to the rate at 1 stays above 10⁻³), or
grows without bound as the substrates scale (10⁹ against 10⁶), or the mirror for
products."""
function _testhelper_degenerate(m; ndraws = 2)
    em = EnzymeRates.compile_mechanism(m)
    fp = collect(EnzymeRates.fitted_params(em))
    cm = m isa EnzymeRates.Mechanism ? m : EnzymeRates._state_mechanism(m, :A)
    subs = [EnzymeRates.name(s) for s in EnzymeRates.substrates(EnzymeRates.reaction(cm))]
    prods = [EnzymeRates.name(p) for p in EnzymeRates.products(EnzymeRates.reaction(cm))]
    mets = sort!(collect(EnzymeRates._concentration_symbols(cm)))
    rng = MersenneTwister(7)
    for _ in 1:ndraws
        θ = exp.(randn(rng, length(fp)))
        p = NamedTuple{(fp..., :Keq, :E_total)}((θ..., exp(randn(rng)), 1.0))
        v(d) = rate_equation(
            em, NamedTuple{Tuple(mets)}(Tuple(get(d, x, 0.0) for x in mets)), p)
        for side in (subs, prods)
            base = v(Dict(x => 1.0 for x in side))
            for x in side
                r = v(Dict(y => (y == x ? 1e-8 : 1.0) for y in side)) / base
                (isfinite(r) && abs(r) < 1e-3) || return true
            end
            big = v(Dict(x => 1e9 for x in side)) / v(Dict(x => 1e6 for x in side))
            big > 10 && return true
        end
    end
    false
end

# Ground truth for the parameter counts the enumeration moves report: the number of
# fitted parameters of the compiled rate equation.
_testhelper_fitted(x) = length(EnzymeRates.fitted_params(x))

# Sorted fitted-parameter count of each child minus its parent's.
function _testhelper_param_deltas(parent, kids)
    base = _testhelper_fitted(parent)
    sort([_testhelper_fitted(k) - base for k in kids])
end

# Whether `m` holds an isomerization step.
_testhelper_holds_iso(m) = any(EnzymeRates.is_iso, Iterators.flatten(EnzymeRates.steps(m)))

# The kinetic groups of `m` whose first step is an isomerization.
_testhelper_iso_groups(m) = [g for g in eachindex(EnzymeRates.steps(m))
                             if EnzymeRates.is_iso(EnzymeRates.steps(m)[g][1])]

# `em` rebuilt on the declared reaction `rxn`. A move's child, a seed and a merged variant
# carry the declared reaction (its atoms and inhibitors), so every fixture is rebuilt on it
# for `==`.
_testhelper_on_reaction(rxn, em) =
    EnzymeRates.Mechanism(rxn, EnzymeRates.steps(EnzymeRates.Mechanism(em)))

# The allosteric mechanism `em` rebuilt on `rxn` at catalytic multiplicity 2 with no
# regulatory sites.
function _testhelper_on_reaction_mult2(rxn, em)
    am = EnzymeRates.AllostericMechanism(em)
    EnzymeRates.AllostericMechanism(rxn, EnzymeRates.steps(am),
        EnzymeRates.cat_allo_states(am), 2, EnzymeRates.RegulatorySite[])
end

# Whether a step binds the metabolite named `x`.
_testhelper_binds(x) = s -> EnzymeRates.bound_metabolite(s) !== nothing &&
                            EnzymeRates.name(EnzymeRates.bound_metabolite(s)) == x

# The kinetic-group index of every step of `m` that binds the metabolite named `x`.
_testhelper_groups_binding(m, x) = [gi for (gi, group) in enumerate(EnzymeRates.steps(m))
                                    for s in group if _testhelper_binds(x)(s)]

# Whether a step is a Theorell–Chance step: it takes up one metabolite and gives off
# another.
_testhelper_is_tc(s) =
    !isempty(EnzymeRates.consumed(s)) && !isempty(EnzymeRates.released(s))

# `gs` with every group holding a step that satisfies one of `preds` turned steady state.
_testhelper_flip_matching(gs, preds) = [
    any(s -> any(p -> p(s), preds), g) ? EnzymeRates._with_equilibrium.(g, false) : g
    for g in gs]

# The children of the four moves `expand_mechanisms` runs on a plain `Mechanism` `m`, in
# its order, before its regulator-type filter.
_testhelper_plain_moves(m, rxn) = append!(
    Union{EnzymeRates.Mechanism, EnzymeRates.AllostericMechanism}[],
    EnzymeRates._expand_re_to_ss(m), EnzymeRates._expand_split_kinetic_group(m),
    EnzymeRates._expand_add_dead_end_regulator(m, rxn),
    EnzymeRates._expand_to_allosteric(m, rxn))

# A kinetic group is catalytic iff it holds a chemistry step (`_is_chemistry`:
# an isomerization, a fused binding or a Theorell–Chance step); otherwise it is
# a binding group. Mirrors the rule `_expand_to_allosteric` uses to decide
# whether a group's `:OnlyA` flip needs a paired regulator to be distinguishable.
_is_catalytic_group(m, g) = any(EnzymeRates._is_chemistry, EnzymeRates.steps(m)[g])

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
    actual(m) = _testhelper_fitted(m)
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
    # One kinetic group holding two bindings of the same metabolite S (context-
    # shared: S binds free E, and S binds A-bound E, with the same K) + a
    # second group with an iso step. Built two ways: reversed outer group
    # order AND swapped inner step order. Canonical-by-construction must
    # collapse them to one struct.
    bind_S = EnzymeRates.Step(
        EnzymeRates.Species(EnzymeRates.Metabolite[], :E),
        EnzymeRates.Species([EnzymeRates.Substrate(:S)], :E_S),
        [EnzymeRates.Substrate(:S)], EnzymeRates.Metabolite[], true)
    bind_S_on_A = EnzymeRates.Step(
        EnzymeRates.Species([EnzymeRates.Substrate(:A)], :E_A),
        EnzymeRates.Species([EnzymeRates.Substrate(:A), EnzymeRates.Substrate(:S)], :E_A_S),
        [EnzymeRates.Substrate(:S)], EnzymeRates.Metabolite[], true)
    iso = EnzymeRates.Step(
        EnzymeRates.Species([EnzymeRates.Substrate(:S)], :E_S),
        EnzymeRates.Species([EnzymeRates.Product(:P)], :E_P),
        EnzymeRates.Metabolite[], EnzymeRates.Metabolite[], false)

    m_orderA = EnzymeRates.Mechanism(rxn, [[bind_S, bind_S_on_A], [iso]])
    m_orderB = EnzymeRates.Mechanism(rxn, [[iso], [bind_S_on_A, bind_S]])
    @test m_orderA == m_orderB
    @test hash(m_orderA) == hash(m_orderB)

    # unique! collapses the duplicate orderings; the mechanisms are not mutated.
    # m_split groups the two binding steps separately, so it is structurally
    # distinct and survives alongside the collapsed orderA/orderB.
    m_split = EnzymeRates.Mechanism(rxn, [[bind_S], [bind_S_on_A], [iso]])
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

@testset "_testhelper_assert_mechanism_invariants: coverage of every reactant" begin
    # POSITIVE: an init mechanism with an unbound declared inhibitor must NOT
    # error — regulators are intentionally excluded from the coverage check
    # (init_mechanisms declares dead-end inhibitors that no step binds yet).
    rxn_inh = @enzyme_reaction begin
        substrates: S[C]
        products:   P[C]
        competitive_inhibitors: R
    end
    for m in EnzymeRates.init_mechanisms(rxn_inh)
        @test _testhelper_assert_mechanism_invariants(m) === nothing
    end

    # NEGATIVE 1: a declared SUBSTRATE that no step binds → error. This half checks
    # that the oracle can fail.
    rxn_unused = @enzyme_reaction begin
        substrates: S[C], T[C]
        products:   P[C2]
    end
    s1 = EnzymeRates.Step(EnzymeRates.Species(EnzymeRates.Metabolite[], :E),
                          EnzymeRates.Species([EnzymeRates.Substrate(:S)], :E_S),
                          [EnzymeRates.Substrate(:S)], EnzymeRates.Metabolite[], true)
    s2 = EnzymeRates.Step(EnzymeRates.Species([EnzymeRates.Substrate(:S)], :E_S),
                          EnzymeRates.Species([EnzymeRates.Product(:P)], :E_P),
                          EnzymeRates.Metabolite[], EnzymeRates.Metabolite[], false)
    s3 = EnzymeRates.Step(EnzymeRates.Species([EnzymeRates.Product(:P)], :E_P),
                          EnzymeRates.Species(EnzymeRates.Metabolite[], :E),
                          EnzymeRates.Metabolite[], [EnzymeRates.Product(:P)], true)
    m_unused = EnzymeRates.Mechanism(rxn_unused, [[s1], [s2], [s3]])
    @test_throws ErrorException _testhelper_assert_mechanism_invariants(m_unused)
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
        EnzymeRates.Metabolite[], EnzymeRates.Metabolite[], false)
    bad_m = EnzymeRates.Mechanism(ldh_rxn, [[bad_iso]])
    @test_throws ErrorException EnzymeRates._assert_atom_conserving(bad_m)
    # The label names the difference the message prints: Lac has 18 fewer C than NADH.
    label = "atoms(to) + atoms(released) − atoms(from) − atoms(consumed) = "
    @test_throws label EnzymeRates._assert_atom_conserving(bad_m)
    @test_throws ":C => -18" EnzymeRates._assert_atom_conserving(bad_m)

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
        [EnzymeRates.Step(E, E_A, [A], EnzymeRates.Metabolite[], true)],
        [EnzymeRates.Step(E_A, E_P_res, EnzymeRates.Metabolite[],
                          EnzymeRates.Metabolite[], false)],
        [EnzymeRates.Step(E_P_res, F, EnzymeRates.Metabolite[], [P], true)],
        [EnzymeRates.Step(F, E_B_res, [B], EnzymeRates.Metabolite[], true)],
        [EnzymeRates.Step(E_B_res, E_Q, EnzymeRates.Metabolite[],
                          EnzymeRates.Metabolite[], true)],
        [EnzymeRates.Step(E_Q, E, EnzymeRates.Metabolite[], [Q], true)],
    ]
    good_m = EnzymeRates.Mechanism(bi_bi_pp_rxn, good_steps)
    @test EnzymeRates._assert_atom_conserving(good_m) === nothing
end


# ═══════════════════════════════════════════════════════════════════════
# 1. Support functions
# ═══════════════════════════════════════════════════════════════════════

# ─── _catalytic_topologies ──────────────────────────────────────────────
@testset "_catalytic_topologies" begin

ter_ter_topos = EnzymeRates._catalytic_topologies(ter_ter_rxn)
uni_bi_topos = EnzymeRates._catalytic_topologies(uni_bi_rxn)

@testset "topology counts" begin
    for (topos, n) in [
        (EnzymeRates._catalytic_topologies(uni_uni_rxn), 1),
        (uni_bi_topos, 3),
        # 9 sequential topologies. A[C]→P[C] and B[N]→Q[N] each leave no
        # covalent residue, so the only ping-pong is degenerate
        # (empty-residue) and is rejected by the admissible-residual rule.
        (EnzymeRates._catalytic_topologies(bi_bi_rxn), 9),
        (EnzymeRates._catalytic_topologies(bi_bi_pp_rxn), 10),
        (ter_ter_topos, 223),
        # 45 = 39 sequential + 6 genuine-residual ping-pong. P[CN] combines
        # A+B, so ping-pong covalent residuals persist on :E; the degenerate
        # empty-residue variants are rejected by the admissible-residual rule.
        (EnzymeRates._catalytic_topologies(ter_bi_rxn), 45),
    ]
        @test length(topos) == n
        @test all(isempty(_connectivity_violations(t)) for t in topos)
        # Every topology has exactly one SS step.
        @test all(count(!s.is_equilibrium for s in t) == 1 for t in topos)
    end
    @test all(EnzymeMechanism(_topo_mech(uni_bi_rxn, t)) isa EnzymeMechanism
              for t in uni_bi_topos)
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
    @test all(EnzymeRates.conformation(sp) === :E
              for t in topos for s in t
              for sp in (EnzymeRates.from_species(s), EnzymeRates.to_species(s)))
    # Surviving ping-pong topologies each carry a genuine (non-empty)
    # covalent intermediate — no empty-residue ping-pong remains.
    pingpong = filter(t -> count(EnzymeRates.is_iso, t) >= 2, topos)
    @test !isempty(pingpong)
    @test all(any(EnzymeRates.has_residual(sp)
                  for s in t
                  for sp in (EnzymeRates.from_species(s), EnzymeRates.to_species(s)))
              for t in pingpong)
end

@testset "isomerization constraints" begin
    topos = ter_ter_topos
    sub_names_set = Set([:A, :B, :D])
    forms(s) = (EnzymeRates.from_species(s), EnzymeRates.to_species(s))

    # C5: at most max(n_subs, n_prods) = 3 metabolites bound on any form.
    @test all(length(EnzymeRates.bound(sp)) <= 3
              for t in topos for s in t for sp in forms(s))

    # C7: every iso source must contain at
    # least one substrate. Use bound list directly (name-parsing is
    # ambiguous with the concat form naming convention).
    @test all(any(b -> EnzymeRates.name(b) ∈ sub_names_set,
                  EnzymeRates.bound(_iso_orient(s)[1]))
              for t in topos for s in t if EnzymeRates.is_iso(s))

    # C8: an iso's product-side (destination) form is built with only
    # products bound, never substrates.
    @test all(EnzymeRates.name(b) ∉ sub_names_set
              for t in topos for s in t if EnzymeRates.is_iso(s)
              for b in EnzymeRates.bound(_iso_orient(s)[2]))
end

_testhelper_bound_names(sp) = Set(EnzymeRates.name(b) for b in EnzymeRates.bound(sp))
# An enzyme form binding exactly the metabolites `names`, with or without a
# covalent residual.
_testhelper_form(names, residual) = sp ->
    _testhelper_bound_names(sp) == Set(names) && EnzymeRates.has_residual(sp) == residual
# Whether `spec` holds an isomerization with one end satisfying `p` and the
# other `q`, whichever way the step is oriented.
_testhelper_iso_between(spec, p, q) = any(spec) do s
    EnzymeRates.is_iso(s) || return false
    f, t = EnzymeRates.from_species(s), EnzymeRates.to_species(s)
    (p(f) && q(t)) || (p(t) && q(f))
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
    @test any(topos) do spec
        _testhelper_iso_between(spec, _testhelper_form([:ATP, :HCO3], false),
                                EnzymeRates.has_residual) &&
        _testhelper_iso_between(spec, _testhelper_form([:OAA], false),
            t -> :Pyr in _testhelper_bound_names(t) && EnzymeRates.has_residual(t))
    end

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
    @test all(count(EnzymeRates.is_iso, t) >= 2 for t in topos)
end

@testset "no admissible catalytic cycle throws" begin
    # Uni-quad: one isomerization would have to give off all four products (C6), and
    # one that gives off fewer leaves a covalent residue after the last substrate has
    # bound.
    uni_quad_rxn = @enzyme_reaction begin
        substrates: S[CNOX]
        products: P[C], Q[N], R[O], T[X]
    end
    @test_throws "no catalytic cycle" EnzymeRates.init_mechanisms(uni_quad_rxn)
    # The product A shares the substrate A's name; reactant atoms are matched by name,
    # so no isomerization balances.
    shared_name_rxn = @enzyme_reaction begin
        substrates: A[C], B[N]
        products: A[CN]
    end
    causes = ["no catalytic cycle", "a product shares a substrate's name"]
    @test_throws causes EnzymeRates.init_mechanisms(shared_name_rxn)
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

# Invalid: {A↔P, B↔P} leaves Q uncovered
@test Set([(:A, :P), (:B, :P)]) ∉ EnzymeRates._competition_patterns(
    Set([:A, :B]), Set([:P, :Q]))

# Every pattern covers all vertices
covers(pat, S, P) = all(any(p -> (s, p) in pat, P) for s in S) &&
                    all(any(s -> (s, p) in pat, S) for p in P)
for (S, P, n) in [
    ([:A, :B], [:P, :Q], 7),              # bi-bi
    ([:A, :B, :C], [:P, :Q, :R], 265),    # ter-ter
    # 2 substrates × 3 products. Bipartite-cover count on K(2,3)
    # by inclusion-exclusion = 25; by symmetry 3 × 2 is also 25.
    ([:A, :B], [:P, :Q, :R], 25),
    ([:A, :B, :D], [:P, :Q], 25),
]
    pats = EnzymeRates._competition_patterns(Set(S), Set(P))
    @test length(pats) == n
    @test all(pat -> covers(pat, S, P), pats)
end
end

# ─── _inhibitor_competition_patterns ────────────────────────────────────
@testset "_inhibitor_competition_patterns" begin
# Uni-uni, no existing inhibitors
pats = EnzymeRates._inhibitor_competition_patterns(
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
    pats = EnzymeRates._inhibitor_competition_patterns(Set(S), Set(P), existing)
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
        E + P ⇌ E(P)
        E + S ⇌ E(S)
        E(S) <--> E(P)
    end
end
mech_uu = EnzymeRates.Mechanism(m_uu)
@test EnzymeRates._forms_where_free(
    mech_uu, EnzymeRates.Reactant, :S) == Set([:E])
@test EnzymeRates._forms_where_free(
    mech_uu, EnzymeRates.Reactant, :P) == Set([:E])

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
@test EnzymeRates._forms_where_free(
    mech_bb, EnzymeRates.Reactant, :B) == Set([:E, :EA])
@test EnzymeRates._forms_where_free(
    mech_bb, EnzymeRates.Reactant, :A) == Set([:E, :EB])
@test EnzymeRates._forms_where_free(
    mech_bb, EnzymeRates.Reactant, :P) == Set([:E, :EQ])
@test EnzymeRates._forms_where_free(
    mech_bb, EnzymeRates.Reactant, :Q) == Set([:E, :EP])
end

@testset "a copy binds at its own sites, apart from its reactant's" begin
    # Ordered bi-bi with A as its own inhibitor at {E, E(Q)}: A binds productively
    # at E, its copy at E and E(Q); Q binds at E and at the copy's form E(A*).
    m = EnzymeRates.Mechanism(@enzyme_mechanism begin
        substrates: A, B
        products: P, Q
        steps: begin
            E + A ⇌ E(A)
            E(A) + B ⇌ E(A, B)
            E(A, B) <--> E(P, Q)
            E(Q) + P ⇌ E(P, Q)
            (E + Q ⇌ E(Q), E(A::Inh) + Q ⇌ E(A::Inh, Q))
            (E + A::Inh ⇌ E(A::Inh), E(Q) + A::Inh ⇌ E(A::Inh, Q))
        end
    end)
    @test EnzymeRates._forms_where_free(
        m, EnzymeRates.Reactant, :A) == Set([:E])
    @test EnzymeRates._forms_where_free(
        m, EnzymeRates.CompetitiveInhibitor, :A) == Set([:E, :EQ])
    @test EnzymeRates._forms_where_free(
        m, EnzymeRates.Reactant, :Q) == Set([:E, :EAinh])
    @test isempty(EnzymeRates._forms_where_free(
        m, EnzymeRates.CompetitiveInhibitor, :Q))
end

@testset "a Theorell–Chance step frees B at E(A) and P at E(Q)" begin
    # E(A) + B → E(Q) + P takes B up at E(A) and gives P off at E(Q); A and Q
    # bind at E.
    m = EnzymeRates.Mechanism(@enzyme_mechanism begin
        substrates: A, B
        products: P, Q
        steps: begin
            E + A <--> E(A)
            E(A) + B <--> E(Q) + P
            E + Q <--> E(Q)
        end
    end)
    @test EnzymeRates._forms_where_free(m, EnzymeRates.Reactant, :A) == Set([:E])
    @test EnzymeRates._forms_where_free(m, EnzymeRates.Reactant, :B) == Set([:EA])
    @test EnzymeRates._forms_where_free(m, EnzymeRates.Reactant, :P) == Set([:EQ])
    @test EnzymeRates._forms_where_free(m, EnzymeRates.Reactant, :Q) == Set([:E])
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
              for steps in result)
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
              for steps in result)
    # 4 unique dead-end forms, 7 competition patterns,
    # all 7 produce distinct dead-end sets → 7 variants
    @test length(result) == 7

    # Each variant adds exactly the dead-end forms whose (substrate, product)
    # pair its competition pattern leaves unforbidden.
    seed_forms = _form_names(topo)
    @test Set(setdiff(_form_names(r), seed_forms) for r in result) == Set([
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
              for steps in result)
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
            E + A ⇌ E(A)
            Estar + B ⇌ Estar(B)
            E + Q ⇌ E(Q)
            Estar + P ⇌ Estar(P)
            E(A) <--> Estar(P)
            Estar(B) ⇌ E(Q)
        end
    end
    topo = _flat_topo(m)
    result =
        EnzymeRates._expand_substrate_product_dead_ends(
            [topo],bi_bi_pp_rxn)
    @test all(isempty(_connectivity_violations(steps))
              for steps in result)
    # 6 dead-end forms (E_A_P, E_A_Q, E_B_Q from
    # E-side + Estar_A_P, Estar_B_P, Estar_B_Q from
    # Estar-side), competition-filtered
    @test length(result) == 7

    # Each variant adds exactly the dead-end forms whose (substrate, product)
    # pair its competition pattern leaves unforbidden.
    seed_forms = _form_names(topo)
    @test Set(setdiff(_form_names(r), seed_forms) for r in result) == Set([
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
        new_forms = setdiff(_form_names(r), seed_forms)
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
    topos = EnzymeRates._catalytic_topologies(
        ter_ter_rxn)
    @test length(topos) == 223
    # Test first (random, most forms) and last topology
    for topo in [topos[1], topos[end]]
        result =
            EnzymeRates._expand_substrate_product_dead_ends(
                [topo], ter_ter_rxn)
        @test all(isempty(_connectivity_violations(steps))
                  for steps in result)
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
    topos = EnzymeRates._catalytic_topologies(ter_ter_rxn)
    _, idx = findmax(length(_form_names(t)) for t in topos)
    result = EnzymeRates._expand_substrate_product_dead_ends([topos[idx]], ter_ter_rxn)
    added = [setdiff(_form_names(steps), _form_names(topos[idx])) for steps in result]
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
end

# ═══════════════════════════════════════════════════════════════════════
# Testsets covering non-enumeration features (AllostericEnzymeMechanism
# accessor identity)
# ═══════════════════════════════════════════════════════════════════════

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
         if pred(EnzymeRates.bound_metabolite(first(EnzymeRates.steps(am)[g])))))
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

# The init mechanisms of uni_uni_rxn, uni_bi_rxn, bi_bi_rxn and bi_bi_pp_rxn.
uu, ub, bb, pp = (EnzymeRates.init_mechanisms(rxn)
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
        @test all(isempty(_connectivity_violations(
            EnzymeRates.steps(m))) for m in specs)
        for s in specs
            if _testhelper_holds_iso(s)
                @test count(st -> !st.is_equilibrium,
                            Iterators.flatten(s.steps)) == 1
            else
                @test count(g -> !EnzymeRates.is_equilibrium(first(g)),
                            s.steps) in (2, 3)
            end
        end
    end
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

@testset "init_mechanisms: the seeds, then their merged and Theorell–Chance variants" begin
    ER = EnzymeRates
    function kind(m)
        st = collect(Iterators.flatten(ER.steps(m)))
        tc = any(s -> !isempty(ER.consumed(s)) && !isempty(ER.released(s)), st)
        fused = any(s -> ER.is_binding(s) && ER._is_chemistry(s), st)
        (tc ? (fused ? :half_tc : :tc) : :merged,
         count(g -> !ER.is_equilibrium(first(g)), ER.steps(m)))
    end
    tally(ms) = Dict(k => count(m -> kind(m) == k, ms) for k in unique(kind.(ms)))
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
    @test tally(pp[63:end]) == Dict((:merged, 2) => 122, (:merged, 3) => 56,
                                    (:tc, 3) => 20, (:half_tc, 3) => 4)
    @test count(ER._degenerate, pp[1:62]) == 7 && !any(ER._degenerate, pp[63:end])
    @test all(m -> ER._assert_emission_rules(m) === nothing, pp)
    # bi_bi_rxn's atoms admit no ping-pong: the 55 sequential seeds and their 184 variants.
    @test length(bb) == 239 && allunique(bb)
    @test all(_testhelper_holds_iso, bb[1:55]) && !any(_testhelper_holds_iso, bb[56:end])
    @test tally(bb[56:end]) == Dict((:merged, 2) => 108, (:merged, 3) => 56, (:tc, 3) => 20)
end

@testset "init_mechanisms on uni-bi: one merged variant per ordered seed" begin
    ER = EnzymeRates
    # Aggregate pin over the uni-bi seed set: the two ordered seeds (P or Q released first)
    # and the random-release seed. An ordered seed's merged base holds S (E + S → E(P, Q)),
    # the inner product's group (into E(P, Q)) and the outer product's group. Singletons
    # fail V; {S, inner} is the lumping twin; {inner, outer} has no reverse maximal rate;
    # {S, outer} is the one variant (n + m − 2 = 1). Its Theorell–Chance base gives none:
    # each single group leaves the TC step RE or closes the cycle, and {TC, outer} turned
    # RE for products closes it as well. The random-release base holds S, Pˣ and Qˣ, each
    # product group with a step into E(P, Q); {S, Pˣ} and {S, Qˣ} pass and are no twins
    # (E(P, Q) has three steps), {Pˣ, Qˣ} has no reverse maximal rate, and ELIM refuses the
    # three-step complex. Each variant fits 2 + 1 + 2 − 1 = 4 (2 + 2 + 1 − 1 = 4 for the
    # random seed's), the seeds' count.
    ordered_p = _testhelper_on_reaction(uni_bi_rxn, @enzyme_mechanism begin
        substrates: S; products: P, Q
        steps: begin
            E + S <--> E(P, Q)
            E(Q) + P ⇌ E(P, Q)
            E + Q <--> E(Q)
        end
    end)
    ordered_q = _testhelper_on_reaction(uni_bi_rxn, @enzyme_mechanism begin
        substrates: S; products: P, Q
        steps: begin
            E + S <--> E(P, Q)
            E(P) + Q ⇌ E(P, Q)
            E + P <--> E(P)
        end
    end)
    random_p = _testhelper_on_reaction(uni_bi_rxn, @enzyme_mechanism begin
        substrates: S; products: P, Q
        steps: begin
            E + S <--> E(P, Q)
            (E(Q) + P <--> E(P, Q), E + P <--> E(P))
            (E(P) + Q ⇌ E(P, Q), E + Q ⇌ E(Q))
        end
    end)
    random_q = _testhelper_on_reaction(uni_bi_rxn, @enzyme_mechanism begin
        substrates: S; products: P, Q
        steps: begin
            E + S <--> E(P, Q)
            (E(Q) + P ⇌ E(P, Q), E + P ⇌ E(P))
            (E(P) + Q <--> E(P, Q), E + Q <--> E(Q))
        end
    end)
    @test length(ub) == 7
    @test all(_testhelper_holds_iso, ub[1:3])
    @test Set(ub[4:7]) == Set([ordered_p, ordered_q, random_p, random_q])
    for seed in ub[1:3]
        @test length(ER._seed_variants(seed)) ==
              (sum(length, ER.steps(seed)) == 4 ? 1 : 2)
    end
    for v in ub[4:7]
        fitted = _testhelper_fitted(v)
        @test (fitted, _testhelper_identifiable_rank(v)) == (4, 4)
    end
end

@testset "init_mechanisms on ter-ter within 150 s" begin
    # Aggregate pin over the whole ter-ter seed set: 35,665 seeds and their 215,190 merged
    # and Theorell–Chance variants. Enumeration, not compilation, dominates the call: it
    # takes ~50 s locally and 70-115 s on CI runners. The budget is 2× the Linux runner.
    t = @elapsed ms = EnzymeRates.init_mechanisms(ter_ter_rxn)
    @test length(ms) == 250855
    @test t < 150
end

@testset "Drops unbound regulators from init Mechanism" begin
    # init_mechanisms produces Mechanisms without dead-end regulators
    # bound. When compiled to EnzymeMechanism, the regulator must NOT
    # appear among the regulators of the lifted reaction
    # (`_testhelper_regulators`) — only the catalytic mechanism is
    # built. After expand_mechanisms adds the dead-end regulator, it
    # should appear.
    init_mechs = EnzymeRates.init_mechanisms(uni_uni_with_reg)
    @test all(isempty(_connectivity_violations(
        EnzymeRates.steps(m))) for m in init_mechs)
    @test !isempty(init_mechs)
    for m in init_mechs
        em = EnzymeRates.compile_mechanism(m)
        @test :I ∉ _testhelper_regulators(em)
    end

    expanded = EnzymeRates.expand_mechanisms(init_mechs, uni_uni_with_reg)
    found_with_reg = false
    for mm in expanded
        em = EnzymeRates.compile_mechanism(mm)
        if :I in _testhelper_regulators(em)
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
    # The two pairs at E(P) and at E(Q) isolate one form joined to the rest by two
    # edges of stored weights −1 and +1, whose 2-cycle sums to 0 (a balanced block):
    # neither flipped group carries flux, so they fail the same way and keep their
    # parent's rank. Children: the S flip plus the three other pairs.
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
    @test length(children) == 4
    @test Set(children) == Set(expected)
    # The two zero-flux pairs are absent, and each has exactly its parent's rank.
    r0 = _testhelper_identifiable_rank(parent)
    for absent in EnzymeRates.Mechanism.([
        (@enzyme_mechanism begin      # P@E and Q@E(P)
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
        (@enzyme_mechanism begin      # Q@E and P@E(Q)
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
        end)
    ])
        @test !(absent in children)
        @test EnzymeRates._re_segment_count(EnzymeRates.steps(absent)) >
              EnzymeRates._re_segment_count(EnzymeRates.steps(parent))
        @test _testhelper_identifiable_rank(absent) == r0
    end
end

@testset "Mechanism — the flip skips the flanks of a qualifying chain" begin
    # Ordered RE seed: E + A ⇌ E(A), E(A) + B ⇌ E(A, B), E(A, B) → E(P, Q) (SS),
    # E(Q) + P ⇌ E(P, Q), E + Q ⇌ E(Q). E(A, B) and E(P, Q) each have one other step, a
    # binding into them, and every group holds one step: a qualifying chain (RSR). Its
    # flanks, the B and P groups, never flip: a flipped flank gives its parent's family
    # plus a phantom (chain lemma). The A and Q flips each cut a segment and carry flux.
    seed = EnzymeRates.Mechanism(@enzyme_mechanism begin
        substrates: A, B; products: P, Q
        steps: begin
            E + A ⇌ E(A)
            E(A) + B ⇌ E(A, B)
            E(A, B) <--> E(P, Q)
            E(Q) + P ⇌ E(P, Q)
            E + Q ⇌ E(Q)
        end
    end)
    a_flip = EnzymeRates.Mechanism(@enzyme_mechanism begin
        substrates: A, B; products: P, Q
        steps: begin
            E + A <--> E(A)
            E(A) + B ⇌ E(A, B)
            E(A, B) <--> E(P, Q)
            E(Q) + P ⇌ E(P, Q)
            E + Q ⇌ E(Q)
        end
    end)
    q_flip = EnzymeRates.Mechanism(@enzyme_mechanism begin
        substrates: A, B; products: P, Q
        steps: begin
            E + A ⇌ E(A)
            E(A) + B ⇌ E(A, B)
            E(A, B) <--> E(P, Q)
            E(Q) + P ⇌ E(P, Q)
            E + Q <--> E(Q)
        end
    end)
    kids = EnzymeRates._expand_re_to_ss(seed)
    @test length(kids) == 2
    @test Set(kids) == Set([a_flip, q_flip])
    # The seed has 5 fitted constants, all identifiable; each child adds one and
    # gains one. Each flank flip, absent, adds a constant and no rank.
    @test _testhelper_fitted(seed) == _testhelper_identifiable_rank(seed) == 5
    for k in kids
        @test _testhelper_fitted(k) == _testhelper_identifiable_rank(k) == 6
    end
    for absent in EnzymeRates.Mechanism.([
        (@enzyme_mechanism begin      # the B flank
            substrates: A, B; products: P, Q
            steps: begin
                E + A ⇌ E(A)
                E(A) + B <--> E(A, B)
                E(A, B) <--> E(P, Q)
                E(Q) + P ⇌ E(P, Q)
                E + Q ⇌ E(Q)
            end
        end),
        (@enzyme_mechanism begin      # the P flank
            substrates: A, B; products: P, Q
            steps: begin
                E + A ⇌ E(A)
                E(A) + B ⇌ E(A, B)
                E(A, B) <--> E(P, Q)
                E(Q) + P <--> E(P, Q)
                E + Q ⇌ E(Q)
            end
        end)
    ])
        @test !(absent in kids)
        @test _testhelper_fitted(absent) == 6
        @test _testhelper_identifiable_rank(absent) == 5
    end
end

@testset "Mechanism — a flank whose group is shared still flips" begin
    # The ordered seed with B's abortive binding at E(Q) in B's group. One RE segment
    # holds every form; the iso is steady state. E(A, B)'s binding shares its group
    # with E(Q) + B ⇌ E(B, Q), so the chain E(A, B) → E(P, Q) does not qualify and
    # no group is a flank. Each flip cuts the segment alone:
    #   A: {E(A), E(A, B)} | {E, E(Q), E(B, Q), E(P, Q)}             → 2 segments
    #   B: {E(A, B)} | {E(B, Q)} | {E, E(A), E(Q), E(P, Q)}           → 3 segments
    #   P: {E(P, Q)} | {E, E(A), E(A, B), E(Q), E(B, Q)}              → 2 segments
    #   Q: {E, E(A), E(A, B)} | {E(Q), E(B, Q), E(P, Q)}              → 2 segments
    # Every segment keeps a bottom form, and each flipped group holds a step of the
    # cycle E → E(A) → E(A, B) → E(P, Q) → E(Q) → E, so it carries flux: four children,
    # each one constant and one rank above the seed's 5.
    seed = EnzymeRates.Mechanism(@enzyme_mechanism begin
        substrates: A, B; products: P, Q
        steps: begin
            E + A ⇌ E(A)
            (E(A) + B ⇌ E(A, B), E(Q) + B ⇌ E(B, Q))
            E(A, B) <--> E(P, Q)
            E(Q) + P ⇌ E(P, Q)
            E + Q ⇌ E(Q)
        end
    end)
    a_flip = EnzymeRates.Mechanism(@enzyme_mechanism begin
        substrates: A, B; products: P, Q
        steps: begin
            E + A <--> E(A)
            (E(A) + B ⇌ E(A, B), E(Q) + B ⇌ E(B, Q))
            E(A, B) <--> E(P, Q)
            E(Q) + P ⇌ E(P, Q)
            E + Q ⇌ E(Q)
        end
    end)
    b_flip = EnzymeRates.Mechanism(@enzyme_mechanism begin
        substrates: A, B; products: P, Q
        steps: begin
            E + A ⇌ E(A)
            (E(A) + B <--> E(A, B), E(Q) + B <--> E(B, Q))
            E(A, B) <--> E(P, Q)
            E(Q) + P ⇌ E(P, Q)
            E + Q ⇌ E(Q)
        end
    end)
    p_flip = EnzymeRates.Mechanism(@enzyme_mechanism begin
        substrates: A, B; products: P, Q
        steps: begin
            E + A ⇌ E(A)
            (E(A) + B ⇌ E(A, B), E(Q) + B ⇌ E(B, Q))
            E(A, B) <--> E(P, Q)
            E(Q) + P <--> E(P, Q)
            E + Q ⇌ E(Q)
        end
    end)
    q_flip = EnzymeRates.Mechanism(@enzyme_mechanism begin
        substrates: A, B; products: P, Q
        steps: begin
            E + A ⇌ E(A)
            (E(A) + B ⇌ E(A, B), E(Q) + B ⇌ E(B, Q))
            E(A, B) <--> E(P, Q)
            E(Q) + P ⇌ E(P, Q)
            E + Q <--> E(Q)
        end
    end)
    @test isempty(EnzymeRates._chain_flank_groups(seed))
    kids = EnzymeRates._expand_re_to_ss(seed)
    @test length(kids) == 4
    @test Set(kids) == Set([a_flip, b_flip, p_flip, q_flip])
    @test _testhelper_fitted(seed) == _testhelper_identifiable_rank(seed) == 5
    for k in kids
        @test _testhelper_fitted(k) == _testhelper_identifiable_rank(k) == 6
    end
end

@testset "AllostericMechanism — a chain flank still flips at multiplicity 2" begin
    # A K-type variant of the ordered seed over two catalytic subunits: the iso and the
    # A binding :OnlyA, the other bindings :EqualAI. Its flank flips stay children: two
    # conformations over two subunits can see a flank's flip. The A-state graph is the
    # plain seed's, so each single flip cuts the one segment in two, keeps a bottom
    # form and carries flux. The parent needs a hyperbolic scheme
    # (`_hyperbolic_catalysis`). With two segments an arborescence toward a root holds
    # one edge, which leaves the other segment, so the scheme fails only when a
    # metabolite scores 2 on one edge or one segment, or scores in one segment and on
    # an edge leaving the other. Each flip passes:
    #   A: {E, E(Q), E(P, Q)} | {E(A), E(A, B)}. A scores on E → E(A) only; B in the
    #      second segment and on E(A, B) → E(P, Q), which leaves it; P and Q in the
    #      first and on E(P, Q) → E(A, B), which leaves it.
    #   Q: the mirror of A.
    #   B: {E, E(A), E(Q), E(P, Q)} | {E(A, B)}. A, P and Q score in the first segment
    #      and on edges leaving it (E(A) → E(A, B), E(P, Q) → E(A, B)); B scores 1 on
    #      E(A) → E(A, B) and in no segment.
    #   P: the mirror of B.
    # Four children, tags kept. The parent is identifiable at 6; every child adds one
    # constant and gains one, the B and P flank flips included, where the plain
    # seed's flank flips gain none.
    seed = EnzymeRates.Mechanism(@enzyme_mechanism begin
        substrates: A, B; products: P, Q
        steps: begin
            E + A ⇌ E(A)
            E(A) + B ⇌ E(A, B)
            E(A, B) <--> E(P, Q)
            E(Q) + P ⇌ E(P, Q)
            E + Q ⇌ E(Q)
        end
    end)
    rxn2 = @enzyme_reaction begin
        substrates: A[C], B[N]
        products: P[C], Q[N]
        oligomeric_state: 2
    end
    am = EnzymeRates.AllostericMechanism(@allosteric_mechanism begin
        substrates: A, B ; products: P, Q ; catalytic_multiplicity: 2
        catalytic_steps: begin
            E + A ⇌ E(A)                :: OnlyA
            E(A) + B ⇌ E(A, B)          :: EqualAI
            E(A, B) <--> E(P, Q)        :: OnlyA
            E(Q) + P ⇌ E(P, Q)          :: EqualAI
            E + Q ⇌ E(Q)                :: EqualAI
        end
    end)
    @test am in EnzymeRates._expand_to_allosteric(seed, rxn2)
    @test isempty(EnzymeRates._chain_flank_groups(am))
    a_flip = EnzymeRates.AllostericMechanism(@allosteric_mechanism begin
        substrates: A, B ; products: P, Q ; catalytic_multiplicity: 2
        catalytic_steps: begin
            E + A <--> E(A)             :: OnlyA
            E(A) + B ⇌ E(A, B)          :: EqualAI
            E(A, B) <--> E(P, Q)        :: OnlyA
            E(Q) + P ⇌ E(P, Q)          :: EqualAI
            E + Q ⇌ E(Q)                :: EqualAI
        end
    end)
    b_flip = EnzymeRates.AllostericMechanism(@allosteric_mechanism begin
        substrates: A, B ; products: P, Q ; catalytic_multiplicity: 2
        catalytic_steps: begin
            E + A ⇌ E(A)                :: OnlyA
            E(A) + B <--> E(A, B)       :: EqualAI
            E(A, B) <--> E(P, Q)        :: OnlyA
            E(Q) + P ⇌ E(P, Q)          :: EqualAI
            E + Q ⇌ E(Q)                :: EqualAI
        end
    end)
    p_flip = EnzymeRates.AllostericMechanism(@allosteric_mechanism begin
        substrates: A, B ; products: P, Q ; catalytic_multiplicity: 2
        catalytic_steps: begin
            E + A ⇌ E(A)                :: OnlyA
            E(A) + B ⇌ E(A, B)          :: EqualAI
            E(A, B) <--> E(P, Q)        :: OnlyA
            E(Q) + P <--> E(P, Q)       :: EqualAI
            E + Q ⇌ E(Q)                :: EqualAI
        end
    end)
    q_flip = EnzymeRates.AllostericMechanism(@allosteric_mechanism begin
        substrates: A, B ; products: P, Q ; catalytic_multiplicity: 2
        catalytic_steps: begin
            E + A ⇌ E(A)                :: OnlyA
            E(A) + B ⇌ E(A, B)          :: EqualAI
            E(A, B) <--> E(P, Q)        :: OnlyA
            E(Q) + P ⇌ E(P, Q)          :: EqualAI
            E + Q <--> E(Q)             :: EqualAI
        end
    end)
    kids = EnzymeRates._expand_re_to_ss(am)
    @test length(kids) == 4
    @test Set(kids) == Set([a_flip, b_flip, p_flip, q_flip])
    @test all(EnzymeRates._hyperbolic_catalysis, kids)
    @test _testhelper_fitted(am) == _testhelper_identifiable_rank(am) == 6
    for k in kids
        @test _testhelper_fitted(k) == _testhelper_identifiable_rank(k) == 7
    end
end

@testset "_chain_flank_groups: the flanks of steady-state chains only" begin
    # Ping-pong seed. E(A) → Estar(P) is steady state, alone in its group, and its two
    # ends each have one other step, a binding into that end, alone in its group: the
    # A and P groups are its flanks. Estar(B) ⇌ E(Q) is a chain of the same shape
    # whose isomerization is at rapid equilibrium, so the B and Q groups are no
    # flanks.
    m = EnzymeRates.Mechanism(@enzyme_mechanism begin
        substrates: A, B
        products: P, Q
        steps: begin
            E + A ⇌ E(A)
            Estar + B ⇌ Estar(B)
            E + Q ⇌ E(Q)
            Estar + P ⇌ Estar(P)
            E(A) <--> Estar(P)
            Estar(B) ⇌ E(Q)
        end
    end)
    flank_steps = Set(EnzymeRates.steps(m)[g]
                      for g in EnzymeRates._chain_flank_groups(m))
    @test length(flank_steps) == 2
    @test all(g -> length(g) == 1 && EnzymeRates.is_binding(only(g)), flank_steps)
    @test Set(EnzymeRates.name(EnzymeRates.bound_metabolite(only(g)))
              for g in flank_steps) == Set([:A, :P])
end

@testset "_chain_flank_groups: a flank binds into its end of the chain" begin
    # Iso uni-uni with P also binding E(S) abortively. Both isomerizations are steady
    # state and alone in their groups, and no chain qualifies:
    #   E(S) → Estar(P): E(S) has two other steps, the S binding into it and the P
    #     binding out of it (E(S) + P ⇌ E(P, S)), so it is no chain end.
    #   Estar → E: each end has exactly one other step, a binding alone in its group,
    #     but out of the end, not into it (Estar + P ⇌ Estar(P), E + S ⇌ E(S)), so
    #     neither is a flank.
    # Without the orientation condition the second isomerization would name the S and P
    # groups, the flanks the first would name without the abortive step.
    m = EnzymeRates.Mechanism(@enzyme_mechanism begin
        substrates: S
        products: P
        steps: begin
            E + S ⇌ E(S)
            E(S) + P ⇌ E(P, S)
            E(S) <--> Estar(P)
            Estar + P ⇌ Estar(P)
            Estar <--> E
        end
    end)
    @test isempty(EnzymeRates._chain_flank_groups(m))
end

@testset "Mechanism — bi-bi ping-pong: 5 RE groups → 3 variants" begin
    # SEED: bi-bi ping-pong with Estar (residual) form. 5 singleton RE
    # groups + 1 SS iso group. E(A) → Estar(P) is a qualifying chain with a
    # steady-state isomerization: E(A) and Estar(P) each have one other step,
    # a binding into that form alone in its group, so its flanks, the A and P
    # groups, never flip (`_chain_flank_groups`). Estar(B) ⇌ E(Q) has the same
    # shape at rapid equilibrium, so the B and Q groups stay units. The B, Q
    # and RE iso groups each flip alone → 3.
    em_seed = @enzyme_mechanism begin
        substrates: A, B
        products: P, Q
        steps: begin
            E + A ⇌ E(A)
            Estar + B ⇌ Estar(B)
            E + Q ⇌ E(Q)
            Estar + P ⇌ Estar(P)
            E(A) <--> Estar(P)
            Estar(B) ⇌ E(Q)
        end
    end
    m = EnzymeRates.Mechanism(em_seed)
    _testhelper_assert_mechanism_invariants(m)

    result = EnzymeRates._expand_re_to_ss(m)

    # 1. count: 5 all-RE groups, of which A and P are chain flanks → 3 variants.
    @test length(result) == 3
    for r in result
        @test r isa EnzymeRates.Mechanism
        _testhelper_assert_mechanism_invariants(r)
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
    @test length(unique(flipped_groups)) == 3

    # 3. preservation
    for r in result
        @test EnzymeRates.reaction(r) == EnzymeRates.reaction(m)
    end

    # 4. The A and P flank flips, absent, keep the seed's rank (chain lemma).
    r0 = _testhelper_identifiable_rank(m)
    for absent in EnzymeRates.Mechanism.([
        (@enzyme_mechanism begin      # the A flank
            substrates: A, B
            products: P, Q
            steps: begin
                E + A <--> E(A)
                Estar + B ⇌ Estar(B)
                E + Q ⇌ E(Q)
                Estar + P ⇌ Estar(P)
                E(A) <--> Estar(P)
                Estar(B) ⇌ E(Q)
            end
        end),
        (@enzyme_mechanism begin      # the P flank
            substrates: A, B
            products: P, Q
            steps: begin
                E + A ⇌ E(A)
                Estar + B ⇌ Estar(B)
                E + Q ⇌ E(Q)
                Estar + P <--> Estar(P)
                E(A) <--> Estar(P)
                Estar(B) ⇌ E(Q)
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
    # E(A, B, D) → E(P, Q, R) (`_chain_flank_groups`) and never flip; the other
    # four each cut the one segment alone → 4 variants.
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
    _testhelper_assert_mechanism_invariants(m)

    result = EnzymeRates._expand_re_to_ss(m)
    @test length(result) == 4
    for r in result
        @test r isa EnzymeRates.Mechanism
        _testhelper_assert_mechanism_invariants(r)
        @test EnzymeRates.compile_mechanism(r) isa EnzymeMechanism
    end
    # The D and P flank flips, absent, keep the seed's rank (chain lemma).
    r0 = _testhelper_identifiable_rank(m)
    for absent in EnzymeRates.Mechanism.([
        (@enzyme_mechanism begin      # the D flank
            substrates: A, B, D
            products: P, Q, R
            steps: begin
                E + A ⇌ E(A)
                E(A) + B ⇌ E(A, B)
                E(A, B) + D <--> E(A, B, D)
                E + R ⇌ E(R)
                E(R) + Q ⇌ E(Q, R)
                E(Q, R) + P ⇌ E(P, Q, R)
                E(A, B, D) <--> E(P, Q, R)
            end
        end),
        (@enzyme_mechanism begin      # the P flank
            substrates: A, B, D
            products: P, Q, R
            steps: begin
                E + A ⇌ E(A)
                E(A) + B ⇌ E(A, B)
                E(A, B) + D ⇌ E(A, B, D)
                E + R ⇌ E(R)
                E(R) + Q ⇌ E(Q, R)
                E(Q, R) + P <--> E(P, Q, R)
                E(A, B, D) <--> E(P, Q, R)
            end
        end)
    ])
        @test !(absent in result)
        @test _testhelper_identifiable_rank(absent) == r0
    end
end

@testset "AllostericMechanism — uni-uni RE→SS flips one group and preserves the tags" begin
    # Each seed is a uni-uni allosteric mechanism whose two binding groups are all-RE and
    # whose isomerization is SS: _expand_re_to_ss fires per all-RE group → 2 variants.
    for (am, deltas) in [
        # :EqualAI groups. Each converted group stays :EqualAI (RE K → SS (kf, kr),
        # shared R/T). Cheap-tag RE→SS adds 1 fitted param per variant (no separate
        # T-state pair).
        (EnzymeRates.AllostericMechanism(@allosteric_mechanism begin
            substrates: S
            products: P
            catalytic_multiplicity: 2
            catalytic_steps: begin
                E + P ⇌ E(P)      :: EqualAI
                E + S ⇌ E(S)      :: EqualAI
                E(S) <--> E(P)    :: EqualAI
            end
        end), [1, 1]),
        # Both binding groups :OnlyA (a balanced substrate/product pair, so the
        # mechanism is Haldane-valid) and :EqualAI catalysis. :OnlyA groups live in the
        # R-state only; the T-state contributes no kf_T/kr_T after RE→SS. The move MUST
        # NOT change R/T-state semantics, so the converted group stays :OnlyA.
        # Connected-component pruning populates the inactive state from free E over ALL
        # surviving steps, so because catalysis is :EqualAI the base mechanism's Q_I
        # ALREADY carries E(S) and E(P) through the inactive catalytic step — even
        # though both bindings are :OnlyA and cannot bind directly in the T-state.
        # Those forms are therefore present in the inactive state of the base AND of
        # every variant, so flipping a binding group RE→SS adds only that group's own SS
        # rate constant (+1) — it does not introduce a new inactive-state form. Both
        # variants → Δ=1.
        (EnzymeRates.AllostericMechanism(@allosteric_mechanism begin
            substrates: S
            products: P
            catalytic_multiplicity: 2
            catalytic_steps: begin
                E + P ⇌ E(P)      :: OnlyA
                E + S ⇌ E(S)      :: OnlyA
                E(S) <--> E(P)    :: EqualAI
            end
        end), [1, 1]),
        # One :NonequalAI group (S-binding), the others :EqualAI. When RE→SS converts a
        # :NonequalAI group, BOTH the R-state K and the T-state K_T must split into
        # (kf, kr) and (kf_T, kr_T): P-binding :EqualAI → +1; S-binding :NonequalAI → +2.
        # SS :NonequalAI adds kf_T, kr_T on top of the R-state pair.
        (EnzymeRates.AllostericMechanism(@allosteric_mechanism begin
            substrates: S
            products: P
            catalytic_multiplicity: 2
            catalytic_steps: begin
                E + P ⇌ E(P)      :: EqualAI
                E + S ⇌ E(S)      :: NonequalAI
                E(S) <--> E(P)    :: EqualAI
            end
        end), [1, 2]),
        # Every group :NonequalAI — a distinguishable mechanism (an all-:EqualAI seed, the
        # first entry above, is a model-space no-op that is never enumerated; the move is
        # still tested on it). :NonequalAI RE→SS adds 2 fitted params per variant (the A-
        # and I-state binding K each split into (kf, kr)).
        (EnzymeRates.AllostericMechanism(@allosteric_mechanism begin
            substrates: S
            products: P
            catalytic_multiplicity: 2
            catalytic_steps: begin
                E + S ⇌ E(S)    :: NonequalAI
                E + P ⇌ E(P)    :: NonequalAI
                E(S) <--> E(P)  :: NonequalAI
            end
        end), [2, 2]),
    ]
        _testhelper_assert_mechanism_invariants(am)
        result = EnzymeRates._expand_re_to_ss(am)
        @test length(result) == 2
        @test _testhelper_param_deltas(am, result) == deltas
        for r in result
            _testhelper_assert_mechanism_invariants(r)
            # Exactly one group newly all-SS; every tag, the multiplicity, the
            # regulatory sites and the reaction preserved.
            n_newly_ss = count(zip(am.cat_steps, r.cat_steps)) do (old, new)
                all(EnzymeRates.is_equilibrium, old) &&
                    !any(EnzymeRates.is_equilibrium, new)
            end
            @test n_newly_ss == 1
            @test r.cat_allo_states == am.cat_allo_states
            @test r.catalytic_multiplicity == am.catalytic_multiplicity
            @test r.regulatory_sites == am.regulatory_sites
            @test EnzymeRates.reaction(r) == EnzymeRates.reaction(am)
        end
    end
end

@testset "Mechanism — a substrate's inhibitor copy never flips; catalytic groups do" begin
    # Ordered bi-bi with A also bound as a competitive-inhibitor copy at E and
    # E(Q), the Q binding mirrored onto the copy forms in the Q group, as the
    # dead-end move builds it. The copy group binds a regulator and is never a
    # unit. A, B, P and Q each divide the one segment alone (the Q flip takes its
    # mirror with it and isolates {E(Q), E(P, Q), E(A::Inh, Q)}), and each flipped
    # group shares a block with the chemistry edge, so all four carry flux. The copy
    # binds at E and E(Q), so E(A, B) and E(P, Q) keep one step each beside the
    # isomerization, a binding into that form alone in its group: B and P are the
    # flanks of a qualifying chain and never flip (`_chain_flank_groups`). Their
    # flips, absent, keep the parent's rank. Children: the A and Q flips.
    m = EnzymeRates.Mechanism(@enzyme_mechanism begin
        substrates: A, B
        products: P, Q
        steps: begin
            E + A ⇌ E(A)
            E(A) + B ⇌ E(A, B)
            E(A, B) <--> E(P, Q)
            E(Q) + P ⇌ E(P, Q)
            (E + Q ⇌ E(Q), E(A::Inh) + Q ⇌ E(A::Inh, Q))
            (E + A::Inh ⇌ E(A::Inh), E(Q) + A::Inh ⇌ E(A::Inh, Q))
        end
    end)
    expected = [
        EnzymeRates.Mechanism(@enzyme_mechanism begin
            substrates: A, B
            products: P, Q
            steps: begin
                E + A <--> E(A)
                E(A) + B ⇌ E(A, B)
                E(A, B) <--> E(P, Q)
                E(Q) + P ⇌ E(P, Q)
                (E + Q ⇌ E(Q), E(A::Inh) + Q ⇌ E(A::Inh, Q))
                (E + A::Inh ⇌ E(A::Inh), E(Q) + A::Inh ⇌ E(A::Inh, Q))
            end
        end),
        EnzymeRates.Mechanism(@enzyme_mechanism begin
            substrates: A, B
            products: P, Q
            steps: begin
                E + A ⇌ E(A)
                E(A) + B ⇌ E(A, B)
                E(A, B) <--> E(P, Q)
                E(Q) + P ⇌ E(P, Q)
                (E + Q <--> E(Q), E(A::Inh) + Q <--> E(A::Inh, Q))
                (E + A::Inh ⇌ E(A::Inh), E(Q) + A::Inh ⇌ E(A::Inh, Q))
            end
        end),
    ]
    kids = EnzymeRates._expand_re_to_ss(m)
    @test length(kids) == 2
    @test Set(kids) == Set(expected)
    r0 = _testhelper_identifiable_rank(m)
    for absent in [
        EnzymeRates.Mechanism(@enzyme_mechanism begin      # the B flank
            substrates: A, B
            products: P, Q
            steps: begin
                E + A ⇌ E(A)
                E(A) + B <--> E(A, B)
                E(A, B) <--> E(P, Q)
                E(Q) + P ⇌ E(P, Q)
                (E + Q ⇌ E(Q), E(A::Inh) + Q ⇌ E(A::Inh, Q))
                (E + A::Inh ⇌ E(A::Inh), E(Q) + A::Inh ⇌ E(A::Inh, Q))
            end
        end),
        EnzymeRates.Mechanism(@enzyme_mechanism begin      # the P flank
            substrates: A, B
            products: P, Q
            steps: begin
                E + A ⇌ E(A)
                E(A) + B ⇌ E(A, B)
                E(A, B) <--> E(P, Q)
                E(Q) + P <--> E(P, Q)
                (E + Q ⇌ E(Q), E(A::Inh) + Q ⇌ E(A::Inh, Q))
                (E + A::Inh ⇌ E(A::Inh), E(Q) + A::Inh ⇌ E(A::Inh, Q))
            end
        end)
    ]
        @test !(absent in kids)
        @test _testhelper_identifiable_rank(absent) == r0
    end
    for r in kids
        _testhelper_assert_mechanism_invariants(r)
        @test EnzymeRates.compile_mechanism(r) isa EnzymeMechanism
    end
    for r in kids, grp in EnzymeRates.steps(r), s in grp
        EnzymeRates.bound_metabolite(s) isa EnzymeRates.Regulator &&
            @test EnzymeRates.is_equilibrium(s)
    end
end

@testset "AllostericMechanism — the inhibitor copy never flips; tags are preserved" begin
    # The same mechanism as an allosteric one with a dead inactive conformation
    # (the chemistry `:OnlyA`, every binding `:EqualAI`) and one catalytic subunit.
    # All four groups A, B, P and Q flip: an allosteric mechanism has no chain
    # flanks (`_chain_flank_groups`). Every child keeps one `:OnlyA` group, the
    # chemistry, and the multiplicity and (empty) regulatory sites.
    am = EnzymeRates.AllostericMechanism(@allosteric_mechanism begin
        substrates: A, B
        products: P, Q
        catalytic_inhibitors: A
        catalytic_multiplicity: 1
        catalytic_steps: begin
            E + A ⇌ E(A)                                            :: EqualAI
            E(A) + B ⇌ E(A, B)                                      :: EqualAI
            E(A, B) <--> E(P, Q)                                    :: OnlyA
            E(Q) + P ⇌ E(P, Q)                                      :: EqualAI
            (E + Q ⇌ E(Q), E(A::Inh) + Q ⇌ E(A::Inh, Q))            :: EqualAI
            (E + A::Inh ⇌ E(A::Inh), E(Q) + A::Inh ⇌ E(A::Inh, Q))  :: EqualAI
        end
    end)
    lift(em) = EnzymeRates.AllostericMechanism(em)
    expected = [
        lift(@allosteric_mechanism begin
            substrates: A, B
            products: P, Q
            catalytic_inhibitors: A
            catalytic_multiplicity: 1
            catalytic_steps: begin
                E + A <--> E(A)                                         :: EqualAI
                E(A) + B ⇌ E(A, B)                                      :: EqualAI
                E(A, B) <--> E(P, Q)                                    :: OnlyA
                E(Q) + P ⇌ E(P, Q)                                      :: EqualAI
                (E + Q ⇌ E(Q), E(A::Inh) + Q ⇌ E(A::Inh, Q))            :: EqualAI
                (E + A::Inh ⇌ E(A::Inh), E(Q) + A::Inh ⇌ E(A::Inh, Q))  :: EqualAI
            end
        end),
        lift(@allosteric_mechanism begin
            substrates: A, B
            products: P, Q
            catalytic_inhibitors: A
            catalytic_multiplicity: 1
            catalytic_steps: begin
                E + A ⇌ E(A)                                            :: EqualAI
                E(A) + B <--> E(A, B)                                   :: EqualAI
                E(A, B) <--> E(P, Q)                                    :: OnlyA
                E(Q) + P ⇌ E(P, Q)                                      :: EqualAI
                (E + Q ⇌ E(Q), E(A::Inh) + Q ⇌ E(A::Inh, Q))            :: EqualAI
                (E + A::Inh ⇌ E(A::Inh), E(Q) + A::Inh ⇌ E(A::Inh, Q))  :: EqualAI
            end
        end),
        lift(@allosteric_mechanism begin
            substrates: A, B
            products: P, Q
            catalytic_inhibitors: A
            catalytic_multiplicity: 1
            catalytic_steps: begin
                E + A ⇌ E(A)                                            :: EqualAI
                E(A) + B ⇌ E(A, B)                                      :: EqualAI
                E(A, B) <--> E(P, Q)                                    :: OnlyA
                E(Q) + P <--> E(P, Q)                                   :: EqualAI
                (E + Q ⇌ E(Q), E(A::Inh) + Q ⇌ E(A::Inh, Q))            :: EqualAI
                (E + A::Inh ⇌ E(A::Inh), E(Q) + A::Inh ⇌ E(A::Inh, Q))  :: EqualAI
            end
        end),
        lift(@allosteric_mechanism begin
            substrates: A, B
            products: P, Q
            catalytic_inhibitors: A
            catalytic_multiplicity: 1
            catalytic_steps: begin
                E + A ⇌ E(A)                                            :: EqualAI
                E(A) + B ⇌ E(A, B)                                      :: EqualAI
                E(A, B) <--> E(P, Q)                                    :: OnlyA
                E(Q) + P ⇌ E(P, Q)                                      :: EqualAI
                (E + Q <--> E(Q), E(A::Inh) + Q <--> E(A::Inh, Q))      :: EqualAI
                (E + A::Inh ⇌ E(A::Inh), E(Q) + A::Inh ⇌ E(A::Inh, Q))  :: EqualAI
            end
        end),
    ]
    kids = EnzymeRates._expand_re_to_ss(am)
    @test length(kids) == 4
    @test Set(kids) == Set(expected)
    for r in kids
        _testhelper_assert_mechanism_invariants(r)
        @test EnzymeRates.compile_mechanism(r) isa EnzymeRates.AllostericEnzymeMechanism
    end
    for r in kids
        @test count(==(:OnlyA), r.cat_allo_states) == 1
        onlya = findfirst(==(:OnlyA), r.cat_allo_states)
        @test EnzymeRates.is_iso(first(r.cat_steps[onlya]))
        @test r.catalytic_multiplicity == am.catalytic_multiplicity
        @test r.regulatory_sites == am.regulatory_sites
        @test EnzymeRates.reaction(r) == EnzymeRates.reaction(am)
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
                                EnzymeRates.consumed(s), EnzymeRates.released(s), false)
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

@testset "Mechanism — all-SS seed: empty (negative)" begin
    # If every catalytic group is already SS, the move has no RE
    # group to fire on → empty.
    m_all_ss = EnzymeRates.Mechanism(@enzyme_mechanism begin
        substrates: S
        products: P
        steps: begin
            E + P <--> E(P)
            E + S <--> E(S)
            E(S) <--> E(P)
        end
    end)
    @test isempty(EnzymeRates._expand_re_to_ss(m_all_ss))
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
    _testhelper_assert_mechanism_invariants(m)

    result = EnzymeRates._expand_split_kinetic_group(m)
    @test !isempty(result)
    base = EnzymeRates._independent_param_count(m)
    for r in result
        @test EnzymeRates._independent_param_count(r) > base
        @test r isa EnzymeRates.Mechanism
        _testhelper_assert_mechanism_invariants(r)
        @test EnzymeRates.compile_mechanism(r) isa EnzymeMechanism
        @test sum(length, EnzymeRates.steps(r)) == sum(length, EnzymeRates.steps(m))
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
    _testhelper_assert_mechanism_invariants(am)

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
        _testhelper_assert_mechanism_invariants(r)
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

@testset "Mechanism — a split part without flux is emitted at rapid equilibrium" begin
    # Ordered bi-bi whose steady-state B group holds the catalytic binding
    # E(A) + B and the abortive binding E(Q) + B. Splitting by context leaves
    # E(Q) + B → E(B, Q) alone: E(B, Q) is a dead end, so the step is a pendant
    # bridge of the segment graph and carries no flux. Its two constants would
    # enter the rate only as their ratio, so the part is emitted at rapid
    # equilibrium, one constant fewer than the raw split and the same family.
    # E(B, Q) lies on no cycle, so the new constant raises the independent count: one
    # child.
    m = EnzymeRates.Mechanism(@enzyme_mechanism begin
        substrates: A, B
        products: P, Q
        steps: begin
            E + A ⇌ E(A)
            (E(A) + B <--> E(A, B), E(Q) + B <--> E(B, Q))
            E(A, B) <--> E(P, Q)
            E(Q) + P ⇌ E(P, Q)
            E + Q ⇌ E(Q)
        end
    end)
    reverted = EnzymeRates.Mechanism(@enzyme_mechanism begin
        substrates: A, B
        products: P, Q
        steps: begin
            E + A ⇌ E(A)
            E(A) + B <--> E(A, B)
            E(Q) + B ⇌ E(B, Q)
            E(A, B) <--> E(P, Q)
            E(Q) + P ⇌ E(P, Q)
            E + Q ⇌ E(Q)
        end
    end)
    raw = EnzymeRates.Mechanism(@enzyme_mechanism begin
        substrates: A, B
        products: P, Q
        steps: begin
            E + A ⇌ E(A)
            E(A) + B <--> E(A, B)
            E(Q) + B <--> E(B, Q)
            E(A, B) <--> E(P, Q)
            E(Q) + P ⇌ E(P, Q)
            E + Q ⇌ E(Q)
        end
    end)
    kids = EnzymeRates._expand_split_kinetic_group(m)
    @test length(kids) == 1
    @test Set(kids) == Set([reverted])
    @test !(raw in kids)
    @test EnzymeRates._independent_param_count(reverted) ==
        EnzymeRates._independent_param_count(m) + 1
    r_rev = _testhelper_identifiable_rank(reverted)
    r_raw = _testhelper_identifiable_rank(raw)
    @test r_rev == r_raw                                          # the same family
    @test _testhelper_fitted(reverted) == _testhelper_fitted(raw) - 1  # one constant fewer
    # Both children keep one phantom of another class: once the abortive step no
    # longer pins kf_B/kr_B, E(A) + B → E(A, B) is a steady-state binding into a
    # form with one exit, and its three constants enter the law through two
    # combinations. The revert removes exactly the zero-flux phantom.
    @test _testhelper_fitted(raw) - r_raw == _testhelper_fitted(reverted) - r_rev + 1
end

@testset "Mechanism — a split leaving a redundant copy part is not emitted" begin
    # Ordered bi-bi with B also bound as a competitive-inhibitor copy at E(A) and
    # E(Q). E(A, B::Inh) has the composition of E(A, B); E(Q, B::Inh) is new, so
    # the placement stands. Splitting the copy group by context would leave
    # E(A) + B::Inh alone, a group that binds only where it duplicates E(A, B), and
    # E(A, B) is formed by the B group alone and left by the chemistry alone, so the
    # gauge exists: the child would be its parent plus a phantom. Every other group
    # holds one step, so the move emits nothing.
    m = EnzymeRates.Mechanism(@enzyme_mechanism begin
        substrates: A, B
        products: P, Q
        steps: begin
            E + A ⇌ E(A)
            E(A) + B ⇌ E(A, B)
            E(A, B) <--> E(P, Q)
            E(Q) + P ⇌ E(P, Q)
            E + Q ⇌ E(Q)
            (E(A) + B::Inh ⇌ E(A, B::Inh), E(Q) + B::Inh ⇌ E(B::Inh, Q))
        end
    end)
    @test isempty(EnzymeRates._expand_split_kinetic_group(m))
    absent = EnzymeRates.Mechanism(@enzyme_mechanism begin
        substrates: A, B
        products: P, Q
        steps: begin
            E + A ⇌ E(A)
            E(A) + B ⇌ E(A, B)
            E(A, B) <--> E(P, Q)
            E(Q) + P ⇌ E(P, Q)
            E + Q ⇌ E(Q)
            E(A) + B::Inh ⇌ E(A, B::Inh)
            E(Q) + B::Inh ⇌ E(B::Inh, Q)
        end
    end)
    @test EnzymeRates._independent_param_count(absent) ==
        EnzymeRates._independent_param_count(m) + 1
    @test _testhelper_identifiable_rank(absent) == _testhelper_identifiable_rank(m)
end

@testset "Mechanism — a split that makes a copy group redundant is not emitted" begin
    # The shared A group with abortive E(A, Q) and A as its own inhibitor at E: the
    # copy duplicates E(A), but the A group also binds at E(Q), so the gauge fails and
    # the mechanism is valid. Splitting the A group alone or the Q group alone gains
    # nothing (the cycle E, E(A), E(A, Q), E(Q) ties the new constant back); splitting
    # both frees one. In that child the A binding at E sits alone in its group and the
    # Q binding at E(A) in its own: every group holds one ratio, the gauge exists, and
    # the copy is redundant (7 fitted, rank 6). Its family is that of the split parent
    # without the copy, so the move emits nothing.
    m = EnzymeRates.Mechanism(@enzyme_mechanism begin
        substrates: A, B
        products: P, Q
        steps: begin
            (E + A ⇌ E(A), E(Q) + A ⇌ E(A, Q))
            E(A) + B ⇌ E(A, B)
            E(A, B) <--> E(P, Q)
            E(Q) + P ⇌ E(P, Q)
            (E + Q ⇌ E(Q), E(A) + Q ⇌ E(A, Q))
            E + A::Inh ⇌ E(A::Inh)
        end
    end)
    absent = EnzymeRates.Mechanism(@enzyme_mechanism begin
        substrates: A, B
        products: P, Q
        steps: begin
            E + A ⇌ E(A)
            E(Q) + A ⇌ E(A, Q)
            E(A) + B ⇌ E(A, B)
            E(A, B) <--> E(P, Q)
            E(Q) + P ⇌ E(P, Q)
            E + Q ⇌ E(Q)
            E(A) + Q ⇌ E(A, Q)
            E + A::Inh ⇌ E(A::Inh)
        end
    end)
    @test EnzymeRates._assert_emission_rules(m) === nothing
    @test isempty(EnzymeRates._expand_split_kinetic_group(m))
    @test !isempty(EnzymeRates._redundant_copy_groups(absent))
    @test EnzymeRates._independent_param_count(absent) ==
        EnzymeRates._independent_param_count(m) + 1
    @test _testhelper_fitted(absent) == 7 && _testhelper_identifiable_rank(absent) == 6
end

@testset "Mechanism — two copies split down to complexes of one composition" begin
    # Ordered bi-bi with Q bound as a competitive-inhibitor copy at E(A) and at E(Q) in two
    # groups, and A as a copy at E(A) and E(Q) in one group. E(Q, A::Inh) has the
    # composition {A, Q} of E(A, Q::Inh), a copy's complex, and no productive form has
    # it, so neither part of the A copy's context bipartition is twin-only; every other
    # group holds one step, and the move emits the pair. The child carries one phantom:
    # with E(Q) + A::Inh alone in its group, its complex and E(A, Q::Inh) enter the law
    # through one coefficient, K_A·K_Q* + K_Q·K_A*, so 9 fitted constants have rank 8.
    # The rule tolerates it: rejecting it would also reject identifiable splits whose
    # other copy is pinned.
    m = EnzymeRates.Mechanism(@enzyme_mechanism begin
        substrates: A, B
        products: P, Q
        steps: begin
            E + A ⇌ E(A)
            E(A) + B ⇌ E(A, B)
            E(A, B) <--> E(P, Q)
            E(Q) + P ⇌ E(P, Q)
            E + Q ⇌ E(Q)
            E(A) + Q::Inh ⇌ E(A, Q::Inh)
            E(Q) + Q::Inh ⇌ E(Q, Q::Inh)
            (E(A) + A::Inh ⇌ E(A, A::Inh), E(Q) + A::Inh ⇌ E(A::Inh, Q))
        end
    end)
    emitted = EnzymeRates.Mechanism(@enzyme_mechanism begin
        substrates: A, B
        products: P, Q
        steps: begin
            E + A ⇌ E(A)
            E(A) + B ⇌ E(A, B)
            E(A, B) <--> E(P, Q)
            E(Q) + P ⇌ E(P, Q)
            E + Q ⇌ E(Q)
            E(A) + Q::Inh ⇌ E(A, Q::Inh)
            E(Q) + Q::Inh ⇌ E(Q, Q::Inh)
            E(A) + A::Inh ⇌ E(A, A::Inh)
            E(Q) + A::Inh ⇌ E(A::Inh, Q)
        end
    end)
    kids = EnzymeRates._expand_split_kinetic_group(m)
    @test length(kids) == 1
    @test Set(kids) == Set([emitted])
    @test _testhelper_fitted(emitted) == _testhelper_identifiable_rank(emitted) + 1
    @test _testhelper_fitted(m) == _testhelper_identifiable_rank(m)
end

end

# ─── _expand_add_dead_end_regulator ────────────────────────────────────
@testset "_expand_add_dead_end_regulator" begin

@testset "Mechanism — bi-bi + I: variants" begin
    for (em_seed, rxn, n_variants, multi_form) in [
        # SEED: bi-bi sequential. The expansion should produce 4 form sets after dedup.
        (@enzyme_mechanism(begin
            substrates: A, B
            products: P, Q
            steps: begin
                E + A ⇌ E(A)
                E(A) + B ⇌ E(A, B)
                E + Q ⇌ E(Q)
                E(Q) + P ⇌ E(P, Q)
                E(A, B) <--> E(P, Q)
            end
        end), @enzyme_reaction(begin
            substrates: A[C], B[N]
            products: P[C], Q[N]
            dead_end_inhibitors: I
        end), 4, false),
        # SEED: bi-bi random. It has 5 eligible forms (E, E_A, E_B, E_P, E_Q);
        # competition patterns × dedup → 9 variants. Some variants bind I at several
        # forms.
        (@enzyme_mechanism(begin
            substrates: A, B
            products: P, Q
            steps: begin
                (E + A ⇌ E(A), E(B) + A ⇌ E(A, B))
                (E + B ⇌ E(B), E(A) + B ⇌ E(A, B))
                (E + P ⇌ E(P), E(Q) + P ⇌ E(P, Q))
                (E + Q ⇌ E(Q), E(P) + Q ⇌ E(P, Q))
                E(A, B) <--> E(P, Q)
            end
        end), @enzyme_reaction(begin
            substrates: A[C], B[N]
            products: P[C], Q[N]
            dead_end_inhibitors: I
        end), 9, true),
        # SEED: bi-bi ping-pong with Estar. It has 4 eligible forms (E, E_A, Estar,
        # E_Q); competition patterns × dedup → 3 variants.
        (@enzyme_mechanism(begin
            substrates: A, B
            products: P, Q
            steps: begin
                E + A ⇌ E(A)
                Estar + B ⇌ Estar(B)
                E + Q ⇌ E(Q)
                Estar + P ⇌ Estar(P)
                E(A) <--> Estar(P)
                Estar(B) ⇌ E(Q)
            end
        end), @enzyme_reaction(begin
            substrates: A[CX], B[N]
            products: P[C], Q[NX]
            dead_end_inhibitors: I
        end), 3, false),
    ]
        m = EnzymeRates.Mechanism(em_seed)
        _testhelper_assert_mechanism_invariants(m)

        result = EnzymeRates._expand_add_dead_end_regulator(m, rxn)

        # 1. count
        @test length(result) == n_variants

        # 2. Δ params: +1 each (one new K_I parameter), measured against ground-truth
        # `fitted_params(compile_mechanism(...))` — the canonical source for exact
        # parameter counts.
        @test _testhelper_param_deltas(m, result) == fill(1, n_variants)

        # 3. compilability
        for r in result
            _testhelper_assert_mechanism_invariants(r)
            @test EnzymeRates.EnzymeMechanism(r) isa EnzymeMechanism
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
    # Ordered Theorell–Chance, all SS: E + A → E(A), E(A) + B → E(Q) + P, E + Q → E(Q).
    # Sites: A is free at E, B at E(A) (taken up by the TC step), P at E(Q) (given off
    # there), Q at E. A foreign inhibitor I competes with one or both substrates and one
    # or both products; a site is skipped when it already holds a competing reactant.
    # {A}×{P}: E, E(Q). {A}×{Q}: E. {A}×{P,Q}: E (E(Q) holds Q). {B}×{P}: E(A), E(Q).
    # {B}×{Q}: E, E(A). {B}×{P,Q}: E, E(A) (E(Q) holds Q). {A,B}×{P}: E, E(Q) (E(A)
    # holds A). {A,B}×{Q}: E. {A,B}×{P,Q}: E. Distinct placements: {E}, {E, E(Q)},
    # {E(A), E(Q)}, {E, E(A)}. Mirrors join two sites: E–E(Q) by Q's binding, E–E(A) by
    # A's, E(A)–E(Q) by the TC step, each in its parent's group, at its parent's flag.
    # A foreign inhibitor has no twins, so the copy rule keeps every placement.
    rxn = @enzyme_reaction begin
        substrates: A[C], B[N]
        products: P[C], Q[N]
        dead_end_inhibitors: I
    end
    withrxn(em) = EnzymeRates.Mechanism(rxn, EnzymeRates.steps(EnzymeRates.Mechanism(em)))
    m = EnzymeRates.Mechanism(rxn, EnzymeRates.steps(EnzymeRates.Mechanism(
        @enzyme_mechanism begin
            substrates: A, B
            products: P, Q
            steps: begin
                E + A <--> E(A)
                E(A) + B <--> E(Q) + P
                E + Q <--> E(Q)
            end
        end)))
    at_E = withrxn(@enzyme_mechanism begin
        substrates: A, B; products: P, Q; regulators: I
        steps: begin
            E + A <--> E(A)
            E(A) + B <--> E(Q) + P
            E + Q <--> E(Q)
            E + I ⇌ E(I)
        end
    end)
    at_E_EQ = withrxn(@enzyme_mechanism begin
        substrates: A, B; products: P, Q; regulators: I
        steps: begin
            E + A <--> E(A)
            E(A) + B <--> E(Q) + P
            (E + Q <--> E(Q), E(I) + Q <--> E(I, Q))
            (E + I ⇌ E(I), E(Q) + I ⇌ E(I, Q))
        end
    end)
    at_EA_EQ = withrxn(@enzyme_mechanism begin
        substrates: A, B; products: P, Q; regulators: I
        steps: begin
            E + A <--> E(A)
            (E(A) + B <--> E(Q) + P, E(A, I) + B <--> E(I, Q) + P)
            E + Q <--> E(Q)
            (E(A) + I ⇌ E(A, I), E(Q) + I ⇌ E(I, Q))
        end
    end)
    at_E_EA = withrxn(@enzyme_mechanism begin
        substrates: A, B; products: P, Q; regulators: I
        steps: begin
            (E + A <--> E(A), E(I) + A <--> E(A, I))
            E(A) + B <--> E(Q) + P
            E + Q <--> E(Q)
            (E + I ⇌ E(I), E(A) + I ⇌ E(A, I))
        end
    end)
    kids = EnzymeRates._expand_add_dead_end_regulator(m, rxn)
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
    _testhelper_assert_mechanism_invariants(m)
    rxn = @enzyme_reaction begin
        substrates: S[C]
        products: P[C]
        dead_end_inhibitors: I, J
    end

    # Step A: 2 eligible regs (I, J), 1 form each → 2 variants total.
    i_or_j_ms = EnzymeRates._expand_add_dead_end_regulator(m, rxn)
    @test length(i_or_j_ms) == 2
    with_i = first(filter(r -> !isempty(_testhelper_groups_binding(r, :I)), i_or_j_ms))

    # Step B: add J on top of the I-bound variant.
    j_ms = EnzymeRates._expand_add_dead_end_regulator(with_i, rxn)
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
    _testhelper_assert_mechanism_invariants(m)
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

@testset "Mechanism — uni-uni: a reactant as its own inhibitor has no placement" begin
    # The only site that carries neither S nor P is free E, and E·S* has the
    # composition of E(S) (E·P* that of E(P)), formed by the S group alone: the
    # gauge exists, and the copy's constant would enter the rate only as
    # K_S + K_S,inh with kcat rescaled. No child is emitted.
    rxn_S = @enzyme_reaction begin
        substrates: S[C]
        products: P[C]
        dead_end_inhibitors: S
    end
    rxn_P = @enzyme_reaction begin
        substrates: S[C]
        products: P[C]
        dead_end_inhibitors: P
    end
    for rxn in (rxn_S, rxn_P)
        m = EnzymeRates.Mechanism(@enzyme_mechanism begin
            substrates: S
            products: P
            steps: begin
                E + S ⇌ E(S)
                E(S) <--> E(P)
                E + P ⇌ E(P)
            end
        end)
        @test isempty(EnzymeRates._expand_add_dead_end_regulator(m, rxn))
    end
end

@testset "Mechanism — ordered bi-bi with A as its own inhibitor: three placements" begin
    # A binds productively at E only, so A* competing with P puts it at {E, E(Q)},
    # with Q puts it at {E} alone, and competing with B puts it at {E(A), E(Q)} or
    # {E, E(A)}. E·A* has the composition of E(A); E(Q)·A* and E(A)·A* are new. The
    # placement at E alone is all-twin, and its gauge exists: E(A) is formed by the
    # A group alone, so the A group's ratio is 1/ρ and the B group's ρ, each over
    # one step. It is redundant and is not emitted. Mirrors: E + Q ⇌ E(Q) onto
    # {E, E(Q)} in the Q group; E + A ⇌ E(A) onto {E, E(A)} in the A group.
    rxn = @enzyme_reaction begin
        substrates: A[C], B[N]
        products: P[C], Q[N]
        dead_end_inhibitors: A
    end
    m = _testhelper_on_reaction(rxn, @enzyme_mechanism begin
        substrates: A, B
        products: P, Q
        steps: begin
            E + A ⇌ E(A)
            E(A) + B ⇌ E(A, B)
            E(A, B) <--> E(P, Q)
            E(Q) + P ⇌ E(P, Q)
            E + Q ⇌ E(Q)
        end
    end)
    at_E_EQ = _testhelper_on_reaction(rxn, @enzyme_mechanism begin
        substrates: A, B
        products: P, Q
        steps: begin
            E + A ⇌ E(A)
            E(A) + B ⇌ E(A, B)
            E(A, B) <--> E(P, Q)
            E(Q) + P ⇌ E(P, Q)
            (E + Q ⇌ E(Q), E(A::Inh) + Q ⇌ E(A::Inh, Q))
            (E + A::Inh ⇌ E(A::Inh), E(Q) + A::Inh ⇌ E(A::Inh, Q))
        end
    end)
    at_EA_EQ = _testhelper_on_reaction(rxn, @enzyme_mechanism begin
        substrates: A, B
        products: P, Q
        steps: begin
            E + A ⇌ E(A)
            E(A) + B ⇌ E(A, B)
            E(A, B) <--> E(P, Q)
            E(Q) + P ⇌ E(P, Q)
            E + Q ⇌ E(Q)
            (E(A) + A::Inh ⇌ E(A, A::Inh), E(Q) + A::Inh ⇌ E(A::Inh, Q))
        end
    end)
    at_E_EA = _testhelper_on_reaction(rxn, @enzyme_mechanism begin
        substrates: A, B
        products: P, Q
        steps: begin
            (E + A ⇌ E(A), E(A::Inh) + A ⇌ E(A, A::Inh))
            E(A) + B ⇌ E(A, B)
            E(A, B) <--> E(P, Q)
            E(Q) + P ⇌ E(P, Q)
            E + Q ⇌ E(Q)
            (E + A::Inh ⇌ E(A::Inh), E(A) + A::Inh ⇌ E(A, A::Inh))
        end
    end)
    at_E = _testhelper_on_reaction(rxn, @enzyme_mechanism begin
        substrates: A, B
        products: P, Q
        steps: begin
            E + A ⇌ E(A)
            E(A) + B ⇌ E(A, B)
            E(A, B) <--> E(P, Q)
            E(Q) + P ⇌ E(P, Q)
            E + Q ⇌ E(Q)
            E + A::Inh ⇌ E(A::Inh)
        end
    end)
    kids = EnzymeRates._expand_add_dead_end_regulator(m, rxn)
    @test length(kids) == 3
    @test Set(kids) == Set([at_E_EQ, at_EA_EQ, at_E_EA])
    # The absent placement adds a constant the data cannot see.
    @test !(at_E in kids)
    # Each child binds the copy as a CompetitiveInhibitor named :A, in one new
    # kinetic group, and adds exactly one fitted constant; every child has full rank.
    base_fitted = _testhelper_fitted(m)
    is_copy(bm) = bm isa EnzymeRates.CompetitiveInhibitor && EnzymeRates.name(bm) === :A
    for r in kids
        @test any(grp -> is_copy(EnzymeRates.bound_metabolite(first(grp))),
                  EnzymeRates.steps(r))
        @test length(EnzymeRates.steps(r)) == length(EnzymeRates.steps(m)) + 1
        fitted = _testhelper_fitted(r)
        @test fitted == base_fitted + 1
        @test fitted == _testhelper_identifiable_rank(r)
    end
end

@testset "Mechanism — shared A group with abortive E(A, Q): the gauge decides" begin
    # A binds at E and, in the same kinetic group, at E(Q) to form the abortive
    # complex E(A, Q). A copy of A at E duplicates E(A) and at E(Q) duplicates
    # E(A, Q), so the placements {E} and {E, E(Q)} are all-twin, and the dwell gauge
    # decides. At {E, E(Q)}, with E + Q ⇌ E(Q) mirrored onto the copy forms in the Q
    # group, both twins come from the one A group and share ρ: every step of the A
    # group has ratio 1/ρ and every step of the Q group ratio 1, so the gauge exists
    # and the copy is redundant (not emitted; 6 fitted, rank 5). At {E} the A group
    # holds E + A ⇌ E(A), ratio 1/ρ, beside E(Q) + A ⇌ E(A, Q), ratio 1: no gauge,
    # since the shared group pins K_A at E(Q), where no copy binds. The child is
    # emitted and has full rank (6 fitted, rank 6). {E(A), E(Q)} and {E, E(A)} create
    # E(A, A*) and are emitted; the second mirrors E + A ⇌ E(A) in the A group.
    rxn = @enzyme_reaction begin
        substrates: A[C], B[N]
        products: P[C], Q[N]
        dead_end_inhibitors: A
    end
    m = _testhelper_on_reaction(rxn, @enzyme_mechanism begin
        substrates: A, B
        products: P, Q
        steps: begin
            (E + A ⇌ E(A), E(Q) + A ⇌ E(A, Q))
            E(A) + B ⇌ E(A, B)
            E(A, B) <--> E(P, Q)
            E(Q) + P ⇌ E(P, Q)
            (E + Q ⇌ E(Q), E(A) + Q ⇌ E(A, Q))
        end
    end)
    at_EA_EQ = _testhelper_on_reaction(rxn, @enzyme_mechanism begin
        substrates: A, B
        products: P, Q
        steps: begin
            (E + A ⇌ E(A), E(Q) + A ⇌ E(A, Q))
            E(A) + B ⇌ E(A, B)
            E(A, B) <--> E(P, Q)
            E(Q) + P ⇌ E(P, Q)
            (E + Q ⇌ E(Q), E(A) + Q ⇌ E(A, Q))
            (E(A) + A::Inh ⇌ E(A, A::Inh), E(Q) + A::Inh ⇌ E(A::Inh, Q))
        end
    end)
    at_E_EA = _testhelper_on_reaction(rxn, @enzyme_mechanism begin
        substrates: A, B
        products: P, Q
        steps: begin
            (E + A ⇌ E(A), E(Q) + A ⇌ E(A, Q), E(A::Inh) + A ⇌ E(A, A::Inh))
            E(A) + B ⇌ E(A, B)
            E(A, B) <--> E(P, Q)
            E(Q) + P ⇌ E(P, Q)
            (E + Q ⇌ E(Q), E(A) + Q ⇌ E(A, Q))
            (E + A::Inh ⇌ E(A::Inh), E(A) + A::Inh ⇌ E(A, A::Inh))
        end
    end)
    at_E = _testhelper_on_reaction(rxn, @enzyme_mechanism begin
        substrates: A, B
        products: P, Q
        steps: begin
            (E + A ⇌ E(A), E(Q) + A ⇌ E(A, Q))
            E(A) + B ⇌ E(A, B)
            E(A, B) <--> E(P, Q)
            E(Q) + P ⇌ E(P, Q)
            (E + Q ⇌ E(Q), E(A) + Q ⇌ E(A, Q))
            E + A::Inh ⇌ E(A::Inh)
        end
    end)
    at_E_EQ = _testhelper_on_reaction(rxn, @enzyme_mechanism begin
        substrates: A, B
        products: P, Q
        steps: begin
            (E + A ⇌ E(A), E(Q) + A ⇌ E(A, Q))
            E(A) + B ⇌ E(A, B)
            E(A, B) <--> E(P, Q)
            E(Q) + P ⇌ E(P, Q)
            (E + Q ⇌ E(Q), E(A) + Q ⇌ E(A, Q), E(A::Inh) + Q ⇌ E(A::Inh, Q))
            (E + A::Inh ⇌ E(A::Inh), E(Q) + A::Inh ⇌ E(A::Inh, Q))
        end
    end)
    kids = EnzymeRates._expand_add_dead_end_regulator(m, rxn)
    @test length(kids) == 3
    @test Set(kids) == Set([at_EA_EQ, at_E_EA, at_E])
    @test !(at_E_EQ in kids)
    @test _testhelper_identifiable_rank(m) == _testhelper_fitted(m)
end

@testset "Mechanism — a second copy may share an older copy's composition" begin
    # Ordered bi-bi with Q as its own inhibitor at {E, E(A)}, the A binding
    # mirrored onto the copy forms. Competition is decided per site: Q* binds a
    # dead-end site of its own, so E(Q*) carries no productive Q. A binds
    # productively at E and E(Q*), B at E(A), P at E(Q), Q at E; Q* binds at E and
    # E(A). A copy of A goes to {E(A), E(Q)} competing with B and P, to {E, E(A)}
    # competing with B and Q, to {E, E(Q), E(Q*)} competing with A and P, to
    # {E, E(Q*)} competing with A and Q, to {E, E(Q)} competing with A, P and Q*, to
    # {E, E(A), E(Q)} competing with B, P and Q*, and to {E} competing with A, Q and
    # Q*. E·A* duplicates E(A), so {E} is all-twin; but the A group that forms E(A)
    # also binds at E(Q*), where no A* binds, so it holds a step of ratio 1/ρ beside
    # one of ratio 1, the gauge fails, and {E} is emitted. E(A)·A* is new, and
    # E(Q)·A* and E(Q*)·A* have the composition of E(A, Q*), a copy's complex, which
    # is no twin source. Where A* binds at E(A), E(A, A*) pins K_A* and the child has
    # full rank. In the child at {E(A), E(Q)} the older Q* group binds at E, a twin of
    # E(Q), and at E(A), whose complex E(A, Q*) has the composition of E(A*, Q); its
    # constant still shows, since the A·Q term is A·Q/(K_A·K_Q*) + A·Q/(K_Q·K_A*).
    # Where A* binds at E and otherwise, if anywhere, only where its complex has the
    # composition of E(A, Q*), in the children at {E, E(Q), E(Q*)}, {E, E(Q*)},
    # {E, E(Q)} and {E}, no site pins either copy: K_A, K_A*, K_Q and K_Q* enter the
    # law only through the A, Q and A·Q coefficients, so the child carries one
    # phantom. The rule tolerates it: it rejects a copy only on a proof.
    rxn = @enzyme_reaction begin
        substrates: A[C], B[N]
        products: P[C], Q[N]
        dead_end_inhibitors: A, Q
    end
    m = _testhelper_on_reaction(rxn, @enzyme_mechanism begin
        substrates: A, B
        products: P, Q
        steps: begin
            (E + A ⇌ E(A), E(Q::Inh) + A ⇌ E(A, Q::Inh))
            E(A) + B ⇌ E(A, B)
            E(A, B) <--> E(P, Q)
            E(Q) + P ⇌ E(P, Q)
            E + Q ⇌ E(Q)
            (E + Q::Inh ⇌ E(Q::Inh), E(A) + Q::Inh ⇌ E(A, Q::Inh))
        end
    end)
    at_EA_EQ = _testhelper_on_reaction(rxn, @enzyme_mechanism begin
        substrates: A, B
        products: P, Q
        steps: begin
            (E + A ⇌ E(A), E(Q::Inh) + A ⇌ E(A, Q::Inh))
            E(A) + B ⇌ E(A, B)
            E(A, B) <--> E(P, Q)
            E(Q) + P ⇌ E(P, Q)
            E + Q ⇌ E(Q)
            (E + Q::Inh ⇌ E(Q::Inh), E(A) + Q::Inh ⇌ E(A, Q::Inh))
            (E(A) + A::Inh ⇌ E(A, A::Inh), E(Q) + A::Inh ⇌ E(A::Inh, Q))
        end
    end)
    at_E_EA = _testhelper_on_reaction(rxn, @enzyme_mechanism begin
        substrates: A, B
        products: P, Q
        steps: begin
            (E + A ⇌ E(A), E(Q::Inh) + A ⇌ E(A, Q::Inh), E(A::Inh) + A ⇌ E(A, A::Inh))
            E(A) + B ⇌ E(A, B)
            E(A, B) <--> E(P, Q)
            E(Q) + P ⇌ E(P, Q)
            E + Q ⇌ E(Q)
            (E + Q::Inh ⇌ E(Q::Inh), E(A) + Q::Inh ⇌ E(A, Q::Inh))
            (E + A::Inh ⇌ E(A::Inh), E(A) + A::Inh ⇌ E(A, A::Inh))
        end
    end)
    at_E_EQ_EQinh = _testhelper_on_reaction(rxn, @enzyme_mechanism begin
        substrates: A, B
        products: P, Q
        steps: begin
            (E + A ⇌ E(A), E(Q::Inh) + A ⇌ E(A, Q::Inh))
            E(A) + B ⇌ E(A, B)
            E(A, B) <--> E(P, Q)
            E(Q) + P ⇌ E(P, Q)
            (E + Q ⇌ E(Q), E(A::Inh) + Q ⇌ E(A::Inh, Q))
            (E + Q::Inh ⇌ E(Q::Inh), E(A) + Q::Inh ⇌ E(A, Q::Inh),
             E(A::Inh) + Q::Inh ⇌ E(A::Inh, Q::Inh))
            (E + A::Inh ⇌ E(A::Inh), E(Q) + A::Inh ⇌ E(A::Inh, Q),
             E(Q::Inh) + A::Inh ⇌ E(A::Inh, Q::Inh))
        end
    end)
    at_E_EQinh = _testhelper_on_reaction(rxn, @enzyme_mechanism begin
        substrates: A, B
        products: P, Q
        steps: begin
            (E + A ⇌ E(A), E(Q::Inh) + A ⇌ E(A, Q::Inh))
            E(A) + B ⇌ E(A, B)
            E(A, B) <--> E(P, Q)
            E(Q) + P ⇌ E(P, Q)
            E + Q ⇌ E(Q)
            (E + Q::Inh ⇌ E(Q::Inh), E(A) + Q::Inh ⇌ E(A, Q::Inh),
             E(A::Inh) + Q::Inh ⇌ E(A::Inh, Q::Inh))
            (E + A::Inh ⇌ E(A::Inh), E(Q::Inh) + A::Inh ⇌ E(A::Inh, Q::Inh))
        end
    end)
    at_E_EQ = _testhelper_on_reaction(rxn, @enzyme_mechanism begin
        substrates: A, B
        products: P, Q
        steps: begin
            (E + A ⇌ E(A), E(Q::Inh) + A ⇌ E(A, Q::Inh))
            E(A) + B ⇌ E(A, B)
            E(A, B) <--> E(P, Q)
            E(Q) + P ⇌ E(P, Q)
            (E + Q ⇌ E(Q), E(A::Inh) + Q ⇌ E(A::Inh, Q))
            (E + Q::Inh ⇌ E(Q::Inh), E(A) + Q::Inh ⇌ E(A, Q::Inh))
            (E + A::Inh ⇌ E(A::Inh), E(Q) + A::Inh ⇌ E(A::Inh, Q))
        end
    end)
    at_E_EA_EQ = _testhelper_on_reaction(rxn, @enzyme_mechanism begin
        substrates: A, B
        products: P, Q
        steps: begin
            (E + A ⇌ E(A), E(Q::Inh) + A ⇌ E(A, Q::Inh), E(A::Inh) + A ⇌ E(A, A::Inh))
            E(A) + B ⇌ E(A, B)
            E(A, B) <--> E(P, Q)
            E(Q) + P ⇌ E(P, Q)
            (E + Q ⇌ E(Q), E(A::Inh) + Q ⇌ E(A::Inh, Q))
            (E + Q::Inh ⇌ E(Q::Inh), E(A) + Q::Inh ⇌ E(A, Q::Inh))
            (E + A::Inh ⇌ E(A::Inh), E(A) + A::Inh ⇌ E(A, A::Inh),
             E(Q) + A::Inh ⇌ E(A::Inh, Q))
        end
    end)
    at_E = _testhelper_on_reaction(rxn, @enzyme_mechanism begin
        substrates: A, B
        products: P, Q
        steps: begin
            (E + A ⇌ E(A), E(Q::Inh) + A ⇌ E(A, Q::Inh))
            E(A) + B ⇌ E(A, B)
            E(A, B) <--> E(P, Q)
            E(Q) + P ⇌ E(P, Q)
            E + Q ⇌ E(Q)
            (E + Q::Inh ⇌ E(Q::Inh), E(A) + Q::Inh ⇌ E(A, Q::Inh))
            E + A::Inh ⇌ E(A::Inh)
        end
    end)
    kids = EnzymeRates._expand_add_dead_end_regulator(m, rxn)
    @test length(kids) == 7
    @test Set(kids) == Set([at_EA_EQ, at_E_EA, at_E_EQ_EQinh, at_E_EQinh, at_E_EQ,
                            at_E_EA_EQ, at_E])
    for r in kids
        @test isempty(EnzymeRates._redundant_copy_groups(r))
        @test EnzymeRates._assert_emission_rules(r) === nothing
    end
    for k in (at_EA_EQ, at_E_EA, at_E_EA_EQ)
        @test _testhelper_identifiable_rank(k) == _testhelper_fitted(k)
    end
    for k in (at_E_EQ_EQinh, at_E_EQinh, at_E_EQ, at_E)
        @test _testhelper_identifiable_rank(k) == _testhelper_fitted(k) - 1
    end
end

@testset "Mechanism — a copy's form hosts a copy that competes with its reactant" begin
    # Ordered bi-bi with A as its own inhibitor at {E, E(Q)}, placing a copy of Q.
    # Competition is decided per site: A* binds a dead-end site of its own, so
    # E(A*) carries no productive A and is a site for a Q* that competes with A. A
    # binds productively at E, B at E(A), P at E(Q), Q at E and E(A*); A* binds at E
    # and E(Q). E(A, B) and E(P, Q) hold every substrate or every product and are no
    # site. Q* goes to {E, E(Q)} competing with A and P, to {E, E(A*)} competing
    # with A and Q, to {E(A), E(Q)} competing with B and P, to {E, E(A), E(A*)}
    # competing with B and Q, to {E, E(A), E(Q)} competing with B, P and A*, to
    # {E, E(A)} competing with B, Q and A*, and to {E} competing with A, Q and A*.
    # E·Q* duplicates E(Q), so {E} is all-twin; but the Q group that forms E(Q) also
    # binds at E(A*), where no Q* binds, so it holds a step of ratio 1/ρ beside one of
    # ratio 1, the gauge fails, and {E} is emitted. E(Q)·Q* is new, and so are
    # E(A)·Q* and E(A*)·Q*: E(A, Q*) and E(A*, Q*) have the composition of E(A*, Q),
    # a copy's complex, which is no twin source. Mirrors: E + A ⇌ E(A) onto E(Q*) in
    # A's group wherever E(A) is a site, E + Q ⇌ E(Q) in Q's group wherever E(Q) is,
    # and E + A* ⇌ E(A*) in A*'s group wherever E(A*) is. Where Q* binds at E(Q),
    # E(Q, Q*) pins K_Q*, and {E(A), E(Q)} has no twin site; those children have full
    # rank. Where Q* binds at E and otherwise, if anywhere, only at E(A) or E(A*), no
    # site pins K_Q* apart from K_Q, so the child carries one phantom, which the rule
    # tolerates: it rejects a copy only on a proof.
    rxn = @enzyme_reaction begin
        substrates: A[C], B[N]
        products: P[C], Q[N]
        dead_end_inhibitors: A, Q
    end
    m = _testhelper_on_reaction(rxn, @enzyme_mechanism begin
        substrates: A, B
        products: P, Q
        steps: begin
            E + A ⇌ E(A)
            E(A) + B ⇌ E(A, B)
            E(A, B) <--> E(P, Q)
            E(Q) + P ⇌ E(P, Q)
            (E + Q ⇌ E(Q), E(A::Inh) + Q ⇌ E(A::Inh, Q))
            (E + A::Inh ⇌ E(A::Inh), E(Q) + A::Inh ⇌ E(A::Inh, Q))
        end
    end)
    at_E_EQ = _testhelper_on_reaction(rxn, @enzyme_mechanism begin
        substrates: A, B
        products: P, Q
        steps: begin
            E + A ⇌ E(A)
            E(A) + B ⇌ E(A, B)
            E(A, B) <--> E(P, Q)
            E(Q) + P ⇌ E(P, Q)
            (E + Q ⇌ E(Q), E(A::Inh) + Q ⇌ E(A::Inh, Q), E(Q::Inh) + Q ⇌ E(Q, Q::Inh))
            (E + A::Inh ⇌ E(A::Inh), E(Q) + A::Inh ⇌ E(A::Inh, Q))
            (E + Q::Inh ⇌ E(Q::Inh), E(Q) + Q::Inh ⇌ E(Q, Q::Inh))
        end
    end)
    at_E_EAinh = _testhelper_on_reaction(rxn, @enzyme_mechanism begin
        substrates: A, B
        products: P, Q
        steps: begin
            E + A ⇌ E(A)
            E(A) + B ⇌ E(A, B)
            E(A, B) <--> E(P, Q)
            E(Q) + P ⇌ E(P, Q)
            (E + Q ⇌ E(Q), E(A::Inh) + Q ⇌ E(A::Inh, Q))
            (E + A::Inh ⇌ E(A::Inh), E(Q) + A::Inh ⇌ E(A::Inh, Q),
             E(Q::Inh) + A::Inh ⇌ E(A::Inh, Q::Inh))
            (E + Q::Inh ⇌ E(Q::Inh), E(A::Inh) + Q::Inh ⇌ E(A::Inh, Q::Inh))
        end
    end)
    at_EA_EQ = _testhelper_on_reaction(rxn, @enzyme_mechanism begin
        substrates: A, B
        products: P, Q
        steps: begin
            E + A ⇌ E(A)
            E(A) + B ⇌ E(A, B)
            E(A, B) <--> E(P, Q)
            E(Q) + P ⇌ E(P, Q)
            (E + Q ⇌ E(Q), E(A::Inh) + Q ⇌ E(A::Inh, Q))
            (E + A::Inh ⇌ E(A::Inh), E(Q) + A::Inh ⇌ E(A::Inh, Q))
            (E(A) + Q::Inh ⇌ E(A, Q::Inh), E(Q) + Q::Inh ⇌ E(Q, Q::Inh))
        end
    end)
    at_E_EA_EAinh = _testhelper_on_reaction(rxn, @enzyme_mechanism begin
        substrates: A, B
        products: P, Q
        steps: begin
            (E + A ⇌ E(A), E(Q::Inh) + A ⇌ E(A, Q::Inh))
            E(A) + B ⇌ E(A, B)
            E(A, B) <--> E(P, Q)
            E(Q) + P ⇌ E(P, Q)
            (E + Q ⇌ E(Q), E(A::Inh) + Q ⇌ E(A::Inh, Q))
            (E + A::Inh ⇌ E(A::Inh), E(Q) + A::Inh ⇌ E(A::Inh, Q),
             E(Q::Inh) + A::Inh ⇌ E(A::Inh, Q::Inh))
            (E + Q::Inh ⇌ E(Q::Inh), E(A) + Q::Inh ⇌ E(A, Q::Inh),
             E(A::Inh) + Q::Inh ⇌ E(A::Inh, Q::Inh))
        end
    end)
    at_E_EA_EQ = _testhelper_on_reaction(rxn, @enzyme_mechanism begin
        substrates: A, B
        products: P, Q
        steps: begin
            (E + A ⇌ E(A), E(Q::Inh) + A ⇌ E(A, Q::Inh))
            E(A) + B ⇌ E(A, B)
            E(A, B) <--> E(P, Q)
            E(Q) + P ⇌ E(P, Q)
            (E + Q ⇌ E(Q), E(A::Inh) + Q ⇌ E(A::Inh, Q), E(Q::Inh) + Q ⇌ E(Q, Q::Inh))
            (E + A::Inh ⇌ E(A::Inh), E(Q) + A::Inh ⇌ E(A::Inh, Q))
            (E + Q::Inh ⇌ E(Q::Inh), E(A) + Q::Inh ⇌ E(A, Q::Inh),
             E(Q) + Q::Inh ⇌ E(Q, Q::Inh))
        end
    end)
    at_E_EA = _testhelper_on_reaction(rxn, @enzyme_mechanism begin
        substrates: A, B
        products: P, Q
        steps: begin
            (E + A ⇌ E(A), E(Q::Inh) + A ⇌ E(A, Q::Inh))
            E(A) + B ⇌ E(A, B)
            E(A, B) <--> E(P, Q)
            E(Q) + P ⇌ E(P, Q)
            (E + Q ⇌ E(Q), E(A::Inh) + Q ⇌ E(A::Inh, Q))
            (E + A::Inh ⇌ E(A::Inh), E(Q) + A::Inh ⇌ E(A::Inh, Q))
            (E + Q::Inh ⇌ E(Q::Inh), E(A) + Q::Inh ⇌ E(A, Q::Inh))
        end
    end)
    at_E = _testhelper_on_reaction(rxn, @enzyme_mechanism begin
        substrates: A, B
        products: P, Q
        steps: begin
            E + A ⇌ E(A)
            E(A) + B ⇌ E(A, B)
            E(A, B) <--> E(P, Q)
            E(Q) + P ⇌ E(P, Q)
            (E + Q ⇌ E(Q), E(A::Inh) + Q ⇌ E(A::Inh, Q))
            (E + A::Inh ⇌ E(A::Inh), E(Q) + A::Inh ⇌ E(A::Inh, Q))
            E + Q::Inh ⇌ E(Q::Inh)
        end
    end)
    kids = EnzymeRates._expand_add_dead_end_regulator(m, rxn)
    @test length(kids) == 7
    @test Set(kids) == Set([at_E_EQ, at_E_EAinh, at_EA_EQ, at_E_EA_EAinh, at_E_EA_EQ,
                            at_E_EA, at_E])
    @test at_E_EAinh in kids
    for r in kids
        @test EnzymeRates._assert_emission_rules(r) === nothing
    end
    for k in (at_E_EQ, at_EA_EQ, at_E_EA_EQ)
        @test _testhelper_identifiable_rank(k) == _testhelper_fitted(k)
    end
    for k in (at_E_EAinh, at_E_EA_EAinh, at_E_EA, at_E)
        @test _testhelper_identifiable_rank(k) == _testhelper_fitted(k) - 1
    end
end

@testset "Mechanism — an inhibitor that competes with A binds the copy's form E(A*)" begin
    # Ordered bi-bi with A as its own inhibitor at {E, E(Q)}, placing a foreign
    # inhibitor I. A* binds a dead-end site of its own, so E(A*) carries no
    # productive A and is a site for an I that competes with A. A binds productively
    # at E, B at E(A), P at E(Q), Q at E and E(A*); A* binds at E and E(Q). I goes
    # to {E, E(Q)} competing with A and P, to {E, E(A*)} competing with A and Q, to
    # {E(A), E(Q)} competing with B and P, to {E, E(A), E(A*)} competing with B and
    # Q, to {E, E(A), E(Q)} competing with B, P and A*, to {E, E(A)} competing with
    # B, Q and A*, and to {E} competing with A, Q and A*. I copies no reactant, so
    # every complex it forms is new: all seven children are emitted, each with full
    # rank. Mirrors: E + A ⇌ E(A) onto E(I) in A's group wherever E(A) is a site,
    # E + Q ⇌ E(Q) in Q's group wherever E(Q) is, and E + A* ⇌ E(A*) in A*'s group
    # wherever E(A*) is.
    rxn = @enzyme_reaction begin
        substrates: A[C], B[N]
        products: P[C], Q[N]
        dead_end_inhibitors: A, I
    end
    m = _testhelper_on_reaction(rxn, @enzyme_mechanism begin
        substrates: A, B
        products: P, Q
        steps: begin
            E + A ⇌ E(A)
            E(A) + B ⇌ E(A, B)
            E(A, B) <--> E(P, Q)
            E(Q) + P ⇌ E(P, Q)
            (E + Q ⇌ E(Q), E(A::Inh) + Q ⇌ E(A::Inh, Q))
            (E + A::Inh ⇌ E(A::Inh), E(Q) + A::Inh ⇌ E(A::Inh, Q))
        end
    end)
    at_E_EQ = _testhelper_on_reaction(rxn, @enzyme_mechanism begin
        substrates: A, B
        products: P, Q
        regulators: I
        steps: begin
            E + A ⇌ E(A)
            E(A) + B ⇌ E(A, B)
            E(A, B) <--> E(P, Q)
            E(Q) + P ⇌ E(P, Q)
            (E + Q ⇌ E(Q), E(A::Inh) + Q ⇌ E(A::Inh, Q), E(I) + Q ⇌ E(I, Q))
            (E + A::Inh ⇌ E(A::Inh), E(Q) + A::Inh ⇌ E(A::Inh, Q))
            (E + I ⇌ E(I), E(Q) + I ⇌ E(I, Q))
        end
    end)
    at_E_EAinh = _testhelper_on_reaction(rxn, @enzyme_mechanism begin
        substrates: A, B
        products: P, Q
        regulators: I
        steps: begin
            E + A ⇌ E(A)
            E(A) + B ⇌ E(A, B)
            E(A, B) <--> E(P, Q)
            E(Q) + P ⇌ E(P, Q)
            (E + Q ⇌ E(Q), E(A::Inh) + Q ⇌ E(A::Inh, Q))
            (E + A::Inh ⇌ E(A::Inh), E(Q) + A::Inh ⇌ E(A::Inh, Q),
             E(I) + A::Inh ⇌ E(A::Inh, I))
            (E + I ⇌ E(I), E(A::Inh) + I ⇌ E(A::Inh, I))
        end
    end)
    at_EA_EQ = _testhelper_on_reaction(rxn, @enzyme_mechanism begin
        substrates: A, B
        products: P, Q
        regulators: I
        steps: begin
            E + A ⇌ E(A)
            E(A) + B ⇌ E(A, B)
            E(A, B) <--> E(P, Q)
            E(Q) + P ⇌ E(P, Q)
            (E + Q ⇌ E(Q), E(A::Inh) + Q ⇌ E(A::Inh, Q))
            (E + A::Inh ⇌ E(A::Inh), E(Q) + A::Inh ⇌ E(A::Inh, Q))
            (E(A) + I ⇌ E(A, I), E(Q) + I ⇌ E(I, Q))
        end
    end)
    at_E_EA_EAinh = _testhelper_on_reaction(rxn, @enzyme_mechanism begin
        substrates: A, B
        products: P, Q
        regulators: I
        steps: begin
            (E + A ⇌ E(A), E(I) + A ⇌ E(A, I))
            E(A) + B ⇌ E(A, B)
            E(A, B) <--> E(P, Q)
            E(Q) + P ⇌ E(P, Q)
            (E + Q ⇌ E(Q), E(A::Inh) + Q ⇌ E(A::Inh, Q))
            (E + A::Inh ⇌ E(A::Inh), E(Q) + A::Inh ⇌ E(A::Inh, Q),
             E(I) + A::Inh ⇌ E(A::Inh, I))
            (E + I ⇌ E(I), E(A) + I ⇌ E(A, I), E(A::Inh) + I ⇌ E(A::Inh, I))
        end
    end)
    at_E_EA_EQ = _testhelper_on_reaction(rxn, @enzyme_mechanism begin
        substrates: A, B
        products: P, Q
        regulators: I
        steps: begin
            (E + A ⇌ E(A), E(I) + A ⇌ E(A, I))
            E(A) + B ⇌ E(A, B)
            E(A, B) <--> E(P, Q)
            E(Q) + P ⇌ E(P, Q)
            (E + Q ⇌ E(Q), E(A::Inh) + Q ⇌ E(A::Inh, Q), E(I) + Q ⇌ E(I, Q))
            (E + A::Inh ⇌ E(A::Inh), E(Q) + A::Inh ⇌ E(A::Inh, Q))
            (E + I ⇌ E(I), E(A) + I ⇌ E(A, I), E(Q) + I ⇌ E(I, Q))
        end
    end)
    at_E_EA = _testhelper_on_reaction(rxn, @enzyme_mechanism begin
        substrates: A, B
        products: P, Q
        regulators: I
        steps: begin
            (E + A ⇌ E(A), E(I) + A ⇌ E(A, I))
            E(A) + B ⇌ E(A, B)
            E(A, B) <--> E(P, Q)
            E(Q) + P ⇌ E(P, Q)
            (E + Q ⇌ E(Q), E(A::Inh) + Q ⇌ E(A::Inh, Q))
            (E + A::Inh ⇌ E(A::Inh), E(Q) + A::Inh ⇌ E(A::Inh, Q))
            (E + I ⇌ E(I), E(A) + I ⇌ E(A, I))
        end
    end)
    at_E = _testhelper_on_reaction(rxn, @enzyme_mechanism begin
        substrates: A, B
        products: P, Q
        regulators: I
        steps: begin
            E + A ⇌ E(A)
            E(A) + B ⇌ E(A, B)
            E(A, B) <--> E(P, Q)
            E(Q) + P ⇌ E(P, Q)
            (E + Q ⇌ E(Q), E(A::Inh) + Q ⇌ E(A::Inh, Q))
            (E + A::Inh ⇌ E(A::Inh), E(Q) + A::Inh ⇌ E(A::Inh, Q))
            E + I ⇌ E(I)
        end
    end)
    kids = EnzymeRates._expand_add_dead_end_regulator(m, rxn)
    @test length(kids) == 7
    @test Set(kids) == Set([at_E_EQ, at_E_EAinh, at_EA_EQ, at_E_EA_EAinh, at_E_EA_EQ,
                            at_E_EA, at_E])
    for k in kids
        @test EnzymeRates._assert_emission_rules(k) === nothing
        @test _testhelper_fitted(k) == _testhelper_fitted(m) + 1
        @test _testhelper_identifiable_rank(k) == _testhelper_fitted(k)
    end
end

@testset "Mechanism — a form that holds only a copy has room for an inhibitor" begin
    # Uni-bi with S as its own inhibitor at {E, E(Q)}. The capacity test counts
    # productive bindings only: E(S*) holds S's copy at a dead-end site of its own
    # and no substrate, so it is a site although S is the only substrate. S binds
    # productively at E, P at E(Q), Q at E and E(S*); S* binds at E and E(Q). A
    # foreign inhibitor I goes to {E, E(Q)} competing with S and P, to {E, E(S*)}
    # competing with S and Q, and to {E} competing with S, Q and S*. The child at
    # {E, E(Q)} mirrors E + Q ⇌ E(Q) onto E(I) in Q's group, the one at {E, E(S*)}
    # mirrors E + S* ⇌ E(S*) onto E(I) in S*'s group. Each child has full rank.
    rxn = @enzyme_reaction begin
        substrates: S[CN]
        products: P[C], Q[N]
        dead_end_inhibitors: S, I
    end
    m = _testhelper_on_reaction(rxn, @enzyme_mechanism begin
        substrates: S
        products: P, Q
        steps: begin
            E + S ⇌ E(S)
            E(S) <--> E(P, Q)
            E(Q) + P ⇌ E(P, Q)
            (E + Q ⇌ E(Q), E(S::Inh) + Q ⇌ E(Q, S::Inh))
            (E + S::Inh ⇌ E(S::Inh), E(Q) + S::Inh ⇌ E(Q, S::Inh))
        end
    end)
    at_E_EQ = _testhelper_on_reaction(rxn, @enzyme_mechanism begin
        substrates: S
        products: P, Q
        regulators: I
        steps: begin
            E + S ⇌ E(S)
            E(S) <--> E(P, Q)
            E(Q) + P ⇌ E(P, Q)
            (E + Q ⇌ E(Q), E(S::Inh) + Q ⇌ E(Q, S::Inh), E(I) + Q ⇌ E(I, Q))
            (E + S::Inh ⇌ E(S::Inh), E(Q) + S::Inh ⇌ E(Q, S::Inh))
            (E + I ⇌ E(I), E(Q) + I ⇌ E(I, Q))
        end
    end)
    at_E_ESinh = _testhelper_on_reaction(rxn, @enzyme_mechanism begin
        substrates: S
        products: P, Q
        regulators: I
        steps: begin
            E + S ⇌ E(S)
            E(S) <--> E(P, Q)
            E(Q) + P ⇌ E(P, Q)
            (E + Q ⇌ E(Q), E(S::Inh) + Q ⇌ E(Q, S::Inh))
            (E + S::Inh ⇌ E(S::Inh), E(Q) + S::Inh ⇌ E(Q, S::Inh),
             E(I) + S::Inh ⇌ E(I, S::Inh))
            (E + I ⇌ E(I), E(S::Inh) + I ⇌ E(I, S::Inh))
        end
    end)
    at_E = _testhelper_on_reaction(rxn, @enzyme_mechanism begin
        substrates: S
        products: P, Q
        regulators: I
        steps: begin
            E + S ⇌ E(S)
            E(S) <--> E(P, Q)
            E(Q) + P ⇌ E(P, Q)
            (E + Q ⇌ E(Q), E(S::Inh) + Q ⇌ E(Q, S::Inh))
            (E + S::Inh ⇌ E(S::Inh), E(Q) + S::Inh ⇌ E(Q, S::Inh))
            E + I ⇌ E(I)
        end
    end)
    kids = EnzymeRates._expand_add_dead_end_regulator(m, rxn)
    @test length(kids) == 3
    @test Set(kids) == Set([at_E_EQ, at_E_ESinh, at_E])
    for k in kids
        @test EnzymeRates._assert_emission_rules(k) === nothing
        @test _testhelper_identifiable_rank(k) == _testhelper_fitted(k)
    end
end

@testset "Mechanism — competing with B targets where B binds, not where B* binds" begin
    # Ordered bi-bi with B as its own inhibitor at E (placed competing with A and Q),
    # placing a foreign inhibitor I. Competition with a reactant targets the forms
    # where it binds productively: A at E, B at E(A), P at E(Q), Q at E; the copy B*
    # binds at E, a target only for competition with B* itself. Capacity leaves out
    # E(A, B) and E(P, Q); E(B*) is never a target, since nothing binds there. I goes
    # to {E, E(Q)} competing with A and P, to {E} competing with A and Q, to
    # {E(A), E(Q)} competing with B and P, to {E, E(A)} competing with B and Q, and to
    # {E, E(A), E(Q)} competing with B, P and B*. In the pattern B and P, free E is a
    # site of B's copy, holds nothing that competes and is not a site of P: targeting
    # by name, B and its copy counting as one competitor, would add it and turn
    # {E(A), E(Q)} into {E, E(A), E(Q)}, so {E(A), E(Q)} comes only from per-site
    # targeting. Mirrors: E + Q ⇌ E(Q) onto E(I) in Q's group wherever E and E(Q) are
    # sites, E + A ⇌ E(A) in A's group wherever E and E(A) are. I is foreign, so each
    # child adds one constant and has full rank (7 fitted, rank 7).
    rxn = @enzyme_reaction begin
        substrates: A[C], B[N]
        products: P[C], Q[N]
        dead_end_inhibitors: B, I
    end
    m = _testhelper_on_reaction(rxn, @enzyme_mechanism begin
        substrates: A, B
        products: P, Q
        steps: begin
            E + A ⇌ E(A)
            E(A) + B ⇌ E(A, B)
            E(A, B) <--> E(P, Q)
            E(Q) + P ⇌ E(P, Q)
            E + Q ⇌ E(Q)
            E + B::Inh ⇌ E(B::Inh)
        end
    end)
    at_E_EQ = _testhelper_on_reaction(rxn, @enzyme_mechanism begin
        substrates: A, B
        products: P, Q
        regulators: I
        steps: begin
            E + A ⇌ E(A)
            E(A) + B ⇌ E(A, B)
            E(A, B) <--> E(P, Q)
            E(Q) + P ⇌ E(P, Q)
            (E + Q ⇌ E(Q), E(I) + Q ⇌ E(I, Q))
            E + B::Inh ⇌ E(B::Inh)
            (E + I ⇌ E(I), E(Q) + I ⇌ E(I, Q))
        end
    end)
    at_E = _testhelper_on_reaction(rxn, @enzyme_mechanism begin
        substrates: A, B
        products: P, Q
        regulators: I
        steps: begin
            E + A ⇌ E(A)
            E(A) + B ⇌ E(A, B)
            E(A, B) <--> E(P, Q)
            E(Q) + P ⇌ E(P, Q)
            E + Q ⇌ E(Q)
            E + B::Inh ⇌ E(B::Inh)
            E + I ⇌ E(I)
        end
    end)
    at_EA_EQ = _testhelper_on_reaction(rxn, @enzyme_mechanism begin
        substrates: A, B
        products: P, Q
        regulators: I
        steps: begin
            E + A ⇌ E(A)
            E(A) + B ⇌ E(A, B)
            E(A, B) <--> E(P, Q)
            E(Q) + P ⇌ E(P, Q)
            E + Q ⇌ E(Q)
            E + B::Inh ⇌ E(B::Inh)
            (E(A) + I ⇌ E(A, I), E(Q) + I ⇌ E(I, Q))
        end
    end)
    at_E_EA = _testhelper_on_reaction(rxn, @enzyme_mechanism begin
        substrates: A, B
        products: P, Q
        regulators: I
        steps: begin
            (E + A ⇌ E(A), E(I) + A ⇌ E(A, I))
            E(A) + B ⇌ E(A, B)
            E(A, B) <--> E(P, Q)
            E(Q) + P ⇌ E(P, Q)
            E + Q ⇌ E(Q)
            E + B::Inh ⇌ E(B::Inh)
            (E + I ⇌ E(I), E(A) + I ⇌ E(A, I))
        end
    end)
    at_E_EA_EQ = _testhelper_on_reaction(rxn, @enzyme_mechanism begin
        substrates: A, B
        products: P, Q
        regulators: I
        steps: begin
            (E + A ⇌ E(A), E(I) + A ⇌ E(A, I))
            E(A) + B ⇌ E(A, B)
            E(A, B) <--> E(P, Q)
            E(Q) + P ⇌ E(P, Q)
            (E + Q ⇌ E(Q), E(I) + Q ⇌ E(I, Q))
            E + B::Inh ⇌ E(B::Inh)
            (E + I ⇌ E(I), E(A) + I ⇌ E(A, I), E(Q) + I ⇌ E(I, Q))
        end
    end)
    kids = EnzymeRates._expand_add_dead_end_regulator(m, rxn)
    @test length(kids) == 5
    @test Set(kids) == Set([at_E_EQ, at_E, at_EA_EQ, at_E_EA, at_E_EA_EQ])
    @test _testhelper_identifiable_rank(m) == _testhelper_fitted(m) == 6
    for k in kids
        @test EnzymeRates._assert_emission_rules(k) === nothing
        @test _testhelper_identifiable_rank(k) == _testhelper_fitted(k) == 7
    end
end

@testset "AllostericMechanism — a copy that is new only in the inactive state is kept" begin
    # Uni-uni whose S binding and chemistry are `:OnlyA`: the inactive conformation
    # binds P but not S. An `:EqualAI` copy of S at E duplicates E(S) in the active
    # state and is the only S-bound form in the inactive one, where its constant is
    # visible; the placement stands and adds an identifiable constant. With S
    # binding `:NonequalAI` the copy duplicates E(S) in both states, the S group
    # alone forms E(S) in each, so the gauge exists in each, and it is skipped.
    rxn = @enzyme_reaction begin
        substrates: S[C]
        products: P[C]
        dead_end_inhibitors: S
        oligomeric_state: 2
    end
    onlya = _testhelper_on_reaction_mult2(rxn, @allosteric_mechanism begin
        substrates: S
        products: P
        catalytic_multiplicity: 2
        catalytic_steps: begin
            E + S ⇌ E(S)      :: OnlyA
            E(S) <--> E(P)    :: OnlyA
            E + P ⇌ E(P)      :: EqualAI
        end
    end)
    kept = _testhelper_on_reaction_mult2(rxn, @allosteric_mechanism begin
        substrates: S
        products: P
        catalytic_inhibitors: S
        catalytic_multiplicity: 2
        catalytic_steps: begin
            E + S ⇌ E(S)                :: OnlyA
            E(S) <--> E(P)              :: OnlyA
            E + P ⇌ E(P)                :: EqualAI
            E + S::Inh ⇌ E(S::Inh)      :: EqualAI
        end
    end)
    kids = EnzymeRates._expand_add_dead_end_regulator(onlya, rxn)
    @test length(kids) == 1
    @test Set(kids) == Set([kept])
    @test _testhelper_identifiable_rank(kept) == _testhelper_identifiable_rank(onlya) + 1
    @test isempty(EnzymeRates._redundant_copy_groups(kept))
    nonequal = _testhelper_on_reaction_mult2(rxn, @allosteric_mechanism begin
        substrates: S
        products: P
        catalytic_multiplicity: 2
        catalytic_steps: begin
            E + S ⇌ E(S)      :: NonequalAI
            E(S) <--> E(P)    :: OnlyA
            E + P ⇌ E(P)      :: EqualAI
        end
    end)
    @test isempty(EnzymeRates._expand_add_dead_end_regulator(nonequal, rxn))
end

@testset "AllostericMechanism — free enzyme binds a copy in every conformation" begin
    # Uni-uni whose every catalytic group is `:OnlyA`: the inactive conformation has no
    # step, yet it holds its free enzyme. An `:EqualAI` copy of S at E duplicates E(S) in
    # the active state and is the only S-bound form in the inactive one, where its
    # constant is visible; the placement stands and adds an identifiable constant.
    rxn = @enzyme_reaction begin
        substrates: S[C]
        products: P[C]
        dead_end_inhibitors: S
        oligomeric_state: 2
    end
    dead_inactive = _testhelper_on_reaction_mult2(rxn, @allosteric_mechanism begin
        substrates: S
        products: P
        catalytic_multiplicity: 2
        catalytic_steps: begin
            E + S ⇌ E(S)      :: OnlyA
            E(S) <--> E(P)    :: OnlyA
            E + P ⇌ E(P)      :: OnlyA
        end
    end)
    kept = _testhelper_on_reaction_mult2(rxn, @allosteric_mechanism begin
        substrates: S
        products: P
        catalytic_inhibitors: S
        catalytic_multiplicity: 2
        catalytic_steps: begin
            E + S ⇌ E(S)                :: OnlyA
            E(S) <--> E(P)              :: OnlyA
            E + P ⇌ E(P)                :: OnlyA
            E + S::Inh ⇌ E(S::Inh)      :: EqualAI
        end
    end)
    kids = EnzymeRates._expand_add_dead_end_regulator(dead_inactive, rxn)
    @test length(kids) == 1
    @test Set(kids) == Set([kept])
    @test isempty(EnzymeRates._redundant_copy_groups(kept))
    @test _testhelper_identifiable_rank(kept) ==
        _testhelper_identifiable_rank(dead_inactive) + 1
end

@testset "AllostericMechanism — no relaxation makes a kept copy redundant" begin
    # The copy of S at E is new only in the inactive state, where S does not bind.
    # Relaxing the S binding to `:NonequalAI` brings E(S) into the inactive state:
    # the copy then duplicates E(S) in both states, formed in each by the S group
    # alone, so the gauge exists in each, its constant shows only beside K_S, and
    # the child is not emitted (one fitted constant above its rank).
    # Relaxing the chemistry while S binds the active state only is rejected by
    # `_onlya_haldane_violation` first, since the Haldane row then has one `:OnlyA`
    # column, and partial catalysis would reject it too: the inactive state would
    # catalyze in part. The P binding and the copy relax. The relaxed copy, with one
    # constant per state, duplicates E(S) in the active state and is new in the
    # inactive one, so it keeps one phantom, its active-state constant; the copy rule
    # keeps such a `:NonequalAI` copy.
    rxn = @enzyme_reaction begin
        substrates: S[C]
        products: P[C]
        dead_end_inhibitors: S
        oligomeric_state: 2
    end
    kept = _testhelper_on_reaction_mult2(rxn, @allosteric_mechanism begin
        substrates: S
        products: P
        catalytic_inhibitors: S
        catalytic_multiplicity: 2
        catalytic_steps: begin
            E + S ⇌ E(S)                :: OnlyA
            E(S) <--> E(P)              :: OnlyA
            E + P ⇌ E(P)                :: EqualAI
            E + S::Inh ⇌ E(S::Inh)      :: EqualAI
        end
    end)
    relaxed_p = _testhelper_on_reaction_mult2(rxn, @allosteric_mechanism begin
        substrates: S
        products: P
        catalytic_inhibitors: S
        catalytic_multiplicity: 2
        catalytic_steps: begin
            E + S ⇌ E(S)                :: OnlyA
            E(S) <--> E(P)              :: OnlyA
            E + P ⇌ E(P)                :: NonequalAI
            E + S::Inh ⇌ E(S::Inh)      :: EqualAI
        end
    end)
    relaxed_copy = _testhelper_on_reaction_mult2(rxn, @allosteric_mechanism begin
        substrates: S
        products: P
        catalytic_inhibitors: S
        catalytic_multiplicity: 2
        catalytic_steps: begin
            E + S ⇌ E(S)                :: OnlyA
            E(S) <--> E(P)              :: OnlyA
            E + P ⇌ E(P)                :: EqualAI
            E + S::Inh ⇌ E(S::Inh)      :: NonequalAI
        end
    end)
    relaxed_s = _testhelper_on_reaction_mult2(rxn, @allosteric_mechanism begin
        substrates: S
        products: P
        catalytic_inhibitors: S
        catalytic_multiplicity: 2
        catalytic_steps: begin
            E + S ⇌ E(S)                :: NonequalAI
            E(S) <--> E(P)              :: OnlyA
            E + P ⇌ E(P)                :: EqualAI
            E + S::Inh ⇌ E(S::Inh)      :: EqualAI
        end
    end)
    kids = EnzymeRates._expand_change_allo_state(kept)
    @test length(kids) == 2
    @test Set(kids) == Set([relaxed_p, relaxed_copy])
    @test !isempty(EnzymeRates._redundant_copy_groups(relaxed_s))
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
        dead_end_inhibitors: S
    end
    m = EnzymeRates.Mechanism(@enzyme_mechanism begin
        substrates: S
        products: P
        steps: begin
            E ⇌ Estar
            Estar + S ⇌ Estar(S)
            Estar(S) <--> E(P)
            E + P ⇌ E(P)
        end
    end)
    @test isempty(EnzymeRates._expand_add_dead_end_regulator(m, rxn))
    absent = EnzymeRates.Mechanism(@enzyme_mechanism begin
        substrates: S
        products: P
        steps: begin
            (E ⇌ Estar, E(S::Inh) ⇌ Estar(S::Inh))
            Estar + S ⇌ Estar(S)
            Estar(S) <--> E(P)
            E + P ⇌ E(P)
            (E + S::Inh ⇌ E(S::Inh), Estar + S::Inh ⇌ Estar(S::Inh))
        end
    end)
    @test _testhelper_identifiable_rank(absent) == _testhelper_identifiable_rank(m)
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
    _testhelper_assert_mechanism_invariants(am)
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
    base_fitted = _testhelper_fitted(am)
    for r in result
        r_fitted = _testhelper_fitted(r)
        @test r_fitted == base_fitted + 1
    end

    # 3. compilability
    for r in result
        @test r isa EnzymeRates.AllostericMechanism
        _testhelper_assert_mechanism_invariants(r)
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

@testset "Mechanism — no regulators: empty (negative)" begin
    m = first(EnzymeRates.init_mechanisms(uni_uni_rxn))
    rxn = @enzyme_reaction begin
        substrates: S[C]
        products: P[C]
    end
    @test isempty(EnzymeRates._expand_add_dead_end_regulator(m, rxn))
end

end

# ═══════════════════════════════════════════════════════════════════════
# 4. Allosteric expansion moves
# ═══════════════════════════════════════════════════════════════════════

# ─── _expand_to_allosteric ─────────────────────────────────────────────
@testset "_expand_to_allosteric" begin

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
    _testhelper_assert_mechanism_invariants(m)
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
    @test _testhelper_param_deltas(m, result) == fill(1, 15)

    # 3. invariants (the deltas above compile every result), and no result is the
    # all-:EqualAI baseline.
    for r in result
        _testhelper_assert_mechanism_invariants(r)
        @test !all(==(:EqualAI), EnzymeRates.cat_allo_states(r))
    end
end

@testset "Mechanism — Bi-bi ping-pong: binding-group :OnlyA variants only" begin
    # SEED: bi-bi ping-pong topology, 6 kinetic groups: 4 bindings
    # (A, B substrate-side; P, Q product-side) + 2 chemical steps (the
    # two half-reaction interconversions E(A)<-->Estar(P) and
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
            Estar + P ⇌ Estar(P)
            E(A) <--> Estar(P)
            Estar(B) ⇌ E(Q)
        end
    end
    m = EnzymeRates.Mechanism(em_seed)
    _testhelper_assert_mechanism_invariants(m)
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
    @test _testhelper_param_deltas(m, result) == fill(1, 15)

    # 3. invariants (the deltas above compile every result), and no result is the
    # all-:EqualAI baseline.
    for r in result
        _testhelper_assert_mechanism_invariants(r)
        @test !all(==(:EqualAI), EnzymeRates.cat_allo_states(r))
    end
end

@testset "Mechanism — uni-uni: binding-group :OnlyA variants only" begin
    # SEED: uni-uni, no declared regulator. 3 kinetic groups: 2 bindings (S, P) + 1
    # catalytic. The catalytic group is never a bare primary :OnlyA (that V-type needs a
    # regulator, none declared), so every child starts from a binding promotion, which
    # alone leaves the Haldane unsatisfiable. Each binding's one-sided :OnlyA is closed
    # over its minimal Haldane completion: every non-empty subset of the 2 bindings
    # (S, P) is set :OnlyA, each with the chemical step also :OnlyA (dead inactive).
    m = _testhelper_on_reaction(uni_uni_allo, @enzyme_mechanism begin
        substrates: S
        products: P
        steps: begin
            E + S ⇌ E(S)
            E(S) <--> E(P)
            E + P ⇌ E(P)
        end
    end)
    only_s = _testhelper_on_reaction_mult2(uni_uni_allo, @allosteric_mechanism begin
        substrates: S
        products: P
        catalytic_multiplicity: 2
        catalytic_steps: begin
            E + S ⇌ E(S)      :: OnlyA
            E(S) <--> E(P)    :: OnlyA
            E + P ⇌ E(P)      :: EqualAI
        end
    end)
    only_p = _testhelper_on_reaction_mult2(uni_uni_allo, @allosteric_mechanism begin
        substrates: S
        products: P
        catalytic_multiplicity: 2
        catalytic_steps: begin
            E + S ⇌ E(S)      :: EqualAI
            E(S) <--> E(P)    :: OnlyA
            E + P ⇌ E(P)      :: OnlyA
        end
    end)
    both = _testhelper_on_reaction_mult2(uni_uni_allo, @allosteric_mechanism begin
        substrates: S
        products: P
        catalytic_multiplicity: 2
        catalytic_steps: begin
            E + S ⇌ E(S)      :: OnlyA
            E(S) <--> E(P)    :: OnlyA
            E + P ⇌ E(P)      :: OnlyA
        end
    end)

    result = EnzymeRates._expand_to_allosteric(m, uni_uni_allo)
    @test length(result) == 3
    @test Set(result) == Set([only_s, only_p, both])

    for r in result
        _testhelper_assert_mechanism_invariants(r)
        @test EnzymeRates.compile_mechanism(r) isa EnzymeRates.AllostericEnzymeMechanism
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
    base = _testhelper_fitted(seed)
    Δ(am) = _testhelper_fitted(am) - base
    outs = EnzymeRates._expand_to_allosteric(seed, uni_uni_allo_2reg)
    vtypes = filter(am -> !isempty(EnzymeRates.regulatory_sites(am)), outs)
    ktypes = filter(am -> isempty(EnzymeRates.regulatory_sites(am)), outs)

    @test length(vtypes) == 4
    @test all(am -> Δ(am) == 2, vtypes)
    # Every V-type has its catalytic group :OnlyA.
    @test all(vtypes) do am
        cat_group = only(g for g in eachindex(EnzymeRates.steps(am))
                         if _is_catalytic_group(am, g))
        EnzymeRates.cat_allo_states(am)[cat_group] == :OnlyA
    end
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

@testset "Mechanism — non-hyperbolic catalytic scheme: multiplicity 1 only" begin
    # The same uni-bi with random SS product release, with multiplicities 1 and 2
    # allowed. One catalytic subunit adds no concentration power of its own, so
    # the multiplicity-1 variants are emitted exactly as for a reaction allowing
    # only multiplicity 1; multiplicity 2 would stack a second source of powers
    # and emits nothing. Three binding groups (S, P, Q), each subset :OnlyA with
    # the chemistry :OnlyA: 2^3 - 1 = 7 K-type children, none rejected by
    # `_onlya_haldane_violation`.
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
    rxn_1_2 = @enzyme_reaction begin
        substrates: S[AB]
        products: P[A], Q[B]
        allowed_catalytic_multiplicities: (1, 2)
    end
    rxn_1 = @enzyme_reaction begin
        substrates: S[AB]
        products: P[A], Q[B]
        allowed_catalytic_multiplicities: (1,)
    end
    children = EnzymeRates._expand_to_allosteric(m, rxn_1_2)
    @test length(children) == 7
    @test all(am -> EnzymeRates.catalytic_multiplicity(am) == 1, children)
    @test Set(children) == Set(EnzymeRates._expand_to_allosteric(m, rxn_1))
end

@testset "Mechanism — catalytic-site abortive complex: no children" begin
    # Ordered SS bi-bi with A binding E(Q) as an abortive complex at its
    # catalytic site: the term rooted at {E(Q), E(A, Q)} holds E → E(A), so the
    # scheme carries A², and a conformational equilibrium on top would stack a
    # second source of concentration powers. The promotion emits nothing.
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
    @test isempty(EnzymeRates._expand_to_allosteric(m, rxn))
end

@testset "Mechanism — substrate declared as an inhibitor: promotion proceeds" begin
    # Ordered SS bi-bi with A also declared as a dead-end inhibitor, written
    # with the ::Inh role tag (src/dsl.jl): a plain `A` occurrence keeps A's
    # declared substrate role, while `A::Inh` binds CompetitiveInhibitor(:A),
    # a distinct species from the catalytic A-bound forms. The inhibitor's
    # own group binds E and E(Q); Q's catalytic binding is mirrored onto the
    # inhibitor-bound form E(A::Inh), sharing its kinetic group, exactly as
    # `_expand_add_dead_end_regulator` builds it.
    parent = EnzymeRates.Mechanism(@enzyme_mechanism begin
        substrates: A, B
        products: P, Q
        steps: begin
            E + A <--> E(A)
            E(A) + B <--> E(A, B)
            (E + Q <--> E(Q), E(A::Inh) + Q <--> E(A::Inh, Q))
            E(Q) + P <--> E(P, Q)
            E(A, B) <--> E(P, Q)
            (E + A::Inh ⇌ E(A::Inh), E(Q) + A::Inh ⇌ E(A::Inh, Q))
        end
    end)
    rxn = @enzyme_reaction begin
        substrates: A[C], B[N]
        products: P[C], Q[N]
        dead_end_inhibitors: A
        oligomeric_state: 2
    end

    # Five binding groups (A, B, Q with its inhibitor-bound mirror, P, and the
    # inhibitor's own group), each subset :OnlyA with the chemistry :OnlyA:
    # 2^5 - 1 = 31 candidate K-type children; _onlya_haldane_violation drops
    # none of them. The reaction declares no allosteric regulator, so no
    # V-type variant is emitted. N = 31. The derived denominator carries A²,
    # but from the inhibitor site, so the catalytic scheme is hyperbolic and
    # the promotion proceeds.
    N = 31
    children = EnzymeRates._expand_to_allosteric(parent, rxn)
    @test length(children) == N
    @test all(EnzymeRates._hyperbolic_catalysis, children)
end
end

# ─── _expand_add_allosteric_regulator ──────────────────────────────────
@testset "_expand_add_allosteric_regulator" begin

@testset "AllostericMechanism — a first allo regulator: 3 variants" begin
    # SEED: uni-uni allosteric with all groups :EqualAI and no allosteric regulator
    # added yet. A declared allosteric regulator x is the only un-added one, so it gets
    # one new site (0 existing reg sites) with each of the 3 non-:EqualAI tags. The
    # :EqualAI branch is gated to "existing site with ≥1 non-:EqualAI ligand", so it
    # does not apply.
    am = EnzymeRates.AllostericMechanism(@allosteric_mechanism begin
        substrates: S; products: P
        catalytic_multiplicity: 2
        catalytic_steps: begin
            E + P ⇌ E(P)    :: EqualAI
            E + S ⇌ E(S)    :: EqualAI
            E(S) <--> E(P)  :: EqualAI
        end
    end)
    _testhelper_assert_mechanism_invariants(am)

    for (rxn, x) in [
        # R is declared only as an allosteric regulator.
        (@enzyme_reaction(begin
            substrates: S[C]
            products: P[C]
            allosteric_regulators: R
            oligomeric_state: 2
        end), :R),
        # S is both the substrate and an allosteric regulator; it still plays its
        # catalytic substrate role in the base mechanism.
        (@enzyme_reaction(begin
            substrates: S[C]
            products: P[C]
            allosteric_regulators: S
            oligomeric_state: 2
        end), :S),
        # P is both the product and an allosteric regulator.
        (@enzyme_reaction(begin
            substrates: S[C]
            products: P[C]
            allosteric_regulators: P
            oligomeric_state: 2
        end), :P),
    ]
        result = EnzymeRates._expand_add_allosteric_regulator(am, rxn)

        # 1. count: 3 non-:EqualAI tags × 1 site option (new site only) = 3.
        @test length(result) == 3

        # 2. Δ params: :OnlyA/:OnlyI add one K_x each (+1); :NonequalAI adds K_x and
        # K_x_T (+2). Sorted [1, 1, 2]. Measured against the actual compiled fitted
        # count.
        @test _testhelper_param_deltas(am, result) == [1, 1, 2]

        # 3. each variant has exactly one regulatory site holding the single ligand x,
        # and the tags across the 3 variants are exactly {:OnlyA, :OnlyI, :NonequalAI}.
        # The added site's multiplicity is am.catalytic_multiplicity.
        @test all(r -> length(r.regulatory_sites) == 1, result)
        @test Set((Symbol[EnzymeRates.name(l) for l in EnzymeRates.ligands(site)],
                   collect(EnzymeRates.allo_states(site)),
                   EnzymeRates.multiplicity(site))
                  for r in result for site in r.regulatory_sites) == Set([
            ([x], [:OnlyA],      am.catalytic_multiplicity),
            ([x], [:OnlyI],      am.catalytic_multiplicity),
            ([x], [:NonequalAI], am.catalytic_multiplicity),
        ])

        # 4. invariants and preservation: catalytic side and cat_allo_states untouched.
        for r in result
            _testhelper_assert_mechanism_invariants(r)
            @test r.catalytic_multiplicity == am.catalytic_multiplicity
            @test r.cat_allo_states == am.cat_allo_states
            @test EnzymeRates.reaction(r) == EnzymeRates.reaction(am)
        end
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
    _testhelper_assert_mechanism_invariants(am)
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
    _testhelper_assert_mechanism_invariants(am)

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
    @test _testhelper_param_deltas(am, result) == [1, 1, 1, 1, 2, 2]

    # 3. structural: every result has :R2 in some regulatory site.
    for r in result
        has_r2 = any(r.regulatory_sites) do site
            any(l -> EnzymeRates.name(l) === :R2,
                EnzymeRates.ligands(site))
        end
        @test has_r2
    end

    # 4. structural: at least one result has :R2 :EqualAI at site 1 (the same site as
    # R1). This is the :EqualAI-at-existing branch, gated on R1 (the existing ligand)
    # being non-:EqualAI.
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

    # 5. property-style: separate new-site (#sites grows to 2) and existing-site
    # (#sites stays at 1) placements.
    new_site_variants = filter(r -> length(r.regulatory_sites) == 2, result)
    existing_site_variants = filter(r -> length(r.regulatory_sites) == 1, result)
    @test length(new_site_variants) == 3   # 3 non-:EqualAI × new site
    # :OnlyA + :NonequalAI appends onto R1's :OnlyA site + 1 :EqualAI; the :OnlyI
    # append is skipped as redundant (disjoint conformations).
    @test length(existing_site_variants) == 3

    # 6. invariants and preservation: catalytic side and cat_allo_states untouched.
    for r in result
        _testhelper_assert_mechanism_invariants(r)
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
    _testhelper_assert_mechanism_invariants(am)

    result = EnzymeRates._expand_add_allosteric_regulator(
        am, uni_uni_allo_2reg)

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
    _testhelper_assert_mechanism_invariants(am)
    @test isempty(EnzymeRates._expand_add_allosteric_regulator(am, rxn))
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
        _testhelper_assert_mechanism_invariants(r)
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
    _testhelper_assert_mechanism_invariants(am)

    fp = collect(EnzymeRates.fitted_params(am))
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
    am_x = EnzymeRates.AllostericMechanism(@allosteric_mechanism begin
        substrates: S; products: P
        allosteric_regulators: X::OnlyA
        catalytic_multiplicity: 2
        catalytic_steps: begin
            E + P ⇌ E(P)    :: EqualAI
            E + S ⇌ E(S)    :: EqualAI
            E(S) <--> E(P)  :: EqualAI
        end
    end)
    am_y = EnzymeRates.AllostericMechanism(@allosteric_mechanism begin
        substrates: S; products: P
        allosteric_regulators: Y::NonequalAI
        catalytic_multiplicity: 2
        catalytic_steps: begin
            E + P ⇌ E(P)    :: EqualAI
            E + S ⇌ E(S)    :: EqualAI
            E(S) <--> E(P)  :: EqualAI
        end
    end)
    # The site both routes reach: X::OnlyA and Y::NonequalAI share one site.
    shared = EnzymeRates.AllostericMechanism(@allosteric_mechanism begin
        substrates: S; products: P
        allosteric_regulators: X::OnlyA, Y::NonequalAI
        catalytic_multiplicity: 2
        catalytic_steps: begin
            E + P ⇌ E(P)    :: EqualAI
            E + S ⇌ E(S)    :: EqualAI
            E(S) <--> E(P)  :: EqualAI
        end
        regulatory_site(multiplicity = 2): begin
            ligands: X, Y
        end
    end)
    # X::OnlyA and Y::NonequalAI on sites of their own: also reached from either parent.
    x_y_apart = EnzymeRates.AllostericMechanism(@allosteric_mechanism begin
        substrates: S; products: P
        allosteric_regulators: X::OnlyA, Y::NonequalAI
        catalytic_multiplicity: 2
        catalytic_steps: begin
            E + P ⇌ E(P)    :: EqualAI
            E + S ⇌ E(S)    :: EqualAI
            E(S) <--> E(P)  :: EqualAI
        end
    end)

    # From am_x, Y joins as :OnlyA or :NonequalAI on a new site or X's site, as :OnlyI on a
    # new site only (an :OnlyI ligand on X's all-:OnlyA site acts on a disjoint
    # conformation), and as :EqualAI on X's site.
    y_onlya_new = EnzymeRates.AllostericMechanism(@allosteric_mechanism begin
        substrates: S; products: P
        allosteric_regulators: X::OnlyA, Y::OnlyA
        catalytic_multiplicity: 2
        catalytic_steps: begin
            E + P ⇌ E(P)    :: EqualAI
            E + S ⇌ E(S)    :: EqualAI
            E(S) <--> E(P)  :: EqualAI
        end
    end)
    y_onlya_at_x = EnzymeRates.AllostericMechanism(@allosteric_mechanism begin
        substrates: S; products: P
        allosteric_regulators: X::OnlyA, Y::OnlyA
        catalytic_multiplicity: 2
        catalytic_steps: begin
            E + P ⇌ E(P)    :: EqualAI
            E + S ⇌ E(S)    :: EqualAI
            E(S) <--> E(P)  :: EqualAI
        end
        regulatory_site(multiplicity = 2): begin
            ligands: X, Y
        end
    end)
    y_onlyi_new = EnzymeRates.AllostericMechanism(@allosteric_mechanism begin
        substrates: S; products: P
        allosteric_regulators: X::OnlyA, Y::OnlyI
        catalytic_multiplicity: 2
        catalytic_steps: begin
            E + P ⇌ E(P)    :: EqualAI
            E + S ⇌ E(S)    :: EqualAI
            E(S) <--> E(P)  :: EqualAI
        end
    end)
    y_equalai_at_x = EnzymeRates.AllostericMechanism(@allosteric_mechanism begin
        substrates: S; products: P
        allosteric_regulators: X::OnlyA, Y::EqualAI
        catalytic_multiplicity: 2
        catalytic_steps: begin
            E + P ⇌ E(P)    :: EqualAI
            E + S ⇌ E(S)    :: EqualAI
            E(S) <--> E(P)  :: EqualAI
        end
        regulatory_site(multiplicity = 2): begin
            ligands: X, Y
        end
    end)
    from_x = EnzymeRates._expand_add_allosteric_regulator(am_x, rxn)
    @test length(from_x) == 6
    @test Set(from_x) == Set([y_onlya_new, y_onlya_at_x, y_onlyi_new, x_y_apart, shared,
                              y_equalai_at_x])

    # From am_y, X joins with each non-:EqualAI tag on a new site or Y's site (Y's
    # :NonequalAI site acts on both conformations), and as :EqualAI on Y's site.
    x_onlyi_new = EnzymeRates.AllostericMechanism(@allosteric_mechanism begin
        substrates: S; products: P
        allosteric_regulators: X::OnlyI, Y::NonequalAI
        catalytic_multiplicity: 2
        catalytic_steps: begin
            E + P ⇌ E(P)    :: EqualAI
            E + S ⇌ E(S)    :: EqualAI
            E(S) <--> E(P)  :: EqualAI
        end
    end)
    x_onlyi_at_y = EnzymeRates.AllostericMechanism(@allosteric_mechanism begin
        substrates: S; products: P
        allosteric_regulators: X::OnlyI, Y::NonequalAI
        catalytic_multiplicity: 2
        catalytic_steps: begin
            E + P ⇌ E(P)    :: EqualAI
            E + S ⇌ E(S)    :: EqualAI
            E(S) <--> E(P)  :: EqualAI
        end
        regulatory_site(multiplicity = 2): begin
            ligands: X, Y
        end
    end)
    x_nonequalai_new = EnzymeRates.AllostericMechanism(@allosteric_mechanism begin
        substrates: S; products: P
        allosteric_regulators: X::NonequalAI, Y::NonequalAI
        catalytic_multiplicity: 2
        catalytic_steps: begin
            E + P ⇌ E(P)    :: EqualAI
            E + S ⇌ E(S)    :: EqualAI
            E(S) <--> E(P)  :: EqualAI
        end
    end)
    x_nonequalai_at_y = EnzymeRates.AllostericMechanism(@allosteric_mechanism begin
        substrates: S; products: P
        allosteric_regulators: X::NonequalAI, Y::NonequalAI
        catalytic_multiplicity: 2
        catalytic_steps: begin
            E + P ⇌ E(P)    :: EqualAI
            E + S ⇌ E(S)    :: EqualAI
            E(S) <--> E(P)  :: EqualAI
        end
        regulatory_site(multiplicity = 2): begin
            ligands: X, Y
        end
    end)
    x_equalai_at_y = EnzymeRates.AllostericMechanism(@allosteric_mechanism begin
        substrates: S; products: P
        allosteric_regulators: X::EqualAI, Y::NonequalAI
        catalytic_multiplicity: 2
        catalytic_steps: begin
            E + P ⇌ E(P)    :: EqualAI
            E + S ⇌ E(S)    :: EqualAI
            E(S) <--> E(P)  :: EqualAI
        end
        regulatory_site(multiplicity = 2): begin
            ligands: X, Y
        end
    end)
    from_y = EnzymeRates._expand_add_allosteric_regulator(am_y, rxn)
    @test length(from_y) == 7
    @test Set(from_y) == Set([x_y_apart, shared, x_onlyi_new, x_onlyi_at_y,
                              x_nonequalai_new, x_nonequalai_at_y, x_equalai_at_y])

    # The shared site comes out of each route once, as one mechanism under == and hash.
    @test count(==(shared), from_x) == 1
    @test count(==(shared), from_y) == 1
    via_x = EnzymeRates._make_am_with_added_reg(am_x, :Y, :NonequalAI, 1)
    via_y = EnzymeRates._make_am_with_added_reg(am_y, :X, :OnlyA, 1)
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
    _testhelper_assert_mechanism_invariants(am)

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
    seed_truth = _testhelper_fitted(am)
    @test _testhelper_fitted(only(r_removal)) == seed_truth + 1

    # 4. compilability + invariants on each variant.
    for r in result
        _testhelper_assert_mechanism_invariants(r)
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
    _testhelper_assert_mechanism_invariants(am)

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
    @test _testhelper_param_deltas(am, result) == [0, 0, 1, 1]

    # 3. invariants; the deltas above compile every result.
    foreach(_testhelper_assert_mechanism_invariants, result)

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
    _testhelper_assert_mechanism_invariants(am)

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
    @test _testhelper_param_deltas(am, result) == [0, 0, 1, 1, 1]

    # 3. invariants; the deltas above compile every result.
    foreach(_testhelper_assert_mechanism_invariants, result)

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
    _testhelper_assert_mechanism_invariants(am)
    @test am.catalytic_multiplicity == 4
    @test [EnzymeRates.multiplicity(s) for s in am.regulatory_sites] == [2]

    # _expand_change_allo_state must preserve both
    # catalytic_multiplicity and per-site multiplicity independently.
    result = EnzymeRates._expand_change_allo_state(am)
    @test !isempty(result)
    for r in result
        _testhelper_assert_mechanism_invariants(r)
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

end

# ─── _site_active_states ─────────────────────────────────────────────────
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

# ─── _expand_merge_regulatory_sites ─────────────────────────────────────
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
    kept = filter(c -> EnzymeRates._respects_reg_type(c, merge_rxn), children)
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
    kept = filter(c -> EnzymeRates._respects_reg_type(c, rxn_aa), kids)
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

@testset "no-op: single-site AllostericMechanism" begin
    single = EnzymeRates.AllostericMechanism(
        merge_rxn, copy(EnzymeRates.steps(base)), cat, 4, [site_a])
    @test isempty(EnzymeRates._expand_merge_regulatory_sites(single))
end

@testset "OnlyA+OnlyI merge is derivation-redundant (premise guard)" begin
# The reason the all-keep OnlyA/OnlyI merge is skipped: merged onto one
# site it evaluates to the same rate as on separate sites. Confirm
# numerically over random parameters and concentrations.
mkm(sites) = EnzymeRates.AllostericMechanism(
    merge_rxn, copy(EnzymeRates.steps(base)), cat, 4, sites)
A = EnzymeRates.AllostericRegulator(:A)
I = EnzymeRates.AllostericRegulator(:I)
separate = mkm([EnzymeRates.RegulatorySite([A], 4, [:OnlyA]),
                EnzymeRates.RegulatorySite([I], 4, [:OnlyI])])
merged = mkm([EnzymeRates.RegulatorySite([A, I], 4, [:OnlyA, :OnlyI])])
fps = EnzymeRates.fitted_params(separate)
fpm = EnzymeRates.fitted_params(merged)
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

# ─── _respects_reg_type ──────────────────────────────────────────────────
@testset "_respects_reg_type" begin

@testset "_respects_reg_type — one- and two-ligand sites" begin
    rxn = @enzyme_reaction begin
        substrates: S[C]
        products: P[C]
        allosteric_regulators: A::Activator, A2::Activator, I::Inhibitor, I2::Inhibitor, U
        oligomeric_state: 2
    end
    core = EnzymeRates.AllostericMechanism(@allosteric_mechanism begin
        substrates: S
        products: P
        catalytic_multiplicity: 2
        catalytic_steps: begin
            E + S ⇌ E(S)   :: EqualAI
            E + P ⇌ E(P)   :: EqualAI
            E(S) <--> E(P) :: EqualAI
        end
    end)
    # Whether `core` with one site holding `ligs` at `states` respects the types `r`
    # declares.
    respects(ligs, states, r = rxn) = EnzymeRates._respects_reg_type(
        EnzymeRates._with(core; sites = [EnzymeRates.RegulatorySite(
            EnzymeRates.AllostericRegulator.(collect(ligs)), 2, collect(states))]), r)

    # An :unspecified ligand passes in every state, alone or beside another
    # :unspecified ligand: U is declared untyped and Nope is not declared at all.
    for state in (:OnlyA, :OnlyI, :EqualAI, :NonequalAI)
        @test respects((:U,), (state,))
        @test respects((:Nope,), (state,))
        @test respects((:U, :Nope), (state, :EqualAI))
    end
    # A CompetitiveInhibitor carrying a reg_type (only reachable by direct
    # construction; the DSL restricts type tags to allosteric_regulators:)
    # is not an AllostericRegulator, so its reg_type is ignored.
    rxn_ci = EnzymeRates.EnzymeReaction(
        [EnzymeRates.ReactantAtoms(EnzymeRates.Substrate(:S), [:C => 1]),
         EnzymeRates.ReactantAtoms(EnzymeRates.Product(:P), [:C => 1])],
        [EnzymeRates.RegulatorMults(
            EnzymeRates.CompetitiveInhibitor(:X), [1], :activator)],
        [1])
    @test respects((:X,), (:OnlyI,), rxn_ci)
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
        @test respects(ligs, states) == ok
    end
end

@testset "_respects_reg_type as a filter" begin
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

    kept = filter(c -> EnzymeRates._respects_reg_type(c, rxn),
        Union{EnzymeRates.Mechanism, EnzymeRates.AllostericMechanism}[
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
    m_aa = first(EnzymeRates.init_mechanisms(rxn_aa))
    cat_aa = Symbol[:EqualAI for _ in 1:length(EnzymeRates.steps(m_aa))]
    site_aa = EnzymeRates.RegulatorySite(
        [EnzymeRates.AllostericRegulator(:A1),
         EnzymeRates.AllostericRegulator(:A2)], 2, [:EqualAI, :EqualAI])
    am_aa = EnzymeRates.AllostericMechanism(
        rxn_aa, copy(EnzymeRates.steps(m_aa)), cat_aa, 2, [site_aa])
    @test !EnzymeRates._respects_reg_type(am_aa, rxn_aa)

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
    @test EnzymeRates._respects_reg_type(am_ai, rxn_ai)
end
end

# ═══════════════════════════════════════════════════════════════════════
# 5. Composition (dedup, expand_mechanisms)
# ═══════════════════════════════════════════════════════════════════════

# ─── canonical by construction ─────────────────────────────────────────

@testset "AllostericMechanism — canonical by construction" begin
# Site permutation with DISTINCT multiplicities. The constructor must permute
# cat_allo_states alongside cat_steps (catalytic side) and produce the same
# regulatory_sites ordering (regulatory side) regardless of input site order.
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
@testset "AllostericMechanism — site multiplicity is part of the dedup key" begin
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
_testhelper_assert_mechanism_invariants(am_m2)
_testhelper_assert_mechanism_invariants(am_m4)
@test am_m2 != am_m4
end
end

# ─── mechanism dedup (unique!) ──────────────────────────────────────────
@testset "Dedup" begin

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
    _testhelper_assert_mechanism_invariants(m_seed)
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
    _testhelper_assert_mechanism_invariants(m)
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
    _testhelper_assert_mechanism_invariants(m)
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
    _testhelper_assert_mechanism_invariants(m)
    base_fitted = _testhelper_fitted(m)
    result = EnzymeRates.expand_mechanisms([m], uni_uni_rxn)
    for child in result
        @test _testhelper_fitted(child) >= base_fitted
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
    # ligand named :R (allosteric regulators live in
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

    # Raw children reproduce expand_mechanisms' moves on a Mechanism without the filter.
    raw_act = _testhelper_plain_moves(m_act, rxn_act)
    raw_plain = _testhelper_plain_moves(m_plain, rxn_plain)

    children_act = EnzymeRates.expand_mechanisms([m_act], rxn_act)
    children_plain = EnzymeRates.expand_mechanisms([m_plain], rxn_plain)

    # A designated activator's raw expansion carries an :OnlyI R child (a
    # V-type variant); the reg_type filter drops it, leaving strictly fewer.
    @test any(m -> has_reg_state(m, :R, :OnlyI), raw_act)
    @test !any(m -> has_reg_state(m, :R, :OnlyI), children_act)
    @test children_act == filter(c -> EnzymeRates._respects_reg_type(c, rxn_act), raw_act)
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
    # contributes nothing: expand_mechanisms equals the reg-type-filtered output of the
    # moves that apply to a Mechanism.
    m = first(EnzymeRates.init_mechanisms(uni_uni_rxn))
    raw = _testhelper_plain_moves(m, uni_uni_rxn)
    baseline = filter(c -> EnzymeRates._respects_reg_type(c, uni_uni_rxn), raw)
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
                            first(EnzymeRates.steps(am)[g])))
        for g in EnzymeRates.kinetic_groups(am)
        if EnzymeRates.cat_allo_states(am)[g] === :OnlyA &&
           !EnzymeRates.is_iso(first(EnzymeRates.steps(am)[g])))

    seen = Set{UInt64}()
    frontier = Union{EnzymeRates.Mechanism, EnzymeRates.AllostericMechanism}[]
    for m in EnzymeRates.init_mechanisms(rxn)
        h = hash(m); h in seen || (push!(seen, h); push!(frontier, m))
    end
    maxdistinct = 0
    gen = 0
    # The first generation reaches the target; three also run the regulator moves.
    while !isempty(frontier) && gen < 3
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

uni = enumerate_all_mechanism(uni_uni_rxn; max_params=8)

@testset "Mechanism — Uni-uni full enumeration" begin
    # enumerate_all_mechanism buckets by ACTUAL fitted-parameter count. At
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
        collect(EnzymeRates.init_mechanisms(bi_bi_rxn)))
    @test !isempty(init)
    steps_of(m) = Iterators.flatten(EnzymeRates.steps(m))
    holds_tc(m) = any(s -> !isempty(EnzymeRates.consumed(s)) &&
                           !isempty(EnzymeRates.released(s)), steps_of(m))
    n_ss(m) = count(g -> !EnzymeRates.is_equilibrium(first(g)), EnzymeRates.steps(m))
    kind(m) = _testhelper_holds_iso(m) ? :seed : holds_tc(m) ? :theorell_chance :
              n_ss(m) == 2 ? :merged_two_groups : :merged_three_groups
    # The groups whose every step binds into a form that has no other step.
    function dead_end_groups(m)
        degree = Dict{EnzymeRates.Species, Int}()
        for s in steps_of(m), sp in (EnzymeRates.from_species(s), EnzymeRates.to_species(s))
            degree[sp] = get(degree, sp, 0) + 1
        end
        count(g -> all(s -> degree[EnzymeRates.to_species(s)] == 1, g),
              EnzymeRates.steps(m))
    end
    counts = Set{Int}()
    tally = Dict{Tuple{Symbol, Int}, Int}()
    for m in init
        n = _testhelper_fitted(m)
        push!(counts, n)
        tally[(kind(m), n)] = get(tally, (kind(m), n), 0) + 1
        kind(m) == :theorell_chance && @test n == 5 + dead_end_groups(m)
    end
    # {5,6,7}: the seeds and their two-group merged variants fit 5, the
    # three-group merged variants 6, and a Theorell–Chance variant 5 plus one
    # for each group the eliminated steps leave holding only dead-end steps.
    @test counts == Set([5, 6, 7])
    @test tally == Dict((:seed, 5) => 55, (:merged_two_groups, 5) => 108,
                        (:merged_three_groups, 6) => 56, (:theorell_chance, 5) => 8,
                        (:theorell_chance, 6) => 8, (:theorell_chance, 7) => 4)
    @test all(m -> kind(m) == :theorell_chance || dead_end_groups(m) == 0, init)
end

@testset "Mechanism — With allosteric regulators" begin
    # A uni-uni allosteric reaction enumerates allosteric mechanisms.
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
    total_with_reg = sum(length(v) for v in values(results))
    total_plain = sum(length(v) for v in values(uni))
    @test total_with_reg > total_plain
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

@testset "zero inactive-state numerator keeps the :NonequalAI binding constants" begin
    # K-type allosteric uni-uni: catalytic step is :OnlyA (so the I-state
    # numerator is zero), but binding steps are :NonequalAI.
    # When the I-state numerator is zero, the binding partition function
    # for :NonequalAI groups must still emit their I-state constants in the
    # I-state denominator so they appear in the rate-equation body and in
    # parameters(m).
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
    @test isempty(first(
        EnzymeRates._state_rate_polys(EnzymeRates.AllostericMechanism(m), :I)))
    params = parameters(m)
    # K_I_E_cS_to_E_c_S and K_I_E_cP_to_E_c_P are referenced in the I-state
    # denominator of the body (the binding partition function for :NonequalAI
    # groups is built whether or not the I-state numerator is zero, since the
    # I-state denominator always appears in the rate's denominator).
    @test :K_I_E_cS_to_E_c_S in params
    @test :K_I_E_cP_to_E_c_P in params
end
end


# ─── Expansion-move helpers ──────────────────────────────────────────────
# Helpers and building blocks of the RE→SS flip move and the kinetic-group
# split move. Every fixture is written out with the mechanism macro so the
# mechanism under test is visible in the testset.

@testset "_independent_param_count matches fitted_params on the test specs" begin
    for spec in MECHANISM_TEST_SPECS
        m = spec.mechanism isa EnzymeRates.AllostericEnzymeMechanism ?
            EnzymeRates.AllostericMechanism(spec.mechanism) :
            EnzymeRates.Mechanism(spec.mechanism)
        @test EnzymeRates._independent_param_count(m) ==
              _testhelper_fitted(m)
    end
end

@testset "_unbalanced_blocks" begin
    # Two parallel edges between the same two vertices, of different weight, close a
    # cycle of nonzero weight: one block, unbalanced.
    @test EnzymeRates._unbalanced_blocks(2, [(1, 2), (1, 2)], [1, 2]) == ([1, 1], Set([1]))
    # A triangle whose third edge is stored against the walk 1 → 2 → 3: along its stored
    # direction 3 → 1 it weighs -3, so the cycle weighs 1 + 2 - 3 = 0 and is balanced.
    @test EnzymeRates._unbalanced_blocks(3, [(1, 2), (2, 3), (3, 1)], [1, 2, -3]) ==
          ([1, 1, 1], Set{Int}())
    # The search reaches vertex 2 from vertex 1 over the edge stored as (2, 1), against
    # its direction, so the potential falls by its weight: (0, -1, 0). The cycle
    # 1 → 2 → 3 → 1 weighs -1 + 1 + 0 = 0, balanced.
    @test EnzymeRates._unbalanced_blocks(3, [(2, 1), (2, 3), (3, 1)], [1, 1, 0]) ==
          ([1, 1, 1], Set{Int}())
    # Two components. Vertices 1 and 2 carry an edge and its reverse of opposite weight
    # (balanced). Vertices 3, 4 and 5 form a triangle whose paths 3 → 4 → 5 and 3 → 5
    # weigh 2 and 1 (unbalanced), with a bridge 5 → 6 hanging off it. Blocks are
    # numbered in the order the search closes them: the pair, the bridge, the triangle.
    @test EnzymeRates._unbalanced_blocks(
        6, [(1, 2), (2, 1), (3, 4), (4, 5), (3, 5), (5, 6)], [1, -1, 1, 1, 1, 4]) ==
          ([1, 1, 3, 3, 3, 2], Set([3]))
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
    @test all(EnzymeRates._flux_carrying_groups(
        EnzymeRates._all_steady_state(EnzymeRates.steps(m)), EnzymeRates.reaction(m)))

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
    fc = EnzymeRates._flux_carrying_groups(
        EnzymeRates._all_steady_state(EnzymeRates.steps(leaf)), EnzymeRates.reaction(leaf))
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
            Estar + P ⇌ Estar(P)
            E(A) <--> Estar(P)
            Estar(B) ⇌ E(Q)
        end
    end)
    @test all(EnzymeRates._flux_carrying_groups(
        EnzymeRates._all_steady_state(EnzymeRates.steps(pingpong)),
        EnzymeRates.reaction(pingpong)))

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
    @test all(EnzymeRates._flux_carrying_groups(
        EnzymeRates._all_steady_state(EnzymeRates.steps(mirror)),
        EnzymeRates.reaction(mirror)))

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
    pfc = EnzymeRates._flux_carrying_groups(
        EnzymeRates._all_steady_state(EnzymeRates.steps(pendant)),
        EnzymeRates.reaction(pendant))
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
    @test all(EnzymeRates._flux_carrying_groups(
        EnzymeRates._all_steady_state(EnzymeRates.steps(pp)), EnzymeRates.reaction(pp)))
    # A Theorell–Chance step consumes B and releases P in one edge, so the
    # mechanism has no isomerization step. The cycle E → E(A) → E(Q) → E has
    # weight 1 + 2 + 1 = 4 = one turnover times the four reactants, so it is
    # unbalanced and every step carries flux.
    tc = EnzymeRates.Mechanism(@enzyme_mechanism begin
        substrates: A, B
        products: P, Q
        steps: begin
            E + A ⇌ E(A)
            E(A) + B <--> E(Q) + P
            E + Q ⇌ E(Q)
        end
    end)
    @test all(EnzymeRates._flux_carrying_groups(
        EnzymeRates._all_steady_state(EnzymeRates.steps(tc)), EnzymeRates.reaction(tc)))
    @test EnzymeRates._flux_carrying_groups(EnzymeRates.steps(tc),
                                            EnzymeRates.reaction(tc)) ==
        [!EnzymeRates.is_equilibrium(first(g)) for g in EnzymeRates.steps(tc)]

    # A merged central complex X with two fused steps into it and no
    # isomerization: E + A ⇌ E(A), E(A) + B → X, E(Q) + P → X, E + Q ⇌ E(Q).
    # Every step lies on the unbalanced cycle E → E(A) → X → E(Q) → E.
    merged = EnzymeRates.Mechanism(@enzyme_mechanism begin
        substrates: A, B
        products: P, Q
        steps: begin
            E + A ⇌ E(A)
            E(A) + B <--> E(P, Q)
            E(Q) + P <--> E(P, Q)
            E + Q ⇌ E(Q)
        end
    end)
    @test all(EnzymeRates._flux_carrying_groups(
        EnzymeRates._all_steady_state(EnzymeRates.steps(merged)),
        EnzymeRates.reaction(merged)))

    # An RE shunt. E + A and E + Q are steady-state while the abortive route
    # E(Q) + A ⇌ E(A, Q) ⇌ E(A) + Q stays at rapid equilibrium. The segment
    # {E} joins the rest by two parallel edges of weight 0 (A uptake +1 against
    # E(A)'s offset +1; Q uptake −1 against E(Q)'s offset −1): a balanced block,
    # so neither steady-state binding carries flux. The turnover runs inside the
    # big segment, where the chemistry step is an unbalanced self-loop of
    # weight 4 (E(A, B) sits at A + B, E(P, Q) at −P − Q). On the all-steady-state
    # graph every group carries flux.
    shunt = EnzymeRates.Mechanism(@enzyme_mechanism begin
        substrates: A, B
        products: P, Q
        steps: begin
            E + A <--> E(A)
            E(Q) + A ⇌ E(A, Q)
            E + Q <--> E(Q)
            E(A) + Q ⇌ E(A, Q)
            E(A) + B ⇌ E(A, B)
            E(A, B) <--> E(P, Q)
            E(Q) + P ⇌ E(P, Q)
        end
    end)
    binder(grp) = EnzymeRates.bound_metabolite(first(grp))
    source(grp) = EnzymeRates.name(EnzymeRates.from_species(first(grp)))
    ss_binding(grp) = binder(grp) !== nothing && !EnzymeRates.is_equilibrium(first(grp))
    flags = EnzymeRates._flux_carrying_groups(EnzymeRates.steps(shunt),
                                              EnzymeRates.reaction(shunt))
    for (g, grp) in enumerate(EnzymeRates.steps(shunt))
        if ss_binding(grp)
            @test source(grp) == :E && !flags[g]
        elseif EnzymeRates.is_iso(first(grp))
            @test flags[g]
        else
            @test EnzymeRates.is_equilibrium(first(grp)) && !flags[g]
        end
    end
    @test count(ss_binding, EnzymeRates.steps(shunt)) == 2
    @test all(EnzymeRates._flux_carrying_groups(
        EnzymeRates._all_steady_state(EnzymeRates.steps(shunt)),
        EnzymeRates.reaction(shunt)))

    # A balanced self-loop: E + A is steady-state, but E and E(A) already share
    # a segment through E ⇌ E(Q) ⇌ E(A, Q) ⇌ E(A), so the step is a self-loop
    # of weight 1 + 0 − 1 = 0 and carries no flux. The chemistry self-loop in
    # the same segment has weight 4 and does.
    loop = EnzymeRates.Mechanism(@enzyme_mechanism begin
        substrates: A, B
        products: P, Q
        steps: begin
            E + A <--> E(A)
            E(Q) + A ⇌ E(A, Q)
            (E + Q ⇌ E(Q), E(A) + Q ⇌ E(A, Q))
            E(A) + B ⇌ E(A, B)
            E(A, B) <--> E(P, Q)
            E(Q) + P ⇌ E(P, Q)
        end
    end)
    lflags = EnzymeRates._flux_carrying_groups(EnzymeRates.steps(loop),
                                               EnzymeRates.reaction(loop))
    for (g, grp) in enumerate(EnzymeRates.steps(loop))
        if ss_binding(grp)
            @test !lflags[g]
        elseif EnzymeRates.is_iso(first(grp))
            @test lflags[g]
        end
    end
    @test count(ss_binding, EnzymeRates.steps(loop)) == 1
    @test EnzymeRates._re_segment_count(EnzymeRates.steps(loop)) == 1
end

@testset "MERGE, ELIM and the seed tests on the ordered bi-bi seed" begin
    ER = EnzymeRates
    rxn = @enzyme_reaction begin
        substrates: A[C], B[N]
        products: P[C], Q[N]
    end
    @test ER._reactant_signs(rxn) == Dict(:A => 1, :B => 1, :P => -1, :Q => -1)
    groups(em) = ER.steps(_testhelper_on_reaction(rxn, em))
    seed = groups(@enzyme_mechanism begin
        substrates: A, B; products: P, Q
        steps: begin
            E + A ⇌ E(A)
            E(A) + B ⇌ E(A, B)
            E(A, B) <--> E(P, Q)
            E(Q) + P ⇌ E(P, Q)
            E + Q ⇌ E(Q)
        end
    end)
    iso = only(s for g in seed for s in g if ER.is_iso(s))
    merged = ER._merge_isomerization(seed, iso)
    # The isomerization's product side E(P, Q) takes over E(A, B)'s steps.
    @test Set(s for g in merged for s in g) == Set(s for g in groups(@enzyme_mechanism begin
        substrates: A, B; products: P, Q
        steps: begin
            E + A ⇌ E(A)
            E(A) + B ⇌ E(P, Q)
            E(Q) + P ⇌ E(P, Q)
            E + Q ⇌ E(Q)
        end
    end) for s in g)
    # The merged cycle E → E(A) → E(P, Q) → E(Q) → E takes up A and B and gives off P and
    # Q, weight 1 + 1 + 1 + 1 = 4, all at rapid equilibrium. The seed's cycle runs
    # through its steady-state isomerization, and its RE steps form a tree.
    @test ER._re_turnover_cycle(merged, rxn)        # every step RE: infinite rate
    @test !ER._re_turnover_cycle(seed, rxn)
    variant(ss) = [[ER._with_equilibrium(s, !(ER.bound_metabolite(s) !== nothing &&
                    ER.name(ER.bound_metabolite(s)) in ss)) for s in g] for g in merged]
    clean(gs) = !ER._re_turnover_cycle(gs, rxn) && ER._has_vmax(gs, rxn, ER.Substrate) &&
        ER._has_vmax(gs, rxn, ER.Product) && !ER._chemistry_equilibrates_both_sides(gs, rxn)
    # Track 4 table 4.3: {A, Pˣ}, {Bˣ, Q}, {Bˣ, Pˣ} clean; {A, Q} fails C (B is not
    # needed); {A, Bˣ} has no forward Vmax; {Pˣ, Q} no reverse Vmax; singletons none.
    # The merged groups are A (E + A ⇌ E(A)), Bˣ (E(A) + B → E(P, Q)), Pˣ
    # (E(Q) + P → E(P, Q)) and Q (E + Q ⇌ E(Q)); ˣ marks a step at the merged complex
    # E(P, Q), the one chemistry node. V for substrates turns A and Bˣ RE and fails
    # exactly when Pˣ and Q are RE too, closing the cycle: the steady-state set lies in
    # {A, Bˣ}. V for products fails exactly when it lies in {Pˣ, Q}. Every singleton
    # lies in one of the two. C fails exactly when Bˣ (giving B off from E(P, Q)) and Pˣ
    # (giving P off) both stay RE: among the pairs, {A, Q} only.
    @test clean(variant([:A, :P])) && clean(variant([:B, :Q])) && clean(variant([:B, :P]))
    @test ER._chemistry_equilibrates_both_sides(variant([:A, :Q]), rxn)
    @test !ER._has_vmax(variant([:A, :B]), rxn, ER.Substrate)
    @test !ER._has_vmax(variant([:P, :Q]), rxn, ER.Product)
    @test all(x -> !clean(variant([x])), (:A, :B, :P, :Q))
    # ELIM of the merged complex: one step E(A) + B → E(Q) + P; all three SS is clean.
    x = ER.to_species(iso)
    tc = ER._eliminate_form(merged, x)
    group_of(pred) = only(g for (g, grp) in enumerate(tc) if any(pred, grp))
    gA, gQ = group_of(_testhelper_binds(:A)), group_of(_testhelper_binds(:Q))
    gTC = group_of(_testhelper_is_tc)
    @test length(tc) == 3 &&
        ER.consumed(only(tc[gTC])) == ER.Metabolite[ER.Substrate(:B)] &&
        ER.released(only(tc[gTC])) == ER.Metabolite[ER.Product(:P)]
    @test Set(s for g in tc for s in g) == Set(s for g in groups(@enzyme_mechanism begin
        substrates: A, B; products: P, Q
        steps: begin
            E + A ⇌ E(A)
            E(A) + B ⇌ E(Q) + P
            E + Q ⇌ E(Q)
        end
    end) for s in g)
    @test ER._re_turnover_cycle(tc, rxn)            # E → E(A) → E(Q) → E, all RE
    ss(which) = [g in which ? ER._with_equilibrium.(grp, false) : grp
                 for (g, grp) in enumerate(tc)]
    # The Theorell–Chance step takes up B and gives off P, so V turns it RE for either
    # side. {A, Q, TC}: V for substrates leaves Q steady state, V for products leaves A,
    # and no RE step is left for C. {A, Q}: the RE Theorell–Chance step gives B off one
    # way and P the other, so C fails. {A, TC}: V for substrates closes the cycle with
    # the RE Q step, no forward Vmax; {Q, TC} is the mirror.
    @test clean(ss([gA, gQ, gTC]))
    @test !clean(ss([gA, gQ])) && !clean(ss([gA, gTC])) && !clean(ss([gQ, gTC]))
    @test ER._chemistry_equilibrates_both_sides(ss([gA, gQ]), rxn)
    @test !ER._has_vmax(ss([gA, gTC]), rxn, ER.Substrate)
    @test !ER._has_vmax(ss([gQ, gTC]), rxn, ER.Product)
    # ELIM needs two bindings into the form: E(A)'s second step binds B out of it.
    ea = ER.to_species(only(s for g in merged for s in g if _testhelper_binds(:A)(s)))
    @test ER._eliminate_form(merged, ea) === nothing
    # The seed itself is not degenerate; the probe agrees on every verdict above.
    @test !ER._degenerate(ER.Mechanism(rxn, seed))
    @test !_testhelper_degenerate(ER.Mechanism(rxn, seed))
    pairs = [(:A, :B), (:A, :P), (:A, :Q), (:B, :P), (:B, :Q), (:P, :Q)]
    candidates = vcat([string(x) => variant([x]) for x in (:A, :B, :P, :Q)],
                      [string(x, y) => variant([x, y]) for (x, y) in pairs],
                      ["TC: A Q TC" => ss([gA, gQ, gTC]), "TC: A Q" => ss([gA, gQ]),
                       "TC: A TC" => ss([gA, gTC]), "TC: Q TC" => ss([gQ, gTC])])
    for (label, gs) in candidates
        m = ER.Mechanism(rxn, gs)
        @test (label, ER._degenerate(m)) == (label, !clean(gs))
        @test (label, _testhelper_degenerate(m)) == (label, !clean(gs))
    end
end

@testset "ELIM refuses the uni-uni merged complex" begin
    ER = EnzymeRates
    groups(em) = ER.steps(_testhelper_on_reaction(uni_uni_rxn, em))
    seed = groups(@enzyme_mechanism begin
        substrates: S; products: P
        steps: begin
            E + S ⇌ E(S)
            E(S) <--> E(P)
            E + P ⇌ E(P)
        end
    end)
    iso = only(s for g in seed for s in g if ER.is_iso(s))
    merged = ER._merge_isomerization(seed, iso)
    # S binds E into the merged complex E(P), and P binds E into it as well.
    @test Set(s for g in merged for s in g) == Set(s for g in groups(@enzyme_mechanism begin
        substrates: S; products: P
        steps: begin
            E + S ⇌ E(P)
            E + P ⇌ E(P)
        end
    end) for s in g)
    # Both of E(P)'s steps start at E, so a Theorell–Chance step E + S → E + P would join
    # E to itself.
    @test ER._eliminate_form(merged, ER.to_species(iso)) === nothing
end

@testset "_re_turnover_cycle: a balanced rapid-equilibrium cycle runs no turnover" begin
    ER = EnzymeRates
    # Random-order binding of A and B at rapid equilibrium closes the cycle
    # E → E(A) → E(A, B) ← E(B) ← E of weight 1 + 1 − 1 − 1 = 0: it runs no net reaction.
    # With the isomerization RE as well, E → E(A) → E(A, B) → E(P, Q) → E(Q) → E has
    # weight 4, one turnover, in the same block.
    m = ER.Mechanism(@enzyme_mechanism begin
        substrates: A, B; products: P, Q
        steps: begin
            E + A ⇌ E(A)
            E + B ⇌ E(B)
            E(A) + B ⇌ E(A, B)
            E(B) + A ⇌ E(A, B)
            E(A, B) <--> E(P, Q)
            E(Q) + P ⇌ E(P, Q)
            E + Q ⇌ E(Q)
        end
    end)
    @test !ER._re_turnover_cycle(ER.steps(m), ER.reaction(m))
    @test ER._re_turnover_cycle([ER._with_equilibrium.(g, true) for g in ER.steps(m)],
                                ER.reaction(m))
end

@testset "C on the ping-pong seed and its single flips" begin
    ER = EnzymeRates
    groups(em) = ER.steps(_testhelper_on_reaction(bi_bi_pp_rxn, em))
    seed = groups(@enzyme_mechanism begin
        substrates: A, B; products: P, Q
        steps: begin
            E + A ⇌ E(A)
            E(A) <--> E(P; residual = A - P)
            E(; residual = A - P) + P ⇌ E(P; residual = A - P)
            E(; residual = A - P) + B ⇌ E(B; residual = A - P)
            E(B; residual = A - P) ⇌ E(Q)
            E + Q ⇌ E(Q)
        end
    end)
    flip(x) = [[_testhelper_binds(x)(s) ? ER._with_equilibrium(s, false) : s for s in g]
               for g in seed]
    # The steady-state isomerization E(A) → E(P; res) cuts the one catalytic cycle, and
    # no step it would turn RE touches a substrate or a product, so both maximal rates
    # exist. The RE isomerization E(B; res) ⇌ E(Q) makes {E(B; res), E(Q)} a chemistry
    # node, left by the RE B binding (traversed away from E(B; res) it gives B off) and
    # by the RE Q binding (away from E(Q) it gives Q off): C fails. At zero products the
    # forms past the steady-state step drain to E and E(A), so B is not needed.
    rxn = bi_bi_pp_rxn
    @test ER._has_vmax(seed, rxn, ER.Substrate) && ER._has_vmax(seed, rxn, ER.Product)
    @test ER._chemistry_equilibrates_both_sides(seed, rxn)
    @test ER._degenerate(ER.Mechanism(rxn, seed))
    @test _testhelper_degenerate(ER.Mechanism(rxn, seed))
    # Flipping B or Q leaves the node one RE exit, Q's or B's: C holds (the B flip is the
    # steady-state half of the SRR ping-pong). Flipping A or P leaves both exits RE: C
    # still fails. V never turns the steady-state isomerization RE (it touches no
    # metabolite), so it still cuts the cycle and every flip keeps both maximal rates.
    for x in (:A, :B, :P, :Q)
        gs = flip(x)
        m = ER.Mechanism(rxn, gs)
        fails_c = x in (:A, :P)
        @test ER._has_vmax(gs, rxn, ER.Substrate) && ER._has_vmax(gs, rxn, ER.Product)
        @test (x, ER._chemistry_equilibrates_both_sides(gs, rxn)) == (x, fails_c)
        @test (x, ER._degenerate(m)) == (x, fails_c)
        @test (x, _testhelper_degenerate(m)) == (x, fails_c)
    end
    # A conformational variant is read on its active-state graph, the seed's: still
    # degenerate.
    am = ER.AllostericMechanism(@allosteric_mechanism begin
        substrates: A, B ; products: P, Q ; catalytic_multiplicity: 1
        catalytic_steps: begin
            E + A ⇌ E(A)                                          :: EqualAI
            E(A) <--> E(P; residual = A - P)                      :: OnlyA
            E(; residual = A - P) + P ⇌ E(P; residual = A - P)    :: EqualAI
            E(; residual = A - P) + B ⇌ E(B; residual = A - P)    :: EqualAI
            E(B; residual = A - P) ⇌ E(Q)                         :: OnlyA
            E + Q ⇌ E(Q)                                          :: EqualAI
        end
    end)
    @test ER._degenerate(am)
    # The probe reads the derived allosteric law. Free E and E(; residual = A - P) lie in
    # one rapid-equilibrium segment of the active state; the law divides each
    # conformation by its free enzyme's weight, so at zero products it does not depend on
    # B (the gate in test/allosteric_ground_truth.jl checks it against mass action).
    @test _testhelper_degenerate(am) == ER._degenerate(am)
end

@testset "_seed_variants: the uni-uni seed has none" begin
    ER = EnzymeRates
    seed = _testhelper_on_reaction(uni_uni_rxn, @enzyme_mechanism begin
        substrates: S; products: P
        steps: begin
            E + S ⇌ E(S)
            E(S) <--> E(P)
            E + P ⇌ E(P)
        end
    end)
    # The merged base holds S (E + S → E(P), fused) and P (E + P ⇌ E(P)), both units.
    # Either group alone at steady state has no maximal rate one way: V turns it rapid
    # equilibrium and closes the turnover cycle E → E(P) → E with the other. {S, P} passes
    # every test but is the lumping twin of the seed (both steps at E(P) steady state, each
    # alone in its group), so it is skipped. Both steps at E(P) start at E, so ELIM gives no
    # Theorell–Chance base.
    @test isempty(ER._seed_variants(seed))
    only_s = _testhelper_on_reaction(uni_uni_rxn, @enzyme_mechanism begin
        substrates: S; products: P
        steps: begin
            E + S <--> E(P)
            E + P ⇌ E(P)
        end
    end)
    only_p = _testhelper_on_reaction(uni_uni_rxn, @enzyme_mechanism begin
        substrates: S; products: P
        steps: begin
            E + S ⇌ E(P)
            E + P <--> E(P)
        end
    end)
    twin = _testhelper_on_reaction(uni_uni_rxn, @enzyme_mechanism begin
        substrates: S; products: P
        steps: begin
            E + S <--> E(P)
            E + P <--> E(P)
        end
    end)
    @test !ER._has_vmax(ER.steps(only_s), uni_uni_rxn, ER.Substrate)
    @test !ER._has_vmax(ER.steps(only_p), uni_uni_rxn, ER.Product)
    @test _testhelper_degenerate(only_s) && _testhelper_degenerate(only_p)
    # The twin is skipped as the seed's twin, not as degenerate.
    @test !ER._degenerate(twin) && !_testhelper_degenerate(twin)
end

@testset "_seed_variants: the undecorated ordered seed" begin
    ER = EnzymeRates
    seed = _testhelper_on_reaction(bi_bi_rxn, @enzyme_mechanism begin
        substrates: A, B; products: P, Q
        steps: begin
            E + A ⇌ E(A)
            E(A) + B ⇌ E(A, B)
            E(A, B) <--> E(P, Q)
            E(Q) + P ⇌ E(P, Q)
            E + Q ⇌ E(Q)
        end
    end)
    # The merged base holds A (E + A ⇌ E(A)), Bˣ (E(A) + B → E(P, Q)), Pˣ
    # (E(Q) + P → E(P, Q)) and Q (E + Q ⇌ E(Q)), every group a unit. Every singleton has no
    # maximal rate one way; among the pairs {A, Q} fails C, {A, Bˣ} has no forward and
    # {Pˣ, Q} no reverse maximal rate (verdicts and probes pinned in "MERGE, ELIM and the
    # seed tests on the ordered bi-bi seed"). {A, Pˣ} and {Bˣ, Q} are variants; {Bˣ, Pˣ} is
    # the seed's lumping twin and is skipped; every triple holds a passing pair. The
    # Theorell–Chance base (A, Q, TC = E(A) + B → E(Q) + P) fails every set of one or two
    # groups and gives {A, Q, TC}. Each variant fits two SS groups and two RE groups
    # (2 + 2 + 1 + 1) or three SS groups (2 + 2 + 2), less the Haldane relation: 5.
    a_p = _testhelper_on_reaction(bi_bi_rxn, @enzyme_mechanism begin
        substrates: A, B; products: P, Q
        steps: begin
            E + A <--> E(A)
            E(A) + B ⇌ E(P, Q)
            E(Q) + P <--> E(P, Q)
            E + Q ⇌ E(Q)
        end
    end)
    b_q = _testhelper_on_reaction(bi_bi_rxn, @enzyme_mechanism begin
        substrates: A, B; products: P, Q
        steps: begin
            E + A ⇌ E(A)
            E(A) + B <--> E(P, Q)
            E(Q) + P ⇌ E(P, Q)
            E + Q <--> E(Q)
        end
    end)
    tc = _testhelper_on_reaction(bi_bi_rxn, @enzyme_mechanism begin
        substrates: A, B; products: P, Q
        steps: begin
            E + A <--> E(A)
            E(A) + B <--> E(Q) + P
            E + Q <--> E(Q)
        end
    end)
    variants = ER._seed_variants(seed)
    @test length(variants) == 3
    @test Set(variants) == Set([a_p, b_q, tc])
    for v in variants
        fitted = _testhelper_fitted(v)
        @test (fitted, _testhelper_identifiable_rank(v)) == (5, 5)
        @test !_testhelper_degenerate(v)
    end
end

@testset "_seed_variants: an ordered seed whose B group is shared" begin
    ER = EnzymeRates
    seed = _testhelper_on_reaction(bi_bi_rxn, @enzyme_mechanism begin
        substrates: A, B; products: P, Q
        steps: begin
            E + A ⇌ E(A)
            (E(A) + B ⇌ E(A, B), E(Q) + B ⇌ E(B, Q))
            E(A, B) <--> E(P, Q)
            E(Q) + P ⇌ E(P, Q)
            E + Q ⇌ E(Q)
        end
    end)
    # The merged base holds A, Bˣ (E(A) + B → E(P, Q) with the abortive E(Q) + B ⇌ E(B, Q)),
    # Pˣ and Q; the abortive step is a dead end, and its group carries flux through the
    # fused step, so every group is a unit. The verdicts are the undecorated seed's: the
    # abortive step lies on no cycle, so it changes no turnover test, and it leaves no node
    # (E(P, Q) is the one node, with exits Bˣ and Pˣ). Singletons fail V; {A, Q} fails C;
    # {A, Bˣ} and {Pˣ, Q} fail V. {A, Pˣ}, {Bˣ, Q} and {Bˣ, Pˣ} are variants: {Bˣ, Pˣ} is no
    # lumping twin, because the Bˣ group holds two steps. Each fits 2 + 2 + 1 + 1 − 1 = 5.
    # ELIM of E(P, Q) leaves the abortive step alone in its group, which carries no flux
    # in the all-steady-state base, so it is no unit and stays RE. As in the undecorated
    # base, every set of one or two of A, Q and TC fails, and {A, Q, TC} fits
    # 2 + 2 + 2 + 1 − 1 = 6.
    b_p = _testhelper_on_reaction(bi_bi_rxn, @enzyme_mechanism begin
        substrates: A, B; products: P, Q
        steps: begin
            E + A ⇌ E(A)
            (E(A) + B <--> E(P, Q), E(Q) + B <--> E(B, Q))
            E(Q) + P <--> E(P, Q)
            E + Q ⇌ E(Q)
        end
    end)
    a_p = _testhelper_on_reaction(bi_bi_rxn, @enzyme_mechanism begin
        substrates: A, B; products: P, Q
        steps: begin
            E + A <--> E(A)
            (E(A) + B ⇌ E(P, Q), E(Q) + B ⇌ E(B, Q))
            E(Q) + P <--> E(P, Q)
            E + Q ⇌ E(Q)
        end
    end)
    b_q = _testhelper_on_reaction(bi_bi_rxn, @enzyme_mechanism begin
        substrates: A, B; products: P, Q
        steps: begin
            E + A ⇌ E(A)
            (E(A) + B <--> E(P, Q), E(Q) + B <--> E(B, Q))
            E(Q) + P ⇌ E(P, Q)
            E + Q <--> E(Q)
        end
    end)
    tc = _testhelper_on_reaction(bi_bi_rxn, @enzyme_mechanism begin
        substrates: A, B; products: P, Q
        steps: begin
            E + A <--> E(A)
            E(Q) + B ⇌ E(B, Q)
            E(A) + B <--> E(Q) + P
            E + Q <--> E(Q)
        end
    end)
    variants = ER._seed_variants(seed)
    @test length(variants) == 4
    @test Set(variants) == Set([b_p, a_p, b_q, tc])
    for (v, n) in ((b_p, 5), (a_p, 5), (b_q, 5), (tc, 6))
        fitted = _testhelper_fitted(v)
        @test (fitted, _testhelper_identifiable_rank(v)) == (n, n)
        @test !_testhelper_degenerate(v)
    end
    # The rejected candidates, each degenerate by the probe as well.
    iso = only(s for g in ER.steps(seed) for s in g if ER.is_iso(s))
    base = [ER._with_equilibrium.(g, true)
            for g in ER._merge_isomerization(ER.steps(seed), iso)]
    tc_base = ER._eliminate_form(base, ER.to_species(iso))
    bA, bQ = _testhelper_binds(:A), _testhelper_binds(:Q)
    rejected = vcat(
        [_testhelper_flip_matching(base, [_testhelper_binds(x)]) for x in (:A, :B, :P, :Q)],
        [_testhelper_flip_matching(base, _testhelper_binds.(xs))
         for xs in ([:A, :Q], [:A, :B], [:P, :Q])],
        [_testhelper_flip_matching(tc_base, ps)
         for ps in ([bA], [bQ], [_testhelper_is_tc], [bA, bQ],
                    [bA, _testhelper_is_tc], [bQ, _testhelper_is_tc])])
    for gs in rejected
        m = ER.Mechanism(bi_bi_rxn, gs)
        @test ER._degenerate(m) && _testhelper_degenerate(m)
    end
end

@testset "_seed_variants: the undecorated ping-pong seed" begin
    ER = EnzymeRates
    seed = _testhelper_on_reaction(bi_bi_pp_rxn, @enzyme_mechanism begin
        substrates: A, B; products: P, Q
        steps: begin
            E + A ⇌ E(A)
            E(A) <--> E(P; residual = A - P)
            E(; residual = A - P) + P ⇌ E(P; residual = A - P)
            E(; residual = A - P) + B ⇌ E(B; residual = A - P)
            E(B; residual = A - P) ⇌ E(Q)
            E + Q ⇌ E(Q)
        end
    end)
    # Merging both isomerizations gives A (E + A → E(P; res), fused), P, B
    # (E(; res) + B → E(Q), fused) and Q, all units, with the merged complexes E(P; res) and
    # E(Q) as C's nodes. Every singleton fails V. {A, Q} (PP-AQ) and {B, P} (PP-BP) pass:
    # each node keeps one RE exit (P from E(P; res) and B from E(Q), or A and Q), and each
    # side's steps turned RE leave the other side's steady-state group on the cycle.
    # {A, P} fails C at E(Q) (RE exits B and Q), {B, Q} at E(P; res) (A and P); {A, B} has
    # no forward and {P, Q} no reverse maximal rate. Every triple holds PP-AQ or PP-BP. No
    # merged complex has both its steps steady state in a variant, so neither is a twin.
    # ELIM of E(P; res) gives TC1 = E + A → E(; res) + P beside B and Q; of E(Q),
    # TC2 = E(; res) + B → E + Q beside A and P. In each Theorell–Chance base a set of one
    # or two groups leaves the TC step RE (C fails) or turns both of one side's groups RE
    # (V fails); all three steady state pass. Each variant fits 2 + 2 + 1 + 1 − 1 or
    # 2 + 2 + 2 − 1 = 5.
    pp_aq = _testhelper_on_reaction(bi_bi_pp_rxn, @enzyme_mechanism begin
        substrates: A, B; products: P, Q
        steps: begin
            E + A <--> E(P; residual = A - P)
            E(; residual = A - P) + P ⇌ E(P; residual = A - P)
            E(; residual = A - P) + B ⇌ E(Q)
            E + Q <--> E(Q)
        end
    end)
    pp_bp = _testhelper_on_reaction(bi_bi_pp_rxn, @enzyme_mechanism begin
        substrates: A, B; products: P, Q
        steps: begin
            E + A ⇌ E(P; residual = A - P)
            E(; residual = A - P) + P <--> E(P; residual = A - P)
            E(; residual = A - P) + B <--> E(Q)
            E + Q ⇌ E(Q)
        end
    end)
    half_tc1 = _testhelper_on_reaction(bi_bi_pp_rxn, @enzyme_mechanism begin
        substrates: A, B; products: P, Q
        steps: begin
            E + A <--> E(; residual = A - P) + P
            E(; residual = A - P) + B <--> E(Q)
            E + Q <--> E(Q)
        end
    end)
    half_tc2 = _testhelper_on_reaction(bi_bi_pp_rxn, @enzyme_mechanism begin
        substrates: A, B; products: P, Q
        steps: begin
            E + A <--> E(P; residual = A - P)
            E(; residual = A - P) + P <--> E(P; residual = A - P)
            E(; residual = A - P) + B <--> E + Q
        end
    end)
    variants = ER._seed_variants(seed)
    @test length(variants) == 4
    @test Set(variants) == Set([pp_aq, pp_bp, half_tc1, half_tc2])
    for v in variants
        fitted = _testhelper_fitted(v)
        @test (fitted, _testhelper_identifiable_rank(v)) == (5, 5)
        @test !_testhelper_degenerate(v)
    end
    # The rejected candidates, each degenerate by the probe as well.
    isos = [s for g in ER.steps(seed) for s in g if ER.is_iso(s)]
    base = [ER._with_equilibrium.(g, true)
            for g in foldl(ER._merge_isomerization, isos; init = ER.steps(seed))]
    tc_bases = [ER._eliminate_form(base, ER.to_species(s)) for s in isos]
    # A Theorell–Chance base holds its TC step and the bindings of the other half's two
    # metabolites, x and y.
    function tc_rejected(gs)
        x, y = [ER.name(ER.bound_metabolite(first(g)))
                for g in gs if ER.is_binding(first(g))]
        bx, by = _testhelper_binds(x), _testhelper_binds(y)
        [_testhelper_flip_matching(gs, ps)
         for ps in ([_testhelper_is_tc], [bx], [by], [_testhelper_is_tc, bx],
                    [_testhelper_is_tc, by], [bx, by])]
    end
    rejected = vcat(
        [_testhelper_flip_matching(base, [_testhelper_binds(x)]) for x in (:A, :B, :P, :Q)],
        [_testhelper_flip_matching(base, _testhelper_binds.(xs))
         for xs in ([:A, :P], [:B, :Q], [:A, :B], [:P, :Q])],
        tc_rejected(tc_bases[1]), tc_rejected(tc_bases[2]))
    for gs in rejected
        m = ER.Mechanism(bi_bi_pp_rxn, gs)
        @test ER._degenerate(m) && _testhelper_degenerate(m)
    end
end

@testset "_seed_variants: the ordered/random seed" begin
    ER = EnzymeRates
    seed = _testhelper_on_reaction(bi_bi_rxn, @enzyme_mechanism begin
        substrates: A, B; products: P, Q
        steps: begin
            E + A ⇌ E(A)
            E(A) + B ⇌ E(A, B)
            E(A, B) <--> E(P, Q)
            (E(Q) + P ⇌ E(P, Q), E + P ⇌ E(P))
            (E(P) + Q ⇌ E(P, Q), E + Q ⇌ E(Q))
        end
    end)
    # The merged base holds A, Bˣ (E(A) + B → E(P, Q)), Pˣ (E(Q) + P → E(P, Q),
    # E + P ⇌ E(P)) and Qˣ (E(P) + Q → E(P, Q), E + Q ⇌ E(Q)), all units. The node E(P, Q)
    # has three exits: Bˣ gives B off, Pˣ's step gives P off and Qˣ's gives Q off. Every
    # singleton fails V.
    # {A, Pˣ} keeps the RE exits Bˣ and Qˣ, {A, Qˣ} keeps Bˣ and Pˣ: both fail C. {A, Bˣ}
    # turned RE with Pˣ and Qˣ closes a turnover cycle (no forward maximal rate), and
    # {Pˣ, Qˣ} likewise backward. {Bˣ, Pˣ} passes: V for products turns Pˣ and Qˣ RE, which
    # closes only the balanced square E → E(P) → E(P, Q) ← E(Q) ← E, and the node keeps one
    # RE exit, Q; {Bˣ, Qˣ} is its mirror. E(P, Q) has three steps, so neither is a twin.
    # Among the triples only {A, Pˣ, Qˣ} has three failing pairs, and it passes: the node
    # keeps the one RE exit Bˣ. Fitted: 1 + 2 + 2 + 1 − 1 = 5 for the pairs, and
    # 2 + 1 + 2 + 2 − 1 = 6 for the triple; the shared P and Q groups satisfy the square's
    # Wegscheider relation by themselves.
    b_p = _testhelper_on_reaction(bi_bi_rxn, @enzyme_mechanism begin
        substrates: A, B; products: P, Q
        steps: begin
            E + A ⇌ E(A)
            E(A) + B <--> E(P, Q)
            (E(Q) + P <--> E(P, Q), E + P <--> E(P))
            (E(P) + Q ⇌ E(P, Q), E + Q ⇌ E(Q))
        end
    end)
    b_q = _testhelper_on_reaction(bi_bi_rxn, @enzyme_mechanism begin
        substrates: A, B; products: P, Q
        steps: begin
            E + A ⇌ E(A)
            E(A) + B <--> E(P, Q)
            (E(Q) + P ⇌ E(P, Q), E + P ⇌ E(P))
            (E(P) + Q <--> E(P, Q), E + Q <--> E(Q))
        end
    end)
    a_p_q = _testhelper_on_reaction(bi_bi_rxn, @enzyme_mechanism begin
        substrates: A, B; products: P, Q
        steps: begin
            E + A <--> E(A)
            E(A) + B ⇌ E(P, Q)
            (E(Q) + P <--> E(P, Q), E + P <--> E(P))
            (E(P) + Q <--> E(P, Q), E + Q <--> E(Q))
        end
    end)
    variants = ER._seed_variants(seed)
    @test length(variants) == 3
    @test Set(variants) == Set([b_p, b_q, a_p_q])
    for (v, n) in ((b_p, 5), (b_q, 5), (a_p_q, 6))
        fitted = _testhelper_fitted(v)
        @test (fitted, _testhelper_identifiable_rank(v)) == (n, n)
        @test !_testhelper_degenerate(v)
    end
    # ELIM refuses the merged complex: it has three steps.
    iso = only(s for g in ER.steps(seed) for s in g if ER.is_iso(s))
    base = [ER._with_equilibrium.(g, true)
            for g in ER._merge_isomerization(ER.steps(seed), iso)]
    @test count(s -> ER.to_species(s) == ER.to_species(iso), Iterators.flatten(base)) == 3
    @test ER._eliminate_form(base, ER.to_species(iso)) === nothing
    # The rejected candidates, each degenerate by the probe as well.
    flip(gs, xs) = [any(s -> any(x -> _testhelper_binds(x)(s), xs), g) ?
                    ER._with_equilibrium.(g, false) : g for g in gs]
    for xs in ([:A], [:B], [:P], [:Q], [:A, :P], [:A, :Q], [:A, :B], [:P, :Q])
        m = ER.Mechanism(bi_bi_rxn, flip(base, xs))
        @test (xs, ER._degenerate(m), _testhelper_degenerate(m)) == (xs, true, true)
    end
end

@testset "_seed_candidate_screen agrees with the group predicates on every candidate" begin
    ER = EnzymeRates
    # The screen reads the turnover, maximal-rate and chemistry tests of a candidate from
    # index arrays. On every subset of the groups of eight bases flipped to steady state it
    # must give the verdict of `_re_turnover_cycle`, `_has_vmax` both ways and
    # `_chemistry_equilibrates_both_sides` on the flipped groups.
    function bases(seed)
        isos = [s for g in ER.steps(seed) for s in g if ER.is_iso(s)]
        merged = [ER._with_equilibrium.(g, true)
                  for g in foldl(ER._merge_isomerization, isos; init = ER.steps(seed))]
        tcs = [ER._eliminate_form(merged, ER.to_species(s)) for s in isos]
        vcat([merged], filter(!isnothing, tcs))
    end
    ordered = _testhelper_on_reaction(bi_bi_rxn, @enzyme_mechanism begin
        substrates: A, B; products: P, Q
        steps: begin
            E + A ⇌ E(A)
            E(A) + B ⇌ E(A, B)
            E(A, B) <--> E(P, Q)
            E(Q) + P ⇌ E(P, Q)
            E + Q ⇌ E(Q)
        end
    end)
    shared_b = _testhelper_on_reaction(bi_bi_rxn, @enzyme_mechanism begin
        substrates: A, B; products: P, Q
        steps: begin
            E + A ⇌ E(A)
            (E(A) + B ⇌ E(A, B), E(Q) + B ⇌ E(B, Q))
            E(A, B) <--> E(P, Q)
            E(Q) + P ⇌ E(P, Q)
            E + Q ⇌ E(Q)
        end
    end)
    ordered_random = _testhelper_on_reaction(bi_bi_rxn, @enzyme_mechanism begin
        substrates: A, B; products: P, Q
        steps: begin
            E + A ⇌ E(A)
            E(A) + B ⇌ E(A, B)
            E(A, B) <--> E(P, Q)
            (E(Q) + P ⇌ E(P, Q), E + P ⇌ E(P))
            (E(P) + Q ⇌ E(P, Q), E + Q ⇌ E(Q))
        end
    end)
    ping_pong = _testhelper_on_reaction(bi_bi_pp_rxn, @enzyme_mechanism begin
        substrates: A, B; products: P, Q
        steps: begin
            E + A ⇌ E(A)
            E(A) <--> E(P; residual = A - P)
            E(; residual = A - P) + P ⇌ E(P; residual = A - P)
            E(; residual = A - P) + B ⇌ E(B; residual = A - P)
            E(B; residual = A - P) ⇌ E(Q)
            E + Q ⇌ E(Q)
        end
    end)
    # Bases: merged and Theorell–Chance for the two ordered seeds, merged only for the
    # ordered/random seed, merged and both halves for the ping-pong seed.
    all_bases = [(ER.reaction(seed), b)
                 for seed in (ordered, shared_b, ordered_random, ping_pong)
                 for b in bases(seed)]
    @test length.(last.(all_bases)) == [4, 3, 4, 4, 4, 4, 3, 3]
    # The screen's verdict and the predicates' verdict on every subset of `base`'s groups.
    function sweep(rxn, base)
        screen = ER._seed_candidate_screen(base, rxn)
        map(0:(2^length(base) - 1)) do bits
            ss = BitVector([isodd(bits >> (g - 1)) for g in eachindex(base)])
            gs = [ss[g] ? ER._with_equilibrium.(grp, false) : grp
                  for (g, grp) in enumerate(base)]
            expected = !ER._re_turnover_cycle(gs, rxn) &&
                ER._has_vmax(gs, rxn, ER.Substrate) && ER._has_vmax(gs, rxn, ER.Product) &&
                !ER._chemistry_equilibrates_both_sides(gs, rxn)
            (ss, screen(gs), expected)
        end
    end
    verdicts = Bool[]
    for (rxn, base) in all_bases, (ss, screened, expected) in sweep(rxn, base)
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
    pp_seeds = filter(_testhelper_holds_iso, ER.init_mechanisms(bi_bi_pp_rxn))
    pp_bases = [b for seed in pp_seeds for b in bases(seed)]
    fused(b, side) = [s for s in Iterators.flatten(b) if ER.is_binding(s) &&
                      ER._is_chemistry(s) && ER.bound_metabolite(s) isa side]
    function three_step_complex(b)
        n = Dict{ER.Species, Int}()
        for s in Iterators.flatten(b), sp in (ER.from_species(s), ER.to_species(s))
            n[sp] = get(n, sp, 0) + 1
        end
        any(s -> n[ER.to_species(s)] >= 3, fused(b, ER.Substrate))
    end
    @test (length(pp_seeds), length(pp_bases)) == (62, 86)
    @test count(b -> !isempty(fused(b, ER.Product)), pp_bases) == 8
    @test count(b -> !isempty(fused(b, ER.Product)) && three_step_complex(b), pp_bases) == 8
    @test count(three_step_complex, pp_bases) == 35 + 8
    outcomes = [o for b in pp_bases for o in sweep(bi_bi_pp_rxn, b)]
    @test length(outcomes) == 1344
    @test isempty([(ss, screened) for (ss, screened, expected) in outcomes
                   if screened != expected])
    @test count(o -> o[3], outcomes) == 484
end

@testset "_productive_twin" begin
    ER = EnzymeRates
    # The form named `nm` among the ends of `m`'s steps.
    form(m, nm) = only(unique(sp for g in ER.steps(m) for s in g
                              for sp in (ER.from_species(s), ER.to_species(s))
                              if ER.name(sp) == nm))
    Ainh, Binh, Sinh, Qinh, I = ER.CompetitiveInhibitor.((:A, :B, :S, :Q, :I))

    # Ordered bi-bi. A copy of A at E has the composition of E(A), its twin; at E(Q)
    # it is new; a copy of B at E is new; a foreign inhibitor is never a twin.
    ordered = ER.Mechanism(@enzyme_mechanism begin
        substrates: A, B
        products: P, Q
        steps: begin
            E + A ⇌ E(A)
            E(A) + B ⇌ E(A, B)
            E(A, B) <--> E(P, Q)
            E(Q) + P ⇌ E(P, Q)
            E + Q ⇌ E(Q)
        end
    end)
    twin = ER._productive_twin(ER.steps(ordered))
    @test twin(form(ordered, :E), Ainh) == form(ordered, :EA)
    @test twin(form(ordered, :EQ), Ainh) === nothing
    @test twin(form(ordered, :E), Binh) === nothing
    @test twin(form(ordered, :E), I) === nothing && twin(form(ordered, :EQ), I) === nothing

    # With the abortive complex E(A, Q) present, the copy at E(Q) is a twin too.
    abortive = ER.Mechanism(@enzyme_mechanism begin
        substrates: A, B
        products: P, Q
        steps: begin
            (E + A ⇌ E(A), E(Q) + A ⇌ E(A, Q))
            E(A) + B ⇌ E(A, B)
            E(A, B) <--> E(P, Q)
            E(Q) + P ⇌ E(P, Q)
            (E + Q ⇌ E(Q), E(A) + Q ⇌ E(A, Q))
        end
    end)
    @test ER._productive_twin(ER.steps(abortive))(form(abortive, :EQ), Ainh) ==
        form(abortive, :EAQ)

    # A copy whose complex already exists is judged against the productive forms only.
    with_copy = ER.Mechanism(@enzyme_mechanism begin
        substrates: A, B
        products: P, Q
        steps: begin
            E + A ⇌ E(A)
            E(A) + B ⇌ E(A, B)
            E(A, B) <--> E(P, Q)
            E(Q) + P ⇌ E(P, Q)
            (E + Q ⇌ E(Q), E(A::Inh) + Q ⇌ E(A::Inh, Q))
            (E + A::Inh ⇌ E(A::Inh), E(Q) + A::Inh ⇌ E(A::Inh, Q))
        end
    end)
    twin_c = ER._productive_twin(ER.steps(with_copy))
    @test twin_c(form(with_copy, :E), Ainh) == form(with_copy, :EA)
    @test twin_c(form(with_copy, :EQ), Ainh) === nothing

    # The offsets key. E ⇌ E* is a rapid-equilibrium isomerization, so E and E*
    # have equal offsets; a copy of S at E has a composition no form has, yet its
    # weight is proportional to E*(S)'s. With the isomerization at steady state the
    # two forms sit in different segments and the copy is not a twin.
    iso = ER.Mechanism(@enzyme_mechanism begin
        substrates: S
        products: P
        steps: begin
            E ⇌ Estar
            Estar + S ⇌ Estar(S)
            Estar(S) <--> E(P)
            E + P ⇌ E(P)
        end
    end)
    estar_s = only(ER.to_species(s) for g in ER.steps(iso) for s in g
                   if ER.bound_metabolite(s) !== nothing &&
                      ER.name(ER.bound_metabolite(s)) == :S)
    @test ER._productive_twin(ER.steps(iso))(form(iso, :E), Sinh) == estar_s
    # Twins are found by offsets; composition only breaks a tie between forms with the
    # complex's offsets. This twin differs from the complex in conformation or residual,
    # so no form has the complex's composition.
    @test (ER.conformation(estar_s), ER.residual(estar_s)) !=
          (ER.conformation(form(iso, :E)), ER.residual(form(iso, :E)))
    iso_ss = ER.Mechanism(@enzyme_mechanism begin
        substrates: S
        products: P
        steps: begin
            E <--> Estar
            Estar + S ⇌ Estar(S)
            Estar(S) <--> E(P)
            E + P ⇌ E(P)
        end
    end)
    @test ER._productive_twin(ER.steps(iso_ss))(form(iso_ss, :E), Sinh) === nothing

    # A form of the complex's composition outside the site's segment is no twin. With
    # A bound at steady state, E(A) has the composition of E·A* but sits in another
    # segment, its weight set by the steady state rather than by K·A·[E]: the two
    # weights are not proportional.
    ss_bound = ER.Mechanism(@enzyme_mechanism begin
        substrates: A, B
        products: P, Q
        steps: begin
            E + A <--> E(A)
            E(A) + B ⇌ E(A, B)
            E(A, B) <--> E(P, Q)
            E(Q) + P ⇌ E(P, Q)
            E + Q ⇌ E(Q)
        end
    end)
    @test ER._productive_twin(ER.steps(ss_bound))(form(ss_bound, :E), Ainh) === nothing

    # Ping-pong with the second chemistry step at rapid equilibrium, an abortive
    # E(B, P; residual) and a copy of P at E. A copy of Q at E(P::Inh) has the new
    # composition {P, Q}, but E(B; residual) ⇌ E(Q) puts E(B, P; residual) at the
    # offsets P + Q of E's segment, the same as E(P::Inh)·Q: a twin. At
    # E(; residual) the copy is not.
    pp = ER.Mechanism(@enzyme_mechanism begin
        substrates: A, B
        products: P, Q
        steps: begin
            E + A ⇌ E(A)
            E(A) <--> E(P; residual = A - P)
            E(P; residual = A - P) ⇌ E(; residual = A - P) + P
            E(; residual = A - P) + B ⇌ E(B; residual = A - P)
            E(B; residual = A - P) ⇌ E(Q)
            E(Q) ⇌ E + Q
            E(B; residual = A - P) + P ⇌ E(B, P; residual = A - P)
            E + P::Inh ⇌ E(P::Inh)
        end
    end)
    copy_form = only(ER.to_species(s) for g in ER.steps(pp) for s in g
                     if ER.bound_metabolite(s) isa ER.CompetitiveInhibitor)
    free_res = only(ER.from_species(s) for g in ER.steps(pp) for s in g
                    if ER.bound_metabolite(s) !== nothing &&
                       ER.name(ER.bound_metabolite(s)) == :B)
    abortive_res = only(ER.to_species(s) for g in ER.steps(pp) for s in g
                        if ER.bound_metabolite(s) !== nothing &&
                           ER.name(ER.bound_metabolite(s)) == :P &&
                           ER.has_residual(ER.to_species(s)) &&
                           length(ER.bound(ER.to_species(s))) == 2)
    twin_pp = ER._productive_twin(ER.steps(pp))
    @test twin_pp(copy_form, Qinh) == abortive_res
    @test twin_pp(free_res, Qinh) === nothing

    # Groups that bind a copy redundantly: only at twin sites, with a consistent gauge.
    @test isempty(ER._redundant_copy_groups(with_copy))

    # Two copies. E(A, Q*) has the composition of E(A*, Q), but a form bound to a
    # competitive inhibitor is no twin source: only productive forms are. So the Q*
    # site at E(A) is not a twin, the Q* site at E duplicates E(Q), and both copy
    # groups hold a site with no twin.
    two_copies = ER.Mechanism(@enzyme_mechanism begin
        substrates: A, B
        products: P, Q
        steps: begin
            (E + A ⇌ E(A), E(Q::Inh) + A ⇌ E(A, Q::Inh))
            E(A) + B ⇌ E(A, B)
            E(A, B) <--> E(P, Q)
            E(Q) + P ⇌ E(P, Q)
            E + Q ⇌ E(Q)
            (E + Q::Inh ⇌ E(Q::Inh), E(A) + Q::Inh ⇌ E(A, Q::Inh))
            (E(A) + A::Inh ⇌ E(A, A::Inh), E(Q) + A::Inh ⇌ E(A::Inh, Q))
        end
    end)
    twin_two = ER._productive_twin(ER.steps(two_copies))
    @test twin_two(form(two_copies, :EA), Qinh) === nothing
    @test twin_two(form(two_copies, :E), Qinh) == form(two_copies, :EQ)
    @test twin_two(form(two_copies, :EQ), Ainh) === nothing
    @test isempty(ER._redundant_copy_groups(two_copies))

    # Conformational states. A copy binds the active state always and the inactive state
    # unless it is `:OnlyA`. It is redundant when every complex has a twin in every state
    # where it binds and the gauge holds over the states together. In each fixture below,
    # the copy E + S* ⇌ E(S*) and every other kinetic group are single steps, and a single
    # step has one rescaling, so the gauge holds in each state where E(S*) has a twin. The
    # states agree too. The only `:EqualAI` groups are E + P ⇌ E(P), which touches neither
    # E(S) nor E(S*) and takes factor 1 at both ends in each state, and an `:EqualAI` copy,
    # which takes 1 at E and its one factor s at E(S*) in each. So a copy here is
    # redundant exactly when E(S*) has a twin in every state where it binds. In the active
    # state that twin is E(S), formed from E by the RE S binding.
    # S binding `:OnlyA`: the inactive state holds E, E(P) and an `:EqualAI` copy's
    # E(S*), and no productive form bound to S. The copy is new there and kept.
    onlya_with_equalai_copy = ER.AllostericMechanism(@allosteric_mechanism begin
        substrates: S
        products: P
        catalytic_inhibitors: S
        catalytic_multiplicity: 2
        catalytic_steps: begin
            E + S ⇌ E(S)                :: OnlyA
            E(S) <--> E(P)              :: OnlyA
            E + P ⇌ E(P)                :: EqualAI
            E + S::Inh ⇌ E(S::Inh)      :: EqualAI
        end
    end)
    @test isempty(ER._redundant_copy_groups(onlya_with_equalai_copy))
    @test _testhelper_fitted(onlya_with_equalai_copy) ==
        _testhelper_identifiable_rank(onlya_with_equalai_copy)
    # An `:OnlyA` copy binds the active state only, where E(S) is its twin: redundant,
    # its constant one fitted parameter above the rank.
    onlya_with_onlya_copy = ER.AllostericMechanism(@allosteric_mechanism begin
        substrates: S
        products: P
        catalytic_inhibitors: S
        catalytic_multiplicity: 2
        catalytic_steps: begin
            E + S ⇌ E(S)                :: OnlyA
            E(S) <--> E(P)              :: OnlyA
            E + P ⇌ E(P)                :: EqualAI
            E + S::Inh ⇌ E(S::Inh)      :: OnlyA
        end
    end)
    @test !isempty(ER._redundant_copy_groups(onlya_with_onlya_copy))
    @test _testhelper_fitted(onlya_with_onlya_copy) ==
        _testhelper_identifiable_rank(onlya_with_onlya_copy) + 1
    # S binding `:NonequalAI`: E(S) exists in both states and is the twin in each, so an
    # `:EqualAI` copy is redundant.
    nonequal_with_equalai_copy = ER.AllostericMechanism(@allosteric_mechanism begin
        substrates: S
        products: P
        catalytic_inhibitors: S
        catalytic_multiplicity: 2
        catalytic_steps: begin
            E + S ⇌ E(S)                :: NonequalAI
            E(S) <--> E(P)              :: OnlyA
            E + P ⇌ E(P)                :: EqualAI
            E + S::Inh ⇌ E(S::Inh)      :: EqualAI
        end
    end)
    @test !isempty(ER._redundant_copy_groups(nonequal_with_equalai_copy))
    @test _testhelper_fitted(nonequal_with_equalai_copy) ==
        _testhelper_identifiable_rank(nonequal_with_equalai_copy) + 1
    # Every catalytic group `:OnlyA`: the inactive state has no catalytic step, yet it
    # holds its free enzyme, so an `:EqualAI` copy binds E there, where E(S*) has no
    # twin: kept. An `:OnlyA` copy binds the active state only and is redundant.
    all_onlya_with_equalai_copy = ER.AllostericMechanism(@allosteric_mechanism begin
        substrates: S
        products: P
        catalytic_inhibitors: S
        catalytic_multiplicity: 2
        catalytic_steps: begin
            E + S ⇌ E(S)                :: OnlyA
            E(S) <--> E(P)              :: OnlyA
            E + P ⇌ E(P)                :: OnlyA
            E + S::Inh ⇌ E(S::Inh)      :: EqualAI
        end
    end)
    @test isempty(ER._redundant_copy_groups(all_onlya_with_equalai_copy))
    @test _testhelper_fitted(all_onlya_with_equalai_copy) ==
        _testhelper_identifiable_rank(all_onlya_with_equalai_copy)
    all_onlya_with_onlya_copy = ER.AllostericMechanism(@allosteric_mechanism begin
        substrates: S
        products: P
        catalytic_inhibitors: S
        catalytic_multiplicity: 2
        catalytic_steps: begin
            E + S ⇌ E(S)                :: OnlyA
            E(S) <--> E(P)              :: OnlyA
            E + P ⇌ E(P)                :: OnlyA
            E + S::Inh ⇌ E(S::Inh)      :: OnlyA
        end
    end)
    @test !isempty(ER._redundant_copy_groups(all_onlya_with_onlya_copy))
    @test _testhelper_fitted(all_onlya_with_onlya_copy) ==
        _testhelper_identifiable_rank(all_onlya_with_onlya_copy) + 1
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
    copy_group(m) = only(g for (g, grp) in enumerate(ER.steps(m))
                         if ER.bound_metabolite(first(grp)) isa ER.CompetitiveInhibitor)

    # B1: ordered bi-bi with A* at E. The twin E(A) is formed
    # by E + A ⇌ E(A), alone in its group, so the A group's ratio is 1/ρ, the B group's
    # ρ and the chemistry's forward factor ρ, each over one step: the gauge exists.
    # 6 fitted, rank 5.
    b1 = ER.Mechanism(@enzyme_mechanism begin
        substrates: A, B
        products: P, Q
        steps: begin
            E + A ⇌ E(A)
            E(A) + B ⇌ E(A, B)
            E(A, B) <--> E(P, Q)
            E(Q) + P ⇌ E(P, Q)
            E + Q ⇌ E(Q)
            E + A::Inh ⇌ E(A::Inh)
        end
    end)
    @test ER._redundant_copy_groups(b1) == [copy_group(b1)]
    @test _testhelper_fitted(b1) == 6 && _testhelper_identifiable_rank(b1) == 5

    # Case 3: the A group also binds at E(Q), forming the abortive E(A, Q). The copy at
    # E duplicates E(A), but the A group then holds a step with ratio 1/ρ (at E) and
    # one with ratio 1 (at E(Q), no copy): no gauge, and the copy's constant is
    # identifiable. 6 fitted, rank 6.
    case3 = ER.Mechanism(@enzyme_mechanism begin
        substrates: A, B
        products: P, Q
        steps: begin
            (E + A ⇌ E(A), E(Q) + A ⇌ E(A, Q))
            E(A) + B ⇌ E(A, B)
            E(A, B) <--> E(P, Q)
            E(Q) + P ⇌ E(P, Q)
            (E + Q ⇌ E(Q), E(A) + Q ⇌ E(A, Q))
            E + A::Inh ⇌ E(A::Inh)
        end
    end)
    g3 = copy_group(case3)
    @test ER._all_twin(ER.steps(case3)[g3], ER._productive_twin(ER.steps(case3)))
    @test isempty(ER._redundant_copy_groups(case3))
    @test _testhelper_fitted(case3) == 6 && _testhelper_identifiable_rank(case3) == 6

    # B14: the same parent with A* at E and E(Q), the Q binding E + Q ⇌ E(Q) mirrored
    # onto the copy forms. Both twins, E(A) and E(A, Q), are formed by the one A
    # group, so they share ρ: the A group's two steps both have ratio 1/ρ, and the Q
    # group's three steps (E → E(Q), E(A) → E(A, Q) and the mirror between the two
    # complexes, both at s) all have ratio 1. The gauge exists. 6 fitted, rank 5.
    b14 = ER.Mechanism(@enzyme_mechanism begin
        substrates: A, B
        products: P, Q
        steps: begin
            (E + A ⇌ E(A), E(Q) + A ⇌ E(A, Q))
            E(A) + B ⇌ E(A, B)
            E(A, B) <--> E(P, Q)
            E(Q) + P ⇌ E(P, Q)
            (E + Q ⇌ E(Q), E(A) + Q ⇌ E(A, Q), E(A::Inh) + Q ⇌ E(A::Inh, Q))
            (E + A::Inh ⇌ E(A::Inh), E(Q) + A::Inh ⇌ E(A::Inh, Q))
        end
    end)
    @test ER._redundant_copy_groups(b14) == [copy_group(b14)]
    @test _testhelper_fitted(b14) == 6 && _testhelper_identifiable_rank(b14) == 5

    # H1: A binds E in one group and E(B) in another, so the twins E(A) and E(A, B)
    # have different ρ. The B group shared by E(A) + B ⇌ E(A, B) (ratio ρ₁/ρ₂) and
    # E(Q) + B ⇌ E(B, Q) (ratio 1) cannot absorb both: no gauge. (The B group at E
    # and its mirror between the complexes have ratio 1.) 7 fitted, rank 7.
    h1 = ER.Mechanism(@enzyme_mechanism begin
        substrates: A, B
        products: P, Q
        steps: begin
            E + A ⇌ E(A)
            E(B) + A ⇌ E(A, B)
            (E + B ⇌ E(B), E(A::Inh) + B ⇌ E(A::Inh, B))
            (E(A) + B ⇌ E(A, B), E(Q) + B ⇌ E(B, Q))
            E(A, B) <--> E(P, Q)
            E(Q) + P ⇌ E(P, Q)
            E + Q ⇌ E(Q)
            (E + A::Inh ⇌ E(A::Inh), E(B) + A::Inh ⇌ E(A::Inh, B))
        end
    end)
    gh = copy_group(h1)
    @test ER._all_twin(ER.steps(h1)[gh], ER._productive_twin(ER.steps(h1)))
    @test isempty(ER._redundant_copy_groups(h1))
    @test _testhelper_fitted(h1) == 7 && _testhelper_identifiable_rank(h1) == 7

    # A step that leaves a complex takes the complex's factor s, never its twin's. B14
    # plus P* at the complex E(A*) and at the twin E(A, Q): the P* group holds
    # E(A*) → E(A*, P*) at s beside E(A, Q) → E(A, P*, Q) at ρ, two ratios, so no gauge,
    # and the copy's constant shows (7 fitted, rank 7, as without the copy). The
    # mechanism is a valid parent.
    p_at_twin = ER.Mechanism(@enzyme_mechanism begin
        substrates: A, B
        products: P, Q
        steps: begin
            (E + A ⇌ E(A), E(Q) + A ⇌ E(A, Q))
            E(A) + B ⇌ E(A, B)
            E(A, B) <--> E(P, Q)
            E(Q) + P ⇌ E(P, Q)
            (E + Q ⇌ E(Q), E(A) + Q ⇌ E(A, Q), E(A::Inh) + Q ⇌ E(A::Inh, Q))
            (E + A::Inh ⇌ E(A::Inh), E(Q) + A::Inh ⇌ E(A::Inh, Q))
            (E(A::Inh) + P::Inh ⇌ E(A::Inh, P::Inh), E(A, Q) + P::Inh ⇌ E(A, P::Inh, Q))
        end
    end)
    @test isempty(ER._redundant_copy_groups(p_at_twin))
    @test ER._assert_emission_rules(p_at_twin) === nothing
    @test _testhelper_fitted(p_at_twin) == 7 &&
          _testhelper_identifiable_rank(p_at_twin) == 7
    # P* at E(A*) and at its own twin E(A) instead: E(A*, P*) pairs with E(A, P*), and a
    # phantom remains (7 fitted, rank 6). The gauge does not reach it: it reads
    # E(A*) → E(A*, P*) at s beside E(A) → E(A, P*) at ρ, and the copy is kept.
    p_at_own_twin = ER.Mechanism(@enzyme_mechanism begin
        substrates: A, B
        products: P, Q
        steps: begin
            (E + A ⇌ E(A), E(Q) + A ⇌ E(A, Q))
            E(A) + B ⇌ E(A, B)
            E(A, B) <--> E(P, Q)
            E(Q) + P ⇌ E(P, Q)
            (E + Q ⇌ E(Q), E(A) + Q ⇌ E(A, Q), E(A::Inh) + Q ⇌ E(A::Inh, Q))
            (E + A::Inh ⇌ E(A::Inh), E(Q) + A::Inh ⇌ E(A::Inh, Q))
            (E(A::Inh) + P::Inh ⇌ E(A::Inh, P::Inh), E(A) + P::Inh ⇌ E(A, P::Inh))
        end
    end)
    @test isempty(ER._redundant_copy_groups(p_at_own_twin))
    @test _testhelper_fitted(p_at_own_twin) == 7 &&
          _testhelper_identifiable_rank(p_at_own_twin) == 6

    # An RE mirror joins two complexes, both at s, so its ratio is 1, as is its parent
    # step's between two sites, whatever the twins' classes. Random-order A in two
    # groups with A* at {E, E(B)}: twins E(A) and E(A, B) of two classes, and
    # E + B ⇌ E(B) mirrored onto the complexes. Redundant: 7 fitted, rank 6.
    two_classes = ER.Mechanism(@enzyme_mechanism begin
        substrates: A, B
        products: P, Q
        steps: begin
            E + A ⇌ E(A)
            E(B) + A ⇌ E(A, B)
            (E + B ⇌ E(B), E(A::Inh) + B ⇌ E(A::Inh, B))
            E(A) + B ⇌ E(A, B)
            E(A, B) <--> E(P, Q)
            E(Q) + P ⇌ E(P, Q)
            E + Q ⇌ E(Q)
            (E + A::Inh ⇌ E(A::Inh), E(B) + A::Inh ⇌ E(A::Inh, B))
        end
    end)
    @test ER._redundant_copy_groups(two_classes) == [copy_group(two_classes)]
    @test _testhelper_fitted(two_classes) == 7 &&
          _testhelper_identifiable_rank(two_classes) == 6
    # B14 with the A bindings and E(A) + Q ⇌ E(A, Q) each in a group of their own: the
    # twins E(A) and E(A, Q) have two classes, and the Q mirror joins their complexes.
    # Redundant: 7 fitted, rank 6. With the Q group at steady state the mirror carries
    # flux between the complexes, its forward constant needs s beside its parent step's
    # 1, and the gauge fails.
    split_b14 = ER.Mechanism(@enzyme_mechanism begin
        substrates: A, B
        products: P, Q
        steps: begin
            E + A ⇌ E(A)
            E(Q) + A ⇌ E(A, Q)
            E(A) + B ⇌ E(A, B)
            E(A, B) <--> E(P, Q)
            E(Q) + P ⇌ E(P, Q)
            (E + Q ⇌ E(Q), E(A::Inh) + Q ⇌ E(A::Inh, Q))
            E(A) + Q ⇌ E(A, Q)
            (E + A::Inh ⇌ E(A::Inh), E(Q) + A::Inh ⇌ E(A::Inh, Q))
        end
    end)
    ss_mirror = ER.Mechanism(@enzyme_mechanism begin
        substrates: A, B
        products: P, Q
        steps: begin
            E + A ⇌ E(A)
            E(Q) + A ⇌ E(A, Q)
            E(A) + B ⇌ E(A, B)
            E(A, B) <--> E(P, Q)
            E(Q) + P ⇌ E(P, Q)
            (E + Q <--> E(Q), E(A::Inh) + Q <--> E(A::Inh, Q))
            E(A) + Q ⇌ E(A, Q)
            (E + A::Inh ⇌ E(A::Inh), E(Q) + A::Inh ⇌ E(A::Inh, Q))
        end
    end)
    @test ER._redundant_copy_groups(split_b14) == [copy_group(split_b14)]
    @test _testhelper_fitted(split_b14) == 7 &&
          _testhelper_identifiable_rank(split_b14) == 6
    @test isempty(ER._redundant_copy_groups(ss_mirror))
    @test _testhelper_fitted(ss_mirror) == 8 &&
          _testhelper_identifiable_rank(ss_mirror) == 6
    # One factor for all of a copy's complexes, whatever their twins' classes: the
    # previous mechanism plus P* at the two complexes E(A*) and E(A*, Q), the Q binding
    # mirrored between the P* complexes. Both P* steps leave a complex at s, one ratio,
    # and the gauge exists: 8 fitted, rank 7.
    p_at_complexes = ER.Mechanism(@enzyme_mechanism begin
        substrates: A, B
        products: P, Q
        steps: begin
            E + A ⇌ E(A)
            E(Q) + A ⇌ E(A, Q)
            E(A) + B ⇌ E(A, B)
            E(A, B) <--> E(P, Q)
            E(Q) + P ⇌ E(P, Q)
            (E + Q ⇌ E(Q), E(A::Inh) + Q ⇌ E(A::Inh, Q),
             E(A::Inh, P::Inh) + Q ⇌ E(A::Inh, P::Inh, Q))
            E(A) + Q ⇌ E(A, Q)
            (E + A::Inh ⇌ E(A::Inh), E(Q) + A::Inh ⇌ E(A::Inh, Q))
            (E(A::Inh) + P::Inh ⇌ E(A::Inh, P::Inh),
             E(A::Inh, Q) + P::Inh ⇌ E(A::Inh, P::Inh, Q))
        end
    end)
    a_group = only(g for (g, grp) in enumerate(ER.steps(p_at_complexes))
                   if ER.bound_metabolite(first(grp)) == ER.CompetitiveInhibitor(:A))
    @test ER._redundant_copy_groups(p_at_complexes) == [a_group]
    @test _testhelper_fitted(p_at_complexes) == 8 &&
          _testhelper_identifiable_rank(p_at_complexes) == 7

    # Case 1, uni-uni with S* at E: the twin E(S) is formed by the S group alone, and
    # the chemistry leaves it alone. 4 fitted, rank 3.
    uni = ER.Mechanism(@enzyme_mechanism begin
        substrates: S
        products: P
        steps: begin
            E + S ⇌ E(S)
            E(S) <--> E(P)
            E + P ⇌ E(P)
            E + S::Inh ⇌ E(S::Inh)
        end
    end)
    @test ER._redundant_copy_groups(uni) == [copy_group(uni)]
    @test _testhelper_fitted(uni) == 4 && _testhelper_identifiable_rank(uni) == 3

    # A fused RE binding forms a twin as a plain one does. Merged uni-uni with A* at E:
    # E + A ⇌ E(P) is fused and rapid equilibrium, E + P → E(P) steady state, so E(P)
    # shares E's segment at the offsets {A: 1}, those of E·A*, and is its twin though
    # no form has the complex's composition. The fused group forms that twin alone, so
    # the twin's class is that group: its ratio is 1/ρ, P's binding takes ρ at E(P),
    # and the copy s at its complex; the gauge exists. The law
    # E_t·k·(Keq·A − P)/(1 + A/K_A + A/K*), k the P binding's rate constant, sees K_A
    # and K* only through 1/K_A + 1/K*: 3 fitted, rank 2, against 2 and 2 without A*.
    fused_re = ER.Mechanism(@enzyme_mechanism begin
        substrates: A
        products: P
        steps: begin
            E + A ⇌ E(P)
            E + P <--> E(P)
            E + A::Inh ⇌ E(A::Inh)
        end
    end)
    fused_parent = ER.Mechanism(@enzyme_mechanism begin
        substrates: A
        products: P
        steps: begin
            E + A ⇌ E(P)
            E + P <--> E(P)
        end
    end)
    copy_g = copy_group(fused_re)
    fused_g = only(k for (k, grp) in enumerate(ER.steps(fused_re))
                   if ER._is_chemistry(first(grp)))
    p_g = only(setdiff(eachindex(ER.steps(fused_re)), (copy_g, fused_g)))
    copy_step = only(ER.steps(fused_re)[copy_g])
    twin = ER._productive_twin(ER.steps(fused_re))
    @test twin(ER.from_species(copy_step), ER.bound_metabolite(copy_step)) ==
        ER.to_species(only(ER.steps(fused_re)[fused_g]))
    @test ER._gauge_rescaling(ER.steps(fused_re), copy_g, twin, identity) ==
        Dict(copy_g => (nothing, copy_g), fused_g => (nothing, fused_g),
             p_g => (nothing, fused_g))
    @test ER._redundant_copy_groups(fused_re) == [copy_g]
    @test _testhelper_fitted(fused_re) == 3 && _testhelper_identifiable_rank(fused_re) == 2
    @test _testhelper_fitted(fused_parent) ==
          _testhelper_identifiable_rank(fused_parent) == 2

    # Ping-pong with the second chemistry step at rapid equilibrium, an abortive
    # E(B, P; residual) and P* at E, then Q* at E(P*). E(P*, Q*) has no form's
    # composition, but it shares E's segment and the offsets P + Q with
    # E(B, P; residual), an offsets twin, its own class. Only E(B; residual) + P ⇌
    # E(B, P; residual) touches the twin, alone in its group: the gauge exists, and
    # Q* adds one phantom to the parent's one (9 fitted, rank 7, against 8 and 7).
    pp = ER.Mechanism(@enzyme_mechanism begin
        substrates: A, B
        products: P, Q
        steps: begin
            E + A ⇌ E(A)
            E(A) <--> E(P; residual = A - P)
            E(P; residual = A - P) ⇌ E(; residual = A - P) + P
            E(; residual = A - P) + B ⇌ E(B; residual = A - P)
            E(B; residual = A - P) ⇌ E(Q)
            E(Q) ⇌ E + Q
            E(B; residual = A - P) + P ⇌ E(B, P; residual = A - P)
            E + P::Inh ⇌ E(P::Inh)
        end
    end)
    pp_q = ER.Mechanism(@enzyme_mechanism begin
        substrates: A, B
        products: P, Q
        steps: begin
            E + A ⇌ E(A)
            E(A) <--> E(P; residual = A - P)
            E(P; residual = A - P) ⇌ E(; residual = A - P) + P
            E(; residual = A - P) + B ⇌ E(B; residual = A - P)
            E(B; residual = A - P) ⇌ E(Q)
            E(Q) ⇌ E + Q
            E(B; residual = A - P) + P ⇌ E(B, P; residual = A - P)
            E + P::Inh ⇌ E(P::Inh)
            E(P::Inh) + Q::Inh ⇌ E(P::Inh, Q::Inh)
        end
    end)
    @test isempty(ER._redundant_copy_groups(pp))
    q_group = only(g for (g, grp) in enumerate(ER.steps(pp_q))
                   if ER.bound_metabolite(first(grp)) == ER.CompetitiveInhibitor(:Q))
    @test ER._redundant_copy_groups(pp_q) == [q_group]
    @test _testhelper_fitted(pp) == 8 && _testhelper_identifiable_rank(pp) == 7
    @test _testhelper_fitted(pp_q) == 9 && _testhelper_identifiable_rank(pp_q) == 7

    # A twin must share its site's segment. A binds E at steady state and E·A*
    # duplicates E(A) by composition only: E(A) lies in another segment, so merging
    # E·A* into it is no reparameterization, and the copy's constant shows (7 fitted,
    # rank 7, one above the mechanism without the copy).
    ss_twin = ER.Mechanism(@enzyme_mechanism begin
        substrates: A, B
        products: P, Q
        steps: begin
            E + A <--> E(A)
            (E + P ⇌ E(P), E(Q) + P ⇌ E(P, Q))
            (E + Q ⇌ E(Q), E(P) + Q ⇌ E(P, Q))
            E(A) + B ⇌ E(A, B)
            E(A, B) <--> E(P, Q)
            E + A::Inh ⇌ E(A::Inh)
        end
    end)
    ss_parent = ER.Mechanism(@enzyme_mechanism begin
        substrates: A, B
        products: P, Q
        steps: begin
            E + A <--> E(A)
            (E + P ⇌ E(P), E(Q) + P ⇌ E(P, Q))
            (E + Q ⇌ E(Q), E(P) + Q ⇌ E(P, Q))
            E(A) + B ⇌ E(A, B)
            E(A, B) <--> E(P, Q)
        end
    end)
    @test isempty(ER._redundant_copy_groups(ss_twin))
    @test _testhelper_fitted(ss_twin) == 7 && _testhelper_identifiable_rank(ss_twin) == 7
    @test _testhelper_fitted(ss_parent) == 6 &&
          _testhelper_identifiable_rank(ss_parent) == 6

    # An `:EqualAI` group has one constant in both states, so it must take one rescaling.
    # A binds `:NonequalAI`, so the twin E(A) has ρ = K_A/(K_A + K*) with a different K_A
    # in each state, while the `:EqualAI` B binding leaves E(A) in both: the active state
    # asks K_B·ρ_A, the inactive K_B·ρ_I. Each state alone has a gauge; the two together
    # have none, and the copy adds an identifiable constant (rank one above the
    # mechanism without it).
    rxn2 = @enzyme_reaction begin
        substrates: A[C], B[N]
        products: P[C], Q[N]
        dead_end_inhibitors: A
        oligomeric_state: 2
    end
    lift2(em) = ER.AllostericMechanism(rxn2, ER.steps(em), ER.cat_allo_states(em), 2,
                                       ER.RegulatorySite[])
    tied = lift2(ER.AllostericMechanism(@allosteric_mechanism begin
        substrates: A, B
        products: P, Q
        catalytic_inhibitors: A
        catalytic_multiplicity: 2
        catalytic_steps: begin
            E + A ⇌ E(A)                :: NonequalAI
            E(A) + B ⇌ E(A, B)          :: EqualAI
            E(A, B) <--> E(P, Q)        :: NonequalAI
            E(Q) + P ⇌ E(P, Q)          :: EqualAI
            E + Q ⇌ E(Q)                :: EqualAI
            E + A::Inh ⇌ E(A::Inh)      :: EqualAI
        end
    end))
    tied_parent = lift2(ER.AllostericMechanism(@allosteric_mechanism begin
        substrates: A, B
        products: P, Q
        catalytic_multiplicity: 2
        catalytic_steps: begin
            E + A ⇌ E(A)                :: NonequalAI
            E(A) + B ⇌ E(A, B)          :: EqualAI
            E(A, B) <--> E(P, Q)        :: NonequalAI
            E(Q) + P ⇌ E(P, Q)          :: EqualAI
            E + Q ⇌ E(Q)                :: EqualAI
        end
    end))
    @test isempty(ER._redundant_copy_groups(tied))
    @test _testhelper_identifiable_rank(tied) ==
        _testhelper_identifiable_rank(tied_parent) + 1
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

    # Ordered SS bi-bi with substrate A binding E(Q) as an abortive complex at
    # rapid equilibrium: E(A, Q) sits one A above E(Q) in its segment, and the
    # term rooted there holds E → E(A), so the equation carries A². A binds its
    # catalytic site, so the power counts even though E(A, Q) is a dead end.
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
    @test !EnzymeRates._hyperbolic_catalysis(dead_end)

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
    # E(Q), as the enumerator groups it. The abortive step went to steady state
    # with its group; the tree toward E(A, Q) holds E → E(A) and E(Q) → E(A, Q),
    # both binding A, so the equation carries A².
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
    @test !EnzymeRates._hyperbolic_catalysis(grouped_dead_end)

    # Ping-pong whose second chemistry step is at rapid equilibrium, with B also
    # bound as an abortive complex on E(Q) in the same kinetic group as its
    # catalytic binding. One segment spans both halves: E sits one B above
    # E(; residual) and E(B, Q) two above, so the segment weight carries B².
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
    @test !EnzymeRates._hyperbolic_catalysis(pingpong_abortive)

    # Ordered SS bi-bi with the abortive complex E(A, Q) reachable from both
    # E(A) + Q and E(Q) + A, each step grouped with its metabolite's catalytic
    # binding as the enumerator writes it. The tree toward E(A, Q) holds
    # E → E(A) and E(Q) → E(A, Q), both binding A, so the equation carries A².
    two_sided_abortive = EnzymeRates.Mechanism(@enzyme_mechanism begin
        substrates: A, B
        products: P, Q
        steps: begin
            (E + A <--> E(A), E(Q) + A <--> E(A, Q))
            E(A) + B <--> E(A, B)
            (E + Q <--> E(Q), E(A) + Q <--> E(A, Q))
            E(Q) + P <--> E(P, Q)
            E(A, B) <--> E(P, Q)
        end
    end)
    @test !EnzymeRates._hyperbolic_catalysis(two_sided_abortive)

    # A competitive-inhibitor square (I on E and on E(A), with A binding the
    # inhibitor-bound form too) leaves the verdict unchanged: an inhibitor binds
    # a site of its own, so every step touching a form that carries it is left
    # out.
    rxn_with_i = @enzyme_reaction begin
        substrates: A[C], B[N]
        products: P[C], Q[N]
        dead_end_inhibitors: I
    end
    ordered_ss = EnzymeRates.Mechanism(@enzyme_mechanism begin
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
    with_i = EnzymeRates._expand_add_dead_end_regulator(ordered_ss, rxn_with_i)
    @test !isempty(with_i)
    @test all(EnzymeRates._hyperbolic_catalysis, with_i)

    # A substrate declared as a dead-end inhibitor binds its inhibitor site under
    # the substrate's own name. Its bindings put A² and A³ into the derived
    # denominator, but those powers come from the inhibitor site, which the
    # predicate leaves out: the four placements on the ordered scheme pass, with
    # degrees 1, 2, 2 and 3. A binds at steady state, so E(A) sits in a segment of
    # its own and E·A* is no twin of it: the copy at E alone is emitted, and it adds
    # a phantom (10 fitted, rank 7, over the parent's 9 and 7).
    rxn_a_inhibits = @enzyme_reaction begin
        substrates: A[C], B[N]
        products: P[C], Q[N]
        dead_end_inhibitors: A
    end
    copy_sites(k) = Set(EnzymeRates.name(EnzymeRates.from_species(s))
                        for grp in EnzymeRates.steps(k) for s in grp
                        if EnzymeRates.bound_metabolite(s) isa
                           EnzymeRates.CompetitiveInhibitor)
    with_a = EnzymeRates._expand_add_dead_end_regulator(ordered_ss, rxn_a_inhibits)
    den_a_degree(k) = maximum(get(Dict(mono), :A, 0) for mono in keys(
        EnzymeRates._raw_symbolic_rate_polys(k, EnzymeRates._step_parameters(k),
            EnzymeRates._build_wegscheider_rename_map(k))[2]))
    @test length(with_a) == 4
    @test Set(copy_sites.(with_a)) ==
        Set([Set([:E]), Set([:E, :EQ]), Set([:EA, :EQ]), Set([:E, :EA])])
    @test all(EnzymeRates._hyperbolic_catalysis, with_a)
    @test sort(den_a_degree.(with_a)) == [1, 2, 2, 3]
    at_e(ks) = only(k for k in ks if copy_sites(k) == Set([:E]))
    @test _testhelper_identifiable_rank(ordered_ss) == 7 &&
          _testhelper_fitted(ordered_ss) == 9
    @test _testhelper_identifiable_rank(at_e(with_a)) == 7 &&
          _testhelper_fitted(at_e(with_a)) == 10

    # Inhibitor-bound forms leave, but catalytic-site powers stay: adding the
    # inhibitor role of A to the abortive-complex scheme keeps it non-hyperbolic.
    # E(A) is formed at steady state, so E·A* has no twin and no placement is
    # redundant: all four are emitted. At {E, E(Q)}, where E(Q)·A* duplicates the
    # abortive E(A, Q), the copy adds an identifiable constant (11 fitted, rank 9, over
    # the parent's 10 and 8); at E alone it adds a phantom (11 fitted, rank 8).
    dead_end_with_a = EnzymeRates._expand_add_dead_end_regulator(dead_end, rxn_a_inhibits)
    @test length(dead_end_with_a) == 4
    @test Set(copy_sites.(dead_end_with_a)) ==
        Set([Set([:E]), Set([:E, :EQ]), Set([:EA, :EQ]), Set([:E, :EA])])
    @test !any(EnzymeRates._hyperbolic_catalysis, dead_end_with_a)
    at_e_eq = only(k for k in dead_end_with_a if copy_sites(k) == Set([:E, :EQ]))
    @test _testhelper_identifiable_rank(dead_end) == 8 && _testhelper_fitted(dead_end) == 10
    @test _testhelper_identifiable_rank(at_e_eq) == 9 && _testhelper_fitted(at_e_eq) == 11
    @test _testhelper_identifiable_rank(at_e(dead_end_with_a)) == 8

    @test !EnzymeRates._requires_hyperbolic_catalysis(ordered)
    @test EnzymeRates._requires_hyperbolic_catalysis(allo_unibi_flip_p)
    # One catalytic subunit: the conformational equilibrium reweights each form
    # without raising any catalytic-site binding to a power, so the rule does
    # not apply.
    allo_unibi_flip_p_one_subunit = EnzymeRates.AllostericMechanism(
        @allosteric_mechanism begin
            substrates: S
            products: P, Q
            catalytic_multiplicity: 1
            catalytic_steps: begin
                E + S ⇌ E(S)                                        :: NonequalAI
                E(S) <--> E(P, Q)                                   :: NonequalAI
                (E + P <--> E(P), E(Q) + P <--> E(P, Q))            :: NonequalAI
                (E + Q ⇌ E(Q), E(P) + Q ⇌ E(P, Q))                  :: NonequalAI
            end
        end)
    @test !EnzymeRates._requires_hyperbolic_catalysis(allo_unibi_flip_p_one_subunit)
end

@testset "_all_reach" begin
    # Segment graph 1 ⇄ 2 ⇄ 3 with edge tuples (source, target, form, mets).
    edges = [(1, 2, 0, Symbol[]), (2, 1, 0, Symbol[]),
             (2, 3, 0, Symbol[]), (3, 2, 0, Symbol[])]
    # Unpinned, every segment reaches every root.
    @test EnzymeRates._all_reach(3, edges, 3, ())
    # Pinning 1 → 2 keeps the path 1 → 2 → 3.
    @test EnzymeRates._all_reach(3, edges, 3, (edges[1],))
    # Pinning 2 → 1 strands 2 (and 1) away from root 3.
    @test !EnzymeRates._all_reach(3, edges, 3, (edges[2],))
    # A pinned 2-cycle (1 → 2 and 2 → 1) reaches no root outside it.
    @test !EnzymeRates._all_reach(3, edges, 3, (edges[1], edges[2]))
    # A root that is a pinned edge's target is fine.
    @test EnzymeRates._all_reach(3, edges, 2, (edges[1], edges[4]))
end

@testset "_hyperbolic_catalysis matches the derived denominator" begin
    # The structural predicate against the exponents of the derived denominator,
    # over every mechanism reachable from the seeds in two expansion levels. The
    # ping-pong-capable reaction's sequential seeds are the bi-bi seeds, and the moves,
    # the predicate and the derivation read metabolite names, never atoms, so its
    # population holds every bi-bi mechanism and the allosteric children of the bi-bi
    # allosteric mechanisms. The reactions declare no inhibitors, so the predicate sees
    # every step and is compared with each mechanism's own derived denominator. For an
    # allosteric mechanism the predicate is compared with the A-state, and the I-state
    # is checked to be hyperbolic whenever the A-state is. Allosteric mechanisms with
    # the same state graph share one derivation per state.
    hyperbolic(p, mets) = all(e <= 1 for mono in keys(p) for (s, e) in mono if s in mets)
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
    n_checked = 0
    n_nonhyperbolic = 0
    for rxn in (unibi, pingpong)
        mets = Symbol[EnzymeRates.name(x) for x in
                      vcat(EnzymeRates.substrates(rxn), EnzymeRates.products(rxn))]
        level = Union{EnzymeRates.Mechanism, EnzymeRates.AllostericMechanism}[
            m for m in EnzymeRates.init_mechanisms(rxn)]
        mechs = copy(level)
        for _ in 1:2
            level = unique!(EnzymeRates.expand_mechanisms(level, rxn))
            append!(mechs, level)
        end
        unique!(mechs)
        derived = Dict{Tuple{Symbol, EnzymeRates.Mechanism}, Bool}()
        state_hyperbolic(m, state) = get!(
            derived, (state, EnzymeRates._state_mechanism(m, state))) do
            hyperbolic(EnzymeRates._state_rate_polys(m, state)[2], mets)
        end
        for m in mechs
            EnzymeRates._eq_complexity(m) <= 337 || continue
            structural = EnzymeRates._hyperbolic_catalysis(m)
            if m isa EnzymeRates.Mechanism
                _, den, _ = EnzymeRates._raw_symbolic_rate_polys(
                    m, EnzymeRates._step_parameters(m),
                    EnzymeRates._build_wegscheider_rename_map(m))
                @test structural == hyperbolic(den, mets)
            else
                hyp_a = state_hyperbolic(m, :A)
                @test structural == hyp_a
                @test state_hyperbolic(m, :I) || !hyp_a
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
    m = EnzymeRates.Mechanism(@enzyme_mechanism begin
        substrates: S
        products: P
        steps: begin
            E + S ⇌ E(S)
            E(S) <--> E(P)
            E + P ⇌ E(P)
        end
    end)
    groups = EnzymeRates.steps(m)
    @test EnzymeRates._re_segment_count(groups) == 1
    g = findfirst(grp -> all(EnzymeRates.is_equilibrium, grp), groups)
    flipped = [i == g ? EnzymeRates._with_equilibrium.(grp, false) : grp
               for (i, grp) in enumerate(groups)]
    @test EnzymeRates._re_segment_count(flipped) == 2

    # Allosteric: measured on its catalytic groups.
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
    am_groups = EnzymeRates.steps(am)
    @test EnzymeRates._re_segment_count(am_groups) == 1
    am_g = findfirst(grp -> all(EnzymeRates.is_equilibrium, grp), am_groups)
    am_flipped = [i == am_g ? EnzymeRates._with_equilibrium.(grp, false) : grp
                  for (i, grp) in enumerate(am_groups)]
    @test EnzymeRates._re_segment_count(am_flipped) == 2

    # Aggregate cross-check over the bi-bi seed set: on every single RE→SS flip the count
    # agrees with the segmentation `_re_segment_extras` produces.
    flips = [[i == g ? EnzymeRates._with_equilibrium.(grp, false) : grp
              for (i, grp) in enumerate(EnzymeRates.steps(seed))]
             for seed in EnzymeRates.init_mechanisms(bi_bi_rxn)
             for g in eachindex(EnzymeRates.steps(seed))
             if all(EnzymeRates.is_equilibrium, EnzymeRates.steps(seed)[g])]
    @test length(flips) == 508
    @test all(gs -> EnzymeRates._re_segment_count(gs) ==
                    length(EnzymeRates._re_segment_extras(gs)[2]), flips)
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

"Whole-group flip of groups `gs` of `m` to steady state (test helper)."
function _testhelper_flip_groups(m, gs)
    groups = [g in gs ? EnzymeRates._with_equilibrium.(grp, false) : grp
              for (g, grp) in enumerate(EnzymeRates.steps(m))]
    EnzymeRates._with(m; groups)
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
    # _context_form: an SS dissociation step whose released metabolite is in
    # neither endpoint's bound list (the Segel ping-pong step shape) is stored
    # as the binding it reverses, F + P → E(A), so its context form is F, the
    # form P binds to.
    ping_pong_step = EnzymeRates.Step(
        EnzymeRates.Species([EnzymeRates.Substrate(:A)], :E),
        EnzymeRates.Species(EnzymeRates.Metabolite[], :F),
        EnzymeRates.Metabolite[], [EnzymeRates.Product(:P)], false)
    @test EnzymeRates._context_form(ping_pong_step) ==
          EnzymeRates.Species(EnzymeRates.Metabolite[], :F) ==
          EnzymeRates.from_species(ping_pong_step)
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
    @test sum(length, EnzymeRates.steps(child)) == sum(length, EnzymeRates.steps(m))
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
    for (gi, grp) in enumerate(EnzymeRates.steps(achild))
        Set(grp) ⊆ Set(bpa[1]) || Set(grp) ⊆ Set(bpa[2]) || continue
        @test EnzymeRates.cat_allo_state(achild, gi) == :NonequalAI
    end
    @test EnzymeRates.catalytic_multiplicity(achild) == 2
    @test EnzymeRates.regulatory_sites(achild) == EnzymeRates.regulatory_sites(am)
end

@testset "_revert_zero_flux_parts" begin
    # Random-order product release with B also bound abortively at E(Q) and at E(P),
    # the two abortive bindings one steady-state group. Neither abortive step carries
    # flux, so with both flags false both parts of the group's context bipartition
    # come back at rapid equilibrium; with both flags true the bipartition is returned
    # as it is.
    m = EnzymeRates.Mechanism(@enzyme_mechanism begin
        substrates: A, B
        products: P, Q
        steps: begin
            E + A ⇌ E(A)
            E(A) + B ⇌ E(A, B)
            E(A, B) <--> E(P, Q)
            E(Q) + P ⇌ E(P, Q)
            E(P) + Q ⇌ E(P, Q)
            E + P ⇌ E(P)
            E + Q ⇌ E(Q)
            (E(Q) + B <--> E(B, Q), E(P) + B <--> E(B, P))
        end
    end)
    groups = EnzymeRates.steps(m)
    g = only(findall(grp -> length(grp) == 2, groups))
    group = groups[g]
    @test !any(EnzymeRates.is_equilibrium, group)
    @test EnzymeRates._flux_carrying_steps(groups, EnzymeRates.reaction(m))[g] == falses(2)
    bp = only(EnzymeRates._context_bipartitions(group))
    parts = EnzymeRates._revert_zero_flux_parts(group, bp, falses(2))
    for (part, raw) in zip(parts, bp)
        @test length(part) == length(raw) == 1
        for (s, r) in zip(part, raw)
            @test EnzymeRates.is_equilibrium(s)
            @test EnzymeRates.from_species(s) == EnzymeRates.from_species(r)
            @test EnzymeRates.to_species(s) == EnzymeRates.to_species(r)
            @test EnzymeRates.consumed(s) == EnzymeRates.consumed(r)
            @test EnzymeRates.released(s) == EnzymeRates.released(r)
        end
    end
    @test EnzymeRates._revert_zero_flux_parts(group, bp, trues(2)) == bp
end

@testset "_partition_independent_count agrees with _independent_param_count" begin
    seeds = EnzymeRates.init_mechanisms(bi_bi_rxn)
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

@testset "_partition_independent_count counts a reverted part as its RE constant" begin
    # The E(Q) + B step of the steady-state B group, relabelled into a new group
    # and counted as a binding K, gives the count of the built child with that
    # step at rapid equilibrium: the cycle basis does not depend on the flags.
    m = EnzymeRates.Mechanism(@enzyme_mechanism begin
        substrates: A, B
        products: P, Q
        steps: begin
            E + A ⇌ E(A)
            (E(A) + B <--> E(A, B), E(Q) + B <--> E(B, Q))
            E(A, B) <--> E(P, Q)
            E(Q) + P ⇌ E(P, Q)
            E + Q ⇌ E(Q)
        end
    end)
    child = EnzymeRates.Mechanism(@enzyme_mechanism begin
        substrates: A, B
        products: P, Q
        steps: begin
            E + A ⇌ E(A)
            E(A) + B <--> E(A, B)
            E(Q) + B ⇌ E(B, Q)
            E(A, B) <--> E(P, Q)
            E(Q) + P ⇌ E(P, Q)
            E + Q ⇌ E(Q)
        end
    end)
    counter = EnzymeRates._partition_independent_count(m)
    flat = EnzymeRates._flat_steps(m)
    ids = [g for (_, g) in flat]
    kinds = [EnzymeRates._count_kind(s) for (s, _) in flat]
    j = only(j for (j, (s, _)) in enumerate(flat)
             if EnzymeRates.name(EnzymeRates.from_species(s)) == :EQ &&
                EnzymeRates.bound_metabolite(s) !== nothing &&
                EnzymeRates.name(EnzymeRates.bound_metabolite(s)) == :B)
    @test kinds[j] == :ss
    base = counter(ids)
    ids[j] = length(EnzymeRates.steps(m)) + 1
    kinds[j] = :binding_K
    @test counter(ids, kinds) == EnzymeRates._independent_param_count(child)
    @test counter(ids, kinds) == base + 1
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
        init = EnzymeRates.init_mechanisms(bi_bi_rxn)
        seeds = filter(_testhelper_holds_iso, init)
        variants = filter(!_testhelper_holds_iso, init)
        @test (length(seeds), length(variants)) == (55, 184)
        @test count(m -> !isempty(EnzymeRates._chain_flank_groups(m)), init) == 8
        @test sum(length(EnzymeRates._expand_re_to_ss(m)) for m in seeds) == 204
        @test sum(length(EnzymeRates._expand_re_to_ss(m)) for m in variants) == 272
    end

    @testset "_expand_re_to_ss: uni-uni emits no flip, both bindings are chain flanks" begin
        # Uni-uni steady-state and rapid-equilibrium laws have the same form: E(S) →
        # E(P) is a qualifying chain whose flanks are the S and P groups, so neither
        # flips (`_chain_flank_groups`). The S flip, absent, keeps the parent's rank; the
        # P flip is its mirror under reversing the reaction.
        m = EnzymeRates.Mechanism(@enzyme_mechanism begin
            substrates: S
            products: P
            steps: begin
                E + S ⇌ E(S)
                E(S) <--> E(P)
                E + P ⇌ E(P)
            end
        end)
        absent = EnzymeRates.Mechanism(@enzyme_mechanism begin
            substrates: S
            products: P
            steps: begin
                E + S <--> E(S)
                E(S) <--> E(P)
                E + P ⇌ E(P)
            end
        end)
        kids = EnzymeRates._expand_re_to_ss(m)
        @test EnzymeRates._chain_flank_groups(m) ==
              Set(g for g in EnzymeRates.kinetic_groups(m)
                  if !EnzymeRates.is_iso(only(m.steps[g])))
        @test isempty(kids)
        @test _testhelper_identifiable_rank(absent) == _testhelper_identifiable_rank(m)
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
            @test EnzymeRates._re_segment_count(EnzymeRates.steps(single)) ==
                  EnzymeRates._re_segment_count(EnzymeRates.steps(m))
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
        # joined through its sibling, so nothing flips one part by itself; C, P, Q
        # and R still flip alone. Of the six pairs, {A1, B1} cuts both binding edges
        # at E and E(C) and leaves {E(A), E(B), E(A, B), E(A, C), E(B, C), E(A, B, C)}
        # as a segment with no substrate-free form, which the constructor rejects.
        # {A1, B2} isolates {E(A), E(A, C)}: its four edges to the rest all have
        # weight 1 inward (A uptake) and −1 outward (B uptake against the +2 offset
        # of E(A, B)), so every cycle through them is balanced and the pair carries
        # no flux; the turnover runs through E(B) inside the big segment. {A2, B1}
        # isolates {E(B), E(B, C)} the same way. Both count as failed and have no
        # extension, since every three-set containing them holds a gaining pair.
        # The other three pairs cut a block that also holds the chemistry edge
        # (weight 4 or 5 against flipped edges of weight 2 or 1; every cycle through
        # it sums to 6) and are emitted.
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
            flip(A2, B2),   # {A2, B2}
        ]
        kids = EnzymeRates._expand_re_to_ss(m)
        @test length(kids) == 7
        @test Set(kids) == Set(expected)
        for c in kids
            ss = count((A1, A2, B1, B2)) do key
                grp = only(grp for grp in EnzymeRates.steps(c)
                           if (binds(grp), sources(grp)) == key)
                !any(EnzymeRates.is_equilibrium, grp)
            end
            @test ss != 1
        end
        # The two zero-flux pairs divide a segment yet keep the parent's rank;
        # flip(A2, B1) is flip(A1, B2) with A and B swapped.
        for absent in (flip(A1, B2), flip(A2, B1))
            @test !(absent in kids)
            @test EnzymeRates._re_segment_count(EnzymeRates.steps(absent)) >
                  EnzymeRates._re_segment_count(EnzymeRates.steps(m))
        end
        @test _testhelper_identifiable_rank(flip(A1, B2)) ==
            _testhelper_identifiable_rank(m)
    end

    @testset "_expand_re_to_ss: ping-pong" begin
        # The rapid-equilibrium subgraph is two trees, {E, E(A), E(Q)} and
        # {E(P; res), E(; res), E(B; res)}, so every RE group is a bridge and
        # flipping any one of them alone raises the segment count. Both
        # isomerizations are steady state, and each of their ends has one other
        # step, a binding into that end alone in its group: E(A) → E(P; res) and
        # E(B; res) → E(Q) are qualifying chains whose flanks are all four binding
        # groups (`_chain_flank_groups`). No group is a unit, so no child. The A flip,
        # absent, keeps the parent's rank; the B, P and Q flips are its images under
        # swapping the two half-reactions and reversing the reaction.
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
        kids = EnzymeRates._expand_re_to_ss(m)
        @test isempty(kids)
        @test _testhelper_identifiable_rank(flipA) == _testhelper_identifiable_rank(m)
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

    @testset "_expand_re_to_ss: a one-subunit allosteric parent keeps every flip" begin
        # The same uni-bi as an allosteric mechanism with one catalytic subunit.
        # The conformational equilibrium adds no concentration power, so the
        # flip move keeps the P and Q flips it drops for two subunits, and emits
        # the same three flips as the plain mechanism.
        allo_one = EnzymeRates.AllostericMechanism(@allosteric_mechanism begin
            substrates: S
            products: P, Q
            catalytic_multiplicity: 1
            catalytic_steps: begin
                E + S ⇌ E(S)                                :: NonequalAI
                E(S) <--> E(P, Q)                           :: NonequalAI
                (E + P ⇌ E(P), E(Q) + P ⇌ E(P, Q))          :: NonequalAI
                (E + Q ⇌ E(Q), E(P) + Q ⇌ E(P, Q))          :: NonequalAI
            end
        end)
        flipS_one = EnzymeRates.AllostericMechanism(@allosteric_mechanism begin
            substrates: S
            products: P, Q
            catalytic_multiplicity: 1
            catalytic_steps: begin
                E + S <--> E(S)                             :: NonequalAI
                E(S) <--> E(P, Q)                           :: NonequalAI
                (E + P ⇌ E(P), E(Q) + P ⇌ E(P, Q))          :: NonequalAI
                (E + Q ⇌ E(Q), E(P) + Q ⇌ E(P, Q))          :: NonequalAI
            end
        end)
        flipP_one = EnzymeRates.AllostericMechanism(@allosteric_mechanism begin
            substrates: S
            products: P, Q
            catalytic_multiplicity: 1
            catalytic_steps: begin
                E + S ⇌ E(S)                                :: NonequalAI
                E(S) <--> E(P, Q)                           :: NonequalAI
                (E + P <--> E(P), E(Q) + P <--> E(P, Q))    :: NonequalAI
                (E + Q ⇌ E(Q), E(P) + Q ⇌ E(P, Q))          :: NonequalAI
            end
        end)
        flipQ_one = EnzymeRates.AllostericMechanism(@allosteric_mechanism begin
            substrates: S
            products: P, Q
            catalytic_multiplicity: 1
            catalytic_steps: begin
                E + S ⇌ E(S)                                :: NonequalAI
                E(S) <--> E(P, Q)                           :: NonequalAI
                (E + P ⇌ E(P), E(Q) + P ⇌ E(P, Q))          :: NonequalAI
                (E + Q <--> E(Q), E(P) + Q <--> E(P, Q))    :: NonequalAI
            end
        end)
        kids_one = EnzymeRates._expand_re_to_ss(allo_one)
        @test length(kids_one) == 3
        @test Set(kids_one) == Set([flipS_one, flipP_one, flipQ_one])
    end

    @testset "_expand_re_to_ss: a zero-flux pair is absent, its supersets reachable" begin
        # Ordered bi-bi with the abortive complex E(A, Q), after both the A group and
        # the Q group were split by context: six rapid-equilibrium binding groups,
        # each holding one step. Singles: E(A) + B and E(Q) + P are bridges of the RE
        # graph (E(A, B) and E(P, Q) are leaves there), so each would flip alone, but
        # they are the flanks of the qualifying chain E(A, B) → E(P, Q) and never flip
        # (`_chain_flank_groups`); the four bindings on the square E–E(A)–E(A, Q)–E(Q)
        # are no bridges. Pairs of square edges all divide the segment.
        # {E + A, E + Q} isolates {E}, joined to the rest by
        # two parallel edges of weight 0 (A uptake +1 against E(A)'s offset +1; Q
        # uptake −1 against E(Q)'s offset −1): a balanced block, so neither flipped
        # group carries flux and the child would be its parent plus two phantoms.
        # {E(A) + Q, E(Q) + A} isolates {E(A, Q)} the same way (weights −1 + 1 and
        # 1 − 1). The other four pairs each cut a block that also holds the chemistry
        # edge, whose weight 3 or 4 against a binding edge of weight 1 makes the block
        # unbalanced. Each of the two zero-flux pairs has no extension: every
        # three-set containing it also contains a gaining pair, so its supersets are
        # reached by a later flip of those children.
        m = EnzymeRates.Mechanism(@enzyme_mechanism begin
            substrates: A, B
            products: P, Q
            steps: begin
                E + A ⇌ E(A)
                E(Q) + A ⇌ E(A, Q)
                E + Q ⇌ E(Q)
                E(A) + Q ⇌ E(A, Q)
                E(A) + B ⇌ E(A, B)
                E(A, B) <--> E(P, Q)
                E(Q) + P ⇌ E(P, Q)
            end
        end)
        expected = [
            EnzymeRates.Mechanism(@enzyme_mechanism begin      # {E + A, E(A) + Q}
                substrates: A, B
                products: P, Q
                steps: begin
                    E + A <--> E(A)
                    E(Q) + A ⇌ E(A, Q)
                    E + Q ⇌ E(Q)
                    E(A) + Q <--> E(A, Q)
                    E(A) + B ⇌ E(A, B)
                    E(A, B) <--> E(P, Q)
                    E(Q) + P ⇌ E(P, Q)
                end
            end),
            EnzymeRates.Mechanism(@enzyme_mechanism begin      # {E + A, E(Q) + A}
                substrates: A, B
                products: P, Q
                steps: begin
                    E + A <--> E(A)
                    E(Q) + A <--> E(A, Q)
                    E + Q ⇌ E(Q)
                    E(A) + Q ⇌ E(A, Q)
                    E(A) + B ⇌ E(A, B)
                    E(A, B) <--> E(P, Q)
                    E(Q) + P ⇌ E(P, Q)
                end
            end),
            EnzymeRates.Mechanism(@enzyme_mechanism begin      # {E + Q, E(A) + Q}
                substrates: A, B
                products: P, Q
                steps: begin
                    E + A ⇌ E(A)
                    E(Q) + A ⇌ E(A, Q)
                    E + Q <--> E(Q)
                    E(A) + Q <--> E(A, Q)
                    E(A) + B ⇌ E(A, B)
                    E(A, B) <--> E(P, Q)
                    E(Q) + P ⇌ E(P, Q)
                end
            end),
            EnzymeRates.Mechanism(@enzyme_mechanism begin      # {E + Q, E(Q) + A}
                substrates: A, B
                products: P, Q
                steps: begin
                    E + A ⇌ E(A)
                    E(Q) + A <--> E(A, Q)
                    E + Q <--> E(Q)
                    E(A) + Q ⇌ E(A, Q)
                    E(A) + B ⇌ E(A, B)
                    E(A, B) <--> E(P, Q)
                    E(Q) + P ⇌ E(P, Q)
                end
            end),
        ]
        kids = EnzymeRates._expand_re_to_ss(m)
        @test length(kids) == 4
        @test Set(kids) == Set(expected)
        # The flank flip of B and the two zero-flux pairs are absent; each divides a
        # segment, and each has exactly its parent's rank. The flank flip of P is the
        # flip of B with the reaction reversed.
        r0 = _testhelper_identifiable_rank(m)
        for absent in (
            EnzymeRates.Mechanism(@enzyme_mechanism begin  # {E(A) + B}
                substrates: A, B
                products: P, Q
                steps: begin
                    E + A ⇌ E(A)
                    E(Q) + A ⇌ E(A, Q)
                    E + Q ⇌ E(Q)
                    E(A) + Q ⇌ E(A, Q)
                    E(A) + B <--> E(A, B)
                    E(A, B) <--> E(P, Q)
                    E(Q) + P ⇌ E(P, Q)
                end
            end),
            EnzymeRates.Mechanism(@enzyme_mechanism begin  # {E + A, E + Q}
                substrates: A, B
                products: P, Q
                steps: begin
                    E + A <--> E(A)
                    E(Q) + A ⇌ E(A, Q)
                    E + Q <--> E(Q)
                    E(A) + Q ⇌ E(A, Q)
                    E(A) + B ⇌ E(A, B)
                    E(A, B) <--> E(P, Q)
                    E(Q) + P ⇌ E(P, Q)
                end
            end),
            EnzymeRates.Mechanism(@enzyme_mechanism begin  # {E(A) + Q, E(Q) + A}
                substrates: A, B
                products: P, Q
                steps: begin
                    E + A ⇌ E(A)
                    E(Q) + A <--> E(A, Q)
                    E + Q ⇌ E(Q)
                    E(A) + Q <--> E(A, Q)
                    E(A) + B ⇌ E(A, B)
                    E(A, B) <--> E(P, Q)
                    E(Q) + P ⇌ E(P, Q)
                end
            end))
            @test !(absent in kids)
            @test EnzymeRates._re_segment_count(EnzymeRates.steps(absent)) >
                  EnzymeRates._re_segment_count(EnzymeRates.steps(m))
            @test _testhelper_identifiable_rank(absent) == r0
        end
    end

    @testset "_expand_re_to_ss: a Theorell–Chance parent flips each binding" begin
        # No isomerization step. Every step lies on the unbalanced cycle
        # E → E(A) → E(Q) → E (weights 1, 2, 1), so both rapid-equilibrium bindings
        # are units; each divides the one segment alone.
        m = EnzymeRates.Mechanism(@enzyme_mechanism begin
            substrates: A, B
            products: P, Q
            steps: begin
                E + A ⇌ E(A)
                E(A) + B <--> E(Q) + P
                E + Q ⇌ E(Q)
            end
        end)
        expected = [
            EnzymeRates.Mechanism(@enzyme_mechanism begin
                substrates: A, B
                products: P, Q
                steps: begin
                    E + A <--> E(A)
                    E(A) + B <--> E(Q) + P
                    E + Q ⇌ E(Q)
                end
            end),
            EnzymeRates.Mechanism(@enzyme_mechanism begin
                substrates: A, B
                products: P, Q
                steps: begin
                    E + A ⇌ E(A)
                    E(A) + B <--> E(Q) + P
                    E + Q <--> E(Q)
                end
            end),
        ]
        kids = EnzymeRates._expand_re_to_ss(m)
        @test length(kids) == 2
        @test Set(kids) == Set(expected)
    end
end

@testset "_expand_split_kinetic_group (context bipartitions)" begin
    @testset "_expand_split_kinetic_group: random-order bi-bi frees indep params" begin
        # The split closure of this seed, which has 5 independent parameters, holds 16
        # structures and reaches the form with every step in its own group and 9
        # independent parameters.
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
        # Aggregate regression pin over the whole bi-bi init set: 102 children from the
        # 55 seeds and 420 from their 184 merged and Theorell–Chance variants, which are
        # 176 from the 108 two-group merged, 236 from the 56 three-group merged and 8
        # from the 20 Theorell–Chance ones.
        init = EnzymeRates.init_mechanisms(bi_bi_rxn)
        n_children(ms) = sum(length(EnzymeRates._expand_split_kinetic_group(m)) for m in ms)
        @test n_children(filter(_testhelper_holds_iso, init)) == 102
        @test n_children(filter(!_testhelper_holds_iso, init)) == 420
    end

    @testset "_expand_split_kinetic_group: a rejected split is the parent's model" begin
        # Level-1 candidates the count test rejects have the parent's identifiable
        # rank; the rendered equation may differ only in which tied name survives. The
        # B, P and Q groups' candidates are images of the A group's under relabeling
        # (A and B swapped with P and Q swapped) and reversing the reaction, so the
        # rank is checked on the A group's.
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
        g = findfirst(grp -> EnzymeRates.bound_metabolite(first(grp)) ==
                             EnzymeRates.Substrate(:A), EnzymeRates.steps(m))
        for bp in EnzymeRates._context_bipartitions(EnzymeRates.steps(m)[g])
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

@testset "expansion moves: ter-ter random-order seed within budget" begin
    # The seed of `init_mechanisms(terter)` with the most steps (55) is the enumeration's
    # worst case: random order, each substrate and product binding at nine forms in one
    # rapid-equilibrium group. The split alone must take under 60 s and the four moves
    # that apply to a plain Mechanism (split, flip, dead end, to-allosteric) under 120 s
    # together; they return 81 children: 12 splits, 6 flips and 63 K-types (the reaction
    # declares no regulator). Measured 26 s for the four moves in a cold focused run, JIT
    # included.
    terter = @enzyme_reaction begin
        substrates: A[C], B[N], C[O]
        products: P[C], Q[N], R[O]
    end
    seed = @enzyme_mechanism begin
        substrates: A, B, C
        products: P, Q, R
        steps: begin
            (E + A ⇌ E(A), E(B) + A ⇌ E(A, B), E(B, C) + A ⇌ E(A, B, C),
             E(B, P) + A ⇌ E(A, B, P), E(C) + A ⇌ E(A, C), E(C, Q) + A ⇌ E(A, C, Q),
             E(P) + A ⇌ E(A, P), E(P, Q) + A ⇌ E(A, P, Q), E(Q) + A ⇌ E(A, Q))
            (E + B ⇌ E(B), E(A) + B ⇌ E(A, B), E(A, C) + B ⇌ E(A, B, C),
             E(A, P) + B ⇌ E(A, B, P), E(C) + B ⇌ E(B, C), E(C, R) + B ⇌ E(B, C, R),
             E(P) + B ⇌ E(B, P), E(P, R) + B ⇌ E(B, P, R), E(R) + B ⇌ E(B, R))
            (E + C ⇌ E(C), E(A) + C ⇌ E(A, C), E(A, B) + C ⇌ E(A, B, C),
             E(A, Q) + C ⇌ E(A, C, Q), E(B) + C ⇌ E(B, C), E(B, R) + C ⇌ E(B, C, R),
             E(Q) + C ⇌ E(C, Q), E(Q, R) + C ⇌ E(C, Q, R), E(R) + C ⇌ E(C, R))
            (E + P ⇌ E(P), E(A) + P ⇌ E(A, P), E(A, B) + P ⇌ E(A, B, P),
             E(A, Q) + P ⇌ E(A, P, Q), E(B) + P ⇌ E(B, P), E(B, R) + P ⇌ E(B, P, R),
             E(Q) + P ⇌ E(P, Q), E(Q, R) + P ⇌ E(P, Q, R), E(R) + P ⇌ E(P, R))
            (E + Q ⇌ E(Q), E(A) + Q ⇌ E(A, Q), E(A, C) + Q ⇌ E(A, C, Q),
             E(A, P) + Q ⇌ E(A, P, Q), E(C) + Q ⇌ E(C, Q), E(C, R) + Q ⇌ E(C, Q, R),
             E(P) + Q ⇌ E(P, Q), E(P, R) + Q ⇌ E(P, Q, R), E(R) + Q ⇌ E(Q, R))
            (E + R ⇌ E(R), E(B) + R ⇌ E(B, R), E(B, C) + R ⇌ E(B, C, R),
             E(B, P) + R ⇌ E(B, P, R), E(C) + R ⇌ E(C, R), E(C, Q) + R ⇌ E(C, Q, R),
             E(P) + R ⇌ E(P, R), E(P, Q) + R ⇌ E(P, Q, R), E(Q) + R ⇌ E(Q, R))
            E(A, B, C) <--> E(P, Q, R)
        end
    end
    worst = EnzymeRates.Mechanism(terter, EnzymeRates.steps(EnzymeRates.Mechanism(seed)))
    @test sum(length, EnzymeRates.steps(worst)) == 55
    t_split = @elapsed split = EnzymeRates._expand_split_kinetic_group(worst)
    t_rest = @elapsed rest = vcat(
        EnzymeRates._expand_re_to_ss(worst),
        EnzymeRates._expand_add_dead_end_regulator(worst, terter),
        EnzymeRates._expand_to_allosteric(worst, terter))
    @test length(split) == 12
    @test t_split < 60
    @test length(rest) == 69
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
        dead_end_inhibitors: A, B, P, Q
    end
    moves(m, rxn) = vcat(EnzymeRates._expand_re_to_ss(m),
                         EnzymeRates._expand_split_kinetic_group(m),
                         EnzymeRates._expand_add_dead_end_regulator(m, rxn))
    iso_kept = true
    function levels(rxn)
        level = unique!(EnzymeRates.init_mechanisms(rxn))
        seen = Set(level)
        out = [level]
        for _ in 1:2
            next = EnzymeRates.Mechanism[]
            for m in level, c in moves(m, rxn)
                iso_kept &= _testhelper_holds_iso(c) == _testhelper_holds_iso(m)
                c in seen && continue
                push!(seen, c); push!(next, c)
            end
            level = next
            push!(out, level)
        end
        out
    end
    obeys_rules(m) = EnzymeRates._assert_emission_rules(m) === nothing
    free(m) = isempty(EnzymeRates._bound_comp_inhibitors(m))
    copies = levels(rxn)
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
    # apart from the seeds: 5,109 mechanisms at level 1 and 97,111 at level 2. Grown
    # from all 264 init mechanisms the levels hold the sums, 264, 6,738 and 128,493.
    @test [count(_testhelper_holds_iso, l) for l in copies] == [62, 1629, 31382]
    @test length.(copies) == [264, 6738, 128493]

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
    # level 1 and 1,237 at level 2. Grown from all 264 init mechanisms the levels
    # hold the sums, 264, 1,018 and 2,371.
    @test [count(m -> free(m) && _testhelper_holds_iso(m), l) for l in copies] ==
          [62, 349, 1134]
    @test [count(m -> free(m) && !_testhelper_holds_iso(m), l) for l in copies] ==
          [202, 669, 1237]
    @test [count(free, l) for l in copies] == [264, 1018, 2371]
end

@testset "expand_mechanisms on a merged uni-uni" begin
    # E + S → E(P) fused (SS), E + P ⇌ E(P) (RE). Flip: the P group cuts E from E(P),
    # carries flux, no bottomless segment: one child, both steps SS. Split: single-step
    # groups, nothing. Dead end: no regulator. To-allosteric at multiplicity 1: the
    # chemistry group (the fused step) is :OnlyA, the binding subsets range over {P}:
    # one K-type child; no regulator, so no V-type.
    rxn = @enzyme_reaction begin
        substrates: S[C]
        products: P[C]
    end
    m = _testhelper_on_reaction(rxn, @enzyme_mechanism begin
        substrates: S
        products: P
        steps: begin
            E + S <--> E(P)
            E + P ⇌ E(P)
        end
    end)
    flipped = _testhelper_on_reaction(rxn, @enzyme_mechanism begin
        substrates: S
        products: P
        steps: begin
            E + S <--> E(P)
            E + P <--> E(P)
        end
    end)
    k_type = EnzymeRates.AllostericMechanism(rxn, EnzymeRates.steps(m),
        [:OnlyA, :OnlyA], 1, EnzymeRates.RegulatorySite[])
    kids = EnzymeRates.expand_mechanisms([m], rxn)
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
end

@testset "_expand_change_allo_state on a merged uni-uni K-type" begin
    # The fused step E + S → E(P) is the chemistry group. Relaxing P to :NonequalAI
    # keeps the chemistry :OnlyA: the inactive conformation binds P with a constant of
    # its own and runs no chemistry, a valid child. The chemistry is not relaxed while
    # P's binding is :OnlyA: that would restore inactive catalysis beside it, a partial
    # :OnlyA catalysis whose one-sided :OnlyA binding also leaves the Haldane cycle
    # unsatisfiable once the chemistry is back in the check graph
    # (`_onlya_haldane_violation`). One child.
    k_type = EnzymeRates.AllostericMechanism(@allosteric_mechanism begin
        substrates: S
        products: P
        catalytic_multiplicity: 1
        catalytic_steps: begin
            E + S <--> E(P)    :: OnlyA
            E + P ⇌ E(P)       :: OnlyA
        end
    end)
    p_relaxed = EnzymeRates.AllostericMechanism(@allosteric_mechanism begin
        substrates: S
        products: P
        catalytic_multiplicity: 1
        catalytic_steps: begin
            E + S <--> E(P)    :: OnlyA
            E + P ⇌ E(P)       :: NonequalAI
        end
    end)
    kids = EnzymeRates._expand_change_allo_state(k_type)
    @test length(kids) == 1
    @test Set(kids) == Set([p_relaxed])
    # The child's rate is k·(S − P/Keq)/(1 + L + P·(1/Kp_A + L/Kp_I)): dividing by
    # 1 + L leaves k/(1 + L) and one P coefficient, so k, Kp_A, Kp_I and L give
    # 4 fitted, rank 2.
    @test _testhelper_fitted(p_relaxed) == 4 &&
          _testhelper_identifiable_rank(p_relaxed) == 2
end

@testset "_expand_to_allosteric on an ordered Theorell–Chance mechanism" begin
    # E(A) + B → E(Q) + P takes up B and gives off P in one step: the chemistry group.
    # Each K-type child sets it :OnlyA and ranges the binding subsets over {A, Q}:
    # {A}, {Q} and {A, Q}. `_onlya_haldane_violation` drops the :OnlyA chemistry group
    # from its cycle graph; the A and Q bindings left form the tree E(A) – E – E(Q),
    # so no cycle row remains and no subset is refused. No regulator, so no V-type:
    # three children.
    m = EnzymeRates.Mechanism(@enzyme_mechanism begin
        substrates: A, B
        products: P, Q
        steps: begin
            E + A ⇌ E(A)
            E(A) + B <--> E(Q) + P
            E + Q ⇌ E(Q)
        end
    end)
    rxn = @enzyme_reaction begin
        substrates: A[C], B[N]
        products: P[C], Q[N]
    end
    only_a = EnzymeRates.AllostericMechanism(@allosteric_mechanism begin
        substrates: A, B ; products: P, Q ; catalytic_multiplicity: 1
        catalytic_steps: begin
            E + A ⇌ E(A)                :: OnlyA
            E(A) + B <--> E(Q) + P      :: OnlyA
            E + Q ⇌ E(Q)                :: EqualAI
        end
    end)
    only_q = EnzymeRates.AllostericMechanism(@allosteric_mechanism begin
        substrates: A, B ; products: P, Q ; catalytic_multiplicity: 1
        catalytic_steps: begin
            E + A ⇌ E(A)                :: EqualAI
            E(A) + B <--> E(Q) + P      :: OnlyA
            E + Q ⇌ E(Q)                :: OnlyA
        end
    end)
    a_and_q = EnzymeRates.AllostericMechanism(@allosteric_mechanism begin
        substrates: A, B ; products: P, Q ; catalytic_multiplicity: 1
        catalytic_steps: begin
            E + A ⇌ E(A)                :: OnlyA
            E(A) + B <--> E(Q) + P      :: OnlyA
            E + Q ⇌ E(Q)                :: OnlyA
        end
    end)
    kids = EnzymeRates._expand_to_allosteric(m, rxn)
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
end

@testset "_expand_to_allosteric on a merged ordered bi-bi at multiplicity 2" begin
    # The ordered bi-bi seed with E(A, B) merged into E(P, Q): B's binding is fused
    # with the chemistry (E(A) + B ⇌ E(P, Q)), and the A and P bindings are steady
    # state. The fused B group is the chemistry group, :OnlyA in every K-type child,
    # and the subsets range over the A, P and Q groups: 2^3 − 1 = 7. Without the
    # :OnlyA chemistry group the bindings form the tree E(A) – E – E(Q) – E(P, Q), so
    # `_onlya_haldane_violation` refuses none. The scheme is hyperbolic, so
    # multiplicity 2 stays open: the RE segments are {E, E(Q)}, where E(Q) carries one
    # Q, and {E(A), E(P, Q)}, where E(P, Q) carries one B; the steady-state A and P
    # bindings join them, and each metabolite scores on one edge only, an edge that
    # leaves the segment carrying that metabolite (`_hyperbolic_catalysis`). No
    # regulator, so no V-type: seven children.
    m2 = EnzymeRates.Mechanism(@enzyme_mechanism begin
        substrates: A, B
        products: P, Q
        steps: begin
            E + A <--> E(A)
            E(A) + B ⇌ E(P, Q)
            E(Q) + P <--> E(P, Q)
            E + Q ⇌ E(Q)
        end
    end)
    rxn2 = @enzyme_reaction begin
        substrates: A[C], B[N]
        products: P[C], Q[N]
        oligomeric_state: 2
    end
    only_a = EnzymeRates.AllostericMechanism(@allosteric_mechanism begin
        substrates: A, B ; products: P, Q ; catalytic_multiplicity: 2
        catalytic_steps: begin
            E + A <--> E(A)             :: OnlyA
            E(A) + B ⇌ E(P, Q)          :: OnlyA
            E(Q) + P <--> E(P, Q)       :: EqualAI
            E + Q ⇌ E(Q)                :: EqualAI
        end
    end)
    only_p = EnzymeRates.AllostericMechanism(@allosteric_mechanism begin
        substrates: A, B ; products: P, Q ; catalytic_multiplicity: 2
        catalytic_steps: begin
            E + A <--> E(A)             :: EqualAI
            E(A) + B ⇌ E(P, Q)          :: OnlyA
            E(Q) + P <--> E(P, Q)       :: OnlyA
            E + Q ⇌ E(Q)                :: EqualAI
        end
    end)
    only_q = EnzymeRates.AllostericMechanism(@allosteric_mechanism begin
        substrates: A, B ; products: P, Q ; catalytic_multiplicity: 2
        catalytic_steps: begin
            E + A <--> E(A)             :: EqualAI
            E(A) + B ⇌ E(P, Q)          :: OnlyA
            E(Q) + P <--> E(P, Q)       :: EqualAI
            E + Q ⇌ E(Q)                :: OnlyA
        end
    end)
    a_and_p = EnzymeRates.AllostericMechanism(@allosteric_mechanism begin
        substrates: A, B ; products: P, Q ; catalytic_multiplicity: 2
        catalytic_steps: begin
            E + A <--> E(A)             :: OnlyA
            E(A) + B ⇌ E(P, Q)          :: OnlyA
            E(Q) + P <--> E(P, Q)       :: OnlyA
            E + Q ⇌ E(Q)                :: EqualAI
        end
    end)
    a_and_q = EnzymeRates.AllostericMechanism(@allosteric_mechanism begin
        substrates: A, B ; products: P, Q ; catalytic_multiplicity: 2
        catalytic_steps: begin
            E + A <--> E(A)             :: OnlyA
            E(A) + B ⇌ E(P, Q)          :: OnlyA
            E(Q) + P <--> E(P, Q)       :: EqualAI
            E + Q ⇌ E(Q)                :: OnlyA
        end
    end)
    p_and_q = EnzymeRates.AllostericMechanism(@allosteric_mechanism begin
        substrates: A, B ; products: P, Q ; catalytic_multiplicity: 2
        catalytic_steps: begin
            E + A <--> E(A)             :: EqualAI
            E(A) + B ⇌ E(P, Q)          :: OnlyA
            E(Q) + P <--> E(P, Q)       :: OnlyA
            E + Q ⇌ E(Q)                :: OnlyA
        end
    end)
    all_three = EnzymeRates.AllostericMechanism(@allosteric_mechanism begin
        substrates: A, B ; products: P, Q ; catalytic_multiplicity: 2
        catalytic_steps: begin
            E + A <--> E(A)             :: OnlyA
            E(A) + B ⇌ E(P, Q)          :: OnlyA
            E(Q) + P <--> E(P, Q)       :: OnlyA
            E + Q ⇌ E(Q)                :: OnlyA
        end
    end)
    @test EnzymeRates._hyperbolic_catalysis(m2)
    kids = EnzymeRates._expand_to_allosteric(m2, rxn2)
    @test all(c -> c isa EnzymeRates.AllostericMechanism, kids)
    @test length(kids) == 7
    @test Set(kids) == Set([only_a, only_p, only_q, a_and_p, a_and_q, p_and_q,
                            all_three])
    # The parent has 5 fitted constants, all identifiable. Each child adds L; at
    # multiplicity 2 each conformation's binding polynomial enters squared, so dividing
    # by 1 + L no longer folds L into the other constants: 6 fitted, rank 6.
    @test _testhelper_fitted(m2) == _testhelper_identifiable_rank(m2) == 5
    for k in kids
        @test _testhelper_fitted(k) == _testhelper_identifiable_rank(k) == 6
    end
end

@testset "_expand_to_allosteric on a merged bi-bi whose B group mixes fused and plain" begin
    # The merged ordered bi-bi with B's abortive binding at E(Q) in B's own group:
    # E(A) + B ⇌ E(P, Q) is fused, E(Q) + B ⇌ E(B, Q) is plain, one RE kinetic group. A
    # group holding a chemistry step is chemistry, so the whole B group, its plain step
    # included, is :OnlyA in every K-type child, and the subsets range over the A, P
    # and Q groups: 2^3 − 1 = 7. Without the :OnlyA B group the check graph keeps the
    # tree E(A) – E – E(Q) – E(P, Q) (E(B, Q) leaves with B's group), so
    # `_onlya_haldane_violation` refuses none. No regulator, so no V-type: seven
    # children.
    m = EnzymeRates.Mechanism(@enzyme_mechanism begin
        substrates: A, B
        products: P, Q
        steps: begin
            E + A <--> E(A)
            (E(A) + B ⇌ E(P, Q), E(Q) + B ⇌ E(B, Q))
            E(Q) + P <--> E(P, Q)
            E + Q ⇌ E(Q)
        end
    end)
    rxn = @enzyme_reaction begin
        substrates: A[C], B[N]
        products: P[C], Q[N]
    end
    only_a = EnzymeRates.AllostericMechanism(@allosteric_mechanism begin
        substrates: A, B ; products: P, Q ; catalytic_multiplicity: 1
        catalytic_steps: begin
            E + A <--> E(A)                              :: OnlyA
            (E(A) + B ⇌ E(P, Q), E(Q) + B ⇌ E(B, Q))     :: OnlyA
            E(Q) + P <--> E(P, Q)                        :: EqualAI
            E + Q ⇌ E(Q)                                 :: EqualAI
        end
    end)
    only_p = EnzymeRates.AllostericMechanism(@allosteric_mechanism begin
        substrates: A, B ; products: P, Q ; catalytic_multiplicity: 1
        catalytic_steps: begin
            E + A <--> E(A)                              :: EqualAI
            (E(A) + B ⇌ E(P, Q), E(Q) + B ⇌ E(B, Q))     :: OnlyA
            E(Q) + P <--> E(P, Q)                        :: OnlyA
            E + Q ⇌ E(Q)                                 :: EqualAI
        end
    end)
    only_q = EnzymeRates.AllostericMechanism(@allosteric_mechanism begin
        substrates: A, B ; products: P, Q ; catalytic_multiplicity: 1
        catalytic_steps: begin
            E + A <--> E(A)                              :: EqualAI
            (E(A) + B ⇌ E(P, Q), E(Q) + B ⇌ E(B, Q))     :: OnlyA
            E(Q) + P <--> E(P, Q)                        :: EqualAI
            E + Q ⇌ E(Q)                                 :: OnlyA
        end
    end)
    a_and_p = EnzymeRates.AllostericMechanism(@allosteric_mechanism begin
        substrates: A, B ; products: P, Q ; catalytic_multiplicity: 1
        catalytic_steps: begin
            E + A <--> E(A)                              :: OnlyA
            (E(A) + B ⇌ E(P, Q), E(Q) + B ⇌ E(B, Q))     :: OnlyA
            E(Q) + P <--> E(P, Q)                        :: OnlyA
            E + Q ⇌ E(Q)                                 :: EqualAI
        end
    end)
    a_and_q = EnzymeRates.AllostericMechanism(@allosteric_mechanism begin
        substrates: A, B ; products: P, Q ; catalytic_multiplicity: 1
        catalytic_steps: begin
            E + A <--> E(A)                              :: OnlyA
            (E(A) + B ⇌ E(P, Q), E(Q) + B ⇌ E(B, Q))     :: OnlyA
            E(Q) + P <--> E(P, Q)                        :: EqualAI
            E + Q ⇌ E(Q)                                 :: OnlyA
        end
    end)
    p_and_q = EnzymeRates.AllostericMechanism(@allosteric_mechanism begin
        substrates: A, B ; products: P, Q ; catalytic_multiplicity: 1
        catalytic_steps: begin
            E + A <--> E(A)                              :: EqualAI
            (E(A) + B ⇌ E(P, Q), E(Q) + B ⇌ E(B, Q))     :: OnlyA
            E(Q) + P <--> E(P, Q)                        :: OnlyA
            E + Q ⇌ E(Q)                                 :: OnlyA
        end
    end)
    all_three = EnzymeRates.AllostericMechanism(@allosteric_mechanism begin
        substrates: A, B ; products: P, Q ; catalytic_multiplicity: 1
        catalytic_steps: begin
            E + A <--> E(A)                              :: OnlyA
            (E(A) + B ⇌ E(P, Q), E(Q) + B ⇌ E(B, Q))     :: OnlyA
            E(Q) + P <--> E(P, Q)                        :: OnlyA
            E + Q ⇌ E(Q)                                 :: OnlyA
        end
    end)
    kids = EnzymeRates._expand_to_allosteric(m, rxn)
    @test length(kids) == 7
    @test Set(kids) == Set([only_a, only_p, only_q, a_and_p, a_and_q, p_and_q,
                            all_three])
    # The parent has 5 fitted constants, all identifiable. Each child adds L, and its
    # dead inactive conformation keeps the bindings it reaches from the free enzyme
    # through :EqualAI groups: A's where A is :EqualAI; Q's where Q is, then P's at
    # E(Q) where P is too. Where it keeps one, L multiplies a binding polynomial that
    # carries A or Q, which no rescaling of the parent's constants absorbs: 6 fitted,
    # rank 6. Where A and Q are both :OnlyA it binds nothing, so L multiplies the free
    # enzyme's terms alone, and 1 + L folds into A's binding rate and Q's dissociation
    # constant: 6 fitted, rank 5.
    @test _testhelper_fitted(m) == _testhelper_identifiable_rank(m) == 5
    for k in (only_a, only_p, only_q, a_and_p, p_and_q)
        @test _testhelper_fitted(k) == _testhelper_identifiable_rank(k) == 6
    end
    for k in (a_and_q, all_three)
        @test _testhelper_fitted(k) == 6 && _testhelper_identifiable_rank(k) == 5
    end
end

@testset "_expand_to_allosteric on a random-order merged bi-bi led by plain bindings" begin
    # The random-order merged bi-bi: A and B bind free E in either order, and the second
    # binding is fused with the chemistry. A's steady-state group holds E + A → E(A) and
    # the fused E(B) + A → E(P, Q); B's rapid-equilibrium group holds E + B ⇌ E(B) and
    # the fused E(A) + B ⇌ E(P, Q). Canonical order sorts a group's steps by source form,
    # E before E(A) and E(B), so each group opens with its plain binding. A group holding
    # a chemistry step is chemistry whatever its first step, so both are :OnlyA in every
    # K-type child, and the subsets range over the Q and P groups: 2^2 − 1 = 3. Without
    # the :OnlyA A and B groups the check graph keeps the path E – E(Q) – E(P, Q), so
    # `_onlya_haldane_violation` refuses none. The reaction allows one catalytic subunit
    # and declares no regulator, so no V-type: three children. Read by its first step,
    # neither group would be chemistry, and the subsets would range over all four groups.
    m = EnzymeRates.Mechanism(@enzyme_mechanism begin
        substrates: A, B
        products: P, Q
        steps: begin
            (E + A <--> E(A), E(B) + A <--> E(P, Q))
            (E + B ⇌ E(B), E(A) + B ⇌ E(P, Q))
            E + Q ⇌ E(Q)
            E(Q) + P <--> E(P, Q)
        end
    end)
    rxn = @enzyme_reaction begin
        substrates: A[C], B[N]
        products: P[C], Q[N]
    end
    chemistry = [g for g in EnzymeRates.steps(m) if any(EnzymeRates._is_chemistry, g)]
    @test length(chemistry) == 2
    @test all(g -> EnzymeRates.is_binding(first(g)) &&
                   !EnzymeRates._is_chemistry(first(g)), chemistry)
    only_q = EnzymeRates.AllostericMechanism(@allosteric_mechanism begin
        substrates: A, B ; products: P, Q ; catalytic_multiplicity: 1
        catalytic_steps: begin
            (E + A <--> E(A), E(B) + A <--> E(P, Q))     :: OnlyA
            (E + B ⇌ E(B), E(A) + B ⇌ E(P, Q))           :: OnlyA
            E + Q ⇌ E(Q)                                 :: OnlyA
            E(Q) + P <--> E(P, Q)                        :: EqualAI
        end
    end)
    only_p = EnzymeRates.AllostericMechanism(@allosteric_mechanism begin
        substrates: A, B ; products: P, Q ; catalytic_multiplicity: 1
        catalytic_steps: begin
            (E + A <--> E(A), E(B) + A <--> E(P, Q))     :: OnlyA
            (E + B ⇌ E(B), E(A) + B ⇌ E(P, Q))           :: OnlyA
            E + Q ⇌ E(Q)                                 :: EqualAI
            E(Q) + P <--> E(P, Q)                        :: OnlyA
        end
    end)
    p_and_q = EnzymeRates.AllostericMechanism(@allosteric_mechanism begin
        substrates: A, B ; products: P, Q ; catalytic_multiplicity: 1
        catalytic_steps: begin
            (E + A <--> E(A), E(B) + A <--> E(P, Q))     :: OnlyA
            (E + B ⇌ E(B), E(A) + B ⇌ E(P, Q))           :: OnlyA
            E + Q ⇌ E(Q)                                 :: OnlyA
            E(Q) + P <--> E(P, Q)                        :: OnlyA
        end
    end)
    kids = EnzymeRates._expand_to_allosteric(m, rxn)
    @test length(kids) == 3
    @test Set(kids) == Set([only_q, only_p, p_and_q])
    # The parent has 5 fitted constants, all identifiable. Each child adds L. Its dead
    # inactive conformation reaches no E(B), since B's group is :OnlyA, so it adds L·E to
    # E's rapid-equilibrium segment {E, E(B), E(Q)}, times 1 + Q/K_Q where Q stays
    # :EqualAI. A's group takes A up at E and at E(B) under one rate constant, so the
    # denominator's A·B term is 2·B/K_B times its A term whatever L is: K_B is pinned and
    # cannot also absorb the 1 + L at E. L shows in every child, the two whose inactive
    # conformation binds nothing included: 6 fitted, rank 6.
    @test _testhelper_fitted(m) == _testhelper_identifiable_rank(m) == 5
    for k in (only_q, only_p, p_and_q)
        @test _testhelper_fitted(k) == _testhelper_identifiable_rank(k) == 6
    end
end

@testset "expand_mechanisms expands a parent whose binding changes conformation" begin
    rxn = @enzyme_reaction begin
        substrates: A[CX], B[N]
        products: P[C], Q[NX]
    end
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

@testset "expand_mechanisms rejects a parent with a zero-flux steady-state group" begin
    # E + A and E + Q at steady state while E(Q) + A ⇌ E(A, Q) ⇌ E(A) + Q stays at
    # rapid equilibrium: {E} joins the rest by two edges of weight 0, so neither
    # group carries flux and their constants enter the rate only as ratios. The
    # flip never emits this child; a hand-written parent is refused.
    rxn = @enzyme_reaction begin
        substrates: A[C], B[N]
        products: P[C], Q[N]
    end
    shunt = EnzymeRates.Mechanism(@enzyme_mechanism begin
        substrates: A, B
        products: P, Q
        steps: begin
            E + A <--> E(A)
            E(Q) + A ⇌ E(A, Q)
            E + Q <--> E(Q)
            E(A) + Q ⇌ E(A, Q)
            E(A) + B ⇌ E(A, B)
            E(A, B) <--> E(P, Q)
            E(Q) + P ⇌ E(P, Q)
        end
    end)
    err = try
        EnzymeRates.expand_mechanisms([shunt], rxn); nothing
    catch e
        e
    end
    @test err isa ErrorException
    @test occursin("carries net flux", err.msg)
    @test occursin("{E_A → EA}", err.msg)
end

@testset "expand_mechanisms rejects a parent whose copy group is redundant" begin
    # A bound as its own competitive inhibitor at E alone: E(A::Inh) has the
    # composition of E(A), formed by the A group alone, so the gauge exists and the
    # copy's constant enters the rate only added to K_A's. The dead-end move never
    # emits this child; a hand-written parent is refused.
    rxn = @enzyme_reaction begin
        substrates: A[C], B[N]
        products: P[C], Q[N]
        dead_end_inhibitors: A
    end
    twin_only = EnzymeRates.Mechanism(@enzyme_mechanism begin
        substrates: A, B
        products: P, Q
        steps: begin
            E + A ⇌ E(A)
            E(A) + B ⇌ E(A, B)
            E(A, B) <--> E(P, Q)
            E(Q) + P ⇌ E(P, Q)
            E + Q ⇌ E(Q)
            E + A::Inh ⇌ E(A::Inh)
        end
    end)
    err = try
        EnzymeRates.expand_mechanisms([twin_only], rxn); nothing
    catch e
        e
    end
    @test err isa ErrorException
    @test occursin("duplicates a productive form", err.msg)
    @test occursin("dwell gauge", err.msg)
    # The same copy placed where it also creates a new complex is a valid parent.
    kept = EnzymeRates.Mechanism(@enzyme_mechanism begin
        substrates: A, B
        products: P, Q
        steps: begin
            E + A ⇌ E(A)
            E(A) + B ⇌ E(A, B)
            E(A, B) <--> E(P, Q)
            E(Q) + P ⇌ E(P, Q)
            (E + Q ⇌ E(Q), E(A::Inh) + Q ⇌ E(A::Inh, Q))
            (E + A::Inh ⇌ E(A::Inh), E(Q) + A::Inh ⇌ E(A::Inh, Q))
        end
    end)
    @test !isempty(EnzymeRates.expand_mechanisms([kept], rxn))
    @test EnzymeRates._assert_emission_rules(kept) === nothing
    # A copy at E alone whose A group also binds at E(Q), forming E(A, Q): every
    # complex duplicates a productive form, but the gauge fails, so it is a valid
    # parent.
    no_gauge = EnzymeRates.Mechanism(@enzyme_mechanism begin
        substrates: A, B
        products: P, Q
        steps: begin
            (E + A ⇌ E(A), E(Q) + A ⇌ E(A, Q))
            E(A) + B ⇌ E(A, B)
            E(A, B) <--> E(P, Q)
            E(Q) + P ⇌ E(P, Q)
            (E + Q ⇌ E(Q), E(A) + Q ⇌ E(A, Q))
            E + A::Inh ⇌ E(A::Inh)
        end
    end)
    @test EnzymeRates._assert_emission_rules(no_gauge) === nothing
    @test !isempty(EnzymeRates.expand_mechanisms([no_gauge], rxn))
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
    @test all(c -> EnzymeRates.reaction(c) == rxn, kids)
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

        key_a = EnzymeRates._rate_eq_dedup_key(rate_equation_string(ldh_m_a))
        @test key_a == EnzymeRates._rate_eq_dedup_key(rate_equation_string(ldh_m_b))

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
        @test key_a != EnzymeRates._rate_eq_dedup_key(rate_equation_string(ldh_m_c))
    end

    @testset "rate_equation_string emits section labels" begin
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
    # (e.g. `K_ELactateNADH_to_ENADH_Lactate` vs the active-state
    # `K_ELactateNAD_to_ENAD_Lactate`), un-lumping a shared constant so the combined
    # solve fails to merge it and over-counts.
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
        rs = first(EnzymeRates.steps(m)[g])
        isempty(EnzymeRates.bound(EnzymeRates.from_species(rs))) &&
            EnzymeRates.bound_metabolite(rs) !== nothing
    end
    n_reproducers = 0
    for m in EnzymeRates.init_mechanisms(rxn)
        has_lumped_bound_group(m) || continue
        pn = _testhelper_fitted(m)
        free_e = Set(g for g in EnzymeRates.kinetic_groups(m)
                     if is_free_e_binding(m, g))
        iso = Set(_testhelper_iso_groups(m))
        children = EnzymeRates._expand_to_allosteric(m, rxn)
        # Every emitted child compiles.
        for c in children
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
        @test _testhelper_fitted(only(clean)) == pn + 1
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

    # The seeds of a reaction with two undesignated allosteric regulators, both required.
    # Their names differ from the `rxn` and `seeds` that other subtests assign.
    rxn_xy = @enzyme_reaction begin
        substrates: A[C6H12O6]
        products:   B[C6H12O6]
        allosteric_regulators: X, Y
        oligomeric_state: 2
    end
    seeds_xy = EnzymeRates.seed_mechanisms(rxn_xy, Set([:X, :Y]), Set{Symbol}())

    @testset "undesignated: skeletons × 2²" begin
        @test !isempty(seeds_xy)
        # Every seed is a fully-regulated, cheap-state AllostericMechanism with X
        # and Y each at their own single-ligand site.
        for s in seeds_xy
            @test s isa EnzymeRates.AllostericMechanism
            @test EnzymeRates._bound_allo_regs(s) == Set([:X, :Y])
            sites = EnzymeRates.regulatory_sites(s)
            @test length(sites) == 2
            @test all(length(EnzymeRates.ligands(si)) == 1 for si in sites)
            @test :NonequalAI ∉ [EnzymeRates.cat_allo_states(s);
                                 [st for si in sites for st in EnzymeRates.allo_states(si)]]
        end
        # Each lineage contributes exactly 2² = 4 state assignments.
        skels = Dict{UInt64, Int}()
        for s in seeds_xy
            k = skeleton_key(s)
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
        seeds_d = EnzymeRates.seed_mechanisms(rxn_d, Set([:X, :Y]), Set{Symbol}())
        @test !isempty(seeds_d)
        # One state assignment per lineage: seed count equals the skeleton count.
        skels_d = Set(skeleton_key(s) for s in seeds_d)
        @test length(seeds_d) == length(skels_d)
        # The same catalytic-allostery skeletons as the undesignated build.
        @test length(seeds_xy) == 4 * length(seeds_d)
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
        n_required = 2
        # A required allosteric regulator adds one dissociation constant; the two
        # required regulators additionally lift the mechanism to two
        # conformations, so a seed carries its base catalytic count + n_required
        # + 1 (L). Compile is slow; sample a few seeds.
        for s in seeds_xy[1:min(3, length(seeds_xy))]
            base = EnzymeRates.Mechanism(
                EnzymeRates.reaction(s), copy(EnzymeRates.steps(s)))
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
        base_count = _testhelper_fitted(inits[1])
        for s in seeds[1:min(2, length(seeds))]
            fp = EnzymeRates.fitted_params(s)
            @test length(fp) == base_count + 1
            @test !(:L in fp)
        end
    end
end

@testset "seed_mechanisms errors when no mechanism binds every required regulator" begin
    # Uni-uni with S as its own competitive inhibitor: the only site that carries
    # neither S nor P is free E, where the copy duplicates E(S). No seed exists,
    # and the beam must say so rather than return nothing.
    rxn = @enzyme_reaction begin
        substrates: S[C]
        products: P[C]
        dead_end_inhibitors: S
    end
    err = try
        EnzymeRates.seed_mechanisms(rxn, Set{Symbol}(), Set([:S])); nothing
    catch e
        e
    end
    @test err isa ErrorException
    @test occursin("required regulator", err.msg)
    @test occursin("competitive inhibitors: S", err.msg)
    @test occursin("allosteric regulators: none", err.msg)
    @test occursin("a uni-uni seed has no such site", err.msg)
    @test occursin("optional_competitive_inhibitors", err.msg)
    # Bi-bi with A as its own inhibitor has placements that create a new complex.
    bibi = @enzyme_reaction begin
        substrates: A[C], B[N]
        products: P[C], Q[N]
        dead_end_inhibitors: A
    end
    seeds = EnzymeRates.seed_mechanisms(bibi, Set{Symbol}(), Set([:A]))
    @test !isempty(seeds)
    @test all(m -> EnzymeRates._assert_emission_rules(m) === nothing, seeds)
    @test all(m -> :A in EnzymeRates._bound_comp_inhibitors(m), seeds)
end

@testset "seed_mechanisms grows the merged and Theorell–Chance seeds" begin
    ER = EnzymeRates
    # A required dead-end inhibitor is placed on every init mechanism, the merged and
    # Theorell–Chance variants included, so some seed carries a fused binding.
    rxn = @enzyme_reaction begin
        substrates: A[C], B[N]
        products: P[C], Q[N]
        dead_end_inhibitors: I
    end
    seeds = ER.seed_mechanisms(rxn, Set{Symbol}(), Set([:I]))
    fused(m) = any(s -> ER.is_binding(s) && ER._is_chemistry(s),
                   Iterators.flatten(ER.steps(m)))
    @test all(m -> :I in ER._bound_comp_inhibitors(m), seeds)
    @test any(fused, seeds)
end

@testset "seed_mechanisms wave-parallel equivalence" begin
    # Inline serial FIFO BFS = the reference the parallel version must match.
    binds_all_required(m, req_allo, req_comp) =
        issubset(req_allo, EnzymeRates._bound_allo_regs(m)) &&
        issubset(req_comp, EnzymeRates._bound_comp_inhibitors(m))
    function serial_seed_reference(rxn, req_allo, req_comp)
        visited = Set{UInt64}()
        queue = Union{EnzymeRates.Mechanism, EnzymeRates.AllostericMechanism}[]
        seeds = Union{EnzymeRates.Mechanism, EnzymeRates.AllostericMechanism}[]
        enq(m) = begin
            h = hash(m)
            h in visited && return
            push!(visited, h); push!(queue, m)
            binds_all_required(m, req_allo, req_comp) && push!(seeds, m)
        end
        for m in EnzymeRates.init_mechanisms(rxn); enq(m); end
        while !isempty(queue)
            m = popfirst!(queue)
            lifted = isempty(req_allo) ? [] : m isa EnzymeRates.Mechanism ?
                EnzymeRates._expand_to_allosteric(m, rxn) :
                EnzymeRates._expand_add_allosteric_regulator(m, rxn)
            for c in [lifted; EnzymeRates._expand_add_dead_end_regulator(m, rxn)]
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
    @test all(m -> binds_all_required(m, req, empty), got)
    @test EnzymeRates.seed_mechanisms(uni_uni_allo_reg, req, empty) == got  # deterministic
end

# Tag of each catalytic kinetic group, keyed by the metabolite it binds
# (`:chem` for the chemical step). Group index is never a stable key: the
# AllostericMechanism constructor canonicalizes group order.
_testhelper_tags_by_bound_metabolite(x) = Dict(
    (bm = EnzymeRates.bound_metabolite(grp[1]);
     bm === nothing ? :chem : EnzymeRates.name(bm)) =>
        EnzymeRates.cat_allo_states(x)[g]
    for (g, grp) in enumerate(EnzymeRates.steps(x)))

# ─── _expand_change_allo_state Haldane filter ──────────────────────────
@testset "_expand_change_allo_state drops the chemical step's relaxation" begin
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

    result = ER._expand_change_allo_state(parent)
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
    kids = ER._expand_change_allo_state(balanced)
    @test length(kids) == 2           # the chemical-step relaxation is dropped
    @test !any(_testhelper_tags_by_bound_metabolite(k)[:chem] == :NonequalAI for k in kids)
    @test Set(_testhelper_tags_by_bound_metabolite(k) for k in kids) == Set([
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
    reg_kids = ER._expand_change_allo_state(reg_parent)
    # 2 surviving cat relaxations + the R-ligand relaxation, which survives.
    @test length(reg_kids) == 3
    @test count(k -> ER.allo_states(only(ER.regulatory_sites(k))) ==
                     [:NonequalAI], reg_kids) == 1
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
    # every :OnlyA child is dead-inactive: all iso steps :OnlyA
    for k in onlya
        @test all(EnzymeRates.cat_allo_states(k)[g] === :OnlyA
                  for g in _testhelper_iso_groups(k))
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
    hasonlyabind(k) = any(EnzymeRates.cat_allo_states(k)[g] === :OnlyA &&
                          EnzymeRates.is_binding(EnzymeRates.steps(k)[g][1])
                          for g in eachindex(EnzymeRates.steps(k)))
    for k in kids
        @test !(hasonlyabind(k) &&
                !all(EnzymeRates.cat_allo_states(k)[g] === :OnlyA
                     for g in _testhelper_iso_groups(k)))
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
    kids = EnzymeRates._expand_change_allo_state(am)
    for k in kids
        tags = EnzymeRates.cat_allo_states(k)
        onlya_iso = any(tags[g] === :OnlyA for g in _testhelper_iso_groups(k))
        live_iso = any(tags[g] !== :OnlyA for g in _testhelper_iso_groups(k))
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
                  for g in _testhelper_iso_groups(k)) for k in kids)
end

@testset "change_allo_state relaxes the bindings of a hand-built partial parent" begin
    # A partial parent: every binding :OnlyA beside a live chemistry step. The
    # enumeration never builds one, since every mechanism it emits is fully dead or fully
    # live in the inactive conformation; this one is written by hand. Each binding
    # relaxation leaves three :OnlyA bindings of both signs on the catalytic cycle, so
    # the Haldane check accepts all four. The chemistry relaxation is not tried while a
    # binding is :OnlyA.
    partial = EnzymeRates.AllostericMechanism(@allosteric_mechanism begin
        substrates: A, B ; products: P, Q ; catalytic_multiplicity: 1
        catalytic_steps: begin
            E + A ⇌ E(A)                :: OnlyA
            E(A) + B ⇌ E(A, B)          :: OnlyA
            E(A, B) <--> E(P, Q)        :: EqualAI
            E(Q) + P ⇌ E(P, Q)          :: OnlyA
            E + Q ⇌ E(Q)                :: OnlyA
        end
    end)
    a_relaxed = EnzymeRates.AllostericMechanism(@allosteric_mechanism begin
        substrates: A, B ; products: P, Q ; catalytic_multiplicity: 1
        catalytic_steps: begin
            E + A ⇌ E(A)                :: NonequalAI
            E(A) + B ⇌ E(A, B)          :: OnlyA
            E(A, B) <--> E(P, Q)        :: EqualAI
            E(Q) + P ⇌ E(P, Q)          :: OnlyA
            E + Q ⇌ E(Q)                :: OnlyA
        end
    end)
    b_relaxed = EnzymeRates.AllostericMechanism(@allosteric_mechanism begin
        substrates: A, B ; products: P, Q ; catalytic_multiplicity: 1
        catalytic_steps: begin
            E + A ⇌ E(A)                :: OnlyA
            E(A) + B ⇌ E(A, B)          :: NonequalAI
            E(A, B) <--> E(P, Q)        :: EqualAI
            E(Q) + P ⇌ E(P, Q)          :: OnlyA
            E + Q ⇌ E(Q)                :: OnlyA
        end
    end)
    p_relaxed = EnzymeRates.AllostericMechanism(@allosteric_mechanism begin
        substrates: A, B ; products: P, Q ; catalytic_multiplicity: 1
        catalytic_steps: begin
            E + A ⇌ E(A)                :: OnlyA
            E(A) + B ⇌ E(A, B)          :: OnlyA
            E(A, B) <--> E(P, Q)        :: EqualAI
            E(Q) + P ⇌ E(P, Q)          :: NonequalAI
            E + Q ⇌ E(Q)                :: OnlyA
        end
    end)
    q_relaxed = EnzymeRates.AllostericMechanism(@allosteric_mechanism begin
        substrates: A, B ; products: P, Q ; catalytic_multiplicity: 1
        catalytic_steps: begin
            E + A ⇌ E(A)                :: OnlyA
            E(A) + B ⇌ E(A, B)          :: OnlyA
            E(A, B) <--> E(P, Q)        :: EqualAI
            E(Q) + P ⇌ E(P, Q)          :: OnlyA
            E + Q ⇌ E(Q)                :: NonequalAI
        end
    end)
    kids = EnzymeRates._expand_change_allo_state(partial)
    @test length(kids) == 4
    @test Set(kids) == Set([a_relaxed, b_relaxed, p_relaxed, q_relaxed])
    # expand_mechanisms passes them on.
    @test issubset(kids,
                   EnzymeRates.expand_mechanisms([partial], EnzymeRates.reaction(partial)))
end
