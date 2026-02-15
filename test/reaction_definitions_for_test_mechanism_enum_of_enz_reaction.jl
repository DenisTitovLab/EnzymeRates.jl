# Test specifications for mechanism enumeration pipeline
# Each reaction is defined with expected counts at each stage

# ── EnumerationTestSpec struct ────────────────────────────────────────────

"""
Data-driven test specification for mechanism enumeration.
Contains expected counts at each pipeline stage.
"""
Base.@kwdef struct EnumerationTestSpec
    name::String
    reaction::Any           # EnzymeReaction instance
    max_forms::Int          # max_forms passed to enumerate_mechanisms

    # Stage counts (verified by tests)
    expected_n_forms::Int                # enumerate_enzyme_forms
    expected_n_catalytic::Int            # _enumerate_only_catalytic_mechanisms
    expected_n_cat_with_act::Int         # after _generate_activator_configs
    expected_n_cat_act_de::Int           # after _enumerate_dead_end_configs

    # RE/SS stage
    skip_ress_test::Bool = false         # skip for slow reactions
    expected_n_total::Int = 0            # enumerate_mechanisms total

    # Performance
    max_enumeration_time::Float64 = Inf  # max seconds; Inf = skip check
end

# ── Independent verification helpers ──────────────────────────────────────

"""
Independently compute the expected dead-end count by enumerating
CouldExist edges from topology forms to unused forms, respecting
dependency chains and max_forms budget.
"""
function _compute_expected_dead_end_count(
    base_spec::EnzymeRates.MechanismSpec,
    forms::Vector{EnzymeRates.EnzymeFormSpec},
    max_forms::Int,
)
    topo = EnzymeRates._spec_to_edges(base_spec, forms)
    topo_forms = Set(Iterators.flatten(topo))
    budget = max_forms - length(topo_forms)
    budget < 0 && return 0

    # Find level-1 dead-end candidates (reachable from topo forms)
    level1 = Set{Int}()
    for fi in topo_forms
        for (j, fj) in enumerate(forms)
            j in topo_forms && continue
            ec, _, _ = EnzymeRates.edge_class(forms[fi], fj)
            ec isa EnzymeRates.CouldExist && push!(level1, j)
        end
    end

    # Find level-2 dead-end candidates (reachable from level-1 only)
    level2 = Set{Int}()
    for ci in level1
        for (j, fj) in enumerate(forms)
            j in topo_forms && continue
            j in level1 && continue
            ec, _, _ = EnzymeRates.edge_class(forms[ci], fj)
            ec isa EnzymeRates.CouldExist && push!(level2, j)
        end
    end

    candidates = sort(collect(union(level1, level2)))

    # Enumerate valid subsets: each form must be reachable from topo
    # or from an already-included dead-end
    count = Ref(1)  # count empty subset (= base itself)
    _count_valid_subsets!(
        count, candidates, 1, Int[], budget, topo_forms, level1, forms,
    )
    count[]
end

function _count_valid_subsets!(
    count, candidates, idx, current, budget, topo_forms, level1, forms,
)
    budget <= 0 && return
    for i in idx:length(candidates)
        c = candidates[i]
        # Check reachability: from topo or from already-included dead-end
        reachable = false
        for fi in topo_forms
            ec, _, _ = EnzymeRates.edge_class(forms[fi], forms[c])
            ec isa EnzymeRates.CouldExist && (reachable = true; break)
        end
        if !reachable
            for fi in current
                ec, _, _ = EnzymeRates.edge_class(forms[fi], forms[c])
                ec isa EnzymeRates.CouldExist &&
                    (reachable = true; break)
            end
        end
        !reachable && continue

        push!(current, c)
        count[] += 1
        _count_valid_subsets!(
            count, candidates, i + 1,
            current, budget - 1, topo_forms, level1, forms,
        )
        pop!(current)
    end
end

"""
Independently compute the expected RE/SS + constraint count by iterating
all RE/SS masks and equivalent group constraint combinations.
"""
function _compute_independent_ress_count(
    spec::EnzymeRates.MechanismSpec,
    forms::Vector{EnzymeRates.EnzymeFormSpec},
)
    edges = EnzymeRates._spec_to_edges(spec, forms)
    n = length(edges)

    # Find equivalent binding step groups
    equiv_groups = EnzymeRates._find_equivalent_groups(edges, forms)

    total = 0
    # Iterate all non-all-RE masks (at least one SS)
    for re_mask in 0:((1 << n) - 2)
        eq_steps = Bool[
            (re_mask >> (i - 1)) & 1 == 1 for i in 1:n
        ]
        # Count valid groups (all steps same type)
        n_valid = 0
        for group in equiv_groups
            first_re = eq_steps[group[1]]
            all(eq_steps[s] == first_re for s in group) &&
                (n_valid += 1)
        end
        # Each valid group: constrained or unconstrained = 2 options
        total += 1 << n_valid
    end
    total
end

# ── Build specifications ──────────────────────────────────────────────────

function build_enumeration_test_specs()
    specs = EnumerationTestSpec[]

    # 1. Uni-Uni: simplest case
    let
        rxn = @enzyme_reaction begin
            substrates: S[C]
            products:   P[C]
        end
        push!(specs, EnumerationTestSpec(
            name = "Uni-Uni",
            reaction = rxn,
            max_forms = 6,
            expected_n_forms = 3,
            expected_n_catalytic = 1,
            expected_n_cat_with_act = 1,
            expected_n_cat_act_de = 1,
            expected_n_total = 7,
            max_enumeration_time = 5.0,
        ))
    end

    # 2. Uni-Uni + 1 dead-end inhibitor
    let
        rxn = @enzyme_reaction begin
            substrates: S[C]
            products:   P[C]
            regulators: I[N]
        end
        push!(specs, EnumerationTestSpec(
            name = "Uni-Uni 1 Regulator",
            reaction = rxn,
            max_forms = 6,
            expected_n_forms = 6,
            expected_n_catalytic = 1,
            expected_n_cat_with_act = 3,
            expected_n_cat_act_de = 10,
            expected_n_total = 1415,
            max_enumeration_time = 5.0,
        ))
    end

    # 3. Uni-Uni + 2 regulators (chain dead-ends)
    let
        rxn = @enzyme_reaction begin
            substrates: S[C]
            products:   P[C]
            regulators: I[N], J[P2]
        end
        push!(specs, EnumerationTestSpec(
            name = "Uni-Uni 2 Regulators",
            reaction = rxn,
            max_forms = 6,
            expected_n_forms = 12,
            expected_n_catalytic = 1,
            expected_n_cat_with_act = 5,
            expected_n_cat_act_de = 79,
            expected_n_total = 7208,
            max_enumeration_time = 10.0,
        ))
    end

    # 4. Bi-Bi (multi-cycle topologies)
    let
        rxn = @enzyme_reaction begin
            substrates: A[C], B[N]
            products:   P[C], Q[N]
        end
        push!(specs, EnumerationTestSpec(
            name = "Bi-Bi",
            reaction = rxn,
            max_forms = 11,
            expected_n_forms = 11,
            expected_n_catalytic = 9,
            expected_n_cat_with_act = 9,
            expected_n_cat_act_de = 9,
            expected_n_total = 2094,
            max_enumeration_time = 5.0,
        ))
    end

    # 5. Bi-Bi Ping Pong (residual forms)
    let
        rxn = @enzyme_reaction begin
            substrates: A[CX], B[N]
            products:   P[C], Q[NX]
        end
        push!(specs, EnumerationTestSpec(
            name = "Bi-Bi Ping Pong",
            reaction = rxn,
            max_forms = 20,
            expected_n_forms = 17,
            expected_n_catalytic = 9,
            expected_n_cat_with_act = 9,
            expected_n_cat_act_de = 9,
            expected_n_total = 2094,
            max_enumeration_time = 5.0,
        ))
    end

    # 6. Bi-Bi + 1 regulator
    let
        rxn = @enzyme_reaction begin
            substrates: A[C], B[N]
            products:   P[C], Q[N]
            regulators: I[P2]
        end
        push!(specs, EnumerationTestSpec(
            name = "Bi-Bi 1 Regulator",
            reaction = rxn,
            max_forms = 6,
            expected_n_forms = 22,
            expected_n_catalytic = 8,
            expected_n_cat_with_act = 8,
            expected_n_cat_act_de = 28,
            expected_n_total = 2206,
            max_enumeration_time = 5.0,
        ))
    end

    # 7. Bi-Bi Ping Pong + 1 regulator
    let
        rxn = @enzyme_reaction begin
            substrates: A[CX], B[N]
            products:   P[C], Q[NX]
            regulators: I[P2]
        end
        push!(specs, EnumerationTestSpec(
            name = "Bi-Bi PP 1 Regulator",
            reaction = rxn,
            max_forms = 6,
            expected_n_forms = 34,
            expected_n_catalytic = 8,
            expected_n_cat_with_act = 8,
            expected_n_cat_act_de = 28,
            expected_n_total = 2206,
            max_enumeration_time = 5.0,
        ))
    end

    # 8. Bi-Bi Ping Pong + 2 regulators
    let
        rxn = @enzyme_reaction begin
            substrates: A[CX], B[N]
            products:   P[C], Q[NX]
            regulators: I[P2], J[Y]
        end
        push!(specs, EnumerationTestSpec(
            name = "Bi-Bi PP 2 Regulators",
            reaction = rxn,
            max_forms = 6,
            expected_n_forms = 68,
            expected_n_catalytic = 8,
            expected_n_cat_with_act = 8,
            expected_n_cat_act_de = 48,
            expected_n_total = 3466,
            max_enumeration_time = 10.0,
        ))
    end

    return specs
end

const ENUMERATION_TEST_SPECS = build_enumeration_test_specs()
