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
Independently compute the expected dead-end count using the thermodynamic
box rule: per topology form, choose which regulator sites bind (2^n_reg
subsets). Multi-reg subsets force box closure. Cartesian product across
topology forms with max_forms budget.
"""
function _compute_expected_dead_end_count(
    base_spec::EnzymeRates.MechanismSpec,
    forms::Vector{EnzymeRates.EnzymeFormSpec},
    max_forms::Int,
)
    topo = EnzymeRates._spec_to_edges(base_spec, forms)
    topo_forms_set = Set(Iterators.flatten(topo))
    topo_forms = sort(collect(topo_forms_set))
    budget = max_forms - length(topo_forms)
    budget < 0 && return 0

    # Per topology form: compute option sizes (number of dead-end forms
    # added by each regulator-subset choice)
    per_form_option_sizes = Vector{Int}[]
    for fi in topo_forms
        reg_sites = EnzymeRates._reg_site_positions(forms[fi])
        n_reg = length(reg_sites)
        seen = Set{Vector{Int}}([Int[]])
        sizes = [0]  # empty option (no dead-end binding) always valid

        for mask in 1:((1 << n_reg) - 1)
            chosen = [reg_sites[k] for k in 1:n_reg
                      if (mask >> (k - 1)) & 1 == 1]
            de_forms = Int[]
            valid = true
            for sub_mask in 1:((1 << length(chosen)) - 1)
                positions = [chosen[k]
                             for k in 1:length(chosen)
                             if (sub_mask >> (k - 1)) & 1 == 1]
                form_idx = EnzymeRates._find_dead_end_form(
                    fi, positions, forms,
                )
                if form_idx === nothing
                    valid = false
                    break
                end
                form_idx in topo_forms_set && continue
                push!(de_forms, form_idx)
            end
            !valid && continue
            sort!(de_forms)
            if de_forms ∉ seen
                push!(seen, de_forms)
                push!(sizes, length(de_forms))
            end
        end
        push!(per_form_option_sizes, sizes)
    end

    # Count Cartesian product combinations respecting budget
    total = Ref(0)
    _count_cartesian!(total, per_form_option_sizes, 1, budget)
    total[]
end

function _count_cartesian!(total, option_sizes, idx, budget)
    if idx > length(option_sizes)
        total[] += 1
        return
    end
    for sz in option_sizes[idx]
        sz > budget && continue
        _count_cartesian!(total, option_sizes, idx + 1, budget - sz)
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
            expected_n_total = 2240,
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
            expected_n_cat_act_de = 34,
            expected_n_total = 6647,
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
