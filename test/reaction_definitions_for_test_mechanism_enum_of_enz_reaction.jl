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
    expected_n_catalytic::Int            # catalytic topologies
    expected_n_cat_with_act::Int         # after activator configs
    expected_n_cat_act_de::Int           # after dead-end configs

    # RE/SS stage
    skip_ress_test::Bool = false         # skip for slow reactions
    expected_n_total::Int = 0            # enumerate_mechanisms total

    # Rate equation smoke test
    test_rate_equation::Bool = true

    # Performance
    max_enumeration_time::Float64 = Inf  # max seconds; Inf = skip check
end

# ── Internal helpers (no EnzymeRates._* calls) ───────────────────────────

"""Extract form-index edge pairs from a MechanismSpec's reactions.

Inlines `_spec_to_edges` logic using only public struct fields.
"""
function _spec_to_form_edges(spec, forms)
    name_to_idx = Dict(f.name => i for (i, f) in enumerate(forms))
    edges = Tuple{Int,Int}[]
    for (lhs, rhs) in spec.reactions
        from_idx = nothing
        to_idx = nothing
        for sym in lhs
            idx = get(name_to_idx, sym, nothing)
            idx !== nothing && (from_idx = idx)
        end
        for sym in rhs
            idx = get(name_to_idx, sym, nothing)
            idx !== nothing && (to_idx = idx)
        end
        if from_idx !== nothing && to_idx !== nothing
            push!(edges, (from_idx, to_idx))
        end
    end
    edges
end

# ── Independent verification helpers ─────────────────────────────────────

"""
Independently compute the expected dead-end count using the thermodynamic
box rule. Per topology form, choose which regulator sites bind (2^n_reg
subsets). Multi-reg subsets force box closure. Cartesian product across
topology forms with max_forms budget.

Uses only public struct fields — no `EnzymeRates._*` calls.
"""
function _compute_expected_dead_end_count(
    base_spec::EnzymeRates.MechanismSpec,
    forms::Vector{EnzymeRates.EnzymeFormSpec},
    max_forms::Int,
)
    topo = _spec_to_form_edges(base_spec, forms)
    topo_forms_set = Set(Iterators.flatten(topo))
    topo_forms = sort(collect(topo_forms_set))
    budget = max_forms - length(topo_forms)
    budget < 0 && return 0

    per_form_option_sizes = Vector{Int}[]
    for fi in topo_forms
        # Reg site positions (inlined from _reg_site_positions)
        reg_sites = [k for k in eachindex(forms[fi].sites)
                     if forms[fi].sites[k].role == :reg &&
                        forms[fi].sites[k].atoms === nothing]
        n_reg = length(reg_sites)
        seen = Set{Vector{Int}}([Int[]])
        sizes = [0]

        for mask in 1:((1 << n_reg) - 1)
            chosen = [reg_sites[k] for k in 1:n_reg
                      if (mask >> (k - 1)) & 1 == 1]
            de_forms = Int[]
            valid = true
            for sub_mask in 1:((1 << length(chosen)) - 1)
                positions = [chosen[k]
                             for k in 1:length(chosen)
                             if (sub_mask >> (k - 1)) & 1 == 1]
                # Find dead-end form (inlined from _find_dead_end_form)
                base = forms[fi]
                form_idx = findfirst(forms) do fj
                    all(1:length(base.sites)) do k
                        if k in positions
                            fj.sites[k].atoms !== nothing
                        else
                            base.sites[k].atoms == fj.sites[k].atoms
                        end
                    end
                end
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

Uses only public struct fields — no `EnzymeRates._*` calls.
"""
function _compute_independent_ress_count(
    spec::EnzymeRates.MechanismSpec,
    forms::Vector{EnzymeRates.EnzymeFormSpec},
)
    edges = _spec_to_form_edges(spec, forms)
    n = length(edges)

    # Find equivalent groups (inlined from _find_equivalent_groups)
    binding_key = Dict{Tuple{Symbol,Int}, Vector{Int}}()
    for (i, (a, b)) in enumerate(edges)
        # Detect binding/release: exactly 1 site changes occupancy
        diff_count = 0
        diff_k = 0
        for k in 1:length(forms[a].sites)
            if forms[a].sites[k].atoms != forms[b].sites[k].atoms
                diff_count += 1
                diff_count == 1 && (diff_k = k)
            end
        end
        diff_count == 1 || continue  # skip isomerization

        k = diff_k
        a_occ = forms[a].sites[k].atoms !== nothing
        site = a_occ ? forms[a].sites[k] : forms[b].sites[k]
        site.role == :prod && continue
        key = (site.metabolite, site.index)
        push!(get!(binding_key, key, Int[]), i)
    end
    equiv_groups = [sort(indices) for (_, indices) in binding_key
                    if length(indices) >= 2]
    sort!(equiv_groups; by=first)

    total = 0
    for re_mask in 0:((1 << n) - 2)
        eq_steps = Bool[
            (re_mask >> (i - 1)) & 1 == 1 for i in 1:n
        ]
        n_valid = 0
        for group in equiv_groups
            first_re = eq_steps[group[1]]
            all(eq_steps[s] == first_re for s in group) &&
                (n_valid += 1)
        end
        total += 1 << n_valid
    end
    total
end

# ── Build specifications ─────────────────────────────────────────────────

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
            max_forms = 15,
            expected_n_forms = 22,
            expected_n_catalytic = 9,
            expected_n_cat_with_act = 27,
            expected_n_cat_act_de = 530,
            skip_ress_test = true,
            expected_n_total = 140203886,
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
            max_forms = 15,
            expected_n_forms = 34,
            expected_n_catalytic = 9,
            expected_n_cat_with_act = 27,
            expected_n_cat_act_de = 530,
            skip_ress_test = true,
            expected_n_total = 140203886,
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
            max_forms = 12,
            expected_n_forms = 68,
            expected_n_catalytic = 9,
            expected_n_cat_with_act = 41,
            expected_n_cat_act_de = 12106,
            skip_ress_test = true,
            expected_n_total = 1303046914,
        ))
    end

    # 9. Bi-Bi budget filtering (max_forms=5 restricts to single cycles)
    let
        rxn = @enzyme_reaction begin
            substrates: A[C], B[N]
            products:   P[C], Q[N]
        end
        push!(specs, EnumerationTestSpec(
            name = "Bi-Bi Budget",
            reaction = rxn,
            max_forms = 5,
            expected_n_forms = 11,
            expected_n_catalytic = 4,
            expected_n_cat_with_act = 4,
            expected_n_cat_act_de = 4,
            expected_n_total = 124,
            max_enumeration_time = 5.0,
        ))
    end

    return specs
end

const ENUMERATION_TEST_SPECS = build_enumeration_test_specs()
