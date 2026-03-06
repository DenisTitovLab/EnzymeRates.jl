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

    # Stage counts (verified by tests)
    expected_n_forms::Int                # enumerate_enzyme_forms
    expected_n_catalytic::Int            # catalytic topologies
    expected_n_cat_de::Int               # after dead-end configs

    # RE/SS stage
    skip_ress_test::Bool = false         # skip for slow reactions
    expected_n_total::Int = 0            # enumerate_mechanisms total

    # Performance
    max_enumeration_time::Float64 = Inf  # max seconds; Inf = skip check
end

# ── Internal helpers (no EnzymeRates._* calls) ───────────────────────────

"""
Compute G (number of RE groups) for given edges and eq_steps via
union-find. Test-local reimplementation for independent verification.
"""
function _compute_re_group_count_test(edges, eq_steps)
    form_indices = collect(Set(Iterators.flatten(edges)))
    parent = Dict(i => i for i in form_indices)
    function find(x)
        while parent[x] != x
            parent[x] = parent[parent[x]]
            x = parent[x]
        end
        x
    end
    for (idx, (a, b)) in enumerate(edges)
        eq_steps[idx] || continue
        ra, rb = find(a), find(b)
        ra != rb && (parent[ra] = rb)
    end
    length(Set(find(i) for i in form_indices))
end

"""
Find the first isomerization edge (multi-site diff) in edge order.
Test-local reimplementation — uses site diffs, not adjacency dict.
"""
function _find_first_isomerization_test(edges, forms)
    for (i, (a, b)) in enumerate(edges)
        ndiff = count(
            k -> forms[a].sites[k].atoms != forms[b].sites[k].atoms,
            eachindex(forms[a].sites))
        ndiff > 1 && return i
    end
    return 1
end

"""
Compute the RE/SS + constraint count using brute-force enumeration
with G ≤ max_re_groups cap.

Baseline: first isomerization edge is always SS. Iterates over
subsets of remaining edges to make SS, keeps masks with 2 ≤ G ≤ 7,
and counts constraint combos for each valid mask.
"""
function _compute_expected_n_total(
    spec::EnzymeRates.MechanismSpec,
    forms::Vector{EnzymeRates.EnzymeFormSpec};
    max_re_groups::Int=7,
)
    edges = spec.edges
    n = length(edges)
    n == 0 && return 0
    iso_idx = _find_first_isomerization_test(edges, forms)
    equiv_groups = _find_equiv_groups(edges, forms)
    other_indices = [i for i in 1:n if i != iso_idx]
    n_other = length(other_indices)
    total = 0
    for ss_mask in 0:(1 << n_other) - 1
        eq_steps = fill(true, n)
        eq_steps[iso_idx] = false
        for (bit, idx) in enumerate(other_indices)
            (ss_mask >> (bit - 1)) & 1 == 1 &&
                (eq_steps[idx] = false)
        end
        G = _compute_re_group_count_test(edges, eq_steps)
        (G < 2 || G > max_re_groups) && continue
        valid_groups = [g for g in equiv_groups
            if all(eq_steps[s] == eq_steps[g[1]] for s in g)]
        total += 1 << length(valid_groups)
    end
    total
end

"""Find equivalent groups from edges: non-product binding edges grouped
by metabolite. Uses only public struct fields."""
function _find_equiv_groups(edges, forms)
    binding_key = Dict{Symbol,Vector{Int}}()
    for (i, (a, b)) in enumerate(edges)
        diff_count = 0
        diff_k = 0
        for k in 1:length(forms[a].sites)
            if forms[a].sites[k].atoms != forms[b].sites[k].atoms
                diff_count += 1
                diff_count == 1 && (diff_k = k)
            end
        end
        diff_count == 1 || continue

        k = diff_k
        a_occ = forms[a].sites[k].atoms !== nothing
        site = a_occ ? forms[a].sites[k] : forms[b].sites[k]
        site.role == :prod && continue
        push!(get!(binding_key, site.metabolite, Int[]), i)
    end
    equiv_groups = [sort(indices) for (_, indices) in binding_key
                    if length(indices) >= 2]
    sort!(equiv_groups; by=first)
end

"""
Compute expected dead-end mechanism count with regulator partitioning.

Sums over all 2^n_reg partitions of regulators into {dead-end, allosteric}.
For each partition, dead-end expansion uses formula: (2^r_de)^n_topo
per catalytic topology, where r_de = number of dead-end regulators.
"""
function _compute_expected_dead_end_count(
    catalytic_specs,
    forms::Vector{EnzymeRates.EnzymeFormSpec},
)
    reg_positions = [k for k in eachindex(forms[1].sites)
                     if forms[1].sites[k].role == :reg]
    n_reg = length(reg_positions)

    total = 0
    for reg_mask in 0:(1 << n_reg) - 1
        n_de = count(i -> (reg_mask >> (i - 1)) & 1 == 0,
            1:n_reg)
        for spec in catalytic_specs
            topo_set = Set(Iterators.flatten(spec.edges))
            n_topo = length(topo_set)
            total += (2^n_de)^n_topo
        end
    end
    total
end

"""
Compute expected total with oligomeric expansion (catalytic_n=N).

For each dead-end spec with k allosteric regulators:
  EM contribution: ress_count (same as without oligomeric)
  OEM contribution: k==0 ? ress_count : ress_count * N^k
Total = sum of (EM + OEM) over all dead-end specs.
"""
function _compute_expected_oligomeric_total(
    de_specs,
    forms::Vector{EnzymeRates.EnzymeFormSpec},
    catalytic_n::Int;
    max_re_groups::Int=7,
)
    total = 0
    for spec in de_specs
        ress = _compute_expected_n_total(spec, forms;
            max_re_groups)
        k = length(spec.allosteric_regulators)
        oem = k == 0 ? ress : ress * catalytic_n^k
        total += ress + oem  # EM + OEM
    end
    total
end

# ── Build specifications ─────────────────────────────────────────────────

function build_enumeration_test_specs()
    specs = EnumerationTestSpec[]

    # 1. Uni-Uni: simplest case
    let
        rxn = @enzyme_reaction begin
            substrates:S[C]
            products:P[C]
        end
        push!(specs, EnumerationTestSpec(
            name="Uni-Uni",
            reaction=rxn,

            expected_n_forms=3,
            expected_n_catalytic=1,
            # No regulators → dead-end = catalytic
            expected_n_cat_de=1,
            # 3 valid RE/SS masks with G cap
            expected_n_total=3,
            max_enumeration_time=5.0,
        ))
    end

    # 2. Uni-Uni + 1 regulator
    let
        rxn = @enzyme_reaction begin
            substrates:S[C]
            products:P[C]
            regulators: R
        end
        push!(specs, EnumerationTestSpec(
            name="Uni-Uni 1 Regulator",
            reaction=rxn,

            expected_n_forms=6,
            expected_n_catalytic=1,
            # 2 partitions: {R dead-end} + {R allosteric}
            # {R de}: (2^1)^3 = 8, {R al}: 1. Total = 9
            expected_n_cat_de=9,
            expected_n_total=341,
            max_enumeration_time=5.0,
        ))
    end

    # 3. Uni-Uni + 2 regulators (chain dead-ends)
    let
        rxn = @enzyme_reaction begin
            substrates:S[C]
            products:P[C]
            regulators: R1, R2
        end
        push!(specs, EnumerationTestSpec(
            name="Uni-Uni 2 Regulators",
            reaction=rxn,

            expected_n_forms=12,
            expected_n_catalytic=1,
            # 4 partitions: (2^2)^3 + 2*(2^1)^3 + 1 = 81
            expected_n_cat_de=81,
            expected_n_total=1246220,
            max_enumeration_time=10.0,
        ))
    end

    # 4. Uni-Bi + 1 regulator
    let
        rxn = @enzyme_reaction begin
            substrates:A[C2]
            products:P1[C], P2[C]
            regulators: R
        end
        push!(specs, EnumerationTestSpec(
            name="Uni-Bi",
            reaction=rxn,

            expected_n_forms=16,
            expected_n_catalytic=3,
            # 2 partitions: {R de} + {R al}
            # {R de}: 2*16 + 32 = 64, {R al}: 3. Total: 67
            expected_n_cat_de=67,
            expected_n_total=92177,
            max_enumeration_time=5.0,
        ))
    end

    # 5. Bi-Bi + 1 regulator (same atoms on both substrates)
    let
        rxn = @enzyme_reaction begin
            substrates:A[C], B[C]
            products:P[C], Q[C]
            regulators: R
        end
        push!(specs, EnumerationTestSpec(
            name="Bi-Bi + 1 Regulator",
            reaction=rxn,

            expected_n_forms=22,
            expected_n_catalytic=9,
            expected_n_cat_de=521,
            skip_ress_test=true,
            expected_n_total=28005686,
        ))
    end

    # 6. Bi-Bi + 1 regulator (different atoms)
    let
        rxn = @enzyme_reaction begin
            substrates:A[C], B[N]
            products:P[C], Q[N]
            regulators: I
        end
        push!(specs, EnumerationTestSpec(
            name="Bi-Bi 1 Regulator",
            reaction=rxn,

            expected_n_forms=22,
            expected_n_catalytic=9,
            expected_n_cat_de=521,
            skip_ress_test=true,
            expected_n_total=28005686,
        ))
    end

    # 7. Bi-Bi Ping Pong (no regulator)
    let
        rxn = @enzyme_reaction begin
            substrates:A[CX], B[N]
            products:P[C], Q[NX]
        end
        push!(specs, EnumerationTestSpec(
            name="Bi-Bi PP",
            reaction=rxn,

            expected_n_forms=17,
            expected_n_catalytic=10,
            # No regulators → dead-end = catalytic
            expected_n_cat_de=10,
            expected_n_total=989,
            max_enumeration_time=10.0,
        ))
    end

    return specs
end

const ENUMERATION_TEST_SPECS = build_enumeration_test_specs()
