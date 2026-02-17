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
Compute the RE/SS + constraint count using a closed-form formula.

For `n` edges and `k` non-overlapping equiv groups of sizes `g₁,...,gₖ`:

    f(n; g₁,...,gₖ) = 2^(n - Σgᵢ) × ∏(2^gᵢ + 2) - 2^k

Derivation:
- `2^n - 1` valid RE/SS masks (all-RE excluded)
- Each equiv group of size `g` is "valid" (all edges same RE/SS) in
  `2` of `2^g` bit patterns (all-0 or all-1)
- Each valid group independently contributes ×2 (constrained or not)
- The product factorizes over free edges and groups
- Subtract `2^k` for the all-RE mask where all groups are valid
"""
function _compute_expected_n_total(
    spec::EnzymeRates.MechanismSpec,
    forms::Vector{EnzymeRates.EnzymeFormSpec},
)
    edges = _spec_to_form_edges(spec, forms)
    n = length(edges)
    equiv_groups = _find_equiv_groups(edges, forms)
    k = length(equiv_groups)
    sum_g = sum(length, equiv_groups; init=0)
    result = 1 << (n - sum_g)
    for group in equiv_groups
        result *= (1 << length(group)) + 2
    end
    result - (1 << k)
end

"""Find equivalent groups from edges: non-product binding edges grouped
by (metabolite, site_index). Uses only public struct fields."""
function _find_equiv_groups(edges, forms)
    binding_key = Dict{Tuple{Symbol,Int},Vector{Int}}()
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
        key = (site.metabolite, site.index)
        push!(get!(binding_key, key, Int[]), i)
    end
    equiv_groups = [sort(indices) for (_, indices) in binding_key
                    if length(indices) >= 2]
    sort!(equiv_groups; by=first)
end

"""
Compute expected dead-end mechanism count assuming independent inhibitor
binding: each inhibitor independently binds or not at each topology form.

A regulator is classified as an inhibitor if its binding site is never
occupied in any topology form; otherwise it is an activator.

Formula per activator config: (2^r_inh)^n_topo
  - r_inh: number of inhibitor regulators (sites never occupied in topo)
  - n_topo: number of forms in the catalytic topology

Uses only public struct fields — no `EnzymeRates._*` calls.
"""
function _compute_expected_dead_end_count(
    activator_specs,
    forms::Vector{EnzymeRates.EnzymeFormSpec},
)
    reg_positions = [k for k in eachindex(forms[1].sites)
                     if forms[1].sites[k].role == :reg]

    total = 0
    for spec in activator_specs
        edges = _spec_to_form_edges(spec, forms)
        topo_set = Set(Iterators.flatten(edges))
        n_topo = length(topo_set)

        r_inh = count(reg_positions) do k
            !any(fi -> forms[fi].sites[k].atoms !== nothing, topo_set)
        end

        total += (2^r_inh)^n_topo
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
            max_forms=100,
            # E, ES, EP
            expected_n_forms=3,
            # E+S <-> ES <-> EP <-> E+P
            expected_n_catalytic=1,
            expected_n_cat_with_act=1,
            expected_n_cat_act_de=1,
            # SE/RE combose can be infered from the number of reaction in each mechanism
            # using the following equation 2^(n_reactions) - 1
            # "-1" is because all RE is not biochemically possible
            # 2^(3)-1 = 7 for RE/SE forms
            expected_n_total=7,
            max_enumeration_time=5.0,
        ))
    end

    # 2. Uni-Uni + 1 regulator
    let
        rxn = @enzyme_reaction begin
            substrates:S[C]
            products:P[C]
            regulators:R[N]
        end
        push!(specs, EnumerationTestSpec(
            name="Uni-Uni 1 Regulator",
            reaction=rxn,
            max_forms=100,
            # n_cat_forms × 2^regulators = 3 × 2^1 = 6
            expected_n_forms=6,
            # E+S <-> ES <-> EP <-> E+P
            expected_n_catalytic=1,
            # 1 catalytic +
            # 1 essential activator: E+A <-> EA+S <-> EAS <-> EAP <-> EA+P <-> E+A
            # 1 non-essential activator:
            #          E+S <-> ES  <-> EP  <-> E+P
            #          ↕       ↕        ↕      ↕
            # E+A <-> EA+S <-> EAS <-> EAP <-> EA+P <-> E+A
            expected_n_cat_with_act=3,
            # Regulator is either inhibitor or activator so dead-end only applies to catalytic form:
            # - catalytic reaction has 3 forms (E, ES, EP)
            # - 3 dead-end complexes with inhibitor bound to only one form
            # - 3 dead-end complexes with inhibitor bound to two forms
            # - 1 dead-end complex with inhibitor bound to all three forms
            # In sum, 1 catalytic + 2 with activator + 7 with dead-ends = 10

            #= More generally, each catalytic form is a slot. At each slot, you independently answer one yes/no question per regulator: "does this regulator form a dead-end here?"
            - 1 regulator → 1 yes/no question → 2 choices per slot
            - 2 regulators → 2 independent yes/no questions → 2²= 4 choices per slot
            - r regulators → 2^r choices per slot
            You make this choice independently at each of the n_cat slots, so you multiply:
            2^r × 2^r × … × 2^r (n_cat times) = (2^r)^n_cat
            For Uni-Uni with 1 regulator: 3 slots, 2 choices each → 2 × 2 × 2 = 8.
            The above include the case with no regulator.
            =#
            expected_n_cat_act_de=10,
            # TODO: SE/RE + equivalence forms:
            # catalytic 2^(3)-1 = 7
            # 1 essential activator 2^(4)-1 = 15
            # 1 non-essential activator 2^(9)-1 = 511
            # 1 inhibitor with 1 dead-end 2^(4)-1 = 15 (3x)
            # 1 inhibitor with 2 dead-end 2^(6)-1 = 63 (3x)
            # 1 inhibitor with 3 dead-end 2^(8)-1 = 255
            expected_n_total=1779,
            max_enumeration_time=5.0,
        ))
    end

    # 3. Uni-Uni + 2 regulators (chain dead-ends)
    let
        rxn = @enzyme_reaction begin
            substrates:S[C]
            products:P[C]
            regulators:R1[N], R2[P]
        end
        push!(specs, EnumerationTestSpec(
            name="Uni-Uni 2 Regulators",
            reaction=rxn,
            max_forms=100,
            # n_cat_forms × 2^regulators = 3 × 2^2 = 12
            expected_n_forms=12,
            expected_n_catalytic=1,
            # With 2 regulators, we get 1 catalytic
            # + 2x2=4 with each activators binding alone (like in Uni-Uni + 1 regulator case)
            # + 2 with one activators essential and one non-essential and vice versa
            # + 2 with both activators essential or non-essential
            expected_n_cat_with_act=9,
            # For each mechanism with no deadend the number of deadend mechanism can be
            # combinatorial calculated as (2^n_regultors)^n_cat_enz_forms.
            # For catalytic = (2^2)^3 = 64
            # For 1 essential activator: (2^1)^4 = 16 (2x since either reg can be inh or act)
            # For 1 non-essential activator: (2^1)^6 = 64 (also 2x)
            # No deadend complex with 2 activators since regulator can be either act and inh.
            # In sum, 64 + 16*2 + 64*2 + 4 = 228
            expected_n_cat_act_de=228,
            expected_n_total=24646535,
            max_enumeration_time=10.0,
        ))
    end

    # 4. Uni-Bi + 1 regulator
    let
        rxn = @enzyme_reaction begin
            substrates:A[C2]
            products:P1[C], P2[C]
            regulators:R[N]
        end
        push!(specs, EnumerationTestSpec(
            name="Uni-Bi",
            reaction=rxn,
            max_forms=100,
            expected_n_forms=16,
            # 2x sequencial product release and 1 random order → 3 catalytic topologies
            # For sequencial release: E+S <-> ES <-> EP1P2 <-> EP1 <-> E+P1 and vice versa
            # For random order:
            # E+S <-> ES <-> EP1P2 <-> EP1 + P2 <-> E+P1
            #                  ↕
            #                EP2 + P1 <-> E+P2
            expected_n_catalytic=3,
            # For sequencial release there are 2 activator mechanisms (as in Uni-Uni case)
            # and for random order there are 4 (2 per catalytic cycle)
            expected_n_cat_with_act=9,
            # For sequential release: (2^1)^4 = 16 (2x)
            # For random order: (2^1)^5 = 32
            # No deadend for activators mechanisms since regulator can be either act and inh.
            # In sum, 16*2 + 32 + 6 (activator mechanisms)= 70
            expected_n_cat_act_de=70,
            expected_n_total=435521,
            max_enumeration_time=5.0,
        ))
    end
    # 5. Bi-Bi + 1 regulator (same atoms on both substrates)
    let
        rxn = @enzyme_reaction begin
            substrates:A[C], B[C]
            products:P[C], Q[C]
            regulators:R[N]
        end
        push!(specs, EnumerationTestSpec(
            name="Bi-Bi + 1 Regulator",
            reaction=rxn,
            max_forms=100,
            expected_n_forms=22,
            expected_n_catalytic=9,
            expected_n_cat_with_act=27,
            expected_n_cat_act_de=530,
            skip_ress_test=true,
            expected_n_total=114684452,
        ))
    end

    # 6. Bi-Bi + 1 regulator (different atoms)
    let
        rxn = @enzyme_reaction begin
            substrates:A[C], B[N]
            products:P[C], Q[N]
            regulators:I[P2]
        end
        push!(specs, EnumerationTestSpec(
            name="Bi-Bi 1 Regulator",
            reaction=rxn,
            max_forms=100,
            expected_n_forms=22,
            expected_n_catalytic=9,
            expected_n_cat_with_act=27,
            expected_n_cat_act_de=530,
            skip_ress_test=true,
            expected_n_total=114684452,
        ))
    end

    # 7. Bi-Bi Ping Pong + 1 regulator
    let
        rxn = @enzyme_reaction begin
            substrates:A[CX], B[N]
            products:P[C], Q[NX]
            regulators:I[P2]
        end
        push!(specs, EnumerationTestSpec(
            name="Bi-Bi PP 1 Regulator",
            reaction=rxn,
            max_forms=100,
            expected_n_forms=34,
            # 9 standard + 16 ping-pong catalytic topologies
            expected_n_catalytic=25,
            expected_n_cat_with_act=75,
            expected_n_cat_act_de=2834,
            skip_ress_test=true,
            expected_n_total=56395770327,
        ))
    end

    # 8. Bi-Bi budget filtering (max_forms=5 restricts to single cycles)
    let
        rxn = @enzyme_reaction begin
            substrates:A[C], B[N]
            products:P[C], Q[N]
        end
        push!(specs, EnumerationTestSpec(
            name="Bi-Bi Budget",
            reaction=rxn,
            max_forms=5,
            expected_n_forms=11,
            expected_n_catalytic=4,
            expected_n_cat_with_act=4,
            expected_n_cat_act_de=4,
            expected_n_total=124,
            max_enumeration_time=5.0,
        ))
    end

    return specs
end

const ENUMERATION_TEST_SPECS = build_enumeration_test_specs()
