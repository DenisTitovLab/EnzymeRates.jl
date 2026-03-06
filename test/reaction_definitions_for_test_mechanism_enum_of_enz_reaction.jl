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
Compute RE group partition (test-local reimplementation).
Returns sorted vector of sorted form-index vectors.
"""
function _compute_re_partition_test(edges, eq_steps)
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
    groups = Dict{Int, Vector{Int}}()
    for i in form_indices
        r = find(i)
        push!(get!(groups, r, Int[]), i)
    end
    sort!([sort!(v) for v in values(groups)])
end

const _TestMONO = Vector{Pair{Symbol,Int}}

"""Increment exponent of `met` in a sorted monomial (test-local)."""
function _add_met_test(mono::_TestMONO, met::Symbol)::_TestMONO
    result = copy(mono)
    idx = findfirst(p -> p.first == met, result)
    if idx !== nothing
        result[idx] = met => result[idx].second + 1
    else
        push!(result, met => 1)
        sort!(result; by=first)
    end
    result
end

"""Multiply two sorted monomials (test-local)."""
function _mono_mul_test(a::_TestMONO, b::_TestMONO)::_TestMONO
    d = Dict{Symbol,Int}(p for p in a)
    for (s, c) in b
        d[s] = get(d, s, 0) + c
    end
    sort!([s => c for (s, c) in d]; by=first)
end

"""Spanning arborescence monomials (test-local)."""
function _spanning_arborescence_monomials_test(
    G::Int, R_conc::Dict{Tuple{Int,Int}, Set{_TestMONO}},
    root::Int,
)
    result = Set{_TestMONO}()
    non_root = [g for g in 1:G if g != root]
    isempty(non_root) && return Set{_TestMONO}([_TestMONO()])
    function enumerate_trees!(idx, current_mono)
        if idx > length(non_root)
            push!(result, current_mono)
            return
        end
        node = non_root[idx]
        for dst in 1:G
            dst == node && continue
            edge_monos = get(R_conc, (node, dst), nothing)
            edge_monos === nothing && continue
            for emono in edge_monos
                enumerate_trees!(
                    idx + 1, _mono_mul_test(current_mono, emono))
            end
        end
    end
    enumerate_trees!(1, _TestMONO())
    result
end

"""
Compute concentration fingerprint (test-local reimplementation).
Uses forms and adj_info to classify edges independently.
"""
function _concentration_fingerprint_test(edges, eq_steps, forms,
                                         adj_info)
    partition = _compute_re_partition_test(edges, eq_steps)
    G = length(partition)
    group_set = [Set(g) for g in partition]
    form_to_group = Dict(
        i => g for (g, grp) in enumerate(partition) for i in grp)

    alpha_conc = Dict{Int, _TestMONO}()
    for (g, group) in enumerate(partition)
        ref = group[1]
        alpha_conc[ref] = _TestMONO()
        queue = [ref]
        while !isempty(queue)
            cur = popfirst!(queue)
            for (idx, (a, b)) in enumerate(edges)
                eq_steps[idx] || continue
                neighbor = if a == cur && b in group_set[g]
                    b
                elseif b == cur && a in group_set[g]
                    a
                else
                    nothing
                end
                (neighbor === nothing ||
                    haskey(alpha_conc, neighbor)) && continue
                info = adj_info[minmax(a, b)]
                parent_mono = alpha_conc[cur]
                child_mono = if info.metabolite !== nothing
                    is_binding =
                        (info.type == :binding && a == cur) ||
                        (info.type == :release && b == cur)
                    is_binding ?
                        _add_met_test(parent_mono,
                            info.metabolite) :
                        copy(parent_mono)
                else
                    copy(parent_mono)
                end
                alpha_conc[neighbor] = child_mono
                push!(queue, neighbor)
            end
        end
    end

    sigma_conc = [Set{_TestMONO}(alpha_conc[i] for i in group)
                  for group in partition]

    R_conc = Dict{Tuple{Int,Int}, Set{_TestMONO}}()
    for (idx, (a, b)) in enumerate(edges)
        eq_steps[idx] && continue
        info = adj_info[minmax(a, b)]
        g1, g2 = form_to_group[a], form_to_group[b]
        g1 == g2 && continue
        fwd_met = (info.type == :binding) ? info.metabolite :
            nothing
        fwd_mono = fwd_met !== nothing ?
            _add_met_test(alpha_conc[a], fwd_met) :
            copy(alpha_conc[a])
        push!(get!(R_conc, (g1, g2), Set{_TestMONO}()), fwd_mono)
        rev_met = (info.type == :release) ? info.metabolite :
            nothing
        rev_mono = rev_met !== nothing ?
            _add_met_test(alpha_conc[b], rev_met) :
            copy(alpha_conc[b])
        push!(get!(R_conc, (g2, g1), Set{_TestMONO}()), rev_mono)
    end

    fingerprint = Set{_TestMONO}()
    for g in 1:G
        D_g = _spanning_arborescence_monomials_test(G, R_conc, g)
        for s in sigma_conc[g], d in D_g
            push!(fingerprint, _mono_mul_test(s, d))
        end
    end
    fingerprint
end

"""
Compute constraint descriptor (test-local reimplementation).
"""
function _constraint_descriptor_test(edges, adj_info, eq_steps,
                                      valid_groups, constraint_mask)
    descriptor = Set{Tuple{Symbol, Symbol}}()
    for (gi, g) in enumerate(valid_groups)
        (constraint_mask >> (gi - 1)) & 1 == 1 || continue
        met = adj_info[minmax(edges[g[1]]...)].metabolite
        mode = eq_steps[g[1]] ? :RE : :SS
        push!(descriptor, (met, mode))
    end
    descriptor
end

"""Map each dead-end edge to its catalytic counterpart index.
Test-local version using only public struct fields."""
function _dead_end_catalytic_map_test(edges, n_cat, forms)
    _strip_regs(fi) = Tuple(
        s.role == :reg ? nothing : s.atoms
        for s in forms[fi].sites)
    _key(a, b) = let sa = _strip_regs(a), sb = _strip_regs(b)
        hash(sa) <= hash(sb) ? (sa, sb) : (sb, sa)
    end
    cat_stripped = Dict{Tuple, Int}()
    for i in 1:n_cat
        a, b = edges[i]
        cat_stripped[_key(a, b)] = i
    end
    [get(cat_stripped, _key(a, b), nothing)
     for (a, b) in @view edges[n_cat+1:end]]
end

"""
Compute the RE/SS + constraint count using brute-force enumeration
with G ≤ max_re_groups cap, with fingerprint-based deduplication.

Baseline: first isomerization edge is always SS. Iterates over
subsets of remaining edges to make SS, keeps masks with 2 ≤ G ≤ 7,
and counts unique (fingerprint, constraint_descriptor) keys.
Dead-end edges inherit RE/SS from their catalytic counterpart.
"""
function _compute_expected_n_total(
    spec::EnzymeRates.MechanismSpec,
    forms::Vector{EnzymeRates.EnzymeFormSpec};
    max_re_groups::Int=7,
)
    edges = spec.edges
    n = length(edges)
    n == 0 && return 0
    n_cat = spec.n_catalytic_edges
    iso_idx = _find_first_isomerization_test(edges, forms)
    de_cat_map = _dead_end_catalytic_map_test(edges, n_cat, forms)
    equiv_groups = _find_equiv_groups(edges, forms, n_cat,
        de_cat_map)

    # Build adjacency for fingerprint computation
    adj_info = EnzymeRates._build_adjacency(forms)

    # Map dead-end edges to catalytic counterparts
    de_cat_map = _dead_end_catalytic_map_test(edges, n_cat, forms)

    other_indices = [i for i in 1:n_cat if i != iso_idx]
    n_other = length(other_indices)

    seen = Set{Tuple{Set{_TestMONO},
                      Set{Tuple{Symbol,Symbol}}}}()
    for ss_mask in 0:(1 << n_other) - 1
        eq_steps = fill(true, n)
        eq_steps[iso_idx] = false
        for (bit, idx) in enumerate(other_indices)
            (ss_mask >> (bit - 1)) & 1 == 1 &&
                (eq_steps[idx] = false)
        end
        # Propagate to dead-end edges
        for (di, cat_idx) in enumerate(de_cat_map)
            cat_idx === nothing && continue
            eq_steps[n_cat + di] = eq_steps[cat_idx]
        end
        G = _compute_re_group_count_test(edges, eq_steps)
        (G < 2 || G > max_re_groups) && continue
        fp = _concentration_fingerprint_test(
            edges, eq_steps, forms, adj_info)
        valid_groups = [g for g in equiv_groups
            if all(eq_steps[s] == eq_steps[g[1]] for s in g)]
        for constraint_mask in 0:(1 << length(valid_groups)) - 1
            desc = _constraint_descriptor_test(
                edges, adj_info, eq_steps, valid_groups,
                constraint_mask)
            push!(seen, (fp, desc))
        end
    end
    length(seen)
end

"""Classify a single edge by site diffs (test-local).
Returns (metabolite, is_single_site_binding) or nothing."""
function _classify_edge_test(edges, forms, i)
    (a, b) = edges[i]
    diff_count = 0
    diff_k = 0
    for k in 1:length(forms[a].sites)
        if forms[a].sites[k].atoms != forms[b].sites[k].atoms
            diff_count += 1
            diff_count == 1 && (diff_k = k)
        end
    end
    diff_count == 1 || return nothing
    k = diff_k
    a_occ = forms[a].sites[k].atoms !== nothing
    site = a_occ ? forms[a].sites[k] : forms[b].sites[k]
    site.role == :prod && return nothing
    return site.metabolite
end

"""Find equivalent groups from catalytic + dead-end edges:
non-product binding edges grouped by metabolite.
Uses only public struct fields."""
function _find_equiv_groups(edges, forms, n_catalytic_edges,
                            de_cat_map=nothing)
    binding_key = Dict{Symbol,Vector{Int}}()
    for i in 1:n_catalytic_edges
        met = _classify_edge_test(edges, forms, i)
        met === nothing && continue
        push!(get!(binding_key, met, Int[]), i)
    end
    if de_cat_map !== nothing
        for (di, cat_idx) in enumerate(de_cat_map)
            cat_idx === nothing && continue
            edge_idx = n_catalytic_edges + di
            met = _classify_edge_test(edges, forms, edge_idx)
            met === nothing && continue
            push!(get!(binding_key, met, Int[]), edge_idx)
        end
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

"""Stirling-based partition multiplicity count (test-local)."""
function _partition_mult_count_test(k, N)
    k == 0 && return 1
    S = zeros(Int, k, k)
    S[1, 1] = 1
    for n in 2:k
        for g in 1:n
            S[n, g] = (g > 1 ? S[n-1, g-1] : 0) + g * S[n-1, g]
        end
    end
    sum(S[k, g] * N^g for g in 1:k)
end

"""
Compute expected total with oligomeric expansion (catalytic_n=N).

For each dead-end spec with k allosteric regulators:
  EM contribution: ress_count
  OEM contribution: ress_count * _partition_mult_count(k, N)
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
        oem = ress * _partition_mult_count_test(k, catalytic_n)
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
            # After fingerprint dedup: all 3 RE/SS masks produce
            # the same concentration monomial set {1, S, P}
            expected_n_total=1,
            max_enumeration_time=5.0,
        ))
    end

    # 2. Uni-Uni + 1 regulator
    #
    # Forms (6): E, ES, EP  +  ER, ESR, EPR
    #   (2×2×2=8 minus 2 excluded SP combos)
    # Adjacency (9 edges):
    #   Catalytic:  (E,ES) bind S, (E,EP) bind P, (ES,EP) iso
    #   R-binding:  (E,ER), (ES,ESR), (EP,EPR)
    #   Inter-de:   (ER,ESR) bind S, (ER,EPR) bind P, (ESR,EPR) iso
    # Catalytic cycle: 1 (E→ES→EP→E), 3 catalytic edges.
    #
    # Dead-end stage (9 specs from 2 regulator partitions):
    #   Partition {R dead-end}: each of 3 catalytic forms independently
    #     can have R bound → 2^3 = 8 dead-end configs
    #   Partition {R allosteric}: 1 (bare catalytic topology)
    #
    # RE/SS stage: only the 3 catalytic edges are candidates for SS.
    # Dead-end substrate/product-binding edges (ER↔ESR, ER↔EPR)
    # inherit RE/SS from their catalytic counterpart (E↔ES, E↔EP).
    # Dead-end regulator-binding edges (E↔ER etc.) remain always RE.
    # Dead-end S-binding edges join equiv groups with their catalytic
    # counterpart, adding constrained variants (K_S_de = K_S).
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
            expected_n_total=19,
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
            expected_n_total=247,
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
            expected_n_total=743,
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
            expected_n_total=48607,
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
            expected_n_total=48607,
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
            expected_n_total=264,
            max_enumeration_time=10.0,
        ))
    end

    return specs
end

const ENUMERATION_TEST_SPECS = build_enumeration_test_specs()
