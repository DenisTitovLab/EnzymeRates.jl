# ABOUTME: Staged pipeline for enumerating enzyme mechanism topologies.
# ABOUTME: Converts EnzymeReaction → MechanismSpec/AllostericMechanismSpec via 8 stages.

# ─── Data Types ─────────────────────────────────────────────

"""
Constraint on kinetic parameters: target parameter equals a linear
combination of source parameters.
Format: `(target_sym, coeff, [(src_sym, src_coeff), ...])`.
"""
const ParamConstraint = Tuple{Symbol, Int, Vector{Tuple{Symbol, Int}}}

# ─── Mechanism Spec Types ──────────────────────────────────────

abstract type AbstractMechanismSpec end

"""Elementary step in canonical binding direction (metabolite on LHS)."""
struct StepSpec
    reactants::Vector{Symbol}   # [:E, :S] or [:EAB]
    products::Vector{Symbol}    # [:ES] or [:EPQ]
    is_equilibrium::Bool
end

Base.:(==)(a::StepSpec, b::StepSpec) =
    a.reactants == b.reactants &&
    a.products == b.products &&
    a.is_equilibrium == b.is_equilibrium

Base.hash(s::StepSpec, h::UInt) =
    hash(s.is_equilibrium,
        hash(s.products, hash(s.reactants, h)))

"""
    MechanismSpec <: AbstractMechanismSpec

Represents a monomeric enzyme mechanism specification in the
staged enumeration pipeline. Steps are `StepSpec` values with
inline form names and equilibrium status.
"""
struct MechanismSpec <: AbstractMechanismSpec
    reaction::Any
    steps::Vector{StepSpec}
    param_constraints::Vector{ParamConstraint}
    param_count::Int
end

"""
    AllostericMechanismSpec <: AbstractMechanismSpec

Represents an allosteric enzyme mechanism built from a
base `MechanismSpec` plus allosteric site and multiplicity info.
`tr_equiv_metabolites` lists metabolites with K_T = K_R.
"""
struct AllostericMechanismSpec <: AbstractMechanismSpec
    base::MechanismSpec
    catalytic_n::Int
    allosteric_reg_sites::Vector{Vector{Symbol}}
    allosteric_multiplicities::Vector{Int}
    tr_equiv_metabolites::Vector{Symbol}
end

# ─── StepSpec Helpers ──────────────────────────────────────────

"""Return the metabolite for a step, or nothing for isomerization."""
step_metabolite(s::StepSpec) =
    length(s.reactants) == 2 ? s.reactants[2] : nothing

"""Return (from_form, to_form) for a step."""
step_forms(s::StepSpec) = (s.reactants[1], s.products[1])

"""Collect all unique form names from steps."""
function all_form_names(spec::MechanismSpec)
    forms = Set{Symbol}()
    for s in spec.steps
        push!(forms, s.reactants[1])
        push!(forms, s.products[1])
    end
    forms
end

function all_form_names(steps::Vector{StepSpec})
    forms = Set{Symbol}()
    for s in steps
        push!(forms, s.reactants[1])
        push!(forms, s.products[1])
    end
    forms
end

# ─── MechanismIterator ────────────────────────────────────────

struct MechanismIterator
    inner::Any
    total::Int
end

Base.eltype(::Type{MechanismIterator}) = AbstractMechanismSpec
Base.IteratorSize(::Type{MechanismIterator}) = Base.HasLength()
Base.length(iter::MechanismIterator) = iter.total
Base.iterate(iter::MechanismIterator, s...) = iterate(iter.inner, s...)

# ─── RE Partition (union-find) ─────────────────────────────────

"""
    _compute_re_partition_from_steps(steps) -> Vector{Vector{Symbol}}

Connected components of enzyme forms linked by RE steps.
"""
function _compute_re_partition_from_steps(
    steps::Vector{StepSpec},
)
    form_names = collect(all_form_names(steps))
    parent = Dict(f => f for f in form_names)
    function find(x)
        while parent[x] != x
            parent[x] = parent[parent[x]]
            x = parent[x]
        end
        x
    end
    for s in steps
        s.is_equilibrium || continue
        a, b = s.reactants[1], s.products[1]
        ra, rb = find(a), find(b)
        ra != rb && (parent[ra] = rb)
    end
    groups = Dict{Symbol, Vector{Symbol}}()
    for f in form_names
        r = find(f)
        push!(get!(groups, r, Symbol[]), f)
    end
    sort!([sort!(v) for v in values(groups)])
end

# ─── Concentration Fingerprint ─────────────────────────────────

"""Increment the exponent of `met` in a sorted monomial."""
function _add_met(mono::MONO, met::Symbol)::MONO
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

"""
Enumerate all spanning arborescences of a G-node directed graph
rooted at `root`, returning the set of concentration monomials.
"""
function _spanning_arborescence_monomials(
    G::Int, R_conc::Dict{Tuple{Int,Int}, Set{MONO}}, root::Int,
)
    result = Set{MONO}()
    non_root = [g for g in 1:G if g != root]
    isempty(non_root) && return Set{MONO}([MONO()])
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
                    idx + 1, _mono_mul(current_mono, emono))
            end
        end
    end
    enumerate_trees!(1, MONO())
    result
end

"""
    _concentration_fingerprint(steps, partition) → Set{MONO}

Concentration monomials in the rate equation denominator.
Uses step-based representation with Symbol-named forms.
Binding direction is always forward (reactant→product) since
steps are in canonical form.
"""
function _concentration_fingerprint(
    steps::Vector{StepSpec},
    partition::Vector{Vector{Symbol}},
)
    G = length(partition)
    group_set = [Set(g) for g in partition]
    form_to_group = Dict(
        f => g
        for (g, grp) in enumerate(partition)
        for f in grp)

    # BFS within each RE component to compute
    # relative concentration monomials
    alpha_conc = Dict{Symbol, MONO}()
    for (g, group) in enumerate(partition)
        ref = group[1]
        alpha_conc[ref] = MONO()
        queue = [ref]
        while !isempty(queue)
            cur = popfirst!(queue)
            for s in steps
                s.is_equilibrium || continue
                from, to = step_forms(s)
                met = step_metabolite(s)
                # Check both directions of the edge
                neighbor = if from == cur &&
                        to in group_set[g]
                    to
                elseif to == cur &&
                        from in group_set[g]
                    from
                else
                    nothing
                end
                (neighbor === nothing ||
                    haskey(alpha_conc, neighbor)) &&
                    continue
                parent_mono = alpha_conc[cur]
                child_mono = if met !== nothing
                    # Forward = binding direction
                    if cur == from
                        _add_met(parent_mono, met)
                    else
                        copy(parent_mono)
                    end
                else
                    copy(parent_mono)
                end
                alpha_conc[neighbor] = child_mono
                push!(queue, neighbor)
            end
        end
    end

    sigma_conc = [Set{MONO}(alpha_conc[f] for f in group)
                  for group in partition]

    # Build SS inter-group transitions
    R_conc = Dict{Tuple{Int,Int}, Set{MONO}}()
    for s in steps
        s.is_equilibrium && continue
        from, to = step_forms(s)
        met = step_metabolite(s)
        g1, g2 = form_to_group[from], form_to_group[to]
        g1 == g2 && continue

        # Forward direction: from→to binds metabolite
        fwd_mono = met !== nothing ?
            _add_met(alpha_conc[from], met) :
            copy(alpha_conc[from])
        push!(get!(R_conc, (g1, g2), Set{MONO}()),
            fwd_mono)

        # Reverse direction: to→from unbinds metabolite
        rev_mono = copy(alpha_conc[to])
        push!(get!(R_conc, (g2, g1), Set{MONO}()),
            rev_mono)
    end

    fingerprint = Set{MONO}()
    for g in 1:G
        D_g = _spanning_arborescence_monomials(
            G, R_conc, g)
        for s in sigma_conc[g], d in D_g
            push!(fingerprint, _mono_mul(s, d))
        end
    end
    fingerprint
end

# ─── Equivalence Groups + Constraints ──────────────────────────

"""
    _constraint_descriptor(steps, valid_groups, constraint_mask)

Constraint descriptor for step-based representation.
"""
function _constraint_descriptor(
    steps::Vector{StepSpec},
    valid_groups::Vector{Vector{Int}},
    constraint_mask::Int,
)
    descriptor = Set{Tuple{Symbol, Symbol}}()
    for (gi, g) in enumerate(valid_groups)
        (constraint_mask >> (gi - 1)) & 1 == 1 || continue
        met = step_metabolite(steps[g[1]])
        mode = steps[g[1]].is_equilibrium ? :RE : :SS
        push!(descriptor, (met, mode))
    end
    descriptor
end

const _DedupKey = Tuple{Set{MONO}, Set{Tuple{Symbol,Symbol}}}

# ─── Set Partitions + Allosteric Helpers ─────────────────────

"""
    _set_partitions(elements::Vector{Symbol})

Enumerate all set partitions (Bell number partitions).
"""
function _set_partitions(elements::Vector{Symbol})
    n = length(elements)
    n == 0 && return [Vector{Symbol}[]]
    n == 1 && return [Vector{Symbol}[elements]]
    result = Vector{Vector{Vector{Symbol}}}()
    for partition in _set_partitions(elements[1:end-1])
        last_elem = elements[end]
        for i in eachindex(partition)
            new_part = [copy(g) for g in partition]
            push!(new_part[i], last_elem)
            push!(result, new_part)
        end
        push!(result, [partition; [Symbol[last_elem]]])
    end
    result
end

"""
    _partition_mult_count(k, N) → Int

Count allosteric multiplicity variants: sum_{g=1}^{k} S(k,g) * N^g.
"""
function _partition_mult_count(k::Int, N::Int)
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

# ─── Stage 1: Catalytic Topologies ───────────────────────────

"""Build a form name Symbol from sorted bound metabolite names."""
function _form_name(
    bound_subs::Vector{Symbol},
    bound_prods::Vector{Symbol},
    has_residual::Bool,
)
    parts = sort!(vcat(bound_subs, bound_prods))
    base = isempty(parts) ? "E" : "E_" * join(parts, "_")
    has_residual && isempty(parts) && (base = "Estar")
    has_residual && !isempty(parts) &&
        (base = "Estar_" * join(parts, "_"))
    Symbol(base)
end

"""Extract atom counts as Dict{Symbol,Int} for a metabolite."""
function _atoms_dict(
    @nospecialize(reaction::EnzymeReaction),
    met::Symbol,
)
    result = Dict{Symbol,Int}()
    for (name, atoms) in substrates(reaction)
        name == met || continue
        for (a, c) in atoms
            result[a] = get(result, a, 0) + c
        end
        return result
    end
    for (name, atoms) in products(reaction)
        name == met || continue
        for (a, c) in atoms
            result[a] = get(result, a, 0) + c
        end
        return result
    end
    result
end

"""Check if product atoms are a subset of accumulated atoms."""
function _can_pingpong(
    accumulated::Dict{Symbol,Int},
    prod_atoms::Dict{Symbol,Int},
)
    for (a, c) in prod_atoms
        get(accumulated, a, 0) < c && return false
    end
    true
end

"""Subtract atom counts: accumulated minus product atoms."""
function _subtract_atoms(
    accumulated::Dict{Symbol,Int},
    prod_atoms::Dict{Symbol,Int},
)
    result = copy(accumulated)
    for (a, c) in prod_atoms
        result[a] -= c
        result[a] == 0 && delete!(result, a)
    end
    result
end

"""Add atom counts: accumulated plus substrate atoms."""
function _add_atoms(
    accumulated::Dict{Symbol,Int},
    sub_atoms::Dict{Symbol,Int},
)
    result = copy(accumulated)
    for (a, c) in sub_atoms
        result[a] = get(result, a, 0) + c
    end
    result
end

"""
    _catalytic_topologies(reaction) -> Vector{MechanismSpec}

Build catalytic cycle topologies by constructive backtracking.
Each topology is a set of steps forming one or more complete
catalytic cycles (E -> ... -> E).
"""
function _catalytic_topologies(
    @nospecialize(reaction::EnzymeReaction),
)
    subs = substrates(reaction)
    prods = products(reaction)
    sub_names = Symbol[s[1] for s in subs]
    prod_names = Symbol[p[1] for p in prods]

    # Precompute atom dicts for each metabolite
    sub_atoms = Dict(
        s => _atoms_dict(reaction, s) for s in sub_names
    )
    prod_atoms = Dict(
        p => _atoms_dict(reaction, p) for p in prod_names
    )

    # Collect all complete catalytic paths as step lists
    all_paths = Vector{Vector{StepSpec}}()

    # Backtracking state:
    # - cur_form: current enzyme form name (Symbol)
    # - acc_atoms: atoms currently on the enzyme
    # - consumed_subs: substrates consumed so far (history)
    # - released_prods: products released so far (history)
    # - on_enzyme_subs: substrates currently bound
    # - on_enzyme_prods: products currently bound
    #     (post-final-isomerize)
    # - has_residual: enzyme carries leftover atoms from
    #     ping-pong
    # - post_final: in product-release phase after final
    #     isomerization
    # - steps: path of StepSpec accumulated so far
    function backtrack!(
        cur_form::Symbol,
        acc_atoms::Dict{Symbol,Int},
        consumed_subs::Vector{Symbol},
        released_prods::Vector{Symbol},
        on_enzyme_subs::Vector{Symbol},
        on_enzyme_prods::Vector{Symbol},
        has_residual::Bool,
        post_final::Bool,
        steps::Vector{StepSpec},
    )
        # Check for complete cycle
        if cur_form == :E && !isempty(steps)
            if Set(consumed_subs) == Set(sub_names) &&
                    Set(released_prods) == Set(prod_names)
                push!(all_paths, copy(steps))
                return
            end
        end

        remaining_subs = [
            s for s in sub_names if s ∉ consumed_subs
        ]
        remaining_prods = [
            p for p in prod_names if p ∉ released_prods
        ]

        if post_final
            # Release any currently bound product
            for p in copy(on_enzyme_prods)
                new_on_prods = filter(!=(p), on_enzyme_prods)
                new_released = [released_prods; p]
                new_form = _form_name(
                    Symbol[], new_on_prods, false
                )
                # Canonical: metabolite on LHS
                step = StepSpec(
                    [new_form, p], [cur_form], true
                )
                push!(steps, step)
                backtrack!(
                    new_form,
                    _subtract_atoms(
                        acc_atoms, prod_atoms[p]
                    ),
                    consumed_subs, new_released,
                    Symbol[], new_on_prods,
                    false, !isempty(new_on_prods),
                    steps
                )
                pop!(steps)
            end
            return
        end

        if isempty(on_enzyme_subs) && !has_residual
            # Free enzyme: bind any remaining substrate
            for s in remaining_subs
                new_on = [on_enzyme_subs; s]
                new_consumed = [consumed_subs; s]
                new_form = _form_name(
                    new_on, Symbol[], false
                )
                step = StepSpec(
                    [cur_form, s], [new_form], true
                )
                push!(steps, step)
                backtrack!(
                    new_form,
                    _add_atoms(acc_atoms, sub_atoms[s]),
                    new_consumed, released_prods,
                    new_on, Symbol[],
                    false, false, steps
                )
                pop!(steps)
            end
        elseif !isempty(on_enzyme_subs) && !has_residual
            # Substrates bound, no residual
            # Option 1: bind another substrate
            for s in remaining_subs
                new_on = [on_enzyme_subs; s]
                new_consumed = [consumed_subs; s]
                new_form = _form_name(
                    new_on, Symbol[], false
                )
                step = StepSpec(
                    [cur_form, s], [new_form], true
                )
                push!(steps, step)
                backtrack!(
                    new_form,
                    _add_atoms(acc_atoms, sub_atoms[s]),
                    new_consumed, released_prods,
                    new_on, Symbol[],
                    false, false, steps
                )
                pop!(steps)
            end
            # Option 2: ping-pong isomerize
            if !isempty(remaining_subs)
                for p in remaining_prods
                    _can_pingpong(
                        acc_atoms, prod_atoms[p]
                    ) || continue
                    residual = _subtract_atoms(
                        acc_atoms, prod_atoms[p]
                    )
                    # Genuine ping-pong needs nonzero
                    # residual
                    isempty(residual) && continue
                    iso_form = _form_name(
                        on_enzyme_subs, [p], true
                    )
                    step = StepSpec(
                        [cur_form], [iso_form], true
                    )
                    push!(steps, step)
                    rel_form = _form_name(
                        Symbol[], Symbol[], true
                    )
                    rel_step = StepSpec(
                        [rel_form, p], [iso_form], true
                    )
                    push!(steps, rel_step)
                    backtrack!(
                        rel_form, residual,
                        consumed_subs,
                        [released_prods; p],
                        Symbol[], Symbol[],
                        true, false, steps
                    )
                    pop!(steps)
                    pop!(steps)
                end
            end
            # Option 3: final isomerize (all subs bound)
            if isempty(remaining_subs)
                all_prod_atoms = reduce(
                    _add_atoms, values(prod_atoms);
                    init=Dict{Symbol,Int}()
                )
                if _can_pingpong(
                    acc_atoms, all_prod_atoms
                )
                    new_form = _form_name(
                        Symbol[],
                        copy(remaining_prods), false
                    )
                    step = StepSpec(
                        [cur_form], [new_form], true
                    )
                    push!(steps, step)
                    backtrack!(
                        new_form, acc_atoms,
                        consumed_subs, released_prods,
                        Symbol[],
                        copy(remaining_prods),
                        false, true, steps
                    )
                    pop!(steps)
                end
            end
        elseif isempty(on_enzyme_subs) && has_residual
            # Residual only (E*): bind any remaining sub
            for s in remaining_subs
                new_on = [s]
                new_consumed = [consumed_subs; s]
                new_form = _form_name(
                    new_on, Symbol[], true
                )
                step = StepSpec(
                    [cur_form, s], [new_form], true
                )
                push!(steps, step)
                backtrack!(
                    new_form,
                    _add_atoms(acc_atoms, sub_atoms[s]),
                    new_consumed, released_prods,
                    new_on, Symbol[],
                    true, false, steps
                )
                pop!(steps)
            end
        elseif !isempty(on_enzyme_subs) && has_residual
            # Residual + substrates bound
            # Option 1: bind another remaining substrate
            for s in remaining_subs
                new_on = [on_enzyme_subs; s]
                new_consumed = [consumed_subs; s]
                new_form = _form_name(
                    new_on, Symbol[], true
                )
                step = StepSpec(
                    [cur_form, s], [new_form], true
                )
                push!(steps, step)
                backtrack!(
                    new_form,
                    _add_atoms(acc_atoms, sub_atoms[s]),
                    new_consumed, released_prods,
                    new_on, Symbol[],
                    true, false, steps
                )
                pop!(steps)
            end
            # Option 2: isomerize to release a product
            for p in remaining_prods
                _can_pingpong(
                    acc_atoms, prod_atoms[p]
                ) || continue
                residual_atoms = _subtract_atoms(
                    acc_atoms, prod_atoms[p]
                )
                has_more = !isempty(residual_atoms)
                if has_more ||
                        !isempty(remaining_subs)
                    # Ping-pong: isomerize + release
                    iso_form = _form_name(
                        on_enzyme_subs, [p], has_more
                    )
                    step = StepSpec(
                        [cur_form], [iso_form], true
                    )
                    push!(steps, step)
                    rel_form = _form_name(
                        Symbol[], Symbol[], has_more
                    )
                    rel_step = StepSpec(
                        [rel_form, p],
                        [iso_form], true
                    )
                    push!(steps, rel_step)
                    backtrack!(
                        rel_form, residual_atoms,
                        consumed_subs,
                        [released_prods; p],
                        Symbol[], Symbol[],
                        has_more, false, steps
                    )
                    pop!(steps)
                    pop!(steps)
                end
                if !has_more &&
                        isempty(remaining_subs)
                    # Final: all subs consumed, all
                    # residual consumed — release all
                    # remaining products
                    new_form = _form_name(
                        Symbol[],
                        copy(remaining_prods), false
                    )
                    step = StepSpec(
                        [cur_form], [new_form], true
                    )
                    push!(steps, step)
                    backtrack!(
                        new_form, acc_atoms,
                        consumed_subs, released_prods,
                        Symbol[],
                        copy(remaining_prods),
                        false, true, steps
                    )
                    pop!(steps)
                end
            end
        end
    end

    backtrack!(
        :E, Dict{Symbol,Int}(), Symbol[], Symbol[],
        Symbol[], Symbol[], false, false, StepSpec[]
    )

    isempty(all_paths) && return MechanismSpec[]

    # Deduplicate paths by their step content (as sets)
    StepKey = Tuple{Vector{Symbol}, Vector{Symbol}}
    _step_key(s::StepSpec) = (
        sort(s.reactants), sort(s.products)
    )::StepKey
    unique_paths = Vector{Vector{StepSpec}}()
    seen_path_keys = Set{Set{StepKey}}()
    for path in all_paths
        key = Set(_step_key(s) for s in path)
        key ∈ seen_path_keys && continue
        push!(seen_path_keys, key)
        push!(unique_paths, path)
    end

    # Partition paths into ping-pong (has Estar forms)
    # and sequential groups. Only combine paths within
    # the same group to avoid biochemically unrealistic
    # mixed topologies.
    function _is_pingpong_path(path)
        for s in path
            for sym in Iterators.flatten(
                (s.reactants, s.products)
            )
                startswith(string(sym), "Estar") &&
                    return true
            end
        end
        false
    end
    pp_paths = [p for p in unique_paths
                if _is_pingpong_path(p)]
    seq_paths = [p for p in unique_paths
                 if !_is_pingpong_path(p)]

    # Build step-set unions incrementally within each
    # group. For each path, union its steps with all
    # existing step-sets and add any new results.
    known_combos = Dict{Set{StepKey},
                        Dict{StepKey, StepSpec}}()
    for group in (pp_paths, seq_paths)
        isempty(group) && continue
        grp_keys = [
            Set(_step_key(s) for s in p)
            for p in group
        ]
        grp_dicts = [
            Dict(_step_key(s) => s for s in p)
            for p in group
        ]
        grp_combos = Dict{Set{StepKey},
                          Dict{StepKey, StepSpec}}()
        for (i, pkeys) in enumerate(grp_keys)
            if !haskey(grp_combos, pkeys)
                grp_combos[pkeys] = copy(grp_dicts[i])
            end
            new_entries = Dict{Set{StepKey},
                               Dict{StepKey, StepSpec}}()
            for (ks, sd) in grp_combos
                merged_keys = union(ks, pkeys)
                merged_keys == ks && continue
                haskey(grp_combos, merged_keys) &&
                    continue
                haskey(new_entries, merged_keys) &&
                    continue
                merged = copy(sd)
                merge!(merged, grp_dicts[i])
                new_entries[merged_keys] = merged
            end
            merge!(grp_combos, new_entries)
        end
        merge!(known_combos, grp_combos)
    end

    # Build MechanismSpec for each topology
    result = MechanismSpec[]
    for step_dict in values(known_combos)
        steps = collect(values(step_dict))
        # Sort steps: binding first, then isomerization
        sort!(steps; by=s -> (
            length(s.reactants) == 1 ? 1 : 0,
            join(sort(s.reactants), "_")
        ))

        n_steps = length(steps)
        # Default RE/SS: first isomerization is SS, rest RE
        iso_idx = findfirst(
            s -> length(s.reactants) == 1, steps
        )
        tagged = [
            StepSpec(
                s.reactants, s.products,
                i != iso_idx
            )
            for (i, s) in enumerate(steps)
        ]

        # Compute param_count
        form_names = Set{Symbol}()
        for s in tagged
            union!(form_names, s.reactants)
            union!(form_names, s.products)
        end
        # Remove metabolite names from form set
        met_names = Set{Symbol}(sub_names)
        union!(met_names, prod_names)
        setdiff!(form_names, met_names)
        n_forms = length(form_names)
        n_independent_cycles = n_steps - n_forms + 1
        n_thermo = n_independent_cycles
        n_re = n_steps - 1  # all except the one SS
        n_ss = 1
        param_count = n_re + 2 * n_ss - n_thermo + 2

        push!(result, MechanismSpec(
            reaction, tagged, ParamConstraint[],
            param_count
        ))
    end
    result
end

# ─── Stage 2: RE/SS Assignment ───────────────────────────────

"""
    _expand_ress_variants(specs, reaction; max_re_groups=7)
        -> Vector{MechanismSpec}

Enumerate all RE/SS assignment combinations for mechanism
steps. Every step can independently be RE or SS. The all-RE
assignment (mask=0) is excluded because at least one step
must be SS for the King-Altman method.
"""
function _expand_ress_variants(
    specs::Vector{MechanismSpec},
    @nospecialize(reaction::EnzymeReaction);
    max_re_groups::Int=7,
)
    result = MechanismSpec[]
    for spec in specs
        n = length(spec.steps)
        n == 0 && continue

        for mask in 1:(1 << n) - 1
            steps = [
                StepSpec(
                    s.reactants, s.products,
                    (mask >> (i - 1)) & 1 == 0,
                )
                for (i, s) in enumerate(spec.steps)
            ]
            # At least one step must be SS
            any(!s.is_equilibrium for s in steps) ||
                continue
            # Check RE-connected groups ≤ max_re_groups
            partition =
                _compute_re_partition_from_steps(steps)
            length(partition) > max_re_groups && continue

            n_re = count(s.is_equilibrium for s in steps)
            n_ss = n - n_re
            n_forms = length(all_form_names(steps))
            n_thermo = n - n_forms + 1
            pc = n_re + 2 * n_ss - n_thermo + 2

            push!(result, MechanismSpec(
                spec.reaction, steps,
                spec.param_constraints, pc,
            ))
        end
    end
    result
end

# ─── Stage 2.5: Substrate/Product Dead-End Expansion ──────────

"""
    _bound_metabolites_at_forms(spec, reaction)

Map each form name to its set of bound metabolite names by
tracing binding steps from :E. Isomerization steps swap
all substrates for all products.
"""
function _bound_metabolites_at_forms(
    spec::MechanismSpec,
    @nospecialize(reaction::EnzymeReaction),
)
    sub_names = Set(s[1] for s in substrates(reaction))
    prod_names = Set(p[1] for p in products(reaction))
    bound = Dict{Symbol, Set{Symbol}}(:E => Set{Symbol}())
    changed = true
    while changed
        changed = false
        for s in spec.steps
            from, to = s.reactants[1], s.products[1]
            met = step_metabolite(s)
            if met !== nothing
                # Binding step: to = from + met
                if haskey(bound, from) &&
                        !haskey(bound, to)
                    bound[to] = union(
                        bound[from], Set([met]))
                    changed = true
                end
                if haskey(bound, to) &&
                        !haskey(bound, from)
                    bound[from] = setdiff(
                        bound[to], Set([met]))
                    changed = true
                end
            else
                # Isomerization: all subs become all prods
                if haskey(bound, from) &&
                        !haskey(bound, to)
                    bound[to] = union(
                        setdiff(bound[from], sub_names),
                        prod_names)
                    changed = true
                end
                if haskey(bound, to) &&
                        !haskey(bound, from)
                    bound[from] = union(
                        setdiff(bound[to], prod_names),
                        sub_names)
                    changed = true
                end
            end
        end
    end
    bound
end

"""
    _dead_end_form_name(base_bound, added_met)

Create form name for a dead-end form: base form's bound
metabolites plus the added metabolite.
"""
function _dead_end_form_name(
    base_bound::Set{Symbol}, added_met::Symbol,
)
    all_mets = sort(collect(
        union(base_bound, Set([added_met]))))
    Symbol("E_" * join(all_mets, "_"))
end

"""
    _expand_substrate_product_dead_ends(specs, reaction)
        -> Vector{MechanismSpec}

For each spec, enumerate substrate/product dead-end
form combinations. A dead-end form is created when a
substrate or product binds to a catalytic form where it
doesn't normally bind, subject to:
- The resulting form is not already a catalytic form
- The resulting form doesn't have all substrates + any
  product, or all products + any substrate
"""
function _expand_substrate_product_dead_ends(
    specs::Vector{MechanismSpec},
    @nospecialize(reaction::EnzymeReaction),
)
    sub_names = Set(s[1] for s in substrates(reaction))
    prod_names = Set(p[1] for p in products(reaction))
    all_mets = union(sub_names, prod_names)

    result = MechanismSpec[]
    for spec in specs
        bound = _bound_metabolites_at_forms(
            spec, reaction)
        cat_forms = all_form_names(spec)

        de_opportunities =
            _substrate_product_dead_end_opportunities(
                bound, cat_forms, sub_names,
                prod_names)

        # Deduplicate: multiple catalytic forms may
        # produce the same dead-end form. Group by
        # dead-end form name.
        de_forms = Dict{Symbol,
            Vector{Tuple{Symbol, Symbol}}}()
        for (f, m) in de_opportunities
            de_name = _dead_end_form_name(bound[f], m)
            push!(get!(de_forms, de_name,
                Tuple{Symbol, Symbol}[]), (f, m))
        end
        de_form_names = sort(collect(keys(de_forms)))
        n_de = length(de_form_names)

        # Enumerate 2^n subsets of dead-end forms
        for mask in 0:(1 << n_de) - 1
            active_de = Set{Symbol}()
            for (j, name) in enumerate(de_form_names)
                if (mask >> (j - 1)) & 1 == 1
                    push!(active_de, name)
                end
            end

            # Build new steps: original + dead-end
            new_steps = copy(spec.steps)

            # Add binding steps for active dead-ends
            for de_name in sort(collect(active_de))
                entries = de_forms[de_name]
                for (cat_form, met) in entries
                    # Binding step: [cat_form, met]
                    #   → [de_name] (always RE)
                    push!(new_steps, StepSpec(
                        [cat_form, met],
                        [de_name], true))
                end
            end

            # Add mirror steps: for each catalytic
            # step, if both endpoints have extensions
            # to active dead-end forms with the same
            # metabolite, add a mirror step
            for s in spec.steps
                from, to = s.reactants[1], s.products[1]
                met = step_metabolite(s)
                for de_met in sort(collect(all_mets))
                    haskey(bound, from) || continue
                    haskey(bound, to) || continue
                    de_met in bound[from] && continue
                    de_met in bound[to] && continue
                    from_de = _dead_end_form_name(
                        bound[from], de_met)
                    to_de = _dead_end_form_name(
                        bound[to], de_met)
                    from_de in active_de || continue
                    to_de in active_de || continue
                    # Mirror step inherits RE/SS
                    if met !== nothing
                        push!(new_steps, StepSpec(
                            [from_de, met],
                            [to_de],
                            s.is_equilibrium))
                    else
                        push!(new_steps, StepSpec(
                            [from_de], [to_de],
                            s.is_equilibrium))
                    end
                end
            end

            # Compute param_count
            n_steps = length(new_steps)
            n_re = count(
                s -> s.is_equilibrium, new_steps)
            n_ss = n_steps - n_re
            n_forms = length(all_form_names(new_steps))
            n_thermo = n_steps - n_forms + 1
            pc = n_re + 2 * n_ss - n_thermo + 2

            push!(result, MechanismSpec(
                spec.reaction, new_steps,
                spec.param_constraints, pc))
        end
    end
    result
end

# ─── Stage 3: Combined Dead-End Expansion ─────────────────────

"""
    _substrate_product_dead_end_opportunities(
        bound, cat_forms, sub_names, prod_names)

Find (form, metabolite) dead-end opportunities for
substrates and products. A dead-end is valid when:
- The form doesn't already bind all substrates or all
  products
- The metabolite isn't already bound at the form
- The resulting form isn't a catalytic form
- The result binds at least one substrate AND at least
  one product (mixed binding required)
- The result doesn't have all substrates or all products
"""
function _substrate_product_dead_end_opportunities(
    bound::Dict{Symbol, Set{Symbol}},
    cat_forms::Set{Symbol},
    sub_names::Set{Symbol},
    prod_names::Set{Symbol},
)
    all_mets = union(sub_names, prod_names)
    opportunities = Tuple{Symbol, Symbol}[]
    for f in sort(collect(cat_forms))
        haskey(bound, f) || continue
        fb = bound[f]
        fb_subs = intersect(fb, sub_names)
        fb_prods = intersect(fb, prod_names)
        # Eligible: neither all subs nor all prods
        (fb_subs == sub_names ||
            fb_prods == prod_names) && continue
        for m in sort(collect(all_mets))
            m in fb && continue
            de_name = _dead_end_form_name(fb, m)
            de_name in cat_forms && continue
            new_bound = union(fb, Set([m]))
            new_subs = intersect(
                new_bound, sub_names)
            new_prods = intersect(
                new_bound, prod_names)
            # Must bind at least one of each type
            if isempty(new_subs) || isempty(new_prods)
                continue
            end
            # Must not bind all of either type
            if new_subs == sub_names ||
                    new_prods == prod_names
                continue
            end
            push!(opportunities, (f, m))
        end
    end
    opportunities
end

"""
    _regulator_dead_end_opportunities(
        bound, cat_forms, sub_names, prod_names,
        dead_end_regs)

Find (form, dummy_met) dead-end opportunities for
regulators. Eligible forms are those where neither all
substrates nor all products are bound. Returns pairs of
(catalytic_form, dummy_metabolite_name).
"""
function _regulator_dead_end_opportunities(
    bound::Dict{Symbol, Set{Symbol}},
    cat_forms::Set{Symbol},
    sub_names::Set{Symbol},
    prod_names::Set{Symbol},
    dead_end_regs::Vector{Symbol},
)
    opportunities = Tuple{Symbol, Symbol}[]
    for f in sort(collect(cat_forms))
        haskey(bound, f) || continue
        fb = bound[f]
        fb_subs = intersect(fb, sub_names)
        fb_prods = intersect(fb, prod_names)
        # Eligible: neither all subs nor all prods
        (fb_subs == sub_names ||
            fb_prods == prod_names) && continue
        for (i, reg) in enumerate(
                sort(dead_end_regs))
            dummy = Symbol(
                string(reg) * "__reg" * string(i))
            push!(opportunities, (f, dummy))
        end
    end
    opportunities
end

"""
    _expand_dead_end(specs, reaction;
        dead_end_regs) -> Vector{MechanismSpec}

Combined dead-end expansion: finds substrate/product
AND regulator dead-end opportunities, then enumerates
the power set over all unique dead-end forms in a
single pass. Dead-end binding steps are always RE.
Mirror steps inherit RE/SS from their catalytic
counterparts.
"""
function _expand_dead_end(
    specs::Vector{MechanismSpec},
    @nospecialize(reaction::EnzymeReaction);
    dead_end_regs::Vector{Symbol}=Symbol[],
    include_substrate_product::Bool=true,
)
    sub_names = Set(
        s[1] for s in substrates(reaction))
    prod_names = Set(
        p[1] for p in products(reaction))
    all_mets = union(sub_names, prod_names)

    result = MechanismSpec[]
    for spec in specs
        bound = _bound_metabolites_at_forms(
            spec, reaction)
        cat_forms = all_form_names(spec)

        sp_opps = if include_substrate_product
            _substrate_product_dead_end_opportunities(
                bound, cat_forms, sub_names,
                prod_names)
        else
            Tuple{Symbol, Symbol}[]
        end

        # Collect regulator dead-end opportunities
        reg_opps = _regulator_dead_end_opportunities(
            bound, cat_forms, sub_names, prod_names,
            dead_end_regs)

        # Merge all opportunities
        all_opps = vcat(sp_opps, reg_opps)

        # Group by dead-end form name. For
        # substrate/product mets, the dead-end form
        # name uses the bound set. For regulators, we
        # extend the bound set with the dummy name.
        de_forms = Dict{Symbol,
            Vector{Tuple{Symbol, Symbol}}}()
        for (f, m) in all_opps
            de_name = _dead_end_form_name(
                bound[f], m)
            push!(get!(de_forms, de_name,
                Tuple{Symbol, Symbol}[]), (f, m))
        end
        de_form_names = sort(collect(
            keys(de_forms)))
        n_de = length(de_form_names)

        # Enumerate 2^n subsets of dead-end forms
        for mask in 0:(1 << n_de) - 1
            active_de = Set{Symbol}()
            for (j, name) in enumerate(
                    de_form_names)
                if (mask >> (j - 1)) & 1 == 1
                    push!(active_de, name)
                end
            end

            new_steps = copy(spec.steps)

            # Add binding steps for active dead-ends
            for de_name in sort(collect(active_de))
                entries = de_forms[de_name]
                for (cat_form, met) in entries
                    push!(new_steps, StepSpec(
                        [cat_form, met],
                        [de_name], true))
                end
            end

            # Add mirror steps: for each catalytic
            # step, check if both endpoints can be
            # extended with the same dead-end met
            # All de mets = union of sub/prod mets
            # and regulator dummy names
            all_de_mets = Set{Symbol}()
            for (_, m) in all_opps
                push!(all_de_mets, m)
            end
            for s in spec.steps
                from = s.reactants[1]
                to = s.products[1]
                met = step_metabolite(s)
                for de_met in sort(
                        collect(all_de_mets))
                    haskey(bound, from) || continue
                    haskey(bound, to) || continue
                    de_met in bound[from] && continue
                    de_met in bound[to] && continue
                    from_de = _dead_end_form_name(
                        bound[from], de_met)
                    to_de = _dead_end_form_name(
                        bound[to], de_met)
                    from_de in active_de || continue
                    to_de in active_de || continue
                    if met !== nothing
                        push!(new_steps, StepSpec(
                            [from_de, met],
                            [to_de],
                            s.is_equilibrium))
                    else
                        push!(new_steps, StepSpec(
                            [from_de], [to_de],
                            s.is_equilibrium))
                    end
                end
            end

            # Compute param_count
            n_steps = length(new_steps)
            n_re = count(
                s -> s.is_equilibrium, new_steps)
            n_ss = n_steps - n_re
            n_forms = length(
                all_form_names(new_steps))
            n_thermo = n_steps - n_forms + 1
            pc = n_re + 2 * n_ss - n_thermo + 2

            push!(result, MechanismSpec(
                spec.reaction, new_steps,
                spec.param_constraints, pc))
        end
    end
    result
end

"""
Count equivalence constraints that make Wegscheider
constraints redundant. Works in equilibrium-constant
space (one dimension per step): cycle rows from the
graph null space, equiv rows as p_a - p_b = 0.
Returns the number of redundant thermo constraints.
"""
function _redundant_thermo_count(
    steps::Vector{StepSpec},
    active_groups::Vector{Vector{Int}},
)
    # Build incidence matrix (forms × steps)
    form_list = Symbol[]
    form_idx = Dict{Symbol, Int}()
    for s in steps
        for f in (s.reactants[1], s.products[1])
            if !haskey(form_idx, f)
                push!(form_list, f)
                form_idx[f] = length(form_list)
            end
        end
    end
    n_forms = length(form_list)
    n_steps = length(steps)
    B = zeros(Int, n_forms, n_steps)
    for (j, s) in enumerate(steps)
        B[form_idx[s.reactants[1]], j] -= 1
        B[form_idx[s.products[1]], j] += 1
    end

    NS = _integer_nullspace(B)
    n_cycles = size(NS, 2)
    n_cycles == 0 && return 0

    # Build equiv rows in step space
    equiv_rows = Int[]
    for g in active_groups
        for j in 2:length(g)
            push!(equiv_rows, g[1], g[j])
        end
    end
    n_equiv = length(equiv_rows) ÷ 2
    n_equiv == 0 && return 0

    # Combined matrix: cycles (transposed) + equiv
    n_total = n_cycles + n_equiv
    M = zeros(Rational{BigInt}, n_total, n_steps)
    for i in 1:n_cycles
        for j in 1:n_steps
            M[i, j] = NS[j, i]
        end
    end
    for k in 1:n_equiv
        a = equiv_rows[2k - 1]
        b = equiv_rows[2k]
        M[n_cycles + k, a] = 1
        M[n_cycles + k, b] = -1
    end

    # Rank via row echelon
    rank_combined = 0
    row = 1
    Mc = copy(M)
    for col in 1:n_steps
        piv = findfirst(
            r -> Mc[r, col] != 0, row:n_total)
        piv === nothing && continue
        piv += row - 1
        Mc[row, :], Mc[piv, :] =
            Mc[piv, :], Mc[row, :]
        Mc[row, :] ./= Mc[row, col]
        for r in 1:n_total
            r != row && Mc[r, col] != 0 &&
                (Mc[r, :] .-= Mc[r, col] .* Mc[row, :])
        end
        rank_combined += 1
        row += 1
    end

    n_cycles + n_equiv - rank_combined
end

# ─── Stage 4: Equivalence Constraint Expansion ────────────────

"""
    _expand_equivalence_constraints(specs, reaction)
        -> Vector{MechanismSpec}

For each spec, enumerate equivalence constraint masks.
Groups steps by (metabolite_name, is_equilibrium). Steps
binding the same metabolite with the same RE/SS status can
share parameters. Dummy regulatory names (:X__reg1)
naturally separate catalytic from regulatory binding.
"""
function _expand_equivalence_constraints(
    specs::Vector{MechanismSpec},
    @nospecialize(reaction::EnzymeReaction),
)
    result = MechanismSpec[]
    for spec in specs
        # Group step indices by (metabolite, RE/SS)
        groups = Dict{
            Tuple{Symbol,Bool}, Vector{Int}}()
        for (i, s) in enumerate(spec.steps)
            met = step_metabolite(s)
            met === nothing && continue
            key = (met, s.is_equilibrium)
            push!(get!(groups, key, Int[]), i)
        end

        # Keep only groups with ≥ 2 steps
        valid_groups = sort!(
            [sort!(g) for (_, g) in groups
             if length(g) >= 2];
            by=first)

        n_groups = length(valid_groups)
        for mask in 0:(1 << n_groups) - 1
            constraints = ParamConstraint[]
            delta = 0
            for (gi, g) in enumerate(valid_groups)
                (mask >> (gi - 1)) & 1 == 1 ||
                    continue
                is_re = spec.steps[g[1]].is_equilibrium
                n_constrained = length(g) - 1
                if is_re
                    delta -= n_constrained
                    for j in 2:length(g)
                        push!(constraints, (
                            Symbol("K$(g[j])"),
                            1,
                            [(Symbol("K$(g[1])"), 1)]
                        ))
                    end
                else
                    delta -= 2 * n_constrained
                    for j in 2:length(g)
                        for sfx in ("f", "r")
                            push!(constraints, (
                                Symbol(
                                    "k$(g[j])$sfx"),
                                1,
                                [(Symbol(
                                    "k$(g[1])$sfx"),
                                  1)]
                            ))
                        end
                    end
                end
            end
            # Equivalence constraints can make
            # Wegscheider constraints redundant
            active = Vector{Int}[]
            for (gi, g) in enumerate(valid_groups)
                (mask >> (gi - 1)) & 1 == 1 &&
                    push!(active, g)
            end
            redundancy = _redundant_thermo_count(
                spec.steps, active)
            push!(result, MechanismSpec(
                spec.reaction, spec.steps,
                constraints,
                spec.param_count + delta +
                    redundancy))
        end
    end
    result
end

# ─── Stage 5: Deduplication ────────────────────────────────────

"""
    _deduplicate(specs, reaction) -> Vector{MechanismSpec}

Deduplicate by (concentration fingerprint, constraint descriptor).
Keeps mechanism with fewest parameters.
"""
function _deduplicate(
    specs::Vector{MechanismSpec},
    @nospecialize(reaction::EnzymeReaction),
)
    isempty(specs) && return specs

    best = Dict{_DedupKey, MechanismSpec}()
    for spec in specs
        steps = spec.steps
        partition = _compute_re_partition_from_steps(steps)
        fp = _concentration_fingerprint(steps, partition)

        # Build valid equivalence groups (same as Stage 4)
        groups = Dict{
            Tuple{Symbol,Bool}, Vector{Int}}()
        for (i, s) in enumerate(steps)
            met = step_metabolite(s)
            met === nothing && continue
            key = (met, s.is_equilibrium)
            push!(get!(groups, key, Int[]), i)
        end
        valid_groups = sort!(
            [sort!(g) for (_, g) in groups
             if length(g) >= 2];
            by=first)

        constraint_mask = _constraints_to_mask(
            spec.param_constraints, valid_groups, steps)
        desc = _constraint_descriptor(
            steps, valid_groups, constraint_mask)

        dedup_key = (fp, desc)
        if !haskey(best, dedup_key) ||
                spec.param_count <
                    best[dedup_key].param_count
            best[dedup_key] = spec
        end
    end
    collect(values(best))
end

"""Reverse-map param_constraints to a bitmask over valid_groups."""
function _constraints_to_mask(constraints, valid_groups,
                              steps::Vector{StepSpec})
    mask = 0
    constrained_step_indices = Set{Int}()
    for (target, _, srcs) in constraints
        m = match(r"[kK](\d+)", string(target))
        m !== nothing && m[1] !== nothing &&
            push!(constrained_step_indices,
                parse(Int, m[1]::SubString))
    end
    for (gi, g) in enumerate(valid_groups)
        if any(idx in constrained_step_indices
               for idx in g[2:end])
            mask |= (1 << (gi - 1))
        end
    end
    mask
end

# ─── Allosteric Expansion ─────────────────────────────────────

"""
    _expand_allosteric(specs, reaction;
        catalytic_n, allosteric_regs)
        -> Vector{AllostericMechanismSpec}

Expand monomeric specs into allosteric (MWC) variants.
"""
function _expand_allosteric(
    specs::Vector{MechanismSpec},
    @nospecialize(reaction::EnzymeReaction);
    catalytic_n::Int=2,
    allosteric_regs::Vector{Symbol}=Symbol[],
)
    result = AllostericMechanismSpec[]
    if isempty(allosteric_regs)
        # No allosteric regulators: catalytic subunits only
        for spec in specs
            push!(result, AllostericMechanismSpec(
                spec, catalytic_n,
                Vector{Symbol}[], Int[], Symbol[]))
        end
        return result
    end
    partitions = _set_partitions(allosteric_regs)
    for spec in specs
        for partition in partitions
            n_groups = length(partition)
            for combo in Iterators.product(
                    ntuple(_ -> 1:catalytic_n, n_groups)...)
                push!(result, AllostericMechanismSpec(
                    spec, catalytic_n, partition,
                    collect(combo), Symbol[]))
            end
        end
    end
    result
end

# ─── T/R Equivalence ─────────────────────────────────────────

"""All metabolites with T-state binding parameters."""
function _collect_t_state_metabolites(
    spec::AllostericMechanismSpec,
)
    t_mets = Symbol[]
    for s in spec.base.steps
        s.is_equilibrium || continue
        met = step_metabolite(s)
        met !== nothing && met ∉ t_mets &&
            push!(t_mets, met)
    end
    for site in spec.allosteric_reg_sites
        for lig in site
            lig ∉ t_mets && push!(t_mets, lig)
        end
    end
    t_mets
end

"""
    _expand_tr_equivalence(specs, reaction)
        -> Vector{AllostericMechanismSpec}

Enumerate T/R parameter equivalence variants. For each metabolite
with a T-state parameter, K_T can equal K_R (fewer params) or be
independent. Produces 2^n variants per input spec where n is the
number of metabolites with T-state binding parameters.
"""
function _expand_tr_equivalence(
    specs::Vector{AllostericMechanismSpec},
    @nospecialize(reaction::EnzymeReaction),
)
    result = AllostericMechanismSpec[]
    for spec in specs
        t_mets = _collect_t_state_metabolites(spec)
        n = length(t_mets)
        for mask in 0:(1 << n) - 1
            equiv = Symbol[t_mets[i] for i in 1:n
                          if ((mask >> (i-1)) & 1) == 1]
            push!(result, AllostericMechanismSpec(
                spec.base, spec.catalytic_n,
                spec.allosteric_reg_sites,
                spec.allosteric_multiplicities,
                equiv))
        end
    end
    result
end

# ─── Post-Allosteric Deduplication ────────────────────────────

"""
    _deduplicate_allosteric(specs, reaction) -> Vector{AllostericMechanismSpec}

Remove T<->R mirror duplicates.
"""
function _deduplicate_allosteric(
    specs::Vector{AllostericMechanismSpec},
    @nospecialize(reaction::EnzymeReaction),
)
    seen = Dict{Any, AllostericMechanismSpec}()
    for spec in specs
        key = _allosteric_canonical_key(spec)
        if !haskey(seen, key) ||
                spec.base.param_count <
                seen[key].base.param_count
            seen[key] = spec
        end
    end
    collect(values(seen))
end

"""Canonical key for allosteric dedup.

Maps complementary TR-equiv sets to the same key so
T/R mirror mechanisms are recognized as duplicates.
"""
function _allosteric_canonical_key(spec::AllostericMechanismSpec)
    base_key = (spec.base.steps,
                spec.base.param_constraints)
    # Sort reg sites and multiplicities together
    pairs = collect(zip(spec.allosteric_reg_sites,
                        spec.allosteric_multiplicities))
    sort!(pairs)
    sorted_sites = [p[1] for p in pairs]
    sorted_mults = [p[2] for p in pairs]
    # Canonical TR-equiv: min of set and its complement
    # so T↔R mirrors map to the same key
    tr = sort(spec.tr_equiv_metabolites)
    all_t = sort(_collect_t_state_metabolites(spec))
    complement = sort(setdiff(all_t, tr))
    canonical_tr = min(tr, complement)
    (base_key, spec.catalytic_n, sorted_sites,
     sorted_mults, canonical_tr)
end

# ─── MechanismSpec → EnzymeMechanism Conversion ──────────────

"""
    compile_mechanism(spec::MechanismSpec)

Convert a `MechanismSpec` to an `EnzymeMechanism`.
"""
function compile_mechanism(spec::MechanismSpec)
    _compile_enzyme_mechanism(spec)
end

"""Construct EnzymeMechanism from MechanismSpec."""
EnzymeMechanism(spec::MechanismSpec) =
    _compile_enzyme_mechanism(spec)

function _compile_enzyme_mechanism(spec::MechanismSpec)
    rxn = spec.reaction

    # Collect metabolite names and atom dicts
    met_atoms = Dict{Symbol,Dict{Symbol,Int}}()
    for (name, atoms) in substrates(rxn)
        met_atoms[name] = Dict{Symbol,Int}(
            a => c for (a, c) in atoms
        )
    end
    for (name, atoms) in products(rxn)
        met_atoms[name] = Dict{Symbol,Int}(
            a => c for (a, c) in atoms
        )
    end
    met_set = Set(keys(met_atoms))
    for r in regulators(rxn)
        push!(met_set, r)
        met_atoms[r] = Dict{Symbol,Int}()
    end

    # Strip __regN suffixes from metabolite names
    function _clean_met(sym::Symbol)
        s = string(sym)
        m = match(r"^(.+)__reg\d+$", s)
        m !== nothing ? Symbol(m.captures[1]) : sym
    end

    # Add suffixed regulator names to met_set/met_atoms
    # so BFS recognizes them as metabolites
    for s in spec.steps
        for sym in Iterators.flatten(
                (s.reactants, s.products))
            clean = _clean_met(sym)
            if clean != sym && clean ∈ met_set
                push!(met_set, sym)
                met_atoms[sym] = met_atoms[clean]
            end
        end
    end

    # Collect form names (all symbols in steps
    # that are not metabolites)
    form_set = Set{Symbol}()
    for s in spec.steps
        for sym in s.reactants
            sym ∉ met_set && push!(form_set, sym)
        end
        for sym in s.products
            sym ∉ met_set && push!(form_set, sym)
        end
    end

    # BFS from :E to compute atoms for each form
    form_atoms = Dict{Symbol,Dict{Symbol,Int}}(
        :E => Dict{Symbol,Int}()
    )
    # Build adjacency: form -> [(neighbor, met_or_nothing,
    #   direction)] where direction = :add if metabolite
    #   binds going from form to neighbor
    adj = Dict{Symbol,
        Vector{Tuple{Symbol,Union{Nothing,Symbol},
                      Symbol}}}()
    for s in spec.steps
        from = s.reactants[1]
        to = s.products[1]
        met = length(s.reactants) == 2 ?
            s.reactants[2] : nothing
        # Canonical: met on LHS means binding
        # from→to. Reverse direction = release.
        if !haskey(adj, from)
            adj[from] = Tuple{
                Symbol,Union{Nothing,Symbol},Symbol
            }[]
        end
        if !haskey(adj, to)
            adj[to] = Tuple{
                Symbol,Union{Nothing,Symbol},Symbol
            }[]
        end
        push!(adj[from], (to, met, :add))
        push!(adj[to], (from, met, :subtract))
    end

    queue = Symbol[:E]
    while !isempty(queue)
        cur = popfirst!(queue)
        cur_atoms = form_atoms[cur]
        haskey(adj, cur) || continue
        for (nbr, met, dir) in adj[cur]
            haskey(form_atoms, nbr) && continue
            nbr_atoms = copy(cur_atoms)
            if met !== nothing
                ma = met_atoms[met]
                if dir == :add
                    for (a, c) in ma
                        nbr_atoms[a] =
                            get(nbr_atoms, a, 0) + c
                    end
                else  # :subtract
                    for (a, c) in ma
                        nbr_atoms[a] =
                            get(nbr_atoms, a, 0) - c
                        nbr_atoms[a] == 0 &&
                            delete!(nbr_atoms, a)
                    end
                end
            end
            # Isomerization: no atom change
            form_atoms[nbr] = nbr_atoms
            push!(queue, nbr)
        end
    end

    # Build enzyme forms tuple sorted by cleaned name
    form_names = sort!(collect(form_set))
    enzymes = Tuple(
        (_clean_met(name), Tuple(sort!(
            [Tuple(p) for p in form_atoms[name]];
            by=first
        )))
        for name in form_names
    )

    species = (
        substrates(rxn), products(rxn),
        regulators(rxn), enzymes
    )

    reactions = Tuple(
        let r = s.reactants, p = s.products
            lhs = Tuple(_clean_met(x) for x in r)
            rhs = Tuple(_clean_met(x) for x in p)
            (lhs, rhs)
        end
        for s in spec.steps
    )

    eq_steps = Tuple(s.is_equilibrium for s in spec.steps)

    constraints = Tuple(
        (t, c, Tuple(Tuple.(f)))
        for (t, c, f) in spec.param_constraints
    )

    EnzymeMechanism(species, reactions, eq_steps,
        constraints)
end

function compile_mechanism(spec::AllostericMechanismSpec)
    cm = compile_mechanism(spec.base)
    cat_mets = metabolites(cm)

    # Build Metabolites tuple (catalytic + regulatory)
    reg_syms = Symbol[]
    for site in spec.allosteric_reg_sites
        for s in site
            s in reg_syms || s in cat_mets || push!(reg_syms, s)
        end
    end
    mets = (cat_mets..., reg_syms...)

    # Build CatSites: (catalytic_metabolites, multiplicity, tr_equiv_mets)
    cat_tr = Tuple(m for m in cat_mets
                   if m in spec.tr_equiv_metabolites)
    cat_sites = (cat_mets, spec.catalytic_n, cat_tr)

    # Build RegSites with TR equivalence info
    reg_sites = Tuple(
        (Tuple(group), mult,
         Tuple(lig for lig in group
               if lig in spec.tr_equiv_metabolites))
        for (group, mult) in zip(
            spec.allosteric_reg_sites,
            spec.allosteric_multiplicities))

    AllostericEnzymeMechanism{mets, typeof(cm), cat_sites, reg_sites}()
end

# ─── Pipeline Orchestration ──────────────────────────────────

"""
    enumerate_mechanisms(reaction; max_re_groups=7, catalytic_n=0)

Enumerate valid mechanism topologies for the given reaction
using a staged pipeline.
"""
function enumerate_mechanisms(
    @nospecialize(reaction::EnzymeReaction);
    max_re_groups::Int=7,
    catalytic_n::Int=0,
)
    catalytic = _catalytic_topologies(reaction)

    roles = regulator_roles(reaction)
    fixed_dead_end = Symbol[
        r[1] for r in roles if r[2] == :dead_end]
    fixed_allosteric = Symbol[
        r[1] for r in roles if r[2] == :allosteric]
    unknown = Symbol[
        r[1] for r in roles if r[2] == :unknown]
    n_unknown = length(unknown)

    all_base = MechanismSpec[]
    all_allosteric = AllostericMechanismSpec[]

    for reg_mask in 0:(1 << n_unknown) - 1
        de_regs = Symbol[fixed_dead_end;
            [unknown[i] for i in 1:n_unknown
             if (reg_mask >> (i - 1)) & 1 == 0]]
        allo_regs = Symbol[fixed_allosteric;
            [unknown[i] for i in 1:n_unknown
             if (reg_mask >> (i - 1)) & 1 == 1]]

        # Stages 2-5: base mechanism pipeline
        base = _expand_ress_variants(
            catalytic, reaction; max_re_groups)
        base = _expand_dead_end(
            base, reaction; dead_end_regs=de_regs,
            include_substrate_product=true)
        base = _expand_equivalence_constraints(
            base, reaction)
        base = _deduplicate(base, reaction)
        append!(all_base, base)

        # Stages 6-8: allosteric expansion
        if !isempty(allo_regs)
            cn = catalytic_n > 0 ? catalytic_n : 1
            allo = _expand_allosteric(
                base, reaction;
                catalytic_n=cn,
                allosteric_regs=allo_regs)
            allo = _expand_tr_equivalence(
                allo, reaction)
            allo = _deduplicate_allosteric(
                allo, reaction)
            append!(all_allosteric, allo)
        end
    end

    total = length(all_base) + length(all_allosteric)
    inner = Iterators.flatten(
        (all_base, all_allosteric))
    MechanismIterator(inner, total)
end
