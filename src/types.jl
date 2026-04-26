using LinearAlgebra: rank

"""Sort species tuples alphabetically by name (first element)."""
_sort_species(t::Tuple) = Tuple(sort(collect(t); by=s -> s[1]))

"""
    EnzymeReaction{Substrates, Products, Regulators, OligomericState}

Singleton type encoding an enzyme reaction specification in type parameters.

- `Substrates`, `Products`: tuple of `(name::Symbol, atoms::Tuple{Vararg{Tuple{Symbol,Int}}})`.
- `Regulators`: tuple of `Symbol` (plain names, no atoms).
- `OligomericState`: number of subunits (Int, default 1).
"""
struct EnzymeReaction{Substrates, Products, Regulators, OligomericState} end

function EnzymeReaction(subs::Tuple, prods::Tuple, regs::Tuple=(); oligomeric_state::Int=1)
    isempty(subs) && error("Substrates must not be empty")
    isempty(prods) && error("Products must not be empty")
    subs_names = [s[1] for s in subs]
    prods_names = [s[1] for s in prods]
    length(subs_names) != length(Set(subs_names)) &&
        error("Duplicate substrate names")
    length(prods_names) != length(Set(prods_names)) &&
        error("Duplicate product names")
    # Normalize regulators to (name, role) pairs
    normalized_regs = if isempty(regs)
        regs
    elseif regs[1] isa Symbol
        Tuple((r, :unknown) for r in regs)
    else
        regs
    end
    for r in normalized_regs
        r isa Tuple{Symbol,Symbol} ||
            error("Regulators must be (Symbol, Symbol) pairs, got $r")
    end
    reg_names = [r[1] for r in normalized_regs]
    length(reg_names) != length(Set(reg_names)) &&
        error("Duplicate regulator names")
    subs = _sort_species(subs)
    prods = _sort_species(prods)
    sorted_regs = Tuple(sort(collect(normalized_regs); by=first))
    EnzymeReaction{subs, prods, sorted_regs, oligomeric_state}()
end

abstract type AbstractEnzymeMechanism end

"""
    EnzymeMechanism{Metabolites, Reactions}

Singleton type encoding an enzyme mechanism.

- `Metabolites`: 3-tuple `(substrates::Tuple{Symbol,...}, products::Tuple{Symbol,...},
  regulators::Tuple{Symbol,...})`. Plain symbol names — no atom content stored.
- `Reactions`: tuple of 4-tuples `(lhs_syms, rhs_syms, is_eq::Bool,
  kinetic_group::Int)`. Steps with identical `kinetic_group` share kinetic
  parameters (one `K` for RE groups, one `k_f` and one `k_r` for SS groups).
"""
struct EnzymeMechanism{Metabolites, Reactions} <: AbstractEnzymeMechanism end

"""Count enzymes, metabolites, atoms, and metabolite names on one side of a reaction."""
function _count_side(side, enzyme_set, enzyme_atoms, met_atoms, step_idx)
    n_enz, n_met, atoms, mets = 0, 0, Dict{Symbol,Int}(), Symbol[]
    for s in side
        if s in enzyme_set
            n_enz += 1
            for (a, c) in enzyme_atoms[s]
                atoms[a] = get(atoms, a, 0) + c
            end
        elseif haskey(met_atoms, s)
            n_met += 1
            push!(mets, s)
            for (a, c) in met_atoms[s]
                atoms[a] = get(atoms, a, 0) + c
            end
        else
            error("Reaction $(step_idx) uses unknown species $(s)")
        end
    end
    (n_enz, n_met, atoms, mets)
end

"""
    EnzymeMechanism(metabolites, reactions) → EnzymeMechanism

Construct an `EnzymeMechanism` from explicit metabolite 3-tuple and
reaction 4-tuples. Step order is preserved (no canonicalization).
Validates structure, enzyme-form connectivity, kinetic-group
composition rules, and stoichiometric feasibility.
"""
function EnzymeMechanism(
    mets::Tuple{Tuple{Vararg{Symbol}}, Tuple{Vararg{Symbol}}, Tuple{Vararg{Symbol}}},
    rxns::Tuple,
)
    subs, prods, regs = mets

    # Sort each species list alphabetically (canonical Metabolites parameter).
    # Step order is preserved per spec §4.1.2.
    subs = Tuple(sort(collect(subs)))
    prods = Tuple(sort(collect(prods)))
    regs = Tuple(sort(collect(regs)))
    mets = (subs, prods, regs)

    isempty(rxns) && error("Reactions tuple must not be empty")
    all(step[3] for step in rxns) &&
        error("At least one SS step required (not all steps can be RE)")

    # Build metabolite set
    met_set = Set{Symbol}()
    for group in (subs, prods, regs)
        for name in group; push!(met_set, name); end
    end

    # Validate each step shape
    for (i, step) in enumerate(rxns)
        length(step) == 4 ||
            error("Step $i must be (lhs, rhs, is_eq, kinetic_group); got $step")
        lhs, rhs, is_eq, gnum = step
        is_eq isa Bool || error("Step $i is_eq must be Bool")
        gnum isa Int || error("Step $i kinetic_group must be Int")
        n_enz_lhs = count(s -> s ∉ met_set, lhs)
        n_enz_rhs = count(s -> s ∉ met_set, rhs)
        n_enz_lhs == 1 ||
            error("Step $i LHS must contain exactly one enzyme form; got $lhs")
        n_enz_rhs == 1 ||
            error("Step $i RHS must contain exactly one enzyme form; got $rhs")
        n_met_lhs = count(s -> s ∈ met_set, lhs)
        n_met_rhs = count(s -> s ∈ met_set, rhs)
        n_met_lhs <= 1 || error("Step $i LHS has more than one metabolite")
        n_met_rhs <= 1 || error("Step $i RHS has more than one metabolite")
    end

    # Canonicalize RE step direction (metabolite on LHS for binding steps).
    rxns = ntuple(length(rxns)) do i
        (lhs, rhs, is_eq, gnum) = rxns[i]
        if !is_eq
            return (lhs, rhs, is_eq, gnum)
        end
        rhs_has_met = any(s in met_set for s in rhs)
        lhs_has_met = any(s in met_set for s in lhs)
        if rhs_has_met && !lhs_has_met
            (rhs, lhs, is_eq, gnum)
        else
            (lhs, rhs, is_eq, gnum)
        end
    end

    # Each substrate / product / regulator must appear in some step
    appears = Set{Symbol}()
    for (lhs, rhs, _, _) in rxns
        for s in lhs; push!(appears, s); end
        for s in rhs; push!(appears, s); end
    end
    for name in vcat(collect(subs), collect(prods), collect(regs))
        name in appears ||
            error("Listed metabolite $name does not appear in any reaction step")
    end

    # Kinetic-group composition rules
    _validate_kinetic_groups(rxns, met_set)

    # Build the singleton type and run remaining checks via accessors
    m = EnzymeMechanism{mets, rxns}()

    # Enzyme-form graph weakly connected
    _validate_enzyme_connectivity(m)

    # Stoichiometric feasibility (rank check)
    _validate_stoichiometry(m)

    m
end

"""
Validate kinetic-group composition: 2+ groups must be all RE binding
or all SS binding (same metabolite); iso steps must be singletons.
"""
function _validate_kinetic_groups(rxns, met_set)
    groups = Dict{Int, Vector{Int}}()
    for (i, step) in enumerate(rxns)
        push!(get!(groups, step[4], Int[]), i)
    end
    for (g, idxs) in groups
        length(idxs) == 1 && continue
        kinds = map(idxs) do i
            lhs, rhs, is_eq, _ = rxns[i]
            mets_in = [s for s in lhs if s in met_set]
            mets_out = [s for s in rhs if s in met_set]
            isempty(mets_in) && isempty(mets_out) &&
                error("Iso step (no metabolite) at index $i must be a " *
                      "singleton kinetic group; found in group $g of size $(length(idxs))")
            length(mets_in) == 1 ||
                error("Step $i has $(length(mets_in)) metabolites on LHS; expected 1")
            (is_eq, mets_in[1])
        end
        first_kind = kinds[1]
        for (i, k) in zip(idxs[2:end], kinds[2:end])
            k[1] == first_kind[1] ||
                error("Kinetic group $g contains both RE and SS binding steps")
            k[2] == first_kind[2] ||
                error("Kinetic group $g binds different metabolites: " *
                      "$(first_kind[2]) and $(k[2])")
        end
    end
end

"""Verify the enzyme-form graph is weakly connected."""
function _validate_enzyme_connectivity(m::EnzymeMechanism)
    enz = enzyme_forms(m)
    isempty(enz) && error("Mechanism has no enzyme forms")
    name_set = Set(enz)
    adj = Dict(n => Set{Symbol}() for n in enz)
    for (lhs, rhs, _, _) in reactions(m)
        e_l = first(s for s in lhs if s in name_set)
        e_r = first(s for s in rhs if s in name_set)
        push!(adj[e_l], e_r)
        push!(adj[e_r], e_l)
    end
    visited = Set{Symbol}()
    queue = [first(enz)]
    while !isempty(queue)
        cur = popfirst!(queue)
        cur in visited && continue
        push!(visited, cur)
        for n in adj[cur]; n in visited || push!(queue, n); end
    end
    visited == name_set ||
        error("Enzyme-form graph not connected; orphan forms: " *
              "$(setdiff(name_set, visited))")
end

"""
Stoichiometric feasibility via `r ∈ col(S)` rank test on the full
stoichiometry matrix. The target vector r has 0 on enzyme rows and
on regulator rows, and ±count(M in subs/prods) on substrate/product rows.
"""
function _validate_stoichiometry(m::EnzymeMechanism)
    S = stoich_matrix(m)
    species = (enzyme_forms(m)..., metabolites(m)...)
    sp_idx = Dict(s => i for (i, s) in enumerate(species))
    r = zeros(Int, length(species))
    for s in substrates(m); r[sp_idx[s]] -= 1; end
    for p in products(m);   r[sp_idx[p]] += 1; end

    rs = Rational.(S)
    rr = Rational.(r)
    rank(rs) == rank(hcat(rs, rr)) ||
        error("Mechanism stoichiometry does not match the declared net reaction. " *
              "Check substrate / product multiplicities and that regulators " *
              "have net zero change. Declared: " *
              "$(_pretty_reaction(substrates(m), products(m)))")
end

_pretty_reaction(subs, prods) =
    "$(join(string.(subs), " + ")) → $(join(string.(prods), " + "))"

"""
    AllostericEnzymeMechanism{Metabolites, CatalyticMech, CatSites, RegSites}

Singleton type for allosteric enzymes (MWC model, always 2 conformations).

- `Metabolites`: tuple of `Symbol` names from `metabolites:` block
- `CatalyticMech`: `EnzymeMechanism` type for one catalytic subunit
- `CatSites`: `(catalytic_metabolites, multiplicity, tr_equiv_mets,
  tr_equiv_cat_steps, r_only_mets, t_only_mets, r_only_cat_steps)`
- `RegSites`: tuple of `((ligand_syms...,), multiplicity, tr_equiv_ligands,
  r_only_ligands, t_only_ligands)` quintuples
"""
struct AllostericEnzymeMechanism{
    Metabolites, CatalyticMech, CatSites, RegSites,
} <: AbstractEnzymeMechanism end

# --- Rate equation mode types ---

"""
    AbstractRateEquationMode

Abstract supertype for rate equation parameterization modes.
Controls which form of the rate equation is used and what parameters are expected.
"""
abstract type AbstractRateEquationMode end

"""
    FullMode <: AbstractRateEquationMode

Full rate equation mode using all 2N microscopic rate constants (k1f, k1r, k2f, k2r, ...).
No thermodynamic constraints applied. Parameters: all k's + E_total.
"""
struct FullMode <: AbstractRateEquationMode end

"""
    ReducedMode <: AbstractRateEquationMode

Rate equation with Haldane-Wegscheider thermodynamic constraints applied.
Dependent parameters are substituted in terms of independent k's and Keq.
Parameters: independent k's + Keq + E_total.
"""
struct ReducedMode <: AbstractRateEquationMode end

"""Singleton instance for full (raw) mode."""
const Full = FullMode()

"""Singleton instance for reduced (Haldane-Wegscheider) mode."""
const Reduced = ReducedMode()

# --- Pretty printing ---

function Base.show(io::IO, ::EnzymeReaction{S,P,R,N}) where {S,P,R,N}
    subs_str = join([string(name) for (name, _) in S], " + ")
    prods_str = join([string(name) for (name, _) in P], " + ")
    print(io, "EnzymeReaction: ", subs_str, " ⇌ ", prods_str)
    if !isempty(R)
        regs_str = join([string(r[1]) for r in R], ", ")
        print(io, " | regulators: ", regs_str)
    end
    N > 1 && print(io, " | oligomeric_state: ", N)
end

# ─── Accessors ─────────────────────────────────────────────────

"""Return substrates as a tuple of `Symbol` names."""
substrates(::EnzymeMechanism{M}) where {M} = M[1]
substrates(::EnzymeReaction{S,P,R,N}) where {S,P,R,N} = S

"""Return products as a tuple of `Symbol` names."""
products(::EnzymeMechanism{M}) where {M} = M[2]
products(::EnzymeReaction{S,P,R,N}) where {S,P,R,N} = P

"""Return regulators as a tuple of `Symbol` names."""
regulators(::EnzymeMechanism{M}) where {M} = M[3]
regulators(::EnzymeReaction{S,P,R,N}) where {S,P,R,N} =
    Tuple(r[1] for r in R)

"""Return regulator (name, role) pairs."""
regulator_roles(::EnzymeReaction{S,P,R,N}) where {S,P,R,N} = R

"""Return oligomeric state (number of subunits)."""
oligomeric_state(::EnzymeReaction{S,P,R,N}) where {S,P,R,N} = N

"""
    metabolites(m::EnzymeMechanism) → Tuple{Symbol,...}

Return distinct metabolite names (substrates ∪ products ∪ regulators) as a tuple
of `Symbol`s in declaration order, deduplicated.
"""
@generated function metabolites(::EnzymeMechanism{M}) where {M}
    seen = Set{Symbol}()
    names = Symbol[]
    for group in M
        for name in group
            if name ∉ seen
                push!(seen, name)
                push!(names, name)
            end
        end
    end
    Tuple(names)
end

"""Return the reactions tuple `((lhs, rhs, is_eq, kinetic_group), ...)`."""
reactions(::EnzymeMechanism{M, R}) where {M, R} = R

"""Return the equilibrium-step flags (`true` = rapid-equilibrium, `false` = steady-state)."""
@generated function equilibrium_steps(::EnzymeMechanism{M, R}) where {M, R}
    Tuple(step[3] for step in R)
end

"""Number of steps in the mechanism."""
n_steps(::EnzymeMechanism{M, R}) where {M, R} = length(R)

"""Kinetic group of step `idx`."""
kinetic_group(::EnzymeMechanism{M, R}, idx::Int) where {M, R} = R[idx][4]

"""Sorted tuple of distinct kinetic group ids."""
@generated function kinetic_groups(::EnzymeMechanism{M, R}) where {M, R}
    Tuple(sort(unique(step[4] for step in R)))
end

"""Indices of steps belonging to kinetic group `G`."""
@generated function steps_in_group(
    ::EnzymeMechanism{M, R}, ::Val{G},
) where {M, R, G}
    Tuple(i for (i, step) in enumerate(R) if step[4] == G)
end
steps_in_group(m::EnzymeMechanism, g::Int) = steps_in_group(m, Val(g))

"""
    enzyme_forms(m::EnzymeMechanism) → Tuple{Symbol,...}

Return distinct enzyme-form names (any symbol appearing in a step that is not a
metabolite) as a tuple of `Symbol`s in step-order, deduplicated.
"""
@generated function enzyme_forms(::EnzymeMechanism{M, R}) where {M, R}
    met_names = Set{Symbol}()
    for group in M; for name in group; push!(met_names, name); end; end
    seen = Set{Symbol}()
    forms = Symbol[]
    for (lhs, rhs, _, _) in R
        for s in lhs; s ∉ met_names && s ∉ seen && (push!(seen, s); push!(forms, s)); end
        for s in rhs; s ∉ met_names && s ∉ seen && (push!(seen, s); push!(forms, s)); end
    end
    Tuple(forms)
end

"""Number of distinct enzyme states."""
n_states(m::EnzymeMechanism) = length(enzyme_forms(m))

"""
    stoich_matrix(m::EnzymeMechanism) → Matrix{Int}

Full stoichiometry matrix. Rows are species in the order
`(enzyme_forms..., metabolites...)` (use `enzyme_row_range(m)` and
`metabolite_row_range(m)` to slice). Columns are step indices.
Positive = produced; negative = consumed in the forward direction.

Enzyme-row columns sum to zero by construction (each step has one
enzyme on each side).
"""
@generated function stoich_matrix(::EnzymeMechanism{M, R}) where {M, R}
    met_names_set = Set{Symbol}()
    for group in M; for name in group; push!(met_names_set, name); end; end

    seen = Set{Symbol}()
    enz = Symbol[]
    for (lhs, rhs, _, _) in R
        for s in lhs; s ∉ met_names_set && s ∉ seen && (push!(seen, s); push!(enz, s)); end
        for s in rhs; s ∉ met_names_set && s ∉ seen && (push!(seen, s); push!(enz, s)); end
    end

    met_seen = Set{Symbol}()
    mets = Symbol[]
    for group in M
        for name in group
            name ∉ met_seen && (push!(met_seen, name); push!(mets, name))
        end
    end

    species = [enz; mets]
    sp_idx = Dict(s => i for (i, s) in enumerate(species))
    S = zeros(Int, length(species), length(R))
    for (j, (lhs, rhs, _, _)) in enumerate(R)
        for s in lhs; S[sp_idx[s], j] -= 1; end
        for s in rhs; S[sp_idx[s], j] += 1; end
    end
    S
end

enzyme_row_range(m::EnzymeMechanism) = 1:n_states(m)
metabolite_row_range(m::EnzymeMechanism) =
    (n_states(m) + 1):(n_states(m) + length(metabolites(m)))

# ─── AllostericEnzymeMechanism Accessors ────────────────────────

"""Delegate structural accessors to the CatalyticMech singleton."""
n_states(::AllostericEnzymeMechanism{M,CM,CS,RS}) where {M,CM,CS,RS} =
    n_states(CM())
n_steps(::AllostericEnzymeMechanism{M,CM,CS,RS}) where {M,CM,CS,RS} =
    n_steps(CM())
equilibrium_steps(::AllostericEnzymeMechanism{M,CM,CS,RS}) where {M,CM,CS,RS} =
    equilibrium_steps(CM())
substrates(::AllostericEnzymeMechanism{M,CM,CS,RS}) where {M,CM,CS,RS} =
    substrates(CM())
products(::AllostericEnzymeMechanism{M,CM,CS,RS}) where {M,CM,CS,RS} =
    products(CM())
@generated function regulators(
    ::AllostericEnzymeMechanism{M,CM,CS,RS},
) where {M,CM,CS,RS}
    ligs = Symbol[]
    for entry in RS
        for lig in entry[1]
            lig in ligs || push!(ligs, lig)
        end
    end
    Tuple(ligs)
end
param_constraints(::AllostericEnzymeMechanism) = ()

"""
    allosteric_regulators(m::AllostericEnzymeMechanism) → Tuple{Tuple{Symbol,Symbol},...}

Return `(ligand, tag)` pairs derived from `RegSites` membership:
ligand listed in `r_only_ligands` → `:OnlyR`,
in `t_only_ligands` → `:OnlyT`,
in `tr_equiv_ligands` → `:EqualRT`,
absent from all three → `:NonequalRT`.
"""
@generated function allosteric_regulators(
    ::AllostericEnzymeMechanism{M,CM,CS,RS},
) where {M,CM,CS,RS}
    pairs = Tuple{Symbol,Symbol}[]
    seen = Set{Symbol}()
    for entry in RS
        ligs, _, tr_equiv, r_only, t_only = entry
        for lig in ligs
            lig in seen && continue
            push!(seen, lig)
            tag = if lig in r_only
                :OnlyR
            elseif lig in t_only
                :OnlyT
            elseif lig in tr_equiv
                :EqualRT
            else
                :NonequalRT
            end
            push!(pairs, (lig, tag))
        end
    end
    Tuple(pairs)
end

"""Return all metabolite names (catalytic + regulatory) from the Metabolites type param."""
metabolites(::AllostericEnzymeMechanism{Mets,CM,CS,RS}) where {Mets,CM,CS,RS} = Mets
