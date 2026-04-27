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

    # Each substrate / product must appear in some step. Regulators are
    # optional bindings (a spec may list them before any expansion move
    # has added their dead-end or allosteric binding steps).
    appears = Set{Symbol}()
    for (lhs, rhs, _, _) in rxns
        for s in lhs; push!(appears, s); end
        for s in rhs; push!(appears, s); end
    end
    for name in vcat(collect(subs), collect(prods))
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
    AllostericEnzymeMechanism{CatalyticMech, CatSites, RegSites}

Singleton type encoding a multi-subunit MWC allosteric enzyme.

- `CatalyticMech`: an `EnzymeMechanism` type (single-subunit catalytic mech).
- `CatSites`: `(multiplicity::Int, group_tags::Tuple{Pair{Int,Symbol}...})`.
  Non-default-only storage; absent groups have tag `:NonequalRT`.
- `RegSites`: tuple of entries `((ligands, multiplicity, ligand_tags),)`.
  One entry per reg site.
"""
struct AllostericEnzymeMechanism{
    CatalyticMech, CatSites, RegSites,
} <: AbstractEnzymeMechanism end

"""
    AllostericEnzymeMechanism(cm, cat_sites, reg_sites)

Build an `AllostericEnzymeMechanism` from a catalytic `EnzymeMechanism`,
a `(multiplicity, group_tags)` pair, and a tuple of reg-site entries.

`group_tags` and ligand-tag entries store only non-default tags; absent
entries default to `:NonequalRT`. `group_tags` is sorted by group id for
canonical type identity. Validates:
  - tag values are one of `:OnlyR`, `:OnlyT`, `:EqualRT`, `:NonequalRT`,
  - `group_tags` reference existing kinetic groups,
  - iso-only kinetic groups are not tagged `:OnlyT`,
  - reg sites have at least one ligand and at least one non-`:EqualRT`
    ligand (single-/all-`:EqualRT` reg sites cancel identically).
"""
function AllostericEnzymeMechanism(
    cm::EnzymeMechanism,
    cat_sites::Tuple{Int, <:Tuple},
    reg_sites::Tuple,
)
    multiplicity, group_tags = cat_sites
    valid_groups = Set(kinetic_groups(cm))
    rxns = reactions(cm)
    cat_mets = Set(metabolites(cm))

    for (g, tag) in group_tags
        g in valid_groups ||
            error("group_tag references non-existent kinetic_group $g")
        tag in (:OnlyR, :OnlyT, :EqualRT, :NonequalRT) ||
            error("Invalid group tag: $tag")
        any_iso = false
        for idx in steps_in_group(cm, g)
            lhs, rhs, is_eq, _ = rxns[idx]
            mets_in = any(s in cat_mets for s in (lhs..., rhs...))
            if !is_eq && !mets_in
                any_iso = true
            end
        end
        any_iso && tag == :OnlyT &&
            error("Iso group $g tagged :OnlyT is forbidden " *
                  "(R-inactive is a relabel)")
    end

    for (i, entry) in enumerate(reg_sites)
        ligands, _, lig_tags = entry
        isempty(ligands) && error("Reg site $i has no ligands")
        tag_map = Dict(lig_tags)
        for (lig, tag) in lig_tags
            tag in (:OnlyR, :OnlyT, :EqualRT, :NonequalRT) ||
                error("Invalid reg-site tag $tag for ligand $lig")
        end
        all_eq = all(get(tag_map, l, :NonequalRT) == :EqualRT for l in ligands)
        all_eq &&
            error("Reg site $i with all `:EqualRT` ligands cancels " *
                  "identically (or single-ligand :EqualRT reg site); at " *
                  "least one ligand must have a non-:EqualRT tag. " *
                  "Ligands: $ligands")
    end

    sorted_tags = Tuple(sort(collect(group_tags); by=first))
    cat_sites_canon = (multiplicity, sorted_tags)

    AllostericEnzymeMechanism{typeof(cm), cat_sites_canon, reg_sites}()
end

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

function Base.show(
    io::IO, m::EnzymeMechanism{Mets, Rxns},
) where {Mets, Rxns}
    _, _, regs = Mets
    enz_set = Set(enzyme_forms(m))

    # Linear mechanism iff each enzyme form appears at most once on either side.
    lhs_counts = Dict{Symbol,Int}()
    rhs_counts = Dict{Symbol,Int}()
    for (lhs, rhs, _, _) in Rxns
        for s in lhs; s in enz_set && (lhs_counts[s] = get(lhs_counts, s, 0) + 1); end
        for s in rhs; s in enz_set && (rhs_counts[s] = get(rhs_counts, s, 0) + 1); end
    end
    is_linear = all(v <= 1 for v in values(lhs_counts)) &&
                all(v <= 1 for v in values(rhs_counts))
    _arrow(is_eq) = is_eq ? " ⇌ " : " <--> "

    if is_linear
        print(io, "EnzymeMechanism: ")
        for (i, (lhs, rhs, is_eq, _)) in enumerate(Rxns)
            i == 1 && print(io, join(lhs, " + "))
            print(io, _arrow(is_eq), join(rhs, " + "))
        end
    else
        print(io, "EnzymeMechanism (", length(Rxns), " steps, ",
              length(enz_set), " enzyme forms):")
        for (lhs, rhs, is_eq, _) in Rxns
            print(io, "\n  ", join(lhs, " + "), _arrow(is_eq),
                      join(rhs, " + "))
        end
    end
    if !isempty(regs)
        print(io, " | regulators: ", join(regs, ", "))
    end
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

catalytic_mechanism(::AllostericEnzymeMechanism{CM}) where {CM} = CM()
catalytic_multiplicity(::AllostericEnzymeMechanism{CM, CS}) where {CM, CS} = CS[1]

function group_tag(::AllostericEnzymeMechanism{CM, CS, RS}, g::Int) where {CM, CS, RS}
    for (k, t) in CS[2]
        k == g && return t
    end
    :NonequalRT
end

step_tag(m::AllostericEnzymeMechanism, idx::Int) =
    group_tag(m, kinetic_group(catalytic_mechanism(m), idx))

substrates(m::AllostericEnzymeMechanism)         = substrates(catalytic_mechanism(m))
products(m::AllostericEnzymeMechanism)           = products(catalytic_mechanism(m))
reactions(m::AllostericEnzymeMechanism)          = reactions(catalytic_mechanism(m))
equilibrium_steps(m::AllostericEnzymeMechanism)  = equilibrium_steps(catalytic_mechanism(m))
n_steps(m::AllostericEnzymeMechanism)            = n_steps(catalytic_mechanism(m))
enzyme_forms(m::AllostericEnzymeMechanism)       = enzyme_forms(catalytic_mechanism(m))
n_states(m::AllostericEnzymeMechanism)           = n_states(catalytic_mechanism(m))
kinetic_group(m::AllostericEnzymeMechanism, i::Int) =
    kinetic_group(catalytic_mechanism(m), i)
kinetic_groups(m::AllostericEnzymeMechanism)     = kinetic_groups(catalytic_mechanism(m))
steps_in_group(m::AllostericEnzymeMechanism, g)  =
    steps_in_group(catalytic_mechanism(m), g)
stoich_matrix(m::AllostericEnzymeMechanism)      = stoich_matrix(catalytic_mechanism(m))
enzyme_row_range(m::AllostericEnzymeMechanism)   = enzyme_row_range(catalytic_mechanism(m))
metabolite_row_range(m::AllostericEnzymeMechanism) =
    metabolite_row_range(catalytic_mechanism(m))

# Returns ONLY reg-site ligands, NOT a union with catalytic_mechanism's
# regulators. Downstream rate-equation code reads `regulators(m)` to find
# dead-end binding K's; including allosteric-only ligands would cause it
# to look up nonexistent K names.
regulators(::AllostericEnzymeMechanism{CM, CS, RS}) where {CM, CS, RS} = begin
    syms = Symbol[]
    seen = Set{Symbol}()
    for entry in RS
        for lig in entry[1]
            lig in seen || (push!(seen, lig); push!(syms, lig))
        end
    end
    Tuple(syms)
end

metabolites(::AllostericEnzymeMechanism{CM, CS, RS}) where {CM, CS, RS} = begin
    cat_mets = metabolites(CM())
    extra = Symbol[]
    seen = Set{Symbol}(cat_mets)
    for entry in RS
        for lig in entry[1]
            lig in seen || (push!(seen, lig); push!(extra, lig))
        end
    end
    (cat_mets..., extra...)
end

allosteric_regulators(::AllostericEnzymeMechanism{CM, CS, RS}) where {CM, CS, RS} = begin
    result = Tuple{Symbol, Symbol}[]
    for (ligands, _, lig_tags) in RS
        tag_map = Dict(lig_tags)
        for lig in ligands
            push!(result, (lig, get(tag_map, lig, :NonequalRT)))
        end
    end
    Tuple(result)
end

catalytic_inhibitors(::AllostericEnzymeMechanism{CM, CS, RS}) where {CM, CS, RS} = begin
    rs_names = Set{Symbol}()
    for (ligs, _, _) in RS
        for l in ligs; push!(rs_names, l); end
    end
    cat_regs = regulators(CM())
    Tuple(r for r in cat_regs if r ∉ rs_names)
end

regulatory_sites(::AllostericEnzymeMechanism{CM, CS, RS}) where {CM, CS, RS} = RS
regulatory_site_ligands(m::AllostericEnzymeMechanism, i::Int)     =
    regulatory_sites(m)[i][1]
regulatory_site_multiplicity(m::AllostericEnzymeMechanism, i::Int) =
    regulatory_sites(m)[i][2]

function regulatory_ligand_tag(m::AllostericEnzymeMechanism, i::Int, lig::Symbol)
    for (k, t) in regulatory_sites(m)[i][3]
        k == lig && return t
    end
    :NonequalRT
end
