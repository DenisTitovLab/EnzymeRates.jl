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

"""Sum element counts across a tuple of `(name, atoms)` pairs.
Returns a Dict{Symbol,Int}. Errors if any species's atoms tuple
is empty (atoms are mandatory) or if any per-atom count is not a
positive Int."""
function _sum_atoms(species::Tuple, side::String)
    totals = Dict{Symbol,Int}()
    for (name, atoms) in species
        isempty(atoms) && error(
            "EnzymeReaction: $side metabolite $name has no declared " *
            "atoms; atoms are mandatory (use `[C…]` bracket syntax in " *
            "@enzyme_reaction or pass non-empty atom tuples to the " *
            "constructor).")
        for (elem, count) in atoms
            count isa Integer && !(count isa Bool) && count > 0 ||
                error(
                "EnzymeReaction: $side metabolite $name has " *
                "non-positive atom count for element $elem ($count); " *
                "atom counts must be positive integers (not Bool).")
            totals[elem] = get(totals, elem, 0) + count
        end
    end
    totals
end

function EnzymeReaction(subs::Tuple, prods::Tuple, regs::Tuple=(); oligomeric_state::Int=1)
    isempty(subs) && error("Substrates must not be empty")
    isempty(prods) && error("Products must not be empty")
    subs_names = [s[1] for s in subs]
    prods_names = [s[1] for s in prods]
    length(subs_names) != length(Set(subs_names)) &&
        error("Duplicate substrate names")
    length(prods_names) != length(Set(prods_names)) &&
        error("Duplicate product names")
    sub_atoms = _sum_atoms(subs, "substrate")
    prod_atoms = _sum_atoms(prods, "product")
    all_elems = union(keys(sub_atoms), keys(prod_atoms))
    for elem in all_elems
        s_count = get(sub_atoms, elem, 0)
        p_count = get(prod_atoms, elem, 0)
        s_count == p_count || error(
            "EnzymeReaction: atom imbalance — element $elem appears " *
            "$s_count time(s) on substrate side and $p_count time(s) " *
            "on product side. Declared atoms must balance.")
    end
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

    # Every substrate, product, AND regulator must appear in some step.
    appears = Set{Symbol}()
    for (lhs, rhs, _, _) in rxns
        for s in lhs; push!(appears, s); end
        for s in rhs; push!(appears, s); end
    end
    for name in vcat(collect(subs), collect(prods), collect(regs))
        name in appears ||
            error("Listed metabolite or regulator $name does not " *
                  "appear in any reaction step")
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

Multi-subunit MWC allosteric enzyme. `CatalyticMech` is an
`EnzymeMechanism` describing one catalytic subunit's cycle.

# Type parameters
- `CatSites`: `(multiplicity::Int, cat_allo_states::Tuple{Symbol...})`.
  `cat_allo_states[g]` is the allosteric state of catalytic kinetic
  group `g` (1-indexed, dense — every group must have an entry).
  Allowed values: `:EqualRT`, `:NonequalRT`, `:OnlyR`. `:OnlyT`
  catalytic groups error during construction (R-state-active
  convention).
- `RegSites`: tuple of regulator-site entries
  `(ligands::Tuple{Symbol...}, multiplicity::Int,
   reg_allo_states::Tuple{Symbol...})` where `reg_allo_states` is
  parallel to `ligands`. Allowed values: all four states
  (`:EqualRT`, `:NonequalRT`, `:OnlyR`, `:OnlyT`).

Constructor validates:
- Catalytic state count matches kinetic-group count.
- Regulator state tuple length matches ligand tuple length at each site.
- No catalytic group has `:OnlyT` state.
- At least one ligand at each reg site is non-`:EqualRT`.
"""
struct AllostericEnzymeMechanism{
    CatalyticMech, CatSites, RegSites,
} <: AbstractEnzymeMechanism end

function AllostericEnzymeMechanism(
    cm::EnzymeMechanism, cat_sites::Tuple, reg_sites::Tuple,
)
    multiplicity, cat_allo_states = cat_sites
    multiplicity isa Int && multiplicity ≥ 1 ||
        error("Catalytic multiplicity must be a positive Int, got $multiplicity")

    n_groups = length(unique(kinetic_group(cm, i) for i in 1:n_steps(cm)))
    # Validate kinetic_group numbers are 1..n_groups consecutive — the
    # cat_allo_states tuple is indexed by group number, so non-consecutive
    # numbering would cause OOB or wrong-state lookup at runtime.
    observed_groups = sort!(unique(kinetic_group(cm, i) for i in 1:n_steps(cm)))
    observed_groups == collect(1:n_groups) ||
        error("Catalytic mechanism kinetic_group numbers must be 1..n " *
              "consecutive; got $observed_groups")
    length(cat_allo_states) == n_groups ||
        error("cat_allo_states length $(length(cat_allo_states)) does not " *
              "match catalytic kinetic-group count $n_groups")
    for (g, st) in enumerate(cat_allo_states)
        st === :OnlyT &&
            error("Catalytic kinetic group $g has state :OnlyT; the " *
                  "R-state is the active state by convention. Relabel " *
                  "your mechanism so the active state is R (use :OnlyR " *
                  "instead).")
        st in (:OnlyR, :EqualRT, :NonequalRT) ||
            error("Catalytic kinetic group $g has unknown allo state $st; " *
                  "must be one of (:OnlyR, :EqualRT, :NonequalRT)")
    end

    for (i, entry) in enumerate(reg_sites)
        ligands, n_reg, reg_allo_states = entry
        ligands isa Tuple && all(l isa Symbol for l in ligands) ||
            error("Reg site $i: ligands must be a Tuple of Symbol")
        length(ligands) >= 1 ||
            error("Reg site $i: must have at least one ligand; got empty " *
                  "ligand tuple")
        n_reg isa Int && n_reg ≥ 1 ||
            error("Reg site $i: multiplicity must be a positive Int")
        length(reg_allo_states) == length(ligands) ||
            error("Reg site $i: reg_allo_states length $(length(reg_allo_states)) " *
                  "does not match ligand count $(length(ligands))")
        for (k, st) in enumerate(reg_allo_states)
            st in (:OnlyR, :OnlyT, :EqualRT, :NonequalRT) ||
                error("Reg site $i, ligand $(ligands[k]): unknown allo state $st")
        end
        # All-:EqualRT site cancels identically — error
        all(st === :EqualRT for st in reg_allo_states) &&
            error("Reg site $i: all ligands are :EqualRT, which produces " *
                  "Q_reg_R == Q_reg_T — no allosteric effect. At least one " *
                  "ligand must be :OnlyR, :OnlyT, or :NonequalRT.")
    end

    AllostericEnzymeMechanism{typeof(cm), cat_sites, reg_sites}()
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
    _arrow(is_eq) = is_eq ? " ⇌ " : " <--> "

    # Walk the steps in source order. Each step shares one enzyme form with
    # the previous step's "outgoing" form; the other side is the new
    # outgoing form. If any step has no shared form, the mechanism is
    # branched and we fall back to multi-line rendering. RE binding
    # canonicalization may put a step's "outgoing" side on the LHS — we
    # detect that case and emit the LHS instead of the RHS.
    chain_segments = String[]
    chain_arrows = String[]
    is_linear = !isempty(Rxns)
    current = nothing
    for (i, (lhs, rhs, is_eq, _)) in enumerate(Rxns)
        e_l = first(s for s in lhs if s in enz_set)
        e_r = first(s for s in rhs if s in enz_set)
        if i == 1
            push!(chain_segments, join(lhs, " + "))
            push!(chain_arrows, _arrow(is_eq))
            push!(chain_segments, join(rhs, " + "))
            current = e_r
        elseif current == e_l
            push!(chain_arrows, _arrow(is_eq))
            push!(chain_segments, join(rhs, " + "))
            current = e_r
        elseif current == e_r
            push!(chain_arrows, _arrow(is_eq))
            push!(chain_segments, join(lhs, " + "))
            current = e_l
        else
            is_linear = false
            break
        end
    end

    if is_linear
        print(io, "EnzymeMechanism: ")
        print(io, chain_segments[1])
        for k in 2:length(chain_segments)
            print(io, chain_arrows[k-1], chain_segments[k])
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

"""Render the catalytic mechanism's steps as multi-line text,
grouping steps that share a kinetic_group with parens and a single
`:: Tag` annotation. Mirrors `@allosteric_mechanism` macro syntax."""
function _format_allo_step_groups(
    io::IO, cm::EnzymeMechanism,
    m::AllostericEnzymeMechanism,
)
    rxns = reactions(cm)
    _arrow(is_eq) = is_eq ? " ⇌ " : " <--> "

    groups_seen = Int[]
    group_to_step_idxs = Dict{Int,Vector{Int}}()
    for (i, step) in enumerate(rxns)
        g = step[4]
        if !haskey(group_to_step_idxs, g)
            push!(groups_seen, g)
            group_to_step_idxs[g] = Int[]
        end
        push!(group_to_step_idxs[g], i)
    end

    for g in groups_seen
        idxs = group_to_step_idxs[g]
        tag = cat_allo_state(m, g)
        if length(idxs) == 1
            (lhs, rhs, is_eq, _) = rxns[idxs[1]]
            print(io, "\n  ", join(lhs, " + "),
                  _arrow(is_eq), join(rhs, " + "),
                  " :: ", tag)
        else
            print(io, "\n  (")
            for (k, i) in enumerate(idxs)
                k > 1 && print(io, ", ")
                (lhs, rhs, is_eq, _) = rxns[i]
                print(io, join(lhs, " + "),
                      _arrow(is_eq), join(rhs, " + "))
            end
            print(io, ") :: ", tag)
        end
    end
end

function Base.show(io::IO, m::AllostericEnzymeMechanism)
    cm = catalytic_mechanism(m)
    print(io, "AllostericEnzymeMechanism (cat_n=",
          catalytic_multiplicity(m))
    rs = regulatory_sites(m)
    if !isempty(rs)
        print(io, ", ", length(rs), " reg sites")
    end
    print(io, "):")
    _format_allo_step_groups(io, cm, m)
    for (i, (ligands, mult, reg_allo_states)) in enumerate(rs)
        print(io, "\n  reg site $i (n=", mult, "): ",
              join(ligands, ", "))
        print(io, " [")
        print(io, join(("$(n)::$(t)"
                        for (n, t) in zip(ligands, reg_allo_states)),
                       ", "))
        print(io, "]")
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

"""Return the allosteric state of catalytic kinetic group `g`."""
function cat_allo_state(::AllostericEnzymeMechanism{CM, CS, RS}, g::Int) where {CM, CS, RS}
    _, states = CS
    return states[g]
end

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
    for (ligands, _, reg_allo_states) in RS
        for (lig, st) in zip(ligands, reg_allo_states)
            push!(result, (lig, st))
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

"""Return the allosteric state of regulator ligand `lig` at site `site_idx`."""
function reg_allo_state(
    ::AllostericEnzymeMechanism{CM, CS, RS}, site_idx::Int, lig::Symbol,
) where {CM, CS, RS}
    ligands, _, states = RS[site_idx]
    idx = findfirst(==(lig), ligands)
    idx === nothing && error("Ligand $lig not at regulatory site $site_idx")
    return states[idx]
end
