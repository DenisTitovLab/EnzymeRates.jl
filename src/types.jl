using Graphs

"""Sort species tuples alphabetically by name (first element)."""
_sort_species(t::Tuple) = Tuple(sort(collect(t); by=s -> s[1]))

"""
Normalize species tuples to 3-element
`(name, atoms, max_sites)`, defaulting max_sites to 1.
"""
function _normalize_species_with_sites(t::Tuple)
    Tuple(_normalize_one_species(s) for s in t)
end

function _normalize_one_species(s)
    n = length(s)
    n == 2 && return (s[1], s[2], 1)
    n == 3 && return s
    error("Species tuple must have 2 or 3 elements, got $n")
end

"""
    EnzymeReaction{Substrates, Products, Regulators}

Singleton type encoding an enzyme reaction specification in type parameters.

Each of `Substrates`, `Products`, `Regulators` is a tuple of
`(name::Symbol, atoms::Tuple{Vararg{Tuple{Symbol,Int}}}, max_sites::Int)`.
For backward compatibility, 2-element `(name, atoms)` tuples are auto-normalized
to 3-element with `max_sites=1`.
"""
struct EnzymeReaction{Substrates, Products, Regulators} end

function EnzymeReaction(subs::Tuple, prods::Tuple, regs::Tuple=())
    isempty(subs) && error("Substrates must not be empty")
    isempty(prods) && error("Products must not be empty")
    subs = _normalize_species_with_sites(subs)
    prods = _normalize_species_with_sites(prods)
    regs = _normalize_species_with_sites(regs)
    for group in (subs, prods, regs), s in group
        s[3] >= 1 || error("max_sites must be ≥ 1, got $(s[3]) for $(s[1])")
    end
    subs_names = [s[1] for s in subs]
    prods_names = [s[1] for s in prods]
    length(subs_names) != length(Set(subs_names)) && error("Duplicate substrate names")
    length(prods_names) != length(Set(prods_names)) && error("Duplicate product names")
    subs = _sort_species(subs)
    prods = _sort_species(prods)
    regs = _sort_species(regs)
    EnzymeReaction{subs, prods, regs}()
end

"""
    EnzymeMechanism{Species,Reactions,EquilibriumSteps}

Singleton type encoding an enzyme mechanism in type parameters.

- `Species`: `(substrates, products, regulators, enzyme_species)`
  where each entry is a tuple of
  `(name::Symbol, atoms::Tuple{Vararg{Tuple{Symbol,Int}}})`.
- `Reactions`: tuple of `(lhs, rhs)` where each side is a
  tuple of species `Symbol`s.
- `EquilibriumSteps`: tuple of `Bool` indicating which steps
  are rapid-equilibrium (`true`) vs steady-state (`false`).
"""
struct EnzymeMechanism{
    Species, Reactions, EquilibriumSteps, ParamConstraints,
} end

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
    EnzymeMechanism(species, reactions, eq_steps, [constraints])

Construct an `EnzymeMechanism` from explicit species, reaction tuples, equilibrium step
flags, and optional parameter constraints.

Each constraint is a tuple
`(target::Symbol, coeff::Int, factors::Tuple{...})`.
"""
function EnzymeMechanism(species::Tuple, reactions::Tuple, eq_steps::Tuple{Vararg{Bool}},
                         constraints::Tuple=())
    length(eq_steps) == length(reactions) ||
        error("eq_steps length must match reactions length")
    all(eq_steps) && !isempty(eq_steps) && error(
        "At least one steady-state step is required " *
        "(not all steps can be rapid-equilibrium)",
    )

    length(species) == 4 ||
        error("species must be " *
              "(substrates, products, regulators, enzymes)")
    subs, prods, regs, enzs = species

    for (label, group) in (("substrate", subs), ("product", prods), ("regulator", regs))
        names = [name for (name, _) in group]
        length(names) != length(Set(names)) && error("Duplicate $label names")
    end

    met_atoms = Dict{Symbol,Dict{Symbol,Int}}()
    for (name, atoms) in Iterators.flatten((subs, prods, regs))
        d = Dict{Symbol,Int}(a => c for (a, c) in atoms)
        if haskey(met_atoms, name)
            met_atoms[name] == d || error("Inconsistent atoms for metabolite $name")
        else
            met_atoms[name] = d
        end
    end

    enzyme_atoms = Dict{Symbol,Dict{Symbol,Int}}()
    for (name, atoms) in enzs
        haskey(enzyme_atoms, name) && error("Duplicate enzyme species $name")
        enzyme_atoms[name] = Dict{Symbol,Int}(a => c for (a, c) in atoms)
    end
    isempty(enzyme_atoms) && error("No enzyme species defined")
    for name in keys(enzyme_atoms)
        haskey(met_atoms, name) && error(
            "Species $name defined as both " *
            "enzyme and metabolite",
        )
    end
    free_enzymes = [name for (name, atoms) in enzs if isempty(atoms)]
    isempty(free_enzymes) && error("No free enzyme form (enzyme with empty atoms) defined")

    n_with = count(s -> !isempty(s[2]), Iterators.flatten((subs, prods, regs)))
    n_total = length(subs) + length(prods) + length(regs)
    0 < n_with < n_total && error(
        "All metabolites must either have atoms " *
        "or all lack atoms; found a mix",
    )
    skip_atom_checks = n_with == 0

    expected = Dict{Symbol,Int}()
    for (name, _) in subs; expected[name] = get(expected, name, 0) - 1; end
    for (name, _) in prods; expected[name] = get(expected, name, 0) + 1; end
    for (name, _) in regs; expected[name] = get(expected, name, 0); end

    isempty(reactions) && error("Reactions tuple must not be empty")
    enzyme_set = Set(keys(enzyme_atoms))

    # Canonical step form for RE steps: normalize direction so metabolite
    # binding steps have metabolite on LHS (→ binding Kd convention).
    # SS steps are left unchanged (kf/kr have no automatic inversion).
    reactions = ntuple(length(reactions)) do i
        lhs, rhs = reactions[i]
        if !eq_steps[i]
            return (lhs, rhs)
        end
        rhs_has_met = any(haskey(met_atoms, s) for s in rhs)
        lhs_has_met = any(haskey(met_atoms, s) for s in lhs)
        if rhs_has_met && !lhs_has_met
            (rhs, lhs)
        else
            (lhs, rhs)
        end
    end

    net = Dict{Symbol,Int}()
    for (i, (lhs, rhs)) in enumerate(reactions)
        lhs_enz, lhs_mets, lhs_atoms, lhs_met_names =
            _count_side(lhs, enzyme_set, enzyme_atoms,
                        met_atoms, i)
        rhs_enz, rhs_mets, rhs_atoms, rhs_met_names =
            _count_side(rhs, enzyme_set, enzyme_atoms,
                        met_atoms, i)
        for s in lhs_met_names; net[s] = get(net, s, 0) - 1; end
        for s in rhs_met_names; net[s] = get(net, s, 0) + 1; end
        lhs_enz == 1 || error("Reaction $i lhs must contain exactly one enzyme form")
        rhs_enz == 1 || error("Reaction $i rhs must contain exactly one enzyme form")
        lhs_mets <= 1 || error("Reaction $i lhs has more than one metabolite")
        rhs_mets <= 1 || error("Reaction $i rhs has more than one metabolite")
        if !skip_atom_checks
            filter!(p -> p.second != 0, lhs_atoms)
            filter!(p -> p.second != 0, rhs_atoms)
            lhs_atoms == rhs_atoms || error("Atomic conservation failed at step $i")
        end
    end

    # Net stoichiometry validation: cycle steps give k× the reaction, but dead-end
    # binding steps add extra consumption that may shift or cancel net contributions.
    # We require: each substrate/product appears in at least one reaction,
    # and no unexpected metabolites appear.
    for (name, coeff) in expected
        if coeff != 0
            label = coeff < 0 ? "Substrate" : "Product"
            haskey(net, name) || error(
                "$label $name does not appear " *
                "in any reaction",
            )
        end
    end
    for (name, _) in net
        haskey(expected, name) || error("Metabolite $name not in species tuple")
    end

    length(reactions) != length(Set(reactions)) && error("Duplicate reactions")

    # Enzyme reachability check
    reached = Set{Symbol}([first(free_enzymes)])
    for _ in reactions, (lhs, rhs) in reactions
        e_l = first(s for s in lhs if s in enzyme_set)
        e_r = first(s for s in rhs if s in enzyme_set)
        e_l ∈ reached && push!(reached, e_r)
        e_r ∈ reached && push!(reached, e_l)
    end
    for name in keys(enzyme_atoms)
        name ∈ reached || error("Enzyme form $name is not reachable from free enzyme")
    end

    # Validate constraints
    if !isempty(constraints)
        valid = Set{Symbol}()
        for (i, re) in enumerate(eq_steps)
            if re
                push!(valid, Symbol("K$i"))
            else
                push!(valid, Symbol("k$(i)f"))
                push!(valid, Symbol("k$(i)r"))
            end
        end
        targets = Symbol[]
        for (target, coeff, factors) in constraints
            target ∈ valid || error(
                "Constraint target $target is not a " *
                "valid parameter of this mechanism",
            )
            coeff > 0 || error(
                "Constraint coefficient must be " *
                "positive, got $coeff for $target",
            )
            push!(targets, target)
            for (sym, _) in factors
                sym ∈ valid || error(
                    "Constraint replacement symbol " *
                    "$sym is not a valid parameter " *
                    "of this mechanism",
                )
                sym == target && error(
                    "Self-referencing constraint: " *
                    "$target = ... $target ...",
                )
            end
        end
        dup = findfirst(
            t -> count(==(t), targets) > 1, targets,
        )
        dup !== nothing && error(
            "Duplicate constraint target: " *
            "$(targets[dup])",
        )
    end

    EnzymeMechanism{species, reactions, eq_steps, constraints}()
end

"""
    OligomericEnzymeMechanism{Metabolites, CatalyticMech, CatalyticN, RegSites, NConf}

Singleton type for multi-site, multi-conformation allosteric enzymes.

- `Metabolites`: tuple of `(name, atoms)` pairs from `metabolites:` block
- `CatalyticMech`: `EnzymeMechanism` type for one catalytic subunit
- `CatalyticN`: number of catalytic sites per enzyme molecule
- `RegSites`: tuple of `((ligand_syms...,), multiplicity)` pairs
- `NConf`: number of conformational states (1 = non-cooperative, 2 = two-state MWC)
"""
struct OligomericEnzymeMechanism{
    Metabolites, CatalyticMech, CatalyticN, RegSites, NConf,
} end

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

function Base.show(io::IO, ::EnzymeReaction{S,P,R}) where {S,P,R}
    subs_str = join([string(name) for (name, _) in S], " + ")
    prods_str = join([string(name) for (name, _) in P], " + ")
    print(io, "EnzymeReaction: ", subs_str, " ⇌ ", prods_str)
    if !isempty(R)
        regs_str = join([string(name) for (name, _) in R], ", ")
        print(io, " | regulators: ", regs_str)
    end
end

function Base.show(
    io::IO,
    m::EnzymeMechanism{Species, Reactions, EqSteps, PC},
) where {Species, Reactions, EqSteps, PC}
    subs, prods, regs, enzs = Species
    enz_names = Set(e[1] for e in enzs)

    # Check if mechanism is linear (each enzyme form appears on LHS and RHS at most once)
    lhs_counts = Dict{Symbol,Int}()
    rhs_counts = Dict{Symbol,Int}()
    for (lhs, rhs) in Reactions
        for s in lhs; s in enz_names && (lhs_counts[s] = get(lhs_counts, s, 0) + 1); end
        for s in rhs; s in enz_names && (rhs_counts[s] = get(rhs_counts, s, 0) + 1); end
    end
    is_linear = (all(v <= 1 for v in values(lhs_counts)) &&
                  all(v <= 1 for v in values(rhs_counts)))

    _arrow(is_eq) = is_eq ? " ⇌ " : " <--> "

    if is_linear
        # Compact chain: E + S ⇌ ES <--> E + P
        parts = String[]
        arrows = String[]
        for (i, (lhs, rhs)) in enumerate(Reactions)
            if i == 1
                push!(parts, join(lhs, " + "))
            end
            push!(arrows, _arrow(EqSteps[i]))
            push!(parts, join(rhs, " + "))
        end
        print(io, "EnzymeMechanism: ")
        for (i, part) in enumerate(parts)
            i > 1 && print(io, arrows[i-1])
            print(io, part)
        end
    else
        # Multi-line for branched mechanisms
        n = length(Reactions)
        ne = length(enzs)
        print(io, "EnzymeMechanism (", n, " steps, ", ne, " enzyme forms):")
        for (i, (lhs, rhs)) in enumerate(Reactions)
            print(io, "\n  ", join(lhs, " + "), _arrow(EqSteps[i]), join(rhs, " + "))
        end
    end
    if !isempty(regs)
        regs_str = join([string(name) for (name, _) in regs], ", ")
        print(io, " | regulators: ", regs_str)
    end
    if !isempty(PC)
        cstrs = [
            _user_constraint_to_string(target, coeff, factors)
            for (target, coeff, factors) in PC
        ]
        print(io, " | constraints: ", join(cstrs, ", "))
    end
end

"""Format a user constraint as a string: target = rhs."""
function _user_constraint_to_string(target::Symbol, coeff::Int, factors)
    parts = String[]
    coeff != 1 && push!(parts, string(coeff))
    for (sym, exp) in factors
        if exp == 1
            push!(parts, string(sym))
        elseif exp == -1
            push!(parts, "1 / $sym")
        else
            push!(parts, "$sym^$exp")
        end
    end
    "$target = $(isempty(parts) ? string(coeff) : join(parts, " * "))"
end

# ─── Accessors ─────────────────────────────────────────────────

"""Return substrates (with stoichiometric multiplicity)."""
substrates(::EnzymeMechanism{Species}) where {Species} = Species[1]
substrates(::EnzymeReaction{S,P,R}) where {S,P,R} = S

"""Return products (with stoichiometric multiplicity)."""
products(::EnzymeMechanism{Species}) where {Species} = Species[2]
products(::EnzymeReaction{S,P,R}) where {S,P,R} = P

"""Return regulators."""
regulators(::EnzymeMechanism{Species}) where {Species} = Species[3]
regulators(::EnzymeReaction{S,P,R}) where {S,P,R} = R

"""Return all enzyme forms as a tuple of (name, atoms)."""
enzyme_forms(::EnzymeMechanism{Species}) where {Species} = Species[4]

"""Compile-time helper: collect unique metabolites from Species type parameter."""
function _unique_metabolites(Species)
    subs, prods, regs = Species[1:3]
    seen = Set{Symbol}()
    mets = Tuple{Symbol,Any}[]
    for group in (subs, prods, regs)
        for (name, atoms) in group
            if name ∉ seen
                push!(seen, name)
                push!(mets, (name, atoms))
            end
        end
    end
    return mets
end

"""Return unique metabolites as a tuple of (name, atoms) — internal use."""
@generated function _metabolites_with_sites(::EnzymeMechanism{Species}) where {Species}
    return Tuple(_unique_metabolites(Species))
end

"""
    metabolites(m::EnzymeMechanism) → Tuple{Symbol,...}

Return distinct metabolite names as a tuple of Symbols.
"""
@generated function metabolites(::EnzymeMechanism{Species}) where {Species}
    mets = _unique_metabolites(Species)
    Tuple(m[1] for m in mets)
end

"""Number of distinct enzyme states."""
n_states(::EnzymeMechanism{Species}) where {Species} = length(Species[4])

"""Number of steps in the mechanism."""
function n_steps(
    ::EnzymeMechanism{Species, Reactions},
) where {Species, Reactions}
    length(Reactions)
end

"""Return the reactions tuple directly."""
reactions(::EnzymeMechanism{Species, R}) where {Species, R} = R

"""Return the equilibrium steps tuple (true = rapid-equilibrium, false = steady-state)."""
equilibrium_steps(::EnzymeMechanism{Sp, Rx, Eq}) where {Sp, Rx, Eq} = Eq

"""Return the parameter constraints tuple."""
param_constraints(::EnzymeMechanism{Sp, Rx, Eq, PC}) where {Sp, Rx, Eq, PC} = PC

"""
Build a directed graph of enzyme-form connectivity.
Returns (graph, enzyme_forms_tuple).
"""
@generated function graph(::EnzymeMechanism{Species, Reactions}) where {Species, Reactions}
    enzs = Species[4]
    enz_names = Tuple(e[1] for e in enzs)
    name_to_idx = Dict(n => i for (i, n) in enumerate(enz_names))
    enz_set = Set(enz_names)
    g = SimpleDiGraph(length(enzs))
    for (lhs, rhs) in Reactions
        e_lhs = first(s for s in lhs if s in enz_set)
        e_rhs = first(s for s in rhs if s in enz_set)
        add_edge!(g, name_to_idx[e_lhs], name_to_idx[e_rhs])
        add_edge!(g, name_to_idx[e_rhs], name_to_idx[e_lhs])
    end
    return g, enzs
end

"""
Stoichiometry matrix: rows = metabolites, columns = steps.
Positive = produced, negative = consumed.
"""
@generated function stoich_matrix(
    ::EnzymeMechanism{Species, Reactions, EqSteps},
) where {Species, Reactions, EqSteps}
    mets = _unique_metabolites(Species)
    met_idx = Dict(m[1] => i for (i, m) in enumerate(mets))
    enz_names = Set(e[1] for e in Species[4])
    S = zeros(Int, length(mets), length(Reactions))
    for (step_j, (lhs, rhs)) in enumerate(Reactions)
        for s in lhs
            s in enz_names && continue
            S[met_idx[s], step_j] -= 1
        end
        for s in rhs
            s in enz_names && continue
            S[met_idx[s], step_j] += 1
        end
    end
    return S
end
