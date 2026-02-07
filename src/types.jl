using Graphs

"""Sort species tuples alphabetically by name (first element)."""
_sort_species(t::Tuple) = Tuple(sort(collect(t); by = s -> s[1]))

"""
    AbstractEnzymeReaction

Abstract supertype for enzyme reactions with different type parameters.
"""
abstract type AbstractEnzymeReaction end

"""
    EnzymeReaction{Substrates, Products, Regulators}

Singleton type encoding an enzyme reaction specification in type parameters.

Each of `Substrates`, `Products`, `Regulators` is a tuple of
`(name::Symbol, atoms::Tuple{Vararg{Tuple{Symbol,Int}}})`.
"""
struct EnzymeReaction{Substrates, Products, Regulators} <: AbstractEnzymeReaction end

function EnzymeReaction(subs::Tuple, prods::Tuple, regs::Tuple=())
    isempty(subs) && error("Substrates must not be empty")
    isempty(prods) && error("Products must not be empty")
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
    AbstractEnzymeMechanism

Abstract supertype for enzyme mechanisms. Used as element type in collections
of mechanisms with different type parameters.
"""
abstract type AbstractEnzymeMechanism end

"""
    EnzymeMechanism{Species,Reactions,EquilibriumSteps}

Singleton type encoding an enzyme mechanism in type parameters.

- `Species`: `(substrates, products, regulators, enzyme_species)` where each entry is a tuple
  of `(name::Symbol, atoms::Tuple{Vararg{Tuple{Symbol,Int}}})`.
- `Reactions`: tuple of `(lhs, rhs)` where each side is a tuple of species `Symbol`s.
- `EquilibriumSteps`: tuple of `Bool` indicating which steps are rapid-equilibrium (`true`)
  vs steady-state (`false`).
"""
struct EnzymeMechanism{Species, Reactions, EquilibriumSteps} <: AbstractEnzymeMechanism end

"""Count enzymes, metabolites, atoms, and metabolite names on one side of a reaction."""
function _count_side(side, enzyme_set, enzyme_atoms, met_atoms, step_idx)
    n_enz, n_met, atoms, mets = 0, 0, Dict{Symbol,Int}(), Symbol[]
    for s in side
        if s in enzyme_set
            n_enz += 1
            for (a, c) in enzyme_atoms[s]; atoms[a] = get(atoms, a, 0) + c; end
        elseif haskey(met_atoms, s)
            n_met += 1; push!(mets, s)
            for (a, c) in met_atoms[s]; atoms[a] = get(atoms, a, 0) + c; end
        else
            error("Reaction $(step_idx) uses unknown species $(s)")
        end
    end
    (n_enz, n_met, atoms, mets)
end

"""
    EnzymeMechanism(species::Tuple, reactions::Tuple, eq_steps::Tuple{Vararg{Bool}})

Construct an `EnzymeMechanism` from explicit species, reaction tuples, and equilibrium step flags.

- `species` must be `(substrates, products, regulators, enzymes)` where each entry is a tuple
  of `(name::Symbol, atoms::Tuple{Vararg{Tuple{Symbol,Int}}})`.
- `reactions` is a tuple of `(lhs, rhs)`; each side is a tuple of symbols.
- `eq_steps` is a tuple of `Bool` of the same length as `reactions`, where `true` marks a
  rapid-equilibrium step and `false` marks a steady-state step.
"""
function EnzymeMechanism(species::Tuple, reactions::Tuple, eq_steps::Tuple{Vararg{Bool}}=ntuple(Returns(false), length(reactions)))
    # 0. Validate eq_steps length
    length(eq_steps) == length(reactions) || error("eq_steps length must match reactions length")
    # At least one SS step required
    all(eq_steps) && !isempty(eq_steps) && error("At least one steady-state step is required (not all steps can be rapid-equilibrium)")

    # 1. Validate species tuple structure
    length(species) == 4 || error("species must be (substrates, products, regulators, enzymes)")
    subs, prods, regs, enzs = species
    subs isa Tuple || error("substrates must be a tuple of (name, atoms)")
    prods isa Tuple || error("products must be a tuple of (name, atoms)")
    regs isa Tuple || error("regulators must be a tuple of (name, atoms)")
    enzs isa Tuple || error("enzymes must be a tuple of (name, atoms)")

    # Check for duplicate names within each metabolite category
    subs_names = [name for (name, _) in subs]
    prods_names = [name for (name, _) in prods]
    regs_names = [name for (name, _) in regs]
    length(subs_names) != length(Set(subs_names)) && error("Duplicate substrate names")
    length(prods_names) != length(Set(prods_names)) && error("Duplicate product names")
    length(regs_names) != length(Set(regs_names)) && error("Duplicate regulator names")

    # 2. Check atom consistency across metabolite definitions
    met_atoms = Dict{Symbol,Dict{Symbol,Int}}()
    function add_met(name, atoms)
        atoms_dict = Dict{Symbol,Int}(a => c for (a, c) in atoms)
        if haskey(met_atoms, name)
            met_atoms[name] == atoms_dict || error("Inconsistent atoms for metabolite $(name)")
        else
            met_atoms[name] = atoms_dict
        end
    end
    for (name, atoms) in subs
        add_met(name, atoms)
    end
    for (name, atoms) in prods
        add_met(name, atoms)
    end
    for (name, atoms) in regs
        add_met(name, atoms)
    end

    # 3. Validate enzyme species (no duplicates, not empty, no overlap with metabolites)
    enzyme_atoms = Dict{Symbol,Dict{Symbol,Int}}()
    for (name, atoms) in enzs
        haskey(enzyme_atoms, name) && error("Duplicate enzyme species $(name)")
        enzyme_atoms[name] = Dict{Symbol,Int}(a => c for (a, c) in atoms)
    end
    isempty(enzyme_atoms) && error("No enzyme species defined")

    for name in keys(enzyme_atoms)
        haskey(met_atoms, name) && error("Species $(name) defined as both enzyme and metabolite")
    end

    # Check free enzyme existence
    free_enzymes = [name for (name, atoms) in enzs if isempty(atoms)]
    isempty(free_enzymes) && error("No free enzyme form (enzyme with empty atoms) defined")

    # Check all-or-nothing atoms for metabolites (substrates + products + regulators)
    n_with = count(s -> !isempty(s[2]), Iterators.flatten((subs, prods, regs)))
    n_total = length(subs) + length(prods) + length(regs)
    if 0 < n_with < n_total
        error("All metabolites must either have atoms or all lack atoms; found a mix")
    end
    skip_atom_checks = n_with == 0

    # 4. Compute expected net stoichiometry from species lists
    expected = Dict{Symbol,Int}()
    for name in subs_names
        expected[name] = get(expected, name, 0) - 1
    end
    for name in prods_names
        expected[name] = get(expected, name, 0) + 1
    end
    for name in regs_names
        expected[name] = get(expected, name, 0) + 0
    end

    # 5. Validate reactions not empty
    isempty(reactions) && error("Reactions tuple must not be empty")

    # 6. Per-reaction validation
    enzyme_set = Set(keys(enzyme_atoms))
    net = Dict{Symbol,Int}()

    for (step_idx, reaction) in enumerate(reactions)
        reaction isa Tuple || error("Reaction $(step_idx) is not a tuple")
        length(reaction) == 2 || error("Reaction $(step_idx) must be (lhs, rhs)")
        lhs, rhs = reaction
        lhs isa Tuple || error("Reaction $(step_idx) lhs must be a tuple")
        rhs isa Tuple || error("Reaction $(step_idx) rhs must be a tuple")

        lhs_enz, lhs_mets, lhs_atoms, lhs_met_names = _count_side(lhs, enzyme_set, enzyme_atoms, met_atoms, step_idx)
        rhs_enz, rhs_mets, rhs_atoms, rhs_met_names = _count_side(rhs, enzyme_set, enzyme_atoms, met_atoms, step_idx)

        for s in lhs_met_names; net[s] = get(net, s, 0) - 1; end
        for s in rhs_met_names; net[s] = get(net, s, 0) + 1; end

        lhs_enz == 1 || error("Reaction $(step_idx) lhs must contain exactly one enzyme form")
        rhs_enz == 1 || error("Reaction $(step_idx) rhs must contain exactly one enzyme form")
        lhs_mets <= 1 || error("Reaction $(step_idx) lhs has more than one metabolite")
        rhs_mets <= 1 || error("Reaction $(step_idx) rhs has more than one metabolite")

        if !skip_atom_checks
            filter!(p -> p.second != 0, lhs_atoms)
            filter!(p -> p.second != 0, rhs_atoms)
            lhs_atoms == rhs_atoms || error("Atomic conservation failed at step $(step_idx)")
        end
    end

    # 7. Overall net stoichiometry validation
    for (name, coeff) in expected
        net_coeff = get(net, name, 0)
        if coeff == 0
            net_coeff == 0 || error("Regulator $(name) has nonzero net stoichiometry")
        else
            net_coeff == 0 && error("Net stoichiometry mismatch for $(name)")
            sign(net_coeff) == sign(coeff) || error("Net stoichiometry mismatch for $(name)")
            abs(net_coeff) % abs(coeff) == 0 || error("Net stoichiometry mismatch for $(name)")
        end
    end
    for (name, _) in net
        haskey(expected, name) || error("Metabolite $(name) not in species tuple")
    end

    # 8. Canonical ordering and normalization
    sorted_species = (_sort_species(subs), _sort_species(prods), _sort_species(regs), _sort_species(enzs))

    # Normalize each reaction side so enzyme symbol comes first
    _norm(side) = Tuple(sort(collect(side); by = s -> s in enzyme_set ? Symbol("") : s))
    rxns = [(_norm(lhs), _norm(rhs)) for (lhs, rhs) in reactions]

    # Check for duplicate reactions (after normalization)
    length(rxns) != length(Set(rxns)) && error("Duplicate reactions")

    # Compute each enzyme form's distance from free enzyme along the reaction pathway
    _enz(side) = first(s for s in side if s in enzyme_set)
    free_enz = first(free_enzymes)
    depth = Dict{Symbol,Int}(free_enz => 0)
    for _ in rxns, r in rxns
        haskey(depth, _enz(r[1])) && !haskey(depth, _enz(r[2])) && (depth[_enz(r[2])] = depth[_enz(r[1])] + 1)
        haskey(depth, _enz(r[2])) && !haskey(depth, _enz(r[1])) && (depth[_enz(r[1])] = depth[_enz(r[2])] + 1)
    end

    # Check all enzyme forms are reachable from free enzyme
    for name in keys(enzyme_atoms)
        haskey(depth, name) || error("Enzyme form $(name) is not reachable from free enzyme")
    end

    # Sort reactions by LHS enzyme depth, then alphabetically by LHS metabolites
    # Sort eq_steps together with reactions
    rxn_pairs = collect(zip(rxns, collect(eq_steps)))
    sort!(rxn_pairs; by = pair -> (depth[_enz(pair[1][1])], sort([s for s in pair[1][1] if s ∉ enzyme_set])))
    sorted_rxns = [p[1] for p in rxn_pairs]
    sorted_eq = Tuple(p[2] for p in rxn_pairs)

    EnzymeMechanism{sorted_species, Tuple(sorted_rxns), sorted_eq}()
end

# --- Rate equation mode types ---

"""
    RateEquationMode

Abstract supertype for rate equation parameterization modes.
Controls which form of the rate equation is used and what parameters are expected.
"""
abstract type RateEquationMode end

"""
    RawMode <: RateEquationMode

Raw rate equation mode using all 2N microscopic rate constants (k1f, k1r, k2f, k2r, ...).
No thermodynamic constraints applied. Parameters: all k's + E_total.
"""
struct RawMode <: RateEquationMode end

"""
    HaldaneWegscheiderMode <: RateEquationMode

Rate equation with Haldane-Wegscheider thermodynamic constraints applied.
Dependent parameters are substituted in terms of independent k's and Keq.
Parameters: independent k's + Keq + E_total.
"""
struct HaldaneWegscheiderMode <: RateEquationMode end

"""Singleton instance for raw mode."""
const Raw = RawMode()

"""Singleton instance for Haldane-Wegscheider mode."""
const HaldaneWegscheider = HaldaneWegscheiderMode()

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

function Base.show(io::IO, ::EnzymeMechanism{Species, Reactions, EqSteps}) where {Species, Reactions, EqSteps}
    subs, prods, regs, enzs = Species
    enz_names = Set(e[1] for e in enzs)

    # Check if mechanism is linear (each enzyme form appears on LHS and RHS at most once)
    lhs_counts = Dict{Symbol,Int}()
    rhs_counts = Dict{Symbol,Int}()
    for (lhs, rhs) in Reactions
        for s in lhs; s in enz_names && (lhs_counts[s] = get(lhs_counts, s, 0) + 1); end
        for s in rhs; s in enz_names && (rhs_counts[s] = get(rhs_counts, s, 0) + 1); end
    end
    is_linear = all(v <= 1 for v in values(lhs_counts)) && all(v <= 1 for v in values(rhs_counts))

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
end

# ─── Accessors ─────────────────────────────────────────────────────────────────

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

"""Return unique metabolites as a tuple of (name, atoms)."""
@generated function metabolites(::EnzymeMechanism{Species}) where {Species}
    return Tuple(_unique_metabolites(Species))
end

"""Number of distinct enzyme states."""
n_states(::EnzymeMechanism{Species}) where {Species} = length(Species[4])

"""Number of steps in the mechanism."""
n_steps(::EnzymeMechanism{Species, Reactions}) where {Species, Reactions} = length(Reactions)

"""Return the reactions tuple directly."""
reactions(::EnzymeMechanism{Species, R}) where {Species, R} = R

"""Return the equilibrium steps tuple (true = rapid-equilibrium, false = steady-state)."""
equilibrium_steps(::EnzymeMechanism{Sp, Rx, Eq}) where {Sp, Rx, Eq} = Eq

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
@generated function stoich_matrix(::EnzymeMechanism{Species, Reactions, EqSteps}) where {Species, Reactions, EqSteps}
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
