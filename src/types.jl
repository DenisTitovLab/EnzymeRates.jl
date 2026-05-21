using LinearAlgebra: rank

# ─── Concrete type hierarchy (spec §5.1–5.7) ──────────────────────────

# §5.1 — Metabolite / Reactant / Regulator hierarchy.
abstract type Metabolite end
abstract type Reactant <: Metabolite end
abstract type Regulator <: Metabolite end

struct Substrate <: Reactant
    name::Symbol
end
struct Product <: Reactant
    name::Symbol
end
struct AllostericRegulator <: Regulator
    name::Symbol
end
struct CompetitiveInhibitor <: Regulator
    name::Symbol
end

name(s::Substrate)            = s.name
name(p::Product)              = p.name
name(a::AllostericRegulator)  = a.name
name(c::CompetitiveInhibitor) = c.name

# §5.2 — Residual: substrates added + products subtracted from the enzyme
# (e.g., a covalent adduct after a ping-pong half-reaction). Empty
# `Residual()` means no covalent residue.
struct Residual
    added::Vector{Substrate}
    subtracted::Vector{Product}
    function Residual(added::Vector{Substrate}, subtracted::Vector{Product})
        new(sort(added; by=name), sort(subtracted; by=name))
    end
end
Residual() = Residual(Substrate[], Product[])

added(r::Residual)        = r.added
subtracted(r::Residual)   = r.subtracted
Base.isempty(r::Residual) = isempty(r.added) && isempty(r.subtracted)
Base.:(==)(a::Residual, b::Residual) =
    a.added == b.added && a.subtracted == b.subtracted
Base.hash(r::Residual, h::UInt) =
    hash(r.subtracted, hash(r.added, hash(:Residual, h)))

# §5.3 — Species: an enzyme form. `bound` is sorted by name; the
# rendered Symbol name reads `:E` / `:E_A_B` / `:Estar...` / `:E_A_B_res...`.
struct Species
    bound::Vector{Metabolite}
    conformation::Symbol
    residual::Residual
    function Species(bound::Vector{<:Metabolite}, conformation::Symbol,
                     residual::Residual)
        new(sort(Vector{Metabolite}(bound); by=name), conformation, residual)
    end
end
Species(bound, conformation::Symbol) = Species(bound, conformation, Residual())

bound(s::Species)        = s.bound
conformation(s::Species) = s.conformation
residual(s::Species)     = s.residual
has_residual(s::Species) = !isempty(s.residual)
Base.:(==)(a::Species, b::Species) =
    a.conformation == b.conformation && a.residual == b.residual &&
    a.bound == b.bound
Base.hash(s::Species, h::UInt) =
    hash(s.bound, hash(s.conformation, hash(s.residual, hash(:Species, h))))

# Render species name deterministically from fields:
#   :<conformation>[_<bound1>_<bound2>...][_res[_+<added>...][_-<subtracted>...]]
# Underscore is the field separator, so metabolite Symbols containing `_`
# produce ambiguous output (Species([Substrate(:A_B)]) and
# Species([Substrate(:A), Substrate(:B)]) both render :E_A_B). Domain
# convention: metabolite Symbols must not contain `_`.
function name(s::Species)
    parts = String[String(s.conformation)]
    for m in s.bound
        push!(parts, String(name(m)))
    end
    if has_residual(s)
        push!(parts, "res")
        for a in added(s.residual)
            push!(parts, "+" * String(name(a)))
        end
        for r in subtracted(s.residual)
            push!(parts, "-" * String(name(r)))
        end
    end
    Symbol(join(parts, "_"))
end

# §5.5 — RegulatorySite: a binding site (possibly multimeric) for one or
# more allosteric ligands. `ligands[i]` and `allo_states[i]` are parallel;
# ordering is meaningful (canonicalize at the call site if needed).
struct RegulatorySite
    ligands::Vector{AllostericRegulator}
    multiplicity::Int
    allo_states::Vector{Symbol}
    function RegulatorySite(ligands::Vector{AllostericRegulator},
                            multiplicity::Int, allo_states::Vector{Symbol})
        length(ligands) == length(allo_states) ||
            error("RegulatorySite: length(ligands)=$(length(ligands)) " *
                  "must equal length(allo_states)=$(length(allo_states))")
        multiplicity ≥ 1 ||
            error("RegulatorySite: multiplicity must be ≥ 1, got $multiplicity")
        for st in allo_states
            st in (:OnlyR, :OnlyT, :EqualRT, :NonequalRT) ||
                error("RegulatorySite: allo state $st must be one of " *
                      ":OnlyR, :OnlyT, :EqualRT, :NonequalRT")
        end
        new(ligands, multiplicity, allo_states)
    end
end

ligands(s::RegulatorySite)      = s.ligands
multiplicity(s::RegulatorySite) = s.multiplicity
allo_states(s::RegulatorySite)  = s.allo_states
Base.:(==)(a::RegulatorySite, b::RegulatorySite) =
    a.ligands == b.ligands && a.multiplicity == b.multiplicity &&
    a.allo_states == b.allo_states
Base.hash(s::RegulatorySite, h::UInt) =
    hash(s.allo_states, hash(s.multiplicity,
        hash(s.ligands, hash(:RegulatorySite, h))))

# §5.4 — Step: one elementary transition. Binding steps carry
# `bound_metabolite`; iso steps carry `nothing`. Binding steps are
# canonicalized to have the metabolite on the `from_species` side. Iso
# steps are canonicalized by lexical species-name order.
struct Step
    from_species::Species
    to_species::Species
    bound_metabolite::Union{Metabolite,Nothing}
    is_equilibrium::Bool
    function Step(from_species::Species, to_species::Species,
                  bound_metabolite::Union{Metabolite,Nothing},
                  is_equilibrium::Bool)
        if bound_metabolite !== nothing
            # Canonical binding direction: metabolite is BOUND in
            # to_species (and free, i.e., not in bound(from_species)),
            # matching the existing "E + S ⇌ E_S" convention.
            in_from = any(m -> m == bound_metabolite, bound(from_species))
            in_to   = any(m -> m == bound_metabolite, bound(to_species))
            if in_from && !in_to
                from_species, to_species = to_species, from_species
            end
        else
            # Iso steps: deterministic direction by lex on name(from_species).
            # Stronger ordering (substrate-then-product-content tiebreak) may
            # be needed once consumers exist; today's invariant is determinism
            # only.
            if string(name(from_species)) > string(name(to_species))
                from_species, to_species = to_species, from_species
            end
        end
        new(from_species, to_species, bound_metabolite, is_equilibrium)
    end
end

from_species(s::Step)     = s.from_species
to_species(s::Step)       = s.to_species
bound_metabolite(s::Step) = s.bound_metabolite
is_equilibrium(s::Step)   = s.is_equilibrium
is_binding(s::Step)       = s.bound_metabolite !== nothing
is_iso(s::Step)           = s.bound_metabolite === nothing
direction(s::Step)        = is_binding(s) ? :binding : :iso
Base.:(==)(a::Step, b::Step) =
    a.from_species == b.from_species && a.to_species == b.to_species &&
    a.bound_metabolite == b.bound_metabolite &&
    a.is_equilibrium == b.is_equilibrium
Base.hash(s::Step, h::UInt) =
    hash(s.is_equilibrium, hash(s.bound_metabolite,
        hash(s.to_species, hash(s.from_species, hash(:Step, h)))))

# §5.6 — Parameter family.
abstract type Parameter end

# Step-bound RE parameters
struct Kd   <: Parameter; step::Step; state::Symbol end
struct Kiso <: Parameter; step::Step; state::Symbol end

# Step-bound SS parameters
struct Kon  <: Parameter; step::Step; state::Symbol end
struct Koff <: Parameter; step::Step; state::Symbol end
struct Kfor <: Parameter; step::Step; state::Symbol end
struct Krev <: Parameter; step::Step; state::Symbol end

# Regulator-site parameter: a single ligand at a single site can appear
# in either the R or T branch of the polynomial.
struct Kreg <: Parameter
    site::RegulatorySite
    ligand::AllostericRegulator
    state::Symbol
end

# Mechanism-level scalars (singletons)
struct Keq   <: Parameter end
struct Etot  <: Parameter end
struct Lallo <: Parameter end

# Step-bound governance: only step-bound subtypes have a step. Kreg /
# Keq / Etot / Lallo intentionally have no `governing_step` method.
const StepBoundParameter = Union{Kd, Kiso, Kon, Koff, Kfor, Krev}
const StatefulParameter  = Union{Kd, Kiso, Kon, Koff, Kfor, Krev, Kreg}

governing_step(p::StepBoundParameter) = p.step
is_t_state(p::StatefulParameter)      = p.state === :T

for T in (:Kd, :Kiso, :Kon, :Koff, :Kfor, :Krev)
    @eval Base.:(==)(a::$T, b::$T) =
        a.step == b.step && a.state === b.state
    @eval Base.hash(p::$T, h::UInt) =
        hash(p.state, hash(p.step, hash($(QuoteNode(T)), h)))
end
Base.:(==)(a::Kreg, b::Kreg) =
    a.site == b.site && a.ligand == b.ligand && a.state === b.state
Base.hash(p::Kreg, h::UInt) =
    hash(p.state, hash(p.ligand, hash(p.site, hash(:Kreg, h))))

# §5.7 — Per-reactant and per-regulator bundling structs. Canonical
# ordering of atoms / multiplicities so two equivalent constructions
# compare equal under `==` / `hash`.
struct ReactantAtoms
    metabolite::Reactant
    atoms::Vector{Pair{Symbol,Int}}
    function ReactantAtoms(metabolite::Reactant,
                           atoms::Vector{<:Pair{Symbol,<:Integer}})
        new(metabolite, sort(Vector{Pair{Symbol,Int}}(atoms); by=first))
    end
end

metabolite(r::ReactantAtoms) = r.metabolite
atoms(r::ReactantAtoms)      = r.atoms
Base.:(==)(a::ReactantAtoms, b::ReactantAtoms) =
    a.metabolite == b.metabolite && a.atoms == b.atoms
Base.hash(r::ReactantAtoms, h::UInt) =
    hash(r.atoms, hash(r.metabolite, hash(:ReactantAtoms, h)))

struct RegulatorMults
    regulator::Regulator
    allowed_multiplicities::Vector{Int}
    function RegulatorMults(regulator::Regulator,
                            allowed_multiplicities::Vector{Int})
        all(m -> m ≥ 1, allowed_multiplicities) ||
            error("RegulatorMults: allowed_multiplicities must all be ≥ 1, " *
                  "got $allowed_multiplicities")
        new(regulator, sort(allowed_multiplicities))
    end
end

regulator(r::RegulatorMults)              = r.regulator
allowed_multiplicities(r::RegulatorMults) = r.allowed_multiplicities
Base.:(==)(a::RegulatorMults, b::RegulatorMults) =
    a.regulator == b.regulator &&
    a.allowed_multiplicities == b.allowed_multiplicities
Base.hash(r::RegulatorMults, h::UInt) =
    hash(r.allowed_multiplicities,
         hash(r.regulator, hash(:RegulatorMults, h)))

# §5.7 — EnzymeReaction: the public concrete reaction descriptor. Holds
# reactants (substrate + product atom payload), regulators (with allowed
# multiplicity sets), and the catalytic multiplicities the enumerator is
# allowed to consider. Canonical ordering is enforced so two equivalent
# constructions compare equal under `==` / `hash`.
struct EnzymeReaction
    reactants::Vector{ReactantAtoms}
    regulators::Vector{RegulatorMults}
    allowed_catalytic_multiplicities::Vector{Int}

    function EnzymeReaction(reactants::Vector{ReactantAtoms},
                            regulators::Vector{RegulatorMults},
                            allowed_catalytic_multiplicities::Vector{Int})
        sorted_reactants = sort(reactants; by = ra -> name(metabolite(ra)))
        sorted_regulators = sort(regulators; by = rm -> name(regulator(rm)))
        sorted_mults = sort(unique(allowed_catalytic_multiplicities))
        all(m -> m ≥ 1, sorted_mults) ||
            error("EnzymeReaction: allowed_catalytic_multiplicities must " *
                  "all be ≥ 1, got $allowed_catalytic_multiplicities")
        new(sorted_reactants, sorted_regulators, sorted_mults)
    end
end

reactants(r::EnzymeReaction) = r.reactants
regulators(r::EnzymeReaction) = r.regulators
allowed_catalytic_multiplicities(r::EnzymeReaction) =
    r.allowed_catalytic_multiplicities
substrates(r::EnzymeReaction) =
    Substrate[metabolite(ra) for ra in r.reactants if metabolite(ra) isa Substrate]
products(r::EnzymeReaction) =
    Product[metabolite(ra) for ra in r.reactants if metabolite(ra) isa Product]

Base.:(==)(a::EnzymeReaction, b::EnzymeReaction) =
    a.reactants == b.reactants && a.regulators == b.regulators &&
    a.allowed_catalytic_multiplicities == b.allowed_catalytic_multiplicities
Base.hash(r::EnzymeReaction, h::UInt) =
    hash(r.allowed_catalytic_multiplicities,
         hash(r.regulators,
              hash(r.reactants, hash(:EnzymeReaction, h))))

function Base.show(io::IO, r::EnzymeReaction)
    subs_str  = join(String.(name.(substrates(r))), " + ")
    prods_str = join(String.(name.(products(r))),   " + ")
    print(io, "EnzymeReaction: ", subs_str, " ⇌ ", prods_str)
    if !isempty(r.regulators)
        regs_str = join(
            (String(name(regulator(rm))) for rm in r.regulators), ", ")
        print(io, " | regulators: ", regs_str)
    end
    mults = r.allowed_catalytic_multiplicities
    if length(mults) == 1 && mults[1] > 1
        print(io, " | oligomeric_state: ", mults[1])
    elseif length(mults) > 1 || (length(mults) == 1 && mults[1] != 1)
        print(io, " | allowed_catalytic_multiplicities: (",
              join(mults, ", "), ")")
    end
end

# §5.8 — Mechanism: groups elementary steps by kinetic group (outer
# vector). All steps within a group share kinetic parameters.
struct Mechanism
    reaction::EnzymeReaction
    steps::Vector{Vector{Step}}
end

reaction(m::Mechanism) = m.reaction
steps(m::Mechanism) = m.steps
kinetic_groups(m::Mechanism) = 1:length(m.steps)
n_steps(m::Mechanism) = sum(length, m.steps; init = 0)
rep_step(m::Mechanism, g::Int) = first(m.steps[g])

Base.:(==)(a::Mechanism, b::Mechanism) =
    a.reaction == b.reaction && a.steps == b.steps
Base.hash(m::Mechanism, h::UInt) =
    hash(m.steps, hash(m.reaction, hash(:Mechanism, h)))

# §5.8 — AllostericMechanism: a multi-subunit MWC enzyme. Each catalytic
# kinetic group carries an allosteric-state tag (`:OnlyR`, `:EqualRT`, or
# `:NonequalRT` — `:OnlyT` is rejected by the R-state-active convention).
const _VALID_CAT_ALLO_STATES = (:OnlyR, :EqualRT, :NonequalRT)

struct AllostericMechanism
    reaction::EnzymeReaction
    cat_steps::Vector{Vector{Step}}
    cat_allo_states::Vector{Symbol}
    catalytic_multiplicity::Int
    regulatory_sites::Vector{RegulatorySite}

    function AllostericMechanism(reaction::EnzymeReaction,
                                 cat_steps::Vector{Vector{Step}},
                                 cat_allo_states::Vector{Symbol},
                                 catalytic_multiplicity::Int,
                                 regulatory_sites::Vector{RegulatorySite})
        length(cat_allo_states) == length(cat_steps) ||
            error("AllostericMechanism: cat_allo_states length " *
                  "$(length(cat_allo_states)) must match cat_steps length " *
                  "$(length(cat_steps))")
        catalytic_multiplicity ≥ 1 ||
            error("AllostericMechanism: catalytic_multiplicity must be ≥ 1, " *
                  "got $catalytic_multiplicity")
        for (g, tag) in enumerate(cat_allo_states)
            tag in _VALID_CAT_ALLO_STATES ||
                error("AllostericMechanism: catalytic group $g has invalid " *
                      "allo state $tag (must be one of " *
                      "$_VALID_CAT_ALLO_STATES); :OnlyT is rejected for " *
                      "catalytic groups (R-state-active convention)")
        end
        new(reaction, cat_steps, cat_allo_states,
            catalytic_multiplicity, regulatory_sites)
    end
end

reaction(m::AllostericMechanism) = m.reaction
steps(m::AllostericMechanism) = m.cat_steps
cat_allo_state(m::AllostericMechanism, g::Int) = m.cat_allo_states[g]
catalytic_multiplicity(m::AllostericMechanism) = m.catalytic_multiplicity
regulatory_sites(m::AllostericMechanism) = m.regulatory_sites
kinetic_groups(m::AllostericMechanism) = 1:length(m.cat_steps)
n_steps(m::AllostericMechanism) = sum(length, m.cat_steps; init = 0)
rep_step(m::AllostericMechanism, g::Int) = first(m.cat_steps[g])

function allosteric_regulators(m::AllostericMechanism)
    seen = AllostericRegulator[]
    for site in m.regulatory_sites, lig in site.ligands
        lig in seen || push!(seen, lig)
    end
    seen
end

function competitive_inhibitors(m::AllostericMechanism)
    CompetitiveInhibitor[regulator(rm) for rm in m.reaction.regulators
                         if regulator(rm) isa CompetitiveInhibitor]
end

Base.:(==)(a::AllostericMechanism, b::AllostericMechanism) =
    a.reaction == b.reaction && a.cat_steps == b.cat_steps &&
    a.cat_allo_states == b.cat_allo_states &&
    a.catalytic_multiplicity == b.catalytic_multiplicity &&
    a.regulatory_sites == b.regulatory_sites
Base.hash(m::AllostericMechanism, h::UInt) =
    hash(m.regulatory_sites,
         hash(m.catalytic_multiplicity,
              hash(m.cat_allo_states,
                   hash(m.cat_steps,
                        hash(m.reaction,
                             hash(:AllostericMechanism, h))))))

# ─── Mechanism ↔ Sig (parametric ↔ non-parametric) conversion ──
#
# Every leaf in `sig` MUST be a valid Julia type-parameter value (isbits,
# Symbol, type, or Tuple of those). `Pair{Symbol,Int}` is NOT valid as a
# type parameter — encode pairs as `Tuple{Symbol,Int}`. Vectors are
# NEVER valid — always wrap in `Tuple(...)`.
#
# One polymorphic `_to_sig` with a method per source type; the matching
# `_*_from_sig` family reconstructs the corresponding type.

_to_sig(s::Substrate)            = (:Substrate, s.name)
_to_sig(p::Product)              = (:Product, p.name)
_to_sig(r::AllostericRegulator)  = (:AllostericRegulator, r.name)
_to_sig(c::CompetitiveInhibitor) = (:CompetitiveInhibitor, c.name)

_to_sig(r::Residual) = (
    Tuple(_to_sig(m) for m in r.added),
    Tuple(_to_sig(m) for m in r.subtracted),
)

_to_sig(s::Species) = (
    Tuple(_to_sig(m) for m in s.bound),
    s.conformation,
    _to_sig(s.residual),
)

_to_sig(s::Step) = (
    _to_sig(s.from_species),
    _to_sig(s.to_species),
    s.bound_metabolite === nothing ? nothing : _to_sig(s.bound_metabolite),
    s.is_equilibrium,
)

_to_sig(ra::ReactantAtoms) = (
    _to_sig(ra.metabolite),
    Tuple((p.first, p.second) for p in ra.atoms),   # Tuple{Symbol,Int}, NOT Pair
)

_to_sig(rm::RegulatorMults) = (
    _to_sig(rm.regulator),
    Tuple(rm.allowed_multiplicities),
)

_to_sig(r::EnzymeReaction) = (
    Tuple(_to_sig(ra) for ra in r.reactants),
    Tuple(_to_sig(rm) for rm in r.regulators),
    Tuple(r.allowed_catalytic_multiplicities),
)

function _metabolite_from_sig(sig::Tuple{Symbol, Symbol})
    kind, nm = sig
    kind === :Substrate            ? Substrate(nm)            :
    kind === :Product              ? Product(nm)              :
    kind === :AllostericRegulator  ? AllostericRegulator(nm)  :
    kind === :CompetitiveInhibitor ? CompetitiveInhibitor(nm) :
    error("Unknown metabolite kind in sig: $kind")
end

function _residual_from_sig(sig::Tuple)
    added_sig, sub_sig = sig
    Residual(
        Substrate[_metabolite_from_sig(t) for t in added_sig],
        Product[_metabolite_from_sig(t)   for t in sub_sig],
    )
end

function _species_from_sig(sig::Tuple)
    bound_sig, conformation, residual_sig = sig
    Species(
        Metabolite[_metabolite_from_sig(t) for t in bound_sig],
        conformation,
        _residual_from_sig(residual_sig),
    )
end

function _step_from_sig(sig::Tuple)
    from_sig, to_sig, met_sig, is_eq = sig
    met = met_sig === nothing ? nothing : _metabolite_from_sig(met_sig)
    Step(_species_from_sig(from_sig), _species_from_sig(to_sig), met, is_eq)
end

function _reactant_atoms_from_sig(sig::Tuple)
    met_sig, atoms_sig = sig
    ReactantAtoms(
        _metabolite_from_sig(met_sig)::Reactant,
        Pair{Symbol,Int}[s => c for (s, c) in atoms_sig],
    )
end

function _regulator_mults_from_sig(sig::Tuple)
    reg_sig, mults_sig = sig
    RegulatorMults(
        _metabolite_from_sig(reg_sig)::Regulator,
        Int[m for m in mults_sig],
    )
end

function _reaction_from_sig(sig::Tuple)
    reactants_sig, regulators_sig, mults_sig = sig
    EnzymeReaction(
        ReactantAtoms[_reactant_atoms_from_sig(t) for t in reactants_sig],
        RegulatorMults[_regulator_mults_from_sig(t) for t in regulators_sig],
        Int[m for m in mults_sig],
    )
end

function _steps_from_sig(sig::Tuple)
    Vector{Step}[Step[_step_from_sig(s) for s in group] for group in sig]
end

_sig_of(m::Mechanism) = (
    _to_sig(m.reaction),
    Tuple(Tuple(_to_sig(s) for s in g) for g in m.steps),
)

function _mechanism_from_sig(sig::Tuple)
    reaction_sig, steps_sig = sig
    Mechanism(_reaction_from_sig(reaction_sig), _steps_from_sig(steps_sig))
end

# ─── Parametric mechanism types ───────────────────────────────────────

"""Sort species tuples alphabetically by name (first element)."""
_sort_species(t::Tuple) = Tuple(sort(collect(t); by=s -> s[1]))

"""
    EnzymeReactionLegacy{Substrates, Products, Regulators, OligomericState}

Singleton type encoding an enzyme reaction specification in type parameters.

- `Substrates`, `Products`: tuple of `(name::Symbol, atoms::Tuple{Vararg{Tuple{Symbol,Int}}})`.
- `Regulators`: tuple of `Symbol` (plain names, no atoms).
- `OligomericState`: number of subunits (Int, default 1).
"""
struct EnzymeReactionLegacy{Substrates, Products, Regulators, OligomericState} end

"""Sum element counts across a tuple of `(name, atoms)` pairs.
Returns a Dict{Symbol,Int}. Errors if any species's atoms tuple
is empty (atoms are mandatory) or if any per-atom count is not a
positive Int."""
function _sum_atoms(species::Tuple, side::String)
    totals = Dict{Symbol,Int}()
    for (name, atoms) in species
        isempty(atoms) && error(
            "EnzymeReactionLegacy: $side metabolite $name has no declared " *
            "atoms; atoms are mandatory (use `[C…]` bracket syntax in " *
            "@enzyme_reaction or pass non-empty atom tuples to the " *
            "constructor).")
        for (elem, count) in atoms
            count isa Integer && !(count isa Bool) && count > 0 ||
                error(
                "EnzymeReactionLegacy: $side metabolite $name has " *
                "non-positive atom count for element $elem ($count); " *
                "atom counts must be positive integers (not Bool).")
            totals[elem] = get(totals, elem, 0) + count
        end
    end
    totals
end

function EnzymeReactionLegacy(subs::Tuple, prods::Tuple, regs::Tuple=(); oligomeric_state::Int=1)
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
            "EnzymeReactionLegacy: atom imbalance — element $elem appears " *
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
    EnzymeReactionLegacy{subs, prods, sorted_regs, oligomeric_state}()
end

# Convert the concrete EnzymeReaction into its parametric Legacy form.
# Used to feed EnzymeReaction values into call sites that still
# dispatch on EnzymeReactionLegacy (init_mechanisms, expand_mechanisms,
# IdentifyRateEquationProblem). Stage 5 / Stage 6 retire the consumers
# of the Legacy form, after which this helper becomes dead code.
function _to_legacy_reaction(r::EnzymeReaction)
    subs = Tuple((name(metabolite(ra)),
                  Tuple((e => c) for (e, c) in atoms(ra)))
                 for ra in r.reactants if metabolite(ra) isa Substrate)
    prods = Tuple((name(metabolite(ra)),
                   Tuple((e => c) for (e, c) in atoms(ra)))
                  for ra in r.reactants if metabolite(ra) isa Product)
    # Legacy atoms shape is Tuple{Tuple{Symbol,Int}...} — convert
    # the Pair tuples above to (Symbol, Int) tuples.
    subs_t = Tuple((nm, Tuple((p.first, p.second) for p in at))
                   for (nm, at) in subs)
    prods_t = Tuple((nm, Tuple((p.first, p.second) for p in at))
                    for (nm, at) in prods)
    regs_t = Tuple((name(regulator(rm)),
                    regulator(rm) isa AllostericRegulator ?
                        :allosteric : :dead_end)
                   for rm in r.regulators)
    mults = r.allowed_catalytic_multiplicities
    length(mults) == 1 ||
        error("_to_legacy_reaction: EnzymeReactionLegacy supports a single " *
              "oligomeric_state; got allowed_catalytic_multiplicities=$mults. " *
              "Enumeration over multiple multiplicities is a Stage 5 capability.")
    EnzymeReactionLegacy(subs_t, prods_t, regs_t; oligomeric_state=mults[1])
end

abstract type AbstractEnzymeMechanism end

"""
    EnzymeMechanism{Sig}

Singleton type encoding an enzyme mechanism. `Sig` is the tuple
`(reaction_sig, steps_sig)` produced by `_sig_of(::Mechanism)`:

- `reaction_sig` encodes the `EnzymeReaction` (reactants, regulators,
  catalytic multiplicities) as a tuple of tuples of `Symbol`s/`Int`s.
- `steps_sig` encodes `Vector{Vector{Step}}` as a nested tuple where the
  outer level is kinetic groups and the inner level is `Step` data
  (from-species, to-species, bound metabolite, is-equilibrium).

The conversion functions are `_sig_of` and `_mechanism_from_sig`; the
exact layout is internal — users construct via `EnzymeMechanism(::Mechanism)`
or the `EnzymeMechanism(metabolites, reactions)` shorthand below.
"""
struct EnzymeMechanism{Sig} <: AbstractEnzymeMechanism end

# Boundary converters: non-parametric Mechanism ↔ parametric
# EnzymeMechanism{Sig}. The new-shape Sig encodes the Mechanism's data
# as `(reaction_sig, steps_sig)` produced by `_sig_of`. The legacy-shape
# Sig `(metabolites_3tuple, rxns_4tuple)` is produced by
# `EnzymeMechanism(metabolites, reactions)` and the enumerator's
# `EnzymeMechanism(spec::MechanismSpec)`; `Mechanism(em)` lifts either
# shape so Stage 3 derivation consumers can walk a `Mechanism` uniformly.
EnzymeMechanism(m::Mechanism) = EnzymeMechanism{_sig_of(m)}()

function Mechanism(em::EnzymeMechanism{Sig}) where {Sig}
    _is_new_sig(Sig) && return _mechanism_from_sig(Sig)
    _mechanism_from_legacy_sig(Sig)
end

"""
Build a `Mechanism` from the legacy `(metabolites_3tuple, rxns_4tuple)`
Sig shape. Enzyme-form Symbols become conformation-only `Species` (empty
`bound`, empty `Residual`) so `name(species)` round-trips to the
original Symbol. The bound-metabolite side info is encoded by ordering:
binding direction `from + met → to` places the metabolite in the
to-species's bound list; release direction `from → to + met` places it
in the from-species's bound list. For opaque symbol forms (where the
mechanism-level Symbol does not encode bound metabolites), the metabolite
goes in the bound list of the side opposite the free-met side in the
original tuple, preserving direction info that the `_legacy_step_tuple`
accessor needs to reconstruct the user's source-order shape.
"""
function _mechanism_from_legacy_sig(Sig::Tuple)
    mets_sig, rxns_sig = Sig
    subs, prods, regs = mets_sig

    reactants = ReactantAtoms[]
    for s in subs
        push!(reactants,
            ReactantAtoms(Substrate(s), Pair{Symbol,Int}[:C => 1]))
    end
    for p in prods
        push!(reactants,
            ReactantAtoms(Product(p), Pair{Symbol,Int}[:C => 1]))
    end
    regulators_vec = RegulatorMults[]
    for r in regs
        push!(regulators_vec,
            RegulatorMults(CompetitiveInhibitor(r), Int[1]))
    end
    reaction = EnzymeReaction(reactants, regulators_vec, Int[1])

    sub_set, prod_set, reg_set =
        Set{Symbol}(subs), Set{Symbol}(prods), Set{Symbol}(regs)
    met_set = sub_set ∪ prod_set ∪ reg_set
    metabolite_of(s::Symbol) =
        s in sub_set  ? Substrate(s) :
        s in prod_set ? Product(s) :
        s in reg_set  ? CompetitiveInhibitor(s) :
        error("Symbol $s is not a declared metabolite or regulator")

    groups = Dict{Int, Vector{Step}}()
    group_order = Int[]
    for (lhs, rhs, is_eq, gnum) in rxns_sig
        e_lhs = first(s for s in lhs if s ∉ met_set)
        e_rhs = first(s for s in rhs if s ∉ met_set)
        m_lhs = Symbol[s for s in lhs if s ∈ met_set]
        m_rhs = Symbol[s for s in rhs if s ∈ met_set]

        bound_met  = nothing
        from_bound = Metabolite[]
        to_bound   = Metabolite[]
        if !isempty(m_lhs)
            bound_met = metabolite_of(m_lhs[1])
            push!(to_bound, bound_met)
        elseif !isempty(m_rhs)
            bound_met = metabolite_of(m_rhs[1])
            push!(from_bound, bound_met)
        end

        from = Species(from_bound, e_lhs, Residual())
        to   = Species(to_bound,   e_rhs, Residual())
        step = Step(from, to, bound_met, is_eq)

        if !haskey(groups, gnum)
            groups[gnum] = Step[]
            push!(group_order, gnum)
        end
        push!(groups[gnum], step)
    end
    step_groups = Vector{Step}[groups[g] for g in group_order]
    Mechanism(reaction, step_groups)
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

    # Build the singleton type and run remaining checks via accessors.
    # The legacy 2-tuple constructor keeps the legacy `(metabolites, rxns)`
    # Sig shape; accessors below dispatch on shape so both this path and
    # the new `EnzymeMechanism(::Mechanism)` constructor work uniformly.
    m = EnzymeMechanism{(mets, rxns)}()

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
# The three type parameters cannot be folded into a single value-tuple `Sig`
# (as EnzymeMechanism{Sig} does): the first slot is a DataType (an
# EnzymeMechanism subtype), and Julia rejects DataTypes inside the
# value-tuple position of a type parameter.
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

function Base.show(io::IO, ::EnzymeReactionLegacy{S,P,R,N}) where {S,P,R,N}
    subs_str = join([string(name) for (name, _) in S], " + ")
    prods_str = join([string(name) for (name, _) in P], " + ")
    print(io, "EnzymeReaction: ", subs_str, " ⇌ ", prods_str)
    if !isempty(R)
        regs_str = join([string(r[1]) for r in R], ", ")
        print(io, " | regulators: ", regs_str)
    end
    N > 1 && print(io, " | oligomeric_state: ", N)
end

function Base.show(io::IO, m::EnzymeMechanism)
    # `show` reads through the accessors so it works for both Sig shapes.
    Rxns = reactions(m)
    regs = regulators(m)
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
#
# Accessors on `EnzymeMechanism{Sig}` dispatch on Sig content to
# preserve compatibility with two construction paths during the
# concrete-types transition:
#
# - Legacy Sig `(metabolites_3tuple, rxns_4tuple)` — produced by
#   `EnzymeMechanism(metabolites, reactions)` and the enumeration
#   internals' `EnzymeMechanism(spec::MechanismSpec)`. The first slot
#   is a 3-tuple of name `Tuple{Symbol,...}`s.
# - New Sig `(reaction_sig, steps_sig)` — produced by
#   `EnzymeMechanism(::Mechanism)` (via `_sig_of`). The first slot is
#   the EnzymeReaction encoding (reactants + regulators + mults), so
#   `Sig[1]` is also a 3-tuple but its first element is a tuple of
#   `((kind::Symbol, name::Symbol), atoms_tuple)` reactants, not a
#   tuple of names.
#
# `_is_new_sig` is the discriminator: it checks whether `Sig[1][1]`
# (the first reactant slot) is a 2-tuple with a `Symbol`-`Symbol` head
# (a metabolite encoding from `_to_sig`) versus a `Tuple{Symbol,...}`
# (a name tuple from the legacy shape).

"""True iff `Sig` is the new shape `(reaction_sig, steps_sig)` produced
by `_sig_of`. False iff it is the legacy `(mets, rxns)` shape produced
by `EnzymeMechanism(metabolites, reactions)`."""
function _is_new_sig(Sig::Tuple)
    length(Sig) == 2 || return false
    reaction_or_mets = Sig[1]
    length(reaction_or_mets) == 3 || return false
    first_slot = reaction_or_mets[1]
    first_slot isa Tuple || return false
    # New-shape reactants: a tuple of (metabolite_sig, atoms_sig) pairs,
    # where each metabolite_sig is `(:Substrate|:Product|..., name)`.
    # Legacy mets: substrates as a `Tuple{Symbol,...}` of bare names.
    isempty(first_slot) && return false
    entry = first_slot[1]
    entry isa Tuple && length(entry) == 2 &&
        entry[1] isa Tuple && length(entry[1]) == 2 &&
        entry[1][1] isa Symbol &&
        entry[1][1] in (:Substrate, :Product,
                        :AllostericRegulator, :CompetitiveInhibitor)
end

"""Render a `Species` as the bare Symbol used by the legacy step shape."""
_species_sym(s::Species) = name(s)

"""Walk the steps of `Mechanism(em)` in flat order, yielding
`(step::Step, kinetic_group::Int)` pairs."""
function _flat_steps(m::Mechanism)
    out = Tuple{Step, Int}[]
    for (g, group) in enumerate(m.steps)
        for s in group
            push!(out, (s, g))
        end
    end
    out
end

"""Build a (lhs_syms, rhs_syms, is_eq, kinetic_group) tuple from a `Step`.

The legacy 4-tuple shape pairs the free metabolite with the enzyme form
on the OPPOSITE side from where it is bound: `E + S ⇌ E_S` becomes
`((:E, :S), (:E_S,), true)`. We read `bound(to_species)` to decide
which side carries the free metabolite — if the metabolite is bound
in the to-species (binding direction), the free metabolite goes on
the from side; if bound in the from-species (release direction), the
free metabolite goes on the to side. Weak Species reps (empty `bound`
on both sides) default to the canonical binding direction
(free metabolite on the from side)."""
function _legacy_step_tuple(step::Step, g::Int)
    e_from = _species_sym(from_species(step))
    e_to   = _species_sym(to_species(step))
    met    = bound_metabolite(step)
    is_eq  = is_equilibrium(step)
    if met === nothing
        # Iso step — no metabolite on either side.
        return ((e_from,), (e_to,), is_eq, g)
    end
    met_name = name(met)
    bound_in_to = any(m -> m == met, bound(to_species(step)))
    bound_in_from = any(m -> m == met, bound(from_species(step)))
    if bound_in_to
        return ((e_from, met_name), (e_to,), is_eq, g)
    elseif bound_in_from
        return ((e_from,), (e_to, met_name), is_eq, g)
    end
    ((e_from, met_name), (e_to,), is_eq, g)
end

"""Return substrates as a tuple of `Symbol` names."""
function substrates(em::EnzymeMechanism{Sig}) where {Sig}
    if _is_new_sig(Sig)
        Tuple(name(metabolite(ra)) for ra in reactants(reaction(Mechanism(em)))
              if metabolite(ra) isa Substrate)
    else
        Sig[1][1]
    end
end
substrates(::EnzymeReactionLegacy{S,P,R,N}) where {S,P,R,N} = S

"""Return products as a tuple of `Symbol` names."""
function products(em::EnzymeMechanism{Sig}) where {Sig}
    if _is_new_sig(Sig)
        Tuple(name(metabolite(ra)) for ra in reactants(reaction(Mechanism(em)))
              if metabolite(ra) isa Product)
    else
        Sig[1][2]
    end
end
products(::EnzymeReactionLegacy{S,P,R,N}) where {S,P,R,N} = P

"""Return regulators as a tuple of `Symbol` names."""
function regulators(em::EnzymeMechanism{Sig}) where {Sig}
    if _is_new_sig(Sig)
        Tuple(name(regulator(rm)) for rm in regulators(reaction(Mechanism(em))))
    else
        Sig[1][3]
    end
end
regulators(::EnzymeReactionLegacy{S,P,R,N}) where {S,P,R,N} =
    Tuple(r[1] for r in R)

"""Return regulator (name, role) pairs."""
regulator_roles(::EnzymeReactionLegacy{S,P,R,N}) where {S,P,R,N} = R

"""Return oligomeric state (number of subunits)."""
oligomeric_state(::EnzymeReactionLegacy{S,P,R,N}) where {S,P,R,N} = N

"""
    metabolites(m::EnzymeMechanism) → Tuple{Symbol,...}

Return distinct metabolite names (substrates ∪ products ∪ regulators) as a tuple
of `Symbol`s in declaration order, deduplicated.
"""
@generated function metabolites(::EnzymeMechanism{Sig}) where {Sig}
    if _is_new_sig(Sig)
        # New-shape: substrates first (from reactants Substrate entries),
        # then products, then regulators.
        names = Symbol[]
        seen = Set{Symbol}()
        for entry in Sig[1][1]
            kind, nm = entry[1]
            kind === :Substrate && nm ∉ seen &&
                (push!(seen, nm); push!(names, nm))
        end
        for entry in Sig[1][1]
            kind, nm = entry[1]
            kind === :Product && nm ∉ seen &&
                (push!(seen, nm); push!(names, nm))
        end
        for entry in Sig[1][2]
            kind, nm = entry[1]
            nm ∉ seen && (push!(seen, nm); push!(names, nm))
        end
        return Tuple(names)
    end
    M = Sig[1]
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
function reactions(em::EnzymeMechanism{Sig}) where {Sig}
    if _is_new_sig(Sig)
        mech = Mechanism(em)
        return Tuple(_legacy_step_tuple(s, g) for (s, g) in _flat_steps(mech))
    end
    Sig[2]
end

"""Return the equilibrium-step flags (`true` = rapid-equilibrium, `false` = steady-state)."""
@generated function equilibrium_steps(::EnzymeMechanism{Sig}) where {Sig}
    if _is_new_sig(Sig)
        return Tuple(step[4] for group in Sig[2] for step in group)
    end
    R = Sig[2]
    Tuple(step[3] for step in R)
end

"""Number of steps in the mechanism."""
function n_steps(::EnzymeMechanism{Sig}) where {Sig}
    _is_new_sig(Sig) ?
        sum(length(group) for group in Sig[2]; init=0) :
        length(Sig[2])
end

"""Kinetic group of step `idx`."""
function kinetic_group(em::EnzymeMechanism{Sig}, idx::Int) where {Sig}
    if _is_new_sig(Sig)
        flat = _flat_steps(Mechanism(em))
        1 ≤ idx ≤ length(flat) ||
            error("kinetic_group: step index $idx out of range 1:$(length(flat))")
        return flat[idx][2]
    end
    Sig[2][idx][4]
end

"""Sorted tuple of distinct kinetic group ids."""
@generated function kinetic_groups(::EnzymeMechanism{Sig}) where {Sig}
    if _is_new_sig(Sig)
        return Tuple(1:length(Sig[2]))
    end
    R = Sig[2]
    Tuple(sort(unique(step[4] for step in R)))
end

"""Indices of steps belonging to kinetic group `G`."""
@generated function steps_in_group(
    ::EnzymeMechanism{Sig}, ::Val{G},
) where {Sig, G}
    if _is_new_sig(Sig)
        idxs = Int[]
        flat = 0
        for (g, group) in enumerate(Sig[2])
            for _ in group
                flat += 1
                g == G && push!(idxs, flat)
            end
        end
        return Tuple(idxs)
    end
    R = Sig[2]
    Tuple(i for (i, step) in enumerate(R) if step[4] == G)
end
steps_in_group(m::EnzymeMechanism, g::Int) = steps_in_group(m, Val(g))

"""
    enzyme_forms(m::EnzymeMechanism) → Tuple{Symbol,...}

Return distinct enzyme-form names (any symbol appearing in a step that is not a
metabolite) as a tuple of `Symbol`s in step-order, deduplicated.
"""
function enzyme_forms(em::EnzymeMechanism{Sig}) where {Sig}
    if _is_new_sig(Sig)
        mech = Mechanism(em)
        met_names = Set(metabolites(em))
        seen = Set{Symbol}()
        forms = Symbol[]
        for (s, _) in _flat_steps(mech)
            for nm in (_species_sym(from_species(s)),
                       _species_sym(to_species(s)))
                nm ∉ met_names && nm ∉ seen &&
                    (push!(seen, nm); push!(forms, nm))
            end
        end
        return Tuple(forms)
    end
    _enzyme_forms_legacy(em)
end

@generated function _enzyme_forms_legacy(::EnzymeMechanism{Sig}) where {Sig}
    M, R = Sig
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
@generated function stoich_matrix(::EnzymeMechanism{Sig}) where {Sig}
    if _is_new_sig(Sig)
        # New shape: walk via the accessors at body-build time and
        # constant-fold the resulting matrix into the generated body.
        em_inst = EnzymeMechanism{Sig}()
        Rxns = reactions(em_inst)
        species = (enzyme_forms(em_inst)..., metabolites(em_inst)...)
        sp_idx = Dict(s => i for (i, s) in enumerate(species))
        S = zeros(Int, length(species), length(Rxns))
        for (j, (lhs, rhs, _, _)) in enumerate(Rxns)
            for s in lhs; S[sp_idx[s], j] -= 1; end
            for s in rhs; S[sp_idx[s], j] += 1; end
        end
        return S
    end
    M, R = Sig
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

# ─── name(p::Parameter, m) chokepoint ─────────────────────────────────
#
# Single chokepoint for parameter Symbol production. Renders positional
# names (:K1, :k1f, :K1_T, :K_<lig>_reg<i>) via Parameter subtype
# dispatch. Routing all parameter-name production through one function
# keeps any name-scheme change a single-function edit.

# Map a Step to the rep-idx used in positional parameter names. The rep
# is the position of the group's first step in the flattened step list,
# matching the CLAUDE.md parameter naming convention.
function _rep_idx_for_step(step::Step,
                           m::Union{Mechanism, AllostericMechanism})
    groups = m isa Mechanism ? m.steps : m.cat_steps
    pos = 0
    for group in groups
        if step in group
            return pos + 1
        end
        pos += length(group)
    end
    error("Step not found in mechanism: $step")
end

_rep_idx_for_step(step::Step, m::EnzymeMechanism) =
    _rep_idx_for_step(step, Mechanism(m))

function _site_idx_of(site::RegulatorySite, m::AllostericMechanism)
    idx = findfirst(==(site), m.regulatory_sites)
    idx === nothing && error("RegulatorySite not found in mechanism")
    return idx
end

# Step-bound parameters map to three rendering rules keyed by type:
# Kd/Kiso → :K{rep},  Kon/Kfor → :k{rep}f,  Koff/Krev → :k{rep}r,
# with optional `_T` suffix for T-state.
_step_param_prefix(::Union{Kd, Kiso})              = 'K'
_step_param_prefix(::Union{Kon, Koff, Kfor, Krev}) = 'k'
_step_param_suffix(::Union{Kd, Kiso})              = ""
_step_param_suffix(::Union{Kon, Kfor})             = "f"
_step_param_suffix(::Union{Koff, Krev})            = "r"

function name(p::StepBoundParameter,
              m::Union{Mechanism, EnzymeMechanism, AllostericMechanism})
    rep = _rep_idx_for_step(p.step, m)
    tag = p.state === :T ? "_T" : ""
    Symbol("$(_step_param_prefix(p))$rep$(_step_param_suffix(p))$tag")
end

# Regulator-site parameter — AllostericMechanism only.
function name(p::Kreg, m::AllostericMechanism)
    site_idx = _site_idx_of(p.site, m)
    lig_name = name(p.ligand)
    p.state === :T ? Symbol("K_$(lig_name)_T_reg$site_idx") :
                     Symbol("K_$(lig_name)_reg$site_idx")
end

# Mechanism-level scalars
name(::Keq,   _) = :Keq
name(::Etot,  _) = :E_total
name(::Lallo, _) = :L
