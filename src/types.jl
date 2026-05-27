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
        new(sort(Vector{Metabolite}(bound);
                 by = m -> (name(m), m isa CompetitiveInhibitor)),
            conformation, residual)
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
        push!(parts, m isa CompetitiveInhibitor ?
                     String(name(m)) * "inh" : String(name(m)))
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
#
# `source_idx` is presentation metadata only — the flat source-order
# position used to render positional parameter names (`:K1`, `:k10f`,
# ...). It is auto-assigned by the `Mechanism` constructor and read
# back by `_rep_idx_for_step`. It is NOT part of Step's structural
# identity (see `==` / `hash` below). A future refactor that moves to
# semantic parameter names (`:K_ATP`) will drop this field; all reads
# flow through `_rep_idx_for_step` so removal is a contained edit.
struct Step
    from_species::Species
    to_species::Species
    bound_metabolite::Union{Metabolite,Nothing}
    is_equilibrium::Bool
    source_idx::Int
    function Step(from_species::Species, to_species::Species,
                  bound_metabolite::Union{Metabolite,Nothing},
                  is_equilibrium::Bool;
                  source_idx::Int = 0)
        if bound_metabolite !== nothing
            # RE binding steps canonicalize to "free + enzyme → enzyme-met"
            # so two structurally-equivalent RE steps written in opposite
            # source directions dedup to the same Step. SS steps are NOT
            # canonicalized: their kf/kr labels are direction-sensitive
            # (analytical formulas reference :kNf as the source-forward
            # rate constant), so swapping would silently flip rate-equation
            # output. See CLAUDE.md "Canonical Step Form" for the
            # invariant.
            in_from = any(m -> m == bound_metabolite, bound(from_species))
            in_to   = any(m -> m == bound_metabolite, bound(to_species))
            if is_equilibrium && in_from && !in_to
                from_species, to_species = to_species, from_species
            end
        elseif is_equilibrium
            # RE iso steps: deterministic direction by lex on name(from_species).
            # SS iso steps are NOT canonicalized — their kf/kr labels are
            # direction-sensitive (analytical formulas reference :kNf as the
            # source-forward rate constant), so swapping would silently flip
            # rate-equation output. See CLAUDE.md "Canonical Step Form".
            if string(name(from_species)) > string(name(to_species))
                from_species, to_species = to_species, from_species
            end
        end
        new(from_species, to_species, bound_metabolite, is_equilibrium,
            source_idx)
    end
end

from_species(s::Step)     = s.from_species
to_species(s::Step)       = s.to_species
bound_metabolite(s::Step) = s.bound_metabolite
is_equilibrium(s::Step)   = s.is_equilibrium
is_binding(s::Step)       = s.bound_metabolite !== nothing
is_iso(s::Step)           = s.bound_metabolite === nothing
direction(s::Step)        = is_binding(s) ? :binding : :iso
# source_idx is presentation metadata for parameter naming, NOT part
# of Step's structural identity. Equality and hash IGNORE it so that
# two Steps with the same physics but different source positions
# compare equal — required by Mechanism dedup.
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
#
# The inner constructor assigns each Step's `source_idx` so the
# chokepoint can render today's source-order positional parameter
# names. If ALL incoming Steps have `source_idx == 0` (the default —
# the natural DSL / enumeration path), the constructor auto-assigns by
# flat position across groups. If ANY incoming Step has a non-zero
# `source_idx`, the constructor preserves it (the caller — e.g., the
# legacy-Sig lift — knows the true source position even after
# regrouping). Mixed mode is rejected to keep the convention
# unambiguous.
struct Mechanism
    reaction::EnzymeReaction
    steps::Vector{Vector{Step}}
    function Mechanism(reaction::EnzymeReaction,
                       steps::Vector{Vector{Step}})
        flat = Step[s for group in steps for s in group]
        any_set = any(s -> s.source_idx != 0, flat)
        all_set = all(s -> s.source_idx != 0, flat)
        any_set && !all_set &&
            error("Mechanism: mix of set and unset source_idx values " *
                  "in steps; pass either all zero (auto-assign) or all " *
                  "non-zero (preserve caller's source positions)")
        if any_set
            new(reaction, steps)
        else
            pos = 0
            renumbered = Vector{Vector{Step}}()
            for group in steps
                new_group = Step[]
                for s in group
                    pos += 1
                    push!(new_group,
                          Step(s.from_species, s.to_species,
                               s.bound_metabolite, s.is_equilibrium;
                               source_idx = pos))
                end
                push!(renumbered, new_group)
            end
            new(reaction, renumbered)
        end
    end
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
        # Mirror the non-allosteric `Mechanism` convention: if all
        # incoming Steps have `source_idx == 0` (the default), assign
        # by flat position across catalytic groups in source order;
        # otherwise preserve the caller's per-Step source positions
        # (reject mixed-mode).
        flat = Step[s for group in cat_steps for s in group]
        any_set = any(s -> s.source_idx != 0, flat)
        all_set = all(s -> s.source_idx != 0, flat)
        any_set && !all_set &&
            error("AllostericMechanism: mix of set and unset source_idx " *
                  "values in cat_steps; pass either all zero (auto-" *
                  "assign) or all non-zero (preserve caller's source " *
                  "positions)")
        if any_set
            new(reaction, cat_steps, cat_allo_states,
                catalytic_multiplicity, regulatory_sites)
        else
            pos = 0
            renumbered = Vector{Vector{Step}}()
            for group in cat_steps
                new_group = Step[]
                for s in group
                    pos += 1
                    push!(new_group,
                          Step(s.from_species, s.to_species,
                               s.bound_metabolite, s.is_equilibrium;
                               source_idx = pos))
                end
                push!(renumbered, new_group)
            end
            new(reaction, renumbered, cat_allo_states,
                catalytic_multiplicity, regulatory_sites)
        end
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
    s.source_idx,
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
    from_sig, to_sig, met_sig, is_eq, source_idx = sig
    met = met_sig === nothing ? nothing : _metabolite_from_sig(met_sig)
    Step(_species_from_sig(from_sig), _species_from_sig(to_sig), met, is_eq;
         source_idx = source_idx)
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
# `EnzymeMechanism(metabolites, reactions)`; `Mechanism(em)` lifts either
# shape so derivation consumers can walk a `Mechanism` uniformly.
# Lift a `Mechanism` to its singleton `EnzymeMechanism` type. The Sig
# stores each Step's data including `source_idx`, which the derivation
# pipeline keys parameter naming on (`_rep_idx_for_step`). To keep
# downstream code's "source_idx == position in `reactions(em)`"
# invariant intact across both Sig shapes — and avoid the index-by-
# `source_idx` mismatch when an enumeration move produces a Mechanism
# whose preserved `source_idx` no longer matches its current flat
# group-major position — renumber `source_idx` to flat-position before
# encoding. Mechanisms whose Steps already satisfy the invariant
# (fresh DSL builds, init_mechanisms output, the Mechanism constructor
# auto-assign path) are unaffected by the renumbering.
EnzymeMechanism(m::Mechanism) =
    EnzymeMechanism{_sig_of(
        _drop_unbound_regulators(_renumber_source_idx(m)))}()

# A regulator declared on the reaction that no step actually binds does
# not belong in the compiled catalytic mechanism's `regulators` list
# (e.g. a dead-end inhibitor before any expansion move binds it). Drop
# such regulators so they neither show up in `regulators(em)` nor get a
# parameter. Substrates/products are never dropped.
function _drop_unbound_regulators(m::Mechanism)
    bound_names = Set{Symbol}()
    for group in m.steps, s in group
        for sp in (from_species(s), to_species(s))
            for b in bound(sp)
                push!(bound_names, name(b))
            end
        end
        bm = bound_metabolite(s)
        bm === nothing || push!(bound_names, name(bm))
    end
    regs = regulators(m.reaction)
    kept = RegulatorMults[rm for rm in regs
                          if name(regulator(rm)) in bound_names]
    length(kept) == length(regs) && return m
    filtered_reaction = EnzymeReaction(
        reactants(m.reaction), kept,
        allowed_catalytic_multiplicities(m.reaction))
    Mechanism(filtered_reaction, m.steps)
end

function _renumber_source_idx(m::Mechanism)
    pos = 0
    new_steps = Vector{Vector{Step}}()
    for group in m.steps
        new_group = Step[]
        for s in group
            pos += 1
            push!(new_group,
                  Step(s.from_species, s.to_species,
                       s.bound_metabolite, s.is_equilibrium;
                       source_idx = pos))
        end
        push!(new_steps, new_group)
    end
    # Bypass the auto-assign branch by passing non-zero `source_idx` via
    # the constructor; the constructor preserves the explicit values.
    Mechanism(m.reaction, new_steps)
end

Mechanism(em::EnzymeMechanism{Sig}) where {Sig} = _mechanism_from_sig(Sig)

"""
Build a `Mechanism` from the legacy `(metabolites_3tuple, rxns_4tuple)`
Sig shape. Each enzyme-form Symbol becomes a conformation-only `Species`
with empty `bound` and empty `Residual`, so `name(species)` round-trips
to the original opaque Symbol exactly (`:E`, `:ES`, `:EAB`, ...). The
`bound_metabolite` field of `Step` carries the metabolite, while
`from_species`/`to_species` are taken straight from the source
`(lhs, rhs)` order. Because the metabolite is in neither bound list,
`Step`'s binding-direction canonicalization no-ops and source-order
direction is preserved verbatim.
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

    # Walk `rxns_sig` in source order and stamp each Step's
    # `source_idx` with its original 1-based source position before
    # grouping. The `Mechanism` constructor preserves these
    # explicit values so `_rep_idx_for_step` reproduces today's
    # source-order rep-idx naming (`:K1`, `:k10f`, ...).
    groups = Dict{Int, Vector{Step}}()
    group_order = Int[]
    for (src_pos, (lhs, rhs, is_eq, gnum)) in enumerate(rxns_sig)
        e_lhs = first(s for s in lhs if s ∉ met_set)
        e_rhs = first(s for s in rhs if s ∉ met_set)
        m_lhs = Symbol[s for s in lhs if s ∈ met_set]
        m_rhs = Symbol[s for s in rhs if s ∈ met_set]
        length(m_lhs) <= 1 ||
            error("_mechanism_from_legacy_sig: step $src_pos has more than " *
                  "one metabolite on LHS ($m_lhs); each Step must bind " *
                  "exactly one metabolite")
        length(m_rhs) <= 1 ||
            error("_mechanism_from_legacy_sig: step $src_pos has more than " *
                  "one metabolite on RHS ($m_rhs); each Step must bind " *
                  "exactly one metabolite")

        bound_met = !isempty(m_lhs) ? metabolite_of(m_lhs[1]) :
                    !isempty(m_rhs) ? metabolite_of(m_rhs[1]) :
                    nothing

        from = Species(Metabolite[], e_lhs, Residual())
        to   = Species(Metabolite[], e_rhs, Residual())
        step = Step(from, to, bound_met, is_eq; source_idx = src_pos)

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

"""
    AllostericMechanism(aem::AllostericEnzymeMechanism)

Lift an `AllostericEnzymeMechanism{CM, CS, RS}` to the non-parametric
`AllostericMechanism` struct. The catalytic mechanism is lifted via
`Mechanism(CM())`; `CS = (multiplicity, cat_allo_states)` provides the
catalytic-side data; each `RS` entry `(ligands, mult, reg_allo_states)`
becomes a `RegulatorySite` whose ligand `Symbol`s are wrapped as
`AllostericRegulator`. Mirrors `Mechanism(::EnzymeMechanism)` — bridges
the parametric ↔ non-parametric boundary so derivation code can walk
`AllostericMechanism` uniformly.
"""
function AllostericMechanism(
    ::AllostericEnzymeMechanism{CM, CS, RS},
) where {CM, CS, RS}
    cm_mech = Mechanism(CM())
    multiplicity, cat_allo_states = CS
    sites = RegulatorySite[]
    for entry in RS
        ligands_syms, mult, reg_allo_states = entry
        ligands_vec = AllostericRegulator[AllostericRegulator(l) for l in ligands_syms]
        push!(sites,
              RegulatorySite(ligands_vec, mult, collect(Symbol, reg_allo_states)))
    end
    AllostericMechanism(cm_mech.reaction, cm_mech.steps,
                        collect(Symbol, cat_allo_states),
                        multiplicity, sites)
end

"""
    AllostericEnzymeMechanism(am::AllostericMechanism)

Lift an `AllostericMechanism` to its singleton `AllostericEnzymeMechanism`
type. The catalytic side becomes an `EnzymeMechanism` lifting through
`Mechanism(am.reaction, am.cat_steps)` so the `source_idx`-renumbering
contract is preserved end-to-end. Catalytic and regulatory allosteric
data are encoded directly into the type parameters.
"""
function AllostericEnzymeMechanism(am::AllostericMechanism)
    cat_mech = Mechanism(am.reaction, am.cat_steps)
    cm = EnzymeMechanism(cat_mech)
    cat_sites = (am.catalytic_multiplicity,
                 Tuple(am.cat_allo_states))
    reg_sites = Tuple(
        (Tuple(Symbol[name(l) for l in site.ligands]),
         site.multiplicity,
         Tuple(site.allo_states))
        for site in am.regulatory_sites)
    AllostericEnzymeMechanism(cm, cat_sites, reg_sites)
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
#   `EnzymeMechanism(metabolites, reactions)`. The first slot
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


"""Return substrates as a tuple of `Symbol` names."""
@generated function substrates(::EnzymeMechanism{Sig}) where {Sig}
    names = Symbol[]
    for entry in Sig[1][1]
        kind, nm = entry[1]
        kind === :Substrate && push!(names, nm)
    end
    return Tuple(names)
end

"""Return products as a tuple of `Symbol` names."""
@generated function products(::EnzymeMechanism{Sig}) where {Sig}
    names = Symbol[]
    for entry in Sig[1][1]
        kind, nm = entry[1]
        kind === :Product && push!(names, nm)
    end
    return Tuple(names)
end

"""Return regulators as a tuple of `Symbol` names."""
@generated function regulators(::EnzymeMechanism{Sig}) where {Sig}
    names = Symbol[]
    for entry in Sig[1][2]
        _kind, nm = entry[1]
        push!(names, nm)
    end
    return Tuple(names)
end

"""
    metabolites(m::EnzymeMechanism) → Tuple{Symbol,...}

Return distinct metabolite names (substrates ∪ products ∪ regulators) as a tuple
of `Symbol`s in declaration order, deduplicated.
"""
@generated function metabolites(::EnzymeMechanism{Sig}) where {Sig}
    # Substrates first (from reactants Substrate entries),
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

"""Synthesize the legacy enzyme-form name from a Species sig
`(bound_sigs, conformation, residual_sig)` — mirrors `name(::Species)`.
Bound order matches `Species.bound`, which is sorted by name in the
struct constructor and preserved by `_to_sig`."""
function _species_name_from_sig(species_sig)
    bound_sigs, conformation, residual_sig = species_sig
    added_sig, sub_sig = residual_sig
    has_res = !(isempty(added_sig) && isempty(sub_sig))
    if isempty(bound_sigs) && !has_res
        return conformation
    end
    parts = String[String(conformation)]
    for b in bound_sigs
        push!(parts, b[1] === :CompetitiveInhibitor ?
                     String(b[2]) * "inh" : String(b[2]))
    end
    if has_res
        push!(parts, "res")
        for a in added_sig
            push!(parts, "+" * String(a[2]))
        end
        for r in sub_sig
            push!(parts, "-" * String(r[2]))
        end
    end
    Symbol(join(parts, "_"))
end

"""Build a (lhs, rhs, is_eq, g) tuple from a Sig step at @generated time.
The bound metabolite's side is set by its bound-list membership: present in
the to-species' bound list ⇒ it binds (lhs); present in the from-species'
list ⇒ it releases (rhs). When it is in neither list, an **SS** step on a
**decomposed** form (at least one side has a non-empty bound list) carrying a
product is a fused catalytic release — the product is produced by the step
(dissociation, rhs), since a retained metabolite would appear in a bound list.
Opaque forms (both bound lists empty — RE steps, or SS steps from RE→SS
enumeration moves) are canonicalized metabolite-on-binding-side, so they bind
via the bound-list-size heuristic; the non-empty-bounds guard keeps the
dissociation rule from mis-releasing those."""
function _step_tuple_from_sig(step_sig, g::Int)
    from_sig, to_sig, met_sig, is_eq, _src = step_sig
    e_from = _species_name_from_sig(from_sig)
    e_to   = _species_name_from_sig(to_sig)
    if met_sig === nothing
        return ((e_from,), (e_to,), is_eq, g)
    end
    met_name = met_sig[2]
    from_bound = from_sig[1]
    to_bound   = to_sig[1]
    bound_in_to   = any(b -> b == met_sig, to_bound)
    bound_in_from = any(b -> b == met_sig, from_bound)
    if bound_in_to
        return ((e_from, met_name), (e_to,), is_eq, g)
    elseif bound_in_from
        return ((e_from,), (e_to, met_name), is_eq, g)
    elseif !is_eq && met_sig[1] === :Product &&
           !(isempty(from_bound) && isempty(to_bound))
        return ((e_from,), (e_to, met_name), is_eq, g)
    elseif length(from_bound) > length(to_bound)
        return ((e_from,), (e_to, met_name), is_eq, g)
    end
    ((e_from, met_name), (e_to,), is_eq, g)
end

"""Return the reactions tuple `((lhs, rhs, is_eq, kinetic_group), ...)`."""
@generated function reactions(::EnzymeMechanism{Sig}) where {Sig}
    tuples = Any[]
    for (g, group) in enumerate(Sig[2])
        for step_sig in group
            push!(tuples, _step_tuple_from_sig(step_sig, g))
        end
    end
    return Tuple(tuples)
end

"""Return the equilibrium-step flags (`true` = rapid-equilibrium, `false` = steady-state)."""
@generated function equilibrium_steps(::EnzymeMechanism{Sig}) where {Sig}
    return Tuple(step[4] for group in Sig[2] for step in group)
end

"""Number of steps in the mechanism."""
function n_steps(::EnzymeMechanism{Sig}) where {Sig}
    sum(length(group) for group in Sig[2]; init=0)
end

"""Kinetic group of step `idx`."""
function kinetic_group(em::EnzymeMechanism{Sig}, idx::Int) where {Sig}
    flat = _flat_steps(Mechanism(em))
    1 ≤ idx ≤ length(flat) ||
        error("kinetic_group: step index $idx out of range 1:$(length(flat))")
    return flat[idx][2]
end

"""Sorted tuple of distinct kinetic group ids."""
@generated function kinetic_groups(::EnzymeMechanism{Sig}) where {Sig}
    return Tuple(1:length(Sig[2]))
end

"""Indices of steps belonging to kinetic group `G`."""
@generated function steps_in_group(
    ::EnzymeMechanism{Sig}, ::Val{G},
) where {Sig, G}
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
steps_in_group(m::EnzymeMechanism, g::Int) = steps_in_group(m, Val(g))

"""
    enzyme_forms(m::EnzymeMechanism) → Tuple{Symbol,...}

Return distinct enzyme-form names (any symbol appearing in a step that is not a
metabolite) as a tuple of `Symbol`s in step-order, deduplicated.
"""
@generated function enzyme_forms(::EnzymeMechanism{Sig}) where {Sig}
    # Collect metabolite names so a Species whose synthesized name
    # coincidentally matches a metabolite (e.g., bare `:S` conformation)
    # is excluded.
    met_names = Set{Symbol}()
    for entry in Sig[1][1]
        push!(met_names, entry[1][2])
    end
    for entry in Sig[1][2]
        push!(met_names, entry[1][2])
    end
    seen = Set{Symbol}()
    forms = Symbol[]
    for group in Sig[2]
        for step_sig in group
            from_sig, to_sig, _, _, _ = step_sig
            for sp_sig in (from_sig, to_sig)
                nm = _species_name_from_sig(sp_sig)
                nm ∉ met_names && nm ∉ seen &&
                    (push!(seen, nm); push!(forms, nm))
            end
        end
    end
    return Tuple(forms)
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
    # Walk via the accessors at body-build time and
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
# is the group's first step's `source_idx` — i.e., the lowest source
# index among the group's steps, matching the CLAUDE.md parameter
# naming convention ("rep step is 6, the lowest-indexed step in that
# group"). The `Mechanism` / `AllostericMechanism` constructors
# guarantee `source_idx` is populated by flat source-order position.
function _rep_idx_for_step(step::Step,
                           m::Union{Mechanism, AllostericMechanism})
    groups = m isa Mechanism ? m.steps : m.cat_steps
    for group in groups
        if step in group
            return first(group).source_idx
        end
    end
    error("Step not found in mechanism: $step")
end

_rep_idx_for_step(step::Step, m::EnzymeMechanism) =
    _rep_idx_for_step(step, Mechanism(m))
_rep_idx_for_step(step::Step, m::AllostericEnzymeMechanism) =
    _rep_idx_for_step(step, AllostericMechanism(m))

# Bridge: parametric AllostericEnzymeMechanism → non-parametric
# AllostericMechanism. Used by chokepoint dispatch (`name(p::Parameter,
# m::AllostericEnzymeMechanism)`) so a single Mechanism-family
# implementation handles both forms. Symmetric with
# `Mechanism(::EnzymeMechanism)`.
_to_mechanism(m::EnzymeMechanism)           = Mechanism(m)
_to_mechanism(m::AllostericEnzymeMechanism) = AllostericMechanism(m)

function _site_idx_of(site::RegulatorySite, m::AllostericMechanism)
    idx = findfirst(==(site), m.regulatory_sites)
    idx === nothing && error("RegulatorySite not found in mechanism")
    return idx
end

# Step-bound parameters map to three rendering rules keyed by type:
# Kd/Kiso → :K{rep},  Kon/Kfor → :k{rep}f,  Koff/Krev → :k{rep}r,
# with optional `_T` suffix for T-state. `_param_symbol` is the shared
# formatter — both the value-context `name(p::P, m)` and the
# type/index-context `name(::Type{P}, idx)` companion delegate here so a
# future parameter-naming refactor (semantic names like :K_ATP) changes
# one function body.
_param_symbol(::Type{<:Union{Kd, Kiso}},        idx::Int) = Symbol("K", idx)
_param_symbol(::Type{<:Union{Kon, Kfor}},       idx::Int) = Symbol("k", idx, "f")
_param_symbol(::Type{<:Union{Koff, Krev}},      idx::Int) = Symbol("k", idx, "r")

_param_symbol(::Type{P}, idx::Int, state::Symbol) where {P<:Parameter} =
    state === :T ? Symbol(_param_symbol(P, idx), "_T") :
                   _param_symbol(P, idx)

# Regulator-site Parameter — keyed by site index + ligand name (not the
# step-bound rep-idx), so it has its own formatter signature.
_param_symbol(::Type{Kreg}, site_idx::Int, lig_name::Symbol, state::Symbol) =
    state === :T ? Symbol("K_", lig_name, "_T_reg", site_idx) :
                   Symbol("K_", lig_name, "_reg",   site_idx)

# Type/index-context chokepoint companion. Used by @generated callers in
# `rate_eq_derivation.jl` where only an integer rep-idx is in scope (no
# Step value to construct a Parameter from). Supports only step-bound
# Parameter types (Kd/Kiso/Kon/Kfor/Koff/Krev) — Kreg names need a
# ligand and site, Keq/Etot/Lallo are stateless, so all four remain
# value-context only.
name(::Type{P}, idx::Int) where {P<:Parameter}                = _param_symbol(P, idx)
name(::Type{P}, idx::Int, state::Symbol) where {P<:Parameter} =
    _param_symbol(P, idx, state)

function name(p::StepBoundParameter,
              m::Union{Mechanism, EnzymeMechanism, AllostericMechanism})
    rep = _rep_idx_for_step(p.step, m)
    _param_symbol(typeof(p), rep, p.state)
end

# Regulator-site parameter — AllostericMechanism only.
function name(p::Kreg, m::AllostericMechanism)
    site_idx = _site_idx_of(p.site, m)
    _param_symbol(Kreg, site_idx, name(p.ligand), p.state)
end

# AllostericEnzymeMechanism dispatches via the AllostericMechanism lift
# so the rate-equation body and Parameter-name production share one
# implementation across parametric / non-parametric forms.
name(p::Kreg, m::AllostericEnzymeMechanism) = name(p, AllostericMechanism(m))
name(p::StepBoundParameter, m::AllostericEnzymeMechanism) =
    name(p, AllostericMechanism(m))

# Mechanism-level scalars
name(::Keq,   _) = :Keq
name(::Etot,  _) = :E_total
name(::Lallo, _) = :L

"""
Enumerate every raw rate-constant Parameter for a non-allosteric
mechanism, in kinetic-group order. Each kinetic group's representative
step (the first step in the group) drives the emit: RE binding → `Kd`,
RE iso → `Kiso`, SS binding → `Kon`+`Koff`, SS iso → `Kfor`+`Krev`. All
parameters carry `state === :None` because non-allosteric mechanisms
have no R/T branches.
"""
function _enumerate_parameters_full(m::Mechanism)
    out = Parameter[]
    for group in m.steps
        rep = first(group)
        if is_equilibrium(rep)
            push!(out, is_binding(rep) ? Kd(rep, :None) : Kiso(rep, :None))
        else
            if is_binding(rep)
                push!(out, Kon(rep, :None))
                push!(out, Koff(rep, :None))
            else
                push!(out, Kfor(rep, :None))
                push!(out, Krev(rep, :None))
            end
        end
    end
    out
end
