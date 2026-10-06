# ABOUTME: Core concrete types (EnzymeReaction, Mechanism, Step, Species,
# ABOUTME: Parameter family) plus the Mechanism↔Sig bridge and name(p,m) chokepoint.

# ─── Concrete type hierarchy ──────────────────────────────────────────

# Metabolite / Reactant / Regulator hierarchy.
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

name(m::Metabolite)::Symbol = m.name

"""
Residual: substrates added + products subtracted from the enzyme
(e.g., a covalent adduct after a ping-pong half-reaction). Empty
`Residual()` means no covalent residue.
"""
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

"""
Sort key for metabolite lists (`Species.bound`, `Step.consumed`,
`Step.released`): by name, a competitive-inhibitor copy after the
metabolite of the same name.
"""
_met_sort_key(m::Metabolite) = (name(m), m isa CompetitiveInhibitor)

"""
A metabolite as written in enzyme-form and parameter names: its name, with an
`inh` marker on a competitive-inhibitor copy, so the copy's forms and constants
stay distinct from those of the same metabolite in its reactant role.
"""
_met_label(m::Metabolite) =
    m isa CompetitiveInhibitor ? String(name(m)) * "inh" : String(name(m))

"""
Species: an enzyme form. `bound` is sorted by `_met_sort_key`; the
rendered Symbol name reads `:E` / `:EATP` / `:Estar...` / `:EATP_res_+P`.
The name is rendered once at construction and stored in `name`; it never
takes part in `==`, `hash` or display.
"""
struct Species
    bound::Vector{Metabolite}
    conformation::Symbol
    residual::Residual
    name::Symbol
    function Species(bound::Vector{<:Metabolite}, conformation::Symbol,
                     residual::Residual)
        sorted = sort(Vector{Metabolite}(bound); by = _met_sort_key)
        new(sorted, conformation, residual,
            _species_name(sorted, conformation, residual))
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
Base.show(io::IO, s::Species) = _show_fields(io, s, (:bound, :conformation, :residual))

# Julia's default display limited to the semantic fields: cached naming data never shows.
function _show_fields(io::IO, x, fields::Tuple)
    show(io, typeof(x))
    print(io, '(')
    recur_io = IOContext(io, :SHOWN_SET => x, :typeinfo => Any)
    for (i, f) in enumerate(fields)
        i > 1 && print(io, ", ")
        show(recur_io, getfield(x, f))
    end
    print(io, ')')
end

"""
Render species name deterministically from fields:
  :<conformation><bound1><bound2>...[_res[_+<added>...][_-<subtracted>...]]
Conformation and bound metabolites are concatenated without separator.
Metabolite Symbols must not contain `_` (domain convention): `:EATP`
unambiguously means E with ATP bound, not a conformation named "EATP".
Examples: `:E`, `:ES`, `:EATP`, `:EstarA_res_+P`.
"""
function _species_name(bound::Vector{Metabolite}, conformation::Symbol,
                       residual::Residual)
    head = String(conformation) * join(_met_label(m) for m in bound)
    isempty(residual) && return Symbol(head)
    Symbol(head, "_res", join("_+" * String(name(a)) for a in added(residual)),
           join("_-" * String(name(p)) for p in subtracted(residual)))
end
name(s::Species) = s.name

"""
RegulatorySite: a binding site (possibly multimeric) for one or
more allosteric ligands. `ligands[i]` and `allo_states[i]` are parallel;
ordering is meaningful (canonicalize at the call site if needed).
"""
struct RegulatorySite
    ligands::Vector{AllostericRegulator}
    multiplicity::Int
    allo_states::Vector{Symbol}
    function RegulatorySite(ligands::Vector{AllostericRegulator},
                            multiplicity::Int, allo_states::Vector{Symbol})
        isempty(ligands) && error("RegulatorySite: needs at least one ligand")
        length(ligands) == length(allo_states) ||
            error("RegulatorySite: length(ligands)=$(length(ligands)) " *
                  "must equal length(allo_states)=$(length(allo_states))")
        multiplicity ≥ 1 ||
            error("RegulatorySite: multiplicity must be ≥ 1, got $multiplicity")
        bad = findfirst(∉(_VALID_REG_ALLO_STATES), allo_states)
        bad === nothing ||
            error("RegulatorySite: allo state $(allo_states[bad]) must be one of " *
                  "$_VALID_REG_ALLO_STATES")
        new(ligands, multiplicity, allo_states)
    end
end

ligands(s::RegulatorySite)      = s.ligands
multiplicity(s::RegulatorySite) = s.multiplicity
allo_states(s::RegulatorySite)  = s.allo_states

"""
Whether `bound_form` is `free` with metabolite `m` added: same residual, and
`bound(bound_form)` equals `bound(free)` plus `m` as a multiset. The
conformation may differ (a binding may change the enzyme's conformation).
"""
_binds_ligand(free::Species, bound_form::Species, m::Metabolite) =
    residual(free) == residual(bound_form) &&
    bound(bound_form) == sort(Metabolite[bound(free)..., m]; by = _met_sort_key)

"""
    Step

One elementary transition between two enzyme `Species`, with its stoichiometry:
going `from_species → to_species` it takes up the metabolites in `consumed` from
solution and gives off those in `released`. Both lists may be empty (an
isomerization) or non-empty (for example a Theorell–Chance step EA + B → EQ + P).
`is_equilibrium` flags a rapid-equilibrium step (`true`) versus a steady-state
step (`false`). A binding (`bound_metabolite`) takes up exactly one metabolite and
gives off none. It is plain when `to_species` holds `from_species`'s metabolites
plus that one with the same residual, the conformation free to change
(`_binds_ligand`: E + A → E(A), E + A → E*(A)), and fused otherwise
(E(A) + B → E(P, Q)). The constructor stores a step that takes up nothing and gives
off exactly one metabolite as the binding it reverses, so every binding is stored
with its metabolite consumed; every other step is oriented by the `Mechanism` /
`AllostericMechanism` constructor. See CLAUDE.md "Canonical Step Form".
"""
struct Step
    from_species::Species
    to_species::Species
    consumed::Vector{Metabolite}
    released::Vector{Metabolite}
    is_equilibrium::Bool
    function Step(from_species::Species, to_species::Species,
                  consumed::Vector{<:Metabolite}, released::Vector{<:Metabolite},
                  is_equilibrium::Bool)
        from_species == to_species && error(
            "Step: both ends are the form $(name(from_species)); a step joins two forms")
        c = sort(Vector{Metabolite}(consumed); by = _met_sort_key)
        r = sort(Vector{Metabolite}(released); by = _met_sort_key)
        shared = intersect(c, r)
        isempty(shared) || error(
            "Step $(name(from_species)) → $(name(to_species)): " *
            "$(join(name.(shared), ", ")) both consumed and released")
        # A step that only gives off one metabolite is stored as the binding it
        # reverses (metabolite consumed), so a binding written in either direction
        # dedups to the same Step; binding constants are named after the step's
        # two sides, not its written direction. See CLAUDE.md "Canonical Step Form".
        if isempty(c) && length(r) == 1
            from_species, to_species, c, r = to_species, from_species, r, c
        end
        new(from_species, to_species, c, r, is_equilibrium)
    end
end

from_species(s::Step)   = s.from_species
to_species(s::Step)     = s.to_species
consumed(s::Step)       = s.consumed
released(s::Step)       = s.released
is_equilibrium(s::Step) = s.is_equilibrium

"""The metabolite a binding step binds: it takes up exactly that metabolite and
gives off none, as a plain or a fused binding (`Step`, `_is_chemistry`); `nothing`
for every other step."""
function bound_metabolite(s::Step)
    length(s.consumed) == 1 && isempty(s.released) ? only(s.consumed) : nothing
end
is_binding(s::Step) = bound_metabolite(s) !== nothing
is_iso(s::Step)     = isempty(s.consumed) && isempty(s.released)

"""Whether `s` is anything but a plain binding: one metabolite taken up, none given
off, and `to_species` holding `from_species`'s metabolites plus that one with the
same residual, the conformation free to change (`_binds_ligand`). Isomerizations,
fused bindings (`E(A) + B → E(P, Q)`), Theorell–Chance steps and steps with several
metabolites on a side are chemistry. The allosteric moves (`_expand_to_allosteric`,
`_partial_onlya_catalysis`, `_expand_change_allo_state`), `_onlya_haldane_violation`
and `show` read it to tell catalysis from binding."""
function _is_chemistry(s::Step)
    m = bound_metabolite(s)
    m === nothing || !_binds_ligand(from_species(s), to_species(s), m)
end

# Parameter family.
abstract type Parameter end

# Step-bound RE parameters
struct Kd   <: Parameter; step::Step; state::Symbol end
struct Kiso <: Parameter; step::Step; state::Symbol end

# Step-bound SS parameters
struct Kon  <: Parameter; step::Step; state::Symbol end
struct Koff <: Parameter; step::Step; state::Symbol end
struct Kfor <: Parameter; step::Step; state::Symbol end
struct Krev <: Parameter; step::Step; state::Symbol end

"""
Regulator-site parameter: a single ligand at a single site can appear
in either the A (active) or I (inactive) branch of the polynomial.
"""
struct Kreg <: Parameter
    site::RegulatorySite
    ligand::AllostericRegulator
    state::Symbol
end

# Mechanism-level scalar (singleton): the MWC coupling constant L
struct Lallo <: Parameter end

# Per-reactant and per-regulator bundling structs. Canonical
# ordering of atoms / multiplicities so two equivalent constructions
# compare equal under `==` / `hash`.
struct ReactantAtoms
    metabolite::Reactant
    atoms::Vector{Pair{Symbol,Int}}
    function ReactantAtoms(metabolite::Reactant,
                           atoms::Vector{<:Pair{Symbol,<:Integer}})
        isempty(atoms) && error(
            "ReactantAtoms: $(name(metabolite)) has no declared atoms; atoms " *
            "are mandatory (use `[C…]` bracket syntax in @enzyme_reaction).")
        bad = findfirst(a -> a.second isa Bool || a.second ≤ 0, atoms)
        bad === nothing ||
            error("ReactantAtoms: $(name(metabolite)) atom count for element " *
                  "$(atoms[bad].first) must be a positive integer, got " *
                  "$(atoms[bad].second).")
        new(metabolite, sort(Vector{Pair{Symbol,Int}}(atoms); by=first))
    end
end

metabolite(r::ReactantAtoms) = r.metabolite
atoms(r::ReactantAtoms)      = r.atoms

struct RegulatorMults
    regulator::Regulator
    allowed_multiplicities::Vector{Int}
    reg_type::Symbol
    function RegulatorMults(regulator::Regulator,
                            allowed_multiplicities::Vector{Int},
                            reg_type::Symbol = :unspecified)
        all(m -> m ≥ 1, allowed_multiplicities) ||
            error("RegulatorMults: allowed_multiplicities must all be ≥ 1, " *
                  "got $allowed_multiplicities")
        reg_type in (:activator, :inhibitor, :unspecified) ||
            error("RegulatorMults: reg_type must be :activator, :inhibitor, or " *
                  ":unspecified, got $reg_type")
        new(regulator, sort(allowed_multiplicities), reg_type)
    end
end

regulator(r::RegulatorMults)              = r.regulator
allowed_multiplicities(r::RegulatorMults) = r.allowed_multiplicities
reg_type(r::RegulatorMults)               = r.reg_type

"""
    EnzymeReaction

The public concrete reaction descriptor: substrates, products, optional
regulators, and the set of catalytic multiplicities the mechanism
enumerator is allowed to consider. It is the entry point to the package —
pair it with rate data in an [`IdentifyRateEquationProblem`](@ref), or pass
it to `EnzymeRates.init_mechanisms` to enumerate candidate mechanisms.

# Fields
- `reactants::Vector{ReactantAtoms}` — substrates and products, each
  carrying its per-atom inventory (used for ping-pong residual bookkeeping
  and atom conservation across steps).
- `regulators::Vector{RegulatorMults}` — competitive inhibitors and
  allosteric regulators, each with its allowed oligomeric multiplicities.
- `allowed_catalytic_multiplicities::Vector{Int}` — the oligomeric states
  the allosteric enumerator may assign to the catalytic core.
- `shared_catalytic_site::Vector{Tuple{Symbol,Symbol}}` — `(substrate,
  product)` pairs that bind the same catalytic site, and so are never both
  bound to it at once.

Construct one with the [`@enzyme_reaction`](@ref) DSL. Each reactant takes
an atom bracket (`S[C]`, `A[C1H1]`); the brackets are load-bearing for
ping-pong and multi-substrate reactions. Reactants and regulators are sorted
by name in the constructor, so two equivalent declarations compare equal
under `==`/`hash`.

```jldoctest
julia> using EnzymeRates

julia> rxn = @enzyme_reaction begin
           substrates: S[C]
           products:   P[C]
       end;

julia> rxn isa EnzymeReaction
true

julia> rxn
EnzymeReaction: S ⇌ P
```
"""
struct EnzymeReaction
    reactants::Vector{ReactantAtoms}
    regulators::Vector{RegulatorMults}
    allowed_catalytic_multiplicities::Vector{Int}
    shared_catalytic_site::Vector{Tuple{Symbol,Symbol}}

    function EnzymeReaction(reactants::Vector{ReactantAtoms},
                            regulators::Vector{RegulatorMults},
                            allowed_catalytic_multiplicities::Vector{Int};
                            shared_catalytic_site = Tuple{Symbol,Symbol}[])
        sorted_reactants = sort(reactants; by = ra -> name(metabolite(ra)))
        sorted_regulators = sort(regulators; by = rm -> name(regulator(rm)))
        sorted_mults = sort(unique(allowed_catalytic_multiplicities))
        all(m -> m ≥ 1, sorted_mults) ||
            error("EnzymeReaction: allowed_catalytic_multiplicities must " *
                  "all be ≥ 1, got $allowed_catalytic_multiplicities")

        sub_names  = Symbol[name(metabolite(ra)) for ra in sorted_reactants
                            if metabolite(ra) isa Substrate]
        prod_names = Symbol[name(metabolite(ra)) for ra in sorted_reactants
                            if metabolite(ra) isa Product]
        isempty(sub_names)  && error("EnzymeReaction: substrates must not be empty")
        isempty(prod_names) && error("EnzymeReaction: products must not be empty")

        reg_keys   = [(name(regulator(rm)), typeof(regulator(rm)))
                      for rm in sorted_regulators]
        allunique(sub_names)  || error("EnzymeReaction: duplicate substrate names")
        allunique(prod_names) || error("EnzymeReaction: duplicate product names")
        allunique(reg_keys)   ||
            error("EnzymeReaction: duplicate regulator of the same kind")

        sub_set  = Set(sub_names)
        prod_set = Set(prod_names)
        normalized_shared = Tuple{Symbol,Symbol}[]
        for pair in shared_catalytic_site
            length(pair) == 2 || error(
                "EnzymeReaction: each shared_catalytic_site entry must name " *
                "exactly two metabolites, got $pair")
            a, b = pair[1], pair[2]
            if a in sub_set && b in prod_set
                push!(normalized_shared, (a, b))
            elseif b in sub_set && a in prod_set
                push!(normalized_shared, (b, a))
            else
                error("EnzymeReaction: shared_catalytic_site pair $pair must " *
                      "name one declared substrate and one declared product")
            end
        end
        sorted_shared = sort(unique(normalized_shared))
        length(sorted_shared) == length(normalized_shared) || error(
            "EnzymeReaction: duplicate shared_catalytic_site pair")

        # A name may be declared once per regulator kind — the same metabolite
        # may bind BOTH as an AllostericRegulator and as a CompetitiveInhibitor
        # (the two roles render to distinct parameter names). A regulator may
        # also share a substrate/product name (so `concs.X` drives every role);
        # cross-category names need not be unique.

        # Atom mass-balance: per-element sum over substrate vs product
        # reactants. Dispatch by metabolite TYPE (not name) so a substrate and
        # product sharing a name route to the correct side.
        sub_atoms  = Dict{Symbol,Int}()
        prod_atoms = Dict{Symbol,Int}()
        for ra in sorted_reactants, (elem, c) in atoms(ra)
            tgt = metabolite(ra) isa Substrate ? sub_atoms : prod_atoms
            tgt[elem] = get(tgt, elem, 0) + c
        end
        for elem in union(keys(sub_atoms), keys(prod_atoms))
            s_c = get(sub_atoms, elem, 0)
            p_c = get(prod_atoms, elem, 0)
            s_c == p_c || error(
                "EnzymeReaction: atom imbalance — element $elem appears " *
                "$s_c time(s) on substrate side and $p_c on product side. " *
                "Declared atoms must balance.")
        end

        new(sorted_reactants, sorted_regulators, sorted_mults, sorted_shared)
    end
end

reactants(r::EnzymeReaction) = r.reactants
regulators(r::EnzymeReaction) = r.regulators
allowed_catalytic_multiplicities(r::EnzymeReaction) =
    r.allowed_catalytic_multiplicities
shared_catalytic_site(r::EnzymeReaction) = r.shared_catalytic_site
substrates(r::EnzymeReaction) =
    Substrate[metabolite(ra) for ra in reactants(r) if metabolite(ra) isa Substrate]
products(r::EnzymeReaction) =
    Product[metabolite(ra) for ra in reactants(r) if metabolite(ra) isa Product]

"""The distinct metabolite names of `rxn`: substrates, then products, then regulators,
each in stored order, keeping a name's first occurrence."""
_metabolite_names(rxn::EnzymeReaction) =
    unique!(Symbol[name.(substrates(rxn)); name.(products(rxn));
                   [name(regulator(rm)) for rm in regulators(rxn)]])

function Base.show(io::IO, r::EnzymeReaction)
    subs_str  = join(String.(name.(substrates(r))), " + ")
    prods_str = join(String.(name.(products(r))),   " + ")
    print(io, "EnzymeReaction: ", subs_str, " ⇌ ", prods_str)
    if !isempty(regulators(r))
        regs_str = join(
            (String(name(regulator(rm))) for rm in regulators(r)), ", ")
        print(io, " | regulators: ", regs_str)
    end
    mults = allowed_catalytic_multiplicities(r)
    if length(mults) == 1 && mults[1] > 1
        print(io, " | oligomeric_state: ", mults[1])
    elseif length(mults) > 1
        print(io, " | allowed_catalytic_multiplicities: (",
              join(mults, ", "), ")")
    end
    if !isempty(shared_catalytic_site(r))
        sites_str = join(
            ("($(p[1]), $(p[2]))" for p in shared_catalytic_site(r)), ", ")
        print(io, " | shared_catalytic_site: ", sites_str)
    end
end

"""
Canonicalize a non-binding step's storage direction to physical-forward, so
`from` is further from product-release / closer to substrate-binding.
Applies to RE AND SS steps — the direction question is identical;
only the parameter count differs. (All binding steps — RE and SS —
are stored with their metabolite consumed by the Step constructor; this
function orients every other step: isomerizations and transformations.)
Reversing a step swaps its forms and its lists, and every tier reads the
two sides symmetrically, so the result does not depend on how the step was
written.
"""
function _canonical_step_direction(s::Step, subs::Set{Symbol}, prods::Set{Symbol},
                                   kind::Dict{Species, Tuple{Bool, Bool}})
    is_binding(s) && return s
    f, t = from_species(s), to_species(s)
    flip() = Step(t, f, released(s), consumed(s), is_equilibrium(s))

    # Tier 1: atom-balance progression over bound plus free metabolites.
    score(sp, free) =
        (count(m -> name(m) in subs, bound(sp)) + count(m -> name(m) in subs, free),
         -count(m -> name(m) in prods, bound(sp)) - count(m -> name(m) in prods, free))
    sf, st = score(f, consumed(s)), score(t, released(s))
    sf > st && return s
    sf < st && return flip()

    # Tier 2: 1-hop (RE+SS) graph context — whether substrates and products
    # enter or leave solution at each form, as `(any substrate, any product)`.
    fk = get(kind, f, (false, false))
    tk = get(kind, t, (false, false))
    fk == (false, true) && tk == (true, false) && return s
    fk == (true, false) && tk == (false, true) && return flip()

    # Tier 3: lex fallback (source-independent).
    string(name(f)) ≤ string(name(t)) ? s : flip()
end

"""
Canonicalize the storage direction (RE + SS) of every non-binding step to
physical-forward for every group. Tier 2 reads each step's free metabolites
at both of its ends, so it sees the same context however the steps were
written. Shared by the `Mechanism` and `AllostericMechanism` constructors so
the Canonical Step Form invariant cannot drift between them.

The table `kind` built here classifies each species by the metabolites that
enter or leave solution at it, as `(any substrate, any product)`: the consumed
metabolites of every step leaving it (its `from_species`) and the released
metabolites of every step arriving at it (its `to_species`). A binding is stored
with its metabolite consumed, so it marks the form the metabolite binds to; a
fused release E(S) → F + P, stored as the binding F + P → E(S), marks F, the form
P leaves at. Reversing a step swaps its forms and its lists together, so the
classification does not depend on how any step was written. Isomerizations carry
no free metabolites and mark nothing. `_canonical_step_direction` Tier 2 reads it
to decide direction for non-binding steps where Tier 1 ties.

Why ALL steps (not just RE): the "substrate-entry / product-exit"
property is a chemistry fact about which forms metabolites enter and
leave at — it does NOT depend on whether the step is rapid-
equilibrium or steady-state. The DSL parses `<-->` as SS and `⇌` as RE;
fixtures like Segel Iso Uni Uni (`E + A <--> EA ⇌ EP <--> F + P, F <--> E`)
use `<-->` throughout, so an RE-only filter would mis-classify both `E`
and `F` as marking nothing and the F⇌E case would fall through to Tier 3 lex.
"""
function _canonicalize_step_directions(reaction::EnzymeReaction,
                                       groups::Vector{Vector{Step}})
    subs  = Set{Symbol}(name(s) for s in substrates(reaction))
    prods = Set{Symbol}(name(s) for s in products(reaction))
    kind = Dict{Species, Tuple{Bool, Bool}}()
    for group in groups, s in group, (form, free) in ((from_species(s), consumed(s)),
                                                      (to_species(s), released(s)))
        has_sub, has_prod = get(kind, form, (false, false))
        kind[form] = (has_sub  || any(m -> name(m) in subs, free),
                      has_prod || any(m -> name(m) in prods, free))
    end
    [[_canonical_step_direction(s, subs, prods, kind)
      for s in group] for group in groups]
end

# Canonical key for a `Step`. Gives the sorts a deterministic ordering so two
# physically-equivalent `Mechanism`s end up with identical step storage.
# Keyed on the species' rendered NAMES (not `hash`): `Substrate`/`Product`
# use Julia's default struct hash, which mixes in the type's `objectid` and so
# is NOT stable across precompile sessions — a hash key would make the canonical
# order (and thus `fitted_params` / the reduced rate equation / dedup keys)
# non-deterministic. Names are derived purely from metabolite names, so the
# ordering is fixed and reproducible.
_met_names(v) = join((String(name(m)) for m in v), "+")
_step_canonical_key(s::Step) =
    (String(name(from_species(s))), String(name(to_species(s))),
     _met_names(consumed(s)), _met_names(released(s)), is_equilibrium(s))

"""
Canonical key for a `RegulatorySite`. Orders the outer site vector so two
`AllostericMechanism`s differing only in site presentation order collapse.
Keyed on ligand names (stable) rather than `hash` for the same reason.
"""
_regulatory_site_canonical_key(site::RegulatorySite) =
    (Tuple(String(name(l)) for l in ligands(site)),
     multiplicity(site),
     Tuple(allo_states(site)))

"""
    _canonical_groups(reaction, groups) -> (groups, perm)

The kinetic groups in Canonical Step Form, as both mechanism constructors store
them. First reject an empty kinetic group. Then orient every step
(`_canonicalize_step_directions`), sort the steps within each group by
`_step_canonical_key`, then sort the groups by the canonical key of each group's
first step; `perm` is that group order (a permutation of 1:length),
which the caller applies to any data parallel to the groups. The inner sort must
run BEFORE the outer one so each group's "first step" key reflects the canonical
inner order. Each key is rendered once per sort. Then enforce the kinetic-group
rules (`_assert_uniform_groups`, `_assert_each_reaction_once`) and reject a
rapid-equilibrium segment with no bottom form (`_bottomless_re_segment`). The
returned groups are fresh vectors; the input is not mutated.
"""
function _canonical_groups(reaction::EnzymeReaction, groups::Vector{Vector{Step}})
    bad = findfirst(isempty, groups)
    bad === nothing ||
        error("Mechanism: kinetic group $bad is empty; a kinetic group holds at least " *
              "one step")
    gs = [g[sortperm(_step_canonical_key.(g))]
          for g in _canonicalize_step_directions(reaction, groups)]
    perm = sortperm([_step_canonical_key(first(g)) for g in gs])
    gs = gs[perm]
    _assert_uniform_groups(gs)
    _assert_each_reaction_once(gs)
    segment = _bottomless_re_segment(gs)
    segment === nothing ||
        error("Mechanism: rapid-equilibrium segment {" *
              join(name.(segment), ", ") * "} has no form free of one side's " *
              "metabolites, so its rate is undefined when that side is absent " *
              "(0/0). Make one of the segment's binding steps steady-state.")
    gs, perm
end

"""
A step side as written in parameter names: its enzyme form followed by its free
metabolites, joined by `_` (a competitive-inhibitor copy carries `inh`, as in
`name(::Species)`), e.g. `"EA_B"` for `E(A) + B`.
"""
_side_label(form::Species, mets::Vector{Metabolite}) =
    join([String(name(form)); _met_label.(mets)], "_")

"""
The two sides of `s` in its stored direction: `from_species` with the consumed
metabolites, then `to_species` with the released ones. Every step constant is
named after them.
"""
_forward_sides(s::Step) = (_side_label(from_species(s), consumed(s)),
                           _side_label(to_species(s), released(s)))

"""
The kind of a step for kinetic grouping: the metabolites it takes up and gives off,
`(consumed, released)`. Steps that share a kinetic group share their constants, so
they must take up and give off the same metabolites: two bindings of `m` share a
kind whether or not either runs chemistry, and every isomerization is one kind.
"""
_step_kind(s::Step) = (consumed(s), released(s))

"""
Error unless the steps of every kinetic group take up and give off the same
metabolites (`_step_kind`) and carry one RE/SS flag: a shared constant means the
same reaction type at the same speed.
"""
function _assert_uniform_groups(steps::Vector{Vector{Step}})
    label(s) = join(_forward_sides(s), " → ") * (is_equilibrium(s) ? " (RE)" : " (SS)")
    for group in steps
        s1 = first(group)
        for s in group
            _step_kind(s) == _step_kind(s1) && is_equilibrium(s) == is_equilibrium(s1) &&
                continue
            error("Mechanism: a kinetic group holds $(label(s1)) and $(label(s)), " *
                  "which differ in kind or in RE/SS; the steps of a kinetic group " *
                  "share their constants, so they must be the same kind of step " *
                  "with the same flag")
        end
    end
end

"""
Error when one reaction — the unordered pair of a step's two sides
(`_forward_sides`) — appears in more than one step. A reaction belongs to one
kinetic group. Two groups holding it with the same RE/SS flag would give one
reaction two sets of constants, which render the same names when the reaction
represents both groups; an RE group and an SS group holding it would make one
reaction both fast and slow. One group holding it twice would count its edge twice
in the derivation.
"""
function _assert_each_reaction_once(steps::Vector{Vector{Step}})
    seen = Dict{Tuple{String, String}, Tuple{Int, Step}}()
    flag(s) = is_equilibrium(s) ? "RE" : "SS"
    for (g, group) in enumerate(steps), s in group
        sides = minmax(_forward_sides(s)...)
        if haskey(seen, sides)
            g0, s0 = seen[sides]
            reaction = join(_forward_sides(s0), " ⇌ ")
            g0 == g && error("Mechanism: kinetic group $g holds the reaction " *
                             "$reaction twice; a mechanism holds each reaction once")
            flags = flag(s0) == flag(s) ? "" :
                " ($(flag(s0)) in group $g0, $(flag(s)) in group $g)"
            error("Mechanism: kinetic groups $g0 and $g both hold the reaction " *
                  "$reaction$flags; a reaction belongs to one kinetic group, whose " *
                  "constants are named after it")
        end
        seen[sides] = (g, s)
    end
end

"""
Rapid-equilibrium segments of `steps` with each form's metabolite exponents
cleared to the lowest in its segment: `(species, segments, extras)`, where
`segments[k]` indexes `species` and `extras[i][x]` is how many more `x` form `i`
carries than the lowest form of its segment. Within a segment, rapid equilibrium
fixes each form's weight relative to any other as a monomial in concentrations,
one factor per metabolite bound or released along the RE path between them;
chemistry steps pass the exponents through unchanged.
"""
function _re_segment_extras(steps::Vector{Vector{Step}})
    species = Species[]
    for group in steps, s in group
        from_species(s) in species || push!(species, from_species(s))
        to_species(s) in species   || push!(species, to_species(s))
    end
    idx = Dict(sp => i for (i, sp) in enumerate(species))
    # RE adjacency with the exponent change of each metabolite from → to.
    adj = [Tuple{Int, Dict{Symbol, Int}}[] for _ in species]
    for group in steps, s in group
        is_equilibrium(s) || continue
        d = Dict{Symbol, Int}()
        for m in consumed(s); d[name(m)] = get(d, name(m), 0) + 1; end
        for m in released(s); d[name(m)] = get(d, name(m), 0) - 1; end
        a, b = idx[from_species(s)], idx[to_species(s)]
        push!(adj[a], (b, d))
        push!(adj[b], (a, Dict(x => -e for (x, e) in d)))
    end
    expo = Dict{Int, Dict{Symbol, Int}}()
    segments = Vector{Int}[]
    for root in eachindex(species)
        haskey(expo, root) && continue
        expo[root] = Dict{Symbol, Int}()
        segment = [root]
        queue = [root]
        while !isempty(queue)
            u = popfirst!(queue)
            for (v, d) in adj[u]
                haskey(expo, v) && continue
                expo[v] = mergewith(+, expo[u], d)
                push!(segment, v); push!(queue, v)
            end
        end
        for x in union(Set{Symbol}(), (keys(expo[i]) for i in segment)...)
            lowest = minimum(get(expo[i], x, 0) for i in segment)
            for i in segment
                expo[i][x] = get(expo[i], x, 0) - lowest
            end
        end
        push!(segments, segment)
    end
    extras = [filter(p -> p.second > 0, expo[i]) for i in eachindex(species)]
    species, segments, extras
end

"""
    _bottomless_re_segment(steps) -> Union{Nothing, Vector{Species}}

The forms of a rapid-equilibrium segment whose weights all vanish at zero
products, or all vanish at zero substrates; `nothing` when no segment does.

Within a segment, rapid equilibrium fixes each form's weight relative to any
other as a monomial in concentrations (one factor per RE binding step on the
path between them). Once the weights are scaled to clear every negative
exponent, a segment normally has a bottom form whose weight is free of one
side's metabolites — the form the others drain to when that side is absent. A
segment with no such form (e.g. `{E(P), E(Q), E(P, Q)}` with `E + P` and
`E + Q` at steady state: weights `Kq·P : Kp·Q : P·Q`) is empty at zero
products, and the rate there depends on how `E(P, Q)` splits between `E(P)`
and `E(Q)` — a ratio of two fast release rates the RE approximation does not
carry. The rate equation is `0/0` exactly where initial-rate data sit. A
segment whose weights vanish only at a corner mixing a substrate and a product
(a mixed abortive complex, or the ping-pong `E` / `E(; residual)` pair) is
harmless: no turnover is possible at that corner.
"""
function _bottomless_re_segment(steps::Vector{Vector{Step}})
    side = Dict{Symbol, Type}()
    for group in steps, s in group
        is_equilibrium(s) || continue
        for m in Iterators.flatten((consumed(s), released(s)))
            side[name(m)] = typeof(m)
        end
    end
    species, segments, extras = _re_segment_extras(steps)
    for segment in segments
        length(segment) < 2 && continue
        for S in (Product, Substrate)
            all(i -> any(x -> get(side, x, Nothing) === S, keys(extras[i])),
                segment) && return species[segment]
        end
    end
    nothing
end

"""
Naming data derived from a mechanism's steps, filled on first use: the free-enzyme
form names, each kinetic group's naming representative and that representative's
rendered sides. It never takes part in `==`, `hash`, the compiled type or display.
"""
mutable struct _NamingCache
    free_enz::Union{Nothing, Set{Symbol}}
    reps::Union{Nothing, Vector{Step}}
    sides::Union{Nothing, Vector{Tuple{String, String}}}
end
_NamingCache() = _NamingCache(nothing, nothing, nothing)

"""
    Mechanism

A non-allosteric enzyme mechanism: a `reaction::EnzymeReaction` plus
`steps::Vector{Vector{Step}}`, where the outer vector is kinetic groups and
each inner vector holds the steps that share that group's kinetic parameters.
The constructor canonicalizes step direction and sorts steps and groups,
so two mechanisms that differ only in how their steps were written collapse to
the same struct. Parameter naming and step ordering derive purely from
structure and flat iteration order. Lift to the singleton derivation type with
`EnzymeRates.compile_mechanism(m)` or `EnzymeMechanism(m)`.
"""
struct Mechanism
    reaction::EnzymeReaction
    steps::Vector{Vector{Step}}
    naming::_NamingCache
    function Mechanism(reaction::EnzymeReaction,
                       steps::Vector{Vector{Step}})
        new(reaction, first(_canonical_groups(reaction, steps)), _NamingCache())
    end
end

steps(m::Mechanism) = m.steps
n_steps(m::Mechanism) = sum(length, m.steps; init = 0)
rep_step(m::Mechanism, g::Int) = first(m.steps[g])
Base.show(io::IO, m::Mechanism) = _show_fields(io, m, (:reaction, :steps))

# AllostericMechanism: a multi-subunit MWC enzyme. Each catalytic
# kinetic group carries an allosteric-state tag (`:OnlyA`, `:EqualAI`, or
# `:NonequalAI` — `:OnlyI` is rejected by the active-state convention).
const _VALID_CAT_ALLO_STATES = (:OnlyA, :EqualAI, :NonequalAI)
const _VALID_REG_ALLO_STATES = (:OnlyA, :OnlyI, :EqualAI, :NonequalAI)

struct AllostericMechanism
    reaction::EnzymeReaction
    cat_steps::Vector{Vector{Step}}
    cat_allo_states::Vector{Symbol}
    catalytic_multiplicity::Int
    regulatory_sites::Vector{RegulatorySite}
    naming::_NamingCache

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
        bad = findfirst(∉(_VALID_CAT_ALLO_STATES), cat_allo_states)
        bad === nothing ||
            error("AllostericMechanism: catalytic group $bad has invalid " *
                  "allo state $(cat_allo_states[bad]) (must be one of " *
                  "$_VALID_CAT_ALLO_STATES); :OnlyI is rejected for " *
                  "catalytic groups (active-state-active convention)")
        # cat_steps and cat_allo_states are parallel — permute both with the
        # same group order. regulatory_sites canonicalizes independently. All
        # come back as fresh vectors, so the caller's inputs are not mutated.
        cat_steps, perm = _canonical_groups(reaction, cat_steps)
        cat_allo_states = cat_allo_states[perm]
        regulatory_sites =
            sort(regulatory_sites; by = _regulatory_site_canonical_key)
        # Detect Kreg name collision: a ligand in two distinct regulatory
        # sites would produce identical rendered Kreg names (no site
        # discriminator). Not enumerated; constructor rejects it.
        ligs = [name(l) for site in regulatory_sites for l in ligands(site)]
        dup = findfirst(i -> ligs[i] in view(ligs, 1:i-1), eachindex(ligs))
        dup === nothing ||
            error("AllostericMechanism: ligand $(ligs[dup]) appears in two " *
                  "distinct regulatory sites; rendered Kreg names would " *
                  "collide. Same-ligand-two-sites is not enumerated.")
        violation = _onlya_haldane_violation(reaction, cat_steps, cat_allo_states)
        violation === nothing ||
            error("AllostericMechanism: $violation")
        new(reaction, cat_steps, cat_allo_states,
            catalytic_multiplicity, regulatory_sites, _NamingCache())
    end
end

# `==` and `hash` by content for the value types whose fields are Vectors or structs: the
# default `==` would compare identity. Every field takes part, in declaration order, except
# a mechanism's `naming` cache. `Species` keeps its own pair: its stored `name` is derived.
for T in (Residual, RegulatorySite, Step, Kd, Kiso, Kon, Koff, Kfor, Krev, Kreg,
          ReactantAtoms, RegulatorMults, EnzymeReaction, Mechanism, AllostericMechanism)
    fs = filter(!=(:naming), fieldnames(T))
    @eval Base.:(==)(a::$T, b::$T) =
        $(reduce((l, r) -> :($l && $r), [:(a.$f == b.$f) for f in fs]))
    @eval Base.hash(x::$T, h::UInt) =
        $(foldl((acc, f) -> :(hash(x.$f, $acc)), fs;
                init = :(hash($(QuoteNode(nameof(T))), h))))
end

reaction(m::Union{Mechanism, AllostericMechanism}) = m.reaction
kinetic_groups(m::Union{Mechanism, AllostericMechanism}) = 1:length(steps(m))
steps(m::AllostericMechanism) = m.cat_steps
cat_allo_state(m::AllostericMechanism, g::Int) = m.cat_allo_states[g]
cat_allo_states(m::AllostericMechanism) = m.cat_allo_states
catalytic_multiplicity(m::AllostericMechanism) = m.catalytic_multiplicity
regulatory_sites(m::AllostericMechanism) = m.regulatory_sites
n_steps(m::AllostericMechanism) = sum(length, m.cat_steps; init = 0)
rep_step(m::AllostericMechanism, g::Int) = first(m.cat_steps[g])

function allosteric_regulators(m::AllostericMechanism)
    seen = AllostericRegulator[]
    for site in regulatory_sites(m), lig in ligands(site)
        lig in seen || push!(seen, lig)
    end
    seen
end

Base.show(io::IO, m::AllostericMechanism) = _show_fields(io, m,
    (:reaction, :cat_steps, :cat_allo_states, :catalytic_multiplicity, :regulatory_sites))

"""
    _with(m; groups) / _with(am; groups, states, sites)

`m` rebuilt through its constructor with the given fields replaced: `groups`, the
kinetic groups; for an `AllostericMechanism` also `states`, the catalytic
allo-state tags parallel to `groups`, and `sites`, the regulatory sites. Fields not
given keep `m`'s. Expansion moves use it to change one part of a mechanism. The
constructor stores fresh group, state and site lists, so a child shares no list with
its parent, only the `Step` and `RegulatorySite` values in them.
"""
_with(m::Mechanism; groups = steps(m)) = Mechanism(reaction(m), groups)
_with(am::AllostericMechanism; groups = steps(am), states = cat_allo_states(am),
      sites = regulatory_sites(am)) =
    AllostericMechanism(reaction(am), groups, states, catalytic_multiplicity(am), sites)

# ─── Mechanism ↔ Sig (parametric ↔ non-parametric) conversion ──
#
# Every leaf in `sig` MUST be a valid Julia type-parameter value (isbits,
# Symbol, type, or Tuple of those). `Pair{Symbol,Int}` is NOT valid as a
# type parameter — encode pairs as `Tuple{Symbol,Int}`. Vectors are
# NEVER valid — always wrap in `Tuple(...)`.
#
# `_to_sig` encodes the leaves of a step, `_sig_of` assembles a mechanism's whole
# tuple from them, and `_mechanism_from_sig` decodes it.
#
# A mechanism's Sig is a tuple type of its own, so a decoder specialized on its
# argument compiles again for every mechanism. `_mechanism_from_sig` and its local
# decoders take their tuples `@nospecialize`, and the decoder indexes and loops over
# the whole mechanism's tuple, since destructuring or mapping a tuple also compiles
# once per tuple type.

"""
Encode a Sig leaf: a metabolite as `(TypeTag, name)` with `TypeTag =
nameof(typeof(m))`, a metabolite list as a tuple of those, a species as
`(bound, conformation, (added, subtracted))` and a step as `(from_species,
to_species, consumed, released, is_equilibrium)`.
"""
_to_sig(m::Metabolite) = (nameof(typeof(m)), name(m))
_to_sig(ms::Vector{<:Metabolite}) = Tuple(_to_sig(m) for m in ms)
_to_sig(s::Species) = (_to_sig(bound(s)), conformation(s),
                       (_to_sig(added(residual(s))), _to_sig(subtracted(residual(s)))))
_to_sig(s::Step) = (_to_sig(from_species(s)), _to_sig(to_species(s)),
                    _to_sig(consumed(s)), _to_sig(released(s)), is_equilibrium(s))

"""
The Sig of `m`: `(reaction_sig, steps_sig)`. `reaction_sig` holds the reactants as
`(leaf, atoms)`, the regulators as `(leaf, allowed_multiplicities)` and the allowed
catalytic multiplicities; `steps_sig` holds one tuple of step leaves per kinetic
group.

A regulator declared on the reaction that no step actually binds does not belong in
the compiled catalytic mechanism's `regulators` list (e.g. a dead-end inhibitor
before any expansion move binds it), so the Sig leaves it out: it neither shows up
in `regulators(em)` nor gets a parameter. Substrates and products are always
encoded. `Mechanism` (the working representation used during enumeration)
intentionally KEEPS unbound regulators — expansion moves bind them later.
"""
function _sig_of(m::Mechanism)
    rxn = reaction(m)
    bound_names = Set{Symbol}(
        name(x) for group in steps(m) for s in group
        for x in Iterators.flatten((bound(from_species(s)), bound(to_species(s)),
                                    consumed(s), released(s))))
    reaction_sig = (
        Tuple((_to_sig(metabolite(ra)), Tuple(Tuple.(atoms(ra)))) for ra in reactants(rxn)),
        Tuple((_to_sig(regulator(rm)), Tuple(allowed_multiplicities(rm)))
              for rm in regulators(rxn) if name(regulator(rm)) in bound_names),
        Tuple(allowed_catalytic_multiplicities(rxn)))
    (reaction_sig, Tuple(Tuple(_to_sig(s) for s in group) for group in steps(m)))
end

# The metabolite type each Sig leaf tag names.
const _SIG_METABOLITE_TYPES = (Substrate = Substrate, Product = Product,
                               AllostericRegulator = AllostericRegulator,
                               CompetitiveInhibitor = CompetitiveInhibitor)

"""The `Mechanism` whose Sig (`_sig_of`) is `sig`."""
function _mechanism_from_sig(@nospecialize(sig::Tuple))
    met(t::Tuple{Symbol, Symbol}) = _SIG_METABOLITE_TYPES[t[1]](t[2])
    mets(@nospecialize(ts::Tuple)) = Metabolite[met(t) for t in ts]
    species(@nospecialize(t::Tuple)) = Species(mets(t[1]), t[2],
        Residual(Substrate[met(x) for x in t[3][1]], Product[met(x) for x in t[3][2]]))
    reactants_sig, regulators_sig, mults_sig = sig[1]
    rxn = EnzymeReaction(
        ReactantAtoms[
            ReactantAtoms(met(r[1]), Pair{Symbol, Int}[a => n for (a, n) in r[2]])
            for r in reactants_sig],
        RegulatorMults[RegulatorMults(met(r[1]), collect(Int, r[2])) for r in regulators_sig],
        collect(Int, mults_sig))
    groups = Vector{Step}[]
    for group_sig in sig[2]
        group = Step[]
        for s in group_sig
            push!(group, Step(species(s[1]), species(s[2]), mets(s[3]), mets(s[4]), s[5]))
        end
        push!(groups, group)
    end
    Mechanism(rxn, groups)
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
  (from-species, to-species, consumed metabolites, released metabolites,
  is-equilibrium).

The conversion functions are `_sig_of` and `_mechanism_from_sig`; the
exact layout is internal — users construct via `EnzymeMechanism(::Mechanism)`.
"""
struct EnzymeMechanism{Sig} <: AbstractEnzymeMechanism end

# Boundary converters: non-parametric Mechanism ↔ parametric
# EnzymeMechanism{Sig}. The Sig encodes the Mechanism's data as
# `(reaction_sig, steps_sig)` produced by `_sig_of`; `Mechanism(em)`
# lifts it back so derivation consumers can walk a `Mechanism`
# uniformly.
"""
Lift a `Mechanism` to its singleton `EnzymeMechanism` type. The Sig
is purely structural — two mechanisms differing only in source order
collapse to the same `EnzymeMechanism` type.
"""
EnzymeMechanism(m::Mechanism) = EnzymeMechanism{_sig_of(m)}()

# Reads `Sig` from the type at run time, so the lift compiles once rather than once
# per mechanism type.
"""
Lift an `EnzymeMechanism` to the `Mechanism` its Sig encodes. The lift is not the inverse
of `EnzymeMechanism(m)`: `shared_catalytic_site`, `RegulatorMults.reg_type` and unbound
regulators are not encoded, on purpose, so that equivalent mechanisms share one compile.
"""
Mechanism(@nospecialize(em::EnzymeMechanism)) =
    _mechanism_from_sig(typeof(em).parameters[1])

"""
    AllostericEnzymeMechanism{CatalyticMech, CatSites, RegSites}

Internal singleton type used by `@generated rate_equation` for
allosteric MWC enzymes. User code constructs and inspects via
`AllostericMechanism`; `compile_mechanism(am)` produces this opaque
fast-path handle. The three type parameters encode the catalytic
mechanism plus per-site allosteric data.

The three type parameters cannot be folded into a single value-tuple `Sig`
(as EnzymeMechanism{Sig} does): the first slot is a DataType (an
EnzymeMechanism subtype), and Julia rejects DataTypes inside the
value-tuple position of a type parameter.
"""
struct AllostericEnzymeMechanism{
    CatalyticMech, CatSites, RegSites,
} <: AbstractEnzymeMechanism end

"""
    AllostericMechanism(aem::AllostericEnzymeMechanism)

Lift an `AllostericEnzymeMechanism{CM, CS, RS}` to the non-parametric
`AllostericMechanism` struct. The catalytic mechanism is lifted via
`Mechanism(CM())`; `CS = (multiplicity, cat_allo_states)` provides the
catalytic-side data; each `RS` entry `(ligands, mult, reg_allo_states)`
becomes a `RegulatorySite` whose ligand `Symbol`s are wrapped as
`AllostericRegulator`. Mirrors `Mechanism(::EnzymeMechanism)` — bridges
the parametric ↔ non-parametric boundary so derivation code can walk
`AllostericMechanism` uniformly. It reads the type parameters at run time,
so it compiles once rather than once per mechanism type.
"""
function AllostericMechanism(@nospecialize(aem::AllostericEnzymeMechanism))
    CM, CS, RS = typeof(aem).parameters
    cm_mech = Mechanism(CM())
    multiplicity, cat_allo_states = CS
    sites = RegulatorySite[
        RegulatorySite(AllostericRegulator[AllostericRegulator(l) for l in ligands_syms],
                       mult, collect(Symbol, reg_allo_states))
        for (ligands_syms, mult, reg_allo_states) in RS]
    AllostericMechanism(reaction(cm_mech), steps(cm_mech),
                        collect(Symbol, cat_allo_states),
                        multiplicity, sites)
end

"""
    AllostericEnzymeMechanism(am::AllostericMechanism)

Lift an `AllostericMechanism` to its singleton `AllostericEnzymeMechanism`
type. The catalytic side becomes an `EnzymeMechanism` lifting through
`Mechanism(am.reaction, am.cat_steps)`. Catalytic and regulatory allosteric
data are encoded directly into the type parameters: `CatSites = (multiplicity,
cat_allo_states)` and one `(ligands, multiplicity, allo_states)` entry per
regulatory site. The `AllostericMechanism` and `RegulatorySite` constructors have
validated every field, so the lift checks nothing further.
"""
function AllostericEnzymeMechanism(am::AllostericMechanism)
    cm = EnzymeMechanism(Mechanism(reaction(am), steps(am)))
    AllostericEnzymeMechanism{
        typeof(cm),
        (catalytic_multiplicity(am), Tuple(cat_allo_states(am))),
        Tuple((Tuple(name.(ligands(s))), multiplicity(s), Tuple(allo_states(s)))
              for s in regulatory_sites(am)),
    }()
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

Full rate equation mode using all 2N microscopic rate constants (k_ES_to_EP, k_EP_to_ES, ...).
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

    # Render the steps as a single chain by walking the enzyme-form graph
    # (stored order is canonical, not chain order). Each step is an edge
    # between its two enzyme forms; start at a path endpoint (a degree-1
    # form) if any, else the free enzyme `:E`, else any form, then follow
    # edges and emit each step's far side. If the walk can't consume every
    # step (the mechanism is branched, or the chain would hide a step's
    # metabolite) → multi-line rendering below.
    _enz_forms(lhs, rhs) = (first(s for s in lhs if s in enz_set),
                            first(s for s in rhs if s in enz_set))
    degree = Dict{Symbol,Int}()
    for (lhs, rhs, _, _) in Rxns
        a, b = _enz_forms(lhs, rhs)
        degree[a] = get(degree, a, 0) + 1
        degree[b] = get(degree, b, 0) + 1
    end
    start = nothing
    for (lhs, rhs, _, _) in Rxns
        a, b = _enz_forms(lhs, rhs)
        degree[a] == 1 && (start = a; break)
        degree[b] == 1 && (start = b; break)
    end
    start === nothing && :E in enz_set && (start = :E)
    start === nothing && !isempty(Rxns) &&
        (start = _enz_forms(Rxns[1][1], Rxns[1][2])[1])

    # A single chain exists iff every enzyme form has degree ≤ 2 (a simple
    # path or cycle); a higher-degree form is a branch point → multi-line.
    chain_segments = String[]
    chain_arrows = String[]
    is_linear = !isempty(Rxns) && all(<=(2), values(degree))
    if is_linear
        subs = Set{Symbol}(substrates(m))
        remaining = collect(Rxns)
        remaining_binds =
            [!_is_chemistry(s) for group in steps(Mechanism(m)) for s in group]
        current = start
        while !isempty(remaining)
            idx = nothing
            if isempty(chain_segments)
                # First step: prefer to leave `current` by binding a substrate,
                # so a reversible cycle renders substrate→product.
                idx = findfirst(remaining) do rxn
                    fs = _enz_forms(rxn[1], rxn[2])
                    current in fs &&
                        any(x -> x in subs, current == fs[1] ? rxn[1] : rxn[2])
                end
            end
            if idx === nothing
                idx = findfirst(
                    rxn -> current in _enz_forms(rxn[1], rxn[2]), remaining)
            end
            idx === nothing && (is_linear = false; break)
            lhs, rhs, is_eq, _ = remaining[idx]
            binds = remaining_binds[idx]
            deleteat!(remaining, idx)
            deleteat!(remaining_binds, idx)
            a, b = _enz_forms(lhs, rhs)
            # Only the first step prints its entry side. A later step may leave
            # its entry-side metabolites unprinted only if it is a plain binding,
            # whose bound form names the metabolite; any other step (e.g. a fused
            # binding or a Theorell–Chance step) needs the multi-line rendering.
            !isempty(chain_segments) && length(current == a ? lhs : rhs) > 1 &&
                !binds && (is_linear = false; break)
            in_side  = current == a ? join(lhs, " + ") : join(rhs, " + ")
            out_side = current == a ? join(rhs, " + ") : join(lhs, " + ")
            isempty(chain_segments) && push!(chain_segments, in_side)
            push!(chain_arrows, _arrow(is_eq))
            push!(chain_segments, out_side)
            current = current == a ? b : a
        end
        !isempty(remaining) && (is_linear = false)
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
# Accessors on `EnzymeMechanism` lift to `Mechanism(em)` and walk the
# concrete `reaction` / `steps` fields. Return shapes (tuples of
# `Symbol`s, the `Matrix{Int}` stoichiometry, ranges) are the contract
# consumed by the @generated rate-equation body builders.

"""Walk the steps of `m` in flat order, yielding
`(step::Step, kinetic_group::Int)` pairs."""
_flat_steps(m::Mechanism) = [(s, g) for (g, group) in enumerate(steps(m)) for s in group]

"""Return substrates as a tuple of `Symbol` names."""
function substrates(em::EnzymeMechanism)
    m = Mechanism(em)
    Tuple(name(s) for s in substrates(reaction(m)))
end

"""Return products as a tuple of `Symbol` names."""
function products(em::EnzymeMechanism)
    m = Mechanism(em)
    Tuple(name(p) for p in products(reaction(m)))
end

"""Return regulators as a tuple of `Symbol` names."""
function regulators(em::EnzymeMechanism)
    m = Mechanism(em)
    Tuple(name(regulator(rm)) for rm in regulators(reaction(m)))
end

"""
    metabolites(m::EnzymeMechanism) → Tuple{Symbol,...}

Return distinct metabolite names (substrates ∪ products ∪ regulators) as a tuple
of `Symbol`s in declaration order, deduplicated. The fitter uses this as the
key set for the per-datapoint concentration `NamedTuple` passed to
[`rate_equation`](@ref).

```jldoctest
julia> using EnzymeRates

julia> m = @enzyme_mechanism begin
           substrates: S
           products: P
           steps: begin
               E + S <--> E(S)
               E(S) <--> E + P
           end
       end;

julia> metabolites(m)
(:S, :P)
```
"""
@generated function metabolites(::EnzymeMechanism{Sig}) where {Sig}
    # Kept `@generated`: `loss!` uses `metabolites(m)` as a compile-time-constant
    # tuple to build the per-datapoint `NamedTuple{MetNames}` concs on the fitting
    # hot path. A runtime body would make that NamedTuple type-unstable and allocate.
    Tuple(_metabolite_names(reaction(Mechanism(EnzymeMechanism{Sig}()))))
end

"""Return the reactions tuple `((lhs, rhs, is_eq, kinetic_group), ...)`.
Each step's `lhs` is its `from_species` name followed by the names of the
metabolites it consumes; `rhs` is its `to_species` name followed by the names
of the metabolites it releases."""
function reactions(em::EnzymeMechanism)
    m = Mechanism(em)
    tuples = Any[]
    for (g, group) in enumerate(steps(m))
        for s in group
            push!(tuples,
                  ((name(from_species(s)), name.(consumed(s))...),
                   (name(to_species(s)), name.(released(s))...),
                   is_equilibrium(s), g))
        end
    end
    return Tuple(tuples)
end

"""Return the equilibrium-step flags (`true` = rapid-equilibrium, `false` = steady-state)."""
function equilibrium_steps(em::EnzymeMechanism)
    m = Mechanism(em)
    Tuple(is_equilibrium(s) for group in steps(m) for s in group)
end

"""Number of steps in the mechanism."""
function n_steps(em::EnzymeMechanism)
    m = Mechanism(em)
    sum(length, steps(m); init=0)
end

"""Kinetic group of step `idx`."""
function kinetic_group(em::EnzymeMechanism, idx::Int)
    flat = _flat_steps(Mechanism(em))
    1 ≤ idx ≤ length(flat) ||
        error("kinetic_group: step index $idx out of range 1:$(length(flat))")
    return flat[idx][2]
end

"""Sorted tuple of distinct kinetic group ids."""
function kinetic_groups(em::EnzymeMechanism)
    m = Mechanism(em)
    Tuple(1:length(steps(m)))
end

"""Indices of steps belonging to kinetic group `g`."""
function steps_in_group(em::EnzymeMechanism, g::Int)
    flat = _flat_steps(Mechanism(em))
    Tuple(i for (i, (_, gid)) in enumerate(flat) if gid == g)
end
steps_in_group(em::EnzymeMechanism, ::Val{G}) where {G} = steps_in_group(em, G)

"""
    enzyme_forms(m::EnzymeMechanism) → Tuple{Symbol,...}

Return distinct enzyme-form names (any symbol appearing in a step that is not a
metabolite) as a tuple of `Symbol`s in step-order, deduplicated.
"""
function enzyme_forms(em::EnzymeMechanism)
    # Collect metabolite names so a Species whose synthesized name
    # coincidentally matches a metabolite (e.g., bare `:S` conformation)
    # is excluded.
    m = Mechanism(em)
    met_names = Set(metabolites(em))
    seen = Set{Symbol}()
    forms = Symbol[]
    for group in steps(m), s in group
        for sp in (from_species(s), to_species(s))
            nm = name(sp)
            nm ∉ met_names && nm ∉ seen &&
                (push!(seen, nm); push!(forms, nm))
        end
    end
    return Tuple(forms)
end

"""Number of distinct enzyme states."""
n_states(m::EnzymeMechanism) = length(enzyme_forms(m))

# ─── AllostericEnzymeMechanism Accessors ────────────────────────

catalytic_mechanism(::AllostericEnzymeMechanism{CM}) where {CM} = CM()
catalytic_multiplicity(::AllostericEnzymeMechanism{CM, CS}) where {CM, CS} = CS[1]

"""Return the allosteric state of catalytic kinetic group `g`."""
function cat_allo_state(::AllostericEnzymeMechanism{CM, CS, RS}, g::Int) where {CM, CS, RS}
    _, states = CS
    return states[g]
end

# These accessors forward to the catalytic_mechanism. Generated en masse
# to avoid 13 lines of boilerplate.
for fn in (:substrates, :products, :reactions, :equilibrium_steps,
           :n_steps, :enzyme_forms, :n_states, :kinetic_groups)
    @eval $fn(m::AllostericEnzymeMechanism) = $fn(catalytic_mechanism(m))
end
kinetic_group(m::AllostericEnzymeMechanism, i::Int) =
    kinetic_group(catalytic_mechanism(m), i)
steps_in_group(m::AllostericEnzymeMechanism, g) =
    steps_in_group(catalytic_mechanism(m), g)

"""
Returns ONLY reg-site ligands, NOT a union with catalytic_mechanism's
regulators. Downstream rate-equation code reads `regulators(m)` to find
dead-end binding K's; including allosteric-only ligands would cause it
to look up nonexistent K names.
"""
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

# `@generated` for the same reason as `metabolites(::EnzymeMechanism)`: `loss!` needs
# the tuple as a compile-time constant.
@generated function metabolites(::AllostericEnzymeMechanism{CM, CS, RS}) where {CM, CS, RS}
    ligs = (l for (site_ligands, _, _) in RS for l in site_ligands)
    Tuple(unique!(Symbol[metabolites(CM())..., ligs...]))
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

regulatory_sites(::AllostericEnzymeMechanism{CM, CS, RS}) where {CM, CS, RS} = RS

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
# Single chokepoint for parameter Symbol production. Renders structural
# names (:K_ES_to_E_S, :k_ES_to_EP, :K_I_ES_to_E_S, :K_Ireg) via Parameter
# subtype dispatch. Routing all parameter-name production through one function
# keeps any name-scheme change a single-function edit.

"""
State token — placed right after the type prefix:
  :A       → "A_"  (allosteric active branch; distinguishes allosteric from non-)
  :I       → "I_"  (allosteric inactive branch: OnlyI or NonequalAI-inactive)
  :EqualAI → ""    (allosteric shared symbol; shared between A and I state)
  :None    → ""    (non-allosteric mechanism)
"""
_state_tag(state::Symbol) =
    state === :A ? "A_" : state === :I ? "I_" : state in (:EqualAI, :None) ? "" :
    error("_state_tag: unexpected Parameter.state $state " *
          "(must be one of :None, :EqualAI, :A, :I)")

"""
The constant of the reaction `a → b`: `prefix`, the state tag, then
`"<a>_to_<b>"` (e.g. `:k_ES_to_EP`, `:K_A_ES_to_E_S`).
"""
_render_reaction(prefix::String, (a, b)::Tuple{String, String}, state::Symbol) =
    Symbol(prefix, _state_tag(state), a, "_to_", b)

"""Each kinetic group's naming representative, in group order; cached per mechanism."""
function _group_reps(m::Union{Mechanism, AllostericMechanism})
    c = m.naming
    reps = c.reps
    reps === nothing || return reps
    fes = _free_enz_set(m)
    c.reps = Step[_group_rep(group, fes) for group in steps(m)]
end

"""Each kinetic group's naming representative's `_forward_sides`, in group order; cached
per mechanism."""
function _group_sides(m::Union{Mechanism, AllostericMechanism})
    c = m.naming
    sides = c.sides
    sides === nothing || return sides
    c.sides = Tuple{String, String}[_forward_sides(rep) for rep in _group_reps(m)]
end

"""Find the kinetic group containing `step`; return its naming rep's `_forward_sides`."""
function _rep_sides(step::Step, m::Union{Mechanism, AllostericMechanism})
    for (g, group) in enumerate(steps(m))
        step in group && return _group_sides(m)[g]
    end
    error("Step not found in mechanism: $step")
end

const _AnyMech = Union{Mechanism, AllostericMechanism}

# Every step constant is named after the reaction of its group's representative:
# a rate constant (`k_`) after the direction it drives, an equilibrium constant
# (`K_`) after the direction whose products-over-reactants it is. A binding's K is
# named in the release direction, so it is a dissociation constant (`K_ES_to_E_S`);
# an isomerization's or other RE step's K in the stored direction (`K_ES_to_EP`).
name(p::Union{Kon, Kfor}, m::_AnyMech) =
    _render_reaction("k_", _rep_sides(p.step, m), p.state)
name(p::Union{Koff, Krev}, m::_AnyMech) =
    _render_reaction("k_", reverse(_rep_sides(p.step, m)), p.state)
name(p::Kiso, m::_AnyMech) =
    _render_reaction("K_", _rep_sides(p.step, m), p.state)
name(p::Kd, m::_AnyMech) =
    _render_reaction("K_", reverse(_rep_sides(p.step, m)), p.state)

"""
Regulator-site parameter: state tag + ligand name + "reg". No site index —
the AllostericMechanism constructor enforces that each ligand appears at most
once across all sites.
"""
name(p::Kreg, ::AllostericMechanism) =
    Symbol("K_", _state_tag(p.state), String(name(p.ligand)), "reg")

# Mechanism-level scalar
name(::Lallo, _) = :L

"""
Enumerate every raw rate-constant Parameter for a non-allosteric
mechanism, in kinetic-group order. Each kinetic group's representative
step (the structurally-primary step, `_group_rep`) drives the emit: RE
binding → `Kd`, RE non-binding step → `Kiso`, SS binding → `Kon`+`Koff`, SS
non-binding step → `Kfor`+`Krev`. All parameters carry `state === :None`
because non-allosteric mechanisms have no A/I branches.
"""
_enumerate_parameters_full(m::Mechanism) =
    Parameter[p for rep in _group_reps(m) for p in _step_params(rep, :None)]

"""
The Parameter(s) governing step `s` with the given allosteric state: the 4-way
switch on `is_equilibrium(s)` × `is_binding(s)`. The walkers over group
representatives (`_enumerate_parameters_full`, `_ss_rate_constant_names`) apply it
to each group's representative step.

Returns 1 element for RE steps (`Kd` or `Kiso`) and 2 elements for SS
steps (`Kon`+`Koff` or `Kfor`+`Krev`). A binding, plain or fused, takes `Kd` /
`Kon`+`Koff`; every other step — an isomerization, a Theorell–Chance step, a step
with two metabolites on one side — takes `Kiso` / `Kfor`+`Krev`.
"""
_step_params(s::Step, state::Symbol) =
    is_equilibrium(s) ? Parameter[(is_binding(s) ? Kd : Kiso)(s, state)] :
    is_binding(s) ? Parameter[Kon(s, state), Koff(s, state)] :
    Parameter[Kfor(s, state), Krev(s, state)]
