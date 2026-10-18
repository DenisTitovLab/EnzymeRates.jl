# ABOUTME: Value types: EnzymeReaction, Step and Species, and the mechanism as a reaction
# ABOUTME: graph in canonical form, with its raw input, checks, identity, Sig and views.

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
Species: an enzyme form. `bound` is sorted by `_met_sort_key`. Its internal name
(`_species_name`: `:E` / `:EATP` / `:Estar...` / `:EB_res_+A_-P`) is rendered once at
construction from the form alone and stored in `name`; the `Step` constructor's errors
name a form by it, and it never takes part in `==`, `hash` or display. Parameter
names write a form with `_form_label`, which also reads the mechanism's reactants.
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
The internal name of a form (`name(::Species)`), rendered from its fields alone:
  :<conformation><bound1><bound2>...[_res[_+<added>...][_-<subtracted>...]]
Conformation and bound metabolites are concatenated without separator, a
competitive-inhibitor copy with an `inh` marker. Metabolite Symbols must not contain
`_` (domain convention): `:EATP` unambiguously means E with ATP bound, not a
conformation named "EATP". Examples: `:E`, `:ES`, `:EATP`, `:EIinh`, `:EstarB_res_+A_-P`.
"""
function _species_name(bound::Vector{Metabolite}, conformation::Symbol,
                       residual::Residual)
    head = String(conformation) *
           join(m isa CompetitiveInhibitor ? String(name(m)) * "inh" : String(name(m))
                for m in bound)
    isempty(residual) && return Symbol(head)
    Symbol(head, "_res", join("_+" * String(name(a)) for a in added(residual)),
           join("_-" * String(name(p)) for p in subtracted(residual)))
end
name(s::Species) = s.name

"""
RegulatorySite: a binding site (possibly multimeric) for one or
more allosteric ligands, each given by its name, once. `ligands[i]` and `allo_states[i]`
are parallel. The constructor sorts the ligands by name and carries each ligand's state
along, so the order the ligands are given in never changes `==`, `hash` or
the compiled type.
"""
struct RegulatorySite
    ligands::Vector{Symbol}
    multiplicity::Int
    allo_states::Vector{Symbol}
    function RegulatorySite(ligands::Vector{Symbol},
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
        dup = findfirst(i -> ligands[i] in view(ligands, 1:(i - 1)), eachindex(ligands))
        dup === nothing || error("RegulatorySite: ligand $(ligands[dup]) appears twice " *
                                 "in one regulatory site")
        p = sortperm(ligands)
        new(ligands[p], multiplicity, allo_states[p])
    end
end

ligands(s::RegulatorySite)      = s.ligands
multiplicity(s::RegulatorySite) = s.multiplicity
allo_states(s::RegulatorySite)  = s.allo_states

"""
Whether `bound_form` is `free` with metabolite `m` added, two metabolites matching
when `same` holds: same residual, and `bound(bound_form)` equals `bound(free)` plus `m`
as a multiset. Both lists are sorted by `_met_sort_key`, so one pass decides it without
allocating. `m` is not specialized on: callers pass an element of a
`Vector{Metabolite}`, and a call whose argument type is abstract would otherwise be
dispatched at run time, which boxes both species.
"""
function _adds_metabolite(free::Species, bound_form::Species, @nospecialize(m::Metabolite),
                          same)
    f, b = bound(free), bound(bound_form)
    residual(free) == residual(bound_form) && length(b) == length(f) + 1 || return false
    j, extra = 1, false
    for x in b
        if j <= length(f) && same(x, f[j])
            j += 1
        elseif !extra && same(x, m)
            extra = true
        else
            return false
        end
    end
    return true
end

"""
Whether `bound_form` is `free` with metabolite `m` added: same residual, and
`bound(bound_form)` equals `bound(free)` plus `m` as a multiset. The
conformation may differ (a binding may change the enzyme's conformation).
"""
_binds_ligand(free::Species, bound_form::Species, m::Metabolite) =
    _adds_metabolite(free, bound_form, m, ==)

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
(E(A) + B → E(P, Q)). The constructor rejects a binding whose forms match a plain
binding by name but not by role (E + P::Mod → E(P)). It stores a step that takes up
nothing and gives off exactly one metabolite as the binding it reverses, so every
binding is stored with its metabolite consumed; every other step is oriented by the
`Mechanism` / `AllostericMechanism` constructor. See CLAUDE.md "Canonical Step Form".
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
        any(in(r), c) && error(
            "Step $(name(from_species)) → $(name(to_species)): " *
            "$(join(name.(intersect(c, r)), ", ")) both consumed and released")
        # A step that only gives off one metabolite is stored as the binding it
        # reverses (metabolite consumed), so a binding written in either direction
        # dedups to the same Step; binding constants are named after the step's
        # two sides, not its written direction. See CLAUDE.md "Canonical Step Form".
        if isempty(c) && length(r) == 1
            from_species, to_species, c, r = to_species, from_species, r, c
        end
        # A binding whose forms differ by the bound metabolite's name keeps every
        # metabolite's role: E + P::Mod → E(P) would carry the modifier copy into a
        # product form and give the product a second route into it.
        length(c) == 1 && isempty(r) &&
            _adds_metabolite(from_species, to_species, c[1],
                             (x, y) -> name(x) === name(y)) &&
            !_binds_ligand(from_species, to_species, c[1]) && error(
                "Step $(name(from_species)) → $(name(to_species)) changes a " *
                "metabolite's role: $(name(to_species)) holds $(name(from_species))'s " *
                "metabolites plus $(name(c[1])) by name, not by role; write a " *
                "modifier copy as X::Mod on both sides")
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
`_expand_change_allo_state`), `_onlya_haldane_violation`, the joint solve's priority
(`_priority_key`) and `_fused_substrate_binding` read it to tell catalysis from
binding."""
function _is_chemistry(s::Step)
    m = bound_metabolite(s)
    m === nothing || !_binds_ligand(from_species(s), to_species(s), m)
end

# Parameter family.
abstract type Parameter end

# Step-bound parameters: a rapid-equilibrium step's equilibrium constant, and a
# steady-state step's forward and reverse rate constants.
struct Krapid <: Parameter; step::Step; state::Symbol end
struct Kfor   <: Parameter; step::Step; state::Symbol end
struct Krev   <: Parameter; step::Step; state::Symbol end

"""
Regulator-site parameter: the constant of the ligand named `ligand` in the A (active) or
I (inactive) branch of the polynomial (`state`). A ligand sits at one site, so its name
and state determine the constant.
"""
struct Kreg <: Parameter
    ligand::Symbol
    state::Symbol
end

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
The names no metabolite takes: `E_total`, `Keq`, `L` and `N`, which name or are kept for
rate-equation constants, and `params`, `concs` and `v`, the local names of a rate
equation's code (`rate_equation_string`). A rate equation binds its constants, its
locals and the concentrations in one scope, so a metabolite of one of these names would
overwrite another or be overwritten.
"""
const _RESERVED_NAMES = (:E_total, :Keq, :L, :N, :params, :concs, :v)

"""
Error when a metabolite name of `names` is reserved (`_RESERVED_NAMES`) or holds `_`:
constant names join form and metabolite names with `_`, so a metabolite name holds none
and no constant takes another's name. `who` names the constructor that checks.
"""
function _check_metabolite_names(who::String, names)
    for n in names
        n in _RESERVED_NAMES &&
            error("$who: $n is a reserved name; E_total, Keq, L, N, params, concs and v " *
                  "are kept for a rate equation's constants and local names")
        occursin('_', String(n)) &&
            error("$who: the metabolite name $n holds `_`; constant names join " *
                  "metabolite names with `_`, so a metabolite name holds none")
    end
end

"""
Error unless the metabolite labels `spellings` (each a label `=>` how the DSL writes it)
are uniquely decodable: no string splits into labels two ways. Enzyme forms and their
constants join the labels of the metabolites they bind without a separator
(`_form_label`), so with `Ac`, `CoA` and `AcCoA` the forms E(Ac, CoA) and E(AcCoA) would
share a name. A reaction's labels are its substrate, product and competitive-inhibitor
names, and `Xmod`, written `X::Mod`, for a reactant `X` also declared a competitive
inhibitor; an allosteric regulator binds no form, so its name takes no part. The
Sardinas–Patterson test decides it, and the error spells both splits.
"""
function _check_uniquely_decodable(spellings::Vector{Pair{String, String}})
    spelled = Dict{String, String}()
    for (label, written) in spellings
        other = get!(spelled, label, written)
        other == written || _spelled_twice_error([other], [written])
    end
    words = sort!(collect(keys(spelled)))
    tail(w, prefix) = w[ncodeunits(prefix) + 1:end]
    # Each entry holds a dangling suffix `d` and two splits: join(long) == join(short) * d.
    queue = [(tail(v, u), [v], [u]) for u in words for v in words
             if u != v && startswith(v, u)]
    seen = Set{String}()
    while !isempty(queue)
        d, long, short = popfirst!(queue)
        d in seen && continue
        push!(seen, d)
        haskey(spelled, d) && _spelled_twice_error(getindex.(Ref(spelled), [short; d]),
                                                   getindex.(Ref(spelled), long))
        for w in words
            startswith(w, d) && push!(queue, (tail(w, d), [short; w], long))
            startswith(d, w) && push!(queue, (tail(d, w), long, [short; w]))
        end
    end
end

"""Error that the label splits `a` and `b`, each a list of DSL spellings, spell one
string."""
function _spelled_twice_error(a::Vector{String}, b::Vector{String})
    spelled = join(replace.(a, "::Mod" => "mod"))
    error("EnzymeReaction: $(join(a, " + ")) and $(join(b, " + ")) both spell $spelled: " *
          "enzyme forms and their constants join the names of the metabolites they bind " *
          "without a separator, so two forms would share a name; rename a metabolite")
end

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
- `allowed_catalytic_multiplicities::Vector{Int}` — the oligomeric states, one
  or more, that the allosteric enumerator may assign to the catalytic core.
- `shared_catalytic_site::Vector{Tuple{Symbol,Symbol}}` — `(substrate,
  product)` pairs that bind the same catalytic site, and so are never both
  bound to it at once.

Construct one with the [`@enzyme_reaction`](@ref) DSL. Each reactant takes
an atom bracket (`S[C]`, `A[C1H1]`); the brackets are load-bearing for
ping-pong and multi-substrate reactions. Reactants and regulators are sorted
by name in the constructor, a name's competitive inhibitor before its
allosteric regulator, so two equivalent declarations compare equal under
`==`/`hash`. The constructor rejects a name listed as both a substrate
and a product: concentrations and constants are keyed by name. For the same
reason it rejects a reactant or regulator named `E_total`, `Keq`, `L`, `N`,
`params`, `concs` or `v`, names a rate equation keeps for its constants and the
local names of its code, and a name holding `_`, since constant names join
names with `_`. It also rejects names that spell one string two ways, such as
`Ac`, `CoA` and `AcCoA`: an enzyme form's name joins the names of the
metabolites it binds without a separator, so E with `Ac` and `CoA` bound and E
with `AcCoA` bound would share a name. A reaction prints as the `@enzyme_reaction`
block that builds it.

```jldoctest
julia> using EnzymeRates

julia> rxn = @enzyme_reaction begin
           substrates: S[C]
           products:   P[C]
       end;

julia> rxn isa EnzymeReaction
true

julia> rxn
@enzyme_reaction begin
    substrates: S[C]
    products: P[C]
end
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
        sorted_regulators = sort(regulators;
            by = rm -> (name(regulator(rm)), regulator(rm) isa AllostericRegulator))
        sorted_mults = sort(unique(allowed_catalytic_multiplicities))
        !isempty(sorted_mults) && all(m -> m ≥ 1, sorted_mults) ||
            error("EnzymeReaction: allowed_catalytic_multiplicities must be a " *
                  "nonempty list of Ints ≥ 1, got $allowed_catalytic_multiplicities")

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
        both_roles = intersect(sub_names, prod_names)
        isempty(both_roles) || error(
            "EnzymeReaction: $(join(both_roles, ", ")) named as both a substrate and " *
            "a product; concentrations and constants are keyed by name")
        _check_metabolite_names("EnzymeReaction", [sub_names; prod_names; first.(reg_keys)])
        reactant_names = [sub_names; prod_names]
        _check_uniquely_decodable(Pair{String, String}[
            [String(n) => String(n) for n in reactant_names];
            [n in reactant_names ? "$(n)mod" => "$n::Mod" : String(n) => String(n)
             for (n, kind) in reg_keys if kind === CompetitiveInhibitor]])

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

"""
Print `r` as the `@enzyme_reaction` text that builds it: each reactant with its atoms as
one formula, the competitive inhibitors and allosteric regulators, the catalytic
multiplicities unless they are the default 1, and the shared catalytic sites. A
regulator whose multiplicities are the catalytic ones prints bare, as a bare name takes
them; any other prints with its own, as in `I(1, 2)`.
"""
function Base.show(io::IO, r::EnzymeReaction)
    mults = allowed_catalytic_multiplicities(r)
    formula(ra) = join("$e" * (c == 1 ? "" : "$c") for (e, c) in atoms(ra))
    tag(rm) = reg_type(rm) === :activator ? "::Activator" :
              reg_type(rm) === :inhibitor ? "::Inhibitor" : ""
    entry(rm) = "$(name(regulator(rm)))" * (allowed_multiplicities(rm) == mults ? "" :
                 "(" * join(allowed_multiplicities(rm), ", ") * ")") * tag(rm)
    line(label, entries) = isempty(entries) ||
                           print(io, "\n    ", label, ": ", join(entries, ", "))
    print(io, "@enzyme_reaction begin")
    for (label, T) in (("substrates", Substrate), ("products", Product))
        line(label, ["$(name(metabolite(ra)))[$(formula(ra))]"
                     for ra in reactants(r) if metabolite(ra) isa T])
    end
    for (label, T) in (("competitive_inhibitors", CompetitiveInhibitor),
                       ("allosteric_regulators", AllostericRegulator))
        line(label, [entry(rm) for rm in regulators(r) if regulator(rm) isa T])
    end
    mults == [1] || print(io, "\n    oligomeric_state: ",
                          length(mults) == 1 ? only(mults) : "(" * join(mults, ", ") * ")")
    line("shared_catalytic_site", ["($s, $p)" for (s, p) in shared_catalytic_site(r)])
    print(io, "\nend")
end

"""Whether `n` names a substrate or product of `rxn`."""
_is_reactant(rxn::EnzymeReaction, n::Symbol) =
    any(ra -> name(metabolite(ra)) === n, reactants(rxn))

# ─── Kinetic groups of Steps ───────────────────────────────────────────
#
# Code that reads a mechanism through its Step views (`steps`) finds the groups' forms
# and rapid-equilibrium segments with these.

"""
The kind of a step for kinetic grouping: the metabolites it takes up and gives off,
`(consumed, released)`. Steps that share a kinetic group share their constants, so
they must take up and give off the same metabolites: two bindings of `m` share a
kind whether or not either runs chemistry, and every isomerization is one kind.
"""
_step_kind(s::Step) = (consumed(s), released(s))

"""The distinct enzyme forms of the steps of `groups`, in the canonical form order
(`_form_index_key`: free E first, then by content)."""
function _forms(groups)
    forms = unique(sp for g in groups for s in g for sp in (from_species(s), to_species(s)))
    forms[sortperm(_form_index_key.(forms))]
end

"""
Rapid-equilibrium segments of `steps` with each form's metabolite exponents
cleared to the lowest in its segment: `(species, segments, extras, idx, segment_of)`,
where `species` is `_forms(steps)`, `segments[k]` indexes `species` in
breadth-first order from its lowest-index form, `extras[i][x]` is how many more
`x` form `i` carries than the lowest form of its segment, `idx[sp]` is the
position of form `sp` in `species` and `segment_of[i]` is the segment holding
form `i`. Segments are numbered by their lowest-index form. Within a segment,
rapid equilibrium fixes each form's weight relative to any other as a monomial
in concentrations, one factor per metabolite bound or released along the RE
path between them; chemistry steps pass the exponents through unchanged.
"""
function _re_segment_extras(steps::Vector{Vector{Step}})
    species = _forms(steps)
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
    segment_of = zeros(Int, length(species))
    for root in eachindex(species)
        segment_of[root] == 0 || continue
        k = length(segments) + 1
        expo[root] = Dict{Symbol, Int}()
        segment_of[root] = k
        segment = [root]
        queue = [root]
        while !isempty(queue)
            u = popfirst!(queue)
            for (v, d) in adj[u]
                segment_of[v] == 0 || continue
                expo[v] = mergewith(+, expo[u], d)
                segment_of[v] = k
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
    species, segments, extras, idx, segment_of
end

# ─── Raw input: a mechanism as the DSL writes it ─────────────────────
#
# The mechanism macros read a block into these plain values (`_raw_mechanism`), before
# any constructor orients, sorts or checks it; `Mechanism(raw)` and
# `AllostericMechanism(raw)` build the mechanism they write.

"""
A metabolite as a form holds it or a step exchanges it: its name, and whether it binds
as a modifier (a name under `modifiers:` written bare, or a reactant's copy `A::Mod`). A
reactant and its modifier copy are two different `Met`s.
"""
struct Met
    name::Symbol
    modifier::Bool
end
Base.isless(x::Met, y::Met) = isless((x.name, x.modifier), (y.name, y.modifier))
Base.hash(x::Met, h::UInt) = hash(x.modifier, hash(x.name, hash(:Met, h)))

"""An enzyme form as written: its conformation label, its bound metabolites and its
residual (the substrates it adds, the products it subtracts). The constructor sorts each
list, so two writings of one form are `==`."""
struct RawForm
    conformation::Symbol
    bound::Vector{Met}
    added::Vector{Symbol}
    subtracted::Vector{Symbol}
    RawForm(conformation::Symbol, bound::Vector{Met},
            added::Vector{Symbol} = Symbol[], subtracted::Vector{Symbol} = Symbol[]) =
        new(conformation, sort(bound), sort(added), sort(subtracted))
end

"""A step as written, in its written direction: going `from → to` it takes up `consumed`
and gives off `released`; `rapid_equilibrium` marks `⇌`. The constructor sorts both
lists; it neither turns the step nor checks it."""
struct RawStep
    from::RawForm
    to::RawForm
    consumed::Vector{Met}
    released::Vector{Met}
    rapid_equilibrium::Bool
    RawStep(from::RawForm, to::RawForm, consumed::Vector{Met}, released::Vector{Met},
            rapid_equilibrium::Bool) =
        new(from, to, sort(consumed), sort(released), rapid_equilibrium)
end

"""
A plain mechanism as written: the declared substrates, products and modifiers (every
name under `modifiers:`, a reactant that binds as a modifier included), the steps in
written order and direction, and `group[j]`, the kinetic group label of step `j`: labels
run 1, 2, … in order of first appearance, one per written statement (a step or a
parenthesized group).
"""
struct RawMechanism
    substrates::Vector{Symbol}
    products::Vector{Symbol}
    modifiers::Vector{Symbol}
    steps::Vector{RawStep}
    group::Vector{Int}
end

"""A regulatory site as written: its ligands, its multiplicity and each ligand's tag."""
struct RawSite
    ligands::Vector{Symbol}
    multiplicity::Int
    allo_states::Vector{Symbol}
end

"""An allosteric mechanism as written: its catalytic part, the tag of each written
group label (`cat_allo_states[l]`), the catalytic multiplicity (1 when the DSL omits it)
and the regulatory sites."""
struct RawAllostericMechanism
    catalytic::RawMechanism
    cat_allo_states::Vector{Symbol}
    catalytic_multiplicity::Int
    sites::Vector{RawSite}
end

# ─── The reaction graph ───────────────────────────────────────────────

"""
Copies of each metabolite-table entry, four bits per entry: entry `i` holds bits
`4(i - 1)` to `4i - 1`. A mechanism therefore holds at most 32 table entries and at most
15 copies of one entry on a form or a step side.
"""
const Counts = UInt128

"""How many copies of table entry `i` `c` counts."""
copies(c::Counts, i::Int) = Int((c >> (4 * (i - 1))) & 0xf)

"""The number of copies `c` counts over every entry."""
function total(c::Counts)
    n = 0
    while c != 0
        n += Int(c & 0xf)
        c >>= 4
    end
    n
end

"""The table entries `c` counts, ascending, each as often as counted."""
function entries(c::Counts)
    out = Int[]
    i = 1
    while c != 0
        for _ in 1:Int(c & 0xf)
            push!(out, i)
        end
        c >>= 4
        i += 1
    end
    out
end

"""
An enzyme form: its conformation label (an index into the mechanism's sorted
`conformations`), the table entries it binds, and its residual: the substrates it adds
and the products it subtracts.
"""
struct Form
    conformation::Int
    bound::Counts
    added::Counts
    subtracted::Counts
end

"""
A step between forms `from` and `to`, stored in its physical-forward direction: going
`from → to` it takes up `consumed` and gives off `released`. `rapid_equilibrium` marks
a rapid-equilibrium step.
"""
struct Edge
    from::Int
    to::Int
    consumed::Counts
    released::Counts
    rapid_equilibrium::Bool
end

"""
    Mechanism

A mechanism as its reaction graph, in canonical form: two writings of one mechanism give
equal fields. `substrates`, `products` and `conformations` are sorted; `table` holds
each metabolite entry once, sorted by (name, modifier), a reactant and its modifier copy
being two entries; `forms` hold free E first, then the other forms in content order;
`edges` run physical-forward, in content order; `group[k]` is edge `k`'s kinetic group,
the groups numbered by their first edge. `Mechanism(raw::RawMechanism)` builds one and
checks it; the field constructor checks nothing and serves canonical form and the Sig
decoder. Lift it to its compiled type with `EnzymeRates.compile_mechanism(m)` or
`EnzymeMechanism(m)`.
"""
struct Mechanism
    substrates::Vector{Symbol}
    products::Vector{Symbol}
    table::Vector{Met}
    conformations::Vector{Symbol}
    forms::Vector{Form}
    edges::Vector{Edge}
    group::Vector{Int}
end

"""
    AllostericMechanism

An MWC enzyme: its catalytic `Mechanism`, one allosteric tag per kinetic group of it
(`cat_allo_states[g]`), the catalytic multiplicity and the regulatory sites, sorted. A
future model with its own derivation wraps the same `Mechanism` in a type of its own.
"""
struct AllostericMechanism
    catalytic::Mechanism
    cat_allo_states::Vector{Symbol}
    catalytic_multiplicity::Int
    regulatory_sites::Vector{RegulatorySite}
    function AllostericMechanism(catalytic::Mechanism, cat_allo_states::Vector{Symbol},
                                 catalytic_multiplicity::Int,
                                 regulatory_sites::Vector{RegulatorySite})
        _check_tags(cat_allo_states, ngroups(catalytic), catalytic_multiplicity)
        sites = sort(regulatory_sites; by = _site_key)
        _check_ligands(sites)
        new(catalytic, copy(cat_allo_states), catalytic_multiplicity, sites)
    end
end

"""A site's place in the order of a mechanism's sites: its ligands, its multiplicity,
then its ligands' tags."""
_site_key(s::RegulatorySite) = (s.ligands, s.multiplicity, s.allo_states)

# Each catalytic kinetic group carries an allosteric-state tag (`:OnlyA`, `:EqualAI` or
# `:NonequalAI`; `:OnlyI` is rejected by the active-state convention); a site ligand
# may also be `:OnlyI`.
const _VALID_CAT_ALLO_STATES = (:OnlyA, :EqualAI, :NonequalAI)
const _VALID_REG_ALLO_STATES = (:OnlyA, :OnlyI, :EqualAI, :NonequalAI)

abstract type AbstractEnzymeMechanism end

"""
    EnzymeMechanism{Sig}

A plain enzyme mechanism as a singleton type, the form `rate_equation` and the fitter
take. `Sig`, its only parameter, is the mechanism encoded as one tuple of integers
(`_sig_of`), so equal mechanisms share one type and the derivation compiles once per
mechanism. `@enzyme_mechanism` returns one; `EnzymeMechanism(m)` or
`EnzymeRates.compile_mechanism(m)` lifts a `Mechanism` to it and
`EnzymeRates.Mechanism(em)` reads the mechanism back. The layout of `Sig` is internal.
"""
struct EnzymeMechanism{Sig} <: AbstractEnzymeMechanism end

"""
    AllostericEnzymeMechanism{Sig}

An allosteric (MWC) enzyme mechanism as a singleton type: `Sig`, its only parameter,
is the mechanism encoded as one tuple of integers (`_sig_of`), its catalytic graph
first. `@allosteric_mechanism` returns one; `AllostericEnzymeMechanism(am)` or
`EnzymeRates.compile_mechanism(am)` lifts an `AllostericMechanism` to it and
`EnzymeRates.AllostericMechanism(aem)` reads the mechanism back.
"""
struct AllostericEnzymeMechanism{Sig} <: AbstractEnzymeMechanism end

# ─── Content order ──────────────────────────────────────────────────────

"""
Whether the entries `a` counts, listed ascending with repeats, come before those of `b`
in lexicographic order, a list coming before every longer list it begins. With the table
sorted by (name, modifier), this is the order of the two metabolite lists by name, then
by modifier.
"""
function _seq_less(a::Counts, b::Counts)
    a == b && return false
    i = trailing_zeros(a ⊻ b) >> 2           # first entry, from 0, whose copies differ
    ca, cb = (a >> (4i)) & 0xf, (b >> (4i)) & 0xf
    # Both lists go on with min(ca, cb) copies of entry i. The one with more copies goes
    # on with entry i; the other goes on with a later entry, or ends.
    ca > cb ? (b >> (4i + 4)) != 0 : (a >> (4i + 4)) == 0
end

"""Forms by label, then bound entries, residual additions and residual subtractions,
each entry list in lexicographic order (`_seq_less`)."""
function Base.isless(f::Form, g::Form)
    f.conformation != g.conformation && return f.conformation < g.conformation
    f.bound != g.bound && return _seq_less(f.bound, g.bound)
    f.added != g.added && return _seq_less(f.added, g.added)
    _seq_less(f.subtracted, g.subtracted)
end

"""Edges by form indices, then consumed and released entries as lists (`_seq_less`), then
the flag, steady state first. Tuples and `argmin` over edges compare by it, never by the
packed numbers."""
function Base.isless(e::Edge, f::Edge)
    e.from != f.from && return e.from < f.from
    e.to != f.to && return e.to < f.to
    e.consumed != f.consumed && return _seq_less(e.consumed, f.consumed)
    e.released != f.released && return _seq_less(e.released, f.released)
    e.rapid_equilibrium < f.rapid_equilibrium
end

"""`e` stored in the opposite direction: forms swapped, consumed and released swapped."""
_reversed(e::Edge) = Edge(e.to, e.from, e.released, e.consumed, e.rapid_equilibrium)

# ─── Canonical form ─────────────────────────────────────────────────────

"""
    _canonicalize(raw::RawMechanism)
        -> (; mechanism, edge_of, flipped, group_of, form_of, misaligned)

The canonical `Mechanism` that `raw` writes, with no representation check. `edge_of[j]`
is the edge written step `j` became, `flipped[j]` whether that edge runs against step
`j` as written, `group_of[l]` the kinetic group written label `l` became,
`form_of[2j - 1]` and `form_of[2j]` the forms step `j`'s two written ends became, and
`misaligned` the lead and member of a tied kinetic group that no orientation aligns, or
no edge. Alignment stops at that group, so a result with a misaligned pair is not
canonical. Errors only where no encoding exists: more than 32 table entries or more than
15 copies of one entry on a form or a step side. The checked constructors build on it.
"""
function _canonicalize(raw::RawMechanism)
    substrates, products = sort(raw.substrates), sort(raw.products)
    table = _metabolite_table(raw)
    conformations = sort!(unique!([f.conformation for s in raw.steps
                                   for f in (s.from, s.to)]))
    n = length(raw.steps)
    ends = Vector{Form}(undef, 2n)
    for (j, s) in enumerate(raw.steps)
        ends[2j - 1] = _form(table, conformations, s.from)
        ends[2j] = _form(table, conformations, s.to)
    end
    forms, form_of = _sorted_forms(ends, conformations)
    written = [Edge(form_of[2j - 1], form_of[2j], _counts(table, s.consumed),
                    _counts(table, s.released), s.rapid_equilibrium)
               for (j, s) in enumerate(raw.steps)]
    substrate_mask, product_mask = _role_masks(table, substrates)
    fwd, misaligned = _orientation(forms, written, raw.group, substrate_mask, product_mask)
    stored = [fwd[j] ? written[j] : _reversed(written[j]) for j in 1:n]
    p = sortperm(stored)
    group, group_of = _number_groups(raw.group[p], maximum(raw.group; init = 0))
    # A search keeps every mechanism it emits, so the table, the labels and the forms,
    # each built in a vector grown past its length, are stored as copies at their length.
    m = Mechanism(substrates, products, copy(table), copy(conformations), copy(forms),
                  stored[p], group)
    (; mechanism = m, edge_of = invperm(p), flipped = .!fwd, group_of, form_of,
     misaligned)
end

"""Every metabolite entry `raw` holds, sorted and once each: bound, consumed and released
metabolites, and the residual names as non-modifier entries."""
function _metabolite_table(raw::RawMechanism)
    seen = Met[]
    for s in raw.steps
        append!(seen, s.consumed, s.released)
        for f in (s.from, s.to)
            append!(seen, f.bound)
            foreach(x -> push!(seen, Met(x, false)), f.added)
            foreach(x -> push!(seen, Met(x, false)), f.subtracted)
        end
    end
    table = unique!(sort!(seen))
    length(table) <= 32 || error("Mechanism: a mechanism holds at most 32 metabolite " *
        "entries (a name and its modifier copy count apart); this one has " *
        "$(length(table))")
    table
end

"""The counts of `mets`, entries of the sorted `table`."""
function _counts(table::Vector{Met}, mets)
    c = Counts(0)
    for x in mets
        i = searchsortedfirst(table, x)
        copies(c, i) == 15 && error("Mechanism: a form or a step side holds $(x.name) " *
            "more than 15 times; a mechanism holds at most 15 copies of one metabolite " *
            "there")
        c += Counts(1) << (4 * (i - 1))
    end
    c
end

"""Form `f` as entries of the sorted `table` and `conformations`."""
_form(table::Vector{Met}, conformations::Vector{Symbol}, f::RawForm) =
    Form(searchsortedfirst(conformations, f.conformation), _counts(table, f.bound),
         _counts(table, (Met(x, false) for x in f.added)),
         _counts(table, (Met(x, false) for x in f.subtracted)))

"""The distinct forms of `ends` in canonical order, free E (label `E`, nothing bound, no
residual) first, and the index each end became."""
function _sorted_forms(ends::Vector{Form}, conformations::Vector{Symbol})
    free_E = Form(something(findfirst(==(:E), conformations), 0), 0, 0, 0)
    forms, form_of = Form[], Vector{Int}(undef, length(ends))
    for t in sortperm(ends; by = f -> (f != free_E, f))
        (isempty(forms) || forms[end] != ends[t]) && push!(forms, ends[t])
        form_of[t] = length(forms)
    end
    forms, form_of
end

"""Masks of the table's substrate entries and of its product entries; a modifier entry
is in neither, so it counts on neither side of a reaction."""
function _role_masks(table::Vector{Met}, substrates::Vector{Symbol})
    substrate_mask = product_mask = Counts(0)
    for (i, x) in enumerate(table)
        x.modifier && continue
        nibble = Counts(0xf) << (4 * (i - 1))
        insorted(x.name, substrates) ? (substrate_mask |= nibble) :
                                       (product_mask |= nibble)
    end
    substrate_mask, product_mask
end

"""
    _orientation(forms, written, labels, substrate_mask, product_mask)
        -> (fwd::BitVector, misaligned::Vector{Edge})

Whether each written edge runs in its physical-forward direction. A step with one
metabolite in all takes it up (tier 0). Otherwise the side holding more substrates,
bound plus free, then fewer products, comes first (tier 1); otherwise the step runs from
a form where only products enter or leave to one where only substrates do (tier 2);
otherwise from the earlier form (tier 3). A modifier entry counts on neither side. Entry
sides are read over every edge as written; reversing an edge swaps its ends and lists
together, so they do not depend on any written direction. Tied groups then turn as one
(`_align_tied_groups!`).
"""
function _orientation(forms::Vector{Form}, written::Vector{Edge}, labels::Vector{Int},
                      substrate_mask::Counts, product_mask::Counts)
    side(c::Counts) = UInt8((c & substrate_mask != 0) + 2 * (c & product_mask != 0))
    enters = zeros(UInt8, length(forms))
    for e in written
        enters[e.from] |= side(e.consumed)
        enters[e.to] |= side(e.released)
    end
    score(f::Int, free::Counts) =
        (total(forms[f].bound & substrate_mask) + total(free & substrate_mask),
         -total(forms[f].bound & product_mask) - total(free & product_mask))
    n = length(written)
    fwd, tier = trues(n), zeros(Int, n)
    for (k, e) in enumerate(written)
        if total(e.consumed) + total(e.released) == 1
            fwd[k] = e.consumed != 0
        elseif (a = score(e.from, e.consumed)) != (b = score(e.to, e.released))
            fwd[k], tier[k] = a > b, 1
        elseif (enters[e.from], enters[e.to]) in ((0x02, 0x01), (0x01, 0x02))
            fwd[k], tier[k] = enters[e.from] == 0x02, 2
        else
            fwd[k], tier[k] = e.from < e.to, 3
        end
    end
    misaligned = _align_tied_groups!(fwd, tier, forms, written, labels)
    fwd, misaligned
end

"""
Turn the steps tiers 2 and 3 decide in each kinetic group (`labels`) with the group's
lead: the tied step that changes conformation first, then by tier, then by its oriented
edge. A group's steps share their constants, so a conformational change and its mirror
at other bound forms must run the same way. Every other tied step that takes up and
gives off the lead's metabolites takes the orientation that does, and, if it changes
conformation, makes the lead's change; a step that fits both ways keeps its own. Returns
the lead and the first step that fits the lead's metabolites but not its change, or no
edge.
"""
function _align_tied_groups!(fwd::BitVector, tier::Vector{Int}, forms::Vector{Form},
                             written::Vector{Edge}, labels::Vector{Int})
    label_pair(e::Edge) = (forms[e.from].conformation, forms[e.to].conformation)
    turns(e::Edge) = label_pair(e)[1] != label_pair(e)[2]
    oriented(k::Int) = fwd[k] ? written[k] : _reversed(written[k])
    for l in unique(labels)
        tied = [k for k in eachindex(written) if labels[k] == l && tier[k] > 1]
        length(tied) < 2 && continue
        lead = oriented(argmin(k -> (!turns(oriented(k)), tier[k], oriented(k)), tied))
        for k in tied
            o = oriented(k)
            same = [x for x in (o, _reversed(o))
                    if (x.consumed, x.released) == (lead.consumed, lead.released)]
            isempty(same) && continue
            i = findfirst(x -> !turns(x) || label_pair(x) == label_pair(lead), same)
            i === nothing && return [lead, same[1]]
            fwd[k] = same[i] == written[k]
        end
    end
    Edge[]
end

"""Group labels in edge order renumbered by first occurrence, and the number each
written label took."""
function _number_groups(labels::Vector{Int}, G::Int)
    group_of, group = zeros(Int, G), similar(labels)
    next = 0
    for (k, l) in enumerate(labels)
        if group_of[l] == 0
            next += 1
            group_of[l] = next
        end
        group[k] = group_of[l]
    end
    group, group_of
end

# ─── Constructors and their checks ──────────────────────────────────────

"""
    Mechanism(raw::RawMechanism)

The mechanism `raw` writes, in canonical form, checked: its declarations
(`_check_declarations`), then its representation (`_check_representation`). Whether its
rate law is defined is the boundary check's question, not the constructor's.
"""
function Mechanism(raw::RawMechanism)
    _check_declarations(raw)
    c = _canonicalize(raw)
    _check_representation(c.mechanism, c.misaligned, raw.modifiers)
    c.mechanism
end

"""
    AllostericMechanism(raw::RawAllostericMechanism)

The allosteric mechanism `raw` writes: its regulatory sites, its tags (each written
group label's), then its catalytic part as `Mechanism(raw.catalytic)` checks it, each
tag moved to the canonical group its label became.
"""
function AllostericMechanism(raw::RawAllostericMechanism)
    sites = [RegulatorySite(s.ligands, s.multiplicity, s.allo_states) for s in raw.sites]
    _check_tags(raw.cat_allo_states, maximum(raw.catalytic.group; init = 0),
                raw.catalytic_multiplicity)
    _check_declarations(raw.catalytic)
    c = _canonicalize(raw.catalytic)
    _check_representation(c.mechanism, c.misaligned, raw.catalytic.modifiers)
    tags = similar(raw.cat_allo_states)
    for (l, g) in enumerate(c.group_of)
        tags[g] = raw.cat_allo_states[l]
    end
    AllostericMechanism(c.mechanism, tags, raw.catalytic_multiplicity, sites)
end

# ─── Checks ─────────────────────────────────────────────────────────────

"""Error unless `cat_allo_states` holds one valid tag for each of `G` kinetic groups and
the catalytic multiplicity is at least 1. No catalytic group is `:OnlyI`: the active
conformation is the one that catalyzes."""
function _check_tags(cat_allo_states::Vector{Symbol}, G::Int, catalytic_multiplicity::Int)
    length(cat_allo_states) == G ||
        error("AllostericMechanism: cat_allo_states length $(length(cat_allo_states)) " *
              "must match the number of kinetic groups $G")
    catalytic_multiplicity ≥ 1 ||
        error("AllostericMechanism: catalytic_multiplicity must be ≥ 1, " *
              "got $catalytic_multiplicity")
    bad = findfirst(∉(_VALID_CAT_ALLO_STATES), cat_allo_states)
    bad === nothing ||
        error("AllostericMechanism: catalytic group $bad has invalid " *
              "allo state $(cat_allo_states[bad]) (must be one of " *
              "$_VALID_CAT_ALLO_STATES); :OnlyI is rejected for " *
              "catalytic groups (active-state-active convention)")
    nothing
end

"""Error when a site ligand's name is reserved or holds `_`, or when a ligand sits on two
sites: each ligand's constant is named after the ligand alone."""
function _check_ligands(sites::Vector{RegulatorySite})
    ligs = [l for s in sites for l in s.ligands]
    _check_metabolite_names("AllostericMechanism", ligs)
    dup = findfirst(i -> ligs[i] in view(ligs, 1:(i - 1)), eachindex(ligs))
    dup === nothing ||
        error("AllostericMechanism: ligand $(ligs[dup]) appears in two " *
              "distinct regulatory sites; rendered Kreg names would " *
              "collide. Same-ligand-two-sites is not enumerated.")
    nothing
end

"""Error unless `raw` declares substrates and products, each name once, no name both a
substrate and a product, and gives each step one kinetic group label from 1 to the
number of groups, every label used."""
function _check_declarations(raw::RawMechanism)
    isempty(raw.substrates) && error("Mechanism: substrates must not be empty")
    isempty(raw.products) && error("Mechanism: products must not be empty")
    for (names, kind) in ((raw.substrates, "substrate"), (raw.products, "product"),
                          (raw.modifiers, "modifier"))
        allunique(names) || error("Mechanism: duplicate $kind names")
    end
    both = findfirst(in(raw.products), raw.substrates)
    both === nothing ||
        error("Mechanism: $(raw.substrates[both]) named as both a substrate and a " *
              "product; concentrations and constants are keyed by name")
    length(raw.group) == length(raw.steps) && all(≥(1), raw.group) ||
        error("Mechanism: each step takes one kinetic group label from 1 to the number " *
              "of groups")
    empty = findfirst(∉(raw.group), 1:maximum(raw.group; init = 0))
    empty === nothing ||
        error("Mechanism: kinetic group $empty is empty; a kinetic group holds at least " *
              "one step")
    nothing
end

"""
Error unless the canonical mechanism `m` represents one: every name declared (in
`modifiers`, the declared modifiers, or among its substrates and products), a modifier
entry only for a declared modifier and a plain entry only for a reactant, residuals that
add substrates and subtract products; no step joining a form to itself, taking up and
giving off one metabolite, or binding a metabolite into a form that holds it in another
role; tied steps that make one conformational change (`misaligned`, from
`_canonicalize`); no two forms that render one name; each kinetic group one kind of step
with one flag; each reaction once.
"""
function _check_representation(m::Mechanism, misaligned::Vector{Edge},
                               modifiers::Vector{Symbol})
    for x in m.table
        _is_reactant(m, x.name) || x.name in modifiers ||
            error("Mechanism: $(x.name) is not declared; declare it under substrates:, " *
                  "products: or modifiers:")
        x.modifier && !(x.name in modifiers) &&
            error("Mechanism: $(x.name) binds as a modifier, but modifiers: does not " *
                  "declare it")
        x.modifier || _is_reactant(m, x.name) ||
            error("Mechanism: $(x.name) is a modifier, so it binds only as a modifier")
    end
    substrate_mask, product_mask = _role_masks(m)
    for f in m.forms, (bad, verb) in ((f.added & product_mask, "adds"),
                                      (f.subtracted & substrate_mask, "subtracts"))
        bad == 0 ||
            error("Mechanism: a residual adds substrates and subtracts products; one " *
                  "$verb $(m.table[first(entries(bad))].name)")
    end
    for e in m.edges
        e.from == e.to &&
            error("Step: both ends are the form $(_form_name(m, e.from)); a step joins " *
                  "two forms")
        shared = _occupied(e.consumed) & _occupied(e.released)
        shared == 0 ||
            error("Step $(_form_name(m, e.from)) → $(_form_name(m, e.to)): " *
                  "$(join([m.table[i].name for i in entries(shared)], ", ")) both " *
                  "consumed and released")
        _binds_by_name(m, e) && !_is_plain_binding(m, e) &&
            error("Step $(_form_name(m, e.from)) → $(_form_name(m, e.to)) changes a " *
                  "metabolite's role: $(_form_name(m, e.to)) holds " *
                  "$(_form_name(m, e.from))'s metabolites plus " *
                  "$(m.table[_bound_entry(e)].name) by name, not by role; write a " *
                  "modifier copy as X::Mod on both sides")
    end
    _check_aligned(m, misaligned)
    _check_distinct_form_names(m)
    groups = _group_members(m)
    _check_uniform_groups(m, groups)
    _check_each_reaction_once(m, groups)
end

"""Error unless `misaligned`, from `_canonicalize`, is empty: it holds the lead and a member
of a tied kinetic group that no orientation aligns."""
function _check_aligned(m::Mechanism, misaligned::Vector{Edge})
    isempty(misaligned) ||
        error("Mechanism: a kinetic group holds $(join(_sides(m, misaligned[1]), " → ")) " *
              "and $(join(_sides(m, misaligned[2]), " → ")), which make different " *
              "conformational changes; the steps of a kinetic group share their " *
              "constants, so they must make the same conformational change")
    nothing
end

"""A count of 1 for every entry `c` counts at least once."""
function _occupied(c::Counts)
    x = c | (c >> 1)
    x |= x >> 2
    x & Counts(0x11111111111111111111111111111111)
end

"""Error when two forms of `m` render one name: every constant is named after the forms
its step joins (`_form_name`)."""
function _check_distinct_form_names(m::Mechanism)
    seen = Dict{String, Int}()
    for f in eachindex(m.forms)
        other = get!(seen, _form_name(m, f), f)
        other == f ||
            error("Mechanism: the enzyme forms $(_form_text(m, other)) and " *
                  "$(_form_text(m, f)) both render the name $(_form_name(m, f)); every " *
                  "constant is named after the forms its step joins, so rename a " *
                  "metabolite to tell the two forms apart")
    end
end

"""Error unless the edges of each kinetic group (`groups`) take up and give off the same
entries and carry one flag: they share their constants."""
function _check_uniform_groups(m::Mechanism, groups::Vector{Vector{Int}})
    label(e::Edge) = join(_sides(m, e), " → ") * (e.rapid_equilibrium ? " (RE)" : " (SS)")
    for members in groups
        e1 = m.edges[first(members)]
        for k in members
            e = m.edges[k]
            (e.consumed, e.released, e.rapid_equilibrium) ==
                (e1.consumed, e1.released, e1.rapid_equilibrium) && continue
            error("Mechanism: a kinetic group holds $(label(e1)) and $(label(e)), " *
                  "which differ in kind or in RE/SS; the steps of a kinetic group " *
                  "share their constants, so they must be the same kind of step " *
                  "with the same flag")
        end
    end
end

"""Error when one reaction, the unordered pair of an edge's two sides, appears in two
edges: a reaction belongs to one kinetic group, whose constants are named after it, and
one group holding it twice would count its edge twice. Reads the groups (`groups`) in
order, each in edge order, and names the reaction as the edge seen first stores it."""
function _check_each_reaction_once(m::Mechanism, groups::Vector{Vector{Int}})
    seen = Dict{NTuple{2, Tuple{Int, Counts}}, Tuple{Int, Int}}()
    flag(e::Edge) = e.rapid_equilibrium ? "RE" : "SS"
    for (g, members) in enumerate(groups), k in members
        e = m.edges[k]
        key = minmax((e.from, e.consumed), (e.to, e.released))
        if haskey(seen, key)
            g0, k0 = seen[key]
            text = join(_sides(m, k0), " ⇌ ")
            g0 == g && error("Mechanism: kinetic group $g holds the reaction $text " *
                             "twice; a mechanism holds each reaction once")
            e0 = m.edges[k0]
            flags = flag(e0) == flag(e) ? "" :
                " ($(flag(e0)) in group $g0, $(flag(e)) in group $g)"
            error("Mechanism: kinetic groups $g0 and $g both hold the reaction " *
                  "$text$flags; a reaction belongs to one kinetic group, whose " *
                  "constants are named after it")
        end
        seen[key] = (g, k)
    end
end

# ─── Identity ───────────────────────────────────────────────────────────

# `hash` by content for every value type, and `==` by content for those whose fields hold
# vectors or structs, whose default `==` would compare identity: each field in declaration
# order, seeded with the type's name. A struct's default hash mixes in its type's identity,
# which changes between precompile sessions; a graph's hash and a raw value's read only
# symbols, integers and flags, so every process agrees on them. An isbits type (`Form`,
# `Edge`) keeps Julia's `==`, which compares its fields. `Species` keeps its own pair: its
# stored `name` is derived. `Met` keeps Julia's `==`, which compares its two fields, and
# has its own `hash`.
for T in (Residual, RegulatorySite, Step, Krapid, Kfor, Krev, Kreg, ReactantAtoms,
          RegulatorMults, EnzymeReaction, RawForm, RawStep, RawMechanism, RawSite,
          RawAllostericMechanism, Form, Edge, Mechanism, AllostericMechanism)
    fs = fieldnames(T)
    isbitstype(T) || @eval Base.:(==)(a::$T, b::$T) =
        $(reduce((l, r) -> :($l && $r), [:(a.$f == b.$f) for f in fs]))
    @eval Base.hash(x::$T, h::UInt) =
        $(foldl((acc, f) -> :(hash(x.$f, $acc)), fs;
                init = :(hash($(QuoteNode(nameof(T))), h))))
end

# ─── Sig ────────────────────────────────────────────────────────────────
#
# A mechanism's Sig, every entry an `Int`:
#   plain       the numbers of substrates, products, table entries, conformation labels,
#               forms and edges; the substrate and product names; each table entry's
#               name and modifier flag; the conformation labels; each form (label,
#               bound, added, subtracted); each edge (from, to, consumed, released,
#               flag); each edge's group. A packed count takes two entries, high half
#               first.
#   name        its number of characters, then each character's code point.
#   allosteric  the plain Sig of its catalytic mechanism; one tag code per group; the
#               catalytic multiplicity; the number of sites, then per site its
#               multiplicity, its number of ligands and per ligand its name and tag code.

"""The allosteric tags, coded by position: catalytic groups and site ligands share it."""
const _TAG_CODES = (:OnlyA, :OnlyI, :EqualAI, :NonequalAI)
_tag_code(t::Symbol) = findfirst(==(t), _TAG_CODES)::Int

function _sig_name!(out::Vector{Int}, s::Symbol)
    str = String(s)
    push!(out, length(str))
    for c in str
        push!(out, Int(c))
    end
    out
end
_sig_counts!(out::Vector{Int}, c::Counts) = push!(out, (c >> 64) % Int, c % Int)

function _sig_mechanism!(out::Vector{Int}, m::Mechanism)
    push!(out, length(m.substrates), length(m.products), length(m.table),
          length(m.conformations), length(m.forms), length(m.edges))
    foreach(s -> _sig_name!(out, s), m.substrates)
    foreach(s -> _sig_name!(out, s), m.products)
    for x in m.table
        _sig_name!(out, x.name)
        push!(out, x.modifier)
    end
    foreach(s -> _sig_name!(out, s), m.conformations)
    for f in m.forms
        push!(out, f.conformation)
        _sig_counts!(out, f.bound)
        _sig_counts!(out, f.added)
        _sig_counts!(out, f.subtracted)
    end
    for e in m.edges
        push!(out, e.from, e.to)
        _sig_counts!(out, e.consumed)
        _sig_counts!(out, e.released)
        push!(out, e.rapid_equilibrium)
    end
    append!(out, m.group)
end

"""`m` as one homogeneous integer tuple, the only type parameter of its compiled type."""
_sig_of(m::Mechanism) = Tuple(_sig_mechanism!(Int[], m))
function _sig_of(am::AllostericMechanism)
    out = _sig_mechanism!(Int[], am.catalytic)
    foreach(t -> push!(out, _tag_code(t)), am.cat_allo_states)
    push!(out, am.catalytic_multiplicity, length(am.regulatory_sites))
    for s in am.regulatory_sites
        push!(out, s.multiplicity, length(s.ligands))
        for (l, t) in zip(s.ligands, s.allo_states)
            _sig_name!(out, l)
            push!(out, _tag_code(t))
        end
    end
    Tuple(out)
end

"""A Sig read from its start: a struct, not captured variables, so nothing is boxed."""
mutable struct _SigReader
    sig::Vector{Int}
    at::Int
end
_next!(r::_SigReader) = r.sig[r.at += 1]
function _next_name!(r::_SigReader)
    n = _next!(r)
    Symbol(String([Char(_next!(r)) for _ in 1:n]))
end
function _next_counts!(r::_SigReader)
    hi = _next!(r) % UInt64
    lo = _next!(r) % UInt64
    (Counts(hi) << 64) | Counts(lo)
end

function _read_mechanism!(r::_SigReader)
    ns = _next!(r); np = _next!(r); nt = _next!(r)
    nc = _next!(r); nf = _next!(r); ne = _next!(r)
    substrates = [_next_name!(r) for _ in 1:ns]
    products = [_next_name!(r) for _ in 1:np]
    table = [Met(_next_name!(r), _next!(r) == 1) for _ in 1:nt]
    conformations = [_next_name!(r) for _ in 1:nc]
    forms = [Form(_next!(r), _next_counts!(r), _next_counts!(r), _next_counts!(r))
             for _ in 1:nf]
    edges = [Edge(_next!(r), _next!(r), _next_counts!(r), _next_counts!(r),
                  _next!(r) == 1) for _ in 1:ne]
    group = [_next!(r) for _ in 1:ne]
    Mechanism(substrates, products, table, conformations, forms, edges, group)
end

function _read_site!(r::_SigReader)
    multiplicity = _next!(r)
    n = _next!(r)
    ligands, states = Symbol[], Symbol[]
    for _ in 1:n
        push!(ligands, _next_name!(r))
        push!(states, _TAG_CODES[_next!(r)])
    end
    RegulatorySite(ligands, multiplicity, states)
end

"""The mechanism whose Sig is `sig`, splatted. Decoding trusts the Sig, which only
`_sig_of` writes, so it canonicalizes and checks nothing. The varargs method compiles
once for every Sig length, as long as callers splat a Sig whose length their code leaves
open."""
_mechanism_from_sig(sig::Int...) = _read_mechanism!(_SigReader(collect(sig), 0))
"""The allosteric mechanism whose Sig is `sig`, splatted; only the composition
constructor's tag and site checks run."""
function _allosteric_mechanism_from_sig(sig::Int...)
    r = _SigReader(collect(sig), 0)
    m = _read_mechanism!(r)
    tags = [_TAG_CODES[_next!(r)] for _ in 1:ngroups(m)]
    multiplicity = _next!(r)
    sites = [_read_site!(r) for _ in 1:_next!(r)]
    AllostericMechanism(m, tags, multiplicity, sites)
end

"""Lift a `Mechanism` to its singleton `EnzymeMechanism` type (`_sig_of`)."""
EnzymeMechanism(m::Mechanism) = EnzymeMechanism{_sig_of(m)}()
"""Lift an `AllostericMechanism` to its singleton `AllostericEnzymeMechanism` type."""
AllostericEnzymeMechanism(am::AllostericMechanism) =
    AllostericEnzymeMechanism{_sig_of(am)}()

# The lifts back read the Sig at run time, so they compile once and no method specializes
# on a Sig but the derivation's own. Each gives back the mechanism its Sig encodes.
"""The `Mechanism` the compiled mechanism `em` encodes."""
Mechanism(@nospecialize(em::EnzymeMechanism)) =
    _mechanism_from_sig(typeof(em).parameters[1]...)
"""The `AllostericMechanism` the compiled mechanism `aem` encodes."""
AllostericMechanism(@nospecialize(aem::AllostericEnzymeMechanism)) =
    _allosteric_mechanism_from_sig(typeof(aem).parameters[1]...)
"""Lift a singleton mechanism to its concrete `Mechanism` or `AllostericMechanism`."""
_concrete(@nospecialize(em::EnzymeMechanism)) = Mechanism(em)
_concrete(@nospecialize(aem::AllostericEnzymeMechanism)) = AllostericMechanism(aem)
catalytic_mechanism(@nospecialize(aem::AllostericEnzymeMechanism)) =
    EnzymeMechanism(AllostericMechanism(aem).catalytic)
catalytic_multiplicity(@nospecialize(aem::AllostericEnzymeMechanism)) =
    AllostericMechanism(aem).catalytic_multiplicity

# ─── Reads ──────────────────────────────────────────────────────────────

substrates(m::Mechanism) = m.substrates
products(m::Mechanism) = m.products
"""Every name bound as a modifier, reactant copies included, sorted."""
modifiers(m::Mechanism) = unique!([x.name for x in m.table if x.modifier])
# The concentration names in the order the rate law reads them: substrates, products,
# then the modifier-only names, each by name; an allosteric mechanism appends each site's
# ligands in site order, once each.
metabolites(m::Mechanism) = unique!([m.substrates; m.products; modifiers(m)])
metabolites(am::AllostericMechanism) =
    unique!([metabolites(am.catalytic);
             [l for s in am.regulatory_sites for l in s.ligands]])
"""Whether `x` names a substrate or a product of `m`."""
_is_reactant(m::Mechanism, x::Symbol) = insorted(x, m.substrates) || insorted(x, m.products)
_role_masks(m::Mechanism) = _role_masks(m.table, m.substrates)

# An allosteric mechanism's names are its catalytic graph's, plus its site ligands.
substrates(am::AllostericMechanism) = substrates(am.catalytic)
products(am::AllostericMechanism) = products(am.catalytic)
modifiers(am::AllostericMechanism) = modifiers(am.catalytic)
_is_reactant(am::AllostericMechanism, x::Symbol) = _is_reactant(am.catalytic, x)
"""The ligands of `am`'s regulatory sites, site by site, each once."""
allosteric_regulators(am::AllostericMechanism) =
    unique(Symbol[l for s in am.regulatory_sites for l in s.ligands])

const _AnyMech = Union{Mechanism, AllostericMechanism}

catalytic(m::Mechanism) = m
catalytic(am::AllostericMechanism) = am.catalytic
cat_allo_states(am::AllostericMechanism) = am.cat_allo_states
cat_allo_state(am::AllostericMechanism, g::Int) = am.cat_allo_states[g]
catalytic_multiplicity(am::AllostericMechanism) = am.catalytic_multiplicity
regulatory_sites(am::AllostericMechanism) = am.regulatory_sites
"""The number of kinetic groups of `m`."""
ngroups(m::Mechanism) = maximum(m.group; init = 0)
"""Each kinetic group's edges, in edge order."""
function _group_members(m::Mechanism)
    members = [Int[] for _ in 1:ngroups(m)]
    for (k, g) in enumerate(m.group)
        push!(members[g], k)
    end
    members
end
_has_residual(f::Form) = f.added != 0 || f.subtracted != 0
_is_binding(e::Edge) = e.released == 0 && total(e.consumed) == 1
"""The entry a binding takes up; 0 for any other edge."""
_bound_entry(e::Edge) = _is_binding(e) ? (trailing_zeros(e.consumed) >> 2) + 1 : 0
"""The entry binding `e` takes up and the bound entries of its two forms; the entry is 0
unless `e` is a binding whose two forms carry one residual."""
function _binding_ends(m::Mechanism, e::Edge)
    f, t = m.forms[e.from], m.forms[e.to]
    i = (f.added, f.subtracted) == (t.added, t.subtracted) ? _bound_entry(e) : 0
    i, f.bound, t.bound
end
"""Whether binding `e` takes up its entry into a form that holds its source's entries plus
that one, with the same residual; the label may change."""
function _is_plain_binding(m::Mechanism, e::Edge)
    i, f, t = _binding_ends(m, e)
    i != 0 && copies(f, i) < 15 && t == f + (Counts(1) << (4 * (i - 1)))
end
"""Whether binding `e` takes up its entry into a form that holds its source's metabolites
plus that one by name, with the same residual: a plain binding, or one that binds a
reactant where its modifier copy sits, or the reverse."""
function _binds_by_name(m::Mechanism, e::Edge)
    i, f, t = _binding_ends(m, e)
    i == 0 && return false
    names(c::Counts) = [m.table[j].name for j in entries(c)]
    names(t) == sort!(push!(names(f), m.table[i].name))
end

"""`m` as raw input: every edge a step in its stored direction, labelled with its group.
Every vector is fresh, so a move that edits the raw value never edits `m`."""
function _raw(m::Mechanism)
    form(f::Int) = (x = m.forms[f];
        RawForm(m.conformations[x.conformation], [m.table[i] for i in entries(x.bound)],
                [m.table[i].name for i in entries(x.added)],
                [m.table[i].name for i in entries(x.subtracted)]))
    steps = [RawStep(form(e.from), form(e.to), [m.table[i] for i in entries(e.consumed)],
                     [m.table[i] for i in entries(e.released)], e.rapid_equilibrium)
             for e in m.edges]
    RawMechanism(copy(m.substrates), copy(m.products), modifiers(m), steps,
                 copy(m.group))
end
_raw(am::AllostericMechanism) = RawAllostericMechanism(_raw(am.catalytic),
    copy(am.cat_allo_states), am.catalytic_multiplicity,
    [RawSite(copy(s.ligands), s.multiplicity, copy(s.allo_states))
     for s in am.regulatory_sites])

# ─── Step views and Step input ─────────────────────────────────────────
#
# Code that walks a mechanism step by step reads its kinetic groups as `Step`s between
# `Species` (`steps`), and code that builds kinetic groups of `Step`s turns them into raw
# input (`_raw_of`), which the constructors canonicalize and check.

"""Table entry `i` of `m` as a `Step` or a `Species` holds it: a modifier entry as its
`CompetitiveInhibitor`, any other as the `Substrate` or `Product` of its name."""
function _view_metabolite(m::Mechanism, i::Int)
    x = m.table[i]
    x.modifier ? CompetitiveInhibitor(x.name) :
        insorted(x.name, m.substrates) ? Substrate(x.name) : Product(x.name)
end

"""Form `f` of `m` as a `Species`."""
function _view_species(m::Mechanism, f::Int)
    x = m.forms[f]
    Species(Metabolite[_view_metabolite(m, i) for i in entries(x.bound)],
            m.conformations[x.conformation],
            Residual(Substrate[Substrate(m.table[i].name) for i in entries(x.added)],
                     Product[Product(m.table[i].name) for i in entries(x.subtracted)]))
end

"""
    steps(m) -> Vector{Vector{Step}}

The kinetic groups of `m`, or of an allosteric mechanism's catalytic graph, as `Step`s:
the groups in order, each group's edges in edge order, every step in its stored
direction. Each call builds fresh vectors.
"""
function steps(m::Mechanism)
    forms = [_view_species(m, f) for f in eachindex(m.forms)]
    mets(c::Counts) = Metabolite[_view_metabolite(m, i) for i in entries(c)]
    [Step[Step(forms[e.from], forms[e.to], mets(e.consumed), mets(e.released),
               e.rapid_equilibrium) for e in m.edges[members]]
     for members in _group_members(m)]
end
steps(am::AllostericMechanism) = steps(am.catalytic)

"""
The raw input that `groups`, kinetic groups of `Step`s over the given substrates and
products, write: each step as stored, labelled with the position of its group, and every
name a step binds as a modifier declared one. Errors on an empty group.
"""
function _raw_of(substrates::Vector{Symbol}, products::Vector{Symbol},
                 groups::Vector{Vector{Step}})
    bad = findfirst(isempty, groups)
    bad === nothing ||
        error("Mechanism: kinetic group $bad is empty; a kinetic group holds at least " *
              "one step")
    met(x::Metabolite) = Met(name(x), x isa CompetitiveInhibitor)
    form(sp::Species) = RawForm(conformation(sp), Met[met(x) for x in bound(sp)],
                                Symbol[name(a) for a in added(residual(sp))],
                                Symbol[name(p) for p in subtracted(residual(sp))])
    steps = RawStep[RawStep(form(from_species(s)), form(to_species(s)),
                            Met[met(x) for x in consumed(s)],
                            Met[met(x) for x in released(s)], is_equilibrium(s))
                    for group in groups for s in group]
    modifiers = sort!(unique!(Symbol[x.name for s in steps
                                     for x in [s.consumed; s.released; s.from.bound;
                                               s.to.bound] if x.modifier]))
    RawMechanism(substrates, products, modifiers, steps,
                 [g for (g, group) in enumerate(groups) for _ in group])
end

"""The raw input that `groups` write over `rxn`'s substrates and products (`_raw_of`)."""
_raw_of(rxn::EnzymeReaction, groups::Vector{Vector{Step}}) =
    _raw_of(name.(substrates(rxn)), name.(products(rxn)), groups)

"""
    Mechanism(rxn::EnzymeReaction, groups::Vector{Vector{Step}})

The mechanism whose kinetic groups are `groups`, over `rxn`'s substrates and products:
`Mechanism(raw)` of the raw input they write (`_raw_of`).
"""
Mechanism(rxn::EnzymeReaction, groups::Vector{Vector{Step}}) =
    Mechanism(_raw_of(rxn, groups))

"""
    AllostericMechanism(rxn, groups, cat_allo_states, catalytic_multiplicity, sites)

The allosteric mechanism whose catalytic kinetic groups are `groups`, over `rxn`'s
substrates and products, `groups[g]` tagged `cat_allo_states[g]`, with the given
catalytic multiplicity and regulatory sites.
"""
AllostericMechanism(rxn::EnzymeReaction, groups::Vector{Vector{Step}},
                    cat_allo_states::Vector{Symbol}, catalytic_multiplicity::Int,
                    sites::Vector{RegulatorySite}) =
    _allosteric_of(_raw_of(rxn, groups), cat_allo_states, catalytic_multiplicity, sites)

"""The allosteric mechanism of catalytic raw input `raw`, its group label `l` tagged
`cat_allo_states[l]`, with the given catalytic multiplicity and regulatory sites."""
_allosteric_of(raw::RawMechanism, cat_allo_states::Vector{Symbol},
               catalytic_multiplicity::Int, sites::Vector{RegulatorySite}) =
    AllostericMechanism(RawAllostericMechanism(raw, cat_allo_states, catalytic_multiplicity,
        RawSite[RawSite(s.ligands, s.multiplicity, s.allo_states) for s in sites]))

"""
    _with(m; groups) / _with(am; groups, states, sites)

`m` built again with the given parts replaced: `groups`, the kinetic groups as `Step`s
over `m`'s substrates and products (`_raw_of`); for an `AllostericMechanism` also
`states`, the tags parallel to `groups`, and `sites`, the regulatory sites. Parts not
given keep `m`'s; an allosteric mechanism given no groups keeps its catalytic graph.
Expansion moves use it to change one part of a mechanism. The child shares with `m`
only values no code changes: the catalytic graph it keeps and the regulatory sites.
"""
_with(m::Mechanism; groups = steps(m)) =
    Mechanism(_raw_of(substrates(m), products(m), groups))
function _with(am::AllostericMechanism; groups = nothing, states = cat_allo_states(am),
               sites = regulatory_sites(am))
    groups === nothing &&
        return AllostericMechanism(am.catalytic, states, catalytic_multiplicity(am), sites)
    _allosteric_of(_raw_of(substrates(am), products(am), groups), states,
                   catalytic_multiplicity(am), sites)
end

"""
    _canonical_order(rxn, groups) -> Vector{Vector{Step}}

`groups`, kinetic groups of `Step`s over `rxn`'s substrates and products, in canonical
form (`_canonicalize`) and read back as `Step`s (`steps`), with no check but the one
orienting needs: errors on a tied group no orientation aligns. The seed variants hold
their bases in it (`_seed_variants`).
"""
function _canonical_order(rxn::EnzymeReaction, groups::Vector{Vector{Step}})
    c = _canonicalize(_raw_of(rxn, groups))
    _check_aligned(c.mechanism, c.misaligned)
    steps(c.mechanism)
end

"""
The reaction of the mechanism `raw` writes: its substrates and products, with
placeholder atoms that balance for any counts (each substrate `length(products)`
carbons, each product `length(substrates)`), and its modifiers as competitive
inhibitors at multiplicity 1. Building it checks `raw`'s names as `EnzymeReaction`
checks a reaction's: none reserved, none holding `_`, and labels that spell every string
one way only. The mechanism macros build it before the mechanism
(`_compiled_mechanism`); a reaction's real atoms are written with `@enzyme_reaction`.
"""
function _placeholder_reaction(raw::RawMechanism)
    n_subs, n_prods = length(raw.substrates), length(raw.products)
    reactant_atoms = ReactantAtoms[
        [ReactantAtoms(Substrate(s), [:C => n_prods]) for s in raw.substrates];
        [ReactantAtoms(Product(p), [:C => n_subs]) for p in raw.products]]
    EnzymeReaction(reactant_atoms,
        RegulatorMults[RegulatorMults(CompetitiveInhibitor(x), [1]) for x in raw.modifiers],
        [1])
end

# ─── Pretty printing ────────────────────────────────────────────────────

function Base.show(io::IO, em::EnzymeMechanism)
    m = Mechanism(em)
    groups = steps(m)
    flat = collect(Iterators.flatten(groups))
    forms = _forms(groups)
    print(io, "EnzymeMechanism (", length(flat), " steps, ", length(forms),
          " enzyme forms):")
    for s in flat
        print(io, "\n  ", _step_text(s, m))
    end
    # The names the steps bind as modifiers, each once.
    mods = modifiers(m)
    isempty(mods) || print(io, " | modifiers: ", join(mods, ", "))
end

# Mirrors `@allosteric_mechanism` syntax: a kinetic group of several steps prints in
# parentheses, and each group carries its allosteric-state tag.
function Base.show(io::IO, aem::AllostericEnzymeMechanism)
    am = AllostericMechanism(aem)
    sites = regulatory_sites(am)
    print(io, "AllostericEnzymeMechanism (cat_n=", catalytic_multiplicity(am))
    isempty(sites) || print(io, ", ", length(sites), " reg sites")
    print(io, "):")
    for (g, group) in enumerate(steps(am))
        body = join(_step_text.(group, Ref(am)), ", ")
        print(io, "\n  ", length(group) == 1 ? body : "($body)", " :: ",
              cat_allo_state(am, g))
    end
    for (i, site) in enumerate(sites)
        ligs = ligands(site)
        print(io, "\n  reg site $i (n=", multiplicity(site), "): ", join(ligs, ", "),
              " [", join(("$l::$t" for (l, t) in zip(ligs, allo_states(site))), ", "), "]")
    end
end

# ─── Accessors ─────────────────────────────────────────────────
#
# A compiled mechanism is read by lifting it to its concrete `Mechanism` or
# `AllostericMechanism` (`Mechanism(em)`, `AllostericMechanism(aem)`). Only
# `metabolites`, `catalytic_mechanism` and `catalytic_multiplicity` are read on the
# compiled type itself.

"""Walk the steps of `m`, or of kinetic `groups`, in flat order, yielding
`(step::Step, kinetic_group::Int)` pairs."""
_flat_steps(groups::Vector{Vector{Step}}) =
    [(s, g) for (g, group) in enumerate(groups) for s in group]
_flat_steps(m::Mechanism) = _flat_steps(steps(m))

"""
    metabolites(m::EnzymeMechanism) → Tuple{Symbol,...}

Return distinct metabolite names as a tuple of `Symbol`s: substrates, then products,
then the names the steps bind as modifiers, each sorted by name, then the ligands of an
allosteric mechanism's regulatory sites, site by site, without duplicates. The fitter
uses this as the key set for the per-datapoint concentration `NamedTuple` passed to
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
@generated function metabolites(::M) where {M <: AbstractEnzymeMechanism}
    # Kept `@generated`: `loss!` uses `metabolites(m)` as a compile-time-constant
    # tuple to build the per-datapoint `NamedTuple{MetNames}` concs on the fitting
    # hot path. A runtime body would make that NamedTuple type-unstable and allocate.
    Tuple(metabolites(_concrete(M())))
end
