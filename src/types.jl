# ABOUTME: Core concrete types (EnzymeReaction, Mechanism, Step, Species, raw input,
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
Species: an enzyme form. `bound` is sorted by `_met_sort_key`. Its internal name
(`_species_name`: `:E` / `:EATP` / `:Estar...` / `:EB_res_+A_-P`) is rendered once at
construction from the form alone and stored in `name`; it orders and identifies forms
inside the package, and never takes part in `==`, `hash` or display. Parameter
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
more allosteric ligands, each given by its name. `ligands[i]` and `allo_states[i]` are
parallel. The constructor sorts the ligands by name and carries each ligand's state
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

"""
Canonicalize a non-binding step's storage direction to physical-forward, so
`from` is further from product-release / closer to substrate-binding.
Applies to RE AND SS steps — the direction question is identical;
only the parameter count differs. (All binding steps — RE and SS —
are stored with their metabolite consumed by the Step constructor; this
function orients every other step: isomerizations and transformations.)
Tier 1 reads atom progression, Tier 2 where substrates enter and products leave, and
Tier 3, for a step neither decides, runs it from the form earlier in the canonical form
order (`_form_index_key`: free E first). Reversing a step swaps its forms and its lists,
and every tier reads the two sides symmetrically, so the result does not depend on how
the step was written. Returns the step in that direction and the tier that decided it
(1, 2 or 3), or `(s, 0)` for a binding.
"""
function _canonical_step_direction(s::Step, entry_sides::Dict{Species, Tuple{Bool, Bool}})
    is_binding(s) && return s, 0
    f, t = from_species(s), to_species(s)

    # Tier 1: atom-balance progression over bound plus free metabolites, read by role:
    # a competitive-inhibitor copy counts on neither side, whatever its name.
    score(sp, free) =
        (count(m -> m isa Substrate, bound(sp)) + count(m -> m isa Substrate, free),
         -count(m -> m isa Product, bound(sp)) - count(m -> m isa Product, free))
    sf, st = score(f, consumed(s)), score(t, released(s))
    sf > st && return s, 1
    sf < st && return _reversed(s), 1

    # Tier 2: 1-hop (RE+SS) graph context — whether substrates and products
    # enter or leave solution at each form, as `(any substrate, any product)`.
    fk = get(entry_sides, f, (false, false))
    tk = get(entry_sides, t, (false, false))
    fk == (false, true) && tk == (true, false) && return s, 2
    fk == (true, false) && tk == (false, true) && return _reversed(s), 2

    # Tier 3: the canonical form order (C9), free E first.
    (_form_index_key(f) <= _form_index_key(t) ? s : _reversed(s)), 3
end

"""`s` stored in the opposite direction: forms swapped, consumed and released swapped."""
_reversed(s::Step) =
    Step(to_species(s), from_species(s), released(s), consumed(s), is_equilibrium(s))

"""
The steps of one kinetic group, `oriented` one by one by `_canonical_step_direction`
with the deciding `tiers`, except that the steps Tier 1 leaves tied turn as one. The
lead is the tied step ranked first by whether it changes conformation (those that do
first), then by deciding tier, then by `_step_index_key` (the canonical form order of
its two forms). Every other tied step takes up and gives off the lead's metabolites
(`_step_kind`), and one that changes conformation runs between the lead's two
conformations in the same order; a step that fits in both orientations keeps its own. A
group's steps share their constants, so a conformational change and its mirror at other
bound forms must be stored in the same physical direction; oriented one by one, Tier 2
can run the change one way and Tier 3 its mirror the other, tying a forward constant to
a reverse one. The same holds for one metabolite exchange at two conformations, which
keeps its conformation and turns by its metabolites alone. A tied step whose metabolites
match the lead's in neither orientation is left for `_assert_uniform_groups` to reject.
Errors when a tied step's metabolites match the lead's but its conformations do not,
since no one orientation then matches the grouping.
"""
function _orient_tied_steps(reaction::EnzymeReaction, oriented::Vector{Step},
                            tiers::Vector{Int})
    tied = findall(>(1), tiers)
    length(tied) < 2 && return oriented
    pair(s) = (conformation(from_species(s)), conformation(to_species(s)))
    turns(s) = pair(s)[1] != pair(s)[2]
    lead = oriented[argmin(i -> (!turns(oriented[i]), tiers[i],
                                 _step_index_key(oriented[i])), tied)]
    map(eachindex(oriented)) do i
        s = oriented[i]
        tiers[i] < 2 && return s
        same_kind = filter(o -> _step_kind(o) == _step_kind(lead), [s, _reversed(s)])
        isempty(same_kind) && return s
        fit = findfirst(o -> !turns(o) || pair(o) == pair(lead), same_kind)
        fit === nothing || return same_kind[fit]
        error("Mechanism: a kinetic group holds " *
              "$(join(_forward_sides(lead, reaction), " → ")) and " *
              "$(join(_forward_sides(first(same_kind), reaction), " → ")), which make " *
              "different conformational changes; the steps of a kinetic group share " *
              "their constants, so they must make the same conformational change")
    end
end

"""
Canonicalize the storage direction (RE + SS) of every non-binding step to
physical-forward for every group, turning the steps of a group that Tier 1 leaves
tied as one (`_orient_tied_steps`). Tier 2 reads each step's free metabolites
at both of its ends, so it sees the same context however the steps were
written. Shared by the `Mechanism` and `AllostericMechanism` constructors so
the Canonical Step Form invariant cannot drift between them.

The table `entry_sides` (`_entry_sides`) classifies each species by the metabolites that
enter or leave solution at it, as `(any substrate, any product)`: the consumed
metabolites of every step leaving it (its `from_species`) and the released
metabolites of every step arriving at it (its `to_species`). A binding is stored
with its metabolite consumed, so it marks the form the metabolite binds to; a
fused release E(S) → F + P, stored as the binding F + P → E(S), marks F, the form
P leaves at. Reversing a step swaps its forms and its lists together, so the
classification does not depend on how any step was written. Isomerizations carry
no free metabolites and mark nothing, and a competitive-inhibitor copy, neither a
substrate nor a product, marks nothing even when it shares a reactant's name.
`_canonical_step_direction` Tier 2 reads it to decide direction for non-binding steps
where Tier 1 ties.

Why ALL steps (not just RE): the "substrate-entry / product-exit"
property is a chemistry fact about which forms metabolites enter and
leave at — it does NOT depend on whether the step is rapid-
equilibrium or steady-state. The DSL parses `<-->` as SS and `⇌` as RE;
fixtures like Segel Iso Uni Uni (`E + A <--> EA ⇌ EP <--> F + P, F <--> E`)
use `<-->` throughout, so an RE-only filter would mis-classify both `E`
and `F` as marking nothing and the F⇌E case would fall through to Tier 3.
"""
function _canonicalize_step_directions(reaction::EnzymeReaction,
                                       groups::Vector{Vector{Step}})
    entry_sides = _entry_sides(groups)
    map(groups) do group
        oriented = [_canonical_step_direction(s, entry_sides) for s in group]
        _orient_tied_steps(reaction, first.(oriented), last.(oriented))
    end
end

"""
Each enzyme form of `groups` with the metabolites that enter or leave solution at it, as
`(any substrate, any product)`: the consumed metabolites of every step leaving it and the
released metabolites of every step arriving at it (`_canonicalize_step_directions`).
"""
function _entry_sides(groups::Vector{Vector{Step}})
    entry_sides = Dict{Species, Tuple{Bool, Bool}}()
    for group in groups, s in group, (form, free) in ((from_species(s), consumed(s)),
                                                      (to_species(s), released(s)))
        has_sub, has_prod = get(entry_sides, form, (false, false))
        entry_sides[form] = (has_sub  || any(m -> m isa Substrate, free),
                             has_prod || any(m -> m isa Product, free))
    end
    entry_sides
end

# Canonical key for a `Step`: its content in the canonical form order
# (`_step_index_key`: each form's `_form_index_key`, free E first, then the metabolites
# it takes up and gives off), then its RE/SS flag. Gives the sorts a deterministic
# ordering so two physically-equivalent `Mechanism`s end up with identical step storage.
# Keyed on content (base names and role flags), not on `hash`: `Substrate`/`Product` use
# Julia's default struct hash, which mixes in the type's `objectid` and so is NOT stable
# across precompile sessions — a hash key would make the canonical order
# non-deterministic, and with it the dedup keys and the rate equation.
_step_canonical_key(s::Step) = (_step_index_key(s)..., is_equilibrium(s))

"""
Canonical key for a `RegulatorySite`. Orders the outer site vector so two
`AllostericMechanism`s differing only in site presentation order collapse.
Keyed on ligand names (stable) rather than `hash` for the same reason.
"""
_regulatory_site_canonical_key(site::RegulatorySite) =
    (Tuple(String(l) for l in ligands(site)),
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
inner order. Each key is computed once per sort. Then reject two distinct forms that
render one name (`_assert_distinct_form_names`) and enforce the kinetic-group rules
(`_assert_uniform_groups`, `_assert_each_reaction_once`). Whether the groups give a
defined rate law is the boundary check's question (`_check_well_posed`), not the
constructors'. The returned groups are fresh vectors; the input is not mutated.
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
    _assert_distinct_form_names(reaction, gs)
    _assert_uniform_groups(reaction, gs)
    _assert_each_reaction_once(reaction, gs)
    gs, perm
end

"""Whether `n` names a substrate or product of `rxn`."""
_is_reactant(rxn::EnzymeReaction, n::Symbol) =
    any(ra -> name(metabolite(ra)) === n, reactants(rxn))

"""Whether `x` is the competitive-inhibitor copy of a substrate or product of `ctx` (a
mechanism, or the reaction of one under construction; `_is_reactant`): a reactant bound
as a modifier."""
_is_reactant_copy(x::Metabolite, ctx) =
    x isa CompetitiveInhibitor && _is_reactant(ctx, name(x))

"""
A metabolite as written in parameter names: its name, with suffix `mod` on a
competitive-inhibitor copy of a substrate or product of `ctx` (a mechanism, or the
reaction of one under construction; `_is_reactant_copy`), so the copy's forms and
constants stay distinct from those of the reactant itself. A competitive inhibitor that
is no reactant binds only as a modifier, so its name stands bare.
"""
_met_label(x::Metabolite, ctx) =
    _is_reactant_copy(x, ctx) ? String(name(x)) * "mod" : String(name(x))

"""
An enzyme form of `ctx` (a mechanism, or the reaction of one under construction) as
written in parameter names: its conformation label, then the labels (`_met_label`) of its
bound metabolites without a separator, then, for a residual, `_res`, `_plus<name>` per
substrate it adds and `_minus<name>` per product it subtracts, each in name order: `E`,
`EAB`, `EPmodQ`, `EB_res_plusA_minusP`.
"""
function _form_label(sp::Species, ctx)
    head = String(conformation(sp)) * join(_met_label(x, ctx) for x in bound(sp))
    r = residual(sp)
    isempty(r) && return head
    head * "_res" * join("_plus" * String(name(a)) for a in added(r)) *
        join("_minus" * String(name(p)) for p in subtracted(r))
end

"""
A step side of `ctx` (a mechanism, or the reaction of one under construction) as
written in parameter names: its enzyme form (`_form_label`) followed by its free
metabolites (`_met_label`), joined by `_`, e.g. `"EA_B"` for `E(A) + B`.
"""
_side_label(form::Species, mets::Vector{Metabolite}, ctx) =
    join([_form_label(form, ctx); [_met_label(x, ctx) for x in mets]], "_")

"""
The two sides of `s` in its stored direction, as parameter names write them for `ctx`
(a mechanism, or the reaction of one under construction): `from_species` with the
consumed metabolites, then `to_species` with the released ones. Every step constant is
named after them.
"""
_forward_sides(s::Step, ctx) =
    (_side_label(from_species(s), consumed(s), ctx),
     _side_label(to_species(s), released(s), ctx))

"""
The kind of a step for kinetic grouping: the metabolites it takes up and gives off,
`(consumed, released)`. Steps that share a kinetic group share their constants, so
they must take up and give off the same metabolites: two bindings of `m` share a
kind whether or not either runs chemistry, and every isomerization is one kind.
"""
_step_kind(s::Step) = (consumed(s), released(s))

"""
Error when two distinct enzyme forms of `steps`, a mechanism of `reaction`, render the
same name. A form's name joins its conformation and bound metabolites without a
separator, so E with NAD and P bound and E with NADP bound are both `ENADP`. Every
constant is named after the forms its step joins (`_form_label`), so the two forms
would share constant names. The forms' internal names (`name(::Species)`), which order
and identify forms inside the package, must differ as well.
"""
function _assert_distinct_form_names(reaction::EnzymeReaction, steps::Vector{Vector{Step}})
    # A form as `@enzyme_mechanism` writes it: `E(NAD, P)`, `E(A::Mod)`, `E(I)`,
    # `E(B; residual = A - P)`.
    function text(sp)
        mets = [_is_reactant_copy(m, reaction) ? "$(name(m))::Mod" : "$(name(m))"
                for m in bound(sp)]
        r = residual(sp)
        terms = join(name.(added(r)), " + ") * join(" - $(name(p))" for p in subtracted(r))
        res = isempty(r) ? "" : "; residual = " * lstrip(terms)
        "$(conformation(sp))(" * join(mets, ", ") * res * ")"
    end
    function check(key, seen, clash)
        for group in steps, s in group, sp in (from_species(s), to_species(s))
            other = get!(seen, key(sp), sp)
            other == sp || error("Mechanism: the enzyme forms $(text(other)) and " *
                                 "$(text(sp)) $(clash(key(sp))), so rename a metabolite " *
                                 "to tell the two forms apart")
        end
    end
    check(sp -> _form_label(sp, reaction), Dict{String, Species}(),
          k -> "both render the name $k; every constant is named after the forms its " *
               "step joins")
    check(name, Dict{Symbol, Species}(),
          k -> "share the internal name $k, which orders and identifies forms inside " *
               "the package")
end

"""
Error unless the steps of every kinetic group take up and give off the same
metabolites (`_step_kind`) and carry one RE/SS flag: a shared constant means the
same reaction type at the same speed.
"""
function _assert_uniform_groups(reaction::EnzymeReaction, steps::Vector{Vector{Step}})
    label(s) = join(_forward_sides(s, reaction), " → ") *
               (is_equilibrium(s) ? " (RE)" : " (SS)")
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
function _assert_each_reaction_once(reaction::EnzymeReaction,
                                    steps::Vector{Vector{Step}})
    seen = Dict{Tuple{String, String}, Tuple{Int, Step}}()
    flag(s) = is_equilibrium(s) ? "RE" : "SS"
    for (g, group) in enumerate(steps), s in group
        sides = minmax(_forward_sides(s, reaction)...)
        if haskey(seen, sides)
            g0, s0 = seen[sides]
            text = join(_forward_sides(s0, reaction), " ⇌ ")
            g0 == g && error("Mechanism: kinetic group $g holds the reaction " *
                             "$text twice; a mechanism holds each reaction once")
            flags = flag(s0) == flag(s) ? "" :
                " ($(flag(s0)) in group $g0, $(flag(s)) in group $g)"
            error("Mechanism: kinetic groups $g0 and $g both hold the reaction " *
                  "$text$flags; a reaction belongs to one kinetic group, whose " *
                  "constants are named after it")
        end
        seen[sides] = (g, s)
    end
end

"""The distinct enzyme forms of the steps of `groups`, in first-seen order, each step's
`from_species` before its `to_species`."""
_forms(groups) =
    unique(sp for g in groups for s in g for sp in (from_species(s), to_species(s)))

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

"""
    _bottomless_re_segment(m, steps) -> Union{Nothing, Vector{Species}}

The forms of a rapid-equilibrium segment of `steps`, kinetic groups over mechanism `m`'s
metabolites, whose weights all vanish at zero products, or all vanish at zero
substrates; `nothing` when no segment does. Each metabolite name takes its side from
`m`'s substrates and products, so a competitive-inhibitor copy, which shares its
reactant's name and concentration, counts on that reactant's side even in a graph that
holds no step taking up or giving off the reactant itself (the inactive conformation's
pruned graph, `_inactive_groups`).

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
function _bottomless_re_segment(m, steps::Vector{Vector{Step}})
    species, segments, extras = _re_segment_extras(steps)
    for segment in segments
        length(segment) < 2 && continue
        for side in (products(m), substrates(m))
            all(i -> any(in(side), keys(extras[i])), segment) && return species[segment]
        end
    end
    nothing
end

"""
    Mechanism

A non-allosteric enzyme mechanism: a `reaction::EnzymeReaction` plus
`steps::Vector{Vector{Step}}`, where the outer vector is kinetic groups and
each inner vector holds the steps that share that group's kinetic parameters.
The constructor canonicalizes step direction and sorts steps and groups,
so two mechanisms that differ only in how their steps were written collapse to
the same struct. `==` and `hash` read the steps, not the reaction, so the same steps
built on two reactions are one mechanism. Step order derives from structure alone, and
each rate constant is named after the member of its kinetic class that the joint solve
keeps (`_class_sides`).
Lift to the singleton derivation type with `EnzymeRates.compile_mechanism(m)` or
`EnzymeMechanism(m)`.
"""
struct Mechanism
    reaction::EnzymeReaction
    steps::Vector{Vector{Step}}
    function Mechanism(reaction::EnzymeReaction,
                       steps::Vector{Vector{Step}})
        new(reaction, first(_canonical_groups(reaction, steps)))
    end
end

steps(m::Mechanism) = m.steps
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
        ligs = [l for site in regulatory_sites for l in ligands(site)]
        # The reaction a mechanism macro builds lists no allosteric regulator, so the
        # ligands' names are checked here.
        _check_metabolite_names("AllostericMechanism", ligs)
        dup = findfirst(i -> ligs[i] in view(ligs, 1:i-1), eachindex(ligs))
        dup === nothing ||
            error("AllostericMechanism: ligand $(ligs[dup]) appears in two " *
                  "distinct regulatory sites; rendered Kreg names would " *
                  "collide. Same-ligand-two-sites is not enumerated.")
        new(reaction, cat_steps, cat_allo_states, catalytic_multiplicity, regulatory_sites)
    end
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

# `==` and `hash` by content for the value types whose fields are Vectors or structs: the
# default `==` would compare identity. Every field takes part, in declaration order, except
# a mechanism's `reaction`: a mechanism is its steps, tags, multiplicity and sites, whatever
# reaction it was built on. `Species` keeps its own pair: its stored `name` is derived.
# `Met` keeps Julia's `==`, which compares its two fields, and has its own `hash`: a
# struct's default hash mixes in its type's identity, which changes between precompile
# sessions.
for T in (Residual, RegulatorySite, Step, Krapid, Kfor, Krev, Kreg,
          ReactantAtoms, RegulatorMults, EnzymeReaction, Mechanism, AllostericMechanism,
          RawForm, RawStep, RawMechanism, RawSite, RawAllostericMechanism)
    fs = filter(!=(:reaction), fieldnames(T))
    @eval Base.:(==)(a::$T, b::$T) =
        $(reduce((l, r) -> :($l && $r), [:(a.$f == b.$f) for f in fs]))
    @eval Base.hash(x::$T, h::UInt) =
        $(foldl((acc, f) -> :(hash(x.$f, $acc)), fs;
                init = :(hash($(QuoteNode(nameof(T))), h))))
end

const _AnyMech = Union{Mechanism, AllostericMechanism}

reaction(m::_AnyMech) = m.reaction

"""The names of `m`'s substrates, sorted."""
substrates(m::_AnyMech) = Symbol[name(x) for x in substrates(reaction(m))]

"""The names of `m`'s products, sorted."""
products(m::_AnyMech) = Symbol[name(x) for x in products(reaction(m))]

"""Whether `n` names a substrate or product of `m`."""
_is_reactant(m::_AnyMech, n::Symbol) = _is_reactant(reaction(m), n)

"""Every name `m`'s steps bind as a modifier, a reactant's modifier copy included,
sorted."""
modifiers(m::_AnyMech) =
    sort!(unique!(Symbol[name(x) for group in steps(m) for s in group
                         for sp in (from_species(s), to_species(s))
                         for x in bound(sp) if x isa CompetitiveInhibitor]))

# The concentration names of a mechanism in the order its rate law reads them: the
# substrates, the products, then the other names its steps bind as modifiers, each by
# name; an allosteric mechanism's site ligands follow, site by site, each once.
metabolites(m::Mechanism) = unique!([substrates(m); products(m); modifiers(m)])
metabolites(am::AllostericMechanism) =
    unique!([substrates(am); products(am); modifiers(am);
             [l for site in regulatory_sites(am) for l in ligands(site)]])

steps(m::AllostericMechanism) = m.cat_steps
cat_allo_state(m::AllostericMechanism, g::Int) = m.cat_allo_states[g]
cat_allo_states(m::AllostericMechanism) = m.cat_allo_states
catalytic_multiplicity(m::AllostericMechanism) = m.catalytic_multiplicity
regulatory_sites(m::AllostericMechanism) = m.regulatory_sites
allosteric_regulators(m::AllostericMechanism) =
    unique(Symbol[l for s in regulatory_sites(m) for l in ligands(s)])

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

"""
    Mechanism(raw::RawMechanism)

The mechanism `raw` writes, built and checked by the `Mechanism` constructor: its
reaction holds `raw`'s names (`_placeholder_reaction`), and its kinetic groups hold
`raw`'s steps by group label (`_step_groups`).
"""
Mechanism(raw::RawMechanism) = Mechanism(_placeholder_reaction(raw), _step_groups(raw))

"""
    AllostericMechanism(raw::RawAllostericMechanism)

The allosteric mechanism `raw` writes, built and checked by the `AllostericMechanism`
constructor: its catalytic part read as `Mechanism(raw.catalytic)` reads it, the tag of
each group label, the catalytic multiplicity and the regulatory sites.
"""
function AllostericMechanism(raw::RawAllostericMechanism)
    cat = raw.catalytic
    reaction, groups = _placeholder_reaction(cat), _step_groups(cat)
    sites = RegulatorySite[RegulatorySite(site.ligands, site.multiplicity, site.allo_states)
                           for site in raw.sites]
    AllostericMechanism(reaction, groups, raw.cat_allo_states, raw.catalytic_multiplicity,
                        sites)
end

"""
The reaction of the mechanism `raw` writes: its substrates and products, with
placeholder atoms that balance for any counts (each substrate `length(products)`
carbons, each product `length(substrates)`), and its modifiers as competitive
inhibitors at multiplicity 1. A reaction's real atoms are written with
`@enzyme_reaction`.
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

"""The steps of `raw`, in written order and direction, as one vector of `Step`s per
group label: a modifier binds as its `CompetitiveInhibitor`, any other metabolite as the
substrate or product of its name."""
function _step_groups(raw::RawMechanism)
    metabolite_of(m::Met) = m.modifier ? CompetitiveInhibitor(m.name) :
                            m.name in raw.substrates ? Substrate(m.name) : Product(m.name)
    species_of(f::RawForm) = Species(Metabolite[metabolite_of(m) for m in f.bound],
        f.conformation, Residual(Substrate[Substrate(a) for a in f.added],
                                 Product[Product(p) for p in f.subtracted]))
    groups = [Step[] for _ in 1:maximum(raw.group; init = 0)]
    for (s, g) in zip(raw.steps, raw.group)
        push!(groups[g], Step(species_of(s.from), species_of(s.to),
                              Metabolite[metabolite_of(m) for m in s.consumed],
                              Metabolite[metabolite_of(m) for m in s.released],
                              s.rapid_equilibrium))
    end
    groups
end

# ─── Mechanism ↔ Sig (parametric ↔ non-parametric) conversion ──
#
# Every leaf in `sig` MUST be a valid Julia type-parameter value (isbits,
# Symbol, or Tuple of those). A type is NOT valid inside a tuple
# (`Val{(Int,)}` throws), and neither is `Pair{Symbol,Int}` — encode pairs as
# `Tuple{Symbol,Int}`. Vectors are NEVER valid — always wrap in `Tuple(...)`.
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

A regulator declared on the reaction that no step binds is left out of the Sig (e.g.
a dead-end inhibitor before any expansion move binds it): the reaction of
`Mechanism(em)` does not list it and it gets no parameter. Substrates and products
are always encoded. `Mechanism` (the working representation used during enumeration)
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
        RegulatorMults[
            RegulatorMults(met(r[1]), collect(Int, r[2])) for r in regulators_sig],
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
becomes a `RegulatorySite` of those ligand names. Mirrors `Mechanism(::EnzymeMechanism)`
— bridges the parametric ↔ non-parametric boundary so derivation code can walk
`AllostericMechanism` uniformly. It reads the type parameters at run time,
so it compiles once rather than once per mechanism type.
"""
function AllostericMechanism(@nospecialize(aem::AllostericEnzymeMechanism))
    CM, CS, RS = typeof(aem).parameters
    cm_mech = Mechanism(CM())
    multiplicity, cat_allo_states = CS
    sites = RegulatorySite[RegulatorySite(collect(Symbol, ligands_syms), mult,
                                          collect(Symbol, reg_allo_states))
                           for (ligands_syms, mult, reg_allo_states) in RS]
    AllostericMechanism(reaction(cm_mech), steps(cm_mech),
                        collect(Symbol, cat_allo_states),
                        multiplicity, sites)
end

"""Lift a singleton mechanism to its concrete `Mechanism` or `AllostericMechanism`."""
_concrete(@nospecialize(m::EnzymeMechanism)) = Mechanism(m)
_concrete(@nospecialize(m::AllostericEnzymeMechanism)) = AllostericMechanism(m)

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
        Tuple((Tuple(ligands(s)), multiplicity(s), Tuple(allo_states(s)))
              for s in regulatory_sites(am)),
    }()
end

# --- Pretty printing ---

"""A step of `ctx` (a mechanism, or the reaction of one under construction) as text: its
entry side, an arrow (`⇌` rapid equilibrium, `<-->` steady state) and its exit side. A
side lists the enzyme form as parameter names write it (`_form_label`), then the
metabolites the step consumes or releases."""
_step_text(s::Step, ctx) =
    join((_form_label(from_species(s), ctx), name.(consumed(s))...), " + ") *
    (is_equilibrium(s) ? " ⇌ " : " <--> ") *
    join((_form_label(to_species(s), ctx), name.(released(s))...), " + ")

function Base.show(io::IO, em::EnzymeMechanism)
    m = Mechanism(em)
    flat = collect(Iterators.flatten(steps(m)))
    forms = _forms(steps(m))
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
# `AllostericMechanism` (`Mechanism(em)`, `AllostericMechanism(aem)`) and walking its
# `reaction` and `steps` fields. Only `metabolites`, `catalytic_mechanism` and
# `catalytic_multiplicity` are read on the compiled type itself.

"""Walk the steps of `m`, or of kinetic `groups`, in flat order, yielding
`(step::Step, kinetic_group::Int)` pairs."""
_flat_steps(groups::Vector{Vector{Step}}) =
    [(s, g) for (g, group) in enumerate(groups) for s in group]
_flat_steps(m::Mechanism) = _flat_steps(steps(m))

"""
    metabolites(m::EnzymeMechanism) → Tuple{Symbol,...}

Return distinct metabolite names as a tuple of `Symbol`s: substrates, then products,
then regulators, each sorted by name, without duplicates. The fitter uses this as the
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
    Tuple(metabolites(Mechanism(EnzymeMechanism{Sig}())))
end

# ─── AllostericEnzymeMechanism Accessors ────────────────────────

catalytic_mechanism(::AllostericEnzymeMechanism{CM}) where {CM} = CM()
catalytic_multiplicity(::AllostericEnzymeMechanism{CM, CS}) where {CM, CS} = CS[1]

# `@generated` for the same reason as `metabolites(::EnzymeMechanism)`: `loss!` needs
# the tuple as a compile-time constant.
@generated function metabolites(::AllostericEnzymeMechanism{CM, CS, RS}) where {CM, CS, RS}
    ligs = (l for (site_ligands, _, _) in RS for l in site_ligands)
    Tuple(unique!(Symbol[metabolites(CM())..., ligs...]))
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

"""The `_forward_sides` of each kinetic class's kept member (`_kept_members`), keyed by
`(group, slot, inactive)`. `name(p, m)` computes them for each call, so a caller that
names many constants of one mechanism computes them once and passes them."""
_class_sides(m::_AnyMech) = Dict{Tuple{Int, Int, Bool}, Tuple{String, String}}(
    k => _forward_sides(s, m) for (k, s) in _kept_members(m))

"""
The sides that name constant `slot` of `step` in `state`: those of the kept member of
its kinetic class, from `sides` (`_class_sides(m)`). An `:EqualAI` step's constant is one
class in both conformations; the inactive copy of an `:OnlyA` step is a class of its own,
named after its own step. A slot the step's group has no class for (the reverse rate of
a rapid-equilibrium step) takes the sides of its slot-1 class, so every parameter of
every step renders.
"""
function _constant_sides(step::Step, slot::Int, state::Symbol, m::_AnyMech, sides)
    g = findfirst(group -> step in group, steps(m))
    g === nothing && error("Step not found in mechanism: $step")
    tag = m isa AllostericMechanism ? cat_allo_state(m, g) : :None
    state === :I && tag === :OnlyA && return _forward_sides(step, m)
    inactive = state === :I && tag === :NonequalAI
    get(() -> sides[(g, 1, inactive)], sides, (g, slot, inactive))
end

# Every step constant is named after the reaction of its class's kept member: a rate
# constant (`k_`) after the direction it drives, an equilibrium constant (`K_`) after
# the direction whose products-over-reactants it is. A binding's K is named in the
# release direction, so it is a dissociation constant (`K_ES_to_E_S`); an
# isomerization's or other RE step's K in the stored direction (`K_ES_to_EP`). `sides`
# holds the kept members' sides (`_class_sides`).
name(p::Kfor, m::_AnyMech, sides = _class_sides(m)) =
    _render_reaction("k_", _constant_sides(p.step, 1, p.state, m, sides), p.state)
name(p::Krev, m::_AnyMech, sides = _class_sides(m)) =
    _render_reaction("k_", reverse(_constant_sides(p.step, 2, p.state, m, sides)),
                     p.state)
function name(p::Krapid, m::_AnyMech, sides = _class_sides(m))
    s = _constant_sides(p.step, 1, p.state, m, sides)
    _render_reaction("K_", is_binding(p.step) ? reverse(s) : s, p.state)
end

"""
Regulator-site parameter: state tag + ligand name + "reg". No site index —
the AllostericMechanism constructor enforces that each ligand appears at most
once across all sites.
"""
name(p::Kreg, ::AllostericMechanism) =
    Symbol("K_", _state_tag(p.state), String(p.ligand), "reg")

"""
The Parameter(s) governing step `s` with the given allosteric state, by
`is_equilibrium(s)`. `_step_names` applies it to every step of one conformation, with
that conformation's tag, `_column_parameter` picks a joint-solve column's, and
`_onlya_haldane_violation` applies it to a group's first step to label the group.

Returns 1 element for an RE step (`Krapid`) and 2 elements for an SS step
(`Kfor`+`Krev`). `name` renders a binding's `Krapid`, plain or fused, as a
dissociation constant and every other step's — an isomerization, a Theorell–Chance
step, a step with two metabolites on one side — as the equilibrium constant of its
stored direction.
"""
_step_constants(s::Step, state::Symbol) =
    is_equilibrium(s) ? Parameter[Krapid(s, state)] :
    Parameter[Kfor(s, state), Krev(s, state)]
