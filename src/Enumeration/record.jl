# ABOUTME: A candidate: the mechanism the enumerator proposes and its own record of each
# ABOUTME: form, the naming rule for the forms it builds, and the one place that builds one.

"""
    _FormRecord

The enumerator's record of one enzyme form: its conformation label (`label`), the names it
binds productively (`productive`) and as modifier copies (`copies`, a pure inhibitor among
them), and its covalent residual, the substrates it adds (`adds`) and the products it
subtracts (`subtracts`). Each list is sorted and holds a name as often as the form does.
"""
struct _FormRecord
    label::Symbol
    productive::Vector{Symbol}
    copies::Vector{Symbol}
    adds::Vector{Symbol}
    subtracts::Vector{Symbol}
end

Base.:(==)(a::_FormRecord, b::_FormRecord) =
    a.label === b.label && a.productive == b.productive && a.copies == b.copies &&
    a.adds == b.adds && a.subtracts == b.subtracts
Base.hash(f::_FormRecord, h::UInt) =
    hash(f.subtracts, hash(f.adds, hash(f.copies, hash(f.productive,
         hash(f.label, hash(:_FormRecord, h))))))

"""
The name of the form whose record is `f`, in a mechanism whose reactants are `reactants`:
its label, then the names it binds in name order, a reactant's modifier copy written
`<name>mod` after its reactant, then any residual as `_res`, `_plus<name>` per substrate it
adds and `_minus<name>` per product it subtracts (`E`, `EAB`, `EPmodQ`,
`EB_res_plusA_minusP`). The reaction's names spell every name one way
(`_check_decodable_names`), so two forms take one name only when their records are equal.
"""
function _form_name(f::_FormRecord, reactants::Vector{Symbol})
    io = IOBuffer()
    print(io, f.label)
    for (x, copy) in sort!([[(x, false) for x in f.productive];
                            [(x, true) for x in f.copies]])
        print(io, x, copy && x in reactants ? "mod" : "")
    end
    if !isempty(f.adds) || !isempty(f.subtracts)
        print(io, "_res")
        foreach(x -> print(io, "_plus", x), f.adds)
        foreach(x -> print(io, "_minus", x), f.subtracts)
    end
    Symbol(String(take!(io)))
end

"""
    CandidateMechanism

A mechanism the enumerator proposes for a reaction (`candidate_mechanisms`). Its field
`mechanism`, a `Mechanism`, is the mechanism to fit and report. Two candidates are `==`
and hash alike exactly when their mechanisms are.
"""
struct CandidateMechanism
    mechanism::Mechanisms.Mechanism
    # The record of each form, at the form's index in `mechanism`'s catalytic scheme. The
    # moves read a form's contents here, never from the mechanism.
    forms::Vector{_FormRecord}
end

# Equal mechanisms hold equal records: the enumerator names each form after its record.
Base.:(==)(a::CandidateMechanism, b::CandidateMechanism) = a.mechanism == b.mechanism
Base.hash(c::CandidateMechanism, h::UInt) = hash(c.mechanism, h)

"""
    _candidate(metabolites, forms, substrates, products, steps, allostery = <plain>)
        -> CandidateMechanism
    _candidate(parent::CandidateMechanism, allosteric_states, oligomeric_state,
               regulatory_sites) -> CandidateMechanism

The candidate whose catalytic scheme has the metabolites `metabolites`, of which
`substrates` and `products` index the reactants, the forms whose records are `forms`, each
named after its record (`_form_name`), and the steps `steps` (`from` and `to` index
`forms`, the lists `metabolites`), with the allostery `allostery`, whose
`allosteric_states[l]` is the state of written kinetic group `l`, built and checked by the
`Mechanism` constructor; a plain mechanism by default. A record given twice names one
form; two different records that take one name are an internal error. The second gives
`parent`'s catalytic scheme, with its records, the allosteric states, oligomeric state and
regulatory sites given. The only place in `Enumeration` that builds a mechanism; the
candidate holds each form's record in the order the constructor
stores the forms.
"""
function _candidate(metabolites::Vector{Symbol}, forms::Vector{_FormRecord},
                    substrates::Vector{Int}, products::Vector{Int},
                    steps::Vector{Mechanisms.Step},
                    allostery::Mechanisms.Allostery =
                        Mechanisms.Allostery(Symbol[], 1, Mechanisms.RegulatorySite[]))
    reactants = metabolites[[substrates; products]]
    names = Symbol[_form_name(f, reactants) for f in forms]
    record = Dict{Symbol, _FormRecord}()
    for (n, f) in zip(names, forms)
        get!(record, n, f) == f ||
            error("_candidate: two different form records are both named $n")
    end
    distinct = unique(names)
    if length(distinct) < length(names)
        at = Dict{Symbol, Int}(n => i for (i, n) in enumerate(distinct))
        steps = Mechanisms.Step[Mechanisms.Step(at[names[e.from]], at[names[e.to]],
                                                e.taken_up, e.given_off,
                                                e.rapid_equilibrium, e.kinetic_group)
                                for e in steps]
    end
    m = Mechanisms.Mechanism(
        Mechanisms.CatalyticScheme(metabolites, distinct, substrates, products, steps),
        allostery)
    CandidateMechanism(m, _FormRecord[record[n] for n in m.catalytic.enzyme_forms])
end
_candidate(parent::CandidateMechanism, allosteric_states::Vector{Symbol},
           oligomeric_state::Int, regulatory_sites::Vector{Mechanisms.RegulatorySite}) =
    CandidateMechanism(
        Mechanisms.Mechanism(parent.mechanism.catalytic,
                             Mechanisms.Allostery(allosteric_states, oligomeric_state,
                                                  regulatory_sites)),
        parent.forms)
