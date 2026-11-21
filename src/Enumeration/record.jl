# ABOUTME: A candidate: the mechanism the enumerator proposes and its own record of each
# ABOUTME: form, with the two places that turn records into a mechanism and back.

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

# Equal mechanisms hold equal records: the enumerator builds each form from its contents.
Base.:(==)(a::CandidateMechanism, b::CandidateMechanism) = a.mechanism == b.mechanism
Base.hash(c::CandidateMechanism, h::UInt) = hash(c.mechanism, h)

"""
    _candidate(substrates, products, forms, steps, allostery = <plain>)
        -> CandidateMechanism
    _candidate(parent::CandidateMechanism, allosteric_states, oligomeric_state,
               regulatory_sites) -> CandidateMechanism

The candidate whose forms have the records `forms` and whose steps are `steps` (`from` and
`to` index `forms`), with the allostery `allostery`, whose `allosteric_states[l]` is the
state of written kinetic group `l`, built and checked by the `Mechanism` constructor; a
plain mechanism by default. The second gives `parent`'s catalytic scheme, with its
records, the allosteric states, oligomeric state and regulatory sites given. The only
place in `Enumeration` that builds a mechanism. The constructors put the forms in
canonical order, so the candidate holds the record of each form as the built mechanism
orders them.
"""
function _candidate(substrates::Vector{Symbol}, products::Vector{Symbol},
                    forms::Vector{_FormRecord}, steps::Vector{Mechanisms.Step},
                    allostery::Mechanisms.Allostery =
                        Mechanisms.Allostery(Symbol[], 1, Mechanisms.RegulatorySite[]))
    m = Mechanisms.Mechanism(
        Mechanisms.CatalyticScheme(substrates, products, _enzyme_form.(forms), steps),
        allostery)
    CandidateMechanism(m, _form_records(m))
end
_candidate(parent::CandidateMechanism, allosteric_states::Vector{Symbol},
           oligomeric_state::Int, regulatory_sites::Vector{Mechanisms.RegulatorySite}) =
    CandidateMechanism(
        Mechanisms.Mechanism(parent.mechanism.catalytic,
                             Mechanisms.Allostery(allosteric_states, oligomeric_state,
                                                  regulatory_sites)),
        parent.forms)

"""The `Mechanisms` form that the record `f` describes."""
_enzyme_form(f::_FormRecord) =
    Mechanisms.EnzymeForm(f.label, [Mechanisms.Metabolite.(f.productive, false);
                                    Mechanisms.Metabolite.(f.copies, true)],
                          f.adds, f.subtracts)

"""
    _form_records(m) -> Vector{_FormRecord}

The record of each form of the catalytic scheme of mechanism `m`, read from the form's
contents, in form order. It makes the record of a mechanism written by hand (and, through
`_candidate`, of one the enumerator builds): the only place in `Enumeration` that reads a
form's contents from a mechanism.
"""
_form_records(m::Mechanisms.Mechanism) =
    _FormRecord[_FormRecord(n.conformation, Symbol[x.name for x in n.bound if !x.modifier],
                            Symbol[x.name for x in n.bound if x.modifier], n.added,
                            n.subtracted)
                for n in m.catalytic.nodes]
