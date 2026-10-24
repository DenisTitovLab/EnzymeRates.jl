# ABOUTME: Verifier helpers on the package side: a mechanism as plain Julia data (the slice's raw
# ABOUTME: input, tags and sites included), labelled expansion, and one-line records for the slice.
const ER = EnzymeRates

"""`m` as plain data `(subs, prods, steps, group, tags, n, sites)`, each step
`(from, to, consumed, released, re)`, a form `(conf, bound, added, subtracted)`, a metabolite
`(name, binds as a competitive inhibitor)`."""
function rawrec(m)
    rxn = ER.reaction(m)
    subs = Symbol[ER.name(x) for x in ER.substrates(rxn)]
    prods = Symbol[ER.name(x) for x in ER.products(rxn)]
    met(x) = (ER.name(x), x isa ER.CompetitiveInhibitor)
    form(sp) = (Symbol(ER.conformation(sp)), Tuple{Symbol, Bool}[met(x) for x in ER.bound(sp)],
                Symbol[ER.name(x) for x in ER.added(ER.residual(sp))],
                Symbol[ER.name(x) for x in ER.subtracted(ER.residual(sp))])
    steps, group = Any[], Int[]
    for (g, grp) in enumerate(ER.steps(m)), s in grp
        push!(steps, (form(ER.from_species(s)), form(ER.to_species(s)),
                      Tuple{Symbol, Bool}[met(x) for x in ER.consumed(s)],
                      Tuple{Symbol, Bool}[met(x) for x in ER.released(s)],
                      ER.is_equilibrium(s)))
        push!(group, g)
    end
    sites = Tuple{Vector{Symbol}, Vector{Symbol}, Int}[]
    tags, n = Symbol[], 1
    if m isa ER.AllostericMechanism
        lig(l) = l isa Symbol ? l : ER.name(l)
        for s in ER.regulatory_sites(m)
            push!(sites, (Symbol[lig(l) for l in ER.ligands(s)],
                          Symbol[Symbol(t) for t in ER.allo_states(s)], ER.multiplicity(s)))
        end
        tags = Symbol[Symbol(t) for t in ER.cat_allo_states(m)]
        n = ER.catalytic_multiplicity(m)
    end
    (subs, prods, steps, group, tags, n, sites)
end

"""The children of `m` with the move that made each, in `expand_mechanisms`' order and with
its filter and check (graph_mechanism_wt src/mechanism_enumeration.jl:1920-1940)."""
function children_labeled(m, rxn)
    ER._assert_emission_rules(m)
    out = Tuple{Any, Symbol}[]
    put!(cs, mv) = foreach(c -> push!(out, (c, mv)), cs)
    put!(ER._expand_re_to_ss(m), :re_to_ss)
    put!(ER._expand_split_kinetic_group(m), :split)
    put!(ER._expand_add_dead_end_regulator(m, rxn), :dead_end)
    if m isa ER.Mechanism
        put!(ER._expand_to_allosteric(m, rxn), :to_allosteric)
    else
        put!(ER._expand_add_allosteric_regulator(m, rxn), :add_allo_reg)
        put!(ER._expand_change_allo_state(m), :change_allo_state)
        put!(ER._expand_merge_regulatory_sites(m), :merge_sites)
    end
    filter!(c -> ER._respects_reg_type(c[1], rxn), out)
    foreach(c -> ER._assert_atom_conserving(c[1]), out)
    out
end

cxA(m) = Int(ER._eq_complexity(m))
cxI(m) = m isa ER.AllostericMechanism ? Int(ER._eq_complexity(ER._inactive_groups(m))) : 0
