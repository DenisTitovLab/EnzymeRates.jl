# ABOUTME: The two-conformation mass-action network at one subunit (formulation 1) for the P5b oracle.
# ABOUTME: Included by validate.jl and two_free.jl after Allo.jl and Oracle.jl.

# ─── Oracle: the two-conformation network at n = 1 ───────────────────────

"""Formulation 1 at one subunit as a plain mechanism: the active graph, the inactive graph
(its kept steps on copies of their forms) and the flip `X ⇌ X_I` (RE, K = `L`) of the free
form `flip`, with the per-step constants of `vals`."""
function combined(md, vals; flip = findfirst(Pop.is_free, md.m.species), L = vals[:L])
    m = md.m
    n0 = length(m.species)
    iE = flip
    jsI = findall(md.keepI)
    formsI = sort!(unique([iE; [f for j in jsI for f in (m.steps[j].from, m.steps[j].to)]]))
    loc = Dict(f => n0 + k for (k, f) in enumerate(formsI))
    tagI(sp) = Pop.PSpecies(sp.bound, Symbol(:I_, sp.conformation), sp.res_added,
                            sp.res_subtracted, Symbol(:I_, sp.name))
    steps = [m.steps;
             [Pop.PStep(loc[m.steps[j].from], loc[m.steps[j].to], m.steps[j].consumed,
                        m.steps[j].released, m.steps[j].re) for j in jsI];
             Pop.PStep(iE, loc[iE], Pop.PMet[], Pop.PMet[], true)]
    k = [[[vals[Allo.colname(m, (:A, j, r))] for r in 1:Allo.nslots(s)] for (j, s) in enumerate(m.steps)];
         [[vals[Allo.colname(m, (:I, j, r))] for r in 1:Allo.nslots(m.steps[j])] for j in jsI];
         [[L]]]
    pm = Pop.PMech(; id = m.id, source = m.source, substrates = m.substrates, products = m.products,
                   inhibitors = m.inhibitors, metabolites = m.metabolites,
                   species = [m.species; tagI.(m.species[formsI])], steps,
                   groups = [[i] for i in eachindex(steps)])
    pm, k
end

