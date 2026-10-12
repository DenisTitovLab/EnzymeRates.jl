# ABOUTME: Prints which package-output fields each converted set's Pop records carry, to pick
# ABOUTME: the records each solve-stage check can use.
include(joinpath(@__DIR__, "lib.jl"))
recs = load_records()
for set in unique(r.set for r in recs)
    rs = [r for r in recs if r.set == set]
    f(g) = count(g, rs)
    @printf("%-18s n=%4d allo=%4d ok=%4d fitted=%4d deps=%4d stepnames=%4d samples=%4d metab=%4d\n",
            set, length(rs), f(r -> S.is_allosteric(r.c.mech)), f(r -> r.p.ok),
            f(r -> !isempty(r.p.fitted_params)), f(r -> !isempty(r.p.dependent_monomials)),
            f(r -> !isempty(r.p.step_names) || !isempty(r.p.step_names_A)),
            f(r -> !isempty(r.p.samples)), f(r -> !isempty(r.p.metabolites)))
end
