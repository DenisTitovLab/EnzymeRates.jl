# ABOUTME: First look at the saved census records: record fields, allosteric multiplicities and
# ABOUTME: sites, class labels per round, and the slice API on one plain and one allosteric case.
const CEN = "/tmp/claude-501/-home-denis-linux--julia-dev-EnzymeRates/12f42995-f5ba-4fda-9510-dd9b87f7d069/scratchpad/census/large"
include(joinpath(CEN, "scripts", "census_lib.jl"))
include(joinpath(CEN, "scripts", "summarize.jl"))

t0 = time()
for name in (:ter_bi_rxn, :ter_ter_rxn, :pyruvate_carboxylase_rxn)
    data = deserialize(joinpath(CEN, "data", "$(name).jls"))
    rounds = data.rounds
    label!(rounds)
    println("== $name: loaded at $(round(time() - t0; digits = 1)) s")
    r1 = first(r for recs in rounds for r in recs if r[:fitted] && get(r, :allo, false))
    println("record keys: ", sort!(collect(keys(r1))))
    for (ri, recs) in enumerate(rounds)
        fit = [r for r in recs if r[:fitted]]
        mults = Dict{Int, Int}(); nsites = Dict{Int, Int}()
        for r in fit
            sm = get(r, :smech, nothing)
            sm === nothing && continue
            mults[sm.multiplicity] = get(mults, sm.multiplicity, 0) + 1
            nsites[length(sm.sites)] = get(nsites, length(sm.sites), 0) + 1
        end
        tags = Dict{Vector{Symbol}, Int}()
        for r in fit
            sm = get(r, :smech, nothing)
            (sm === nothing || !S.is_allosteric(sm)) && continue
            t = sort!(unique(sm.tags))
            tags[t] = get(tags, t, 0) + 1
        end
        println("round $(ri - 1): fitted $(length(fit)); multiplicity $(mults); sites $(nsites); tag sets $(tags)")
        println("   no smech: ", count(r -> get(r, :smech, nothing) === nothing, fit),
                "; conv_err ", count(r -> haskey(r, :conv_err), fit),
                "; law_err ", count(r -> haskey(r, :law_err), fit))
    end
end

# The slice API on one allosteric record
data = deserialize(joinpath(CEN, "data", "ter_bi_rxn.jls"))
rounds = data.rounds
byid = Dict(r[:id] => r for recs in rounds for r in recs)
for id in ((1, 8), (1, 1157), (0, 5))
    r = byid[id]
    m = r[:smech]
    println("-- id $id: ", mech_text(m))
    println("   forms: ", [S.form_text(m, i) for i in eachindex(m.forms)])
    println("   mets: ", m.mets, " subs ", m.substrates, " prods ", m.products)
    println("   tags: ", m.tags, " group: ", m.group)
    println("   fitted: ", S.fitted_params(m))
    println("   constant_names: ", S.constant_names(m))
    S.is_allosteric(m) && println("   inactive names: ", S.constant_names(m; inactive = true))
    S.is_allosteric(m) && println("   inactive_mask: ", S.inactive_mask(m))
    cols, A, b = S.constraint_system(m)
    println("   constraint system: ", size(A), " cols ", length(cols))
    sol = S.solve(m)
    L = S.rate_law(m, sol)
    println("   law vars: ", L.vars.names)
    println("   admit: ", S.admit(m))
    println("   canonical order: ", S.canonical(L).order)
end
