# ABOUTME: Face oracle on large random samples of fitted candidates (plain and n = 1 allosteric),
# ABOUTME: tabulated against the census's structural face verdicts to estimate misses and false alarms.
include(joinpath(@__DIR__, "vlib.jl"))
include(joinpath(@__DIR__, "moracle.jl"))

rng = Xoshiro(919191)
census_face(r) = any(!=(:none), r[:face_kind]) || any(r[:dvan])
for (name, k) in (("uni_bi_rxn", 0), ("bi_bi_rxn", 300), ("bi_bi_pp_rxn", 300),
                  ("bi_bi_rxn.r2_100000", 400), ("bi_bi_pp_rxn.r2_1200", 400))
    d, mechs = loadrun(name)
    fitted = [r[:gid] for r in d.recs if r[:status] === :ok && r[:fitted] &&
              (!r[:allo] || get(r, :allo_n, 1) == 1)]
    pick = k == 0 ? fitted : fitted[sort!(randperm(rng, length(fitted))[1:min(k, length(fitted))])]
    tab = Dict{Tuple{Bool, Bool}, Int}()
    t = @elapsed for g in pick
        r = d.recs[g]
        ft = try
            face_test(mstruct(mechs[g]), rng)
        catch e
            println("  error #$g: ", errmsg(e)); continue
        end
        bad = any(f -> f.fails, ft.faces)
        key = (census_face(r), bad)
        tab[key] = get(tab, key, 0) + 1
        key[1] != key[2] && println("  disagreement #$g census $(key[1]) oracle $(key[2]) ",
            join(("$(f.face): pattern $(round(f.pattern; sigdigits = 2)) gap $(round(f.gap; sigdigits = 2))"
                  for f in ft.faces if f.fails), "; "), " | ", r[:text])
    end
    @printf("[%s] %d candidates (%.0f s): (census flag, oracle fails) %s\n", name, length(pick), t,
            join(("$(k[1])/$(k[2]) $v" for (k, v) in sort!(collect(tab))), ", "))
    d = mechs = nothing
    GC.gc()
end
