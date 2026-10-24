# ABOUTME: Re-checks the comparator-in-variant direction at fixed L for chosen variants with more
# ABOUTME: restarts and three L values, to tell a real gap at fixed L from an optimizer miss.
# usage: deep.jl <data dir> <rxn> <gid>...
push!(ARGS, "", "")          # contain.jl reads three ARGS at load; main() does not run
include(joinpath(@__DIR__, "contain.jl"))
rxn = ARGS[2]
for gs in ARGS[3:(end - 2)]
    gid = parse(Int, gs)
    r = recs(rxn)[gid]
    subs, prods = roles(rxn)
    sm = parse_neutral(r[10], subs, prods)
    pm = parse_neutral(r[10], subs, prods; plain = true)
    f = Target(blaw(S.rate_law(sm)), filter(!=(:L), collect(S.fitted_params(sm))))
    g = Target(blaw(S.law(pm)), collect(S.fitted_params(pm)))
    for Lv in (0.01, 0.1, 1.0, 10.0), trial in 1:3
        rng = Xoshiro(1000 * trial + gid)
        e = contained(g, f, [:L], rng; restarts = 30, fixval = Lv)
        @printf("%s %d L=%g trial %d: comparator in variant at fixed L, rms log error %.1e\n",
                rxn, gid, Lv, trial, e)
        flush(stdout)
    end
end
