# ABOUTME: Counts the bi-bi and ping-pong rounds (package pipeline only, no facts) and times
# ABOUTME: them, to size the census runs; prints per-move child counts of round 1.
include(joinpath(@__DIR__, "lib.jl"))
for name in (:bi_bi_rxn, :bi_bi_pp_rxn)
    rxn = getfield(Main, name)
    run = Run(rxn, string(name))
    t = @elapsed f0 = round0!(run)
    @printf("%s r0: fresh=%d ok=%d fitted=%d (%.1f s)\n", name, f0[:fresh], f0[:ok], f0[:fitted], t)
    t = @elapsed f1 = expand_round!(run, frontier(run, 0), 1)
    @printf("   r1: parents=%d raw=%d generated=%d seen=%d fresh=%d ok=%d fitted=%d eqdup=%d cx=%d np=%d err=%d expfail=%d (%.1f s; expand %.1f s)\n",
            f1[:parents], f1[:raw_children], f1[:generated], f1[:seen_skip], f1[:fresh], f1[:ok],
            f1[:fitted], f1[:eqhash_dup], f1[:cx_skip], f1[:param_skip], f1[:error],
            f1[:expand_failures], t, f1[:t_expand])
    mv = Dict{Symbol, Int}()
    for r in run.recs
        r[:round] == 1 || continue
        mv[r[:move]] = get(mv, r[:move], 0) + 1
    end
    println("   r1 fresh by move: ", sort!(collect(mv)))
    mvf = Dict{Symbol, Int}()
    for r in run.recs
        r[:round] == 1 && get(r, :fitted, false) || continue
        mvf[r[:move]] = get(mvf, r[:move], 0) + 1
    end
    println("   r1 fitted by move: ", sort!(collect(mvf)))
    errs = unique(first(get(r, :err, ""), 150) for r in run.recs if r[:status] === :error)
    foreach(e -> println("   error: ", e), errs[1:min(end, 3)])
    flush(stdout)
end
@printf("maxrss %.0f MB\n", Sys.maxrss() / 2^20)
