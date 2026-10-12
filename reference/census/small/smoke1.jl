# ABOUTME: Smoke test of the census pipeline (package side only): seeds, base tier, one
# ABOUTME: expansion round per reaction, funnel counts and timings; checks expand_labeled.
include(joinpath(@__DIR__, "lib.jl"))
t_load = time()
for name in (:uni_uni_rxn, :uni_uni_with_reg, :uni_uni_allo, :uni_uni_allo_reg,
             :uni_uni_allo_2reg, :uni_uni_reg_and_inhibitor, :uni_bi_rxn, :bi_bi_rxn,
             :bi_bi_pp_rxn)
    rxn = getfield(Main, name)
    run = Run(rxn, string(name))
    t = @elapsed f0 = round0!(run)
    @printf("%-28s r0: %s seeds=%d deg=%d flipadd=%d fresh=%d ok=%d fitted=%d cx=%d np=%d err=%d (%.1f s)\n",
            name, f0[:seed_source], f0[:seeds], f0[:degenerate_seeds], f0[:flip_children_added],
            f0[:fresh], f0[:ok], f0[:fitted], f0[:cx_skip], f0[:param_skip], f0[:error], t)
    flush(stdout)
    # expand_labeled equals expand_mechanisms on up to 20 parents
    par = frontier(run, 0)
    bad = 0
    for g in par[1:min(end, 20)]
        a = [c for (c, _) in expand_labeled(run.mechs[g], rxn)]
        b = ER.expand_mechanisms(Union{ER.Mechanism, ER.AllostericMechanism}[run.mechs[g]], rxn)
        a == b || (bad += 1)
    end
    println("   expand_labeled != expand_mechanisms on $bad of $(min(length(par), 20)) parents")
    name in (:bi_bi_rxn, :bi_bi_pp_rxn) && continue
    t = @elapsed f1 = expand_round!(run, par, 1)
    @printf("   r1: parents=%d raw=%d generated=%d seen=%d fresh=%d ok=%d fitted=%d eqdup=%d cx=%d np=%d err=%d expfail=%d (%.1f s)\n",
            f1[:parents], f1[:raw_children], f1[:generated], f1[:seen_skip], f1[:fresh], f1[:ok],
            f1[:fitted], f1[:eqhash_dup], f1[:cx_skip], f1[:param_skip], f1[:error],
            f1[:expand_failures], t)
    errs = unique(get(r, :err, "") for r in run.recs if r[:status] === :error)
    foreach(e -> println("   error: ", first(e, 200)), errs[1:min(end, 3)])
    flush(stdout)
end
@printf("total %.0f s; maxrss %.0f MB\n", time() - t_load, Sys.maxrss() / 2^20)
