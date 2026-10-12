# ABOUTME: Smoke test of the verification library on a handful of census candidates: structure,
# ABOUTME: admission, rank, same-law, growth, face oracle, :OnlyA family and law-vs-oracle.
include(joinpath(@__DIR__, "vlib.jl"))
include(joinpath(@__DIR__, "moracle.jl"))

function smoke(name, gids)
    d, mechs = loadrun(name)
    recs = d.recs
    rng = Xoshiro(7)
    for g in gids
        r = recs[g]; m = mechs[g]
        println("\n#$g ($name) ", r[:text])
        ms = mstruct(m)
        println("  forms $(length(ms.forms)) steps $(length(ms.steps)) groups $(ms.ngroups) frees $(ms.frees) allo $(ms.allo)")
        println("  my_admit $(my_admit(ms))  census admit: ", r[:admit] === nothing ? "-" : first(r[:admit], 90))
        t = @elapsed pl = plaw(m)
        t2 = @elapsed rk = cs_rank(pl, rng)
        @printf("  np %d (census %d) rank %d (census %d) kept %.1e dropped %.1e  [compile+plaw %.2f s, rank %.2f s]\n",
                length(pl.names), r[:np], rk.rank, r[:rank_pkg], rk.kept, rk.dropped, t, t2)
        gr = growth(pl, ms, rng)
        @printf("  growth slopes fwd %.3g rev %.3g\n", gr[:fwd], gr[:rev])
        t = @elapsed ft = face_test(ms, rng)
        @printf("  face test (%.2f s): positive pattern %.1e\n", t, ft.positive_pattern)
        for f in ft.faces
            @printf("    %-16s exact %s interior %s conv %.1e pattern %.1e ipattern %.1e gap %.1e %s\n", f.face,
                    string(round.(f.exact; sigdigits = 6)), string(round.(f.interior; sigdigits = 6)),
                    f.conv, f.pattern, f.ipattern, f.gap, f.fails ? "FAILS" : "")
        end
        if ms.allo
            t = @elapsed ot = onlya_test(ms, rng)
            println("  onlya test ($(round(t; digits = 2)) s): family $(ot.family) nfree $(ot.nfree)")
            for p in ot.points
                @printf("    %-16s pkg %.6g fam %s spread %.1e vs_pkg %.1e\n", p.point, p.pkg,
                        string(round.(p.fam; sigdigits = 6)), p.spread, p.vs_pkg)
            end
        else
            t = @elapsed lo = law_vs_oracle(m, ms, pl, rng)
            println("  law vs oracle ($(round(t; digits = 2)) s):")
            for p in lo
                @printf("    %-16s law %.6g interior %.6g exact %s law/int %.1e spread %.1e law/exact %.1e\n",
                        p.point, p.law, p.interior, string(round.(p.exact; sigdigits = 6)),
                        p.law_vs_interior, p.exact_spread, p.law_vs_exact)
            end
        end
    end
end

smoke("uni_bi_rxn", [1, 9, 13, 50, 213])
smoke("bi_bi_pp_rxn", [3048])
smoke("bi_bi_rxn.r2_100000", [14463, 17966, 17968])
