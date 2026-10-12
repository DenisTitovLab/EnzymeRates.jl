# ABOUTME: Summarizes the P5b validation results (data/results_<pop>.jls) per population.
# ABOUTME: Stdlib only; prints counts, value errors, ranks, normalization cases and text checks.
using Serialization, Printf, Statistics
include(joinpath(@__DIR__, "..", "common", "Population.jl"))

pops = isempty(ARGS) ? ["specs", "handwritten", "allo_uniuni", "allo_bibi", "expand_bibi"] : ARGS
mx(v) = isempty(v) ? NaN : maximum(v)
tally(v) = join(("$k=$c" for (k, c) in sort!([(k, count(==(k), v)) for k in unique(v)]; by = string ∘ first)), ", ")
for name in pops
    path = joinpath(@__DIR__, "data", "results_$(name).jls")
    isfile(path) || continue
    rs = deserialize(path)
    ok = filter(r -> !haskey(r, :error), rs)
    println("== $name: $(length(rs)) allosteric, errors $(length(rs) - length(ok))")
    for r in rs
        haskey(r, :error) && println("   ERROR $(r[:id]): $(first(r[:error], 300))")
    end
    println("   inactive graph differs from the package's: ", count(r -> !r[:keepI_match], ok))
    println("   non-integer exponent in the solve: ", count(r -> r[:nonint], ok),
            "; dependent RE constants written through independents in the law: ",
            count(r -> r[:re_subst] > 0, ok), " mechanisms")
    d = [r[:count] - r[:pkg_count] for r in ok]
    println("   count mine - pkg: ", tally(d))
    println("   (a) fitted names outside the declared constants: mine ",
            count(r -> !isempty(r[:undeclared_mine]), ok), ", package ",
            count(r -> !isempty(r[:undeclared_pkg]), ok), " (",
            sum(r -> length(r[:undeclared_pkg]), ok; init = 0), " names)")
    println("   (b) package: opposite folds ", count(r -> r[:pkg_opp], ok), ", folded in one state ",
            "but used by the other ", count(r -> r[:pkg_one], ok), "; mine: a fold or dependent in ",
            "the law or fitted list ", count(r -> r[:dep_in_law] > 0, ok), ", opposite twin folds ",
            count(r -> r[:opp_mine] > 0, ok), "; my folds: ", sum(r -> r[:folds], ok; init = 0),
            " in ", count(r -> r[:folds] > 0, ok), " mechanisms (collapses ",
            sum(r -> r[:collapses], ok; init = 0), ")")
    # (c) ranks
    rk = filter(r -> haskey(r, :rank_mine), ok)
    if !isempty(rk)
        println("   (c) ranks computed on $(length(rk)): mine count == rank ",
                count(r -> r[:rank_mine][1] == r[:count], rk), ", package count == rank ",
                count(r -> r[:rank_pkg][1] == r[:pkg_count], rk))
        gaps = [(r[:rank_mine][2], r[:rank_mine][3]) for r in rk]
        @printf("       rank gap: smallest kept pivot %.1e, largest dropped %.1e (relative)\n",
                minimum(first.(gaps)), maximum(last.(gaps)))
        for r in rk
            (r[:rank_mine][1] != r[:count] || r[:rank_pkg][1] != r[:pkg_count]) &&
                println("       rank != count: $(r[:id]) mine $(r[:count]) rank $(r[:rank_mine][1]); ",
                        "pkg $(r[:pkg_count]) rank $(r[:rank_pkg][1])  $(first(r[:source], 60))")
        end
        dif = filter(r -> r[:count] != r[:pkg_count], rk)
        println("       count differs on $(length(dif)): mine == its rank ",
                count(r -> r[:rank_mine][1] == r[:count], dif), ", package == its rank ",
                count(r -> r[:rank_pkg][1] == r[:pkg_count], dif))
    end
    # (d) values
    same = filter(r -> r[:count] == r[:pkg_count] && r[:cons_err] < 1e-10, ok)
    @printf("   (d) same model (equal count, package point satisfies my constraints): %d of %d; ",
            length(same), length(ok))
    @printf("my v vs package sample v max %.1e, vs package law in BigFloat max %.1e\n",
            mx([r[:v_vs_sample] for r in same]), mx([r[:v_vs_pkgbig] for r in same]))
    eqc_bad = filter(r -> r[:count] == r[:pkg_count] && r[:cons_err] >= 1e-10, ok)
    println("       equal count but the package point violates my constraints: ", length(eqc_bad),
            isempty(eqc_bad) ? "" : " e.g. " * join([r[:id] for r in eqc_bad[1:min(5, end)]], ", "))
    @printf("       package law rebuilt from its exported pieces vs its sample v: max %.1e\n",
            mx([r[:pkg_rebuild] for r in ok]))
    @printf("   (d1) my point in the package's names: name spread max %.1e; package bare point vs mine max %.1e; v max %.1e\n",
            mx([r[:restrict_spread] for r in ok]), mx([r[:restrict_bare] for r in ok]),
            mx([r[:restrict_v] for r in ok]))
    bad = filter(r -> r[:restrict_v] > 1e-20 || r[:restrict_spread] > 1e-20, ok)
    println("       package law does not contain mine: ", length(bad),
            isempty(bad) ? "" : " e.g. " * join([r[:id] for r in bad[1:min(5, end)]], ", "))
    @printf("   detailed balance |v(Keq)|/|v(Keq/10)|: mine max %.1e, package max %.1e; n = 1 mass-action oracle vs mine max %.1e\n",
            mx([r[:db_mine] for r in ok]), mx([r[:db_pkg] for r in ok]), mx([r[:oracle_n1] for r in ok]))
    println("   normalization case package: ", tally([r[:norm_pkg] for r in ok]),
            " | mine: ", tally([r[:norm_mine] for r in ok]),
            " | differ: ", count(r -> r[:norm_pkg] != r[:norm_mine], ok))
    println("   free-enzyme weight kind (A, I) package: ", tally([r[:dkind_pkg] for r in ok]))
    println("   free-enzyme weight kind (A, I) mine: ", tally([r[:dkind_mine] for r in ok]))
    println("   (f) printed-law artifacts: mine ", count(r -> r[:artifacts_mine] > 0, ok),
            " mechanisms, package ", count(r -> r[:artifacts_pkg] > 0, ok), " (",
            sum(r -> r[:artifacts_pkg], ok; init = 0), " occurrences)")
    println("       state-symmetric polynomial pairs: mine ", sum(r -> r[:sym_pairs], ok; init = 0),
            ", printed identically ", sum(r -> r[:sym_same_text], ok; init = 0),
            "; package ", sum(r -> r[:pkg_sym_pairs], ok; init = 0), ", same term order ",
            sum(r -> r[:pkg_sym_same_order], ok; init = 0))
    ts = [r[:derive_s] for r in ok]
    @printf("   derive time median %.1f ms, max %.1f ms (first call includes compilation)\n",
            1e3 * median(ts), 1e3 * maximum(ts))
end
