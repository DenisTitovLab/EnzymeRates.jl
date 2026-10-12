# ABOUTME: Runs the package's RE->SS flip, kinetic-group split and dead-end moves on a parent set
# ABOUTME: and saves parents, children (structure only, as Pop.PMech), timings and counts.
# Usage: julia --project=<er_oracle> pkg_moves.jl <set>   (init_bibi_pp, expand_bibi,
#        allo_bibi, p1_enum_inh); init_bibi_pp also writes the set dead1 (level-1 dead ends).
using EnzymeRates, Random, Serialization, Printf
const S = "/tmp/claude-501/-home-denis-linux--julia-dev-EnzymeRates/12f42995-f5ba-4fda-9510-dd9b87f7d069/scratchpad"
include(joinpath(S, "proto/common/Population.jl"))
const OUT = joinpath(S, "slice/enum/data")
const SET = ARGS[1]

include(joinpath(@__DIR__, "pkg_lib.jl"))

_errtext(e) = first(sprint(showerror, e), 400)

"""`f()` and its wall time: the minimum of three runs when one run is under 50 ms."""
function timed(f)
    t = @elapsed r = f()
    t < 0.05 && for _ in 1:2
        t = min(t, @elapsed f())
    end
    r, t
end

count_of(m) = try ER._independent_param_count(m) catch; -1 end

"""Every move on every parent: the children (as Pop.PMech), the time, the error if any,
and the package's independent-parameter count of the parent and of each split child."""
function run_moves(parents, dead_rxns, label)
    out = NamedTuple[]
    kids = Dict(:flip => Any[], :split => Any[], :dead => Any[])
    for (i, m) in enumerate(parents)
        rec = Dict{Symbol, Any}(:index => i, :parent => to_pmech(m, "$label/$i"),
                                :count => count_of(m))
        for (mv, f) in ((:flip, () -> ER._expand_re_to_ss(m)),
                        (:split, () -> ER._expand_split_kinetic_group(m)),
                        (:dead, () -> ER._expand_add_dead_end_regulator(m, dead_rxns[i])))
            try
                cs, t = timed(f)
                rec[mv] = (children = [to_pmech(c, "$label/$i/$mv/$j") for (j, c) in enumerate(cs)],
                           time = t, error = "",
                           counts = mv === :split ? count_of.(cs) : Int[])
                append!(kids[mv], cs)
            catch e
                rec[mv] = (children = Pop.PMech[], time = NaN, error = _errtext(e), counts = Int[])
            end
        end
        push!(out, (; (k => v for (k, v) in rec)...))
        i % 50 == 0 && (@printf("%s: %d/%d parents, %.0f s\n", label, i, length(parents),
                                time() - T0); flush(stdout))
    end
    out, kids
end

const T0 = time()
pick(xs, n, seed) = length(xs) <= n ? collect(1:length(xs)) :
    sort!(randperm(Random.Xoshiro(seed), length(xs))[1:n])

function main()
    if SET == "init_bibi_pp" || SET == "expand_bibi"
        ps = Pop.load_population(SET)
        parents = [to_er(p, bi_bi_pp_rxn) for p in ps]
        recs, kids = run_moves(parents, fill(bi_bi_pp_inh, length(parents)), SET)
        serialize(joinpath(OUT, "pkg_$SET.jls"), (; set = SET, indices = collect(1:length(ps)),
                                                  records = recs))
        if SET == "init_bibi_pp"
            dead = unique(kids[:dead])
            sel = pick(dead, 300, 20261010)
            @printf("dead1: %d distinct level-1 dead-end children, %d picked\n",
                    length(dead), length(sel))
            recs2, _ = run_moves(dead[sel], fill(bi_bi_pp_inh, length(sel)), "dead1")
            serialize(joinpath(OUT, "pkg_dead1.jls"), (; set = "dead1", indices = sel,
                                                       records = recs2))
        end
    elseif SET == "allo_bibi"
        ps = Pop.load_population(SET)
        sel = pick(ps, 300, 20261011)
        recs, _ = run_moves([to_er(ps[i], bi_bi_allo_rxn) for i in sel],
                            fill(bi_bi_allo_inh, length(sel)), SET)
        serialize(joinpath(OUT, "pkg_$SET.jls"), (; set = SET, indices = sel, records = recs))
    elseif SET == "p1_enum_inh"
        ps = deserialize(joinpath(S, "proto/skeptic_P1-engine/data/enum_inh.jls"))
        sel = pick(ps, 300, 20261012)
        rxn(p) = startswith(p.id, "uni_uni") ? uni_uni_with_reg : bi_bi_with_inh
        recs, _ = run_moves([to_er(ps[i], rxn(ps[i])) for i in sel], [rxn(ps[i]) for i in sel],
                            SET)
        serialize(joinpath(OUT, "pkg_$SET.jls"), (; set = SET, indices = sel, records = recs))
    else
        error("unknown set $SET")
    end
    @printf("%s done in %.0f s; maxrss %.0f MB\n", SET, time() - T0, Sys.maxrss() / 2^20)
end
main()
