# ABOUTME: Today's bugs on the integrated pipeline, from the saved runs: :OnlyA verdicts on mixed
# ABOUTME: groups and second free forms, kcat errors, the slowest laws, same-model regroupings.
include(joinpath(@__DIR__, "lib.jl"))

load(prefix) = [x for l in readdir(IDATA) if startswith(l, prefix)
                for x in deserialize(joinpath(IDATA, l))]
pipe, bench, keys_ = load("pipe_"), load("bench_"), load("keys_")
byid = Dict((x.set, x.id) => x for x in pipe)

println("== admission (decision 6): rejected per set")
for set in sort!(unique(x.set for x in pipe if x.allo))
    v = filter(x -> x.set == set && x.allo, pipe)
    rej = count(x -> !x.ok, v)
    rej > 0 && @printf("  %-18s %4d of %4d allosteric rejected\n", set, rej, length(v))
end

println("== hand-built cases (verdict: message)")
for x in sort(filter(x -> startswith(x.set, "p5a_adv"), pipe); by = x -> x.id)
    occursin(r"^AD(1|2|3b|3d|5|8|11|14|15)", x.id) || continue
    @printf("  %-30s %s\n", x.id, x.ok ? "admitted; kcat: " * (x.kcat isa String ?
            x.kcat : "$(x.kcat) candidates") : "rejected: " * first(x.msg, 150))
end

println("== kcat on the mechanisms where the package says \"no kcat components\"")
for id in ("handwritten/3", "handwritten/6", "handwritten/10", "handwritten/20")
    x = byid[("handwritten", id)]
    println("  ", id, ": ", x.kcat isa String ? x.kcat : "$(x.kcat) candidates")
end
println("== kcat errors per set (derived laws)")
for set in sort!(unique(x.set for x in pipe))
    v = filter(x -> x.set == set && x.ok && x.kcat isa String, pipe)
    isempty(v) && continue
    kinds = Dict{String, Int}()
    foreach(x -> (k = replace(first(x.kcat, 64), r"with [A-Za-z0-9_, ]+" => "with …");
                  kinds[k] = get(kinds, k, 0) + 1), v)
    println("  ", set, ": ", join(("$n × $k" for (k, n) in kinds), "; "))
end

println("== speed of the largest enumerated laws and the slowest laws")
for id in ("expand_bibi/626", "expand_bibi/642", "expand_bibi/658")
    x = only(filter(x -> x.id == id, bench))
    @printf("  %-16s %6.1f ns, %d allocations, %d ops\n", id, x.t1n_ns, x.allocs, x.ops)
end
for x in first(sort(bench; by = x -> -x.t1n_ns), 6)
    @printf("  slowest: %-18s %-20s %6.1f ns, ops %d, D terms %d\n", x.set, x.id, x.t1n_ns,
            x.ops, x.nD)
end

println("== same-model regroupings (P2's 663 single-group splits of the 806 plain laws)")
p2sets = ("specs", "init_uniuni", "init_unibi", "init_bibi", "init_bibi_pp", "expand_bibi")
v = filter(x -> x.set in p2sets && !x.allo, keys_)
@printf("  splits %d, same model %d: key changes %d, fitted-name-set changes %d \
        (package: 1824 / 663 / 327 / 312)\n", sum(x.p2 for x in v), sum(x.p2same for x in v),
        sum(x.p2key for x in v), sum(x.p2names for x in v))
s(f, w) = sum(getfield(x, f) for x in w; init = 0)
@printf("  all records: same-model splits %d of %d; key changes %d; fitted-name changes %d; \
        solve changes %d\n", s(:nsame, keys_), s(:nsplit, keys_), s(:keyfail, keys_),
        s(:namechange, keys_), s(:solvefail, keys_))
