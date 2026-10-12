# ABOUTME: Summarizes saved pipeline runs: rejection messages, stage-time totals and medians,
# ABOUTME: sizes, normalizer rules, kcat outcomes; usage: summarize_pipe.jl labels...
include(joinpath(@__DIR__, "lib.jl"))
rows = reduce(vcat, [deserialize(joinpath(IDATA, "pipe_$l.jls")) for l in ARGS])
ok = filter(x -> x.ok, rows)
println(length(rows), " records, ", length(ok), " derived")
msgs = Dict{String, Int}()
foreach(x -> (k = replace(first(x.msg, 70), r"\{[^}]*\}" => "{…}"); msgs[k] = get(msgs, k, 0) + 1),
        filter(x -> !x.ok, rows))
foreach(((k, v),) -> println("  rejected ", v, " × ", k), sort!(collect(msgs); by = last, rev = true))
for f in (:t_derive, :t_solve, :t_law, :t_key, :t_print, :t_code, :t_kcat, :t_oracle)
    v = filter(!isnan, [1e3 * get(x, f, NaN) for x in ok])
    @printf("  %-9s total %8.1f s; %s ms\n", f, sum(v) / 1e3, q3(v))
end
for f in (:nD, :nQ, :ops, :nfit, :printed)
    println("  ", rpad(f, 9), q3(Float64[getfield(x, f) for x in ok]))
end
println("  trees    ", q3(Float64[sum(x.trees) for x in ok]))
println("  rref     ", q3(Float64[prod(x.rref) for x in ok]), " (rows x cols)")
rules = Dict{Any, Int}()
foreach(x -> (rules[x.rules] = get(rules, x.rules, 0) + 1), ok)
println("  rules: ", rules)
kc = Dict{String, Int}()
foreach(x -> (k = x.kcat isa String ? first(x.kcat, 60) : "candidates";
              kc[k] = get(kc, k, 0) + 1), ok)
println("  kcat: ", kc)
println("  oracle worst ", maximum(x -> x.oracle, ok; init = 0.0), "; half exponents ",
        count(x -> x.halves, ok), "; BigInt laws ", count(x -> x.big, ok))
