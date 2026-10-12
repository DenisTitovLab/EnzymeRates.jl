# ABOUTME: Times the package's _raw_symbolic_rate_polys on every plain mechanism of the probe populations.
# ABOUTME: Rebuilds the exporter's members (same order), checks them against the Pop records, writes pkg_times.jls.
empty!(ARGS); push!(ARGS, "specs")          # the exporter's prelude loads the spec fixtures for "specs"
const EXPORT = joinpath(@__DIR__, "..", "common", "export", "export_pop.jl")
let src = read(EXPORT, String)
    cut = findfirst("const MEMBERS = members(POP)", src)
    include_string(Main, src[1:first(cut) - 1], EXPORT)
end
const POPS = ["specs", "init_uniuni", "init_unibi", "init_bibi", "init_bibi_pp", "expand_bibi",
              "handwritten"]
topoly(p) = Dict{Vector{Pair{Symbol, Int}}, Rational{Int}}(collect(k) => c for (k, c) in p)
function best_ns(f; reps = 3)
    f()
    minimum(begin t = time_ns(); f(); Float64(time_ns() - t) end for _ in 1:reps)
end
times = Dict{String, NamedTuple}()
for pop in POPS
    recs = Dict(r.id => r for r in Pop.load_population(pop))
    ms = members(pop)
    nplain = 0
    for (i, (m, src)) in enumerate(ms)
        lm = ER._concrete(ER.compile_mechanism(m))
        lm isa ER.AllostericMechanism && continue
        id = "$pop/$i"
        rec = recs[id]
        (rec.source == src && !rec.allosteric) || error("member mismatch at $id")
        sp = ER._step_parameters(lm)
        num, den, _ = ER._raw_symbolic_rate_polys(lm, sp, Dict{Symbol, Symbol}())
        (topoly(num), topoly(den)) == (topoly(rec.raw_num), topoly(rec.raw_den)) || error("polys differ at $id")
        t_raw = best_ns(() -> ER._raw_symbolic_rate_polys(lm, ER._step_parameters(lm),
                                                          Dict{Symbol, Symbol}()))
        t_default = best_ns(() -> ER._raw_symbolic_rate_polys(lm))
        t_dep = best_ns(() -> ER._dependent_param_exprs(lm))
        times[id] = (; t_raw, t_default, t_dep)
        nplain += 1
    end
    println(pop, ": timed ", nplain, " plain mechanisms"); flush(stdout)
end
serialize(joinpath(@__DIR__, "pkg_times.jls"), times)
t(f) = sort!([getfield(v, f) for v in values(times)])
for f in (:t_raw, :t_default, :t_dep)
    v = t(f)
    println(f, ": n=", length(v), " median ", round(v[(end + 1) ÷ 2] / 1e3; digits = 1), " µs, max ",
            round(v[end] / 1e3; digits = 1), " µs, total ", round(sum(v) / 1e9; digits = 3), " s")
end
