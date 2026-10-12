# ABOUTME: With the package (er_oracle project): _kcat_forward at each plain record's first
# ABOUTME: sample, and the package's derivation times, for the 806 shared plain mechanisms.
empty!(ARGS); push!(ARGS, "specs")          # the exporter's prelude loads the spec fixtures
const SCR = "/tmp/claude-501/-home-denis-linux--julia-dev-EnzymeRates/12f42995-f5ba-4fda-9510-dd9b87f7d069/scratchpad"
const EXPORT = joinpath(SCR, "proto", "common", "export", "export_pop.jl")
let src = read(EXPORT, String)
    cut = findfirst("const MEMBERS = members(POP)", src)
    include_string(Main, src[1:first(cut) - 1], EXPORT)
end
const POPS = ["specs", "init_uniuni", "init_unibi", "init_bibi", "init_bibi_pp", "expand_bibi",
              "handwritten"]
best_ns(f; reps = 3) = (f(); minimum(begin t = time_ns(); f(); Float64(time_ns() - t) end
                                     for _ in 1:reps))
out = Dict{String, NamedTuple}()
for pop in POPS
    recs = Dict(r.id => r for r in Pop.load_population(pop))
    n = 0
    for (i, (m, src)) in enumerate(members(pop))
        em = ER.compile_mechanism(m)
        lm = ER._concrete(em)
        lm isa ER.AllostericMechanism && continue
        rec = recs["$pop/$i"]
        (rec.source == src && !rec.allosteric) || error("member mismatch at $pop/$i")
        params = (; (k => v for (k, v) in rec.samples[1].params)...)
        kcat, msg = try
            Base.invokelatest(ER._kcat_forward, em, params), ""
        catch e
            NaN, sprint(showerror, e)
        end
        t_polys = best_ns(() -> ER._raw_symbolic_rate_polys(lm))
        t_dep = best_ns(() -> ER._dependent_param_exprs(lm))
        out[rec.id] = (; kcat, msg, t_polys, t_dep)
        n += 1
    end
    println(pop, ": ", n, " plain mechanisms"); flush(stdout)
end
serialize(joinpath(SCR, "slice", "engine", "data", "pkg_kcat_times.jls"), out)
println("done: ", length(out))
