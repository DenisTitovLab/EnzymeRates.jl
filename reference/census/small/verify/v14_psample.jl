# ABOUTME: Re-tests every same-law verdict of the census's P sample (p_sample.log) on the package laws
# ABOUTME: with this verifier's Levenberg–Marquardt family fit, both ways.
include(joinpath(@__DIR__, "vlib.jl"))
rng = Xoshiro(141414)

entries = Tuple{String, Symbol, Int, String, Int}[]       # run, sub-class, gid, verdict, matched gid
let cur = ""
    for ln in eachline(joinpath(CEN, "logs", "p_sample.log"))
        m = match(r"^== (\S+): ", ln)
        m === nothing || (cur = m.captures[1]; continue)
        m = match(r"^\s+\[(\w+)\] #(\d+) .*→ (same law as a smaller candidate|contains a smaller candidate's law, not equal|no smaller candidate with the same coefficient support)(?: #(\d+))?", ln)
        m === nothing && continue
        push!(entries, (cur, Symbol(m.captures[1]), parse(Int, m.captures[2]), m.captures[3],
                        m.captures[4] === nothing ? 0 : parse(Int, m.captures[4])))
    end
end
println("parsed $(length(entries)) verdicts")
tab = Dict{Tuple{Symbol, String, String}, Int}()
cache = Dict{String, Any}()
for (name, sub, g, verdict, h) in entries
    h == 0 && (k = (sub, verdict, "not tested (no comparator named)"); tab[k] = get(tab, k, 0) + 1; continue)
    d, mechs = get!(() -> loadrun(name), cache, name)
    a, b = plaw(mechs[g]), plaw(mechs[h])
    ab, ba = same_family(a, b, rng)
    mine = ab < 1e-8 && ba < 1e-8 ? "same" : ba < 1e-8 ? "P contains comparator" :
           ab < 1e-8 ? "P inside comparator" : "different"
    k = (sub, verdict, mine)
    tab[k] = get(tab, k, 0) + 1
    @printf("  %s [%s] #%d (np %d rank %d) vs #%d (np %d rank %d): census '%s' | package laws: P-in-comparator %.1e, comparator-in-P %.1e → %s\n",
            name, sub, g, d.recs[g][:np], d.recs[g][:rank_pkg], h, d.recs[h][:np], d.recs[h][:rank_pkg],
            verdict, ab, ba, mine)
    flush(stdout)
end
println("\n(sub-class, census verdict, package-law verdict):")
for (k, v) in sort!(collect(tab); by = first)
    println("  $(k[1]) | $(k[2]) | $(k[3]): $v")
end
