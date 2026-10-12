# ABOUTME: Follow-up B1: two real examples per R_D class (DSL-style steps with tags), the structural classes
# ABOUTME: of P5a's hand-built cases against their (i)/(ii) verdicts, and the cost of the structural check.
include(joinpath(@__DIR__, "b1_lib.jl"))
include(joinpath(S, "proto/p5a/limits_lib.jl"))      # retag() for P5a's rejected uni-uni twins

pops = Dict(name => load_allo(name) for name in POPS)
allm = [m for name in POPS for m in pops[name]]

# ── examples: smallest mechanisms of each class, phantom with one kinetic door first ──
function pick(pred, n)
    c = sort(filter(pred, allm); by = m -> (length(m.steps), m.id))
    c[1:min(end, n)]
end
println("######## R_D examples")
ex = Dict(
    "no_island" => pick(m -> :OnlyA in m.group_tags && analyze(m, :D).cls == "no_island", 2),
    "emptied" => pick(m -> (a = analyze(m, :D); a.cls == "emptied" && a.nkdoor == 1 && !isempty(a.forced)), 1) ∪
                 pick(m -> (a = analyze(m, :D); a.cls == "emptied" && a.vanish_cycle), 1),
    "phantom" => pick(m -> (a = analyze(m, :D); a.cls == "phantom" && a.nkdoor == 1), 1) ∪
                 pick(m -> (a = analyze(m, :D); a.cls == "phantom" && a.nkdoor == 2 && a.nadoor >= 1), 1),
    "undefined" => pick(m -> (a = analyze(m, :D); a.cls == "undefined" && a.nadoor == 0), 1) ∪
                   pick(m -> (a = analyze(m, :D); a.cls == "undefined" && a.nadoor >= 1 && a.nisl == 1), 1))
for c in ("no_island", "emptied", "phantom", "undefined")
    println("\n== class $c")
    for m in ex[c]; print_mech(stdout, m, :D); end
end

# ── P5a's hand-built cases (uni-uni grid, ping-pong traps): class vs P5a (i)/(ii) verdict ──
println("\n######## P5a hand-built cases (cases.jls): structural class vs P5a verdict")
cases = [retag(m) for m in deserialize(joinpath(S, "proto/p5a/data/cases.jls"))]
recs = Dict((r[:id], r[:conv]) => r for r in deserialize(joinpath(S, "proto/p5a/data/results_cases.jls")))
for (rd, conv) in ((:iu, :i), (:D, :ii))
    tab = Dict{Tuple{String, String}, Int}()
    for m in cases
        a = analyze(m, rd)
        r = get(recs, (m.id, conv), nothing)
        adm = r === nothing ? "?" : (r[:admitted] ? "adm" : "REJ")
        v = verdict(r)
        tab[(a.cls, "$adm $v")] = get(tab, (a.cls, "$adm $v"), 0) + 1
    end
    println("  R_$rd vs ($conv) [(i) ≡ (iu) for single-step groups]:")
    for ((c, v), n) in sort!(collect(tab); by = x -> (x[1][1], x[1][2]))
        @printf("    %-13s %-16s %3d\n", c, v, n)
    end
end
println("\n  hand-built non-grid cases under R_D / R_iu:")
for m in cases
    startswith(m.id, "uni/") && continue
    aD, aU = analyze(m, :D), analyze(m, :iu)
    vs = join(["($c) " * verdict(get(recs, (m.id, c), nothing)) for c in (:i, :ii)], " ")
    println("    ", rpad(m.id, 18), " R_D=", rpad(aD.cls, 13), " R_iu=", rpad(aU.cls, 13), " ", vs)
end
for id in ("pp_err1", "pp_dead")
    m = only(filter(x -> x.id == id, cases))
    print_mech(stdout, m, :D)
end

# ── cost of the structural check ──
println("\n######## cost")
onlyA = filter(m -> :OnlyA in m.group_tags, allm)
for m in onlyA; analyze(m, :D); end
for rd in (:D, :iu)
    t = minimum(@elapsed(for m in onlyA; analyze(m, rd); end) for _ in 1:5)
    @printf("  R_%-2s: %d mechanisms in %.3f s = %.1f µs per mechanism (mean steps %.1f, forms %.1f)\n", rd,
            length(onlyA), t, 1e6t / length(onlyA), sum(length(m.steps) for m in onlyA) / length(onlyA),
            sum(length(m.species) for m in onlyA) / length(onlyA))
end
big = sort(onlyA; by = m -> -length(m.steps))[1:20]
t = minimum(@elapsed(for _ in 1:100, m in big; analyze(m, :D); end) for _ in 1:3)
@printf("  20 largest (%d–%d steps): %.1f µs per mechanism\n", length(big[end].steps), length(big[1].steps), 1e6t / 2000)
