# ABOUTME: Allosteric check (b): admission verdicts against follow-up B1's structural class and
# ABOUTME: full-direction limits, P5a's and the skeptic's (iu) verdicts, and the package's guard.
include(joinpath(@__DIR__, "lib.jl"))
module B1
import Main.Pop
include(joinpath(Main.SCRATCH, "followups", "B1", "b1_lib.jl"))
end

"""Follow-up B1's class under per-step classes with free kinetic ratios (its crosscheck2
`cls_free` of `analyze(p, :iu)`); admission there is `no_island` or `emptied`."""
function b1_class(p)
    a = B1.analyze(p, :iu)
    a.cls in ("free_forms", "infeasible", "no_island") && return a.cls
    isl = [c for c in 1:a.ncomp if c != a.main]
    inmain = [a.clus[c] == a.clus[a.main] for c in isl]
    any(k -> inmain[k] && a.isl_sign[k] == -1, eachindex(isl)) && return "sequestering"
    any(!, inmain) && return "kinetic_island"
    any(==(2), a.isl_sign) && return "ambiguous"
    "emptied"
end

"""B1's full-direction records of population `name`, keyed by (id, conv)."""
function b1_full(name)
    d = Dict{Tuple{String, Symbol}, Dict{Symbol, Any}}()
    dir = joinpath(SCRATCH, "followups", "B1", "data")
    for f in readdir(dir)
        startswith(f, "full_$(name).jls") || startswith(f, "full_$(name)_") || continue
        endswith(f, "_ids.jls") && continue
        foreach(r -> d[(r[:id], r[:conv])] = r, deserialize(joinpath(dir, f)))
    end
    d
end

"""The skeptic's adversarial-case records keyed by (id, conv)."""
adv = Dict((r[:id], r[:conv]) => r for f in ("adv_results", "adv_results3", "adv_results4")
           for r in deserialize(joinpath(PROTO, "skeptic_P5a-onlyA", "data", f * ".jls")))

kind(msg) = msg === nothing ? "admit" : startswith(msg, "the conformation flips") ? "flip" :
    occursin("still join", msg) ? "G1 chord" : occursin("close a cycle", msg) ? "G1 cycle" :
    occursin("reach the free enzyme", msg) ? "islands" : "other"

const POPS = ("specs", "handwritten", "allo_uniuni", "allo_bibi", "expand_bibi",
              "allo_uniuni_inh", "allo_unibi", "allo_bibi_pp", "allo_terbi")
records = [r for r in load_records() if S.is_allosteric(r.c.mech)]
tab = Dict{Any, Int}()
bump!(k) = (tab[k] = get(tab, k, 0) + 1)
disagree = String[]
full = Dict(n => b1_full(n) for n in POPS)
p5a = Dict(n => B1.load_records(n) for n in POPS)
for r in records
    m, p = r.c.mech, r.p
    msg = S.admit(m)
    mine = kind(msg)
    if r.set in POPS
        cls = b1_class(p)
        b1ok = cls in ("no_island", "emptied")
        fv = B1.verdict(get(full[r.set], (p.id, :iu), nothing))
        pv = B1.verdict(get(p5a[r.set], (p.id, :iu), get(p5a[r.set], (p.id, :i), nothing)))
        bump!(("population", mine, cls, fv))
        (mine == "admit") == b1ok || push!(disagree, "$(p.id): mine $mine, B1 $cls")
        mine == "admit" && fv ∉ ("U=D", "—") && push!(disagree, "$(p.id): admitted, B1 full $fv")
        mine != "admit" && fv == "U=D" && push!(disagree, "$(p.id): rejected, B1 full U=D")
        bump!(("P5a 16 dirs", mine, pv))
    else
        cls = b1_class(p)
        sk = get(adv, (p.id, :iu), nothing)
        skv = sk === nothing ? "—" : !haskey(sk, :cone_v) ? "empty" :
              sk[:cone_spread] >= 1e-20 ? "NU" : sk[:cone_vs_del] isa Number &&
              sk[:cone_vs_del] < 1e-20 ? "U=D" : "U≠D"
        bump!(("hand-built", r.set, mine, cls, skv, p.ok ? "pkg admits" : "pkg rejects"))
        println(rpad("$(r.set) $(p.id)", 52), rpad(mine, 10), rpad(cls, 16), rpad(skv, 7),
                p.ok ? "package admits" : "package rejects: " * first(string(p.error), 60))
    end
end
println("\n== tallies (key => mechanisms)")
foreach(((k, v),) -> println("  ", k, " => ", v), sort!(collect(tab); by = string))
println("\ndisagreements: ", length(disagree))
foreach(x -> println("  ", x), disagree)
