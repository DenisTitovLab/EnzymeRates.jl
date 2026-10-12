# ABOUTME: Follow-up B1 cross-check: structural class (readings D, iu, ii, i) × P5a's numeric cone verdict
# ABOUTME: per population; lists every disagreement with the hypotheses; serializes the joined table.
include(joinpath(@__DIR__, "b1_lib.jl"))
const CLASSES = ["no_island", "emptied", "phantom", "undefined", "sequestering", "ambiguous", "infeasible", "free_forms"]
const VERDS = ["U=D", "U=D+seq", "U≠D", "U≠D+seq", "NU", "NU+seq", "empty", "err", "—"]
"""Hypothesis check: the verdict a class predicts (R_iu: no_island/emptied ⇔ U=D; R_D also
phantom ⇔ U≠D and undefined ⇔ NU)."""
function agrees(cls, v, rd)
    v == "—" && return true
    ok_ud = cls in ("no_island", "emptied")
    rd in (:iu, :i) && return ok_ud == (v == "U=D")
    cls in ("no_island", "emptied") && return v == "U=D"
    cls == "phantom" && return startswith(v, "U≠D")
    cls == "undefined" && return startswith(v, "NU")
    cls == "sequestering" && return endswith(v, "+seq")
    false
end
single_step_onlyA(m) = all(g -> m.group_tags[g] !== :OnlyA || length(m.groups[g]) == 1, eachindex(m.groups))
rows = Dict{Symbol, Any}[]
pairs = [(:D, :ii), (:iu, :iu), (:ii, :ii), (:i, :i)]
for name in POPS
    pop = load_allo(name)
    recs = load_records(name)
    println("\n######## $name: ", length(pop), " allosteric")
    for (rd, conv) in pairs
        tab = Dict{Tuple{String, String}, Int}()
        bad = String[]
        for m in pop
            a = analyze(m, rd)
            r = get(recs, (m.id, conv), nothing)
            # (iu) records exist only for P5a's populations and the skeptic's; use (i) where identical
            src = conv
            if r === nothing && conv === :iu && single_step_onlyA(m)
                r = get(recs, (m.id, :i), nothing); src = :i
            end
            v = :OnlyA in m.group_tags ? verdict(r) : "U=D"   # no blocked step: the limit is the model
            tab[(a.cls, v)] = get(tab, (a.cls, v), 0) + 1
            agrees(a.cls, v, rd) || push!(bad, "$(m.id)[$(a.cls)→$v]")
            push!(rows, Dict(:pop => name, :id => m.id, :reading => rd, :conv => conv, :src => src,
                             :cls => a.cls, :verdict => v, :nkdoor => a.nkdoor, :nadoor => a.nadoor,
                             :nisl => a.nisl, :forced => length(a.forced), :agree => agrees(a.cls, v, rd)))
        end
        vs = [v for v in VERDS if any(k -> k[2] == v, keys(tab))]
        @printf("  reading %-3s vs conv (%s):  %s\n", rd, conv, join(rpad.(vs, 8)))
        for c in CLASSES
            any(k -> k[1] == c, keys(tab)) || continue
            @printf("    %-13s %s\n", c, join([rpad(string(get(tab, (c, v), 0)), 8) for v in vs]))
        end
        println("    disagreements with the hypothesis: ", length(bad))
        isempty(bad) || println("      ", join(bad[1:min(end, 40)], " "), length(bad) > 40 ? " …" : "")
    end
end
serialize(joinpath(@__DIR__, "crosscheck_rows.jls"), rows)
