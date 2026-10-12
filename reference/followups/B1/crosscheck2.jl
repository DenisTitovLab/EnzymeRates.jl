# ABOUTME: Follow-up B1 final cross-check: structural class (R_D, R_iu) × numeric verdict, against P5a's
# ABOUTME: 16-direction records and the full-direction re-evaluation (ii, iu, iiu, Dk); lists disagreements.
include(joinpath(@__DIR__, "b1_lib.jl"))
const CLASSES = ["no_island", "emptied", "phantom", "undefined", "kinetic_island", "sequestering", "ambiguous", "infeasible", "free_forms"]
const VERDS = ["U=D", "U=D+seq", "U≠D", "U≠D+seq", "NU", "NU+seq", "empty", "err", "—"]

"""Full-direction records of population `name` keyed by (id, conv)."""
function load_full(name)
    d = Dict{Tuple{String, Symbol}, Dict{Symbol, Any}}()
    dir = joinpath(@__DIR__, "data")
    for f in readdir(dir)
        (startswith(f, "full_$(name).jls") || startswith(f, "full_$(name)_")) || continue
        endswith(f, "_ids.jls") && continue
        for r in deserialize(joinpath(dir, f)); d[(r[:id], r[:conv])] = r; end
    end
    for f in readdir(dir)                       # Dka: only mechanisms with kinetic cycle doors
        startswith(f, "dka_") || continue
        for r in deserialize(joinpath(dir, f)); startswith(r[:id], name * "/") && (d[(r[:id], :Dka)] = r); end
    end
    d
end

"""Class under free kinetic ratios (P5a's cone): an island outside the main affinity cluster has a
free mass order (not unique, can sequester): "kinetic_island"; otherwise as the kept-ratio class."""
function cls_free(a)
    a.cls in ("free_forms", "infeasible", "no_island") && return a.cls
    isl = [c for c in 1:a.ncomp if c != a.main]
    inmain = [a.clus[c] == a.clus[a.main] for c in isl]
    any(k -> inmain[k] && a.isl_sign[k] == -1, eachindex(isl)) && return "sequestering"
    any(!, inmain) && return "kinetic_island"
    any(==(2), a.isl_sign) && return "ambiguous"
    "emptied"
end
klass(a, sem) = sem === :kept ? a.cls : cls_free(a)

"""Strict agreement with the hypotheses: no_island/emptied ⇔ U=D; phantom ⇒ U≠D without
sequestration; undefined and kinetic_island ⇒ NU; sequestering ⇒ +seq."""
function agrees(cls, v)
    v == "—" && return true
    cls in ("no_island", "emptied") && return v == "U=D"
    v == "U=D" && return false
    cls == "phantom" && return v == "U≠D"
    cls in ("undefined", "kinetic_island") && return startswith(v, "NU")
    cls == "sequestering" && return endswith(v, "+seq")
    false
end
single_step_onlyA(m) = all(g -> m.group_tags[g] !== :OnlyA || length(m.groups[g]) == 1, eachindex(m.groups))

# (reading, ratio semantics, convention, source)
const PAIRS = [(:D, :kept, :ii, :p5a), (:D, :kept, :Dk, :full), (:D, :kept, :Dka, :full),
               (:D, :free, :ii, :p5a), (:D, :free, :ii, :full), (:D, :free, :iiu, :full),
               (:iu, :kept, :iu, :p5a), (:iu, :free, :iu, :p5a), (:iu, :free, :iu, :full)]
total = Dict{NTuple{4, Symbol}, Dict{Tuple{String, String}, Int}}()
disagree = Dict{NTuple{4, Symbol}, Vector{String}}()
classes = Dict{Tuple{String, Symbol}, Dict{String, Int}}()
kd = Dict{Tuple{Symbol, String}, Vector{Int}}()     # class => [0, 1, ≥2] kinetic doors
flags = Dict{Tuple{String, Symbol}, Vector{Int}}()  # forced, vanishing cycle, cyc_single_form
for name in POPS
    pop = load_allo(name)
    p5a = load_records(name); full = load_full(name)
    for rd in (:D, :iu)
        cnt = get!(classes, (name, rd), Dict{String, Int}())
        fl = get!(flags, (name, rd), zeros(Int, 4))
        for m in pop
            a = analyze(m, rd)
            cnt[a.cls] = get(cnt, a.cls, 0) + 1
            cf = "free:" * cls_free(a); cnt[cf] = get(cnt, cf, 0) + 1
            k = get!(kd, (rd, a.cls), zeros(Int, 3)); k[min(a.nkdoor, 2) + 1] += 1
            fl[1] += !isempty(a.forced); fl[2] += a.vanish_cycle; fl[3] += a.cyc_single_form
            fl[4] += a.cyc_hang
        end
    end
    println("\n######## $name: ", length(pop), " allosteric (", count(m -> :OnlyA in m.group_tags, pop), " with :OnlyA)")
    for (rd, sem, conv, src) in PAIRS
        recs = src === :p5a ? p5a : full
        tab = get!(total, (rd, sem, conv, src), Dict{Tuple{String, String}, Int}())
        bad = get!(disagree, (rd, sem, conv, src), String[])
        loc = Dict{Tuple{String, String}, Int}()
        for m in pop
            a = analyze(m, rd)
            r = get(recs, (m.id, conv), nothing)
            conv === :Dka && r === nothing && (r = get(recs, (m.id, :Dk), nothing))   # no cycle door: Dka ≡ Dk
            if r === nothing && src === :p5a && conv === :iu && single_step_onlyA(m)
                r = get(recs, (m.id, :i), nothing)        # (i) ≡ (iu) when every :OnlyA group has one step
            end
            v = :OnlyA in m.group_tags ? verdict(r) : "U=D"
            c = klass(a, sem)
            loc[(c, v)] = get(loc, (c, v), 0) + 1
            tab[(c, v)] = get(tab, (c, v), 0) + 1
            agrees(c, v) || push!(bad, "$(m.id)[$c→$v]")
        end
        vs = [v for v in VERDS if any(k -> k[2] == v, keys(loc))]
        @printf("  R_%-3s %-4s vs %-4s (%s): %s\n", rd, sem, conv, src, join(rpad.(vs, 8)))
        for c in CLASSES
            any(k -> k[1] == c, keys(loc)) || continue
            @printf("    %-13s %s\n", c, join([rpad(string(get(loc, (c, v), 0)), 8) for v in vs]))
        end
    end
end
println("\n######## classes per population")
for rd in (:D, :iu)
    println("  R_$rd:")
    for name in POPS
        c = classes[(name, rd)]; f = flags[(name, rd)]
        @printf("    %-16s %s | free: %s | forced-ratio kinetic doors %d, vanishing cycles %d, cycle with one attach form %d, hanging cycles %d\n",
                name, join(("$k=$(c[k])" for k in CLASSES if haskey(c, k)), " "),
                join(("$k=$(c["free:" * k])" for k in CLASSES if haskey(c, "free:" * k)), " "), f...)
    end
    println("  kinetic doors per mechanism (0 / 1 / ≥2) by class:")
    for c in CLASSES
        haskey(kd, (rd, c)) && @printf("    %-13s %s\n", c, join(kd[(rd, c)], " / "))
    end
end
println("\n######## totals over all populations")
for (rd, sem, conv, src) in PAIRS
    tab = total[(rd, sem, conv, src)]
    vs = [v for v in VERDS if any(k -> k[2] == v, keys(tab))]
    @printf("  R_%-3s %-4s vs %-4s (%s): %s\n", rd, sem, conv, src, join(rpad.(vs, 8)))
    for c in CLASSES
        any(k -> k[1] == c, keys(tab)) || continue
        @printf("    %-13s %s\n", c, join([rpad(string(get(tab, (c, v), 0)), 8) for v in vs]))
    end
    bad = disagree[(rd, sem, conv, src)]
    println("    disagreements: ", length(bad))
    isempty(bad) || println("      ", join(bad, " "))
end
