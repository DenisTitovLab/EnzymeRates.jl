# ABOUTME: Prints the census tables from the analyze.jl rows: kcat classes per set, the zero-kcat
# ABOUTME: analysis, sup > kcat counts with classes and examples, gauge failures; usage: labels...
using Serialization, Printf, Statistics
const KDATA = "/tmp/claude-501/-home-denis-linux--julia-dev-EnzymeRates/12f42995-f5ba-4fda-9510-dd9b87f7d069/scratchpad/round4/kcat/data"
const NEX = parse(Int, get(ENV, "NEX", "3"))

rows = [x for l in ARGS for x in deserialize(joinpath(KDATA, "rows_$l.jls"))]
pop(set) = startswith(set, "init_terter") ? "init_terter" : set
sets = unique(pop(x.set) for x in rows)
ok = filter(x -> x.ok, rows)
showex(x) = (println("    ", x.id, " (", x.set, "/", x.index, ")"); foreach(l -> println("      ", l), x.lines))

println("== Task 1: kcat class per set (core = the slice's kcat, dead-end branches removed; full = whole law)")
@printf("%-18s %7s %7s | %7s %6s %6s %6s %6s | full: %7s %6s %6s %6s %6s\n", "set", "records",
        "derived", "finite", "zero", "infin", "nouniq", "other", "finite", "zero", "infin", "nouniq", "other")
classes = (:finite, :zero, :infinite, :nounique, :other)
for set in [sets; "ALL"]
    xs = set == "ALL" ? rows : filter(x -> pop(x.set) == set, rows)
    ys = filter(x -> x.ok, xs)
    c = [count(x -> x.cls === k, ys) for k in classes]
    f = [count(x -> x.fcls === k, ys) for k in classes]
    @printf("%-18s %7d %7d | %7d %6d %6d %6d %6d |       %7d %6d %6d %6d %6d\n", set, length(xs),
            length(ys), c..., f...)
end
others = filter(x -> x.cls === :other || x.fcls === :other, ok)
isempty(others) || println("  other messages: ", unique([x.cls === :other ? x.msg : x.fmsg for x in others]))
notok = filter(x -> !x.ok, rows)
println("  not derivable: ", length(notok), "; reasons: ",
        sort!(collect(Dict(r => count(x -> first(x.why, 40) == r, notok) for r in unique(first(x.why, 40) for x in notok))); by = last, rev = true))
cross = Dict{Tuple{Symbol, Symbol}, Int}()
foreach(x -> cross[(x.cls, x.fcls)] = get(cross, (x.cls, x.fcls), 0) + 1, ok)
println("  core class => full-law class: ", sort!(collect(cross); by = string))

println("\n== Task 2: kcat = 0 laws")
zs = filter(x -> x.cls === :zero, ok)
for set in [sets; "ALL"]
    xs = set == "ALL" ? zs : filter(x -> pop(x.set) == set, zs)
    isempty(xs) && continue
    z2 = [p.z2 for x in xs for p in x.pts]
    pos = count(x -> all(p -> p.z2.npos == 20, x.pts), xs)
    sym = count(x -> x.zp !== nothing && x.zp.qkeeps && x.zp.positive, xs)
    nsym = count(x -> x.zp !== nothing, xs)
    # inhibiting substrates: numeric (all three points agree) and symbolic (plain laws)
    drv(x) = [all(p -> p.z2.drives[j], x.pts) for j in 1:x.nsubs]
    one = count(x -> count(drv(x)) == 1, xs)
    alltog = count(x -> all(p -> p.z2.drives[end], x.pts), xs)
    agree = count(x -> x.zp === nothing || drv(x) == (x.zp.gap .< 0), xs)
    agreec = count(x -> x.zpc === nothing || drv(x) == (x.zpc.gap .< 0), xs)
    traps = [t for x in xs for p in x.pts for t in p.z2.traps]
    @printf("%-18s laws %4d | v>0 at 20 random points x 3 params: %4d | Q(P=0) nonzero & coeffs > 0: %4d of %4d plain\n",
            set, length(xs), pos, sym, nsym)
    @printf("%-18s   one substrate drives v->0: %4d; all together -> 0: %4d; numeric = symbolic (full law) %4d, (core law) %4d\n",
            "", one, alltog, agree, agreec)
    @printf("%-18s   traps %5d: hold X %5d, hold product/residual %5d, exits only X-release or product-binding %5d; share >= 0.999: %5d; oracle worst %.1e\n",
            "", length(traps), count(t -> t.holdsX, traps), count(t -> t.holdsP, traps),
            count(t -> isempty(t.other), traps), count(t -> t.share >= 0.999, traps),
            maximum(p.orc for p in z2))
    kinds = Dict{String, Int}()
    for x in xs, (j, p) in enumerate(x.pts), t in p.z2.traps
        k = !isempty(t.other) ? "other exits" : !t.holdsX ? "trap without X" :
            t.holdsP ? "X + product/residual, " * (isempty(t.blocked) ? "no product-binding exit" : "left only by X release or product binding") :
            "X without product (e.g. second X copy), " * (isempty(t.blocked) ? "left only by X release" : "X release or product binding")
        kinds[k] = get(kinds, k, 0) + 1
    end
    println("                     trap kinds (laws x points x inhibiting substrates): ", sort!(collect(kinds); by = last, rev = true))
    nsub = Dict{String, Int}()
    foreach(x -> (k = join([String(t.X) for t in x.pts[1].z2.traps], "+"); nsub[k] = get(nsub, k, 0) + 1), xs)
    println("                     inhibiting substrate(s) (point 1): ", sort!(collect(nsub)))
    pat = Dict{String, Int}()
    foreach(t -> pat[join(t.forms, " + ")] = get(pat, join(t.forms, " + "), 0) + 1, traps)
    println("                     trap form sets (count over laws x points), top: ",
            first(sort!(collect(pat); by = last, rev = true), 8))
end
bad = filter(x -> any(p -> any(t -> !isempty(t.other) || !t.holdsX, p.z2.traps) || p.z2.npos < 20, x.pts), zs)
println("  laws with a trap exit other than X release / product binding, a trap without X, or v <= 0 somewhere: ", length(bad))
for x in bad[1:min(end, 4)]
    showex(x)
    for p in x.pts, t in p.z2.traps
        println("      X=", t.X, " trap ", t.forms, " share ", round(t.share; sigdigits = 4), " other exits ", t.other, " npos ", p.z2.npos)
    end
end
println("  examples:")
for set in sets
    xs = filter(x -> pop(x.set) == set, zs)
    for x in xs[1:min(end, NEX)]
        showex(x)
        p = x.pts[1].z2
        println("      slopes (log10 v(1e8)/v(1e7)) per substrate then all: ", round.(p.slopes; digits = 2),
                "; v at 20 random points in [", @sprintf("%.3g", p.vmin), ", ", @sprintf("%.3g", p.vmax), "]")
        x.zp === nothing || println("      symbolic gap deg_X(S*Q) - deg_X(D) per substrate (full law): ", x.zp.gap, ", core law: ", x.zpc.gap)
        for t in p.traps
            println("      X = ", t.X, ": trap ", t.forms, " (share ", @sprintf("%.6f", t.share), "), X-release exits ", t.nx,
                    ", exits needing a product: ", t.blocked)
        end
    end
end

println("\n== Task 3: sup v/E_total at zero products (box 1e-4..1e4) against kcat, finite kcat only")
fin = filter(x -> x.cls === :finite, ok)
label(x) = x.allo ? "MWC n=$(x.n)" :
           x.feat.ss_nt_sub > 0 && x.feat.ss_nt_prod > 0 ? "SS cycles: substrate orders + product forms" :
           x.feat.ss_nt_sub > 0 ? "SS cycles: substrate-binding orders" :
           x.feat.ss_nt_prod > 0 ? "SS cycles: product-holding forms only" :
           x.feat.nt > 0 ? "non-turnover cycles all RE" : "no non-turnover cycle"
for set in [sets; "ALL"]
    xs = set == "ALL" ? fin : filter(x -> pop(x.set) == set, fin)
    isempty(xs) && continue
    over = filter(x -> x.ratio > 1.01, xs)
    overs = filter(x -> x.ratios > 1.01, xs)
    @printf("%-18s finite %5d  sup > 1.01 kcat: %5d (%.1f%%), substrates only (modifiers 0): %5d; max ratio %.3f; sup < 0.5 kcat %5d (of them full-law limit < 0.5 kcat: %d)\n",
            set, length(xs), length(over), 100 * length(over) / length(xs), length(overs),
            maximum(x -> x.ratio, xs), count(x -> x.ratio < 0.5, xs),
            count(x -> x.ratio < 0.5 && any(p -> p.kf < 0.5 * p.kc, x.pts), xs))
end
println("  by class (finite laws: total / sup > 1.01 kcat / substrates only):")
for k in sort!(unique(label.(fin)))
    xs = filter(x -> label(x) == k, fin)
    @printf("    %-46s %6d / %5d / %5d\n", k, length(xs), count(x -> x.ratio > 1.01, xs),
            count(x -> x.ratios > 1.01, xs))
end
println("  sup > 1.01 kcat: with modifiers in the mechanism ", count(x -> x.ratio > 1.01 && !isempty(x.mods), fin),
        ", without ", count(x -> x.ratio > 1.01 && isempty(x.mods), fin))
println("  full-law limit kf vs core kcat kc among finite laws: kf < 0.99 kc at some point: ",
        count(x -> any(p -> p.kf < 0.99 * p.kc, x.pts), fin), "; kf > 1.01 kc: ",
        count(x -> any(p -> p.kf > 1.01 * p.kc, x.pts), fin), "; full class: ",
        sort!(collect(Dict(k => count(x -> x.fcls === k, fin) for k in unique(x.fcls for x in fin)))))
println("  ratio distribution (finite): ", [q => round(quantile([x.ratio for x in fin], q); digits = 4) for q in (0.01, 0.1, 0.5, 0.9, 0.99)])
over = sort(filter(x -> x.ratio > 1.01, fin); by = x -> -x.ratio)
println("  examples (largest ratios, then one per class):")
picked = over[1:min(end, NEX)]
for k in unique(label.(over))
    y = findfirst(x -> label(x) == k && !(x in picked), over)
    y === nothing || push!(picked, over[y])
end
for x in picked
    showex(x)
    i = argmax([p.sup / p.kc for p in x.pts])
    p = x.pts[i]
    println("      class: ", label(x), "; sup ", @sprintf("%.4g", p.sup), " kcat ", @sprintf("%.4g", p.kc),
            " ratio ", @sprintf("%.3f", p.sup / p.kc), " (substrates only ", @sprintf("%.3f", p.sups / p.kc),
            "), full-law limit ", @sprintf("%.4g", p.kf), "; at (concentration index => value) ",
            sort!([k => round(v; sigdigits = 3) for (k, v) in p.at]))
end

println("\n== Task 4: reference-rate gauge (50 random concentration vectors x 3 parameter points per law)")
for set in [sets; "ALL"]
    xs = set == "ALL" ? ok : filter(x -> pop(x.set) == set, ok)
    nn = count(x -> any(p -> p.nnan > 0, x.pts), xs)
    ni = count(x -> any(p -> p.ninf > 0, x.pts), xs)
    nz = count(x -> any(p -> !(0 < p.maxabs < Inf), x.pts), xs)
    nz0 = count(x -> any(p -> !(0 < p.maxabs0 < Inf), x.pts), xs)
    @printf("%-18s laws %6d  NaN %4d  Inf %4d  max|v| not finite>0: %4d; over the 10 product-free vectors: %4d; min|v|/max|v| median %.2e\n",
            set, length(xs), nn, ni, nz, nz0, median([p.minabs / p.maxabs for x in xs for p in x.pts]))
end
fails = filter(x -> any(p -> p.nnan > 0 || p.ninf > 0 || !(0 < p.maxabs < Inf), x.pts), ok)
for x in fails[1:min(end, 5)]
    showex(x)
    println("      ", [(p.nnan, p.ninf, p.maxabs, p.maxabs0) for p in x.pts])
end
