# ABOUTME: Summarizes the unit-rule census rows (data/<pop>_*.jls from census.jl): rule kinds,
# ABOUTME: Cleland and order-variant agreement, units, value errors, artifacts, sizes; args: pops.
using Serialization, Printf, Statistics
include(joinpath(@__DIR__, "ulib.jl"))   # the row types (Norm, QQ) must be defined to read

const DATA = joinpath(@__DIR__, "data")

function load(pop)
    files = sort(filter(f -> startswith(f, pop * "_") && endswith(f, ".jls"), readdir(DATA)))
    rows = reduce(vcat, [deserialize(joinpath(DATA, f)) for f in files]; init = Any[])
    files, rows
end

cnt(xs) = sort!(collect(pairs(Dict(k => count(==(k), xs) for k in unique(xs)))); by = string)
pct(a, b) = b == 0 ? "-" : @sprintf("%.1f%%", 100a / b)
fmax(xs) = isempty(xs) ? NaN : maximum(xs)

for pop in ARGS
    files, R = load(pop)
    println("\n================ ", pop, ": ", length(R), " rows from ", length(files), " files")
    bad = filter(r -> !r.ok, R)
    println("errors: ", length(bad))
    for (msg, n) in cnt([first(r.err, 140) for r in bad])
        println("  ", n, " × ", msg)
    end
    ok = filter(r -> r.ok, R)
    isempty(ok) && continue
    println("allosteric laws: ", count(r -> r.allo, ok), "; n > 1: ", count(r -> r.n > 1, ok))
    C = [(r.id, t, x) for r in ok for (t, x) in r.rules]
    println("conformations: ", length(C), " (", join(("$k $n" for (k, n) in cnt([c[2] for c in C])), ", "), ")")
    println("spec rule kinds: ", join(("$k $n" for (k, n) in cnt([c[3].spec_kind for c in C])), ", "))
    println("slice rule kinds: ", join(("$k $n" for (k, n) in cnt([c[3].slice_kind for c in C])), ", "))
    NC = filter(c -> !c[3].constant, C)
    println("no constant term: ", length(NC), " conformations (", pct(length(NC), length(C)), ")")
    cle = filter(c -> c[3].cleland, NC)
    println("  Cleland's monomials all present: ", length(cle))
    diffC = filter(c -> c[3].same_cleland !== true, cle)
    println("  spec N == Cleland N: ", length(cle) - length(diffC), "; differs: ", length(diffC))
    for c in first(diffC, 12)
        x = c[3]
        println("    ", first(c[1], 70), " [", c[2], "] spec ", x.spec_kind, " λ: ", x.spec_lambda,
                " | Cleland λ: ", x.cleland_lambda)
    end
    sdiff = filter(c -> c[3].same_slice !== true, NC)
    println("  slice N == spec N: ", length(NC) - length(sdiff), "; differs: ", length(sdiff),
            " (", join(("$k $n" for (k, n) in cnt([(c[3].slice_kind, c[3].spec_kind) for c in sdiff])), ", "), ")")
    for o in ORDERS
        d = count(c -> c[3].same_order[o] !== true, NC)
        dk = cnt([c[3].kind_order[o] for c in NC if c[3].same_order[o] !== true])
        ce = count(c -> c[3].cle_order[o] === true, cle)
        println(@sprintf("  order %-9s: N differs from spec rule in %d (%s); equals Cleland in %d of %d",
                         o, d, join(("$k $n" for (k, n) in dk), ", "), ce, length(cle)))
    end
    frac = filter(c -> c[3].spec_kind === :virtual_frac, C)
    rat = filter(c -> c[3].spec_kind === :ratio, C)
    for (label, xs) in (("fractional λ", frac), ("ratio form", rat))
        println(label, ": ", length(xs))
        for c in first(xs, 6)
            println("    ", first(c[1], 70), " [", c[2], "] λ: ", c[3].spec_lambda)
        end
    end
    println("units (spec rule): full fail ", count(r -> !r.units_full, ok), ", coefficient fail ",
            count(r -> !r.units_coef, ok))
    for r in first(filter(r -> !(r.units_full && r.units_coef), ok), 5)
        println("    ", first(r.id, 60), ": ", first(r.units_note, 200))
    end
    println("units (slice printer): full fail ", count(r -> !r.units_sfull, ok),
            ", coefficient fail ", count(r -> !r.units_scoef, ok))
    for r in first(filter(r -> !(r.units_sfull && r.units_scoef), ok), 5)
        println("    ", first(r.id, 60), ": ", first(r.units_snote, 200))
    end
    for f in (:e_full, :e_coef, :e_full2, :e_coef2, :e_sfull, :e_scoef, :e_ref, :e_ref2)
        xs = [getfield(r, f) for r in ok]
        @printf("value %-8s max %.1e, over 1e-100: %d\n", f, fmax(xs), count(>(1e-100), xs))
    end
    println("reference-line exponent denominators: ", join(("$k: $n" for (k, n) in cnt([r.refden for r in ok])), ", "))
    println("plain numerators: one concentration monomial ", count(r -> !r.allo && r.nq == 1, ok),
            ", more ", count(r -> !r.allo && r.nq > 1, ok))
    println("artifacts: spec layout ", sum(r -> r.art, ok), ", slice layout ", sum(r -> r.art2, ok))
    a2 = filter(r -> r.art2 > 0, ok)
    isempty(a2) || println("  slice-layout artifacts in: ", join(first([r.id for r in a2], 8), ", "))
    println("sqrt in full law: ", count(r -> r.sqrt, ok), "; ⁻ in coefficient names: ",
            count(r -> r.ratio_names, ok), "; N printed apart: ", count(r -> r.Napart, ok))
    for f in (:len_full, :len_coef, :len_full2, :len_sfull, :len_scoef)
        xs = [getfield(r, f) for r in ok]
        @printf("length %-9s median %6.0f  max %7d\n", f, median(xs), maximum(xs))
    end
    s = [r.secs for r in R]
    @printf("seconds per law: median %.3f, max %.2f, total %.0f\n", median(s), maximum(s), sum(s))
end
