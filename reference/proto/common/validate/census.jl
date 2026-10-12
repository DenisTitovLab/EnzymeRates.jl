# ABOUTME: Census of every population: counts, failures, Pass-2 folds, name partitions and costs.
# ABOUTME: Stdlib only; reads data/*.jls through Population.jl.
include(joinpath(@__DIR__, "..", "Population.jl"))
using Printf, Statistics

q(v, p) = isempty(v) ? NaN : quantile(v, p)
regnames(m) = Set(Symbol("K_", st, String(l), "reg") for site in m.sites
                  for l in site.ligands for st in ("", "A_", "I_"))

for name in Pop.populations()
    pop = Pop.load_population(name)
    plain = [m for m in pop if !m.allosteric]
    allo = [m for m in pop if m.allosteric]
    println("== $name: $(length(pop)) mechanisms (plain $(length(plain)), allosteric ",
            "$(length(allo))); not ok $(count(m -> !m.ok, pop))")
    stages = Dict{Symbol, Int}()
    for m in pop, (k, _) in m.stage_errors; stages[k] = get(stages, k, 0) + 1; end
    isempty(stages) || println("   stage errors: ", join(("$k=$v" for (k, v) in stages), ", "))
    for m in pop
        m.ok || println("   NOT OK $(m.id) $(first(m.source, 50)): $(first(m.error, 200))")
    end
    nk = count(m -> any(s -> !isnan(s.kcat), m.samples), pop)
    println("   kcat computed on $nk; kcat threw on ",
            count(m -> any(e -> first(e) === :kcat, m.stage_errors), pop))
    mism = count(m -> m.ok && m.independent_param_count != length(m.fitted_params), pop)
    println("   _independent_param_count != length(fitted_params): $mism")
    if !isempty(plain)
        p2 = count(m -> !isempty(m.pass2_rename), plain)
        p2c = count(m -> m.pass2_changes_polys, plain)
        rc = count(m -> Pop.re_cycle_closure(m)[1] > 0, plain)
        closure = [Pop.re_cycle_closure(m) for m in plain]
        rc_open = count(c -> c[2] > 0, closure)
        bad = 0
        for m in plain
            m.ok || continue
            names = Set(n for ns in m.step_names for n in ns)
            fitted, dep = Set(m.fitted_params), Set(first.(m.dependents))
            absorbed = Set(first.(m.pass2_rename))
            parts = (fitted, dep, absorbed)
            ok = union(parts...) == names &&
                 sum(length, parts) == length(union(parts...))
            ok || (bad += 1; bad <= 5 && println("   partition broken: $(m.id) ",
                   "fitted∖names=$(setdiff(fitted, names)) names∖all=",
                   "$(setdiff(names, union(parts...)))"))
        end
        println("   plain: Pass-2 fold present $p2, fold changes the polynomials $p2c, ",
                "RE steps closing a cycle $rc (with a cycle the grouping leaves open: ",
                "$rc_open; open RE cycles in total $(sum(last, closure))), names not ",
                "partitioned into fitted/dependent/absorbed $bad")
        for m in plain
            isempty(m.pass2_rename) && continue
            println("     Pass-2 $(m.id): ", join(("$a=>$b" for (a, b) in m.pass2_rename), ", "),
                    m.pass2_changes_polys ? "  (changes polynomials)" : "")
        end
    end
    if !isempty(allo)
        phantom = 0; renamed = 0; reversed = 0; p2A = 0; p2I = 0; opp = 0
        for m in allo
            m.ok || continue
            declared = union(Set(n for ns in [m.step_names_A; m.step_names_I] for n in ns),
                             regnames(m), Set([:L]))
            any(p -> p ∉ declared, m.fitted_params) && (phantom += 1)
            if !isempty(m.step_names_I_pkg)
                any(j -> m.inactive_kept[j] && m.step_names_I_pkg[j] != m.step_names_I[j],
                    eachindex(m.steps)) && (renamed += 1)
                any(m.inactive_reversed) && (reversed += 1)
            end
            isempty(m.pass2_rename_A) || (p2A += 1)
            isempty(m.pass2_rename_I) || (p2I += 1)
        end
        println("   allosteric: fitted name outside the declared constants $phantom; ",
                "inactive step named differently by the package's pruned graph $renamed; ",
                "pruned graph re-orients a step $reversed; Pass-2 fold in A $p2A, in I $p2I")
        tags = Dict{Symbol, Int}()
        for m in allo, t in unique(m.group_tags); tags[t] = get(tags, t, 0) + 1; end
        println("   tags (mechanisms carrying): ", join(("$k=$v" for (k, v) in tags), ", "),
                "; multiplicities ", join(sort(unique(m.cat_multiplicity for m in allo)), ","),
                "; with sites $(count(m -> !isempty(m.sites), allo))")
    end
    ok = [m for m in pop if m.ok]
    ctor = [m.ctor_ns for m in ok if !isnan(m.ctor_ns)]
    ev = [m.eval_ns for m in ok if !isnan(m.eval_ns)]
    fc = [m.first_call_s for m in ok if !isnan(m.first_call_s)]
    @printf("   ctor_ns median %.0f p90 %.0f max %.0f | eval_ns median %.1f p90 %.1f max %.1f \
            (>120 ns: %d; allocating: %d) | first call s median %.3f max %.2f\n",
            q(ctor, 0.5), q(ctor, 0.9), isempty(ctor) ? NaN : maximum(ctor), q(ev, 0.5),
            q(ev, 0.9), isempty(ev) ? NaN : maximum(ev), count(>(120), ev),
            count(m -> m.eval_allocs > 0, ok), q(fc, 0.5), isempty(fc) ? NaN : maximum(fc))
    nf = [length(m.species) for m in pop]; ns = [length(m.steps) for m in pop]
    nfit = [length(m.fitted_params) for m in ok]
    @printf("   forms %d–%d (median %.0f), steps %d–%d, fitted %d–%d, eq_complexity max %d\n",
            minimum(nf), maximum(nf), median(nf), minimum(ns), maximum(ns),
            isempty(nfit) ? 0 : minimum(nfit), isempty(nfit) ? 0 : maximum(nfit),
            maximum(m.eq_complexity for m in pop))
    fam = Dict{Symbol, Int}()
    for m in pop, f in m.family; fam[f] = get(fam, f, 0) + 1; end
    println("   families: ", join(("$k=$v" for (k, v) in sort!(collect(fam); by = first)), ", "))
    println("   distinct eq_hash among ok: $(length(unique(m.eq_hash for m in ok)))")
end
