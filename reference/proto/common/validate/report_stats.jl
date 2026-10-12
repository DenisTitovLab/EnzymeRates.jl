# ABOUTME: Extra statistics for the Foundation report: slow evaluations, phantom names, folds.
# ABOUTME: Stdlib only.
include(joinpath(@__DIR__, "..", "Population.jl"))
using Printf, Statistics
all = [m for n in Pop.populations() for m in Pop.load_population(n)]
println("total mechanisms: ", length(all), " (plain ", count(m -> !m.allosteric, all),
        ", allosteric ", count(m -> m.allosteric, all), ")")
for m in sort([m for m in all if m.eval_ns > 60]; by = m -> -m.eval_ns)[1:min(8, end)]
    @printf("slow eval %-20s %6.1f ns  eq_complexity %4d  den terms %d  fitted %d  %s\n", m.id,
            m.eval_ns, m.eq_complexity, length(m.raw_den), length(m.fitted_params), m.family)
end
ec = [(m.eq_complexity, m.eval_ns) for m in all if !isnan(m.eval_ns)]
println("corr(eq_complexity, eval_ns) = ", round(cor(first.(ec), last.(ec)); digits = 3))
regnames(m) = Set(Symbol("K_", st, String(l), "reg") for site in m.sites
                  for l in site.ligands for st in ("", "A_", "I_"))
k = 0
for m in all
    m.allosteric || continue
    declared = union(Set(n for ns in [m.step_names_A; m.step_names_I] for n in ns), regnames(m),
                     Set([:L]))
    ph = [p for p in m.fitted_params if p ∉ declared]
    isempty(ph) && continue
    global k += 1
    k <= 3 || continue
    println("phantom fitted names in $(m.id): $ph")
    for (j, s) in enumerate(m.steps)
        m.inactive_kept[j] && m.step_names_I_pkg[j] != m.step_names_I[j] &&
            println("    step $j $(Pop.step_text(m, s)) [$(m.group_tags[Pop.group_of_step(m)[j]])]: ",
                    "whole-mechanism I name $(m.step_names_I[j]), package I name $(m.step_names_I_pkg[j])")
    end
end
for m in all
    m.allosteric && (!isempty(m.pass2_rename_A) || !isempty(m.pass2_rename_I)) || continue
    opp = [(a, b) for (a, b) in m.pass2_rename_A for (c, d) in m.pass2_rename_I
           if replace(String(a), "A_" => "") == replace(String(d), "I_" => "") &&
              replace(String(b), "A_" => "") == replace(String(c), "I_" => "")]
    isempty(opp) || println("opposite-direction folds in $(m.id): A $(m.pass2_rename_A) I $(m.pass2_rename_I)")
end
println("allosteric with a Pass-2 fold in exactly one state: ",
        count(m -> m.allosteric && xor(isempty(m.pass2_rename_A), isempty(m.pass2_rename_I)), all))
fc = [m.first_call_s for m in all if !isnan(m.first_call_s)]
@printf("first rate_equation call: median %.3f s, p99 %.3f s, max %.3f s; total %.0f s\n",
        median(fc), quantile(fc, 0.99), maximum(fc), sum(fc))
ct = [m.ctor_ns for m in all]
@printf("constructor _with: plain median %.1f µs, allosteric median %.1f µs\n",
        median([m.ctor_ns for m in all if !m.allosteric]) / 1e3,
        median([m.ctor_ns for m in all if m.allosteric]) / 1e3)
kfail = [(m.id, m.source) for m in all if any(e -> first(e) === :kcat, m.stage_errors)]
println("kcat threw on: ", kfail)
