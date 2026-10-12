# ABOUTME: Counts plain SS steps whose two constants never reach the rate polynomials.
# ABOUTME: Stdlib only; such a step adds one fitted direction (its overall speed) that v ignores.
include(joinpath(@__DIR__, "..", "Population.jl"))
syms(p) = Set(s for (mono, _) in p for (s, _) in mono)
for name in Pop.populations()
    plain = [m for m in Pop.load_population(name) if !m.allosteric && m.ok]
    isempty(plain) && continue
    hits = 0; ex = String[]
    for m in plain
        used = union(syms(m.pkg_num), syms(m.pkg_den))
        silent = [j for (j, s) in enumerate(m.steps)
                  if !s.re && all(n -> get(Dict(m.pass2_rename), n, n) ∉ used, m.step_names[j])]
        isempty(silent) && continue
        hits += 1
        fitted = [n for j in silent for n in m.step_names[j] if n in m.fitted_params]
        length(ex) < 4 && push!(ex, "$(m.id) $([Pop.step_text(m, m.steps[j]) for j in silent]) " *
                                     "fitted among them: $fitted  ($(first(m.source, 50)))")
    end
    println("$name: plain $(length(plain)); with an SS step absent from the rate polynomials $hits")
    foreach(e -> println("    ", e), ex)
end
