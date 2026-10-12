# ABOUTME: Measures how the RE spanning-tree (BFS) order changes the Full-mode and Reduced laws, and how
# ABOUTME: often group or dependent substitution changes the concentration lowest-terms shift.
const C = joinpath(@__DIR__, "..", "..", "common")
include(joinpath(C, "Population.jl")); include(joinpath(C, "Poly.jl")); include(joinpath(@__DIR__, "EngineLean.jl"))
using Printf
const POPS = ["specs", "init_uniuni", "init_unibi", "init_bibi", "init_bibi_pp", "expand_bibi", "handwritten"]
const PI = Poly.P{Int}
"""The law in package names: steps → group names (Pass-2 folded); with `dep`, dependents substituted too."""
function named(m, V, sv, num, den; reduced)
    rename = Dict(m.pass2_rename)
    sub = Dict(v => Poly.mono(Poly.var!(V, get(rename, m.step_names[j][i], m.step_names[j][i])) => 1)
               for (j, vs) in enumerate(sv) for (i, v) in enumerate(vs))
    dep = Dict(Poly.var!(V, k) => Poly.mono((Poly.var!(V, s) => Int(e) for (s, e) in mono)...)
               for (k, mono) in m.dependent_monomials)
    f(p) = (q = Poly.subst(p, v -> get(sub, v, nothing)); reduced ? Poly.subst(q, v -> get(dep, v, nothing)) : q)
    n2, d2 = f(num), f(den)
    shift = Poly.lowest_shift((n2, d2), [V.index[x] for x in m.metabolites])
    n2, d2, !isempty(shift)
end
for pop in POPS
    ms = [m for m in Pop.load_population(pop) if !m.allosteric]
    full_dep = red_dep = shift_g = shift_r = 0; ex = String[]
    for m in ms
        V, sv, num, den, _ = EngineLean.bare_rate_law(m)
        _, sv2, num2, den2, _ = EngineLean.bare_rate_law(m; order = reverse(eachindex(m.steps)))
        # both laws share V (same names) only if built on one table: rebuild the second on V's names
        V2, _, _, _, _ = EngineLean.bare_rate_law(m; order = reverse(eachindex(m.steps)))
        conv(p) = Poly.subst(p, v -> Poly.mono(Poly.var!(V, V2.names[v]) => 1))
        num2, den2 = conv(num2), conv(den2)
        for reduced in (false, true)
            a, b, sa = named(m, V, sv, num, den; reduced)
            c, d, _ = named(m, V, sv, num2, den2; reduced)
            differ = a * d != c * b
            reduced ? (red_dep += differ) : (full_dep += differ)
            reduced ? (shift_r += sa) : (shift_g += sa)
            differ && !reduced && length(ex) < 3 && push!(ex, m.id)
        end
    end
    @printf("%-13s n=%3d  BFS order changes the law: Full-mode (groups) %3d, Reduced %d | substitution moves the lowest-terms shift: groups %d, groups+dependents %d   e.g. %s\n",
            pop, length(ms), full_dep, red_dep, shift_g, shift_r, join(ex, ", "))
end
