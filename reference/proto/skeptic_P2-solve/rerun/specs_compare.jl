# ABOUTME: Side-by-side dependent choices (mine vs package) on the plain specs, in the package's names.
# ABOUTME: Throwaway P2 probe script; stdlib only.
include(joinpath("/tmp/claude-501/-home-denis-linux--julia-dev-EnzymeRates/12f42995-f5ba-4fda-9510-dd9b87f7d069/scratchpad/proto/common", "Thermo.jl"))
include(joinpath(@__DIR__, "census_helpers.jl"))
const T = Thermo
pkgname(m, (g, r)) = m.step_names[m.groups[g][1]][r]
same = 0
for p in ("specs", "handwritten"), m in Pop.load_population(p)
    m.allosteric && continue
    sol = T.solve(m)
    mine = sort([string(pkgname(m, my_slot(sol, d.col)), d.label === :haldane ? "(H)" : "(W)", d.fold ? "=fold" : "")
                 for d in sol.deps if d.label !== :kinetic_group])
    pk = sort([string(n, any(q -> q.first === :Keq, mono) ? "(H)" : "(W)") for (n, mono) in m.dependent_monomials])
    append!(pk, [string(a, "(W)=fold") for (a, _) in m.pass2_rename]); sort!(pk)
    ms = Set(replace.(mine, "=fold" => "")); ps = Set(replace.(pk, "=fold" => ""))
    ms == ps && (global same += 1)
    println(rpad(m.source[1:min(end, 58)], 60), ms == ps ? "SAME " : "DIFF ", "mine: ", join(mine, ", "),
            ms == ps ? "" : "   | pkg: " * join(pk, ", "))
end
println("same: ", same)
