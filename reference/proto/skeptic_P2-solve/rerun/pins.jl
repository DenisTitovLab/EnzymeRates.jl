# ABOUTME: My solve on the three fixtures of the package testset "rate_equation_string marks only the
# ABOUTME: ties folded into v": a binding-K fold, an isomerization-K tie, and a constant pinned to Keq.
include(joinpath("/tmp/claude-501/-home-denis-linux--julia-dev-EnzymeRates/12f42995-f5ba-4fda-9510-dd9b87f7d069/scratchpad/proto/common", "Thermo.jl"))
const T = Thermo
M(n, r) = Pop.PMet(n, r)
sp(conf, b...) = Pop.PSpecies(collect(b), conf, Symbol[], Symbol[], Symbol(conf, join(string.(getfield.(collect(b), :name)))))
mk(id, subs, prods, inh, species, steps, groups) = Pop.PMech(; id, source = id, substrates = subs, products = prods,
    inhibitors = inh, metabolites = [subs; prods; inh], species, steps, groups)
A = M(:A, :S); B = M(:B, :S); P = M(:P, :P); S = M(:S, :S); I = M(:I, :Inh); X = Pop.PMet[]
folded = mk("folded", [:A, :B], [:P], Symbol[],
    [sp(:E), sp(:E, A), sp(:E, B), sp(:E, A, B), sp(:E, P)],
    [Pop.PStep(1, 2, [A], X, true), Pop.PStep(3, 4, [A], X, true), Pop.PStep(1, 3, [B], X, true),
     Pop.PStep(2, 4, [B], X, true), Pop.PStep(4, 5, X, X, false), Pop.PStep(1, 5, [P], X, true)],
    [[1, 2], [3], [4], [5], [6]])
iso_tie = mk("iso_tie", [:S], [:P], Symbol[],
    [sp(:E), sp(:Estar), sp(:E, S), sp(:Estar, S), sp(:E, P)],
    [Pop.PStep(1, 3, [S], X, true), Pop.PStep(2, 4, [S], X, true), Pop.PStep(1, 2, X, X, true),
     Pop.PStep(3, 4, X, X, true), Pop.PStep(3, 5, X, X, false), Pop.PStep(1, 5, [P], X, true)],
    [[1, 2], [3], [4], [5], [6]])
keq_tie = mk("keq_tie", [:S], [:P], [:I],
    [sp(:E), sp(:Estar), sp(:E, I)],
    [Pop.PStep(1, 2, [S], [P], true), Pop.PStep(1, 3, [I], X, false), Pop.PStep(2, 3, [I], X, false)],
    [[1], [2, 3]])
for m in (folded, iso_tie, keq_tie)
    sol = T.solve(m)
    println("== ", m.id, ": independents ", T.independent_names(sol))
    for d in sol.deps
        println("   ", T.colname(sol, d.col), " = Keq^", d.keq, " ", join([string(T.colname(sol, f), "^", e) for (f, e) in d.terms], " "),
                "   [", d.label, d.fold ? ", fold" : "", isempty(d.terms) ? ", pinned" : "", "]")
    end
end
