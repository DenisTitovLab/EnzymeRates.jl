# ABOUTME: Checks that rule-A normalization (divide by the constant term) reproduces the package's pinned
# ABOUTME: expected_factored_num/denom strings of MECHANISM_TEST_SPECS exactly; P4 probe, stdlib only.
include(joinpath(@__DIR__, "P4.jl"))
using .Poly: Mono, mono, P, term
const PI = P{Int}
const RF = P4.RF
const FIX = "/tmp/claude-501/-home-denis-linux--julia-dev-EnzymeRates/12f42995-f5ba-4fda-9510-dd9b87f7d069/scratchpad/er_oracle/test/mechanism_definitions_for_test_enzyme_derivation.jl"
ast = Meta.parseall(read(FIX, String))
walk(f, ex) = (f(ex); ex isa Expr && foreach(a -> walk(f, a), ex.args))
strs = Dict{String, Dict{Symbol, String}}()
walk(ast) do ex
    (ex isa Expr && ex.head === :call && ex.args[1] === :_testhelper_MechanismTestSpec) || return
    kw = Dict(a.args[1] => a.args[2] for a in ex.args[2:end] if a isa Expr && a.head === :kw)
    haskey(kw, :expected_factored_denom) || return
    # string concatenations (`"a" * "b"`) are evaluated
    strs[kw[:name]] = Dict(k => eval(kw[k]) for k in (:expected_factored_num, :expected_factored_denom))
end
specs = Pop.load_population("specs")
for (name, d) in sort(collect(strs); by = first)
    m = only(filter(mm -> mm.source == "spec:" * name, specs))
    L = P4.law(m)
    # every package name (fitted or dependent, any step slot) as the law's monomial
    env = Dict{Symbol, RF}(L.vars.names[v] => RF(term(PI, mono(v => 1))) for v in values(L.conc))
    for (j, ns) in enumerate(m.step_names), (r, n) in enumerate(ns)
        env[n] = RF(term(PI, L.stepmono[j][r]))
    end
    env[:Keq] = RF(term(PI, mono(L.keq => 1)))
    sub(ex) = ex isa Symbol ? get(env, ex, ex) : ex isa Expr ? Expr(ex.head, map(sub, ex.args)...) : ex
    ev(s) = eval(sub(Meta.parse(s)))
    D = P4.cmap(L, L.den)
    c0 = RF(D[Mono()])
    okd = P4.same(RF(L.den) / c0, ev(d[:expected_factored_denom]))
    okn = P4.same(RF(L.num) / c0, ev(d[:expected_factored_num]))
    println(rpad(name, 42), "constant term a monomial: ", length(D[Mono()].t) == 1,
            "; D/c0 == expected denominator: ", okd, "; N/c0 == expected numerator: ", okn)
end
