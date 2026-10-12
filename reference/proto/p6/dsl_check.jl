# ABOUTME: Lowers every MECHANISM_TEST_SPECS fixture, as written in the test source, through the P6
# ABOUTME: DSL prototype and checks it equals the GM of the package's exported spec mechanism.
include(joinpath(@__DIR__, "harness.jl"))
include(joinpath(@__DIR__, "DslLower.jl"))
using Printf

const FIXTURES = joinpath(@__DIR__, "..", "..", "er_oracle", "test",
                          "mechanism_definitions_for_test_enzyme_derivation.jl")
const MACROS = Symbol.(("@_testhelper_enzyme_mechanism_src", "@_testhelper_allosteric_mechanism_src",
                        "@enzyme_mechanism", "@allosteric_mechanism"))
bodies(x, acc = Expr[]) = (x isa Expr && (x.head === :macrocall && x.args[1] in MACROS ?
    push!(acc, x.args[end]) : foreach(a -> bodies(a, acc), x.args)); acc)
specnames(x, acc = String[]) = (x isa Expr && (x.head === :call && x.args[1] === :_testhelper_MechanismTestSpec ?
    append!(acc, [a.args[2] for a in x.args if a isa Expr && a.head === :kw && a.args[1] === :name]) :
    foreach(a -> specnames(a, acc), x.args)); acc)

specs = Dict(replace(m.source, "spec:" => "") => m for m in Pop.load_population("specs"))
ast = Meta.parseall(read(FIXTURES, String))
nequal = 0; npaired = 0; unpaired = String[]; differ = String[]; lowered = Dict{String, GM.GMech}()
for st in ast.args
    names, bs = specnames(st), bodies(st)
    isempty(names) && continue
    if length(names) != 1 || length(bs) != 1
        push!(unpaired, "$(names) with $(length(bs)) macro bodies"); continue
    end
    global npaired += 1
    raw = DslLower.lower(only(bs))
    g = _testhelper_gm(raw)
    lowered[only(names)] = g
    ref = _testhelper_gm(_testhelper_toraw(specs[only(names)]))
    if g == ref && hash(g) == hash(ref)
        global nequal += 1
    else
        push!(differ, only(names))
    end
end
@printf("spec fixtures paired with one macro body: %d of %d specs; equal to the exported spec: %d\n",
        npaired, length(specs), nequal)
println("unpaired: ", unpaired)
println("different: ", differ)
refs = Dict(k => _testhelper_gm(_testhelper_toraw(m)) for (k, m) in specs)
cross = [(a, b) for (a, g) in lowered, (b, r) in refs if a != b && g == r]
println("lowered fixtures equal to a different spec: ", cross)
ex = lowered["Segel Ping Pong Bi Bi"]
println("Segel Ping Pong Bi Bi, canonical edges (group, stored a < b):")
for (k, e) in enumerate(ex.edges)
    println("  g", ex.group[k], "  ", GM.step_text(ex, e))
end
