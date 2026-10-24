# ABOUTME: Exports unit-rule census populations from the package (census env) as data/<pop>.tsv:
# ABOUTME: id, substrates, products, neutral text, today's printed law (specs). Args: <pop>...
using EnzymeRates, Random, Printf, LinearAlgebra, Test
const ER = EnzymeRates
const CEN = "/tmp/claude-501/-home-denis-linux--julia-dev-EnzymeRates/12f42995-f5ba-4fda-9510-dd9b87f7d069/scratchpad/phase3/census"
include(joinpath(CEN, "census_common.jl"))
include(joinpath(CEN, "neutral_steps.jl"))
const WT = dirname(dirname(pathof(EnzymeRates)))
println("EnzymeRates ", pkgversion(EnzymeRates), " from ", WT, " at ",
        readchomp(`git -C $WT rev-parse --short HEAD`))

# The census reactions (phase3/census/corpus.jl:22-52) and ter_ter_rxn
# (test/test_mechanism_enumeration.jl:252-255).
const RXN = Dict(
    "bibi" => """
@enzyme_reaction begin
    substrates: A[C], B[N]
    products: P[C], Q[N]
    competitive_inhibitors: A, I
end""",
    "pingpong" => """
@enzyme_reaction begin
    substrates: A[CX], B[N]
    products: P[C], Q[NX]
    allosteric_regulators: R
    oligomeric_state: 2
    shared_catalytic_site: (A, P)
end""",
    "terter" => """
@enzyme_reaction begin
    substrates: A[C], B[N], D[X]
    products: P[C], Q[N], R[X]
end""")

"""The 42 specs: the worktree's spec file evaluated in a fresh module with corpus.jl's
stand-ins for the test-only source macros (phase3/census/corpus.jl:194-233)."""
function load_specs()
    path = joinpath(WT, "test", "mechanism_definitions_for_test_enzyme_derivation.jl")
    mod = Module(:CensusSpecs)
    Core.eval(mod, :(using EnzymeRates, LinearAlgebra, Random, Test))
    Core.eval(mod, quote
        macro _testhelper_enzyme_mechanism_src(block)
            esc(Expr(:tuple, Expr(:macrocall, Symbol("@enzyme_mechanism"), __source__,
                                  block), nothing))
        end
        macro _testhelper_allosteric_mechanism_src(block)
            esc(Expr(:tuple, Expr(:macrocall, Symbol("@allosteric_mechanism"), __source__,
                                  block), nothing, nothing))
        end
        macro _testhelper_with_raw(call)
            esc(Expr(:tuple, call, nothing))
        end
    end)
    srcmacro(ex) = ex isa Expr && (ex.head === :macro && ex.args[1].args[1] in
        (:_testhelper_enzyme_mechanism_src, :_testhelper_allosteric_mechanism_src,
         :_testhelper_with_raw) ||
        ex.head === :macrocall && ex.args[1] isa GlobalRef &&
        ex.args[1].name === Symbol("@doc") && srcmacro(ex.args[end]))
    skipped = 0
    for ex in Meta.parseall(read(path, String); filename = path).args
        ex isa LineNumberNode && continue
        srcmacro(ex) && continue
        try
            Core.eval(mod, ex)
        catch
            skipped += 1
        end
    end
    println("spec file: $skipped statements skipped")
    [("spec:" * s.name, s.mechanism) for s in Core.eval(mod, :MECHANISM_TEST_SPECS)]
end

function write_pop(name, items)
    out = joinpath(@__DIR__, "data", name * ".tsv")
    seen = Set{String}()
    open(out, "w") do io
        for (id, m) in items
            t = census_text(m)
            key = id === nothing ? mech_id(t) : id
            key in seen && continue
            push!(seen, key)
            s, p, _ = census_roles(m)
            # today's printed law, for the specs only (each costs a derivation)
            today = name == "spec" ? ER.rate_equation_string(m) : ""
            write_row(io, key, join(s, ","), join(p, ","), t, today)
        end
    end
    @printf("wrote %s: %d rows; maxrss %.0f MB\n", out, length(seen), Sys.maxrss() / 2^20)
end

for pop in ARGS
t0 = time()
if pop == "spec"
    write_pop(pop, load_specs())
elseif pop == "terter"
    rxn = Core.eval(Main, Meta.parse(RXN[pop]))
    t = @elapsed ms = ER.init_mechanisms(rxn)
    @printf("init_mechanisms(ter_ter_rxn): %d in %.1f s; maxrss %.0f MB\n", length(ms), t,
            Sys.maxrss() / 2^20)
    pick = sort(randperm(Xoshiro(20261014), length(ms))[1:2000])
    write_pop(pop, [(nothing, ms[i]) for i in pick])
else
    rxn = Core.eval(Main, Meta.parse(RXN[pop]))
    init = ER.init_mechanisms(rxn)
    items = Any[(nothing, m) for m in init]
    errs = 0
    for p in init
        kids = try
            ER.expand_mechanisms([p], rxn)
        catch e
            errs += 1
            continue
        end
        append!(items, [(nothing, c) for c in kids])
    end
    @printf("%s: init %d, child rows %d, parent errors %d\n", pop, length(init),
            length(items) - length(init), errs)
    write_pop(pop, items)
end
@printf("%s done in %.0f s\n", pop, time() - t0)
end
