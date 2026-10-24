# ABOUTME: Exports the R5.6/R5.7 inputs from the package as census text: the 42 specs, the six
# ABOUTME: size-gate laws and the largest ter-bi seeds (small), or a stratified ter-ter set (terter).
using EnzymeRates, Random, Printf, SHA, LinearAlgebra, Test, Serialization
const ER = EnzymeRates
const CEN = "/tmp/claude-501/-home-denis-linux--julia-dev-EnzymeRates/12f42995-f5ba-4fda-9510-dd9b87f7d069/scratchpad/phase3/census"
include(joinpath(CEN, "census_common.jl"))
include(joinpath(CEN, "neutral_steps.jl"))
const WT = dirname(dirname(pathof(EnzymeRates)))
println("EnzymeRates ", pkgversion(EnzymeRates), " at ",
        readchomp(`git -C $WT rev-parse --short HEAD`), "; Julia ", VERSION)

cxA(m) = Int(ER._eq_complexity(census_concrete(m)))
cxI(m) = (c = census_concrete(m); c isa ER.AllostericMechanism ?
          Int(ER._eq_complexity(ER._inactive_groups(c))) : 0)

"""Seconds per call of the package's `rate_equation`, the gate helper's way
(test_rate_eq_derivation.jl:209-230: minimum over 5 batches of 10,000 calls)."""
function pkg_time(em, rng)
    lu() = 10.0^(2 * rand(rng) - 1)
    concs = (; (x => lu() for x in ER.metabolites(em))...)
    params = (; (n => lu() for n in ER.fitted_params(em))..., Keq = 10.0, E_total = 1.0)
    ER.rate_equation(em, concs, params)
    allocs = @allocated ER.rate_equation(em, concs, params)
    best, acc = Inf, 0.0
    for _ in 1:5
        acc = 0.0
        t = @elapsed for _ in 1:10_000
            acc += ER.rate_equation(em, concs, params)
        end
        best = min(best, t / 10_000)
    end
    best, allocs
end

function row!(io, set, label, m; timed = false, rng = Xoshiro(1))
    t = census_text(m)
    s, p, _ = census_roles(m)
    pkgD, ns, al = "", "", ""
    if timed
        em = census_compiled(m)
        pkgD = try
            string(length(ER._raw_symbolic_rate_polys(typeof(em))[2]))
        catch
            ""
        end
        tt, a = pkg_time(em, rng)
        ns, al = @sprintf("%.2f", 1e9 * tt), string(a)
    end
    write_row(io, set, label, mech_id(t), join(s, ","), join(p, ","), string(cxA(m)),
              string(cxI(m)), pkgD, ns, al, t)
end

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
    [(s.name, s.mechanism) for s in Core.eval(mod, :MECHANISM_TEST_SPECS)]
end

"""The six laws of the size-scaled gate (test_rate_eq_derivation.jl:1253-1302)."""
function gate_laws()
    path = joinpath(WT, "test", "test_rate_eq_derivation.jl")
    mod = Module(:GateLaws)
    Core.eval(mod, :(using EnzymeRates))
    names = (:_testhelper_largest_enumerated_mechanisms, :_testhelper_large_ter_bi_mechanisms)
    for ex in Meta.parseall(read(path, String); filename = path).args
        ex isa Expr && ex.head === :const && ex.args[1].args[1] in names &&
            Core.eval(mod, ex)
    end
    [("bibi396_$i", m) for (i, m) in enumerate(Core.eval(mod, names[1]))] ∪
        [("terbi_$i", m) for (i, m) in enumerate(Core.eval(mod, names[2]))]
end

const TERBI = """
@enzyme_reaction begin
    substrates: A[C], B[N], D[X]
    products: P[CN], Q[X]
end"""
const TERTER = """
@enzyme_reaction begin
    substrates: A[C], B[N], D[X]
    products: P[C], Q[N], R[X]
end"""

t0 = time()
mode = ARGS[1]
out = joinpath(@__DIR__, "data", "x5_$mode.tsv")
open(out, "w") do io
    if mode == "small"
        rng = Xoshiro(20261014)
        for (n, m) in load_specs()
            row!(io, "spec", n, m; timed = true, rng)
        end
        for (n, m) in gate_laws()
            row!(io, "gate", n, m; timed = true, rng)
        end
        rxn = Core.eval(Main, Meta.parse(TERBI))
        ms = ER.init_mechanisms(rxn)
        cx = cxA.(ms)
        top = sortperm(cx; rev = true, alg = MergeSort)
        top = [i for i in top if cx[i] <= 5000][1:30]
        @printf("ter-bi init %d; V×τ max %d, top30 within 5000 from %d to %d\n", length(ms),
                maximum(cx), cx[top[1]], cx[top[end]])
        for i in top
            row!(io, "terbi_top", string(i), ms[i]; timed = true, rng)
        end
    else
        rxn = Core.eval(Main, Meta.parse(TERTER))
        t = @elapsed ms = ER.init_mechanisms(rxn)
        cx = cxA.(ms)
        @printf("ter-ter init %d in %.0f s; V×τ computed by %.0f s; maxrss %.0f MB\n",
                length(ms), t, time() - t0, Sys.maxrss() / 2^20)
        fitpick = Any[]
        bins = [(1, 50, 150), (51, 337, 150), (338, 1000, 150), (1001, 2000, 200),
                (2001, 5000, 300)]
        rng = Xoshiro(20261014)
        for (lo, hi, k) in bins
            idx = findall(x -> lo <= x <= hi, cx)
            pick = length(idx) <= k ? idx : idx[sort(randperm(rng, length(idx))[1:k])]
            @printf("bin [%d, %d]: %d mechanisms, %d picked\n", lo, hi, length(idx), length(pick))
            foreach(i -> row!(io, "terter_init", string(i), ms[i]), pick)
            # two mechanisms per bin for fittime.jl, with their census ids
            for i in pick[1:2]
                push!(fitpick, ("terter_init", string(i), mech_id(census_text(ms[i])), ms[i]))
            end
        end
        top = [i for i in sortperm(cx; rev = true, alg = MergeSort) if cx[i] <= 5000][1:40]
        foreach(i -> row!(io, "terter_top", string(i), ms[i]), top)
        @printf("over 5000: %d; top40 within 5000 from %d to %d\n", count(>(5000), cx),
                cx[top[1]], cx[top[end]])
        small = findall(<=(337), cx)
        mid = findall(x -> 1000 < x <= 2500, cx)
        parents = [ms[small[randperm(rng, length(small))[1:15]]];
                   ms[mid[randperm(rng, length(mid))[1:10]]]]
        ms = nothing
        GC.gc()
        kids = Any[]
        for p in parents
            append!(kids, ER.expand_mechanisms([p], rxn))
        end
        unique!(kids)
        ok = [c for c in kids if cxA(c) + cxI(c) <= 5000]
        pick = length(ok) <= 400 ? ok : ok[sort(randperm(rng, length(ok))[1:400])]
        @printf("round 1 from 25 parents: %d children, %d allosteric, %d with V×τ sum <= 5000, %d picked\n",
                length(kids), count(c -> c isa ER.AllostericMechanism, kids), length(ok),
                length(pick))
        foreach(c -> row!(io, "terter_r1", "", c), pick)
        for c in filter(c -> c isa ER.AllostericMechanism, pick)[1:2]
            push!(fitpick, ("terter_r1", "", mech_id(census_text(c)), c))
        end
        serialize(joinpath(@__DIR__, "data", "fitpick.jls"), fitpick)
    end
end
@printf("wrote %s in %.0f s; maxrss %.0f MB\n", out, time() - t0, Sys.maxrss() / 2^20)
