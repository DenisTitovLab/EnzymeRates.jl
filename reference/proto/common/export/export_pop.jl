# ABOUTME: Builds one mechanism population with EnzymeRates and serializes it as Pop.PMech records.
# ABOUTME: Usage: julia --project=<er_oracle> export_pop.jl <population> [chunk nchunks]

using EnzymeRates, Random, Serialization, Printf, Statistics
include(joinpath(@__DIR__, "..", "Population.jl"))
include(joinpath(@__DIR__, "export_lib.jl"))

const POP = ARGS[1]
const CHUNK = length(ARGS) >= 2 ? parse(Int, ARGS[2]) : 1
const NCHUNKS = length(ARGS) >= 3 ? parse(Int, ARGS[3]) : 1
const TESTDIR = joinpath(dirname(dirname(pathof(EnzymeRates))), "test")
const T0 = time()

"""Top-level `const` statements of a test file whose names are in `names`, in file order."""
function test_consts(file, names)
    ast = Meta.parseall(read(joinpath(TESTDIR, file), String))
    [ex for ex in ast.args if ex isa Expr && ex.head === :const &&
     ex.args[1] isa Expr && ex.args[1].args[1] in names]
end

# The reactions as test/test_mechanism_enumeration.jl defines them.
for ex in test_consts("test_mechanism_enumeration.jl",
                      (:uni_uni_rxn, :uni_bi_rxn, :bi_bi_rxn, :bi_bi_pp_rxn, :uni_uni_allo_reg))
    Core.eval(Main, ex)
end
# The bi-bi reaction at catalytic multiplicity 2 (the `bi_bi_allo_rxn` of the enumeration
# tests and the reaction of the `_expand_to_allosteric` aggregate pin).
const bi_bi_allo_rxn = @enzyme_reaction begin
    substrates: A[C], B[N]
    products: P[C], Q[N]
    oligomeric_state: 2
end

POP == "specs" &&
    include(joinpath(TESTDIR, "mechanism_definitions_for_test_enzyme_derivation.jl"))

"""Deterministic subsample of `n` of `xs`, kept in emission order."""
subsample(xs, n, seed) = length(xs) <= n ? xs :
    xs[sort!(randperm(Random.Xoshiro(seed), length(xs))[1:n])]

const NOTES = String[]
note(s) = (push!(NOTES, s); println(s))

"""`expand_mechanisms` run parent by parent (the same children in the same order as one
call on all parents), recording a parent whose expansion throws instead of aborting."""
function expand_each(parents, rxn, label)
    out = Union{ER.Mechanism, ER.AllostericMechanism}[]
    for (i, p) in enumerate(parents)
        try
            append!(out, ER.expand_mechanisms([p], rxn))
        catch e
            note("$label: expand_mechanisms failed on parent $i: $(_errtext(e))")
        end
    end
    out
end

"""Handwritten mechanisms of the BigFloat mass-action testsets and the collapse tests."""
function handwritten_members()
    sandbox = Module(:Handwritten)
    Core.eval(sandbox, :(using EnzymeRates, Random))
    out = Tuple{Any, String}[]
    file = "test_rate_eq_derivation.jl"
    for ex in test_consts(file, (:_testhelper_fused_and_plain_binding,
                                 :_testhelper_fused_cases, :_testhelper_tc_cases))
        Core.eval(sandbox, ex)
    end
    fused = Base.invokelatest(getglobal, sandbox, :_testhelper_fused_cases)
    for (i, em) in enumerate(fused)
        push!(out, (em, "handwritten:$file:_testhelper_fused_cases[$i]"))
    end
    tc = Base.invokelatest(getglobal, sandbox, :_testhelper_tc_cases)
    for k in keys(tc)
        push!(out, (tc[k], "handwritten:$file:_testhelper_tc_cases.$k"))
    end
    # Mechanisms assigned inside the named testsets.
    wanted = ("the free-enzyme weight is free E's weight in the reduced denominator",
              "a fused release keeps its catalytic rate constant fitted",
              "reversing written steps changes nothing",
              "a constant thermodynamics pins to 1 is the dependent constant 1.0")
    ismech(x) = x isa Expr && x.head === :macrocall &&
                x.args[1] in (Symbol("@enzyme_mechanism"), Symbol("@allosteric_mechanism"))
    function walk(ex, title)
        ex isa Expr || return
        if ex.head === :macrocall && ex.args[1] === Symbol("@testset")
            t = findfirst(a -> a isa String, ex.args)
            title = t === nothing ? title : ex.args[t]
        end
        if title in wanted && ex.head === :(=) && ex.args[1] isa Symbol && ismech(ex.args[2])
            push!(out, (Core.eval(sandbox, ex.args[2]),
                        "handwritten:$file:testset \"$title\": $(ex.args[1])"))
            return
        end
        foreach(a -> walk(a, title), ex.args)
    end
    walk(Meta.parseall(read(joinpath(TESTDIR, file), String)), "")
    # test_allosteric_collapse.jl: its helper definitions, called with the testsets' tags.
    cfile = "test_allosteric_collapse.jl"
    mod = only(ex for ex in Meta.parseall(read(joinpath(TESTDIR, cfile), String)).args
               if ex isa Expr && ex.head === :module)
    for ex in mod.args[3].args
        ex isa Expr || continue
        ex.head === :using && continue
        ex.head === :macrocall && ex.args[1] === Symbol("@testset") && continue
        Core.eval(sandbox, ex)
    end
    uni = Base.invokelatest(getglobal, sandbox, :_testhelper_uni)
    ro = Base.invokelatest(getglobal, sandbox, :_testhelper_ro_bibi)
    st3 = fill(:EqualAI, 9); st3[3] = :NonequalAI
    st23 = fill(:EqualAI, 9); st23[2] = :NonequalAI; st23[3] = :NonequalAI
    calls = [
        ("uni([:EqualAI,:EqualAI,:EqualAI])", () -> uni([:EqualAI, :EqualAI, :EqualAI])),
        ("uni([:NonequalAI,:EqualAI,:EqualAI])", () -> uni([:NonequalAI, :EqualAI, :EqualAI])),
        ("uni([:NonequalAI,:EqualAI,:NonequalAI])",
         () -> uni([:NonequalAI, :EqualAI, :NonequalAI])),
        ("uni([:NonequalAI,:NonequalAI,:EqualAI])",
         () -> uni([:NonequalAI, :NonequalAI, :EqualAI])),
        ("ro_bibi(st[3]=:NonequalAI)", () -> ro(st3)),
        ("uni([:NonequalAI,:EqualAI,:EqualAI]; re=(false,false,false))",
         () -> uni([:NonequalAI, :EqualAI, :EqualAI]; re = (false, false, false))),
        ("uni([:NonequalAI,:NonequalAI,:EqualAI]; re=(false,false,false))",
         () -> uni([:NonequalAI, :NonequalAI, :EqualAI]; re = (false, false, false))),
        ("uni([:NonequalAI,:EqualAI,:NonequalAI]; re=(true,false,false))",
         () -> uni([:NonequalAI, :EqualAI, :NonequalAI]; re = (true, false, false))),
        ("ro_bibi(st[2]=st[3]=:NonequalAI)", () -> ro(st23)),
        ("uni([:NonequalAI,:OnlyA,:EqualAI])", () -> uni([:NonequalAI, :OnlyA, :EqualAI])),
    ]
    for (label, f) in calls
        push!(out, (Base.invokelatest(f), "handwritten:$cfile:$label"))
    end
    # Lift compiled values to concrete mechanisms; keep the first of structural duplicates.
    lifted = Tuple{Any, String}[]
    for (x, src) in out
        m = x isa EnzymeRates.AbstractEnzymeMechanism ? ER._concrete(x) : x
        dup = findfirst(((y, _),) -> typeof(y) == typeof(m) && y == m, lifted)
        if dup === nothing
            push!(lifted, (m, src))
        else
            note("skipped (duplicate of $(lifted[dup][2])): $src")
        end
    end
    lifted
end

function members(pop)
    if pop == "specs"
        return [(ER._concrete(s.mechanism), "spec:" * s.name) for s in MECHANISM_TEST_SPECS]
    elseif pop in ("init_uniuni", "init_unibi", "init_bibi", "init_bibi_pp")
        rxn = Dict("init_uniuni" => uni_uni_rxn, "init_unibi" => uni_bi_rxn,
                   "init_bibi" => bi_bi_rxn, "init_bibi_pp" => bi_bi_pp_rxn)[pop]
        raw = ER.init_mechanisms(rxn)
        ms = unique(raw)
        note("$pop: init_mechanisms returned $(length(raw)), $(length(ms)) distinct")
        return [(m, "$pop[$i]") for (i, m) in enumerate(ms)]
    elseif pop == "expand_bibi"
        init = unique(ER.init_mechanisms(bi_bi_pp_rxn))
        pick = sort!(randperm(Random.Xoshiro(20261009), length(init))[1:60])
        children = unique(expand_each(init[pick], bi_bi_pp_rxn, "expand_bibi"))
        note("expand_bibi: parents init_bibi_pp$(pick); $(length(children)) distinct children")
        kept = subsample(children, 2000, 20261010)
        note("expand_bibi: kept $(length(kept)) (cap 2000)")
        return [(m, "expand_bibi[$i]") for (i, m) in enumerate(kept)]
    elseif pop == "allo_uniuni"
        seeds = unique(ER.init_mechanisms(uni_uni_allo_reg))
        l1 = unique(filter(c -> c isa ER.AllostericMechanism,
                           expand_each(seeds, uni_uni_allo_reg, "allo_uniuni l1")))
        l2 = unique(filter(c -> c isa ER.AllostericMechanism && !(c in l1),
                           expand_each(l1, uni_uni_allo_reg, "allo_uniuni l2")))
        note("allo_uniuni: $(length(seeds)) seeds, $(length(l1)) allosteric children, " *
             "$(length(l2)) new allosteric grandchildren")
        l2 = subsample(l2, 1000, 20261011)
        return [[(m, "allo_uniuni:child[$i]") for (i, m) in enumerate(l1)];
                [(m, "allo_uniuni:grandchild[$i]") for (i, m) in enumerate(l2)]]
    elseif pop == "allo_bibi"
        seeds = unique(ER.init_mechanisms(bi_bi_allo_rxn))
        l1 = [c for m in seeds for c in ER._expand_to_allosteric(m, bi_bi_allo_rxn)]
        u1 = unique(l1)
        note("allo_bibi: $(length(seeds)) seeds, $(length(l1)) _expand_to_allosteric " *
             "children, $(length(u1)) distinct")
        l2 = unique(filter(c -> !(c in u1), expand_each(u1, bi_bi_allo_rxn, "allo_bibi l2")))
        note("allo_bibi: $(length(l2)) distinct expand_mechanisms children not in level 1 " *
             "($(count(c -> c isa ER.AllostericMechanism, l2)) allosteric)")
        l2 = subsample(l2, 500, 20261012)
        return [[(m, "allo_bibi:to_allosteric[$i]") for (i, m) in enumerate(u1)];
                [(m, "allo_bibi:expand[$i]") for (i, m) in enumerate(l2)]]
    elseif pop == "handwritten"
        return handwritten_members()
    end
    error("unknown population $pop")
end

const MEMBERS = members(POP)
const RANGE = let n = length(MEMBERS), per = cld(n, NCHUNKS)
    ((CHUNK - 1) * per + 1):min(CHUNK * per, n)
end
const POPSEED = Dict("specs" => 1, "init_uniuni" => 2, "init_unibi" => 3, "init_bibi" => 4,
                     "init_bibi_pp" => 5, "expand_bibi" => 6, "allo_uniuni" => 7,
                     "allo_bibi" => 8, "handwritten" => 9)[POP]
note("$POP chunk $CHUNK/$NCHUNKS: members $(RANGE) of $(length(MEMBERS)); " *
     "enumeration took $(round(time() - T0; digits = 1)) s")

# kcat on every spec and handwritten mechanism, on 100 of each enumerated population.
const KCAT_IDX = POP in ("specs", "handwritten") ? Set(1:length(MEMBERS)) :
    Set(subsample(collect(1:length(MEMBERS)), 100, 1000 + POPSEED))

out = Pop.PMech[]
const T1 = time()
for i in RANGE
    m, src = MEMBERS[i]
    t = time()
    pm = try
        to_pmech(m, "$POP/$i", src; seed = 100_000 * POPSEED + i, do_kcat = i in KCAT_IDX)
    catch e
        note("FATAL $POP/$i ($src): $(_errtext(e))")
        rethrow()
    end
    push!(out, pm)
    dt = time() - t
    msg = isempty(pm.stage_errors) ? "" :
          "  errors: " * join(("$(k): $(first(v, 160))" for (k, v) in pm.stage_errors), " | ")
    @printf("%5d %-10s ok=%d %6.2fs  %s%s\n", i, pm.allosteric ? "allo" : "plain", pm.ok,
            dt, first(src, 70), msg)
    flush(stdout)
end
const ELAPSED = time() - T1
part = NCHUNKS == 1 ? POP : "$POP.part$CHUNK"
path = Pop.save_population(part, out)

fam = Dict{Symbol, Int}()
for pm in out, f in pm.family; fam[f] = get(fam, f, 0) + 1; end
open(joinpath(Pop.DATA_DIR, "$part.summary.txt"), "w") do io
    println(io, "population $POP chunk $CHUNK/$NCHUNKS: $(length(out)) mechanisms ",
            "(members $(RANGE) of $(length(MEMBERS)))")
    println(io, "export time $(round(ELAPSED; digits = 1)) s; ",
            "Julia $(VERSION); EnzymeRates at $(pkgdir(EnzymeRates))")
    println(io, "plain $(count(!(pm.allosteric) for pm in out)), ",
            "allosteric $(count(pm.allosteric for pm in out)), ",
            "not ok $(count(!(pm.ok) for pm in out))")
    for n in NOTES; println(io, "note: ", n); end
    println(io, "families: ", join(("$k=$v" for (k, v) in sort!(collect(fam); by = first)), ", "))
    for pm in out, (k, v) in pm.stage_errors
        println(io, "error ", pm.id, " (", pm.source, ") ", k, ": ", replace(v, '\n' => ' '))
    end
end
@printf("wrote %s (%d mechanisms) in %.1f s; maxrss %.0f MB\n", path, length(out), ELAPSED,
        Sys.maxrss() / 2^20)
