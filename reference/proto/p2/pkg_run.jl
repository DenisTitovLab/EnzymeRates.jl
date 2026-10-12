# ABOUTME: Package-side P2 measurements (needs EnzymeRates at the oracle worktree): time of
# ABOUTME: _dependent_param_exprs per plain mechanism, and the package on same-model regroupings.
using EnzymeRates, Random, Serialization, Printf, Statistics
const ER = EnzymeRates
include(joinpath(@__DIR__, "..", "common", "Population.jl"))
include(joinpath(@__DIR__, "..", "common", "Thermo.jl"))
const TESTDIR = joinpath(dirname(dirname(pathof(EnzymeRates))), "test")

function test_consts(file, names)
    ast = Meta.parseall(read(joinpath(TESTDIR, file), String))
    [ex for ex in ast.args if ex isa Expr && ex.head === :const &&
     ex.args[1] isa Expr && ex.args[1].args[1] in names]
end
for ex in test_consts("test_mechanism_enumeration.jl", (:uni_uni_rxn, :uni_bi_rxn, :bi_bi_rxn, :bi_bi_pp_rxn))
    Core.eval(Main, ex)
end
include(joinpath(TESTDIR, "mechanism_definitions_for_test_enzyme_derivation.jl"))

function expand_each(parents, rxn)
    out = Union{ER.Mechanism, ER.AllostericMechanism}[]
    for p in parents; append!(out, ER.expand_mechanisms([p], rxn)); end
    out
end
function members(pop)
    pop == "specs" && return [ER._concrete(s.mechanism) for s in MECHANISM_TEST_SPECS]
    if pop in ("init_uniuni", "init_unibi", "init_bibi", "init_bibi_pp")
        rxn = Dict("init_uniuni" => uni_uni_rxn, "init_unibi" => uni_bi_rxn,
                   "init_bibi" => bi_bi_rxn, "init_bibi_pp" => bi_bi_pp_rxn)[pop]
        return unique(ER.init_mechanisms(rxn))
    end
    @assert pop == "expand_bibi"
    init = unique(ER.init_mechanisms(bi_bi_pp_rxn))
    pick = sort!(randperm(Random.Xoshiro(20261009), length(init))[1:60])
    unique(expand_each(init[pick], bi_bi_pp_rxn))          # 700 < cap 2000, so no subsample
end

regrouped(pm, groups) = Pop.PMech(; id = pm.id, source = pm.source, substrates = pm.substrates,
    products = pm.products, inhibitors = pm.inhibitors, metabolites = pm.metabolites,
    species = pm.species, steps = pm.steps, groups = groups)

tp = Float64[]; tp_warm = Float64[]; tm = Float64[]
nmis = 0; nchk = 0
split_same_pkg = 0; split_same_mine = 0; split_agree = 0
split_fitted_diff = 0; split_key_diff = 0; split_count_diff = 0; split_err = 0
examples = String[]
for pop in ("specs", "init_uniuni", "init_unibi", "init_bibi", "init_bibi_pp", "expand_bibi")
    ms = members(pop)
    pms = Pop.load_population(pop)
    @assert length(ms) == length(pms) "$pop: $(length(ms)) vs $(length(pms))"
    for (i, m) in enumerate(ms)
        pm = pms[i]
        pm.allosteric && continue
        em = ER.compile_mechanism(m); lm = ER._concrete(em)
        s0 = ER.rate_equation_string(em)
        global nchk += 1
        ER._rate_eq_dedup_key(s0) == pm.eq_hash || (global nmis += 1; println("MISALIGNED ", pm.id); continue)
        # timing: cold naming caches (fresh rebuilds), then a warm call
        fresh = [ER._with(lm) for _ in 1:6]
        ER._dependent_param_exprs(fresh[1])
        push!(tp, minimum(@elapsed(ER._dependent_param_exprs(f)) for f in fresh[2:end]))
        push!(tp_warm, minimum(@elapsed(ER._dependent_param_exprs(lm)) for _ in 1:5))
        Thermo.solve(pm)
        push!(tm, minimum(@elapsed(Thermo.solve(pm)) for _ in 1:5))
        # same-model splits
        n0 = ER._independent_param_count(lm)
        _, ind0 = ER._dependent_param_exprs(lm)
        key0 = ER._rate_eq_dedup_key(s0)
        groups = ER.steps(lm)
        nmine = length(Thermo.solve(pm).indep)
        for (gi, g) in enumerate(groups)
            length(g) < 2 && continue
            # my solve on the same split (package group gi = PMech group gi)
            pg = [[gg for (k, gg) in enumerate(pm.groups) if k != gi]; [[j] for j in pm.groups[gi]]]
            mine_same = length(Thermo.solve(regrouped(pm, pg)).indep) == nmine
            mine_same && (global split_same_mine += 1)
            sm = try
                ER.Mechanism(ER.reaction(lm), [[gg for (k, gg) in enumerate(groups) if k != gi];
                                                [[s] for s in g]])
            catch e
                global split_err += 1; println("SPLIT-ERR ", pm.id, " ", sprint(showerror, e)[1:min(end, 200)]); continue
            end
            n1 = ER._independent_param_count(sm)
            (n1 == n0) == mine_same || (global split_count_diff += 1;
                println("COUNT-DISAGREE ", pm.id, " group ", gi, " pkg ", n0, "->", n1, " mine same=", mine_same))
            n1 == n0 || continue
            global split_same_pkg += 1
            _, ind1 = ER._dependent_param_exprs(sm)
            s1 = ER.rate_equation_string(ER.compile_mechanism(sm))
            fd = Set(ind1) != Set(ind0); kd = ER._rate_eq_dedup_key(s1) != key0
            fd && (global split_fitted_diff += 1)
            kd && (global split_key_diff += 1)
            (fd || kd) || (global split_agree += 1)
            if (fd || kd) && length(examples) < 6
                push!(examples, string("== ", pm.id, " split group ", gi, ": fitted differ=", fd, " key differ=", kd,
                    "\n-- grouped:\n", s0, "\n-- split:\n", s1))
            end
        end
    end
    println(pop, " done")
end
@printf("aligned %d of %d plain mechanisms (misaligned %d)\n", nchk - nmis, nchk, nmis)
@printf("_dependent_param_exprs cold: median %.1f us, p90 %.1f us, max %.1f us, total %.2f s\n",
        1e6median(tp), 1e6quantile(tp, 0.9), 1e6maximum(tp), sum(tp))
@printf("_dependent_param_exprs warm: median %.1f us, max %.1f us\n", 1e6median(tp_warm), 1e6maximum(tp_warm))
@printf("Thermo.solve (same process): median %.1f us, p90 %.1f us, max %.1f us, total %.3f s\n",
        1e6median(tm), 1e6quantile(tm, 0.9), 1e6maximum(tm), sum(tm))
@printf("ratio of medians package/mine: %.1f\n", median(tp) / median(tm))
println("same-model splits: package keeps count in $split_same_pkg, mine in $split_same_mine, count verdicts disagree $split_count_diff, split errors $split_err")
println("package on same-model splits: identical fitted set and dedup key $split_agree; fitted set differs $split_fitted_diff; dedup key differs $split_key_diff")
foreach(println, examples)
