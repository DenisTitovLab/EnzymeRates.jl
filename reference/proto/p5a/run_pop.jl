# ABOUTME: Light :OnlyA limit study over whole allosteric populations (P5a): cone feasibility,
# ABOUTME: uniqueness, sequestration and agreement with deletion for conventions (i) and (ii).
include(joinpath(@__DIR__, "limits_lib.jl"))

const NAME = ARGS[1]
pop = filter(m -> m.allosteric && :OnlyA in m.group_tags, Pop.load_population(NAME))
length(ARGS) > 1 && (pop = pop[parse(Int, ARGS[2]):min(end, parse(Int, ARGS[3]))])

"""Cone limit study of one mechanism (no ε sweep): up to `ncone` directions, deletion, the
RREF path of the package-like order and its reversal, and one oracle check."""
function light(m, conv; seed = 11, ncone = 16)
    rng = Random.Xoshiro(seed)
    untie = conv === :iu
    conv = untie ? :i : conv
    J = MWC.joint_solve(m; order = :pkg, untie)
    hat = MWC.random_hat(J, rng)
    keq, L = logu(rng), logu(rng)
    rA, rI = draw_reg(m, rng)
    concs = Dict(x => logu(rng) for x in m.metabolites)
    aint = MWC.intended(J, conv); a = MWC.realized(J, aint)
    partner, wrong = MWC.path_report(J, aint, a)
    r = Dict{Symbol, Any}(:id => m.id, :conv => untie ? :iu : conv, :npartner => length(partner),
                          :nwrong => length(wrong))
    v1 = v_path(J, hat, a, big(1), keq, L, rA, rI, concs)
    S = sym_setup(J, hat, keq, L, rA, rI, concs)
    srref = sym_eval(S, a)
    r[:v_rref] = srref.v
    θ = MWC.all_values(J, hat, keq)
    key(Jx, i) = (Jx.cols[i].st, Jx.cols[i].j, Jx.cols[i].r)
    θk = Dict(key(J, i) => v for (i, v) in θ)
    Jr = MWC.joint_solve(m; order = :rev, untie)
    ar = MWC.realized(Jr, MWC.intended(Jr, conv))
    r[:v_rref_rev] = sym_limit(Jr, ar, Dict(i => θk[key(Jr, i)] for i in Jr.indep), keq, L,
                               rA, rI, concs).v
    r[:order_spread] = reldiff(r[:v_rref_rev], srref.v, 0)
    r[:v_del] = v_deletion(J, hat, keq, L, rA, rI, concs)
    dirs, d = MWC.cone_directions(J, conv)
    r[:cone_dim], r[:cone_n] = d, length(dirs)
    if !isempty(dirs)
        pick = unique([1; round.(Int, range(1, length(dirs); length = min(ncone, length(dirs))))])
        res = [sym_eval(S, dirs[k]) for k in pick]
        r[:cone_v] = res[1].v
        r[:cone_spread] = maximum(reldiff(s.v, res[1].v, 0) for s in res)
        r[:cone_seq] = any(s -> s.ordI[2] !== nothing && s.ordI[3] !== nothing &&
                                s.ordI[2] < s.ordI[3], res)
        r[:cone_vs_del] = maximum(reldiff(s.v, r[:v_del], 0) for s in res)
        r[:zero_over_zero] = any(s -> s.ordI[2] !== nothing && s.ordI[2][1] > 0, res)
        e = setprecision(() -> BigFloat(10)^(-max(6, 30 ÷ max(1, Int(maximum(abs, dirs[1]))))),
                         BigFloat, MWC.PREC)
        r[:oracle_err] = reldiff(v_path(J, hat, dirs[1], e, keq, L, rA, rI, concs), res[1].v, v1)
    end
    r
end

recs = Dict{Symbol, Any}[]
t0 = time()
const CONVS = (:i, :ii, :iu)
for (k, m) in enumerate(pop), conv in CONVS
    r = try
        light(m, conv)
    catch e
        Dict{Symbol, Any}(:id => m.id, :conv => conv, :error => first(sprint(showerror, e), 200))
    end
    r[:tags] = m.group_tags
    push!(recs, r)
    k % 50 == 0 && conv === :iu && (@printf("%d/%d %.0fs\n", k, length(pop), time() - t0); flush(stdout))
end
tag = length(ARGS) > 1 ? "_$(ARGS[2])" : ""
serialize(joinpath(@__DIR__, "data", "pop_$(NAME)$(tag).jls"), recs)
for conv in CONVS
    rs = filter(r -> r[:conv] === conv, recs)
    err = count(r -> haskey(r, :error), rs)
    ok = filter(r -> !haskey(r, :error), rs)
    feas = filter(r -> haskey(r, :cone_v), ok)
    uniq = filter(r -> r[:cone_spread] < 1e-20, feas)
    @printf("%s conv (%s): %d mechanisms, errors %d; cone empty %d; unique %d (= deletion %d); non-unique %d; sequestration possible %d; naive 0/0 %d; RREF order-dependent %d; oracle err max %.1e\n",
            NAME, conv, length(rs), err, length(ok) - length(feas), length(uniq),
            count(r -> r[:cone_vs_del] < 1e-20, uniq), length(feas) - length(uniq),
            count(r -> r[:cone_seq], feas), count(r -> r[:zero_over_zero], feas),
            count(r -> r[:order_spread] > 1e-12, ok),
            maximum((r[:oracle_err] for r in feas); init = 0.0))
    for r in filter(r -> haskey(r, :error), rs)[1:min(end, 5)]
        println("   error ", r[:id], ": ", r[:error])
    end
end
