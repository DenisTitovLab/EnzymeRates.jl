# ABOUTME: Skeptic P5a: for the cone records with the largest oracle-vs-symbolic gap, checks that the
# ABOUTME: oracle converges to the symbolic limit as ε shrinks (the gap must fall like ε).
include(joinpath(@__DIR__, "limits_lib.jl"))
using Printf
for (name, conv) in (("allo_terbi", :iu), ("allo_bibi_pp", :i), ("allo_unibi", :i))
    pop = Dict(m.id => m for m in deserialize(joinpath(@__DIR__, "data", "$name.jls")))
    recs = filter(r -> r[:conv] === conv && haskey(r, :oracle_err) && !isnan(r[:oracle_err]),
                  deserialize(joinpath(@__DIR__, "data", "res_$name.jls")))
    sort!(recs; by = r -> -r[:oracle_err])
    for r in recs[1:2]
        m = pop[r[:id]]
        untie = conv === :iu
        rng = Random.Xoshiro(11)
        J = MWC.joint_solve(m; order = :pkg, untie)
        hat = MWC.random_hat(J, rng); keq, L = logu(rng), logu(rng)
        rA, rI = draw_reg(m, rng)
        concs = Dict(x => logu(rng) for x in m.metabolites)
        dirs, _ = MWC.cone_directions(J, :i)
        a = dirs[1]
        S = sym_setup(J, hat, keq, L, rA, rI, concs)
        vs = sym_eval(S, a).v
        v1 = v_path(J, hat, zeros(MWC.Q, length(J.cols)), big(1), keq, L, rA, rI, concs)
        errs = [(k, reldiff(v_path(J, hat, a, MWC.big1024(10)^(-k), keq, L, rA, rI, concs), vs, v1)) for k in (3, 6, 9, 12)]
        @printf("  %-16s %-2s max|a|=%d  symbolic %.6e  |oracle−symbolic| at ε=1e-3,1e-6,1e-9,1e-12: %s\n", r[:id], conv,
                Int(maximum(abs, a)), Float64(vs), join([@sprintf("%.1e", e) for (_, e) in errs], ", "))
    end
end
