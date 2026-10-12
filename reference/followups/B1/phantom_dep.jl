# ABOUTME: Follow-up B1: shows that a phantom's (Dk) limit depends on the blocked bridge step's inactive ratio
# ABOUTME: kf_I/kr_I and not on its speed: scales the hat values of that step's two constants. Usage: julia phantom_dep.jl id...
include(joinpath(@__DIR__, "full_eval.jl"))
for id in ARGS
    name = split(id, "/")[1]
    m = only(filter(x -> x.id == id, load_allo(String(name))))
    rng = Random.Xoshiro(11)
    J = MWC.joint_solve(m; order = :pkg, untie = true)
    hat = MWC.random_hat(J, rng); keq, L = logu(rng), logu(rng)
    rA, rI = draw_reg(m, rng)
    concs = Dict(x => logu(rng) for x in m.metabolites)
    br = bridge_kinetic_doors(m)
    extra = [(J.idx[(:I, j, 1)], J.idx[(:I, j, 2)]) for j in br]
    dirs, _ = cone_dirs(J, :ii; extra)
    vdel = v_deletion(J, hat, keq, L, rA, rI, concs)
    lim(h) = sym_eval(sym_setup(J, h, keq, L, rA, rI, concs), dirs[1]).v
    v0 = lim(hat)
    println("== $id: bridge kinetic doors ", [step_dsl(m, m.steps[j]) for j in br], "; $(length(dirs)) directions")
    @printf("   deletion %.6e | limit %.6e (rel to deletion %.1e)\n", Float64(vdel), Float64(v0), reldiff(v0, vdel, 0))
    for (lab, f1, f2) in (("kf_I ×10", 10, 1), ("kr_I ×10", 1, 10), ("both ×10", 10, 10))
        h = copy(hat)
        for j in br
            for (r, f) in ((1, f1), (2, f2))
                i = J.idx[(:I, j, r)]
                haskey(h, i) || error("column $(J.cols[i].name) is not independent")
                h[i] *= f
            end
        end
        v = lim(h)
        @printf("   %-9s limit %.6e (rel change %.1e)\n", lab, Float64(v), reldiff(v, v0, 0))
    end
end
