# ABOUTME: Follow-up B1: enumerates P5a's own cone directions for one mechanism and convention, prints the
# ABOUTME: moved exponents per step and the leading-face limit of every direction. Usage: julia cone_probe.jl pop id conv
include(joinpath(@__DIR__, "b1_lib.jl"))
include(joinpath(S, "proto/p5a/limits_lib.jl"))
name, id, conv = ARGS[1], ARGS[2], Symbol(ARGS[3])
m = only(filter(x -> x.id == id, load_allo(name)))
untie = conv === :iu; c = untie ? :i : conv === :iiu ? :ii : conv
untie |= conv === :iiu
J = MWC.joint_solve(m; order = :pkg, untie)
dirs, d = MWC.cone_directions(J, c)
println("$id conv $conv: cone dim $d, $(length(dirs)) directions")
sig = MWC.cone_signs(J, c)
O = sort!(collect(keys(sig)))
println("moved columns: ", join(["$(J.cols[i].j).$(J.cols[i].r)$(J.cs[:I][J.cols[i].j].re ? "K" : "k")($(sig[i]))" for i in O], " "))
rng = Random.Xoshiro(11)
hat = MWC.random_hat(J, rng); keq, L = logu(rng), logu(rng)
rA, rI = draw_reg(m, rng)
concs = Dict(x => logu(rng) for x in m.metabolites)
Sx = sym_setup(J, hat, keq, L, rA, rI, concs)
vdel = v_deletion(J, hat, keq, L, rA, rI, concs)
vals = BigFloat[]
for (k, a) in enumerate(dirs)
    s = sym_eval(Sx, a)
    push!(vals, s.v)
    k <= parse(Int, get(ARGS, 4, "12")) || continue
    seq = s.ordI[2] !== nothing && s.ordI[3] !== nothing && s.ordI[2] < s.ordI[3]
    @printf("  dir %3d: %s | v=%.6e  vs_del=%.1e seq=%s ordI=%s\n", k, join([string(Int(a[i])) for i in O], " "),
            Float64(s.v), reldiff(s.v, vdel, 0), seq, [o === nothing ? "-" : string(o[1]) for o in s.ordI])
end
spread = maximum(reldiff(v, vals[1], 0) for v in vals)
@printf("ALL %d directions: spread %.2e; max vs deletion %.2e; distinct values %d\n", length(vals), spread,
        maximum(reldiff(v, vdel, 0) for v in vals), length(unique(round.(Float64.(vals); sigdigits = 10))))
