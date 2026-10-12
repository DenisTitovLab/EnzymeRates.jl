# ABOUTME: Follow-up B1: evaluates P5a's leading-face limit along a hand-built cone direction (exponents per
# ABOUTME: step slot), checking C·a = 0 and the sign constraints. Usage: julia hand_dir.jl pop id conv "j:e1[,e2] ..."
include(joinpath(@__DIR__, "full_eval.jl"))
name, id, conv = ARGS[1], ARGS[2], Symbol(ARGS[3])
m = only(filter(x -> x.id == id, load_allo(name)))
untie = conv in (:iu, :iiu, :Dk); c = conv === :iu ? :i : :ii
J = MWC.joint_solve(m; order = :pkg, untie)
for (j, s) in enumerate(m.steps)
    cs = J.cs[:I][j]
    println("  step $j: ", step_dsl(m, s), "  I:", cs.re ? "RE" : "SS", " canonical ", form_text(m.species[cs.a]), " -> ",
            form_text(m.species[cs.b]), " flipped=", cs.flipped, " tag=", m.group_tags[Pop.group_of_step(m)[j]])
end
sig = MWC.cone_signs(J, c)
rng = Random.Xoshiro(11)
hat = MWC.random_hat(J, rng); keq, L = logu(rng), logu(rng)
rA, rI = draw_reg(m, rng)
concs = Dict(x => logu(rng) for x in m.metabolites)
Sx = sym_setup(J, hat, keq, L, rA, rI, concs)
vdel = v_deletion(J, hat, keq, L, rA, rI, concs)
for spec in ARGS[4:end]
    a = zeros(MWC.Q, length(J.cols))
    for item in split(spec)
        j, es = split(item, ":")
        for (r, e) in enumerate(split(es, ","))
            a[J.idx[(:I, parse(Int, j), r)]] = parse(Int, e)
        end
    end
    okC = all(iszero, J.C * a)
    oks = all(i -> (s = get(sig, i, 0); s == 0 ? a[i] == 0 : sign(a[i]) == s), eachindex(a))
    s = sym_eval(Sx, a)
    seq = s.ordI[2] !== nothing && s.ordI[3] !== nothing && s.ordI[2] < s.ordI[3]
    @printf("  dir [%s]: C·a=0 %s, signs ok %s | v=%.6e (deletion %.6e, rel %.1e) seq=%s ordI=%s\n", spec, okC, oks,
            Float64(s.v), Float64(vdel), reldiff(s.v, vdel, 0), seq, [o === nothing ? "-" : string(o[1]) for o in s.ordI])
end
