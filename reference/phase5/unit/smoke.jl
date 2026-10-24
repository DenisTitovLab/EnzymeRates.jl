# ABOUTME: First run of the unit-rule research code: identifier checks for coefficient names,
# ABOUTME: the 42 specs through every rule and check, and a few printed laws to read.
include(joinpath(@__DIR__, "ulib.jl"))
include(joinpath(@__DIR__, "pops.jl"))

for n in ("d_P_B⁻¹", "d_A²_B", "n_B⁻¹", "dA_S", "nI_0", "d_A⁻²_B", "N_A", "x⁻", "d_P_B⁻¹⁰")
    ex = try Meta.parse("$n = 2") catch e e end
    println(rpad(n, 12), " isidentifier=", Base.isidentifier(n), "  parse -> ", repr(ex))
end
println("parse(\"a * d_P_B⁻¹ * P / B\") -> ", repr(Meta.parse("a * d_P_B⁻¹ * P / B")))

P = population("spec")
println(length(P), " specs")
rng = Xoshiro(1)
for (id, text, subs, prods) in P
    try
        m = parse_neutral(text, subs, prods)
        r = analyze(m, rng)
        rs = join(("$(t):$(x.spec_kind)/$(x.slice_kind)/cle=$(x.cleland)/sameC=$(x.same_cleland)/sameS=$(x.same_slice)" for (t, x) in r.rules), " ")
        @printf("%-58s ef=%.0e ec=%.0e ef2=%.0e ec2=%.0e es=%.0e esc=%.0e ref=%.0e/%.0e U=%s/%s sU=%s/%s art=%d/%d %s\n",
                first(id, 58), r.e_full, r.e_coef, r.e_full2, r.e_coef2, r.e_sfull, r.e_scoef,
                r.e_ref, r.e_ref2, r.units_full, r.units_coef, r.units_sfull, r.units_scoef,
                r.art, r.art2, rs)
        (r.units_full && r.units_coef) || println("   units: ", r.units_note)
    catch err
        println(first(id, 58), " ERROR ", first(sprint(showerror, err), 400))
    end
end
for (id, text, subs, prods) in P
    occursin(r"Ping Pong|MWC Dimer|^spec:Uni-Uni$", id) || continue
    m = parse_neutral(text, subs, prods)
    r = analyze(m, Xoshiro(2))
    println("\n=== ", id, "\n", r.full, "\n--- coefficient form\n", r.coef)
end
