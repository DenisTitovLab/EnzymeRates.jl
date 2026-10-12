# ABOUTME: Validates the two-conformation oracle (P5a): against the package where no group is
# ABOUTME: :OnlyA, explicit vs factorized networks, and the deletion oracle on every :OnlyA population.
include(joinpath(@__DIR__, "pkgvals.jl"))
using Printf

const Λ = MWC.big1024(10)^40
relerr(a, b) = abs(a - b) / max(abs(b), MWC.big1024(1e-300))

pops = Dict(name => filter(m -> m.allosteric, Pop.load_population(name))
            for name in ("specs", "handwritten", "allo_uniuni", "allo_bibi", "expand_bibi"))

"""Oracle inputs of a package sample: per-step constants, deletion activity, regulators."""
function inputs(m, s)
    vals = pkg_values(m, s)
    kA = pkg_step_constants(m, vals, :A)
    kI = pkg_step_constants(m, vals, :I)
    tags = MWC.step_tags(m)
    activeI = [tags[j] !== :OnlyA && kI[j] !== nothing for j in eachindex(m.steps)]
    kIf = [kI[j] === nothing ? [big(1), big(1)] : kI[j] for j in eachindex(m.steps)]
    keepI = MWC.reachable(m, [tags[j] !== :OnlyA for j in eachindex(m.steps)])
    regA, regI = pkg_regK(m, vals)
    concs = Dict(k => MWC.big1024(v) for (k, v) in s.concs)
    vals, kA, kIf, activeI, keepI, regA, regI, concs
end

println("== V1: no-:OnlyA mechanisms, factorized oracle at the mechanism's n vs package")
worst_f, worst_p = 0.0, 0.0
for name in ("specs", "handwritten", "allo_uniuni"), m in pops[name]
    :OnlyA in m.group_tags && continue
    ef, ep = 0.0, 0.0
    for s in m.samples
        vals, kA, kI, act, keep, rA, rI, c = inputs(m, s)
        reI = [st.re for st in m.steps]
        v = MWC.v_factorized(m, kA, kI, reI, act, keep, vals[:L], rA, rI, c;
                             n = m.cat_multiplicity, Λ = Λ)
        vp = pkg_poly_v(m, vals, Dict(s.concs), rA, rI)
        ef = max(ef, Float64(relerr(v, MWC.big1024(s.v)))); ep = max(ep, Float64(relerr(v, vp)))
    end
    global worst_f = max(worst_f, ef); global worst_p = max(worst_p, ep)
    @printf("  %-16s n=%d sites=%d  vs package v %.1e   vs exact package polys %.1e\n",
            m.id, m.cat_multiplicity, length(m.sites), ef, ep)
end
@printf("  worst: vs package Float64 v %.2e, vs exact polys %.2e\n", worst_f, worst_p)

println("== V2: explicit 2-protomer network vs factorized formula (n = 2)")
for id in ("specs/27", "handwritten/25", "handwritten/34", "allo_uniuni/29", "specs/28")
    m = only(filter(x -> x.id == id, reduce(vcat, values(pops))))
    e = 0.0
    for s in m.samples[1:2]
        vals, kA, kI, act, keep, rA, rI, c = inputs(m, s)
        reI = [st.re for st in m.steps]
        vf = MWC.v_factorized(m, kA, kI, reI, act, keep, vals[:L], rA, rI, c; n = 2, Λ = Λ)
        t = @elapsed ve = MWC.v_explicit(m, kA, kI, reI, act, keep, vals[:L], rA, rI, c;
                                         nprot = 2, Λ = Λ)
        e = max(e, Float64(relerr(ve, vf)))
    end
    @printf("  %-16s tags=%s sites=%d  explicit vs factorized %.1e\n", id, m.group_tags,
            length(m.sites), e)
end

println("== V3: explicit 1-protomer network (regulators explicit) vs factorized at n = 1")
for name in ("specs", "handwritten", "allo_uniuni"), m in pops[name]
    nreg = isempty(m.sites) ? 1 : length(MWC.reg_states(m))
    nreg * length(m.species) * 2 > 400 && continue
    e = 0.0
    for s in m.samples[1:2]
        vals, kA, kI, act, keep, rA, rI, c = inputs(m, s)
        reI = [st.re for st in m.steps]
        vf = MWC.v_factorized(m, kA, kI, reI, act, keep, vals[:L], rA, rI, c; n = 1, Λ = Λ)
        ve = MWC.v_explicit(m, kA, kI, reI, act, keep, vals[:L], rA, rI, c; nprot = 1, Λ = Λ)
        e = max(e, Float64(relerr(ve, vf)))
    end
    e > 1e-30 && @printf("  %-16s regstates=%d  explicit vs factorized %.1e\n", m.id, nreg, e)
end
println("  (only mechanisms above 1e-30 are listed)")
nv3 = sum(count(m -> (isempty(m.sites) ? 1 : length(MWC.reg_states(m))) * length(m.species) * 2 <= 400, pops[n]) for n in ("specs", "handwritten", "allo_uniuni"))
println("  checked: $nv3 mechanisms")

println("== V4: deletion oracle (the package's semantics) vs package v on :OnlyA mechanisms")
for name in ("specs", "handwritten", "allo_uniuni", "expand_bibi", "allo_bibi")
    worst = 0.0; nbad = 0; nm = 0; bad = String[]
    t = @elapsed for m in pops[name]
        :OnlyA in m.group_tags || continue
        nm += 1
        e = 0.0
        for s in m.samples
            vals, kA, kI, act, keep, rA, rI, c = inputs(m, s)
            reI = [st.re for st in m.steps]
            v = MWC.v_factorized(m, kA, kI, reI, act, keep, vals[:L], rA, rI, c;
                                 n = m.cat_multiplicity, Λ = Λ)
            e = max(e, Float64(relerr(v, MWC.big1024(s.v))))
        end
        worst = max(worst, e)
        e > 1e-10 && (nbad += 1; length(bad) < 8 && push!(bad, @sprintf("%s(%.1e)", m.id, e)))
    end
    @printf("  %-12s :OnlyA mechanisms %4d  worst rel %.1e  above 1e-10: %d %s  (%.1f s)\n",
            name, nm, worst, nbad, join(bad, " "), t)
end
