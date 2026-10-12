# ABOUTME: Package-side checks of the guard probe's claims: the package's own law vs the
# ABOUTME: independent oracle on the face, enumerator emission, NaN/zero laws in the fitter.
using EnzymeRates, Serialization, Random, Printf, LinearAlgebra
const ER = EnzymeRates
include(joinpath(@__DIR__, "voracle.jl"))
const R = "/tmp/claude-501/-home-denis-linux--julia-dev-EnzymeRates/12f42995-f5ba-4fda-9510-dd9b87f7d069/scratchpad/round4/guard"

errline(e) = replace(first(sprint(showerror, e), 220), '\n' => ' ')

# ─── per-step mass-action constants of a package mechanism (plain) ──────────────────────
evalx(x::Symbol, P) = P[x]
evalx(x::Real, P) = big(x)
evalx(e::Expr, P) = getfield(Base, e.args[1])((evalx(a, P) for a in e.args[2:end])...)

"""Mass-action steps of `mech` with the constants the package's law uses for `params`
(fitted values, Keq): dependents from `_dependent_param_exprs`, absorbed Wegscheider
names through `_build_wegscheider_rename_map`, RE ratios as `_re_weight_ratio` reads them."""
function mass_action(mech::ER.Mechanism, params)
    dep, _ = ER._dependent_param_exprs(mech)
    rename = ER._build_wegscheider_rename_map(mech)
    P = Dict{Symbol, BigFloat}(k => big(v) for (k, v) in pairs(params))
    for (k, e) in dep
        P[k] = evalx(e, P)
    end
    look(k) = haskey(P, k) ? P[k] : P[rename[k]]
    map(ER._flat_steps(mech)) do (s, g)
        rep = first(ER.steps(mech)[g])
        f, t = string(ER.name(ER.from_species(s))), string(ER.name(ER.to_species(s)))
        c, r = Symbol[ER.name(x) for x in ER.consumed(s)], Symbol[ER.name(x) for x in ER.released(s)]
        if ER.is_equilibrium(s)
            K = look(ER.name(ER.Krapid(rep, :None), mech))
            kf, kr = ER.is_binding(s) ? (big(1), K) : (K, big(1))
        else
            kf, kr = look(ER.name(ER.Kfor(rep, :None), mech)), look(ER.name(ER.Krev(rep, :None), mech))
        end
        MStep(f, t, c, r, ER.is_equilibrium(s), kf, kr)
    end
end

"""Package law vs the oracle with the same constants: at a positive point (Λ → ∞), on the
faces `faces` (law value, interior limit, exact-face range over speed patterns)."""
function compare(label, em, faces, ref, rng; npar = 4)
    mech = ER.Mechanism(em)
    names = ER.fitted_params(em)
    mets = collect(ER.metabolites(em))
    println("\n=== $label  (fitted params: $(length(names)))")
    for ip in 1:npar
        params = (; (x => 10.0^(2rand(rng) - 1) for x in names)..., Keq = 10.0^(2rand(rng) - 1),
                  E_total = 1.0)
        steps = mass_action(mech, params)
        base = Dict{Symbol, BigFloat}(x => big(10.0^(2rand(rng) - 1)) for x in mets)
        nt(c) = (; (x => Float64(c[x]) for x in mets)...)
        vlaw = ER.rate_equation(em, nt(base), params)
        vor = rate_at(steps, base, big(10)^40, ones(BigFloat, length(steps)), ref)[1]
        @printf("  [param %d] positive point: law %.10g  oracle(Λ=1e40) %s  rel.diff %.2e\n", ip,
                vlaw, fmt(vor), abs(vlaw - Float64(vor)) / abs(vlaw))
        for zs in faces
            c0 = copy(base); foreach(z -> c0[z] = big(0), zs)
            lawface = ER.rate_equation(em, nt(c0), params)
            il = interior_limits(steps, base, zs, ref, rng)
            fl = face_limits(steps, base, zs, ref, rng)
            fv = fl.vals
            @printf("    face %-9s law %-12s interior limit %-12s exact face [%s, %s] (nterm %d) max|face-law|/|law| %s\n",
                    join(zs, ","), @sprintf("%.6g", lawface), join(unique(fmt.([r[end] for r in il])), "/"),
                    fmt(minimum(fv)), fmt(maximum(fv)), fl.nterm,
                    isnan(lawface) ? "law NaN" : fmt(maximum(abs.(fv .- lawface)) / abs(lawface)))
        end
    end
end

const RXN = @enzyme_reaction begin
    substrates: A[C], B[N]
    products: P[C], Q[N]
end

# ─── 1. the flagged enumerator children (probe's specs), through the public expand ───────
"""A package mechanism from the probe's plain-data spec, on reaction `rxn`."""
function mybuild(sp, rxn)
    mods = Set(x[1] for s in sp.steps for x in [s[3]; s[4]; s[1][2]; s[2][2]] if x[2])
    met(x) = x[2] ? ER.CompetitiveInhibitor(x[1]) :
             x[1] in sp.substrates ? ER.Substrate(x[1]) : ER.Product(x[1])
    form(f) = ER.Species(ER.Metabolite[met(x) for x in f[2]], f[1],
                         ER.Residual(ER.Substrate[ER.Substrate(x) for x in f[3]],
                                     ER.Product[ER.Product(x) for x in f[4]]))
    stp(s) = ER.Step(form(s[1]), form(s[2]), ER.Metabolite[met(x) for x in s[3]],
                     ER.Metabolite[met(x) for x in s[4]], s[5])
    groups = [ER.Step[stp(sp.steps[j]) for j in findall(==(g), sp.group)] for g in 1:maximum(sp.group)]
    isempty(sp.tags) && return ER.Mechanism(rxn, groups)
    sites = [ER.RegulatorySite(ER.AllostericRegulator[ER.AllostericRegulator(l) for l in s[1]],
                               s[3], collect(Symbol, s[2])) for s in sp.sites]
    ER.AllostericMechanism(rxn, groups, collect(Symbol, sp.tags), sp.multiplicity, sites)
end
dummy_rxn(sp) = begin
    mods = sort(unique([x[1] for s in sp.steps for x in [s[3]; s[4]; s[1][2]; s[2][2]] if x[2]]))
    ER.EnzymeReaction(
        [[ER.ReactantAtoms(ER.Substrate(s), [:C => length(sp.products)]) for s in sp.substrates];
         [ER.ReactantAtoms(ER.Product(p), [:C => length(sp.substrates)]) for p in sp.products]],
        [ER.RegulatorMults(ER.CompetitiveInhibitor(r), [1], :unspecified) for r in mods], [1])
end

function check_flagged()
    fl = deserialize(joinpath(R, "data", "enum_flagged.jls"))
    println("flagged children in the probe's file: $(length(fl))")
    seeds = ER.init_mechanisms(RXN)
    nexp = 0; nmove = 0; nseed = 0; nbibi = 0; ndeg_child = 0; nerr = 0
    for (i, r) in enumerate(fl)
        sp = r.parent
        bibi = isempty(sp.tags) && sp.substrates == [:A, :B] && sp.products == [:P, :Q] &&
               !any(s -> any(x -> x[2], [s[3]; s[4]; s[1][2]; s[2][2]]), sp.steps)
        rxn = bibi ? RXN : dummy_rxn(sp)
        try
            p, c = mybuild(r.parent, rxn), mybuild(r.child, rxn)
            inmove = any(==(c), ER._expand_re_to_ss(p))
            inexp = any(==(c), ER.expand_mechanisms(typeof(p)[p], rxn))
            isseed = bibi && any(==(p), seeds)
            dc = ER._degenerate(c)
            nmove += inmove; nexp += inexp; nseed += isseed; nbibi += bibi; ndeg_child += dc
            @printf("  %2d %-6s %-5s plain-bibi=%-5s parent∈init_mechanisms=%-5s in _expand_re_to_ss=%-5s in expand_mechanisms=%-5s child degenerate=%-5s %s\n",
                    i, r.level, r.move, bibi, isseed, inmove, inexp, dc, first(r.face, 50))
        catch e
            nerr += 1
            println("  $i error: ", errline(e))
        end
    end
    println("summary: in _expand_re_to_ss $nmove, in expand_mechanisms $nexp, parent is an init_mechanisms seed $nseed of $nbibi plain bi-bi, child degenerate $ndeg_child, errors $nerr")
end

# ─── 2. the mechanisms themselves, written in the DSL ─────────────────────────────────────
const D_parent = @enzyme_mechanism begin
    substrates: A, B
    products: P, Q
    steps: begin
        (E + A <--> E(A), E(B) + A <--> E(P, Q))
        (E + B ⇌ E(B), E(A) + B ⇌ E(P, Q))
        E + Q ⇌ E(Q)
        E(B) + Q ⇌ E(B, Q)
        E(Q) + B ⇌ E(B, Q)
        E(Q) + P <--> E(P, Q)
    end
end
const D_child = @enzyme_mechanism begin
    substrates: A, B
    products: P, Q
    steps: begin
        (E + A <--> E(A), E(B) + A <--> E(P, Q))
        (E + B <--> E(B), E(A) + B <--> E(P, Q))
        E + Q ⇌ E(Q)
        E(B) + Q ⇌ E(B, Q)
        E(Q) + B ⇌ E(B, Q)
        E(Q) + P <--> E(P, Q)
    end
end
const C_ldh_like = @enzyme_mechanism begin
    substrates: A, B
    products: P, Q
    steps: begin
        E + A ⇌ E(A)
        E(A) + B <--> E(A, B)
        E(A, B) <--> E(P, Q)
        E(P) + Q ⇌ E(P, Q)
        E(Q) + P ⇌ E(P, Q)
        E + P <--> E(P)
        E + Q ⇌ E(Q)
    end
end
const LDH_plain = @enzyme_mechanism begin
    substrates: NADH, Pyruvate
    products: Lactate, NAD
    steps: begin
        (E + Lactate <--> E(Lactate), E(NADH) + Lactate <--> E(Lactate, NADH))
        (E + NAD ⇌ E(NAD), E(Lactate) + NAD ⇌ E(Lactate, NAD), E(Pyruvate) + NAD ⇌ E(NAD, Pyruvate))
        (E + NADH <--> E(NADH), E(Lactate) + NADH <--> E(Lactate, NADH),
         E(Pyruvate) + NADH <--> E(NADH, Pyruvate))
        (E + Pyruvate ⇌ E(Pyruvate), E(NAD) + Pyruvate ⇌ E(NAD, Pyruvate),
         E(NADH) + Pyruvate ⇌ E(NADH, Pyruvate))
        E(NAD) + Lactate ⇌ E(Lactate, NAD)
        E(NADH, Pyruvate) <--> E(Lactate, NAD)
    end
end
const LDH_spec40 = @allosteric_mechanism begin
    substrates: NADH, Pyruvate
    products: Lactate, NAD
    catalytic_multiplicity: 4
    catalytic_steps: begin
        (E + Lactate <--> E(Lactate),
         E(NADH) + Lactate <--> E(Lactate, NADH)) :: NonequalAI
        (E + NAD ⇌ E(NAD),
         E(Lactate) + NAD ⇌ E(Lactate, NAD),
         E(Pyruvate) + NAD ⇌ E(NAD, Pyruvate)) :: EqualAI
        (E + NADH <--> E(NADH),
         E(Lactate) + NADH <--> E(Lactate, NADH),
         E(Pyruvate) + NADH <--> E(NADH, Pyruvate)) :: OnlyA
        (E + Pyruvate ⇌ E(Pyruvate),
         E(NAD) + Pyruvate ⇌ E(NAD, Pyruvate),
         E(NADH) + Pyruvate ⇌ E(NADH, Pyruvate)) :: EqualAI
        E(NAD) + Lactate ⇌ E(Lactate, NAD) :: EqualAI
        E(NADH, Pyruvate) <--> E(Lactate, NAD) :: OnlyA
    end
end
const Spec14 = @enzyme_mechanism begin
    substrates: S1, S2
    products: P
    steps: begin
        E + S1 <--> E(S1)
        E + S2 ⇌ E(S2)
        E(S1) + S2 ⇌ E(S1, S2)
        E(S2) + S1 ⇌ E(S1, S2)
        E(S1, S2) <--> E(P)
        E(P) ⇌ E + P
    end
end
const G_dead_end = @enzyme_mechanism begin
    substrates: S
    products: P, Q
    steps: begin
        E + S ⇌ E(S)
        E(S) <--> E(P, Q)
        E(Q) + P ⇌ E(P, Q)
        E + Q <--> E(Q)
        E(P) + Q <--> E(P, Q)
    end
end

function main()
    rng = Xoshiro(77)
    check_flagged()

    println("\n## D: hand-written parent/child, through the package moves")
    p = ER.Mechanism(RXN, ER.steps(ER.Mechanism(D_parent)))
    c = ER.Mechanism(RXN, ER.steps(ER.Mechanism(D_child)))
    println("  parent ∈ init_mechanisms(A+B→P+Q): ", any(==(p), ER.init_mechanisms(RXN)),
            "; _degenerate(parent) = ", ER._degenerate(p), "; _degenerate(child) = ", ER._degenerate(c))
    println("  child ∈ expand_mechanisms([parent]): ", any(==(c), ER.expand_mechanisms([p], RXN)))

    compare("D enumerator flip child 1 (package law)", D_child, [[:P, :Q], [:Q]], :A, rng)
    compare("C hand LDH-like trap (package law)", C_ldh_like, [[:P, :Q], [:Q]], :A, rng)
    compare("LDH spec 40, active conformation as a plain mechanism (package law)", LDH_plain,
            [[:Lactate, :NAD], [:NAD]], :NADH, rng)
    compare("Spec 14 'Numerator: redundant SS-bind' (package law)", Spec14, [[:S1, :S2]], :S1, rng)
    compare("G hand ss_dead_end_at_face (package law)", G_dead_end, [[:P, :Q]], :S, rng; npar = 2)

    println("\n## LDH spec 40 as written (allosteric, n = 4): rate_equation on the forward face")
    names = ER.fitted_params(LDH_spec40)
    for _ in 1:3
        params = (; (x => 10.0^(2rand(rng) - 1) for x in names)..., Keq = 10.0^(2rand(rng) - 1), E_total = 1.0)
        conc = (; NADH = 10.0^(2rand(rng) - 1), Pyruvate = 10.0^(2rand(rng) - 1), Lactate = 0.0, NAD = 0.0)
        @printf("  v(Lactate = NAD = 0) = %.6g\n", ER.rate_equation(LDH_spec40, conc, params))
    end
end

main()
