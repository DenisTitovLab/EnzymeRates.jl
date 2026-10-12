# ABOUTME: Oracle-only census of the face problem on enumerated plain mechanisms and their
# ABOUTME: children, on full faces and on one-reactant-zero faces, with the package's constants.
include(joinpath(@__DIR__, "vpkglib.jl"))
const SP = "/tmp/claude-501/-home-denis-linux--julia-dev-EnzymeRates/12f42995-f5ba-4fda-9510-dd9b87f7d069/scratchpad"
include(joinpath(SP, "proto/common/Population.jl"))
include(joinpath(SP, "slice/src/Slice.jl"))
const S = Slice
setprecision(BigFloat, 512)

"""The package mechanism of a slice mechanism (plain), on the DSL's placeholder reaction."""
function to_pkg(m::S.Mechanism)
    (subs, prods, steps, group), _ = S.raw(m)
    mods = sort(unique([x.name for s in steps for x in [s.consumed; s.released; s.from.bound; s.to.bound]
                        if x.modifier]))
    met(x) = x.modifier ? ER.CompetitiveInhibitor(x.name) :
             x.name in subs ? ER.Substrate(x.name) : ER.Product(x.name)
    form(f) = ER.Species(ER.Metabolite[met(x) for x in f.bound], f.conf,
                         ER.Residual(ER.Substrate[ER.Substrate(x) for x in f.added],
                                     ER.Product[ER.Product(x) for x in f.subtracted]))
    stp(s) = ER.Step(form(s.from), form(s.to), ER.Metabolite[met(x) for x in s.consumed],
                     ER.Metabolite[met(x) for x in s.released], s.re)
    rxn = ER.EnzymeReaction(
        [[ER.ReactantAtoms(ER.Substrate(s), [:C => length(prods)]) for s in subs];
         [ER.ReactantAtoms(ER.Product(p), [:C => length(subs)]) for p in prods]],
        [ER.RegulatorMults(ER.CompetitiveInhibitor(r), [1], :unspecified) for r in mods], [1])
    ER.Mechanism(rxn, [ER.Step[stp(steps[j]) for j in findall(==(g), group)] for g in 1:maximum(group)])
end

"""Per face (each reactant alone at zero, each side at zero): the largest relative gap
between the exact-face RE limit (3 speed patterns, Λ = 1e30) and the limit from inside
(Λ = 1e40 at 1e-15 of the face metabolites), or :nonunique."""
function face_gaps(mech, rng)
    rxn = ER.reaction(mech)
    subs = [ER.name(x) for x in ER.substrates(rxn)]
    prods = [ER.name(x) for x in ER.products(rxn)]
    _, indep = ER._dependent_param_exprs(mech)
    params = (; (x => 10.0^(2rand(rng) - 1) for x in indep)..., Keq = 10.0^(2rand(rng) - 1), E_total = 1.0)
    steps = mass_action(mech, params)
    allm = unique!([x for s in steps for x in [s.cons; s.rel]])
    base = Dict{Symbol, BigFloat}(x => big(10.0^(2rand(rng) - 1)) for x in allm)
    ref = subs[1]
    scale = abs(rate_at(steps, base, big(10)^30, ones(BigFloat, length(steps)), ref)[1])
    faces = unique!([[[p] for p in prods]; [prods]; [[s] for s in subs]; [subs]])
    out = Dict{Vector{Symbol}, Any}()
    for zs in faces
        c = copy(base); foreach(z -> c[z] = big(10)^-15 * lu10(rng, -1, 1), zs)
        vin = rate_at(steps, c, big(10)^40, ones(BigFloat, length(steps)), ref)[1]
        c0 = copy(base); foreach(z -> c0[z] = big(0), zs)
        gap = big(0); nonu = false
        for k in 1:3
            sp = k == 1 ? ones(BigFloat, length(steps)) : [lu10(rng, -3, 3) for _ in steps]
            v, _, _, nt = rate_at(steps, c0, big(10)^30, sp, ref)
            nt > 1 && (nonu = true; break)
            gap = max(gap, abs(v - vin) / max(abs(vin), big(1e-6) * scale, big(1e-300)))
        end
        out[zs] = nonu ? :nonunique : gap
    end
    out, subs, prods
end

const PLAIN_ENUM = ["init_uniuni", "init_unibi", "init_bibi", "init_bibi_pp", "expand_bibi", "p1_enum_inh"]

function run(label, ms, rng, flagged)
    cnt = Dict{String, Int}(); ex = Dict{String, Vector{Any}}()
    nflag_hit = 0; nflag_seen = 0
    for m in ms
        S.is_allosteric(m) && continue
        k, pm = try
            pm = to_pkg(m)
            g, subs, prods = face_gaps(pm, rng)
            full = [prods, subs]
            bad(z) = g[z] === :nonunique || g[z] > 1e-6
            kf = any(z -> bad(z), full) ? "full face" : any(z -> bad(z), keys(g)) ? "one-reactant face only" : "ok"
            (kf, pm)
        catch e
            ("error: " * first(errline(e), 50), nothing)
        end
        cnt[k] = get(cnt, k, 0) + 1
        length(get!(ex, k, Any[])) < 3 && push!(ex[k], m)
        if pm !== nothing && any(==(pm), flagged)
            nflag_seen += 1; k == "ok" || (nflag_hit += 1)
        end
    end
    println("  $label: ", sort!(collect(cnt)), "  (probe-flagged among them: $nflag_seen, of which the oracle finds a face gap in $nflag_hit)")
    for k in ("one-reactant face only", "full face")
        haskey(ex, k) || continue
        for m in ex[k][1:min(2, end)]
            println("    [$k] ", join([S.step_text(m, e) for e in m.edges], "; "))
        end
    end
end

function main()
    rng = Xoshiro(2026)
    fl = deserialize(joinpath(R, "data", "enum_flagged.jls"))
    flagged = [mybuild(r.child, dummy_rxn(r.child)) for r in fl if isempty(r.child.tags)]
    parents = [c.mech for set in PLAIN_ENUM for c in deserialize(joinpath(SP, "slice/core/data/$set.jls"))]
    t0 = time()
    run("enumerated plain parents ($(length(parents)))", parents, rng, flagged)
    @printf("  (%.0f s)\n", time() - t0)
    kids = S.Mechanism[]
    for m in parents, f in (S.flip_children, S.split_children)
        try append!(kids, f(m)) catch end
    end
    unique!(kids)
    run("their flip and split children ($(length(kids)))", kids, rng, flagged)
    @printf("  (%.0f s)\n", time() - t0)
    sample = kids[randperm(rng, length(kids))[1:min(600, length(kids))]]
    kids2 = S.Mechanism[]
    for m in sample, f in (S.flip_children, S.split_children)
        try append!(kids2, f(m)) catch end
    end
    unique!(kids2)
    run("flip and split children of 600 random children ($(length(kids2)))", kids2, rng, flagged)
    @printf("  (%.0f s)\n", time() - t0)
end

main()
