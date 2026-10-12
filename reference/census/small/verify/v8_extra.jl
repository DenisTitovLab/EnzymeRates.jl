# ABOUTME: Extra checks: the structural cause of D pairs read from the package types (tags differing only
# ABOUTME: outside the inactive state), and faces with a modifier at zero for the inhibitor reactions.
include(joinpath(@__DIR__, "vlib.jl"))
include(joinpath(@__DIR__, "moracle.jl"))

"""Groups with a step on the forms the free form reaches over non-:OnlyA steps."""
function igroups(m)
    ms = mstruct(m)
    kept = [ms.tags[s.g] !== :OnlyA for s in ms.steps]
    r = reach_from(ms, ms.frees[1], kept)
    Set(s.g for (k, s) in enumerate(ms.steps) if kept[k] && r[s.a] && r[s.b])
end

function dcause(name)
    d, mechs = loadrun(name)
    first_key = Dict{String, Int}()
    n = 0; tagonly = 0; same_steps = 0
    moves = Dict{Tuple{Symbol, Symbol}, Int}()
    for r in d.recs
        (r[:status] === :ok && r[:fitted]) || continue
        k = r[:key]
        haskey(first_key, k) || (first_key[k] = r[:gid]; continue)
        x, y = mechs[r[:gid]], mechs[first_key[k]]
        n += 1
        mv = (r[:move], d.recs[first_key[k]][:move])
        moves[mv] = get(moves, mv, 0) + 1
        (x isa ER.AllostericMechanism && y isa ER.AllostericMechanism) || continue
        ER.steps(x) == ER.steps(y) || continue
        same_steps += 1
        ER.catalytic_multiplicity(x) == ER.catalytic_multiplicity(y) &&
            ER.regulatory_sites(x) == ER.regulatory_sites(y) || continue
        tx, ty = ER.cat_allo_states(x), ER.cat_allo_states(y)
        gi = union(igroups(x), igroups(y))
        all(g -> tx[g] == ty[g], gi) && (tagonly += 1)
    end
    @printf("[%s] D pairs %d: same catalytic steps %d; tags differ only on groups outside both inactive states %d (%.1f%%)\n",
            name, n, same_steps, tagonly, 100tagonly / max(n, 1))
    println("   (this move, earlier move): ", join(("$(k[1])←$(k[2]) $v" for (k, v) in sort!(collect(moves); by = last, rev = true)), ", "))
end

for name in ["uni_bi_rxn", "bi_bi_rxn", "bi_bi_pp_rxn", "bi_bi_rxn.r2_100000", "bi_bi_pp_rxn.r2_1200"]
    dcause(name)
    GC.gc()
end

# ── faces with a modifier at zero (dead-end / competitive inhibitors), n = 1 oracle ─────
function modifier_faces(name)
    d, mechs = loadrun(name)
    rng = Xoshiro(818181)
    nfail = 0; ntest = 0
    for r in d.recs
        (r[:status] === :ok && r[:fitted]) || continue
        ms = mstruct(mechs[r[:gid]])
        mods = [x for x in ms.mets if !(x in ms.subs) && !(x in ms.prods)]
        isempty(mods) && continue
        ntest += 1
        A = aconsts(ms, rng)
        I = ms.allo ? iconsts(ms, A, rng) : nothing
        lnkI = I === nothing ? nothing : I.lnk(I.draw(rng))
        L = big(lu10(rng, -1, 1))
        conc = Dict{Symbol, BigFloat}(x => big(lu10(rng, -1, 1)) for x in ms.mets)
        σs = [[big(lu10(rng, -3, 3)) for _ in 1:(length(ms.steps) + 1)] for _ in 1:2]
        bad = String[]
        for zs in [[mods]; [[x] for x in mods if length(mods) > 1]]
            # the modifier at zero alone, and also with every product at zero
            for extra in (Symbol[], ms.prods)
                c0 = copy(conc); foreach(z -> c0[z] = big(0), [zs; extra])
                ce = copy(conc); foreach(z -> ce[z] = big(10)^-25 * big(lu10(rng, -1, 1)), zs)
                foreach(z -> ce[z] = big(0), extra)
                ex = [oracle_rate(ms, A, I, lnkI, L, c0, big(10)^40, σ) for σ in σs]
                it = oracle_rate(ms, A, I, lnkI, L, ce, big(10)^60, σs[1])
                (rel(ex[1], ex[2]) > 1e-8 || rel(ex[1], it) > 1e-8) &&
                    push!(bad, join(zs, ",") * (isempty(extra) ? "" : "+products") * "=0")
            end
        end
        isempty(bad) || (nfail += 1; println("   #$(r[:gid]) fails at ", join(bad, "; "), ": ", r[:text]))
    end
    println("[$name] modifier faces: $ntest candidates with a modifier step tested, $nfail fail")
end
modifier_faces("uni_uni_with_reg")
modifier_faces("uni_uni_reg_and_inhibitor")
