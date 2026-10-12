# ABOUTME: Recomputes census table rows end to end with a driver of its own: sampled seeds, base
# ABOUTME: tier, batch stage, expansion rounds; then fitted/allosteric, W_admit and D per round.
# usage: v5_e2e.jl <rxn> <N1> <N2>
const CEN = "/tmp/claude-501/-home-denis-linux--julia-dev-EnzymeRates/12f42995-f5ba-4fda-9510-dd9b87f7d069/scratchpad/census/large"
using EnzymeRates, Random, Serialization, Printf
const ER = EnzymeRates

# The reactions, as the census reads them: the test file's `const` definitions.
let ast = Meta.parseall(read(joinpath(dirname(dirname(pathof(EnzymeRates))), "test",
                                      "test_mechanism_enumeration.jl"), String))
    for ex in ast.args
        ex isa Expr && ex.head === :const && ex.args[1] isa Expr &&
            ex.args[1].args[1] in (:ter_bi_rxn, :ter_ter_rxn, :pyruvate_carboxylase_rxn) &&
            Core.eval(Main, ex)
    end
end

const SEED = 20261011

"""The batch stage of `_process_batch` on `mechs`: structural `seen`, the complexity
filter (337), the parameter cap (20), `rate_equation_string` and its dedup key; one
representative per key not in `memo` is fitted. Returns one record per fresh mechanism."""
function batch!(mechs, seen, memo)
    fresh = Any[]
    nseen = 0
    for m in mechs
        h = hash(m)
        h in seen ? (nseen += 1) : (push!(seen, h); push!(fresh, m))
    end
    recs = map(fresh) do m
        try
            ER._eq_complexity(m) > 337 && return (mech = m, status = :cx, eq = UInt64(0))
            em = ER.compile_mechanism(m)
            length(ER.fitted_params(em)) > 20 && return (mech = m, status = :cap, eq = UInt64(0))
            (mech = m, status = :ok, eq = ER._rate_eq_dedup_key(ER.rate_equation_string(em)))
        catch e
            (mech = m, status = :fail, eq = UInt64(0))
        end
    end
    fitted = falses(length(recs))
    inmemo = 0
    firsts = Dict{UInt64, Int}()
    for (i, r) in enumerate(recs)
        r.status === :ok || continue
        r.eq in memo && (inmemo += 1; continue)
        haskey(firsts, r.eq) || (firsts[r.eq] = i; fitted[i] = true)
    end
    union!(memo, keys(firsts))
    recs, fitted, (nseen = nseen, inmemo = inmemo)
end

"""The components of the species joined by steps outside `:OnlyA` groups that no chain of
`:OnlyA` plain bindings (unbound to bound) reaches from the free enzyme's component: the
island rule, read off the package mechanism."""
function lost_components(am::ER.AllostericMechanism)
    gs, tags = ER.steps(am), ER.cat_allo_states(am)
    sp = unique!([s for g in gs for st in g for s in (ER.from_species(st), ER.to_species(st))])
    at = Dict(s => i for (i, s) in enumerate(sp))
    parent = collect(eachindex(sp))
    top(i) = (while parent[i] != i; i = parent[i]; end; i)
    for (g, grp) in enumerate(gs), st in grp
        tags[g] === :OnlyA && continue
        parent[top(at[ER.from_species(st)])] = top(at[ER.to_species(st)])
    end
    comp = top.(eachindex(sp))
    arcs = [(comp[at[ER.from_species(st)]], comp[at[ER.to_species(st)]])
            for (g, grp) in enumerate(gs) if tags[g] === :OnlyA
            for st in grp if !ER._is_chemistry(st)]
    free = findall(s -> isempty(ER.bound(s)) && isempty(ER.residual(s)), sp)
    length(free) == 1 || return -1
    reached = Set([comp[free[1]]])
    grown = true
    while grown
        grown = false
        for (a, b) in arcs
            a in reached && !(b in reached) && (push!(reached, b); grown = true)
        end
    end
    length(setdiff(Set(comp), reached))
end

rxn_name = Symbol(ARGS[1])
N1, N2 = parse(Int, ARGS[2]), parse(Int, ARGS[3])
rxn = getfield(Main, rxn_name)
T0 = time()
seeds_all = unique!(collect(ER.init_mechanisms(rxn)))
n_all = length(seeds_all)
idx = n_all > 2000 ? sort!(randperm(Xoshiro(SEED), n_all)[1:2000]) : collect(1:n_all)
seeds = seeds_all[idx]
seeds_all = nothing
GC.gc()
ndeg = count(ER._degenerate, seeds)
base, bfail = ER._base_tier(seeds, rxn)
@printf("%s: %d seeds, sample %d, degenerate %d, base %d (flip children %d, failures %d) [%.0f s]\n",
        rxn_name, n_all, length(seeds), ndeg, length(base), length(base) - (length(seeds) - ndeg),
        length(bfail), time() - T0)
# The degenerate seeds the screen removes: growth of today's law (`rate_equation` in BigFloat)
# with the substrates (products at zero) and the products (substrates at zero) scaled by
# t = 1e10, 1e20, 1e30 along the diagonal and each axis.
setprecision(BigFloat, 256)
function grows(em, rng)
    mets = collect(ER.metabolites(em))
    names = collect(ER.fitted_params(em))
    p = (; (n => big(10.0)^(2rand(rng) - 1) for n in names)..., Keq = big(10.0)^(2rand(rng) - 1),
         E_total = big(1.0))
    base = Dict(x => big(10.0)^(2rand(rng) - 1) for x in mets)
    subs = [ER.name(x) for x in ER.substrates(rxn)]
    prods = [ER.name(x) for x in ER.products(rxn)]
    out = String[]
    for (side, grow, zero) in (("fwd", subs, prods), ("rev", prods, subs))
        for w in [[ones(length(grow))]; [[i == j ? 1.0 : 0.0 for j in eachindex(grow)] for i in eachindex(grow)]]
            vs = map((big(10.0)^10, big(10.0)^20, big(10.0)^30)) do t
                c = copy(base)
                foreach(z -> c[z] = big(0.0), zero)
                foreach(((x, e),) -> c[x] = base[x] * t^e, zip(grow, w))
                Base.invokelatest(ER.rate_equation, em, (; (x => c[x] for x in mets)...), p)
            end
            abs(vs[3]) > 1.5 * abs(vs[2]) && abs(vs[2]) > 1.5 * abs(vs[1]) &&
                push!(out, "$side w=$w")
        end
    end
    out
end
let rng = Xoshiro(SEED + 9), ng = 0, nd = 0, ex = String[]
    for m in seeds
        ER._degenerate(m) || continue
        nd += 1
        g = try
            grows(ER.compile_mechanism(m), rng)
        catch e
            ["error: " * first(sprint(showerror, e), 80)]
        end
        isempty(g) || (ng += 1; length(ex) < 3 && push!(ex, join(g, "; ")))
    end
    @printf("degenerate seeds of the sample: %d; today's law grows without bound in %d %s [%.0f s]\n",
            nd, ng, ex, time() - T0)
end
seen, memo = Set{UInt64}(), Set{UInt64}()
rounds = Any[]
recs, fitted, c = batch!(base, seen, memo)
push!(rounds, (recs = recs, fitted = fitted, c = c, parents = Int[]))
for (rnd, NP) in ((1, N1), (2, N2))
    NP > 0 || break
    prev = rounds[end]
    pool = [r.mech for r in prev.recs if r.status === :ok]
    pick = length(pool) > NP ? sort!(randperm(Xoshiro(SEED + rnd), length(pool))[1:NP]) :
           collect(eachindex(pool))
    kids, owner = Any[], Dict{Any, Int}()
    nthrow, nemit = 0, 0
    for (pi, p) in enumerate(pool[pick])
        ks = try
            ER.expand_mechanisms(Union{ER.Mechanism, ER.AllostericMechanism}[p], rxn)
        catch
            nthrow += 1
            continue
        end
        nemit += length(ks)
        for k in ks
            haskey(owner, k) && continue
            owner[k] = pi
            push!(kids, k)
        end
    end
    rr, ff, cc = batch!(kids, seen, memo)
    @printf("round %d: %d parents of %d, %d threw, %d emitted, %d distinct [%.0f s]\n", rnd,
            length(pick), length(pool), nthrow, nemit, length(kids), time() - T0)
    push!(rounds, (recs = rr, fitted = ff, c = cc, parents = pick))
end

out = Any[]
for (ri, R) in enumerate(rounds)
    st = [r.status for r in R.recs]
    nf = count(R.fitted)
    fm = [r.mech for (r, f) in zip(R.recs, R.fitted) if f]
    nallo = count(m -> m isa ER.AllostericMechanism, fm)
    lost = [m isa ER.AllostericMechanism ? lost_components(m) : 0 for m in fm]
    @printf("ROUND %d: fresh %d (seen %d); complexity skip %d; cap %d; fail %d; compiled %d; eq_hash in memo %d, in-batch duplicate %d; FITTED %d (allosteric %d); island rule (own, package mechanism): %d with a lost component (%s), %d free-form errors\n",
            ri - 1, length(R.recs), R.c.nseen, count(==(:cx), st), count(==(:cap), st),
            count(==(:fail), st), count(==(:ok), st), R.c.inmemo,
            count(==(:ok), st) - R.c.inmemo - nf, nf, nallo, count(>(0), lost),
            string(round(100 * count(>(0), lost) / max(nf, 1); digits = 1), "%"),
            count(<(0), lost))
    push!(out, (fitted = fm, lost = lost, hashes = [hash(m) for m in fm]))
end
serialize(joinpath(CEN, "verify", "data", "v5_$(rxn_name).jls"), out)
@printf("total %.0f s, maxrss %.0f MB\n", time() - T0, Sys.maxrss() / 2^20)
