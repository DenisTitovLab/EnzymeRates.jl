# ABOUTME: R5.1 GC workloads for one law-type variant (VARIANT = U, A, B or C): a derivation loop and
# ABOUTME: an identify-like search (key cache, fits, beam per key), 3 rounds each, full GC between.
include(joinpath(@__DIR__, "lib.jl"))
const V = get(ENV, "VARIANT", "C")
const S = load_variant(V)
const TINY = get(ENV, "TINY", "0") == "1"
const NFIT = TINY ? 2 : 10
t0 = time()

pick(set) = group(set) in ("spec", "pingpong", "bibi", "twosites", "uniuni")
pop, sets = S.Mechanism[], String[]
nbad = 0
for r in census_rows()
    pick(r.set) || continue
    try
        push!(pop, slice_mech(S, r.net))
        push!(sets, r.set)
    catch
        global nbad += 1
    end
end
# identify seeds: small mechanisms, a third each from bi-bi init, ping-pong init and the
# allosteric ping-pong children, so the beam's members have flip and split children
seeds = reduce(vcat, [first([m for (m, s) in zip(pop, sets) if s == set && length(m.edges) <= 9 &&
                                 (set != "pingpong:L1" || S.is_allosteric(m))], TINY ? 2 : 20)
                      for set in ("bibi:init", "pingpong:init", "pingpong:L1")])
const POP = let io = IOBuffer(); serialize(io, TINY ? pop[1:40] : pop); take!(io) end
const SEEDS = let io = IOBuffer(); serialize(io, seeds); take!(io) end
println("julia ", VERSION, " variant ", V, ": population ", length(pop), " (",
        count(S.is_allosteric, pop), " allosteric; ", nbad, " rows the slice cannot build), ",
        "identify seeds ", length(seeds)); flush(stdout)

"""One derivation pass over the population read back from its bytes, reading the law's
fields in the loop as callers do (the crash probe's `L.D` read among them). The loop has no
try/catch, as the crash probe's loop had none; the population derives without error."""
function derive_round(bytes)
    tot = zeros(Int, 9)
    for m in deserialize(IOBuffer(bytes))
        sol = S.solve(m)
        L = S.rate_law(m, sol)
        nD = L isa S.Law ? S.nterms(L.D) : S.nterms(L.active.D)
        nsig = L isa S.Law ? sum(S.nterms, L.sigma; init = 0) :
                             sum(S.nterms, L.active.sigma; init = 0)
        nsite = 0
        if L isa S.AlloLaw
            for (WA, WI, k) in L.sites
                nsite += S.nterms(WA) + S.nterms(WI) + k
            end
            L.inactive === nothing || (nsite += S.nterms(L.inactive.D))
            nsite += S.nterms(L.fA) + S.nterms(L.fI)
        end
        s = S.rate_equation_string(m, L, sol)
        c = S.coefficient_form(m, L)
        key = S.canonical(L)
        code = S.body(S.program(m, sol, L))
        v5 = class5(m, L)
        tot .+= (1, nD, S.nterms(L.H), nsig, nsite, length(s) + length(c),
                 length(key.order), length(string(code)), isempty(v5) ? 0 : 1)
    end
    tot
end

"""Pass 1 of identify on `ms`: the complexity filter, then law, key, printed law and fitted
names; `(m, key, order, nchars, nfitted)` per kept candidate, `nothing` when skipped."""
function pass1(ms)
    map(ms) do m
        try
            S.complexity(m) > 337 && return nothing
            sol = S.solve(m)
            L = S.rate_law(m, sol)
            key = S.canonical(L)
            (m = m, key = key.key, order = key.order,
             n = length(S.rate_equation_string(m, L, sol)), np = length(S.fitted_params(m)))
        catch
            nothing
        end
    end
end

"""A fit's stand-in: the loss of the compiled law over seeded rows and parameter draws."""
function fit_loss(m, key)
    em = S.compile(m)
    ps, cs = S.parameters(m), S.metabolites(m)
    rng = Xoshiro(hash(key))
    loss = 0.0
    for _ in 1:20
        params = NamedTuple{Tuple(ps)}(Tuple(p === :E_total ? 1.0 : exp(randn(rng)) for p in ps))
        for _ in 1:25
            concs = NamedTuple{Tuple(cs)}(Tuple(exp(randn(rng)) for _ in cs))
            v = S.rate_equation(em, concs, params)
            loss += isfinite(v) ? log1p(abs(v)) : 1e3
        end
    end
    round(loss; digits = 6)
end

"""One identify-like search from the seeds read back from their bytes: pass 1, a cache of
fits by key, a beam of the three best keys whose members all expand (flip and split), pass 1
on the children, fits of new keys, and a cross-validation list deduplicated by key."""
function identify_round(bytes)
    cache, members, cv = Dict{String, Float64}(), Dict{String, Vector{Any}}(), String[]
    nfit = 0
    function absorb!(rows)
        for r in rows
            r === nothing && continue
            push!(get!(members, r.key, Any[]), r.m)
            haskey(cache, r.key) && continue
            nfit < NFIT || continue
            cache[r.key] = fit_loss(r.m, r.key)
            nfit += 1
        end
    end
    absorb!(pass1(deserialize(IOBuffer(bytes))))
    beam = first(sort!(collect(keys(cache)); by = k -> (cache[k], k)), 3)
    kids, xerr = Any[], String[]
    for k in beam, m in members[k]
        try
            append!(kids, S.flip_children(m))
            append!(kids, S.split_children(m))
        catch e
            push!(xerr, first(sprint(showerror, e), 100))
        end
    end
    absorb!(pass1(unique(kids)))
    for k in sort!(collect(keys(cache)); by = k -> (cache[k], k))
        k in cv || push!(cv, k)
    end
    (length(kids), length(members), nfit, round(sum(values(cache)); digits = 3), length(cv),
     sum(length, values(members)), unique(xerr))
end

"""Three rounds of `f` on `bytes`, each result compared as text with the first, a full
collection after each. The result is text at once, so no loop variable is a Union."""
function rounds(label, f, bytes)
    ref = ""
    for r in 1:3
        t = @elapsed (x = repr(f(bytes)))
        isempty(ref) && (ref = x)
        x == ref || error("$label round $r differs: $x vs $ref")
        println(label, " round ", r, " ", x, " ", round(t; digits = 1), " s"); flush(stdout)
        GC.gc(true)
        println(label, " full gc ", r, " ok"); flush(stdout)
    end
end
rounds("derivation", derive_round, POP)
rounds("identify", identify_round, SEEDS)
println("finished in ", round(Int, time() - t0), " s")
