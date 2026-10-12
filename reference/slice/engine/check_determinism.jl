# ABOUTME: Engine-stage check (c): bare and final polynomials, printed laws and kcat candidates
# ABOUTME: under rewritten inputs (step order, directions, group labels) and across processes.
include(joinpath(@__DIR__, "lib.jl"))

"""Everything the engine stage derives from `m`, as comparable values."""
function derived(m)
    b = S.bare_law(m)
    sol = S.solve(m)
    L = S.law(m, sol)
    k = try
        S.kcat_terms(m, sol)[2]
    catch err
        sprint(showerror, err)
    end
    (b.num, b.den, b.sigma, b.trees, b.cols, L.vars.names, L.H, L.Q, L.D, L.sigma, L.trees,
     L.shift, S.rate_equation_string(m, L, sol), S.coefficient_form(m, L), k, S.complexity(m))
end

"""`raw` with its steps shuffled, each reversed with probability 1/2, and its group labels
permuted."""
function rewrite(raw, rng)
    n, G = length(raw.steps), maximum(raw.group)
    perm, relabel = randperm(rng, n), randperm(rng, G)
    flip(s) = rand(rng, Bool) ? S.RawStep(s.to, s.from, s.released, s.consumed, s.re) : s
    (; raw.substrates, raw.products, steps = [flip(raw.steps[j]) for j in perm],
     group = [relabel[raw.group[j]] for j in perm], tags = raw.tags, raw.multiplicity, raw.sites)
end

records = [r for r in load_records() if !S.is_allosteric(r.c.mech)]
path = joinpath(@__DIR__, "data", "derived.jls")
first_run = !isfile(path)
rng = Xoshiro(7)
mine = Dict{Tuple{String, Int}, Any}()
nvar, neq, bad, t0 = 0, 0, String[], time()
for r in records
    m = r.c.mech
    d = try
        derived(m)
    catch err
        sprint(showerror, err)
    end
    mine[(r.set, r.c.index)] = d
    first_run || continue
    for _ in 1:3
        m2 = first(S.canonicalize(rewrite(S.pop_raw(r.p), rng)))
        global nvar += 1
        m2 == m || (push!(bad, "$(r.set) $(r.p.id): rewritten mechanism differs"); continue)
        d2 = try
            derived(m2)
        catch err
            sprint(showerror, err)
        end
        d2 == d ? (global neq += 1) : push!(bad, "$(r.set) $(r.p.id): derived values differ")
    end
end
if first_run
    serialize(path, mine)
    @printf("(c) %d plain records, %d rewritten inputs (step order, direction, group labels): \
            derived values identical in %d; differences %d (%.0f s)\n", length(records), nvar,
            neq, length(bad), time() - t0)
    foreach(println, first(bad, 20))
else
    old = deserialize(path)
    same = count(k -> old[k] == mine[k], keys(mine))
    @printf("(c) across processes: %d of %d records derive identical values (%.0f s)\n", same,
            length(mine), time() - t0)
end
