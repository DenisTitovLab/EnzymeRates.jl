# ABOUTME: K1 key invariance per record: rewritten inputs, same-model regroupings, metabolite
# ABOUTME: renamings, and a second process; usage: check_keys.jl label sets… (needs pipe rows)
include(joinpath(@__DIR__, "lib.jl"))

const LABEL, SETS = ARGS[1], ARGS[2:end]

"""The core's rewritten input: steps reordered and randomly reversed, group labels permuted,
every list shuffled."""
function variant(r, rng)
    rf(f) = S.RawForm(f.conf, shuffle(rng, f.bound), shuffle(rng, f.added),
                      shuffle(rng, f.subtracted))
    n, G = length(r.steps), maximum(r.group; init = 0)
    perm, lab = randperm(rng, n), randperm(rng, G)
    steps = map(perm) do j
        s = r.steps[j]
        f, t, c, x = rf(s.from), rf(s.to), shuffle(rng, s.consumed), shuffle(rng, s.released)
        rand(rng, Bool) ? S.RawStep(t, f, x, c, s.re) : S.RawStep(f, t, c, x, s.re)
    end
    tags = similar(r.tags)
    for g in eachindex(r.tags); tags[lab[g]] = r.tags[g]; end
    sites = shuffle(rng, [(q = randperm(rng, length(s.ligands));
                           S.Site(s.ligands[q], s.states[q], s.multiplicity)) for s in r.sites])
    first(S.canonicalize(shuffle(rng, r.substrates), shuffle(rng, r.products), steps,
                         [lab[r.group[j]] for j in perm]; tags, r.multiplicity, sites))
end

"""`m` with its steps regrouped by `group` (one label per edge), tags following `tags`."""
function regrouped(m, group, tags)
    (subs, prods, steps, _), kw = S.raw(m)
    first(S.canonicalize(subs, prods, steps, group; tags, kw.multiplicity, kw.sites))
end

"""`m` with every metabolite name `x` renamed `f(x)` (forms, steps, roles and sites)."""
function renamed(m, f)
    (subs, prods, steps, group), kw = S.raw(m)
    met(x) = S.Met(f(x.name), x.modifier)
    form(q) = S.RawForm(q.conf, met.(q.bound), f.(q.added), f.(q.subtracted))
    sites = [S.Site(f.(s.ligands), s.states, s.multiplicity) for s in kw.sites]
    S.Mechanism(f.(subs), f.(prods), [S.RawStep(form(s.from), form(s.to), met.(s.consumed),
                                                met.(s.released), s.re) for s in steps],
                group; kw.tags, kw.multiplicity, sites)
end

"""A renaming of `m`'s metabolites that reverses the sorted order of the substrates, of the
products and of the other names, so every position changes where a role has two names."""
function reverser(m)
    out = Dict{Symbol, Symbol}()
    for (tag, xs) in (("s", m.substrates), ("p", m.products),
                      ("o", setdiff(S.metabolites(m), [m.substrates; m.products])))
        foreach(((i, x),) -> out[x] = Symbol(tag, length(xs) + 1 - i), enumerate(xs))
    end
    x -> out[x]
end

rank_of(m) = (s = S.constraint_system(m); rank_exact([s[2] s[3]]))
key_of(m) = S.canonical(S.rate_law(m)).key
solve_sig(sol) = (sol.indep, [(d.col, d.terms, d.keq) for d in sol.deps])

rng = Xoshiro(hash(LABEL) % UInt32)
out = NamedTuple[]
t0 = time()
pipe = Dict((x.set, x.index) => x for l in readdir(IDATA) if startswith(l, "pipe_")
            for x in deserialize(joinpath(IDATA, l)) if x.set in SETS && x.ok)
for set in SETS
    for r in records_of(set)
        haskey(pipe, (set, r.c.index)) || continue
        # init_terter: every 50th record (the compile checks recompute every chosen key)
        startswith(set, "init_terter") && r.c.index % 50 != 0 && continue
        m, key = r.c.mech, pipe[(set, r.c.index)].key
        sol, fit = S.solve(m), S.fitted_params(m)
        second = key_of(m) == key                              # another process's key
        raw = S.pop_raw(r.p)
        nvar, badvar = 0, 0
        for _ in 1:2
            v = variant(raw, rng)
            nvar += 1
            badvar += !(v == m && key_of(v) == key)
        end
        # same-model regroupings: each group of two or more split into singletons (P2's
        # single-group splits) and into a random bipartition; same model when the rank of the
        # constraint system is unchanged (for a plain mechanism: the independent count)
        base = S.is_allosteric(m) ? rank_of(m) : length(sol.indep)
        same_model(m2) = base == (S.is_allosteric(m) ? rank_of(m2) : length(S.solve(m2).indep))
        nsplit, nsame, keyfail, namechange, solvefail = 0, 0, 0, 0, 0
        p2, p2same, p2key, p2names = 0, 0, 0, 0
        G = S.ngroups(m)
        for g in 1:G
            mem = S.members(m, g)
            length(mem) < 2 && continue
            for (kind, moved) in ((:singletons, [[k] for k in mem[2:end]]),
                                  (:bipartition, (q = shuffle(rng, mem);
                                                  [q[1:rand(rng, 1:length(q) - 1)]])))
                group, tags = Int.(m.group), copy(m.tags)
                for (i, ks) in enumerate(moved)
                    foreach(k -> group[k] = G + i, ks)
                    S.is_allosteric(m) && push!(tags, m.tags[g])
                end
                m2 = regrouped(m, group, tags)
                nsplit += 1
                p2 += kind === :singletons
                same_model(m2) || continue
                nsame += 1
                k2, f2 = key_of(m2), S.fitted_params(m2)
                keyfail += k2 != key
                namechange += f2 != fit
                solvefail += solve_sig(S.solve(m2)) != solve_sig(sol)
                if kind === :singletons && length(S.solve(m2).indep) == length(sol.indep)
                    p2same += 1
                    p2key += k2 != key
                    p2names += Set(f2) != Set(fit)
                end
            end
        end
        same_order = key_of(renamed(m, x -> Symbol("z", x))) == key
        rev = renamed(m, reverser(m))
        rev_key, rev_solve = key_of(rev) == key, solve_sig(S.solve(rev)) == solve_sig(sol)
        push!(out, (set = r.set, index = r.c.index, id = r.p.id, allo = S.is_allosteric(m),
                    second, nvar, badvar, nsplit, nsame, keyfail, namechange, solvefail,
                    p2, p2same, p2key, p2names, same_order, rev_key, rev_solve))
    end
    @printf("%-20s done at %.0f s\n", set, time() - t0)
    flush(stdout)
end
serialize(joinpath(IDATA, "keys_$LABEL.jls"), out)
s(f) = sum(getfield(x, f) for x in out; init = 0)
@printf("%s: %d records; second process differs %d; variants %d, differ %d; same-model splits \
        %d of %d, key changes %d, fitted-name changes %d, solve changes %d; P2 splits %d \
        (same model %d): key changes %d, fitted-name-set changes %d; order-preserving renaming \
        changes the key in %d; order-reversing renaming changes it in %d (solve in %d); \
        %.0f s\n", LABEL, length(out), count(x -> !x.second, out), s(:nvar), s(:badvar),
        s(:nsame), s(:nsplit), s(:keyfail), s(:namechange), s(:solvefail), s(:p2), s(:p2same),
        s(:p2key), s(:p2names), count(x -> !x.same_order, out), count(x -> !x.rev_key, out),
        count(x -> !x.rev_solve, out), time() - t0)
