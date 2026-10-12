# ABOUTME: Solve-stage check (c): identical solves under rewritten inputs (order, direction,
# ABOUTME: group labels), row-space-preserving row changes, and same-model regroupings.
include(joinpath(@__DIR__, "lib.jl"))

records = [r for r in load_records() if (try S.solve(r.c.mech); true catch; false end)]
@printf("records solved: %d\n", length(records))

signature(sol) = (sol.cols, sol.indep, [(d.col, d.terms, d.keq) for d in sol.deps], sol.fold)
labels(sol) = [d.label for d in sol.deps]

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

rank_of(m) = (s = S.constraint_system(m); rank_exact([s[2] s[3]]))

rng = Xoshiro(20261011)
nvar, badvar, badvarnames = 0, String[], 0
nrow, badrow = 0, String[]
split_tried, split_same, split_fail, split_label, split_names = 0, 0, String[], 0, 0
merge_tried, merge_same, merge_fail, merge_names = 0, 0, String[], 0
plainallo, name_examples = [0, 0], String[]
for r in records
    m = r.c.mech
    sol = S.solve(m)
    sig, fit = signature(sol), S.fitted_params(m)
    # (i) rewritten inputs
    raw = S.pop_raw(r.p)
    for _ in 1:5
        global nvar += 1
        v = variant(raw, rng)
        sv = S.solve(v)
        signature(sv) == sig && labels(sv) == labels(sol) ||
            push!(badvar, "$(r.set) $(r.p.id)")
        S.fitted_params(v) == fit || (global badvarnames += 1)
    end
    # (ii) row-space-preserving row changes: shuffle, then add random row multiples
    cols, A, b = S.constraint_system(m)
    n = length(cols)
    M0 = Rational{BigInt}.([A b])
    p0 = S.rref!(M0, n)
    for _ in 1:3
        global nrow += 1
        B = [A b][shuffle(rng, 1:size(A, 1)), :]
        for _ in 1:size(B, 1)
            i, j = rand(rng, axes(B, 1)), rand(rng, axes(B, 1))
            i == j || (B[i, :] .+= rand(rng, -2:2) .* B[j, :])
        end
        M1 = Rational{BigInt}.(B)
        p1 = S.rref!(M1, n)
        p1 == p0 && M1[1:length(p1), :] == M0[1:length(p0), :] ||
            push!(badrow, "$(r.set) $(r.p.id)")
    end
    # (iii) same-model regroupings: splits of each group into singletons and into a random
    # bipartition; merges along folds between groups with one tag
    base = rank_of(m)
    G = S.ngroups(m)
    for g in 1:G
        mem = S.members(m, g)
        length(mem) < 2 && continue
        parts = [[[k] for k in mem[2:end]],
                 (q = shuffle(rng, mem); [q[1:rand(rng, 1:length(q) - 1)]])]
        for moved in parts
            global split_tried += 1
            group, tags = Int.(m.group), copy(m.tags)
            for (i, ks) in enumerate(moved)
                foreach(k -> group[k] = G + i, ks)
                S.is_allosteric(m) && push!(tags, m.tags[g])
            end
            m2 = regrouped(m, group, tags)
            rank_of(m2) == base || continue
            global split_same += 1
            s2 = S.solve(m2)
            signature(s2) == sig || push!(split_fail, "$(r.set) $(r.p.id) group $g")
            labels(s2) == labels(sol) || (global split_label += 1)
            plainallo[S.is_allosteric(m) ? 2 : 1] += 1
            f2 = S.fitted_params(m2)
            f2 == fit && continue
            global split_names += 1
            length(name_examples) < 6 && push!(name_examples,
                "$(r.set) $(r.p.id) group $g: $(setdiff(fit, f2)) -> $(setdiff(f2, fit))")
        end
    end
    for (d, i) in sol.fold
        cd, ci = sol.cols[d], sol.cols[i]
        g, h = m.group[cd.edge], m.group[ci.edge]
        (g == h || cd.inactive != ci.inactive || cd.slot != ci.slot) && continue
        S.is_allosteric(m) && m.tags[g] !== m.tags[h] && continue
        global merge_tried += 1
        group = [x == h ? Int(g) : Int(x) for x in m.group]
        group = [x > h ? x - 1 : x for x in group]
        tags = S.is_allosteric(m) ? deleteat!(copy(m.tags), h) : Symbol[]
        m2 = try
            regrouped(m, group, tags)
        catch
            continue
        end
        rank_of(m2) == base || continue
        global merge_same += 1
        signature(S.solve(m2)) == sig || push!(merge_fail, "$(r.set) $(r.p.id)")
        S.fitted_params(m2) == fit || (global merge_names += 1)
    end
end
@printf("(i) rewritten inputs: %d, solve or labels differ %d, fitted names differ %d\n",
        nvar, length(badvar), badvarnames)
@printf("(ii) row-space-preserving row changes: %d, RREF differs %d\n", nrow, length(badrow))
@printf("(iii) splits: %d tried, %d keep the row space; solve differs %d; labels change \
        in %d; fitted names change in %d\n", split_tried, split_same, length(split_fail),
        split_label, split_names)
@printf("     merges along folds: %d built, %d keep the row space; solve differs %d; \
        fitted names change in %d\n", merge_tried, merge_same, length(merge_fail),
        merge_names)
@printf("     same-model splits: %d of plain mechanisms, %d of allosteric ones\n",
        plainallo...)
println("     fitted-name changes under same-model splits, examples:")
foreach(x -> println("       ", x), name_examples)
foreach(x -> println("    ", x), first([badvar; badrow; split_fail; merge_fail], 20))
