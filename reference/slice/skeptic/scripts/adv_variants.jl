# ABOUTME: Rewrites each adversarial mechanism (step order, directions, group labels, list order)
# ABOUTME: and checks that the mechanism, its key and its compiled type do not change.
include(joinpath(@__DIR__, "advlib.jl"))
include(joinpath(@__DIR__, "adv_cases.jl"))

"""`p`'s raw input rewritten: steps permuted and randomly reversed, group labels permuted."""
function variant(p, rng)
    r = S.pop_raw(p)
    rf(f) = S.RawForm(f.conf, shuffle(rng, f.bound), shuffle(rng, f.added),
                      shuffle(rng, f.subtracted))
    n, G = length(r.steps), maximum(r.group)
    perm, lab = randperm(rng, n), randperm(rng, G)
    steps = map(perm) do j
        s = r.steps[j]
        f, t, c, x = rf(s.from), rf(s.to), shuffle(rng, s.consumed), shuffle(rng, s.released)
        rand(rng, Bool) ? S.RawStep(t, f, x, c, s.re) : S.RawStep(f, t, c, x, s.re)
    end
    tags = similar(r.tags)
    for g in eachindex(r.tags); tags[lab[g]] = r.tags[g]; end
    first(S.canonicalize(shuffle(rng, r.substrates), shuffle(rng, r.products), steps,
                         [lab[r.group[j]] for j in perm]; tags, r.multiplicity,
                         sites = shuffle(rng, r.sites)))
end

function main()
    rng = Xoshiro(23)
    for cs in CASES
        p = try cs.thunk() catch; continue end
        m = try S.Converted(p, 0).mech catch; continue end
        key = try S.canonical(S.rate_law(m)).key catch; continue end
        bad = 0
        for _ in 1:5
            v = variant(p, rng)
            bad += !(v == m && S.compile(v) === S.compile(m) && S.canonical(S.rate_law(v)).key == key)
        end
        @printf("%-30s 5 rewrites, %d differ\n", cs.name, bad)
    end
end

main()
