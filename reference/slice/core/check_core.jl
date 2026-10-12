# ABOUTME: Checks the core stage on every converted set: equality classes against the package,
# ABOUTME: equivalent variants, mutations, orientation, names, Sig round trip, boundary checks.
const PROTO = "/tmp/claude-501/-home-denis-linux--julia-dev-EnzymeRates/12f42995-f5ba-4fda-9510-dd9b87f7d069/scratchpad/proto"
const DATA = "/tmp/claude-501/-home-denis-linux--julia-dev-EnzymeRates/12f42995-f5ba-4fda-9510-dd9b87f7d069/scratchpad/slice/core/data"
include(joinpath(PROTO, "common", "Population.jl"))
include(joinpath(@__DIR__, "..", "src", "Slice.jl"))
include(joinpath(@__DIR__, "retag.jl"))
using Serialization, Printf, Random, Statistics
const S = Slice

manifest = [split(l, '\t') for l in readlines(joinpath(DATA, "MANIFEST.tsv"))[2:end]]
sets = [(name = String(r[1]), pops = deserialize(r[4]),
         conv = deserialize(joinpath(DATA, r[1] * ".jls"))) for r in manifest]
# Package-built records: the probes exported them through the package's own constructor.
const HANDBUILT = ("p1_direct", "p2_hand", "p2_fuzz_dump", "p2_fuzz_dump_flip")
pkgbuilt(set) = !(set.name in HANDBUILT)
records = [(set = s.name, p = s.name == "p5a_cases" ? retag(s.pops[c.index]) : s.pops[c.index],
            c = c) for s in sets for c in s.conv]
@printf("records: %d converted in %d sets\n", length(records), length(sets))

# ─── 1. Orientation against the package's stored direction ──────────────────────────────
tiers, bad_orient, nsteps = zeros(Int, 4), String[], 0
for r in records
    fwd, tier = S.orientation(r.c.mech)
    foreach(t -> tiers[t + 1] += 1, tier)
    pkgbuilt((name = r.set,)) || continue
    for j in eachindex(r.p.steps)
        global nsteps += 1
        fwd[r.c.edge_of[j]] == !r.c.flipped[j] ||
            push!(bad_orient, "$(r.set) $(r.p.id) step $j: $(Pop.step_text(r.p, r.p.steps[j]))")
    end
end
@printf("orientation: %d package steps, %d differ from the package's stored direction\n",
        nsteps, length(bad_orient))
@printf("  edges decided by tier 0/1/2/3 (all records): %s\n", tiers)
foreach(x -> println("    ", x), first(bad_orient, 10))

# ─── 2. Equality and hash classes against the package's structural key ──────────────────
pkgkey(m) = join([join(m.substrates, ","), join(m.products, ","),
                  string(m.allosteric, m.cat_multiplicity), join(m.group_tags, ","),
                  [string(Pop.step_text(m, s), " @", g)
                   for (s, g) in zip(m.steps, Pop.group_of_step(m))]...,
                  [string(s.ligands, s.states, s.multiplicity) for s in m.sites]...], " | ")
# Rejected P5a grid cases: the package's structure key cannot see their restored tags.
pk = [r for r in records if pkgbuilt((name = r.set,)) && (r.p.ok || r.set != "p5a_cases")]
bykey = Dict{String, Vector{Int}}()
bymech = Dict{S.Mechanism, Vector{Int}}()
for (i, r) in enumerate(pk)
    push!(get!(bykey, pkgkey(r.p), Int[]), i)
    push!(get!(bymech, r.c.mech, Int[]), i)
end
split_ = count(v -> length(unique(pk[i].c.mech for i in v)) > 1, values(bykey))
merge_ = count(v -> length(unique(pkgkey(pk[i].p) for i in v)) > 1, values(bymech))
hashes = unique(hash(m) for m in keys(bymech))
@printf("classes: %d package-built records, %d package classes, %d mechanism classes, \
        %d distinct hashes; false splits %d, false merges %d\n", length(pk), length(bykey),
        length(bymech), length(hashes), split_, merge_)

# ─── 3. Equivalent variants and the Sig ─────────────────────────────────────────────────
function variant(r, rng; flipall = false)
    sh(v) = rng === nothing ? copy(v) : shuffle(rng, v)
    rf(f) = S.RawForm(f.conf, sh(f.bound), sh(f.added), sh(f.subtracted))
    n, G = length(r.steps), maximum(r.group; init = 0)
    perm = rng === nothing ? collect(1:n) : randperm(rng, n)
    lab = rng === nothing ? collect(1:G) : randperm(rng, G)
    steps = map(perm) do j
        s = r.steps[j]
        f, t, c, x = rf(s.from), rf(s.to), sh(s.consumed), sh(s.released)
        flipall || (rng !== nothing && rand(rng, Bool)) ? S.RawStep(t, f, x, c, s.re) :
                                                          S.RawStep(f, t, c, x, s.re)
    end
    tags = similar(r.tags)
    for g in eachindex(r.tags); tags[lab[g]] = r.tags[g]; end
    sites = sh([(q = rng === nothing ? collect(eachindex(s.ligands)) :
                     randperm(rng, length(s.ligands));
                 S.Site(s.ligands[q], s.states[q], s.multiplicity)) for s in r.sites])
    S.canonicalize(sh(r.substrates), sh(r.products), steps, [lab[r.group[j]] for j in perm];
                   tags, r.multiplicity, sites)
end
rng = Xoshiro(20261010)
nvar, badvar, badsig, badtype = 0, String[], String[], 0
for r in records
    m = r.c.mech
    raw = S.pop_raw(r.p)
    for k in 1:6
        v = first(k == 1 ? variant(raw, nothing; flipall = true) : variant(raw, rng))
        global nvar += 1
        v == m && hash(v) == hash(m) || push!(badvar, "$(r.set) $(r.p.id) variant $k")
        k == 2 && (typeof(S.compile(v)) === typeof(S.compile(m)) || (global badtype += 1))
    end
    s = S.sig(m)
    d = S.decode(s...)
    d == m && hash(d) == hash(m) && S.sig(d) == s && S.mechanism(S.compile(m)) == m ||
        push!(badsig, "$(r.set) $(r.p.id)")
end
@printf("variants: %d built, %d unequal or with another hash; compiled types differing \
        between equal variants: %d\n", nvar, length(badvar), badtype)
foreach(x -> println("    ", x), first(badvar, 10))
siglen = [length(S.sig(r.c.mech)) for r in records]
@printf("Sig: %d of %d round trips fail; length median %d, max %d Ints\n", length(badsig),
        length(records), median(siglen), maximum(siglen))

# ─── 4. Mutations: never equal to the original, never sharing its hash ─────────────────
function mutate(raw, kind, rng)
    steps, group, tags = copy(raw.steps), copy(raw.group), copy(raw.tags)
    sites, mult, subs, prods = copy(raw.sites), raw.multiplicity, raw.substrates, raw.products
    G = maximum(group)
    if kind === :re_flip
        g = rand(rng, 1:G)
        steps = [group[j] == g ? S.RawStep(s.from, s.to, s.consumed, s.released, !s.re) : s
                 for (j, s) in enumerate(steps)]
    elseif kind === :split
        big = [g for g in 1:G if count(==(g), group) > 1]
        isempty(big) && return nothing
        g = rand(rng, big)
        group[rand(rng, findall(==(g), group))] = G + 1
        isempty(tags) || push!(tags, tags[g])
    elseif kind === :merge
        G > 1 || return nothing
        g, h = minmax(randperm(rng, G)[1:2]...)
        group = [x == h ? g : x > h ? x - 1 : x for x in group]
        isempty(tags) || deleteat!(tags, h)
    elseif kind === :tag
        isempty(tags) && return nothing
        g = rand(rng, 1:G)
        tags[g] = rand(rng, filter(!=(tags[g]), collect(S.TAGS)))
    elseif kind === :multiplicity
        isempty(tags) && return nothing
        mult += 1
    elseif kind === :site
        isempty(sites) && return nothing
        i = rand(rng, eachindex(sites))
        s = sites[i]
        st = copy(s.states)
        st[1] = rand(rng, filter(!=(st[1]), collect(S.SITE_STATES)))
        sites[i] = S.Site(s.ligands, st, s.multiplicity)
    elseif kind === :modifier
        j = rand(rng, eachindex(steps))
        s = steps[j]
        flip(x) = S.Met(x.name, !x.modifier)
        isempty(s.consumed) && return nothing
        steps[j] = S.RawStep(s.from, s.to, [flip(s.consumed[1]); s.consumed[2:end]],
                             s.released, s.re)
    elseif kind === :conf
        c = rand(rng, unique(f.conf for s in steps for f in (s.from, s.to)))
        rn(f) = S.RawForm(f.conf === c ? :Zq : f.conf, f.bound, f.added, f.subtracted)
        steps = [S.RawStep(rn(s.from), rn(s.to), s.consumed, s.released, s.re) for s in steps]
    elseif kind === :swap
        subs, prods = copy(subs), copy(prods)
        i, k = rand(rng, eachindex(subs)), rand(rng, eachindex(prods))
        subs[i], prods[k] = prods[k], subs[i]
    elseif kind === :drop
        length(steps) > 1 || return nothing
        j = rand(rng, eachindex(steps))
        deleteat!(steps, j)
        g = group[j]
        deleteat!(group, j)
        if !(g in group)
            group = [x > g ? x - 1 : x for x in group]
            isempty(tags) || deleteat!(tags, g)
        end
    end
    S.canonicalize(subs, prods, steps, group; tags, multiplicity = mult, sites)
end
const KINDS = (:re_flip, :split, :merge, :tag, :multiplicity, :site, :modifier, :conf,
               :swap, :drop)
mstats = Dict(k => [0, 0, 0, 0] for k in KINDS)    # tried, rejected, differ, equal or same hash
for r in records, kind in KINDS
    raw = S.pop_raw(r.p)
    st = mstats[kind]
    res = try
        mutate(raw, kind, rng)
    catch
        st[1] += 1; st[2] += 1; continue
    end
    res === nothing && continue
    st[1] += 1
    m2 = first(res)
    m2 != r.c.mech && hash(m2) != hash(r.c.mech) ? (st[3] += 1) : (st[4] += 1)
end
println("mutations (tried, rejected by the constructor, differ, equal or same hash):")
for k in KINDS; @printf("    %-13s %s\n", k, mstats[k]); end

# ─── 5. Names against the package's (decided changes applied to its names) ──────────────
"""The package's name with the decided changes applied, and which rules changed it."""
function transform(name::Symbol, m::S.Mechanism)
    s0 = String(name)
    s = replace(s0, "_+" => "_plus", "_-" => "_minus")
    rules = s == s0 ? String[] : ["residual"]
    for x in S.inhibitors(m)
        dual = S.is_reactant(m, x)
        t = replace(s, string(x, "inh") => dual ? string(x, "mod") : string(x))
        t == s || push!(rules, dual ? "dual-role modifier" : "modifier-only")
        s = t
    end
    Symbol(s), unique!(rules)
end
nnames, changed, mism, untied = 0, Dict{String, Int}(), String[], 0
for r in records
    p, c, m = r.p, r.c, r.c.mech
    allo = S.is_allosteric(m)
    (allo ? !isempty(p.step_names_A) : !isempty(p.step_names)) || continue
    p.ok || continue                # rejected P5a cases: names come from a NonequalAI twin
    states = allo ? ((false, p.step_names_A), (true, p.step_names_I)) :
                    ((false, p.step_names),)
    for (inactive, pkg) in states
        mine = S.constant_names(m; inactive)
        for j in eachindex(p.steps), (a, b) in zip(mine[c.edge_of[j]], pkg[j])
            global nnames += 1
            t, rules = transform(b, m)
            foreach(k -> changed[k] = get(changed, k, 0) + 1, rules)
            a == t && continue
            g = m.group[c.edge_of[j]]
            if inactive && m.tags[g] === :OnlyA && length(S.members(m, g)) > 1
                global untied += 1
            else
                push!(mism, "$(r.set) $(p.id) step $j $(inactive ? "I" : "A"): mine $a, \
                             package $b")
            end
        end
    end
end
@printf("names: %d compared (per step and state); changed by a decided rule: %s; untied \
         inactive :OnlyA members %d; unexplained %d\n", nnames, sort!(collect(changed)),
        untied, length(mism))
foreach(x -> println("    ", x), first(mism, 20))
regc, regbad, regbare = 0, String[], 0
for r in records
    p, m = r.p, r.c.mech
    S.is_allosteric(m) && p.ok || continue
    mine = Set(n for (_, ns) in S.regulator_names(m) for n in ns)
    for f in p.fitted_params
        occursin(r"reg$", String(f)) || continue
        global regc += 1
        lig = Symbol(replace(String(f), r"^K_(A_|I_)?" => "", r"reg$" => ""))
        st = only(st for s in m.sites for (l, st) in zip(s.ligands, s.states) if l === lig)
        want = st === :EqualAI ? Symbol("K_", lig, "reg") : f
        want == f || (global regbare += 1)
        want in mine || push!(regbad, "$(p.id) $f")
    end
end
@printf("regulator names: %d fitted package names, %d renamed bare (:EqualAI), %d without \
         a match\n", regc, regbare, length(regbad))
foreach(x -> println("    ", x), first(regbad, 10))

# ─── 6. Derived lists against the package ───────────────────────────────────────────────
metorder, metset, inhset, regset = 0, 0, 0, 0
for r in records
    p, m = r.p, r.c.mech
    pkgbuilt((name = r.set,)) || continue
    Set(S.metabolites(m)) == Set(p.metabolites) || (global metset += 1)
    S.metabolites(m) == p.metabolites || (global metorder += 1)
    Set(S.inhibitors(m)) == Set(p.inhibitors) || (global inhset += 1)
    Set(S.regulators(m)) == Set(p.regulators) || (global regset += 1)
end
@printf("derived lists vs the package: metabolites differ as sets %d, in order only %d; \
        inhibitors differ %d; regulators differ %d\n", metset, metorder - metset, inhset,
        regset)

# ─── 7. Boundary checks: every package-built record passes ──────────────────────────────
bfail = String[]
for r in records
    pkgbuilt((name = r.set,)) || continue
    try
        S.check_boundary(r.c.mech)
    catch e
        push!(bfail, "$(r.set) $(r.p.id): " * first(sprint(showerror, e), 160))
    end
end
@printf("boundary checks on package-built records: %d fail\n", length(bfail))
foreach(x -> println("    ", x), first(bfail, 10))
const ATOMS = Dict(
    "uniuni" => Dict(:S => [:C => 1], :P => [:C => 1]),
    "unibi" => Dict(:S => [:A => 1, :B => 1], :P => [:A => 1], :Q => [:B => 1]),
    "bibi" => Dict(:A => [:C => 1], :B => [:N => 1], :P => [:C => 1], :Q => [:N => 1]),
    "bibi_pp" => Dict(:A => [:C => 1, :X => 1], :B => [:N => 1], :P => [:C => 1],
                      :Q => [:N => 1, :X => 1]),
    "terbi" => Dict(:A => [:C => 1], :B => [:N => 1], :D => [:X => 1],
                    :P => [:C => 1, :N => 1], :Q => [:X => 1]))
const RXN = Dict("init_uniuni" => "uniuni", "allo_uniuni" => "uniuni",
                 "allo_uniuni_inh" => "uniuni", "init_unibi" => "unibi",
                 "allo_unibi" => "unibi", "init_bibi" => "bibi", "allo_bibi" => "bibi",
                 "init_bibi_pp" => "bibi_pp", "expand_bibi" => "bibi_pp",
                 "allo_bibi_pp" => "bibi_pp", "allo_terbi" => "terbi")
natoms, afail = 0, String[]
for r in records
    rx = r.set == "p1_enum_inh" ? (startswith(r.p.id, "uni") ? "uniuni" : "bibi") :
         get(RXN, r.set, nothing)
    rx === nothing && continue
    global natoms += 1
    try
        S.check_atoms(r.c.mech, ATOMS[rx])
    catch e
        push!(afail, "$(r.set) $(r.p.id): " * first(sprint(showerror, e), 160))
    end
end
@printf("atom conservation on %d enumerated records: %d fail\n", natoms, length(afail))
foreach(x -> println("    ", x), first(afail, 10))

# ─── 7b. Rejections: each structural and boundary error fires on a hand-built input ───────
F(conf, ms...; add = Symbol[], sub = Symbol[]) = S.RawForm(conf, collect(S.Met, ms), add, sub)
X(x, mod = false) = S.Met(x, mod)
st(f, t, c, r, re) = S.RawStep(f, t, collect(S.Met, c), collect(S.Met, r), re)
uni = [st(F(:E), F(:E, X(:S)), [X(:S)], [], true), st(F(:E, X(:S)), F(:E, X(:P)), [], [], false),
       st(F(:E, X(:P)), F(:E), [], [X(:P)], true)]
ctor(steps; subs = [:S], prods = [:P], kw...) =
    first(S.canonicalize(subs, prods, steps, collect(eachindex(steps)); kw...))
bottom = [st(F(:E), F(:E, X(:S)), [X(:S)], [], false),
          st(F(:E, X(:S)), F(:E, X(:P), X(:Q)), [], [], false),
          st(F(:E, X(:P), X(:Q)), F(:E, X(:P)), [], [X(:Q)], true),
          st(F(:E, X(:P), X(:Q)), F(:E, X(:Q)), [], [X(:P)], true),
          st(F(:E, X(:P)), F(:E), [], [X(:P)], false), st(F(:E, X(:Q)), F(:E), [], [X(:Q)], false)]
nadp = [st(F(:E), F(:E, X(:NAD)), [X(:NAD)], [], true),
        st(F(:E, X(:NAD)), F(:E, X(:NAD), X(:P)), [X(:P)], [], true),
        st(F(:E), F(:E, X(:NADP)), [X(:NADP)], [], true),
        st(F(:E, X(:NAD), X(:P)), F(:E, X(:NADP)), [], [], false)]
cases = [
    ("33 metabolite entries", "at most 32 metabolite entries", () -> ctor(
        [st(F(:E), F(:E, X(Symbol("I", i), true)), [X(Symbol("I", i), true)], [], true)
         for i in 1:33]; subs = [:S], prods = [:P])),
    ("16 copies on a side", "more than 15 times", () -> ctor(
        [st(F(:E), F(:E, fill(X(:S), 16)...), fill(X(:S), 16), [], true)])),
    ("non-reactant without the modifier flag", "binds only as a modifier",
     () -> ctor([uni; st(F(:E), F(:E, X(:I)), [X(:I)], [], true)])),
    ("residual subtracting a substrate", "a residual adds substrates",
     () -> ctor([uni; st(F(:E), F(:E; sub = [:S]), [], [], false)])),
    ("self-loop", "joins a form to itself", () -> ctor([uni; st(F(:E), F(:E), [], [], false)])),
    ("one metabolite on both sides", "takes up and gives off one metabolite",
     () -> ctor([uni; st(F(:E), F(:F), [X(:S)], [X(:S)], false)])),
    ("one reaction twice", "one reaction", () -> ctor([uni; uni[2]])),
    ("disconnected steps", "one connected graph",
     () -> ctor([uni; st(F(:F), F(:F, X(:S)), [X(:S)], [], true)])),
    ("sites on a plain mechanism", "a plain mechanism",
     () -> ctor(uni; sites = [S.Site([:R], [:OnlyA], 1)])),
    ("tags of the wrong count", "tags each kinetic group", () -> ctor(uni; tags = [:OnlyA])),
    ("group of two kinds", "differ in kind", () -> first(S.canonicalize([:S], [:P],
        [uni; st(F(:E), F(:F), [], [], false)], [1, 2, 3, 1]))),
    ("ligand at two sites", "one regulatory site", () -> ctor(uni;
        tags = [:EqualAI, :EqualAI, :EqualAI],
        sites = [S.Site([:R], [:OnlyA], 1), S.Site([:R, :T], [:OnlyI, :OnlyA], 2)])),
    ("bottomless RE segment", "has no form free of one side", () -> S.check_boundary(
        ctor(bottom; prods = [:P, :Q]))),
    ("role change E + P::Mod → E(P)", "changes a metabolite's role", () -> S.check_boundary(
        ctor([uni; st(F(:E), F(:E, X(:P)), [X(:P, true)], [], true)]))),
    ("colliding form names E(NAD, P), E(NADP)", "both render the name ENADP",
     () -> S.check_boundary(ctor(nadp; subs = [:NAD, :NADP], prods = [:P]))),
    ("atoms not conserved", "does not conserve atoms",
     () -> S.check_atoms(ctor(uni), Dict(:S => [:C => 1], :P => [:N => 1]))),
]
negfail = 0
for (what, msg, f) in cases
    err = try
        f(); "no error"
    catch e
        sprint(showerror, e)
    end
    ok = occursin(msg, err)
    global negfail += !ok
    @printf("    rejection %-42s %s\n", what, ok ? "fires" : "MISSING: " * first(err, 120))
end
@printf("rejections: %d of %d fire as expected\n", length(cases) - negfail, length(cases))

# ─── 8. Costs ───────────────────────────────────────────────────────────────────────────
function time_ns_per(f, k; reps = 11)
    f()
    ts = map(1:reps) do _
        t0 = time_ns()
        for _ in 1:k; Base.donotdelete(f()); end
        (time_ns() - t0) / k
    end
    median(ts)
end
costs = Dict(k => Float64[] for k in (:ctor, :orient, :names, :sig, :decode, :hash, :eq))
for r in records[1:3:end]
    raw, m = S.pop_raw(r.p), r.c.mech
    s = S.sig(m)
    push!(costs[:ctor], time_ns_per(() -> S.canonicalize(raw), 50))
    push!(costs[:orient], time_ns_per(() -> S.orient(m), 50))
    push!(costs[:names], time_ns_per(() -> S.constant_names(m), 20))
    push!(costs[:sig], time_ns_per(() -> S.sig(m), 50))
    push!(costs[:decode], time_ns_per(() -> S.decode(s...), 50))
    push!(costs[:hash], time_ns_per(() -> hash(m), 200))
    m2 = deepcopy(m)
    push!(costs[:eq], time_ns_per(() -> m == m2, 200))
end
for (k, v) in sort!(collect(costs); by = first)
    @printf("cost %-7s median %8.0f ns, p90 %8.0f ns, max %8.0f ns (%d mechanisms)\n", k,
            median(v), quantile(v, 0.9), maximum(v), length(v))
end
