# ABOUTME: Structural :OnlyA island classification (follow-up B1): components of the inactive step graph
# ABOUTME: without blocked steps, affinity/kinetic doors per reading, mass orders, classes; stdlib only.
const S = "/tmp/claude-501/-home-denis-linux--julia-dev-EnzymeRates/12f42995-f5ba-4fda-9510-dd9b87f7d069/scratchpad"
isdefined(Main, :Pop) || include(joinpath(S, "proto/common/Population.jl"))
using Serialization, Printf

# ─── Populations and P5a records ───────────────────────────────────────────

const POPS = ["allo_bibi", "expand_bibi", "allo_uniuni", "specs", "handwritten",
              "allo_uniuni_inh", "allo_unibi", "allo_bibi_pp", "allo_terbi"]
const SKEPTIC = Set(["allo_uniuni_inh", "allo_unibi", "allo_bibi_pp", "allo_terbi"])

"""Allosteric mechanisms of population `name` (common data or the skeptic's exports)."""
load_allo(name) = filter(m -> m.allosteric,
    name in SKEPTIC ? deserialize(joinpath(S, "proto/skeptic_P5a-onlyA/data/$name.jls"))::Vector{Pop.PMech} :
                      Pop.load_population(name))

"""P5a / skeptic limit-study records relevant to population `name`, keyed by (id, conv)."""
function load_records(name)
    files = name in SKEPTIC ? ["proto/skeptic_P5a-onlyA/data/res_$name.jls"] :
            name in ("allo_bibi", "expand_bibi") ? ["proto/p5a/data/pop_$name.jls"] :
            name == "allo_uniuni" ? ["proto/p5a/data/results_allo_uniuni.jls"] :
            ["proto/p5a/data/results_cases.jls"]
    d = Dict{Tuple{String, Symbol}, Dict{Symbol, Any}}()
    for f in files, r in deserialize(joinpath(S, f))
        startswith(r[:id], name * "/") && (d[(r[:id], r[:conv])] = r)
    end
    d
end

"""P5a's numeric verdict of one record: U=D, U≠D, NU (+seq when some direction sequesters),
empty (no feasible direction), err."""
function verdict(r)
    r === nothing && return "—"
    haskey(r, :error) && return "err"
    haskey(r, :cone_v) || return "empty"
    s = r[:cone_seq] ? "+seq" : ""
    r[:cone_spread] >= 1e-20 ? "NU" * s : r[:cone_vs_del] >= 1e-20 ? "U≠D" * s : "U=D" * s
end

# ─── Step classes and door types ───────────────────────────────────────────

mkey(x::Pop.PMet) = (x.name, x.role)

""":plain (one metabolite consumed, none released, `to` = `from` plus it, same residual and
conformation), :fused (any other binding), :iso (nothing taken up or given off), :tc (one taken
up, one given off), :other."""
function step_class(m::Pop.PMech, s::Pop.PStep)
    nc, nr = length(s.consumed), length(s.released)
    nc == 0 && nr == 0 && return :iso
    if nc == 1 && nr == 0
        a, b = m.species[s.from], m.species[s.to]
        plain = a.conformation === b.conformation &&
                sort(mkey.([a.bound; s.consumed])) == sort(mkey.(b.bound)) &&
                sort(a.res_added) == sort(b.res_added) && sort(a.res_subtracted) == sort(b.res_subtracted)
        return plain ? :plain : :fused
    end
    nc == 1 && nr == 1 && return :tc
    :other
end

"""Whether each step is an affinity door (true) or a kinetic door (false) when blocked, under
`reading`: `:D` (an RE plain binding is affinity, everything else kinetic), `:iu` (every plain
binding, RE or SS, is affinity), and the group-tied variants P5a used, `:ii` (as `:D`, but a group
holding any non-plain step is kinetic as a whole) and `:i` (as `:iu`, same group rule)."""
function affinity_flags(m::Pop.PMech, reading::Symbol)
    cls = [step_class(m, s) for s in m.steps]
    plain = cls .=== :plain
    if reading in (:ii, :i)
        g = Pop.group_of_step(m)
        chem = [any(!plain[i] for i in m.groups[g[j]]) for j in eachindex(m.steps)]
        reading === :ii && return [plain[j] && !chem[j] && m.steps[j].re for j in eachindex(m.steps)]
        return [plain[j] && !chem[j] for j in eachindex(m.steps)]
    end
    reading === :D && return [plain[j] && m.steps[j].re for j in eachindex(m.steps)]
    reading === :iu && return plain
    error("unknown reading $reading")
end

blocked_flags(m::Pop.PMech) = (g = Pop.group_of_step(m);
    [!isempty(m.group_tags) && m.group_tags[g[j]] === :OnlyA for j in eachindex(m.steps)])

# ─── Components, doors, orders ─────────────────────────────────────────────

"""Union-find labels 1..nc of the forms over the steps with `keep[j]`."""
function components(n::Int, steps, keep)
    p = collect(1:n)
    root(x) = p[x] == x ? x : (p[x] = root(p[x]))
    for (j, s) in enumerate(steps)
        keep[j] && (p[root(s.from)] = root(s.to))
    end
    r = [root(i) for i in 1:n]
    lab = Dict(v => k for (k, v) in enumerate(unique(r)))
    [lab[x] for x in r], length(lab)
end

"""Sign algebra of mass orders relative to main: 0 finite, 1 vanishing, -1 diverging,
2 ambiguous (depends on unequal exponents), 3 undefined (depends on vanishing rates)."""
function addsign(a, b)
    (a == 3 || b == 3) && return 3
    (a == 2 || b == 2) && return 2
    a == 0 && return b
    b == 0 && return a
    a == b ? a : 2
end

"""Structural analysis of one allosteric mechanism under `reading`.

Components: forms joined by unblocked steps; `main` holds the flip form. Doors: blocked steps
joining two components (affinity: stored `from` = unbound side U, `to` = bound side B).
Mass orders relative to main (exponent of ε, every affinity door's K_I ∝ 1/ε):
- inside an affinity cluster (components joined by affinity doors) orders follow RE equilibrium:
  B = U + 1 (unit potentials `φ`); a component below main (directed U→B path from main) has
  order > 0 (sign 1), above main (path to main) < 0 (sign -1), neither: sign 2 (ambiguous: the
  order is a difference of independent exponents), both: a directed cycle (infeasible);
- clusters hang off each other through kinetic doors (the cluster graph). A cluster reached through
  a bridge of that graph sits at equilibrium through the door's kept ratio: its door form takes the
  parent door form's sign. A cluster on a cycle of the cluster graph has sign 3 (undefined: its mass
  depends on how fast each door vanishes);
- a kinetic door inside one cluster ("self-loop") whose sides differ in φ cannot keep its ratio:
  thermodynamics (the cycle through the affinity doors) forces it (`forced`).
Class precedence: free_forms, infeasible (affinity chord or directed affinity cycle), sequestering
(an island of sign -1), undefined (3), ambiguous (2), phantom (0), emptied (all islands 1), no_island."""
function analyze(m::Pop.PMech, reading::Symbol)
    blocked = blocked_flags(m)
    aff = affinity_flags(m, reading)
    free = findall(Pop.is_free, m.species)
    nfree = length(free)
    comp, nc = components(length(m.species), m.steps, .!blocked)
    main = nfree == 0 ? 0 : comp[free[1]]
    U(j) = comp[m.steps[j].from]; B(j) = comp[m.steps[j].to]
    doors = [j for j in eachindex(m.steps) if blocked[j] && U(j) != B(j)]
    chords = [j for j in eachindex(m.steps) if blocked[j] && U(j) == B(j)]
    aff_chords = [j for j in chords if aff[j]]
    adoors = [j for j in doors if aff[j]]
    kdoors = [j for j in doors if !aff[j]]
    islands = [c for c in 1:nc if c != main]
    succ = [Int[] for _ in 1:nc]; pred = [Int[] for _ in 1:nc]
    for j in adoors; push!(succ[U(j)], B(j)); push!(pred[B(j)], U(j)); end
    function reach(src, nbrs)
        seen = falses(nc); seen[src] = true; q = [src]
        while !isempty(q)
            u = popfirst!(q)
            for v in nbrs[u]; seen[v] || (seen[v] = true; push!(q, v)); end
        end
        seen
    end
    clus, ncl = components(nc, [(from = U(j), to = B(j)) for j in adoors], trues(length(adoors)))
    # unit potentials inside each cluster; `unbalanced` when an affinity cycle needs unequal exponents
    φ = fill(typemin(Int), nc); unbalanced = false
    for r in 1:nc
        φ[r] == typemin(Int) || continue
        φ[r] = 0; q = [r]
        while !isempty(q)
            u = popfirst!(q)
            for v in succ[u]; φ[v] == typemin(Int) && (φ[v] = φ[u] + 1; push!(q, v)); end
            for v in pred[u]; φ[v] == typemin(Int) && (φ[v] = φ[u] - 1; push!(q, v)); end
        end
    end
    unbalanced = any(j -> φ[B(j)] - φ[U(j)] != 1, adoors)
    directed_cycle = any(j -> reach(B(j), succ)[U(j)], adoors)
    function rel_signs(ref)
        below, above = reach(ref, succ), reach(ref, pred)
        [clus[c] != clus[ref] ? 4 : c == ref ? 0 : below[c] && above[c] ? 2 :
         below[c] ? 1 : above[c] ? -1 : 2 for c in 1:nc]
    end
    kself = [j for j in kdoors if clus[U(j)] == clus[B(j)]]
    kedge = [j for j in kdoors if clus[U(j)] != clus[B(j)]]
    forced = [j for j in kself if φ[U(j)] != φ[B(j)]]
    sign = fill(3, nc)
    cyc_main = false; cyc_hang = false; bridges = Int[]; vanish_cycle = false; cyc_single_form = false
    if main != 0
        mainc = clus[main]
        rs = rel_signs(main)
        for c in 1:nc; clus[c] == mainc && (sign[c] = rs[c]); end
        cadj = [Tuple{Int, Int}[] for _ in 1:ncl]
        for j in kedge
            a, b = clus[U(j)], clus[B(j)]
            push!(cadj[a], (b, j)); push!(cadj[b], (a, j))
        end
        tin = zeros(Int, ncl); low = zeros(Int, ncl); timer = Ref(0); bridge = Set{Int}()
        function dfs(u, pe)
            timer[] += 1; tin[u] = low[u] = timer[]
            for (v, j) in cadj[u]
                j == pe && continue
                if tin[v] == 0
                    dfs(v, j); low[u] = min(low[u], low[v]); low[v] > tin[u] && push!(bridge, j)
                else
                    low[u] = min(low[u], tin[v])
                end
            end
        end
        dfs(mainc, 0)
        blk, nb = components(ncl, [(from = clus[U(j)], to = clus[B(j)]) for j in kedge],
                             [!(j in bridge) for j in kedge])
        blksize = [count(==(b), blk) for b in 1:nb]
        # the kinetic cycle block through the main cluster: its attachment points in the main cluster.
        # With kept ratios every island of the block has an order between the least and the largest
        # attachment order (steady state = weighted mean of its neighbours), so it vanishes when every
        # attachment vanishes ("vanishing cycle") and is rate-dependent otherwise.
        inblk(j) = !(j in bridge) && blk[clus[U(j)]] == blk[mainc]
        attach = [(clus[U(j)] == mainc ? U(j) : B(j), clus[U(j)] == mainc ? m.steps[j].from : m.steps[j].to)
                  for j in kedge if inblk(j) && (clus[U(j)] == mainc || clus[B(j)] == mainc)]
        vanishing = !isempty(attach) && all(a -> sign[a[1]] == 1, attach)
        single_form = length(unique(last.(attach))) == 1
        done = falses(ncl); done[mainc] = true; q = [mainc]
        while !isempty(q)
            u = popfirst!(q)
            for (v, j) in cadj[u]
                done[v] && continue
                done[v] = true; push!(q, v)
                cu, cv = clus[U(j)] == u ? (U(j), B(j)) : (B(j), U(j))
                if j in bridge && blksize[blk[v]] == 1
                    rsv = rel_signs(cv)
                    for c in 1:nc; clus[c] == v && (sign[c] = addsign(sign[cu], rsv[c])); end
                elseif blk[v] == blk[mainc] && vanishing
                    vanish_cycle = true
                    doorforms = unique([clus[U(k)] == v ? U(k) : B(k) for k in kedge
                                        if inblk(k) && (clus[U(k)] == v || clus[B(k)] == v)])
                    rels = [rel_signs(f) for f in doorforms]
                    for c in 1:nc
                        clus[c] == v || continue
                        sign[c] = all(r -> r[c] in (0, 1), rels) ? 1 : 2
                    end
                else
                    blk[v] == blk[mainc] ? (cyc_main = true) : (cyc_hang = true)
                    blk[v] == blk[mainc] && single_form && (cyc_single_form = true)
                    for c in 1:nc; clus[c] == v && (sign[c] = 3); end
                end
            end
        end
        bridges = collect(bridge)
    end
    isl_sign = [sign[c] for c in islands]
    cls = nfree != 1 ? "free_forms" :
          !isempty(aff_chords) || directed_cycle ? "infeasible" :
          isempty(islands) ? "no_island" :
          any(==(-1), isl_sign) ? "sequestering" :
          any(==(3), isl_sign) ? "undefined" :
          any(==(2), isl_sign) ? "ambiguous" :
          any(==(0), isl_sign) ? "phantom" : "emptied"
    # the task's literal tests, without order propagation: an island holds an unbound side of an
    # affinity door (sequestering) / every island holds a bound side and no unbound side (emptied)
    lit_seq = any(j -> U(j) != main, adoors)
    lit_emptied = all(c -> any(j -> B(j) == c, adoors) && !any(j -> U(j) == c, adoors), islands)
    (cls = cls, nfree = nfree, ncomp = nc, nisl = length(islands), nadoor = length(adoors),
     nkdoor = length(kdoors), nkself = length(kself), nkedge = length(kedge), forced = forced,
     aff_chords = aff_chords, kchords = [j for j in chords if !aff[j]], directed_cycle = directed_cycle,
     unbalanced = unbalanced, cyc_main = cyc_main, cyc_hang = cyc_hang, isl_sign = isl_sign,
     vanish_cycle = vanish_cycle, cyc_single_form = cyc_single_form,
     comp = comp, main = main, doors = doors, aff = aff, bridges = bridges, clus = clus,
     lit_seq = lit_seq, lit_emptied = lit_emptied)
end

# ─── Printing ──────────────────────────────────────────────────────────────

metname(x::Pop.PMet) = string(x.name, x.role === :Inh ? "::Inh" : "")
function form_text(sp::Pop.PSpecies)
    res = isempty(sp.res_added) && isempty(sp.res_subtracted) ? "" :
          "; residual = " * join([string.(sp.res_added); "-" .* string.(sp.res_subtracted)], " ")
    string(sp.conformation, "(", join(metname.(sp.bound), ", "), res, ")")
end
function step_dsl(m::Pop.PMech, s::Pop.PStep)
    side(f, ms) = join([form_text(m.species[f]); metname.(ms)], " + ")
    side(s.from, s.consumed) * (s.re ? " ⇌ " : " <--> ") * side(s.to, s.released)
end

"""Print every step with its group tag, its door type under `reading` and its component."""
function print_mech(io, m::Pop.PMech, reading::Symbol; a = analyze(m, reading))
    g = Pop.group_of_step(m)
    println(io, "  -- ", m.id, "  n=", m.cat_multiplicity, "  class[", reading, "]=", a.cls,
            "  comps=", a.ncomp, " main=", a.main, " island signs=", a.isl_sign)
    for (j, s) in enumerate(m.steps)
        door = j in a.doors ? (a.aff[j] ? "  AFFINITY door" : "  KINETIC door") :
               m.group_tags[g[j]] === :OnlyA ? "  (blocked chord)" : ""
        @printf(io, "     g%-2d %-10s %-62s c%d-c%d%s\n", g[j], "[" * string(m.group_tags[g[j]]) * "]",
                step_dsl(m, s), a.comp[s.from], a.comp[s.to], door)
    end
end
