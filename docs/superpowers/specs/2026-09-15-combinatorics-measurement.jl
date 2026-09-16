# ABOUTME: Measures per-move child counts, depth dependence, and flip/split closure sizes
# ABOUTME: on a bi-bi reaction with one allosteric regulator and one competitive inhibitor.

using EnzymeRates
using Statistics: median
const ER = EnzymeRates

say(x...) = (println(x...); flush(stdout))
const T0 = time()
elapsed() = string(round(time() - T0; digits = 1), " s elapsed\n")

# ─── mechanism printing ───────────────────────────────────────────────────────

function res_str(sp)
    ER.has_residual(sp) || return ""
    s = join([string(ER.name(a)) for a in ER.added(ER.residual(sp))], " + ")
    for r in ER.subtracted(ER.residual(sp))
        s *= " - " * string(ER.name(r))
    end
    "; residual = " * s
end

function sp_str(sp)
    b = String[]
    for x in ER.bound(sp)
        push!(b, x isa ER.CompetitiveInhibitor ? string(ER.name(x), "::Inh") :
                 x isa ER.AllostericRegulator ? string(ER.name(x), "::Allo") :
                 string(ER.name(x)))
    end
    inner = join(b, ", ") * res_str(sp)
    string(ER.conformation(sp), isempty(inner) ? "" : "(" * inner * ")")
end

function step_str(s)
    bm = ER.bound_metabolite(s)
    arrow = ER.is_equilibrium(s) ? "⇌" : "<-->"
    f, t = sp_str(ER.from_species(s)), sp_str(ER.to_species(s))
    bm === nothing && return "$f $arrow $t"
    if any(x -> x == bm, ER.bound(ER.to_species(s)))
        "$f + $(ER.name(bm)) $arrow $t"
    else
        "$f $arrow $t + $(ER.name(bm))"
    end
end

"Grouped step lists plus, for an allosteric mechanism, its site summary."
function mech_lines(m)
    tags = m isa ER.AllostericMechanism ? ER.cat_allo_states(m) : nothing
    out = String[]
    for (g, grp) in enumerate(ER.steps(m))
        body = length(grp) == 1 ? step_str(grp[1]) :
               "(" * join(step_str.(grp), ", ") * ")"
        tags === nothing || (body *= " :: " * string(tags[g]))
        push!(out, body)
    end
    if m isa ER.AllostericMechanism
        push!(out, "catalytic_multiplicity: $(ER.catalytic_multiplicity(m))")
        for site in ER.regulatory_sites(m)
            ls = join([string(ER.name(l), "::", st) for (l, st) in
                       zip(ER.ligands(site), ER.allo_states(site))], ", ")
            push!(out, "regulatory_site(multiplicity = " *
                       "$(ER.multiplicity(site))): $ls")
        end
    end
    out
end

sig(m) = join(mech_lines(m), " | ")

function describe(m)
    np = try
        string(ER._independent_param_count(m))
    catch e
        "err:" * first(sprint(showerror, e), 60)
    end
    say("    ", m isa ER.AllostericMechanism ? "AllostericMechanism" : "Mechanism",
        ", ", ER.n_steps(m), " steps, ", length(ER.steps(m)), " groups, ",
        ER._re_segment_count(m), " RE segments, ", np, " indep params")
    for l in mech_lines(m)
        say("      ", l)
    end
end

# ─── the reaction ─────────────────────────────────────────────────────────────

rxn = @enzyme_reaction begin
    substrates: A[C], B[N]
    products: P[C], Q[N]
    allosteric_regulators: R(1)
    competitive_inhibitors: I
end

# ─── 1. init_mechanisms ───────────────────────────────────────────────────────

say("## 1. `init_mechanisms`\n")

topos = ER._catalytic_topologies(rxn)
pats = ER._competition_patterns(Set([:A, :B]), Set([:P, :Q]))
inits = ER.init_mechanisms(rxn)

"`init_mechanisms` restricted to a single catalytic topology."
function init_for_topo(topo)
    ms = ER.Mechanism[]
    for (st, gr) in ER._expand_substrate_product_dead_ends([topo], rxn)
        mst, mgr = ER._apply_equivalence_grouping(st, gr)
        push!(ms, ER.Mechanism(rxn, ER._to_group_list(mst, mgr)))
    end
    ms
end

per_topo = [length(init_for_topo(t)) for t in topos]

say("| quantity | value |")
say("|---|---|")
say("| catalytic topologies | ", length(topos), " |")
say("| dead-end (competition) patterns | ", length(pats), " |")
say("| `init_mechanisms(rxn)` | ", length(inits), " |")
say("| seeds per topology: min | ", minimum(per_topo), " |")
say("| seeds per topology: median | ", median(per_topo), " |")
say("| seeds per topology: max | ", maximum(per_topo), " |")
say("| sum over topologies == `init_mechanisms` | ",
    sum(per_topo) == length(inits), " |")
say()
ord = sortperm(collect(zip(per_topo, [join(step_str.(t), " | ") for t in topos])))
say("Topology with the FEWEST init mechanisms (", per_topo[ord[1]], "):")
for s in topos[ord[1]]
    say("      ", step_str(s))
end
say("\nTopology with the MOST init mechanisms (", per_topo[ord[end]], "):")
for s in topos[ord[end]]
    say("      ", step_str(s))
end
say("\n", elapsed())

# ─── 2. seed_mechanisms ───────────────────────────────────────────────────────

say("## 2. `seed_mechanisms(rxn, Set([:R]), Set([:I]))`\n")
seeds = ER.seed_mechanisms(rxn, Set([:R]), Set([:I]))
sort!(seeds; by = sig)
n_allo = count(m -> m isa ER.AllostericMechanism, seeds)
say("| quantity | value |")
say("|---|---|")
say("| seeds | ", length(seeds), " |")
say("| allosteric seeds | ", n_allo, " |")
say("| plain `Mechanism` seeds | ", length(seeds) - n_allo, " |")
say("| steps per seed: min / median / max | ",
    minimum(ER.n_steps.(seeds)), " / ", median(ER.n_steps.(seeds)), " / ",
    maximum(ER.n_steps.(seeds)), " |")
say("\n", elapsed())

# ─── the seven moves ──────────────────────────────────────────────────────────

# `scope` is `:full` when the move is cheap enough to run over all 19_904 seeds
# and `:sample` when it is not (see the report's Method note).
moves = [
    ("_expand_re_to_ss", m -> ER._expand_re_to_ss(m), :full),
    ("_expand_split_kinetic_group",
     m -> ER._expand_split_kinetic_group(m), :sample),
    ("_expand_add_dead_end_regulator",
     m -> ER._expand_add_dead_end_regulator(m, rxn), :full),
    ("_expand_to_allosteric", m -> ER._expand_to_allosteric(m, rxn), :full),
    ("_expand_add_allosteric_regulator",
     m -> ER._expand_add_allosteric_regulator(m, rxn), :full),
    ("_expand_change_allo_state", m -> ER._expand_change_allo_state(m), :sample),
    ("_expand_merge_regulatory_sites",
     m -> ER._expand_merge_regulatory_sites(m), :full),
]

const SEED_SAMPLE_MAX = 500
const DEPTH1_SAMPLE_MAX = 200

"Every `k`-th element, at most `n`, from a vector already in canonical order."
function subsample(v, n)
    k = max(1, cld(length(v), n))
    s = v[1:k:end]
    (length(s) > n ? s[1:n] : s, k)
end

seed_sample, seed_k = subsample(seeds, SEED_SAMPLE_MAX)
sample_pos = Set(1:seed_k:length(seeds))

say("## 3. Children per parent, over the seed set of item 2\n")
say("Counts are raw per-move output; `expand_mechanisms` applies ",
    "`_filter_by_reg_type`\nto the union, so its child count is ≤ the row sums ",
    "below. Scope `all` = all ", length(seeds), "\nseeds; scope `sample` = every ",
    seed_k, "-th seed in canonical order (", length(seed_sample), " parents).\n")
say("| move | scope | parents | parents with ≥1 child | children | min | ",
    "median | max |")
say("|---|---|---|---|---|---|---|---|")

# Children of the sampled seeds under every move — the depth-1 parent pool.
kids = Union{ER.Mechanism, ER.AllostericMechanism}[]
extremes = Tuple{String, Any, Any, Int, Int}[]
for (nm, f, scope) in moves
    parents = scope === :full ? seeds : seed_sample
    counts = Int[]
    for (i, m) in enumerate(parents)
        cs = f(m)
        push!(counts, length(cs))
        collect_it = scope === :sample || i in sample_pos
        collect_it && append!(kids, cs)
    end
    i_lo, i_hi = argmin(counts), argmax(counts)
    say("| `", nm, "` | ", scope === :full ? "all" : "sample", " | ",
        length(parents), " | ", count(>(0), counts), " | ", sum(counts), " | ",
        minimum(counts), " | ", median(counts), " | ", maximum(counts), " |")
    push!(extremes, (nm, parents[i_lo], parents[i_hi], counts[i_lo], counts[i_hi]))
end
say()
for (nm, lo, hi, nlo, nhi) in extremes
    if nhi == 0
        say("**", nm, "** emits no children on any parent in scope.\n")
        continue
    end
    say("**", nm, "** — parent at the minimum (", nlo, " children):")
    describe(lo)
    say("**", nm, "** — parent at the maximum (", nhi, " children):")
    describe(hi)
    say()
end

# ─── 4. one level down ────────────────────────────────────────────────────────

say(elapsed())
say("## 4. Children per parent, one level down\n")
unique!(kids)
sort!(kids; by = sig)
say("depth-1 children (all moves over the ", length(seed_sample),
    " sampled seeds, deduped): ", length(kids))
depth1_sample, d1_k = subsample(kids, DEPTH1_SAMPLE_MAX)
say("parents: every ", d1_k, "-th in canonical order — ", length(depth1_sample),
    " mechanisms\n")
say("| move | parents with ≥1 child | children | min | median | max |")
say("|---|---|---|---|---|---|")
extremes2 = Tuple{String, Any, Any, Int, Int}[]
for (nm, f, _) in moves
    counts = [length(f(m)) for m in depth1_sample]
    i_lo, i_hi = argmin(counts), argmax(counts)
    say("| `", nm, "` | ", count(>(0), counts), " | ", sum(counts), " | ",
        minimum(counts), " | ", median(counts), " | ", maximum(counts), " |")
    push!(extremes2, (nm, depth1_sample[i_lo], depth1_sample[i_hi],
                      counts[i_lo], counts[i_hi]))
end
say()
for (nm, lo, hi, nlo, nhi) in extremes2
    if nhi == 0
        say("**", nm, "** emits no children on any parent in scope.\n")
        continue
    end
    say("**", nm, "** — parent at the minimum (", nlo, " children):")
    describe(lo)
    say("**", nm, "** — parent at the maximum (", nhi, " children):")
    describe(hi)
    say()
end

# ─── 5. closures ──────────────────────────────────────────────────────────────

say(elapsed())
say("## 5. Flip and split closures\n")

const CLOSURE_CAP = 20_000
const CLOSURE_SECONDS = 120.0

"""
Breadth-first closure of `seed` under repeated application of `gen` alone.
Returns `(nodes, stop)` where `stop` is `""`, `"cap"` or `"time"`.
"""
function closure(seed, gen)
    seen = Set{UInt64}([hash(seed)])
    nodes = Any[seed]
    queue = Any[seed]
    t0 = time()
    stop = ""
    while !isempty(queue) && isempty(stop)
        m = popfirst!(queue)
        for c in gen(m)
            h = hash(c)
            h in seen && continue
            push!(seen, h)
            push!(nodes, c)
            push!(queue, c)
        end
        length(nodes) >= CLOSURE_CAP && (stop = "cap")
        time() - t0 > CLOSURE_SECONDS && isempty(stop) && (stop = "time")
    end
    (nodes, stop)
end

"The RE segmentation of `m`: the partition of its enzyme forms into RE segments."
function seg_key(m)
    cm = m isa ER.AllostericMechanism ? ER._state_mechanism(m, :A) : m
    species, groups, _ = ER._compute_re_groups(cm)
    Tuple(sort([Tuple(sort([string(ER.name(species[i])) for i in g]))
                for g in groups]))
end

byorder = sortperm(collect(zip([-ER.n_steps(m) for m in seeds], sig.(seeds))))
# The three seeds with the most steps, and — because the split closure of a
# largest seed does not finish inside the budget — the three with the fewest.
probes = vcat([("most-$i", seeds[byorder[i]]) for i in 1:3],
              [("fewest-$i", seeds[byorder[end + 1 - i]]) for i in 1:3])

say("Closure caps: ", CLOSURE_CAP, " visited nodes or ", CLOSURE_SECONDS,
    " s per closure.\n")
say("| seed | steps | groups | RE segments | indep params | flip closure | ",
    "distinct RE segmentations | split closure |")
say("|---|---|---|---|---|---|---|---|")
for (lbl, m) in probes
    fl, fstop = closure(m, ER._expand_re_to_ss)
    sp, sstop = closure(m, ER._expand_split_kinetic_group)
    nsegm = length(unique(seg_key.(fl)))
    say("| ", lbl, " | ", ER.n_steps(m), " | ", length(ER.steps(m)), " | ",
        ER._re_segment_count(m), " | ", ER._independent_param_count(m), " | ",
        length(fl), isempty(fstop) ? "" : " ($fstop hit)", " | ", nsegm, " | ",
        length(sp), isempty(sstop) ? "" : " ($sstop hit)", " |")
end
say()
for (lbl, m) in probes
    say("seed ", lbl, ":")
    describe(m)
    say()
end
say(elapsed())
say("DONE")
