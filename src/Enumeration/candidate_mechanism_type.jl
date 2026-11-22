# ABOUTME: CandidateMechanism, a mechanism with the enumerator's record of each form, the
# ABOUTME: structural tests a candidate must pass, and the helpers moves and seeds share.

# ─── The candidate and its form records ─────────────────────────────────

"""
    _FormRecord

The enumerator's record of one enzyme form: its conformation label (`label`), the names it
binds productively (`productive`) and as modifier copies (`copies`, a pure inhibitor among
them), and its covalent residual, the substrates it adds (`adds`) and the products it
subtracts (`subtracts`). Each list is sorted and holds a name as often as the form does.
"""
struct _FormRecord
    label::Symbol
    productive::Vector{Symbol}
    copies::Vector{Symbol}
    adds::Vector{Symbol}
    subtracts::Vector{Symbol}
end

Base.:(==)(a::_FormRecord, b::_FormRecord) =
    a.label === b.label && a.productive == b.productive && a.copies == b.copies &&
    a.adds == b.adds && a.subtracts == b.subtracts
Base.hash(f::_FormRecord, h::UInt) =
    hash(f.subtracts, hash(f.adds, hash(f.copies, hash(f.productive,
         hash(f.label, hash(:_FormRecord, h))))))

"""
The name of the form whose record is `f`, in a mechanism whose reactants are `reactants`:
its label, then the names it binds in name order, a reactant's modifier copy written
`<name>mod` after its reactant, then any residual as `_res`, `_plus<name>` per substrate it
adds and `_minus<name>` per product it subtracts (`E`, `EAB`, `EPmodQ`,
`EB_res_plusA_minusP`). The reaction's names spell every name one way
(`_check_decodable_names`), so two forms take one name only when their records are equal.
"""
function _form_name(f::_FormRecord, reactants::Vector{Symbol})
    io = IOBuffer()
    print(io, f.label)
    for (x, copy) in sort!([[(x, false) for x in f.productive];
                            [(x, true) for x in f.copies]])
        print(io, x, copy && x in reactants ? "mod" : "")
    end
    if !isempty(f.adds) || !isempty(f.subtracts)
        print(io, "_res")
        foreach(x -> print(io, "_plus", x), f.adds)
        foreach(x -> print(io, "_minus", x), f.subtracts)
    end
    Symbol(String(take!(io)))
end

"""
    CandidateMechanism

A mechanism the enumerator proposes for a reaction (`candidate_mechanisms`). Its field
`mechanism`, a `Mechanism`, is the mechanism to fit and report. Two candidates are `==`
and hash alike exactly when their mechanisms are.
"""
struct CandidateMechanism
    mechanism::Mechanisms.Mechanism
    # The record of each form, at the form's index in `mechanism`'s catalytic scheme. The
    # moves read a form's contents here, never from the mechanism.
    forms::Vector{_FormRecord}
end

# Equal mechanisms hold equal records: the enumerator names each form after its record.
Base.:(==)(a::CandidateMechanism, b::CandidateMechanism) = a.mechanism == b.mechanism
Base.hash(c::CandidateMechanism, h::UInt) = hash(c.mechanism, h)

"""
    _candidate(metabolites, forms, substrates, products, steps, allostery = <plain>)
        -> CandidateMechanism
    _candidate(parent::CandidateMechanism, allosteric_states, oligomeric_state,
               regulatory_sites) -> CandidateMechanism

The candidate whose catalytic scheme has the metabolites `metabolites`, of which
`substrates` and `products` index the reactants, the forms whose records are `forms`, each
named after its record (`_form_name`), and the steps `steps` (`from` and `to` index
`forms`, the lists `metabolites`), with the allostery `allostery`, whose
`allosteric_states[l]` is the state of written kinetic group `l`, built and checked by the
`Mechanism` constructor; a plain mechanism by default. A record given twice names one
form; two different records that take one name are an internal error. The second gives
`parent`'s catalytic scheme, with its records, the allosteric states, oligomeric state and
regulatory sites given. The only place in `Enumeration` that builds a mechanism; the
candidate holds each form's record in the order the constructor
stores the forms.
"""
function _candidate(metabolites::Vector{Symbol}, forms::Vector{_FormRecord},
                    substrates::Vector{Int}, products::Vector{Int},
                    steps::Vector{Mechanisms.Step},
                    allostery::Mechanisms.Allostery =
                        Mechanisms.Allostery(Symbol[], 1, Mechanisms.RegulatorySite[]))
    reactants = metabolites[[substrates; products]]
    names = Symbol[_form_name(f, reactants) for f in forms]
    record = Dict{Symbol, _FormRecord}()
    for (n, f) in zip(names, forms)
        get!(record, n, f) == f ||
            error("_candidate: two different form records are both named $n")
    end
    distinct = unique(names)
    if length(distinct) < length(names)
        at = Dict{Symbol, Int}(n => i for (i, n) in enumerate(distinct))
        steps = Mechanisms.Step[Mechanisms.Step(at[names[e.from]], at[names[e.to]],
                                                e.taken_up, e.given_off,
                                                e.rapid_equilibrium, e.kinetic_group)
                                for e in steps]
    end
    m = Mechanisms.Mechanism(
        Mechanisms.CatalyticScheme(metabolites, distinct, substrates, products, steps),
        allostery)
    CandidateMechanism(m, _FormRecord[record[n] for n in m.catalytic.enzyme_forms])
end
_candidate(parent::CandidateMechanism, allosteric_states::Vector{Symbol},
           oligomeric_state::Int, regulatory_sites::Vector{Mechanisms.RegulatorySite}) =
    CandidateMechanism(
        Mechanisms.Mechanism(parent.mechanism.catalytic,
                             Mechanisms.Allostery(allosteric_states, oligomeric_state,
                                                  regulatory_sites)),
        parent.forms)

# ─── The structural tests a candidate must pass ─────────────────────────

"""Whether the rate of graph `m`, edge `k` at rapid equilibrium where `re[k]` holds, has a
maximum as the metabolites `names` (its substrates or its products, by index) grow: turning
rapid equilibrium every steady-state edge that takes up or gives off one of them leaves no
rapid-equilibrium turnover cycle (condition V). Only the edges `on` flags count, so an
allosteric mechanism's test can read its active law (`_degenerate`)."""
function _has_vmax(m::Mechanisms.CatalyticScheme, re::BitVector, names::Vector{Int},
                   on::BitVector = trues(length(m.steps)))
    touches = BitVector([_exchanges(m, k, names) for k in eachindex(m.steps)])
    !Mechanisms.has_turnover_cycle(m, (re .| touches) .& on)
end

"""Whether the metabolite list `s` holds a metabolite of `names`, both indices into a
scheme's metabolites."""
_holds(s::Vector{Int}, names::Vector{Int}) = any(in(names), s)

"""Whether step `k` of `m` takes up or gives off a metabolite of `names` (`_holds`)."""
_exchanges(m::Mechanisms.CatalyticScheme, k::Int, names::Vector{Int}) =
    _holds(m.steps[k].taken_up, names) || _holds(m.steps[k].given_off, names)

"""Whether step `k` of `m` takes up a metabolite of one of `subs` and `prods` (a reactant
of one side of the reaction) and gives off one of the other, as a Theorell–Chance step
does."""
function _crosses_sides(m::Mechanisms.CatalyticScheme, k::Int, subs::Vector{Int},
                        prods::Vector{Int})
    t, g = m.steps[k].taken_up, m.steps[k].given_off
    _holds(t, subs) && _holds(g, prods) || _holds(t, prods) && _holds(g, subs)
end

"""
Whether some chemistry node of the catalytic graph `m` of candidate `c`, edge `k` at rapid
equilibrium where `re[k]` holds, is left both by a rapid-equilibrium edge that releases a
substrate and by one that releases a product (condition C). A chemistry node is a merged
complex (a form entered by a fused binding of a substrate) with every form joined to it by
RE isomerizations, a set of forms joined to each other by RE isomerizations, or an RE
Theorell–Chance step. An RE edge leaves a node through a node form F and releases M when
traversing it away from F gives M off: F is its `to` form and M is consumed, or F is its
`from` form and M is released. An RE Theorell–Chance node leaves through its own edge in
both directions. Then the chemistry sits in rapid equilibrium with both sides and some
reactant is not needed at zero products or substrates. Only the edges `on` flags count
(`_degenerate`).
"""
function _chemistry_equilibrates_both_sides(
        c::CandidateMechanism, re::BitVector,
        on::BitVector = trues(length(c.mechanism.catalytic.steps)))
    m = c.mechanism.catalytic
    re = re .& on
    subs, prods = m.substrates, m.products
    releases(node, side) = any(eachindex(m.steps)) do k
        (; from, to) = m.steps[k]
        re[k] &&
            (to in node && !(from in node) && _holds(m.steps[k].taken_up, side) ||
             from in node && !(to in node) && _holds(m.steps[k].given_off, side))
    end
    # The forms joined by RE isomerizations, one node per segment; a merged complex
    # outside every RE isomerization is a node of its own.
    isos = BitVector([re[k] && _is_iso(m.steps[k]) for k in eachindex(m.steps)])
    s = Mechanisms.rapid_equilibrium_segments(m, isos, isos)
    nodes = [s.segments; [[x] for x in _merged_complexes(c, subs, on)
                          if s.segment_of[x] == 0]]
    any(node -> releases(node, subs) && releases(node, prods), nodes) ||
        any(k -> re[k] && _crosses_sides(m, k, subs, prods), eachindex(m.steps))
end

"""Whether the mechanism of candidate `c` is degenerate by the structural conditions V and
C: no maximal rate in one direction (`_has_vmax`), or chemistry in rapid equilibrium with
both sides (`_chemistry_equilibrates_both_sides`). Read on the catalytic graph with its own
flags, over the steps the active conformation runs
(`steps_running_in_allosteric_conformations`): every edge of a plain mechanism, and the
edges its tags keep of an allosteric one. The beam expands a degenerate seed without
fitting it (`_base_tier`)."""
function _degenerate(c::CandidateMechanism)
    cm = c.mechanism.catalytic
    re = _rapid_equilibrium_flags(cm)
    on = Mechanisms.steps_running_in_allosteric_conformations(c.mechanism).active
    !(_has_vmax(cm, re, cm.substrates, on) &&
      _has_vmax(cm, re, cm.products, on)) ||
        _chemistry_equilibrates_both_sides(c, re, on)
end

"""
    _kept_forms_keep_a_binding(m) -> Bool

Whether, in each conformation of `m`, every form but free E that the conformation's law
keeps (a form of a step it runs, `steps_running_in_allosteric_conformations`) and that
some binding enters (a step that takes up one metabolite and gives off none) is entered by a
binding the conformation does not exclude: one outside the `:OnlyI` groups in the active
conformation, outside the `:OnlyA` groups in the inactive one. True for a plain mechanism.
The guard admits a mechanism that fails this when its law is the limit, so a user may write
one; the enumerator does not build it.
"""
function _kept_forms_keep_a_binding(am::Mechanisms.Mechanism)
    # A plain mechanism excludes no binding.
    _is_allosteric(am) || return true
    steps, n = am.catalytic.steps, length(am.catalytic.enzyme_forms)
    runs = Mechanisms.steps_running_in_allosteric_conformations(am)
    for (law, excluded_state) in ((runs.active, :OnlyI), (runs.inactive, :OnlyA))
        # Per form: some binding enters it; some binding the conformation runs does. A
        # step with one metabolite is stored as its uptake, so it enters its `to` form.
        entered, made = falses(n), falses(n)
        for (k, e) in enumerate(steps)
            _is_binding(e) || continue
            entered[e.to] = true
            am.allosteric.allosteric_states[e.kinetic_group] === excluded_state ||
                (made[e.to] = true)
        end
        for (k, e) in enumerate(steps), f in (e.from, e.to)
            law[k] && f != 1 && entered[f] && !made[f] && return false
        end
    end
    true
end

"""
Whether the enumerator may only give this mechanism a catalytic scheme that
passes `_hyperbolic_catalysis`. True when a conformational equilibrium sits over
a catalytic site with more than one subunit: the equilibrium then raises every
catalytic-site binding to the power of the oligomeric state. Over one subunit, as for every
plain mechanism, it only reweights each enzyme form and adds no power of its own.
"""
_requires_hyperbolic_catalysis(m::Mechanisms.Mechanism) = m.allosteric.oligomeric_state > 1

"""Whether candidate `c` passes `_hyperbolic_catalysis` where its mechanism
`_requires_hyperbolic_catalysis`; true where it does not."""
_hyperbolic_where_required(c::CandidateMechanism) =
    !_requires_hyperbolic_catalysis(c.mechanism) || _hyperbolic_catalysis(c)

"""
    _hyperbolic_catalysis(c::CandidateMechanism, on::BitVector) -> Bool
    _hyperbolic_catalysis(c::CandidateMechanism) -> Bool

Whether the King–Altman denominator of the catalytic scheme `m` of candidate `c` has degree
at most 1 in every substrate and product concentration. Every binding of a substrate or
product at its catalytic site counts, abortive complexes included. Edges touching a form
that holds a modifier are left out (a mask of the catalytic-site edges): a modifier binds a
site of its own by definition, so the powers its binding adds, including those of a
reactant's modifier copy, are a separate source that conformational mechanisms keep. Only
the edges `on` flags count, and `rapid_equilibrium_segments` leaves out the forms they do
not keep, so an allosteric mechanism's test can read each conformation's law. Without `on`,
the test reads each conformation's law (`steps_running_in_allosteric_conformations`), once
when the two run the same steps, as both conformations of a plain mechanism run every step:
a conformation that runs part of the catalytic graph can carry a power the whole graph
lacks, as a ping-pong's abortive complexes put A² into the inactive law once that
conformation stops one half-reaction.

The equation is a sum over rapid-equilibrium (RE) segments and spanning
arborescences of the segment graph toward each segment. A denominator term is
the root segment's weight times the weight of every tree edge, so the exponent
of `X` in a term is the most `X` any form of the root segment carries beyond
the segment's bottom form, plus, per tree edge, one if the step binds `X` in
the tree direction and the count of `X` the edge's source form carries beyond
its bottom (`rapid_equilibrium_segments`). The degree exceeds 1 exactly when one
segment or one edge scores 2 or more, or the root scores 1 and some
arborescence toward it holds a scoring edge, or some arborescence holds two
scoring edges. An arborescence toward `S` containing given edges exists iff
every segment still reaches `S` once each given edge's source keeps that edge
as its only way out (`_all_reach`).
"""
function _hyperbolic_catalysis(c::CandidateMechanism, on::BitVector)
    m = c.mechanism.catalytic
    on_site = on .& BitVector([isempty(c.forms[e.from].copies) &&
                               isempty(c.forms[e.to].copies) for e in m.steps])
    s = Mechanisms.rapid_equilibrium_segments(m, on_site)
    # Directed segment-graph edges: source segment, target segment, source form, the
    # metabolites (their rows of the offsets, their indices) bound in that direction.
    edges = Tuple{Int, Int, Int, Vector{Int}}[]
    for k in s.steps
        (; from, to) = m.steps[k]
        a, b = s.segment_of[from], s.segment_of[to]
        push!(edges, (a, b, from, m.steps[k].taken_up))
        push!(edges, (b, a, to, m.steps[k].given_off))
    end
    n = length(s.segments)
    for x in [m.substrates; m.products]
        score(e) = count(==(x), e[4]) + s.offsets[x, e[3]]
        carrying = [e for e in edges if score(e) > 0]
        any(e -> score(e) > 1, carrying) && return false
        root_score(k) = maximum((s.offsets[x, f] for f in s.segments[k]); init=0)
        any(k -> root_score(k) > 1, 1:n) && return false
        roots = [k for k in 1:n if root_score(k) == 1]
        for e in carrying, k in roots
            k != e[1] && _all_reach(n, edges, k, (e,)) && return false
        end
        for (p, e1) in enumerate(carrying), e2 in carrying[p + 1:end]
            e1[1] == e2[1] && continue
            any(k -> k != e1[1] && k != e2[1] && _all_reach(n, edges, k, (e1, e2)),
                1:n) && return false
        end
    end
    true
end

function _hyperbolic_catalysis(c::CandidateMechanism)
    runs = Mechanisms.steps_running_in_allosteric_conformations(c.mechanism)
    all(mask -> _hyperbolic_catalysis(c, mask), unique!([runs.active, runs.inactive]))
end

"""
Whether every segment of the segment graph reaches `root` when each edge in
`fixed` is its source segment's only way out. A digraph has a spanning
arborescence toward `root` iff every vertex reaches `root`, and with the fixed
edges as their sources' only exits every such arborescence contains them.
"""
function _all_reach(n::Int, edges, root::Int, fixed)
    pinned = Dict(e[1] => e[2] for e in fixed)
    into = [Int[] for _ in 1:n]
    for (u, v, _, _) in edges
        get(pinned, u, v) == v && push!(into[v], u)
    end
    seen = falses(n)
    seen[root] = true
    queue = [root]
    while !isempty(queue)
        v = popfirst!(queue)
        for u in into[v]
            seen[u] || (seen[u] = true; push!(queue, u))
        end
    end
    all(seen)
end

"""
A class of the dwell gauge (`_gauge_rescaling`): a kinetic group's class (`group > 0`)
or a twin form's class of its own (`form > 0`), in one conformational state (`state` 1,
or 2 for the inactive one) or in both (0). Group and form indices count apart, so group
2 and form 2 are two classes.
"""
struct _GaugeClass
    group::Int
    form::Int
    state::Int
end

"""The class of the factor 1, which every form outside a copy's complexes and twins
takes."""
const _NO_RESCALING = _GaugeClass(0, 0, 0)

"""The gauge class of group `h` or form `t` (`_GaugeClass`) in a mechanism with one
conformational state."""
_one_state_class(h::Int, t::Int) = _GaugeClass(h, t, 1)

"""
    _productive_twin(c, mask) -> (site, ligand) -> Int

The productive form of the sub-graph of the catalytic graph `m` of candidate `c` whose edges
`mask` keeps, a form that holds no modifier, whose weight is proportional to that of the
complex a modifier named `ligand` forms at form `site` by a rapid-equilibrium step, or 0
when none is: a form of the site's rapid-equilibrium segment with the complex's offsets
(`rapid_equilibrium_segments`; the modifier counts as its name, so a reactant's copy as the
reactant), the one of the complex's composition when several are, else the first in form
order. Within a segment rapid equilibrium fixes each weight relative to another as a
constant times the metabolites bound along the route between them, so equal offsets mean
proportional weights, and steady-state data see only their sum. The offsets reach a copy
across a rapid-equilibrium isomerization (a ping-pong second chemistry step, a
conformational isomer), where no form shares the complex's composition. A form of the same
composition in another segment, formed by a steady-state binding, is no twin: its weight is
not a constant multiple of the complex's, merging the two is no reparameterization, and the
copy's constant may stay visible.

Only productive forms are twin sources, so the complex itself never is one. A complex
that duplicates only a copy's complex shares its weight with that copy's constant, and
the sites that pin either copy may keep both constants separable, so such a match never
makes a twin. The forms are `m`'s own, so the twins of a mechanism the dead-end move
builds are found on that mechanism (`_dead_end_child`).
"""
function _productive_twin(c::CandidateMechanism, mask::BitVector)
    m = c.mechanism.catalytic
    s = Mechanisms.rapid_equilibrium_segments(m, mask)
    # A form's bound names, sorted, its conformation and its residual.
    composition(f) = (sort!([c.forms[f].productive; c.forms[f].copies]), c.forms[f].label,
                      (c.forms[f].adds, c.forms[f].subtracts))
    function (site::Int, ligand::Symbol)
        # The ligand's row of the offsets, its index among the metabolites.
        n = findfirst(==(ligand), m.metabolites)
        n === nothing && return 0
        twins = [f for f in s.segments[s.segment_of[site]]
                 if isempty(c.forms[f].copies) &&
                    all(x -> s.offsets[x, f] == s.offsets[x, site] + (x == n),
                        axes(s.offsets, 1))]
        isempty(twins) && return 0
        bound_names, label, residual = composition(site)
        wanted = (sort!(push!(bound_names, ligand)), label, residual)
        k = findfirst(t -> composition(t) == wanted, twins)
        k === nothing ? minimum(twins) : twins[k]
    end
end

"""The name of the modifier copy that step `k` of candidate `c` binds (a binding,
`_is_binding`, that takes up a copy, `_copy_change`), or `nothing` when the step binds no
modifier copy."""
function _bound_copy(c::CandidateMechanism, k::Int)
    cm = c.mechanism.catalytic
    e = cm.steps[k]
    _is_binding(e) || return nothing
    x = cm.metabolites[only(e.taken_up)]
    _copy_change(c, k, x) > 0 ? x : nothing
end

"""
Whether the edges `ks` of one kinetic group of candidate `c` bind a modifier, and only at
sites whose complex has a productive twin (`twin`, a `_productive_twin` of one state's
graph).
"""
function _all_twin(c::CandidateMechanism, ks::Vector{Int}, twin)
    ligand = _bound_copy(c, first(ks))
    steps = c.mechanism.catalytic.steps
    ligand !== nothing && all(k -> twin(steps[k].from, ligand) != 0, ks)
end

"""
    _gauge_rescaling(c, members, g, twin, label)
        -> Union{Dict{Int, Tuple{_GaugeClass, _GaugeClass}}, Nothing}

The dwell gauge of the modifier copy bound by kinetic group `g` of the catalytic graph `m`
of candidate `c`, in one conformational state's graph, whose groups hold the edges `members`
lists: the rescaling each nonempty kinetic group needs, or `nothing` when the gauge does not
exist. It does not exist when group `g` binds no modifier or a complex it forms has no
productive twin (`_all_twin`, with `twin` a `_productive_twin` of that graph), or when some
group would need two rescalings. The gauge gives every form a factor σ: ρ_T to each twin T,
one factor s to every complex of the copy and 1 to every other form. Multiplying a
rapid-equilibrium group's K by σ(from)/σ(to), and a steady-state group's forward constant by
σ(from) and its reverse by σ(to), divides each form's weight by its factor and leaves every
flux as it was. A twin's weight is a constant multiple of its complexes', so the factors can
keep each twin's total with its complexes, and the rate law, unchanged while K* moves: the
copy's constant is a phantom. A kinetic group shares its constants, so every group, the
copy's own included, must have one ratio σ(from)/σ(to) over its edges (rapid equilibrium),
or one σ(from) and one σ(to) (steady state). The copy's group, at 1 on each site, holds all
its complexes at the one factor s, and a complex never takes its twin's factor: its weight
moves opposite to the twin's. A rapid-equilibrium mirror edge, between two complexes,
therefore has ratio 1, as its parent edge between two sites has, whatever the twins; a
steady-state mirror needs s beside its parent edge's 1, and the gauge fails. When each site
binds the copied ligand at rapid equilibrium, the gauge is the merge of each complex into
its twin, the twin's binding group taking K_h + K* in place of K_h, and the family is that
of the mechanism without the copy.

The ρ are generic numbers, equal only where the structure makes them equal: twins
formed from their sites by RE bindings of the copied ligand in one kinetic group h, each
the twin of one complex, share one factor (the merge gives it as ρ = K_h/(K_h + K*));
every other twin has its own. `label(h, t)` names the classes (`_GaugeClass`): it
receives the group of a shared class (h) or of the copy's complexes (g) with `t = 0`, or
`h = 0` and the twin `t` of a class of its own, so the caller decides which classes are
one number. A group's rescaling is the pair of its ends' classes, `_NO_RESCALING` standing
for the factor 1, and two of them for a rapid-equilibrium group whose two ends scale
alike. The test is bookkeeping over the edges, with no parameters and no numerics.
"""
function _gauge_rescaling(c::CandidateMechanism, members::Vector{Vector{Int}}, g::Int,
                          twin, label)
    m = c.mechanism.catalytic
    _all_twin(c, members[g], twin) || return nothing
    ligand = _bound_copy(c, first(members[g]))::Symbol
    twin_of = Dict(m.steps[k].to => twin(m.steps[k].from, ligand)
                   for k in members[g])
    # The rapid-equilibrium bindings of the ligand itself, not of a modifier copy.
    binding_group = Dict{Tuple{Int, Int}, Int}()
    for (h, ks) in enumerate(members), k in ks
        e = m.steps[k]
        e.rapid_equilibrium && _is_binding(e) &&
            m.metabolites[only(e.taken_up)] === ligand && _copy_change(c, k, ligand) == 0 &&
            (binding_group[(e.from, e.to)] = h)
    end
    class = Dict{Int, _GaugeClass}(complex => label(g, 0) for complex in keys(twin_of))
    for k in members[g]
        (; from, to) = m.steps[k]
        t = twin_of[to]
        haskey(class, t) && continue
        h = count(==(t), values(twin_of)) == 1 ? get(binding_group, (from, t), 0) : 0
        class[t] = h > 0 ? label(h, 0) : label(0, t)
    end
    σ(f) = get(class, f, _NO_RESCALING)
    rescale(k) = (a = σ(m.steps[k].from); b = σ(m.steps[k].to);
                  m.steps[k].rapid_equilibrium && a == b ?
                      (_NO_RESCALING, _NO_RESCALING) : (a, b))
    rescaling = Dict{Int, Tuple{_GaugeClass, _GaugeClass}}()
    for (h, ks) in enumerate(members)
        isempty(ks) && continue
        r = rescale(first(ks))
        all(k -> rescale(k) == r, ks) || return nothing
        rescaling[h] = r
    end
    rescaling
end

"""
    _redundant_copy_groups(c) -> Vector{Int}

Kinetic groups of the mechanism `m` of candidate `c` that bind a modifier redundantly: in
every conformational state where the copy binds, every complex has a productive twin, a form
whose weight is proportional to the complex's (`_productive_twin`, `_all_twin`), and the
dwell gauge exists over the states together (`_gauge_rescaling`). The gauge then moves the
copy's constant along a direction the rate law cannot see, so the copy adds a phantom, and
the mechanism's family is its family without the copy (`_gauge_rescaling`). A copy group
whose complexes all have twins but whose gauge fails is not redundant: a group that forms or
leaves a twin also binds where no copy does (a shared group that pins the existing binding),
a group holds a step that leaves a complex beside one that leaves another form (a second
copy bound at the complex and at a twin), or a mirror is steady-state. The gauge is the
proof that the copy's constant is invisible, and without it the constant may be
identifiable. Two copies that are twins only of each other, neither pinned, pass the test
and carry a phantom: the rule admits them.

An allosteric copy binds the active state unless its tag is `:OnlyI` and the inactive
state unless it is `:OnlyA`. Each state's graph is its law's
(`steps_running_in_allosteric_conformations`): the edges its tags keep between forms they
join to free E. The copy binds in a state only at the sites that
graph keeps; in a state where it binds nothing, every factor of that state is 1. A copy
that binds in neither state is not redundant: its constant enters neither law and is not
fitted. The two states form one system.
A group tagged `:EqualAI` has one set of constants in both, so it must take the same
rescaling in each. A shared class is one number in both states only when its binding
group and the copy are both `:EqualAI` (one K_h, one K*). The copy's complexes share one
factor in both states when the copy is `:EqualAI` (one K*). Every other class is a number
per state.
"""
function _redundant_copy_groups(c::CandidateMechanism)
    am = c.mechanism
    cm = am.catalytic
    members = _kinetic_groups(cm)
    any(ks -> _bound_copy(c, first(ks)) !== nothing, members) || return Int[]
    # A plain mechanism has one state; an allosteric one has two, which form one system.
    if !_is_allosteric(am)
        twin = _productive_twin(c, trues(length(cm.steps)))
        return [g for g in eachindex(members)
                if _gauge_rescaling(c, members, g, twin, _one_state_class) !== nothing]
    end
    equal(h) = am.allosteric.allosteric_states[h] === :EqualAI
    runs = Mechanisms.steps_running_in_allosteric_conformations(am)
    masks = (runs.active, runs.inactive)
    kept = [[filter(k -> mask[k], ks) for ks in members] for mask in masks]
    twins = [_productive_twin(c, mask) for mask in masks]
    # `s` indexes the conformational state: 1 active, 2 inactive.
    unscaled(s) = Dict(h => (_NO_RESCALING, _NO_RESCALING)
                       for h in eachindex(kept[s]) if !isempty(kept[s][h]))
    filter(collect(eachindex(members))) do g
        all(s -> isempty(kept[s][g]), 1:2) && return false
        one_number(h) = h > 0 && equal(h) && equal(g)
        label(state) = (h, t) -> _GaugeClass(h, t, one_number(h) ? 0 : state)
        rescaling(s) = isempty(kept[s][g]) ? unscaled(s) :
                       _gauge_rescaling(c, kept[s], g, twins[s], label(s))
        active = rescaling(1)
        active === nothing && return false
        inactive = rescaling(2)
        inactive === nothing && return false
        all(h -> !equal(h) || !haskey(inactive, h) || active[h] == inactive[h],
            keys(active))
    end
end

# ─── Atom-conservation validation ────────────────────────────

"""
Net atom multiset of the reactants named in `plus` minus those named in `minus`, with
zero counts dropped. A name may repeat; a name that names no reactant carries no atoms.
"""
function _net_atoms(reaction::Mechanisms.EnzymeReaction, plus, minus)
    acc = Dict{Symbol,Int}()
    for ra in [reaction.substrates; reaction.products]
        k = count(==(ra.name), plus) - count(==(ra.name), minus)
        for (a, c) in ra.atoms
            acc[a] = get(acc, a, 0) + k * c
        end
    end
    filter!(kv -> kv.second != 0, acc)
end

"""
    _assert_emission_rules(c::CandidateMechanism)

The two rules every mechanism the moves emit satisfies, checked on a parent candidate before
it is expanded. Every steady-state kinetic group holds an edge that carries net flux
(`_flux_carrying_groups` on the catalytic graph, an allosteric mechanism's active state). No
kinetic group binds a competitive inhibitor redundantly (`_redundant_copy_groups`). A parent
must obey both because a flip tests only the groups it flips. The split, the dead-end move,
and `_expand_change_allo_state` filter their children; the other moves preserve both rules.
"""
function _assert_emission_rules(c::CandidateMechanism)
    cm = c.mechanism.catalytic
    members = _kinetic_groups(cm)
    # A group as its steps' DSL lines (`step_text`).
    label(g) = join((Mechanisms.step_text(cm, k) for k in members[g]), ", ")
    flux = _flux_carrying_groups(cm, _rapid_equilibrium_flags(cm))
    for (g, ks) in enumerate(members)
        cm.steps[first(ks)].rapid_equilibrium || flux[g] || error(
            "candidate_mechanisms: steady-state kinetic group {" * label(g) * "} has no " *
            "step that carries net flux, so its two constants enter the rate only as " *
            "their ratio; write the group at rapid equilibrium")
    end
    for g in _redundant_copy_groups(c)
        error("candidate_mechanisms: kinetic group {" * label(g) * "} binds a " *
              "competitive inhibitor only where the complex duplicates a productive " *
              "form, and a dwell gauge absorbs its constant into the existing binding's; " *
              "bind the inhibitor where it forms a new complex, or drop it")
    end
    nothing
end

"""
    _assert_atom_conserving(c::CandidateMechanism, rxn::EnzymeReaction)

Assert every edge of the catalytic graph of an enumerated candidate `c` conserves the atoms
that reaction `rxn` gives its metabolites: the atoms of its source form plus those it takes
up equal the atoms of its target form plus those it gives off, a form carrying the atoms of
its bound metabolites and of its residual's added substrates, less those of its residual's
subtracted products, each read from the reaction's inventory by name (`_net_atoms`). Errors
on the first edge that does not, group by group, naming its forms and the metabolites it
exchanges. This is the enumeration-path guardrail against atom-non-conserving covalent
intermediates; it is intentionally NOT a `Mechanism` constructor check, since a hand-written
`@enzyme_mechanism` block carries no atoms and may fold steps.
"""
function _assert_atom_conserving(c::CandidateMechanism, rxn::Mechanisms.EnzymeReaction)
    cm = c.mechanism.catalytic
    names(v::Vector{Int}) = cm.metabolites[v]
    # The names counted on one side of an edge: `form`'s bound metabolites, productive or
    # copies, and residual additions, the residual subtractions of the `other` form, and
    # the free metabolites.
    side(form::Int, other::Int, free::Vector{Int}) =
        [c.forms[form].productive; c.forms[form].copies; c.forms[form].adds;
         c.forms[other].subtracts; names(free)]
    for ks in _kinetic_groups(cm), k in ks
        f, t = cm.steps[k].from, cm.steps[k].to
        taken, given = cm.steps[k].taken_up, cm.steps[k].given_off
        diff = _net_atoms(rxn, side(t, f, given), side(f, t, taken))
        isempty(diff) && continue
        error(_atom_conserving_message(cm, f, t, names(taken), names(given), diff))
    end
    nothing
end

function _atom_conserving_message(cm, f::Int, t::Int, taken, given, diff)
    from, to = cm.enzyme_forms[f], cm.enzyme_forms[t]
    "atom-non-conserving step $from → $to (consumed $taken, released $given): " *
    "atoms(to) + atoms(released) − atoms(from) − atoms(consumed) = $diff"
end

# ─── Regulator-Type Filter ────────────────────────────────────

"""
    _respects_reg_type(m, rxn::EnzymeReaction) -> Bool

Whether every regulatory ligand of `m` respects the type `rxn` declares for it, given the
other ligands at its site. A plain mechanism has no regulatory sites and passes. A ligand's
type is the `effect` of its allosteric `Regulator` entry in `rxn`, or `:unspecified`
when `rxn` has no such entry (the name is absent, or declared only as a competitive
inhibitor). An `:unspecified` ligand always passes. A designated activator is
never `:OnlyI` and a designated inhibitor is never `:OnlyA`; the type-matching pure
state and `:NonequalAI` are always allowed. `:EqualAI` is an antagonist state and is
rejected only when another ligand at the site is a same-type designated effector, since
it would counteract that ligand's declared direction.
"""
function _respects_reg_type(am::Mechanisms.Mechanism, rxn::Mechanisms.EnzymeReaction)
    declared = Dict(reg.name => reg.effect
                    for reg in rxn.regulators if reg.kind === :allosteric)
    all(am.allosteric.regulatory_sites) do site
        types = [get(declared, l, :unspecified) for l in site.ligands]
        all(zip(types, site.allosteric_states)) do (t, state)
            t === :unspecified ||
                !(t === :activator && state === :OnlyI ||
                  t === :inhibitor && state === :OnlyA ||
                  state === :EqualAI && count(==(t), types) > 1)
        end
    end
end

# ─── Limits, kinetic groups and plain bindings ──────────────────────────

"""The candidates of `cands` whose law has at most `eq_complexity_filter` terms and at
most `max_param_count` fitted constants, in order. A limit of `typemax(Int)` limits
nothing, so its count is not computed: the fitted count solves the mechanism."""
_within_size_limits(cands::Vector{CandidateMechanism}, max_param_count::Int,
                    eq_complexity_filter::Int) =
    filter(cands) do c
        (eq_complexity_filter == typemax(Int) ||
         Derivation.rate_equation_term_count(c.mechanism) <= eq_complexity_filter) &&
            (max_param_count == typemax(Int) ||
             Derivation.fitted_parameter_count(c.mechanism) <= max_param_count)
    end

"""Whether the enumerator keeps `m`: the guard (`Derivation.check_valid_rate_equation`)
accepts it, and no conformation keeps a form that only bindings it excludes enter
(`_kept_forms_keep_a_binding`). The test the enumerator's exits and its flip, split and
tagging searches apply to a candidate."""
_enumerator_keeps(m) =
    Derivation.check_valid_rate_equation(Bool, m) && _kept_forms_keep_a_binding(m)

"""The number of copies of metabolite `i` (an index into a scheme's `metabolites`) in the
list `v`."""
_copies(v::Vector{Int}, i::Int) = count(==(i), v)

"""Whether `m` is allosteric: its kinetic groups carry allosteric states."""
_is_allosteric(m::Mechanisms.Mechanism) = !isempty(m.allosteric.allosteric_states)

"""The rapid-equilibrium flag of each step of `m`, in step order."""
_rapid_equilibrium_flags(m::Mechanisms.CatalyticScheme) =
    BitVector([e.rapid_equilibrium for e in m.steps])

"""Each kinetic group of graph `m`, its steps in step order."""
function _kinetic_groups(m::Mechanisms.CatalyticScheme)
    members = [Int[] for _ in 1:maximum(e -> e.kinetic_group, m.steps; init = 0)]
    for (k, e) in enumerate(m.steps)
        push!(members[e.kinetic_group], k)
    end
    members
end

"""How many more modifier copies of `name` the record of candidate `c` holds at the second
form of step `k` than at its first: positive when the step takes copies of `name` up,
negative when it gives them off, zero when it exchanges `name` only productively. A copy
binds a site of its own and takes no part in the chemistry, so it enters or leaves a form
only with the step that takes it up or gives it off."""
function _copy_change(c::CandidateMechanism, k::Int, name::Symbol)
    e = c.mechanism.catalytic.steps[k]
    count(==(name), c.forms[e.to].copies) - count(==(name), c.forms[e.from].copies)
end

"""Whether step `k` of candidate `c` is a plain binding: it takes up one metabolite and
gives off none (`_is_binding`), and the record of its second form is that of its first
with the metabolite added, among the copies when the step binds a copy (`_copy_change`),
the residual the same and the label free to change. Every other step is chemistry: an
isomerization, a fused binding (the last substrate binding straight into the
product-bound form), a Theorell–Chance step or an exchange."""
function _is_plain_binding(c::CandidateMechanism, k::Int)
    cm = c.mechanism.catalytic
    e = cm.steps[k]
    _is_binding(e) || return false
    f, t, x = c.forms[e.from], c.forms[e.to], cm.metabolites[only(e.taken_up)]
    f.adds == t.adds && f.subtracts == t.subtracts || return false
    _copy_change(c, k, x) > 0 ?
        f.productive == t.productive && _one_more(t.copies, f.copies, x) :
        f.copies == t.copies && _one_more(t.productive, f.productive, x)
end

"""Whether the sorted list `t` is the sorted list `f` with one more `x`."""
function _one_more(t::Vector{Symbol}, f::Vector{Symbol}, x::Symbol)
    length(t) == length(f) + 1 || return false
    j, skipped = 1, false
    for y in t
        if j <= length(f) && y === f[j]
            j += 1
        elseif !skipped && y === x
            skipped = true
        else
            return false
        end
    end
    skipped
end

# ─── Steps, subsets, merged complexes and bound regulators ──────────────

"""Whether step `s` takes up and gives off nothing: an isomerization."""
_is_iso(s::Mechanisms.Step) = isempty(s.taken_up) && isempty(s.given_off)

"""Whether step `s` takes up one metabolite and gives off none: a binding, plain or fused
(`_is_plain_binding`)."""
_is_binding(s::Mechanisms.Step) = isempty(s.given_off) && length(s.taken_up) == 1

"""All subsets of `v` in binary-counting order (bit `i - 1` selects `v[i]`), empty first."""
_subsets(v::AbstractVector) =
    [v[[isodd(mask >> (i - 1)) for i in eachindex(v)]] for mask in 0:(1 << length(v)) - 1]

"""The merged complexes of the catalytic graph of candidate `c`, each once: the forms that a
fused binding of a substrate enters, a binding that is not plain (`_is_binding`,
`_is_plain_binding`) and takes up a metabolite of `subs`, indices into the scheme's
metabolites. Only the edges `on` flags count."""
function _merged_complexes(c::CandidateMechanism, subs::Vector{Int},
                           on::BitVector = trues(length(
                               c.mechanism.catalytic.steps)))
    steps = c.mechanism.catalytic.steps
    unique(steps[k].to for k in eachindex(steps)
           if on[k] && _is_binding(steps[k]) && !_is_plain_binding(c, k) &&
              _holds(steps[k].taken_up, subs))
end

"""Names of the allosteric regulators bound at `m`'s regulatory sites (empty for a plain
mechanism)."""
_bound_allo_regs(m::Mechanisms.Mechanism) =
    Set{Symbol}(lig for site in m.allosteric.regulatory_sites for lig in site.ligands)

"""Names of the competitive inhibitors some step of candidate `c` takes up as a modifier
copy (`_copy_change`)."""
function _bound_comp_inhibitors(c::CandidateMechanism)
    cm = c.mechanism.catalytic
    Set{Symbol}(cm.metabolites[x] for k in eachindex(cm.steps) for x in cm.steps[k].taken_up
                if _copy_change(c, k, cm.metabolites[x]) > 0)
end

# ─── Kinetic groups that carry flux ─────────────────────────────────────

"""Whether each kinetic group of graph `m` holds some step `k` with `pred(k)`."""
function _groups_with(pred, m::Mechanisms.CatalyticScheme)
    out = falses(maximum(e -> e.kinetic_group, m.steps; init = 0))
    for k in eachindex(m.steps)
        pred(k) && (out[m.steps[k].kinetic_group] = true)
    end
    out
end

"""
    _flux_carrying_edges(m, re) -> BitVector

Whether each edge of graph `m`, edge `k` at rapid equilibrium where `re[k]` holds,
carries net steady-state flux for generic parameter values. Rapid-equilibrium edges hold
equilibrium mass and are flagged `false`; the rule that uses these flags concerns
steady-state groups only.

Contract each rapid-equilibrium segment to one vertex (`rapid_equilibrium_segments`) and
make every steady-state edge an edge between its two forms' segments, weighted by the
substrates it takes up minus the products (`_uptake_weight`), plus its first end's offsets
minus its second end's (`from`, `to`; a modifier copy counts as its reactant). Around any
cycle the offsets telescope and the uptakes sum to the net turnover times the reactant
count, so a cycle has weight zero exactly when it runs no net reaction. A self-loop (both
forms in one segment) carries flux iff its weight is nonzero. Any other edge carries flux
iff its biconnected block (`_unbalanced_blocks`) holds a cycle of nonzero weight: every edge
of such a block lies on one (join the edge to the cycle by two disjoint paths; one of the
two resulting cycles is unbalanced), while in a balanced block detailed balance holds along
every cycle and each edge's flux vanishes. The zero-flux verdict holds for any parameters
and grouping; the flux-carrying verdict for one-way steps, which every step the enumerator
emits is. The test reads each edge's metabolite lists, so it holds for fused and
Theorell–Chance steps as well.
"""
function _flux_carrying_edges(m::Mechanisms.CatalyticScheme, re::BitVector)
    s = Mechanisms.rapid_equilibrium_segments(m, trues(length(m.steps)), re)
    rho = _name_signs(m)
    offset(f) = sum(rho[n] * s.offsets[n, f] for n in eachindex(rho); init = 0)
    flags = falses(length(m.steps))
    edges = Tuple{Int, Int}[]; weights = Int[]; owner = Int[]
    for k in eachindex(m.steps)
        re[k] && continue
        (; from, to) = m.steps[k]
        weight = _uptake_weight(m, k) + offset(from) - offset(to)
        a, b = s.segment_of[from], s.segment_of[to]
        if a == b
            flags[k] = weight != 0
        else
            push!(edges, (a, b)); push!(weights, weight); push!(owner, k)
        end
    end
    block, unbalanced = _unbalanced_blocks(length(s.segments), edges, weights)
    for (j, k) in enumerate(owner)
        flags[k] = block[j] in unbalanced
    end
    flags
end

"""
    _flux_carrying_groups(m, re) -> BitVector

One flag per kinetic group of graph `m`, edge `k` at rapid equilibrium where `re[k]`
holds: whether some edge of the group carries net steady-state flux
(`_flux_carrying_edges`). A steady-state group with no such edge exposes only the ratio
of its two constants, so the moves never emit one; a zero-flux edge inside a group that
also holds a flux-carrying edge costs nothing, because the group's shared constants are
pinned by the edge that carries flux.
"""
function _flux_carrying_groups(m::Mechanisms.CatalyticScheme, re::BitVector)
    flux = _flux_carrying_edges(m, re)
    _groups_with(k -> flux[k], m)
end

"""Each metabolite's net uptake sign, in the order of `m.metabolites` (the rows of the
offsets of `rapid_equilibrium_segments`): +1 for a substrate, −1 for a product, 0 for a
modifier."""
_name_signs(m::Mechanisms.CatalyticScheme) =
    Int[n in m.substrates ? 1 : n in m.products ? -1 : 0 for n in eachindex(m.metabolites)]

"""The substrates step `k` of `m` takes up net, less the products it takes up net."""
function _uptake_weight(m::Mechanisms.CatalyticScheme, k::Int)
    t, g = m.steps[k].taken_up, m.steps[k].given_off
    net(x) = _copies(t, x) - _copies(g, x)
    sum(net, m.substrates; init = 0) -
        sum(net, m.products; init = 0)
end

"""
    _unbalanced_blocks(nv, edges, weights) -> (block, unbalanced)

The biconnected block of every edge of an undirected multigraph without self-loops,
`edges[e] = (u, v)` with vertices `1:nv` (Tarjan's edge-stack algorithm; a bridge is a
block of its own), and the set of blocks holding a cycle of nonzero weight, edge `e`
weighing `weights[e]` along its stored direction. The search gives each vertex a
potential that rises by each tree edge's weight along its stored direction. An edge's
defect, its weight minus the potential difference of its ends, is zero on a tree edge
and, on any other edge, the weight of the cycle it closes with the tree path between its
ends, a cycle of its own block. The tree edges of a block span it, so these cycles
generate the block's cycles and weight adds over them: a block holds a cycle of nonzero
weight exactly when one of its edges has a nonzero defect.
"""
function _unbalanced_blocks(nv::Int, edges::Vector{Tuple{Int, Int}}, weights::Vector{Int})
    adj = [Int[] for _ in 1:nv]
    for (e, (u, v)) in enumerate(edges)
        push!(adj[u], e); push!(adj[v], e)
    end
    disc = zeros(Int, nv); low = zeros(Int, nv); phi = zeros(Int, nv)
    block = zeros(Int, length(edges))
    stack = Int[]; clock = Ref(0); nblocks = Ref(0)
    # It takes itself as an argument because a closure that names itself is boxed.
    function visit(self, u, parent_edge)
        clock[] += 1; disc[u] = low[u] = clock[]
        for e in adj[u]
            e == parent_edge && continue
            w = edges[e][1] == u ? edges[e][2] : edges[e][1]
            if disc[w] == 0
                phi[w] = phi[u] + (edges[e][1] == u ? weights[e] : -weights[e])
                push!(stack, e)
                self(self, w, e)
                low[u] = min(low[u], low[w])
                if low[w] >= disc[u]
                    nblocks[] += 1
                    while true
                        x = pop!(stack); block[x] = nblocks[]
                        x == e && break
                    end
                end
            elseif disc[w] < disc[u]
                push!(stack, e)
                low[u] = min(low[u], disc[w])
            end
        end
    end
    for v in 1:nv
        disc[v] == 0 && visit(visit, v, 0)
    end
    block, Set{Int}(block[e] for (e, (u, v)) in enumerate(edges)
                    if weights[e] != phi[v] - phi[u])
end
