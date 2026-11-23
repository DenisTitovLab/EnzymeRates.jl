# ABOUTME: Mechanism, a reaction graph in canonical form and its allostery, with their
# ABOUTME: checks, identity and text, and what the graph gives: steps, segments, cycles.

# ─── Regulatory sites ───────────────────────────────────────────────────

"""The allosteric states a kinetic group or a regulatory site's ligand takes; the
constructors accept no other (`_check_tags`, `RegulatorySite`)."""
const _ALLOSTERIC_STATES = (:OnlyA, :OnlyI, :EqualAI, :NonequalAI)

"""
RegulatorySite: a binding site (possibly multimeric) for one or
more allosteric ligands, each given by its name, once. `ligands[i]` and
`allosteric_states[i]` are parallel. The constructor sorts the ligands by name and
carries each ligand's state along, so the order the ligands are given in never changes
`==`, `hash` or the compiled rate law.
"""
struct RegulatorySite
    ligands::Vector{Symbol}
    multiplicity::Int
    allosteric_states::Vector{Symbol}
    function RegulatorySite(ligands::Vector{Symbol},
                            multiplicity::Int, allosteric_states::Vector{Symbol})
        isempty(ligands) && error("RegulatorySite: needs at least one ligand")
        length(ligands) == length(allosteric_states) ||
            error("RegulatorySite: length(ligands)=$(length(ligands)) " *
                  "must equal length(allosteric_states)=$(length(allosteric_states))")
        multiplicity ≥ 1 ||
            error("RegulatorySite: multiplicity must be ≥ 1, got $multiplicity")
        bad = findfirst(∉(_ALLOSTERIC_STATES), allosteric_states)
        bad === nothing ||
            error("RegulatorySite: allosteric state $(allosteric_states[bad]) " *
                  "must be one of $_ALLOSTERIC_STATES")
        dup = findfirst(i -> ligands[i] in view(ligands, 1:(i - 1)), eachindex(ligands))
        dup === nothing ||
            error("RegulatorySite: allosteric regulator `$(ligands[dup])` is listed " *
                  "twice; an allosteric regulator sits at exactly one regulatory site.")
        p = sortperm(ligands)
        new(ligands[p], multiplicity, allosteric_states[p])
    end
end

# ─── The reaction graph ─────────────────────────────────────────────────

"""
    Step(from, to, taken_up, given_off, rapid_equilibrium, kinetic_group)

A step between the enzyme forms `from` and `to` of its catalytic scheme, indices into its
`enzyme_forms`: going `from → to` it takes up the metabolites `taken_up` and gives off
`given_off`, indices into its `metabolites`, each list sorted and holding a metabolite as
often as the step moves it. `rapid_equilibrium` marks a rapid-equilibrium step, and
`kinetic_group` is the step's kinetic group, whose steps share their constants.
"""
struct Step
    from::Int
    to::Int
    taken_up::Vector{Int}
    given_off::Vector{Int}
    rapid_equilibrium::Bool
    kinetic_group::Int
end

"""
    CatalyticScheme(metabolites, enzyme_forms, substrates, products, steps)

The catalytic scheme of a mechanism: a reaction graph. `metabolites` names the metabolites
its steps take up or give off, and `substrates` and `products` index them; every other
metabolite is a modifier. `enzyme_forms` names its enzyme forms, the free enzyme `E` among
them; a name is opaque, and only `E` means something. Each step of `steps` joins two forms
and takes up and gives off metabolites by index (`Step`); its `kinetic_group` labels its
kinetic group, the labels running 1, 2, … with each used. The lists may come in any order
and a step may be written in either direction.

The constructor checks the representation and stores the scheme in canonical form, all
but its kinetic group labels: the metabolites as substrates, then products, then
modifiers, each by name; `E` first and the other forms by name; each step in the direction
the graph gives it (`_orientation`), the steps sorted. Each step keeps the kinetic group
label it is written with. Only [`Mechanism`](@ref), which reads the labels to give each
group its allosteric state, numbers the groups canonically, by their first step, so two
schemes of one graph can differ (`!=`) where the mechanisms built on them are equal.

It throws an `ArgumentError` when no form is named `E`, and an `ErrorException` when:
- an index is out of range, or the substrates or the products are empty;
- a name is listed twice, a substrate is also a product, or a modifier is also a reactant;
- a metabolite's name is one `_check_metabolite_names` rejects: `E`, a reserved name such
  as `Keq`, a name holding `_`, or one no DSL text can write;
- a form's name is one no DSL text can write, or a form is named like a metabolite;
- the kinetic group labels do not run from 1 with each used;
- a step joins a form to itself, or takes up and gives off one metabolite;
- a form is not joined to `E`;
- a kinetic group mixes kinds of step, or a reaction appears twice;
- a metabolite takes part in no step.
"""
struct CatalyticScheme
    metabolites::Vector{Symbol}
    enzyme_forms::Vector{Symbol}
    substrates::Vector{Int}
    products::Vector{Int}
    steps::Vector{Step}
    function CatalyticScheme(metabolites::Vector{Symbol}, enzyme_forms::Vector{Symbol},
                             substrates::Vector{Int}, products::Vector{Int},
                             steps::Vector{Step})
        _check_written(metabolites, enzyme_forms, substrates, products, steps)
        c = _canonical(metabolites, enzyme_forms, substrates, products, steps)
        _check_graph(c)
        c
    end
    # The scheme with these fields, unchecked and as given: the canonical scheme is built
    # through it, and a mechanism renumbers its kinetic groups through it.
    global _unchecked_scheme(metabolites::Vector{Symbol}, enzyme_forms::Vector{Symbol},
                             substrates::Vector{Int}, products::Vector{Int},
                             steps::Vector{Step}) =
        new(metabolites, enzyme_forms, substrates, products, steps)
end

# ─── Canonical form ─────────────────────────────────────────────────────

"""
    _canonical(metabolites, enzyme_forms, substrates, products, steps) -> CatalyticScheme

The canonical scheme that the arguments write, in the order `CatalyticScheme` gives,
unchecked but for what it needs: indices in range and a form named `E` (`_check_written`).
Each step's lists are sorted, and the steps sorted by their ends, then their lists, then
their flag.
"""
function _canonical(metabolites::Vector{Symbol}, enzyme_forms::Vector{Symbol},
                    substrates::Vector{Int}, products::Vector{Int}, steps::Vector{Step})
    role(i) = i in substrates ? 1 : i in products ? 2 : 3
    met_order = sortperm(eachindex(metabolites); by = i -> (role(i), metabolites[i]))
    met_of = invperm(met_order)
    form_order = sortperm(eachindex(enzyme_forms);
                          by = f -> (enzyme_forms[f] !== :E, enzyme_forms[f]))
    form_of = invperm(form_order)
    ns, np = length(substrates), length(products)
    written = [Step(form_of[e.from], form_of[e.to], sort!(met_of[e.taken_up]),
                    sort!(met_of[e.given_off]), e.rapid_equilibrium, e.kinetic_group)
               for e in steps]
    fwd = _orientation(written, length(enzyme_forms), ns, np)
    stored = [fwd[k] ? e : _reversed(e) for (k, e) in enumerate(written)]
    sort!(stored; by = e -> (e.from, e.to, e.taken_up, e.given_off, e.rapid_equilibrium))
    _unchecked_scheme(metabolites[met_order], enzyme_forms[form_order], collect(1:ns),
                      collect((ns + 1):(ns + np)), stored)
end

_reversed(e::Step) =
    Step(e.to, e.from, e.given_off, e.taken_up, e.rapid_equilibrium, e.kinetic_group)

_is_iso(e::Step) = isempty(e.taken_up) && isempty(e.given_off)

"""
    _orientation(written, nf, ns, np) -> BitVector

Whether to store each step of `written` as written, from the graph alone: its `nf` forms
stored in canonical order, `E` first and the others by name, and its metabolites the
`ns` substrates, then the `np` products, then the modifiers. The first rule that decides
gives a step's direction:

1. A step that only takes up, or only gives off, metabolites runs as the uptake.
2. A step that takes up some and gives off others runs the way that takes up more
   substrates and gives off more products, net.
3. Any other step, an isomerization or one that trades one modifier for another, runs
   the way the shortest cycle through it that turns the reaction over runs
   (`_turnover_directions!`).
4. With no such cycle, or one each way of equal length, it runs away from `E`, from the
   form fewer steps from it.
5. With both ends equally far, it runs from the form stored first.

The steps of a kinetic group share their constants, so they turn together, after the
member nearest `E` among those rules 1 to 3 decide, or among all of them when they decide
none (`_turn_groups_together!`).
"""
function _orientation(written::Vector{Step}, nf::Int, ns::Int, np::Int)
    net(list::Vector{Int}) = count(<=(ns), list) - count(i -> ns < i <= ns + np, list)
    # Rules 1 and 2: positive to store a step as written, negative to turn it, zero when
    # neither decides, as for an isomerization.
    verdict = [isempty(e.taken_up) || isempty(e.given_off) ?
               length(e.taken_up) - length(e.given_off) :
               net(e.taken_up) - net(e.given_off) for e in written]
    at = _steps_at_forms(written, nf)
    # The first substrate, metabolite 1, counts turnover.
    _turnover_directions!(verdict, written, at, 1)
    distance = _distance_from_free_E(written, at)
    fwd, decided = trues(length(written)), trues(length(written))
    for (k, e) in enumerate(written)
        decided[k] = verdict[k] != 0
        a, b = distance[e.from], distance[e.to]
        fwd[k] = decided[k] ? verdict[k] > 0 : a == b ? e.from < e.to : a < b
    end
    _turn_groups_together!(fwd, decided, written, distance)
    fwd
end

"""The steps of `written` at each of its `nf` forms: those that touch form `u` are
`at.steps[at.start[u]:(at.start[u + 1] - 1)]`, a step joining a form to itself twice."""
function _steps_at_forms(written::Vector{Step}, nf::Int)
    start = zeros(Int, nf + 1)
    for e in written
        start[e.from + 1] += 1
        start[e.to + 1] += 1
    end
    start[1] = 1
    for u in 1:nf
        start[u + 1] += start[u]
    end
    steps, filled = Vector{Int}(undef, 2 * length(written)), start[1:nf]
    for (k, e) in enumerate(written), u in (e.from, e.to)
        steps[filled[u]] = k
        filled[u] += 1
    end
    (; start, steps)
end

_other_end(e::Step, u::Int) = e.from == u ? e.to : e.from

"""How many steps of `written` separate each of its forms from `E`, form 1, along
steps taken either way (`at`, from `_steps_at_forms`); `typemax(Int)` for a form no
steps join to it."""
function _distance_from_free_E(written::Vector{Step}, at)
    distance = fill(typemax(Int), length(at.start) - 1)
    isempty(distance) && return distance
    distance[1] = 0
    queue, head = [1], 0
    while head < length(queue)
        u = queue[head += 1]
        for i in at.start[u]:(at.start[u + 1] - 1)
            v = _other_end(written[at.steps[i]], u)
            distance[v] == typemax(Int) || continue
            distance[v] = distance[u] + 1
            push!(queue, v)
        end
    end
    distance
end

"""
    _turnover_directions!(verdict, written, at, s1) -> verdict

Decide each step `k` of `written` (`at` lists its steps by form, `_steps_at_forms`) that
`verdict` leaves at `0`: `verdict[k]` becomes `1` when the shortest closed walk that takes
step `k` once, as written, and takes up one of the first substrate, metabolite `s1`, net
is shorter than any that takes the step the other way, `-1` when it
is longer, and stays `0` when there is none either way or both are equally long. The walk
may pass a form twice, so a step on a cycle that turns nothing over finds a walk each way
of one length, round a turnover cycle in either direction, and is decided by neither.
"""
function _turnover_directions!(verdict::Vector{Int}, written::Vector{Step}, at, s1::Int)
    any(iszero, verdict) || return verdict
    nf = length(at.start) - 1
    uptake = [count(==(s1), e.taken_up) - count(==(s1), e.given_off) for e in written]
    reached, forms = falses(nf), Int[]
    # The search runs breadth first over (form, net uptake of s1) and keeps the net uptake
    # within ±span, as many as there are forms, so that it ends: a walk round a turnover
    # cycle takes up more on every turn. Each queued state carries its walk's length.
    span = nf
    seen = falses(nf, 2span + 1)
    queue = NTuple{3, Int}[]
    for k in eachindex(written)
        verdict[k] == 0 || continue
        from, to = written[k].from, written[k].to
        # A step whose two ends the other steps do not join lies on no cycle. Looking
        # first keeps the search below from visiting every (form, net uptake) pair for
        # such a step.
        fill!(reached, false)
        empty!(forms)
        reached[to] = true
        push!(forms, to)
        head = 0
        while head < length(forms) && !reached[from]
            u = forms[head += 1]
            for i in at.start[u]:(at.start[u + 1] - 1)
                at.steps[i] == k && continue
                v = _other_end(written[at.steps[i]], u)
                reached[v] || (reached[v] = true; push!(forms, v))
            end
        end
        reached[from] || continue
        # The walk back from `to` to `from`, without step `k`, takes up `1 - uptake[k]`
        # net when step `k` as written turns the reaction over, and `-1 - uptake[k]` when
        # step `k` taken the other way does.
        ahead, behind = 1 - uptake[k], -1 - uptake[k]
        abs(ahead) <= span && abs(behind) <= span || continue
        fill!(seen, false)
        empty!(queue)
        seen[to, span + 1] = true
        push!(queue, (to, 0, 0))
        a = b = typemax(Int)
        head = 0
        # Once the search reaches a state as far out as the shortest walk found, every
        # walk of that length has been found.
        while head < length(queue) && queue[head + 1][3] < min(a, b)
            u, n, d = queue[head += 1]
            for i in at.start[u]:(at.start[u + 1] - 1)
                j = at.steps[i]
                j == k && continue
                v = _other_end(written[j], u)
                x = written[j].from == u ? uptake[j] : -uptake[j]
                (abs(n + x) > span || seen[v, n + x + span + 1]) && continue
                seen[v, n + x + span + 1] = true
                v == from && n + x == ahead && (a = d + 1)
                v == from && n + x == behind && (b = d + 1)
                push!(queue, (v, n + x, d + 1))
            end
        end
        verdict[k] = a == b ? 0 : a < b ? 1 : -1
    end
    verdict
end

"""
Turn the steps of each kinetic group (each written step's `kinetic_group` label) with the
group's lead. The lead is its member nearest `E` among those `decided` marks, or among all
of them when none is marked: nearest by the distances of its two ends from `E`
(`distance`), the nearer end first, then by its two forms, the one stored first first. A
member that takes up or gives off metabolites takes the direction in which it takes up
what the lead takes up; an isomerization keeps the direction it was written in relative to
the lead. `fwd[k]` says whether written step `k` is stored as written.
"""
function _turn_groups_together!(fwd::BitVector, decided::BitVector, written::Vector{Step},
                                distance::Vector{Int})
    near(k::Int) = (e = written[k];
                    (minmax(distance[e.from], distance[e.to])..., minmax(e.from, e.to)...))
    for g in 1:maximum(e -> e.kinetic_group, written; init = 0)
        in_group(k::Int) = written[k].kinetic_group == g
        count(in_group, eachindex(written)) < 2 && continue
        any_decided = any(k -> in_group(k) && decided[k], eachindex(written))
        lead = 0
        for k in eachindex(written)
            in_group(k) && (decided[k] || !any_decided) || continue
            (lead == 0 || near(k) < near(lead)) && (lead = k)
        end
        taken_up = fwd[lead] ? written[lead].taken_up : written[lead].given_off
        for k in eachindex(written)
            in_group(k) || continue
            e = written[k]
            if _is_iso(e)
                fwd[k] = fwd[lead]
            elseif e.taken_up == taken_up || e.given_off == taken_up
                fwd[k] = e.taken_up == taken_up
            end
        end
    end
    fwd
end

# ─── Checks ─────────────────────────────────────────────────────────────

"""
Error unless the scheme as written can be stored:
- every index in range, and `substrates` and `products` not empty;
- each metabolite name once, a substrate's not a product's and a modifier's not a
  reactant's, and each one a name a metabolite can take (`_check_metabolite_names`);
- the form names distinct, names a text can write (`_check_identifiers`) and no
  metabolite's;
- a form named `E` (an `ArgumentError`);
- the kinetic group labels running from 1 with each used.
"""
function _check_written(metabolites::Vector{Symbol}, enzyme_forms::Vector{Symbol},
                        substrates::Vector{Int}, products::Vector{Int},
                        steps::Vector{Step})
    nm, nf = length(metabolites), length(enzyme_forms)
    for (k, e) in enumerate(steps)
        e.from in 1:nf && e.to in 1:nf ||
            error("Mechanism: step $k joins forms $(e.from) and $(e.to), but there are " *
                  "$nf forms; a step's ends index the enzyme forms")
        all(in(1:nm), e.taken_up) && all(in(1:nm), e.given_off) ||
            error("Mechanism: step $k takes up $(e.taken_up) and gives off " *
                  "$(e.given_off), but there are $nm metabolites; a step's lists index " *
                  "the metabolites")
    end
    for (list, kind) in ((substrates, "substrate"), (products, "product"))
        isempty(list) && error("Mechanism: $(kind)s must not be empty")
        all(in(1:nm), list) ||
            error("Mechanism: $(kind)s $list, but there are $nm metabolites; the " *
                  "substrates and products index the metabolites")
        allunique(metabolites[list]) || error("Mechanism: duplicate $kind names")
    end
    both = findfirst(in(metabolites[products]), metabolites[substrates])
    both === nothing ||
        error("Mechanism: $(metabolites[substrates[both]]) named as both a substrate " *
              "and a product; concentrations and constants are keyed by name")
    modifiers = metabolites[[i for i in 1:nm if !(i in substrates || i in products)]]
    allunique(modifiers) || error("Mechanism: duplicate modifier names")
    reactant = findfirst(in(metabolites[[substrates; products]]), modifiers)
    reactant === nothing ||
        error("Mechanism: $(modifiers[reactant]) is declared both as a reactant and as a " *
              "modifier; a step that binds a reactant as a modifier simply takes it up, " *
              "as `ES + S ⇌ ESS` does")
    _check_metabolite_names("Mechanism", metabolites)
    _check_identifiers("Mechanism", "enzyme form", enzyme_forms)
    dup = findfirst(i -> enzyme_forms[i] in view(enzyme_forms, 1:(i - 1)),
                    eachindex(enzyme_forms))
    dup === nothing || error("Mechanism: two enzyme forms are named $(enzyme_forms[dup])")
    shared = findfirst(in(metabolites), enzyme_forms)
    shared === nothing ||
        error("Mechanism: $(enzyme_forms[shared]) names both an enzyme form and a " *
              "metabolite; a text could not tell them apart")
    :E in enzyme_forms ||
        throw(ArgumentError("Mechanism: no enzyme form is named `E`; the free enzyme is " *
                            "the form named `E`"))
    labels = Int[e.kinetic_group for e in steps]
    all(≥(1), labels) ||
        error("Mechanism: each step takes one kinetic group label from 1 to the number " *
              "of groups")
    empty = findfirst(∉(labels), 1:maximum(labels; init = 0))
    empty === nothing ||
        error("Mechanism: kinetic group $empty is empty; a kinetic group holds at least " *
              "one step")
    nothing
end

"""
Error unless the canonical scheme `c` is a graph a mechanism can have: the checks of
`CatalyticScheme` that read its steps.
"""
function _check_graph(c::CatalyticScheme)
    for e in c.steps
        e.from == e.to &&
            error("Step: both ends are the form $(c.enzyme_forms[e.from]); a step joins " *
                  "two forms")
        shared = unique!([i for i in e.taken_up if i in e.given_off])
        isempty(shared) ||
            error("Step $(c.enzyme_forms[e.from]) → $(c.enzyme_forms[e.to]): " *
                  "$(join(c.metabolites[shared], ", ")) both consumed and released")
    end
    cut_off = findall(!, _joined_to_free_E(c, falses(length(c.steps))))
    isempty(cut_off) ||
        error("Mechanism: the steps do not form one connected graph; no steps lead from " *
              "free E to " * join(c.enzyme_forms[cut_off], ", "))
    # The checks name each kinetic group by the number a mechanism gives it, by its first
    # step; the scheme keeps the written labels, which `Mechanism` reads.
    written = Int[e.kinetic_group for e in c.steps]
    labels, _ = _number_groups(written, maximum(written; init = 0))
    groups = [findall(==(g), labels) for g in 1:maximum(labels; init = 0)]
    _check_uniform_groups(c, groups)
    _check_each_reaction_once(c, groups)
    used = falses(length(c.metabolites))
    for e in c.steps, i in Iterators.flatten((e.taken_up, e.given_off))
        used[i] = true
    end
    k = findfirst(!, used)
    k === nothing && return nothing
    k in c.substrates || k in c.products ||
        error("Mechanism: modifier $(c.metabolites[k]) is declared, but no step binds it")
    error("Mechanism: $(k in c.substrates ? "substrate" : "product") " *
          "$(c.metabolites[k]) is declared, but no step uses it")
end

"""Error unless the steps of each kinetic group (`groups`) of `c` take up and give off the
same metabolites and carry one flag: they share their constants."""
function _check_uniform_groups(c::CatalyticScheme, groups::Vector{Vector{Int}})
    for members in groups
        e1 = c.steps[first(members)]
        for k in members
            e = c.steps[k]
            (e.taken_up, e.given_off, e.rapid_equilibrium) ==
                (e1.taken_up, e1.given_off, e1.rapid_equilibrium) && continue
            error("Mechanism: a kinetic group holds $(step_text(c, first(members))) and " *
                  "$(step_text(c, k)), which differ in kind or in RE/SS; the steps of a " *
                  "kinetic group share their constants, so they must be the same kind of " *
                  "step with the same flag")
        end
    end
end

"""Error when one reaction, the unordered pair of a step's two sides, appears in two
steps: a reaction belongs to one kinetic group, whose constants are named after it, and
one group holding it twice would count its step twice. Reads the groups (`groups`) in
order, each in step order, and names the reaction as the step seen first stores it."""
function _check_each_reaction_once(c::CatalyticScheme, groups::Vector{Vector{Int}})
    seen = Dict{NTuple{2, Tuple{Int, Vector{Int}}}, Tuple{Int, Int}}()
    flag(e::Step) = e.rapid_equilibrium ? "RE" : "SS"
    for (g, members) in enumerate(groups), k in members
        e = c.steps[k]
        key = minmax((e.from, e.taken_up), (e.to, e.given_off))
        if haskey(seen, key)
            g0, k0 = seen[key]
            e0 = c.steps[k0]
            text = step_text(c, k0)
            g0 == g && error("Mechanism: kinetic group $g holds the reaction $text " *
                             "twice; a mechanism holds each reaction once")
            flags = flag(e0) == flag(e) ? "" :
                " ($(flag(e0)) in group $g0, $(flag(e)) in group $g)"
            error("Mechanism: kinetic groups $g0 and $g both hold the reaction " *
                  "$text$flags; a reaction belongs to one kinetic group, whose " *
                  "constants are named after it")
        end
        seen[key] = (g, k)
    end
end

"""
    Allostery(allosteric_states, oligomeric_state, regulatory_sites)

What a mechanism adds to its catalytic scheme as an MWC enzyme, whose subunits switch
together between an active and an inactive conformation: `allosteric_states[g]`, the
allosteric state of kinetic group `g`; the `oligomeric_state`, the number of subunits that
switch together; and the `regulatory_sites`. A plain mechanism has no allosteric state,
oligomeric state 1 and no regulatory site: `Allostery(Symbol[], 1, RegulatorySite[])`.
[`Mechanism`](@ref) checks it against its catalytic scheme.

A kinetic group's allosteric state says which conformations run its steps, and with which
constants:
- `:EqualAI`: both, with the same constants;
- `:NonequalAI`: both, each with its own constants;
- `:OnlyA`: the active one alone; the inactive conformation excludes the steps;
- `:OnlyI`: the inactive one alone; the active conformation excludes them.

An exclusion reaches beyond its group. The conformations interconvert at free E alone, so
a form that the steps a conformation does not exclude no longer join to free E holds no
enzyme of that conformation, and a step between two such forms does not run in it,
whatever its state. With `E + S ⇌ ES` and `EP ⇌ E + P` both `:OnlyA`, the inactive
conformation has neither `ES` nor `EP`, and `ES <--> EP` does not run in it even as
`:EqualAI`.
`steps_running_in_allosteric_conformations` gives the steps each conformation runs.
"""
struct Allostery
    allosteric_states::Vector{Symbol}
    oligomeric_state::Int
    regulatory_sites::Vector{RegulatorySite}
    # Typed arguments, so that no `Bool` converts to an oligomeric state.
    Allostery(allosteric_states::Vector{Symbol}, oligomeric_state::Int,
              regulatory_sites::Vector{RegulatorySite}) =
        new(allosteric_states, oligomeric_state, regulatory_sites)
end

"""
    Mechanism(catalytic::CatalyticScheme, allosteric::Allostery)

An enzyme mechanism: its catalytic scheme (`catalytic`), a reaction graph, and its
allostery (`allosteric`). It is plain when `allosteric.allosteric_states` is empty, and
allosteric (MWC) otherwise. `@enzyme_mechanism` builds one from its text, and `show`
prints that text.

`allosteric.allosteric_states[l]` is the allosteric state of the steps of `catalytic`
whose kinetic group label is `l`; `Mechanisms.Allostery` says what a state means for the
steps each conformation runs. The constructor stores the mechanism in canonical form,
so two writings of one mechanism are `==` and hash alike: the kinetic groups numbered by
their first step, the states moved with them, and the regulatory sites sorted. A plain
mechanism and an allosteric one with the same catalytic scheme differ. It throws an
`ErrorException` when:
- an allosteric mechanism does not give each kinetic group one known state, or has an
  oligomeric state below 1;
- a plain one has an oligomeric state other than 1 or a regulatory site;
- a ligand cannot name a metabolite, names an enzyme form or sits at two sites.
Whether its rate law is defined is `check_valid_rate_equation`'s question, not the
constructor's.
"""
struct Mechanism
    catalytic::CatalyticScheme
    allosteric::Allostery
    function Mechanism(catalytic::CatalyticScheme, allosteric::Allostery)
        states = allosteric.allosteric_states
        if isempty(states)
            n = length(allosteric.regulatory_sites)
            allosteric.oligomeric_state == 1 && n == 0 ||
                error("Mechanism: a mechanism with no allosteric states has oligomeric " *
                      "state 1 and no regulatory site; got oligomeric state " *
                      "$(allosteric.oligomeric_state) and $n regulatory " *
                      (n == 1 ? "site" : "sites"))
        else
            _check_tags(states, maximum(e -> e.kinetic_group, catalytic.steps; init = 0),
                        allosteric.oligomeric_state)
        end
        sites = sort(allosteric.regulatory_sites; by = _site_key)
        _check_ligands(sites, catalytic)
        _labels_in_order(catalytic.steps) &&
            return new(catalytic,
                       Allostery(copy(states), allosteric.oligomeric_state, sites))
        group, group_of = _number_groups(Int[e.kinetic_group for e in catalytic.steps],
                                         maximum(e -> e.kinetic_group, catalytic.steps))
        steps = [Step(e.from, e.to, e.taken_up, e.given_off, e.rapid_equilibrium, group[k])
                 for (k, e) in enumerate(catalytic.steps)]
        new(_unchecked_scheme(catalytic.metabolites, catalytic.enzyme_forms,
                              catalytic.substrates, catalytic.products, steps),
            Allostery(isempty(states) ? Symbol[] : _tags_by_group(group_of, states),
                      allosteric.oligomeric_state, sites))
    end
end

"""Whether the kinetic group labels of `steps` are numbered by their first step: the first
step of each group comes after the first step of each group with a smaller label."""
function _labels_in_order(steps::Vector{Step})
    seen = 0
    for e in steps
        e.kinetic_group > seen + 1 && return false
        seen = max(seen, e.kinetic_group)
    end
    true
end

"""A site's place in the order of a mechanism's sites: its ligands, its multiplicity, then
its ligands' states."""
_site_key(s::RegulatorySite) = (s.ligands, s.multiplicity, s.allosteric_states)

"""Group labels in step order renumbered by first occurrence, and the number each
written label took."""
function _number_groups(labels::Vector{Int}, G::Int)
    group_of, group = zeros(Int, G), similar(labels)
    next = 0
    for (k, l) in enumerate(labels)
        if group_of[l] == 0
            next += 1
            group_of[l] = next
        end
        group[k] = group_of[l]
    end
    group, group_of
end

"""The states by canonical group: `allosteric_states[l]`, the state of written label `l`,
moved to group `group_of[l]`."""
function _tags_by_group(group_of::Vector{Int}, allosteric_states::Vector{Symbol})
    tags = similar(allosteric_states)
    for (l, g) in enumerate(group_of)
        tags[g] = allosteric_states[l]
    end
    tags
end

"""Error unless `allosteric_states` holds one valid state (`_ALLOSTERIC_STATES`) for each
of `G` kinetic groups and the oligomeric state is at least 1."""
function _check_tags(allosteric_states::Vector{Symbol}, G::Int, oligomeric_state::Int)
    length(allosteric_states) == G ||
        error("Mechanism: allosteric_states length " *
              "$(length(allosteric_states)) must match the number of kinetic groups $G")
    oligomeric_state ≥ 1 ||
        error("Mechanism: oligomeric_state must be ≥ 1, " *
              "got $oligomeric_state")
    bad = findfirst(∉(_ALLOSTERIC_STATES), allosteric_states)
    bad === nothing ||
        error("Mechanism: catalytic group $bad has invalid allosteric state " *
              "$(allosteric_states[bad]) (must be one of $_ALLOSTERIC_STATES)")
    nothing
end

"""Error when a site ligand cannot name a metabolite (`_check_metabolite_names`) or names
an enzyme form of `catalytic`, which a text would read as the form, or when a ligand sits
on two sites: each ligand's constant is named after the ligand alone."""
function _check_ligands(sites::Vector{RegulatorySite}, catalytic::CatalyticScheme)
    ligs = Symbol[l for s in sites for l in s.ligands]
    _check_metabolite_names("Mechanism", ligs)
    for l in ligs
        l in catalytic.enzyme_forms &&
            error("Mechanism: ligand $l is also the name of an enzyme form; a text could " *
                  "not tell them apart")
    end
    dup = findfirst(i -> ligs[i] in view(ligs, 1:(i - 1)), eachindex(ligs))
    dup === nothing ||
        error(
            "Mechanism: allosteric regulator `$(ligs[dup])` is listed at two regulatory " *
            "sites; an allosteric regulator sits at exactly one regulatory site.")
    nothing
end

# ─── Identity ───────────────────────────────────────────────────────────

# `==` and `hash` by content for every value type whose fields hold vectors or structs,
# whose default `==` would compare identity: each field in declaration order, `hash`
# seeded with the type's name. A struct's default hash mixes in its type's identity, which
# changes between precompile sessions; a value's hash reads only symbols, integers and
# flags, so every process agrees on it.
for T in (RegulatorySite, Reactant, Regulator, EnzymeReaction, Step, CatalyticScheme,
          Allostery, Mechanism)
    fs = fieldnames(T)
    @eval Base.:(==)(a::$T, b::$T) =
        $(reduce((l, r) -> :($l && $r), [:(a.$f == b.$f) for f in fs]))
    @eval Base.hash(x::$T, h::UInt) =
        $(foldl((acc, f) -> :(hash(x.$f, $acc)), fs;
                init = :(hash($(QuoteNode(nameof(T))), h))))
end

# ─── Reads ──────────────────────────────────────────────────────────────

"""
    metabolites(m) -> Vector{Symbol}

The names of the metabolites of mechanism `m`, each once: those of its catalytic scheme,
its substrates, then its products, then its modifiers, each sorted by name, then the
ligands of an allosteric mechanism's regulatory sites that are not among them, site by
site. These are the concentrations its rate law reads, in the order `metabolites(rate)` of
its compiled law
([`compile_rate_equation`](@ref EnzymeRates.Derivation.compile_rate_equation)) gives
them as a tuple.

```jldoctest
julia> using EnzymeRates

julia> m = @enzyme_mechanism begin
           substrates: S
           products: P
           steps: begin
               E + S <--> ES
               ES <--> E + P
           end
       end;

julia> metabolites(m)
2-element Vector{Symbol}:
 :S
 :P
```
"""
metabolites(m::Mechanism) =
    unique!([m.catalytic.metabolites;
             [l for s in m.allosteric.regulatory_sites for l in s.ligands]])

# ─── The steps each conformation runs, and the kinetic groups ───────────

"""Each kinetic group's steps of `steps`, in step order, by the steps' `kinetic_group`
labels; each list stored at its length."""
function _groups(steps::Vector{Step})
    G = maximum(e -> e.kinetic_group, steps; init = 0)
    members = [Vector{Int}(undef, count(e -> e.kinetic_group == g, steps)) for g in 1:G]
    filled = zeros(Int, G)
    for (k, e) in enumerate(steps)
        g = e.kinetic_group
        members[g][filled[g] += 1] = k
    end
    members
end

"""
    steps_running_in_allosteric_conformations(m) -> (; active, inactive)

The steps of mechanism `m` that run in its active and in its inactive conformation, each a
`BitVector` with one flag per step of its catalytic scheme. A step runs in a conformation
unless its kinetic group's allosteric state excludes it there or the conformation's
exclusions cut it off from free E ([`Allostery`](@ref) has the rule). Every step of a
plain mechanism runs in both.
"""
steps_running_in_allosteric_conformations(m::Mechanism) =
    (; active = _steps_in_conformation(m, :OnlyI),
       inactive = _steps_in_conformation(m, :OnlyA))

"""The steps of `m` that run in the conformation that excludes the kinetic groups in
allosteric state `state` (`steps_running_in_allosteric_conformations`)."""
function _steps_in_conformation(m::Mechanism, state::Symbol)
    c, states = m.catalytic, m.allosteric.allosteric_states
    isempty(states) && return trues(length(c.steps))
    excluded = (e -> states[e.kinetic_group] === state).(c.steps)
    joined = _joined_to_free_E(c, excluded)
    (k -> !excluded[k] && joined[c.steps[k].from]).(eachindex(c.steps))
end

"""Which forms of `m` the steps not in `excluded` join to free E, form 1."""
function _joined_to_free_E(m::CatalyticScheme, excluded::BitVector)
    reached = falses(length(m.enzyme_forms))
    reached[1] = true
    grew = true
    while grew
        grew = false
        for (k, e) in enumerate(m.steps)
            (excluded[k] || reached[e.from] == reached[e.to]) && continue
            reached[e.from] = reached[e.to] = true
            grew = true
        end
    end
    reached
end

# ─── The mechanism as macro text ────────────────────────────────────────

"""
Print `m` as the `@enzyme_mechanism` text that builds it: its substrates, products and
modifiers (the metabolites that are neither), each list in stored order and an empty one
left out; an allosteric mechanism's oligomeric state; its steps (`_print_steps`), each
line ending in its kinetic group's allosteric state when it has one; then one
`regulatory_site(n):` line per regulatory site. The text is a function of `m` alone, so
equal mechanisms print one text. `DSL.parse_mechanism` reads the text of a mechanism the
macro accepts back as that mechanism: one that is well posed (`check_valid_rate_equation`).
"""
function Base.show(io::IO, m::Mechanism)
    c, a = m.catalytic, m.allosteric
    print(io, "@enzyme_mechanism begin")
    _print_line(io, "substrates", c.metabolites[c.substrates])
    _print_line(io, "products", c.metabolites[c.products])
    modifiers = [i for i in eachindex(c.metabolites)
                 if !(i in c.substrates || i in c.products)]
    _print_line(io, "modifiers", c.metabolites[modifiers])
    isempty(a.allosteric_states) ||
        _print_line(io, "oligomeric_state", [a.oligomeric_state])
    _print_steps(io, c, a.allosteric_states)
    for s in a.regulatory_sites
        _print_line(io, "regulatory_site($(s.multiplicity))",
                    ["$l::$t" for (l, t) in zip(s.ligands, s.allosteric_states)])
    end
    print(io, "\nend")
    nothing
end

"""Print the line `label: item, item, ...` of a mechanism's block, or nothing when
`items` is empty."""
function _print_line(io::IO, label::String, items::Vector)
    isempty(items) || print(io, "\n    ", label, ": ", join(items, ", "))
    nothing
end

"""
Print the `steps:` block of `m`: one line per kinetic group in group order, a group of
one step as the step (`step_text`), a larger group as the tuple of its steps in step
order, and, when `tags` holds one tag per group, each line ending in its group's tag.
"""
function _print_steps(io::IO, m::CatalyticScheme, tags::Vector{Symbol})
    print(io, "\n    steps: begin")
    for (g, members) in enumerate(_groups(m.steps))
        text = join((step_text(m, k) for k in members), ", ")
        print(io, "\n        ", length(members) == 1 ? text : "($text)")
        isempty(tags) || print(io, " :: ", tags[g])
    end
    print(io, "\n    end")
    nothing
end

"""Step `k` of `m` as a mechanism macro writes it, in its stored direction: each side is
its form's name and the metabolites the step takes up or gives off there, joined by `⇌`
(rapid equilibrium) or `<-->` (steady state): `E + S ⇌ ES`, `EA + B <--> EQ + P`."""
function step_text(m::CatalyticScheme, k::Int)
    e = m.steps[k]
    side(f::Int, held::Vector{Int}) = join([m.enzyme_forms[f]; m.metabolites[held]], " + ")
    side(e.from, e.taken_up) * (e.rapid_equilibrium ? " ⇌ " : " <--> ") *
        side(e.to, e.given_off)
end

# ─── Segments and cycles of a step mask ─────────────────────────────────

# Sub-graphs: a mask says which steps of a mechanism a sub-graph keeps, such as an
# allosteric mechanism's inactive conformation. Nothing here constructs a mechanism, so no
# check meant for the mechanisms a user writes fires on a sub-graph.

"""
    _present_species(m, on) -> BitVector

The forms of `m` that the steps `on` keeps touch, one flag per form: a sub-graph's law,
segments and network cover these.
"""
function _present_species(m::CatalyticScheme, on::BitVector)
    present = falses(length(m.enzyme_forms))
    for (k, e) in enumerate(m.steps)
        on[k] && (present[e.from] = present[e.to] = true)
    end
    present
end

"""
    rapid_equilibrium_segments(m, on, re = <m's rapid-equilibrium flags>) -> (; segment_of,
        segments, offsets, steps)

The rapid-equilibrium segments of the steps of `m` that `on` keeps, joined by the kept
steps that `re` flags. `segments[s]` lists segment `s`'s forms breadth-first from its first
form, each form's steps in the order of the groups, and the segments are numbered by their
first forms; `segment_of[f]` is form `f`'s segment, 0 for a form no kept step touches
(`_present_species`). `offsets[n, f]` is how many more of metabolite `n` (an index into
`m.metabolites`) form `f` carries than the form of its segment that carries fewest of it.
Within a segment, rapid equilibrium fixes each form's weight relative to any other as a
monomial in concentrations, one factor per metabolite taken up or given off along the path
between them; chemistry passes the counts through unchanged.

`steps` lists the kept steps that `re` leaves at steady state and that join two segments,
in step order: the edges of the steady-state graph between the segments.
"""
function rapid_equilibrium_segments(m::CatalyticScheme, on::BitVector,
                                    re::BitVector = BitVector(s.rapid_equilibrium
                                                              for s in m.steps))
    nf = length(m.enzyme_forms)
    # The kept rapid-equilibrium steps at each form, group by group: `(k, 1)` at the form
    # step `k` leaves in its stored direction, `(k, -1)` at the form it leaves in reverse.
    out = [Tuple{Int, Int}[] for _ in 1:nf]
    for k in reduce(vcat, _groups(m.steps); init = Int[])
        on[k] && re[k] || continue
        push!(out[m.steps[k].from], (k, 1))
        push!(out[m.steps[k].to], (k, -1))
    end
    present = _present_species(m, on)
    segment_of, segments = zeros(Int, nf), Vector{Int}[]
    offsets = zeros(Int, length(m.metabolites), nf)
    for root in 1:nf
        present[root] && segment_of[root] == 0 || continue
        members = [root]
        push!(segments, members)
        segment_of[root] = length(segments)
        head = 0
        while head < length(members)
            head += 1
            u = members[head]
            for (k, direction) in out[u]
                e = m.steps[k]
                v = direction == 1 ? e.to : e.from
                segment_of[v] == 0 || continue
                segment_of[v] = length(segments)
                push!(members, v)
                for n in axes(offsets, 1)
                    offsets[n, v] = offsets[n, u]
                end
                for x in e.taken_up
                    offsets[x, v] += direction
                end
                for x in e.given_off
                    offsets[x, v] -= direction
                end
            end
        end
    end
    for segment in segments, n in axes(offsets, 1)
        low = minimum(f -> offsets[n, f], segment)
        for f in segment
            offsets[n, f] -= low
        end
    end
    steps = [k for (k, e) in enumerate(m.steps)
             if on[k] && !re[k] && segment_of[e.from] != segment_of[e.to]]
    (; segment_of, segments, offsets, steps)
end

"""Whether the steps `k` with `on[k]` close a cycle of nonzero weight, step `k` joining
forms `from[k]` and `to[k]` with weight `weight[k]`: giving every form a potential that
rises by each step's weight along it then fails somewhere, which a weighted union-find in
the buffers `parent` and `offset` (one entry per form) detects."""
function has_nonzero_weight_cycle(from, to, weight, on, parent, offset)
    parent .= eachindex(parent); fill!(offset, 0)
    function root(x)
        p = 0
        while parent[x] != x
            p += offset[x]; x = parent[x]
        end
        x, p
    end
    for k in eachindex(from)
        on[k] || continue
        ra, pa = root(from[k]); rb, pb = root(to[k])
        if ra == rb
            pb == pa + weight[k] || return true
        else
            parent[rb] = ra; offset[rb] = pa + weight[k] - pb
        end
    end
    false
end

"""
    has_turnover_cycle(m::CatalyticScheme, on::BitVector) -> Bool

Whether a cycle of the steps of `m` that `on` keeps turns the reaction over: no potential
of the forms rises along each step by its uptake of the first substrate, which a weighted
union-find finds (`has_nonzero_weight_cycle`). A cycle of a
mechanism whose steps conserve atoms runs the reaction a whole number of times, so its
uptake of the first substrate is zero only when it turns nothing over; a hand-written
cycle that does not is read by that uptake alone. Without a turnover cycle the rate is zero
everywhere, and a turnover cycle of rapid-equilibrium steps alone makes it infinite;
`check_valid_rate_equation` rejects both.
"""
function has_turnover_cycle(m::CatalyticScheme, on::BitVector)
    # The first substrate, whose net uptake is the rate.
    s1 = m.substrates[1]
    uptake(list::Vector{Int}) = count(==(s1), list)
    n = length(m.enzyme_forms)
    has_nonzero_weight_cycle([e.from for e in m.steps], [e.to for e in m.steps],
                             [uptake(e.taken_up) - uptake(e.given_off) for e in m.steps],
                             on, zeros(Int, n), zeros(Int, n))
end
