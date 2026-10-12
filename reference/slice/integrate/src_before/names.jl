# ABOUTME: The naming chokepoint: every form, side and constant name of a mechanism,
# ABOUTME: rendered from its full graph, and the DSL text of forms and steps for messages.

"""Table entry `i` as names write it: its name, with `mod` on a reactant's modifier copy."""
label(m::Mechanism, i) = string(m.mets[i].name, is_dual(m, m.mets[i]) ? "mod" : "")

"""The labels of the entries counted in `p`, each as often as counted."""
labels(m::Mechanism, p::Pack) =
    [label(m, i) for i in eachindex(m.mets) for _ in 1:count_at(p, i)]

"""Form `i`'s name: conformation, bound labels, then any residual as `_res_plusA_minusP`
(`E`, `EAB`, `EAmod`, `EB_res_plusA_minusP`)."""
function form_name(m::Mechanism, i)
    f = m.forms[i]
    s = string(m.confs[f.conf], labels(m, f.bound)...)
    f.res == 0 && return s
    res = [x.name for x in mets_of(m, f.res)]
    string(s, "_res", ("_plus$x" for x in res if x in m.substrates)...,
           ("_minus$x" for x in res if x in m.products)...)
end

"""A step side: form `i`'s name and the labels of the free metabolites `p`, joined by
`_`."""
side_name(m::Mechanism, i, p::Pack) = join([form_name(m, i); labels(m, p)], "_")

"""Today's free-enzyme forms: no rapid-equilibrium step takes a metabolite up into them,
and they hold no metabolite or some step takes one up into them."""
function free_forms(m::Mechanism, fwd)
    into, re_into = falses(length(m.forms)), falses(length(m.forms))
    for k in eachindex(m.edges)
        e = oriented(m, k, fwd)
        e.cons == 0 && continue
        into[e.b] = true
        re_into[e.b] |= e.re
    end
    [!re_into[i] && (m.forms[i].bound == 0 || into[i]) for i in eachindex(m.forms)]
end

"""
    naming_reps(m, fwd = orient(m)) -> rep[g]

Each kinetic group's naming representative, today's rule: the member of lowest priority
(a free-enzyme metabolite step at rapid equilibrium -1, at steady state 0, a metabolite
step between bound forms 10, an isomerization 20), ties broken by the rendered names of
its forms in physical-forward order.
"""
function naming_reps(m::Mechanism, fwd = orient(m))
    free = free_forms(m, fwd)
    function key(k)
        e = oriented(m, k, fwd)
        iso, atfree = e.cons == 0 && e.rel == 0, free[e.a] || free[e.b]
        (e.re && !iso && atfree ? -1 : iso ? 20 : atfree ? 0 : 10,
         form_name(m, e.a), form_name(m, e.b))
    end
    [argmin(key, members(m, g)) for g in 1:ngroups(m)]
end

"""The state tag of group `g`'s constants: none in a plain mechanism or an `:EqualAI` group,
else `A_` or, in the inactive state, `I_`."""
state_tag(m::Mechanism, g, inactive::Bool) =
    !is_allosteric(m) || m.tags[g] === :EqualAI ? "" : inactive ? "I_" : "A_"

"""
    constant_names(m; inactive = false) -> names[k]

The names of edge `k`'s constants in its physical-forward direction: `[K]` at rapid
equilibrium, `[kf, kr]` at steady state. An edge renders its kinetic group's naming
representative (`naming_reps`), except that in the inactive state of an allosteric
mechanism an `:OnlyA` edge, untied there, renders itself. A rate constant is named after
the direction it drives; a binding's `K` after its release, as a dissociation constant;
any other `K` after the physical-forward direction.
"""
function constant_names(m::Mechanism; inactive::Bool = false)
    fwd = orient(m)
    reps = naming_reps(m, fwd)
    map(eachindex(m.edges)) do k
        g = m.group[k]
        own = inactive && is_allosteric(m) && m.tags[g] === :OnlyA
        e = oriented(m, own ? k : reps[g], fwd)
        a, b = side_name(m, e.a, e.cons), side_name(m, e.b, e.rel)
        t = state_tag(m, g, inactive)
        e.re ? [Symbol("K_", t, is_binding(e) ? "$(b)_to_$a" : "$(a)_to_$b")] :
               [Symbol("k_", t, a, "_to_", b), Symbol("k_", t, b, "_to_", a)]
    end
end

"""Each site ligand with the names of its constants: `K_A_Rreg` (`:OnlyA`), `K_I_Rreg`
(`:OnlyI`), both (`:NonequalAI`), or the bare `K_Rreg` (`:EqualAI`)."""
regulator_names(m::Mechanism) =
    [l => (st === :EqualAI ? [Symbol("K_", l, "reg")] :
           [Symbol("K_", t, "_", l, "reg") for t in ("A", "I")
            if st === :NonequalAI || st === Symbol("Only", t)])
     for s in m.sites for (l, st) in zip(s.ligands, s.states)]

"""Error when two forms render one name, or two constants (each group's, each inactive
class's, each regulator's, `L`, `Keq`, `E_total`) render one name."""
function check_names(m::Mechanism)
    fs = [form_name(m, i) for i in eachindex(m.forms)]
    for i in eachindex(fs), j in 1:i - 1
        fs[i] == fs[j] && error("the enzyme forms $(form_text(m, j)) and " *
                                "$(form_text(m, i)) both render the name $(fs[i]); " *
                                "rename a metabolite to tell them apart")
    end
    reps, active = naming_reps(m), constant_names(m)
    names = [n for k in reps for n in active[k]]
    if is_allosteric(m)
        inactive = constant_names(m; inactive = true)
        append!(names, (n for k in eachindex(m.edges) for n in inactive[k]
                        if m.tags[m.group[k]] === :OnlyA ||
                           (m.tags[m.group[k]] === :NonequalAI && k in reps)))
        append!(names, (n for (_, ns) in regulator_names(m) for n in ns))
        push!(names, :L)
    end
    push!(names, :Keq, :E_total)
    i = findfirst(j -> names[j] in view(names, 1:j - 1), eachindex(names))
    i === nothing || error("two constants render the name $(names[i])")
end

"""Metabolite `x` as the DSL writes it: a reactant's modifier copy as `A::Mod`."""
met_text(m::Mechanism, x::Met) = string(x.name, is_dual(m, x) ? "::Mod" : "")

"""Form `i` as the DSL writes it: `E`, `E(A, P::Mod)`, `E(B; residual = A - P)`."""
function form_text(m::Mechanism, i)
    f = m.forms[i]
    ms = [met_text(m, x) for x in mets_of(m, f.bound)]
    res = [x.name for x in mets_of(m, f.res)]
    added, lost = join(filter(in(m.substrates), res), " + "),
                  join(" - $x" for x in res if x in m.products)
    r = isempty(res) ? "" : "; residual = " * lstrip(added * lost)
    conf = string(m.confs[f.conf])
    isempty(ms) && isempty(r) ? conf : "$conf($(join(ms, ", "))$r)"
end

"""Edge `e` as the DSL writes it, `a` side first."""
function step_text(m::Mechanism, e::Edge)
    side(i, p) = join([form_text(m, i); [met_text(m, x) for x in mets_of(m, p)]], " + ")
    string(side(e.a, e.cons), e.re ? " ⇌ " : " <--> ", side(e.b, e.rel))
end
