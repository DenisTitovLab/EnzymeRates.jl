# ABOUTME: Admission of allosteric tags (decision 6): one flip form, the exact cone test G1
# ABOUTME: on per-step classes, and the inactive-island rule.

"""
    admit(m) -> nothing or the rejection message

Whether the inactive conformation of `m` is the limit of its active one as the `:OnlyA`
steps lose their inactive affinity (a plain binding: `K_I → ∞`, or `kon_I → 0` with `koff_I`
fixed) or speed (any other step: both rates → 0, the ratio left to thermodynamics), so that
deletion plus pruning (`inactive_mask`) derives it; `nothing` for a plain or an admitted
mechanism. The tests run on the components of the forms that the steps outside `:OnlyA`
groups join, each `:OnlyA` plain binding an arc from its unbound to its bound side:
- one flip form: every free form lies in the first one's rapid-equilibrium segment of those
  steps, so the choice among them only rescales `L`;
- G1, the exact cone test: no arc stays inside a component and the arcs close no cycle.
  Thermodynamics makes each arc's equilibrium-ratio exponent a difference of component
  potentials and only arcs carry a sign, so a direction exists exactly then;
- islands: the flip form's component reaches every component along the arcs, so each
  island's weight vanishes; any other island's would depend on how its steps vanish.
"""
function admit(m::Mechanism)
    is_allosteric(m) || return nothing
    is_free(m.forms[1]) || return "an allosteric mechanism flips its conformation at a " *
        "free form (no bound metabolite, no residual); this one has none"
    kept = BitVector(m.tags[g] !== :OnlyA for g in m.group)
    segment = first(re_forest(m, kept))
    other = findfirst(i -> is_free(m.forms[i]) && segment[i] != segment[1],
                      eachindex(m.forms))
    other === nothing || return "the conformation flips at one free form, " *
        "$(form_text(m, 1)); the free form $(form_text(m, other)) is not joined to it by " *
        "rapid-equilibrium steps outside :OnlyA groups, so the choice would change the " *
        "model"
    comp = collect(eachindex(m.forms))
    top(i) = (while comp[i] != i; i = comp[i]; end; i)
    foreach(k -> kept[k] && (comp[top(m.edges[k].a)] = top(m.edges[k].b)), eachindex(kept))
    comp .= top.(eachindex(comp))
    fwd = orient(m)
    arcs = [(comp[o.a], comp[o.b], k) for (k, o) in
            ((k, oriented(m, k, fwd)) for k in findall(.!kept)) if is_plain_binding(m, o)]
    texts(js) = join((step_text(m, m.edges[arcs[j][3]]) for j in js), ", ")
    chord = findfirst(a -> a[1] == a[2], arcs)
    chord === nothing || return "in the inactive conformation the steps outside :OnlyA " *
        "groups still join the two sides of the :OnlyA binding $(texts(chord)), so " *
        "thermodynamics keeps its affinity finite"
    indeg = zeros(Int, length(comp))
    foreach(a -> indeg[a[2]] += 1, arcs)
    todo = [c for c in unique(comp) if indeg[c] == 0]
    while !isempty(todo)
        u = pop!(todo)
        for (x, y, _) in arcs
            x == u && (indeg[y] -= 1) == 0 && push!(todo, y)
        end
    end
    if any(>(0), indeg)
        back(c) = findfirst(a -> a[2] == c && indeg[a[1]] > 0, arcs)
        c = findfirst(>(0), indeg)
        for _ in comp                                   # afterwards c lies on a cycle
            c = arcs[back(c)][1]
        end
        loop = [back(c)]
        while arcs[loop[end]][1] != c
            push!(loop, back(arcs[loop[end]][1]))
        end
        return "the :OnlyA bindings $(texts(reverse!(loop))) close a cycle in the " *
            "inactive conformation, so thermodynamics forbids their affinities to vanish"
    end
    seen, todo = falses(length(comp)), [comp[1]]
    seen[comp[1]] = true
    while !isempty(todo)
        u = pop!(todo)
        for (x, y, _) in arcs
            x == u && !seen[y] && (seen[y] = true; push!(todo, y))
        end
    end
    lost = findfirst(c -> !seen[c], comp)
    lost === nothing || return "in the inactive conformation the forms " *
        join((form_text(m, i) for i in eachindex(comp) if comp[i] == comp[lost]), ", ") *
        " reach the free enzyme only through :OnlyA steps that are not bindings toward " *
        "them, so their weight depends on how fast those steps vanish"
    nothing
end
