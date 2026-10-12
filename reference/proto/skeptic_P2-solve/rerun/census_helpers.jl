# ABOUTME: Shared helpers of the P2 census scripts: package slots, my slots, class checks.
# ABOUTME: Included after Thermo.jl; stdlib only.
"""(group index, stored slot) of a package parameter name: the first step that renders it."""
function pkg_slot(m, name)
    for (gi, g) in enumerate(m.groups), j in g, (r, n) in enumerate(m.step_names[j])
        n === name && return (gi, r)
    end
    nothing
end

"""(group, stored slot) of one of my columns: slot in the stored orientation."""
function my_slot(sol, i)
    c = sol.cols[i]; s = sol.cs[c.j]
    g = findfirst(gr -> c.j in gr, sol.m.groups)
    (g, s.re ? 1 : (s.flipped ? 3 - c.r : c.r))
end

"""Whether class (group g, stored slot r) is a chemistry member's constant against its physical
forward direction, or an RE chemistry constant."""
function textbook_class(sol, cl)
    g, r = cl
    any(sol.m.groups[g]) do j
        s = sol.cs[j]
        s.kind >= 2 || return false
        s.re && return true
        i = findfirst(c -> c.j == j && c.r == (s.flipped ? 3 - r : r), sol.cols)
        sol.cols[i].prio[4] == 1
    end
end

"""Short description of a class: the most chemical member's kind, RE/SS, against/forward."""
function class_desc(sol, cl)
    g, r = cl
    j = sol.m.groups[g][argmax([sol.cs[j].kind for j in sol.m.groups[g]])]
    s = sol.cs[j]
    kind = ("bind", "chem", "iso")[s.kind]
    s.re && return "RE $(kind) K"
    i = findfirst(c -> c.j == j && c.r == (s.flipped ? 3 - r : r), sol.cols)
    "SS $(kind) $(sol.cols[i].prio[4] == 1 ? "reverse" : "forward")"
end

