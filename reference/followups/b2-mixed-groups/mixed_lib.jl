# ABOUTME: Step classification shared by the mixed-group scripts: plain binding, fused binding,
# ABOUTME: isomerization, Theorell-Chance; mixed groups and DSL-style printing of a mechanism.

metkey(x::Pop.PMet) = (x.name, x.role)

"""Step class: :plain (one metabolite consumed, none released, `to` = `from` plus it with
the same residual), :fused (a binding that is not plain), :iso, :tc (one taken up, one given
off), :other."""
function step_class(m::Pop.PMech, s::Pop.PStep)
    nc, nr = length(s.consumed), length(s.released)
    nc == 0 && nr == 0 && return :iso
    if nc == 1 && nr == 0
        a, b = m.species[s.from], m.species[s.to]
        plain = sort(metkey.([a.bound; s.consumed])) == sort(metkey.(b.bound)) &&
                sort(a.res_added) == sort(b.res_added) &&
                sort(a.res_subtracted) == sort(b.res_subtracted)
        return plain ? :plain : :fused
    end
    nc == 1 && nr == 1 && return :tc
    :other
end

"""Whether a form holds both a substrate and a product (a dead-end complex)."""
deadend(sp::Pop.PSpecies) =
    any(x -> x.role === :S, sp.bound) && any(x -> x.role === :P, sp.bound)

"""Groups holding at least one plain binding and at least one chemistry step."""
mixed_groups(m::Pop.PMech) =
    [g for (g, grp) in enumerate(m.groups)
     if any(j -> step_class(m, m.steps[j]) === :plain, grp) &&
        any(j -> step_class(m, m.steps[j]) !== :plain, grp)]

group_text(m, g) = join([Pop.step_text(m, m.steps[j]) for j in m.groups[g]], ", ")

"""Print every kinetic group of `m` as a parenthesized tuple, with its tag when allosteric."""
function print_mech(m::Pop.PMech)
    println("-- ", m.id, "  (", m.source, ")  family ", m.family,
            m.allosteric ? "  n = $(m.cat_multiplicity)" : "")
    mg = Set(mixed_groups(m))
    for (g, grp) in enumerate(m.groups)
        tag = m.allosteric ? "  [$(m.group_tags[g])]" : ""
        mark = g in mg ? "   <- mixed" : ""
        txt = length(grp) == 1 ? group_text(m, g) : "(" * group_text(m, g) * ")"
        println("     ", txt, tag, mark)
    end
end
