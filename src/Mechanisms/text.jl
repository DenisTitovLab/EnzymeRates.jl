# ABOUTME: Prints a mechanism as the macro text that builds it: its names, then one line
# ABOUTME: per kinetic group, and for an allosteric one its states and regulatory sites.

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
