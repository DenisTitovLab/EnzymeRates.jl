# ABOUTME: Prints a mechanism as the macro text that builds it: its names, then one line
# ABOUTME: per kinetic group, and for an allosteric one its tags and regulatory sites.

"""
Print `m` as the `@enzyme_mechanism` text that builds it: its substrates, products and
modifiers (the names its steps bind as modifiers), each list sorted and an empty one left
out, then its steps (`_print_steps`). The text is a function of `m` alone, so equal
mechanisms print one text. `DSL.parse_mechanism` reads the text of a mechanism the macros
accept back as that mechanism: one that is well posed (`check_valid_rate_equation`) and
whose names spell every string one way (`check_decodable_names`).
"""
function Base.show(io::IO, m::Mechanism)
    print(io, "@enzyme_mechanism begin")
    _print_names(io, m)
    _print_steps(io, m, Symbol[])
    print(io, "\nend")
    nothing
end

"""
Print `am` as the `@allosteric_mechanism` text that builds it: its catalytic graph's
names, its allosteric regulators (left out when it has no site), its catalytic
multiplicity, its steps with their tags (`_print_steps`), then one `regulatory_site(n):`
line per site.
"""
function Base.show(io::IO, am::AllostericMechanism)
    print(io, "@allosteric_mechanism begin")
    _print_names(io, am.catalytic)
    _print_line(io, "allosteric_regulators", sort(allosteric_regulators(am)))
    _print_line(io, "catalytic_multiplicity", [am.catalytic_multiplicity])
    _print_steps(io, am.catalytic, am.mwc_tags)
    for s in am.regulatory_sites
        _print_line(io, "regulatory_site($(s.multiplicity))",
                    ["$l::$t" for (l, t) in zip(s.ligands, s.mwc_tags)])
    end
    print(io, "\nend")
    nothing
end

# A compiled mechanism prints its mechanism's text. Each method compiles once, whatever
# the Sig: it reads the Sig at run time (`Mechanism(em)`), and `@nospecializeinfer` keeps
# a caller compiled for one Sig, such as `sprint`, from inferring a show for that Sig.
Base.@nospecializeinfer Base.show(io::IO, @nospecialize(em::EnzymeMechanism)) =
    show(io, Mechanism(em))
Base.@nospecializeinfer Base.show(io::IO, @nospecialize(aem::AllostericEnzymeMechanism)) =
    show(io, AllostericMechanism(aem))

"""Print the line `label: item, item, ...` of a mechanism's block, or nothing when
`items` is empty."""
function _print_line(io::IO, label::String, items::Vector)
    isempty(items) || print(io, "\n    ", label, ": ", join(items, ", "))
    nothing
end

"""Print the substrates:, products: and modifiers: lines of `m`."""
function _print_names(io::IO, m::Mechanism)
    _print_line(io, "substrates", m.substrates)
    _print_line(io, "products", m.products)
    _print_line(io, "modifiers", modifiers(m))
end

"""
Print the `steps:` block of `m`: one line per kinetic group in group order, a group of
one step as the step (`step_text`), a larger group as the tuple of its steps in edge
order, and, when `tags` holds one tag per group, each line ending in its group's tag.
"""
function _print_steps(io::IO, m::Mechanism, tags::Vector{Symbol})
    print(io, "\n    steps: begin")
    for (g, members) in enumerate(_group_members(m))
        text = join((step_text(m, k) for k in members), ", ")
        print(io, "\n        ", length(members) == 1 ? text : "($text)")
        isempty(tags) || print(io, " :: ", tags[g])
    end
    print(io, "\n    end")
    nothing
end
