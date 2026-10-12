# ABOUTME: Confirms with the package that its own moves emit the flagged children: runs the
# ABOUTME: package's RE→SS flip / group split on each parent and looks for the child.
using EnzymeRates, Serialization, Printf
const ER = EnzymeRates
const R = "/tmp/claude-501/-home-denis-linux--julia-dev-EnzymeRates/12f42995-f5ba-4fda-9510-dd9b87f7d069/scratchpad/round4/guard"
include(joinpath(R, "pkg_build.jl"))
fl = deserialize(joinpath(R, "data", "enum_flagged.jls"))
found = 0
for r in fl
    p, c = build(r.parent), build(r.child)
    kids = r.move === :flip ? ER._expand_re_to_ss(p) : ER._expand_split_kinetic_group(p)
    hit = any(==(c), kids)
    hit && (global found += 1)
    s = ER.rate_equation_string(c)
    @printf("%-7s %-6s emitted by the package: %-5s derives: yes  %s\n", r.level, r.move, hit, first(r.face, 70))
end
println("flagged children emitted by the package's own moves: $found of $(length(fl))")
