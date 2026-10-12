# ABOUTME: Prints the package's allosteric rate-equation strings for a few spec mechanisms,
# ABOUTME: to see today's MWC printed form (allosteric stage, exploration only).
const SCRATCH = "/tmp/claude-501/-home-denis-linux--julia-dev-EnzymeRates/12f42995-f5ba-4fda-9510-dd9b87f7d069/scratchpad"
include(joinpath(SCRATCH, "proto", "common", "Population.jl"))
specs = Pop.load_population("specs")
for p in specs
    p.allosteric || continue
    println("== ", p.id, "  ", p.source, "  n=", p.cat_multiplicity, " tags=", p.group_tags,
            " sites=", [(s.ligands, s.states, s.multiplicity) for s in p.sites])
    for s in p.steps
        println("   ", s.from, " -> ", s.to, " cons=", [x.name for x in s.consumed],
                " rel=", [x.name for x in s.released], " re=", s.re)
    end
    println(p.rate_equation_string)
end
