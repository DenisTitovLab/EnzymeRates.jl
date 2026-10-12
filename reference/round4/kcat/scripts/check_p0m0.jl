# ABOUTME: Counts laws whose rate at zero products and zero modifiers (the substrate-only sup search
# ABOUTME: points) is not finite and positive anywhere, per label.
using Serialization, Printf
const KDATA = "/tmp/claude-501/-home-denis-linux--julia-dev-EnzymeRates/12f42995-f5ba-4fda-9510-dd9b87f7d069/scratchpad/round4/kcat/data"
for l in ("core", "expand_terter", "terter_sample", "terter_nounique")
    rows = filter(x -> x.ok, deserialize(joinpath(KDATA, "rows_$l.jls")))
    bad = filter(x -> any(p -> !(isfinite(p.sups) && p.sups > 0), x.pts), rows)
    withmods = count(x -> !isempty(x.mods), rows)
    @printf("%-16s laws %5d (with modifiers %4d): substrate-only sup (products 0, modifiers 0) not finite > 0: %d %s\n",
            l, length(rows), withmods, length(bad), join([x.id for x in bad], ", "))
end
