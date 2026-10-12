# ABOUTME: Today's package verdict on each spec: build it with the EnzymeRates constructors (as
# ABOUTME: the DSL does) and derive its reduced rate equation; records each error message.
# Usage: julia --project=<EnzymeRates checkout> pkg_verdict.jl <mode>
using EnzymeRates, Serialization, Printf
const ER = EnzymeRates
const R = "/tmp/claude-501/-home-denis-linux--julia-dev-EnzymeRates/12f42995-f5ba-4fda-9510-dd9b87f7d069/scratchpad/round4/guard"
const MODE = isempty(ARGS) ? "smoke" : ARGS[1]

errmsg(e) = replace(first(sprint(showerror, e), 400), '\n' => ' ')

include(joinpath(R, "pkg_build.jl"))

function main()
    specs = deserialize(joinpath(R, "data", "specs_$MODE.jls"))
    out = Vector{Any}(undef, length(specs))
    t0 = time()
    for (i, r) in enumerate(specs)
        t1 = time_ns()
        m, ctor = try
            build(r.spec), nothing
        catch e
            nothing, errmsg(e)
        end
        t2 = time_ns()
        derive, vline = nothing, ""
        if m !== nothing
            derive = try
                s = ER.rate_equation_string(m)
                l = last(split(s, "\n"))
                vline = first(l, 300)
                nothing
            catch e
                errmsg(e)
            end
        end
        t3 = time_ns()
        out[i] = (key = r.key, ctor = ctor, derive = derive, vline = vline,
                  t_ctor = (t2 - t1) / 1e9, t_derive = (t3 - t2) / 1e9)
        i % 500 == 0 && (@printf("%d/%d %.0f s\n", i, length(specs), time() - t0); flush(stdout))
    end
    serialize(joinpath(R, "data", "pkg_$MODE.jls"), out)
    @printf("done %d in %.0f s\n", length(out), time() - t0)
end

main()
