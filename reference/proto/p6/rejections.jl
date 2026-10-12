# ABOUTME: Prints why the P6 constructor rejected the :inh and :roles mutations of equality.jl.
# ABOUTME: Re-runs those two mutations only and tallies the error messages.
include(joinpath(@__DIR__, "harness.jl"))
src = read(joinpath(@__DIR__, "equality.jl"), String)
include_string(Main, src[findfirst("function mutate", src)[1]:findfirst("kinds = (", src)[1]-1])
ms = _testhelper_load_all()
msgs = Dict{String, Int}()
for (i, m) in enumerate(ms)
    raw = _testhelper_toraw(m; rng = Xoshiro(77i), flip = bitrand(Xoshiro(78i), length(m.steps)))
    for k in (:inh, :roles)
        r = mutate(raw, k, m)
        r === nothing && continue
        try _testhelper_gm(r) catch e
            key = string(k, ": ", sprint(showerror, e), " (", m.id, ")")
            key = k === :roles ? string(k, ": ", sprint(showerror, e)) : key
            msgs[key] = get(msgs, key, 0) + 1
        end
    end
end
for (k, v) in msgs; println(v, " × ", k); end
