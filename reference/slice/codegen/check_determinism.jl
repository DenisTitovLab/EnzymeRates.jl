# ABOUTME: Codegen determinism: the SHA-256 of every record's generated body, saved per process;
# ABOUTME: a second process compares against the first run's file.
include(joinpath(@__DIR__, "lib.jl"))
using SHA
const OUT = joinpath(@__DIR__, "data", "body_hashes.jls")
hashes = Dict{Tuple{String, Int}, String}()
t0 = time()
for r in load_records()
    b = try
        S.body(S.program(r.c.mech))
    catch err
        println("$(r.set) $(r.p.id): ", sprint(showerror, err))
        continue
    end
    hashes[(r.set, r.c.index)] = bytes2hex(sha256(string(b)))
end
@printf("%d bodies in %.0f s\n", length(hashes), time() - t0)
if isfile(OUT)
    old = deserialize(OUT)
    diff = count(k -> get(old, k, "") != hashes[k], keys(hashes))
    println("against the first process: $(length(old)) bodies, $diff differ")
else
    serialize(OUT, hashes)
    println("saved; run again in a new process to compare")
end
