# ABOUTME: The slice's expansion driver (expand) on the 60 sampled ter-ter parents against the
# ABOUTME: package's plain children of the same parents, and admission of its allosteric ones.
include(joinpath(@__DIR__, "lib.jl"))

explog = read(joinpath(@__DIR__, "logs", "export_terter_expand.log"), String)
picks = parse.(Int, split(match(r"parents \[([0-9, ]+)\]", explog)[1], ", "))
failed = [parse(Int, x[1]) for x in eachmatch(r"init_terter/(\d+):\s+0 children .*failed",
                                              explog)]
init = Dict{Int, Any}()
for k in 1:11, (i, p) in enumerate(deserialize(joinpath(IDATA, "pop_init_terter.part$k.jls")))
    j = (k - 1) * 25_000 + i
    j in picks && (init[j] = S.Converted(p, j).mech)
end
atoms = Dict(:A => [:C => 1], :B => [:N => 1], :D => [:X => 1], :P => [:C => 1],
             :Q => [:N => 1], :R => [:X => 1])
mine, per = Set{S.Mechanism}(), Dict{Int, Int}()
t = @elapsed for j in picks
    cs = S.expand(init[j], Symbol[], atoms)
    per[j] = length(cs)
    j in failed || union!(mine, cs)
end
pkg = [S.Converted(p, i).mech for (i, p) in
       enumerate(deserialize(joinpath(IDATA, "pop_expand_terter.jls")))]
plain = Set(filter(!S.is_allosteric, pkg))
allo = filter(S.is_allosteric, pkg)
@printf("%d parents (%d where the package's expansion threw: %s, slice children %s); slice \
        children over the others %d, package plain children %d; equal sets: %s; only slice \
        %d, only package %d; slice expansion %.2f s in all\n", length(picks), length(failed),
        failed, [per[j] for j in failed], length(mine), length(plain), mine == plain,
        length(setdiff(mine, plain)), length(setdiff(plain, mine)), t)
rej = [S.admit(m) for m in allo]
@printf("package's allosteric children: %d; admission rejects %d (%s)\n", length(allo),
        count(!isnothing, rej), join(unique(first(x, 60) for x in rej if x !== nothing), " | "))
