# ABOUTME: Tests of the routing the integration added: one-argument printing, kcat(m, params)
# ABOUTME: for plain and allosteric mechanisms, derive's artifacts, and the expand driver.
include(joinpath(@__DIR__, "lib.jl"))
using Test

rng = Xoshiro(11)
specs = records_of("specs")
@testset "printing and kcat route plain and allosteric mechanisms" begin
    for r in specs
        m = r.c.mech
        S.admit(m) === nothing || continue
        d = S.derive(m)
        @test S.rate_equation_string(m) == d.printed
        @test S.coefficient_form(m) == d.coefficients
        @test d.law isa (S.is_allosteric(m) ? S.AlloLaw : S.Law)
        d.kcat isa String && continue
        params = (; (n => Oracle.logu(rng) for n in S.parameters(m))...)
        vars, terms = d.kcat
        x = [vars.conc[i] ? big(1) : BigFloat(params[n]) for (i, n) in enumerate(vars.names)]
        @test S.kcat(m, params) ≈ maximum(S.evaluate(a, x) / S.evaluate(b, x) for (a, b) in terms)
    end
end

@testset "derive rejects with admission's message" begin
    r = first(filter(r -> S.admit(r.c.mech) !== nothing, records_of("p5a_adv_cases")))
    msg = S.admit(r.c.mech)
    @test_throws ErrorException(msg) S.derive(r.c.mech)
end

@testset "expand: the three moves' children, checked and admitted" begin
    atoms = Dict(:A => [:C => 1, :X => 1], :B => [:N => 1], :P => [:C => 1],
                 :Q => [:N => 1, :X => 1])
    for r in records_of("init_bibi_pp")[1:40]
        m = r.c.mech
        cs = S.expand(m, [:I, :A, :Q], atoms)
        moves = unique!([S.flip_children(m); S.split_children(m);
                         S.dead_end_children(m, [:I, :A, :Q])])
        @test Set(cs) == Set(moves) && allunique(cs)
    end
end
