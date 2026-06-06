# ABOUTME: Extracts ```julia code blocks from README.md and runs them in
# ABOUTME: one REPL session, skipping blocks tagged # README-SKIP-IN-TEST.

using Test
using EnzymeRates
using Random

@testset "README runs" begin
    md = replace(read(joinpath(@__DIR__, "..", "README.md"), String), "\r\n" => "\n")
    blocks = String[]
    for m in eachmatch(r"```julia\n(.*?)\n```"s, md)
        block = m.captures[1]
        startswith(strip(block), "# README-SKIP-IN-TEST") && continue
        push!(blocks, block)
    end
    @test !isempty(blocks)

    script = join(blocks, "\n\n")
    sandbox = Module()
    Core.eval(sandbox, :(using EnzymeRates, Random))
    Core.eval(sandbox, Meta.parse("begin\n$script\nend"))
end
