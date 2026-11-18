# ABOUTME: Rewrites the Mermaid module map of docs/src/developer.md and the module table of
# ABOUTME: .claude/CLAUDE.md from src/; run it after a change test_boundaries.jl reports.
include("module_map.jl")

for (block, path) in (_testhelper_MAP_BLOCK => ("docs", "src", "developer.md"),
                      _testhelper_TABLE_BLOCK => (".claude", "CLAUDE.md"))
    file = normpath(joinpath(@__DIR__, "..", path...))
    write(file, _testhelper_replace_generated_block(
        read(file, String), block, block.generate(_testhelper_SRC)))
end
