# ABOUTME: Rewrites from src/ the module map and file table of docs/src/developer.md and
# ABOUTME: the module table of .claude/CLAUDE.md; run it when test_boundaries.jl asks.
include("module_map.jl")

for (block, path) in (_testhelper_MAP_BLOCK => ("docs", "src", "developer.md"),
                      _testhelper_FILE_TABLE_BLOCK => ("docs", "src", "developer.md"),
                      _testhelper_TABLE_BLOCK => (".claude", "CLAUDE.md"))
    file = normpath(joinpath(@__DIR__, "..", path...))
    write(file, _testhelper_replace_generated_block(
        read(file, String), block, block.generate(_testhelper_SRC)))
end
