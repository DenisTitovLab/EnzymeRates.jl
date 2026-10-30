# ABOUTME: Rewrites the Mermaid module map in docs/src/developer.md from src/; run it with
# ABOUTME: `julia test/update_module_map.jl` after a change that test_boundaries.jl reports.
include("module_map.jl")

const _testhelper_DEVELOPER_PAGE = normpath(joinpath(@__DIR__, "..", "docs", "src",
    "developer.md"))

write(_testhelper_DEVELOPER_PAGE, _testhelper_replace_module_map(
    read(_testhelper_DEVELOPER_PAGE, String), _testhelper_module_map(_testhelper_SRC)))
