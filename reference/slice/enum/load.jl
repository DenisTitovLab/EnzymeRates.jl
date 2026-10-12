# ABOUTME: Loads Pop and a Slice module of only the files the enumeration stage needs (core,
# ABOUTME: thermo.jl and enum.jl), so that other stages' work in progress cannot break a run.
const S = "/tmp/claude-501/-home-denis-linux--julia-dev-EnzymeRates/12f42995-f5ba-4fda-9510-dd9b87f7d069/scratchpad"
include(joinpath(S, "proto/common/Population.jl"))
module Slice
const SRC = "/tmp/claude-501/-home-denis-linux--julia-dev-EnzymeRates/12f42995-f5ba-4fda-9510-dd9b87f7d069/scratchpad/slice/src"
for file in ("types.jl", "names.jl", "poly.jl", "convert.jl", "thermo.jl", "enum.jl")
    include(joinpath(SRC, file))
end
end
using Serialization, Printf, Statistics, SHA
"""SHA-256 (first 12 hex digits) of a source file, to record the version a run used."""
_testhelper_src_sha(f) = bytes2hex(open(sha256, joinpath(Slice.SRC, f)))[1:12]
