# ABOUTME: Builds law-type variant F (L1's layout as the verifiers asked): C's Int-only engine, E's mutable
# ABOUTME: Poly, mutable Law and AlloLaw with const fields, and the flip-only inactive law stored as (0, 1).
import os, sys
SP = "/tmp/claude-501/-home-denis-linux--julia-dev-EnzymeRates/12f42995-f5ba-4fda-9510-dd9b87f7d069/scratchpad"
SRC = f"{SP}/slice/src"
C = f"{SP}/phase45/research5/law/variants/C"
E = f"{SP}/phase45/research5/verify/variantE"
OUT = f"{SP}/phase45/research5/synth/variantF"

def read(path):
    return open(path, encoding="utf-8").read()

def patch(text, pairs, label):
    for old, new in pairs:
        n = text.count(old)
        if n != 1:
            sys.exit(f"{label}: patch matches {n} times: {old[:70]!r}")
        text = text.replace(old, new)
    return text

LAW_OLD = """struct Law{T <: Integer}
    vars::Vars
    H::Poly{T}
    Q::Poly{T}
    D::Poly{T}
    sigma::Vector{Poly{T}}
    trees::Vector{Poly{T}}
    shift::Mono
end"""
LAW_NEW = """mutable struct Law{T <: Integer}
    const vars::Vars
    const H::Poly{T}
    const Q::Poly{T}
    const D::Poly{T}
    const sigma::Vector{Poly{T}}
    const trees::Vector{Poly{T}}
    const shift::Mono
end"""
ALLO_OLD = """struct AlloLaw{T <: Integer}
    vars::Vars
    H::Poly{T}
    active::Law{T}
    inactive::Union{Nothing, Law{T}}
    fA::Poly{T}
    fI::Poly{T}
    sites::Vector{Tuple{Poly{T}, Poly{T}, Int}}
    n::Int
end"""
ALLO_NEW = """mutable struct AlloLaw{T <: Integer}
    const vars::Vars
    const H::Poly{T}
    const active::Law{T}
    const inactive::Law{T}
    const fA::Poly{T}
    const fI::Poly{T}
    const sites::Vector{Tuple{Poly{T}, Poly{T}, Int}}
    const n::Int
end

\"\"\"The inactive law of a conformation that is the flip form alone: `Q = 0`, `D = 1`, no
segment, the flip form's weight 1.\"\"\"
flip_only_law(LA::Law{T}) where {T} =
    Law(LA.vars, LA.H, zero(Poly{T}), one(Poly{T}), Poly{T}[], [one(Poly{T})], Mono())
\"\"\"Whether `L` is `flip_only_law`'s law.\"\"\"
flip_only(L::Law) = isempty(L.sigma)"""

files = {
    "poly.jl": read(f"{E}/poly.jl"),
    "engine.jl": patch(read(f"{C}/engine.jl"), [(LAW_OLD, LAW_NEW)], "engine.jl"),
    "allosteric.jl": patch(read(f"{C}/allosteric.jl"), [
        (ALLO_OLD, ALLO_NEW),
        ("    LI = any(imask) ? law(m, sol, imask; inactive = true) : nothing\n",
         "    LI = any(imask) ? law(m, sol, imask; inactive = true) : flip_only_law(LA)\n"),
        ("cross_weights(flip(LA), LI === nothing ? one(P) : flip(LI))",
         "cross_weights(flip(LA), flip(LI))"),
        ("    A.inactive === nothing ? (zero(Poly{T}), one(Poly{T})) : (A.inactive.Q, A.inactive.D)",
         "    (A.inactive.Q, A.inactive.D)"),
    ], "allosteric.jl"),
    "codegen.jl": patch(read(f"{SRC}/codegen.jl"), [
        ("    QI, DI = A.inactive === nothing ? (zero(Poly{T}), one(Poly{T})) :",
         "    QI, DI = flip_only(A.inactive) ? (zero(Poly{T}), one(Poly{T})) :"),
    ], "codegen.jl"),
    "print_allo.jl": patch(read(f"{SRC}/print_allo.jl"), [
        ("    A.inactive === nothing || push!(cs, conformation(m, A.inactive, \"I\"))",
         "    flip_only(A.inactive) || push!(cs, conformation(m, A.inactive, \"I\"))"),
    ], "print_allo.jl"),
}
for name, text in files.items():
    lines = text.split("\n")
    assert lines[0].startswith("# ABOUTME:") and lines[1].startswith("# ABOUTME:"), name
    lines[0] = "# ABOUTME: Variant F (synth/mkF.py) of the slice file below. " + lines[0][len("# ABOUTME: "):]
    open(f"{OUT}/{name}", "w", encoding="utf-8").write("\n".join(lines))

ORDER = ["types.jl", "names.jl", "poly.jl", "convert.jl", "thermo.jl", "api.jl", "engine.jl",
         "normalize.jl", "print.jl", "allosteric.jl", "admission.jl", "print_allo.jl",
         "codegen.jl", "enum.jl", "key.jl"]
open(f"{OUT}/SliceF.jl", "w", encoding="utf-8").write(
    "# ABOUTME: Law-type variant F of the slice as module SliceF: this directory's patched files\n"
    "# ABOUTME: (made by ../mkF.py) in place of the slice's originals.\n"
    "module SliceF\n"
    f"for file in ({', '.join(repr(f).replace(chr(39), chr(34)) for f in ORDER)})\n"
    "    local p = joinpath(@__DIR__, file)\n"
    f"    include(isfile(p) ? p : joinpath(\"{SRC}\", file))\n"
    "end\n"
    "end # module\n")
print("variant F written:", sorted(files))
