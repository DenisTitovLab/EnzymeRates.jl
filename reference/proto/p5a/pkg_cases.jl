# ABOUTME: Builds the P5a hand-written allosteric cases with EnzymeRates: the guard's verdict for
# ABOUTME: each, and Pop.PMech records (package outputs) for the admitted ones and tag twins of the rest.
using EnzymeRates, Random, Serialization, Printf, Statistics
include(joinpath(@__DIR__, "..", "common", "Population.jl"))
include(joinpath(@__DIR__, "..", "common", "export", "export_lib.jl"))

const SB = Module(:P5aCases)
Core.eval(SB, :(using EnzymeRates))

"""Construct an allosteric mechanism from DSL text; `nothing` and the message on failure."""
function build(txt)
    try
        Core.eval(SB, Meta.parse(txt)), ""
    catch e
        nothing, first(sprint(showerror, e isa LoadError ? e.error : e), 300)
    end
end

dsl(subs, prods, lines; n = 1) = """
@allosteric_mechanism begin
    substrates: $(join(subs, ", "))
    products: $(join(prods, ", "))
    catalytic_multiplicity: $n
    catalytic_steps: begin
        $(join(lines, "\n        "))
    end
end"""

cases = Tuple{String, String, String}[]     # (id, source/tags, dsl)
const T3 = (:OnlyA, :EqualAI, :NonequalAI)
for sre in (true, false), pre in (true, false), tS in T3, tC in T3, tP in T3
    kind = string(sre ? "RE" : "SS", "/", pre ? "RE" : "SS")
    lines = ["E + S $(sre ? "⇌" : "<-->") E(S) :: $tS",
             "E(S) <--> E(P) :: $tC",
             "E + P $(pre ? "⇌" : "<-->") E(P) :: $tP"]
    push!(cases, ("uni/$kind/$(tS),$(tC),$(tP)", "uni:$kind:$tS,$tC,$tP",
                  dsl(["S"], ["P"], lines)))
end
push!(cases, ("familyA_cut", "gate:metabD:OnlyA,EqualAI,OnlyA,EqualAI", dsl(["S", "B"], ["P"], [
    "E + S <--> E(S) :: OnlyA", "E(S) + B ⇌ E(S, B) :: EqualAI",
    "E(S, B) <--> E(P) :: OnlyA", "E + P ⇌ E(P) :: EqualAI"])))
push!(cases, ("familyA_cut_RE", "gate:metabD-RE-S", dsl(["S", "B"], ["P"], [
    "E + S ⇌ E(S) :: OnlyA", "E(S) + B ⇌ E(S, B) :: EqualAI",
    "E(S, B) <--> E(P) :: OnlyA", "E + P ⇌ E(P) :: EqualAI"])))
push!(cases, ("multiA", "gate:multiA", dsl(["A", "B"], ["P"], [
    "E + A ⇌ E(A) :: OnlyA", "E(A) + B ⇌ E(A, B) :: OnlyA",
    "E(A, B) <--> E(P) :: OnlyA", "E + P ⇌ E(P) :: EqualAI"])))
const PP = ["E + ATP ⇌ E(ATP) :: EqualAI",
            "E(ATP) <--> E(F16BP; residual = ATP - F16BP) :: OnlyA",
            "E(; residual = ATP - F16BP) + F16BP ⇌ E(F16BP; residual = ATP - F16BP) :: EqualAI",
            "E(; residual = ATP - F16BP) + F6P ⇌ E(F6P; residual = ATP - F16BP) :: OnlyA",
            "E(F6P; residual = ATP - F16BP) ⇌ E(ADP) :: EqualAI",
            "E + ADP ⇌ E(ADP) :: EqualAI"]
push!(cases, ("pp_err1", "gate:err1", dsl(["ATP", "F6P"], ["ADP", "F16BP"], PP)))
push!(cases, ("pp_dead", "gate:dead-inactive", dsl(["ATP", "F6P"], ["ADP", "F16BP"],
    [PP[1:4]; "E(F6P; residual = ATP - F16BP) ⇌ E(ADP) :: OnlyA"; PP[6]])))
push!(cases, ("pp_err1_ss", "err1 with SS F6P binding", dsl(["ATP", "F6P"], ["ADP", "F16BP"],
    [PP[1:3]; "E(; residual = ATP - F16BP) + F6P <--> E(F6P; residual = ATP - F16BP) :: OnlyA"; PP[5:6]])))
push!(cases, ("pp_re_chem", "gate:RE segment holding E and F", dsl(["A", "B"], ["P", "Q"], [
    "E + A ⇌ E(A) :: EqualAI", "E(A) <--> E(P; residual = A - P) :: OnlyA",
    "E(; residual = A - P) + P ⇌ E(P; residual = A - P) :: EqualAI",
    "E(; residual = A - P) + B ⇌ E(B; residual = A - P) :: EqualAI",
    "E(B; residual = A - P) ⇌ E(Q) :: OnlyA", "E + Q ⇌ E(Q) :: EqualAI"])))
push!(cases, ("pp_nonequal", "gate:NonequalAI ping-pong", dsl(["A", "B"], ["P", "Q"], [
    "E + A <--> E(A) :: EqualAI", "E(A) <--> E(; residual = A - P) + P :: NonequalAI",
    "E(; residual = A - P) + B <--> E(B; residual = A - P) :: EqualAI",
    "E(B; residual = A - P) <--> E + Q :: NonequalAI"])))

out = Pop.PMech[]
for (i, (id, src, txt)) in enumerate(cases)
    m, err = build(txt)
    if m !== nothing
        p = to_pmech(ER.AllostericMechanism(m), id, src; seed = 1000 + i, do_bench = false)
        push!(out, p)
        @printf("%-34s admitted  fitted=%d ok=%s %s\n", id, length(p.fitted_params), p.ok,
                isempty(p.stage_errors) ? "" : string(p.stage_errors))
    else
        # the structure of the same steps with every group :NonequalAI; the stdlib side
        # restores the requested tags from `source`
        tw, err2 = build(replace(txt, r":: (OnlyA|EqualAI|NonequalAI)" => ":: NonequalAI"))
        tw === nothing && (println("$id: twin failed too: $err2"); continue)
        p = to_pmech(ER.AllostericMechanism(tw), id, src; seed = 1000 + i, do_bench = false,
                     do_kcat = false)
        d = Dict(f => getfield(p, f) for f in fieldnames(Pop.PMech))
        d[:ok] = false; d[:error] = "rejected: " * err
        push!(out, Pop.PMech(; d...))
        @printf("%-34s REJECTED  %s\n", id, first(err, 110))
    end
end
mkpath(joinpath(@__DIR__, "data"))
serialize(joinpath(@__DIR__, "data", "cases.jls"), out)
println("wrote $(length(out)) cases")
