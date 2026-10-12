# ABOUTME: Two residual-free free forms in an allosteric mechanism: does the form that carries L
# ABOUTME: change the model? Synthetic mechanisms, each choice checked against the n = 1 oracle.
using Random, Printf
include(joinpath(@__DIR__, "Allo.jl"))
include(joinpath(@__DIR__, "..", "common", "Oracle.jl"))

met(x, r) = Pop.PMet(x, r)
sp(name, conf, bound...) = Pop.PSpecies([met(b, r) for (b, r) in bound], conf, Symbol[], Symbol[], name)

"""An allosteric `PMech` of `species` and one-step groups `steps` (all `:NonequalAI`, n = 1)
with names `<K|kf|kr>_<A|I>_<j>`; `first_free` picks which free form comes first."""
function mech(id, subs, prods, species, steps)
    nm(st, j, s) = s.re ? [Symbol("K_", st, "_", j)] : [Symbol("kf_", st, "_", j), Symbol("kr_", st, "_", j)]
    Pop.PMech(; id, source = id, substrates = subs, products = prods, metabolites = [subs; prods],
              species, steps, groups = [[j] for j in eachindex(steps)], allosteric = true,
              cat_multiplicity = 1, group_tags = fill(:NonequalAI, length(steps)),
              step_names_A = [nm("A", j, s) for (j, s) in enumerate(steps)],
              step_names_I = [nm("I", j, s) for (j, s) in enumerate(steps)])
end

include(joinpath(@__DIR__, "validate_lib.jl"))   # combined(md, vals): the n = 1 two-conformation network

"""The law of `m` (L on its first free form) against the mass-action network whose flip is at
the free form `f` instead: fits that network's L' to the law at one concentration point and
reports the largest mismatch over five others (0 means the choice only reparametrizes L)."""
function compare(m, f, rng; label)
    md = Allo.derive(m)
    setprecision(BigFloat, 512) do
        fit = Dict{Symbol, BigFloat}(n => Oracle.logu(rng) for n in md.names.fitted)
        fit[:Keq] = Oracle.logu(rng)
        vals = Allo.all_values(md, fit)
        cs = [Oracle.concs_at_ratio(m, rng, fit[:Keq] * 10 * Oracle.logu(rng)) for _ in 1:6]
        law(c) = Allo.rate(md, vals, c)
        own = maximum(abs(Oracle.steady_state_v(combined(md, vals)..., c) - law(c)) / abs(law(c)) for c in cs)
        alt(c, L) = Oracle.steady_state_v(combined(md, vals; flip = f, L)..., c)
        lo, hi = big(10)^-30, big(10)^30                       # bisection in log L'
        g(L) = alt(cs[1], L) - law(cs[1])
        sg = sign(g(lo))
        sign(g(hi)) == sg && return @printf("%s: no L' reproduces the law at the first point\n", label)
        for _ in 1:300
            mid = sqrt(lo * hi); sign(g(mid)) == sg ? (lo = mid) : (hi = mid)
        end
        Lp = sqrt(lo * hi)
        mism = maximum(abs(alt(c, Lp) - law(c)) / abs(law(c)) for c in cs[2:end])
        @printf("%s: law vs its own network %.1e; flip at %s instead: best L'/L = %.6g, mismatch at 5 other points %.2e\n",
                label, Float64(own), m.species[f].name, Float64(Lp / fit[:L]), Float64(mism))
    end
end

A, B, P, Q, S = (:A, :S), (:B, :S), (:P, :P), (:Q, :P), (:S, :S)
rng = Xoshiro(5)
# 1. Ping-pong written with a conformation label: E and F are both free, in different SS segments.
pp = mech("pingpong_EF", [:A, :B], [:P, :Q],
          [sp(:E, :E), sp(:EA, :E, A), sp(:F, :F), sp(:FB, :F, B)],
          [Pop.PStep(1, 2, [met(A...)], Pop.PMet[], true),          # E + A ⇌ EA
           Pop.PStep(2, 3, Pop.PMet[], [met(P...)], false),         # EA <--> F + P
           Pop.PStep(3, 4, [met(B...)], Pop.PMet[], true),          # F + B ⇌ FB
           Pop.PStep(4, 1, Pop.PMet[], [met(Q...)], false)])        # FB <--> E + Q
compare(pp, 3, rng; label = "ping-pong E/F in two SS segments")
# 2. Two free forms in one RE segment: E ⇌ F rapid, S binds E.
re = mech("re_EF", [:S], [:P],
          [sp(:E, :E), sp(:ES, :E, S), sp(:EP, :E, P), sp(:F, :F)],
          [Pop.PStep(1, 4, Pop.PMet[], Pop.PMet[], true),           # E ⇌ F
           Pop.PStep(1, 2, [met(S...)], Pop.PMet[], true),          # E + S ⇌ ES
           Pop.PStep(2, 3, Pop.PMet[], Pop.PMet[], false),          # ES <--> EP
           Pop.PStep(1, 3, [met(P...)], Pop.PMet[], true)])         # E + P ⇌ EP
compare(re, 4, rng; label = "E ⇌ F in one RE segment")
# 3. Two free forms in two SS segments joined by a metabolite-free SS isomerization E <--> F.
iso = mech("ss_iso_EF", [:S], [:P],
           [sp(:E, :E), sp(:ES, :E, S), sp(:EP, :E, P), sp(:F, :F)],
           [Pop.PStep(1, 4, Pop.PMet[], Pop.PMet[], false),         # E <--> F
            Pop.PStep(1, 2, [met(S...)], Pop.PMet[], true),         # E + S ⇌ ES
            Pop.PStep(2, 3, Pop.PMet[], Pop.PMet[], false),         # ES <--> EP
            Pop.PStep(1, 3, [met(P...)], Pop.PMet[], true)])        # E + P ⇌ EP
compare(iso, 4, rng; label = "E <--> F as an SS dead-end isomerization")
