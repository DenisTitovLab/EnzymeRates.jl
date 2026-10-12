# ABOUTME: Critic check of P5b's L-carrier rule: does moving the flip between two free forms joined
# ABOUTME: by a flux-carrying metabolite-free SS isomerization (Cleland iso) only rescale L?
using Random, Printf
const P5B = joinpath(@__DIR__, "..", "p5b")
include(joinpath(P5B, "Allo.jl"))
include(joinpath(P5B, "..", "common", "Oracle.jl"))

met(x, r) = Pop.PMet(x, r)
sp(name, conf, bound...) = Pop.PSpecies([met(b, r) for (b, r) in bound], conf, Symbol[], Symbol[], name)

# Copied from p5b/two_free.jl (one-step :NonequalAI groups, n = 1).
function mech(id, subs, prods, species, steps)
    nm(st, j, s) = s.re ? [Symbol("K_", st, "_", j)] : [Symbol("kf_", st, "_", j), Symbol("kr_", st, "_", j)]
    Pop.PMech(; id, source = id, substrates = subs, products = prods, metabolites = [subs; prods],
              species, steps, groups = [[j] for j in eachindex(steps)], allosteric = true,
              cat_multiplicity = 1, group_tags = fill(:NonequalAI, length(steps)),
              step_names_A = [nm("A", j, s) for (j, s) in enumerate(steps)],
              step_names_I = [nm("I", j, s) for (j, s) in enumerate(steps)])
end

include(joinpath(P5B, "validate_lib.jl"))

# Copied from p5b/two_free.jl: fit L' of the flip-at-f network at one point, report the
# largest mismatch at five other points (0 means the choice only reparametrizes L).
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
        lo, hi = big(10)^-30, big(10)^30
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

S, P = (:S, :S), (:P, :P)
for seed in (5, 6)
    rng = Xoshiro(seed)
    # Control: P5b's case 3 (F a dead end of an SS isomerization).
    dead = mech("ss_iso_dead_end", [:S], [:P],
                [sp(:E, :E), sp(:ES, :E, S), sp(:EP, :E, P), sp(:F, :F)],
                [Pop.PStep(1, 4, Pop.PMet[], Pop.PMet[], false),
                 Pop.PStep(1, 2, [met(S...)], Pop.PMet[], true),
                 Pop.PStep(2, 3, Pop.PMet[], Pop.PMet[], false),
                 Pop.PStep(1, 3, [met(P...)], Pop.PMet[], true)])
    compare(dead, 4, rng; label = "seed $seed control: dead-end SS iso E <--> F")
    # Cleland iso uni-uni: E + S ⇌ ES; ES <--> F + P; F <--> E (the isomerization carries all flux).
    iso = mech("cleland_iso", [:S], [:P],
               [sp(:E, :E), sp(:ES, :E, S), sp(:F, :F)],
               [Pop.PStep(1, 2, [met(S...)], Pop.PMet[], true),
                Pop.PStep(2, 3, Pop.PMet[], [met(P...)], false),
                Pop.PStep(3, 1, Pop.PMet[], Pop.PMet[], false)])
    compare(iso, 3, rng; label = "seed $seed Cleland iso: flux-carrying SS iso F <--> E")
end
