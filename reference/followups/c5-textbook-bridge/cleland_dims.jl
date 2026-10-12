# ABOUTME: C5 follow-up: dimension (Keq fixed) of each Cleland textbook formula with its Haldane relations
# ABOUTME: imposed, as the rank of ∂v/∂log(free constants) over random concentrations; stdlib only.
include("/tmp/claude-501/-home-denis-linux--julia-dev-EnzymeRates/12f42995-f5ba-4fda-9510-dd9b87f7d069/scratchpad/proto/skeptic_P1-engine/SkLib.jl")
include("/tmp/claude-501/-home-denis-linux--julia-dev-EnzymeRates/12f42995-f5ba-4fda-9510-dd9b87f7d069/scratchpad/proto/p4/bridge_text.jl")
using Random, Printf
logu(rng) = exp(big(log(10.0)) * (4 * big(rand(rng)) - 2))

function dim(rate, build, nfree, mets, rng)
    setprecision(BigFloat, 512) do
        Keq = logu(rng); x0 = [logu(rng) for _ in 1:nfree]
        pts = [(; (m => logu(rng) for m in mets)...) for _ in 1:(3 * nfree + 6)]
        h = big(10)^-30
        J = [(rate(build(x0 .* [i == j ? exp(h) : big(1) for i in 1:nfree], Keq), c) -
              rate(build(x0 .* [i == j ? exp(-h) : big(1) for i in 1:nfree], Keq), c)) / (2h)
             for c in pts, j in 1:nfree]
        Sk.numeric_rank(J; tol = big(10)^-20)
    end
end
rng = Random.Xoshiro(17)
fams = [
    ("uni-uni: V1, Ka, Kp free; V2 from V1*Kp/(V2*Ka) = Keq", uniuni, 3, (:A, :P),
     (x, Keq) -> (V1 = x[1], Ka = x[2], Kp = x[3], V2 = x[1] * x[3] / (Keq * x[2]))),
    ("ordered: 8 free; V2 from H1, Kip from H2", ordered, 8, (:A, :B, :P, :Q),
     (x, Keq) -> (V2 = x[1] * x[4] * x[8] / (Keq * x[6] * x[3]);
                  (V1 = x[1], V2 = V2, Ka = x[2], Kb = x[3], Kp = x[4], Kq = x[5], Kia = x[6],
                   Kib = x[7], Kiq = x[8], Kip = Keq * x[2] * x[7] / ((x[1] / V2)^2 * x[5]), Keq = Keq))),
    ("ping-pong: 7 free; V2 from (V1/V2)^2*Kp*Kq/(Ka*Kb) = Keq", pingpong, 7, (:A, :B, :P, :Q),
     (x, Keq) -> (V1 = x[1], V2 = x[1] * sqrt(x[4] * x[5] / (Keq * x[2] * x[3])), Ka = x[2], Kb = x[3],
                  Kp = x[4], Kq = x[5], Kia = x[6], Kiq = x[7], Keq = Keq)),
    ("RE random: 7 free; V2 from V1*Kip*Kq/(V2*Kia*Kb) = Keq", re_random, 7, (:A, :B, :P, :Q),
     (x, Keq) -> (V1 = x[1], Kia = x[2], Kib = x[3], Kb = x[4], Kip = x[5], Kq = x[6], Kiq = x[7],
                  V2 = x[1] * x[5] * x[6] / (Keq * x[2] * x[4]))),
]
for (label, rate, n, mets, build) in fams
    @printf("%-62s dimension %d\n", label, dim(rate, build, n, mets, rng))
end
