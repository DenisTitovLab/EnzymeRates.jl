# ABOUTME: Runs a small identify_rate_equation search (the package test's smoke settings, NW
# ABOUTME: Distributed workers) so its compiled code can be scanned for JuliaLICM hoists.
using TestEnv; TestEnv.activate()
using Distributed
const NW = parse(Int, get(ENV, "NW", "0"))
NW > 0 && addprocs(NW; exeflags = ["--project=$(Base.active_project())", "--heap-size-hint=700M"])
@everywhere using EnzymeRates, Random, OptimizationCMAEvolutionStrategy
using EnzymeRates, Random, OptimizationCMAEvolutionStrategy

rxn = @enzyme_reaction begin
    substrates: S[C]
    products: P[C]
    competitive_inhibitors: R
end
base = first(EnzymeRates.init_mechanisms(rxn))
tags = [EnzymeRates.bound_metabolite(first(g)) isa EnzymeRates.Reactant ? :OnlyA : :NonequalAI
        for g in EnzymeRates.steps(base)]
site = EnzymeRates.RegulatorySite([EnzymeRates.AllostericRegulator(:R)], 1, [:OnlyI])
truth = EnzymeRates.AllostericEnzymeMechanism(EnzymeRates.AllostericMechanism(
    EnzymeRates.reaction(base), copy(EnzymeRates.steps(base)), tags, 1, [site]))
params = (K_A_EP_to_E_P = 1.0, K_A_ES_to_E_S = 0.5, k_A_ES_to_EP = 5.0, K_I_Rreg = 2.0,
          L = 0.1, Keq = 2.0, E_total = 1.0)
Random.seed!(42)
rows = [(g, 0.1 + 9.9rand(), 0.1 + 9.9rand(), 0.1 + 9.9rand()) for g in 1:5 for _ in 1:10]
data = (group = ["G$(r[1])" for r in rows], S = [r[2] for r in rows], P = [r[3] for r in rows],
        R = [r[4] for r in rows],
        Rate = [rate_equation(truth, (S = r[2], P = r[3], R = r[4]), params) for r in rows])
prob = IdentifyRateEquationProblem(rxn, data; Keq = 2.0)
for round in 1:parse(Int, get(ENV, "ROUNDS", "1"))
    res = redirect_stdout(devnull) do
        identify_rate_equation(prob; min_beam_width = 1, loss_rel_threshold = 1.0,
            loss_abs_threshold = 0.0, max_param_count = 8, n_cv_candidates = 1,
            save_dir = mktempdir(), optimizer = CMAEvolutionStrategyOpt(),
            n_restarts = 1, maxtime = 1.0, show_progress = false)
    end
    println("round ", round, ": ", size(res.cv_results, 1), " cv rows"); flush(stdout)
    GC.gc(true)
    NW > 0 && @everywhere GC.gc(true)
    println("full gc ", round, " ok"); flush(stdout)
end
println("finished")
