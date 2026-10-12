# ABOUTME: Rebuilds as Pop.PMech the hand-built adversarial mechanisms that exist only as code:
# ABOUTME: Skeptic P1's D1–D5 (adv_direct.jl) and Skeptic P2's five `hand:` cases (adv/*.jl).

"""Skeptic P1's direct mechanisms D1–D5, as `adv_direct.jl` defines them."""
function _testhelper_p1_direct()
    F(name, mets...; conf = :E) = Pop.PSpecies([Pop.PMet(x, r) for (x, r) in mets], conf,
                                               Symbol[], Symbol[], name)
    M(x, r = :S) = Pop.PMet(x, r)
    function mech(id, subs, prods, forms, steps; inh = Symbol[])
        idx = Dict(f.name => i for (i, f) in enumerate(forms))
        st = [Pop.PStep(idx[a], idx[b], c, r, re) for (a, b, c, r, re) in steps]
        Pop.PMech(; id, source = id, substrates = subs, products = prods, inhibitors = inh,
                  metabolites = [subs; prods; inh], species = forms, steps = st,
                  groups = [[j] for j in eachindex(st)])
    end
    [mech("D1_free_form_not_root", [:S], [:P],
          [F(:E), F(:ES, (:S, :S)), F(:D; conf = :D)],
          [(:E, :ES, [M(:S)], [], false), (:ES, :E, [], [M(:P, :P)], false),
           (:D, :E, [], [], true)]),
     mech("D2_bridge_between_reactive_parts", [:S], [:P],
          [F(:E), F(:ES, (:S, :S)), F(:F; conf = :F), F(:FS, (:S, :S); conf = :F)],
          [(:E, :ES, [M(:S)], [], false), (:ES, :E, [], [M(:P, :P)], false),
           (:F, :FS, [M(:S)], [], false), (:FS, :F, [], [M(:P, :P)], false),
           (:E, :F, [], [], false)]),
     mech("D3_single_segment_core_plus_ss_leaf", [:S], [:P],
          [F(:E), F(:ES, (:S, :S)), F(:EP, (:P, :P)), F(:EI, (:I, :Inh))],
          [(:E, :ES, [M(:S)], [], true), (:ES, :EP, [], [], false),
           (:EP, :E, [], [M(:P, :P)], true), (:E, :EI, [M(:I, :Inh)], [], false)];
          inh = [:I]),
     mech("D4_nonreactive_square_between_reactive_parts", [:S], [:P],
          [F(:E), F(:ES, (:S, :S)), F(:F; conf = :F), F(:FS, (:S, :S); conf = :F),
           F(:G; conf = :G), F(:H; conf = :H)],
          [(:E, :ES, [M(:S)], [], false), (:ES, :E, [], [M(:P, :P)], false),
           (:F, :FS, [M(:S)], [], false), (:FS, :F, [], [M(:P, :P)], false),
           (:E, :G, [], [], false), (:G, :F, [], [], false), (:E, :H, [], [], false),
           (:H, :F, [], [], false)]),
     mech("D5_parallel_rank_21", [:S], [:P],
          [F(:E); F(:ES, (:S, :S)); [F(Symbol("EI", k), (Symbol("I", k), :Inh)) for k in 1:20];
           [F(Symbol("ESI", k), (Symbol("I", k), :Inh), (:S, :S)) for k in 1:20]],
          [(:E, :ES, [M(:S)], [], false); (:ES, :E, [], [M(:P, :P)], false);
           [(:E, Symbol("EI", k), [M(Symbol("I", k), :Inh)], [], true) for k in 1:20];
           [(:ES, Symbol("ESI", k), [M(Symbol("I", k), :Inh)], [], true) for k in 1:20];
           [(Symbol("EI", k), Symbol("ESI", k), [M(:S)], [], false) for k in 1:20]];
          inh = [Symbol("I", k) for k in 1:20])]
end

"""Skeptic P2's hand-built mechanisms (`adversarial.jl`, `followup.jl`, `last_checks.jl`)."""
function _testhelper_p2_hand()
    M(n, r) = Pop.PMet(n, r)
    S_, P_, I_, J_, NONE = M(:S, :S), M(:P, :P), M(:I, :Inh), M(:J, :Inh), Pop.PMet[]
    sp(conf, b...) = Pop.PSpecies(collect(b), conf, Symbol[], Symbol[],
                                  Symbol(conf, join(string.(getfield.(collect(b), :name)))))
    mk(id, species, steps, groups; inh = Symbol[]) = Pop.PMech(; id, source = id,
        substrates = [:S], products = [:P], inhibitors = inh, metabolites = [:S; :P; inh],
        species, steps, groups)
    [mk("hand:parallel_identical_steps", [sp(:E), sp(:E, S_), sp(:E, P_)],
        [Pop.PStep(1, 2, [S_], NONE, true), Pop.PStep(2, 3, NONE, NONE, false),
         Pop.PStep(2, 3, NONE, NONE, false), Pop.PStep(1, 3, [P_], NONE, true)],
        [[1], [2], [3], [4]]),
     mk("hand:group_different_conformational_changes",
        [sp(:E), sp(:E, S_), sp(:E, P_), sp(:F), sp(:G)],
        [Pop.PStep(1, 2, [S_], NONE, false), Pop.PStep(2, 3, NONE, NONE, false),
         Pop.PStep(1, 3, [P_], NONE, false), Pop.PStep(1, 4, NONE, NONE, false),
         Pop.PStep(4, 5, NONE, NONE, false), Pop.PStep(5, 1, NONE, NONE, false)],
        [[1], [2], [3], [4, 5], [6]]),
     mk("hand:half_keq", [sp(:E), sp(:F), sp(:E, S_), sp(:F, S_)],
        [Pop.PStep(3, 4, NONE, NONE, false), Pop.PStep(1, 2, NONE, NONE, false),
         Pop.PStep(4, 1, NONE, [P_], false), Pop.PStep(2, 3, [S_], NONE, false)],
        [[1, 2], [3], [4]]),
     mk("hand:re_dead_end_chain", [sp(:E), sp(:E, S_), sp(:E, P_), sp(:E, I_), sp(:F, I_)],
        [Pop.PStep(1, 2, [S_], NONE, true), Pop.PStep(2, 3, NONE, NONE, false),
         Pop.PStep(1, 3, [P_], NONE, true), Pop.PStep(1, 4, [I_], NONE, true),
         Pop.PStep(4, 5, NONE, NONE, true)], [[1], [2], [3], [4], [5]]; inh = [:I]),
     Pop.PMech(; id = "hand:crossed_reps", source = "hand:crossed_reps", substrates = [:S],
               products = [:P], inhibitors = [:I, :J], metabolites = [:S, :P, :I, :J],
               species = [sp(:E), sp(:E, S_), sp(:E, P_), sp(:F), sp(:E, J_),
                          sp(:F, J_, S_), sp(:F, I_), sp(:E, I_, S_)],
               steps = [Pop.PStep(1, 2, [S_], NONE, true), Pop.PStep(2, 3, NONE, NONE, false),
                        Pop.PStep(1, 3, [P_], NONE, true), Pop.PStep(1, 4, NONE, NONE, true),
                        Pop.PStep(1, 5, [J_], NONE, true), Pop.PStep(4, 7, [I_], NONE, true),
                        Pop.PStep(5, 6, [S_], NONE, false), Pop.PStep(7, 8, [S_], NONE, false)],
               groups = [[1], [2], [3], [4], [5], [6], [7, 8]])]
end
