# ABOUTME: The guard's error messages: the cause of each fault, naming each step at fault by
# ABOUTME: its DSL line, and the fix, printed when the fixed variant passes the guard.

# ─── Guard messages ────────────────────────────────────────────────────

"""
    _guard_message(mech) -> String

The error for the guard's verdict on `mech`:
`Mechanism: <cause, naming each step by its line>. <fix>.`, the
cause read under the reading the verdict was judged under (`_guard_with_reading`). Class
2's fix is printed when each step's variant leaves no fast turnover cycle
(`_every_step_cures`); every other fix only when its variant passes `_guard`. When two
constants share a name (`_SHARED_NAME`), the error names them and asks for a renaming;
for a cycle that does not balance a metabolite (`_UNBALANCED_CYCLE`), it names the
cycle's steps and its net change (`_unbalanced_cycle`) and the slips that make one.
An error when the guard admits `mech`.
"""
function _guard_message(mech::Mechanisms.Mechanism)
    (; verdict, reading) = _guard_with_reading(mech)
    verdict.fault === _NO_FAULT && error("_guard_message: verdict $verdict holds no fault")
    verdict.fault === _SHARED_NAME &&
        return "Mechanism: constants of this mechanism's rate equation share a name (" *
               join(_shared_constant_names(mech), ", ") * "). Rename a metabolite or an " *
               "enzyme form so that each constant has a name of its own."
    if verdict.fault === _UNBALANCED_CYCLE
        m = mech.catalytic
        (; steps, change) = _unbalanced_cycle(m)
        lines = join([Mechanisms.step_text(m, k) for k in steps], ", ", " and ")
        # The metabolites on one side of the net change, each with its count when above 1.
        function side(sgn)
            listed = [(abs(c) == 1 ? "" : "$(abs(c)) ") * String(m.metabolites[x])
                      for (x, c) in enumerate(change) if sgn * c > 0]
            isempty(listed) ? "nothing" : join(listed, ", ", " and ")
        end
        return "Mechanism: the steps $lines close a cycle that, net, takes up $(side(1)) " *
               "and gives off $(side(-1)). Round every cycle, the substrates go in and " *
               "the products come out as often as the reaction runs, and each modifier " *
               "comes off as often as it binds. Check whether one of these steps leaves " *
               "out a metabolite, writes one on the wrong side or names the wrong form."
    end
    verdict.fault === _NO_CYCLE &&
        return "Mechanism: no cycle of steps turns substrates into products, so the rate " *
               "is zero everywhere."
    verdict.fault === _RE_TURNOVER && return _re_turnover_message(mech, reading)
    # Only an allosteric mechanism excludes steps, so only its verdict is admission's.
    verdict.fault in _ADMISSION_FAULTS && return _admission_message(mech, verdict, reading)
    _face_message(mech, verdict)
end

"""Whether making each step of `cycle`, transitions of `mech`'s ε-network `net` under
`reading`, steady-state with its kinetic group leaves no fast cycle that turns the
reaction over: then class 2's fix holds whichever step of the cycle is picked."""
function _every_step_cures(mech::Mechanisms.Mechanism, net::_Network, cycle::Vector{Int},
                           reading::_Reading)
    m = mech.catalytic
    groups = _kinetic_groups(m)
    for t in cycle
        re = _re_flags(m)
        re[groups[m.steps[net.edge[t]].kinetic_group]] .= false
        isempty(_re_turnover_transitions(_eps_network(mech, re, reading))) || return false
    end
    true
end

"""Class 2's message, the excluded steps switched off as `reading` says: the steps of the
first fast cycle that turns the reaction over, each excluded step among them marked with
the conformation where it loses a direction (its affinity, when it takes up one
metabolite and loses that uptake), and the fix where each of those steps made
steady-state leaves no fast cycle that turns the reaction over (`_every_step_cures`)."""
function _re_turnover_message(mech::Mechanisms.Mechanism, reading::_Reading)
    m = mech.catalytic
    net = _eps_network(mech, _re_flags(m), reading)
    cycle = _re_turnover_transitions(net)
    c = net.conformation[cycle[1]]
    state = c == 2 ? "inactive" : "active"
    lines = String[]
    for t in cycle
        k = net.edge[t]
        line = Mechanisms.step_text(m, k)
        if _excluded_steps(mech, c)[k]
            e = m.steps[k]
            affinity = reading[k, c] === _LOSES_FORWARD && length(e.taken_up) == 1 &&
                       isempty(e.given_off)
            line *= affinity ? " (whose affinity vanishes in the $state conformation)" :
                    " (whose rate one way vanishes in the $state conformation)"
        end
        push!(lines, line)
    end
    # Only an allosteric mechanism's network has an inactive copy.
    place = c == 2 ? " in the inactive conformation" : ""
    msg = "Mechanism: the rapid-equilibrium steps $(join(lines, ", ", " and ")) close a " *
          "cycle$place that turns the reaction over, so the mechanism has no finite rate."
    fix = " Make one step of the cycle steady-state."
    _every_step_cures(mech, net, cycle, reading) ? msg * fix : msg
end

"""The suffix a message writes after the form of node `u` of the two-copy network `net`:
` (active)` or ` (inactive)`."""
_conformation_suffix(net::_Network, u::Int) =
    _copy_of(net, u) == 1 ? " (active)" : " (inactive)"

"""Where a class-3 fault of `mech` lies: nothing for a plain mechanism; the conformation
`c` of an allosteric one, nothing when `c` is 0."""
_place(mech::Mechanisms.Mechanism, c::Int) =
    !_is_allosteric(mech) ? "" : c == 2 ? " in the inactive conformation" :
    c == 1 ? " in the active conformation" : ""

"""
The message of a class-3 verdict `v` on `mech`: the face, the forms and steps at fault,
and a fix only where its variant passes the guard (`_steady_state_passes`): for a trap,
the kinetic group of the step that feeds it or joins it to its segment made steady-state;
for a bottomless segment, the kinetic group of the first of its rapid-equilibrium bindings
whose variant passes. A verdict in no one conformation labels each form with its own
(`_conformation_suffix`).
"""
function _face_message(mech::Mechanisms.Mechanism, v::_Verdict)
    m = mech.catalytic
    net = _law_network(mech)
    face = _faces(m)[v.face]
    s = _face_structure(net, face)
    at = _face_phrase(m, face, _place(mech, v.conformation))
    # Forms that lie in different conformations name each form's conformation.
    label(u) = v.conformation == 0 ? _conformation_suffix(net, u) : ""
    forms(nodes) =
        "{" * join([String(m.enzyme_forms[_form_of(net, u)]) * label(u) for u in nodes],
                   ", ") * "}"
    lines(edges) = join([Mechanisms.step_text(m, k) for k in edges], ", ", " and ")
    if v.fault === _NONUNIQUE
        C = _away_class(s)
        return "Mechanism: $(at)the forms $(forms(C)) keep enzyme away from free E: " *
               "every step that leaves them takes up a zeroed metabolite " *
               "($(lines(_exits(net, s, C)))), so the steady state is not unique."
    elseif v.fault === _BOTTOMLESS
        binding(k) = length(m.steps[k].taken_up) + length(m.steps[k].given_off) == 1
        bindings = unique!([net.edge[t] for t in eachindex(net.from)
                            if net.fast[t] && net.edge[t] != 0 &&
                               s.segment[net.from[t]] == s.bottomless &&
                               binding(net.edge[t])])
        k = findfirst(k -> _steady_state_passes(mech, k), bindings)
        fix = k === nothing ? "" :
              " Make the kinetic group of $(Mechanisms.step_text(m, bindings[k])) " *
              "steady-state."
        return "Mechanism: $(at)the rapid-equilibrium segment " *
               "$(forms(s.segments[s.bottomless])) has no form whose weight keeps the " *
               "segment's lowest power of each zeroed metabolite, so the forms' shares " *
               "depend on how those concentrations vanish and the rate is 0/0 there." * fix
    elseif v.fault === _TRAP
        trap, edge, joined = _trap_feeder(net, s)
        fix = _steady_state_passes(mech, edge) ?
              " Make that step's kinetic group steady-state." : ""
        joined &&
            return "Mechanism: $(at)the forms $(forms(s.fast_classes[trap])) hold " *
                   "enzyme, but the rapid-equilibrium step " *
                   "$(Mechanisms.step_text(m, edge)) that joins them to the rest of " *
                   "their segment takes up a zeroed metabolite, so the rate law " *
                   "misweights them." * fix
        return "Mechanism: $(at)the rapid-equilibrium step " *
               "$(Mechanisms.step_text(m, edge)) fills the forms " *
               "$(forms(s.fast_classes[trap])), which no rapid-equilibrium step running " *
               "there empties, so the rate law misweights them." * fix
    elseif v.fault === _SPLIT
        split = _split_parts(net, s)
        return "Mechanism: $(at)the forms $(join(forms.(split.parts), " and ")) of one " *
               "rapid-equilibrium segment are joined only by steps that take up a " *
               "zeroed metabolite ($(lines(split.joins))), so the rate law cannot weight " *
               "one against the other."
    elseif v.fault === _LAW_MULTI || v.fault === _MULTILEVEL
        why = v.limit === _LIMIT_LAW_NAN ?
                  "the rate law has no value (0/0), while the steady state does" :
              v.limit === _LIMIT_SPEED ?
                  "the steady state depends on how fast the rapid-equilibrium steps " *
                  "run, which the rate law does not carry" :
                  "the steady state's fast-equilibrium limit differs from the rate law"
        return "Mechanism: $(at)$why."
    end
    error("_face_message: no message for $(v.fault)")
end
