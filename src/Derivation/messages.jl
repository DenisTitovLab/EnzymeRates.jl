# ABOUTME: The guard's error messages: the cause of each fault, naming each step at fault by
# ABOUTME: its DSL line, and the fix, printed when the fixed variant passes the guard.

# ─── Guard messages ────────────────────────────────────────────────────

"""
    _guard_message(mech, v) -> String

The error for the guard's verdict `v` on `mech`:
`<Mechanism|AllostericMechanism>: <cause, naming each step by its line>. <fix>.` Class 2's
fix is printed when each step's variant leaves no fast turnover cycle
(`_every_step_cures`); every other fix only when its variant passes `_guard`.
"""
function _guard_message(mech::Mechanisms._AnyMech, v::_Verdict)
    who = mech isa Mechanisms.AllostericMechanism ? "AllostericMechanism" : "Mechanism"
    v.fault === _NO_CYCLE &&
        return "$who: no cycle of steps turns substrates into products, so the rate is " *
               "zero everywhere."
    v.fault === _RE_TURNOVER && return _re_turnover_message(who, mech)
    mech isa Mechanisms.AllostericMechanism && v.fault in _ADMISSION_FAULTS &&
        return _admission_message(mech, v)
    v.fault === _NO_FAULT && error("_guard_message: verdict $v holds no fault")
    _face_message(who, mech, v)
end

"""Whether making each step of `cycle`, transitions of `mech`'s ε-network `net`,
steady-state with its kinetic group leaves no fast cycle that turns the reaction over: then
class 2's fix holds whichever step of the cycle is picked."""
function _every_step_cures(mech::Mechanisms._AnyMech, net::_Network, cycle::Vector{Int})
    m = Mechanisms.catalytic(mech)
    for t in cycle
        re = Mechanisms.rapid_equilibrium_steps(m)
        re[Mechanisms.members(m, Mechanisms.group(m, net.edge[t]))] .= false
        isempty(_re_turnover_transitions(_eps_network(mech, re))) || return false
    end
    true
end

"""Class 2's message: the steps of the first fast cycle that turns the reaction over,
each tagged binding among them marked with the conformation where its affinity
vanishes, and the fix where each of those steps made steady-state leaves no fast cycle
that turns the reaction over (`_every_step_cures`)."""
function _re_turnover_message(who::String, mech::Mechanisms._AnyMech)
    m = Mechanisms.catalytic(mech)
    net = _eps_network(mech)
    cycle = _re_turnover_transitions(net)
    c = net.conformation[cycle[1]]
    state = c == 2 ? "inactive" : "active"
    lines = String[]
    for t in cycle
        k = net.edge[t]
        line = Mechanisms._step_text(m, k)
        mech isa Mechanisms.AllostericMechanism && Mechanisms._excludes(mech, k, c) &&
            (line *= " (whose affinity vanishes in the $state conformation)")
        push!(lines, line)
    end
    place = mech isa Mechanisms.AllostericMechanism && c == 2 ?
            " in the inactive conformation" : ""
    msg = "$who: the rapid-equilibrium steps $(join(lines, ", ", " and ")) close a " *
          "cycle$place that turns the reaction over, so the mechanism has no finite rate."
    fix = " Make one step of the cycle steady-state."
    _every_step_cures(mech, net, cycle) ? msg * fix : msg
end

"""The suffix a message writes after the form of node `u` of the two-copy network `net`:
` (active)` or ` (inactive)`."""
_conformation_suffix(net::_Network, u::Int) =
    _copy_of(net, u) == 1 ? " (active)" : " (inactive)"

"""Where a class-3 fault of `mech` lies: nothing for a plain mechanism; the conformation
`c` of an allosteric one, nothing when `c` is 0."""
_place(::Mechanisms.Mechanism, ::Int) = ""
_place(::Mechanisms.AllostericMechanism, c::Int) =
    c == 2 ? " in the inactive conformation" : c == 1 ? " in the active conformation" : ""

"""
The message of a class-3 verdict `v` on `mech`: the face, the forms and steps at fault,
and a fix only where its variant passes the guard (`_steady_state_passes`): for a trap,
the kinetic group of the step that feeds it or joins it to its segment made steady-state;
for a bottomless segment, the group of the first of its rapid-equilibrium bindings whose
variant passes. A verdict in no one conformation labels each form with its own
(`_conformation_suffix`).
"""
function _face_message(who::String, mech::Mechanisms._AnyMech, v::_Verdict)
    m = Mechanisms.catalytic(mech)
    net = _law_network(mech)
    face = _faces(m)[v.face]
    s = _face_structure(net, face)
    at = _face_phrase(m, face, _place(mech, v.conformation))
    # Forms that lie in different conformations name each form's conformation.
    label(u) = v.conformation == 0 ? _conformation_suffix(net, u) : ""
    forms(nodes) =
        "{" * join([Mechanisms._form_text(m, _form_of(net, u)) * label(u) for u in nodes],
                   ", ") * "}"
    lines(edges) = join([Mechanisms._step_text(m, k) for k in edges], ", ", " and ")
    if v.fault === _NONUNIQUE
        C = _away_class(s)
        return "$who: $(at)the forms $(forms(C)) keep enzyme away from free E: every " *
               "step that leaves them takes up a zeroed metabolite " *
               "($(lines(_exits(net, s, C)))), so the steady state is not unique."
    elseif v.fault === _BOTTOMLESS
        bindings = unique!([net.edge[t] for t in eachindex(net.from)
                            if net.fast[t] && net.edge[t] != 0 &&
                               s.segment[net.from[t]] == s.bottomless &&
                               Mechanisms.is_binding(m, net.edge[t])])
        k = findfirst(k -> _steady_state_passes(mech, k), bindings)
        fix = k === nothing ? "" :
              " Make the group of $(Mechanisms._step_text(m, bindings[k])) steady-state."
        return "$who: $(at)the rapid-equilibrium segment " *
               "$(forms(s.segments[s.bottomless])) has no form whose weight keeps the " *
               "segment's lowest power of each zeroed metabolite, so the forms' shares " *
               "depend on how those concentrations vanish and the rate is 0/0 there." * fix
    elseif v.fault === _TRAP
        trap, edge, joined = _trap_feeder(net, s)
        fix = _steady_state_passes(mech, edge) ?
              " Make that step's kinetic group steady-state." : ""
        joined &&
            return "$who: $(at)the forms $(forms(s.fast_classes[trap])) hold enzyme, but " *
                   "the rapid-equilibrium step $(Mechanisms._step_text(m, edge)) that " *
                   "joins them to the rest of their segment takes up a zeroed " *
                   "metabolite, so the rate law misweights them." * fix
        return "$who: $(at)the rapid-equilibrium step $(Mechanisms._step_text(m, edge)) " *
               "fills the forms $(forms(s.fast_classes[trap])), which no " *
               "rapid-equilibrium step running there empties, so the rate law " *
               "misweights them." * fix
    elseif v.fault === _SPLIT
        split = _split_parts(net, s)
        return "$who: $(at)the forms $(join(forms.(split.parts), " and ")) of one " *
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
        return "$who: $(at)$why."
    end
    error("_face_message: no message for $(v.fault)")
end
