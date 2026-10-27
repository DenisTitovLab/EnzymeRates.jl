# ABOUTME: The derivation memo: the mechanism a task last asked about, with its guard
# ABOUTME: verdict, joint solve and law, so consecutive questions about it derive once.

"""
The derivation memo of a task (`_memo`): the concrete mechanism `m` last asked about, and
its guard verdict (`_guard`), joint solve (`_solve`) and law (`_law`), each filled on
first need. Identify's PASS 1 asks one mechanism several questions in a row
(`check_valid_rate_equation`, `fitted_params`, `rate_equation_string`,
`rate_equation_key`); through the memo they run the guard, the solve and the law once.
The memo holds one mechanism, replaced by the next one asked about (compared by `===`),
and lives in the task's local storage, so tasks share nothing. The compiled types'
derivation (`_derivation_input`, `rate_equation`) never reads it.
"""
mutable struct _Memo
    m::Any
    verdict::Union{Nothing, _Verdict}
    solve::Union{Nothing, _Solve}
    law::Union{Nothing, _Law, _AlloLaw}
end

"""The task's memo (`_Memo`), emptied and given to `m` unless it holds `m`."""
function _memo(m::Mechanisms._AnyMech)
    c = get!(() -> _Memo(nothing, nothing, nothing, nothing), task_local_storage(),
             _Memo)::_Memo
    c.m === m || (c.m = m; c.verdict = nothing; c.solve = nothing; c.law = nothing)
    c
end

"""`m`'s guard verdict (`_guard`), through the memo."""
function _memo_verdict(m::Mechanisms._AnyMech)
    v = _memo(m).verdict
    v === nothing || return v
    _memo(m).verdict = _guard(m)
end

"""`m`'s joint solve (`_solve`), through the memo."""
function _memo_solve(m::Mechanisms._AnyMech)
    s = _memo(m).solve
    s === nothing || return s
    _memo(m).solve = _solve(m)
end

"""`m`'s law (`_law`) from its joint solve, through the memo: a `_Law` for a plain
mechanism, an `_AlloLaw` for an allosteric one."""
function _memo_law(m::Mechanisms._AnyMech)
    L = _memo(m).law
    L === nothing || return L::(m isa Mechanisms.Mechanism ? _Law : _AlloLaw)
    _memo(m).law = _law(m, _memo_solve(m))
end

"""`f` called on `m`, its joint solve `s` and its law `L`, both through the memo: the
memo's counterpart of `_with_law`."""
_with_memo_law(f, m::Mechanisms._AnyMech) = f(m, _memo_solve(m), _memo_law(m))
