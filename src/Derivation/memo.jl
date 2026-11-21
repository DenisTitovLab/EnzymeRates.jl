# ABOUTME: The derivation memo: the mechanism a task last asked about, with its guard
# ABOUTME: verdict, joint solve and law, so consecutive questions about it derive once.

"""
The derivation memo of a task (`_memo`): the concrete mechanism `m` last asked about, and
its guard verdict (`_guard`), joint solve (`_solve`) and law (`_law`), each filled on
first need. Identify's PASS 1 asks one mechanism three questions in a row (`parameters`,
`rate_equation_string`, `rate_equation_dedup_key`); through the memo they run the solve
and the law once.
The memo holds one mechanism, replaced by the next one asked about (compared by `===`),
and lives in the task's local storage, so tasks share nothing. A compiled rate law's
generated methods (`compile_rate_equation`, `_derivation_input`) never read it.
"""
mutable struct _Memo
    m::Any
    verdict::Union{Nothing, _Verdict}
    solve::Union{Nothing, _Solve}
    law::Union{Nothing, _Law, _AlloLaw}
end

"""The task's memo (`_Memo`), emptied and given to `m` unless it holds `m`."""
function _memo(m::_AnyMech)
    c = get!(() -> _Memo(nothing, nothing, nothing, nothing), task_local_storage(),
             _Memo)::_Memo
    c.m === m || (c.m = m; c.verdict = nothing; c.solve = nothing; c.law = nothing)
    c
end

"""`m`'s guard verdict (`_guard`), through the memo."""
function _memo_verdict(m::_AnyMech)
    v = _memo(m).verdict
    v === nothing || return v
    _memo(m).verdict = _guard(m)
end

"""`m`'s joint solve (`_solve`), through the memo."""
function _memo_solve(m::_AnyMech)
    s = _memo(m).solve
    s === nothing || return s
    _memo(m).solve = _solve(m)
end

"""`m`'s law (`_law`) from its joint solve, through the memo: a `_Law` for a plain
mechanism, an `_AlloLaw` for an allosteric one."""
function _memo_law(m::_AnyMech)
    L = _memo(m).law
    L === nothing || return L::(m isa Mechanisms.Mechanism ? _Law : _AlloLaw)
    _memo(m).law = _law(m, _memo_solve(m))
end

"""`f` called on `m`, its joint solve `s` and its law `L`, both through the memo: the
memo's counterpart of `_with_law`."""
_with_memo_law(f, m::_AnyMech) = f(m, _memo_solve(m), _memo_law(m))
