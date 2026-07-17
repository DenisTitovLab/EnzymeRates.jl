# Prototype scripts — solve-then-limit MWC derivation

Scratch validation scripts referenced by the design spec
(`docs/superpowers/specs/2026-07-16-mwc-solve-then-limit-derivation-design.md`) and
plan (`docs/superpowers/plans/2026-07-16-mwc-solve-then-limit-derivation.md`).
They are NOT part of the package; they run standalone against `EnzymeRates` to
validate the redesign's claims against an independent mass-action oracle.

Run: `cd <repo> && julia --project=. docs/superpowers/prototypes/2026-07-16-mwc-solve-then-limit/<script>.jl`

## Core validation
- `solve_then_limit.jl` — implements solve-then-limit for n=1 uni-uni; matches the ground truth on all constructable tag combos.
- `n2_spotcheck.jl` — n=2 Family A vs an explicit 2-protomer network.
- `attack.jl`, `attack_claimB_strong.jl` — adversarial checks (note: attack.jl's "cross-weight 8-37% off" used a bogus formula; the shipped cross-weight is correct — see the spec's corrected motivation).
- `brief_confirm_trap.jl`, `brief_enum_check.jl`, `brief_all_constructable.jl` — the correctness boundary equals the constructor guard.

## Normalization (the crux)
- `redesign_norm_vs_skip.jl`, `redesign_pingpong.jl` — per-state-normalized vs un-normalized combine; ping-pong `oracle_Eflip`.
- `ldh_koff_ambiguity.jl`, `eq_consistency.jl` — the `:OnlyA` SS limit convention (only `k_I→0` is thermodynamically legal).
- `final_check.jl`, `pp_eqcheck.jl`, `limit_vs_deletion.jl` — supporting checks.
- `compare_pingpong.jl`, `probe_mech.jl` — the shipped `rate_equation` vs the ping-pong oracle (3e-16; the deferred `:NonequalAI` ping-pong case).

## Ping-pong :OnlyA kcat bug root cause
- `RC_gt.jl` — n=1 two-conformation ground truth for the two PFK ping-pong `:OnlyA` mechanisms (both valid).
- `RC_probe.jl` — the inactive-graph fragmentation (two empty-bound forms).
- `RC_dfree.jl` — the `d_free_I` values (0 for one, product-bearing for the other).
See `docs/superpowers/findings/2026-07-16-pingpong-onlya-kcat-bug.md`.
