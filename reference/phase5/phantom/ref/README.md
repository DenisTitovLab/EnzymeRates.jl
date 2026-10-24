# Reference material for the reaction-graph rewrite

This branch archives the prototypes, probe reports and decision records behind
`docs/superpowers/specs/2026-10-12-reaction-graph-rewrite-design.md`. It is reference only:
none of this code is package code, and the branch is never merged. Scripts refer to the
session scratchpad they ran in; the generated populations (about 400 MB of serialized
mechanisms and laws) are left out because the tests regenerate them from
`init_mechanisms` and `expand_mechanisms`.

## Layout

| Path | What it holds | Spec sections |
|---|---|---|
| `decisions/` | The decision page (`graph-redesign-decisions.html`, cards A1–E11 with evidence and options), Denis's choices and notes (`cards/*.json`), the replies on the cards (`replies/*.json`) and the decisions made in chat (`chat_decisions_ledger.md`) | Appendix A |
| `slice/src/` | The vertical-slice prototype of the whole design: value type, names, polynomials, joint solve, parameter API, engine, normalization, printing, kcat, allosteric law, admission, code generation, enumeration moves, key, pipeline | 4–13 |
| `slice/INTERFACES.md` | The slice's interface contract, stage by stage, with the skeptic's and critic's notes | 4–13 |
| `slice/reports/` | Stage reports: core, solve, engine, allosteric, codegen, enumeration, integration, skeptic, critic | Appendix B |
| `slice/<stage>/` | Each stage's check scripts and oracles | 14 |
| `probes/` | The derisking probes' reports (0–12; `12.md` is the critic) | Appendix B |
| `explore/` | The first exploration of the package (DSL, enumeration, consumers, tests, docs, bugs, MWC, code generation, literature) | 1 |
| `followups/` | Follow-up investigations behind cards B1, B2, B5, C4, C5 | 7, 9 |
| `proto/` | The probes' throwaway code (engine, joint solve, limits, normalization, mechanism type, code generation, key) | — |
| `round4/reports/` | Verification of the round-4 answers: kcat census, well-posedness guard, flip at E, their verifiers and critic, and the random-order peak search | 7, 8, 12 |
| `round4/kcat/`, `round4/guard/`, `round4/flipE/` | Their scripts; `round4/guard/lib.jl` holds the RE-speed oracle and the face check, `round4/flipE/scripts/lib.jl` the two-conformation mass-action oracle, `round4/flipE/prototype_flip_at_E.diff` the flip-at-E change to PR #86's branch | 7, 8, 14 |
| `census/` | The enumeration census (small and large reactions, fit-cost calibration), its verifiers and critic | 13 |
| `peaksearch/` | The adversarial peak search on random-order bi-bi and variants | 12 |

## Provenance

Produced between 2026-10-09 and 2026-10-12 in one Claude Code session with Denis, by
workflows of probe, verifier and critic agents. Each report states what it checked, its
counts and what it did not check; later reports correct earlier ones where they disagree
(for example `round4/reports/verify_*.md` and `census/reports/critic.md`).
