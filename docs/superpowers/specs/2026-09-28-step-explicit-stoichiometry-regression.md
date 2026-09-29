# Steps with explicit stoichiometry: regression check

Date: 2026-09-29. This is Verification 1 of `2026-09-28-step-explicit-stoichiometry-design.md`:
every mechanism the enumerator emits must keep its fitted names, rate equations and place in the
enumeration, except the fixtures with fused steps.

## Method

A scratch script (not committed) records, for each mechanism, `fitted_params`, the Reduced and
Full `rate_equation_string`, any error, the stored step and group order (`okey`) and the time
to derive all three. A mechanism is keyed by its steps written with explicit consumed and
released lists, a key that both versions of the code produce. The baseline ran on ab7f0cf,
before the refactor; the comparison run on f7e6394 derived the baseline's mechanisms. Both runs
enumerated independently, and their per-level counts are compared.

Population: `init_mechanisms` plus two levels of the re_to_ss, split and dead_end moves on
reactions R1–R6 of the findings document, every mechanism with V×τ ≤ 337, at most 6,000 derived
per reaction; plus the 42 `MECHANISM_TEST_SPECS`. R6 derives all 1,831 mechanisms of levels 0–1
and a fixed-seed sample of 4,169 at level 2.

## Results

| Set | Distinct per level (both runs) | Derived | fitted | reduced | full | error | okey |
|---|---|---|---|---|---|---|---|
| R1 | 1, 2, 1 | 4 | 0 | 0 | 0 | 0 | 0 |
| R2 | 1, 3, 3 | 7 | 0 | 0 | 0 | 0 | 0 |
| R3 | 1, 4, 6 | 11 | 0 | 0 | 0 | 0 | 0 |
| R4 | 62, 369, 1388 | 1,819 | 0 | 0 | 0 | 0 | 0 |
| R5 | 62, 719, 3934 | 4,715 | 0 | 2 | 2 | 0 | 0 |
| R6 | 62, 1769, 28304 | 6,000 | 0 | 0 | 0 | 0 | 0 |
| specs | – | 42 | 12 | 12 | 12 | 4 | – |

The difference columns count mechanisms whose field differs. Both runs enumerated without
errors, and neither run lacks a mechanism of the other.

## Differences and verdicts

### R5: two fixed bugs

Two level-2 dead_end children of ping-pong seeds differ in their Reduced and Full strings;
their fitted names and step order are identical. Each carries a second Q binding onto a
substrate-bound form, E(A) + Q → E(A, Q), in the kinetic group of E + Q → E(Q); the Q group is
SS and the B group RE in one, the reverse in the other. The A binding group joins E(A, Q) to
E(Q) at rapid equilibrium, so E(A, Q) → E(A) + Q is a second Q-release route. The baseline's
reaction-cut search skipped every step onto a complex holding both a substrate and a product,
took the one-step cut {E + Q → E(Q)}, and so dropped that route.

The oracle is brute-force mass action over every enzyme form (BigFloat, rapid-equilibrium steps
run at rate scale 1e60), at three random points per mechanism. It was compared with both runs'
Full strings, with every constant drawn independently (the Full law holds for any values), and
with their Reduced strings, with the fitted constants drawn and the dependents computed by each
string's own constraint lines.

| Mechanism | Refactor: largest relative error | Baseline: relative error |
|---|---|---|
| Q group RE, B group SS | 1.9e-59 | 0.098 to 0.82 |
| Q group SS, B group RE | 5.3e-60 | 0.35 to 0.69 |

The refactor's law matches mass action to the precision the finite rapid-equilibrium rate
allows; the baseline's does not. Both are fixed bugs.

### R6: eight more cases outside the sample

A comparison of the two numerator polynomials over all 30,135 R6 mechanisms, run during the
implementation, found eight more of the same shape (level-2 dead_end children; the inhibitor
copy is A, B, P or Q). None falls in the level-2 sample, so the table shows no R6 difference.
Re-enumerated and derived on f7e6394, all eight, and the two R5 mechanisms, match the oracle to
1.3e-58 or better in both Full and Reduced. The baseline's numerator disagreed with the oracle
at all 24 points drawn for the eight.

### specs: twelve renamed fixtures

Twelve fixtures write chemistry and product release as one step. Such a step is a
transformation, named by its two forms:

| Spec | Renames |
|---|---|
| Uni-Uni | kon_P_ES → k_ES_to_E, koff_P_ES → k_E_to_ES |
| RE Uni-Uni | kon_P_EA → k_EA_to_E, koff_P_EA → k_E_to_EA |
| Segel Ordered Uni Bi | kon_P_EA → k_EA_to_EQ, koff_P_EA → k_EQ_to_EA |
| Segel Ordered Bi Bi, RE Ordered Bi-Bi | kon_P_EAB → k_EAB_to_EQ, koff_P_EAB → k_EQ_to_EAB |
| Segel Ping Pong Bi Bi | kon_P_EA → k_EA_to_F, koff_P_EA → k_F_to_EA, kon_Q_FB → k_FB_to_E, koff_Q_FB → k_E_to_FB |
| Segel Ordered Ter Bi | kon_P_EABC → k_EABC_to_EQ, koff_P_EABC → k_EQ_to_EABC |
| Segel Ordered Ter Ter | kon_P_EABC → k_EABC_to_EQR, koff_P_EABC → k_EQR_to_EABC |
| Segel Bi Uni Uni Uni PP Ter Bi | kon_P_EAB → k_EAB_to_F, koff_P_EAB → k_F_to_EAB, kon_Q_FC → k_FC_to_E, koff_Q_FC → k_E_to_FC |
| Segel Bi Uni Uni Bi PP Ter Ter | kon_P_EAB → k_EAB_to_F, koff_P_EAB → k_F_to_EAB, kon_Q_FC → k_FC_to_ER, koff_Q_FC → k_ER_to_FC |
| Segel Bi Bi Uni Uni PP Ter Ter | kon_P_EAB → k_EAB_to_FQ, koff_P_EAB → k_FQ_to_EAB, kon_R_FC → k_FC_to_E, koff_R_FC → k_E_to_FC |
| Segel Hexa Uni Ping Pong | kon_P_EA → k_EA_to_F, koff_P_EA → k_F_to_EA, kon_Q_FB → k_FB_to_G, koff_Q_FB → k_G_to_FB, kon_R_GC → k_GC_to_E, koff_R_GC → k_E_to_GC |

For each, the renamed baseline fitted set equals the refactor's, and the two Reduced laws agree
at five random points to a relative 7.6e-15 or better. No other spec's names changed.

### specs: four error strings

The script asks every spec for a Full string, which `AllostericEnzymeMechanism` does not
provide, so all 12 allosteric specs record the same `MethodError` in both runs. The message
prints the mechanism's type, whose step encoding the refactor changes to (from, to, consumed
tuple, released tuple, flag), as design section 3 states. In four specs the message's
300-character truncation reaches the first step, so the text differs. Their fitted names and
Reduced strings are identical. This is an artefact of the script, not a change in behaviour.

## Timing

Median derivation time per mechanism, in seconds (`compile_mechanism`, `fitted_params` and both
strings):

| Set | Baseline | Refactor | Change |
|---|---|---|---|
| R1 | 0.126 | 0.107 | −15% |
| R3 | 0.120 | 0.083 | −31% |
| R4 | 0.570 | 0.437 | −23% |
| R5 | 0.663 | 0.490 | −26% |
| R6 | 1.255 | 0.906 | −28% |
| specs | 0.230 | 0.177 | −23% |

R2's median is 0.001 s in both runs. Total derivation time over all 12,598 mechanisms fell from
12,779 s to 9,617 s. The design's limit is a 20% rise; no set rose. The two runs ran hours
apart on the same machine, so part of the fall may be load; its sign is the same in every set.

Enumeration is slower. Levels 0–2 of R6 took 7.3 s in the baseline run and 8.6 s in the
refactor run (+18%; R5 1.1 s to 1.2 s, R4 0.4 s in both). A warm measurement at 9a404e0 gave
R6 7.32 s to 8.27 s (+13%); there the parameter-name collision guard cost 0.37 s and
step-direction canonicalization 0.31 s.

## Conclusion

Every enumerated mechanism compared keeps its fitted names, both rate equations, its step order
and its place in the enumeration, except two R5 mechanisms whose baseline law was wrong. The
fixtures with fused steps change only their parameter names. Derivation is not slower.
