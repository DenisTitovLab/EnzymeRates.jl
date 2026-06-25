# kcat rescaling at measured concentrations (future enhancement)

Status: idea, deferred to a separate PR.

## Context

`rescale_parameter_values` rescales a mechanism's SS rate constants so that
`_kcat_forward(m, params) ≈ scale_k_to_kcat`. The default `scale_k_to_kcat = 1.0`
is a *relative* normalization anchor — its absolute value is arbitrary, so the
analytic `_kcat_forward` (the extrapolated productive turnover at substrate
saturation) is fine for it.

When a user supplies a **real measured kcat** instead of `1.0`, the rescaling
target becomes a physical number, and the model's kcat must be computed the same
way the experimentalist obtained theirs. The analytic saturating-limit
`_kcat_forward` does not generally match a measured value for enzymes with
**substrate inhibition** or other **non-monotonic kinetics**: the analytic value
extrapolates to a saturating productive plateau that the real `v([S])` curve never
reaches (it peaks at finite substrate and declines).

## Proposal

Let `rescale_parameter_values` accept the **concentrations of the metabolites at
which the kcat was measured**, and anchor to the model rate *there* rather than to
the analytic saturating limit:

- rescale the SS rate constants so that `rate_equation(m, measured_concs, params)`
  equals the provided kcat, instead of solving `_kcat_forward(m, params) = kcat`.

Evaluating `rate_equation` at the actual assay concentrations is unambiguous and
needs no saturating-limit extraction, so it avoids every substrate-inhibition /
non-monotonicity edge case the analytic path has to reason about (dead-end
substrate inhibition, V-type allosteric self-inhibition, product-activator
corners). It also matches how kcat is measured in practice.

## Scope

Separate PR. Independent of the current allosteric fixes (A: dead-I-state
undefined params; B: multi-pattern `_kcat_forward`; C: failure-row mechanism
recording). Those keep the analytic `_kcat_forward` for the relative
`scale_k_to_kcat = 1.0` default; this enhancement adds a measured-concentration
path for absolute kcat values.
