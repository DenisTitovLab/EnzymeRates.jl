#!/bin/bash
# ABOUTME: Package-side skeptic runs under one hold of the shared Julia lock: the adversarial
# ABOUTME: export, then the probe's own pkg_run.jl (copied) to confirm its package numbers.
E=/tmp/claude-501/-home-denis-linux--julia-dev-EnzymeRates/12f42995-f5ba-4fda-9510-dd9b87f7d069/scratchpad/er_oracle
t0=$(date +%s)
julia --startup-file=no --project=$E --heap-size-hint=2500M /tmp/claude-501/-home-denis-linux--julia-dev-EnzymeRates/12f42995-f5ba-4fda-9510-dd9b87f7d069/scratchpad/proto/skeptic_P2-solve/adv/export_adv.jl > /tmp/claude-501/-home-denis-linux--julia-dev-EnzymeRates/12f42995-f5ba-4fda-9510-dd9b87f7d069/scratchpad/proto/skeptic_P2-solve/logs/export_adv.log 2>&1
echo "export_adv exit $? $(($(date +%s)-t0)) s"
t0=$(date +%s)
julia --startup-file=no --project=$E --heap-size-hint=2500M /tmp/claude-501/-home-denis-linux--julia-dev-EnzymeRates/12f42995-f5ba-4fda-9510-dd9b87f7d069/scratchpad/proto/skeptic_P2-solve/rerun/pkg_run.jl > /tmp/claude-501/-home-denis-linux--julia-dev-EnzymeRates/12f42995-f5ba-4fda-9510-dd9b87f7d069/scratchpad/proto/skeptic_P2-solve/logs/rerun_pkg_run.log 2>&1
echo "pkg_run exit $? $(($(date +%s)-t0)) s"
