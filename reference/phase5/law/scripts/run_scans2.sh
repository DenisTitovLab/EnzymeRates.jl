#!/bin/bash
# ABOUTME: Runs the remaining R5.1 remark scans with fewer lock turns: scan_all.jl runs several
# ABOUTME: variants per Julia process (1.12: C, B, A; 1.13: C, B, A, U), plus the 1.13 reproducer.
S=/tmp/claude-501/-home-denis-linux--julia-dev-EnzymeRates/12f42995-f5ba-4fda-9510-dd9b87f7d069/scratchpad
L=$S/phase45/research5/law
J=$S/phase45/jl_research.sh
export JULIA_LLVM_ARGS="-pass-remarks=julia-licm -pass-remarks-analysis=julia-licm"
VARIANTS=C,B,A $J $L/logs/remarks_scanall_CBA_112.log 600 --startup-file=no --pkgimages=no --heap-size-hint=1500M $L/scripts/scan_all.jl
$J $L/logs/remarks_repro_113.log 300 +1.13 --startup-file=no --pkgimages=no $S/phase3/gcprobe/repro_upstream.jl
VARIANTS=C,B,A,U $J $L/logs/remarks_scanall_CBAU_113.log 600 +1.13 --startup-file=no --pkgimages=no --heap-size-hint=1500M $L/scripts/scan_all.jl
echo "scans2 done" > $L/logs/run_scans2.done
