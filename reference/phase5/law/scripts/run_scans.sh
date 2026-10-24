#!/bin/bash
# ABOUTME: Runs the R5.1 remark scans one Julia at a time through the research lock: the upstream
# ABOUTME: reproducer as a positive control, then scan.jl per variant, on Julia 1.12 and 1.13.
S=/tmp/claude-501/-home-denis-linux--julia-dev-EnzymeRates/12f42995-f5ba-4fda-9510-dd9b87f7d069/scratchpad
R=$S/phase45/research5/law
J=$S/phase45/jl_research.sh
export JULIA_LLVM_ARGS="-pass-remarks=julia-licm -pass-remarks-analysis=julia-licm"
for ver in ${VERSIONS:-112 113}; do
    ch=""; [ "$ver" = 113 ] && ch="+1.13"
    [ -n "$CONTROL" ] && $J $R/logs/remarks_repro_$ver.log 300 $ch --startup-file=no --pkgimages=no $S/phase3/gcprobe/repro_upstream.jl
    for v in ${VARIANTS:-U A B C}; do
        VARIANT=$v $J $R/logs/remarks_scan_${v}_$ver.log 600 $ch --startup-file=no --pkgimages=no --heap-size-hint=1500M $R/scripts/scan.jl
    done
done
echo "all scans done" > $R/logs/run_scans.done
