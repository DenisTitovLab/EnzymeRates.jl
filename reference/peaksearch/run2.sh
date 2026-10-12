#!/bin/bash
# ABOUTME: Run 2 of the peak search: full adversarial search, focused parallel-route search,
# ABOUTME: targeted guesses and census, profiles of the best points, RE self-check.
D=/tmp/claude-501/-home-denis-linux--julia-dev-EnzymeRates/12f42995-f5ba-4fda-9510-dd9b87f7d069/scratchpad
L=$D/proto/julia.lock
P=$D/peaksearch
J="julia --project=$P --heap-size-hint=1500M"
echo "full start $(date +%T)"
timeout -s KILL 5400 flock -w 14400 $L $J $P/search.jl full 2000000 > $P/full.log 2>&1
echo "full exit $? $(date +%T)"
timeout -s KILL 3600 flock -w 14400 $L $J $P/search.jl d1focus 5000000 d1,d2 zero -6,-3,0,3,6 > $P/d1focus.log 2>&1
echo "d1focus exit $? $(date +%T)"
timeout -s KILL 3600 flock -w 14400 $L $J $P/guesses.jl > $P/guesses.log 2>&1
echo "guesses exit $? $(date +%T)"
timeout -s KILL 3600 flock -w 14400 $L $J $P/profiles.jl results_full.jls > $P/profiles.log 2>&1
echo "profiles exit $? $(date +%T)"
timeout -s KILL 3600 flock -w 14400 $L $J $P/selfcheck_re.jl > $P/selfcheck_re.log 2>&1
echo "selfcheck_re exit $? $(date +%T)"
