#!/usr/bin/env bash
# ABOUTME: Test-integrity gate enforced at every stage closeout.
# ABOUTME: Hard checks: no file deletion, no @testset count drop,
# ABOUTME: no @test_skip/@test_broken, no commented-out @test lines.
# ABOUTME: Soft check: WARNs about any modified @test/@test_throws line
# ABOUTME: for human review (weakened assertions can't be caught statically).
# Per spec §2 (NON-NEGOTIABLE). Run from repo root.
set -uo pipefail

BASE_REF="${1:-main}"
fail=0

# Hard Check 1: no test file deleted (renames OK; deletions forbidden).
deleted=$(git diff --name-status "$BASE_REF"..HEAD -- test/ 2>/dev/null \
    | awk '$1 == "D" { print $2 }')
if [ -n "$deleted" ]; then
    echo "FAIL [Check 1]: test file(s) deleted vs $BASE_REF:"
    echo "$deleted" | sed 's/^/    /'
    fail=1
fi

# Hard Check 2: @testset count never decreases — EXCEPT for testsets
# documented in docs/superpowers/refactor-deleted-tests.md per spec §2.1
# narrow exception (tests of deleted helpers; conditions are strict).
n_base=0
for f in $(git ls-tree -r "$BASE_REF" --name-only test/ 2>/dev/null | grep "\.jl$"); do
    c=$(git show "$BASE_REF":"$f" 2>/dev/null | grep -c "@testset" || true)
    n_base=$((n_base + c))
done
n_head=$(grep -rh "@testset" test/ 2>/dev/null | grep -c "@testset" || true)
# Each "### test_..." heading in the deleted-tests log = one permitted deletion.
n_documented=$(grep -c "^### test_" docs/superpowers/refactor-deleted-tests.md 2>/dev/null || echo 0)
n_base_adj=$((n_base - n_documented))
echo "@testset count: $BASE_REF=$n_base, documented deletions=$n_documented, adjusted base=$n_base_adj, HEAD=$n_head"
if [ "$n_head" -lt "$n_base_adj" ]; then
    echo "FAIL [Check 2]: @testset count decreased beyond documented deletions"
    echo "    Add an entry to docs/superpowers/refactor-deleted-tests.md per spec §2.1,"
    echo "    OR restore the deleted testset(s) and re-apply as mechanical adaptation."
    fail=1
fi

# Hard Check 3: no @test_skip or @test_broken added.
forbidden=$(git diff "$BASE_REF"..HEAD -- test/ | grep "^+" | grep -v "^+++" \
    | grep -E "@test_skip|@test_broken" || true)
if [ -n "$forbidden" ]; then
    echo "FAIL [Check 3]: @test_skip or @test_broken added vs $BASE_REF:"
    echo "$forbidden" | sed 's/^/    /'
    fail=1
fi

# Hard Check 4: no @test lines commented out.
commented=$(git diff "$BASE_REF"..HEAD -- test/ | grep "^+" | grep -v "^+++" \
    | grep -E "^\+\s*#+\s*@test\b" || true)
if [ -n "$commented" ]; then
    echo "FAIL [Check 4]: @test line(s) commented out vs $BASE_REF:"
    echo "$commented" | sed 's/^/    /'
    fail=1
fi

# Soft Check 5: WARN on any modified @test/@test_throws line.
# Spec §2 forbids changing hardcoded values / weakening assertions.
# Static detection is unreliable — we emit a WARN list for human review.
modified_tests=$(git diff "$BASE_REF"..HEAD -U0 -- test/ \
    | grep -E "^[+-]\s*@test(_throws)?\b" \
    | grep -v "^---" | grep -v "^+++" || true)
if [ -n "$modified_tests" ]; then
    echo ""
    echo "WARN [Check 5]: @test/@test_throws lines modified vs $BASE_REF."
    echo "  Human review REQUIRED — spec §2 forbids weakening assertions"
    echo "  (e.g., '== 3' → 'isa Int', tolerance relaxation, swapped operators):"
    echo "$modified_tests" | sed 's/^/    /'
    echo "  If any change is a weakening, REVERT and re-apply as mechanical only."
fi

if [ "$fail" -eq 0 ]; then
    echo "Test-integrity hard checks PASSED vs $BASE_REF."
fi
exit "$fail"
