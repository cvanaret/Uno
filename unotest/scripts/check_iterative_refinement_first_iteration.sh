#!/usr/bin/env bash
# Checks that, on the given instance, iterative refinement fails at the very first iteration and the pivoting tolerance
# of the linear solver is increased before the first iterate is accepted.
# Usage: check_iterative_refinement_first_iteration.sh <uno_ampl> <instance.nl> <linear_solver>
set -uo pipefail
UNO_AMPL="$1"; INSTANCE="$2"; LINEAR_SOLVER="$3"
LOG="$(mktemp)"
trap 'rm -f "$LOG"' EXIT
"$UNO_AMPL" "$INSTANCE" preset=ipopt linear_solver="$LINEAR_SOLVER" max_iterations=2 logger=DEBUG > "$LOG" 2>&1

increase_line=$(grep -n -m1 "Increasing the ${LINEAR_SOLVER} pivoting tolerance" "$LOG" | cut -d: -f1)
# window of the first outer iteration in the DEBUG log
first_iteration_line=$(grep -n -m1 '^### Outer iteration 1$' "$LOG" | cut -d: -f1)
second_iteration_line=$(grep -n -m1 '^### Outer iteration 2$' "$LOG" | cut -d: -f1)
second_iteration_line=${second_iteration_line:-$(wc -l < "$LOG")}
if [ -z "$increase_line" ]; then
   echo "FAIL: the ${LINEAR_SOLVER} pivoting tolerance was not increased"; tail -30 "$LOG"; exit 1
fi
if [ -z "$first_iteration_line" ] || [ "$increase_line" -lt "$first_iteration_line" ] || [ "$increase_line" -gt "$second_iteration_line" ]; then
   echo "FAIL: the pivoting tolerance was not increased during the first iteration"; tail -30 "$LOG"; exit 1
fi
grep -m1 "Iterative refinement failed" "$LOG"
grep -m1 "Increasing the ${LINEAR_SOLVER} pivoting tolerance" "$LOG"
echo "PASS"
