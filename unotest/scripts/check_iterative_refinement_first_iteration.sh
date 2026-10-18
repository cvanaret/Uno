#!/usr/bin/env bash
# Checks that, on the given instance, iterative refinement fails at the very first iteration and the pivoting tolerance
# of the linear solver is increased before the first iterate is accepted.
# Usage: check_iterative_refinement_first_iteration.sh <uno_ampl> <instance.nl> <linear_solver>
# The check relies on the following DEBUG log messages (keep in sync with the code):
#  - "### Outer iteration <k>"                          (Uno.cpp)
#  - "Iterative refinement failed"                      (EQPSolver.cpp)
#  - "Increasing the <solver> pivoting tolerance"       (MA27Solver.cpp, MA57Solver.cpp, MUMPSSolver.cpp)
set -uo pipefail
if [ "$#" -ne 3 ]; then
   echo "Usage: $0 <uno_ampl> <instance.nl> <linear_solver>"; exit 2
fi
UNO_AMPL="$1"; INSTANCE="$2"; LINEAR_SOLVER="$3"
LOG="$(mktemp)"
trap 'rm -f "$LOG"' EXIT

fail() {
   echo "FAIL: $1"; tail -30 "$LOG"; exit 1
}

"$UNO_AMPL" "$INSTANCE" preset=ipopt linear_solver="$LINEAR_SOLVER" max_iterations=2 logger=DEBUG > "$LOG" 2>&1
status=$?
# uno_ampl returns a nonzero status only upon an exception; >= 128 means it was killed by a signal (crash, sanitizer abort)
if [ "$status" -ne 0 ]; then
   fail "uno_ampl exited with status ${status}"
fi
if grep -q -e "ERROR: AddressSanitizer" -e "runtime error:" "$LOG"; then
   fail "sanitizer error"
fi

# window of the first outer iteration in the DEBUG log
first_iteration_line=$(grep -n -m1 '^### Outer iteration 1$' "$LOG" | cut -d: -f1)
second_iteration_line=$(grep -n -m1 '^### Outer iteration 2$' "$LOG" | cut -d: -f1)
second_iteration_line=${second_iteration_line:-$(wc -l < "$LOG")}
if [ -z "$first_iteration_line" ]; then
   fail "the first outer iteration was not found in the log"
fi
failure_line=$(grep -n -m1 "Iterative refinement failed" "$LOG" | cut -d: -f1)
increase_line=$(grep -n -m1 "Increasing the ${LINEAR_SOLVER} pivoting tolerance" "$LOG" | cut -d: -f1)
if [ -z "$failure_line" ]; then
   fail "iterative refinement never failed"
fi
if [ -z "$increase_line" ]; then
   fail "the ${LINEAR_SOLVER} pivoting tolerance was not increased"
fi
if [ "$increase_line" -lt "$first_iteration_line" ] || [ "$increase_line" -gt "$second_iteration_line" ]; then
   fail "the pivoting tolerance was not increased during the first iteration"
fi
if [ "$failure_line" -gt "$increase_line" ]; then
   fail "the pivoting tolerance was increased before iterative refinement failed"
fi
sed -n "${failure_line}p;${increase_line}p" "$LOG"
echo "PASS"
