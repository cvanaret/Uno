#!/usr/bin/env bash
# Checks that, on the given instance, iterative refinement fails at the very first iteration and that the quality of the
# factorization is improved before the first iterate is accepted:
#  - mode "pivoting" (linear_system_scaling = none): the pivoting tolerance of the linear solver is increased;
#  - mode "scaling" (linear_system_scaling = ruiz, on demand): the linear system scaling is switched on and the matrix is
#    refactorized, before any increase of the pivoting tolerance.
# Usage: check_iterative_refinement_first_iteration.sh <uno_ampl> <instance.nl> <linear_solver> [pivoting|scaling] [options...]
# The check relies on the following DEBUG log messages (keep in sync with the code):
#  - "### Outer iteration <k>"                          (Uno.cpp)
#  - "Iterative refinement failed"                      (EQPSolver.cpp)
#  - "Refactorizing the augmented matrix with scaling"  (EQPSolver.cpp)
#  - "Increasing the <solver> pivoting tolerance"       (MA27Solver.cpp, MA57Solver.cpp, MUMPSSolver.cpp)
set -uo pipefail
if [ "$#" -lt 3 ]; then
   echo "Usage: $0 <uno_ampl> <instance.nl> <linear_solver> [pivoting|scaling] [options...]"; exit 2
fi
UNO_AMPL="$1"; INSTANCE="$2"; LINEAR_SOLVER="$3"; MODE="${4:-pivoting}"
shift $(( $# < 4 ? $# : 4 ))
case "$MODE" in
   pivoting) SCALING_OPTIONS=(linear_system_scaling=none) ;;
   scaling) SCALING_OPTIONS=(linear_system_scaling=ruiz linear_system_scaling_on_demand=true) ;;
   *) echo "Unknown mode ${MODE} (pivoting|scaling)"; exit 2 ;;
esac
LOG="$(mktemp)"
trap 'rm -f "$LOG"' EXIT

fail() {
   echo "FAIL: $1"; tail -30 "$LOG"; exit 1
}

"$UNO_AMPL" "$INSTANCE" preset=ipopt linear_solver="$LINEAR_SOLVER" max_iterations=2 logger=DEBUG "${SCALING_OPTIONS[@]}" "$@" \
   > "$LOG" 2>&1
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
in_first_iteration() {
   [ "$1" -ge "$first_iteration_line" ] && [ "$1" -le "$second_iteration_line" ]
}
failure_line=$(grep -n -m1 "Iterative refinement failed" "$LOG" | cut -d: -f1)
scaling_line=$(grep -n -m1 "Refactorizing the augmented matrix with scaling" "$LOG" | cut -d: -f1)
increase_line=$(grep -n -m1 "Increasing the ${LINEAR_SOLVER} pivoting tolerance" "$LOG" | cut -d: -f1)
if [ -z "$failure_line" ]; then
   fail "iterative refinement never failed"
fi

if [ "$MODE" = "pivoting" ]; then
   if [ -n "$scaling_line" ]; then
      fail "the linear systems were scaled with linear_system_scaling = none"
   fi
   if [ -z "$increase_line" ]; then
      fail "the ${LINEAR_SOLVER} pivoting tolerance was not increased"
   fi
   if ! in_first_iteration "$increase_line"; then
      fail "the pivoting tolerance was not increased during the first iteration"
   fi
   if [ "$failure_line" -gt "$increase_line" ]; then
      fail "the pivoting tolerance was increased before iterative refinement failed"
   fi
   sed -n "${failure_line}p;${increase_line}p" "$LOG"
else
   if [ -z "$scaling_line" ]; then
      fail "the linear system scaling was not switched on"
   fi
   if ! in_first_iteration "$scaling_line"; then
      fail "the linear system scaling was not switched on during the first iteration"
   fi
   if [ "$failure_line" -gt "$scaling_line" ]; then
      fail "the linear system scaling was switched on before iterative refinement failed"
   fi
   # the scaling comes first; the pivoting tolerance may be increased afterwards (if refinement fails again)
   if [ -n "$increase_line" ] && [ "$increase_line" -lt "$scaling_line" ]; then
      fail "the pivoting tolerance was increased before the linear system scaling was switched on"
   fi
   sed -n "${failure_line}p;${scaling_line}p" "$LOG"
fi
echo "PASS"
