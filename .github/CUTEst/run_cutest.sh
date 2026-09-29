#!/usr/bin/env bash
# Usage: run_cutest.sh <label> <problem>

set -u
label=$1; problem=$2
export TERM=xterm GFORTRAN_UNBUFFERED_PRECONNECTED=y # unbuffer the Fortran driver's stdout/stderr
# 2>&1 so stderr (std::terminate, gfortran backtrace) is captured; tee prints everything live
stdbuf -oL -eL "$CUTEST/bin/runcutest" -A "$MYARCH" -p uno -D "$problem" -o 1 2>&1 | tee "$label.log"
status=${PIPESTATUS[0]}

if [ "$status" -ne 0 ]; then
  echo "::error::[$label] runcutest exited with status $status"; exit 1
fi
# runcutest swallows the solver's exit code, so detect crashes from the output
if grep -Eq 'terminate called|Program received signal|Segmentation fault|Backtrace for this error|could not be found' "$label.log"; then
  echo "::error::[$label] Uno crashed on $problem"; exit 1
fi