#!/usr/bin/env bash
# Usage: benchmark_instance.sh <problem> <SIF file (.SIF or .SIF.bz2)>
#
# Solves one problem in SIF format with Uno through `runcutest -p uno`, in its
# own directory under $WORK_DIR, and writes:
#   $LOG_DIR/<problem>.log            full output (SIF decoding, compilation, Uno)
#   $ROWS_DIR/<problem>.tsv           <problem> \t solved|failed \t <seconds>
#
# The time is Uno's own "CPU time (secs)" when Uno reports it (it excludes the
# SIF decoding and compilation), otherwise the wall-clock time of the whole run.
#
# Environment: CUTEST, MYARCH, OPTION_FILE (absolute path), TIMELIMIT,
#              WORK_DIR, LOG_DIR, ROWS_DIR

set -u
problem=$1
sif_file=$2
dir="$WORK_DIR/$problem"
log="$LOG_DIR/$problem.log"
rm -rf "$dir"
mkdir -p "$dir"
cd "$dir" || exit 1

# the CUTEst driver reads ./uno.opt; options are applied in order, so the
# appended time limit overrides any time_limit of the option file
cp "$OPTION_FILE" uno.opt
printf '\ntime_limit %s\n' "$TIMELIMIT" >> uno.opt

# sifdecoder looks for <problem>.SIF in the current directory first: copying it
# here makes the problem's collection explicit (some names exist in several
# collections, with different contents). Very large problems may be bzip2-compressed.
case "$sif_file" in
  *.bz2) bunzip2 -kc "$sif_file" > "$problem.SIF" ;;
  *)     cp "$sif_file" "$problem.SIF" ;;
esac

echo "[$(date '+%Y-%m-%d %H:%M:%S')] Solving ${problem}"
export TERM=xterm GFORTRAN_UNBUFFERED_PRECONNECTED=y
start=$(date +%s.%N)
# Uno enforces time_limit itself; the outer `timeout` also covers SIF decoding
# and compilation, which can take a while for the largest problems
timeout -k 30 $((TIMELIMIT + 600)) "$CUTEST/bin/runcutest" -A "$MYARCH" -p uno -D "$problem" > "$log" 2>&1
status=$?
end=$(date +%s.%N)

secs=$(awk '/CPU time \(secs\)/ {value = $NF} END {if (value != "") printf "%.2f", value}' "$log")
if [ -z "$secs" ]; then
  secs=$(awk "BEGIN{printf \"%.2f\", ${end} - ${start}}")
fi

# runcutest swallows the solver's exit code, so a crash is detected from the
# output; "Success" is Uno's optimization status, as in the Mittelmann benchmark
if [ "$status" -eq 0 ] \
   && grep -Eq '^[[:space:]]*Optimization status:[[:space:]]+Success' "$log" \
   && ! grep -Eq 'terminate called|Program received signal|Segmentation fault|Backtrace for this error' "$log"; then
  result=solved
else
  result=failed
fi
printf '%s\t%s\t%s\n' "$problem" "$result" "$secs" > "$ROWS_DIR/$problem.tsv"
echo "${problem}: ${result} in ${secs}s"

# the decoded sources and the executable take a lot of space over hundreds of problems
cd "$WORK_DIR" && rm -rf "$dir"
exit 0
