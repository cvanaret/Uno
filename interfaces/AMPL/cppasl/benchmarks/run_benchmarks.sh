#!/bin/sh
# Generates the Mittelmann qcqp instances and compares cppasl with ASL2.
# usage: benchmarks/run_benchmarks.sh <build directory> <scratch directory>
set -e
BUILD=${1:-build}; OUT=${2:-/tmp/cppasl_instances}; mkdir -p "$OUT"
G="$BUILD/generate_qcqp_nl"; C="$BUILD/compare_with_asl2"
#      name        n    ml  mq  pl    pq sd sq   sp   plf  pqf  seed
while read name n ml mq pl pq sd sq sp plf pqf seed; do
   [ -f "$OUT/$name.nl" ] || "$G" $n $ml $mq $pl $pq $sd $sq $sp $plf $pqf $seed text "$OUT/$name.nl"
   "$C" "$OUT/$name.nl" 3
done <<INSTANCES
qcqp500_A     500  10  0  100   10 1 .01 .2  .1  .2 1
qcqp1500_B   1500 500  0  10000  8 1 .01 .01 .01 .2 2
qcqp750_nc    750  10  0  100   10 0 .01 .2  .1  .2 3
INSTANCES
