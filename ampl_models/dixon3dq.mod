# AMPL Model by Hande Y. Benson
#
# Copyright (C) 2001 Princeton University
# All Rights Reserved
#
# Permission to use, copy, modify, and distribute this software and
# its documentation for any purpose and without fee is hereby
# granted, provided that the above copyright notice appear in all
# copies and that the copyright notice and this
# permission notice appear in all supporting documentation.                     

#   Source: problem 156 (p. 51) in
#   A.R. Buckley,
#   "Test functions for unconstrained minimization",
#   TR 1989CS-3, Mathematics, statistics and computing centre,
#   Dalhousie University, Halifax (CDN), 1989.

#   SIF input: Ph. Toint, Dec 1989.

#   classification QUR2-AN-V-0

param n := 10;

var x {1..n} := -1.0;

# rvdb comment: the sum should start at 1.
minimize obj: (x[1]-1.0)^2 + sum {j in 2..n-1} (x[j]-x[j+1])^2 + (x[n]-1.0)^2;

solve;
display x;
display obj;
