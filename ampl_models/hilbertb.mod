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

#   Source: problem 19 (p. 59) in
#   A.R. Buckley,
#   "Test functions for unconstrained minimization",
#   TR 1989CS-3, Mathematics, statistics and computing centre,
#   Dalhousie University, Halifax (CDN), 1989.

#   SIF input: Ph. Toint, Dec 1989.

#   classification QUR2-AN-V-0

param N:=50;
param D:=5.0;
var x{1..N} := -3.0;

minimize f:
	sum {i in 1..N} (sum {j in 1..i-1} x[i]*x[j]/(i+j-1) + (x[i]^2)*(D+1/(4*i-2)));

solve; display f; display x;
