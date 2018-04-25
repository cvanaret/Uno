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

#   Source: problem P1 in
#   W.W. Hager,
#   "Multiplier Methods for Nonlinear Optimal Control",
#   SIAM J. on Numerical Analysis 27(4): 1061-1080, 1990.

#   SIF input: Ph. Toint, March 1991.

#   classification SLR2-AN-V-V

param N:=5000;
var x{0..N} := 0.0;
var u{1..N} := 0.0;

minimize f:
	0.5*x[N]^2 + sum {i in 1..N} (u[i]^2)/(2*N);
subject to cons1{i in 1..N}:
	(N-0.5)*x[i] + (-N-0.5)*x[i-1] - u[i] = 0;
subject to cons2:
	x[0] = 1.0;

solve; display f; display x;
