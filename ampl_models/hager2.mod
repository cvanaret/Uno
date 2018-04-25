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

#   Source: problem P2 in
#   W.W. Hager,
#   "Multiplier Methods for Nonlinear Optimal Control",
#   SIAM J. on Numercal Analysis 27(4): 1061-1080, 1990.

#   SIF input: Ph. Toint, March 1991.

#   classification OLR2-AN-V-V

param n := 5000;

param h := 1/n;

var x{0..n} := 0.0;
var u{1..n} := 0.0;

minimize f:
	sum {i in 1..n} h*(x[i-1]^2 + x[i-1]*x[i] + x[i]^2)/6 +
	sum {i in 1..n} h*(u[i])^2/4 ;

subject to cons1{i in 1..n}:
	(n-0.25)*x[i] - (n+0.25)*x[i-1] - u[i] = 0;

fix x[0] := 1.0;
solve; 

