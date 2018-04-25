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

#   Source: Problem 121 (p. 99) in
#   A.R. Buckley,
#   "Test functions for unconstrained minimization",
#   TR 1989CS-3, Mathematics, statistics and computing centre,
#   Dalhousie University, Halifax (CDN), 1989.

#   SIF input: Ph. Toint, Dec 1989.

#   classification OXR2-MN-V-0

param N := 1001;
param lambda := -3.4;
param h := 1/(N+1);

var x{i in 0..N+1} := if (i==0 || i==N+1) then 0
		      else -0.1*h*i^2;

minimize f:
	2*lambda*h*(exp(x[1])-exp(x[0]))/(x[1]-x[0])
	+ sum {i in 1..N} 2*x[i]^2/h
	- sum {i in 1..N} 2*x[i]*x[i-1]/h
	+ sum {i in 1..N} 2*lambda*h*(exp(x[i+1])-exp(x[i]))/(x[i+1]-x[i])
;

fix x[0] := 0;
fix x[N+1] := 0;

#option loqo_options "verbose=2 timing=1 convex";	
solve;
display f;
display x;
