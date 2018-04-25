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

#   Source:
#   S. Lyle and N.K. Nichols,
#   "Numerical Methods for Optimal Control Problems with State Constraints",
#   Numerical Analysis Report 8/91, Dept of Mathematics, 
#   University of Reading, UK.

#   SIF input: Nick Gould, July 1991.

#   classification OOR2-MN-V-V

param n := 5000;
param a := 0.07716;
param pi := 4*atan(1);

var x{0..n} <= 0.5, >= -0.5, := 0.0;
var u{0..n} <= 1.0, >= 0.0, := 0.0;

minimize f:
	sum {i in 1..n} (-u[i]*(x[i]-cos(2*pi*i/n))^2/(2*n) - u[i-1]*(x[i-1]-cos(2*pi*(i-1)/n))^2/(2*n));
subject to cons1{i in 1..n}:
	x[i]*u[i]/(2*a) + x[i-1]*u[i-1]/(2*a) + x[i]*n - x[i-1]*n + u[i]*cos(2*pi*i/n)/(-2*a) + u[i-1]*cos(2*pi*(i-1)/n)/(-2*a) = 0;

fix x[0] := 0.25;
solve;

