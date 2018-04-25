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

#   classification LLR2-MN-V-V

param N:=5000;
param A:=0.07716;
param H:=1/N;
param pi := 3.1415;

var x1{0..N};
var x2{0..N} <= 0.125, >= -0.125;
var u{0..N} >= -1.0, <= 1.0;

minimize f:
	sum {i in 1..N} (-0.5*H*cos(2*pi*i*H)*x1[i] - 0.5*H*cos(2*pi*(i-1)*H)*x1[i] + H*(u[i]+u[i-1])/(8*pi^2));
subject to cons1{i in 1..N}:
	(x1[i]-x1[i-1])/H - 0.5*(x2[i]+x2[i-1]) = 0;
subject to cons2{i in 1..N}:
	(x2[i]-x2[i-1])/H - 0.5*(u[i]+u[i-1]) = 0;
subject to cons3:
	x1[0] = 0.0;
subject to cons4:
	x2[0] = 0.0;

solve; display f; display x1, x2, u;
