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

#   classification QOR2-AN-V-V

param n := 5000;
param h := 1/n;

var x{1..n-1};
var y{1..n,1..2};

minimize f:
	5*h*(0.5*y[1,1]^2+0.5*y[1,2]^2+x[1]^2 + sum {t in 2..n-1} (y[t,1]^2 + y[t,2]^2 + x[t]^2) + 0.5*y[n,1]^2 + 0.5*y[n,2]^2);
subject to cons1{t in 1..n-1}:
	-y[t+1,1] + (1+5*h)*y[t,1] - 5*h*y[t,2] + 5*h*x[t] - 5*h*y[t,1]*y[t,2]^2= 0;
subject to cons2{t in 1..n-1}:
	y[t,2] - y[t+1,2] + 5*h*y[t,1] = 0;
fix y[1,1] := 0.0;
fix y[1,2] := 1.0;

solve;
