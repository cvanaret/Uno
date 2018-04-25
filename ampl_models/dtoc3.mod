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

#   classification QLR2-AN-V-V

param n := 5000;
param s := 1/n;

var x{1..n-1};
var y{1..n,1..2};

minimize f:
	0.5*s*(sum{t in 1..n-1} (2*y[t+1,1]^2 + y[t+1,2]^2 + 6*x[t]^2));
subject to cons1{t in 1..n-1}:
	y[t,1] - y[t+1,1] + s*y[t,2] = 0;
subject to cons2{t in 1..n-1}:
	y[t,2] - y[t+1,2] - s*y[t,1] + s*x[t] = 0;

fix y[1,1] := 15.0;
fix y[1,2] := 5.0;
solve;
