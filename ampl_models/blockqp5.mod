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

param n := 1000;
param b := 5;
var x{1..n} >= -1, <= 1, := 0.99;
var y{1..n} >= -1, <= 1, := -0.99;
var z{1..b} >= 0, <= 1, := 0.5;

minimize f:
	sum {i in 1..n} (i/n)*(x[i]*y[i]) + sum {j in 1..b} 0.5*z[j]^2
;

subject to cons1:
	sum {i in 1..n} (x[i] + y[i]) + sum {j in 1..b} z[j] >= b+1;
subject to cons2{i in 1..n}:
	x[i] + y[i] + sum {j in 1..b} z[j] = b;

solve;
