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

#   classification OOR2-MN-V-V

param n := 5000;
param a := 0.07716;
param pi := 4*atan(1);

var x{0..n} <= 0.5, >= -0.5, := 0.25;
var u{0..n} <= 1.0, >= 0.0, := 0.25;
var us{i in 0..n} = u[i]/n;

minimize f:
	sum {i in 1..n} (-us[i]*(x[i]-cos(2*pi*i/n))^2/(2*n) - us[i-1]*(x[i-1]-cos(2*pi*(i-1)/n))^2/(2*n));
subject to cons1{i in 1..n}:
	x[i]*us[i]/(2*a) + x[i-1]*us[i-1]/(2*a) + x[i]*n - x[i-1]*n + us[i]*cos(2*pi*i/n)/(-2*a) + us[i-1]*cos(2*pi*(i-1)/n)/(-2*a) = 0;
subject to cons2:
	(x[0] - x[n])*n = 0;
fix x[0] := 0.25;
solve;

