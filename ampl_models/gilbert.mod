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
#   J.Ch. Gilbert,
#   "On the Realization of the Wolfe Conditions in Reduced Quasi-Newton
#   Methods for Equality Constrained Optimization",
#   RR-2127, INRIA (F), 1993.

#   SIF input: Ph. Toint, April 1994

#   classification QQR2-AN-V-1

param n := 1000;
var x{i in 1..n} := (-1.0)^(i+1)*10.0;
minimize f:
	sum {i in 1..n} ((n+1-i)*x[i]/n-1.0)^2/2;
subject to cons1:
	(sum {i in 1..n} x[i]^2 - 1.0)/2 = 0.0;
solve;
display f;
display x;
