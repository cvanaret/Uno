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
#   G. Li,
#   "The secant/finite difference algorithm for solving sparse
#   nonlinear systems of equations",
#   SIAM Journal on Optimization, (to appear), 1990.

#   SIF input: Ph. Toint, Aug 1990.

#   classification SUR2-AN-V-0

param N:=10000;
var x{1..N} := 4.0;

minimize f:
	sum {i in 1..N} 4*(-x[1]+x[i]^2)^2 + sum {i in 1..N} (x[i]-1.0)^2;

solve; display f; display x;
