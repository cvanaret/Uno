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

#   Source: problem 7 in
#   J.J. More', B.S. Garbow and K.E. Hillstrom,
#   "Testing Unconstrained Optimization Software",
#   ACM Transactions on Mathematical Software, vol. 7(1), pp. 17-41, 1981.

#   See also Buckley#12 (p. 58)
#   SIF input: Ph. Toint, Dec 1989.

#   classification  SUR2-AN-3-0

param x_init{1..3};
var x{i in 1..3} := x_init[i];
var theta = if (x[1] > 0) then atan(x[2]/x[1])/(2*3.1415) else 
	    if (x[1] < 0) then atan(x[2]/x[1])/(2*3.1415) + 0.5 else 0.0;

minimize f:
	(10*(x[3]-10*theta))^2 + (10*(sqrt(x[1]^2+x[2]^2)-1))^2 + x[3]^2;

data;
param x_init := 1 -1 2 0 3 0;

solve;
display f;
display x;
