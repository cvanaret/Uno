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

#   Source: an example problem (p. 83) in
#   J.E. Dennis and R.B. Schnabel,
#   "Numerical Methods for Unconstrained Optimization and Nonlinear
#   Equations",
#   Prentice-Hall, Englewood Cliffs, 1983.

#   SIF input: Ph. Toint, Nov 1990.

#   classification SUR2-AN-3-0

var x{1..3} := 10.0;

minimize f:
	(x[1]^2+x[2]^3-x[3]^4)^2
	+(2*x[1]*x[2]*x[3])^2
	+(2*x[1]*x[2]-3*x[2]*x[3]+x[1]*x[3])^2
;

solve;
display f;
display x;
