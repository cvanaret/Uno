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
#   A.A. Brown and M. Bartholomew-Biggs,
#   "Some effective methods for unconstrained optimization based on
#   the solution of ordinary differential equations",
#   Technical Report 178, Numerical Optimization Centre, Hatfield
#   Polytechnic, (Hatfield, UK), 1987.

#   SIF input: Ph. Toint, June 1990.

#   classification OUR2-AN-2-0

param p := 10000;
var x1 := 0.86;
var x2 := 0.72;

minimize f:
	- 2*(x1-1)^2
	+ p*(-0.02 + (x2-x1^2)^2/p + (x1-1)^2)^2
;

solve;
display f;
display x1, x2;

