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

#   Source: an example problem (p. 107) in
#   J.E. Dennis and R.B. Schnabel,
#   "Numerical Methods for Unconstrained Optimization and Nonlinear
#   Equations",
#   Prentice-Hall, Englewood Cliffs, 1983.

#   SIF input: Ph. Toint, Nov 1990.

#   classification SUR2-AY-2-0

var x1 := 2;
var x2 := 0;

minimize f:
	(2*(x1+x2)^2+(x1-x2)^2-8)^2	
	+ (5*x1^2+(x2-3)^2-9)^2
;

solve;
display f;
display x1, x2;
