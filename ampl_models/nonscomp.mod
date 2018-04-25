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
#   M. Lescrenier,
#   "Towards the use of supercomputers for large scale nonlinear
#   partially separable optimization",
#   PhD Thesis, FUNDP (Namur, B), 1989.

#   SIF input: Ph. Toint, May 1990.

#   classification SBR2-AN-V-0

param N := 10000;
var x{i in 1..N} := 3, >= if (i mod 3 == 0) then 1 else -100, <= 100;

minimize f:
	(x[1] - 1)^2
	+ sum {i in 2..N} 4*(x[i]-x[i-1]^2)^2
;

solve;
display f;
display x;
