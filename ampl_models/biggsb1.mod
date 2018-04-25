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
#   M. Batholomew-Biggs and F.G. Hernandez,
#   "Some improvements to the subroutine OPALQP for dealing with large
#    problems",
#   Numerical Optimization Centre, Hatfield, 1992.

#   SIF input: Ph Toint, April 1992.

#   classification QBR2-AN-V-V

param N := 1000;
var x{i in 1..N};

minimize f:
	(x[1]-1)^2 + sum {i in 1..N-1} (x[i+1]-x[i])^2 + (1-x[N])^2;

subject to cons1{i in 1..N-1}:
	0.0 <= x[i] <= 0.9;

solve;
