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

#   Source: a variant of problem 4 in
#   P.T. Boggs and J.W. Tolle,
#   "A strategy for global convergence in a sequential 
#    quadratic programming algorithm",
#   SINUM 26(3), pp. 600-623, 1989.

#   The original problem seems to be unbounded.  The contribution of
#   x3 in the first constraint has been squared instead of cubed.

#   The problem is not convex.
 
#   SIF input: Ph. Toint, June 1993.

#   classification QQR2-AN-3-2

param xinit{1..3};
var x{i in 1..3} := xinit[i];

minimize f:
	x[1]-x[2]+x[2]^3;
subject to cons1:
	-25+x[1]^2+x[2]^2+x[3]^2 = 0;
subject to cons2:
	x[1]+x[2]+x[3]-1 = 0;

data;
param xinit:=
1	4.0382
2	-2.9470
3	-0.09115;

solve; display f; display x;
