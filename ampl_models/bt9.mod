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

#   Source: problem 9 in
#   P.T. Boggs and J.W. Tolle,
#   "A strategy for global convergence in a sequential 
#    quadratic programming algorithm",
#   SINUM 26(3), pp. 600-623, 1989.


#   The problem as stated in the paper seems to contain a typo.
#   In order to make the problem bounded below and the second constraint
#   feasible at the proposed solution, the sign of x2 in the second constraint
#   has been set to - instead of +.

#   The problem is not convex.

#   SIF input: Ph. Toint, June 1993.

#   classification LOR2-AN-4-2

var x{1..4} := 2.0;

minimize f:
	-x[1];
subject to cons1:
	x[2]-x[1]^3-x[3]^2 = 0;
subject to cons2:
	-x[2]+x[1]^2-x[4]^2 = 0;

solve; display f; display x;
