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
#   P.T. Boggs and J.W. Tolle,
#   "A strategy for global convergence in a sequential 
#    quadratic programming algorithm",
#   SINUM 26(3), pp. 600-623, 1989.

#   The problem is not convex.

#   SIF input: Ph. Toint, June 1993.

#   classification OQR2-AN-5-3

param x_init{1..5};
var x{i in 1..5} := x_init[i];

minimize f:
	100*(x[2]-x[1]^2)^2 + (x[1]-1.0)^2;
subject to cons1:
	x[1]*x[2] - x[3]^2 = 1.0;
subject to cons2:
	x[2]^2 - x[4]^2 + x[1] = 0.0;
subject to cons3:
	x[5]^2 + x[1] = 0.5;

data;
param x_init:=
1	-2.0
2	1.0
3	1.0
4	1.0
5	1.0;

solve;
display f;
display x;
