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

#   Source: problem 5 in
#   P.T. Boggs and J.W. Tolle,
#   "A strategy for global convergence in a sequential 
#    quadratic programming algorithm",
#   SINUM 26(3), pp. 600-623, 1989.

#   The problem as stated in the paper seems to contain a typo.
#   The sign of the x3 squared term in the first constraint has been
#   set to + instead of - in order to ensdure that the problem is 
#   bounded below and the optimal point stated recovered.

#   The problem is not convex.

#   SIF input: Ph. Toint, June 1993.

#   classification QQR2-AN-3-2

param N;
param x_init{1..N};
var x{i in 1..N} := x_init[i];

minimize f:
	1000 - x[1]^2 - x[3]^2 - 2*x[2]^2 - x[1]*x[2] - x[1]*x[3];

subject to cons1:
	-25 + x[1]^2 + x[2]^2 + x[3]^2 = 0;

subject to cons2:
	8*x[1]+14*x[2]+7*x[3] - 56 = 0;

data;
param N:= 3;
param x_init:=
1	2.0
2	2.0
3	2.0;

solve; display x; 
