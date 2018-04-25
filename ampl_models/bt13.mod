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

#   Source: problem 13 in
#   P.T. Boggs and J.W. Tolle,
#   "A strategy for global convergence in a sequential 
#    quadratic programming algorithm",
#   SINUM 26(3), pp. 600-623, 1989.

#   The problem has been modified by adding a lower bound of 0.0 on x5
#   in order to make it bounded below.

#   The problem is not convex.

#   SIF input: Ph. Toint, June 1993.

#   classification LQR2-AY-5-1

param N;
param x_init{1..N};
var x{i in 1..N}:= x_init[i];

minimize f:
	x[5];

subject to cons1:
	x[1]^2+(x[1]-2*x[2])^2 + (x[2]-3*x[3])^2+(x[3]-4*x[4])^2 
	- x[5]^2= 0;

subject to cons2:
	x[5] >= 0;

data;

param N := 5;
param x_init:=
1	1.0
2	2.0
3	3.0
4	3.0
5	228.0;

solve;
display x;
