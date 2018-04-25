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

#   Source: problem 4 starting point D in
#   M.H. Wright:
#   "Numerical Methods for Nonlinearly Constrained Optimization",
#   Ph.D. Thesis, Stanford (USA), 1976.

#   SIF input: Ph. Toint, March 1991.
#              correction by Ph. Shott, January, 1995.

#   classification OQR2-AN-5-3

param x_init{1..5};
var x{i in 1..5} := x_init[i];

minimize f:
	x[1]^2
	+ (x[1]-x[2])^2
	+ (x[2]-x[3])^3
	+ (x[3]-x[4])^4
	+ (x[4]-x[5])^4
;
subject to cons1:
	x[2]^2+x[3]^2+x[1]-3*sqrt(2)-2 = 0;
subject to cons2:
	-x[3]^2+x[2]+x[4]-2*sqrt(2)+2 = 0;
subject to cons3:
	x[1]*x[5]-2 = 0;

data;
param x_init:=
1	-1
2	2
3	1
4	-2
5	-2;

solve;
display f;
display x;
