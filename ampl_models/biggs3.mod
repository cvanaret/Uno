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

#   Source: Problem 152 in
#   A.R. Buckley,
#   "Test functions for unconstrained minimization",
#   TR 1989CS-3, Mathematics, statistics and computing centre,
#   Dalhousie University, Halifax (CDN), 1989.

#   SIF input: Ph. Toint, Dec 1989.

#   classification SXR2-AN-6-0

param N:=6;
param M:=13;
param xinit{1..N};
var x{i in 1..N}:=xinit[i];

minimize f:
	sum {i in 1..M} (-exp(-0.1*i)+5*exp(-i)-3*exp(-0.4*i)
	+ x[3]*exp(-0.1*i*x[1]) - x[4]*exp(-0.1*i*x[2]) +
	x[6]*exp(-0.1*i*x[5]))^2;
subject to cons1:
	x[3] = 1;
subject to cons2:
	x[5] = 4;
subject to cons3:
	x[6] = 3;
data;
param xinit:=
1	1
2	2
3	1
4	1
5	4
6	3;

solve;
display f;
display x;

