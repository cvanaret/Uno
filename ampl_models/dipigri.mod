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
#   G. Di Pillo and L. Grippo,
#   "An new augmented Lagrangian function for inequality constraints
#   in nonlinear programming problems",
#   JOTA, vol. 36, pp. 495-519, 1982.

#   SIF input: Ph. Toint, June 1990.

#   classification OOR2-AN-7-4

param x_init{1..7};
var x{i in 1..7} := x_init[i];

minimize f:
	(x[1]-10)^2+5*(x[2]-12)^2+x[3]^4+3*(x[4]-11)^2+
	10*x[5]^6 + 7*x[6]^2+x[7]^4-4*x[6]*x[7] -10*x[6]-8*x[7]
;
subject to cons1:
	2*x[1]^2+3*x[2]^4+4*x[4]^2+x[3]+5*x[5]-127.0 <= 0;
subject to cons2:
	10*x[3]^2+7*x[1]+3*x[2]+x[4]-x[5]-282.0 <= 0;
subject to cons3:
	x[2]^2+6*x[6]^2+23*x[1]-8*x[7]-196.0 <= 0;
subject to cons4:
	4*x[1]^2+x[2]^2-3*x[1]*x[2]+2*x[3]^2+5*x[6]-11*x[7] <=
	0;

data;
param x_init:=
1	1
2	2
3	0
4	4
5	0
6	1
7	1;

solve;
display f;
display x;
