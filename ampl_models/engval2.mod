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

#   Source: problem 15 in
#   A.R. Buckley,
#   "Test functions for unconstrained minimization",
#   TR 1989CS-3, Mathematics, statistics and computing centre,
#   Dalhousie University, Halifax (CDN), 1989.

#   SIF input: Ph. Toint, Dec 1989.

#   classification SUR2-AN-3-0

param x_init{1..3};
var x{i in 1..3} := x_init[i];

minimize f:
	(x[1]^2+x[2]^2+x[3]^2-1)^2
	+ (x[1]^2+x[2]^2+(x[3]-2)^2-1)^2
	+ (x[1]+x[2]+x[3]-1)^2
	+ (x[1]+x[2]-x[3]+1)^2
	+ (3*x[2]^2+x[1]^3+(5*x[3]-x[1]+1)^2-36)^2
;

data;
param x_init:= 1 1 2 2 3 0;

solve;
display f;
display x;
