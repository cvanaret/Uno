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
#   Zhengua Lina and Yong Li,
#   "A Modified Scaling Newton-Type Method for Nonlinear Programming"
#   Department of Mathematics, Jilin University, Changchun, China, 1994.

#   SIF input: Ph. Toint, January 1994.

#   classification OQR2-AN-3-2

param xinit{1..3};
var x{i in 1..3} := xinit[i], >= 0;

minimize f:
	x[1]^3 - 6*x[1]^2 + 11*x[1] + x[2] + x[3];
subject to cons1:
	4 <= x[1]^2 + x[2]^2 + x[3]^2 <= 10;
subject to cons2:
	x[3] <= 5;

data;
param xinit:= 1 0.1 2 0.1 3 3.0;

solve; display f; display x;
