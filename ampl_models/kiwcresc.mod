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
#   K.C. Kiwiel,
#   "Methods of Descent for Nondifferentiable Optimization"
#   Lectures Notes in Mathematics 1133, Springer Verlag, 1985.

#   SIF input: Ph. Toint, Nov 1993.

#   classification  LQR2-AN-3-2

param xinit{1..2};
var x{i in 1..2} := xinit[i];
var u;
minimize f:
	u;
subject to cons1:
	-u+x[2]-1+x[1]^2+(x[2]-1)^2 <= 0;
subject to cons2:
	-u+x[2]+1-x[1]^2-(x[2]-1)^2 <= 0;

data;
param xinit:= 1 -1.5 2 2.0;

solve; display f; display x;
