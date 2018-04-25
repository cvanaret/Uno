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

#   Source: problem 27 in
#   D.H. Himmelblau,
#   "Applied nonlinear programming",
#   McGraw-Hill, New-York, 1972.


#   See Buckley#77 (p. 62)

#   SIF input: Ph. Toint, Dec 1989.

#   classification OUR2-AN-2-0

param xinit{1..2};
var x{i in 1..2} := xinit[i];

minimize f:
	(x[1]*x[2]*(1-x[1])*(1-x[2]-x[1]*(1-x[1]^5)))^2;

data;
param xinit:= 1 -1.2 2 1.0;

solve; display f; display x;
