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

#   Source: problem 30 in
#   D.H. Himmelblau,
#   "Applied nonlinear programming",
#   McGraw-Hill, New-York, 1972.

#   See Buckley#88 (p. 65)

#   SIF input: Ph. Toint, Dec 1989.

#   classification NQR2-AY-3-3

option presolve 0;
param xinit{1..3};
var x{i in 1..3}:= xinit[i];
minimize f: 0;
subject to cons1:
	( -x[3]+0.25*(x[1]+x[2])^2 ) = 0;
subject to cons2:
	(-x[1]+1.0) = 0;
subject to cons3:
	(-x[2]+1.0) = 0;

data;
param xinit:= 1 -1.2 2 2.0 3 2.0;
solve;
display f;
display x; 
