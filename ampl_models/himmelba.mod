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

#   Source: problem 25 in
#   D.H. Himmelblau,
#   "Applied nonlinear programming",
#   McGraw-Hill, New-York, 1972.

#   See Buckley#215 (p. 61)

#   SIF input: Ph. Toint, Dec 1989.

#   classification NLR2-AN-2-2

option presolve 0;
param xinit{1..2};
var x{i in 1..2} := xinit[i];
minimize f: 0;
subject to cons1:
	 (4*x[1]-20) = 0;
subject to cons2:
	 (x[2]-6) = 0;

data;
param xinit:= 1 8.0 2 9.0;

solve;
display f;
display x;
