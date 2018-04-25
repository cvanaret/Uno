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
#   M.J.D. Powell,
#   " A hybrid method for nonlinear equations",
#   In P. Rabinowitz(ed.) "Numerical Methods for Nonlinear Algebraic
#   Equations", Gordon and Breach, 1970.

#   See also Buckley#217 (p.84.)

#   classification NOR2-AN-2-2

param xinit{1..2};
var x{i in 1..2} := xinit[i];

minimize f: 0;
subject to cons1:
	(x[1]^2) = 0;
subject to cons2:
	(10*x[1]/(x[1]+0.1)+2*x[2]^2) = 0;

data;
param xinit:= 1 3 2 1;

solve;
display f;
display x;

