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

#   Source: problem 29 in
#   D.H. Himmelblau,
#   "Applied nonlinear programming",
#   McGraw-Hill, New-York, 1972.

#   See Buckley#86 (p. 64)

#   SIF input: Ph. Toint, Dec 1989.

#   classification NQR2-AN-2-2

var x{1..2} := 1.0;
minimize f: 0;
subject to cons1:
	(12*x[2]-1.0+x[1]^2) = 0;
subject to cons2:
	(84*x[1]+2324*x[2]-681.0+49*x[1]^2+49*x[2]^2) = 0;

solve;
display f;
display x;

