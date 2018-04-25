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

#   Source: problem 33 in
#   D.H. Himmelblau,
#   "Applied nonlinear programming",
#   McGraw-Hill, New-York, 1972.

#   See Buckley#87 (p. 67)

#   SIF input: Ph. Toint, Dec 1989.

#   classification OUR2-AN-2-0

var x{1..2} := 0.5;

minimize f:
	exp(-x[1]-x[2])*(2*x[1]^2+3*x[2]^2);

solve; display f; display x;
