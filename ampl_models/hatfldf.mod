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
#   "The OPTIMA user manual (issue No.8, p. 47)",
#   Numerical Optimization Centre, Hatfield Polytechnic (UK), 1989.

#   SIF input: Ph. Toint, May 1990.

#   classification NOR2-AN-3-3

var x{1..3} := 0.1;

minimize f: 0;
subject to cons1:
	(x[1]-0.032+x[2]*exp(x[3])) = 0;
subject to cons2:
	(x[1]-0.056+x[2]*exp(2*x[3])) = 0;
subject to cons3:
	(x[1]-0.099+x[2]*exp(3*x[3])) = 0;

solve; display f; display x;
