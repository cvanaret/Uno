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
#   a contribution to fullfill the LANCELOT academic licence agreement,
#   inspired by
#   Y. Lin and J. Pang,
#   "Iterative methods for large convex quadratic programs: a survey",
#   SIAM Journal on Control and Optimization 25, pp.383-411, 1987.

#   SIF input: J. Judice, University of Coimbra, January 1995.
#              condensed by Ph. Toint, January 1995.

#   classification QBR2-AN-V-0

param N:=1000;
var x{1..N} >= 0;

minimize f:
	sum {i in 1..N} 6*x[i]^2 - 3*x[1]+x[2]+x[0.5*N-1]-3*x[N/2] + 4*x[N/2+1] + sum {i in N/2+3..N} x[i] + sum {i in 1..N-2}
(-4*x[i]*x[i+1]+x[i]*x[i+2]);

solve; display f; 
#display x;
