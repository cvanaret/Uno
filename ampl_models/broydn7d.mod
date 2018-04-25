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
#   Ph.L. Toint,
#   "Some numerical results using a sparse matrix updating formula in
#   unconstrained optimization",
#   Mathematics of Computation, vol. 32(114), pp. 839-852, 1978.

#   See also Buckley#84
#   SIF input: Ph. Toint, Dec 1989.

#   classification OUR2-AN-V-0

param N:=1000;

var x{1..N} := 1.0;

minimize f:
	abs(-2*x[2]+1+(3-2*x[1])*x[1])^(7/3) +
	sum {i in 2..N-1} abs(1-x[i-1]-2*x[i+1]+(3-2*x[i])*x[i])^(7/3)+
	abs(-x[N-1]+1 +(3-2*x[N])*x[N])^(7/3) +
	sum {i in 1..N/2} abs(x[i]+x[i+N/2])^(7/3);

solve;
display f;
display x;
