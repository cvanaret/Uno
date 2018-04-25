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
#   L.C.W. Dixon and Z. Maany,
#   "A family of test problems with sparse Hessians for unconstrained
#   optimization",
#   TR 206, Numerical Optimization Centre, Hatfield Polytechnic, 1988.

#   See also Buckley#221 (p. 49)
#   SIF input: Ph. Toint, Dec 1989.
#              correction by Ph. Shott, January 1995.

#   classification OUR2-AN-V-0

param M:=1000;
param N:=3*M;

param alpha:=1.0;
param beta:=0.125;
param gamma:=0.125;
param delta := 0.125;

param K{1..4} := 0;
var x{1..N}:=2.0;

minimize f:
	1.0 + sum {i in 1..N} alpha*x[i]^2*(i/N)^K[1] +
	sum {i in 1..N-1} beta*x[i]^2*(x[i+1]+x[i+1]^2)^2*(i/N)^K[2] +
	sum {i in 1..2*M} gamma*x[i]^2*x[i+M]^4*(i/N)^K[3] +
	sum {i in 1..M} delta*x[i]*x[i+2*M]*(i/N)^K[4];

solve;
display f;
display x;
