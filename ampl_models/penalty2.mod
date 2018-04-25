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

#   Source:  Problem 24 in
#   J.J. More', B.S. Garbow and K.E. Hillstrom,
#   "Testing Unconstrained Optimization Software",
#   ACM Transactions on Mathematical Software, vol. 7(1), pp. 17-41, 1981.

#   See also Buckley#112 (p. 80)

#   SIF input: Ph. Toint, Dec 1989.

#   classification SUR2-AN-V-0

param N:=100;
param M := 2*N;

param a := 10^-5;
var x{1..N} := 1/2;
param y{i in 1..M} := exp(i/10)  + exp( (i-1)/10 );

minimize f:
	(x[1]-0.2)^2 + sum {i in 2..N} a*(exp(x[i]/10)+exp(x[i-1]/10)-y[i])^2 + sum {i in N+1..2*N-1}
	a*(exp(x[i-N+1]/10)-exp(-1/10))^2 + ( sum {j in 1..N} (N-j+1)*x[j]^2 - 1)^2;

solve;
display f;
display x;
