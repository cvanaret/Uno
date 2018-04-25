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

#   Source:  Problem 23 in
#   J.J. More', B.S. Garbow and K.E. Hillstrom,
#   "Testing Unconstrained Optimization Software",
#   ACM Transactions on Mathematical Software, vol. 7(1), pp. 17-41, 1981.

#   See also Buckley #181 (p. 79).

#   SIF input: Ph. Toint, Dec 1989.

#   classification SUR2-AN-V-0

param N:=1000;
param M:=N+1;

var x{i in 1..N} := i;

param a := 10^-5;

minimize f:
	sum {i in 1..N} a*(x[i]-1)^2 + ( sum {j in 1..N} x[j]^2 - 1/4 )^2;

solve;
display f;
display x;
