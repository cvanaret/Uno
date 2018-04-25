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

#   Source:  problem 21 in
#   J.J. More', B.S. Garbow and K.E. Hillstrom,
#   "Testing Unconstrained Optimization Software",
#   ACM Transactions on Mathematical Software, vol. 7(1), pp. 17-41, 1981.

#   SIF input: Ph. Toint, Dec 1989.

#   classification SUR2-AN-V-0

param N:=10000;

var x{i in 1..N} := if (i mod 2 == 1) then -1.2 else 1;

minimize f:
	sum {i in 1..N/2} ( 100*(x[2*i]-x[2*i-1]^2)^2 + (x[2*i-1]-1)^2 );

solve;
display f;
display x;
