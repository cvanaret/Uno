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

#   Source: Problem 27 in
#   J.J. More', B.S. Garbow and K.E. Hillstrom,
#   "Testing Unconstrained Optimization Software",
#   ACM Transactions on Mathematical Software, vol. 7(1), pp. 17-41, 1981.

#   See also Buckley#79
#   SIF input: Ph. Toint, Dec 1989.

#   classification SUR2-AN-V-0

param N := 10;

var x{1..N} := 1/2;

minimize f:
	sum {i in 1..N-1} (x[i]+sum {j in 1..N} x[j] - (N+1) )^2 + (prod {j in 1..N} x[j] - 1)^2;

solve;
display f;
display x;

