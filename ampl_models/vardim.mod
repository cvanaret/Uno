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

#   Source:  problem 25 in
#   J.J. More', B.S. Garbow and K.E. Hillstrom,
#   "Testing Unconstrained Optimization Software",
#   ACM Transactions on Mathematical Software, vol. 7(1), pp. 17-41, 1981.

#   See also Buckley#72 (p.98).

#   SIF input: Ph. Toint, Dec 1989.

#   classification  OUR2-AN-V-0

param N := 100;
var x{i in 1..N} := 1-i/N;

minimize f:
	sum {i in 1..N} (x[i]-1)^2 + (sum {i in 1..N} i*x[i] - N*(N+1)/2)^2 + (sum {i in 1..N} i*x[i] - N*(N+1)/2)^4;

solve; display f; display x;
