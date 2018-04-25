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

#   Source:  problem 14 in
#   J.J. More', B.S. Garbow and K.E. Hillstrom,
#   "Testing Unconstrained Optimization Software",
#   ACM Transactions on Mathematical Software, vol. 7(1), pp. 17-41, 1981.

#   See also Toint#27, Buckley#17 (p. 101), Conn, Gould, Toint#7

#   SIF input: Ph. Toint, Dec 1989.

#   classification SUR2-AN-V-0

var x{i in 1..10000} := if (i mod 2 == 1) then -3 else -1;

minimize f:
	sum {i in 1..2500} (100*(x[4*i-2]-x[4*i-3]^2)^2 + (1-x[4*i-3])^2 + 90*(x[4*i]-x[4*i-1]^2)^2 + (1-x[4*i-1])^2 +
	10*(x[4*i-2]+x[4*i]-2)^2 + 0.1*(x[4*i-2]-x[4*i])^2 );

solve;
display f;
display x;

