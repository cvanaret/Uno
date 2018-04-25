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

#   Source:  problem 28 in
#   J.J. More', B.S. Garbow and K.E. Hillstrom,
#   "Testing Unconstrained Optimization Software",
#   ACM Transactions on Mathematical Software, vol. 7(1), pp. 17-41, 1981.

#   See also Buckley#17 (p. 75).

#   SIF input: Ph. Toint, Dec 1989 and Nick Gould, Oct 1992.

#   classification SUR2-MN-V-0

param N:=5000;

param h:=1/(N+1);
param t{i in 1..N} := i*h;

var x{i in 0..N+1} := if (1 <= i <= N) then t[i]*(t[i]-1);

minimize f:
	sum {i in 1..N} (2*x[i]-x[i-1]-x[i+1]+h^2*(x[i]+t[i]+1)^3/2)^2;

subject to cons1: x[0] = 0;
subject to cons2: x[N+1] = 0;

solve;
display f;
display x;
