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

#   Source:  problem 30 in
#   J.J. More', B.S. Garbow and K.E. Hillstrom,
#   "Testing Unconstrained Optimization Software",
#   ACM Transactions on Mathematical Software, vol. 7(1), pp. 17-41, 1981.

#   See also Toint#17 and Buckley#78.
#   SIF input: Ph. Toint, Dec 1989.

#   classification NOR2-AN-V-V

param N:=10000;
param kappa1:=2.0;
param kappa2:=1.0;
var x{1..N} := -1.0;

minimize f: 0;
subject to cons1:
	(-2*x[2]+kappa2+(3-kappa1*x[1])*x[1]) = 0;
subject to cons2{i in 2..N-1}:
	(-x[i-1]-2*x[i+1]+kappa2+(3-kappa1*x[i])*x[i]) = 0;
subject to cons3:
	(-x[N-1]+kappa2+(3-kappa1*x[N])*x[N]) = 0;

solve; display f; display x;
