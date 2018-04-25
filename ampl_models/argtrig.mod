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

#   Source:  Problem 26 in
#   J.J. More', B.S. Garbow and K.E. Hillstrom,
#   "Testing Unconstrained Optimization Software",
#   ACM Transactions on Mathematical Software, vol. 7(1), pp. 17-41, 1981.

#   SIF input: Ph. Toint, Dec 1989.

#   classification NOR2-AN-V-V

param N:=100;
var x{1..N} := 1/N;

minimize f: 0;
subject to cons{i in 1..N}:
	(i*(cos(x[i])+sin(x[i])) + sum {j in 1..N} cos(x[j]) - (N+i) ) = 0;

solve; display f; display x;
