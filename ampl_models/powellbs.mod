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

#   Source:  Problem 3 in
#   J.J. More', B.S. Garbow and K.E. Hillstrom,
#   "Testing Unconstrained Optimization Software",
#   ACM Transactions on Mathematical Software, vol. 7(1), pp. 17-41, 1981.

#   See also  Toint#34, Buckley#22 (p. 82).

#   SIF input: Ph. Toint, Dec 1989.

#   classification NOR2-AN-2-2

param N:=2;
param xinit{1..N};
var x{i in 1..N}:=xinit[i];

minimize f:0;
subject to cons{i in 1..N-1}:
	(-1.0+10000*x[i]*x[i+1]) = 0;
subject to cons2{i in 1..N-1}:
	(-1.0001+exp(-x[i])+exp(-x[i+1])) = 0;

data;
param xinit:= 1 0.0 2 1.0;

solve; display f; display x;
