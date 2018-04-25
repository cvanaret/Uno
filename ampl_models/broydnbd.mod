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

#   Source: problem 31 in
#   J.J. More', B.S. Garbow and K.E. Hillstrom,
#   "Testing Unconstrained Optimization Software",
#   ACM Transactions on Mathematical Software, vol. 7(1), pp. 17-41, 1981.

#   See also Buckley#73 and Toint#18
#   SIF input: Ph. Toint, Dec 1989.

#   classification NOR2-AN-V-V

param N:=5000;

param ml := 5;
param mu := 1;

set J{i in 1..N} := {j in 1..N : j != i && max(1,i-ml) <= j <= min(N,i+mu) };

var x{1..N} := -1;

minimize f: 0;
subject to cons{i in 1..N}:
	( x[i]*(2+5*x[i]^2) + 1 - sum {j in J[i]} x[j]*(1+x[j]) ) = 0;

solve;
display f;
display x;
