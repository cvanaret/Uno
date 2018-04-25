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

#   Source:  Problem 29 in
#   J.J. More', B.S. Garbow and K.E. Hillstrom,
#   "Testing Unconstrained Optimization Software",
#   ACM Transactions on Mathematical Software, vol. 7(1), pp. 17-41, 1981.

#   SIF input: Ph. Toint, Feb 1990.

#   classification NOR2-AN-V-V

param N:=100;

param h := 1/(N+1);
param t{i in 1..N} := i*h;
var x{i in 0..N+1} := if (1 <= i <= N) then t[i]*(t[i]-1);

minimize f: 0;
subject to cons{i in 1..N}:
	 ( x[i]+h*( (1-t[i])*sum {j in 1..i} t[j]*(x[j]+t[j]+1)^3 + t[i]*sum {j in i+1..N}
	(1-t[j])*(x[j]+t[j]+1)^3)/2 ) = 0;

fix x[0] := 0;
fix x[N+1] := 0;

solve;
display f;
display x;
