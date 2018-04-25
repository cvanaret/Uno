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

#   Source: problem 2 in
#   J.J. More', B.S. Garbow and K.E. Hillstrom,
#   "Testing Unconstrained Optimization Software",
#   ACM Transactions on Mathematical Software, vol. 7(1), pp. 17-41, 1981.

#   See also Toint#33, Buckley#24
#   SIF input: Ph. Toint, Dec 1989.

#   classification SUR2-AN-V-0

param n := 5000;
param ngs := n-1;

var x{1..n};

minimize f:
	sum {i in 1..ngs} ((5.0-x[i+1])*x[i+1]^2+x[i]-2*x[i+1]-13.0)^2 +
	sum {i in 1..ngs} ((1.0+x[i+1])*x[i+1]^2+x[i]-14*x[i+1]-29.0)^2;

let x[1] := 0.5;
let x[2] := -2.0;
let {i in 1..n: i>2} x[i] := 0.0;

solve;
display f;
display x;
