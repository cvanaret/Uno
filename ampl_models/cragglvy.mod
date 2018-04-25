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

#   Source:  problem 32 in
#   Ph. L. Toint,
#   "Test problems for partially separable optimization and results
#   for the routine PSPMIN",
#   Report 83/4, Department of Mathematics, FUNDP (Namur, B), 1983.

#   See  also Buckley#18
#   SIF input: Ph. Toint, Dec 1989.

#   classification OUR2-AY-V-0

param m := 2499;
param n := 2*m+2;

var x{i in 1..n} := if (i==1) then 1.0 else 2.0;

minimize f:
	sum {i in 1..m} (
	(exp(x[2*i-1])-x[2*i])^4 +
	100*(x[2*i]-x[2*i+1])^6 +
	(tan(x[2*i+1]-x[2*i+2])+x[2*i+1]-x[2*i+2])^4 +
	(x[2*i-1])^8 +
	(x[2*i+2]-1.0)^2
	);

solve;
display f;
display x;
