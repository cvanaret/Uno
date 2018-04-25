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

#   Source: problem 46 in
#   Ph.L. Toint,
#   "Test problems for partially separable optimization and results
#   for the routine PSPMIN",
#   Report 83/4, Department of Mathematics, FUNDP (Namur, B), 1983.

#   SIF input: Ph. Toint, Dec 1989.

#   classification OBR2-AN-V-0

param N := 500;

var x{1..N} >= -5, <= 5, := 0.5;

minimize f:
	sum {i in 1..N-1} (x[i]+x[i+1])*0.0001*exp(-x[i]*x[i+1])/N
	+ 100*(sum {i in 1..N} x[i] - 1)^2;

solve;
display f;
display x;
