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

#   Source: problem 47 in
#   Ph.L. Toint,
#   "Test problems for partially separable optimization and results
#   for the routine PSPMIN",
#   Report 83/4, Department of Mathematics, FUNDP (Namur, B), 1983.

#   SIF input: Ph. Toint, Dec 1989.

#   classification SBR2-AY-4-0

param N:=4;
param NGS := N-2;

var x{1..N} := 3;

minimize f:
	sum {i in 1..NGS} sqrt(x[i]^2+(x[i+1]-x[i+2])^2+1) ;
subject to cons1:
	x[1] <= -1;

solve;
display f;
display x;
