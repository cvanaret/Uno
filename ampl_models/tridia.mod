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

#   Source: problem 8 in
#   Ph.L. Toint,
#   "Test problems for partially separable optimization and results
#   for the routine PSPMIN",
#   Report 83/4, Department of Mathematics, FUNDP (Namur, B), 1983.

#   See also Buckley#40 (p.96)

#   SIF input: Ph. Toint, Dec 1989.

#   classification QUR2-AN-V-0

param N;
param alpha;
param beta;
param gamma;
param delta;

var x{1..N}:=1.0;

minimize f:
	gamma*(x[1]*delta-1.0)^2 + sum {i in 2..N} i*(-beta*x[i-1]+alpha*x[i])^2;

data;
param N:=10000;
param alpha:=2.0;
param beta:=1.0;
param gamma:=1.0;
param delta:=1.0;

solve; display x;

