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

#   Source: E. Polak and A. L. Tits,
#   "A recursive quadratic programming algorithm for semi-infinite
#    optimization problems",
#   Appl. Math. Optim. 8, 1982, pp 325-349.

#   SIF input: Nick Gould, February, 1994.

#   classification LLR2-AN-2-V

param M:=500;
param lower := 0.0;
param upper := 1.0;
param diff := upper-lower;
param h := diff/M;
var u;
var x;

minimize f:
	u;
subject to cons1{i in 0..M}:
	-(i*h+lower)+(i*h+lower)^2+u+ ( (i*h+lower)-3*(i*h+lower)^2 + 1)*x >= 0;

solve; display f; display x;
