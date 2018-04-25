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

#   Source: a very simple box-constrained quadratic

#   SIF input: Nick Gould, March 1992

#   classification QBR2-AN-2-0

param xinit{1..2};
var x{i in 1..2} := xinit[i];

minimize f:
	x[2]+(-x[1]+x[2])^2+(2*x[1]+x[2])^2;
subject to cons1:
	0.0 <= x[2] <= 0.5;

data;
param xinit:= 1 10 2 1;

solve;
display f;
display x;
	
