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

#   Source:
#   "The OPTIMA user manual (issue No.8, p. 26)",
#   Numerical Optimization Centre, Hatfield Polytechnic (UK), 1989.

#   SIF input: Ph. Toint, May 1990.

#   classification SBR2-AN-25-0

param N:=4;
var x{1..N} := 0.9;

minimize f:
	(x[1]-1)^2 + sum {i in 2..N-1} (x[i+1]-x[i]^2)^2 + (x[N]-1)^2;
subject to cons1{i in 1..N-1}:
	0.0 <= x[i] <= 10.0;
solve; display f; display x;
