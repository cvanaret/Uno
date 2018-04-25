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

#   Source: a generalization of problem 20 in
#   M. J. D. Powell
#   "On the quadratic progamming algorithm of Goldfarb and Idnani",
#   Mathematical Programmimg Study 25 (1985) 46-61.

#   SIF input: Nick Gould, August 1994.

#   classification QLR2-AN-V-V

param N:=1000;
var x{i in 1..N}:=if ( (i mod 2) == 1) then 0 else -0.5-i;

minimize f:
	0.5*sum {i in 1..N} x[i]^2;
subject to cons1{k in 1..N-1}:
	x[k+1]-x[k] >= -0.5+(-1)^k*k;
subject to cons2:
	x[1]-x[N] >= N-0.5;
solve;
