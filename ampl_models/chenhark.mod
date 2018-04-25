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
#   B. Chen and P. T. Harker,
#   SIMAX 14 (1993) 1168-1190

#   SDIF input: Nick Gould, November 1993.

#   classification QBR2-AN-V-V

param n:=1000;
param nfree := 500;
param ndegen := 200;
param x_p{i in -1..n+2} := if (i <= 0 || i > nfree) then 0.0 else 1.0;

var x{i in 1..n} >= 0, := 0.5;

minimize f:
	sum {i in 2..n-1} 0.5*(x[i+1]+x[i-1] - 2*x[i])^2 +
	0.5*x[1]^2 +
	0.5*(2*x[1] - x[2])^2 +
	0.5*(2*x[n] - x[n-1])^2 +
	0.5*(x[n])^2 +
	sum {i in 1..nfree+ndegen} x[i]*(-6*x_p[i] + 4*x_p[i+1] + 
		4*x_p[i-1] -
		x_p[i+2] - x_p[i-2]) +
	sum {i in nfree+ndegen+1..n} x[i]*(-6*x_p[i] + 4*x_p[i+1] + 4*x_p[i-1] -
                x_p[i+2] - x_p[i-2] + 1)
	;

solve;
