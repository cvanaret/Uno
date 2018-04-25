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

#   classification OBR2-MN-V-0

param N := 120;
param M := 10;

var x{1..N} := 0.0;

minimize f:
	sum {i in 1..M} ((i/M)*(x[i]*x[i+1])^4-10*i*x[i])
	+ sum {i in M+1..N-1} (4*x[i]^2+2*x[N]^2+x[i]*x[N]-10*i*x[i]);
subject to cons{i in 1..M}:
	0.0 <= x[i] <= 10.0;

solve;
display f;
display x;
