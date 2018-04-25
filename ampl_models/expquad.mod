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

#   classification OBR2-AN-V-V

param n := 120;
param m := 10;

var x{1..n} := 0.0;

minimize f:
	sum {i in 1..m} exp(0.1*i*m*x[i]*x[i+1]) +
	sum {i in m+1..n-1} (4.0*x[i]*x[i]+2.0*x[n]*x[n] + x[i]*x[n]) +
	sum {i in 1..n} (-10.0*i*x[i]);
subject to cons{i in 1..m}:
	0.0 <= x[i] <= 10.0;

solve;

	
