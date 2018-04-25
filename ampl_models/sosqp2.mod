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

#   classification QLR2-AN-V-V

param N:=10000;
var x{1..N} <= 1, >= -1;
var y{1..N} <= 1, >= -1;

minimize f:
	sum {i in 1..N} x[i]*y[i];
subject to cons1{i in 1..N}:
	i*x[i]-y[i] = i;
subject to cons2:
	sum {i in 1..N} (x[i]+y[i]) = N/2;

solve; display f; display x, y;
