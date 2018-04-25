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

param p := 2000;
param k := 2;
var x{1..p+k};
minimize f:
	sum {i in 1..p+k} 0.5*(x[i]-sin(i/(p+k)))^2;
subject to cons1{i in 1..p}:
	x[i]-2*x[i+1]+x[i+2] >= 0;
subject to cons2:
	x[1] >= 0.08;
subject to cons3{i in p+1..p+k}:
	x[i] = 0.0;

solve; display f; display x;
