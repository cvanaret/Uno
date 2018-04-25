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

#   Source:  problem 4 in A. Zecevic, "Contribution to methods 
#            of external penalty functions - algorithm MNBS"
#            Advanced Business School, Belgrade, 
#            (whatever is left of) Yugoslavia.

#   SIF input: Nick Gould, April 1993.

#   classification QQR2-AN-2-2

param xinit{1..2};
var x{i in 1..2} := xinit[i];

minimize f:
	6*x[1]^2+x[2]^2-60*x[1]-8*x[2]+166;
subject to cons1:
	x[1]*x[2]-x[1]-x[2] <= 0;
subject to cons2:
	-x[1]-x[2] <= -3.0;
subject to cons3:
	0 <= x[1] <= 10;
subject to cons4:
	0 <= x[2] <= 10;

data;
param xinit := 1 0.1 2 -0.1;
solve; display f; display x;
