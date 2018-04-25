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
#   a contribution to fullfill the LANCELOT academic licence agreement.

#   SIF input: Li-zhi Liao, Dept. of Mathematics, 
#              Hong Kong Baptist College, May 1994.

#   classification OQR2-AN-3-2

param xinit{1..3};
var x{i in 1..3} := xinit[i], >= 0.0;

minimize f:
	11*x[1]+x[3]+x[1]^3-6*x[1]^2;
subject to cons1:	
	-x[1]^2-x[2]^2+x[3]^2 >= 0;
subject to cons2:
	-4+x[1]^2+x[2]^2+x[3]^2 >= 0;
subject to cons3:
	x[3] <= 5.0;

data;
param xinit:= 1 0.0 2 0.0 3 3.0;

solve; display f; display x;
