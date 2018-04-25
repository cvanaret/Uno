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
#   E. Polak, D.H. Mayne and J.E. Higgins,
#   "Superlinearly convergent algorithm for min-max problems"
#   JOTA 69, pp. 407-439, 1991.

#   SIF input: Ph. Toint, Nov 1993.

#   classification  LQR2-AN-3-3

param xinit{1..2};
var x{i in 1..2} := xinit[i];
var u;
minimize f:
	u;
subject to cons1:
	-u-x[1]-1+2*x[1]^2+2*x[2]^2 <=0;
subject to cons2:
	-u-0.01+0.01*x[1]^2+0.01*x[2]^2 <= 0;
subject to cons3:
	-u-100000.0+100000.0*(x[1]-2)^2+x[2]^2 <= 0;

data;
param xinit:= 1 0.9 2 0.1;

solve; display f; display x;
