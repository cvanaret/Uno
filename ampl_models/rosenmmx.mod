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
#   M.M. Makela,
#   "Nonsmooth optimization",
#   Ph.D. thesis, Jyvaskyla University, 1990

#   SIF input: Ph. Toint, Nov 1993.

#   classification  LQR2-AN-5-4

var x{1..4} := 0.0;
var u := 0.0;
minimize f:
	u;
subject to cons1:
	-u-5*x[1]-5*x[2]-21*x[3]+7*x[4]+x[1]^2+x[2]^2+2*x[3]^2+x[4]^2 <= 0;
subject to cons2:
	-u+5*x[1]-15*x[2]-11*x[3]-3*x[4]-80+11*sum {i in 1..4} x[i]^2 +x[3]^2<= 0;
subject to cons3:
	-u-15*x[1]-5*x[2]-21*x[3]-3*x[4]-100+11*x[1]^2+21*x[2]^2+12*x[3]^2+21*x[4]^2 <= 0;
subject to cons4:
	-u+15*x[1]-15*x[2]-21*x[3]-3*x[4]-50+11*(x[1]^2+x[2]^2)+12*x[3]^2+x[4]^2 <= 0;

solve;
display f;
display x;

