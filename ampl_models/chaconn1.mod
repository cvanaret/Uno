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
#   C. Charalambous and A.R. Conn,
#   "An efficient method to solve the minmax problem directly",
#   SINUM 15, pp. 162-187, 1978.

#   SIF input: Ph. Toint, Nov 1993.

#   classification  LOR2-AY-3-3

param xinit{1..2};
var x{i in 1..2} := xinit[i];
var u;

minimize f:
	u;
subject to cons1:
	-u+x[1]^2+x[2]^4 <= 0;
subject to cons2:
	-u+(2-x[1])^2+(2-x[2])^2 <= 0;
subject to cons3:
	-u+2*exp(x[2]-x[1]) <= 0;

data;
param xinit:= 1 1.0 2 -0.1;

solve; display f; display x;
