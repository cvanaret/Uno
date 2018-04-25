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

#   classification  LOR2-AN-5-4

var x{1..4} := 0;
var u;

minimize f:
	u;
subject to cons1:
	(x[1]-(x[4]+1)^4)^2 + (x[2]-(x[1]-(x[4]+1)^4)^4)^2 +
	2*x[3]^2 + x[4]^2 + 5*(x[4]+1)^4 + 5*(x[1]-(x[4]+1)^4)^4 
	-u-5*x[1]-5*x[2]-21*x[3]+7*x[4] <= 0;
subject to cons2:
	11*(x[1]-(x[4]+1)^4)^2 + 11*(x[2]-(x[1]-(x[4]+1)^4)^4)^2 
	+ 12*x[3]^2 + 11*x[4]^2 - 5*(x[4]+1)^4 + 
	15*(x[1]-(x[4]+1)^4)^4-u+5*x[1]-15*x[2]-11*x[3]-3*x[4]-80
	<= 0; 
subject to cons3:
	11*(x[1]-(x[4]+1)^4)^2 + 21*(x[2]-(x[1]-(x[4]+1)^4)^4)^2 
	+ 12*x[3]^2 + 21*x[4]^2 + 15*(x[4]+1)^4
	+5*(x[1]-(x[4]+1)^4)^4-u-15*x[1]-5*x[2]-21*x[3]-3*x[4]-100
	<= 0; 
subject to cons4:
	11*(x[1]-(x[4]+1)^4)^2 + 11*(x[2]-(x[1]-(x[4]+1)^4)^4)^2 
	+ 12*x[3]^2 + x[4]^2 - 15*(x[4]+1)^4 + 15*(x[1]-(x[4]+1)^4)^4 
	-u+15*x[1]-15*x[2]-21*x[3]-3*x[4]-50 <= 0;

solve;
display f;
display x;
