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

#   classification  LOR2-AN-3-2

var x{1..2} := 0.1;
var u;

minimize f:
	u;
subject to cons1:
	-u+3*x[1]^2+50*(x[1]-x[2]^4-1)^2 <= 0;
subject to cons2:
	-u+3*x[1]^2+50*(x[1]-x[2]^4+1)^2 <= 0;

solve;
display f;
display x;
