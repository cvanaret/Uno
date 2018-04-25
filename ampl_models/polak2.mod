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

#   classification  LOR2-AN-11-2

var x{i in 1..10} := if (i == 1) then 100 else 0.1;
var u;

minimize f:
	u;
subject to cons1:
	exp(1.0D-8*x[1]^2+(x[2]+2)^2+x[3]^2+4*x[4]^2+sum {i in
	5..10} x[i]^2 ) -u <= 0;
subject to cons2:
	exp(1.0D-8*x[1]^2+(x[2]-2)^2+x[3]^2+4*x[4]^2+sum {i in
	5..10} x[i]^2 ) -u <= 0;

solve;
display f;
display x;
