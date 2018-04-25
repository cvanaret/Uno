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
#   A.R. Conn, N. Gould and Ph.L. Toint,
#   "The LANCELOT User's Manual",
#   Dept of Maths, FUNDP, 1991.

#   SIF input:  Ph. Toint, Jan 1991.

#   classification OOR2-AN-2-1

var x;
var y;

minimize f:
	exp(x-2*y);

subject to cons1:
	sin(-x+y-1)=0;
subject to cons2:
	-2 <= x <= 2;
subject to cons3:
	-1.5 <= y <= 1.5;

solve;
display f;
display x, y;
