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

#   SIF input: Ph. Toint, Jan 1991.

#   classification LLR2-AN-2-2

var x >= 0;
var y;
minimize f:
	x;
subject to cons1:
	x+2*y-2 = 0;
subject to cons2:
	2*x+y-2 = 0;

solve;
display f;
display x,y;
