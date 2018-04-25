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
#   R. Fletcher
#   "Practical Methods of Optimization",
#   second edition, Wiley, 1987.

#   SIF input: Ph. Toint, March 1994.

#   classification QOR2-AN-4-4

var x{1..4} := 1;

minimize f:
	x[1]*x[2];
subject to cons1:
	(x[1]*x[3]+x[2]*x[4])^2/(x[1]^2+x[2]^2) - x[3]^2 - x[4]^2 +
	1= 0;
subject to cons2:
	x[1]-x[3]-1 >= 0;
subject to cons3:
	x[2]-x[4]-1 >= 0;
subject to cons4:
	x[3]-x[4] >= 0;
subject to cons5:
	x[4] >= 1;

solve;
display f;
display x;
