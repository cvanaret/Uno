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
#   M. Batholomew-Biggs and F.G. Hernandez,
#   "Some improvements to the subroutine OPALQP for dealing with large
#    problems",
#   Numerical Optimization Centre, Hatfield, 1992.

#   SIF input: Ph Toint, April 1992.

#   classification QLR2-AN-4-7

var x{1..4} >= 0.0, <= 5.0, := 0.0;

minimize f:
	(-x[1]*x[3]-x[2]*x[4]);
subject to cons1:
	0 <= x[1]+x[2] -2.5 <= 5.0;
subject to cons2:
	0 <= x[1]+x[3] -2.5 <= 5.0;
subject to cons3:
	0 <= x[1]+x[4] -2.5 <= 5.0;
subject to cons4:
	0 <= x[2]+x[3] -2.0 <= 5.0;
subject to cons5:
	0 <= x[2]+x[4] -2.0 <= 5.0;
subject to cons6:
	0 <= x[3]+x[4] -1.5 <= 5.0;
subject to cons7:
	x[1]+x[2]+x[3]+x[4]-5.0 >= 0;

solve; display f; display x;
