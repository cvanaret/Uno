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
#   N. Gould: private communication.

#   SIF input: Nick Gould, June 1990.

#   classification OOR2-AY-4-1

var x{1..4};

minimize f:
	x[3]-1 +
	x[1]^2+
	x[2]^2 + (x[3]+x[4])^2 +
	sin(x[3])^2 + x[1]^2*x[2]^2 + x[4]-3 +
	sin(x[3])^2 +
	(x[4]-1)^2 +
	(x[2]^2)^2+
	(x[3]^2 + (x[4]+x[1])^2)^2 +
	(x[1]-4 + sin(x[4])^2 + x[2]^2*x[3]^2)^2 +
	sin(x[4])^4;
subject to cons1:
	x[2] >= 1;
subject to cons2:
	-1D+10 <= x[3] <= 1;
subject to cons3:
	x[4] = 2;
subject to cons4:
	x[1]^2+x[2]^2 -1 <= 0;
solve;
display f;
display x;
