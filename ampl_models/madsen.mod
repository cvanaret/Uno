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
#   K. Madsen
#   "An algorithm for minmax solution of overdetermined systems of non-linear
#   equations"
#   JIMA, vol.16, pp. 321-328, 1975.

#   SIF input: Ph. Toint, April 1992.

#   classification LOR2-AN-3-6

var x1 := 3;
var x2 := 1;
var u := 1;

minimize f:
	u;
subject to cons1:
	u-x1^2-x2^2-x1*x2 >= 0;
subject to cons2:
	u+x1^2+x2^2+x1*x2 >= 0;
subject to cons3:
	u-sin(x1) >= 0;
subject to cons4:
	u+sin(x1) >= 0;
subject to cons5:
	u-cos(x2) >= 0;
subject to cons6:
	u+cos(x2) >= 0;

solve;
display f;
display x1,x2,u;
