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
#   Ph. Toint, private communication,

#   SIF input: Ph. Toint, Dec 1989.

#   classification OUR2-AY-2-0

param hlength:=30;
param cslope:=100;

var x1 := -5;
var x2 := -7;

minimize f:
	sin(7*x1)^2*cos(7*x2)^2*hlength +
cslope*sqrt(0.01+(x1-x2)^2) + cslope*sqrt(0.01+x1^2)
;

solve;
display f;
display x1, x2;
