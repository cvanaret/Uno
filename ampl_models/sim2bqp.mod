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

#   Source: another very simple box-constrained quadratic

#   classification QBR2-AN-2-0

var x1 := 10;
var x2 := 1, >= 0, <= 0.5;

minimize f:
	x2
	+ (-x1+x2)^2
	+ (x1+x2)^2
	;

solve;
display f;
display x1, x2;
