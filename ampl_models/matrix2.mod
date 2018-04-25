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

#   classification  QOR2-AY-6-2

var x11 >= 0, := 1;
var x12 := 1;
var x22 >= 0, := 1;
var y11 <= 0, := 1;
var y12 := 1;
var y22 <= 0, := 1;

minimize f:
	(x11-y11)^2 + 2*(x12-y12)^2 + (x22-y22)^2;
subject to cons1:
	x11*x22-x12^2 >= 0;
subject to cons2:
	y11*y22-y12^2<= 0;

solve;
display f;
display x11,x12,x22,y11,y12,y22;
