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

#   classification QQR2-AN-2-1

var x >= 0, := 10.0;
var y >= 0, := 10.0;

minimize f:
	(x-1)^2;
subject to cons1:
	-1+(x-1)^2+(y-10)^2 = 0;

solve;
display f;
display x;
display y;
