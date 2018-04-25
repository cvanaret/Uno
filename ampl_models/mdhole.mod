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
#   Ph. Toint, private communication.

#   SIF input: Ph. Toint, March 1991.

#   classification OBR2-AN-2-0

var x:=10.0, >= 0;
var y:=10.0;

minimize f:
	((-y+sin(x))^2)/0.01 + x;

solve; display f; display x;
