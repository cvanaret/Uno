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
#   Ph. Toint, private communication, 1991.

#   SDIF input: Ph. Toint, June 1993.

#   classification OBR2-AN-2-0

var x:=-1.2, >= 0;
var y:=1.0, >= 0;

minimize f:
	log(1+10000*(y-x^2)^2+(1-x)^2);

solve; display f; display x;
