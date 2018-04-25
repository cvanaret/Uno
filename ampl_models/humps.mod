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
#   Ph. Toint, private communication, 1997.

#   SDIF input: Ph. Toint, May 1997.

#   classification SUR2-AN-2-0

param zeta := 20.0;
var x:=-506.0;
var y:=-506.2;

minimize f:
	0.05*(x^2+y^2) + (sin(zeta*x)*sin(zeta*y))^2;

solve; display f; display x;
