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
#   a problem designed by Ph. Toint for experimenting with feasibility
#   issues in barrier approaches to nonlinear inequality constraints.

#   SIF input: Ph.L. Toint, September 93.

#   classification LOR2-AN-2-2

param tip := 0.0001;
var x:=1.0;
var y:=5.0;
minimize f:
	x;
subject to cons1:
	-y+sin(x) <= 0;
subject to cons2:
	-y+tip*x+sin(x) >= 0;

solve; display f; display x; display y;
