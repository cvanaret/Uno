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
#   R.  Wommersley and R. Fletcher,
#   "An algorithm for composite nonsmooth optimization problems"
#   JOTA, vol.48, pp.493-523, 1986

#   SIF input: Ph. Toint, April 1992.

#   classification LOR2-AN-3-3

var x{1..2} := 2.0;
var u := 1.0;

minimize f:
	u;
subject to cons1:
	u-x[1]^4-x[2]^2 >= 0;
subject to cons2:
	u-(2.0-x[1])^2-(2.0-x[2])^2 >= 0;
subject to cons3:
	u-2*exp(x[2]-x[1]) >= 0;

solve; display f; display x;
