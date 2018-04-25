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
#   R.S. Wommersley and R. Fletcher,
#   "An algorithm for composite nonsmooth optimization problems",
#   JOTA,, vol.48, pp. 493-523, 1986.

#   SIF input: Ph. Toint, April 1992.

#   classification LOR2-AN-3-3

param xinit{1..2};
var x{i in 1..2} := xinit[i];
var u := 7.5;

minimize f:
	u;
subject to cons1:
	u-0.5*x[1]-x[2]^2-5*x[1]/(x[1]+0.1) >= 0;
subject to cons2:
	u+0.5*x[1]-x[2]^2-5*x[1]/(x[1]+0.1) >= 0;
subject to cons3:
	u+0.5*x[1]+x[2]^2+5*x[1]/(x[1]+0.1) >= 0;

data;
param xinit:= 1 3.0 2 1.0;

solve; display f; display x;
